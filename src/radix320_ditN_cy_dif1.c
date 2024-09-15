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
#include "radix64.h"

#define RADIX 320	// Use #define rather than const int to ensure it's really a compile-time const in the C sense
#define ODD_RADIX 5	// ODD_RADIX = [radix >> trailz(radix)]

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

	#include "sse2_macro_gcc64.h"

  // For Mersenne-mod need (16 [SSE2] or 64 [AVX]) + (4 [HIACC] or 40 [LOACC]) added slots for half_arr lookup tables.
  // Max = (40 [SSE2]; 132 [AVX]), add to (half_arr_offset320[=0x508] + RADIX) to get required value of radix320_creals_in_local_store:
  #ifdef USE_AVX512	// RADIX/8 = 0x18 fewer carry slots than AVX:
	const int half_arr_offset320 = 0x508 + (RADIX>>3);	// RADIX/8 vec_dbl slots for carries in AVX-512 mode
	// (half_arr_offset320 + RADIX) + 132 (=0x84) and round up to nearest multiple of 4:
	const int radix320_creals_in_local_store = (0x508 + (RADIX>>3) + RADIX + 132 + 3)&(~0x3);
  #elif defined(USE_AVX)
	const int half_arr_offset320 = 0x508 + (RADIX>>2);	// RADIX/4 vec_dbl slots for carries in AVX mode
	// (half_arr_offset320 + RADIX) + 132 (=0x84) and round up to nearest multiple of 4:
	const int radix320_creals_in_local_store = (0x508 + (RADIX>>2) + RADIX + 132 + 3)&(~0x3);
  #else
	const int half_arr_offset320 = 0x508 + (RADIX>>1);	// RADIX/2 vec_dbl slots for carries in SSE2 mode
	// (half_arr_offset320 + RADIX) = 40 and round up to nearest multiple of 4:
	const int radix320_creals_in_local_store = (0x508 + (RADIX>>1) + RADIX + 40 + 3)&(~0x3);
  #endif

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
		vec_dbl *r00;
		vec_dbl *half_arr;
	#else
		double *r00;
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

/***************/

int radix320_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-320 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-320 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix320_ditN_cy_dif1";
  #ifdef USE_SSE2
	static int thr_id = 0;	// Master thread gets this special id
  #endif
  #if !defined(MULTITHREAD) && defined(USE_SSE2)
	const int pfetch_dist = PFETCH_DIST;
  #endif
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
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
	int j2,jp,ntmp;
  #endif
  #ifndef MULTITHREAD
	int jstart, ilo,ihi,hi_neg;;
  #endif
  #ifdef USE_SSE2
	int nbytes;
  #endif
  #ifndef MULTITHREAD
	int k0,k1,k2,k3,k4,nshift;
  #endif
  #if !defined(MULTITHREAD) && defined(USE_SSE2)
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 20|40|80 for avx512,avx,sse, respectively:
	// Note: Not able to get away with as large a setting of incr for radix 320 as we did for 160.
	int incr;
  #ifdef USE_AVX512
	const int incr_long = 10,incr_med = 5,incr_short = 5;
  #elif defined(USE_AVX2)
	const int incr_long = 20,incr_med =10,incr_short = 4;
  #else
	const int incr_long = 20,incr_med =10,incr_short = 4;
  #endif
  // Have no specialized HIACC carry macro in USE_AVX512 and ARMv8 SIMD, so use "goes to 11" in LOACC mode via an incr_hiacc = 1:
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

	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,pg0,ph0,pi0,pj0, nsave = 0;
	static int poff[RADIX>>2];	// Store mults of p4 offset for loop control
#ifndef MULTITHREAD
// Shared DIF+DIT:
  #ifndef USE_SSE2
	double rt,it;
  #endif
	static int t_offsets[64];
	static int plo[16];
  #ifndef USE_SSE2
	static int phi[20];
  #endif
	// Need storage for 4 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 4*9 elts:
	static int dit_p20_cperms[36];
	static int dif_o_offsets[RADIX], dit_i_offsets[RADIX];
	// Consts containing [in little-endian hex-char fashion] the 6 [A-F] 16-subperms defined in the output-ordering commentary below:
	const uint64 dif_16_operms[6] = {
		0x01326745cdfe98abull,	// [A]
		0x98abfedc32104576ull,	// [B]
		0x32104576fedcab89ull,	// [C]
		0x6745321098abfedcull,	// [D]
		0xcdfe98ab67453210ull,	// [E]
		0xfedcab8945761023ull	// [F]
	};
	// Consts containing [in little-endian hex-char fashion] the 9 [A-I] 16-subperms defined in the Input-ordering commentary below:
	const uint64 dit_16_iperms[9] = {
		0x01327654fedcba98ull,	// [A]
		0xfedcba9876543210ull,	// [B]
		0x32105467dcef98abull,	// [C]
		0x98abefcd10236745ull,	// [D]
		0x10236745efcdab89ull,	// [E]
		0x67452301ab89cdfeull,	// [F]
		0xab89cdfe23014576ull,	// [G]
		0xcdfe89ba45760132ull,	// [H]
		0x4576013289bafedcull	// [I]
	};
	uint64 i64;
#endif
	static double radix_inv, n2inv;
#if defined(USE_SSE2) || !defined(MULTITHREAD)
	const double cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
				 cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
				 s2  =  0.95105651629515357211,	/*  sin(u) */
				 ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
				 ss2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
#endif
	double scale, dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
#ifndef MULTITHREAD
	double *addr;
#endif
	struct complex t[RADIX];
#if !defined(MULTITHREAD) && !defined(USE_SSE2)
	struct complex *tptr;
#endif
#ifndef MULTITHREAD
	int *itmp;	// Pointer into the bjmodn array
#endif
#if !defined(MULTITHREAD) && defined(USE_AVX) && !defined(USE_AVX512)
	int *itm2;	// Pointer into the bjmodn array
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
#ifndef MULTITHREAD
	int col,co2,co3;
  #ifdef USE_AVX512
	double t0,t1,t2,t3;
   #ifdef CARRY_16_WAY
	static struct uint32x16 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
   #else
	static struct uint32x8  *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
   #endif
  #elif defined(USE_AVX)
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
  #endif
#endif // !MULTITHREAD

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
	double *add0,*add1,*add2,*add3;	/* Addresses into array sections */
  #endif

  #ifndef MULTITHREAD
	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
  #endif
  #ifndef USE_AVX512
	const double crnd = 3.0*0x4000000*0x2000000;
  #endif
  #ifndef USE_AVX
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
  #endif
	vec_dbl
	#ifndef MULTITHREAD
		*va0,*va1,*va2,*va3,*va4, *vb0,*vb1,*vb2,*vb3,*vb4,
	  #if defined(USE_AVX2) && !defined(USE_AVX512)
		*wa0,*wa1,*wa2,*wa3,*wa4, *wb0,*wb1,*wb2,*wb3,*wb4,
	  #endif
	#endif
		*tmp,*tm2	// Non-static utility ptrs
	#ifndef USE_SSE2
		,*tm0	// Non-static utility ptrs
	#endif
	#if !defined(MULTITHREAD) && defined(USE_SSE2)
		,*tm1	// Non-static utility ptrs
	#endif
		;
	static vec_dbl *two,
		*ycc1,*yss1,*ycc2,*yss2,*yss3,	// radiy-5 DFT trig consts
		*max_err, *sse2_rnd, *half_arr,
		*r00,	// Head of RADIX*vec_cmplx-sized local store #1
	  #if !defined(MULTITHREAD) && defined(USE_SSE2)
		*s1p00,	// Head of RADIX*vec_cmplx-sized local store #2
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
	static task_control_t   task_control = {NULL, (void*)cy320_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	int m,m2;
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	int bjmodn[RADIX];
	double temp,frac,cy[RADIX];

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

#ifndef MULTITHREAD
	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=-1;
#endif
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
		// consisting of radix320_creals_in_local_store dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix320_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix320_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = sc_ptr;	r00   = tmp;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x280;					// Head of RADIX*vec_cmplx-sized local store #2
	  #if !defined(MULTITHREAD) && defined(USE_SSE2)
						s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
	  #endif
		tmp += 0x280;
		two   = tmp + 0x0;
		// DFT-64 roots all def'd locally in the DFT-64 functions, only need the radix-5 roots here:
		ycc1  = tmp + 0x1;	// radix-5 DFT trig consts
		ycc2  = tmp + 0x2;
		yss1  = tmp + 0x3;
		yss2  = tmp + 0x4;
		yss3  = tmp + 0x5;
		tmp += 0x6;	// sc_ptr += 0x506; add a pad so the tmp += 2 below gives half_arr = 0x508 + (RADIX/vec_dbl-length)
	  #ifdef USE_AVX512
		cy = tmp;		tmp += 0x28;	// RADIX/8 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;	// sc_ptr += 0x(508 + 28) = 0x538; This is where the value of half_arr_offset320 comes from
	  #elif defined(USE_AVX)
		cy = tmp;		tmp += 0x50;	// RADIX/4 vec_dbl slots
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(508 + 50) = 0x560; This is where the value of half_arr_offset320 comes from
		half_arr= tmp + 0x02;	// This table needs 20 vec_dbl in both avx and sse2 mode
	  #else
		cy = tmp;		tmp += 0xa0;	// RADIX/2 vec_dbl slots
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(508 + a0) = 0x5b0; This is where the value of half_arr_offset320 comes from
		half_arr= tmp + 0x02;
	  #endif
		ASSERT((radix320_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r00) + (20 << L2_SZ_VD), "radix320_creals_in_local_store checksum failed!");
		/* These remain fixed: */
		VEC_DBL_INIT(two  , 2.0  );
		VEC_DBL_INIT(ycc1, cc1  );	// radix-5 DFT trig consts
		VEC_DBL_INIT(ycc2, cc2  );
		VEC_DBL_INIT(yss1, s2   );
		VEC_DBL_INIT(yss2, ss1  );
		VEC_DBL_INIT(yss3, ss2  );
		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		sse2_rnd->d0 = base[0]; sse2_rnd->d1 = baseinv[1]; sse2_rnd->d2 = wts_mult[1]; sse2_rnd->d3 = inv_mult[0];
	  #else
		VEC_DBL_INIT(sse2_rnd, crnd);
	  #endif

	// Init-mode calls to these functions which maintain an internal local-alloc static store:
		thr_id = -1;	// Use this special thread id for any macro-required thread-local-data inits...
		SSE2_RADIX_64_DIF(CY_THREADS,thr_id, 0,0,0,0,0,0);
		SSE2_RADIX_64_DIT(CY_THREADS,thr_id, 0,0,0,0,0);
		thr_id = 0;	// ...then revert to 0.

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

		pini = NDIVR/CY_THREADS;
		/*   constant index offsets for array load/stores are here.	*/
		p1 = NDIVR;
		p2 = p1 + NDIVR;		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );
		p3 = p2 + NDIVR;		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );
		p4 = p3 + NDIVR;		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );
		p5 = p4 + NDIVR;		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );
		p6 = p5 + NDIVR;		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );
		p7 = p6 + NDIVR;		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );
		p8 = p7 + NDIVR;		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );
		p9 = p8 + NDIVR;		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );
		pa = p9 + NDIVR;		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );
		pb = pa + NDIVR;		pa += ( (pa >> DAT_BITS) << PAD_BITS );
		pc = pb + NDIVR;		pb += ( (pb >> DAT_BITS) << PAD_BITS );
		pd = pc + NDIVR;		pc += ( (pc >> DAT_BITS) << PAD_BITS );
		pe = pd + NDIVR;		pd += ( (pd >> DAT_BITS) << PAD_BITS );
		pf = pe + NDIVR;		pe += ( (pe >> DAT_BITS) << PAD_BITS );
		NDIVR <<= 4;			pf += ( (pf >> DAT_BITS) << PAD_BITS );
		p10 = NDIVR;
		p20 = p10 + NDIVR;		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p30 = p20 + NDIVR;		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p40 = p30 + NDIVR;		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
		p50 = p40 + NDIVR;		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p60 = p50 + NDIVR;		p50 += ( (p50 >> DAT_BITS) << PAD_BITS );
		p70 = p60 + NDIVR;		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		p80 = p70 + NDIVR;		p70 += ( (p70 >> DAT_BITS) << PAD_BITS );
		p90 = p80 + NDIVR;		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		pa0 = p90 + NDIVR;		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pb0 = pa0 + NDIVR;		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pc0 = pb0 + NDIVR;		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		pd0 = pc0 + NDIVR;		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pe0 = pd0 + NDIVR;		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		pf0 = pe0 + NDIVR;		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		pg0 = pf0 + NDIVR;		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
		ph0 = pg0 + NDIVR;		pg0 += ( (pg0 >> DAT_BITS) << PAD_BITS );
		pi0 = ph0 + NDIVR;		ph0 += ( (ph0 >> DAT_BITS) << PAD_BITS );
		pj0 = pi0 + NDIVR;		pi0 += ( (pi0 >> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			pj0 += ( (pj0 >> DAT_BITS) << PAD_BITS );
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
		poff[0x40+0] = pg0; poff[0x40+1] = pg0+p4; poff[0x40+2] = pg0+p8; poff[0x40+3] = pg0+pc;
		poff[0x44+0] = ph0; poff[0x44+1] = ph0+p4; poff[0x44+2] = ph0+p8; poff[0x44+3] = ph0+pc;
		poff[0x48+0] = pi0; poff[0x48+1] = pi0+p4; poff[0x48+2] = pi0+p8; poff[0x48+3] = pi0+pc;
		poff[0x4c+0] = pj0; poff[0x4c+1] = pj0+p4; poff[0x4c+2] = pj0+p8; poff[0x4c+3] = pj0+pc;

	#ifndef MULTITHREAD
	// Shared:
		// Set array offsets for radix-64 DFT in/outputs:
		// t_offsets w.r.to: t-array, same for all 5 DFTs:
		t_offsets[0x00] = 0x00<<1;	t_offsets[0x10] = 0x10<<1;	t_offsets[0x20] = 0x20<<1;	t_offsets[0x30] = 0x30<<1;
		t_offsets[0x01] = 0x01<<1;	t_offsets[0x11] = 0x11<<1;	t_offsets[0x21] = 0x21<<1;	t_offsets[0x31] = 0x31<<1;
		t_offsets[0x02] = 0x02<<1;	t_offsets[0x12] = 0x12<<1;	t_offsets[0x22] = 0x22<<1;	t_offsets[0x32] = 0x32<<1;
		t_offsets[0x03] = 0x03<<1;	t_offsets[0x13] = 0x13<<1;	t_offsets[0x23] = 0x23<<1;	t_offsets[0x33] = 0x33<<1;
		t_offsets[0x04] = 0x04<<1;	t_offsets[0x14] = 0x14<<1;	t_offsets[0x24] = 0x24<<1;	t_offsets[0x34] = 0x34<<1;
		t_offsets[0x05] = 0x05<<1;	t_offsets[0x15] = 0x15<<1;	t_offsets[0x25] = 0x25<<1;	t_offsets[0x35] = 0x35<<1;
		t_offsets[0x06] = 0x06<<1;	t_offsets[0x16] = 0x16<<1;	t_offsets[0x26] = 0x26<<1;	t_offsets[0x36] = 0x36<<1;
		t_offsets[0x07] = 0x07<<1;	t_offsets[0x17] = 0x17<<1;	t_offsets[0x27] = 0x27<<1;	t_offsets[0x37] = 0x37<<1;
		t_offsets[0x08] = 0x08<<1;	t_offsets[0x18] = 0x18<<1;	t_offsets[0x28] = 0x28<<1;	t_offsets[0x38] = 0x38<<1;
		t_offsets[0x09] = 0x09<<1;	t_offsets[0x19] = 0x19<<1;	t_offsets[0x29] = 0x29<<1;	t_offsets[0x39] = 0x39<<1;
		t_offsets[0x0a] = 0x0a<<1;	t_offsets[0x1a] = 0x1a<<1;	t_offsets[0x2a] = 0x2a<<1;	t_offsets[0x3a] = 0x3a<<1;
		t_offsets[0x0b] = 0x0b<<1;	t_offsets[0x1b] = 0x1b<<1;	t_offsets[0x2b] = 0x2b<<1;	t_offsets[0x3b] = 0x3b<<1;
		t_offsets[0x0c] = 0x0c<<1;	t_offsets[0x1c] = 0x1c<<1;	t_offsets[0x2c] = 0x2c<<1;	t_offsets[0x3c] = 0x3c<<1;
		t_offsets[0x0d] = 0x0d<<1;	t_offsets[0x1d] = 0x1d<<1;	t_offsets[0x2d] = 0x2d<<1;	t_offsets[0x3d] = 0x3d<<1;
		t_offsets[0x0e] = 0x0e<<1;	t_offsets[0x1e] = 0x1e<<1;	t_offsets[0x2e] = 0x2e<<1;	t_offsets[0x3e] = 0x3e<<1;
		t_offsets[0x0f] = 0x0f<<1;	t_offsets[0x1f] = 0x1f<<1;	t_offsets[0x2f] = 0x2f<<1;	t_offsets[0x3f] = 0x3f<<1;

		/*** NB: Delay the plo[]-analog of the phi[] SIMD fixed-addressing to *after* the inits below
		because DIF o_offset and DIT i_offset-array inits need the p-multiplied versions of plo,phi ***/
		plo[0x0] =  0; plo[0x1] = p1; plo[0x2] = p2; plo[0x3] = p3;
		plo[0x4] = p4; plo[0x5] = p5; plo[0x6] = p6; plo[0x7] = p7;
		plo[0x8] = p8; plo[0x9] = p9; plo[0xa] = pa; plo[0xb] = pb;
		plo[0xc] = pc; plo[0xd] = pd; plo[0xe] = pe; plo[0xf] = pf;

	   #ifdef USE_SSE2
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p1*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:
		//phi[0x00]= 0x00<<1; phi[0x01]= 0x10<<1; phi[0x02]= 0x20<<1; phi[0x03]= 0x30<<1;
		//phi[0x04]= 0x40<<1; phi[0x05]= 0x50<<1; phi[0x06]= 0x60<<1; phi[0x07]= 0x70<<1;
		//phi[0x08]= 0x80<<1; phi[0x09]= 0x90<<1; phi[0x0a]= 0xa0<<1; phi[0x0b]= 0xb0<<1;
		//phi[0x0c]= 0xc0<<1; phi[0x0d]= 0xd0<<1; phi[0x0e]= 0xe0<<1; phi[0x0f]= 0xf0<<1;
		//phi[0x10]=0x100<<1; phi[0x11]=0x110<<1; phi[0x12]=0x120<<1; phi[0x13]=0x130<<1;
	   #else
		phi[0x00]=   0; phi[0x01]= p10; phi[0x02]= p20; phi[0x03]= p30;
		phi[0x04]= p40; phi[0x05]= p50; phi[0x06]= p60; phi[0x07]= p70;
		phi[0x08]= p80; phi[0x09]= p90; phi[0x0a]= pa0; phi[0x0b]= pb0;
		phi[0x0c]= pc0; phi[0x0d]= pd0; phi[0x0e]= pe0; phi[0x0f]= pf0;
		phi[0x10]= pg0; phi[0x11]= ph0; phi[0x12]= pi0; phi[0x13]= pj0;
	   #endif	// sse2?

	/*** DIF Output permutation: */
	// Operm 1:
		i64 = dif_16_operms[0];	// [A] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;	// p-offset indices encoded in little-endian hex-char fashion:
			dif_o_offsets[     l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dif_16_operms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x10+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x20+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dif_16_operms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x30+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
	// Operm 2:
		i64 = dif_16_operms[4];	// [E] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x40+l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x50+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dif_16_operms[5];	// [F] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x60+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dif_16_operms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x70+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
	// Operm 3:
		i64 = dif_16_operms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x80+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dif_16_operms[5];	// [F] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x90+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dif_16_operms[4];	// [E] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0xa0+l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0xb0+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
	// Operm 4:
		i64 = dif_16_operms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0xc0+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dif_16_operms[4];	// [E] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0xd0+l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dif_16_operms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0xe0+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dif_16_operms[5];	// [F] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0xf0+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
	// Operm 5:
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x100+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dif_16_operms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x110+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dif_16_operms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x120+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dif_16_operms[4];	// [E] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x130+l] = plo[(i64 >> nshift)&0xf];
		}

	/*** DIT indexing stuff: ***/
		// Init storage for 4 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 4*9:
		l = 0;
	  #ifdef USE_SSE2
		dit_p20_cperms[l++] = 0x000<<1; dit_p20_cperms[l++] = 0x100<<1; dit_p20_cperms[l++] = 0x0c0<<1; dit_p20_cperms[l++] = 0x080<<1; dit_p20_cperms[l++] = 0x040<<1; dit_p20_cperms[l++] = 0x000<<1; dit_p20_cperms[l++] = 0x100<<1; dit_p20_cperms[l++] = 0x0c0<<1; dit_p20_cperms[l++] = 0x080<<1;
		dit_p20_cperms[l++] = 0x010<<1; dit_p20_cperms[l++] = 0x110<<1; dit_p20_cperms[l++] = 0x0d0<<1; dit_p20_cperms[l++] = 0x090<<1; dit_p20_cperms[l++] = 0x050<<1; dit_p20_cperms[l++] = 0x010<<1; dit_p20_cperms[l++] = 0x110<<1; dit_p20_cperms[l++] = 0x0d0<<1; dit_p20_cperms[l++] = 0x090<<1;
		dit_p20_cperms[l++] = 0x020<<1; dit_p20_cperms[l++] = 0x120<<1; dit_p20_cperms[l++] = 0x0e0<<1; dit_p20_cperms[l++] = 0x0a0<<1; dit_p20_cperms[l++] = 0x060<<1; dit_p20_cperms[l++] = 0x020<<1; dit_p20_cperms[l++] = 0x120<<1; dit_p20_cperms[l++] = 0x0e0<<1; dit_p20_cperms[l++] = 0x0a0<<1;
		dit_p20_cperms[l++] = 0x030<<1; dit_p20_cperms[l++] = 0x130<<1; dit_p20_cperms[l++] = 0x0f0<<1; dit_p20_cperms[l++] = 0x0b0<<1; dit_p20_cperms[l++] = 0x070<<1; dit_p20_cperms[l++] = 0x030<<1; dit_p20_cperms[l++] = 0x130<<1; dit_p20_cperms[l++] = 0x0f0<<1; dit_p20_cperms[l++] = 0x0b0<<1;
	  #else
		dit_p20_cperms[l++] =   0; dit_p20_cperms[l++] = pg0; dit_p20_cperms[l++] = pc0; dit_p20_cperms[l++] = p80; dit_p20_cperms[l++] = p40; dit_p20_cperms[l++] =   0; dit_p20_cperms[l++] = pg0; dit_p20_cperms[l++] = pc0; dit_p20_cperms[l++] = p80;
		dit_p20_cperms[l++] = p10; dit_p20_cperms[l++] = ph0; dit_p20_cperms[l++] = pd0; dit_p20_cperms[l++] = p90; dit_p20_cperms[l++] = p50; dit_p20_cperms[l++] = p10; dit_p20_cperms[l++] = ph0; dit_p20_cperms[l++] = pd0; dit_p20_cperms[l++] = p90;
		dit_p20_cperms[l++] = p20; dit_p20_cperms[l++] = pi0; dit_p20_cperms[l++] = pe0; dit_p20_cperms[l++] = pa0; dit_p20_cperms[l++] = p60; dit_p20_cperms[l++] = p20; dit_p20_cperms[l++] = pi0; dit_p20_cperms[l++] = pe0; dit_p20_cperms[l++] = pa0;
		dit_p20_cperms[l++] = p30; dit_p20_cperms[l++] = pj0; dit_p20_cperms[l++] = pf0; dit_p20_cperms[l++] = pb0; dit_p20_cperms[l++] = p70; dit_p20_cperms[l++] = p30; dit_p20_cperms[l++] = pj0; dit_p20_cperms[l++] = pf0; dit_p20_cperms[l++] = pb0;
	  #endif
	// Iperm 1:
		i64 = dit_16_iperms[0];	// [A] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;	// p-offset indices encoded in little-endian hex-char fashion:
			dit_i_offsets[     l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dit_16_iperms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0x10+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dit_16_iperms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0x20+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dit_16_iperms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0x30+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
	// Iperm 2:
		i64 = dit_16_iperms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0x40+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dit_16_iperms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0x50+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dit_16_iperms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0x60+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dit_16_iperms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0x70+l] = plo[(i64 >> nshift)&0xf];
		}
	// Iperm 3:
		i64 = dit_16_iperms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0x80+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dit_16_iperms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0x90+l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dit_16_iperms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0xa0+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dit_16_iperms[4];	// [E] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0xb0+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
	// Iperm 4:
		i64 = dit_16_iperms[5];	// [F] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0xc0+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dit_16_iperms[6];	// [G] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0xd0+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dit_16_iperms[5];	// [F] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0xe0+l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dit_16_iperms[6];	// [G] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0xf0+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
	// Iperm 5:
		i64 = dit_16_iperms[7];	// [H] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0x100+l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dit_16_iperms[8];	// [i] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0x110+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dit_16_iperms[8];	// [I] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0x120+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dit_16_iperms[8];	// [I] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0x130+l] = plo[(i64 >> nshift)&0xf] + p20;
		}

	   #ifdef USE_SSE2
		/*** Delayed part of SIMD fixed-offsets init: ***/
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p1*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:
		plo[0x0] = 0x0<<1; plo[0x1] = 0x1<<1; plo[0x2] = 0x2<<1; plo[0x3] = 0x3<<1;
		plo[0x4] = 0x4<<1; plo[0x5] = 0x5<<1; plo[0x6] = 0x6<<1; plo[0x7] = 0x7<<1;
		plo[0x8] = 0x8<<1; plo[0x9] = 0x9<<1; plo[0xa] = 0xa<<1; plo[0xb] = 0xb<<1;
		plo[0xc] = 0xc<<1; plo[0xd] = 0xd<<1; plo[0xe] = 0xe<<1; plo[0xf] = 0xf<<1;
	   #endif	// sse2?

	#endif

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
			tdat[ithread].r00 = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (vec_dbl *)((intptr_t)tdat[ithread].r00 + ((intptr_t)half_arr - (intptr_t)r00));
		#else	// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
			tdat[ithread].r00      = (double *)base;
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


/*...The radix-320 final DIT pass is here.	*/

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
		ASSERT(tdat[ithread].r00 == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
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
		#include "radix320_main_carry_loop.h"

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
		ASSERT(0x0 == cy320_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
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
				jt = j1 + poff[l];	// poff[] = p0,4,8,...
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

void radix320_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-320 complex DIF FFT pass on the data in the length-N real vector A.
*/
	int j,j1,j2,jp,jt,l,ilo,ihi,hi_neg,k0,k1,k2,k3,k4,nshift;
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,pg0,ph0,pi0,pj0, first_entry=TRUE;
	static int plo[16],phi[20];
	const double cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
				cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
				s2  =  0.95105651629515357211,	/*  sin(u) */
				ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
				ss2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	static int t_offsets[64], o_offsets[RADIX];
	// Consts containing [in little-endian hex-char fashion] the 6 [A-F] 16-subperms defined in the output-ordering commentary below:
	const uint64 dif_16_operms[6] = {
		0x01326745cdfe98abull,	// [A]
		0x98abfedc32104576ull,	// [B]
		0x32104576fedcab89ull,	// [C]
		0x6745321098abfedcull,	// [D]
		0xcdfe98ab67453210ull,	// [E]
		0xfedcab8945761023ull	// [F]
	};
	uint64 i64;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex t[RADIX], *tptr;
	double rt,it;

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
		p1 = NDIVR;
		p2 = p1 + NDIVR;		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );
		p3 = p2 + NDIVR;		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );
		p4 = p3 + NDIVR;		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );
		p5 = p4 + NDIVR;		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );
		p6 = p5 + NDIVR;		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );
		p7 = p6 + NDIVR;		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );
		p8 = p7 + NDIVR;		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );
		p9 = p8 + NDIVR;		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );
		pa = p9 + NDIVR;		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );
		pb = pa + NDIVR;		pa += ( (pa >> DAT_BITS) << PAD_BITS );
		pc = pb + NDIVR;		pb += ( (pb >> DAT_BITS) << PAD_BITS );
		pd = pc + NDIVR;		pc += ( (pc >> DAT_BITS) << PAD_BITS );
		pe = pd + NDIVR;		pd += ( (pd >> DAT_BITS) << PAD_BITS );
		pf = pe + NDIVR;		pe += ( (pe >> DAT_BITS) << PAD_BITS );
		NDIVR <<= 4;			pf += ( (pf >> DAT_BITS) << PAD_BITS );
		p10 = NDIVR;
		p20 = p10 + NDIVR;		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p30 = p20 + NDIVR;		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p40 = p30 + NDIVR;		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
		p50 = p40 + NDIVR;		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p60 = p50 + NDIVR;		p50 += ( (p50 >> DAT_BITS) << PAD_BITS );
		p70 = p60 + NDIVR;		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		p80 = p70 + NDIVR;		p70 += ( (p70 >> DAT_BITS) << PAD_BITS );
		p90 = p80 + NDIVR;		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		pa0 = p90 + NDIVR;		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pb0 = pa0 + NDIVR;		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pc0 = pb0 + NDIVR;		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		pd0 = pc0 + NDIVR;		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pe0 = pd0 + NDIVR;		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		pf0 = pe0 + NDIVR;		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		pg0 = pf0 + NDIVR;		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
		ph0 = pg0 + NDIVR;		pg0 += ( (pg0 >> DAT_BITS) << PAD_BITS );
		pi0 = ph0 + NDIVR;		ph0 += ( (ph0 >> DAT_BITS) << PAD_BITS );
		pj0 = pi0 + NDIVR;		pi0 += ( (pi0 >> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			pj0 += ( (pj0 >> DAT_BITS) << PAD_BITS );

		plo[0x0] =  0; plo[0x1] = p1; plo[0x2] = p2; plo[0x3] = p3;
		plo[0x4] = p4; plo[0x5] = p5; plo[0x6] = p6; plo[0x7] = p7;
		plo[0x8] = p8; plo[0x9] = p9; plo[0xa] = pa; plo[0xb] = pb;
		plo[0xc] = pc; plo[0xd] = pd; plo[0xe] = pe; plo[0xf] = pf;

		phi[0x00]=   0; phi[0x01]= p10; phi[0x02]= p20; phi[0x03]= p30;
		phi[0x04]= p40; phi[0x05]= p50; phi[0x06]= p60; phi[0x07]= p70;
		phi[0x08]= p80; phi[0x09]= p90; phi[0x0a]= pa0; phi[0x0b]= pb0;
		phi[0x0c]= pc0; phi[0x0d]= pd0; phi[0x0e]= pe0; phi[0x0f]= pf0;
		phi[0x10]= pg0; phi[0x11]= ph0; phi[0x12]= pi0; phi[0x13]= pj0;

		t_offsets[0x00] = 0x00<<1;	t_offsets[0x10] = 0x10<<1;	t_offsets[0x20] = 0x20<<1;	t_offsets[0x30] = 0x30<<1;
		t_offsets[0x01] = 0x01<<1;	t_offsets[0x11] = 0x11<<1;	t_offsets[0x21] = 0x21<<1;	t_offsets[0x31] = 0x31<<1;
		t_offsets[0x02] = 0x02<<1;	t_offsets[0x12] = 0x12<<1;	t_offsets[0x22] = 0x22<<1;	t_offsets[0x32] = 0x32<<1;
		t_offsets[0x03] = 0x03<<1;	t_offsets[0x13] = 0x13<<1;	t_offsets[0x23] = 0x23<<1;	t_offsets[0x33] = 0x33<<1;
		t_offsets[0x04] = 0x04<<1;	t_offsets[0x14] = 0x14<<1;	t_offsets[0x24] = 0x24<<1;	t_offsets[0x34] = 0x34<<1;
		t_offsets[0x05] = 0x05<<1;	t_offsets[0x15] = 0x15<<1;	t_offsets[0x25] = 0x25<<1;	t_offsets[0x35] = 0x35<<1;
		t_offsets[0x06] = 0x06<<1;	t_offsets[0x16] = 0x16<<1;	t_offsets[0x26] = 0x26<<1;	t_offsets[0x36] = 0x36<<1;
		t_offsets[0x07] = 0x07<<1;	t_offsets[0x17] = 0x17<<1;	t_offsets[0x27] = 0x27<<1;	t_offsets[0x37] = 0x37<<1;
		t_offsets[0x08] = 0x08<<1;	t_offsets[0x18] = 0x18<<1;	t_offsets[0x28] = 0x28<<1;	t_offsets[0x38] = 0x38<<1;
		t_offsets[0x09] = 0x09<<1;	t_offsets[0x19] = 0x19<<1;	t_offsets[0x29] = 0x29<<1;	t_offsets[0x39] = 0x39<<1;
		t_offsets[0x0a] = 0x0a<<1;	t_offsets[0x1a] = 0x1a<<1;	t_offsets[0x2a] = 0x2a<<1;	t_offsets[0x3a] = 0x3a<<1;
		t_offsets[0x0b] = 0x0b<<1;	t_offsets[0x1b] = 0x1b<<1;	t_offsets[0x2b] = 0x2b<<1;	t_offsets[0x3b] = 0x3b<<1;
		t_offsets[0x0c] = 0x0c<<1;	t_offsets[0x1c] = 0x1c<<1;	t_offsets[0x2c] = 0x2c<<1;	t_offsets[0x3c] = 0x3c<<1;
		t_offsets[0x0d] = 0x0d<<1;	t_offsets[0x1d] = 0x1d<<1;	t_offsets[0x2d] = 0x2d<<1;	t_offsets[0x3d] = 0x3d<<1;
		t_offsets[0x0e] = 0x0e<<1;	t_offsets[0x1e] = 0x1e<<1;	t_offsets[0x2e] = 0x2e<<1;	t_offsets[0x3e] = 0x3e<<1;
		t_offsets[0x0f] = 0x0f<<1;	t_offsets[0x1f] = 0x1f<<1;	t_offsets[0x2f] = 0x2f<<1;	t_offsets[0x3f] = 0x3f<<1;
	/*** Output permutation: */
	#if 0
	/* Use these only for initial run via test_fft_radices, used to compute the needed output-perm: */
		o_offsets[0x00] =  0	;	o_offsets[0x10] =  0+p10;	o_offsets[0x20] =  0+p20;	o_offsets[0x30] =  0+p30;
		o_offsets[0x01] = p1	;	o_offsets[0x11] = p1+p10;	o_offsets[0x21] = p1+p20;	o_offsets[0x31] = p1+p30;
		o_offsets[0x02] = p2	;	o_offsets[0x12] = p2+p10;	o_offsets[0x22] = p2+p20;	o_offsets[0x32] = p2+p30;
		o_offsets[0x03] = p3	;	o_offsets[0x13] = p3+p10;	o_offsets[0x23] = p3+p20;	o_offsets[0x33] = p3+p30;
		o_offsets[0x04] = p4	;	o_offsets[0x14] = p4+p10;	o_offsets[0x24] = p4+p20;	o_offsets[0x34] = p4+p30;
		o_offsets[0x05] = p5	;	o_offsets[0x15] = p5+p10;	o_offsets[0x25] = p5+p20;	o_offsets[0x35] = p5+p30;
		o_offsets[0x06] = p6	;	o_offsets[0x16] = p6+p10;	o_offsets[0x26] = p6+p20;	o_offsets[0x36] = p6+p30;
		o_offsets[0x07] = p7	;	o_offsets[0x17] = p7+p10;	o_offsets[0x27] = p7+p20;	o_offsets[0x37] = p7+p30;
		o_offsets[0x08] = p8	;	o_offsets[0x18] = p8+p10;	o_offsets[0x28] = p8+p20;	o_offsets[0x38] = p8+p30;
		o_offsets[0x09] = p9	;	o_offsets[0x19] = p9+p10;	o_offsets[0x29] = p9+p20;	o_offsets[0x39] = p9+p30;
		o_offsets[0x0a] = pa	;	o_offsets[0x1a] = pa+p10;	o_offsets[0x2a] = pa+p20;	o_offsets[0x3a] = pa+p30;
		o_offsets[0x0b] = pb	;	o_offsets[0x1b] = pb+p10;	o_offsets[0x2b] = pb+p20;	o_offsets[0x3b] = pb+p30;
		o_offsets[0x0c] = pc	;	o_offsets[0x1c] = pc+p10;	o_offsets[0x2c] = pc+p20;	o_offsets[0x3c] = pc+p30;
		o_offsets[0x0d] = pd	;	o_offsets[0x1d] = pd+p10;	o_offsets[0x2d] = pd+p20;	o_offsets[0x3d] = pd+p30;
		o_offsets[0x0e] = pe	;	o_offsets[0x1e] = pe+p10;	o_offsets[0x2e] = pe+p20;	o_offsets[0x3e] = pe+p30;
		o_offsets[0x0f] = pf	;	o_offsets[0x1f] = pf+p10;	o_offsets[0x2f] = pf+p20;	o_offsets[0x3f] = pf+p30;
		for(l = 0; l < 64; l++) {
			o_offsets[0x040+l] = o_offsets[l];
			o_offsets[0x080+l] = o_offsets[l];
			o_offsets[0x0c0+l] = o_offsets[l];
			o_offsets[0x100+l] = o_offsets[l];
		}
	//
	#else
	// Extract the 6 [A-F] 16-subperms defined in the output-ordering commentary below:
	/* Operm 1:
		 0,1,3,2,6,7,4,5,c,d,f,e,9,8,a,b + p00  =  [A] + p00
		,9,8,a,b,f,e,d,c,3,2,1,0,4,5,7,6 + p10  =  [B] + p10
		,3,2,1,0,4,5,7,6,f,e,d,c,a,b,8,9 + p30  =  [C] + p30
		,6,7,4,5,3,2,1,0,9,8,a,b,f,e,d,c + p20  =  [D] + p20
	*/
		i64 = dif_16_operms[0];	// [A] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;	// p-offset indices encoded in little-endian hex-char fashion:
			o_offsets[     l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dif_16_operms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x10+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x20+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dif_16_operms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x30+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
	/* Operm 2:
		,c,d,f,e,9,8,a,b,6,7,4,5,3,2,1,0 + pg0  =  [E] + p00
		,3,2,1,0,4,5,7,6,f,e,d,c,a,b,8,9 + ph0  =  [C] + p10
		,f,e,d,c,a,b,8,9,4,5,7,6,1,0,2,3 + pj0  =  [F] + p30 (mod p40)
		,9,8,a,b,f,e,d,c,3,2,1,0,4,5,7,6 + pi0  =  [B] + p20
	*/
		i64 = dif_16_operms[4];	// [E] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x40+l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x50+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dif_16_operms[5];	// [F] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x60+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dif_16_operms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x70+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
	/* Operm 3:
		,6,7,4,5,3,2,1,0,9,8,a,b,f,e,d,c + pe0  =  [D] + p20
		,f,e,d,c,a,b,8,9,4,5,7,6,1,0,2,3 + pf0  =  [F] + p30
		,c,d,f,e,9,8,a,b,6,7,4,5,3,2,1,0 + pc0  =  [E] + p00 (mod p40)
		,3,2,1,0,4,5,7,6,f,e,d,c,a,b,8,9 + pd0  =  [C] + p10
	*/
		i64 = dif_16_operms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x80+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dif_16_operms[5];	// [F] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x90+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dif_16_operms[4];	// [E] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0xa0+l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0xb0+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
	/* Operm 4:
		,9,8,a,b,f,e,d,c,3,2,1,0,4,5,7,6 + p90  =  [B] + p10
		,c,d,f,e,9,8,a,b,6,7,4,5,3,2,1,0 + p80  =  [E] + p00
		,6,7,4,5,3,2,1,0,9,8,a,b,f,e,d,c + pa0  =  [D] + p20 (mod p40)
		,f,e,d,c,a,b,8,9,4,5,7,6,1,0,2,3 + pb0  =  [F] + p30
	*/
		i64 = dif_16_operms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0xc0+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dif_16_operms[4];	// [E] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0xd0+l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dif_16_operms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0xe0+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dif_16_operms[5];	// [F] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0xf0+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
	/* Operm 5:
		,3,2,1,0,4,5,7,6,f,e,d,c,a,b,8,9 + p70  =  [C] + p30
		,6,7,4,5,3,2,1,0,9,8,a,b,f,e,d,c + p60  =  [D] + p20
		,9,8,a,b,f,e,d,c,3,2,1,0,4,5,7,6 + p50  =  [B] + p10 (mod p40)
		,c,d,f,e,9,8,a,b,6,7,4,5,3,2,1,0 + p40  =  [E] + p00,
	*/
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x100+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dif_16_operms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x110+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dif_16_operms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x120+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dif_16_operms[4];	// [E] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x130+l] = plo[(i64 >> nshift)&0xf];
		}
	#endif
	}

/*...The radix-320 pass is here.	*/

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
	Twiddleless version arranges 64 sets of radix-5 DFT inputs as follows: 0 in upper left corner,
	decrement 64 horizontally and 5 vertically (mod 320):

	bc-code to generate indices:
	obase=16;
	n=320;
	a=0;b=n-64;c=b-64;d=c-64;e=d-64;
	i = 64;
	while(i) { print a+(a<0)*n,",",b+(b<0)*n,",",c+(c<0)*n,",",d+(d<0)*n,",",e+(e<0)*n,"\n"; a-=5;b-=5;c-=5;d-=5;e-=5; i-=1; }

	DIF/DIT input-scramble array = [
		000,100,0c0,080,040   000,100,0c0,080,040 + p0		0f0,0b0,070,030,130   0f0,0b0,070,030,130 + p0		0a0,060,020,120,0e0   0a0,060,020,120,0e0 + p0		050,010,110,0d0,090   050,010,110,0d0,090 + p0
		13b,0fb,0bb,07b,03b   130,0f0,0b0,070,030 + pb		0eb,0ab,06b,02b,12b   0e0,0a0,060,020,120 + pb		09b,05b,01b,11b,0db   090,050,010,110,0d0 + pb		04b,00b,10b,0cb,08b   040,000,100,0c0,080 + pb
		136,0f6,0b6,076,036   130,0f0,0b0,070,030 + p6		0e6,0a6,066,026,126   0e0,0a0,060,020,120 + p6		096,056,016,116,0d6   090,050,010,110,0d0 + p6		046,006,106,0c6,086   040,000,100,0c0,080 + p6
		131,0f1,0b1,071,031   130,0f0,0b0,070,030 + p1		0e1,0a1,061,021,121   0e0,0a0,060,020,120 + p1		091,051,011,111,0d1   090,050,010,110,0d0 + p1		041,001,101,0c1,081   040,000,100,0c0,080 + p1
		12c,0ec,0ac,06c,02c   120,0e0,0a0,060,020 + pc		0dc,09c,05c,01c,11c   0d0,090,050,010,110 + pc		08c,04c,00c,10c,0cc   080,040,000,100,0c0 + pc		03c,13c,0fc,0bc,07c   030,130,0f0,0b0,070 + pc
		127,0e7,0a7,067,027   120,0e0,0a0,060,020 + p7		0d7,097,057,017,117   0d0,090,050,010,110 + p7		087,047,007,107,0c7   080,040,000,100,0c0 + p7		037,137,0f7,0b7,077   030,130,0f0,0b0,070 + p7
		122,0e2,0a2,062,022   120,0e0,0a0,060,020 + p2		0d2,092,052,012,112   0d0,090,050,010,110 + p2		082,042,002,102,0c2   080,040,000,100,0c0 + p2		032,132,0f2,0b2,072   030,130,0f0,0b0,070 + p2
		11d,0dd,09d,05d,01d   110,0d0,090,050,010 + pd		0cd,08d,04d,00d,10d   0c0,080,040,000,100 + pd		07d,03d,13d,0fd,0bd   070,030,130,0f0,0b0 + pd		02d,12d,0ed,0ad,06d   020,120,0e0,0a0,060 + pd
		118,0d8,098,058,018 = 110,0d0,090,050,010 + p8		0c8,088,048,008,108 = 0c0,080,040,000,100 + p8		078,038,138,0f8,0b8 = 070,030,130,0f0,0b0 + p8		028,128,0e8,0a8,068 = 020,120,0e0,0a0,060 + p8
		113,0d3,093,053,013   110,0d0,090,050,010 + p3		0c3,083,043,003,103   0c0,080,040,000,100 + p3		073,033,133,0f3,0b3   070,030,130,0f0,0b0 + p3		023,123,0e3,0a3,063   020,120,0e0,0a0,060 + p3
		10e,0ce,08e,04e,00e   100,0c0,080,040,000 + pe		0be,07e,03e,13e,0fe   0b0,070,030,130,0f0 + pe		06e,02e,12e,0ee,0ae   060,020,120,0e0,0a0 + pe		01e,11e,0de,09e,05e   010,110,0d0,090,050 + pe
		109,0c9,089,049,009   100,0c0,080,040,000 + p9		0b9,079,039,139,0f9   0b0,070,030,130,0f0 + p9		069,029,129,0e9,0a9   060,020,120,0e0,0a0 + p9		019,119,0d9,099,059   010,110,0d0,090,050 + p9
		104,0c4,084,044,004   100,0c0,080,040,000 + p4		0b4,074,034,134,0f4   0b0,070,030,130,0f0 + p4		064,024,124,0e4,0a4   060,020,120,0e0,0a0 + p4		014,114,0d4,094,054   010,110,0d0,090,050 + p4
		0ff,0bf,07f,03f,13f   0f0,0b0,070,030,130 + pf		0af,06f,02f,12f,0ef   0a0,060,020,120,0e0 + pf		05f,01f,11f,0df,09f   050,010,110,0d0,090 + pf		00f,10f,0cf,08f,04f   000,100,0c0,080,040 + pf
		0fa,0ba,07a,03a,13a   0f0,0b0,070,030,130 + pa		0aa,06a,02a,12a,0ea   0a0,060,020,120,0e0 + pa		05a,01a,11a,0da,09a   050,010,110,0d0,090 + pa		00a,10a,0ca,08a,04a   000,100,0c0,080,040 + pa
		0f5,0b5,075,035,135   0f0,0b0,070,030,130 + p5		0a5,065,025,125,0e5   0a0,060,020,120,0e0 + p5		055,015,115,0d5,095   050,010,110,0d0,090 + p5		005,105,0c5,085,045   000,100,0c0,080,040 + p5
					[cont. in col 2]									[cont. in col 3]									[cont. in col 4]

	If we prestore the length-5 perms appearing left of the + above we get 12 such perms added to the p0-f offsets,
	which are really just [0-4]-element leftward-cshifts of the 4 base-perms

	[0] = 000,100,0c0,080,040;	[1] = 130,0f0,0b0,070,030;	[2] = 120,0e0,0a0,060,020;	[3] = 110,0d0,090,050,010.

	Thus in terms of these base-perms and the lcshift count:

		[0]<<0 + p0		[1]<<1 + p0		[2]<<2 + p0		[1]<<3 + p0
		[1]<<0 + pb		[2]<<1 + pb		[3]<<2 + pb		[2]<<4 + pb
		[1]<<0 + p6		[2]<<1 + p6		[3]<<2 + p6		[2]<<4 + p6
		[1]<<0 + p1		[2]<<1 + p1		[3]<<2 + p1		[2]<<4 + p1
		[2]<<0 + pc		[3]<<1 + pc		[3]<<3 + pc		[2]<<4 + pc
		[2]<<0 + p7		[3]<<1 + p7		[3]<<3 + p7		[2]<<4 + p7
		[2]<<0 + p2		[3]<<1 + p2		[0]<<3 + p2		[1]<<4 + p2
		[3]<<0 + pd		[0]<<2 + pd		[0]<<3 + pd		[1]<<4 + pd
		[3]<<0 + p8		[0]<<2 + p8		[0]<<3 + p8		[1]<<4 + p8
		[3]<<0 + p3		[0]<<2 + p3		[0]<<3 + p3		[1]<<4 + p3
		[0]<<1 + pe		[1]<<2 + pe		[0]<<3 + pe		[1]<<4 + pe
		[0]<<1 + p9		[1]<<2 + p9		[1]<<3 + p9		[0]<<4 + p9
		[0]<<1 + p4		[1]<<2 + p4		[1]<<3 + p4		[0]<<4 + p4
		[1]<<1 + pf		[2]<<2 + pf		[1]<<3 + pf		[0]<<0 + pf
		[1]<<1 + pa		[2]<<2 + pa		[1]<<3 + pa		[0]<<0 + pa
		[1]<<1 + p5		[2]<<2 + p5		[1]<<3 + p5		[0]<<0 + p5
			[cont.2]		[cont.3]		[cont.4]
	Shift counts by row [0-63]:
		 Row		Shift
		 0- 9	0
		10-22	1
		23-35	2
		36-48	3
		49-60	4	[only 12 shift-4s]
		61-63	0
	Thus specification of each complete 5-perm needs:

		2 bits for the base-perm [0-3];
		3 bits for the base-perm lcshift count [0-2];	<*** trouble - now need 9 bits per perm, no longer fits in a byte ***
		4 bits for the p0-f addend in terms of plo array index.

	Thus 64 x 8 bits for the complete set of 64 3-perms.
	Or ... just do some fancy-indexing footwork to directly implement the horiz/vert index decrementing, as below:
	*/
	//...gather the needed data (320 64-bit complex) and do 64 radix-5 transforms - We want unit-strides in the radix64-DFT macros, so use large output strides here:
		tptr = t;
		ilo = ihi = 0;	// ints tracking the low and high nibbles of the p* offset at head  of each trio
		k0 =   0;	k1 = pg0;	k2 = pc0;	k3 = p80;	k4 = p40;
		for(l = 0; l < 64; l++) {
			jp = plo[ilo];	// Same for each elt of the current-offset trio
//printf("Set %2u: k0-4 = %2X,%2X,%2X,%2X,%2X, ilo = %2X\n",l, k0>>1,k1>>1,k2>>1,k3>>1,k4>>1, ilo);
			jt = j1 + jp; jp += j2;
//printf("Set %2u: a[re,im]0-4 = %1.0f,%1.0f[%u,%u]\n",l, a[jt+k0],a[jp+k0],ref[jt+k0/2],ref[jp+k0/2]);//, a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4]);
			RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],
				tptr->re,tptr->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x100)->re,(tptr+0x100)->im,
				rt,it);
			tptr++;
			ilo -= 5;
			// If low nibble underflows need to add 16 to it, decrement high nibble ihi, and recompute k1-4.
			// Use branch here because dependent-op count > 10:
			if(ilo < 0) {
				ilo += 16;
				ihi--;			hi_neg = ihi< 0;	ihi+= ((-hi_neg) & 20);	// 1st high-part offset
				k1 = ihi- 4;	hi_neg = k1 < 0;	k1 += ((-hi_neg) & 20);	// high part of "decrement 64 horizontally" to get 2nd offset
				k2 = k1 - 4;	hi_neg = k2 < 0;	k2 += ((-hi_neg) & 20);	// high part of "decrement 64 horizontally" to get 3rd offset
				k3 = k2 - 4;	hi_neg = k3 < 0;	k3 += ((-hi_neg) & 20);	// high part of "decrement 64 horizontally" to get 4th offset
				k4 = k3 - 4;	hi_neg = k4 < 0;	k4 += ((-hi_neg) & 20);	// high part of "decrement 64 horizontally" to get 5th offset
				// Use distinct ihi and k0 = phi[ihi] for 1st high-part offset because basic-index ihi must survive from one loop pass to next:
				k0 = phi[ihi];	k1 = phi[k1];	k2 = phi[k2];	k3 = phi[k3];	k4 = phi[k4];
			}
		}
//exit(0);
	/*...and now do 5 radix-64 transforms.
	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF]
	to properly permute the outputs to the required ordering. Required output index permutation =
	[
		   0,  1,  3,  2,  6,  7,  4,  5,  c,  d,  f,  e,  9,  8,  a,  b
		, 19, 18, 1a, 1b, 1f, 1e, 1d, 1c, 13, 12, 11, 10, 14, 15, 17, 16
		, 33, 32, 31, 30, 34, 35, 37, 36, 3f, 3e, 3d, 3c, 3a, 3b, 38, 39
		, 26, 27, 24, 25, 23, 22, 21, 20, 29, 28, 2a, 2b, 2f, 2e, 2d, 2c
		,10c,10d,10f,10e,109,108,10a,10b,106,107,104,105,103,102,101,100
		,113,112,111,110,114,115,117,116,11f,11e,11d,11c,11a,11b,118,119
		,13f,13e,13d,13c,13a,13b,138,139,134,135,137,136,131,130,132,133
		,129,128,12a,12b,12f,12e,12d,12c,123,122,121,120,124,125,127,126
		, e6, e7, e4, e5, e3, e2, e1, e0, e9, e8, ea, eb, ef, ee, ed, ec
		, ff, fe, fd, fc, fa, fb, f8, f9, f4, f5, f7, f6, f1, f0, f2, f3
		, cc, cd, cf, ce, c9, c8, ca, cb, c6, c7, c4, c5, c3, c2, c1, c0
		, d3, d2, d1, d0, d4, d5, d7, d6, df, de, dd, dc, da, db, d8, d9
		, 99, 98, 9a, 9b, 9f, 9e, 9d, 9c, 93, 92, 91, 90, 94, 95, 97, 96
		, 8c, 8d, 8f, 8e, 89, 88, 8a, 8b, 86, 87, 84, 85, 83, 82, 81, 80
		, a6, a7, a4, a5, a3, a2, a1, a0, a9, a8, aa, ab, af, ae, ad, ac
		, bf, be, bd, bc, ba, bb, b8, b9, b4, b5, b7, b6, b1, b0, b2, b3
		, 73, 72, 71, 70, 74, 75, 77, 76, 7f, 7e, 7d, 7c, 7a, 7b, 78, 79
		, 66, 67, 64, 65, 63, 62, 61, 60, 69, 68, 6a, 6b, 6f, 6e, 6d, 6c
		, 59, 58, 5a, 5b, 5f, 5e, 5d, 5c, 53, 52, 51, 50, 54, 55, 57, 56
		, 4c, 4d, 4f, 4e, 49, 48, 4a, 4b, 46, 47, 44, 45, 43, 42, 41, 40
	]
	Which in terms of our 5 radix-64 DFTs is
		 0,1,3,2,6,7,4,5,c,d,f,e,9,8,a,b + p00  =  [A] + p00
		,9,8,a,b,f,e,d,c,3,2,1,0,4,5,7,6 + p10  =  [B] + p10
		,3,2,1,0,4,5,7,6,f,e,d,c,a,b,8,9 + p30  =  [C] + p30
		,6,7,4,5,3,2,1,0,9,8,a,b,f,e,d,c + p20  =  [D] + p20

		,c,d,f,e,9,8,a,b,6,7,4,5,3,2,1,0 + pg0  =  [E] + pg0
		,3,2,1,0,4,5,7,6,f,e,d,c,a,b,8,9 + ph0  =  [C] + ph0
		,f,e,d,c,a,b,8,9,4,5,7,6,1,0,2,3 + pj0  =  [F] + pj0
		,9,8,a,b,f,e,d,c,3,2,1,0,4,5,7,6 + pi0  =  [B] + pi0

		,6,7,4,5,3,2,1,0,9,8,a,b,f,e,d,c + pe0  =  [D] + pe0
		,f,e,d,c,a,b,8,9,4,5,7,6,1,0,2,3 + pf0  =  [F] + pf0
		,c,d,f,e,9,8,a,b,6,7,4,5,3,2,1,0 + pc0  =  [E] + pc0
		,3,2,1,0,4,5,7,6,f,e,d,c,a,b,8,9 + pd0  =  [C] + pd0

		,9,8,a,b,f,e,d,c,3,2,1,0,4,5,7,6 + p90  =  [B] + p90
		,c,d,f,e,9,8,a,b,6,7,4,5,3,2,1,0 + p80  =  [E] + p80
		,6,7,4,5,3,2,1,0,9,8,a,b,f,e,d,c + pa0  =  [D] + pa0
		,f,e,d,c,a,b,8,9,4,5,7,6,1,0,2,3 + pb0  =  [F] + pb0

		,3,2,1,0,4,5,7,6,f,e,d,c,a,b,8,9 + p70  =  [C] + p70
		,6,7,4,5,3,2,1,0,9,8,a,b,f,e,d,c + p60  =  [D] + p60
		,9,8,a,b,f,e,d,c,3,2,1,0,4,5,7,6 + p50  =  [B] + p50
		,c,d,f,e,9,8,a,b,6,7,4,5,3,2,1,0 + p40  =  [E] + p40,
	i.e. 6 mod-16 index patterns - the unique A, 4xB, 4xC, 4xD, 4xE, 3xF .

	Thus we can encode the needed o_offset 64-perms in (mod p40) fashion in terms of the above 6 [A-F] 16-subperms as

		Block 1: a[]		Block 2: a[]+pg0		Block 3: a[]+pc0		Block 2: a[]+p80		Block 3: a[]+p40
		[A] + p00			[E] + p00				[D] + p20				[B] + p10				[C] + p30
		[B] + p10			[C] + p10				[F] + p30				[E] + p00				[D] + p20
		[C] + p30			[F] + p30				[E] + p00				[D] + p20				[B] + p10
		[D] + p20			[B] + p20				[C] + p10				[F] + p30				[E] + p00
	*/
		//	NOTE that RADIX_64_DIF outputs are IN-ORDER rather than BR:
		jt = j1    ;	RADIX_64_DIF((double *)(t+0x000),t_offsets,1, (a+jt),(o_offsets+0x000),RE_IM_STRIDE);	// Inputs in t[00-63]
		jt = j1+pg0;	RADIX_64_DIF((double *)(t+0x040),t_offsets,1, (a+jt),(o_offsets+0x040),RE_IM_STRIDE);	// Inputs in t[64-127]
		jt = j1+pc0;	RADIX_64_DIF((double *)(t+0x080),t_offsets,1, (a+jt),(o_offsets+0x080),RE_IM_STRIDE);	// Inputs in t[128-191]
		jt = j1+p80;	RADIX_64_DIF((double *)(t+0x0c0),t_offsets,1, (a+jt),(o_offsets+0x0c0),RE_IM_STRIDE);	// Inputs in t[192-255]
		jt = j1+p40;	RADIX_64_DIF((double *)(t+0x100),t_offsets,1, (a+jt),(o_offsets+0x100),RE_IM_STRIDE);	// Inputs in t[256-319]
	}
}

/***************/

void radix320_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-320 complex DIT FFT pass on the data in the length-N real vector A.
*/
	int j,j1,j2,jp,jt,k,l,k0,k1,k2,k3,k4,nshift, lm64,lne0,ic;
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,pg0,ph0,pi0,pj0, first_entry=TRUE;
	static int plo[16];
  #ifndef USE_SSE2
	//static int phi[20];
  #endif
	// Need storage for 4 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 4*9 elts:
	static int dit_p20_cperms[36];
	const double cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
				cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
				s2  =  0.95105651629515357211,	/*  sin(u) */
				ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
				ss2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	static int t_offsets[64], i_offsets[RADIX];
	// Consts containing [in little-endian hex-char fashion] the 9 [A-I] 16-subperms defined in the Input-ordering commentary below:
	const uint64 dit_16_iperms[9] = {
		0x01327654fedcba98ull,	// [A]
		0xfedcba9876543210ull,	// [B]
		0x32105467dcef98abull,	// [C]
		0x98abefcd10236745ull,	// [D]
		0x10236745efcdab89ull,	// [E]
		0x67452301ab89cdfeull,	// [F]
		0xab89cdfe23014576ull,	// [G]
		0xcdfe89ba45760132ull,	// [H]
		0x4576013289bafedcull	// [I]
	};
	uint64 i64;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex t[RADIX], *tptr;
	double rt,it;

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
		p1 = NDIVR;
		p2 = p1 + NDIVR;		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );
		p3 = p2 + NDIVR;		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );
		p4 = p3 + NDIVR;		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );
		p5 = p4 + NDIVR;		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );
		p6 = p5 + NDIVR;		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );
		p7 = p6 + NDIVR;		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );
		p8 = p7 + NDIVR;		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );
		p9 = p8 + NDIVR;		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );
		pa = p9 + NDIVR;		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );
		pb = pa + NDIVR;		pa += ( (pa >> DAT_BITS) << PAD_BITS );
		pc = pb + NDIVR;		pb += ( (pb >> DAT_BITS) << PAD_BITS );
		pd = pc + NDIVR;		pc += ( (pc >> DAT_BITS) << PAD_BITS );
		pe = pd + NDIVR;		pd += ( (pd >> DAT_BITS) << PAD_BITS );
		pf = pe + NDIVR;		pe += ( (pe >> DAT_BITS) << PAD_BITS );
		NDIVR <<= 4;			pf += ( (pf >> DAT_BITS) << PAD_BITS );
		p10 = NDIVR;
		p20 = p10 + NDIVR;		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p30 = p20 + NDIVR;		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p40 = p30 + NDIVR;		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
		p50 = p40 + NDIVR;		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p60 = p50 + NDIVR;		p50 += ( (p50 >> DAT_BITS) << PAD_BITS );
		p70 = p60 + NDIVR;		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		p80 = p70 + NDIVR;		p70 += ( (p70 >> DAT_BITS) << PAD_BITS );
		p90 = p80 + NDIVR;		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		pa0 = p90 + NDIVR;		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pb0 = pa0 + NDIVR;		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pc0 = pb0 + NDIVR;		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		pd0 = pc0 + NDIVR;		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pe0 = pd0 + NDIVR;		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		pf0 = pe0 + NDIVR;		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		pg0 = pf0 + NDIVR;		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
		ph0 = pg0 + NDIVR;		pg0 += ( (pg0 >> DAT_BITS) << PAD_BITS );
		pi0 = ph0 + NDIVR;		ph0 += ( (ph0 >> DAT_BITS) << PAD_BITS );
		pj0 = pi0 + NDIVR;		pi0 += ( (pi0 >> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			pj0 += ( (pj0 >> DAT_BITS) << PAD_BITS );

		plo[0x0] =  0; plo[0x1] = p1; plo[0x2] = p2; plo[0x3] = p3;
		plo[0x4] = p4; plo[0x5] = p5; plo[0x6] = p6; plo[0x7] = p7;
		plo[0x8] = p8; plo[0x9] = p9; plo[0xa] = pa; plo[0xb] = pb;
		plo[0xc] = pc; plo[0xd] = pd; plo[0xe] = pe; plo[0xf] = pf;

	#ifndef USE_SSE2
		/* phi[0x00]=   0; phi[0x01]= p10; phi[0x02]= p20; phi[0x03]= p30; */
		/* phi[0x04]= p40; phi[0x05]= p50; phi[0x06]= p60; phi[0x07]= p70; */
		/* phi[0x08]= p80; phi[0x09]= p90; phi[0x0a]= pa0; phi[0x0b]= pb0; */
		/* phi[0x0c]= pc0; phi[0x0d]= pd0; phi[0x0e]= pe0; phi[0x0f]= pf0; */
		/* phi[0x10]= pg0; phi[0x11]= ph0; phi[0x12]= pi0; phi[0x13]= pj0; */
	#endif

	// Set array offsets for radix-64 outputs.
		// t_offsets w.r.to: t-array, same for all 5 DFTs:
		t_offsets[0x00] = 0x00<<1;	t_offsets[0x10] = 0x10<<1;	t_offsets[0x20] = 0x20<<1;	t_offsets[0x30] = 0x30<<1;
		t_offsets[0x01] = 0x01<<1;	t_offsets[0x11] = 0x11<<1;	t_offsets[0x21] = 0x21<<1;	t_offsets[0x31] = 0x31<<1;
		t_offsets[0x02] = 0x02<<1;	t_offsets[0x12] = 0x12<<1;	t_offsets[0x22] = 0x22<<1;	t_offsets[0x32] = 0x32<<1;
		t_offsets[0x03] = 0x03<<1;	t_offsets[0x13] = 0x13<<1;	t_offsets[0x23] = 0x23<<1;	t_offsets[0x33] = 0x33<<1;
		t_offsets[0x04] = 0x04<<1;	t_offsets[0x14] = 0x14<<1;	t_offsets[0x24] = 0x24<<1;	t_offsets[0x34] = 0x34<<1;
		t_offsets[0x05] = 0x05<<1;	t_offsets[0x15] = 0x15<<1;	t_offsets[0x25] = 0x25<<1;	t_offsets[0x35] = 0x35<<1;
		t_offsets[0x06] = 0x06<<1;	t_offsets[0x16] = 0x16<<1;	t_offsets[0x26] = 0x26<<1;	t_offsets[0x36] = 0x36<<1;
		t_offsets[0x07] = 0x07<<1;	t_offsets[0x17] = 0x17<<1;	t_offsets[0x27] = 0x27<<1;	t_offsets[0x37] = 0x37<<1;
		t_offsets[0x08] = 0x08<<1;	t_offsets[0x18] = 0x18<<1;	t_offsets[0x28] = 0x28<<1;	t_offsets[0x38] = 0x38<<1;
		t_offsets[0x09] = 0x09<<1;	t_offsets[0x19] = 0x19<<1;	t_offsets[0x29] = 0x29<<1;	t_offsets[0x39] = 0x39<<1;
		t_offsets[0x0a] = 0x0a<<1;	t_offsets[0x1a] = 0x1a<<1;	t_offsets[0x2a] = 0x2a<<1;	t_offsets[0x3a] = 0x3a<<1;
		t_offsets[0x0b] = 0x0b<<1;	t_offsets[0x1b] = 0x1b<<1;	t_offsets[0x2b] = 0x2b<<1;	t_offsets[0x3b] = 0x3b<<1;
		t_offsets[0x0c] = 0x0c<<1;	t_offsets[0x1c] = 0x1c<<1;	t_offsets[0x2c] = 0x2c<<1;	t_offsets[0x3c] = 0x3c<<1;
		t_offsets[0x0d] = 0x0d<<1;	t_offsets[0x1d] = 0x1d<<1;	t_offsets[0x2d] = 0x2d<<1;	t_offsets[0x3d] = 0x3d<<1;
		t_offsets[0x0e] = 0x0e<<1;	t_offsets[0x1e] = 0x1e<<1;	t_offsets[0x2e] = 0x2e<<1;	t_offsets[0x3e] = 0x3e<<1;
		t_offsets[0x0f] = 0x0f<<1;	t_offsets[0x1f] = 0x1f<<1;	t_offsets[0x2f] = 0x2f<<1;	t_offsets[0x3f] = 0x3f<<1;

		// Init storage for 4 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 4*9:
		l = 0;
		dit_p20_cperms[l++] =   0; dit_p20_cperms[l++] = pg0; dit_p20_cperms[l++] = pc0; dit_p20_cperms[l++] = p80; dit_p20_cperms[l++] = p40; dit_p20_cperms[l++] =   0; dit_p20_cperms[l++] = pg0; dit_p20_cperms[l++] = pc0; dit_p20_cperms[l++] = p80;
		dit_p20_cperms[l++] = p10; dit_p20_cperms[l++] = ph0; dit_p20_cperms[l++] = pd0; dit_p20_cperms[l++] = p90; dit_p20_cperms[l++] = p50; dit_p20_cperms[l++] = p10; dit_p20_cperms[l++] = ph0; dit_p20_cperms[l++] = pd0; dit_p20_cperms[l++] = p90;
		dit_p20_cperms[l++] = p20; dit_p20_cperms[l++] = pi0; dit_p20_cperms[l++] = pe0; dit_p20_cperms[l++] = pa0; dit_p20_cperms[l++] = p60; dit_p20_cperms[l++] = p20; dit_p20_cperms[l++] = pi0; dit_p20_cperms[l++] = pe0; dit_p20_cperms[l++] = pa0;
		dit_p20_cperms[l++] = p30; dit_p20_cperms[l++] = pj0; dit_p20_cperms[l++] = pf0; dit_p20_cperms[l++] = pb0; dit_p20_cperms[l++] = p70; dit_p20_cperms[l++] = p30; dit_p20_cperms[l++] = pj0; dit_p20_cperms[l++] = pf0; dit_p20_cperms[l++] = pb0;

	/* Iperm 1:
		0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00 = [A] + p00
		f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10 = [B] + p10
		f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p30 = [B] + p30
		f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p20 = [B] + p20
	*/
		i64 = dit_16_iperms[0];	// [A] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;	// p-offset indices encoded in little-endian hex-char fashion:
			i_offsets[     l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dit_16_iperms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0x10+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dit_16_iperms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0x20+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dit_16_iperms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0x30+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
	/* Iperm 2:
		3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p70 = [C] + p30
		3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p60 = [C] + p20
		3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p50 = [C] + p10 (mod p40)
		3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p40 = [C] + p00
	*/
		i64 = dit_16_iperms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0x40+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dit_16_iperms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0x50+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dit_16_iperms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0x60+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dit_16_iperms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0x70+l] = plo[(i64 >> nshift)&0xf];
		}
	/* Iperm 3:
		9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p90 = [D] + p10
		9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p80 = [D] + p00
		9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + pa0 = [D] + p20 (mod p40)
		1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pb0 = [E] + p30
	*/
		i64 = dit_16_iperms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0x80+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dit_16_iperms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0x90+l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dit_16_iperms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0xa0+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dit_16_iperms[4];	// [E] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0xb0+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
	/* Iperm 4:
		6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pe0 = [F] + p20
		a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + pf0 = [G] + p30
		6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pc0 = [F] + p00 (mod p40)
		a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + pd0 = [G] + p10
	*/
		i64 = dit_16_iperms[5];	// [F] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0xc0+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dit_16_iperms[6];	// [G] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0xd0+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dit_16_iperms[5];	// [F] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0xe0+l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dit_16_iperms[6];	// [G] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0xf0+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
	/* Iperm 5:
		c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + pg0 = [H] + p00
		4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + ph0 = [I] + p10
		4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + pj0 = [I] + p30 (mod p40)
		4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + pi0 = [I] + p20
	*/
		i64 = dit_16_iperms[7];	// [H] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0x100+l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dit_16_iperms[8];	// [i] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0x110+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dit_16_iperms[8];	// [I] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0x120+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dit_16_iperms[8];	// [I] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			i_offsets[0x130+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
	}

/*...The radix-320 pass is here.	*/

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

	//...gather the needed data (320 64-bit complex) and do 3 radix-64 transforms:
	/*
	Twiddleless version uses same linear-index-vector-form permutation as in DIF -
	Remember, inputs to DIT are bit-reversed, so using output of test_fft_radix(), we have

	Combined DIT input-scramble array =
		  0,  1,  3,  2,  7,  6,  5,  4,  f,  e,  d,  c,  b,  a,  9,  8 = 0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00 = [A] + p00
		 1f, 1e, 1d, 1c, 1b, 1a, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10 = f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10 = [B] + p10
		 3f, 3e, 3d, 3c, 3b, 3a, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30 = f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p30 = [B] + p30
		 2f, 2e, 2d, 2c, 2b, 2a, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20 = f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p20 = [B] + p20
		 73, 72, 71, 70, 75, 74, 76, 77, 7d, 7c, 7e, 7f, 79, 78, 7a, 7b = 3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p70 = [C] + p70
		 63, 62, 61, 60, 65, 64, 66, 67, 6d, 6c, 6e, 6f, 69, 68, 6a, 6b = 3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p60 = [C] + p60
		 53, 52, 51, 50, 55, 54, 56, 57, 5d, 5c, 5e, 5f, 59, 58, 5a, 5b = 3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p50 = [C] + p50
		 43, 42, 41, 40, 45, 44, 46, 47, 4d, 4c, 4e, 4f, 49, 48, 4a, 4b = 3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p40 = [C] + p40
		 99, 98, 9a, 9b, 9e, 9f, 9c, 9d, 91, 90, 92, 93, 96, 97, 94, 95 = 9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p90 = [D] + p90
		 89, 88, 8a, 8b, 8e, 8f, 8c, 8d, 81, 80, 82, 83, 86, 87, 84, 85 = 9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p80 = [D] + p80
		 a9, a8, aa, ab, ae, af, ac, ad, a1, a0, a2, a3, a6, a7, a4, a5 = 9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + pa0 = [D] + pa0
		 b1, b0, b2, b3, b6, b7, b4, b5, be, bf, bc, bd, ba, bb, b8, b9 = 1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pb0 = [E] + pb0
		 e6, e7, e4, e5, e2, e3, e0, e1, ea, eb, e8, e9, ec, ed, ef, ee = 6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pe0 = [F] + pe0
		 fa, fb, f8, f9, fc, fd, ff, fe, f2, f3, f0, f1, f4, f5, f7, f6 = a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + pf0 = [G] + pf0
		 c6, c7, c4, c5, c2, c3, c0, c1, ca, cb, c8, c9, cc, cd, cf, ce = 6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pc0 = [F] + pc0
		 da, db, d8, d9, dc, dd, df, de, d2, d3, d0, d1, d4, d5, d7, d6 = a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + pd0 = [G] + pd0
		10c,10d,10f,10e,108,109,10b,10a,104,105,107,106,100,101,103,102 = c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + pg0 = [H] + pg0
		114,115,117,116,110,111,113,112,118,119,11b,11a,11f,11e,11d,11c = 4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + ph0 = [I] + ph0
		134,135,137,136,130,131,133,132,138,139,13b,13a,13f,13e,13d,13c = 4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + pj0 = [I] + pj0
		124,125,127,126,120,121,123,122,128,129,12b,12a,12f,12e,12d,12c = 4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + pi0 = [I] + pi0
	*/
		jt = j1    ;	RADIX_64_DIT(a+jt,(i_offsets+0x000),RE_IM_STRIDE, (double *)(t+0x000),t_offsets,1);	// Outputs in t[00-63]
		jt = j1+p40;	RADIX_64_DIT(a+jt,(i_offsets+0x040),RE_IM_STRIDE, (double *)(t+0x040),t_offsets,1);	// Outputs in t[64-127]
		jt = j1+p80;	RADIX_64_DIT(a+jt,(i_offsets+0x080),RE_IM_STRIDE, (double *)(t+0x080),t_offsets,1);	// Outputs in t[64-127]
		jt = j1+pc0;	RADIX_64_DIT(a+jt,(i_offsets+0x0c0),RE_IM_STRIDE, (double *)(t+0x0c0),t_offsets,1);	// Outputs in t[64-127]
		jt = j1+pg0;	RADIX_64_DIT(a+jt,(i_offsets+0x100),RE_IM_STRIDE, (double *)(t+0x100),t_offsets,1);	// Outputs in t[64-127]

	/*...and now do 64 radix-5 transforms, with the columns of t*[r,i] output pairs in the above 5x radix-64 set now acting as input rows.
	Since our first-look oindex ordering was +p0x[0,40,80,c0,100] for each radix-5 and incrementing += p1 between those DFTs,
	arrange resulting mismatched-data-sorted index permutation into 5 vertical 64-entry columns to get needed oindex patterns.
	Indexing in hex for clarity and using [0-4]<<{0-4} notation in the rightmost column to flag reusable 5-perms, as left-circular
	(0-4)-elt shifts of basic perms [0] := 0,100,c0,80,40, [1] := 10,110,d0,90,50, [2] := 20,120,e0,a0,60, [3] := 30,130,f0,b0,70:
		Required output index permutation = [
		  0,100, c0, 80, 40		  0,100, c0, 80, 40 + p0		[0]<<0 + p0
		 ff, bf, 7f, 3f,13f		 f0, b0, 70, 30,130 + pf		[3]<<2 + pf
		 be, 7e, 3e,13e, fe		 b0, 70, 30,130, f0 + pe		[3]<<3 + pe
		 7d, 3d,13d, fd, bd		 70, 30,130, f0, b0 + pd		[3]<<4 + pd
		 3c,13c, fc, bc, 7c		 30,130, f0, b0, 70 + pc		[3]<<0 + pc
		13b, fb, bb, 7b, 3b		130, f0, b0, 70, 30 + pb		[3]<<1 + pb
		 fa, ba, 7a, 3a,13a		 f0, b0, 70, 30,130 + pa		[3]<<2 + pa
		 b9, 79, 39,139, f9		 b0, 70, 30,130, f0 + p9		[3]<<3 + p9
		 78, 38,138, f8, b8		 70, 30,130, f0, b0 + p8		[3]<<4 + p8
		 37,137, f7, b7, 77		 30,130, f0, b0, 70 + p7		[3]<<0 + p7
		136, f6, b6, 76, 36		130, f0, b0, 70, 30 + p6		[3]<<1 + p6
		 f5, b5, 75, 35,135		 f0, b0, 70, 30,130 + p5		[3]<<2 + p5
		 b4, 74, 34,134, f4		 b0, 70, 30,130, f0 + p4		[3]<<3 + p4
		 73, 33,133, f3, b3		 70, 30,130, f0, b0 + p3		[3]<<4 + p3
		 32,132, f2, b2, 72		 30,130, f0, b0, 70 + p2		[3]<<0 + p2
		131, f1, b1, 71, 31		130, f0, b0, 70, 30 + p1		[3]<<1 + p1
		 f0, b0, 70, 30,130	=	 f0, b0, 70, 30,130 + p0	=	[3]<<2 + p0
		 af, 6f, 2f,12f, ef		 a0, 60, 20,120, e0 + pf		[2]<<3 + pf
		 6e, 2e,12e, ee, ae		 60, 20,120, e0, a0 + pe		[2]<<4 + pe
		 2d,12d, ed, ad, 6d		 20,120, e0, a0, 60 + pd		[2]<<0 + pd
		12c, ec, ac, 6c, 2c		120, e0, a0, 60, 20 + pc		[2]<<1 + pc
		 eb, ab, 6b, 2b,12b		 e0, a0, 60, 20,120 + pb		[2]<<2 + pb
		 aa, 6a, 2a,12a, ea		 a0, 60, 20,120, e0 + pa		[2]<<3 + pa
		 69, 29,129, e9, a9		 60, 20,120, e0, a0 + p9		[2]<<4 + p9
		 28,128, e8, a8, 68		 20,120, e0, a0, 60 + p8		[2]<<0 + p8
		127, e7, a7, 67, 27		120, e0, a0, 60, 20 + p7		[2]<<1 + p7
		 e6, a6, 66, 26,126		 e0, a0, 60, 20,120 + p6		[2]<<2 + p6
		 a5, 65, 25,125, e5		 a0, 60, 20,120, e0 + p5		[2]<<3 + p5
		 64, 24,124, e4, a4		 60, 20,120, e0, a0 + p4		[2]<<4 + p4
		 23,123, e3, a3, 63		 20,120, e0, a0, 60 + p3		[2]<<0 + p3
		122, e2, a2, 62, 22		120, e0, a0, 60, 20 + p2		[2]<<1 + p2
		 e1, a1, 61, 21,121		 e0, a0, 60, 20,120 + p1		[2]<<2 + p1
		 a0, 60, 20,120, e0		 a0, 60, 20,120, e0 + p0		[2]<<3 + p0
		 5f, 1f,11f, df, 9f		 50, 10,110, d0, 90 + pf		[1]<<4 + pf
		 1e,11e, de, 9e, 5e		 10,110, d0, 90, 50 + pe		[1]<<0 + pe
		11d, dd, 9d, 5d, 1d		110, d0, 90, 50, 10 + pd		[1]<<1 + pd
		 dc, 9c, 5c, 1c,11c		 d0, 90, 50, 10,110 + pc		[1]<<2 + pc
		 9b, 5b, 1b,11b, db		 90, 50, 10,110, d0 + pb		[1]<<3 + pb
		 5a, 1a,11a, da, 9a		 50, 10,110, d0, 90 + pa		[1]<<4 + pa
		 19,119, d9, 99, 59		 10,110, d0, 90, 50 + p9		[1]<<0 + p9
		118, d8, 98, 58, 18		110, d0, 90, 50, 10 + p8		[1]<<1 + p8
		 d7, 97, 57, 17,117		 d0, 90, 50, 10,110 + p7		[1]<<2 + p7
		 96, 56, 16,116, d6		 90, 50, 10,110, d0 + p6		[1]<<3 + p6
		 55, 15,115, d5, 95		 50, 10,110, d0, 90 + p5		[1]<<4 + p5
		 14,114, d4, 94, 54		 10,110, d0, 90, 50 + p4		[1]<<0 + p4
		113, d3, 93, 53, 13		110, d0, 90, 50, 10 + p3		[1]<<1 + p3
		 d2, 92, 52, 12,112		 d0, 90, 50, 10,110 + p2		[1]<<2 + p2
		 91, 51, 11,111, d1		 90, 50, 10,110, d0 + p1		[1]<<3 + p1
		 50, 10,110, d0, 90	=	 50, 10,110, d0, 90 + p0	=	[1]<<4 + p0
		  f,10f, cf, 8f, 4f		  0,100, c0, 80, 40 + pf		[0]<<0 + pf
		10e, ce, 8e, 4e,  e		100, c0, 80, 40,  0 + pe		[0]<<1 + pe
		 cd, 8d, 4d,  d,10d		 c0, 80, 40,  0,100 + pd		[0]<<2 + pd
		 8c, 4c,  c,10c, cc		 80, 40,  0,100, c0 + pc		[0]<<3 + pc
		 4b,  b,10b, cb, 8b		 40,  0,100, c0, 80 + pb		[0]<<4 + pb
		  a,10a, ca, 8a, 4a		  0,100, c0, 80, 40 + pa		[0]<<0 + pa
		109, c9, 89, 49,  9		100, c0, 80, 40,  0 + p9		[0]<<1 + p9
		 c8, 88, 48,  8,108		 c0, 80, 40,  0,100 + p8		[0]<<2 + p8
		 87, 47,  7,107, c7		 80, 40,  0,100, c0 + p7		[0]<<3 + p7
		 46,  6,106, c6, 86		 40,  0,100, c0, 80 + p6		[0]<<4 + p6
		  5,105, c5, 85, 45		  0,100, c0, 80, 40 + p5		[0]<<0 + p5
		104, c4, 84, 44,  4		100, c0, 80, 40,  0 + p4		[0]<<1 + p4
		 c3, 83, 43,  3,103		 c0, 80, 40,  0,100 + p3		[0]<<2 + p3
		 82, 42,  2,102, c2		 80, 40,  0,100, c0 + p2		[0]<<3 + p2
		 41,  1,101, c1, 81		 40,  0,100, c0, 80 + p1		[0]<<4 + p1

		Leftward cshift count satisfies	(i + (i!=0)) % 5
		Which-5-perm-index satisfies	((64-i)>>2) & (-(i!=0))
		p-offset index satisfies		(64-i) % 16
	*/
		tptr = t;
		for(l = 0; l < 64; l++) {
		#if 0	// Initial try, used to compute needed output-permutation:
			k = plo[l&15] + (phi[l>>4]);
			k0 =   0;	k1 = p40;	k2 = p80;	k3 = pc0;	k4 = pg0;
		#else
			lm64 = 64-l, lne0 = (l != 0);
			k  = (lm64>>4) & (-lne0);	k += (k<<3);// Which-5-perm-index, *= 9
			ic = ((l + lne0) % 5) + k;			// Leftward cshift count plus 9-elts offset per increment in k
//printf("Set %2u: perm[%u] << %u + p%X, ic = %u",l,k/9,(l + lne0) % 5,lm64 & 15, ic);
			k = plo[lm64 & 15];						// Re-use k for p-offset index, plo[(64-i) % 16]
			k0 = dit_p20_cperms[ic]; k1 = dit_p20_cperms[ic+1]; k2 = dit_p20_cperms[ic+2]; k3 = dit_p20_cperms[ic+3]; k4 = dit_p20_cperms[ic+4];
		#endif
//printf(" ... k0-4 = %2X,%2X,%2X,%2X,%2X\n",k0>>1,k1>>1,k2>>1,k3>>1,k4>>1);
			jt = j1+k; jp = j2+k;
			RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,
				tptr->re,tptr->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0x100)->re,(tptr+0x100)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],
				rt,it);
			tptr++;
		}
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy320_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr;
	#ifdef USE_SSE2
		const int pfetch_dist = PFETCH_DIST;
	#endif
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
			,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,pg0,ph0,pi0,pj0;
	// Shared DIF+DIT:
	#ifndef USE_SSE2
		double rt,it;
		double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	#endif
		int poff[RADIX>>2];	// Store [RADIX/4] mults of p04 offset for loop control
		int ilo,ihi,hi_neg,k0,k1,k2,k3,k4,nshift;
		int plo[16], t_offsets[64];
	#ifndef USE_SSE2
		int phi[20];
	#endif
		uint64 i0,i1,i2,i3;
	// DIF:
		int dif_o_offsets[RADIX];
		// Consts containing [in little-endian hex-char fashion] the 6 [A-F] 16-subperms defined in the output-ordering commentary below:
		const uint64 dif_16_operms[6] = {
			0x01326745cdfe98abull,	// [A]
			0x98abfedc32104576ull,	// [B]
			0x32104576fedcab89ull,	// [C]
			0x6745321098abfedcull,	// [D]
			0xcdfe98ab67453210ull,	// [E]
			0xfedcab8945761023ull	// [F]
		};
	// DIT:
		// Need storage for 4 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 4*9 elts:
		int dit_p20_cperms[36];
		int dit_i_offsets[RADIX];
		// Consts containing [in little-endian hex-char fashion] the 9 [A-I] 16-subperms defined in the Input-ordering commentary below:
		const uint64 dit_16_iperms[9] = {
			0x01327654fedcba98ull,	// [A]
			0xfedcba9876543210ull,	// [B]
			0x32105467dcef98abull,	// [C]
			0x98abefcd10236745ull,	// [D]
			0x10236745efcdab89ull,	// [E]
			0x67452301ab89cdfeull,	// [F]
			0xab89cdfe23014576ull,	// [G]
			0xcdfe89ba45760132ull,	// [H]
			0x4576013289bafedcull	// [I]
		};

		int j,j1,jt,l;
	#ifndef USE_SSE2
		int j2,jp,ntmp;
	#endif
	#ifdef USE_SSE2
		// incr = Carry-chain wts-multipliers recurrence length, which must divide
		// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 20|40|80 for avx512,avx,sse, respectively:
		// Note: Not able to get away with as large a setting of incr for radix 320 as we did for 160.
		int incr;
	  #ifdef USE_AVX512
		const int incr_long = 10,incr_med = 5,incr_short = 5;
	  #elif defined(USE_AVX2)
		const int incr_long = 20,incr_med =10,incr_short = 4;
	  #else
		const int incr_long = 20,incr_med =10,incr_short = 4;
	  #endif
	  // Have no specialized HIACC carry macro in USE_AVX512 and ARMv8 SIMD, so use "goes to 11" in LOACC mode via an incr_hiacc = 1:
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
	#endif // USE_SSE2

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
		int *itmp;	// Pointer into the bjmodn array
	  #if defined(USE_AVX) && !defined(USE_AVX512)
		int *itm2;	// Pointer into the bjmodn array
	  #endif
	  #ifndef USE_AVX
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	  #endif
		double *add0,*add1,*add2,*add3;
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *tmp,*tm1,*tm2,	// Non-static utility ptrs
			*va0,*va1,*va2,*va3,*va4, *vb0,*vb1,*vb2,*vb3,*vb4
		  #if defined(USE_AVX2) && !defined(USE_AVX512)
			,*wa0,*wa1,*wa2,*wa3,*wa4, *wb0,*wb1,*wb2,*wb3,*wb4
		  #endif
			;
		vec_dbl
		  #if defined(USE_AVX2) && !defined(USE_AVX512)
			*two,
		  #endif
			*ycc1,/* *yss1,*ycc2,*yss2,*yss3, */	// radiy-5 DFT trig consts
			*max_err,
		  #ifndef USE_AVX512
			*sse2_rnd,
		  #endif
			*half_arr,
			*r00,	// Head of RADIX*vec_cmplx-sized local store #1
			*s1p00,	// Head of RADIX*vec_cmplx-sized local store #2
			*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx
	  #ifndef USE_AVX512
		double dtmp;
	  #endif
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else

		const double cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
				cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
				s2  =  0.95105651629515357211,	/*  sin(u) */
				ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
				ss2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
		double *base, *baseinv;
		int p0123[4];
		const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int m,m2;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		int bjmodn[RADIX];	// Thread only carries a base datum here, must alloc a local array for remaining values
		double *cy = thread_arg->cy, temp,frac;
		struct complex t[RADIX], *tptr;
		int *itmp;	// Pointer into the bjmodn array

	#endif

	// int data:
	#ifdef USE_SSE2
		int thr_id = thread_arg->tid;
	#endif
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
		p1 = NDIVR;
		p2 = p1 + NDIVR;		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );
		p3 = p2 + NDIVR;		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );
		p4 = p3 + NDIVR;		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );
		p5 = p4 + NDIVR;		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );
		p6 = p5 + NDIVR;		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );
		p7 = p6 + NDIVR;		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );
		p8 = p7 + NDIVR;		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );
		p9 = p8 + NDIVR;		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );
		pa = p9 + NDIVR;		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );
		pb = pa + NDIVR;		pa += ( (pa >> DAT_BITS) << PAD_BITS );
		pc = pb + NDIVR;		pb += ( (pb >> DAT_BITS) << PAD_BITS );
		pd = pc + NDIVR;		pc += ( (pc >> DAT_BITS) << PAD_BITS );
		pe = pd + NDIVR;		pd += ( (pd >> DAT_BITS) << PAD_BITS );
		pf = pe + NDIVR;		pe += ( (pe >> DAT_BITS) << PAD_BITS );
		NDIVR <<= 4;			pf += ( (pf >> DAT_BITS) << PAD_BITS );
		p10 = NDIVR;
		p20 = p10 + NDIVR;		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p30 = p20 + NDIVR;		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p40 = p30 + NDIVR;		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
		p50 = p40 + NDIVR;		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p60 = p50 + NDIVR;		p50 += ( (p50 >> DAT_BITS) << PAD_BITS );
		p70 = p60 + NDIVR;		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		p80 = p70 + NDIVR;		p70 += ( (p70 >> DAT_BITS) << PAD_BITS );
		p90 = p80 + NDIVR;		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		pa0 = p90 + NDIVR;		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pb0 = pa0 + NDIVR;		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pc0 = pb0 + NDIVR;		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		pd0 = pc0 + NDIVR;		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pe0 = pd0 + NDIVR;		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		pf0 = pe0 + NDIVR;		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		pg0 = pf0 + NDIVR;		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
		ph0 = pg0 + NDIVR;		pg0 += ( (pg0 >> DAT_BITS) << PAD_BITS );
		pi0 = ph0 + NDIVR;		ph0 += ( (ph0 >> DAT_BITS) << PAD_BITS );
		pj0 = pi0 + NDIVR;		pi0 += ( (pi0 >> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			pj0 += ( (pj0 >> DAT_BITS) << PAD_BITS );
	#ifndef USE_SSE2
		p0123[0] = 0; p0123[1] = p1; p0123[2] = p2; p0123[3] = p3;
	#endif
		poff[     0] =   0; poff[0x04+0] = p10; poff[0x08+0] = p20; poff[0x0c+0] = p30;
		poff[0x10+0] = p40; poff[0x14+0] = p50; poff[0x18+0] = p60; poff[0x1c+0] = p70;
		poff[0x20+0] = p80; poff[0x24+0] = p90; poff[0x28+0] = pa0; poff[0x2c+0] = pb0;
		poff[0x30+0] = pc0; poff[0x34+0] = pd0; poff[0x38+0] = pe0; poff[0x3c+0] = pf0;
		poff[0x40+0] = pg0; poff[0x44+0] = ph0; poff[0x48+0] = pi0; poff[0x4c+0] = pj0;
		for(l = 0; l < 80; l += 4) {
			poff[l+1] = poff[l] + p4; poff[l+2] = poff[l] + p8; poff[l+3] = poff[l] + pc;
		}
	// Shared:
		// Set array offsets for radix-64 DFT in/outputs:
		// t_offsets w.r.to: t-array, same for all 5 DFTs:
		for(l = 0; l < 64; l += 4) {
			t_offsets[l] = l+l; t_offsets[l+1] = (l+1)<<1; t_offsets[l+2] = (l+2)<<1; t_offsets[l+3] = (l+3)<<1;
		}

		plo[0x0] =  0; plo[0x1] = p1; plo[0x2] = p2; plo[0x3] = p3;
		plo[0x4] = p4; plo[0x5] = p5; plo[0x6] = p6; plo[0x7] = p7;
		plo[0x8] = p8; plo[0x9] = p9; plo[0xa] = pa; plo[0xb] = pb;
		plo[0xc] = pc; plo[0xd] = pd; plo[0xe] = pe; plo[0xf] = pf;

		/*** NB: Delay the plo[]-analog of the phi[] SIMD fixed-addressing to *after* the inits below
		because DIF o_offset and DIT i_offset-array inits need the p-multiplied versions of plo,phi ***/
	   #ifdef USE_SSE2
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p1*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:
		/* phi[0x00]= 0x00<<1; phi[0x01]= 0x10<<1; phi[0x02]= 0x20<<1; phi[0x03]= 0x30<<1; */
		/* phi[0x04]= 0x40<<1; phi[0x05]= 0x50<<1; phi[0x06]= 0x60<<1; phi[0x07]= 0x70<<1; */
		/* phi[0x08]= 0x80<<1; phi[0x09]= 0x90<<1; phi[0x0a]= 0xa0<<1; phi[0x0b]= 0xb0<<1; */
		/* phi[0x0c]= 0xc0<<1; phi[0x0d]= 0xd0<<1; phi[0x0e]= 0xe0<<1; phi[0x0f]= 0xf0<<1; */
		/* phi[0x10]=0x100<<1; phi[0x11]=0x110<<1; phi[0x12]=0x120<<1; phi[0x13]=0x130<<1; */
	   #else
		phi[0x00]=   0; phi[0x01]= p10; phi[0x02]= p20; phi[0x03]= p30;
		phi[0x04]= p40; phi[0x05]= p50; phi[0x06]= p60; phi[0x07]= p70;
		phi[0x08]= p80; phi[0x09]= p90; phi[0x0a]= pa0; phi[0x0b]= pb0;
		phi[0x0c]= pc0; phi[0x0d]= pd0; phi[0x0e]= pe0; phi[0x0f]= pf0;
		phi[0x10]= pg0; phi[0x11]= ph0; phi[0x12]= pi0; phi[0x13]= pj0;
	   #endif	// sse2?

	/*** DIF Output permutation: */
	// Operm 1:
		i0 = dif_16_operms[0]; i1 = dif_16_operms[1]; i2 = dif_16_operms[2]; i3 = dif_16_operms[3];	// [A-D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;	// p-offset indices encoded in little-endian hex-char fashion:
			dif_o_offsets[     l] = plo[(i0 >> nshift)&0xf];
			dif_o_offsets[0x10+l] = plo[(i1 >> nshift)&0xf] + p10;
			dif_o_offsets[0x20+l] = plo[(i2 >> nshift)&0xf] + p30;
			dif_o_offsets[0x30+l] = plo[(i3 >> nshift)&0xf] + p20;
		}
	// Operm 2:
		i0 = dif_16_operms[4]; i1 = dif_16_operms[2]; i2 = dif_16_operms[5]; i3 = dif_16_operms[1];	// [E,C,F,B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x40+l] = plo[(i0 >> nshift)&0xf];
			dif_o_offsets[0x50+l] = plo[(i1 >> nshift)&0xf] + p10;
			dif_o_offsets[0x60+l] = plo[(i2 >> nshift)&0xf] + p30;
			dif_o_offsets[0x70+l] = plo[(i3 >> nshift)&0xf] + p20;
		}
	// Operm 3:
		i0 = dif_16_operms[3]; i1 = dif_16_operms[5]; i2 = dif_16_operms[4]; i3 = dif_16_operms[2];	// [D,F,E,C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x80+l] = plo[(i0 >> nshift)&0xf] + p20;
			dif_o_offsets[0x90+l] = plo[(i1 >> nshift)&0xf] + p30;
			dif_o_offsets[0xa0+l] = plo[(i2 >> nshift)&0xf];
			dif_o_offsets[0xb0+l] = plo[(i3 >> nshift)&0xf] + p10;
		}
	// Operm 4:
		i0 = dif_16_operms[1]; i1 = dif_16_operms[4]; i2 = dif_16_operms[3]; i3 = dif_16_operms[5];	// [B,E,D,F] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0xc0+l] = plo[(i0 >> nshift)&0xf] + p10;
			dif_o_offsets[0xd0+l] = plo[(i1 >> nshift)&0xf];
			dif_o_offsets[0xe0+l] = plo[(i2 >> nshift)&0xf] + p20;
			dif_o_offsets[0xf0+l] = plo[(i3 >> nshift)&0xf] + p30;
		}
	// Operm 5:
		i0 = dif_16_operms[2]; i1 = dif_16_operms[3]; i2 = dif_16_operms[1]; i3 = dif_16_operms[4];	// [C,D,B,E] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x100+l] = plo[(i0 >> nshift)&0xf] + p30;
			dif_o_offsets[0x110+l] = plo[(i1 >> nshift)&0xf] + p20;
			dif_o_offsets[0x120+l] = plo[(i2 >> nshift)&0xf] + p10;
			dif_o_offsets[0x130+l] = plo[(i3 >> nshift)&0xf];
		}

	/*** DIT indexing stuff: ***/
		// Init storage for 4 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 4*9:
		l = 0;
	  #ifdef USE_SSE2
		dit_p20_cperms[l++] = 0x000<<1; dit_p20_cperms[l++] = 0x100<<1; dit_p20_cperms[l++] = 0x0c0<<1; dit_p20_cperms[l++] = 0x080<<1; dit_p20_cperms[l++] = 0x040<<1; dit_p20_cperms[l++] = 0x000<<1; dit_p20_cperms[l++] = 0x100<<1; dit_p20_cperms[l++] = 0x0c0<<1; dit_p20_cperms[l++] = 0x080<<1;
		dit_p20_cperms[l++] = 0x010<<1; dit_p20_cperms[l++] = 0x110<<1; dit_p20_cperms[l++] = 0x0d0<<1; dit_p20_cperms[l++] = 0x090<<1; dit_p20_cperms[l++] = 0x050<<1; dit_p20_cperms[l++] = 0x010<<1; dit_p20_cperms[l++] = 0x110<<1; dit_p20_cperms[l++] = 0x0d0<<1; dit_p20_cperms[l++] = 0x090<<1;
		dit_p20_cperms[l++] = 0x020<<1; dit_p20_cperms[l++] = 0x120<<1; dit_p20_cperms[l++] = 0x0e0<<1; dit_p20_cperms[l++] = 0x0a0<<1; dit_p20_cperms[l++] = 0x060<<1; dit_p20_cperms[l++] = 0x020<<1; dit_p20_cperms[l++] = 0x120<<1; dit_p20_cperms[l++] = 0x0e0<<1; dit_p20_cperms[l++] = 0x0a0<<1;
		dit_p20_cperms[l++] = 0x030<<1; dit_p20_cperms[l++] = 0x130<<1; dit_p20_cperms[l++] = 0x0f0<<1; dit_p20_cperms[l++] = 0x0b0<<1; dit_p20_cperms[l++] = 0x070<<1; dit_p20_cperms[l++] = 0x030<<1; dit_p20_cperms[l++] = 0x130<<1; dit_p20_cperms[l++] = 0x0f0<<1; dit_p20_cperms[l++] = 0x0b0<<1;
	  #else
		dit_p20_cperms[l++] =   0; dit_p20_cperms[l++] = pg0; dit_p20_cperms[l++] = pc0; dit_p20_cperms[l++] = p80; dit_p20_cperms[l++] = p40; dit_p20_cperms[l++] =   0; dit_p20_cperms[l++] = pg0; dit_p20_cperms[l++] = pc0; dit_p20_cperms[l++] = p80;
		dit_p20_cperms[l++] = p10; dit_p20_cperms[l++] = ph0; dit_p20_cperms[l++] = pd0; dit_p20_cperms[l++] = p90; dit_p20_cperms[l++] = p50; dit_p20_cperms[l++] = p10; dit_p20_cperms[l++] = ph0; dit_p20_cperms[l++] = pd0; dit_p20_cperms[l++] = p90;
		dit_p20_cperms[l++] = p20; dit_p20_cperms[l++] = pi0; dit_p20_cperms[l++] = pe0; dit_p20_cperms[l++] = pa0; dit_p20_cperms[l++] = p60; dit_p20_cperms[l++] = p20; dit_p20_cperms[l++] = pi0; dit_p20_cperms[l++] = pe0; dit_p20_cperms[l++] = pa0;
		dit_p20_cperms[l++] = p30; dit_p20_cperms[l++] = pj0; dit_p20_cperms[l++] = pf0; dit_p20_cperms[l++] = pb0; dit_p20_cperms[l++] = p70; dit_p20_cperms[l++] = p30; dit_p20_cperms[l++] = pj0; dit_p20_cperms[l++] = pf0; dit_p20_cperms[l++] = pb0;
	  #endif
	// Iperm 1:
		i0 = dit_16_iperms[0]; i1 = dit_16_iperms[1]; i2 = dit_16_iperms[1]; i3 = dit_16_iperms[1];	// [A,B,B,B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;	// p-offset indices encoded in little-endian hex-char fashion:
			dit_i_offsets[     l] = plo[(i0 >> nshift)&0xf];
			dit_i_offsets[0x10+l] = plo[(i1 >> nshift)&0xf] + p10;
			dit_i_offsets[0x20+l] = plo[(i2 >> nshift)&0xf] + p30;
			dit_i_offsets[0x30+l] = plo[(i3 >> nshift)&0xf] + p20;
		}
	// Iperm 2:
		i0 = dit_16_iperms[2]; i1 = dit_16_iperms[2]; i2 = dit_16_iperms[2]; i3 = dit_16_iperms[2];	// [C,C,C,C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0x40+l] = plo[(i0 >> nshift)&0xf] + p30;
			dit_i_offsets[0x50+l] = plo[(i1 >> nshift)&0xf] + p20;
			dit_i_offsets[0x60+l] = plo[(i2 >> nshift)&0xf] + p10;
			dit_i_offsets[0x70+l] = plo[(i3 >> nshift)&0xf];
		}
	// Iperm 3:
		i0 = dit_16_iperms[3]; i1 = dit_16_iperms[3]; i2 = dit_16_iperms[3]; i3 = dit_16_iperms[4];	// [D,D,D,E] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0x80+l] = plo[(i0 >> nshift)&0xf] + p10;
			dit_i_offsets[0x90+l] = plo[(i1 >> nshift)&0xf];
			dit_i_offsets[0xa0+l] = plo[(i2 >> nshift)&0xf] + p20;
			dit_i_offsets[0xb0+l] = plo[(i3 >> nshift)&0xf] + p30;
		}
	// Iperm 4:
		i0 = dit_16_iperms[5]; i1 = dit_16_iperms[6]; i2 = dit_16_iperms[5]; i3 = dit_16_iperms[6];	// [F,G,F,G] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0xc0+l] = plo[(i0 >> nshift)&0xf] + p20;
			dit_i_offsets[0xd0+l] = plo[(i1 >> nshift)&0xf] + p30;
			dit_i_offsets[0xe0+l] = plo[(i2 >> nshift)&0xf];
			dit_i_offsets[0xf0+l] = plo[(i3 >> nshift)&0xf] + p10;
		}
	// Iperm 5:
		i0 = dit_16_iperms[7]; i1 = dit_16_iperms[8]; i2 = dit_16_iperms[8]; i3 = dit_16_iperms[8];	// [H,I,I,I] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dit_i_offsets[0x100+l] = plo[(i0 >> nshift)&0xf];
			dit_i_offsets[0x110+l] = plo[(i1 >> nshift)&0xf] + p10;
			dit_i_offsets[0x120+l] = plo[(i2 >> nshift)&0xf] + p30;
			dit_i_offsets[0x130+l] = plo[(i3 >> nshift)&0xf] + p20;
		}

	#ifdef USE_SSE2
		/*** Delayed part of SIMD fixed-offsets init: ***/
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p1*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:
		plo[0x0] = 0x0<<1; plo[0x1] = 0x1<<1; plo[0x2] = 0x2<<1; plo[0x3] = 0x3<<1;
		plo[0x4] = 0x4<<1; plo[0x5] = 0x5<<1; plo[0x6] = 0x6<<1; plo[0x7] = 0x7<<1;
		plo[0x8] = 0x8<<1; plo[0x9] = 0x9<<1; plo[0xa] = 0xa<<1; plo[0xb] = 0xb<<1;
		plo[0xc] = 0xc<<1; plo[0xd] = 0xd<<1; plo[0xe] = 0xe<<1; plo[0xf] = 0xf<<1;
	#endif

	#ifdef USE_SSE2
		tmp	= r00 = thread_arg->r00;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x280;	s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
		tmp += 0x280;
	  #if defined(USE_AVX2) && !defined(USE_AVX512)
		two   = tmp + 0x0;
	  #endif
		// DFT-64 roots all def'd locally in the DFT-64 functions, only need the radix-5 roots here:
		ycc1  = tmp + 0x1;	// radix-5 DFT trig consts
		//ycc2  = tmp + 0x2;
		//yss1  = tmp + 0x3;
		//yss2  = tmp + 0x4;
		//yss3  = tmp + 0x5;
		tmp += 0x6;	// sc_ptr += 0x506; add a pad so the tmp += 2 below gives half_arr = 0x508 + (RADIX/vec_dbl-length)
	  #ifdef USE_AVX512
		cy = tmp;		tmp += 0x28;	// RADIX/8 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		//sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;	// sc_ptr += 0x(508 + 28) = 0x538; This is where the value of half_arr_offset320 comes from
	  #elif defined(USE_AVX)
		cy = tmp;		tmp += 0x50;	// RADIX/4 vec_dbl slots
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(508 + 50) = 0x560; This is where the value of half_arr_offset320 comes from
		half_arr= tmp + 0x02;	// This table needs 20 vec_dbl in both avx and sse2 mode
	  #else
		cy = tmp;		tmp += 0xa0;	// RADIX/2 vec_dbl slots
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(508 + a0) = 0x5b0; This is where the value of half_arr_offset320 comes from
		half_arr= tmp + 0x02;
	  #endif

		ASSERT((r00 == thread_arg->r00), "thread-local memcheck failed!");
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

		sign_mask = (uint64*)(r00 + radix320_creals_in_local_store);
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
		base    = (double *)thread_arg->r00  ;
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
		#include "radix320_main_carry_loop.h"

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
#undef ODD_RADIX
#undef PFETCH_DIST
