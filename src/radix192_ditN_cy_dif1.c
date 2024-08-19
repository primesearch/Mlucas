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

#define RADIX 192	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

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

  // For Mersenne-mod need (16 [SSE2] or 64 [AVX]) + (4 [HIACC] or 40 [LOACC]) added slots for half_arr lookup tables.
  // Max = (40 [SSE2]; 132 [AVX]), add to (half_arr_offset192 + RADIX) to get required value of radix192_creals_in_local_store:
  #ifdef USE_AVX512	// RADIX/8 = 0x18 fewer carry slots than AVX:
	const int half_arr_offset192 = 0x304 + (RADIX>>3);	// RADIX/8 vec_dbl slots for carries in AVX-512 mode
	// (half_arr_offset192 + RADIX) + 132 (=0x84) and round up to nearest multiple of 4:
/*	const int radix192_creals_in_local_store = (half_arr_offset192 + RADIX + 132 + 3)&(~0x3);	<*** GCC [not Clang] gives 'error: initializer element is not constant'... */
	const int radix192_creals_in_local_store = (0x304 + (RADIX>>3) + RADIX + 132 + 3)&(~0x3); // <*** ...here is that workaround for that stupidity.
  #elif defined(USE_AVX)
	const int half_arr_offset192 = 0x304 + (RADIX>>2);	// RADIX/4 vec_dbl slots for carries in AVX mode
	// (half_arr_offset192 + RADIX) + 132 (=0x84) and round up to nearest multiple of 4:
/*	const int radix192_creals_in_local_store = (half_arr_offset192 + RADIX + 132 + 3)&(~0x3);	<*** GCC [not Clang] gives 'error: initializer element is not constant'... */
	const int radix192_creals_in_local_store = (0x304 + (RADIX>>2) + RADIX + 132 + 3)&(~0x3); // <*** ...here is that workaround for that stupidity.
  #else
	const int half_arr_offset192 = 0x304 + (RADIX>>1);	// RADIX/2 vec_dbl slots for carries in SSE2 mode
	// (half_arr_offset192 + RADIX) = 40 and round up to nearest multiple of 4:
/*	const int radix192_creals_in_local_store = (half_arr_offset192 + RADIX + 40 + 3)&(~0x3); see above note in AVX decl. */
	const int radix192_creals_in_local_store = (0x304 + (RADIX>>1) + RADIX + 40 + 3)&(~0x3); // <*** ...here is that workaround for that stupidity.
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

/**************/

int radix192_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-192 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-192 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix192_ditN_cy_dif1";
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
	int j2,ntmp;
  #endif
  #ifndef MULTITHREAD
	int jp,jstart;
  #endif
  #ifdef USE_SSE2
	int nbytes;
  #endif
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 12|24|48 for avx512,avx,sse, respectively:
	int incr;
  #ifdef USE_AVX512
	const int incr_long = 12,incr_med = 6,incr_short = 4;
  #else
	const int incr_long = 12,incr_med = 8,incr_short = 4;
  #endif
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

	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0, nsave = 0;
  #ifndef MULTITHREAD
// Shared DIF+DIT:
	int ilo,ihi,idx,jdx,hi_neg,k0,k1,k2, *iptr;
	int nshift;
	static int plo[16],phi[12], t_offsets[64];
	uint64 i64;
// DIF:
	static int dif_o_offsets[RADIX];
	// Consts containing [in little-endian hex-char fashion] the 4 [A-D] 16-subperms defined in the output-ordering commentary below:
	const uint64 dif_16_operms[4] = {
		0x01235476ab98fecdull,	// [A]
		0x54762310fecd98baull,	// [B]
		0xab98fecd54762310ull,	// [C]
		0xfecd98ba23107645ull	// [D]
	};
// DIT:
	// To-be-inited array used to support circular (0-2)-element shifts of basic 3-perm p10-multiple patterns:
	static int p10_3perm[4][5] = { {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0} };
	static int dit_i_offsets[RADIX];
  #endif

	static int poff[RADIX>>2];	// Store mults of p4 offset for loop control
	static double radix_inv, n2inv;
#if defined(USE_SSE2) || !defined(MULTITHREAD)
	const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675;	/* sin(twopi/3)		*/
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
	double *add0,*add1,*add2,*add3;
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
	vec_dbl	// Non-static utility ptrs
  #ifndef MULTITHREAD
		*va0,*va1,*va2, *vb0,*vb1,*vb2,
		*vc0,*vc1,*vc2, *vd0,*vd1,*vd2,
  #endif
		*tmp,
  #ifndef MULTITHREAD
		*tm1,
  #endif
		*tm2;
	static vec_dbl *cc0,*ss0, *max_err, *sse2_rnd, *half_arr,
		*r00,
  #ifndef MULTITHREAD
		*r40,*r80,	// Head of RADIX*vec_cmplx-sized local store #1
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
	static task_control_t   task_control = {NULL, (void*)cy192_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	int m,m2;
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	int bjmodn[RADIX];
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
		if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {	// Only care about this divisibility property for LOACC carry modes:
			ASSERT(0 == ((RADIX/i) % incr),"Carry-chain wts-multipliers recurrence length must divide RADIX/[n-wayness of carry macro]!");
		}
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
		// consisting of radix192_creals_in_local_store dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix192_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix192_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = sc_ptr;	r00   = tmp;	// Head of RADIX*vec_cmplx-sized local store #1
	  #ifndef MULTITHREAD
		r40 = tmp + 0x080;
		r80 = tmp + 0x100;
	  #endif
		tmp += 0x180;					// Head of RADIX*vec_cmplx-sized local store #2
	  #ifndef MULTITHREAD
						s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
	  #endif
		tmp += 0x180;
		// DFT-64 roots all def'd locally in the DFT-64 functions, only need the radix-3 roots here:
		cc0	= tmp;
		ss0	= tmp + 1;
		tmp += 0x2;	// sc_ptr += 0x302
	  #ifdef USE_AVX512
		cy = tmp;		tmp += 0x18;	// RADIX/8 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;	// sc_ptr += 0x31c; This is where the value of half_arr_offset192 comes from
	  #elif defined(USE_AVX)
		cy = tmp;		tmp += 0x30;	// RADIX/4 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x334; This is where the value of half_arr_offset192 comes from
		half_arr= tmp + 0x02;	// This table needs 96*SZ_VD bytes in avx mode
	  #else
		cy = tmp;		tmp += 0x60;	// RADIX/2 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x364; This is where the value of half_arr_offset192 comes from
		half_arr= tmp + 0x02;	// This table needs 32*SZ_VD bytes in sse2 mode
	  #endif
//		ASSERT(half_arr_offset == (uint32)(half_arr-sc_ptr), "half_arr_offset mismatches actual!");
		ASSERT((radix192_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r00) + (20 << L2_SZ_VD), "radix192_creals_in_local_store checksum failed!");

		/* Roots of 1 for radix-3 DFTs: cc0 = (cc1+cc2+cc3)/3 - 1; subtract 1 from Nussbaumer's definition in order to ease in-place computation */
		VEC_DBL_INIT(cc0, c3m1);
		VEC_DBL_INIT(ss0, s   );
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
		nbytes = (intptr_t)cy - (intptr_t)cc0;	// #bytes in above sincos block of data
		tmp = cc0;
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
		NDIVR >>= 4;			pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
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

	#ifndef MULTITHREAD
	// Shared:
		plo[0x0] =  0; plo[0x1] = p1; plo[0x2] = p2; plo[0x3] = p3;
		plo[0x4] = p4; plo[0x5] = p5; plo[0x6] = p6; plo[0x7] = p7;
		plo[0x8] = p8; plo[0x9] = p9; plo[0xa] = pa; plo[0xb] = pb;
		plo[0xc] = pc; plo[0xd] = pd; plo[0xe] = pe; plo[0xf] = pf;

		phi[0x0] =   0; phi[0x1] = p10; phi[0x2] = p20; phi[0x3] = p30;
		phi[0x4] = p40; phi[0x5] = p50; phi[0x6] = p60; phi[0x7] = p70;
		phi[0x8] = p80; phi[0x9] = p90; phi[0xa] = pa0; phi[0xb] = pb0;

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

	// DIF:
	// Extract the 4 [A-D] 16-subperms defined in the output-ordering commentary below:
	/* Operm 1:
		0,1,2,3,5,4,7,6,a,b,9,8,f,e,c,d + p00  =  [A] + p00
		5,4,7,6,2,3,1,0,f,e,c,d,9,8,b,a + p10  =  [B] + p10
		a,b,9,8,f,e,c,d,5,4,7,6,2,3,1,0 + p20  =  [C] + p20
		f,e,c,d,9,8,b,a,2,3,1,0,7,6,4,5 + p30  =  [D] + p30
	*/
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
			dif_o_offsets[0x20+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dif_16_operms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x30+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
	/* Operm 2:
		5,4,7,6,2,3,1,0,f,e,c,d,9,8,b,a + p90  =  [B] + p10
		a,b,9,8,f,e,c,d,5,4,7,6,2,3,1,0 + p80  =  [C] + p00 (mod p40)
		f,e,c,d,9,8,b,a,2,3,1,0,7,6,4,5 + pb0  =  [D] + p30
		5,4,7,6,2,3,1,0,f,e,c,d,9,8,b,a + pa0  =  [B] + p20
	*/
		i64 = dif_16_operms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x40+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x50+l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dif_16_operms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x60+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dif_16_operms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x70+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
	/* Operm 3:
		a,b,9,8,f,e,c,d,5,4,7,6,2,3,1,0 + p60  =  [C] + p20
		f,e,c,d,9,8,b,a,2,3,1,0,7,6,4,5 + p70  =  [D] + p30 (mod p40)
		5,4,7,6,2,3,1,0,f,e,c,d,9,8,b,a + p50  =  [B] + p10
		a,b,9,8,f,e,c,d,5,4,7,6,2,3,1,0 + p40  =  [C] + p00
	*/
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x80+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dif_16_operms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x90+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dif_16_operms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0xa0+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0xb0+l] = plo[(i64 >> nshift)&0xf];
		}

	// DIT:
	/* Iperm 1:
		0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00
		f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10
		f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p30
		f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p20
	*/																		// swap [0x2*] and [0x3*] cols:
		dit_i_offsets[0x00] =  0;	dit_i_offsets[0x10] = pf+p10;	dit_i_offsets[0x30] = pf+p20;	dit_i_offsets[0x20] = pf+p30;
		dit_i_offsets[0x01] = p1;	dit_i_offsets[0x11] = pe+p10;	dit_i_offsets[0x31] = pe+p20;	dit_i_offsets[0x21] = pe+p30;
		dit_i_offsets[0x02] = p3;	dit_i_offsets[0x12] = pd+p10;	dit_i_offsets[0x32] = pd+p20;	dit_i_offsets[0x22] = pd+p30;
		dit_i_offsets[0x03] = p2;	dit_i_offsets[0x13] = pc+p10;	dit_i_offsets[0x33] = pc+p20;	dit_i_offsets[0x23] = pc+p30;
		dit_i_offsets[0x04] = p7;	dit_i_offsets[0x14] = pb+p10;	dit_i_offsets[0x34] = pb+p20;	dit_i_offsets[0x24] = pb+p30;
		dit_i_offsets[0x05] = p6;	dit_i_offsets[0x15] = pa+p10;	dit_i_offsets[0x35] = pa+p20;	dit_i_offsets[0x25] = pa+p30;
		dit_i_offsets[0x06] = p5;	dit_i_offsets[0x16] = p9+p10;	dit_i_offsets[0x36] = p9+p20;	dit_i_offsets[0x26] = p9+p30;
		dit_i_offsets[0x07] = p4;	dit_i_offsets[0x17] = p8+p10;	dit_i_offsets[0x37] = p8+p20;	dit_i_offsets[0x27] = p8+p30;
		dit_i_offsets[0x08] = pf;	dit_i_offsets[0x18] = p7+p10;	dit_i_offsets[0x38] = p7+p20;	dit_i_offsets[0x28] = p7+p30;
		dit_i_offsets[0x09] = pe;	dit_i_offsets[0x19] = p6+p10;	dit_i_offsets[0x39] = p6+p20;	dit_i_offsets[0x29] = p6+p30;
		dit_i_offsets[0x0a] = pd;	dit_i_offsets[0x1a] = p5+p10;	dit_i_offsets[0x3a] = p5+p20;	dit_i_offsets[0x2a] = p5+p30;
		dit_i_offsets[0x0b] = pc;	dit_i_offsets[0x1b] = p4+p10;	dit_i_offsets[0x3b] = p4+p20;	dit_i_offsets[0x2b] = p4+p30;
		dit_i_offsets[0x0c] = pb;	dit_i_offsets[0x1c] = p3+p10;	dit_i_offsets[0x3c] = p3+p20;	dit_i_offsets[0x2c] = p3+p30;
		dit_i_offsets[0x0d] = pa;	dit_i_offsets[0x1d] = p2+p10;	dit_i_offsets[0x3d] = p2+p20;	dit_i_offsets[0x2d] = p2+p30;
		dit_i_offsets[0x0e] = p9;	dit_i_offsets[0x1e] = p1+p10;	dit_i_offsets[0x3e] = p1+p20;	dit_i_offsets[0x2e] = p1+p30;
		dit_i_offsets[0x0f] = p8;	dit_i_offsets[0x1f] =  0+p10;	dit_i_offsets[0x3f] =  0+p20;	dit_i_offsets[0x2f] =  0+p30;
	/* Iperm 2:
		5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p10
		5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p00 (mod p40)
		5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p20
		9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p30
	*/			// swap [0x4*] and [0x5*] cols:
		dit_i_offsets[0x50] = p5;	dit_i_offsets[0x40] = p5+p10;	dit_i_offsets[0x60] = p5+p20;	dit_i_offsets[0x70] = p9+p30;
		dit_i_offsets[0x51] = p4;	dit_i_offsets[0x41] = p4+p10;	dit_i_offsets[0x61] = p4+p20;	dit_i_offsets[0x71] = p8+p30;
		dit_i_offsets[0x52] = p6;	dit_i_offsets[0x42] = p6+p10;	dit_i_offsets[0x62] = p6+p20;	dit_i_offsets[0x72] = pa+p30;
		dit_i_offsets[0x53] = p7;	dit_i_offsets[0x43] = p7+p10;	dit_i_offsets[0x63] = p7+p20;	dit_i_offsets[0x73] = pb+p30;
		dit_i_offsets[0x54] = p1;	dit_i_offsets[0x44] = p1+p10;	dit_i_offsets[0x64] = p1+p20;	dit_i_offsets[0x74] = pe+p30;
		dit_i_offsets[0x55] =  0;	dit_i_offsets[0x45] =  0+p10;	dit_i_offsets[0x65] =  0+p20;	dit_i_offsets[0x75] = pf+p30;
		dit_i_offsets[0x56] = p2;	dit_i_offsets[0x46] = p2+p10;	dit_i_offsets[0x66] = p2+p20;	dit_i_offsets[0x76] = pc+p30;
		dit_i_offsets[0x57] = p3;	dit_i_offsets[0x47] = p3+p10;	dit_i_offsets[0x67] = p3+p20;	dit_i_offsets[0x77] = pd+p30;
		dit_i_offsets[0x58] = p9;	dit_i_offsets[0x48] = p9+p10;	dit_i_offsets[0x68] = p9+p20;	dit_i_offsets[0x78] = p1+p30;
		dit_i_offsets[0x59] = p8;	dit_i_offsets[0x49] = p8+p10;	dit_i_offsets[0x69] = p8+p20;	dit_i_offsets[0x79] =  0+p30;
		dit_i_offsets[0x5a] = pa;	dit_i_offsets[0x4a] = pa+p10;	dit_i_offsets[0x6a] = pa+p20;	dit_i_offsets[0x7a] = p2+p30;
		dit_i_offsets[0x5b] = pb;	dit_i_offsets[0x4b] = pb+p10;	dit_i_offsets[0x6b] = pb+p20;	dit_i_offsets[0x7b] = p3+p30;
		dit_i_offsets[0x5c] = pe;	dit_i_offsets[0x4c] = pe+p10;	dit_i_offsets[0x6c] = pe+p20;	dit_i_offsets[0x7c] = p6+p30;
		dit_i_offsets[0x5d] = pf;	dit_i_offsets[0x4d] = pf+p10;	dit_i_offsets[0x6d] = pf+p20;	dit_i_offsets[0x7d] = p7+p30;
		dit_i_offsets[0x5e] = pc;	dit_i_offsets[0x4e] = pc+p10;	dit_i_offsets[0x6e] = pc+p20;	dit_i_offsets[0x7e] = p4+p30;
		dit_i_offsets[0x5f] = pd;	dit_i_offsets[0x4f] = pd+p10;	dit_i_offsets[0x6f] = pd+p20;	dit_i_offsets[0x7f] = p5+p30;
	/* Iperm 3:
		a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p20
		2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p30 (mod p40)
		a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p00
		2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p10
	*/		// swap [0x8*] <=> [0xa*] colpair, and [0x9*] <=> [0xb*] pair:
		dit_i_offsets[0xa0] = pa;	dit_i_offsets[0xb0] = p2+p10;	dit_i_offsets[0x80] = pa+p20;	dit_i_offsets[0x90] = p2+p30;
		dit_i_offsets[0xa1] = pb;	dit_i_offsets[0xb1] = p3+p10;	dit_i_offsets[0x81] = pb+p20;	dit_i_offsets[0x91] = p3+p30;
		dit_i_offsets[0xa2] = p8;	dit_i_offsets[0xb2] =  0+p10;	dit_i_offsets[0x82] = p8+p20;	dit_i_offsets[0x92] =  0+p30;
		dit_i_offsets[0xa3] = p9;	dit_i_offsets[0xb3] = p1+p10;	dit_i_offsets[0x83] = p9+p20;	dit_i_offsets[0x93] = p1+p30;
		dit_i_offsets[0xa4] = pc;	dit_i_offsets[0xb4] = p4+p10;	dit_i_offsets[0x84] = pc+p20;	dit_i_offsets[0x94] = p4+p30;
		dit_i_offsets[0xa5] = pd;	dit_i_offsets[0xb5] = p5+p10;	dit_i_offsets[0x85] = pd+p20;	dit_i_offsets[0x95] = p5+p30;
		dit_i_offsets[0xa6] = pf;	dit_i_offsets[0xb6] = p7+p10;	dit_i_offsets[0x86] = pf+p20;	dit_i_offsets[0x96] = p7+p30;
		dit_i_offsets[0xa7] = pe;	dit_i_offsets[0xb7] = p6+p10;	dit_i_offsets[0x87] = pe+p20;	dit_i_offsets[0x97] = p6+p30;
		dit_i_offsets[0xa8] = p2;	dit_i_offsets[0xb8] = pc+p10;	dit_i_offsets[0x88] = p2+p20;	dit_i_offsets[0x98] = pc+p30;
		dit_i_offsets[0xa9] = p3;	dit_i_offsets[0xb9] = pd+p10;	dit_i_offsets[0x89] = p3+p20;	dit_i_offsets[0x99] = pd+p30;
		dit_i_offsets[0xaa] =  0;	dit_i_offsets[0xba] = pf+p10;	dit_i_offsets[0x8a] =  0+p20;	dit_i_offsets[0x9a] = pf+p30;
		dit_i_offsets[0xab] = p1;	dit_i_offsets[0xbb] = pe+p10;	dit_i_offsets[0x8b] = p1+p20;	dit_i_offsets[0x9b] = pe+p30;
		dit_i_offsets[0xac] = p4;	dit_i_offsets[0xbc] = p8+p10;	dit_i_offsets[0x8c] = p4+p20;	dit_i_offsets[0x9c] = p8+p30;
		dit_i_offsets[0xad] = p5;	dit_i_offsets[0xbd] = p9+p10;	dit_i_offsets[0x8d] = p5+p20;	dit_i_offsets[0x9d] = p9+p30;
		dit_i_offsets[0xae] = p7;	dit_i_offsets[0xbe] = pb+p10;	dit_i_offsets[0x8e] = p7+p20;	dit_i_offsets[0x9e] = pb+p30;
		dit_i_offsets[0xaf] = p6;	dit_i_offsets[0xbf] = pa+p10;	dit_i_offsets[0x8f] = p6+p20;	dit_i_offsets[0x9f] = pa+p30;

	// Init elts of array used to support circular (0-2)-element shifts of the basic 3-perm p10-multiple patterns
	// [A] = 00,40,80, [B] = 10,50,90, [C] = 20,60,a0, [D] = 30,70,b0:
	  #ifndef USE_SSE2
		l = 0;	p10_3perm[0][l++] =   0; p10_3perm[0][l++] = p40; p10_3perm[0][l++] = p80; p10_3perm[0][l++] =   0; p10_3perm[0][l++] = p40;
		l = 0;	p10_3perm[1][l++] = p10; p10_3perm[1][l++] = p50; p10_3perm[1][l++] = p90; p10_3perm[1][l++] = p10; p10_3perm[1][l++] = p50;
		l = 0;	p10_3perm[2][l++] = p20; p10_3perm[2][l++] = p60; p10_3perm[2][l++] = pa0; p10_3perm[2][l++] = p20; p10_3perm[2][l++] = p60;
		l = 0;	p10_3perm[3][l++] = p30; p10_3perm[3][l++] = p70; p10_3perm[3][l++] = pb0; p10_3perm[3][l++] = p30; p10_3perm[3][l++] = p70;
	  #else	// Double the RHS immediately here to avoid awkward loop over 2D array:
		l = 0;	p10_3perm[0][l++] =   0; p10_3perm[0][l++] =0x80; p10_3perm[0][l++]=0x100; p10_3perm[0][l++] =   0; p10_3perm[0][l++] =0x80;
		l = 0;	p10_3perm[1][l++] =0x20; p10_3perm[1][l++] =0xa0; p10_3perm[1][l++]=0x120; p10_3perm[1][l++] =0x20; p10_3perm[1][l++] =0xa0;
		l = 0;	p10_3perm[2][l++] =0x40; p10_3perm[2][l++] =0xc0; p10_3perm[2][l++]=0x140; p10_3perm[2][l++] =0x40; p10_3perm[2][l++] =0xc0;
		l = 0;	p10_3perm[3][l++] =0x60; p10_3perm[3][l++] =0xe0; p10_3perm[3][l++]=0x160; p10_3perm[3][l++] =0x60; p10_3perm[3][l++] =0xe0;
	  #endif
	#endif // !MULTITHREAD

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

/*...The radix-192 final DIT pass is here.	*/

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
		#include "radix192_main_carry_loop.h"

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
		ASSERT(0x0 == cy192_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
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
	//	printf("Iter = %d, prp_mult = %4.1f, maxerr = %20.15f, a[0] = %20.15f\n",iter,prp_mult,maxerr,a[0]);
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

void radix192_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-192 complex DIF FFT pass on the data in the length-N real vector A.
*/
	int j,j1,j2,jp,jt,l,ilo,ihi,hi_neg,k0,k1,k2,nshift;
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0, first_entry=TRUE;
	static int plo[16],phi[12];
	const double c3m1= -1.50000000000000000000, s = 0.86602540378443864675;	// cos(twopi/3)-1, sin(twopi/3)
	static int i_offsets[64], o_offsets[RADIX];
	// Consts containing [in little-endian hex-char fashion] the 4 [A-D] 16-subperms defined in the output-ordering commentary below:
	const uint64 dif_16_operms[4] = {
		0x01235476ab98fecdull,	// [A]
		0x54762310fecd98baull,	// [B]
		0xab98fecd54762310ull,	// [C]
		0xfecd98ba23107645ull	// [D]
	};
	uint64 i64;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex t[RADIX], *tptr;
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
		NDIVR >>= 4;			pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );

		plo[0x0] =  0; plo[0x1] = p1; plo[0x2] = p2; plo[0x3] = p3;
		plo[0x4] = p4; plo[0x5] = p5; plo[0x6] = p6; plo[0x7] = p7;
		plo[0x8] = p8; plo[0x9] = p9; plo[0xa] = pa; plo[0xb] = pb;
		plo[0xc] = pc; plo[0xd] = pd; plo[0xe] = pe; plo[0xf] = pf;

		phi[0x0] =   0; phi[0x1] = p10; phi[0x2] = p20; phi[0x3] = p30;
		phi[0x4] = p40; phi[0x5] = p50; phi[0x6] = p60; phi[0x7] = p70;
		phi[0x8] = p80; phi[0x9] = p90; phi[0xa] = pa0; phi[0xb] = pb0;

		i_offsets[0x00] = 0x00<<1;	i_offsets[0x10] = 0x10<<1;	i_offsets[0x20] = 0x20<<1;	i_offsets[0x30] = 0x30<<1;
		i_offsets[0x01] = 0x01<<1;	i_offsets[0x11] = 0x11<<1;	i_offsets[0x21] = 0x21<<1;	i_offsets[0x31] = 0x31<<1;
		i_offsets[0x02] = 0x02<<1;	i_offsets[0x12] = 0x12<<1;	i_offsets[0x22] = 0x22<<1;	i_offsets[0x32] = 0x32<<1;
		i_offsets[0x03] = 0x03<<1;	i_offsets[0x13] = 0x13<<1;	i_offsets[0x23] = 0x23<<1;	i_offsets[0x33] = 0x33<<1;
		i_offsets[0x04] = 0x04<<1;	i_offsets[0x14] = 0x14<<1;	i_offsets[0x24] = 0x24<<1;	i_offsets[0x34] = 0x34<<1;
		i_offsets[0x05] = 0x05<<1;	i_offsets[0x15] = 0x15<<1;	i_offsets[0x25] = 0x25<<1;	i_offsets[0x35] = 0x35<<1;
		i_offsets[0x06] = 0x06<<1;	i_offsets[0x16] = 0x16<<1;	i_offsets[0x26] = 0x26<<1;	i_offsets[0x36] = 0x36<<1;
		i_offsets[0x07] = 0x07<<1;	i_offsets[0x17] = 0x17<<1;	i_offsets[0x27] = 0x27<<1;	i_offsets[0x37] = 0x37<<1;
		i_offsets[0x08] = 0x08<<1;	i_offsets[0x18] = 0x18<<1;	i_offsets[0x28] = 0x28<<1;	i_offsets[0x38] = 0x38<<1;
		i_offsets[0x09] = 0x09<<1;	i_offsets[0x19] = 0x19<<1;	i_offsets[0x29] = 0x29<<1;	i_offsets[0x39] = 0x39<<1;
		i_offsets[0x0a] = 0x0a<<1;	i_offsets[0x1a] = 0x1a<<1;	i_offsets[0x2a] = 0x2a<<1;	i_offsets[0x3a] = 0x3a<<1;
		i_offsets[0x0b] = 0x0b<<1;	i_offsets[0x1b] = 0x1b<<1;	i_offsets[0x2b] = 0x2b<<1;	i_offsets[0x3b] = 0x3b<<1;
		i_offsets[0x0c] = 0x0c<<1;	i_offsets[0x1c] = 0x1c<<1;	i_offsets[0x2c] = 0x2c<<1;	i_offsets[0x3c] = 0x3c<<1;
		i_offsets[0x0d] = 0x0d<<1;	i_offsets[0x1d] = 0x1d<<1;	i_offsets[0x2d] = 0x2d<<1;	i_offsets[0x3d] = 0x3d<<1;
		i_offsets[0x0e] = 0x0e<<1;	i_offsets[0x1e] = 0x1e<<1;	i_offsets[0x2e] = 0x2e<<1;	i_offsets[0x3e] = 0x3e<<1;
		i_offsets[0x0f] = 0x0f<<1;	i_offsets[0x1f] = 0x1f<<1;	i_offsets[0x2f] = 0x2f<<1;	i_offsets[0x3f] = 0x3f<<1;
	/*
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
	*/
	// Extract the 4 [A-D] 16-subperms defined in the output-ordering commentary below:
	/* Operm 1:
		0,1,2,3,5,4,7,6,a,b,9,8,f,e,c,d + p00  =  [A] + p00
		5,4,7,6,2,3,1,0,f,e,c,d,9,8,b,a + p10  =  [B] + p10
		a,b,9,8,f,e,c,d,5,4,7,6,2,3,1,0 + p20  =  [C] + p20
		f,e,c,d,9,8,b,a,2,3,1,0,7,6,4,5 + p30  =  [D] + p30
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
			o_offsets[0x20+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dif_16_operms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x30+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
	/* Operm 2:
		5,4,7,6,2,3,1,0,f,e,c,d,9,8,b,a + p90  =  [B] + p10
		a,b,9,8,f,e,c,d,5,4,7,6,2,3,1,0 + p80  =  [C] + p00 (mod p40)
		f,e,c,d,9,8,b,a,2,3,1,0,7,6,4,5 + pb0  =  [D] + p30
		5,4,7,6,2,3,1,0,f,e,c,d,9,8,b,a + pa0  =  [B] + p20
	*/
		i64 = dif_16_operms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x40+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x50+l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dif_16_operms[3];	// [D] subperm
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
		a,b,9,8,f,e,c,d,5,4,7,6,2,3,1,0 + p60  =  [C] + p20
		f,e,c,d,9,8,b,a,2,3,1,0,7,6,4,5 + p70  =  [D] + p30 (mod p40)
		5,4,7,6,2,3,1,0,f,e,c,d,9,8,b,a + p50  =  [B] + p10
		a,b,9,8,f,e,c,d,5,4,7,6,2,3,1,0 + p40  =  [C] + p00
	*/
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x80+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dif_16_operms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0x90+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dif_16_operms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0xa0+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			o_offsets[0xb0+l] = plo[(i64 >> nshift)&0xf];
		}
	}

/*...The radix-192 pass is here.	*/

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
	Twiddleless version arranges 64 sets of radix-3 DFT inputs as follows: 0 in upper left corner,
	decrement 64 horizontally and 3 vertically (mod 192):

	DIF/DIT input-scramble array = [
		00,80,40   00,80,40 + p0	90,50,10   90,50,10 + p0	60,20,a0   60,20,a0 + p0	30,b0,70   30,b0,70 + p0
		bd,7d,3d   b0,70,30 + pd	8d,4d,0d   80,40,00 + pd	5d,1d,9d   50,10,90 + pd	2d,ad,6d   20,a0,60 + pd
		ba,7a,3a   b0,70,30 + pa	8a,4a,0a   80,40,00 + pa	5a,1a,9a   50,10,90 + pa	2a,aa,6a   20,a0,60 + pa
		b7,77,37   b0,70,30 + p7	87,47,07   80,40,00 + p7	57,17,97   50,10,90 + p7	27,a7,67   20,a0,60 + p7
		b4,74,34   b0,70,30 + p4	84,44,04   80,40,00 + p4	54,14,94   50,10,90 + p4	24,a4,64   20,a0,60 + p4
		b1,71,31   b0,70,30 + p1	81,41,01   80,40,00 + p1	51,11,91   50,10,90 + p1	21,a1,61   20,a0,60 + p1
		ae,6e,2e   a0,60,20 + pe	7e,3e,be   70,30,b0 + pe	4e,0e,8e   40,00,80 + pe	1e,9e,5e   10,90,50 + pe
		ab,6b,2b   a0,60,20 + pb	7b,3b,bb   70,30,b0 + pb	4b,0b,8b   40,00,80 + pb	1b,9b,5b   10,90,50 + pb
		a8,68,28 = a0,60,20 + p8	78,38,b8 = 70,30,b0 + p8	48,08,88 = 40,00,80 + p8	18,98,58 = 10,90,50 + p8
		a5,65,25   a0,60,20 + p5	75,35,b5   70,30,b0 + p5	45,05,85   40,00,80 + p5	15,95,55   10,90,50 + p5
		a2,62,22   a0,60,20 + p2	72,32,b2   70,30,b0 + p2	42,02,82   40,00,80 + p2	12,92,52   10,90,50 + p2
		9f,5f,1f   90,50,10 + pf	6f,2f,af   60,20,a0 + pf	3f,bf,7f   30,b0,70 + pf	0f,8f,4f   00,80,40 + pf
		9c,5c,1c   90,50,10 + pc	6c,2c,ac   60,20,a0 + pc	3c,bc,7c   30,b0,70 + pc	0c,8c,4c   00,80,40 + pc
		99,59,19   90,50,10 + p9	69,29,a9   60,20,a0 + p9	39,b9,79   30,b0,70 + p9	09,89,49   00,80,40 + p9
		96,56,16   90,50,10 + p6	66,26,a6   60,20,a0 + p6	36,b6,76   30,b0,70 + p6	06,86,46   00,80,40 + p6
		93,53,13   90,50,10 + p3	63,23,a3   60,20,a0 + p3	33,b3,73   30,b0,70 + p3	03,83,43   00,80,40 + p3
		[cont. in col 2]			[cont. in col 3]			[cont. in col 4]
	If we prestore the length-3 perms appearing left of the + above we get 12 such perms added to the p0-f offsets,
	which are really just [0-2]-element leftward-cshifts of the 4 base-perms

		[0] = 00,80,40;		[1] = b0,70,30;		[2] = a0,60,20;		[3] = 90,50,10.

	Thus in terms of these base-perms and the lcshift count:

		[0]<<0 + p0	[3]<<0 + p0	[2]<<1 + p0	[1]<<2 + p0
		[1]<<0 + pd	[0]<<1 + pd	[3]<<1 + pd	[2]<<2 + pd
		[1]<<0 + pa	[0]<<1 + pa	[3]<<1 + pa	[2]<<2 + pa
		[1]<<0 + p7	[0]<<1 + p7	[3]<<1 + p7	[2]<<2 + p7
		[1]<<0 + p4	[0]<<1 + p4	[3]<<1 + p4	[2]<<2 + p4
		[1]<<0 + p1	[0]<<1 + p1	[3]<<1 + p1	[2]<<2 + p1
		[2]<<0 + pe	[1]<<1 + pe	[0]<<2 + pe	[1]<<2 + pe
		[2]<<0 + pb	[1]<<1 + pb	[0]<<2 + pb	[1]<<2 + pb
		[2]<<0 + p8	[1]<<1 + p8	[0]<<2 + p8	[1]<<2 + p8
		[2]<<0 + p5	[1]<<1 + p5	[0]<<2 + p5	[1]<<2 + p5
		[2]<<0 + p2	[1]<<1 + p2	[0]<<2 + p2	[1]<<2 + p2
		[3]<<0 + pf	[2]<<1 + pf	[1]<<2 + pf	[0]<<0 + pf
		[3]<<0 + pc	[2]<<1 + pc	[1]<<2 + pc	[0]<<0 + pc
		[3]<<0 + p9	[2]<<1 + p9	[1]<<2 + p9	[0]<<0 + p9
		[3]<<0 + p6	[2]<<1 + p6	[1]<<2 + p6	[0]<<0 + p6
		[3]<<0 + p3	[2]<<1 + p3	[1]<<2 + p3	[0]<<0 + p3
		[cont.2]	[cont.3]	[cont.4]

	Thus specification of each complete 3-perm needs:

		2 bits for the base-perm [0-3];
		2 bits for the base-perm lcshift count [0-2];
		4 bits for the p0-f addend in terms of plo array index.

	Thus 64 x 8 bits for the complete set of 64 3-perms.
	Or ... just do some fancy-indexing footwork to directly implement the horiz/vert index decrementing, as below:
	*/
	//...gather the needed data (192 64-bit complex) and do 16 radix-3 transforms - We want unit-strides in the radix32-DFT macros, so use large output strides here:
		tptr = t;
		ilo = ihi = 0;	// ints tracking the low and high nibbles of the p* offset at head  of each trio
		k0 =   0;	k1 = p80;	k2 = p40;
		for(l = 0; l < 64; l++) {
			jp = plo[ilo];	// Same for each elt of the current-offset trio
			jt = j1 + jp; jp += j2;
			RADIX_03_DFT(s,c3m1,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],
				t00,t01,t02,t03,t04,t05,
				tptr->re,tptr->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x80)->re,(tptr+0x80)->im
			);	tptr++;
			ilo -= 3;
			// If low nibble underflows need to add 16 to it, decrement high nibble ihi, and recompute k1,k2.
			// Use branch here because dependent-op count > 10:
			if(ilo < 0) {
				ilo += 16;
				ihi--;			hi_neg = ihi< 0;	ihi+= ((-hi_neg) & 12);	// 1st high-part offset
				k1 = ihi- 4;	hi_neg = k1 < 0;	k1 += ((-hi_neg) & 12);	// high part of "decrement 64 horizontally" to get 2nd offset
				k2 = k1 - 4;	hi_neg = k2 < 0;	k2 += ((-hi_neg) & 12);	// high part of "decrement 64 horizontally" to get 3rd offset
				k0 = phi[ihi];	k1 = phi[k1];	k2 = phi[k2];	// Use distinct ihi and k0 = phi[ihi] for 1st high-part offset because basic-index ihi must survive from one loop pass to next
			}
		}

	/*...and now do 3 radix-64 transforms.
	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF]
	to properly permute the outputs to the required ordering, which in terms of our 3 radix-64 DFTs is

		00,01,02,03,05,04,07,06,0a,0b,09,08,0f,0e,0c,0d  =  0,1,2,3,5,4,7,6,a,b,9,8,f,e,c,d + p00  =  [A] + p00
		15,14,17,16,12,13,11,10,1f,1e,1c,1d,19,18,1b,1a  =  5,4,7,6,2,3,1,0,f,e,c,d,9,8,b,a + p10  =  [B] + p10
		2a,2b,29,28,2f,2e,2c,2d,25,24,27,26,22,23,21,20  =  a,b,9,8,f,e,c,d,5,4,7,6,2,3,1,0 + p20  =  [C] + p20
		3f,3e,3c,3d,39,38,3b,3a,32,33,31,30,37,36,34,35  =  f,e,c,d,9,8,b,a,2,3,1,0,7,6,4,5 + p30  =  [D] + p30

		95,94,97,96,92,93,91,90,9f,9e,9c,9d,99,98,9b,9a  =  5,4,7,6,2,3,1,0,f,e,c,d,9,8,b,a + p90  =  [B] + p90
		8a,8b,89,88,8f,8e,8c,8d,85,84,87,86,82,83,81,80  =  a,b,9,8,f,e,c,d,5,4,7,6,2,3,1,0 + p80  =  [C] + p80
		bf,be,bc,bd,b9,b8,bb,ba,b2,b3,b1,b0,b7,b6,b4,b5  =  f,e,c,d,9,8,b,a,2,3,1,0,7,6,4,5 + pb0  =  [D] + pb0
		a5,a4,a7,a6,a2,a3,a1,a0,af,ae,ac,ad,a9,a8,ab,aa  =  5,4,7,6,2,3,1,0,f,e,c,d,9,8,b,a + pa0  =  [B] + pa0

		6a,6b,69,68,6f,6e,6c,6d,65,64,67,66,62,63,61,60  =  a,b,9,8,f,e,c,d,5,4,7,6,2,3,1,0 + p60  =  [C] + p60
		7f,7e,7c,7d,79,78,7b,7a,72,73,71,70,77,76,74,75  =  f,e,c,d,9,8,b,a,2,3,1,0,7,6,4,5 + p70  =  [D] + p70
		55,54,57,56,52,53,51,50,5f,5e,5c,5d,59,58,5b,5a  =  5,4,7,6,2,3,1,0,f,e,c,d,9,8,b,a + p50  =  [B] + p50
		4a,4b,49,48,4f,4e,4c,4d,45,44,47,46,42,43,41,40  =  a,b,9,8,f,e,c,d,5,4,7,6,2,3,1,0 + p40  =  [C] + p40

	Thus we can encode the needed o_offset 64-perms in (mod p40) fashion in terms of the above 4 [A-D] 16-subperms as

		Block 1: a[]		Block 2: a[]+p80		Block 3: a[]+p40
		[A] + p00			[B] + p10				[C] + p20
		[B] + p10			[C] + p00				[D] + p30
		[C] + p20			[D] + p30				[B] + p10
		[D] + p30			[B] + p20				[C] + p00
	*/
		//	NOTE that RADIX_64_DIF outputs are IN-ORDER rather than BR:
		jt = j1    ;	RADIX_64_DIF((double *)(t+0x00),i_offsets,1, (a+jt),(o_offsets+0x00),RE_IM_STRIDE);	// Inputs in t[00-63]
		jt = j1+p80;	RADIX_64_DIF((double *)(t+0x40),i_offsets,1, (a+jt),(o_offsets+0x40),RE_IM_STRIDE);	// Inputs in t[64-127]
		jt = j1+p40;	RADIX_64_DIF((double *)(t+0x80),i_offsets,1, (a+jt),(o_offsets+0x80),RE_IM_STRIDE);	// Inputs in t[128-191]
	}
}

/***************/

void radix192_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-192 complex DIT FFT pass on the data in the length-N real vector A.
*/
	int j,j1,j2,jp,jt,l,ilo,ihi,idx,jdx,/* lo_neg,hi_neg, */k0,k1,k2, *iptr;
	// In order to preserve length-2 numeric-index-offset property here, use hex for digit:
	static int NDIVR,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0, first_entry=TRUE;
	static int plo[16];
	//static int phi[12];
	// To-be-inited array used to support circular (0-2)-element shifts of basic 3-perm p10-multiple patterns:
	static int p10_3perm[4][5] = { {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0} };
	const double c3m1= -1.50000000000000000000, s = 0.86602540378443864675;	// cos(twopi/3)-1, sin(twopi/3)
	static int o_offsets[64], i_offsets[RADIX];
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex t[RADIX], *tptr;
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
		NDIVR >>= 4;			pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );

		plo[0x0] =  0; plo[0x1] = p1; plo[0x2] = p2; plo[0x3] = p3;
		plo[0x4] = p4; plo[0x5] = p5; plo[0x6] = p6; plo[0x7] = p7;
		plo[0x8] = p8; plo[0x9] = p9; plo[0xa] = pa; plo[0xb] = pb;
		plo[0xc] = pc; plo[0xd] = pd; plo[0xe] = pe; plo[0xf] = pf;

	#if 0
		phi[0x0] =   0; phi[0x1] = p10; phi[0x2] = p20; phi[0x3] = p30;
		phi[0x4] = p40; phi[0x5] = p50; phi[0x6] = p60; phi[0x7] = p70;
		phi[0x8] = p80; phi[0x9] = p90; phi[0xa] = pa0; phi[0xb] = pb0;
	#endif

	/* Iperm 1:
		0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00
		f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10
		f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p30
		f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p20
	*/															// swap [0x2*] and [0x3*] cols:
		i_offsets[0x00] =  0;	i_offsets[0x10] = pf+p10;	i_offsets[0x30] = pf+p20;	i_offsets[0x20] = pf+p30;
		i_offsets[0x01] = p1;	i_offsets[0x11] = pe+p10;	i_offsets[0x31] = pe+p20;	i_offsets[0x21] = pe+p30;
		i_offsets[0x02] = p3;	i_offsets[0x12] = pd+p10;	i_offsets[0x32] = pd+p20;	i_offsets[0x22] = pd+p30;
		i_offsets[0x03] = p2;	i_offsets[0x13] = pc+p10;	i_offsets[0x33] = pc+p20;	i_offsets[0x23] = pc+p30;
		i_offsets[0x04] = p7;	i_offsets[0x14] = pb+p10;	i_offsets[0x34] = pb+p20;	i_offsets[0x24] = pb+p30;
		i_offsets[0x05] = p6;	i_offsets[0x15] = pa+p10;	i_offsets[0x35] = pa+p20;	i_offsets[0x25] = pa+p30;
		i_offsets[0x06] = p5;	i_offsets[0x16] = p9+p10;	i_offsets[0x36] = p9+p20;	i_offsets[0x26] = p9+p30;
		i_offsets[0x07] = p4;	i_offsets[0x17] = p8+p10;	i_offsets[0x37] = p8+p20;	i_offsets[0x27] = p8+p30;
		i_offsets[0x08] = pf;	i_offsets[0x18] = p7+p10;	i_offsets[0x38] = p7+p20;	i_offsets[0x28] = p7+p30;
		i_offsets[0x09] = pe;	i_offsets[0x19] = p6+p10;	i_offsets[0x39] = p6+p20;	i_offsets[0x29] = p6+p30;
		i_offsets[0x0a] = pd;	i_offsets[0x1a] = p5+p10;	i_offsets[0x3a] = p5+p20;	i_offsets[0x2a] = p5+p30;
		i_offsets[0x0b] = pc;	i_offsets[0x1b] = p4+p10;	i_offsets[0x3b] = p4+p20;	i_offsets[0x2b] = p4+p30;
		i_offsets[0x0c] = pb;	i_offsets[0x1c] = p3+p10;	i_offsets[0x3c] = p3+p20;	i_offsets[0x2c] = p3+p30;
		i_offsets[0x0d] = pa;	i_offsets[0x1d] = p2+p10;	i_offsets[0x3d] = p2+p20;	i_offsets[0x2d] = p2+p30;
		i_offsets[0x0e] = p9;	i_offsets[0x1e] = p1+p10;	i_offsets[0x3e] = p1+p20;	i_offsets[0x2e] = p1+p30;
		i_offsets[0x0f] = p8;	i_offsets[0x1f] =  0+p10;	i_offsets[0x3f] =  0+p20;	i_offsets[0x2f] =  0+p30;
	/* Iperm 2:
		5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p10
		5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p00 (mod p40)
		5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p20
		9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p30
	*/		// swap [0x4*] and [0x5*] cols:
		i_offsets[0x50] = p5;	i_offsets[0x40] = p5+p10;	i_offsets[0x60] = p5+p20;	i_offsets[0x70] = p9+p30;
		i_offsets[0x51] = p4;	i_offsets[0x41] = p4+p10;	i_offsets[0x61] = p4+p20;	i_offsets[0x71] = p8+p30;
		i_offsets[0x52] = p6;	i_offsets[0x42] = p6+p10;	i_offsets[0x62] = p6+p20;	i_offsets[0x72] = pa+p30;
		i_offsets[0x53] = p7;	i_offsets[0x43] = p7+p10;	i_offsets[0x63] = p7+p20;	i_offsets[0x73] = pb+p30;
		i_offsets[0x54] = p1;	i_offsets[0x44] = p1+p10;	i_offsets[0x64] = p1+p20;	i_offsets[0x74] = pe+p30;
		i_offsets[0x55] =  0;	i_offsets[0x45] =  0+p10;	i_offsets[0x65] =  0+p20;	i_offsets[0x75] = pf+p30;
		i_offsets[0x56] = p2;	i_offsets[0x46] = p2+p10;	i_offsets[0x66] = p2+p20;	i_offsets[0x76] = pc+p30;
		i_offsets[0x57] = p3;	i_offsets[0x47] = p3+p10;	i_offsets[0x67] = p3+p20;	i_offsets[0x77] = pd+p30;
		i_offsets[0x58] = p9;	i_offsets[0x48] = p9+p10;	i_offsets[0x68] = p9+p20;	i_offsets[0x78] = p1+p30;
		i_offsets[0x59] = p8;	i_offsets[0x49] = p8+p10;	i_offsets[0x69] = p8+p20;	i_offsets[0x79] =  0+p30;
		i_offsets[0x5a] = pa;	i_offsets[0x4a] = pa+p10;	i_offsets[0x6a] = pa+p20;	i_offsets[0x7a] = p2+p30;
		i_offsets[0x5b] = pb;	i_offsets[0x4b] = pb+p10;	i_offsets[0x6b] = pb+p20;	i_offsets[0x7b] = p3+p30;
		i_offsets[0x5c] = pe;	i_offsets[0x4c] = pe+p10;	i_offsets[0x6c] = pe+p20;	i_offsets[0x7c] = p6+p30;
		i_offsets[0x5d] = pf;	i_offsets[0x4d] = pf+p10;	i_offsets[0x6d] = pf+p20;	i_offsets[0x7d] = p7+p30;
		i_offsets[0x5e] = pc;	i_offsets[0x4e] = pc+p10;	i_offsets[0x6e] = pc+p20;	i_offsets[0x7e] = p4+p30;
		i_offsets[0x5f] = pd;	i_offsets[0x4f] = pd+p10;	i_offsets[0x6f] = pd+p20;	i_offsets[0x7f] = p5+p30;
	/* Iperm 3:
		a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p20
		2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p30 (mod p40)
		a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p00
		2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p10
	*/		// swap [0x8*] <=> [0xa*] colpair, and [0x9*] <=> [0xb*] pair:
		i_offsets[0xa0] = pa;	i_offsets[0xb0] = p2+p10;	i_offsets[0x80] = pa+p20;	i_offsets[0x90] = p2+p30;
		i_offsets[0xa1] = pb;	i_offsets[0xb1] = p3+p10;	i_offsets[0x81] = pb+p20;	i_offsets[0x91] = p3+p30;
		i_offsets[0xa2] = p8;	i_offsets[0xb2] =  0+p10;	i_offsets[0x82] = p8+p20;	i_offsets[0x92] =  0+p30;
		i_offsets[0xa3] = p9;	i_offsets[0xb3] = p1+p10;	i_offsets[0x83] = p9+p20;	i_offsets[0x93] = p1+p30;
		i_offsets[0xa4] = pc;	i_offsets[0xb4] = p4+p10;	i_offsets[0x84] = pc+p20;	i_offsets[0x94] = p4+p30;
		i_offsets[0xa5] = pd;	i_offsets[0xb5] = p5+p10;	i_offsets[0x85] = pd+p20;	i_offsets[0x95] = p5+p30;
		i_offsets[0xa6] = pf;	i_offsets[0xb6] = p7+p10;	i_offsets[0x86] = pf+p20;	i_offsets[0x96] = p7+p30;
		i_offsets[0xa7] = pe;	i_offsets[0xb7] = p6+p10;	i_offsets[0x87] = pe+p20;	i_offsets[0x97] = p6+p30;
		i_offsets[0xa8] = p2;	i_offsets[0xb8] = pc+p10;	i_offsets[0x88] = p2+p20;	i_offsets[0x98] = pc+p30;
		i_offsets[0xa9] = p3;	i_offsets[0xb9] = pd+p10;	i_offsets[0x89] = p3+p20;	i_offsets[0x99] = pd+p30;
		i_offsets[0xaa] =  0;	i_offsets[0xba] = pf+p10;	i_offsets[0x8a] =  0+p20;	i_offsets[0x9a] = pf+p30;
		i_offsets[0xab] = p1;	i_offsets[0xbb] = pe+p10;	i_offsets[0x8b] = p1+p20;	i_offsets[0x9b] = pe+p30;
		i_offsets[0xac] = p4;	i_offsets[0xbc] = p8+p10;	i_offsets[0x8c] = p4+p20;	i_offsets[0x9c] = p8+p30;
		i_offsets[0xad] = p5;	i_offsets[0xbd] = p9+p10;	i_offsets[0x8d] = p5+p20;	i_offsets[0x9d] = p9+p30;
		i_offsets[0xae] = p7;	i_offsets[0xbe] = pb+p10;	i_offsets[0x8e] = p7+p20;	i_offsets[0x9e] = pb+p30;
		i_offsets[0xaf] = p6;	i_offsets[0xbf] = pa+p10;	i_offsets[0x8f] = p6+p20;	i_offsets[0x9f] = pa+p30;

		o_offsets[0x00] = 0x00<<1;	o_offsets[0x10] = 0x10<<1;	o_offsets[0x20] = 0x20<<1;	o_offsets[0x30] = 0x30<<1;
		o_offsets[0x01] = 0x01<<1;	o_offsets[0x11] = 0x11<<1;	o_offsets[0x21] = 0x21<<1;	o_offsets[0x31] = 0x31<<1;
		o_offsets[0x02] = 0x02<<1;	o_offsets[0x12] = 0x12<<1;	o_offsets[0x22] = 0x22<<1;	o_offsets[0x32] = 0x32<<1;
		o_offsets[0x03] = 0x03<<1;	o_offsets[0x13] = 0x13<<1;	o_offsets[0x23] = 0x23<<1;	o_offsets[0x33] = 0x33<<1;
		o_offsets[0x04] = 0x04<<1;	o_offsets[0x14] = 0x14<<1;	o_offsets[0x24] = 0x24<<1;	o_offsets[0x34] = 0x34<<1;
		o_offsets[0x05] = 0x05<<1;	o_offsets[0x15] = 0x15<<1;	o_offsets[0x25] = 0x25<<1;	o_offsets[0x35] = 0x35<<1;
		o_offsets[0x06] = 0x06<<1;	o_offsets[0x16] = 0x16<<1;	o_offsets[0x26] = 0x26<<1;	o_offsets[0x36] = 0x36<<1;
		o_offsets[0x07] = 0x07<<1;	o_offsets[0x17] = 0x17<<1;	o_offsets[0x27] = 0x27<<1;	o_offsets[0x37] = 0x37<<1;
		o_offsets[0x08] = 0x08<<1;	o_offsets[0x18] = 0x18<<1;	o_offsets[0x28] = 0x28<<1;	o_offsets[0x38] = 0x38<<1;
		o_offsets[0x09] = 0x09<<1;	o_offsets[0x19] = 0x19<<1;	o_offsets[0x29] = 0x29<<1;	o_offsets[0x39] = 0x39<<1;
		o_offsets[0x0a] = 0x0a<<1;	o_offsets[0x1a] = 0x1a<<1;	o_offsets[0x2a] = 0x2a<<1;	o_offsets[0x3a] = 0x3a<<1;
		o_offsets[0x0b] = 0x0b<<1;	o_offsets[0x1b] = 0x1b<<1;	o_offsets[0x2b] = 0x2b<<1;	o_offsets[0x3b] = 0x3b<<1;
		o_offsets[0x0c] = 0x0c<<1;	o_offsets[0x1c] = 0x1c<<1;	o_offsets[0x2c] = 0x2c<<1;	o_offsets[0x3c] = 0x3c<<1;
		o_offsets[0x0d] = 0x0d<<1;	o_offsets[0x1d] = 0x1d<<1;	o_offsets[0x2d] = 0x2d<<1;	o_offsets[0x3d] = 0x3d<<1;
		o_offsets[0x0e] = 0x0e<<1;	o_offsets[0x1e] = 0x1e<<1;	o_offsets[0x2e] = 0x2e<<1;	o_offsets[0x3e] = 0x3e<<1;
		o_offsets[0x0f] = 0x0f<<1;	o_offsets[0x1f] = 0x1f<<1;	o_offsets[0x2f] = 0x2f<<1;	o_offsets[0x3f] = 0x3f<<1;

	// Init elts of array used to support circular (0-2)-element shifts of the basic 3-perm p10-multiple patterns
	// [A] = 00,40,80, [B] = 10,50,90, [C] = 20,60,a0, [D] = 30,70,b0:
		l = 0;
		p10_3perm[0][l++] =   0; p10_3perm[0][l++] = p40; p10_3perm[0][l++] = p80; p10_3perm[0][l++] =   0; p10_3perm[0][l++] = p40;
		l = 0;
		p10_3perm[1][l++] = p10; p10_3perm[1][l++] = p50; p10_3perm[1][l++] = p90; p10_3perm[1][l++] = p10; p10_3perm[1][l++] = p50;
		l = 0;
		p10_3perm[2][l++] = p20; p10_3perm[2][l++] = p60; p10_3perm[2][l++] = pa0; p10_3perm[2][l++] = p20; p10_3perm[2][l++] = p60;
		l = 0;
		p10_3perm[3][l++] = p30; p10_3perm[3][l++] = p70; p10_3perm[3][l++] = pb0; p10_3perm[3][l++] = p30; p10_3perm[3][l++] = p70;
	}

/*...The radix-192 pass is here.	*/

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
	Remember, inputs to DIT are bit-reversed, so using output of test_fft_radix(), we have

	Combined DIT input-scramble array = [
		00,01,03,02,07,06,05,04,0f,0e,0d,0c,0b,0a,09,08  =  0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00  =  [A] + p00
		1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10  =  f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10  =  [B] + p10
		3f,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30  =  f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p30  =  [B] + p30
		2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20  =  f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p20  =  [B] + p20

		95,94,96,97,91,90,92,93,99,98,9a,9b,9e,9f,9c,9d  =  5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p90  =  [C] + p90
		85,84,86,87,81,80,82,83,89,88,8a,8b,8e,8f,8c,8d  =  5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p80  =  [C] + p80
		a5,a4,a6,a7,a1,a0,a2,a3,a9,a8,aa,ab,ae,af,ac,ad  =  5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + pa0  =  [C] + pa0
		b9,b8,ba,bb,be,bf,bc,bd,b1,b0,b2,b3,b6,b7,b4,b5  =  9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + pb0  =  [D] + pb0

		6a,6b,68,69,6c,6d,6f,6e,62,63,60,61,64,65,67,66  =  a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p60  =  [E] + p60
		72,73,70,71,74,75,77,76,7c,7d,7f,7e,78,79,7b,7a  =  2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p70  =  [F] + p70
		4a,4b,48,49,4c,4d,4f,4e,42,43,40,41,44,45,47,46  =  a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p40  =  [E] + p40
		52,53,50,51,54,55,57,56,5c,5d,5f,5e,58,59,5b,5a  =  2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p50  =  [F] + p50
	*/
	//...gather the needed data (192 64-bit complex) and do 3 radix-64 transforms:
		jt = j1    ;	RADIX_64_DIT(a+jt,(i_offsets+0x00),RE_IM_STRIDE, (double *)(t+0x00),o_offsets,1);	// Outputs in t[00-63]
		jt = j1+p80;	RADIX_64_DIT(a+jt,(i_offsets+0x40),RE_IM_STRIDE, (double *)(t+0x40),o_offsets,1);	// Outputs in t[64-127]
		jt = j1+p40;	RADIX_64_DIT(a+jt,(i_offsets+0x80),RE_IM_STRIDE, (double *)(t+0x80),o_offsets,1);	// Outputs in t[64-127]

	//...and now do 64 radix-3 transforms:
	/*
	Use the inverse supercalafragalistic Ancient Chinese Secret index-munging formula [iSACSIMPF]
	to properly permute the outputs of the 64 radix-3 DFTs.
	Since our first-look oindex ordering was +p0x[0,40,80] for each radix-3 and incrementing += p1 between those DFTs,
	arrange resulting mismatched-data-sorted index permutation into 3 vertical 64-entry columns to get needed oindex patterns.

	Note:
	o the lo part of the p-offset [i.e. index into plo[]] results from init = 0 and decrementing 1 (mod 16) each loop pass;
		[These occurrences are marked with --...*] Note that based on the above data, when ilo-- underflows, we don't meed
		to check if ihi-- - 8 also underflows, but see no obvious way to take advantage of this.
	o the hi part of the p-offset [i.e. index into phi[]] results from init = 0 and decrementing 8 (mod 12) each loop pass.

	Indexing in hex for clarity and using [A-D]0-2 notation in the rightmost column to flag reusable circular
	(0-2)-element shifts of the basic patterns [A] = 00,40,80, [B] = 10,50,90, [C] = 20,60,a0, [D] = 30,70,b0:

		00,40,80		00,40,80 + p0		A0	<*** This is the only [A-D] shift index which is "out of place",
		3f,7f,bf	*--	30,70,b0 + pf*		D0		i.e. does not fit into the increment-mod-3 pattern of the rest.
		7e,be,3e		70,b0,30 + pe		D1
		bd,3d,7d		b0,30,70 + pd		D2
		3c,7c,bc		30,70,b0 + pc		D0
		7b,bb,3b		70,b0,30 + pb		D1
		ba,3a,7a		b0,30,70 + pa		D2
		39,79,b9		30,70,b0 + p9		D0
		78,b8,38		70,b0,30 + p8		D1
		b7,37,77		b0,30,70 + p7		D2
		36,76,b6		30,70,b0 + p6		D0
		75,b5,35		70,b0,30 + p5		D1
		b4,34,74		b0,30,70 + p4		D2
		33,73,b3		30,70,b0 + p3		D0
		72,b2,32		70,b0,30 + p2		D1
		b1,31,71		b0,30,70 + p1		D2
		30,70,b0		30,70,b0 + p0		D0
		6f,af,2f	*--	60,a0,20 + pf*		C1	<*** Whenever lo index underflows, decr. the A-D subpattern idx (mod 4)
		ae,2e,6e		a0,20,60 + pe		C2
		2d,6d,ad		20,60,a0 + pd		C0
		6c,ac,2c		60,a0,20 + pc		C1
		ab,2b,6b		a0,20,60 + pb		C2
		2a,6a,aa		20,60,a0 + pa		C0
		69,a9,29		60,a0,20 + p9		C1
		a8,28,68		a0,20,60 + p8		C2
		27,67,a7	=	20,60,a0 + p7		C0
		66,a6,26		60,a0,20 + p6		C1
		a5,25,65		a0,20,60 + p5		C2
		24,64,a4		20,60,a0 + p4		C0
		63,a3,23		60,a0,20 + p3		C1
		a2,22,62		a0,20,60 + p2		C2
		21,61,a1		20,60,a0 + p1		C0
		60,a0,20		60,a0,20 + p0		C1
		9f,1f,5f	*--	90,10,50 + pf*		B2
		1e,5e,9e		10,50,90 + pe		B0
		5d,9d,1d		50,90,10 + pd		B1
		9c,1c,5c		90,10,50 + pc		B2
		1b,5b,9b		10,50,90 + pb		B0
		5a,9a,1a		50,90,10 + pa		B1
		99,19,59		90,10,50 + p9		B2
		18,58,98		10,50,90 + p8		B0
		57,97,17		50,90,10 + p7		B1
		96,16,56		90,10,50 + p6		B2
		15,55,95		10,50,90 + p5		B0
		54,94,14		50,90,10 + p4		B1
		93,13,53		90,10,50 + p3		B2
		12,52,92		10,50,90 + p2		B0
		51,91,11		50,90,10 + p1		B1
		90,10,50		90,10,50 + p0		B2
		0f,4f,8f	*--	00,40,80 + pf*		A0
		4e,8e,0e		40,80,00 + pe		A1
		8d,0d,4d		80,00,40 + pd		A2
		0c,4c,8c		00,40,80 + pc		A0
		4b,8b,0b		40,80,00 + pb		A1
		8a,0a,4a		80,00,40 + pa		A2
		09,49,89		00,40,80 + p9		A0
		48,88,08		40,80,00 + p8		A1
		87,07,47		80,00,40 + p7		A2
		06,46,86		00,40,80 + p6		A0
		45,85,05		40,80,00 + p5		A1
		84,04,44		80,00,40 + p4		A2
		03,43,83		00,40,80 + p3		A0
		42,82,02		40,80,00 + p2		A1
		81,01,41		80,00,40 + p1		A2

	Can generate idx 0,1,2,... pattern on-the-fly sans explicit (mod 3)ing by initing idx = 2, but setting the k0-2
	offset values corr. to A0 [rather than Aidx] on loop entry to take care of the leading A0, then at end of each
	loop pass doing
					idx++;
					idx += (idx == 3);	// When idx hits 3, bump up to 4...
					idx &= 3;	// ...and reduce (mod 4), which is cheaper than (mod 3)
	*/
		tptr = t;
		ilo = ihi = 0;	// ints tracking the low and high nibbles of the p* offset at head  of each trio
		idx = 2; jdx = 0; iptr = &p10_3perm[jdx][0];	// jdx is [A-D] index, idx is 0-2, so this formally corr. to A2...
		// ... but set the k0-2 offset values corr. to A0 [rather than A2] on loop entry to take care of leading A0:
		k0 =   0;	k1 = p40;	k2 = p80;
		for(l = 0; l < 64; l++) {
			jp = plo[ilo];	// Same for each elt of the current-offset trio
			jt = j1 + jp; jp += j2;
			RADIX_03_DFT(s,c3m1,
				tptr->re,tptr->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x80)->re,(tptr+0x80)->im,
				t00,t01,t02,t03,t04,t05,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2]
			);	tptr++;
		#if 0
			// If low nibble underflows need to add 16 to it and decrement high nibble ihi:
			ilo--;			lo_neg = (ilo < 0);	ilo+= ((-lo_neg) & 16);
			ihi -= 8+lo_neg;hi_neg = ihi< 0;	ihi+= ((-hi_neg) & 12);	// 1st high-part offset
			k1 = ihi- 8;	hi_neg = k1 < 0;	k1 += ((-hi_neg) & 12);	// high part of "decrement 64 horizontally" to get 2nd offset
			k2 = k1 - 8;	hi_neg = k2 < 0;	k2 += ((-hi_neg) & 12);	// high part of "decrement 64 horizontally" to get 3rd offset
			k0 = phi[ihi];	k1 = phi[k1];	k2 = phi[k2];	// Use distinct ihi and k0 = phi[ihi] for 1st high-part offset because basic-index ihi must survive from one loop pass to next
		#else
			// Incr [A-D]-perm cshift counter (mod 3):
			idx++;
			idx += (idx == 3);	// When idx hits 3, bump up to 4...
			idx &= 3;	// ...and reduce (mod 4), which is cheaper than (mod 3)
			// If low nibble underflows, += 16 and decrement [A-D] index (mod 4), which ==> jdx += 3 (mod 4):
			ilo--;
			if(ilo < 0) {	// This branch only taken every 16th loop exec
				ilo += 16;
				jdx = (jdx + 3)&3;
				iptr = &p10_3perm[jdx][0];
			}
			k0 = *(iptr+idx);	k1 = *(iptr+idx+1);	k2 = *(iptr+idx+2);
		#endif
		}
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy192_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr;
	#ifdef USE_SSE2
		const int pfetch_dist = PFETCH_DIST;
	#endif
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
			,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0;
		int poff[RADIX>>2];	// Store mults of p4 offset for loop control
	  #ifndef USE_SSE2
		double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	  #endif
	// Shared DIF+DIT:
		int ilo,ihi,idx,jdx,hi_neg,k0,k1,k2,nshift, *iptr;
		int plo[16], t_offsets[64];
	  #ifndef USE_SSE2
		int phi[12];
	  #endif
		uint64 i64;
	// DIF:
		int dif_o_offsets[RADIX];
		// Consts containing [in little-endian hex-char fashion] the 4 [A-D] 16-subperms defined in the output-ordering commentary below:
		const uint64 dif_16_operms[4] = {
			0x01235476ab98fecdull,	// [A]
			0x54762310fecd98baull,	// [B]
			0xab98fecd54762310ull,	// [C]
			0xfecd98ba23107645ull	// [D]
		};
	// DIT:
		// To-be-inited array used to support circular (0-2)-element shifts of basic 3-perm p10-multiple patterns:
		int p10_3perm[4][5] = { {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0} };
		int dit_i_offsets[RADIX];

		int j,j1,jt,jp,l;
	#ifndef USE_SSE2
		int j2,ntmp;
	#else
		// incr = Carry-chain wts-multipliers recurrence length;
		//incr must divide RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 12|24|48 for avx512,avx,sse, respectively:
		int incr;
	  #ifdef USE_AVX512
		const int incr_long = 12,incr_med = 6,incr_short = 4;
	  #else
		const int incr_long = 12,incr_med = 8,incr_short = 4;
	  #endif
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
		vec_dbl *cc0,/* *ss0, */ *max_err,
		  #ifndef USE_AVX512
			*sse2_rnd,
		  #endif
			*half_arr,
			*r00,*r40,*r80,	// Head of RADIX*vec_cmplx-sized local store #1
			*s1p00,	// Head of RADIX*vec_cmplx-sized local store #2
			*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx
		vec_dbl *tmp,*tm1,*tm2,	// Non-static utility ptrs
			*va0,*va1,*va2, *vb0,*vb1,*vb2,
			*vc0,*vc1,*vc2, *vd0,*vd1,*vd2;
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
		int m,m2;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		int bjmodn[RADIX];	// Thread only carries a base datum here, must alloc a local array for remaining values
		double *cy = thread_arg->cy, temp,frac,
			t00,t01,t02,t03,t04,t05;
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
		NDIVR >>= 4;			pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
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

	// Shared:
		plo[0x0] =  0; plo[0x1] = p1; plo[0x2] = p2; plo[0x3] = p3;
		plo[0x4] = p4; plo[0x5] = p5; plo[0x6] = p6; plo[0x7] = p7;
		plo[0x8] = p8; plo[0x9] = p9; plo[0xa] = pa; plo[0xb] = pb;
		plo[0xc] = pc; plo[0xd] = pd; plo[0xe] = pe; plo[0xf] = pf;

	#ifndef USE_SSE2
		phi[0x0] =   0; phi[0x1] = p10; phi[0x2] = p20; phi[0x3] = p30;
		phi[0x4] = p40; phi[0x5] = p50; phi[0x6] = p60; phi[0x7] = p70;
		phi[0x8] = p80; phi[0x9] = p90; phi[0xa] = pa0; phi[0xb] = pb0;
	#endif

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

	// DIF:
	// Extract the 4 [A-D] 16-subperms defined in the output-ordering commentary below:
	/* Operm 1:
		0,1,2,3,5,4,7,6,a,b,9,8,f,e,c,d + p00  =  [A] + p00
		5,4,7,6,2,3,1,0,f,e,c,d,9,8,b,a + p10  =  [B] + p10
		a,b,9,8,f,e,c,d,5,4,7,6,2,3,1,0 + p20  =  [C] + p20
		f,e,c,d,9,8,b,a,2,3,1,0,7,6,4,5 + p30  =  [D] + p30
	*/
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
			dif_o_offsets[0x20+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dif_16_operms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x30+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
	/* Operm 2:
		5,4,7,6,2,3,1,0,f,e,c,d,9,8,b,a + p90  =  [B] + p10
		a,b,9,8,f,e,c,d,5,4,7,6,2,3,1,0 + p80  =  [C] + p00 (mod p40)
		f,e,c,d,9,8,b,a,2,3,1,0,7,6,4,5 + pb0  =  [D] + p30
		5,4,7,6,2,3,1,0,f,e,c,d,9,8,b,a + pa0  =  [B] + p20
	*/
		i64 = dif_16_operms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x40+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x50+l] = plo[(i64 >> nshift)&0xf];
		}
		i64 = dif_16_operms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x60+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dif_16_operms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x70+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
	/* Operm 3:
		a,b,9,8,f,e,c,d,5,4,7,6,2,3,1,0 + p60  =  [C] + p20
		f,e,c,d,9,8,b,a,2,3,1,0,7,6,4,5 + p70  =  [D] + p30 (mod p40)
		5,4,7,6,2,3,1,0,f,e,c,d,9,8,b,a + p50  =  [B] + p10
		a,b,9,8,f,e,c,d,5,4,7,6,2,3,1,0 + p40  =  [C] + p00
	*/
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x80+l] = plo[(i64 >> nshift)&0xf] + p20;
		}
		i64 = dif_16_operms[3];	// [D] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0x90+l] = plo[(i64 >> nshift)&0xf] + p30;
		}
		i64 = dif_16_operms[1];	// [B] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0xa0+l] = plo[(i64 >> nshift)&0xf] + p10;
		}
		i64 = dif_16_operms[2];	// [C] subperm
		for(l = 0; l < 16; l++) {
			nshift = (15 - l)<<2;
			dif_o_offsets[0xb0+l] = plo[(i64 >> nshift)&0xf];
		}

	// DIT:
	/* Iperm 1:
		0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00
		f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10
		f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p30
		f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p20
	*/																		// swap [0x2*] and [0x3*] cols:
		dit_i_offsets[0x00] =  0;	dit_i_offsets[0x10] = pf+p10;	dit_i_offsets[0x30] = pf+p20;	dit_i_offsets[0x20] = pf+p30;
		dit_i_offsets[0x01] = p1;	dit_i_offsets[0x11] = pe+p10;	dit_i_offsets[0x31] = pe+p20;	dit_i_offsets[0x21] = pe+p30;
		dit_i_offsets[0x02] = p3;	dit_i_offsets[0x12] = pd+p10;	dit_i_offsets[0x32] = pd+p20;	dit_i_offsets[0x22] = pd+p30;
		dit_i_offsets[0x03] = p2;	dit_i_offsets[0x13] = pc+p10;	dit_i_offsets[0x33] = pc+p20;	dit_i_offsets[0x23] = pc+p30;
		dit_i_offsets[0x04] = p7;	dit_i_offsets[0x14] = pb+p10;	dit_i_offsets[0x34] = pb+p20;	dit_i_offsets[0x24] = pb+p30;
		dit_i_offsets[0x05] = p6;	dit_i_offsets[0x15] = pa+p10;	dit_i_offsets[0x35] = pa+p20;	dit_i_offsets[0x25] = pa+p30;
		dit_i_offsets[0x06] = p5;	dit_i_offsets[0x16] = p9+p10;	dit_i_offsets[0x36] = p9+p20;	dit_i_offsets[0x26] = p9+p30;
		dit_i_offsets[0x07] = p4;	dit_i_offsets[0x17] = p8+p10;	dit_i_offsets[0x37] = p8+p20;	dit_i_offsets[0x27] = p8+p30;
		dit_i_offsets[0x08] = pf;	dit_i_offsets[0x18] = p7+p10;	dit_i_offsets[0x38] = p7+p20;	dit_i_offsets[0x28] = p7+p30;
		dit_i_offsets[0x09] = pe;	dit_i_offsets[0x19] = p6+p10;	dit_i_offsets[0x39] = p6+p20;	dit_i_offsets[0x29] = p6+p30;
		dit_i_offsets[0x0a] = pd;	dit_i_offsets[0x1a] = p5+p10;	dit_i_offsets[0x3a] = p5+p20;	dit_i_offsets[0x2a] = p5+p30;
		dit_i_offsets[0x0b] = pc;	dit_i_offsets[0x1b] = p4+p10;	dit_i_offsets[0x3b] = p4+p20;	dit_i_offsets[0x2b] = p4+p30;
		dit_i_offsets[0x0c] = pb;	dit_i_offsets[0x1c] = p3+p10;	dit_i_offsets[0x3c] = p3+p20;	dit_i_offsets[0x2c] = p3+p30;
		dit_i_offsets[0x0d] = pa;	dit_i_offsets[0x1d] = p2+p10;	dit_i_offsets[0x3d] = p2+p20;	dit_i_offsets[0x2d] = p2+p30;
		dit_i_offsets[0x0e] = p9;	dit_i_offsets[0x1e] = p1+p10;	dit_i_offsets[0x3e] = p1+p20;	dit_i_offsets[0x2e] = p1+p30;
		dit_i_offsets[0x0f] = p8;	dit_i_offsets[0x1f] =  0+p10;	dit_i_offsets[0x3f] =  0+p20;	dit_i_offsets[0x2f] =  0+p30;
	/* Iperm 2:
		5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p10
		5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p00 (mod p40)
		5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p20
		9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p30
	*/			// swap [0x4*] and [0x5*] cols:
		dit_i_offsets[0x50] = p5;	dit_i_offsets[0x40] = p5+p10;	dit_i_offsets[0x60] = p5+p20;	dit_i_offsets[0x70] = p9+p30;
		dit_i_offsets[0x51] = p4;	dit_i_offsets[0x41] = p4+p10;	dit_i_offsets[0x61] = p4+p20;	dit_i_offsets[0x71] = p8+p30;
		dit_i_offsets[0x52] = p6;	dit_i_offsets[0x42] = p6+p10;	dit_i_offsets[0x62] = p6+p20;	dit_i_offsets[0x72] = pa+p30;
		dit_i_offsets[0x53] = p7;	dit_i_offsets[0x43] = p7+p10;	dit_i_offsets[0x63] = p7+p20;	dit_i_offsets[0x73] = pb+p30;
		dit_i_offsets[0x54] = p1;	dit_i_offsets[0x44] = p1+p10;	dit_i_offsets[0x64] = p1+p20;	dit_i_offsets[0x74] = pe+p30;
		dit_i_offsets[0x55] =  0;	dit_i_offsets[0x45] =  0+p10;	dit_i_offsets[0x65] =  0+p20;	dit_i_offsets[0x75] = pf+p30;
		dit_i_offsets[0x56] = p2;	dit_i_offsets[0x46] = p2+p10;	dit_i_offsets[0x66] = p2+p20;	dit_i_offsets[0x76] = pc+p30;
		dit_i_offsets[0x57] = p3;	dit_i_offsets[0x47] = p3+p10;	dit_i_offsets[0x67] = p3+p20;	dit_i_offsets[0x77] = pd+p30;
		dit_i_offsets[0x58] = p9;	dit_i_offsets[0x48] = p9+p10;	dit_i_offsets[0x68] = p9+p20;	dit_i_offsets[0x78] = p1+p30;
		dit_i_offsets[0x59] = p8;	dit_i_offsets[0x49] = p8+p10;	dit_i_offsets[0x69] = p8+p20;	dit_i_offsets[0x79] =  0+p30;
		dit_i_offsets[0x5a] = pa;	dit_i_offsets[0x4a] = pa+p10;	dit_i_offsets[0x6a] = pa+p20;	dit_i_offsets[0x7a] = p2+p30;
		dit_i_offsets[0x5b] = pb;	dit_i_offsets[0x4b] = pb+p10;	dit_i_offsets[0x6b] = pb+p20;	dit_i_offsets[0x7b] = p3+p30;
		dit_i_offsets[0x5c] = pe;	dit_i_offsets[0x4c] = pe+p10;	dit_i_offsets[0x6c] = pe+p20;	dit_i_offsets[0x7c] = p6+p30;
		dit_i_offsets[0x5d] = pf;	dit_i_offsets[0x4d] = pf+p10;	dit_i_offsets[0x6d] = pf+p20;	dit_i_offsets[0x7d] = p7+p30;
		dit_i_offsets[0x5e] = pc;	dit_i_offsets[0x4e] = pc+p10;	dit_i_offsets[0x6e] = pc+p20;	dit_i_offsets[0x7e] = p4+p30;
		dit_i_offsets[0x5f] = pd;	dit_i_offsets[0x4f] = pd+p10;	dit_i_offsets[0x6f] = pd+p20;	dit_i_offsets[0x7f] = p5+p30;
	/* Iperm 3:
		a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p20
		2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p30 (mod p40)
		a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p00
		2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p10
	*/		// swap [0x8*] <=> [0xa*] colpair, and [0x9*] <=> [0xb*] pair:
		dit_i_offsets[0xa0] = pa;	dit_i_offsets[0xb0] = p2+p10;	dit_i_offsets[0x80] = pa+p20;	dit_i_offsets[0x90] = p2+p30;
		dit_i_offsets[0xa1] = pb;	dit_i_offsets[0xb1] = p3+p10;	dit_i_offsets[0x81] = pb+p20;	dit_i_offsets[0x91] = p3+p30;
		dit_i_offsets[0xa2] = p8;	dit_i_offsets[0xb2] =  0+p10;	dit_i_offsets[0x82] = p8+p20;	dit_i_offsets[0x92] =  0+p30;
		dit_i_offsets[0xa3] = p9;	dit_i_offsets[0xb3] = p1+p10;	dit_i_offsets[0x83] = p9+p20;	dit_i_offsets[0x93] = p1+p30;
		dit_i_offsets[0xa4] = pc;	dit_i_offsets[0xb4] = p4+p10;	dit_i_offsets[0x84] = pc+p20;	dit_i_offsets[0x94] = p4+p30;
		dit_i_offsets[0xa5] = pd;	dit_i_offsets[0xb5] = p5+p10;	dit_i_offsets[0x85] = pd+p20;	dit_i_offsets[0x95] = p5+p30;
		dit_i_offsets[0xa6] = pf;	dit_i_offsets[0xb6] = p7+p10;	dit_i_offsets[0x86] = pf+p20;	dit_i_offsets[0x96] = p7+p30;
		dit_i_offsets[0xa7] = pe;	dit_i_offsets[0xb7] = p6+p10;	dit_i_offsets[0x87] = pe+p20;	dit_i_offsets[0x97] = p6+p30;
		dit_i_offsets[0xa8] = p2;	dit_i_offsets[0xb8] = pc+p10;	dit_i_offsets[0x88] = p2+p20;	dit_i_offsets[0x98] = pc+p30;
		dit_i_offsets[0xa9] = p3;	dit_i_offsets[0xb9] = pd+p10;	dit_i_offsets[0x89] = p3+p20;	dit_i_offsets[0x99] = pd+p30;
		dit_i_offsets[0xaa] =  0;	dit_i_offsets[0xba] = pf+p10;	dit_i_offsets[0x8a] =  0+p20;	dit_i_offsets[0x9a] = pf+p30;
		dit_i_offsets[0xab] = p1;	dit_i_offsets[0xbb] = pe+p10;	dit_i_offsets[0x8b] = p1+p20;	dit_i_offsets[0x9b] = pe+p30;
		dit_i_offsets[0xac] = p4;	dit_i_offsets[0xbc] = p8+p10;	dit_i_offsets[0x8c] = p4+p20;	dit_i_offsets[0x9c] = p8+p30;
		dit_i_offsets[0xad] = p5;	dit_i_offsets[0xbd] = p9+p10;	dit_i_offsets[0x8d] = p5+p20;	dit_i_offsets[0x9d] = p9+p30;
		dit_i_offsets[0xae] = p7;	dit_i_offsets[0xbe] = pb+p10;	dit_i_offsets[0x8e] = p7+p20;	dit_i_offsets[0x9e] = pb+p30;
		dit_i_offsets[0xaf] = p6;	dit_i_offsets[0xbf] = pa+p10;	dit_i_offsets[0x8f] = p6+p20;	dit_i_offsets[0x9f] = pa+p30;

	// Init elts of array used to support circular (0-2)-element shifts of the basic 3-perm p10-multiple patterns
	// [A] = 00,40,80, [B] = 10,50,90, [C] = 20,60,a0, [D] = 30,70,b0:
	#ifndef USE_SSE2
		l = 0;	p10_3perm[0][l++] =   0; p10_3perm[0][l++] = p40; p10_3perm[0][l++] = p80; p10_3perm[0][l++] =   0; p10_3perm[0][l++] = p40;
		l = 0;	p10_3perm[1][l++] = p10; p10_3perm[1][l++] = p50; p10_3perm[1][l++] = p90; p10_3perm[1][l++] = p10; p10_3perm[1][l++] = p50;
		l = 0;	p10_3perm[2][l++] = p20; p10_3perm[2][l++] = p60; p10_3perm[2][l++] = pa0; p10_3perm[2][l++] = p20; p10_3perm[2][l++] = p60;
		l = 0;	p10_3perm[3][l++] = p30; p10_3perm[3][l++] = p70; p10_3perm[3][l++] = pb0; p10_3perm[3][l++] = p30; p10_3perm[3][l++] = p70;
	#else	// Double the RHS immediately here to avoid awkward loop over 2D array:
		l = 0;	p10_3perm[0][l++] =   0; p10_3perm[0][l++] =0x80; p10_3perm[0][l++]=0x100; p10_3perm[0][l++] =   0; p10_3perm[0][l++] =0x80;
		l = 0;	p10_3perm[1][l++] =0x20; p10_3perm[1][l++] =0xa0; p10_3perm[1][l++]=0x120; p10_3perm[1][l++] =0x20; p10_3perm[1][l++] =0xa0;
		l = 0;	p10_3perm[2][l++] =0x40; p10_3perm[2][l++] =0xc0; p10_3perm[2][l++]=0x140; p10_3perm[2][l++] =0x40; p10_3perm[2][l++] =0xc0;
		l = 0;	p10_3perm[3][l++] =0x60; p10_3perm[3][l++] =0xe0; p10_3perm[3][l++]=0x160; p10_3perm[3][l++] =0x60; p10_3perm[3][l++] =0xe0;
	#endif

	#ifdef USE_SSE2
		tmp = r00 = thread_arg->r00;	// Head of RADIX*vec_cmplx-sized local store #1
		r40 = tmp + 0x080;
		r80 = tmp + 0x100;
		tmp += 0x180;	s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
		tmp += 0x180;
		// DFT-64 roots all def'd locally in the DFT-64 functions, only need the radix-3 roots here:
		cc0	= tmp;
		//ss0	= tmp + 1;
		tmp += 0x2;	// r00 += 0x302
	  #ifdef USE_AVX512
		cy = tmp;		tmp += 0x18;	// RADIX/8 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		//sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;	// sc_ptr += 0x31c; This is where the value of half_arr_offset192 comes from
	  #elif defined(USE_AVX)
		cy = tmp;		tmp += 0x30;	// RADIX/4 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x334; This is where the value of half_arr_offset192 comes from
		half_arr= tmp + 0x02;	// This table needs 68*SZ_VD bytes in avx mode
	  #else
		cy = tmp;		tmp += 0x60;	// RADIX/2 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x364; This is where the value of half_arr_offset192 comes from
		half_arr= tmp + 0x02;	// This table needs 20*SZ_VD bytes in sse2 mode
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

		sign_mask = (uint64*)(r00 + radix192_creals_in_local_store);
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
		base    = (double *)thread_arg->r00    ;
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
		#include "radix192_main_carry_loop.h"

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
