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
#include "radix1024.h"

#define RADIX 1024	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

#define EPS 1e-10

#ifndef RAD1
	#define RAD1	64	// Sets the first-pass (twiddleless DFTs) radix of the 2-pass decomposition of 1024
#endif
#if (RAD1 != 32) && (RAD1 != 64)
	#error Currently only pass-1 radices 32 and 64 supported!
#endif

#ifndef PFETCH_DIST
  #ifdef USE_AVX512
	#define PFETCH_DIST	64	// Feb 2017: Test on KNL point to this as best
  #elif defined(USE_AVX)
	#define PFETCH_DIST	32	// This seems to work best on my Haswell, even though 64 bytes seems more logical in AVX mode
  #else
	#define PFETCH_DIST	32
  #endif
#endif

// SIMD+SSE2 code only available for GCC build:
#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC)

	#include "sse2_macro_gcc64.h"

  // For Mersenne-mod need (16 [SSE2] or 64 [AVX]) + (4 [HIACC] or 40 [LOACC]) added slots for half_arr lookup tables.
  // Max = (40 [SSE2]; 132 [AVX]),
  // For Fermat-mod we use RADIX*4 = 4096 [note there is no Fermat-mod LOACC option for this power-of-2 DFT] more
  // slots in AVX mode for the compact negacyclic-roots chained-multiply scheme. Add larger of the 2 numbers -
  // 4096 for AVX, 40 for SSE2 - to (half_arr_offset1024 + RADIX) to get SIMD value of radix1024_creals_in_local_store:
  #ifdef USE_AVX512
	const int half_arr_offset1024 = 0x1905;	// 0x100 = 2*(RADIX/8) fewer cy-slots than in AVX mode
	const int radix1024_creals_in_local_store = 0x2d08;	// (half_arr_offset1024 + 5*RADIX) and round up to nearest multiple of 4
	#include "radix1024_avx_negadwt_consts.h"
  #elif defined(USE_AVX)
	const int half_arr_offset1024 = 0x1a05;	// + RADIX = 0x1e05; Used for thread local-storage-integrity checking
	const int radix1024_creals_in_local_store = 0x2e08;	// (half_arr_offset1024 + 5*RADIX) and round up to nearest multiple of 4
	#include "radix1024_avx_negadwt_consts.h"
  #else
	const int half_arr_offset1024 = 0x1c05;	// + RADIX = 0x2005; Used for thread local-storage-integrity checking
	const int radix1024_creals_in_local_store = 0x2030;	// (half_arr_offset1024 + RADIX) + 0x28 and round up to nearest multiple of 4
  #endif

#elif defined(USE_SSE2)

	#error SIMD build only supported for GCC-compatible compiler under *nix/macOS!

#endif	/* USE_SSE2 */

#ifdef USE_PTHREAD

	// Use non-pooled simple spawn/rejoin thread-team model
	#include "threadpool.h"

	struct cy_thread_data_t{
	// int data - if needed, pad to yield an even number of these:
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
		struct complex *rn0;
		struct complex *rn1;
	#ifdef USE_SSE2
		vec_dbl *r00;
		vec_dbl *half_arr;
	#else
		double *r00;
		double *half_arr;
	#endif
		uint64*sm_ptr;
		uint32 bjmodnini;
		int bjmodn0;
	// For large radix0 use thread-local arrays for DWT indices/carries - only caveat is these must be SIMD-aligned:
	// Since GCC uses __BIGGEST_ALIGNMENT__ = 16 on x86, which is too small to be useful for avx data,
	// we are forced to resort to fugly hackage - add pad slots to a garbage-named struct-internal array along with
	// a pointer-to-be-inited-at-runtime, when we set ptr to the lowest-index array element having the desired alginment:
		double *cy_r,*cy_i;
	#ifdef USE_AVX512
		double cy_dat[2*RADIX+8] __attribute__ ((__aligned__(8)));
	#else
		double cy_dat[2*RADIX+4] __attribute__ ((__aligned__(8)));	// Enforce min-alignment of 8 bytes in 32-bit builds.
	#endif
	};

#endif

/***************/

int radix1024_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-1024 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-1024 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix16_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
#if !defined(MULTITHREAD) || defined(USE_SSE2)
	#include "radix1024_twiddles.h"
#endif
	const char func[] = "radix1024_ditN_cy_dif1";
	static int thr_id = 0;	// Master thread gets this special id
	const int pfetch_dist = PFETCH_DIST;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
  #ifdef USE_AVX512
	const int jhi_wrap_mers = 15;
	const int jhi_wrap_ferm = 15;
  #else
	const int jhi_wrap_mers =  7;
	const int jhi_wrap_ferm = 15;	// For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
  #endif
	// FMA-based DFT needs the tangent:
#if defined(USE_AVX2) && !defined(USE_IMCI512)
	static double tan = 0.41421356237309504879;
#endif
	int NDIVR,i,incr,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,ntmp,outer,nbytes;
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 64|128|256 for avx512,avx,sse, respectively:
	const int incr_long = 16,incr_med = 8,incr_short = 4;
	// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
	if(USE_SHORT_CY_CHAIN == 0)
		incr = incr_long;
	else if(USE_SHORT_CY_CHAIN == 1)
		incr = incr_med;
	else
		incr = incr_short;
	int k1,k2;
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
	static double wts_mult[2], inv_mult[2];	// Const wts-multiplier and 2*(its multiplicative inverse)
	double rt,it, wt_re,wt_im, wi_re,wi_im;	// Fermat-mod weights stuff, used in both scalar and AVX mode
	static uint32 bjmodnini;
	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,nm1,
		p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
		p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p1a,p1b,p1c,p1d,p1e,p1f,
		p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,
		p100,p110,p120,p130,p140,p150,p160,p170,p180,p190,p1a0,p1b0,p1c0,p1d0,p1e0,p1f0,p200, nsave = 0;
	static int poff[RADIX>>2];
// DIF:
	static int dif_i_offsets[64], dif_po_br[16];
// DIT:
	int ju,jv;
	static int dit_i_offsets[64], dit_poffs[16], dit_po_br[64];
// Shared by both:
	static int dif_o_offsets[64], o_offsets[64];

	static double radix_inv, n2inv;
	double scale, dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex t[RADIX], *tptr;
	const double *addr,*addi;
	int *itmp,*itm2;	// Pointer into the bjmodn array
	int err;
	static int first_entry=TRUE;

	int n_div_nwt;

#ifdef USE_SSE2

	int idx_offset,idx_incr;
	static int cslots_in_local_store;
	static vec_dbl *one = 0x0, one_dat;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_nm1;
	uint64 tmp64;

  #ifndef PFETCH_DIST
  #ifdef USE_AVX512
	#define PFETCH_DIST	64	// Feb 2017: Test on KNL point to this as best
  #elif defined(USE_AVX)
	#define PFETCH_DIST	32	// This seems to work best on my Haswell, even though 64 bytes seems more logical in AVX mode
  #else
	#define PFETCH_DIST	32
  #endif
#endif

#ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else
	double *add0,*add1,*add2,*add3,*add4,*add5,*add6,*add7,*add8,*add9,*adda,*addb,*addc,*addd,*adde,*addf;
  #endif	// MULTITHREAD

	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
	const double crnd = 3.0*0x4000000*0x2000000;
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	static vec_dbl *max_err, *sse2_rnd, *half_arr, *isrt2,*cc0,*ss0,	// Need DFT-16 roots explicitly
		*twid00,// ptr to head of 64 sets of twiddles shared by the 2nd-half DIF and DIT DFT macros
		*r00,	// Head of RADIX*vec_cmplx-sized local store #1
		*s1p00,	// Head of RADIX*vec_cmplx-sized local store #2
		*vd00,	// Head of 64 x vec_cmplx loal store for DFT-64 intermediates
		*cy_r,*cy_i;	// Need RADIX slots for sse2 carries, RADIX/2 for avx
  #ifdef USE_AVX
	static vec_dbl *base_negacyclic_root;
  #endif

	vec_dbl *tmp,*tm0,*tm1,*tm2;	// Non-static utility ptrs

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
	static task_control_t   task_control = {NULL, (void*)cy1024_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	int m,m2;
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	int bjmodn[RADIX];
	double temp,frac,
		cy_r[RADIX],cy_i[RADIX];

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static int *_bjmodnini = 0x0, *_bjmodn[RADIX];
	static double *_cy_r[RADIX],*_cy_i[RADIX];
	if(!_jhi) {
		_cy_r[0] = 0x0;	// First of these used as an "already inited consts?" sentinel, must init = 0x0 at same time do so for non-array static ptrs
	}

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=-1;
	// Jan 2018: To support PRP-testing, read the LR-modpow-scalar-multiply-needed bit for the current iteration from the global array:
	double prp_mult = 1.0;
	// v18: If use residue shift in context of Pépin test, need prp_mult = 2 whenever the 'shift = 2*shift + random[0,1]' update gets a 1-bit in the random slot
	if((TEST_TYPE == TEST_TYPE_PRIMALITY && MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	|| (TEST_TYPE & 0xfffffffe) == TEST_TYPE_PRP) {	// Mask off low bit to lump together PRP and PRP-C tests
		i = (iter-1) % ITERS_BETWEEN_CHECKPOINTS;	// Bit we need to read...iter-counter is unit-offset w.r.to iter-interval, hence the -1
		if((BASE_MULTIPLIER_BITS[i>>6] >> (i&63)) & 1)
			prp_mult = PRP_BASE;
	}

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n >> 10;	ndivr_inv = (double)RADIX/n;
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

		nm1   = n-1;

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
		else
		{
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
			tdat[ithread].rn0 = rn0;
			tdat[ithread].rn1 = rn1;

		// This array pointer must be set based on vec_dbl-sized alignment at runtime for each thread:
			for(l = 0; l < RE_IM_STRIDE; l++) {
				if( ((intptr_t)&tdat[ithread].cy_dat[l] & SZ_VDM1) == 0 ) {
					tdat[ithread].cy_r = &tdat[ithread].cy_dat[l];
					tdat[ithread].cy_i = tdat[ithread].cy_r + RADIX;
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

		// Use vector-double type size (16 bytes for SSE2, 32 for AVX) to alloc a block of local storage
		cslots_in_local_store = radix1024_creals_in_local_store + (20+RADIX/2)/2;	// Just add enough int64 space for both cases, plus some
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix1024_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = sc_ptr;	r00   = tmp;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x800;	s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
		tmp += 0x800;	vd00  = tmp;	// Head of 64 x vec_cmplx = 128 vec_dbl local store for DFT-64 intermediates
		tmp += 0x80;	// sc_ptr += 0x1080
		// DFT-16 roots needed explicitly:
		isrt2  = tmp + 0x00;
		cc0    = tmp + 0x01;
		ss0    = tmp + 0x02;
		tmp += 0x3;
		// ptrs to head of 64 sets (2*15 = 30 vec_dbl data each) of non-unity twiddles shared by the 2nd-half DIF and DIT DFT macros:
		twid00  = tmp + 0x00;
		tmp += 0x780;	// sc_ptr += 0x1083 + 4*0x1e0 = 0x1803
	  #ifdef USE_AVX512
		cy_r = tmp;	cy_i = tmp+0x080;	tmp += 0x100;	// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
		cy_r = tmp;	cy_i = tmp+0x100;	tmp += 0x200;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x1803 + 0x202 = 0x1a05... This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 96 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy_r = tmp;	cy_i = tmp+0x200;	tmp += 0x400;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x1803 + 0x402 = 0x1c05... This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 32 x 16 bytes for Mersenne-mod, 2 for Fermat-mod */
	  #endif

		ASSERT(half_arr_offset1024 == (uint32)(half_arr-sc_ptr), "half_arr_offset mismatches actual!");
		ASSERT((radix1024_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r00) + (20 << L2_SZ_VD), "radix1024_creals_in_local_store checksum failed!");

		/* These remain fixed: */
		VEC_DBL_INIT(isrt2,ISRT2);
		VEC_DBL_INIT(cc0  ,  c16);
		VEC_DBL_INIT(ss0  ,  s16);
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

		// ptrs to 16 sets (30 vec_dbl data each) of non-unity twiddles shared by the 2nd-half DIF and DIT DFT macros.
		// Since we copied the init-blocks here from the code below in which the twiddle-sets appear in BR order, init same way:

		// Use loop-based init for codecompactness and extensibility to larger pow2 radices:
	  #if defined(USE_AVX2) && !defined(USE_IMCI512)
		// The first 2 sets (twid00/20) are processed in non-FMA fashion by both DIF/DIT macros, so init in non-FMA fashion:
		for(l = 0; l < 2; l++) {
	  #else
		for(l = 0; l < 64; l++) {
	  #endif
			j = reverse(l,6)<<1;	// twid00-offsets are processed in BR order
			tmp = twid00 + (j<<4)-j;	// Twid-offsets are multiples of 30 vec_dbl
			addr = DFT1024_TWIDDLES[l]; addi = addr+1;	// Pointer to required row of 2-D twiddles array
			VEC_DBL_INIT(tmp,*(addr+0x00)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x00)); tmp++;
			VEC_DBL_INIT(tmp,*(addr+0x02)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x02)); tmp++;
			VEC_DBL_INIT(tmp,*(addr+0x04)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x04)); tmp++;
			VEC_DBL_INIT(tmp,*(addr+0x06)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x06)); tmp++;
			VEC_DBL_INIT(tmp,*(addr+0x08)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x08)); tmp++;
			VEC_DBL_INIT(tmp,*(addr+0x0a)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x0a)); tmp++;
			VEC_DBL_INIT(tmp,*(addr+0x0c)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x0c)); tmp++;
			VEC_DBL_INIT(tmp,*(addr+0x0e)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x0e)); tmp++;
			VEC_DBL_INIT(tmp,*(addr+0x10)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x10)); tmp++;
			VEC_DBL_INIT(tmp,*(addr+0x12)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x12)); tmp++;
			VEC_DBL_INIT(tmp,*(addr+0x14)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x14)); tmp++;
			VEC_DBL_INIT(tmp,*(addr+0x16)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x16)); tmp++;
			VEC_DBL_INIT(tmp,*(addr+0x18)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x18)); tmp++;
			VEC_DBL_INIT(tmp,*(addr+0x1a)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x1a)); tmp++;
			VEC_DBL_INIT(tmp,*(addr+0x1c)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x1c)); tmp++;
		}

		// The remaining 14 sets are inited differently depending on whether SIMD+FMA is used:
	  #if defined(USE_AVX2) && !defined(USE_IMCI512)

		// Precompute the FMA-modified twiddles for the 2nd-pass radix-16 DFTs:
		#ifdef USE_FMA
			#error USE_FMA flag not supported in SIMD mode - to use FMA under AVX2/FMA3, define *only* USE_AVX2!
		#endif

		#include "radix16_dif_dit_pass_asm.h"	// Need this for FMA_TWIDDLE_FIDDLE macro

		// Init the vec_dbl const 1.0:
		one = &one_dat;	VEC_DBL_INIT(one, 1.0);
		for(l = 2; l < 64; l++) {
			j = reverse(l,6)<<1;	// twid00-offsets are processed in BR order
			tmp = twid00 + (j<<4)-j;	// Twid-offsets are multiples of 30 vec_dbl
			addr = DFT1024_TWIDDLES[l]; addi = addr+1;	// Pointer to required row of 2-D twiddles array
			FMA_TWIDDLE_FIDDLE(
				*(addr+0x00),*(addi+0x00), *(addr+0x02),*(addi+0x02), *(addr+0x04),*(addi+0x04), *(addr+0x06),*(addi+0x06), *(addr+0x08),*(addi+0x08), *(addr+0x0a),*(addi+0x0a), *(addr+0x0c),*(addi+0x0c), *(addr+0x0e),*(addi+0x0e), *(addr+0x10),*(addi+0x10), *(addr+0x12),*(addi+0x12), *(addr+0x14),*(addi+0x14), *(addr+0x16),*(addi+0x16), *(addr+0x18),*(addi+0x18), *(addr+0x1a),*(addi+0x1a), *(addr+0x1c),*(addi+0x1c),
				c16,tan,
				tmp
			)
		}

	#endif

		// Propagate the above consts to the remaining threads:
		nbytes = (intptr_t)cy_r - (intptr_t)isrt2;	// #bytes in 1st of above block of consts
		tmp = isrt2;
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

		if(TRANSFORM_TYPE == RIGHT_ANGLE)
		{
			/* In Fermat-mod mode, use first 2 SIMD-sized slots for base and 1/base: */
			VEC_DBL_INIT(tmp, base   [0]);	++tmp;
			VEC_DBL_INIT(tmp, baseinv[0]);	++tmp;
			/* [+2] slot is for [scale,scale] */

			// Propagate the above consts to the remaining threads:
			nbytes = 2 << L2_SZ_VD;
			tmp = half_arr;
			tm2 = tmp + cslots_in_local_store;
			for(ithread = 1; ithread < CY_THREADS; ++ithread) {
				memcpy(tm2, tmp, nbytes);
				tmp = tm2;		tm2 += cslots_in_local_store;
			}

		#ifdef USE_AVX

			base_negacyclic_root = half_arr + RADIX;

			/*
			The pattern of the negacyclic-DWT-weights ("nDWTs") applied to the RADIX complex outputs of the final-DIT-pass is thus:
			The n DWTs multiplying each set of RADIX DIT DFT outputs are the product of a single complex-root "base multiplier"
			rbase (separately computed for each batch of DFT outputs), said "base root" multiplies the (0 - RADIX-1)st [4*RADIX]th
			roots of unity, i.e.

				 rbase * (j*I*Pi/2)/RADIX, for j = 0, ..., RADIX-1 .

			See the radix28 version of this routine for additional details.
			*/
			#if 0
			// Simple qfloat-based loop to crank out the roots which populate the radix1024_avx_negadwt_consts table:
				struct qfloat qc,qs,qx,qy,qt,qn,qmul;
				qt = qfdiv(QPIHALF, i64_to_q((int64)RADIX));	// Pi/2/RADIX
				qc = qfcos(qt);	qs = qfsin(qt);
				qx = QONE;		qy = QZRO;
				for(j = 0; j < RADIX; j++) {
					printf("j = %3u: cos = %#16" PRIX64 "\n",j,qfdbl_as_uint64(qx));
					// Up-multiply the complex exponential:
					qn = qfmul(qx, qc); qt = qfmul(qy, qs); qmul = qfsub(qn, qt);	// Store qxnew in qmul for now.
					qn = qfmul(qx, qs); qt = qfmul(qy, qc); qy   = qfadd(qn, qt); qx = qmul;
				}
				exit(0);
			#endif

		  #ifdef USE_AVX512	// 8-way-double analog of AVX inits below:

			tmp = base_negacyclic_root + 2*RADIX;	// First 2*RADIX slots reserved for RADIX/8 copies of the Re/Im parts of the 8 base multipliers
			tm2 = tmp + RADIX/4 - 1;
			// First elt-pair needs special handling - have the 1.0 in avx_negadwt_consts[0] but the sine term buggers things
			tmp->d0 = 1.0;	(tmp+1)->d0 = 0.0;
			tmp64 = radix1024_avx_negadwt_consts[1];	tmp->d1 = tm2->d7 = *(double *)&tmp64;
			tmp64 = radix1024_avx_negadwt_consts[2];	tmp->d2 = tm2->d6 = *(double *)&tmp64;
			tmp64 = radix1024_avx_negadwt_consts[3];	tmp->d3 = tm2->d5 = *(double *)&tmp64;
			tmp64 = radix1024_avx_negadwt_consts[4];	tmp->d4 = tm2->d4 = *(double *)&tmp64;
			tmp64 = radix1024_avx_negadwt_consts[5];	tmp->d5 = tm2->d3 = *(double *)&tmp64;
			tmp64 = radix1024_avx_negadwt_consts[6];	tmp->d6 = tm2->d2 = *(double *)&tmp64;
			tmp64 = radix1024_avx_negadwt_consts[7];	tmp->d7 = tm2->d1 = *(double *)&tmp64;	tmp += 2;
			for(j = 8; j < RADIX; j += 8) {
				tmp64 = radix1024_avx_negadwt_consts[j+0];	tmp->d0 = tm2->d0 = *(double *)&tmp64;	tm2 -= 2;
				tmp64 = radix1024_avx_negadwt_consts[j+1];	tmp->d1 = tm2->d7 = *(double *)&tmp64;
				tmp64 = radix1024_avx_negadwt_consts[j+2];	tmp->d2 = tm2->d6 = *(double *)&tmp64;
				tmp64 = radix1024_avx_negadwt_consts[j+3];	tmp->d3 = tm2->d5 = *(double *)&tmp64;
				tmp64 = radix1024_avx_negadwt_consts[j+4];	tmp->d4 = tm2->d4 = *(double *)&tmp64;
				tmp64 = radix1024_avx_negadwt_consts[j+5];	tmp->d5 = tm2->d3 = *(double *)&tmp64;
				tmp64 = radix1024_avx_negadwt_consts[j+6];	tmp->d6 = tm2->d2 = *(double *)&tmp64;
				tmp64 = radix1024_avx_negadwt_consts[j+7];	tmp->d7 = tm2->d1 = *(double *)&tmp64;	tmp += 2;
			}
			tmp = base_negacyclic_root + RADIX*2;	// reset to point to start of above block

		  #else

			tmp = base_negacyclic_root + RADIX*2;	// First 2*RADIX slots reserved for RADIX/4 copies of the Re/Im parts of the 4 base multipliers
			tm2 = tmp + RADIX/2 - 1;
			// First elt-pair needs special handling - have the 1.0 in avx_negadwt_consts[0] but the sine term buggers things
			tmp->d0 = 1.0;	(tmp+1)->d0 = 0.0;
			tmp64 = radix1024_avx_negadwt_consts[1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(  1*I*Pi/(2*RADIX)) = sin((RADIX-  1)*I*Pi/(2*RADIX)) */
			tmp64 = radix1024_avx_negadwt_consts[2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(  2*I*Pi/(2*RADIX)) = sin((RADIX-  2)*I*Pi/(2*RADIX)) */
			tmp64 = radix1024_avx_negadwt_consts[3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(  3*I*Pi/(2*RADIX)) = sin((RADIX-  3)*I*Pi/(2*RADIX)) */	tmp += 2;
			for(j = 4; j < RADIX; j += 4) {
				tmp64 = radix1024_avx_negadwt_consts[j+0];	tmp->d0 = tm2->d0 = *(double *)&tmp64;	tm2 -= 2;
				tmp64 = radix1024_avx_negadwt_consts[j+1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;
				tmp64 = radix1024_avx_negadwt_consts[j+2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;
				tmp64 = radix1024_avx_negadwt_consts[j+3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	tmp += 2;
			}

			tmp = base_negacyclic_root + RADIX*2;	// reset to point to start of above block

		  #endif

			nbytes = RADIX << (L2_SZ_VD-1);	// RADIX*SZ_VD/2; 7 AVX-register-sized complex data

			// Propagate the above consts to the remaining threads:
			tm2 = tmp + cslots_in_local_store;
			for(ithread = 1; ithread < CY_THREADS; ++ithread) {
				memcpy(tm2, tmp, nbytes);
				tmp = tm2;		tm2 += cslots_in_local_store;
			}

		#endif	// AVX?
		}
		else
		{
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

		#else	// USE_SSE2

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

		}	// TRANSFORM_TYPE toggle

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

		sse_nm1 = sse_sw + RE_IM_STRIDE;
		tmp64 = (uint64)nm1;
		tmp64 = tmp64 + (tmp64 << 32);
		for(i = 0; i < RE_IM_STRIDE; ++i) {
			*(sse_nm1 +i) = tmp64;
		}

		nbytes = 4 << L2_SZ_VD;

	  #ifdef USE_AVX512
	   #ifdef CARRY_16_WAY
		n_minus_sil   = (struct uint32x16*)sse_nm1 + 1;
		n_minus_silp1 = (struct uint32x16*)sse_nm1 + 2;
		sinwt         = (struct uint32x16*)sse_nm1 + 3;
		sinwtm1       = (struct uint32x16*)sse_nm1 + 4;
		nbytes += 256;
	   #else
		n_minus_sil   = (struct uint32x8 *)sse_nm1 + 1;
		n_minus_silp1 = (struct uint32x8 *)sse_nm1 + 2;
		sinwt         = (struct uint32x8 *)sse_nm1 + 3;
		sinwtm1       = (struct uint32x8 *)sse_nm1 + 4;
		nbytes += 128;
	   #endif
	  #elif defined(USE_AVX)
		n_minus_sil   = (struct uint32x4 *)sse_nm1 + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_nm1 + 2;
		sinwt         = (struct uint32x4 *)sse_nm1 + 3;
		sinwtm1       = (struct uint32x4 *)sse_nm1 + 4;
		nbytes += 64;
	  #endif

		// Propagate the above consts to the remaining threads:
		tmp = (vec_dbl *)sm_ptr;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	  #ifdef USE_AVX
		bjmodn = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn = (int*)(sse_nm1 + RE_IM_STRIDE);
	  #endif

	#endif	// USE_SSE2

		/*   constant index offsets for load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
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
		p10 = pf + NDIVR;		pf += ( (pf >> DAT_BITS) << PAD_BITS );
		p11 = p10 + NDIVR;		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p12 = p11 + NDIVR;		p11 += ( (p11 >> DAT_BITS) << PAD_BITS );
		p13 = p12 + NDIVR;		p12 += ( (p12 >> DAT_BITS) << PAD_BITS );
		p14 = p13 + NDIVR;		p13 += ( (p13 >> DAT_BITS) << PAD_BITS );
		p15 = p14 + NDIVR;		p14 += ( (p14 >> DAT_BITS) << PAD_BITS );
		p16 = p15 + NDIVR;		p15 += ( (p15 >> DAT_BITS) << PAD_BITS );
		p17 = p16 + NDIVR;		p16 += ( (p16 >> DAT_BITS) << PAD_BITS );
		p18 = p17 + NDIVR;		p17 += ( (p17 >> DAT_BITS) << PAD_BITS );
		p19 = p18 + NDIVR;		p18 += ( (p18 >> DAT_BITS) << PAD_BITS );
		p1a = p19 + NDIVR;		p19 += ( (p19 >> DAT_BITS) << PAD_BITS );
		p1b = p1a + NDIVR;		p1a += ( (p1a >> DAT_BITS) << PAD_BITS );
		p1c = p1b + NDIVR;		p1b += ( (p1b >> DAT_BITS) << PAD_BITS );
		p1d = p1c + NDIVR;		p1c += ( (p1c >> DAT_BITS) << PAD_BITS );
		p1e = p1d + NDIVR;		p1d += ( (p1d >> DAT_BITS) << PAD_BITS );
		p1f = p1e + NDIVR;		p1e += ( (p1e >> DAT_BITS) << PAD_BITS );
								p1f += ( (p1f >> DAT_BITS) << PAD_BITS );
		NDIVR <<= 4;	p20 = NDIVR << 1;	// NDIVR holds unpadded p10
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
		p100= pf0 + NDIVR;		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
		p110= p100 + NDIVR;		p100+= ( (p100>> DAT_BITS) << PAD_BITS );
		p120= p110 + NDIVR;		p110+= ( (p110>> DAT_BITS) << PAD_BITS );
		p130= p120 + NDIVR;		p120+= ( (p120>> DAT_BITS) << PAD_BITS );
		p140= p130 + NDIVR;		p130+= ( (p130>> DAT_BITS) << PAD_BITS );
		p150= p140 + NDIVR;		p140+= ( (p140>> DAT_BITS) << PAD_BITS );
		p160= p150 + NDIVR;		p150+= ( (p150>> DAT_BITS) << PAD_BITS );
		p170= p160 + NDIVR;		p160+= ( (p160>> DAT_BITS) << PAD_BITS );
		p180= p170 + NDIVR;		p170+= ( (p170>> DAT_BITS) << PAD_BITS );
		p190= p180 + NDIVR;		p180+= ( (p180>> DAT_BITS) << PAD_BITS );
		p1a0= p190 + NDIVR;		p190+= ( (p190>> DAT_BITS) << PAD_BITS );
		p1b0= p1a0 + NDIVR;		p1a0+= ( (p1a0>> DAT_BITS) << PAD_BITS );
		p1c0= p1b0 + NDIVR;		p1b0+= ( (p1b0>> DAT_BITS) << PAD_BITS );
		p1d0= p1c0 + NDIVR;		p1c0+= ( (p1c0>> DAT_BITS) << PAD_BITS );
		p1e0= p1d0 + NDIVR;		p1d0+= ( (p1d0>> DAT_BITS) << PAD_BITS );
		p1f0= p1e0 + NDIVR;		p1e0+= ( (p1e0>> DAT_BITS) << PAD_BITS );
		p200= p1f0 + NDIVR;		p1f0+= ( (p1f0>> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			p200+= ( (p200>> DAT_BITS) << PAD_BITS );
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
			poff[192+l] = poff[l] + p100+p200;
		}

	// DIF:
		// Set array offsets for radix-64 DFT in/outputs:
		// set 1 is w.r.to: a-array:										// set 2 w.r.to t-array - These are shared by DIT:
		dif_i_offsets[0x00] = 0   ;		dif_i_offsets[0x20] = 0   +p200;	o_offsets[0x00] = 0x00<<1;	o_offsets[0x20] = 0x20<<1;
		dif_i_offsets[0x01] = p10 ;		dif_i_offsets[0x21] = p10 +p200;	o_offsets[0x01] = 0x01<<1;	o_offsets[0x21] = 0x21<<1;
		dif_i_offsets[0x02] = p20 ;		dif_i_offsets[0x22] = p20 +p200;	o_offsets[0x02] = 0x02<<1;	o_offsets[0x22] = 0x22<<1;
		dif_i_offsets[0x03] = p30 ;		dif_i_offsets[0x23] = p30 +p200;	o_offsets[0x03] = 0x03<<1;	o_offsets[0x23] = 0x23<<1;
		dif_i_offsets[0x04] = p40 ;		dif_i_offsets[0x24] = p40 +p200;	o_offsets[0x04] = 0x04<<1;	o_offsets[0x24] = 0x24<<1;
		dif_i_offsets[0x05] = p50 ;		dif_i_offsets[0x25] = p50 +p200;	o_offsets[0x05] = 0x05<<1;	o_offsets[0x25] = 0x25<<1;
		dif_i_offsets[0x06] = p60 ;		dif_i_offsets[0x26] = p60 +p200;	o_offsets[0x06] = 0x06<<1;	o_offsets[0x26] = 0x26<<1;
		dif_i_offsets[0x07] = p70 ;		dif_i_offsets[0x27] = p70 +p200;	o_offsets[0x07] = 0x07<<1;	o_offsets[0x27] = 0x27<<1;
		dif_i_offsets[0x08] = p80 ;		dif_i_offsets[0x28] = p80 +p200;	o_offsets[0x08] = 0x08<<1;	o_offsets[0x28] = 0x28<<1;
		dif_i_offsets[0x09] = p90 ;		dif_i_offsets[0x29] = p90 +p200;	o_offsets[0x09] = 0x09<<1;	o_offsets[0x29] = 0x29<<1;
		dif_i_offsets[0x0a] = pa0 ;		dif_i_offsets[0x2a] = pa0 +p200;	o_offsets[0x0a] = 0x0a<<1;	o_offsets[0x2a] = 0x2a<<1;
		dif_i_offsets[0x0b] = pb0 ;		dif_i_offsets[0x2b] = pb0 +p200;	o_offsets[0x0b] = 0x0b<<1;	o_offsets[0x2b] = 0x2b<<1;
		dif_i_offsets[0x0c] = pc0 ;		dif_i_offsets[0x2c] = pc0 +p200;	o_offsets[0x0c] = 0x0c<<1;	o_offsets[0x2c] = 0x2c<<1;
		dif_i_offsets[0x0d] = pd0 ;		dif_i_offsets[0x2d] = pd0 +p200;	o_offsets[0x0d] = 0x0d<<1;	o_offsets[0x2d] = 0x2d<<1;
		dif_i_offsets[0x0e] = pe0 ;		dif_i_offsets[0x2e] = pe0 +p200;	o_offsets[0x0e] = 0x0e<<1;	o_offsets[0x2e] = 0x2e<<1;
		dif_i_offsets[0x0f] = pf0 ;		dif_i_offsets[0x2f] = pf0 +p200;	o_offsets[0x0f] = 0x0f<<1;	o_offsets[0x2f] = 0x2f<<1;
		dif_i_offsets[0x10] = p100;		dif_i_offsets[0x30] = p100+p200;	o_offsets[0x10] = 0x10<<1;	o_offsets[0x30] = 0x30<<1;
		dif_i_offsets[0x11] = p110;		dif_i_offsets[0x31] = p110+p200;	o_offsets[0x11] = 0x11<<1;	o_offsets[0x31] = 0x31<<1;
		dif_i_offsets[0x12] = p120;		dif_i_offsets[0x32] = p120+p200;	o_offsets[0x12] = 0x12<<1;	o_offsets[0x32] = 0x32<<1;
		dif_i_offsets[0x13] = p130;		dif_i_offsets[0x33] = p130+p200;	o_offsets[0x13] = 0x13<<1;	o_offsets[0x33] = 0x33<<1;
		dif_i_offsets[0x14] = p140;		dif_i_offsets[0x34] = p140+p200;	o_offsets[0x14] = 0x14<<1;	o_offsets[0x34] = 0x34<<1;
		dif_i_offsets[0x15] = p150;		dif_i_offsets[0x35] = p150+p200;	o_offsets[0x15] = 0x15<<1;	o_offsets[0x35] = 0x35<<1;
		dif_i_offsets[0x16] = p160;		dif_i_offsets[0x36] = p160+p200;	o_offsets[0x16] = 0x16<<1;	o_offsets[0x36] = 0x36<<1;
		dif_i_offsets[0x17] = p170;		dif_i_offsets[0x37] = p170+p200;	o_offsets[0x17] = 0x17<<1;	o_offsets[0x37] = 0x37<<1;
		dif_i_offsets[0x18] = p180;		dif_i_offsets[0x38] = p180+p200;	o_offsets[0x18] = 0x18<<1;	o_offsets[0x38] = 0x38<<1;
		dif_i_offsets[0x19] = p190;		dif_i_offsets[0x39] = p190+p200;	o_offsets[0x19] = 0x19<<1;	o_offsets[0x39] = 0x39<<1;
		dif_i_offsets[0x1a] = p1a0;		dif_i_offsets[0x3a] = p1a0+p200;	o_offsets[0x1a] = 0x1a<<1;	o_offsets[0x3a] = 0x3a<<1;
		dif_i_offsets[0x1b] = p1b0;		dif_i_offsets[0x3b] = p1b0+p200;	o_offsets[0x1b] = 0x1b<<1;	o_offsets[0x3b] = 0x3b<<1;
		dif_i_offsets[0x1c] = p1c0;		dif_i_offsets[0x3c] = p1c0+p200;	o_offsets[0x1c] = 0x1c<<1;	o_offsets[0x3c] = 0x3c<<1;
		dif_i_offsets[0x1d] = p1d0;		dif_i_offsets[0x3d] = p1d0+p200;	o_offsets[0x1d] = 0x1d<<1;	o_offsets[0x3d] = 0x3d<<1;
		dif_i_offsets[0x1e] = p1e0;		dif_i_offsets[0x3e] = p1e0+p200;	o_offsets[0x1e] = 0x1e<<1;	o_offsets[0x3e] = 0x3e<<1;
		dif_i_offsets[0x1f] = p1f0;		dif_i_offsets[0x3f] = p1f0+p200;	o_offsets[0x1f] = 0x1f<<1;	o_offsets[0x3f] = 0x3f<<1;
	#ifdef USE_SSE2
		// o_offs 2x,4x,... what might be expected (= 2 vec_dbl per output) due to cast-to-double-pointer of DIT-64 output arg
		for(l = 0; l < 64; l++) {
			dif_o_offsets[l] = o_offsets[l] << (L2_SZ_VD - 3);	// 2x for sse2, 4x for avx, etc
		}
	#endif

		dif_po_br[0x0] =  0; dif_po_br[0x1] = p8; dif_po_br[0x2] = p4; dif_po_br[0x3] = pc;
		dif_po_br[0x4] = p2; dif_po_br[0x5] = pa; dif_po_br[0x6] = p6; dif_po_br[0x7] = pe;
		dif_po_br[0x8] = p1; dif_po_br[0x9] = p9; dif_po_br[0xa] = p5; dif_po_br[0xb] = pd;
		dif_po_br[0xc] = p3; dif_po_br[0xd] = pb; dif_po_br[0xe] = p7; dif_po_br[0xf] = pf;

	// DIT:
		// Set array offsets for radix-64 DFT in/outputs:
		// set 1 is w.r.to: a-array:
		dit_i_offsets[0x00] =  0;		dit_i_offsets[0x20] =  0+p20;
		dit_i_offsets[0x01] = p1;		dit_i_offsets[0x21] = p1+p20;
		dit_i_offsets[0x02] = p2;		dit_i_offsets[0x22] = p2+p20;
		dit_i_offsets[0x03] = p3;		dit_i_offsets[0x23] = p3+p20;
		dit_i_offsets[0x04] = p4;		dit_i_offsets[0x24] = p4+p20;
		dit_i_offsets[0x05] = p5;		dit_i_offsets[0x25] = p5+p20;
		dit_i_offsets[0x06] = p6;		dit_i_offsets[0x26] = p6+p20;
		dit_i_offsets[0x07] = p7;		dit_i_offsets[0x27] = p7+p20;
		dit_i_offsets[0x08] = p8;		dit_i_offsets[0x28] = p8+p20;
		dit_i_offsets[0x09] = p9;		dit_i_offsets[0x29] = p9+p20;
		dit_i_offsets[0x0a] = pa;		dit_i_offsets[0x2a] = pa+p20;
		dit_i_offsets[0x0b] = pb;		dit_i_offsets[0x2b] = pb+p20;
		dit_i_offsets[0x0c] = pc;		dit_i_offsets[0x2c] = pc+p20;
		dit_i_offsets[0x0d] = pd;		dit_i_offsets[0x2d] = pd+p20;
		dit_i_offsets[0x0e] = pe;		dit_i_offsets[0x2e] = pe+p20;
		dit_i_offsets[0x0f] = pf;		dit_i_offsets[0x2f] = pf+p20;
		dit_i_offsets[0x10] = p10;		dit_i_offsets[0x30] = p10+p20;
		dit_i_offsets[0x11] = p11;		dit_i_offsets[0x31] = p11+p20;
		dit_i_offsets[0x12] = p12;		dit_i_offsets[0x32] = p12+p20;
		dit_i_offsets[0x13] = p13;		dit_i_offsets[0x33] = p13+p20;
		dit_i_offsets[0x14] = p14;		dit_i_offsets[0x34] = p14+p20;
		dit_i_offsets[0x15] = p15;		dit_i_offsets[0x35] = p15+p20;
		dit_i_offsets[0x16] = p16;		dit_i_offsets[0x36] = p16+p20;
		dit_i_offsets[0x17] = p17;		dit_i_offsets[0x37] = p17+p20;
		dit_i_offsets[0x18] = p18;		dit_i_offsets[0x38] = p18+p20;
		dit_i_offsets[0x19] = p19;		dit_i_offsets[0x39] = p19+p20;
		dit_i_offsets[0x1a] = p1a;		dit_i_offsets[0x3a] = p1a+p20;
		dit_i_offsets[0x1b] = p1b;		dit_i_offsets[0x3b] = p1b+p20;
		dit_i_offsets[0x1c] = p1c;		dit_i_offsets[0x3c] = p1c+p20;
		dit_i_offsets[0x1d] = p1d;		dit_i_offsets[0x3d] = p1d+p20;
		dit_i_offsets[0x1e] = p1e;		dit_i_offsets[0x3e] = p1e+p20;
		dit_i_offsets[0x1f] = p1f;		dit_i_offsets[0x3f] = p1f+p20;

		dit_po_br[4*0x0] =  0; dit_po_br[4*0x1] = p8; dit_po_br[4*0x2] = p4; dit_po_br[4*0x3] = pc;
		dit_po_br[4*0x4] = p2; dit_po_br[4*0x5] = pa; dit_po_br[4*0x6] = p6; dit_po_br[4*0x7] = pe;
		dit_po_br[4*0x8] = p1; dit_po_br[4*0x9] = p9; dit_po_br[4*0xa] = p5; dit_po_br[4*0xb] = pd;
		dit_po_br[4*0xc] = p3; dit_po_br[4*0xd] = pb; dit_po_br[4*0xe] = p7; dit_po_br[4*0xf] = pf;
		// Each of the foregoing 16 indices is head of a (i0,i0+p20,i0+p10,i0+p30) quartet:
		for(l = 0; l < 16; l++) {
			j = l << 2;
			dit_po_br[j+1] = dit_po_br[j] + p20;
			dit_po_br[j+2] = dit_po_br[j] + p10;
			dit_po_br[j+3] = dit_po_br[j] + p30;
		}
		dit_poffs[0x0] =           0; dit_poffs[0x1] =         p40; dit_poffs[0x2] =         p80; dit_poffs[0x3] =         pc0;
		dit_poffs[0x4] =        p100; dit_poffs[0x5] =        p140; dit_poffs[0x6] =        p180; dit_poffs[0x7] =        p1c0;
		dit_poffs[0x8] =        p200; dit_poffs[0x9] =  p40 + p200; dit_poffs[0xa] =  p80 + p200; dit_poffs[0xb] =  pc0 + p200;
		dit_poffs[0xc] = p100 + p200; dit_poffs[0xd] = p140 + p200; dit_poffs[0xe] = p180 + p200; dit_poffs[0xf] = p1c0 + p200;

		if(_cy_r[0])	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;
			for(i = 0; i < RADIX; i++) {
				free((void *)_bjmodn[i]); _bjmodn[i] = 0x0;
				free((void *)  _cy_r[i]);   _cy_r[i] = 0x0;
				free((void *)  _cy_i[i]);   _cy_i[i] = 0x0;
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
			_cy_r[i]	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r[i]== 0x0);
			_cy_i[i]	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i[i]== 0x0);
		}

		ASSERT(ptr_prod == 0, "ERROR: unable to allocate one or more auxiliary arrays!");

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
			i.e. the one that n2/radix-separated FFT outputs need:
			*/
			_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"ERROR: unable to allocate array _bjmodnini.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
			_bjmodnini[0] = 0;
			_bjmodnini[1] = 0;
			for(j=0; j < NDIVR/CY_THREADS; j++)
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
			for(j=0; j < NDIVR; j++)
			{
				bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
			}
			ASSERT(_bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");
			if(CY_THREADS > 1)
			{
				for(ithread = 1; ithread < CY_THREADS; ithread++)
				{
					_i[ithread] = ((uint32)(sw - _bjmodnini[ithread]) >> 31);
				}
			}

			// Include 0-thread here ... bjmodn terms all 0 for that, but need jhi computed for all threads:
			j = _bjmodnini[CY_THREADS];
			for(ithread = 0; ithread < CY_THREADS; ithread++)
			{
				_bjmodn[0][ithread] = _bjmodnini[ithread];
				for(i = 1; i < RADIX; i++) {
					MOD_ADD32(_bjmodn[i-1][ithread], j, n, _bjmodn[i][ithread]);
				}
			}
		}

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
		  // For pow2 FFT lengths, bjmodn only needed for mers-mod:
		  if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		  {
			tdat[ithread].bjmodnini = _bjmodnini[CY_THREADS];
			tdat[ithread].bjmodn0 = _bjmodnini[ithread];
		  }
		#ifdef USE_SSE2
			tdat[ithread].r00      = __r0 + ithread*cslots_in_local_store;
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

/*...The radix-1024 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(i = 0; i < RADIX; i++) {
			_cy_r[i][ithread] = 0;
			_cy_i[i][ithread] = 0;
		}
	}
  #if 0	//ndef USE_SSE2	*** v20: Non-SIMD builds now also support shifted-residue
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy_r[0][0] = -2;
	}
  #endif
	*fracmax=0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	// init inverse-weight scale factor = 2/n for normal carry pass, 1 for wrapper pass

for(outer=0; outer <= 1; outer++)
{
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		/*
		Moved this inside the outer-loop, so on cleanup pass can use it to reset _col,_co2,_co3 starting values,
		then simply overwrite it with 1 prior to starting the k-loop.
		*/
		_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (_i[0] = 1).	*/

		khi = n_div_nwt/CY_THREADS;
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_jstart[ithread] = ithread*NDIVR/CY_THREADS;
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + jhi_wrap_mers;	/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + nwt-1;

			_col[ithread] = ithread*(khi*RADIX);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
			_co2[ithread] = (n>>nwt_bits)-1+RADIX - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
			_co3[ithread] = _co2[ithread]-RADIX;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
		}
	}
	else
	{
		khi = 1;
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_jstart[ithread] = ithread*NDIVR/CY_THREADS;
			/*
			For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
			*/
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + jhi_wrap_ferm;	/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + n_div_nwt/CY_THREADS;
		}
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
		ASSERT(tdat[ithread].rn0 == rn0, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].rn1 == rn1, "thread-local memcheck fail!");
	#ifdef USE_SSE2
		ASSERT(tdat[ithread].r00 == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].r00;
		ASSERT(((tmp + 0x1080)->d0 == ISRT2 && (tmp + 0x1080)->d1 == ISRT2), "thread-local memcheck failed!");
		tmp = tdat[ithread].half_arr;
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		ASSERT(((tmp-1)->d0 == base[0] && (tmp-1)->d1 == baseinv[1] && (tmp-1)->d2 == wts_mult[1] && (tmp-1)->d3 == inv_mult[0]), "thread-local memcheck failed!");
	  #else
		ASSERT(((tmp-1)->d0 == crnd && (tmp-1)->d1 == crnd), "thread-local memcheck failed!");
	  #endif
	#endif

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX512
			/* No-Op */
		#elif defined(USE_AVX)
			// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
			dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#elif defined(USE_SSE2)
			dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#endif
			/* init carries	*/
			for(i = 0; i < RADIX; i++) {
				tdat[ithread].cy_r[i] = _cy_r[i][ithread];
			}
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
		#ifdef USE_AVX512
			/* No-Op */
		#elif defined(USE_SSE2)
			// This is slightly different for power-of-2 DFTs: Here, scale is in the +2 slot, base & baseinv remain fixed in 0,+1 slots:
			dtmp = tmp->d0 * (tmp+1)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = tmp->d1 * (tmp+1)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#endif
			// scale gets set immediately prior to calling carry macro, hence no use checking it here.
			/* init carries	*/
			for(i = 0; i < RADIX; i++) {
				tdat[ithread].cy_r[i] = _cy_r[i][ithread];
				tdat[ithread].cy_i[i] = _cy_i[i][ithread];
			}
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

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			col = _col[ithread];
			co2 = _co2[ithread];
			co3 = _co3[ithread];

			for(l = 0; l < RADIX; l++) {
				bjmodn[l] = _bjmodn[l][ithread];
			}
			/* init carries	*/
		#ifdef USE_AVX512
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 8, ++tmp) {
				tmp->d0 = _cy_r[l  ][ithread];
				tmp->d1 = _cy_r[l+1][ithread];
				tmp->d2 = _cy_r[l+2][ithread];
				tmp->d3 = _cy_r[l+3][ithread];
				tmp->d4 = _cy_r[l+4][ithread];
				tmp->d5 = _cy_r[l+5][ithread];
				tmp->d6 = _cy_r[l+6][ithread];
				tmp->d7 = _cy_r[l+7][ithread];
			}
		#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 4, ++tmp) {
				tmp->d0 = _cy_r[l  ][ithread];
				tmp->d1 = _cy_r[l+1][ithread];
				tmp->d2 = _cy_r[l+2][ithread];
				tmp->d3 = _cy_r[l+3][ithread];
			}
		#elif defined(USE_SSE2)
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 2, ++tmp) {
				tmp->d0 = _cy_r[l  ][ithread];
				tmp->d1 = _cy_r[l+1][ithread];
			}
		#else
			for(l = 0; l < RADIX; l++) {
				cy_r[l] = _cy_r[l][ithread];
			}
		#endif
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
			/* init carries	*/
		#ifdef USE_AVX512
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 8, ++tmp, ++tm2) {
				tmp->d0 = _cy_r[l  ][ithread];		tm2->d0 = _cy_i[l  ][ithread];
				tmp->d1 = _cy_r[l+1][ithread];		tm2->d1 = _cy_i[l+1][ithread];
				tmp->d2 = _cy_r[l+2][ithread];		tm2->d2 = _cy_i[l+2][ithread];
				tmp->d3 = _cy_r[l+3][ithread];		tm2->d3 = _cy_i[l+3][ithread];
				tmp->d4 = _cy_r[l+4][ithread];		tm2->d4 = _cy_i[l+4][ithread];
				tmp->d5 = _cy_r[l+5][ithread];		tm2->d5 = _cy_i[l+5][ithread];
				tmp->d6 = _cy_r[l+6][ithread];		tm2->d6 = _cy_i[l+6][ithread];
				tmp->d7 = _cy_r[l+7][ithread];		tm2->d7 = _cy_i[l+7][ithread];
			}
		#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 4, ++tmp, ++tm2) {
				tmp->d0 = _cy_r[l  ][ithread];		tm2->d0 = _cy_i[l  ][ithread];
				tmp->d1 = _cy_r[l+1][ithread];		tm2->d1 = _cy_i[l+1][ithread];
				tmp->d2 = _cy_r[l+2][ithread];		tm2->d2 = _cy_i[l+2][ithread];
				tmp->d3 = _cy_r[l+3][ithread];		tm2->d3 = _cy_i[l+3][ithread];
			}
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			tmp = cy_r;
			for(l = 0; l < RADIX; l++, ++tmp) {
				// This relies on the cy_R,i sections of the SIMD data being contiguous, i.e.
				// step-thru the cy_r data via the tmp-pointer takes us seamlessly into the cy_i:
				tmp->d0 = _cy_r[l][ithread];		tmp->d1 = _cy_i[l][ithread];
			}
		#else
			for(l = 0; l < RADIX; l++) {
				cy_r[l] = _cy_r[l][ithread];		cy_i[l] = _cy_i[l][ithread];
			}
		#endif
		}

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix1024_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX512
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 8, ++tmp) {
				_cy_r[l  ][ithread] = tmp->d0;
				_cy_r[l+1][ithread] = tmp->d1;
				_cy_r[l+2][ithread] = tmp->d2;
				_cy_r[l+3][ithread] = tmp->d3;
				_cy_r[l+4][ithread] = tmp->d4;
				_cy_r[l+5][ithread] = tmp->d5;
				_cy_r[l+6][ithread] = tmp->d6;
				_cy_r[l+7][ithread] = tmp->d7;
			}
			if(full_pass) {
				t0 = MAX(max_err->d0,max_err->d1);
				t1 = MAX(max_err->d2,max_err->d3);
				t2 = MAX(max_err->d4,max_err->d5);
				t3 = MAX(max_err->d6,max_err->d7);
				maxerr = MAX( MAX(t0,t1), MAX(t2,t3) );
			}
		#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 4, ++tmp) {
				_cy_r[l  ][ithread] = tmp->d0;
				_cy_r[l+1][ithread] = tmp->d1;
				_cy_r[l+2][ithread] = tmp->d2;
				_cy_r[l+3][ithread] = tmp->d3;
			}
			if(full_pass) maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 2, ++tmp) {
				_cy_r[l  ][ithread] = tmp->d0;
				_cy_r[l+1][ithread] = tmp->d1;
			}
			if(full_pass) maxerr = MAX(max_err->d0,max_err->d1);
		#else
			for(l = 0; l < RADIX; l++) {
				_cy_r[l][ithread] = cy_r[l];
			}
		#endif
		}
		else
		{
		#ifdef USE_AVX512
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 8, ++tmp, ++tm2) {
				_cy_r[l  ][ithread] = tmp->d0;		_cy_i[l  ][ithread] = tm2->d0;
				_cy_r[l+1][ithread] = tmp->d1;		_cy_i[l+1][ithread] = tm2->d1;
				_cy_r[l+2][ithread] = tmp->d2;		_cy_i[l+2][ithread] = tm2->d2;
				_cy_r[l+3][ithread] = tmp->d3;		_cy_i[l+3][ithread] = tm2->d3;
				_cy_r[l+4][ithread] = tmp->d4;		_cy_i[l+4][ithread] = tm2->d4;
				_cy_r[l+5][ithread] = tmp->d5;		_cy_i[l+5][ithread] = tm2->d5;
				_cy_r[l+6][ithread] = tmp->d6;		_cy_i[l+6][ithread] = tm2->d6;
				_cy_r[l+7][ithread] = tmp->d7;		_cy_i[l+7][ithread] = tm2->d7;
			}
			if(full_pass) {
				t0 = MAX(max_err->d0,max_err->d1);
				t1 = MAX(max_err->d2,max_err->d3);
				t2 = MAX(max_err->d4,max_err->d5);
				t3 = MAX(max_err->d6,max_err->d7);
				maxerr = MAX( MAX(t0,t1), MAX(t2,t3) );
			}
		#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 4, ++tmp, ++tm2) {
				_cy_r[l  ][ithread] = tmp->d0;		_cy_i[l  ][ithread] = tm2->d0;
				_cy_r[l+1][ithread] = tmp->d1;		_cy_i[l+1][ithread] = tm2->d1;
				_cy_r[l+2][ithread] = tmp->d2;		_cy_i[l+2][ithread] = tm2->d2;
				_cy_r[l+3][ithread] = tmp->d3;		_cy_i[l+3][ithread] = tm2->d3;
			}
			if(full_pass) maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			tmp = cy_r;
			for(l = 0; l < RADIX; l++, ++tmp) {
				// This relies on the cy_R,i sections of the SIMD data being contiguous, i.e.
				// step-thru the cy_r data via the tmp-pointer takes us seamlessly into the cy_i:
				_cy_r[l][ithread] = tmp->d0;		_cy_i[l][ithread] = tmp->d1;
			}
			if(full_pass) maxerr = MAX(max_err->d0,max_err->d1);
		#else
			for(l = 0; l < RADIX; l++) {
				_cy_r[l][ithread] = cy_r[l];		_cy_i[l][ithread] = cy_i[l];
			}
		#endif
		}

  #endif	// #ifdef USE_PTHREAD

	}	/******* END OF PARALLEL FOR-LOOP ********/

#ifdef USE_PTHREAD	// End of threadpool-based dispatch: Add a small wait-loop to ensure all threads complete

  #if 0//def OS_TYPE_MACOSX

	/*** Main execution thread executes remaining chunks in serial fashion (but in || with the pool threads): ***/
	for(j = 0; j < main_work_units; ++j)
	{
	//	printf("adding main task %d\n",j + pool_work_units);
		ASSERT(0x0 == cy1024_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
	ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
	ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(0 == mlucas_nanosleep(&ns_time), "nanosleep fail!");
	}
//	printf("%s end  ; #tasks = %d, #free_tasks = %d\n",func, tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		if(maxerr < tdat[ithread].maxerr) {
			maxerr = tdat[ithread].maxerr;
		}
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			for(l = 0; l < RADIX; l++) {
				_cy_r[l][ithread] = tdat[ithread].cy_r[l];
			}
		}
		else
		{
			for(l = 0; l < RADIX; l++) {
				_cy_r[l][ithread] = tdat[ithread].cy_r[l];
				_cy_i[l][ithread] = tdat[ithread].cy_i[l];
			}
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
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		for(l = 0; l < RADIX; l++) {
			t[l].re = _cy_r[l][CY_THREADS - 1];
		}
		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			for(l = 0; l < RADIX; l++) {
				_cy_r[l][ithread] = _cy_r[l][ithread-1];
			}
		}
		_cy_r[0][0] =+t[RADIX-1].re;	/* ...The wraparound carry is here: */
		for(l = 1; l < RADIX; l++) {
			_cy_r[l][0] = t[l-1].re;
		}
	}
	else
	{
		j = CY_THREADS - 1;
		for(l = 0; l < RADIX; l++) {
			t[l].re = _cy_r[l][j];		t[l].im = _cy_i[l][j];
		}
		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			for(l = 0; l < RADIX; l++) {
				_cy_r[l][ithread] = _cy_r[l][ithread-1];	_cy_i[l][ithread] = _cy_i[l][ithread-1];
			}
		}
		_cy_r[0][0] =-t[RADIX-1].im;	_cy_i[0][0] =+t[RADIX-1].re;	/* ...The 2 Mo"bius carries are here: */
		for(l = 1; l < RADIX; l++) {
			_cy_r[l][0] = t[l-1].re;	_cy_i[l][0] = t[l-1].im;
		}
	}

	full_pass = 0;
	scale = prp_mult = 1;

	/*
	For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
	*/
	if((MODULUS_TYPE == MODULUS_TYPE_GENFFTMUL) || (TRANSFORM_TYPE == RIGHT_ANGLE))
		j_jhi = jhi_wrap_ferm;
	else
		j_jhi = jhi_wrap_mers;

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(j = ithread*pini; j <= ithread*pini + j_jhi; j++)
		{
			// Generate padded version of j, since prepadding pini is thread-count unsafe:
			j1 = j + ( (j >> DAT_BITS) << PAD_BITS );
			for(ntmp = 0; ntmp < RADIX>>2; ntmp++) {
				jt = j1 + poff[ntmp];
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
			dtmp += fabs(_cy_r[l][ithread]) + fabs(_cy_i[l][ithread]);
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

/****************/

void radix1024_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-1024 complex DIF FFT pass on the data in the length-N real vector A.
!
!   The data are stored in a 1-D zero-offset array, with 2^PAD_BITS 8-byte padding elements inserted
!   between every block of 2^DAT_BITS contiguous data. The array padding is to prevent data being accessed
!   in strides that are large powers of two and thus to minimize cache thrashing
!   (in cache-based microprocessor architectures) or bank conflicts (in supercomputers.)
!
!   See the documentation in radix64_dif_pass for further details on storage, indexing and DFT decomposition.
*/
	int i,j,j1,j2,jt,jp;
	static int NDIVR,first_entry=TRUE,
			p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
		p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p1a,p1b,p1c,p1d,p1e,p1f,
  #if RAD1 == 32
		p20,p40,p60,p80,pa0,pc0,pe0,p100,p120,p140,p160,p180,p1a0,p1c0,p1e0,p200;
	static int i_offsets[32], o_offsets[32], p_offsets[32], q_offsets[32];
  #else
		p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,p100,p110,p120,p130,p140,p150,p160,p170,p180,p190,p1a0,p1b0,p1c0,p1d0,p1e0,p1f0,p200;
	static int i_offsets[64], o_offsets[64];
	static int po_br[16];
	// We prefer pointer-based array-element access, because that allows our radix16 DFT-with-twiddles
	// to look the same in terms of array-element arglists:
	const double *addr,*addi;
	struct complex *tptr;
	#include "radix1024_twiddles.h"
  #endif
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex t[RADIX];

// Use to auto-gen the twiddles:
//	print_pow2_twiddles(1024, 64,16);
//	exit(0);

	// New runlength?
	if(!first_entry && (n >> 10) != NDIVR)
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		ASSERT((double *)t == &(t[0].re), "Unexpected value for Tmp-array-start pointer!");
		first_entry=FALSE;
		NDIVR = n >> 10;
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
		p10 = pf + NDIVR;		pf += ( (pf >> DAT_BITS) << PAD_BITS );
		p11 = p10 + NDIVR;		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p12 = p11 + NDIVR;		p11 += ( (p11 >> DAT_BITS) << PAD_BITS );
		p13 = p12 + NDIVR;		p12 += ( (p12 >> DAT_BITS) << PAD_BITS );
		p14 = p13 + NDIVR;		p13 += ( (p13 >> DAT_BITS) << PAD_BITS );
		p15 = p14 + NDIVR;		p14 += ( (p14 >> DAT_BITS) << PAD_BITS );
		p16 = p15 + NDIVR;		p15 += ( (p15 >> DAT_BITS) << PAD_BITS );
		p17 = p16 + NDIVR;		p16 += ( (p16 >> DAT_BITS) << PAD_BITS );
		p18 = p17 + NDIVR;		p17 += ( (p17 >> DAT_BITS) << PAD_BITS );
		p19 = p18 + NDIVR;		p18 += ( (p18 >> DAT_BITS) << PAD_BITS );
		p1a = p19 + NDIVR;		p19 += ( (p19 >> DAT_BITS) << PAD_BITS );
		p1b = p1a + NDIVR;		p1a += ( (p1a >> DAT_BITS) << PAD_BITS );
		p1c = p1b + NDIVR;		p1b += ( (p1b >> DAT_BITS) << PAD_BITS );
		p1d = p1c + NDIVR;		p1c += ( (p1c >> DAT_BITS) << PAD_BITS );
		p1e = p1d + NDIVR;		p1d += ( (p1d >> DAT_BITS) << PAD_BITS );
		p1f = p1e + NDIVR;		p1e += ( (p1e >> DAT_BITS) << PAD_BITS );
								p1f += ( (p1f >> DAT_BITS) << PAD_BITS );
	#if RAD1 == 32
		NDIVR <<= 5;	p20 = NDIVR;	// NDIVR holds unpadded p20
		p40 = p20 + NDIVR;		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p60 = p40 + NDIVR;		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p80 = p60 + NDIVR;		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		pa0 = p80 + NDIVR;		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		pc0 = pa0 + NDIVR;		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pe0 = pc0 + NDIVR;		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		p100= pe0 + NDIVR;		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		p120= p100+ NDIVR;		p100+= ( (p100>> DAT_BITS) << PAD_BITS );
		p140= p120+ NDIVR;		p120+= ( (p120>> DAT_BITS) << PAD_BITS );
		p160= p140+ NDIVR;		p140+= ( (p140>> DAT_BITS) << PAD_BITS );
		p180= p160+ NDIVR;		p160+= ( (p160>> DAT_BITS) << PAD_BITS );
		p1a0= p180+ NDIVR;		p180+= ( (p180>> DAT_BITS) << PAD_BITS );
		p1c0= p1a0+ NDIVR;		p1a0+= ( (p1a0>> DAT_BITS) << PAD_BITS );
		p1e0= p1c0+ NDIVR;		p1c0+= ( (p1c0>> DAT_BITS) << PAD_BITS );
		p200= p1e0+ NDIVR;		p1e0+= ( (p1e0>> DAT_BITS) << PAD_BITS );
		NDIVR >>= 5;			p200+= ( (p200>> DAT_BITS) << PAD_BITS );
	#else
		NDIVR <<= 4;	p20 = NDIVR << 1;	// NDIVR holds unpadded p10
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
		p100= pf0 + NDIVR;		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
		p110= p100 + NDIVR;		p100+= ( (p100>> DAT_BITS) << PAD_BITS );
		p120= p110 + NDIVR;		p110+= ( (p110>> DAT_BITS) << PAD_BITS );
		p130= p120 + NDIVR;		p120+= ( (p120>> DAT_BITS) << PAD_BITS );
		p140= p130 + NDIVR;		p130+= ( (p130>> DAT_BITS) << PAD_BITS );
		p150= p140 + NDIVR;		p140+= ( (p140>> DAT_BITS) << PAD_BITS );
		p160= p150 + NDIVR;		p150+= ( (p150>> DAT_BITS) << PAD_BITS );
		p170= p160 + NDIVR;		p160+= ( (p160>> DAT_BITS) << PAD_BITS );
		p180= p170 + NDIVR;		p170+= ( (p170>> DAT_BITS) << PAD_BITS );
		p190= p180 + NDIVR;		p180+= ( (p180>> DAT_BITS) << PAD_BITS );
		p1a0= p190 + NDIVR;		p190+= ( (p190>> DAT_BITS) << PAD_BITS );
		p1b0= p1a0 + NDIVR;		p1a0+= ( (p1a0>> DAT_BITS) << PAD_BITS );
		p1c0= p1b0 + NDIVR;		p1b0+= ( (p1b0>> DAT_BITS) << PAD_BITS );
		p1d0= p1c0 + NDIVR;		p1c0+= ( (p1c0>> DAT_BITS) << PAD_BITS );
		p1e0= p1d0 + NDIVR;		p1d0+= ( (p1d0>> DAT_BITS) << PAD_BITS );
		p1f0= p1e0 + NDIVR;		p1e0+= ( (p1e0>> DAT_BITS) << PAD_BITS );
		p200= p1f0 + NDIVR;		p1f0+= ( (p1f0>> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			p200+= ( (p200>> DAT_BITS) << PAD_BITS );
	#endif

	#if RAD1 == 32
		// Set array offsets for radix-32 DFT in/outputs:
		// set 1 is w.r.to: a-array:								// set 2 w.r.to t-array:
		i_offsets[0x00] = 0   ;		i_offsets[0x10] = 0   +p200;	o_offsets[0x00] = 0x00<<1;	o_offsets[0x10] = 0x10<<1;
		i_offsets[0x01] = p20 ;		i_offsets[0x11] = p20 +p200;	o_offsets[0x01] = 0x01<<1;	o_offsets[0x11] = 0x11<<1;
		i_offsets[0x02] = p40 ;		i_offsets[0x12] = p40 +p200;	o_offsets[0x02] = 0x02<<1;	o_offsets[0x12] = 0x12<<1;
		i_offsets[0x03] = p60 ;		i_offsets[0x13] = p60 +p200;	o_offsets[0x03] = 0x03<<1;	o_offsets[0x13] = 0x13<<1;
		i_offsets[0x04] = p80 ;		i_offsets[0x14] = p80 +p200;	o_offsets[0x04] = 0x04<<1;	o_offsets[0x14] = 0x14<<1;
		i_offsets[0x05] = pa0 ;		i_offsets[0x15] = pa0 +p200;	o_offsets[0x05] = 0x05<<1;	o_offsets[0x15] = 0x15<<1;
		i_offsets[0x06] = pc0 ;		i_offsets[0x16] = pc0 +p200;	o_offsets[0x06] = 0x06<<1;	o_offsets[0x16] = 0x16<<1;
		i_offsets[0x07] = pe0 ;		i_offsets[0x17] = pe0 +p200;	o_offsets[0x07] = 0x07<<1;	o_offsets[0x17] = 0x17<<1;
		i_offsets[0x08] = p100;		i_offsets[0x18] = p100+p200;	o_offsets[0x08] = 0x08<<1;	o_offsets[0x18] = 0x18<<1;
		i_offsets[0x09] = p120;		i_offsets[0x19] = p120+p200;	o_offsets[0x09] = 0x09<<1;	o_offsets[0x19] = 0x19<<1;
		i_offsets[0x0a] = p140;		i_offsets[0x1a] = p140+p200;	o_offsets[0x0a] = 0x0a<<1;	o_offsets[0x1a] = 0x1a<<1;
		i_offsets[0x0b] = p160;		i_offsets[0x1b] = p160+p200;	o_offsets[0x0b] = 0x0b<<1;	o_offsets[0x1b] = 0x1b<<1;
		i_offsets[0x0c] = p180;		i_offsets[0x1c] = p180+p200;	o_offsets[0x0c] = 0x0c<<1;	o_offsets[0x1c] = 0x1c<<1;
		i_offsets[0x0d] = p1a0;		i_offsets[0x1d] = p1a0+p200;	o_offsets[0x0d] = 0x0d<<1;	o_offsets[0x1d] = 0x1d<<1;
		i_offsets[0x0e] = p1c0;		i_offsets[0x1e] = p1c0+p200;	o_offsets[0x0e] = 0x0e<<1;	o_offsets[0x1e] = 0x1e<<1;
		i_offsets[0x0f] = p1e0;		i_offsets[0x1f] = p1e0+p200;	o_offsets[0x0f] = 0x0f<<1;	o_offsets[0x1f] = 0x1f<<1;

		// 2nd-pass (with-twiddles) I/O index offsets reverse things:
		// set 1 is w.r.to: t-array:								// set 2 w.r.to a-array:
		p_offsets[0x00] = 0x000;	p_offsets[0x01] = 0x200;		q_offsets[0x00] =  0;		q_offsets[0x10] = p10;
		p_offsets[0x02] = 0x100;	p_offsets[0x03] = 0x300;		q_offsets[0x01] = p1;		q_offsets[0x11] = p11;
		p_offsets[0x04] = 0x080;	p_offsets[0x05] = 0x280;		q_offsets[0x02] = p2;		q_offsets[0x12] = p12;
		p_offsets[0x06] = 0x180;	p_offsets[0x07] = 0x380;		q_offsets[0x03] = p3;		q_offsets[0x13] = p13;
		p_offsets[0x08] = 0x040;	p_offsets[0x09] = 0x240;		q_offsets[0x04] = p4;		q_offsets[0x14] = p14;
		p_offsets[0x0a] = 0x140;	p_offsets[0x0b] = 0x340;		q_offsets[0x05] = p5;		q_offsets[0x15] = p15;
		p_offsets[0x0c] = 0x0c0;	p_offsets[0x0d] = 0x2c0;		q_offsets[0x06] = p6;		q_offsets[0x16] = p16;
		p_offsets[0x0e] = 0x1c0;	p_offsets[0x0f] = 0x3c0;		q_offsets[0x07] = p7;		q_offsets[0x17] = p17;
		p_offsets[0x10] = 0x020;	p_offsets[0x11] = 0x220;		q_offsets[0x08] = p8;		q_offsets[0x18] = p18;
		p_offsets[0x12] = 0x120;	p_offsets[0x13] = 0x320;		q_offsets[0x09] = p9;		q_offsets[0x19] = p19;
		p_offsets[0x14] = 0x0a0;	p_offsets[0x15] = 0x2a0;		q_offsets[0x0a] = pa;		q_offsets[0x1a] = p1a;
		p_offsets[0x16] = 0x1a0;	p_offsets[0x17] = 0x3a0;		q_offsets[0x0b] = pb;		q_offsets[0x1b] = p1b;
		p_offsets[0x18] = 0x060;	p_offsets[0x19] = 0x260;		q_offsets[0x0c] = pc;		q_offsets[0x1c] = p1c;
		p_offsets[0x1a] = 0x160;	p_offsets[0x1b] = 0x360;		q_offsets[0x0d] = pd;		q_offsets[0x1d] = p1d;
		p_offsets[0x1c] = 0x0e0;	p_offsets[0x1d] = 0x2e0;		q_offsets[0x0e] = pe;		q_offsets[0x1e] = p1e;
		p_offsets[0x1e] = 0x1e0;	p_offsets[0x1f] = 0x3e0;		q_offsets[0x0f] = pf;		q_offsets[0x1f] = p1f;
		// Need real-array offsets due to pointer-cast-to-double:
		for(j=0;j<32;j++){
			p_offsets[j] <<= 1;
		//	printf("offsets[%2d] = %8d,%8d,%8d,%8d\n",j,i_offsets[j],o_offsets[j],p_offsets[j],q_offsets[j]);
		}
	//	exit(0);

	#else

		// Set array offsets for radix-64 DFT in/outputs:
		// set 1 is w.r.to: a-array:								// set 2 w.r.to t-array:
		i_offsets[0x00] = 0   ;		i_offsets[0x20] = 0   +p200;	o_offsets[0x00] = 0x00<<1;	o_offsets[0x20] = 0x20<<1;
		i_offsets[0x01] = p10 ;		i_offsets[0x21] = p10 +p200;	o_offsets[0x01] = 0x01<<1;	o_offsets[0x21] = 0x21<<1;
		i_offsets[0x02] = p20 ;		i_offsets[0x22] = p20 +p200;	o_offsets[0x02] = 0x02<<1;	o_offsets[0x22] = 0x22<<1;
		i_offsets[0x03] = p30 ;		i_offsets[0x23] = p30 +p200;	o_offsets[0x03] = 0x03<<1;	o_offsets[0x23] = 0x23<<1;
		i_offsets[0x04] = p40 ;		i_offsets[0x24] = p40 +p200;	o_offsets[0x04] = 0x04<<1;	o_offsets[0x24] = 0x24<<1;
		i_offsets[0x05] = p50 ;		i_offsets[0x25] = p50 +p200;	o_offsets[0x05] = 0x05<<1;	o_offsets[0x25] = 0x25<<1;
		i_offsets[0x06] = p60 ;		i_offsets[0x26] = p60 +p200;	o_offsets[0x06] = 0x06<<1;	o_offsets[0x26] = 0x26<<1;
		i_offsets[0x07] = p70 ;		i_offsets[0x27] = p70 +p200;	o_offsets[0x07] = 0x07<<1;	o_offsets[0x27] = 0x27<<1;
		i_offsets[0x08] = p80 ;		i_offsets[0x28] = p80 +p200;	o_offsets[0x08] = 0x08<<1;	o_offsets[0x28] = 0x28<<1;
		i_offsets[0x09] = p90 ;		i_offsets[0x29] = p90 +p200;	o_offsets[0x09] = 0x09<<1;	o_offsets[0x29] = 0x29<<1;
		i_offsets[0x0a] = pa0 ;		i_offsets[0x2a] = pa0 +p200;	o_offsets[0x0a] = 0x0a<<1;	o_offsets[0x2a] = 0x2a<<1;
		i_offsets[0x0b] = pb0 ;		i_offsets[0x2b] = pb0 +p200;	o_offsets[0x0b] = 0x0b<<1;	o_offsets[0x2b] = 0x2b<<1;
		i_offsets[0x0c] = pc0 ;		i_offsets[0x2c] = pc0 +p200;	o_offsets[0x0c] = 0x0c<<1;	o_offsets[0x2c] = 0x2c<<1;
		i_offsets[0x0d] = pd0 ;		i_offsets[0x2d] = pd0 +p200;	o_offsets[0x0d] = 0x0d<<1;	o_offsets[0x2d] = 0x2d<<1;
		i_offsets[0x0e] = pe0 ;		i_offsets[0x2e] = pe0 +p200;	o_offsets[0x0e] = 0x0e<<1;	o_offsets[0x2e] = 0x2e<<1;
		i_offsets[0x0f] = pf0 ;		i_offsets[0x2f] = pf0 +p200;	o_offsets[0x0f] = 0x0f<<1;	o_offsets[0x2f] = 0x2f<<1;
		i_offsets[0x10] = p100;		i_offsets[0x30] = p100+p200;	o_offsets[0x10] = 0x10<<1;	o_offsets[0x30] = 0x30<<1;
		i_offsets[0x11] = p110;		i_offsets[0x31] = p110+p200;	o_offsets[0x11] = 0x11<<1;	o_offsets[0x31] = 0x31<<1;
		i_offsets[0x12] = p120;		i_offsets[0x32] = p120+p200;	o_offsets[0x12] = 0x12<<1;	o_offsets[0x32] = 0x32<<1;
		i_offsets[0x13] = p130;		i_offsets[0x33] = p130+p200;	o_offsets[0x13] = 0x13<<1;	o_offsets[0x33] = 0x33<<1;
		i_offsets[0x14] = p140;		i_offsets[0x34] = p140+p200;	o_offsets[0x14] = 0x14<<1;	o_offsets[0x34] = 0x34<<1;
		i_offsets[0x15] = p150;		i_offsets[0x35] = p150+p200;	o_offsets[0x15] = 0x15<<1;	o_offsets[0x35] = 0x35<<1;
		i_offsets[0x16] = p160;		i_offsets[0x36] = p160+p200;	o_offsets[0x16] = 0x16<<1;	o_offsets[0x36] = 0x36<<1;
		i_offsets[0x17] = p170;		i_offsets[0x37] = p170+p200;	o_offsets[0x17] = 0x17<<1;	o_offsets[0x37] = 0x37<<1;
		i_offsets[0x18] = p180;		i_offsets[0x38] = p180+p200;	o_offsets[0x18] = 0x18<<1;	o_offsets[0x38] = 0x38<<1;
		i_offsets[0x19] = p190;		i_offsets[0x39] = p190+p200;	o_offsets[0x19] = 0x19<<1;	o_offsets[0x39] = 0x39<<1;
		i_offsets[0x1a] = p1a0;		i_offsets[0x3a] = p1a0+p200;	o_offsets[0x1a] = 0x1a<<1;	o_offsets[0x3a] = 0x3a<<1;
		i_offsets[0x1b] = p1b0;		i_offsets[0x3b] = p1b0+p200;	o_offsets[0x1b] = 0x1b<<1;	o_offsets[0x3b] = 0x3b<<1;
		i_offsets[0x1c] = p1c0;		i_offsets[0x3c] = p1c0+p200;	o_offsets[0x1c] = 0x1c<<1;	o_offsets[0x3c] = 0x3c<<1;
		i_offsets[0x1d] = p1d0;		i_offsets[0x3d] = p1d0+p200;	o_offsets[0x1d] = 0x1d<<1;	o_offsets[0x3d] = 0x3d<<1;
		i_offsets[0x1e] = p1e0;		i_offsets[0x3e] = p1e0+p200;	o_offsets[0x1e] = 0x1e<<1;	o_offsets[0x3e] = 0x3e<<1;
		i_offsets[0x1f] = p1f0;		i_offsets[0x3f] = p1f0+p200;	o_offsets[0x1f] = 0x1f<<1;	o_offsets[0x3f] = 0x3f<<1;

		po_br[0x0] = 0; po_br[0x1] = p8; po_br[0x2] = p4; po_br[0x3] = pc; po_br[0x4] = p2; po_br[0x5] = pa; po_br[0x6] = p6; po_br[0x7] = pe; po_br[0x8] = p1; po_br[0x9] = p9; po_br[0xa] = p5; po_br[0xb] = pd; po_br[0xc] = p3; po_br[0xd] = pb; po_br[0xe] = p7; po_br[0xf] = pf;
	#endif
	}

/*...The radix-1024 pass is here.	*/

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
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;

	#if RAD1 == 32	// 32 x 32 2-pass DFT version:

	// Gather the needed data and do 16 twiddleless length-32 subtransforms, with p-offsets in br16 order: 084c2a6e195d3b7f:
	// NOTE that RADIX_32_DIF outputs are IN-ORDER rather than BR:
		/*...Block 00: */	jt = j1      ;	jp = 0;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 10: */	jt = jt + p10;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 01: */	jt = j1 +  p8;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 11: */	jt = jt + p10;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 02: */	jt = j1 +  p4;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 12: */	jt = jt + p10;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 03: */	jt = j1 +  pc;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 13: */	jt = jt + p10;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 04: */	jt = j1 +  p2;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 14: */	jt = jt + p10;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 05: */	jt = j1 +  pa;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 15: */	jt = jt + p10;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 06: */	jt = j1 +  p6;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 16: */	jt = jt + p10;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 07: */	jt = j1 +  pe;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 17: */	jt = jt + p10;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 08: */	jt = j1 +  p1;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 18: */	jt = jt + p10;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 09: */	jt = j1 +  p9;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 19: */	jt = jt + p10;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 0a: */	jt = j1 +  p5;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 1a: */	jt = jt + p10;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 0b: */	jt = j1 +  pd;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 1b: */	jt = jt + p10;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 0c: */	jt = j1 +  p3;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 1c: */	jt = jt + p10;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 0d: */	jt = j1 +  pb;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 1d: */	jt = jt + p10;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 0e: */	jt = j1 +  p7;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 1e: */	jt = jt + p10;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 0f: */	jt = j1 +  pf;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 1f: */	jt = jt + p10;	jp+=32;		RADIX_32_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);

	/*...and now do 32 radix-32 subtransforms, including the internal twiddles (define F := e^(I*2pi/1024) and use F^n = -F^(n-0x200) for "easy mod" of larger powers):

			  0     1      2      3      4      5      6      7      8      9      a      b      c      d      e      f     10     11     12     13     14     15     16     17     18     19     1a     1b     1c     1d     1e     1f

Block  0:  {   1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1}
Block  1:  {   1, F^  1, F^  2, F^  3, F^  4, F^  5, F^  6, F^  7, F^  8, F^  9, F^  a, F^  b, F^  c, F^  d, F^  e, F^  f, F^ 10, F^ 11, F^ 12, F^ 13, F^ 14, F^ 15, F^ 16, F^ 17, F^ 18, F^ 19, F^ 1a, F^ 1b, F^ 1c, F^ 1d, F^ 1e, F^ 1f}
Block  2:  {   1, F^  2, F^  4, F^  6, F^  8, F^  a, F^  c, F^  e, F^ 10, F^ 12, F^ 14, F^ 16, F^ 18, F^ 1a, F^ 1c, F^ 1e, F^ 20, F^ 22, F^ 24, F^ 26, F^ 28, F^ 2a, F^ 2c, F^ 2e, F^ 30, F^ 32, F^ 34, F^ 36, F^ 38, F^ 3a, F^ 3c, F^ 3e}
Block  3:  {   1, F^  3, F^  6, F^  9, F^  c, F^  f, F^ 12, F^ 15, F^ 18, F^ 1b, F^ 1e, F^ 21, F^ 24, F^ 27, F^ 2a, F^ 2d, F^ 30, F^ 33, F^ 36, F^ 39, F^ 3c, F^ 3f, F^ 42, F^ 45, F^ 48, F^ 4b, F^ 4e, F^ 51, F^ 54, F^ 57, F^ 5a, F^ 5d}
Block  4:  {   1, F^  4, F^  8, F^  c, F^ 10, F^ 14, F^ 18, F^ 1c, F^ 20, F^ 24, F^ 28, F^ 2c, F^ 30, F^ 34, F^ 38, F^ 3c, F^ 40, F^ 44, F^ 48, F^ 4c, F^ 50, F^ 54, F^ 58, F^ 5c, F^ 60, F^ 64, F^ 68, F^ 6c, F^ 70, F^ 74, F^ 78, F^ 7c}
Block  5:  {   1, F^  5, F^  a, F^  f, F^ 14, F^ 19, F^ 1e, F^ 23, F^ 28, F^ 2d, F^ 32, F^ 37, F^ 3c, F^ 41, F^ 46, F^ 4b, F^ 50, F^ 55, F^ 5a, F^ 5f, F^ 64, F^ 69, F^ 6e, F^ 73, F^ 78, F^ 7d, F^ 82, F^ 87, F^ 8c, F^ 91, F^ 96, F^ 9b}
Block  6:  {   1, F^  6, F^  c, F^ 12, F^ 18, F^ 1e, F^ 24, F^ 2a, F^ 30, F^ 36, F^ 3c, F^ 42, F^ 48, F^ 4e, F^ 54, F^ 5a, F^ 60, F^ 66, F^ 6c, F^ 72, F^ 78, F^ 7e, F^ 84, F^ 8a, F^ 90, F^ 96, F^ 9c, F^ a2, F^ a8, F^ ae, F^ b4, F^ ba}
Block  7:  {   1, F^  7, F^  e, F^ 15, F^ 1c, F^ 23, F^ 2a, F^ 31, F^ 38, F^ 3f, F^ 46, F^ 4d, F^ 54, F^ 5b, F^ 62, F^ 69, F^ 70, F^ 77, F^ 7e, F^ 85, F^ 8c, F^ 93, F^ 9a, F^ a1, F^ a8, F^ af, F^ b6, F^ bd, F^ c4, F^ cb, F^ d2, F^ d9}
Block  8:  {   1, F^  8, F^ 10, F^ 18, F^ 20, F^ 28, F^ 30, F^ 38, F^ 40, F^ 48, F^ 50, F^ 58, F^ 60, F^ 68, F^ 70, F^ 78, F^ 80, F^ 88, F^ 90, F^ 98, F^ a0, F^ a8, F^ b0, F^ b8, F^ c0, F^ c8, F^ d0, F^ d8, F^ e0, F^ e8, F^ f0, F^ f8}
Block  9:  {   1, F^  9, F^ 12, F^ 1b, F^ 24, F^ 2d, F^ 36, F^ 3f, F^ 48, F^ 51, F^ 5a, F^ 63, F^ 6c, F^ 75, F^ 7e, F^ 87, F^ 90, F^ 99, F^ a2, F^ ab, F^ b4, F^ bd, F^ c6, F^ cf, F^ d8, F^ e1, F^ ea, F^ f3, F^ fc, F^105, F^10e, F^117}
Block  a:  {   1, F^  a, F^ 14, F^ 1e, F^ 28, F^ 32, F^ 3c, F^ 46, F^ 50, F^ 5a, F^ 64, F^ 6e, F^ 78, F^ 82, F^ 8c, F^ 96, F^ a0, F^ aa, F^ b4, F^ be, F^ c8, F^ d2, F^ dc, F^ e6, F^ f0, F^ fa, F^104, F^10e, F^118, F^122, F^12c, F^136}
Block  b:  {   1, F^  b, F^ 16, F^ 21, F^ 2c, F^ 37, F^ 42, F^ 4d, F^ 58, F^ 63, F^ 6e, F^ 79, F^ 84, F^ 8f, F^ 9a, F^ a5, F^ b0, F^ bb, F^ c6, F^ d1, F^ dc, F^ e7, F^ f2, F^ fd, F^108, F^113, F^11e, F^129, F^134, F^13f, F^14a, F^155}
Block  c:  {   1, F^  c, F^ 18, F^ 24, F^ 30, F^ 3c, F^ 48, F^ 54, F^ 60, F^ 6c, F^ 78, F^ 84, F^ 90, F^ 9c, F^ a8, F^ b4, F^ c0, F^ cc, F^ d8, F^ e4, F^ f0, F^ fc, F^108, F^114, F^120, F^12c, F^138, F^144, F^150, F^15c, F^168, F^174}
Block  d:  {   1, F^  d, F^ 1a, F^ 27, F^ 34, F^ 41, F^ 4e, F^ 5b, F^ 68, F^ 75, F^ 82, F^ 8f, F^ 9c, F^ a9, F^ b6, F^ c3, F^ d0, F^ dd, F^ ea, F^ f7, F^104, F^111, F^11e, F^12b, F^138, F^145, F^152, F^15f, F^16c, F^179, F^186, F^193}
Block  e:  {   1, F^  e, F^ 1c, F^ 2a, F^ 38, F^ 46, F^ 54, F^ 62, F^ 70, F^ 7e, F^ 8c, F^ 9a, F^ a8, F^ b6, F^ c4, F^ d2, F^ e0, F^ ee, F^ fc, F^10a, F^118, F^126, F^134, F^142, F^150, F^15e, F^16c, F^17a, F^188, F^196, F^1a4, F^1b2}
Block  f:  {   1, F^  f, F^ 1e, F^ 2d, F^ 3c, F^ 4b, F^ 5a, F^ 69, F^ 78, F^ 87, F^ 96, F^ a5, F^ b4, F^ c3, F^ d2, F^ e1, F^ f0, F^ ff, F^10e, F^11d, F^12c, F^13b, F^14a, F^159, F^168, F^177, F^186, F^195, F^1a4, F^1b3, F^1c2, F^1d1}
Block 10:  {   1, F^ 10, F^ 20, F^ 30, F^ 40, F^ 50, F^ 60, F^ 70, F^ 80, F^ 90, F^ a0, F^ b0, F^ c0, F^ d0, F^ e0, F^ f0, F^100, F^110, F^120, F^130, F^140, F^150, F^160, F^170, F^180, F^190, F^1a0, F^1b0, F^1c0, F^1d0, F^1e0, F^1f0}
Block 11:  {   1, F^ 11, F^ 22, F^ 33, F^ 44, F^ 55, F^ 66, F^ 77, F^ 88, F^ 99, F^ aa, F^ bb, F^ cc, F^ dd, F^ ee, F^ ff, F^110, F^121, F^132, F^143, F^154, F^165, F^176, F^187, F^198, F^1a9, F^1ba, F^1cb, F^1dc, F^1ed, F^1fe,-F^ 0f}
Block 12:  {   1, F^ 12, F^ 24, F^ 36, F^ 48, F^ 5a, F^ 6c, F^ 7e, F^ 90, F^ a2, F^ b4, F^ c6, F^ d8, F^ ea, F^ fc, F^10e, F^120, F^132, F^144, F^156, F^168, F^17a, F^18c, F^19e, F^1b0, F^1c2, F^1d4, F^1e6, F^1f8,-F^ 0a,-F^ 1c,-F^ 2e}
Block 13:  {   1, F^ 13, F^ 26, F^ 39, F^ 4c, F^ 5f, F^ 72, F^ 85, F^ 98, F^ ab, F^ be, F^ d1, F^ e4, F^ f7, F^10a, F^11d, F^130, F^143, F^156, F^169, F^17c, F^18f, F^1a2, F^1b5, F^1c8, F^1db, F^1ee,-F^ 01,-F^ 14,-F^ 27,-F^ 3a,-F^ 4d}
Block 14:  {   1, F^ 14, F^ 28, F^ 3c, F^ 50, F^ 64, F^ 78, F^ 8c, F^ a0, F^ b4, F^ c8, F^ dc, F^ f0, F^104, F^118, F^12c, F^140, F^154, F^168, F^17c, F^190, F^1a4, F^1b8, F^1cc, F^1e0, F^1f4,-F^ 08,-F^ 1c,-F^ 30,-F^ 44,-F^ 58,-F^ 6c}
Block 15:  {   1, F^ 15, F^ 2a, F^ 3f, F^ 54, F^ 69, F^ 7e, F^ 93, F^ a8, F^ bd, F^ d2, F^ e7, F^ fc, F^111, F^126, F^13b, F^150, F^165, F^17a, F^18f, F^1a4, F^1b9, F^1ce, F^1e3, F^1f8,-F^ 0d,-F^ 22,-F^ 37,-F^ 4c,-F^ 61,-F^ 76,-F^ 8b}
Block 16:  {   1, F^ 16, F^ 2c, F^ 42, F^ 58, F^ 6e, F^ 84, F^ 9a, F^ b0, F^ c6, F^ dc, F^ f2, F^108, F^11e, F^134, F^14a, F^160, F^176, F^18c, F^1a2, F^1b8, F^1ce, F^1e4, F^1fa,-F^ 10,-F^ 26,-F^ 3c,-F^ 52,-F^ 68,-F^ 7e,-F^ 94,-F^ aa}
Block 17:  {   1, F^ 17, F^ 2e, F^ 45, F^ 5c, F^ 73, F^ 8a, F^ a1, F^ b8, F^ cf, F^ e6, F^ fd, F^114, F^12b, F^142, F^159, F^170, F^187, F^19e, F^1b5, F^1cc, F^1e3, F^1fa,-F^ 11,-F^ 28,-F^ 3f,-F^ 56,-F^ 6d,-F^ 84,-F^ 9b,-F^ b2,-F^ c9}
Block 18:  {   1, F^ 18, F^ 30, F^ 48, F^ 60, F^ 78, F^ 90, F^ a8, F^ c0, F^ d8, F^ f0, F^108, F^120, F^138, F^150, F^168, F^180, F^198, F^1b0, F^1c8, F^1e0, F^1f8,-F^ 10,-F^ 28,-F^ 40,-F^ 58,-F^ 70,-F^ 88,-F^ a0,-F^ b8,-F^ d0,-F^ e8}
Block 19:  {   1, F^ 19, F^ 32, F^ 4b, F^ 64, F^ 7d, F^ 96, F^ af, F^ c8, F^ e1, F^ fa, F^113, F^12c, F^145, F^15e, F^177, F^190, F^1a9, F^1c2, F^1db, F^1f4,-F^ 0d,-F^ 26,-F^ 3f,-F^ 58,-F^ 71,-F^ 8a,-F^ a3,-F^ bc,-F^ d5,-F^ ee,-F^107}
Block 1a:  {   1, F^ 1a, F^ 34, F^ 4e, F^ 68, F^ 82, F^ 9c, F^ b6, F^ d0, F^ ea, F^104, F^11e, F^138, F^152, F^16c, F^186, F^1a0, F^1ba, F^1d4, F^1ee,-F^ 08,-F^ 22,-F^ 3c,-F^ 56,-F^ 70,-F^ 8a,-F^ a4,-F^ be,-F^ d8,-F^ f2,-F^10c,-F^126}
Block 1b:  {   1, F^ 1b, F^ 36, F^ 51, F^ 6c, F^ 87, F^ a2, F^ bd, F^ d8, F^ f3, F^10e, F^129, F^144, F^15f, F^17a, F^195, F^1b0, F^1cb, F^1e6,-F^ 01,-F^ 1c,-F^ 37,-F^ 52,-F^ 6d,-F^ 88,-F^ a3,-F^ be,-F^ d9,-F^ f4,-F^10f,-F^12a,-F^145}
Block 1c:  {   1, F^ 1c, F^ 38, F^ 54, F^ 70, F^ 8c, F^ a8, F^ c4, F^ e0, F^ fc, F^118, F^134, F^150, F^16c, F^188, F^1a4, F^1c0, F^1dc, F^1f8,-F^ 14,-F^ 30,-F^ 4c,-F^ 68,-F^ 84,-F^ a0,-F^ bc,-F^ d8,-F^ f4,-F^110,-F^12c,-F^148,-F^164}
Block 1d:  {   1, F^ 1d, F^ 3a, F^ 57, F^ 74, F^ 91, F^ ae, F^ cb, F^ e8, F^105, F^122, F^13f, F^15c, F^179, F^196, F^1b3, F^1d0, F^1ed,-F^ 0a,-F^ 27,-F^ 44,-F^ 61,-F^ 7e,-F^ 9b,-F^ b8,-F^ d5,-F^ f2,-F^10f,-F^12c,-F^149,-F^166,-F^183}
Block 1e:  {   1, F^ 1e, F^ 3c, F^ 5a, F^ 78, F^ 96, F^ b4, F^ d2, F^ f0, F^10e, F^12c, F^14a, F^168, F^186, F^1a4, F^1c2, F^1e0, F^1fe,-F^ 1c,-F^ 3a,-F^ 58,-F^ 76,-F^ 94,-F^ b2,-F^ d0,-F^ ee,-F^10c,-F^12a,-F^148,-F^166,-F^184,-F^1a2}
Block 1f:  {   1, F^ 1f, F^ 3e, F^ 5d, F^ 7c, F^ 9b, F^ ba, F^ d9, F^ f8, F^117, F^136, F^155, F^174, F^193, F^1b2, F^1d1, F^1f0,-F^ 0f,-F^ 2e,-F^ 4d,-F^ 6c,-F^ 8b,-F^ aa,-F^ c9,-F^ e8,-F^107,-F^126,-F^145,-F^164,-F^183,-F^1a2,-F^1c1} .

	Now further reduce remaining powers:
		~ denotes complex conjugation, i.e. ~F^n := ( Re(E),-Im(F)),
		* denotes interchange of real and imaginary part, i.e. *F^n := ( Im(F), Re(F)),
		- denotes negation, i.e. -F^n := (-Re(F),-Im(F)),
	and any combination of these operators is evaluated right-to-left, e.g.
		-~*F^n = -~(*F^n) = -~( Im(F), Re(F)) = -( Im(F),-Re(F)) = (-Im(F),+Re(F)) .
	Note that the - and ~ operators commute, as do the - and *, but ~ and * anticommute, i.e. ~*F^n = -*~F^n) ,
	and the values of the individual exponentials are, in terms of the sincos parameters defined in this module:

		F^ 80 = exp(i* 1*twopi/8) = isrt2*( 1    , 1    )
		F^100 = exp(i* 2*twopi/8) =       ( 0    , 1    ) = I
		F^180 = exp(i* 3*twopi/8) = isrt2*(-1    , 1    ) = *~e^80

		F^{ 81, 82, ... , fe, ff} =  *F^{7f,7e, ... , 2, 1}
		F^{101,102, ... ,17e,17f} = I.F^{ 1, 2, ... ,7e,7f} = *~F^{ 1, 2, ... ,7e,7f}	(I.{} denotes complex multiply by I)
		F^{181,182, ... ,1fe,1ff} = -~F^{7f,7e, ... , 2, 1} .

	Exploiting these symmetries allows our 32 x 32 twiddles matrix to be expressed in terms of the powers F^1-7f, 80,
	the imaginary constant I and isrt2 (via F^80 = isrt2*(1+I)) as follows - Note the left halves of the even rows
	(counting the top row 	as "row 0") are identical to the twiddles-sets of the radix-512 DFT.
	(We omit the F^ here to save space, and just list the applicable complex-arithmetic operation and the power of F):

	 0  1  2  3  4   5   6   7   8    9    a    b    c    d    e    f   10   11   12   13   14   15   16   17   18   19   1a   1b   1c   1d   1e   1f

00	{1, 1, 1, 1, 1,  1,  1,  1,  1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1}
01	{1, 1, 2, 3, 4,  5,  6,  7,  8,   9,   a,   b,   c,   d,   e,   f,  10,  11,  12,  13,  14,  15,  16,  17,  18,  19,  1a,  1b,  1c,  1d,  1e,  1f}
02	{1, 2, 4, 6, 8,  a,  c,  e, 10,  12,  14,  16,  18,  1a,  1c,  1e,  20,  22,  24,  26,  28,  2a,  2c,  2e,  30,  32,  34,  36,  38,  3a,  3c,  3e}
03	{1, 3, 6, 9, c,  f, 12, 15, 18,  1b,  1e,  21,  24,  27,  2a,  2d,  30,  33,  36,  39,  3c,  3f,  42,  45,  48,  4b,  4e,  51,  54,  57,  5a,  5d}
04	{1, 4, 8, c,10, 14, 18, 1c, 20,  24,  28,  2c,  30,  34,  38,  3c,  40,  44,  48,  4c,  50,  54,  58,  5c,  60,  64,  68,  6c,  70,  74,  78,  7c}
05	{1, 5, a, f,14, 19, 1e, 23, 28,  2d,  32,  37,  3c,  41,  46,  4b,  50,  55,  5a,  5f,  64,  69,  6e,  73,  78,  7d, *7e, *79, *74, *6f, *6a, *65}
06	{1, 6, c,12,18, 1e, 24, 2a, 30,  36,  3c,  42,  48,  4e,  54,  5a,  60,  66,  6c,  72,  78,  7e, *7c, *76, *70, *6a, *64, *5e, *58, *52, *4c, *46}
07	{1, 7, e,15,1c, 23, 2a, 31, 38,  3f,  46,  4d,  54,  5b,  62,  69,  70,  77,  7e, *7b, *74, *6d, *66, *5f, *58, *51, *4a, *43, *3c, *35, *2e, *27}
08	{1, 8,10,18,20, 28, 30, 38, 40,  48,  50,  58,  60,  68,  70,  78,  80, *78, *70, *68, *60, *58, *50, *48, *40, *38, *30, *28, *20, *18, *10, * 8}
09	{1, 9,12,1b,24, 2d, 36, 3f, 48,  51,  5a,  63,  6c,  75,  7e, *79, *70, *67, *5e, *55, *4c, *43, *3a, *31, *28, *1f, *16, * d, * 4,*~ 5,*~ e,*~17}
0a	{1, a,14,1e,28, 32, 3c, 46, 50,  5a,  64,  6e,  78, *7e, *74, *6a, *60, *56, *4c, *42, *38, *2e, *24, *1a, *10, * 6,*~ 4,*~ e,*~18,*~22,*~2c,*~36}
0b	{1, b,16,21,2c, 37, 42, 4d, 58,  63,  6e,  79, *7c, *71, *66, *5b, *50, *45, *3a, *2f, *24, *19, * e, * 3,*~ 8,*~13,*~1e,*~29,*~34,*~3f,*~4a,*~55}
0c	{1, c,18,24,30, 3c, 48, 54, 60,  6c,  78, *7c, *70, *64, *58, *4c, *40, *34, *28, *1c, *10, * 4,*~ 8,*~14,*~20,*~2c,*~38,*~44,*~50,*~5c,*~68,*~74}
 d:	{1, d,1a,27,34, 41, 4e, 5b, 68,  75,  82,  8f,  9c,  a9,  b6,  c3,  d0,  dd,  ea,  f7, 104, 111, 11e, 12b, 138, 145, 152, 15f, 16c, 179, 186, 193}
0d	{1, d,1a,27,34, 41, 4e, 5b, 68,  75, *7e, *71, *64, *57, *4a, *3d, *30, *23, *16, * 9,*~ 4,*~11,*~1e,*~2b,*~38,*~45,*~52,*~5f,*~6c,*~79,-~7a,-~6d}
0e	{1, e,1c,2a,38, 46, 54, 62, 70,  7e, *74, *66, *58, *4a, *3c, *2e, *20, *12, * 4,*~ a,*~18,*~26,*~34,*~42,*~50,*~5e,*~6c,*~7a,-~78,-~6a,-~5c,-~4e}
0f	{1, f,1e,2d,3c, 4b, 5a, 69, 78, *79, *6a, *5b, *4c, *3d, *2e, *1f, *10, * 1,*~ e,*~1d,*~2c,*~3b,*~4a,*~59,*~68,*~77,-~7a,-~6b,-~5c,-~4d,-~3e,-~2f}
10	{1,10,20,30,40, 50, 60, 70, 80, *70, *60, *50, *40, *30, *20, *10,I.{},*~10,*~20,*~30,*~40,*~50,*~60,*~70,*~80,-~70,-~60,-~50,-~40,-~30,-~20,-~10}
11	{1,11,22,33,44, 55, 66, 77,*78, *67, *56, *45, *34, *23, *12, * 1,*~10,*~21,*~32,*~43,*~54,*~65,*~76,-~79,-~68,-~57,-~46,-~35,-~24,-~13,-~ 2, - f}
12	{1,12,24,36,48, 5a, 6c, 7e,*70, *5e, *4c, *3a, *28, *16, * 4,*~ e,*~20,*~32,*~44,*~56,*~68,*~7a,-~74,-~62,-~50,-~3e,-~2c,-~1a,-~ 8, - a, -1c, -2e}
13	{1,13,26,39,4c, 5f, 72,*7b,*68, *55, *42, *2f, *1c, * 9,*~ a,*~1d,*~30,*~43,*~56,*~69,*~7c,-~71,-~5e,-~4b,-~38,-~25,-~12, - 1, -14, -27, -3a, -4d}
14	{1,14,28,3c,50, 64, 78,*74,*60, *4c, *38, *24, *10,*~ 4,*~18,*~2c,*~40,*~54,*~68,*~7c,-~70,-~5c,-~48,-~34,-~20,-~ c, - 8, -1c, -30, -44, -58, -6c}
15	{1,15,2a,3f,54, 69, 7e,*6d,*58, *43, *2e, *19, * 4,*~11,*~26,*~3b,*~50,*~65,*~7a,-~71,-~5c,-~47,-~32,-~1d,-~ 8, - d, -22, -37, -4c, -61, -76,-*75}
16	{1,16,2c,42,58, 6e,*7c,*66,*50, *3a, *24, * e,*~ 8,*~1e,*~34,*~4a,*~60,*~76,-~74,-~5e,-~48,-~32,-~1c,-~ 6, -10, -26, -3c, -52, -68, -7e,-*6c,-*56}
17	{1,17,2e,45,5c, 73,*76,*5f,*48, *31, *1a, * 3,*~14,*~2b,*~42,*~59,*~70,-~79,-~62,-~4b,-~34,-~1d,-~ 6, -11, -28, -3f, -56, -6d,-*7c,-*65,-*4e,-*37}
18	{1,18,30,48,60, 78,*70,*58,*40, *28, *10,*~ 8,*~20,*~38,*~50,*~68,*~80,-~68,-~50,-~38,-~20,-~ 8, -10, -28, -40, -58, -70,-*78,-*60,-*48,-*30,-*18}
19	{1,19,32,4b,64, 7d,*6a,*51,*38, *1f, * 6,*~13,*~2c,*~45,*~5e,*~77,-~70,-~57,-~3e,-~25,-~ c, - d, -26, -3f, -58, -71,-*76,-*5d,-*44,-*2b,-*12,~* 7}
1a	{1,1a,34,4e,68,*7e,*64,*4a,*30, *16,*~ 4,*~1e,*~38,*~52,*~6c,-~7a,-~60,-~46,-~2c,-~12, - 8, -22, -3c, -56, -70,-*76,-*5c,-*42,-*28,-* e,~* c,~*26}
1b	{1,1b,36,51,6c,*79,*5e,*43,*28, * d,*~ e,*~29,*~44,*~5f,*~7a,-~6b,-~50,-~35,-~1a, - 1, -1c, -37, -52, -6d,-*78,-*5d,-*42,-*27,-* c,~* f,~*2a,~*45}
1c	{1,1c,38,54,70,*74,*58,*3c,*20, * 4,*~18,*~34,*~50,*~6c,-~78,-~5c,-~40,-~24,-~ 8, -14, -30, -4c, -68,-*7c,-*60,-*44,-*28,-* c,~*10,~*2c,~*48,~*64}
1d	{1,1d,3a,57,74,*6f,*52,*35,*18,*~ 5,*~22,*~3f,*~5c,*~79,-~6a,-~4d,-~30,-~13, - a, -27, -44, -61, -7e,-*65,-*48,-*2b,-* e,~* f,~*2c,~*49,~*66, ~7d}
1e	{1,1e,3c,5a,78,*6a,*4c,*2e,*10,*~ e,*~2c,*~4a,*~68,-~7a,-~5c,-~3e,-~20,-~ 2, -1c, -3a, -58, -76,-*6c,-*4e,-*30,-*12,~* c,~*2a,~*48,~*66, ~7c, ~5e}
1f	{1,1f,3e,5d,7c,*65,*46,*27,* 8,*~17,*~36,*~55,*~74,-~6d,-~4e,-~2f,-~10, - f, -2e, -4d, -6c,-*75,-*56,-*37,-*18,~* 7,~*26,~*45,~*64, ~7d, ~5e, ~3f} .

	Or with rows (can do by pasting already-BRed-index col at left and sorting) and columns rearranged in BR32 {0+8+4+c+2+a+6+e+1+9+5+d+3+b+7+f+}
	(where + denotes 'preceding index + 0x10'; pairwise-swap cols 1/8,2/4,3/c,5/a,7/e,b/d, then paste corr. '+' col to right of each) order:

	 0   10   8   18  4   14    c   1c  2   12    a   1a   6   16    e   1e  1   11    9   19   5   15    d   1d  3   13    b   1b   7   17    f   1f

00	{1,   1,  1,   1, 1,   1,   1,   1, 1,   1,   1,   1,  1,   1,   1,   1, 1,   1,   1,   1,  1,   1,   1,   1, 1,   1,   1,   1,  1,   1,   1,   1}
10	{1,I.{}, 80,*~80,40,*~40, *40,-~40,20,*~20, *60,-~60, 60,*~60, *20,-~20,10,*~10, *70,-~70, 50,*~50, *30,-~30,30,*~30, *50,-~50, 70,*~70, *10,-~10}
08	{1,  80, 40, *40,20, *60,  60, *20,10, *70,  50, *30, 30, *50,  70, *10, 8, *78,  48, *38, 28, *58,  68, *18,18, *68,  58, *28, 38, *48,  78, * 8}
18	{1,*~80,*40, -40,60,-~20,*~20,-*60,30,-~50, *10, -70,*70, -10,*~50,-*30,18,-~68, *28, -58, 78,-~ 8,*~38,-*48,48,-~38,*~ 8,-*78,*58, -28,*~68,-*18}
04	{1,  40, 20,  60,10,  50,  30,  70, 8,  48,  28,  68, 18,  58,  38,  78, 4,  44,  24,  64, 14,  54,  34,  74, c,  4c,  2c,  6c, 1c,  5c,  3c,  7c}
14	{1,*~40,*60,-~20,50,-~70, *10, -30,28,*~68, *38, - 8, 78,-~48,*~18, -58,14,*~54, *4c,-~ c, 64,-~5c,*~ 4, -44,3c,*~7c, *24, -1c,*74,-~34,*~2c, -6c}
0c	{1, *40, 60,*~20,30, *10, *70,*~50,18, *28,  78,*~38, 48,*~ 8, *58,*~68, c, *34,  6c,*~2c, 3c, * 4, *64,*~5c,24, *1c, *7c,*~44, 54,*~14, *4c,*~74}
1c	{1,-~40,*20,-*60,70, -30,*~50,~*10,38,-~ 8,*~18,-*28,*58, -68,-~78,~*48,1c,-~24, * 4,-*44,*74, -4c,*~6c,~*2c,54, -14,*~34,-* c,*3c,-*7c,-~5c,~*64}
02	{1,  20, 10,  30, 8,  28,  18,  38, 4,  24,  14,  34,  c,  2c,  1c,  3c, 2,  22,  12,  32,  a,  2a,  1a,  3a, 6,  26,  16,  36,  e,  2e,  1e,  3e}
12	{1,*~20,*70,-~50,48,*~68, *28,-~ 8,24,*~44, *4c,-~2c, 6c,-~74, * 4, -1c,12,*~32, *5e,-~3e, 5a,*~7a, *16, - a,36,*~56, *3a,-~1a, 7e,-~62,*~ e, -2e}
0a	{1, *60, 50, *10,28, *38,  78,*~18,14, *4c,  64,*~ 4, 3c, *24, *74,*~2c, a, *56,  5a, * 6, 32, *2e, *7e,*~22,1e, *42,  6e,*~ e, 46, *1a, *6a,*~36}
1a	{1,-~60,*30, -70,68, - 8,*~38,-*28,34,-~2c,*~ 4,-*5c,*64, -3c,*~6c,~* c,1a,-~46, *16,-*76,*7e, -22,*~52,-* e,4e,-~12,*~1e,-*42,*4a, -56,-~7a,~*26}
06	{1,  60, 30, *70,18,  78,  48, *58, c,  6c,  3c, *64, 24, *7c,  54, *4c, 6,  66,  36, *6a, 1e,  7e,  4e, *52,12,  72,  42, *5e, 2a, *76,  5a, *46}
16	{1,*~60,*50, -10,58,-~48,*~ 8, -68,2c,-~74, *24, -3c,*7c,-~1c,*~34,-*6c,16,*~76, *3a, -26, 6e,-~32,*~1e, -7e,42,-~5e, * e, -52,*66,-~ 6,*~4a,-*56}
0e	{1, *20, 70,*~50,38,*~18, *58,-~78,1c, * 4, *74,*~6c, 54,*~34, *3c,-~5c, e, *12,  7e,*~5e, 46,*~26, *4a,-~6a,2a,*~ a, *66,*~7a, 62,*~42, *2e,-~4e}
1e	{1,-~20,*10,-*30,78, -58,*~68,~*48,3c, -1c,*~2c,~* c,*4c,-*6c,-~5c, ~7c,1e,-~ 2,*~ e,-*12,*6a, -76,-~7a,~*66,5a, -3a,*~4a,~*2a,*2e,-*4e,-~3e, ~5e}
01	{1,  10,  8,  18, 4,  14,   c,  1c, 2,  12,   a,  1a,  6,  16,   e,  1e, 1,  11,   9,  19,  5,  15,   d,  1d, 3,  13,   b,  1b,  7,  17,   f,  1f}
11	{1,*~10,*78,-~68,44,*~54, *34,-~24,22,*~32, *56,-~46, 66,*~76, *12,-~ 2,11,*~21, *67,-~57, 55,*~65, *23,-~13,33,*~43, *45,-~35, 77,-~79, * 1, - f}
09	{1, *70, 48, *28,24, *4c,  6c, * 4,12, *5e,  5a, *16, 36, *3a,  7e,*~ e, 9, *67,  51, *1f, 2d, *43,  75,*~ 5,1b, *55,  63, * d, 3f, *31, *79,*~17}
19	{1,-~70,*38, -58,64,-~ c,*~2c,-*44,32,-~3e, * 6,-*76,*6a, -26,*~5e,-*12,19,-~57, *1f, -71, 7d, - d,*~45,-*2b,4b,-~25,*~13,-*5d,*51, -3f,*~77,~* 7}
05	{1,  50, 28,  78,14,  64,  3c, *74, a,  5a,  32, *7e, 1e,  6e,  46, *6a, 5,  55,  2d,  7d, 19,  69,  41, *6f, f,  5f,  37, *79, 23,  73,  4b, *65}
15	{1,*~50,*58,-~ 8,54,-~5c, * 4, -4c,2a,*~7a, *2e, -22, 7e,-~32,*~26, -76,15,*~65, *43, - d, 69,-~47,*~11, -61,3f,-~71, *19, -37,*6d,-~1d,*~3b,-*75}
0d	{1, *30, 68,*~38,34,*~ 4, *64,*~6c,1a, *16, *7e,*~52, 4e,*~1e, *4a,-~7a, d, *23,  75,*~45, 41,*~11, *57,*~79,27, * 9, *71,*~5f, 5b,*~2b, *3d,-~6d}
1d	{1,-~30,*18,-*48,74, -44,*~5c,~*2c,3a, - a,*~22,-* e,*52, -7e,-~6a,~*66,1d,-~13,*~ 5,-*2b,*6f, -61,*~79,~*49,57, -27,*~3f,~* f,*35,-*65,-~4d, ~7d}
03	{1,  30, 18,  48, c,  3c,  24,  54, 6,  36,  1e,  4e, 12,  42,  2a,  5a, 3,  33,  1b,  4b,  f,  3f,  27,  57, 9,  39,  21,  51, 15,  45,  2d,  5d}
13	{1,*~30,*68,-~38,4c,*~7c, *1c, -14,26,*~56, *42,-~12, 72,-~5e,*~ a, -3a,13,*~43, *55,-~25, 5f,-~71, * 9, -27,39,*~69, *2f, - 1,*7b,-~4b,*~1d, -4d}
0b	{1, *50, 58,*~ 8,2c, *24, *7c,*~34,16, *3a,  6e,*~1e, 42, * e, *66,*~4a, b, *45,  63,*~13, 37, *19, *71,*~3f,21, *2f,  79,*~29, 4d, * 3, *5b,*~55}
1b	{1,-~50,*28,-*78,6c, -1c,*~44,-* c,36,-~1a,*~ e,-*42,*5e, -52,*~7a,~*2a,1b,-~35, * d,-*5d,*79, -37,*~5f,~* f,51, - 1,*~29,-*27,*43, -6d,-~6b,~*45}
07	{1,  70, 38, *58,1c, *74,  54, *3c, e,  7e,  46, *4a, 2a, *66,  62, *2e, 7,  77,  3f, *51, 23, *6d,  5b, *35,15, *7b,  4d, *43, 31, *5f,  69, *27}
17	{1,*~70,*48, -28,5c,-~34,*~14,-*7c,2e,-~62, *1a, -56,*76,-~ 6,*~42,-*4e,17,-~79, *31, -3f, 73,-~1d,*~2b,-*65,45,-~4b, * 3, -6d,*5f, -11,*~59,-*37}
0f	{1, *10, 78,*~68,3c,*~2c, *4c,-~5c,1e,*~ e, *6a,-~7a, 5a,*~4a, *2e,-~3e, f, * 1, *79,*~77, 4b,*~3b, *3d,-~4d,2d,*~1d, *5b,-~6b, 69,*~59, *1f,-~2f}
1f	{1,-~10,* 8,-*18,7c, -6c,*~74,~*64,3e, -2e,*~36,~*26,*46,-*56,-~4e, ~5e,1f, - f,*~17,~* 7,*65,-*75,-~6d, ~7d,5d, -4d,*~55,~*45,*27,-*37,-~2f, ~3f} .

Express powers which are multiples of 2,4,8,16,32 resp. in terms of E,D,C,B,A := radix-512,256,128,64,32 fundamental roots

	 0   10    8    18    4    14      c     1c    2     12      a     1a     6     16      e     1e    1     11      9     19     5     15      d     1d    3     13      b     1b     7     17      f     1f
00
10	{1, I.{}, A^4,*~A^4,A^2 ,*~A^2 , *A^2 ,-~A^2 ,A^1 ,*~A^1 , *A^3 ,-~A^3 , A^3 ,*~A^3 , *A^1 ,-~A^1 ,B^1 ,*~B^1 , *B^7 ,-~B^7 , B^5 ,*~B^5 , *B^3 ,-~B^3 ,B^3 ,*~B^3 , *B^5 ,-~B^5 , B^7 ,*~B^7 , *B^1 ,-~B^1 }
08	{1,  A^4, A^2, *A^2,A^1 , *A^3 ,  A^3 , *A^1 ,B^1 , *B^7 ,  B^5 , *B^3 , B^3 , *B^5 ,  B^7 , *B^1 ,C^1 , *C^f ,  C^9 , *C^7 , C^5 , *C^b ,  C^d , *C^3 ,C^3 , *C^d ,  C^b , *C^5 , C^7 , *C^9 ,  C^f , *C^1 }
18	{1,*~A^4,*A^2, -A^2,A^3 ,-~A^1 ,*~A^1 ,-*A^3 ,B^3 ,-~B^5 , *B^1 , -B^7 ,*B^7 , -B^1 ,*~B^5 ,-*B^3 ,C^3 ,-~C^d , *C^5 , -C^b , C^f ,-~C^1 ,*~C^7 ,-*C^9 ,C^9 ,-~C^7 ,*~C^1 ,-*C^f ,*C^b , -C^5 ,*~C^d ,-*C^3 }
04	{1,  A^2, A^1,  A^3,B^1 ,  B^5 ,  B^3 ,  B^7 ,C^1 ,  C^9 ,  C^5 ,  C^d , C^3 ,  C^b ,  C^7 ,  C^f ,D^01,  D^11,  D^09,  D^19, D^05,  D^15,  D^0d,  D^1d,D^03,  D^13,  D^0b,  D^1b, D^07,  D^17,  D^0f,  D^1f}
14	{1,*~A^2,*A^3,-~A^1,B^5 ,-~B^7 , *B^1 , -B^3 ,C^5 ,*~C^d , *C^7 , -C^1 , C^f ,-~C^9 ,*~C^3 , -C^b ,D^05,*~D^15, *D^13,-~D^03, D^19,-~D^17,*~D^01, -D^11,D^0f,*~D^1f, *D^09, -D^07,*D^1d,-~D^0d,*~D^0b, -D^1b}
0c	{1, *A^2, A^3,*~A^1,B^3 , *B^1 , *B^7 ,*~B^5 ,C^3 , *C^5 ,  C^f ,*~C^7 , C^9 ,*~C^1 , *C^b ,*~C^d ,D^03, *D^0d,  D^1b,*~D^0b, D^0f, *D^01, *D^19,*~D^17,D^09, *D^07, *D^1f,*~D^11, D^15,*~D^05, *D^13,*~D^1d}
1c	{1,-~A^2,*A^1,-*A^3,B^7 , -B^3 ,*~B^5 ,~*B^1 ,C^7 ,-~C^1 ,*~C^3 ,-*C^5 ,*C^b , -C^d ,-~C^f ,~*C^9 ,D^07,-~D^09, *D^01,-*D^11,*D^1d, -D^13,*~D^1b,~*D^0b,D^15, -D^05,*~D^0d,-*D^03,*D^0f,-*D^1f,-~D^17,~*D^19}
02	{1,  A^1, B^1,  B^3,C^1 ,  C^5 ,  C^3 ,  C^7 ,D^01,  D^09,  D^05,  D^0d, D^03,  D^0b,  D^07,  D^0f,E^01,  E^11,  E^09,  E^19, E^05,  E^15,  E^0d,  E^1d,E^03,  E^13,  E^0b,  E^1b, E^07,  E^17,  E^0f,  E^1f}
12	{1,*~A^1,*B^7,-~B^5,C^9 ,*~C^d , *C^5 ,-~C^1 ,D^09,*~D^11, *D^13,-~D^0b, D^1b,-~D^1d, *D^01, -D^07,E^09,*~E^19, *E^2f,-~E^1f, E^2d,*~E^3d, *E^0b, -E^05,E^1b,*~E^2b, *E^1d,-~E^0d, E^3f,-~E^31,*~E^07, -E^17}
0a	{1, *A^3, B^5, *B^1,C^5 , *C^7 ,  C^f ,*~C^3 ,D^05, *D^13,  D^19,*~D^01, D^0f, *D^09, *D^1d,*~D^0b,E^05, *E^2b,  E^2d, *E^03, E^19, *E^17, *E^3f,*~E^11,E^0f, *E^21,  E^37,*~E^07, E^23, *E^0d, *E^35,*~E^1b}
1a	{1,-~A^3,*B^3, -B^7,C^d , -C^1 ,*~C^7 ,-*C^5 ,D^0d,-~D^0b,*~D^01,-*D^17,*D^19, -D^0f,*~D^1b,~*D^03,E^0d,-~E^23, *E^0b,-*E^3b,*E^3f, -E^11,*~E^29,-*E^07,E^27,-~E^09,*~E^0f,-*E^21,*E^25, -E^2b,-~E^3d,~*E^13}
06	{1,  A^3, B^3, *B^7,C^3 ,  C^f ,  C^9 , *C^b ,D^03,  D^1b,  D^0f, *D^19, D^09, *D^1f,  D^15, *D^13,E^03,  E^33,  E^1b, *E^35, E^0f,  E^3f,  E^27, *E^29,E^09,  E^39,  E^21, *E^2f, E^15, *E^3b,  E^2d, *E^23}
16	{1,*~A^3,*B^5, -B^1,C^b ,-~C^9 ,*~C^1 , -C^d ,D^0b,-~D^1d, *D^09, -D^0f,*D^1f,-~D^07,*~D^0d,-*D^1b,E^0b,*~E^3b, *E^1d, -E^13, E^37,-~E^19,*~E^0f, -E^3f,E^21,-~E^2f, *E^07, -E^29,*E^33,-~E^03,*~E^25,-*E^2b}
0e	{1, *A^1, B^7,*~B^5,C^7 ,*~C^3 , *C^b ,-~C^f ,D^07, *D^01, *D^1d,*~D^1b, D^15,*~D^0d, *D^0f,-~D^17,E^07, *E^09,  E^3f,*~E^2f, E^23,*~E^13, *E^25,-~E^35,E^15,*~E^05, *E^33,*~E^3d, E^31,*~E^21, *E^17,-~E^27}
1e	{1,-~A^1,*B^1,-*B^3,C^f , -C^b ,*~C^d ,~*C^9 ,D^0f, -D^07,*~D^0b,~*D^03,*D^13,-*D^1b,-~D^17, ~D^1f,E^0f,-~E^01,*~E^07,-*E^09,*E^35, -E^3b,-~E^3d,~*E^33,E^2d, -E^1d,*~E^25,~*E^15,*E^17,-*E^27,-~E^1f, ~E^2f}
01	{1,  B^1, C^1,  C^3,D^01,  D^05,  D^03,  D^07,E^01,  E^09,  E^05,  E^0d, E^03,  E^0b,  E^07,  E^0f,F^01,  F^11,  F^09,  F^19, F^05,  F^15,  F^0d,  F^1d,F^03,  F^13,  F^0b,  F^1b, F^07,  F^17,  F^0f,  F^1f}
11	{1,*~B^1,*C^f,-~C^d,D^11,*~D^15, *D^0d,-~D^09,E^11,*~E^19, *E^2b,-~E^23, E^33,*~E^3b, *E^09,-~E^01,F^11,*~F^21, *F^67,-~F^57, F^55,*~F^65, *F^23,-~F^13,F^33,*~F^43, *F^45,-~F^35, F^77,-~F^79, *F^01, -F^0f}
09	{1, *B^7, C^9, *C^5,D^09, *D^13,  D^1b, *D^01,E^09, *E^2f,  E^2d, *E^0b, E^1b, *E^1d,  E^3f,*~E^07,F^09, *F^67,  F^51, *F^1f, F^2d, *F^43,  F^75,*~F^05,F^1b, *F^55,  F^63, *F^0d, F^3f, *F^31, *F^79,*~F^17}
19	{1,-~B^7,*C^7, -C^b,D^19,-~D^03,*~D^0b,-*D^11,E^19,-~E^1f, *E^03,-*E^3b,*E^35, -E^13,*~E^2f,-*E^09,F^19,-~F^57, *F^1f, -F^71, F^7d, -F^0d,*~F^45,-*F^2b,F^4b,-~F^25,*~F^13,-*F^5d,*F^51, -F^3f,*~F^77,~*F^ 7}
05	{1,  B^5, C^5,  C^f,D^05,  D^19,  D^0f, *D^1d,E^05,  E^2d,  E^19, *E^3f, E^0f,  E^37,  E^23, *E^35,F^05,  F^55,  F^2d,  F^7d, F^19,  F^69,  F^41, *F^6f,F^0f,  F^5f,  F^37, *F^79, F^23,  F^73,  F^4b, *F^65}
15	{1,*~B^5,*C^b,-~C^1,D^15,-~D^17, *D^01, -D^13,E^15,*~E^3d, *E^17, -E^11, E^3f,-~E^19,*~E^13, -E^3b,F^15,*~F^65, *F^43, -F^0d, F^69,-~F^47,*~F^11, -F^61,F^3f,-~F^71, *F^19, -F^37,*F^6d,-~F^1d,*~F^3b,-*F^75}
0d	{1, *B^3, C^d,*~C^7,D^0d,*~D^01, *D^19,*~D^1b,E^0d, *E^0b, *E^3f,*~E^29, E^27,*~E^0f, *E^25,-~E^3d,F^0d, *F^23,  F^75,*~F^45, F^41,*~F^11, *F^57,*~F^79,F^27, *F^09, *F^71,*~F^5f, F^5b,*~F^2b, *F^3d,-~F^6d}
1d	{1,-~B^3,*C^3,-*C^9,D^1d, -D^11,*~D^17,~*D^0b,E^1d, -E^05,*~E^11,-*E^07,*E^29, -E^3f,-~E^35,~*E^33,F^1d,-~F^13,*~F^05,-*F^2b,*F^6f, -F^61,*~F^79,~*F^49,F^57, -F^27,*~F^3f,~*F^0f,*F^35,-*F^65,-~F^4d, ~F^7d}
03	{1,  B^3, C^3,  C^9,D^03,  D^0f,  D^09,  D^15,E^03,  E^1b,  E^0f,  E^27, E^09,  E^21,  E^15,  E^2d,F^03,  F^33,  F^1b,  F^4b, F^0f,  F^3f,  F^27,  F^57,F^09,  F^39,  F^21,  F^51, F^15,  F^45,  F^2d,  F^5d}
13	{1,*~B^3,*C^d,-~C^7,D^13,*~D^1f, *D^07, -D^05,E^13,*~E^2b, *E^21,-~E^09, E^39,-~E^2f,*~E^05, -E^1d,F^13,*~F^43, *F^55,-~F^25, F^5f,-~F^71, *F^09, -F^27,F^39,*~F^69, *F^2f, -F^01,*F^7b,-~F^4b,*~F^1d, -F^4d}
0b	{1, *B^5, C^b,*~C^1,D^0b, *D^09, *D^1f,*~D^0d,E^0b, *E^1d,  E^37,*~E^0f, E^21, *E^07, *E^33,*~E^25,F^0b, *F^45,  F^63,*~F^13, F^37, *F^19, *F^71,*~F^3f,F^21, *F^2f,  F^79,*~F^29, F^4d, *F^03, *F^5b,*~F^55}
1b	{1,-~B^5,*C^5,-*C^f,D^1b, -D^07,*~D^11,-*D^03,E^1b,-~E^0d,*~E^07,-*E^21,*E^2f, -E^29,*~E^3d,~*E^15,F^1b,-~F^35, *F^0d,-*F^5d,*F^79, -F^37,*~F^5f,~*F^0f,F^51, -F^01,*~F^29,-*F^27,*F^43, -F^6d,-~F^6b,~*F^45}
07	{1,  B^7, C^7, *C^b,D^07, *D^1d,  D^15, *D^0f,E^07,  E^3f,  E^23, *E^25, E^15, *E^33,  E^31, *E^17,F^07,  F^77,  F^3f, *F^51, F^23, *F^6d,  F^5b, *F^35,F^15, *F^7b,  F^4d, *F^43, F^31, *F^5f,  F^69, *F^27}
17	{1,*~B^7,*C^9, -C^5,D^17,-~D^0d,*~D^05,-*D^1f,E^17,-~E^31, *E^0d, -E^2b,*E^3b,-~E^03,*~E^21,-*E^27,F^17,-~F^79, *F^31, -F^3f, F^73,-~F^1d,*~F^2b,-*F^65,F^45,-~F^4b, *F^03, -F^6d,*F^5f, -F^11,*~F^59,-*F^37}
0f	{1, *B^1, C^f,*~C^d,D^0f,*~D^0b, *D^13,-~D^17,E^0f,*~E^07, *E^35,-~E^3d, E^2d,*~E^25, *E^17,-~E^1f,F^0f, *F^01, *F^79,*~F^77, F^4b,*~F^3b, *F^3d,-~F^4d,F^2d,*~F^1d, *F^5b,-~F^6b, F^69,*~F^59, *F^1f,-~F^2f}
1f	{1,-~B^1,*C^1,-*C^3,D^1f, -D^1b,*~D^1d,~*D^19,E^1f, -E^17,*~E^1b,~*E^13,*E^23,-*E^2b,-~E^27, ~E^2f,F^1f, -F^0f,*~F^17,~*F^07,*F^65,-*F^75,-~F^6d, ~F^7d,F^5d, -F^4d,*~F^55,~*F^45,*F^27,-*F^37,-~F^2f, ~F^3f} .

Now re-express these in terms of the sincos consts defined above, using the following operator -> sincos replacements:

	~*	s,-c
	-*	-s,-c
	*~	-s,c
	*	s,c
	-~	-c,s
	-	-c,-s
	~	c,-s
		c,s

	 0      10        8            18         4             14           c             1c          2             12           a             1a          6             16           e             1e          1              11            9              19           5              15            d              1d           3              13            b              1b           7              17            f              1f
	--  ---------  ---------  -----------  ----------  ------------  -----------  ------------  ----------  ------------  -----------  ------------  ----------  ------------  -----------  ------------  -----------  -------------  ------------  -------------  -----------  -------------  ------------  -------------  -----------  -------------  ------------  -------------  -----------  -------------  ------------  -------------
00
10	{1,  0,1     ,[ 1,1]IRT2,[ 1, 1]IRT2 ,[c,s]16    ,[-s, c]16    ,[ s,c]16    ,[-c, s]16    ,[c,s]32_1  ,[-s, c]32_1  ,[ s,c]32_3  ,[-c, s]32_3  ,[c,s]32_3  ,[-s, c]32_3  ,[ s,c]32_1  ,[-c, s]32_1  ,[c,s]64_1   ,[-s, c]64_1   ,[ s,c]64_7   ,[-c, s]64_7   ,[c,s]64_5   ,[-s, c]64_5   ,[ s,c]64_3   ,[-c, s]64_3   ,[c,s]64_3   ,[-s, c]64_3   ,[ s,c]64_5   ,[-c, s]64_5   ,[c,s]64_7   ,[-s, c]64_7   ,[ s,c]64_1   ,[-c, s]64_1   }
08	{1,[ 1,1]IRT2,[c,s]16   ,[ s, c]16   ,[c,s]32_1  ,[ s, c]32_3  ,[ c,s]32_3  ,[ s, c]32_1  ,[c,s]64_1  ,[ s, c]64_7  ,[ c,s]64_5  ,[ s, c]64_3  ,[c,s]64_3  ,[ s, c]64_5  ,[ c,s]64_7  ,[ s, c]64_1  ,[c,s]128_1  ,[ s, c]128_f  ,[ c,s]128_9  ,[ s, c]128_7  ,[c,s]128_5  ,[ s, c]128_b  ,[ c,s]128_d  ,[ s, c]128_3  ,[c,s]128_3  ,[ s, c]128_d  ,[ c,s]128_b  ,[ s, c]128_5  ,[c,s]128_7  ,[ s, c]128_9  ,[ c,s]128_f  ,[ s, c]128_1  }
18	{1,[-1,1]IRT2,[s,c]16   ,[-c,-s]16   ,[c,s]32_3  ,[-c, s]32_1  ,[-s,c]32_1  ,[-s,-c]32_3  ,[c,s]64_3  ,[-c, s]64_5  ,[ s,c]64_1  ,[-c,-s]64_7  ,[s,c]64_7  ,[-c,-s]64_1  ,[-s,c]64_5  ,[-s,-c]64_3  ,[c,s]128_3  ,[-c, s]128_d  ,[ s,c]128_5  ,[-c,-s]128_b  ,[c,s]128_f  ,[-c, s]128_1  ,[-s,c]128_7  ,[-s,-c]128_9  ,[c,s]128_9  ,[-c, s]128_7  ,[-s,c]128_1  ,[-s,-c]128_f  ,[s,c]128_b  ,[-c,-s]128_5  ,[-s,c]128_d  ,[-s,-c]128_3  }
04	{1,[ c,s]16  ,[c,s]32_1 ,[ c, s]32_3 ,[c,s]64_1  ,[ c, s]64_5  ,[ c,s]64_3  ,[ c, s]64_7  ,[c,s]128_1 ,[ c, s]128_9 ,[ c,s]128_5 ,[ c, s]128_d ,[c,s]128_3 ,[ c, s]128_b ,[ c,s]128_7 ,[ c, s]128_f ,[c,s]256_01 ,[ c, s]256_11 ,[ c,s]256_09 ,[ c, s]256_19 ,[c,s]256_05 ,[ c, s]256_15 ,[ c,s]256_0d ,[ c, s]256_1d ,[c,s]256_03 ,[ c, s]256_13 ,[ c,s]256_0b ,[ c, s]256_1b ,[c,s]256_07 ,[ c, s]256_17 ,[ c,s]256_0f ,[ c, s]256_1f }
14	{1,[-s,c]16  ,[s,c]32_3 ,[-c, s]32_1 ,[c,s]64_5  ,[-c, s]64_7  ,[ s,c]64_1  ,[-c,-s]64_3  ,[c,s]128_5 ,[-s, c]128_d ,[ s,c]128_7 ,[-c,-s]128_1 ,[c,s]128_f ,[-c, s]128_9 ,[-s,c]128_3 ,[-c,-s]128_b ,[c,s]256_05 ,[-s, c]256_15 ,[ s,c]256_13 ,[-c, s]256_03 ,[c,s]256_19 ,[-c, s]256_17 ,[-s,c]256_01 ,[-c,-s]256_11 ,[c,s]256_0f ,[-s, c]256_1f ,[ s,c]256_09 ,[-c,-s]256_07 ,[s,c]256_1d ,[-c, s]256_0d ,[-s,c]256_0b ,[-c,-s]256_1b }
0c	{1,[ s,c]16  ,[c,s]32_3 ,[-s, c]32_1 ,[c,s]64_3  ,[ s, c]64_1  ,[ s,c]64_7  ,[-s, c]64_5  ,[c,s]128_3 ,[ s, c]128_5 ,[ c,s]128_f ,[-s, c]128_7 ,[c,s]128_9 ,[-s, c]128_1 ,[ s,c]128_b ,[-s, c]128_d ,[c,s]256_03 ,[ s, c]256_0d ,[ c,s]256_1b ,[-s, c]256_0b ,[c,s]256_0f ,[ s, c]256_01 ,[ s,c]256_19 ,[-s, c]256_17 ,[c,s]256_09 ,[ s, c]256_07 ,[ s,c]256_1f ,[-s, c]256_11 ,[c,s]256_15 ,[-s, c]256_05 ,[ s,c]256_13 ,[-s, c]256_1d }
1c	{1,[-c,s]16  ,[s,c]32_1 ,[-s,-c]32_3 ,[c,s]64_7  ,[-c,-s]64_3  ,[-s,c]64_5  ,[ s,-c]64_1  ,[c,s]128_7 ,[-c, s]128_1 ,[-s,c]128_3 ,[-s,-c]128_5 ,[s,c]128_b ,[-c,-s]128_d ,[-c,s]128_f ,[ s,-c]128_9 ,[c,s]256_07 ,[-c, s]256_09 ,[ s,c]256_01 ,[-s,-c]256_11 ,[s,c]256_1d ,[-c,-s]256_13 ,[-s,c]256_1b ,[ s,-c]256_0b ,[c,s]256_15 ,[-c,-s]256_05 ,[-s,c]256_0d ,[-s,-c]256_03 ,[s,c]256_0f ,[-s,-c]256_1f ,[-c,s]256_17 ,[ s,-c]256_19 }
02	{1,[ c,s]32_1,[c,s]64_1 ,[ c, s]64_3 ,[c,s]128_1 ,[ c, s]128_5 ,[ c,s]128_3 ,[ c, s]128_7 ,[c,s]256_01,[ c, s]256_09,[ c,s]256_05,[ c, s]256_0d,[c,s]256_03,[ c, s]256_0b,[ c,s]256_07,[ c, s]256_0f,[c,s]512_01 ,[ c, s]512_11 ,[ c,s]512_09 ,[ c, s]512_19 ,[c,s]512_05 ,[ c, s]512_15 ,[ c,s]512_0d ,[ c, s]512_1d ,[c,s]512_03 ,[ c, s]512_13 ,[ c,s]512_0b ,[ c, s]512_1b ,[c,s]512_07 ,[ c, s]512_17 ,[ c,s]512_0f ,[ c, s]512_1f }
12	{1,[-s,c]32_1,[s,c]64_7 ,[-c, s]64_5 ,[c,s]128_9 ,[-s, c]128_d ,[ s,c]128_5 ,[-c, s]128_1 ,[c,s]256_09,[-s, c]256_11,[ s,c]256_13,[-c, s]256_0b,[c,s]256_1b,[-c, s]256_1d,[ s,c]256_01,[-c,-s]256_07,[c,s]512_09 ,[-s, c]512_19 ,[ s,c]512_2f ,[-c, s]512_1f ,[c,s]512_2d ,[-s, c]512_3d ,[ s,c]512_0b ,[-c,-s]512_05 ,[c,s]512_1b ,[-s, c]512_2b ,[ s,c]512_1d ,[-c, s]512_0d ,[c,s]512_3f ,[-c, s]512_31 ,[-s,c]512_07 ,[-c,-s]512_17 }
0a	{1,[ s,c]32_3,[c,s]64_5 ,[ s, c]64_1 ,[c,s]128_5 ,[ s, c]128_7 ,[ c,s]128_f ,[-s, c]128_3 ,[c,s]256_05,[ s, c]256_13,[ c,s]256_19,[-s, c]256_01,[c,s]256_0f,[ s, c]256_09,[ s,c]256_1d,[-s, c]256_0b,[c,s]512_05 ,[ s, c]512_2b ,[ c,s]512_2d ,[ s, c]512_03 ,[c,s]512_19 ,[ s, c]512_17 ,[ s,c]512_3f ,[-s, c]512_11 ,[c,s]512_0f ,[ s, c]512_21 ,[ c,s]512_37 ,[-s, c]512_07 ,[c,s]512_23 ,[ s, c]512_0d ,[ s,c]512_35 ,[-s, c]512_1b }
1a	{1,[-c,s]32_3,[s,c]64_3 ,[-c,-s]64_7 ,[c,s]128_d ,[-c,-s]128_1 ,[-s,c]128_7 ,[-s,-c]128_5 ,[c,s]256_0d,[-c, s]256_0b,[-s,c]256_01,[-s,-c]256_17,[s,c]256_19,[-c,-s]256_0f,[-s,c]256_1b,[ s,-c]256_03,[c,s]512_0d ,[-c, s]512_23 ,[ s,c]512_0b ,[-s,-c]512_3b ,[s,c]512_3f ,[-c,-s]512_11 ,[-s,c]512_29 ,[-s,-c]512_07 ,[c,s]512_27 ,[-c, s]512_09 ,[-s,c]512_0f ,[-s,-c]512_21 ,[s,c]512_25 ,[-c,-s]512_2b ,[-c,s]512_3d ,[ s,-c]512_13 }
06	{1,[ c,s]32_3,[c,s]64_3 ,[ s, c]64_7 ,[c,s]128_3 ,[ c, s]128_f ,[ c,s]128_9 ,[ s, c]128_b ,[c,s]256_03,[ c, s]256_1b,[ c,s]256_0f,[ s, c]256_19,[c,s]256_09,[ s, c]256_1f,[ c,s]256_15,[ s, c]256_13,[c,s]512_03 ,[ c, s]512_33 ,[ c,s]512_1b ,[ s, c]512_35 ,[c,s]512_0f ,[ c, s]512_3f ,[ c,s]512_27 ,[ s, c]512_29 ,[c,s]512_09 ,[ c, s]512_39 ,[ c,s]512_21 ,[ s, c]512_2f ,[c,s]512_15 ,[ s, c]512_3b ,[ c,s]512_2d ,[ s, c]512_23 }
16	{1,[-s,c]32_3,[s,c]64_5 ,[-c,-s]64_1 ,[c,s]128_b ,[-c, s]128_9 ,[-s,c]128_1 ,[-c,-s]128_d ,[c,s]256_0b,[-c, s]256_1d,[ s,c]256_09,[-c,-s]256_0f,[s,c]256_1f,[-c, s]256_07,[-s,c]256_0d,[-s,-c]256_1b,[c,s]512_0b ,[-s, c]512_3b ,[ s,c]512_1d ,[-c,-s]512_13 ,[c,s]512_37 ,[-c, s]512_19 ,[-s,c]512_0f ,[-c,-s]512_3f ,[c,s]512_21 ,[-c, s]512_2f ,[ s,c]512_07 ,[-c,-s]512_29 ,[s,c]512_33 ,[-c, s]512_03 ,[-s,c]512_25 ,[-s,-c]512_2b }
0e	{1,[ s,c]32_1,[c,s]64_7 ,[-s, c]64_5 ,[c,s]128_7 ,[-s, c]128_3 ,[ s,c]128_b ,[-c, s]128_f ,[c,s]256_07,[ s, c]256_01,[ s,c]256_1d,[-s, c]256_1b,[c,s]256_15,[-s, c]256_0d,[ s,c]256_0f,[-c, s]256_17,[c,s]512_07 ,[ s, c]512_09 ,[ c,s]512_3f ,[-s, c]512_2f ,[c,s]512_23 ,[-s, c]512_13 ,[ s,c]512_25 ,[-c, s]512_35 ,[c,s]512_15 ,[-s, c]512_05 ,[ s,c]512_33 ,[-s, c]512_3d ,[c,s]512_31 ,[-s, c]512_21 ,[ s,c]512_17 ,[-c, s]512_27 }
1e	{1,[-c,s]32_1,[s,c]64_1 ,[-s,-c]64_3 ,[c,s]128_f ,[-c,-s]128_b ,[-s,c]128_d ,[ s,-c]128_9 ,[c,s]256_0f,[-c,-s]256_07,[-s,c]256_0b,[ s,-c]256_03,[s,c]256_13,[-s,-c]256_1b,[-c,s]256_17,[ c,-s]256_1f,[c,s]512_0f ,[-c, s]512_01 ,[-s,c]512_07 ,[-s,-c]512_09 ,[s,c]512_35 ,[-c,-s]512_3b ,[-c,s]512_3d ,[ s,-c]512_33 ,[c,s]512_2d ,[-c,-s]512_1d ,[-s,c]512_25 ,[ s,-c]512_15 ,[s,c]512_17 ,[-s,-c]512_27 ,[-c,s]512_1f ,[ c,-s]512_2f }
01	{1,[ c,s]64_1,[c,s]128_1,[ c, s]128_3,[c,s]256_01,[ c, s]256_05,[ c,s]256_03,[ c, s]256_07,[c,s]512_01,[ c, s]512_09,[ c,s]512_05,[ c, s]512_0d,[c,s]512_03,[ c, s]512_0b,[ c,s]512_07,[ c, s]512_0f,[c,s]1024_01,[ c, s]1024_11,[ c,s]1024_09,[ c, s]1024_19,[c,s]1024_05,[ c, s]1024_15,[ c,s]1024_0d,[ c, s]1024_1d,[c,s]1024_03,[ c, s]1024_13,[ c,s]1024_0b,[ c, s]1024_1b,[c,s]1024_07,[ c, s]1024_17,[ c,s]1024_0f,[ c, s]1024_1f}
11	{1,[-s,c]64_1,[s,c]128_f,[-c, s]128_d,[c,s]256_11,[-s, c]256_15,[ s,c]256_0d,[-c, s]256_09,[c,s]512_11,[-s, c]512_19,[ s,c]512_2b,[-c, s]512_23,[c,s]512_33,[-s, c]512_3b,[ s,c]512_09,[-c, s]512_01,[c,s]1024_11,[-s, c]1024_21,[ s,c]1024_67,[-c, s]1024_57,[c,s]1024_55,[-s, c]1024_65,[ s,c]1024_23,[-c, s]1024_13,[c,s]1024_33,[-s, c]1024_43,[ s,c]1024_45,[-c, s]1024_35,[c,s]1024_77,[-c, s]1024_79,[ s,c]1024_01,[-c,-s]1024_0f}
09	{1,[ s,c]64_7,[c,s]128_9,[ s, c]128_5,[c,s]256_09,[ s, c]256_13,[ c,s]256_1b,[ s, c]256_01,[c,s]512_09,[ s, c]512_2f,[ c,s]512_2d,[ s, c]512_0b,[c,s]512_1b,[ s, c]512_1d,[ c,s]512_3f,[-s, c]512_07,[c,s]1024_09,[ s, c]1024_67,[ c,s]1024_51,[ s, c]1024_1f,[c,s]1024_2d,[ s, c]1024_43,[ c,s]1024_75,[-s, c]1024_05,[c,s]1024_1b,[ s, c]1024_55,[ c,s]1024_63,[ s, c]1024_0d,[c,s]1024_3f,[ s, c]1024_31,[ s,c]1024_79,[-s, c]1024_17}
19	{1,[-c,s]64_7,[s,c]128_7,[-c,-s]128_b,[c,s]256_19,[-c, s]256_03,[-s,c]256_0b,[-s,-c]256_11,[c,s]512_19,[-c, s]512_1f,[ s,c]512_03,[-s,-c]512_3b,[s,c]512_35,[-c,-s]512_13,[-s,c]512_2f,[-s,-c]512_09,[c,s]1024_19,[-c, s]1024_57,[ s,c]1024_1f,[-c,-s]1024_71,[c,s]1024_7d,[-c,-s]1024_0d,[-s,c]1024_45,[-s,-c]1024_2b,[c,s]1024_4b,[-c, s]1024_25,[-s,c]1024_13,[-s,-c]1024_5d,[s,c]1024_51,[-c,-s]1024_3f,[-s,c]1024_77,[ s,-c]1024_07}
05	{1,[ c,s]64_5,[c,s]128_5,[ c, s]128_f,[c,s]256_05,[ c, s]256_19,[ c,s]256_0f,[ s, c]256_1d,[c,s]512_05,[ c, s]512_2d,[ c,s]512_19,[ s, c]512_3f,[c,s]512_0f,[ c, s]512_37,[ c,s]512_23,[ s, c]512_35,[c,s]1024_05,[ c, s]1024_55,[ c,s]1024_2d,[ c, s]1024_7d,[c,s]1024_19,[ c, s]1024_69,[ c,s]1024_41,[ s, c]1024_6f,[c,s]1024_0f,[ c, s]1024_5f,[ c,s]1024_37,[ s, c]1024_79,[c,s]1024_23,[ c, s]1024_73,[ c,s]1024_4b,[ s, c]1024_65}
15	{1,[-s,c]64_5,[s,c]128_b,[-c, s]128_1,[c,s]256_15,[-c, s]256_17,[ s,c]256_01,[-c,-s]256_13,[c,s]512_15,[-s, c]512_3d,[ s,c]512_17,[-c,-s]512_11,[c,s]512_3f,[-c, s]512_19,[-s,c]512_13,[-c,-s]512_3b,[c,s]1024_15,[-s, c]1024_65,[ s,c]1024_43,[-c,-s]1024_0d,[c,s]1024_69,[-c, s]1024_47,[-s,c]1024_11,[-c,-s]1024_61,[c,s]1024_3f,[-c, s]1024_71,[ s,c]1024_19,[-c,-s]1024_37,[s,c]1024_6d,[-c, s]1024_1d,[-s,c]1024_3b,[-s,-c]1024_75}
0d	{1,[ s,c]64_3,[c,s]128_d,[-s, c]128_7,[c,s]256_0d,[-s, c]256_01,[ s,c]256_19,[-s, c]256_1b,[c,s]512_0d,[ s, c]512_0b,[ s,c]512_3f,[-s, c]512_29,[c,s]512_27,[-s, c]512_0f,[ s,c]512_25,[-c, s]512_3d,[c,s]1024_0d,[ s, c]1024_23,[ c,s]1024_75,[-s, c]1024_45,[c,s]1024_41,[-s, c]1024_11,[ s,c]1024_57,[-s, c]1024_79,[c,s]1024_27,[ s, c]1024_09,[ s,c]1024_71,[-s, c]1024_5f,[c,s]1024_5b,[-s, c]1024_2b,[ s,c]1024_3d,[-c, s]1024_6d}
1d	{1,[-c,s]64_3,[s,c]128_3,[-s,-c]128_9,[c,s]256_1d,[-c,-s]256_11,[-s,c]256_17,[ s,-c]256_0b,[c,s]512_1d,[-c,-s]512_05,[-s,c]512_11,[-s,-c]512_07,[s,c]512_29,[-c,-s]512_3f,[-c,s]512_35,[ s,-c]512_33,[c,s]1024_1d,[-c, s]1024_13,[-s,c]1024_05,[-s,-c]1024_2b,[s,c]1024_6f,[-c,-s]1024_61,[-s,c]1024_79,[ s,-c]1024_49,[c,s]1024_57,[-c,-s]1024_27,[-s,c]1024_3f,[ s,-c]1024_0f,[s,c]1024_35,[-s,-c]1024_65,[-c,s]1024_4d,[ c,-s]1024_7d}
03	{1,[ c,s]64_3,[c,s]128_3,[ c, s]128_9,[c,s]256_03,[ c, s]256_0f,[ c,s]256_09,[ c, s]256_15,[c,s]512_03,[ c, s]512_1b,[ c,s]512_0f,[ c, s]512_27,[c,s]512_09,[ c, s]512_21,[ c,s]512_15,[ c, s]512_2d,[c,s]1024_03,[ c, s]1024_33,[ c,s]1024_1b,[ c, s]1024_4b,[c,s]1024_0f,[ c, s]1024_3f,[ c,s]1024_27,[ c, s]1024_57,[c,s]1024_09,[ c, s]1024_39,[ c,s]1024_21,[ c, s]1024_51,[c,s]1024_15,[ c, s]1024_45,[ c,s]1024_2d,[ c, s]1024_5d}
13	{1,[-s,c]64_3,[s,c]128_d,[-c, s]128_7,[c,s]256_13,[-s, c]256_1f,[ s,c]256_07,[-c,-s]256_05,[c,s]512_13,[-s, c]512_2b,[ s,c]512_21,[-c, s]512_09,[c,s]512_39,[-c, s]512_2f,[-s,c]512_05,[-c,-s]512_1d,[c,s]1024_13,[-s, c]1024_43,[ s,c]1024_55,[-c, s]1024_25,[c,s]1024_5f,[-c, s]1024_71,[ s,c]1024_09,[-c,-s]1024_27,[c,s]1024_39,[-s, c]1024_69,[ s,c]1024_2f,[-c,-s]1024_01,[s,c]1024_7b,[-c, s]1024_4b,[-s,c]1024_1d,[-c,-s]1024_4d}
0b	{1,[ s,c]64_5,[c,s]128_b,[-s, c]128_1,[c,s]256_0b,[ s, c]256_09,[ s,c]256_1f,[-s, c]256_0d,[c,s]512_0b,[ s, c]512_1d,[ c,s]512_37,[-s, c]512_0f,[c,s]512_21,[ s, c]512_07,[ s,c]512_33,[-s, c]512_25,[c,s]1024_0b,[ s, c]1024_45,[ c,s]1024_63,[-s, c]1024_13,[c,s]1024_37,[ s, c]1024_19,[ s,c]1024_71,[-s, c]1024_3f,[c,s]1024_21,[ s, c]1024_2f,[ c,s]1024_79,[-s, c]1024_29,[c,s]1024_4d,[ s, c]1024_03,[ s,c]1024_5b,[-s, c]1024_55}
1b	{1,[-c,s]64_5,[s,c]128_5,[-s,-c]128_f,[c,s]256_1b,[-c,-s]256_07,[-s,c]256_11,[-s,-c]256_03,[c,s]512_1b,[-c, s]512_0d,[-s,c]512_07,[-s,-c]512_21,[s,c]512_2f,[-c,-s]512_29,[-s,c]512_3d,[ s,-c]512_15,[c,s]1024_1b,[-c, s]1024_35,[ s,c]1024_0d,[-s,-c]1024_5d,[s,c]1024_79,[-c,-s]1024_37,[-s,c]1024_5f,[ s,-c]1024_0f,[c,s]1024_51,[-c,-s]1024_01,[-s,c]1024_29,[-s,-c]1024_27,[s,c]1024_43,[-c,-s]1024_6d,[-c,s]1024_6b,[ s,-c]1024_45}
07	{1,[ c,s]64_7,[c,s]128_7,[ s, c]128_b,[c,s]256_07,[ s, c]256_1d,[ c,s]256_15,[ s, c]256_0f,[c,s]512_07,[ c, s]512_3f,[ c,s]512_23,[ s, c]512_25,[c,s]512_15,[ s, c]512_33,[ c,s]512_31,[ s, c]512_17,[c,s]1024_07,[ c, s]1024_77,[ c,s]1024_3f,[ s, c]1024_51,[c,s]1024_23,[ s, c]1024_6d,[ c,s]1024_5b,[ s, c]1024_35,[c,s]1024_15,[ s, c]1024_7b,[ c,s]1024_4d,[ s, c]1024_43,[c,s]1024_31,[ s, c]1024_5f,[ c,s]1024_69,[ s, c]1024_27}
17	{1,[-s,c]64_7,[s,c]128_9,[-c,-s]128_5,[c,s]256_17,[-c, s]256_0d,[-s,c]256_05,[-s,-c]256_1f,[c,s]512_17,[-c, s]512_31,[ s,c]512_0d,[-c,-s]512_2b,[s,c]512_3b,[-c, s]512_03,[-s,c]512_21,[-s,-c]512_27,[c,s]1024_17,[-c, s]1024_79,[ s,c]1024_31,[-c,-s]1024_3f,[c,s]1024_73,[-c, s]1024_1d,[-s,c]1024_2b,[-s,-c]1024_65,[c,s]1024_45,[-c, s]1024_4b,[ s,c]1024_03,[-c,-s]1024_6d,[s,c]1024_5f,[-c,-s]1024_11,[-s,c]1024_59,[-s,-c]1024_37}
0f	{1,[ s,c]64_1,[c,s]128_f,[-s, c]128_d,[c,s]256_0f,[-s, c]256_0b,[ s,c]256_13,[-c, s]256_17,[c,s]512_0f,[-s, c]512_07,[ s,c]512_35,[-c, s]512_3d,[c,s]512_2d,[-s, c]512_25,[ s,c]512_17,[-c, s]512_1f,[c,s]1024_0f,[ s, c]1024_01,[ s,c]1024_79,[-s, c]1024_77,[c,s]1024_4b,[-s, c]1024_3b,[ s,c]1024_3d,[-c, s]1024_4d,[c,s]1024_2d,[-s, c]1024_1d,[ s,c]1024_5b,[-c, s]1024_6b,[c,s]1024_69,[-s, c]1024_59,[ s,c]1024_1f,[-c, s]1024_2f}
1f	{1,[-c,s]64_1,[s,c]128_1,[-s,-c]128_3,[c,s]256_1f,[-c,-s]256_1b,[-s,c]256_1d,[ s,-c]256_19,[c,s]512_1f,[-c,-s]512_17,[-s,c]512_1b,[ s,-c]512_13,[s,c]512_23,[-s,-c]512_2b,[-c,s]512_27,[ c,-s]512_2f,[c,s]1024_1f,[-c,-s]1024_0f,[-s,c]1024_17,[ s,-c]1024_07,[s,c]1024_65,[-s,-c]1024_75,[-c,s]1024_6d,[ c,-s]1024_7d,[c,s]1024_5d,[-c,-s]1024_4d,[-s,c]1024_55,[ s,-c]1024_45,[s,c]1024_27,[-s,-c]1024_37,[-c,s]1024_2f,[ c,-s]1024_3f} .

	Only the last 31 inputs to each of the radix-32 transforms ("Blocks") 1 through 31 are multiplied by non-unity twiddles.
	For DIF we process both the blocks, and the twiddles within each block, in bit-reversed order.
	One can see from the data below that aside from using a twiddleless DIF for Block 0 there is
	little to be gained from trying to exploit other "special form" twiddles such as I and isrt2*[+-1,+-1].
	Thus our radix-32-DIF-with-twiddles macro uses generic complex MUL for the 31 non-unity twiddles of each invocation.
	*/
// First 15 (of 31 nontrivial) of each twiddles set must match those of the analogous radix-16 DFT call in the radix-512 DIF.

	/* Block 0: has all-unity twiddles: */		jt = j1;
		RADIX_32_DIF(
			(double *)(t+ 0),p_offsets,1, (a+jt),q_offsets,RE_IM_STRIDE
		);
	/* Block 10: */		jt = j1 + p20;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+ 1),p_offsets, (a+jt),q_offsets,
			0,1,ISRT2,ISRT2,-ISRT2,ISRT2,c16,s16,-s16,c16,s16,c16,-c16,s16,c32_1,s32_1,-s32_1,c32_1,s32_3,c32_3,-c32_3,s32_3,c32_3,s32_3,-s32_3,c32_3,s32_1,c32_1,-c32_1,s32_1,
			c64_1,s64_1,-s64_1,c64_1,s64_7,c64_7,-c64_7,s64_7,c64_5,s64_5,-s64_5,c64_5,s64_3,c64_3,-c64_3,s64_3,c64_3,s64_3,-s64_3,c64_3,s64_5,c64_5,-c64_5,s64_5,c64_7,s64_7,-s64_7,c64_7,s64_1,c64_1,-c64_1,s64_1
		);
	/* Block 8: */		jt = j1 + p40;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+ 2),p_offsets, (a+jt),q_offsets,
			ISRT2,ISRT2,c16,s16,s16,c16,c32_1,s32_1,s32_3,c32_3,c32_3,s32_3,s32_1,c32_1,c64_1,s64_1,s64_7,c64_7,c64_5,s64_5,s64_3,c64_3,c64_3,s64_3,s64_5,c64_5,c64_7,s64_7,s64_1,c64_1,
			c128_1,s128_1,s128_f,c128_f,c128_9,s128_9,s128_7,c128_7,c128_5,s128_5,s128_b,c128_b,c128_d,s128_d,s128_3,c128_3,c128_3,s128_3,s128_d,c128_d,c128_b,s128_b,s128_5,c128_5,c128_7,s128_7,s128_9,c128_9,c128_f,s128_f,s128_1,c128_1
		);
	/* Block 18: */		jt = j1 + p60;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+ 3),p_offsets, (a+jt),q_offsets,
			-ISRT2,ISRT2,s16,c16,-c16,-s16,c32_3,s32_3,-c32_1,s32_1,-s32_1,c32_1,-s32_3,-c32_3,c64_3,s64_3,-c64_5,s64_5,s64_1,c64_1,-c64_7,-s64_7,s64_7,c64_7,-c64_1,-s64_1,-s64_5,c64_5,-s64_3,-c64_3,
			c128_3,s128_3,-c128_d,s128_d,s128_5,c128_5,-c128_b,-s128_b,c128_f,s128_f,-c128_1,s128_1,-s128_7,c128_7,-s128_9,-c128_9,c128_9,s128_9,-c128_7,s128_7,-s128_1,c128_1,-s128_f,-c128_f,s128_b,c128_b,-c128_5,-s128_5,-s128_d,c128_d,-s128_3,-c128_3
		);
	/* Block 4: */		jt = j1 + p80;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+ 4),p_offsets, (a+jt),q_offsets,
			c16,s16,c32_1,s32_1,c32_3,s32_3,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7,c128_1,s128_1,c128_9,s128_9,c128_5,s128_5,c128_d,s128_d,c128_3,s128_3,c128_b,s128_b,c128_7,s128_7,c128_f,s128_f,
			c256_01,s256_01,c256_11,s256_11,c256_09,s256_09,c256_19,s256_19,c256_05,s256_05,c256_15,s256_15,c256_0d,s256_0d,c256_1d,s256_1d,c256_03,s256_03,c256_13,s256_13,c256_0b,s256_0b,c256_1b,s256_1b,c256_07,s256_07,c256_17,s256_17,c256_0f,s256_0f,c256_1f,s256_1f
		);
	/* Block 14: */		jt = j1 + pa0;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+ 5),p_offsets, (a+jt),q_offsets,
			-s16,c16,s32_3,c32_3,-c32_1,s32_1,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3,c128_5,s128_5,-s128_d,c128_d,s128_7,c128_7,-c128_1,-s128_1,c128_f,s128_f,-c128_9,s128_9,-s128_3,c128_3,-c128_b,-s128_b,
			c256_05,s256_05,-s256_15,c256_15,s256_13,c256_13,-c256_03,s256_03,c256_19,s256_19,-c256_17,s256_17,-s256_01,c256_01,-c256_11,-s256_11,c256_0f,s256_0f,-s256_1f,c256_1f,s256_09,c256_09,-c256_07,-s256_07,s256_1d,c256_1d,-c256_0d,s256_0d,-s256_0b,c256_0b,-c256_1b,-s256_1b
		);
	/* Block c: */		jt = j1 + pc0;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+ 6),p_offsets, (a+jt),q_offsets,
			s16,c16,c32_3,s32_3,-s32_1,c32_1,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5,c128_3,s128_3,s128_5,c128_5,c128_f,s128_f,-s128_7,c128_7,c128_9,s128_9,-s128_1,c128_1,s128_b,c128_b,-s128_d,c128_d,
			c256_03,s256_03,s256_0d,c256_0d,c256_1b,s256_1b,-s256_0b,c256_0b,c256_0f,s256_0f,s256_01,c256_01,s256_19,c256_19,-s256_17,c256_17,c256_09,s256_09,s256_07,c256_07,s256_1f,c256_1f,-s256_11,c256_11,c256_15,s256_15,-s256_05,c256_05,s256_13,c256_13,-s256_1d,c256_1d
		);
	/* Block 1c: */		jt = j1 + pe0;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+ 7),p_offsets, (a+jt),q_offsets,
			-c16,s16,s32_1,c32_1,-s32_3,-c32_3,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1,c128_7,s128_7,-c128_1,s128_1,-s128_3,c128_3,-s128_5,-c128_5,s128_b,c128_b,-c128_d,-s128_d,-c128_f,s128_f,s128_9,-c128_9,
			c256_07,s256_07,-c256_09,s256_09,s256_01,c256_01,-s256_11,-c256_11,s256_1d,c256_1d,-c256_13,-s256_13,-s256_1b,c256_1b,s256_0b,-c256_0b,c256_15,s256_15,-c256_05,-s256_05,-s256_0d,c256_0d,-s256_03,-c256_03,s256_0f,c256_0f,-s256_1f,-c256_1f,-c256_17,s256_17,s256_19,-c256_19
		);
	/* Block 2: */		jt = j1 + p100;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+ 8),p_offsets, (a+jt),q_offsets,
			c32_1,s32_1,c64_1,s64_1,c64_3,s64_3,c128_1,s128_1,c128_5,s128_5,c128_3,s128_3,c128_7,s128_7,c256_01,s256_01,c256_09,s256_09,c256_05,s256_05,c256_0d,s256_0d,c256_03,s256_03,c256_0b,s256_0b,c256_07,s256_07,c256_0f,s256_0f,
			c512_01,s512_01,c512_11,s512_11,c512_09,s512_09,c512_19,s512_19,c512_05,s512_05,c512_15,s512_15,c512_0d,s512_0d,c512_1d,s512_1d,c512_03,s512_03,c512_13,s512_13,c512_0b,s512_0b,c512_1b,s512_1b,c512_07,s512_07,c512_17,s512_17,c512_0f,s512_0f,c512_1f,s512_1f
		);
	/* Block 12: */		jt = j1 + p120;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+ 9),p_offsets, (a+jt),q_offsets,
			-s32_1,c32_1,s64_7,c64_7,-c64_5,s64_5,c128_9,s128_9,-s128_d,c128_d,s128_5,c128_5,-c128_1,s128_1,c256_09,s256_09,-s256_11,c256_11,s256_13,c256_13,-c256_0b,s256_0b,c256_1b,s256_1b,-c256_1d,s256_1d,s256_01,c256_01,-c256_07,-s256_07,
			c512_09,s512_09,-s512_19,c512_19,s512_2f,c512_2f,-c512_1f,s512_1f,c512_2d,s512_2d,-s512_3d,c512_3d,s512_0b,c512_0b,-c512_05,-s512_05,c512_1b,s512_1b,-s512_2b,c512_2b,s512_1d,c512_1d,-c512_0d,s512_0d,c512_3f,s512_3f,-c512_31,s512_31,-s512_07,c512_07,-c512_17,-s512_17
		);
	/* Block a: */		jt = j1 + p140;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+10),p_offsets, (a+jt),q_offsets,
			s32_3,c32_3,c64_5,s64_5,s64_1,c64_1,c128_5,s128_5,s128_7,c128_7,c128_f,s128_f,-s128_3,c128_3,c256_05,s256_05,s256_13,c256_13,c256_19,s256_19,-s256_01,c256_01,c256_0f,s256_0f,s256_09,c256_09,s256_1d,c256_1d,-s256_0b,c256_0b,
			c512_05,s512_05,s512_2b,c512_2b,c512_2d,s512_2d,s512_03,c512_03,c512_19,s512_19,s512_17,c512_17,s512_3f,c512_3f,-s512_11,c512_11,c512_0f,s512_0f,s512_21,c512_21,c512_37,s512_37,-s512_07,c512_07,c512_23,s512_23,s512_0d,c512_0d,s512_35,c512_35,-s512_1b,c512_1b
		);
	/* Block 1a: */		jt = j1 + p160;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+11),p_offsets, (a+jt),q_offsets,
			-c32_3,s32_3,s64_3,c64_3,-c64_7,-s64_7,c128_d,s128_d,-c128_1,-s128_1,-s128_7,c128_7,-s128_5,-c128_5,c256_0d,s256_0d,-c256_0b,s256_0b,-s256_01,c256_01,-s256_17,-c256_17,s256_19,c256_19,-c256_0f,-s256_0f,-s256_1b,c256_1b,s256_03,-c256_03,
			c512_0d,s512_0d,-c512_23,s512_23,s512_0b,c512_0b,-s512_3b,-c512_3b,s512_3f,c512_3f,-c512_11,-s512_11,-s512_29,c512_29,-s512_07,-c512_07,c512_27,s512_27,-c512_09,s512_09,-s512_0f,c512_0f,-s512_21,-c512_21,s512_25,c512_25,-c512_2b,-s512_2b,-c512_3d,s512_3d,s512_13,-c512_13
		);
	/* Block 6: */		jt = j1 + p180;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+12),p_offsets, (a+jt),q_offsets,
			c32_3,s32_3,c64_3,s64_3,s64_7,c64_7,c128_3,s128_3,c128_f,s128_f,c128_9,s128_9,s128_b,c128_b,c256_03,s256_03,c256_1b,s256_1b,c256_0f,s256_0f,s256_19,c256_19,c256_09,s256_09,s256_1f,c256_1f,c256_15,s256_15,s256_13,c256_13,
			c512_03,s512_03,c512_33,s512_33,c512_1b,s512_1b,s512_35,c512_35,c512_0f,s512_0f,c512_3f,s512_3f,c512_27,s512_27,s512_29,c512_29,c512_09,s512_09,c512_39,s512_39,c512_21,s512_21,s512_2f,c512_2f,c512_15,s512_15,s512_3b,c512_3b,c512_2d,s512_2d,s512_23,c512_23
		);
	/* Block 16: */		jt = j1 + p1a0;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+13),p_offsets, (a+jt),q_offsets,
			-s32_3,c32_3,s64_5,c64_5,-c64_1,-s64_1,c128_b,s128_b,-c128_9,s128_9,-s128_1,c128_1,-c128_d,-s128_d,c256_0b,s256_0b,-c256_1d,s256_1d,s256_09,c256_09,-c256_0f,-s256_0f,s256_1f,c256_1f,-c256_07,s256_07,-s256_0d,c256_0d,-s256_1b,-c256_1b,
			c512_0b,s512_0b,-s512_3b,c512_3b,s512_1d,c512_1d,-c512_13,-s512_13,c512_37,s512_37,-c512_19,s512_19,-s512_0f,c512_0f,-c512_3f,-s512_3f,c512_21,s512_21,-c512_2f,s512_2f,s512_07,c512_07,-c512_29,-s512_29,s512_33,c512_33,-c512_03,s512_03,-s512_25,c512_25,-s512_2b,-c512_2b
		);
	/* Block e: */		jt = j1 + p1c0;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+14),p_offsets, (a+jt),q_offsets,
			s32_1,c32_1,c64_7,s64_7,-s64_5,c64_5,c128_7,s128_7,-s128_3,c128_3,s128_b,c128_b,-c128_f,s128_f,c256_07,s256_07,s256_01,c256_01,s256_1d,c256_1d,-s256_1b,c256_1b,c256_15,s256_15,-s256_0d,c256_0d,s256_0f,c256_0f,-c256_17,s256_17,
			c512_07,s512_07,s512_09,c512_09,c512_3f,s512_3f,-s512_2f,c512_2f,c512_23,s512_23,-s512_13,c512_13,s512_25,c512_25,-c512_35,s512_35,c512_15,s512_15,-s512_05,c512_05,s512_33,c512_33,-s512_3d,c512_3d,c512_31,s512_31,-s512_21,c512_21,s512_17,c512_17,-c512_27,s512_27
		);
	/* Block 1e: */		jt = j1 + p1e0;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+15),p_offsets, (a+jt),q_offsets,
			-c32_1,s32_1,s64_1,c64_1,-s64_3,-c64_3,c128_f,s128_f,-c128_b,-s128_b,-s128_d,c128_d,s128_9,-c128_9,c256_0f,s256_0f,-c256_07,-s256_07,-s256_0b,c256_0b,s256_03,-c256_03,s256_13,c256_13,-s256_1b,-c256_1b,-c256_17,s256_17,c256_1f,-s256_1f,
			c512_0f,s512_0f,-c512_01,s512_01,-s512_07,c512_07,-s512_09,-c512_09,s512_35,c512_35,-c512_3b,-s512_3b,-c512_3d,s512_3d,s512_33,-c512_33,c512_2d,s512_2d,-c512_1d,-s512_1d,-s512_25,c512_25,s512_15,-c512_15,s512_17,c512_17,-s512_27,-c512_27,-c512_1f,s512_1f,c512_2f,-s512_2f
		);

	/***************************** ODD-ORDER TWIDDLES ROWS: *****************************/
	j1 += p200;

	/* Block  1: */		jt = j1;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+16),p_offsets, (a+jt),q_offsets,
			c64_1,s64_1,c128_1,s128_1,c128_3,s128_3,c256_01,s256_01,c256_05,s256_05,c256_03,s256_03,c256_07,s256_07,c512_01,s512_01,c512_09,s512_09,c512_05,s512_05,c512_0d,s512_0d,c512_03,s512_03,c512_0b,s512_0b,c512_07,s512_07,c512_0f,s512_0f,
			c1024_01,s1024_01,c1024_11,s1024_11,c1024_09,s1024_09,c1024_19,s1024_19,c1024_05,s1024_05,c1024_15,s1024_15,c1024_0d,s1024_0d,c1024_1d,s1024_1d,c1024_03,s1024_03,c1024_13,s1024_13,c1024_0b,s1024_0b,c1024_1b,s1024_1b,c1024_07,s1024_07,c1024_17,s1024_17,c1024_0f,s1024_0f,c1024_1f,s1024_1f
		);
	/* Block 11: */		jt = j1 + p20;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+17),p_offsets, (a+jt),q_offsets,
			-s64_1,c64_1,s128_f,c128_f,-c128_d,s128_d,c256_11,s256_11,-s256_15,c256_15,s256_0d,c256_0d,-c256_09,s256_09,c512_11,s512_11,-s512_19,c512_19,s512_2b,c512_2b,-c512_23,s512_23,c512_33,s512_33,-s512_3b,c512_3b,s512_09,c512_09,-c512_01,s512_01,
			c1024_11,s1024_11,-s1024_21,c1024_21,s1024_67,c1024_67,-c1024_57,s1024_57,c1024_55,s1024_55,-s1024_65,c1024_65,s1024_23,c1024_23,-c1024_13,s1024_13,c1024_33,s1024_33,-s1024_43,c1024_43,s1024_45,c1024_45,-c1024_35,s1024_35,c1024_77,s1024_77,-c1024_79,s1024_79,s1024_01,c1024_01,-c1024_0f,-s1024_0f
		);
	/* Block  9: */		jt = j1 + p40;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+18),p_offsets, (a+jt),q_offsets,
			s64_7,c64_7,c128_9,s128_9,s128_5,c128_5,c256_09,s256_09,s256_13,c256_13,c256_1b,s256_1b,s256_01,c256_01,c512_09,s512_09,s512_2f,c512_2f,c512_2d,s512_2d,s512_0b,c512_0b,c512_1b,s512_1b,s512_1d,c512_1d,c512_3f,s512_3f,-s512_07,c512_07,
			c1024_09,s1024_09,s1024_67,c1024_67,c1024_51,s1024_51,s1024_1f,c1024_1f,c1024_2d,s1024_2d,s1024_43,c1024_43,c1024_75,s1024_75,-s1024_05,c1024_05,c1024_1b,s1024_1b,s1024_55,c1024_55,c1024_63,s1024_63,s1024_0d,c1024_0d,c1024_3f,s1024_3f,s1024_31,c1024_31,s1024_79,c1024_79,-s1024_17,c1024_17
		);
	/* Block 19: */		jt = j1 + p60;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+19),p_offsets, (a+jt),q_offsets,
			-c64_7,s64_7,s128_7,c128_7,-c128_b,-s128_b,c256_19,s256_19,-c256_03,s256_03,-s256_0b,c256_0b,-s256_11,-c256_11,c512_19,s512_19,-c512_1f,s512_1f,s512_03,c512_03,-s512_3b,-c512_3b,s512_35,c512_35,-c512_13,-s512_13,-s512_2f,c512_2f,-s512_09,-c512_09,
			c1024_19,s1024_19,-c1024_57,s1024_57,s1024_1f,c1024_1f,-c1024_71,-s1024_71,c1024_7d,s1024_7d,-c1024_0d,-s1024_0d,-s1024_45,c1024_45,-s1024_2b,-c1024_2b,c1024_4b,s1024_4b,-c1024_25,s1024_25,-s1024_13,c1024_13,-s1024_5d,-c1024_5d,s1024_51,c1024_51,-c1024_3f,-s1024_3f,-s1024_77,c1024_77,s1024_07,-c1024_07
		);
	/* Block  5: */		jt = j1 + p80;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+20),p_offsets, (a+jt),q_offsets,
			c64_5,s64_5,c128_5,s128_5,c128_f,s128_f,c256_05,s256_05,c256_19,s256_19,c256_0f,s256_0f,s256_1d,c256_1d,c512_05,s512_05,c512_2d,s512_2d,c512_19,s512_19,s512_3f,c512_3f,c512_0f,s512_0f,c512_37,s512_37,c512_23,s512_23,s512_35,c512_35,
			c1024_05,s1024_05,c1024_55,s1024_55,c1024_2d,s1024_2d,c1024_7d,s1024_7d,c1024_19,s1024_19,c1024_69,s1024_69,c1024_41,s1024_41,s1024_6f,c1024_6f,c1024_0f,s1024_0f,c1024_5f,s1024_5f,c1024_37,s1024_37,s1024_79,c1024_79,c1024_23,s1024_23,c1024_73,s1024_73,c1024_4b,s1024_4b,s1024_65,c1024_65
		);
	/* Block 15: */		jt = j1 + pa0;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+21),p_offsets, (a+jt),q_offsets,
			-s64_5,c64_5,s128_b,c128_b,-c128_1,s128_1,c256_15,s256_15,-c256_17,s256_17,s256_01,c256_01,-c256_13,-s256_13,c512_15,s512_15,-s512_3d,c512_3d,s512_17,c512_17,-c512_11,-s512_11,c512_3f,s512_3f,-c512_19,s512_19,-s512_13,c512_13,-c512_3b,-s512_3b,
			c1024_15,s1024_15,-s1024_65,c1024_65,s1024_43,c1024_43,-c1024_0d,-s1024_0d,c1024_69,s1024_69,-c1024_47,s1024_47,-s1024_11,c1024_11,-c1024_61,-s1024_61,c1024_3f,s1024_3f,-c1024_71,s1024_71,s1024_19,c1024_19,-c1024_37,-s1024_37,s1024_6d,c1024_6d,-c1024_1d,s1024_1d,-s1024_3b,c1024_3b,-s1024_75,-c1024_75
		);
	/* Block  d: */		jt = j1 + pc0;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+22),p_offsets, (a+jt),q_offsets,
			s64_3,c64_3,c128_d,s128_d,-s128_7,c128_7,c256_0d,s256_0d,-s256_01,c256_01,s256_19,c256_19,-s256_1b,c256_1b,c512_0d,s512_0d,s512_0b,c512_0b,s512_3f,c512_3f,-s512_29,c512_29,c512_27,s512_27,-s512_0f,c512_0f,s512_25,c512_25,-c512_3d,s512_3d,
			c1024_0d,s1024_0d,s1024_23,c1024_23,c1024_75,s1024_75,-s1024_45,c1024_45,c1024_41,s1024_41,-s1024_11,c1024_11,s1024_57,c1024_57,-s1024_79,c1024_79,c1024_27,s1024_27,s1024_09,c1024_09,s1024_71,c1024_71,-s1024_5f,c1024_5f,c1024_5b,s1024_5b,-s1024_2b,c1024_2b,s1024_3d,c1024_3d,-c1024_6d,s1024_6d
		);
	/* Block 1d: */		jt = j1 + pe0;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+23),p_offsets, (a+jt),q_offsets,
			-c64_3,s64_3,s128_3,c128_3,-s128_9,-c128_9,c256_1d,s256_1d,-c256_11,-s256_11,-s256_17,c256_17,s256_0b,-c256_0b,c512_1d,s512_1d,-c512_05,-s512_05,-s512_11,c512_11,-s512_07,-c512_07,s512_29,c512_29,-c512_3f,-s512_3f,-c512_35,s512_35,s512_33,-c512_33,
			c1024_1d,s1024_1d,-c1024_13,s1024_13,-s1024_05,c1024_05,-s1024_2b,-c1024_2b,s1024_6f,c1024_6f,-c1024_61,-s1024_61,-s1024_79,c1024_79,s1024_49,-c1024_49,c1024_57,s1024_57,-c1024_27,-s1024_27,-s1024_3f,c1024_3f,s1024_0f,-c1024_0f,s1024_35,c1024_35,-s1024_65,-c1024_65,-c1024_4d,s1024_4d,c1024_7d,-s1024_7d
		);
	/* Block  3: */		jt = j1 + p100;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+24),p_offsets, (a+jt),q_offsets,
			c64_3,s64_3,c128_3,s128_3,c128_9,s128_9,c256_03,s256_03,c256_0f,s256_0f,c256_09,s256_09,c256_15,s256_15,c512_03,s512_03,c512_1b,s512_1b,c512_0f,s512_0f,c512_27,s512_27,c512_09,s512_09,c512_21,s512_21,c512_15,s512_15,c512_2d,s512_2d,
			c1024_03,s1024_03,c1024_33,s1024_33,c1024_1b,s1024_1b,c1024_4b,s1024_4b,c1024_0f,s1024_0f,c1024_3f,s1024_3f,c1024_27,s1024_27,c1024_57,s1024_57,c1024_09,s1024_09,c1024_39,s1024_39,c1024_21,s1024_21,c1024_51,s1024_51,c1024_15,s1024_15,c1024_45,s1024_45,c1024_2d,s1024_2d,c1024_5d,s1024_5d
		);
	/* Block 13: */		jt = j1 + p120;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+25),p_offsets, (a+jt),q_offsets,
			-s64_3,c64_3,s128_d,c128_d,-c128_7,s128_7,c256_13,s256_13,-s256_1f,c256_1f,s256_07,c256_07,-c256_05,-s256_05,c512_13,s512_13,-s512_2b,c512_2b,s512_21,c512_21,-c512_09,s512_09,c512_39,s512_39,-c512_2f,s512_2f,-s512_05,c512_05,-c512_1d,-s512_1d,
			c1024_13,s1024_13,-s1024_43,c1024_43,s1024_55,c1024_55,-c1024_25,s1024_25,c1024_5f,s1024_5f,-c1024_71,s1024_71,s1024_09,c1024_09,-c1024_27,-s1024_27,c1024_39,s1024_39,-s1024_69,c1024_69,s1024_2f,c1024_2f,-c1024_01,-s1024_01,s1024_7b,c1024_7b,-c1024_4b,s1024_4b,-s1024_1d,c1024_1d,-c1024_4d,-s1024_4d
		);
	/* Block  b: */		jt = j1 + p140;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+26),p_offsets, (a+jt),q_offsets,
			s64_5,c64_5,c128_b,s128_b,-s128_1,c128_1,c256_0b,s256_0b,s256_09,c256_09,s256_1f,c256_1f,-s256_0d,c256_0d,c512_0b,s512_0b,s512_1d,c512_1d,c512_37,s512_37,-s512_0f,c512_0f,c512_21,s512_21,s512_07,c512_07,s512_33,c512_33,-s512_25,c512_25,
			c1024_0b,s1024_0b,s1024_45,c1024_45,c1024_63,s1024_63,-s1024_13,c1024_13,c1024_37,s1024_37,s1024_19,c1024_19,s1024_71,c1024_71,-s1024_3f,c1024_3f,c1024_21,s1024_21,s1024_2f,c1024_2f,c1024_79,s1024_79,-s1024_29,c1024_29,c1024_4d,s1024_4d,s1024_03,c1024_03,s1024_5b,c1024_5b,-s1024_55,c1024_55
		);
	/* Block 1b: */		jt = j1 + p160;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+27),p_offsets, (a+jt),q_offsets,
			-c64_5,s64_5,s128_5,c128_5,-s128_f,-c128_f,c256_1b,s256_1b,-c256_07,-s256_07,-s256_11,c256_11,-s256_03,-c256_03,c512_1b,s512_1b,-c512_0d,s512_0d,-s512_07,c512_07,-s512_21,-c512_21,s512_2f,c512_2f,-c512_29,-s512_29,-s512_3d,c512_3d,s512_15,-c512_15,
			c1024_1b,s1024_1b,-c1024_35,s1024_35,s1024_0d,c1024_0d,-s1024_5d,-c1024_5d,s1024_79,c1024_79,-c1024_37,-s1024_37,-s1024_5f,c1024_5f,s1024_0f,-c1024_0f,c1024_51,s1024_51,-c1024_01,-s1024_01,-s1024_29,c1024_29,-s1024_27,-c1024_27,s1024_43,c1024_43,-c1024_6d,-s1024_6d,-c1024_6b,s1024_6b,s1024_45,-c1024_45
		);
	/* Block  7: */		jt = j1 + p180;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+28),p_offsets, (a+jt),q_offsets,
			c64_7,s64_7,c128_7,s128_7,s128_b,c128_b,c256_07,s256_07,s256_1d,c256_1d,c256_15,s256_15,s256_0f,c256_0f,c512_07,s512_07,c512_3f,s512_3f,c512_23,s512_23,s512_25,c512_25,c512_15,s512_15,s512_33,c512_33,c512_31,s512_31,s512_17,c512_17,
			c1024_07,s1024_07,c1024_77,s1024_77,c1024_3f,s1024_3f,s1024_51,c1024_51,c1024_23,s1024_23,s1024_6d,c1024_6d,c1024_5b,s1024_5b,s1024_35,c1024_35,c1024_15,s1024_15,s1024_7b,c1024_7b,c1024_4d,s1024_4d,s1024_43,c1024_43,c1024_31,s1024_31,s1024_5f,c1024_5f,c1024_69,s1024_69,s1024_27,c1024_27
		);
	/* Block 17: */		jt = j1 + p1a0;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+29),p_offsets, (a+jt),q_offsets,
			-s64_7,c64_7,s128_9,c128_9,-c128_5,-s128_5,c256_17,s256_17,-c256_0d,s256_0d,-s256_05,c256_05,-s256_1f,-c256_1f,c512_17,s512_17,-c512_31,s512_31,s512_0d,c512_0d,-c512_2b,-s512_2b,s512_3b,c512_3b,-c512_03,s512_03,-s512_21,c512_21,-s512_27,-c512_27,
			c1024_17,s1024_17,-c1024_79,s1024_79,s1024_31,c1024_31,-c1024_3f,-s1024_3f,c1024_73,s1024_73,-c1024_1d,s1024_1d,-s1024_2b,c1024_2b,-s1024_65,-c1024_65,c1024_45,s1024_45,-c1024_4b,s1024_4b,s1024_03,c1024_03,-c1024_6d,-s1024_6d,s1024_5f,c1024_5f,-c1024_11,-s1024_11,-s1024_59,c1024_59,-s1024_37,-c1024_37
		);
	/* Block  f: */		jt = j1 + p1c0;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+30),p_offsets, (a+jt),q_offsets,
			s64_1,c64_1,c128_f,s128_f,-s128_d,c128_d,c256_0f,s256_0f,-s256_0b,c256_0b,s256_13,c256_13,-c256_17,s256_17,c512_0f,s512_0f,-s512_07,c512_07,s512_35,c512_35,-c512_3d,s512_3d,c512_2d,s512_2d,-s512_25,c512_25,s512_17,c512_17,-c512_1f,s512_1f,
			c1024_0f,s1024_0f,s1024_01,c1024_01,s1024_79,c1024_79,-s1024_77,c1024_77,c1024_4b,s1024_4b,-s1024_3b,c1024_3b,s1024_3d,c1024_3d,-c1024_4d,s1024_4d,c1024_2d,s1024_2d,-s1024_1d,c1024_1d,s1024_5b,c1024_5b,-c1024_6b,s1024_6b,c1024_69,s1024_69,-s1024_59,c1024_59,s1024_1f,c1024_1f,-c1024_2f,s1024_2f
		);
	/* Block 1f: */		jt = j1 + p1e0;
		RADIX_32_DIF_TWIDDLE_OOP(
			(double *)(t+31),p_offsets, (a+jt),q_offsets,
			-c64_1,s64_1,s128_1,c128_1,-s128_3,-c128_3,c256_1f,s256_1f,-c256_1b,-s256_1b,-s256_1d,c256_1d,s256_19,-c256_19,c512_1f,s512_1f,-c512_17,-s512_17,-s512_1b,c512_1b,s512_13,-c512_13,s512_23,c512_23,-s512_2b,-c512_2b,-c512_27,s512_27,c512_2f,-s512_2f,
			c1024_1f,s1024_1f,-c1024_0f,-s1024_0f,-s1024_17,c1024_17,s1024_07,-c1024_07,s1024_65,c1024_65,-s1024_75,-c1024_75,-c1024_6d,s1024_6d,c1024_7d,-s1024_7d,c1024_5d,s1024_5d,-c1024_4d,-s1024_4d,-s1024_55,c1024_55,s1024_45,-c1024_45,s1024_27,c1024_27,-s1024_37,-c1024_37,-c1024_2f,s1024_2f,c1024_3f,-s1024_3f
		);

	j1 -= p200;

	#else	// 64 x 16 version:

	/* Jan 2014: Having just spent a week implementing the SSE2/AVX code version of the radix-256 dft/carry routine, it is
	clearly preferable to structure larger-radix DFTs to re-use as much already-written-and-debugged code as possible.
	For n = 1024 that means factoring the radix as 64*16 rather than 32*32, that is, preferring a opening salvo of
	16 twiddleless DFTS [already coded for SIMD] followed by 64 radix-16 with-twiddles [again, already coded] instead of
	the 32x32 scheme, which needs new code for the radix-32-with-twiddles. */

	// Gather the needed data and do 16 twiddleless length-64 subtransforms, with p-offsets in br16 order: 084c2a6e195d3b7f:

		for(i = 0, jp = 0; i < 16; i++, jp += 64) {
			jt = j1 + po_br[i];	// po_br[] = p[084c2a6e195d3b7f]
			//	NOTE that RADIX_64_DIF outputs are IN-ORDER rather than BR:
			RADIX_64_DIF((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		}

	/*...and now do 64 radix-16 subtransforms, including the internal twiddles.

	Instead of again tediously working out the twiddles matrix for the 64x16 scheme, spent roughly the same amount of time
	coding up a utility function br.c::print_pow2_twiddles() which automates the process, and which gives the twiddles
	for our 64 x 16 impl of radix-1024 DFT in already-row/col-bit-reversed form.
	*/
	// Block 0: has all-unity twiddles
	tptr = t;
	jt = j1;	jp = j2;
	// Twiddleless DIF bit-reverses its outputs, so a_p* terms appear in BR-order [swap index pairs 1/8,2/4,3/c,5/a,7/e.b/d]:
	RADIX_16_DIF(
		tptr->re,tptr->im,(tptr+0x200)->re,(tptr+0x200)->im,(tptr+0x100)->re,(tptr+0x100)->im,(tptr+0x300)->re,(tptr+0x300)->im,(tptr+0x080)->re,(tptr+0x080)->im,(tptr+0x280)->re,(tptr+0x280)->im,(tptr+0x180)->re,(tptr+0x180)->im,(tptr+0x380)->re,(tptr+0x380)->im,(tptr+0x040)->re,(tptr+0x040)->im,(tptr+0x240)->re,(tptr+0x240)->im,(tptr+0x140)->re,(tptr+0x140)->im,(tptr+0x340)->re,(tptr+0x340)->im,(tptr+0x0c0)->re,(tptr+0x0c0)->im,(tptr+0x2c0)->re,(tptr+0x2c0)->im,(tptr+0x1c0)->re,(tptr+0x1c0)->im,(tptr+0x3c0)->re,(tptr+0x3c0)->im,
		a[jt],a[jp],a[jt+p1],a[jp+p1],a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],a[jt+p4],a[jp+p4],a[jt+p5],a[jp+p5],a[jt+p6],a[jp+p6],a[jt+p7],a[jp+p7],a[jt+p8],a[jp+p8],a[jt+p9],a[jp+p9],a[jt+pa],a[jp+pa],a[jt+pb],a[jp+pb],a[jt+pc],a[jp+pc],a[jt+pd],a[jp+pd],a[jt+pe],a[jp+pe],a[jt+pf],a[jp+pf],
		c16,s16
	);	tptr++;

	// Remaining 63 sets of macro calls done in loop:
	for(i = 1; i < 64; i++) {
		jt = j1 + i_offsets[i]; jp = j2 + i_offsets[i];	// poffs[] = p10,p20,...,p3f0
		addr = DFT1024_TWIDDLES[i]; addi = addr+1;	// Pointer to required row of 2-D twiddles array
		RADIX_16_DIF_TWIDDLE_OOP(
			tptr->re,tptr->im,(tptr+0x200)->re,(tptr+0x200)->im,(tptr+0x100)->re,(tptr+0x100)->im,(tptr+0x300)->re,(tptr+0x300)->im,(tptr+0x080)->re,(tptr+0x080)->im,(tptr+0x280)->re,(tptr+0x280)->im,(tptr+0x180)->re,(tptr+0x180)->im,(tptr+0x380)->re,(tptr+0x380)->im,(tptr+0x040)->re,(tptr+0x040)->im,(tptr+0x240)->re,(tptr+0x240)->im,(tptr+0x140)->re,(tptr+0x140)->im,(tptr+0x340)->re,(tptr+0x340)->im,(tptr+0x0c0)->re,(tptr+0x0c0)->im,(tptr+0x2c0)->re,(tptr+0x2c0)->im,(tptr+0x1c0)->re,(tptr+0x1c0)->im,(tptr+0x3c0)->re,(tptr+0x3c0)->im,
			a[jt],a[jp],a[jt+p1],a[jp+p1],a[jt+p2],a[jp+p2],a[jt+p3],a[jp+p3],a[jt+p4],a[jp+p4],a[jt+p5],a[jp+p5],a[jt+p6],a[jp+p6],a[jt+p7],a[jp+p7],a[jt+p8],a[jp+p8],a[jt+p9],a[jp+p9],a[jt+pa],a[jp+pa],a[jt+pb],a[jp+pb],a[jt+pc],a[jp+pc],a[jt+pd],a[jp+pd],a[jt+pe],a[jp+pe],a[jt+pf],a[jp+pf],
			*(addr+0x00),*(addi+0x00), *(addr+0x02),*(addi+0x02), *(addr+0x04),*(addi+0x04), *(addr+0x06),*(addi+0x06), *(addr+0x08),*(addi+0x08), *(addr+0x0a),*(addi+0x0a), *(addr+0x0c),*(addi+0x0c), *(addr+0x0e),*(addi+0x0e), *(addr+0x10),*(addi+0x10), *(addr+0x12),*(addi+0x12), *(addr+0x14),*(addi+0x14), *(addr+0x16),*(addi+0x16), *(addr+0x18),*(addi+0x18), *(addr+0x1a),*(addi+0x1a), *(addr+0x1c),*(addi+0x1c),
			c16,s16
		);	tptr++;
	}

	#endif
	}
}

/**************/

void radix1024_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-1024 complex inverse DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix64_dif_pass for further details on storage, indexing and DFT decomposition.
*/
	int i,j,j1,j2,jt,jp;
	static int NDIVR,first_entry=TRUE,
			p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
		p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p1a,p1b,p1c,p1d,p1e,p1f,
  #if RAD1 == 32
		p20,p40,p60,p80,pa0,pc0,pe0,p100,p120,p140,p160,p180,p1a0,p1c0,p1e0,p200;
	static int i_offsets[32], o_offsets[32], p_offsets[32], q_offsets[32];
  #else
		p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,p100,p110,p120,p130,p140,p150,p160,p170,p180,p190,p1a0,p1b0,p1c0,p1d0,p1e0,p1f0,p200;
	int ju,jv;
	static int i_offsets[64], o_offsets[64];
	static int poffs[16],po_br[64];
	// We prefer pointer-based array-element access, because that allows our radix16 DFT-with-twiddles
	// to look the same in terms of array-element arglists:
	const double *addr,*addi;
	struct complex *tptr;
	#include "radix1024_twiddles.h"
  #endif
	struct complex t[RADIX];

	// New runlength?
	if(!first_entry && (n >> 10) != NDIVR)
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		ASSERT((double *)t == &(t[0].re), "Unexpected value for Tmp-array-start pointer!");
		first_entry=FALSE;
		NDIVR = n >> 10;
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
		p10 = pf + NDIVR;		pf += ( (pf >> DAT_BITS) << PAD_BITS );
		p11 = p10 + NDIVR;		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p12 = p11 + NDIVR;		p11 += ( (p11 >> DAT_BITS) << PAD_BITS );
		p13 = p12 + NDIVR;		p12 += ( (p12 >> DAT_BITS) << PAD_BITS );
		p14 = p13 + NDIVR;		p13 += ( (p13 >> DAT_BITS) << PAD_BITS );
		p15 = p14 + NDIVR;		p14 += ( (p14 >> DAT_BITS) << PAD_BITS );
		p16 = p15 + NDIVR;		p15 += ( (p15 >> DAT_BITS) << PAD_BITS );
		p17 = p16 + NDIVR;		p16 += ( (p16 >> DAT_BITS) << PAD_BITS );
		p18 = p17 + NDIVR;		p17 += ( (p17 >> DAT_BITS) << PAD_BITS );
		p19 = p18 + NDIVR;		p18 += ( (p18 >> DAT_BITS) << PAD_BITS );
		p1a = p19 + NDIVR;		p19 += ( (p19 >> DAT_BITS) << PAD_BITS );
		p1b = p1a + NDIVR;		p1a += ( (p1a >> DAT_BITS) << PAD_BITS );
		p1c = p1b + NDIVR;		p1b += ( (p1b >> DAT_BITS) << PAD_BITS );
		p1d = p1c + NDIVR;		p1c += ( (p1c >> DAT_BITS) << PAD_BITS );
		p1e = p1d + NDIVR;		p1d += ( (p1d >> DAT_BITS) << PAD_BITS );
		p1f = p1e + NDIVR;		p1e += ( (p1e >> DAT_BITS) << PAD_BITS );
								p1f += ( (p1f >> DAT_BITS) << PAD_BITS );
	#if RAD1 == 32
		NDIVR <<= 5;	p20 = NDIVR;	// NDIVR holds unpadded p20
		p40 = p20 + NDIVR;		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p60 = p40 + NDIVR;		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p80 = p60 + NDIVR;		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		pa0 = p80 + NDIVR;		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		pc0 = pa0 + NDIVR;		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pe0 = pc0 + NDIVR;		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		p100= pe0 + NDIVR;		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		p120= p100+ NDIVR;		p100+= ( (p100>> DAT_BITS) << PAD_BITS );
		p140= p120+ NDIVR;		p120+= ( (p120>> DAT_BITS) << PAD_BITS );
		p160= p140+ NDIVR;		p140+= ( (p140>> DAT_BITS) << PAD_BITS );
		p180= p160+ NDIVR;		p160+= ( (p160>> DAT_BITS) << PAD_BITS );
		p1a0= p180+ NDIVR;		p180+= ( (p180>> DAT_BITS) << PAD_BITS );
		p1c0= p1a0+ NDIVR;		p1a0+= ( (p1a0>> DAT_BITS) << PAD_BITS );
		p1e0= p1c0+ NDIVR;		p1c0+= ( (p1c0>> DAT_BITS) << PAD_BITS );
		p200= p1e0+ NDIVR;		p1e0+= ( (p1e0>> DAT_BITS) << PAD_BITS );
		NDIVR >>= 5;			p200+= ( (p200>> DAT_BITS) << PAD_BITS );
	#else
		NDIVR <<= 4;	p20 = NDIVR << 1;	// NDIVR holds unpadded p10
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
		p100= pf0 + NDIVR;		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
		p110= p100 + NDIVR;		p100+= ( (p100>> DAT_BITS) << PAD_BITS );
		p120= p110 + NDIVR;		p110+= ( (p110>> DAT_BITS) << PAD_BITS );
		p130= p120 + NDIVR;		p120+= ( (p120>> DAT_BITS) << PAD_BITS );
		p140= p130 + NDIVR;		p130+= ( (p130>> DAT_BITS) << PAD_BITS );
		p150= p140 + NDIVR;		p140+= ( (p140>> DAT_BITS) << PAD_BITS );
		p160= p150 + NDIVR;		p150+= ( (p150>> DAT_BITS) << PAD_BITS );
		p170= p160 + NDIVR;		p160+= ( (p160>> DAT_BITS) << PAD_BITS );
		p180= p170 + NDIVR;		p170+= ( (p170>> DAT_BITS) << PAD_BITS );
		p190= p180 + NDIVR;		p180+= ( (p180>> DAT_BITS) << PAD_BITS );
		p1a0= p190 + NDIVR;		p190+= ( (p190>> DAT_BITS) << PAD_BITS );
		p1b0= p1a0 + NDIVR;		p1a0+= ( (p1a0>> DAT_BITS) << PAD_BITS );
		p1c0= p1b0 + NDIVR;		p1b0+= ( (p1b0>> DAT_BITS) << PAD_BITS );
		p1d0= p1c0 + NDIVR;		p1c0+= ( (p1c0>> DAT_BITS) << PAD_BITS );
		p1e0= p1d0 + NDIVR;		p1d0+= ( (p1d0>> DAT_BITS) << PAD_BITS );
		p1f0= p1e0 + NDIVR;		p1e0+= ( (p1e0>> DAT_BITS) << PAD_BITS );
		p200= p1f0 + NDIVR;		p1f0+= ( (p1f0>> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			p200+= ( (p200>> DAT_BITS) << PAD_BITS );
	#endif

	#if RAD1 == 32
		// Set array offsets for radix-32 DFT in/outputs:
		// set 1 is w.r.to: t-array:								// set 2 w.r.to a-array:
		i_offsets[0x00] =  0;		i_offsets[0x10] = p10;			o_offsets[0x00] = 0x00<<1;	o_offsets[0x10] = 0x10<<1;
		i_offsets[0x01] = p1;		i_offsets[0x11] = p11;			o_offsets[0x01] = 0x01<<1;	o_offsets[0x11] = 0x11<<1;
		i_offsets[0x02] = p2;		i_offsets[0x12] = p12;			o_offsets[0x02] = 0x02<<1;	o_offsets[0x12] = 0x12<<1;
		i_offsets[0x03] = p3;		i_offsets[0x13] = p13;			o_offsets[0x03] = 0x03<<1;	o_offsets[0x13] = 0x13<<1;
		i_offsets[0x04] = p4;		i_offsets[0x14] = p14;			o_offsets[0x04] = 0x04<<1;	o_offsets[0x14] = 0x14<<1;
		i_offsets[0x05] = p5;		i_offsets[0x15] = p15;			o_offsets[0x05] = 0x05<<1;	o_offsets[0x15] = 0x15<<1;
		i_offsets[0x06] = p6;		i_offsets[0x16] = p16;			o_offsets[0x06] = 0x06<<1;	o_offsets[0x16] = 0x16<<1;
		i_offsets[0x07] = p7;		i_offsets[0x17] = p17;			o_offsets[0x07] = 0x07<<1;	o_offsets[0x17] = 0x17<<1;
		i_offsets[0x08] = p8;		i_offsets[0x18] = p18;			o_offsets[0x08] = 0x08<<1;	o_offsets[0x18] = 0x18<<1;
		i_offsets[0x09] = p9;		i_offsets[0x19] = p19;			o_offsets[0x09] = 0x09<<1;	o_offsets[0x19] = 0x19<<1;
		i_offsets[0x0a] = pa;		i_offsets[0x1a] = p1a;			o_offsets[0x0a] = 0x0a<<1;	o_offsets[0x1a] = 0x1a<<1;
		i_offsets[0x0b] = pb;		i_offsets[0x1b] = p1b;			o_offsets[0x0b] = 0x0b<<1;	o_offsets[0x1b] = 0x1b<<1;
		i_offsets[0x0c] = pc;		i_offsets[0x1c] = p1c;			o_offsets[0x0c] = 0x0c<<1;	o_offsets[0x1c] = 0x1c<<1;
		i_offsets[0x0d] = pd;		i_offsets[0x1d] = p1d;			o_offsets[0x0d] = 0x0d<<1;	o_offsets[0x1d] = 0x1d<<1;
		i_offsets[0x0e] = pe;		i_offsets[0x1e] = p1e;			o_offsets[0x0e] = 0x0e<<1;	o_offsets[0x1e] = 0x1e<<1;
		i_offsets[0x0f] = pf;		i_offsets[0x1f] = p1f;			o_offsets[0x0f] = 0x0f<<1;	o_offsets[0x1f] = 0x1f<<1;

		// 2nd-pass (with-twiddles) I/O index offsets reverse things:
		// set 1 is w.r.to: t-array:								// set 2 w.r.to a-array:
		p_offsets[0x00] = 0x000;	p_offsets[0x10] = 0x200;		q_offsets[0x00] = 0   ;		q_offsets[0x10] = 0   +p200;
		p_offsets[0x01] = 0x020;	p_offsets[0x11] = 0x220;		q_offsets[0x01] = p20 ;		q_offsets[0x11] = p20 +p200;
		p_offsets[0x02] = 0x040;	p_offsets[0x12] = 0x240;		q_offsets[0x02] = p40 ;		q_offsets[0x12] = p40 +p200;
		p_offsets[0x03] = 0x060;	p_offsets[0x13] = 0x260;		q_offsets[0x03] = p60 ;		q_offsets[0x13] = p60 +p200;
		p_offsets[0x04] = 0x080;	p_offsets[0x14] = 0x280;		q_offsets[0x04] = p80 ;		q_offsets[0x14] = p80 +p200;
		p_offsets[0x05] = 0x0a0;	p_offsets[0x15] = 0x2a0;		q_offsets[0x05] = pa0 ;		q_offsets[0x15] = pa0 +p200;
		p_offsets[0x06] = 0x0c0;	p_offsets[0x16] = 0x2c0;		q_offsets[0x06] = pc0 ;		q_offsets[0x16] = pc0 +p200;
		p_offsets[0x07] = 0x0e0;	p_offsets[0x17] = 0x2e0;		q_offsets[0x07] = pe0 ;		q_offsets[0x17] = pe0 +p200;
		p_offsets[0x08] = 0x100;	p_offsets[0x18] = 0x300;		q_offsets[0x08] = p100;		q_offsets[0x18] = p100+p200;
		p_offsets[0x09] = 0x120;	p_offsets[0x19] = 0x320;		q_offsets[0x09] = p120;		q_offsets[0x19] = p120+p200;
		p_offsets[0x0a] = 0x140;	p_offsets[0x1a] = 0x340;		q_offsets[0x0a] = p140;		q_offsets[0x1a] = p140+p200;
		p_offsets[0x0b] = 0x160;	p_offsets[0x1b] = 0x360;		q_offsets[0x0b] = p160;		q_offsets[0x1b] = p160+p200;
		p_offsets[0x0c] = 0x180;	p_offsets[0x1c] = 0x380;		q_offsets[0x0c] = p180;		q_offsets[0x1c] = p180+p200;
		p_offsets[0x0d] = 0x1a0;	p_offsets[0x1d] = 0x3a0;		q_offsets[0x0d] = p1a0;		q_offsets[0x1d] = p1a0+p200;
		p_offsets[0x0e] = 0x1c0;	p_offsets[0x1e] = 0x3c0;		q_offsets[0x0e] = p1c0;		q_offsets[0x1e] = p1c0+p200;
		p_offsets[0x0f] = 0x1e0;	p_offsets[0x1f] = 0x3e0;		q_offsets[0x0f] = p1e0;		q_offsets[0x1f] = p1e0+p200;
		// Need real-array offsets due to pointer-cast-to-double:
		for(j=0;j<32;j++){
			p_offsets[j] <<= 1;
		}

	#else

		// Set array offsets for radix-64 DFT in/outputs:
		// set 1 is w.r.to: a-array:								// set 2 w.r.to t-array:
		i_offsets[0x00] =  0;		i_offsets[0x20] =  0+p20;		o_offsets[0x00] = 0x00<<1;	o_offsets[0x20] = 0x20<<1;
		i_offsets[0x01] = p1;		i_offsets[0x21] = p1+p20;		o_offsets[0x01] = 0x01<<1;	o_offsets[0x21] = 0x21<<1;
		i_offsets[0x02] = p2;		i_offsets[0x22] = p2+p20;		o_offsets[0x02] = 0x02<<1;	o_offsets[0x22] = 0x22<<1;
		i_offsets[0x03] = p3;		i_offsets[0x23] = p3+p20;		o_offsets[0x03] = 0x03<<1;	o_offsets[0x23] = 0x23<<1;
		i_offsets[0x04] = p4;		i_offsets[0x24] = p4+p20;		o_offsets[0x04] = 0x04<<1;	o_offsets[0x24] = 0x24<<1;
		i_offsets[0x05] = p5;		i_offsets[0x25] = p5+p20;		o_offsets[0x05] = 0x05<<1;	o_offsets[0x25] = 0x25<<1;
		i_offsets[0x06] = p6;		i_offsets[0x26] = p6+p20;		o_offsets[0x06] = 0x06<<1;	o_offsets[0x26] = 0x26<<1;
		i_offsets[0x07] = p7;		i_offsets[0x27] = p7+p20;		o_offsets[0x07] = 0x07<<1;	o_offsets[0x27] = 0x27<<1;
		i_offsets[0x08] = p8;		i_offsets[0x28] = p8+p20;		o_offsets[0x08] = 0x08<<1;	o_offsets[0x28] = 0x28<<1;
		i_offsets[0x09] = p9;		i_offsets[0x29] = p9+p20;		o_offsets[0x09] = 0x09<<1;	o_offsets[0x29] = 0x29<<1;
		i_offsets[0x0a] = pa;		i_offsets[0x2a] = pa+p20;		o_offsets[0x0a] = 0x0a<<1;	o_offsets[0x2a] = 0x2a<<1;
		i_offsets[0x0b] = pb;		i_offsets[0x2b] = pb+p20;		o_offsets[0x0b] = 0x0b<<1;	o_offsets[0x2b] = 0x2b<<1;
		i_offsets[0x0c] = pc;		i_offsets[0x2c] = pc+p20;		o_offsets[0x0c] = 0x0c<<1;	o_offsets[0x2c] = 0x2c<<1;
		i_offsets[0x0d] = pd;		i_offsets[0x2d] = pd+p20;		o_offsets[0x0d] = 0x0d<<1;	o_offsets[0x2d] = 0x2d<<1;
		i_offsets[0x0e] = pe;		i_offsets[0x2e] = pe+p20;		o_offsets[0x0e] = 0x0e<<1;	o_offsets[0x2e] = 0x2e<<1;
		i_offsets[0x0f] = pf;		i_offsets[0x2f] = pf+p20;		o_offsets[0x0f] = 0x0f<<1;	o_offsets[0x2f] = 0x2f<<1;
		i_offsets[0x10] = p10;		i_offsets[0x30] = p10+p20;		o_offsets[0x10] = 0x10<<1;	o_offsets[0x30] = 0x30<<1;
		i_offsets[0x11] = p11;		i_offsets[0x31] = p11+p20;		o_offsets[0x11] = 0x11<<1;	o_offsets[0x31] = 0x31<<1;
		i_offsets[0x12] = p12;		i_offsets[0x32] = p12+p20;		o_offsets[0x12] = 0x12<<1;	o_offsets[0x32] = 0x32<<1;
		i_offsets[0x13] = p13;		i_offsets[0x33] = p13+p20;		o_offsets[0x13] = 0x13<<1;	o_offsets[0x33] = 0x33<<1;
		i_offsets[0x14] = p14;		i_offsets[0x34] = p14+p20;		o_offsets[0x14] = 0x14<<1;	o_offsets[0x34] = 0x34<<1;
		i_offsets[0x15] = p15;		i_offsets[0x35] = p15+p20;		o_offsets[0x15] = 0x15<<1;	o_offsets[0x35] = 0x35<<1;
		i_offsets[0x16] = p16;		i_offsets[0x36] = p16+p20;		o_offsets[0x16] = 0x16<<1;	o_offsets[0x36] = 0x36<<1;
		i_offsets[0x17] = p17;		i_offsets[0x37] = p17+p20;		o_offsets[0x17] = 0x17<<1;	o_offsets[0x37] = 0x37<<1;
		i_offsets[0x18] = p18;		i_offsets[0x38] = p18+p20;		o_offsets[0x18] = 0x18<<1;	o_offsets[0x38] = 0x38<<1;
		i_offsets[0x19] = p19;		i_offsets[0x39] = p19+p20;		o_offsets[0x19] = 0x19<<1;	o_offsets[0x39] = 0x39<<1;
		i_offsets[0x1a] = p1a;		i_offsets[0x3a] = p1a+p20;		o_offsets[0x1a] = 0x1a<<1;	o_offsets[0x3a] = 0x3a<<1;
		i_offsets[0x1b] = p1b;		i_offsets[0x3b] = p1b+p20;		o_offsets[0x1b] = 0x1b<<1;	o_offsets[0x3b] = 0x3b<<1;
		i_offsets[0x1c] = p1c;		i_offsets[0x3c] = p1c+p20;		o_offsets[0x1c] = 0x1c<<1;	o_offsets[0x3c] = 0x3c<<1;
		i_offsets[0x1d] = p1d;		i_offsets[0x3d] = p1d+p20;		o_offsets[0x1d] = 0x1d<<1;	o_offsets[0x3d] = 0x3d<<1;
		i_offsets[0x1e] = p1e;		i_offsets[0x3e] = p1e+p20;		o_offsets[0x1e] = 0x1e<<1;	o_offsets[0x3e] = 0x3e<<1;
		i_offsets[0x1f] = p1f;		i_offsets[0x3f] = p1f+p20;		o_offsets[0x1f] = 0x1f<<1;	o_offsets[0x3f] = 0x3f<<1;

		po_br[4*0x0] =   0; po_br[4*0x1] = p8; po_br[4*0x2] = p4; po_br[4*0x3] = pc; po_br[4*0x4] = p2; po_br[4*0x5] = pa; po_br[4*0x6] = p6; po_br[4*0x7] = pe; po_br[4*0x8] = p1; po_br[4*0x9] = p9; po_br[4*0xa] = p5; po_br[4*0xb] = pd; po_br[4*0xc] = p3; po_br[4*0xd] = pb; po_br[4*0xe] = p7; po_br[4*0xf] = pf;
		// Each of the foregoing 16 indices is head of a (i0,i0+p20,i0+p10,i0+p30) quartet:
		for(i = 0; i < 16; i++) {
			j = i << 2;
			po_br[j+1] = po_br[j] + p20;
			po_br[j+2] = po_br[j] + p10;
			po_br[j+3] = po_br[j] + p30;
		}
		poffs[0x0] = 0; poffs[0x1] = p40; poffs[0x2] = p80; poffs[0x3] = pc0; poffs[0x4] = p100; poffs[0x5] = p140; poffs[0x6] = p180; poffs[0x7] = p1c0; poffs[0x8] = p200; poffs[0x9] = p40 + p200; poffs[0xa] = p80 + p200; poffs[0xb] = pc0 + p200; poffs[0xc] = p100 + p200; poffs[0xd] = p140 + p200; poffs[0xe] = p180 + p200; poffs[0xf] = p1c0 + p200;
	#endif
	}

/*...The radix-1024 pass is here.	*/

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
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;

	#if RAD1 == 32	// 32 x 32 2-pass DFT version:

	// Gather the needed data and do 32 twiddleless length-32 subtransforms, with p-offsets in-order:
		/*...Block 00: */	jt = j1      ;	jp = 0;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 01: */	jt = j1 + p20;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 02: */	jt = j1 + p40;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 03: */	jt = j1 + p60;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 04: */	jt = j1 + p80;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 05: */	jt = j1 + pa0;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 06: */	jt = j1 + pc0;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 07: */	jt = j1 + pe0;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 08: */	jt = j1 +p100;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 09: */	jt = j1 +p120;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 0a: */	jt = j1 +p140;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 0b: */	jt = j1 +p160;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 0c: */	jt = j1 +p180;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 0d: */	jt = j1 +p1a0;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 0e: */	jt = j1 +p1c0;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 0f: */	jt = j1 +p1e0;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		j1 += p200;
		/*...Block 10: */	jt = j1      ;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 11: */	jt = j1 + p20;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 12: */	jt = j1 + p40;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 13: */	jt = j1 + p60;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 14: */	jt = j1 + p80;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 15: */	jt = j1 + pa0;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 16: */	jt = j1 + pc0;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 17: */	jt = j1 + pe0;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 18: */	jt = j1 +p100;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 19: */	jt = j1 +p120;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 1a: */	jt = j1 +p140;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 1b: */	jt = j1 +p160;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 1c: */	jt = j1 +p180;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 1d: */	jt = j1 +p1a0;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 1e: */	jt = j1 +p1c0;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		/*...Block 1f: */	jt = j1 +p1e0;	jp+=32;		RADIX_32_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		j1 -= p200;

	/*...and now do 32 radix-32 subtransforms, including the internal twiddle factors - we use the same positive-power
	roots as in the DIF here, just fiddle with signs within the macro to effect the conjugate-multiplies. Twiddles occur
	in the same order here as DIF, but the in-and-output-index offsets are BRed: j1 + p[0+8+4+c+2+a+6+e+1+9+5+d+3+b+7+f+].
	*/
	/* Block 00 */		jt = j1;		i = 0;
		RADIX_32_DIT            (
			(double *)(t+i),p_offsets,1, (a+jt),q_offsets,RE_IM_STRIDE
		);
	/* Block 10 */		jt += p10;		i += 0x10;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			0,1,ISRT2,ISRT2,-ISRT2,ISRT2,c16,s16,-s16,c16,s16,c16,-c16,s16,c32_1,s32_1,-s32_1,c32_1,s32_3,c32_3,-c32_3,s32_3,c32_3,s32_3,-s32_3,c32_3,s32_1,c32_1,-c32_1,s32_1,
			c64_1,s64_1,-s64_1,c64_1,s64_7,c64_7,-c64_7,s64_7,c64_5,s64_5,-s64_5,c64_5,s64_3,c64_3,-c64_3,s64_3,c64_3,s64_3,-s64_3,c64_3,s64_5,c64_5,-c64_5,s64_5,c64_7,s64_7,-s64_7,c64_7,s64_1,c64_1,-c64_1,s64_1
		);
	/* Block 08 */		jt = j1 + p8;	i = 0x8;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			ISRT2,ISRT2,c16,s16,s16,c16,c32_1,s32_1,s32_3,c32_3,c32_3,s32_3,s32_1,c32_1,c64_1,s64_1,s64_7,c64_7,c64_5,s64_5,s64_3,c64_3,c64_3,s64_3,s64_5,c64_5,c64_7,s64_7,s64_1,c64_1,
			c128_1,s128_1,s128_f,c128_f,c128_9,s128_9,s128_7,c128_7,c128_5,s128_5,s128_b,c128_b,c128_d,s128_d,s128_3,c128_3,c128_3,s128_3,s128_d,c128_d,c128_b,s128_b,s128_5,c128_5,c128_7,s128_7,s128_9,c128_9,c128_f,s128_f,s128_1,c128_1
		);
	/* Block 18 */		jt += p10;		i += 0x10;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			-ISRT2,ISRT2,s16,c16,-c16,-s16,c32_3,s32_3,-c32_1,s32_1,-s32_1,c32_1,-s32_3,-c32_3,c64_3,s64_3,-c64_5,s64_5,s64_1,c64_1,-c64_7,-s64_7,s64_7,c64_7,-c64_1,-s64_1,-s64_5,c64_5,-s64_3,-c64_3,
			c128_3,s128_3,-c128_d,s128_d,s128_5,c128_5,-c128_b,-s128_b,c128_f,s128_f,-c128_1,s128_1,-s128_7,c128_7,-s128_9,-c128_9,c128_9,s128_9,-c128_7,s128_7,-s128_1,c128_1,-s128_f,-c128_f,s128_b,c128_b,-c128_5,-s128_5,-s128_d,c128_d,-s128_3,-c128_3
		);
	/* Block 04 */		jt = j1 + p4;	i = 0x4;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			c16,s16,c32_1,s32_1,c32_3,s32_3,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7,c128_1,s128_1,c128_9,s128_9,c128_5,s128_5,c128_d,s128_d,c128_3,s128_3,c128_b,s128_b,c128_7,s128_7,c128_f,s128_f,
			c256_01,s256_01,c256_11,s256_11,c256_09,s256_09,c256_19,s256_19,c256_05,s256_05,c256_15,s256_15,c256_0d,s256_0d,c256_1d,s256_1d,c256_03,s256_03,c256_13,s256_13,c256_0b,s256_0b,c256_1b,s256_1b,c256_07,s256_07,c256_17,s256_17,c256_0f,s256_0f,c256_1f,s256_1f
		);
	/* Block 14 */		jt += p10;		i += 0x10;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			-s16,c16,s32_3,c32_3,-c32_1,s32_1,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3,c128_5,s128_5,-s128_d,c128_d,s128_7,c128_7,-c128_1,-s128_1,c128_f,s128_f,-c128_9,s128_9,-s128_3,c128_3,-c128_b,-s128_b,
			c256_05,s256_05,-s256_15,c256_15,s256_13,c256_13,-c256_03,s256_03,c256_19,s256_19,-c256_17,s256_17,-s256_01,c256_01,-c256_11,-s256_11,c256_0f,s256_0f,-s256_1f,c256_1f,s256_09,c256_09,-c256_07,-s256_07,s256_1d,c256_1d,-c256_0d,s256_0d,-s256_0b,c256_0b,-c256_1b,-s256_1b
		);
	/* Block 0c */		jt = j1 + pc;	i = 0xc;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			s16,c16,c32_3,s32_3,-s32_1,c32_1,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5,c128_3,s128_3,s128_5,c128_5,c128_f,s128_f,-s128_7,c128_7,c128_9,s128_9,-s128_1,c128_1,s128_b,c128_b,-s128_d,c128_d,
			c256_03,s256_03,s256_0d,c256_0d,c256_1b,s256_1b,-s256_0b,c256_0b,c256_0f,s256_0f,s256_01,c256_01,s256_19,c256_19,-s256_17,c256_17,c256_09,s256_09,s256_07,c256_07,s256_1f,c256_1f,-s256_11,c256_11,c256_15,s256_15,-s256_05,c256_05,s256_13,c256_13,-s256_1d,c256_1d
		);
	/* Block 1c */		jt += p10;		i += 0x10;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			-c16,s16,s32_1,c32_1,-s32_3,-c32_3,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1,c128_7,s128_7,-c128_1,s128_1,-s128_3,c128_3,-s128_5,-c128_5,s128_b,c128_b,-c128_d,-s128_d,-c128_f,s128_f,s128_9,-c128_9,
			c256_07,s256_07,-c256_09,s256_09,s256_01,c256_01,-s256_11,-c256_11,s256_1d,c256_1d,-c256_13,-s256_13,-s256_1b,c256_1b,s256_0b,-c256_0b,c256_15,s256_15,-c256_05,-s256_05,-s256_0d,c256_0d,-s256_03,-c256_03,s256_0f,c256_0f,-s256_1f,-c256_1f,-c256_17,s256_17,s256_19,-c256_19
		);
	/* Block 02 */		jt = j1 + p2;	i = 0x2;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			c32_1,s32_1,c64_1,s64_1,c64_3,s64_3,c128_1,s128_1,c128_5,s128_5,c128_3,s128_3,c128_7,s128_7,c256_01,s256_01,c256_09,s256_09,c256_05,s256_05,c256_0d,s256_0d,c256_03,s256_03,c256_0b,s256_0b,c256_07,s256_07,c256_0f,s256_0f,
			c512_01,s512_01,c512_11,s512_11,c512_09,s512_09,c512_19,s512_19,c512_05,s512_05,c512_15,s512_15,c512_0d,s512_0d,c512_1d,s512_1d,c512_03,s512_03,c512_13,s512_13,c512_0b,s512_0b,c512_1b,s512_1b,c512_07,s512_07,c512_17,s512_17,c512_0f,s512_0f,c512_1f,s512_1f
		);
	/* Block 12 */		jt += p10;		i += 0x10;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			-s32_1,c32_1,s64_7,c64_7,-c64_5,s64_5,c128_9,s128_9,-s128_d,c128_d,s128_5,c128_5,-c128_1,s128_1,c256_09,s256_09,-s256_11,c256_11,s256_13,c256_13,-c256_0b,s256_0b,c256_1b,s256_1b,-c256_1d,s256_1d,s256_01,c256_01,-c256_07,-s256_07,
			c512_09,s512_09,-s512_19,c512_19,s512_2f,c512_2f,-c512_1f,s512_1f,c512_2d,s512_2d,-s512_3d,c512_3d,s512_0b,c512_0b,-c512_05,-s512_05,c512_1b,s512_1b,-s512_2b,c512_2b,s512_1d,c512_1d,-c512_0d,s512_0d,c512_3f,s512_3f,-c512_31,s512_31,-s512_07,c512_07,-c512_17,-s512_17
		);
	/* Block 0a */		jt = j1 + pa;	i = 0xa;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			s32_3,c32_3,c64_5,s64_5,s64_1,c64_1,c128_5,s128_5,s128_7,c128_7,c128_f,s128_f,-s128_3,c128_3,c256_05,s256_05,s256_13,c256_13,c256_19,s256_19,-s256_01,c256_01,c256_0f,s256_0f,s256_09,c256_09,s256_1d,c256_1d,-s256_0b,c256_0b,
			c512_05,s512_05,s512_2b,c512_2b,c512_2d,s512_2d,s512_03,c512_03,c512_19,s512_19,s512_17,c512_17,s512_3f,c512_3f,-s512_11,c512_11,c512_0f,s512_0f,s512_21,c512_21,c512_37,s512_37,-s512_07,c512_07,c512_23,s512_23,s512_0d,c512_0d,s512_35,c512_35,-s512_1b,c512_1b
		);
	/* Block 1a */		jt += p10;		i += 0x10;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			-c32_3,s32_3,s64_3,c64_3,-c64_7,-s64_7,c128_d,s128_d,-c128_1,-s128_1,-s128_7,c128_7,-s128_5,-c128_5,c256_0d,s256_0d,-c256_0b,s256_0b,-s256_01,c256_01,-s256_17,-c256_17,s256_19,c256_19,-c256_0f,-s256_0f,-s256_1b,c256_1b,s256_03,-c256_03,
			c512_0d,s512_0d,-c512_23,s512_23,s512_0b,c512_0b,-s512_3b,-c512_3b,s512_3f,c512_3f,-c512_11,-s512_11,-s512_29,c512_29,-s512_07,-c512_07,c512_27,s512_27,-c512_09,s512_09,-s512_0f,c512_0f,-s512_21,-c512_21,s512_25,c512_25,-c512_2b,-s512_2b,-c512_3d,s512_3d,s512_13,-c512_13
		);
	/* Block 06 */		jt = j1 + p6;	i = 0x6;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			c32_3,s32_3,c64_3,s64_3,s64_7,c64_7,c128_3,s128_3,c128_f,s128_f,c128_9,s128_9,s128_b,c128_b,c256_03,s256_03,c256_1b,s256_1b,c256_0f,s256_0f,s256_19,c256_19,c256_09,s256_09,s256_1f,c256_1f,c256_15,s256_15,s256_13,c256_13,
			c512_03,s512_03,c512_33,s512_33,c512_1b,s512_1b,s512_35,c512_35,c512_0f,s512_0f,c512_3f,s512_3f,c512_27,s512_27,s512_29,c512_29,c512_09,s512_09,c512_39,s512_39,c512_21,s512_21,s512_2f,c512_2f,c512_15,s512_15,s512_3b,c512_3b,c512_2d,s512_2d,s512_23,c512_23
		);
	/* Block 16 */		jt += p10;		i += 0x10;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			-s32_3,c32_3,s64_5,c64_5,-c64_1,-s64_1,c128_b,s128_b,-c128_9,s128_9,-s128_1,c128_1,-c128_d,-s128_d,c256_0b,s256_0b,-c256_1d,s256_1d,s256_09,c256_09,-c256_0f,-s256_0f,s256_1f,c256_1f,-c256_07,s256_07,-s256_0d,c256_0d,-s256_1b,-c256_1b,
			c512_0b,s512_0b,-s512_3b,c512_3b,s512_1d,c512_1d,-c512_13,-s512_13,c512_37,s512_37,-c512_19,s512_19,-s512_0f,c512_0f,-c512_3f,-s512_3f,c512_21,s512_21,-c512_2f,s512_2f,s512_07,c512_07,-c512_29,-s512_29,s512_33,c512_33,-c512_03,s512_03,-s512_25,c512_25,-s512_2b,-c512_2b
		);
	/* Block 0e */		jt = j1 + pe;	i = 0xe;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			s32_1,c32_1,c64_7,s64_7,-s64_5,c64_5,c128_7,s128_7,-s128_3,c128_3,s128_b,c128_b,-c128_f,s128_f,c256_07,s256_07,s256_01,c256_01,s256_1d,c256_1d,-s256_1b,c256_1b,c256_15,s256_15,-s256_0d,c256_0d,s256_0f,c256_0f,-c256_17,s256_17,
			c512_07,s512_07,s512_09,c512_09,c512_3f,s512_3f,-s512_2f,c512_2f,c512_23,s512_23,-s512_13,c512_13,s512_25,c512_25,-c512_35,s512_35,c512_15,s512_15,-s512_05,c512_05,s512_33,c512_33,-s512_3d,c512_3d,c512_31,s512_31,-s512_21,c512_21,s512_17,c512_17,-c512_27,s512_27
		);
	/* Block 1e */		jt += p10;		i += 0x10;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			-c32_1,s32_1,s64_1,c64_1,-s64_3,-c64_3,c128_f,s128_f,-c128_b,-s128_b,-s128_d,c128_d,s128_9,-c128_9,c256_0f,s256_0f,-c256_07,-s256_07,-s256_0b,c256_0b,s256_03,-c256_03,s256_13,c256_13,-s256_1b,-c256_1b,-c256_17,s256_17,c256_1f,-s256_1f,
			c512_0f,s512_0f,-c512_01,s512_01,-s512_07,c512_07,-s512_09,-c512_09,s512_35,c512_35,-c512_3b,-s512_3b,-c512_3d,s512_3d,s512_33,-c512_33,c512_2d,s512_2d,-c512_1d,-s512_1d,-s512_25,c512_25,s512_15,-c512_15,s512_17,c512_17,-s512_27,-c512_27,-c512_1f,s512_1f,c512_2f,-s512_2f
		);
	/* Block 01 */		jt = j1 + p1;	i = 0x1;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			c64_1,s64_1,c128_1,s128_1,c128_3,s128_3,c256_01,s256_01,c256_05,s256_05,c256_03,s256_03,c256_07,s256_07,c512_01,s512_01,c512_09,s512_09,c512_05,s512_05,c512_0d,s512_0d,c512_03,s512_03,c512_0b,s512_0b,c512_07,s512_07,c512_0f,s512_0f,
			c1024_01,s1024_01,c1024_11,s1024_11,c1024_09,s1024_09,c1024_19,s1024_19,c1024_05,s1024_05,c1024_15,s1024_15,c1024_0d,s1024_0d,c1024_1d,s1024_1d,c1024_03,s1024_03,c1024_13,s1024_13,c1024_0b,s1024_0b,c1024_1b,s1024_1b,c1024_07,s1024_07,c1024_17,s1024_17,c1024_0f,s1024_0f,c1024_1f,s1024_1f
		);
	/* Block 11 */		jt += p10;		i += 0x10;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			-s64_1,c64_1,s128_f,c128_f,-c128_d,s128_d,c256_11,s256_11,-s256_15,c256_15,s256_0d,c256_0d,-c256_09,s256_09,c512_11,s512_11,-s512_19,c512_19,s512_2b,c512_2b,-c512_23,s512_23,c512_33,s512_33,-s512_3b,c512_3b,s512_09,c512_09,-c512_01,s512_01,
			c1024_11,s1024_11,-s1024_21,c1024_21,s1024_67,c1024_67,-c1024_57,s1024_57,c1024_55,s1024_55,-s1024_65,c1024_65,s1024_23,c1024_23,-c1024_13,s1024_13,c1024_33,s1024_33,-s1024_43,c1024_43,s1024_45,c1024_45,-c1024_35,s1024_35,c1024_77,s1024_77,-c1024_79,s1024_79,s1024_01,c1024_01,-c1024_0f,-s1024_0f
		);
	/* Block 09 */		jt = j1 + p9;	i = 0x9;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			s64_7,c64_7,c128_9,s128_9,s128_5,c128_5,c256_09,s256_09,s256_13,c256_13,c256_1b,s256_1b,s256_01,c256_01,c512_09,s512_09,s512_2f,c512_2f,c512_2d,s512_2d,s512_0b,c512_0b,c512_1b,s512_1b,s512_1d,c512_1d,c512_3f,s512_3f,-s512_07,c512_07,
			c1024_09,s1024_09,s1024_67,c1024_67,c1024_51,s1024_51,s1024_1f,c1024_1f,c1024_2d,s1024_2d,s1024_43,c1024_43,c1024_75,s1024_75,-s1024_05,c1024_05,c1024_1b,s1024_1b,s1024_55,c1024_55,c1024_63,s1024_63,s1024_0d,c1024_0d,c1024_3f,s1024_3f,s1024_31,c1024_31,s1024_79,c1024_79,-s1024_17,c1024_17
		);
	/* Block 19 */		jt += p10;		i += 0x10;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			-c64_7,s64_7,s128_7,c128_7,-c128_b,-s128_b,c256_19,s256_19,-c256_03,s256_03,-s256_0b,c256_0b,-s256_11,-c256_11,c512_19,s512_19,-c512_1f,s512_1f,s512_03,c512_03,-s512_3b,-c512_3b,s512_35,c512_35,-c512_13,-s512_13,-s512_2f,c512_2f,-s512_09,-c512_09,
			c1024_19,s1024_19,-c1024_57,s1024_57,s1024_1f,c1024_1f,-c1024_71,-s1024_71,c1024_7d,s1024_7d,-c1024_0d,-s1024_0d,-s1024_45,c1024_45,-s1024_2b,-c1024_2b,c1024_4b,s1024_4b,-c1024_25,s1024_25,-s1024_13,c1024_13,-s1024_5d,-c1024_5d,s1024_51,c1024_51,-c1024_3f,-s1024_3f,-s1024_77,c1024_77,s1024_07,-c1024_07
		);
	/* Block 05 */		jt = j1 + p5;	i = 0x5;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			c64_5,s64_5,c128_5,s128_5,c128_f,s128_f,c256_05,s256_05,c256_19,s256_19,c256_0f,s256_0f,s256_1d,c256_1d,c512_05,s512_05,c512_2d,s512_2d,c512_19,s512_19,s512_3f,c512_3f,c512_0f,s512_0f,c512_37,s512_37,c512_23,s512_23,s512_35,c512_35,
			c1024_05,s1024_05,c1024_55,s1024_55,c1024_2d,s1024_2d,c1024_7d,s1024_7d,c1024_19,s1024_19,c1024_69,s1024_69,c1024_41,s1024_41,s1024_6f,c1024_6f,c1024_0f,s1024_0f,c1024_5f,s1024_5f,c1024_37,s1024_37,s1024_79,c1024_79,c1024_23,s1024_23,c1024_73,s1024_73,c1024_4b,s1024_4b,s1024_65,c1024_65
		);
	/* Block 15 */		jt += p10;		i += 0x10;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			-s64_5,c64_5,s128_b,c128_b,-c128_1,s128_1,c256_15,s256_15,-c256_17,s256_17,s256_01,c256_01,-c256_13,-s256_13,c512_15,s512_15,-s512_3d,c512_3d,s512_17,c512_17,-c512_11,-s512_11,c512_3f,s512_3f,-c512_19,s512_19,-s512_13,c512_13,-c512_3b,-s512_3b,
			c1024_15,s1024_15,-s1024_65,c1024_65,s1024_43,c1024_43,-c1024_0d,-s1024_0d,c1024_69,s1024_69,-c1024_47,s1024_47,-s1024_11,c1024_11,-c1024_61,-s1024_61,c1024_3f,s1024_3f,-c1024_71,s1024_71,s1024_19,c1024_19,-c1024_37,-s1024_37,s1024_6d,c1024_6d,-c1024_1d,s1024_1d,-s1024_3b,c1024_3b,-s1024_75,-c1024_75
		);
	/* Block 0d */		jt = j1 + pd;	i = 0xd;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			s64_3,c64_3,c128_d,s128_d,-s128_7,c128_7,c256_0d,s256_0d,-s256_01,c256_01,s256_19,c256_19,-s256_1b,c256_1b,c512_0d,s512_0d,s512_0b,c512_0b,s512_3f,c512_3f,-s512_29,c512_29,c512_27,s512_27,-s512_0f,c512_0f,s512_25,c512_25,-c512_3d,s512_3d,
			c1024_0d,s1024_0d,s1024_23,c1024_23,c1024_75,s1024_75,-s1024_45,c1024_45,c1024_41,s1024_41,-s1024_11,c1024_11,s1024_57,c1024_57,-s1024_79,c1024_79,c1024_27,s1024_27,s1024_09,c1024_09,s1024_71,c1024_71,-s1024_5f,c1024_5f,c1024_5b,s1024_5b,-s1024_2b,c1024_2b,s1024_3d,c1024_3d,-c1024_6d,s1024_6d
		);
	/* Block 1d */		jt += p10;		i += 0x10;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			-c64_3,s64_3,s128_3,c128_3,-s128_9,-c128_9,c256_1d,s256_1d,-c256_11,-s256_11,-s256_17,c256_17,s256_0b,-c256_0b,c512_1d,s512_1d,-c512_05,-s512_05,-s512_11,c512_11,-s512_07,-c512_07,s512_29,c512_29,-c512_3f,-s512_3f,-c512_35,s512_35,s512_33,-c512_33,
			c1024_1d,s1024_1d,-c1024_13,s1024_13,-s1024_05,c1024_05,-s1024_2b,-c1024_2b,s1024_6f,c1024_6f,-c1024_61,-s1024_61,-s1024_79,c1024_79,s1024_49,-c1024_49,c1024_57,s1024_57,-c1024_27,-s1024_27,-s1024_3f,c1024_3f,s1024_0f,-c1024_0f,s1024_35,c1024_35,-s1024_65,-c1024_65,-c1024_4d,s1024_4d,c1024_7d,-s1024_7d
		);
	/* Block 03 */		jt = j1 + p3;	i = 0x3;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			c64_3,s64_3,c128_3,s128_3,c128_9,s128_9,c256_03,s256_03,c256_0f,s256_0f,c256_09,s256_09,c256_15,s256_15,c512_03,s512_03,c512_1b,s512_1b,c512_0f,s512_0f,c512_27,s512_27,c512_09,s512_09,c512_21,s512_21,c512_15,s512_15,c512_2d,s512_2d,
			c1024_03,s1024_03,c1024_33,s1024_33,c1024_1b,s1024_1b,c1024_4b,s1024_4b,c1024_0f,s1024_0f,c1024_3f,s1024_3f,c1024_27,s1024_27,c1024_57,s1024_57,c1024_09,s1024_09,c1024_39,s1024_39,c1024_21,s1024_21,c1024_51,s1024_51,c1024_15,s1024_15,c1024_45,s1024_45,c1024_2d,s1024_2d,c1024_5d,s1024_5d
		);
	/* Block 13 */		jt += p10;		i += 0x10;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			-s64_3,c64_3,s128_d,c128_d,-c128_7,s128_7,c256_13,s256_13,-s256_1f,c256_1f,s256_07,c256_07,-c256_05,-s256_05,c512_13,s512_13,-s512_2b,c512_2b,s512_21,c512_21,-c512_09,s512_09,c512_39,s512_39,-c512_2f,s512_2f,-s512_05,c512_05,-c512_1d,-s512_1d,
			c1024_13,s1024_13,-s1024_43,c1024_43,s1024_55,c1024_55,-c1024_25,s1024_25,c1024_5f,s1024_5f,-c1024_71,s1024_71,s1024_09,c1024_09,-c1024_27,-s1024_27,c1024_39,s1024_39,-s1024_69,c1024_69,s1024_2f,c1024_2f,-c1024_01,-s1024_01,s1024_7b,c1024_7b,-c1024_4b,s1024_4b,-s1024_1d,c1024_1d,-c1024_4d,-s1024_4d
		);
	/* Block 0b */		jt = j1 + pb;	i = 0xb;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			s64_5,c64_5,c128_b,s128_b,-s128_1,c128_1,c256_0b,s256_0b,s256_09,c256_09,s256_1f,c256_1f,-s256_0d,c256_0d,c512_0b,s512_0b,s512_1d,c512_1d,c512_37,s512_37,-s512_0f,c512_0f,c512_21,s512_21,s512_07,c512_07,s512_33,c512_33,-s512_25,c512_25,
			c1024_0b,s1024_0b,s1024_45,c1024_45,c1024_63,s1024_63,-s1024_13,c1024_13,c1024_37,s1024_37,s1024_19,c1024_19,s1024_71,c1024_71,-s1024_3f,c1024_3f,c1024_21,s1024_21,s1024_2f,c1024_2f,c1024_79,s1024_79,-s1024_29,c1024_29,c1024_4d,s1024_4d,s1024_03,c1024_03,s1024_5b,c1024_5b,-s1024_55,c1024_55
		);
	/* Block 1b */		jt += p10;		i += 0x10;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			-c64_5,s64_5,s128_5,c128_5,-s128_f,-c128_f,c256_1b,s256_1b,-c256_07,-s256_07,-s256_11,c256_11,-s256_03,-c256_03,c512_1b,s512_1b,-c512_0d,s512_0d,-s512_07,c512_07,-s512_21,-c512_21,s512_2f,c512_2f,-c512_29,-s512_29,-s512_3d,c512_3d,s512_15,-c512_15,
			c1024_1b,s1024_1b,-c1024_35,s1024_35,s1024_0d,c1024_0d,-s1024_5d,-c1024_5d,s1024_79,c1024_79,-c1024_37,-s1024_37,-s1024_5f,c1024_5f,s1024_0f,-c1024_0f,c1024_51,s1024_51,-c1024_01,-s1024_01,-s1024_29,c1024_29,-s1024_27,-c1024_27,s1024_43,c1024_43,-c1024_6d,-s1024_6d,-c1024_6b,s1024_6b,s1024_45,-c1024_45
		);
	/* Block 07 */		jt = j1 + p7;	i = 0x7;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			c64_7,s64_7,c128_7,s128_7,s128_b,c128_b,c256_07,s256_07,s256_1d,c256_1d,c256_15,s256_15,s256_0f,c256_0f,c512_07,s512_07,c512_3f,s512_3f,c512_23,s512_23,s512_25,c512_25,c512_15,s512_15,s512_33,c512_33,c512_31,s512_31,s512_17,c512_17,
			c1024_07,s1024_07,c1024_77,s1024_77,c1024_3f,s1024_3f,s1024_51,c1024_51,c1024_23,s1024_23,s1024_6d,c1024_6d,c1024_5b,s1024_5b,s1024_35,c1024_35,c1024_15,s1024_15,s1024_7b,c1024_7b,c1024_4d,s1024_4d,s1024_43,c1024_43,c1024_31,s1024_31,s1024_5f,c1024_5f,c1024_69,s1024_69,s1024_27,c1024_27
		);
	/* Block 17 */		jt += p10;		i += 0x10;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			-s64_7,c64_7,s128_9,c128_9,-c128_5,-s128_5,c256_17,s256_17,-c256_0d,s256_0d,-s256_05,c256_05,-s256_1f,-c256_1f,c512_17,s512_17,-c512_31,s512_31,s512_0d,c512_0d,-c512_2b,-s512_2b,s512_3b,c512_3b,-c512_03,s512_03,-s512_21,c512_21,-s512_27,-c512_27,
			c1024_17,s1024_17,-c1024_79,s1024_79,s1024_31,c1024_31,-c1024_3f,-s1024_3f,c1024_73,s1024_73,-c1024_1d,s1024_1d,-s1024_2b,c1024_2b,-s1024_65,-c1024_65,c1024_45,s1024_45,-c1024_4b,s1024_4b,s1024_03,c1024_03,-c1024_6d,-s1024_6d,s1024_5f,c1024_5f,-c1024_11,-s1024_11,-s1024_59,c1024_59,-s1024_37,-c1024_37
		);
	/* Block 0f */		jt = j1 + pf;	i = 0xf;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			s64_1,c64_1,c128_f,s128_f,-s128_d,c128_d,c256_0f,s256_0f,-s256_0b,c256_0b,s256_13,c256_13,-c256_17,s256_17,c512_0f,s512_0f,-s512_07,c512_07,s512_35,c512_35,-c512_3d,s512_3d,c512_2d,s512_2d,-s512_25,c512_25,s512_17,c512_17,-c512_1f,s512_1f,
			c1024_0f,s1024_0f,s1024_01,c1024_01,s1024_79,c1024_79,-s1024_77,c1024_77,c1024_4b,s1024_4b,-s1024_3b,c1024_3b,s1024_3d,c1024_3d,-c1024_4d,s1024_4d,c1024_2d,s1024_2d,-s1024_1d,c1024_1d,s1024_5b,c1024_5b,-c1024_6b,s1024_6b,c1024_69,s1024_69,-s1024_59,c1024_59,s1024_1f,c1024_1f,-c1024_2f,s1024_2f
		);
	/* Block 1f */		jt += p10;		i += 0x10;
		RADIX_32_DIT_TWIDDLE_OOP(
			(double *)(t+i),p_offsets, (a+jt),q_offsets,
			-c64_1,s64_1,s128_1,c128_1,-s128_3,-c128_3,c256_1f,s256_1f,-c256_1b,-s256_1b,-s256_1d,c256_1d,s256_19,-c256_19,c512_1f,s512_1f,-c512_17,-s512_17,-s512_1b,c512_1b,s512_13,-c512_13,s512_23,c512_23,-s512_2b,-c512_2b,-c512_27,s512_27,c512_2f,-s512_2f,
			c1024_1f,s1024_1f,-c1024_0f,-s1024_0f,-s1024_17,c1024_17,s1024_07,-c1024_07,s1024_65,c1024_65,-s1024_75,-c1024_75,-c1024_6d,s1024_6d,c1024_7d,-s1024_7d,c1024_5d,s1024_5d,-c1024_4d,-s1024_4d,-s1024_55,c1024_55,s1024_45,-c1024_45,s1024_27,c1024_27,-s1024_37,-c1024_37,-c1024_2f,s1024_2f,c1024_3f,-s1024_3f
		);

	#else	// 64 x 16 version:

	// Gather the needed data and do 16 twiddleless length-64 subtransforms, with p-offsets in-order:

		for(i = 0, jp = 0; i < 16; i++, jp += 64) {
			jt = j1 + poffs[i];	// poffs[] = p40,p80,...,p3c0
			RADIX_64_DIT((a+jt),i_offsets,RE_IM_STRIDE, (double *)(t+jp),o_offsets,1);
		}

	/*...and now do 64 radix-16 subtransforms, including the internal twiddle factors - we use the same positive-power
	roots as in the DIF here, just fiddle with signs within the macro to effect the conjugate-multiplies. Twiddles occur
	in the same order here as DIF, but the in-and-output-index offsets are BRed: j1 + p[0,8,4,c,2,a,6,e,1,9,5,d,3,b,7,f],
	with each of the foregoing 16 indices being head of a (i0,i0+p20,i0+p10,i0+p30) quartet:
	*/
		// Block 0: has all-unity twiddles
		tptr = t;
		jt = j1;	jp = j2;
		ju = jt+p200;	jv = jp+p200;
		// Twiddleless DIF bit-reverses its outputs, so a_p* terms appear in BR-order [swap index pairs 1/8,2/4,3/c,5/a,7/e.b/d]:
		RADIX_16_DIT(
			tptr->re,tptr->im,(tptr+0x040)->re,(tptr+0x040)->im,(tptr+0x080)->re,(tptr+0x080)->im,(tptr+0x0c0)->re,(tptr+0x0c0)->im,(tptr+0x100)->re,(tptr+0x100)->im,(tptr+0x140)->re,(tptr+0x140)->im,(tptr+0x180)->re,(tptr+0x180)->im,(tptr+0x1c0)->re,(tptr+0x1c0)->im,(tptr+0x200)->re,(tptr+0x200)->im,(tptr+0x240)->re,(tptr+0x240)->im,(tptr+0x280)->re,(tptr+0x280)->im,(tptr+0x2c0)->re,(tptr+0x2c0)->im,(tptr+0x300)->re,(tptr+0x300)->im,(tptr+0x340)->re,(tptr+0x340)->im,(tptr+0x380)->re,(tptr+0x380)->im,(tptr+0x3c0)->re,(tptr+0x3c0)->im,
			a[jt     ],a[jp     ],a[jt+ p40],a[jp+ p40],a[jt+ p80],a[jp+ p80],a[jt+ pc0],a[jp+ pc0],a[jt+p100],a[jp+p100],a[jt+p140],a[jp+p140],a[jt+p180],a[jp+p180],a[jt+p1c0],a[jp+p1c0],a[ju     ],a[jv     ],a[ju+ p40],a[jv+ p40],a[ju+ p80],a[jv+ p80],a[ju+ pc0],a[jv+ pc0],a[ju+p100],a[jv+p100],a[ju+p140],a[jv+p140],a[ju+p180],a[jv+p180],a[ju+p1c0],a[jv+p1c0],
			c16,s16
		);

		// Remaining 63 sets of macro calls done in loop:
		for(i = 1; i < 64; i++) {
			tptr = t + reverse(i,6);
			jt = j1 + po_br[i]; jp = j2 + po_br[i];	// po_br[] = p[084c2a6e195d3b7f]
			ju = jt+p200;	jv = jp+p200;
			addr = DFT1024_TWIDDLES[i]; addi = addr+1;	// Pointer to required row of 2-D twiddles array
			RADIX_16_DIT_TWIDDLE_OOP(
				tptr->re,tptr->im,(tptr+0x040)->re,(tptr+0x040)->im,(tptr+0x080)->re,(tptr+0x080)->im,(tptr+0x0c0)->re,(tptr+0x0c0)->im,(tptr+0x100)->re,(tptr+0x100)->im,(tptr+0x140)->re,(tptr+0x140)->im,(tptr+0x180)->re,(tptr+0x180)->im,(tptr+0x1c0)->re,(tptr+0x1c0)->im,(tptr+0x200)->re,(tptr+0x200)->im,(tptr+0x240)->re,(tptr+0x240)->im,(tptr+0x280)->re,(tptr+0x280)->im,(tptr+0x2c0)->re,(tptr+0x2c0)->im,(tptr+0x300)->re,(tptr+0x300)->im,(tptr+0x340)->re,(tptr+0x340)->im,(tptr+0x380)->re,(tptr+0x380)->im,(tptr+0x3c0)->re,(tptr+0x3c0)->im,
				a[jt     ],a[jp     ],a[jt+ p40],a[jp+ p40],a[jt+ p80],a[jp+ p80],a[jt+ pc0],a[jp+ pc0],a[jt+p100],a[jp+p100],a[jt+p140],a[jp+p140],a[jt+p180],a[jp+p180],a[jt+p1c0],a[jp+p1c0],a[ju     ],a[jv     ],a[ju+ p40],a[jv+ p40],a[ju+ p80],a[jv+ p80],a[ju+ pc0],a[jv+ pc0],a[ju+p100],a[jv+p100],a[ju+p140],a[jv+p140],a[ju+p180],a[jv+p180],a[ju+p1c0],a[jv+p1c0],
				*(addr+0x00),*(addi+0x00), *(addr+0x02),*(addi+0x02), *(addr+0x04),*(addi+0x04), *(addr+0x06),*(addi+0x06), *(addr+0x08),*(addi+0x08), *(addr+0x0a),*(addi+0x0a), *(addr+0x0c),*(addi+0x0c), *(addr+0x0e),*(addi+0x0e), *(addr+0x10),*(addi+0x10), *(addr+0x12),*(addi+0x12), *(addr+0x14),*(addi+0x14), *(addr+0x16),*(addi+0x16), *(addr+0x18),*(addi+0x18), *(addr+0x1a),*(addi+0x1a), *(addr+0x1c),*(addi+0x1c),
				c16,s16
			);
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
	cy1024_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
#if 0
	#error This file has active debug code!
	#ifdef USE_AVX512
	const char dbg_fname[] = "dbg2.txt";
	#elif defined(USE_AVX)
	const char dbg_fname[] = "dbg1.txt";
	#endif
	FILE *o_file = 0x0;
#endif
	const char func[] = "radix1024_ditN_cy_dif1";
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr,*addi;
		struct complex *tptr;
		const int pfetch_dist = PFETCH_DIST;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
			p11,p12,p13,p14,p15,p16,p17,p18,p19,p1a,p1b,p1c,p1d,p1e,p1f,
			p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,
			p100,p110,p120,p130,p140,p150,p160,p170,p180,p190,p1a0,p1b0,p1c0,p1d0,p1e0,p1f0,p200;
		int poff[RADIX>>2];
	// DIF:
		int dif_i_offsets[64], dif_po_br[16];
	// DIT:
		int ju,jv;
		int dit_i_offsets[64], dit_poffs[16], dit_po_br[64];
	// Shared by both:
		int dif_o_offsets[64], o_offsets[64];

		int incr,j,j1,j2,jt,jp,k,l,ntmp;
		// incr = Carry-chain wts-multipliers recurrence length, which must divide
		// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 64|128|256 for avx512,avx,sse, respectively:
		const int incr_long = 16,incr_med = 8,incr_short = 4;
		// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
		if(USE_SHORT_CY_CHAIN == 0)
			incr = incr_long;
		else if(USE_SHORT_CY_CHAIN == 1)
			incr = incr_med;
		else
			incr = incr_short;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
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
	#endif
		double rt,it, wt_re,wt_im, wi_re,wi_im;	// Fermat-mod weights stuff, used in both scalar and AVX mode
		int k1,k2;

	#ifdef USE_SSE2

		const double crnd = 3.0*0x4000000*0x2000000;
		double *add0,*add1,*add2,*add3,*add4,*add5,*add6,*add7,*add8,*add9,*adda,*addb,*addc,*addd,*adde,*addf;
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *tmp,*tm0,*tm1,*tm2;	// utility ptrs
		int *itmp,*itm2;			// Pointer into the bjmodn array
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		vec_dbl *max_err, *sse2_rnd, *half_arr, *isrt2,*cc0,*ss0,	// Need DFT-16 roots explicitly
			*twid00,// ptr to head of 64 sets of twiddles shared by the 2nd-half DIF and DIT DFT macros
			*r00,	// Head of RADIX*vec_cmplx-sized local store #1
			*s1p00,	// Head of RADIX*vec_cmplx-sized local store #2
			*vd00,	// Head of 64 x vec_cmplx loal store for DFT-64 intermediates
			*cy_r,*cy_i;	// Need RADIX slots for sse2 carries, RADIX/2 for avx
	  #ifdef USE_AVX
		vec_dbl *base_negacyclic_root;
	  #endif

		/* These are used in conjunction with the langth-odd_radix arrays in the USE_SCALAR_CARRY #define below;
		In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
		int idx_offset, idx_incr;
		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_nm1;

	#else

		double *base, *baseinv;
		int p0123[4];
		const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int m,m2;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		int bjmodn[RADIX];	// Thread only carries a base datum here, must alloc a local array for remaining values
		double *cy_r = thread_arg->cy_r,*cy_i = thread_arg->cy_i, temp,frac;
		// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
		// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
		// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
		struct complex t[RADIX];
		#include "radix1024_twiddles.h"
		int *itmp;	// Pointer into the bjmodn array

	#endif

	// int data:
		int thr_id = thread_arg->tid;
		int iter = thread_arg->iter;
		int NDIVR = thread_arg->ndivr;
		int n = NDIVR*RADIX, nm1 = n-1;
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
		int sw  = thread_arg->sw, bw = n - sw;
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
		struct complex *rn0 = thread_arg->rn0;
		struct complex *rn1 = thread_arg->rn1;

		/*   constant index offsets for load/stores are here.	*/
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
		p10 = pf + NDIVR;		pf += ( (pf >> DAT_BITS) << PAD_BITS );
		p11 = p10 + NDIVR;		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p12 = p11 + NDIVR;		p11 += ( (p11 >> DAT_BITS) << PAD_BITS );
		p13 = p12 + NDIVR;		p12 += ( (p12 >> DAT_BITS) << PAD_BITS );
		p14 = p13 + NDIVR;		p13 += ( (p13 >> DAT_BITS) << PAD_BITS );
		p15 = p14 + NDIVR;		p14 += ( (p14 >> DAT_BITS) << PAD_BITS );
		p16 = p15 + NDIVR;		p15 += ( (p15 >> DAT_BITS) << PAD_BITS );
		p17 = p16 + NDIVR;		p16 += ( (p16 >> DAT_BITS) << PAD_BITS );
		p18 = p17 + NDIVR;		p17 += ( (p17 >> DAT_BITS) << PAD_BITS );
		p19 = p18 + NDIVR;		p18 += ( (p18 >> DAT_BITS) << PAD_BITS );
		p1a = p19 + NDIVR;		p19 += ( (p19 >> DAT_BITS) << PAD_BITS );
		p1b = p1a + NDIVR;		p1a += ( (p1a >> DAT_BITS) << PAD_BITS );
		p1c = p1b + NDIVR;		p1b += ( (p1b >> DAT_BITS) << PAD_BITS );
		p1d = p1c + NDIVR;		p1c += ( (p1c >> DAT_BITS) << PAD_BITS );
		p1e = p1d + NDIVR;		p1d += ( (p1d >> DAT_BITS) << PAD_BITS );
		p1f = p1e + NDIVR;		p1e += ( (p1e >> DAT_BITS) << PAD_BITS );
								p1f += ( (p1f >> DAT_BITS) << PAD_BITS );
		NDIVR <<= 4;	p20 = NDIVR << 1;	// NDIVR holds unpadded p10
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
		p100= pf0 + NDIVR;		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
		p110= p100 + NDIVR;		p100+= ( (p100>> DAT_BITS) << PAD_BITS );
		p120= p110 + NDIVR;		p110+= ( (p110>> DAT_BITS) << PAD_BITS );
		p130= p120 + NDIVR;		p120+= ( (p120>> DAT_BITS) << PAD_BITS );
		p140= p130 + NDIVR;		p130+= ( (p130>> DAT_BITS) << PAD_BITS );
		p150= p140 + NDIVR;		p140+= ( (p140>> DAT_BITS) << PAD_BITS );
		p160= p150 + NDIVR;		p150+= ( (p150>> DAT_BITS) << PAD_BITS );
		p170= p160 + NDIVR;		p160+= ( (p160>> DAT_BITS) << PAD_BITS );
		p180= p170 + NDIVR;		p170+= ( (p170>> DAT_BITS) << PAD_BITS );
		p190= p180 + NDIVR;		p180+= ( (p180>> DAT_BITS) << PAD_BITS );
		p1a0= p190 + NDIVR;		p190+= ( (p190>> DAT_BITS) << PAD_BITS );
		p1b0= p1a0 + NDIVR;		p1a0+= ( (p1a0>> DAT_BITS) << PAD_BITS );
		p1c0= p1b0 + NDIVR;		p1b0+= ( (p1b0>> DAT_BITS) << PAD_BITS );
		p1d0= p1c0 + NDIVR;		p1c0+= ( (p1c0>> DAT_BITS) << PAD_BITS );
		p1e0= p1d0 + NDIVR;		p1d0+= ( (p1d0>> DAT_BITS) << PAD_BITS );
		p1f0= p1e0 + NDIVR;		p1e0+= ( (p1e0>> DAT_BITS) << PAD_BITS );
		p200= p1f0 + NDIVR;		p1f0+= ( (p1f0>> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			p200+= ( (p200>> DAT_BITS) << PAD_BITS );
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
			poff[192+l] = poff[l] + p100+p200;
		}

	// DIF:
		// Set array offsets for radix-64 DFT in/outputs:
		// set 1 is w.r.to: a-array:										// set 2 w.r.to t-array - These are shared by DIT:
		dif_i_offsets[0x00] = 0   ;		dif_i_offsets[0x20] = 0   +p200;	o_offsets[0x00] = 0x00<<1;	o_offsets[0x20] = 0x20<<1;
		dif_i_offsets[0x01] = p10 ;		dif_i_offsets[0x21] = p10 +p200;	o_offsets[0x01] = 0x01<<1;	o_offsets[0x21] = 0x21<<1;
		dif_i_offsets[0x02] = p20 ;		dif_i_offsets[0x22] = p20 +p200;	o_offsets[0x02] = 0x02<<1;	o_offsets[0x22] = 0x22<<1;
		dif_i_offsets[0x03] = p30 ;		dif_i_offsets[0x23] = p30 +p200;	o_offsets[0x03] = 0x03<<1;	o_offsets[0x23] = 0x23<<1;
		dif_i_offsets[0x04] = p40 ;		dif_i_offsets[0x24] = p40 +p200;	o_offsets[0x04] = 0x04<<1;	o_offsets[0x24] = 0x24<<1;
		dif_i_offsets[0x05] = p50 ;		dif_i_offsets[0x25] = p50 +p200;	o_offsets[0x05] = 0x05<<1;	o_offsets[0x25] = 0x25<<1;
		dif_i_offsets[0x06] = p60 ;		dif_i_offsets[0x26] = p60 +p200;	o_offsets[0x06] = 0x06<<1;	o_offsets[0x26] = 0x26<<1;
		dif_i_offsets[0x07] = p70 ;		dif_i_offsets[0x27] = p70 +p200;	o_offsets[0x07] = 0x07<<1;	o_offsets[0x27] = 0x27<<1;
		dif_i_offsets[0x08] = p80 ;		dif_i_offsets[0x28] = p80 +p200;	o_offsets[0x08] = 0x08<<1;	o_offsets[0x28] = 0x28<<1;
		dif_i_offsets[0x09] = p90 ;		dif_i_offsets[0x29] = p90 +p200;	o_offsets[0x09] = 0x09<<1;	o_offsets[0x29] = 0x29<<1;
		dif_i_offsets[0x0a] = pa0 ;		dif_i_offsets[0x2a] = pa0 +p200;	o_offsets[0x0a] = 0x0a<<1;	o_offsets[0x2a] = 0x2a<<1;
		dif_i_offsets[0x0b] = pb0 ;		dif_i_offsets[0x2b] = pb0 +p200;	o_offsets[0x0b] = 0x0b<<1;	o_offsets[0x2b] = 0x2b<<1;
		dif_i_offsets[0x0c] = pc0 ;		dif_i_offsets[0x2c] = pc0 +p200;	o_offsets[0x0c] = 0x0c<<1;	o_offsets[0x2c] = 0x2c<<1;
		dif_i_offsets[0x0d] = pd0 ;		dif_i_offsets[0x2d] = pd0 +p200;	o_offsets[0x0d] = 0x0d<<1;	o_offsets[0x2d] = 0x2d<<1;
		dif_i_offsets[0x0e] = pe0 ;		dif_i_offsets[0x2e] = pe0 +p200;	o_offsets[0x0e] = 0x0e<<1;	o_offsets[0x2e] = 0x2e<<1;
		dif_i_offsets[0x0f] = pf0 ;		dif_i_offsets[0x2f] = pf0 +p200;	o_offsets[0x0f] = 0x0f<<1;	o_offsets[0x2f] = 0x2f<<1;
		dif_i_offsets[0x10] = p100;		dif_i_offsets[0x30] = p100+p200;	o_offsets[0x10] = 0x10<<1;	o_offsets[0x30] = 0x30<<1;
		dif_i_offsets[0x11] = p110;		dif_i_offsets[0x31] = p110+p200;	o_offsets[0x11] = 0x11<<1;	o_offsets[0x31] = 0x31<<1;
		dif_i_offsets[0x12] = p120;		dif_i_offsets[0x32] = p120+p200;	o_offsets[0x12] = 0x12<<1;	o_offsets[0x32] = 0x32<<1;
		dif_i_offsets[0x13] = p130;		dif_i_offsets[0x33] = p130+p200;	o_offsets[0x13] = 0x13<<1;	o_offsets[0x33] = 0x33<<1;
		dif_i_offsets[0x14] = p140;		dif_i_offsets[0x34] = p140+p200;	o_offsets[0x14] = 0x14<<1;	o_offsets[0x34] = 0x34<<1;
		dif_i_offsets[0x15] = p150;		dif_i_offsets[0x35] = p150+p200;	o_offsets[0x15] = 0x15<<1;	o_offsets[0x35] = 0x35<<1;
		dif_i_offsets[0x16] = p160;		dif_i_offsets[0x36] = p160+p200;	o_offsets[0x16] = 0x16<<1;	o_offsets[0x36] = 0x36<<1;
		dif_i_offsets[0x17] = p170;		dif_i_offsets[0x37] = p170+p200;	o_offsets[0x17] = 0x17<<1;	o_offsets[0x37] = 0x37<<1;
		dif_i_offsets[0x18] = p180;		dif_i_offsets[0x38] = p180+p200;	o_offsets[0x18] = 0x18<<1;	o_offsets[0x38] = 0x38<<1;
		dif_i_offsets[0x19] = p190;		dif_i_offsets[0x39] = p190+p200;	o_offsets[0x19] = 0x19<<1;	o_offsets[0x39] = 0x39<<1;
		dif_i_offsets[0x1a] = p1a0;		dif_i_offsets[0x3a] = p1a0+p200;	o_offsets[0x1a] = 0x1a<<1;	o_offsets[0x3a] = 0x3a<<1;
		dif_i_offsets[0x1b] = p1b0;		dif_i_offsets[0x3b] = p1b0+p200;	o_offsets[0x1b] = 0x1b<<1;	o_offsets[0x3b] = 0x3b<<1;
		dif_i_offsets[0x1c] = p1c0;		dif_i_offsets[0x3c] = p1c0+p200;	o_offsets[0x1c] = 0x1c<<1;	o_offsets[0x3c] = 0x3c<<1;
		dif_i_offsets[0x1d] = p1d0;		dif_i_offsets[0x3d] = p1d0+p200;	o_offsets[0x1d] = 0x1d<<1;	o_offsets[0x3d] = 0x3d<<1;
		dif_i_offsets[0x1e] = p1e0;		dif_i_offsets[0x3e] = p1e0+p200;	o_offsets[0x1e] = 0x1e<<1;	o_offsets[0x3e] = 0x3e<<1;
		dif_i_offsets[0x1f] = p1f0;		dif_i_offsets[0x3f] = p1f0+p200;	o_offsets[0x1f] = 0x1f<<1;	o_offsets[0x3f] = 0x3f<<1;
	#ifdef USE_SSE2
		// o_offs 2x,4x,... what might be expected (= 2 vec_dbl per output) due to cast-to-double-pointer of DIT-64 output arg
		for(l = 0; l < 64; l++) {
			dif_o_offsets[l] = o_offsets[l] << (L2_SZ_VD - 3);	// 2x for sse2, 4x for avx, etc
		}
	#endif

		dif_po_br[0x0] =  0; dif_po_br[0x1] = p8; dif_po_br[0x2] = p4; dif_po_br[0x3] = pc;
		dif_po_br[0x4] = p2; dif_po_br[0x5] = pa; dif_po_br[0x6] = p6; dif_po_br[0x7] = pe;
		dif_po_br[0x8] = p1; dif_po_br[0x9] = p9; dif_po_br[0xa] = p5; dif_po_br[0xb] = pd;
		dif_po_br[0xc] = p3; dif_po_br[0xd] = pb; dif_po_br[0xe] = p7; dif_po_br[0xf] = pf;

	// DIT:
		// Set array offsets for radix-64 DFT in/outputs:
		// set 1 is w.r.to: a-array:
		dit_i_offsets[0x00] =  0;		dit_i_offsets[0x20] =  0+p20;
		dit_i_offsets[0x01] = p1;		dit_i_offsets[0x21] = p1+p20;
		dit_i_offsets[0x02] = p2;		dit_i_offsets[0x22] = p2+p20;
		dit_i_offsets[0x03] = p3;		dit_i_offsets[0x23] = p3+p20;
		dit_i_offsets[0x04] = p4;		dit_i_offsets[0x24] = p4+p20;
		dit_i_offsets[0x05] = p5;		dit_i_offsets[0x25] = p5+p20;
		dit_i_offsets[0x06] = p6;		dit_i_offsets[0x26] = p6+p20;
		dit_i_offsets[0x07] = p7;		dit_i_offsets[0x27] = p7+p20;
		dit_i_offsets[0x08] = p8;		dit_i_offsets[0x28] = p8+p20;
		dit_i_offsets[0x09] = p9;		dit_i_offsets[0x29] = p9+p20;
		dit_i_offsets[0x0a] = pa;		dit_i_offsets[0x2a] = pa+p20;
		dit_i_offsets[0x0b] = pb;		dit_i_offsets[0x2b] = pb+p20;
		dit_i_offsets[0x0c] = pc;		dit_i_offsets[0x2c] = pc+p20;
		dit_i_offsets[0x0d] = pd;		dit_i_offsets[0x2d] = pd+p20;
		dit_i_offsets[0x0e] = pe;		dit_i_offsets[0x2e] = pe+p20;
		dit_i_offsets[0x0f] = pf;		dit_i_offsets[0x2f] = pf+p20;
		dit_i_offsets[0x10] = p10;		dit_i_offsets[0x30] = p10+p20;
		dit_i_offsets[0x11] = p11;		dit_i_offsets[0x31] = p11+p20;
		dit_i_offsets[0x12] = p12;		dit_i_offsets[0x32] = p12+p20;
		dit_i_offsets[0x13] = p13;		dit_i_offsets[0x33] = p13+p20;
		dit_i_offsets[0x14] = p14;		dit_i_offsets[0x34] = p14+p20;
		dit_i_offsets[0x15] = p15;		dit_i_offsets[0x35] = p15+p20;
		dit_i_offsets[0x16] = p16;		dit_i_offsets[0x36] = p16+p20;
		dit_i_offsets[0x17] = p17;		dit_i_offsets[0x37] = p17+p20;
		dit_i_offsets[0x18] = p18;		dit_i_offsets[0x38] = p18+p20;
		dit_i_offsets[0x19] = p19;		dit_i_offsets[0x39] = p19+p20;
		dit_i_offsets[0x1a] = p1a;		dit_i_offsets[0x3a] = p1a+p20;
		dit_i_offsets[0x1b] = p1b;		dit_i_offsets[0x3b] = p1b+p20;
		dit_i_offsets[0x1c] = p1c;		dit_i_offsets[0x3c] = p1c+p20;
		dit_i_offsets[0x1d] = p1d;		dit_i_offsets[0x3d] = p1d+p20;
		dit_i_offsets[0x1e] = p1e;		dit_i_offsets[0x3e] = p1e+p20;
		dit_i_offsets[0x1f] = p1f;		dit_i_offsets[0x3f] = p1f+p20;

		dit_po_br[4*0x0] =  0; dit_po_br[4*0x1] = p8; dit_po_br[4*0x2] = p4; dit_po_br[4*0x3] = pc;
		dit_po_br[4*0x4] = p2; dit_po_br[4*0x5] = pa; dit_po_br[4*0x6] = p6; dit_po_br[4*0x7] = pe;
		dit_po_br[4*0x8] = p1; dit_po_br[4*0x9] = p9; dit_po_br[4*0xa] = p5; dit_po_br[4*0xb] = pd;
		dit_po_br[4*0xc] = p3; dit_po_br[4*0xd] = pb; dit_po_br[4*0xe] = p7; dit_po_br[4*0xf] = pf;
		// Each of the foregoing 16 indices is head of a (i0,i0+p20,i0+p10,i0+p30) quartet:
		for(l = 0; l < 16; l++) {
			j = l << 2;
			dit_po_br[j+1] = dit_po_br[j] + p20;
			dit_po_br[j+2] = dit_po_br[j] + p10;
			dit_po_br[j+3] = dit_po_br[j] + p30;
		}
		dit_poffs[0x0] =           0; dit_poffs[0x1] =         p40; dit_poffs[0x2] =         p80; dit_poffs[0x3] =         pc0;
		dit_poffs[0x4] =        p100; dit_poffs[0x5] =        p140; dit_poffs[0x6] =        p180; dit_poffs[0x7] =        p1c0;
		dit_poffs[0x8] =        p200; dit_poffs[0x9] =  p40 + p200; dit_poffs[0xa] =  p80 + p200; dit_poffs[0xb] =  pc0 + p200;
		dit_poffs[0xc] = p100 + p200; dit_poffs[0xd] = p140 + p200; dit_poffs[0xe] = p180 + p200; dit_poffs[0xf] = p1c0 + p200;

	#ifdef USE_SSE2
		tmp = r00 = thread_arg->r00;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x800;	s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
		tmp += 0x800;	vd00  = tmp;	// Head of 64 x vec_cmplx = 128 vec_dbl local store for DFT-64 intermediates
		tmp += 0x80;	// sc_ptr += 0x1080
		// DFT-16 roots needed explicitly:
		isrt2  = tmp + 0x00;
		cc0    = tmp + 0x01;
		ss0    = tmp + 0x02;
		tmp += 0x3;
		// ptrs to 64 sets (2*15 = 30 vec_dbl data each) of non-unity twiddles shared by the 2nd-half DIF and DIT DFT macros:
		twid00  = tmp;
		tmp += 0x780;	// sc_ptr += 0x1083 + 4*0x1e0 = 0x1803
	  #ifdef USE_AVX512
		cy_r = tmp;	cy_i = tmp+0x080;	tmp += 0x100;	// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod
	  #elif defined(USE_AVX)
		cy_r = tmp;	cy_i = tmp+0x100;	tmp += 0x200;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x1803 + 0x202 = 0x1a05... This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 68 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod
	  #else
		cy_r = tmp;	cy_i = tmp+0x200;	tmp += 0x400;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// 0x6c5 +2 = 0x6c7 = 1735 complex
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, 2 for Fermat-mod */
	  #endif

		ASSERT((r00 == thread_arg->r00), "thread-local memcheck failed!");
		ASSERT((half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
		ASSERT((isrt2->d0 == ISRT2 && isrt2->d1 == ISRT2), "thread-local memcheck failed!");
	  #ifndef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts:
		ASSERT((sse2_rnd->d0 == crnd && sse2_rnd->d1 == crnd), "thread-local memcheck failed!");
	  #endif

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
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
		} else {
		#ifdef USE_AVX512
			/* No-Op */
		#else
			dtmp = (half_arr)->d0 * (half_arr+1)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (half_arr)->d1 * (half_arr+1)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#endif
		}

		VEC_DBL_INIT(max_err, 0.0);

		sign_mask = (uint64*)(r00 + radix1024_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (#doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		sse_sw  = sse_bw    + RE_IM_STRIDE;
		sse_nm1 = sse_sw    + RE_IM_STRIDE;

	  #ifdef USE_AVX512
	   #ifdef CARRY_16_WAY
		n_minus_sil   = (struct uint32x16*)sse_nm1 + 1;
		n_minus_silp1 = (struct uint32x16*)sse_nm1 + 2;
		sinwt         = (struct uint32x16*)sse_nm1 + 3;
		sinwtm1       = (struct uint32x16*)sse_nm1 + 4;
	   #else
		n_minus_sil   = (struct uint32x8 *)sse_nm1 + 1;
		n_minus_silp1 = (struct uint32x8 *)sse_nm1 + 2;
		sinwt         = (struct uint32x8 *)sse_nm1 + 3;
		sinwtm1       = (struct uint32x8 *)sse_nm1 + 4;
	   #endif
	  #elif defined(USE_AVX)
		n_minus_sil   = (struct uint32x4 *)sse_nm1 + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_nm1 + 2;
		sinwt         = (struct uint32x4 *)sse_nm1 + 3;
		sinwtm1       = (struct uint32x4 *)sse_nm1 + 4;
	  #endif
	  #ifdef USE_AVX
		bjmodn = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn = (int*)(sse_nm1 + RE_IM_STRIDE);
	  #endif

	#else

		// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
		base    = (double *)thread_arg->r00     ;
		baseinv = (double *)thread_arg->half_arr;

	#endif	// USE_SSE2 ?

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Init DWT-indices: */
			uint32 bjmodnini = thread_arg->bjmodnini;
			bjmodn[0] = thread_arg->bjmodn0;
			for(l = 1; l < RADIX; l++) {	// must use e.g. l for loop idx here as i is used for dwt indexing
				MOD_ADD32(bjmodn[l-1], bjmodnini, n, bjmodn[l]);
			}

			/* init carries	*/
			addr = thread_arg->cy_r;
		#ifdef USE_AVX512
			tmp = cy_r;
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
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 4, ++tmp) {
				tmp->d0 = *(addr+l  );
				tmp->d1 = *(addr+l+1);
				tmp->d2 = *(addr+l+2);
				tmp->d3 = *(addr+l+3);
			}
		#elif defined(USE_SSE2)
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 2, ++tmp) {
				tmp->d0 = *(addr+l  );
				tmp->d1 = *(addr+l+1);
			}
		#elif 0	// No-op in scalar case, since carry pattern matches that of thread data
			for(l = 0; l < RADIX; l++) {
				cy_r[l] = *(addr+l);
			}
		#endif
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
			/* init carries	*/
			addr = thread_arg->cy_r;	addi = thread_arg->cy_i;
		#ifdef USE_AVX512
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 8, ++tmp, ++tm2) {
				tmp->d0 = *(addr+l  );		tm2->d0 = *(addi+l  );
				tmp->d1 = *(addr+l+1);		tm2->d1 = *(addi+l+1);
				tmp->d2 = *(addr+l+2);		tm2->d2 = *(addi+l+2);
				tmp->d3 = *(addr+l+3);		tm2->d3 = *(addi+l+3);
				tmp->d4 = *(addr+l+4);		tm2->d4 = *(addi+l+4);
				tmp->d5 = *(addr+l+5);		tm2->d5 = *(addi+l+5);
				tmp->d6 = *(addr+l+6);		tm2->d6 = *(addi+l+6);
				tmp->d7 = *(addr+l+7);		tm2->d7 = *(addi+l+7);
			}
		#elif defined(USE_AVX)
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 4, ++tmp, ++tm2) {
				tmp->d0 = *(addr+l  );		tm2->d0 = *(addi+l  );
				tmp->d1 = *(addr+l+1);		tm2->d1 = *(addi+l+1);
				tmp->d2 = *(addr+l+2);		tm2->d2 = *(addi+l+2);
				tmp->d3 = *(addr+l+3);		tm2->d3 = *(addi+l+3);
			}
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			tmp = cy_r;
			for(l = 0; l < RADIX; l++, ++tmp) {
				// This relies on the cy_R,i sections of the SIMD data being contiguous, i.e.
				// step-thru the cy_r data via the tmp-pointer takes us seamlessly into the cy_i:
				tmp->d0 = *(addr+l  );		tmp->d1 = *(addi+l  );
			}
		#elif 0	// No-op in scalar case, since carry pattern matches that of thread data
			for(l = 0; l < RADIX; l++) {
				cy_r[l] = *(addr+l);		cy_i[l] = *(addi+l);
			}
		#endif
		}

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix1024_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			addr = thread_arg->cy_r;
		#ifdef USE_AVX512
			tmp = cy_r;
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
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 4, ++tmp) {
				*(addr+l  ) = tmp->d0;
				*(addr+l+1) = tmp->d1;
				*(addr+l+2) = tmp->d2;
				*(addr+l+3) = tmp->d3;
			}
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			tmp = cy_r;
			for(l = 0; l < RADIX; l += 2, ++tmp) {
				*(addr+l  ) = tmp->d0;
				*(addr+l+1) = tmp->d1;
			}
			maxerr = MAX(max_err->d0,max_err->d1);
		#elif 0	// No_op in scalar case, since carry pattern matches that of thread data
			for(l = 0; l < RADIX; l++) {
				*(addr+l) = cy_r[l];
			}
		#endif
		}
		else
		{
			addr = thread_arg->cy_r;	addi = thread_arg->cy_i;
		#ifdef USE_AVX512
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 8, ++tmp, ++tm2) {
				*(addr+l  ) = tmp->d0;		*(addi+l  ) = tm2->d0;
				*(addr+l+1) = tmp->d1;		*(addi+l+1) = tm2->d1;
				*(addr+l+2) = tmp->d2;		*(addi+l+2) = tm2->d2;
				*(addr+l+3) = tmp->d3;		*(addi+l+3) = tm2->d3;
				*(addr+l+4) = tmp->d4;		*(addi+l+4) = tm2->d4;
				*(addr+l+5) = tmp->d5;		*(addi+l+5) = tm2->d5;
				*(addr+l+6) = tmp->d6;		*(addi+l+6) = tm2->d6;
				*(addr+l+7) = tmp->d7;		*(addi+l+7) = tm2->d7;
			}
			t0 = MAX(max_err->d0,max_err->d1);
			t1 = MAX(max_err->d2,max_err->d3);
			t2 = MAX(max_err->d4,max_err->d5);
			t3 = MAX(max_err->d6,max_err->d7);
			maxerr = MAX( MAX(t0,t1), MAX(t2,t3) );
		#elif defined(USE_AVX)
			tmp = cy_r;	tm2 = cy_i;
			for(l = 0; l < RADIX; l += 4, ++tmp, ++tm2) {
				*(addr+l  ) = tmp->d0;		*(addi+l  ) = tm2->d0;
				*(addr+l+1) = tmp->d1;		*(addi+l+1) = tm2->d1;
				*(addr+l+2) = tmp->d2;		*(addi+l+2) = tm2->d2;
				*(addr+l+3) = tmp->d3;		*(addi+l+3) = tm2->d3;
			}
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			tmp = cy_r;
			for(l = 0; l < RADIX; l++, ++tmp) {
				// This relies on the cy_R,i sections of the SIMD data being contiguous, i.e.
				// step-thru the cy_r data via the tmp-pointer takes us seamlessly into the cy_i:
				*(addr+l  ) = tmp->d0;		*(addi+l  ) = tmp->d1;
			}
			maxerr = MAX(max_err->d0,max_err->d1);
		#elif 0	// No_op in scalar case, since carry pattern matches that of thread data
			for(l = 0; l < RADIX; l++) {
				*(addr+l) = cy_r[l];		*(addi+l) = cy_i[l];
			}
		#endif
		}

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
