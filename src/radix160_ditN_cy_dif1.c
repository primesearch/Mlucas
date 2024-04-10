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
#include "radix32.h"

#define RADIX 160	// Use #define rather than const int to ensure it's really a compile-time const in the C sense
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

// See the radix28 version of this routine for details about the
// small-cyclic-array indexing scheme used in the fermat_carry_norm_errcheckB macros.

#ifdef USE_SSE2

	#include "sse2_macro_gcc64.h"

  // For Mersenne-mod need (16 [SSE2] or 64 [AVX]) + (4 [HIACC] or 40 [LOACC]) added slots for half_arr lookup tables.
  // Max = (40 [SSE2]; 132 [AVX]), add to (half_arr_offset160[=0x292] + RADIX) to get required value of radix160_creals_in_local_store:
  #ifdef USE_AVX512	// RADIX/8 = 0x14 fewer carry slots than AVX:
	const int half_arr_offset160 = 0x2a6;	// + RADIX = 0x2a6 + 0xa0 = 0x346; Used for thread local-storage-integrity checking
	const int radix160_creals_in_local_store = 0x3cc;	// += 132 (=0x84) and round up to nearest multiple of 4
  #elif defined(USE_AVX)
	const int half_arr_offset160 = 0x2ba;	// + RADIX = 0x2ba + 0xa0 = 0x35a; Used for thread local-storage-integrity checking
	const int radix160_creals_in_local_store = 0x3e0;	// += 132 (=0x84) and round up to nearest multiple of 4
  #else
	const int half_arr_offset160 = 0x2e2;	// + RADIX = 0x2e2 + 0xa0 = 0x382; Used for thread local-storage-integrity checking
	const int radix160_creals_in_local_store = 0x3ac;	// += 40 (=0x28) and round up to nearest multiple of 4
  #endif

	#include "radix32_ditN_cy_dif1_asm.h"

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

int radix160_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-160 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-160 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix160_ditN_cy_dif1";
	const int pfetch_dist = PFETCH_DIST;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
	static double wts_mult[2], inv_mult[2];	// Const wts-multiplier and 2*(its multiplicative inverse)
	double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	// Cleanup loop assumes carryins propagate at most 4 words up, but need at least 1 vec_cmplx
	// (2 vec_dbl)'s worth of doubles in wraparound step, hence AVX-512 needs value bumped up:
  #ifdef USE_AVX512
	const int jhi_wrap = 15;
  #else
	const int jhi_wrap =  7;
  #endif
	int NDIVR,i,incr,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,ntmp,outer,nbytes;
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 10|20|40 for avx512,avx,sse, respectively:
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

	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90, nsave = 0;
	static int poff[RADIX>>2];	// Store mults of p4 offset for loop control
#ifndef MULTITHREAD
// Shared DIF+DIT:
	double rt,it;
	static int t_offsets[32];
	// Need storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 2*9 elts:
	static int dif_offsets[RADIX], dif_p20_cperms[18], dif_p20_lo_offset[32], dif_phi[ODD_RADIX];
	static int dit_offsets[RADIX], dit_p20_cperms[18], dit_p20_lo_offset[32], dit_phi[ODD_RADIX];
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
	double *addr;
	struct complex t[RADIX], *tptr;
	int *itmp,*itm2;	// Pointer into the bjmodn array
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
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

	// Uint64 bitmaps for alternate "rounded the other way" copies of sqrt2,isrt2. Default round-to-nearest versions
	// (SQRT2, ISRT2) end in ...3BCD. Since we round these down as ...3BCC90... --> ..3BCC, append _dn to varnames:
	const uint64 sqrt2_dn = 0x3FF6A09E667F3BCCull, isrt2_dn = 0x3FE6A09E667F3BCCull;
	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
	const double crnd = 3.0*0x4000000*0x2000000;
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	vec_dbl
	#ifndef MULTITHREAD
		*va0,*va1,*va2,*va3,*va4, *vb0,*vb1,*vb2,*vb3,*vb4,
		*wa0,*wa1,*wa2,*wa3,*wa4, *wb0,*wb1,*wb2,*wb3,*wb4,
	#endif
		*tmp,*tm0,*tm1,*tm2;	// Non-static utility ptrs
	static vec_dbl *two,*one,*sqrt2,*isrt2,*xcc1,*xss1,*xcc2,*xss2,*xcc3,*xss3,	// radix-32 DFT trig consts
		*ycc1,*yss1,*ycc2,*yss2,*yss3,	// radiy-5 DFT trig consts
		*max_err, *sse2_rnd, *half_arr,
		*r00,	// Head of RADIX*vec_cmplx-sized local store #1
		*s1p00,	// Head of RADIX*vec_cmplx-sized local store #2
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
	static task_control_t   task_control = {NULL, (void*)cy160_process_chunk, NULL, 0x0};

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

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=-1;
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
		// consisting of radix160_creals_in_local_store dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix160_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix160_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = sc_ptr;	r00   = tmp;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x140;	s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
		tmp += 0x140;
		two   = tmp + 0x0;	// AVX+ versions of Radix-32 DFT macros assume consts 2.0,1.0,sqrt2,isrt2 laid out thusly
		one   = tmp + 0x1;
		sqrt2 = tmp + 0x2;
		isrt2 = tmp + 0x3;
		xcc2  = tmp + 0x4;	// Radix-32 DFT macros assume roots stored in this [8th, 16th, 32nd_1,3] order
		xss2  = tmp + 0x5;
		xcc1  = tmp + 0x6;
		xss1  = tmp + 0x7;
		xcc3  = tmp + 0x8;
		xss3  = tmp + 0x9;
		ycc1  = tmp + 0xa;	// radix-5 DFT trig consts
		ycc2  = tmp + 0xb;
		yss1  = tmp + 0xc;
		yss2  = tmp + 0xd;
		yss3  = tmp + 0xe;
		tmp += 0x10;	// sc_ptr += 0x290
	  #ifdef USE_AVX512
		cy = tmp;		tmp += 0x14;	// RADIX/8 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
		cy = tmp;		tmp += 0x28;	// RADIX/4 vec_dbl slots
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(290 + 28 + 2) = 0x2ba; This is where the value of half_arr_offset160 comes from
		half_arr= tmp + 0x02;	// This table needs 20 vec_dbl in both avx and sse2 mode
	  #else
		cy = tmp;		tmp += 0x50;	// RADIX/2 vec_dbl slots
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(290 + 50 + 2) = 0x2e2; This is where the value of half_arr_offset160 comes from
		half_arr= tmp + 0x02;
	  #endif
		ASSERT((radix160_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r00) + (20 << L2_SZ_VD), "radix160_creals_in_local_store checksum failed!");
		/* These remain fixed: */
		VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );
	  #if 1
		// 2 unnamed slots for alternate "rounded the other way" copies of sqrt2,isrt2:
		dtmp = *(double *)&sqrt2_dn;	VEC_DBL_INIT(sqrt2, dtmp);
		dtmp = *(double *)&isrt2_dn;	VEC_DBL_INIT(isrt2, dtmp);
	  #else
		VEC_DBL_INIT(sqrt2, SQRT2);		VEC_DBL_INIT(isrt2, ISRT2);
	  #endif
		VEC_DBL_INIT(xcc2, c16  );	VEC_DBL_INIT(xss2, s16  );
		VEC_DBL_INIT(xcc1, c32_1);	VEC_DBL_INIT(xss1, s32_1);
		VEC_DBL_INIT(xcc3, c32_3);	VEC_DBL_INIT(xss3, s32_3);
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
		NDIVR >>= 4;			p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
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

	#ifndef MULTITHREAD
	// Shared:
		// Set array offsets for radix-32 DFT in/outputs:
		// t_offsets w.r.to: t-array, same for all 5 DFTs:
		t_offsets[0x00] = 0x00<<1;	t_offsets[0x10] = 0x10<<1;
		t_offsets[0x01] = 0x01<<1;	t_offsets[0x11] = 0x11<<1;
		t_offsets[0x02] = 0x02<<1;	t_offsets[0x12] = 0x12<<1;
		t_offsets[0x03] = 0x03<<1;	t_offsets[0x13] = 0x13<<1;
		t_offsets[0x04] = 0x04<<1;	t_offsets[0x14] = 0x14<<1;
		t_offsets[0x05] = 0x05<<1;	t_offsets[0x15] = 0x15<<1;
		t_offsets[0x06] = 0x06<<1;	t_offsets[0x16] = 0x16<<1;
		t_offsets[0x07] = 0x07<<1;	t_offsets[0x17] = 0x17<<1;
		t_offsets[0x08] = 0x08<<1;	t_offsets[0x18] = 0x18<<1;
		t_offsets[0x09] = 0x09<<1;	t_offsets[0x19] = 0x19<<1;
		t_offsets[0x0a] = 0x0a<<1;	t_offsets[0x1a] = 0x1a<<1;
		t_offsets[0x0b] = 0x0b<<1;	t_offsets[0x1b] = 0x1b<<1;
		t_offsets[0x0c] = 0x0c<<1;	t_offsets[0x1c] = 0x1c<<1;
		t_offsets[0x0d] = 0x0d<<1;	t_offsets[0x1d] = 0x1d<<1;
		t_offsets[0x0e] = 0x0e<<1;	t_offsets[0x1e] = 0x1e<<1;
		t_offsets[0x0f] = 0x0f<<1;	t_offsets[0x1f] = 0x1f<<1;

	/*** DIF indexing stuff: ***/

		dif_phi[0] =   0;		dit_phi[0] =   0;
		dif_phi[1] = p80;		dit_phi[1] = p60;
		dif_phi[2] = p60;		dit_phi[2] = p20;
		dif_phi[3] = p40;		dit_phi[3] = p80;
		dif_phi[4] = p20;		dit_phi[4] = p40;

		// Init storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 2*9

	   #ifdef USE_SSE2
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p10*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:
		l = 0;
		dif_p20_cperms[l++] = 0x00<<1; dif_p20_cperms[l++] = 0x80<<1; dif_p20_cperms[l++] = 0x60<<1; dif_p20_cperms[l++] = 0x40<<1; dif_p20_cperms[l++] = 0x20<<1; dif_p20_cperms[l++] = 0x00<<1; dif_p20_cperms[l++] = 0x80<<1; dif_p20_cperms[l++] = 0x60<<1; dif_p20_cperms[l++] = 0x40<<1;
		dif_p20_cperms[l++] = 0x90<<1; dif_p20_cperms[l++] = 0x70<<1; dif_p20_cperms[l++] = 0x50<<1; dif_p20_cperms[l++] = 0x30<<1; dif_p20_cperms[l++] = 0x10<<1; dif_p20_cperms[l++] = 0x90<<1; dif_p20_cperms[l++] = 0x70<<1; dif_p20_cperms[l++] = 0x50<<1; dif_p20_cperms[l++] = 0x30<<1;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-5 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-$ index, bits <4:30> for the p0-f. In SIMD mode again replace p[0-f] with 0x[0-f]<<1 for use with
		// contig-local-mem, thus replace 'p' prefixes with 0x and 3-left-shifts with << 4 to account for the << 1:
		l = 0;
		dif_p20_lo_offset[l++] = ((0x0 << 4) + 0);
		dif_p20_lo_offset[l++] = ((0xb << 4) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x6 << 4) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x1 << 4) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xc << 4) + 1);
		dif_p20_lo_offset[l++] = ((0x7 << 4) + 1);
		dif_p20_lo_offset[l++] = ((0x2 << 4) + 1);
		dif_p20_lo_offset[l++] = ((0xd << 4) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x8 << 4) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x3 << 4) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xe << 4) + 2);
		dif_p20_lo_offset[l++] = ((0x9 << 4) + 2);
		dif_p20_lo_offset[l++] = ((0x4 << 4) + 2);
		dif_p20_lo_offset[l++] = ((0xf << 4) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xa << 4) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x5 << 4) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x0 << 4) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xb << 4) + 3);
		dif_p20_lo_offset[l++] = ((0x6 << 4) + 3);
		dif_p20_lo_offset[l++] = ((0x1 << 4) + 3);
		dif_p20_lo_offset[l++] = ((0xc << 4) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x7 << 4) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x2 << 4) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xd << 4) + 4);
		dif_p20_lo_offset[l++] = ((0x8 << 4) + 4);
		dif_p20_lo_offset[l++] = ((0x3 << 4) + 4);
		dif_p20_lo_offset[l++] = ((0xe << 4) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x9 << 4) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x4 << 4) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xf << 4) + 0);
		dif_p20_lo_offset[l++] = ((0xa << 4) + 0);
		dif_p20_lo_offset[l++] = ((0x5 << 4) + 0);

	   #else

		// Init storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 2*9
		l = 0;
		dif_p20_cperms[l++] =   0; dif_p20_cperms[l++] = p80; dif_p20_cperms[l++] = p60; dif_p20_cperms[l++] = p40; dif_p20_cperms[l++] = p20; dif_p20_cperms[l++] =   0; dif_p20_cperms[l++] = p80; dif_p20_cperms[l++] = p60; dif_p20_cperms[l++] = p40;
		dif_p20_cperms[l++] = p90; dif_p20_cperms[l++] = p70; dif_p20_cperms[l++] = p50; dif_p20_cperms[l++] = p30; dif_p20_cperms[l++] = p10; dif_p20_cperms[l++] = p90; dif_p20_cperms[l++] = p70; dif_p20_cperms[l++] = p50; dif_p20_cperms[l++] = p30;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-5 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-4 index, bits <4:30> for the p0-f:
		l = 0;
		dif_p20_lo_offset[l++] = (( 0 << 3) + 0);
		dif_p20_lo_offset[l++] = ((pb << 3) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p6 << 3) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p1 << 3) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pc << 3) + 1);
		dif_p20_lo_offset[l++] = ((p7 << 3) + 1);
		dif_p20_lo_offset[l++] = ((p2 << 3) + 1);
		dif_p20_lo_offset[l++] = ((pd << 3) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p8 << 3) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p3 << 3) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pe << 3) + 2);
		dif_p20_lo_offset[l++] = ((p9 << 3) + 2);
		dif_p20_lo_offset[l++] = ((p4 << 3) + 2);
		dif_p20_lo_offset[l++] = ((pf << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pa << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p5 << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = (( 0 << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pb << 3) + 3);
		dif_p20_lo_offset[l++] = ((p6 << 3) + 3);
		dif_p20_lo_offset[l++] = ((p1 << 3) + 3);
		dif_p20_lo_offset[l++] = ((pc << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p7 << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p2 << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pd << 3) + 4);
		dif_p20_lo_offset[l++] = ((p8 << 3) + 4);
		dif_p20_lo_offset[l++] = ((p3 << 3) + 4);
		dif_p20_lo_offset[l++] = ((pe << 3) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p9 << 3) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p4 << 3) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pf << 3) + 0);
		dif_p20_lo_offset[l++] = ((pa << 3) + 0);
		dif_p20_lo_offset[l++] = ((p5 << 3) + 0);

	   #endif	// sse2?

	// dif_offsets are w.r.to a-array, need 7 distinct sets of these, one for each DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		// Set 0: [0,1,3,2,6,7,4,5,c,d,f,e,9,8,a,b + p00],[9,8,a,b,f,e,d,c,3,2,1,0,4,5,7,6 + p10]
		l = 0;
		dif_offsets[0x00] =  0;		dif_offsets[0x10] = p9+p10;
		dif_offsets[0x01] = p1;		dif_offsets[0x11] = p8+p10;
		dif_offsets[0x02] = p3;		dif_offsets[0x12] = pa+p10;
		dif_offsets[0x03] = p2;		dif_offsets[0x13] = pb+p10;
		dif_offsets[0x04] = p6;		dif_offsets[0x14] = pf+p10;
		dif_offsets[0x05] = p7;		dif_offsets[0x15] = pe+p10;
		dif_offsets[0x06] = p4;		dif_offsets[0x16] = pd+p10;
		dif_offsets[0x07] = p5;		dif_offsets[0x17] = pc+p10;
		dif_offsets[0x08] = pc;		dif_offsets[0x18] = p3+p10;
		dif_offsets[0x09] = pd;		dif_offsets[0x19] = p2+p10;
		dif_offsets[0x0a] = pf;		dif_offsets[0x1a] = p1+p10;
		dif_offsets[0x0b] = pe;		dif_offsets[0x1b] =    p10;
		dif_offsets[0x0c] = p9;		dif_offsets[0x1c] = p4+p10;
		dif_offsets[0x0d] = p8;		dif_offsets[0x1d] = p5+p10;
		dif_offsets[0x0e] = pa;		dif_offsets[0x1e] = p7+p10;
		dif_offsets[0x0f] = pb;		dif_offsets[0x1f] = p6+p10;
		// Set 1: [6,7,4,5,3,2,1,0,9,8,a,b,f,e,d,c + p80],[f,e,d,c,a,b,8,9,4,5,7,6,1,0,2,3 + p90]
		l += 32;
		dif_offsets[l+0x00] = p6;		dif_offsets[l+0x10] = pf+p10;
		dif_offsets[l+0x01] = p7;		dif_offsets[l+0x11] = pe+p10;
		dif_offsets[l+0x02] = p4;		dif_offsets[l+0x12] = pd+p10;
		dif_offsets[l+0x03] = p5;		dif_offsets[l+0x13] = pc+p10;
		dif_offsets[l+0x04] = p3;		dif_offsets[l+0x14] = pa+p10;
		dif_offsets[l+0x05] = p2;		dif_offsets[l+0x15] = pb+p10;
		dif_offsets[l+0x06] = p1;		dif_offsets[l+0x16] = p8+p10;
		dif_offsets[l+0x07] =  0;		dif_offsets[l+0x17] = p9+p10;
		dif_offsets[l+0x08] = p9;		dif_offsets[l+0x18] = p4+p10;
		dif_offsets[l+0x09] = p8;		dif_offsets[l+0x19] = p5+p10;
		dif_offsets[l+0x0a] = pa;		dif_offsets[l+0x1a] = p7+p10;
		dif_offsets[l+0x0b] = pb;		dif_offsets[l+0x1b] = p6+p10;
		dif_offsets[l+0x0c] = pf;		dif_offsets[l+0x1c] = p1+p10;
		dif_offsets[l+0x0d] = pe;		dif_offsets[l+0x1d] =    p10;
		dif_offsets[l+0x0e] = pd;		dif_offsets[l+0x1e] = p2+p10;
		dif_offsets[l+0x0f] = pc;		dif_offsets[l+0x1f] = p3+p10;
		// Set 2: [3,2,1,0,4,5,7,6,f,e,d,c,a,b,8,9 + p70],[6,7,4,5,3,2,1,0,9,8,a,b,f,e,d,c + p60]
		l += 32;
		dif_offsets[l+0x00] = p3+p10;		dif_offsets[l+0x10] = p6;
		dif_offsets[l+0x01] = p2+p10;		dif_offsets[l+0x11] = p7;
		dif_offsets[l+0x02] = p1+p10;		dif_offsets[l+0x12] = p4;
		dif_offsets[l+0x03] =    p10;		dif_offsets[l+0x13] = p5;
		dif_offsets[l+0x04] = p4+p10;		dif_offsets[l+0x14] = p3;
		dif_offsets[l+0x05] = p5+p10;		dif_offsets[l+0x15] = p2;
		dif_offsets[l+0x06] = p7+p10;		dif_offsets[l+0x16] = p1;
		dif_offsets[l+0x07] = p6+p10;		dif_offsets[l+0x17] =  0;
		dif_offsets[l+0x08] = pf+p10;		dif_offsets[l+0x18] = p9;
		dif_offsets[l+0x09] = pe+p10;		dif_offsets[l+0x19] = p8;
		dif_offsets[l+0x0a] = pd+p10;		dif_offsets[l+0x1a] = pa;
		dif_offsets[l+0x0b] = pc+p10;		dif_offsets[l+0x1b] = pb;
		dif_offsets[l+0x0c] = pa+p10;		dif_offsets[l+0x1c] = pf;
		dif_offsets[l+0x0d] = pb+p10;		dif_offsets[l+0x1d] = pe;
		dif_offsets[l+0x0e] = p8+p10;		dif_offsets[l+0x1e] = pd;
		dif_offsets[l+0x0f] = p9+p10;		dif_offsets[l+0x1f] = pc;
		// Set 3: [c,d,f,e,9,8,a,b,6,7,4,5,3,2,1,0 + p40],[3,2,1,0,4,5,7,6,f,e,d,c,a,b,8,9 + p50]
		l += 32;
		dif_offsets[l+0x00] = pc;		dif_offsets[l+0x10] = p3+p10;
		dif_offsets[l+0x01] = pd;		dif_offsets[l+0x11] = p2+p10;
		dif_offsets[l+0x02] = pf;		dif_offsets[l+0x12] = p1+p10;
		dif_offsets[l+0x03] = pe;		dif_offsets[l+0x13] =    p10;
		dif_offsets[l+0x04] = p9;		dif_offsets[l+0x14] = p4+p10;
		dif_offsets[l+0x05] = p8;		dif_offsets[l+0x15] = p5+p10;
		dif_offsets[l+0x06] = pa;		dif_offsets[l+0x16] = p7+p10;
		dif_offsets[l+0x07] = pb;		dif_offsets[l+0x17] = p6+p10;
		dif_offsets[l+0x08] = p6;		dif_offsets[l+0x18] = pf+p10;
		dif_offsets[l+0x09] = p7;		dif_offsets[l+0x19] = pe+p10;
		dif_offsets[l+0x0a] = p4;		dif_offsets[l+0x1a] = pd+p10;
		dif_offsets[l+0x0b] = p5;		dif_offsets[l+0x1b] = pc+p10;
		dif_offsets[l+0x0c] = p3;		dif_offsets[l+0x1c] = pa+p10;
		dif_offsets[l+0x0d] = p2;		dif_offsets[l+0x1d] = pb+p10;
		dif_offsets[l+0x0e] = p1;		dif_offsets[l+0x1e] = p8+p10;
		dif_offsets[l+0x0f] =  0;		dif_offsets[l+0x1f] = p9+p10;
		// Set 4: [9,8,a,b,f,e,d,c,3,2,1,0,4,5,7,6 + p30],[c,d,f,e,9,8,a,b,6,7,4,5,3,2,1,0 + p20]
		l += 32;
		dif_offsets[l+0x00] = p9+p10;		dif_offsets[l+0x10] = pc;
		dif_offsets[l+0x01] = p8+p10;		dif_offsets[l+0x11] = pd;
		dif_offsets[l+0x02] = pa+p10;		dif_offsets[l+0x12] = pf;
		dif_offsets[l+0x03] = pb+p10;		dif_offsets[l+0x13] = pe;
		dif_offsets[l+0x04] = pf+p10;		dif_offsets[l+0x14] = p9;
		dif_offsets[l+0x05] = pe+p10;		dif_offsets[l+0x15] = p8;
		dif_offsets[l+0x06] = pd+p10;		dif_offsets[l+0x16] = pa;
		dif_offsets[l+0x07] = pc+p10;		dif_offsets[l+0x17] = pb;
		dif_offsets[l+0x08] = p3+p10;		dif_offsets[l+0x18] = p6;
		dif_offsets[l+0x09] = p2+p10;		dif_offsets[l+0x19] = p7;
		dif_offsets[l+0x0a] = p1+p10;		dif_offsets[l+0x1a] = p4;
		dif_offsets[l+0x0b] =    p10;		dif_offsets[l+0x1b] = p5;
		dif_offsets[l+0x0c] = p4+p10;		dif_offsets[l+0x1c] = p3;
		dif_offsets[l+0x0d] = p5+p10;		dif_offsets[l+0x1d] = p2;
		dif_offsets[l+0x0e] = p7+p10;		dif_offsets[l+0x1e] = p1;
		dif_offsets[l+0x0f] = p6+p10;		dif_offsets[l+0x1f] =  0;
	  #ifdef USE_SSE2
		// IN SIMD mode preshift all the above offsets << 3 to turn into double=array pointer offsets:
		for(l = 0; l < RADIX; l++) {
			dif_offsets[l] <<= 3;
		}
	  #endif

	/*** DIT indexing stuff: ***/
		// Init storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 2*9

	   #ifdef USE_SSE2
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p10*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:
		l = 0;
		dit_p20_cperms[l++] = 0x00<<1; dit_p20_cperms[l++] = 0x60<<1; dit_p20_cperms[l++] = 0x20<<1; dit_p20_cperms[l++] = 0x80<<1; dit_p20_cperms[l++] = 0x40<<1; dit_p20_cperms[l++] = 0x00<<1; dit_p20_cperms[l++] = 0x60<<1; dit_p20_cperms[l++] = 0x20<<1; dit_p20_cperms[l++] = 0x80<<1;
		dit_p20_cperms[l++] = 0x50<<1; dit_p20_cperms[l++] = 0x10<<1; dit_p20_cperms[l++] = 0x70<<1; dit_p20_cperms[l++] = 0x30<<1; dit_p20_cperms[l++] = 0x90<<1; dit_p20_cperms[l++] = 0x50<<1; dit_p20_cperms[l++] = 0x10<<1; dit_p20_cperms[l++] = 0x70<<1; dit_p20_cperms[l++] = 0x30<<1;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-5 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-4 index, bits <4:30> for the p0-f. In SIMD mode again replace p[0-f] with 0x[0-f]<<1 for use with
		// contig-local-mem, thus replace 'p' prefixes with 0x and 3-left-shifts with << 4 to account for the << 1:
		l = 0;
		dit_p20_lo_offset[l++] = ((0x0 << 4) + 0);
		dit_p20_lo_offset[l++] = ((0xf << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xe << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xd << 4) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xc << 4) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xb << 4) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xa << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x9 << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x8 << 4) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x7 << 4) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x6 << 4) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x5 << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x4 << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x3 << 4) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x2 << 4) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x1 << 4) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x0 << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xf << 4) + 0);
		dit_p20_lo_offset[l++] = ((0xe << 4) + 1);
		dit_p20_lo_offset[l++] = ((0xd << 4) + 2);
		dit_p20_lo_offset[l++] = ((0xc << 4) + 3);
		dit_p20_lo_offset[l++] = ((0xb << 4) + 4);
		dit_p20_lo_offset[l++] = ((0xa << 4) + 0);
		dit_p20_lo_offset[l++] = ((0x9 << 4) + 1);
		dit_p20_lo_offset[l++] = ((0x8 << 4) + 2);
		dit_p20_lo_offset[l++] = ((0x7 << 4) + 3);
		dit_p20_lo_offset[l++] = ((0x6 << 4) + 4);
		dit_p20_lo_offset[l++] = ((0x5 << 4) + 0);
		dit_p20_lo_offset[l++] = ((0x4 << 4) + 1);
		dit_p20_lo_offset[l++] = ((0x3 << 4) + 2);
		dit_p20_lo_offset[l++] = ((0x2 << 4) + 3);
		dit_p20_lo_offset[l++] = ((0x1 << 4) + 4);

	   #else

		// Init storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 2*9
		l = 0;
		dit_p20_cperms[l++] =   0; dit_p20_cperms[l++] = p60; dit_p20_cperms[l++] = p20; dit_p20_cperms[l++] = p80; dit_p20_cperms[l++] = p40; dit_p20_cperms[l++] =   0; dit_p20_cperms[l++] = p60; dit_p20_cperms[l++] = p20; dit_p20_cperms[l++] = p80;
		dit_p20_cperms[l++] = p50; dit_p20_cperms[l++] = p10; dit_p20_cperms[l++] = p70; dit_p20_cperms[l++] = p30; dit_p20_cperms[l++] = p90; dit_p20_cperms[l++] = p50; dit_p20_cperms[l++] = p10; dit_p20_cperms[l++] = p70; dit_p20_cperms[l++] = p30;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-5 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-4 index, bits <4:30> for the p0-f:
		l = 0;
		dit_p20_lo_offset[l++] = (( 0 << 3) + 0);
		dit_p20_lo_offset[l++] = ((pf << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pe << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pd << 3) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pc << 3) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pb << 3) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pa << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p9 << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p8 << 3) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p7 << 3) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p6 << 3) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p5 << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p4 << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p3 << 3) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p2 << 3) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p1 << 3) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = (( 0 << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pf << 3) + 0);
		dit_p20_lo_offset[l++] = ((pe << 3) + 1);
		dit_p20_lo_offset[l++] = ((pd << 3) + 2);
		dit_p20_lo_offset[l++] = ((pc << 3) + 3);
		dit_p20_lo_offset[l++] = ((pb << 3) + 4);
		dit_p20_lo_offset[l++] = ((pa << 3) + 0);
		dit_p20_lo_offset[l++] = ((p9 << 3) + 1);
		dit_p20_lo_offset[l++] = ((p8 << 3) + 2);
		dit_p20_lo_offset[l++] = ((p7 << 3) + 3);
		dit_p20_lo_offset[l++] = ((p6 << 3) + 4);
		dit_p20_lo_offset[l++] = ((p5 << 3) + 0);
		dit_p20_lo_offset[l++] = ((p4 << 3) + 1);
		dit_p20_lo_offset[l++] = ((p3 << 3) + 2);
		dit_p20_lo_offset[l++] = ((p2 << 3) + 3);
		dit_p20_lo_offset[l++] = ((p1 << 3) + 4);

	   #endif	// sse2?

	// dit_offsets are w.r.to a-array, need 5 distinct sets of these, one for each radix-32 DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		// Set 0: [0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00],[f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10]
		l = 0;
		dit_offsets[0x00] =  0;		dit_offsets[0x10] = pf+p10;
		dit_offsets[0x01] = p1;		dit_offsets[0x11] = pe+p10;
		dit_offsets[0x02] = p3;		dit_offsets[0x12] = pd+p10;
		dit_offsets[0x03] = p2;		dit_offsets[0x13] = pc+p10;
		dit_offsets[0x04] = p7;		dit_offsets[0x14] = pb+p10;
		dit_offsets[0x05] = p6;		dit_offsets[0x15] = pa+p10;
		dit_offsets[0x06] = p5;		dit_offsets[0x16] = p9+p10;
		dit_offsets[0x07] = p4;		dit_offsets[0x17] = p8+p10;
		dit_offsets[0x08] = pf;		dit_offsets[0x18] = p7+p10;
		dit_offsets[0x09] = pe;		dit_offsets[0x19] = p6+p10;
		dit_offsets[0x0a] = pd;		dit_offsets[0x1a] = p5+p10;
		dit_offsets[0x0b] = pc;		dit_offsets[0x1b] = p4+p10;
		dit_offsets[0x0c] = pb;		dit_offsets[0x1c] = p3+p10;
		dit_offsets[0x0d] = pa;		dit_offsets[0x1d] = p2+p10;
		dit_offsets[0x0e] = p9;		dit_offsets[0x1e] = p1+p10;
		dit_offsets[0x0f] = p8;		dit_offsets[0x1f] =    p10;
		// Set 1: [3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p70],[3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p60] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = p3+p10;		dit_offsets[l+0x10] = p3;
		dit_offsets[l+0x01] = p2+p10;		dit_offsets[l+0x11] = p2;
		dit_offsets[l+0x02] = p1+p10;		dit_offsets[l+0x12] = p1;
		dit_offsets[l+0x03] =    p10;		dit_offsets[l+0x13] =  0;
		dit_offsets[l+0x04] = p5+p10;		dit_offsets[l+0x14] = p5;
		dit_offsets[l+0x05] = p4+p10;		dit_offsets[l+0x15] = p4;
		dit_offsets[l+0x06] = p6+p10;		dit_offsets[l+0x16] = p6;
		dit_offsets[l+0x07] = p7+p10;		dit_offsets[l+0x17] = p7;
		dit_offsets[l+0x08] = pd+p10;		dit_offsets[l+0x18] = pd;
		dit_offsets[l+0x09] = pc+p10;		dit_offsets[l+0x19] = pc;
		dit_offsets[l+0x0a] = pe+p10;		dit_offsets[l+0x1a] = pe;
		dit_offsets[l+0x0b] = pf+p10;		dit_offsets[l+0x1b] = pf;
		dit_offsets[l+0x0c] = p9+p10;		dit_offsets[l+0x1c] = p9;
		dit_offsets[l+0x0d] = p8+p10;		dit_offsets[l+0x1d] = p8;
		dit_offsets[l+0x0e] = pa+p10;		dit_offsets[l+0x1e] = pa;
		dit_offsets[l+0x0f] = pb+p10;		dit_offsets[l+0x1f] = pb;
		// Set 2: [9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p30],[9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p20] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = p9+p10;		dit_offsets[l+0x10] = p9;
		dit_offsets[l+0x01] = p8+p10;		dit_offsets[l+0x11] = p8;
		dit_offsets[l+0x02] = pa+p10;		dit_offsets[l+0x12] = pa;
		dit_offsets[l+0x03] = pb+p10;		dit_offsets[l+0x13] = pb;
		dit_offsets[l+0x04] = pe+p10;		dit_offsets[l+0x14] = pe;
		dit_offsets[l+0x05] = pf+p10;		dit_offsets[l+0x15] = pf;
		dit_offsets[l+0x06] = pc+p10;		dit_offsets[l+0x16] = pc;
		dit_offsets[l+0x07] = pd+p10;		dit_offsets[l+0x17] = pd;
		dit_offsets[l+0x08] = p1+p10;		dit_offsets[l+0x18] = p1;
		dit_offsets[l+0x09] =    p10;		dit_offsets[l+0x19] =  0;
		dit_offsets[l+0x0a] = p2+p10;		dit_offsets[l+0x1a] = p2;
		dit_offsets[l+0x0b] = p3+p10;		dit_offsets[l+0x1b] = p3;
		dit_offsets[l+0x0c] = p6+p10;		dit_offsets[l+0x1c] = p6;
		dit_offsets[l+0x0d] = p7+p10;		dit_offsets[l+0x1d] = p7;
		dit_offsets[l+0x0e] = p4+p10;		dit_offsets[l+0x1e] = p4;
		dit_offsets[l+0x0f] = p5+p10;		dit_offsets[l+0x1f] = p5;
		// Set 3: [6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + p80],[a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p90] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = p6;		dit_offsets[l+0x10] = pa+p10;
		dit_offsets[l+0x01] = p7;		dit_offsets[l+0x11] = pb+p10;
		dit_offsets[l+0x02] = p4;		dit_offsets[l+0x12] = p8+p10;
		dit_offsets[l+0x03] = p5;		dit_offsets[l+0x13] = p9+p10;
		dit_offsets[l+0x04] = p2;		dit_offsets[l+0x14] = pc+p10;
		dit_offsets[l+0x05] = p3;		dit_offsets[l+0x15] = pd+p10;
		dit_offsets[l+0x06] =  0;		dit_offsets[l+0x16] = pf+p10;
		dit_offsets[l+0x07] = p1;		dit_offsets[l+0x17] = pe+p10;
		dit_offsets[l+0x08] = pa;		dit_offsets[l+0x18] = p2+p10;
		dit_offsets[l+0x09] = pb;		dit_offsets[l+0x19] = p3+p10;
		dit_offsets[l+0x0a] = p8;		dit_offsets[l+0x1a] =    p10;
		dit_offsets[l+0x0b] = p9;		dit_offsets[l+0x1b] = p1+p10;
		dit_offsets[l+0x0c] = pc;		dit_offsets[l+0x1c] = p4+p10;
		dit_offsets[l+0x0d] = pd;		dit_offsets[l+0x1d] = p5+p10;
		dit_offsets[l+0x0e] = pf;		dit_offsets[l+0x1e] = p7+p10;
		dit_offsets[l+0x0f] = pe;		dit_offsets[l+0x1f] = p6+p10;
		// Set 4: [c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p40],[4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + p50]] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = pc;		dit_offsets[l+0x10] = p4+p10;
		dit_offsets[l+0x01] = pd;		dit_offsets[l+0x11] = p5+p10;
		dit_offsets[l+0x02] = pf;		dit_offsets[l+0x12] = p7+p10;
		dit_offsets[l+0x03] = pe;		dit_offsets[l+0x13] = p6+p10;
		dit_offsets[l+0x04] = p8;		dit_offsets[l+0x14] =    p10;
		dit_offsets[l+0x05] = p9;		dit_offsets[l+0x15] = p1+p10;
		dit_offsets[l+0x06] = pb;		dit_offsets[l+0x16] = p3+p10;
		dit_offsets[l+0x07] = pa;		dit_offsets[l+0x17] = p2+p10;
		dit_offsets[l+0x08] = p4;		dit_offsets[l+0x18] = p8+p10;
		dit_offsets[l+0x09] = p5;		dit_offsets[l+0x19] = p9+p10;
		dit_offsets[l+0x0a] = p7;		dit_offsets[l+0x1a] = pb+p10;
		dit_offsets[l+0x0b] = p6;		dit_offsets[l+0x1b] = pa+p10;
		dit_offsets[l+0x0c] =  0;		dit_offsets[l+0x1c] = pf+p10;
		dit_offsets[l+0x0d] = p1;		dit_offsets[l+0x1d] = pe+p10;
		dit_offsets[l+0x0e] = p3;		dit_offsets[l+0x1e] = pd+p10;
		dit_offsets[l+0x0f] = p2;		dit_offsets[l+0x1f] = pc+p10;
	  #ifdef USE_SSE2
		// IN SIMD mode preshift all the above offsets << 3 to turn into double=array pointer offsets:
		for(l = 0; l < RADIX; l++) {
			dit_offsets[l] <<= 3;
		}
	  #endif

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

/*...The radix-160 final DIT pass is here.	*/

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
		#include "radix160_main_carry_loop.h"

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
		ASSERT(0x0 == cy160_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
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

/****************/

void radix160_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-160 complex DIF FFT pass on the data in the length-N real vector A.
*/
	int i,j,j1,j2,jt,jp;
	// p-indexing is hexadecimal here:
	static int NDIVR,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90, first_entry=TRUE;
	static int t_offsets[32], dif_offsets[RADIX];
	// Need storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 2*9 elts:
	static int dif_p20_cperms[18], dif_p20_lo_offset[32];
	const double cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
				cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
				s2  =  0.95105651629515357211,	/*  sin(u) */
				ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
				ss2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it;
	struct complex t[RADIX], *tptr;

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
		NDIVR >>= 4;			p90 += ( (p90 >> DAT_BITS) << PAD_BITS );

		// Set array offsets for radix-32 DFT in/outputs:
		// t_offsets w.r.to: t-array, same for all 5 DFTs:
		t_offsets[0x00] = 0x00<<1;	t_offsets[0x10] = 0x10<<1;
		t_offsets[0x01] = 0x01<<1;	t_offsets[0x11] = 0x11<<1;
		t_offsets[0x02] = 0x02<<1;	t_offsets[0x12] = 0x12<<1;
		t_offsets[0x03] = 0x03<<1;	t_offsets[0x13] = 0x13<<1;
		t_offsets[0x04] = 0x04<<1;	t_offsets[0x14] = 0x14<<1;
		t_offsets[0x05] = 0x05<<1;	t_offsets[0x15] = 0x15<<1;
		t_offsets[0x06] = 0x06<<1;	t_offsets[0x16] = 0x16<<1;
		t_offsets[0x07] = 0x07<<1;	t_offsets[0x17] = 0x17<<1;
		t_offsets[0x08] = 0x08<<1;	t_offsets[0x18] = 0x18<<1;
		t_offsets[0x09] = 0x09<<1;	t_offsets[0x19] = 0x19<<1;
		t_offsets[0x0a] = 0x0a<<1;	t_offsets[0x1a] = 0x1a<<1;
		t_offsets[0x0b] = 0x0b<<1;	t_offsets[0x1b] = 0x1b<<1;
		t_offsets[0x0c] = 0x0c<<1;	t_offsets[0x1c] = 0x1c<<1;
		t_offsets[0x0d] = 0x0d<<1;	t_offsets[0x1d] = 0x1d<<1;
		t_offsets[0x0e] = 0x0e<<1;	t_offsets[0x1e] = 0x1e<<1;
		t_offsets[0x0f] = 0x0f<<1;	t_offsets[0x1f] = 0x1f<<1;

		// Init storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 2*9
		i = 0;
		dif_p20_cperms[i++] =   0; dif_p20_cperms[i++] = p80; dif_p20_cperms[i++] = p60; dif_p20_cperms[i++] = p40; dif_p20_cperms[i++] = p20; dif_p20_cperms[i++] =   0; dif_p20_cperms[i++] = p80; dif_p20_cperms[i++] = p60; dif_p20_cperms[i++] = p40;
		dif_p20_cperms[i++] = p90; dif_p20_cperms[i++] = p70; dif_p20_cperms[i++] = p50; dif_p20_cperms[i++] = p30; dif_p20_cperms[i++] = p10; dif_p20_cperms[i++] = p90; dif_p20_cperms[i++] = p70; dif_p20_cperms[i++] = p50; dif_p20_cperms[i++] = p30;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-5 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-4 index, bits <4:30> for the p0-f:
		i = 0;
		dif_p20_lo_offset[i++] = (( 0 << 3) + 0);
		dif_p20_lo_offset[i++] = ((pb << 3) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p6 << 3) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p1 << 3) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pc << 3) + 1);
		dif_p20_lo_offset[i++] = ((p7 << 3) + 1);
		dif_p20_lo_offset[i++] = ((p2 << 3) + 1);
		dif_p20_lo_offset[i++] = ((pd << 3) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p8 << 3) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p3 << 3) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pe << 3) + 2);
		dif_p20_lo_offset[i++] = ((p9 << 3) + 2);
		dif_p20_lo_offset[i++] = ((p4 << 3) + 2);
		dif_p20_lo_offset[i++] = ((pf << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pa << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p5 << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = (( 0 << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pb << 3) + 3);
		dif_p20_lo_offset[i++] = ((p6 << 3) + 3);
		dif_p20_lo_offset[i++] = ((p1 << 3) + 3);
		dif_p20_lo_offset[i++] = ((pc << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p7 << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p2 << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pd << 3) + 4);
		dif_p20_lo_offset[i++] = ((p8 << 3) + 4);
		dif_p20_lo_offset[i++] = ((p3 << 3) + 4);
		dif_p20_lo_offset[i++] = ((pe << 3) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p9 << 3) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p4 << 3) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pf << 3) + 0);
		dif_p20_lo_offset[i++] = ((pa << 3) + 0);
		dif_p20_lo_offset[i++] = ((p5 << 3) + 0);

	// dif_offsets are w.r.to a-array, need 7 distinct sets of these, one for each DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		// Set 0: [0,1,3,2,6,7,4,5,c,d,f,e,9,8,a,b + p00],[9,8,a,b,f,e,d,c,3,2,1,0,4,5,7,6 + p10]
		i = 0;
		dif_offsets[0x00] =  0;		dif_offsets[0x10] = p9+p10;
		dif_offsets[0x01] = p1;		dif_offsets[0x11] = p8+p10;
		dif_offsets[0x02] = p3;		dif_offsets[0x12] = pa+p10;
		dif_offsets[0x03] = p2;		dif_offsets[0x13] = pb+p10;
		dif_offsets[0x04] = p6;		dif_offsets[0x14] = pf+p10;
		dif_offsets[0x05] = p7;		dif_offsets[0x15] = pe+p10;
		dif_offsets[0x06] = p4;		dif_offsets[0x16] = pd+p10;
		dif_offsets[0x07] = p5;		dif_offsets[0x17] = pc+p10;
		dif_offsets[0x08] = pc;		dif_offsets[0x18] = p3+p10;
		dif_offsets[0x09] = pd;		dif_offsets[0x19] = p2+p10;
		dif_offsets[0x0a] = pf;		dif_offsets[0x1a] = p1+p10;
		dif_offsets[0x0b] = pe;		dif_offsets[0x1b] =    p10;
		dif_offsets[0x0c] = p9;		dif_offsets[0x1c] = p4+p10;
		dif_offsets[0x0d] = p8;		dif_offsets[0x1d] = p5+p10;
		dif_offsets[0x0e] = pa;		dif_offsets[0x1e] = p7+p10;
		dif_offsets[0x0f] = pb;		dif_offsets[0x1f] = p6+p10;
		// Set 1: [6,7,4,5,3,2,1,0,9,8,a,b,f,e,d,c + p80],[f,e,d,c,a,b,8,9,4,5,7,6,1,0,2,3 + p90]
		i += 32;
		dif_offsets[i+0x00] = p6;		dif_offsets[i+0x10] = pf+p10;
		dif_offsets[i+0x01] = p7;		dif_offsets[i+0x11] = pe+p10;
		dif_offsets[i+0x02] = p4;		dif_offsets[i+0x12] = pd+p10;
		dif_offsets[i+0x03] = p5;		dif_offsets[i+0x13] = pc+p10;
		dif_offsets[i+0x04] = p3;		dif_offsets[i+0x14] = pa+p10;
		dif_offsets[i+0x05] = p2;		dif_offsets[i+0x15] = pb+p10;
		dif_offsets[i+0x06] = p1;		dif_offsets[i+0x16] = p8+p10;
		dif_offsets[i+0x07] =  0;		dif_offsets[i+0x17] = p9+p10;
		dif_offsets[i+0x08] = p9;		dif_offsets[i+0x18] = p4+p10;
		dif_offsets[i+0x09] = p8;		dif_offsets[i+0x19] = p5+p10;
		dif_offsets[i+0x0a] = pa;		dif_offsets[i+0x1a] = p7+p10;
		dif_offsets[i+0x0b] = pb;		dif_offsets[i+0x1b] = p6+p10;
		dif_offsets[i+0x0c] = pf;		dif_offsets[i+0x1c] = p1+p10;
		dif_offsets[i+0x0d] = pe;		dif_offsets[i+0x1d] =    p10;
		dif_offsets[i+0x0e] = pd;		dif_offsets[i+0x1e] = p2+p10;
		dif_offsets[i+0x0f] = pc;		dif_offsets[i+0x1f] = p3+p10;
		// Set 2: [3,2,1,0,4,5,7,6,f,e,d,c,a,b,8,9 + p70],[6,7,4,5,3,2,1,0,9,8,a,b,f,e,d,c + p60]
		i += 32;
		dif_offsets[i+0x00] = p3+p10;		dif_offsets[i+0x10] = p6;
		dif_offsets[i+0x01] = p2+p10;		dif_offsets[i+0x11] = p7;
		dif_offsets[i+0x02] = p1+p10;		dif_offsets[i+0x12] = p4;
		dif_offsets[i+0x03] =    p10;		dif_offsets[i+0x13] = p5;
		dif_offsets[i+0x04] = p4+p10;		dif_offsets[i+0x14] = p3;
		dif_offsets[i+0x05] = p5+p10;		dif_offsets[i+0x15] = p2;
		dif_offsets[i+0x06] = p7+p10;		dif_offsets[i+0x16] = p1;
		dif_offsets[i+0x07] = p6+p10;		dif_offsets[i+0x17] =  0;
		dif_offsets[i+0x08] = pf+p10;		dif_offsets[i+0x18] = p9;
		dif_offsets[i+0x09] = pe+p10;		dif_offsets[i+0x19] = p8;
		dif_offsets[i+0x0a] = pd+p10;		dif_offsets[i+0x1a] = pa;
		dif_offsets[i+0x0b] = pc+p10;		dif_offsets[i+0x1b] = pb;
		dif_offsets[i+0x0c] = pa+p10;		dif_offsets[i+0x1c] = pf;
		dif_offsets[i+0x0d] = pb+p10;		dif_offsets[i+0x1d] = pe;
		dif_offsets[i+0x0e] = p8+p10;		dif_offsets[i+0x1e] = pd;
		dif_offsets[i+0x0f] = p9+p10;		dif_offsets[i+0x1f] = pc;
		// Set 3: [c,d,f,e,9,8,a,b,6,7,4,5,3,2,1,0 + p40],[3,2,1,0,4,5,7,6,f,e,d,c,a,b,8,9 + p50]
		i += 32;
		dif_offsets[i+0x00] = pc;		dif_offsets[i+0x10] = p3+p10;
		dif_offsets[i+0x01] = pd;		dif_offsets[i+0x11] = p2+p10;
		dif_offsets[i+0x02] = pf;		dif_offsets[i+0x12] = p1+p10;
		dif_offsets[i+0x03] = pe;		dif_offsets[i+0x13] =    p10;
		dif_offsets[i+0x04] = p9;		dif_offsets[i+0x14] = p4+p10;
		dif_offsets[i+0x05] = p8;		dif_offsets[i+0x15] = p5+p10;
		dif_offsets[i+0x06] = pa;		dif_offsets[i+0x16] = p7+p10;
		dif_offsets[i+0x07] = pb;		dif_offsets[i+0x17] = p6+p10;
		dif_offsets[i+0x08] = p6;		dif_offsets[i+0x18] = pf+p10;
		dif_offsets[i+0x09] = p7;		dif_offsets[i+0x19] = pe+p10;
		dif_offsets[i+0x0a] = p4;		dif_offsets[i+0x1a] = pd+p10;
		dif_offsets[i+0x0b] = p5;		dif_offsets[i+0x1b] = pc+p10;
		dif_offsets[i+0x0c] = p3;		dif_offsets[i+0x1c] = pa+p10;
		dif_offsets[i+0x0d] = p2;		dif_offsets[i+0x1d] = pb+p10;
		dif_offsets[i+0x0e] = p1;		dif_offsets[i+0x1e] = p8+p10;
		dif_offsets[i+0x0f] =  0;		dif_offsets[i+0x1f] = p9+p10;
		// Set 4: [9,8,a,b,f,e,d,c,3,2,1,0,4,5,7,6 + p30],[c,d,f,e,9,8,a,b,6,7,4,5,3,2,1,0 + p20]
		i += 32;
		dif_offsets[i+0x00] = p9+p10;		dif_offsets[i+0x10] = pc;
		dif_offsets[i+0x01] = p8+p10;		dif_offsets[i+0x11] = pd;
		dif_offsets[i+0x02] = pa+p10;		dif_offsets[i+0x12] = pf;
		dif_offsets[i+0x03] = pb+p10;		dif_offsets[i+0x13] = pe;
		dif_offsets[i+0x04] = pf+p10;		dif_offsets[i+0x14] = p9;
		dif_offsets[i+0x05] = pe+p10;		dif_offsets[i+0x15] = p8;
		dif_offsets[i+0x06] = pd+p10;		dif_offsets[i+0x16] = pa;
		dif_offsets[i+0x07] = pc+p10;		dif_offsets[i+0x17] = pb;
		dif_offsets[i+0x08] = p3+p10;		dif_offsets[i+0x18] = p6;
		dif_offsets[i+0x09] = p2+p10;		dif_offsets[i+0x19] = p7;
		dif_offsets[i+0x0a] = p1+p10;		dif_offsets[i+0x1a] = p4;
		dif_offsets[i+0x0b] =    p10;		dif_offsets[i+0x1b] = p5;
		dif_offsets[i+0x0c] = p4+p10;		dif_offsets[i+0x1c] = p3;
		dif_offsets[i+0x0d] = p5+p10;		dif_offsets[i+0x1d] = p2;
		dif_offsets[i+0x0e] = p7+p10;		dif_offsets[i+0x1e] = p1;
		dif_offsets[i+0x0f] = p6+p10;		dif_offsets[i+0x1f] =  0;
	}

/*...The radix-160 pass is here.	*/

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

	//...gather the needed data (160 64-bit complex) and do 32 radix-5 transforms:
	/*
	Twiddleless version arranges 32 sets of radix-5 DFT inputs as follows: 0 in upper left corner,
	decrement (mod 160 = 0xa0) 32 (= 0x20) horizontally and 5 vertically.

	DIF/DIT input-scramble array =

	Indexing in hex for clarity and using [evn|odd]0-4 notation in the rightmost column to flag reusable
	5-perms [in fact simple circular (0-4)-element shifts of evn = [00,80,60,40,20] and odd = [90,70,50,30,10]:

		00,80,60,40,20		00,80,60,40,20 + p0		[evn0] + p0
		9b,7b,5b,3b,1b		90,70,50,30,10 + pb		[odd0] + pb
		96,76,56,36,16		90,70,50,30,10 + p6		[odd0] + p6
		91,71,51,31,11		90,70,50,30,10 + p1		[odd0] + p1
		8c,6c,4c,2c,0c		80,60,40,20,00 + pc		[evn1] + pc
		87,67,47,27,07		80,60,40,20,00 + p7		[evn1] + p7
		82,62,42,22,02		80,60,40,20,00 + p2		[evn1] + p2
		7d,5d,3d,1d,9d		70,50,30,10,90 + pd		[odd1] + pd
		78,58,38,18,98		70,50,30,10,90 + p8		[odd1] + p8
		73,53,33,13,93		70,50,30,10,90 + p3		[odd1] + p3
		6e,4e,2e,0e,8e		60,40,20,00,80 + pe		[evn2] + pe
		69,49,29,09,89		60,40,20,00,80 + p9		[evn2] + p9
		64,44,24,04,84		60,40,20,00,80 + p4		[evn2] + p4
		5f,3f,1f,9f,7f		50,30,10,90,70 + pf		[odd2] + pf
		5a,3a,1a,9a,7a		50,30,10,90,70 + pa		[odd2] + pa
		55,35,15,95,75		50,30,10,90,70 + p5		[odd2] + p5
		50,30,10,90,70	=	50,30,10,90,70 + p0	=	[odd2] + p0	<<< p0-5 pattern repeats here
		4b,2b,0b,8b,6b		40,20,00,80,60 + pb		[evn3] + pb
		46,26,06,86,66		40,20,00,80,60 + p6		[evn3] + p6
		41,21,01,81,61		40,20,00,80,60 + p1		[evn3] + p1
		3c,1c,9c,7c,5c		30,10,90,70,50 + pc		[odd3] + pc
		37,17,97,77,57		30,10,90,70,50 + p7		[odd3] + p7
		32,12,92,72,52		30,10,90,70,50 + p2		[odd3] + p2
		2d,0d,8d,6d,4d		20,00,80,60,40 + pd		[evn4] + pd
		28,08,88,68,48		20,00,80,60,40 + p8		[evn4] + p8
		23,03,83,63,43		20,00,80,60,40 + p3		[evn4] + p3
		1e,9e,7e,5e,3e		10,90,70,50,30 + pe		[odd4] + pe
		19,99,79,59,39		10,90,70,50,30 + p9		[odd4] + p9
		14,94,74,54,34		10,90,70,50,30 + p4		[odd4] + p4
		0f,8f,6f,4f,2f		00,80,60,40,20 + pf		[evn0] + pf
		0a,8a,6a,4a,2a		00,80,60,40,20 + pa		[evn0] + pa
		05,85,65,45,25		00,80,60,40,20 + p5		[evn0] + p5
	*/
		tptr = t;
		for(i = 0; i < 32; i++) {
			int k = dif_p20_lo_offset[i];
			// Extract index (in [0-4]) into circ-shift array used for high parts of p-mults. The [0-4] value is
			// in low 3 bits of k; the "which length-9 half of the dif_p20_cperms array?" selector is via (k < 0):
			int ic = ((-(k < 0)) & 9)	// +/- sign on k puts us into lower/upper half of the cshift array (base index 0/9)
						+ (k & 0x7);	// ...and low 3 bits give the element index w.r.to the array-half in question.
			int k0 = dif_p20_cperms[ic], k1 = dif_p20_cperms[ic+1], k2 = dif_p20_cperms[ic+2], k3 = dif_p20_cperms[ic+3], k4 = dif_p20_cperms[ic+4];
			// Extract Low part, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-7 DFTs:
			k = (k & 0x7fffffff) >> 3;
			jt = j1+k; jp = j2+k;
			RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],
				tptr->re,tptr->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x80)->re,(tptr+0x80)->im,
				rt,it);
			tptr++;
		}

	/*...and now do 5 radix-32 transforms;
	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF]
	to properly permute the outputs of the radix-32 DFTs to the required ordering, which in terms of our p-offsets is

		00,01,03,02,06,07,04,05,0c,0d,0f,0e,09,08,0a,0b,19,18,1a,1b,1f,1e,1d,1c,13,12,11,10,14,15,17,16
		86,87,84,85,83,82,81,80,89,88,8a,8b,8f,8e,8d,8c,9f,9e,9d,9c,9a,9b,98,99,94,95,97,96,91,90,92,93
		73,72,71,70,74,75,77,76,7f,7e,7d,7c,7a,7b,78,79,66,67,64,65,63,62,61,60,69,68,6a,6b,6f,6e,6d,6c
		4c,4d,4f,4e,49,48,4a,4b,46,47,44,45,43,42,41,40,53,52,51,50,54,55,57,56,5f,5e,5d,5c,5a,5b,58,59
		39,38,3a,3b,3f,3e,3d,3c,33,32,31,30,34,35,37,36,2c,2d,2f,2e,29,28,2a,2b,26,27,24,25,23,22,21,20

		[0,1,3,2,6,7,4,5,c,d,f,e,9,8,a,b + p00],[9,8,a,b,f,e,d,c,3,2,1,0,4,5,7,6 + p10]
		[6,7,4,5,3,2,1,0,9,8,a,b,f,e,d,c + p80],[f,e,d,c,a,b,8,9,4,5,7,6,1,0,2,3 + p90]
	=	[3,2,1,0,4,5,7,6,f,e,d,c,a,b,8,9 + p70],[6,7,4,5,3,2,1,0,9,8,a,b,f,e,d,c + p60]
		[c,d,f,e,9,8,a,b,6,7,4,5,3,2,1,0 + p40],[3,2,1,0,4,5,7,6,f,e,d,c,a,b,8,9 + p50]
		[9,8,a,b,f,e,d,c,3,2,1,0,4,5,7,6 + p30],[c,d,f,e,9,8,a,b,6,7,4,5,3,2,1,0 + p20]
	*/
		tptr = t;
		jt = j1    ; RADIX_32_DIF((double *)tptr,t_offsets,1, (a+jt),dif_offsets     ,RE_IM_STRIDE);	tptr += 32;
		jt = j1+p80; RADIX_32_DIF((double *)tptr,t_offsets,1, (a+jt),dif_offsets+0x20,RE_IM_STRIDE);	tptr += 32;
		jt = j1+p60; RADIX_32_DIF((double *)tptr,t_offsets,1, (a+jt),dif_offsets+0x40,RE_IM_STRIDE);	tptr += 32;
		jt = j1+p40; RADIX_32_DIF((double *)tptr,t_offsets,1, (a+jt),dif_offsets+0x60,RE_IM_STRIDE);	tptr += 32;
		jt = j1+p20; RADIX_32_DIF((double *)tptr,t_offsets,1, (a+jt),dif_offsets+0x80,RE_IM_STRIDE);
	}
}

/***************/

void radix160_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-160 complex DIT FFT pass on the data in the length-N real vector A.
*/
	int i,j,j1,j2,jt,jp;
	// p-indexing is hexadecimal here:
	static int NDIVR,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90, first_entry=TRUE;
	static int t_offsets[32], dit_offsets[RADIX];
	// Need storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 2*9 elts:
	static int dit_p20_cperms[18], dit_p20_lo_offset[32];
	const double cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
				cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
				s2  =  0.95105651629515357211,	/*  sin(u) */
				ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
				ss2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it;
	struct complex t[RADIX], *tptr;

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
		NDIVR >>= 4;			p90 += ( (p90 >> DAT_BITS) << PAD_BITS );

	// Set array offsets for radix-32 outputs.
		// t_offsets w.r.to: t-array, same for all 5 DFTs:
		t_offsets[0x00] = 0x00<<1;	t_offsets[0x10] = 0x10<<1;
		t_offsets[0x01] = 0x01<<1;	t_offsets[0x11] = 0x11<<1;
		t_offsets[0x02] = 0x02<<1;	t_offsets[0x12] = 0x12<<1;
		t_offsets[0x03] = 0x03<<1;	t_offsets[0x13] = 0x13<<1;
		t_offsets[0x04] = 0x04<<1;	t_offsets[0x14] = 0x14<<1;
		t_offsets[0x05] = 0x05<<1;	t_offsets[0x15] = 0x15<<1;
		t_offsets[0x06] = 0x06<<1;	t_offsets[0x16] = 0x16<<1;
		t_offsets[0x07] = 0x07<<1;	t_offsets[0x17] = 0x17<<1;
		t_offsets[0x08] = 0x08<<1;	t_offsets[0x18] = 0x18<<1;
		t_offsets[0x09] = 0x09<<1;	t_offsets[0x19] = 0x19<<1;
		t_offsets[0x0a] = 0x0a<<1;	t_offsets[0x1a] = 0x1a<<1;
		t_offsets[0x0b] = 0x0b<<1;	t_offsets[0x1b] = 0x1b<<1;
		t_offsets[0x0c] = 0x0c<<1;	t_offsets[0x1c] = 0x1c<<1;
		t_offsets[0x0d] = 0x0d<<1;	t_offsets[0x1d] = 0x1d<<1;
		t_offsets[0x0e] = 0x0e<<1;	t_offsets[0x1e] = 0x1e<<1;
		t_offsets[0x0f] = 0x0f<<1;	t_offsets[0x1f] = 0x1f<<1;

		// Init storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 2*9
		i = 0;
		dit_p20_cperms[i++] =   0; dit_p20_cperms[i++] = p60; dit_p20_cperms[i++] = p20; dit_p20_cperms[i++] = p80; dit_p20_cperms[i++] = p40; dit_p20_cperms[i++] =   0; dit_p20_cperms[i++] = p60; dit_p20_cperms[i++] = p20; dit_p20_cperms[i++] = p80;
		dit_p20_cperms[i++] = p50; dit_p20_cperms[i++] = p10; dit_p20_cperms[i++] = p70; dit_p20_cperms[i++] = p30; dit_p20_cperms[i++] = p90; dit_p20_cperms[i++] = p50; dit_p20_cperms[i++] = p10; dit_p20_cperms[i++] = p70; dit_p20_cperms[i++] = p30;

		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-5 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-4 index, bits <4:30> for the p0-f:
		i = 0;
		dit_p20_lo_offset[i++] = (( 0 << 3) + 0);
		dit_p20_lo_offset[i++] = ((pf << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pe << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pd << 3) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pc << 3) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pb << 3) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pa << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p9 << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p8 << 3) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p7 << 3) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p6 << 3) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p5 << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p4 << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p3 << 3) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p2 << 3) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p1 << 3) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = (( 0 << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pf << 3) + 0);
		dit_p20_lo_offset[i++] = ((pe << 3) + 1);
		dit_p20_lo_offset[i++] = ((pd << 3) + 2);
		dit_p20_lo_offset[i++] = ((pc << 3) + 3);
		dit_p20_lo_offset[i++] = ((pb << 3) + 4);
		dit_p20_lo_offset[i++] = ((pa << 3) + 0);
		dit_p20_lo_offset[i++] = ((p9 << 3) + 1);
		dit_p20_lo_offset[i++] = ((p8 << 3) + 2);
		dit_p20_lo_offset[i++] = ((p7 << 3) + 3);
		dit_p20_lo_offset[i++] = ((p6 << 3) + 4);
		dit_p20_lo_offset[i++] = ((p5 << 3) + 0);
		dit_p20_lo_offset[i++] = ((p4 << 3) + 1);
		dit_p20_lo_offset[i++] = ((p3 << 3) + 2);
		dit_p20_lo_offset[i++] = ((p2 << 3) + 3);
		dit_p20_lo_offset[i++] = ((p1 << 3) + 4);

	// dit_offsets are w.r.to a-array, need 5 distinct sets of these, one for each radix-32 DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		// Set 0: [0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00],[f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10]
		i = 0;
		dit_offsets[0x00] =  0;		dit_offsets[0x10] = pf+p10;
		dit_offsets[0x01] = p1;		dit_offsets[0x11] = pe+p10;
		dit_offsets[0x02] = p3;		dit_offsets[0x12] = pd+p10;
		dit_offsets[0x03] = p2;		dit_offsets[0x13] = pc+p10;
		dit_offsets[0x04] = p7;		dit_offsets[0x14] = pb+p10;
		dit_offsets[0x05] = p6;		dit_offsets[0x15] = pa+p10;
		dit_offsets[0x06] = p5;		dit_offsets[0x16] = p9+p10;
		dit_offsets[0x07] = p4;		dit_offsets[0x17] = p8+p10;
		dit_offsets[0x08] = pf;		dit_offsets[0x18] = p7+p10;
		dit_offsets[0x09] = pe;		dit_offsets[0x19] = p6+p10;
		dit_offsets[0x0a] = pd;		dit_offsets[0x1a] = p5+p10;
		dit_offsets[0x0b] = pc;		dit_offsets[0x1b] = p4+p10;
		dit_offsets[0x0c] = pb;		dit_offsets[0x1c] = p3+p10;
		dit_offsets[0x0d] = pa;		dit_offsets[0x1d] = p2+p10;
		dit_offsets[0x0e] = p9;		dit_offsets[0x1e] = p1+p10;
		dit_offsets[0x0f] = p8;		dit_offsets[0x1f] =    p10;
		// Set 1: [3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p70],[3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p60] (mod p20):
		i += 32;
		dit_offsets[i+0x00] = p3+p10;		dit_offsets[i+0x10] = p3;
		dit_offsets[i+0x01] = p2+p10;		dit_offsets[i+0x11] = p2;
		dit_offsets[i+0x02] = p1+p10;		dit_offsets[i+0x12] = p1;
		dit_offsets[i+0x03] =    p10;		dit_offsets[i+0x13] =  0;
		dit_offsets[i+0x04] = p5+p10;		dit_offsets[i+0x14] = p5;
		dit_offsets[i+0x05] = p4+p10;		dit_offsets[i+0x15] = p4;
		dit_offsets[i+0x06] = p6+p10;		dit_offsets[i+0x16] = p6;
		dit_offsets[i+0x07] = p7+p10;		dit_offsets[i+0x17] = p7;
		dit_offsets[i+0x08] = pd+p10;		dit_offsets[i+0x18] = pd;
		dit_offsets[i+0x09] = pc+p10;		dit_offsets[i+0x19] = pc;
		dit_offsets[i+0x0a] = pe+p10;		dit_offsets[i+0x1a] = pe;
		dit_offsets[i+0x0b] = pf+p10;		dit_offsets[i+0x1b] = pf;
		dit_offsets[i+0x0c] = p9+p10;		dit_offsets[i+0x1c] = p9;
		dit_offsets[i+0x0d] = p8+p10;		dit_offsets[i+0x1d] = p8;
		dit_offsets[i+0x0e] = pa+p10;		dit_offsets[i+0x1e] = pa;
		dit_offsets[i+0x0f] = pb+p10;		dit_offsets[i+0x1f] = pb;
		// Set 2: [9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p30],[9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p20] (mod p20):
		i += 32;
		dit_offsets[i+0x00] = p9+p10;		dit_offsets[i+0x10] = p9;
		dit_offsets[i+0x01] = p8+p10;		dit_offsets[i+0x11] = p8;
		dit_offsets[i+0x02] = pa+p10;		dit_offsets[i+0x12] = pa;
		dit_offsets[i+0x03] = pb+p10;		dit_offsets[i+0x13] = pb;
		dit_offsets[i+0x04] = pe+p10;		dit_offsets[i+0x14] = pe;
		dit_offsets[i+0x05] = pf+p10;		dit_offsets[i+0x15] = pf;
		dit_offsets[i+0x06] = pc+p10;		dit_offsets[i+0x16] = pc;
		dit_offsets[i+0x07] = pd+p10;		dit_offsets[i+0x17] = pd;
		dit_offsets[i+0x08] = p1+p10;		dit_offsets[i+0x18] = p1;
		dit_offsets[i+0x09] =    p10;		dit_offsets[i+0x19] =  0;
		dit_offsets[i+0x0a] = p2+p10;		dit_offsets[i+0x1a] = p2;
		dit_offsets[i+0x0b] = p3+p10;		dit_offsets[i+0x1b] = p3;
		dit_offsets[i+0x0c] = p6+p10;		dit_offsets[i+0x1c] = p6;
		dit_offsets[i+0x0d] = p7+p10;		dit_offsets[i+0x1d] = p7;
		dit_offsets[i+0x0e] = p4+p10;		dit_offsets[i+0x1e] = p4;
		dit_offsets[i+0x0f] = p5+p10;		dit_offsets[i+0x1f] = p5;
		// Set 3: [6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + p80],[a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p90] (mod p20):
		i += 32;
		dit_offsets[i+0x00] = p6;		dit_offsets[i+0x10] = pa+p10;
		dit_offsets[i+0x01] = p7;		dit_offsets[i+0x11] = pb+p10;
		dit_offsets[i+0x02] = p4;		dit_offsets[i+0x12] = p8+p10;
		dit_offsets[i+0x03] = p5;		dit_offsets[i+0x13] = p9+p10;
		dit_offsets[i+0x04] = p2;		dit_offsets[i+0x14] = pc+p10;
		dit_offsets[i+0x05] = p3;		dit_offsets[i+0x15] = pd+p10;
		dit_offsets[i+0x06] =  0;		dit_offsets[i+0x16] = pf+p10;
		dit_offsets[i+0x07] = p1;		dit_offsets[i+0x17] = pe+p10;
		dit_offsets[i+0x08] = pa;		dit_offsets[i+0x18] = p2+p10;
		dit_offsets[i+0x09] = pb;		dit_offsets[i+0x19] = p3+p10;
		dit_offsets[i+0x0a] = p8;		dit_offsets[i+0x1a] =    p10;
		dit_offsets[i+0x0b] = p9;		dit_offsets[i+0x1b] = p1+p10;
		dit_offsets[i+0x0c] = pc;		dit_offsets[i+0x1c] = p4+p10;
		dit_offsets[i+0x0d] = pd;		dit_offsets[i+0x1d] = p5+p10;
		dit_offsets[i+0x0e] = pf;		dit_offsets[i+0x1e] = p7+p10;
		dit_offsets[i+0x0f] = pe;		dit_offsets[i+0x1f] = p6+p10;
		// Set 4: [c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p40],[4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + p50]] (mod p20):
		i += 32;
		dit_offsets[i+0x00] = pc;		dit_offsets[i+0x10] = p4+p10;
		dit_offsets[i+0x01] = pd;		dit_offsets[i+0x11] = p5+p10;
		dit_offsets[i+0x02] = pf;		dit_offsets[i+0x12] = p7+p10;
		dit_offsets[i+0x03] = pe;		dit_offsets[i+0x13] = p6+p10;
		dit_offsets[i+0x04] = p8;		dit_offsets[i+0x14] =    p10;
		dit_offsets[i+0x05] = p9;		dit_offsets[i+0x15] = p1+p10;
		dit_offsets[i+0x06] = pb;		dit_offsets[i+0x16] = p3+p10;
		dit_offsets[i+0x07] = pa;		dit_offsets[i+0x17] = p2+p10;
		dit_offsets[i+0x08] = p4;		dit_offsets[i+0x18] = p8+p10;
		dit_offsets[i+0x09] = p5;		dit_offsets[i+0x19] = p9+p10;
		dit_offsets[i+0x0a] = p7;		dit_offsets[i+0x1a] = pb+p10;
		dit_offsets[i+0x0b] = p6;		dit_offsets[i+0x1b] = pa+p10;
		dit_offsets[i+0x0c] =  0;		dit_offsets[i+0x1c] = pf+p10;
		dit_offsets[i+0x0d] = p1;		dit_offsets[i+0x1d] = pe+p10;
		dit_offsets[i+0x0e] = p3;		dit_offsets[i+0x1e] = pd+p10;
		dit_offsets[i+0x0f] = p2;		dit_offsets[i+0x1f] = pc+p10;
	}

/*...The radix-160 pass is here.	*/

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
	store the 5 index-offset 32-tets going into the radix-32 DFTs in the dit_offset array.

	(***NOTE*** The following set of permutations can be auto-generated by running test_dft_radix()
	for the radix in question in TEST_TYPE = 1 [DIT] mode, skipping the actual DIT-pass step initially):
	Combined DIT input-scramble array:

	Combined DIT input-scramble array = [
		[0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00],[f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10]
		[3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p70],[3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p60]
		[9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p30],[9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p20]
		[6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + p80],[a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p90]
		[c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p40],[4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + p50]
	*/
	//...gather the needed data (160 64-bit complex) and do 5 radix-32 transforms:
		tptr = t;
		jt = j1    ; RADIX_32_DIT((a+jt),dit_offsets     ,RE_IM_STRIDE, (double *)tptr,t_offsets,1);	tptr += 32;
		jt = j1+p60; RADIX_32_DIT((a+jt),dit_offsets+0x20,RE_IM_STRIDE, (double *)tptr,t_offsets,1);	tptr += 32;
		jt = j1+p20; RADIX_32_DIT((a+jt),dit_offsets+0x40,RE_IM_STRIDE, (double *)tptr,t_offsets,1);	tptr += 32;
		jt = j1+p80; RADIX_32_DIT((a+jt),dit_offsets+0x60,RE_IM_STRIDE, (double *)tptr,t_offsets,1);	tptr += 32;
		jt = j1+p40; RADIX_32_DIT((a+jt),dit_offsets+0x80,RE_IM_STRIDE, (double *)tptr,t_offsets,1);

	/*...and now do 32 radix-5 transforms, with the columns of t*[r,i] output pairs in the above 5x radix-32 set now acting as input rows.
	Since our first-look oindex ordering was +p0x[0,20,40,60,80] for each radix-5 and incrementing += p1 between those DFTs,
	arrange resulting mismatched-data-sorted index permutation into 5 vertical 32-entry columns to get needed oindex patterns.
	Indexing in hex for clarity and using [evn|odd]0-4 notation in the rightmost column to flag reusable 5-perms
	[in fact simple circular (0-4)-element shifts of the basic patterns evn := 00,60,20,80,40 and odd := 50,10,70,30,90:

		00,60,20,80,40		00,60,20,80,40 + p0		[evn0] + p0
		5f,1f,7f,3f,9f		50,10,70,30,90 + pf		[odd0] + pf
		1e,7e,3e,9e,5e		10,70,30,90,50 + pe		[odd1] + pe
		7d,3d,9d,5d,1d		70,30,90,50,10 + pd		[odd2] + pd
		3c,9c,5c,1c,7c		30,90,50,10,70 + pc		[odd3] + pc
		9b,5b,1b,7b,3b		90,50,10,70,30 + pb		[odd4] + pb
		5a,1a,7a,3a,9a		50,10,70,30,90 + pa		[odd0] + pa
		19,79,39,99,59		10,70,30,90,50 + p9		[odd1] + p9
		78,38,98,58,18		70,30,90,50,10 + p8		[odd2] + p8
		37,97,57,17,77		30,90,50,10,70 + p7		[odd3] + p7
		96,56,16,76,36		90,50,10,70,30 + p6		[odd4] + p6
		55,15,75,35,95		50,10,70,30,90 + p5		[odd0] + p5
		14,74,34,94,54		10,70,30,90,50 + p4		[odd1] + p4
		73,33,93,53,13		70,30,90,50,10 + p3		[odd2] + p3
		32,92,52,12,72		30,90,50,10,70 + p2		[odd3] + p2
		91,51,11,71,31		90,50,10,70,30 + p1		[odd4] + p1
		50,10,70,30,90	=	50,10,70,30,90 + p0	=	[odd0] + p0
		0f,6f,2f,8f,4f		00,60,20,80,40 + pf		[evn0] + pf
		6e,2e,8e,4e,0e		60,20,80,40,00 + pe		[evn1] + pe
		2d,8d,4d,0d,6d		20,80,40,00,60 + pd		[evn2] + pd
		8c,4c,0c,6c,2c		80,40,00,60,20 + pc		[evn3] + pc
		4b,0b,6b,2b,8b		40,00,60,20,80 + pb		[evn4] + pb
		0a,6a,2a,8a,4a		00,60,20,80,40 + pa		[evn0] + pa
		69,29,89,49,09		60,20,80,40,00 + p9		[evn1] + p9
		28,88,48,08,68		20,80,40,00,60 + p8		[evn2] + p8
		87,47,07,67,27		80,40,00,60,20 + p7		[evn3] + p7
		46,06,66,26,86		40,00,60,20,80 + p6		[evn4] + p6
		05,65,25,85,45		00,60,20,80,40 + p5		[evn0] + p5
		64,24,84,44,04		60,20,80,40,00 + p4		[evn1] + p4
		23,83,43,03,63		20,80,40,00,60 + p3		[evn2] + p3
		82,42,02,62,22		80,40,00,60,20 + p2		[evn3] + p2
		41,01,61,21,81		40,00,60,20,80 + p1		[evn4] + p1
	*/
		tptr = t;
		for(i = 0; i < 32; i++) {
			int k = dit_p20_lo_offset[i];
			// Extract index (in [0-4]) into circ-shift array used for high parts of p-mults. The [0-4] value is
			// in low 3 bits of k; the "which length-9 half of the dit_p20_cperms array?" selector is via (k < 0):
			int ic = ((-(k < 0)) & 9)	// +/- sign on k puts us into lower/upper half of the cshift array (base index 0/9)
						+ (k & 0x7);	// ...and low 3 bits give the element index w.r.to the array-half in question.
			int k0 = dit_p20_cperms[ic], k1 = dit_p20_cperms[ic+1], k2 = dit_p20_cperms[ic+2], k3 = dit_p20_cperms[ic+3], k4 = dit_p20_cperms[ic+4];
			// Extract Low part, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-5 DFTs:
			k = (k & 0x7fffffff) >> 3;
			jt = j1+k; jp = j2+k;
			RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,
				tptr->re,tptr->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x80)->re,(tptr+0x80)->im,
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
	cy160_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr;
		const int pfetch_dist = PFETCH_DIST;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
					,p10,p20,p30,p40,p50,p60,p70,p80,p90;
	// Shared DIF+DIT:
		double rt,it;
		double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
		int poff[RADIX>>2];	// Store [RADIX/4] mults of p04 offset for loop control
		int t_offsets[32];
		// Need storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 2*9 elts:
		int dif_offsets[RADIX], dif_p20_cperms[18], dif_p20_lo_offset[32], dif_phi[ODD_RADIX];
		int dit_offsets[RADIX], dit_p20_cperms[18], dit_p20_lo_offset[32], dit_phi[ODD_RADIX];

		int incr,j,j1,j2,jt,jp,k,l,ntmp;
		// incr = Carry-chain wts-multipliers recurrence length, which must divide
		// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 10|20|40 for avx512,avx,sse, respectively:
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
		int *itmp,*itm2;	// Pointer into the bjmodn array
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		double *add0, *add1, *add2, *add3;
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *tmp,*tm1,*tm2,	// Non-static utility ptrs
			*va0,*va1,*va2,*va3,*va4, *vb0,*vb1,*vb2,*vb3,*vb4,
			*wa0,*wa1,*wa2,*wa3,*wa4, *wb0,*wb1,*wb2,*wb3,*wb4;
		vec_dbl *two,*one,*sqrt2,*isrt2,*xcc1,*xss1,*xcc2,*xss2,*xcc3,*xss3,	// radix-32 DFT trig consts
			*ycc1,*yss1,*ycc2,*yss2,*yss3,	// radiy-5 DFT trig consts
			*max_err, *sse2_rnd, *half_arr,
			*r00,	// Head of RADIX*vec_cmplx-sized local store #1
			*s1p00,	// Head of RADIX*vec_cmplx-sized local store #2
			*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx
		double dtmp;
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
		int thr_id = thread_arg->tid;
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
		NDIVR >>= 4;			p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
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

	// Shared:
		// Set array offsets for radix-32 DFT in/outputs:
		// t_offsets w.r.to: t-array, same for all 5 DFTs:
		t_offsets[0x00] = 0x00<<1;	t_offsets[0x10] = 0x10<<1;
		t_offsets[0x01] = 0x01<<1;	t_offsets[0x11] = 0x11<<1;
		t_offsets[0x02] = 0x02<<1;	t_offsets[0x12] = 0x12<<1;
		t_offsets[0x03] = 0x03<<1;	t_offsets[0x13] = 0x13<<1;
		t_offsets[0x04] = 0x04<<1;	t_offsets[0x14] = 0x14<<1;
		t_offsets[0x05] = 0x05<<1;	t_offsets[0x15] = 0x15<<1;
		t_offsets[0x06] = 0x06<<1;	t_offsets[0x16] = 0x16<<1;
		t_offsets[0x07] = 0x07<<1;	t_offsets[0x17] = 0x17<<1;
		t_offsets[0x08] = 0x08<<1;	t_offsets[0x18] = 0x18<<1;
		t_offsets[0x09] = 0x09<<1;	t_offsets[0x19] = 0x19<<1;
		t_offsets[0x0a] = 0x0a<<1;	t_offsets[0x1a] = 0x1a<<1;
		t_offsets[0x0b] = 0x0b<<1;	t_offsets[0x1b] = 0x1b<<1;
		t_offsets[0x0c] = 0x0c<<1;	t_offsets[0x1c] = 0x1c<<1;
		t_offsets[0x0d] = 0x0d<<1;	t_offsets[0x1d] = 0x1d<<1;
		t_offsets[0x0e] = 0x0e<<1;	t_offsets[0x1e] = 0x1e<<1;
		t_offsets[0x0f] = 0x0f<<1;	t_offsets[0x1f] = 0x1f<<1;

	/*** DIF indexing stuff: ***/

		dif_phi[0] =   0;		dit_phi[0] =   0;
		dif_phi[1] = p80;		dit_phi[1] = p60;
		dif_phi[2] = p60;		dit_phi[2] = p20;
		dif_phi[3] = p40;		dit_phi[3] = p80;
		dif_phi[4] = p20;		dit_phi[4] = p40;

		// Init storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 2*9

	   #ifdef USE_SSE2
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p10*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:
		l = 0;
		dif_p20_cperms[l++] = 0x00<<1; dif_p20_cperms[l++] = 0x80<<1; dif_p20_cperms[l++] = 0x60<<1; dif_p20_cperms[l++] = 0x40<<1; dif_p20_cperms[l++] = 0x20<<1; dif_p20_cperms[l++] = 0x00<<1; dif_p20_cperms[l++] = 0x80<<1; dif_p20_cperms[l++] = 0x60<<1; dif_p20_cperms[l++] = 0x40<<1;
		dif_p20_cperms[l++] = 0x90<<1; dif_p20_cperms[l++] = 0x70<<1; dif_p20_cperms[l++] = 0x50<<1; dif_p20_cperms[l++] = 0x30<<1; dif_p20_cperms[l++] = 0x10<<1; dif_p20_cperms[l++] = 0x90<<1; dif_p20_cperms[l++] = 0x70<<1; dif_p20_cperms[l++] = 0x50<<1; dif_p20_cperms[l++] = 0x30<<1;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-5 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-$ index, bits <4:30> for the p0-f. In SIMD mode again replace p[0-f] with 0x[0-f]<<1 for use with
		// contig-local-mem, thus replace 'p' prefixes with 0x and 3-left-shifts with << 4 to account for the << 1:
		l = 0;
		dif_p20_lo_offset[l++] = ((0x0 << 4) + 0);
		dif_p20_lo_offset[l++] = ((0xb << 4) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x6 << 4) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x1 << 4) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xc << 4) + 1);
		dif_p20_lo_offset[l++] = ((0x7 << 4) + 1);
		dif_p20_lo_offset[l++] = ((0x2 << 4) + 1);
		dif_p20_lo_offset[l++] = ((0xd << 4) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x8 << 4) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x3 << 4) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xe << 4) + 2);
		dif_p20_lo_offset[l++] = ((0x9 << 4) + 2);
		dif_p20_lo_offset[l++] = ((0x4 << 4) + 2);
		dif_p20_lo_offset[l++] = ((0xf << 4) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xa << 4) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x5 << 4) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x0 << 4) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xb << 4) + 3);
		dif_p20_lo_offset[l++] = ((0x6 << 4) + 3);
		dif_p20_lo_offset[l++] = ((0x1 << 4) + 3);
		dif_p20_lo_offset[l++] = ((0xc << 4) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x7 << 4) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x2 << 4) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xd << 4) + 4);
		dif_p20_lo_offset[l++] = ((0x8 << 4) + 4);
		dif_p20_lo_offset[l++] = ((0x3 << 4) + 4);
		dif_p20_lo_offset[l++] = ((0xe << 4) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x9 << 4) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x4 << 4) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xf << 4) + 0);
		dif_p20_lo_offset[l++] = ((0xa << 4) + 0);
		dif_p20_lo_offset[l++] = ((0x5 << 4) + 0);

	   #else

		// Init storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 2*9
		l = 0;
		dif_p20_cperms[l++] =   0; dif_p20_cperms[l++] = p80; dif_p20_cperms[l++] = p60; dif_p20_cperms[l++] = p40; dif_p20_cperms[l++] = p20; dif_p20_cperms[l++] =   0; dif_p20_cperms[l++] = p80; dif_p20_cperms[l++] = p60; dif_p20_cperms[l++] = p40;
		dif_p20_cperms[l++] = p90; dif_p20_cperms[l++] = p70; dif_p20_cperms[l++] = p50; dif_p20_cperms[l++] = p30; dif_p20_cperms[l++] = p10; dif_p20_cperms[l++] = p90; dif_p20_cperms[l++] = p70; dif_p20_cperms[l++] = p50; dif_p20_cperms[l++] = p30;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-5 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-4 index, bits <4:30> for the p0-f:
		l = 0;
		dif_p20_lo_offset[l++] = (( 0 << 3) + 0);
		dif_p20_lo_offset[l++] = ((pb << 3) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p6 << 3) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p1 << 3) + 0) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pc << 3) + 1);
		dif_p20_lo_offset[l++] = ((p7 << 3) + 1);
		dif_p20_lo_offset[l++] = ((p2 << 3) + 1);
		dif_p20_lo_offset[l++] = ((pd << 3) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p8 << 3) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p3 << 3) + 1) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pe << 3) + 2);
		dif_p20_lo_offset[l++] = ((p9 << 3) + 2);
		dif_p20_lo_offset[l++] = ((p4 << 3) + 2);
		dif_p20_lo_offset[l++] = ((pf << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pa << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p5 << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = (( 0 << 3) + 2) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pb << 3) + 3);
		dif_p20_lo_offset[l++] = ((p6 << 3) + 3);
		dif_p20_lo_offset[l++] = ((p1 << 3) + 3);
		dif_p20_lo_offset[l++] = ((pc << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p7 << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p2 << 3) + 3) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pd << 3) + 4);
		dif_p20_lo_offset[l++] = ((p8 << 3) + 4);
		dif_p20_lo_offset[l++] = ((p3 << 3) + 4);
		dif_p20_lo_offset[l++] = ((pe << 3) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p9 << 3) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p4 << 3) + 4) + ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pf << 3) + 0);
		dif_p20_lo_offset[l++] = ((pa << 3) + 0);
		dif_p20_lo_offset[l++] = ((p5 << 3) + 0);

	   #endif	// sse2?

	// dif_offsets are w.r.to a-array, need 7 distinct sets of these, one for each DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		// Set 0: [0,1,3,2,6,7,4,5,c,d,f,e,9,8,a,b + p00],[9,8,a,b,f,e,d,c,3,2,1,0,4,5,7,6 + p10]
		l = 0;
		dif_offsets[0x00] =  0;		dif_offsets[0x10] = p9+p10;
		dif_offsets[0x01] = p1;		dif_offsets[0x11] = p8+p10;
		dif_offsets[0x02] = p3;		dif_offsets[0x12] = pa+p10;
		dif_offsets[0x03] = p2;		dif_offsets[0x13] = pb+p10;
		dif_offsets[0x04] = p6;		dif_offsets[0x14] = pf+p10;
		dif_offsets[0x05] = p7;		dif_offsets[0x15] = pe+p10;
		dif_offsets[0x06] = p4;		dif_offsets[0x16] = pd+p10;
		dif_offsets[0x07] = p5;		dif_offsets[0x17] = pc+p10;
		dif_offsets[0x08] = pc;		dif_offsets[0x18] = p3+p10;
		dif_offsets[0x09] = pd;		dif_offsets[0x19] = p2+p10;
		dif_offsets[0x0a] = pf;		dif_offsets[0x1a] = p1+p10;
		dif_offsets[0x0b] = pe;		dif_offsets[0x1b] =    p10;
		dif_offsets[0x0c] = p9;		dif_offsets[0x1c] = p4+p10;
		dif_offsets[0x0d] = p8;		dif_offsets[0x1d] = p5+p10;
		dif_offsets[0x0e] = pa;		dif_offsets[0x1e] = p7+p10;
		dif_offsets[0x0f] = pb;		dif_offsets[0x1f] = p6+p10;
		// Set 1: [6,7,4,5,3,2,1,0,9,8,a,b,f,e,d,c + p80],[f,e,d,c,a,b,8,9,4,5,7,6,1,0,2,3 + p90]
		l += 32;
		dif_offsets[l+0x00] = p6;		dif_offsets[l+0x10] = pf+p10;
		dif_offsets[l+0x01] = p7;		dif_offsets[l+0x11] = pe+p10;
		dif_offsets[l+0x02] = p4;		dif_offsets[l+0x12] = pd+p10;
		dif_offsets[l+0x03] = p5;		dif_offsets[l+0x13] = pc+p10;
		dif_offsets[l+0x04] = p3;		dif_offsets[l+0x14] = pa+p10;
		dif_offsets[l+0x05] = p2;		dif_offsets[l+0x15] = pb+p10;
		dif_offsets[l+0x06] = p1;		dif_offsets[l+0x16] = p8+p10;
		dif_offsets[l+0x07] =  0;		dif_offsets[l+0x17] = p9+p10;
		dif_offsets[l+0x08] = p9;		dif_offsets[l+0x18] = p4+p10;
		dif_offsets[l+0x09] = p8;		dif_offsets[l+0x19] = p5+p10;
		dif_offsets[l+0x0a] = pa;		dif_offsets[l+0x1a] = p7+p10;
		dif_offsets[l+0x0b] = pb;		dif_offsets[l+0x1b] = p6+p10;
		dif_offsets[l+0x0c] = pf;		dif_offsets[l+0x1c] = p1+p10;
		dif_offsets[l+0x0d] = pe;		dif_offsets[l+0x1d] =    p10;
		dif_offsets[l+0x0e] = pd;		dif_offsets[l+0x1e] = p2+p10;
		dif_offsets[l+0x0f] = pc;		dif_offsets[l+0x1f] = p3+p10;
		// Set 2: [3,2,1,0,4,5,7,6,f,e,d,c,a,b,8,9 + p70],[6,7,4,5,3,2,1,0,9,8,a,b,f,e,d,c + p60]
		l += 32;
		dif_offsets[l+0x00] = p3+p10;		dif_offsets[l+0x10] = p6;
		dif_offsets[l+0x01] = p2+p10;		dif_offsets[l+0x11] = p7;
		dif_offsets[l+0x02] = p1+p10;		dif_offsets[l+0x12] = p4;
		dif_offsets[l+0x03] =    p10;		dif_offsets[l+0x13] = p5;
		dif_offsets[l+0x04] = p4+p10;		dif_offsets[l+0x14] = p3;
		dif_offsets[l+0x05] = p5+p10;		dif_offsets[l+0x15] = p2;
		dif_offsets[l+0x06] = p7+p10;		dif_offsets[l+0x16] = p1;
		dif_offsets[l+0x07] = p6+p10;		dif_offsets[l+0x17] =  0;
		dif_offsets[l+0x08] = pf+p10;		dif_offsets[l+0x18] = p9;
		dif_offsets[l+0x09] = pe+p10;		dif_offsets[l+0x19] = p8;
		dif_offsets[l+0x0a] = pd+p10;		dif_offsets[l+0x1a] = pa;
		dif_offsets[l+0x0b] = pc+p10;		dif_offsets[l+0x1b] = pb;
		dif_offsets[l+0x0c] = pa+p10;		dif_offsets[l+0x1c] = pf;
		dif_offsets[l+0x0d] = pb+p10;		dif_offsets[l+0x1d] = pe;
		dif_offsets[l+0x0e] = p8+p10;		dif_offsets[l+0x1e] = pd;
		dif_offsets[l+0x0f] = p9+p10;		dif_offsets[l+0x1f] = pc;
		// Set 3: [c,d,f,e,9,8,a,b,6,7,4,5,3,2,1,0 + p40],[3,2,1,0,4,5,7,6,f,e,d,c,a,b,8,9 + p50]
		l += 32;
		dif_offsets[l+0x00] = pc;		dif_offsets[l+0x10] = p3+p10;
		dif_offsets[l+0x01] = pd;		dif_offsets[l+0x11] = p2+p10;
		dif_offsets[l+0x02] = pf;		dif_offsets[l+0x12] = p1+p10;
		dif_offsets[l+0x03] = pe;		dif_offsets[l+0x13] =    p10;
		dif_offsets[l+0x04] = p9;		dif_offsets[l+0x14] = p4+p10;
		dif_offsets[l+0x05] = p8;		dif_offsets[l+0x15] = p5+p10;
		dif_offsets[l+0x06] = pa;		dif_offsets[l+0x16] = p7+p10;
		dif_offsets[l+0x07] = pb;		dif_offsets[l+0x17] = p6+p10;
		dif_offsets[l+0x08] = p6;		dif_offsets[l+0x18] = pf+p10;
		dif_offsets[l+0x09] = p7;		dif_offsets[l+0x19] = pe+p10;
		dif_offsets[l+0x0a] = p4;		dif_offsets[l+0x1a] = pd+p10;
		dif_offsets[l+0x0b] = p5;		dif_offsets[l+0x1b] = pc+p10;
		dif_offsets[l+0x0c] = p3;		dif_offsets[l+0x1c] = pa+p10;
		dif_offsets[l+0x0d] = p2;		dif_offsets[l+0x1d] = pb+p10;
		dif_offsets[l+0x0e] = p1;		dif_offsets[l+0x1e] = p8+p10;
		dif_offsets[l+0x0f] =  0;		dif_offsets[l+0x1f] = p9+p10;
		// Set 4: [9,8,a,b,f,e,d,c,3,2,1,0,4,5,7,6 + p30],[c,d,f,e,9,8,a,b,6,7,4,5,3,2,1,0 + p20]
		l += 32;
		dif_offsets[l+0x00] = p9+p10;		dif_offsets[l+0x10] = pc;
		dif_offsets[l+0x01] = p8+p10;		dif_offsets[l+0x11] = pd;
		dif_offsets[l+0x02] = pa+p10;		dif_offsets[l+0x12] = pf;
		dif_offsets[l+0x03] = pb+p10;		dif_offsets[l+0x13] = pe;
		dif_offsets[l+0x04] = pf+p10;		dif_offsets[l+0x14] = p9;
		dif_offsets[l+0x05] = pe+p10;		dif_offsets[l+0x15] = p8;
		dif_offsets[l+0x06] = pd+p10;		dif_offsets[l+0x16] = pa;
		dif_offsets[l+0x07] = pc+p10;		dif_offsets[l+0x17] = pb;
		dif_offsets[l+0x08] = p3+p10;		dif_offsets[l+0x18] = p6;
		dif_offsets[l+0x09] = p2+p10;		dif_offsets[l+0x19] = p7;
		dif_offsets[l+0x0a] = p1+p10;		dif_offsets[l+0x1a] = p4;
		dif_offsets[l+0x0b] =    p10;		dif_offsets[l+0x1b] = p5;
		dif_offsets[l+0x0c] = p4+p10;		dif_offsets[l+0x1c] = p3;
		dif_offsets[l+0x0d] = p5+p10;		dif_offsets[l+0x1d] = p2;
		dif_offsets[l+0x0e] = p7+p10;		dif_offsets[l+0x1e] = p1;
		dif_offsets[l+0x0f] = p6+p10;		dif_offsets[l+0x1f] =  0;
	  #ifdef USE_SSE2
		// IN SIMD mode preshift all the above offsets << 3 to turn into double=array pointer offsets:
		for(l = 0; l < RADIX; l++) {
			dif_offsets[l] <<= 3;
		}
	  #endif

	/*** DIT indexing stuff: ***/
		// Init storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 2*9

	   #ifdef USE_SSE2
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p10*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:
		l = 0;
		dit_p20_cperms[l++] = 0x00<<1; dit_p20_cperms[l++] = 0x60<<1; dit_p20_cperms[l++] = 0x20<<1; dit_p20_cperms[l++] = 0x80<<1; dit_p20_cperms[l++] = 0x40<<1; dit_p20_cperms[l++] = 0x00<<1; dit_p20_cperms[l++] = 0x60<<1; dit_p20_cperms[l++] = 0x20<<1; dit_p20_cperms[l++] = 0x80<<1;
		dit_p20_cperms[l++] = 0x50<<1; dit_p20_cperms[l++] = 0x10<<1; dit_p20_cperms[l++] = 0x70<<1; dit_p20_cperms[l++] = 0x30<<1; dit_p20_cperms[l++] = 0x90<<1; dit_p20_cperms[l++] = 0x50<<1; dit_p20_cperms[l++] = 0x10<<1; dit_p20_cperms[l++] = 0x70<<1; dit_p20_cperms[l++] = 0x30<<1;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-5 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-4 index, bits <4:30> for the p0-f. In SIMD mode again replace p[0-f] with 0x[0-f]<<1 for use with
		// contig-local-mem, thus replace 'p' prefixes with 0x and 3-left-shifts with << 4 to account for the << 1:
		l = 0;
		dit_p20_lo_offset[l++] = ((0x0 << 4) + 0);
		dit_p20_lo_offset[l++] = ((0xf << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xe << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xd << 4) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xc << 4) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xb << 4) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xa << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x9 << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x8 << 4) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x7 << 4) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x6 << 4) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x5 << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x4 << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x3 << 4) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x2 << 4) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x1 << 4) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x0 << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xf << 4) + 0);
		dit_p20_lo_offset[l++] = ((0xe << 4) + 1);
		dit_p20_lo_offset[l++] = ((0xd << 4) + 2);
		dit_p20_lo_offset[l++] = ((0xc << 4) + 3);
		dit_p20_lo_offset[l++] = ((0xb << 4) + 4);
		dit_p20_lo_offset[l++] = ((0xa << 4) + 0);
		dit_p20_lo_offset[l++] = ((0x9 << 4) + 1);
		dit_p20_lo_offset[l++] = ((0x8 << 4) + 2);
		dit_p20_lo_offset[l++] = ((0x7 << 4) + 3);
		dit_p20_lo_offset[l++] = ((0x6 << 4) + 4);
		dit_p20_lo_offset[l++] = ((0x5 << 4) + 0);
		dit_p20_lo_offset[l++] = ((0x4 << 4) + 1);
		dit_p20_lo_offset[l++] = ((0x3 << 4) + 2);
		dit_p20_lo_offset[l++] = ((0x2 << 4) + 3);
		dit_p20_lo_offset[l++] = ((0x1 << 4) + 4);

	   #else

		// Init storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,4] that means 2*9
		l = 0;
		dit_p20_cperms[l++] =   0; dit_p20_cperms[l++] = p60; dit_p20_cperms[l++] = p20; dit_p20_cperms[l++] = p80; dit_p20_cperms[l++] = p40; dit_p20_cperms[l++] =   0; dit_p20_cperms[l++] = p60; dit_p20_cperms[l++] = p20; dit_p20_cperms[l++] = p80;
		dit_p20_cperms[l++] = p50; dit_p20_cperms[l++] = p10; dit_p20_cperms[l++] = p70; dit_p20_cperms[l++] = p30; dit_p20_cperms[l++] = p90; dit_p20_cperms[l++] = p50; dit_p20_cperms[l++] = p10; dit_p20_cperms[l++] = p70; dit_p20_cperms[l++] = p30;
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-5 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 3 bits for
		// the 0-4 index, bits <4:30> for the p0-f:
		l = 0;
		dit_p20_lo_offset[l++] = (( 0 << 3) + 0);
		dit_p20_lo_offset[l++] = ((pf << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pe << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pd << 3) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pc << 3) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pb << 3) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pa << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p9 << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p8 << 3) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p7 << 3) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p6 << 3) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p5 << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p4 << 3) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p3 << 3) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p2 << 3) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p1 << 3) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = (( 0 << 3) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pf << 3) + 0);
		dit_p20_lo_offset[l++] = ((pe << 3) + 1);
		dit_p20_lo_offset[l++] = ((pd << 3) + 2);
		dit_p20_lo_offset[l++] = ((pc << 3) + 3);
		dit_p20_lo_offset[l++] = ((pb << 3) + 4);
		dit_p20_lo_offset[l++] = ((pa << 3) + 0);
		dit_p20_lo_offset[l++] = ((p9 << 3) + 1);
		dit_p20_lo_offset[l++] = ((p8 << 3) + 2);
		dit_p20_lo_offset[l++] = ((p7 << 3) + 3);
		dit_p20_lo_offset[l++] = ((p6 << 3) + 4);
		dit_p20_lo_offset[l++] = ((p5 << 3) + 0);
		dit_p20_lo_offset[l++] = ((p4 << 3) + 1);
		dit_p20_lo_offset[l++] = ((p3 << 3) + 2);
		dit_p20_lo_offset[l++] = ((p2 << 3) + 3);
		dit_p20_lo_offset[l++] = ((p1 << 3) + 4);

	   #endif	// sse2?

	// dit_offsets are w.r.to a-array, need 5 distinct sets of these, one for each radix-32 DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		// Set 0: [0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00],[f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10]
		l = 0;
		dit_offsets[0x00] =  0;		dit_offsets[0x10] = pf+p10;
		dit_offsets[0x01] = p1;		dit_offsets[0x11] = pe+p10;
		dit_offsets[0x02] = p3;		dit_offsets[0x12] = pd+p10;
		dit_offsets[0x03] = p2;		dit_offsets[0x13] = pc+p10;
		dit_offsets[0x04] = p7;		dit_offsets[0x14] = pb+p10;
		dit_offsets[0x05] = p6;		dit_offsets[0x15] = pa+p10;
		dit_offsets[0x06] = p5;		dit_offsets[0x16] = p9+p10;
		dit_offsets[0x07] = p4;		dit_offsets[0x17] = p8+p10;
		dit_offsets[0x08] = pf;		dit_offsets[0x18] = p7+p10;
		dit_offsets[0x09] = pe;		dit_offsets[0x19] = p6+p10;
		dit_offsets[0x0a] = pd;		dit_offsets[0x1a] = p5+p10;
		dit_offsets[0x0b] = pc;		dit_offsets[0x1b] = p4+p10;
		dit_offsets[0x0c] = pb;		dit_offsets[0x1c] = p3+p10;
		dit_offsets[0x0d] = pa;		dit_offsets[0x1d] = p2+p10;
		dit_offsets[0x0e] = p9;		dit_offsets[0x1e] = p1+p10;
		dit_offsets[0x0f] = p8;		dit_offsets[0x1f] =    p10;
		// Set 1: [3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p70],[3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p60] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = p3+p10;		dit_offsets[l+0x10] = p3;
		dit_offsets[l+0x01] = p2+p10;		dit_offsets[l+0x11] = p2;
		dit_offsets[l+0x02] = p1+p10;		dit_offsets[l+0x12] = p1;
		dit_offsets[l+0x03] =    p10;		dit_offsets[l+0x13] =  0;
		dit_offsets[l+0x04] = p5+p10;		dit_offsets[l+0x14] = p5;
		dit_offsets[l+0x05] = p4+p10;		dit_offsets[l+0x15] = p4;
		dit_offsets[l+0x06] = p6+p10;		dit_offsets[l+0x16] = p6;
		dit_offsets[l+0x07] = p7+p10;		dit_offsets[l+0x17] = p7;
		dit_offsets[l+0x08] = pd+p10;		dit_offsets[l+0x18] = pd;
		dit_offsets[l+0x09] = pc+p10;		dit_offsets[l+0x19] = pc;
		dit_offsets[l+0x0a] = pe+p10;		dit_offsets[l+0x1a] = pe;
		dit_offsets[l+0x0b] = pf+p10;		dit_offsets[l+0x1b] = pf;
		dit_offsets[l+0x0c] = p9+p10;		dit_offsets[l+0x1c] = p9;
		dit_offsets[l+0x0d] = p8+p10;		dit_offsets[l+0x1d] = p8;
		dit_offsets[l+0x0e] = pa+p10;		dit_offsets[l+0x1e] = pa;
		dit_offsets[l+0x0f] = pb+p10;		dit_offsets[l+0x1f] = pb;
		// Set 2: [9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p30],[9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p20] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = p9+p10;		dit_offsets[l+0x10] = p9;
		dit_offsets[l+0x01] = p8+p10;		dit_offsets[l+0x11] = p8;
		dit_offsets[l+0x02] = pa+p10;		dit_offsets[l+0x12] = pa;
		dit_offsets[l+0x03] = pb+p10;		dit_offsets[l+0x13] = pb;
		dit_offsets[l+0x04] = pe+p10;		dit_offsets[l+0x14] = pe;
		dit_offsets[l+0x05] = pf+p10;		dit_offsets[l+0x15] = pf;
		dit_offsets[l+0x06] = pc+p10;		dit_offsets[l+0x16] = pc;
		dit_offsets[l+0x07] = pd+p10;		dit_offsets[l+0x17] = pd;
		dit_offsets[l+0x08] = p1+p10;		dit_offsets[l+0x18] = p1;
		dit_offsets[l+0x09] =    p10;		dit_offsets[l+0x19] =  0;
		dit_offsets[l+0x0a] = p2+p10;		dit_offsets[l+0x1a] = p2;
		dit_offsets[l+0x0b] = p3+p10;		dit_offsets[l+0x1b] = p3;
		dit_offsets[l+0x0c] = p6+p10;		dit_offsets[l+0x1c] = p6;
		dit_offsets[l+0x0d] = p7+p10;		dit_offsets[l+0x1d] = p7;
		dit_offsets[l+0x0e] = p4+p10;		dit_offsets[l+0x1e] = p4;
		dit_offsets[l+0x0f] = p5+p10;		dit_offsets[l+0x1f] = p5;
		// Set 3: [6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + p80],[a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p90] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = p6;		dit_offsets[l+0x10] = pa+p10;
		dit_offsets[l+0x01] = p7;		dit_offsets[l+0x11] = pb+p10;
		dit_offsets[l+0x02] = p4;		dit_offsets[l+0x12] = p8+p10;
		dit_offsets[l+0x03] = p5;		dit_offsets[l+0x13] = p9+p10;
		dit_offsets[l+0x04] = p2;		dit_offsets[l+0x14] = pc+p10;
		dit_offsets[l+0x05] = p3;		dit_offsets[l+0x15] = pd+p10;
		dit_offsets[l+0x06] =  0;		dit_offsets[l+0x16] = pf+p10;
		dit_offsets[l+0x07] = p1;		dit_offsets[l+0x17] = pe+p10;
		dit_offsets[l+0x08] = pa;		dit_offsets[l+0x18] = p2+p10;
		dit_offsets[l+0x09] = pb;		dit_offsets[l+0x19] = p3+p10;
		dit_offsets[l+0x0a] = p8;		dit_offsets[l+0x1a] =    p10;
		dit_offsets[l+0x0b] = p9;		dit_offsets[l+0x1b] = p1+p10;
		dit_offsets[l+0x0c] = pc;		dit_offsets[l+0x1c] = p4+p10;
		dit_offsets[l+0x0d] = pd;		dit_offsets[l+0x1d] = p5+p10;
		dit_offsets[l+0x0e] = pf;		dit_offsets[l+0x1e] = p7+p10;
		dit_offsets[l+0x0f] = pe;		dit_offsets[l+0x1f] = p6+p10;
		// Set 4: [c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p40],[4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + p50]] (mod p20):
		l += 32;
		dit_offsets[l+0x00] = pc;		dit_offsets[l+0x10] = p4+p10;
		dit_offsets[l+0x01] = pd;		dit_offsets[l+0x11] = p5+p10;
		dit_offsets[l+0x02] = pf;		dit_offsets[l+0x12] = p7+p10;
		dit_offsets[l+0x03] = pe;		dit_offsets[l+0x13] = p6+p10;
		dit_offsets[l+0x04] = p8;		dit_offsets[l+0x14] =    p10;
		dit_offsets[l+0x05] = p9;		dit_offsets[l+0x15] = p1+p10;
		dit_offsets[l+0x06] = pb;		dit_offsets[l+0x16] = p3+p10;
		dit_offsets[l+0x07] = pa;		dit_offsets[l+0x17] = p2+p10;
		dit_offsets[l+0x08] = p4;		dit_offsets[l+0x18] = p8+p10;
		dit_offsets[l+0x09] = p5;		dit_offsets[l+0x19] = p9+p10;
		dit_offsets[l+0x0a] = p7;		dit_offsets[l+0x1a] = pb+p10;
		dit_offsets[l+0x0b] = p6;		dit_offsets[l+0x1b] = pa+p10;
		dit_offsets[l+0x0c] =  0;		dit_offsets[l+0x1c] = pf+p10;
		dit_offsets[l+0x0d] = p1;		dit_offsets[l+0x1d] = pe+p10;
		dit_offsets[l+0x0e] = p3;		dit_offsets[l+0x1e] = pd+p10;
		dit_offsets[l+0x0f] = p2;		dit_offsets[l+0x1f] = pc+p10;
	  #ifdef USE_SSE2
		// IN SIMD mode preshift all the above offsets << 3 to turn into double=array pointer offsets:
		for(l = 0; l < RADIX; l++) {
			dit_offsets[l] <<= 3;
		}
	  #endif

	#ifdef USE_SSE2
		tmp	= r00 = thread_arg->r00;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x140;	s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
		tmp += 0x140;
		two   = tmp + 0x0;	// AVX+ versions of Radix-32 DFT macros assume consts 2.0,1.0,sqrt2,isrt2 laid out thusly
		one   = tmp + 0x1;
		sqrt2 = tmp + 0x2;
		isrt2 = tmp + 0x3;
		xcc2  = tmp + 0x4;	// Radix-32 DFT macros assume roots stored in this [8th, 16th, 32nd_1,3] order
		xss2  = tmp + 0x5;
		xcc1  = tmp + 0x6;
		xss1  = tmp + 0x7;
		xcc3  = tmp + 0x8;
		xss3  = tmp + 0x9;
		ycc1  = tmp + 0xa;	// radix-5 DFT trig consts
		ycc2  = tmp + 0xb;
		yss1  = tmp + 0xc;
		yss2  = tmp + 0xd;
		yss3  = tmp + 0xe;
		tmp += 0x10;	// sc_ptr += 0x290
	  #ifdef USE_AVX512
		cy = tmp;		tmp += 0x14;	// RADIX/8 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
		cy = tmp;		tmp += 0x28;	// RADIX/4 vec_dbl slots
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(290 + 28 + 2) = 0x2ba; This is where the value of half_arr_offset160 comes from
		half_arr= tmp + 0x02;	// This table needs 20 vec_dbl in both avx and sse2 mode
	  #else
		cy = tmp;		tmp += 0x50;	// RADIX/2 vec_dbl slots
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(290 + 50 + 2) = 0x2e2; This is where the value of half_arr_offset160 comes from
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

		sign_mask = (uint64*)(r00 + radix160_creals_in_local_store);
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
		#include "radix160_main_carry_loop.h"

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
