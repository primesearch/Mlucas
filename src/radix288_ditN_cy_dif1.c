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

#define RADIX 288	// Use #define rather than const int to ensure it's really a compile-time const in the C sense
#define ODD_RADIX 9	// ODD_RADIX = [radix >> trailz(radix)]

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
  // Max = (40 [SSE2]; 132 [AVX]), add to (half_arr_offset288 + RADIX) to get required value of radix288_creals_in_local_store:
  #ifdef USE_AVX512	// RADIX/8 = 36 = 0x24 fewer carry slots than AVX:
	const int half_arr_offset288 = 0x4b8;	// + RADIX = 0x4b8 + 0x120 = 0x5d8; Used for thread local-storage-integrity checking
	const int radix288_creals_in_local_store = 0x65c;	// += 132 (=0x84) and round up to nearest multiple of 4
  #elif defined(USE_AVX)
	const int half_arr_offset288 = 0x4dc;	// + RADIX = 0x4dc + 0x120 = 0x5fc; Used for thread local-storage-integrity checking
	const int radix288_creals_in_local_store = 0x684;	// += 132 (=0x84) and round up to nearest multiple of 4
  #else
	const int half_arr_offset288 = 0x524;	// + RADIX = 0x524 + 0x120 = 0x644; Used for thread local-storage-integrity checking
	const int radix288_creals_in_local_store = 0x66c;	// += 40 (=0x28) and round up to nearest multiple of 4
  #endif

	#include "sse2_macro_gcc64.h"
	#include "radix09_sse_macro.h"
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

int radix288_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-288 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-288 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix288_ditN_cy_dif1";
	const int pfetch_dist = PFETCH_DIST;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
	static double wts_mult[2], inv_mult[2];	// Const wts-multiplier and 2*(its multiplicative inverse)
	double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
  #ifdef USE_AVX512
	const int jhi_wrap = 15;
  #else
	const int jhi_wrap =  7;
  #endif
	int NDIVR,i,incr,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,ntmp,outer,nbytes;
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 18|36|72 for avx512,avx,sse, respectively:
  #ifdef USE_AVX512	// Oddly, FMA-based builds appear to need slightly shorter chains for radices 144,288
	const int incr_long =  9,incr_med = 6,incr_short = 3;
  #elif defined(USE_AVX2)
	const int incr_long = 12,incr_med = 6,incr_short = 4;
  #else
	const int incr_long = 18,incr_med = 9,incr_short = 4;
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
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,p100,p110, nsave = 0;
	static int poff[RADIX>>2];	// Store mults of p4 offset for loop control
#ifndef MULTITHREAD
// Shared DIF+DIT:
	double rt,it,re;
	static int t_offsets[32];
	// Need storage for 2 circular-shifts perms of a basic 9-vector, with shift count in [0,8] that means 2*17 elts:
	static int dif_offsets[RADIX], dif_p20_cperms[34], dif_p20_lo_offset[32], dif_phi[ODD_RADIX];
	static int dit_offsets[RADIX], dit_p20_cperms[34], dit_p20_lo_offset[32], dit_phi[ODD_RADIX];
#endif
	static double radix_inv, n2inv;
#if defined(USE_SSE2) || !defined(MULTITHREAD)
	const double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
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
  #ifdef USE_AVX2
	// Due to GCC macro argc limit of 30, to enable 16-register data-doubled version of the radix-9 macros need 2 length-9 ptr arrays:
	vec_dbl *rad9_iptr[9], *rad9_optr[9];
  #endif

#ifdef USE_SSE2

	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else
	double *add0,*add1,*add2,*add3;	/* Addresses into array sections */
  #endif

	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
	const double crnd = 3.0*0x4000000*0x2000000;
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	vec_dbl
	#ifndef MULTITHREAD
		*va0,*va1,*va2,*va3,*va4,*va5,*va6,*va7,*va8,
		*vb0,*vb1,*vb2,*vb3,*vb4,*vb5,*vb6,*vb7,*vb8,
	#endif
		*tmp,*tm0,*tm1,*tm2;	// Non-static utility ptrs
	static vec_dbl *two,*one,*sqrt2,*isrt2,*xcc1,*xss1,*xcc2,*xss2,*xcc3,*xss3,	// radix-32 DFT trig consts
									*ycc1,*yss1,*ycc2,*yss2,*ycc3m1,*yss3,*ycc4,*yss4,	// radix-9 DFT trig consts
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
	static task_control_t   task_control = {NULL, (void*)cy288_process_chunk, NULL, 0x0};

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
		// consisting of radix288_creals_in_local_store dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix288_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix288_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = sc_ptr;	r00   = tmp;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x240;	s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
		tmp	+= 0x243;	// Extra 3 slots here for two,one below - added those late, too lazy to rejigger all the existing offsets following
		two     = tmp - 3;	// AVX+ versions of various DFT macros assume consts 2.0,1.0,isrt2 laid out thusly
		one     = tmp - 2;
		sqrt2	= tmp - 1;
		isrt2   = tmp + 0x0;
		xcc2	= tmp + 0x1;	// Radix-32 DFT macros assume roots stored in this [8th, 16th, 32nd_1,3] order
		xss2	= tmp + 0x2;
		xcc1	= tmp + 0x3;
		xss1	= tmp + 0x4;
		xcc3	= tmp + 0x5;
		xss3	= tmp + 0x6;
		// Roots for radix-9 DFTs:
		ycc1    = tmp + 0x7;
		yss1    = tmp + 0x8;
		ycc2    = tmp + 0x9;
		yss2    = tmp + 0xa;
		ycc3m1  = tmp + 0xb;
		yss3    = tmp + 0xc;
		ycc4    = tmp + 0xd;
		yss4    = tmp + 0xe;
		tmp += 0xf;	// sc_ptr += 0x492
	  #ifdef USE_AVX512
		cy = tmp;		tmp += 0x24;	// RADIX/8 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
		cy = tmp;		tmp += 0x48;	// RADIX/4 vec_dbl slots
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(492 + 48 + 2) = 0x4dc; This is where the value of half_arr_offset288 comes from
		half_arr= tmp + 0x02;	// This table needs 20 vec_dbl in both avx and sse2 mode
	  #else
		cy = tmp;		tmp += 0x90;	// RADIX/2 vec_dbl slots
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(492 + 90 + 2) = 0x524; This is where the value of half_arr_offset288 comes from
		half_arr= tmp + 0x02;
	  #endif
		ASSERT((radix288_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r00) + (20 << L2_SZ_VD), "radix288_creals_in_local_store checksum failed!");
		/* These remain fixed: */
		VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );
		VEC_DBL_INIT(sqrt2, SQRT2);	VEC_DBL_INIT(isrt2, ISRT2);
		VEC_DBL_INIT(xcc2, c16  );	VEC_DBL_INIT(xss2, s16  );	// radix-32 DFT trig consts
		VEC_DBL_INIT(xcc1, c32_1);	VEC_DBL_INIT(xss1, s32_1);
		VEC_DBL_INIT(xcc3, c32_3);	VEC_DBL_INIT(xss3, s32_3);
		VEC_DBL_INIT(ycc1  , c	);		VEC_DBL_INIT(yss1, s );	// radix-5 DFT trig consts
		VEC_DBL_INIT(ycc2  , c2  );		VEC_DBL_INIT(yss2, s2);
		VEC_DBL_INIT(ycc3m1, c3m1);		VEC_DBL_INIT(yss3, s3);
		VEC_DBL_INIT(ycc4  , c4  );		VEC_DBL_INIT(yss4, s4);
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
		nbytes += 64;;
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
		p1 = NDIVR;		p10 = NDIVR<<4;		p100 = NDIVR<<8;
		p2 = p1 + p1;	p20 = p10 + p10;	p110 = p100 + p10;
		p3 = p2 + p1;	p30 = p20 + p10;
		p4 = p3 + p1;	p40 = p30 + p10;
		p5 = p4 + p1;	p50 = p40 + p10;
		p6 = p5 + p1;	p60 = p50 + p10;
		p7 = p6 + p1;	p70 = p60 + p10;
		p8 = p7 + p1;	p80 = p70 + p10;
		p9 = p8 + p1;	p90 = p80 + p10;
		pa = p9 + p1;	pa0 = p90 + p10;
		pb = pa + p1;	pb0 = pa0 + p10;
		pc = pb + p1;	pc0 = pb0 + p10;
		pd = pc + p1;	pd0 = pc0 + p10;
		pe = pd + p1;	pe0 = pd0 + p10;
		pf = pe + p1;	pf0 = pe0 + p10;

		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );		p50 += ( (p50 >> DAT_BITS) << PAD_BITS );
		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );		p70 += ( (p70 >> DAT_BITS) << PAD_BITS );
		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pa += ( (pa >> DAT_BITS) << PAD_BITS );		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pb += ( (pb >> DAT_BITS) << PAD_BITS );		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		pc += ( (pc >> DAT_BITS) << PAD_BITS );		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pd += ( (pd >> DAT_BITS) << PAD_BITS );		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		pe += ( (pe >> DAT_BITS) << PAD_BITS );		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		pf += ( (pf >> DAT_BITS) << PAD_BITS );		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
													p100 += ( (p100 >> DAT_BITS) << PAD_BITS );
													p110 += ( (p110 >> DAT_BITS) << PAD_BITS );
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
		poff[0x40+0] =p100; poff[0x40+1] =p100+p4; poff[0x40+2] =p100+p8; poff[0x40+3] =p100+pc;
		poff[0x44+0] =p110; poff[0x44+1] =p110+p4; poff[0x44+2] =p110+p8; poff[0x44+3] =p110+pc;

	#ifndef MULTITHREAD

		l = 0;
		dif_phi[l++] =   0;
		dif_phi[l++] = p100;
		dif_phi[l++] = pa0;
		dif_phi[l++] = p40;
		dif_phi[l++] = pe0;
		dif_phi[l++] = p80;
		dif_phi[l++] = p20;
		dif_phi[l++] = pc0;
		dif_phi[l++] = p60;
		l = 0;
		dit_phi[l++] =   0;
		dit_phi[l++] = p20;
		dit_phi[l++] = p40;
		dit_phi[l++] = p80;
		dit_phi[l++] = pa0;
		dit_phi[l++] = p60;
		dit_phi[l++] =p100;
		dit_phi[l++] = pc0;
		dit_phi[l++] = pe0;

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

	   #ifdef USE_SSE2
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p10*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:

		// Init storage for 2 circular-shifts perms of a basic 9-vector, with shift count in [0,8] that means 2*17
		l = 0;
		// Even multiples of p10 cshift array:
		dif_p20_cperms[l++] = 0;
		dif_p20_cperms[l++] = 0x100<<1;
		dif_p20_cperms[l++] = 0xe0<<1;
		dif_p20_cperms[l++] = 0xc0<<1;
		dif_p20_cperms[l++] = 0xa0<<1;
		dif_p20_cperms[l++] = 0x80<<1;
		dif_p20_cperms[l++] = 0x60<<1;
		dif_p20_cperms[l++] = 0x40<<1;
		dif_p20_cperms[l++] = 0x20<<1;
		while(l < 2*ODD_RADIX-1) {
			dif_p20_cperms[l] = dif_p20_cperms[l - ODD_RADIX]; ++l;
		}
		// Odd multiples of p10 cshift array:
		dif_p20_cperms[l++] = 0x110<<1;
		dif_p20_cperms[l++] = 0xf0<<1;
		dif_p20_cperms[l++] = 0xd0<<1;
		dif_p20_cperms[l++] = 0xb0<<1;
		dif_p20_cperms[l++] = 0x90<<1;
		dif_p20_cperms[l++] = 0x70<<1;
		dif_p20_cperms[l++] = 0x50<<1;
		dif_p20_cperms[l++] = 0x30<<1;
		dif_p20_cperms[l++] = 0x10<<1;
		while(l < 4*ODD_RADIX-2) {
			dif_p20_cperms[l] = dif_p20_cperms[l - ODD_RADIX]; ++l;
		}
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-9 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-8; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-8 index, bits <4:30> for the p0-f. In SIMD mode again replace p[0-f] with 0x[0-f]<<1 for use with
		// contig-local-mem, thus replace 'p' prefixes with 0x and 4-left-shifts with << 5 to account for the << 1:
		l = 0;
		dif_p20_lo_offset[l++] = ((  0 << 5) + 0);
		dif_p20_lo_offset[l++] = ((0x7 << 5) + 0)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xe << 5) + 1);
		dif_p20_lo_offset[l++] = ((0x5 << 5) + 1);
		dif_p20_lo_offset[l++] = ((0xc << 5) + 1)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x3 << 5) + 1)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xa << 5) + 2);
		dif_p20_lo_offset[l++] = ((0x1 << 5) + 2);
		dif_p20_lo_offset[l++] = ((0x8 << 5) + 2)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xf << 5) + 3);
		dif_p20_lo_offset[l++] = ((0x6 << 5) + 3);
		dif_p20_lo_offset[l++] = ((0xd << 5) + 3)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x4 << 5) + 3)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xb << 5) + 4);
		dif_p20_lo_offset[l++] = ((0x2 << 5) + 4);
		dif_p20_lo_offset[l++] = ((0x9 << 5) + 4)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((  0 << 5) + 4)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x7 << 5) + 5);
		dif_p20_lo_offset[l++] = ((0xe << 5) + 5)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x5 << 5) + 5)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xc << 5) + 6);
		dif_p20_lo_offset[l++] = ((0x3 << 5) + 6);
		dif_p20_lo_offset[l++] = ((0xa << 5) + 6)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x1 << 5) + 6)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x8 << 5) + 7);
		dif_p20_lo_offset[l++] = ((0xf << 5) + 7)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x6 << 5) + 7)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xd << 5) + 8);
		dif_p20_lo_offset[l++] = ((0x4 << 5) + 8);
		dif_p20_lo_offset[l++] = ((0xb << 5) + 8)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x2 << 5) + 8)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x9 << 5) + 0);

	   #else

		// Init storage for 2 circular-shifts perms of a basic 9-vector, with shift count in [0,8] that means 2*17
		l = 0;
		// Even multiples of p10 cshift array:
		dif_p20_cperms[l++] = 0;
		dif_p20_cperms[l++] = p100;
		dif_p20_cperms[l++] = pe0;
		dif_p20_cperms[l++] = pc0;
		dif_p20_cperms[l++] = pa0;
		dif_p20_cperms[l++] = p80;
		dif_p20_cperms[l++] = p60;
		dif_p20_cperms[l++] = p40;
		dif_p20_cperms[l++] = p20;
		while(l < 2*ODD_RADIX-1) {
			dif_p20_cperms[l] = dif_p20_cperms[l - ODD_RADIX]; ++l;
		}
		// Odd multiples of p10 cshift array:
		dif_p20_cperms[l++] = p110;
		dif_p20_cperms[l++] = pf0;
		dif_p20_cperms[l++] = pd0;
		dif_p20_cperms[l++] = pb0;
		dif_p20_cperms[l++] = p90;
		dif_p20_cperms[l++] = p70;
		dif_p20_cperms[l++] = p50;
		dif_p20_cperms[l++] = p30;
		dif_p20_cperms[l++] = p10;
		while(l < 4*ODD_RADIX-2) {
			dif_p20_cperms[l] = dif_p20_cperms[l - ODD_RADIX]; ++l;
		}
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-9 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-8; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-8 index, bits <4:30> for the p0-f:
		l = 0;
		dif_p20_lo_offset[l++] = (( 0 << 4) + 0);
		dif_p20_lo_offset[l++] = ((p7 << 4) + 0)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pe << 4) + 1);
		dif_p20_lo_offset[l++] = ((p5 << 4) + 1);
		dif_p20_lo_offset[l++] = ((pc << 4) + 1)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p3 << 4) + 1)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pa << 4) + 2);
		dif_p20_lo_offset[l++] = ((p1 << 4) + 2);
		dif_p20_lo_offset[l++] = ((p8 << 4) + 2)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pf << 4) + 3);
		dif_p20_lo_offset[l++] = ((p6 << 4) + 3);
		dif_p20_lo_offset[l++] = ((pd << 4) + 3)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p4 << 4) + 3)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pb << 4) + 4);
		dif_p20_lo_offset[l++] = ((p2 << 4) + 4);
		dif_p20_lo_offset[l++] = ((p9 << 4) + 4)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = (( 0 << 4) + 4)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p7 << 4) + 5);
		dif_p20_lo_offset[l++] = ((pe << 4) + 5)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p5 << 4) + 5)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pc << 4) + 6);
		dif_p20_lo_offset[l++] = ((p3 << 4) + 6);
		dif_p20_lo_offset[l++] = ((pa << 4) + 6)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p1 << 4) + 6)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p8 << 4) + 7);
		dif_p20_lo_offset[l++] = ((pf << 4) + 7)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p6 << 4) + 7)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pd << 4) + 8);
		dif_p20_lo_offset[l++] = ((p4 << 4) + 8);
		dif_p20_lo_offset[l++] = ((pb << 4) + 8)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p2 << 4) + 8)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p9 << 4) + 0);

	   #endif	// sse2?

	// dif_offsets are w.r.to a-array, need 9 distinct sets of these, one for each DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		// Set 0: [0,1,3,2,7,6,5,4,e,f,c,d,a,b,8,9 + p00],[c,d,f,e,8,9,b,a,5,4,6,7,1,0,2,3 + p10]
		l = 0;
		dif_offsets[0x00] =  0;		dif_offsets[0x10] = pc+p10;
		dif_offsets[0x01] = p1;		dif_offsets[0x11] = pd+p10;
		dif_offsets[0x02] = p3;		dif_offsets[0x12] = pf+p10;
		dif_offsets[0x03] = p2;		dif_offsets[0x13] = pe+p10;
		dif_offsets[0x04] = p7;		dif_offsets[0x14] = p8+p10;
		dif_offsets[0x05] = p6;		dif_offsets[0x15] = p9+p10;
		dif_offsets[0x06] = p5;		dif_offsets[0x16] = pb+p10;
		dif_offsets[0x07] = p4;		dif_offsets[0x17] = pa+p10;
		dif_offsets[0x08] = pe;		dif_offsets[0x18] = p5+p10;
		dif_offsets[0x09] = pf;		dif_offsets[0x19] = p4+p10;
		dif_offsets[0x0a] = pc;		dif_offsets[0x1a] = p6+p10;
		dif_offsets[0x0b] = pd;		dif_offsets[0x1b] = p7+p10;
		dif_offsets[0x0c] = pa;		dif_offsets[0x1c] = p1+p10;
		dif_offsets[0x0d] = pb;		dif_offsets[0x1d] =    p10;
		dif_offsets[0x0e] = p8;		dif_offsets[0x1e] = p2+p10;
		dif_offsets[0x0f] = p9;		dif_offsets[0x1f] = p3+p10;
		// Set 1: [3,2,1,0,5,4,6,7,c,d,f,e,8,9,b,a +p100],[f,e,d,c,b,a,9,8,6,7,4,5,2,3,0,1 +p110]
		l += 32;
		dif_offsets[l+0x00] = p3;		dif_offsets[l+0x10] = pf+p10;
		dif_offsets[l+0x01] = p2;		dif_offsets[l+0x11] = pe+p10;
		dif_offsets[l+0x02] = p1;		dif_offsets[l+0x12] = pd+p10;
		dif_offsets[l+0x03] =  0;		dif_offsets[l+0x13] = pc+p10;
		dif_offsets[l+0x04] = p5;		dif_offsets[l+0x14] = pb+p10;
		dif_offsets[l+0x05] = p4;		dif_offsets[l+0x15] = pa+p10;
		dif_offsets[l+0x06] = p6;		dif_offsets[l+0x16] = p9+p10;
		dif_offsets[l+0x07] = p7;		dif_offsets[l+0x17] = p8+p10;
		dif_offsets[l+0x08] = pc;		dif_offsets[l+0x18] = p6+p10;
		dif_offsets[l+0x09] = pd;		dif_offsets[l+0x19] = p7+p10;
		dif_offsets[l+0x0a] = pf;		dif_offsets[l+0x1a] = p4+p10;
		dif_offsets[l+0x0b] = pe;		dif_offsets[l+0x1b] = p5+p10;
		dif_offsets[l+0x0c] = p8;		dif_offsets[l+0x1c] = p2+p10;
		dif_offsets[l+0x0d] = p9;		dif_offsets[l+0x1d] = p3+p10;
		dif_offsets[l+0x0e] = pb;		dif_offsets[l+0x1e] =    p10;
		dif_offsets[l+0x0f] = pa;		dif_offsets[l+0x1f] = p1+p10;
		// Set 2: [1,0,2,3,6,7,4,5,f,e,d,c,b,a,9,8 + pb0],[3,2,1,0,5,4,6,7,c,d,f,e,8,9,b,a + pa0]
		l += 32;
		dif_offsets[l+0x00] = p1+p10;		dif_offsets[l+0x10] = p3;
		dif_offsets[l+0x01] =    p10;		dif_offsets[l+0x11] = p2;
		dif_offsets[l+0x02] = p2+p10;		dif_offsets[l+0x12] = p1;
		dif_offsets[l+0x03] = p3+p10;		dif_offsets[l+0x13] =  0;
		dif_offsets[l+0x04] = p6+p10;		dif_offsets[l+0x14] = p5;
		dif_offsets[l+0x05] = p7+p10;		dif_offsets[l+0x15] = p4;
		dif_offsets[l+0x06] = p4+p10;		dif_offsets[l+0x16] = p6;
		dif_offsets[l+0x07] = p5+p10;		dif_offsets[l+0x17] = p7;
		dif_offsets[l+0x08] = pf+p10;		dif_offsets[l+0x18] = pc;
		dif_offsets[l+0x09] = pe+p10;		dif_offsets[l+0x19] = pd;
		dif_offsets[l+0x0a] = pd+p10;		dif_offsets[l+0x1a] = pf;
		dif_offsets[l+0x0b] = pc+p10;		dif_offsets[l+0x1b] = pe;
		dif_offsets[l+0x0c] = pb+p10;		dif_offsets[l+0x1c] = p8;
		dif_offsets[l+0x0d] = pa+p10;		dif_offsets[l+0x1d] = p9;
		dif_offsets[l+0x0e] = p9+p10;		dif_offsets[l+0x1e] = pb;
		dif_offsets[l+0x0f] = p8+p10;		dif_offsets[l+0x1f] = pa;
		// Set 3: [a,b,8,9,c,d,f,e,3,2,1,0,5,4,6,7 + p40],[1,0,2,3,6,7,4,5,f,e,d,c,b,a,9,8 + p50]
		l += 32;
		dif_offsets[l+0x00] = pa;		dif_offsets[l+0x10] = p1+p10;
		dif_offsets[l+0x01] = pb;		dif_offsets[l+0x11] =    p10;
		dif_offsets[l+0x02] = p8;		dif_offsets[l+0x12] = p2+p10;
		dif_offsets[l+0x03] = p9;		dif_offsets[l+0x13] = p3+p10;
		dif_offsets[l+0x04] = pc;		dif_offsets[l+0x14] = p6+p10;
		dif_offsets[l+0x05] = pd;		dif_offsets[l+0x15] = p7+p10;
		dif_offsets[l+0x06] = pf;		dif_offsets[l+0x16] = p4+p10;
		dif_offsets[l+0x07] = pe;		dif_offsets[l+0x17] = p5+p10;
		dif_offsets[l+0x08] = p3;		dif_offsets[l+0x18] = pf+p10;
		dif_offsets[l+0x09] = p2;		dif_offsets[l+0x19] = pe+p10;
		dif_offsets[l+0x0a] = p1;		dif_offsets[l+0x1a] = pd+p10;
		dif_offsets[l+0x0b] =  0;		dif_offsets[l+0x1b] = pc+p10;
		dif_offsets[l+0x0c] = p5;		dif_offsets[l+0x1c] = pb+p10;
		dif_offsets[l+0x0d] = p4;		dif_offsets[l+0x1d] = pa+p10;
		dif_offsets[l+0x0e] = p6;		dif_offsets[l+0x1e] = p9+p10;
		dif_offsets[l+0x0f] = p7;		dif_offsets[l+0x1f] = p8+p10;
		// Set 4: [8,9,b,a,f,e,d,c,1,0,2,3,6,7,4,5 + pf0],[a,b,8,9,c,d,f,e,3,2,1,0,5,4,6,7 + pe0]
		l += 32;
		dif_offsets[l+0x00] = p8+p10;		dif_offsets[l+0x10] = pa;
		dif_offsets[l+0x01] = p9+p10;		dif_offsets[l+0x11] = pb;
		dif_offsets[l+0x02] = pb+p10;		dif_offsets[l+0x12] = p8;
		dif_offsets[l+0x03] = pa+p10;		dif_offsets[l+0x13] = p9;
		dif_offsets[l+0x04] = pf+p10;		dif_offsets[l+0x14] = pc;
		dif_offsets[l+0x05] = pe+p10;		dif_offsets[l+0x15] = pd;
		dif_offsets[l+0x06] = pd+p10;		dif_offsets[l+0x16] = pf;
		dif_offsets[l+0x07] = pc+p10;		dif_offsets[l+0x17] = pe;
		dif_offsets[l+0x08] = p1+p10;		dif_offsets[l+0x18] = p3;
		dif_offsets[l+0x09] =    p10;		dif_offsets[l+0x19] = p2;
		dif_offsets[l+0x0a] = p2+p10;		dif_offsets[l+0x1a] = p1;
		dif_offsets[l+0x0b] = p3+p10;		dif_offsets[l+0x1b] =  0;
		dif_offsets[l+0x0c] = p6+p10;		dif_offsets[l+0x1c] = p5;
		dif_offsets[l+0x0d] = p7+p10;		dif_offsets[l+0x1d] = p4;
		dif_offsets[l+0x0e] = p4+p10;		dif_offsets[l+0x1e] = p6;
		dif_offsets[l+0x0f] = p5+p10;		dif_offsets[l+0x1f] = p7;
		// Set 5: [7,6,5,4,3,2,1,0,a,b,8,9,c,d,f,e + p80],[8,9,b,a,f,e,d,c,1,0,2,3,6,7,4,5 + p90]
		l += 32;
		dif_offsets[l+0x00] = p7;		dif_offsets[l+0x10] = p8+p10;
		dif_offsets[l+0x01] = p6;		dif_offsets[l+0x11] = p9+p10;
		dif_offsets[l+0x02] = p5;		dif_offsets[l+0x12] = pb+p10;
		dif_offsets[l+0x03] = p4;		dif_offsets[l+0x13] = pa+p10;
		dif_offsets[l+0x04] = p3;		dif_offsets[l+0x14] = pf+p10;
		dif_offsets[l+0x05] = p2;		dif_offsets[l+0x15] = pe+p10;
		dif_offsets[l+0x06] = p1;		dif_offsets[l+0x16] = pd+p10;
		dif_offsets[l+0x07] =  0;		dif_offsets[l+0x17] = pc+p10;
		dif_offsets[l+0x08] = pa;		dif_offsets[l+0x18] = p1+p10;
		dif_offsets[l+0x09] = pb;		dif_offsets[l+0x19] =    p10;
		dif_offsets[l+0x0a] = p8;		dif_offsets[l+0x1a] = p2+p10;
		dif_offsets[l+0x0b] = p9;		dif_offsets[l+0x1b] = p3+p10;
		dif_offsets[l+0x0c] = pc;		dif_offsets[l+0x1c] = p6+p10;
		dif_offsets[l+0x0d] = pd;		dif_offsets[l+0x1d] = p7+p10;
		dif_offsets[l+0x0e] = pf;		dif_offsets[l+0x1e] = p4+p10;
		dif_offsets[l+0x0f] = pe;		dif_offsets[l+0x1f] = p5+p10;
		// Set 6: [5,4,6,7,1,0,2,3,8,9,b,a,f,e,d,c + p30],[7,6,5,4,3,2,1,0,a,b,8,9,c,d,f,e + p20]
		l += 32;
		dif_offsets[l+0x00] = p5+p10;		dif_offsets[l+0x10] = p7;
		dif_offsets[l+0x01] = p4+p10;		dif_offsets[l+0x11] = p6;
		dif_offsets[l+0x02] = p6+p10;		dif_offsets[l+0x12] = p5;
		dif_offsets[l+0x03] = p7+p10;		dif_offsets[l+0x13] = p4;
		dif_offsets[l+0x04] = p1+p10;		dif_offsets[l+0x14] = p3;
		dif_offsets[l+0x05] =    p10;		dif_offsets[l+0x15] = p2;
		dif_offsets[l+0x06] = p2+p10;		dif_offsets[l+0x16] = p1;
		dif_offsets[l+0x07] = p3+p10;		dif_offsets[l+0x17] =  0;
		dif_offsets[l+0x08] = p8+p10;		dif_offsets[l+0x18] = pa;
		dif_offsets[l+0x09] = p9+p10;		dif_offsets[l+0x19] = pb;
		dif_offsets[l+0x0a] = pb+p10;		dif_offsets[l+0x1a] = p8;
		dif_offsets[l+0x0b] = pa+p10;		dif_offsets[l+0x1b] = p9;
		dif_offsets[l+0x0c] = pf+p10;		dif_offsets[l+0x1c] = pc;
		dif_offsets[l+0x0d] = pe+p10;		dif_offsets[l+0x1d] = pd;
		dif_offsets[l+0x0e] = pd+p10;		dif_offsets[l+0x1e] = pf;
		dif_offsets[l+0x0f] = pc+p10;		dif_offsets[l+0x1f] = pe;
		// Set 7: [e,f,c,d,a,b,8,9,7,6,5,4,3,2,1,0 + pc0],[5,4,6,7,1,0,2,3,8,9,b,a,f,e,d,c + pd0]
		l += 32;
		dif_offsets[l+0x00] = pe;		dif_offsets[l+0x10] = p5+p10;
		dif_offsets[l+0x01] = pf;		dif_offsets[l+0x11] = p4+p10;
		dif_offsets[l+0x02] = pc;		dif_offsets[l+0x12] = p6+p10;
		dif_offsets[l+0x03] = pd;		dif_offsets[l+0x13] = p7+p10;
		dif_offsets[l+0x04] = pa;		dif_offsets[l+0x14] = p1+p10;
		dif_offsets[l+0x05] = pb;		dif_offsets[l+0x15] =    p10;
		dif_offsets[l+0x06] = p8;		dif_offsets[l+0x16] = p2+p10;
		dif_offsets[l+0x07] = p9;		dif_offsets[l+0x17] = p3+p10;
		dif_offsets[l+0x08] = p7;		dif_offsets[l+0x18] = p8+p10;
		dif_offsets[l+0x09] = p6;		dif_offsets[l+0x19] = p9+p10;
		dif_offsets[l+0x0a] = p5;		dif_offsets[l+0x1a] = pb+p10;
		dif_offsets[l+0x0b] = p4;		dif_offsets[l+0x1b] = pa+p10;
		dif_offsets[l+0x0c] = p3;		dif_offsets[l+0x1c] = pf+p10;
		dif_offsets[l+0x0d] = p2;		dif_offsets[l+0x1d] = pe+p10;
		dif_offsets[l+0x0e] = p1;		dif_offsets[l+0x1e] = pd+p10;
		dif_offsets[l+0x0f] =  0;		dif_offsets[l+0x1f] = pc+p10;
		// Set 8: [c,d,f,e,8,9,b,a,5,4,6,7,1,0,2,3 + p70],[e,f,c,d,a,b,8,9,7,6,5,4,3,2,1,0 + p60]
		l += 32;
		dif_offsets[l+0x00] = pc+p10;		dif_offsets[l+0x10] = pe;
		dif_offsets[l+0x01] = pd+p10;		dif_offsets[l+0x11] = pf;
		dif_offsets[l+0x02] = pf+p10;		dif_offsets[l+0x12] = pc;
		dif_offsets[l+0x03] = pe+p10;		dif_offsets[l+0x13] = pd;
		dif_offsets[l+0x04] = p8+p10;		dif_offsets[l+0x14] = pa;
		dif_offsets[l+0x05] = p9+p10;		dif_offsets[l+0x15] = pb;
		dif_offsets[l+0x06] = pb+p10;		dif_offsets[l+0x16] = p8;
		dif_offsets[l+0x07] = pa+p10;		dif_offsets[l+0x17] = p9;
		dif_offsets[l+0x08] = p5+p10;		dif_offsets[l+0x18] = p7;
		dif_offsets[l+0x09] = p4+p10;		dif_offsets[l+0x19] = p6;
		dif_offsets[l+0x0a] = p6+p10;		dif_offsets[l+0x1a] = p5;
		dif_offsets[l+0x0b] = p7+p10;		dif_offsets[l+0x1b] = p4;
		dif_offsets[l+0x0c] = p1+p10;		dif_offsets[l+0x1c] = p3;
		dif_offsets[l+0x0d] =    p10;		dif_offsets[l+0x1d] = p2;
		dif_offsets[l+0x0e] = p2+p10;		dif_offsets[l+0x1e] = p1;
		dif_offsets[l+0x0f] = p3+p10;		dif_offsets[l+0x1f] =  0;
	  #ifdef USE_SSE2
		// IN SIMD mode preshift all the above offsets << 3 to turn into double=array pointer offsets:
		for(l = 0; l < RADIX; l++) {
			dif_offsets[l] <<= 3;
		}
	  #endif

	/*** DIT indexing stuff: ***/

	   #ifdef USE_SSE2
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p10*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:
		l = 0;
		// Init storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,8] that means 2*17
		// Even multiples of p10 cshift array: evn0 := 0x10*[00,02,04,06,08,0a,0c,0e,10]
		dit_p20_cperms[l++] =    0;
		dit_p20_cperms[l++] = 0x20<<1;
		dit_p20_cperms[l++] = 0x40<<1;
		dit_p20_cperms[l++] = 0x60<<1;
		dit_p20_cperms[l++] = 0x80<<1;
		dit_p20_cperms[l++] = 0xa0<<1;
		dit_p20_cperms[l++] = 0xc0<<1;
		dit_p20_cperms[l++] = 0xe0<<1;
		dit_p20_cperms[l++] = 0x100<<1;
		while(l < 2*ODD_RADIX-1) {
			dit_p20_cperms[l] = dit_p20_cperms[l - ODD_RADIX]; ++l;
		}
		// Odd multiples of p10 cshift array: odd0 := 0x10*[03,05,07,09,0b,0d,0f,11,01]
		dit_p20_cperms[l++] = 0x30<<1;
		dit_p20_cperms[l++] = 0x50<<1;
		dit_p20_cperms[l++] = 0x70<<1;
		dit_p20_cperms[l++] = 0x90<<1;
		dit_p20_cperms[l++] = 0xb0<<1;
		dit_p20_cperms[l++] = 0xd0<<1;
		dit_p20_cperms[l++] = 0xf0<<1;
		dit_p20_cperms[l++] = 0x110<<1;
		dit_p20_cperms[l++] = 0x10<<1;
		while(l < 4*ODD_RADIX-2) {
			dit_p20_cperms[l] = dit_p20_cperms[l - ODD_RADIX]; ++l;
		}

		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-9 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-8 index, bits <4:30> for the p0-f. In SIMD mode again replace p[0-f] with 0x[0-f]<<1 for use with
		// contig-local-mem, thus replace 'p' prefixes with 0x and 4-left-shifts with << 5 to account for the << 1:
		l = 0;
		dit_p20_lo_offset[l++] = ((  0 << 5) + 0);
		dit_p20_lo_offset[l++] = ((0xf << 5) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xe << 5) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xd << 5) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xc << 5) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xb << 5) + 8) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xa << 5) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x9 << 5) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x8 << 5) + 5) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x7 << 5) + 7) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x6 << 5) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x5 << 5) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x4 << 5) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x3 << 5) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x2 << 5) + 8) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x1 << 5) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((  0 << 5) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xf << 5) + 6);
		dit_p20_lo_offset[l++] = ((0xe << 5) + 8);
		dit_p20_lo_offset[l++] = ((0xd << 5) + 1);
		dit_p20_lo_offset[l++] = ((0xc << 5) + 3);
		dit_p20_lo_offset[l++] = ((0xb << 5) + 5);
		dit_p20_lo_offset[l++] = ((0xa << 5) + 7);
		dit_p20_lo_offset[l++] = ((0x9 << 5) + 0);
		dit_p20_lo_offset[l++] = ((0x8 << 5) + 2);
		dit_p20_lo_offset[l++] = ((0x7 << 5) + 4);
		dit_p20_lo_offset[l++] = ((0x6 << 5) + 6);
		dit_p20_lo_offset[l++] = ((0x5 << 5) + 8);
		dit_p20_lo_offset[l++] = ((0x4 << 5) + 1);
		dit_p20_lo_offset[l++] = ((0x3 << 5) + 3);
		dit_p20_lo_offset[l++] = ((0x2 << 5) + 5);
		dit_p20_lo_offset[l++] = ((0x1 << 5) + 7);

	   #else

		l = 0;
		// Init storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,8] that means 2*17
		// Even multiples of p10 cshift array: evn0 := p10*[00,02,04,06,08,0a,0c,0e,10]
		dit_p20_cperms[l++] =   0;
		dit_p20_cperms[l++] = p20;
		dit_p20_cperms[l++] = p40;
		dit_p20_cperms[l++] = p60;
		dit_p20_cperms[l++] = p80;
		dit_p20_cperms[l++] = pa0;
		dit_p20_cperms[l++] = pc0;
		dit_p20_cperms[l++] = pe0;
		dit_p20_cperms[l++] = p100;
		while(l < 2*ODD_RADIX-1) {
			dit_p20_cperms[l] = dit_p20_cperms[l - ODD_RADIX]; ++l;
		}
		// Odd multiples of p10 cshift array: odd0 := p10*[03,05,07,09,0b,0d,0f,11,01]
		dit_p20_cperms[l++] = p30;
		dit_p20_cperms[l++] = p50;
		dit_p20_cperms[l++] = p70;
		dit_p20_cperms[l++] = p90;
		dit_p20_cperms[l++] = pb0;
		dit_p20_cperms[l++] = pd0;
		dit_p20_cperms[l++] = pf0;
		dit_p20_cperms[l++] = p110;
		dit_p20_cperms[l++] = p10;
		while(l < 4*ODD_RADIX-2) {
			dit_p20_cperms[l] = dit_p20_cperms[l - ODD_RADIX]; ++l;
		}

		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-9 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-4 index, bits <4:30> for the p0-f:
		l = 0;
		dit_p20_lo_offset[l++] = (( 0 << 4) + 0);
		dit_p20_lo_offset[l++] = ((pf << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pe << 4) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pd << 4) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pc << 4) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pb << 4) + 8) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pa << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p9 << 4) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p8 << 4) + 5) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p7 << 4) + 7) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p6 << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p5 << 4) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p4 << 4) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p3 << 4) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p2 << 4) + 8) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p1 << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = (( 0 << 4) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pf << 4) + 6);
		dit_p20_lo_offset[l++] = ((pe << 4) + 8);
		dit_p20_lo_offset[l++] = ((pd << 4) + 1);
		dit_p20_lo_offset[l++] = ((pc << 4) + 3);
		dit_p20_lo_offset[l++] = ((pb << 4) + 5);
		dit_p20_lo_offset[l++] = ((pa << 4) + 7);
		dit_p20_lo_offset[l++] = ((p9 << 4) + 0);
		dit_p20_lo_offset[l++] = ((p8 << 4) + 2);
		dit_p20_lo_offset[l++] = ((p7 << 4) + 4);
		dit_p20_lo_offset[l++] = ((p6 << 4) + 6);
		dit_p20_lo_offset[l++] = ((p5 << 4) + 8);
		dit_p20_lo_offset[l++] = ((p4 << 4) + 1);
		dit_p20_lo_offset[l++] = ((p3 << 4) + 3);
		dit_p20_lo_offset[l++] = ((p2 << 4) + 5);
		dit_p20_lo_offset[l++] = ((p1 << 4) + 7);

	   #endif	// sse2?

	// dit_offsets are w.r.to a-array, need 9 distinct sets of these, one for each radix-32 DFT.
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
		// Set 1: [5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p30],[5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p20]
		l += 32;
		dit_offsets[l+0x00] = p5+p10;		dit_offsets[l+0x10] = p5;
		dit_offsets[l+0x01] = p4+p10;		dit_offsets[l+0x11] = p4;
		dit_offsets[l+0x02] = p6+p10;		dit_offsets[l+0x12] = p6;
		dit_offsets[l+0x03] = p7+p10;		dit_offsets[l+0x13] = p7;
		dit_offsets[l+0x04] = p1+p10;		dit_offsets[l+0x14] = p1;
		dit_offsets[l+0x05] =    p10;		dit_offsets[l+0x15] =  0;
		dit_offsets[l+0x06] = p2+p10;		dit_offsets[l+0x16] = p2;
		dit_offsets[l+0x07] = p3+p10;		dit_offsets[l+0x17] = p3;
		dit_offsets[l+0x08] = p9+p10;		dit_offsets[l+0x18] = p9;
		dit_offsets[l+0x09] = p8+p10;		dit_offsets[l+0x19] = p8;
		dit_offsets[l+0x0a] = pa+p10;		dit_offsets[l+0x1a] = pa;
		dit_offsets[l+0x0b] = pb+p10;		dit_offsets[l+0x1b] = pb;
		dit_offsets[l+0x0c] = pe+p10;		dit_offsets[l+0x1c] = pe;
		dit_offsets[l+0x0d] = pf+p10;		dit_offsets[l+0x1d] = pf;
		dit_offsets[l+0x0e] = pc+p10;		dit_offsets[l+0x1e] = pc;
		dit_offsets[l+0x0f] = pd+p10;		dit_offsets[l+0x1f] = pd;
		// Set 2: [a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p40],[2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p50]
		l += 32;
		dit_offsets[l+0x00] = pa;		dit_offsets[l+0x10] = p2+p10;
		dit_offsets[l+0x01] = pb;		dit_offsets[l+0x11] = p3+p10;
		dit_offsets[l+0x02] = p8;		dit_offsets[l+0x12] =    p10;
		dit_offsets[l+0x03] = p9;		dit_offsets[l+0x13] = p1+p10;
		dit_offsets[l+0x04] = pc;		dit_offsets[l+0x14] = p4+p10;
		dit_offsets[l+0x05] = pd;		dit_offsets[l+0x15] = p5+p10;
		dit_offsets[l+0x06] = pf;		dit_offsets[l+0x16] = p7+p10;
		dit_offsets[l+0x07] = pe;		dit_offsets[l+0x17] = p6+p10;
		dit_offsets[l+0x08] = p2;		dit_offsets[l+0x18] = pc+p10;
		dit_offsets[l+0x09] = p3;		dit_offsets[l+0x19] = pd+p10;
		dit_offsets[l+0x0a] =  0;		dit_offsets[l+0x1a] = pf+p10;
		dit_offsets[l+0x0b] = p1;		dit_offsets[l+0x1b] = pe+p10;
		dit_offsets[l+0x0c] = p4;		dit_offsets[l+0x1c] = p8+p10;
		dit_offsets[l+0x0d] = p5;		dit_offsets[l+0x1d] = p9+p10;
		dit_offsets[l+0x0e] = p7;		dit_offsets[l+0x1e] = pb+p10;
		dit_offsets[l+0x0f] = p6;		dit_offsets[l+0x1f] = pa+p10;
		// Set 3: [7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + p80],[b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p90]
		l += 32;
		dit_offsets[l+0x00] = p7;		dit_offsets[l+0x10] = pb+p10;
		dit_offsets[l+0x01] = p6;		dit_offsets[l+0x11] = pa+p10;
		dit_offsets[l+0x02] = p5;		dit_offsets[l+0x12] = p9+p10;
		dit_offsets[l+0x03] = p4;		dit_offsets[l+0x13] = p8+p10;
		dit_offsets[l+0x04] = p3;		dit_offsets[l+0x14] = pd+p10;
		dit_offsets[l+0x05] = p2;		dit_offsets[l+0x15] = pc+p10;
		dit_offsets[l+0x06] = p1;		dit_offsets[l+0x16] = pe+p10;
		dit_offsets[l+0x07] =  0;		dit_offsets[l+0x17] = pf+p10;
		dit_offsets[l+0x08] = pb;		dit_offsets[l+0x18] = p3+p10;
		dit_offsets[l+0x09] = pa;		dit_offsets[l+0x19] = p2+p10;
		dit_offsets[l+0x0a] = p9;		dit_offsets[l+0x1a] = p1+p10;
		dit_offsets[l+0x0b] = p8;		dit_offsets[l+0x1b] =    p10;
		dit_offsets[l+0x0c] = pd;		dit_offsets[l+0x1c] = p5+p10;
		dit_offsets[l+0x0d] = pc;		dit_offsets[l+0x1d] = p4+p10;
		dit_offsets[l+0x0e] = pe;		dit_offsets[l+0x1e] = p6+p10;
		dit_offsets[l+0x0f] = pf;		dit_offsets[l+0x1f] = p7+p10;
		// Set 4: [1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pb0],[1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pa0]
		l += 32;
		dit_offsets[l+0x00] = p1+p10;		dit_offsets[l+0x10] = p1;
		dit_offsets[l+0x01] =    p10;		dit_offsets[l+0x11] =  0;
		dit_offsets[l+0x02] = p2+p10;		dit_offsets[l+0x12] = p2;
		dit_offsets[l+0x03] = p3+p10;		dit_offsets[l+0x13] = p3;
		dit_offsets[l+0x04] = p6+p10;		dit_offsets[l+0x14] = p6;
		dit_offsets[l+0x05] = p7+p10;		dit_offsets[l+0x15] = p7;
		dit_offsets[l+0x06] = p4+p10;		dit_offsets[l+0x16] = p4;
		dit_offsets[l+0x07] = p5+p10;		dit_offsets[l+0x17] = p5;
		dit_offsets[l+0x08] = pe+p10;		dit_offsets[l+0x18] = pe;
		dit_offsets[l+0x09] = pf+p10;		dit_offsets[l+0x19] = pf;
		dit_offsets[l+0x0a] = pc+p10;		dit_offsets[l+0x1a] = pc;
		dit_offsets[l+0x0b] = pd+p10;		dit_offsets[l+0x1b] = pd;
		dit_offsets[l+0x0c] = pa+p10;		dit_offsets[l+0x1c] = pa;
		dit_offsets[l+0x0d] = pb+p10;		dit_offsets[l+0x1d] = pb;
		dit_offsets[l+0x0e] = p8+p10;		dit_offsets[l+0x1e] = p8;
		dit_offsets[l+0x0f] = p9+p10;		dit_offsets[l+0x1f] = p9;
		// Set 5: [c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p70],[c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p60]
		l += 32;
		dit_offsets[l+0x00] = pc+p10;		dit_offsets[l+0x10] = pc;
		dit_offsets[l+0x01] = pd+p10;		dit_offsets[l+0x11] = pd;
		dit_offsets[l+0x02] = pf+p10;		dit_offsets[l+0x12] = pf;
		dit_offsets[l+0x03] = pe+p10;		dit_offsets[l+0x13] = pe;
		dit_offsets[l+0x04] = p8+p10;		dit_offsets[l+0x14] = p8;
		dit_offsets[l+0x05] = p9+p10;		dit_offsets[l+0x15] = p9;
		dit_offsets[l+0x06] = pb+p10;		dit_offsets[l+0x16] = pb;
		dit_offsets[l+0x07] = pa+p10;		dit_offsets[l+0x17] = pa;
		dit_offsets[l+0x08] = p4+p10;		dit_offsets[l+0x18] = p4;
		dit_offsets[l+0x09] = p5+p10;		dit_offsets[l+0x19] = p5;
		dit_offsets[l+0x0a] = p7+p10;		dit_offsets[l+0x1a] = p7;
		dit_offsets[l+0x0b] = p6+p10;		dit_offsets[l+0x1b] = p6;
		dit_offsets[l+0x0c] =    p10;		dit_offsets[l+0x1c] =  0;
		dit_offsets[l+0x0d] = p1+p10;		dit_offsets[l+0x1d] = p1;
		dit_offsets[l+0x0e] = p3+p10;		dit_offsets[l+0x1e] = p3;
		dit_offsets[l+0x0f] = p2+p10;		dit_offsets[l+0x1f] = p2;
		// Set 6: [3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b +p100],[d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 +p110]
		l += 32;
		dit_offsets[l+0x00] = p3;		dit_offsets[l+0x10] = pd+p10;
		dit_offsets[l+0x01] = p2;		dit_offsets[l+0x11] = pc+p10;
		dit_offsets[l+0x02] = p1;		dit_offsets[l+0x12] = pe+p10;
		dit_offsets[l+0x03] =  0;		dit_offsets[l+0x13] = pf+p10;
		dit_offsets[l+0x04] = p5;		dit_offsets[l+0x14] = p9+p10;
		dit_offsets[l+0x05] = p4;		dit_offsets[l+0x15] = p8+p10;
		dit_offsets[l+0x06] = p6;		dit_offsets[l+0x16] = pa+p10;
		dit_offsets[l+0x07] = p7;		dit_offsets[l+0x17] = pb+p10;
		dit_offsets[l+0x08] = pd;		dit_offsets[l+0x18] = p5+p10;
		dit_offsets[l+0x09] = pc;		dit_offsets[l+0x19] = p4+p10;
		dit_offsets[l+0x0a] = pe;		dit_offsets[l+0x1a] = p6+p10;
		dit_offsets[l+0x0b] = pf;		dit_offsets[l+0x1b] = p7+p10;
		dit_offsets[l+0x0c] = p9;		dit_offsets[l+0x1c] = p1+p10;
		dit_offsets[l+0x0d] = p8;		dit_offsets[l+0x1d] =    p10;
		dit_offsets[l+0x0e] = pa;		dit_offsets[l+0x1e] = p2+p10;
		dit_offsets[l+0x0f] = pb;		dit_offsets[l+0x1f] = p3+p10;
		// Set 7: [e,f,c,d,a,b,8,9,6,7,4,5,2,3,0,1 + pc0],[6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pd0]
		l += 32;
		dit_offsets[l+0x00] = pe;		dit_offsets[l+0x10] = p6+p10;
		dit_offsets[l+0x01] = pf;		dit_offsets[l+0x11] = p7+p10;
		dit_offsets[l+0x02] = pc;		dit_offsets[l+0x12] = p4+p10;
		dit_offsets[l+0x03] = pd;		dit_offsets[l+0x13] = p5+p10;
		dit_offsets[l+0x04] = pa;		dit_offsets[l+0x14] = p2+p10;
		dit_offsets[l+0x05] = pb;		dit_offsets[l+0x15] = p3+p10;
		dit_offsets[l+0x06] = p8;		dit_offsets[l+0x16] =    p10;
		dit_offsets[l+0x07] = p9;		dit_offsets[l+0x17] = p1+p10;
		dit_offsets[l+0x08] = p6;		dit_offsets[l+0x18] = pa+p10;
		dit_offsets[l+0x09] = p7;		dit_offsets[l+0x19] = pb+p10;
		dit_offsets[l+0x0a] = p4;		dit_offsets[l+0x1a] = p8+p10;
		dit_offsets[l+0x0b] = p5;		dit_offsets[l+0x1b] = p9+p10;
		dit_offsets[l+0x0c] = p2;		dit_offsets[l+0x1c] = pc+p10;
		dit_offsets[l+0x0d] = p3;		dit_offsets[l+0x1d] = pd+p10;
		dit_offsets[l+0x0e] =  0;		dit_offsets[l+0x1e] = pf+p10;
		dit_offsets[l+0x0f] = p1;		dit_offsets[l+0x1f] = pe+p10;
		// Set 8: [8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + pf0],[8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + pe0]
		l += 32;
		dit_offsets[l+0x00] = p8+p10;		dit_offsets[l+0x10] = p8;
		dit_offsets[l+0x01] = p9+p10;		dit_offsets[l+0x11] = p9;
		dit_offsets[l+0x02] = pb+p10;		dit_offsets[l+0x12] = pb;
		dit_offsets[l+0x03] = pa+p10;		dit_offsets[l+0x13] = pa;
		dit_offsets[l+0x04] = pf+p10;		dit_offsets[l+0x14] = pf;
		dit_offsets[l+0x05] = pe+p10;		dit_offsets[l+0x15] = pe;
		dit_offsets[l+0x06] = pd+p10;		dit_offsets[l+0x16] = pd;
		dit_offsets[l+0x07] = pc+p10;		dit_offsets[l+0x17] = pc;
		dit_offsets[l+0x08] =    p10;		dit_offsets[l+0x18] =  0;
		dit_offsets[l+0x09] = p1+p10;		dit_offsets[l+0x19] = p1;
		dit_offsets[l+0x0a] = p3+p10;		dit_offsets[l+0x1a] = p3;
		dit_offsets[l+0x0b] = p2+p10;		dit_offsets[l+0x1b] = p2;
		dit_offsets[l+0x0c] = p7+p10;		dit_offsets[l+0x1c] = p7;
		dit_offsets[l+0x0d] = p6+p10;		dit_offsets[l+0x1d] = p6;
		dit_offsets[l+0x0e] = p5+p10;		dit_offsets[l+0x1e] = p5;
		dit_offsets[l+0x0f] = p4+p10;		dit_offsets[l+0x1f] = p4;
	  #ifdef USE_SSE2
		// IN SIMD mode preshift all the above offsets << 3 to turn into double=array pointer offsets:
		for(l = 0; l < RADIX; l++) {
			dit_offsets[l] <<= 3;
		}
	  #endif

	#endif	// #ifndef MULTITHREAD

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


/*...The radix-288 final DIT pass is here.	*/

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
		#include "radix288_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
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
		ASSERT(0x0 == cy288_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
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
	//	printf("Iter = %d, prp_mult = %4.1f, maxerr = %20.15f\n",iter,prp_mult,maxerr);
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

void radix288_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!...Subroutine to perform an initial radix-288 complex DIF FFT pass on the data in the length-N real vector A.
!   See the documentation in radix9_dif_pass for details on the radix-9 subtransform.
*/
	const double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
	int i,j,j1,j2,jt,jp;
	// p-indexing is hexadecimal here:
	static int NDIVR,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,p100,p110, first_entry=TRUE;
	static int t_offsets[32], dif_offsets[RADIX];
	// Need storage for 2 circular-shifts perms of a basic 9-vector, with shift count in [0,8] that means 2*17 elts:
	static int dif_p20_cperms[34], dif_p20_lo_offset[32], dif_phi[ODD_RADIX];
	double rt,it,re;
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
		p1 = NDIVR;		p10 = NDIVR<<4;		p100 = NDIVR<<8;
		p2 = p1 + p1;	p20 = p10 + p10;	p110 = p100 + p10;
		p3 = p2 + p1;	p30 = p20 + p10;
		p4 = p3 + p1;	p40 = p30 + p10;
		p5 = p4 + p1;	p50 = p40 + p10;
		p6 = p5 + p1;	p60 = p50 + p10;
		p7 = p6 + p1;	p70 = p60 + p10;
		p8 = p7 + p1;	p80 = p70 + p10;
		p9 = p8 + p1;	p90 = p80 + p10;
		pa = p9 + p1;	pa0 = p90 + p10;
		pb = pa + p1;	pb0 = pa0 + p10;
		pc = pb + p1;	pc0 = pb0 + p10;
		pd = pc + p1;	pd0 = pc0 + p10;
		pe = pd + p1;	pe0 = pd0 + p10;
		pf = pe + p1;	pf0 = pe0 + p10;

		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );		p50 += ( (p50 >> DAT_BITS) << PAD_BITS );
		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );		p70 += ( (p70 >> DAT_BITS) << PAD_BITS );
		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pa += ( (pa >> DAT_BITS) << PAD_BITS );		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pb += ( (pb >> DAT_BITS) << PAD_BITS );		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		pc += ( (pc >> DAT_BITS) << PAD_BITS );		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pd += ( (pd >> DAT_BITS) << PAD_BITS );		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		pe += ( (pe >> DAT_BITS) << PAD_BITS );		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		pf += ( (pf >> DAT_BITS) << PAD_BITS );		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
													p100 += ( (p100 >> DAT_BITS) << PAD_BITS );
													p110 += ( (p110 >> DAT_BITS) << PAD_BITS );
		i = 0;				// p-offsets for outputs have 'big' part of form p[0,16,10,4,14,8,2,12,6]0,
		dif_phi[i++] =   0;	// i.e. multiple-of-16 part gets decremented by 6 (mod 16) each term of the
		dif_phi[i++] = p100;// sequence, aside from the initial pair of terms, 0,16.
		dif_phi[i++] = pa0;
		dif_phi[i++] = p40;
		dif_phi[i++] = pe0;
		dif_phi[i++] = p80;
		dif_phi[i++] = p20;
		dif_phi[i++] = pc0;
		dif_phi[i++] = p60;

		// Set array offsets for radix-32 DFT in/outputs:
		// t_offsets w.r.to: t-array, same for all 9 32-DFTs:
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

		// Init storage for 2 circular-shifts perms of a basic 9-vector, with shift count in [0,8] that means 2*17 = 34:
		i = 0;
		// Even multiples of p10 cshift array:
		dif_p20_cperms[i++] = 0;
		dif_p20_cperms[i++] = p100;
		dif_p20_cperms[i++] = pe0;
		dif_p20_cperms[i++] = pc0;
		dif_p20_cperms[i++] = pa0;
		dif_p20_cperms[i++] = p80;
		dif_p20_cperms[i++] = p60;
		dif_p20_cperms[i++] = p40;
		dif_p20_cperms[i++] = p20;
		while(i < 2*ODD_RADIX-1) {
			dif_p20_cperms[i] = dif_p20_cperms[i - ODD_RADIX]; ++i;
		}
		// Odd multiples of p10 cshift array:
		dif_p20_cperms[i++] = p110;
		dif_p20_cperms[i++] = pf0;
		dif_p20_cperms[i++] = pd0;
		dif_p20_cperms[i++] = pb0;
		dif_p20_cperms[i++] = p90;
		dif_p20_cperms[i++] = p70;
		dif_p20_cperms[i++] = p50;
		dif_p20_cperms[i++] = p30;
		dif_p20_cperms[i++] = p10;
		while(i < 4*ODD_RADIX-2) {
			dif_p20_cperms[i] = dif_p20_cperms[i - ODD_RADIX]; ++i;
		}
		// Low parts, i.e. (mod p20) of the p-index offsets in the below circ-perm-indexing scheme for the radix-9 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-8; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-8 index, bits <4:30> for the p0-f:
		i = 0;
		dif_p20_lo_offset[i++] = (( 0 << 4) + 0);
		dif_p20_lo_offset[i++] = ((p7 << 4) + 0)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pe << 4) + 1);
		dif_p20_lo_offset[i++] = ((p5 << 4) + 1);
		dif_p20_lo_offset[i++] = ((pc << 4) + 1)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p3 << 4) + 1)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pa << 4) + 2);
		dif_p20_lo_offset[i++] = ((p1 << 4) + 2);
		dif_p20_lo_offset[i++] = ((p8 << 4) + 2)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pf << 4) + 3);
		dif_p20_lo_offset[i++] = ((p6 << 4) + 3);
		dif_p20_lo_offset[i++] = ((pd << 4) + 3)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p4 << 4) + 3)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pb << 4) + 4);
		dif_p20_lo_offset[i++] = ((p2 << 4) + 4);
		dif_p20_lo_offset[i++] = ((p9 << 4) + 4)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = (( 0 << 4) + 4)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p7 << 4) + 5);
		dif_p20_lo_offset[i++] = ((pe << 4) + 5)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p5 << 4) + 5)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pc << 4) + 6);
		dif_p20_lo_offset[i++] = ((p3 << 4) + 6);
		dif_p20_lo_offset[i++] = ((pa << 4) + 6)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p1 << 4) + 6)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p8 << 4) + 7);
		dif_p20_lo_offset[i++] = ((pf << 4) + 7)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p6 << 4) + 7)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pd << 4) + 8);
		dif_p20_lo_offset[i++] = ((p4 << 4) + 8);
		dif_p20_lo_offset[i++] = ((pb << 4) + 8)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p2 << 4) + 8)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p9 << 4) + 0);
	#if 0	// Monotone initial-o-indices used for extracting the needed operm:
		i = 0;
		dif_offsets[0x00] =  0;
		dif_offsets[0x01] = p1;
		dif_offsets[0x02] = p2;
		dif_offsets[0x03] = p3;
		dif_offsets[0x04] = p4;
		dif_offsets[0x05] = p5;
		dif_offsets[0x06] = p6;
		dif_offsets[0x07] = p7;
		dif_offsets[0x08] = p8;
		dif_offsets[0x09] = p9;
		dif_offsets[0x0a] = pa;
		dif_offsets[0x0b] = pb;
		dif_offsets[0x0c] = pc;
		dif_offsets[0x0d] = pd;
		dif_offsets[0x0e] = pe;
		dif_offsets[0x0f] = pf;
		for(i = 0; i < 16; i++) { dif_offsets[i+0x010] = dif_offsets[i] + p10; }
		// And now make 8 copies of the initial 32 elts:
		for(i = 0; i < 32; i++) { dif_offsets[i+0x020] = dif_offsets[i]; }
		for(i = 0; i < 32; i++) { dif_offsets[i+0x040] = dif_offsets[i]; }
		for(i = 0; i < 32; i++) { dif_offsets[i+0x060] = dif_offsets[i]; }
		for(i = 0; i < 32; i++) { dif_offsets[i+0x080] = dif_offsets[i]; }
		for(i = 0; i < 32; i++) { dif_offsets[i+0x0a0] = dif_offsets[i]; }
		for(i = 0; i < 32; i++) { dif_offsets[i+0x0c0] = dif_offsets[i]; }
		for(i = 0; i < 32; i++) { dif_offsets[i+0x0e0] = dif_offsets[i]; }
		for(i = 0; i < 32; i++) { dif_offsets[i+0x100] = dif_offsets[i]; }
	#else
	// dif_offsets are w.r.to a-array, need 9 distinct sets of these, one for each DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		// Set 0: [0,1,3,2,7,6,5,4,e,f,c,d,a,b,8,9 + p00],[c,d,f,e,8,9,b,a,5,4,6,7,1,0,2,3 + p10]
		i = 0;
		dif_offsets[0x00] =  0;		dif_offsets[0x10] = pc+p10;
		dif_offsets[0x01] = p1;		dif_offsets[0x11] = pd+p10;
		dif_offsets[0x02] = p3;		dif_offsets[0x12] = pf+p10;
		dif_offsets[0x03] = p2;		dif_offsets[0x13] = pe+p10;
		dif_offsets[0x04] = p7;		dif_offsets[0x14] = p8+p10;
		dif_offsets[0x05] = p6;		dif_offsets[0x15] = p9+p10;
		dif_offsets[0x06] = p5;		dif_offsets[0x16] = pb+p10;
		dif_offsets[0x07] = p4;		dif_offsets[0x17] = pa+p10;
		dif_offsets[0x08] = pe;		dif_offsets[0x18] = p5+p10;
		dif_offsets[0x09] = pf;		dif_offsets[0x19] = p4+p10;
		dif_offsets[0x0a] = pc;		dif_offsets[0x1a] = p6+p10;
		dif_offsets[0x0b] = pd;		dif_offsets[0x1b] = p7+p10;
		dif_offsets[0x0c] = pa;		dif_offsets[0x1c] = p1+p10;
		dif_offsets[0x0d] = pb;		dif_offsets[0x1d] =    p10;
		dif_offsets[0x0e] = p8;		dif_offsets[0x1e] = p2+p10;
		dif_offsets[0x0f] = p9;		dif_offsets[0x1f] = p3+p10;
		// Set 1: [3,2,1,0,5,4,6,7,c,d,f,e,8,9,b,a +p100],[f,e,d,c,b,a,9,8,6,7,4,5,2,3,0,1 +p110]
		i += 32;
		dif_offsets[i+0x00] = p3;		dif_offsets[i+0x10] = pf+p10;
		dif_offsets[i+0x01] = p2;		dif_offsets[i+0x11] = pe+p10;
		dif_offsets[i+0x02] = p1;		dif_offsets[i+0x12] = pd+p10;
		dif_offsets[i+0x03] =  0;		dif_offsets[i+0x13] = pc+p10;
		dif_offsets[i+0x04] = p5;		dif_offsets[i+0x14] = pb+p10;
		dif_offsets[i+0x05] = p4;		dif_offsets[i+0x15] = pa+p10;
		dif_offsets[i+0x06] = p6;		dif_offsets[i+0x16] = p9+p10;
		dif_offsets[i+0x07] = p7;		dif_offsets[i+0x17] = p8+p10;
		dif_offsets[i+0x08] = pc;		dif_offsets[i+0x18] = p6+p10;
		dif_offsets[i+0x09] = pd;		dif_offsets[i+0x19] = p7+p10;
		dif_offsets[i+0x0a] = pf;		dif_offsets[i+0x1a] = p4+p10;
		dif_offsets[i+0x0b] = pe;		dif_offsets[i+0x1b] = p5+p10;
		dif_offsets[i+0x0c] = p8;		dif_offsets[i+0x1c] = p2+p10;
		dif_offsets[i+0x0d] = p9;		dif_offsets[i+0x1d] = p3+p10;
		dif_offsets[i+0x0e] = pb;		dif_offsets[i+0x1e] =    p10;
		dif_offsets[i+0x0f] = pa;		dif_offsets[i+0x1f] = p1+p10;
		// Set 2: [1,0,2,3,6,7,4,5,f,e,d,c,b,a,9,8 + pb0],[3,2,1,0,5,4,6,7,c,d,f,e,8,9,b,a + pa0]
		i += 32;
		dif_offsets[i+0x00] = p1+p10;		dif_offsets[i+0x10] = p3;
		dif_offsets[i+0x01] =    p10;		dif_offsets[i+0x11] = p2;
		dif_offsets[i+0x02] = p2+p10;		dif_offsets[i+0x12] = p1;
		dif_offsets[i+0x03] = p3+p10;		dif_offsets[i+0x13] =  0;
		dif_offsets[i+0x04] = p6+p10;		dif_offsets[i+0x14] = p5;
		dif_offsets[i+0x05] = p7+p10;		dif_offsets[i+0x15] = p4;
		dif_offsets[i+0x06] = p4+p10;		dif_offsets[i+0x16] = p6;
		dif_offsets[i+0x07] = p5+p10;		dif_offsets[i+0x17] = p7;
		dif_offsets[i+0x08] = pf+p10;		dif_offsets[i+0x18] = pc;
		dif_offsets[i+0x09] = pe+p10;		dif_offsets[i+0x19] = pd;
		dif_offsets[i+0x0a] = pd+p10;		dif_offsets[i+0x1a] = pf;
		dif_offsets[i+0x0b] = pc+p10;		dif_offsets[i+0x1b] = pe;
		dif_offsets[i+0x0c] = pb+p10;		dif_offsets[i+0x1c] = p8;
		dif_offsets[i+0x0d] = pa+p10;		dif_offsets[i+0x1d] = p9;
		dif_offsets[i+0x0e] = p9+p10;		dif_offsets[i+0x1e] = pb;
		dif_offsets[i+0x0f] = p8+p10;		dif_offsets[i+0x1f] = pa;
		// Set 3: [a,b,8,9,c,d,f,e,3,2,1,0,5,4,6,7 + p40],[1,0,2,3,6,7,4,5,f,e,d,c,b,a,9,8 + p50]
		i += 32;
		dif_offsets[i+0x00] = pa;		dif_offsets[i+0x10] = p1+p10;
		dif_offsets[i+0x01] = pb;		dif_offsets[i+0x11] =    p10;
		dif_offsets[i+0x02] = p8;		dif_offsets[i+0x12] = p2+p10;
		dif_offsets[i+0x03] = p9;		dif_offsets[i+0x13] = p3+p10;
		dif_offsets[i+0x04] = pc;		dif_offsets[i+0x14] = p6+p10;
		dif_offsets[i+0x05] = pd;		dif_offsets[i+0x15] = p7+p10;
		dif_offsets[i+0x06] = pf;		dif_offsets[i+0x16] = p4+p10;
		dif_offsets[i+0x07] = pe;		dif_offsets[i+0x17] = p5+p10;
		dif_offsets[i+0x08] = p3;		dif_offsets[i+0x18] = pf+p10;
		dif_offsets[i+0x09] = p2;		dif_offsets[i+0x19] = pe+p10;
		dif_offsets[i+0x0a] = p1;		dif_offsets[i+0x1a] = pd+p10;
		dif_offsets[i+0x0b] =  0;		dif_offsets[i+0x1b] = pc+p10;
		dif_offsets[i+0x0c] = p5;		dif_offsets[i+0x1c] = pb+p10;
		dif_offsets[i+0x0d] = p4;		dif_offsets[i+0x1d] = pa+p10;
		dif_offsets[i+0x0e] = p6;		dif_offsets[i+0x1e] = p9+p10;
		dif_offsets[i+0x0f] = p7;		dif_offsets[i+0x1f] = p8+p10;
		// Set 4: [8,9,b,a,f,e,d,c,1,0,2,3,6,7,4,5 + pf0],[a,b,8,9,c,d,f,e,3,2,1,0,5,4,6,7 + pe0]
		i += 32;
		dif_offsets[i+0x00] = p8+p10;		dif_offsets[i+0x10] = pa;
		dif_offsets[i+0x01] = p9+p10;		dif_offsets[i+0x11] = pb;
		dif_offsets[i+0x02] = pb+p10;		dif_offsets[i+0x12] = p8;
		dif_offsets[i+0x03] = pa+p10;		dif_offsets[i+0x13] = p9;
		dif_offsets[i+0x04] = pf+p10;		dif_offsets[i+0x14] = pc;
		dif_offsets[i+0x05] = pe+p10;		dif_offsets[i+0x15] = pd;
		dif_offsets[i+0x06] = pd+p10;		dif_offsets[i+0x16] = pf;
		dif_offsets[i+0x07] = pc+p10;		dif_offsets[i+0x17] = pe;
		dif_offsets[i+0x08] = p1+p10;		dif_offsets[i+0x18] = p3;
		dif_offsets[i+0x09] =    p10;		dif_offsets[i+0x19] = p2;
		dif_offsets[i+0x0a] = p2+p10;		dif_offsets[i+0x1a] = p1;
		dif_offsets[i+0x0b] = p3+p10;		dif_offsets[i+0x1b] =  0;
		dif_offsets[i+0x0c] = p6+p10;		dif_offsets[i+0x1c] = p5;
		dif_offsets[i+0x0d] = p7+p10;		dif_offsets[i+0x1d] = p4;
		dif_offsets[i+0x0e] = p4+p10;		dif_offsets[i+0x1e] = p6;
		dif_offsets[i+0x0f] = p5+p10;		dif_offsets[i+0x1f] = p7;
		// Set 5: [7,6,5,4,3,2,1,0,a,b,8,9,c,d,f,e + p80],[8,9,b,a,f,e,d,c,1,0,2,3,6,7,4,5 + p90]
		i += 32;
		dif_offsets[i+0x00] = p7;		dif_offsets[i+0x10] = p8+p10;
		dif_offsets[i+0x01] = p6;		dif_offsets[i+0x11] = p9+p10;
		dif_offsets[i+0x02] = p5;		dif_offsets[i+0x12] = pb+p10;
		dif_offsets[i+0x03] = p4;		dif_offsets[i+0x13] = pa+p10;
		dif_offsets[i+0x04] = p3;		dif_offsets[i+0x14] = pf+p10;
		dif_offsets[i+0x05] = p2;		dif_offsets[i+0x15] = pe+p10;
		dif_offsets[i+0x06] = p1;		dif_offsets[i+0x16] = pd+p10;
		dif_offsets[i+0x07] =  0;		dif_offsets[i+0x17] = pc+p10;
		dif_offsets[i+0x08] = pa;		dif_offsets[i+0x18] = p1+p10;
		dif_offsets[i+0x09] = pb;		dif_offsets[i+0x19] =    p10;
		dif_offsets[i+0x0a] = p8;		dif_offsets[i+0x1a] = p2+p10;
		dif_offsets[i+0x0b] = p9;		dif_offsets[i+0x1b] = p3+p10;
		dif_offsets[i+0x0c] = pc;		dif_offsets[i+0x1c] = p6+p10;
		dif_offsets[i+0x0d] = pd;		dif_offsets[i+0x1d] = p7+p10;
		dif_offsets[i+0x0e] = pf;		dif_offsets[i+0x1e] = p4+p10;
		dif_offsets[i+0x0f] = pe;		dif_offsets[i+0x1f] = p5+p10;
		// Set 6: [5,4,6,7,1,0,2,3,8,9,b,a,f,e,d,c + p30],[7,6,5,4,3,2,1,0,a,b,8,9,c,d,f,e + p20]
		i += 32;
		dif_offsets[i+0x00] = p5+p10;		dif_offsets[i+0x10] = p7;
		dif_offsets[i+0x01] = p4+p10;		dif_offsets[i+0x11] = p6;
		dif_offsets[i+0x02] = p6+p10;		dif_offsets[i+0x12] = p5;
		dif_offsets[i+0x03] = p7+p10;		dif_offsets[i+0x13] = p4;
		dif_offsets[i+0x04] = p1+p10;		dif_offsets[i+0x14] = p3;
		dif_offsets[i+0x05] =    p10;		dif_offsets[i+0x15] = p2;
		dif_offsets[i+0x06] = p2+p10;		dif_offsets[i+0x16] = p1;
		dif_offsets[i+0x07] = p3+p10;		dif_offsets[i+0x17] =  0;
		dif_offsets[i+0x08] = p8+p10;		dif_offsets[i+0x18] = pa;
		dif_offsets[i+0x09] = p9+p10;		dif_offsets[i+0x19] = pb;
		dif_offsets[i+0x0a] = pb+p10;		dif_offsets[i+0x1a] = p8;
		dif_offsets[i+0x0b] = pa+p10;		dif_offsets[i+0x1b] = p9;
		dif_offsets[i+0x0c] = pf+p10;		dif_offsets[i+0x1c] = pc;
		dif_offsets[i+0x0d] = pe+p10;		dif_offsets[i+0x1d] = pd;
		dif_offsets[i+0x0e] = pd+p10;		dif_offsets[i+0x1e] = pf;
		dif_offsets[i+0x0f] = pc+p10;		dif_offsets[i+0x1f] = pe;
		// Set 7: [e,f,c,d,a,b,8,9,7,6,5,4,3,2,1,0 + pc0],[5,4,6,7,1,0,2,3,8,9,b,a,f,e,d,c + pd0]
		i += 32;
		dif_offsets[i+0x00] = pe;		dif_offsets[i+0x10] = p5+p10;
		dif_offsets[i+0x01] = pf;		dif_offsets[i+0x11] = p4+p10;
		dif_offsets[i+0x02] = pc;		dif_offsets[i+0x12] = p6+p10;
		dif_offsets[i+0x03] = pd;		dif_offsets[i+0x13] = p7+p10;
		dif_offsets[i+0x04] = pa;		dif_offsets[i+0x14] = p1+p10;
		dif_offsets[i+0x05] = pb;		dif_offsets[i+0x15] =    p10;
		dif_offsets[i+0x06] = p8;		dif_offsets[i+0x16] = p2+p10;
		dif_offsets[i+0x07] = p9;		dif_offsets[i+0x17] = p3+p10;
		dif_offsets[i+0x08] = p7;		dif_offsets[i+0x18] = p8+p10;
		dif_offsets[i+0x09] = p6;		dif_offsets[i+0x19] = p9+p10;
		dif_offsets[i+0x0a] = p5;		dif_offsets[i+0x1a] = pb+p10;
		dif_offsets[i+0x0b] = p4;		dif_offsets[i+0x1b] = pa+p10;
		dif_offsets[i+0x0c] = p3;		dif_offsets[i+0x1c] = pf+p10;
		dif_offsets[i+0x0d] = p2;		dif_offsets[i+0x1d] = pe+p10;
		dif_offsets[i+0x0e] = p1;		dif_offsets[i+0x1e] = pd+p10;
		dif_offsets[i+0x0f] =  0;		dif_offsets[i+0x1f] = pc+p10;
		// Set 8: [c,d,f,e,8,9,b,a,5,4,6,7,1,0,2,3 + p70],[e,f,c,d,a,b,8,9,7,6,5,4,3,2,1,0 + p60]
		i += 32;
		dif_offsets[i+0x00] = pc+p10;		dif_offsets[i+0x10] = pe;
		dif_offsets[i+0x01] = pd+p10;		dif_offsets[i+0x11] = pf;
		dif_offsets[i+0x02] = pf+p10;		dif_offsets[i+0x12] = pc;
		dif_offsets[i+0x03] = pe+p10;		dif_offsets[i+0x13] = pd;
		dif_offsets[i+0x04] = p8+p10;		dif_offsets[i+0x14] = pa;
		dif_offsets[i+0x05] = p9+p10;		dif_offsets[i+0x15] = pb;
		dif_offsets[i+0x06] = pb+p10;		dif_offsets[i+0x16] = p8;
		dif_offsets[i+0x07] = pa+p10;		dif_offsets[i+0x17] = p9;
		dif_offsets[i+0x08] = p5+p10;		dif_offsets[i+0x18] = p7;
		dif_offsets[i+0x09] = p4+p10;		dif_offsets[i+0x19] = p6;
		dif_offsets[i+0x0a] = p6+p10;		dif_offsets[i+0x1a] = p5;
		dif_offsets[i+0x0b] = p7+p10;		dif_offsets[i+0x1b] = p4;
		dif_offsets[i+0x0c] = p1+p10;		dif_offsets[i+0x1c] = p3;
		dif_offsets[i+0x0d] =    p10;		dif_offsets[i+0x1d] = p2;
		dif_offsets[i+0x0e] = p2+p10;		dif_offsets[i+0x1e] = p1;
		dif_offsets[i+0x0f] = p3+p10;		dif_offsets[i+0x1f] =  0;
	#endif
	}

/*...The radix-288 pass is here.	*/

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

	//...gather the needed data (288 64-bit complex) and do 32 radix-9 transforms:
	/*
	Twiddleless version arranges 32 sets of radix-9 DFT inputs as follows: 0 in upper left corner,
	decrement 32 (= 0x20) horizontally and 9 vertically, all indexing done (mod 288 = 0x120).

	Indexing in hex for clarity and using [evn|odd]0-8 notation in the rightmost column to flag reusable
	9-perms [in fact simple circular (0-8)-element shifts of evn = [000,100,0e0,0c0,0a0,080,060,040,020]
	and odd = [110,0f0,0d0,0b0,090,070,050,030,010]:

	DIF/DIT input-scramble array:		[000,100,...,020 = basic even-offset array] row	leftward-circ-shift count
										[110,0f0,...,010 = basic  odd-offset array] idx	of basic even|odd arrays:
	000,100,0e0,0c0,0a0,080,060,040,020 = 000,100,0e0,0c0,0a0,080,060,040,020 + p0	0	[evn0] + p0
	117,0f7,0d7,0b7,097,077,057,037,017 = 110,0f0,0d0,0b0,090,070,050,030,010 + p7	1	[odd0] + p7
	10e,0ee,0ce,0ae,08e,06e,04e,02e,00e = 100,0e0,0c0,0a0,080,060,040,020,000 + pe	2	[evn1] + pe
	105,0e5,0c5,0a5,085,065,045,025,005 = 100,0e0,0c0,0a0,080,060,040,020,000 + p5	3	[evn1] + p5
	0fc,0dc,0bc,09c,07c,05c,03c,01c,11c = 0f0,0d0,0b0,090,070,050,030,010,110 + pc	4	[odd1] + pc
	0f3,0d3,0b3,093,073,053,033,013,113 = 0f0,0d0,0b0,090,070,050,030,010,110 + p3	5	[odd1] + p3
	0ea,0ca,0aa,08a,06a,04a,02a,00a,10a = 0e0,0c0,0a0,080,060,040,020,000,100 + pa	6	[evn2] + pa
	0e1,0c1,0a1,081,061,041,021,001,101 = 0e0,0c0,0a0,080,060,040,020,000,100 + p1	7	[evn2] + p1
	0d8,0b8,098,078,058,038,018,118,0f8 = 0d0,0b0,090,070,050,030,010,110,0f0 + p8	8	[odd2] + p8
	0cf,0af,08f,06f,04f,02f,00f,10f,0ef = 0c0,0a0,080,060,040,020,000,100,0e0 + pf	9	[evn3] + pf
	0c6,0a6,086,066,046,026,006,106,0e6 = 0c0,0a0,080,060,040,020,000,100,0e0 + p6	a	[evn3] + p6
	0bd,09d,07d,05d,03d,01d,11d,0fd,0dd = 0b0,090,070,050,030,010,110,0f0,0d0 + pd	b	[odd3] + pd
	0b4,094,074,054,034,014,114,0f4,0d4 = 0b0,090,070,050,030,010,110,0f0,0d0 + p4	c	[odd3] + p4
	0ab,08b,06b,04b,02b,00b,10b,0eb,0cb = 0a0,080,060,040,020,000,100,0e0,0c0 + pb	d	[evn4] + pb
	0a2,082,062,042,022,002,102,0e2,0c2 = 0a0,080,060,040,020,000,100,0e0,0c0 + p2	e	[evn4] + p2
	099,079,059,039,019,119,0f9,0d9,0b9 = 090,070,050,030,010,110,0f0,0d0,0b0 + p9	f	[odd4] + p9
	090,070,050,030,010,110,0f0,0d0,0b0 = 090,070,050,030,010,110,0f0,0d0,0b0 + p0	10	[odd4] + p0	<<< p0-f pattern repeats here
	087,067,047,027,007,107,0e7,0c7,0a7 = 080,060,040,020,000,100,0e0,0c0,0a0 + p7	11	[evn5] + p7
	07e,05e,03e,01e,11e,0fe,0de,0be,09e = 070,050,030,010,110,0f0,0d0,0b0,090 + pe	12	[odd5] + pe
	075,055,035,015,115,0f5,0d5,0b5,095 = 070,050,030,010,110,0f0,0d0,0b0,090 + p5	13	[odd5] + p5
	06c,04c,02c,00c,10c,0ec,0cc,0ac,08c = 060,040,020,000,100,0e0,0c0,0a0,080 + pc	14	[evn6] + pc
	063,043,023,003,103,0e3,0c3,0a3,083 = 060,040,020,000,100,0e0,0c0,0a0,080 + p3	15	[evn6] + p3
	05a,03a,01a,11a,0fa,0da,0ba,09a,07a = 050,030,010,110,0f0,0d0,0b0,090,070 + pa	16	[odd6] + pa
	051,031,011,111,0f1,0d1,0b1,091,071 = 050,030,010,110,0f0,0d0,0b0,090,070 + p1	17	[odd6] + p1
	048,028,008,108,0e8,0c8,0a8,088,068 = 040,020,000,100,0e0,0c0,0a0,080,060 + p8	18	[evn7] + p8
	03f,01f,11f,0ff,0df,0bf,09f,07f,05f = 030,010,110,0f0,0d0,0b0,090,070,050 + pf	19	[odd7] + pf
	036,016,116,0f6,0d6,0b6,096,076,056 = 030,010,110,0f0,0d0,0b0,090,070,050 + p6	1a	[odd7] + p6
	02d,00d,10d,0ed,0cd,0ad,08d,06d,04d = 020,000,100,0e0,0c0,0a0,080,060,040 + pd	1b	[evn8] + pd
	024,004,104,0e4,0c4,0a4,084,064,044 = 020,000,100,0e0,0c0,0a0,080,060,040 + p4	1c	[evn8] + p4
	01b,11b,0fb,0db,0bb,09b,07b,05b,03b = 010,110,0f0,0d0,0b0,090,070,050,030 + pb	1d	[odd8] + pb
	012,112,0f2,0d2,0b2,092,072,052,032 = 010,110,0f0,0d0,0b0,090,070,050,030 + p2	1e	[odd8] + p2
	009,109,0e9,0c9,0a9,089,069,049,029 = 000,100,0e0,0c0,0a0,080,060,040,020 + p9	1f	[evn0] + p9
	*/
		tptr = t;
		for(i = 0; i < 32; i++) {
			int k = dif_p20_lo_offset[i];
			// Extract index (in [0-8]) into circ-shift array used for high parts of p-mults. The [0-8] value is
			// in low 4 bits of k; the "which length-17 half of the dif_p20_cperms array?" selector is via (k < 0):
			int ic = ((-(k < 0)) & 17)	// +/- sign on k puts us into lower/upper half of the cshift array (base index 0/17)
						+ (k & 0xf);	// ...and low 4 bits give the element index w.r.to the array-half in question.
			int k0 = dif_p20_cperms[ic], k1 = dif_p20_cperms[ic+1], k2 = dif_p20_cperms[ic+2], k3 = dif_p20_cperms[ic+3], k4 = dif_p20_cperms[ic+4], k5 = dif_p20_cperms[ic+5], k6 = dif_p20_cperms[ic+6], k7 = dif_p20_cperms[ic+7], k8 = dif_p20_cperms[ic+8];
			// Extract Low part, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-9 DFTs:
			k = (k & 0x7fffffff) >> 4;
			jt = j1+k; jp = j2+k;
			RADIX_09_DIF(
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],
				tptr->re,tptr->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x100)->re,(tptr+0x100)->im,
				rt,it,re
			);	tptr++;
		}

	/*...and now do 9 radix-32 transforms;
	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF]
	to properly permute the outputs of the radix-32 DFTs to the required ordering, which in terms of our p-offsets is

	000,001,003,002,007,006,005,004,00e,00f,00c,00d,00a,00b,008,009,01c,01d,01f,01e,018,019,01b,01a,015,014,016,017,011,010,012,013
	103,102,101,100,105,104,106,107,10c,10d,10f,10e,108,109,10b,10a,11f,11e,11d,11c,11b,11a,119,118,116,117,114,115,112,113,110,111
	0b1,0b0,0b2,0b3,0b6,0b7,0b4,0b5,0bf,0be,0bd,0bc,0bb,0ba,0b9,0b8,0a3,0a2,0a1,0a0,0a5,0a4,0a6,0a7,0ac,0ad,0af,0ae,0a8,0a9,0ab,0aa
	04a,04b,048,049,04c,04d,04f,04e,043,042,041,040,045,044,046,047,051,050,052,053,056,057,054,055,05f,05e,05d,05c,05b,05a,059,058
	0f8,0f9,0fb,0fa,0ff,0fe,0fd,0fc,0f1,0f0,0f2,0f3,0f6,0f7,0f4,0f5,0ea,0eb,0e8,0e9,0ec,0ed,0ef,0ee,0e3,0e2,0e1,0e0,0e5,0e4,0e6,0e7
	087,086,085,084,083,082,081,080,08a,08b,088,089,08c,08d,08f,08e,098,099,09b,09a,09f,09e,09d,09c,091,090,092,093,096,097,094,095
	035,034,036,037,031,030,032,033,038,039,03b,03a,03f,03e,03d,03c,027,026,025,024,023,022,021,020,02a,02b,028,029,02c,02d,02f,02e
	0ce,0cf,0cc,0cd,0ca,0cb,0c8,0c9,0c7,0c6,0c5,0c4,0c3,0c2,0c1,0c0,0d5,0d4,0d6,0d7,0d1,0d0,0d2,0d3,0d8,0d9,0db,0da,0df,0de,0dd,0dc
	07c,07d,07f,07e,078,079,07b,07a,075,074,076,077,071,070,072,073,06e,06f,06c,06d,06a,06b,068,069,067,066,065,064,063,062,061,060
=
	[0,1,3,2,7,6,5,4,e,f,c,d,a,b,8,9 + p00],[c,d,f,e,8,9,b,a,5,4,6,7,1,0,2,3 + p10]
	[3,2,1,0,5,4,6,7,c,d,f,e,8,9,b,a +p100],[f,e,d,c,b,a,9,8,6,7,4,5,2,3,0,1 +p110]
	[1,0,2,3,6,7,4,5,f,e,d,c,b,a,9,8 + pb0],[3,2,1,0,5,4,6,7,c,d,f,e,8,9,b,a + pa0]
	[a,b,8,9,c,d,f,e,3,2,1,0,5,4,6,7 + p40],[1,0,2,3,6,7,4,5,f,e,d,c,b,a,9,8 + p50]
	[8,9,b,a,f,e,d,c,1,0,2,3,6,7,4,5 + pf0],[a,b,8,9,c,d,f,e,3,2,1,0,5,4,6,7 + pe0]
	[7,6,5,4,3,2,1,0,a,b,8,9,c,d,f,e + p80],[8,9,b,a,f,e,d,c,1,0,2,3,6,7,4,5 + p90]
	[5,4,6,7,1,0,2,3,8,9,b,a,f,e,d,c + p30],[7,6,5,4,3,2,1,0,a,b,8,9,c,d,f,e + p20]
	[e,f,c,d,a,b,8,9,7,6,5,4,3,2,1,0 + pc0],[5,4,6,7,1,0,2,3,8,9,b,a,f,e,d,c + pd0]
	[c,d,f,e,8,9,b,a,5,4,6,7,1,0,2,3 + p70],[e,f,c,d,a,b,8,9,7,6,5,4,3,2,1,0 + p60]
	*/
		tptr = t;
		for(i = 0; i < ODD_RADIX; i++) {
			jt = j1+dif_phi[i]; RADIX_32_DIF((double *)tptr,t_offsets,1, (a+jt),dif_offsets+(i<<5),RE_IM_STRIDE);	tptr += 32;
		}
	}
}

/***************/

void radix288_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-288 complex DIT FFT pass on the data in the length-N real vector A.
*/
	const double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
	int i,j,j1,j2,jt,jp;
	// p-indexing is hexadecimal here:
	static int NDIVR,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,p100,p110, first_entry=TRUE;
	static int t_offsets[32], dit_offsets[RADIX];
	// Need storage for 2 circular-shifts perms of a basic 9-vector, with shift count in [0,8] that means 2*17 elts:
	static int dit_p20_cperms[34], dit_p20_lo_offset[32], dit_phi[ODD_RADIX];
	double rt,it,re;
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
		p1 = NDIVR;		p10 = NDIVR<<4;		p100 = NDIVR<<8;
		p2 = p1 + p1;	p20 = p10 + p10;	p110 = p100 + p10;
		p3 = p2 + p1;	p30 = p20 + p10;
		p4 = p3 + p1;	p40 = p30 + p10;
		p5 = p4 + p1;	p50 = p40 + p10;
		p6 = p5 + p1;	p60 = p50 + p10;
		p7 = p6 + p1;	p70 = p60 + p10;
		p8 = p7 + p1;	p80 = p70 + p10;
		p9 = p8 + p1;	p90 = p80 + p10;
		pa = p9 + p1;	pa0 = p90 + p10;
		pb = pa + p1;	pb0 = pa0 + p10;
		pc = pb + p1;	pc0 = pb0 + p10;
		pd = pc + p1;	pd0 = pc0 + p10;
		pe = pd + p1;	pe0 = pd0 + p10;
		pf = pe + p1;	pf0 = pe0 + p10;

		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );		p50 += ( (p50 >> DAT_BITS) << PAD_BITS );
		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );		p70 += ( (p70 >> DAT_BITS) << PAD_BITS );
		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pa += ( (pa >> DAT_BITS) << PAD_BITS );		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pb += ( (pb >> DAT_BITS) << PAD_BITS );		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		pc += ( (pc >> DAT_BITS) << PAD_BITS );		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pd += ( (pd >> DAT_BITS) << PAD_BITS );		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		pe += ( (pe >> DAT_BITS) << PAD_BITS );		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		pf += ( (pf >> DAT_BITS) << PAD_BITS );		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
													p100 += ( (p100 >> DAT_BITS) << PAD_BITS );
													p110 += ( (p110 >> DAT_BITS) << PAD_BITS );
		i = 0;
		dit_phi[i++] =   0;
		dit_phi[i++] = p20;
		dit_phi[i++] = p40;
		dit_phi[i++] = p80;
		dit_phi[i++] = pa0;
		dit_phi[i++] = p60;
		dit_phi[i++] =p100;
		dit_phi[i++] = pc0;
		dit_phi[i++] = pe0;

		// Set array offsets for radix-32 DFT in/outputs:
		// t_offsets w.r.to: t-array, same for all 9 DFTs:
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

		// Init storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,8] that means 2*17
		i = 0;
		// Even multiples of p10 cshift array: evn0 := p10*[00,02,04,06,08,0a,0c,0e,10]
		dit_p20_cperms[i++] =   0;
		dit_p20_cperms[i++] = p20;
		dit_p20_cperms[i++] = p40;
		dit_p20_cperms[i++] = p60;
		dit_p20_cperms[i++] = p80;
		dit_p20_cperms[i++] = pa0;
		dit_p20_cperms[i++] = pc0;
		dit_p20_cperms[i++] = pe0;
		dit_p20_cperms[i++] = p100;
		while(i < 2*ODD_RADIX-1) {
			dit_p20_cperms[i] = dit_p20_cperms[i - ODD_RADIX]; ++i;
		}
		// Odd multiples of p10 cshift array: odd0 := p10*[03,05,07,09,0b,0d,0f,11,01]
		dit_p20_cperms[i++] = p30;
		dit_p20_cperms[i++] = p50;
		dit_p20_cperms[i++] = p70;
		dit_p20_cperms[i++] = p90;
		dit_p20_cperms[i++] = pb0;
		dit_p20_cperms[i++] = pd0;
		dit_p20_cperms[i++] = pf0;
		dit_p20_cperms[i++] = p110;
		dit_p20_cperms[i++] = p10;
		while(i < 4*ODD_RADIX-2) {
			dit_p20_cperms[i] = dit_p20_cperms[i - ODD_RADIX]; ++i;
		}

		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-9 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-8 index, bits <4:30> for the p0-f:
		i = 0;
		dit_p20_lo_offset[i++] = (( 0 << 4) + 0);
		dit_p20_lo_offset[i++] = ((pf << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pe << 4) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pd << 4) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pc << 4) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pb << 4) + 8) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pa << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p9 << 4) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p8 << 4) + 5) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p7 << 4) + 7) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p6 << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p5 << 4) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p4 << 4) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p3 << 4) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p2 << 4) + 8) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p1 << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = (( 0 << 4) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pf << 4) + 6);
		dit_p20_lo_offset[i++] = ((pe << 4) + 8);
		dit_p20_lo_offset[i++] = ((pd << 4) + 1);
		dit_p20_lo_offset[i++] = ((pc << 4) + 3);
		dit_p20_lo_offset[i++] = ((pb << 4) + 5);
		dit_p20_lo_offset[i++] = ((pa << 4) + 7);
		dit_p20_lo_offset[i++] = ((p9 << 4) + 0);
		dit_p20_lo_offset[i++] = ((p8 << 4) + 2);
		dit_p20_lo_offset[i++] = ((p7 << 4) + 4);
		dit_p20_lo_offset[i++] = ((p6 << 4) + 6);
		dit_p20_lo_offset[i++] = ((p5 << 4) + 8);
		dit_p20_lo_offset[i++] = ((p4 << 4) + 1);
		dit_p20_lo_offset[i++] = ((p3 << 4) + 3);
		dit_p20_lo_offset[i++] = ((p2 << 4) + 5);
		dit_p20_lo_offset[i++] = ((p1 << 4) + 7);

	// dit_offsets are w.r.to a-array, need 9 distinct sets of these, one for each radix-32 DFT.
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
		// Set 1: [5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p30],[5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p20]
		i += 32;
		dit_offsets[i+0x00] = p5+p10;		dit_offsets[i+0x10] = p5;
		dit_offsets[i+0x01] = p4+p10;		dit_offsets[i+0x11] = p4;
		dit_offsets[i+0x02] = p6+p10;		dit_offsets[i+0x12] = p6;
		dit_offsets[i+0x03] = p7+p10;		dit_offsets[i+0x13] = p7;
		dit_offsets[i+0x04] = p1+p10;		dit_offsets[i+0x14] = p1;
		dit_offsets[i+0x05] =    p10;		dit_offsets[i+0x15] =  0;
		dit_offsets[i+0x06] = p2+p10;		dit_offsets[i+0x16] = p2;
		dit_offsets[i+0x07] = p3+p10;		dit_offsets[i+0x17] = p3;
		dit_offsets[i+0x08] = p9+p10;		dit_offsets[i+0x18] = p9;
		dit_offsets[i+0x09] = p8+p10;		dit_offsets[i+0x19] = p8;
		dit_offsets[i+0x0a] = pa+p10;		dit_offsets[i+0x1a] = pa;
		dit_offsets[i+0x0b] = pb+p10;		dit_offsets[i+0x1b] = pb;
		dit_offsets[i+0x0c] = pe+p10;		dit_offsets[i+0x1c] = pe;
		dit_offsets[i+0x0d] = pf+p10;		dit_offsets[i+0x1d] = pf;
		dit_offsets[i+0x0e] = pc+p10;		dit_offsets[i+0x1e] = pc;
		dit_offsets[i+0x0f] = pd+p10;		dit_offsets[i+0x1f] = pd;
		// Set 2: [a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p40],[2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p50]
		i += 32;
		dit_offsets[i+0x00] = pa;		dit_offsets[i+0x10] = p2+p10;
		dit_offsets[i+0x01] = pb;		dit_offsets[i+0x11] = p3+p10;
		dit_offsets[i+0x02] = p8;		dit_offsets[i+0x12] =    p10;
		dit_offsets[i+0x03] = p9;		dit_offsets[i+0x13] = p1+p10;
		dit_offsets[i+0x04] = pc;		dit_offsets[i+0x14] = p4+p10;
		dit_offsets[i+0x05] = pd;		dit_offsets[i+0x15] = p5+p10;
		dit_offsets[i+0x06] = pf;		dit_offsets[i+0x16] = p7+p10;
		dit_offsets[i+0x07] = pe;		dit_offsets[i+0x17] = p6+p10;
		dit_offsets[i+0x08] = p2;		dit_offsets[i+0x18] = pc+p10;
		dit_offsets[i+0x09] = p3;		dit_offsets[i+0x19] = pd+p10;
		dit_offsets[i+0x0a] =  0;		dit_offsets[i+0x1a] = pf+p10;
		dit_offsets[i+0x0b] = p1;		dit_offsets[i+0x1b] = pe+p10;
		dit_offsets[i+0x0c] = p4;		dit_offsets[i+0x1c] = p8+p10;
		dit_offsets[i+0x0d] = p5;		dit_offsets[i+0x1d] = p9+p10;
		dit_offsets[i+0x0e] = p7;		dit_offsets[i+0x1e] = pb+p10;
		dit_offsets[i+0x0f] = p6;		dit_offsets[i+0x1f] = pa+p10;
		// Set 3: [7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + p80],[b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p90]
		i += 32;
		dit_offsets[i+0x00] = p7;		dit_offsets[i+0x10] = pb+p10;
		dit_offsets[i+0x01] = p6;		dit_offsets[i+0x11] = pa+p10;
		dit_offsets[i+0x02] = p5;		dit_offsets[i+0x12] = p9+p10;
		dit_offsets[i+0x03] = p4;		dit_offsets[i+0x13] = p8+p10;
		dit_offsets[i+0x04] = p3;		dit_offsets[i+0x14] = pd+p10;
		dit_offsets[i+0x05] = p2;		dit_offsets[i+0x15] = pc+p10;
		dit_offsets[i+0x06] = p1;		dit_offsets[i+0x16] = pe+p10;
		dit_offsets[i+0x07] =  0;		dit_offsets[i+0x17] = pf+p10;
		dit_offsets[i+0x08] = pb;		dit_offsets[i+0x18] = p3+p10;
		dit_offsets[i+0x09] = pa;		dit_offsets[i+0x19] = p2+p10;
		dit_offsets[i+0x0a] = p9;		dit_offsets[i+0x1a] = p1+p10;
		dit_offsets[i+0x0b] = p8;		dit_offsets[i+0x1b] =    p10;
		dit_offsets[i+0x0c] = pd;		dit_offsets[i+0x1c] = p5+p10;
		dit_offsets[i+0x0d] = pc;		dit_offsets[i+0x1d] = p4+p10;
		dit_offsets[i+0x0e] = pe;		dit_offsets[i+0x1e] = p6+p10;
		dit_offsets[i+0x0f] = pf;		dit_offsets[i+0x1f] = p7+p10;
		// Set 4: [1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pb0],[1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pa0]
		i += 32;
		dit_offsets[i+0x00] = p1+p10;		dit_offsets[i+0x10] = p1;
		dit_offsets[i+0x01] =    p10;		dit_offsets[i+0x11] =  0;
		dit_offsets[i+0x02] = p2+p10;		dit_offsets[i+0x12] = p2;
		dit_offsets[i+0x03] = p3+p10;		dit_offsets[i+0x13] = p3;
		dit_offsets[i+0x04] = p6+p10;		dit_offsets[i+0x14] = p6;
		dit_offsets[i+0x05] = p7+p10;		dit_offsets[i+0x15] = p7;
		dit_offsets[i+0x06] = p4+p10;		dit_offsets[i+0x16] = p4;
		dit_offsets[i+0x07] = p5+p10;		dit_offsets[i+0x17] = p5;
		dit_offsets[i+0x08] = pe+p10;		dit_offsets[i+0x18] = pe;
		dit_offsets[i+0x09] = pf+p10;		dit_offsets[i+0x19] = pf;
		dit_offsets[i+0x0a] = pc+p10;		dit_offsets[i+0x1a] = pc;
		dit_offsets[i+0x0b] = pd+p10;		dit_offsets[i+0x1b] = pd;
		dit_offsets[i+0x0c] = pa+p10;		dit_offsets[i+0x1c] = pa;
		dit_offsets[i+0x0d] = pb+p10;		dit_offsets[i+0x1d] = pb;
		dit_offsets[i+0x0e] = p8+p10;		dit_offsets[i+0x1e] = p8;
		dit_offsets[i+0x0f] = p9+p10;		dit_offsets[i+0x1f] = p9;
		// Set 5: [c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p70],[c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p60]
		i += 32;
		dit_offsets[i+0x00] = pc+p10;		dit_offsets[i+0x10] = pc;
		dit_offsets[i+0x01] = pd+p10;		dit_offsets[i+0x11] = pd;
		dit_offsets[i+0x02] = pf+p10;		dit_offsets[i+0x12] = pf;
		dit_offsets[i+0x03] = pe+p10;		dit_offsets[i+0x13] = pe;
		dit_offsets[i+0x04] = p8+p10;		dit_offsets[i+0x14] = p8;
		dit_offsets[i+0x05] = p9+p10;		dit_offsets[i+0x15] = p9;
		dit_offsets[i+0x06] = pb+p10;		dit_offsets[i+0x16] = pb;
		dit_offsets[i+0x07] = pa+p10;		dit_offsets[i+0x17] = pa;
		dit_offsets[i+0x08] = p4+p10;		dit_offsets[i+0x18] = p4;
		dit_offsets[i+0x09] = p5+p10;		dit_offsets[i+0x19] = p5;
		dit_offsets[i+0x0a] = p7+p10;		dit_offsets[i+0x1a] = p7;
		dit_offsets[i+0x0b] = p6+p10;		dit_offsets[i+0x1b] = p6;
		dit_offsets[i+0x0c] =    p10;		dit_offsets[i+0x1c] =  0;
		dit_offsets[i+0x0d] = p1+p10;		dit_offsets[i+0x1d] = p1;
		dit_offsets[i+0x0e] = p3+p10;		dit_offsets[i+0x1e] = p3;
		dit_offsets[i+0x0f] = p2+p10;		dit_offsets[i+0x1f] = p2;
		// Set 6: [3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b +p100],[d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 +p110]
		i += 32;
		dit_offsets[i+0x00] = p3;		dit_offsets[i+0x10] = pd+p10;
		dit_offsets[i+0x01] = p2;		dit_offsets[i+0x11] = pc+p10;
		dit_offsets[i+0x02] = p1;		dit_offsets[i+0x12] = pe+p10;
		dit_offsets[i+0x03] =  0;		dit_offsets[i+0x13] = pf+p10;
		dit_offsets[i+0x04] = p5;		dit_offsets[i+0x14] = p9+p10;
		dit_offsets[i+0x05] = p4;		dit_offsets[i+0x15] = p8+p10;
		dit_offsets[i+0x06] = p6;		dit_offsets[i+0x16] = pa+p10;
		dit_offsets[i+0x07] = p7;		dit_offsets[i+0x17] = pb+p10;
		dit_offsets[i+0x08] = pd;		dit_offsets[i+0x18] = p5+p10;
		dit_offsets[i+0x09] = pc;		dit_offsets[i+0x19] = p4+p10;
		dit_offsets[i+0x0a] = pe;		dit_offsets[i+0x1a] = p6+p10;
		dit_offsets[i+0x0b] = pf;		dit_offsets[i+0x1b] = p7+p10;
		dit_offsets[i+0x0c] = p9;		dit_offsets[i+0x1c] = p1+p10;
		dit_offsets[i+0x0d] = p8;		dit_offsets[i+0x1d] =    p10;
		dit_offsets[i+0x0e] = pa;		dit_offsets[i+0x1e] = p2+p10;
		dit_offsets[i+0x0f] = pb;		dit_offsets[i+0x1f] = p3+p10;
		// Set 7: [e,f,c,d,a,b,8,9,6,7,4,5,2,3,0,1 + pc0],[6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pd0]
		i += 32;
		dit_offsets[i+0x00] = pe;		dit_offsets[i+0x10] = p6+p10;
		dit_offsets[i+0x01] = pf;		dit_offsets[i+0x11] = p7+p10;
		dit_offsets[i+0x02] = pc;		dit_offsets[i+0x12] = p4+p10;
		dit_offsets[i+0x03] = pd;		dit_offsets[i+0x13] = p5+p10;
		dit_offsets[i+0x04] = pa;		dit_offsets[i+0x14] = p2+p10;
		dit_offsets[i+0x05] = pb;		dit_offsets[i+0x15] = p3+p10;
		dit_offsets[i+0x06] = p8;		dit_offsets[i+0x16] =    p10;
		dit_offsets[i+0x07] = p9;		dit_offsets[i+0x17] = p1+p10;
		dit_offsets[i+0x08] = p6;		dit_offsets[i+0x18] = pa+p10;
		dit_offsets[i+0x09] = p7;		dit_offsets[i+0x19] = pb+p10;
		dit_offsets[i+0x0a] = p4;		dit_offsets[i+0x1a] = p8+p10;
		dit_offsets[i+0x0b] = p5;		dit_offsets[i+0x1b] = p9+p10;
		dit_offsets[i+0x0c] = p2;		dit_offsets[i+0x1c] = pc+p10;
		dit_offsets[i+0x0d] = p3;		dit_offsets[i+0x1d] = pd+p10;
		dit_offsets[i+0x0e] =  0;		dit_offsets[i+0x1e] = pf+p10;
		dit_offsets[i+0x0f] = p1;		dit_offsets[i+0x1f] = pe+p10;
		// Set 8: [8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + pf0],[8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + pe0]
		i += 32;
		dit_offsets[i+0x00] = p8+p10;		dit_offsets[i+0x10] = p8;
		dit_offsets[i+0x01] = p9+p10;		dit_offsets[i+0x11] = p9;
		dit_offsets[i+0x02] = pb+p10;		dit_offsets[i+0x12] = pb;
		dit_offsets[i+0x03] = pa+p10;		dit_offsets[i+0x13] = pa;
		dit_offsets[i+0x04] = pf+p10;		dit_offsets[i+0x14] = pf;
		dit_offsets[i+0x05] = pe+p10;		dit_offsets[i+0x15] = pe;
		dit_offsets[i+0x06] = pd+p10;		dit_offsets[i+0x16] = pd;
		dit_offsets[i+0x07] = pc+p10;		dit_offsets[i+0x17] = pc;
		dit_offsets[i+0x08] =    p10;		dit_offsets[i+0x18] =  0;
		dit_offsets[i+0x09] = p1+p10;		dit_offsets[i+0x19] = p1;
		dit_offsets[i+0x0a] = p3+p10;		dit_offsets[i+0x1a] = p3;
		dit_offsets[i+0x0b] = p2+p10;		dit_offsets[i+0x1b] = p2;
		dit_offsets[i+0x0c] = p7+p10;		dit_offsets[i+0x1c] = p7;
		dit_offsets[i+0x0d] = p6+p10;		dit_offsets[i+0x1d] = p6;
		dit_offsets[i+0x0e] = p5+p10;		dit_offsets[i+0x1e] = p5;
		dit_offsets[i+0x0f] = p4+p10;		dit_offsets[i+0x1f] = p4;
	}

/*...The radix-288 pass is here.	*/

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
	/*
	Twiddleless version uses same linear-index-vector-form permutation as in DIF -
	Remember, inputs to DIT are bit-reversed, so use output of test_fft_radix() with
	TTYPE=0 to auto-generate needed input-index permutation:

	Combined DIT input-scramble array =
		000,001,003,002,007,006,005,004,00f,00e,00d,00c,00b,00a,009,008,01f,01e,01d,01c,01b,01a,019,018,017,016,015,014,013,012,011,010,
		035,034,036,037,031,030,032,033,039,038,03a,03b,03e,03f,03c,03d,025,024,026,027,021,020,022,023,029,028,02a,02b,02e,02f,02c,02d,
		04a,04b,048,049,04c,04d,04f,04e,042,043,040,041,044,045,047,046,052,053,050,051,054,055,057,056,05c,05d,05f,05e,058,059,05b,05a,
		087,086,085,084,083,082,081,080,08b,08a,089,088,08d,08c,08e,08f,09b,09a,099,098,09d,09c,09e,09f,093,092,091,090,095,094,096,097,
		0b1,0b0,0b2,0b3,0b6,0b7,0b4,0b5,0be,0bf,0bc,0bd,0ba,0bb,0b8,0b9,0a1,0a0,0a2,0a3,0a6,0a7,0a4,0a5,0ae,0af,0ac,0ad,0aa,0ab,0a8,0a9,
		07c,07d,07f,07e,078,079,07b,07a,074,075,077,076,070,071,073,072,06c,06d,06f,06e,068,069,06b,06a,064,065,067,066,060,061,063,062,
		103,102,101,100,105,104,106,107,10d,10c,10e,10f,109,108,10a,10b,11d,11c,11e,11f,119,118,11a,11b,115,114,116,117,111,110,112,113,
		0ce,0cf,0cc,0cd,0ca,0cb,0c8,0c9,0c6,0c7,0c4,0c5,0c2,0c3,0c0,0c1,0d6,0d7,0d4,0d5,0d2,0d3,0d0,0d1,0da,0db,0d8,0d9,0dc,0dd,0df,0de,
		0f8,0f9,0fb,0fa,0ff,0fe,0fd,0fc,0f0,0f1,0f3,0f2,0f7,0f6,0f5,0f4,0e8,0e9,0eb,0ea,0ef,0ee,0ed,0ec,0e0,0e1,0e3,0e2,0e7,0e6,0e5,0e4,
	=
		[0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00],[f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10]
		[5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p30],[5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p20]
		[a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p40],[2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p50]
		[7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + p80],[b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p90]
		[1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pb0],[1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pa0]
		[c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p70],[c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p60]
		[3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b +p100],[d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 +p110]
		[e,f,c,d,a,b,8,9,6,7,4,5,2,3,0,1 + pc0],[6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pd0]
		[8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + pf0],[8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + pe0]
	*/
	//...gather the needed data (288 64-bit complex) and do 9 radix-32 transforms:
		tptr = t;
		for(i = 0; i < ODD_RADIX; i++) {
			jt = j1+dit_phi[i]; RADIX_32_DIT((a+jt),dit_offsets+(i<<5),RE_IM_STRIDE, (double *)tptr,t_offsets,1);	tptr += 32;
		}

	/*...and now do 32 radix-9 transforms, with the columns of t*[r,i] output pairs in the above 9x radix-32 set now acting as input rows.
	Since our first-look oindex ordering was +p0-8 for each radix-9 and incrementing += p9 between those DFTs,
	simply break resulting o-index perm generated by the test_fft_radix code into 9-element rows.
	Indexing in hex for clarity and using [evn|odd]0-8 notation in the rightmost column to flag reusable 9-perms, in fact
	simple circular <<(0-8) shifts of the basic patterns:
		evn0 := p10*[00,02,04,06,08,0a,0c,0e,10]
		odd0 := p10*[03,05,07,09,0b,0d,0f,11,01]

	Required output index permutation =
		000,020,040,060,080,0a0,0c0,0e0,100 + p0		[evn0] + p0
		03f,05f,07f,09f,0bf,0df,0ff,11f,01f + pf		[odd0] + pf
		07e,09e,0be,0de,0fe,11e,01e,03e,05e + pe		[odd2] + pe
		0bd,0dd,0fd,11d,01d,03d,05d,07d,09d + pd		[odd4] + pd
		0fc,11c,01c,03c,05c,07c,09c,0bc,0dc + pc		[odd6] + pc
		01b,03b,05b,07b,09b,0bb,0db,0fb,11b + pb		[odd8] + pb
		05a,07a,09a,0ba,0da,0fa,11a,01a,03a + pa		[odd1] + pa
		099,0b9,0d9,0f9,119,019,039,059,079 + p9		[odd3] + p9
		0d8,0f8,118,018,038,058,078,098,0b8 + p8		[odd5] + p8
		117,017,037,057,077,097,0b7,0d7,0f7 + p7		[odd7] + p7
		036,056,076,096,0b6,0d6,0f6,116,016 + p6		[odd0] + p6
		075,095,0b5,0d5,0f5,115,015,035,055 + p5		[odd2] + p5
		0b4,0d4,0f4,114,014,034,054,074,094 + p4		[odd4] + p4
		0f3,113,013,033,053,073,093,0b3,0d3 + p3		[odd6] + p3
		012,032,052,072,092,0b2,0d2,0f2,112 + p2		[odd8] + p2
		051,071,091,0b1,0d1,0f1,111,011,031 + p1		[odd1] + p1
		090,0b0,0d0,0f0,110,010,030,050,070 + p0	=	[odd3] + p0
		0cf,0ef,10f,00f,02f,04f,06f,08f,0af + pf		[evn6] + pf
		10e,00e,02e,04e,06e,08e,0ae,0ce,0ee + pe		[evn8] + pe
		02d,04d,06d,08d,0ad,0cd,0ed,10d,00d + pd		[evn1] + pd
		06c,08c,0ac,0cc,0ec,10c,00c,02c,04c + pc		[evn3] + pc
		0ab,0cb,0eb,10b,00b,02b,04b,06b,08b + pb		[evn5] + pb
		0ea,10a,00a,02a,04a,06a,08a,0aa,0ca + pa		[evn7] + pa
		009,029,049,069,089,0a9,0c9,0e9,109 + p9		[evn0] + p9
		048,068,088,0a8,0c8,0e8,108,008,028 + p8		[evn2] + p8
		087,0a7,0c7,0e7,107,007,027,047,067 + p7		[evn4] + p7
		0c6,0e6,106,006,026,046,066,086,0a6 + p6		[evn6] + p6
		105,005,025,045,065,085,0a5,0c5,0e5 + p5		[evn8] + p5
		024,044,064,084,0a4,0c4,0e4,104,004 + p4		[evn1] + p4
		063,083,0a3,0c3,0e3,103,003,023,043 + p3		[evn3] + p3
		0a2,0c2,0e2,102,002,022,042,062,082 + p2		[evn5] + p2
		0e1,101,001,021,041,061,081,0a1,0c1 + p1		[evn7] + p1
	*/
		tptr = t;
		for(i = 0; i < 32; i++) {
		#if 0	// First-look linear indexing used to allow test_fft_radix to extract the needed operm:
			int k = i*18;
			int k0=0,k1=2,k2=4,k3=6,k4=8,k5=10,k6=12,k7=14,k8=16;
		#else
			int k = dit_p20_lo_offset[i];
			// Extract index (in [0-8]) into circ-shift array used for high parts of p-mults. The [0-8] value is
			// in low 4 bits of k; the "which length-17 half of the dit_p20_cperms array?" selector is via (k < 0):
			int ic = ((-(k < 0)) & 17)	// +/- sign on k puts us into lower/upper half of the cshift array (base index 0/17)
						+ (k & 0xf);	// ...and low 4 bits give the element index w.r.to the array-half in question.
			int k0 = dit_p20_cperms[ic], k1 = dit_p20_cperms[ic+1], k2 = dit_p20_cperms[ic+2], k3 = dit_p20_cperms[ic+3], k4 = dit_p20_cperms[ic+4], k5 = dit_p20_cperms[ic+5], k6 = dit_p20_cperms[ic+6], k7 = dit_p20_cperms[ic+7], k8 = dit_p20_cperms[ic+8];
			// Extract Low part, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-9 DFTs:
			k = (k & 0x7fffffff) >> 4;
		#endif
			jt = j1+k; jp = j2+k;
			RADIX_09_DIT(
				tptr->re,tptr->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x100)->re,(tptr+0x100)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],
				rt,it,re
			);	tptr++;
		}
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy288_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr;
		const int pfetch_dist = PFETCH_DIST;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
				,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,p100,p110;
	// Shared DIF+DIT:
		double rt,it,re;
		int poff[RADIX>>2];	// Store [RADIX/4] mults of p04 offset for loop control
		int t_offsets[32];
		// Need storage for 2 circular-shifts perms of a basic 9-vector, with shift count in [0,8] that means 2*17 elts:
		int dif_offsets[RADIX], dif_p20_cperms[34], dif_p20_lo_offset[32], dif_phi[ODD_RADIX];
		int dit_offsets[RADIX], dit_p20_cperms[34], dit_p20_lo_offset[32], dit_phi[ODD_RADIX];

		int incr,j,j1,j2,jt,jp,k,l,ntmp;
		// incr = Carry-chain wts-multipliers recurrence length, which must divide
		// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 18|36|72 for avx512,avx,sse, respectively:
	  #ifdef USE_AVX512	// Oddly, FMA-based builds appear to need slightly shorter chains for radices 144,288
		const int incr_long =  9,incr_med = 6,incr_short = 3;
	  #elif defined(USE_AVX2)
		const int incr_long = 12,incr_med = 6,incr_short = 4;
	  #else
		const int incr_long = 18,incr_med = 9,incr_short = 4;
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
	#ifdef USE_AVX2
		// Due to GCC macro argc limit of 30, to enable 16-register data-doubled version of the radix-9 macros need 2 length-9 ptr arrays:
		vec_dbl *rad9_iptr[9], *rad9_optr[9];
	#endif
		double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode

	#ifdef USE_SSE2

		const double crnd = 3.0*0x4000000*0x2000000;
		int *itmp,*itm2;	// Pointer into the bjmodn array
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		double *add0, *add1, *add2, *add3;
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *tmp,*tm0,*tm1,*tm2,	// Non-static utility ptrs
			*va0,*va1,*va2,*va3,*va4,*va5,*va6,*va7,*va8,
			*vb0,*vb1,*vb2,*vb3,*vb4,*vb5,*vb6,*vb7,*vb8,
				*two,*one,*sqrt2,*isrt2,*xcc1,*xss1,*xcc2,*xss2,*xcc3,*xss3,				// radix-32 DFT trig consts
								*ycc1,*yss1,*ycc2,*yss2,*ycc3m1,*yss3,*ycc4,*yss4,	// radiy-9 DFT trig consts
			*max_err, *sse2_rnd, *half_arr,
			*r00,	// Head of RADIX*vec_cmplx-sized local store #1
			*s1p00,	// Head of RADIX*vec_cmplx-sized local store #2
			*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx
		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else

	const double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
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
		p1 = NDIVR;		p10 = NDIVR<<4;		p100 = NDIVR<<8;
		p2 = p1 + p1;	p20 = p10 + p10;	p110 = p100 + p10;
		p3 = p2 + p1;	p30 = p20 + p10;
		p4 = p3 + p1;	p40 = p30 + p10;
		p5 = p4 + p1;	p50 = p40 + p10;
		p6 = p5 + p1;	p60 = p50 + p10;
		p7 = p6 + p1;	p70 = p60 + p10;
		p8 = p7 + p1;	p80 = p70 + p10;
		p9 = p8 + p1;	p90 = p80 + p10;
		pa = p9 + p1;	pa0 = p90 + p10;
		pb = pa + p1;	pb0 = pa0 + p10;
		pc = pb + p1;	pc0 = pb0 + p10;
		pd = pc + p1;	pd0 = pc0 + p10;
		pe = pd + p1;	pe0 = pd0 + p10;
		pf = pe + p1;	pf0 = pe0 + p10;

		p1 += ( (p1 >> DAT_BITS) << PAD_BITS );		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p2 += ( (p2 >> DAT_BITS) << PAD_BITS );		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p3 += ( (p3 >> DAT_BITS) << PAD_BITS );		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
		p4 += ( (p4 >> DAT_BITS) << PAD_BITS );		p40 += ( (p40 >> DAT_BITS) << PAD_BITS );
		p5 += ( (p5 >> DAT_BITS) << PAD_BITS );		p50 += ( (p50 >> DAT_BITS) << PAD_BITS );
		p6 += ( (p6 >> DAT_BITS) << PAD_BITS );		p60 += ( (p60 >> DAT_BITS) << PAD_BITS );
		p7 += ( (p7 >> DAT_BITS) << PAD_BITS );		p70 += ( (p70 >> DAT_BITS) << PAD_BITS );
		p8 += ( (p8 >> DAT_BITS) << PAD_BITS );		p80 += ( (p80 >> DAT_BITS) << PAD_BITS );
		p9 += ( (p9 >> DAT_BITS) << PAD_BITS );		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pa += ( (pa >> DAT_BITS) << PAD_BITS );		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pb += ( (pb >> DAT_BITS) << PAD_BITS );		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		pc += ( (pc >> DAT_BITS) << PAD_BITS );		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pd += ( (pd >> DAT_BITS) << PAD_BITS );		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		pe += ( (pe >> DAT_BITS) << PAD_BITS );		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		pf += ( (pf >> DAT_BITS) << PAD_BITS );		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
													p100 += ( (p100 >> DAT_BITS) << PAD_BITS );
													p110 += ( (p110 >> DAT_BITS) << PAD_BITS );
		l = 0;
		dif_phi[l++] =   0;
		dif_phi[l++] = p100;
		dif_phi[l++] = pa0;
		dif_phi[l++] = p40;
		dif_phi[l++] = pe0;
		dif_phi[l++] = p80;
		dif_phi[l++] = p20;
		dif_phi[l++] = pc0;
		dif_phi[l++] = p60;
		l = 0;
		dit_phi[l++] =   0;
		dit_phi[l++] = p20;
		dit_phi[l++] = p40;
		dit_phi[l++] = p80;
		dit_phi[l++] = pa0;
		dit_phi[l++] = p60;
		dit_phi[l++] =p100;
		dit_phi[l++] = pc0;
		dit_phi[l++] = pe0;
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
		poff[0x40+0] =p100; poff[0x40+1] =p100+p4; poff[0x40+2] =p100+p8; poff[0x40+3] =p100+pc;
		poff[0x44+0] =p110; poff[0x44+1] =p110+p4; poff[0x44+2] =p110+p8; poff[0x44+3] =p110+pc;

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

	   #ifdef USE_SSE2
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p10*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:

		// Init storage for 2 circular-shifts perms of a basic 9-vector, with shift count in [0,8] that means 2*17
		l = 0;
		// Even multiples of p10 cshift array:
		dif_p20_cperms[l++] = 0;
		dif_p20_cperms[l++] = 0x100<<1;
		dif_p20_cperms[l++] = 0xe0<<1;
		dif_p20_cperms[l++] = 0xc0<<1;
		dif_p20_cperms[l++] = 0xa0<<1;
		dif_p20_cperms[l++] = 0x80<<1;
		dif_p20_cperms[l++] = 0x60<<1;
		dif_p20_cperms[l++] = 0x40<<1;
		dif_p20_cperms[l++] = 0x20<<1;
		while(l < 2*ODD_RADIX-1) {
			dif_p20_cperms[l] = dif_p20_cperms[l - ODD_RADIX]; ++l;
		}
		// Odd multiples of p10 cshift array:
		dif_p20_cperms[l++] = 0x110<<1;
		dif_p20_cperms[l++] = 0xf0<<1;
		dif_p20_cperms[l++] = 0xd0<<1;
		dif_p20_cperms[l++] = 0xb0<<1;
		dif_p20_cperms[l++] = 0x90<<1;
		dif_p20_cperms[l++] = 0x70<<1;
		dif_p20_cperms[l++] = 0x50<<1;
		dif_p20_cperms[l++] = 0x30<<1;
		dif_p20_cperms[l++] = 0x10<<1;
		while(l < 4*ODD_RADIX-2) {
			dif_p20_cperms[l] = dif_p20_cperms[l - ODD_RADIX]; ++l;
		}
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-9 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-8; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-8 index, bits <4:30> for the p0-f. In SIMD mode again replace p[0-f] with 0x[0-f]<<1 for use with
		// contig-local-mem, thus replace 'p' prefixes with 0x and 4-left-shifts with << 5 to account for the << 1:
		l = 0;
		dif_p20_lo_offset[l++] = ((  0 << 5) + 0);
		dif_p20_lo_offset[l++] = ((0x7 << 5) + 0)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xe << 5) + 1);
		dif_p20_lo_offset[l++] = ((0x5 << 5) + 1);
		dif_p20_lo_offset[l++] = ((0xc << 5) + 1)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x3 << 5) + 1)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xa << 5) + 2);
		dif_p20_lo_offset[l++] = ((0x1 << 5) + 2);
		dif_p20_lo_offset[l++] = ((0x8 << 5) + 2)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xf << 5) + 3);
		dif_p20_lo_offset[l++] = ((0x6 << 5) + 3);
		dif_p20_lo_offset[l++] = ((0xd << 5) + 3)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x4 << 5) + 3)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xb << 5) + 4);
		dif_p20_lo_offset[l++] = ((0x2 << 5) + 4);
		dif_p20_lo_offset[l++] = ((0x9 << 5) + 4)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((  0 << 5) + 4)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x7 << 5) + 5);
		dif_p20_lo_offset[l++] = ((0xe << 5) + 5)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x5 << 5) + 5)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xc << 5) + 6);
		dif_p20_lo_offset[l++] = ((0x3 << 5) + 6);
		dif_p20_lo_offset[l++] = ((0xa << 5) + 6)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x1 << 5) + 6)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x8 << 5) + 7);
		dif_p20_lo_offset[l++] = ((0xf << 5) + 7)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x6 << 5) + 7)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xd << 5) + 8);
		dif_p20_lo_offset[l++] = ((0x4 << 5) + 8);
		dif_p20_lo_offset[l++] = ((0xb << 5) + 8)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x2 << 5) + 8)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x9 << 5) + 0);

	   #else

		// Init storage for 2 circular-shifts perms of a basic 9-vector, with shift count in [0,8] that means 2*17
		l = 0;
		// Even multiples of p10 cshift array:
		dif_p20_cperms[l++] = 0;
		dif_p20_cperms[l++] = p100;
		dif_p20_cperms[l++] = pe0;
		dif_p20_cperms[l++] = pc0;
		dif_p20_cperms[l++] = pa0;
		dif_p20_cperms[l++] = p80;
		dif_p20_cperms[l++] = p60;
		dif_p20_cperms[l++] = p40;
		dif_p20_cperms[l++] = p20;
		while(l < 2*ODD_RADIX-1) {
			dif_p20_cperms[l] = dif_p20_cperms[l - ODD_RADIX]; ++l;
		}
		// Odd multiples of p10 cshift array:
		dif_p20_cperms[l++] = p110;
		dif_p20_cperms[l++] = pf0;
		dif_p20_cperms[l++] = pd0;
		dif_p20_cperms[l++] = pb0;
		dif_p20_cperms[l++] = p90;
		dif_p20_cperms[l++] = p70;
		dif_p20_cperms[l++] = p50;
		dif_p20_cperms[l++] = p30;
		dif_p20_cperms[l++] = p10;
		while(l < 4*ODD_RADIX-2) {
			dif_p20_cperms[l] = dif_p20_cperms[l - ODD_RADIX]; ++l;
		}
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-9 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-8; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-8 index, bits <4:30> for the p0-f:
		l = 0;
		dif_p20_lo_offset[l++] = (( 0 << 4) + 0);
		dif_p20_lo_offset[l++] = ((p7 << 4) + 0)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pe << 4) + 1);
		dif_p20_lo_offset[l++] = ((p5 << 4) + 1);
		dif_p20_lo_offset[l++] = ((pc << 4) + 1)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p3 << 4) + 1)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pa << 4) + 2);
		dif_p20_lo_offset[l++] = ((p1 << 4) + 2);
		dif_p20_lo_offset[l++] = ((p8 << 4) + 2)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pf << 4) + 3);
		dif_p20_lo_offset[l++] = ((p6 << 4) + 3);
		dif_p20_lo_offset[l++] = ((pd << 4) + 3)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p4 << 4) + 3)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pb << 4) + 4);
		dif_p20_lo_offset[l++] = ((p2 << 4) + 4);
		dif_p20_lo_offset[l++] = ((p9 << 4) + 4)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = (( 0 << 4) + 4)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p7 << 4) + 5);
		dif_p20_lo_offset[l++] = ((pe << 4) + 5)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p5 << 4) + 5)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pc << 4) + 6);
		dif_p20_lo_offset[l++] = ((p3 << 4) + 6);
		dif_p20_lo_offset[l++] = ((pa << 4) + 6)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p1 << 4) + 6)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p8 << 4) + 7);
		dif_p20_lo_offset[l++] = ((pf << 4) + 7)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p6 << 4) + 7)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pd << 4) + 8);
		dif_p20_lo_offset[l++] = ((p4 << 4) + 8);
		dif_p20_lo_offset[l++] = ((pb << 4) + 8)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p2 << 4) + 8)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p9 << 4) + 0);

	   #endif	// sse2?

	// dif_offsets are w.r.to a-array, need 9 distinct sets of these, one for each DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		// Set 0: [0,1,3,2,7,6,5,4,e,f,c,d,a,b,8,9 + p00],[c,d,f,e,8,9,b,a,5,4,6,7,1,0,2,3 + p10]
		l = 0;
		dif_offsets[0x00] =  0;		dif_offsets[0x10] = pc+p10;
		dif_offsets[0x01] = p1;		dif_offsets[0x11] = pd+p10;
		dif_offsets[0x02] = p3;		dif_offsets[0x12] = pf+p10;
		dif_offsets[0x03] = p2;		dif_offsets[0x13] = pe+p10;
		dif_offsets[0x04] = p7;		dif_offsets[0x14] = p8+p10;
		dif_offsets[0x05] = p6;		dif_offsets[0x15] = p9+p10;
		dif_offsets[0x06] = p5;		dif_offsets[0x16] = pb+p10;
		dif_offsets[0x07] = p4;		dif_offsets[0x17] = pa+p10;
		dif_offsets[0x08] = pe;		dif_offsets[0x18] = p5+p10;
		dif_offsets[0x09] = pf;		dif_offsets[0x19] = p4+p10;
		dif_offsets[0x0a] = pc;		dif_offsets[0x1a] = p6+p10;
		dif_offsets[0x0b] = pd;		dif_offsets[0x1b] = p7+p10;
		dif_offsets[0x0c] = pa;		dif_offsets[0x1c] = p1+p10;
		dif_offsets[0x0d] = pb;		dif_offsets[0x1d] =    p10;
		dif_offsets[0x0e] = p8;		dif_offsets[0x1e] = p2+p10;
		dif_offsets[0x0f] = p9;		dif_offsets[0x1f] = p3+p10;
		// Set 1: [3,2,1,0,5,4,6,7,c,d,f,e,8,9,b,a +p100],[f,e,d,c,b,a,9,8,6,7,4,5,2,3,0,1 +p110]
		l += 32;
		dif_offsets[l+0x00] = p3;		dif_offsets[l+0x10] = pf+p10;
		dif_offsets[l+0x01] = p2;		dif_offsets[l+0x11] = pe+p10;
		dif_offsets[l+0x02] = p1;		dif_offsets[l+0x12] = pd+p10;
		dif_offsets[l+0x03] =  0;		dif_offsets[l+0x13] = pc+p10;
		dif_offsets[l+0x04] = p5;		dif_offsets[l+0x14] = pb+p10;
		dif_offsets[l+0x05] = p4;		dif_offsets[l+0x15] = pa+p10;
		dif_offsets[l+0x06] = p6;		dif_offsets[l+0x16] = p9+p10;
		dif_offsets[l+0x07] = p7;		dif_offsets[l+0x17] = p8+p10;
		dif_offsets[l+0x08] = pc;		dif_offsets[l+0x18] = p6+p10;
		dif_offsets[l+0x09] = pd;		dif_offsets[l+0x19] = p7+p10;
		dif_offsets[l+0x0a] = pf;		dif_offsets[l+0x1a] = p4+p10;
		dif_offsets[l+0x0b] = pe;		dif_offsets[l+0x1b] = p5+p10;
		dif_offsets[l+0x0c] = p8;		dif_offsets[l+0x1c] = p2+p10;
		dif_offsets[l+0x0d] = p9;		dif_offsets[l+0x1d] = p3+p10;
		dif_offsets[l+0x0e] = pb;		dif_offsets[l+0x1e] =    p10;
		dif_offsets[l+0x0f] = pa;		dif_offsets[l+0x1f] = p1+p10;
		// Set 2: [1,0,2,3,6,7,4,5,f,e,d,c,b,a,9,8 + pb0],[3,2,1,0,5,4,6,7,c,d,f,e,8,9,b,a + pa0]
		l += 32;
		dif_offsets[l+0x00] = p1+p10;		dif_offsets[l+0x10] = p3;
		dif_offsets[l+0x01] =    p10;		dif_offsets[l+0x11] = p2;
		dif_offsets[l+0x02] = p2+p10;		dif_offsets[l+0x12] = p1;
		dif_offsets[l+0x03] = p3+p10;		dif_offsets[l+0x13] =  0;
		dif_offsets[l+0x04] = p6+p10;		dif_offsets[l+0x14] = p5;
		dif_offsets[l+0x05] = p7+p10;		dif_offsets[l+0x15] = p4;
		dif_offsets[l+0x06] = p4+p10;		dif_offsets[l+0x16] = p6;
		dif_offsets[l+0x07] = p5+p10;		dif_offsets[l+0x17] = p7;
		dif_offsets[l+0x08] = pf+p10;		dif_offsets[l+0x18] = pc;
		dif_offsets[l+0x09] = pe+p10;		dif_offsets[l+0x19] = pd;
		dif_offsets[l+0x0a] = pd+p10;		dif_offsets[l+0x1a] = pf;
		dif_offsets[l+0x0b] = pc+p10;		dif_offsets[l+0x1b] = pe;
		dif_offsets[l+0x0c] = pb+p10;		dif_offsets[l+0x1c] = p8;
		dif_offsets[l+0x0d] = pa+p10;		dif_offsets[l+0x1d] = p9;
		dif_offsets[l+0x0e] = p9+p10;		dif_offsets[l+0x1e] = pb;
		dif_offsets[l+0x0f] = p8+p10;		dif_offsets[l+0x1f] = pa;
		// Set 3: [a,b,8,9,c,d,f,e,3,2,1,0,5,4,6,7 + p40],[1,0,2,3,6,7,4,5,f,e,d,c,b,a,9,8 + p50]
		l += 32;
		dif_offsets[l+0x00] = pa;		dif_offsets[l+0x10] = p1+p10;
		dif_offsets[l+0x01] = pb;		dif_offsets[l+0x11] =    p10;
		dif_offsets[l+0x02] = p8;		dif_offsets[l+0x12] = p2+p10;
		dif_offsets[l+0x03] = p9;		dif_offsets[l+0x13] = p3+p10;
		dif_offsets[l+0x04] = pc;		dif_offsets[l+0x14] = p6+p10;
		dif_offsets[l+0x05] = pd;		dif_offsets[l+0x15] = p7+p10;
		dif_offsets[l+0x06] = pf;		dif_offsets[l+0x16] = p4+p10;
		dif_offsets[l+0x07] = pe;		dif_offsets[l+0x17] = p5+p10;
		dif_offsets[l+0x08] = p3;		dif_offsets[l+0x18] = pf+p10;
		dif_offsets[l+0x09] = p2;		dif_offsets[l+0x19] = pe+p10;
		dif_offsets[l+0x0a] = p1;		dif_offsets[l+0x1a] = pd+p10;
		dif_offsets[l+0x0b] =  0;		dif_offsets[l+0x1b] = pc+p10;
		dif_offsets[l+0x0c] = p5;		dif_offsets[l+0x1c] = pb+p10;
		dif_offsets[l+0x0d] = p4;		dif_offsets[l+0x1d] = pa+p10;
		dif_offsets[l+0x0e] = p6;		dif_offsets[l+0x1e] = p9+p10;
		dif_offsets[l+0x0f] = p7;		dif_offsets[l+0x1f] = p8+p10;
		// Set 4: [8,9,b,a,f,e,d,c,1,0,2,3,6,7,4,5 + pf0],[a,b,8,9,c,d,f,e,3,2,1,0,5,4,6,7 + pe0]
		l += 32;
		dif_offsets[l+0x00] = p8+p10;		dif_offsets[l+0x10] = pa;
		dif_offsets[l+0x01] = p9+p10;		dif_offsets[l+0x11] = pb;
		dif_offsets[l+0x02] = pb+p10;		dif_offsets[l+0x12] = p8;
		dif_offsets[l+0x03] = pa+p10;		dif_offsets[l+0x13] = p9;
		dif_offsets[l+0x04] = pf+p10;		dif_offsets[l+0x14] = pc;
		dif_offsets[l+0x05] = pe+p10;		dif_offsets[l+0x15] = pd;
		dif_offsets[l+0x06] = pd+p10;		dif_offsets[l+0x16] = pf;
		dif_offsets[l+0x07] = pc+p10;		dif_offsets[l+0x17] = pe;
		dif_offsets[l+0x08] = p1+p10;		dif_offsets[l+0x18] = p3;
		dif_offsets[l+0x09] =    p10;		dif_offsets[l+0x19] = p2;
		dif_offsets[l+0x0a] = p2+p10;		dif_offsets[l+0x1a] = p1;
		dif_offsets[l+0x0b] = p3+p10;		dif_offsets[l+0x1b] =  0;
		dif_offsets[l+0x0c] = p6+p10;		dif_offsets[l+0x1c] = p5;
		dif_offsets[l+0x0d] = p7+p10;		dif_offsets[l+0x1d] = p4;
		dif_offsets[l+0x0e] = p4+p10;		dif_offsets[l+0x1e] = p6;
		dif_offsets[l+0x0f] = p5+p10;		dif_offsets[l+0x1f] = p7;
		// Set 5: [7,6,5,4,3,2,1,0,a,b,8,9,c,d,f,e + p80],[8,9,b,a,f,e,d,c,1,0,2,3,6,7,4,5 + p90]
		l += 32;
		dif_offsets[l+0x00] = p7;		dif_offsets[l+0x10] = p8+p10;
		dif_offsets[l+0x01] = p6;		dif_offsets[l+0x11] = p9+p10;
		dif_offsets[l+0x02] = p5;		dif_offsets[l+0x12] = pb+p10;
		dif_offsets[l+0x03] = p4;		dif_offsets[l+0x13] = pa+p10;
		dif_offsets[l+0x04] = p3;		dif_offsets[l+0x14] = pf+p10;
		dif_offsets[l+0x05] = p2;		dif_offsets[l+0x15] = pe+p10;
		dif_offsets[l+0x06] = p1;		dif_offsets[l+0x16] = pd+p10;
		dif_offsets[l+0x07] =  0;		dif_offsets[l+0x17] = pc+p10;
		dif_offsets[l+0x08] = pa;		dif_offsets[l+0x18] = p1+p10;
		dif_offsets[l+0x09] = pb;		dif_offsets[l+0x19] =    p10;
		dif_offsets[l+0x0a] = p8;		dif_offsets[l+0x1a] = p2+p10;
		dif_offsets[l+0x0b] = p9;		dif_offsets[l+0x1b] = p3+p10;
		dif_offsets[l+0x0c] = pc;		dif_offsets[l+0x1c] = p6+p10;
		dif_offsets[l+0x0d] = pd;		dif_offsets[l+0x1d] = p7+p10;
		dif_offsets[l+0x0e] = pf;		dif_offsets[l+0x1e] = p4+p10;
		dif_offsets[l+0x0f] = pe;		dif_offsets[l+0x1f] = p5+p10;
		// Set 6: [5,4,6,7,1,0,2,3,8,9,b,a,f,e,d,c + p30],[7,6,5,4,3,2,1,0,a,b,8,9,c,d,f,e + p20]
		l += 32;
		dif_offsets[l+0x00] = p5+p10;		dif_offsets[l+0x10] = p7;
		dif_offsets[l+0x01] = p4+p10;		dif_offsets[l+0x11] = p6;
		dif_offsets[l+0x02] = p6+p10;		dif_offsets[l+0x12] = p5;
		dif_offsets[l+0x03] = p7+p10;		dif_offsets[l+0x13] = p4;
		dif_offsets[l+0x04] = p1+p10;		dif_offsets[l+0x14] = p3;
		dif_offsets[l+0x05] =    p10;		dif_offsets[l+0x15] = p2;
		dif_offsets[l+0x06] = p2+p10;		dif_offsets[l+0x16] = p1;
		dif_offsets[l+0x07] = p3+p10;		dif_offsets[l+0x17] =  0;
		dif_offsets[l+0x08] = p8+p10;		dif_offsets[l+0x18] = pa;
		dif_offsets[l+0x09] = p9+p10;		dif_offsets[l+0x19] = pb;
		dif_offsets[l+0x0a] = pb+p10;		dif_offsets[l+0x1a] = p8;
		dif_offsets[l+0x0b] = pa+p10;		dif_offsets[l+0x1b] = p9;
		dif_offsets[l+0x0c] = pf+p10;		dif_offsets[l+0x1c] = pc;
		dif_offsets[l+0x0d] = pe+p10;		dif_offsets[l+0x1d] = pd;
		dif_offsets[l+0x0e] = pd+p10;		dif_offsets[l+0x1e] = pf;
		dif_offsets[l+0x0f] = pc+p10;		dif_offsets[l+0x1f] = pe;
		// Set 7: [e,f,c,d,a,b,8,9,7,6,5,4,3,2,1,0 + pc0],[5,4,6,7,1,0,2,3,8,9,b,a,f,e,d,c + pd0]
		l += 32;
		dif_offsets[l+0x00] = pe;		dif_offsets[l+0x10] = p5+p10;
		dif_offsets[l+0x01] = pf;		dif_offsets[l+0x11] = p4+p10;
		dif_offsets[l+0x02] = pc;		dif_offsets[l+0x12] = p6+p10;
		dif_offsets[l+0x03] = pd;		dif_offsets[l+0x13] = p7+p10;
		dif_offsets[l+0x04] = pa;		dif_offsets[l+0x14] = p1+p10;
		dif_offsets[l+0x05] = pb;		dif_offsets[l+0x15] =    p10;
		dif_offsets[l+0x06] = p8;		dif_offsets[l+0x16] = p2+p10;
		dif_offsets[l+0x07] = p9;		dif_offsets[l+0x17] = p3+p10;
		dif_offsets[l+0x08] = p7;		dif_offsets[l+0x18] = p8+p10;
		dif_offsets[l+0x09] = p6;		dif_offsets[l+0x19] = p9+p10;
		dif_offsets[l+0x0a] = p5;		dif_offsets[l+0x1a] = pb+p10;
		dif_offsets[l+0x0b] = p4;		dif_offsets[l+0x1b] = pa+p10;
		dif_offsets[l+0x0c] = p3;		dif_offsets[l+0x1c] = pf+p10;
		dif_offsets[l+0x0d] = p2;		dif_offsets[l+0x1d] = pe+p10;
		dif_offsets[l+0x0e] = p1;		dif_offsets[l+0x1e] = pd+p10;
		dif_offsets[l+0x0f] =  0;		dif_offsets[l+0x1f] = pc+p10;
		// Set 8: [c,d,f,e,8,9,b,a,5,4,6,7,1,0,2,3 + p70],[e,f,c,d,a,b,8,9,7,6,5,4,3,2,1,0 + p60]
		l += 32;
		dif_offsets[l+0x00] = pc+p10;		dif_offsets[l+0x10] = pe;
		dif_offsets[l+0x01] = pd+p10;		dif_offsets[l+0x11] = pf;
		dif_offsets[l+0x02] = pf+p10;		dif_offsets[l+0x12] = pc;
		dif_offsets[l+0x03] = pe+p10;		dif_offsets[l+0x13] = pd;
		dif_offsets[l+0x04] = p8+p10;		dif_offsets[l+0x14] = pa;
		dif_offsets[l+0x05] = p9+p10;		dif_offsets[l+0x15] = pb;
		dif_offsets[l+0x06] = pb+p10;		dif_offsets[l+0x16] = p8;
		dif_offsets[l+0x07] = pa+p10;		dif_offsets[l+0x17] = p9;
		dif_offsets[l+0x08] = p5+p10;		dif_offsets[l+0x18] = p7;
		dif_offsets[l+0x09] = p4+p10;		dif_offsets[l+0x19] = p6;
		dif_offsets[l+0x0a] = p6+p10;		dif_offsets[l+0x1a] = p5;
		dif_offsets[l+0x0b] = p7+p10;		dif_offsets[l+0x1b] = p4;
		dif_offsets[l+0x0c] = p1+p10;		dif_offsets[l+0x1c] = p3;
		dif_offsets[l+0x0d] =    p10;		dif_offsets[l+0x1d] = p2;
		dif_offsets[l+0x0e] = p2+p10;		dif_offsets[l+0x1e] = p1;
		dif_offsets[l+0x0f] = p3+p10;		dif_offsets[l+0x1f] =  0;
	  #ifdef USE_SSE2
		// IN SIMD mode preshift all the above offsets << 3 to turn into double=array pointer offsets:
		for(l = 0; l < RADIX; l++) {
			dif_offsets[l] <<= 3;
		}
	  #endif

	/*** DIT indexing stuff: ***/

	   #ifdef USE_SSE2
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p10*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:
		l = 0;
		// Init storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,8] that means 2*17
		// Even multiples of p10 cshift array: evn0 := 0x10*[00,02,04,06,08,0a,0c,0e,10]
		dit_p20_cperms[l++] =    0;
		dit_p20_cperms[l++] = 0x20<<1;
		dit_p20_cperms[l++] = 0x40<<1;
		dit_p20_cperms[l++] = 0x60<<1;
		dit_p20_cperms[l++] = 0x80<<1;
		dit_p20_cperms[l++] = 0xa0<<1;
		dit_p20_cperms[l++] = 0xc0<<1;
		dit_p20_cperms[l++] = 0xe0<<1;
		dit_p20_cperms[l++] = 0x100<<1;
		while(l < 2*ODD_RADIX-1) {
			dit_p20_cperms[l] = dit_p20_cperms[l - ODD_RADIX]; ++l;
		}
		// Odd multiples of p10 cshift array: odd0 := 0x10*[03,05,07,09,0b,0d,0f,11,01]
		dit_p20_cperms[l++] = 0x30<<1;
		dit_p20_cperms[l++] = 0x50<<1;
		dit_p20_cperms[l++] = 0x70<<1;
		dit_p20_cperms[l++] = 0x90<<1;
		dit_p20_cperms[l++] = 0xb0<<1;
		dit_p20_cperms[l++] = 0xd0<<1;
		dit_p20_cperms[l++] = 0xf0<<1;
		dit_p20_cperms[l++] = 0x110<<1;
		dit_p20_cperms[l++] = 0x10<<1;
		while(l < 4*ODD_RADIX-2) {
			dit_p20_cperms[l] = dit_p20_cperms[l - ODD_RADIX]; ++l;
		}

		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-9 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-8 index, bits <4:30> for the p0-f. In SIMD mode again replace p[0-f] with 0x[0-f]<<1 for use with
		// contig-local-mem, thus replace 'p' prefixes with 0x and 4-left-shifts with << 5 to account for the << 1:
		l = 0;
		dit_p20_lo_offset[l++] = ((  0 << 5) + 0);
		dit_p20_lo_offset[l++] = ((0xf << 5) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xe << 5) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xd << 5) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xc << 5) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xb << 5) + 8) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xa << 5) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x9 << 5) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x8 << 5) + 5) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x7 << 5) + 7) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x6 << 5) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x5 << 5) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x4 << 5) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x3 << 5) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x2 << 5) + 8) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x1 << 5) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((  0 << 5) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xf << 5) + 6);
		dit_p20_lo_offset[l++] = ((0xe << 5) + 8);
		dit_p20_lo_offset[l++] = ((0xd << 5) + 1);
		dit_p20_lo_offset[l++] = ((0xc << 5) + 3);
		dit_p20_lo_offset[l++] = ((0xb << 5) + 5);
		dit_p20_lo_offset[l++] = ((0xa << 5) + 7);
		dit_p20_lo_offset[l++] = ((0x9 << 5) + 0);
		dit_p20_lo_offset[l++] = ((0x8 << 5) + 2);
		dit_p20_lo_offset[l++] = ((0x7 << 5) + 4);
		dit_p20_lo_offset[l++] = ((0x6 << 5) + 6);
		dit_p20_lo_offset[l++] = ((0x5 << 5) + 8);
		dit_p20_lo_offset[l++] = ((0x4 << 5) + 1);
		dit_p20_lo_offset[l++] = ((0x3 << 5) + 3);
		dit_p20_lo_offset[l++] = ((0x2 << 5) + 5);
		dit_p20_lo_offset[l++] = ((0x1 << 5) + 7);

	   #else

		l = 0;
		// Init storage for 2 circular-shifts perms of a basic 5-vector, with shift count in [0,8] that means 2*17
		// Even multiples of p10 cshift array: evn0 := p10*[00,02,04,06,08,0a,0c,0e,10]
		dit_p20_cperms[l++] =   0;
		dit_p20_cperms[l++] = p20;
		dit_p20_cperms[l++] = p40;
		dit_p20_cperms[l++] = p60;
		dit_p20_cperms[l++] = p80;
		dit_p20_cperms[l++] = pa0;
		dit_p20_cperms[l++] = pc0;
		dit_p20_cperms[l++] = pe0;
		dit_p20_cperms[l++] = p100;
		while(l < 2*ODD_RADIX-1) {
			dit_p20_cperms[l] = dit_p20_cperms[l - ODD_RADIX]; ++l;
		}
		// Odd multiples of p10 cshift array: odd0 := p10*[03,05,07,09,0b,0d,0f,11,01]
		dit_p20_cperms[l++] = p30;
		dit_p20_cperms[l++] = p50;
		dit_p20_cperms[l++] = p70;
		dit_p20_cperms[l++] = p90;
		dit_p20_cperms[l++] = pb0;
		dit_p20_cperms[l++] = pd0;
		dit_p20_cperms[l++] = pf0;
		dit_p20_cperms[l++] = p110;
		dit_p20_cperms[l++] = p10;
		while(l < 4*ODD_RADIX-2) {
			dit_p20_cperms[l] = dit_p20_cperms[l - ODD_RADIX]; ++l;
		}

		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-9 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-4 index, bits <4:30> for the p0-f:
		l = 0;
		dit_p20_lo_offset[l++] = (( 0 << 4) + 0);
		dit_p20_lo_offset[l++] = ((pf << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pe << 4) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pd << 4) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pc << 4) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pb << 4) + 8) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pa << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p9 << 4) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p8 << 4) + 5) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p7 << 4) + 7) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p6 << 4) + 0) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p5 << 4) + 2) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p4 << 4) + 4) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p3 << 4) + 6) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p2 << 4) + 8) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p1 << 4) + 1) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = (( 0 << 4) + 3) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pf << 4) + 6);
		dit_p20_lo_offset[l++] = ((pe << 4) + 8);
		dit_p20_lo_offset[l++] = ((pd << 4) + 1);
		dit_p20_lo_offset[l++] = ((pc << 4) + 3);
		dit_p20_lo_offset[l++] = ((pb << 4) + 5);
		dit_p20_lo_offset[l++] = ((pa << 4) + 7);
		dit_p20_lo_offset[l++] = ((p9 << 4) + 0);
		dit_p20_lo_offset[l++] = ((p8 << 4) + 2);
		dit_p20_lo_offset[l++] = ((p7 << 4) + 4);
		dit_p20_lo_offset[l++] = ((p6 << 4) + 6);
		dit_p20_lo_offset[l++] = ((p5 << 4) + 8);
		dit_p20_lo_offset[l++] = ((p4 << 4) + 1);
		dit_p20_lo_offset[l++] = ((p3 << 4) + 3);
		dit_p20_lo_offset[l++] = ((p2 << 4) + 5);
		dit_p20_lo_offset[l++] = ((p1 << 4) + 7);

	   #endif	// sse2?

	// dit_offsets are w.r.to a-array, need 9 distinct sets of these, one for each radix-32 DFT.
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
		// Set 1: [5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p30],[5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p20]
		l += 32;
		dit_offsets[l+0x00] = p5+p10;		dit_offsets[l+0x10] = p5;
		dit_offsets[l+0x01] = p4+p10;		dit_offsets[l+0x11] = p4;
		dit_offsets[l+0x02] = p6+p10;		dit_offsets[l+0x12] = p6;
		dit_offsets[l+0x03] = p7+p10;		dit_offsets[l+0x13] = p7;
		dit_offsets[l+0x04] = p1+p10;		dit_offsets[l+0x14] = p1;
		dit_offsets[l+0x05] =    p10;		dit_offsets[l+0x15] =  0;
		dit_offsets[l+0x06] = p2+p10;		dit_offsets[l+0x16] = p2;
		dit_offsets[l+0x07] = p3+p10;		dit_offsets[l+0x17] = p3;
		dit_offsets[l+0x08] = p9+p10;		dit_offsets[l+0x18] = p9;
		dit_offsets[l+0x09] = p8+p10;		dit_offsets[l+0x19] = p8;
		dit_offsets[l+0x0a] = pa+p10;		dit_offsets[l+0x1a] = pa;
		dit_offsets[l+0x0b] = pb+p10;		dit_offsets[l+0x1b] = pb;
		dit_offsets[l+0x0c] = pe+p10;		dit_offsets[l+0x1c] = pe;
		dit_offsets[l+0x0d] = pf+p10;		dit_offsets[l+0x1d] = pf;
		dit_offsets[l+0x0e] = pc+p10;		dit_offsets[l+0x1e] = pc;
		dit_offsets[l+0x0f] = pd+p10;		dit_offsets[l+0x1f] = pd;
		// Set 2: [a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p40],[2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p50]
		l += 32;
		dit_offsets[l+0x00] = pa;		dit_offsets[l+0x10] = p2+p10;
		dit_offsets[l+0x01] = pb;		dit_offsets[l+0x11] = p3+p10;
		dit_offsets[l+0x02] = p8;		dit_offsets[l+0x12] =    p10;
		dit_offsets[l+0x03] = p9;		dit_offsets[l+0x13] = p1+p10;
		dit_offsets[l+0x04] = pc;		dit_offsets[l+0x14] = p4+p10;
		dit_offsets[l+0x05] = pd;		dit_offsets[l+0x15] = p5+p10;
		dit_offsets[l+0x06] = pf;		dit_offsets[l+0x16] = p7+p10;
		dit_offsets[l+0x07] = pe;		dit_offsets[l+0x17] = p6+p10;
		dit_offsets[l+0x08] = p2;		dit_offsets[l+0x18] = pc+p10;
		dit_offsets[l+0x09] = p3;		dit_offsets[l+0x19] = pd+p10;
		dit_offsets[l+0x0a] =  0;		dit_offsets[l+0x1a] = pf+p10;
		dit_offsets[l+0x0b] = p1;		dit_offsets[l+0x1b] = pe+p10;
		dit_offsets[l+0x0c] = p4;		dit_offsets[l+0x1c] = p8+p10;
		dit_offsets[l+0x0d] = p5;		dit_offsets[l+0x1d] = p9+p10;
		dit_offsets[l+0x0e] = p7;		dit_offsets[l+0x1e] = pb+p10;
		dit_offsets[l+0x0f] = p6;		dit_offsets[l+0x1f] = pa+p10;
		// Set 3: [7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + p80],[b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p90]
		l += 32;
		dit_offsets[l+0x00] = p7;		dit_offsets[l+0x10] = pb+p10;
		dit_offsets[l+0x01] = p6;		dit_offsets[l+0x11] = pa+p10;
		dit_offsets[l+0x02] = p5;		dit_offsets[l+0x12] = p9+p10;
		dit_offsets[l+0x03] = p4;		dit_offsets[l+0x13] = p8+p10;
		dit_offsets[l+0x04] = p3;		dit_offsets[l+0x14] = pd+p10;
		dit_offsets[l+0x05] = p2;		dit_offsets[l+0x15] = pc+p10;
		dit_offsets[l+0x06] = p1;		dit_offsets[l+0x16] = pe+p10;
		dit_offsets[l+0x07] =  0;		dit_offsets[l+0x17] = pf+p10;
		dit_offsets[l+0x08] = pb;		dit_offsets[l+0x18] = p3+p10;
		dit_offsets[l+0x09] = pa;		dit_offsets[l+0x19] = p2+p10;
		dit_offsets[l+0x0a] = p9;		dit_offsets[l+0x1a] = p1+p10;
		dit_offsets[l+0x0b] = p8;		dit_offsets[l+0x1b] =    p10;
		dit_offsets[l+0x0c] = pd;		dit_offsets[l+0x1c] = p5+p10;
		dit_offsets[l+0x0d] = pc;		dit_offsets[l+0x1d] = p4+p10;
		dit_offsets[l+0x0e] = pe;		dit_offsets[l+0x1e] = p6+p10;
		dit_offsets[l+0x0f] = pf;		dit_offsets[l+0x1f] = p7+p10;
		// Set 4: [1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pb0],[1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pa0]
		l += 32;
		dit_offsets[l+0x00] = p1+p10;		dit_offsets[l+0x10] = p1;
		dit_offsets[l+0x01] =    p10;		dit_offsets[l+0x11] =  0;
		dit_offsets[l+0x02] = p2+p10;		dit_offsets[l+0x12] = p2;
		dit_offsets[l+0x03] = p3+p10;		dit_offsets[l+0x13] = p3;
		dit_offsets[l+0x04] = p6+p10;		dit_offsets[l+0x14] = p6;
		dit_offsets[l+0x05] = p7+p10;		dit_offsets[l+0x15] = p7;
		dit_offsets[l+0x06] = p4+p10;		dit_offsets[l+0x16] = p4;
		dit_offsets[l+0x07] = p5+p10;		dit_offsets[l+0x17] = p5;
		dit_offsets[l+0x08] = pe+p10;		dit_offsets[l+0x18] = pe;
		dit_offsets[l+0x09] = pf+p10;		dit_offsets[l+0x19] = pf;
		dit_offsets[l+0x0a] = pc+p10;		dit_offsets[l+0x1a] = pc;
		dit_offsets[l+0x0b] = pd+p10;		dit_offsets[l+0x1b] = pd;
		dit_offsets[l+0x0c] = pa+p10;		dit_offsets[l+0x1c] = pa;
		dit_offsets[l+0x0d] = pb+p10;		dit_offsets[l+0x1d] = pb;
		dit_offsets[l+0x0e] = p8+p10;		dit_offsets[l+0x1e] = p8;
		dit_offsets[l+0x0f] = p9+p10;		dit_offsets[l+0x1f] = p9;
		// Set 5: [c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p70],[c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p60]
		l += 32;
		dit_offsets[l+0x00] = pc+p10;		dit_offsets[l+0x10] = pc;
		dit_offsets[l+0x01] = pd+p10;		dit_offsets[l+0x11] = pd;
		dit_offsets[l+0x02] = pf+p10;		dit_offsets[l+0x12] = pf;
		dit_offsets[l+0x03] = pe+p10;		dit_offsets[l+0x13] = pe;
		dit_offsets[l+0x04] = p8+p10;		dit_offsets[l+0x14] = p8;
		dit_offsets[l+0x05] = p9+p10;		dit_offsets[l+0x15] = p9;
		dit_offsets[l+0x06] = pb+p10;		dit_offsets[l+0x16] = pb;
		dit_offsets[l+0x07] = pa+p10;		dit_offsets[l+0x17] = pa;
		dit_offsets[l+0x08] = p4+p10;		dit_offsets[l+0x18] = p4;
		dit_offsets[l+0x09] = p5+p10;		dit_offsets[l+0x19] = p5;
		dit_offsets[l+0x0a] = p7+p10;		dit_offsets[l+0x1a] = p7;
		dit_offsets[l+0x0b] = p6+p10;		dit_offsets[l+0x1b] = p6;
		dit_offsets[l+0x0c] =    p10;		dit_offsets[l+0x1c] =  0;
		dit_offsets[l+0x0d] = p1+p10;		dit_offsets[l+0x1d] = p1;
		dit_offsets[l+0x0e] = p3+p10;		dit_offsets[l+0x1e] = p3;
		dit_offsets[l+0x0f] = p2+p10;		dit_offsets[l+0x1f] = p2;
		// Set 6: [3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b +p100],[d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 +p110]
		l += 32;
		dit_offsets[l+0x00] = p3;		dit_offsets[l+0x10] = pd+p10;
		dit_offsets[l+0x01] = p2;		dit_offsets[l+0x11] = pc+p10;
		dit_offsets[l+0x02] = p1;		dit_offsets[l+0x12] = pe+p10;
		dit_offsets[l+0x03] =  0;		dit_offsets[l+0x13] = pf+p10;
		dit_offsets[l+0x04] = p5;		dit_offsets[l+0x14] = p9+p10;
		dit_offsets[l+0x05] = p4;		dit_offsets[l+0x15] = p8+p10;
		dit_offsets[l+0x06] = p6;		dit_offsets[l+0x16] = pa+p10;
		dit_offsets[l+0x07] = p7;		dit_offsets[l+0x17] = pb+p10;
		dit_offsets[l+0x08] = pd;		dit_offsets[l+0x18] = p5+p10;
		dit_offsets[l+0x09] = pc;		dit_offsets[l+0x19] = p4+p10;
		dit_offsets[l+0x0a] = pe;		dit_offsets[l+0x1a] = p6+p10;
		dit_offsets[l+0x0b] = pf;		dit_offsets[l+0x1b] = p7+p10;
		dit_offsets[l+0x0c] = p9;		dit_offsets[l+0x1c] = p1+p10;
		dit_offsets[l+0x0d] = p8;		dit_offsets[l+0x1d] =    p10;
		dit_offsets[l+0x0e] = pa;		dit_offsets[l+0x1e] = p2+p10;
		dit_offsets[l+0x0f] = pb;		dit_offsets[l+0x1f] = p3+p10;
		// Set 7: [e,f,c,d,a,b,8,9,6,7,4,5,2,3,0,1 + pc0],[6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pd0]
		l += 32;
		dit_offsets[l+0x00] = pe;		dit_offsets[l+0x10] = p6+p10;
		dit_offsets[l+0x01] = pf;		dit_offsets[l+0x11] = p7+p10;
		dit_offsets[l+0x02] = pc;		dit_offsets[l+0x12] = p4+p10;
		dit_offsets[l+0x03] = pd;		dit_offsets[l+0x13] = p5+p10;
		dit_offsets[l+0x04] = pa;		dit_offsets[l+0x14] = p2+p10;
		dit_offsets[l+0x05] = pb;		dit_offsets[l+0x15] = p3+p10;
		dit_offsets[l+0x06] = p8;		dit_offsets[l+0x16] =    p10;
		dit_offsets[l+0x07] = p9;		dit_offsets[l+0x17] = p1+p10;
		dit_offsets[l+0x08] = p6;		dit_offsets[l+0x18] = pa+p10;
		dit_offsets[l+0x09] = p7;		dit_offsets[l+0x19] = pb+p10;
		dit_offsets[l+0x0a] = p4;		dit_offsets[l+0x1a] = p8+p10;
		dit_offsets[l+0x0b] = p5;		dit_offsets[l+0x1b] = p9+p10;
		dit_offsets[l+0x0c] = p2;		dit_offsets[l+0x1c] = pc+p10;
		dit_offsets[l+0x0d] = p3;		dit_offsets[l+0x1d] = pd+p10;
		dit_offsets[l+0x0e] =  0;		dit_offsets[l+0x1e] = pf+p10;
		dit_offsets[l+0x0f] = p1;		dit_offsets[l+0x1f] = pe+p10;
		// Set 8: [8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + pf0],[8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + pe0]
		l += 32;
		dit_offsets[l+0x00] = p8+p10;		dit_offsets[l+0x10] = p8;
		dit_offsets[l+0x01] = p9+p10;		dit_offsets[l+0x11] = p9;
		dit_offsets[l+0x02] = pb+p10;		dit_offsets[l+0x12] = pb;
		dit_offsets[l+0x03] = pa+p10;		dit_offsets[l+0x13] = pa;
		dit_offsets[l+0x04] = pf+p10;		dit_offsets[l+0x14] = pf;
		dit_offsets[l+0x05] = pe+p10;		dit_offsets[l+0x15] = pe;
		dit_offsets[l+0x06] = pd+p10;		dit_offsets[l+0x16] = pd;
		dit_offsets[l+0x07] = pc+p10;		dit_offsets[l+0x17] = pc;
		dit_offsets[l+0x08] =    p10;		dit_offsets[l+0x18] =  0;
		dit_offsets[l+0x09] = p1+p10;		dit_offsets[l+0x19] = p1;
		dit_offsets[l+0x0a] = p3+p10;		dit_offsets[l+0x1a] = p3;
		dit_offsets[l+0x0b] = p2+p10;		dit_offsets[l+0x1b] = p2;
		dit_offsets[l+0x0c] = p7+p10;		dit_offsets[l+0x1c] = p7;
		dit_offsets[l+0x0d] = p6+p10;		dit_offsets[l+0x1d] = p6;
		dit_offsets[l+0x0e] = p5+p10;		dit_offsets[l+0x1e] = p5;
		dit_offsets[l+0x0f] = p4+p10;		dit_offsets[l+0x1f] = p4;
	  #ifdef USE_SSE2
		// IN SIMD mode preshift all the above offsets << 3 to turn into double=array pointer offsets:
		for(l = 0; l < RADIX; l++) {
			dit_offsets[l] <<= 3;
		}
	  #endif

	#ifdef USE_SSE2
		tmp	= r00 = thread_arg->r00;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x240;	s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
		tmp	+= 0x243;	// Extra 3 slots here for two,one below - added those late, too lazy to rejigger all the existing offsets following
		two     = tmp - 3;	// AVX+ versions of various DFT macros assume consts 2.0,1.0,isrt2 laid out thusly
		one     = tmp - 2;
		sqrt2	= tmp - 1;
		isrt2   = tmp + 0x0;
		xcc2	= tmp + 0x1;	// Radix-32 DFT macros assume roots stored in this [8th, 16th, 32nd_1,3] order
		xss2	= tmp + 0x2;
		xcc1	= tmp + 0x3;
		xss1	= tmp + 0x4;
		xcc3	= tmp + 0x5;
		xss3	= tmp + 0x6;
		// Roots for radix-9 DFTs:
		ycc1    = tmp + 0x7;
		yss1    = tmp + 0x8;
		ycc2    = tmp + 0x9;
		yss2    = tmp + 0xa;
		ycc3m1  = tmp + 0xb;
		yss3    = tmp + 0xc;
		ycc4    = tmp + 0xd;
		yss4    = tmp + 0xe;
		tmp += 0xf;	// sc_ptr += 0x492
	  #ifdef USE_AVX512
		cy = tmp;		tmp += 0x24;	// RADIX/8 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
		cy = tmp;		tmp += 0x48;	// RADIX/4 vec_dbl slots
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(492 + 48 + 2) = 0x4dc; This is where the value of half_arr_offset288 comes from
		half_arr= tmp + 0x02;	// This table needs 20 vec_dbl in both avx and sse2 mode
	  #else
		cy = tmp;		tmp += 0x90;	// RADIX/2 vec_dbl slots
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(492 + 90 + 2) = 0x524; This is where the value of half_arr_offset288 comes from
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

		sign_mask = (uint64*)(r00 + radix288_creals_in_local_store);
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
		#include "radix288_main_carry_loop.h"

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
