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
#include "radix4032.h"	// 4032-element dif64_oidx_lo byte-array

#define RADIX 4032	// Use #define rather than const int to ensure it's really a compile-time const in the C sense
#define ODD_RADIX 63	// ODD_RADIX = [radix >> trailz(radix)]

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

  // In hex, RADIX = 0xfc0, *4 = 0x3f00
  // For Mersenne-mod we need max(4*ODD_RADIX, (16 [SSE2] or 64 [AVX]) + 4) added slots for the half_arr lookup tables.
  // 4*ODD_RADIX = 252 here, trumps both elements of the (20,68), thus use extra 252 = 0xfc slots in SSE2 mode.
  // In AVX mode.this gets supplemented by the extra storage we use for chained computation of the negacyclic-weights.
  // Add relevant number (half_arr_offset4032 + RADIX) to get required value of radix4032_creals_in_local_store:
  #ifdef USE_AVX512	// 0x3f0 fewer carry slots than AVX:
	const int half_arr_offset4032 = 0x42f2;	// + RADIX = 0x52b2; Used for thread local-storage-integrity checking
	const int radix4032_creals_in_local_store = 0x53c0;	// AVX+LOACC: (half_arr_offset4032 + RADIX) + 0xfc, and round up to nearest multiple of 16
  #elif defined(USE_AVX)
	const int half_arr_offset4032 = 0x46e2;	// + RADIX = 0x56a2; Used for thread local-storage-integrity checking
   #if HIACC	// Need extra 4*RADIX = 0xfc0 vec_dbl slots in this mode
	const int radix4032_creals_in_local_store = 0x95b0;	// AVX+HIACC: (half_arr_offset4032 + 5*RADIX) + 0xfc, and round up to nearest multiple of 16
   #else
	const int radix4032_creals_in_local_store = 0x57b0;	// AVX+LOACC: (half_arr_offset4032 + RADIX) + 0xfc, and round up to nearest multiple of 16
   #endif
  #else
	const int half_arr_offset4032 = 0x4ec2;	// + RADIX = 0x5e82; Used for thread local-storage-integrity checking
	const int radix4032_creals_in_local_store = 0x5f80;	// (half_arr_offset4032 + RADIX) + 0xfc, and round up to nearest multiple of 16
  #endif

  #ifdef USE_AVX
	#include "radix4032_avx_negadwt_consts.h"
  #endif
	#include "sse2_macro_gcc64.h"
	#include "radix09_sse_macro.h"

#endif	// USE_SSE2

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
		int wts_idx_inc2;
		int icycle[ODD_RADIX];
	#ifdef USE_SSE2
		int jcycle[ODD_RADIX];
	#endif
	#ifdef USE_AVX
		int kcycle[ODD_RADIX];
		int lcycle[ODD_RADIX];
	#endif
	#ifdef USE_AVX512
		int mcycle[ODD_RADIX];
		int ncycle[ODD_RADIX];
		int ocycle[ODD_RADIX];
		int pcycle[ODD_RADIX];
	#endif

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

int radix4032_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-4032 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-4032 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix4032_ditN_cy_dif1";
	static int thr_id = 0;	// Master thread gets this special id
	const int pfetch_dist = PFETCH_DIST;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
	static double wts_mult[2], inv_mult[2];	// Const wts-multiplier and 2*(its multiplicative inverse)
	double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	int NDIVR,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,ntmp,outer,nbytes;
	uint32 nwt16 = nwt << L2_SZ_VD;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code

	// Need these both in scalar mode and to ease the SSE2-array init...dimension = ODD_RADIX;
	// In order to ease the ptr-access for the || routine, lump these 4*ODD_RADIX doubles together with copies of
	// the 4 in the passed-in bs[2] and bsinv[2] arrays [and used in this 4-double form by the mersenne-mod carry macros]
	// into a single foo_array[4*(ODD_RADIX+1)], then convert what used to be disparate ODD_RADIX-sized arrays to pointers.
	static double foo_array[(ODD_RADIX+1)<<2], *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr, bs,bsinv;

	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4, nsave = 0;
	static int poff[RADIX>>2];
#ifndef MULTITHREAD
// Shared by DIF and DIT:
	int kk, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf;
	const uint8 *iptr;
	static int io_offsets[64], t_offsets[64];
	static int plo[64], phi[ODD_RADIX], toff[ODD_RADIX];
	const uint8 dft_p40_cperms[128] = {	// Pad byte-array to next-higher 8-byte multiple
		0,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01,
		0,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01
	},	dft_phi[64] = {	// Only 63 entries, but again pad out to nearest 8-multiple
		0,0x02,0x01,0x08,0x07,0x06,0x05,0x04,0x03,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09
	};
// DIF:
	/*** Now done via #include "radix4032.h" @top of this file ***/	// 4032-element dif64_oidx_lo byte-array
// DIT:
	// Low parts [p0-3f] of output-index perms - need 64 bytes (one init row below) for each radix-64 DFT:
	const uint8 dit64_iidx_lo[64] = {	// Keep only first row...CF. comments in radix4032_dit_pass1().
		0x00,0x01,0x03,0x02,0x07,0x06,0x05,0x04,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x3f,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20
	};
#endif
	static double radix_inv, n2inv;
	double scale, dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	double *addr, *addi;
	struct complex t[RADIX], *tptr;
	int *itmp,*itm2;	// Pointer into the bjmodn array
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	static int n_div_nwt;
	int col,co2,co3,m,m2;
  #ifdef USE_AVX512
	double t0,t1,t2,t3;
	static struct uint32x8 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #elif defined(USE_AVX)
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
  #endif
	double rt,it;
	// indices into weights arrays (mod NWT):
	static int ii[ODD_RADIX] = {
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1
	};
	/* These are used in conjunction with the langth-ODD_RADIX arrays in the USE_SCALAR_CARRY code flow;
	In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
	static int idx_offset, idx_incr, wts_idx_incr = 0, wts_idx_inc2 = 0
		,icycle[ODD_RADIX],ic_idx;
#ifdef USE_SSE2
	static int jcycle[ODD_RADIX],jc_idx;
  #ifdef USE_AVX
	static int kcycle[ODD_RADIX],kc_idx;
	static int lcycle[ODD_RADIX],lc_idx;
  #endif
  #ifdef USE_AVX512
	static int mcycle[ODD_RADIX],mc_idx;
	static int ncycle[ODD_RADIX],nc_idx;
	static int ocycle[ODD_RADIX],oc_idx;
	static int pcycle[ODD_RADIX],pc_idx;
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

	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
	const double crnd = 3.0*0x4000000*0x2000000;
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	vec_dbl
		*tmp,*tm0,*tm1,*tm2;	// Non-static utility ptrs
	static vec_dbl *max_err, *sse2_rnd, *half_arr,
		*r00,	// Head of RADIX*vec_cmplx-sized local store #1
		*s1p00,	// Head of RADIX*vec_cmplx-sized local store #2
		*cy_r,*cy_i;	// Need RADIX total vec_dbl slots for sse2 carries, RADIX/2 for avx
  #ifdef USE_AVX
	static vec_dbl *base_negacyclic_root;
  #endif

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
	static task_control_t   task_control = {NULL, (void*)cy4032_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
  #if PFETCH
	double *addp;
  #endif
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

	foo_array[0] = base[0];
	foo_array[1] = base[1];
	foo_array[2] = baseinv[0];
	foo_array[3] = baseinv[1];
	wt_arr    = foo_array + 4;
	wtinv_arr = wt_arr    + ODD_RADIX;
	bs_arr    = wtinv_arr + ODD_RADIX;
	bsinv_arr = bs_arr    + ODD_RADIX;

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
	NDIVR   = n/RADIX;	ndivr_inv = (double)RADIX/n;
	n_div_nwt = NDIVR >> nwt_bits;
//	printf("CY-step: n_div_nwt = %u\n",n_div_nwt);
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

	  #ifdef USE_AVX	// AVX LOACC: Make CARRY_8_WAY default here:
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
		// consisting of 128*2 vec_dbl and (8+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix4032_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix4032_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		// RADIX = 0xfc0, each of the r00 and s1p00 local-store blocks has 2*RADIX vec_dbl elts:
		tmp = sc_ptr;	r00   = tmp;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x1f80;	s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
		tmp += 0x1f80;
		// DFT-63 and -64 roots all def'd locally in the DFT-63,64 functions
	  #ifdef USE_AVX512
		cy_r = tmp;	cy_i = tmp+0x1f8;	tmp += 2*0x1f8;	// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
		cy_r = tmp;	cy_i = tmp+0x3f0;	tmp += 2*0x3f0;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// 0x3f00 += 0x7e0 + 2 => sc_ptr += 0x46e2
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 68 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy_r = tmp;	cy_i = tmp+0x7e0;	tmp += 2*0x7e0;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// 0x3f00 += 0xfc0 + 2 => sc_ptr += 0x4ec2
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*ODD_RADIX] x 16 for Fermat-mod */
	  #endif
		ASSERT(half_arr_offset4032 == (uint32)(half_arr-sc_ptr), "half_arr_offset mismatches actual!");
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
			j = (1<<(2*(L2_SZ_VD-2))) + 4;	// 16+4 for sse2, 64+4 for avx
		} else {
			j = ODD_RADIX<<2;				// 4*ODD_RADIX
		}
		ASSERT((radix4032_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r00) + (j << L2_SZ_VD), "radix4032_creals_in_local_store checksum failed!");

		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
		VEC_DBL_INIT(sse2_rnd, crnd);

	// Init-mode calls to these functions which maintain an internal local-alloc static store:
		thr_id = -1;	// Use this special thread id for any macro-required thread-local-data inits...
		SSE2_RADIX_63_DIF(CY_THREADS,thr_id, 0,0,0,0);
		SSE2_RADIX_63_DIT(CY_THREADS,thr_id, 0,0,0,0);
		SSE2_RADIX_64_DIF(CY_THREADS,thr_id, 0,0,0,0,0,0);
		SSE2_RADIX_64_DIT(CY_THREADS,thr_id, 0,0,0,0,0);
		thr_id = 0;	// ...then revert to 0.

		// Propagate the above consts to the remaining threads:
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
	#ifdef USE_AVX

		base_negacyclic_root = half_arr + RADIX;

	  #if HIACC
		/*
		The pattern of the negacyclic-DWT-weights ("nDWTs") applied to the RADIX complex outputs of the final-DIT-pass is like so:
		The nDWTs multiplying each set of RADIX DIT DFT outputs are simply the product of a single complex-root "base multiplier" rbase
		(separately computed for each batch of DFT outputs), which "base root" multiplies the (0 - RADIX-1)st [4*RADIX]th roots of unity,
		i.e.
			 rbase * (j*I*Pi/2)/RADIX, for j = 0, ..., RADIX-1 .

		See the radix28 version of this routine for additional details.
		*/
		#if 0
		// Simple qfloat-based loop to crank out the roots which populate the radix4032_avx_negadwt_consts table:
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

		tmp = base_negacyclic_root + RADIX*2;	// First 2*RADIX slots reserved for RADIX/4 copies of the Re/Im parts of the 4 base multipliers
		tm2 = tmp + RADIX/2 - 1;
		// First elt-pair needs special handling - have the 1.0 in avx_negadwt_consts[0] but the sine term buggers things
		tmp->d0 = 1.0;	(tmp+1)->d0 = 0.0;
		tmp64 = radix4032_avx_negadwt_consts[1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(  1*I*Pi/(2*RADIX)) = sin((RADIX-  1)*I*Pi/(2*RADIX)) */
		tmp64 = radix4032_avx_negadwt_consts[2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(  2*I*Pi/(2*RADIX)) = sin((RADIX-  2)*I*Pi/(2*RADIX)) */
		tmp64 = radix4032_avx_negadwt_consts[3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(  3*I*Pi/(2*RADIX)) = sin((RADIX-  3)*I*Pi/(2*RADIX)) */	tmp += 2;
		for(j = 4; j < RADIX; j += 4) {
			tmp64 = radix4032_avx_negadwt_consts[j+0];	tmp->d0 = tm2->d0 = *(double *)&tmp64;	tm2 -= 2;
			tmp64 = radix4032_avx_negadwt_consts[j+1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;
			tmp64 = radix4032_avx_negadwt_consts[j+2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;
			tmp64 = radix4032_avx_negadwt_consts[j+3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	tmp += 2;
		}
		tmp = base_negacyclic_root + RADIX*2;	// reset to point to start of above block
		nbytes = RADIX*SZ_VD/2;	// RADIX/4 AVX-register-sized complex data

	  #else	// HIACC = false:

		// lower-precision version, which yields slightly more roundoff error, but is simpler and more storage-compact.
	   #ifdef USE_AVX512

		// Init exp(j*I*Pi/2/RADIX), for j = 0-7:
		tmp = base_negacyclic_root + 16;	// First 16 slots reserved for Re/Im parts of the 8 base multipliers
		tmp->d0 = 1.0;
		tmp64 = radix4032_avx_negadwt_consts[      1];	tmp->d1 = *(double *)&tmp64;	// cos(01*I*Pi/(2*RADIX))
		tmp64 = radix4032_avx_negadwt_consts[      2];	tmp->d2 = *(double *)&tmp64;	// cos(02*I*Pi/(2*RADIX))
		tmp64 = radix4032_avx_negadwt_consts[      3];	tmp->d3 = *(double *)&tmp64;	// cos(03*I*Pi/(2*RADIX))
		tmp64 = radix4032_avx_negadwt_consts[      4];	tmp->d4 = *(double *)&tmp64;	// cos(04*I*Pi/(2*RADIX))
		tmp64 = radix4032_avx_negadwt_consts[      5];	tmp->d5 = *(double *)&tmp64;	// cos(05*I*Pi/(2*RADIX))
		tmp64 = radix4032_avx_negadwt_consts[      6];	tmp->d6 = *(double *)&tmp64;	// cos(06*I*Pi/(2*RADIX))
		tmp64 = radix4032_avx_negadwt_consts[      7];	tmp->d7 = *(double *)&tmp64;	// cos(07*I*Pi/(2*RADIX))
		++tmp;
		tmp->d0 = 0.0;
		tmp64 = radix4032_avx_negadwt_consts[RADIX-1];	tmp->d1 = *(double *)&tmp64;	// sin(01*I*Pi/(2*RADIX))
		tmp64 = radix4032_avx_negadwt_consts[RADIX-2];	tmp->d2 = *(double *)&tmp64;	// sin(02*I*Pi/(2*RADIX))
		tmp64 = radix4032_avx_negadwt_consts[RADIX-3];	tmp->d3 = *(double *)&tmp64;	// sin(03*I*Pi/(2*RADIX))
		tmp64 = radix4032_avx_negadwt_consts[RADIX-4];	tmp->d4 = *(double *)&tmp64;	// sin(04*I*Pi/(2*RADIX))
		tmp64 = radix4032_avx_negadwt_consts[RADIX-5];	tmp->d5 = *(double *)&tmp64;	// sin(05*I*Pi/(2*RADIX))
		tmp64 = radix4032_avx_negadwt_consts[RADIX-6];	tmp->d6 = *(double *)&tmp64;	// sin(06*I*Pi/(2*RADIX))
		tmp64 = radix4032_avx_negadwt_consts[RADIX-7];	tmp->d7 = *(double *)&tmp64;	// sin(07*I*Pi/(2*RADIX))
		++tmp;	// 0x480(base_negacyclic_root)
		tmp64 = radix4032_avx_negadwt_consts[      8];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// cos(08*I*Pi/(2*RADIX))
		++tmp;	// 0x4c0(base_negacyclic_root)
		tmp64 = radix4032_avx_negadwt_consts[RADIX-8];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// sin(08*I*Pi/(2*RADIX))
		tmp = base_negacyclic_root + 16;	// reset to point to start of above block

	   #elif defined(USE_AVX)

		// Init exp(j*I*Pi/2/RADIX), for j = 0-3:
		tmp = base_negacyclic_root + 8;	// First 8 slots reserved for Re/Im parts of the 4 base multipliers
		tmp->d0 = 1.0;
		tmp64 = radix4032_avx_negadwt_consts[      1];	tmp->d1 = *(double *)&tmp64;	// cos(01*I*Pi/(2*RADIX))
		tmp64 = radix4032_avx_negadwt_consts[      2];	tmp->d2 = *(double *)&tmp64;	// cos(02*I*Pi/(2*RADIX))
		tmp64 = radix4032_avx_negadwt_consts[      3];	tmp->d3 = *(double *)&tmp64;	// cos(03*I*Pi/(2*RADIX))

		(++tmp)->d0 = 0.0;
		tmp64 = radix4032_avx_negadwt_consts[RADIX-1];	tmp->d1 = *(double *)&tmp64;	// sin(01*I*Pi/(2*RADIX))
		tmp64 = radix4032_avx_negadwt_consts[RADIX-2];	tmp->d2 = *(double *)&tmp64;	// sin(02*I*Pi/(2*RADIX))
		tmp64 = radix4032_avx_negadwt_consts[RADIX-3];	tmp->d3 = *(double *)&tmp64;	// sin(03*I*Pi/(2*RADIX))
		++tmp;
		tmp64 = radix4032_avx_negadwt_consts[      4];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// cos(04*I*Pi/(2*RADIX))
		++tmp;
		tmp64 = radix4032_avx_negadwt_consts[RADIX-4];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// sin(04*I*Pi/(2*RADIX))
		tmp = base_negacyclic_root + 8;	// reset to point to start of above block

	   #endif

		nbytes = 4*SZ_VD;	// 2 SIMD-complex data

	  #endif	// HIACC toggle

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
		// Each lookup-catagory in the 'mini-tables' used in AVX mode balloons from 16x32-bytes to 64x64-bytes,
		// so switch to an opmask-based scheme which starts with e.g. a broadcast constant and onditional doubling
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
	n_minus_sil   = (struct uint32x8 *)sse_n + 1;
	n_minus_silp1 = (struct uint32x8 *)sse_n + 2;
	sinwt         = (struct uint32x8 *)sse_n + 3;
	sinwtm1       = (struct uint32x8 *)sse_n + 4;
	nbytes += 128;
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

	/*********** Defer the per-thread local-mem-block copy until after added wts-index precomputation below ************/
	#endif	/* USE_SSE2 */

		/*   constant index offsets for array load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		// Legacy code needs 3 lowest nonzero fixed-index p[] terms:
		p1 = 1*NDIVR;	 p1 += ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = 2*NDIVR;	 p2 += ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = 3*NDIVR;	 p3 += ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = 4*NDIVR;	 p4 += ( (p4 >> DAT_BITS) << PAD_BITS );
	#ifndef USE_SSE2
		p0123[0] = 0; p0123[1] = p1; p0123[2] = p2; p0123[3] = p3;
	#endif
		for(l = 0; l < (RADIX>>2); l++) {
			poff[l] = (l<<2)*NDIVR;	// Corr. to every 4th plo[] term
			poff[l] += ( (poff[l] >> DAT_BITS) << PAD_BITS );
		}

	#ifndef MULTITHREAD
		/* Constant padded-array index offsets for array load/stores are here. */
		for(l = 0; l < 64; l++)
		{
			plo[l] = l*NDIVR;
			plo[l] += ( (plo[l] >> DAT_BITS) << PAD_BITS );
			t_offsets[l] = ((l<<6)-l)<<1;	// 2*(l*63), extra 2x is due to cast-to-double of t-array in call to RADIX_64_DIF|T
		}
		for(l = 0; l < ODD_RADIX; l++)
		{
			phi[l] = (l<<6)*NDIVR;
			phi[l] += ( (phi[l] >> DAT_BITS) << PAD_BITS );
			toff[l] = l+l;
		}
	#endif

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

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/radix-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"ERROR: unable to allocate array _bjmodnini in %s.\n", func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		_bjmodnini[0] = 0;
		_bjmodnini[1] = 0;

		/* For Fermat-mod, since 'adjacent' words are actually stride-2 separated
		in terms of the floating residue array, block boundaries have half the i-index
		(e.g. as in sw*i and bw*i) value they do in the Mersenne-mod case:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			jhi = NDIVR/CY_THREADS;
		}
		else
		{
			jhi = NDIVR/CY_THREADS/2;
		}

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

		// In non-power-of-2-runlength case, both Mersenne and Fermat-mod share these next 2 loops:
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

			// Every (ODD_RADIX)th bjmodn initializer needs to be forced-to-bigword in fermat-mod DWT case:
			if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
			{
				/* Start this value off at N in Fermat-mod case, so (bjmodn >= sw) check in
				fermat_carry_norm_errcheck (cf. carry.h) yields a bigword (i == 1) for j= 0:
				*/
				for(i = 0; i < RADIX; i += ODD_RADIX) {
					_bjmodn[i][ithread] = n;
				}
			}
		}

		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		{
			/* For Fermat-mod, IBDWT access patterns repeat with period NWT = {odd part of radix0},
			so for even radix0 values still only need [radix0 >> trailz(radix0)] bjmodn and ii's:
			*/
			/* indices into IBDWT weights arrays (mod NWT) is here: */
			ii[0]= 0;
			ii[1]= (SW_DIV_N*(NDIVR >> 1)) % nwt;	// nwt *not* a power of 2, must use library-mod!
			for(i = 2; i < ODD_RADIX; i++) {
				MOD_ADD32(ii[i-1], ii[1], nwt, ii[i]);
			}

			/* Find the circular-index-shift (cf. the head-of-file comments of radix28_ditN_cy_dif1.c) by searching bjmodn01 ... bjmodn[nwt] for the one == bw: */
			for(i = 1; i < ODD_RADIX; i++) {
				if( _bjmodn[i][0] == bw ) {
					wts_idx_incr = i;
					break;
				};
			}	//	printf("wts_idx_incr = %u\n",wts_idx_incr);
			ASSERT(wts_idx_incr != 0, "wts_idx_incr init failed!");

		#ifdef USE_SSE2
			wts_idx_inc2 = wts_idx_incr << (2*L2_SZ_VD - 3);	/* In the SIMD version, use icycle0-6 as actual address
							offsets, so wts_idx_incr includes a *sizeof(vec_dbl) for the array-of-vector-doubles indexing, and another
							doubling|quadrupling|... to reflect the fact that the SIMD version of the loop is equivalent to 2|4|... scalar
							loop executions, i.e. corresponds to [#doubles in each vec_dbl] scalar-code increments of the icycle indices. */
			wts_idx_inc2 %= nwt16;	/* Need an extra mod since [2|4|...]*wts_idx_incr may be >= nwt */
		#endif
			/* Subtract nwt from the increments to ease fast-mod */
			wts_idx_incr -= nwt;
		#ifdef USE_SSE2
			wts_idx_inc2 -= nwt16;
		#endif

			for(i = 0; i < ODD_RADIX; i++) {
				/* Need this both in scalar mode and to ease the SSE2-array init */
				j = _bjmodn[i][0] > sw;	bs_arr[i] = base[j];	bsinv_arr[i] = baseinv[j];
				wt_arr[i] = wt0[ii[i]];	// inverse wts must be reinited on each pass, since these have a *scale multiplier
				/* Give icycle indices their proper starting values: */
				icycle[i] = i;
			}
			/*
				Ex: wts_idx_incr = 8, cols have contents if various *cycle[] arrays:
				i	j	k	l	m	n	o	p
				0	8	1	9	2	10	3	11
				1	9	2	10	3	11	4	12
				2	10	3	11	4	12	5	13	<*** For base[]/base_inv[] vectors, access is rowwise, e.g. will see [i-l]cycleA = 2,10,3,11
				3	11	4	12	5	13	6	14
				4	12	5	13	6	14	7	0	<*** AVX: successive calls to 4-way cy-macro incr row by 4, e.g. [i-l]cycleA on call[0]: 0,8,1,9, call[1]: 4,12,5,13
				5	13	6	14	7	0	8	1
				6	14	7	0	8	1	9	2	<AVX-512: successive calls to 8-way cy-macro incr row by 8
				7	0	8	1	9	2	10	3
				8	1	9	2	10	3	11	4
				9	2	10	3	11	4	12	5
				10	3	11	4	12	5	13	6
				11	4	12	5	13	6	14	7
				12	5	13	6	14	7	0	8
				13	6	14	7	0	8	1	9
				14	7	0	8	1	9	2	10
			*/
		#ifdef USE_SSE2
			tmp = half_arr;
			for(i = 0; i < ODD_RADIX; i++, tmp++) {
																									tmp->d0 = wt_arr[icycle[i]];
			/* Now set the imaginary parts to the values corresponding to the 2nd of each pair of scalar-mode loop passes.
			Use this sequence for mod-add, as it is faster than general-mod '% nwt': */
				jcycle[i] = icycle[i] + wts_idx_incr;	jcycle[i] += ( (-(jcycle[i] < 0)) & nwt);	tmp->d1 = wt_arr[jcycle[i]];
		  #ifdef USE_AVX
				kcycle[i] = jcycle[i] + wts_idx_incr;	kcycle[i] += ( (-(kcycle[i] < 0)) & nwt);	tmp->d2 = wt_arr[kcycle[i]];
				lcycle[i] = kcycle[i] + wts_idx_incr;	lcycle[i] += ( (-(lcycle[i] < 0)) & nwt);	tmp->d3 = wt_arr[lcycle[i]];
		  #endif
		  #ifdef USE_AVX512
				mcycle[i] = lcycle[i] + wts_idx_incr;	mcycle[i] += ( (-(mcycle[i] < 0)) & nwt);	tmp->d4 = wt_arr[mcycle[i]];
				ncycle[i] = mcycle[i] + wts_idx_incr;	ncycle[i] += ( (-(ncycle[i] < 0)) & nwt);	tmp->d5 = wt_arr[ncycle[i]];
				ocycle[i] = ncycle[i] + wts_idx_incr;	ocycle[i] += ( (-(ocycle[i] < 0)) & nwt);	tmp->d6 = wt_arr[ocycle[i]];
				pcycle[i] = ocycle[i] + wts_idx_incr;	pcycle[i] += ( (-(pcycle[i] < 0)) & nwt);	tmp->d7 = wt_arr[pcycle[i]];
		  #endif
			}

			// Propagate the above wts-consts to the remaining threads:
			nbytes = ODD_RADIX*SZ_VD;
			tmp = half_arr;
			tm2 = tmp + cslots_in_local_store;
			for(ithread = 1; ithread < CY_THREADS; ++ithread) {
				memcpy(tm2, tmp, nbytes);
				tmp = tm2;		tm2 += cslots_in_local_store;
			}

			tmp = half_arr + ODD_RADIX*2;	/* Put the base-mini-arrays right after the weights */

		  #ifdef USE_AVX512

			// Each transposed-data octet in the AVX-512 carry macro needs linearly incrementing bs_arr data (mod ODD_RADIX);
			// Need all [ODD_RADIX] possible such length-8 index subsequences, which will be accessed via their head element
			// by the [ijklmnop]cycle* index octets in the respective carry-macro call:
			tm2 = tmp + ODD_RADIX;
			for(i = 0; i < ODD_RADIX; i++, tmp++, tm2++) {
				tmp->d0 = bs_arr   [i];	tmp->d1 = bs_arr   [(i+1)%ODD_RADIX];	tmp->d2 = bs_arr   [(i+2)%ODD_RADIX];	tmp->d3 = bs_arr   [(i+3)%ODD_RADIX];	tmp->d4 = bs_arr   [(i+4)%ODD_RADIX];	tmp->d5 = bs_arr   [(i+5)%ODD_RADIX];	tmp->d6 = bs_arr   [(i+6)%ODD_RADIX];	tmp->d7 = bs_arr   [(i+7)%ODD_RADIX];
				tm2->d0 = bsinv_arr[i];	tm2->d1 = bsinv_arr[(i+1)%ODD_RADIX];	tm2->d2 = bsinv_arr[(i+2)%ODD_RADIX];	tm2->d3 = bsinv_arr[(i+3)%ODD_RADIX];	tm2->d4 = bsinv_arr[(i+4)%ODD_RADIX];	tm2->d5 = bsinv_arr[(i+5)%ODD_RADIX];	tm2->d6 = bsinv_arr[(i+6)%ODD_RADIX];	tm2->d7 = bsinv_arr[(i+7)%ODD_RADIX];
			}

		  #elif defined(USE_AVX)

			// Each transposed-data quartet in the AVX carry macro needs linearly incrementing bs_arr data (mod ODD_RADIX);
			// Need all [ODD_RADIX] possible such length-4 index subsequences, which will be accessed via their head element
			// by the [ijkl]cycle* index quartets in the respective carry-macro call:
			tm2 = tmp + ODD_RADIX;
			for(i = 0; i < ODD_RADIX; i++, tmp++, tm2++) {
				tmp->d0 = bs_arr   [i];	tmp->d1 = bs_arr   [(i+1)%ODD_RADIX];	tmp->d2 = bs_arr   [(i+2)%ODD_RADIX];	tmp->d3 = bs_arr   [(i+3)%ODD_RADIX];
				tm2->d0 = bsinv_arr[i];	tm2->d1 = bsinv_arr[(i+1)%ODD_RADIX];	tm2->d2 = bsinv_arr[(i+2)%ODD_RADIX];	tm2->d3 = bsinv_arr[(i+3)%ODD_RADIX];
			}

		  #else

			/* In SSE2 mode, because we apply doubled weights to data arranged as [a.re,b.re,...],[a.im,b.im,...] but apply
			doubled base multipliers to shuffled data [a.re,a.im],[b.re,b.im],... (i.e. shuffled to yield same data layout as
			in the scalar case), the weights need to have disparate real and imag parts, but the base/baseinv terms do not: */
			for(i = 0; i < ODD_RADIX; i++) {
				VEC_DBL_INIT(tmp, bs_arr[i]);	++tmp;
			}
			for(i = 0; i < ODD_RADIX; i++) {
				VEC_DBL_INIT(tmp, bsinv_arr[i]);	++tmp;
			}

		  #endif

			// Propagate the above consts to the remaining threads:
			nbytes <<= 1;	// [base+binv] ==> 2x as many consts as [wts], since the wtinv data done each pass of outer-loop
			tmp = half_arr + ODD_RADIX*2;	/* Put the base-mini-arrays right after the weights */
			tm2 = tmp + cslots_in_local_store;
			for(ithread = 1; ithread < CY_THREADS; ++ithread) {
				memcpy(tm2, tmp, nbytes);
				tmp = tm2;		tm2 += cslots_in_local_store;
			}

			for(i = 0; i < ODD_RADIX; i++) {
				icycle[i] <<= L2_SZ_VD;	jcycle[i] <<= L2_SZ_VD;
			#ifdef USE_AVX
				kcycle[i] <<= L2_SZ_VD;	lcycle[i] <<= L2_SZ_VD;
			#endif
			#ifdef USE_AVX512
				mcycle[i] <<= L2_SZ_VD;	ncycle[i] <<= L2_SZ_VD;	ocycle[i] <<= L2_SZ_VD;	pcycle[i] <<= L2_SZ_VD;
			#endif
			}

		#endif	// USE_SSE2 ?
		}	// MODULUS_TYPE_FERMAT ?

	#ifdef USE_PTHREAD

		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			tdat[ithread].bjmodnini = _bjmodnini[CY_THREADS];
			tdat[ithread].bjmodn0 = _bjmodnini[ithread];
		#ifdef USE_SSE2
			tdat[ithread].r00      = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (vec_dbl *)((intptr_t)tdat[ithread].r00 + ((intptr_t)half_arr - (intptr_t)r00));
		#else
			// In scalar mode use these 2 ptrs to pass wts_idx_incr and the base/baseinv/etc array-ptrs:
			tdat[ithread].r00      = (double *)foo_array;
			tdat[ithread].half_arr = (double *)&wts_idx_incr;
		#endif	// USE_SSE2
		}

		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		{
			// These inits must occur just once, in loop-simulated full-pass mode,
			// in order to get the array-index-offset values of the icycle/jcycle indices right:
			for(i = 0; i < ODD_RADIX; i++) {
				tdat[0].icycle[i] = icycle[i];
			  #ifdef USE_SSE2
				tdat[0].wts_idx_inc2 = wts_idx_inc2;
				tdat[0].jcycle[i] = jcycle[i];
			  #endif
			  #ifdef USE_AVX
				tdat[0].kcycle[i] = kcycle[i];
				tdat[0].lcycle[i] = lcycle[i];
			  #endif
			  #ifdef USE_AVX512
				tdat[0].mcycle[i] = mcycle[i];
				tdat[0].ncycle[i] = ncycle[i];
				tdat[0].ocycle[i] = ocycle[i];
				tdat[0].pcycle[i] = pcycle[i];
			  #endif
			}
			// For remaining threads, simulate the loop-evolution of the above indices.
			// Note that the non-thread-associated *cycle[] arry data will get changed fom their above-inited
			// values in the loop here, but that's OK because in || mode only the thread-associated values matter:
			for(ithread = 1; ithread < CY_THREADS; ++ithread) {
				jstart = 0;
				jhi = NDIVR/CY_THREADS;	// Earlier setting = NDIVR/CY_THREADS/2 was for simulating bjmodn evolution, need 2x that here
				// Get value of (negative) increment resulting from (jhi-jstart)/stride execs of *cycle[] += wts_idx_inc* (mod nwt*):
			#ifndef USE_SSE2
				j = ((int64)wts_idx_incr * ( (jhi-jstart)>>(L2_SZ_VD-2) ) % nwt16);	// []>>(L2_SZ_VD-2) is fast subst. for []/stride
			#else
				j = ((int64)wts_idx_inc2 * ( (jhi-jstart)>>(L2_SZ_VD-2) ) % nwt16);	// Cast wts_idx_inc* to signed 64-bit to avoid
						// overflow of product; further compute (jhi-jstart)/stride prior to multiply to gain more bits-to-spare.
			#endif
				// khi = 1 for Fermat-mod, thus no outer loop needed here
				for(i = 0; i < ODD_RADIX; i++) {
				#ifndef USE_SSE2	// Scalar-double mode uses non-pointerized icycle values:
					icycle[i] += j;		icycle[i] += ( (-(int)((uint32)icycle[i] >> 31)) & nwt);
				#else
					icycle[i] += j;		icycle[i] += ( (-(int)((uint32)icycle[i] >> 31)) & nwt16);
					jcycle[i] += j;		jcycle[i] += ( (-(int)((uint32)jcycle[i] >> 31)) & nwt16);
				  #ifdef USE_AVX
					kcycle[i] += j;		kcycle[i] += ( (-(int)((uint32)kcycle[i] >> 31)) & nwt16);
					lcycle[i] += j;		lcycle[i] += ( (-(int)((uint32)lcycle[i] >> 31)) & nwt16);
				  #endif
				  #ifdef USE_AVX512
					mcycle[i] += j;		mcycle[i] += ( (-(int)((uint32)mcycle[i] >> 31)) & nwt16);
					ncycle[i] += j;		ncycle[i] += ( (-(int)((uint32)ncycle[i] >> 31)) & nwt16);
					ocycle[i] += j;		ocycle[i] += ( (-(int)((uint32)ocycle[i] >> 31)) & nwt16);
					pcycle[i] += j;		pcycle[i] += ( (-(int)((uint32)pcycle[i] >> 31)) & nwt16);
				  #endif
				#endif
				}
				for(i = 0; i < ODD_RADIX; i++) {
					tdat[ithread].icycle[i] = icycle[i];
				  #ifdef USE_SSE2
					tdat[ithread].wts_idx_inc2 = wts_idx_inc2;
					tdat[ithread].jcycle[i] = jcycle[i];
				  #endif
				  #ifdef USE_AVX
					tdat[ithread].kcycle[i] = kcycle[i];
					tdat[ithread].lcycle[i] = lcycle[i];
				  #endif
				  #ifdef USE_AVX512
					tdat[ithread].mcycle[i] = mcycle[i];
					tdat[ithread].ncycle[i] = ncycle[i];
					tdat[ithread].ocycle[i] = ocycle[i];
					tdat[ithread].pcycle[i] = pcycle[i];
				  #endif
				}
			}
			// Restore the original loop-start values of the cycle arrays, since we use these for init of inv-wts below:
			for(i = 0; i < ODD_RADIX; i++) {
				icycle[i] = tdat[0].icycle[i];
			  #ifdef USE_SSE2
				jcycle[i] = tdat[0].jcycle[i];
			  #endif
			  #ifdef USE_AVX
				kcycle[i] = tdat[0].kcycle[i];
				lcycle[i] = tdat[0].lcycle[i];
			  #endif
			  #ifdef USE_AVX512
				mcycle[i] = tdat[0].mcycle[i];
				ncycle[i] = tdat[0].ncycle[i];
				ocycle[i] = tdat[0].ocycle[i];
				pcycle[i] = tdat[0].pcycle[i];
			  #endif
			}
		}

	#endif	// USE_PTHREAD

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

/*...The radix-4032 final DIT pass is here.	*/

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
		_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (_i[0] = 1).	*/

		khi = n_div_nwt/CY_THREADS;
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_jstart[ithread] = ithread*NDIVR/CY_THREADS;
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + nwt-1;

			_col[ithread] = ithread*(khi*RADIX);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
			_co2[ithread] = (n>>nwt_bits)-1+RADIX - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
			_co3[ithread] = _co2[ithread]-RADIX;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
		}
	}
	else
	{
		_i[0] = 0;		/* Pointer to the BASE and BASEINV arrays. If n divides p, lowest-order digit is always a smallword (_i[0] = 0).	*/
		khi = 1;
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_jstart[ithread] = ithread*NDIVR/CY_THREADS;
			/*
			For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
			*/
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + 15;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + n_div_nwt/CY_THREADS;
		}

		// Now that full_pass-dependent scale factor known, init inverse weights tiny-table used for Fermat-mod
		for(i = 0; i < ODD_RADIX; i++) {
			wtinv_arr[i] = scale*wt1[ii[i]];
		}

	// In threaded mode, the master *cycle[] values are unmodified during main loop exec; only the thread-associated
	// copies of these index arrays get modified. In non-threaded mode we must separately store copies of the masters
	// in order to solve the save/restore issue. We start from the (static, unmodified during loop) ii[]-index values:
	#ifndef MULTITHREAD
		for(i = 0; i < ODD_RADIX; i++) {
			/* Reinit *cycle indices their proper starting values - recall in SIMD mode these all are ( << L2_SZ_VD): */
			icycle[i] = i;
		#ifdef USE_SSE2
			jcycle[i] = icycle[i] + wts_idx_incr;	jcycle[i] += ( (-(jcycle[i] < 0)) & nwt);
		  #ifdef USE_AVX
			kcycle[i] = jcycle[i] + wts_idx_incr;	kcycle[i] += ( (-(kcycle[i] < 0)) & nwt);
			lcycle[i] = kcycle[i] + wts_idx_incr;	lcycle[i] += ( (-(lcycle[i] < 0)) & nwt);
		   #ifdef USE_AVX512
			mcycle[i] = lcycle[i] + wts_idx_incr;	mcycle[i] += ( (-(mcycle[i] < 0)) & nwt);
			ncycle[i] = mcycle[i] + wts_idx_incr;	ncycle[i] += ( (-(ncycle[i] < 0)) & nwt);
			ocycle[i] = ncycle[i] + wts_idx_incr;	ocycle[i] += ( (-(ocycle[i] < 0)) & nwt);
			pcycle[i] = ocycle[i] + wts_idx_incr;	pcycle[i] += ( (-(pcycle[i] < 0)) & nwt);
			mcycle[i] <<= L2_SZ_VD;	ncycle[i] <<= L2_SZ_VD;	ocycle[i] <<= L2_SZ_VD;	pcycle[i] <<= L2_SZ_VD;
		   #endif
			kcycle[i] <<= L2_SZ_VD;	lcycle[i] <<= L2_SZ_VD;
		  #endif
			icycle[i] <<= L2_SZ_VD;	jcycle[i] <<= L2_SZ_VD;
		#endif
		}
	#endif

	#ifdef USE_SSE2
		// Remember: *cycle[] entries all << L2_SZ_VD here - must left-shift-on-the-fly before using:
		tm2 = half_arr + ODD_RADIX;
		for(i = 0; i < ODD_RADIX; i++, tm2++) {
			tm2->d0 = wtinv_arr[icycle[i] >> L2_SZ_VD];
			tm2->d1 = wtinv_arr[jcycle[i] >> L2_SZ_VD];
		#ifdef USE_AVX
			tm2->d2 = wtinv_arr[kcycle[i] >> L2_SZ_VD];
			tm2->d3 = wtinv_arr[lcycle[i] >> L2_SZ_VD];
		#endif
		#ifdef USE_AVX512
			tm2->d4 = wtinv_arr[mcycle[i] >> L2_SZ_VD];
			tm2->d5 = wtinv_arr[ncycle[i] >> L2_SZ_VD];
			tm2->d6 = wtinv_arr[ocycle[i] >> L2_SZ_VD];
			tm2->d7 = wtinv_arr[pcycle[i] >> L2_SZ_VD];
		#endif
		}

		// Propagate the above inv-wts to the remaining threads - surrounding consts are unchanged:
		nbytes = ODD_RADIX*SZ_VD;
		tmp = half_arr + ODD_RADIX;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	#endif	// USE_SSE2 ?
	}	// MODULUS_TYPE_FERMAT ?

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
		ASSERT(tdat[ithread].wts_idx_inc2 == wts_idx_inc2, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].r00 == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].half_arr;
		ASSERT(((tmp-1)->d0 == crnd && (tmp-1)->d1 == crnd), "thread-local memcheck failed!");
	#endif
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX
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
		#ifdef USE_SSE2
			dtmp = (tmp)->d0 * (tmp+ODD_RADIX)->d0;	ASSERT(fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp)->d1 * (tmp+ODD_RADIX)->d1;	ASSERT(fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
		#endif
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
		#include "radix4032_main_carry_loop.h"

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
		ASSERT(0x0 == cy4032_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
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
	//	printf("Iter = %d, cy0,1 =  %20.15f, %20.15f, maxerr = %20.15f\n",iter,_cy_r[0][0],_cy_r[1][0],maxerr);
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
	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		j_jhi =15;
	}
	else
	{
		j_jhi = 7;
	}

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

void radix4032_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-4032 = 63x64 complex DIF FFT pass on the data in the length-N real vector A.
*/
	int k,l, j,j1,j2,jt;
	static int NDIVR,first_entry=TRUE;
	struct complex t[RADIX], *tptr;
	// Need storage for circular-shifts perms of a basic 63-vector, with shift count in [0,63] that means 2*63 elts:
	const uint8 *iptr, dif_p40_cperms[128] = {	// Pad byte-array to next-higher 8-byte multiple
		0,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01,
		0,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01
	},	dif_phi[64] = {	// Only 63 entries, but again pad out to nearest 8-multiple
		0,0x02,0x01,0x08,0x07,0x06,0x05,0x04,0x03,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09
	};
	static int t_offsets[64], io_offsets[64];
	static int plo[64], phi[ODD_RADIX], toff[ODD_RADIX];
	/*** Now done via #include "radix4032.h" @top of this file***/	// 4032-element dif64_oidx_lo byte-array

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR = n/RADIX;
		/* Constant padded-array index offsets for array load/stores are here. */
		for(l = 0; l < 64; l++)
		{
			plo[l] = l*NDIVR;
			plo[l] += ( (plo[l] >> DAT_BITS) << PAD_BITS );
			t_offsets[l] = ((l<<6)-l)<<1;	// 2*(l*63), extra 2x is due to cast-to-double of t-array in call to RADIX_64_DIF
		//	o_offsets[l] = plo[l];	// Initial pass (which we use to auto-gen the needed operm) values of o_offsets
		}
		for(l = 0; l < ODD_RADIX; l++)
		{
			phi[l] = (l<<6)*NDIVR;
			phi[l] += ( (phi[l] >> DAT_BITS) << PAD_BITS );
			toff[l] = l+l;
		}
	}

/*...The radix-4032 pass is here.	*/

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

	//...gather the needed data (4032 64-bit complex, i.e 8064 64-bit reals) and do 64 radix-63 transforms...
	/*
	Twiddleless version arranges 64 sets of radix-63 DIF inputs as follows: 0 in upper left corner,
	decrement 64 horizontally and 63 vertically, all resulting indexing modulo 4032 ( = 0xfc0):
	(we can auto-generate these by compiling test_fft_radix.c with -DTTYPE=0 -DRADIX=4032, running
	the resulting executable and snarfing the first set of index-outputs, "DIF/DIT input-scramble array").

	DIF/DIT input-scramble array = [in hex]
		000,f80,f40,f00,ec0,e80,e40,e00,dc0,d80,d40,d00,cc0,c80,c40,c00,bc0,b80,b40,b00,ac0,a80,a40,a00,9c0,980,940,900,8c0,880,840,800,7c0,780,740,700,6c0,680,640,600,5c0,580,540,500,4c0,480,440,400,3c0,380,340,300,2c0,280,240,200,1c0,180,140,100,0c0,080,040,
		f81,f41,f01,ec1,e81,e41,e01,dc1,d81,d41,d01,cc1,c81,c41,c01,bc1,b81,b41,b01,ac1,a81,a41,a01,9c1,981,941,901,8c1,881,841,801,7c1,781,741,701,6c1,681,641,601,5c1,581,541,501,4c1,481,441,401,3c1,381,341,301,2c1,281,241,201,1c1,181,141,101,0c1,081,041,001,
		f42,f02,ec2,e82,e42,e02,dc2,d82,d42,d02,cc2,c82,c42,c02,bc2,b82,b42,b02,ac2,a82,a42,a02,9c2,982,942,902,8c2,882,842,802,7c2,782,742,702,6c2,682,642,602,5c2,582,542,502,4c2,482,442,402,3c2,382,342,302,2c2,282,242,202,1c2,182,142,102,0c2,082,042,002,f82,
		f03,ec3,e83,e43,e03,dc3,d83,d43,d03,cc3,c83,c43,c03,bc3,b83,b43,b03,ac3,a83,a43,a03,9c3,983,943,903,8c3,883,843,803,7c3,783,743,703,6c3,683,643,603,5c3,583,543,503,4c3,483,443,403,3c3,383,343,303,2c3,283,243,203,1c3,183,143,103,0c3,083,043,003,f83,f43,
		ec4,e84,e44,e04,dc4,d84,d44,d04,cc4,c84,c44,c04,bc4,b84,b44,b04,ac4,a84,a44,a04,9c4,984,944,904,8c4,884,844,804,7c4,784,744,704,6c4,684,644,604,5c4,584,544,504,4c4,484,444,404,3c4,384,344,304,2c4,284,244,204,1c4,184,144,104,0c4,084,044,004,f84,f44,f04,
		e85,e45,e05,dc5,d85,d45,d05,cc5,c85,c45,c05,bc5,b85,b45,b05,ac5,a85,a45,a05,9c5,985,945,905,8c5,885,845,805,7c5,785,745,705,6c5,685,645,605,5c5,585,545,505,4c5,485,445,405,3c5,385,345,305,2c5,285,245,205,1c5,185,145,105,0c5,085,045,005,f85,f45,f05,ec5,
		e46,e06,dc6,d86,d46,d06,cc6,c86,c46,c06,bc6,b86,b46,b06,ac6,a86,a46,a06,9c6,986,946,906,8c6,886,846,806,7c6,786,746,706,6c6,686,646,606,5c6,586,546,506,4c6,486,446,406,3c6,386,346,306,2c6,286,246,206,1c6,186,146,106,0c6,086,046,006,f86,f46,f06,ec6,e86,
		e07,dc7,d87,d47,d07,cc7,c87,c47,c07,bc7,b87,b47,b07,ac7,a87,a47,a07,9c7,987,947,907,8c7,887,847,807,7c7,787,747,707,6c7,687,647,607,5c7,587,547,507,4c7,487,447,407,3c7,387,347,307,2c7,287,247,207,1c7,187,147,107,0c7,087,047,007,f87,f47,f07,ec7,e87,e47,
		dc8,d88,d48,d08,cc8,c88,c48,c08,bc8,b88,b48,b08,ac8,a88,a48,a08,9c8,988,948,908,8c8,888,848,808,7c8,788,748,708,6c8,688,648,608,5c8,588,548,508,4c8,488,448,408,3c8,388,348,308,2c8,288,248,208,1c8,188,148,108,0c8,088,048,008,f88,f48,f08,ec8,e88,e48,e08,
		d89,d49,d09,cc9,c89,c49,c09,bc9,b89,b49,b09,ac9,a89,a49,a09,9c9,989,949,909,8c9,889,849,809,7c9,789,749,709,6c9,689,649,609,5c9,589,549,509,4c9,489,449,409,3c9,389,349,309,2c9,289,249,209,1c9,189,149,109,0c9,089,049,009,f89,f49,f09,ec9,e89,e49,e09,dc9,
		d4a,d0a,cca,c8a,c4a,c0a,bca,b8a,b4a,b0a,aca,a8a,a4a,a0a,9ca,98a,94a,90a,8ca,88a,84a,80a,7ca,78a,74a,70a,6ca,68a,64a,60a,5ca,58a,54a,50a,4ca,48a,44a,40a,3ca,38a,34a,30a,2ca,28a,24a,20a,1ca,18a,14a,10a,0ca,08a,04a,00a,f8a,f4a,f0a,eca,e8a,e4a,e0a,dca,d8a,
		d0b,ccb,c8b,c4b,c0b,bcb,b8b,b4b,b0b,acb,a8b,a4b,a0b,9cb,98b,94b,90b,8cb,88b,84b,80b,7cb,78b,74b,70b,6cb,68b,64b,60b,5cb,58b,54b,50b,4cb,48b,44b,40b,3cb,38b,34b,30b,2cb,28b,24b,20b,1cb,18b,14b,10b,0cb,08b,04b,00b,f8b,f4b,f0b,ecb,e8b,e4b,e0b,dcb,d8b,d4b,
		ccc,c8c,c4c,c0c,bcc,b8c,b4c,b0c,acc,a8c,a4c,a0c,9cc,98c,94c,90c,8cc,88c,84c,80c,7cc,78c,74c,70c,6cc,68c,64c,60c,5cc,58c,54c,50c,4cc,48c,44c,40c,3cc,38c,34c,30c,2cc,28c,24c,20c,1cc,18c,14c,10c,0cc,08c,04c,00c,f8c,f4c,f0c,ecc,e8c,e4c,e0c,dcc,d8c,d4c,d0c,
		c8d,c4d,c0d,bcd,b8d,b4d,b0d,acd,a8d,a4d,a0d,9cd,98d,94d,90d,8cd,88d,84d,80d,7cd,78d,74d,70d,6cd,68d,64d,60d,5cd,58d,54d,50d,4cd,48d,44d,40d,3cd,38d,34d,30d,2cd,28d,24d,20d,1cd,18d,14d,10d,0cd,08d,04d,00d,f8d,f4d,f0d,ecd,e8d,e4d,e0d,dcd,d8d,d4d,d0d,ccd,
		c4e,c0e,bce,b8e,b4e,b0e,ace,a8e,a4e,a0e,9ce,98e,94e,90e,8ce,88e,84e,80e,7ce,78e,74e,70e,6ce,68e,64e,60e,5ce,58e,54e,50e,4ce,48e,44e,40e,3ce,38e,34e,30e,2ce,28e,24e,20e,1ce,18e,14e,10e,0ce,08e,04e,00e,f8e,f4e,f0e,ece,e8e,e4e,e0e,dce,d8e,d4e,d0e,cce,c8e,
		c0f,bcf,b8f,b4f,b0f,acf,a8f,a4f,a0f,9cf,98f,94f,90f,8cf,88f,84f,80f,7cf,78f,74f,70f,6cf,68f,64f,60f,5cf,58f,54f,50f,4cf,48f,44f,40f,3cf,38f,34f,30f,2cf,28f,24f,20f,1cf,18f,14f,10f,0cf,08f,04f,00f,f8f,f4f,f0f,ecf,e8f,e4f,e0f,dcf,d8f,d4f,d0f,ccf,c8f,c4f,
		bd0,b90,b50,b10,ad0,a90,a50,a10,9d0,990,950,910,8d0,890,850,810,7d0,790,750,710,6d0,690,650,610,5d0,590,550,510,4d0,490,450,410,3d0,390,350,310,2d0,290,250,210,1d0,190,150,110,0d0,090,050,010,f90,f50,f10,ed0,e90,e50,e10,dd0,d90,d50,d10,cd0,c90,c50,c10,
		b91,b51,b11,ad1,a91,a51,a11,9d1,991,951,911,8d1,891,851,811,7d1,791,751,711,6d1,691,651,611,5d1,591,551,511,4d1,491,451,411,3d1,391,351,311,2d1,291,251,211,1d1,191,151,111,0d1,091,051,011,f91,f51,f11,ed1,e91,e51,e11,dd1,d91,d51,d11,cd1,c91,c51,c11,bd1,
		b52,b12,ad2,a92,a52,a12,9d2,992,952,912,8d2,892,852,812,7d2,792,752,712,6d2,692,652,612,5d2,592,552,512,4d2,492,452,412,3d2,392,352,312,2d2,292,252,212,1d2,192,152,112,0d2,092,052,012,f92,f52,f12,ed2,e92,e52,e12,dd2,d92,d52,d12,cd2,c92,c52,c12,bd2,b92,
		b13,ad3,a93,a53,a13,9d3,993,953,913,8d3,893,853,813,7d3,793,753,713,6d3,693,653,613,5d3,593,553,513,4d3,493,453,413,3d3,393,353,313,2d3,293,253,213,1d3,193,153,113,0d3,093,053,013,f93,f53,f13,ed3,e93,e53,e13,dd3,d93,d53,d13,cd3,c93,c53,c13,bd3,b93,b53,
		ad4,a94,a54,a14,9d4,994,954,914,8d4,894,854,814,7d4,794,754,714,6d4,694,654,614,5d4,594,554,514,4d4,494,454,414,3d4,394,354,314,2d4,294,254,214,1d4,194,154,114,0d4,094,054,014,f94,f54,f14,ed4,e94,e54,e14,dd4,d94,d54,d14,cd4,c94,c54,c14,bd4,b94,b54,b14,
		a95,a55,a15,9d5,995,955,915,8d5,895,855,815,7d5,795,755,715,6d5,695,655,615,5d5,595,555,515,4d5,495,455,415,3d5,395,355,315,2d5,295,255,215,1d5,195,155,115,0d5,095,055,015,f95,f55,f15,ed5,e95,e55,e15,dd5,d95,d55,d15,cd5,c95,c55,c15,bd5,b95,b55,b15,ad5,
		a56,a16,9d6,996,956,916,8d6,896,856,816,7d6,796,756,716,6d6,696,656,616,5d6,596,556,516,4d6,496,456,416,3d6,396,356,316,2d6,296,256,216,1d6,196,156,116,0d6,096,056,016,f96,f56,f16,ed6,e96,e56,e16,dd6,d96,d56,d16,cd6,c96,c56,c16,bd6,b96,b56,b16,ad6,a96,
		a17,9d7,997,957,917,8d7,897,857,817,7d7,797,757,717,6d7,697,657,617,5d7,597,557,517,4d7,497,457,417,3d7,397,357,317,2d7,297,257,217,1d7,197,157,117,0d7,097,057,017,f97,f57,f17,ed7,e97,e57,e17,dd7,d97,d57,d17,cd7,c97,c57,c17,bd7,b97,b57,b17,ad7,a97,a57,
		9d8,998,958,918,8d8,898,858,818,7d8,798,758,718,6d8,698,658,618,5d8,598,558,518,4d8,498,458,418,3d8,398,358,318,2d8,298,258,218,1d8,198,158,118,0d8,098,058,018,f98,f58,f18,ed8,e98,e58,e18,dd8,d98,d58,d18,cd8,c98,c58,c18,bd8,b98,b58,b18,ad8,a98,a58,a18,
		999,959,919,8d9,899,859,819,7d9,799,759,719,6d9,699,659,619,5d9,599,559,519,4d9,499,459,419,3d9,399,359,319,2d9,299,259,219,1d9,199,159,119,0d9,099,059,019,f99,f59,f19,ed9,e99,e59,e19,dd9,d99,d59,d19,cd9,c99,c59,c19,bd9,b99,b59,b19,ad9,a99,a59,a19,9d9,
		95a,91a,8da,89a,85a,81a,7da,79a,75a,71a,6da,69a,65a,61a,5da,59a,55a,51a,4da,49a,45a,41a,3da,39a,35a,31a,2da,29a,25a,21a,1da,19a,15a,11a,0da,09a,05a,01a,f9a,f5a,f1a,eda,e9a,e5a,e1a,dda,d9a,d5a,d1a,cda,c9a,c5a,c1a,bda,b9a,b5a,b1a,ada,a9a,a5a,a1a,9da,99a,
		91b,8db,89b,85b,81b,7db,79b,75b,71b,6db,69b,65b,61b,5db,59b,55b,51b,4db,49b,45b,41b,3db,39b,35b,31b,2db,29b,25b,21b,1db,19b,15b,11b,0db,09b,05b,01b,f9b,f5b,f1b,edb,e9b,e5b,e1b,ddb,d9b,d5b,d1b,cdb,c9b,c5b,c1b,bdb,b9b,b5b,b1b,adb,a9b,a5b,a1b,9db,99b,95b,
		8dc,89c,85c,81c,7dc,79c,75c,71c,6dc,69c,65c,61c,5dc,59c,55c,51c,4dc,49c,45c,41c,3dc,39c,35c,31c,2dc,29c,25c,21c,1dc,19c,15c,11c,0dc,09c,05c,01c,f9c,f5c,f1c,edc,e9c,e5c,e1c,ddc,d9c,d5c,d1c,cdc,c9c,c5c,c1c,bdc,b9c,b5c,b1c,adc,a9c,a5c,a1c,9dc,99c,95c,91c,
		89d,85d,81d,7dd,79d,75d,71d,6dd,69d,65d,61d,5dd,59d,55d,51d,4dd,49d,45d,41d,3dd,39d,35d,31d,2dd,29d,25d,21d,1dd,19d,15d,11d,0dd,09d,05d,01d,f9d,f5d,f1d,edd,e9d,e5d,e1d,ddd,d9d,d5d,d1d,cdd,c9d,c5d,c1d,bdd,b9d,b5d,b1d,add,a9d,a5d,a1d,9dd,99d,95d,91d,8dd,
		85e,81e,7de,79e,75e,71e,6de,69e,65e,61e,5de,59e,55e,51e,4de,49e,45e,41e,3de,39e,35e,31e,2de,29e,25e,21e,1de,19e,15e,11e,0de,09e,05e,01e,f9e,f5e,f1e,ede,e9e,e5e,e1e,dde,d9e,d5e,d1e,cde,c9e,c5e,c1e,bde,b9e,b5e,b1e,ade,a9e,a5e,a1e,9de,99e,95e,91e,8de,89e,
		81f,7df,79f,75f,71f,6df,69f,65f,61f,5df,59f,55f,51f,4df,49f,45f,41f,3df,39f,35f,31f,2df,29f,25f,21f,1df,19f,15f,11f,0df,09f,05f,01f,f9f,f5f,f1f,edf,e9f,e5f,e1f,ddf,d9f,d5f,d1f,cdf,c9f,c5f,c1f,bdf,b9f,b5f,b1f,adf,a9f,a5f,a1f,9df,99f,95f,91f,8df,89f,85f,
		7e0,7a0,760,720,6e0,6a0,660,620,5e0,5a0,560,520,4e0,4a0,460,420,3e0,3a0,360,320,2e0,2a0,260,220,1e0,1a0,160,120,0e0,0a0,060,020,fa0,f60,f20,ee0,ea0,e60,e20,de0,da0,d60,d20,ce0,ca0,c60,c20,be0,ba0,b60,b20,ae0,aa0,a60,a20,9e0,9a0,960,920,8e0,8a0,860,820,
		7a1,761,721,6e1,6a1,661,621,5e1,5a1,561,521,4e1,4a1,461,421,3e1,3a1,361,321,2e1,2a1,261,221,1e1,1a1,161,121,0e1,0a1,061,021,fa1,f61,f21,ee1,ea1,e61,e21,de1,da1,d61,d21,ce1,ca1,c61,c21,be1,ba1,b61,b21,ae1,aa1,a61,a21,9e1,9a1,961,921,8e1,8a1,861,821,7e1,
		762,722,6e2,6a2,662,622,5e2,5a2,562,522,4e2,4a2,462,422,3e2,3a2,362,322,2e2,2a2,262,222,1e2,1a2,162,122,0e2,0a2,062,022,fa2,f62,f22,ee2,ea2,e62,e22,de2,da2,d62,d22,ce2,ca2,c62,c22,be2,ba2,b62,b22,ae2,aa2,a62,a22,9e2,9a2,962,922,8e2,8a2,862,822,7e2,7a2,
		723,6e3,6a3,663,623,5e3,5a3,563,523,4e3,4a3,463,423,3e3,3a3,363,323,2e3,2a3,263,223,1e3,1a3,163,123,0e3,0a3,063,023,fa3,f63,f23,ee3,ea3,e63,e23,de3,da3,d63,d23,ce3,ca3,c63,c23,be3,ba3,b63,b23,ae3,aa3,a63,a23,9e3,9a3,963,923,8e3,8a3,863,823,7e3,7a3,763,
		6e4,6a4,664,624,5e4,5a4,564,524,4e4,4a4,464,424,3e4,3a4,364,324,2e4,2a4,264,224,1e4,1a4,164,124,0e4,0a4,064,024,fa4,f64,f24,ee4,ea4,e64,e24,de4,da4,d64,d24,ce4,ca4,c64,c24,be4,ba4,b64,b24,ae4,aa4,a64,a24,9e4,9a4,964,924,8e4,8a4,864,824,7e4,7a4,764,724,
		6a5,665,625,5e5,5a5,565,525,4e5,4a5,465,425,3e5,3a5,365,325,2e5,2a5,265,225,1e5,1a5,165,125,0e5,0a5,065,025,fa5,f65,f25,ee5,ea5,e65,e25,de5,da5,d65,d25,ce5,ca5,c65,c25,be5,ba5,b65,b25,ae5,aa5,a65,a25,9e5,9a5,965,925,8e5,8a5,865,825,7e5,7a5,765,725,6e5,
		666,626,5e6,5a6,566,526,4e6,4a6,466,426,3e6,3a6,366,326,2e6,2a6,266,226,1e6,1a6,166,126,0e6,0a6,066,026,fa6,f66,f26,ee6,ea6,e66,e26,de6,da6,d66,d26,ce6,ca6,c66,c26,be6,ba6,b66,b26,ae6,aa6,a66,a26,9e6,9a6,966,926,8e6,8a6,866,826,7e6,7a6,766,726,6e6,6a6,
		627,5e7,5a7,567,527,4e7,4a7,467,427,3e7,3a7,367,327,2e7,2a7,267,227,1e7,1a7,167,127,0e7,0a7,067,027,fa7,f67,f27,ee7,ea7,e67,e27,de7,da7,d67,d27,ce7,ca7,c67,c27,be7,ba7,b67,b27,ae7,aa7,a67,a27,9e7,9a7,967,927,8e7,8a7,867,827,7e7,7a7,767,727,6e7,6a7,667,
		5e8,5a8,568,528,4e8,4a8,468,428,3e8,3a8,368,328,2e8,2a8,268,228,1e8,1a8,168,128,0e8,0a8,068,028,fa8,f68,f28,ee8,ea8,e68,e28,de8,da8,d68,d28,ce8,ca8,c68,c28,be8,ba8,b68,b28,ae8,aa8,a68,a28,9e8,9a8,968,928,8e8,8a8,868,828,7e8,7a8,768,728,6e8,6a8,668,628,
		5a9,569,529,4e9,4a9,469,429,3e9,3a9,369,329,2e9,2a9,269,229,1e9,1a9,169,129,0e9,0a9,069,029,fa9,f69,f29,ee9,ea9,e69,e29,de9,da9,d69,d29,ce9,ca9,c69,c29,be9,ba9,b69,b29,ae9,aa9,a69,a29,9e9,9a9,969,929,8e9,8a9,869,829,7e9,7a9,769,729,6e9,6a9,669,629,5e9,
		56a,52a,4ea,4aa,46a,42a,3ea,3aa,36a,32a,2ea,2aa,26a,22a,1ea,1aa,16a,12a,0ea,0aa,06a,02a,faa,f6a,f2a,eea,eaa,e6a,e2a,dea,daa,d6a,d2a,cea,caa,c6a,c2a,bea,baa,b6a,b2a,aea,aaa,a6a,a2a,9ea,9aa,96a,92a,8ea,8aa,86a,82a,7ea,7aa,76a,72a,6ea,6aa,66a,62a,5ea,5aa,
		52b,4eb,4ab,46b,42b,3eb,3ab,36b,32b,2eb,2ab,26b,22b,1eb,1ab,16b,12b,0eb,0ab,06b,02b,fab,f6b,f2b,eeb,eab,e6b,e2b,deb,dab,d6b,d2b,ceb,cab,c6b,c2b,beb,bab,b6b,b2b,aeb,aab,a6b,a2b,9eb,9ab,96b,92b,8eb,8ab,86b,82b,7eb,7ab,76b,72b,6eb,6ab,66b,62b,5eb,5ab,56b,
		4ec,4ac,46c,42c,3ec,3ac,36c,32c,2ec,2ac,26c,22c,1ec,1ac,16c,12c,0ec,0ac,06c,02c,fac,f6c,f2c,eec,eac,e6c,e2c,dec,dac,d6c,d2c,cec,cac,c6c,c2c,bec,bac,b6c,b2c,aec,aac,a6c,a2c,9ec,9ac,96c,92c,8ec,8ac,86c,82c,7ec,7ac,76c,72c,6ec,6ac,66c,62c,5ec,5ac,56c,52c,
		4ad,46d,42d,3ed,3ad,36d,32d,2ed,2ad,26d,22d,1ed,1ad,16d,12d,0ed,0ad,06d,02d,fad,f6d,f2d,eed,ead,e6d,e2d,ded,dad,d6d,d2d,ced,cad,c6d,c2d,bed,bad,b6d,b2d,aed,aad,a6d,a2d,9ed,9ad,96d,92d,8ed,8ad,86d,82d,7ed,7ad,76d,72d,6ed,6ad,66d,62d,5ed,5ad,56d,52d,4ed,
		46e,42e,3ee,3ae,36e,32e,2ee,2ae,26e,22e,1ee,1ae,16e,12e,0ee,0ae,06e,02e,fae,f6e,f2e,eee,eae,e6e,e2e,dee,dae,d6e,d2e,cee,cae,c6e,c2e,bee,bae,b6e,b2e,aee,aae,a6e,a2e,9ee,9ae,96e,92e,8ee,8ae,86e,82e,7ee,7ae,76e,72e,6ee,6ae,66e,62e,5ee,5ae,56e,52e,4ee,4ae,
		42f,3ef,3af,36f,32f,2ef,2af,26f,22f,1ef,1af,16f,12f,0ef,0af,06f,02f,faf,f6f,f2f,eef,eaf,e6f,e2f,def,daf,d6f,d2f,cef,caf,c6f,c2f,bef,baf,b6f,b2f,aef,aaf,a6f,a2f,9ef,9af,96f,92f,8ef,8af,86f,82f,7ef,7af,76f,72f,6ef,6af,66f,62f,5ef,5af,56f,52f,4ef,4af,46f,
		3f0,3b0,370,330,2f0,2b0,270,230,1f0,1b0,170,130,0f0,0b0,070,030,fb0,f70,f30,ef0,eb0,e70,e30,df0,db0,d70,d30,cf0,cb0,c70,c30,bf0,bb0,b70,b30,af0,ab0,a70,a30,9f0,9b0,970,930,8f0,8b0,870,830,7f0,7b0,770,730,6f0,6b0,670,630,5f0,5b0,570,530,4f0,4b0,470,430,
		3b1,371,331,2f1,2b1,271,231,1f1,1b1,171,131,0f1,0b1,071,031,fb1,f71,f31,ef1,eb1,e71,e31,df1,db1,d71,d31,cf1,cb1,c71,c31,bf1,bb1,b71,b31,af1,ab1,a71,a31,9f1,9b1,971,931,8f1,8b1,871,831,7f1,7b1,771,731,6f1,6b1,671,631,5f1,5b1,571,531,4f1,4b1,471,431,3f1,
		372,332,2f2,2b2,272,232,1f2,1b2,172,132,0f2,0b2,072,032,fb2,f72,f32,ef2,eb2,e72,e32,df2,db2,d72,d32,cf2,cb2,c72,c32,bf2,bb2,b72,b32,af2,ab2,a72,a32,9f2,9b2,972,932,8f2,8b2,872,832,7f2,7b2,772,732,6f2,6b2,672,632,5f2,5b2,572,532,4f2,4b2,472,432,3f2,3b2,
		333,2f3,2b3,273,233,1f3,1b3,173,133,0f3,0b3,073,033,fb3,f73,f33,ef3,eb3,e73,e33,df3,db3,d73,d33,cf3,cb3,c73,c33,bf3,bb3,b73,b33,af3,ab3,a73,a33,9f3,9b3,973,933,8f3,8b3,873,833,7f3,7b3,773,733,6f3,6b3,673,633,5f3,5b3,573,533,4f3,4b3,473,433,3f3,3b3,373,
		2f4,2b4,274,234,1f4,1b4,174,134,0f4,0b4,074,034,fb4,f74,f34,ef4,eb4,e74,e34,df4,db4,d74,d34,cf4,cb4,c74,c34,bf4,bb4,b74,b34,af4,ab4,a74,a34,9f4,9b4,974,934,8f4,8b4,874,834,7f4,7b4,774,734,6f4,6b4,674,634,5f4,5b4,574,534,4f4,4b4,474,434,3f4,3b4,374,334,
		2b5,275,235,1f5,1b5,175,135,0f5,0b5,075,035,fb5,f75,f35,ef5,eb5,e75,e35,df5,db5,d75,d35,cf5,cb5,c75,c35,bf5,bb5,b75,b35,af5,ab5,a75,a35,9f5,9b5,975,935,8f5,8b5,875,835,7f5,7b5,775,735,6f5,6b5,675,635,5f5,5b5,575,535,4f5,4b5,475,435,3f5,3b5,375,335,2f5,
		276,236,1f6,1b6,176,136,0f6,0b6,076,036,fb6,f76,f36,ef6,eb6,e76,e36,df6,db6,d76,d36,cf6,cb6,c76,c36,bf6,bb6,b76,b36,af6,ab6,a76,a36,9f6,9b6,976,936,8f6,8b6,876,836,7f6,7b6,776,736,6f6,6b6,676,636,5f6,5b6,576,536,4f6,4b6,476,436,3f6,3b6,376,336,2f6,2b6,
		237,1f7,1b7,177,137,0f7,0b7,077,037,fb7,f77,f37,ef7,eb7,e77,e37,df7,db7,d77,d37,cf7,cb7,c77,c37,bf7,bb7,b77,b37,af7,ab7,a77,a37,9f7,9b7,977,937,8f7,8b7,877,837,7f7,7b7,777,737,6f7,6b7,677,637,5f7,5b7,577,537,4f7,4b7,477,437,3f7,3b7,377,337,2f7,2b7,277,
		1f8,1b8,178,138,0f8,0b8,078,038,fb8,f78,f38,ef8,eb8,e78,e38,df8,db8,d78,d38,cf8,cb8,c78,c38,bf8,bb8,b78,b38,af8,ab8,a78,a38,9f8,9b8,978,938,8f8,8b8,878,838,7f8,7b8,778,738,6f8,6b8,678,638,5f8,5b8,578,538,4f8,4b8,478,438,3f8,3b8,378,338,2f8,2b8,278,238,
		1b9,179,139,0f9,0b9,079,039,fb9,f79,f39,ef9,eb9,e79,e39,df9,db9,d79,d39,cf9,cb9,c79,c39,bf9,bb9,b79,b39,af9,ab9,a79,a39,9f9,9b9,979,939,8f9,8b9,879,839,7f9,7b9,779,739,6f9,6b9,679,639,5f9,5b9,579,539,4f9,4b9,479,439,3f9,3b9,379,339,2f9,2b9,279,239,1f9,
		17a,13a,0fa,0ba,07a,03a,fba,f7a,f3a,efa,eba,e7a,e3a,dfa,dba,d7a,d3a,cfa,cba,c7a,c3a,bfa,bba,b7a,b3a,afa,aba,a7a,a3a,9fa,9ba,97a,93a,8fa,8ba,87a,83a,7fa,7ba,77a,73a,6fa,6ba,67a,63a,5fa,5ba,57a,53a,4fa,4ba,47a,43a,3fa,3ba,37a,33a,2fa,2ba,27a,23a,1fa,1ba,
		13b,0fb,0bb,07b,03b,fbb,f7b,f3b,efb,ebb,e7b,e3b,dfb,dbb,d7b,d3b,cfb,cbb,c7b,c3b,bfb,bbb,b7b,b3b,afb,abb,a7b,a3b,9fb,9bb,97b,93b,8fb,8bb,87b,83b,7fb,7bb,77b,73b,6fb,6bb,67b,63b,5fb,5bb,57b,53b,4fb,4bb,47b,43b,3fb,3bb,37b,33b,2fb,2bb,27b,23b,1fb,1bb,17b,
		0fc,0bc,07c,03c,fbc,f7c,f3c,efc,ebc,e7c,e3c,dfc,dbc,d7c,d3c,cfc,cbc,c7c,c3c,bfc,bbc,b7c,b3c,afc,abc,a7c,a3c,9fc,9bc,97c,93c,8fc,8bc,87c,83c,7fc,7bc,77c,73c,6fc,6bc,67c,63c,5fc,5bc,57c,53c,4fc,4bc,47c,43c,3fc,3bc,37c,33c,2fc,2bc,27c,23c,1fc,1bc,17c,13c,
		0bd,07d,03d,fbd,f7d,f3d,efd,ebd,e7d,e3d,dfd,dbd,d7d,d3d,cfd,cbd,c7d,c3d,bfd,bbd,b7d,b3d,afd,abd,a7d,a3d,9fd,9bd,97d,93d,8fd,8bd,87d,83d,7fd,7bd,77d,73d,6fd,6bd,67d,63d,5fd,5bd,57d,53d,4fd,4bd,47d,43d,3fd,3bd,37d,33d,2fd,2bd,27d,23d,1fd,1bd,17d,13d,0fd,
		07e,03e,fbe,f7e,f3e,efe,ebe,e7e,e3e,dfe,dbe,d7e,d3e,cfe,cbe,c7e,c3e,bfe,bbe,b7e,b3e,afe,abe,a7e,a3e,9fe,9be,97e,93e,8fe,8be,87e,83e,7fe,7be,77e,73e,6fe,6be,67e,63e,5fe,5be,57e,53e,4fe,4be,47e,43e,3fe,3be,37e,33e,2fe,2be,27e,23e,1fe,1be,17e,13e,0fe,0be,
		03f,fbf,f7f,f3f,eff,ebf,e7f,e3f,dff,dbf,d7f,d3f,cff,cbf,c7f,c3f,bff,bbf,b7f,b3f,aff,abf,a7f,a3f,9ff,9bf,97f,93f,8ff,8bf,87f,83f,7ff,7bf,77f,73f,6ff,6bf,67f,63f,5ff,5bf,57f,53f,4ff,4bf,47f,43f,3ff,3bf,37f,33f,2ff,2bf,27f,23f,1ff,1bf,17f,13f,0ff,0bf,07f,

	If we subtract each row's common (mod 64) low p-part from all terms of the row
	as we do in the implementation to reduce the number of index offsets needing to be stored,
	we decrement 64 horizontally and 64 vertically, all (mod 4032), and where we elide the trailing 0 on all indices:

		00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04 + p00
		f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00 + p01
		f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8 + p02
		f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4 + p03
		ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0 + p04
		e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec + p05
		e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8 + p06
		e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4 + p07
		dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0 + p08
		d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc + p09
		d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8 + p0a
		d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4 + p0b
		cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0 + p0c
		c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc + p0d
		c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8 + p0e
		c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4 + p0f
		bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0 + p10
		b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc + p11
		b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8 + p12
		b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4 + p13
		ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0 + p14
		a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac + p15
		a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8 + p16
		a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4 + p17
		9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0 + p18
		98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c + p19
		94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98 + p1a
		90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94 + p1b
		8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90 + p1c
		88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c + p1d
		84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88 + p1e
		80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84 + p1f
		7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80 + p20
		78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c + p21
		74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78 + p22
		70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74 + p23
		6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70 + p24
		68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c + p25
		64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68 + p26
		60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64 + p27
		5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60 + p28
		58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c + p29
		54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58 + p2a
		50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54 + p2b
		4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50 + p2c
		48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c + p2d
		44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48 + p2e
		40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44 + p2f
		3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40 + p30
		38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c + p31
		34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38 + p32
		30,2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34 + p33
		2c,28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30 + p34
		28,24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c + p35
		24,20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28 + p36
		20,1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24 + p37
		1c,18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20 + p38
		18,14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c + p39
		14,10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18 + p3a
		10,0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14 + p3b
		0c,08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10 + p3c
		08,04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c + p3d
		04,00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08 + p3e
		00,f8,f4,f0,ec,e8,e4,e0,dc,d8,d4,d0,cc,c8,c4,c0,bc,b8,b4,b0,ac,a8,a4,a0,9c,98,94,90,8c,88,84,80,7c,78,74,70,6c,68,64,60,5c,58,54,50,4c,48,44,40,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04 + p3f

	Further dividing the above by 4 to obtain thing in terms of indices into our phi[] (multiples of p40) array:

		00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01 + p00
		3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00 + p01
		3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e + p02
		3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d + p03
		3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c + p04
		3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b + p05
		39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a + p06
		38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39 + p07
		37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38 + p08
		36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37 + p09
		35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36 + p0a
		34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35 + p0b
		33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34 + p0c
		32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33 + p0d
		31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32 + p0e
		30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31 + p0f
		2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30 + p10
		2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f + p11
		2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e + p12
		2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d + p13
		2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c + p14
		2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b + p15
		29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a + p16
		28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29 + p17
		27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28 + p18
		26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27 + p19
		25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26 + p1a
		24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25 + p1b
		23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24 + p1c
		22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23 + p1d
		21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22 + p1e
		20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21 + p1f
		1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20 + p20
		1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f + p21
		1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e + p22
		1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d + p23
		1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c + p24
		1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b + p25
		19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a + p26
		18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19 + p27
		17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18 + p28
		16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17 + p29
		15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16 + p2a
		14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15 + p2b
		13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14 + p2c
		12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13 + p2d
		11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12 + p2e
		10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11 + p2f
		0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10 + p30
		0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f + p31
		0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e + p32
		0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d + p33
		0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c + p34
		0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b + p35
		09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a + p36
		08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09 + p37
		07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08 + p38
		06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07 + p39
		05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06 + p3a
		04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05 + p3b
		03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04 + p3c
		02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03 + p3d
		01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02 + p3e
		00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01 + p3f

	In order to make the result amenable to loop-based execution, we need to encode the indices
	to the left & right of the + in easily-computable-index fashion. This scheme has 2 ingredients:
	[Use i as loop idx in comments, although use 'l' in actual code]

	[0] RHS is simply i = 0,..,3f, i.e.plo[i]

	[1] Note that the (/= p40, where 40 = hex) row data are simply cicrular-shift perms of the basic (row 0) set.

	For row i (counting from 0 to 63), the leftward-shift count needing to be applied = i.
	*/
		tptr = t;
		for(k = 0; k < 64; ++k)
		{
			iptr = dif_p40_cperms + k;
			for(l = 0; l < ODD_RADIX; l++)
			{
				io_offsets[l] = phi[iptr[l]];
			}
			RADIX_63_DIF(
				a+j1+plo[k], io_offsets, RE_IM_STRIDE,
				(double *)tptr, toff, 1
			);
			tptr += ODD_RADIX;
		}
	//...and now do 63 radix-64 transforms:
	/*
		O-perm in hex form - generated this with a view toward Impl 2-table-style streamlined indexing
		to match other large-radix DFTs in the compact-obj-code framework.
		The required output permutation is, grouping things into blocks of 64, one for each DFT output set [in hex]:

		row:						o-Indices:
		00	000,001,002,003,004,005,006,007,008,009,00a,00b,00c,00d,00e,00f,010,011,012,013,014,015,016,017,018,019,01a,01b,01c,01d,01e,01f,020,021,022,023,024,025,026,027,028,029,02a,02b,02c,02d,02e,02f,030,031,032,033,034,035,036,037,038,039,03a,03b,03c,03d,03e,03f,
		01	095,094,097,096,093,092,090,091,09d,09c,09f,09e,09b,09a,098,099,08d,08c,08f,08e,08b,08a,088,089,083,082,080,081,087,086,084,085,0b5,0b4,0b7,0b6,0b3,0b2,0b0,0b1,0bd,0bc,0bf,0be,0bb,0ba,0b8,0b9,0ad,0ac,0af,0ae,0ab,0aa,0a8,0a9,0a3,0a2,0a0,0a1,0a7,0a6,0a4,0a5,
		02	06a,06b,069,068,06e,06f,06d,06c,066,067,065,064,061,060,063,062,07a,07b,079,078,07e,07f,07d,07c,076,077,075,074,071,070,073,072,05a,05b,059,058,05e,05f,05d,05c,056,057,055,054,051,050,053,052,046,047,045,044,041,040,043,042,04e,04f,04d,04c,049,048,04b,04a,
		03	207,206,204,205,200,201,202,203,20f,20e,20c,20d,208,209,20a,20b,217,216,214,215,210,211,212,213,21f,21e,21c,21d,218,219,21a,21b,227,226,224,225,220,221,222,223,22f,22e,22c,22d,228,229,22a,22b,237,236,234,235,230,231,232,233,23f,23e,23c,23d,238,239,23a,23b,
		04	1f1,1f0,1f3,1f2,1f5,1f4,1f7,1f6,1f9,1f8,1fb,1fa,1fd,1fc,1ff,1fe,1e9,1e8,1eb,1ea,1ed,1ec,1ef,1ee,1e5,1e4,1e7,1e6,1e3,1e2,1e0,1e1,1c9,1c8,1cb,1ca,1cd,1cc,1cf,1ce,1c5,1c4,1c7,1c6,1c3,1c2,1c0,1c1,1d9,1d8,1db,1da,1dd,1dc,1df,1de,1d5,1d4,1d7,1d6,1d3,1d2,1d0,1d1,
		05	19c,19d,19e,19f,19a,19b,199,198,192,193,191,190,196,197,195,194,182,183,181,180,186,187,185,184,18a,18b,189,188,18e,18f,18d,18c,1bc,1bd,1be,1bf,1ba,1bb,1b9,1b8,1b2,1b3,1b1,1b0,1b6,1b7,1b5,1b4,1a2,1a3,1a1,1a0,1a6,1a7,1a5,1a4,1aa,1ab,1a9,1a8,1ae,1af,1ad,1ac,
		06	163,162,160,161,167,166,164,165,16b,16a,168,169,16f,16e,16c,16d,173,172,170,171,177,176,174,175,17b,17a,178,179,17f,17e,17c,17d,153,152,150,151,157,156,154,155,15b,15a,158,159,15f,15e,15c,15d,14b,14a,148,149,14f,14e,14c,14d,147,146,144,145,140,141,142,143,
		07	10e,10f,10d,10c,109,108,10b,10a,101,100,103,102,105,104,107,106,11e,11f,11d,11c,119,118,11b,11a,111,110,113,112,115,114,117,116,12e,12f,12d,12c,129,128,12b,12a,121,120,123,122,125,124,127,126,13e,13f,13d,13c,139,138,13b,13a,131,130,133,132,135,134,137,136,
		08	0f8,0f9,0fa,0fb,0fc,0fd,0fe,0ff,0f4,0f5,0f6,0f7,0f2,0f3,0f1,0f0,0e4,0e5,0e6,0e7,0e2,0e3,0e1,0e0,0ec,0ed,0ee,0ef,0ea,0eb,0e9,0e8,0c4,0c5,0c6,0c7,0c2,0c3,0c1,0c0,0cc,0cd,0ce,0cf,0ca,0cb,0c9,0c8,0d4,0d5,0d6,0d7,0d2,0d3,0d1,0d0,0dc,0dd,0de,0df,0da,0db,0d9,0d8,
		09	f9f,f9e,f9c,f9d,f98,f99,f9a,f9b,f90,f91,f92,f93,f94,f95,f96,f97,f80,f81,f82,f83,f84,f85,f86,f87,f88,f89,f8a,f8b,f8c,f8d,f8e,f8f,fbf,fbe,fbc,fbd,fb8,fb9,fba,fbb,fb0,fb1,fb2,fb3,fb4,fb5,fb6,fb7,fa0,fa1,fa2,fa3,fa4,fa5,fa6,fa7,fa8,fa9,faa,fab,fac,fad,fae,faf,
		0a	f65,f64,f67,f66,f63,f62,f60,f61,f6d,f6c,f6f,f6e,f6b,f6a,f68,f69,f75,f74,f77,f76,f73,f72,f70,f71,f7d,f7c,f7f,f7e,f7b,f7a,f78,f79,f55,f54,f57,f56,f53,f52,f50,f51,f5d,f5c,f5f,f5e,f5b,f5a,f58,f59,f4d,f4c,f4f,f4e,f4b,f4a,f48,f49,f43,f42,f40,f41,f47,f46,f44,f45,
		0b	f0a,f0b,f09,f08,f0e,f0f,f0d,f0c,f06,f07,f05,f04,f01,f00,f03,f02,f1a,f1b,f19,f18,f1e,f1f,f1d,f1c,f16,f17,f15,f14,f11,f10,f13,f12,f2a,f2b,f29,f28,f2e,f2f,f2d,f2c,f26,f27,f25,f24,f21,f20,f23,f22,f3a,f3b,f39,f38,f3e,f3f,f3d,f3c,f36,f37,f35,f34,f31,f30,f33,f32,
		0c	efb,efa,ef8,ef9,eff,efe,efc,efd,ef7,ef6,ef4,ef5,ef0,ef1,ef2,ef3,ee7,ee6,ee4,ee5,ee0,ee1,ee2,ee3,eef,eee,eec,eed,ee8,ee9,eea,eeb,ec7,ec6,ec4,ec5,ec0,ec1,ec2,ec3,ecf,ece,ecc,ecd,ec8,ec9,eca,ecb,ed7,ed6,ed4,ed5,ed0,ed1,ed2,ed3,edf,ede,edc,edd,ed8,ed9,eda,edb,
		0d	e91,e90,e93,e92,e95,e94,e97,e96,e99,e98,e9b,e9a,e9d,e9c,e9f,e9e,e89,e88,e8b,e8a,e8d,e8c,e8f,e8e,e85,e84,e87,e86,e83,e82,e80,e81,eb1,eb0,eb3,eb2,eb5,eb4,eb7,eb6,eb9,eb8,ebb,eba,ebd,ebc,ebf,ebe,ea9,ea8,eab,eaa,ead,eac,eaf,eae,ea5,ea4,ea7,ea6,ea3,ea2,ea0,ea1,
		0e	e6c,e6d,e6e,e6f,e6a,e6b,e69,e68,e62,e63,e61,e60,e66,e67,e65,e64,e7c,e7d,e7e,e7f,e7a,e7b,e79,e78,e72,e73,e71,e70,e76,e77,e75,e74,e5c,e5d,e5e,e5f,e5a,e5b,e59,e58,e52,e53,e51,e50,e56,e57,e55,e54,e42,e43,e41,e40,e46,e47,e45,e44,e4a,e4b,e49,e48,e4e,e4f,e4d,e4c,
		0f	e03,e02,e00,e01,e07,e06,e04,e05,e0b,e0a,e08,e09,e0f,e0e,e0c,e0d,e13,e12,e10,e11,e17,e16,e14,e15,e1b,e1a,e18,e19,e1f,e1e,e1c,e1d,e23,e22,e20,e21,e27,e26,e24,e25,e2b,e2a,e28,e29,e2f,e2e,e2c,e2d,e33,e32,e30,e31,e37,e36,e34,e35,e3b,e3a,e38,e39,e3f,e3e,e3c,e3d,
		10	df6,df7,df5,df4,df1,df0,df3,df2,dfe,dff,dfd,dfc,df9,df8,dfb,dfa,dee,def,ded,dec,de9,de8,deb,dea,de1,de0,de3,de2,de5,de4,de7,de6,dce,dcf,dcd,dcc,dc9,dc8,dcb,dca,dc1,dc0,dc3,dc2,dc5,dc4,dc7,dc6,dde,ddf,ddd,ddc,dd9,dd8,ddb,dda,dd1,dd0,dd3,dd2,dd5,dd4,dd7,dd6,
		11	d98,d99,d9a,d9b,d9c,d9d,d9e,d9f,d94,d95,d96,d97,d92,d93,d91,d90,d84,d85,d86,d87,d82,d83,d81,d80,d8c,d8d,d8e,d8f,d8a,d8b,d89,d88,db8,db9,dba,dbb,dbc,dbd,dbe,dbf,db4,db5,db6,db7,db2,db3,db1,db0,da4,da5,da6,da7,da2,da3,da1,da0,dac,dad,dae,daf,daa,dab,da9,da8,
		12	d6f,d6e,d6c,d6d,d68,d69,d6a,d6b,d60,d61,d62,d63,d64,d65,d66,d67,d7f,d7e,d7c,d7d,d78,d79,d7a,d7b,d70,d71,d72,d73,d74,d75,d76,d77,d5f,d5e,d5c,d5d,d58,d59,d5a,d5b,d50,d51,d52,d53,d54,d55,d56,d57,d40,d41,d42,d43,d44,d45,d46,d47,d48,d49,d4a,d4b,d4c,d4d,d4e,d4f,
		13	d05,d04,d07,d06,d03,d02,d00,d01,d0d,d0c,d0f,d0e,d0b,d0a,d08,d09,d15,d14,d17,d16,d13,d12,d10,d11,d1d,d1c,d1f,d1e,d1b,d1a,d18,d19,d25,d24,d27,d26,d23,d22,d20,d21,d2d,d2c,d2f,d2e,d2b,d2a,d28,d29,d35,d34,d37,d36,d33,d32,d30,d31,d3d,d3c,d3f,d3e,d3b,d3a,d38,d39,
		14	cf2,cf3,cf1,cf0,cf6,cf7,cf5,cf4,cfa,cfb,cf9,cf8,cfe,cff,cfd,cfc,cea,ceb,ce9,ce8,cee,cef,ced,cec,ce6,ce7,ce5,ce4,ce1,ce0,ce3,ce2,cca,ccb,cc9,cc8,cce,ccf,ccd,ccc,cc6,cc7,cc5,cc4,cc1,cc0,cc3,cc2,cda,cdb,cd9,cd8,cde,cdf,cdd,cdc,cd6,cd7,cd5,cd4,cd1,cd0,cd3,cd2,
		15	c9b,c9a,c98,c99,c9f,c9e,c9c,c9d,c97,c96,c94,c95,c90,c91,c92,c93,c87,c86,c84,c85,c80,c81,c82,c83,c8f,c8e,c8c,c8d,c88,c89,c8a,c8b,cbb,cba,cb8,cb9,cbf,cbe,cbc,cbd,cb7,cb6,cb4,cb5,cb0,cb1,cb2,cb3,ca7,ca6,ca4,ca5,ca0,ca1,ca2,ca3,caf,cae,cac,cad,ca8,ca9,caa,cab,
		16	c61,c60,c63,c62,c65,c64,c67,c66,c69,c68,c6b,c6a,c6d,c6c,c6f,c6e,c71,c70,c73,c72,c75,c74,c77,c76,c79,c78,c7b,c7a,c7d,c7c,c7f,c7e,c51,c50,c53,c52,c55,c54,c57,c56,c59,c58,c5b,c5a,c5d,c5c,c5f,c5e,c49,c48,c4b,c4a,c4d,c4c,c4f,c4e,c45,c44,c47,c46,c43,c42,c40,c41,
		17	c0c,c0d,c0e,c0f,c0a,c0b,c09,c08,c02,c03,c01,c00,c06,c07,c05,c04,c1c,c1d,c1e,c1f,c1a,c1b,c19,c18,c12,c13,c11,c10,c16,c17,c15,c14,c2c,c2d,c2e,c2f,c2a,c2b,c29,c28,c22,c23,c21,c20,c26,c27,c25,c24,c3c,c3d,c3e,c3f,c3a,c3b,c39,c38,c32,c33,c31,c30,c36,c37,c35,c34,
		18	bfd,bfc,bff,bfe,bfb,bfa,bf8,bf9,bf3,bf2,bf0,bf1,bf7,bf6,bf4,bf5,be3,be2,be0,be1,be7,be6,be4,be5,beb,bea,be8,be9,bef,bee,bec,bed,bc3,bc2,bc0,bc1,bc7,bc6,bc4,bc5,bcb,bca,bc8,bc9,bcf,bce,bcc,bcd,bd3,bd2,bd0,bd1,bd7,bd6,bd4,bd5,bdb,bda,bd8,bd9,bdf,bde,bdc,bdd,
		19	b96,b97,b95,b94,b91,b90,b93,b92,b9e,b9f,b9d,b9c,b99,b98,b9b,b9a,b8e,b8f,b8d,b8c,b89,b88,b8b,b8a,b81,b80,b83,b82,b85,b84,b87,b86,bb6,bb7,bb5,bb4,bb1,bb0,bb3,bb2,bbe,bbf,bbd,bbc,bb9,bb8,bbb,bba,bae,baf,bad,bac,ba9,ba8,bab,baa,ba1,ba0,ba3,ba2,ba5,ba4,ba7,ba6,
		1a	b68,b69,b6a,b6b,b6c,b6d,b6e,b6f,b64,b65,b66,b67,b62,b63,b61,b60,b78,b79,b7a,b7b,b7c,b7d,b7e,b7f,b74,b75,b76,b77,b72,b73,b71,b70,b58,b59,b5a,b5b,b5c,b5d,b5e,b5f,b54,b55,b56,b57,b52,b53,b51,b50,b44,b45,b46,b47,b42,b43,b41,b40,b4c,b4d,b4e,b4f,b4a,b4b,b49,b48,
		1b	b0f,b0e,b0c,b0d,b08,b09,b0a,b0b,b00,b01,b02,b03,b04,b05,b06,b07,b1f,b1e,b1c,b1d,b18,b19,b1a,b1b,b10,b11,b12,b13,b14,b15,b16,b17,b2f,b2e,b2c,b2d,b28,b29,b2a,b2b,b20,b21,b22,b23,b24,b25,b26,b27,b3f,b3e,b3c,b3d,b38,b39,b3a,b3b,b30,b31,b32,b33,b34,b35,b36,b37,
		1c	af9,af8,afb,afa,afd,afc,aff,afe,af5,af4,af7,af6,af3,af2,af0,af1,ae5,ae4,ae7,ae6,ae3,ae2,ae0,ae1,aed,aec,aef,aee,aeb,aea,ae8,ae9,ac5,ac4,ac7,ac6,ac3,ac2,ac0,ac1,acd,acc,acf,ace,acb,aca,ac8,ac9,ad5,ad4,ad7,ad6,ad3,ad2,ad0,ad1,add,adc,adf,ade,adb,ada,ad8,ad9,
		1d	a92,a93,a91,a90,a96,a97,a95,a94,a9a,a9b,a99,a98,a9e,a9f,a9d,a9c,a8a,a8b,a89,a88,a8e,a8f,a8d,a8c,a86,a87,a85,a84,a81,a80,a83,a82,ab2,ab3,ab1,ab0,ab6,ab7,ab5,ab4,aba,abb,ab9,ab8,abe,abf,abd,abc,aaa,aab,aa9,aa8,aae,aaf,aad,aac,aa6,aa7,aa5,aa4,aa1,aa0,aa3,aa2,
		1e	a6b,a6a,a68,a69,a6f,a6e,a6c,a6d,a67,a66,a64,a65,a60,a61,a62,a63,a7b,a7a,a78,a79,a7f,a7e,a7c,a7d,a77,a76,a74,a75,a70,a71,a72,a73,a5b,a5a,a58,a59,a5f,a5e,a5c,a5d,a57,a56,a54,a55,a50,a51,a52,a53,a47,a46,a44,a45,a40,a41,a42,a43,a4f,a4e,a4c,a4d,a48,a49,a4a,a4b,
		1f	a01,a00,a03,a02,a05,a04,a07,a06,a09,a08,a0b,a0a,a0d,a0c,a0f,a0e,a11,a10,a13,a12,a15,a14,a17,a16,a19,a18,a1b,a1a,a1d,a1c,a1f,a1e,a21,a20,a23,a22,a25,a24,a27,a26,a29,a28,a2b,a2a,a2d,a2c,a2f,a2e,a31,a30,a33,a32,a35,a34,a37,a36,a39,a38,a3b,a3a,a3d,a3c,a3f,a3e,
		20	9f4,9f5,9f6,9f7,9f2,9f3,9f1,9f0,9fc,9fd,9fe,9ff,9fa,9fb,9f9,9f8,9ec,9ed,9ee,9ef,9ea,9eb,9e9,9e8,9e2,9e3,9e1,9e0,9e6,9e7,9e5,9e4,9cc,9cd,9ce,9cf,9ca,9cb,9c9,9c8,9c2,9c3,9c1,9c0,9c6,9c7,9c5,9c4,9dc,9dd,9de,9df,9da,9db,9d9,9d8,9d2,9d3,9d1,9d0,9d6,9d7,9d5,9d4,
		21	99d,99c,99f,99e,99b,99a,998,999,993,992,990,991,997,996,994,995,983,982,980,981,987,986,984,985,98b,98a,988,989,98f,98e,98c,98d,9bd,9bc,9bf,9be,9bb,9ba,9b8,9b9,9b3,9b2,9b0,9b1,9b7,9b6,9b4,9b5,9a3,9a2,9a0,9a1,9a7,9a6,9a4,9a5,9ab,9aa,9a8,9a9,9af,9ae,9ac,9ad,
		22	966,967,965,964,961,960,963,962,96e,96f,96d,96c,969,968,96b,96a,976,977,975,974,971,970,973,972,97e,97f,97d,97c,979,978,97b,97a,956,957,955,954,951,950,953,952,95e,95f,95d,95c,959,958,95b,95a,94e,94f,94d,94c,949,948,94b,94a,941,940,943,942,945,944,947,946,
		23	908,909,90a,90b,90c,90d,90e,90f,904,905,906,907,902,903,901,900,918,919,91a,91b,91c,91d,91e,91f,914,915,916,917,912,913,911,910,928,929,92a,92b,92c,92d,92e,92f,924,925,926,927,922,923,921,920,938,939,93a,93b,93c,93d,93e,93f,934,935,936,937,932,933,931,930,
		24	8f7,8f6,8f4,8f5,8f0,8f1,8f2,8f3,8ff,8fe,8fc,8fd,8f8,8f9,8fa,8fb,8ef,8ee,8ec,8ed,8e8,8e9,8ea,8eb,8e0,8e1,8e2,8e3,8e4,8e5,8e6,8e7,8cf,8ce,8cc,8cd,8c8,8c9,8ca,8cb,8c0,8c1,8c2,8c3,8c4,8c5,8c6,8c7,8df,8de,8dc,8dd,8d8,8d9,8da,8db,8d0,8d1,8d2,8d3,8d4,8d5,8d6,8d7,
		25	899,898,89b,89a,89d,89c,89f,89e,895,894,897,896,893,892,890,891,885,884,887,886,883,882,880,881,88d,88c,88f,88e,88b,88a,888,889,8b9,8b8,8bb,8ba,8bd,8bc,8bf,8be,8b5,8b4,8b7,8b6,8b3,8b2,8b0,8b1,8a5,8a4,8a7,8a6,8a3,8a2,8a0,8a1,8ad,8ac,8af,8ae,8ab,8aa,8a8,8a9,
		26	862,863,861,860,866,867,865,864,86a,86b,869,868,86e,86f,86d,86c,872,873,871,870,876,877,875,874,87a,87b,879,878,87e,87f,87d,87c,852,853,851,850,856,857,855,854,85a,85b,859,858,85e,85f,85d,85c,84a,84b,849,848,84e,84f,84d,84c,846,847,845,844,841,840,843,842,
		27	80b,80a,808,809,80f,80e,80c,80d,807,806,804,805,800,801,802,803,81b,81a,818,819,81f,81e,81c,81d,817,816,814,815,810,811,812,813,82b,82a,828,829,82f,82e,82c,82d,827,826,824,825,820,821,822,823,83b,83a,838,839,83f,83e,83c,83d,837,836,834,835,830,831,832,833,
		28	7fe,7ff,7fd,7fc,7f9,7f8,7fb,7fa,7f1,7f0,7f3,7f2,7f5,7f4,7f7,7f6,7e1,7e0,7e3,7e2,7e5,7e4,7e7,7e6,7e9,7e8,7eb,7ea,7ed,7ec,7ef,7ee,7c1,7c0,7c3,7c2,7c5,7c4,7c7,7c6,7c9,7c8,7cb,7ca,7cd,7cc,7cf,7ce,7d1,7d0,7d3,7d2,7d5,7d4,7d7,7d6,7d9,7d8,7db,7da,7dd,7dc,7df,7de,
		29	794,795,796,797,792,793,791,790,79c,79d,79e,79f,79a,79b,799,798,78c,78d,78e,78f,78a,78b,789,788,782,783,781,780,786,787,785,784,7b4,7b5,7b6,7b7,7b2,7b3,7b1,7b0,7bc,7bd,7be,7bf,7ba,7bb,7b9,7b8,7ac,7ad,7ae,7af,7aa,7ab,7a9,7a8,7a2,7a3,7a1,7a0,7a6,7a7,7a5,7a4,
		2a	76d,76c,76f,76e,76b,76a,768,769,763,762,760,761,767,766,764,765,77d,77c,77f,77e,77b,77a,778,779,773,772,770,771,777,776,774,775,75d,75c,75f,75e,75b,75a,758,759,753,752,750,751,757,756,754,755,743,742,740,741,747,746,744,745,74b,74a,748,749,74f,74e,74c,74d,
		2b	706,707,705,704,701,700,703,702,70e,70f,70d,70c,709,708,70b,70a,716,717,715,714,711,710,713,712,71e,71f,71d,71c,719,718,71b,71a,726,727,725,724,721,720,723,722,72e,72f,72d,72c,729,728,72b,72a,736,737,735,734,731,730,733,732,73e,73f,73d,73c,739,738,73b,73a,
		2c	6f0,6f1,6f2,6f3,6f4,6f5,6f6,6f7,6f8,6f9,6fa,6fb,6fc,6fd,6fe,6ff,6e8,6e9,6ea,6eb,6ec,6ed,6ee,6ef,6e4,6e5,6e6,6e7,6e2,6e3,6e1,6e0,6c8,6c9,6ca,6cb,6cc,6cd,6ce,6cf,6c4,6c5,6c6,6c7,6c2,6c3,6c1,6c0,6d8,6d9,6da,6db,6dc,6dd,6de,6df,6d4,6d5,6d6,6d7,6d2,6d3,6d1,6d0,
		2d	697,696,694,695,690,691,692,693,69f,69e,69c,69d,698,699,69a,69b,68f,68e,68c,68d,688,689,68a,68b,680,681,682,683,684,685,686,687,6b7,6b6,6b4,6b5,6b0,6b1,6b2,6b3,6bf,6be,6bc,6bd,6b8,6b9,6ba,6bb,6af,6ae,6ac,6ad,6a8,6a9,6aa,6ab,6a0,6a1,6a2,6a3,6a4,6a5,6a6,6a7,
		2e	669,668,66b,66a,66d,66c,66f,66e,665,664,667,666,663,662,660,661,679,678,67b,67a,67d,67c,67f,67e,675,674,677,676,673,672,670,671,659,658,65b,65a,65d,65c,65f,65e,655,654,657,656,653,652,650,651,645,644,647,646,643,642,640,641,64d,64c,64f,64e,64b,64a,648,649,
		2f	602,603,601,600,606,607,605,604,60a,60b,609,608,60e,60f,60d,60c,612,613,611,610,616,617,615,614,61a,61b,619,618,61e,61f,61d,61c,622,623,621,620,626,627,625,624,62a,62b,629,628,62e,62f,62d,62c,632,633,631,630,636,637,635,634,63a,63b,639,638,63e,63f,63d,63c,
		30	5f3,5f2,5f0,5f1,5f7,5f6,5f4,5f5,5fb,5fa,5f8,5f9,5ff,5fe,5fc,5fd,5eb,5ea,5e8,5e9,5ef,5ee,5ec,5ed,5e7,5e6,5e4,5e5,5e0,5e1,5e2,5e3,5cb,5ca,5c8,5c9,5cf,5ce,5cc,5cd,5c7,5c6,5c4,5c5,5c0,5c1,5c2,5c3,5db,5da,5d8,5d9,5df,5de,5dc,5dd,5d7,5d6,5d4,5d5,5d0,5d1,5d2,5d3,
		31	59e,59f,59d,59c,599,598,59b,59a,591,590,593,592,595,594,597,596,581,580,583,582,585,584,587,586,589,588,58b,58a,58d,58c,58f,58e,5be,5bf,5bd,5bc,5b9,5b8,5bb,5ba,5b1,5b0,5b3,5b2,5b5,5b4,5b7,5b6,5a1,5a0,5a3,5a2,5a5,5a4,5a7,5a6,5a9,5a8,5ab,5aa,5ad,5ac,5af,5ae,
		32	564,565,566,567,562,563,561,560,56c,56d,56e,56f,56a,56b,569,568,574,575,576,577,572,573,571,570,57c,57d,57e,57f,57a,57b,579,578,554,555,556,557,552,553,551,550,55c,55d,55e,55f,55a,55b,559,558,54c,54d,54e,54f,54a,54b,549,548,542,543,541,540,546,547,545,544,
		33	50d,50c,50f,50e,50b,50a,508,509,503,502,500,501,507,506,504,505,51d,51c,51f,51e,51b,51a,518,519,513,512,510,511,517,516,514,515,52d,52c,52f,52e,52b,52a,528,529,523,522,520,521,527,526,524,525,53d,53c,53f,53e,53b,53a,538,539,533,532,530,531,537,536,534,535,
		34	4fa,4fb,4f9,4f8,4fe,4ff,4fd,4fc,4f6,4f7,4f5,4f4,4f1,4f0,4f3,4f2,4e6,4e7,4e5,4e4,4e1,4e0,4e3,4e2,4ee,4ef,4ed,4ec,4e9,4e8,4eb,4ea,4c6,4c7,4c5,4c4,4c1,4c0,4c3,4c2,4ce,4cf,4cd,4cc,4c9,4c8,4cb,4ca,4d6,4d7,4d5,4d4,4d1,4d0,4d3,4d2,4de,4df,4dd,4dc,4d9,4d8,4db,4da,
		35	490,491,492,493,494,495,496,497,498,499,49a,49b,49c,49d,49e,49f,488,489,48a,48b,48c,48d,48e,48f,484,485,486,487,482,483,481,480,4b0,4b1,4b2,4b3,4b4,4b5,4b6,4b7,4b8,4b9,4ba,4bb,4bc,4bd,4be,4bf,4a8,4a9,4aa,4ab,4ac,4ad,4ae,4af,4a4,4a5,4a6,4a7,4a2,4a3,4a1,4a0,
		36	467,466,464,465,460,461,462,463,46f,46e,46c,46d,468,469,46a,46b,477,476,474,475,470,471,472,473,47f,47e,47c,47d,478,479,47a,47b,457,456,454,455,450,451,452,453,45f,45e,45c,45d,458,459,45a,45b,44f,44e,44c,44d,448,449,44a,44b,440,441,442,443,444,445,446,447,
		37	409,408,40b,40a,40d,40c,40f,40e,405,404,407,406,403,402,400,401,419,418,41b,41a,41d,41c,41f,41e,415,414,417,416,413,412,410,411,429,428,42b,42a,42d,42c,42f,42e,425,424,427,426,423,422,420,421,439,438,43b,43a,43d,43c,43f,43e,435,434,437,436,433,432,430,431,
		38	3fc,3fd,3fe,3ff,3fa,3fb,3f9,3f8,3f2,3f3,3f1,3f0,3f6,3f7,3f5,3f4,3e2,3e3,3e1,3e0,3e6,3e7,3e5,3e4,3ea,3eb,3e9,3e8,3ee,3ef,3ed,3ec,3c2,3c3,3c1,3c0,3c6,3c7,3c5,3c4,3ca,3cb,3c9,3c8,3ce,3cf,3cd,3cc,3d2,3d3,3d1,3d0,3d6,3d7,3d5,3d4,3da,3db,3d9,3d8,3de,3df,3dd,3dc,
		39	393,392,390,391,397,396,394,395,39b,39a,398,399,39f,39e,39c,39d,38b,38a,388,389,38f,38e,38c,38d,387,386,384,385,380,381,382,383,3b3,3b2,3b0,3b1,3b7,3b6,3b4,3b5,3bb,3ba,3b8,3b9,3bf,3be,3bc,3bd,3ab,3aa,3a8,3a9,3af,3ae,3ac,3ad,3a7,3a6,3a4,3a5,3a0,3a1,3a2,3a3,
		3a	36e,36f,36d,36c,369,368,36b,36a,361,360,363,362,365,364,367,366,37e,37f,37d,37c,379,378,37b,37a,371,370,373,372,375,374,377,376,35e,35f,35d,35c,359,358,35b,35a,351,350,353,352,355,354,357,356,341,340,343,342,345,344,347,346,349,348,34b,34a,34d,34c,34f,34e,
		3b	304,305,306,307,302,303,301,300,30c,30d,30e,30f,30a,30b,309,308,314,315,316,317,312,313,311,310,31c,31d,31e,31f,31a,31b,319,318,324,325,326,327,322,323,321,320,32c,32d,32e,32f,32a,32b,329,328,334,335,336,337,332,333,331,330,33c,33d,33e,33f,33a,33b,339,338,
		3c	2f5,2f4,2f7,2f6,2f3,2f2,2f0,2f1,2fd,2fc,2ff,2fe,2fb,2fa,2f8,2f9,2ed,2ec,2ef,2ee,2eb,2ea,2e8,2e9,2e3,2e2,2e0,2e1,2e7,2e6,2e4,2e5,2cd,2cc,2cf,2ce,2cb,2ca,2c8,2c9,2c3,2c2,2c0,2c1,2c7,2c6,2c4,2c5,2dd,2dc,2df,2de,2db,2da,2d8,2d9,2d3,2d2,2d0,2d1,2d7,2d6,2d4,2d5,
		3d	29a,29b,299,298,29e,29f,29d,29c,296,297,295,294,291,290,293,292,286,287,285,284,281,280,283,282,28e,28f,28d,28c,289,288,28b,28a,2ba,2bb,2b9,2b8,2be,2bf,2bd,2bc,2b6,2b7,2b5,2b4,2b1,2b0,2b3,2b2,2a6,2a7,2a5,2a4,2a1,2a0,2a3,2a2,2ae,2af,2ad,2ac,2a9,2a8,2ab,2aa,
		3e	260,261,262,263,264,265,266,267,268,269,26a,26b,26c,26d,26e,26f,270,271,272,273,274,275,276,277,278,279,27a,27b,27c,27d,27e,27f,250,251,252,253,254,255,256,257,258,259,25a,25b,25c,25d,25e,25f,248,249,24a,24b,24c,24d,24e,24f,244,245,246,247,242,243,241,240,

	Separating the above into (mod 64) - indices into plo[] array - and (/= 64) - indices into phi[] (multiples of p40) array -
	Note the phi index pattern is identical to radix-1008 DIF, except here phi has 64-multiples of plo, rather than 16-multiples:

		row:						o-Indices:
		00	p[00,01,02,03,04,05,06,07,08,09,0a,0b,0c,0d,0e,0f,10,11,12,13,14,15,16,17,18,19,1a,1b,1c,1d,1e,1f,20,21,22,23,24,25,26,27,28,29,2a,2b,2c,2d,2e,2f,30,31,32,33,34,35,36,37,38,39,3a,3b,3c,3d,3e,3f] + phi[00]
		01	p[15,14,17,16,13,12,10,11,1d,1c,1f,1e,1b,1a,18,19,0d,0c,0f,0e,0b,0a,08,09,03,02,00,01,07,06,04,05,35,34,37,36,33,32,30,31,3d,3c,3f,3e,3b,3a,38,39,2d,2c,2f,2e,2b,2a,28,29,23,22,20,21,27,26,24,25] + phi[02]
		02	p[2a,2b,29,28,2e,2f,2d,2c,26,27,25,24,21,20,23,22,3a,3b,39,38,3e,3f,3d,3c,36,37,35,34,31,30,33,32,1a,1b,19,18,1e,1f,1d,1c,16,17,15,14,11,10,13,12,06,07,05,04,01,00,03,02,0e,0f,0d,0c,09,08,0b,0a] + phi[01]
		03	p[07,06,04,05,00,01,02,03,0f,0e,0c,0d,08,09,0a,0b,17,16,14,15,10,11,12,13,1f,1e,1c,1d,18,19,1a,1b,27,26,24,25,20,21,22,23,2f,2e,2c,2d,28,29,2a,2b,37,36,34,35,30,31,32,33,3f,3e,3c,3d,38,39,3a,3b] + phi[08]
		04	p[31,30,33,32,35,34,37,36,39,38,3b,3a,3d,3c,3f,3e,29,28,2b,2a,2d,2c,2f,2e,25,24,27,26,23,22,20,21,09,08,0b,0a,0d,0c,0f,0e,05,04,07,06,03,02,00,01,19,18,1b,1a,1d,1c,1f,1e,15,14,17,16,13,12,10,11] + phi[07]
		05	p[1c,1d,1e,1f,1a,1b,19,18,12,13,11,10,16,17,15,14,02,03,01,00,06,07,05,04,0a,0b,09,08,0e,0f,0d,0c,3c,3d,3e,3f,3a,3b,39,38,32,33,31,30,36,37,35,34,22,23,21,20,26,27,25,24,2a,2b,29,28,2e,2f,2d,2c] + phi[06]
		06	p[23,22,20,21,27,26,24,25,2b,2a,28,29,2f,2e,2c,2d,33,32,30,31,37,36,34,35,3b,3a,38,39,3f,3e,3c,3d,13,12,10,11,17,16,14,15,1b,1a,18,19,1f,1e,1c,1d,0b,0a,08,09,0f,0e,0c,0d,07,06,04,05,00,01,02,03] + phi[05]
		07	p[0e,0f,0d,0c,09,08,0b,0a,01,00,03,02,05,04,07,06,1e,1f,1d,1c,19,18,1b,1a,11,10,13,12,15,14,17,16,2e,2f,2d,2c,29,28,2b,2a,21,20,23,22,25,24,27,26,3e,3f,3d,3c,39,38,3b,3a,31,30,33,32,35,34,37,36] + phi[04]
		08	p[38,39,3a,3b,3c,3d,3e,3f,34,35,36,37,32,33,31,30,24,25,26,27,22,23,21,20,2c,2d,2e,2f,2a,2b,29,28,04,05,06,07,02,03,01,00,0c,0d,0e,0f,0a,0b,09,08,14,15,16,17,12,13,11,10,1c,1d,1e,1f,1a,1b,19,18] + phi[03]
		09	p[1f,1e,1c,1d,18,19,1a,1b,10,11,12,13,14,15,16,17,00,01,02,03,04,05,06,07,08,09,0a,0b,0c,0d,0e,0f,3f,3e,3c,3d,38,39,3a,3b,30,31,32,33,34,35,36,37,20,21,22,23,24,25,26,27,28,29,2a,2b,2c,2d,2e,2f] + phi[3e]
		0a	p[25,24,27,26,23,22,20,21,2d,2c,2f,2e,2b,2a,28,29,35,34,37,36,33,32,30,31,3d,3c,3f,3e,3b,3a,38,39,15,14,17,16,13,12,10,11,1d,1c,1f,1e,1b,1a,18,19,0d,0c,0f,0e,0b,0a,08,09,03,02,00,01,07,06,04,05] + phi[3d]
		0b	p[0a,0b,09,08,0e,0f,0d,0c,06,07,05,04,01,00,03,02,1a,1b,19,18,1e,1f,1d,1c,16,17,15,14,11,10,13,12,2a,2b,29,28,2e,2f,2d,2c,26,27,25,24,21,20,23,22,3a,3b,39,38,3e,3f,3d,3c,36,37,35,34,31,30,33,32] + phi[3c]
		0c	p[3b,3a,38,39,3f,3e,3c,3d,37,36,34,35,30,31,32,33,27,26,24,25,20,21,22,23,2f,2e,2c,2d,28,29,2a,2b,07,06,04,05,00,01,02,03,0f,0e,0c,0d,08,09,0a,0b,17,16,14,15,10,11,12,13,1f,1e,1c,1d,18,19,1a,1b] + phi[3b]
		0d	p[11,10,13,12,15,14,17,16,19,18,1b,1a,1d,1c,1f,1e,09,08,0b,0a,0d,0c,0f,0e,05,04,07,06,03,02,00,01,31,30,33,32,35,34,37,36,39,38,3b,3a,3d,3c,3f,3e,29,28,2b,2a,2d,2c,2f,2e,25,24,27,26,23,22,20,21] + phi[3a]
		0e	p[2c,2d,2e,2f,2a,2b,29,28,22,23,21,20,26,27,25,24,3c,3d,3e,3f,3a,3b,39,38,32,33,31,30,36,37,35,34,1c,1d,1e,1f,1a,1b,19,18,12,13,11,10,16,17,15,14,02,03,01,00,06,07,05,04,0a,0b,09,08,0e,0f,0d,0c] + phi[39]
		0f	p[03,02,00,01,07,06,04,05,0b,0a,08,09,0f,0e,0c,0d,13,12,10,11,17,16,14,15,1b,1a,18,19,1f,1e,1c,1d,23,22,20,21,27,26,24,25,2b,2a,28,29,2f,2e,2c,2d,33,32,30,31,37,36,34,35,3b,3a,38,39,3f,3e,3c,3d] + phi[38]
		10	p[36,37,35,34,31,30,33,32,3e,3f,3d,3c,39,38,3b,3a,2e,2f,2d,2c,29,28,2b,2a,21,20,23,22,25,24,27,26,0e,0f,0d,0c,09,08,0b,0a,01,00,03,02,05,04,07,06,1e,1f,1d,1c,19,18,1b,1a,11,10,13,12,15,14,17,16] + phi[37]
		11	p[18,19,1a,1b,1c,1d,1e,1f,14,15,16,17,12,13,11,10,04,05,06,07,02,03,01,00,0c,0d,0e,0f,0a,0b,09,08,38,39,3a,3b,3c,3d,3e,3f,34,35,36,37,32,33,31,30,24,25,26,27,22,23,21,20,2c,2d,2e,2f,2a,2b,29,28] + phi[36]
		12	p[2f,2e,2c,2d,28,29,2a,2b,20,21,22,23,24,25,26,27,3f,3e,3c,3d,38,39,3a,3b,30,31,32,33,34,35,36,37,1f,1e,1c,1d,18,19,1a,1b,10,11,12,13,14,15,16,17,00,01,02,03,04,05,06,07,08,09,0a,0b,0c,0d,0e,0f] + phi[35]
		13	p[05,04,07,06,03,02,00,01,0d,0c,0f,0e,0b,0a,08,09,15,14,17,16,13,12,10,11,1d,1c,1f,1e,1b,1a,18,19,25,24,27,26,23,22,20,21,2d,2c,2f,2e,2b,2a,28,29,35,34,37,36,33,32,30,31,3d,3c,3f,3e,3b,3a,38,39] + phi[34]
		14	p[32,33,31,30,36,37,35,34,3a,3b,39,38,3e,3f,3d,3c,2a,2b,29,28,2e,2f,2d,2c,26,27,25,24,21,20,23,22,0a,0b,09,08,0e,0f,0d,0c,06,07,05,04,01,00,03,02,1a,1b,19,18,1e,1f,1d,1c,16,17,15,14,11,10,13,12] + phi[33]
		15	p[1b,1a,18,19,1f,1e,1c,1d,17,16,14,15,10,11,12,13,07,06,04,05,00,01,02,03,0f,0e,0c,0d,08,09,0a,0b,3b,3a,38,39,3f,3e,3c,3d,37,36,34,35,30,31,32,33,27,26,24,25,20,21,22,23,2f,2e,2c,2d,28,29,2a,2b] + phi[32]
		16	p[21,20,23,22,25,24,27,26,29,28,2b,2a,2d,2c,2f,2e,31,30,33,32,35,34,37,36,39,38,3b,3a,3d,3c,3f,3e,11,10,13,12,15,14,17,16,19,18,1b,1a,1d,1c,1f,1e,09,08,0b,0a,0d,0c,0f,0e,05,04,07,06,03,02,00,01] + phi[31]
		17	p[0c,0d,0e,0f,0a,0b,09,08,02,03,01,00,06,07,05,04,1c,1d,1e,1f,1a,1b,19,18,12,13,11,10,16,17,15,14,2c,2d,2e,2f,2a,2b,29,28,22,23,21,20,26,27,25,24,3c,3d,3e,3f,3a,3b,39,38,32,33,31,30,36,37,35,34] + phi[30]
		18	p[3d,3c,3f,3e,3b,3a,38,39,33,32,30,31,37,36,34,35,23,22,20,21,27,26,24,25,2b,2a,28,29,2f,2e,2c,2d,03,02,00,01,07,06,04,05,0b,0a,08,09,0f,0e,0c,0d,13,12,10,11,17,16,14,15,1b,1a,18,19,1f,1e,1c,1d] + phi[2f]
		19	p[16,17,15,14,11,10,13,12,1e,1f,1d,1c,19,18,1b,1a,0e,0f,0d,0c,09,08,0b,0a,01,00,03,02,05,04,07,06,36,37,35,34,31,30,33,32,3e,3f,3d,3c,39,38,3b,3a,2e,2f,2d,2c,29,28,2b,2a,21,20,23,22,25,24,27,26] + phi[2e]
		1a	p[28,29,2a,2b,2c,2d,2e,2f,24,25,26,27,22,23,21,20,38,39,3a,3b,3c,3d,3e,3f,34,35,36,37,32,33,31,30,18,19,1a,1b,1c,1d,1e,1f,14,15,16,17,12,13,11,10,04,05,06,07,02,03,01,00,0c,0d,0e,0f,0a,0b,09,08] + phi[2d]
		1b	p[0f,0e,0c,0d,08,09,0a,0b,00,01,02,03,04,05,06,07,1f,1e,1c,1d,18,19,1a,1b,10,11,12,13,14,15,16,17,2f,2e,2c,2d,28,29,2a,2b,20,21,22,23,24,25,26,27,3f,3e,3c,3d,38,39,3a,3b,30,31,32,33,34,35,36,37] + phi[2c]
		1c	p[39,38,3b,3a,3d,3c,3f,3e,35,34,37,36,33,32,30,31,25,24,27,26,23,22,20,21,2d,2c,2f,2e,2b,2a,28,29,05,04,07,06,03,02,00,01,0d,0c,0f,0e,0b,0a,08,09,15,14,17,16,13,12,10,11,1d,1c,1f,1e,1b,1a,18,19] + phi[2b]
		1d	p[12,13,11,10,16,17,15,14,1a,1b,19,18,1e,1f,1d,1c,0a,0b,09,08,0e,0f,0d,0c,06,07,05,04,01,00,03,02,32,33,31,30,36,37,35,34,3a,3b,39,38,3e,3f,3d,3c,2a,2b,29,28,2e,2f,2d,2c,26,27,25,24,21,20,23,22] + phi[2a]
		1e	p[2b,2a,28,29,2f,2e,2c,2d,27,26,24,25,20,21,22,23,3b,3a,38,39,3f,3e,3c,3d,37,36,34,35,30,31,32,33,1b,1a,18,19,1f,1e,1c,1d,17,16,14,15,10,11,12,13,07,06,04,05,00,01,02,03,0f,0e,0c,0d,08,09,0a,0b] + phi[29]
		1f	p[01,00,03,02,05,04,07,06,09,08,0b,0a,0d,0c,0f,0e,11,10,13,12,15,14,17,16,19,18,1b,1a,1d,1c,1f,1e,21,20,23,22,25,24,27,26,29,28,2b,2a,2d,2c,2f,2e,31,30,33,32,35,34,37,36,39,38,3b,3a,3d,3c,3f,3e] + phi[28]
		20	p[34,35,36,37,32,33,31,30,3c,3d,3e,3f,3a,3b,39,38,2c,2d,2e,2f,2a,2b,29,28,22,23,21,20,26,27,25,24,0c,0d,0e,0f,0a,0b,09,08,02,03,01,00,06,07,05,04,1c,1d,1e,1f,1a,1b,19,18,12,13,11,10,16,17,15,14] + phi[27]
		21	p[1d,1c,1f,1e,1b,1a,18,19,13,12,10,11,17,16,14,15,03,02,00,01,07,06,04,05,0b,0a,08,09,0f,0e,0c,0d,3d,3c,3f,3e,3b,3a,38,39,33,32,30,31,37,36,34,35,23,22,20,21,27,26,24,25,2b,2a,28,29,2f,2e,2c,2d] + phi[26]
		22	p[26,27,25,24,21,20,23,22,2e,2f,2d,2c,29,28,2b,2a,36,37,35,34,31,30,33,32,3e,3f,3d,3c,39,38,3b,3a,16,17,15,14,11,10,13,12,1e,1f,1d,1c,19,18,1b,1a,0e,0f,0d,0c,09,08,0b,0a,01,00,03,02,05,04,07,06] + phi[25]
		23	p[08,09,0a,0b,0c,0d,0e,0f,04,05,06,07,02,03,01,00,18,19,1a,1b,1c,1d,1e,1f,14,15,16,17,12,13,11,10,28,29,2a,2b,2c,2d,2e,2f,24,25,26,27,22,23,21,20,38,39,3a,3b,3c,3d,3e,3f,34,35,36,37,32,33,31,30] + phi[24]
		24	p[37,36,34,35,30,31,32,33,3f,3e,3c,3d,38,39,3a,3b,2f,2e,2c,2d,28,29,2a,2b,20,21,22,23,24,25,26,27,0f,0e,0c,0d,08,09,0a,0b,00,01,02,03,04,05,06,07,1f,1e,1c,1d,18,19,1a,1b,10,11,12,13,14,15,16,17] + phi[23]
		25	p[19,18,1b,1a,1d,1c,1f,1e,15,14,17,16,13,12,10,11,05,04,07,06,03,02,00,01,0d,0c,0f,0e,0b,0a,08,09,39,38,3b,3a,3d,3c,3f,3e,35,34,37,36,33,32,30,31,25,24,27,26,23,22,20,21,2d,2c,2f,2e,2b,2a,28,29] + phi[22]
		26	p[22,23,21,20,26,27,25,24,2a,2b,29,28,2e,2f,2d,2c,32,33,31,30,36,37,35,34,3a,3b,39,38,3e,3f,3d,3c,12,13,11,10,16,17,15,14,1a,1b,19,18,1e,1f,1d,1c,0a,0b,09,08,0e,0f,0d,0c,06,07,05,04,01,00,03,02] + phi[21]
		27	p[0b,0a,08,09,0f,0e,0c,0d,07,06,04,05,00,01,02,03,1b,1a,18,19,1f,1e,1c,1d,17,16,14,15,10,11,12,13,2b,2a,28,29,2f,2e,2c,2d,27,26,24,25,20,21,22,23,3b,3a,38,39,3f,3e,3c,3d,37,36,34,35,30,31,32,33] + phi[20]
		28	p[3e,3f,3d,3c,39,38,3b,3a,31,30,33,32,35,34,37,36,21,20,23,22,25,24,27,26,29,28,2b,2a,2d,2c,2f,2e,01,00,03,02,05,04,07,06,09,08,0b,0a,0d,0c,0f,0e,11,10,13,12,15,14,17,16,19,18,1b,1a,1d,1c,1f,1e] + phi[1f]
		29	p[14,15,16,17,12,13,11,10,1c,1d,1e,1f,1a,1b,19,18,0c,0d,0e,0f,0a,0b,09,08,02,03,01,00,06,07,05,04,34,35,36,37,32,33,31,30,3c,3d,3e,3f,3a,3b,39,38,2c,2d,2e,2f,2a,2b,29,28,22,23,21,20,26,27,25,24] + phi[1e]
		2a	p[2d,2c,2f,2e,2b,2a,28,29,23,22,20,21,27,26,24,25,3d,3c,3f,3e,3b,3a,38,39,33,32,30,31,37,36,34,35,1d,1c,1f,1e,1b,1a,18,19,13,12,10,11,17,16,14,15,03,02,00,01,07,06,04,05,0b,0a,08,09,0f,0e,0c,0d] + phi[1d]
		2b	p[06,07,05,04,01,00,03,02,0e,0f,0d,0c,09,08,0b,0a,16,17,15,14,11,10,13,12,1e,1f,1d,1c,19,18,1b,1a,26,27,25,24,21,20,23,22,2e,2f,2d,2c,29,28,2b,2a,36,37,35,34,31,30,33,32,3e,3f,3d,3c,39,38,3b,3a] + phi[1c]
		2c	p[30,31,32,33,34,35,36,37,38,39,3a,3b,3c,3d,3e,3f,28,29,2a,2b,2c,2d,2e,2f,24,25,26,27,22,23,21,20,08,09,0a,0b,0c,0d,0e,0f,04,05,06,07,02,03,01,00,18,19,1a,1b,1c,1d,1e,1f,14,15,16,17,12,13,11,10] + phi[1b]
		2d	p[17,16,14,15,10,11,12,13,1f,1e,1c,1d,18,19,1a,1b,0f,0e,0c,0d,08,09,0a,0b,00,01,02,03,04,05,06,07,37,36,34,35,30,31,32,33,3f,3e,3c,3d,38,39,3a,3b,2f,2e,2c,2d,28,29,2a,2b,20,21,22,23,24,25,26,27] + phi[1a]
		2e	p[29,28,2b,2a,2d,2c,2f,2e,25,24,27,26,23,22,20,21,39,38,3b,3a,3d,3c,3f,3e,35,34,37,36,33,32,30,31,19,18,1b,1a,1d,1c,1f,1e,15,14,17,16,13,12,10,11,05,04,07,06,03,02,00,01,0d,0c,0f,0e,0b,0a,08,09] + phi[19]
		2f	p[02,03,01,00,06,07,05,04,0a,0b,09,08,0e,0f,0d,0c,12,13,11,10,16,17,15,14,1a,1b,19,18,1e,1f,1d,1c,22,23,21,20,26,27,25,24,2a,2b,29,28,2e,2f,2d,2c,32,33,31,30,36,37,35,34,3a,3b,39,38,3e,3f,3d,3c] + phi[18]
		30	p[33,32,30,31,37,36,34,35,3b,3a,38,39,3f,3e,3c,3d,2b,2a,28,29,2f,2e,2c,2d,27,26,24,25,20,21,22,23,0b,0a,08,09,0f,0e,0c,0d,07,06,04,05,00,01,02,03,1b,1a,18,19,1f,1e,1c,1d,17,16,14,15,10,11,12,13] + phi[17]
		31	p[1e,1f,1d,1c,19,18,1b,1a,11,10,13,12,15,14,17,16,01,00,03,02,05,04,07,06,09,08,0b,0a,0d,0c,0f,0e,3e,3f,3d,3c,39,38,3b,3a,31,30,33,32,35,34,37,36,21,20,23,22,25,24,27,26,29,28,2b,2a,2d,2c,2f,2e] + phi[16]
		32	p[24,25,26,27,22,23,21,20,2c,2d,2e,2f,2a,2b,29,28,34,35,36,37,32,33,31,30,3c,3d,3e,3f,3a,3b,39,38,14,15,16,17,12,13,11,10,1c,1d,1e,1f,1a,1b,19,18,0c,0d,0e,0f,0a,0b,09,08,02,03,01,00,06,07,05,04] + phi[15]
		33	p[0d,0c,0f,0e,0b,0a,08,09,03,02,00,01,07,06,04,05,1d,1c,1f,1e,1b,1a,18,19,13,12,10,11,17,16,14,15,2d,2c,2f,2e,2b,2a,28,29,23,22,20,21,27,26,24,25,3d,3c,3f,3e,3b,3a,38,39,33,32,30,31,37,36,34,35] + phi[14]
		34	p[3a,3b,39,38,3e,3f,3d,3c,36,37,35,34,31,30,33,32,26,27,25,24,21,20,23,22,2e,2f,2d,2c,29,28,2b,2a,06,07,05,04,01,00,03,02,0e,0f,0d,0c,09,08,0b,0a,16,17,15,14,11,10,13,12,1e,1f,1d,1c,19,18,1b,1a] + phi[13]
		35	p[10,11,12,13,14,15,16,17,18,19,1a,1b,1c,1d,1e,1f,08,09,0a,0b,0c,0d,0e,0f,04,05,06,07,02,03,01,00,30,31,32,33,34,35,36,37,38,39,3a,3b,3c,3d,3e,3f,28,29,2a,2b,2c,2d,2e,2f,24,25,26,27,22,23,21,20] + phi[12]
		36	p[27,26,24,25,20,21,22,23,2f,2e,2c,2d,28,29,2a,2b,37,36,34,35,30,31,32,33,3f,3e,3c,3d,38,39,3a,3b,17,16,14,15,10,11,12,13,1f,1e,1c,1d,18,19,1a,1b,0f,0e,0c,0d,08,09,0a,0b,00,01,02,03,04,05,06,07] + phi[11]
		37	p[09,08,0b,0a,0d,0c,0f,0e,05,04,07,06,03,02,00,01,19,18,1b,1a,1d,1c,1f,1e,15,14,17,16,13,12,10,11,29,28,2b,2a,2d,2c,2f,2e,25,24,27,26,23,22,20,21,39,38,3b,3a,3d,3c,3f,3e,35,34,37,36,33,32,30,31] + phi[10]
		38	p[3c,3d,3e,3f,3a,3b,39,38,32,33,31,30,36,37,35,34,22,23,21,20,26,27,25,24,2a,2b,29,28,2e,2f,2d,2c,02,03,01,00,06,07,05,04,0a,0b,09,08,0e,0f,0d,0c,12,13,11,10,16,17,15,14,1a,1b,19,18,1e,1f,1d,1c] + phi[0f]
		39	p[13,12,10,11,17,16,14,15,1b,1a,18,19,1f,1e,1c,1d,0b,0a,08,09,0f,0e,0c,0d,07,06,04,05,00,01,02,03,33,32,30,31,37,36,34,35,3b,3a,38,39,3f,3e,3c,3d,2b,2a,28,29,2f,2e,2c,2d,27,26,24,25,20,21,22,23] + phi[0e]
		3a	p[2e,2f,2d,2c,29,28,2b,2a,21,20,23,22,25,24,27,26,3e,3f,3d,3c,39,38,3b,3a,31,30,33,32,35,34,37,36,1e,1f,1d,1c,19,18,1b,1a,11,10,13,12,15,14,17,16,01,00,03,02,05,04,07,06,09,08,0b,0a,0d,0c,0f,0e] + phi[0d]
		3b	p[04,05,06,07,02,03,01,00,0c,0d,0e,0f,0a,0b,09,08,14,15,16,17,12,13,11,10,1c,1d,1e,1f,1a,1b,19,18,24,25,26,27,22,23,21,20,2c,2d,2e,2f,2a,2b,29,28,34,35,36,37,32,33,31,30,3c,3d,3e,3f,3a,3b,39,38] + phi[0c]
		3c	p[35,34,37,36,33,32,30,31,3d,3c,3f,3e,3b,3a,38,39,2d,2c,2f,2e,2b,2a,28,29,23,22,20,21,27,26,24,25,0d,0c,0f,0e,0b,0a,08,09,03,02,00,01,07,06,04,05,1d,1c,1f,1e,1b,1a,18,19,13,12,10,11,17,16,14,15] + phi[0b]
		3d	p[1a,1b,19,18,1e,1f,1d,1c,16,17,15,14,11,10,13,12,06,07,05,04,01,00,03,02,0e,0f,0d,0c,09,08,0b,0a,3a,3b,39,38,3e,3f,3d,3c,36,37,35,34,31,30,33,32,26,27,25,24,21,20,23,22,2e,2f,2d,2c,29,28,2b,2a] + phi[0a]
		3e	p[20,21,22,23,24,25,26,27,28,29,2a,2b,2c,2d,2e,2f,30,31,32,33,34,35,36,37,38,39,3a,3b,3c,3d,3e,3f,10,11,12,13,14,15,16,17,18,19,1a,1b,1c,1d,1e,1f,08,09,0a,0b,0c,0d,0e,0f,04,05,06,07,02,03,01,00] + phi[09]

		While the phi-indices may have a not-wildly-badly-exploitable pattern, also simplest just to use a byte-array.
		*/
		tptr = t;
		for(k = 0; k < ODD_RADIX; ++k) {
		/* Initial-trial code for auto-gen of operm:
			jt = j1 + phi[k];
			RADIX_64_DIF((double *)tptr,t_offsets,1, (a+jt),o_offsets,RE_IM_STRIDE);
		*/
			iptr = dif64_oidx_lo + (k<<6);
			for(l = 0; l < 64; l++)
			{
				io_offsets[l] = plo[*iptr++];
			}
			jt = j1 + phi[dif_phi[k]];
			RADIX_64_DIF((double *)tptr,t_offsets,1, (a+jt),io_offsets,RE_IM_STRIDE);
			tptr++;
		}
	}
}

/***************/

void radix4032_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-4032 = 63x16 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix4032_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
*/
	int k,l, j,j1,j2,jt;
	static int NDIVR,first_entry=TRUE;
	struct complex t[RADIX], *tptr;
	// Need storage for circular-shifts perms of a basic 63-vector, with shift count in [0,15] that means 63+15 elts:
	const uint8 *iptr, dit_p40_cperms[128] = {	// Pad byte-array to next-higher 8-byte multiple
		0,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01,
		0,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01
	},	dit_phi[64] = {	// Only 63 entries, but again pad out to nearest 8-multiple
		0,0x02,0x01,0x08,0x07,0x06,0x05,0x04,0x03,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09
	};
	static int io_offsets[64], t_offsets[64];
	static int plo[64], phi[ODD_RADIX], toff[ODD_RADIX];
	// Low parts [p0-3f] of output-index perms - need 64 bytes (one init row below) for each radix-64 DFT:
	/**** NOTE that this is just a column-permuted version of dif64_oidx_lo[], 1st index row gives the column permutation ****/
	const uint8 dit64_iidx_lo[64] = {	// Keep only first row...
		0x00,0x01,0x03,0x02,0x07,0x06,0x05,0x04,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x3f,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20
	};	// ...and use it to access same array as used for DIF-input scramble - the main carry-loop routine thus needs just one full-sized byte-array:
	/*** Now done via #include "radix4032.h" @top of this file ***/	// 4032-element dif64_oidx_lo byte-array

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR = n/RADIX;
		/* Constant padded-array index offsets for array load/stores are here. */
		for(l = 0; l < 64; l++)
		{
			plo[l] = l*NDIVR;
			plo[l] += ( (plo[l] >> DAT_BITS) << PAD_BITS );
			t_offsets[l] = ((l<<6)-l)<<1;	// 2*(l*63), extra 2x is due to cast-to-double of t-array in call to RADIX_64_DIT
		//	i_offsets[l] = plo[l];	// Initial pass (which we use to auto-gen the needed operm) values of i_offsets
		}
		for(l = 0; l < ODD_RADIX; l++)
		{
			phi[l] = (l<<6)*NDIVR;
			phi[l] += ( (phi[l] >> DAT_BITS) << PAD_BITS );
			toff[l] = l+l;
		}
	}

/*...The radix-4032 pass is here.	*/

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
	Gather the needed data (4032 64-bit complex) and do 63 radix-16 transforms:

	Twiddleless version uses same linear-index-vector-form permutation as in DIF -
	Remember, inputs to DIT are bit-reversed, so use output of test_fft_radix() with
	TTYPE=0 to auto-generate needed input-index permutation:

	Combined DIT input-scramble array = [in hex]
		000,001,003,002,007,006,005,004,00f,00e,00d,00c,00b,00a,009,008,01f,01e,01d,01c,01b,01a,019,018,017,016,015,014,013,012,011,010,03f,03e,03d,03c,03b,03a,039,038,037,036,035,034,033,032,031,030,02f,02e,02d,02c,02b,02a,029,028,027,026,025,024,023,022,021,020,
		095,094,096,097,091,090,092,093,099,098,09a,09b,09e,09f,09c,09d,085,084,086,087,081,080,082,083,089,088,08a,08b,08e,08f,08c,08d,0a5,0a4,0a6,0a7,0a1,0a0,0a2,0a3,0a9,0a8,0aa,0ab,0ae,0af,0ac,0ad,0b9,0b8,0ba,0bb,0be,0bf,0bc,0bd,0b1,0b0,0b2,0b3,0b6,0b7,0b4,0b5,
		06a,06b,068,069,06c,06d,06f,06e,062,063,060,061,064,065,067,066,072,073,070,071,074,075,077,076,07c,07d,07f,07e,078,079,07b,07a,04a,04b,048,049,04c,04d,04f,04e,042,043,040,041,044,045,047,046,052,053,050,051,054,055,057,056,05c,05d,05f,05e,058,059,05b,05a,
		207,206,205,204,203,202,201,200,20b,20a,209,208,20d,20c,20e,20f,21b,21a,219,218,21d,21c,21e,21f,213,212,211,210,215,214,216,217,23b,23a,239,238,23d,23c,23e,23f,233,232,231,230,235,234,236,237,22b,22a,229,228,22d,22c,22e,22f,223,222,221,220,225,224,226,227,
		1f1,1f0,1f2,1f3,1f6,1f7,1f4,1f5,1fe,1ff,1fc,1fd,1fa,1fb,1f8,1f9,1e1,1e0,1e2,1e3,1e6,1e7,1e4,1e5,1ee,1ef,1ec,1ed,1ea,1eb,1e8,1e9,1d1,1d0,1d2,1d3,1d6,1d7,1d4,1d5,1de,1df,1dc,1dd,1da,1db,1d8,1d9,1c1,1c0,1c2,1c3,1c6,1c7,1c4,1c5,1ce,1cf,1cc,1cd,1ca,1cb,1c8,1c9,
		19c,19d,19f,19e,198,199,19b,19a,194,195,197,196,190,191,193,192,18c,18d,18f,18e,188,189,18b,18a,184,185,187,186,180,181,183,182,1ac,1ad,1af,1ae,1a8,1a9,1ab,1aa,1a4,1a5,1a7,1a6,1a0,1a1,1a3,1a2,1b4,1b5,1b7,1b6,1b0,1b1,1b3,1b2,1b8,1b9,1bb,1ba,1bf,1be,1bd,1bc,
		163,162,161,160,165,164,166,167,16d,16c,16e,16f,169,168,16a,16b,17d,17c,17e,17f,179,178,17a,17b,175,174,176,177,171,170,172,173,143,142,141,140,145,144,146,147,14d,14c,14e,14f,149,148,14a,14b,15d,15c,15e,15f,159,158,15a,15b,155,154,156,157,151,150,152,153,
		10e,10f,10c,10d,10a,10b,108,109,106,107,104,105,102,103,100,101,116,117,114,115,112,113,110,111,11a,11b,118,119,11c,11d,11f,11e,136,137,134,135,132,133,130,131,13a,13b,138,139,13c,13d,13f,13e,126,127,124,125,122,123,120,121,12a,12b,128,129,12c,12d,12f,12e,
		0f8,0f9,0fb,0fa,0ff,0fe,0fd,0fc,0f0,0f1,0f3,0f2,0f7,0f6,0f5,0f4,0e8,0e9,0eb,0ea,0ef,0ee,0ed,0ec,0e0,0e1,0e3,0e2,0e7,0e6,0e5,0e4,0d8,0d9,0db,0da,0df,0de,0dd,0dc,0d0,0d1,0d3,0d2,0d7,0d6,0d5,0d4,0c8,0c9,0cb,0ca,0cf,0ce,0cd,0cc,0c0,0c1,0c3,0c2,0c7,0c6,0c5,0c4,
		f9f,f9e,f9d,f9c,f9b,f9a,f99,f98,f97,f96,f95,f94,f93,f92,f91,f90,f8f,f8e,f8d,f8c,f8b,f8a,f89,f88,f87,f86,f85,f84,f83,f82,f81,f80,faf,fae,fad,fac,fab,faa,fa9,fa8,fa7,fa6,fa5,fa4,fa3,fa2,fa1,fa0,fb7,fb6,fb5,fb4,fb3,fb2,fb1,fb0,fbb,fba,fb9,fb8,fbd,fbc,fbe,fbf,
		f65,f64,f66,f67,f61,f60,f62,f63,f69,f68,f6a,f6b,f6e,f6f,f6c,f6d,f79,f78,f7a,f7b,f7e,f7f,f7c,f7d,f71,f70,f72,f73,f76,f77,f74,f75,f45,f44,f46,f47,f41,f40,f42,f43,f49,f48,f4a,f4b,f4e,f4f,f4c,f4d,f59,f58,f5a,f5b,f5e,f5f,f5c,f5d,f51,f50,f52,f53,f56,f57,f54,f55,
		f0a,f0b,f08,f09,f0c,f0d,f0f,f0e,f02,f03,f00,f01,f04,f05,f07,f06,f12,f13,f10,f11,f14,f15,f17,f16,f1c,f1d,f1f,f1e,f18,f19,f1b,f1a,f32,f33,f30,f31,f34,f35,f37,f36,f3c,f3d,f3f,f3e,f38,f39,f3b,f3a,f22,f23,f20,f21,f24,f25,f27,f26,f2c,f2d,f2f,f2e,f28,f29,f2b,f2a,
		efb,efa,ef9,ef8,efd,efc,efe,eff,ef3,ef2,ef1,ef0,ef5,ef4,ef6,ef7,eeb,eea,ee9,ee8,eed,eec,eee,eef,ee3,ee2,ee1,ee0,ee5,ee4,ee6,ee7,edb,eda,ed9,ed8,edd,edc,ede,edf,ed3,ed2,ed1,ed0,ed5,ed4,ed6,ed7,ecb,eca,ec9,ec8,ecd,ecc,ece,ecf,ec3,ec2,ec1,ec0,ec5,ec4,ec6,ec7,
		e91,e90,e92,e93,e96,e97,e94,e95,e9e,e9f,e9c,e9d,e9a,e9b,e98,e99,e81,e80,e82,e83,e86,e87,e84,e85,e8e,e8f,e8c,e8d,e8a,e8b,e88,e89,ea1,ea0,ea2,ea3,ea6,ea7,ea4,ea5,eae,eaf,eac,ead,eaa,eab,ea8,ea9,ebe,ebf,ebc,ebd,eba,ebb,eb8,eb9,eb6,eb7,eb4,eb5,eb2,eb3,eb0,eb1,
		e6c,e6d,e6f,e6e,e68,e69,e6b,e6a,e64,e65,e67,e66,e60,e61,e63,e62,e74,e75,e77,e76,e70,e71,e73,e72,e78,e79,e7b,e7a,e7f,e7e,e7d,e7c,e4c,e4d,e4f,e4e,e48,e49,e4b,e4a,e44,e45,e47,e46,e40,e41,e43,e42,e54,e55,e57,e56,e50,e51,e53,e52,e58,e59,e5b,e5a,e5f,e5e,e5d,e5c,
		e03,e02,e01,e00,e05,e04,e06,e07,e0d,e0c,e0e,e0f,e09,e08,e0a,e0b,e1d,e1c,e1e,e1f,e19,e18,e1a,e1b,e15,e14,e16,e17,e11,e10,e12,e13,e3d,e3c,e3e,e3f,e39,e38,e3a,e3b,e35,e34,e36,e37,e31,e30,e32,e33,e2d,e2c,e2e,e2f,e29,e28,e2a,e2b,e25,e24,e26,e27,e21,e20,e22,e23,
		df6,df7,df4,df5,df2,df3,df0,df1,dfa,dfb,df8,df9,dfc,dfd,dff,dfe,de6,de7,de4,de5,de2,de3,de0,de1,dea,deb,de8,de9,dec,ded,def,dee,dd6,dd7,dd4,dd5,dd2,dd3,dd0,dd1,dda,ddb,dd8,dd9,ddc,ddd,ddf,dde,dc6,dc7,dc4,dc5,dc2,dc3,dc0,dc1,dca,dcb,dc8,dc9,dcc,dcd,dcf,dce,
		d98,d99,d9b,d9a,d9f,d9e,d9d,d9c,d90,d91,d93,d92,d97,d96,d95,d94,d88,d89,d8b,d8a,d8f,d8e,d8d,d8c,d80,d81,d83,d82,d87,d86,d85,d84,da8,da9,dab,daa,daf,dae,dad,dac,da0,da1,da3,da2,da7,da6,da5,da4,db0,db1,db3,db2,db7,db6,db5,db4,dbf,dbe,dbd,dbc,dbb,dba,db9,db8,
		d6f,d6e,d6d,d6c,d6b,d6a,d69,d68,d67,d66,d65,d64,d63,d62,d61,d60,d77,d76,d75,d74,d73,d72,d71,d70,d7b,d7a,d79,d78,d7d,d7c,d7e,d7f,d4f,d4e,d4d,d4c,d4b,d4a,d49,d48,d47,d46,d45,d44,d43,d42,d41,d40,d57,d56,d55,d54,d53,d52,d51,d50,d5b,d5a,d59,d58,d5d,d5c,d5e,d5f,
		d05,d04,d06,d07,d01,d00,d02,d03,d09,d08,d0a,d0b,d0e,d0f,d0c,d0d,d19,d18,d1a,d1b,d1e,d1f,d1c,d1d,d11,d10,d12,d13,d16,d17,d14,d15,d39,d38,d3a,d3b,d3e,d3f,d3c,d3d,d31,d30,d32,d33,d36,d37,d34,d35,d29,d28,d2a,d2b,d2e,d2f,d2c,d2d,d21,d20,d22,d23,d26,d27,d24,d25,
		cf2,cf3,cf0,cf1,cf4,cf5,cf7,cf6,cfc,cfd,cff,cfe,cf8,cf9,cfb,cfa,ce2,ce3,ce0,ce1,ce4,ce5,ce7,ce6,cec,ced,cef,cee,ce8,ce9,ceb,cea,cd2,cd3,cd0,cd1,cd4,cd5,cd7,cd6,cdc,cdd,cdf,cde,cd8,cd9,cdb,cda,cc2,cc3,cc0,cc1,cc4,cc5,cc7,cc6,ccc,ccd,ccf,cce,cc8,cc9,ccb,cca,
		c9b,c9a,c99,c98,c9d,c9c,c9e,c9f,c93,c92,c91,c90,c95,c94,c96,c97,c8b,c8a,c89,c88,c8d,c8c,c8e,c8f,c83,c82,c81,c80,c85,c84,c86,c87,cab,caa,ca9,ca8,cad,cac,cae,caf,ca3,ca2,ca1,ca0,ca5,ca4,ca6,ca7,cb3,cb2,cb1,cb0,cb5,cb4,cb6,cb7,cbd,cbc,cbe,cbf,cb9,cb8,cba,cbb,
		c61,c60,c62,c63,c66,c67,c64,c65,c6e,c6f,c6c,c6d,c6a,c6b,c68,c69,c7e,c7f,c7c,c7d,c7a,c7b,c78,c79,c76,c77,c74,c75,c72,c73,c70,c71,c41,c40,c42,c43,c46,c47,c44,c45,c4e,c4f,c4c,c4d,c4a,c4b,c48,c49,c5e,c5f,c5c,c5d,c5a,c5b,c58,c59,c56,c57,c54,c55,c52,c53,c50,c51,
		c0c,c0d,c0f,c0e,c08,c09,c0b,c0a,c04,c05,c07,c06,c00,c01,c03,c02,c14,c15,c17,c16,c10,c11,c13,c12,c18,c19,c1b,c1a,c1f,c1e,c1d,c1c,c34,c35,c37,c36,c30,c31,c33,c32,c38,c39,c3b,c3a,c3f,c3e,c3d,c3c,c24,c25,c27,c26,c20,c21,c23,c22,c28,c29,c2b,c2a,c2f,c2e,c2d,c2c,
		bfd,bfc,bfe,bff,bf9,bf8,bfa,bfb,bf5,bf4,bf6,bf7,bf1,bf0,bf2,bf3,bed,bec,bee,bef,be9,be8,bea,beb,be5,be4,be6,be7,be1,be0,be2,be3,bdd,bdc,bde,bdf,bd9,bd8,bda,bdb,bd5,bd4,bd6,bd7,bd1,bd0,bd2,bd3,bcd,bcc,bce,bcf,bc9,bc8,bca,bcb,bc5,bc4,bc6,bc7,bc1,bc0,bc2,bc3,
		b96,b97,b94,b95,b92,b93,b90,b91,b9a,b9b,b98,b99,b9c,b9d,b9f,b9e,b86,b87,b84,b85,b82,b83,b80,b81,b8a,b8b,b88,b89,b8c,b8d,b8f,b8e,ba6,ba7,ba4,ba5,ba2,ba3,ba0,ba1,baa,bab,ba8,ba9,bac,bad,baf,bae,bba,bbb,bb8,bb9,bbc,bbd,bbf,bbe,bb2,bb3,bb0,bb1,bb4,bb5,bb7,bb6,
		b68,b69,b6b,b6a,b6f,b6e,b6d,b6c,b60,b61,b63,b62,b67,b66,b65,b64,b70,b71,b73,b72,b77,b76,b75,b74,b7f,b7e,b7d,b7c,b7b,b7a,b79,b78,b48,b49,b4b,b4a,b4f,b4e,b4d,b4c,b40,b41,b43,b42,b47,b46,b45,b44,b50,b51,b53,b52,b57,b56,b55,b54,b5f,b5e,b5d,b5c,b5b,b5a,b59,b58,
		b0f,b0e,b0d,b0c,b0b,b0a,b09,b08,b07,b06,b05,b04,b03,b02,b01,b00,b17,b16,b15,b14,b13,b12,b11,b10,b1b,b1a,b19,b18,b1d,b1c,b1e,b1f,b37,b36,b35,b34,b33,b32,b31,b30,b3b,b3a,b39,b38,b3d,b3c,b3e,b3f,b27,b26,b25,b24,b23,b22,b21,b20,b2b,b2a,b29,b28,b2d,b2c,b2e,b2f,
		af9,af8,afa,afb,afe,aff,afc,afd,af1,af0,af2,af3,af6,af7,af4,af5,ae9,ae8,aea,aeb,aee,aef,aec,aed,ae1,ae0,ae2,ae3,ae6,ae7,ae4,ae5,ad9,ad8,ada,adb,ade,adf,adc,add,ad1,ad0,ad2,ad3,ad6,ad7,ad4,ad5,ac9,ac8,aca,acb,ace,acf,acc,acd,ac1,ac0,ac2,ac3,ac6,ac7,ac4,ac5,
		a92,a93,a90,a91,a94,a95,a97,a96,a9c,a9d,a9f,a9e,a98,a99,a9b,a9a,a82,a83,a80,a81,a84,a85,a87,a86,a8c,a8d,a8f,a8e,a88,a89,a8b,a8a,aa2,aa3,aa0,aa1,aa4,aa5,aa7,aa6,aac,aad,aaf,aae,aa8,aa9,aab,aaa,abc,abd,abf,abe,ab8,ab9,abb,aba,ab4,ab5,ab7,ab6,ab0,ab1,ab3,ab2,
		a6b,a6a,a69,a68,a6d,a6c,a6e,a6f,a63,a62,a61,a60,a65,a64,a66,a67,a73,a72,a71,a70,a75,a74,a76,a77,a7d,a7c,a7e,a7f,a79,a78,a7a,a7b,a4b,a4a,a49,a48,a4d,a4c,a4e,a4f,a43,a42,a41,a40,a45,a44,a46,a47,a53,a52,a51,a50,a55,a54,a56,a57,a5d,a5c,a5e,a5f,a59,a58,a5a,a5b,
		a01,a00,a02,a03,a06,a07,a04,a05,a0e,a0f,a0c,a0d,a0a,a0b,a08,a09,a1e,a1f,a1c,a1d,a1a,a1b,a18,a19,a16,a17,a14,a15,a12,a13,a10,a11,a3e,a3f,a3c,a3d,a3a,a3b,a38,a39,a36,a37,a34,a35,a32,a33,a30,a31,a2e,a2f,a2c,a2d,a2a,a2b,a28,a29,a26,a27,a24,a25,a22,a23,a20,a21,
		9f4,9f5,9f7,9f6,9f0,9f1,9f3,9f2,9f8,9f9,9fb,9fa,9ff,9fe,9fd,9fc,9e4,9e5,9e7,9e6,9e0,9e1,9e3,9e2,9e8,9e9,9eb,9ea,9ef,9ee,9ed,9ec,9d4,9d5,9d7,9d6,9d0,9d1,9d3,9d2,9d8,9d9,9db,9da,9df,9de,9dd,9dc,9c4,9c5,9c7,9c6,9c0,9c1,9c3,9c2,9c8,9c9,9cb,9ca,9cf,9ce,9cd,9cc,
		99d,99c,99e,99f,999,998,99a,99b,995,994,996,997,991,990,992,993,98d,98c,98e,98f,989,988,98a,98b,985,984,986,987,981,980,982,983,9ad,9ac,9ae,9af,9a9,9a8,9aa,9ab,9a5,9a4,9a6,9a7,9a1,9a0,9a2,9a3,9b5,9b4,9b6,9b7,9b1,9b0,9b2,9b3,9b9,9b8,9ba,9bb,9be,9bf,9bc,9bd,
		966,967,964,965,962,963,960,961,96a,96b,968,969,96c,96d,96f,96e,97a,97b,978,979,97c,97d,97f,97e,972,973,970,971,974,975,977,976,946,947,944,945,942,943,940,941,94a,94b,948,949,94c,94d,94f,94e,95a,95b,958,959,95c,95d,95f,95e,952,953,950,951,954,955,957,956,
		908,909,90b,90a,90f,90e,90d,90c,900,901,903,902,907,906,905,904,910,911,913,912,917,916,915,914,91f,91e,91d,91c,91b,91a,919,918,930,931,933,932,937,936,935,934,93f,93e,93d,93c,93b,93a,939,938,920,921,923,922,927,926,925,924,92f,92e,92d,92c,92b,92a,929,928,
		8f7,8f6,8f5,8f4,8f3,8f2,8f1,8f0,8fb,8fa,8f9,8f8,8fd,8fc,8fe,8ff,8e7,8e6,8e5,8e4,8e3,8e2,8e1,8e0,8eb,8ea,8e9,8e8,8ed,8ec,8ee,8ef,8d7,8d6,8d5,8d4,8d3,8d2,8d1,8d0,8db,8da,8d9,8d8,8dd,8dc,8de,8df,8c7,8c6,8c5,8c4,8c3,8c2,8c1,8c0,8cb,8ca,8c9,8c8,8cd,8cc,8ce,8cf,
		899,898,89a,89b,89e,89f,89c,89d,891,890,892,893,896,897,894,895,889,888,88a,88b,88e,88f,88c,88d,881,880,882,883,886,887,884,885,8a9,8a8,8aa,8ab,8ae,8af,8ac,8ad,8a1,8a0,8a2,8a3,8a6,8a7,8a4,8a5,8b1,8b0,8b2,8b3,8b6,8b7,8b4,8b5,8be,8bf,8bc,8bd,8ba,8bb,8b8,8b9,
		862,863,860,861,864,865,867,866,86c,86d,86f,86e,868,869,86b,86a,87c,87d,87f,87e,878,879,87b,87a,874,875,877,876,870,871,873,872,842,843,840,841,844,845,847,846,84c,84d,84f,84e,848,849,84b,84a,85c,85d,85f,85e,858,859,85b,85a,854,855,857,856,850,851,853,852,
		80b,80a,809,808,80d,80c,80e,80f,803,802,801,800,805,804,806,807,813,812,811,810,815,814,816,817,81d,81c,81e,81f,819,818,81a,81b,833,832,831,830,835,834,836,837,83d,83c,83e,83f,839,838,83a,83b,823,822,821,820,825,824,826,827,82d,82c,82e,82f,829,828,82a,82b,
		7fe,7ff,7fc,7fd,7fa,7fb,7f8,7f9,7f6,7f7,7f4,7f5,7f2,7f3,7f0,7f1,7ee,7ef,7ec,7ed,7ea,7eb,7e8,7e9,7e6,7e7,7e4,7e5,7e2,7e3,7e0,7e1,7de,7df,7dc,7dd,7da,7db,7d8,7d9,7d6,7d7,7d4,7d5,7d2,7d3,7d0,7d1,7ce,7cf,7cc,7cd,7ca,7cb,7c8,7c9,7c6,7c7,7c4,7c5,7c2,7c3,7c0,7c1,
		794,795,797,796,790,791,793,792,798,799,79b,79a,79f,79e,79d,79c,784,785,787,786,780,781,783,782,788,789,78b,78a,78f,78e,78d,78c,7a4,7a5,7a7,7a6,7a0,7a1,7a3,7a2,7a8,7a9,7ab,7aa,7af,7ae,7ad,7ac,7b8,7b9,7bb,7ba,7bf,7be,7bd,7bc,7b0,7b1,7b3,7b2,7b7,7b6,7b5,7b4,
		76d,76c,76e,76f,769,768,76a,76b,765,764,766,767,761,760,762,763,775,774,776,777,771,770,772,773,779,778,77a,77b,77e,77f,77c,77d,74d,74c,74e,74f,749,748,74a,74b,745,744,746,747,741,740,742,743,755,754,756,757,751,750,752,753,759,758,75a,75b,75e,75f,75c,75d,
		706,707,704,705,702,703,700,701,70a,70b,708,709,70c,70d,70f,70e,71a,71b,718,719,71c,71d,71f,71e,712,713,710,711,714,715,717,716,73a,73b,738,739,73c,73d,73f,73e,732,733,730,731,734,735,737,736,72a,72b,728,729,72c,72d,72f,72e,722,723,720,721,724,725,727,726,
		6f0,6f1,6f3,6f2,6f7,6f6,6f5,6f4,6ff,6fe,6fd,6fc,6fb,6fa,6f9,6f8,6e0,6e1,6e3,6e2,6e7,6e6,6e5,6e4,6ef,6ee,6ed,6ec,6eb,6ea,6e9,6e8,6d0,6d1,6d3,6d2,6d7,6d6,6d5,6d4,6df,6de,6dd,6dc,6db,6da,6d9,6d8,6c0,6c1,6c3,6c2,6c7,6c6,6c5,6c4,6cf,6ce,6cd,6cc,6cb,6ca,6c9,6c8,
		697,696,695,694,693,692,691,690,69b,69a,699,698,69d,69c,69e,69f,687,686,685,684,683,682,681,680,68b,68a,689,688,68d,68c,68e,68f,6a7,6a6,6a5,6a4,6a3,6a2,6a1,6a0,6ab,6aa,6a9,6a8,6ad,6ac,6ae,6af,6bb,6ba,6b9,6b8,6bd,6bc,6be,6bf,6b3,6b2,6b1,6b0,6b5,6b4,6b6,6b7,
		669,668,66a,66b,66e,66f,66c,66d,661,660,662,663,666,667,664,665,671,670,672,673,676,677,674,675,67e,67f,67c,67d,67a,67b,678,679,649,648,64a,64b,64e,64f,64c,64d,641,640,642,643,646,647,644,645,651,650,652,653,656,657,654,655,65e,65f,65c,65d,65a,65b,658,659,
		602,603,600,601,604,605,607,606,60c,60d,60f,60e,608,609,60b,60a,61c,61d,61f,61e,618,619,61b,61a,614,615,617,616,610,611,613,612,63c,63d,63f,63e,638,639,63b,63a,634,635,637,636,630,631,633,632,62c,62d,62f,62e,628,629,62b,62a,624,625,627,626,620,621,623,622,
		5f3,5f2,5f1,5f0,5f5,5f4,5f6,5f7,5fd,5fc,5fe,5ff,5f9,5f8,5fa,5fb,5e3,5e2,5e1,5e0,5e5,5e4,5e6,5e7,5ed,5ec,5ee,5ef,5e9,5e8,5ea,5eb,5d3,5d2,5d1,5d0,5d5,5d4,5d6,5d7,5dd,5dc,5de,5df,5d9,5d8,5da,5db,5c3,5c2,5c1,5c0,5c5,5c4,5c6,5c7,5cd,5cc,5ce,5cf,5c9,5c8,5ca,5cb,
		59e,59f,59c,59d,59a,59b,598,599,596,597,594,595,592,593,590,591,58e,58f,58c,58d,58a,58b,588,589,586,587,584,585,582,583,580,581,5ae,5af,5ac,5ad,5aa,5ab,5a8,5a9,5a6,5a7,5a4,5a5,5a2,5a3,5a0,5a1,5b6,5b7,5b4,5b5,5b2,5b3,5b0,5b1,5ba,5bb,5b8,5b9,5bc,5bd,5bf,5be,
		564,565,567,566,560,561,563,562,568,569,56b,56a,56f,56e,56d,56c,578,579,57b,57a,57f,57e,57d,57c,570,571,573,572,577,576,575,574,544,545,547,546,540,541,543,542,548,549,54b,54a,54f,54e,54d,54c,558,559,55b,55a,55f,55e,55d,55c,550,551,553,552,557,556,555,554,
		50d,50c,50e,50f,509,508,50a,50b,505,504,506,507,501,500,502,503,515,514,516,517,511,510,512,513,519,518,51a,51b,51e,51f,51c,51d,535,534,536,537,531,530,532,533,539,538,53a,53b,53e,53f,53c,53d,525,524,526,527,521,520,522,523,529,528,52a,52b,52e,52f,52c,52d,
		4fa,4fb,4f8,4f9,4fc,4fd,4ff,4fe,4f2,4f3,4f0,4f1,4f4,4f5,4f7,4f6,4ea,4eb,4e8,4e9,4ec,4ed,4ef,4ee,4e2,4e3,4e0,4e1,4e4,4e5,4e7,4e6,4da,4db,4d8,4d9,4dc,4dd,4df,4de,4d2,4d3,4d0,4d1,4d4,4d5,4d7,4d6,4ca,4cb,4c8,4c9,4cc,4cd,4cf,4ce,4c2,4c3,4c0,4c1,4c4,4c5,4c7,4c6,
		490,491,493,492,497,496,495,494,49f,49e,49d,49c,49b,49a,499,498,480,481,483,482,487,486,485,484,48f,48e,48d,48c,48b,48a,489,488,4a0,4a1,4a3,4a2,4a7,4a6,4a5,4a4,4af,4ae,4ad,4ac,4ab,4aa,4a9,4a8,4bf,4be,4bd,4bc,4bb,4ba,4b9,4b8,4b7,4b6,4b5,4b4,4b3,4b2,4b1,4b0,
		467,466,465,464,463,462,461,460,46b,46a,469,468,46d,46c,46e,46f,47b,47a,479,478,47d,47c,47e,47f,473,472,471,470,475,474,476,477,447,446,445,444,443,442,441,440,44b,44a,449,448,44d,44c,44e,44f,45b,45a,459,458,45d,45c,45e,45f,453,452,451,450,455,454,456,457,
		409,408,40a,40b,40e,40f,40c,40d,401,400,402,403,406,407,404,405,411,410,412,413,416,417,414,415,41e,41f,41c,41d,41a,41b,418,419,431,430,432,433,436,437,434,435,43e,43f,43c,43d,43a,43b,438,439,421,420,422,423,426,427,424,425,42e,42f,42c,42d,42a,42b,428,429,
		3fc,3fd,3ff,3fe,3f8,3f9,3fb,3fa,3f4,3f5,3f7,3f6,3f0,3f1,3f3,3f2,3ec,3ed,3ef,3ee,3e8,3e9,3eb,3ea,3e4,3e5,3e7,3e6,3e0,3e1,3e3,3e2,3dc,3dd,3df,3de,3d8,3d9,3db,3da,3d4,3d5,3d7,3d6,3d0,3d1,3d3,3d2,3cc,3cd,3cf,3ce,3c8,3c9,3cb,3ca,3c4,3c5,3c7,3c6,3c0,3c1,3c3,3c2,
		393,392,391,390,395,394,396,397,39d,39c,39e,39f,399,398,39a,39b,383,382,381,380,385,384,386,387,38d,38c,38e,38f,389,388,38a,38b,3a3,3a2,3a1,3a0,3a5,3a4,3a6,3a7,3ad,3ac,3ae,3af,3a9,3a8,3aa,3ab,3bd,3bc,3be,3bf,3b9,3b8,3ba,3bb,3b5,3b4,3b6,3b7,3b1,3b0,3b2,3b3,
		36e,36f,36c,36d,36a,36b,368,369,366,367,364,365,362,363,360,361,376,377,374,375,372,373,370,371,37a,37b,378,379,37c,37d,37f,37e,34e,34f,34c,34d,34a,34b,348,349,346,347,344,345,342,343,340,341,356,357,354,355,352,353,350,351,35a,35b,358,359,35c,35d,35f,35e,
		304,305,307,306,300,301,303,302,308,309,30b,30a,30f,30e,30d,30c,318,319,31b,31a,31f,31e,31d,31c,310,311,313,312,317,316,315,314,338,339,33b,33a,33f,33e,33d,33c,330,331,333,332,337,336,335,334,328,329,32b,32a,32f,32e,32d,32c,320,321,323,322,327,326,325,324,
		2f5,2f4,2f6,2f7,2f1,2f0,2f2,2f3,2f9,2f8,2fa,2fb,2fe,2ff,2fc,2fd,2e5,2e4,2e6,2e7,2e1,2e0,2e2,2e3,2e9,2e8,2ea,2eb,2ee,2ef,2ec,2ed,2d5,2d4,2d6,2d7,2d1,2d0,2d2,2d3,2d9,2d8,2da,2db,2de,2df,2dc,2dd,2c5,2c4,2c6,2c7,2c1,2c0,2c2,2c3,2c9,2c8,2ca,2cb,2ce,2cf,2cc,2cd,
		29a,29b,298,299,29c,29d,29f,29e,292,293,290,291,294,295,297,296,28a,28b,288,289,28c,28d,28f,28e,282,283,280,281,284,285,287,286,2aa,2ab,2a8,2a9,2ac,2ad,2af,2ae,2a2,2a3,2a0,2a1,2a4,2a5,2a7,2a6,2b2,2b3,2b0,2b1,2b4,2b5,2b7,2b6,2bc,2bd,2bf,2be,2b8,2b9,2bb,2ba,
		260,261,263,262,267,266,265,264,26f,26e,26d,26c,26b,26a,269,268,27f,27e,27d,27c,27b,27a,279,278,277,276,275,274,273,272,271,270,240,241,243,242,247,246,245,244,24f,24e,24d,24c,24b,24a,249,248,25f,25e,25d,25c,25b,25a,259,258,257,256,255,254,253,252,251,250,

	Separating the above into (mod 64) - indices into plo[] array - and (/= 64) - indices into phi[] (multiples of p40) array -
	Note the phi index pattern is identical to radix-4032 DIF, except here phi has 64-multiples of plo, rather than 16-multiples:

		row:						o-Indices:
		00	p[00,01,03,02,07,06,05,04,0f,0e,0d,0c,0b,0a,09,08,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,3f,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20] + phi[00]
		01	p[15,14,16,17,11,10,12,13,19,18,1a,1b,1e,1f,1c,1d,05,04,06,07,01,00,02,03,09,08,0a,0b,0e,0f,0c,0d,25,24,26,27,21,20,22,23,29,28,2a,2b,2e,2f,2c,2d,39,38,3a,3b,3e,3f,3c,3d,31,30,32,33,36,37,34,35] + phi[02]
		02	p[2a,2b,28,29,2c,2d,2f,2e,22,23,20,21,24,25,27,26,32,33,30,31,34,35,37,36,3c,3d,3f,3e,38,39,3b,3a,0a,0b,08,09,0c,0d,0f,0e,02,03,00,01,04,05,07,06,12,13,10,11,14,15,17,16,1c,1d,1f,1e,18,19,1b,1a] + phi[01]
		03	p[07,06,05,04,03,02,01,00,0b,0a,09,08,0d,0c,0e,0f,1b,1a,19,18,1d,1c,1e,1f,13,12,11,10,15,14,16,17,3b,3a,39,38,3d,3c,3e,3f,33,32,31,30,35,34,36,37,2b,2a,29,28,2d,2c,2e,2f,23,22,21,20,25,24,26,27] + phi[08]
		04	p[31,30,32,33,36,37,34,35,3e,3f,3c,3d,3a,3b,38,39,21,20,22,23,26,27,24,25,2e,2f,2c,2d,2a,2b,28,29,11,10,12,13,16,17,14,15,1e,1f,1c,1d,1a,1b,18,19,01,00,02,03,06,07,04,05,0e,0f,0c,0d,0a,0b,08,09] + phi[07]
		05	p[1c,1d,1f,1e,18,19,1b,1a,14,15,17,16,10,11,13,12,0c,0d,0f,0e,08,09,0b,0a,04,05,07,06,00,01,03,02,2c,2d,2f,2e,28,29,2b,2a,24,25,27,26,20,21,23,22,34,35,37,36,30,31,33,32,38,39,3b,3a,3f,3e,3d,3c] + phi[06]
		06	p[23,22,21,20,25,24,26,27,2d,2c,2e,2f,29,28,2a,2b,3d,3c,3e,3f,39,38,3a,3b,35,34,36,37,31,30,32,33,03,02,01,00,05,04,06,07,0d,0c,0e,0f,09,08,0a,0b,1d,1c,1e,1f,19,18,1a,1b,15,14,16,17,11,10,12,13] + phi[05]
		07	p[0e,0f,0c,0d,0a,0b,08,09,06,07,04,05,02,03,00,01,16,17,14,15,12,13,10,11,1a,1b,18,19,1c,1d,1f,1e,36,37,34,35,32,33,30,31,3a,3b,38,39,3c,3d,3f,3e,26,27,24,25,22,23,20,21,2a,2b,28,29,2c,2d,2f,2e] + phi[04]
		08	p[38,39,3b,3a,3f,3e,3d,3c,30,31,33,32,37,36,35,34,28,29,2b,2a,2f,2e,2d,2c,20,21,23,22,27,26,25,24,18,19,1b,1a,1f,1e,1d,1c,10,11,13,12,17,16,15,14,08,09,0b,0a,0f,0e,0d,0c,00,01,03,02,07,06,05,04] + phi[03]
		09	p[1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,37,36,35,34,33,32,31,30,3b,3a,39,38,3d,3c,3e,3f] + phi[3e]
		0a	p[25,24,26,27,21,20,22,23,29,28,2a,2b,2e,2f,2c,2d,39,38,3a,3b,3e,3f,3c,3d,31,30,32,33,36,37,34,35,05,04,06,07,01,00,02,03,09,08,0a,0b,0e,0f,0c,0d,19,18,1a,1b,1e,1f,1c,1d,11,10,12,13,16,17,14,15] + phi[3d]
		0b	p[0a,0b,08,09,0c,0d,0f,0e,02,03,00,01,04,05,07,06,12,13,10,11,14,15,17,16,1c,1d,1f,1e,18,19,1b,1a,32,33,30,31,34,35,37,36,3c,3d,3f,3e,38,39,3b,3a,22,23,20,21,24,25,27,26,2c,2d,2f,2e,28,29,2b,2a] + phi[3c]
		0c	p[3b,3a,39,38,3d,3c,3e,3f,33,32,31,30,35,34,36,37,2b,2a,29,28,2d,2c,2e,2f,23,22,21,20,25,24,26,27,1b,1a,19,18,1d,1c,1e,1f,13,12,11,10,15,14,16,17,0b,0a,09,08,0d,0c,0e,0f,03,02,01,00,05,04,06,07] + phi[3b]
		0d	p[11,10,12,13,16,17,14,15,1e,1f,1c,1d,1a,1b,18,19,01,00,02,03,06,07,04,05,0e,0f,0c,0d,0a,0b,08,09,21,20,22,23,26,27,24,25,2e,2f,2c,2d,2a,2b,28,29,3e,3f,3c,3d,3a,3b,38,39,36,37,34,35,32,33,30,31] + phi[3a]
		0e	p[2c,2d,2f,2e,28,29,2b,2a,24,25,27,26,20,21,23,22,34,35,37,36,30,31,33,32,38,39,3b,3a,3f,3e,3d,3c,0c,0d,0f,0e,08,09,0b,0a,04,05,07,06,00,01,03,02,14,15,17,16,10,11,13,12,18,19,1b,1a,1f,1e,1d,1c] + phi[39]
		0f	p[03,02,01,00,05,04,06,07,0d,0c,0e,0f,09,08,0a,0b,1d,1c,1e,1f,19,18,1a,1b,15,14,16,17,11,10,12,13,3d,3c,3e,3f,39,38,3a,3b,35,34,36,37,31,30,32,33,2d,2c,2e,2f,29,28,2a,2b,25,24,26,27,21,20,22,23] + phi[38]
		10	p[36,37,34,35,32,33,30,31,3a,3b,38,39,3c,3d,3f,3e,26,27,24,25,22,23,20,21,2a,2b,28,29,2c,2d,2f,2e,16,17,14,15,12,13,10,11,1a,1b,18,19,1c,1d,1f,1e,06,07,04,05,02,03,00,01,0a,0b,08,09,0c,0d,0f,0e] + phi[37]
		11	p[18,19,1b,1a,1f,1e,1d,1c,10,11,13,12,17,16,15,14,08,09,0b,0a,0f,0e,0d,0c,00,01,03,02,07,06,05,04,28,29,2b,2a,2f,2e,2d,2c,20,21,23,22,27,26,25,24,30,31,33,32,37,36,35,34,3f,3e,3d,3c,3b,3a,39,38] + phi[36]
		12	p[2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,37,36,35,34,33,32,31,30,3b,3a,39,38,3d,3c,3e,3f,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,17,16,15,14,13,12,11,10,1b,1a,19,18,1d,1c,1e,1f] + phi[35]
		13	p[05,04,06,07,01,00,02,03,09,08,0a,0b,0e,0f,0c,0d,19,18,1a,1b,1e,1f,1c,1d,11,10,12,13,16,17,14,15,39,38,3a,3b,3e,3f,3c,3d,31,30,32,33,36,37,34,35,29,28,2a,2b,2e,2f,2c,2d,21,20,22,23,26,27,24,25] + phi[34]
		14	p[32,33,30,31,34,35,37,36,3c,3d,3f,3e,38,39,3b,3a,22,23,20,21,24,25,27,26,2c,2d,2f,2e,28,29,2b,2a,12,13,10,11,14,15,17,16,1c,1d,1f,1e,18,19,1b,1a,02,03,00,01,04,05,07,06,0c,0d,0f,0e,08,09,0b,0a] + phi[33]
		15	p[1b,1a,19,18,1d,1c,1e,1f,13,12,11,10,15,14,16,17,0b,0a,09,08,0d,0c,0e,0f,03,02,01,00,05,04,06,07,2b,2a,29,28,2d,2c,2e,2f,23,22,21,20,25,24,26,27,33,32,31,30,35,34,36,37,3d,3c,3e,3f,39,38,3a,3b] + phi[32]
		16	p[21,20,22,23,26,27,24,25,2e,2f,2c,2d,2a,2b,28,29,3e,3f,3c,3d,3a,3b,38,39,36,37,34,35,32,33,30,31,01,00,02,03,06,07,04,05,0e,0f,0c,0d,0a,0b,08,09,1e,1f,1c,1d,1a,1b,18,19,16,17,14,15,12,13,10,11] + phi[31]
		17	p[0c,0d,0f,0e,08,09,0b,0a,04,05,07,06,00,01,03,02,14,15,17,16,10,11,13,12,18,19,1b,1a,1f,1e,1d,1c,34,35,37,36,30,31,33,32,38,39,3b,3a,3f,3e,3d,3c,24,25,27,26,20,21,23,22,28,29,2b,2a,2f,2e,2d,2c] + phi[30]
		18	p[3d,3c,3e,3f,39,38,3a,3b,35,34,36,37,31,30,32,33,2d,2c,2e,2f,29,28,2a,2b,25,24,26,27,21,20,22,23,1d,1c,1e,1f,19,18,1a,1b,15,14,16,17,11,10,12,13,0d,0c,0e,0f,09,08,0a,0b,05,04,06,07,01,00,02,03] + phi[2f]
		19	p[16,17,14,15,12,13,10,11,1a,1b,18,19,1c,1d,1f,1e,06,07,04,05,02,03,00,01,0a,0b,08,09,0c,0d,0f,0e,26,27,24,25,22,23,20,21,2a,2b,28,29,2c,2d,2f,2e,3a,3b,38,39,3c,3d,3f,3e,32,33,30,31,34,35,37,36] + phi[2e]
		1a	p[28,29,2b,2a,2f,2e,2d,2c,20,21,23,22,27,26,25,24,30,31,33,32,37,36,35,34,3f,3e,3d,3c,3b,3a,39,38,08,09,0b,0a,0f,0e,0d,0c,00,01,03,02,07,06,05,04,10,11,13,12,17,16,15,14,1f,1e,1d,1c,1b,1a,19,18] + phi[2d]
		1b	p[0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,17,16,15,14,13,12,11,10,1b,1a,19,18,1d,1c,1e,1f,37,36,35,34,33,32,31,30,3b,3a,39,38,3d,3c,3e,3f,27,26,25,24,23,22,21,20,2b,2a,29,28,2d,2c,2e,2f] + phi[2c]
		1c	p[39,38,3a,3b,3e,3f,3c,3d,31,30,32,33,36,37,34,35,29,28,2a,2b,2e,2f,2c,2d,21,20,22,23,26,27,24,25,19,18,1a,1b,1e,1f,1c,1d,11,10,12,13,16,17,14,15,09,08,0a,0b,0e,0f,0c,0d,01,00,02,03,06,07,04,05] + phi[2b]
		1d	p[12,13,10,11,14,15,17,16,1c,1d,1f,1e,18,19,1b,1a,02,03,00,01,04,05,07,06,0c,0d,0f,0e,08,09,0b,0a,22,23,20,21,24,25,27,26,2c,2d,2f,2e,28,29,2b,2a,3c,3d,3f,3e,38,39,3b,3a,34,35,37,36,30,31,33,32] + phi[2a]
		1e	p[2b,2a,29,28,2d,2c,2e,2f,23,22,21,20,25,24,26,27,33,32,31,30,35,34,36,37,3d,3c,3e,3f,39,38,3a,3b,0b,0a,09,08,0d,0c,0e,0f,03,02,01,00,05,04,06,07,13,12,11,10,15,14,16,17,1d,1c,1e,1f,19,18,1a,1b] + phi[29]
		1f	p[01,00,02,03,06,07,04,05,0e,0f,0c,0d,0a,0b,08,09,1e,1f,1c,1d,1a,1b,18,19,16,17,14,15,12,13,10,11,3e,3f,3c,3d,3a,3b,38,39,36,37,34,35,32,33,30,31,2e,2f,2c,2d,2a,2b,28,29,26,27,24,25,22,23,20,21] + phi[28]
		20	p[34,35,37,36,30,31,33,32,38,39,3b,3a,3f,3e,3d,3c,24,25,27,26,20,21,23,22,28,29,2b,2a,2f,2e,2d,2c,14,15,17,16,10,11,13,12,18,19,1b,1a,1f,1e,1d,1c,04,05,07,06,00,01,03,02,08,09,0b,0a,0f,0e,0d,0c] + phi[27]
		21	p[1d,1c,1e,1f,19,18,1a,1b,15,14,16,17,11,10,12,13,0d,0c,0e,0f,09,08,0a,0b,05,04,06,07,01,00,02,03,2d,2c,2e,2f,29,28,2a,2b,25,24,26,27,21,20,22,23,35,34,36,37,31,30,32,33,39,38,3a,3b,3e,3f,3c,3d] + phi[26]
		22	p[26,27,24,25,22,23,20,21,2a,2b,28,29,2c,2d,2f,2e,3a,3b,38,39,3c,3d,3f,3e,32,33,30,31,34,35,37,36,06,07,04,05,02,03,00,01,0a,0b,08,09,0c,0d,0f,0e,1a,1b,18,19,1c,1d,1f,1e,12,13,10,11,14,15,17,16] + phi[25]
		23	p[08,09,0b,0a,0f,0e,0d,0c,00,01,03,02,07,06,05,04,10,11,13,12,17,16,15,14,1f,1e,1d,1c,1b,1a,19,18,30,31,33,32,37,36,35,34,3f,3e,3d,3c,3b,3a,39,38,20,21,23,22,27,26,25,24,2f,2e,2d,2c,2b,2a,29,28] + phi[24]
		24	p[37,36,35,34,33,32,31,30,3b,3a,39,38,3d,3c,3e,3f,27,26,25,24,23,22,21,20,2b,2a,29,28,2d,2c,2e,2f,17,16,15,14,13,12,11,10,1b,1a,19,18,1d,1c,1e,1f,07,06,05,04,03,02,01,00,0b,0a,09,08,0d,0c,0e,0f] + phi[23]
		25	p[19,18,1a,1b,1e,1f,1c,1d,11,10,12,13,16,17,14,15,09,08,0a,0b,0e,0f,0c,0d,01,00,02,03,06,07,04,05,29,28,2a,2b,2e,2f,2c,2d,21,20,22,23,26,27,24,25,31,30,32,33,36,37,34,35,3e,3f,3c,3d,3a,3b,38,39] + phi[22]
		26	p[22,23,20,21,24,25,27,26,2c,2d,2f,2e,28,29,2b,2a,3c,3d,3f,3e,38,39,3b,3a,34,35,37,36,30,31,33,32,02,03,00,01,04,05,07,06,0c,0d,0f,0e,08,09,0b,0a,1c,1d,1f,1e,18,19,1b,1a,14,15,17,16,10,11,13,12] + phi[21]
		27	p[0b,0a,09,08,0d,0c,0e,0f,03,02,01,00,05,04,06,07,13,12,11,10,15,14,16,17,1d,1c,1e,1f,19,18,1a,1b,33,32,31,30,35,34,36,37,3d,3c,3e,3f,39,38,3a,3b,23,22,21,20,25,24,26,27,2d,2c,2e,2f,29,28,2a,2b] + phi[20]
		28	p[3e,3f,3c,3d,3a,3b,38,39,36,37,34,35,32,33,30,31,2e,2f,2c,2d,2a,2b,28,29,26,27,24,25,22,23,20,21,1e,1f,1c,1d,1a,1b,18,19,16,17,14,15,12,13,10,11,0e,0f,0c,0d,0a,0b,08,09,06,07,04,05,02,03,00,01] + phi[1f]
		29	p[14,15,17,16,10,11,13,12,18,19,1b,1a,1f,1e,1d,1c,04,05,07,06,00,01,03,02,08,09,0b,0a,0f,0e,0d,0c,24,25,27,26,20,21,23,22,28,29,2b,2a,2f,2e,2d,2c,38,39,3b,3a,3f,3e,3d,3c,30,31,33,32,37,36,35,34] + phi[1e]
		2a	p[2d,2c,2e,2f,29,28,2a,2b,25,24,26,27,21,20,22,23,35,34,36,37,31,30,32,33,39,38,3a,3b,3e,3f,3c,3d,0d,0c,0e,0f,09,08,0a,0b,05,04,06,07,01,00,02,03,15,14,16,17,11,10,12,13,19,18,1a,1b,1e,1f,1c,1d] + phi[1d]
		2b	p[06,07,04,05,02,03,00,01,0a,0b,08,09,0c,0d,0f,0e,1a,1b,18,19,1c,1d,1f,1e,12,13,10,11,14,15,17,16,3a,3b,38,39,3c,3d,3f,3e,32,33,30,31,34,35,37,36,2a,2b,28,29,2c,2d,2f,2e,22,23,20,21,24,25,27,26] + phi[1c]
		2c	p[30,31,33,32,37,36,35,34,3f,3e,3d,3c,3b,3a,39,38,20,21,23,22,27,26,25,24,2f,2e,2d,2c,2b,2a,29,28,10,11,13,12,17,16,15,14,1f,1e,1d,1c,1b,1a,19,18,00,01,03,02,07,06,05,04,0f,0e,0d,0c,0b,0a,09,08] + phi[1b]
		2d	p[17,16,15,14,13,12,11,10,1b,1a,19,18,1d,1c,1e,1f,07,06,05,04,03,02,01,00,0b,0a,09,08,0d,0c,0e,0f,27,26,25,24,23,22,21,20,2b,2a,29,28,2d,2c,2e,2f,3b,3a,39,38,3d,3c,3e,3f,33,32,31,30,35,34,36,37] + phi[1a]
		2e	p[29,28,2a,2b,2e,2f,2c,2d,21,20,22,23,26,27,24,25,31,30,32,33,36,37,34,35,3e,3f,3c,3d,3a,3b,38,39,09,08,0a,0b,0e,0f,0c,0d,01,00,02,03,06,07,04,05,11,10,12,13,16,17,14,15,1e,1f,1c,1d,1a,1b,18,19] + phi[19]
		2f	p[02,03,00,01,04,05,07,06,0c,0d,0f,0e,08,09,0b,0a,1c,1d,1f,1e,18,19,1b,1a,14,15,17,16,10,11,13,12,3c,3d,3f,3e,38,39,3b,3a,34,35,37,36,30,31,33,32,2c,2d,2f,2e,28,29,2b,2a,24,25,27,26,20,21,23,22] + phi[18]
		30	p[33,32,31,30,35,34,36,37,3d,3c,3e,3f,39,38,3a,3b,23,22,21,20,25,24,26,27,2d,2c,2e,2f,29,28,2a,2b,13,12,11,10,15,14,16,17,1d,1c,1e,1f,19,18,1a,1b,03,02,01,00,05,04,06,07,0d,0c,0e,0f,09,08,0a,0b] + phi[17]
		31	p[1e,1f,1c,1d,1a,1b,18,19,16,17,14,15,12,13,10,11,0e,0f,0c,0d,0a,0b,08,09,06,07,04,05,02,03,00,01,2e,2f,2c,2d,2a,2b,28,29,26,27,24,25,22,23,20,21,36,37,34,35,32,33,30,31,3a,3b,38,39,3c,3d,3f,3e] + phi[16]
		32	p[24,25,27,26,20,21,23,22,28,29,2b,2a,2f,2e,2d,2c,38,39,3b,3a,3f,3e,3d,3c,30,31,33,32,37,36,35,34,04,05,07,06,00,01,03,02,08,09,0b,0a,0f,0e,0d,0c,18,19,1b,1a,1f,1e,1d,1c,10,11,13,12,17,16,15,14] + phi[15]
		33	p[0d,0c,0e,0f,09,08,0a,0b,05,04,06,07,01,00,02,03,15,14,16,17,11,10,12,13,19,18,1a,1b,1e,1f,1c,1d,35,34,36,37,31,30,32,33,39,38,3a,3b,3e,3f,3c,3d,25,24,26,27,21,20,22,23,29,28,2a,2b,2e,2f,2c,2d] + phi[14]
		34	p[3a,3b,38,39,3c,3d,3f,3e,32,33,30,31,34,35,37,36,2a,2b,28,29,2c,2d,2f,2e,22,23,20,21,24,25,27,26,1a,1b,18,19,1c,1d,1f,1e,12,13,10,11,14,15,17,16,0a,0b,08,09,0c,0d,0f,0e,02,03,00,01,04,05,07,06] + phi[13]
		35	p[10,11,13,12,17,16,15,14,1f,1e,1d,1c,1b,1a,19,18,00,01,03,02,07,06,05,04,0f,0e,0d,0c,0b,0a,09,08,20,21,23,22,27,26,25,24,2f,2e,2d,2c,2b,2a,29,28,3f,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30] + phi[12]
		36	p[27,26,25,24,23,22,21,20,2b,2a,29,28,2d,2c,2e,2f,3b,3a,39,38,3d,3c,3e,3f,33,32,31,30,35,34,36,37,07,06,05,04,03,02,01,00,0b,0a,09,08,0d,0c,0e,0f,1b,1a,19,18,1d,1c,1e,1f,13,12,11,10,15,14,16,17] + phi[11]
		37	p[09,08,0a,0b,0e,0f,0c,0d,01,00,02,03,06,07,04,05,11,10,12,13,16,17,14,15,1e,1f,1c,1d,1a,1b,18,19,31,30,32,33,36,37,34,35,3e,3f,3c,3d,3a,3b,38,39,21,20,22,23,26,27,24,25,2e,2f,2c,2d,2a,2b,28,29] + phi[10]
		38	p[3c,3d,3f,3e,38,39,3b,3a,34,35,37,36,30,31,33,32,2c,2d,2f,2e,28,29,2b,2a,24,25,27,26,20,21,23,22,1c,1d,1f,1e,18,19,1b,1a,14,15,17,16,10,11,13,12,0c,0d,0f,0e,08,09,0b,0a,04,05,07,06,00,01,03,02] + phi[0f]
		39	p[13,12,11,10,15,14,16,17,1d,1c,1e,1f,19,18,1a,1b,03,02,01,00,05,04,06,07,0d,0c,0e,0f,09,08,0a,0b,23,22,21,20,25,24,26,27,2d,2c,2e,2f,29,28,2a,2b,3d,3c,3e,3f,39,38,3a,3b,35,34,36,37,31,30,32,33] + phi[0e]
		3a	p[2e,2f,2c,2d,2a,2b,28,29,26,27,24,25,22,23,20,21,36,37,34,35,32,33,30,31,3a,3b,38,39,3c,3d,3f,3e,0e,0f,0c,0d,0a,0b,08,09,06,07,04,05,02,03,00,01,16,17,14,15,12,13,10,11,1a,1b,18,19,1c,1d,1f,1e] + phi[0d]
		3b	p[04,05,07,06,00,01,03,02,08,09,0b,0a,0f,0e,0d,0c,18,19,1b,1a,1f,1e,1d,1c,10,11,13,12,17,16,15,14,38,39,3b,3a,3f,3e,3d,3c,30,31,33,32,37,36,35,34,28,29,2b,2a,2f,2e,2d,2c,20,21,23,22,27,26,25,24] + phi[0c]
		3c	p[35,34,36,37,31,30,32,33,39,38,3a,3b,3e,3f,3c,3d,25,24,26,27,21,20,22,23,29,28,2a,2b,2e,2f,2c,2d,15,14,16,17,11,10,12,13,19,18,1a,1b,1e,1f,1c,1d,05,04,06,07,01,00,02,03,09,08,0a,0b,0e,0f,0c,0d] + phi[0b]
		3d	p[1a,1b,18,19,1c,1d,1f,1e,12,13,10,11,14,15,17,16,0a,0b,08,09,0c,0d,0f,0e,02,03,00,01,04,05,07,06,2a,2b,28,29,2c,2d,2f,2e,22,23,20,21,24,25,27,26,32,33,30,31,34,35,37,36,3c,3d,3f,3e,38,39,3b,3a] + phi[0a]
		3e	p[20,21,23,22,27,26,25,24,2f,2e,2d,2c,2b,2a,29,28,3f,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,00,01,03,02,07,06,05,04,0f,0e,0d,0c,0b,0a,09,08,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10] + phi[09]
	*/
	/*...gather the needed data (4032 64-bit complex, i.e. 2016 64-bit reals) and do 63 radix-64 transforms...*/
		tptr = t;
		for(k = 0; k < ODD_RADIX; ++k) {
		//	iptr = dit64_iidx_lo + (k<<6);	Original code using full-length dit64_iidx_lo[] array
			iptr = dif64_oidx_lo + (k<<6);
			for(l = 0; l < 64; l++)
			{
			//	io_offsets[l] = plo[*iptr++];	Original code using full-length dit64_iidx_lo[] array
				io_offsets[l] = plo[*(iptr+dit64_iidx_lo[l])];
			}
			jt = j1 + phi[dit_phi[k]];
			RADIX_64_DIT(a+jt,io_offsets,RE_IM_STRIDE, (double *)tptr,t_offsets,1);
			tptr++;
		}
	/*...and now do 64 radix-63 transforms. The required output permutation is as follows:
	int operm[RADIX] = {
		000,f80,f40,f00,ec0,e80,e40,e00,dc0,d80,d40,d00,cc0,c80,c40,c00,bc0,b80,b40,b00,ac0,a80,a40,a00,9c0,980,940,900,8c0,880,840,800,7c0,780,740,700,6c0,680,640,600,5c0,580,540,500,4c0,480,440,400,3c0,380,340,300,2c0,280,240,200,1c0,180,140,100,0c0,080,040,
		03f,fbf,f7f,f3f,eff,ebf,e7f,e3f,dff,dbf,d7f,d3f,cff,cbf,c7f,c3f,bff,bbf,b7f,b3f,aff,abf,a7f,a3f,9ff,9bf,97f,93f,8ff,8bf,87f,83f,7ff,7bf,77f,73f,6ff,6bf,67f,63f,5ff,5bf,57f,53f,4ff,4bf,47f,43f,3ff,3bf,37f,33f,2ff,2bf,27f,23f,1ff,1bf,17f,13f,0ff,0bf,07f,
		07e,03e,fbe,f7e,f3e,efe,ebe,e7e,e3e,dfe,dbe,d7e,d3e,cfe,cbe,c7e,c3e,bfe,bbe,b7e,b3e,afe,abe,a7e,a3e,9fe,9be,97e,93e,8fe,8be,87e,83e,7fe,7be,77e,73e,6fe,6be,67e,63e,5fe,5be,57e,53e,4fe,4be,47e,43e,3fe,3be,37e,33e,2fe,2be,27e,23e,1fe,1be,17e,13e,0fe,0be,
		0bd,07d,03d,fbd,f7d,f3d,efd,ebd,e7d,e3d,dfd,dbd,d7d,d3d,cfd,cbd,c7d,c3d,bfd,bbd,b7d,b3d,afd,abd,a7d,a3d,9fd,9bd,97d,93d,8fd,8bd,87d,83d,7fd,7bd,77d,73d,6fd,6bd,67d,63d,5fd,5bd,57d,53d,4fd,4bd,47d,43d,3fd,3bd,37d,33d,2fd,2bd,27d,23d,1fd,1bd,17d,13d,0fd,
		0fc,0bc,07c,03c,fbc,f7c,f3c,efc,ebc,e7c,e3c,dfc,dbc,d7c,d3c,cfc,cbc,c7c,c3c,bfc,bbc,b7c,b3c,afc,abc,a7c,a3c,9fc,9bc,97c,93c,8fc,8bc,87c,83c,7fc,7bc,77c,73c,6fc,6bc,67c,63c,5fc,5bc,57c,53c,4fc,4bc,47c,43c,3fc,3bc,37c,33c,2fc,2bc,27c,23c,1fc,1bc,17c,13c,
		13b,0fb,0bb,07b,03b,fbb,f7b,f3b,efb,ebb,e7b,e3b,dfb,dbb,d7b,d3b,cfb,cbb,c7b,c3b,bfb,bbb,b7b,b3b,afb,abb,a7b,a3b,9fb,9bb,97b,93b,8fb,8bb,87b,83b,7fb,7bb,77b,73b,6fb,6bb,67b,63b,5fb,5bb,57b,53b,4fb,4bb,47b,43b,3fb,3bb,37b,33b,2fb,2bb,27b,23b,1fb,1bb,17b,
		17a,13a,0fa,0ba,07a,03a,fba,f7a,f3a,efa,eba,e7a,e3a,dfa,dba,d7a,d3a,cfa,cba,c7a,c3a,bfa,bba,b7a,b3a,afa,aba,a7a,a3a,9fa,9ba,97a,93a,8fa,8ba,87a,83a,7fa,7ba,77a,73a,6fa,6ba,67a,63a,5fa,5ba,57a,53a,4fa,4ba,47a,43a,3fa,3ba,37a,33a,2fa,2ba,27a,23a,1fa,1ba,
		1b9,179,139,0f9,0b9,079,039,fb9,f79,f39,ef9,eb9,e79,e39,df9,db9,d79,d39,cf9,cb9,c79,c39,bf9,bb9,b79,b39,af9,ab9,a79,a39,9f9,9b9,979,939,8f9,8b9,879,839,7f9,7b9,779,739,6f9,6b9,679,639,5f9,5b9,579,539,4f9,4b9,479,439,3f9,3b9,379,339,2f9,2b9,279,239,1f9,
		1f8,1b8,178,138,0f8,0b8,078,038,fb8,f78,f38,ef8,eb8,e78,e38,df8,db8,d78,d38,cf8,cb8,c78,c38,bf8,bb8,b78,b38,af8,ab8,a78,a38,9f8,9b8,978,938,8f8,8b8,878,838,7f8,7b8,778,738,6f8,6b8,678,638,5f8,5b8,578,538,4f8,4b8,478,438,3f8,3b8,378,338,2f8,2b8,278,238,
		237,1f7,1b7,177,137,0f7,0b7,077,037,fb7,f77,f37,ef7,eb7,e77,e37,df7,db7,d77,d37,cf7,cb7,c77,c37,bf7,bb7,b77,b37,af7,ab7,a77,a37,9f7,9b7,977,937,8f7,8b7,877,837,7f7,7b7,777,737,6f7,6b7,677,637,5f7,5b7,577,537,4f7,4b7,477,437,3f7,3b7,377,337,2f7,2b7,277,
		276,236,1f6,1b6,176,136,0f6,0b6,076,036,fb6,f76,f36,ef6,eb6,e76,e36,df6,db6,d76,d36,cf6,cb6,c76,c36,bf6,bb6,b76,b36,af6,ab6,a76,a36,9f6,9b6,976,936,8f6,8b6,876,836,7f6,7b6,776,736,6f6,6b6,676,636,5f6,5b6,576,536,4f6,4b6,476,436,3f6,3b6,376,336,2f6,2b6,
		2b5,275,235,1f5,1b5,175,135,0f5,0b5,075,035,fb5,f75,f35,ef5,eb5,e75,e35,df5,db5,d75,d35,cf5,cb5,c75,c35,bf5,bb5,b75,b35,af5,ab5,a75,a35,9f5,9b5,975,935,8f5,8b5,875,835,7f5,7b5,775,735,6f5,6b5,675,635,5f5,5b5,575,535,4f5,4b5,475,435,3f5,3b5,375,335,2f5,
		2f4,2b4,274,234,1f4,1b4,174,134,0f4,0b4,074,034,fb4,f74,f34,ef4,eb4,e74,e34,df4,db4,d74,d34,cf4,cb4,c74,c34,bf4,bb4,b74,b34,af4,ab4,a74,a34,9f4,9b4,974,934,8f4,8b4,874,834,7f4,7b4,774,734,6f4,6b4,674,634,5f4,5b4,574,534,4f4,4b4,474,434,3f4,3b4,374,334,
		333,2f3,2b3,273,233,1f3,1b3,173,133,0f3,0b3,073,033,fb3,f73,f33,ef3,eb3,e73,e33,df3,db3,d73,d33,cf3,cb3,c73,c33,bf3,bb3,b73,b33,af3,ab3,a73,a33,9f3,9b3,973,933,8f3,8b3,873,833,7f3,7b3,773,733,6f3,6b3,673,633,5f3,5b3,573,533,4f3,4b3,473,433,3f3,3b3,373,
		372,332,2f2,2b2,272,232,1f2,1b2,172,132,0f2,0b2,072,032,fb2,f72,f32,ef2,eb2,e72,e32,df2,db2,d72,d32,cf2,cb2,c72,c32,bf2,bb2,b72,b32,af2,ab2,a72,a32,9f2,9b2,972,932,8f2,8b2,872,832,7f2,7b2,772,732,6f2,6b2,672,632,5f2,5b2,572,532,4f2,4b2,472,432,3f2,3b2,
		3b1,371,331,2f1,2b1,271,231,1f1,1b1,171,131,0f1,0b1,071,031,fb1,f71,f31,ef1,eb1,e71,e31,df1,db1,d71,d31,cf1,cb1,c71,c31,bf1,bb1,b71,b31,af1,ab1,a71,a31,9f1,9b1,971,931,8f1,8b1,871,831,7f1,7b1,771,731,6f1,6b1,671,631,5f1,5b1,571,531,4f1,4b1,471,431,3f1,
		3f0,3b0,370,330,2f0,2b0,270,230,1f0,1b0,170,130,0f0,0b0,070,030,fb0,f70,f30,ef0,eb0,e70,e30,df0,db0,d70,d30,cf0,cb0,c70,c30,bf0,bb0,b70,b30,af0,ab0,a70,a30,9f0,9b0,970,930,8f0,8b0,870,830,7f0,7b0,770,730,6f0,6b0,670,630,5f0,5b0,570,530,4f0,4b0,470,430,
		42f,3ef,3af,36f,32f,2ef,2af,26f,22f,1ef,1af,16f,12f,0ef,0af,06f,02f,faf,f6f,f2f,eef,eaf,e6f,e2f,def,daf,d6f,d2f,cef,caf,c6f,c2f,bef,baf,b6f,b2f,aef,aaf,a6f,a2f,9ef,9af,96f,92f,8ef,8af,86f,82f,7ef,7af,76f,72f,6ef,6af,66f,62f,5ef,5af,56f,52f,4ef,4af,46f,
		46e,42e,3ee,3ae,36e,32e,2ee,2ae,26e,22e,1ee,1ae,16e,12e,0ee,0ae,06e,02e,fae,f6e,f2e,eee,eae,e6e,e2e,dee,dae,d6e,d2e,cee,cae,c6e,c2e,bee,bae,b6e,b2e,aee,aae,a6e,a2e,9ee,9ae,96e,92e,8ee,8ae,86e,82e,7ee,7ae,76e,72e,6ee,6ae,66e,62e,5ee,5ae,56e,52e,4ee,4ae,
		4ad,46d,42d,3ed,3ad,36d,32d,2ed,2ad,26d,22d,1ed,1ad,16d,12d,0ed,0ad,06d,02d,fad,f6d,f2d,eed,ead,e6d,e2d,ded,dad,d6d,d2d,ced,cad,c6d,c2d,bed,bad,b6d,b2d,aed,aad,a6d,a2d,9ed,9ad,96d,92d,8ed,8ad,86d,82d,7ed,7ad,76d,72d,6ed,6ad,66d,62d,5ed,5ad,56d,52d,4ed,
		4ec,4ac,46c,42c,3ec,3ac,36c,32c,2ec,2ac,26c,22c,1ec,1ac,16c,12c,0ec,0ac,06c,02c,fac,f6c,f2c,eec,eac,e6c,e2c,dec,dac,d6c,d2c,cec,cac,c6c,c2c,bec,bac,b6c,b2c,aec,aac,a6c,a2c,9ec,9ac,96c,92c,8ec,8ac,86c,82c,7ec,7ac,76c,72c,6ec,6ac,66c,62c,5ec,5ac,56c,52c,
		52b,4eb,4ab,46b,42b,3eb,3ab,36b,32b,2eb,2ab,26b,22b,1eb,1ab,16b,12b,0eb,0ab,06b,02b,fab,f6b,f2b,eeb,eab,e6b,e2b,deb,dab,d6b,d2b,ceb,cab,c6b,c2b,beb,bab,b6b,b2b,aeb,aab,a6b,a2b,9eb,9ab,96b,92b,8eb,8ab,86b,82b,7eb,7ab,76b,72b,6eb,6ab,66b,62b,5eb,5ab,56b,
		56a,52a,4ea,4aa,46a,42a,3ea,3aa,36a,32a,2ea,2aa,26a,22a,1ea,1aa,16a,12a,0ea,0aa,06a,02a,faa,f6a,f2a,eea,eaa,e6a,e2a,dea,daa,d6a,d2a,cea,caa,c6a,c2a,bea,baa,b6a,b2a,aea,aaa,a6a,a2a,9ea,9aa,96a,92a,8ea,8aa,86a,82a,7ea,7aa,76a,72a,6ea,6aa,66a,62a,5ea,5aa,
		5a9,569,529,4e9,4a9,469,429,3e9,3a9,369,329,2e9,2a9,269,229,1e9,1a9,169,129,0e9,0a9,069,029,fa9,f69,f29,ee9,ea9,e69,e29,de9,da9,d69,d29,ce9,ca9,c69,c29,be9,ba9,b69,b29,ae9,aa9,a69,a29,9e9,9a9,969,929,8e9,8a9,869,829,7e9,7a9,769,729,6e9,6a9,669,629,5e9,
		5e8,5a8,568,528,4e8,4a8,468,428,3e8,3a8,368,328,2e8,2a8,268,228,1e8,1a8,168,128,0e8,0a8,068,028,fa8,f68,f28,ee8,ea8,e68,e28,de8,da8,d68,d28,ce8,ca8,c68,c28,be8,ba8,b68,b28,ae8,aa8,a68,a28,9e8,9a8,968,928,8e8,8a8,868,828,7e8,7a8,768,728,6e8,6a8,668,628,
		627,5e7,5a7,567,527,4e7,4a7,467,427,3e7,3a7,367,327,2e7,2a7,267,227,1e7,1a7,167,127,0e7,0a7,067,027,fa7,f67,f27,ee7,ea7,e67,e27,de7,da7,d67,d27,ce7,ca7,c67,c27,be7,ba7,b67,b27,ae7,aa7,a67,a27,9e7,9a7,967,927,8e7,8a7,867,827,7e7,7a7,767,727,6e7,6a7,667,
		666,626,5e6,5a6,566,526,4e6,4a6,466,426,3e6,3a6,366,326,2e6,2a6,266,226,1e6,1a6,166,126,0e6,0a6,066,026,fa6,f66,f26,ee6,ea6,e66,e26,de6,da6,d66,d26,ce6,ca6,c66,c26,be6,ba6,b66,b26,ae6,aa6,a66,a26,9e6,9a6,966,926,8e6,8a6,866,826,7e6,7a6,766,726,6e6,6a6,
		6a5,665,625,5e5,5a5,565,525,4e5,4a5,465,425,3e5,3a5,365,325,2e5,2a5,265,225,1e5,1a5,165,125,0e5,0a5,065,025,fa5,f65,f25,ee5,ea5,e65,e25,de5,da5,d65,d25,ce5,ca5,c65,c25,be5,ba5,b65,b25,ae5,aa5,a65,a25,9e5,9a5,965,925,8e5,8a5,865,825,7e5,7a5,765,725,6e5,
		6e4,6a4,664,624,5e4,5a4,564,524,4e4,4a4,464,424,3e4,3a4,364,324,2e4,2a4,264,224,1e4,1a4,164,124,0e4,0a4,064,024,fa4,f64,f24,ee4,ea4,e64,e24,de4,da4,d64,d24,ce4,ca4,c64,c24,be4,ba4,b64,b24,ae4,aa4,a64,a24,9e4,9a4,964,924,8e4,8a4,864,824,7e4,7a4,764,724,
		723,6e3,6a3,663,623,5e3,5a3,563,523,4e3,4a3,463,423,3e3,3a3,363,323,2e3,2a3,263,223,1e3,1a3,163,123,0e3,0a3,063,023,fa3,f63,f23,ee3,ea3,e63,e23,de3,da3,d63,d23,ce3,ca3,c63,c23,be3,ba3,b63,b23,ae3,aa3,a63,a23,9e3,9a3,963,923,8e3,8a3,863,823,7e3,7a3,763,
		762,722,6e2,6a2,662,622,5e2,5a2,562,522,4e2,4a2,462,422,3e2,3a2,362,322,2e2,2a2,262,222,1e2,1a2,162,122,0e2,0a2,062,022,fa2,f62,f22,ee2,ea2,e62,e22,de2,da2,d62,d22,ce2,ca2,c62,c22,be2,ba2,b62,b22,ae2,aa2,a62,a22,9e2,9a2,962,922,8e2,8a2,862,822,7e2,7a2,
		7a1,761,721,6e1,6a1,661,621,5e1,5a1,561,521,4e1,4a1,461,421,3e1,3a1,361,321,2e1,2a1,261,221,1e1,1a1,161,121,0e1,0a1,061,021,fa1,f61,f21,ee1,ea1,e61,e21,de1,da1,d61,d21,ce1,ca1,c61,c21,be1,ba1,b61,b21,ae1,aa1,a61,a21,9e1,9a1,961,921,8e1,8a1,861,821,7e1,
		7e0,7a0,760,720,6e0,6a0,660,620,5e0,5a0,560,520,4e0,4a0,460,420,3e0,3a0,360,320,2e0,2a0,260,220,1e0,1a0,160,120,0e0,0a0,060,020,fa0,f60,f20,ee0,ea0,e60,e20,de0,da0,d60,d20,ce0,ca0,c60,c20,be0,ba0,b60,b20,ae0,aa0,a60,a20,9e0,9a0,960,920,8e0,8a0,860,820,
		81f,7df,79f,75f,71f,6df,69f,65f,61f,5df,59f,55f,51f,4df,49f,45f,41f,3df,39f,35f,31f,2df,29f,25f,21f,1df,19f,15f,11f,0df,09f,05f,01f,f9f,f5f,f1f,edf,e9f,e5f,e1f,ddf,d9f,d5f,d1f,cdf,c9f,c5f,c1f,bdf,b9f,b5f,b1f,adf,a9f,a5f,a1f,9df,99f,95f,91f,8df,89f,85f,
		85e,81e,7de,79e,75e,71e,6de,69e,65e,61e,5de,59e,55e,51e,4de,49e,45e,41e,3de,39e,35e,31e,2de,29e,25e,21e,1de,19e,15e,11e,0de,09e,05e,01e,f9e,f5e,f1e,ede,e9e,e5e,e1e,dde,d9e,d5e,d1e,cde,c9e,c5e,c1e,bde,b9e,b5e,b1e,ade,a9e,a5e,a1e,9de,99e,95e,91e,8de,89e,
		89d,85d,81d,7dd,79d,75d,71d,6dd,69d,65d,61d,5dd,59d,55d,51d,4dd,49d,45d,41d,3dd,39d,35d,31d,2dd,29d,25d,21d,1dd,19d,15d,11d,0dd,09d,05d,01d,f9d,f5d,f1d,edd,e9d,e5d,e1d,ddd,d9d,d5d,d1d,cdd,c9d,c5d,c1d,bdd,b9d,b5d,b1d,add,a9d,a5d,a1d,9dd,99d,95d,91d,8dd,
		8dc,89c,85c,81c,7dc,79c,75c,71c,6dc,69c,65c,61c,5dc,59c,55c,51c,4dc,49c,45c,41c,3dc,39c,35c,31c,2dc,29c,25c,21c,1dc,19c,15c,11c,0dc,09c,05c,01c,f9c,f5c,f1c,edc,e9c,e5c,e1c,ddc,d9c,d5c,d1c,cdc,c9c,c5c,c1c,bdc,b9c,b5c,b1c,adc,a9c,a5c,a1c,9dc,99c,95c,91c,
		91b,8db,89b,85b,81b,7db,79b,75b,71b,6db,69b,65b,61b,5db,59b,55b,51b,4db,49b,45b,41b,3db,39b,35b,31b,2db,29b,25b,21b,1db,19b,15b,11b,0db,09b,05b,01b,f9b,f5b,f1b,edb,e9b,e5b,e1b,ddb,d9b,d5b,d1b,cdb,c9b,c5b,c1b,bdb,b9b,b5b,b1b,adb,a9b,a5b,a1b,9db,99b,95b,
		95a,91a,8da,89a,85a,81a,7da,79a,75a,71a,6da,69a,65a,61a,5da,59a,55a,51a,4da,49a,45a,41a,3da,39a,35a,31a,2da,29a,25a,21a,1da,19a,15a,11a,0da,09a,05a,01a,f9a,f5a,f1a,eda,e9a,e5a,e1a,dda,d9a,d5a,d1a,cda,c9a,c5a,c1a,bda,b9a,b5a,b1a,ada,a9a,a5a,a1a,9da,99a,
		999,959,919,8d9,899,859,819,7d9,799,759,719,6d9,699,659,619,5d9,599,559,519,4d9,499,459,419,3d9,399,359,319,2d9,299,259,219,1d9,199,159,119,0d9,099,059,019,f99,f59,f19,ed9,e99,e59,e19,dd9,d99,d59,d19,cd9,c99,c59,c19,bd9,b99,b59,b19,ad9,a99,a59,a19,9d9,
		9d8,998,958,918,8d8,898,858,818,7d8,798,758,718,6d8,698,658,618,5d8,598,558,518,4d8,498,458,418,3d8,398,358,318,2d8,298,258,218,1d8,198,158,118,0d8,098,058,018,f98,f58,f18,ed8,e98,e58,e18,dd8,d98,d58,d18,cd8,c98,c58,c18,bd8,b98,b58,b18,ad8,a98,a58,a18,
		a17,9d7,997,957,917,8d7,897,857,817,7d7,797,757,717,6d7,697,657,617,5d7,597,557,517,4d7,497,457,417,3d7,397,357,317,2d7,297,257,217,1d7,197,157,117,0d7,097,057,017,f97,f57,f17,ed7,e97,e57,e17,dd7,d97,d57,d17,cd7,c97,c57,c17,bd7,b97,b57,b17,ad7,a97,a57,
		a56,a16,9d6,996,956,916,8d6,896,856,816,7d6,796,756,716,6d6,696,656,616,5d6,596,556,516,4d6,496,456,416,3d6,396,356,316,2d6,296,256,216,1d6,196,156,116,0d6,096,056,016,f96,f56,f16,ed6,e96,e56,e16,dd6,d96,d56,d16,cd6,c96,c56,c16,bd6,b96,b56,b16,ad6,a96,
		a95,a55,a15,9d5,995,955,915,8d5,895,855,815,7d5,795,755,715,6d5,695,655,615,5d5,595,555,515,4d5,495,455,415,3d5,395,355,315,2d5,295,255,215,1d5,195,155,115,0d5,095,055,015,f95,f55,f15,ed5,e95,e55,e15,dd5,d95,d55,d15,cd5,c95,c55,c15,bd5,b95,b55,b15,ad5,
		ad4,a94,a54,a14,9d4,994,954,914,8d4,894,854,814,7d4,794,754,714,6d4,694,654,614,5d4,594,554,514,4d4,494,454,414,3d4,394,354,314,2d4,294,254,214,1d4,194,154,114,0d4,094,054,014,f94,f54,f14,ed4,e94,e54,e14,dd4,d94,d54,d14,cd4,c94,c54,c14,bd4,b94,b54,b14,
		b13,ad3,a93,a53,a13,9d3,993,953,913,8d3,893,853,813,7d3,793,753,713,6d3,693,653,613,5d3,593,553,513,4d3,493,453,413,3d3,393,353,313,2d3,293,253,213,1d3,193,153,113,0d3,093,053,013,f93,f53,f13,ed3,e93,e53,e13,dd3,d93,d53,d13,cd3,c93,c53,c13,bd3,b93,b53,
		b52,b12,ad2,a92,a52,a12,9d2,992,952,912,8d2,892,852,812,7d2,792,752,712,6d2,692,652,612,5d2,592,552,512,4d2,492,452,412,3d2,392,352,312,2d2,292,252,212,1d2,192,152,112,0d2,092,052,012,f92,f52,f12,ed2,e92,e52,e12,dd2,d92,d52,d12,cd2,c92,c52,c12,bd2,b92,
		b91,b51,b11,ad1,a91,a51,a11,9d1,991,951,911,8d1,891,851,811,7d1,791,751,711,6d1,691,651,611,5d1,591,551,511,4d1,491,451,411,3d1,391,351,311,2d1,291,251,211,1d1,191,151,111,0d1,091,051,011,f91,f51,f11,ed1,e91,e51,e11,dd1,d91,d51,d11,cd1,c91,c51,c11,bd1,
		bd0,b90,b50,b10,ad0,a90,a50,a10,9d0,990,950,910,8d0,890,850,810,7d0,790,750,710,6d0,690,650,610,5d0,590,550,510,4d0,490,450,410,3d0,390,350,310,2d0,290,250,210,1d0,190,150,110,0d0,090,050,010,f90,f50,f10,ed0,e90,e50,e10,dd0,d90,d50,d10,cd0,c90,c50,c10,
		c0f,bcf,b8f,b4f,b0f,acf,a8f,a4f,a0f,9cf,98f,94f,90f,8cf,88f,84f,80f,7cf,78f,74f,70f,6cf,68f,64f,60f,5cf,58f,54f,50f,4cf,48f,44f,40f,3cf,38f,34f,30f,2cf,28f,24f,20f,1cf,18f,14f,10f,0cf,08f,04f,00f,f8f,f4f,f0f,ecf,e8f,e4f,e0f,dcf,d8f,d4f,d0f,ccf,c8f,c4f,
		c4e,c0e,bce,b8e,b4e,b0e,ace,a8e,a4e,a0e,9ce,98e,94e,90e,8ce,88e,84e,80e,7ce,78e,74e,70e,6ce,68e,64e,60e,5ce,58e,54e,50e,4ce,48e,44e,40e,3ce,38e,34e,30e,2ce,28e,24e,20e,1ce,18e,14e,10e,0ce,08e,04e,00e,f8e,f4e,f0e,ece,e8e,e4e,e0e,dce,d8e,d4e,d0e,cce,c8e,
		c8d,c4d,c0d,bcd,b8d,b4d,b0d,acd,a8d,a4d,a0d,9cd,98d,94d,90d,8cd,88d,84d,80d,7cd,78d,74d,70d,6cd,68d,64d,60d,5cd,58d,54d,50d,4cd,48d,44d,40d,3cd,38d,34d,30d,2cd,28d,24d,20d,1cd,18d,14d,10d,0cd,08d,04d,00d,f8d,f4d,f0d,ecd,e8d,e4d,e0d,dcd,d8d,d4d,d0d,ccd,
		ccc,c8c,c4c,c0c,bcc,b8c,b4c,b0c,acc,a8c,a4c,a0c,9cc,98c,94c,90c,8cc,88c,84c,80c,7cc,78c,74c,70c,6cc,68c,64c,60c,5cc,58c,54c,50c,4cc,48c,44c,40c,3cc,38c,34c,30c,2cc,28c,24c,20c,1cc,18c,14c,10c,0cc,08c,04c,00c,f8c,f4c,f0c,ecc,e8c,e4c,e0c,dcc,d8c,d4c,d0c,
		d0b,ccb,c8b,c4b,c0b,bcb,b8b,b4b,b0b,acb,a8b,a4b,a0b,9cb,98b,94b,90b,8cb,88b,84b,80b,7cb,78b,74b,70b,6cb,68b,64b,60b,5cb,58b,54b,50b,4cb,48b,44b,40b,3cb,38b,34b,30b,2cb,28b,24b,20b,1cb,18b,14b,10b,0cb,08b,04b,00b,f8b,f4b,f0b,ecb,e8b,e4b,e0b,dcb,d8b,d4b,
		d4a,d0a,cca,c8a,c4a,c0a,bca,b8a,b4a,b0a,aca,a8a,a4a,a0a,9ca,98a,94a,90a,8ca,88a,84a,80a,7ca,78a,74a,70a,6ca,68a,64a,60a,5ca,58a,54a,50a,4ca,48a,44a,40a,3ca,38a,34a,30a,2ca,28a,24a,20a,1ca,18a,14a,10a,0ca,08a,04a,00a,f8a,f4a,f0a,eca,e8a,e4a,e0a,dca,d8a,
		d89,d49,d09,cc9,c89,c49,c09,bc9,b89,b49,b09,ac9,a89,a49,a09,9c9,989,949,909,8c9,889,849,809,7c9,789,749,709,6c9,689,649,609,5c9,589,549,509,4c9,489,449,409,3c9,389,349,309,2c9,289,249,209,1c9,189,149,109,0c9,089,049,009,f89,f49,f09,ec9,e89,e49,e09,dc9,
		dc8,d88,d48,d08,cc8,c88,c48,c08,bc8,b88,b48,b08,ac8,a88,a48,a08,9c8,988,948,908,8c8,888,848,808,7c8,788,748,708,6c8,688,648,608,5c8,588,548,508,4c8,488,448,408,3c8,388,348,308,2c8,288,248,208,1c8,188,148,108,0c8,088,048,008,f88,f48,f08,ec8,e88,e48,e08,
		e07,dc7,d87,d47,d07,cc7,c87,c47,c07,bc7,b87,b47,b07,ac7,a87,a47,a07,9c7,987,947,907,8c7,887,847,807,7c7,787,747,707,6c7,687,647,607,5c7,587,547,507,4c7,487,447,407,3c7,387,347,307,2c7,287,247,207,1c7,187,147,107,0c7,087,047,007,f87,f47,f07,ec7,e87,e47,
		e46,e06,dc6,d86,d46,d06,cc6,c86,c46,c06,bc6,b86,b46,b06,ac6,a86,a46,a06,9c6,986,946,906,8c6,886,846,806,7c6,786,746,706,6c6,686,646,606,5c6,586,546,506,4c6,486,446,406,3c6,386,346,306,2c6,286,246,206,1c6,186,146,106,0c6,086,046,006,f86,f46,f06,ec6,e86,
		e85,e45,e05,dc5,d85,d45,d05,cc5,c85,c45,c05,bc5,b85,b45,b05,ac5,a85,a45,a05,9c5,985,945,905,8c5,885,845,805,7c5,785,745,705,6c5,685,645,605,5c5,585,545,505,4c5,485,445,405,3c5,385,345,305,2c5,285,245,205,1c5,185,145,105,0c5,085,045,005,f85,f45,f05,ec5,
		ec4,e84,e44,e04,dc4,d84,d44,d04,cc4,c84,c44,c04,bc4,b84,b44,b04,ac4,a84,a44,a04,9c4,984,944,904,8c4,884,844,804,7c4,784,744,704,6c4,684,644,604,5c4,584,544,504,4c4,484,444,404,3c4,384,344,304,2c4,284,244,204,1c4,184,144,104,0c4,084,044,004,f84,f44,f04,
		f03,ec3,e83,e43,e03,dc3,d83,d43,d03,cc3,c83,c43,c03,bc3,b83,b43,b03,ac3,a83,a43,a03,9c3,983,943,903,8c3,883,843,803,7c3,783,743,703,6c3,683,643,603,5c3,583,543,503,4c3,483,443,403,3c3,383,343,303,2c3,283,243,203,1c3,183,143,103,0c3,083,043,003,f83,f43,
		f42,f02,ec2,e82,e42,e02,dc2,d82,d42,d02,cc2,c82,c42,c02,bc2,b82,b42,b02,ac2,a82,a42,a02,9c2,982,942,902,8c2,882,842,802,7c2,782,742,702,6c2,682,642,602,5c2,582,542,502,4c2,482,442,402,3c2,382,342,302,2c2,282,242,202,1c2,182,142,102,0c2,082,042,002,f82,
		f81,f41,f01,ec1,e81,e41,e01,dc1,d81,d41,d01,cc1,c81,c41,c01,bc1,b81,b41,b01,ac1,a81,a41,a01,9c1,981,941,901,8c1,881,841,801,7c1,781,741,701,6c1,681,641,601,5c1,581,541,501,4c1,481,441,401,3c1,381,341,301,2c1,281,241,201,1c1,181,141,101,0c1,081,041,001
	};
	Separating the above into (mod 64) - indices into plo[] array - and (/= 64) - indices into phi[] (multiples of p40) array -
		row	[00,3e,...,01]-perm lshift count
		00	00	phi[00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01] + plo[00]
		01	00	phi[00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01] + plo[3f]
		02	3e	phi[01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02] + plo[3e]
		03	3d	phi[02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03] + plo[3d]
		04	3c	phi[03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04] + plo[3c]
		05	3b	phi[04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05] + plo[3b]
		06	3a	phi[05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06] + plo[3a]
		07	39	phi[06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07] + plo[39]
		08	38	phi[07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08] + plo[38]
		09	37	phi[08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09] + plo[37]
		0a	36	phi[09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a] + plo[36]
		0b	35	phi[0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b] + plo[35]
		0c	34	phi[0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c] + plo[34]
		0d	33	phi[0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d] + plo[33]
		0e	32	phi[0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e] + plo[32]
		0f	31	phi[0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f] + plo[31]
		10	30	phi[0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10] + plo[30]
		11	2f	phi[10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11] + plo[2f]
		12	2e	phi[11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12] + plo[2e]
		13	2d	phi[12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13] + plo[2d]
		14	2c	phi[13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14] + plo[2c]
		15	2b	phi[14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15] + plo[2b]
		16	2a	phi[15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16] + plo[2a]
		17	29	phi[16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17] + plo[29]
		18	28	phi[17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18] + plo[28]
		19	27	phi[18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19] + plo[27]
		1a	26	phi[19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a] + plo[26]
		1b	25	phi[1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b] + plo[25]
		1c	24	phi[1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c] + plo[24]
		1d	23	phi[1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d] + plo[23]
		1e	22	phi[1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e] + plo[22]
		1f	21	phi[1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f] + plo[21]
		20	20	phi[1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20] + plo[20]
		21	1f	phi[20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21] + plo[1f]
		22	1e	phi[21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22] + plo[1e]
		23	1d	phi[22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23] + plo[1d]
		24	1c	phi[23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24] + plo[1c]
		25	1b	phi[24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25] + plo[1b]
		26	1a	phi[25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26] + plo[1a]
		27	19	phi[26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27] + plo[19]
		28	18	phi[27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28] + plo[18]
		29	17	phi[28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29] + plo[17]
		2a	16	phi[29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a] + plo[16]
		2b	15	phi[2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b] + plo[15]
		2c	14	phi[2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c] + plo[14]
		2d	13	phi[2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d] + plo[13]
		2e	12	phi[2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e] + plo[12]
		2f	11	phi[2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f] + plo[11]
		30	10	phi[2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30] + plo[10]
		31	0f	phi[30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31] + plo[0f]
		32	0e	phi[31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32] + plo[0e]
		33	0d	phi[32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33] + plo[0d]
		34	0c	phi[33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34] + plo[0c]
		35	0b	phi[34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35] + plo[0b]
		36	0a	phi[35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36] + plo[0a]
		37	09	phi[36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37] + plo[09]
		38	08	phi[37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38] + plo[08]
		39	07	phi[38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39] + plo[07]
		3a	06	phi[39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a] + plo[06]
		3b	05	phi[3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b] + plo[05]
		3c	04	phi[3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c] + plo[04]
		3d	03	phi[3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d] + plo[03]
		3e	02	phi[3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e] + plo[02]
		3f	01	phi[3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00] + plo[01]

	In order to make the result amenable to loop-based execution, we need to encode the
	plo-and-phi-indices in easily-computable-index fashion. This scheme has 2 ingredients:
	[Use i as loop idx in comments, although use 'k' in actual code]

	[0] plo index is simply 0,3f,3e,..,1, i.e. plo[(64-i) & (-(i > 0))]

	[1] Note that the (/= p40, where 40 = hex) row data are simply circular-shift perms of the basic (row 0) set;
	For row i (counting from 0 to 63), leftward-shift count needing to be applied = 0,0,3e,3d,..,1, same as plo idx in [0],
	except that we also need to zero the i=1 datum.
*/
//printf("Required output index permutation = [");
		tptr = t;
		for(k = 0; k < 64; ++k)
		{
			jt = (64-k) & (-(k > 0));
			iptr = dit_p40_cperms + (jt & (-(k != 1)));
			for(l = 0; l < ODD_RADIX; l++)
			{
			/*
				int i = k*ODD_RADIX + l;	// Initial run to obtain required operm uses simple in-order p0,1,..., indexing
			//if(l==0)printf("[%2x]",operm[i]&63);
			//printf("%2x,",operm[i]>>6);
				io_offsets[l] = phi[i>>6] + plo[i&63];
			*/
				io_offsets[l] = phi[iptr[l]];
			}
			RADIX_63_DIT(
				(double *)tptr, toff, 1,
			//	a+j1, io_offsets	// Initial run to obtain required operm uses simple in-order p0,1,..., indexing
				a+j1+plo[jt], io_offsets, RE_IM_STRIDE
			);
			tptr += ODD_RADIX;
		}
//printf("]\n");exit(0);
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy4032_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr,*addi;
		struct complex *tptr;
		const int pfetch_dist = PFETCH_DIST;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p1,p2,p3,p4;
		int poff[RADIX>>2];
	// DFT stuff:
	// Shared by DIF and DIT:
		int kk, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf;
		const uint8 *iptr;
		int plo[64], phi[ODD_RADIX], toff[ODD_RADIX], io_offsets[64], t_offsets[64];
		const uint8 dft_p40_cperms[128] = {	// Pad byte-array to next-higher 8-byte multiple
			0,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01,
			0,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01
		},	dft_phi[64] = {	// Only 63 entries, but again pad out to nearest 8-multiple
			0,0x02,0x01,0x08,0x07,0x06,0x05,0x04,0x03,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09
		};
	// DIF:
		/*** Now done via #include "radix4032.h" @top of this file ***/	// 4032-element dif64_oidx_lo byte-array
	// DIT:
		// Low parts [p0-3f] of output-index perms - need 64 bytes (one init row below) for each radix-64 DFT:
		const uint8 dit64_iidx_lo[64] = {	// Keep only first row...CF. comments in radix4032_dit_pass1().
			0x00,0x01,0x03,0x02,0x07,0x06,0x05,0x04,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x3f,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20
		};

		int j,j1,j2,jt,jp,k,l,ntmp;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
	#ifdef USE_AVX512
		double t0,t1,t2,t3;
		struct uint32x8 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#elif defined(USE_AVX)
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#else
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	#endif
		double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
		double rt,it;

	#ifdef USE_SSE2

		const double crnd = 3.0*0x4000000*0x2000000;
		double *add0,*add1,*add2,*add3;
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *tmp,*tm0,*tm1,*tm2;	// utility ptrs
		int *itmp,*itm2;	// Pointer into the bjmodn array
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		vec_dbl *max_err, *sse2_rnd, *half_arr,
			*r00,	// Head of RADIX*vec_cmplx-sized local store #1
			*s1p00,	// Head of RADIX*vec_cmplx-sized local store #2
			*cy_r,*cy_i;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx
	  #ifdef USE_AVX
		vec_dbl *base_negacyclic_root;
	  #endif

		/* These are used in conjunction with the langth-ODD_RADIX arrays in the USE_SCALAR_CARRY #define below;
		In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
		int idx_offset, idx_incr;
		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else

		double *base, *baseinv, *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr, bs,bsinv, wts_idx_incr;
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
		struct complex *rn0 = thread_arg->rn0;
		struct complex *rn1 = thread_arg->rn1;

		/*   constant index offsets for array load/stores are here.	*/
		// Legacy code needs 3 lowest nonzero fixed-index p[] terms:
		p1 = 1*NDIVR;	 p1 += ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = 2*NDIVR;	 p2 += ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = 3*NDIVR;	 p3 += ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = 4*NDIVR;	 p4 += ( (p4 >> DAT_BITS) << PAD_BITS );
	#ifndef USE_SSE2
		p0123[0] = 0; p0123[1] = p1; p0123[2] = p2; p0123[3] = p3;
	#endif
		for(l = 0; l < (RADIX>>2); l++) {
			poff[l] = (l<<2)*NDIVR;	// Corr. to every 4th plo[] term
			poff[l] += ( (poff[l] >> DAT_BITS) << PAD_BITS );
		}

		/* Constant padded-array index offsets for array load/stores are here. */
		for(l = 0; l < 64; l++)
		{
			plo[l] = l*NDIVR;
			plo[l] += ( (plo[l] >> DAT_BITS) << PAD_BITS );
			t_offsets[l] = ((l<<6)-l)<<1;	// 2*(l*63), extra 2x is due to cast-to-double of t-array in call to RADIX_64_DIF|T
		}
		for(l = 0; l < ODD_RADIX; l++)
		{
			phi[l] = (l<<6)*NDIVR;
			phi[l] += ( (phi[l] >> DAT_BITS) << PAD_BITS );
			toff[l] = l+l;
		}

	#ifdef USE_SSE2
		uint32 nwt16 = nwt << L2_SZ_VD;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code

		r00 = tmp = thread_arg->r00;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x1f80;	s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
		tmp += 0x1f80;
		// DFT-63 and -64 roots all def'd locally in the DFT-63,64 functions
	  #ifdef USE_AVX512
		cy_r = tmp;	cy_i = tmp+0x1f8;	tmp += 2*0x1f8;	// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod in AVX mode
	  #elif defined(USE_AVX)
		cy_r = tmp;	cy_i = tmp+0x3f0;	tmp += 2*0x3f0;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod in AVX mode
	  #else
		cy_r = tmp;	cy_i = tmp+0x7e0;	tmp += 2*0x7e0;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	  #endif

		ASSERT((r00 == thread_arg->r00), "thread-local memcheck failed!");
		ASSERT((half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
		ASSERT((sse2_rnd->d0 == crnd && sse2_rnd->d1 == crnd), "thread-local memcheck failed!");
		tmp = half_arr;
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
	  #ifdef USE_AVX
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif
	} else {
		dtmp = (tmp)->d0 * (tmp+ODD_RADIX)->d0;	ASSERT(fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp)->d1 * (tmp+ODD_RADIX)->d1;	ASSERT(fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
	}

		VEC_DBL_INIT(max_err, 0.0);

		sign_mask = (uint64*)(r00 + radix4032_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (#doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		sse_sw  = sse_bw    + RE_IM_STRIDE;
		sse_n   = sse_sw    + RE_IM_STRIDE;
	  #ifdef USE_AVX512
		n_minus_sil   = (struct uint32x8 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x8 *)sse_n + 2;
		sinwt         = (struct uint32x8 *)sse_n + 3;
		sinwtm1       = (struct uint32x8 *)sse_n + 4;
		bjmodn = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #elif defined(USE_AVX)
		n_minus_sil   = (struct uint32x4 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_n + 2;
		sinwt         = (struct uint32x4 *)sse_n + 3;
		sinwtm1       = (struct uint32x4 *)sse_n + 4;
		bjmodn = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn = (int*)(sse_n + RE_IM_STRIDE);
	  #endif

	#else

		// In scalar mode use these 2 ptrs to pass wts_idx_incr and the base/baseinv/etc array-ptrs:
		wts_idx_incr = *(int *)thread_arg->half_arr;
		base      = (double *)thread_arg->r00;
		baseinv   = base + 2;
		wt_arr    = base + 4;
		wtinv_arr = wt_arr    + ODD_RADIX;
		bs_arr    = wtinv_arr + ODD_RADIX;
		bsinv_arr = bs_arr    + ODD_RADIX;

	#endif	// USE_SSE2 ?

		// Can't simply use thread-associated values of these *cycle index arrays here, since
		// thread values must be ***read-only*** so as to retain the proper first-init values
		// on each entry to this thread-task. Instead use the bjmodn data storage block - which
		// is otherwise unused in Fermat-Mod mode - for local storage of these cycle tables.
		/*** Pointer-inits: ***/
		int *icycle = bjmodn,ic_idx;
	#ifdef USE_SSE2
		int wts_idx_inc2 = thread_arg->wts_idx_inc2;
		int *jcycle = icycle + ODD_RADIX,jc_idx;
	  #ifdef USE_AVX
		int *kcycle = jcycle + ODD_RADIX,kc_idx;
		int *lcycle = kcycle + ODD_RADIX,lc_idx;
	  #endif
	  #ifdef USE_AVX512
		int *mcycle = lcycle + ODD_RADIX,mc_idx;
		int *ncycle = mcycle + ODD_RADIX,nc_idx;
		int *ocycle = ncycle + ODD_RADIX,oc_idx;
		int *pcycle = ocycle + ODD_RADIX,pc_idx;
	  #endif
	#endif
		/*** Value-inits: ***/
		for(j = 0; j < ODD_RADIX; j++) {
			icycle[j] = thread_arg->icycle[j];
		#ifdef USE_SSE2
			jcycle[j] = thread_arg->jcycle[j];
		  #ifdef USE_AVX
			kcycle[j] = thread_arg->kcycle[j];
			lcycle[j] = thread_arg->lcycle[j];
		  #endif
		  #ifdef USE_AVX512
			mcycle[j] = thread_arg->mcycle[j];
			ncycle[j] = thread_arg->ncycle[j];
			ocycle[j] = thread_arg->ocycle[j];
			pcycle[j] = thread_arg->pcycle[j];
		  #endif
		#endif
		}

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
		#elif 0	// No_op in scalar case, since carry pattern matches that of thread data
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
		#elif 0	// No_op in scalar case, since carry pattern matches that of thread data
			for(l = 0; l < RADIX; l++) {
				cy_r[l] = *(addr+l);		cy_i[l] = *(addi+l);
			}
		#endif
		}

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix4032_main_carry_loop.h"

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
#undef ODD_RADIX
#undef PFETCH_DIST
