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
#include "radix16.h"

#define RADIX 208	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

#define EPS 1e-10

// Using a flag-value (rather than 'defined or not?' scheme allows us to override the enclosed default at compile time:
#ifdef USE_AVX2
  #ifndef DFT_13_FMA	// To force non-default value of toggle, use -DDFT_13_FMA=1 at compile time
	#define DFT_13_FMA	0	// Toggle between 'naive' but lower-ROE DFT [36 ADD, 168 FMA (incl 12 MUL), 268 memref] and lower-opcount
							// but higher-ROE van Buskirk-style 'tangent' DFT [120 ADD, 78 FMA (incl 34 MUL), 168 memref].
							// Toggle only respected for build modes where FMA instructions exist. With FMA, 'naive' DFT has
							// a comparable opcount to VB one, since it yields a target-rich environment for FMA.
  #endif
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

#ifdef USE_SSE2

  // For Mersenne-mod need (16 [SSE2] or 64 [AVX]) + (4 [HIACC] or 40 [LOACC]) added slots for half_arr lookup tables.
  // Max = (40 [SSE2]; 132 [AVX]), add to (half_arr_offset208 + RADIX) to get required value of radix208_creals_in_local_store:
  #ifdef USE_AVX512	// RADIX/8 = 0x1a fewer carry slots than AVX:
	const int half_arr_offset208 = 0x376;	// + RADIX = += 0xd0 = 0x446; Used for thread local-storage-integrity checking
	const int radix208_creals_in_local_store = 0x4cc;	// (half_arr_offset + RADIX) + 132 (=0x84) and round up to nearest multiple of 4
  #elif defined(USE_AVX)
	const int half_arr_offset208 = 0x390;	// + RADIX = += 0xd0 = 0x460; Used for thread local-storage-integrity checking
	const int radix208_creals_in_local_store = 0x4e4;	// (half_arr_offset + RADIX) + 132 (=0x84) and round up to nearest multiple of 4
  #else
	const int half_arr_offset208 = 0x3c4;	// + RADIX = 0x494; Used for thread local-storage-integrity checking
	const int radix208_creals_in_local_store = 0x4bc;	// (half_arr_offset + RADIX) = 40 and round up to nearest multiple of 4
  #endif

	#include "sse2_macro_gcc64.h"
	#include "radix13_sse_macro.h"

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

int radix208_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-208 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-208 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix208_ditN_cy_dif1";
#if !defined(MULTITHREAD) && defined(USE_SSE2)
	const int pfetch_dist = PFETCH_DIST;
#endif
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
#if DFT_13_FMA || defined(USE_ARM_V8_SIMD)
  #warning Using FMA-heavy / lo-add 13-DFT
const double cc1=  0.88545602565320989590,	/* Real part of exp(i*2*pi/13), the radix-13 fundamental sincos datum	*/
			ss1 =  0.46472317204376854565,	/* Imag part of exp(i*2*pi/13).	*/
			cc2 =  0.56806474673115580252,	/* cos(2u)	*/
			ss2 =  0.82298386589365639457,	/* sin(2u)	*/
			cc3 =  0.12053668025532305336,	/* cos(3u)	*/
			ss3 =  0.99270887409805399279,	/* sin(3u)	*/
			cc4 = -0.35460488704253562594,	/* cos(4u)	*/
			ss4 =  0.93501624268541482344,	/* sin(4u)	*/
			cc5 = -0.74851074817110109861,	/* cos(5u)	*/
			ss5 =  0.66312265824079520240,	/* sin(5u)	*/
			cc6 = -0.97094181742605202714,	/* cos(6u)	*/
			ss6 =  0.23931566428755776718;	/* sin(6u)	*/
#else	// Consts for van Buskirk-style tangent DFT in this header:
  #warning Using FMA-lite / hi-add 13-DFT
	#include "radix13.h"
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
	int j2,jp,ntmp;
  #endif
  #ifndef MULTITHREAD
	int jstart;
  #endif
  #ifdef USE_SSE2
	int nbytes;
  #endif
  #if !defined(MULTITHREAD) && defined(USE_SSE2)
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 13|26|52 for avx512,avx,sse, respectively.
	// But fixed-incr too restrictive here, so 'divide 13|26|52 into pieces' via increment-array whose elts sum to 13|26|52:
	const int *incr,*inc_arr;
  #ifdef USE_AVX512	// Have no specialized HIACC carry macro in AVX-512 and ARMv8 SIMD, so these get an "goes to 11" in LOACC mode via an incr_hiacc[] array:
	const int incr_long[] = {13}, incr_med[] = {6,7}, incr_short[] = {4,5,4}, incr_hiacc[] = {3,2,3,2,3};
  #elif defined(USE_AVX)
	const int incr_long[] = {13,13}, incr_med[] = {9,8,9}, incr_short[] = {4,5,4,4,5,4}, incr_hiacc[] = {0};
  #elif defined(USE_ARM_V8_SIMD)
	const int incr_long[] = {13,13,13,13}, incr_med[] = {9,8,9,9,8,9}, incr_short[] = {4,5,4,4,5,4,4,5,4,4,5,4}, incr_hiacc[] = {3,2,3,2,3,3,2,3,2,3,3,2,3,2,3,3,2,3,2,3};
  #else
	const int incr_long[] = {13,13,13,13}, incr_med[] = {9,8,9,9,8,9}, incr_short[] = {4,5,4,4,5,4,4,5,4,4,5,4}, incr_hiacc[] = {0};
  #endif
	// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
	if(USE_SHORT_CY_CHAIN == 0)
		inc_arr = incr_long;
	else if(USE_SHORT_CY_CHAIN == 1)
		inc_arr = incr_med;
	else if(USE_SHORT_CY_CHAIN == 2)
		inc_arr = incr_short;
	else
		inc_arr = incr_hiacc;
  #endif // !MULTITHREAD && USE_SSE2

	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
								,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0, nsave = 0;
	static int poff[RADIX>>2];	// Store [RADIX/4] mults of p04 offset for loop control
#ifndef MULTITHREAD
// DIF:
	const uint64 dif16_oidx_lo[13] = {
		0x01326745dcef89baull, 0x23015467fedcab89ull, 0x98abfedc23015467ull, 0x7654230198abfedcull, 0xcdfe98ab76542301ull,
		0x10237654cdfe98abull, 0xba98cdfe10237654ull, 0x45761023ba98cdfeull, 0xefcdba9845761023ull, 0x32104576efcdba98ull,
		0x89baefcd32104576ull, 0x6745321089baefcdull, 0xdcef89ba67453210ull
	};
	// circ-shift count of basic array needed for stage 1:
	const int dif_ncshft[16] = {0x0,0x1,0x2,0x3,0x4,0x5,0x5,0x6,0x7,0x8,0x9,0x9,0xa,0xb,0xc,0x0},
		// dif_pcshft has 13 base elts followed by repeat of first 12 of those to support circ-shift perms
		dif_pcshft[25] = {0x0,0xc,0xb,0xa,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2,0x1,0x0,0xc,0xb,0xa,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2}
	  #ifdef USE_SSE2
		// Jan 2021: need 2nd version of dif_pcshft[] with elts pre-shifted by <<(L2_SZ_VD+5) to support simple additive-address-offset computation in revised 11-DFT macro:
		,dif_pcshft2[25] = {0x0,0xc<<(L2_SZ_VD+5),0xb<<(L2_SZ_VD+5),0xa<<(L2_SZ_VD+5),0x9<<(L2_SZ_VD+5),0x8<<(L2_SZ_VD+5),0x7<<(L2_SZ_VD+5),0x6<<(L2_SZ_VD+5),0x5<<(L2_SZ_VD+5),0x4<<(L2_SZ_VD+5),0x3<<(L2_SZ_VD+5),0x2<<(L2_SZ_VD+5),0x1<<(L2_SZ_VD+5),0x0<<(L2_SZ_VD+5),0xc<<(L2_SZ_VD+5),0xb<<(L2_SZ_VD+5),0xa<<(L2_SZ_VD+5),0x9<<(L2_SZ_VD+5),0x8<<(L2_SZ_VD+5),0x7<<(L2_SZ_VD+5),0x6<<(L2_SZ_VD+5),0x5<<(L2_SZ_VD+5),0x4<<(L2_SZ_VD+5),0x3<<(L2_SZ_VD+5),0x2<<(L2_SZ_VD+5)}
	  #endif
		;
// DIT:
	const uint64 dit16_iidx_lo[13] = {
		0x01327654fedcba98ull, 0x76543210ba98dcefull, 0xba98dcef32105467ull, 0x32105467dcef98abull, 0xdcef98ab54671023ull,
		0x98abefcd10236745ull, 0x10236745efcdab89ull, 0xefcdab8967452301ull, 0x67452301ab89cdfeull, 0x23014576cdfe89baull,
		0xcdfe89ba45760132ull, 0x4576013289bafedcull, 0x89bafedc01327654ull
	};
	// circ-shift count of basic array needed for stage 2:
	const int dit_ncshft[16] = {0x0,0xb,0xc,0x0,0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xa,0xb,0xc}
	  #ifndef USE_SSE2
		// dit_pcshft has 13 base elts followed by repeat of first 12 of those to support circ-shift perms
		,dit_pcshft[25] = {0x0,0x9,0x5,0x1,0xa,0x6,0x2,0xb,0x7,0x3,0xc,0x8,0x4,0x0,0x9,0x5,0x1,0xa,0x6,0x2,0xb,0x7,0x3,0xc,0x8}
	  #else
		// Jan 2021: need 2nd version of dit_pcshft[] with elts pre-shifted by <<(L2_SZ_VD+5) to support simple additive-address-offset computation in revised SIMD 11-DFT macro:
		,dit_pcshft2[25] = {0x0,0x9<<(L2_SZ_VD+5),0x5<<(L2_SZ_VD+5),0x1<<(L2_SZ_VD+5),0xa<<(L2_SZ_VD+5),0x6<<(L2_SZ_VD+5),0x2<<(L2_SZ_VD+5),0xb<<(L2_SZ_VD+5),0x7<<(L2_SZ_VD+5),0x3<<(L2_SZ_VD+5),0xc<<(L2_SZ_VD+5),0x8<<(L2_SZ_VD+5),0x4<<(L2_SZ_VD+5),0x0<<(L2_SZ_VD+5),0x9<<(L2_SZ_VD+5),0x5<<(L2_SZ_VD+5),0x1<<(L2_SZ_VD+5),0xa<<(L2_SZ_VD+5),0x6<<(L2_SZ_VD+5),0x2<<(L2_SZ_VD+5),0xb<<(L2_SZ_VD+5),0x7<<(L2_SZ_VD+5),0x3<<(L2_SZ_VD+5),0xc<<(L2_SZ_VD+5),0x8<<(L2_SZ_VD+5),}
	  #endif
		;
// Shared by both DIF+DIT:
	static int plo[16],phi[13];
	uint64 i64;
	const int *iptr;
	int kk;
  #ifndef USE_SSE2
	int k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf;
  #else
	int po_kperm[16],*po_ptr = &(po_kperm[0]);
  #endif
  #ifdef USE_SSE2
	const uint32 dft13_offs[13] = {0,0x20<<L2_SZ_VD,0x40<<L2_SZ_VD,0x60<<L2_SZ_VD,0x80<<L2_SZ_VD,0xa0<<L2_SZ_VD,0xc0<<L2_SZ_VD,0xe0<<L2_SZ_VD,0x100<<L2_SZ_VD,0x120<<L2_SZ_VD,0x140<<L2_SZ_VD,0x160<<L2_SZ_VD,0x180<<L2_SZ_VD}, *dft13_offptr = &(dft13_offs[0]);
  #endif
#endif
	static double radix_inv, n2inv;
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
		*tmp,*tm2	// Non-static utility ptrs
	  #ifndef MULTITHREAD
		,*tm1	// Non-static utility ptrs
	  #endif
		;
	static vec_dbl *two,*one,*sqrt2,*isrt2, *rad13_const, *max_err, *sse2_rnd, *half_arr, *cc0, *ss0,	// rad13_const needs 18*sizeof(vec_dbl) bytes
		*r00	// Head of RADIX*vec_cmplx-sized local store #1
	  #ifndef MULTITHREAD
		,*s1p00	// Head of RADIX*vec_cmplx-sized local store #2
		,*cy	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx
	  #endif
		;
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
	static task_control_t   task_control = {NULL, (void*)cy208_process_chunk, NULL, 0x0};

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
		// consisting of 88 dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix208_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix208_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = sc_ptr;	r00   = tmp;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x1a0;					// Head of RADIX*vec_cmplx-sized local store #2
	  #ifndef MULTITHREAD
						s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
	  #endif
		tmp += 0x1a0;
		two   = tmp + 0;	// AVX+ versions of radix-8,16,32 twiddleless-DFT macros need consts [2,1,sqrt2,isrt2] quartet laid out thusly
		one   = tmp + 1;
		sqrt2 = tmp + 2;
		isrt2 = tmp + 3;
		cc0   = tmp + 4;
		ss0   = tmp + 5;
	//	one   = tmp + 6;	Unnamed 1.0 slot to be used in radix-13
	//	two   = tmp + 7;	Unnamed 2.0 slot to be used in radix-13
		rad13_const = tmp + 0x08;	// Needs 17 vec_dbl slots
		tmp += 0x1a;	// Need 8 + 17 = 25 vec_dbl slots for DFT sincos; round up nearest even
	// sc_ptr += 0x36a
	  #ifdef USE_AVX512
	   #ifndef MULTITHREAD
		cy = tmp;						// RADIX/8 vec_dbl slots for carry sub-array
	   #endif
						tmp += 0x1a;	// RADIX/8 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
	   #ifndef MULTITHREAD
		cy = tmp;						// RADIX/4 vec_dbl slots for carry sub-array
	   #endif
						tmp += 0x34;	// RADIX/4 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x390; This is where the value of half_arr_offset208 comes from
		half_arr= tmp + 0x02;	// This table needs 96*SZ_VD bytes in avx mode
	  #else
	   #ifndef MULTITHREAD
		cy = tmp;						// RADIX/2 vec_dbl slots for carry sub-array
	   #endif
						tmp += 0x68;	// RADIX/2 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x3c4; This is where the value of half_arr_offset208 comes from
		half_arr= tmp + 0x02;	// This table needs 32*SZ_VD bytes in sse2 mode
	  #endif
		ASSERT((radix208_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r00) + (20 << L2_SZ_VD), "radix208_creals_in_local_store checksum failed!");
		/* These remain fixed: */
		VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );
	  #if 1
		// 2 unnamed slots for alternate "rounded the other way" copies of sqrt2,isrt2:
		dtmp = u64_to_f64(sqrt2_dn);	VEC_DBL_INIT(sqrt2, dtmp);
		dtmp = u64_to_f64(isrt2_dn);	VEC_DBL_INIT(isrt2, dtmp);
	  #else
		VEC_DBL_INIT(sqrt2, SQRT2);		VEC_DBL_INIT(isrt2, ISRT2);
	  #endif
		VEC_DBL_INIT(cc0  ,  c16);
		VEC_DBL_INIT(ss0  ,  s16);
		tmp = rad13_const-2;		/* __cc pointer offsets: */
		VEC_DBL_INIT(tmp,  1.0);	++tmp;	/*	-0x020 = 1.0 */
		VEC_DBL_INIT(tmp,  2.0);	++tmp;	/*	-0x010 = 2.0 */
	  #if DFT_13_FMA || defined(USE_ARM_V8_SIMD)	// FMA version (like radix-11/FMA) based on simple radix-13 DFT implementation with lots of opportunities for FMAs
		VEC_DBL_INIT(tmp, cc1 );	++tmp;	/*	0x000 = cc1 */
		VEC_DBL_INIT(tmp, cc2 );	++tmp;	/*	0x010 = cc2 */
		VEC_DBL_INIT(tmp, cc3 );	++tmp;	/*	0x020 = cc3 */
		VEC_DBL_INIT(tmp, cc4 );	++tmp;	/*	0x030 = cc4 */
		VEC_DBL_INIT(tmp, cc5 );	++tmp;	/*	0x040 = cc5 */
		VEC_DBL_INIT(tmp, cc6 );	++tmp;	/*	0x050 = cc6 */
		VEC_DBL_INIT(tmp, ss1 );	++tmp;	/*	0x060 = ss1 */
		VEC_DBL_INIT(tmp, ss2 );	++tmp;	/*	0x070 = ss2 */
		VEC_DBL_INIT(tmp, ss3 );	++tmp;	/*	0x080 = ss3 */
		VEC_DBL_INIT(tmp, ss4 );	++tmp;	/*	0x090 = ss4 */
		VEC_DBL_INIT(tmp, ss5 );	++tmp;	/*	0x0a0 = ss5 */
		VEC_DBL_INIT(tmp, ss6 );	++tmp;	/*	0x0b0 = ss6 */
		VEC_DBL_INIT(tmp, 0.0 );	++tmp;	/*	0x0c0 = 0.0 */// Upper 5 slots unused here; init = 0
		VEC_DBL_INIT(tmp, 0.0 );	++tmp;	/*	0x0d0 = 0.0 */
		VEC_DBL_INIT(tmp, 0.0 );	++tmp;	/*	0x0e0 = 0.0 */
		VEC_DBL_INIT(tmp, 0.0 );	++tmp;	/*	0x0f0 = 0.0 */
		VEC_DBL_INIT(tmp, 0.0 );	++tmp;	/*	0x100 = 0.0 */
	  #else
		VEC_DBL_INIT(tmp,  DC1);	++tmp;	/*	0x000 =  DC1 */
		VEC_DBL_INIT(tmp,  DC3);	++tmp;	/*	0x010 =  DC3 */
		VEC_DBL_INIT(tmp,  DC4);	++tmp;	/*	0x020 =  DC4 */
		VEC_DBL_INIT(tmp,  DS1);	++tmp;	/*	0x030 =  DS1 */
		VEC_DBL_INIT(tmp,  DS2);	++tmp;	/*	0x040 =  DS2 */
		VEC_DBL_INIT(tmp,  DS3);	++tmp;	/*	0x050 =  DS3 */
		VEC_DBL_INIT(tmp,  DS4);	++tmp;	/*	0x060 =  DS4 */
		VEC_DBL_INIT(tmp,  DS5);	++tmp;	/*	0x070 =  DS5 */
		VEC_DBL_INIT(tmp, DC23);	++tmp;	/*	0x080 = DC23 */
		VEC_DBL_INIT(tmp, DC54);	++tmp;	/*	0x090 = DC54 */
		VEC_DBL_INIT(tmp, DC65);	++tmp;	/*	0x0a0 = DC65 */
		VEC_DBL_INIT(tmp, DS63);	++tmp;	/*	0x0b0 = DS63 */
		VEC_DBL_INIT(tmp, DS74);	++tmp;	/*	0x0c0 = DS74 */
		VEC_DBL_INIT(tmp, DS85);	++tmp;	/*	0x0d0 = DS85 */
		VEC_DBL_INIT(tmp, DS93);	++tmp;	/*	0x0e0 = DS93 */
		VEC_DBL_INIT(tmp, DSa4);	++tmp;	/*	0x0f0 = DSa4 */
		VEC_DBL_INIT(tmp, DSb5);	++tmp;	/*	0x100 = DSb5 */
	  #endif
		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		sse2_rnd->d0 = base[0]; sse2_rnd->d1 = baseinv[1]; sse2_rnd->d2 = wts_mult[1]; sse2_rnd->d3 = inv_mult[0];
	  #else
		VEC_DBL_INIT(sse2_rnd, crnd);
	  #endif

		// Propagate the above consts to the remaining threads:
		nbytes = (intptr_t)tmp - (intptr_t)two;	// #bytes in above sincos block of data
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
		NDIVR >>= 4;			pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
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

	#ifndef MULTITHREAD
	// Shared by both DIF+DIT:
		plo[0x0] =  0; plo[0x1] = p1; plo[0x2] = p2; plo[0x3] = p3;
		plo[0x4] = p4; plo[0x5] = p5; plo[0x6] = p6; plo[0x7] = p7;
		plo[0x8] = p8; plo[0x9] = p9; plo[0xa] = pa; plo[0xb] = pb;
		plo[0xc] = pc; plo[0xd] = pd; plo[0xe] = pe; plo[0xf] = pf;
		l = 0;
		phi[l++] =   0; phi[l++] = p10; phi[l++] = p20; phi[l++] = p30; phi[l++] = p40;
		phi[l++] = p50; phi[l++] = p60; phi[l++] = p70; phi[l++] = p80; phi[l++] = p90;
		phi[l++] = pa0; phi[l++] = pb0; phi[l++] = pc0;
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

/*...The radix-208 final DIT pass is here.	*/

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
			_jhi[ithread] = _jstart[ithread] + jhi_wrap;	/* Cleanup loop assumes carryins propagate at most 4 words up. */
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
		#include "radix208_main_carry_loop.h"

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
		ASSERT(0x0 == cy208_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
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

void radix208_dif_pass1(double a[], int n)
{
#include "radix13.h"	// In these wrappers we always use the tangent-DFT
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-208 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix13_dif_pass for details on the radix-13 subtransform.
*/
	int l,j,j1,j2,jt,jp;
	const int *iptr;
	// p-indexing is hexadecimal here:
	static int NDIVR,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0, first_entry=TRUE;
	const uint64 dif16_oidx_lo[13] = {
		0x01326745dcef89baull, 0x23015467fedcab89ull, 0x98abfedc23015467ull, 0x7654230198abfedcull, 0xcdfe98ab76542301ull,
		0x10237654cdfe98abull, 0xba98cdfe10237654ull, 0x45761023ba98cdfeull, 0xefcdba9845761023ull, 0x32104576efcdba98ull,
		0x89baefcd32104576ull, 0x6745321089baefcdull, 0xdcef89ba67453210ull
	};
	static int plo[16],phi[13];
	// circ-shift count of basic array needed for stage 1:
	const int dif_ncshft[16] = {0x0,0x1,0x2,0x3,0x4,0x5,0x5,0x6,0x7,0x8,0x9,0x9,0xa,0xb,0xc,0x0},
		// dif_pcshft has 13 base elts followed by repeat of first 12 of those to support circ-shift perms
		dif_pcshft[25] = {0x0,0xc,0xb,0xa,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2,0x1,0x0,0xc,0xb,0xa,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2};
	uint64 i64;
	int k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf;
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
		pa0 = p90 + NDIVR;		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pb0 = pa0 + NDIVR;		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pc0 = pb0 + NDIVR;		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );

		plo[0x0] =  0; plo[0x1] = p1; plo[0x2] = p2; plo[0x3] = p3;
		plo[0x4] = p4; plo[0x5] = p5; plo[0x6] = p6; plo[0x7] = p7;
		plo[0x8] = p8; plo[0x9] = p9; plo[0xa] = pa; plo[0xb] = pb;
		plo[0xc] = pc; plo[0xd] = pd; plo[0xe] = pe; plo[0xf] = pf;
		l = 0;
		phi[l++] =   0; phi[l++] = p10; phi[l++] = p20; phi[l++] = p30; phi[l++] = p40;
		phi[l++] = p50; phi[l++] = p60; phi[l++] = p70; phi[l++] = p80; phi[l++] = p90;
		phi[l++] = pa0; phi[l++] = pb0; phi[l++] = pc0;
	}

/*...The radix-208 pass is here.	*/

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
	Twiddleless version arranges 16 sets of radix-13 DFT inputs as follows:
	0 in upper left corner, decrement 16 horizontally and 13 vertically (mod 208).
	Use hex here to match p-indexing ...Can auto-generate these by running test_fft_radix with TTYPE = 0;
	note this input-offset pattern is shared by DIF and DIT, but DIT layers a generalized bit-reversal atop it:

	DIF/DIT input-scramble array =              [vvv 00,...,10 vvv = basic offset array:]		leftward-circ-shift count of basic array
		00,c0,b0,a0,90,80,70,60,50,40,30,20,10   00,c0,b0,a0,90,80,70,60,50,40,30,20,10 + p0	0
		c3,b3,a3,93,83,73,63,53,43,33,23,13,03   c0,b0,a0,90,80,70,60,50,40,30,20,10,00 + p3	1
		b6,a6,96,86,76,66,56,46,36,26,16,06,c6   b0,a0,90,80,70,60,50,40,30,20,10,00,c0 + p6	2
		a9,99,89,79,69,59,49,39,29,19,09,c9,b9   a0,90,80,70,60,50,40,30,20,10,00,c0,b0 + p9	3
		9c,8c,7c,6c,5c,4c,3c,2c,1c,0c,cc,bc,ac   90,80,70,60,50,40,30,20,10,00,c0,b0,a0 + pc	4
		8f,7f,6f,5f,4f,3f,2f,1f,0f,cf,bf,af,9f   80,70,60,50,40,30,20,10,00,c0,b0,a0,90 + pf	5
		82,72,62,52,42,32,22,12,02,c2,b2,a2,92   80,70,60,50,40,30,20,10,00,c0,b0,a0,90 + p2	5
		75,65,55,45,35,25,15,05,c5,b5,a5,95,85   70,60,50,40,30,20,10,00,c0,b0,a0,90,80 + p5	6
		68,58,48,38,28,18,08,c8,b8,a8,98,88,78 = 60,50,40,30,20,10,00,c0,b0,a0,90,80,70 + p8	7
		5b,4b,3b,2b,1b,0b,cb,bb,ab,9b,8b,7b,6b   50,40,30,20,10,00,c0,b0,a0,90,80,70,60 + pb	8
		4e,3e,2e,1e,0e,ce,be,ae,9e,8e,7e,6e,5e   40,30,20,10,00,c0,b0,a0,90,80,70,60,50 + pe	9
		41,31,21,11,01,c1,b1,a1,91,81,71,61,51   40,30,20,10,00,c0,b0,a0,90,80,70,60,50 + p1	9
		34,24,14,04,c4,b4,a4,94,84,74,64,54,44   30,20,10,00,c0,b0,a0,90,80,70,60,50,40 + p4	a
		27,17,07,c7,b7,a7,97,87,77,67,57,47,37   20,10,00,c0,b0,a0,90,80,70,60,50,40,30 + p7	b
		1a,0a,ca,ba,aa,9a,8a,7a,6a,5a,4a,3a,2a   10,00,c0,b0,a0,90,80,70,60,50,40,30,20 + pa	c
		0d,cd,bd,ad,9d,8d,7d,6d,5d,4d,3d,2d,1d   00,c0,b0,a0,90,80,70,60,50,40,30,20,10 + pd	0
	*/
	/*...gather the needed data (208 64-bit complex) and do 16 radix-13 transforms...*/
		tptr = t;
		for(l = 0; l < 16; l++) {
			iptr = dif_pcshft + dif_ncshft[l];
			// Hi-part of p-offset indices:
			k0 = phi[*iptr];
			k1 = phi[*(iptr+0x1)];
			k2 = phi[*(iptr+0x2)];
			k3 = phi[*(iptr+0x3)];
			k4 = phi[*(iptr+0x4)];
			k5 = phi[*(iptr+0x5)];
			k6 = phi[*(iptr+0x6)];
			k7 = phi[*(iptr+0x7)];
			k8 = phi[*(iptr+0x8)];
			k9 = phi[*(iptr+0x9)];
			ka = phi[*(iptr+0xa)];
			kb = phi[*(iptr+0xb)];
			kc = phi[*(iptr+0xc)];
			jp = plo[((l<<1)+l) & 0xf];	// Low-part offset = p[3*l (mod 16)] ...
			jt = j1 + jp; jp += j2;		// ... = p0,3,6,9,c,f,2,5,...
			RADIX_13_DFT(
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],
				tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im
			);	tptr++;
		}
	/*...and now do 13 radix-16 transforms.
	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF]
	to properly permute the outputs of the radix-16 DFTs to the required ordering, which in terms of our p-offsets is

		00,01,03,02,06,07,04,05,0d,0c,0e,0f,08,09,0b,0a   0,1,3,2,6,7,4,5,d,c,e,f,8,9,b,a
		c2,c3,c0,c1,c5,c4,c6,c7,cf,ce,cd,cc,ca,cb,c8,c9   2,3,0,1,5,4,6,7,f,e,d,c,a,b,8,9 + pc0
		b9,b8,ba,bb,bf,be,bd,bc,b2,b3,b0,b1,b5,b4,b6,b7   9,8,a,b,f,e,d,c,2,3,0,1,5,4,6,7 + pb0
		a7,a6,a5,a4,a2,a3,a0,a1,a9,a8,aa,ab,af,ae,ad,ac   7,6,5,4,2,3,0,1,9,8,a,b,f,e,d,c + pa0
		9c,9d,9f,9e,99,98,9a,9b,97,96,95,94,92,93,90,91   c,d,f,e,9,8,a,b,7,6,5,4,2,3,0,1 + p90
		81,80,82,83,87,86,85,84,8c,8d,8f,8e,89,88,8a,8b   1,0,2,3,7,6,5,4,c,d,f,e,9,8,a,b + p80
		7b,7a,79,78,7c,7d,7f,7e,71,70,72,73,77,76,75,74 = b,a,9,8,c,d,f,e,1,0,2,3,7,6,5,4 + p70
		64,65,67,66,61,60,62,63,6b,6a,69,68,6c,6d,6f,6e   4,5,7,6,1,0,2,3,b,a,9,8,c,d,f,e + p60
		5e,5f,5c,5d,5b,5a,59,58,54,55,57,56,51,50,52,53   e,f,c,d,b,a,9,8,4,5,7,6,1,0,2,3 + p50
		43,42,41,40,44,45,47,46,4e,4f,4c,4d,4b,4a,49,48   3,2,1,0,4,5,7,6,e,f,c,d,b,a,9,8 + p40
		38,39,3b,3a,3e,3f,3c,3d,33,32,31,30,34,35,37,36   8,9,b,a,e,f,c,d,3,2,1,0,4,5,7,6 + p30
		26,27,24,25,23,22,21,20,28,29,2b,2a,2e,2f,2c,2d   6,7,4,5,3,2,1,0,8,9,b,a,e,f,c,d + p20
		1d,1c,1e,1f,18,19,1b,1a,16,17,14,15,13,12,11,10   d,c,e,f,8,9,b,a,6,7,4,5,3,2,1,0 + p10

	In compact-obj-code mode we define the following "submotifs" for the various perms of the
	index quartets 0-3,4-7,8-b,c-f appearing to left of the '+' in the above rcol expressions:

		a0: 0,1,3,2		b0: 4,5,7,6		c0: 8,9,b,a		d0: c,d,f,e
		a1: 1,0,2,3		b1: 5,4,6,7		c1: 9,8,a,b		d1: d,c,e,f
		a2: 2,3,0,1		b2: 6,7,4,5		c2: a,b,8,9		d2: e,f,c,d
		a3: 3,2,1,0		b3: 7,6,5,4		c3: b,a,9,8		d3: f,e,d,c

	In terms of which the output ordering 16-tets become

		[a0],[b2],[d1],[c0]
		[a2],[b1],[d3],[c2] + pc0
		[c1],[d3],[a2],[b1] + pb0
		[b3],[a2],[c1],[d3] + pa0
		[d0],[c1],[b3],[a2] + p90
		[a1],[b3],[d0],[c1] + p80
		[c3],[d0],[a1],[b3] + p70
		[b0],[a1],[c3],[d0] + p60
		[d2],[c3],[b0],[a1] + p50
		[a3],[b0],[d2],[c3] + p40
		[c0],[d2],[a3],[b0] + p30
		[b2],[a3],[c0],[d2] + p20
		[d1],[c0],[b2],[a3] + p10 .

	We can encode each [a-d][0-3] motif index via 4 bits, thus each output row needs 4 x 4 = 16 bits ...
	but much easier, albeit less "compression efficient" to directly convert each 16-perm above
	into the hex digits of a little-endian uint64.
	*/
	//...and now do 13 radix-16 transforms:
		tptr = t;
		for(l = 0; l < 13; l++) {
			i64 = dif16_oidx_lo[l];
			// p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];
			k1 = plo[(i64 >> 56)&0xf];
			k2 = plo[(i64 >> 52)&0xf];
			k3 = plo[(i64 >> 48)&0xf];
			k4 = plo[(i64 >> 44)&0xf];
			k5 = plo[(i64 >> 40)&0xf];
			k6 = plo[(i64 >> 36)&0xf];
			k7 = plo[(i64 >> 32)&0xf];
			k8 = plo[(i64 >> 28)&0xf];
			k9 = plo[(i64 >> 24)&0xf];
			ka = plo[(i64 >> 20)&0xf];
			kb = plo[(i64 >> 16)&0xf];
			kc = plo[(i64 >> 12)&0xf];
			kd = plo[(i64 >>  8)&0xf];
			ke = plo[(i64 >>  4)&0xf];
			kf = plo[(i64      )&0xf];
			jp = phi[dif_pcshft[l]];	// = p0,pc0,pb0,...,p10
			jt = j1 + jp; jp += j2;
			RADIX_16_DIF(
				tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],a[jt+kf],a[jp+kf],
				c16,s16);	tptr += 0x10;
		}
	}
}

/***************/

void radix208_dit_pass1(double a[], int n)
{
#include "radix13.h"	// In these wrappers we always use the tangent-DFT
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-208 complex DIT FFT pass on the data in the length-N real vector A.
*/
	int l,j,j1,j2,jt,jp;
	const int *iptr;
	// p-indexing is hexadecimal here:
	static int NDIVR,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0, first_entry=TRUE;
	const uint64 dit16_iidx_lo[13] = {
		0x01327654fedcba98ull, 0x76543210ba98dcefull, 0xba98dcef32105467ull, 0x32105467dcef98abull, 0xdcef98ab54671023ull,
		0x98abefcd10236745ull, 0x10236745efcdab89ull, 0xefcdab8967452301ull, 0x67452301ab89cdfeull, 0x23014576cdfe89baull,
		0xcdfe89ba45760132ull, 0x4576013289bafedcull, 0x89bafedc01327654ull
	};
	static int plo[16],phi[13];
	// circ-shift count of basic array needed for stage 2:
	const int dit_ncshft[16] = {0x0,0xb,0xc,0x0,0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xa,0xb,0xc},
		// dit_pcshft has 13 base elts followed by repeat of first 12 of those to support circ-shift perms
		dit_pcshft[25] = {0x0,0x9,0x5,0x1,0xa,0x6,0x2,0xb,0x7,0x3,0xc,0x8,0x4,0x0,0x9,0x5,0x1,0xa,0x6,0x2,0xb,0x7,0x3,0xc,0x8};
	uint64 i64;
	int kk, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf;
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
		pa0 = p90 + NDIVR;		p90 += ( (p90 >> DAT_BITS) << PAD_BITS );
		pb0 = pa0 + NDIVR;		pa0 += ( (pa0 >> DAT_BITS) << PAD_BITS );
		pc0 = pb0 + NDIVR;		pb0 += ( (pb0 >> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );

		plo[0x0] =  0; plo[0x1] = p1; plo[0x2] = p2; plo[0x3] = p3;
		plo[0x4] = p4; plo[0x5] = p5; plo[0x6] = p6; plo[0x7] = p7;
		plo[0x8] = p8; plo[0x9] = p9; plo[0xa] = pa; plo[0xb] = pb;
		plo[0xc] = pc; plo[0xd] = pd; plo[0xe] = pe; plo[0xf] = pf;
		l = 0;
		phi[l++] =   0; phi[l++] = p10; phi[l++] = p20; phi[l++] = p30; phi[l++] = p40;
		phi[l++] = p50; phi[l++] = p60; phi[l++] = p70; phi[l++] = p80; phi[l++] = p90;
		phi[l++] = pa0; phi[l++] = pb0; phi[l++] = pc0;

		// Stage 1 [the 16 radix-13 DFTs] needs circ-perms of this, so pad basic 16-elt array with copy of first 15 elts:
		l = 0;
	}

/*...The radix-208 pass is here.	*/

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

	/*...gather the needed data (208 64-bit complex) and do 13 radix-16 transforms:

	Twiddleless version uses same linear-index-vector-form permutation as in DIF -
	Remember, inputs to DIT are bit-reversed, so use output of test_fft_radix() with
	TTYPE=0 to auto-generate needed input-index permutation:

	Combined DIT input-scramble array =
		00,01,03,02,07,06,05,04,0f,0e,0d,0c,0b,0a,09,08   0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00
		a7,a6,a5,a4,a3,a2,a1,a0,ab,aa,a9,a8,ad,ac,ae,af   7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + pa0
		7b,7a,79,78,7d,7c,7e,7f,73,72,71,70,75,74,76,77   b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p70
		43,42,41,40,45,44,46,47,4d,4c,4e,4f,49,48,4a,4b   3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p40
		1d,1c,1e,1f,19,18,1a,1b,15,14,16,17,11,10,12,13   d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + p10
		b9,b8,ba,bb,be,bf,bc,bd,b1,b0,b2,b3,b6,b7,b4,b5   9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + pb0
		81,80,82,83,86,87,84,85,8e,8f,8c,8d,8a,8b,88,89 = 1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + p80
		5e,5f,5c,5d,5a,5b,58,59,56,57,54,55,52,53,50,51   e,f,c,d,a,b,8,9,6,7,4,5,2,3,0,1 + p50
		26,27,24,25,22,23,20,21,2a,2b,28,29,2c,2d,2f,2e   6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + p20
		c2,c3,c0,c1,c4,c5,c7,c6,cc,cd,cf,ce,c8,c9,cb,ca   2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + pc0
		9c,9d,9f,9e,98,99,9b,9a,94,95,97,96,90,91,93,92   c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p90
		64,65,67,66,60,61,63,62,68,69,6b,6a,6f,6e,6d,6c   4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + p60
		38,39,3b,3a,3f,3e,3d,3c,30,31,33,32,37,36,35,34   8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + p30

	...which has the same pattern of p*0 offsets as the outputs of the DIF, but different row-internal permutation patterns.

	In order to support a compact-object-code version of the above simply encode each 16-perm as a hex-char string:
	*NOTE* this means we must extract each p-offset in little-endian fashion, e.g. low 4 bits have rightmost p-offset above.
	*/
		tptr = t;
		kk = 0;
		for(l = 0; l < 13; l++) {
			i64 = dit16_iidx_lo[l];
			// p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];
			k1 = plo[(i64 >> 56)&0xf];
			k2 = plo[(i64 >> 52)&0xf];
			k3 = plo[(i64 >> 48)&0xf];
			k4 = plo[(i64 >> 44)&0xf];
			k5 = plo[(i64 >> 40)&0xf];
			k6 = plo[(i64 >> 36)&0xf];
			k7 = plo[(i64 >> 32)&0xf];
			k8 = plo[(i64 >> 28)&0xf];
			k9 = plo[(i64 >> 24)&0xf];
			ka = plo[(i64 >> 20)&0xf];
			kb = plo[(i64 >> 16)&0xf];
			kc = plo[(i64 >> 12)&0xf];
			kd = plo[(i64 >>  8)&0xf];
			ke = plo[(i64 >>  4)&0xf];
			kf = plo[(i64      )&0xf];
			jp = phi[kk];	// = p10*[0,a,7,4,1,b,8,5,2,c,9,6,3], start idx = 0 and decr 3 (mod 13) each loop
			jt = j1 + jp; jp += j2;
			RADIX_16_DIT(
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],a[jt+kf],a[jp+kf],
				tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,(tptr+0xf)->re,(tptr+0xf)->im,
				c16,s16
			);	tptr += 0x10;
			kk -= 3; kk += ((-(kk < 0)) & 13);
		}
	/*...and now do 16 radix-13 transforms:
	Since our first-look oindex ordering was +p0x[0,10,20,30,40,50,60,70,80,90,a0,b0,c0] for each radix-13 and incrementing += p1 between those DFTs,
	arrange resulting mismatched-data-sorted index permutation into 13 vertical 16-entry columns to get needed oindex patterns.
	Indexing in hex for clarity and using [evn|odd]0-6 notation in the rightmost column to flag reusable 7-perms
	[in fact simple circular (0-6)-element shifts of a basic pattern] of (0-6)*p20 and p10 + (0-6)*p20:

												[vvv 00,...,40 vvv = basic offset array:]		leftward-circ-shift count of basic array
		00,90,50,10,a0,60,20,b0,70,30,c0,80,40   00,90,50,10,a0,60,20,b0,70,30,c0,80,40 + p0	0
		8f,4f,0f,9f,5f,1f,af,6f,2f,bf,7f,3f,cf   80,40,00,90,50,10,a0,60,20,b0,70,30,c0 + pf	b
		4e,0e,9e,5e,1e,ae,6e,2e,be,7e,3e,ce,8e   40,00,90,50,10,a0,60,20,b0,70,30,c0,80 + pe	c
		0d,9d,5d,1d,ad,6d,2d,bd,7d,3d,cd,8d,4d   00,90,50,10,a0,60,20,b0,70,30,c0,80,40 + pd	0
		9c,5c,1c,ac,6c,2c,bc,7c,3c,cc,8c,4c,0c   90,50,10,a0,60,20,b0,70,30,c0,80,40,00 + pc	1
		5b,1b,ab,6b,2b,bb,7b,3b,cb,8b,4b,0b,9b   50,10,a0,60,20,b0,70,30,c0,80,40,00,90 + pb	2
		1a,aa,6a,2a,ba,7a,3a,ca,8a,4a,0a,9a,5a   10,a0,60,20,b0,70,30,c0,80,40,00,90,50 + pa	3
		a9,69,29,b9,79,39,c9,89,49,09,99,59,19   a0,60,20,b0,70,30,c0,80,40,00,90,50,10 + p9	4
		68,28,b8,78,38,c8,88,48,08,98,58,18,a8 = 60,20,b0,70,30,c0,80,40,00,90,50,10,a0 + p8	5
		27,b7,77,37,c7,87,47,07,97,57,17,a7,67   20,b0,70,30,c0,80,40,00,90,50,10,a0,60 + p7	6
		b6,76,36,c6,86,46,06,96,56,16,a6,66,26   b0,70,30,c0,80,40,00,90,50,10,a0,60,20 + p6	7
		75,35,c5,85,45,05,95,55,15,a5,65,25,b5   70,30,c0,80,40,00,90,50,10,a0,60,20,b0 + p5	8
		34,c4,84,44,04,94,54,14,a4,64,24,b4,74   30,c0,80,40,00,90,50,10,a0,60,20,b0,70 + p4	9
		c3,83,43,03,93,53,13,a3,63,23,b3,73,33   c0,80,40,00,90,50,10,a0,60,20,b0,70,30 + p3	a
		82,42,02,92,52,12,a2,62,22,b2,72,32,c2   80,40,00,90,50,10,a0,60,20,b0,70,30,c0 + p2	b
		41,01,91,51,11,a1,61,21,b1,71,31,c1,81   40,00,90,50,10,a0,60,20,b0,70,30,c0,80 + p1	c
	*/
		tptr = t;
		for(l = 0; l < 16; l++) {
			iptr = dit_pcshft + dit_ncshft[l];
			// Hi-part of p-offset indices:
			k0 = phi[*iptr];
			k1 = phi[*(iptr+0x1)];
			k2 = phi[*(iptr+0x2)];
			k3 = phi[*(iptr+0x3)];
			k4 = phi[*(iptr+0x4)];
			k5 = phi[*(iptr+0x5)];
			k6 = phi[*(iptr+0x6)];
			k7 = phi[*(iptr+0x7)];
			k8 = phi[*(iptr+0x8)];
			k9 = phi[*(iptr+0x9)];
			ka = phi[*(iptr+0xa)];
			kb = phi[*(iptr+0xb)];
			kc = phi[*(iptr+0xc)];
			jp = plo[(16 - l)&0xf];	// Low-part offset = p0,f,e,...,2,1
			jt = j1 + jp; jp += j2;
			RADIX_13_DFT(
				tptr->re,tptr->im,(tptr+0x10)->re,(tptr+0x10)->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x30)->re,(tptr+0x30)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x50)->re,(tptr+0x50)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x70)->re,(tptr+0x70)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0x90)->re,(tptr+0x90)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xb0)->re,(tptr+0xb0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc]
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
	cy208_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr;
	#ifdef USE_SSE2
		const int pfetch_dist = PFETCH_DIST;
	#endif
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
					,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0;
		int poff[RADIX>>2];	// Store [RADIX/4] mults of p04 offset for loop control
	#ifndef USE_SSE2
		double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	#endif
	// DIF:
		const uint64 dif16_oidx_lo[13] = {
			0x01326745dcef89baull, 0x23015467fedcab89ull, 0x98abfedc23015467ull, 0x7654230198abfedcull, 0xcdfe98ab76542301ull,
			0x10237654cdfe98abull, 0xba98cdfe10237654ull, 0x45761023ba98cdfeull, 0xefcdba9845761023ull, 0x32104576efcdba98ull,
			0x89baefcd32104576ull, 0x6745321089baefcdull, 0xdcef89ba67453210ull
		};
		// circ-shift count of basic array needed for stage 1:
		const int dif_ncshft[16] = {0x0,0x1,0x2,0x3,0x4,0x5,0x5,0x6,0x7,0x8,0x9,0x9,0xa,0xb,0xc,0x0},
			// dif_pcshft has 13 base elts followed by repeat of first 12 of those to support circ-shift perms
			dif_pcshft[25] = {0x0,0xc,0xb,0xa,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2,0x1,0x0,0xc,0xb,0xa,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2}
		  #ifdef USE_SSE2
			// Jan 2021: need 2nd version of dif_pcshft[] with elts pre-shifted by <<(L2_SZ_VD+5) to support simple additive-address-offset computation in revised 11-DFT macro:
			,dif_pcshft2[25] = {0x0,0xc<<(L2_SZ_VD+5),0xb<<(L2_SZ_VD+5),0xa<<(L2_SZ_VD+5),0x9<<(L2_SZ_VD+5),0x8<<(L2_SZ_VD+5),0x7<<(L2_SZ_VD+5),0x6<<(L2_SZ_VD+5),0x5<<(L2_SZ_VD+5),0x4<<(L2_SZ_VD+5),0x3<<(L2_SZ_VD+5),0x2<<(L2_SZ_VD+5),0x1<<(L2_SZ_VD+5),0x0<<(L2_SZ_VD+5),0xc<<(L2_SZ_VD+5),0xb<<(L2_SZ_VD+5),0xa<<(L2_SZ_VD+5),0x9<<(L2_SZ_VD+5),0x8<<(L2_SZ_VD+5),0x7<<(L2_SZ_VD+5),0x6<<(L2_SZ_VD+5),0x5<<(L2_SZ_VD+5),0x4<<(L2_SZ_VD+5),0x3<<(L2_SZ_VD+5),0x2<<(L2_SZ_VD+5)}
		  #endif
			;
	// DIT:
		const uint64 dit16_iidx_lo[13] = {
			0x01327654fedcba98ull, 0x76543210ba98dcefull, 0xba98dcef32105467ull, 0x32105467dcef98abull, 0xdcef98ab54671023ull,
			0x98abefcd10236745ull, 0x10236745efcdab89ull, 0xefcdab8967452301ull, 0x67452301ab89cdfeull, 0x23014576cdfe89baull,
			0xcdfe89ba45760132ull, 0x4576013289bafedcull, 0x89bafedc01327654ull
		};
		// circ-shift count of basic array needed for stage 2:
		const int dit_ncshft[16] = {0x0,0xb,0xc,0x0,0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xa,0xb,0xc}
		  #ifndef USE_SSE2
			// dit_pcshft has 13 base elts followed by repeat of first 12 of those to support circ-shift perms
			,dit_pcshft[25] = {0x0,0x9,0x5,0x1,0xa,0x6,0x2,0xb,0x7,0x3,0xc,0x8,0x4,0x0,0x9,0x5,0x1,0xa,0x6,0x2,0xb,0x7,0x3,0xc,0x8}
		  #else
			// Jan 2021: need 2nd version of dit_pcshft[] with elts pre-shifted by <<(L2_SZ_VD+5) to support simple additive-address-offset computation in revised SIMD 11-DFT macro:
			,dit_pcshft2[25] = {0x0,0x9<<(L2_SZ_VD+5),0x5<<(L2_SZ_VD+5),0x1<<(L2_SZ_VD+5),0xa<<(L2_SZ_VD+5),0x6<<(L2_SZ_VD+5),0x2<<(L2_SZ_VD+5),0xb<<(L2_SZ_VD+5),0x7<<(L2_SZ_VD+5),0x3<<(L2_SZ_VD+5),0xc<<(L2_SZ_VD+5),0x8<<(L2_SZ_VD+5),0x4<<(L2_SZ_VD+5),0x0<<(L2_SZ_VD+5),0x9<<(L2_SZ_VD+5),0x5<<(L2_SZ_VD+5),0x1<<(L2_SZ_VD+5),0xa<<(L2_SZ_VD+5),0x6<<(L2_SZ_VD+5),0x2<<(L2_SZ_VD+5),0xb<<(L2_SZ_VD+5),0x7<<(L2_SZ_VD+5),0x3<<(L2_SZ_VD+5),0xc<<(L2_SZ_VD+5),0x8<<(L2_SZ_VD+5),}
		  #endif
			;
	// Shared by both DIF+DIT:
		int plo[16],phi[13];
		uint64 i64;
		const int *iptr;
		int kk;
	#ifndef USE_SSE2
		int k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf;
	#else
		int po_kperm[16],*po_ptr = &(po_kperm[0]);
	#endif
	#ifdef USE_SSE2
		const uint32 dft13_offs[13] = {0,0x20<<L2_SZ_VD,0x40<<L2_SZ_VD,0x60<<L2_SZ_VD,0x80<<L2_SZ_VD,0xa0<<L2_SZ_VD,0xc0<<L2_SZ_VD,0xe0<<L2_SZ_VD,0x100<<L2_SZ_VD,0x120<<L2_SZ_VD,0x140<<L2_SZ_VD,0x160<<L2_SZ_VD,0x180<<L2_SZ_VD}, *dft13_offptr = &(dft13_offs[0]);
	#endif
		int j,j1,l;
	#ifndef USE_SSE2
		int j2,jt,jp,ntmp;
	#endif
	#ifdef USE_SSE2
		// incr = Carry-chain wts-multipliers recurrence length, which must divide
		// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 13|26|52 for avx512,avx,sse, respectively.
		// But fixed-incr too restrictive here, so 'divide 13|26|52 into pieces' via increment-array whose elts sum to 13|26|52:
		const int *incr,*inc_arr;
	  #ifdef USE_AVX512	// Have no specialized HIACC carry macro in AVX-512 and ARMv8 SIMD, so these get an "goes to 11" in LOACC mode via an incr_hiacc[] array:
		const int incr_long[] = {13}, incr_med[] = {6,7}, incr_short[] = {4,5,4}, incr_hiacc[] = {3,2,3,2,3};
	  #elif defined(USE_AVX)
		const int incr_long[] = {13,13}, incr_med[] = {9,8,9}, incr_short[] = {4,5,4,4,5,4}, incr_hiacc[] = {0};
	  #elif defined(USE_ARM_V8_SIMD)
		const int incr_long[] = {13,13,13,13}, incr_med[] = {9,8,9,9,8,9}, incr_short[] = {4,5,4,4,5,4,4,5,4,4,5,4}, incr_hiacc[] = {3,2,3,2,3,3,2,3,2,3,3,2,3,2,3,3,2,3,2,3};
	  #else
		const int incr_long[] = {13,13,13,13}, incr_med[] = {9,8,9,9,8,9}, incr_short[] = {4,5,4,4,5,4,4,5,4,4,5,4}, incr_hiacc[] = {0};
	  #endif
		// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
		if(USE_SHORT_CY_CHAIN == 0)
			inc_arr = incr_long;
		else if(USE_SHORT_CY_CHAIN == 1)
			inc_arr = incr_med;
		else if(USE_SHORT_CY_CHAIN == 2)
			inc_arr = incr_short;
		else
			inc_arr = incr_hiacc;
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
		double *add0,*add1,*add2,*add3;	/* Addresses into array sections */
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *tmp,*tm1,*tm2;	// Non-static utility ptrs
		vec_dbl *two,/* *one,*sqrt2, */*isrt2, *rad13_const, *max_err, *half_arr, /* *cc0, *ss0, */	/* rad13_const needs 18*16 bytes allocated */
	  #ifndef USE_AVX512
			*sse2_rnd,
	  #endif
			*r00,	// Head of RADIX*vec_cmplx-sized local store #1
			*s1p00,	// Head of RADIX*vec_cmplx-sized local store #2
			*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx
	  #ifndef USE_AVX512
		double dtmp;
	  #endif
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else

	  #ifdef DFT_13_FMA
		const double cc1=  0.88545602565320989590,	/* Real part of exp(i*2*pi/13), the radix-13 fundamental sincos datum	*/
				ss1 =  0.46472317204376854565,	/* Imag part of exp(i*2*pi/13).	*/
				cc2 =  0.56806474673115580252,	/* cos(2u)	*/
				ss2 =  0.82298386589365639457,	/* sin(2u)	*/
				cc3 =  0.12053668025532305336,	/* cos(3u)	*/
				ss3 =  0.99270887409805399279,	/* sin(3u)	*/
				cc4 = -0.35460488704253562594,	/* cos(4u)	*/
				ss4 =  0.93501624268541482344,	/* sin(4u)	*/
				cc5 = -0.74851074817110109861,	/* cos(5u)	*/
				ss5 =  0.66312265824079520240,	/* sin(5u)	*/
				cc6 = -0.97094181742605202714,	/* cos(6u)	*/
				ss6 =  0.23931566428755776718;	/* sin(6u)	*/
	  #else	// Consts for van Buskirk-style tangent DFT in this header:
		#include "radix13.h"
	  #endif
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
		NDIVR >>= 4;			pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
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

	// Shared by both DIF+DIT:
		plo[0x0] =  0; plo[0x1] = p1; plo[0x2] = p2; plo[0x3] = p3;
		plo[0x4] = p4; plo[0x5] = p5; plo[0x6] = p6; plo[0x7] = p7;
		plo[0x8] = p8; plo[0x9] = p9; plo[0xa] = pa; plo[0xb] = pb;
		plo[0xc] = pc; plo[0xd] = pd; plo[0xe] = pe; plo[0xf] = pf;
		l = 0;
		phi[l++] =   0; phi[l++] = p10; phi[l++] = p20; phi[l++] = p30; phi[l++] = p40;
		phi[l++] = p50; phi[l++] = p60; phi[l++] = p70; phi[l++] = p80; phi[l++] = p90;
		phi[l++] = pa0; phi[l++] = pb0; phi[l++] = pc0;

	#ifdef USE_SSE2
		tmp	= r00 = thread_arg->r00;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x1a0;	s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
		tmp += 0x1a0;
		two   = tmp + 0;	// AVX+ versions of radix-8,16,32 twiddleless-DFT macros need consts [2,1,sqrt2,isrt2] quartet laid out thusly
		//one   = tmp + 1;
		//sqrt2 = tmp + 2;
		isrt2 = tmp + 3;
		//cc0   = tmp + 4;
		//ss0   = tmp + 5;
	//	one   = tmp + 6;	Unnamed 1.0 slot to be used in radix-13
	//	two   = tmp + 7;	Unnamed 2.0 slot to be used in radix-13
		rad13_const = tmp + 0x08;	// Needs 17 vec_dbl slots
		tmp += 0x1a;	// Need 8 + 17 = 25 vec_dbl slots for DFT sincos; round up nearest even
	// sc_ptr += 0x36a
	  #ifdef USE_AVX512
		cy = tmp;		tmp += 0x1a;	// RADIX/8 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		//sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
		cy = tmp;		tmp += 0x34;	// RADIX/4 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x390; This is where the value of half_arr_offset208 comes from
		half_arr= tmp + 0x02;	// This table needs 68*SZ_VD bytes in avx mode
	  #else
		cy = tmp;		tmp += 0x68;	// RADIX/2 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x3c4; This is where the value of half_arr_offset208 comes from
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

		sign_mask = (uint64*)(r00 + radix208_creals_in_local_store);
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
		#include "radix208_main_carry_loop.h"

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
