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

#define RADIX 44	// = 0x2c; Use #define rather than const int to ensure it's really a compile-time const in the C sense

#define EPS 1e-10

// Using a flag-value (rather than 'defined or not?' scheme allows us to override the enclosed default at compile time:
#ifdef USE_AVX2
  #ifndef DFT_11_FMA	// To force cyclic-5x5-subconvo-based DFT [160 ADD, 44 MUL] over [10 ADD, 140 FMA]-DFT, use -DDFT_11_FMA=0 at compile time
	#define DFT_11_FMA	1	// Toggle between the 'naive' but good ROE properties DFT and the low-total-arithmetic-opcount
							// but high-ROE van Buskirk-style 'tangent' DFT. Toggle only respected for build modes - that
							// means ones where FMA instructions exist - which allow both options. With FMA, the 'naive'
							// DFT has a comparable opcount to the VB one, since it yields a target-rich environment for FMA.
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
  // Max = (40 [SSE2]; 132 [AVX]), add to (half_arr_offset44 + RADIX) to get value of radix44_creals_in_local_store.
  #ifdef USE_AVX512	// RADIX/8 = 6 fewer carry slots than AVX:
	const int half_arr_offset44 = 0xd0;	// + RADIX = 0xfc; Used for thread local-storage-integrity checking
	const int radix44_creals_in_local_store = 0x198;	// (half_arr_offset44 + RADIX) + 132(=0x84) and round up to nearest multiple of 8
  #elif defined(USE_AVX)
	const int half_arr_offset44 = 0xd5;	// + RADIX = 0x101; Used for thread local-storage-integrity checking
	const int radix44_creals_in_local_store = 0x188;	// (half_arr_offset44 + RADIX) + 132(=0x84) and round up to nearest multiple of 8
  #else
	const int half_arr_offset44 = 0xe0;	// + RADIX = 0x10c; Used for thread local-storage-integrity checking
	const int radix44_creals_in_local_store = 0x134;	// (half_arr_offset44 + RADIX) + 40 and round up to nearest multiple of 4
  #endif

	#include "sse2_macro_gcc64.h"
	#include "radix11_sse_macro.h"

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

/****************/

int radix44_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-44 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-44 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix44_ditN_cy_dif1";
  #if !defined(MULTITHREAD) && defined(USE_SSE2)
	const int pfetch_dist = PFETCH_DIST;
  #endif
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
#if !defined(MULTITHREAD) && defined(USE_SSE2)
	const uint32 dft11_offs[11] = {0,0x2<<L2_SZ_VD,0x4<<L2_SZ_VD,0x6<<L2_SZ_VD,0x8<<L2_SZ_VD,0xa<<L2_SZ_VD,0xc<<L2_SZ_VD,0xe<<L2_SZ_VD,0x10<<L2_SZ_VD,0x12<<L2_SZ_VD,0x14<<L2_SZ_VD}, *dft11_offptr = &(dft11_offs[0]), *ui32_ptr;
	int i0,i1,i2,i3;
#endif
  #ifndef MULTITHREAD
	static uint32 p0123[4];
  #endif
	static double wts_mult[2], inv_mult[2];	// Const wts-multiplier and 2*(its multiplicative inverse)
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
  #endif
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
	int jstart,ntmp;
  #endif
  #ifdef USE_SSE2
	int nbytes;
  #endif
	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40, nsave = 0;
	static int poff[RADIX>>2];	// Store mults of p04 offset for loop control
// FMA-based SIMD or (scalar-double) + (LO_ADD = 1 in masterdefs.h) use these sincos constants:
#if (defined(USE_AVX2) && DFT_11_FMA) || defined(USE_ARM_V8_SIMD) || (!defined(USE_SSE2) && (LO_ADD != 0))
	#warning Using FMA-heavy lo-add 11-DFT
  #if !defined(MULTITHREAD) || (defined(USE_AVX2) && DFT_11_FMA) || defined(USE_ARM_V8_SIMD)
	// FMA based on simple radix-11 DFT implementation, same as LO_ADD - more accurate, and with FMA, faster as well
	const double cc1 =  0.84125353283118116886,	/* Real part of exp(i*2*pi/11), the radix-11 fundamental sincos datum	*/
			ss1 =  0.54064081745559758210,	/* Imag part of exp(i*2*pi/11).	*/
			cc2 =  0.41541501300188642553,	/* cos(2u)	*/
			ss2 =  0.90963199535451837140,	/* sin(2u)	*/
			cc3 = -0.14231483827328514043,	/* cos(3u)	*/
			ss3 =  0.98982144188093273237,	/* sin(3u)	*/
			cc4 = -0.65486073394528506404,	/* cos(4u)	*/
			ss4 =  0.75574957435425828378,	/* sin(4u)	*/
			cc5 = -0.95949297361449738988,	/* cos(5u)	*/
			ss5 =  0.28173255684142969773;	/* sin(5u)	*/
  #endif
#else	// AVX, SSE2 and non-LO_ADD (cf. masterdefs.h) scalar-double builds all use these:
	#warning LO_ADD = 0 defined at compile time ... using low-mul Nussbaumer-style 11-DFT
	// 5 Dec 2014: Change const -> static to enable wrong-way-rounding sensitivity experiments:
	const double a0 = 2.31329240211767848235, /* a0 = (   cq0      -  cq3+  cq2-  cq4)		*/
			a1 =  1.88745388228838373902, /* a1 = (         cq1-  cq3+  cq2-  cq4)		*/
			a2 = -1.41435370755978245393, /* a2 = (-2*cq0-2*cq1+3*cq3-2*cq2+3*cq4)/5	*/
			a3 =  0.08670737584270518028, /* a3 = (-  cq0+  cq1-  cq3+  cq2      )		*/
			a4 = -0.73047075949850706917, /* a4 = (-  cq0+  cq1-  cq3      +  cq4)		*/
			a5 =  0.38639279888589610480, /* a5 = ( 3*cq0-2*cq1+3*cq3-2*cq2-2*cq4)/5	*/
			a6 =  0.51254589567199992361, /* a6 = (            -  cq3+  cq2      )		*/
			a7 =  1.07027574694717148957, /* a7 = (         cq1-  cq3            )		*/
			a8 = -0.55486073394528506404, /* a8 = (-  cq0-  cq1+4*cq3-  cq2-  cq4)/5	*/
			a9 = -1.10000000000000000000, /* a9 = (   cq0+  cq1+  cq3+  cq2+  cq4)/5 - 1*/
			/* Sine terms same form as cos, but cq1 => -sq1 and no -1 on b9 term: */
			b0 =  0.49298012814084233296, /* b0 = (   sq0      -  sq3+  sq2-  sq4)		*/
			b1 = -0.95729268466927362054, /* b1 = (      -  sq1-  sq3+  sq2-  sq4)		*/
			b2 =  0.37415717312460801167, /* b2 = (-2*sq0+2*sq1+3*sq3-2*sq2+3*sq4)/5	*/
			b3 = -1.21620094528344150491, /* b3 = (-  sq0-  sq1-  sq3+  sq2      )		*/
			b4 = -1.92428983032294453955, /* b4 = (-  sq0-  sq1-  sq3      +  sq4)		*/
			b5 =  0.63306543373877589604, /* b5 = ( 3*sq0+2*sq1+3*sq3-2*sq2-2*sq4)/5	*/
			b6 =  0.23407186752667444859, /* b6 = (            -  sq3+  sq2      )		*/
			b7 = -1.66538156970877665518, /* b7 = (      -  sq1-  sq3            )		*/
			b8 =  0.42408709531871829886, /* b8 = (-  sq0+  sq1+4*sq3-  sq2-  sq4)/5	*/
			b9 =  0.33166247903553998491; /* b9 = (   sq0-  sq1+  sq3+  sq2+  sq4)/5	*/
#endif
	static double radix_inv, n2inv;
	double scale, dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex t[RADIX];
#ifndef MULTITHREAD
  #ifndef USE_SSE2
	struct complex *tptr;
  #endif
	int *itmp;	// Pointer into the bjmodn array
  #if defined(USE_AVX) && !defined(USE_AVX512)
	int *itm2;	// Pointer into the bjmodn array
  #endif
#endif
	int err;
	static int first_entry=TRUE;

#if !defined(MULTITHREAD) && defined(USE_SSE2)
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = --|--|11 for avx512,avx,sse, respectively.
	// But even for SSE2, fixed-incr too restrictive here, so 'divide 11 into pieces' via increment-array whose elts sum to 11.
	// For AVX and AVX-512 we just need a single-element array whose 0-term = [0|1] for [AVX|AVX-512] in the HIACC case:
  #ifndef USE_AVX
	const int *incr;
  #endif
	const int *inc_arr;
  #ifdef USE_AVX512	// Have no specialized HIACC carry macro in AVX-512 and ARMv8 SIMD, so these get an "goes to 11" in LOACC mode via an incr_hiacc[] array:
	const int incr_long[] = {1}, incr_med[] = {1}, incr_short[] = {1}, incr_hiacc[] = {1};
  #elif defined(USE_AVX)
	const int incr_long[] = {1}, incr_med[] = {1}, incr_short[] = {1}, incr_hiacc[] = {0};
  #elif defined(USE_ARM_V8_SIMD)
	const int incr_long[] = {6,5}, incr_med[] = {4,3,4}, incr_short[] = {2,2,3,2,2}, incr_hiacc[] = {1,1,1,1,1,1,1,1,1,1,1};
  #else
	const int incr_long[] = {6,5}, incr_med[] = {4,3,4}, incr_short[] = {2,2,3,2,2}, incr_hiacc[] = {0};
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
#endif

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
  #ifndef MULTITHREAD
	int col,co2,co3;
   #ifdef USE_AVX512
	double t0,t1,t2,t3;
	static struct uint32x8 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
   #elif defined(USE_AVX)
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
   #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
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
	double *addr, *add0,*add1,*add2,*add3;	/* Addresses into array sections */
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
	vec_dbl *tmp,*tm2;	// Non-static utility ptrs
  #ifndef MULTITHREAD
	vec_dbl *tm1;	// Non-static utility ptrs
  #endif
	static vec_dbl *two,*one,*five, *ua0,*ua1,*ua2,*ua3,*ua4,*ua5,*ua6,*ua7,*ua8,*ua9, *ub0,*ub1,*ub2,*ub3,*ub4,*ub5,*ub6,*ub7,*ub8,*ub9, *max_err, *sse2_rnd, *half_arr,
		*r00,
	  #ifndef MULTITHREAD
		*s1p00,
	  #endif
		*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx

#endif

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
	static task_control_t   task_control = {NULL, (void*)cy44_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	int m,m2;
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	double *addr;
	int bjmodn[RADIX];
	double rt,it,temp,frac,cy[RADIX];

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_bjmodnini = 0x0,*_bjmodn[RADIX];
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_cy[RADIX];
  #if defined(USE_AVX512) && !defined(USE_THREADS)
	WARN(HERE, "radix44_ditN_cy_dif1: AVX-512 support requires multithreaded build (-DUSE_THREADS); Skipping this leading radix.", "", 1); return(ERR_RADIX0_UNAVAILABLE);
  #endif

	if(!_jhi) {
		_cy[0] = 0x0;	// First of these used as an "already inited consts?" sentinel, must init = 0x0 at same time do so for non-array static ptrs
	}

  #ifdef USE_IMCI512
	WARN(HERE, "radix36_ditN_cy_dif1: No k10m/IMCI-512 support; Skipping this leading radix.", "", 1); return(ERR_RADIX0_UNAVAILABLE);
  #endif

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

	  #ifdef USE_AVX	// AVX LOACC: Make CARRY_8_WAY default here:
		i = 8;
	  #elif defined(USE_SSE2)	// AVX and SSE2 modes use 4-way carry macros
		i = 4;
	  #else	// Scalar-double mode:
		i = 1;
	  #endif
	  {	// Use {} since qt separately def'd elsewhere
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
	  }

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
		// consisting of radix44_creals_in_local_store dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix44_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix44_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x1f) == 0, "sm_ptr not 32-byte aligned!");

	/* Use low 88x2 16-byte slots of sc_arr for temporaries, next 21 for the constants needed by the radix-11 DFT,
	next RADIX/2 = 22 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	#ifdef USE_PTHREAD
		__r0 = sc_ptr;
	#endif
		r00 = sc_ptr;			tmp = r00 + 0x58;
		r00 = r00 + 0x00;
	#ifndef MULTITHREAD
								s1p00 = tmp + 0x00;
	#endif
		tmp += 0x58;	// sc_ptr += 0xb0
		two     = tmp + 0;	// AVX+ versions of various DFT macros need consts [2,1] pair laid out thusly
		one     = tmp + 1;
	//	two     = tmp + 2; Unnamed ptr, used in AVX2 mode to hold 2.0 (and *five holds 1.0) for the radix-11 DFT code
		five    = tmp + 3;
		ua0     = tmp + 4;
		ua1     = tmp + 5;
		ua2     = tmp + 6;
		ua3     = tmp + 7;
		ua4     = tmp + 8;
		ua5     = tmp + 9;
		ua6     = tmp + 0xa;
		ua7     = tmp + 0xb;
		ua8     = tmp + 0xc;
		ua9     = tmp + 0xd;
		ub0     = tmp + 0xe;
		ub1     = tmp + 0xf;
		ub2     = tmp + 0x10;
		ub3     = tmp + 0x11;
		ub4     = tmp + 0x12;
		ub5     = tmp + 0x13;
		ub6     = tmp + 0x14;
		ub7     = tmp + 0x15;
		ub8     = tmp + 0x16;
		ub9     = tmp + 0x17;
		tmp += 0x18;	// sc_ptr += 0xc8
	#ifdef USE_AVX512
		cy = tmp;		tmp += 6;	// sc_ptr += 0xce [0xc8 + RADIX/8 (rounded up) vec_dbl slots for carry sub-array]
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;	// sc_ptr += 0xd0; This is where the value of half_arr_offset44 comes from
	#elif defined(USE_AVX)
		cy = tmp;		tmp += 0x0b;	// sc_ptr += 0xd3 [0xc8 + RADIX/4 vec_dbl slots for carry sub-array]
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0xd5; This is where the value of half_arr_offset44 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy = tmp;		tmp += 0x16;	// sc_ptr += 0xde [0xc8 + RADIX/2 vec_dbl slots for carry sub-array]
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0xe0; This is where the value of half_arr_offset44 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*odd_radix] x 16 for Fermat-mod */
	  #endif

		ASSERT((radix44_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r00) + (20 << L2_SZ_VD), "radix44_creals_in_local_store checksum failed!");
	  #if (defined(USE_AVX2) && DFT_11_FMA) || defined(USE_ARM_V8_SIMD)
	  	/* no-op */
	  #else
		//======= Test the accuracy of the radix-11 trig consts: =========
		struct qfloat qtheta,qt,cq0,cq1,cq2,cq3,cq4,sq1,sq2,sq3,sq4;
		//struct qfloat qs,qfifth,sq0;
		qtheta = qfdiv(Q2PI, i64_to_q((uint64)11));
		qt =          qtheta ;	cq0 = qfcos(qt);	//sq0 = qfsin(qt);
		qt = qfadd(qt,qtheta);	cq1 = qfcos(qt);	sq1 = qfsin(qt);
		qt = qfadd(qt,qtheta);	cq2 = qfcos(qt);	sq2 = qfsin(qt);
		qt = qfadd(qt,qtheta);	cq3 = qfcos(qt);	sq3 = qfsin(qt);
		qt = qfadd(qt,qtheta);	cq4 = qfcos(qt);	sq4 = qfsin(qt);
		//qfifth = qf_rational_quotient((uint64)1, (uint64)5);
	/*
		Sensitivity analysis of the various radix-11 DFT params. Here is the baseline, with all consts rounded-to-nearest:
		./Mlucas -fftlen 352 -iters 10000 -radset 0
		Using complex FFT radices        44        16        16        16

		M7036339 Roundoff warning on iteration      335, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration      632, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration      702, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     1138, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     1597, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     2187, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     2457, maxerr =   0.437500000000
		M7036339 Roundoff warning on iteration     3322, maxerr =   0.460937500000
		M7036339 Roundoff warning on iteration     4207, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     4414, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     4427, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     4436, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     4918, maxerr =   0.437500000000
		M7036339 Roundoff warning on iteration     5438, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     5652, maxerr =   0.421875000000
		M7036339 Roundoff warning on iteration     6806, maxerr =   0.437500000000
		M7036339 Roundoff warning on iteration     7030, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     7379, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     8087, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     8575, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     8660, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     8984, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     9245, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     9371, maxerr =   0.437500000000
		M7036339 Roundoff warning on iteration     9412, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     9498, maxerr =   0.437500000000
		M7036339 Roundoff warning on iteration     9590, maxerr =   0.406250000000
		M7036339 Roundoff warning on iteration     9686, maxerr =   0.406250000000

		28 warnings, 21@0.40625, 1@0.421875, 5@0.4375, 1 > 0.45

		Results of wrong-way-rounding each param in turn and running the same test - we only detail the 'improved' cases:

			#ROE					Summary							Note
			----	---------------------------------------------	----
		a0	3		AvgMaxErr = 0.283830276. MaxErr = 0.457031250	!
		a1	4		AvgMaxErr = 0.279624450. MaxErr = 0.414062500	!!
		a3	19		AvgMaxErr = 0.294811107. MaxErr = 0.437500000
		b0	22		AvgMaxErr = 0.293350167. MaxErr = 0.437500000
		b1	22		AvgMaxErr = 0.294640647. MaxErr = 0.437500000

		Now try all 5 of the above together ... that is actually worse, so there are complex interactions here.
		Try pairwise combos, denoting a0,a1,a3,b0,b1 by indices 0-4, resp., listing results only for combos
		which behave synergistically, i.e. where te result is better than for either term individually:

		0,2		8	AvgMaxErr = 0.283457815. MaxErr = 0.421875000
		0,4		2	AvgMaxErr = 0.284313121. MaxErr = 0.414062500
		1,4		3	AvgMaxErr = 0.280305692. MaxErr = 0.406250000	<< Still a tad worse than 1 [=a1] alone

		024		2	AvgMaxErr = 0.283626017. MaxErr = 0.406250000
		124		5	AvgMaxErr = 0.280089419. MaxErr = 0.406250000	<< #warns greater than 024 but AME is less.

		0124 is again awful.

		Try both of the above triplets on the next-larger prime, 7036343:
		024		2	AvgMaxErr = 0.283818714. MaxErr = 0.406250000
		124		4	AvgMaxErr = 0.279806160. MaxErr = 0.437500000
		..and the next-next-larger-prime, 7036361:
		024		8	AvgMaxErr = 0.283602988. MaxErr = 0.437500000	<< hmm, tough call - but blowup of 024 err stats here shows how noisy
		124		3	AvgMaxErr = 0.278681073. MaxErr = 0.437500000	<< these data are, only consistent trend is in the AME, 124 better there.
		Lastly, try a prime 0.5% larger than above, p = 7071521:
		024		109 warns (18@0.4375, 3 > 0.4375) thru iter 8705 (12.5 warn/Kiter), at whch point get fatal err of 0.484375
		124		 82 warns (18@0.4375, 2 > 0.4375) thru iter 7123 (11.5 warn/Kiter), at whch point get fatal err of 0.500000

	*/
		// Propagate to the VEC_DBL_INITs below just the ones we finally settle on of the full set in the following #if 0 block:
	   #if 0
		//================================================================
		// *** Cosine terms: ***
		qt = qfsub(qfadd(cq0,cq2),qfadd(cq3,cq4));	dtmp = qfdbl(qt);	ASSERT(dtmp == a0, "a0");	/* a0 = (   cq0      -  cq3+  cq2-  cq4) */
		qt = qfsub(qfadd(cq1,cq2),qfadd(cq3,cq4));	dtmp = qfdbl(qt);	ASSERT(dtmp == a1, "a1");	/* a1 = (         cq1-  cq3+  cq2-  cq4) */
		qt = qfsub(qfadd(cq1,cq2),qfadd(cq0,cq3));	dtmp = qfdbl(qt);	ASSERT(dtmp == a3, "a3");	/* a3 = (-  cq0+  cq1-  cq3+  cq2      ) */
		qt = qfsub(qfadd(cq1,cq4),qfadd(cq0,cq3));	dtmp = qfdbl(qt);	ASSERT(dtmp == a4, "a4");	/* a4 = (-  cq0+  cq1-  cq3      +  cq4) */
		qt = qfsub(cq2,cq3);						dtmp = qfdbl(qt);	ASSERT(dtmp == a6, "a6");	/* a6 = (            -  cq3+  cq2      ) */
		qt = qfsub(cq1,cq3);						dtmp = qfdbl(qt);	ASSERT(dtmp == a7, "a7");	/* a7 = (         cq1-  cq3            ) */
		qt = qfmul(qfifth, qfsub( qfmul_pow2(cq3,2), qfadd(qfadd(cq0,cq1),qfadd(cq2,cq4)) ));	dtmp = qfdbl(qt);	ASSERT(dtmp == a8, "a8");	/* a8 = (-  cq0-  cq1+4*cq3-  cq2-  cq4)/5 */
		qt = qfsub( qfmul(qfifth, qfadd( cq0 , qfadd(qfadd(cq1,cq2),qfadd(cq3,cq4)) )), QONE);	dtmp = qfdbl(qt);	ASSERT(dtmp == a9, "a9");	/* a9 = (   cq0+  cq1+  cq3+  cq2+  cq4)/5 - 1 */
		qs = qfadd(qfadd(cq0,cq1), cq2);	qs = qfmul_pow2(qs,1);	// 2*(cq0+cq1+cq2)
		qt = qfadd(cq3,cq4);	qt = qfadd(qt, qfmul_pow2(qt,1));	// 3*(cq3+cq4)
		qt = qfmul(qfifth, qfsub(qt,qs));	dtmp = qfdbl(qt);	ASSERT(dtmp == a2, "a2");	/* a2 = (-2*cq0-2*cq1+3*cq3-2*cq2+3*cq4)/5 */
		qs = qfadd(qfadd(cq4,cq1), cq2);	qs = qfmul_pow2(qs,1);	// 2*(cq4+cq1+cq2)
		qt = qfadd(cq3,cq0);	qt = qfadd(qt, qfmul_pow2(qt,1));	// 3*(cq3+cq0)
		qt = qfmul(qfifth, qfsub(qt,qs));	dtmp = qfdbl(qt);	ASSERT(dtmp == a5, "a5");	/* a5 = ( 3*cq0-2*cq1+3*cq3-2*cq2-2*cq4)/5 */
		// *** Sine terms: ***
		qt = qfsub(qfadd(sq0,sq2),qfadd(sq3,sq4));	dtmp = qfdbl(qt);	ASSERT(dtmp == b0, "b0");	/* b0 = (   sq0      -  sq3+  sq2-  sq4) */
		qt = qfsub(qfsub(sq2,sq1),qfadd(sq3,sq4));	dtmp = qfdbl(qt);	ASSERT(dtmp == b1, "b1");	/* b1 = (        -sq1-  sq3+  sq2-  sq4) */
		qt = qfsub(qfsub(sq2,sq1),qfadd(sq0,sq3));	dtmp = qfdbl(qt);	ASSERT(dtmp == b3, "b3");	/* b3 = (-  sq0-  sq1-  sq3+  sq2      ) */
		qt = qfsub(qfsub(sq4,sq1),qfadd(sq0,sq3));	dtmp = qfdbl(qt);	ASSERT(dtmp == b4, "b4");	/* b4 = (-  sq0-  sq1-  sq3      +  sq4) */
		qt = qfsub(sq2,sq3);						dtmp = qfdbl(qt);	ASSERT(dtmp == b6, "b6");	/* b6 = (            -  sq3+  sq2      ) */
		qt = qfneg(qfadd(sq1,sq3));					dtmp = qfdbl(qt);	ASSERT(dtmp == b7, "b7");	/* b7 = (        -sq1-  sq3            ) */
		qt = qfmul(qfifth, qfsub( qfmul_pow2(sq3,2), qfadd(qfsub(sq0,sq1),qfadd(sq2,sq4)) ));	dtmp = qfdbl(qt);	ASSERT(dtmp == b8, "b8");	/* b8 = (-  sq0+  sq1+4*sq3-  sq2-  sq4)/5 */
		qt = qfmul(qfifth, qfadd( sq0 , qfadd(qfsub(sq2,sq1),qfadd(sq3,sq4)) ));	dtmp = qfdbl(qt);	ASSERT(dtmp == b9, "b9");	/* b9 = (   sq0-  sq1+  sq3+  sq2+  sq4)/5 - 1 */
		qs = qfadd(qfsub(sq0,sq1), sq2);	qs = qfmul_pow2(qs,1);	// 2*(sq0-sq1+sq2)
		qt = qfadd(sq3,sq4);	qt = qfadd(qt, qfmul_pow2(qt,1));	// 3*(sq3+sq4)
		qt = qfmul(qfifth, qfsub(qt,qs));	dtmp = qfdbl(qt);	ASSERT(dtmp == b2, "b2");	/* b2 = (-2*sq0+2*sq1+3*sq3-2*sq2+3*sq4)/5 */
		qs = qfadd(qfsub(sq4,sq1), sq2);	qs = qfmul_pow2(qs,1);	// 2*(sq4-sq1+sq2)
		qt = qfadd(sq3,sq0);	qt = qfadd(qt, qfmul_pow2(qt,1));	// 3*(sq3+sq0)
		qt = qfmul(qfifth, qfsub(qt,qs));	dtmp = qfdbl(qt);	ASSERT(dtmp == b5, "b5");	/* b5 = ( 3*sq0+2*sq1+3*sq3-2*sq2-2*sq4)/5 */
		//================================================================
	   #endif
	  #endif
		/* These remain fixed: */
		VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );
	  #if (defined(USE_AVX2) && DFT_11_FMA) || defined(USE_ARM_V8_SIMD)	// FMA version based on simple radix-11 DFT implementation, same as LO_ADD
		tmp = five-1;
		VEC_DBL_INIT(tmp ,2.0);		VEC_DBL_INIT(five,1.0);	// *five-1,*five used for [2.0,1.0] here
		VEC_DBL_INIT(ua0 ,cc1);		VEC_DBL_INIT(ub0 ,0.0);	// upper 10 slots unused here; init = 0
		VEC_DBL_INIT(ua1 ,cc2);		VEC_DBL_INIT(ub1 ,0.0);
		VEC_DBL_INIT(ua2 ,cc3);		VEC_DBL_INIT(ub2 ,0.0);
		VEC_DBL_INIT(ua3 ,cc4);		VEC_DBL_INIT(ub3 ,0.0);
		VEC_DBL_INIT(ua4 ,cc5);		VEC_DBL_INIT(ub4 ,0.0);
		VEC_DBL_INIT(ua5 ,ss1);		VEC_DBL_INIT(ub5 ,0.0);
		VEC_DBL_INIT(ua6 ,ss2);		VEC_DBL_INIT(ub6 ,0.0);
		VEC_DBL_INIT(ua7 ,ss3);		VEC_DBL_INIT(ub7 ,0.0);
		VEC_DBL_INIT(ua8 ,ss4);		VEC_DBL_INIT(ub8 ,0.0);
		VEC_DBL_INIT(ua9 ,ss5);		VEC_DBL_INIT(ub9 ,0.0);
	  #else
		tmp = five-1;
		VEC_DBL_INIT(tmp ,2.0);		VEC_DBL_INIT(five,5.0);	// *five-1 used for 2.0 in SSE2 mode [The fused SSE2_RADIX44_DIT_NOTWIDDLE macro needs the 2.0 there]
		VEC_DBL_INIT(ua0 , a0);		VEC_DBL_INIT(ub0 , b0);
		VEC_DBL_INIT(ua1 , a1);		VEC_DBL_INIT(ub1 , b1);
		VEC_DBL_INIT(ua2 , a2);		VEC_DBL_INIT(ub2 , b2);
		VEC_DBL_INIT(ua3 , a3);		VEC_DBL_INIT(ub3 , b3);
		VEC_DBL_INIT(ua4 , a4);		VEC_DBL_INIT(ub4 , b4);
		VEC_DBL_INIT(ua5 , a5);		VEC_DBL_INIT(ub5 , b5);
		VEC_DBL_INIT(ua6 , a6);		VEC_DBL_INIT(ub6 , b6);
		VEC_DBL_INIT(ua7 , a7);		VEC_DBL_INIT(ub7 , b7);
		VEC_DBL_INIT(ua8 , a8);		VEC_DBL_INIT(ub8 , b8);
		VEC_DBL_INIT(ua9 , a9);		VEC_DBL_INIT(ub9 ,-b9);	/* Flip sign to simplify code re-use in radix-11 SSE2 macro */
		qt = qfsub(qfadd(cq1,cq2),qfadd(cq3,cq4));	dtmp = qfdbl_wrong_way_rnd(qt);	VEC_DBL_INIT(ua1,dtmp);	// wrong-way-rounding of a1
		qt = qfsub(qfadd(cq1,cq2),qfadd(cq0,cq3));	dtmp = qfdbl_wrong_way_rnd(qt);	VEC_DBL_INIT(ua3,dtmp);	// wrong-way-rounding of a3
		qt = qfsub(qfsub(sq2,sq1),qfadd(sq3,sq4));	dtmp = qfdbl_wrong_way_rnd(qt);	VEC_DBL_INIT(ub1,dtmp);	// wrong-way-rounding of b1
	  #endif
		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		sse2_rnd->d0 = base[0]; sse2_rnd->d1 = baseinv[1]; sse2_rnd->d2 = wts_mult[1]; sse2_rnd->d3 = inv_mult[0];
	  #else
		VEC_DBL_INIT(sse2_rnd, crnd);
	  #endif

		// Propagate the above consts to the remaining threads:
		nbytes = (intptr_t)cy - (intptr_t)two;
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
	  #ifndef MULTITHREAD
		n_minus_sil   = (struct uint32x8 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x8 *)sse_n + 2;
		sinwt         = (struct uint32x8 *)sse_n + 3;
		sinwtm1       = (struct uint32x8 *)sse_n + 4;
	  #endif
		nbytes += 128;
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
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;
		p28 = p24 + p04;
		p32 = p28 + p04;
		p36 = p32 + p04;
		p40 = p36 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p36 = p36 + ( (p36 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );

	#ifndef MULTITHREAD
		p0123[0] = 0; p0123[1] = p01; p0123[2] = p02; p0123[3] = p03;
	#endif

		poff[0x0] =   0; poff[0x1] = p04; poff[0x2] = p08; poff[0x3] = p12;
		poff[0x4] = p16; poff[0x5] = p20; poff[0x6] = p24; poff[0x7] = p28;
		poff[0x8] = p32; poff[0x9] = p36; poff[0xa] = p40;

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

/*...The radix-44 final DIT pass is here.	*/

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
		#include "radix44_main_carry_loop.h"

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
		ASSERT(0x0 == cy44_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
	ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
	ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(0 == mlucas_nanosleep(&ns_time), "nanosleep fail!");
	}
//	printf("radix44_ditN_cy_dif1 end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

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
				jt = j1 + poff[l];	// poff[] = p04,p08,...,p56
				a[jt    ] *= radix_inv;
				a[jt+p01] *= radix_inv;
				a[jt+p02] *= radix_inv;
				a[jt+p03] *= radix_inv;
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

void radix44_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-44 complex DIF FFT pass on the data in the length-N real vector A.
*/
	int j,jt,jp,j1,j2;
	static int NDIVR,first_entry=TRUE,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40;
#if LO_ADD
	const double cc1 =  0.84125353283118116886,	/* Real part of exp(i*2*pi/11), the radix-11 fundamental sincos datum	*/
			ss1 =  0.54064081745559758210,	/* Imag part of exp(i*2*pi/11).	*/
			cc2 =  0.41541501300188642553,	/* cos(2u)	*/
			ss2 =  0.90963199535451837140,	/* sin(2u)	*/
			cc3 = -0.14231483827328514043,	/* cos(3u)	*/
			ss3 =  0.98982144188093273237,	/* sin(3u)	*/
			cc4 = -0.65486073394528506404,	/* cos(4u)	*/
			ss4 =  0.75574957435425828378,	/* sin(4u)	*/
			cc5 = -0.95949297361449738988,	/* cos(5u)	*/
			ss5 =  0.28173255684142969773;	/* sin(5u)	*/
#else
	const double a0 =  2.31329240211767848235, /* a0 = (   cq0      -  cq3+  cq2-  cq4)		*/
			a1 =  1.88745388228838373902, /* a1 = (         cq1-  cq3+  cq2-  cq4)		*/
			a2 = -1.41435370755978245393, /* a2 = (-2*cq0-2*cq1+3*cq3-2*cq2+3*cq4)/5	*/
			a3 =  0.08670737584270518028, /* a3 = (-  cq0+  cq1-  cq3+  cq2      )		*/
			a4 = -0.73047075949850706917, /* a4 = (-  cq0+  cq1-  cq3      +  cq4)		*/
			a5 =  0.38639279888589610480, /* a5 = ( 3*cq0-2*cq1+3*cq3-2*cq2-2*cq4)/5	*/
			a6 =  0.51254589567199992361, /* a6 = (            -  cq3+  cq2      )		*/
			a7 =  1.07027574694717148957, /* a7 = (         cq1-  cq3            )		*/
			a8 = -0.55486073394528506404, /* a8 = (-  cq0-  cq1+4*cq3-  cq2-  cq4)/5	*/
			a9 = -1.10000000000000000000, /* a9 = (   cq0+  cq1+  cq3+  cq2+  cq4)/5 - 1*/

			b0 =  0.49298012814084233296, /* b0 = (   sq0      -  sq3+  sq2-  sq4)		*/
			b1 = -0.95729268466927362054, /* b1 = (      -  sq1-  sq3+  sq2-  sq4)		*/
			b2 =  0.37415717312460801167, /* b2 = (-2*sq0+2*sq1+3*sq3-2*sq2+3*sq4)/5	*/
			b3 = -1.21620094528344150491, /* b3 = (-  sq0-  sq1-  sq3+  sq2      )		*/
			b4 = -1.92428983032294453955, /* b4 = (-  sq0-  sq1-  sq3      +  sq4)		*/
			b5 =  0.63306543373877589604, /* b5 = ( 3*sq0+2*sq1+3*sq3-2*sq2-2*sq4)/5	*/
			b6 =  0.23407186752667444859, /* b6 = (            -  sq3+  sq2      )		*/
			b7 = -1.66538156970877665518, /* b7 = (      -  sq1-  sq3            )		*/
			b8 =  0.42408709531871829886, /* b8 = (-  sq0+  sq1+4*sq3-  sq2-  sq4)/5	*/
			b9 =  0.33166247903553998491; /* b9 = (   sq0-  sq1+  sq3+  sq2+  sq4)/5	*/
#endif
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

		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;
		p28 = p24 + p04;
		p32 = p28 + p04;
		p36 = p32 + p04;
		p40 = p36 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p36 = p36 + ( (p36 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-44 pass is here.	*/

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
	/*...gather the needed data (44 64-bit complex, i.e. 88 64-bit reals) and do 4 radix-11 transforms...*/
	/* Twiddleless version arranges 4 sets of radix-11 DFT inputs as follows:
		0 in upper left corner, decrement 4 horizontally and 11 vertically (mod 44):

		The required output permutation is i*12 (mod 44) for i = 0-43:
			Base-10 input-index offsets:					Hex input-index offsets:
			00,40,36,32,28,24,20,16,12,08,04			00,28,24,20,1c,18,14,10,0c,08,04
			33,29,25,21,17,13,09,05,01,41,37			21,1d,19,15,11,0d,09,05,01,29,25
			22,18,14,10,06,02,42,38,34,30,26			16,12,0e,0a,06,02,2a,26,22,1e,1a
			11,07,03,43,39,35,31,27,23,19,15			0b,07,03,2b,27,23,1f,1b,17,13,0f

	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF] to properly permute the radix-4 DFT outputs.
	*/
		tptr = t;
	#if LO_ADD
		jt = j1    ; jp = j2    ;	RADIX_11_DFT_BASIC(a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	tptr += 11;
		jt = j1+p01; jp = j2+p01;	RADIX_11_DFT_BASIC(a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	tptr += 11;
		jt = j1+p02; jp = j2+p02;	RADIX_11_DFT_BASIC(a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	tptr += 11;
		jt = j1+p03; jp = j2+p03;	RADIX_11_DFT_BASIC(a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);
	#else
		jt = j1    ; jp = j2    ;	RADIX_11_DFT(a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);	tptr += 11;
		jt = j1+p01; jp = j2+p01;	RADIX_11_DFT(a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);	tptr += 11;
		jt = j1+p02; jp = j2+p02;	RADIX_11_DFT(a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);	tptr += 11;
		jt = j1+p03; jp = j2+p03;	RADIX_11_DFT(a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);
	#endif
		/*...and now do 11 radix-4 transforms...*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p40; jp = j2+p40;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p36; jp = j2+p36;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	tptr++;
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	tptr++;
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	tptr++;
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	tptr++;
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	tptr++;
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIF(tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
		/* Totals: 4*radix11 + 11*radix04 = 4*(152 FADD, 44 FMUL)	+ 11*(16 FADD, 0 FMUL) = 784 FADD, 176 FMUL	*/
	}
}
/*
Unpermuted DIF outputs in left column, and the output index where each needs to go in []
in order to match the bit-reversed matrix-multiply DFT outputs in the right set of columns:

   0 [   0] 236.0000000000   189.0000000000  (   0)  236.0000000000   189.0000000000
   1 [   1] -32.0000000000   -27.0000000000  (  22)  -32.0000000000   -27.0000000000
   2 [   2] -10.0000000000    -5.0000000000  (  11)  -10.0000000000    -5.0000000000
   3 [   3] -18.0000000000   -25.0000000000  (  33)  -18.0000000000   -25.0000000000

   4 [  41]  11.1512215466    24.0242163196  (   1)   -3.5575051052    -7.5004478474
   5 [  40]  -3.6545658269    18.0551079573  (  23)    3.1041298533    -7.8963424647
   6 [  43]   7.4820065980   -12.0231995596  (  12)   14.9388279007   -12.3752085181
   7 [  42]  16.9375513781    23.8501660484  (  34)    8.8463140932   -16.0678583967

   8 [  38]  39.0936733850   -32.3447548001  (   2)   -9.4839628566     6.8049779883
   9 [  39]  -0.5216089174    15.2315827674  (  24)   29.2378241518    33.2630685362
  10 [  37] -48.7357676843    10.4812779061  (  13)   22.0526961791    -6.2994315330
  11 [  36]  -3.7958386977    14.7194188005  (  35)  -12.8767451770    -6.2427936444

  12 [  32] -10.8100095734    21.3984549366  (   3)  -19.1813630782   -13.9272285822
  13 [  33] -19.1224413133   -21.5236610012  (  25)  -10.4370837863   -26.2672878420
  14 [  34]   1.1660991398    28.0446190599  (  14)   28.5057795698   -23.4048812091
  15 [  35]  -1.1477307972    -3.9728027190  (  36)   13.2967408741    10.4272914314

  16 [  31]  -4.2498272043    10.6726404285  (   4)  -15.5870869025   -44.7914305956
  17 [  30]  -0.9626434619     0.5229653520  (  26)   -1.0423591098    40.9542963604
  18 [  28] -25.0400889167   -18.4141433297  (  15)    6.7261054808    -4.8184504585
  19 [  29] -28.7125623109   -32.1535503331  (  37)   17.0043942516    -9.8950027594

  20 [  25]  -4.1708965255     8.4490527718  (   5)  -24.4303401948   -18.0065547227
  21 [  24]  12.2156779504   -12.9675770627  (  27)    1.4536458478   -23.1558811910
  22 [  27] -14.0980220633   -24.0478577359  (  16)  -33.9004676525    34.2766694898
  23 [  26]  12.0864190834    -8.4745070923  (  38)    6.2198098725    -2.6049527558

  24 [  22] -33.9004676525    34.2766694898  (   6)   12.2156779504   -12.9675770627
  25 [  23]   6.2198098725    -2.6049527558  (  28)   -4.1708965255     8.4490527718
  26 [  21]   1.4536458478   -23.1558811910  (  17)   12.0864190834    -8.4745070923
  27 [  20] -24.4303401948   -18.0065547227  (  39)  -14.0980220633   -24.0478577359

  28 [  16] -15.5870869025   -44.7914305956  (   7)  -25.0400889167   -18.4141433297
  29 [  17]  -1.0423591098    40.9542963604  (  29)  -28.7125623109   -32.1535503331
  30 [  18]   6.7261054808    -4.8184504585  (  18)   -0.9626434619     0.5229653520
  31 [  19]  17.0043942516    -9.8950027594  (  40)   -4.2498272043    10.6726404285

  32 [  15]  13.2967408741    10.4272914314  (   8)  -10.8100095734    21.3984549366
  33 [  14]  28.5057795698   -23.4048812091  (  30)  -19.1224413133   -21.5236610012
  34 [  12] -19.1813630782   -13.9272285822  (  19)    1.1660991398    28.0446190599
  35 [  13] -10.4370837863   -26.2672878420  (  41)   -1.1477307972    -3.9728027190

  36 [   9]  29.2378241518    33.2630685362  (   9)   -3.7958386977    14.7194188005
  37 [   8]  -9.4839628566     6.8049779883  (  31)  -48.7357676843    10.4812779061
  38 [  11] -12.8767451770    -6.2427936444  (  20)   39.0936733850   -32.3447548001
  39 [  10]  22.0526961791    -6.2994315330  (  42)   -0.5216089174    15.2315827674

  40 [   6]  14.9388279007   -12.3752085181  (  10)   -3.6545658269    18.0551079573
  41 [   7]   8.8463140932   -16.0678583967  (  32)   11.1512215466    24.0242163196
  42 [   5]   3.1041298533    -7.8963424647  (  21)   16.9375513781    23.8501660484
  43 [   4]  -3.5575051052    -7.5004478474  (  43)    7.4820065980   -12.0231995596
*/
/***************/

void radix44_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-44 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix44_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
*/
	int j,jt,jp,j1,j2;
	static int NDIVR,first_entry=TRUE,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40;
#if LO_ADD
	const double cc1 =  0.84125353283118116886,	/* Real part of exp(i*2*pi/11), the radix-11 fundamental sincos datum	*/
			ss1 =  0.54064081745559758210,	/* Imag part of exp(i*2*pi/11).	*/
			cc2 =  0.41541501300188642553,	/* cos(2u)	*/
			ss2 =  0.90963199535451837140,	/* sin(2u)	*/
			cc3 = -0.14231483827328514043,	/* cos(3u)	*/
			ss3 =  0.98982144188093273237,	/* sin(3u)	*/
			cc4 = -0.65486073394528506404,	/* cos(4u)	*/
			ss4 =  0.75574957435425828378,	/* sin(4u)	*/
			cc5 = -0.95949297361449738988,	/* cos(5u)	*/
			ss5 =  0.28173255684142969773;	/* sin(5u)	*/
#else
	const double a0 =  2.31329240211767848235, /* a0 = (   cq0      -  cq3+  cq2-  cq4)		*/
			a1 =  1.88745388228838373902, /* a1 = (         cq1-  cq3+  cq2-  cq4)		*/
			a2 = -1.41435370755978245393, /* a2 = (-2*cq0-2*cq1+3*cq3-2*cq2+3*cq4)/5	*/
			a3 =  0.08670737584270518028, /* a3 = (-  cq0+  cq1-  cq3+  cq2      )		*/
			a4 = -0.73047075949850706917, /* a4 = (-  cq0+  cq1-  cq3      +  cq4)		*/
			a5 =  0.38639279888589610480, /* a5 = ( 3*cq0-2*cq1+3*cq3-2*cq2-2*cq4)/5	*/
			a6 =  0.51254589567199992361, /* a6 = (            -  cq3+  cq2      )		*/
			a7 =  1.07027574694717148957, /* a7 = (         cq1-  cq3            )		*/
			a8 = -0.55486073394528506404, /* a8 = (-  cq0-  cq1+4*cq3-  cq2-  cq4)/5	*/
			a9 = -1.10000000000000000000, /* a9 = (   cq0+  cq1+  cq3+  cq2+  cq4)/5 - 1*/

			b0 =  0.49298012814084233296, /* b0 = (   sq0      -  sq3+  sq2-  sq4)		*/
			b1 = -0.95729268466927362054, /* b1 = (      -  sq1-  sq3+  sq2-  sq4)		*/
			b2 =  0.37415717312460801167, /* b2 = (-2*sq0+2*sq1+3*sq3-2*sq2+3*sq4)/5	*/
			b3 = -1.21620094528344150491, /* b3 = (-  sq0-  sq1-  sq3+  sq2      )		*/
			b4 = -1.92428983032294453955, /* b4 = (-  sq0-  sq1-  sq3      +  sq4)		*/
			b5 =  0.63306543373877589604, /* b5 = ( 3*sq0+2*sq1+3*sq3-2*sq2-2*sq4)/5	*/
			b6 =  0.23407186752667444859, /* b6 = (            -  sq3+  sq2      )		*/
			b7 = -1.66538156970877665518, /* b7 = (      -  sq1-  sq3            )		*/
			b8 =  0.42408709531871829886, /* b8 = (-  sq0+  sq1+4*sq3-  sq2-  sq4)/5	*/
			b9 =  0.33166247903553998491; /* b9 = (   sq0-  sq1+  sq3+  sq2+  sq4)/5	*/
#endif
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

		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;
		p28 = p24 + p04;
		p32 = p28 + p04;
		p36 = p32 + p04;
		p40 = p36 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p36 = p36 + ( (p36 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-44 pass is here.	*/

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
		Twiddleless version requires us to swap inputs as follows:
		indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43
			  -> 0,40,36,32,28,24,20,16,12, 8, 4|33,29,25,21,17,13, 9, 5, 1,41,37|22,18,14,10, 6, 2,42,38,34,30,26|11, 7, 3,43,39,35,31,27,23,19,15

		I.e. start out with first quartet of indices {0,11,22,33}, permute those according to
		(0,11,22,33}*43%44 = {0,33,22,11), then each is head of a length-11 list of indices with decrement 4.

		Remember, inputs to DIT are bit-reversed, so
		a[ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43] contain
		x[ 0,22,11,33, 1,23,12,34, 2,24,13,35, 3,25,14,36, 4,26,15,37, 5,27,16,38, 6,28,17,39, 7,29,18,40, 8,30,19,41, 9,31,20,42,10,32,21,43], which get swapped to
	(Look at first nontrivial one, i.e. x[1 -> 40] ... in terms of a[] this translates to a[4 -> 31])
		x[ 0,22,33,11,40,18,29, 7,36,14,25, 3,32,10,21,43,28, 6,17,39,24, 2,13,35,20,42, 9,31,16,38, 5,27,12,34, 1,23, 8,30,41,19, 4,26,37,15], which means the a-indices get swapped as
		a[ 0, 1, 3, 2,31,30,29,28,15,14,13,12,41,40,42,43,25,24,26,27, 9, 8,10,11,38,39,36,37,22,23,20,21, 6, 7, 4, 5,32,33,35,34,16,17,19,18].
	*/
	/*...gather the needed data (44 64-bit complex, i.e. 88 64-bit reals) and do 11 radix-4 transforms...*/
					/*                                   inputs                                   */ /*              outputs              */
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p40; jp = j2+p40;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p36; jp = j2+p36;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);	tptr++;
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02], tptr->re,tptr->im,(tptr+11)->re,(tptr+11)->im,(tptr+22)->re,(tptr+22)->im,(tptr+33)->re,(tptr+33)->im, rt,it);

		/*...and now do 4 radix-11 transforms.
		The required output permutation is i*12 (mod 44) for i = 0-43:
			Base-10 output-index offsets:					Hex output-index offsets:
			00,12,24,36,04,16,28,40,08,20,32			00,0c,18,24,04,10,1c,28,08,14,20
			11,23,35,03,15,27,39,07,19,31,43			0b,17,23,03,0f,1b,27,07,13,1f,2b
			22,34,02,14,26,38,06,18,30,42,10			16,22,02,0e,1a,26,06,12,1e,2a,0a
			33,01,13,25,37,05,17,29,41,09,21			21,01,0d,19,25,05,11,1d,29,09,15
		*/
		tptr = t;
	#if LO_ADD
		jt = j1    ; jp = j2    ;	RADIX_11_DFT_BASIC(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	tptr += 11;
		jt = j1+p03; jp = j2+p03;	RADIX_11_DFT_BASIC(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	tptr += 11;
		jt = j1+p02; jp = j2+p02;	RADIX_11_DFT_BASIC(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);	tptr += 11;
		jt = j1+p01; jp = j2+p01;	RADIX_11_DFT_BASIC(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20], cc1,cc2,cc3,cc4,cc5,ss1,ss2,ss3,ss4,ss5);
	#else
		jt = j1    ; jp = j2    ;	RADIX_11_DFT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32], a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);	tptr += 11;
		jt = j1+p03; jp = j2+p03;	RADIX_11_DFT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40], a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);	tptr += 11;
		jt = j1+p02; jp = j2+p02;	RADIX_11_DFT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a[jt+p20],a[jp+p20],a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08], a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);	tptr += 11;
		jt = j1+p01; jp = j2+p01;	RADIX_11_DFT(tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im, a[jt+p32],a[jp+p32],a[jt    ],a[jp    ],a[jt+p12],a[jp+p12],a[jt+p24],a[jp+p24],a[jt+p36],a[jp+p36],a[jt+p04],a[jp+p04],a[jt+p16],a[jp+p16],a[jt+p28],a[jp+p28],a[jt+p40],a[jp+p40],a[jt+p08],a[jp+p08],a[jt+p20],a[jp+p20], a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9);
	#endif
	}
}
/*
Unpermuted DIF outputs in left column, and the output index where each needs to go in []
in order to match the matrix-multiply DFT outputs in the right set of columns:

   0 [00]   236.0000000000  -189.0000000000     236.0000000000  -189.0000000000
   1 [11]   -10.0000000000     5.0000000000      -3.5575051052     7.5004478474
   2 [22]   -32.0000000000    27.0000000000      -9.4839628566    -6.8049779883
   3 [33]   -18.0000000000    25.0000000000     -19.1813630782    13.9272285822

   4 [12]    14.9388279007    12.3752085181     -15.5870869025    44.7914305956
   5 [23]     3.1041298533     7.8963424647     -24.4303401948    18.0065547227
   6 [34]     8.8463140932    16.0678583967      12.2156779504    12.9675770627
   7 [01]    -3.5575051052     7.5004478474     -25.0400889167    18.4141433297

   8 [24]    29.2378241518   -33.2630685362     -10.8100095734   -21.3984549366
   9 [35]   -12.8767451770     6.2427936444      -3.7958386977   -14.7194188005
  10 [02]    -9.4839628566    -6.8049779883      -3.6545658269   -18.0551079573
  11 [13]    22.0526961791     6.2994315330     -10.0000000000     5.0000000000

  12 [36]    13.2967408741   -10.4272914314      14.9388279007    12.3752085181
  13 [03]   -19.1813630782    13.9272285822      22.0526961791     6.2994315330
  14 [14]    28.5057795698    23.4048812091      28.5057795698    23.4048812091
  15 [25]   -10.4370837863    26.2672878420       6.7261054808     4.8184504585

  16 [04]   -15.5870869025    44.7914305956     -33.9004676525   -34.2766694898
  17 [15]     6.7261054808     4.8184504585      12.0864190834     8.4745070923
  18 [26]    -1.0423591098   -40.9542963604      -0.9626434619    -0.5229653520
  19 [37]    17.0043942516     9.8950027594       1.1660991398   -28.0446190599

  20 [16]   -33.9004676525   -34.2766694898      39.0936733850    32.3447548001
  21 [27]     1.4536458478    23.1558811910      16.9375513781   -23.8501660484
  22 [38]     6.2198098725     2.6049527558     -32.0000000000    27.0000000000
  23 [05]   -24.4303401948    18.0065547227       3.1041298533     7.8963424647

  24 [28]    -4.1708965255    -8.4490527718      29.2378241518   -33.2630685362
  25 [39]   -14.0980220633    24.0478577359     -10.4370837863    26.2672878420
  26 [06]    12.2156779504    12.9675770627      -1.0423591098   -40.9542963604
  27 [17]    12.0864190834     8.4745070923       1.4536458478    23.1558811910

  28 [40]    -4.2498272043   -10.6726404285      -4.1708965255    -8.4490527718
  29 [07]   -25.0400889167    18.4141433297     -28.7125623109    32.1535503331
  30 [18]    -0.9626434619    -0.5229653520     -19.1224413133    21.5236610012
  31 [29]   -28.7125623109    32.1535503331     -48.7357676843   -10.4812779061

  32 [08]   -10.8100095734   -21.3984549366      11.1512215466   -24.0242163196
  33 [19]     1.1660991398   -28.0446190599     -18.0000000000    25.0000000000
  34 [30]   -19.1224413133    21.5236610012       8.8463140932    16.0678583967
  35 [41]    -1.1477307972     3.9728027190     -12.8767451770     6.2427936444

  36 [20]    39.0936733850    32.3447548001      13.2967408741   -10.4272914314
  37 [31]   -48.7357676843   -10.4812779061      17.0043942516     9.8950027594
  38 [42]    -0.5216089174   -15.2315827674       6.2198098725     2.6049527558
  39 [09]    -3.7958386977   -14.7194188005     -14.0980220633    24.0478577359

  40 [32]    11.1512215466   -24.0242163196      -4.2498272043   -10.6726404285
  41 [43]     7.4820065980    12.0231995596      -1.1477307972     3.9728027190
  42 [10]    -3.6545658269   -18.0551079573      -0.5216089174   -15.2315827674
  43 [21]    16.9375513781   -23.8501660484       7.4820065980    12.0231995596

Now because we output things such that the 0,1,2,3 index offsets are the *columns* of of 4 x 9 set of out-array-indices,
this means that the output permutation translates (in terms of of 4 radix-11 macro outputs above) translates into the following:

	00,12,24,36,04,16,28,40,08,20,32
	11,23,35,03,15,27,39,07,19,31,43
	22,34,02,14,26,38,06,18,30,42,10
	33,01,13,25,37,05,17,29,41,09,21
*/

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy44_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr;
	#ifdef USE_SSE2
		const int pfetch_dist = PFETCH_DIST;
	#endif
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40;
		int poff[RADIX>>2];
		int j,j1,l,ntmp;
	#ifndef USE_SSE2
		int j2,jt,jp;
	#endif
	#ifdef USE_SSE2
	  #ifndef USE_AVX
		const int *incr;
	  #endif
		const int *inc_arr;
	  #ifdef USE_AVX512	// Have no specialized HIACC carry macro in AVX-512 and ARMv8 SIMD, so these get an "goes to 11" in LOACC mode via an incr_hiacc[] array:
		const int incr_long[] = {1}, incr_med[] = {1}, incr_short[] = {1}, incr_hiacc[] = {1};
	  #elif defined(USE_AVX)
		const int incr_long[] = {1}, incr_med[] = {1}, incr_short[] = {1}, incr_hiacc[] = {0};
	  #elif defined(USE_ARM_V8_SIMD)
		const int incr_long[] = {6,5}, incr_med[] = {4,3,4}, incr_short[] = {2,2,3,2,2}, incr_hiacc[] = {1,1,1,1,1,1,1,1,1,1,1};
	  #else
		const int incr_long[] = {6,5}, incr_med[] = {4,3,4}, incr_short[] = {2,2,3,2,2}, incr_hiacc[] = {0};
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
	#endif

	#ifndef USE_AVX
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
	#endif
	#ifdef USE_AVX512
		double t0,t1,t2,t3;
		struct uint32x8 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#elif defined(USE_AVX)
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#else
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	#endif
	#ifndef USE_SSE2
		double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	#endif
		uint32 p0123[4];

	#ifdef USE_SSE2
		const uint32 dft11_offs[11] = {0,0x2<<L2_SZ_VD,0x4<<L2_SZ_VD,0x6<<L2_SZ_VD,0x8<<L2_SZ_VD,0xa<<L2_SZ_VD,0xc<<L2_SZ_VD,0xe<<L2_SZ_VD,0x10<<L2_SZ_VD,0x12<<L2_SZ_VD,0x14<<L2_SZ_VD}, *dft11_offptr = &(dft11_offs[0]), *ui32_ptr;
		int i0,i1,i2,i3;
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
		vec_dbl /* *two,*one,*five, */ *ua0,/* *ua1,*ua2,*ua3,*ua4,*ua5,*ua6,*ua7,*ua8,*ua9, *ub0,*ub1,*ub2,*ub3,*ub4,*ub5,*ub6,*ub7,*ub8,*ub9, */ *max_err,
		  #ifndef USE_AVX512
			*sse2_rnd,
		  #endif
			*half_arr,
			*r00,*s1p00,
			*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx
		vec_dbl *tmp,*tm1;	// Non-static utility ptrs
	  #ifndef USE_AVX512
		vec_dbl *tm2;	// Non-static utility ptrs
		double dtmp;
	  #endif
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else
	  // FMA-based SIMD or (scalar-double) + (LO_ADD = 1 in masterdefs.h) use these sincos constants:
	  #if (defined(USE_AVX2) && DFT_11_FMA) || defined(USE_ARM_V8_SIMD) || (!defined(USE_SSE2) && (LO_ADD != 0))
		// FMA based on simple radix-11 DFT implementation, same as LO_ADD - more accurate, and with FMA, faster as well
		const double cc1 =  0.84125353283118116886,	/* Real part of exp(i*2*pi/11), the radix-11 fundamental sincos datum	*/
			ss1 =  0.54064081745559758210,	/* Imag part of exp(i*2*pi/11).	*/
			cc2 =  0.41541501300188642553,	/* cos(2u)	*/
			ss2 =  0.90963199535451837140,	/* sin(2u)	*/
			cc3 = -0.14231483827328514043,	/* cos(3u)	*/
			ss3 =  0.98982144188093273237,	/* sin(3u)	*/
			cc4 = -0.65486073394528506404,	/* cos(4u)	*/
			ss4 =  0.75574957435425828378,	/* sin(4u)	*/
			cc5 = -0.95949297361449738988,	/* cos(5u)	*/
			ss5 =  0.28173255684142969773;	/* sin(5u)	*/
	  #else	// AVX, SSE2 and non-LO_ADD (cf. masterdefs.h) scalar-double builds all use these:
		// 5 Dec 2014: Change const -> static to enable wrong-way-rounding sensitivity experiments:
		const double a0 = 2.31329240211767848235, /* a0 = (   cq0      -  cq3+  cq2-  cq4)		*/
			a1 =  1.88745388228838373902, /* a1 = (         cq1-  cq3+  cq2-  cq4)		*/
			a2 = -1.41435370755978245393, /* a2 = (-2*cq0-2*cq1+3*cq3-2*cq2+3*cq4)/5	*/
			a3 =  0.08670737584270518028, /* a3 = (-  cq0+  cq1-  cq3+  cq2      )		*/
			a4 = -0.73047075949850706917, /* a4 = (-  cq0+  cq1-  cq3      +  cq4)		*/
			a5 =  0.38639279888589610480, /* a5 = ( 3*cq0-2*cq1+3*cq3-2*cq2-2*cq4)/5	*/
			a6 =  0.51254589567199992361, /* a6 = (            -  cq3+  cq2      )		*/
			a7 =  1.07027574694717148957, /* a7 = (         cq1-  cq3            )		*/
			a8 = -0.55486073394528506404, /* a8 = (-  cq0-  cq1+4*cq3-  cq2-  cq4)/5	*/
			a9 = -1.10000000000000000000, /* a9 = (   cq0+  cq1+  cq3+  cq2+  cq4)/5 - 1*/
			/* Sine terms same form as cos, but cq1 => -sq1 and no -1 on b9 term: */
			b0 =  0.49298012814084233296, /* b0 = (   sq0      -  sq3+  sq2-  sq4)		*/
			b1 = -0.95729268466927362054, /* b1 = (      -  sq1-  sq3+  sq2-  sq4)		*/
			b2 =  0.37415717312460801167, /* b2 = (-2*sq0+2*sq1+3*sq3-2*sq2+3*sq4)/5	*/
			b3 = -1.21620094528344150491, /* b3 = (-  sq0-  sq1-  sq3+  sq2      )		*/
			b4 = -1.92428983032294453955, /* b4 = (-  sq0-  sq1-  sq3      +  sq4)		*/
			b5 =  0.63306543373877589604, /* b5 = ( 3*sq0+2*sq1+3*sq3-2*sq2-2*sq4)/5	*/
			b6 =  0.23407186752667444859, /* b6 = (            -  sq3+  sq2      )		*/
			b7 = -1.66538156970877665518, /* b7 = (      -  sq1-  sq3            )		*/
			b8 =  0.42408709531871829886, /* b8 = (-  sq0+  sq1+4*sq3-  sq2-  sq4)/5	*/
			b9 =  0.33166247903553998491; /* b9 = (   sq0-  sq1+  sq3+  sq2+  sq4)/5	*/
	  #endif
		double *base, *baseinv;
		const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int m,m2;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		int bjmodn[RADIX];	// Thread only carries a base datum here, must alloc a local array for remaining values
		double *cy = thread_arg->cy, rt,it,temp,frac;
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
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;
		p28 = p24 + p04;
		p32 = p28 + p04;
		p36 = p32 + p04;
		p40 = p36 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p36 = p36 + ( (p36 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );

		p0123[0] = 0; p0123[1] = p01; p0123[2] = p02; p0123[3] = p03;

		poff[0x0] =   0; poff[0x1] = p04; poff[0x2] = p08; poff[0x3] = p12;
		poff[0x4] = p16; poff[0x5] = p20; poff[0x6] = p24; poff[0x7] = p28;
		poff[0x8] = p32; poff[0x9] = p36; poff[0xa] = p40;

	#ifdef USE_SSE2
		r00 = thread_arg->r00;	tmp = r00 + 0x58;
		r00 = r00 + 0x00;		s1p00 = tmp + 0x00;
		tmp += 0x58;	// sc_ptr += 0xb0
		//two     = tmp + 0;	// AVX+ versions of various DFT macros need consts [2,1] pair laid out thusly
		//one     = tmp + 1;
	//	two     = tmp + 2; Unnamed ptr, used in AVX2 mode to hold 2.0 (and *five holds 1.0) for the radix-11 DFT code
		//five    = tmp + 3;
		ua0     = tmp + 4;
		//ua1     = tmp + 5;
		//ua2     = tmp + 6;
		//ua3     = tmp + 7;
		//ua4     = tmp + 8;
		//ua5     = tmp + 9;
		//ua6     = tmp + 0xa;
		//ua7     = tmp + 0xb;
		//ua8     = tmp + 0xc;
		//ua9     = tmp + 0xd;
		//ub0     = tmp + 0xe;
		//ub1     = tmp + 0xf;
		//ub2     = tmp + 0x10;
		//ub3     = tmp + 0x11;
		//ub4     = tmp + 0x12;
		//ub5     = tmp + 0x13;
		//ub6     = tmp + 0x14;
		//ub7     = tmp + 0x15;
		//ub8     = tmp + 0x16;
		//ub9     = tmp + 0x17;
		tmp += 0x18;	// sc_ptr += 0xde
	#ifdef USE_AVX512
		cy = tmp;		tmp += 6;		// RADIX/8 and round up
		max_err = tmp + 0x00;
		//sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	#elif defined(USE_AVX)
		cy = tmp;		tmp += 0x0b;	// sc_ptr += 0xe9
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0xeb; This is where the value of half_arr_offset44 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy = tmp;		tmp += 0x16;	// sc_ptr += 0xf4
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0xf6; This is where the value of half_arr_offset44 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*odd_radix] x 16 for Fermat-mod */
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

		sign_mask = (uint64*)(r00 + radix44_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (  #doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		sse_sw  = sse_bw    + RE_IM_STRIDE;
		sse_n   = sse_sw    + RE_IM_STRIDE;
	  #ifdef USE_AVX512
		n_minus_sil   = (struct uint32x8 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x8 *)sse_n + 2;
		sinwt         = (struct uint32x8 *)sse_n + 3;
		sinwtm1       = (struct uint32x8 *)sse_n + 4;
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
		#include "radix44_main_carry_loop.h"

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
