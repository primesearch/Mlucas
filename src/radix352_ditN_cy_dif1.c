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

#define RADIX		352	// Use #define rather than const int to ensure it's really a compile-time const in the C sense
#define ODD_RADIX	11	// ODD_RADIX = [radix >> trailz(radix)]

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
  // RADIX = 352 = 0x160:
  // For Mersenne-mod need (16 [SSE2] or 64 [AVX]) + (4 [HIACC] or 40 [LOACC]) added slots for half_arr lookup tables.
  // Max = (40 [SSE2]; 132 [AVX]), add to (half_arr_offset352 + RADIX) to get required value of radix352_creals_in_local_store:
  #ifdef USE_AVX512	// RADIX/8 = 44 = 0x2c fewer carry slots than AVX:
	const int half_arr_offset352 = 0x5ce;	// + RADIX = 0x72e; Used for thread local-storage-integrity checking
	const int radix352_creals_in_local_store = 0x7b4;	// (half_arr_offset + RADIX) + 132 (=0x84) and round up to nearest multiple of 4
  #elif defined(USE_AVX)
	const int half_arr_offset352 = 0x5fa;	// + RADIX = 0x75a; Used for thread local-storage-integrity checking
	const int radix352_creals_in_local_store = 0x7e0;	// (half_arr_offset + RADIX) + 132 (=0x84) and round up to nearest multiple of 4
  #else
	const int half_arr_offset352 = 0x652;	// + RADIX = 0x7b2; Used for thread local-storage-integrity checking
	const int radix352_creals_in_local_store = 0x7dc;	// (half_arr_offset + RADIX) + 40 (=0x28) and round up to nearest multiple of 4
  #endif

	#include "sse2_macro_gcc64.h"
	#include "radix11_sse_macro.h"
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

/****************/

int radix352_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-352 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-352 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix352_ditN_cy_dif1";
#if !defined(MULTITHREAD) && defined(USE_SSE2)
	const int pfetch_dist = PFETCH_DIST;
#endif
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
#if defined(USE_SSE2) && !((defined(USE_AVX2) && DFT_11_FMA) || defined(USE_ARM_V8_SIMD))
	struct qfloat qtheta,cq0,cq1,cq2,cq3,cq4,/* sq0, */sq1,sq2,sq3,sq4;
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
	int j2,jp,ntmp;
#endif
#ifndef MULTITHREAD
	int jstart, *ip;
#endif
#ifdef USE_SSE2
	int nbytes;
#endif
#if !defined(MULTITHREAD) && defined(USE_SSE2)
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 22|44|88 for avx512,avx,sse, respectively.
	// But fixed-incr too restrictive here, so 'divide 22|44|88 into pieces' via increment-array whose elts sum to 22|44|88:
	const int *incr,*inc_arr;
  #ifdef USE_AVX512	// Have no specialized HIACC carry macro in AVX-512 and ARMv8 SIMD, so these get an "goes to 11" in LOACC mode via an incr_hiacc[] array:
	const int incr_long[] = {11,11}, incr_med[] = {6,5,6,5}, incr_short[] = {4,3,4,4,3,4}, incr_hiacc[] = {2,2,2,2,2,2,2,2,2,2,2};
  #elif defined(USE_AVX)
	const int incr_long[] = {11,11,11,11}, incr_med[] = {6,5,6,5,6,5,6,5}, incr_short[] = {4,4,4,4,4,4,4,4,4,4,4}, incr_hiacc[] = {0};
  #elif defined(USE_ARM_V8_SIMD)
	const int incr_long[] = {11,11,11,11,11,11,11,11}, incr_med[] = {6,5,6,5,6,5,6,5,6,5,6,5,6,5,6,5}, incr_short[] = {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4}, incr_hiacc[] = {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
  #else
	const int incr_long[] = {11,11,11,11,11,11,11,11}, incr_med[] = {6,5,6,5,6,5,6,5,6,5,6,5,6,5,6,5}, incr_short[] = {4,3,4,4,3,4,4,3,4,4,3,4,4,3,4,4,3,4,4,3,4,4,3,4}, incr_hiacc[] = {0};
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
#endif // !MULTITHEAD && USE_SSE2

	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,p100,p110,p120,p130,p140,p150, nsave = 0;
	static int poff[RADIX>>2];	// Store mults of p4 offset for loop control
#ifndef MULTITHREAD
// Shared DIF+DIT:
  #ifndef USE_SSE2
	static int t_offsets[32];
  #endif
	// Need storage for 2 circular-shifts perms of a basic 11-vector, with shift count in [0,10] that means 2*21 elts:
	static int dif_offsets[RADIX], dif_p20_cperms[42], dif_p20_lo_offset[32], dif_phi[ODD_RADIX];
	static int dit_offsets[RADIX], dit_p20_cperms[42], dit_p20_lo_offset[32], dit_phi[ODD_RADIX];
  #ifdef USE_SSE2
	const uint32 dft11_offs[11] = {0,0x40<<L2_SZ_VD,0x80<<L2_SZ_VD,0xc0<<L2_SZ_VD,0x100<<L2_SZ_VD,0x140<<L2_SZ_VD,0x180<<L2_SZ_VD,0x1c0<<L2_SZ_VD,0x200<<L2_SZ_VD,0x240<<L2_SZ_VD,0x280<<L2_SZ_VD}, *dft11_offptr = &(dft11_offs[0]);
  #endif
#endif
	static double radix_inv, n2inv;
// FMA-based SIMD or (scalar-double) + (LO_ADD = 1 in masterdefs.h) use these sincos constants:
#if (defined(USE_AVX2) && DFT_11_FMA) || defined(USE_ARM_V8_SIMD) || (!defined(USE_SSE2) && (LO_ADD != 0))
	#warning Using FMA-heavy lo-add 11-DFT
  #if !defined(MULTITHREAD) || (defined(USE_AVX2) && DFT_11_FMA) || defined(USE_ARM_V8_SIMD)
	const double
		a1 =  0.84125353283118116886,	/* Real part of exp(i*2*pi/11), the radix-11 fundamental sincos datum	*/
		b1 =  0.54064081745559758210,	/* Imag part of exp(i*2*pi/11).	*/
		a2 =  0.41541501300188642553,	/* cos(2u)	*/
		b2 =  0.90963199535451837140,	/* sin(2u)	*/
		a3 = -0.14231483827328514043,	/* cos(3u)	*/
		b3 =  0.98982144188093273237,	/* sin(3u)	*/
		a4 = -0.65486073394528506404,	/* cos(4u)	*/
		b4 =  0.75574957435425828378,	/* sin(4u)	*/
		a5 = -0.95949297361449738988,	/* cos(5u)	*/
		b5 =  0.28173255684142969773;	/* sin(5u)	*/
  #endif
#else
	#warning LO_ADD = 0 defined at compile time ... Using FMA-lite / hi-add 11-DFT
	const double
		a0 =  2.31329240211767848235, /* a0 = (   cq0      -  cq3+  cq2-  cq4)		*/
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
#endif	// LO_ADD ?
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
	vec_dbl *tmp,*tm2;	// Non-static utility ptrs
  #ifndef MULTITHREAD
	vec_dbl *tm1;	// Non-static utility ptrs
  #endif
	static vec_dbl *two,*one,*sqrt2,*isrt2, *cc1,*ss1,*cc2,*ss2,*cc3,*ss3, *five,
		*ua0,*ua1,*ua2,*ua3,*ua4,*ua5,*ua6,*ua7,*ua8,*ua9,
		*ub0,*ub1,*ub2,*ub3,*ub4,*ub5,*ub6,*ub7,*ub8,*ub9,
		*max_err, *sse2_rnd, *half_arr,
		*r00,	// Head of RADIX*vec_cmplx-sized local store #1
	  #ifndef MULTITHREAD
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
	static task_control_t   task_control = {NULL, (void*)cy352_process_chunk, NULL, 0x0};

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
		// consisting of radix352_creals_in_local_store dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix352_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix352_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = sc_ptr;	r00   = tmp;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x2c0;					// Head of RADIX*vec_cmplx-sized local store #2
	  #ifndef MULTITHREAD
						s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
	  #endif
		tmp	+= 0x2c0;
		two     = tmp;	// AVX+ versions of various DFT macros assume consts 2.0,1.0,isrt2 laid out thusly
		one     = tmp + 0x1;
		sqrt2	= tmp + 0x2;
		isrt2   = tmp + 0x3;
		cc2	    = tmp + 0x4;	// Radix-32 DFT macros assume roots stored in this [8th, 16th, 32nd_1,3] order
		ss2	    = tmp + 0x5;
		cc1	    = tmp + 0x6;
		ss1	    = tmp + 0x7;
		cc3	    = tmp + 0x8;
		ss3	    = tmp + 0x9;
	//	two     = tmp + 0xa; Unnamed ptr, used in AVX2 mode to hold 2.0 (and *five holds 1.0) for the radix-11 DFT code
		five    = tmp + 0xb;
		ua0     = tmp + 0xc;
		ua1     = tmp + 0xd;
		ua2     = tmp + 0xe;
		ua3     = tmp + 0xf;
		ua4     = tmp + 0x10;
		ua5     = tmp + 0x11;
		ua6     = tmp + 0x12;
		ua7     = tmp + 0x13;
		ua8     = tmp + 0x14;
		ua9     = tmp + 0x15;
		ub0     = tmp + 0x16;
		ub1     = tmp + 0x17;
		ub2     = tmp + 0x18;
		ub3     = tmp + 0x19;
		ub4     = tmp + 0x1a;
		ub5     = tmp + 0x1b;
		ub6     = tmp + 0x1c;
		ub7     = tmp + 0x1d;
		ub8     = tmp + 0x1e;
		ub9     = tmp + 0x1f;
		tmp += 0x20;	// sc_ptr += 0x5a0
	  #ifdef USE_AVX512
		cy = tmp;		tmp += 0x2c;	// RADIX/8 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(5a0 + 2c + 2) = 0x5ce; This is where the value of half_arr_offset352 comes from
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
		cy = tmp;		tmp += 0x58;	// RADIX/4 vec_dbl slots
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(5a0 + 58 + 2) = 0x5fa; This is where the value of half_arr_offset352 comes from
		half_arr= tmp + 0x02;	// This table needs 20 vec_dbl in both avx and sse2 mode
	  #else
		cy = tmp;		tmp += 0xb0;	// RADIX/2 vec_dbl slots
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(5a0 + b0 + 2) = 0x652; This is where the value of half_arr_offset352 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*odd_radix] x 16 for Fermat-mod */
	  #endif

	  #if (defined(USE_AVX2) && DFT_11_FMA) || defined(USE_ARM_V8_SIMD)
	  	/* no-op */
	  #else
		/*======= Cf. the radix-44 routine for details on the wrong-way-rounding of a
		selected subset of DFT-11 consts in order to improve overall accuracy: =========*/
		qtheta = qfdiv(Q2PI, i64_to_q((uint64)11));
		qt =          qtheta ;	cq0 = qfcos(qt);	//sq0 = qfsin(qt);
		qt = qfadd(qt,qtheta);	cq1 = qfcos(qt);	sq1 = qfsin(qt);
		qt = qfadd(qt,qtheta);	cq2 = qfcos(qt);	sq2 = qfsin(qt);
		qt = qfadd(qt,qtheta);	cq3 = qfcos(qt);	sq3 = qfsin(qt);
		qt = qfadd(qt,qtheta);	cq4 = qfcos(qt);	sq4 = qfsin(qt);
		//================================================================
	  #endif
		ASSERT((radix352_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r00) + (20 << L2_SZ_VD), "radix352_creals_in_local_store checksum failed!");
		/* These remain fixed: */
		VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );
	  #if 1
		// 2 unnamed slots for alternate "rounded the other way" copies of sqrt2,isrt2:
		dtmp = u64_to_f64(sqrt2_dn);	VEC_DBL_INIT(sqrt2, dtmp);
		dtmp = u64_to_f64(isrt2_dn);	VEC_DBL_INIT(isrt2, dtmp);
	  #else
		VEC_DBL_INIT(sqrt2, SQRT2);		VEC_DBL_INIT(isrt2, ISRT2);
	  #endif
		VEC_DBL_INIT(cc2, c16  );	VEC_DBL_INIT(ss2, s16  );	// radix-32 DFT trig consts
		VEC_DBL_INIT(cc1, c32_1);	VEC_DBL_INIT(ss1, s32_1);
		VEC_DBL_INIT(cc3, c32_3);	VEC_DBL_INIT(ss3, s32_3);
	  #if (defined(USE_AVX2) && DFT_11_FMA) || defined(USE_ARM_V8_SIMD)	// FMA version based on simple radix-11 DFT implementation, same as LO_ADD
		tmp = five-1;
		VEC_DBL_INIT(tmp ,2.0);		VEC_DBL_INIT(five,1.0);	// *five-1,*five used for [2.0,1.0] here
		VEC_DBL_INIT(ua0 , a1);		VEC_DBL_INIT(ub0 ,0.0);	// upper 10 slots unused here; init = 0
		VEC_DBL_INIT(ua1 , a2);		VEC_DBL_INIT(ub1 ,0.0);
		VEC_DBL_INIT(ua2 , a3);		VEC_DBL_INIT(ub2 ,0.0);
		VEC_DBL_INIT(ua3 , a4);		VEC_DBL_INIT(ub3 ,0.0);
		VEC_DBL_INIT(ua4 , a5);		VEC_DBL_INIT(ub4 ,0.0);
		VEC_DBL_INIT(ua5 , b1);		VEC_DBL_INIT(ub5 ,0.0);
		VEC_DBL_INIT(ua6 , b2);		VEC_DBL_INIT(ub6 ,0.0);
		VEC_DBL_INIT(ua7 , b3);		VEC_DBL_INIT(ub7 ,0.0);
		VEC_DBL_INIT(ua8 , b4);		VEC_DBL_INIT(ub8 ,0.0);
		VEC_DBL_INIT(ua9 , b5);		VEC_DBL_INIT(ub9 ,0.0);
	  #else
		tmp = five-1;
		VEC_DBL_INIT(tmp ,2.0);		VEC_DBL_INIT(five,5.0);	// *five-1,*five used for [2.0,5.0] here
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
		p1 = NDIVR;		p10 = NDIVR<<4;		p100 = NDIVR<<8;
		p2 = p1 + p1;	p20 = p10 + p10;	p110 = p100 + p10;
		p3 = p2 + p1;	p30 = p20 + p10;	p120 = p110 + p10;
		p4 = p3 + p1;	p40 = p30 + p10;	p130 = p120 + p10;
		p5 = p4 + p1;	p50 = p40 + p10;	p140 = p130 + p10;
		p6 = p5 + p1;	p60 = p50 + p10;	p150 = p140 + p10;
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
													p120 += ( (p120 >> DAT_BITS) << PAD_BITS );
													p130 += ( (p130 >> DAT_BITS) << PAD_BITS );
													p140 += ( (p140 >> DAT_BITS) << PAD_BITS );
													p150 += ( (p150 >> DAT_BITS) << PAD_BITS );
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
		poff[0x40+0] =p100; poff[0x40+1] =p100+p4; poff[0x40+2] =p100+p8; poff[0x40+3] =p100+pc;
		poff[0x44+0] =p110; poff[0x44+1] =p110+p4; poff[0x44+2] =p110+p8; poff[0x44+3] =p110+pc;
		poff[0x48+0] =p120; poff[0x48+1] =p120+p4; poff[0x48+2] =p120+p8; poff[0x48+3] =p120+pc;
		poff[0x4c+0] =p130; poff[0x4c+1] =p130+p4; poff[0x4c+2] =p130+p8; poff[0x4c+3] =p130+pc;
		poff[0x50+0] =p140; poff[0x50+1] =p140+p4; poff[0x50+2] =p140+p8; poff[0x50+3] =p140+pc;
		poff[0x54+0] =p150; poff[0x54+1] =p150+p4; poff[0x54+2] =p150+p8; poff[0x54+3] =p150+pc;

	#ifndef MULTITHREAD

		l = 0;
		dif_phi[l++] =   0;
		dif_phi[l++] = p140;
		dif_phi[l++] = p120;
		dif_phi[l++] = p100;
		dif_phi[l++] = pe0;
		dif_phi[l++] = pc0;
		dif_phi[l++] = pa0;
		dif_phi[l++] = p80;
		dif_phi[l++] = p60;
		dif_phi[l++] = p40;
		dif_phi[l++] = p20;
		l = 0;
		dit_phi[l++] =   0;
		dit_phi[l++] = p20;
		dit_phi[l++] = p40;
		dit_phi[l++] = p60;
		dit_phi[l++] = p80;
		dit_phi[l++] = pa0;
		dit_phi[l++] = pc0;
		dit_phi[l++] = pe0;
		dit_phi[l++] =p100;
		dit_phi[l++] =p120;
		dit_phi[l++] =p140;

	// Shared:
	  #ifndef USE_SSE2
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
	  #endif

	/*** DIF indexing stuff: ***/

		// Init storage for 2 circular-shifts perms of a basic 11-vector, with shift count in [0,10] that means 2*22 = 42 elts:
	   #ifdef USE_SSE2
		/* Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		replacing p10*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode.
		But since the resulting offset-indices will be getting added to a vec_dbl* which implies a pointer-arithmetic
		addend left-shifted [L2_SZ_VD] bits, which will be 4,5,6 for sse2,avx,avx-512 and since x86 only supports
		an LEA instruction with a max-shift of 3, to avoid the resulting performance penalty, simply inline the needed
		operand-shift in this init stage, thus allowing the ensuing pointer arithmetic to be done via simple ADD:
		*/
		i = 0; j = L2_SZ_VD+1;
		// Even multiples of p10 cshift array:
		dif_p20_cperms[i++] = 0;
		dif_p20_cperms[i++] = 0x140<<j;
		dif_p20_cperms[i++] = 0x120<<j;
		dif_p20_cperms[i++] = 0x100<<j;
		dif_p20_cperms[i++] = 0xe0<<j;
		dif_p20_cperms[i++] = 0xc0<<j;
		dif_p20_cperms[i++] = 0xa0<<j;
		dif_p20_cperms[i++] = 0x80<<j;
		dif_p20_cperms[i++] = 0x60<<j;
		dif_p20_cperms[i++] = 0x40<<j;
		dif_p20_cperms[i++] = 0x20<<j;
		while(i < 2*ODD_RADIX-1) {
			dif_p20_cperms[i] = dif_p20_cperms[i - ODD_RADIX]; ++i;
		}
		// Odd multiples of p10 cshift array:
		dif_p20_cperms[i++] = 0x150<<j;
		dif_p20_cperms[i++] = 0x130<<j;
		dif_p20_cperms[i++] = 0x110<<j;
		dif_p20_cperms[i++] = 0xf0<<j;
		dif_p20_cperms[i++] = 0xd0<<j;
		dif_p20_cperms[i++] = 0xb0<<j;
		dif_p20_cperms[i++] = 0x90<<j;
		dif_p20_cperms[i++] = 0x70<<j;
		dif_p20_cperms[i++] = 0x50<<j;
		dif_p20_cperms[i++] = 0x30<<j;
		dif_p20_cperms[i++] = 0x10<<j;	// *** This last shift had somehow gotten munged from <<j to <<1 in v19, but only in this 0thr-init section...since 0thr-build is deprecated, nobody noticed.
		while(i < 4*ODD_RADIX-2) {
			dif_p20_cperms[i] = dif_p20_cperms[i - ODD_RADIX]; ++i;
		}
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-11 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-10; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-11 index, bits <4:30> for the p0-f:
		i = 0;
		dif_p20_lo_offset[i++] = ((  0 << 5) + 0 );
		dif_p20_lo_offset[i++] = ((0x5 << 5) + 0 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((0xa << 5) + 1 );
		dif_p20_lo_offset[i++] = ((0xf << 5) + 1 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((0x4 << 5) + 1 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((0x9 << 5) + 2 );
		dif_p20_lo_offset[i++] = ((0xe << 5) + 2 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((0x3 << 5) + 2 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((0x8 << 5) + 3 );
		dif_p20_lo_offset[i++] = ((0xd << 5) + 3 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((0x2 << 5) + 3 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((0x7 << 5) + 4 );
		dif_p20_lo_offset[i++] = ((0xc << 5) + 4 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((0x1 << 5) + 4 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((0x6 << 5) + 5 );
		dif_p20_lo_offset[i++] = ((0xb << 5) + 5 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((  0 << 5) + 5 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((0x5 << 5) + 6 );
		dif_p20_lo_offset[i++] = ((0xa << 5) + 6 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((0xf << 5) + 7 );
		dif_p20_lo_offset[i++] = ((0x4 << 5) + 7 );
		dif_p20_lo_offset[i++] = ((0x9 << 5) + 7 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((0xe << 5) + 8 );
		dif_p20_lo_offset[i++] = ((0x3 << 5) + 8 );
		dif_p20_lo_offset[i++] = ((0x8 << 5) + 8 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((0xd << 5) + 9 );
		dif_p20_lo_offset[i++] = ((0x2 << 5) + 9 );
		dif_p20_lo_offset[i++] = ((0x7 << 5) + 9 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((0xc << 5) + 10);
		dif_p20_lo_offset[i++] = ((0x1 << 5) + 10);
		dif_p20_lo_offset[i++] = ((0x6 << 5) + 10)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((0xb << 5) + 0 );

	   #else

		// Init storage for 2 circular-shifts perms of a basic 11-vector, with shift count in [0,10] that means 2*22 = 42 elts:
		i = 0;
		// Even multiples of p10 cshift array:
		dif_p20_cperms[i++] = 0;
		dif_p20_cperms[i++] = p140;
		dif_p20_cperms[i++] = p120;
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
		dif_p20_cperms[i++] = p150;
		dif_p20_cperms[i++] = p130;
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
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-11 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-10; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-11 index, bits <4:30> for the p0-f:
		i = 0;
		dif_p20_lo_offset[i++] = (( 0 << 4) + 0 );
		dif_p20_lo_offset[i++] = ((p5 << 4) + 0 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pa << 4) + 1 );
		dif_p20_lo_offset[i++] = ((pf << 4) + 1 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p4 << 4) + 1 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p9 << 4) + 2 );
		dif_p20_lo_offset[i++] = ((pe << 4) + 2 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p3 << 4) + 2 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p8 << 4) + 3 );
		dif_p20_lo_offset[i++] = ((pd << 4) + 3 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p2 << 4) + 3 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p7 << 4) + 4 );
		dif_p20_lo_offset[i++] = ((pc << 4) + 4 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p1 << 4) + 4 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p6 << 4) + 5 );
		dif_p20_lo_offset[i++] = ((pb << 4) + 5 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = (( 0 << 4) + 5 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p5 << 4) + 6 );
		dif_p20_lo_offset[i++] = ((pa << 4) + 6 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pf << 4) + 7 );
		dif_p20_lo_offset[i++] = ((p4 << 4) + 7 );
		dif_p20_lo_offset[i++] = ((p9 << 4) + 7 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pe << 4) + 8 );
		dif_p20_lo_offset[i++] = ((p3 << 4) + 8 );
		dif_p20_lo_offset[i++] = ((p8 << 4) + 8 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pd << 4) + 9 );
		dif_p20_lo_offset[i++] = ((p2 << 4) + 9 );
		dif_p20_lo_offset[i++] = ((p7 << 4) + 9 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pc << 4) + 10);
		dif_p20_lo_offset[i++] = ((p1 << 4) + 10);
		dif_p20_lo_offset[i++] = ((p6 << 4) + 10)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pb << 4) + 0 );

	   #endif	// sse2?

	// dif_offsets are w.r.to a-array, need 11 distinct sets of these, one for each DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		i = 0; ip = dif_offsets;	// In the following, the 'big' += p[>= 20] offsets already in dif_phi[], only care about += p0|p10:
		// Set 0: [0,1,2,3,5,4,7,6,b,a,8,9,e,f,d,c + p00],[7,6,4,5,1,0,3,2,d,c,f,e,a,b,9,8 + p10]
		*ip++= 0;*ip++=p1;*ip++=p2;*ip++=p3;*ip++=p5;*ip++=p4;*ip++=p7;*ip++=p6;*ip++=pb;*ip++=pa;*ip++=p8;*ip++=p9;*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16; *ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16;
		// Set 1: [8,9,a,b,d,c,f,e,7,6,4,5,1,0,3,2 +p140],[f,e,c,d,9,8,b,a,3,2,0,1,6,7,5,4 +p150]
		*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16; *ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd;*ip++=p9;*ip++=p8;*ip++=pb;*ip++=pa;*ip++=p3;*ip++=p2;*ip++= 0;*ip++=p1;*ip++=p6;*ip++=p7;*ip++=p5;*ip++=p4; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16;
		// Set 2: [4,5,6,7,3,2,0,1,f,e,c,d,9,8,b,a +p130],[8,9,a,b,d,c,f,e,7,6,4,5,1,0,3,2 +p120]
		*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=p3;*ip++=p2;*ip++= 0;*ip++=p1;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd;*ip++=p9;*ip++=p8;*ip++=pb;*ip++=pa; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16; *ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16;
		// Set 3: [2,3,1,0,7,6,4,5,8,9,a,b,d,c,f,e +p100],[4,5,6,7,3,2,0,1,f,e,c,d,9,8,b,a +p110]
		*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16; *ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=p3;*ip++=p2;*ip++= 0;*ip++=p1;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd;*ip++=p9;*ip++=p8;*ip++=pb;*ip++=pa; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16;
		// Set 4: [a,b,9,8,f,e,c,d,4,5,6,7,3,2,0,1 + pf0],[2,3,1,0,7,6,4,5,8,9,a,b,d,c,f,e + pe0]
		*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=p3;*ip++=p2;*ip++= 0;*ip++=p1; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16; *ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16;
		// Set 5: [e,f,d,c,8,9,a,b,2,3,1,0,7,6,4,5 + pc0],[a,b,9,8,f,e,c,d,4,5,6,7,3,2,0,1 + pd0]
		*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16; *ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=p3;*ip++=p2;*ip++= 0;*ip++=p1; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16;
		// Set 6: [1,0,3,2,4,5,6,7,a,b,9,8,f,e,c,d + pb0],[e,f,d,c,8,9,a,b,2,3,1,0,7,6,4,5 + pa0]
		*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16; *ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16;
		// Set 7: [5,4,7,6,2,3,1,0,e,f,d,c,8,9,a,b + p80],[1,0,3,2,4,5,6,7,a,b,9,8,f,e,c,d + p90]
		*ip++=p5;*ip++=p4;*ip++=p7;*ip++=p6;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16; *ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16;
		// Set 8: [d,c,f,e,a,b,9,8,1,0,3,2,4,5,6,7 + p70],[5,4,7,6,2,3,1,0,e,f,d,c,8,9,a,b + p60]
		*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16; *ip++=p5;*ip++=p4;*ip++=p7;*ip++=p6;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16;
		// Set 9: [b,a,8,9,e,f,d,c,5,4,7,6,2,3,1,0 + p40],[d,c,f,e,a,b,9,8,1,0,3,2,4,5,6,7 + p50]
		*ip++=pb;*ip++=pa;*ip++=p8;*ip++=p9;*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p5;*ip++=p4;*ip++=p7;*ip++=p6;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16; *ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16;
		// Set A: [7,6,4,5,1,0,3,2,d,c,f,e,a,b,9,8 + p30],[b,a,8,9,e,f,d,c,5,4,7,6,2,3,1,0 + p20]
		*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16; *ip++=pb;*ip++=pa;*ip++=p8;*ip++=p9;*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p5;*ip++=p4;*ip++=p7;*ip++=p6;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16;
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
		// Init storage for 2 circular-shifts perms of a basic 11-vector, with shift count in [0,10] that means 2*22 = 42 elts:
		i = 0; j = L2_SZ_VD+1;
		// Even multiples of p10 cshift array:
		dit_p20_cperms[i++] = 0;
		dit_p20_cperms[i++] = 0x140<<j;
		dit_p20_cperms[i++] = 0x120<<j;
		dit_p20_cperms[i++] = 0x100<<j;
		dit_p20_cperms[i++] = 0xe0<<j;
		dit_p20_cperms[i++] = 0xc0<<j;
		dit_p20_cperms[i++] = 0xa0<<j;
		dit_p20_cperms[i++] = 0x80<<j;
		dit_p20_cperms[i++] = 0x60<<j;
		dit_p20_cperms[i++] = 0x40<<j;
		dit_p20_cperms[i++] = 0x20<<j;
		while(i < 2*ODD_RADIX-1) {
			dit_p20_cperms[i] = dit_p20_cperms[i - ODD_RADIX]; ++i;
		}
		// Odd multiples of p10 cshift array:
		dit_p20_cperms[i++] = 0x150<<j;
		dit_p20_cperms[i++] = 0x130<<j;
		dit_p20_cperms[i++] = 0x110<<j;
		dit_p20_cperms[i++] = 0xf0<<j;
		dit_p20_cperms[i++] = 0xd0<<j;
		dit_p20_cperms[i++] = 0xb0<<j;
		dit_p20_cperms[i++] = 0x90<<j;
		dit_p20_cperms[i++] = 0x70<<j;
		dit_p20_cperms[i++] = 0x50<<j;
		dit_p20_cperms[i++] = 0x30<<j;
		dit_p20_cperms[i++] = 0x10<<j;
		while(i < 4*ODD_RADIX-2) {
			dit_p20_cperms[i] = dit_p20_cperms[i - ODD_RADIX]; ++i;
		}

		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-11 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-8 index, bits <4:30> for the p0-f:
		i = 0;
		dit_p20_lo_offset[i++] = ((  0 << 5) + 0 );
		dit_p20_lo_offset[i++] = ((0xf << 5) + 1 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((0xe << 5) + 2 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((0xd << 5) + 3 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((0xc << 5) + 4 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((0xb << 5) + 5 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((0xa << 5) + 6 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((0x9 << 5) + 7 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((0x8 << 5) + 8 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((0x7 << 5) + 9 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((0x6 << 5) + 10) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((0x5 << 5) + 0 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((0x4 << 5) + 1 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((0x3 << 5) + 2 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((0x2 << 5) + 3 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((0x1 << 5) + 4 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((  0 << 5) + 5 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((0xf << 5) + 7 );
		dit_p20_lo_offset[i++] = ((0xe << 5) + 8 );
		dit_p20_lo_offset[i++] = ((0xd << 5) + 9 );
		dit_p20_lo_offset[i++] = ((0xc << 5) + 10);
		dit_p20_lo_offset[i++] = ((0xb << 5) + 0 );
		dit_p20_lo_offset[i++] = ((0xa << 5) + 1 );
		dit_p20_lo_offset[i++] = ((0x9 << 5) + 2 );
		dit_p20_lo_offset[i++] = ((0x8 << 5) + 3 );
		dit_p20_lo_offset[i++] = ((0x7 << 5) + 4 );
		dit_p20_lo_offset[i++] = ((0x6 << 5) + 5 );
		dit_p20_lo_offset[i++] = ((0x5 << 5) + 6 );
		dit_p20_lo_offset[i++] = ((0x4 << 5) + 7 );
		dit_p20_lo_offset[i++] = ((0x3 << 5) + 8 );
		dit_p20_lo_offset[i++] = ((0x2 << 5) + 9 );
		dit_p20_lo_offset[i++] = ((0x1 << 5) + 10);

	   #else

		// Init storage for 2 circular-shifts perms of a basic 11-vector, with shift count in [0,10] that means 2*22 = 42 elts:
		i = 0;
		// Even multiples of p10 cshift array:
		dit_p20_cperms[i++] = 0;
		dit_p20_cperms[i++] = p140;
		dit_p20_cperms[i++] = p120;
		dit_p20_cperms[i++] = p100;
		dit_p20_cperms[i++] = pe0;
		dit_p20_cperms[i++] = pc0;
		dit_p20_cperms[i++] = pa0;
		dit_p20_cperms[i++] = p80;
		dit_p20_cperms[i++] = p60;
		dit_p20_cperms[i++] = p40;
		dit_p20_cperms[i++] = p20;
		while(i < 2*ODD_RADIX-1) {
			dit_p20_cperms[i] = dit_p20_cperms[i - ODD_RADIX]; ++i;
		}
		// Odd multiples of p10 cshift array:
		dit_p20_cperms[i++] = p150;
		dit_p20_cperms[i++] = p130;
		dit_p20_cperms[i++] = p110;
		dit_p20_cperms[i++] = pf0;
		dit_p20_cperms[i++] = pd0;
		dit_p20_cperms[i++] = pb0;
		dit_p20_cperms[i++] = p90;
		dit_p20_cperms[i++] = p70;
		dit_p20_cperms[i++] = p50;
		dit_p20_cperms[i++] = p30;
		dit_p20_cperms[i++] = p10;
		while(i < 4*ODD_RADIX-2) {
			dit_p20_cperms[i] = dit_p20_cperms[i - ODD_RADIX]; ++i;
		}

		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-11 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-8 index, bits <4:30> for the p0-f:
		i = 0;
		dit_p20_lo_offset[i++] = (( 0 << 4) + 0 );
		dit_p20_lo_offset[i++] = ((pf << 4) + 1 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pe << 4) + 2 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pd << 4) + 3 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pc << 4) + 4 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pb << 4) + 5 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pa << 4) + 6 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p9 << 4) + 7 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p8 << 4) + 8 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p7 << 4) + 9 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p6 << 4) + 10) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p5 << 4) + 0 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p4 << 4) + 1 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p3 << 4) + 2 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p2 << 4) + 3 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p1 << 4) + 4 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = (( 0 << 4) + 5 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pf << 4) + 7 );
		dit_p20_lo_offset[i++] = ((pe << 4) + 8 );
		dit_p20_lo_offset[i++] = ((pd << 4) + 9 );
		dit_p20_lo_offset[i++] = ((pc << 4) + 10);
		dit_p20_lo_offset[i++] = ((pb << 4) + 0 );
		dit_p20_lo_offset[i++] = ((pa << 4) + 1 );
		dit_p20_lo_offset[i++] = ((p9 << 4) + 2 );
		dit_p20_lo_offset[i++] = ((p8 << 4) + 3 );
		dit_p20_lo_offset[i++] = ((p7 << 4) + 4 );
		dit_p20_lo_offset[i++] = ((p6 << 4) + 5 );
		dit_p20_lo_offset[i++] = ((p5 << 4) + 6 );
		dit_p20_lo_offset[i++] = ((p4 << 4) + 7 );
		dit_p20_lo_offset[i++] = ((p3 << 4) + 8 );
		dit_p20_lo_offset[i++] = ((p2 << 4) + 9 );
		dit_p20_lo_offset[i++] = ((p1 << 4) + 10);

	   #endif	// sse2?

	// dit_offsets are w.r.to a-array, need 11 distinct sets of these, one for each radix-32 DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		i = 0; ip = dit_offsets;	// In the following, the 'big' += p[>= 20] offsets already in dif_phi[], only care about += p0|p10:
		// Set 0: [0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00],[f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10]
		*ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2;*ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4;*ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc;*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16; *ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc;*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8;*ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4;*ip++=p3;*ip++=p2;*ip++=p1;*ip++= 0; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16;
		// Set 1: [7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + p30],[7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + p20]
		*ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4;*ip++=p3;*ip++=p2;*ip++=p1;*ip++= 0;*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8;*ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16; *ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4;*ip++=p3;*ip++=p2;*ip++=p1;*ip++= 0;*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8;*ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16;
		// Set 2: [b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p40],[3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p50]
		*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8;*ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf;*ip++=p3;*ip++=p2;*ip++=p1;*ip++= 0;*ip++=p5;*ip++=p4;*ip++=p6;*ip++=p7; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16; *ip++=p3;*ip++=p2;*ip++=p1;*ip++= 0;*ip++=p5;*ip++=p4;*ip++=p6;*ip++=p7;*ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf;*ip++=p9;*ip++=p8;*ip++=pa;*ip++=pb; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16;
		// Set 3: [d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + p70],[d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + p60]
		*ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf;*ip++=p9;*ip++=p8;*ip++=pa;*ip++=pb;*ip++=p5;*ip++=p4;*ip++=p6;*ip++=p7;*ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16; *ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf;*ip++=p9;*ip++=p8;*ip++=pa;*ip++=pb;*ip++=p5;*ip++=p4;*ip++=p6;*ip++=p7;*ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16;
		// Set 4: [5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p80],[9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p90]
		*ip++=p5;*ip++=p4;*ip++=p6;*ip++=p7;*ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3;*ip++=p9;*ip++=p8;*ip++=pa;*ip++=pb;*ip++=pe;*ip++=pf;*ip++=pc;*ip++=pd; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16; *ip++=p9;*ip++=p8;*ip++=pa;*ip++=pb;*ip++=pe;*ip++=pf;*ip++=pc;*ip++=pd;*ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3;*ip++=p6;*ip++=p7;*ip++=p4;*ip++=p5; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16;
		// Set 5: [1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pb0],[1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pa0]
		*ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3;*ip++=p6;*ip++=p7;*ip++=p4;*ip++=p5;*ip++=pe;*ip++=pf;*ip++=pc;*ip++=pd;*ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16; *ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3;*ip++=p6;*ip++=p7;*ip++=p4;*ip++=p5;*ip++=pe;*ip++=pf;*ip++=pc;*ip++=pd;*ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16;
		// Set 6: [e,f,c,d,a,b,8,9,6,7,4,5,2,3,0,1 + pc0],[6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pd0]
		*ip++=pe;*ip++=pf;*ip++=pc;*ip++=pd;*ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9;*ip++=p6;*ip++=p7;*ip++=p4;*ip++=p5;*ip++=p2;*ip++=p3;*ip++= 0;*ip++=p1; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16; *ip++=p6;*ip++=p7;*ip++=p4;*ip++=p5;*ip++=p2;*ip++=p3;*ip++= 0;*ip++=p1;*ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9;*ip++=pc;*ip++=pd;*ip++=pf;*ip++=pe; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16;
		// Set 7: [a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + pf0],[a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + pe0]
		*ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9;*ip++=pc;*ip++=pd;*ip++=pf;*ip++=pe;*ip++=p2;*ip++=p3;*ip++= 0;*ip++=p1;*ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16; *ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9;*ip++=pc;*ip++=pd;*ip++=pf;*ip++=pe;*ip++=p2;*ip++=p3;*ip++= 0;*ip++=p1;*ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16;
		// Set 8: [2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a +p100],[c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 +p110]
		*ip++=p2;*ip++=p3;*ip++= 0;*ip++=p1;*ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6;*ip++=pc;*ip++=pd;*ip++=pf;*ip++=pe;*ip++=p8;*ip++=p9;*ip++=pb;*ip++=pa; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16; *ip++=pc;*ip++=pd;*ip++=pf;*ip++=pe;*ip++=p8;*ip++=p9;*ip++=pb;*ip++=pa;*ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6;*ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16;
		// Set 9: [4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c +p130],[4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c +p120]
		*ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6;*ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2;*ip++=p8;*ip++=p9;*ip++=pb;*ip++=pa;*ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16; *ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6;*ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2;*ip++=p8;*ip++=p9;*ip++=pb;*ip++=pa;*ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16;
		// Set A: [8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 +p140],[0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 +p150]
		*ip++=p8;*ip++=p9;*ip++=pb;*ip++=pa;*ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc;*ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2;*ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16; *ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2;*ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4;*ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc;*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16;
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


/*...The radix-352 final DIT pass is here.	*/

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
		#include "radix352_main_carry_loop.h"

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
		ASSERT(0x0 == cy352_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
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

void radix352_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-352 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix11_dif_pass for details on the radix-11 subtransform.
*/
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
	int i,j,j1,j2,jt,jp, *ip;
	// p-indexing is hexadecimal here:
	static int NDIVR,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,p100,p110,p120,p130,p140,p150, first_entry=TRUE;
	static int t_offsets[32], dif_offsets[RADIX];
	// Need storage for 2 circular-shifts perms of a basic 11-vector, with shift count in [0,10] that means 2*22 elts:
	static int dif_p20_cperms[42], dif_p20_lo_offset[32], dif_phi[ODD_RADIX];
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
		p3 = p2 + p1;	p30 = p20 + p10;	p120 = p110 + p10;
		p4 = p3 + p1;	p40 = p30 + p10;	p130 = p120 + p10;
		p5 = p4 + p1;	p50 = p40 + p10;	p140 = p130 + p10;
		p6 = p5 + p1;	p60 = p50 + p10;	p150 = p140 + p10;
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
													p120 += ( (p120 >> DAT_BITS) << PAD_BITS );
													p130 += ( (p130 >> DAT_BITS) << PAD_BITS );
													p140 += ( (p140 >> DAT_BITS) << PAD_BITS );
													p150 += ( (p150 >> DAT_BITS) << PAD_BITS );
	#if 0	// Monotone initial-o-indices used for extracting the needed operm:
		i = 0;
		dif_phi[i++] =   0;
		dif_phi[i++] = p20;
		dif_phi[i++] = p40;
		dif_phi[i++] = p60;
		dif_phi[i++] = p80;
		dif_phi[i++] = pa0;
		dif_phi[i++] = pc0;
		dif_phi[i++] = pe0;
		dif_phi[i++] = p100;
		dif_phi[i++] = p120;
		dif_phi[i++] = p140;
	#else
		i = 0;				// p-offsets for outputs have 'big' part of form p[0,16,10,4,14,8,2,12,6]0,
		dif_phi[i++] =   0;	// i.e. multiple-of-16 part gets decremented by 6 (mod 16) each term of the
		dif_phi[i++] = p140;// sequence, aside from the initial pair of terms, 0,16.
		dif_phi[i++] = p120;
		dif_phi[i++] = p100;
		dif_phi[i++] = pe0;
		dif_phi[i++] = pc0;
		dif_phi[i++] = pa0;
		dif_phi[i++] = p80;
		dif_phi[i++] = p60;
		dif_phi[i++] = p40;
		dif_phi[i++] = p20;
	#endif
		// Set array offsets for radix-32 DFT in/outputs:
		// t_offsets w.r.to: t-array, same for all 11 32-DFTs:
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

		// Init storage for 2 circular-shifts perms of a basic 11-vector, with shift count in [0,10] that means 2*22 = 42 elts:
		i = 0;
		// Even multiples of p10 cshift array:
		dif_p20_cperms[i++] = 0;
		dif_p20_cperms[i++] = p140;
		dif_p20_cperms[i++] = p120;
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
		dif_p20_cperms[i++] = p150;
		dif_p20_cperms[i++] = p130;
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
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-11 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-10; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-11 index, bits <4:30> for the p0-f:
		i = 0;
		dif_p20_lo_offset[i++] = (( 0 << 4) + 0 );
		dif_p20_lo_offset[i++] = ((p5 << 4) + 0 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pa << 4) + 1 );
		dif_p20_lo_offset[i++] = ((pf << 4) + 1 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p4 << 4) + 1 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p9 << 4) + 2 );
		dif_p20_lo_offset[i++] = ((pe << 4) + 2 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p3 << 4) + 2 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p8 << 4) + 3 );
		dif_p20_lo_offset[i++] = ((pd << 4) + 3 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p2 << 4) + 3 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p7 << 4) + 4 );
		dif_p20_lo_offset[i++] = ((pc << 4) + 4 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p1 << 4) + 4 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p6 << 4) + 5 );
		dif_p20_lo_offset[i++] = ((pb << 4) + 5 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = (( 0 << 4) + 5 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((p5 << 4) + 6 );
		dif_p20_lo_offset[i++] = ((pa << 4) + 6 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pf << 4) + 7 );
		dif_p20_lo_offset[i++] = ((p4 << 4) + 7 );
		dif_p20_lo_offset[i++] = ((p9 << 4) + 7 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pe << 4) + 8 );
		dif_p20_lo_offset[i++] = ((p3 << 4) + 8 );
		dif_p20_lo_offset[i++] = ((p8 << 4) + 8 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pd << 4) + 9 );
		dif_p20_lo_offset[i++] = ((p2 << 4) + 9 );
		dif_p20_lo_offset[i++] = ((p7 << 4) + 9 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pc << 4) + 10);
		dif_p20_lo_offset[i++] = ((p1 << 4) + 10);
		dif_p20_lo_offset[i++] = ((p6 << 4) + 10)+ ((uint32)1 << 31);
		dif_p20_lo_offset[i++] = ((pb << 4) + 0 );

	#if 0	// Monotone initial-o-indices used for extracting the needed operm:
		ip = dif_offsets;
		*ip++=0;*ip++=p1;*ip++=p2;*ip++=p3;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=pc;*ip++=pd;*ip++=pe;*ip++=pf;
		for(i = 0; i < 16; i++) { *ip++ = dif_offsets[i] + p10; }
		// And now make (ODD_RADIX-1) copies of the initial 32 elts:
		for(i = 0; i < 32; i++) { *ip++ = dif_offsets[i]; }
		for(i = 0; i < 32; i++) { *ip++ = dif_offsets[i]; }
		for(i = 0; i < 32; i++) { *ip++ = dif_offsets[i]; }
		for(i = 0; i < 32; i++) { *ip++ = dif_offsets[i]; }
		for(i = 0; i < 32; i++) { *ip++ = dif_offsets[i]; }
		for(i = 0; i < 32; i++) { *ip++ = dif_offsets[i]; }
		for(i = 0; i < 32; i++) { *ip++ = dif_offsets[i]; }
		for(i = 0; i < 32; i++) { *ip++ = dif_offsets[i]; }
		for(i = 0; i < 32; i++) { *ip++ = dif_offsets[i]; }
		for(i = 0; i < 32; i++) { *ip++ = dif_offsets[i]; }
	#else
	// dif_offsets are w.r.to a-array, need 11 distinct sets of these, one for each DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		i = 0; ip = dif_offsets;	// In the following, the 'big' += p[>= 20] offsets already in dif_phi[], only care about += p0|p10:
		// Set 0: [0,1,2,3,5,4,7,6,b,a,8,9,e,f,d,c + p00],[7,6,4,5,1,0,3,2,d,c,f,e,a,b,9,8 + p10]
		*ip++= 0;*ip++=p1;*ip++=p2;*ip++=p3;*ip++=p5;*ip++=p4;*ip++=p7;*ip++=p6;*ip++=pb;*ip++=pa;*ip++=p8;*ip++=p9;*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16; *ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16;
		// Set 1: [8,9,a,b,d,c,f,e,7,6,4,5,1,0,3,2 +p140],[f,e,c,d,9,8,b,a,3,2,0,1,6,7,5,4 +p150]
		*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16; *ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd;*ip++=p9;*ip++=p8;*ip++=pb;*ip++=pa;*ip++=p3;*ip++=p2;*ip++= 0;*ip++=p1;*ip++=p6;*ip++=p7;*ip++=p5;*ip++=p4; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16;
		// Set 2: [4,5,6,7,3,2,0,1,f,e,c,d,9,8,b,a +p130],[8,9,a,b,d,c,f,e,7,6,4,5,1,0,3,2 +p120]
		*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=p3;*ip++=p2;*ip++= 0;*ip++=p1;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd;*ip++=p9;*ip++=p8;*ip++=pb;*ip++=pa; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16; *ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16;
		// Set 3: [2,3,1,0,7,6,4,5,8,9,a,b,d,c,f,e +p100],[4,5,6,7,3,2,0,1,f,e,c,d,9,8,b,a +p110]
		*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16; *ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=p3;*ip++=p2;*ip++= 0;*ip++=p1;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd;*ip++=p9;*ip++=p8;*ip++=pb;*ip++=pa; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16;
		// Set 4: [a,b,9,8,f,e,c,d,4,5,6,7,3,2,0,1 + pf0],[2,3,1,0,7,6,4,5,8,9,a,b,d,c,f,e + pe0]
		*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=p3;*ip++=p2;*ip++= 0;*ip++=p1; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16; *ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16;
		// Set 5: [e,f,d,c,8,9,a,b,2,3,1,0,7,6,4,5 + pc0],[a,b,9,8,f,e,c,d,4,5,6,7,3,2,0,1 + pd0]
		*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16; *ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=p3;*ip++=p2;*ip++= 0;*ip++=p1; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16;
		// Set 6: [1,0,3,2,4,5,6,7,a,b,9,8,f,e,c,d + pb0],[e,f,d,c,8,9,a,b,2,3,1,0,7,6,4,5 + pa0]
		*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16; *ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16;
		// Set 7: [5,4,7,6,2,3,1,0,e,f,d,c,8,9,a,b + p80],[1,0,3,2,4,5,6,7,a,b,9,8,f,e,c,d + p90]
		*ip++=p5;*ip++=p4;*ip++=p7;*ip++=p6;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16; *ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16;
		// Set 8: [d,c,f,e,a,b,9,8,1,0,3,2,4,5,6,7 + p70],[5,4,7,6,2,3,1,0,e,f,d,c,8,9,a,b + p60]
		*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16; *ip++=p5;*ip++=p4;*ip++=p7;*ip++=p6;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16;
		// Set 9: [b,a,8,9,e,f,d,c,5,4,7,6,2,3,1,0 + p40],[d,c,f,e,a,b,9,8,1,0,3,2,4,5,6,7 + p50]
		*ip++=pb;*ip++=pa;*ip++=p8;*ip++=p9;*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p5;*ip++=p4;*ip++=p7;*ip++=p6;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16; *ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16;
		// Set A: [7,6,4,5,1,0,3,2,d,c,f,e,a,b,9,8 + p30],[b,a,8,9,e,f,d,c,5,4,7,6,2,3,1,0 + p20]
		*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8; for(j = i; j < i+16; j++) { dif_offsets[j] += p10; }; i += 16; *ip++=pb;*ip++=pa;*ip++=p8;*ip++=p9;*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p5;*ip++=p4;*ip++=p7;*ip++=p6;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0; for(j = i; j < i+16; j++) { dif_offsets[j] +=   0; }; i += 16;
	#endif
	}

/*...The radix-352 pass is here.	*/

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
	Twiddleless version arranges 32 sets of radix-11 DFT inputs as follows: 0 in upper left corner,
	decrement 32 (= 0x20) horizontally and 11 vertically, all indexing done (mod 352 = 0x160).

	Indexing in hex for clarity and using [evn|odd]0-10 notation in the rightmost column to flag reusable
	11-perms [in fact simple circular (0-10)-element shifts of evn = [000,140,120,100,0e0,0c0,0a0,080,060,040,020]
	and odd = [150,130,110,0f0,0d0,0b0,090,070,050,030,010].
		...Can auto-generate these by running test_fft_radix with TTYPE = 0;
	note this input-offset pattern is shared by DIF and DIT, but DIT layers a generalized bit-reversal atop it:

			DIF/DIT input-scramble array:			[000,140,...,020 = basic even-offset array]			row	leftward-circ-shift count
													[150,130,...,010 = basic  odd-offset array] 		idx	of basic even|odd arrays:
	  0,140,120,100, e0, c0, a0, 80, 60, 40, 20		000,140,120,100,0e0,0c0,0a0,080,060,040,020 + p0	0	[evn0 ] + p0
	155,135,115, f5, d5, b5, 95, 75, 55, 35, 15		150,130,110,0f0,0d0,0b0,090,070,050,030,010 + p5	1	[odd0 ] + p5
	14a,12a,10a, ea, ca, aa, 8a, 6a, 4a, 2a,  a		140,120,100,0e0,0c0,0a0,080,060,040,020,000 + pa	2	[evn1 ] + pa
	13f,11f, ff, df, bf, 9f, 7f, 5f, 3f, 1f,15f		130,110,0f0,0d0,0b0,090,070,050,030,010,150 + pf	3	[odd1 ] + pf
	134,114, f4, d4, b4, 94, 74, 54, 34, 14,154		130,110,0f0,0d0,0b0,090,070,050,030,010,150 + p4	4	[odd1 ] + p4
	129,109, e9, c9, a9, 89, 69, 49, 29,  9,149		120,100,0e0,0c0,0a0,080,060,040,020,000,140 + p9	5	[evn2 ] + p9
	11e, fe, de, be, 9e, 7e, 5e, 3e, 1e,15e,13e		110,0f0,0d0,0b0,090,070,050,030,010,150,130 + pe	6	[odd2 ] + pe
	113, f3, d3, b3, 93, 73, 53, 33, 13,153,133		110,0f0,0d0,0b0,090,070,050,030,010,150,130 + p3	7	[odd2 ] + p3
	108, e8, c8, a8, 88, 68, 48, 28,  8,148,128		100,0e0,0c0,0a0,080,060,040,020,000,140,120 + p8	8	[evn3 ] + p8
	 fd, dd, bd, 9d, 7d, 5d, 3d, 1d,15d,13d,11d		0f0,0d0,0b0,090,070,050,030,010,150,130,110 + pd	9	[odd3 ] + pd
	 f2, d2, b2, 92, 72, 52, 32, 12,152,132,112		0f0,0d0,0b0,090,070,050,030,010,150,130,110 + p2	a	[odd3 ] + p2
	 e7, c7, a7, 87, 67, 47, 27,  7,147,127,107		0e0,0c0,0a0,080,060,040,020,000,140,120,100 + p7	b	[evn4 ] + p7
	 dc, bc, 9c, 7c, 5c, 3c, 1c,15c,13c,11c, fc		0d0,0b0,090,070,050,030,010,150,130,110,0f0 + pc	c	[odd4 ] + pc
	 d1, b1, 91, 71, 51, 31, 11,151,131,111, f1		0d0,0b0,090,070,050,030,010,150,130,110,0f0 + p1	d	[odd4 ] + p1
	 c6, a6, 86, 66, 46, 26,  6,146,126,106, e6		0c0,0a0,080,060,040,020,000,140,120,100,0e0 + p6	e	[evn5 ] + p6
	 bb, 9b, 7b, 5b, 3b, 1b,15b,13b,11b, fb, db		0b0,090,070,050,030,010,150,130,110,0f0,0d0 + pb	f	[odd5 ] + pb
	 b0, 90, 70, 50, 30, 10,150,130,110, f0, d0	=	0b0,090,070,050,030,010,150,130,110,0f0,0d0 + p0	10	[odd5 ] + p0 <<< p0-f pattern repeats here
	 a5, 85, 65, 45, 25,  5,145,125,105, e5, c5		0a0,080,060,040,020,000,140,120,100,0e0,0c0 + p5	11	[evn6 ] + p5
	 9a, 7a, 5a, 3a, 1a,15a,13a,11a, fa, da, ba		090,070,050,030,010,150,130,110,0f0,0d0,0b0 + pa	12	[odd6 ] + pa
	 8f, 6f, 4f, 2f,  f,14f,12f,10f, ef, cf, af		080,060,040,020,000,140,120,100,0e0,0c0,0a0 + pf	13	[evn7 ] + pf
	 84, 64, 44, 24,  4,144,124,104, e4, c4, a4		080,060,040,020,000,140,120,100,0e0,0c0,0a0 + p4	14	[evn7 ] + p4
	 79, 59, 39, 19,159,139,119, f9, d9, b9, 99		070,050,030,010,150,130,110,0f0,0d0,0b0,090 + p9	15	[odd7 ] + p9
	 6e, 4e, 2e,  e,14e,12e,10e, ee, ce, ae, 8e		060,040,020,000,140,120,100,0e0,0c0,0a0,080 + pe	16	[evn8 ] + pe
	 63, 43, 23,  3,143,123,103, e3, c3, a3, 83		060,040,020,000,140,120,100,0e0,0c0,0a0,080 + p3	17	[evn8 ] + p3
	 58, 38, 18,158,138,118, f8, d8, b8, 98, 78		050,030,010,150,130,110,0f0,0d0,0b0,090,070 + p8	18	[odd8 ] + p8
	 4d, 2d,  d,14d,12d,10d, ed, cd, ad, 8d, 6d		040,020,000,140,120,100,0e0,0c0,0a0,080,060 + pd	19	[evn9 ] + pd
	 42, 22,  2,142,122,102, e2, c2, a2, 82, 62		040,020,000,140,120,100,0e0,0c0,0a0,080,060 + p2	1a	[evn9 ] + p2
	 37, 17,157,137,117, f7, d7, b7, 97, 77, 57		030,010,150,130,110,0f0,0d0,0b0,090,070,050 + p7	1b	[odd9 ] + p7
	 2c,  c,14c,12c,10c, ec, cc, ac, 8c, 6c, 4c		020,000,140,120,100,0e0,0c0,0a0,080,060,040 + pc	1c	[evn10] + pc
	 21,  1,141,121,101, e1, c1, a1, 81, 61, 41		020,000,140,120,100,0e0,0c0,0a0,080,060,040 + p1	1d	[evn10] + p1
	 16,156,136,116, f6, d6, b6, 96, 76, 56, 36		010,150,130,110,0f0,0d0,0b0,090,070,050,030 + p6	1e	[odd10] + p6
	  b,14b,12b,10b, eb, cb, ab, 8b, 6b, 4b, 2b		000,140,120,100,0e0,0c0,0a0,080,060,040,020 + pb	1f	[evn0 ] + pb
	*/
	/*...gather the needed data (352 64-bit complex) and do 32 radix-11 transforms...*/
		tptr = t;
		for(i = 0; i < 32; i++) {
			int k = dif_p20_lo_offset[i];
			// Extract index (in [0-10]) into circ-shift array used for high parts of p-mults. The [0-10] value is
			// in low 4 bits of k; the "which length-21 half of the dif_p20_cperms array?" selector is via (k < 0):
			int ic = ((-(k < 0)) & 21)	// +/- sign on k puts us into lower/upper half of the cshift array (base index 0/21)
						+ (k & 0xf);	// ...and low 4 bits give the element index w.r.to the array-half in question.
			int k0 = dif_p20_cperms[ic], k1 = dif_p20_cperms[ic+1], k2 = dif_p20_cperms[ic+2], k3 = dif_p20_cperms[ic+3], k4 = dif_p20_cperms[ic+4], k5 = dif_p20_cperms[ic+5], k6 = dif_p20_cperms[ic+6], k7 = dif_p20_cperms[ic+7], k8 = dif_p20_cperms[ic+8], k9 = dif_p20_cperms[ic+9], ka = dif_p20_cperms[ic+10];
			// Extract Low part, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-11 DFTs:
			k = (k & 0x7fffffff) >> 4;
			jt = j1+k; jp = j2+k;
			RADIX_11_DFT(
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],
				tptr->re,tptr->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x100)->re,(tptr+0x100)->im,(tptr+0x120)->re,(tptr+0x120)->im,(tptr+0x140)->re,(tptr+0x140)->im,
				a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9
			);	tptr++;
		}

	/*...and now do 11 radix-32 transforms.
	Required output index permutation =
	  0,  1,  2,  3,  5,  4,  7,  6,  b,  a,  8,  9,  e,  f,  d,  c, 17, 16, 14, 15, 11, 10, 13, 12, 1d, 1c, 1f, 1e, 1a, 1b, 19, 18,
	148,149,14a,14b,14d,14c,14f,14e,147,146,144,145,141,140,143,142,15f,15e,15c,15d,159,158,15b,15a,153,152,150,151,156,157,155,154,
	134,135,136,137,133,132,130,131,13f,13e,13c,13d,139,138,13b,13a,128,129,12a,12b,12d,12c,12f,12e,127,126,124,125,121,120,123,122,
	102,103,101,100,107,106,104,105,108,109,10a,10b,10d,10c,10f,10e,114,115,116,117,113,112,110,111,11f,11e,11c,11d,119,118,11b,11a,
	 fa, fb, f9, f8, ff, fe, fc, fd, f4, f5, f6, f7, f3, f2, f0, f1, e2, e3, e1, e0, e7, e6, e4, e5, e8, e9, ea, eb, ed, ec, ef, ee,
	 ce, cf, cd, cc, c8, c9, ca, cb, c2, c3, c1, c0, c7, c6, c4, c5, da, db, d9, d8, df, de, dc, dd, d4, d5, d6, d7, d3, d2, d0, d1,
	 b1, b0, b3, b2, b4, b5, b6, b7, ba, bb, b9, b8, bf, be, bc, bd, ae, af, ad, ac, a8, a9, aa, ab, a2, a3, a1, a0, a7, a6, a4, a5,
	 85, 84, 87, 86, 82, 83, 81, 80, 8e, 8f, 8d, 8c, 88, 89, 8a, 8b, 91, 90, 93, 92, 94, 95, 96, 97, 9a, 9b, 99, 98, 9f, 9e, 9c, 9d,
	 7d, 7c, 7f, 7e, 7a, 7b, 79, 78, 71, 70, 73, 72, 74, 75, 76, 77, 65, 64, 67, 66, 62, 63, 61, 60, 6e, 6f, 6d, 6c, 68, 69, 6a, 6b,
	 4b, 4a, 48, 49, 4e, 4f, 4d, 4c, 45, 44, 47, 46, 42, 43, 41, 40, 5d, 5c, 5f, 5e, 5a, 5b, 59, 58, 51, 50, 53, 52, 54, 55, 56, 57,
	 37, 36, 34, 35, 31, 30, 33, 32, 3d, 3c, 3f, 3e, 3a, 3b, 39, 38, 2b, 2a, 28, 29, 2e, 2f, 2d, 2c, 25, 24, 27, 26, 22, 23, 21, 20,
	=
	[0,1,2,3,5,4,7,6,b,a,8,9,e,f,d,c + p00],[7,6,4,5,1,0,3,2,d,c,f,e,a,b,9,8 + p10]
	[8,9,a,b,d,c,f,e,7,6,4,5,1,0,3,2 +p140],[f,e,c,d,9,8,b,a,3,2,0,1,6,7,5,4 +p150]
	[4,5,6,7,3,2,0,1,f,e,c,d,9,8,b,a +p130],[8,9,a,b,d,c,f,e,7,6,4,5,1,0,3,2 +p120]
	[2,3,1,0,7,6,4,5,8,9,a,b,d,c,f,e +p100],[4,5,6,7,3,2,0,1,f,e,c,d,9,8,b,a +p110]
	[a,b,9,8,f,e,c,d,4,5,6,7,3,2,0,1 + pf0],[2,3,1,0,7,6,4,5,8,9,a,b,d,c,f,e + pe0]
	[e,f,d,c,8,9,a,b,2,3,1,0,7,6,4,5 + pc0],[a,b,9,8,f,e,c,d,4,5,6,7,3,2,0,1 + pd0]
	[1,0,3,2,4,5,6,7,a,b,9,8,f,e,c,d + pb0],[e,f,d,c,8,9,a,b,2,3,1,0,7,6,4,5 + pa0]
	[5,4,7,6,2,3,1,0,e,f,d,c,8,9,a,b + p80],[1,0,3,2,4,5,6,7,a,b,9,8,f,e,c,d + p90]
	[d,c,f,e,a,b,9,8,1,0,3,2,4,5,6,7 + p70],[5,4,7,6,2,3,1,0,e,f,d,c,8,9,a,b + p60]
	[b,a,8,9,e,f,d,c,5,4,7,6,2,3,1,0 + p40],[d,c,f,e,a,b,9,8,1,0,3,2,4,5,6,7 + p50]
	[7,6,4,5,1,0,3,2,d,c,f,e,a,b,9,8 + p30],[b,a,8,9,e,f,d,c,5,4,7,6,2,3,1,0 + p20]
	*/
		tptr = t;
		for(i = 0; i < ODD_RADIX; i++) {
			jt = j1+dif_phi[i]; RADIX_32_DIF((double *)tptr,t_offsets,1, (a+jt),dif_offsets+(i<<5),RE_IM_STRIDE);	tptr += 32;
		}
	}
}

/***************/

void radix352_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-352 complex DIT FFT pass on the data in the length-N real vector A.
*/
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
	int i,j,j1,j2,jt,jp, *ip;
	// p-indexing is hexadecimal here:
	static int NDIVR,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
		,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,p100,p110,p120,p130,p140,p150, first_entry=TRUE;
	static int t_offsets[32], dit_offsets[RADIX];
	// Need storage for 2 circular-shifts perms of a basic 11-vector, with shift count in [0,10] that means 2*21 elts:
	static int dit_p20_cperms[42], dit_p20_lo_offset[32], dit_phi[ODD_RADIX];
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
		p3 = p2 + p1;	p30 = p20 + p10;	p120 = p110 + p10;
		p4 = p3 + p1;	p40 = p30 + p10;	p130 = p120 + p10;
		p5 = p4 + p1;	p50 = p40 + p10;	p140 = p130 + p10;
		p6 = p5 + p1;	p60 = p50 + p10;	p150 = p140 + p10;
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
													p120 += ( (p120 >> DAT_BITS) << PAD_BITS );
													p130 += ( (p130 >> DAT_BITS) << PAD_BITS );
													p140 += ( (p140 >> DAT_BITS) << PAD_BITS );
													p150 += ( (p150 >> DAT_BITS) << PAD_BITS );
		i = 0;
		dit_phi[i++] =   0;
		dit_phi[i++] = p20;
		dit_phi[i++] = p40;
		dit_phi[i++] = p60;
		dit_phi[i++] = p80;
		dit_phi[i++] = pa0;
		dit_phi[i++] = pc0;
		dit_phi[i++] = pe0;
		dit_phi[i++] =p100;
		dit_phi[i++] =p120;
		dit_phi[i++] =p140;

		// Set array offsets for radix-32 DFT in/outputs:
		// t_offsets w.r.to: t-array, same for all 11 DFTs:
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

		// Init storage for 2 circular-shifts perms of a basic 11-vector, with shift count in [0,10] that means 2*22 = 42 elts:
		i = 0;
		// Even multiples of p10 cshift array:
		dit_p20_cperms[i++] = 0;
		dit_p20_cperms[i++] = p140;
		dit_p20_cperms[i++] = p120;
		dit_p20_cperms[i++] = p100;
		dit_p20_cperms[i++] = pe0;
		dit_p20_cperms[i++] = pc0;
		dit_p20_cperms[i++] = pa0;
		dit_p20_cperms[i++] = p80;
		dit_p20_cperms[i++] = p60;
		dit_p20_cperms[i++] = p40;
		dit_p20_cperms[i++] = p20;
		while(i < 2*ODD_RADIX-1) {
			dit_p20_cperms[i] = dit_p20_cperms[i - ODD_RADIX]; ++i;
		}
		// Odd multiples of p10 cshift array:
		dit_p20_cperms[i++] = p150;
		dit_p20_cperms[i++] = p130;
		dit_p20_cperms[i++] = p110;
		dit_p20_cperms[i++] = pf0;
		dit_p20_cperms[i++] = pd0;
		dit_p20_cperms[i++] = pb0;
		dit_p20_cperms[i++] = p90;
		dit_p20_cperms[i++] = p70;
		dit_p20_cperms[i++] = p50;
		dit_p20_cperms[i++] = p30;
		dit_p20_cperms[i++] = p10;
		while(i < 4*ODD_RADIX-2) {
			dit_p20_cperms[i] = dit_p20_cperms[i - ODD_RADIX]; ++i;
		}

		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-11 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-8 index, bits <4:30> for the p0-f:
		i = 0;
		dit_p20_lo_offset[i++] = (( 0 << 4) + 0 );
		dit_p20_lo_offset[i++] = ((pf << 4) + 1 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pe << 4) + 2 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pd << 4) + 3 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pc << 4) + 4 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pb << 4) + 5 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pa << 4) + 6 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p9 << 4) + 7 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p8 << 4) + 8 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p7 << 4) + 9 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p6 << 4) + 10) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p5 << 4) + 0 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p4 << 4) + 1 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p3 << 4) + 2 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p2 << 4) + 3 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((p1 << 4) + 4 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = (( 0 << 4) + 5 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[i++] = ((pf << 4) + 7 );
		dit_p20_lo_offset[i++] = ((pe << 4) + 8 );
		dit_p20_lo_offset[i++] = ((pd << 4) + 9 );
		dit_p20_lo_offset[i++] = ((pc << 4) + 10);
		dit_p20_lo_offset[i++] = ((pb << 4) + 0 );
		dit_p20_lo_offset[i++] = ((pa << 4) + 1 );
		dit_p20_lo_offset[i++] = ((p9 << 4) + 2 );
		dit_p20_lo_offset[i++] = ((p8 << 4) + 3 );
		dit_p20_lo_offset[i++] = ((p7 << 4) + 4 );
		dit_p20_lo_offset[i++] = ((p6 << 4) + 5 );
		dit_p20_lo_offset[i++] = ((p5 << 4) + 6 );
		dit_p20_lo_offset[i++] = ((p4 << 4) + 7 );
		dit_p20_lo_offset[i++] = ((p3 << 4) + 8 );
		dit_p20_lo_offset[i++] = ((p2 << 4) + 9 );
		dit_p20_lo_offset[i++] = ((p1 << 4) + 10);

	// dit_offsets are w.r.to a-array, need 11 distinct sets of these, one for each radix-32 DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		i = 0; ip = dit_offsets;	// In the following, the 'big' += p[>= 20] offsets already in dif_phi[], only care about += p0|p10:
		// Set 0: [0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00],[f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10]
		*ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2;*ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4;*ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc;*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16; *ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc;*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8;*ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4;*ip++=p3;*ip++=p2;*ip++=p1;*ip++= 0; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16;
		// Set 1: [7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + p30],[7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + p20]
		*ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4;*ip++=p3;*ip++=p2;*ip++=p1;*ip++= 0;*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8;*ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16; *ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4;*ip++=p3;*ip++=p2;*ip++=p1;*ip++= 0;*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8;*ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16;
		// Set 2: [b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p40],[3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p50]
		*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8;*ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf;*ip++=p3;*ip++=p2;*ip++=p1;*ip++= 0;*ip++=p5;*ip++=p4;*ip++=p6;*ip++=p7; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16; *ip++=p3;*ip++=p2;*ip++=p1;*ip++= 0;*ip++=p5;*ip++=p4;*ip++=p6;*ip++=p7;*ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf;*ip++=p9;*ip++=p8;*ip++=pa;*ip++=pb; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16;
		// Set 3: [d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + p70],[d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + p60]
		*ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf;*ip++=p9;*ip++=p8;*ip++=pa;*ip++=pb;*ip++=p5;*ip++=p4;*ip++=p6;*ip++=p7;*ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16; *ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf;*ip++=p9;*ip++=p8;*ip++=pa;*ip++=pb;*ip++=p5;*ip++=p4;*ip++=p6;*ip++=p7;*ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16;
		// Set 4: [5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p80],[9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p90]
		*ip++=p5;*ip++=p4;*ip++=p6;*ip++=p7;*ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3;*ip++=p9;*ip++=p8;*ip++=pa;*ip++=pb;*ip++=pe;*ip++=pf;*ip++=pc;*ip++=pd; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16; *ip++=p9;*ip++=p8;*ip++=pa;*ip++=pb;*ip++=pe;*ip++=pf;*ip++=pc;*ip++=pd;*ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3;*ip++=p6;*ip++=p7;*ip++=p4;*ip++=p5; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16;
		// Set 5: [1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pb0],[1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pa0]
		*ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3;*ip++=p6;*ip++=p7;*ip++=p4;*ip++=p5;*ip++=pe;*ip++=pf;*ip++=pc;*ip++=pd;*ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16; *ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3;*ip++=p6;*ip++=p7;*ip++=p4;*ip++=p5;*ip++=pe;*ip++=pf;*ip++=pc;*ip++=pd;*ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16;
		// Set 6: [e,f,c,d,a,b,8,9,6,7,4,5,2,3,0,1 + pc0],[6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pd0]
		*ip++=pe;*ip++=pf;*ip++=pc;*ip++=pd;*ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9;*ip++=p6;*ip++=p7;*ip++=p4;*ip++=p5;*ip++=p2;*ip++=p3;*ip++= 0;*ip++=p1; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16; *ip++=p6;*ip++=p7;*ip++=p4;*ip++=p5;*ip++=p2;*ip++=p3;*ip++= 0;*ip++=p1;*ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9;*ip++=pc;*ip++=pd;*ip++=pf;*ip++=pe; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16;
		// Set 7: [a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + pf0],[a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + pe0]
		*ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9;*ip++=pc;*ip++=pd;*ip++=pf;*ip++=pe;*ip++=p2;*ip++=p3;*ip++= 0;*ip++=p1;*ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16; *ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9;*ip++=pc;*ip++=pd;*ip++=pf;*ip++=pe;*ip++=p2;*ip++=p3;*ip++= 0;*ip++=p1;*ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16;
		// Set 8: [2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a +p100],[c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 +p110]
		*ip++=p2;*ip++=p3;*ip++= 0;*ip++=p1;*ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6;*ip++=pc;*ip++=pd;*ip++=pf;*ip++=pe;*ip++=p8;*ip++=p9;*ip++=pb;*ip++=pa; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16; *ip++=pc;*ip++=pd;*ip++=pf;*ip++=pe;*ip++=p8;*ip++=p9;*ip++=pb;*ip++=pa;*ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6;*ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16;
		// Set 9: [4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c +p130],[4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c +p120]
		*ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6;*ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2;*ip++=p8;*ip++=p9;*ip++=pb;*ip++=pa;*ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16; *ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6;*ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2;*ip++=p8;*ip++=p9;*ip++=pb;*ip++=pa;*ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16;
		// Set A: [8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 +p140],[0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 +p150]
		*ip++=p8;*ip++=p9;*ip++=pb;*ip++=pa;*ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc;*ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2;*ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4; for(j = i; j < i+16; j++) { dit_offsets[j] +=   0; }; i += 16; *ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2;*ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4;*ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc;*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8; for(j = i; j < i+16; j++) { dit_offsets[j] += p10; }; i += 16;
	}

/*...The radix-352 pass is here.	*/

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

	/*...gather the needed data (352 64-bit complex) and do 11 radix-32 transforms.

	Twiddleless version uses same linear-index-vector-form permutation as in DIF -
	Remember, inputs to DIT are bit-reversed, so use output of test_fft_radix() with
	TTYPE=0 to auto-generate needed input-index permutation:

	Combined DIT input-scramble array =
	  0,  1,  3,  2,  7,  6,  5,  4,  f,  e,  d,  c,  b,  a,  9,  8, 1f, 1e, 1d, 1c, 1b, 1a, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10,
	 37, 36, 35, 34, 33, 32, 31, 30, 3b, 3a, 39, 38, 3d, 3c, 3e, 3f, 27, 26, 25, 24, 23, 22, 21, 20, 2b, 2a, 29, 28, 2d, 2c, 2e, 2f,
	 4b, 4a, 49, 48, 4d, 4c, 4e, 4f, 43, 42, 41, 40, 45, 44, 46, 47, 53, 52, 51, 50, 55, 54, 56, 57, 5d, 5c, 5e, 5f, 59, 58, 5a, 5b,
	 7d, 7c, 7e, 7f, 79, 78, 7a, 7b, 75, 74, 76, 77, 71, 70, 72, 73, 6d, 6c, 6e, 6f, 69, 68, 6a, 6b, 65, 64, 66, 67, 61, 60, 62, 63,
	 85, 84, 86, 87, 81, 80, 82, 83, 89, 88, 8a, 8b, 8e, 8f, 8c, 8d, 99, 98, 9a, 9b, 9e, 9f, 9c, 9d, 91, 90, 92, 93, 96, 97, 94, 95,
	 b1, b0, b2, b3, b6, b7, b4, b5, be, bf, bc, bd, ba, bb, b8, b9, a1, a0, a2, a3, a6, a7, a4, a5, ae, af, ac, ad, aa, ab, a8, a9,
	 ce, cf, cc, cd, ca, cb, c8, c9, c6, c7, c4, c5, c2, c3, c0, c1, d6, d7, d4, d5, d2, d3, d0, d1, da, db, d8, d9, dc, dd, df, de,
	 fa, fb, f8, f9, fc, fd, ff, fe, f2, f3, f0, f1, f4, f5, f7, f6, ea, eb, e8, e9, ec, ed, ef, ee, e2, e3, e0, e1, e4, e5, e7, e6,
	102,103,100,101,104,105,107,106,10c,10d,10f,10e,108,109,10b,10a,11c,11d,11f,11e,118,119,11b,11a,114,115,117,116,110,111,113,112,
	134,135,137,136,130,131,133,132,138,139,13b,13a,13f,13e,13d,13c,124,125,127,126,120,121,123,122,128,129,12b,12a,12f,12e,12d,12c,
	148,149,14b,14a,14f,14e,14d,14c,140,141,143,142,147,146,145,144,150,151,153,152,157,156,155,154,15f,15e,15d,15c,15b,15a,159,158,
	=
	[0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00],[f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10]
	[7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + p30],[7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + p20]
	[b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p40],[3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p50]
	[d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + p70],[d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + p60]
	[5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p80],[9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p90]
	[1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pb0],[1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pa0]
	[e,f,c,d,a,b,8,9,6,7,4,5,2,3,0,1 + pc0],[6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pd0]
	[a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + pf0],[a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + pe0]
	[2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a +p100],[c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 +p110]
	[4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c +p130],[4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c +p120]
	[8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 +p140],[0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 +p150]
	These offset patterns are encoded by the combination of the dit_phi[] (p20-multiples) and dit_offsets (mod-p20 parts) arrays.
	*/
		tptr = t;
		for(i = 0; i < ODD_RADIX; i++) {
			jt = j1+dit_phi[i]; RADIX_32_DIT((a+jt),dit_offsets+(i<<5),RE_IM_STRIDE, (double *)tptr,t_offsets,1);	tptr += 32;
		}

	/*...and now do 32 radix-11 transforms, with the columns of t*[r,i] output pairs in the above 11x radix-32 set now acting as input rows.
	Since our first-look oindex ordering was +p0-a for each radix-11 and incrementing += p11 between those DFTs,
	simply break resulting o-index perm generated by the test_fft_radix code into 11-element rows.
	Indexing in hex for clarity and using [evn|odd]0-a notation in the rightmost column to flag reusable 11-perms, in fact
	simple circular <<(0-a) shifts of the basic patterns:

	Required output index permutation =
			DIF/DIT input-scramble array:			[000,140,...,020 = basic even-offset array]			row	leftward-circ-shift count
													[150,130,...,010 = basic  odd-offset array] 		idx	of basic even|odd arrays:
	  0,140,120,100, e0, c0, a0, 80, 60, 40, 20		  0,140,120,100, e0, c0, a0, 80, 60, 40, 20 + p0	0	[evn0 ] + p0
	13f,11f, ff, df, bf, 9f, 7f, 5f, 3f, 1f,15f		130,110, f0, d0, b0, 90, 70, 50, 30, 10,150 + pf	1	[odd1 ] + pf
	11e, fe, de, be, 9e, 7e, 5e, 3e, 1e,15e,13e		110, f0, d0, b0, 90, 70, 50, 30, 10,150,130 + pe	2	[odd2 ] + pe
	 fd, dd, bd, 9d, 7d, 5d, 3d, 1d,15d,13d,11d		 f0, d0, b0, 90, 70, 50, 30, 10,150,130,110 + pd	3	[odd3 ] + pd
	 dc, bc, 9c, 7c, 5c, 3c, 1c,15c,13c,11c, fc		 d0, b0, 90, 70, 50, 30, 10,150,130,110, f0 + pc	4	[odd4 ] + pc
	 bb, 9b, 7b, 5b, 3b, 1b,15b,13b,11b, fb, db		 b0, 90, 70, 50, 30, 10,150,130,110, f0, d0 + pb	5	[odd5 ] + pb
	 9a, 7a, 5a, 3a, 1a,15a,13a,11a, fa, da, ba		 90, 70, 50, 30, 10,150,130,110, f0, d0, b0 + pa	6	[odd6 ] + pa
	 79, 59, 39, 19,159,139,119, f9, d9, b9, 99		 70, 50, 30, 10,150,130,110, f0, d0, b0, 90 + p9	7	[odd7 ] + p9
	 58, 38, 18,158,138,118, f8, d8, b8, 98, 78		 50, 30, 10,150,130,110, f0, d0, b0, 90, 70 + p8	8	[odd8 ] + p8
	 37, 17,157,137,117, f7, d7, b7, 97, 77, 57		 30, 10,150,130,110, f0, d0, b0, 90, 70, 50 + p7	9	[odd9 ] + p7
	 16,156,136,116, f6, d6, b6, 96, 76, 56, 36		 10,150,130,110, f0, d0, b0, 90, 70, 50, 30 + p6	a	[odd10] + p6
	155,135,115, f5, d5, b5, 95, 75, 55, 35, 15		150,130,110, f0, d0, b0, 90, 70, 50, 30, 10 + p5	b	[odd0 ] + p5
	134,114, f4, d4, b4, 94, 74, 54, 34, 14,154		130,110, f0, d0, b0, 90, 70, 50, 30, 10,150 + p4	c	[odd1 ] + p4
	113, f3, d3, b3, 93, 73, 53, 33, 13,153,133		110, f0, d0, b0, 90, 70, 50, 30, 10,150,130 + p3	d	[odd2 ] + p3
	 f2, d2, b2, 92, 72, 52, 32, 12,152,132,112		 f0, d0, b0, 90, 70, 50, 30, 10,150,130,110 + p2	e	[odd3 ] + p2
	 d1, b1, 91, 71, 51, 31, 11,151,131,111, f1		 d0, b0, 90, 70, 50, 30, 10,150,130,110, f0 + p1	f	[odd4 ] + p1
	 b0, 90, 70, 50, 30, 10,150,130,110, f0, d0		 b0, 90, 70, 50, 30, 10,150,130,110, f0, d0 + p0	10	[odd5 ] + p0 <<< p0-f pattern repeats here
	 8f, 6f, 4f, 2f,  f,14f,12f,10f, ef, cf, af		 80, 60, 40, 20,  0,140,120,100, e0, c0, a0 + pf	11	[evn7 ] + pf
	 6e, 4e, 2e,  e,14e,12e,10e, ee, ce, ae, 8e		 60, 40, 20,  0,140,120,100, e0, c0, a0, 80 + pe	12	[evn8 ] + pe
	 4d, 2d,  d,14d,12d,10d, ed, cd, ad, 8d, 6d		 40, 20,  0,140,120,100, e0, c0, a0, 80, 60 + pd	13	[evn9 ] + pd
	 2c,  c,14c,12c,10c, ec, cc, ac, 8c, 6c, 4c		 20,  0,140,120,100, e0, c0, a0, 80, 60, 40 + pc	14	[evn10] + pc
	  b,14b,12b,10b, eb, cb, ab, 8b, 6b, 4b, 2b		  0,140,120,100, e0, c0, a0, 80, 60, 40, 20 + pb	15	[evn0 ] + pb
	14a,12a,10a, ea, ca, aa, 8a, 6a, 4a, 2a,  a		140,120,100, e0, c0, a0, 80, 60, 40, 20,  0 + pa	16	[evn1 ] + pa
	129,109, e9, c9, a9, 89, 69, 49, 29,  9,149		120,100, e0, c0, a0, 80, 60, 40, 20,  0,140 + p9	17	[evn2 ] + p9
	108, e8, c8, a8, 88, 68, 48, 28,  8,148,128		100, e0, c0, a0, 80, 60, 40, 20,  0,140,120 + p8	18	[evn3 ] + p8
	 e7, c7, a7, 87, 67, 47, 27,  7,147,127,107		 e0, c0, a0, 80, 60, 40, 20,  0,140,120,100 + p7	19	[evn4 ] + p7
	 c6, a6, 86, 66, 46, 26,  6,146,126,106, e6		 c0, a0, 80, 60, 40, 20,  0,140,120,100, e0 + p6	1a	[evn5 ] + p6
	 a5, 85, 65, 45, 25,  5,145,125,105, e5, c5		 a0, 80, 60, 40, 20,  0,140,120,100, e0, c0 + p5	1b	[evn6 ] + p5
	 84, 64, 44, 24,  4,144,124,104, e4, c4, a4		 80, 60, 40, 20,  0,140,120,100, e0, c0, a0 + p4	1c	[evn7 ] + p4
	 63, 43, 23,  3,143,123,103, e3, c3, a3, 83		 60, 40, 20,  0,140,120,100, e0, c0, a0, 80 + p3	1d	[evn8 ] + p3
	 42, 22,  2,142,122,102, e2, c2, a2, 82, 62		 40, 20,  0,140,120,100, e0, c0, a0, 80, 60 + p2	1e	[evn9 ] + p2
	 21,  1,141,121,101, e1, c1, a1, 81, 61, 41		 20,  0,140,120,100, e0, c0, a0, 80, 60, 40 + p1	1f	[evn10] + p1
	*/
		tptr = t;
		for(i = 0; i < 32; i++) {
		#if 0	// First-look linear indexing used to allow test_fft_radix to extract the needed operm:
			int k = i*22;
			int k0=0,k1=2,k2=4,k3=6,k4=8,k5=10,k6=12,k7=14,k8=16,k9=18,ka=20;
		#else
			int k = dit_p20_lo_offset[i];
			// Extract index (in [0-10]) into circ-shift array used for high parts of p-mults. The [0-10] value is
			// in low 4 bits of k; the "which length-21 half of the dit_p20_cperms array?" selector is via (k < 0):
			int ic = ((-(k < 0)) & 21)	// +/- sign on k puts us into lower/upper half of the cshift array (base index 0/21)
						+ (k & 0xf);	// ...and low 4 bits give the element index w.r.to the array-half in question.
			int k0 = dit_p20_cperms[ic], k1 = dit_p20_cperms[ic+1], k2 = dit_p20_cperms[ic+2], k3 = dit_p20_cperms[ic+3], k4 = dit_p20_cperms[ic+4], k5 = dit_p20_cperms[ic+5], k6 = dit_p20_cperms[ic+6], k7 = dit_p20_cperms[ic+7], k8 = dit_p20_cperms[ic+8], k9 = dit_p20_cperms[ic+9], ka = dit_p20_cperms[ic+10];
			// Extract Low part, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-11 DFTs:
			k = (k & 0x7fffffff) >> 4;
		#endif
			jt = j1+k; jp = j2+k;
			RADIX_11_DFT(
				tptr->re,tptr->im,(tptr+0x20)->re,(tptr+0x20)->im,(tptr+0x40)->re,(tptr+0x40)->im,(tptr+0x60)->re,(tptr+0x60)->im,(tptr+0x80)->re,(tptr+0x80)->im,(tptr+0xa0)->re,(tptr+0xa0)->im,(tptr+0xc0)->re,(tptr+0xc0)->im,(tptr+0xe0)->re,(tptr+0xe0)->im,(tptr+0x100)->re,(tptr+0x100)->im,(tptr+0x120)->re,(tptr+0x120)->im,(tptr+0x140)->re,(tptr+0x140)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],
				a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b0,b1,b2,b3,b4,b5,b6,b7,b8,b9
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
	cy352_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr;
	#ifdef USE_SSE2
		const int pfetch_dist = PFETCH_DIST;
	#endif
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf
				,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,p100,p110,p120,p130,p140,p150;
	// Shared DIF+DIT:
		int poff[RADIX>>2];	// Store [RADIX/4] mults of p04 offset for loop control
	#ifndef USE_SSE2
		int t_offsets[32];
	#endif
		// Need storage for 2 circular-shifts perms of a basic 11-vector, with shift count in [0,10] that means 2*21 elts:
		int dif_offsets[RADIX], dif_p20_cperms[42], dif_p20_lo_offset[32], dif_phi[ODD_RADIX];
		int dit_offsets[RADIX], dit_p20_cperms[42], dit_p20_lo_offset[32], dit_phi[ODD_RADIX];
	#ifdef USE_SSE2
		const uint32 dft11_offs[11] = {0,0x40<<L2_SZ_VD,0x80<<L2_SZ_VD,0xc0<<L2_SZ_VD,0x100<<L2_SZ_VD,0x140<<L2_SZ_VD,0x180<<L2_SZ_VD,0x1c0<<L2_SZ_VD,0x200<<L2_SZ_VD,0x240<<L2_SZ_VD,0x280<<L2_SZ_VD}, *dft11_offptr = &(dft11_offs[0]);
	#endif
		int j,j1,l, *ip;
	#ifndef USE_SSE2
		int j2,jt,jp,ntmp;
	#endif
		// incr = Carry-chain wts-multipliers recurrence length, which must divide
		// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 22|44|88 for avx512,avx,sse, respectively.
		// But fixed-incr too restrictive here, so 'divide 22|44|88 into pieces' via increment-array whose elts sum to 22|44|88:
	#ifdef USE_SSE2
		const int *incr,*inc_arr;
	  #ifdef USE_AVX512	// Have no specialized HIACC carry macro in AVX-512 and ARMv8 SIMD, so these get an "goes to 11" in LOACC mode via an incr_hiacc[] array:
		const int incr_long[] = {11,11}, incr_med[] = {6,5,6,5}, incr_short[] = {4,3,4,4,3,4}, incr_hiacc[] = {2,2,2,2,2,2,2,2,2,2,2};
	  #elif defined(USE_AVX)
		const int incr_long[] = {11,11,11,11}, incr_med[] = {6,5,6,5,6,5,6,5}, incr_short[] = {4,4,4,4,4,4,4,4,4,4,4}, incr_hiacc[] = {0};
	  #elif defined(USE_ARM_V8_SIMD)
		const int incr_long[] = {11,11,11,11,11,11,11,11}, incr_med[] = {6,5,6,5,6,5,6,5,6,5,6,5,6,5,6,5}, incr_short[] = {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4}, incr_hiacc[] = {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
	  #else
		const int incr_long[] = {11,11,11,11,11,11,11,11}, incr_med[] = {6,5,6,5,6,5,6,5,6,5,6,5,6,5,6,5}, incr_short[] = {4,3,4,4,3,4,4,3,4,4,3,4,4,3,4,4,3,4,4,3,4,4,3,4}, incr_hiacc[] = {0};
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
	#endif // !USE_SSE2

	#ifndef USE_AVX
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
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
	#endif
	#ifndef USE_SSE2
		double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
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
		double *add0, *add1, *add2, *add3;
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *tmp,*tm1,*tm2,	// Non-static utility ptrs
			/* *two,*one,*sqrt2, */*isrt2, /* *cc1,*ss1,*cc2,*ss2,*cc3,*ss3, *five, */
			*ua0,//*ua1,*ua2,*ua3,*ua4,*ua5,*ua6,*ua7,*ua8,*ua9,
			//*ub0,*ub1,*ub2,*ub3,*ub4,*ub5,*ub6,*ub7,*ub8,*ub9,
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

	  // FMA-based SIMD or (scalar-double) + (LO_ADD = 1 in masterdefs.h) use these sincos constants:
	  #if (defined(USE_AVX2) && DFT_11_FMA) || defined(USE_ARM_V8_SIMD) || (!defined(USE_SSE2) && (LO_ADD != 0))
		const double
			a1 =  0.84125353283118116886,	/* Real part of exp(i*2*pi/11), the radix-11 fundamental sincos datum	*/
			b1 =  0.54064081745559758210,	/* Imag part of exp(i*2*pi/11).	*/
			a2 =  0.41541501300188642553,	/* cos(2u)	*/
			b2 =  0.90963199535451837140,	/* sin(2u)	*/
			a3 = -0.14231483827328514043,	/* cos(3u)	*/
			b3 =  0.98982144188093273237,	/* sin(3u)	*/
			a4 = -0.65486073394528506404,	/* cos(4u)	*/
			b4 =  0.75574957435425828378,	/* sin(4u)	*/
			a5 = -0.95949297361449738988,	/* cos(5u)	*/
			b5 =  0.28173255684142969773;	/* sin(5u)	*/
	  #else
		const double
			a0 =  2.31329240211767848235, /* a0 = (   cq0      -  cq3+  cq2-  cq4)		*/
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
	  #endif	// LO_ADD ?
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
		p1 = NDIVR;		p10 = NDIVR<<4;		p100 = NDIVR<<8;
		p2 = p1 + p1;	p20 = p10 + p10;	p110 = p100 + p10;
		p3 = p2 + p1;	p30 = p20 + p10;	p120 = p110 + p10;
		p4 = p3 + p1;	p40 = p30 + p10;	p130 = p120 + p10;
		p5 = p4 + p1;	p50 = p40 + p10;	p140 = p130 + p10;
		p6 = p5 + p1;	p60 = p50 + p10;	p150 = p140 + p10;
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
													p120 += ( (p120 >> DAT_BITS) << PAD_BITS );
													p130 += ( (p130 >> DAT_BITS) << PAD_BITS );
													p140 += ( (p140 >> DAT_BITS) << PAD_BITS );
													p150 += ( (p150 >> DAT_BITS) << PAD_BITS );
		l = 0;
		dif_phi[l++] =   0;
		dif_phi[l++] = p140;
		dif_phi[l++] = p120;
		dif_phi[l++] = p100;
		dif_phi[l++] = pe0;
		dif_phi[l++] = pc0;
		dif_phi[l++] = pa0;
		dif_phi[l++] = p80;
		dif_phi[l++] = p60;
		dif_phi[l++] = p40;
		dif_phi[l++] = p20;
		l = 0;
		dit_phi[l++] =   0;
		dit_phi[l++] = p20;
		dit_phi[l++] = p40;
		dit_phi[l++] = p60;
		dit_phi[l++] = p80;
		dit_phi[l++] = pa0;
		dit_phi[l++] = pc0;
		dit_phi[l++] = pe0;
		dit_phi[l++] =p100;
		dit_phi[l++] =p120;
		dit_phi[l++] =p140;
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
		poff[0x48+0] =p120; poff[0x48+1] =p120+p4; poff[0x48+2] =p120+p8; poff[0x48+3] =p120+pc;
		poff[0x4c+0] =p130; poff[0x4c+1] =p130+p4; poff[0x4c+2] =p130+p8; poff[0x4c+3] =p130+pc;
		poff[0x50+0] =p140; poff[0x50+1] =p140+p4; poff[0x50+2] =p140+p8; poff[0x50+3] =p140+pc;
		poff[0x54+0] =p150; poff[0x54+1] =p150+p4; poff[0x54+2] =p150+p8; poff[0x54+3] =p150+pc;

	// Shared:
	#ifndef USE_SSE2
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
	#endif

	/*** DIF indexing stuff: ***/

	   #ifdef USE_SSE2
		// Since SIMD code stores DIT-outs into contig-local mem rather than back into large-strided main-array locs,
		// replacing p10*[] with []<<1 gives vec_dbl-complex stride analogs of the p-mults used here in scalar-double mode:

		// Init storage for 2 circular-shifts perms of a basic 11-vector, with shift count in [0,10] that means 2*22 = 42 elts:
		l = 0; j = L2_SZ_VD+1;
		// Even multiples of p10 cshift array:
		dif_p20_cperms[l++] = 0;
		dif_p20_cperms[l++] = 0x140<<j;
		dif_p20_cperms[l++] = 0x120<<j;
		dif_p20_cperms[l++] = 0x100<<j;
		dif_p20_cperms[l++] = 0xe0<<j;
		dif_p20_cperms[l++] = 0xc0<<j;
		dif_p20_cperms[l++] = 0xa0<<j;
		dif_p20_cperms[l++] = 0x80<<j;
		dif_p20_cperms[l++] = 0x60<<j;
		dif_p20_cperms[l++] = 0x40<<j;
		dif_p20_cperms[l++] = 0x20<<j;
		while(l < 2*ODD_RADIX-1) {
			dif_p20_cperms[l] = dif_p20_cperms[l - ODD_RADIX]; ++l;
		}
		// Odd multiples of p10 cshift array:
		dif_p20_cperms[l++] = 0x150<<j;
		dif_p20_cperms[l++] = 0x130<<j;
		dif_p20_cperms[l++] = 0x110<<j;
		dif_p20_cperms[l++] = 0xf0<<j;
		dif_p20_cperms[l++] = 0xd0<<j;
		dif_p20_cperms[l++] = 0xb0<<j;
		dif_p20_cperms[l++] = 0x90<<j;
		dif_p20_cperms[l++] = 0x70<<j;
		dif_p20_cperms[l++] = 0x50<<j;
		dif_p20_cperms[l++] = 0x30<<j;
		dif_p20_cperms[l++] = 0x10<<j;
		while(l < 4*ODD_RADIX-2) {
			dif_p20_cperms[l] = dif_p20_cperms[l - ODD_RADIX]; ++l;
		}
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-11 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-10; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-11 index, bits <4:30> for the p0-f:
		l = 0;
		dif_p20_lo_offset[l++] = ((  0 << 5) + 0 );
		dif_p20_lo_offset[l++] = ((0x5 << 5) + 0 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xa << 5) + 1 );
		dif_p20_lo_offset[l++] = ((0xf << 5) + 1 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x4 << 5) + 1 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x9 << 5) + 2 );
		dif_p20_lo_offset[l++] = ((0xe << 5) + 2 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x3 << 5) + 2 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x8 << 5) + 3 );
		dif_p20_lo_offset[l++] = ((0xd << 5) + 3 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x2 << 5) + 3 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x7 << 5) + 4 );
		dif_p20_lo_offset[l++] = ((0xc << 5) + 4 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x1 << 5) + 4 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x6 << 5) + 5 );
		dif_p20_lo_offset[l++] = ((0xb << 5) + 5 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((  0 << 5) + 5 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0x5 << 5) + 6 );
		dif_p20_lo_offset[l++] = ((0xa << 5) + 6 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xf << 5) + 7 );
		dif_p20_lo_offset[l++] = ((0x4 << 5) + 7 );
		dif_p20_lo_offset[l++] = ((0x9 << 5) + 7 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xe << 5) + 8 );
		dif_p20_lo_offset[l++] = ((0x3 << 5) + 8 );
		dif_p20_lo_offset[l++] = ((0x8 << 5) + 8 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xd << 5) + 9 );
		dif_p20_lo_offset[l++] = ((0x2 << 5) + 9 );
		dif_p20_lo_offset[l++] = ((0x7 << 5) + 9 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xc << 5) + 10);
		dif_p20_lo_offset[l++] = ((0x1 << 5) + 10);
		dif_p20_lo_offset[l++] = ((0x6 << 5) + 10)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((0xb << 5) + 0 );

	   #else

		// Init storage for 2 circular-shifts perms of a basic 11-vector, with shift count in [0,10] that means 2*22 = 42 elts:
		l = 0;
		// Even multiples of p10 cshift array:
		dif_p20_cperms[l++] = 0;
		dif_p20_cperms[l++] = p140;
		dif_p20_cperms[l++] = p120;
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
		dif_p20_cperms[l++] = p150;
		dif_p20_cperms[l++] = p130;
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
		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-11 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-10; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-11 index, bits <4:30> for the p0-f:
		l = 0;
		dif_p20_lo_offset[l++] = (( 0 << 4) + 0 );
		dif_p20_lo_offset[l++] = ((p5 << 4) + 0 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pa << 4) + 1 );
		dif_p20_lo_offset[l++] = ((pf << 4) + 1 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p4 << 4) + 1 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p9 << 4) + 2 );
		dif_p20_lo_offset[l++] = ((pe << 4) + 2 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p3 << 4) + 2 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p8 << 4) + 3 );
		dif_p20_lo_offset[l++] = ((pd << 4) + 3 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p2 << 4) + 3 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p7 << 4) + 4 );
		dif_p20_lo_offset[l++] = ((pc << 4) + 4 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p1 << 4) + 4 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p6 << 4) + 5 );
		dif_p20_lo_offset[l++] = ((pb << 4) + 5 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = (( 0 << 4) + 5 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((p5 << 4) + 6 );
		dif_p20_lo_offset[l++] = ((pa << 4) + 6 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pf << 4) + 7 );
		dif_p20_lo_offset[l++] = ((p4 << 4) + 7 );
		dif_p20_lo_offset[l++] = ((p9 << 4) + 7 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pe << 4) + 8 );
		dif_p20_lo_offset[l++] = ((p3 << 4) + 8 );
		dif_p20_lo_offset[l++] = ((p8 << 4) + 8 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pd << 4) + 9 );
		dif_p20_lo_offset[l++] = ((p2 << 4) + 9 );
		dif_p20_lo_offset[l++] = ((p7 << 4) + 9 )+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pc << 4) + 10);
		dif_p20_lo_offset[l++] = ((p1 << 4) + 10);
		dif_p20_lo_offset[l++] = ((p6 << 4) + 10)+ ((uint32)1 << 31);
		dif_p20_lo_offset[l++] = ((pb << 4) + 0 );

	   #endif	// sse2?

	// dif_offsets are w.r.to a-array, need 11 distinct sets of these, one for each DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		l = 0; ip = dif_offsets;	// In the following, the 'big' += p[>= 20] offsets already in dif_phi[], only care about += p0|p10:
		// Set 0: [0,1,2,3,5,4,7,6,b,a,8,9,e,f,d,c + p00],[7,6,4,5,1,0,3,2,d,c,f,e,a,b,9,8 + p10]
		*ip++= 0;*ip++=p1;*ip++=p2;*ip++=p3;*ip++=p5;*ip++=p4;*ip++=p7;*ip++=p6;*ip++=pb;*ip++=pa;*ip++=p8;*ip++=p9;*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc; for(j = l; j < l+16; j++) { dif_offsets[j] +=   0; }; l += 16; *ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8; for(j = l; j < l+16; j++) { dif_offsets[j] += p10; }; l += 16;
		// Set 1: [8,9,a,b,d,c,f,e,7,6,4,5,1,0,3,2 +p140],[f,e,c,d,9,8,b,a,3,2,0,1,6,7,5,4 +p150]
		*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2; for(j = l; j < l+16; j++) { dif_offsets[j] +=   0; }; l += 16; *ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd;*ip++=p9;*ip++=p8;*ip++=pb;*ip++=pa;*ip++=p3;*ip++=p2;*ip++= 0;*ip++=p1;*ip++=p6;*ip++=p7;*ip++=p5;*ip++=p4; for(j = l; j < l+16; j++) { dif_offsets[j] += p10; }; l += 16;
		// Set 2: [4,5,6,7,3,2,0,1,f,e,c,d,9,8,b,a +p130],[8,9,a,b,d,c,f,e,7,6,4,5,1,0,3,2 +p120]
		*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=p3;*ip++=p2;*ip++= 0;*ip++=p1;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd;*ip++=p9;*ip++=p8;*ip++=pb;*ip++=pa; for(j = l; j < l+16; j++) { dif_offsets[j] += p10; }; l += 16; *ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2; for(j = l; j < l+16; j++) { dif_offsets[j] +=   0; }; l += 16;
		// Set 3: [2,3,1,0,7,6,4,5,8,9,a,b,d,c,f,e +p100],[4,5,6,7,3,2,0,1,f,e,c,d,9,8,b,a +p110]
		*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe; for(j = l; j < l+16; j++) { dif_offsets[j] +=   0; }; l += 16; *ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=p3;*ip++=p2;*ip++= 0;*ip++=p1;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd;*ip++=p9;*ip++=p8;*ip++=pb;*ip++=pa; for(j = l; j < l+16; j++) { dif_offsets[j] += p10; }; l += 16;
		// Set 4: [a,b,9,8,f,e,c,d,4,5,6,7,3,2,0,1 + pf0],[2,3,1,0,7,6,4,5,8,9,a,b,d,c,f,e + pe0]
		*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=p3;*ip++=p2;*ip++= 0;*ip++=p1; for(j = l; j < l+16; j++) { dif_offsets[j] += p10; }; l += 16; *ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe; for(j = l; j < l+16; j++) { dif_offsets[j] +=   0; }; l += 16;
		// Set 5: [e,f,d,c,8,9,a,b,2,3,1,0,7,6,4,5 + pc0],[a,b,9,8,f,e,c,d,4,5,6,7,3,2,0,1 + pd0]
		*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5; for(j = l; j < l+16; j++) { dif_offsets[j] +=   0; }; l += 16; *ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=p3;*ip++=p2;*ip++= 0;*ip++=p1; for(j = l; j < l+16; j++) { dif_offsets[j] += p10; }; l += 16;
		// Set 6: [1,0,3,2,4,5,6,7,a,b,9,8,f,e,c,d + pb0],[e,f,d,c,8,9,a,b,2,3,1,0,7,6,4,5 + pa0]
		*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd; for(j = l; j < l+16; j++) { dif_offsets[j] += p10; }; l += 16; *ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5; for(j = l; j < l+16; j++) { dif_offsets[j] +=   0; }; l += 16;
		// Set 7: [5,4,7,6,2,3,1,0,e,f,d,c,8,9,a,b + p80],[1,0,3,2,4,5,6,7,a,b,9,8,f,e,c,d + p90]
		*ip++=p5;*ip++=p4;*ip++=p7;*ip++=p6;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb; for(j = l; j < l+16; j++) { dif_offsets[j] +=   0; }; l += 16; *ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=pf;*ip++=pe;*ip++=pc;*ip++=pd; for(j = l; j < l+16; j++) { dif_offsets[j] += p10; }; l += 16;
		// Set 8: [d,c,f,e,a,b,9,8,1,0,3,2,4,5,6,7 + p70],[5,4,7,6,2,3,1,0,e,f,d,c,8,9,a,b + p60]
		*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7; for(j = l; j < l+16; j++) { dif_offsets[j] += p10; }; l += 16; *ip++=p5;*ip++=p4;*ip++=p7;*ip++=p6;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0;*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p8;*ip++=p9;*ip++=pa;*ip++=pb; for(j = l; j < l+16; j++) { dif_offsets[j] +=   0; }; l += 16;
		// Set 9: [b,a,8,9,e,f,d,c,5,4,7,6,2,3,1,0 + p40],[d,c,f,e,a,b,9,8,1,0,3,2,4,5,6,7 + p50]
		*ip++=pb;*ip++=pa;*ip++=p8;*ip++=p9;*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p5;*ip++=p4;*ip++=p7;*ip++=p6;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0; for(j = l; j < l+16; j++) { dif_offsets[j] +=   0; }; l += 16; *ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=p4;*ip++=p5;*ip++=p6;*ip++=p7; for(j = l; j < l+16; j++) { dif_offsets[j] += p10; }; l += 16;
		// Set A: [7,6,4,5,1,0,3,2,d,c,f,e,a,b,9,8 + p30],[b,a,8,9,e,f,d,c,5,4,7,6,2,3,1,0 + p20]
		*ip++=p7;*ip++=p6;*ip++=p4;*ip++=p5;*ip++=p1;*ip++= 0;*ip++=p3;*ip++=p2;*ip++=pd;*ip++=pc;*ip++=pf;*ip++=pe;*ip++=pa;*ip++=pb;*ip++=p9;*ip++=p8; for(j = l; j < l+16; j++) { dif_offsets[j] += p10; }; l += 16; *ip++=pb;*ip++=pa;*ip++=p8;*ip++=p9;*ip++=pe;*ip++=pf;*ip++=pd;*ip++=pc;*ip++=p5;*ip++=p4;*ip++=p7;*ip++=p6;*ip++=p2;*ip++=p3;*ip++=p1;*ip++= 0; for(j = l; j < l+16; j++) { dif_offsets[j] +=   0; }; l += 16;
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
		// Init storage for 2 circular-shifts perms of a basic 11-vector, with shift count in [0,10] that means 2*22 = 42 elts:
		l = 0; j = L2_SZ_VD+1;
		// Even multiples of p10 cshift array:
		dit_p20_cperms[l++] = 0;
		dit_p20_cperms[l++] = 0x140<<j;
		dit_p20_cperms[l++] = 0x120<<j;
		dit_p20_cperms[l++] = 0x100<<j;
		dit_p20_cperms[l++] = 0xe0<<j;
		dit_p20_cperms[l++] = 0xc0<<j;
		dit_p20_cperms[l++] = 0xa0<<j;
		dit_p20_cperms[l++] = 0x80<<j;
		dit_p20_cperms[l++] = 0x60<<j;
		dit_p20_cperms[l++] = 0x40<<j;
		dit_p20_cperms[l++] = 0x20<<j;
		while(l < 2*ODD_RADIX-1) {
			dit_p20_cperms[l] = dit_p20_cperms[l - ODD_RADIX]; ++l;
		}
		// Odd multiples of p10 cshift array:
		dit_p20_cperms[l++] = 0x150<<j;
		dit_p20_cperms[l++] = 0x130<<j;
		dit_p20_cperms[l++] = 0x110<<j;
		dit_p20_cperms[l++] = 0xf0<<j;
		dit_p20_cperms[l++] = 0xd0<<j;
		dit_p20_cperms[l++] = 0xb0<<j;
		dit_p20_cperms[l++] = 0x90<<j;
		dit_p20_cperms[l++] = 0x70<<j;
		dit_p20_cperms[l++] = 0x50<<j;
		dit_p20_cperms[l++] = 0x30<<j;
		dit_p20_cperms[l++] = 0x10<<j;
		while(l < 4*ODD_RADIX-2) {
			dit_p20_cperms[l] = dit_p20_cperms[l - ODD_RADIX]; ++l;
		}

		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-11 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-8 index, bits <4:30> for the p0-f:
		l = 0;
		dit_p20_lo_offset[l++] = ((  0 << 5) + 0 );
		dit_p20_lo_offset[l++] = ((0xf << 5) + 1 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xe << 5) + 2 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xd << 5) + 3 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xc << 5) + 4 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xb << 5) + 5 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xa << 5) + 6 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x9 << 5) + 7 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x8 << 5) + 8 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x7 << 5) + 9 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x6 << 5) + 10) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x5 << 5) + 0 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x4 << 5) + 1 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x3 << 5) + 2 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x2 << 5) + 3 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0x1 << 5) + 4 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((  0 << 5) + 5 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((0xf << 5) + 7 );
		dit_p20_lo_offset[l++] = ((0xe << 5) + 8 );
		dit_p20_lo_offset[l++] = ((0xd << 5) + 9 );
		dit_p20_lo_offset[l++] = ((0xc << 5) + 10);
		dit_p20_lo_offset[l++] = ((0xb << 5) + 0 );
		dit_p20_lo_offset[l++] = ((0xa << 5) + 1 );
		dit_p20_lo_offset[l++] = ((0x9 << 5) + 2 );
		dit_p20_lo_offset[l++] = ((0x8 << 5) + 3 );
		dit_p20_lo_offset[l++] = ((0x7 << 5) + 4 );
		dit_p20_lo_offset[l++] = ((0x6 << 5) + 5 );
		dit_p20_lo_offset[l++] = ((0x5 << 5) + 6 );
		dit_p20_lo_offset[l++] = ((0x4 << 5) + 7 );
		dit_p20_lo_offset[l++] = ((0x3 << 5) + 8 );
		dit_p20_lo_offset[l++] = ((0x2 << 5) + 9 );
		dit_p20_lo_offset[l++] = ((0x1 << 5) + 10);

	   #else

		// Init storage for 2 circular-shifts perms of a basic 11-vector, with shift count in [0,10] that means 2*22 = 42 elts:
		l = 0;
		// Even multiples of p10 cshift array:
		dit_p20_cperms[l++] = 0;
		dit_p20_cperms[l++] = p140;
		dit_p20_cperms[l++] = p120;
		dit_p20_cperms[l++] = p100;
		dit_p20_cperms[l++] = pe0;
		dit_p20_cperms[l++] = pc0;
		dit_p20_cperms[l++] = pa0;
		dit_p20_cperms[l++] = p80;
		dit_p20_cperms[l++] = p60;
		dit_p20_cperms[l++] = p40;
		dit_p20_cperms[l++] = p20;
		while(l < 2*ODD_RADIX-1) {
			dit_p20_cperms[l] = dit_p20_cperms[l - ODD_RADIX]; ++l;
		}
		// Odd multiples of p10 cshift array:
		dit_p20_cperms[l++] = p150;
		dit_p20_cperms[l++] = p130;
		dit_p20_cperms[l++] = p110;
		dit_p20_cperms[l++] = pf0;
		dit_p20_cperms[l++] = pd0;
		dit_p20_cperms[l++] = pb0;
		dit_p20_cperms[l++] = p90;
		dit_p20_cperms[l++] = p70;
		dit_p20_cperms[l++] = p50;
		dit_p20_cperms[l++] = p30;
		dit_p20_cperms[l++] = p10;
		while(l < 4*ODD_RADIX-2) {
			dit_p20_cperms[l] = dit_p20_cperms[l - ODD_RADIX]; ++l;
		}

		// Low parts, i.e. (mod p20) of the p-index offsets in the above circ-perm-indexing scheme for the radix-11 DFTs.
		// Each elt of form p[0-f] + [evn|odd]0-4; use high-bit-toggle to encode the [evn|odd] selector, low 4 bits for
		// the 0-8 index, bits <4:30> for the p0-f:
		l = 0;
		dit_p20_lo_offset[l++] = (( 0 << 4) + 0 );
		dit_p20_lo_offset[l++] = ((pf << 4) + 1 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pe << 4) + 2 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pd << 4) + 3 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pc << 4) + 4 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pb << 4) + 5 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pa << 4) + 6 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p9 << 4) + 7 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p8 << 4) + 8 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p7 << 4) + 9 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p6 << 4) + 10) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p5 << 4) + 0 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p4 << 4) + 1 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p3 << 4) + 2 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p2 << 4) + 3 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((p1 << 4) + 4 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = (( 0 << 4) + 5 ) + ((uint32)1 << 31);
		dit_p20_lo_offset[l++] = ((pf << 4) + 7 );
		dit_p20_lo_offset[l++] = ((pe << 4) + 8 );
		dit_p20_lo_offset[l++] = ((pd << 4) + 9 );
		dit_p20_lo_offset[l++] = ((pc << 4) + 10);
		dit_p20_lo_offset[l++] = ((pb << 4) + 0 );
		dit_p20_lo_offset[l++] = ((pa << 4) + 1 );
		dit_p20_lo_offset[l++] = ((p9 << 4) + 2 );
		dit_p20_lo_offset[l++] = ((p8 << 4) + 3 );
		dit_p20_lo_offset[l++] = ((p7 << 4) + 4 );
		dit_p20_lo_offset[l++] = ((p6 << 4) + 5 );
		dit_p20_lo_offset[l++] = ((p5 << 4) + 6 );
		dit_p20_lo_offset[l++] = ((p4 << 4) + 7 );
		dit_p20_lo_offset[l++] = ((p3 << 4) + 8 );
		dit_p20_lo_offset[l++] = ((p2 << 4) + 9 );
		dit_p20_lo_offset[l++] = ((p1 << 4) + 10);

	   #endif	// sse2?

	// dit_offsets are w.r.to a-array, need 11 distinct sets of these, one for each radix-32 DFT.
	// NB: We could trivially include the full p10 multiples here, but prefer to do things (mod p20)
	// firstly for aesthetic reasons - all array elements in [0,p20) - and secondly to provide
	// for the possibility of a streamlined smaller-sub-array-based encoding at a later date.
		l = 0; ip = dit_offsets;	// In the following, the 'big' += p[>= 20] offsets already in dif_phi[], only care about += p0|p10:
		// Set 0: [0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p00],[f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p10]
		*ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2;*ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4;*ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc;*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8; for(j = l; j < l+16; j++) { dit_offsets[j] +=   0; }; l += 16; *ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc;*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8;*ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4;*ip++=p3;*ip++=p2;*ip++=p1;*ip++= 0; for(j = l; j < l+16; j++) { dit_offsets[j] += p10; }; l += 16;
		// Set 1: [7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + p30],[7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + p20]
		*ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4;*ip++=p3;*ip++=p2;*ip++=p1;*ip++= 0;*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8;*ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf; for(j = l; j < l+16; j++) { dit_offsets[j] += p10; }; l += 16; *ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4;*ip++=p3;*ip++=p2;*ip++=p1;*ip++= 0;*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8;*ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf; for(j = l; j < l+16; j++) { dit_offsets[j] +=   0; }; l += 16;
		// Set 2: [b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p40],[3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p50]
		*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8;*ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf;*ip++=p3;*ip++=p2;*ip++=p1;*ip++= 0;*ip++=p5;*ip++=p4;*ip++=p6;*ip++=p7; for(j = l; j < l+16; j++) { dit_offsets[j] +=   0; }; l += 16; *ip++=p3;*ip++=p2;*ip++=p1;*ip++= 0;*ip++=p5;*ip++=p4;*ip++=p6;*ip++=p7;*ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf;*ip++=p9;*ip++=p8;*ip++=pa;*ip++=pb; for(j = l; j < l+16; j++) { dit_offsets[j] += p10; }; l += 16;
		// Set 3: [d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + p70],[d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + p60]
		*ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf;*ip++=p9;*ip++=p8;*ip++=pa;*ip++=pb;*ip++=p5;*ip++=p4;*ip++=p6;*ip++=p7;*ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3; for(j = l; j < l+16; j++) { dit_offsets[j] += p10; }; l += 16; *ip++=pd;*ip++=pc;*ip++=pe;*ip++=pf;*ip++=p9;*ip++=p8;*ip++=pa;*ip++=pb;*ip++=p5;*ip++=p4;*ip++=p6;*ip++=p7;*ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3; for(j = l; j < l+16; j++) { dit_offsets[j] +=   0; }; l += 16;
		// Set 4: [5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p80],[9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p90]
		*ip++=p5;*ip++=p4;*ip++=p6;*ip++=p7;*ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3;*ip++=p9;*ip++=p8;*ip++=pa;*ip++=pb;*ip++=pe;*ip++=pf;*ip++=pc;*ip++=pd; for(j = l; j < l+16; j++) { dit_offsets[j] +=   0; }; l += 16; *ip++=p9;*ip++=p8;*ip++=pa;*ip++=pb;*ip++=pe;*ip++=pf;*ip++=pc;*ip++=pd;*ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3;*ip++=p6;*ip++=p7;*ip++=p4;*ip++=p5; for(j = l; j < l+16; j++) { dit_offsets[j] += p10; }; l += 16;
		// Set 5: [1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pb0],[1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + pa0]
		*ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3;*ip++=p6;*ip++=p7;*ip++=p4;*ip++=p5;*ip++=pe;*ip++=pf;*ip++=pc;*ip++=pd;*ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9; for(j = l; j < l+16; j++) { dit_offsets[j] += p10; }; l += 16; *ip++=p1;*ip++= 0;*ip++=p2;*ip++=p3;*ip++=p6;*ip++=p7;*ip++=p4;*ip++=p5;*ip++=pe;*ip++=pf;*ip++=pc;*ip++=pd;*ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9; for(j = l; j < l+16; j++) { dit_offsets[j] +=   0; }; l += 16;
		// Set 6: [e,f,c,d,a,b,8,9,6,7,4,5,2,3,0,1 + pc0],[6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + pd0]
		*ip++=pe;*ip++=pf;*ip++=pc;*ip++=pd;*ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9;*ip++=p6;*ip++=p7;*ip++=p4;*ip++=p5;*ip++=p2;*ip++=p3;*ip++= 0;*ip++=p1; for(j = l; j < l+16; j++) { dit_offsets[j] +=   0; }; l += 16; *ip++=p6;*ip++=p7;*ip++=p4;*ip++=p5;*ip++=p2;*ip++=p3;*ip++= 0;*ip++=p1;*ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9;*ip++=pc;*ip++=pd;*ip++=pf;*ip++=pe; for(j = l; j < l+16; j++) { dit_offsets[j] += p10; }; l += 16;
		// Set 7: [a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + pf0],[a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + pe0]
		*ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9;*ip++=pc;*ip++=pd;*ip++=pf;*ip++=pe;*ip++=p2;*ip++=p3;*ip++= 0;*ip++=p1;*ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6; for(j = l; j < l+16; j++) { dit_offsets[j] += p10; }; l += 16; *ip++=pa;*ip++=pb;*ip++=p8;*ip++=p9;*ip++=pc;*ip++=pd;*ip++=pf;*ip++=pe;*ip++=p2;*ip++=p3;*ip++= 0;*ip++=p1;*ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6; for(j = l; j < l+16; j++) { dit_offsets[j] +=   0; }; l += 16;
		// Set 8: [2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a +p100],[c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 +p110]
		*ip++=p2;*ip++=p3;*ip++= 0;*ip++=p1;*ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6;*ip++=pc;*ip++=pd;*ip++=pf;*ip++=pe;*ip++=p8;*ip++=p9;*ip++=pb;*ip++=pa; for(j = l; j < l+16; j++) { dit_offsets[j] +=   0; }; l += 16; *ip++=pc;*ip++=pd;*ip++=pf;*ip++=pe;*ip++=p8;*ip++=p9;*ip++=pb;*ip++=pa;*ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6;*ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2; for(j = l; j < l+16; j++) { dit_offsets[j] += p10; }; l += 16;
		// Set 9: [4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c +p130],[4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c +p120]
		*ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6;*ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2;*ip++=p8;*ip++=p9;*ip++=pb;*ip++=pa;*ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc; for(j = l; j < l+16; j++) { dit_offsets[j] += p10; }; l += 16; *ip++=p4;*ip++=p5;*ip++=p7;*ip++=p6;*ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2;*ip++=p8;*ip++=p9;*ip++=pb;*ip++=pa;*ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc; for(j = l; j < l+16; j++) { dit_offsets[j] +=   0; }; l += 16;
		// Set A: [8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 +p140],[0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 +p150]
		*ip++=p8;*ip++=p9;*ip++=pb;*ip++=pa;*ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc;*ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2;*ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4; for(j = l; j < l+16; j++) { dit_offsets[j] +=   0; }; l += 16; *ip++= 0;*ip++=p1;*ip++=p3;*ip++=p2;*ip++=p7;*ip++=p6;*ip++=p5;*ip++=p4;*ip++=pf;*ip++=pe;*ip++=pd;*ip++=pc;*ip++=pb;*ip++=pa;*ip++=p9;*ip++=p8; for(j = l; j < l+16; j++) { dit_offsets[j] += p10; }; l += 16;
	  #ifdef USE_SSE2
		// IN SIMD mode preshift all the above offsets << 3 to turn into double=array pointer offsets:
		for(l = 0; l < RADIX; l++) {
			dit_offsets[l] <<= 3;
		}
	  #endif

	#ifdef USE_SSE2
		tmp	= r00 = thread_arg->r00;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x2c0;	s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
		tmp += 0x2c0;
		//two     = tmp;	// AVX+ versions of various DFT macros assume consts 2.0,1.0,isrt2 laid out thusly
		//one     = tmp + 0x1;
		//sqrt2	= tmp + 0x2;
		isrt2   = tmp + 0x3;
		//cc2	    = tmp + 0x4;	// Radix-32 DFT macros assume roots stored in this [8th, 16th, 32nd_1,3] order
		//ss2	    = tmp + 0x5;
		//cc1	    = tmp + 0x6;
		//ss1	    = tmp + 0x7;
		//cc3	    = tmp + 0x8;
		//ss3	    = tmp + 0x9;
	//	two     = tmp + 0xa; Unnamed ptr, used in AVX2 mode to hold 2.0 (and *five holds 1.0) for the radix-11 DFT code
		//five    = tmp + 0xb;
		ua0     = tmp + 0xc;
		//ua1     = tmp + 0xd;
		//ua2     = tmp + 0xe;
		//ua3     = tmp + 0xf;
		//ua4     = tmp + 0x10;
		//ua5     = tmp + 0x11;
		//ua6     = tmp + 0x12;
		//ua7     = tmp + 0x13;
		//ua8     = tmp + 0x14;
		//ua9     = tmp + 0x15;
		//ub0     = tmp + 0x16;
		//ub1     = tmp + 0x17;
		//ub2     = tmp + 0x18;
		//ub3     = tmp + 0x19;
		//ub4     = tmp + 0x1a;
		//ub5     = tmp + 0x1b;
		//ub6     = tmp + 0x1c;
		//ub7     = tmp + 0x1d;
		//ub8     = tmp + 0x1e;
		//ub9     = tmp + 0x1f;
		tmp += 0x20;	// sc_ptr += 0x5a0
	  #ifdef USE_AVX512
		cy = tmp;		tmp += 0x2c;	// RADIX/8 vec_dbl slots for carry sub-array
		max_err = tmp + 0x00;
		//sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(5a0 + 2c + 2) = 0x5ce; This is where the value of half_arr_offset352 comes from
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
		cy = tmp;		tmp += 0x58;	// RADIX/4 vec_dbl slots
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(5a0 + 58 + 2) = 0x5fa; This is where the value of half_arr_offset352 comes from
		half_arr= tmp + 0x02;	// This table needs 20 vec_dbl in both avx and sse2 mode
	  #else
		cy = tmp;		tmp += 0xb0;	// RADIX/2 vec_dbl slots
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x(5a0 + b0 + 2) = 0x652; This is where the value of half_arr_offset352 comes from
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

		sign_mask = (uint64*)(r00 + radix352_creals_in_local_store);
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
		#include "radix352_main_carry_loop.h"

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
