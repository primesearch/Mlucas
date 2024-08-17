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

#define RADIX 60	// Use #define rather than const int to ensure it's really a compile-time const in the C sense
#define ODD_RADIX 15	// ODD_RADIX = [radix >> trailz(radix)]

#define EPS 1e-10

#ifndef COMPACT_OBJ	// Toggle for parametrized-loop-DFT compact-object code scheme
  #ifdef USE_AVX
	#define COMPACT_OBJ	1
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
  // Max = (40 [SSE2]; 132 [AVX]),
  // For Fermat-mod in AVX mode we need RADIX*4 = 240 [if HIACC] or 12 [if not] slots for the compact
  // negacyclic-roots chained-multiply scheme.
  // Add larger number in each case - i.e. max(132,240) = 240 if AVX, 40 if SSE2
  // to (half_arr_offset60 + RADIX) to get required value of radix60_creals_in_local_store:
  #ifdef USE_AVX512	// 2*floor(RADIX/8) = 14 fewer carry slots than AVX:
	const int half_arr_offset60 = 326;	// + RADIX = 386;  Used for thread local-storage-integrity checking
	const int radix60_creals_in_local_store = 632;	// AVX: 386 + 240 and round up to nearest multiple of 8
  #elif defined(USE_AVX)
	const int half_arr_offset60 = 340;	// + RADIX = 400;  Used for thread local-storage-integrity checking
	const int radix60_creals_in_local_store = 640;	// AVX: 400 + 240
  #else
	const int half_arr_offset60 = 370;	// + RADIX = 430; Used for thread local-storage-integrity checking
	const int radix60_creals_in_local_store = 472;	// (half_arr_offset60 + RADIX) + 40 and round up to nearest multiple of 4
  #endif
  #ifdef USE_AVX
	const uint64 radix60_avx_negadwt_consts[RADIX] = {	// 8 entries per line ==> RADIX/8 lines:
		0x3FF0000000000000ull,0x3FEFFD315BBF4275ull,0x3FEFF4C5ED12E61Dull,0x3FEFE6BF2E2660AFull,0x3FEFD31F94F867C6ull,0x3FEFB9EA92EC689Bull,0x3FEF9B24942FE45Cull,0x3FEF76D2FEF3CC4Bull,
		0x3FEF4CFC327A0080ull,0x3FEF1DA785F71BCEull,0x3FEEE8DD4748BF15ull,0x3FEEAEA6B98095C0ull,0x3FEE6F0E134454FFull,0x3FEE2A1E7D02FE9Full,0x3FEDDFE40EFFB805ull,0x3FED906BCF328D46ull,
		0x3FED3BC3AEFF7F95ull,0x3FECE1FA88C445BBull,0x3FEC83201D3D2C6Dull,0x3FEC1F4510C18B95ull,0x3FEBB67AE8584CAAull,0x3FEB48D406A50540ull,0x3FEAD663A8AE2FDCull,0x3FEA5F3DE27D13F2ull,
		0x3FE9E3779B97F4A8ull,0x3FE963268B572492ull,0x3FE8DE613515A328ull,0x3FE8553EE43DEF13ull,0x3FE7C7D7A833BEC2ull,0x3FE73644501B56CDull,0x3FE6A09E667F3BCDull,0x3FE607002CD5031Dull,
		0x3FE5698496E20BD8ull,0x3FE4C8474600EEEEull,0x3FE4236484487ABEull,0x3FE37AF93F9513EAull,0x3FE2CF2304755A5Eull,0x3FE21FFFF8FAF674ull,0x3FE16DAED770771Dull,0x3FE0B84EE8F52E9Dull,
		0x3FE0000000000000ull,0x3FDE89C4E59427B1ull,0x3FDD0E2E2B44DE01ull,0x3FDB8D7E6A56D476ull,0x3FDA07F921061AD1ull,0x3FD87DE2A6AEA963ull,0x3FD6EF801FCED33Cull,0x3FD55D1771E5BAB9ull,
		0x3FD3C6EF372FE950ull,0x3FD22D4EB2443163ull,0x3FD0907DC1930690ull,0x3FCDE189A594FBCCull,0x3FCA9CD9AC4258F6ull,0x3FC7537E63143E2Eull,0x3FC4060B67A85375ull,0x3FC0B5150F6DA2D1ull,
		0x3FBAC2609B3C576Cull,0x3FB415E532398E49ull,0x3FAACBC748EFC90Eull,0x3F9ACE214390CA91ull
	};
  #endif

	#include "sse2_macro_gcc64.h"
	#include "radix15_sse_macro.h"

#endif	/* USE_SSE2 */

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

int radix60_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-60 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-60 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix60_ditN_cy_dif1";
  #if !defined(MULTITHREAD) && defined(USE_SSE2)
	const int pfetch_dist = PFETCH_DIST;
  #endif
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
  #if !defined(MULTITHREAD) && (!defined(USE_SSE2) || (COMPACT_OBJ != 0))
	static uint32 p0123[4];
  #endif
	static double wts_mult[2], inv_mult[2];	// Const wts-multiplier and 2*(its multiplicative inverse)
  #ifndef MULTITHREAD
	#if !defined(USE_SSE2) || defined(USE_AVX)
	double wt_re,wt_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	#endif
	#ifndef USE_SSE2
	double wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	#endif
  #endif
	// Cleanup loop assumes carryins propagate at most 4 words up, but need at least 1 vec_cmplx
	// (2 vec_dbl)'s worth of doubles in wraparound step, hence AVX-512 needs mers-value bumped up:
  #ifdef USE_AVX512
	const int jhi_wrap_mers = 15;
	const int jhi_wrap_ferm = 15;
  #else
	const int jhi_wrap_mers =  7;
	const int jhi_wrap_ferm = 15;	// For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
  #endif
	int NDIVR,i,j,j1,jt,jstart,jhi,full_pass,khi,l,outer;
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	int j2,jp;
  #endif
  #ifndef MULTITHREAD
	int ntmp;
  #endif
  #ifdef USE_SSE2
	int nbytes;
  #endif
#if !defined(MULTITHREAD) && defined(USE_SSE2)
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = --|--|15 for avx512,avx,sse, respectively.
	// For the present case we structure the carry-macro calls such that this divisibility requirement is only for 128-bit SIMD builds;
	// note fixed-incr too restrictive here, so in sse2 case 'divide 15 into pieces' via increment-array whose elts sum to 15:
  #ifndef USE_AVX
	const int *incr;
  #endif
	const int *inc_arr = 0x0;	// Fermat-mod only uses inc_arr in AVX+ mode, init = 0 to quiet 'uninit' warnings.
	const int incr_long[] = {15}, incr_med[] = {8,7}, incr_short[] = {5,5,5};
  // Have no specialized HIACC carry macro in USE_AVX512 and ARMv8 SIMD, use nonzero incr-value to differenetiate vs AVX/AVX2:
  #if defined(USE_AVX512) || defined(USE_ARM_V8_SIMD)
	const int incr_hiacc[] = {3,3,3,3,3};	// Don't actually use the elements in AVX-512 mode; values are for ARMv8
  #else
	const int incr_hiacc[] = {0};
  #endif
	// Fermat-mod: For AVX+, define carry-subchain length in terms of 2^nfold subchains:
  #ifdef USE_AVX512
	/*** No avx-512 support for radix-60 Fermat-mod, so just declare these with 0-values ***/
	int nfold = 0;
	const int nexec_long[] = {0}, nexec_med[] = {0}, nexec_short[] = {0}, nexec_hiacc[] = {0};
  #elif defined(USE_AVX)
	// For nfold > 4,  RADIX/4 not divisible by 2^nfold, so use a more-general inner-loop scheme which can handle that:
	int nfold = USE_SHORT_CY_CHAIN + 0;
	const int nexec_long[] = {15}, nexec_med[] = {8,7}, nexec_short[] = {5,5,5}, nexec_hiacc[] = {3,3,3,3,3};
  #endif
	// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
		if(USE_SHORT_CY_CHAIN == 0)
			inc_arr = incr_long;
		else if(USE_SHORT_CY_CHAIN == 1)
			inc_arr = incr_med;
		else if(USE_SHORT_CY_CHAIN == 2)
			inc_arr = incr_short;
		else
			inc_arr = incr_hiacc;
	} else {	// MODULUS_TYPE_FERMAT:
	#ifdef USE_AVX
		if(USE_SHORT_CY_CHAIN == 0)
			inc_arr = nexec_long;
		else if(USE_SHORT_CY_CHAIN == 1)
			inc_arr = nexec_med;
		else if(USE_SHORT_CY_CHAIN == 2)
			inc_arr = nexec_short;
		else
			inc_arr = nexec_hiacc;
	#endif
	}
#endif // MULTITHREAD && USE_SSE2

  #ifdef USE_SSE2
	uint32 nwt16 = nwt << L2_SZ_VD;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code
  #endif
  #ifndef MULTITHREAD
	uint32 k1,k2;
   #if defined(USE_SSE2) && !defined(USE_AVX512)
	uint32 k3,k4;
   #endif
  #endif
	/* Need these both in scalar mode and to ease the SSE2-array init...dimension = ODD_RADIX;
	In order to ease the ptr-access for the || routine, lump these 4*ODD_RADIX doubles together with copies of
	the 4 in the passed-in bs[2] and bsinv[2] arrays [and used in this 4-double form by the mersenne-mod carry macros]
	into a single foo_array[4*(ODD_RADIX+1)], then convert what used to be disparate ODD_RADIX-sized arrays to pointers.

	NOTE: Use literal value of ODD_RADIX in dimensioning in order to allow us to explicitly make static ...
	some compilers [e.g. gcc] are OK with foo_array[(ODD_RADIX+1)<<2], others [e.g. icc] will complain like
	         "error: a variable length array cannot have static storage duration"
	even though ODD_RADIX has been declared a const: */
	static double foo_array[(ODD_RADIX+1)<<2], *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr;
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	static double bs,bsinv;
  #endif

	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40,p44,p48,p52,p56, nsave = 0;
	static int poff[RADIX>>2];
	static double radix_inv, n2inv;
  #if !defined(MULTITHREAD) || defined(USE_SSE2)
	const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/
					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
  #endif
	double scale,dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
  #ifndef MULTITHREAD
	double *addr;
   #ifndef USE_SSE2
	double *addi;
   #endif
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
	static int n_div_nwt;
  #ifndef MULTITHREAD
	int col,co2,co3;
   #ifndef USE_SSE2
	int m,m2;
   #endif
   #ifdef USE_AVX512
	double t0,t1,t2,t3;
	static struct uint32x8 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
   #elif defined(USE_AVX)
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
   #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
	#if !defined(USE_SSE2) || defined(USE_AVX)
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	#endif
   #endif
   #if !defined(USE_SSE2) || defined(USE_AVX)
	double rt,it;
   #endif
  #endif
	static int ii[ODD_RADIX] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};	/* indices into weights arrays (mod NWT) */
	/* These are used in conjunction with the langth-ODD_RADIX arrays in the USE_SCALAR_CARRY code flow;
	In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
  #if !defined(MULTITHREAD) && defined(USE_SSE2) && !defined(USE_AVX)
	static int idx_offset, idx_incr;
  #endif
	static int wts_idx_incr = 0
	  #ifdef USE_SSE2
		,wts_idx_inc2 = 0
	  #endif
		,icycle[ODD_RADIX];
  #ifndef MULTITHREAD
	static int ic_idx;
  #endif
#ifdef USE_SSE2
   #if !defined(MULTITHREAD) && defined(USE_AVX)
	int i0,i1,i2,i3;
   #endif
	static int jcycle[ODD_RADIX];
   #ifndef MULTITHREAD
	static int jc_idx;
   #endif
  #ifdef USE_AVX
	static int kcycle[ODD_RADIX];
	static int lcycle[ODD_RADIX];
   #ifndef MULTITHREAD
	static int kc_idx, lc_idx;
   #endif
  #endif
  #ifdef USE_AVX512
	static int mcycle[ODD_RADIX];
	static int ncycle[ODD_RADIX];
	static int ocycle[ODD_RADIX];
	static int pcycle[ODD_RADIX];
   #ifndef MULTITHREAD
	static int mc_idx, nc_idx, oc_idx, pc_idx;
   #endif
  #endif

	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else
	double *add0, *add1, *add2, *add3;
  #endif	// MULTITHREAD

  #ifndef MULTITHREAD
	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
  #endif
  #if !defined(MULTIHREAD) && !defined(USE_AVX512)
	const double crnd = 3.0*0x4000000*0x2000000;
  #endif
	static vec_dbl *max_err, *sse2_rnd, *half_arr, *two, *sse2_c3m1, *sse2_s, *sse2_cn1, *sse2_cn2, *sse2_ss3, *sse2_sn1, *sse2_sn2,
	  #ifndef MULTITHREAD
		*s1p00,
	   #if !COMPACT_OBJ
		*s1p01,*s1p02,*s1p03,*s1p04,*s1p05,*s1p06,*s1p07,*s1p08,*s1p09,*s1p0a,*s1p0b,*s1p0c,*s1p0d,*s1p0e,
		*s1p10,*s1p11,*s1p12,*s1p13,*s1p14,*s1p15,*s1p16,*s1p17,*s1p18,*s1p19,*s1p1a,*s1p1b,*s1p1c,*s1p1d,*s1p1e,
		*s1p20,*s1p21,*s1p22,*s1p23,*s1p24,*s1p25,*s1p26,*s1p27,*s1p28,*s1p29,*s1p2a,*s1p2b,*s1p2c,*s1p2d,*s1p2e,
		*s1p30,*s1p31,*s1p32,*s1p33,*s1p34,*s1p35,*s1p36,*s1p37,*s1p38,*s1p39,*s1p3a,*s1p3b,*s1p3c,*s1p3d,*s1p3e,
	   #endif
	  #endif
		*r00,
	  #ifndef MULTITHREAD
	   #if !COMPACT_OBJ
		*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09,*r0a,*r0b,*r0c,*r0d,*r0e,
		*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r1a,*r1b,*r1c,*r1d,*r1e,
		*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r2a,*r2b,*r2c,*r2d,*r2e,
		*r30,*r31,*r32,*r33,*r34,*r35,*r36,*r37,*r38,*r39,*r3a,*r3b,*r3c,*r3d,*r3e,
	   #endif
		*x00,*x01,*x02,*x03,*x04,*x05,*x06,*x07,*x08,*x09,*x0a,*x0b,*x0c,*x0d,*x0e,
		*y00,*y01,*y02,*y03,*y04,*y05,*y06,*y07,*y08,*y09,*y0a,*y0b,*y0c,*y0d,*y0e,
	  #endif
		*cy_r	// Need RADIX slots for sse2 carries, RADIX/2 for avx
	  #if !defined(MULTITHREAD) && defined(USE_AVX)
		,*cy_i	// Need RADIX slots for sse2 carries, RADIX/2 for avx
	  #endif
		;
  #ifdef USE_AVX
	static vec_dbl *base_negacyclic_root;
  #endif

	vec_dbl *tmp,*tm2;	// Non-static utility ptrs
  #if !defined(MULTITHREAD) && defined(USE_AVX)
	vec_dbl *tm0;	// Non-static utility ptrs
  #endif
  #ifndef MULTITHREAD
	vec_dbl *tm1;	// Non-static utility ptrs
  #endif
  #ifndef USE_AVX
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
  #endif

#endif	// USE_SSE2

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
	static task_control_t   task_control = {NULL, (void*)cy60_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
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

  #if defined(USE_AVX512) && !defined(USE_THREADS)
	WARN(HERE, "radix60_ditN_cy_dif1: AVX-512 support requires multithreaded build (-DUSE_THREADS); Skipping this leading radix.", "", 1); return(ERR_RADIX0_UNAVAILABLE);
  #endif
  #ifdef USE_IMCI512
	WARN(HERE, "radix60_ditN_cy_dif1: No k10m/IMCI-512 support; Skipping this leading radix.", "", 1); return(ERR_RADIX0_UNAVAILABLE);
  #endif
  #ifdef USE_AVX512
	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		WARN(HERE, "radix60_ditN_cy_dif1: No AVX-512 support for Fermat-mod; Skipping this leading radix.", "", 1); return(ERR_RADIX0_UNAVAILABLE);
	}
  #endif

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

  #ifndef MULTITHREAD
	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=-1;
  #endif
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
	//	printf("0: wt*inv-1 = %15.8e\n",fabs(wts_mult[0]*inv_mult[0] - 1.));
		ASSERT(fabs(wts_mult[0]*inv_mult[0] - 1.0) < EPS, "wts_mults fail accuracy check!");
		//curr have w, 2/w, separate-mul-by-1-or-0.5 gives [w,w/2] and [1/w,2/w] for i = 0,1, resp:
		wts_mult[1] = 0.5*wts_mult[0];
		inv_mult[1] = 2.0*inv_mult[0];
		ASSERT(fabs(wts_mult[1]*inv_mult[1] - 1.0) < EPS, "wts_mults fail accuracy check!");
	//	printf("1: wt*inv-1 = %15.8e\n",fabs(wts_mult[1]*inv_mult[1] - 1.));

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
		// consisting of 128*2 vec_dbl and (8+RADIX/2) uint64 element slots per thread, the latter of which
		// provide thread-local storage for int-data and tables
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix60_creals_in_local_store + (((12+RADIX/2)/2 + ODD_RADIX + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix60_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 2*RADIX vector-double-sized slots of sc_arr for s1p* temporaries, next 2*RADIX slots for r* temps,
	next RADIX slots for x and y in-place DFT temps, next 7 for the complex root combs needed for the radix-3 and -5 sub-DFTs,
	next 30[avx] | 60[sse2] for the doubled carry pairs, next 2 for ROE and RND_CONST, next 5*RADIX[avx] | 20[sse2] for the
	half_arr table lookup stuff, plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		r00 = sc_ptr;			tmp = r00 + 0x78;
		r00 = r00 + 0x00;
	  #ifndef MULTITHREAD
								s1p00 = tmp + 0x00;
	   #if !COMPACT_OBJ
		r01 = r00 + 0x02;		s1p01 = tmp + 0x02;
		r02 = r00 + 0x04;		s1p02 = tmp + 0x04;
		r03 = r00 + 0x06;		s1p03 = tmp + 0x06;
		r04 = r00 + 0x08;		s1p04 = tmp + 0x08;
		r05 = r00 + 0x0a;		s1p05 = tmp + 0x0a;
		r06 = r00 + 0x0c;		s1p06 = tmp + 0x0c;
		r07 = r00 + 0x0e;		s1p07 = tmp + 0x0e;
		r08 = r00 + 0x10;		s1p08 = tmp + 0x10;
		r09 = r00 + 0x12;		s1p09 = tmp + 0x12;
		r0a = r00 + 0x14;		s1p0a = tmp + 0x14;
		r0b = r00 + 0x16;		s1p0b = tmp + 0x16;
		r0c = r00 + 0x18;		s1p0c = tmp + 0x18;
		r0d = r00 + 0x1a;		s1p0d = tmp + 0x1a;
		r0e = r00 + 0x1c;		s1p0e = tmp + 0x1c;
		r10 = r00 + 0x1e;		s1p10 = tmp + 0x1e;
		r11 = r00 + 0x20;		s1p11 = tmp + 0x20;
		r12 = r00 + 0x22;		s1p12 = tmp + 0x22;
		r13 = r00 + 0x24;		s1p13 = tmp + 0x24;
		r14 = r00 + 0x26;		s1p14 = tmp + 0x26;
		r15 = r00 + 0x28;		s1p15 = tmp + 0x28;
		r16 = r00 + 0x2a;		s1p16 = tmp + 0x2a;
		r17 = r00 + 0x2c;		s1p17 = tmp + 0x2c;
		r18 = r00 + 0x2e;		s1p18 = tmp + 0x2e;
		r19 = r00 + 0x30;		s1p19 = tmp + 0x30;
		r1a = r00 + 0x32;		s1p1a = tmp + 0x32;
		r1b = r00 + 0x34;		s1p1b = tmp + 0x34;
		r1c = r00 + 0x36;		s1p1c = tmp + 0x36;
		r1d = r00 + 0x38;		s1p1d = tmp + 0x38;
		r1e = r00 + 0x3a;		s1p1e = tmp + 0x3a;
		r20 = r00 + 0x3c;		s1p20 = tmp + 0x3c;
		r21 = r00 + 0x3e;		s1p21 = tmp + 0x3e;
		r22 = r00 + 0x40;		s1p22 = tmp + 0x40;
		r23 = r00 + 0x42;		s1p23 = tmp + 0x42;
		r24 = r00 + 0x44;		s1p24 = tmp + 0x44;
		r25 = r00 + 0x46;		s1p25 = tmp + 0x46;
		r26 = r00 + 0x48;		s1p26 = tmp + 0x48;
		r27 = r00 + 0x4a;		s1p27 = tmp + 0x4a;
		r28 = r00 + 0x4c;		s1p28 = tmp + 0x4c;
		r29 = r00 + 0x4e;		s1p29 = tmp + 0x4e;
		r2a = r00 + 0x50;		s1p2a = tmp + 0x50;
		r2b = r00 + 0x52;		s1p2b = tmp + 0x52;
		r2c = r00 + 0x54;		s1p2c = tmp + 0x54;
		r2d = r00 + 0x56;		s1p2d = tmp + 0x56;
		r2e = r00 + 0x58;		s1p2e = tmp + 0x58;
		r30 = r00 + 0x5a;		s1p30 = tmp + 0x5a;
		r31 = r00 + 0x5c;		s1p31 = tmp + 0x5c;
		r32 = r00 + 0x5e;		s1p32 = tmp + 0x5e;
		r33 = r00 + 0x60;		s1p33 = tmp + 0x60;
		r34 = r00 + 0x62;		s1p34 = tmp + 0x62;
		r35 = r00 + 0x64;		s1p35 = tmp + 0x64;
		r36 = r00 + 0x66;		s1p36 = tmp + 0x66;
		r37 = r00 + 0x68;		s1p37 = tmp + 0x68;
		r38 = r00 + 0x6a;		s1p38 = tmp + 0x6a;
		r39 = r00 + 0x6c;		s1p39 = tmp + 0x6c;
		r3a = r00 + 0x6e;		s1p3a = tmp + 0x6e;
		r3b = r00 + 0x70;		s1p3b = tmp + 0x70;
		r3c = r00 + 0x72;		s1p3c = tmp + 0x72;
		r3d = r00 + 0x74;		s1p3d = tmp + 0x74;
		r3e = r00 + 0x76;		s1p3e = tmp + 0x76;
	   #endif
	  #endif
		tmp += 0x78;	// sc_ptr + 240
	  #ifndef MULTITHREAD
		x00   = tmp + 0x00;
		x01   = tmp + 0x02;
		x02   = tmp + 0x04;
		x03   = tmp + 0x06;
		x04   = tmp + 0x08;
		x05   = tmp + 0x0a;
		x06   = tmp + 0x0c;
		x07   = tmp + 0x0e;
		x08   = tmp + 0x10;
		x09   = tmp + 0x12;
		x0a   = tmp + 0x14;
		x0b   = tmp + 0x16;
		x0c   = tmp + 0x18;
		x0d   = tmp + 0x1a;
		x0e   = tmp + 0x1c;
	  #endif
		tmp += 0x1e;	// sc_ptr + 270
	  #ifndef MULTITHREAD
		y00   = tmp + 0x00;
		y01   = tmp + 0x02;
		y02   = tmp + 0x04;
		y03   = tmp + 0x06;
		y04   = tmp + 0x08;
		y05   = tmp + 0x0a;
		y06   = tmp + 0x0c;
		y07   = tmp + 0x0e;
		y08   = tmp + 0x10;
		y09   = tmp + 0x12;
		y0a   = tmp + 0x14;
		y0b   = tmp + 0x16;
		y0c   = tmp + 0x18;
		y0d   = tmp + 0x1a;
		y0e   = tmp + 0x1c;
	  #endif
		tmp += 0x1e;	// +30= 300
		sse2_c3m1 = tmp + 0x00;
		sse2_s    = tmp + 0x01;
		sse2_cn1  = tmp + 0x02;
		sse2_cn2  = tmp + 0x03;
		sse2_ss3  = tmp + 0x04;
		sse2_sn1  = tmp + 0x05;
		sse2_sn2  = tmp + 0x06;
		two       = tmp + 0x07;
		tmp += 0x08;	// +8 = 308
	#ifdef USE_AVX512
		cy_r = tmp;										// RADIX/8 (round up) vec_dbl slots for each of cy_r and cy_i carry sub-arrays
	  #ifndef MULTITHREAD
					cy_i = tmp+0x8;						// RADIX/8 (round up) vec_dbl slots for each of cy_r and cy_i carry sub-arrays
	  #endif
										tmp += 2*0x8;	// RADIX/8 (round up) vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	#elif defined(USE_AVX)
		cy_r = tmp;										// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays;  +30 = 338
	  #ifndef MULTITHREAD
					cy_i = tmp+0xf;						// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays;  +30 = 338
	  #endif
										tmp += 2*0xf;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays;  +30 = 338
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 30 + 2 = 340
		// This is where the value of half_arr_offset60 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	#else
		cy_r = tmp;	/* cy_i = tmp+0x1e; */	tmp += 2*0x1e;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays; +60 = 368 complex
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 60 + 2 = 370 complex
		// This is where the value of half_arr_offset60 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*ODD_RADIX] x 16 for Fermat-mod */
	#endif

		ASSERT(half_arr_offset60 == (uint32)(half_arr-sc_ptr), "half_arr_offset60 mismatches actual!");
		ASSERT((radix60_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r00) + (20 << L2_SZ_VD), "radix60_creals_in_local_store checksum failed!");
		/* These remain fixed: */
		VEC_DBL_INIT(sse2_c3m1, c3m1);
		VEC_DBL_INIT(sse2_s   , s   );
		VEC_DBL_INIT(sse2_cn1 , cn1 );
		VEC_DBL_INIT(sse2_cn2 , cn2 );
		VEC_DBL_INIT(sse2_ss3 , ss3 );
		VEC_DBL_INIT(sse2_sn1 , sn1 );
		VEC_DBL_INIT(sse2_sn2 , sn2 );
		VEC_DBL_INIT(two      , 2.0 );

		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		sse2_rnd->d0 = base[0]; sse2_rnd->d1 = baseinv[1]; sse2_rnd->d2 = wts_mult[1]; sse2_rnd->d3 = inv_mult[0];
	  #else
		VEC_DBL_INIT(sse2_rnd, crnd);
	  #endif

		// Propagate the above consts to the remaining threads:
		nbytes = (intptr_t)cy_r - (intptr_t)sse2_c3m1;	// #bytes in 1st of above block of consts
		tmp = sse2_c3m1;
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
	#ifdef USE_AVX

		base_negacyclic_root = half_arr + RADIX;

	  #ifdef HIACC
		#ifdef USE_AVX512
		  #error Fermat-mod HIACC mode only available in AVX/AVX2 builds!
		#endif
		/*
		The pattern of the negacyclic-DWT-weights ("nDWTs") applied to the RADIX complex outputs of the final-DIT-pass is like so:
		The nDWTs multiplying each set of RADIX DIT DFT outputs are simply the product of a single complex-root "base multiplier" rbase
		(separately computed for each batch of DFT outputs), which "base root" multiplies the (0 - RADIX-1)st [4*RADIX]th roots of unity,
		i.e.
			 rbase * (j*I*Pi/2)/RADIX, for j = 0, ..., RADIX-1 .

		See the radix28 version of this routine for additional details.
		*/
		#if 0
		// Simple qfloat-based loop to crank out the roots which populate the radix60_avx_negadwt_consts table:
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
		tmp64 = radix60_avx_negadwt_consts[1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(  1*I*Pi/(2*RADIX)) = sin((RADIX-  1)*I*Pi/(2*RADIX)) */
		tmp64 = radix60_avx_negadwt_consts[2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(  2*I*Pi/(2*RADIX)) = sin((RADIX-  2)*I*Pi/(2*RADIX)) */
		tmp64 = radix60_avx_negadwt_consts[3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(  3*I*Pi/(2*RADIX)) = sin((RADIX-  3)*I*Pi/(2*RADIX)) */	tmp += 2;
		for(j = 4; j < RADIX; j += 4) {
			tmp64 = radix60_avx_negadwt_consts[j+0];	tmp->d0 = tm2->d0 = *(double *)&tmp64;	tm2 -= 2;
			tmp64 = radix60_avx_negadwt_consts[j+1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;
			tmp64 = radix60_avx_negadwt_consts[j+2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;
			tmp64 = radix60_avx_negadwt_consts[j+3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	tmp += 2;
		}
		tmp = base_negacyclic_root + RADIX*2;	// reset to point to start of above block
		nbytes = RADIX*SZ_VD/2;	// RADIX/4 AVX-register-sized complex data

	  #else	// HIACC = false:

		// lower-precision version, which yields slightly more roundoff error, but is simpler and more storage-compact.
	   #ifdef USE_AVX512

		// Init exp(j*I*Pi/2/RADIX), for j = 0-7:
		tmp = base_negacyclic_root + 16;	// First 16 slots reserved for Re/Im parts of the 8 base multipliers
		tmp->d0 = 1.0;
		tmp64 = radix60_avx_negadwt_consts[      1];	tmp->d1 = *(double *)&tmp64;	// cos(01*I*Pi/(2*RADIX))
		tmp64 = radix60_avx_negadwt_consts[      2];	tmp->d2 = *(double *)&tmp64;	// cos(02*I*Pi/(2*RADIX))
		tmp64 = radix60_avx_negadwt_consts[      3];	tmp->d3 = *(double *)&tmp64;	// cos(03*I*Pi/(2*RADIX))
		tmp64 = radix60_avx_negadwt_consts[      4];	tmp->d4 = *(double *)&tmp64;	// cos(04*I*Pi/(2*RADIX))
		tmp64 = radix60_avx_negadwt_consts[      5];	tmp->d5 = *(double *)&tmp64;	// cos(05*I*Pi/(2*RADIX))
		tmp64 = radix60_avx_negadwt_consts[      6];	tmp->d6 = *(double *)&tmp64;	// cos(06*I*Pi/(2*RADIX))
		tmp64 = radix60_avx_negadwt_consts[      7];	tmp->d7 = *(double *)&tmp64;	// cos(07*I*Pi/(2*RADIX))
		++tmp;
		tmp->d0 = 0.0;
		tmp64 = radix60_avx_negadwt_consts[RADIX-1];	tmp->d1 = *(double *)&tmp64;	// sin(01*I*Pi/(2*RADIX))
		tmp64 = radix60_avx_negadwt_consts[RADIX-2];	tmp->d2 = *(double *)&tmp64;	// sin(02*I*Pi/(2*RADIX))
		tmp64 = radix60_avx_negadwt_consts[RADIX-3];	tmp->d3 = *(double *)&tmp64;	// sin(03*I*Pi/(2*RADIX))
		tmp64 = radix60_avx_negadwt_consts[RADIX-4];	tmp->d4 = *(double *)&tmp64;	// sin(04*I*Pi/(2*RADIX))
		tmp64 = radix60_avx_negadwt_consts[RADIX-5];	tmp->d5 = *(double *)&tmp64;	// sin(05*I*Pi/(2*RADIX))
		tmp64 = radix60_avx_negadwt_consts[RADIX-6];	tmp->d6 = *(double *)&tmp64;	// sin(06*I*Pi/(2*RADIX))
		tmp64 = radix60_avx_negadwt_consts[RADIX-7];	tmp->d7 = *(double *)&tmp64;	// sin(07*I*Pi/(2*RADIX))
		++tmp;	// 0x480(base_negacyclic_root)
		tmp64 = radix60_avx_negadwt_consts[      8];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// cos(08*I*Pi/(2*RADIX))
		++tmp;	// 0x4c0(base_negacyclic_root)
		tmp64 = radix60_avx_negadwt_consts[RADIX-8];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// sin(08*I*Pi/(2*RADIX))
		tmp = base_negacyclic_root + 16;	// reset to point to start of above block

	   #elif defined(USE_AVX)

		// Init exp(j*I*Pi/2/RADIX), for j = 0-3:
		tmp = base_negacyclic_root + 8;	// First 8 slots reserved for Re/Im parts of the 4 base multipliers
		tmp->d0 = 1.0;
		tmp64 = radix60_avx_negadwt_consts[      1];	tmp->d1 = *(double *)&tmp64;	// cos(01*I*Pi/(2*RADIX))
		tmp64 = radix60_avx_negadwt_consts[      2];	tmp->d2 = *(double *)&tmp64;	// cos(02*I*Pi/(2*RADIX))
		tmp64 = radix60_avx_negadwt_consts[      3];	tmp->d3 = *(double *)&tmp64;	// cos(03*I*Pi/(2*RADIX))

		(++tmp)->d0 = 0.0;
		tmp64 = radix60_avx_negadwt_consts[RADIX-1];	tmp->d1 = *(double *)&tmp64;	// sin(01*I*Pi/(2*RADIX))
		tmp64 = radix60_avx_negadwt_consts[RADIX-2];	tmp->d2 = *(double *)&tmp64;	// sin(02*I*Pi/(2*RADIX))
		tmp64 = radix60_avx_negadwt_consts[RADIX-3];	tmp->d3 = *(double *)&tmp64;	// sin(03*I*Pi/(2*RADIX))
		++tmp;
		tmp64 = radix60_avx_negadwt_consts[      4];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// cos(04*I*Pi/(2*RADIX))
		++tmp;
		tmp64 = radix60_avx_negadwt_consts[RADIX-4];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// sin(04*I*Pi/(2*RADIX))
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
	}		/************************************************************************/
	else	/*                MODULUS_TYPE_MERSENNE:                                */
	{		/************************************************************************/
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

	/*********** Defer the per-thread local-mem-block copy until after added wts-index precomputation below ************/
	#endif	/* USE_SSE2 */

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
		p44 = p40 + p04;
		p48 = p44 + p04;
		p52 = p48 + p04;
		p56 = p52 + p04;

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
		p44 = p44 + ( (p44 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p52 = p52 + ( (p52 >> DAT_BITS) << PAD_BITS );
		p56 = p56 + ( (p56 >> DAT_BITS) << PAD_BITS );
	#if !defined(MULTITHREAD) && (!defined(USE_SSE2) || (COMPACT_OBJ != 0))
		p0123[0] = 0; p0123[1] = p01; p0123[2] = p02; p0123[3] = p03;
	#endif
		poff[0x0] =   0; poff[0x1] = p04; poff[0x2] = p08; poff[0x3] = p12;
		poff[0x4] = p16; poff[0x5] = p20; poff[0x6] = p24; poff[0x7] = p28;
		poff[0x8] = p32; poff[0x9] = p36; poff[0xa] = p40; poff[0xb] = p44;
		poff[0xc] = p48; poff[0xd] = p52; poff[0xe] = p56;

		if(_cy_r[0])	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			ASSERT(0 != _i, "free(_i) but ptr = 0x0!");
			for(i= 0; i < RADIX; i++) {
				ASSERT(0 != _bjmodn[i], "free(_bjmodn[i]) but ptr = 0x0!");
				ASSERT(0 !=   _cy_r[i], "free(_cy_r[i]) but ptr = 0x0!");
				ASSERT(0 !=   _cy_i[i], "free(_cy_i[i]) but ptr = 0x0!");
			}
			ASSERT(0 != _jstart, "free(_jstart) but ptr = 0x0!");
			ASSERT(0 != _jhi, "free(_jhi) but ptr = 0x0!");
			ASSERT(0 != _col, "free(_col) but ptr = 0x0!");
			ASSERT(0 != _co2, "free(_co2) but ptr = 0x0!");
			ASSERT(0 != _co3, "free(_co3) but ptr = 0x0!");
			ASSERT(0 != _bjmodnini, "free(_bjmodnini) but ptr = 0x0!");

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

		ASSERT(ptr_prod == 0, "ERROR: unable to allocate one or more auxiliary arrays.");

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
			}
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
			tmp->d0 = bs_arr[0x0];	tmp->d1 = bs_arr[0x1];	tmp->d2 = bs_arr[0x2];	tmp->d3 = bs_arr[0x3];	tmp->d4 = bs_arr[0x4];	tmp->d5 = bs_arr[0x5];	tmp->d6 = bs_arr[0x6];	tmp->d7 = bs_arr[0x7];	++tmp;
			tmp->d0 = bs_arr[0x1];	tmp->d1 = bs_arr[0x2];	tmp->d2 = bs_arr[0x3];	tmp->d3 = bs_arr[0x4];	tmp->d4 = bs_arr[0x5];	tmp->d5 = bs_arr[0x6];	tmp->d6 = bs_arr[0x7];	tmp->d7 = bs_arr[0x8];	++tmp;
			tmp->d0 = bs_arr[0x2];	tmp->d1 = bs_arr[0x3];	tmp->d2 = bs_arr[0x4];	tmp->d3 = bs_arr[0x5];	tmp->d4 = bs_arr[0x6];	tmp->d5 = bs_arr[0x7];	tmp->d6 = bs_arr[0x8];	tmp->d7 = bs_arr[0x9];	++tmp;
			tmp->d0 = bs_arr[0x3];	tmp->d1 = bs_arr[0x4];	tmp->d2 = bs_arr[0x5];	tmp->d3 = bs_arr[0x6];	tmp->d4 = bs_arr[0x7];	tmp->d5 = bs_arr[0x8];	tmp->d6 = bs_arr[0x9];	tmp->d7 = bs_arr[0xa];	++tmp;
			tmp->d0 = bs_arr[0x4];	tmp->d1 = bs_arr[0x5];	tmp->d2 = bs_arr[0x6];	tmp->d3 = bs_arr[0x7];	tmp->d4 = bs_arr[0x8];	tmp->d5 = bs_arr[0x9];	tmp->d6 = bs_arr[0xa];	tmp->d7 = bs_arr[0xb];	++tmp;
			tmp->d0 = bs_arr[0x5];	tmp->d1 = bs_arr[0x6];	tmp->d2 = bs_arr[0x7];	tmp->d3 = bs_arr[0x8];	tmp->d4 = bs_arr[0x9];	tmp->d5 = bs_arr[0xa];	tmp->d6 = bs_arr[0xb];	tmp->d7 = bs_arr[0xc];	++tmp;
			tmp->d0 = bs_arr[0x6];	tmp->d1 = bs_arr[0x7];	tmp->d2 = bs_arr[0x8];	tmp->d3 = bs_arr[0x9];	tmp->d4 = bs_arr[0xa];	tmp->d5 = bs_arr[0xb];	tmp->d6 = bs_arr[0xc];	tmp->d7 = bs_arr[0xd];	++tmp;
			tmp->d0 = bs_arr[0x7];	tmp->d1 = bs_arr[0x8];	tmp->d2 = bs_arr[0x9];	tmp->d3 = bs_arr[0xa];	tmp->d4 = bs_arr[0xb];	tmp->d5 = bs_arr[0xc];	tmp->d6 = bs_arr[0xd];	tmp->d7 = bs_arr[0xe];	++tmp;
			tmp->d0 = bs_arr[0x8];	tmp->d1 = bs_arr[0x9];	tmp->d2 = bs_arr[0xa];	tmp->d3 = bs_arr[0xb];	tmp->d4 = bs_arr[0xc];	tmp->d5 = bs_arr[0xd];	tmp->d6 = bs_arr[0xe];	tmp->d7 = bs_arr[0x0];	++tmp;
			tmp->d0 = bs_arr[0x9];	tmp->d1 = bs_arr[0xa];	tmp->d2 = bs_arr[0xb];	tmp->d3 = bs_arr[0xc];	tmp->d4 = bs_arr[0xd];	tmp->d5 = bs_arr[0xe];	tmp->d6 = bs_arr[0x0];	tmp->d7 = bs_arr[0x1];	++tmp;
			tmp->d0 = bs_arr[0xa];	tmp->d1 = bs_arr[0xb];	tmp->d2 = bs_arr[0xc];	tmp->d3 = bs_arr[0xd];	tmp->d4 = bs_arr[0xe];	tmp->d5 = bs_arr[0x0];	tmp->d6 = bs_arr[0x1];	tmp->d7 = bs_arr[0x2];	++tmp;
			tmp->d0 = bs_arr[0xb];	tmp->d1 = bs_arr[0xc];	tmp->d2 = bs_arr[0xd];	tmp->d3 = bs_arr[0xe];	tmp->d4 = bs_arr[0x0];	tmp->d5 = bs_arr[0x1];	tmp->d6 = bs_arr[0x2];	tmp->d7 = bs_arr[0x3];	++tmp;
			tmp->d0 = bs_arr[0xc];	tmp->d1 = bs_arr[0xd];	tmp->d2 = bs_arr[0xe];	tmp->d3 = bs_arr[0x0];	tmp->d4 = bs_arr[0x1];	tmp->d5 = bs_arr[0x2];	tmp->d6 = bs_arr[0x3];	tmp->d7 = bs_arr[0x4];	++tmp;
			tmp->d0 = bs_arr[0xd];	tmp->d1 = bs_arr[0xe];	tmp->d2 = bs_arr[0x0];	tmp->d3 = bs_arr[0x1];	tmp->d4 = bs_arr[0x2];	tmp->d5 = bs_arr[0x3];	tmp->d6 = bs_arr[0x4];	tmp->d7 = bs_arr[0x5];	++tmp;
			tmp->d0 = bs_arr[0xe];	tmp->d1 = bs_arr[0x0];	tmp->d2 = bs_arr[0x1];	tmp->d3 = bs_arr[0x2];	tmp->d4 = bs_arr[0x3];	tmp->d5 = bs_arr[0x4];	tmp->d6 = bs_arr[0x5];	tmp->d7 = bs_arr[0x6];	++tmp;

			tmp->d0 = bsinv_arr[0x0];	tmp->d1 = bsinv_arr[0x1];	tmp->d2 = bsinv_arr[0x2];	tmp->d3 = bsinv_arr[0x3];	tmp->d4 = bsinv_arr[0x4];	tmp->d5 = bsinv_arr[0x5];	tmp->d6 = bsinv_arr[0x6];	tmp->d7 = bsinv_arr[0x7];	++tmp;
			tmp->d0 = bsinv_arr[0x1];	tmp->d1 = bsinv_arr[0x2];	tmp->d2 = bsinv_arr[0x3];	tmp->d3 = bsinv_arr[0x4];	tmp->d4 = bsinv_arr[0x5];	tmp->d5 = bsinv_arr[0x6];	tmp->d6 = bsinv_arr[0x7];	tmp->d7 = bsinv_arr[0x8];	++tmp;
			tmp->d0 = bsinv_arr[0x2];	tmp->d1 = bsinv_arr[0x3];	tmp->d2 = bsinv_arr[0x4];	tmp->d3 = bsinv_arr[0x5];	tmp->d4 = bsinv_arr[0x6];	tmp->d5 = bsinv_arr[0x7];	tmp->d6 = bsinv_arr[0x8];	tmp->d7 = bsinv_arr[0x9];	++tmp;
			tmp->d0 = bsinv_arr[0x3];	tmp->d1 = bsinv_arr[0x4];	tmp->d2 = bsinv_arr[0x5];	tmp->d3 = bsinv_arr[0x6];	tmp->d4 = bsinv_arr[0x7];	tmp->d5 = bsinv_arr[0x8];	tmp->d6 = bsinv_arr[0x9];	tmp->d7 = bsinv_arr[0xa];	++tmp;
			tmp->d0 = bsinv_arr[0x4];	tmp->d1 = bsinv_arr[0x5];	tmp->d2 = bsinv_arr[0x6];	tmp->d3 = bsinv_arr[0x7];	tmp->d4 = bsinv_arr[0x8];	tmp->d5 = bsinv_arr[0x9];	tmp->d6 = bsinv_arr[0xa];	tmp->d7 = bsinv_arr[0xb];	++tmp;
			tmp->d0 = bsinv_arr[0x5];	tmp->d1 = bsinv_arr[0x6];	tmp->d2 = bsinv_arr[0x7];	tmp->d3 = bsinv_arr[0x8];	tmp->d4 = bsinv_arr[0x9];	tmp->d5 = bsinv_arr[0xa];	tmp->d6 = bsinv_arr[0xb];	tmp->d7 = bsinv_arr[0xc];	++tmp;
			tmp->d0 = bsinv_arr[0x6];	tmp->d1 = bsinv_arr[0x7];	tmp->d2 = bsinv_arr[0x8];	tmp->d3 = bsinv_arr[0x9];	tmp->d4 = bsinv_arr[0xa];	tmp->d5 = bsinv_arr[0xb];	tmp->d6 = bsinv_arr[0xc];	tmp->d7 = bsinv_arr[0xd];	++tmp;
			tmp->d0 = bsinv_arr[0x7];	tmp->d1 = bsinv_arr[0x8];	tmp->d2 = bsinv_arr[0x9];	tmp->d3 = bsinv_arr[0xa];	tmp->d4 = bsinv_arr[0xb];	tmp->d5 = bsinv_arr[0xc];	tmp->d6 = bsinv_arr[0xd];	tmp->d7 = bsinv_arr[0xe];	++tmp;
			tmp->d0 = bsinv_arr[0x8];	tmp->d1 = bsinv_arr[0x9];	tmp->d2 = bsinv_arr[0xa];	tmp->d3 = bsinv_arr[0xb];	tmp->d4 = bsinv_arr[0xc];	tmp->d5 = bsinv_arr[0xd];	tmp->d6 = bsinv_arr[0xe];	tmp->d7 = bsinv_arr[0x0];	++tmp;
			tmp->d0 = bsinv_arr[0x9];	tmp->d1 = bsinv_arr[0xa];	tmp->d2 = bsinv_arr[0xb];	tmp->d3 = bsinv_arr[0xc];	tmp->d4 = bsinv_arr[0xd];	tmp->d5 = bsinv_arr[0xe];	tmp->d6 = bsinv_arr[0x0];	tmp->d7 = bsinv_arr[0x1];	++tmp;
			tmp->d0 = bsinv_arr[0xa];	tmp->d1 = bsinv_arr[0xb];	tmp->d2 = bsinv_arr[0xc];	tmp->d3 = bsinv_arr[0xd];	tmp->d4 = bsinv_arr[0xe];	tmp->d5 = bsinv_arr[0x0];	tmp->d6 = bsinv_arr[0x1];	tmp->d7 = bsinv_arr[0x2];	++tmp;
			tmp->d0 = bsinv_arr[0xb];	tmp->d1 = bsinv_arr[0xc];	tmp->d2 = bsinv_arr[0xd];	tmp->d3 = bsinv_arr[0xe];	tmp->d4 = bsinv_arr[0x0];	tmp->d5 = bsinv_arr[0x1];	tmp->d6 = bsinv_arr[0x2];	tmp->d7 = bsinv_arr[0x3];	++tmp;
			tmp->d0 = bsinv_arr[0xc];	tmp->d1 = bsinv_arr[0xd];	tmp->d2 = bsinv_arr[0xe];	tmp->d3 = bsinv_arr[0x0];	tmp->d4 = bsinv_arr[0x1];	tmp->d5 = bsinv_arr[0x2];	tmp->d6 = bsinv_arr[0x3];	tmp->d7 = bsinv_arr[0x4];	++tmp;
			tmp->d0 = bsinv_arr[0xd];	tmp->d1 = bsinv_arr[0xe];	tmp->d2 = bsinv_arr[0x0];	tmp->d3 = bsinv_arr[0x1];	tmp->d4 = bsinv_arr[0x2];	tmp->d5 = bsinv_arr[0x3];	tmp->d6 = bsinv_arr[0x4];	tmp->d7 = bsinv_arr[0x5];	++tmp;
			tmp->d0 = bsinv_arr[0xe];	tmp->d1 = bsinv_arr[0x0];	tmp->d2 = bsinv_arr[0x1];	tmp->d3 = bsinv_arr[0x2];	tmp->d4 = bsinv_arr[0x3];	tmp->d5 = bsinv_arr[0x4];	tmp->d6 = bsinv_arr[0x5];	tmp->d7 = bsinv_arr[0x6];	++tmp;

		  #elif defined(USE_AVX)

			// Each transposed-data quartet in the AVX carry macro needs linearly incrementing bs_arr data (mod ODD_RADIX);
			// Need all [ODD_RADIX] possible such length-4 index subsequences, which will be accessed via their head element
			// by the [ijkl]cycle* index quartets in the respective carry-macro call:
			tmp->d0 = bs_arr[0x0];	tmp->d1 = bs_arr[0x1];	tmp->d2 = bs_arr[0x2];	tmp->d3 = bs_arr[0x3];	++tmp;
			tmp->d0 = bs_arr[0x1];	tmp->d1 = bs_arr[0x2];	tmp->d2 = bs_arr[0x3];	tmp->d3 = bs_arr[0x4];	++tmp;
			tmp->d0 = bs_arr[0x2];	tmp->d1 = bs_arr[0x3];	tmp->d2 = bs_arr[0x4];	tmp->d3 = bs_arr[0x5];	++tmp;
			tmp->d0 = bs_arr[0x3];	tmp->d1 = bs_arr[0x4];	tmp->d2 = bs_arr[0x5];	tmp->d3 = bs_arr[0x6];	++tmp;
			tmp->d0 = bs_arr[0x4];	tmp->d1 = bs_arr[0x5];	tmp->d2 = bs_arr[0x6];	tmp->d3 = bs_arr[0x7];	++tmp;
			tmp->d0 = bs_arr[0x5];	tmp->d1 = bs_arr[0x6];	tmp->d2 = bs_arr[0x7];	tmp->d3 = bs_arr[0x8];	++tmp;
			tmp->d0 = bs_arr[0x6];	tmp->d1 = bs_arr[0x7];	tmp->d2 = bs_arr[0x8];	tmp->d3 = bs_arr[0x9];	++tmp;
			tmp->d0 = bs_arr[0x7];	tmp->d1 = bs_arr[0x8];	tmp->d2 = bs_arr[0x9];	tmp->d3 = bs_arr[0xa];	++tmp;
			tmp->d0 = bs_arr[0x8];	tmp->d1 = bs_arr[0x9];	tmp->d2 = bs_arr[0xa];	tmp->d3 = bs_arr[0xb];	++tmp;
			tmp->d0 = bs_arr[0x9];	tmp->d1 = bs_arr[0xa];	tmp->d2 = bs_arr[0xb];	tmp->d3 = bs_arr[0xc];	++tmp;
			tmp->d0 = bs_arr[0xa];	tmp->d1 = bs_arr[0xb];	tmp->d2 = bs_arr[0xc];	tmp->d3 = bs_arr[0xd];	++tmp;
			tmp->d0 = bs_arr[0xb];	tmp->d1 = bs_arr[0xc];	tmp->d2 = bs_arr[0xd];	tmp->d3 = bs_arr[0xe];	++tmp;
			tmp->d0 = bs_arr[0xc];	tmp->d1 = bs_arr[0xd];	tmp->d2 = bs_arr[0xe];	tmp->d3 = bs_arr[0x0];	++tmp;
			tmp->d0 = bs_arr[0xd];	tmp->d1 = bs_arr[0xe];	tmp->d2 = bs_arr[0x0];	tmp->d3 = bs_arr[0x1];	++tmp;
			tmp->d0 = bs_arr[0xe];	tmp->d1 = bs_arr[0x0];	tmp->d2 = bs_arr[0x1];	tmp->d3 = bs_arr[0x2];	++tmp;

			tmp->d0 = bsinv_arr[0x0];	tmp->d1 = bsinv_arr[0x1];	tmp->d2 = bsinv_arr[0x2];	tmp->d3 = bsinv_arr[0x3];	++tmp;
			tmp->d0 = bsinv_arr[0x1];	tmp->d1 = bsinv_arr[0x2];	tmp->d2 = bsinv_arr[0x3];	tmp->d3 = bsinv_arr[0x4];	++tmp;
			tmp->d0 = bsinv_arr[0x2];	tmp->d1 = bsinv_arr[0x3];	tmp->d2 = bsinv_arr[0x4];	tmp->d3 = bsinv_arr[0x5];	++tmp;
			tmp->d0 = bsinv_arr[0x3];	tmp->d1 = bsinv_arr[0x4];	tmp->d2 = bsinv_arr[0x5];	tmp->d3 = bsinv_arr[0x6];	++tmp;
			tmp->d0 = bsinv_arr[0x4];	tmp->d1 = bsinv_arr[0x5];	tmp->d2 = bsinv_arr[0x6];	tmp->d3 = bsinv_arr[0x7];	++tmp;
			tmp->d0 = bsinv_arr[0x5];	tmp->d1 = bsinv_arr[0x6];	tmp->d2 = bsinv_arr[0x7];	tmp->d3 = bsinv_arr[0x8];	++tmp;
			tmp->d0 = bsinv_arr[0x6];	tmp->d1 = bsinv_arr[0x7];	tmp->d2 = bsinv_arr[0x8];	tmp->d3 = bsinv_arr[0x9];	++tmp;
			tmp->d0 = bsinv_arr[0x7];	tmp->d1 = bsinv_arr[0x8];	tmp->d2 = bsinv_arr[0x9];	tmp->d3 = bsinv_arr[0xa];	++tmp;
			tmp->d0 = bsinv_arr[0x8];	tmp->d1 = bsinv_arr[0x9];	tmp->d2 = bsinv_arr[0xa];	tmp->d3 = bsinv_arr[0xb];	++tmp;
			tmp->d0 = bsinv_arr[0x9];	tmp->d1 = bsinv_arr[0xa];	tmp->d2 = bsinv_arr[0xb];	tmp->d3 = bsinv_arr[0xc];	++tmp;
			tmp->d0 = bsinv_arr[0xa];	tmp->d1 = bsinv_arr[0xb];	tmp->d2 = bsinv_arr[0xc];	tmp->d3 = bsinv_arr[0xd];	++tmp;
			tmp->d0 = bsinv_arr[0xb];	tmp->d1 = bsinv_arr[0xc];	tmp->d2 = bsinv_arr[0xd];	tmp->d3 = bsinv_arr[0xe];	++tmp;
			tmp->d0 = bsinv_arr[0xc];	tmp->d1 = bsinv_arr[0xd];	tmp->d2 = bsinv_arr[0xe];	tmp->d3 = bsinv_arr[0x0];	++tmp;
			tmp->d0 = bsinv_arr[0xd];	tmp->d1 = bsinv_arr[0xe];	tmp->d2 = bsinv_arr[0x0];	tmp->d3 = bsinv_arr[0x1];	++tmp;
			tmp->d0 = bsinv_arr[0xe];	tmp->d1 = bsinv_arr[0x0];	tmp->d2 = bsinv_arr[0x1];	tmp->d3 = bsinv_arr[0x2];	++tmp;

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
			tdat[ithread].r00    = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (vec_dbl *)((intptr_t)tdat[ithread].r00 + ((intptr_t)half_arr - (intptr_t)r00));
		#else
			// In scalar mode use these 2 ptrs to pass wts_idx_incr and the base/baseinv/etc array-ptrs:
			tdat[ithread].r00    = (double *)foo_array;
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
				j = ((int64)wts_idx_incr * ( (jhi-jstart)>>(L2_SZ_VD-2) ) % nwt  );	// []>>(L2_SZ_VD-2) is fast subst. for []/stride
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

/*...The radix-60 final DIT pass is here.	*/

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
				_jhi[ithread] = _jstart[ithread] + jhi_wrap_mers;
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
				_jhi[ithread] = _jstart[ithread] + jhi_wrap_ferm;
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
		tmp = tdat[ithread].r00;
		ASSERT(((tmp + 300)->d0 == c3m1 && (tmp + 300)->d1 == c3m1), "thread-local memcheck failed!");
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
		#include "radix60_main_carry_loop.h"

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
		ASSERT(0x0 == cy60_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
	ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
	ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(0 == mlucas_nanosleep(&ns_time), "nanosleep fail!");
	}
//	printf("radix32_ditN_cy_dif1 end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

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

	// For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large:
	if(TRANSFORM_TYPE == RIGHT_ANGLE)
		j_jhi = jhi_wrap_ferm;
	else
		j_jhi = jhi_wrap_mers;

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

void radix60_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-60 complex DIF FFT pass on the data in the length-N real vector A.
*/
	int j,jt,jp,j1,j2;
	static int NDIVR,first_entry=TRUE,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40,p44,p48,p52,p56;
	static double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/
					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it,
		t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,
		t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,
		t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,
		t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei;

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
		p44 = p40 + p04;
		p48 = p44 + p04;
		p52 = p48 + p04;
		p56 = p52 + p04;

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
		p44 = p44 + ( (p44 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p52 = p52 + ( (p52 >> DAT_BITS) << PAD_BITS );
		p56 = p56 + ( (p56 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-60 pass is here.	*/

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

		/*...gather the needed data (60 64-bit complex, i.e.120 64-bit reals) and do 4 radix-15 transforms...*/
	/*
	Twiddleless version arranges 4 sets of radix-15 DFT inputs as follows:
	0 in upper left corner, decrement 4 horizontally and 15 vertically, indexing modulo 60:

		RADIX_15_DFT(00,56,52,48,44,40,36,32,28,24,20,16,12,08,04)
		RADIX_15_DFT(45,41,37,33,29,25,21,17,13,09,05,01,57,53,49)
		RADIX_15_DFT(30,26,22,18,14,10,06,02,58,54,50,46,42,38,34)
		RADIX_15_DFT(15,11,07,03,59,55,51,47,43,39,35,31,27,23,19)

	If we subtract costant offset 1,2,3 from the last 3 rows as we do in the implementation to reduce the number
	of index offsets needing to be stored, we decrement 4 horizontally and 16 vertically.
	*/
	#if 0
		jt = j1    ; jp = j2    ;	RADIX_15_DIF(a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,rt,it);
		jt = j1+p01; jp = j2+p01;	RADIX_15_DIF(a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,rt,it);
		jt = j1+p02; jp = j2+p02;	RADIX_15_DIF(a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,rt,it);
		jt = j1+p03; jp = j2+p03;	RADIX_15_DIF(a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei,rt,it);
	#else
		jt = j1    ; jp = j2    ;	RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04], t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei);
		jt = j1+p01; jp = j2+p01;	RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48], t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei);
		jt = j1+p02; jp = j2+p02;	RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32], t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei);
		jt = j1+p03; jp = j2+p03;	RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16], t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei);
	#endif
		/*...and now do 15 radix-4 transforms.
		The required output permutation is:

			[ 0, 1, 2, 3
			, 9, 8,11,10
			, 6, 7, 5, 4
			,57,56,59,58
			,54,55,53,52
			,48,49,50,51
			,46,47,45,44
			,40,41,42,43
			,39,38,36,37
			,32,33,34,35
			,31,30,28,29
			,25,24,27,26
			,23,22,20,21
			,17,16,19,18
			,14,15,13,12].
		*/
												/*          inputs                    */ /*                                      outputs                              */
		jt = j1    ; jp = j2    ;	RADIX_04_DIF(t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIF(t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIF(t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
		jt = j1+p56; jp = j2+p56;	RADIX_04_DIF(t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
		jt = j1+p52; jp = j2+p52;	RADIX_04_DIF(t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
		jt = j1+p48; jp = j2+p48;	RADIX_04_DIF(t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
		jt = j1+p44; jp = j2+p44;	RADIX_04_DIF(t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
		jt = j1+p40; jp = j2+p40;	RADIX_04_DIF(t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
		jt = j1+p36; jp = j2+p36;	RADIX_04_DIF(t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIF(t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIF(t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIF(t0br,t0bi,t1br,t1bi,t2br,t2bi,t3br,t3bi,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIF(t0cr,t0ci,t1cr,t1ci,t2cr,t2ci,t3cr,t3ci,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIF(t0dr,t0di,t1dr,t1di,t2dr,t2di,t3dr,t3di,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIF(t0er,t0ei,t1er,t1ei,t2er,t2ei,t3er,t3ei,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
		/* Totals: 4*radix15 + 15*radix04 = 4*(162 FADD, 50 FMUL) + 15*(16 FADD, 0 FMUL) = 888 FADD, 200 FMUL	*/
	}
}

/**************/

void radix60_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-60 complex DIT FFT pass on the data in the length-N real vector A.
*/
	int j,jt,jp,j1,j2;
	static int NDIVR,first_entry=TRUE,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40,p44,p48,p52,p56;
	static double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/
					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it,
		t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,
		t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,
		t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,
		t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei;

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
		p44 = p40 + p04;
		p48 = p44 + p04;
		p52 = p48 + p04;
		p56 = p52 + p04;

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
		p44 = p44 + ( (p44 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p52 = p52 + ( (p52 >> DAT_BITS) << PAD_BITS );
		p56 = p56 + ( (p56 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-60 pass is here.	*/

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
		Twiddleless version requires us to swap inputs as follows:
		indices x[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59
			  ->x[0 56,52,48,44,40,36,32,28,24,20,16,12,08,04|45,41,37,33,29,25,21,17,13,09,05,01,57,53,49|30,26,22,18,14,10,06,02,58,54,50,46,42,38,34|15,11,07,03,59,55,51,47,43,39,35,31,27,23,19

		I.e. start out with first quartet of indices {0,15,30,45}, permute those according to
		(0,15,30,45}*59%60 = {0,45,30,15), then each is head of a length-15 list of indices with decrement 4.

		Remember, inputs to DIT are bit-reversed, so
		a[ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59] contain
		x[ 0 30 15 45 05 35 20 50 10 40 25 55 01 31 16 46 06 36 21 51 11 41 26 56 02 32 17 47 07 37 22 52 12 42 27 57 03 33 18 48 08 38 23 53 13 43 28 58 04 34 19 49 09 39 24 54 14 44 29 59], which get swapped to
	(Look at first nontrivial one, i.e. x[1 -> 56] ... in terms of a[] this translates to a[12 -> 23])
	-->	x[ 0 30 45 15 40 10 25 55 20 50 05 35 56 26 41 11 36 06 21 51 16 46 01 31 52 22 37 07 32 02 17 47 12 42 57 27 48 18 33 03 28 58 13 43 08 38 53 23 44 14 29 59 24 54 09 39 04 34 49 19]
		a[ 0  1  3  2| 9  8 10 11| 6  7  4  5|23 22 21 20|17 16 18 19|14 15 12 13|31 30 29 28|25 24 26 27|32 33 35 34|39 38 37 36|46 47 44 45|40 41 43 42|57 56 58 59|54 55 52 53|48 49 51 50].
	*/
		/*...gather the needed data (60 64-bit complex, i.e.120 64-bit reals) and do 15 radix-4 transforms...*/
												/*                                   inputs                                   */ /*              outputs              */
		jt = j1    ; jp = j2    ;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,rt,it);
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,rt,it);
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,rt,it);
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,rt,it);
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,rt,it);
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,rt,it);
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,rt,it);
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,rt,it);
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,rt,it);
		jt = j1+p36; jp = j2+p36;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,rt,it);
		jt = j1+p44; jp = j2+p44;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,rt,it);
		jt = j1+p40; jp = j2+p40;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0br,t0bi,t1br,t1bi,t2br,t2bi,t3br,t3bi,rt,it);
		jt = j1+p56; jp = j2+p56;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t0cr,t0ci,t1cr,t1ci,t2cr,t2ci,t3cr,t3ci,rt,it);
		jt = j1+p52; jp = j2+p52;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t0dr,t0di,t1dr,t1di,t2dr,t2di,t3dr,t3di,rt,it);
		jt = j1+p48; jp = j2+p48;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0er,t0ei,t1er,t1ei,t2er,t2ei,t3er,t3ei,rt,it);

		/*...and now do 4 radix-15 transforms.
		The required output permutation is:

			[ 0,44,28,12,56,40,24, 8,52,36,20, 4,48,32,16
			 15,59,43,27,11,55,39,23, 7,51,35,19, 3,47,31
			 30,14,58,42,26,10,54,38,22, 6,50,34,18, 2,46
			 45,29,13,57,41,25, 9,53,37,21, 5,49,33,17, 1]
		*/
												/*                                                  inputs                                                                                          */ /*                                                                    temps                                                                        */ /* outputs: --> */
	#if 0
		jt = j1    ; jp = j2    ;	RADIX_15_DIT_C(t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a[jt    ],a[jp    ],a[jt+p44],a[jp+p44],a[jt+p28],a[jp+p28],a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],rt,it);
		jt = j1+p03; jp = j2+p03;	RADIX_15_DIT_C(t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p44],a[jp+p44],a[jt+p28],a[jp+p28],rt,it);
		jt = j1+p02; jp = j2+p02;	RADIX_15_DIT_C(t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a[jt+p28],a[jp+p28],a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p44],a[jp+p44],rt,it);
		jt = j1+p01; jp = j2+p01;	RADIX_15_DIT_C(t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a[jt+p44],a[jp+p44],a[jt+p28],a[jp+p28],a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],rt,it);
	#else
		jt = j1    ; jp = j2    ;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei, a[jt    ],a[jp    ],a[jt+p44],a[jp+p44],a[jt+p28],a[jp+p28],a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16]);
		jt = j1+p03; jp = j2+p03;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei, a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p44],a[jp+p44],a[jt+p28],a[jp+p28]);
		jt = j1+p02; jp = j2+p02;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei, a[jt+p28],a[jp+p28],a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p44],a[jp+p44]);
		jt = j1+p01; jp = j2+p01;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2, t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei, a[jt+p44],a[jp+p44],a[jt+p28],a[jp+p28],a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ]);
	#endif
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy60_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr,*addi;
	#ifndef USE_SSE2
		struct complex *tptr;
	#endif
	#ifdef USE_SSE2
		const int pfetch_dist = PFETCH_DIST;
	#endif
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40,p44,p48,p52,p56;
		int poff[RADIX>>2];
		int j,j1,l,ntmp;
	#ifndef USE_SSE2
		int j2;
	#else // USE_SSE2
		// incr = Carry-chain wts-multipliers recurrence length, which must divide
		// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = --|--|15 for avx512,avx,sse, respectively.
		// For the present case we structure the carry-macro calls such that this divisibility requirement is only for 128-bit SIMD builds;
		// note fixed-incr too restrictive here, so in sse2 case 'divide 15 into pieces' via increment-array whose elts sum to 15:
	  #ifndef USE_AVX
		const int *incr;
	  #endif
		const int *inc_arr = 0x0;	// Fermat-mod only uses inc_arr in AVX+ mode, init = 0 to quiet 'uninit' warnings.
		const int incr_long[] = {15}, incr_med[] = {8,7}, incr_short[] = {5,5,5};
	  // Have no specialized HIACC carry macro in USE_AVX512 and ARMv8 SIMD, use nonzero incr-value to differenetiate vs AVX/AVX2:
	  #if defined(USE_AVX512) || defined(USE_ARM_V8_SIMD)
		const int incr_hiacc[] = {3,3,3,3,3};	// Don't actually use the elements in AVX-512 mode; values are for ARMv8
	  #else
		const int incr_hiacc[] = {0};
	  #endif
		// Fermat-mod: For AVX+, define carry-subchain length in terms of 2^nfold subchains:
	  #ifdef USE_AVX512
		/*** No avx-512 support for radix-60 Fermat-mod, so just declare these with 0-values ***/
		int nfold = 0;
		const int nexec_long[] = {0}, nexec_med[] = {0}, nexec_short[] = {0}, nexec_hiacc[] = {0};
	  #elif defined(USE_AVX)
		// For nfold > 4,  RADIX/4 not divisible by 2^nfold, so use a more-general inner-loop scheme which can handle that:
		int nfold = USE_SHORT_CY_CHAIN + 0;
		const int nexec_long[] = {15}, nexec_med[] = {8,7}, nexec_short[] = {5,5,5}, nexec_hiacc[] = {3,3,3,3,3};
	  #endif
		// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
			if(USE_SHORT_CY_CHAIN == 0)
				inc_arr = incr_long;
			else if(USE_SHORT_CY_CHAIN == 1)
				inc_arr = incr_med;
			else if(USE_SHORT_CY_CHAIN == 2)
				inc_arr = incr_short;
			else
				inc_arr = incr_hiacc;
		} else {	// MODULUS_TYPE_FERMAT:
		#ifdef USE_AVX
			if(USE_SHORT_CY_CHAIN == 0)
				inc_arr = nexec_long;
			else if(USE_SHORT_CY_CHAIN == 1)
				inc_arr = nexec_med;
			else if(USE_SHORT_CY_CHAIN == 2)
				inc_arr = nexec_short;
			else
				inc_arr = nexec_hiacc;
		#endif
		}
	#endif // USE_SSE2

		uint32 k1,k2;
	#if defined(USE_SSE2) && !defined(USE_AVX512)
		uint32 k3,k4;
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
	#if !defined(USE_SSE2) || defined(USE_AVX)
		double wt_re,wt_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	#endif
	#ifndef USE_SSE2
		double wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	#endif
	#if !defined(USE_SSE2) || defined(USE_AVX)
		double rt,it;
	#endif
	#if !defined(USE_SSE2) || (COMPACT_OBJ != 0)
		uint32 p0123[4];
	#endif
	#ifdef USE_SSE2
	  #ifdef USE_AVX
		int i0,i1,i2,i3;
	  #endif
		const double c3m1= -1.50000000000000000000;	/* cos(twopi/3)-1	*/
	  #ifndef USE_AVX512
		const double crnd = 3.0*0x4000000*0x2000000;
	  #endif
		double *add0, *add1, *add2, *add3;
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *tmp,*tm1,*tm2;	// utility ptrs
	  #ifdef USE_AVX
		vec_dbl *tm0;	// utility ptrs
	  #endif
		int *itmp;			// Pointer into the bjmodn array
	  #if defined(USE_AVX) && !defined(USE_AVX512)
		int *itm2;			// Pointer into the bjmodn array
	  #endif
	  #ifndef USE_AVX
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	  #endif
		vec_dbl *max_err,
		  #ifndef USE_AVX512
			*sse2_rnd,
		  #endif
			*half_arr,
		  #ifndef USE_ARM_V8_SIMD // TODO that looks strange that this is needed to get rid of unused 'two'
			*two,
		  #endif
			*sse2_c3m1, /* *sse2_s, */ *sse2_cn1, /* *sse2_cn2, *sse2_ss3, *sse2_sn1, *sse2_sn2, */
			*s1p00,
		  #if !COMPACT_OBJ
			*s1p01,*s1p02,*s1p03,*s1p04,*s1p05,*s1p06,*s1p07,*s1p08,*s1p09,*s1p0a,*s1p0b,*s1p0c,*s1p0d,*s1p0e,
			*s1p10,*s1p11,*s1p12,*s1p13,*s1p14,*s1p15,*s1p16,*s1p17,*s1p18,*s1p19,*s1p1a,*s1p1b,*s1p1c,*s1p1d,*s1p1e,
			*s1p20,*s1p21,*s1p22,*s1p23,*s1p24,*s1p25,*s1p26,*s1p27,*s1p28,*s1p29,*s1p2a,*s1p2b,*s1p2c,*s1p2d,*s1p2e,
			*s1p30,*s1p31,*s1p32,*s1p33,*s1p34,*s1p35,*s1p36,*s1p37,*s1p38,*s1p39,*s1p3a,*s1p3b,*s1p3c,*s1p3d,*s1p3e,
		  #endif
			*r00,
		  #if !COMPACT_OBJ
			*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09,*r0a,*r0b,*r0c,*r0d,*r0e,
			*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r1a,*r1b,*r1c,*r1d,*r1e,
			*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r2a,*r2b,*r2c,*r2d,*r2e,
			*r30,*r31,*r32,*r33,*r34,*r35,*r36,*r37,*r38,*r39,*r3a,*r3b,*r3c,*r3d,*r3e,
		  #endif
			*x00,*x01,*x02,*x03,*x04,*x05,*x06,*x07,*x08,*x09,*x0a,*x0b,*x0c,*x0d,*x0e,
			*y00,*y01,*y02,*y03,*y04,*y05,*y06,*y07,*y08,*y09,*y0a,*y0b,*y0c,*y0d,*y0e,
			*cy_r	// Need RADIX slots for sse2 carries, RADIX/2 for avx
		  #ifdef USE_AVX
			,*cy_i	// Need RADIX slots for sse2 carries, RADIX/2 for avx
		  #endif
			;
	  #ifdef USE_AVX
		vec_dbl *base_negacyclic_root;
	  #endif

	  #ifndef USE_AVX
		/* These are used in conjunction with the langth-ODD_RADIX arrays in the USE_SCALAR_CARRY #define below;
		In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
		int idx_offset, idx_incr;
	  #endif
		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else

		const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
						s   =  0.86602540378443864675,	/* sin(twopi/3)		*/
						cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
						cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
						ss3 =  0.95105651629515357210,	/*  sin(u) */
						sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
						sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
		double *base, *baseinv, *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr, bs,bsinv;
		int wts_idx_incr;
		const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int jt,jp,m,m2;
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
		struct complex *rn0 = thread_arg->rn0;
		struct complex *rn1 = thread_arg->rn1;

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
		p44 = p40 + p04;
		p48 = p44 + p04;
		p52 = p48 + p04;
		p56 = p52 + p04;

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
		p44 = p44 + ( (p44 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p52 = p52 + ( (p52 >> DAT_BITS) << PAD_BITS );
		p56 = p56 + ( (p56 >> DAT_BITS) << PAD_BITS );
	#if !defined(USE_SSE2) || (COMPACT_OBJ != 0)
		p0123[0] = 0; p0123[1] = p01; p0123[2] = p02; p0123[3] = p03;
	#endif
		poff[0x0] =   0; poff[0x1] = p04; poff[0x2] = p08; poff[0x3] = p12;
		poff[0x4] = p16; poff[0x5] = p20; poff[0x6] = p24; poff[0x7] = p28;
		poff[0x8] = p32; poff[0x9] = p36; poff[0xa] = p40; poff[0xb] = p44;
		poff[0xc] = p48; poff[0xd] = p52; poff[0xe] = p56;

	#ifdef USE_SSE2
		uint32 nwt16 = nwt << L2_SZ_VD;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code

		r00 = thread_arg->r00;tmp = r00 + 0x78;
		r00 = r00 + 0x00;		s1p00 = tmp + 0x00;
	  #if !COMPACT_OBJ
		r01 = r00 + 0x02;		s1p01 = tmp + 0x02;
		r02 = r00 + 0x04;		s1p02 = tmp + 0x04;
		r03 = r00 + 0x06;		s1p03 = tmp + 0x06;
		r04 = r00 + 0x08;		s1p04 = tmp + 0x08;
		r05 = r00 + 0x0a;		s1p05 = tmp + 0x0a;
		r06 = r00 + 0x0c;		s1p06 = tmp + 0x0c;
		r07 = r00 + 0x0e;		s1p07 = tmp + 0x0e;
		r08 = r00 + 0x10;		s1p08 = tmp + 0x10;
		r09 = r00 + 0x12;		s1p09 = tmp + 0x12;
		r0a = r00 + 0x14;		s1p0a = tmp + 0x14;
		r0b = r00 + 0x16;		s1p0b = tmp + 0x16;
		r0c = r00 + 0x18;		s1p0c = tmp + 0x18;
		r0d = r00 + 0x1a;		s1p0d = tmp + 0x1a;
		r0e = r00 + 0x1c;		s1p0e = tmp + 0x1c;
		r10 = r00 + 0x1e;		s1p10 = tmp + 0x1e;
		r11 = r00 + 0x20;		s1p11 = tmp + 0x20;
		r12 = r00 + 0x22;		s1p12 = tmp + 0x22;
		r13 = r00 + 0x24;		s1p13 = tmp + 0x24;
		r14 = r00 + 0x26;		s1p14 = tmp + 0x26;
		r15 = r00 + 0x28;		s1p15 = tmp + 0x28;
		r16 = r00 + 0x2a;		s1p16 = tmp + 0x2a;
		r17 = r00 + 0x2c;		s1p17 = tmp + 0x2c;
		r18 = r00 + 0x2e;		s1p18 = tmp + 0x2e;
		r19 = r00 + 0x30;		s1p19 = tmp + 0x30;
		r1a = r00 + 0x32;		s1p1a = tmp + 0x32;
		r1b = r00 + 0x34;		s1p1b = tmp + 0x34;
		r1c = r00 + 0x36;		s1p1c = tmp + 0x36;
		r1d = r00 + 0x38;		s1p1d = tmp + 0x38;
		r1e = r00 + 0x3a;		s1p1e = tmp + 0x3a;
		r20 = r00 + 0x3c;		s1p20 = tmp + 0x3c;
		r21 = r00 + 0x3e;		s1p21 = tmp + 0x3e;
		r22 = r00 + 0x40;		s1p22 = tmp + 0x40;
		r23 = r00 + 0x42;		s1p23 = tmp + 0x42;
		r24 = r00 + 0x44;		s1p24 = tmp + 0x44;
		r25 = r00 + 0x46;		s1p25 = tmp + 0x46;
		r26 = r00 + 0x48;		s1p26 = tmp + 0x48;
		r27 = r00 + 0x4a;		s1p27 = tmp + 0x4a;
		r28 = r00 + 0x4c;		s1p28 = tmp + 0x4c;
		r29 = r00 + 0x4e;		s1p29 = tmp + 0x4e;
		r2a = r00 + 0x50;		s1p2a = tmp + 0x50;
		r2b = r00 + 0x52;		s1p2b = tmp + 0x52;
		r2c = r00 + 0x54;		s1p2c = tmp + 0x54;
		r2d = r00 + 0x56;		s1p2d = tmp + 0x56;
		r2e = r00 + 0x58;		s1p2e = tmp + 0x58;
		r30 = r00 + 0x5a;		s1p30 = tmp + 0x5a;
		r31 = r00 + 0x5c;		s1p31 = tmp + 0x5c;
		r32 = r00 + 0x5e;		s1p32 = tmp + 0x5e;
		r33 = r00 + 0x60;		s1p33 = tmp + 0x60;
		r34 = r00 + 0x62;		s1p34 = tmp + 0x62;
		r35 = r00 + 0x64;		s1p35 = tmp + 0x64;
		r36 = r00 + 0x66;		s1p36 = tmp + 0x66;
		r37 = r00 + 0x68;		s1p37 = tmp + 0x68;
		r38 = r00 + 0x6a;		s1p38 = tmp + 0x6a;
		r39 = r00 + 0x6c;		s1p39 = tmp + 0x6c;
		r3a = r00 + 0x6e;		s1p3a = tmp + 0x6e;
		r3b = r00 + 0x70;		s1p3b = tmp + 0x70;
		r3c = r00 + 0x72;		s1p3c = tmp + 0x72;
		r3d = r00 + 0x74;		s1p3d = tmp + 0x74;
		r3e = r00 + 0x76;		s1p3e = tmp + 0x76;
	  #endif
		tmp += 0x78;
		x00   = tmp + 0x00;
		x01   = tmp + 0x02;
		x02   = tmp + 0x04;
		x03   = tmp + 0x06;
		x04   = tmp + 0x08;
		x05   = tmp + 0x0a;
		x06   = tmp + 0x0c;
		x07   = tmp + 0x0e;
		x08   = tmp + 0x10;
		x09   = tmp + 0x12;
		x0a   = tmp + 0x14;
		x0b   = tmp + 0x16;
		x0c   = tmp + 0x18;
		x0d   = tmp + 0x1a;
		x0e   = tmp + 0x1c;
		tmp += 0x1e;
		y00   = tmp + 0x00;
		y01   = tmp + 0x02;
		y02   = tmp + 0x04;
		y03   = tmp + 0x06;
		y04   = tmp + 0x08;
		y05   = tmp + 0x0a;
		y06   = tmp + 0x0c;
		y07   = tmp + 0x0e;
		y08   = tmp + 0x10;
		y09   = tmp + 0x12;
		y0a   = tmp + 0x14;
		y0b   = tmp + 0x16;
		y0c   = tmp + 0x18;
		y0d   = tmp + 0x1a;
		y0e   = tmp + 0x1c;
		tmp += 0x1e;

		ASSERT((tmp->d0 == tmp->d1) && (tmp->d0 == c3m1), "thread-local memcheck failed!");
		sse2_c3m1 = tmp + 0x00;
		//sse2_s    = tmp + 0x01;
		sse2_cn1  = tmp + 0x02;
		//sse2_cn2  = tmp + 0x03;
		//sse2_ss3  = tmp + 0x04;
		//sse2_sn1  = tmp + 0x05;
		//sse2_sn2  = tmp + 0x06;
	#ifndef USE_ARM_V8_SIMD // TODO that looks strange that this is needed to get rid of unused 'two'
		two       = tmp + 0x07;
	#endif
		tmp += 0x08;
	   	tm2 = tmp + 0x1e;	/* Use these as handy temps */

	#ifdef USE_AVX512
		cy_r = tmp;	cy_i = tmp+0x8;	tmp += 2*0x8;	//  RADIX/8 (and round up) vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		//sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod
	#elif defined(USE_AVX)
		cy_r = tmp;	cy_i = tmp+0xf;	tmp += 2*0xf;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 340 complex
		// This is where the value of half_arr_offset60 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
								// +20 = 390 complex, round up to nearby multiple of 4
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod
	  #else
		cy_r = tmp;	/* cy_i = tmp+0x1e; */	tmp += 2*0x1e;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 370 complex
		// This is where the value of half_arr_offset60 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*ODD_RADIX] x 16 for Fermat-mod */
								// +20 = 390 complex, round up to nearby multiple of 4
	  #endif

		ASSERT((r00 == thread_arg->r00), "thread-local memcheck failed!");
		ASSERT((half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
	  #ifndef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts:
		ASSERT((sse2_rnd->d0 == crnd && sse2_rnd->d1 == crnd), "thread-local memcheck failed!");
	  #endif
		tmp = half_arr;
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
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
		dtmp = (tmp)->d0 * (tmp+ODD_RADIX)->d0;	ASSERT(fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp)->d1 * (tmp+ODD_RADIX)->d1;	ASSERT(fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
	}

		VEC_DBL_INIT(max_err, 0.0);

		sign_mask = (uint64*)(r00 + radix60_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (#doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
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
		#include "radix60_main_carry_loop.h"

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
