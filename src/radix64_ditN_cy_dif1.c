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

#define RADIX 64	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

#define EPS 1e-10

#ifdef USE_SSE2	// This toggle only supported for SSE2,AVX,AVX2 ... all other SIMD modes are hardcoded == 0:
  #if !defined(USE_ARM_V8_SIMD) && !defined(USE_AVX512)
	#define USE_SCALAR_DFT_MACRO	0	// Seems a tad slower in all of SSE2,AVX,AVX2, so set = 0
  #endif
#else
	#define USE_SCALAR_DFT_MACRO	1
	#warning USE_SCALAR_DFT_MACRO = 1
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
  // For Fermat-mod we use RADIX*4 = 256 [note there is no Fermat-mod LOACC option for this power-of-2 DFT] more
  // slots in AVX mode for the compact negacyclic-roots chained-multiply scheme. Add larger of the 2 numbers -
  // 256 for AVX, 40 for SSE2 - to (half_arr_offset64 + RADIX) to get AVX value of radix64_creals_in_local_store:
  #if USE_SCALAR_DFT_MACRO
	#define OFF	0
  #else
	#define OFF 0x32	// Extra alloc for sincos data
  #endif

  #ifdef USE_AVX512
	const int half_arr_offset64 = 0x100 + OFF + (RADIX>>3);	// RADIX/8 vec_dbl slots for carries in AVX mode
	const int radix64_creals_in_local_store = (0x100 + OFF + (RADIX>>3) + RADIX + 256 + 3)&(~0x3); // May not need this many slots, but just use AVX-alloc for now
  #elif defined(USE_AVX)
	const int half_arr_offset64 = 0x100 + OFF + (RADIX>>2);	// RADIX/4 vec_dbl slots for carries in AVX mode
	// (half_arr_offset64 + RADIX) + 256 and round up to nearest multiple of 8
/*	const int radix64_creals_in_local_store = (half_arr_offset64 + RADIX + 256 + 3)&(~0x3);	<*** GCC [not Clang] gives 'error: initializer element is not constant'... */
	const int radix64_creals_in_local_store = (0x100 + OFF + (RADIX>>2) + RADIX + 256 + 3)&(~0x3); // <*** ...here is that workaround for that stupidity.
  #else
	const int half_arr_offset64 = 0x100 + OFF + (RADIX>>1);	// RADIX/2 vec_dbl slots for carries in SSE2 mode
	// (half_arr_offset64 + RADIX) + 40 and round up to nearest multiple of 8
/*	const int radix64_creals_in_local_store = (half_arr_offset64 + RADIX +  40 + 3)&(~0x3);	<*** GCC [not Clang] gives 'error: initializer element is not constant'... */
	const int radix64_creals_in_local_store = (0x100 + OFF + (RADIX>>1) + RADIX +  40 + 3)&(~0x3); // <*** ...here is that workaround for that stupidity.
  #endif

  #ifdef USE_AVX
	const uint64 radix64_avx_negadwt_consts[RADIX] = {	// 8 entries per line ==> RADIX/8 lines:
		0x3FF0000000000000ull,0x3FEFFD886084CD0Dull,0x3FEFF621E3796D7Eull,0x3FEFE9CDAD01883Aull,0x3FEFD88DA3D12526ull,0x3FEFC26470E19FD3ull,0x3FEFA7557F08A517ull,0x3FEF8764FA714BA9ull,
		0x3FEF6297CFF75CB0ull,0x3FEF38F3AC64E589ull,0x3FEF0A7EFB9230D7ull,0x3FEED740E7684963ull,0x3FEE9F4156C62DDAull,0x3FEE6288EC48E112ull,0x3FEE212104F686E5ull,0x3FEDDB13B6CCC23Cull,
		0x3FED906BCF328D46ull,0x3FED4134D14DC93Aull,0x3FECED7AF43CC773ull,0x3FEC954B213411F5ull,0x3FEC38B2F180BDB1ull,0x3FEBD7C0AC6F952Aull,0x3FEB728345196E3Eull,0x3FEB090A58150200ull,
		0x3FEA9B66290EA1A3ull,0x3FEA29A7A0462782ull,0x3FE9B3E047F38741ull,0x3FE93A22499263FBull,0x3FE8BC806B151741ull,0x3FE83B0E0BFF976Eull,0x3FE7B5DF226AAFAFull,0x3FE72D0837EFFF96ull,
		0x3FE6A09E667F3BCDull,0x3FE610B7551D2CDFull,0x3FE57D69348CECA0ull,0x3FE4E6CABBE3E5E9ull,0x3FE44CF325091DD6ull,0x3FE3AFFA292050B9ull,0x3FE30FF7FCE17035ull,0x3FE26D054CDD12DFull,
		0x3FE1C73B39AE68C8ull,0x3FE11EB3541B4B23ull,0x3FE073879922FFEEull,0x3FDF8BA4DBF89ABAull,0x3FDE2B5D3806F63Bull,0x3FDCC66E9931C45Eull,0x3FDB5D1009E15CC0ull,0x3FD9EF7943A8ED8Aull,
		0x3FD87DE2A6AEA963ull,0x3FD7088530FA459Full,0x3FD58F9A75AB1FDDull,0x3FD4135C94176601ull,0x3FD294062ED59F06ull,0x3FD111D262B1F677ull,0x3FCF19F97B215F1Bull,0x3FCC0B826A7E4F63ull,
		0x3FC8F8B83C69A60Bull,0x3FC5E214448B3FC6ull,0x3FC2C8106E8E613Aull,0x3FBF564E56A9730Eull,0x3FB917A6BC29B42Cull,0x3FB2D52092CE19F6ull,0x3FA91F65F10DD814ull,0x3F992155F7A3667Eull
	};
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

int radix64_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-64 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-64 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix16_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix64_ditN_cy_dif1";
#if USE_SCALAR_DFT_MACRO && defined(USE_SSE2)
	static int thr_id = 0;	// Master thread gets this special id
#endif
#if USE_SCALAR_DFT_MACRO && defined(USE_SSE2) && !defined(MULTITHREAD)
	static int dft_offsets[RADIX], c_offsets[RADIX];
#endif
#if !defined(MULTITHREAD) && defined(USE_SSE2)
	const int pfetch_dist = PFETCH_DIST;
#endif
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
  #ifdef USE_AVX512
	const int jhi_wrap_mers = 15;
	const int jhi_wrap_ferm = 15;
  #else
	const int jhi_wrap_mers =  7;
	const int jhi_wrap_ferm = 15;	// For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
  #endif
	int NDIVR,i,j,j1,jt,full_pass,khi,l,outer;
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	int j2,jp;
  #endif
  #ifndef MULTITHREAD
	int jstart,jhi;
  #endif
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	int ntmp;
  #endif
  #ifdef USE_SSE2
	int nbytes;
  #endif
#if !defined(MULTITHREAD) && defined(USE_SSE2) && !defined(USE_AVX)
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 4|8|16 for avx512,avx,sse, respectively:
	int incr;
	const int incr_long = 16,incr_med = 8,incr_short = 4;
	// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
	if(USE_SHORT_CY_CHAIN == 0)
		incr = incr_long;
	else if(USE_SHORT_CY_CHAIN == 1)
		incr = incr_med;
	else
		incr = incr_short;
#endif
#ifndef MULTITHREAD
	int k1;
  #if !defined(USE_SSE2) || defined(USE_AVX)
	int k2;
  #endif
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
	static double wts_mult[2], inv_mult[2];	// Const wts-multiplier and 2*(its multiplicative inverse)
#if !defined(MULTITHREAD) && (!defined(USE_SSE2) || defined(USE_AVX))
	double rt,it, wt_re,wt_im;	// Fermat-mod weights stuff, used in both scalar and AVX mode
#endif
#if !defined(MULTITHREAD) && !defined(USE_SSE2)
	double wi_re,wi_im;	// Fermat-mod weights stuff, used in both scalar and AVX mode
#endif
	static uint32 bjmodnini;
	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30,p38, nsave = 0;
  #if !defined(MULTITHREAD) || defined(USE_SSE2)
	static uint32 nm1;
  #endif
	static int poff[RADIX>>2];	// Store mults of p-offsets for loop-controlled DFT macro calls
  #ifndef MULTITHREAD
   #if defined(USE_SSE2) && !USE_SCALAR_DFT_MACRO
	static int po_lin[8];	// Store mults of p-offsets for loop-controlled DFT macro calls
   #endif
   #ifndef USE_SSE2
	static int po_br[8];	// Store mults of p-offsets for loop-controlled DFT macro calls
   #endif
  #endif
	static double radix_inv, n2inv;
	double scale, dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex t[RADIX];
#if !defined(MULTITHREAD) && !defined(USE_SSE2)
	struct complex *tptr;
#endif
#ifndef MULTITHREAD
	int *itmp;	// Pointer into the bjmodn array
  #if defined(USE_AVX) && !defined(USE_AVX512)
	int *itm2;	// Pointer into the bjmodn array
  #endif
#endif
	int err;
	static int first_entry=TRUE;

	int n_div_nwt;

#ifdef USE_SSE2

	// Uint64 bitmaps for alternate "rounded the other way" copies of sqrt2,isrt2. Default round-to-nearest versions
	// (SQRT2, ISRT2) end in ...3BCD. Since we round these down as ...3BCC90... --> ..3BCC, append _dn to varnames:
	const uint64 sqrt2_dn = 0x3FF6A09E667F3BCCull, isrt2_dn = 0x3FE6A09E667F3BCCull;
  #if !defined(MULTITHREAD) && !defined(USE_AVX)
	int idx_offset,idx_incr;
  #endif
	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_nm1;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else
	double *addr, *add0,*add1,*add2,*add3,*add4,*add5,*add6,*add7;
  #endif	// MULTITHREAD

  #ifndef MULTITHREAD
	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
  #endif
  #ifndef USE_AVX512
	const double crnd = 3.0*0x4000000*0x2000000;
  #endif
  #ifndef USE_AVX
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
  #endif
	static vec_dbl *max_err, *sse2_rnd, *half_arr,
	#if !USE_SCALAR_DFT_MACRO
		*cc0, *ss0, *two,*one,*sqrt2,
		 *isrt2, *cc1, *ss1, *cc2, *ss2, *cc3, *ss3, *cc4, *ss4, *cc5, *ss5, *cc6, *ss6, *cc7, *ss7,
		*nisrt2,*ncc1,*nss1,*ncc2,*nss2,*ncc3,*nss3,*ncc4,*nss4,*ncc5,*nss5,*ncc6,*nss6,*ncc7,*nss7,	// each non-unity root now needs a negated counterpart
	  #ifndef MULTITHREAD
		// Array of 7*14 SIMD twiddles-pointers for the reduced-#args version of SSE2_RADIX8_DI[F,T]_TWIDDLE_OOP
		*w[98], **twid_ptrs,
	  #endif
	#endif // !USE_SCALAR_DFT_MACRO
		// This routine dates from the benighted early days when I treated r's as pointers-to-real-SIMD...
		*r00,
	  #ifndef MULTITHREAD
		*r02,*r04,*r06,*r08,*r0A,*r0C,*r0E,*r10,*r20,*r30,*r40,*r50,*r60,*r70,
		//*r12,*r14,*r16,*r18,*r1A,*r1C,*r1E,
		//*r22,*r24,*r26,*r28,*r2A,*r2C,*r2E,*r32,*r34,*r36,*r38,*r3A,*r3C,*r3E,
		//*r42,*r44,*r46,*r48,*r4A,*r4C,*r4E,*r52,*r54,*r56,*r58,*r5A,*r5C,*r5E,
		//*r62,*r64,*r66,*r68,*r6A,*r6C,*r6E,*r72,*r74,*r76,*r78,*r7A,*r7C,*r7E,
	  #endif
	  #ifndef MULTITHREAD
		// ...and s's as pointers-to-complex-SIMD; thus the r-indices run 2x faster than the s-ones:
		*s1p00,*s1p01,*s1p02,*s1p03,*s1p04,*s1p05,*s1p06,*s1p07,*s1p08,*s1p10,*s1p18,*s1p20,*s1p28,*s1p30,*s1p38,
		//*s1p09,*s1p0a,*s1p0b,*s1p0c,*s1p0d,*s1p0e,*s1p0f,
		//*s1p11,*s1p12,*s1p13,*s1p14,*s1p15,*s1p16,*s1p17,*s1p19,*s1p1a,*s1p1b,*s1p1c,*s1p1d,*s1p1e,*s1p1f,
		//*s1p21,*s1p22,*s1p23,*s1p24,*s1p25,*s1p26,*s1p27,*s1p29,*s1p2a,*s1p2b,*s1p2c,*s1p2d,*s1p2e,*s1p2f,
		//*s1p31,*s1p32,*s1p33,*s1p34,*s1p35,*s1p36,*s1p37,*s1p39,*s1p3a,*s1p3b,*s1p3c,*s1p3d,*s1p3e,*s1p3f,
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
  #ifndef MULTITHREAD
   #ifdef USE_AVX
	vec_dbl *tm0;	// Non-static utility ptrs
   #endif
	vec_dbl *tm1;	// Non-static utility ptrs
  #endif

#else
  #ifndef MULTITHREAD
	static int p0123[4];
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
	static task_control_t   task_control = {NULL, (void*)cy64_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// In scalar-double mode use a const-double array of 7 sets of 14 doubles:
	const double DFT64_TWIDDLES[7][14] = {
		{ 0.0,1.0,ISRT2,ISRT2,-ISRT2,ISRT2,c16,s16,-s16,c16,s16,c16,-c16,s16 },
		{ ISRT2,ISRT2,c16,s16,s16,c16,c32_1,s32_1,s32_3,c32_3,c32_3,s32_3,s32_1,c32_1 },
		{ -ISRT2,ISRT2,s16,c16,-c16,-s16,c32_3,s32_3,-c32_1,s32_1,-s32_1,c32_1,-s32_3,-c32_3 },
		{ c16,s16,c32_1,s32_1,c32_3,s32_3,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7 },
		{ -s16,c16,s32_3,c32_3,-c32_1,s32_1,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3 },
		{ s16,c16,c32_3,s32_3,-s32_1,c32_1,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5 },
		{ -c16,s16,s32_1,c32_1,-s32_3,-c32_3,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1 }
	};

	// Vars needed in scalar mode only:
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	int m,m2;
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	double *addr,*addi;
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

	  #if !defined(MULTITHREAD) || defined(USE_SSE2)
		nm1   = n-1;
	  #endif

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
		// consisting of 128 vec_dbl and ([8 if SSE2, 16 if AVX] + RADIX/2) uint64 element slots per thread
		cslots_in_local_store = radix64_creals_in_local_store + (20+RADIX/2)/2;	// Just add enough int64 space for both cases, plus some
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix64_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 64 vec_ddl-sized slots of sc_arr for temporaries, next 7 for the nontrivial complex 16th roots,
	next 32 for the vector carries, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = sc_ptr;
		r00 = tmp + 0x00;
	  #ifndef MULTITHREAD
		r02 = tmp + 0x02;
		r04 = tmp + 0x04;
		r06 = tmp + 0x06;
		r08 = tmp + 0x08;
		r0A = tmp + 0x0a;
		r0C = tmp + 0x0c;
		r0E = tmp + 0x0e;
		r10 = tmp + 0x10;
		//r12 = tmp + 0x12;
		//r14 = tmp + 0x14;
		//r16 = tmp + 0x16;
		//r18 = tmp + 0x18;
		//r1A = tmp + 0x1a;
		//r1C = tmp + 0x1c;
		//r1E = tmp + 0x1e;
		r20 = tmp + 0x20;
		//r22 = tmp + 0x22;
		//r24 = tmp + 0x24;
		//r26 = tmp + 0x26;
		//r28 = tmp + 0x28;
		//r2A = tmp + 0x2a;
		//r2C = tmp + 0x2c;
		//r2E = tmp + 0x2e;
		r30 = tmp + 0x30;
		//r32 = tmp + 0x32;
		//r34 = tmp + 0x34;
		//r36 = tmp + 0x36;
		//r38 = tmp + 0x38;
		//r3A = tmp + 0x3a;
		//r3C = tmp + 0x3c;
		//r3E = tmp + 0x3e;
		r40 = tmp + 0x40;
		//r42 = tmp + 0x42;
		//r44 = tmp + 0x44;
		//r46 = tmp + 0x46;
		//r48 = tmp + 0x48;
		//r4A = tmp + 0x4a;
		//r4C = tmp + 0x4c;
		//r4E = tmp + 0x4e;
		r50 = tmp + 0x50;
		//r52 = tmp + 0x52;
		//r54 = tmp + 0x54;
		//r56 = tmp + 0x56;
		//r58 = tmp + 0x58;
		//r5A = tmp + 0x5a;
		//r5C = tmp + 0x5c;
		//r5E = tmp + 0x5e;
		r60 = tmp + 0x60;
		//r62 = tmp + 0x62;
		//r64 = tmp + 0x64;
		//r66 = tmp + 0x66;
		//r68 = tmp + 0x68;
		//r6A = tmp + 0x6a;
		//r6C = tmp + 0x6c;
		//r6E = tmp + 0x6e;
		r70 = tmp + 0x70;
		//r72 = tmp + 0x72;
		//r74 = tmp + 0x74;
		//r76 = tmp + 0x76;
		//r78 = tmp + 0x78;
		//r7A = tmp + 0x7a;
		//r7C = tmp + 0x7c;
		//r7E = tmp + 0x7e;
	  #endif
		tmp += 0x80;
	  #ifndef MULTITHREAD
		s1p00 = tmp + 0x00;
		s1p01 = tmp + 0x02;
		s1p02 = tmp + 0x04;
		s1p03 = tmp + 0x06;
		s1p04 = tmp + 0x08;
		s1p05 = tmp + 0x0a;
		s1p06 = tmp + 0x0c;
		s1p07 = tmp + 0x0e;
		s1p08 = tmp + 0x10;
		//s1p09 = tmp + 0x12;
		//s1p0a = tmp + 0x14;
		//s1p0b = tmp + 0x16;
		//s1p0c = tmp + 0x18;
		//s1p0d = tmp + 0x1a;
		//s1p0e = tmp + 0x1c;
		//s1p0f = tmp + 0x1e;
		s1p10 = tmp + 0x20;
		//s1p11 = tmp + 0x22;
		//s1p12 = tmp + 0x24;
		//s1p13 = tmp + 0x26;
		//s1p14 = tmp + 0x28;
		//s1p15 = tmp + 0x2a;
		//s1p16 = tmp + 0x2c;
		//s1p17 = tmp + 0x2e;
		s1p18 = tmp + 0x30;
		//s1p19 = tmp + 0x32;
		//s1p1a = tmp + 0x34;
		//s1p1b = tmp + 0x36;
		//s1p1c = tmp + 0x38;
		//s1p1d = tmp + 0x3a;
		//s1p1e = tmp + 0x3c;
		//s1p1f = tmp + 0x3e;
		s1p20 = tmp + 0x40;
		//s1p21 = tmp + 0x42;
		//s1p22 = tmp + 0x44;
		//s1p23 = tmp + 0x46;
		//s1p24 = tmp + 0x48;
		//s1p25 = tmp + 0x4a;
		//s1p26 = tmp + 0x4c;
		//s1p27 = tmp + 0x4e;
		s1p28 = tmp + 0x50;
		//s1p29 = tmp + 0x52;
		//s1p2a = tmp + 0x54;
		//s1p2b = tmp + 0x56;
		//s1p2c = tmp + 0x58;
		//s1p2d = tmp + 0x5a;
		//s1p2e = tmp + 0x5c;
		//s1p2f = tmp + 0x5e;
		s1p30 = tmp + 0x60;
		//s1p31 = tmp + 0x62;
		//s1p32 = tmp + 0x64;
		//s1p33 = tmp + 0x66;
		//s1p34 = tmp + 0x68;
		//s1p35 = tmp + 0x6a;
		//s1p36 = tmp + 0x6c;
		//s1p37 = tmp + 0x6e;
		s1p38 = tmp + 0x70;
		//s1p39 = tmp + 0x72;
		//s1p3a = tmp + 0x74;
		//s1p3b = tmp + 0x76;
		//s1p3c = tmp + 0x78;
		//s1p3d = tmp + 0x7a;
		//s1p3e = tmp + 0x7c;
		//s1p3f = tmp + 0x7e;
	  #endif
		tmp += 0x80;
	  #if !USE_SCALAR_DFT_MACRO
		two     = tmp + 0;	// AVX+ versions of various DFT macros need consts [2,1,sqrt2,isrt2] quartet laid out thusly
		one     = tmp + 1;
		sqrt2	= tmp + 2;
	//	isrt2   = tmp + 3;	Unnamed slot, since previous layout below already has an iart2 pointer
		tmp += 4;
		// Each non-unity root now needs a negated counterpart:
		/* Stupidity: Since a truly general-purpose [in the sense that it can be used for our radix-128 internal-twiddles]
		radix-8 DFT-with-twiddles macro needs 8 in-addresses [corr. to the 8 real parts of the input data], 8 o-addresses,
		and 7 each of cosine and sine data [which cannot be assumed to occur in fixed-stride pairs - cf. our usage of
		SSE2_RADIX8_DIT_TWIDDLE_OOP() below], that hits the GCC hard limit of 30-operands for ASM macros, but we still
		need one more operand for the ISRT2 pointer. Only easy workaround I found for this is to stick a vector-ISRT2 copy
		in between each +-[cc,ss] vector-data pair, thus any time we need a vector-isrt2 for the radix-8 internal twiddles
		we get it at (vec_dbl*)cc-1.
		/Stupidity */
		nisrt2	= tmp + 0x00;	// For the +- isrt2 pair put the - datum first, thus cc0 satisfies
		 isrt2	= tmp + 0x01;	// the same "cc-1 gets you isrt2" property as do the other +-[cc,ss] pairs.
		 cc0	= tmp + 0x02;
		 ss0	= tmp + 0x03;
	// [copy isrt2]	= tmp + 0x04;
		 cc1	= tmp + 0x05;
		 ss1	= tmp + 0x06;
	// [copy isrt2]	= tmp + 0x07;
		ncc1	= tmp + 0x08;
		nss1	= tmp + 0x09;
	// [copy isrt2]	= tmp + 0x0a;
		 cc2	= tmp + 0x0b;
		 ss2	= tmp + 0x0c;
	// [copy isrt2]	= tmp + 0x0d;
		ncc2	= tmp + 0x0e;
		nss2	= tmp + 0x0f;
	// [copy isrt2]	= tmp + 0x10;
		 cc3	= tmp + 0x11;
		 ss3	= tmp + 0x12;
	// [copy isrt2]	= tmp + 0x13;
		ncc3	= tmp + 0x14;
		nss3	= tmp + 0x15;
	// [copy isrt2]	= tmp + 0x16;
		 cc4	= tmp + 0x17;
		 ss4	= tmp + 0x18;
	// [copy isrt2]	= tmp + 0x19;
		ncc4	= tmp + 0x1a;
		nss4	= tmp + 0x1b;
	// [copy isrt2]	= tmp + 0x1c;
		 cc5	= tmp + 0x1d;
		 ss5	= tmp + 0x1e;
	// [copy isrt2]	= tmp + 0x1f;
		ncc5	= tmp + 0x20;
		nss5	= tmp + 0x21;
	// [copy isrt2]	= tmp + 0x22;
		 cc6	= tmp + 0x23;
		 ss6	= tmp + 0x24;
	// [copy isrt2]	= tmp + 0x25;
		ncc6	= tmp + 0x26;
		nss6	= tmp + 0x27;
	// [copy isrt2]	= tmp + 0x28;
		 cc7	= tmp + 0x29;
		 ss7	= tmp + 0x2a;
	// [copy isrt2]	= tmp + 0x2b;
		ncc7	= tmp + 0x2c;
		nss7	= tmp + 0x2d;
		tmp += 0x2e;	// sc_ptr += 0x132
	   #ifndef MULTITHREAD
		// Init elements of array of 7*14 SIMD twiddles-pointers for the reduced-#args version of SSE2_RADIX8_DI[F,T]_TWIDDLE_OOP
		i = 0;
		w[i++]=ss0;w[i++]=cc0;w[i++]=isrt2;w[i++]=isrt2;w[i++]=nisrt2;w[i++]=isrt2;w[i++]=cc4;w[i++]=ss4;w[i++]=nss4;w[i++]=cc4;w[i++]=ss4;w[i++]=cc4;w[i++]=ncc4;w[i++]=ss4;
		w[i++]=isrt2;w[i++]=isrt2;w[i++]=cc4;w[i++]=ss4;w[i++]=ss4;w[i++]=cc4;w[i++]=cc2;w[i++]=ss2;w[i++]=ss6;w[i++]=cc6;w[i++]=cc6;w[i++]=ss6;w[i++]=ss2;w[i++]=cc2;
		w[i++]=nisrt2;w[i++]=isrt2;w[i++]=ss4;w[i++]=cc4;w[i++]=ncc4;w[i++]=nss4;w[i++]=cc6;w[i++]=ss6;w[i++]=ncc2;w[i++]=ss2;w[i++]=nss2;w[i++]=cc2;w[i++]=nss6;w[i++]=ncc6;
		w[i++]=cc4;w[i++]=ss4;w[i++]=cc2;w[i++]=ss2;w[i++]=cc6;w[i++]=ss6;w[i++]=cc1;w[i++]=ss1;w[i++]=cc5;w[i++]=ss5;w[i++]=cc3;w[i++]=ss3;w[i++]=cc7;w[i++]=ss7;
		w[i++]=nss4;w[i++]=cc4;w[i++]=ss6;w[i++]=cc6;w[i++]=ncc2;w[i++]=ss2;w[i++]=cc5;w[i++]=ss5;w[i++]=ncc7;w[i++]=ss7;w[i++]=ss1;w[i++]=cc1;w[i++]=ncc3;w[i++]=nss3;
		w[i++]=ss4;w[i++]=cc4;w[i++]=cc6;w[i++]=ss6;w[i++]=nss2;w[i++]=cc2;w[i++]=cc3;w[i++]=ss3;w[i++]=ss1;w[i++]=cc1;w[i++]=ss7;w[i++]=cc7;w[i++]=nss5;w[i++]=cc5;
		w[i++]=ncc4;w[i++]=ss4;w[i++]=ss2;w[i++]=cc2;w[i++]=nss6;w[i++]=ncc6;w[i++]=cc7;w[i++]=ss7;w[i++]=ncc3;w[i++]=nss3;w[i++]=nss5;w[i++]=cc5;w[i++]=ss1;w[i++]=ncc1;
	   #endif
	  #endif // !USE_SCALAR_DFT_MACRO
	  #ifdef USE_AVX512
		cy_r = tmp;										// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
	   #ifndef MULTITHREAD
					cy_i = tmp+0x08;					// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
	   #endif
										tmp += 2*0x08;	// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
		cy_r = tmp;										// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
	   #ifndef MULTITHREAD
					cy_i = tmp+0x10;					// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
	   #endif
										tmp += 2*0x10;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 340 vec_dbl
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 96 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy_r = tmp;	/* cy_i = tmp+0x20; */	tmp += 2*0x20;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 372 complex
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 32 x 16 bytes for Mersenne-mod, 2 for Fermat-mod */
	  #endif
//		ASSERT(half_arr_offset == (uint32)(half_arr-sc_ptr), "half_arr_offset mismatches actual!");
		ASSERT((radix64_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r00) + (20 << L2_SZ_VD), "radix64_creals_in_local_store checksum failed!");

	  #if !USE_SCALAR_DFT_MACRO
		/* These remain fixed: */
		VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );
		tmp = sqrt2+1;
		// 2 unnamed slots for alternate "rounded the other way" copies of sqrt2,isrt2:
		dtmp = *(double *)&sqrt2_dn;	VEC_DBL_INIT(sqrt2, dtmp);
		dtmp = *(double *)&isrt2_dn;	VEC_DBL_INIT(isrt2, dtmp);
		VEC_DBL_INIT(nisrt2,-dtmp);
		VEC_DBL_INIT( isrt2, dtmp);									// Copies of +ISRT2 needed for 30-asm-macro-operand-GCC-limit workaround:
		VEC_DBL_INIT( cc0,   1.0);		VEC_DBL_INIT( ss0,   0.0);	//	tmp =  cc0-1; ASSERT(tmp->d0 == ISRT2 && tmp->d1 == ISRT2, "tmp->d0,1 != ISRT2");	Disable to allow "round down" variant
		VEC_DBL_INIT( cc1, c64_1);		VEC_DBL_INIT( ss1, s64_1);		tmp =  cc1-1; VEC_DBL_INIT(tmp, dtmp);
		VEC_DBL_INIT( cc2, c32_1);		VEC_DBL_INIT( ss2, s32_1);		tmp =  cc2-1; VEC_DBL_INIT(tmp, dtmp);
		VEC_DBL_INIT( cc3, c64_3);		VEC_DBL_INIT( ss3, s64_3);		tmp =  cc3-1; VEC_DBL_INIT(tmp, dtmp);
		VEC_DBL_INIT( cc4, c16  );		VEC_DBL_INIT( ss4, s16  );		tmp =  cc4-1; VEC_DBL_INIT(tmp, dtmp);
		VEC_DBL_INIT( cc5, c64_5);		VEC_DBL_INIT( ss5, s64_5);		tmp =  cc5-1; VEC_DBL_INIT(tmp, dtmp);
		VEC_DBL_INIT( cc6, c32_3);		VEC_DBL_INIT( ss6, s32_3);		tmp =  cc6-1; VEC_DBL_INIT(tmp, dtmp);
		VEC_DBL_INIT( cc7, c64_7);		VEC_DBL_INIT( ss7, s64_7);		tmp =  cc7-1; VEC_DBL_INIT(tmp, dtmp);
		VEC_DBL_INIT(ncc1,-c64_1);		VEC_DBL_INIT(nss1,-s64_1);		tmp = ncc1-1; VEC_DBL_INIT(tmp, dtmp);
		VEC_DBL_INIT(ncc2,-c32_1);		VEC_DBL_INIT(nss2,-s32_1);		tmp = ncc2-1; VEC_DBL_INIT(tmp, dtmp);
		VEC_DBL_INIT(ncc3,-c64_3);		VEC_DBL_INIT(nss3,-s64_3);		tmp = ncc3-1; VEC_DBL_INIT(tmp, dtmp);
		VEC_DBL_INIT(ncc4,-c16  );		VEC_DBL_INIT(nss4,-s16  );		tmp = ncc4-1; VEC_DBL_INIT(tmp, dtmp);
		VEC_DBL_INIT(ncc5,-c64_5);		VEC_DBL_INIT(nss5,-s64_5);		tmp = ncc5-1; VEC_DBL_INIT(tmp, dtmp);
		VEC_DBL_INIT(ncc6,-c32_3);		VEC_DBL_INIT(nss6,-s32_3);		tmp = ncc6-1; VEC_DBL_INIT(tmp, dtmp);
		VEC_DBL_INIT(ncc7,-c64_7);		VEC_DBL_INIT(nss7,-s64_7);		tmp = ncc7-1; VEC_DBL_INIT(tmp, dtmp);
	  #else
	// Init-mode calls to these functions which maintain an internal local-alloc static store:
		thr_id = -1;	// Use this special thread id for any macro-required thread-local-data inits...
		SSE2_RADIX_64_DIF(CY_THREADS,thr_id, 0,0,0,0,0,0);
		SSE2_RADIX_64_DIT(CY_THREADS,thr_id, 0,0,0,0,0);
		thr_id = 0;	// ...then revert to 0.
	  #endif

		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		sse2_rnd->d0 = base[0]; sse2_rnd->d1 = baseinv[1]; sse2_rnd->d2 = wts_mult[1]; sse2_rnd->d3 = inv_mult[0];
	  #else
		VEC_DBL_INIT(sse2_rnd, crnd);
	  #endif

	  #if !USE_SCALAR_DFT_MACRO
		// Propagate the above consts to the remaining threads:
		nbytes = (intptr_t)cy_r - (intptr_t)two;	// #bytes in 1st of above block of consts
		tmp = two;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}
	  #endif
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
			The pattern of the negacyclic-DWT-weights ("nDWTs") applied to the RADIX complex outputs of the final-DIT-pass is like so:
			The nDWTs multiplying each set of RADIX DIT DFT outputs are simply the product of a single complex-root "base multiplier" rbase
			(separately computed for each batch of DFT outputs), which "base root" multiplies the (0 - RADIX-1)st [4*RADIX]th roots of unity,
			i.e.
				 rbase * (j*I*Pi/2)/RADIX, for j = 0, ..., RADIX-1 .

			See the radix28 version of this routine for additional details.
			*/
			#if 0
			// Simple qfloat-based loop to crank out the roots which populate the radix64_avx_negadwt_consts table:
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
			tmp64 = radix64_avx_negadwt_consts[1];	tmp->d1 = tm2->d7 = *(double *)&tmp64;
			tmp64 = radix64_avx_negadwt_consts[2];	tmp->d2 = tm2->d6 = *(double *)&tmp64;
			tmp64 = radix64_avx_negadwt_consts[3];	tmp->d3 = tm2->d5 = *(double *)&tmp64;
			tmp64 = radix64_avx_negadwt_consts[4];	tmp->d4 = tm2->d4 = *(double *)&tmp64;
			tmp64 = radix64_avx_negadwt_consts[5];	tmp->d5 = tm2->d3 = *(double *)&tmp64;
			tmp64 = radix64_avx_negadwt_consts[6];	tmp->d6 = tm2->d2 = *(double *)&tmp64;
			tmp64 = radix64_avx_negadwt_consts[7];	tmp->d7 = tm2->d1 = *(double *)&tmp64;	tmp += 2;
			for(j = 8; j < RADIX; j += 8) {
				tmp64 = radix64_avx_negadwt_consts[j+0];	tmp->d0 = tm2->d0 = *(double *)&tmp64;	tm2 -= 2;
				tmp64 = radix64_avx_negadwt_consts[j+1];	tmp->d1 = tm2->d7 = *(double *)&tmp64;
				tmp64 = radix64_avx_negadwt_consts[j+2];	tmp->d2 = tm2->d6 = *(double *)&tmp64;
				tmp64 = radix64_avx_negadwt_consts[j+3];	tmp->d3 = tm2->d5 = *(double *)&tmp64;
				tmp64 = radix64_avx_negadwt_consts[j+4];	tmp->d4 = tm2->d4 = *(double *)&tmp64;
				tmp64 = radix64_avx_negadwt_consts[j+5];	tmp->d5 = tm2->d3 = *(double *)&tmp64;
				tmp64 = radix64_avx_negadwt_consts[j+6];	tmp->d6 = tm2->d2 = *(double *)&tmp64;
				tmp64 = radix64_avx_negadwt_consts[j+7];	tmp->d7 = tm2->d1 = *(double *)&tmp64;	tmp += 2;
			}
			tmp = base_negacyclic_root + RADIX*2;	// reset to point to start of above block

		  #else

			tmp = base_negacyclic_root + RADIX*2;	// First 2*RADIX slots reserved for RADIX/4 copies of the Re/Im parts of the 4 base multipliers
			tm2 = tmp + RADIX/2 - 1;
			// First elt-pair needs special handling - have the 1.0 in avx_negadwt_consts[0] but the sine term buggers things
			tmp->d0 = 1.0;	(tmp+1)->d0 = 0.0;
			tmp64 = radix64_avx_negadwt_consts[1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;				/* cos(  1*I*Pi/(2*RADIX)) = sin((RADIX-  1)*I*Pi/(2*RADIX)) */
			tmp64 = radix64_avx_negadwt_consts[2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;				/* cos(  2*I*Pi/(2*RADIX)) = sin((RADIX-  2)*I*Pi/(2*RADIX)) */
			tmp64 = radix64_avx_negadwt_consts[3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	tmp += 2;	/* cos(  3*I*Pi/(2*RADIX)) = sin((RADIX-  3)*I*Pi/(2*RADIX)) */
			for(j = 4; j < RADIX; j += 4) {
				tmp64 = radix64_avx_negadwt_consts[j+0];	tmp->d0 = tm2->d0 = *(double *)&tmp64;	tm2 -= 2;
				tmp64 = radix64_avx_negadwt_consts[j+1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;
				tmp64 = radix64_avx_negadwt_consts[j+2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;
				tmp64 = radix64_avx_negadwt_consts[j+3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	tmp += 2;
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
		#ifndef MULTITHREAD
		n_minus_sil   = (struct uint32x16*)sse_nm1 + 1;
		n_minus_silp1 = (struct uint32x16*)sse_nm1 + 2;
		sinwt         = (struct uint32x16*)sse_nm1 + 3;
		sinwtm1       = (struct uint32x16*)sse_nm1 + 4;
		#endif
		nbytes += 256;
	   #else
		#ifndef MULTITHREAD
		n_minus_sil   = (struct uint32x8 *)sse_nm1 + 1;
		n_minus_silp1 = (struct uint32x8 *)sse_nm1 + 2;
		sinwt         = (struct uint32x8 *)sse_nm1 + 3;
		sinwtm1       = (struct uint32x8 *)sse_nm1 + 4;
		#endif
		nbytes += 128;
	   #endif
	  #elif defined(USE_AVX)
	   #ifndef MULTITHREAD
		n_minus_sil   = (struct uint32x4 *)sse_nm1 + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_nm1 + 2;
		sinwt         = (struct uint32x4 *)sse_nm1 + 3;
		sinwtm1       = (struct uint32x4 *)sse_nm1 + 4;
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
	   #ifdef USE_AVX
		bjmodn = (int*)(sinwtm1 + RE_IM_STRIDE);
	   #else
		bjmodn = (int*)(sse_nm1 + RE_IM_STRIDE);
	   #endif
	  #endif

	#endif	// USE_SSE2

		/*   constant index offsets for load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		p01 = NDIVR;
		p02 = p01 +p01;
		p03 = p02 +p01;
		p04 = p03 +p01;
		p05 = p04 +p01;
		p06 = p05 +p01;
		p07 = p06 +p01;
		p08 = p07 +p01;
		p10 = p08 +p08;
		p18 = p10 +p08;
		p20 = p18 +p08;
		p28 = p20 +p08;
		p30 = p28 +p08;
		p38 = p30 +p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p30 = p30 + ( (p30 >> DAT_BITS) << PAD_BITS );
		p38 = p38 + ( (p38 >> DAT_BITS) << PAD_BITS );
	#if !defined(USE_SSE2) && !defined(MULTITHREAD)
		p0123[0] = 0; p0123[1] = p01; p0123[2] = p02; p0123[3] = p03;
	#endif
		poff[0x0] =   0; poff[0x1] = p04    ; poff[0x2] = p08; poff[0x3] = p04+p08;
		poff[0x4] = p10; poff[0x5] = p04+p10; poff[0x6] = p18; poff[0x7] = p04+p18;
		poff[0x8] = p20; poff[0x9] = p04+p20; poff[0xa] = p28; poff[0xb] = p04+p28;
		poff[0xc] = p30; poff[0xd] = p04+p30; poff[0xe] = p38; poff[0xf] = p04+p38;
	#if defined(USE_SSE2) && !USE_SCALAR_DFT_MACRO && !defined(MULTITHREAD)
		// _lin = linear p-multiples (but padding means we can't assume e.g. p02 = 2*p01):
		po_lin[0] =   0; po_lin[1] = p01; po_lin[2] = p02; po_lin[3] = p03;
		po_lin[4] = p04; po_lin[5] = p05; po_lin[6] = p06; po_lin[7] = p07;
	#endif
	#if !defined(USE_SSE2) && !defined(MULTITHREAD)
		// _br  = bit-reversed p-multiples:
		po_br[0] =   0; po_br[1] = p04; po_br[2] = p02; po_br[3] = p06;
		po_br[4] = p01; po_br[5] = p05; po_br[6] = p03; po_br[7] = p07;
	#endif

	#if USE_SCALAR_DFT_MACRO && defined(USE_SSE2) && !defined(MULTITHREAD)
		for(l = 0; l < RADIX; l++) {
			c_offsets[l] = (l<<1);
		}
		// Set array offsets for in/outputs - low parts in fisrt 16 slots, high parts in next 16:
		dft_offsets[0x0] =   0; dft_offsets[0x1] = p01; dft_offsets[0x2] = p02; dft_offsets[0x3] = p03; dft_offsets[0x4] = p04; dft_offsets[0x5] = p05; dft_offsets[0x6] = p06; dft_offsets[0x7] = p07;
		// Now add above low parts to the 7 nonzero high parts:
		l = 0x08;	jp = p08;
		dft_offsets[l+0x0] = jp+  0; dft_offsets[l+0x1] = jp+p01; dft_offsets[l+0x2] = jp+p02; dft_offsets[l+0x3] = jp+p03; dft_offsets[l+0x4] = jp+p04; dft_offsets[l+0x5] = jp+p05; dft_offsets[l+0x6] = jp+p06; dft_offsets[l+0x7] = jp+p07;
		l = 0x10;	jp = p10;
		dft_offsets[l+0x0] = jp+  0; dft_offsets[l+0x1] = jp+p01; dft_offsets[l+0x2] = jp+p02; dft_offsets[l+0x3] = jp+p03; dft_offsets[l+0x4] = jp+p04; dft_offsets[l+0x5] = jp+p05; dft_offsets[l+0x6] = jp+p06; dft_offsets[l+0x7] = jp+p07;
		l = 0x18;	jp = p18;
		dft_offsets[l+0x0] = jp+  0; dft_offsets[l+0x1] = jp+p01; dft_offsets[l+0x2] = jp+p02; dft_offsets[l+0x3] = jp+p03; dft_offsets[l+0x4] = jp+p04; dft_offsets[l+0x5] = jp+p05; dft_offsets[l+0x6] = jp+p06; dft_offsets[l+0x7] = jp+p07;
		l = 0x20;	jp = p20;
		dft_offsets[l+0x0] = jp+  0; dft_offsets[l+0x1] = jp+p01; dft_offsets[l+0x2] = jp+p02; dft_offsets[l+0x3] = jp+p03; dft_offsets[l+0x4] = jp+p04; dft_offsets[l+0x5] = jp+p05; dft_offsets[l+0x6] = jp+p06; dft_offsets[l+0x7] = jp+p07;
		l = 0x28;	jp = p28;
		dft_offsets[l+0x0] = jp+  0; dft_offsets[l+0x1] = jp+p01; dft_offsets[l+0x2] = jp+p02; dft_offsets[l+0x3] = jp+p03; dft_offsets[l+0x4] = jp+p04; dft_offsets[l+0x5] = jp+p05; dft_offsets[l+0x6] = jp+p06; dft_offsets[l+0x7] = jp+p07;
		l = 0x30;	jp = p30;
		dft_offsets[l+0x0] = jp+  0; dft_offsets[l+0x1] = jp+p01; dft_offsets[l+0x2] = jp+p02; dft_offsets[l+0x3] = jp+p03; dft_offsets[l+0x4] = jp+p04; dft_offsets[l+0x5] = jp+p05; dft_offsets[l+0x6] = jp+p06; dft_offsets[l+0x7] = jp+p07;
		l = 0x38;	jp = p38;
		dft_offsets[l+0x0] = jp+  0; dft_offsets[l+0x1] = jp+p01; dft_offsets[l+0x2] = jp+p02; dft_offsets[l+0x3] = jp+p03; dft_offsets[l+0x4] = jp+p04; dft_offsets[l+0x5] = jp+p05; dft_offsets[l+0x6] = jp+p06; dft_offsets[l+0x7] = jp+p07;
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

/*...The radix-64 final DIT pass is here.	*/

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
		#include "radix64_main_carry_loop.h"

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
		ASSERT(0x0 == cy64_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
	ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
	ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(0 == mlucas_nanosleep(&ns_time), "nanosleep fail!");
	}
//	printf("radix64_ditN_cy_dif1 end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

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
			for(l = 0; l < RADIX>>2; l++) {
				jt = j1 + poff[l];
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

/***************/

void radix64_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-64 complex DIF FFT pass on the data in the length-N real vector A.
!
!   The data are stored in a 1-D zero-offset array, with 2^PAD_BITS 8-byte padding elements inserted
!   between every block of 2^DAT_BITS contiguous data. The array padding is to prevent data being accessed
!   in strides that are large powers of two and thus to minimize cache thrashing
!   (in cache-based microprocessor architectures) or bank conflicts (in supercomputers.)
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1;
#if !USE_SCALAR_DFT_MACRO
	int j2,jt,jp;
#endif
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30,p38, first_entry=TRUE;
#if USE_SCALAR_DFT_MACRO
	static int i_offsets[64];
#else
	double
		 t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F
		,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F
		,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F
		,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t4A,t4B,t4C,t4D,t4E,t4F
		,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t5A,t5B,t5C,t5D,t5E,t5F
		,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t6A,t6B,t6C,t6D,t6E,t6F
		,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t7A,t7B,t7C,t7D,t7E,t7F;
#endif

	if(!first_entry && (n >> 6) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR = n >> 6;

		p01 = NDIVR;
		p02 = p01 +p01;
		p03 = p02 +p01;
		p04 = p03 +p01;
		p05 = p04 +p01;
		p06 = p05 +p01;
		p07 = p06 +p01;
		p08 = p07 +p01;
		p10 = p08 +p08;
		p18 = p10 +p08;
		p20 = p18 +p08;
		p28 = p20 +p08;
		p30 = p28 +p08;
		p38 = p30 +p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p30 = p30 + ( (p30 >> DAT_BITS) << PAD_BITS );
		p38 = p38 + ( (p38 >> DAT_BITS) << PAD_BITS );

	#if USE_SCALAR_DFT_MACRO
		// Set array offsets for in/outputs:
		i_offsets[0x00] = 0  ;
		i_offsets[0x01] = p01;
		i_offsets[0x02] = p02;
		i_offsets[0x03] = p03;
		i_offsets[0x04] = p04;
		i_offsets[0x05] = p05;
		i_offsets[0x06] = p06;
		i_offsets[0x07] = p07;
		for(j = 0; j < 8; ++j) { i_offsets[j+0x08] = i_offsets[j] + p08; }
		for(j = 0; j < 8; ++j) { i_offsets[j+0x10] = i_offsets[j] + p10; }
		for(j = 0; j < 8; ++j) { i_offsets[j+0x18] = i_offsets[j] + p18; }
		for(j = 0; j < 8; ++j) { i_offsets[j+0x20] = i_offsets[j] + p20; }
		for(j = 0; j < 8; ++j) { i_offsets[j+0x28] = i_offsets[j] + p28; }
		for(j = 0; j < 8; ++j) { i_offsets[j+0x30] = i_offsets[j] + p30; }
		for(j = 0; j < 8; ++j) { i_offsets[j+0x38] = i_offsets[j] + p38; }
	#endif
	}

/*...The radix-64 pass is here.	*/

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
	#if !USE_SCALAR_DFT_MACRO
		j2 = j1+RE_IM_STRIDE;
	#endif

	#if USE_SCALAR_DFT_MACRO

		RADIX_64_DIF(
			(a+j1),i_offsets,RE_IM_STRIDE,
			(a+j1),i_offsets,RE_IM_STRIDE
		);

	#else

	/* Gather the needed data (64 64-bit complex, i.e. 128 64-bit reals) and do 8 twiddleless length-8 subtransforms: */
		/*...Block 0: */
		jt = j1;	jp = j2;
		RADIX_08_DIF_OOP(
			a[jt],a[jp],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],
			t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F
	    );
		/*...Block 1: */
		jt = j1 + p04;	jp = j2 + p04;
		RADIX_08_DIF_OOP(
			a[jt],a[jp],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],
			t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F
	    );
		/*...Block 2: */
		jt = j1 + p02;	jp = j2 + p02;
		RADIX_08_DIF_OOP(
			a[jt],a[jp],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],
			t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F
	    );
		/*...Block 3: */
		jt = j1 + p06;	jp = j2 + p06;
		RADIX_08_DIF_OOP(
			a[jt],a[jp],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],
			t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F
	    );
		/*...Block 4: */
		jt = j1 + p01;	jp = j2 + p01;
		RADIX_08_DIF_OOP(
			a[jt],a[jp],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],
			t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t4A,t4B,t4C,t4D,t4E,t4F
	    );
		/*...Block 5: */
		jt = j1 + p05;	jp = j2 + p05;
		RADIX_08_DIF_OOP(
			a[jt],a[jp],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],
			t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t5A,t5B,t5C,t5D,t5E,t5F
	    );
		/*...Block 6: */
		jt = j1 + p03;	jp = j2 + p03;
		RADIX_08_DIF_OOP(
			a[jt],a[jp],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],
			t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t6A,t6B,t6C,t6D,t6E,t6F
	    );
		/*...Block 7: */
		jt = j1 + p07;	jp = j2 + p07;
		RADIX_08_DIF_OOP(
			a[jt],a[jp],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],
			t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t7A,t7B,t7C,t7D,t7E,t7F
	    );

	/*...and now do eight radix-8 subtransforms, including the internal twiddle factors:

		Block 0: twiddles = {   1,   1,   1,   1,   1,   1,   1,   1}
		Block 1: twiddles = {   1,E^ 1,E^ 2,E^ 3,E^ 4,E^ 5,E^ 6,E^ 7}
		Block 2: twiddles = {   1,E^ 2,E^ 4,E^ 6,E^ 8,E^10,E^12,E^14}
		Block 3: twiddles = {   1,E^ 3,E^ 6,E^ 9,E^12,E^15,E^18,E^21}
		Block 4: twiddles = {   1,E^ 4,E^ 8,E^12,E^16,E^20,E^24,E^28}
		Block 5: twiddles = {   1,E^ 5,E^10,E^15,E^20,E^25,E^30,E^35}
		Block 6: twiddles = {   1,E^ 6,E^12,E^18,E^24,E^30,E^36,E^42}
		Block 7: twiddles = {   1,E^ 7,E^14,E^21,E^28,E^35,E^42,E^49}

	where ~ denotes complex conjugation, * denotes interchange of real and imaginary part, - denotes negation
	(i.e. ~E^n := ( Re(E),-Im(E)), *E^n := ( Im(E), Re(E)), -E^n := (-Re(E),-Im(E)), and any sequence of these
	operators is evaluated right-to-left, e.g. -~*E^n = -~(*E^n) = -~( Im(E), Re(E)) = -( Im(E),-Re(E)) = (-Im(E),+Re(E));
	note that the - and ~ operators commute, as do the - and *, but ~ and * anticommute, i.e. ~*E^n = -*~E^n) ,
	and the values of the individual exponentials are, in terms of the sincos parameters defined in this module:

		E^ 1 = exp(i* 1*twopi/64) =       ( c64_1, s64_1)
		E^ 2 = exp(i* 2*twopi/64) =       ( c32_1, s32_1)
		E^ 3 = exp(i* 3*twopi/64) =       ( c64_3, s64_3)
		E^ 4 = exp(i* 4*twopi/64) =       ( c16  , s16  )
		E^ 5 = exp(i* 5*twopi/64) =       ( c64_5, s64_5)
		E^ 6 = exp(i* 6*twopi/64) =       ( c32_3, s32_3)
		E^ 7 = exp(i* 7*twopi/64) =       ( c64_7, s64_7)
		E^ 8 = exp(i* 8*twopi/64) = isrt2*( 1    , 1    )
		E^16 = exp(i*16*twopi/64) =       ( 0    , 1    ) = I
		E^24 = exp(i*24*twopi/64) = isrt2*(-1    , 1    ) = *~E^ 8

		E^{ 9,10,11,12,13,14,15}  = *E^{ 7, 6, 5, 4, 3, 2, 1}

		E^{17,18,19,20,21,22,23}  =I.E^{ 1, 2, 3, 4, 5, 6, 7} =*~E^{ 1, 2, 3, 4, 5, 6, 7} (I.{} denotes complex multiply by I)

		E^{25,26,27,28,29,30,31}  =-~E^{ 7, 6, 5, 4, 3, 2, 1}

	and using that E^n = -E^(n-32), we get

		E^{33,34,35,36,37,38,39}  = -E^{ 1, 2, 3, 4, 5, 6, 7}

		E^42 =-E^10 = -*E^ 6

		E^49 =-E^17 = -*~E^ 1 = ~*E^ 1 .

	Thus, exploiting these symmetries allows our 8x8 twiddles matrix to be expressed in terms of the powers E^1-7, 1
	and the imaginary constant I as:

		Block 1: twiddles = {   1,     1,     1,     1,     1,     1,     1,     1}
		Block 2: twiddles = {   1,  E^ 1,  E^ 2,  E^ 3,  E^ 4,  E^ 5,  E^ 6,  E^ 7}
		Block 3: twiddles = {   1,  E^ 2,  E^ 4,  E^ 6,  E^ 8, *E^ 6, *E^ 4, *E^ 2}
		Block 4: twiddles = {   1,  E^ 3,  E^ 6, *E^ 7, *E^ 4, *E^ 1,*~E^ 2,*~E^ 5}
		Block 5: twiddles = {   1,  E^ 4,  E^ 8, *E^ 4,  I.{},*~E^ 4,*~E^ 8,-~E^ 4}
		Block 6: twiddles = {   1,  E^ 5, *E^ 6, *E^ 1,*~E^ 4,-~E^ 7,-~E^ 2, -E^ 3}
		Block 7: twiddles = {   1,  E^ 6, *E^ 4,*~E^ 2,*~E^ 8,-~E^ 2, -E^ 4,-*E^ 6}
		Block 8: twiddles = {   1,  E^ 7, *E^ 2,*~E^ 5,-~E^ 4, -E^ 3,-*E^ 6,~*E^ 1} ,

	and only the last 7 inputs to each of the radix-8 transforms 2 through 8 are multiplied by non-unity twiddles.
	For DIF we process both the blocks, and the twiddles within each block, in bit-reversed order.
	One can see from the data below that aside from using a twiddleless DIF for Block 0 there is
	little to be gained from trying to exploit other "special form" twiddles such as I and isrt2*[+-1,+-1],
	as only 5 of the remaining 7*8 = 56 non-unity twiddles have such a special form (one = I, four of the isrt2-form).
	Thus our radix-8-DIF-with-twiddles macro uses generic complex MUL for the 7 non-unity twiddles.
	*/
		/* Block 0: t*0,1 twiddles = {   1,     1,     1,     1,     1,     1,     1,     1}: */
		jt = j1;	jp = j2;
		/* 0-index block has all-unity twiddles: Remember, the twiddleless DIT also bit-reverses its outputs, so a_p* terms appear in BR-order here [swap index pairs 1/4 and 3/6]: */
		RADIX_08_DIF_OOP(
			t00,t01,t40,t41,t20,t21,t60,t61,t10,t11,t50,t51,t30,t31,t70,t71,
			a[jt],a[jp],a[jt+p04],a[jp+p04],a[jt+p02],a[jp+p02],a[jt+p06],a[jp+p06],a[jt+p01],a[jp+p01],a[jt+p05],a[jp+p05],a[jt+p03],a[jp+p03],a[jt+p07],a[jp+p07]
		);
	/* Here is the with-twiddles version of the above macro call:
		RADIX_08_DIF_TWIDDLE_OOP(
			t00,t01,t10,t11,t20,t21,t30,t31,t40,t41,t50,t51,t60,t61,t70,t71,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			1.0,0.0, 1.0,0.0, 1.0,0.0, 1.0,0.0, 1.0,0.0, 1.0,0.0, 1.0,0.0
		);
	*/
		/* Block 4: t*8,9, twiddles = {   1,  E^ 4,  E^ 8, *E^ 4,  I.{},*~E^ 4,*~E^ 8,-~E^ 4}
							BR order: {   1,  I.{},  E^ 8,*~E^ 8,  E^ 4,*~E^ 4, *E^ 4,-~E^ 4}
			In terms of [Re,Im] pairs:{[1,0],[0,1],isrt2*[1,1],isrt2*[-1,1],[c4,s4],[-s4,c4],[s4,c4],[-c4,s4]}
		*/
		jt = j1 + p08;	jp = j2 + p08;
		RADIX_08_DIF_TWIDDLE_OOP(
			t08,t09,t18,t19,t28,t29,t38,t39,t48,t49,t58,t59,t68,t69,t78,t79,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			0.0,1.0,ISRT2,ISRT2,-ISRT2,ISRT2,c16  ,s16  ,-s16  ,c16  ,s16  ,c16  ,-c16  ,s16
		);
		/* Block 2: t*4,5, twiddles = {   1,  E^ 2,  E^ 4,  E^ 6,  E^ 8, *E^ 6, *E^ 4, *E^ 2}
							BR order: {   1,  E^ 8,  E^ 4, *E^ 4,  E^ 2, *E^ 6,  E^ 6, *E^ 2}
			In terms of [Re,Im] pairs:{[1,0],isrt2*[1,1],[c4,s4],[s4,c4],[c2,s2],[s6,c6],[c6,s6],[s2,c2]}
		*/
		jt = j1 + p10;	jp = j2 + p10;
		RADIX_08_DIF_TWIDDLE_OOP(
			t04,t05,t14,t15,t24,t25,t34,t35,t44,t45,t54,t55,t64,t65,t74,t75,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			ISRT2,ISRT2,c16  ,s16  ,s16  ,c16  ,c32_1,s32_1,s32_3,c32_3,c32_3,s32_3,s32_1,c32_1
		);
		/* Block 6: t*C,D, twiddles = {   1,  E^ 6, *E^ 4,*~E^ 2,*~E^ 8,-~E^ 2, -E^ 4,-*E^ 6}
							BR order: {   1,*~E^ 8, *E^ 4, -E^ 4,  E^ 6,-~E^ 2,*~E^ 2,-*E^ 6}
			In terms of [Re,Im] pairs:{[1,0],isrt2*[-1,1],[s4,c4],[-c4,-s4],[c6,s6],[-c2,s2],[-s2,c2],[-s6,-c6]}
		*/
		jt = j1 + p18;	jp = j2 + p18;
		RADIX_08_DIF_TWIDDLE_OOP(
			t0C,t0D,t1C,t1D,t2C,t2D,t3C,t3D,t4C,t4D,t5C,t5D,t6C,t6D,t7C,t7D,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			-ISRT2,ISRT2,s16  ,c16  ,-c16  ,-s16  ,c32_3,s32_3,-c32_1,s32_1,-s32_1,c32_1,-s32_3,-c32_3
		);
		/* Block 1: t*2,3, twiddles = {   1,  E^ 1,  E^ 2,  E^ 3,  E^ 4,  E^ 5,  E^ 6,  E^ 7}
							BR order: {   1,  E^ 4,  E^ 2,  E^ 6,  E^ 1,  E^ 5,  E^ 3,  E^ 7}
			In terms of [Re,Im] pairs:{[1,0],[c4,s4],[c2,s2],[c6,s6],[c1,s1],[c5,s5],[c3,s3],[c7,s7]}
		*/
		jt = j1 + p20;	jp = j2 + p20;
		RADIX_08_DIF_TWIDDLE_OOP(
			t02,t03,t12,t13,t22,t23,t32,t33,t42,t43,t52,t53,t62,t63,t72,t73,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			c16  ,s16  ,c32_1,s32_1,c32_3,s32_3,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7
		);
		/* Block 5: t*A,B, twiddles = {   1,  E^ 5, *E^ 6, *E^ 1,*~E^ 4,-~E^ 7,-~E^ 2, -E^ 3}
							BR order: {   1,*~E^ 4, *E^ 6,-~E^ 2,  E^ 5,-~E^ 7, *E^ 1, -E^ 3}
			In terms of [Re,Im] pairs:{[1,0],[-s4,c4],[s6,c6],[-c2,s2],[c5,s5],[-c7,s7],[s1,c1],[-c3,-s3]}
		*/
		jt = j1 + p28;	jp = j2 + p28;
		RADIX_08_DIF_TWIDDLE_OOP(
			t0A,t0B,t1A,t1B,t2A,t2B,t3A,t3B,t4A,t4B,t5A,t5B,t6A,t6B,t7A,t7B,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			-s16  ,c16  ,s32_3,c32_3,-c32_1,s32_1,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3
		);
		/* Block 3: t*6,7, twiddles = {   1,  E^ 3,  E^ 6, *E^ 7, *E^ 4, *E^ 1,*~E^ 2,*~E^ 5}
							BR order: {   1, *E^ 4,  E^ 6,*~E^ 2,  E^ 3, *E^ 1, *E^ 7,*~E^ 5}
			In terms of [Re,Im] pairs:{[1,0],[s4,c4],[c6,s6],[-s2,c2],[c3,s3],[s1,c1],[s7,c7],[-s5,c5]}
		*/
		jt = j1 + p30;	jp = j2 + p30;
		RADIX_08_DIF_TWIDDLE_OOP(
			t06,t07,t16,t17,t26,t27,t36,t37,t46,t47,t56,t57,t66,t67,t76,t77,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			s16  ,c16  ,c32_3,s32_3,-s32_1,c32_1,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5
		);
		/* Block 7: t*E,F, twiddles = {   1,  E^ 7, *E^ 2,*~E^ 5,-~E^ 4, -E^ 3,-*E^ 6,~*E^ 1} ,
							BR order: {   1,-~E^ 4, *E^ 2,-*E^ 6,  E^ 7, -E^ 3,*~E^ 5,~*E^ 1} ,
			In terms of [Re,Im] pairs:{[1,0],[-c4,s4],[s2,c2],[-s6,-c6],[c7,s7],[-c3,-s3],[-s5,c5],[s1,-c1]}
		*/
		jt = j1 + p38;	jp = j2 + p38;
		RADIX_08_DIF_TWIDDLE_OOP(
			t0E,t0F,t1E,t1F,t2E,t2F,t3E,t3F,t4E,t4F,t5E,t5F,t6E,t6F,t7E,t7F,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			-c16  ,s16  ,s32_1,c32_1,-s32_3,-c32_3,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1
		);
	#endif	// USE_SCALAR_DFT_MACRO ?
	}
}

/**************/

void radix64_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-64 complex inverse DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details.
*/
	int j,j1;
#if !USE_SCALAR_DFT_MACRO
	int j2,jt,jp;
#endif
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30,p38, first_entry=TRUE;
#if USE_SCALAR_DFT_MACRO
	static int i_offsets[64];
#else
	double
		 t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F
		,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F
		,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F
		,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t4A,t4B,t4C,t4D,t4E,t4F
		,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t5A,t5B,t5C,t5D,t5E,t5F
		,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t6A,t6B,t6C,t6D,t6E,t6F
		,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t7A,t7B,t7C,t7D,t7E,t7F;
#endif

	if(!first_entry && (n >> 6) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR = n >> 6;

		p01 = NDIVR;
		p02 = p01 +p01;
		p03 = p02 +p01;
		p04 = p03 +p01;
		p05 = p04 +p01;
		p06 = p05 +p01;
		p07 = p06 +p01;
		p08 = p07 +p01;
		p10 = p08 +p08;
		p18 = p10 +p08;
		p20 = p18 +p08;
		p28 = p20 +p08;
		p30 = p28 +p08;
		p38 = p30 +p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p30 = p30 + ( (p30 >> DAT_BITS) << PAD_BITS );
		p38 = p38 + ( (p38 >> DAT_BITS) << PAD_BITS );

	#if USE_SCALAR_DFT_MACRO
		// Set array offsets for radix-64 DFT in/outputs:
		// set 1 is w.r.to: a-array:								// set 2 w.r.to t-array:
		i_offsets[0x00] = 0  ;
		i_offsets[0x01] = p01;
		i_offsets[0x02] = p02;
		i_offsets[0x03] = p03;
		i_offsets[0x04] = p04;
		i_offsets[0x05] = p05;
		i_offsets[0x06] = p06;
		i_offsets[0x07] = p07;
		for(j = 0; j < 8; ++j) { i_offsets[j+0x08] = i_offsets[j] + p08; }
		for(j = 0; j < 8; ++j) { i_offsets[j+0x10] = i_offsets[j] + p10; }
		for(j = 0; j < 8; ++j) { i_offsets[j+0x18] = i_offsets[j] + p18; }
		for(j = 0; j < 8; ++j) { i_offsets[j+0x20] = i_offsets[j] + p20; }
		for(j = 0; j < 8; ++j) { i_offsets[j+0x28] = i_offsets[j] + p28; }
		for(j = 0; j < 8; ++j) { i_offsets[j+0x30] = i_offsets[j] + p30; }
		for(j = 0; j < 8; ++j) { i_offsets[j+0x38] = i_offsets[j] + p38; }
	#endif
	}

/*...The radix-64 pass is here.	*/

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
	#if !USE_SCALAR_DFT_MACRO
		j2 = j1+RE_IM_STRIDE;
	#endif

	#if USE_SCALAR_DFT_MACRO

		RADIX_64_DIT(
			(a+j1),i_offsets,RE_IM_STRIDE,
			(a+j1),i_offsets,RE_IM_STRIDE
		);

	#else

	/* Gather the needed data (64 64-bit complex, i.e. 128 64-bit reals) and do 8 twiddleless length-8 subtransforms: */
		/*...Block 0: */
		jt = j1;	jp = j2;
		RADIX_08_DIT_OOP(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F
	    );
		/*...Block 1: */
		jt = j1 + p08;	jp = j2 + p08;
		RADIX_08_DIT_OOP(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F
	    );
		/*...Block 2: */
		jt = j1 + p10;	jp = j2 + p10;
		RADIX_08_DIT_OOP(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F
	    );
		/*...Block 3: */
		jt = j1 + p18;	jp = j2 + p18;
		RADIX_08_DIT_OOP(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F
	    );
		/*...Block 4: */
		jt = j1 + p20;	jp = j2 + p20;
		RADIX_08_DIT_OOP(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t4A,t4B,t4C,t4D,t4E,t4F
	    );
		/*...Block 5: */
		jt = j1 + p28;	jp = j2 + p28;
		RADIX_08_DIT_OOP(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t5A,t5B,t5C,t5D,t5E,t5F
	    );
		/*...Block 6: */
		jt = j1 + p30;	jp = j2 + p30;
		RADIX_08_DIT_OOP(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t6A,t6B,t6C,t6D,t6E,t6F
	    );
		/*...Block 7: */
		jt = j1 + p38;	jp = j2 + p38;
		RADIX_08_DIT_OOP(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t7A,t7B,t7C,t7D,t7E,t7F
	    );

	/*...and now do eight radix-8 subtransforms, including the internal twiddle factors - we use the same positive-power
	roots as in the DIF here, just fiddle with signs along the way to effect the conjugate-multiplies.
	The twiddles occur in the same order here as for DIF, but now the output-index offsets are BRed: j1 + p[0,4,2,6,1,5,3,7],
	as are the index offsets of each sets of complex outputs in the A-array: [jt,jp] + p08*[0,4,2,6,1,5,3,7]:
	*/
		/* Block 0: t*0,1 twiddles = {   1,     1,     1,     1,     1,     1,     1,     1}: */
		jt = j1;	jp = j2;
		/* 0-index block has all-unity twiddles: Remember, the twiddleless DIT also bit-reverses its outputs, so a_p* terms appear in-order here: */
		RADIX_08_DIT_OOP(
			t00,t01,t10,t11,t20,t21,t30,t31,t40,t41,t50,t51,t60,t61,t70,t71,
			a[jt],a[jp],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38]
		);
	/* Here is the with-twiddles version of the above macro call:
		RADIX_08_DIT_TWIDDLE_OOP(
			t00,t01,t10,t11,t20,t21,t30,t31,t40,t41,t50,t51,t60,t61,t70,t71,
			a[jt],a[jp],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt+p30],a[jp+p30],a[jt+p08],a[jp+p08],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p38],a[jp+p38],
			1.0,0.0, 1.0,0.0, 1.0,0.0, 1.0,0.0, 1.0,0.0, 1.0,0.0, 1.0,0.0
		);
	*/
		/* Block 4: t*8,9, twiddles = ~{   1,  E^ 4,  E^ 8, *E^ 4,  I.{},*~E^ 4,*~E^ 8,-~E^ 4}
							BR order: ~{   1,  I.{},  E^ 8,*~E^ 8,  E^ 4,*~E^ 4, *E^ 4,-~E^ 4}
			In terms of [Re,Im] pairs:~{[1,0],[0,1],isrt2*[1,1],isrt2*[-1,1],[c4,s4],[-s4,c4],[s4,c4],[-c4,s4]}
		*/
		jt = j1 + p04;	jp = j2 + p04;
		RADIX_08_DIT_TWIDDLE_OOP(
			t08,t09,t18,t19,t28,t29,t38,t39,t48,t49,t58,t59,t68,t69,t78,t79,
			a[jt],a[jp],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt+p30],a[jp+p30],a[jt+p08],a[jp+p08],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p38],a[jp+p38],
			0.0,1.0,ISRT2,ISRT2,-ISRT2,ISRT2,c16  ,s16  ,-s16  ,c16  ,s16  ,c16  ,-c16  ,s16
		);
		/* Block 2: t*4,5, twiddles = ~{   1,  E^ 2,  E^ 4,  E^ 6,  E^ 8, *E^ 6, *E^ 4, *E^ 2}
							BR order: ~{   1,  E^ 8,  E^ 4, *E^ 4,  E^ 2, *E^ 6,  E^ 6, *E^ 2}
			In terms of [Re,Im] pairs:~{[1,0],isrt2*[1,1],[c4,s4],[s4,c4],[c2,s2],[s6,c6],[c6,s6],[s2,c2]}
		*/
		jt = j1 + p02;	jp = j2 + p02;
		RADIX_08_DIT_TWIDDLE_OOP(
			t04,t05,t14,t15,t24,t25,t34,t35,t44,t45,t54,t55,t64,t65,t74,t75,
			a[jt],a[jp],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt+p30],a[jp+p30],a[jt+p08],a[jp+p08],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p38],a[jp+p38],
			ISRT2,ISRT2,c16  ,s16  ,s16  ,c16  ,c32_1,s32_1,s32_3,c32_3,c32_3,s32_3,s32_1,c32_1
		);
		/* Block 6: t*C,D, twiddles = ~{   1,  E^ 6, *E^ 4,*~E^ 2,*~E^ 8,-~E^ 2, -E^ 4,-*E^ 6}
							BR order: ~{   1,*~E^ 8, *E^ 4, -E^ 4,  E^ 6,-~E^ 2,*~E^ 2,-*E^ 6}
			In terms of [Re,Im] pairs:~{[1,0],isrt2*[-1,1],[s4,c4],[-c4,-s4],[c6,s6],[-c2,s2],[-s2,c2],[-s6,-c6]}
		*/
		jt = j1 + p06;	jp = j2 + p06;
		RADIX_08_DIT_TWIDDLE_OOP(
			t0C,t0D,t1C,t1D,t2C,t2D,t3C,t3D,t4C,t4D,t5C,t5D,t6C,t6D,t7C,t7D,
			a[jt],a[jp],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt+p30],a[jp+p30],a[jt+p08],a[jp+p08],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p38],a[jp+p38],
			-ISRT2,ISRT2,s16  ,c16  ,-c16  ,-s16  ,c32_3,s32_3,-c32_1,s32_1,-s32_1,c32_1,-s32_3,-c32_3
		);
		/* Block 1: t*2,3, twiddles = ~{   1,  E^ 1,  E^ 2,  E^ 3,  E^ 4,  E^ 5,  E^ 6,  E^ 7}
							BR order: ~{   1,  E^ 4,  E^ 2,  E^ 6,  E^ 1,  E^ 5,  E^ 3,  E^ 7}
			In terms of [Re,Im] pairs:~{[1,0],[c4,s4],[c2,s2],[c6,s6],[c1,s1],[c5,s5],[c3,s3],[c7,s7]}
		*/
		jt = j1 + p01;	jp = j2 + p01;
		RADIX_08_DIT_TWIDDLE_OOP(
			t02,t03,t12,t13,t22,t23,t32,t33,t42,t43,t52,t53,t62,t63,t72,t73,
			a[jt],a[jp],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt+p30],a[jp+p30],a[jt+p08],a[jp+p08],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p38],a[jp+p38],
			c16  ,s16  ,c32_1,s32_1,c32_3,s32_3,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7
		);
		/* Block 5: t*A,B, twiddles = ~{   1,  E^ 5, *E^ 6, *E^ 1,*~E^ 4,-~E^ 7,-~E^ 2, -E^ 3}
							BR order: ~{   1,*~E^ 4, *E^ 6,-~E^ 2,  E^ 5,-~E^ 7, *E^ 1, -E^ 3}
			In terms of [Re,Im] pairs:~{[1,0],[-s4,c4],[s6,c6],[-c2,s2],[c5,s5],[-c7,s7],[s1,c1],[-c3,-s3]}
		*/
		jt = j1 + p05;	jp = j2 + p05;
		RADIX_08_DIT_TWIDDLE_OOP(
			t0A,t0B,t1A,t1B,t2A,t2B,t3A,t3B,t4A,t4B,t5A,t5B,t6A,t6B,t7A,t7B,
			a[jt],a[jp],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt+p30],a[jp+p30],a[jt+p08],a[jp+p08],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p38],a[jp+p38],
			-s16  ,c16  ,s32_3,c32_3,-c32_1,s32_1,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3
		);
		/* Block 3: t*6,7, twiddles = ~{   1,  E^ 3,  E^ 6, *E^ 7, *E^ 4, *E^ 1,*~E^ 2,*~E^ 5}
							BR order: ~{   1, *E^ 4,  E^ 6,*~E^ 2,  E^ 3, *E^ 1, *E^ 7,*~E^ 5}
			In terms of [Re,Im] pairs:~{[1,0],[s4,c4],[c6,s6],[-s2,c2],[c3,s3],[s1,c1],[s7,c7],[-s5,c5]}
		*/
		jt = j1 + p03;	jp = j2 + p03;
		RADIX_08_DIT_TWIDDLE_OOP(
			t06,t07,t16,t17,t26,t27,t36,t37,t46,t47,t56,t57,t66,t67,t76,t77,
			a[jt],a[jp],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt+p30],a[jp+p30],a[jt+p08],a[jp+p08],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p38],a[jp+p38],
			s16  ,c16  ,c32_3,s32_3,-s32_1,c32_1,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5
		);
		/* Block 7: t*E,F, twiddles = ~{   1,  E^ 7, *E^ 2,*~E^ 5,-~E^ 4, -E^ 3,-*E^ 6,~*E^ 1} ,
							BR order: ~{   1,-~E^ 4, *E^ 2,-*E^ 6,  E^ 7, -E^ 3,*~E^ 5,~*E^ 1} ,
			In terms of [Re,Im] pairs:~{[1,0],[-c4,s4],[s2,c2],[-s6,-c6],[c7,s7],[-c3,-s3],[-s5,c5],[s1,-c1]}
		*/
		jt = j1 + p07;	jp = j2 + p07;
		RADIX_08_DIT_TWIDDLE_OOP(
			t0E,t0F,t1E,t1F,t2E,t2F,t3E,t3F,t4E,t4F,t5E,t5F,t6E,t6F,t7E,t7F,
			a[jt],a[jp],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt+p30],a[jp+p30],a[jt+p08],a[jp+p08],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p38],a[jp+p38],
			-c16  ,s16  ,s32_1,c32_1,-s32_3,-c32_3,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1
		);
	#endif	// USE_SCALAR_DFT_MACRO ?
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy64_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr,*addi;
	#ifndef USE_SSE2
		struct complex *tptr;
	#else
		const int pfetch_dist = PFETCH_DIST;
	#endif
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30,p38;
		int poff[RADIX>>2];	// Store mults of p-offsets for loop-controlled DFT macro calls
	#if defined(USE_SSE2) && !USE_SCALAR_DFT_MACRO
		int po_lin[8];		// Store mults of p-offsets for loop-controlled DFT macro calls
	#endif
	#ifndef USE_SSE2
		int po_br[8];		// Store mults of p-offsets for loop-controlled DFT macro calls
	#endif
	#if defined(USE_SSE2) && USE_SCALAR_DFT_MACRO
		static int dft_offsets[RADIX], c_offsets[RADIX];
	#endif
		int j,j1,l;
	#ifndef USE_SSE2
		int j2,jt,jp;
	#endif
	#if defined(USE_SSE2) && !defined(USE_AVX)
		int incr;
		// incr = Carry-chain wts-multipliers recurrence length, which must divide
		// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 4|8|16 for avx512,avx,sse, respectively:
		const int incr_long = 16,incr_med = 8,incr_short = 4;
		// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
		if(USE_SHORT_CY_CHAIN == 0)
			incr = incr_long;
		else if(USE_SHORT_CY_CHAIN == 1)
			incr = incr_med;
		else
			incr = incr_short;
	#endif
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
	#if !defined(USE_SSE2) || defined(USE_AVX)
		double rt,it, wt_re,wt_im;	// Fermat-mod weights stuff, used in both scalar and AVX mode
	#endif
	#ifndef USE_SSE2
		double wi_re,wi_im;	// Fermat-mod weights stuff, used in both scalar and AVX mode
	#endif
		int k1;
	#if !defined(USE_SSE2) || defined(USE_AVX)
		int k2;
	#endif

	#ifdef USE_SSE2

	  #ifndef USE_AVX512
		const double crnd = 3.0*0x4000000*0x2000000;
	  #endif
		double *add0,*add1,*add2,*add3,*add4,*add5,*add6,*add7;
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *tmp,*tm1,*tm2;	// utility ptrs
	  #ifdef USE_AVX
		vec_dbl *tm0;
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
		#if !USE_SCALAR_DFT_MACRO
			*cc0, *ss0, *two,/* *one,*sqrt2, */
			 *isrt2, *cc1, *ss1, *cc2, *ss2, *cc3, *ss3, *cc4, *ss4, *cc5, *ss5, *cc6, *ss6, *cc7, *ss7,
			*nisrt2,*ncc1,*ncc2,*nss2,*ncc3,*nss3,*ncc4,*nss4,/* *ncc5, */*nss5,*ncc6,*nss6,*ncc7,/* *nss7, */	// each non-unity root now needs a negated counterpart
		#endif
			// This routine dates from the benighted early days when I treated r's as pointers-to-real-SIMD...
			*r00,*r02,*r04,*r06,*r08,*r0A,*r0C,*r0E,*r10,/* *r12,*r14,*r16,*r18,*r1A,*r1C,*r1E, */
			*r20,/* *r22,*r24,*r26,*r28,*r2A,*r2C,*r2E, */*r30,/* *r32,*r34,*r36,*r38,*r3A,*r3C,*r3E, */
			*r40,/* *r42,*r44,*r46,*r48,*r4A,*r4C,*r4E, */*r50,/* *r52,*r54,*r56,*r58,*r5A,*r5C,*r5E, */
			*r60,/* *r62,*r64,*r66,*r68,*r6A,*r6C,*r6E, */*r70,/* *r72,*r74,*r76,*r78,*r7A,*r7C,*r7E, */
			// ...and s's as pointers-to-complex-SIMD; thus the r-indices run 2x faster than the s-ones:
			*s1p00,*s1p01,*s1p02,*s1p03,*s1p04,*s1p05,*s1p06,*s1p07,*s1p08,/* *s1p09,*s1p0a,*s1p0b,*s1p0c,*s1p0d,*s1p0e,*s1p0f, */
			*s1p10,/* *s1p11,*s1p12,*s1p13,*s1p14,*s1p15,*s1p16,*s1p17,*/ *s1p18,/* *s1p19,*s1p1a,*s1p1b,*s1p1c,*s1p1d,*s1p1e,*s1p1f, */
			*s1p20,/* *s1p21,*s1p22,*s1p23,*s1p24,*s1p25,*s1p26,*s1p27,*/ *s1p28,/* *s1p29,*s1p2a,*s1p2b,*s1p2c,*s1p2d,*s1p2e,*s1p2f, */
			*s1p30,/* *s1p31,*s1p32,*s1p33,*s1p34,*s1p35,*s1p36,*s1p37,*/ *s1p38,/* *s1p39,*s1p3a,*s1p3b,*s1p3c,*s1p3d,*s1p3e,*s1p3f, */
			*cy_r	// Need RADIX slots for sse2 carries, RADIX/2 for avx
		#ifdef USE_AVX
			,*cy_i	// Need RADIX slots for sse2 carries, RADIX/2 for avx
		#endif
			;
	  #ifdef USE_AVX
		vec_dbl *base_negacyclic_root;
	  #endif

	  #ifndef USE_AVX
		/* These are used in conjunction with the langth-odd_radix arrays in the USE_SCALAR_CARRY #define below;
		In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
		int idx_offset, idx_incr;
	  #endif
		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_nm1;

	#else

		double *base, *baseinv;
		int p0123[4];
		const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int m,m2,ntmp;
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
	#ifndef USE_SSE2
		int nm1 = n-1;
	#endif
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
	#ifndef USE_SSE2
		int bw = n - sw;
	#endif
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

	// Save the needed (i.e. non-unity) twiddles shared by the radix-64 DIF and DIT:
	#ifndef USE_SSE2
		// In scalar-double mode use a const-double array of 7 sets of 14 doubles:
		const double DFT64_TWIDDLES[7][14] = {
			{ 0.0,1.0,ISRT2,ISRT2,-ISRT2,ISRT2,c16,s16,-s16,c16,s16,c16,-c16,s16 },
			{ ISRT2,ISRT2,c16,s16,s16,c16,c32_1,s32_1,s32_3,c32_3,c32_3,s32_3,s32_1,c32_1 },
			{ -ISRT2,ISRT2,s16,c16,-c16,-s16,c32_3,s32_3,-c32_1,s32_1,-s32_1,c32_1,-s32_3,-c32_3 },
			{ c16,s16,c32_1,s32_1,c32_3,s32_3,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7 },
			{ -s16,c16,s32_3,c32_3,-c32_1,s32_1,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3 },
			{ s16,c16,c32_3,s32_3,-s32_1,c32_1,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5 },
			{ -c16,s16,s32_1,c32_1,-s32_3,-c32_3,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1 }
		};
	#endif

		/*   constant index offsets for load/stores are here.	*/
		p01 = NDIVR;
		p02 = p01 +p01;
		p03 = p02 +p01;
		p04 = p03 +p01;
		p05 = p04 +p01;
		p06 = p05 +p01;
		p07 = p06 +p01;
		p08 = p07 +p01;
		p10 = p08 +p08;
		p18 = p10 +p08;
		p20 = p18 +p08;
		p28 = p20 +p08;
		p30 = p28 +p08;
		p38 = p30 +p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p30 = p30 + ( (p30 >> DAT_BITS) << PAD_BITS );
		p38 = p38 + ( (p38 >> DAT_BITS) << PAD_BITS );
	#ifndef USE_SSE2
		p0123[0] = 0; p0123[1] = p01; p0123[2] = p02; p0123[3] = p03;
	#endif
		poff[0x0] =   0; poff[0x1] = p04    ; poff[0x2] = p08; poff[0x3] = p04+p08;
		poff[0x4] = p10; poff[0x5] = p04+p10; poff[0x6] = p18; poff[0x7] = p04+p18;
		poff[0x8] = p20; poff[0x9] = p04+p20; poff[0xa] = p28; poff[0xb] = p04+p28;
		poff[0xc] = p30; poff[0xd] = p04+p30; poff[0xe] = p38; poff[0xf] = p04+p38;
	#if defined(USE_SSE2) && !USE_SCALAR_DFT_MACRO
		// _lin = linear p-multiples (but padding means we can't assume e.g. p02 = 2*p01):
		po_lin[0] =   0; po_lin[1] = p01; po_lin[2] = p02; po_lin[3] = p03;
		po_lin[4] = p04; po_lin[5] = p05; po_lin[6] = p06; po_lin[7] = p07;
	#endif
	#ifndef USE_SSE2
		// _br  = bit-reversed p-multiples:
		po_br[0] =   0; po_br[1] = p04; po_br[2] = p02; po_br[3] = p06;
		po_br[4] = p01; po_br[5] = p05; po_br[6] = p03; po_br[7] = p07;
	#endif

	#if defined(USE_SSE2) && USE_SCALAR_DFT_MACRO
		for(l = 0; l < RADIX; l++) {
			c_offsets[l] = (l<<1);
		}
		// Set array offsets for in/outputs - low parts in fisrt 16 slots, high parts in next 16:
		dft_offsets[0x0] =   0; dft_offsets[0x1] = p01; dft_offsets[0x2] = p02; dft_offsets[0x3] = p03; dft_offsets[0x4] = p04; dft_offsets[0x5] = p05; dft_offsets[0x6] = p06; dft_offsets[0x7] = p07;
		// Now add above low parts to the 7 nonzero high parts:
		l = 0x08;	jp = p08;
		dft_offsets[l+0x0] = jp+  0; dft_offsets[l+0x1] = jp+p01; dft_offsets[l+0x2] = jp+p02; dft_offsets[l+0x3] = jp+p03; dft_offsets[l+0x4] = jp+p04; dft_offsets[l+0x5] = jp+p05; dft_offsets[l+0x6] = jp+p06; dft_offsets[l+0x7] = jp+p07;
		l = 0x10;	jp = p10;
		dft_offsets[l+0x0] = jp+  0; dft_offsets[l+0x1] = jp+p01; dft_offsets[l+0x2] = jp+p02; dft_offsets[l+0x3] = jp+p03; dft_offsets[l+0x4] = jp+p04; dft_offsets[l+0x5] = jp+p05; dft_offsets[l+0x6] = jp+p06; dft_offsets[l+0x7] = jp+p07;
		l = 0x18;	jp = p18;
		dft_offsets[l+0x0] = jp+  0; dft_offsets[l+0x1] = jp+p01; dft_offsets[l+0x2] = jp+p02; dft_offsets[l+0x3] = jp+p03; dft_offsets[l+0x4] = jp+p04; dft_offsets[l+0x5] = jp+p05; dft_offsets[l+0x6] = jp+p06; dft_offsets[l+0x7] = jp+p07;
		l = 0x20;	jp = p20;
		dft_offsets[l+0x0] = jp+  0; dft_offsets[l+0x1] = jp+p01; dft_offsets[l+0x2] = jp+p02; dft_offsets[l+0x3] = jp+p03; dft_offsets[l+0x4] = jp+p04; dft_offsets[l+0x5] = jp+p05; dft_offsets[l+0x6] = jp+p06; dft_offsets[l+0x7] = jp+p07;
		l = 0x28;	jp = p28;
		dft_offsets[l+0x0] = jp+  0; dft_offsets[l+0x1] = jp+p01; dft_offsets[l+0x2] = jp+p02; dft_offsets[l+0x3] = jp+p03; dft_offsets[l+0x4] = jp+p04; dft_offsets[l+0x5] = jp+p05; dft_offsets[l+0x6] = jp+p06; dft_offsets[l+0x7] = jp+p07;
		l = 0x30;	jp = p30;
		dft_offsets[l+0x0] = jp+  0; dft_offsets[l+0x1] = jp+p01; dft_offsets[l+0x2] = jp+p02; dft_offsets[l+0x3] = jp+p03; dft_offsets[l+0x4] = jp+p04; dft_offsets[l+0x5] = jp+p05; dft_offsets[l+0x6] = jp+p06; dft_offsets[l+0x7] = jp+p07;
		l = 0x38;	jp = p38;
		dft_offsets[l+0x0] = jp+  0; dft_offsets[l+0x1] = jp+p01; dft_offsets[l+0x2] = jp+p02; dft_offsets[l+0x3] = jp+p03; dft_offsets[l+0x4] = jp+p04; dft_offsets[l+0x5] = jp+p05; dft_offsets[l+0x6] = jp+p06; dft_offsets[l+0x7] = jp+p07;
	#endif

	#ifdef USE_SSE2
		r00 = thread_arg->r00;	// declared above
		tmp = r00;
		r00 = tmp + 0x00;		r40 = tmp + 0x40;
		r02 = tmp + 0x02;		//r42 = tmp + 0x42;
		r04 = tmp + 0x04;		//r44 = tmp + 0x44;
		r06 = tmp + 0x06;		//r46 = tmp + 0x46;
		r08 = tmp + 0x08;		//r48 = tmp + 0x48;
		r0A = tmp + 0x0a;		//r4A = tmp + 0x4a;
		r0C = tmp + 0x0c;		//r4C = tmp + 0x4c;
		r0E = tmp + 0x0e;		//r4E = tmp + 0x4e;
		r10 = tmp + 0x10;		r50 = tmp + 0x50;
		//r12 = tmp + 0x12;		r52 = tmp + 0x52;
		//r14 = tmp + 0x14;		r54 = tmp + 0x54;
		//r16 = tmp + 0x16;		r56 = tmp + 0x56;
		//r18 = tmp + 0x18;		r58 = tmp + 0x58;
		//r1A = tmp + 0x1a;		r5A = tmp + 0x5a;
		//r1C = tmp + 0x1c;		r5C = tmp + 0x5c;
		//r1E = tmp + 0x1e;		r5E = tmp + 0x5e;
		r20 = tmp + 0x20;		r60 = tmp + 0x60;
		//r22 = tmp + 0x22;		r62 = tmp + 0x62;
		//r24 = tmp + 0x24;		r64 = tmp + 0x64;
		//r26 = tmp + 0x26;		r66 = tmp + 0x66;
		//r28 = tmp + 0x28;		r68 = tmp + 0x68;
		//r2A = tmp + 0x2a;		r6A = tmp + 0x6a;
		//r2C = tmp + 0x2c;		r6C = tmp + 0x6c;
		//r2E = tmp + 0x2e;		r6E = tmp + 0x6e;
		r30 = tmp + 0x30;		r70 = tmp + 0x70;
		//r32 = tmp + 0x32;		r72 = tmp + 0x72;
		//r34 = tmp + 0x34;		r74 = tmp + 0x74;
		//r36 = tmp + 0x36;		r76 = tmp + 0x76;
		//r38 = tmp + 0x38;		r78 = tmp + 0x78;
		//r3A = tmp + 0x3a;		r7A = tmp + 0x7a;
		//r3C = tmp + 0x3c;		r7C = tmp + 0x7c;
		//r3E = tmp + 0x3e;		r7E = tmp + 0x7e;
		tmp += 0x80;
		s1p00 = tmp + 0x00;	s1p20 = tmp + 0x40;
		s1p01 = tmp + 0x02;	//s1p21 = tmp + 0x42;
		s1p02 = tmp + 0x04;	//s1p22 = tmp + 0x44;
		s1p03 = tmp + 0x06;	//s1p23 = tmp + 0x46;
		s1p04 = tmp + 0x08;	//s1p24 = tmp + 0x48;
		s1p05 = tmp + 0x0a;	//s1p25 = tmp + 0x4a;
		s1p06 = tmp + 0x0c;	//s1p26 = tmp + 0x4c;
		s1p07 = tmp + 0x0e;	//s1p27 = tmp + 0x4e;
		s1p08 = tmp + 0x10;	s1p28 = tmp + 0x50;
		//s1p09 = tmp + 0x12;	s1p29 = tmp + 0x52;
		//s1p0a = tmp + 0x14;	s1p2a = tmp + 0x54;
		//s1p0b = tmp + 0x16;	s1p2b = tmp + 0x56;
		//s1p0c = tmp + 0x18;	s1p2c = tmp + 0x58;
		//s1p0d = tmp + 0x1a;	s1p2d = tmp + 0x5a;
		//s1p0e = tmp + 0x1c;	s1p2e = tmp + 0x5c;
		//s1p0f = tmp + 0x1e;	s1p2f = tmp + 0x5e;
		s1p10 = tmp + 0x20;	s1p30 = tmp + 0x60;
		//s1p11 = tmp + 0x22;	s1p31 = tmp + 0x62;
		//s1p12 = tmp + 0x24;	s1p32 = tmp + 0x64;
		//s1p13 = tmp + 0x26;	s1p33 = tmp + 0x66;
		//s1p14 = tmp + 0x28;	s1p34 = tmp + 0x68;
		//s1p15 = tmp + 0x2a;	s1p35 = tmp + 0x6a;
		//s1p16 = tmp + 0x2c;	s1p36 = tmp + 0x6c;
		//s1p17 = tmp + 0x2e;	s1p37 = tmp + 0x6e;
		s1p18 = tmp + 0x30;	s1p38 = tmp + 0x70;
		//s1p19 = tmp + 0x32;	s1p39 = tmp + 0x72;
		//s1p1a = tmp + 0x34;	s1p3a = tmp + 0x74;
		//s1p1b = tmp + 0x36;	s1p3b = tmp + 0x76;
		//s1p1c = tmp + 0x38;	s1p3c = tmp + 0x78;
		//s1p1d = tmp + 0x3a;	s1p3d = tmp + 0x7a;
		//s1p1e = tmp + 0x3c;	s1p3e = tmp + 0x7c;
		//s1p1f = tmp + 0x3e;	s1p3f = tmp + 0x7e;
		tmp += 0x80;
	  #if !USE_SCALAR_DFT_MACRO
		// To support FMA versions of the radix-8 macros used to build radix-64 we insert a standalone copy of the [2,1,sqrt2,isrt2] quartet:
		two     = tmp + 0;	// AVX+ versions of various DFT macros assume consts 2.0,1.0,isrt2 laid out thusly
		//one     = tmp + 1;
		//sqrt2	= tmp + 2;
	//	isrt2   = tmp + 3;	Unnamed slot, since previous layout below already has an iart2 pointer
		tmp += 4;
		nisrt2	= tmp + 0x00;	// For the +- isrt2 pair put the - datum first, thus cc0 satisfies
		 isrt2	= tmp + 0x01;	// the same "cc-1 gets you isrt2" property as do the other +-[cc,ss] pairs.
/*===================
29 Jan 2021: -O3 -march=knl build on KNL hits SIGILL, illegal instruction exception on above pointer-add:
	...
   0x0000000000723126 <+1718>:	vmovdqa %xmm5,0x690(%rsp)
   0x000000000072312f <+1727>:	vmovq  %rsi,%xmm5
   0x0000000000723134 <+1732>:	lea    0x4540(%rbp),%rsi
   0x000000000072313b <+1739>:	vpinsrq $0x1,%rcx,%xmm5,%xmm13
   0x0000000000723141 <+1745>:	vmovq  %rsi,%xmm5
   0x0000000000723146 <+1750>:	lea    0x4580(%rbp),%rsi
=> 0x000000000072314d <+1757>:	vmovdqa64 %xmm16,%xmm6
   0x0000000000723153 <+1763>:	vpinsrq $0x1,%rsi,%xmm5,%xmm5

There are other vmovdqa in the disassembly, but this is the only 64-suffixed one.
This family of instructions is listed as "Move Aligned Packed Integer Values"
Only the 512-bit ZMM-operand version is available in AVX512F; the XMM,YMM forms need AVX512VL.
Workaround: Compiled just this file with -O2, rest with usual -O3.
===================*/
		 cc0	= tmp + 0x02;
		 ss0	= tmp + 0x03;
	// [copy isrt2]	= tmp + 0x04;
		 cc1	= tmp + 0x05;
		 ss1	= tmp + 0x06;
	// [copy isrt2]	= tmp + 0x07;
	    ncc1	= tmp + 0x08;
		//nss1	= tmp + 0x09;
	// [copy isrt2]	= tmp + 0x0a;
		 cc2	= tmp + 0x0b;
		 ss2	= tmp + 0x0c;
	// [copy isrt2]	= tmp + 0x0d;
		ncc2	= tmp + 0x0e;
		nss2	= tmp + 0x0f;
	// [copy isrt2]	= tmp + 0x10;
		 cc3	= tmp + 0x11;
		 ss3	= tmp + 0x12;
	// [copy isrt2]	= tmp + 0x13;
		ncc3	= tmp + 0x14;
		nss3	= tmp + 0x15;
	// [copy isrt2]	= tmp + 0x16;
		 cc4	= tmp + 0x17;
		 ss4	= tmp + 0x18;
	// [copy isrt2]	= tmp + 0x19;
		ncc4	= tmp + 0x1a;
		nss4	= tmp + 0x1b;
	// [copy isrt2]	= tmp + 0x1c;
		 cc5	= tmp + 0x1d;
		 ss5	= tmp + 0x1e;
	// [copy isrt2]	= tmp + 0x1f;
		//ncc5	= tmp + 0x20;
		nss5	= tmp + 0x21;
	// [copy isrt2]	= tmp + 0x22;
		 cc6	= tmp + 0x23;
		 ss6	= tmp + 0x24;
	// [copy isrt2]	= tmp + 0x25;
		ncc6	= tmp + 0x26;
		nss6	= tmp + 0x27;
	// [copy isrt2]	= tmp + 0x28;
		 cc7	= tmp + 0x29;
		 ss7	= tmp + 0x2a;
	// [copy isrt2]	= tmp + 0x2b;
		ncc7	= tmp + 0x2c;
		//nss7	= tmp + 0x2d;
		tmp += 0x2e;	// sc_ptr += 0x132
		// Array of 7*14 SIMD twiddles-pointers for the reduced-#args version of SSE2_RADIX8_DI[F,T]_TWIDDLE_OOP
		vec_dbl *w[98] = {
			ss0,cc0,isrt2,isrt2,nisrt2,isrt2,cc4,ss4,nss4,cc4,ss4,cc4,ncc4,ss4,
			isrt2,isrt2,cc4,ss4,ss4,cc4,cc2,ss2,ss6,cc6,cc6,ss6,ss2,cc2,
			nisrt2,isrt2,ss4,cc4,ncc4,nss4,cc6,ss6,ncc2,ss2,nss2,cc2,nss6,ncc6,
			cc4,ss4,cc2,ss2,cc6,ss6,cc1,ss1,cc5,ss5,cc3,ss3,cc7,ss7,
			nss4,cc4,ss6,cc6,ncc2,ss2,cc5,ss5,ncc7,ss7,ss1,cc1,ncc3,nss3,
			ss4,cc4,cc6,ss6,nss2,cc2,cc3,ss3,ss1,cc1,ss7,cc7,nss5,cc5,
			ncc4,ss4,ss2,cc2,nss6,ncc6,cc7,ss7,ncc3,nss3,nss5,cc5,ss1,ncc1
		}, **twid_ptrs;
	  #endif
	  #ifdef USE_AVX512
		cy_r = tmp;	cy_i = tmp+0x08;	tmp += 2*0x08;	// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		//sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod
	  #elif defined(USE_AVX)
		cy_r = tmp;	cy_i = tmp+0x10;	tmp += 2*0x10;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 340 vec_dbl
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 68 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod
	  #else
		cy_r = tmp;	/* cy_i = tmp+0x20; */	tmp += 2*0x20;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 372 complex
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, 2 for Fermat-mod */
	  #endif

		ASSERT((r00 == thread_arg->r00), "thread-local memcheck failed!");
		ASSERT((half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
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

		sign_mask = (uint64*)(r00 + radix64_creals_in_local_store);
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
		#include "radix64_main_carry_loop.h"

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
