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
#include "radix128.h"

#define RADIX 128	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

#define EPS 1e-10

#ifndef COMPACT_OBJ	// Toggle for parametrized-loop-DFT compact-object code scheme
	#define COMPACT_OBJ	1	// This radix large enough that we want small-obj to be default in all build modes
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
  // For Fermat-mod we use RADIX*4 = 512 [note there is no Fermat-mod LOACC option for this power-of-2 DFT] more
  // slots in AVX mode for the compact negacyclic-roots chained-multiply scheme. Add larger of the 2 numbers -
  // 512 for AVX, 40 for SSE2 - to (half_arr_offset128 + RADIX) to get AVX value of radix128_creals_in_local_store:
  #if COMPACT_OBJ
	#define NEXTRA	0xf0
  #else
	#define NEXTRA	0
  #endif

  #ifdef USE_AVX512
	const int half_arr_offset128 = 0x288 + NEXTRA;	// 0x20 = 2*(RADIX/8) fewer cy-slots than in AVX mode
	const int radix128_creals_in_local_store = 0x508 + NEXTRA;	// May not need this many slots, but just use AVX-alloc for now
  #elif defined(USE_AVX)
	const int half_arr_offset128 = 0x2a8 + NEXTRA;	// Used for thread local-storage-integrity checking
	// GCC issues bogus 'error: initializer element is not constant' if try to eval this vis the summation in the RHS comment, so literalize it:
	const int radix128_creals_in_local_store = 0x528 + NEXTRA;	// half_arr_offset128 + RADIX + 0x200
  #else
	const int half_arr_offset128 = 0x2e8 + NEXTRA;	// Used for thread local-storage-integrity checking
	const int radix128_creals_in_local_store = 0x390 + NEXTRA;	// half_arr_offset128 + RADIX + 0x28
  #endif

  #ifdef USE_AVX
	const uint64 radix128_avx_negadwt_consts[RADIX] = {	// 8 entries per line ==> RADIX/8 lines:
		0x3FF0000000000000ull,0x3FEFFF62169B92DBull,0x3FEFFD886084CD0Dull,0x3FEFFA72EFFEF75Dull,0x3FEFF621E3796D7Eull,0x3FEFF095658E71ADull,0x3FEFE9CDAD01883Aull,0x3FEFE1CAFCBD5B09ull,
		0x3FEFD88DA3D12526ull,0x3FEFCE15FD6DA67Bull,0x3FEFC26470E19FD3ull,0x3FEFB5797195D741ull,0x3FEFA7557F08A517ull,0x3FEF97F924C9099Bull,0x3FEF8764FA714BA9ull,0x3FEF7599A3A12077ull,
		0x3FEF6297CFF75CB0ull,0x3FEF4E603B0B2F2Dull,0x3FEF38F3AC64E589ull,0x3FEF2252F7763ADAull,0x3FEF0A7EFB9230D7ull,0x3FEEF178A3E473C2ull,0x3FEED740E7684963ull,0x3FEEBBD8C8DF0B74ull,
		0x3FEE9F4156C62DDAull,0x3FEE817BAB4CD10Dull,0x3FEE6288EC48E112ull,0x3FEE426A4B2BC17Eull,0x3FEE212104F686E5ull,0x3FEDFEAE622DBE2Bull,0x3FEDDB13B6CCC23Cull,0x3FEDB6526238A09Bull,
		0x3FED906BCF328D46ull,0x3FED696173C9E68Bull,0x3FED4134D14DC93Aull,0x3FED17E7743E35DCull,0x3FECED7AF43CC773ull,0x3FECC1F0F3FCFC5Cull,0x3FEC954B213411F5ull,0x3FEC678B3488739Bull,
		0x3FEC38B2F180BDB1ull,0x3FEC08C426725549ull,0x3FEBD7C0AC6F952Aull,0x3FEBA5AA673590D2ull,0x3FEB728345196E3Eull,0x3FEB3E4D3EF55712ull,0x3FEB090A58150200ull,0x3FEAD2BC9E21D511ull,
		0x3FEA9B66290EA1A3ull,0x3FEA63091B02FAE2ull,0x3FEA29A7A0462782ull,0x3FE9EF43EF29AF94ull,0x3FE9B3E047F38741ull,0x3FE9777EF4C7D742ull,0x3FE93A22499263FBull,0x3FE8FBCCA3EF940Dull,
		0x3FE8BC806B151741ull,0x3FE87C400FBA2EBFull,0x3FE83B0E0BFF976Eull,0x3FE7F8ECE3571771ull,0x3FE7B5DF226AAFAFull,0x3FE771E75F037261ull,0x3FE72D0837EFFF96ull,0x3FE6E74454EAA8AFull,
		0x3FE6A09E667F3BCDull,0x3FE6591925F0783Dull,0x3FE610B7551D2CDFull,0x3FE5C77BBE65018Cull,0x3FE57D69348CECA0ull,0x3FE5328292A35596ull,0x3FE4E6CABBE3E5E9ull,0x3FE49A449B9B0939ull,
		0x3FE44CF325091DD6ull,0x3FE3FED9534556D4ull,0x3FE3AFFA292050B9ull,0x3FE36058B10659F3ull,0x3FE30FF7FCE17035ull,0x3FE2BEDB25FAF3EAull,0x3FE26D054CDD12DFull,0x3FE21A799933EB59ull,
		0x3FE1C73B39AE68C8ull,0x3FE1734D63DEDB49ull,0x3FE11EB3541B4B23ull,0x3FE0C9704D5D898Full,0x3FE073879922FFEEull,0x3FE01CFC874C3EB7ull,0x3FDF8BA4DBF89ABAull,0x3FDEDC1952EF78D6ull,
		0x3FDE2B5D3806F63Bull,0x3FDD79775B86E389ull,0x3FDCC66E9931C45Eull,0x3FDC1249D8011EE7ull,0x3FDB5D1009E15CC0ull,0x3FDAA6C82B6D3FCAull,0x3FD9EF7943A8ED8Aull,0x3FD9372A63BC93D7ull,
		0x3FD87DE2A6AEA963ull,0x3FD7C3A9311DCCE7ull,0x3FD7088530FA459Full,0x3FD64C7DDD3F27C6ull,0x3FD58F9A75AB1FDDull,0x3FD4D1E24278E76Aull,0x3FD4135C94176601ull,0x3FD35410C2E18152ull,
		0x3FD294062ED59F06ull,0x3FD1D3443F4CDB3Eull,0x3FD111D262B1F677ull,0x3FD04FB80E37FDAEull,0x3FCF19F97B215F1Bull,0x3FCD934FE5454311ull,0x3FCC0B826A7E4F63ull,0x3FCA82A025B00451ull,
		0x3FC8F8B83C69A60Bull,0x3FC76DD9DE50BF31ull,0x3FC5E214448B3FC6ull,0x3FC45576B1293E5Aull,0x3FC2C8106E8E613Aull,0x3FC139F0CEDAF577ull,0x3FBF564E56A9730Eull,0x3FBC3785C79EC2D5ull,
		0x3FB917A6BC29B42Cull,0x3FB5F6D00A9AA419ull,0x3FB2D52092CE19F6ull,0x3FAF656E79F820E0ull,0x3FA91F65F10DD814ull,0x3FA2D865759455CDull,0x3F992155F7A3667Eull,0x3F8921D1FCDEC784ull
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

int radix128_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-128 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-128 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix16_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
#if !defined(MULTITHREAD) || defined(USE_SSE2)
	#include "radix128_twiddles.h"
#endif
	const char func[] = "radix128_ditN_cy_dif1";
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
	int j2,jp,ntmp;
  #endif
  #ifndef MULTITHREAD
	int jstart,jhi;
  #endif
  #ifdef USE_SSE2
	int nbytes;
  #endif
  #if !defined(MULTITHREAD) && defined(USE_SSE2)
	int incr;
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 64|128|256 for avx512,avx,sse, respectively:
	const int incr_long = 8,incr_short = 4;
	if(USE_SHORT_CY_CHAIN)	// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
		incr = incr_short;
	else
		incr = incr_long;
  #endif
  #ifndef MULTITHREAD
	#if !defined(USE_SSE2) || defined(USE_AVX)
	int k1,k2;
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
  #endif
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
	static uint32 bw,sw,p01,p02,p03,p04,p05,p06,p07,p08,p09,p0a,p0b,p0c,p0d,p0e,p0f,p10,p18,p20,p28,p30,p38,p40,p48,p50,p58,p60,p68,p70,p78, nsave = 0;
  #if !defined(MULTITHREAD) || defined(USE_SSE2)
	static uint32 nm1;
  #endif
	static int poff[RADIX>>2];	// Store mults of p-offsets for loop-controlled DFT macro calls
  #if 0 && !defined(MULTITHREAD)
	static int po_lin[16];		// Store mults of p-offsets for loop-controlled DFT macro calls
  #endif
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	static int po_br[16];		// Store mults of p-offsets for loop-controlled DFT macro calls
  #endif
	static double radix_inv, n2inv;
	double scale, dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
  #ifdef USE_SSE2
	const double *addr, *addi;
  #endif
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

//	int i0,i1,m0,m1,m3;	/* m2 already def'd for regular carry sequence */
	double *add0,*add1,*add2,*add3/* ,*add4,*add5,*add6,*add7 */;	/* Addresses into array sections */

  #endif
	// Uint64 bitmaps for alternate "rounded the other way" copies of sqrt2,isrt2. Default round-to-nearest versions
	// (SQRT2, ISRT2) end in ...3BCD. Since we round these down as ...3BCC90... --> ..3BCC, append _dn to varnames:
	const uint64 sqrt2_dn = 0x3FF6A09E667F3BCCull, isrt2_dn = 0x3FE6A09E667F3BCCull;
  #if !defined(MULTITHREAD) || !defined(USE_SSE2)
	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
  #endif
  #ifndef USE_AVX512
	const double crnd = 3.0*0x4000000*0x2000000;
  #endif
  #ifndef USE_AVX
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
  #endif
	static vec_dbl *max_err, *sse2_rnd, *half_arr, *two,*one,*sqrt2,
	#ifndef MULTITHREAD
		*isrt2,
	#endif
	#if !COMPACT_OBJ
		*cc0,*ss0,
	#endif
	#if COMPACT_OBJ
		// ptrs to 16 sets of twiddles shared by the 2nd-half DIF and DIT DFT macros:
		*twid0,
	#else
		*nisrt2,	// each non-unity root now needs a negated counterpart:
		*cc128_1,*ss128_1,*ncc128_1,*nss128_1,*cc64_1,*ss64_1,*ncc64_1,*nss64_1,*cc128_3,*ss128_3,*ncc128_3,*nss128_3,
		*cc32_1,*ss32_1,*ncc32_1,*nss32_1,*cc128_5,*ss128_5,*ncc128_5,*nss128_5,*cc64_3,*ss64_3,*ncc64_3,*nss64_3,*cc128_7,*ss128_7,*ncc128_7,*nss128_7,
		*cc16,*ss16,*ncc16,*nss16,*cc128_9,*ss128_9,*ncc128_9,*nss128_9,*cc64_5,*ss64_5,*ncc64_5,*nss64_5,*cc128_b,*ss128_b,*ncc128_b,*nss128_b,
		*cc32_3,*ss32_3,*ncc32_3,*nss32_3,*cc128_d,*ss128_d,*ncc128_d,*nss128_d,*cc64_7,*ss64_7,*ncc64_7,*nss64_7,*cc128_f,*ss128_f,*ncc128_f,*nss128_f,
	#endif
		*r00,//*r10,*r20,*r30,*r40,*r50,*r60,*r70,
	#ifndef MULTITHREAD
		*s1p00,
	#endif
		//*s1p10,*s1p20,*s1p30,*s1p40,*s1p50,*s1p60,*s1p70,
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
	static task_control_t   task_control = {NULL, (void*)cy128_process_chunk, NULL, 0x0};

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
		cslots_in_local_store = radix128_creals_in_local_store + (20+RADIX/2)/2;	// Just add enough int64 space for both cases, plus some
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix128_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = sc_ptr;
		r00 = tmp + 0x00;	//r40 = tmp + 0x80;
		//r10 = tmp + 0x20;	r50 = tmp + 0xa0;
		//r20 = tmp + 0x40;	r60 = tmp + 0xc0;
		//r30 = tmp + 0x60;	r70 = tmp + 0xe0;
		tmp += 0x100;
	  #ifndef MULTITHREAD
		s1p00 = tmp + 0x00;	//s1p40 = tmp + 0x80;
	  #endif
		//s1p10 = tmp + 0x20;	s1p50 = tmp + 0xa0;
		//s1p20 = tmp + 0x40;	s1p60 = tmp + 0xc0;
		//s1p30 = tmp + 0x60;	s1p70 = tmp + 0xe0;
		tmp += 0x100;
		two     = tmp + 0;	// AVX+ versions of radix-8,16,32 twiddleless-DFT macros need consts [2,1,sqrt2,isrt2] quartet laid out thusly
		one     = tmp + 1;
		sqrt2   = tmp + 2;
		// Next 3 slots unnamed in non-compact-obj-code mode, since that layout below already has iart2,c16,s16 pointers
	  #ifndef MULTITHREAD
		isrt2   = tmp + 3;
	  #endif
		//cc0		= tmp + 4;	// ...and the radix-16 DFT macrod need the base 16-root [re,im] pair right above isrt2.
		//ss0		= tmp + 5;
		tmp += 6;
		// Each non-unity root now needs a negated counterpart:
		/* Stupidity: Since a truly general-purpose [in the sense that it can be used for our radix-128 internal-twiddles]
		radix-8 DFT-with-twiddles macro needs 8 in-addresses [corr. to the 8 real parts of the input data], 8 o-addresses,
		and 7 each of cosine and sine data [which cannot be assumed to occur in fixed-stride pairs - cf. our usage of
		SSE2_RADIX8_DIT_TWIDDLE_OOP() below], that hits the GCC hard limit of 30-operands for ASM macros, but we still
		need one more operand for the ISRT2 pointer. Only easy workaround I found for this is to stick a vector-ISRT2 copy
		in between each +-[cc,ss] vector-data pair, thus any time we need a vector-isrt2 for the radix-8 internal twiddles
		we get it at (vec_dbl*)cc-1.	/Stupidity */
	  #if COMPACT_OBJ
	   #ifndef MULTITHREAD
		isrt2   = tmp - 3;	// Alternate layout below already has [isrt2,cc0,ss0] pointer-triplet in a different place, leave these 3 slots unnamed there
	   #endif
		//cc0		= tmp - 2;
		//ss0		= tmp - 1;
		// ptrs to 16 sets (3*7 = 21 vec_dbl data each) of non-unity twiddles shared by the 2nd-half DIF and DIT DFT macros:
		twid0  = tmp + 0x00;
		tmp += 0x150;	// += 16*21 = 0x150 => sc_ptr + 0x356 (alternate version below incurs 0x60 elts, thus diff = 0xf0
	  #else
		// For the +- isrt2 pair put the - datum first, thus cc0 satisfies the same "cc-1 gets you isrt2" property as do the other +-[cc,ss] pairs.
									// [copy isrt2]	= tmp + 0x2f;
		nisrt2  = tmp + 0x00;				 cc16   = tmp + 0x30;
		 isrt2  = tmp + 0x01;				 ss16   = tmp + 0x31;// [copy isrt2]	= tmp + 0x32;
		 cc0    = tmp + 0x02;				ncc16   = tmp + 0x32;
		 ss0    = tmp + 0x03;				nss16   = tmp + 0x33;
// [copy isrt2]	= tmp + 0x04;		// [copy isrt2]	= tmp + 0x34;
		 cc128_1= tmp + 0x05;				 cc128_9= tmp + 0x35;
		 ss128_1= tmp + 0x06;				 ss128_9= tmp + 0x36;
// [copy isrt2]	= tmp + 0x07;		// [copy isrt2] = tmp + 0x37;
		ncc128_1= tmp + 0x08;				ncc128_9= tmp + 0x38;
		nss128_1= tmp + 0x09;				nss128_9= tmp + 0x39;
// [copy isrt2]	= tmp + 0x0a;		// [copy isrt2]	= tmp + 0x3a;
		 cc64_1	= tmp + 0x0b;				 cc64_5	= tmp + 0x3b;
		 ss64_1	= tmp + 0x0c;				 ss64_5	= tmp + 0x3c;
// [copy isrt2]	= tmp + 0x0d;		// [copy isrt2]	= tmp + 0x3d;
		ncc64_1	= tmp + 0x0e;				ncc64_5	= tmp + 0x3e;
		nss64_1	= tmp + 0x0f;				nss64_5	= tmp + 0x3f;
// [copy isrt2]	= tmp + 0x10;		// [copy isrt2]	= tmp + 0x40;
		 cc128_3= tmp + 0x11;				 cc128_b= tmp + 0x41;
		 ss128_3= tmp + 0x12;				 ss128_b= tmp + 0x42;
// [copy isrt2]	= tmp + 0x13;		// [copy isrt2] = tmp + 0x43;
		ncc128_3= tmp + 0x14;				ncc128_b= tmp + 0x44;
		nss128_3= tmp + 0x15;				nss128_b= tmp + 0x45;
// [copy isrt2]	= tmp + 0x16;		// [copy isrt2]	= tmp + 0x46;
		 cc32_1	= tmp + 0x17;				 cc32_3	= tmp + 0x47;
		 ss32_1	= tmp + 0x18;				 ss32_3	= tmp + 0x48;
// [copy isrt2]	= tmp + 0x19;		// [copy isrt2]	= tmp + 0x49;
		ncc32_1	= tmp + 0x1a;				ncc32_3	= tmp + 0x4a;
		nss32_1	= tmp + 0x1b;				nss32_3	= tmp + 0x4b;
// [copy isrt2]	= tmp + 0x1c;		// [copy isrt2]	= tmp + 0x4c;
		 cc128_5= tmp + 0x1d;				 cc128_d= tmp + 0x4d;
		 ss128_5= tmp + 0x1e;				 ss128_d= tmp + 0x4e;
// [copy isrt2]	= tmp + 0x1f;		// [copy isrt2] = tmp + 0x4f;
		ncc128_5= tmp + 0x20;				ncc128_d= tmp + 0x50;
		nss128_5= tmp + 0x21;				nss128_d= tmp + 0x51;
// [copy isrt2]	= tmp + 0x22;		// [copy isrt2]	= tmp + 0x52;
		 cc64_3	= tmp + 0x23;				 cc64_7	= tmp + 0x53;
		 ss64_3	= tmp + 0x24;				 ss64_7	= tmp + 0x54;
// [copy isrt2]	= tmp + 0x25;		// [copy isrt2]	= tmp + 0x55;
		ncc64_3	= tmp + 0x26;				ncc64_7	= tmp + 0x56;
		nss64_3	= tmp + 0x27;				nss64_7	= tmp + 0x57;
// [copy isrt2]	= tmp + 0x28;		// [copy isrt2]	= tmp + 0x58;
		 cc128_7= tmp + 0x29;				 cc128_f= tmp + 0x59;
		 ss128_7= tmp + 0x2a;				 ss128_f= tmp + 0x5a;
// [copy isrt2]	= tmp + 0x2b;		// [copy isrt2] = tmp + 0x5b;
		ncc128_7= tmp + 0x2c;				ncc128_f= tmp + 0x5c;
		nss128_7= tmp + 0x2d;				nss128_f= tmp + 0x5d;
// [copy isrt2]	= tmp + 0x2e;		// [copy isrt2]	= tmp + 0x5e;
		tmp += 0x60;	// Add pad slot to make offset even; sc_ptr += 0x266
	  #endif
	  #ifdef USE_AVX512
		cy_r = tmp;										// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
	   #ifndef MULTITHREAD
					cy_i = tmp+0x10;					// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
	   #endif
										tmp += 2*0x10;	// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
		cy_r = tmp;										// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
	   #ifndef MULTITHREAD
					cy_i = tmp+0x20;					// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
	   #endif
										tmp += 2*0x20;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 0x3f8 vec_dbl if COMPACT_OBJ = true, 0x2a8 if not.
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 96 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy_r = tmp;	/* cy_i = tmp+0x40; */	tmp += 2*0x40;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 0x438 vec_dbl
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 32 x 16 bytes for Mersenne-mod, 2 for Fermat-mod */
	  #endif
//		ASSERT(half_arr_offset == (uint32)(half_arr-sc_ptr), "half_arr_offset mismatches actual!");
		ASSERT((radix128_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r00) + (20 << L2_SZ_VD), "radix128_creals_in_local_store checksum failed!");

		/* These remain fixed: */
		VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one  , 1.0  );
		tmp = sqrt2+1;
	  #if 1
		// Use alternate "rounded the other way" copies of sqrt2,isrt2 for greater accuracy:
		dtmp = *(double *)&sqrt2_dn;	VEC_DBL_INIT(sqrt2, dtmp);
		dtmp = *(double *)&isrt2_dn;	VEC_DBL_INIT(tmp  , dtmp);	// This copy of isrt2 only gets assigned the same-name ptr in compact-obj_code mode
	  #else
		VEC_DBL_INIT(sqrt2, SQRT2);		VEC_DBL_INIT(tmp  , ISRT2);	// This copy of isrt2 only gets assigned the same-name ptr in compact-obj_code mode
	  #endif
		VEC_DBL_INIT(tmp+1, c16  );	VEC_DBL_INIT(tmp+2, s16  );	// Ditto for these 2 ptrs of the needed [isrt2,c16,s16] triplet

	  #if COMPACT_OBJ
		// Use loop-based init for code compactness and extensibility to larger pow2 radices:
		for(l = 0; l < 16; l++) {
			j = reverse(l,4);	// twid0-offsets are processed in BR16 order
			tmp = twid0 + (j<<4)+(j<<2)+j;	// Twid-offsets are multiples of 21 vec_dbl
			addr = DFT128_TWIDDLES[l]; addi = addr+1;	// Pointer to required row of 2-D twiddles array
			VEC_DBL_INIT(tmp, ISRT2); tmp++;	VEC_DBL_INIT(tmp,*(addr+0x00)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x00)); tmp++;
			VEC_DBL_INIT(tmp, ISRT2); tmp++;	VEC_DBL_INIT(tmp,*(addr+0x02)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x02)); tmp++;
			VEC_DBL_INIT(tmp, ISRT2); tmp++;	VEC_DBL_INIT(tmp,*(addr+0x04)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x04)); tmp++;
			VEC_DBL_INIT(tmp, ISRT2); tmp++;	VEC_DBL_INIT(tmp,*(addr+0x06)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x06)); tmp++;
			VEC_DBL_INIT(tmp, ISRT2); tmp++;	VEC_DBL_INIT(tmp,*(addr+0x08)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x08)); tmp++;
			VEC_DBL_INIT(tmp, ISRT2); tmp++;	VEC_DBL_INIT(tmp,*(addr+0x0a)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x0a)); tmp++;
			VEC_DBL_INIT(tmp, ISRT2); tmp++;	VEC_DBL_INIT(tmp,*(addr+0x0c)); tmp++;	VEC_DBL_INIT(tmp,*(addi+0x0c)); tmp++;
		}
	  #else
		VEC_DBL_INIT(nisrt2,-ISRT2);
		VEC_DBL_INIT( isrt2, ISRT2);										// Copies of +ISRT2 needed for 30-asm-macro-operand-GCC-limit workaround:
		VEC_DBL_INIT( cc0    ,    1.0);		VEC_DBL_INIT( ss0    ,    0.0);		tmp =  ss0    +1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc128_1, c128_1);		VEC_DBL_INIT( ss128_1, s128_1);		tmp =  ss128_1+1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc128_1,-c128_1);		VEC_DBL_INIT(nss128_1,-s128_1);		tmp = nss128_1+1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc64_1 , c64_1 );		VEC_DBL_INIT( ss64_1 , s64_1 );		tmp =  ss64_1 +1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc64_1 ,-c64_1 );		VEC_DBL_INIT(nss64_1 ,-s64_1 );		tmp = nss64_1 +1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc128_3, c128_3);		VEC_DBL_INIT( ss128_3, s128_3);		tmp =  ss128_3+1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc128_3,-c128_3);		VEC_DBL_INIT(nss128_3,-s128_3);		tmp = nss128_3+1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc32_1 , c32_1 );		VEC_DBL_INIT( ss32_1 , s32_1 );		tmp =  ss32_1 +1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc32_1 ,-c32_1 );		VEC_DBL_INIT(nss32_1 ,-s32_1 );		tmp = nss32_1 +1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc128_5, c128_5);		VEC_DBL_INIT( ss128_5, s128_5);		tmp =  ss128_5+1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc128_5,-c128_5);		VEC_DBL_INIT(nss128_5,-s128_5);		tmp = nss128_5+1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc64_3 , c64_3 );		VEC_DBL_INIT( ss64_3 , s64_3 );		tmp =  ss64_3 +1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc64_3 ,-c64_3 );		VEC_DBL_INIT(nss64_3 ,-s64_3 );		tmp = nss64_3 +1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc128_7, c128_7);		VEC_DBL_INIT( ss128_7, s128_7);		tmp =  ss128_7+1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc128_7,-c128_7);		VEC_DBL_INIT(nss128_7,-s128_7);		tmp = nss128_7+1; VEC_DBL_INIT(tmp, ISRT2);
		tmp =  cc16-1; VEC_DBL_INIT(tmp, ISRT2);	// Special init for start of rcol sincos data, needed for isrt/cc16/ss16 triplet expected by radix-16 DFT macros
		VEC_DBL_INIT( cc16   , c16   );		VEC_DBL_INIT( ss16   , s16   );		tmp =  ss16   +1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc16   ,-c16   );		VEC_DBL_INIT(nss16   ,-s16   );		tmp = nss16   +1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc128_9, c128_9);		VEC_DBL_INIT( ss128_9, s128_9);		tmp =  ss128_9+1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc128_9,-c128_9);		VEC_DBL_INIT(nss128_9,-s128_9);		tmp = nss128_9+1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc64_5 , c64_5 );		VEC_DBL_INIT( ss64_5 , s64_5 );		tmp =  ss64_5 +1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc64_5 ,-c64_5 );		VEC_DBL_INIT(nss64_5 ,-s64_5 );		tmp = nss64_5 +1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc128_b, c128_b);		VEC_DBL_INIT( ss128_b, s128_b);		tmp =  ss128_b+1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc128_b,-c128_b);		VEC_DBL_INIT(nss128_b,-s128_b);		tmp = nss128_b+1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc32_3 , c32_3 );		VEC_DBL_INIT( ss32_3 , s32_3 );		tmp =  ss32_3 +1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc32_3 ,-c32_3 );		VEC_DBL_INIT(nss32_3 ,-s32_3 );		tmp = nss32_3 +1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc128_d, c128_d);		VEC_DBL_INIT( ss128_d, s128_d);		tmp =  ss128_d+1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc128_d,-c128_d);		VEC_DBL_INIT(nss128_d,-s128_d);		tmp = nss128_d+1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc64_7 , c64_7 );		VEC_DBL_INIT( ss64_7 , s64_7 );		tmp =  ss64_7 +1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc64_7 ,-c64_7 );		VEC_DBL_INIT(nss64_7 ,-s64_7 );		tmp = nss64_7 +1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc128_f, c128_f);		VEC_DBL_INIT( ss128_f, s128_f);		tmp =  ss128_f+1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc128_f,-c128_f);		VEC_DBL_INIT(nss128_f,-s128_f);		tmp = nss128_f+1; VEC_DBL_INIT(tmp, ISRT2);
	  #endif

		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		sse2_rnd->d0 = base[0]; sse2_rnd->d1 = baseinv[1]; sse2_rnd->d2 = wts_mult[1]; sse2_rnd->d3 = inv_mult[0];
	  #else
		VEC_DBL_INIT(sse2_rnd, crnd);
	  #endif

		// Propagate the above consts to the remaining threads:
		nbytes = (intptr_t)cy_r - (intptr_t)two;	// #bytes in 1st of above block of consts
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
			// Simple qfloat-based loop to crank out the roots which populate the radix128_avx_negadwt_consts table:
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
			tmp64 = radix128_avx_negadwt_consts[1];	tmp->d1 = tm2->d7 = *(double *)&tmp64;
			tmp64 = radix128_avx_negadwt_consts[2];	tmp->d2 = tm2->d6 = *(double *)&tmp64;
			tmp64 = radix128_avx_negadwt_consts[3];	tmp->d3 = tm2->d5 = *(double *)&tmp64;
			tmp64 = radix128_avx_negadwt_consts[4];	tmp->d4 = tm2->d4 = *(double *)&tmp64;
			tmp64 = radix128_avx_negadwt_consts[5];	tmp->d5 = tm2->d3 = *(double *)&tmp64;
			tmp64 = radix128_avx_negadwt_consts[6];	tmp->d6 = tm2->d2 = *(double *)&tmp64;
			tmp64 = radix128_avx_negadwt_consts[7];	tmp->d7 = tm2->d1 = *(double *)&tmp64;	tmp += 2;
			for(j = 8; j < RADIX; j += 8) {
				tmp64 = radix128_avx_negadwt_consts[j+0];	tmp->d0 = tm2->d0 = *(double *)&tmp64;	tm2 -= 2;
				tmp64 = radix128_avx_negadwt_consts[j+1];	tmp->d1 = tm2->d7 = *(double *)&tmp64;
				tmp64 = radix128_avx_negadwt_consts[j+2];	tmp->d2 = tm2->d6 = *(double *)&tmp64;
				tmp64 = radix128_avx_negadwt_consts[j+3];	tmp->d3 = tm2->d5 = *(double *)&tmp64;
				tmp64 = radix128_avx_negadwt_consts[j+4];	tmp->d4 = tm2->d4 = *(double *)&tmp64;
				tmp64 = radix128_avx_negadwt_consts[j+5];	tmp->d5 = tm2->d3 = *(double *)&tmp64;
				tmp64 = radix128_avx_negadwt_consts[j+6];	tmp->d6 = tm2->d2 = *(double *)&tmp64;
				tmp64 = radix128_avx_negadwt_consts[j+7];	tmp->d7 = tm2->d1 = *(double *)&tmp64;	tmp += 2;
			}
			tmp = base_negacyclic_root + RADIX*2;	// reset to point to start of above block

		  #else

			tmp = base_negacyclic_root + RADIX*2;	// First 2*RADIX slots reserved for RADIX/4 copies of the Re/Im parts of the 4 base multipliers
			tm2 = tmp + RADIX/2 - 1;
			// First elt-pair needs special handling - have the 1.0 in avx_negadwt_consts[0] but the sine term buggers things
			tmp->d0 = 1.0;	(tmp+1)->d0 = 0.0;
			tmp64 = radix128_avx_negadwt_consts[1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(  1*I*Pi/(2*RADIX)) = sin((RADIX-  1)*I*Pi/(2*RADIX)) */
			tmp64 = radix128_avx_negadwt_consts[2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(  2*I*Pi/(2*RADIX)) = sin((RADIX-  2)*I*Pi/(2*RADIX)) */
			tmp64 = radix128_avx_negadwt_consts[3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(  3*I*Pi/(2*RADIX)) = sin((RADIX-  3)*I*Pi/(2*RADIX)) */	tmp += 2;
			for(j = 4; j < RADIX; j += 4) {
				tmp64 = radix128_avx_negadwt_consts[j+0];	tmp->d0 = tm2->d0 = *(double *)&tmp64;	tm2 -= 2;
				tmp64 = radix128_avx_negadwt_consts[j+1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;
				tmp64 = radix128_avx_negadwt_consts[j+2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;
				tmp64 = radix128_avx_negadwt_consts[j+3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	tmp += 2;
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
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p0a = p09 + p01;
		p0b = p0a + p01;
		p0c = p0b + p01;
		p0d = p0c + p01;
		p0e = p0d + p01;
		p0f = p0e + p01;
		p10 = p08 + p08;
		p18 = p10 + p08;
		p20 = p18 + p08;
		p28 = p20 + p08;
		p30 = p28 + p08;
		p38 = p30 + p08;
		p40 = p38 + p08;
		p48 = p40 + p08;
		p50 = p48 + p08;
		p58 = p50 + p08;
		p60 = p58 + p08;
		p68 = p60 + p08;
		p70 = p68 + p08;
		p78 = p70 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p0a = p0a + ( (p0a >> DAT_BITS) << PAD_BITS );
		p0b = p0b + ( (p0b >> DAT_BITS) << PAD_BITS );
		p0c = p0c + ( (p0c >> DAT_BITS) << PAD_BITS );
		p0d = p0d + ( (p0d >> DAT_BITS) << PAD_BITS );
		p0e = p0e + ( (p0e >> DAT_BITS) << PAD_BITS );
		p0f = p0f + ( (p0f >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p30 = p30 + ( (p30 >> DAT_BITS) << PAD_BITS );
		p38 = p38 + ( (p38 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p50 = p50 + ( (p50 >> DAT_BITS) << PAD_BITS );
		p58 = p58 + ( (p58 >> DAT_BITS) << PAD_BITS );
		p60 = p60 + ( (p60 >> DAT_BITS) << PAD_BITS );
		p68 = p68 + ( (p68 >> DAT_BITS) << PAD_BITS );
		p70 = p70 + ( (p70 >> DAT_BITS) << PAD_BITS );
		p78 = p78 + ( (p78 >> DAT_BITS) << PAD_BITS );
	#if !defined(USE_SSE2) && !defined(MULTITHREAD)
		p0123[0] = 0; p0123[1] = p01; p0123[2] = p02; p0123[3] = p03;
	#endif
		poff[0x00] =   0; poff[0x01] = p04    ; poff[0x02] = p08; poff[0x03] = p04+p08;
		poff[0x04] = p10; poff[0x05] = p04+p10; poff[0x06] = p18; poff[0x07] = p04+p18;
		poff[0x08] = p20; poff[0x09] = p04+p20; poff[0x0a] = p28; poff[0x0b] = p04+p28;
		poff[0x0c] = p30; poff[0x0d] = p04+p30; poff[0x0e] = p38; poff[0x0f] = p04+p38;
		poff[0x10] = p40; poff[0x11] = p04+p40; poff[0x12] = p48; poff[0x13] = p04+p48;
		poff[0x14] = p50; poff[0x15] = p04+p50; poff[0x16] = p58; poff[0x17] = p04+p58;
		poff[0x18] = p60; poff[0x19] = p04+p60; poff[0x1a] = p68; poff[0x1b] = p04+p68;
		poff[0x1c] = p70; poff[0x1d] = p04+p70; poff[0x1e] = p78; poff[0x1f] = p04+p78;
	#if 0 && !defined(MULTITHREAD)
		// _lin = linear p-multiples (but padding means we can't assume e.g. p02 = 2*p01):
		po_lin[0x0] =   0; po_lin[0x1] = p01; po_lin[0x2] = p02; po_lin[0x3] = p03;
		po_lin[0x4] = p04; po_lin[0x5] = p05; po_lin[0x6] = p06; po_lin[0x7] = p07;
		po_lin[0x8] = p08; po_lin[0x9] = p09; po_lin[0xa] = p0a; po_lin[0xb] = p0b;
		po_lin[0xc] = p0c; po_lin[0xd] = p0d; po_lin[0xe] = p0e; po_lin[0xf] = p0f;
	#endif
	#if !defined(USE_SSE2) && !defined(MULTITHREAD)
		// _br  = bit-reversed p-multiples:
		po_br[0x0] =   0; po_br[0x1] = p08; po_br[0x2] = p04; po_br[0x3] = p0c;
		po_br[0x4] = p02; po_br[0x5] = p0a; po_br[0x6] = p06; po_br[0x7] = p0e;
		po_br[0x8] = p01; po_br[0x9] = p09; po_br[0xa] = p05; po_br[0xb] = p0d;
		po_br[0xc] = p03; po_br[0xd] = p0b; po_br[0xe] = p07; po_br[0xf] = p0f;
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

/*...The radix-128 final DIT pass is here.	*/

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
		#include "radix128_main_carry_loop.h"

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
		ASSERT(0x0 == cy128_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
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

void radix128_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-128 complex DIF FFT pass on the data in the length-N real vector A.
!
!   The data are stored in a 1-D zero-offset array, with 2^PAD_BITS 8-byte padding elements inserted
!   between every block of 2^DAT_BITS contiguous data. The array padding is to prevent data being accessed
!   in strides that are large powers of two and thus to minimize cache thrashing
!   (in cache-based microprocessor architectures) or bank conflicts (in supercomputers.)
!
!   See the documentation in radix64_dif_pass for further details on storage, indexing and DFT decomposition.
*/
	int j,j1,j2,jt,jp;
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p09,p0a,p0b,p0c,p0d,p0e,p0f,p10,p18,p20,p28,p30,p38,p40,p48,p50,p58,p60,p68,p70,p78, first_entry=TRUE;
#if USE_SCALAR_DFT_MACRO
	static int i_offsets[RADIX];
#else
	// We name the even-order roots here the same as we do in the radix-64 code - that allows us to take the
	// radix-8 DFT calls using all-even-order-radix-128 roots more or less straight from the radix-64 code.
	double
		t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,t0fr,t0fi,
		t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,t1fr,t1fi,
		t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,t2fr,t2fi,
		t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei,t3fr,t3fi,
		t40r,t40i,t41r,t41i,t42r,t42i,t43r,t43i,t44r,t44i,t45r,t45i,t46r,t46i,t47r,t47i,t48r,t48i,t49r,t49i,t4ar,t4ai,t4br,t4bi,t4cr,t4ci,t4dr,t4di,t4er,t4ei,t4fr,t4fi,
		t50r,t50i,t51r,t51i,t52r,t52i,t53r,t53i,t54r,t54i,t55r,t55i,t56r,t56i,t57r,t57i,t58r,t58i,t59r,t59i,t5ar,t5ai,t5br,t5bi,t5cr,t5ci,t5dr,t5di,t5er,t5ei,t5fr,t5fi,
		t60r,t60i,t61r,t61i,t62r,t62i,t63r,t63i,t64r,t64i,t65r,t65i,t66r,t66i,t67r,t67i,t68r,t68i,t69r,t69i,t6ar,t6ai,t6br,t6bi,t6cr,t6ci,t6dr,t6di,t6er,t6ei,t6fr,t6fi,
		t70r,t70i,t71r,t71i,t72r,t72i,t73r,t73i,t74r,t74i,t75r,t75i,t76r,t76i,t77r,t77i,t78r,t78i,t79r,t79i,t7ar,t7ai,t7br,t7bi,t7cr,t7ci,t7dr,t7di,t7er,t7ei,t7fr,t7fi;
#endif

	if(!first_entry && (n >> 7) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR = n >> 7;

		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p0a = p09 + p01;
		p0b = p0a + p01;
		p0c = p0b + p01;
		p0d = p0c + p01;
		p0e = p0d + p01;
		p0f = p0e + p01;
		p10 = p08 + p08;
		p18 = p10 + p08;
		p20 = p18 + p08;
		p28 = p20 + p08;
		p30 = p28 + p08;
		p38 = p30 + p08;
		p40 = p38 + p08;
		p48 = p40 + p08;
		p50 = p48 + p08;
		p58 = p50 + p08;
		p60 = p58 + p08;
		p68 = p60 + p08;
		p70 = p68 + p08;
		p78 = p70 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p0a = p0a + ( (p0a >> DAT_BITS) << PAD_BITS );
		p0b = p0b + ( (p0b >> DAT_BITS) << PAD_BITS );
		p0c = p0c + ( (p0c >> DAT_BITS) << PAD_BITS );
		p0d = p0d + ( (p0d >> DAT_BITS) << PAD_BITS );
		p0e = p0e + ( (p0e >> DAT_BITS) << PAD_BITS );
		p0f = p0f + ( (p0f >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p30 = p30 + ( (p30 >> DAT_BITS) << PAD_BITS );
		p38 = p38 + ( (p38 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p50 = p50 + ( (p50 >> DAT_BITS) << PAD_BITS );
		p58 = p58 + ( (p58 >> DAT_BITS) << PAD_BITS );
		p60 = p60 + ( (p60 >> DAT_BITS) << PAD_BITS );
		p68 = p68 + ( (p68 >> DAT_BITS) << PAD_BITS );
		p70 = p70 + ( (p70 >> DAT_BITS) << PAD_BITS );
		p78 = p78 + ( (p78 >> DAT_BITS) << PAD_BITS );

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
		i_offsets[0x08] = p08;
		i_offsets[0x09] = p09;
		i_offsets[0x0a] = p0a;
		i_offsets[0x0b] = p0b;
		i_offsets[0x0c] = p0c;
		i_offsets[0x0d] = p0d;
		i_offsets[0x0e] = p0e;
		i_offsets[0x0f] = p0f;
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x10] = i_offsets[j] + p10; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x20] = i_offsets[j] + p20; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x30] = i_offsets[j] + p30; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x40] = i_offsets[j] + p40; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x50] = i_offsets[j] + p50; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x60] = i_offsets[j] + p60; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x70] = i_offsets[j] + p70; }
	#endif
	}

/*...The radix-128 pass is here.	*/

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

	#if USE_SCALAR_DFT_MACRO

		RADIX_128_DIF(
			(a+j1),i_offsets,RE_IM_STRIDE,
			(a+j1),i_offsets,RE_IM_STRIDE
		);

	#else

	// Gather the needed data and do 8 twiddleless length-16 subtransforms, with p-offsets in br8 order: 04261537:
	// NOTE that unlike the RADIX_08_DIF_OOP() macro used for pass 1 of the radix-64 DFT, RADIX_16_DIF outputs are IN-ORDER rather than BR:
		/*...Block 0: */		jt = j1;	jp = j2;
		RADIX_16_DIF(
			a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],a[jt+p40],a[jp+p40],a[jt+p48],a[jp+p48],a[jt+p50],a[jp+p50],a[jt+p58],a[jp+p58],a[jt+p60],a[jp+p60],a[jt+p68],a[jp+p68],a[jt+p70],a[jp+p70],a[jt+p78],a[jp+p78],
			t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,t0fr,t0fi,
			c16,s16	// These = (c16,s16) typically def'd for use in the radix-16 DFT
		);
		/*...Block 1: */		jt = j1 + p04;	jp = j2 + p04;
		RADIX_16_DIF(
			a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],a[jt+p40],a[jp+p40],a[jt+p48],a[jp+p48],a[jt+p50],a[jp+p50],a[jt+p58],a[jp+p58],a[jt+p60],a[jp+p60],a[jt+p68],a[jp+p68],a[jt+p70],a[jp+p70],a[jt+p78],a[jp+p78],
			t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,t1fr,t1fi,
			c16,s16
		);
		/*...Block 2: */		jt = j1 + p02;	jp = j2 + p02;
		RADIX_16_DIF(
			a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],a[jt+p40],a[jp+p40],a[jt+p48],a[jp+p48],a[jt+p50],a[jp+p50],a[jt+p58],a[jp+p58],a[jt+p60],a[jp+p60],a[jt+p68],a[jp+p68],a[jt+p70],a[jp+p70],a[jt+p78],a[jp+p78],
			t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,t2fr,t2fi,
			c16,s16
		);
		/*...Block 3: */		jt = j1 + p06;	jp = j2 + p06;
		RADIX_16_DIF(
			a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],a[jt+p40],a[jp+p40],a[jt+p48],a[jp+p48],a[jt+p50],a[jp+p50],a[jt+p58],a[jp+p58],a[jt+p60],a[jp+p60],a[jt+p68],a[jp+p68],a[jt+p70],a[jp+p70],a[jt+p78],a[jp+p78],
			t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei,t3fr,t3fi,
			c16,s16
		);
		/*...Block 4: */		jt = j1 + p01;	jp = j2 + p01;
		RADIX_16_DIF(
			a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],a[jt+p40],a[jp+p40],a[jt+p48],a[jp+p48],a[jt+p50],a[jp+p50],a[jt+p58],a[jp+p58],a[jt+p60],a[jp+p60],a[jt+p68],a[jp+p68],a[jt+p70],a[jp+p70],a[jt+p78],a[jp+p78],
			t40r,t40i,t41r,t41i,t42r,t42i,t43r,t43i,t44r,t44i,t45r,t45i,t46r,t46i,t47r,t47i,t48r,t48i,t49r,t49i,t4ar,t4ai,t4br,t4bi,t4cr,t4ci,t4dr,t4di,t4er,t4ei,t4fr,t4fi,
			c16,s16
		);
		/*...Block 5: */		jt = j1 + p05;	jp = j2 + p05;
		RADIX_16_DIF(
			a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],a[jt+p40],a[jp+p40],a[jt+p48],a[jp+p48],a[jt+p50],a[jp+p50],a[jt+p58],a[jp+p58],a[jt+p60],a[jp+p60],a[jt+p68],a[jp+p68],a[jt+p70],a[jp+p70],a[jt+p78],a[jp+p78],
			t50r,t50i,t51r,t51i,t52r,t52i,t53r,t53i,t54r,t54i,t55r,t55i,t56r,t56i,t57r,t57i,t58r,t58i,t59r,t59i,t5ar,t5ai,t5br,t5bi,t5cr,t5ci,t5dr,t5di,t5er,t5ei,t5fr,t5fi,
			c16,s16
		);
		/*...Block 6: */		jt = j1 + p03;	jp = j2 + p03;
		RADIX_16_DIF(
			a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],a[jt+p40],a[jp+p40],a[jt+p48],a[jp+p48],a[jt+p50],a[jp+p50],a[jt+p58],a[jp+p58],a[jt+p60],a[jp+p60],a[jt+p68],a[jp+p68],a[jt+p70],a[jp+p70],a[jt+p78],a[jp+p78],
			t60r,t60i,t61r,t61i,t62r,t62i,t63r,t63i,t64r,t64i,t65r,t65i,t66r,t66i,t67r,t67i,t68r,t68i,t69r,t69i,t6ar,t6ai,t6br,t6bi,t6cr,t6ci,t6dr,t6di,t6er,t6ei,t6fr,t6fi,
			c16,s16
		);
		/*...Block 7: */		jt = j1 + p07;	jp = j2 + p07;
		RADIX_16_DIF(
			a[jt    ],a[jp    ],a[jt+p08],a[jp+p08],a[jt+p10],a[jp+p10],a[jt+p18],a[jp+p18],a[jt+p20],a[jp+p20],a[jt+p28],a[jp+p28],a[jt+p30],a[jp+p30],a[jt+p38],a[jp+p38],a[jt+p40],a[jp+p40],a[jt+p48],a[jp+p48],a[jt+p50],a[jp+p50],a[jt+p58],a[jp+p58],a[jt+p60],a[jp+p60],a[jt+p68],a[jp+p68],a[jt+p70],a[jp+p70],a[jt+p78],a[jp+p78],
			t70r,t70i,t71r,t71i,t72r,t72i,t73r,t73i,t74r,t74i,t75r,t75i,t76r,t76i,t77r,t77i,t78r,t78i,t79r,t79i,t7ar,t7ai,t7br,t7bi,t7cr,t7ci,t7dr,t7di,t7er,t7ei,t7fr,t7fi,
			c16,s16
		);

	/*...and now do 16 radix-8 subtransforms, including the internal twiddle factors:

		Block 0: twiddles = {   1,     1,    1,   1,   1,   1,   1,   1}
		Block 1: twiddles = {   1,  E^ 1, E^ 2,E^ 3,E^ 4,E^ 5,E^ 6,E^ 7}
		Block 2: twiddles = {   1,  E^ 2, E^ 4,E^ 6,E^ 8,E^ a,E^ c,E^ e}
		Block 3: twiddles = {   1,  E^ 3, E^ 6,E^ 9,E^ c,E^ f,E^12,E^15}
		Block 4: twiddles = {   1,  E^ 4, E^ 8,E^ c,E^10,E^14,E^18,E^1c}
		Block 5: twiddles = {   1,  E^ 5, E^ a,E^ f,E^14,E^19,E^1e,E^23}
		Block 6: twiddles = {   1,  E^ 6, E^ c,E^12,E^18,E^1e,E^24,E^2a}
		Block 7: twiddles = {   1,  E^ 7, E^ e,E^15,E^1c,E^23,E^2a,E^31}
		Block 8: twiddles = {   1,  E^ 8, E^10,E^18,E^20,E^28,E^30,E^38}
		Block 9: twiddles = {   1,  E^ 9, E^12,E^1b,E^24,E^2d,E^36,E^3f}
		Block a: twiddles = {   1,  E^ a, E^14,E^1e,E^28,E^32,E^3c,E^46}
		Block b: twiddles = {   1,  E^ b, E^16,E^21,E^2c,E^37,E^42,E^4d}
		Block c: twiddles = {   1,  E^ c, E^18,E^24,E^30,E^3c,E^48,E^54}
		Block d: twiddles = {   1,  E^ d, E^1a,E^27,E^34,E^41,E^4e,E^5b}
		Block e: twiddles = {   1,  E^ e, E^1c,E^2a,E^38,E^46,E^54,E^62}
		Block f: twiddles = {   1,  E^ f, E^1e,E^2d,E^3c,E^4b,E^5a,E^69}, highest power: 0x69 = 105 = 7*15, checks.

	where ~ denotes complex conjugation, * denotes interchange of real and imaginary part, - denotes negation
	(i.e. ~E^n := ( Re(E),-Im(E)), *E^n := ( Im(E), Re(E)), -E^n := (-Re(E),-Im(E)), and any sequence of these
	operators is evaluated right-to-left, e.g. -~*E^n = -~(*E^n) = -~( Im(E), Re(E)) = -( Im(E),-Re(E)) = (-Im(E),+Re(E));
	note that the - and ~ operators commute, as do the - and *, but ~ and * anticommute, i.e. ~*E^n = -*~E^n) ,
	and the values of the individual exponentials are, in terms of the sincos parameters defined in this module:

		E^ 1 = exp(i* 1*twopi/128) =       ( c128_1, s128_1)
		E^ 2 = exp(i* 2*twopi/128) =       ( c128_2, s128_2) =       ( c64_1, s64_1)
		E^ 3 = exp(i* 3*twopi/128) =       ( c128_3, s128_3)
		E^ 4 = exp(i* 4*twopi/128) =       ( c128_4, s128_4) =       ( c32_1, s32_1)
		E^ 5 = exp(i* 5*twopi/128) =       ( c128_5, s128_5)
		E^ 6 = exp(i* 6*twopi/128) =       ( c128_6, s128_6) =       ( c64_3, s64_3)
		E^ 7 = exp(i* 7*twopi/128) =       ( c128_7, s128_7)
		E^ 8 = exp(i* 8*twopi/128) =       ( c128_8, s128_8) =       ( c16, s16)
		E^ 9 = exp(i* 9*twopi/128) =       ( c128_9, s128_9)
		E^ a = exp(i* a*twopi/128) =       ( c128_a, s128_a) =       ( c64_5, s64_5)
		E^ b = exp(i* b*twopi/128) =       ( c128_b, s128_b)
		E^ c = exp(i* c*twopi/128) =       ( c128_c, s128_c) =       ( c32_3, s32_3)
		E^ d = exp(i* d*twopi/128) =       ( c128_d, s128_d)
		E^ e = exp(i* e*twopi/128) =       ( c128_e, s128_e) =       ( c64_7, s64_7)
		E^ f = exp(i* f*twopi/128) =       ( c128_f, s128_f)
		E^10 = exp(i*10*twopi/128) = isrt2*( 1    , 1    )
		E^20 = exp(i*20*twopi/128) =       ( 0    , 1    ) = I
		E^30 = exp(i*30*twopi/128) = isrt2*(-1    , 1    ) = *~E^10

		E^{11,12,13,14,15,16,17,18,19,1a,1b,1c,1d,1e,1f}
	=  *E^{ f, e, d, c, b, a, 9, 8, 7, 6, 5, 4, 3, 2, 1}

		E^{21,22,23,24,25,26,27,28,29,2a,2b,2c,2d,2e,2f} = I.E^{ 1, 2, 3, 4, 5, 6, 7, 8, 9, a, b, c, d, e, f}
	= *~E^{ 1, 2, 3, 4, 5, 6, 7, 8, 9, a, b, c, d, e, f}	(I.{} denotes complex multiply by I)

		E^{31,32,33,34,35,36,37,38,39,3a,3b,3c,3d,3e,3f}
	= -~E^{ f, e, d, c, b, a, 9, 8, 7, 6, 5, 4, 3, 2, 1}

	and using that E^n = -E^(n-0x40), we get e.g.

		E^69 = -E^29 = -*~E^9 = ~*E^9 = ~[s(9),c(9)] = [+s(9),-c(9)] .

	Thus, exploiting these symmetries allows our 8x8 twiddles matrix to be expressed in terms of the powers E^1-f, 1,
	the imaginary constant I and isrt2 (via E^10 = isrt2*(1+I)) as follows - Note the even rows (counting the top row
	as "row 0") are identical to the twiddles-sets of the radix-64 DFT.

		Block 0: twiddles = {1,     1,     1,     1,     1,     1,     1,     1}
		Block 1: twiddles = {1,  E^ 1,  E^ 2,  E^ 3,  E^ 4,  E^ 5,  E^ 6,  E^ 7}
		Block 2: twiddles = {1,  E^ 2,  E^ 4,  E^ 6,  E^ 8,  E^ a,  E^ c,  E^ e}
		Block 3: twiddles = {1,  E^ 3,  E^ 6,  E^ 9,  E^ c,  E^ f, *E^ e, *E^ b}
		Block 4: twiddles = {1,  E^ 4,  E^ 8,  E^ c,  E^10, *E^ c, *E^ 8, *E^ 4}
		Block 5: twiddles = {1,  E^ 5,  E^ a,  E^ f, *E^ c, *E^ 7, *E^ 2,*~E^ 3}
		Block 6: twiddles = {1,  E^ 6,  E^ c, *E^ e, *E^ 8, *E^ 2,*~E^ 4,*~E^ a}
		Block 7: twiddles = {1,  E^ 7,  E^ e, *E^ b, *E^ 4,*~E^ 3,*~E^ a,-~E^ f}
		Block 8: twiddles = {1,  E^ 8,  E^10, *E^ 8,  I.{},*~E^ 8,*~E^10,-~E^ 8}
		Block 9: twiddles = {1,  E^ 9, *E^ e, *E^ 5,*~E^ 4,*~E^ d,-~E^ a,-~E^ 1}
		Block a: twiddles = {1,  E^ a, *E^ c, *E^ 2,*~E^ 8,-~E^ e,-~E^ 4, -E^ 6}
		Block b: twiddles = {1,  E^ b, *E^ a,*~E^ 1,*~E^ c,-~E^ 9, -E^ 2, -E^ d}
		Block c: twiddles = {1,  E^ c, *E^ 8,*~E^ 4,*~E^10,-~E^ 4, -E^ 8,-*E^ c}
		Block d: twiddles = {1,  E^ d, *E^ 6,*~E^ 7,-~E^ c, -E^ 1, -E^ e,-*E^ 5}
		Block e: twiddles = {1,  E^ e, *E^ 4,*~E^ a,-~E^ 8, -E^ 6,-*E^ c,~*E^ 2}
		Block f: twiddles = {1,  E^ f, *E^ 2,*~E^ d,-~E^ 4, -E^ b,-*E^ 6,~*E^ 9}

	Or with the columns arranged in BR {04261537} order:

		Block 0: twiddles = {1,     1,     1,     1,     1,     1,     1,     1}
		Block 1: twiddles = {1,  E^ 4,  E^ 2,  E^ 6,  E^ 1,  E^ 5,  E^ 3,  E^ 7}
		Block 2: twiddles = {1,  E^ 8,  E^ 4,  E^ c,  E^ 2,  E^ a,  E^ 6,  E^ e}
		Block 3: twiddles = {1,  E^ c,  E^ 6, *E^ e,  E^ 3,  E^ f,  E^ 9, *E^ b}
		Block 4: twiddles = {1,  E^10,  E^ 8, *E^ 8,  E^ 4, *E^ c,  E^ c, *E^ 4}
		Block 5: twiddles = {1, *E^ c,  E^ a, *E^ 2,  E^ 5, *E^ 7,  E^ f,*~E^ 3}
		Block 6: twiddles = {1, *E^ 8,  E^ c,*~E^ 4,  E^ 6, *E^ 2, *E^ e,*~E^ a}
		Block 7: twiddles = {1, *E^ 4,  E^ e,*~E^ a,  E^ 7,*~E^ 3, *E^ b,-~E^ f}
		Block 8: twiddles = {1,  I.{},  E^10,*~E^10,  E^ 8,*~E^ 8, *E^ 8,-~E^ 8}
		Block 9: twiddles = {1,*~E^ 4, *E^ e,-~E^ a,  E^ 9,*~E^ d, *E^ 5,-~E^ 1}
		Block a: twiddles = {1,*~E^ 8, *E^ c,-~E^ 4,  E^ a,-~E^ e, *E^ 2, -E^ 6}
		Block b: twiddles = {1,*~E^ c, *E^ a, -E^ 2,  E^ b,-~E^ 9,*~E^ 1, -E^ d}
		Block c: twiddles = {1,*~E^10, *E^ 8, -E^ 8,  E^ c,-~E^ 4,*~E^ 4,-*E^ c}
		Block d: twiddles = {1,-~E^ c, *E^ 6, -E^ e,  E^ d, -E^ 1,*~E^ 7,-*E^ 5}
		Block e: twiddles = {1,-~E^ 8, *E^ 4,-*E^ c,  E^ e, -E^ 6,*~E^ a,~*E^ 2}
		Block f: twiddles = {1,-~E^ 4, *E^ 2,-*E^ 6,  E^ f, -E^ b,*~E^ d,~*E^ 9}

	Only the last 7 inputs to each of the radix-8 transforms 2 through 8 are multiplied by non-unity twiddles.
	For DIF we process both the blocks, and the twiddles within each block, in bit-reversed order.
	One can see from the data below that aside from using a twiddleless DIF for Block 0 there is
	little to be gained from trying to exploit other "special form" twiddles such as I and isrt2*[+-1,+-1],
	as only 5 of the remaining 15*8 = 120 non-unity twiddles have such a special form (one = I, four of the isrt2-form).
	Thus our radix-8-DIF-with-twiddles macro uses generic complex MUL for the 7 non-unity twiddles of each invocation.
	*/
/**** FOR THE FIRST 8 BLOCKS - EVEN-ORDER POWERS - EXPRESS THINGS IN TERMS OF D :+ E^2, THE RADIX-64 FUNDAMENTAL ROOT ****/
		// Block 0: t*0ri, has all-unity twiddles
		jt = j1;	jp = j2;
		// Remember, the twiddleless DIT also bit-reverses its outputs, so a_p* terms appear in BR-order here [swap index pairs 1/4 and 3/6]:
		RADIX_08_DIF_OOP(
			t00r,t00i,t40r,t40i,t20r,t20i,t60r,t60i,t10r,t10i,t50r,t50i,t30r,t30i,t70r,t70i,
			a[jt],a[jp],a[jt+p04],a[jp+p04],a[jt+p02],a[jp+p02],a[jt+p06],a[jp+p06],a[jt+p01],a[jp+p01],a[jt+p05],a[jp+p05],a[jt+p03],a[jp+p03],a[jt+p07],a[jp+p07]
		);
		// Block 8: t*1ri, BR twiddle row = {   1,  I.{},  D^ 8,*~D^ 8,  D^ 4,*~D^ 4, *D^ 4,-~D^ 4}
		jt = j1 + p08;	jp = j2 + p08;
		RADIX_08_DIF_TWIDDLE_OOP(
			t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,t41r,t41i,t51r,t51i,t61r,t61i,t71r,t71i,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			0,1,ISRT2,ISRT2,-ISRT2,ISRT2,c16,s16,-s16,c16,s16,c16,-c16,s16
		);
		// Block 4: t*2ri, BR twiddle row = {   1,  D^ 8,  D^ 4, *D^ 4,  D^ 2, *D^ 6,  D^ 6, *D^ 2}
		jt = j1 + p10;	jp = j2 + p10;
		RADIX_08_DIF_TWIDDLE_OOP(
			t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,t42r,t42i,t52r,t52i,t62r,t62i,t72r,t72i,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			ISRT2,ISRT2,c16,s16,s16,c16,c32_1,s32_1,s32_3,c32_3,c32_3,s32_3,s32_1,c32_1
		);
		// Block c: t*3ri, BR twiddle row = {   1,*~D^ 8, *D^ 4, -D^ 4,  D^ 6,-~D^ 2,*~D^ 2,-*D^ 6}
		jt = j1 + p18;	jp = j2 + p18;
		RADIX_08_DIF_TWIDDLE_OOP(
			t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,t43r,t43i,t53r,t53i,t63r,t63i,t73r,t73i,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			-ISRT2,ISRT2,s16,c16,-c16,-s16,c32_3,s32_3,-c32_1,s32_1,-s32_1,c32_1,-s32_3,-c32_3
		);
		// Block 2: t*4ri, BR twiddle row = {   1,  D^ 4,  D^ 2,  D^ 6,  D^ 1,  D^ 5,  D^ 3,  D^ 7}
		jt = j1 + p20;	jp = j2 + p20;
		RADIX_08_DIF_TWIDDLE_OOP(
			t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,t44r,t44i,t54r,t54i,t64r,t64i,t74r,t74i,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			c16,s16,c32_1,s32_1,c32_3,s32_3,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7
		);
		// Block a: t*5ri, BR twiddle row = {   1,*~D^ 4, *D^ 6,-~D^ 2,  D^ 5,-~D^ 7, *D^ 1, -D^ 3}
		jt = j1 + p28;	jp = j2 + p28;
		RADIX_08_DIF_TWIDDLE_OOP(
			t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,t45r,t45i,t55r,t55i,t65r,t65i,t75r,t75i,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			-s16,c16,s32_3,c32_3,-c32_1,s32_1,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3
		);
		// Block 6: t*6ri, BR twiddle row = {   1, *D^ 4,  D^ 6,*~D^ 2,  D^ 3, *D^ 1, *D^ 7,*~D^ 5}
		jt = j1 + p30;	jp = j2 + p30;
		RADIX_08_DIF_TWIDDLE_OOP(
			t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,t46r,t46i,t56r,t56i,t66r,t66i,t76r,t76i,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			s16,c16,c32_3,s32_3,-s32_1,c32_1,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5
		);
		// Block e: t*7ri, BR twiddle row = {   1,-~D^ 4, *D^ 2,-*D^ 6,  D^ 7, -D^ 3,*~D^ 5,~*D^ 1} ,
		jt = j1 + p38;	jp = j2 + p38;
		RADIX_08_DIF_TWIDDLE_OOP(
			t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,t47r,t47i,t57r,t57i,t67r,t67i,t77r,t77i,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			-c16,s16,s32_1,c32_1,-s32_3,-c32_3,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1
		);

	/***************************** ODD-ORDER TWIDDLES ROWS: *****************************/

		// Block 1: BR twiddle row = {1,  E^ 4,  E^ 2,  E^ 6,  E^ 1,  E^ 5,  E^ 3,  E^ 7}
		jt = j1 + p40;	jp = j2 + p40;
		RADIX_08_DIF_TWIDDLE_OOP(
			t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,t48r,t48i,t58r,t58i,t68r,t68i,t78r,t78i,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			c32_1,s32_1, c64_1,s64_1, c64_3,s64_3, c128_1,s128_1, c128_5,s128_5, c128_3,s128_3, c128_7,s128_7
		);
		// Block 9: BR twiddle row = {1,*~E^ 4, *E^ e,-~E^ a,  E^ 9,*~E^ d, *E^ 5,-~E^ 1}
		jt = j1 + p48;	jp = j2 + p48;
		RADIX_08_DIF_TWIDDLE_OOP(
			t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,t49r,t49i,t59r,t59i,t69r,t69i,t79r,t79i,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			-s32_1,c32_1, s64_7,c64_7, -c64_5,s64_5, c128_9,s128_9, -s128_d,c128_d, s128_5,c128_5, -c128_1,s128_1
		);
		// Block 5: BR twiddle row = {1, *E^ c,  E^ a, *E^ 2,  E^ 5, *E^ 7,  E^ f,*~E^ 3}
		jt = j1 + p50;	jp = j2 + p50;
		RADIX_08_DIF_TWIDDLE_OOP(
			t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,t4ar,t4ai,t5ar,t5ai,t6ar,t6ai,t7ar,t7ai,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			s32_3,c32_3, c64_5,s64_5, s64_1,c64_1, c128_5,s128_5, s128_7,c128_7, c128_f,s128_f, -s128_3,c128_3
		);
		// Block d: BR twiddle row = {1,-~E^ c, *E^ 6, -E^ e,  E^ d, -E^ 1,*~E^ 7,-*E^ 5}
		jt = j1 + p58;	jp = j2 + p58;
		RADIX_08_DIF_TWIDDLE_OOP(
			t0br,t0bi,t1br,t1bi,t2br,t2bi,t3br,t3bi,t4br,t4bi,t5br,t5bi,t6br,t6bi,t7br,t7bi,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			-c32_3,s32_3, s64_3,c64_3, -c64_7,-s64_7, c128_d,s128_d, -c128_1,-s128_1, -s128_7,c128_7, -s128_5,-c128_5
		);
		// Block 3: BR twiddle row = {1,  E^ c,  E^ 6, *E^ e,  E^ 3,  E^ f,  E^ 9, *E^ b}
		jt = j1 + p60;	jp = j2 + p60;
		RADIX_08_DIF_TWIDDLE_OOP(
			t0cr,t0ci,t1cr,t1ci,t2cr,t2ci,t3cr,t3ci,t4cr,t4ci,t5cr,t5ci,t6cr,t6ci,t7cr,t7ci,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			c32_3,s32_3, c64_3,s64_3, s64_7,c64_7, c128_3,s128_3, c128_f,s128_f, c128_9,s128_9, s128_b,c128_b
		);
		// Block b: BR twiddle row = {1,*~E^ c, *E^ a, -E^ 2,  E^ b,-~E^ 9,*~E^ 1, -E^ d}
		jt = j1 + p68;	jp = j2 + p68;
		RADIX_08_DIF_TWIDDLE_OOP(
			t0dr,t0di,t1dr,t1di,t2dr,t2di,t3dr,t3di,t4dr,t4di,t5dr,t5di,t6dr,t6di,t7dr,t7di,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			-s32_3,c32_3, s64_5,c64_5, -c64_1,-s64_1, c128_b,s128_b, -c128_9,s128_9, -s128_1,c128_1, -c128_d,-s128_d
		);
		// Block 7: BR twiddle row = {1, *E^ 4,  E^ e,*~E^ a,  E^ 7,*~E^ 3, *E^ b,-~E^ f}
		jt = j1 + p70;	jp = j2 + p70;
		RADIX_08_DIF_TWIDDLE_OOP(
			t0er,t0ei,t1er,t1ei,t2er,t2ei,t3er,t3ei,t4er,t4ei,t5er,t5ei,t6er,t6ei,t7er,t7ei,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			s32_1,c32_1, c64_7,s64_7, -s64_5,c64_5, c128_7,s128_7, -s128_3,c128_3, s128_b,c128_b, -c128_f,s128_f
		);
		// Block f: BR twiddle row = {1,-~E^ 4, *E^ 2,-*E^ 6,  E^ f, -E^ b,*~E^ d,~*E^ 9}
		jt = j1 + p78;	jp = j2 + p78;
		RADIX_08_DIF_TWIDDLE_OOP(
			t0fr,t0fi,t1fr,t1fi,t2fr,t2fi,t3fr,t3fi,t4fr,t4fi,t5fr,t5fi,t6fr,t6fi,t7fr,t7fi,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			-c32_1,s32_1, s64_1,c64_1, -s64_3,-c64_3, c128_f,s128_f, -c128_b,-s128_b, -s128_d,c128_d, s128_9,-c128_9
		);
	#endif	// USE_SCALAR_DFT_MACRO ?
	}
}

/**************/

void radix128_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-128 complex inverse DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix128_dif_pass for further details on storage, indexing and DFT decomposition.
*/
	int j,j1,j2,jt,jp;
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p09,p0a,p0b,p0c,p0d,p0e,p0f,p10,p18,p20,p28,p30,p38,p40,p48,p50,p58,p60,p68,p70,p78, first_entry=TRUE;
#if USE_SCALAR_DFT_MACRO
	static int i_offsets[RADIX];
#else
	// We name the even-order roots here the same as we do in the radix-64 code - that allows us to take the
	// radix-8 DFT calls using all-even-order-radix-128 roots more or less straight from the radix-64 code.
	double
		t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,t0fr,t0fi,
		t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,t1fr,t1fi,
		t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,t2fr,t2fi,
		t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei,t3fr,t3fi,
		t40r,t40i,t41r,t41i,t42r,t42i,t43r,t43i,t44r,t44i,t45r,t45i,t46r,t46i,t47r,t47i,t48r,t48i,t49r,t49i,t4ar,t4ai,t4br,t4bi,t4cr,t4ci,t4dr,t4di,t4er,t4ei,t4fr,t4fi,
		t50r,t50i,t51r,t51i,t52r,t52i,t53r,t53i,t54r,t54i,t55r,t55i,t56r,t56i,t57r,t57i,t58r,t58i,t59r,t59i,t5ar,t5ai,t5br,t5bi,t5cr,t5ci,t5dr,t5di,t5er,t5ei,t5fr,t5fi,
		t60r,t60i,t61r,t61i,t62r,t62i,t63r,t63i,t64r,t64i,t65r,t65i,t66r,t66i,t67r,t67i,t68r,t68i,t69r,t69i,t6ar,t6ai,t6br,t6bi,t6cr,t6ci,t6dr,t6di,t6er,t6ei,t6fr,t6fi,
		t70r,t70i,t71r,t71i,t72r,t72i,t73r,t73i,t74r,t74i,t75r,t75i,t76r,t76i,t77r,t77i,t78r,t78i,t79r,t79i,t7ar,t7ai,t7br,t7bi,t7cr,t7ci,t7dr,t7di,t7er,t7ei,t7fr,t7fi;
#endif

	if(!first_entry && (n >> 7) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR = n >> 7;

		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p0a = p09 + p01;
		p0b = p0a + p01;
		p0c = p0b + p01;
		p0d = p0c + p01;
		p0e = p0d + p01;
		p0f = p0e + p01;
		p10 = p08 + p08;
		p18 = p10 + p08;
		p20 = p18 + p08;
		p28 = p20 + p08;
		p30 = p28 + p08;
		p38 = p30 + p08;
		p40 = p38 + p08;
		p48 = p40 + p08;
		p50 = p48 + p08;
		p58 = p50 + p08;
		p60 = p58 + p08;
		p68 = p60 + p08;
		p70 = p68 + p08;
		p78 = p70 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p0a = p0a + ( (p0a >> DAT_BITS) << PAD_BITS );
		p0b = p0b + ( (p0b >> DAT_BITS) << PAD_BITS );
		p0c = p0c + ( (p0c >> DAT_BITS) << PAD_BITS );
		p0d = p0d + ( (p0d >> DAT_BITS) << PAD_BITS );
		p0e = p0e + ( (p0e >> DAT_BITS) << PAD_BITS );
		p0f = p0f + ( (p0f >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p30 = p30 + ( (p30 >> DAT_BITS) << PAD_BITS );
		p38 = p38 + ( (p38 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p50 = p50 + ( (p50 >> DAT_BITS) << PAD_BITS );
		p58 = p58 + ( (p58 >> DAT_BITS) << PAD_BITS );
		p60 = p60 + ( (p60 >> DAT_BITS) << PAD_BITS );
		p68 = p68 + ( (p68 >> DAT_BITS) << PAD_BITS );
		p70 = p70 + ( (p70 >> DAT_BITS) << PAD_BITS );
		p78 = p78 + ( (p78 >> DAT_BITS) << PAD_BITS );

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
		i_offsets[0x08] = p08;
		i_offsets[0x09] = p09;
		i_offsets[0x0a] = p0a;
		i_offsets[0x0b] = p0b;
		i_offsets[0x0c] = p0c;
		i_offsets[0x0d] = p0d;
		i_offsets[0x0e] = p0e;
		i_offsets[0x0f] = p0f;
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x10] = i_offsets[j] + p10; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x20] = i_offsets[j] + p20; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x30] = i_offsets[j] + p30; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x40] = i_offsets[j] + p40; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x50] = i_offsets[j] + p50; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x60] = i_offsets[j] + p60; }
		for(j = 0; j < 0x10; ++j) { i_offsets[j+0x70] = i_offsets[j] + p70; }
	#endif
	}

/*...The radix-128 pass is here.	*/

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

	#if USE_SCALAR_DFT_MACRO

		RADIX_128_DIT(
			(a+j1),i_offsets,RE_IM_STRIDE,
			(a+j1),i_offsets,RE_IM_STRIDE
		);

	#else

	// Gather the needed data and do 8 twiddleless length-16 subtransforms:
		/*...Block 0: */
		jt = j1;	jp = j2;
		RADIX_16_DIT(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
			t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,t0fr,t0fi,
			c16,s16	// These = (c16,s16) typically def'd for use in the radix-16 DFT
		);
		/*...Block 1: */
		jt = j1 + p10;	jp = j2 + p10;
		RADIX_16_DIT(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
			t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,t1fr,t1fi,
			c16,s16
		);
		/*...Block 2: */
		jt = j1 + p20;	jp = j2 + p20;
		RADIX_16_DIT(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
			t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,t2fr,t2fi,
			c16,s16
		);
		/*...Block 3: */
		jt = j1 + p30;	jp = j2 + p30;
		RADIX_16_DIT(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
			t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei,t3fr,t3fi,
			c16,s16
		);
		/*...Block 4: */
		jt = j1 + p40;	jp = j2 + p40;
		RADIX_16_DIT(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
			t40r,t40i,t41r,t41i,t42r,t42i,t43r,t43i,t44r,t44i,t45r,t45i,t46r,t46i,t47r,t47i,t48r,t48i,t49r,t49i,t4ar,t4ai,t4br,t4bi,t4cr,t4ci,t4dr,t4di,t4er,t4ei,t4fr,t4fi,
			c16,s16
		);
		/*...Block 5: */
		jt = j1 + p50;	jp = j2 + p50;
		RADIX_16_DIT(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
			t50r,t50i,t51r,t51i,t52r,t52i,t53r,t53i,t54r,t54i,t55r,t55i,t56r,t56i,t57r,t57i,t58r,t58i,t59r,t59i,t5ar,t5ai,t5br,t5bi,t5cr,t5ci,t5dr,t5di,t5er,t5ei,t5fr,t5fi,
			c16,s16
		);
		/*...Block 6: */
		jt = j1 + p60;	jp = j2 + p60;
		RADIX_16_DIT(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
			t60r,t60i,t61r,t61i,t62r,t62i,t63r,t63i,t64r,t64i,t65r,t65i,t66r,t66i,t67r,t67i,t68r,t68i,t69r,t69i,t6ar,t6ai,t6br,t6bi,t6cr,t6ci,t6dr,t6di,t6er,t6ei,t6fr,t6fi,
			c16,s16
		);
		/*...Block 7: */
		jt = j1 + p70;	jp = j2 + p70;
		RADIX_16_DIT(
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p08],a[jp+p08],a[jt+p09],a[jp+p09],a[jt+p0a],a[jp+p0a],a[jt+p0b],a[jp+p0b],a[jt+p0c],a[jp+p0c],a[jt+p0d],a[jp+p0d],a[jt+p0e],a[jp+p0e],a[jt+p0f],a[jp+p0f],
			t70r,t70i,t71r,t71i,t72r,t72i,t73r,t73i,t74r,t74i,t75r,t75i,t76r,t76i,t77r,t77i,t78r,t78i,t79r,t79i,t7ar,t7ai,t7br,t7bi,t7cr,t7ci,t7dr,t7di,t7er,t7ei,t7fr,t7fi,
			c16,s16
		);

	/*...and now do 16 radix-8 subtransforms, including the internal twiddle factors - we use the same positive-power
	roots as in the DIF here, just fiddle with signs within the macro to effect the conjugate-multiplies. Twiddles occur
	in the same order here as DIF, but the in-and-output-index offsets are BRed: j1 + p[0,8,4,c,2,a,6,e,1,9,5,d,3,b,7,f],
	as are the index offsets of each sets of complex outputs in the A-array: [jt,jp] + p10*[0,4,2,6,1,5,3,7]:
	*/
		/* Block 0: t*0,1 twiddles = {   1,     1,     1,     1,     1,     1,     1,     1}: */
		jt = j1;	jp = j2;
		/* 0-index block has all-unity twiddles: Remember, the twiddleless DIT also bit-reverses its outputs, so a_p* terms appear in-order here: */
		RADIX_08_DIT_OOP(
			t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,t40r,t40i,t50r,t50i,t60r,t60i,t70r,t70i,
			a[jt],a[jp],a[jt+p10],a[jp+p10],a[jt+p20],a[jp+p20],a[jt+p30],a[jp+p30],a[jt+p40],a[jp+p40],a[jt+p50],a[jp+p50],a[jt+p60],a[jp+p60],a[jt+p70],a[jp+p70]
		);
		// Block 8: t*8ri, BR twiddle row = {   1,  I.{},  D^ 8,*~D^ 8,  D^ 4,*~D^ 4, *D^ 4,-~D^ 4}
		jt = j1 + p08;	jp = j2 + p08;
		RADIX_08_DIT_TWIDDLE_OOP(
			t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,t48r,t48i,t58r,t58i,t68r,t68i,t78r,t78i,
			a[jt],a[jp],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],
			0,1,ISRT2,ISRT2,-ISRT2,ISRT2,c16,s16,-s16,c16,s16,c16,-c16,s16
		);
		// Block 4: t*4ri, BR twiddle row = {   1,  D^ 8,  D^ 4, *D^ 4,  D^ 2, *D^ 6,  D^ 6, *D^ 2}
		jt = j1 + p04;	jp = j2 + p04;
		RADIX_08_DIT_TWIDDLE_OOP(
			t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,t44r,t44i,t54r,t54i,t64r,t64i,t74r,t74i,
			a[jt],a[jp],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],
			ISRT2,ISRT2,c16,s16,s16,c16,c32_1,s32_1,s32_3,c32_3,c32_3,s32_3,s32_1,c32_1
		);
		// Block c: t*cri, BR twiddle row = {   1,*~D^ 8, *D^ 4, -D^ 4,  D^ 6,-~D^ 2,*~D^ 2,-*D^ 6}
		jt = j1 + p0c;	jp = j2 + p0c;
		RADIX_08_DIT_TWIDDLE_OOP(
			t0cr,t0ci,t1cr,t1ci,t2cr,t2ci,t3cr,t3ci,t4cr,t4ci,t5cr,t5ci,t6cr,t6ci,t7cr,t7ci,
			a[jt],a[jp],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],
			-ISRT2,ISRT2,s16,c16,-c16,-s16,c32_3,s32_3,-c32_1,s32_1,-s32_1,c32_1,-s32_3,-c32_3
		);
		// Block 2: t*2ri, BR twiddle row = {   1,  D^ 4,  D^ 2,  D^ 6,  D^ 1,  D^ 5,  D^ 3,  D^ 7}
		jt = j1 + p02;	jp = j2 + p02;
		RADIX_08_DIT_TWIDDLE_OOP(
			t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,t42r,t42i,t52r,t52i,t62r,t62i,t72r,t72i,
			a[jt],a[jp],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],
			c16,s16,c32_1,s32_1,c32_3,s32_3,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7
		);
		// Block a: t*ari, BR twiddle row = {   1,*~D^ 4, *D^ 6,-~D^ 2,  D^ 5,-~D^ 7, *D^ 1, -D^ 3}
		jt = j1 + p0a;	jp = j2 + p0a;
		RADIX_08_DIT_TWIDDLE_OOP(
			t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,t4ar,t4ai,t5ar,t5ai,t6ar,t6ai,t7ar,t7ai,
			a[jt],a[jp],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],
			-s16,c16,s32_3,c32_3,-c32_1,s32_1,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3
		);
		// Block 6: t*6ri, BR twiddle row = {   1, *D^ 4,  D^ 6,*~D^ 2,  D^ 3, *D^ 1, *D^ 7,*~D^ 5}
		jt = j1 + p06;	jp = j2 + p06;
		RADIX_08_DIT_TWIDDLE_OOP(
			t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,t46r,t46i,t56r,t56i,t66r,t66i,t76r,t76i,
			a[jt],a[jp],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],
			s16,c16,c32_3,s32_3,-s32_1,c32_1,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5
		);
		// Block e: t*eri, BR twiddle row = {   1,-~D^ 4, *D^ 2,-*D^ 6,  D^ 7, -D^ 3,*~D^ 5,~*D^ 1} ,
		jt = j1 + p0e;	jp = j2 + p0e;
		RADIX_08_DIT_TWIDDLE_OOP(
			t0er,t0ei,t1er,t1ei,t2er,t2ei,t3er,t3ei,t4er,t4ei,t5er,t5ei,t6er,t6ei,t7er,t7ei,
			a[jt],a[jp],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],
			-c16,s16,s32_1,c32_1,-s32_3,-c32_3,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1
		);

	/***************************** ODD-ORDER TWIDDLES ROWS: *****************************/

		// Block 1: BR twiddle row = {1,  E^ 4,  E^ 2,  E^ 6,  E^ 1,  E^ 5,  E^ 3,  E^ 7}
		jt = j1 + p01;	jp = j2 + p01;
		RADIX_08_DIT_TWIDDLE_OOP(
			t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,t41r,t41i,t51r,t51i,t61r,t61i,t71r,t71i,
			a[jt],a[jp],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],
			c32_1,s32_1, c64_1,s64_1, c64_3,s64_3, c128_1,s128_1, c128_5,s128_5, c128_3,s128_3, c128_7,s128_7
		);
		// Block 9: BR twiddle row = {1,*~E^ 4, *E^ e,-~E^ a,  E^ 9,*~E^ d, *E^ 5,-~E^ 1}
		jt = j1 + p09;	jp = j2 + p09;
		RADIX_08_DIT_TWIDDLE_OOP(
			t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,t49r,t49i,t59r,t59i,t69r,t69i,t79r,t79i,
			a[jt],a[jp],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],
			-s32_1,c32_1, s64_7,c64_7, -c64_5,s64_5, c128_9,s128_9, -s128_d,c128_d, s128_5,c128_5, -c128_1,s128_1
		);
		// Block 5: BR twiddle row = {1, *E^ c,  E^ a, *E^ 2,  E^ 5, *E^ 7,  E^ f,*~E^ 3}
		jt = j1 + p05;	jp = j2 + p05;
		RADIX_08_DIT_TWIDDLE_OOP(
			t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,t45r,t45i,t55r,t55i,t65r,t65i,t75r,t75i,
			a[jt],a[jp],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],
			s32_3,c32_3, c64_5,s64_5, s64_1,c64_1, c128_5,s128_5, s128_7,c128_7, c128_f,s128_f, -s128_3,c128_3
		);
		// Block d: BR twiddle row = {1,-~E^ c, *E^ 6, -E^ e,  E^ d, -E^ 1,*~E^ 7,-*E^ 5}
		jt = j1 + p0d;	jp = j2 + p0d;
		RADIX_08_DIT_TWIDDLE_OOP(
			t0dr,t0di,t1dr,t1di,t2dr,t2di,t3dr,t3di,t4dr,t4di,t5dr,t5di,t6dr,t6di,t7dr,t7di,
			a[jt],a[jp],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],
			-c32_3,s32_3, s64_3,c64_3, -c64_7,-s64_7, c128_d,s128_d, -c128_1,-s128_1, -s128_7,c128_7, -s128_5,-c128_5
		);
		// Block 3: BR twiddle row = {1,  E^ c,  E^ 6, *E^ e,  E^ 3,  E^ f,  E^ 9, *E^ b}
		jt = j1 + p03;	jp = j2 + p03;
		RADIX_08_DIT_TWIDDLE_OOP(
			t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,t43r,t43i,t53r,t53i,t63r,t63i,t73r,t73i,
			a[jt],a[jp],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],
			c32_3,s32_3, c64_3,s64_3, s64_7,c64_7, c128_3,s128_3, c128_f,s128_f, c128_9,s128_9, s128_b,c128_b
		);
		// Block b: BR twiddle row = {1,*~E^ c, *E^ a, -E^ 2,  E^ b,-~E^ 9,*~E^ 1, -E^ d}
		jt = j1 + p0b;	jp = j2 + p0b;
		RADIX_08_DIT_TWIDDLE_OOP(
			t0br,t0bi,t1br,t1bi,t2br,t2bi,t3br,t3bi,t4br,t4bi,t5br,t5bi,t6br,t6bi,t7br,t7bi,
			a[jt],a[jp],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],
			-s32_3,c32_3, s64_5,c64_5, -c64_1,-s64_1, c128_b,s128_b, -c128_9,s128_9, -s128_1,c128_1, -c128_d,-s128_d
		);
		// Block 7: BR twiddle row = {1, *E^ 4,  E^ e,*~E^ a,  E^ 7,*~E^ 3, *E^ b,-~E^ f}
		jt = j1 + p07;	jp = j2 + p07;
		RADIX_08_DIT_TWIDDLE_OOP(
			t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,t47r,t47i,t57r,t57i,t67r,t67i,t77r,t77i,
			a[jt],a[jp],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],
			s32_1,c32_1, c64_7,s64_7, -s64_5,c64_5, c128_7,s128_7, -s128_3,c128_3, s128_b,c128_b, -c128_f,s128_f
		);
		// Block f: BR twiddle row = {1,-~E^ 4, *E^ 2,-*E^ 6,  E^ f, -E^ b,*~E^ d,~*E^ 9}
		jt = j1 + p0f;	jp = j2 + p0f;
		RADIX_08_DIT_TWIDDLE_OOP(
			t0fr,t0fi,t1fr,t1fi,t2fr,t2fi,t3fr,t3fi,t4fr,t4fi,t5fr,t5fi,t6fr,t6fi,t7fr,t7fi,
			a[jt],a[jp],a[jt+p40],a[jp+p40],a[jt+p20],a[jp+p20],a[jt+p60],a[jp+p60],a[jt+p10],a[jp+p10],a[jt+p50],a[jp+p50],a[jt+p30],a[jp+p30],a[jt+p70],a[jp+p70],
			-c32_1,s32_1, s64_1,c64_1, -s64_3,-c64_3, c128_f,s128_f, -c128_b,-s128_b, -s128_d,c128_d, s128_9,-c128_9
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
	cy128_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
	#ifdef USE_SSE2
		double *addr,*addi;
	#else
		struct complex *tptr;
	#endif
	#ifdef USE_SSE2
		const int pfetch_dist = PFETCH_DIST;
	#endif
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p01,p02,p03,p04,p05,p06,p07,p08,p09,p0a,p0b,p0c,p0d,p0e,p0f,p10,p18,p20,p28,p30,p38,p40,p48,p50,p58,p60,p68,p70,p78;
		int poff[RADIX>>2];	// Store mults of p-offsets for loop-controlled DFT macro calls
	#if 0
		int po_lin[16];		// Store mults of p-offsets for loop-controlled DFT macro calls
	#endif
	#ifndef USE_SSE2
		int po_br[16];		// Store mults of p-offsets for loop-controlled DFT macro calls
	#endif
		int j,j1,l;
	#ifndef USE_SSE2
		int j2,jt,jp;
	#endif
	#ifdef USE_SSE2
		int incr;
		// incr = Carry-chain wts-multipliers recurrence length, which must divide
		// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 64|128|256 for avx512,avx,sse, respectively:
		const int incr_long = 8,incr_short = 4;
		if(USE_SHORT_CY_CHAIN)	// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
			incr = incr_short;
		else
			incr = incr_long;
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
	#if !defined(USE_SSE2) || defined(USE_AVX)
		int k1,k2;
	#endif

	#ifdef USE_SSE2

	  #ifndef USE_AVX512
		const double crnd = 3.0*0x4000000*0x2000000;
	  #endif
		double *add0,*add1,*add2,*add3/* ,*add4,*add5,*add6,*add7 */;	/* Addresses into array sections */
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
			*half_arr, *two,/* *one,*sqrt2, */*isrt2,/* *cc0,*ss0, */
	  #if COMPACT_OBJ
			// ptrs to 16 sets of twiddles shared by the 2nd-half DIF and DIT DFT macros:
			*twid0,
	  #else
		#error COMPACT_OBJ must be defined nonzero!
	  #endif
			*r00,//*r10,*r20,*r30,*r40,*r50,*r60,*r70,
			*s1p00,//*s1p10,*s1p20,*s1p30,*s1p40,*s1p50,*s1p60,*s1p70,
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

		int p0123[4];
		double *base, *baseinv;
		const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int m,m2,ntmp;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		int bjmodn[RADIX];	// Thread only carries a base datum here, must alloc a local array for remaining values
		double *cy_r = thread_arg->cy_r,*cy_i = thread_arg->cy_i, temp,frac;
		// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
		// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
		// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
		struct complex t[RADIX];
		#include "radix128_twiddles.h"
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

		/*   constant index offsets for load/stores are here.	*/
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p0a = p09 + p01;
		p0b = p0a + p01;
		p0c = p0b + p01;
		p0d = p0c + p01;
		p0e = p0d + p01;
		p0f = p0e + p01;
		p10 = p08 + p08;
		p18 = p10 + p08;
		p20 = p18 + p08;
		p28 = p20 + p08;
		p30 = p28 + p08;
		p38 = p30 + p08;
		p40 = p38 + p08;
		p48 = p40 + p08;
		p50 = p48 + p08;
		p58 = p50 + p08;
		p60 = p58 + p08;
		p68 = p60 + p08;
		p70 = p68 + p08;
		p78 = p70 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p0a = p0a + ( (p0a >> DAT_BITS) << PAD_BITS );
		p0b = p0b + ( (p0b >> DAT_BITS) << PAD_BITS );
		p0c = p0c + ( (p0c >> DAT_BITS) << PAD_BITS );
		p0d = p0d + ( (p0d >> DAT_BITS) << PAD_BITS );
		p0e = p0e + ( (p0e >> DAT_BITS) << PAD_BITS );
		p0f = p0f + ( (p0f >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p30 = p30 + ( (p30 >> DAT_BITS) << PAD_BITS );
		p38 = p38 + ( (p38 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p50 = p50 + ( (p50 >> DAT_BITS) << PAD_BITS );
		p58 = p58 + ( (p58 >> DAT_BITS) << PAD_BITS );
		p60 = p60 + ( (p60 >> DAT_BITS) << PAD_BITS );
		p68 = p68 + ( (p68 >> DAT_BITS) << PAD_BITS );
		p70 = p70 + ( (p70 >> DAT_BITS) << PAD_BITS );
		p78 = p78 + ( (p78 >> DAT_BITS) << PAD_BITS );
	#ifndef USE_SSE2
		p0123[0] = 0; p0123[1] = p01; p0123[2] = p02; p0123[3] = p03;
	#endif
		poff[0x00] =   0; poff[0x01] = p04    ; poff[0x02] = p08; poff[0x03] = p04+p08;
		poff[0x04] = p10; poff[0x05] = p04+p10; poff[0x06] = p18; poff[0x07] = p04+p18;
		poff[0x08] = p20; poff[0x09] = p04+p20; poff[0x0a] = p28; poff[0x0b] = p04+p28;
		poff[0x0c] = p30; poff[0x0d] = p04+p30; poff[0x0e] = p38; poff[0x0f] = p04+p38;
		poff[0x10] = p40; poff[0x11] = p04+p40; poff[0x12] = p48; poff[0x13] = p04+p48;
		poff[0x14] = p50; poff[0x15] = p04+p50; poff[0x16] = p58; poff[0x17] = p04+p58;
		poff[0x18] = p60; poff[0x19] = p04+p60; poff[0x1a] = p68; poff[0x1b] = p04+p68;
		poff[0x1c] = p70; poff[0x1d] = p04+p70; poff[0x1e] = p78; poff[0x1f] = p04+p78;
	#if 0
		// _lin = linear p-multiples (but padding means we can't assume e.g. p02 = 2*p01):
		po_lin[0x0] =   0; po_lin[0x1] = p01; po_lin[0x2] = p02; po_lin[0x3] = p03;
		po_lin[0x4] = p04; po_lin[0x5] = p05; po_lin[0x6] = p06; po_lin[0x7] = p07;
		po_lin[0x8] = p08; po_lin[0x9] = p09; po_lin[0xa] = p0a; po_lin[0xb] = p0b;
		po_lin[0xc] = p0c; po_lin[0xd] = p0d; po_lin[0xe] = p0e; po_lin[0xf] = p0f;
	#endif
	#ifndef USE_SSE2
		// _br  = bit-reversed p-multiples:
		po_br[0x0] =   0; po_br[0x1] = p08; po_br[0x2] = p04; po_br[0x3] = p0c;
		po_br[0x4] = p02; po_br[0x5] = p0a; po_br[0x6] = p06; po_br[0x7] = p0e;
		po_br[0x8] = p01; po_br[0x9] = p09; po_br[0xa] = p05; po_br[0xb] = p0d;
		po_br[0xc] = p03; po_br[0xd] = p0b; po_br[0xe] = p07; po_br[0xf] = p0f;
	#endif

	#ifdef USE_SSE2
		tmp = thread_arg->r00;	// declared above
		r00 = tmp + 0x00;	//r40 = tmp + 0x80;
		//r10 = tmp + 0x20;	r50 = tmp + 0xa0;
		//r20 = tmp + 0x40;	r60 = tmp + 0xc0;
		//r30 = tmp + 0x60;	r70 = tmp + 0xe0;
		tmp += 0x100;
		s1p00 = tmp + 0x00;	//s1p40 = tmp + 0x80;
		//s1p10 = tmp + 0x20;	s1p50 = tmp + 0xa0;
		//s1p20 = tmp + 0x40;	s1p60 = tmp + 0xc0;
		//s1p30 = tmp + 0x60;	s1p70 = tmp + 0xe0;
		tmp += 0x100;
		two     = tmp + 0;	// AVX+ versions of radix-8,16,32 twiddleless-DFT macros need consts [2,1,sqrt2,isrt2] quartet laid out thusly
		//one     = tmp + 1;
		//sqrt2   = tmp + 2;
		// Next 3 slots unnamed in non-compact-obj-code mode, since that layout below already has iart2,c16,s16 pointers
		isrt2   = tmp + 3;
		//cc0		= tmp + 4;	// ...and the radix-16 DFT macrod need the base 16-root [re,im] pair right above isrt2.
		//ss0		= tmp + 5;
		tmp += 6;
		// Each non-unity root now needs a negated counterpart:
		/* Stupidity: Since a truly general-purpose [in the sense that it can be used for our radix-128 internal-twiddles]
		radix-8 DFT-with-twiddles macro needs 8 in-addresses [corr. to the 8 real parts of the input data], 8 o-addresses,
		and 7 each of cosine and sine data [which cannot be assumed to occur in fixed-stride pairs - cf. our usage of
		SSE2_RADIX8_DIT_TWIDDLE_OOP() below], that hits the GCC hard limit of 30-operands for ASM macros, but we still
		need one more operand for the ISRT2 pointer. Only easy workaround I found for this is to stick a vector-ISRT2 copy
		in between each +-[cc,ss] vector-data pair, thus any time we need a vector-isrt2 for the radix-8 internal twiddles
		we get it at (vec_dbl*)cc-1.	/Stupidity */
	  #if COMPACT_OBJ
		isrt2   = tmp - 3;	// Alternate layout below already has [isrt2,cc0,ss0] pointer-triplet in a different place, leave these 3 slots unnamed there
		//cc0		= tmp - 2;
		//ss0		= tmp - 1;
		// ptrs to 16 sets (3*7 = 21 vec_dbl data each) of non-unity twiddles shared by the 2nd-half DIF and DIT DFT macros:
		twid0  = tmp + 0x00;
		tmp += 0x150;	// += 16*21 = 0x150 => sc_ptr + 0x356 (alternate version below incurs 0x60 elts, thus diff = 0xf0)
	  #else
		#error COMPACT_OBJ must be defined nonzero!
	  #endif
	  #ifdef USE_AVX512
		cy_r = tmp;	cy_i = tmp+0x10;	tmp += 2*0x10;	// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		//sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod
	  #elif defined(USE_AVX)
		cy_r = tmp;	cy_i = tmp+0x20;	tmp += 2*0x20;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 0x(2a8 + diff) vec_dbl
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 68 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod
	  #else
		cy_r = tmp;	/* cy_i = tmp+0x40; */	tmp += 2*0x40;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 0x(2e8 + diff)
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, 2 for Fermat-mod */
	  #endif

		ASSERT((r00 == thread_arg->r00), "thread-local memcheck failed!");
		ASSERT((half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
		ASSERT((two->d0 == 2.0 && two->d1 == 2.0), "thread-local memcheck failed!");
		// Must make this check 'fuzzy' to allow for wrong-way-round experiments:
		ASSERT((fabs(isrt2->d0 - ISRT2) < EPS && fabs(isrt2->d1 - ISRT2) < EPS), "thread-local memcheck failed!");
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

		sign_mask = (uint64*)(r00 + radix128_creals_in_local_store);
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
		#ifdef USE_SSE2
			addr = thread_arg->cy_r;
		#endif
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
		#ifdef USE_SSE2
			addr = thread_arg->cy_r;	addi = thread_arg->cy_i;
		#endif
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
		#include "radix128_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_SSE2
			addr = thread_arg->cy_r;
		#endif
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
		#ifdef USE_SSE2
			addr = thread_arg->cy_r;	addi = thread_arg->cy_i;
		#endif
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
