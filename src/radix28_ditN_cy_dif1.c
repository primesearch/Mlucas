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

#define RADIX 28	// 0x1c; Use #define rather than const int to ensure it's really a compile-time const in the C sense
#define ODD_RADIX 7	// ODD_RADIX = [radix >> trailz(radix)]

#define EPS 1e-10

#ifndef PFETCH_DIST
  #ifdef USE_AVX
	#define PFETCH_DIST	32	// This seems to work best on my Haswell, even though 64 bytes seems more logical in AVX mode
  #else
	#define PFETCH_DIST	32
  #endif
#endif

#ifdef USE_SSE2

  // For Mersenne-mod need (16 [SSE2] or 64 [AVX]) + (4 [HIACC] or 40 [LOACC]) added slots for half_arr lookup tables.
  // Max = (40 [SSE2]; 96 [AVX]).
  // For Fermat-mod in AVX mode we need RADIX*4 = 0x70 [if HIACC] or 0xc [if not] slots for the compact
  // negacyclic-roots chained-multiply scheme.
  // Add larger number in each case - i.e. max(0x50,0x70) = 0x70 if AVX, 0x28 if SSE2
  // to (half_arr_offset28 + RADIX) = (82-RADIX/2) + 16 = 90 to get required value of radix28_creals_in_local_store:
  #ifdef USE_AVX
	const int half_arr_offset28 = 0x56;	// + RADIX = 0x72; Used for thread local-storage-integrity checking
	const int radix28_creals_in_local_store = 0xf8;	// AVX: 0x72 + 0x84(=132) and round up to nearest multiple of 8
  #else
	const int half_arr_offset28 = 0x64;	// + RADIX = 0x80; Used for thread local-storage-integrity checking
	const int radix28_creals_in_local_store = 0xa8;	// SSE2: 0x80 + 0x28 and round up to nearest multiple of 8
  #endif

	/*
	Here more details about the small-cyclic-array indexing scheme used in the fermat_carry_norm_errcheckB macros.
	Here are sample data for:

	F24 using a length-7*2^17 transform, complex radices	F24 using a length-7*2^18 transform, complex radices
	= 28,16,32,32 [only the leading 28 matters here]:		= 28,32,32,32:

		J = 0:												J = 0:
		ii0 = 0, bjmodn00 = 917504, i = 1					ii0 = 0, bjmodn00 = 1835008, i = 1
		ii1 = 6, bjmodn01 = 131072, i = 0					ii1 = 5, bjmodn01 =  524288, i = 0
		ii2 = 5, bjmodn02 = 262144, i = 0					ii2 = 3, bjmodn02 = 1048576, i = 0
		ii3 = 4, bjmodn03 = 393216, i = 0					ii3 = 1, bjmodn03 = 1572864, i = 1
		ii4 = 3, bjmodn04 = 524288, i = 0					ii4 = 6, bjmodn04 =  262144, i = 0
		ii5 = 2, bjmodn05 = 655360, i = 0					ii5 = 4, bjmodn05 =  786432, i = 0
		ii6 = 1, bjmodn06 = 786432, i = 1					ii6 = 2, bjmodn06 = 1310720, i = 0
		J = 2:												J = 2:
		ii0 = 5, bjmodn00 = 262144, i = 0					ii0 = 5, bjmodn00 =  524288, i = 0
		ii1 = 4, bjmodn01 = 393216, i = 0					ii1 = 3, bjmodn01 = 1048576, i = 0
		ii2 = 3, bjmodn02 = 524288, i = 0					ii2 = 1, bjmodn02 = 1572864, i = 1
		ii3 = 2, bjmodn03 = 655360, i = 0					ii3 = 6, bjmodn03 =  262144, i = 0
		ii4 = 1, bjmodn04 = 786432, i = 1					ii4 = 4, bjmodn04 =  786432, i = 0
		ii5 = 0, bjmodn05 = 917504, i = 1					ii5 = 2, bjmodn05 = 1310720, i = 0
		ii6 = 6, bjmodn06 = 131072, i = 0					ii6 = 0, bjmodn06 = 1835008, i = 1
		J = 4:												J = 4:
		ii0 = 3, bjmodn00 = 524288, i = 0					ii0 = 3, bjmodn00 = 1048576, i = 0
		ii1 = 2, bjmodn01 = 655360, i = 0					ii1 = 1, bjmodn01 = 1572864, i = 1
		ii2 = 1, bjmodn02 = 786432, i = 1					ii2 = 6, bjmodn02 =  262144, i = 0
		ii3 = 0, bjmodn03 = 917504, i = 1					ii3 = 4, bjmodn03 =  786432, i = 0
		ii4 = 6, bjmodn04 = 131072, i = 0					ii4 = 2, bjmodn04 = 1310720, i = 0
		ii5 = 5, bjmodn05 = 262144, i = 0					ii5 = 0, bjmodn05 = 1835008, i = 1
		ii6 = 4, bjmodn06 = 393216, i = 0					ii6 = 5, bjmodn06 =  524288, i = 0
		J = 6:												J = 6:
		ii0 = 1, bjmodn00 = 786432, i = 1					ii0 = 1, bjmodn00 = 1572864, i = 1
		ii1 = 0, bjmodn01 = 917504, i = 1					ii1 = 6, bjmodn01 =  262144, i = 0
		ii2 = 6, bjmodn02 = 131072, i = 0					ii2 = 4, bjmodn02 =  786432, i = 0
		ii3 = 5, bjmodn03 = 262144, i = 0					ii3 = 2, bjmodn03 = 1310720, i = 0
		ii4 = 4, bjmodn04 = 393216, i = 0					ii4 = 0, bjmodn04 = 1835008, i = 1
		ii5 = 3, bjmodn05 = 524288, i = 0					ii5 = 5, bjmodn05 =  524288, i = 0
		ii6 = 2, bjmodn06 = 655360, i = 0					ii6 = 3, bjmodn06 = 1048576, i = 0
		J = 8:												J = 8:
		ii0 = 6, bjmodn00 = 131072, i = 0					ii0 = 6, bjmodn00 =  262144, i = 0
		ii1 = 5, bjmodn01 = 262144, i = 0					ii1 = 4, bjmodn01 =  786432, i = 0
		ii2 = 4, bjmodn02 = 393216, i = 0					ii2 = 2, bjmodn02 = 1310720, i = 0
		ii3 = 3, bjmodn03 = 524288, i = 0					ii3 = 0, bjmodn03 = 1835008, i = 1
		ii4 = 2, bjmodn04 = 655360, i = 0					ii4 = 5, bjmodn04 =  524288, i = 0
		ii5 = 1, bjmodn05 = 786432, i = 1					ii5 = 3, bjmodn05 = 1048576, i = 0
		ii6 = 0, bjmodn06 = 917504, i = 1					ii6 = 1, bjmodn06 = 1572864, i = 1
		J = 10:												J = 10:
		ii0 = 4, bjmodn00 = 393216, i = 0					ii0 = 4, bjmodn00 =  786432, i = 0
		ii1 = 3, bjmodn01 = 524288, i = 0					ii1 = 2, bjmodn01 = 1310720, i = 0
		ii2 = 2, bjmodn02 = 655360, i = 0					ii2 = 0, bjmodn02 = 1835008, i = 1
		ii3 = 1, bjmodn03 = 786432, i = 1					ii3 = 5, bjmodn03 =  524288, i = 0
		ii4 = 0, bjmodn04 = 917504, i = 1					ii4 = 3, bjmodn04 = 1048576, i = 0
		ii5 = 6, bjmodn05 = 131072, i = 0					ii5 = 1, bjmodn05 = 1572864, i = 1
		ii6 = 5, bjmodn06 = 262144, i = 0					ii6 = 6, bjmodn06 =  262144, i = 0
		J = 12:												J = 12:
		ii0 = 2, bjmodn00 = 655360, i = 0					ii0 = 2, bjmodn00 = 1310720, i = 0
		ii1 = 1, bjmodn01 = 786432, i = 1					ii1 = 0, bjmodn01 = 1835008, i = 1
		ii2 = 0, bjmodn02 = 917504, i = 1					ii2 = 5, bjmodn02 =  524288, i = 0
		ii3 = 6, bjmodn03 = 131072, i = 0					ii3 = 3, bjmodn03 = 1048576, i = 0
		ii4 = 5, bjmodn04 = 262144, i = 0					ii4 = 1, bjmodn04 = 1572864, i = 1
		ii5 = 4, bjmodn05 = 393216, i = 0					ii5 = 6, bjmodn05 =  262144, i = 0
		ii6 = 3, bjmodn06 = 524288, i = 0					ii6 = 4, bjmodn06 =  786432, i = 0
`	...And now (i.e. every [nwt]th pass) repeat the j=0 pattern: ...
		J = 14:												J = 14:
		ii0 = 0, bjmodn00 = 917504, i = 1					ii0 = 0, bjmodn00 = 1835008, i = 1
		ii1 = 6, bjmodn01 = 131072, i = 0					ii1 = 5, bjmodn01 =  524288, i = 0
		ii2 = 5, bjmodn02 = 262144, i = 0					ii2 = 3, bjmodn02 = 1048576, i = 0
		ii3 = 4, bjmodn03 = 393216, i = 0					ii3 = 1, bjmodn03 = 1572864, i = 1
		ii4 = 3, bjmodn04 = 524288, i = 0					ii4 = 6, bjmodn04 =  262144, i = 0
		ii5 = 2, bjmodn05 = 655360, i = 0					ii5 = 4, bjmodn05 =  786432, i = 0
		ii6 = 1, bjmodn06 = 786432, i = 1					ii6 = 2, bjmodn06 = 1310720, i = 0

	For the F24 case, the cyclical per-loop-pass index-pattern shift = 2; for F25 it = 1. How to compute it in advance?

	The per-loop-pass update of the bjmodn terms is this:

		i = (bjmodn > sw);					//       i = 1 if a bigword,   0 if a smallword
		bjmodn -= sw;						// result >= 0 if a bigword, < 0 if a smallword
		bjmodn += ( ((int)i-1) & n);		//       add 0 if a bigword,   N if a smallword

	So need to take bjmodn00 start value,which = n, subtract sw, see how many bjmodn's further on we need to go to get
	the resulting bjmodn00-sw value.

	Take advantage of this as follows:
	Init length-[nwt] integer arrays ii_arr = (initial values of ii0-ii6), i_arr = (initial values of (bjmodn00-06 > sw)),
	treat these as circular arrays, for j = 0 starting index into these arrays = 0,
	on each loop execution we advance the starting index in these arrays by wts_idx_incr := (bw*radix0/n) places.
	These ii and i-array values are used as lookups into suitably initialized weights and base arrays and arrays of their inverses.
	In fact if we use the initial patterns of the ii and i-indices (i.e. those corresponding to the j=0 initial loop interation)
	to init correspondingly permuted weights and base arrays, we can replace actual in-loop use of separate ii and i-arrays
	with a simple set of [nwt] indices which start with values [0,1,2,...,nwt-1] and get incremented by wts_idx_incr (mod nwt)
	on each pass.

Example 1: p=2^24, n = 917504 = 7*2^17, bjmodn00 - sw = bw = p%n = 262144 which corr. to bjmodn02 in the initial-bjmodn-values list,
	==> per-pass circular shift amount of the wts-array [and related] index = 2 := wts_idx_incr.

Example 2: p=2^25, n = 1835008 = 7*2^18, bjmodn00 - sw = bw = p%n = 524288 which corr. to bjmodn01 in the initial-bjmodn-values list,
	==> per-pass circular shift amount of the wts-array [and related] index = 1 := wts_idx_incr.

	The 2 indexing subarrays that come into play are:

	[1] ii: Define SW_DIV_N = sw*nwt/n = (n-bw)*nwt/n = 655360*7/917504 = [5*2^17]*7/[7*2^17] = nwt-wts_idx_incr = 5
		We need (nwt) starting values of ii, with the (j)th of these defined by

		For the above example, with radix0 = 28, we get

		ii[j] = j*(SW_DIV_N*(n/[FFT radix used for the carry step]/2) % nwt
			  = j*[nwt-wts_idx_incr]*[n/28/2] % nwt
			  = j*5*2^14 % 7
			  = j*20 % 7
			  = -j%7
			  = 0,6,5,4,3,2,1 .

		Latter is order in which access wt[] and wtinv[] arrays, so simply init length-nwt local versions of these arrays
		containing the weights ordered according to that index pattern.

	[2] i : See any of the fused final-iFFT-radix/carry/initial-FFT-radix pass routines to see how the
		parameter bjmodnini is inited...define it as

			bjmodnini := bw*(n/2)/[FFT radix used for the carry step] % n .

		For the above example, with radix0 = 28, we get

			bjmodnini = [2*2^17]*[7*2^16]/28 % n
			= 2^32 % n
			= 131072 ,

		with the (j)th value of this

			bjmodn[j] := j*bjmodini % n,

		and we start this value at n (rather that 0) for j = 0, so resulting index i into base and baseinv arrays comes out right:

			i [j] = (bjmodn[j] > sw) .

		For the above example, with radix0 = 28, we get

			i [j] = [1, (j*131072) for j=1...6] > sw = 1,0,0,0,0,0,1,

		so init length-nwt local versions of base[] and baseinv[] arrays which contain the entries ordered via this index pattern.
	*/

	#include "sse2_macro_gcc64.h"

	#if (OS_BITS == 64)
	  #ifdef USE_AVX
		#define GCC_ASM_FULL_INLINE  0	// 0 to use small-macros to assemble radix-28 DFTs, 1 to inline fuse macros as a few big blobs of asm
		#define USE_64BIT_ASM_STYLE  1
	  #else
		#define GCC_ASM_FULL_INLINE  0
		#define USE_64BIT_ASM_STYLE  1
	  #endif
	#else	// 32-bit mode
		#define GCC_ASM_FULL_INLINE  1
	#endif

	#if GCC_ASM_FULL_INLINE
	  #ifdef USE_ARM_V8_SIMD
		#error ARMv8 build supports only small-macro-based 28-DFTs!
	  #endif
		#include "radix28_ditN_cy_dif1_gcc64.h"
	#endif

#endif

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
		vec_dbl *s1p00r;
		vec_dbl *half_arr;
	#else
		double *s1p00r;
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
		double cy_dat[2*RADIX+4] __attribute__ ((__aligned__(8)));	// Enforce min-alignment of 8 bytes in 32-bit builds.
	};

#endif

/**************/

int radix28_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-28 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-28 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix7/8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix28_ditN_cy_dif1";
  #if !defined(MULTITHREAD) && defined(USE_SSE2)
	const int pfetch_dist = PFETCH_DIST;
  #endif
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
	static double wts_mult[2], inv_mult[2];	// Const wts-multiplier and 2*(its multiplicative inverse)
  #if !defined(MULTITHREAD) && (!defined(USE_SSE2) || (defined(USE_AVX) && !defined(USE_AVX512)))
	double wt_re,wt_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
  #endif
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	double wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
  #endif
	int NDIVR,i,j,j1,jt,jstart,jhi,full_pass,khi,l,outer;
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	int j2;
  #endif
  #ifndef MULTITHREAD
	int ntmp;
  #endif
  #ifdef USE_SSE2
	int nbytes;
  #endif
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	int jp;
  #endif
  #if !defined(MULTITHREAD) && defined(USE_SSE2) && !defined(USE_AVX)
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = --|--|7 for avx512,avx,sse, respectively (This radix unavailable
	// for avx-512 builds; for avx/avx2 use old-style non-parametrized carry-macro loop, 3 x 8-way follwed by 1 x 4-way, sum = 28).
	// But fixed-incr too restrictive here, so in sse2 case 'divide 7 into pieces' via increment-array whose elts sum to 7:
	const int *incr,*inc_arr;
  // recurrence length 7 is 'long' for radix-56, so make long = medium for radix-28:
   #ifdef USE_ARM_V8_SIMD	// Have no specialized HIACC carry macro in ARMv8 SIMD, so this gets an "goes to 11" in LOACC mode via an incr_hiacc[] array:
	const int incr_long[] = {7}, incr_med[] = {7}, incr_short[] = {3,4}, incr_hiacc[] = {2,3,2};
   #else
	const int incr_long[] = {7}, incr_med[] = {7}, incr_short[] = {3,4}, incr_hiacc[] = {0};
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
  #endif // !MULTITHREAD && USE_SSE2 && !USE_AVX

	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p08,p12,p16,p20,p24, nsave = 0;
	static int poff[RADIX>>2];	// Store [RADIX/4] mults of p04 offset for loop control
  #if !defined(MULTITHREAD) && !defined(USE_SSE2)
	static int p0123[4];
  #endif
	static double radix_inv, n2inv;
#ifdef USE_SSE2
	uint32 nwt16 = nwt << L2_SZ_VD;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code
#endif
	// Need these both in scalar mode and to ease the SSE2-array init...dimension = ODD_RADIX;
	// In order to ease the ptr-access for the || routine, lump these 4*ODD_RADIX doubles together with copies of
	// the 4 in the passed-in bs[2] and bsinv[2] arrays [and used in this 4-double form by the mersenne-mod carry macros]
	// into a single foo_array[4*(ODD_RADIX+1)], then convert what used to be disparate ODD_RADIX-sized arrays to pointers.
	//
	// NOTE: Use literal value of ODD_RADIX in dimensioning in order to allow us to explicitly make static ...
	// some compilers [e.g. gcc] are OK with foo_array[(ODD_RADIX+1)<<2], others [e.g. icc] will complain like
	//          "error: a variable length array cannot have static storage duration"
	// even though ODD_RADIX has been declared a const:
	static double foo_array[(ODD_RADIX+1)<<2], *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr;
#if !defined(USE_SSE2) && !defined(MULTITHREAD)
	static double bs,bsinv;
#endif

// FMA-based SIMD or (scalar-double) + (LO_ADD = 1 in masterdefs.h) means non-Nussbaumer radix-7, uses these sincos constants:
#if (defined(USE_AVX2) || defined(USE_ARM_V8_SIMD)) || (!defined(MULTITHREAD) && !defined(USE_SSE2) && (LO_ADD != 0))
	const double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
#elif (defined(USE_SSE2) && !(defined(USE_AVX2) || defined(USE_ARM_V8_SIMD))) || (!defined(USE_SSE2) && !LO_ADD)	// Low-MUL/high-ADD Nussbaumer radix-7 needs these trig combos:
	#warning LO_ADD = 0 defined at compile time ... using low-mul Nussbaumer-style 7-DFT. [NOTE: SSE2-SIMD build *requires* this]
	// SSE2 version assumes LO_ADD = 0, i.e. the low-mul Nussbaumer-style DFT implementation:
	const double	cx0 =-0.16666666666666666667,	/* (cc1+cc2+cc3)/3 */
				 	cx1 = 1.52445866976115265675, 	/*  cc1-cc3		*/
				 	cx2 = 0.67844793394610472196, 	/*  cc2-cc3		*/
				 	cx3 = 0.73430220123575245957,	/* (cc1+cc2-2*cc3)/3	*/
				/* Switch the sign of ss3 in these: */
				 	sx0 = 0.44095855184409843174,	/* (ss1+ss2-ss3)/3	*/
				 	sx1 = 1.21571522158558792920, 	/*  ss1+ss3		*/
				 	sx2 = 1.40881165129938172752, 	/*  ss2+ss3		*/
				 	sx3 = 0.87484229096165655224;	/* (ss1+ss2+2*ss3)/3	*/
#endif
	double scale,dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex t[RADIX];
  #ifndef MULTITHREAD
	double *addr;
   #ifndef USE_SSE2
	double *addi;
	double t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13;
   #endif
	int *itmp;	// Pointer into the bjmodn array
   #if defined(USE_AVX) && !defined(USE_AVX512)
	int *itm2;	// Pointer into the bjmodn array
   #endif
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
   #ifdef USE_AVX
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
   #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
   #endif
   #if !defined(USE_SSE2) || defined(USE_AVX) && !defined(USE_AVX512)
	double rt,it;
   #endif
  #endif
	static int ii[ODD_RADIX] = {-1,-1,-1,-1,-1,-1,-1};	/* indices into weights arrays (mod NWT) */
	/* These are used in conjunction with the langth-ODD_RADIX arrays in the USE_SCALAR_CARRY code flow;
	In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
  #if !defined(MULTITHREAD) && defined(USE_SSE2) && !defined(USE_AVX)
	static int idx_offset, idx_incr;
  #endif
	static int wts_idx_incr = 0;
  #ifdef USE_SSE2
	static int wts_idx_inc2 = 0;
  #endif
	static int icycle[ODD_RADIX];
  #if !defined(MULTITHREAD) && !defined(USE_AVX512)
	static int ic_idx;
  #endif
#ifdef USE_SSE2
	static int jcycle[ODD_RADIX];
  #if !defined(MULTITHREAD) && !defined(USE_AVX512)
	static int jc_idx;
  #endif
  #ifdef USE_AVX
	static int kcycle[ODD_RADIX];
	static int lcycle[ODD_RADIX];
   #if !defined(MULTITHREAD) && !defined(USE_AVX512)
	static int kc_idx, lc_idx;
   #endif
  #endif

	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else

	double *add0, *add1, *add2, *add3;	/* Addresses into array sections */

  #endif

  #ifndef MULTITHREAD
	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
  #endif
	const double crnd = 3.0*0x4000000*0x2000000;
	static vec_dbl *one,*two, *cc0, *ss0, *cc1, *ss1, *cc2, *ss2, *cc3, *ss3, *max_err, *sse2_rnd, *half_arr
	,*s1p00r
  #if !defined(MULTITHREAD) && !GCC_ASM_FULL_INLINE
	,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r
	,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r
  #endif
		,*cy_r	// Need RADIX slots for sse2 carries, RADIX/2 for avx
  #if !defined(MULTITHREAD) && defined(USE_AVX)
		,*cy_i	// Need RADIX slots for sse2 carries, RADIX/2 for avx
  #endif
	;
  #ifdef USE_AVX
	static vec_dbl *base_negacyclic_root;
  #endif
	vec_dbl *tmp,*tm2;	// Non-static utility ptrs
  #ifndef MULTITHREAD
   #if defined(USE_AVX) && !defined(USE_AVX512)
	vec_dbl *tm0;	// Non-static utility ptrs
   #endif
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
	static task_control_t   task_control = {NULL, (void*)cy28_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
  #if PFETCH
	double *addp;
  #endif
	int bjmodn[RADIX];
	double re,im,temp,frac
		,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r
		,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i
		,cy_r[RADIX],cy_i[RADIX];

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

  #ifdef USE_AVX512	// This includes IMCI512 builds
	WARN(HERE, "radix28_ditN_cy_dif1: No AVX-512 support; Skipping this leading radix.", "", 1); return(ERR_RADIX0_UNAVAILABLE);
  #endif

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
		ASSERT(LO_ADD,"LO_ADD");
		psave = p;	nsave = n;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)RADIX));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;	/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

	  #ifdef USE_AVX	// AVX LOACC: Make CARRY_8_WAY the default for this mode
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
		// consisting of radix28_creals_in_local_store vec_dbl and (8+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix28_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix28_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 56 16-byte slots of sc_arr for temporaries, next 8 for the nontrivial complex 16th roots,
	next 28 for the doubled carry pairs, next 2 for ROE and RND_CONST, next RADIX for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array. We also use a few more slots in AVX mode
	for the compact negacyclic-roots chained-multiply scheme, which drastically reduces accesses to the rn0,rn1 tables.
	*/
	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = sc_ptr;
		s1p00r = tmp + 0x00;
	  #if !defined(MULTITHREAD) && !GCC_ASM_FULL_INLINE
		s1p01r = tmp + 0x02;
		s1p02r = tmp + 0x04;
		s1p03r = tmp + 0x06;
		s1p04r = tmp + 0x08;
		s1p05r = tmp + 0x0a;
		s1p06r = tmp + 0x0c;
		s1p07r = tmp + 0x0e;
		s1p08r = tmp + 0x10;
		s1p09r = tmp + 0x12;
		s1p10r = tmp + 0x14;
		s1p11r = tmp + 0x16;
		s1p12r = tmp + 0x18;
		s1p13r = tmp + 0x1a;
		s1p14r = tmp + 0x1c;
		s1p15r = tmp + 0x1e;
		s1p16r = tmp + 0x20;
		s1p17r = tmp + 0x22;
		s1p18r = tmp + 0x24;
		s1p19r = tmp + 0x26;
		s1p20r = tmp + 0x28;
		s1p21r = tmp + 0x2a;
		s1p22r = tmp + 0x2c;
		s1p23r = tmp + 0x2e;
		s1p24r = tmp + 0x30;
		s1p25r = tmp + 0x32;
		s1p26r = tmp + 0x34;
		s1p27r = tmp + 0x36;
	  #endif
		tmp += 0x38;
		two   = tmp + 0x0;	// FMA2 versions of various DFT macros assume consts 2.0,1.0 laid out thusly
		one   = tmp + 0x1;
		cc0   = tmp + 0x2;
		ss0   = tmp + 0x3;
		cc1   = tmp + 0x4;
		ss1   = tmp + 0x5;
		cc2   = tmp + 0x6;
		ss2   = tmp + 0x7;
		cc3   = tmp + 0x8;
		ss3   = tmp + 0x9;	/* Pad with extra 4 slots for scratch storage needed by SSE2_RADIX_07_DFT macro here */
		tmp += 0xe;
	  #ifdef USE_AVX
		cy_r = tmp;
	   #ifndef MULTITHREAD
		cy_i = tmp+0x7;
	   #endif
		tmp += 2*0x7;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays;
		max_err = tmp + 0x00;	// +2 for max_err, sse2_rnd
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x58; This is where the value of half_arr_offset56 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy_r = tmp;	/* cy_i = tmp+0xe; */	tmp += 2*0xe;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays; +60 = 368 complex
		max_err = tmp + 0x00;					// ^^^ sc_ptr += 0x126
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x66; This is where the value of half_arr_offset56 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*ODD_RADIX] x 16 for Fermat-mod */
	  #endif

		/* These remain fixed: */
		VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );	// 1.0,2.0 needed for FMA support
	  #if defined(USE_AVX2) || defined(USE_ARM_V8_SIMD)
		// AVX2 (i.e. FMA) means non-Nussbaumer radix-7, uses these sincos constants:
		VEC_DBL_INIT(cc0, uc1  );	VEC_DBL_INIT(ss0, us1);
		VEC_DBL_INIT(cc1, uc2  );	VEC_DBL_INIT(ss1, us2);
		VEC_DBL_INIT(cc2, uc3  );	VEC_DBL_INIT(ss2, us3);
		VEC_DBL_INIT(cc3, 0.0  );	VEC_DBL_INIT(ss3, 0.0);	// Unused in non-Nussbaumer mode
	  #else
		/* SSE2 version assumes LO_ADD = 0, i.e. the low-mul Nussbaumer-style DFT implementation: */
		/* cc0 = (cc1+cc2+cc3)/3 - 1; subtract 1 from Nussbaumer's definition in order to ease in-place computation */
		VEC_DBL_INIT(cc0, cx0-1);	VEC_DBL_INIT(ss0, sx0);
		VEC_DBL_INIT(cc1, cx1  );	VEC_DBL_INIT(ss1, sx1);
		VEC_DBL_INIT(cc2, cx2  );	VEC_DBL_INIT(ss2, sx2);
		VEC_DBL_INIT(cc3, cx3  );	VEC_DBL_INIT(ss3, sx3);
	  #endif
		VEC_DBL_INIT(sse2_rnd, crnd);		/* SSE2 math = 53-mantissa-bit IEEE double-float: */

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
	#ifdef USE_AVX

		base_negacyclic_root = half_arr + RADIX;

	  #if HIACC
		/*
		The pattern of the negacyclic-DWT-weights ("nDWTs") applied to the RADIX complex outputs of the final-DIT-pass is like so:
		The nDWTs multiplying each set of RADIX DIT DFT outputs are simply the product of a single complex-root "base multiplier" rbase
		(separately computed for each batch of DFT outputs), which "base root" multiplies the (0 - RADIX-1)st [4*RADIX]th roots of unity,
		i.e.
			 rbase * (j*I*Pi/2)/RADIX, for j = 0, ..., RADIX-1 .

		As applied to the 4-way-SIMD data in AVX mode, we will in fact have 4 such base roots at a time in a pair of YMM registers
		(real parts in one YMM, imaginary parts in another. The above [4*28]th roots will be precomputed and stored in 14 YMM-sized
		local-data slots (7 for the Re parts, 7 for the Im) like so, using the above j-index as a shorthand:

			Slot  0: Re[j =  0 + 0, 1, 2, 3]
			Slot  1: Im[j =  0 + 0, 1, 2, 3]

			Slot  2: Re[j =  4 + 0, 1, 2, 3]
			Slot  3: Im[j =  4 + 0, 1, 2, 3]

			Slot  4: Re[j =  8 + 0, 1, 2, 3]
			Slot  5: Im[j =  8 + 0, 1, 2, 3]

			Slot  6: Re[j = 12 + 0, 1, 2, 3]
			Slot  7: Im[j = 12 + 0, 1, 2, 3]

			Slot  8: Re[j = 16 + 0, 1, 2, 3]
			Slot  9: Im[j = 16 + 0, 1, 2, 3]

			Slot 10: Re[j = 20 + 0, 1, 2, 3]
			Slot 11: Im[j = 20 + 0, 1, 2, 3]

			Slot 12: Re[j = 24 + 0, 1, 2, 3]
			Slot 13: Im[j = 24 + 0, 1, 2, 3]

		Prior to performing the normalize-and-carry step on each set of 28 AVX-complex ( = 4*28 double-complex) DIT DFT outputs,
		we compute the 4 base multipliers needed for that set of data:

			Re[rbase 0, 1, 2, 3]
			Im[rbase 0, 1, 2, 3]

		and then do a complex multiply of that quartet of complex-double data with each of the above 7 precomputed AVX-complex
		constants, storing the results in another set of local-mem slots and/or YMM registers, as desired.
		*/
		// Init exp(j*I*Pi/2)/28, for j = 0-27:
		tmp = base_negacyclic_root + RADIX*2;	// First 56 = 2*RADIX slots reserved for RADIX/4 copies of the Re/Im parts of the 4 base multipliers
		tm2 = tmp + RADIX/2 - 1;
										tmp->d0 = 1.0;	(tmp+1)->d0 = 0.0;
		tmp64 = 0x3FEFF31CCABE6E4Cull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(01*I*Pi/56) = sin(27*I*Pi/56) */
		tmp64 = 0x3FEFCC7D8C64135Full;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(02*I*Pi/56) = sin(26*I*Pi/56) */
		tmp64 = 0x3FEF8C4160D38565ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(03*I*Pi/56) = sin(25*I*Pi/56) */
		tmp += 2;
		tmp64 = 0x3FEF329C0558E969ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(04*I*Pi/56) = sin(24*I*Pi/56) */	tm2 -= 2;
		tmp64 = 0x3FEEBFD5AEFD405Aull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(05*I*Pi/56) = sin(23*I*Pi/56) */
		tmp64 = 0x3FEE344AD05D3F86ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(06*I*Pi/56) = sin(22*I*Pi/56) */
		tmp64 = 0x3FED906BCF328D46ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(07*I*Pi/56) = sin(21*I*Pi/56) */
		tmp += 2;
		tmp64 = 0x3FECD4BCA9CB5C71ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(08*I*Pi/56) = sin(20*I*Pi/56) */	tm2 -= 2;
		tmp64 = 0x3FEC01D48CB95263ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(09*I*Pi/56) = sin(19*I*Pi/56) */
		tmp64 = 0x3FEB185D590D5A44ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(10*I*Pi/56) = sin(18*I*Pi/56) */
		tmp64 = 0x3FEA19131B8279C4ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(11*I*Pi/56) = sin(17*I*Pi/56) */
		tmp += 2;
		tmp64 = 0x3FE904C37505DE4Bull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(12*I*Pi/56) = sin(16*I*Pi/56) */	tm2 -= 2;
		tmp64 = 0x3FE7DC4CF5162385ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(13*I*Pi/56) = sin(15*I*Pi/56) */
		tmp64 = 0x3FE6A09E667F3BCDull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(14*I*Pi/56) = sin(14*I*Pi/56) */
		tmp64 = 0x3FE552B60F035F34ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(15*I*Pi/56) = sin(13*I*Pi/56) */
		tmp += 2;
		tmp64 = 0x3FE3F3A0E28BEDD1ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(16*I*Pi/56) = sin(12*I*Pi/56) */	tm2 -= 2;
		tmp64 = 0x3FE28479AA873CFEull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(17*I*Pi/56) = sin(11*I*Pi/56) */
		tmp64 = 0x3FE106682221CD8Aull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(18*I*Pi/56) = sin(10*I*Pi/56) */
		tmp64 = 0x3FDEF5401024C4F4ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(19*I*Pi/56) = sin(09*I*Pi/56) */
		tmp += 2;
		tmp64 = 0x3FDBC4C04D71ABC1ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(20*I*Pi/56) = sin(08*I*Pi/56) */	tm2 -= 2;
		tmp64 = 0x3FD87DE2A6AEA963ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(21*I*Pi/56) = sin(07*I*Pi/56) */
		tmp64 = 0x3FD5234ACA69A9FEull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(22*I*Pi/56) = sin(06*I*Pi/56) */
		tmp64 = 0x3FD1B7AC4AFC3C02ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(23*I*Pi/56) = sin(05*I*Pi/56) */
		tmp += 2;
		tmp64 = 0x3FCC7B90E3024582ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(24*I*Pi/56) = sin(04*I*Pi/56) */	tm2 -= 2;
		tmp64 = 0x3FC570D80B7C3350ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(25*I*Pi/56) = sin(03*I*Pi/56) */
		tmp64 = 0x3FBCA9B4332D6F61ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(26*I*Pi/56) = sin(02*I*Pi/56) */
		tmp64 = 0x3FACB544024FC940ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(27*I*Pi/56) = sin(01*I*Pi/56) */
		tmp = base_negacyclic_root + RADIX*2;	// reset to point to start of above block
		nbytes = RADIX*SZ_VD/2;	// 7 AVX-register-sized complex data

	  #else	// HIACC = false ==> lower-precision version:

		/* Alternatively, we
		could save on local storage by precomputing just

			Slot 0: Re[j =  0 + 0, 1, 2, 3]
			Slot 1: Im[j =  0 + 0, 1, 2, 3]

		and prior to performing the normalize-and-carry step on each set of 28 AVX-complex ( = 4*28 double-complex) DIT DFT outputs,
		cmultiplying this AVX-complex datum by the 4 base multipliers needed for the first set of 4 AVX-complex DFT outputs:

			Re[rbase 0, 1, 2, 3]
			Im[rbase 0, 1, 2, 3] .

		After processing each such set of 4 AVX-complex DFT outputs, we get ready for the next set of outputs by multiplying our
		current set of 4 AXV-complex roots by the following complex "up-multiply" constant:

			Slot 2: Re[j = 4, 4, 4, 4]
			Slot 3: Im[j = 4, 4, 4, 4] .

		This will have slightly more roundoff error, but is simpler and more storage-compact.
		*/
		// Init exp(j*I*Pi/2)/28, for j = 0-3:
		tmp = base_negacyclic_root + 8;	// First 8 slots reserved for Re/Im parts of the 4 base multipliers
		tmp->d0 = 1.0;
		tmp64 = 0x3FEFF31CCABE6E4Cull;	tmp->d1 = *(double *)&tmp64;	// cos(01*I*Pi/56)
		tmp64 = 0x3FEFCC7D8C64135Full;	tmp->d2 = *(double *)&tmp64;	// cos(02*I*Pi/56)
		tmp64 = 0x3FEF8C4160D38565ull;	tmp->d3 = *(double *)&tmp64;	// cos(03*I*Pi/56)
		(++tmp)->d0 = 0.0;
		tmp64 = 0x3FC570D80B7C3350ull;	tmp->d3 = *(double *)&tmp64;	// sin(03*I*Pi/120)
		tmp64 = 0x3FBCA9B4332D6F61ull;	tmp->d2 = *(double *)&tmp64;	// sin(02*I*Pi/120)
		tmp64 = 0x3FACB544024FC940ull;	tmp->d1 = *(double *)&tmp64;	// sin(01*I*Pi/120)
		++tmp;
		tmp64 = 0x3FEF329C0558E969ull;	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// cos(04*I*Pi/56)
		++tmp;
		tmp64 = 0x3FCC7B90E3024582ull;	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// sin(04*I*Pi/56)
		tmp = base_negacyclic_root + 8;	// reset to point to start of above block
		nbytes = 4*SZ_VD;	// 2 AVX-register-sized complex data

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
	#ifdef USE_AVX

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

#ifdef USE_AVX
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

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );

	#if !defined(MULTITHREAD) && !defined(USE_SSE2)
		p0123[0] = 0; p0123[1] = p01; p0123[2] = p02; p0123[3] = p03;
	#endif

		poff[0] =   0; poff[1] = p04; poff[2] = p08; poff[3] = p12; poff[4] = p16; poff[5] = p20; poff[6] = p24;

		ASSERT(p01+p01 == p02, "p01+p01 != p02");
		ASSERT(p02+p02 == p04, "p02+p02 != p04");
		ASSERT(p04+p04 == p08, "p04+p04 != p08");
		ASSERT(p08+p04 == p12, "p08+p04 != p12");
		ASSERT(p12+p04 == p16, "p12+p04 != p16");
		ASSERT(p16+p04 == p20, "p16+p04 != p20");
		ASSERT(p20+p04 == p24, "p20+p04 != p24");

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

		  #ifdef USE_AVX

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
			}

		#endif	/* USE_SSE2 */
		}

	#ifdef USE_PTHREAD

		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			tdat[ithread].bjmodnini = _bjmodnini[CY_THREADS];
			tdat[ithread].bjmodn0 = _bjmodnini[ithread];
		#ifdef USE_SSE2
			tdat[ithread].s1p00r = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (vec_dbl *)((intptr_t)tdat[ithread].s1p00r + ((intptr_t)half_arr - (intptr_t)s1p00r));
		#else
			// In scalar mode use these 2 ptrs to pass wts_idx_incr and the base/baseinv/etc array-ptrs:
			tdat[ithread].s1p00r     = (double *)foo_array;
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

/*...The radix-28 final DIT pass is here.	*/

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
		}

		// Propagate the above inv-wts to the remaining threads - surrounding consts are unchanged:
		nbytes = ODD_RADIX*SZ_VD;
		tmp = half_arr + ODD_RADIX;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	#endif	/* USE_SSE2 */
	}	// 	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)

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

#endif	// USE_SSE2

	/* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
	if(!full_pass)
	{
		khi = 1;
	//	fprintf(dbg_file, "radix28_ditN_cy_dif1: Cleanup Pass:\n");
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
		ASSERT(tdat[ithread].s1p00r == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].s1p00r;
		ASSERT(((tmp + 0x38)->d0 == 2.0 && (tmp + 0x38)->d1 == 2.0), "thread-local memcheck failed!");
		ASSERT(((tmp + half_arr_offset28-1)->d0 == crnd && (tmp + half_arr_offset28-1)->d1 == crnd), "thread-local memcheck failed!");
	#endif
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX
			// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
			dtmp = (tmp + half_arr_offset28+40)->d0 * (tmp + half_arr_offset28+56)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp + half_arr_offset28+40)->d1 * (tmp + half_arr_offset28+56)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#elif defined(USE_SSE2)
			dtmp = (tmp + half_arr_offset28+10)->d0 * (tmp + half_arr_offset28+14)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp + half_arr_offset28+10)->d1 * (tmp + half_arr_offset28+14)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#endif
			/* init carries	*/
			for(i = 0; i < RADIX; i++) {
				tdat[ithread].cy_r[i] = _cy_r[i][ithread];
			}
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
		#ifdef USE_SSE2
			dtmp = (tmp + half_arr_offset28)->d0 * (tmp + half_arr_offset28+ODD_RADIX)->d0;	ASSERT(fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp + half_arr_offset28)->d1 * (tmp + half_arr_offset28+ODD_RADIX)->d1;	ASSERT(fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
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
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
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
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
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
		#include "radix28_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
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
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
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
		ASSERT(0x0 == cy28_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
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
		j_jhi =15;
	else
		j_jhi = 7;

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

/***************/

void radix28_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-28 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix7_dif_pass for details on the radix-15 subtransforms.
*/
	int j,j1,j2;
	static int n28,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27, first_entry=TRUE;
	static double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	double re,im,rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i;

	if(!first_entry && (n/28) != n28)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n28=n/28;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n28;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p10 = p09 + p01;
		p11 = p10 + p01;
		p12 = p11 + p01;
		p13 = p12 + p01;
		p14 = p13 + p01;
		p15 = p14 + p01;
		p16 = p15 + p01;
		p17 = p16 + p01;
		p18 = p17 + p01;
		p19 = p18 + p01;
		p20 = p19 + p01;
		p21 = p20 + p01;
		p22 = p21 + p01;
		p23 = p22 + p01;
		p24 = p23 + p01;
		p25 = p24 + p01;
		p26 = p25 + p01;
		p27 = p26 + p01;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p17 = p17 + ( (p17 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p19 = p19 + ( (p19 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p21 = p21 + ( (p21 >> DAT_BITS) << PAD_BITS );
		p22 = p22 + ( (p22 >> DAT_BITS) << PAD_BITS );
		p23 = p23 + ( (p23 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p25 = p25 + ( (p25 >> DAT_BITS) << PAD_BITS );
		p26 = p26 + ( (p26 >> DAT_BITS) << PAD_BITS );
		p27 = p27 + ( (p27 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-28 pass is here.	*/

	for(j = 0; j < n28; j += 2)
	{
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 += ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

	/*...gather the needed data (28 64-bit complex, i.e. 56 64-bit reals) and do 4 radix-7 transforms...*/

	/*
		Twiddleless version requires us to swap inputs as follows:
		indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27
			  -> 0,21,14, 7,24,17,10, 3,20,13, 6,27,16, 9, 2,23,12, 5,26,19, 8, 1,22,15, 4,25,18,11
		I.e. start out with first septet of indices {0,4,8,12,16,20,24}, permute those according to
		{0,4,8,12,16,20,24}*27%28 = {0,24,20,16,12,8,4}, then each is head of a length-4 list of indices with decrement 7.
	*/
						 /*                                                                      inputs                                                                         */ /*               intermediates              */ /*                                                  outputs                                                  */
		RADIX_07_DFT(a[j1    ],a[j2    ],a[j1+p24],a[j2+p24],a[j1+p20],a[j2+p20],a[j1+p16],a[j2+p16],a[j1+p12],a[j2+p12],a[j1+p08],a[j2+p08],a[j1+p04],a[j2+p04],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a[j1+p21],a[j2+p21],a[j1+p17],a[j2+p17],a[j1+p13],a[j2+p13],a[j1+p09],a[j2+p09],a[j1+p05],a[j2+p05],a[j1+p01],a[j2+p01],a[j1+p25],a[j2+p25],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a[j1+p14],a[j2+p14],a[j1+p10],a[j2+p10],a[j1+p06],a[j2+p06],a[j1+p02],a[j2+p02],a[j1+p26],a[j2+p26],a[j1+p22],a[j2+p22],a[j1+p18],a[j2+p18],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a[j1+p07],a[j2+p07],a[j1+p03],a[j2+p03],a[j1+p27],a[j2+p27],a[j1+p23],a[j2+p23],a[j1+p19],a[j2+p19],a[j1+p15],a[j2+p15],a[j1+p11],a[j2+p11],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);

	/*...and now do 7 radix-4 transforms...*/
						 /*                          inputs                           */ /*                                      outputs                                      */
		RADIX_04_DIF(a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],rt,it);
		RADIX_04_DIF(a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,a[j1+p25],a[j2+p25],a[j1+p24],a[j2+p24],a[j1+p27],a[j2+p27],a[j1+p26],a[j2+p26],rt,it);
		RADIX_04_DIF(a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a[j1+p22],a[j2+p22],a[j1+p23],a[j2+p23],a[j1+p21],a[j2+p21],a[j1+p20],a[j2+p20],rt,it);
		RADIX_04_DIF(a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p18],a[j2+p18],a[j1+p19],a[j2+p19],rt,it);
		RADIX_04_DIF(a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],rt,it);
		RADIX_04_DIF(a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p11],a[j2+p11],a[j1+p10],a[j2+p10],rt,it);
		RADIX_04_DIF(a1p24r,a1p24i,a1p25r,a1p25i,a1p26r,a1p26i,a1p27r,a1p27i,a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],rt,it);
	}
}

/***************/

void radix28_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-28 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2;
	static int n28,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27, first_entry=TRUE;
	static double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	double re,im,rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i;

	if(!first_entry && (n/28) != n28)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n28=n/28;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n28;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p10 = p09 + p01;
		p11 = p10 + p01;
		p12 = p11 + p01;
		p13 = p12 + p01;
		p14 = p13 + p01;
		p15 = p14 + p01;
		p16 = p15 + p01;
		p17 = p16 + p01;
		p18 = p17 + p01;
		p19 = p18 + p01;
		p20 = p19 + p01;
		p21 = p20 + p01;
		p22 = p21 + p01;
		p23 = p22 + p01;
		p24 = p23 + p01;
		p25 = p24 + p01;
		p26 = p25 + p01;
		p27 = p26 + p01;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p17 = p17 + ( (p17 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p19 = p19 + ( (p19 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p21 = p21 + ( (p21 >> DAT_BITS) << PAD_BITS );
		p22 = p22 + ( (p22 >> DAT_BITS) << PAD_BITS );
		p23 = p23 + ( (p23 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p25 = p25 + ( (p25 >> DAT_BITS) << PAD_BITS );
		p26 = p26 + ( (p26 >> DAT_BITS) << PAD_BITS );
		p27 = p27 + ( (p27 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-28 pass is here.	*/

	for(j=0; j < n28; j += 2)
	{
	#ifdef USE_AVX
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
		indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27
			  -> 0,24,20,16,12, 8, 4,21,17,13, 9, 5, 1,25,14,10, 6, 2,26,22,18, 7, 3,27,23,19,15,11

		I.e. start out with first quartet of indices {0,7,14,21}, permute those according to
		  {0,7,14,21}*27%28 = {0,21,14,7}, then each is head of a length-7 list of indices with decrement 4

		Remember, inputs to DIT are bit-reversed, so
		a[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27] contain
		x[0,14, 7,21, 1,15, 8,22, 2,16, 9,23, 3,17,10,24, 4,18,11,25, 5,19,12,26, 6,20,13,27], which get swapped to
	(Look at first nontrivial one, i.e. x[1 -> 24] ... in terms of a[] this translates to a[4 -> 15])
		x[0,14,21, 7,24,10,17, 3,20, 6,13,27,16, 2, 9,23,12,26, 5,19, 8,22, 1,15, 4,18,25,11], which means the a-indices get swapped as
		a[0, 1, 3, 2,15,14,13,12,25,24,26,27, 9, 8,10,11,22,23,20,21, 6, 7, 4, 5,16,17,19,18].
	*/
	/* EWM: 10/18/04: We swap the odd-index outputs of each of the radix-4 DIT transforms (1<=>3, 5<=>7, etc.) so that the indexing
					  of all the radix-7 transforms (really just the 2nd and 4th of these) winds up being in-place. This allows us
					  to properly re-use the ajp1 variables in the carry-pass version of this routine.
					  (NOTE: Didn't bother to do the swapping in the macro-less version of this routine - would only need to
					   do the swapping there if for some reason we wanted the macro-less DIT code in the carry routines.)
	*/
	/*...gather the needed data (28 64-bit complex, i.e. 56 64-bit reals) and do 7 radix-4 transforms...*/
					/*                                    inputs                                  */ /*                         outputs                   */
		RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a1p00r,a1p00i,a1p03r,a1p03i,a1p02r,a1p02i,a1p01r,a1p01i,rt,it);
		RADIX_04_DIT(a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p13],a[j2+p13],a[j1+p12],a[j2+p12],a1p04r,a1p04i,a1p07r,a1p07i,a1p06r,a1p06i,a1p05r,a1p05i,rt,it);
		RADIX_04_DIT(a[j1+p25],a[j2+p25],a[j1+p24],a[j2+p24],a[j1+p26],a[j2+p26],a[j1+p27],a[j2+p27],a1p08r,a1p08i,a1p11r,a1p11i,a1p10r,a1p10i,a1p09r,a1p09i,rt,it);
		RADIX_04_DIT(a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a1p12r,a1p12i,a1p15r,a1p15i,a1p14r,a1p14i,a1p13r,a1p13i,rt,it);
		RADIX_04_DIT(a[j1+p22],a[j2+p22],a[j1+p23],a[j2+p23],a[j1+p20],a[j2+p20],a[j1+p21],a[j2+p21],a1p16r,a1p16i,a1p19r,a1p19i,a1p18r,a1p18i,a1p17r,a1p17i,rt,it);
		RADIX_04_DIT(a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05],a1p20r,a1p20i,a1p23r,a1p23i,a1p22r,a1p22i,a1p21r,a1p21i,rt,it);
		RADIX_04_DIT(a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],a1p24r,a1p24i,a1p27r,a1p27i,a1p26r,a1p26i,a1p25r,a1p25i,rt,it);

	/*...and now do 4 radix-4 transforms...*/
					/*                                              inputs                                          */ /*               intermediates              */ /*                                                                     outputs                                                           */
		RADIX_07_DFT(a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a[j1    ],a[j2    ],a[j1+p08],a[j2+p08],a[j1+p16],a[j2+p16],a[j1+p24],a[j2+p24],a[j1+p04],a[j2+p04],a[j1+p12],a[j2+p12],a[j1+p20],a[j2+p20],uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a[j1+p07],a[j2+p07],a[j1+p15],a[j2+p15],a[j1+p23],a[j2+p23],a[j1+p03],a[j2+p03],a[j1+p11],a[j2+p11],a[j1+p19],a[j2+p19],a[j1+p27],a[j2+p27],uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a[j1+p14],a[j2+p14],a[j1+p22],a[j2+p22],a[j1+p02],a[j2+p02],a[j1+p10],a[j2+p10],a[j1+p18],a[j2+p18],a[j1+p26],a[j2+p26],a[j1+p06],a[j2+p06],uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a[j1+p21],a[j2+p21],a[j1+p01],a[j2+p01],a[j1+p09],a[j2+p09],a[j1+p17],a[j2+p17],a[j1+p25],a[j2+p25],a[j1+p05],a[j2+p05],a[j1+p13],a[j2+p13],uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy28_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
	  // LO_ADD = 1 in masterdefs.h means non-Nussbaumer radix-7, uses these sincos constants:
	  #if !defined(USE_SSE2) && (LO_ADD != 0)
		const double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
						us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
						uc2 =-.22252093395631440426,	 /* cos(2u)	*/
						us2 = .97492791218182360702,	 /* sin(2u)	*/
						uc3 =-.90096886790241912622,	 /* cos(3u)	*/
						us3 = .43388373911755812050;	 /* sin(3u)	*/
	  #elif !defined(USE_SSE2) && !LO_ADD	// Low-MUL/high-ADD Nussbaumer radix-7 needs these trig combos:
		const double	cx0 =-0.16666666666666666667,	/* (cc1+cc2+cc3)/3 */
						cx1 = 1.52445866976115265675, 	/*  cc1-cc3		*/
						cx2 = 0.67844793394610472196, 	/*  cc2-cc3		*/
						cx3 = 0.73430220123575245957,	/* (cc1+cc2-2*cc3)/3	*/
					/* Switch the sign of ss3 in these: */
						sx0 = 0.44095855184409843174,	/* (ss1+ss2-ss3)/3	*/
						sx1 = 1.21571522158558792920, 	/*  ss1+ss3		*/
						sx2 = 1.40881165129938172752, 	/*  ss2+ss3		*/
						sx3 = 0.87484229096165655224;	/* (ss1+ss2+2*ss3)/3	*/
	  #elif defined(USE_SSE2) && !defined(USE_AVX2) && !defined(USE_ARM_V8_SIMD)
		const double	sx3 = 0.87484229096165655224;	/* (ss1+ss2+2*ss3)/3	*/
	  #endif
	  #ifdef USE_SSE2
		const int pfetch_dist = PFETCH_DIST;
	  #endif
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p01,p02,p03,p04,p08,p12,p16,p20,p24;
		int poff[RADIX>>2];	// Store [RADIX/4] mults of p04 offset for loop control
	  #ifndef USE_SSE2
		int p0123[4];
	  #endif
		int j,j1,l,ntmp;
	  #ifndef USE_SSE2
		int j2;
	  #endif
	  #ifndef USE_AVX
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
	  #endif
	#if defined(USE_SSE2) && !defined(USE_AVX)
		// incr = Carry-chain wts-multipliers recurrence length, which must divide
		// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = --|--|7 for avx512,avx,sse, respectively (This radix unavailable
		// for avx-512 builds; for avx/avx2 use old-style non-parametrized carry-macro loop, 3 x 8-way follwed by 1 x 4-way, sum = 28).
		// But fixed-incr too restrictive here, so in sse2 case 'divide 7 into pieces' via increment-array whose elts sum to 7:
		const int *incr,*inc_arr;
	  // recurrence length 7 is 'long' for radix-56, so make long = medium for radix-28:
	  #ifdef USE_ARM_V8_SIMD	// Have no specialized HIACC carry macro in ARMv8 SIMD, so this gets an "goes to 11" in LOACC mode via an incr_hiacc[] array:
		const int incr_long[] = {7}, incr_med[] = {7}, incr_short[] = {3,4}, incr_hiacc[] = {2,3,2};
	  #else
		const int incr_long[] = {7}, incr_med[] = {7}, incr_short[] = {3,4}, incr_hiacc[] = {0};
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
	#endif // USE_SSE2 && !USE_AVX

	#ifdef USE_AVX
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#else
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	#endif
		double *addr,*addi;
	#if !defined(USE_SSE2) || (defined(USE_AVX) && !defined(USE_AVX512))
		double wt_re,wt_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	#endif
	#ifndef USE_SSE2
		double wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	#endif
	#if !defined(USE_SSE2) || (defined(USE_AVX) && !defined(USE_AVX512))
		double rt,it;
	#endif

	#ifdef USE_SSE2

		const double crnd = 3.0*0x4000000*0x2000000;
		double *add0, *add1, *add2, *add3;
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *tmp,*tm1,*tm2;	// utility ptrs
	  #if defined(USE_AVX) && !defined(USE_AVX512)
		vec_dbl *tm0;	// utility ptrs
	  #endif
		int *itmp;			// Pointer into the bjmodn array
	  #if defined(USE_AVX) && !defined(USE_AVX512)
		int *itm2;			// Pointer into the bjmodn array
	  #endif
	  #ifndef USE_AVX
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	  #endif
		vec_dbl /* *one, */*two, *cc0/*, *ss0, *cc1, *ss1, *cc2, *ss2, *cc3 */, *ss3, *max_err, *sse2_rnd, *half_arr,
			*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,
			*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r,
			*cy_r	// Need RADIX slots for sse2 carries, RADIX/2 for avx
	  #ifdef USE_AVX
			,*cy_i	// Need RADIX slots for sse2 carries, RADIX/2 for avx
	  #endif
			;
	  #if defined(USE_AVX) && !defined(USE_AVX512)
		vec_dbl *base_negacyclic_root;
	  #endif

		/* These are used in conjunction with the langth-ODD_RADIX arrays in the USE_SCALAR_CARRY #define below;
		In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
	  #ifndef USE_AVX
		int idx_offset, idx_incr;
	  #endif
		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else

		double *base, *baseinv, *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr, bs,bsinv;
		int wts_idx_incr;
		const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int jt,jp,m,m2;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	  #if PFETCH
		double *addp;
	  #endif
		int bjmodn[RADIX];	// Thread only carries a base datum here, must alloc a local array for remaining values
		double *cy_r = thread_arg->cy_r,*cy_i = thread_arg->cy_i, temp,frac;
		// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
		// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
		// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
		double re,im,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13	// tmps for scalar-double DFT macros
			,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r
			,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i;
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
	  #ifndef USE_AVX512
		struct complex *rn0 = thread_arg->rn0;
		struct complex *rn1 = thread_arg->rn1;
	  #endif

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

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );

	#ifndef USE_SSE2
		p0123[0] = 0; p0123[1] = p01; p0123[2] = p02; p0123[3] = p03;
	#endif

		poff[0] =   0; poff[1] = p04; poff[2] = p08; poff[3] = p12; poff[4] = p16; poff[5] = p20; poff[6] = p24;

	#ifdef USE_SSE2
		uint32 nwt16 = nwt << L2_SZ_VD;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code

		tmp = thread_arg->s1p00r;
		s1p00r = tmp + 0x00;
		s1p01r = tmp + 0x02;
		s1p02r = tmp + 0x04;
		s1p03r = tmp + 0x06;
		s1p04r = tmp + 0x08;
		s1p05r = tmp + 0x0a;
		s1p06r = tmp + 0x0c;
		s1p07r = tmp + 0x0e;
		s1p08r = tmp + 0x10;
		s1p09r = tmp + 0x12;
		s1p10r = tmp + 0x14;
		s1p11r = tmp + 0x16;
		s1p12r = tmp + 0x18;
		s1p13r = tmp + 0x1a;
		s1p14r = tmp + 0x1c;
		s1p15r = tmp + 0x1e;
		s1p16r = tmp + 0x20;
		s1p17r = tmp + 0x22;
		s1p18r = tmp + 0x24;
		s1p19r = tmp + 0x26;
		s1p20r = tmp + 0x28;
		s1p21r = tmp + 0x2a;
		s1p22r = tmp + 0x2c;
		s1p23r = tmp + 0x2e;
		s1p24r = tmp + 0x30;
		s1p25r = tmp + 0x32;
		s1p26r = tmp + 0x34;
		s1p27r = tmp + 0x36;
		tmp += 0x38;
		two   = tmp + 0x0;	// FMA2 versions of various DFT macros assume consts 2.0,1.0 laid out thusly
		//one   = tmp + 0x1;
		cc0   = tmp + 0x2;
		//ss0   = tmp + 0x3;
		//cc1   = tmp + 0x4;
		//ss1   = tmp + 0x5;
		//cc2   = tmp + 0x6;
		//ss2   = tmp + 0x7;
		//cc3   = tmp + 0x8;
		ss3   = tmp + 0x9;	/* Pad with extra 4 slots for scratch storage needed by SSE2_RADIX_07_DFT macro here */
		tmp += 0xe;
	  #ifdef USE_AVX
		cy_r = tmp;	cy_i = tmp+0x7;	tmp += 2*0x7;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 2 = 0x58 [avx] or 0x66 [sse2]; This is where the value of half_arr_offset28 comes from
		half_arr= tmp + 0x02;	/* This table needs 20x16 bytes */
	   #if defined(USE_AVX) && !defined(USE_AVX512)
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod
	   #endif
	  #else
		cy_r = tmp;	/* cy_i = tmp+0xe; */	tmp += 2*0xe;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 2 = 0x58 [avx] or 0x66 [sse2]; This is where the value of half_arr_offset28 comes from
		half_arr= tmp + 0x02;	/* This table needs 20x16 bytes */
	  #endif

		ASSERT((two->d0 == 2.0 && two->d1 == 2.0), "thread-local memcheck failed!");
	  #if defined(USE_AVX2) || defined(USE_ARM_V8_SIMD)
		// AVX2 (i.e. FMA)means non-Nussbaumer radix-7, uses these sincos constants:
		ASSERT((ss3->d0 == 0.0 && ss3->d1 == 0.0), "thread-local memcheck failed!");
	  #else
		/* SSE2 version assumes LO_ADD = 0, i.e. the low-mul Nussbaumer-style DFT implementation: */
		ASSERT((ss3->d0 == sx3 && ss3->d1 == sx3), "thread-local memcheck failed!");
	  #endif
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

		sign_mask = (uint64*)(s1p00r + radix28_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (  #doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		sse_sw  = sse_bw    + RE_IM_STRIDE;
		sse_n   = sse_sw    + RE_IM_STRIDE;
	  #ifdef USE_AVX
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
		base      = (double *)thread_arg->s1p00r;
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
		int *icycle = bjmodn;
	#ifndef USE_AVX512
		int ic_idx;
	#endif
	#ifdef USE_SSE2
		int wts_idx_inc2 = thread_arg->wts_idx_inc2;
		int *jcycle = icycle + ODD_RADIX;
	  #ifndef USE_AVX512
		int jc_idx;
	  #endif
	  #ifdef USE_AVX
		int *kcycle = jcycle + ODD_RADIX;
		int *lcycle = kcycle + ODD_RADIX;
	   #ifndef USE_AVX512
		int kc_idx, lc_idx;
	   #endif
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
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
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
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
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
		#include "radix28_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			addr = thread_arg->cy_r;
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
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
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
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
