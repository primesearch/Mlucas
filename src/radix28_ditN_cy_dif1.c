/*******************************************************************************
*                                                                              *
*   (C) 1997-2017 by Ernst W. Mayer.                                           *
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
  #ifdef USE_AVX512
	#define PFETCH_DIST	64	// Feb 2017: Test on KNL point to this as best
  #elif defined(USE_AVX)
	#define PFETCH_DIST	32	// This seems to work best on my Haswell, even though 64 bytes seems more logical in AVX mode
  #else
	#define PFETCH_DIST	32
  #endif
#endif

/* Use for toggling higher-accuracy version of the twiddles computation */
//#define HIACC 0	<*** prefer to set via compile-time flag; default is FALSE [= LOACC]

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

	#include "sse2_macro.h"

	#if defined(COMPILER_TYPE_MSVC)

		#define GCC_ASM_FULL_INLINE  0	// Only small-macro form available under MSVC

	#else	/* GCC-style inline ASM: */

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
		  #if OS_BITS == 32
			#include "radix28_ditN_cy_dif1_gcc32.h"
		  #else
			#include "radix28_ditN_cy_dif1_gcc64.h"
		  #endif
		#endif

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

		int icycle0;
		int icycle1;
		int icycle2;
		int icycle3;
		int icycle4;
		int icycle5;
		int icycle6;

	#ifdef USE_SSE2
		int jcycle0;
		int jcycle1;
		int jcycle2;
		int jcycle3;
		int jcycle4;
		int jcycle5;
		int jcycle6;
	  #ifdef USE_AVX
		int kcycle0;
		int kcycle1;
		int kcycle2;
		int kcycle3;
		int kcycle4;
		int kcycle5;
		int kcycle6;

		int lcycle0;
		int lcycle1;
		int lcycle2;
		int lcycle3;
		int lcycle4;
		int lcycle5;
		int lcycle6;
	  #endif
	#endif

	// double data:
		double maxerr;
		double scale;

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
	#ifdef LOACC
		double *wts_mult, *inv_mult;
	#endif
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

		int bjmodn00;
		int bjmodn01;
		int bjmodn02;
		int bjmodn03;
		int bjmodn04;
		int bjmodn05;
		int bjmodn06;
		int bjmodn07;
		int bjmodn08;
		int bjmodn09;
		int bjmodn10;
		int bjmodn11;
		int bjmodn12;
		int bjmodn13;
		int bjmodn14;
		int bjmodn15;
		int bjmodn16;
		int bjmodn17;
		int bjmodn18;
		int bjmodn19;
		int bjmodn20;
		int bjmodn21;
		int bjmodn22;
		int bjmodn23;
		int bjmodn24;
		int bjmodn25;
		int bjmodn26;
		int bjmodn27;

		// Pad to make size a multiple of 64 bytes:
		double dpad0;
		double dpad1;
		double dpad2;

		/* carries: */
		double cy_r00;
		double cy_r01;
		double cy_r02;
		double cy_r03;
		double cy_r04;
		double cy_r05;
		double cy_r06;
		double cy_r07;
		double cy_r08;
		double cy_r09;
		double cy_r10;
		double cy_r11;
		double cy_r12;
		double cy_r13;
		double cy_r14;
		double cy_r15;
		double cy_r16;
		double cy_r17;
		double cy_r18;
		double cy_r19;
		double cy_r20;
		double cy_r21;
		double cy_r22;
		double cy_r23;
		double cy_r24;
		double cy_r25;
		double cy_r26;
		double cy_r27;

		double cy_i00;
		double cy_i01;
		double cy_i02;
		double cy_i03;
		double cy_i04;
		double cy_i05;
		double cy_i06;
		double cy_i07;
		double cy_i08;
		double cy_i09;
		double cy_i10;
		double cy_i11;
		double cy_i12;
		double cy_i13;
		double cy_i14;
		double cy_i15;
		double cy_i16;
		double cy_i17;
		double cy_i18;
		double cy_i19;
		double cy_i20;
		double cy_i21;
		double cy_i22;
		double cy_i23;
		double cy_i24;
		double cy_i25;
		double cy_i26;
		double cy_i27;
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
	const int pfetch_dist = PFETCH_DIST;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
#ifdef USE_SSE2
	const int sz_vd = sizeof(vec_dbl);
	// lg(sizeof(vec_dbl)):
  #ifdef USE_AVX
	const int l2_sz_vd = 5;
  #else
	const int l2_sz_vd = 4;
  #endif
#endif
  #ifdef LOACC
	static double wts_mult[2], inv_mult[2];	// Const wts-multiplier and 2*(its multiplicative inverse)
  #endif
	double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	int NDIVR,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,ntmp,outer,nbytes;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p08,p12,p16,p20,p24, nsave = 0;
	static int poff[RADIX>>2];	// Store [RADIX/4] mults of p04 offset for loop control
	static double radix_inv, n2inv;
#ifdef USE_SSE2
	uint32 nwt16 = nwt << l2_sz_vd;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code
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
	static double foo_array[(7+1)<<2], *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr, bs,bsinv;

// AVX2 (i.e. FMA) or (scalar-double) + (LO_ADD = 1 in masterdefs.h) means non-Nussbaumer radix-7, uses these sincos constants:
#if defined(USE_AVX2) || (!defined(USE_SSE2) && defined(LO_ADD))

	const double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/

#elif defined(USE_SSE2) || !defined (LO_ADD)	// Low-MUL/high-ADD Nussbaumer radix-7 needs these trig combos:

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
#else
	#error Unhandled combination of preprocessor flags!	// Just in case I ever 'decide' to leave some holes in the above PP logic tree
#endif
	double *addr;
	double scale
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55;
	double maxerr = 0.0;
	int *itmp;	// Pointer into the bjmodn array
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	static int n_div_nwt;
	int col,co2,co3,m,m2;
  #ifdef USE_AVX
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
  #endif
	double rt,it;
	int k1,k2;
	int ii0,ii1,ii2,ii3,ii4,ii5,ii6;	/* indices into weights arrays (mod NWT) */
	/* These are used in conjunction with the langth-7 arrays in the USE_SCALAR_CARRY code flow;
	In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
	static int idx_offset, idx_incr, wts_idx_incr = 0, wts_idx_inc2 = 0
		,icycle0,icycle1,icycle2,icycle3,icycle4,icycle5,icycle6;

#ifdef USE_SSE2

	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else

	double *add0, *add1, *add2, *add3;	/* Addresses into array sections */

  #endif

	const double crnd = 3.0*0x4000000*0x2000000;
	static vec_dbl *one,*two, *cc0, *ss0, *cc1, *ss1, *cc2, *ss2, *cc3, *ss3, *max_err, *sse2_rnd, *half_arr
	,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r
	,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r;
  #ifdef USE_AVX
	static vec_dbl *base_negacyclic_root;
  #endif
	double dtmp;
	vec_dbl *tmp,*tm1,*tm2;	// Non-static utility ptrs
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	static int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27;
	static vec_dbl *cy_r00,*cy_r04,*cy_r08,*cy_r12,*cy_r16,*cy_r20,*cy_r24,*cy_i00,*cy_i04,*cy_i08,*cy_i12,*cy_i16,*cy_i20,*cy_i24;
  #ifndef USE_AVX
	static vec_dbl *cy_r02,*cy_r06,*cy_r10,*cy_r14,*cy_r18,*cy_r22,*cy_r26,*cy_i02,*cy_i06,*cy_i10,*cy_i14,*cy_i18,*cy_i22,*cy_i26;
  #endif

	/* In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
	static int jcycle0, jcycle1, jcycle2, jcycle3, jcycle4, jcycle5, jcycle6;
  #ifdef USE_AVX
	static int kcycle0, kcycle1, kcycle2, kcycle3, kcycle4, kcycle5, kcycle6;
	static int lcycle0, lcycle1, lcycle2, lcycle3, lcycle4, lcycle5, lcycle6;
  #endif

#endif

#ifdef MULTITHREAD

	static struct cy_thread_data_t *tdat = 0x0;
	// Threadpool-based dispatch stuff:
	static int main_work_units = 0, pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)cy28_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
  #if PFETCH
	double *addp;
  #endif
	int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27;
	double re,im,temp,frac
		,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r
		,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i
		,cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27
		,cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_bjmodn00 = 0x0,*_bjmodn01 = 0x0,*_bjmodn02 = 0x0,*_bjmodn03 = 0x0,*_bjmodn04 = 0x0,*_bjmodn05 = 0x0,*_bjmodn06 = 0x0,*_bjmodn07 = 0x0,*_bjmodn08 = 0x0,*_bjmodn09 = 0x0,*_bjmodn10 = 0x0,*_bjmodn11 = 0x0,*_bjmodn12 = 0x0,*_bjmodn13 = 0x0,*_bjmodn14 = 0x0,*_bjmodn15 = 0x0,*_bjmodn16 = 0x0,*_bjmodn17 = 0x0,*_bjmodn18 = 0x0,*_bjmodn19 = 0x0,*_bjmodn20 = 0x0,*_bjmodn21 = 0x0,*_bjmodn22 = 0x0,*_bjmodn23 = 0x0,*_bjmodn24 = 0x0,*_bjmodn25 = 0x0,*_bjmodn26 = 0x0,*_bjmodn27 = 0x0;
	static int *_bjmodnini = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double
	*_cy_r00 = 0x0,*_cy_r01 = 0x0,*_cy_r02 = 0x0,*_cy_r03 = 0x0,*_cy_r04 = 0x0,*_cy_r05 = 0x0,*_cy_r06 = 0x0,*_cy_r07 = 0x0,*_cy_r08 = 0x0,*_cy_r09 = 0x0,*_cy_r10 = 0x0,*_cy_r11 = 0x0,*_cy_r12 = 0x0,*_cy_r13 = 0x0,*_cy_r14 = 0x0,*_cy_r15 = 0x0,*_cy_r16 = 0x0,*_cy_r17 = 0x0,*_cy_r18 = 0x0,*_cy_r19 = 0x0,*_cy_r20 = 0x0,*_cy_r21 = 0x0,*_cy_r22 = 0x0,*_cy_r23 = 0x0,*_cy_r24 = 0x0,*_cy_r25 = 0x0,*_cy_r26 = 0x0,*_cy_r27 = 0x0,
	*_cy_i00 = 0x0,*_cy_i01 = 0x0,*_cy_i02 = 0x0,*_cy_i03 = 0x0,*_cy_i04 = 0x0,*_cy_i05 = 0x0,*_cy_i06 = 0x0,*_cy_i07 = 0x0,*_cy_i08 = 0x0,*_cy_i09 = 0x0,*_cy_i10 = 0x0,*_cy_i11 = 0x0,*_cy_i12 = 0x0,*_cy_i13 = 0x0,*_cy_i14 = 0x0,*_cy_i15 = 0x0,*_cy_i16 = 0x0,*_cy_i17 = 0x0,*_cy_i18 = 0x0,*_cy_i19 = 0x0,*_cy_i20 = 0x0,*_cy_i21 = 0x0,*_cy_i22 = 0x0,*_cy_i23 = 0x0,*_cy_i24 = 0x0,*_cy_i25 = 0x0,*_cy_i26 = 0x0,*_cy_i27 = 0x0;

  #ifdef USE_AVX512
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

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=ii0=ii1=ii2=ii3=ii4=ii5=ii6=-1;

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/RADIX;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/%d in %s.\n", iter,RADIX,func);
		if(INTERACT)fprintf(stderr,"%s",cbuf);
		fp = mlucas_fopen(   OFILE,"a");
		fq = mlucas_fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;
		err=ERR_CARRY;
		return(err);
	}

	if(p != psave || n != nsave) {	/* Exponent or array length change triggers re-init */
		first_entry=TRUE;
		/* To-do: Support #thread change here! */
	}

/*...initialize things upon first entry: */

	if(first_entry)
	{
		ASSERT(HERE, LO_ADD,"LO_ADD");
		psave = p;	nsave = n;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)RADIX));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;	/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

	#ifdef LOACC

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
		ASSERT(HERE,fabs(wts_mult[0]*inv_mult[0] - 1.0) < EPS, "wts_mults fail accuracy check!");
		//curr have w, 2/w, separate-mul-by-1-or-0.5 gives [w,w/2] and [1/w,2/w] for i = 0,1, resp:
		wts_mult[1] = 0.5*wts_mult[0];
		inv_mult[1] = 2.0*inv_mult[0];
		ASSERT(HERE,fabs(wts_mult[1]*inv_mult[1] - 1.0) < EPS, "wts_mults fail accuracy check!");

	#endif

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
			if(NDIVR    %CY_THREADS != 0) { WARN(HERE, "NDIVR    %CY_THREADS != 0", "", 1); return(ERR_ASSERT); }
			if(n_div_nwt%CY_THREADS != 0) { WARN(HERE, "n_div_nwt%CY_THREADS != 0", "", 1); return(ERR_ASSERT); }
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
					ASSERT(HERE, 0x0 != (tpool = threadpool_init(pool_work_units, MAX_THREADS, pool_work_units, &thread_control)), "threadpool_init failed!");
					printf("radix%d_ditN_cy_dif1: Init threadpool of %d threads\n", RADIX, pool_work_units);
				} else {
					main_work_units = 1;
					printf("radix%d_ditN_cy_dif1: CY_THREADS = 1: Using main execution thread, no threadpool needed.\n", RADIX);
				}
	
			#else
	
				pool_work_units = CY_THREADS;
				ASSERT(HERE, 0x0 != (tpool = threadpool_init(CY_THREADS, MAX_THREADS, CY_THREADS, &thread_control)), "threadpool_init failed!");
	
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
			tdat[ithread].arrdat = a;			/* Main data array */
			tdat[ithread].wt0 = wt0;
			tdat[ithread].wt1 = wt1;
		#ifdef LOACC
			tdat[ithread].wts_mult = wts_mult;
			tdat[ithread].inv_mult = inv_mult;
		#endif
			tdat[ithread].si  = si;
			tdat[ithread].rn0 = rn0;
			tdat[ithread].rn1 = rn1;
		}
	#endif

	#ifdef USE_SSE2

		ASSERT(HERE, ((long)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((long)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use vector-double type size (16 bytes for SSE2, 32 for AVX) to alloc a block of local storage
		// consisting of 128*2 vec_dbl and (8+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix28_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((long)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix28_creals_in_local_store);
		ASSERT(HERE, ((long)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

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
		cy_r00	= tmp + 0x0;
		cy_r04	= tmp + 0x1;
		cy_r08	= tmp + 0x2;
		cy_r12	= tmp + 0x3;
		cy_r16	= tmp + 0x4;
		cy_r20	= tmp + 0x5;
		cy_r24	= tmp + 0x6;
		cy_i00	= tmp + 0x7;
		cy_i04	= tmp + 0x8;
		cy_i08	= tmp + 0x9;
		cy_i12	= tmp + 0xa;
		cy_i16	= tmp + 0xb;
		cy_i20	= tmp + 0xc;
		cy_i24	= tmp + 0xd;
		tmp += 0xe;
	  #else
								tm2 = tmp + 0xe;
		cy_r00	= tmp + 0x0;	cy_i00	= tm2 + 0x0;
		cy_r02	= tmp + 0x1;	cy_i02	= tm2 + 0x1;
		cy_r04	= tmp + 0x2;	cy_i04	= tm2 + 0x2;
		cy_r06	= tmp + 0x3;	cy_i06	= tm2 + 0x3;
		cy_r08	= tmp + 0x4;	cy_i08	= tm2 + 0x4;
		cy_r10	= tmp + 0x5;	cy_i10	= tm2 + 0x5;
		cy_r12	= tmp + 0x6;	cy_i12	= tm2 + 0x6;
		cy_r14	= tmp + 0x7;	cy_i14	= tm2 + 0x7;
		cy_r16	= tmp + 0x8;	cy_i16	= tm2 + 0x8;
		cy_r18	= tmp + 0x9;	cy_i18	= tm2 + 0x9;
		cy_r20	= tmp + 0xa;	cy_i20	= tm2 + 0xa;
		cy_r22	= tmp + 0xb;	cy_i22	= tm2 + 0xb;
		cy_r24	= tmp + 0xc;	cy_i24	= tm2 + 0xc;
		cy_r26	= tmp + 0xd;	cy_i26	= tm2 + 0xd;
		tmp += 0x1c;
	  #endif
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 2 = 0x58 [avx] or 0x66 [sse2]; This is where the value of half_arr_offset28 comes from
		half_arr= tmp + 0x02;	/* This table needs 20x16 bytes */

		/* These remain fixed: */
		VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one, 1.0  );	// 1.0,2.0 needed for FMA support
	  #ifdef USE_AVX2
		// AVX2 (i.e. FMA)means non-Nussbaumer radix-7, uses these sincos constants:
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
		nbytes = (long)cy_r00 - (long)two;	// #bytes in above block of data, allowing for 'holes' and assuming only that cc0 and sse2_rnd bookend the block
		tmp = two;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}
		nbytes = sz_vd;	// sse2_rnd is a solo (in the SIMD-vector) datum
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
		The pattern of the negacyclic-DWT-weights ("nDWTs") applied to the 28 complex outputs of the final-DIT-pass is like so:
		The nDWTs multiplying each set of 28 DIT DFT outputs are simply the product of a single complex-root "base multiplier" rbase
		(separately computed for each batch of DFT outputs), which "base root" multiplies the first 28 [4*28]th roots of unity, i.e.

			 rbase * (j*I*Pi/2)/28, for j = 0, ..., 27 .

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
		nbytes = RADIX*sz_vd/2;	// 7 AVX-register-sized complex data

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
		nbytes = 4*sz_vd;	// 2 AVX-register-sized complex data

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
	  #ifdef LOACC
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
		nbytes = 96 << l2_sz_vd;
	  #else
		nbytes = 64 << l2_sz_vd;
	  #endif

	#else

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
	  #ifdef LOACC
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
		nbytes = 24 << l2_sz_vd;
	  #else
		nbytes = 16 << l2_sz_vd;
	  #endif

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

	nbytes = 4 << l2_sz_vd;

#ifdef USE_AVX
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

	#ifdef USE_AVX
		bjmodn00 = (int*)(sinwtm1 + RE_IM_STRIDE);
	#else
		bjmodn00 = (int*)(sse_n   + RE_IM_STRIDE);
	#endif
		bjmodn01 = bjmodn00 +  1;
		bjmodn02 = bjmodn00 +  2;
		bjmodn03 = bjmodn00 +  3;
		bjmodn04 = bjmodn00 +  4;
		bjmodn05 = bjmodn00 +  5;
		bjmodn06 = bjmodn00 +  6;
		bjmodn07 = bjmodn00 +  7;
		bjmodn08 = bjmodn00 +  8;
		bjmodn09 = bjmodn00 +  9;
		bjmodn10 = bjmodn00 + 10;
		bjmodn11 = bjmodn00 + 11;
		bjmodn12 = bjmodn00 + 12;
		bjmodn13 = bjmodn00 + 13;
		bjmodn14 = bjmodn00 + 14;
		bjmodn15 = bjmodn00 + 15;
		bjmodn16 = bjmodn00 + 16;
		bjmodn17 = bjmodn00 + 17;
		bjmodn18 = bjmodn00 + 18;
		bjmodn19 = bjmodn00 + 19;
		bjmodn20 = bjmodn00 + 20;
		bjmodn21 = bjmodn00 + 21;
		bjmodn22 = bjmodn00 + 22;
		bjmodn23 = bjmodn00 + 23;
		bjmodn24 = bjmodn00 + 24;
		bjmodn25 = bjmodn00 + 25;
		bjmodn26 = bjmodn00 + 26;
		bjmodn27 = bjmodn00 + 27;

	/*********** Defer the per-thread local-mem-block copy until after added wts-index precomputation below ************/
	#endif	/* USE_SSE2 */

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
		#ifdef USE_SSE2
			tdat[ithread].s1p00r = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (long)tdat[ithread].s1p00r + ((long)half_arr - (long)s1p00r);
		#else
			// In scalar mode use these 2 ptrs to pass wts_idx_incr and the base/baseinv/etc array-ptrs:
			tdat[ithread].s1p00r   = (double *)foo_array;
			tdat[ithread].half_arr = (double *)&wts_idx_incr;
		#endif	// USE_SSE2
		}
	#endif

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

		poff[0] =   0; poff[1] = p04; poff[2] = p08; poff[3] = p12; poff[4] = p16; poff[5] = p20; poff[6] = p24;

		ASSERT(HERE, p01+p01 == p02, "p01+p01 != p02");
		ASSERT(HERE, p02+p02 == p04, "p02+p02 != p04");
		ASSERT(HERE, p04+p04 == p08, "p04+p04 != p08");
		ASSERT(HERE, p08+p04 == p12, "p08+p04 != p12");
		ASSERT(HERE, p12+p04 == p16, "p12+p04 != p16");
		ASSERT(HERE, p16+p04 == p20, "p16+p04 != p20");
		ASSERT(HERE, p20+p04 == p24, "p20+p04 != p24");

		if(_cy_r00)	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;

			free((void *)_bjmodn00); _bjmodn00 = 0x0;
			free((void *)_bjmodn01); _bjmodn01 = 0x0;
			free((void *)_bjmodn02); _bjmodn02 = 0x0;
			free((void *)_bjmodn03); _bjmodn03 = 0x0;
			free((void *)_bjmodn04); _bjmodn04 = 0x0;
			free((void *)_bjmodn05); _bjmodn05 = 0x0;
			free((void *)_bjmodn06); _bjmodn06 = 0x0;
			free((void *)_bjmodn07); _bjmodn07 = 0x0;
			free((void *)_bjmodn08); _bjmodn08 = 0x0;
			free((void *)_bjmodn09); _bjmodn09 = 0x0;
			free((void *)_bjmodn10); _bjmodn10 = 0x0;
			free((void *)_bjmodn11); _bjmodn11 = 0x0;
			free((void *)_bjmodn12); _bjmodn12 = 0x0;
			free((void *)_bjmodn13); _bjmodn13 = 0x0;
			free((void *)_bjmodn14); _bjmodn14 = 0x0;
			free((void *)_bjmodn15); _bjmodn15 = 0x0;
			free((void *)_bjmodn16); _bjmodn16 = 0x0;
			free((void *)_bjmodn17); _bjmodn17 = 0x0;
			free((void *)_bjmodn18); _bjmodn18 = 0x0;
			free((void *)_bjmodn19); _bjmodn19 = 0x0;
			free((void *)_bjmodn20); _bjmodn20 = 0x0;
			free((void *)_bjmodn21); _bjmodn21 = 0x0;
			free((void *)_bjmodn22); _bjmodn22 = 0x0;
			free((void *)_bjmodn23); _bjmodn23 = 0x0;
			free((void *)_bjmodn24); _bjmodn24 = 0x0;
			free((void *)_bjmodn25); _bjmodn25 = 0x0;
			free((void *)_bjmodn26); _bjmodn26 = 0x0;
			free((void *)_bjmodn27); _bjmodn27 = 0x0;

			free((void *)_cy_r00); _cy_r00 = 0x0;		free((void *)_cy_i00); _cy_i00 = 0x0;
			free((void *)_cy_r01); _cy_r01 = 0x0;		free((void *)_cy_i01); _cy_i01 = 0x0;
			free((void *)_cy_r02); _cy_r02 = 0x0;		free((void *)_cy_i02); _cy_i02 = 0x0;
			free((void *)_cy_r03); _cy_r03 = 0x0;		free((void *)_cy_i03); _cy_i03 = 0x0;
			free((void *)_cy_r04); _cy_r04 = 0x0;		free((void *)_cy_i04); _cy_i04 = 0x0;
			free((void *)_cy_r05); _cy_r05 = 0x0;		free((void *)_cy_i05); _cy_i05 = 0x0;
			free((void *)_cy_r06); _cy_r06 = 0x0;		free((void *)_cy_i06); _cy_i06 = 0x0;
			free((void *)_cy_r07); _cy_r07 = 0x0;		free((void *)_cy_i07); _cy_i07 = 0x0;
			free((void *)_cy_r08); _cy_r08 = 0x0;		free((void *)_cy_i08); _cy_i08 = 0x0;
			free((void *)_cy_r09); _cy_r09 = 0x0;		free((void *)_cy_i09); _cy_i09 = 0x0;
			free((void *)_cy_r10); _cy_r10 = 0x0;		free((void *)_cy_i10); _cy_i10 = 0x0;
			free((void *)_cy_r11); _cy_r11 = 0x0;		free((void *)_cy_i11); _cy_i11 = 0x0;
			free((void *)_cy_r12); _cy_r12 = 0x0;		free((void *)_cy_i12); _cy_i12 = 0x0;
			free((void *)_cy_r13); _cy_r13 = 0x0;		free((void *)_cy_i13); _cy_i13 = 0x0;
			free((void *)_cy_r14); _cy_r14 = 0x0;		free((void *)_cy_i14); _cy_i14 = 0x0;
			free((void *)_cy_r15); _cy_r15 = 0x0;		free((void *)_cy_i15); _cy_i15 = 0x0;
			free((void *)_cy_r16); _cy_r16 = 0x0;		free((void *)_cy_i16); _cy_i16 = 0x0;
			free((void *)_cy_r17); _cy_r17 = 0x0;		free((void *)_cy_i17); _cy_i17 = 0x0;
			free((void *)_cy_r18); _cy_r18 = 0x0;		free((void *)_cy_i18); _cy_i18 = 0x0;
			free((void *)_cy_r19); _cy_r19 = 0x0;		free((void *)_cy_i19); _cy_i19 = 0x0;
			free((void *)_cy_r20); _cy_r20 = 0x0;		free((void *)_cy_i20); _cy_i20 = 0x0;
			free((void *)_cy_r21); _cy_r21 = 0x0;		free((void *)_cy_i21); _cy_i21 = 0x0;
			free((void *)_cy_r22); _cy_r22 = 0x0;		free((void *)_cy_i22); _cy_i22 = 0x0;
			free((void *)_cy_r23); _cy_r23 = 0x0;		free((void *)_cy_i23); _cy_i23 = 0x0;
			free((void *)_cy_r24); _cy_r24 = 0x0;		free((void *)_cy_i24); _cy_i24 = 0x0;
			free((void *)_cy_r25); _cy_r25 = 0x0;		free((void *)_cy_i25); _cy_i25 = 0x0;
			free((void *)_cy_r26); _cy_r26 = 0x0;		free((void *)_cy_i26); _cy_i26 = 0x0;
			free((void *)_cy_r27); _cy_r27 = 0x0;		free((void *)_cy_i27); _cy_i27 = 0x0;

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
		_bjmodn00	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn00== 0x0);
		_bjmodn01	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn01== 0x0);
		_bjmodn02	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn02== 0x0);
		_bjmodn03	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn03== 0x0);
		_bjmodn04	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn04== 0x0);
		_bjmodn05	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn05== 0x0);
		_bjmodn06	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn06== 0x0);
		_bjmodn07	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn07== 0x0);
		_bjmodn08	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn08== 0x0);
		_bjmodn09	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn09== 0x0);
		_bjmodn10	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn10== 0x0);
		_bjmodn11	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn11== 0x0);
		_bjmodn12	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn12== 0x0);
		_bjmodn13	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn13== 0x0);
		_bjmodn14	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn14== 0x0);
		_bjmodn15	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn15== 0x0);
		_bjmodn16	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn16== 0x0);
		_bjmodn17	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn17== 0x0);
		_bjmodn18	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn18== 0x0);
		_bjmodn19	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn19== 0x0);
		_bjmodn20	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn20== 0x0);
		_bjmodn21	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn21== 0x0);
		_bjmodn22	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn22== 0x0);
		_bjmodn23	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn23== 0x0);
		_bjmodn24	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn24== 0x0);
		_bjmodn25	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn25== 0x0);
		_bjmodn26	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn26== 0x0);
		_bjmodn27	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn27== 0x0);
		_jstart  	= (int *)malloc(j);	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(j);	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(j);	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(j);	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(j);	ptr_prod += (uint32)(_co3     == 0x0);

		j = CY_THREADS*sizeof(double);
		_cy_r00	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r00== 0x0);
		_cy_r01	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r01== 0x0);
		_cy_r02	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r02== 0x0);
		_cy_r03	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r03== 0x0);
		_cy_r04	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r04== 0x0);
		_cy_r05	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r05== 0x0);
		_cy_r06	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r06== 0x0);
		_cy_r07	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r07== 0x0);
		_cy_r08	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r08== 0x0);
		_cy_r09	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r09== 0x0);
		_cy_r10	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r10== 0x0);
		_cy_r11	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r11== 0x0);
		_cy_r12	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r12== 0x0);
		_cy_r13	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r13== 0x0);
		_cy_r14	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r14== 0x0);
		_cy_r15	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r15== 0x0);
		_cy_r16	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r16== 0x0);
		_cy_r17	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r17== 0x0);
		_cy_r18	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r18== 0x0);
		_cy_r19	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r19== 0x0);
		_cy_r20	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r20== 0x0);
		_cy_r21	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r21== 0x0);
		_cy_r22	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r22== 0x0);
		_cy_r23	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r23== 0x0);
		_cy_r24	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r24== 0x0);
		_cy_r25	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r25== 0x0);
		_cy_r26	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r26== 0x0);
		_cy_r27	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r27== 0x0);

		_cy_i00	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i00== 0x0);
		_cy_i01	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i01== 0x0);
		_cy_i02	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i02== 0x0);
		_cy_i03	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i03== 0x0);
		_cy_i04	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i04== 0x0);
		_cy_i05	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i05== 0x0);
		_cy_i06	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i06== 0x0);
		_cy_i07	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i07== 0x0);
		_cy_i08	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i08== 0x0);
		_cy_i09	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i09== 0x0);
		_cy_i10	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i10== 0x0);
		_cy_i11	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i11== 0x0);
		_cy_i12	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i12== 0x0);
		_cy_i13	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i13== 0x0);
		_cy_i14	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i14== 0x0);
		_cy_i15	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i15== 0x0);
		_cy_i16	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i16== 0x0);
		_cy_i17	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i17== 0x0);
		_cy_i18	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i18== 0x0);
		_cy_i19	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i19== 0x0);
		_cy_i20	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i20== 0x0);
		_cy_i21	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i21== 0x0);
		_cy_i22	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i22== 0x0);
		_cy_i23	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i23== 0x0);
		_cy_i24	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i24== 0x0);
		_cy_i25	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i25== 0x0);
		_cy_i26	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i26== 0x0);
		_cy_i27	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i27== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays in radix28_ditN_cy_dif1.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/radix-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
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
		ASSERT(HERE, _bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");

		first_entry=FALSE;
	}	/* endif(first_entry) */

/*...The radix-28 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_cy_r00[ithread] = 0;	_cy_i00[ithread] = 0;
		_cy_r01[ithread] = 0;	_cy_i01[ithread] = 0;
		_cy_r02[ithread] = 0;	_cy_i02[ithread] = 0;
		_cy_r03[ithread] = 0;	_cy_i03[ithread] = 0;
		_cy_r04[ithread] = 0;	_cy_i04[ithread] = 0;
		_cy_r05[ithread] = 0;	_cy_i05[ithread] = 0;
		_cy_r06[ithread] = 0;	_cy_i06[ithread] = 0;
		_cy_r07[ithread] = 0;	_cy_i07[ithread] = 0;
		_cy_r08[ithread] = 0;	_cy_i08[ithread] = 0;
		_cy_r09[ithread] = 0;	_cy_i09[ithread] = 0;
		_cy_r10[ithread] = 0;	_cy_i10[ithread] = 0;
		_cy_r11[ithread] = 0;	_cy_i11[ithread] = 0;
		_cy_r12[ithread] = 0;	_cy_i12[ithread] = 0;
		_cy_r13[ithread] = 0;	_cy_i13[ithread] = 0;
		_cy_r14[ithread] = 0;	_cy_i14[ithread] = 0;
		_cy_r15[ithread] = 0;	_cy_i15[ithread] = 0;
		_cy_r16[ithread] = 0;	_cy_i16[ithread] = 0;
		_cy_r17[ithread] = 0;	_cy_i17[ithread] = 0;
		_cy_r18[ithread] = 0;	_cy_i18[ithread] = 0;
		_cy_r19[ithread] = 0;	_cy_i19[ithread] = 0;
		_cy_r20[ithread] = 0;	_cy_i20[ithread] = 0;
		_cy_r21[ithread] = 0;	_cy_i21[ithread] = 0;
		_cy_r22[ithread] = 0;	_cy_i22[ithread] = 0;
		_cy_r23[ithread] = 0;	_cy_i23[ithread] = 0;
		_cy_r24[ithread] = 0;	_cy_i24[ithread] = 0;
		_cy_r25[ithread] = 0;	_cy_i25[ithread] = 0;
		_cy_r26[ithread] = 0;	_cy_i26[ithread] = 0;
		_cy_r27[ithread] = 0;	_cy_i27[ithread] = 0;
	}
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy_r00[      0] = -2;
	}

	*fracmax=0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

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

		/* For Fermat-mod, IBDWT access patterns repeat with period NWT = {odd part of radix0},
		so for even radix0 values still only need [radix0 >> trailz(radix0)] bjmodn and ii's:
		*/
		/* indices into IBDWT weights arrays (mod NWT) is here: */
		ii0= 0;
		ii1= (SW_DIV_N*(NDIVR >> 1)) % nwt;	// nwt *not* a power of 2, must use library-mod!
		MOD_ADD32(ii1,ii1,nwt,ii2);
		MOD_ADD32(ii2,ii1,nwt,ii3);
		MOD_ADD32(ii3,ii1,nwt,ii4);
		MOD_ADD32(ii4,ii1,nwt,ii5);
		MOD_ADD32(ii5,ii1,nwt,ii6);
	}

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
//	printf("_bjmodnini[CY_THREADS] = %u\n",_bjmodnini[CY_THREADS]);
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_bjmodn00[ithread] = _bjmodnini[ithread];
		MOD_ADD32(_bjmodn00[ithread], j, n, _bjmodn01[ithread]);
		MOD_ADD32(_bjmodn01[ithread], j, n, _bjmodn02[ithread]);
		MOD_ADD32(_bjmodn02[ithread], j, n, _bjmodn03[ithread]);
		MOD_ADD32(_bjmodn03[ithread], j, n, _bjmodn04[ithread]);
		MOD_ADD32(_bjmodn04[ithread], j, n, _bjmodn05[ithread]);
		MOD_ADD32(_bjmodn05[ithread], j, n, _bjmodn06[ithread]);
		MOD_ADD32(_bjmodn06[ithread], j, n, _bjmodn07[ithread]);
		MOD_ADD32(_bjmodn07[ithread], j, n, _bjmodn08[ithread]);
		MOD_ADD32(_bjmodn08[ithread], j, n, _bjmodn09[ithread]);
		MOD_ADD32(_bjmodn09[ithread], j, n, _bjmodn10[ithread]);
		MOD_ADD32(_bjmodn10[ithread], j, n, _bjmodn11[ithread]);
		MOD_ADD32(_bjmodn11[ithread], j, n, _bjmodn12[ithread]);
		MOD_ADD32(_bjmodn12[ithread], j, n, _bjmodn13[ithread]);
		MOD_ADD32(_bjmodn13[ithread], j, n, _bjmodn14[ithread]);
		MOD_ADD32(_bjmodn14[ithread], j, n, _bjmodn15[ithread]);
		MOD_ADD32(_bjmodn15[ithread], j, n, _bjmodn16[ithread]);
		MOD_ADD32(_bjmodn16[ithread], j, n, _bjmodn17[ithread]);
		MOD_ADD32(_bjmodn17[ithread], j, n, _bjmodn18[ithread]);
		MOD_ADD32(_bjmodn18[ithread], j, n, _bjmodn19[ithread]);
		MOD_ADD32(_bjmodn19[ithread], j, n, _bjmodn20[ithread]);
		MOD_ADD32(_bjmodn20[ithread], j, n, _bjmodn21[ithread]);
		MOD_ADD32(_bjmodn21[ithread], j, n, _bjmodn22[ithread]);
		MOD_ADD32(_bjmodn22[ithread], j, n, _bjmodn23[ithread]);
		MOD_ADD32(_bjmodn23[ithread], j, n, _bjmodn24[ithread]);
		MOD_ADD32(_bjmodn24[ithread], j, n, _bjmodn25[ithread]);
		MOD_ADD32(_bjmodn25[ithread], j, n, _bjmodn26[ithread]);
		MOD_ADD32(_bjmodn26[ithread], j, n, _bjmodn27[ithread]);

		// Every (ODD_RADIX)th bjmodn initializer needs to be forced-to-bigword in fermat-mod DWT case:
		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		{
			/* Start this value off at N in Fermat-mod case, so (bjmodn >= sw) check in
			fermat_carry_norm_errcheck (cf. carry.h) yields a bigword (i == 1) for j= 0:
			*/
			_bjmodn00[ithread] = n;
			_bjmodn07[ithread] = n;
			_bjmodn14[ithread] = n;
			_bjmodn21[ithread] = n;
		}
	}

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		/* Find the circular-index-shift described in the head-of-file comments of radix56_ditN_cy_dif1.c, by serching bjmodn01 ... bjmodn[nwt] for the one == bw: */
		if( _bjmodn01[0] == bw ) { wts_idx_incr = 1; };
		if( _bjmodn02[0] == bw ) { wts_idx_incr = 2; };
		if( _bjmodn03[0] == bw ) { wts_idx_incr = 3; };
		if( _bjmodn04[0] == bw ) { wts_idx_incr = 4; };
		if( _bjmodn05[0] == bw ) { wts_idx_incr = 5; };
		if( _bjmodn06[0] == bw ) { wts_idx_incr = 6; };

	#ifdef USE_SSE2
		wts_idx_inc2 = wts_idx_incr << (2*l2_sz_vd - 3);	/* In the SIMD version, use icycle0-6 as actual address
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

		/* Need this both in scalar mode and to ease the SSE2-array init */
		j = _bjmodn00[0] > sw;		bs_arr[0] = base[j];	bsinv_arr[0] = baseinv[j];
		j = _bjmodn01[0] > sw;		bs_arr[1] = base[j];	bsinv_arr[1] = baseinv[j];
		j = _bjmodn02[0] > sw;		bs_arr[2] = base[j];	bsinv_arr[2] = baseinv[j];
		j = _bjmodn03[0] > sw;		bs_arr[3] = base[j];	bsinv_arr[3] = baseinv[j];
		j = _bjmodn04[0] > sw;		bs_arr[4] = base[j];	bsinv_arr[4] = baseinv[j];
		j = _bjmodn05[0] > sw;		bs_arr[5] = base[j];	bsinv_arr[5] = baseinv[j];
		j = _bjmodn06[0] > sw;		bs_arr[6] = base[j];	bsinv_arr[6] = baseinv[j];

		/* Give icycle indices their proper starting values: */
		icycle0 = 0x0;
		icycle1 = 0x1;
		icycle2 = 0x2;
		icycle3 = 0x3;
		icycle4 = 0x4;
		icycle5 = 0x5;
		icycle6 = 0x6;

		/* Need this both in scalar mode and to ease the SSE2-array init */
		wt_arr[0] = wt0[ii0];	wtinv_arr[0] = scale*wt1[ii0];
		wt_arr[1] = wt0[ii1];	wtinv_arr[1] = scale*wt1[ii1];
		wt_arr[2] = wt0[ii2];	wtinv_arr[2] = scale*wt1[ii2];
		wt_arr[3] = wt0[ii3];	wtinv_arr[3] = scale*wt1[ii3];
		wt_arr[4] = wt0[ii4];	wtinv_arr[4] = scale*wt1[ii4];
		wt_arr[5] = wt0[ii5];	wtinv_arr[5] = scale*wt1[ii5];
		wt_arr[6] = wt0[ii6];	wtinv_arr[6] = scale*wt1[ii6];

	#ifdef USE_SSE2

		tmp = half_arr;
		tmp->d0 = wt_arr[icycle0];	++tmp;
		tmp->d0 = wt_arr[icycle1];	++tmp;
		tmp->d0 = wt_arr[icycle2];	++tmp;
		tmp->d0 = wt_arr[icycle3];	++tmp;
		tmp->d0 = wt_arr[icycle4];	++tmp;
		tmp->d0 = wt_arr[icycle5];	++tmp;
		tmp->d0 = wt_arr[icycle6];	++tmp;

		tmp->d0 = wtinv_arr[icycle0];	++tmp;
		tmp->d0 = wtinv_arr[icycle1];	++tmp;
		tmp->d0 = wtinv_arr[icycle2];	++tmp;
		tmp->d0 = wtinv_arr[icycle3];	++tmp;
		tmp->d0 = wtinv_arr[icycle4];	++tmp;
		tmp->d0 = wtinv_arr[icycle5];	++tmp;
		tmp->d0 = wtinv_arr[icycle6];	++tmp;

		/* Now set the imaginary parts to the values corresponding to the 2nd of each pair of scalar-mode loop passes.
		Use this sequence for mod-add, as it is faster than general-mod '% nwt'. The reason we do not use the MOD_ADD32
		macor here is that wts_idx_incr is precomputed constant, so we have also pre-subtracted the modulus nwt from it:
		*/
		jcycle0 = icycle0 + wts_idx_incr;		jcycle0 += ( (-(jcycle0 < 0)) & nwt);
		jcycle1 = icycle1 + wts_idx_incr;		jcycle1 += ( (-(jcycle1 < 0)) & nwt);
		jcycle2 = icycle2 + wts_idx_incr;		jcycle2 += ( (-(jcycle2 < 0)) & nwt);
		jcycle3 = icycle3 + wts_idx_incr;		jcycle3 += ( (-(jcycle3 < 0)) & nwt);
		jcycle4 = icycle4 + wts_idx_incr;		jcycle4 += ( (-(jcycle4 < 0)) & nwt);
		jcycle5 = icycle5 + wts_idx_incr;		jcycle5 += ( (-(jcycle5 < 0)) & nwt);
		jcycle6 = icycle6 + wts_idx_incr;		jcycle6 += ( (-(jcycle6 < 0)) & nwt);

		tmp = half_arr;
		tmp->d1 = wt_arr[jcycle0];	++tmp;
		tmp->d1 = wt_arr[jcycle1];	++tmp;
		tmp->d1 = wt_arr[jcycle2];	++tmp;
		tmp->d1 = wt_arr[jcycle3];	++tmp;
		tmp->d1 = wt_arr[jcycle4];	++tmp;
		tmp->d1 = wt_arr[jcycle5];	++tmp;
		tmp->d1 = wt_arr[jcycle6];	++tmp;

		tmp->d1 = wtinv_arr[jcycle0];	++tmp;
		tmp->d1 = wtinv_arr[jcycle1];	++tmp;
		tmp->d1 = wtinv_arr[jcycle2];	++tmp;
		tmp->d1 = wtinv_arr[jcycle3];	++tmp;
		tmp->d1 = wtinv_arr[jcycle4];	++tmp;
		tmp->d1 = wtinv_arr[jcycle5];	++tmp;
		tmp->d1 = wtinv_arr[jcycle6];	++tmp;

	  #ifdef USE_AVX
		kcycle0 = jcycle0 + wts_idx_incr;		kcycle0 += ( (-(kcycle0 < 0)) & nwt);
		kcycle1 = jcycle1 + wts_idx_incr;		kcycle1 += ( (-(kcycle1 < 0)) & nwt);
		kcycle2 = jcycle2 + wts_idx_incr;		kcycle2 += ( (-(kcycle2 < 0)) & nwt);
		kcycle3 = jcycle3 + wts_idx_incr;		kcycle3 += ( (-(kcycle3 < 0)) & nwt);
		kcycle4 = jcycle4 + wts_idx_incr;		kcycle4 += ( (-(kcycle4 < 0)) & nwt);
		kcycle5 = jcycle5 + wts_idx_incr;		kcycle5 += ( (-(kcycle5 < 0)) & nwt);
		kcycle6 = jcycle6 + wts_idx_incr;		kcycle6 += ( (-(kcycle6 < 0)) & nwt);

		tmp = half_arr;
		tmp->d2 = wt_arr[kcycle0];	++tmp;
		tmp->d2 = wt_arr[kcycle1];	++tmp;
		tmp->d2 = wt_arr[kcycle2];	++tmp;
		tmp->d2 = wt_arr[kcycle3];	++tmp;
		tmp->d2 = wt_arr[kcycle4];	++tmp;
		tmp->d2 = wt_arr[kcycle5];	++tmp;
		tmp->d2 = wt_arr[kcycle6];	++tmp;

		tmp->d2 = wtinv_arr[kcycle0];	++tmp;
		tmp->d2 = wtinv_arr[kcycle1];	++tmp;
		tmp->d2 = wtinv_arr[kcycle2];	++tmp;
		tmp->d2 = wtinv_arr[kcycle3];	++tmp;
		tmp->d2 = wtinv_arr[kcycle4];	++tmp;
		tmp->d2 = wtinv_arr[kcycle5];	++tmp;
		tmp->d2 = wtinv_arr[kcycle6];	++tmp;

		lcycle0 = kcycle0 + wts_idx_incr;		lcycle0 += ( (-(lcycle0 < 0)) & nwt);
		lcycle1 = kcycle1 + wts_idx_incr;		lcycle1 += ( (-(lcycle1 < 0)) & nwt);
		lcycle2 = kcycle2 + wts_idx_incr;		lcycle2 += ( (-(lcycle2 < 0)) & nwt);
		lcycle3 = kcycle3 + wts_idx_incr;		lcycle3 += ( (-(lcycle3 < 0)) & nwt);
		lcycle4 = kcycle4 + wts_idx_incr;		lcycle4 += ( (-(lcycle4 < 0)) & nwt);
		lcycle5 = kcycle5 + wts_idx_incr;		lcycle5 += ( (-(lcycle5 < 0)) & nwt);
		lcycle6 = kcycle6 + wts_idx_incr;		lcycle6 += ( (-(lcycle6 < 0)) & nwt);

		tmp = half_arr;
		tmp->d3 = wt_arr[lcycle0];	++tmp;
		tmp->d3 = wt_arr[lcycle1];	++tmp;
		tmp->d3 = wt_arr[lcycle2];	++tmp;
		tmp->d3 = wt_arr[lcycle3];	++tmp;
		tmp->d3 = wt_arr[lcycle4];	++tmp;
		tmp->d3 = wt_arr[lcycle5];	++tmp;
		tmp->d3 = wt_arr[lcycle6];	++tmp;

		tmp->d3 = wtinv_arr[lcycle0];	++tmp;
		tmp->d3 = wtinv_arr[lcycle1];	++tmp;
		tmp->d3 = wtinv_arr[lcycle2];	++tmp;
		tmp->d3 = wtinv_arr[lcycle3];	++tmp;
		tmp->d3 = wtinv_arr[lcycle4];	++tmp;
		tmp->d3 = wtinv_arr[lcycle5];	++tmp;
		tmp->d3 = wtinv_arr[lcycle6];	++tmp;
	  #endif

		tmp = half_arr + ODD_RADIX*2;	/* Put the base-mini-arrays right after the weights */

	  #ifdef USE_AVX
		// Each transposed-data quartet in the AVX carry macro needs linearly incrementing bs_arr data (mod ODD_RADIX);
		// Need all [ODD_RADIX] possible such length-4 index subsequences, which will be accessed via their head element
		// by the [ijkl]cycle* index quartets in the respective carry-macro call:
		tmp->d0 = bs_arr[0];	tmp->d1 = bs_arr[1];	tmp->d2 = bs_arr[2];	tmp->d3 = bs_arr[3];	++tmp;
		tmp->d0 = bs_arr[1];	tmp->d1 = bs_arr[2];	tmp->d2 = bs_arr[3];	tmp->d3 = bs_arr[4];	++tmp;
		tmp->d0 = bs_arr[2];	tmp->d1 = bs_arr[3];	tmp->d2 = bs_arr[4];	tmp->d3 = bs_arr[5];	++tmp;
		tmp->d0 = bs_arr[3];	tmp->d1 = bs_arr[4];	tmp->d2 = bs_arr[5];	tmp->d3 = bs_arr[6];	++tmp;
		tmp->d0 = bs_arr[4];	tmp->d1 = bs_arr[5];	tmp->d2 = bs_arr[6];	tmp->d3 = bs_arr[0];	++tmp;
		tmp->d0 = bs_arr[5];	tmp->d1 = bs_arr[6];	tmp->d2 = bs_arr[0];	tmp->d3 = bs_arr[1];	++tmp;
		tmp->d0 = bs_arr[6];	tmp->d1 = bs_arr[0];	tmp->d2 = bs_arr[1];	tmp->d3 = bs_arr[2];	++tmp;

		tmp->d0 = bsinv_arr[0];	tmp->d1 = bsinv_arr[1];	tmp->d2 = bsinv_arr[2];	tmp->d3 = bsinv_arr[3];	++tmp;
		tmp->d0 = bsinv_arr[1];	tmp->d1 = bsinv_arr[2];	tmp->d2 = bsinv_arr[3];	tmp->d3 = bsinv_arr[4];	++tmp;
		tmp->d0 = bsinv_arr[2];	tmp->d1 = bsinv_arr[3];	tmp->d2 = bsinv_arr[4];	tmp->d3 = bsinv_arr[5];	++tmp;
		tmp->d0 = bsinv_arr[3];	tmp->d1 = bsinv_arr[4];	tmp->d2 = bsinv_arr[5];	tmp->d3 = bsinv_arr[6];	++tmp;
		tmp->d0 = bsinv_arr[4];	tmp->d1 = bsinv_arr[5];	tmp->d2 = bsinv_arr[6];	tmp->d3 = bsinv_arr[0];	++tmp;
		tmp->d0 = bsinv_arr[5];	tmp->d1 = bsinv_arr[6];	tmp->d2 = bsinv_arr[0];	tmp->d3 = bsinv_arr[1];	++tmp;
		tmp->d0 = bsinv_arr[6];	tmp->d1 = bsinv_arr[0];	tmp->d2 = bsinv_arr[1];	tmp->d3 = bsinv_arr[2];	++tmp;

	  #else

		/* In SSE2 mode, because we apply doubled weights to data arranged as [a.re,b.re,...],[a.im,b.im,...] but apply
		doubled base multipliers to shuffled data [a.re,a.im],[b.re,b.im],... (i.e. shuffled to yield same data layout as
		in the scalar case), the weights need to have disparate real and imag parts, but the base/baseinv terms do not: */
		VEC_DBL_INIT(tmp, bs_arr[0]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[1]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[2]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[3]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[4]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[5]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[6]);	++tmp;

		VEC_DBL_INIT(tmp, bsinv_arr[0]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[1]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[2]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[3]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[4]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[5]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[6]);	++tmp;

	  #endif

		// Propagate the above consts to the remaining threads:
		nbytes = RADIX*sz_vd;
		tmp = half_arr;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

		icycle0 <<= l2_sz_vd;		jcycle0 <<= l2_sz_vd;
		icycle1 <<= l2_sz_vd;		jcycle1 <<= l2_sz_vd;
		icycle2 <<= l2_sz_vd;		jcycle2 <<= l2_sz_vd;
		icycle3 <<= l2_sz_vd;		jcycle3 <<= l2_sz_vd;
		icycle4 <<= l2_sz_vd;		jcycle4 <<= l2_sz_vd;
		icycle5 <<= l2_sz_vd;		jcycle5 <<= l2_sz_vd;
		icycle6 <<= l2_sz_vd;		jcycle6 <<= l2_sz_vd;
	  #ifdef USE_AVX
		kcycle0 <<= l2_sz_vd;		lcycle0 <<= l2_sz_vd;
		kcycle1 <<= l2_sz_vd;		lcycle1 <<= l2_sz_vd;
		kcycle2 <<= l2_sz_vd;		lcycle2 <<= l2_sz_vd;
		kcycle3 <<= l2_sz_vd;		lcycle3 <<= l2_sz_vd;
		kcycle4 <<= l2_sz_vd;		lcycle4 <<= l2_sz_vd;
		kcycle5 <<= l2_sz_vd;		lcycle5 <<= l2_sz_vd;
		kcycle6 <<= l2_sz_vd;		lcycle6 <<= l2_sz_vd;
	  #endif

	#endif	/* USE_SSE2 */
	}	// if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)

/*	fprintf(stderr, "radix28_ditN_cy_dif1: wts_idx_incr = %d\n", wts_idx_incr);*/

  #ifdef USE_PTHREAD

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart = _jstart[ithread];
		jhi    = _jhi[ithread];
	}
	// These inits must occur just once, in full-pass mode,
	// in order to get the array-index-offset values of the icycle/jcycle indices right:
	else if(full_pass)	// Fermat-mod
	{
		tdat[0].icycle0 = icycle0;
		tdat[0].icycle1 = icycle1;
		tdat[0].icycle2 = icycle2;
		tdat[0].icycle3 = icycle3;
		tdat[0].icycle4 = icycle4;
		tdat[0].icycle5 = icycle5;
		tdat[0].icycle6 = icycle6;

	#ifdef USE_SSE2
		tdat[0].wts_idx_inc2 = wts_idx_inc2;
		tdat[0].jcycle0 = jcycle0;
		tdat[0].jcycle1 = jcycle1;
		tdat[0].jcycle2 = jcycle2;
		tdat[0].jcycle3 = jcycle3;
		tdat[0].jcycle4 = jcycle4;
		tdat[0].jcycle5 = jcycle5;
		tdat[0].jcycle6 = jcycle6;

	  #ifdef USE_AVX
		tdat[0].kcycle0 = kcycle0;
		tdat[0].kcycle1 = kcycle1;
		tdat[0].kcycle2 = kcycle2;
		tdat[0].kcycle3 = kcycle3;
		tdat[0].kcycle4 = kcycle4;
		tdat[0].kcycle5 = kcycle5;
		tdat[0].kcycle6 = kcycle6;

		tdat[0].lcycle0 = lcycle0;
		tdat[0].lcycle1 = lcycle1;
		tdat[0].lcycle2 = lcycle2;
		tdat[0].lcycle3 = lcycle3;
		tdat[0].lcycle4 = lcycle4;
		tdat[0].lcycle5 = lcycle5;
		tdat[0].lcycle6 = lcycle6;
	  #endif
	#endif

		// For remaining threads, simulate the loop-evolution of the above indices:
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			jstart = _jstart[ithread];
			jhi    = _jhi[ithread];
			for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
			{
				for(j = jstart; j < jhi; j += stride)
				{
				#ifndef USE_SSE2	// Scalar-double mode uses non-pointerized icycle values:

					icycle0 += wts_idx_incr;		icycle0 += ( (-(int)((uint32)icycle0 >> 31)) & nwt);
					icycle1 += wts_idx_incr;		icycle1 += ( (-(int)((uint32)icycle1 >> 31)) & nwt);
					icycle2 += wts_idx_incr;		icycle2 += ( (-(int)((uint32)icycle2 >> 31)) & nwt);
					icycle3 += wts_idx_incr;		icycle3 += ( (-(int)((uint32)icycle3 >> 31)) & nwt);
					icycle4 += wts_idx_incr;		icycle4 += ( (-(int)((uint32)icycle4 >> 31)) & nwt);
					icycle5 += wts_idx_incr;		icycle5 += ( (-(int)((uint32)icycle5 >> 31)) & nwt);
					icycle6 += wts_idx_incr;		icycle6 += ( (-(int)((uint32)icycle6 >> 31)) & nwt);

				#else

					icycle0 += wts_idx_inc2;		icycle0 += ( (-(int)((uint32)icycle0 >> 31)) & nwt16);
					icycle1 += wts_idx_inc2;		icycle1 += ( (-(int)((uint32)icycle1 >> 31)) & nwt16);
					icycle2 += wts_idx_inc2;		icycle2 += ( (-(int)((uint32)icycle2 >> 31)) & nwt16);
					icycle3 += wts_idx_inc2;		icycle3 += ( (-(int)((uint32)icycle3 >> 31)) & nwt16);
					icycle4 += wts_idx_inc2;		icycle4 += ( (-(int)((uint32)icycle4 >> 31)) & nwt16);
					icycle5 += wts_idx_inc2;		icycle5 += ( (-(int)((uint32)icycle5 >> 31)) & nwt16);
					icycle6 += wts_idx_inc2;		icycle6 += ( (-(int)((uint32)icycle6 >> 31)) & nwt16);

					jcycle0 += wts_idx_inc2;		jcycle0 += ( (-(int)((uint32)jcycle0 >> 31)) & nwt16);
					jcycle1 += wts_idx_inc2;		jcycle1 += ( (-(int)((uint32)jcycle1 >> 31)) & nwt16);
					jcycle2 += wts_idx_inc2;		jcycle2 += ( (-(int)((uint32)jcycle2 >> 31)) & nwt16);
					jcycle3 += wts_idx_inc2;		jcycle3 += ( (-(int)((uint32)jcycle3 >> 31)) & nwt16);
					jcycle4 += wts_idx_inc2;		jcycle4 += ( (-(int)((uint32)jcycle4 >> 31)) & nwt16);
					jcycle5 += wts_idx_inc2;		jcycle5 += ( (-(int)((uint32)jcycle5 >> 31)) & nwt16);
					jcycle6 += wts_idx_inc2;		jcycle6 += ( (-(int)((uint32)jcycle6 >> 31)) & nwt16);
				  #ifdef USE_AVX
					kcycle0 += wts_idx_inc2;		kcycle0 += ( (-(int)((uint32)kcycle0 >> 31)) & nwt16);
					kcycle1 += wts_idx_inc2;		kcycle1 += ( (-(int)((uint32)kcycle1 >> 31)) & nwt16);
					kcycle2 += wts_idx_inc2;		kcycle2 += ( (-(int)((uint32)kcycle2 >> 31)) & nwt16);
					kcycle3 += wts_idx_inc2;		kcycle3 += ( (-(int)((uint32)kcycle3 >> 31)) & nwt16);
					kcycle4 += wts_idx_inc2;		kcycle4 += ( (-(int)((uint32)kcycle4 >> 31)) & nwt16);
					kcycle5 += wts_idx_inc2;		kcycle5 += ( (-(int)((uint32)kcycle5 >> 31)) & nwt16);
					kcycle6 += wts_idx_inc2;		kcycle6 += ( (-(int)((uint32)kcycle6 >> 31)) & nwt16);

					lcycle0 += wts_idx_inc2;		lcycle0 += ( (-(int)((uint32)lcycle0 >> 31)) & nwt16);
					lcycle1 += wts_idx_inc2;		lcycle1 += ( (-(int)((uint32)lcycle1 >> 31)) & nwt16);
					lcycle2 += wts_idx_inc2;		lcycle2 += ( (-(int)((uint32)lcycle2 >> 31)) & nwt16);
					lcycle3 += wts_idx_inc2;		lcycle3 += ( (-(int)((uint32)lcycle3 >> 31)) & nwt16);
					lcycle4 += wts_idx_inc2;		lcycle4 += ( (-(int)((uint32)lcycle4 >> 31)) & nwt16);
					lcycle5 += wts_idx_inc2;		lcycle5 += ( (-(int)((uint32)lcycle5 >> 31)) & nwt16);
					lcycle6 += wts_idx_inc2;		lcycle6 += ( (-(int)((uint32)lcycle6 >> 31)) & nwt16);
				  #endif
				#endif
				}
			}
			tdat[ithread].icycle0 = icycle0;
			tdat[ithread].icycle1 = icycle1;
			tdat[ithread].icycle2 = icycle2;
			tdat[ithread].icycle3 = icycle3;
			tdat[ithread].icycle4 = icycle4;
			tdat[ithread].icycle5 = icycle5;
			tdat[ithread].icycle6 = icycle6;

		#ifdef USE_SSE2
			tdat[ithread].wts_idx_inc2 = wts_idx_inc2;
			tdat[ithread].jcycle0 = jcycle0;
			tdat[ithread].jcycle1 = jcycle1;
			tdat[ithread].jcycle2 = jcycle2;
			tdat[ithread].jcycle3 = jcycle3;
			tdat[ithread].jcycle4 = jcycle4;
			tdat[ithread].jcycle5 = jcycle5;
			tdat[ithread].jcycle6 = jcycle6;
		  #ifdef USE_AVX
			tdat[ithread].kcycle0 = kcycle0;
			tdat[ithread].kcycle1 = kcycle1;
			tdat[ithread].kcycle2 = kcycle2;
			tdat[ithread].kcycle3 = kcycle3;
			tdat[ithread].kcycle4 = kcycle4;
			tdat[ithread].kcycle5 = kcycle5;
			tdat[ithread].kcycle6 = kcycle6;

			tdat[ithread].lcycle0 = lcycle0;
			tdat[ithread].lcycle1 = lcycle1;
			tdat[ithread].lcycle2 = lcycle2;
			tdat[ithread].lcycle3 = lcycle3;
			tdat[ithread].lcycle4 = lcycle4;
			tdat[ithread].lcycle5 = lcycle5;
			tdat[ithread].lcycle6 = lcycle6;
		  #endif
		#endif
		}
	}

  #endif

#ifdef USE_SSE2

	tmp = max_err;	VEC_DBL_INIT(tmp, 0.0);
	tm2 = tmp + cslots_in_local_store;
	for(ithread = 1; ithread < CY_THREADS; ++ithread) {
		memcpy(tm2, tmp, sz_vd);
		tmp = tm2;		tm2 += cslots_in_local_store;
	}

#endif	// USE_PTHREAD

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
	// int data:
		ASSERT(HERE, tdat[ithread].tid == ithread, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].ndivr == NDIVR, "thread-local memcheck fail!");

		tdat[ithread].khi    = khi;
		tdat[ithread].i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		tdat[ithread].jstart = _jstart[ithread];
		tdat[ithread].jhi    = _jhi[ithread];

		tdat[ithread].col = _col[ithread];
		tdat[ithread].co2 = _co2[ithread];
		tdat[ithread].co3 = _co3[ithread];
		ASSERT(HERE, tdat[ithread].sw  == sw, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].nwt == nwt, "thread-local memcheck fail!");

	// double data:
		tdat[ithread].maxerr = 0.0;
		tdat[ithread].scale = scale;

	// pointer data:
		ASSERT(HERE, tdat[ithread].arrdat == a, "thread-local memcheck fail!");			/* Main data array */
		ASSERT(HERE, tdat[ithread].wt0 == wt0, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].wt1 == wt1, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].si  == si, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].rn0 == rn0, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].rn1 == rn1, "thread-local memcheck fail!");
	#ifdef USE_SSE2
		ASSERT(HERE, tdat[ithread].wts_idx_inc2 == wts_idx_inc2, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].s1p00r == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].s1p00r;
		ASSERT(HERE, ((tmp + 0x38)->d0 == 2.0 && (tmp + 0x38)->d1 == 2.0), "thread-local memcheck failed!");
		ASSERT(HERE, ((tmp + half_arr_offset28-1)->d0 == crnd && (tmp + half_arr_offset28-1)->d1 == crnd), "thread-local memcheck failed!");
	#endif
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX
			// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
			dtmp = (tmp + half_arr_offset28+40)->d0 * (tmp + half_arr_offset28+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp + half_arr_offset28+40)->d1 * (tmp + half_arr_offset28+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#elif defined(USE_SSE2)
			dtmp = (tmp + half_arr_offset28+10)->d0 * (tmp + half_arr_offset28+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp + half_arr_offset28+10)->d1 * (tmp + half_arr_offset28+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#endif
			tdat[ithread].bjmodn00 = _bjmodn00[ithread];
			tdat[ithread].bjmodn01 = _bjmodn01[ithread];
			tdat[ithread].bjmodn02 = _bjmodn02[ithread];
			tdat[ithread].bjmodn03 = _bjmodn03[ithread];
			tdat[ithread].bjmodn04 = _bjmodn04[ithread];
			tdat[ithread].bjmodn05 = _bjmodn05[ithread];
			tdat[ithread].bjmodn06 = _bjmodn06[ithread];
			tdat[ithread].bjmodn07 = _bjmodn07[ithread];
			tdat[ithread].bjmodn08 = _bjmodn08[ithread];
			tdat[ithread].bjmodn09 = _bjmodn09[ithread];
			tdat[ithread].bjmodn10 = _bjmodn10[ithread];
			tdat[ithread].bjmodn11 = _bjmodn11[ithread];
			tdat[ithread].bjmodn12 = _bjmodn12[ithread];
			tdat[ithread].bjmodn13 = _bjmodn13[ithread];
			tdat[ithread].bjmodn14 = _bjmodn14[ithread];
			tdat[ithread].bjmodn15 = _bjmodn15[ithread];
			tdat[ithread].bjmodn16 = _bjmodn16[ithread];
			tdat[ithread].bjmodn17 = _bjmodn17[ithread];
			tdat[ithread].bjmodn18 = _bjmodn18[ithread];
			tdat[ithread].bjmodn19 = _bjmodn19[ithread];
			tdat[ithread].bjmodn20 = _bjmodn20[ithread];
			tdat[ithread].bjmodn21 = _bjmodn21[ithread];
			tdat[ithread].bjmodn22 = _bjmodn22[ithread];
			tdat[ithread].bjmodn23 = _bjmodn23[ithread];
			tdat[ithread].bjmodn24 = _bjmodn24[ithread];
			tdat[ithread].bjmodn25 = _bjmodn25[ithread];
			tdat[ithread].bjmodn26 = _bjmodn26[ithread];
			tdat[ithread].bjmodn27 = _bjmodn27[ithread];
			/* init carries	*/
			tdat[ithread].cy_r00 = _cy_r00[ithread];
			tdat[ithread].cy_r01 = _cy_r01[ithread];
			tdat[ithread].cy_r02 = _cy_r02[ithread];
			tdat[ithread].cy_r03 = _cy_r03[ithread];
			tdat[ithread].cy_r04 = _cy_r04[ithread];
			tdat[ithread].cy_r05 = _cy_r05[ithread];
			tdat[ithread].cy_r06 = _cy_r06[ithread];
			tdat[ithread].cy_r07 = _cy_r07[ithread];
			tdat[ithread].cy_r08 = _cy_r08[ithread];
			tdat[ithread].cy_r09 = _cy_r09[ithread];
			tdat[ithread].cy_r10 = _cy_r10[ithread];
			tdat[ithread].cy_r11 = _cy_r11[ithread];
			tdat[ithread].cy_r12 = _cy_r12[ithread];
			tdat[ithread].cy_r13 = _cy_r13[ithread];
			tdat[ithread].cy_r14 = _cy_r14[ithread];
			tdat[ithread].cy_r15 = _cy_r15[ithread];
			tdat[ithread].cy_r16 = _cy_r16[ithread];
			tdat[ithread].cy_r17 = _cy_r17[ithread];
			tdat[ithread].cy_r18 = _cy_r18[ithread];
			tdat[ithread].cy_r19 = _cy_r19[ithread];
			tdat[ithread].cy_r20 = _cy_r20[ithread];
			tdat[ithread].cy_r21 = _cy_r21[ithread];
			tdat[ithread].cy_r22 = _cy_r22[ithread];
			tdat[ithread].cy_r23 = _cy_r23[ithread];
			tdat[ithread].cy_r24 = _cy_r24[ithread];
			tdat[ithread].cy_r25 = _cy_r25[ithread];
			tdat[ithread].cy_r26 = _cy_r26[ithread];
			tdat[ithread].cy_r27 = _cy_r27[ithread];
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
		#ifdef USE_SSE2
			dtmp = (tmp + half_arr_offset28)->d0 * (tmp + half_arr_offset28+ODD_RADIX)->d0;	ASSERT(HERE, fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp + half_arr_offset28)->d1 * (tmp + half_arr_offset28+ODD_RADIX)->d1;	ASSERT(HERE, fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
		#endif
			/* init carries	*/
			tdat[ithread].cy_r00 = _cy_r00[ithread];	tdat[ithread].cy_i00 = _cy_i00[ithread];
			tdat[ithread].cy_r01 = _cy_r01[ithread];	tdat[ithread].cy_i01 = _cy_i01[ithread];
			tdat[ithread].cy_r02 = _cy_r02[ithread];	tdat[ithread].cy_i02 = _cy_i02[ithread];
			tdat[ithread].cy_r03 = _cy_r03[ithread];	tdat[ithread].cy_i03 = _cy_i03[ithread];
			tdat[ithread].cy_r04 = _cy_r04[ithread];	tdat[ithread].cy_i04 = _cy_i04[ithread];
			tdat[ithread].cy_r05 = _cy_r05[ithread];	tdat[ithread].cy_i05 = _cy_i05[ithread];
			tdat[ithread].cy_r06 = _cy_r06[ithread];	tdat[ithread].cy_i06 = _cy_i06[ithread];
			tdat[ithread].cy_r07 = _cy_r07[ithread];	tdat[ithread].cy_i07 = _cy_i07[ithread];
			tdat[ithread].cy_r08 = _cy_r08[ithread];	tdat[ithread].cy_i08 = _cy_i08[ithread];
			tdat[ithread].cy_r09 = _cy_r09[ithread];	tdat[ithread].cy_i09 = _cy_i09[ithread];
			tdat[ithread].cy_r10 = _cy_r10[ithread];	tdat[ithread].cy_i10 = _cy_i10[ithread];
			tdat[ithread].cy_r11 = _cy_r11[ithread];	tdat[ithread].cy_i11 = _cy_i11[ithread];
			tdat[ithread].cy_r12 = _cy_r12[ithread];	tdat[ithread].cy_i12 = _cy_i12[ithread];
			tdat[ithread].cy_r13 = _cy_r13[ithread];	tdat[ithread].cy_i13 = _cy_i13[ithread];
			tdat[ithread].cy_r14 = _cy_r14[ithread];	tdat[ithread].cy_i14 = _cy_i14[ithread];
			tdat[ithread].cy_r15 = _cy_r15[ithread];	tdat[ithread].cy_i15 = _cy_i15[ithread];
			tdat[ithread].cy_r16 = _cy_r16[ithread];	tdat[ithread].cy_i16 = _cy_i16[ithread];
			tdat[ithread].cy_r17 = _cy_r17[ithread];	tdat[ithread].cy_i17 = _cy_i17[ithread];
			tdat[ithread].cy_r18 = _cy_r18[ithread];	tdat[ithread].cy_i18 = _cy_i18[ithread];
			tdat[ithread].cy_r19 = _cy_r19[ithread];	tdat[ithread].cy_i19 = _cy_i19[ithread];
			tdat[ithread].cy_r20 = _cy_r20[ithread];	tdat[ithread].cy_i20 = _cy_i20[ithread];
			tdat[ithread].cy_r21 = _cy_r21[ithread];	tdat[ithread].cy_i21 = _cy_i21[ithread];
			tdat[ithread].cy_r22 = _cy_r22[ithread];	tdat[ithread].cy_i22 = _cy_i22[ithread];
			tdat[ithread].cy_r23 = _cy_r23[ithread];	tdat[ithread].cy_i23 = _cy_i23[ithread];
			tdat[ithread].cy_r24 = _cy_r24[ithread];	tdat[ithread].cy_i24 = _cy_i24[ithread];
			tdat[ithread].cy_r25 = _cy_r25[ithread];	tdat[ithread].cy_i25 = _cy_i25[ithread];
			tdat[ithread].cy_r26 = _cy_r26[ithread];	tdat[ithread].cy_i26 = _cy_i26[ithread];
			tdat[ithread].cy_r27 = _cy_r27[ithread];	tdat[ithread].cy_i27 = _cy_i27[ithread];
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

	#ifdef USE_SSE2
		*bjmodn00 = _bjmodn00[ithread];
		*bjmodn01 = _bjmodn01[ithread];
		*bjmodn02 = _bjmodn02[ithread];
		*bjmodn03 = _bjmodn03[ithread];
		*bjmodn04 = _bjmodn04[ithread];
		*bjmodn05 = _bjmodn05[ithread];
		*bjmodn06 = _bjmodn06[ithread];
		*bjmodn07 = _bjmodn07[ithread];
		*bjmodn08 = _bjmodn08[ithread];
		*bjmodn09 = _bjmodn09[ithread];
		*bjmodn10 = _bjmodn10[ithread];
		*bjmodn11 = _bjmodn11[ithread];
		*bjmodn12 = _bjmodn12[ithread];
		*bjmodn13 = _bjmodn13[ithread];
		*bjmodn14 = _bjmodn14[ithread];
		*bjmodn15 = _bjmodn15[ithread];
		*bjmodn16 = _bjmodn16[ithread];
		*bjmodn17 = _bjmodn17[ithread];
		*bjmodn18 = _bjmodn18[ithread];
		*bjmodn19 = _bjmodn19[ithread];
		*bjmodn20 = _bjmodn20[ithread];
		*bjmodn21 = _bjmodn21[ithread];
		*bjmodn22 = _bjmodn22[ithread];
		*bjmodn23 = _bjmodn23[ithread];
		*bjmodn24 = _bjmodn24[ithread];
		*bjmodn25 = _bjmodn25[ithread];
		*bjmodn26 = _bjmodn26[ithread];
		*bjmodn27 = _bjmodn27[ithread];
	#else
		bjmodn00 = _bjmodn00[ithread];
		bjmodn01 = _bjmodn01[ithread];
		bjmodn02 = _bjmodn02[ithread];
		bjmodn03 = _bjmodn03[ithread];
		bjmodn04 = _bjmodn04[ithread];
		bjmodn05 = _bjmodn05[ithread];
		bjmodn06 = _bjmodn06[ithread];
		bjmodn07 = _bjmodn07[ithread];
		bjmodn08 = _bjmodn08[ithread];
		bjmodn09 = _bjmodn09[ithread];
		bjmodn10 = _bjmodn10[ithread];
		bjmodn11 = _bjmodn11[ithread];
		bjmodn12 = _bjmodn12[ithread];
		bjmodn13 = _bjmodn13[ithread];
		bjmodn14 = _bjmodn14[ithread];
		bjmodn15 = _bjmodn15[ithread];
		bjmodn16 = _bjmodn16[ithread];
		bjmodn17 = _bjmodn17[ithread];
		bjmodn18 = _bjmodn18[ithread];
		bjmodn19 = _bjmodn19[ithread];
		bjmodn20 = _bjmodn20[ithread];
		bjmodn21 = _bjmodn21[ithread];
		bjmodn22 = _bjmodn22[ithread];
		bjmodn23 = _bjmodn23[ithread];
		bjmodn24 = _bjmodn24[ithread];
		bjmodn25 = _bjmodn25[ithread];
		bjmodn26 = _bjmodn26[ithread];
		bjmodn27 = _bjmodn27[ithread];
	#endif

		i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		jstart = _jstart[ithread];
		jhi    = _jhi[ithread];

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			col = _col[ithread];
			co2 = _co2[ithread];
			co3 = _co3[ithread];

			/* init carries	*/
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
			cy_r00->d0 = _cy_r00[ithread];	cy_r00->d1 = _cy_r01[ithread];	cy_r00->d2 = _cy_r02[ithread];	cy_r00->d3 = _cy_r03[ithread];
			cy_r04->d0 = _cy_r04[ithread];	cy_r04->d1 = _cy_r05[ithread];	cy_r04->d2 = _cy_r06[ithread];	cy_r04->d3 = _cy_r07[ithread];
			cy_r08->d0 = _cy_r08[ithread];	cy_r08->d1 = _cy_r09[ithread];	cy_r08->d2 = _cy_r10[ithread];	cy_r08->d3 = _cy_r11[ithread];
			cy_r12->d0 = _cy_r12[ithread];	cy_r12->d1 = _cy_r13[ithread];	cy_r12->d2 = _cy_r14[ithread];	cy_r12->d3 = _cy_r15[ithread];
			cy_r16->d0 = _cy_r16[ithread];	cy_r16->d1 = _cy_r17[ithread];	cy_r16->d2 = _cy_r18[ithread];	cy_r16->d3 = _cy_r19[ithread];
			cy_r20->d0 = _cy_r20[ithread];	cy_r20->d1 = _cy_r21[ithread];	cy_r20->d2 = _cy_r22[ithread];	cy_r20->d3 = _cy_r23[ithread];
			cy_r24->d0 = _cy_r24[ithread];	cy_r24->d1 = _cy_r25[ithread];	cy_r24->d2 = _cy_r26[ithread];	cy_r24->d3 = _cy_r27[ithread];
		#elif defined(USE_SSE2)
			cy_r00->d0 = _cy_r00[ithread];	cy_r00->d1 = _cy_r01[ithread];
			cy_r02->d0 = _cy_r02[ithread];	cy_r02->d1 = _cy_r03[ithread];
			cy_r04->d0 = _cy_r04[ithread];	cy_r04->d1 = _cy_r05[ithread];
			cy_r06->d0 = _cy_r06[ithread];	cy_r06->d1 = _cy_r07[ithread];
			cy_r08->d0 = _cy_r08[ithread];	cy_r08->d1 = _cy_r09[ithread];
			cy_r10->d0 = _cy_r10[ithread];	cy_r10->d1 = _cy_r11[ithread];
			cy_r12->d0 = _cy_r12[ithread];	cy_r12->d1 = _cy_r13[ithread];
			cy_r14->d0 = _cy_r14[ithread];	cy_r14->d1 = _cy_r15[ithread];
			cy_r16->d0 = _cy_r16[ithread];	cy_r16->d1 = _cy_r17[ithread];
			cy_r18->d0 = _cy_r18[ithread];	cy_r18->d1 = _cy_r19[ithread];
			cy_r20->d0 = _cy_r20[ithread];	cy_r20->d1 = _cy_r21[ithread];
			cy_r22->d0 = _cy_r22[ithread];	cy_r22->d1 = _cy_r23[ithread];
			cy_r24->d0 = _cy_r24[ithread];	cy_r24->d1 = _cy_r25[ithread];
			cy_r26->d0 = _cy_r26[ithread];	cy_r26->d1 = _cy_r27[ithread];
		#else
			cy_r00 = _cy_r00[ithread];
			cy_r01 = _cy_r01[ithread];
			cy_r02 = _cy_r02[ithread];
			cy_r03 = _cy_r03[ithread];
			cy_r04 = _cy_r04[ithread];
			cy_r05 = _cy_r05[ithread];
			cy_r06 = _cy_r06[ithread];
			cy_r07 = _cy_r07[ithread];
			cy_r08 = _cy_r08[ithread];
			cy_r09 = _cy_r09[ithread];
			cy_r10 = _cy_r10[ithread];
			cy_r11 = _cy_r11[ithread];
			cy_r12 = _cy_r12[ithread];
			cy_r13 = _cy_r13[ithread];
			cy_r14 = _cy_r14[ithread];
			cy_r15 = _cy_r15[ithread];
			cy_r16 = _cy_r16[ithread];
			cy_r17 = _cy_r17[ithread];
			cy_r18 = _cy_r18[ithread];
			cy_r19 = _cy_r19[ithread];
			cy_r20 = _cy_r20[ithread];
			cy_r21 = _cy_r21[ithread];
			cy_r22 = _cy_r22[ithread];
			cy_r23 = _cy_r23[ithread];
			cy_r24 = _cy_r24[ithread];
			cy_r25 = _cy_r25[ithread];
			cy_r26 = _cy_r26[ithread];
			cy_r27 = _cy_r27[ithread];
		#endif
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
		#ifdef USE_AVX
			cy_r00->d0 = _cy_r00[ithread];	cy_r00->d1 = _cy_r01[ithread];	cy_r00->d2 = _cy_r02[ithread];	cy_r00->d3 = _cy_r03[ithread];
			cy_r04->d0 = _cy_r04[ithread];	cy_r04->d1 = _cy_r05[ithread];	cy_r04->d2 = _cy_r06[ithread];	cy_r04->d3 = _cy_r07[ithread];
			cy_r08->d0 = _cy_r08[ithread];	cy_r08->d1 = _cy_r09[ithread];	cy_r08->d2 = _cy_r10[ithread];	cy_r08->d3 = _cy_r11[ithread];
			cy_r12->d0 = _cy_r12[ithread];	cy_r12->d1 = _cy_r13[ithread];	cy_r12->d2 = _cy_r14[ithread];	cy_r12->d3 = _cy_r15[ithread];
			cy_r16->d0 = _cy_r16[ithread];	cy_r16->d1 = _cy_r17[ithread];	cy_r16->d2 = _cy_r18[ithread];	cy_r16->d3 = _cy_r19[ithread];
			cy_r20->d0 = _cy_r20[ithread];	cy_r20->d1 = _cy_r21[ithread];	cy_r20->d2 = _cy_r22[ithread];	cy_r20->d3 = _cy_r23[ithread];
			cy_r24->d0 = _cy_r24[ithread];	cy_r24->d1 = _cy_r25[ithread];	cy_r24->d2 = _cy_r26[ithread];	cy_r24->d3 = _cy_r27[ithread];

			cy_i00->d0 = _cy_i00[ithread];	cy_i00->d1 = _cy_i01[ithread];	cy_i00->d2 = _cy_i02[ithread];	cy_i00->d3 = _cy_i03[ithread];
			cy_i04->d0 = _cy_i04[ithread];	cy_i04->d1 = _cy_i05[ithread];	cy_i04->d2 = _cy_i06[ithread];	cy_i04->d3 = _cy_i07[ithread];
			cy_i08->d0 = _cy_i08[ithread];	cy_i08->d1 = _cy_i09[ithread];	cy_i08->d2 = _cy_i10[ithread];	cy_i08->d3 = _cy_i11[ithread];
			cy_i12->d0 = _cy_i12[ithread];	cy_i12->d1 = _cy_i13[ithread];	cy_i12->d2 = _cy_i14[ithread];	cy_i12->d3 = _cy_i15[ithread];
			cy_i16->d0 = _cy_i16[ithread];	cy_i16->d1 = _cy_i17[ithread];	cy_i16->d2 = _cy_i18[ithread];	cy_i16->d3 = _cy_i19[ithread];
			cy_i20->d0 = _cy_i20[ithread];	cy_i20->d1 = _cy_i21[ithread];	cy_i20->d2 = _cy_i22[ithread];	cy_i20->d3 = _cy_i23[ithread];
			cy_i24->d0 = _cy_i24[ithread];	cy_i24->d1 = _cy_i25[ithread];	cy_i24->d2 = _cy_i26[ithread];	cy_i24->d3 = _cy_i27[ithread];
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			cy_r00->d0 = _cy_r00[ithread];	cy_r00->d1 = _cy_i00[ithread];
			cy_r02->d0 = _cy_r01[ithread];	cy_r02->d1 = _cy_i01[ithread];
			cy_r04->d0 = _cy_r02[ithread];	cy_r04->d1 = _cy_i02[ithread];
			cy_r06->d0 = _cy_r03[ithread];	cy_r06->d1 = _cy_i03[ithread];
			cy_r08->d0 = _cy_r04[ithread];	cy_r08->d1 = _cy_i04[ithread];
			cy_r10->d0 = _cy_r05[ithread];	cy_r10->d1 = _cy_i05[ithread];
			cy_r12->d0 = _cy_r06[ithread];	cy_r12->d1 = _cy_i06[ithread];
			cy_r14->d0 = _cy_r07[ithread];	cy_r14->d1 = _cy_i07[ithread];
			cy_r16->d0 = _cy_r08[ithread];	cy_r16->d1 = _cy_i08[ithread];
			cy_r18->d0 = _cy_r09[ithread];	cy_r18->d1 = _cy_i09[ithread];
			cy_r20->d0 = _cy_r10[ithread];	cy_r20->d1 = _cy_i10[ithread];
			cy_r22->d0 = _cy_r11[ithread];	cy_r22->d1 = _cy_i11[ithread];
			cy_r24->d0 = _cy_r12[ithread];	cy_r24->d1 = _cy_i12[ithread];
			cy_r26->d0 = _cy_r13[ithread];	cy_r26->d1 = _cy_i13[ithread];
			cy_i00->d0 = _cy_r14[ithread];	cy_i00->d1 = _cy_i14[ithread];
			cy_i02->d0 = _cy_r15[ithread];	cy_i02->d1 = _cy_i15[ithread];
			cy_i04->d0 = _cy_r16[ithread];	cy_i04->d1 = _cy_i16[ithread];
			cy_i06->d0 = _cy_r17[ithread];	cy_i06->d1 = _cy_i17[ithread];
			cy_i08->d0 = _cy_r18[ithread];	cy_i08->d1 = _cy_i18[ithread];
			cy_i10->d0 = _cy_r19[ithread];	cy_i10->d1 = _cy_i19[ithread];
			cy_i12->d0 = _cy_r20[ithread];	cy_i12->d1 = _cy_i20[ithread];
			cy_i14->d0 = _cy_r21[ithread];	cy_i14->d1 = _cy_i21[ithread];
			cy_i16->d0 = _cy_r22[ithread];	cy_i16->d1 = _cy_i22[ithread];
			cy_i18->d0 = _cy_r23[ithread];	cy_i18->d1 = _cy_i23[ithread];
			cy_i20->d0 = _cy_r24[ithread];	cy_i20->d1 = _cy_i24[ithread];
			cy_i22->d0 = _cy_r25[ithread];	cy_i22->d1 = _cy_i25[ithread];
			cy_i24->d0 = _cy_r26[ithread];	cy_i24->d1 = _cy_i26[ithread];
			cy_i26->d0 = _cy_r27[ithread];	cy_i26->d1 = _cy_i27[ithread];
		#else
			cy_r00 = _cy_r00[ithread];		cy_i00 = _cy_i00[ithread];
			cy_r01 = _cy_r01[ithread];		cy_i01 = _cy_i01[ithread];
			cy_r02 = _cy_r02[ithread];		cy_i02 = _cy_i02[ithread];
			cy_r03 = _cy_r03[ithread];		cy_i03 = _cy_i03[ithread];
			cy_r04 = _cy_r04[ithread];		cy_i04 = _cy_i04[ithread];
			cy_r05 = _cy_r05[ithread];		cy_i05 = _cy_i05[ithread];
			cy_r06 = _cy_r06[ithread];		cy_i06 = _cy_i06[ithread];
			cy_r07 = _cy_r07[ithread];		cy_i07 = _cy_i07[ithread];
			cy_r08 = _cy_r08[ithread];		cy_i08 = _cy_i08[ithread];
			cy_r09 = _cy_r09[ithread];		cy_i09 = _cy_i09[ithread];
			cy_r10 = _cy_r10[ithread];		cy_i10 = _cy_i10[ithread];
			cy_r11 = _cy_r11[ithread];		cy_i11 = _cy_i11[ithread];
			cy_r12 = _cy_r12[ithread];		cy_i12 = _cy_i12[ithread];
			cy_r13 = _cy_r13[ithread];		cy_i13 = _cy_i13[ithread];
			cy_r14 = _cy_r14[ithread];		cy_i14 = _cy_i14[ithread];
			cy_r15 = _cy_r15[ithread];		cy_i15 = _cy_i15[ithread];
			cy_r16 = _cy_r16[ithread];		cy_i16 = _cy_i16[ithread];
			cy_r17 = _cy_r17[ithread];		cy_i17 = _cy_i17[ithread];
			cy_r18 = _cy_r18[ithread];		cy_i18 = _cy_i18[ithread];
			cy_r19 = _cy_r19[ithread];		cy_i19 = _cy_i19[ithread];
			cy_r20 = _cy_r20[ithread];		cy_i20 = _cy_i20[ithread];
			cy_r21 = _cy_r21[ithread];		cy_i21 = _cy_i21[ithread];
			cy_r22 = _cy_r22[ithread];		cy_i22 = _cy_i22[ithread];
			cy_r23 = _cy_r23[ithread];		cy_i23 = _cy_i23[ithread];
			cy_r24 = _cy_r24[ithread];		cy_i24 = _cy_i24[ithread];
			cy_r25 = _cy_r25[ithread];		cy_i25 = _cy_i25[ithread];
			cy_r26 = _cy_r26[ithread];		cy_i26 = _cy_i26[ithread];
			cy_r27 = _cy_r27[ithread];		cy_i27 = _cy_i27[ithread];
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
		#ifdef USE_AVX
			_cy_r00[ithread] = cy_r00->d0;	_cy_r01[ithread] = cy_r00->d1;	_cy_r02[ithread] = cy_r00->d2;	_cy_r03[ithread] = cy_r00->d3;
			_cy_r04[ithread] = cy_r04->d0;	_cy_r05[ithread] = cy_r04->d1;	_cy_r06[ithread] = cy_r04->d2;	_cy_r07[ithread] = cy_r04->d3;
			_cy_r08[ithread] = cy_r08->d0;	_cy_r09[ithread] = cy_r08->d1;	_cy_r10[ithread] = cy_r08->d2;	_cy_r11[ithread] = cy_r08->d3;
			_cy_r12[ithread] = cy_r12->d0;	_cy_r13[ithread] = cy_r12->d1;	_cy_r14[ithread] = cy_r12->d2;	_cy_r15[ithread] = cy_r12->d3;
			_cy_r16[ithread] = cy_r16->d0;	_cy_r17[ithread] = cy_r16->d1;	_cy_r18[ithread] = cy_r16->d2;	_cy_r19[ithread] = cy_r16->d3;
			_cy_r20[ithread] = cy_r20->d0;	_cy_r21[ithread] = cy_r20->d1;	_cy_r22[ithread] = cy_r20->d2;	_cy_r23[ithread] = cy_r20->d3;
			_cy_r24[ithread] = cy_r24->d0;	_cy_r25[ithread] = cy_r24->d1;	_cy_r26[ithread] = cy_r24->d2;	_cy_r27[ithread] = cy_r24->d3;
			if(full_pass) maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			_cy_r00[ithread] = cy_r00->d0;	_cy_r01[ithread] = cy_r00->d1;
			_cy_r02[ithread] = cy_r02->d0;	_cy_r03[ithread] = cy_r02->d1;
			_cy_r04[ithread] = cy_r04->d0;	_cy_r05[ithread] = cy_r04->d1;
			_cy_r06[ithread] = cy_r06->d0;	_cy_r07[ithread] = cy_r06->d1;
			_cy_r08[ithread] = cy_r08->d0;	_cy_r09[ithread] = cy_r08->d1;
			_cy_r10[ithread] = cy_r10->d0;	_cy_r11[ithread] = cy_r10->d1;
			_cy_r12[ithread] = cy_r12->d0;	_cy_r13[ithread] = cy_r12->d1;
			_cy_r14[ithread] = cy_r14->d0;	_cy_r15[ithread] = cy_r14->d1;
			_cy_r16[ithread] = cy_r16->d0;	_cy_r17[ithread] = cy_r16->d1;
			_cy_r18[ithread] = cy_r18->d0;	_cy_r19[ithread] = cy_r18->d1;
			_cy_r20[ithread] = cy_r20->d0;	_cy_r21[ithread] = cy_r20->d1;
			_cy_r22[ithread] = cy_r22->d0;	_cy_r23[ithread] = cy_r22->d1;
			_cy_r24[ithread] = cy_r24->d0;	_cy_r25[ithread] = cy_r24->d1;
			_cy_r26[ithread] = cy_r26->d0;	_cy_r27[ithread] = cy_r26->d1;
			if(full_pass) maxerr = MAX(max_err->d0,max_err->d1);
		#else
			_cy_r00[ithread] = cy_r00;
			_cy_r01[ithread] = cy_r01;
			_cy_r02[ithread] = cy_r02;
			_cy_r03[ithread] = cy_r03;
			_cy_r04[ithread] = cy_r04;
			_cy_r05[ithread] = cy_r05;
			_cy_r06[ithread] = cy_r06;
			_cy_r07[ithread] = cy_r07;
			_cy_r08[ithread] = cy_r08;
			_cy_r09[ithread] = cy_r09;
			_cy_r10[ithread] = cy_r10;
			_cy_r11[ithread] = cy_r11;
			_cy_r12[ithread] = cy_r12;
			_cy_r13[ithread] = cy_r13;
			_cy_r14[ithread] = cy_r14;
			_cy_r15[ithread] = cy_r15;
			_cy_r16[ithread] = cy_r16;
			_cy_r17[ithread] = cy_r17;
			_cy_r18[ithread] = cy_r18;
			_cy_r19[ithread] = cy_r19;
			_cy_r20[ithread] = cy_r20;
			_cy_r21[ithread] = cy_r21;
			_cy_r22[ithread] = cy_r22;
			_cy_r23[ithread] = cy_r23;
			_cy_r24[ithread] = cy_r24;
			_cy_r25[ithread] = cy_r25;
			_cy_r26[ithread] = cy_r26;
			_cy_r27[ithread] = cy_r27;
		#endif
		}
		else
		{
		#ifdef USE_AVX
			_cy_r00[ithread] = cy_r00->d0;	_cy_r01[ithread] = cy_r00->d1;	_cy_r02[ithread] = cy_r00->d2;	_cy_r03[ithread] = cy_r00->d3;
			_cy_r04[ithread] = cy_r04->d0;	_cy_r05[ithread] = cy_r04->d1;	_cy_r06[ithread] = cy_r04->d2;	_cy_r07[ithread] = cy_r04->d3;
			_cy_r08[ithread] = cy_r08->d0;	_cy_r09[ithread] = cy_r08->d1;	_cy_r10[ithread] = cy_r08->d2;	_cy_r11[ithread] = cy_r08->d3;
			_cy_r12[ithread] = cy_r12->d0;	_cy_r13[ithread] = cy_r12->d1;	_cy_r14[ithread] = cy_r12->d2;	_cy_r15[ithread] = cy_r12->d3;
			_cy_r16[ithread] = cy_r16->d0;	_cy_r17[ithread] = cy_r16->d1;	_cy_r18[ithread] = cy_r16->d2;	_cy_r19[ithread] = cy_r16->d3;
			_cy_r20[ithread] = cy_r20->d0;	_cy_r21[ithread] = cy_r20->d1;	_cy_r22[ithread] = cy_r20->d2;	_cy_r23[ithread] = cy_r20->d3;
			_cy_r24[ithread] = cy_r24->d0;	_cy_r25[ithread] = cy_r24->d1;	_cy_r26[ithread] = cy_r24->d2;	_cy_r27[ithread] = cy_r24->d3;

			_cy_i00[ithread] = cy_i00->d0;	_cy_i01[ithread] = cy_i00->d1;	_cy_i02[ithread] = cy_i00->d2;	_cy_i03[ithread] = cy_i00->d3;
			_cy_i04[ithread] = cy_i04->d0;	_cy_i05[ithread] = cy_i04->d1;	_cy_i06[ithread] = cy_i04->d2;	_cy_i07[ithread] = cy_i04->d3;
			_cy_i08[ithread] = cy_i08->d0;	_cy_i09[ithread] = cy_i08->d1;	_cy_i10[ithread] = cy_i08->d2;	_cy_i11[ithread] = cy_i08->d3;
			_cy_i12[ithread] = cy_i12->d0;	_cy_i13[ithread] = cy_i12->d1;	_cy_i14[ithread] = cy_i12->d2;	_cy_i15[ithread] = cy_i12->d3;
			_cy_i16[ithread] = cy_i16->d0;	_cy_i17[ithread] = cy_i16->d1;	_cy_i18[ithread] = cy_i16->d2;	_cy_i19[ithread] = cy_i16->d3;
			_cy_i20[ithread] = cy_i20->d0;	_cy_i21[ithread] = cy_i20->d1;	_cy_i22[ithread] = cy_i20->d2;	_cy_i23[ithread] = cy_i20->d3;
			_cy_i24[ithread] = cy_i24->d0;	_cy_i25[ithread] = cy_i24->d1;	_cy_i26[ithread] = cy_i24->d2;	_cy_i27[ithread] = cy_i24->d3;
			if(full_pass) maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			_cy_r00[ithread] = cy_r00->d0;	_cy_i00[ithread] = cy_r00->d1;
			_cy_r01[ithread] = cy_r02->d0;	_cy_i01[ithread] = cy_r02->d1;
			_cy_r02[ithread] = cy_r04->d0;	_cy_i02[ithread] = cy_r04->d1;
			_cy_r03[ithread] = cy_r06->d0;	_cy_i03[ithread] = cy_r06->d1;
			_cy_r04[ithread] = cy_r08->d0;	_cy_i04[ithread] = cy_r08->d1;
			_cy_r05[ithread] = cy_r10->d0;	_cy_i05[ithread] = cy_r10->d1;
			_cy_r06[ithread] = cy_r12->d0;	_cy_i06[ithread] = cy_r12->d1;
			_cy_r07[ithread] = cy_r14->d0;	_cy_i07[ithread] = cy_r14->d1;
			_cy_r08[ithread] = cy_r16->d0;	_cy_i08[ithread] = cy_r16->d1;
			_cy_r09[ithread] = cy_r18->d0;	_cy_i09[ithread] = cy_r18->d1;
			_cy_r10[ithread] = cy_r20->d0;	_cy_i10[ithread] = cy_r20->d1;
			_cy_r11[ithread] = cy_r22->d0;	_cy_i11[ithread] = cy_r22->d1;
			_cy_r12[ithread] = cy_r24->d0;	_cy_i12[ithread] = cy_r24->d1;
			_cy_r13[ithread] = cy_r26->d0;	_cy_i13[ithread] = cy_r26->d1;
			_cy_r14[ithread] = cy_i00->d0;	_cy_i14[ithread] = cy_i00->d1;
			_cy_r15[ithread] = cy_i02->d0;	_cy_i15[ithread] = cy_i02->d1;
			_cy_r16[ithread] = cy_i04->d0;	_cy_i16[ithread] = cy_i04->d1;
			_cy_r17[ithread] = cy_i06->d0;	_cy_i17[ithread] = cy_i06->d1;
			_cy_r18[ithread] = cy_i08->d0;	_cy_i18[ithread] = cy_i08->d1;
			_cy_r19[ithread] = cy_i10->d0;	_cy_i19[ithread] = cy_i10->d1;
			_cy_r20[ithread] = cy_i12->d0;	_cy_i20[ithread] = cy_i12->d1;
			_cy_r21[ithread] = cy_i14->d0;	_cy_i21[ithread] = cy_i14->d1;
			_cy_r22[ithread] = cy_i16->d0;	_cy_i22[ithread] = cy_i16->d1;
			_cy_r23[ithread] = cy_i18->d0;	_cy_i23[ithread] = cy_i18->d1;
			_cy_r24[ithread] = cy_i20->d0;	_cy_i24[ithread] = cy_i20->d1;
			_cy_r25[ithread] = cy_i22->d0;	_cy_i25[ithread] = cy_i22->d1;
			_cy_r26[ithread] = cy_i24->d0;	_cy_i26[ithread] = cy_i24->d1;
			_cy_r27[ithread] = cy_i26->d0;	_cy_i27[ithread] = cy_i26->d1;
			if(full_pass) maxerr = MAX(max_err->d0,max_err->d1);
		#else
			_cy_r00[ithread] = cy_r00;	_cy_i00[ithread] = cy_i00;
			_cy_r01[ithread] = cy_r01;	_cy_i01[ithread] = cy_i01;
			_cy_r02[ithread] = cy_r02;	_cy_i02[ithread] = cy_i02;
			_cy_r03[ithread] = cy_r03;	_cy_i03[ithread] = cy_i03;
			_cy_r04[ithread] = cy_r04;	_cy_i04[ithread] = cy_i04;
			_cy_r05[ithread] = cy_r05;	_cy_i05[ithread] = cy_i05;
			_cy_r06[ithread] = cy_r06;	_cy_i06[ithread] = cy_i06;
			_cy_r07[ithread] = cy_r07;	_cy_i07[ithread] = cy_i07;
			_cy_r08[ithread] = cy_r08;	_cy_i08[ithread] = cy_i08;
			_cy_r09[ithread] = cy_r09;	_cy_i09[ithread] = cy_i09;
			_cy_r10[ithread] = cy_r10;	_cy_i10[ithread] = cy_i10;
			_cy_r11[ithread] = cy_r11;	_cy_i11[ithread] = cy_i11;
			_cy_r12[ithread] = cy_r12;	_cy_i12[ithread] = cy_i12;
			_cy_r13[ithread] = cy_r13;	_cy_i13[ithread] = cy_i13;
			_cy_r14[ithread] = cy_r14;	_cy_i14[ithread] = cy_i14;
			_cy_r15[ithread] = cy_r15;	_cy_i15[ithread] = cy_i15;
			_cy_r16[ithread] = cy_r16;	_cy_i16[ithread] = cy_i16;
			_cy_r17[ithread] = cy_r17;	_cy_i17[ithread] = cy_i17;
			_cy_r18[ithread] = cy_r18;	_cy_i18[ithread] = cy_i18;
			_cy_r19[ithread] = cy_r19;	_cy_i19[ithread] = cy_i19;
			_cy_r20[ithread] = cy_r20;	_cy_i20[ithread] = cy_i20;
			_cy_r21[ithread] = cy_r21;	_cy_i21[ithread] = cy_i21;
			_cy_r22[ithread] = cy_r22;	_cy_i22[ithread] = cy_i22;
			_cy_r23[ithread] = cy_r23;	_cy_i23[ithread] = cy_i23;
			_cy_r24[ithread] = cy_r24;	_cy_i24[ithread] = cy_i24;
			_cy_r25[ithread] = cy_r25;	_cy_i25[ithread] = cy_i25;
			_cy_r26[ithread] = cy_r26;	_cy_i26[ithread] = cy_i26;
			_cy_r27[ithread] = cy_r27;	_cy_i27[ithread] = cy_i27;
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
		ASSERT(HERE, 0x0 == cy28_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
	ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
	ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(HERE, 0 == nanosleep(&ns_time, 0x0), "nanosleep fail!");
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
			_cy_r00[ithread] = tdat[ithread].cy_r00;
			_cy_r01[ithread] = tdat[ithread].cy_r01;
			_cy_r02[ithread] = tdat[ithread].cy_r02;
			_cy_r03[ithread] = tdat[ithread].cy_r03;
			_cy_r04[ithread] = tdat[ithread].cy_r04;
			_cy_r05[ithread] = tdat[ithread].cy_r05;
			_cy_r06[ithread] = tdat[ithread].cy_r06;
			_cy_r07[ithread] = tdat[ithread].cy_r07;
			_cy_r08[ithread] = tdat[ithread].cy_r08;
			_cy_r09[ithread] = tdat[ithread].cy_r09;
			_cy_r10[ithread] = tdat[ithread].cy_r10;
			_cy_r11[ithread] = tdat[ithread].cy_r11;
			_cy_r12[ithread] = tdat[ithread].cy_r12;
			_cy_r13[ithread] = tdat[ithread].cy_r13;
			_cy_r14[ithread] = tdat[ithread].cy_r14;
			_cy_r15[ithread] = tdat[ithread].cy_r15;
			_cy_r16[ithread] = tdat[ithread].cy_r16;
			_cy_r17[ithread] = tdat[ithread].cy_r17;
			_cy_r18[ithread] = tdat[ithread].cy_r18;
			_cy_r19[ithread] = tdat[ithread].cy_r19;
			_cy_r20[ithread] = tdat[ithread].cy_r20;
			_cy_r21[ithread] = tdat[ithread].cy_r21;
			_cy_r22[ithread] = tdat[ithread].cy_r22;
			_cy_r23[ithread] = tdat[ithread].cy_r23;
			_cy_r24[ithread] = tdat[ithread].cy_r24;
			_cy_r25[ithread] = tdat[ithread].cy_r25;
			_cy_r26[ithread] = tdat[ithread].cy_r26;
			_cy_r27[ithread] = tdat[ithread].cy_r27;
		}
		else
		{
			_cy_r00[ithread] = tdat[ithread].cy_r00;	_cy_i00[ithread] = tdat[ithread].cy_i00;
			_cy_r01[ithread] = tdat[ithread].cy_r01;	_cy_i01[ithread] = tdat[ithread].cy_i01;
			_cy_r02[ithread] = tdat[ithread].cy_r02;	_cy_i02[ithread] = tdat[ithread].cy_i02;
			_cy_r03[ithread] = tdat[ithread].cy_r03;	_cy_i03[ithread] = tdat[ithread].cy_i03;
			_cy_r04[ithread] = tdat[ithread].cy_r04;	_cy_i04[ithread] = tdat[ithread].cy_i04;
			_cy_r05[ithread] = tdat[ithread].cy_r05;	_cy_i05[ithread] = tdat[ithread].cy_i05;
			_cy_r06[ithread] = tdat[ithread].cy_r06;	_cy_i06[ithread] = tdat[ithread].cy_i06;
			_cy_r07[ithread] = tdat[ithread].cy_r07;	_cy_i07[ithread] = tdat[ithread].cy_i07;
			_cy_r08[ithread] = tdat[ithread].cy_r08;	_cy_i08[ithread] = tdat[ithread].cy_i08;
			_cy_r09[ithread] = tdat[ithread].cy_r09;	_cy_i09[ithread] = tdat[ithread].cy_i09;
			_cy_r10[ithread] = tdat[ithread].cy_r10;	_cy_i10[ithread] = tdat[ithread].cy_i10;
			_cy_r11[ithread] = tdat[ithread].cy_r11;	_cy_i11[ithread] = tdat[ithread].cy_i11;
			_cy_r12[ithread] = tdat[ithread].cy_r12;	_cy_i12[ithread] = tdat[ithread].cy_i12;
			_cy_r13[ithread] = tdat[ithread].cy_r13;	_cy_i13[ithread] = tdat[ithread].cy_i13;
			_cy_r14[ithread] = tdat[ithread].cy_r14;	_cy_i14[ithread] = tdat[ithread].cy_i14;
			_cy_r15[ithread] = tdat[ithread].cy_r15;	_cy_i15[ithread] = tdat[ithread].cy_i15;
			_cy_r16[ithread] = tdat[ithread].cy_r16;	_cy_i16[ithread] = tdat[ithread].cy_i16;
			_cy_r17[ithread] = tdat[ithread].cy_r17;	_cy_i17[ithread] = tdat[ithread].cy_i17;
			_cy_r18[ithread] = tdat[ithread].cy_r18;	_cy_i18[ithread] = tdat[ithread].cy_i18;
			_cy_r19[ithread] = tdat[ithread].cy_r19;	_cy_i19[ithread] = tdat[ithread].cy_i19;
			_cy_r20[ithread] = tdat[ithread].cy_r20;	_cy_i20[ithread] = tdat[ithread].cy_i20;
			_cy_r21[ithread] = tdat[ithread].cy_r21;	_cy_i21[ithread] = tdat[ithread].cy_i21;
			_cy_r22[ithread] = tdat[ithread].cy_r22;	_cy_i22[ithread] = tdat[ithread].cy_i22;
			_cy_r23[ithread] = tdat[ithread].cy_r23;	_cy_i23[ithread] = tdat[ithread].cy_i23;
			_cy_r24[ithread] = tdat[ithread].cy_r24;	_cy_i24[ithread] = tdat[ithread].cy_i24;
			_cy_r25[ithread] = tdat[ithread].cy_r25;	_cy_i25[ithread] = tdat[ithread].cy_i25;
			_cy_r26[ithread] = tdat[ithread].cy_r26;	_cy_i26[ithread] = tdat[ithread].cy_i26;
			_cy_r27[ithread] = tdat[ithread].cy_r27;	_cy_i27[ithread] = tdat[ithread].cy_i27;
		}
	}
#endif

	if(full_pass) {
//		printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

	/*   Wraparound carry cleanup loop is here: ***
	!
	!   (1) Invert the radix-28 forward DIF FFT of the first block of 28 complex elements in A and unweight;
	!   (2) Propagate cleanup carries among the real and imaginary parts of the 28 outputs of (1);
	!   (3) Reweight and perform a radix-28 forward DIF FFT on the result of (2);
	!   (4) If any of the exit carries from (2) are nonzero, advance to the next 28 elements and repeat (1-4).
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t00= _cy_r00[CY_THREADS - 1];
		t02= _cy_r01[CY_THREADS - 1];
		t04= _cy_r02[CY_THREADS - 1];
		t06= _cy_r03[CY_THREADS - 1];
		t08= _cy_r04[CY_THREADS - 1];
		t10= _cy_r05[CY_THREADS - 1];
		t12= _cy_r06[CY_THREADS - 1];
		t14= _cy_r07[CY_THREADS - 1];
		t16= _cy_r08[CY_THREADS - 1];
		t18= _cy_r09[CY_THREADS - 1];
		t20= _cy_r10[CY_THREADS - 1];
		t22= _cy_r11[CY_THREADS - 1];
		t24= _cy_r12[CY_THREADS - 1];
		t26= _cy_r13[CY_THREADS - 1];
		t28= _cy_r14[CY_THREADS - 1];
		t30= _cy_r15[CY_THREADS - 1];
		t32= _cy_r16[CY_THREADS - 1];
		t34= _cy_r17[CY_THREADS - 1];
		t36= _cy_r18[CY_THREADS - 1];
		t38= _cy_r19[CY_THREADS - 1];
		t40= _cy_r20[CY_THREADS - 1];
		t42= _cy_r21[CY_THREADS - 1];
		t44= _cy_r22[CY_THREADS - 1];
		t46= _cy_r23[CY_THREADS - 1];
		t48= _cy_r24[CY_THREADS - 1];
		t50= _cy_r25[CY_THREADS - 1];
		t52= _cy_r26[CY_THREADS - 1];
		t54= _cy_r27[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1,"radix28_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
			_cy_r00[ithread] = _cy_r00[ithread-1];
			_cy_r01[ithread] = _cy_r01[ithread-1];
			_cy_r02[ithread] = _cy_r02[ithread-1];
			_cy_r03[ithread] = _cy_r03[ithread-1];
			_cy_r04[ithread] = _cy_r04[ithread-1];
			_cy_r05[ithread] = _cy_r05[ithread-1];
			_cy_r06[ithread] = _cy_r06[ithread-1];
			_cy_r07[ithread] = _cy_r07[ithread-1];
			_cy_r08[ithread] = _cy_r08[ithread-1];
			_cy_r09[ithread] = _cy_r09[ithread-1];
			_cy_r10[ithread] = _cy_r10[ithread-1];
			_cy_r11[ithread] = _cy_r11[ithread-1];
			_cy_r12[ithread] = _cy_r12[ithread-1];
			_cy_r13[ithread] = _cy_r13[ithread-1];
			_cy_r14[ithread] = _cy_r14[ithread-1];
			_cy_r15[ithread] = _cy_r15[ithread-1];
			_cy_r16[ithread] = _cy_r16[ithread-1];
			_cy_r17[ithread] = _cy_r17[ithread-1];
			_cy_r18[ithread] = _cy_r18[ithread-1];
			_cy_r19[ithread] = _cy_r19[ithread-1];
			_cy_r20[ithread] = _cy_r20[ithread-1];
			_cy_r21[ithread] = _cy_r21[ithread-1];
			_cy_r22[ithread] = _cy_r22[ithread-1];
			_cy_r23[ithread] = _cy_r23[ithread-1];
			_cy_r24[ithread] = _cy_r24[ithread-1];
			_cy_r25[ithread] = _cy_r25[ithread-1];
			_cy_r26[ithread] = _cy_r26[ithread-1];
			_cy_r27[ithread] = _cy_r27[ithread-1];
		}

		_cy_r00[0] =+t54;	/* ...The wraparound carry is here: */
		_cy_r01[0] = t00;
		_cy_r02[0] = t02;
		_cy_r03[0] = t04;
		_cy_r04[0] = t06;
		_cy_r05[0] = t08;
		_cy_r06[0] = t10;
		_cy_r07[0] = t12;
		_cy_r08[0] = t14;
		_cy_r09[0] = t16;
		_cy_r10[0] = t18;
		_cy_r11[0] = t20;
		_cy_r12[0] = t22;
		_cy_r13[0] = t24;
		_cy_r14[0] = t26;
		_cy_r15[0] = t28;
		_cy_r16[0] = t30;
		_cy_r17[0] = t32;
		_cy_r18[0] = t34;
		_cy_r19[0] = t36;
		_cy_r20[0] = t38;
		_cy_r21[0] = t40;
		_cy_r22[0] = t42;
		_cy_r23[0] = t44;
		_cy_r24[0] = t46;
		_cy_r25[0] = t48;
		_cy_r26[0] = t50;
		_cy_r27[0] = t52;
	}
	else
	{
		t00= _cy_r00[CY_THREADS - 1];	t01= _cy_i00[CY_THREADS - 1];
		t02= _cy_r01[CY_THREADS - 1];	t03= _cy_i01[CY_THREADS - 1];
		t04= _cy_r02[CY_THREADS - 1];	t05= _cy_i02[CY_THREADS - 1];
		t06= _cy_r03[CY_THREADS - 1];	t07= _cy_i03[CY_THREADS - 1];
		t08= _cy_r04[CY_THREADS - 1];	t09= _cy_i04[CY_THREADS - 1];
		t10= _cy_r05[CY_THREADS - 1];	t11= _cy_i05[CY_THREADS - 1];
		t12= _cy_r06[CY_THREADS - 1];	t13= _cy_i06[CY_THREADS - 1];
		t14= _cy_r07[CY_THREADS - 1];	t15= _cy_i07[CY_THREADS - 1];
		t16= _cy_r08[CY_THREADS - 1];	t17= _cy_i08[CY_THREADS - 1];
		t18= _cy_r09[CY_THREADS - 1];	t19= _cy_i09[CY_THREADS - 1];
		t20= _cy_r10[CY_THREADS - 1];	t21= _cy_i10[CY_THREADS - 1];
		t22= _cy_r11[CY_THREADS - 1];	t23= _cy_i11[CY_THREADS - 1];
		t24= _cy_r12[CY_THREADS - 1];	t25= _cy_i12[CY_THREADS - 1];
		t26= _cy_r13[CY_THREADS - 1];	t27= _cy_i13[CY_THREADS - 1];
		t28= _cy_r14[CY_THREADS - 1];	t29= _cy_i14[CY_THREADS - 1];
		t30= _cy_r15[CY_THREADS - 1];	t31= _cy_i15[CY_THREADS - 1];
		t32= _cy_r16[CY_THREADS - 1];	t33= _cy_i16[CY_THREADS - 1];
		t34= _cy_r17[CY_THREADS - 1];	t35= _cy_i17[CY_THREADS - 1];
		t36= _cy_r18[CY_THREADS - 1];	t37= _cy_i18[CY_THREADS - 1];
		t38= _cy_r19[CY_THREADS - 1];	t39= _cy_i19[CY_THREADS - 1];
		t40= _cy_r20[CY_THREADS - 1];	t41= _cy_i20[CY_THREADS - 1];
		t42= _cy_r21[CY_THREADS - 1];	t43= _cy_i21[CY_THREADS - 1];
		t44= _cy_r22[CY_THREADS - 1];	t45= _cy_i22[CY_THREADS - 1];
		t46= _cy_r23[CY_THREADS - 1];	t47= _cy_i23[CY_THREADS - 1];
		t48= _cy_r24[CY_THREADS - 1];	t49= _cy_i24[CY_THREADS - 1];
		t50= _cy_r25[CY_THREADS - 1];	t51= _cy_i25[CY_THREADS - 1];
		t52= _cy_r26[CY_THREADS - 1];	t53= _cy_i26[CY_THREADS - 1];
		t54= _cy_r27[CY_THREADS - 1];	t55= _cy_i27[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1,"radix28_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
			_cy_r00[ithread] = _cy_r00[ithread-1];		_cy_i00[ithread] = _cy_i00[ithread-1];
			_cy_r01[ithread] = _cy_r01[ithread-1];		_cy_i01[ithread] = _cy_i01[ithread-1];
			_cy_r02[ithread] = _cy_r02[ithread-1];		_cy_i02[ithread] = _cy_i02[ithread-1];
			_cy_r03[ithread] = _cy_r03[ithread-1];		_cy_i03[ithread] = _cy_i03[ithread-1];
			_cy_r04[ithread] = _cy_r04[ithread-1];		_cy_i04[ithread] = _cy_i04[ithread-1];
			_cy_r05[ithread] = _cy_r05[ithread-1];		_cy_i05[ithread] = _cy_i05[ithread-1];
			_cy_r06[ithread] = _cy_r06[ithread-1];		_cy_i06[ithread] = _cy_i06[ithread-1];
			_cy_r07[ithread] = _cy_r07[ithread-1];		_cy_i07[ithread] = _cy_i07[ithread-1];
			_cy_r08[ithread] = _cy_r08[ithread-1];		_cy_i08[ithread] = _cy_i08[ithread-1];
			_cy_r09[ithread] = _cy_r09[ithread-1];		_cy_i09[ithread] = _cy_i09[ithread-1];
			_cy_r10[ithread] = _cy_r10[ithread-1];		_cy_i10[ithread] = _cy_i10[ithread-1];
			_cy_r11[ithread] = _cy_r11[ithread-1];		_cy_i11[ithread] = _cy_i11[ithread-1];
			_cy_r12[ithread] = _cy_r12[ithread-1];		_cy_i12[ithread] = _cy_i12[ithread-1];
			_cy_r13[ithread] = _cy_r13[ithread-1];		_cy_i13[ithread] = _cy_i13[ithread-1];
			_cy_r14[ithread] = _cy_r14[ithread-1];		_cy_i14[ithread] = _cy_i14[ithread-1];
			_cy_r15[ithread] = _cy_r15[ithread-1];		_cy_i15[ithread] = _cy_i15[ithread-1];
			_cy_r16[ithread] = _cy_r16[ithread-1];		_cy_i16[ithread] = _cy_i16[ithread-1];
			_cy_r17[ithread] = _cy_r17[ithread-1];		_cy_i17[ithread] = _cy_i17[ithread-1];
			_cy_r18[ithread] = _cy_r18[ithread-1];		_cy_i18[ithread] = _cy_i18[ithread-1];
			_cy_r19[ithread] = _cy_r19[ithread-1];		_cy_i19[ithread] = _cy_i19[ithread-1];
			_cy_r20[ithread] = _cy_r20[ithread-1];		_cy_i20[ithread] = _cy_i20[ithread-1];
			_cy_r21[ithread] = _cy_r21[ithread-1];		_cy_i21[ithread] = _cy_i21[ithread-1];
			_cy_r22[ithread] = _cy_r22[ithread-1];		_cy_i22[ithread] = _cy_i22[ithread-1];
			_cy_r23[ithread] = _cy_r23[ithread-1];		_cy_i23[ithread] = _cy_i23[ithread-1];
			_cy_r24[ithread] = _cy_r24[ithread-1];		_cy_i24[ithread] = _cy_i24[ithread-1];
			_cy_r25[ithread] = _cy_r25[ithread-1];		_cy_i25[ithread] = _cy_i25[ithread-1];
			_cy_r26[ithread] = _cy_r26[ithread-1];		_cy_i26[ithread] = _cy_i26[ithread-1];
			_cy_r27[ithread] = _cy_r27[ithread-1];		_cy_i27[ithread] = _cy_i27[ithread-1];
		}

		_cy_r00[0] =-t55;	_cy_i00[0] =+t54;	/* ...The 2 Mo"bius carries are here: */
		_cy_r01[0] = t00;	_cy_i01[0] = t01;
		_cy_r02[0] = t02;	_cy_i02[0] = t03;
		_cy_r03[0] = t04;	_cy_i03[0] = t05;
		_cy_r04[0] = t06;	_cy_i04[0] = t07;
		_cy_r05[0] = t08;	_cy_i05[0] = t09;
		_cy_r06[0] = t10;	_cy_i06[0] = t11;
		_cy_r07[0] = t12;	_cy_i07[0] = t13;
		_cy_r08[0] = t14;	_cy_i08[0] = t15;
		_cy_r09[0] = t16;	_cy_i09[0] = t17;
		_cy_r10[0] = t18;	_cy_i10[0] = t19;
		_cy_r11[0] = t20;	_cy_i11[0] = t21;
		_cy_r12[0] = t22;	_cy_i12[0] = t23;
		_cy_r13[0] = t24;	_cy_i13[0] = t25;
		_cy_r14[0] = t26;	_cy_i14[0] = t27;
		_cy_r15[0] = t28;	_cy_i15[0] = t29;
		_cy_r16[0] = t30;	_cy_i16[0] = t31;
		_cy_r17[0] = t32;	_cy_i17[0] = t33;
		_cy_r18[0] = t34;	_cy_i18[0] = t35;
		_cy_r19[0] = t36;	_cy_i19[0] = t37;
		_cy_r20[0] = t38;	_cy_i20[0] = t39;
		_cy_r21[0] = t40;	_cy_i21[0] = t41;
		_cy_r22[0] = t42;	_cy_i22[0] = t43;
		_cy_r23[0] = t44;	_cy_i23[0] = t45;
		_cy_r24[0] = t46;	_cy_i24[0] = t47;
		_cy_r25[0] = t48;	_cy_i25[0] = t49;
		_cy_r26[0] = t50;	_cy_i26[0] = t51;
		_cy_r27[0] = t52;	_cy_i27[0] = t53;
	}

	full_pass = 0;
	scale = 1;

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
			j2 = j + ( (j >> DAT_BITS) << PAD_BITS );
			jt = j2;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j2 + p04;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j2 + p08;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j2 + p12;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j2 + p16;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j2 + p20;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j2 + p24;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
		}
	}
}	/* endfor(outer) */

	t00 = 0;
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		t00 += fabs(_cy_r00[0])+fabs(_cy_r01[0])+fabs(_cy_r02[0])+fabs(_cy_r03[0])+fabs(_cy_r04[0])+fabs(_cy_r05[0])+fabs(_cy_r06[0])+fabs(_cy_r07[0])+fabs(_cy_r08[0])+fabs(_cy_r09[0])+fabs(_cy_r10[0])+fabs(_cy_r11[0])+fabs(_cy_r12[0])+fabs(_cy_r13[0])+fabs(_cy_r14[0])+fabs(_cy_r15[0])+fabs(_cy_r16[0])+fabs(_cy_r17[0])+fabs(_cy_r18[0])+fabs(_cy_r19[0])+fabs(_cy_r20[0])+fabs(_cy_r21[0])+fabs(_cy_r22[0])+fabs(_cy_r23[0])+fabs(_cy_r24[0])+fabs(_cy_r25[0])+fabs(_cy_r26[0])+fabs(_cy_r27[0]);
		t00 += fabs(_cy_i00[0])+fabs(_cy_i01[0])+fabs(_cy_i02[0])+fabs(_cy_i03[0])+fabs(_cy_i04[0])+fabs(_cy_i05[0])+fabs(_cy_i06[0])+fabs(_cy_i07[0])+fabs(_cy_i08[0])+fabs(_cy_i09[0])+fabs(_cy_i10[0])+fabs(_cy_i11[0])+fabs(_cy_i12[0])+fabs(_cy_i13[0])+fabs(_cy_i14[0])+fabs(_cy_i15[0])+fabs(_cy_i16[0])+fabs(_cy_i17[0])+fabs(_cy_i18[0])+fabs(_cy_i19[0])+fabs(_cy_i20[0])+fabs(_cy_i21[0])+fabs(_cy_i22[0])+fabs(_cy_i23[0])+fabs(_cy_i24[0])+fabs(_cy_i25[0])+fabs(_cy_i26[0])+fabs(_cy_i27[0]);

		*fracmax = maxerr;
	}

	if(t00 != 0.0)
	{
		sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in %s - input wordsize may be too small.\n",iter,func);
		if(INTERACT)fprintf(stderr,"%s",cbuf);
		fp = mlucas_fopen(   OFILE,"a");
		fq = mlucas_fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;
		err=ERR_CARRY;
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
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;

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
!...Acronym: DIF = Decimation In Time
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
		const char func[] = "radix28_ditN_cy_dif1";
		const int pfetch_dist = PFETCH_DIST;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		int j,j1,j2,k,l;
		uint32 p01,p02,p03,p04,p08,p12,p16,p20,p24;
		int poff[RADIX>>2];	// Store [RADIX/4] mults of p04 offset for loop control
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
	#ifdef USE_AVX
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#else
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	#endif
		double *addr;
		double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
		double rt,it;
		int k1,k2;

	#ifdef USE_SSE2

	  #ifdef USE_AVX
		const int l2_sz_vd = 5;
	  #else
		const int l2_sz_vd = 4;
	  #endif
		const double crnd = 3.0*0x4000000*0x2000000;
		int *itmp;	// Pointer into the bjmodn array
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		const double sx3 = 0.87484229096165655224;
		double *add0, *add1, *add2, *add3;
		vec_dbl *one,*two, *cc0, *ss0, *cc1, *ss1, *cc2, *ss2, *cc3, *ss3, *max_err, *sse2_rnd, *half_arr, *tmp,*tm1,*tm2,
			*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,
			*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r;
		vec_dbl *cy_r00,*cy_r04,*cy_r08,*cy_r12,*cy_r16,*cy_r20,*cy_r24,*cy_i00,*cy_i04,*cy_i08,*cy_i12,*cy_i16,*cy_i20,*cy_i24;
	  #ifndef USE_AVX
		vec_dbl *cy_r02,*cy_r06,*cy_r10,*cy_r14,*cy_r18,*cy_r22,*cy_r26,*cy_i02,*cy_i06,*cy_i10,*cy_i14,*cy_i18,*cy_i22,*cy_i26;
	  #else
		vec_dbl *base_negacyclic_root;
	  #endif
		int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27;
		/* These are used in conjunction with the langth-7 arrays in the USE_SCALAR_CARRY #define below;
		In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
		int idx_offset, idx_incr;
		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else

	  // LO_ADD = 1 in masterdefs.h means non-Nussbaumer radix-7, uses these sincos constants:
	  #ifdef LO_ADD
		const double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
						us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
						uc2 =-.22252093395631440426,	 /* cos(2u)	*/
						us2 = .97492791218182360702,	 /* sin(2u)	*/
						uc3 =-.90096886790241912622,	 /* cos(3u)	*/
						us3 = .43388373911755812050;	 /* sin(3u)	*/
	  #else	// Low-MUL/high-ADD Nussbaumer radix-7 needs these trig combos:
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
		double *base, *baseinv, *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr, bs,bsinv;
		int wts_idx_incr;
		const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int jt,jp,m,m2,ntmp;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	  #if PFETCH
		double *addp;
	  #endif
		double re,im,temp,frac
			,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13
			,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r
			,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i
			,cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27
			,cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27;
		int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27;

	#endif

		struct cy_thread_data_t* thread_arg = targ;

	// int data:
		int iter = thread_arg->iter;
		int NDIVR = thread_arg->ndivr;
		int n = NDIVR*RADIX;
		int khi    = thread_arg->khi;
		int i      = thread_arg->i;	/* Pointer to the BASE and BASEINV arrays.	*/
		int jstart = thread_arg->jstart;
		int jhi    = thread_arg->jhi;
		int col = thread_arg->col;
		int co2 = thread_arg->co2;
		int co3 = thread_arg->co3;
		int sw  = thread_arg->sw;
		int nwt = thread_arg->nwt;
		int icycle0 = thread_arg->icycle0;
		int icycle1 = thread_arg->icycle1;
		int icycle2 = thread_arg->icycle2;
		int icycle3 = thread_arg->icycle3;
		int icycle4 = thread_arg->icycle4;
		int icycle5 = thread_arg->icycle5;
		int icycle6 = thread_arg->icycle6;
	#ifdef USE_SSE2
		int wts_idx_inc2 = thread_arg->wts_idx_inc2;
		int jcycle0 = thread_arg->jcycle0;
		int jcycle1 = thread_arg->jcycle1;
		int jcycle2 = thread_arg->jcycle2;
		int jcycle3 = thread_arg->jcycle3;
		int jcycle4 = thread_arg->jcycle4;
		int jcycle5 = thread_arg->jcycle5;
		int jcycle6 = thread_arg->jcycle6;
	  #ifdef USE_AVX
		int kcycle0 = thread_arg->kcycle0;
		int kcycle1 = thread_arg->kcycle1;
		int kcycle2 = thread_arg->kcycle2;
		int kcycle3 = thread_arg->kcycle3;
		int kcycle4 = thread_arg->kcycle4;
		int kcycle5 = thread_arg->kcycle5;
		int kcycle6 = thread_arg->kcycle6;

		int lcycle0 = thread_arg->lcycle0;
		int lcycle1 = thread_arg->lcycle1;
		int lcycle2 = thread_arg->lcycle2;
		int lcycle3 = thread_arg->lcycle3;
		int lcycle4 = thread_arg->lcycle4;
		int lcycle5 = thread_arg->lcycle5;
		int lcycle6 = thread_arg->lcycle6;
	  #endif
	#endif

	// double data:
		double maxerr = thread_arg->maxerr;
		double scale = thread_arg->scale;	int full_pass = scale < 0.5;

	// pointer data:
		double *a = thread_arg->arrdat;
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
	#ifdef LOACC
		double *wts_mult = thread_arg->wts_mult;	// Const Intra-block wts-multiplier...
		double *inv_mult = thread_arg->inv_mult;	// ...and 2*(its multiplicative inverse).
		ASSERT(HERE,fabs(wts_mult[0]*inv_mult[0] - 1.0) < EPS, "wts_mults fail accuracy check!");
		ASSERT(HERE,fabs(wts_mult[1]*inv_mult[1] - 1.0) < EPS, "wts_mults fail accuracy check!");
	#endif
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

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );

		poff[0] =   0; poff[1] = p04; poff[2] = p08; poff[3] = p12; poff[4] = p16; poff[5] = p20; poff[6] = p24;

	#ifdef USE_SSE2
		uint32 nwt16 = nwt << l2_sz_vd;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code

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
		cy_r00	= tmp + 0x0;
		cy_r04	= tmp + 0x1;
		cy_r08	= tmp + 0x2;
		cy_r12	= tmp + 0x3;
		cy_r16	= tmp + 0x4;
		cy_r20	= tmp + 0x5;
		cy_r24	= tmp + 0x6;
		cy_i00	= tmp + 0x7;
		cy_i04	= tmp + 0x8;
		cy_i08	= tmp + 0x9;
		cy_i12	= tmp + 0xa;
		cy_i16	= tmp + 0xb;
		cy_i20	= tmp + 0xc;
		cy_i24	= tmp + 0xd;
		tmp += 0xe;
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 2 = 0x58 [avx] or 0x66 [sse2]; This is where the value of half_arr_offset28 comes from
		half_arr= tmp + 0x02;	/* This table needs 20x16 bytes */
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod
	  #else
								tm2 = tmp + 0xe;
		cy_r00	= tmp + 0x0;	cy_i00	= tm2 + 0x0;
		cy_r02	= tmp + 0x1;	cy_i02	= tm2 + 0x1;
		cy_r04	= tmp + 0x2;	cy_i04	= tm2 + 0x2;
		cy_r06	= tmp + 0x3;	cy_i06	= tm2 + 0x3;
		cy_r08	= tmp + 0x4;	cy_i08	= tm2 + 0x4;
		cy_r10	= tmp + 0x5;	cy_i10	= tm2 + 0x5;
		cy_r12	= tmp + 0x6;	cy_i12	= tm2 + 0x6;
		cy_r14	= tmp + 0x7;	cy_i14	= tm2 + 0x7;
		cy_r16	= tmp + 0x8;	cy_i16	= tm2 + 0x8;
		cy_r18	= tmp + 0x9;	cy_i18	= tm2 + 0x9;
		cy_r20	= tmp + 0xa;	cy_i20	= tm2 + 0xa;
		cy_r22	= tmp + 0xb;	cy_i22	= tm2 + 0xb;
		cy_r24	= tmp + 0xc;	cy_i24	= tm2 + 0xc;
		cy_r26	= tmp + 0xd;	cy_i26	= tm2 + 0xd;
		tmp += 0x1c;
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 2 = 0x58 [avx] or 0x66 [sse2]; This is where the value of half_arr_offset28 comes from
		half_arr= tmp + 0x02;	/* This table needs 20x16 bytes */
	  #endif

		ASSERT(HERE, (two->d0 == 2.0 && two->d1 == 2.0), "thread-local memcheck failed!");
	  #ifdef USE_AVX2
		// AVX2 (i.e. FMA)means non-Nussbaumer radix-7, uses these sincos constants:
		ASSERT(HERE, (ss3->d0 == 0.0 && ss3->d1 == 0.0), "thread-local memcheck failed!");
	  #else
		/* SSE2 version assumes LO_ADD = 0, i.e. the low-mul Nussbaumer-style DFT implementation: */
		ASSERT(HERE, (ss3->d0 == sx3 && ss3->d1 == sx3), "thread-local memcheck failed!");
	  #endif
		ASSERT(HERE, (sse2_rnd->d0 == crnd && sse2_rnd->d1 == crnd), "thread-local memcheck failed!");
		tmp = half_arr;
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
	  #ifdef USE_AVX
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif
	} else {
		dtmp = (tmp)->d0 * (tmp+ODD_RADIX)->d0;	ASSERT(HERE, fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp)->d1 * (tmp+ODD_RADIX)->d1;	ASSERT(HERE, fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
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

		bjmodn00 = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn00 = (int*)(sse_n + RE_IM_STRIDE);
	  #endif
		bjmodn01 = bjmodn00 +  1;
		bjmodn02 = bjmodn00 +  2;
		bjmodn03 = bjmodn00 +  3;
		bjmodn04 = bjmodn00 +  4;
		bjmodn05 = bjmodn00 +  5;
		bjmodn06 = bjmodn00 +  6;
		bjmodn07 = bjmodn00 +  7;
		bjmodn08 = bjmodn00 +  8;
		bjmodn09 = bjmodn00 +  9;
		bjmodn10 = bjmodn00 + 10;
		bjmodn11 = bjmodn00 + 11;
		bjmodn12 = bjmodn00 + 12;
		bjmodn13 = bjmodn00 + 13;
		bjmodn14 = bjmodn00 + 14;
		bjmodn15 = bjmodn00 + 15;
		bjmodn16 = bjmodn00 + 16;
		bjmodn17 = bjmodn00 + 17;
		bjmodn18 = bjmodn00 + 18;
		bjmodn19 = bjmodn00 + 19;
		bjmodn20 = bjmodn00 + 20;
		bjmodn21 = bjmodn00 + 21;
		bjmodn22 = bjmodn00 + 22;
		bjmodn23 = bjmodn00 + 23;
		bjmodn24 = bjmodn00 + 24;
		bjmodn25 = bjmodn00 + 25;
		bjmodn26 = bjmodn00 + 26;
		bjmodn27 = bjmodn00 + 27;
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

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Init DWT-indices: */							/* init carries	*/
		#ifdef USE_AVX
			*bjmodn00 = thread_arg->bjmodn00;			cy_r00->d0 = thread_arg->cy_r00;
			*bjmodn01 = thread_arg->bjmodn01;			cy_r00->d1 = thread_arg->cy_r01;
			*bjmodn02 = thread_arg->bjmodn02;			cy_r00->d2 = thread_arg->cy_r02;
			*bjmodn03 = thread_arg->bjmodn03;			cy_r00->d3 = thread_arg->cy_r03;
			*bjmodn04 = thread_arg->bjmodn04;			cy_r04->d0 = thread_arg->cy_r04;
			*bjmodn05 = thread_arg->bjmodn05;			cy_r04->d1 = thread_arg->cy_r05;
			*bjmodn06 = thread_arg->bjmodn06;			cy_r04->d2 = thread_arg->cy_r06;
			*bjmodn07 = thread_arg->bjmodn07;			cy_r04->d3 = thread_arg->cy_r07;
			*bjmodn08 = thread_arg->bjmodn08;			cy_r08->d0 = thread_arg->cy_r08;
			*bjmodn09 = thread_arg->bjmodn09;			cy_r08->d1 = thread_arg->cy_r09;
			*bjmodn10 = thread_arg->bjmodn10;			cy_r08->d2 = thread_arg->cy_r10;
			*bjmodn11 = thread_arg->bjmodn11;			cy_r08->d3 = thread_arg->cy_r11;
			*bjmodn12 = thread_arg->bjmodn12;			cy_r12->d0 = thread_arg->cy_r12;
			*bjmodn13 = thread_arg->bjmodn13;			cy_r12->d1 = thread_arg->cy_r13;
			*bjmodn14 = thread_arg->bjmodn14;			cy_r12->d2 = thread_arg->cy_r14;
			*bjmodn15 = thread_arg->bjmodn15;			cy_r12->d3 = thread_arg->cy_r15;
			*bjmodn16 = thread_arg->bjmodn16;			cy_r16->d0 = thread_arg->cy_r16;
			*bjmodn17 = thread_arg->bjmodn17;			cy_r16->d1 = thread_arg->cy_r17;
			*bjmodn18 = thread_arg->bjmodn18;			cy_r16->d2 = thread_arg->cy_r18;
			*bjmodn19 = thread_arg->bjmodn19;			cy_r16->d3 = thread_arg->cy_r19;
			*bjmodn20 = thread_arg->bjmodn20;			cy_r20->d0 = thread_arg->cy_r20;
			*bjmodn21 = thread_arg->bjmodn21;			cy_r20->d1 = thread_arg->cy_r21;
			*bjmodn22 = thread_arg->bjmodn22;			cy_r20->d2 = thread_arg->cy_r22;
			*bjmodn23 = thread_arg->bjmodn23;			cy_r20->d3 = thread_arg->cy_r23;
			*bjmodn24 = thread_arg->bjmodn24;			cy_r24->d0 = thread_arg->cy_r24;
			*bjmodn25 = thread_arg->bjmodn25;			cy_r24->d1 = thread_arg->cy_r25;
			*bjmodn26 = thread_arg->bjmodn26;			cy_r24->d2 = thread_arg->cy_r26;
			*bjmodn27 = thread_arg->bjmodn27;			cy_r24->d3 = thread_arg->cy_r27;

		#elif defined(USE_SSE2)

			*bjmodn00 = thread_arg->bjmodn00;			cy_r00->d0 = thread_arg->cy_r00;
			*bjmodn01 = thread_arg->bjmodn01;			cy_r00->d1 = thread_arg->cy_r01;
			*bjmodn02 = thread_arg->bjmodn02;			cy_r02->d0 = thread_arg->cy_r02;
			*bjmodn03 = thread_arg->bjmodn03;			cy_r02->d1 = thread_arg->cy_r03;
			*bjmodn04 = thread_arg->bjmodn04;			cy_r04->d0 = thread_arg->cy_r04;
			*bjmodn05 = thread_arg->bjmodn05;			cy_r04->d1 = thread_arg->cy_r05;
			*bjmodn06 = thread_arg->bjmodn06;			cy_r06->d0 = thread_arg->cy_r06;
			*bjmodn07 = thread_arg->bjmodn07;			cy_r06->d1 = thread_arg->cy_r07;
			*bjmodn08 = thread_arg->bjmodn08;			cy_r08->d0 = thread_arg->cy_r08;
			*bjmodn09 = thread_arg->bjmodn09;			cy_r08->d1 = thread_arg->cy_r09;
			*bjmodn10 = thread_arg->bjmodn10;			cy_r10->d0 = thread_arg->cy_r10;
			*bjmodn11 = thread_arg->bjmodn11;			cy_r10->d1 = thread_arg->cy_r11;
			*bjmodn12 = thread_arg->bjmodn12;			cy_r12->d0 = thread_arg->cy_r12;
			*bjmodn13 = thread_arg->bjmodn13;			cy_r12->d1 = thread_arg->cy_r13;
			*bjmodn14 = thread_arg->bjmodn14;			cy_r14->d0 = thread_arg->cy_r14;
			*bjmodn15 = thread_arg->bjmodn15;			cy_r14->d1 = thread_arg->cy_r15;
			*bjmodn16 = thread_arg->bjmodn16;			cy_r16->d0 = thread_arg->cy_r16;
			*bjmodn17 = thread_arg->bjmodn17;			cy_r16->d1 = thread_arg->cy_r17;
			*bjmodn18 = thread_arg->bjmodn18;			cy_r18->d0 = thread_arg->cy_r18;
			*bjmodn19 = thread_arg->bjmodn19;			cy_r18->d1 = thread_arg->cy_r19;
			*bjmodn20 = thread_arg->bjmodn20;			cy_r20->d0 = thread_arg->cy_r20;
			*bjmodn21 = thread_arg->bjmodn21;			cy_r20->d1 = thread_arg->cy_r21;
			*bjmodn22 = thread_arg->bjmodn22;			cy_r22->d0 = thread_arg->cy_r22;
			*bjmodn23 = thread_arg->bjmodn23;			cy_r22->d1 = thread_arg->cy_r23;
			*bjmodn24 = thread_arg->bjmodn24;			cy_r24->d0 = thread_arg->cy_r24;
			*bjmodn25 = thread_arg->bjmodn25;			cy_r24->d1 = thread_arg->cy_r25;
			*bjmodn26 = thread_arg->bjmodn26;			cy_r26->d0 = thread_arg->cy_r26;
			*bjmodn27 = thread_arg->bjmodn27;			cy_r26->d1 = thread_arg->cy_r27;

		#else

			bjmodn00 = thread_arg->bjmodn00;			cy_r00 = thread_arg->cy_r00;
			bjmodn01 = thread_arg->bjmodn01;			cy_r01 = thread_arg->cy_r01;
			bjmodn02 = thread_arg->bjmodn02;			cy_r02 = thread_arg->cy_r02;
			bjmodn03 = thread_arg->bjmodn03;			cy_r03 = thread_arg->cy_r03;
			bjmodn04 = thread_arg->bjmodn04;			cy_r04 = thread_arg->cy_r04;
			bjmodn05 = thread_arg->bjmodn05;			cy_r05 = thread_arg->cy_r05;
			bjmodn06 = thread_arg->bjmodn06;			cy_r06 = thread_arg->cy_r06;
			bjmodn07 = thread_arg->bjmodn07;			cy_r07 = thread_arg->cy_r07;
			bjmodn08 = thread_arg->bjmodn08;			cy_r08 = thread_arg->cy_r08;
			bjmodn09 = thread_arg->bjmodn09;			cy_r09 = thread_arg->cy_r09;
			bjmodn10 = thread_arg->bjmodn10;			cy_r10 = thread_arg->cy_r10;
			bjmodn11 = thread_arg->bjmodn11;			cy_r11 = thread_arg->cy_r11;
			bjmodn12 = thread_arg->bjmodn12;			cy_r12 = thread_arg->cy_r12;
			bjmodn13 = thread_arg->bjmodn13;			cy_r13 = thread_arg->cy_r13;
			bjmodn14 = thread_arg->bjmodn14;			cy_r14 = thread_arg->cy_r14;
			bjmodn15 = thread_arg->bjmodn15;			cy_r15 = thread_arg->cy_r15;
			bjmodn16 = thread_arg->bjmodn16;			cy_r16 = thread_arg->cy_r16;
			bjmodn17 = thread_arg->bjmodn17;			cy_r17 = thread_arg->cy_r17;
			bjmodn18 = thread_arg->bjmodn18;			cy_r18 = thread_arg->cy_r18;
			bjmodn19 = thread_arg->bjmodn19;			cy_r19 = thread_arg->cy_r19;
			bjmodn20 = thread_arg->bjmodn20;			cy_r20 = thread_arg->cy_r20;
			bjmodn21 = thread_arg->bjmodn21;			cy_r21 = thread_arg->cy_r21;
			bjmodn22 = thread_arg->bjmodn22;			cy_r22 = thread_arg->cy_r22;
			bjmodn23 = thread_arg->bjmodn23;			cy_r23 = thread_arg->cy_r23;
			bjmodn24 = thread_arg->bjmodn24;			cy_r24 = thread_arg->cy_r24;
			bjmodn25 = thread_arg->bjmodn25;			cy_r25 = thread_arg->cy_r25;
			bjmodn26 = thread_arg->bjmodn26;			cy_r26 = thread_arg->cy_r26;
			bjmodn27 = thread_arg->bjmodn27;			cy_r27 = thread_arg->cy_r27;

		#endif
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
			/* init carries	*/
		#ifdef USE_AVX
			cy_r00->d0 = thread_arg->cy_r00;	cy_i00->d0 = thread_arg->cy_i00;
			cy_r00->d1 = thread_arg->cy_r01;	cy_i00->d1 = thread_arg->cy_i01;
			cy_r00->d2 = thread_arg->cy_r02;	cy_i00->d2 = thread_arg->cy_i02;
			cy_r00->d3 = thread_arg->cy_r03;	cy_i00->d3 = thread_arg->cy_i03;
			cy_r04->d0 = thread_arg->cy_r04;	cy_i04->d0 = thread_arg->cy_i04;
			cy_r04->d1 = thread_arg->cy_r05;	cy_i04->d1 = thread_arg->cy_i05;
			cy_r04->d2 = thread_arg->cy_r06;	cy_i04->d2 = thread_arg->cy_i06;
			cy_r04->d3 = thread_arg->cy_r07;	cy_i04->d3 = thread_arg->cy_i07;
			cy_r08->d0 = thread_arg->cy_r08;	cy_i08->d0 = thread_arg->cy_i08;
			cy_r08->d1 = thread_arg->cy_r09;	cy_i08->d1 = thread_arg->cy_i09;
			cy_r08->d2 = thread_arg->cy_r10;	cy_i08->d2 = thread_arg->cy_i10;
			cy_r08->d3 = thread_arg->cy_r11;	cy_i08->d3 = thread_arg->cy_i11;
			cy_r12->d0 = thread_arg->cy_r12;	cy_i12->d0 = thread_arg->cy_i12;
			cy_r12->d1 = thread_arg->cy_r13;	cy_i12->d1 = thread_arg->cy_i13;
			cy_r12->d2 = thread_arg->cy_r14;	cy_i12->d2 = thread_arg->cy_i14;
			cy_r12->d3 = thread_arg->cy_r15;	cy_i12->d3 = thread_arg->cy_i15;
			cy_r16->d0 = thread_arg->cy_r16;	cy_i16->d0 = thread_arg->cy_i16;
			cy_r16->d1 = thread_arg->cy_r17;	cy_i16->d1 = thread_arg->cy_i17;
			cy_r16->d2 = thread_arg->cy_r18;	cy_i16->d2 = thread_arg->cy_i18;
			cy_r16->d3 = thread_arg->cy_r19;	cy_i16->d3 = thread_arg->cy_i19;
			cy_r20->d0 = thread_arg->cy_r20;	cy_i20->d0 = thread_arg->cy_i20;
			cy_r20->d1 = thread_arg->cy_r21;	cy_i20->d1 = thread_arg->cy_i21;
			cy_r20->d2 = thread_arg->cy_r22;	cy_i20->d2 = thread_arg->cy_i22;
			cy_r20->d3 = thread_arg->cy_r23;	cy_i20->d3 = thread_arg->cy_i23;
			cy_r24->d0 = thread_arg->cy_r24;	cy_i24->d0 = thread_arg->cy_i24;
			cy_r24->d1 = thread_arg->cy_r25;	cy_i24->d1 = thread_arg->cy_i25;
			cy_r24->d2 = thread_arg->cy_r26;	cy_i24->d2 = thread_arg->cy_i26;
			cy_r24->d3 = thread_arg->cy_r27;	cy_i24->d3 = thread_arg->cy_i27;

		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			cy_r00->d0 = thread_arg->cy_r00;	cy_r00->d1 = thread_arg->cy_i00;
			cy_r02->d0 = thread_arg->cy_r01;	cy_r02->d1 = thread_arg->cy_i01;
			cy_r04->d0 = thread_arg->cy_r02;	cy_r04->d1 = thread_arg->cy_i02;
			cy_r06->d0 = thread_arg->cy_r03;	cy_r06->d1 = thread_arg->cy_i03;
			cy_r08->d0 = thread_arg->cy_r04;	cy_r08->d1 = thread_arg->cy_i04;
			cy_r10->d0 = thread_arg->cy_r05;	cy_r10->d1 = thread_arg->cy_i05;
			cy_r12->d0 = thread_arg->cy_r06;	cy_r12->d1 = thread_arg->cy_i06;
			cy_r14->d0 = thread_arg->cy_r07;	cy_r14->d1 = thread_arg->cy_i07;
			cy_r16->d0 = thread_arg->cy_r08;	cy_r16->d1 = thread_arg->cy_i08;
			cy_r18->d0 = thread_arg->cy_r09;	cy_r18->d1 = thread_arg->cy_i09;
			cy_r20->d0 = thread_arg->cy_r10;	cy_r20->d1 = thread_arg->cy_i10;
			cy_r22->d0 = thread_arg->cy_r11;	cy_r22->d1 = thread_arg->cy_i11;
			cy_r24->d0 = thread_arg->cy_r12;	cy_r24->d1 = thread_arg->cy_i12;
			cy_r26->d0 = thread_arg->cy_r13;	cy_r26->d1 = thread_arg->cy_i13;
			cy_i00->d0 = thread_arg->cy_r14;	cy_i00->d1 = thread_arg->cy_i14;
			cy_i02->d0 = thread_arg->cy_r15;	cy_i02->d1 = thread_arg->cy_i15;
			cy_i04->d0 = thread_arg->cy_r16;	cy_i04->d1 = thread_arg->cy_i16;
			cy_i06->d0 = thread_arg->cy_r17;	cy_i06->d1 = thread_arg->cy_i17;
			cy_i08->d0 = thread_arg->cy_r18;	cy_i08->d1 = thread_arg->cy_i18;
			cy_i10->d0 = thread_arg->cy_r19;	cy_i10->d1 = thread_arg->cy_i19;
			cy_i12->d0 = thread_arg->cy_r20;	cy_i12->d1 = thread_arg->cy_i20;
			cy_i14->d0 = thread_arg->cy_r21;	cy_i14->d1 = thread_arg->cy_i21;
			cy_i16->d0 = thread_arg->cy_r22;	cy_i16->d1 = thread_arg->cy_i22;
			cy_i18->d0 = thread_arg->cy_r23;	cy_i18->d1 = thread_arg->cy_i23;
			cy_i20->d0 = thread_arg->cy_r24;	cy_i20->d1 = thread_arg->cy_i24;
			cy_i22->d0 = thread_arg->cy_r25;	cy_i22->d1 = thread_arg->cy_i25;
			cy_i24->d0 = thread_arg->cy_r26;	cy_i24->d1 = thread_arg->cy_i26;
			cy_i26->d0 = thread_arg->cy_r27;	cy_i26->d1 = thread_arg->cy_i27;

		#else

			cy_r00 = thread_arg->cy_r00;		cy_i00 = thread_arg->cy_i00;
			cy_r01 = thread_arg->cy_r01;		cy_i01 = thread_arg->cy_i01;
			cy_r02 = thread_arg->cy_r02;		cy_i02 = thread_arg->cy_i02;
			cy_r03 = thread_arg->cy_r03;		cy_i03 = thread_arg->cy_i03;
			cy_r04 = thread_arg->cy_r04;		cy_i04 = thread_arg->cy_i04;
			cy_r05 = thread_arg->cy_r05;		cy_i05 = thread_arg->cy_i05;
			cy_r06 = thread_arg->cy_r06;		cy_i06 = thread_arg->cy_i06;
			cy_r07 = thread_arg->cy_r07;		cy_i07 = thread_arg->cy_i07;
			cy_r08 = thread_arg->cy_r08;		cy_i08 = thread_arg->cy_i08;
			cy_r09 = thread_arg->cy_r09;		cy_i09 = thread_arg->cy_i09;
			cy_r10 = thread_arg->cy_r10;		cy_i10 = thread_arg->cy_i10;
			cy_r11 = thread_arg->cy_r11;		cy_i11 = thread_arg->cy_i11;
			cy_r12 = thread_arg->cy_r12;		cy_i12 = thread_arg->cy_i12;
			cy_r13 = thread_arg->cy_r13;		cy_i13 = thread_arg->cy_i13;
			cy_r14 = thread_arg->cy_r14;		cy_i14 = thread_arg->cy_i14;
			cy_r15 = thread_arg->cy_r15;		cy_i15 = thread_arg->cy_i15;
			cy_r16 = thread_arg->cy_r16;		cy_i16 = thread_arg->cy_i16;
			cy_r17 = thread_arg->cy_r17;		cy_i17 = thread_arg->cy_i17;
			cy_r18 = thread_arg->cy_r18;		cy_i18 = thread_arg->cy_i18;
			cy_r19 = thread_arg->cy_r19;		cy_i19 = thread_arg->cy_i19;
			cy_r20 = thread_arg->cy_r20;		cy_i20 = thread_arg->cy_i20;
			cy_r21 = thread_arg->cy_r21;		cy_i21 = thread_arg->cy_i21;
			cy_r22 = thread_arg->cy_r22;		cy_i22 = thread_arg->cy_i22;
			cy_r23 = thread_arg->cy_r23;		cy_i23 = thread_arg->cy_i23;
			cy_r24 = thread_arg->cy_r24;		cy_i24 = thread_arg->cy_i24;
			cy_r25 = thread_arg->cy_r25;		cy_i25 = thread_arg->cy_i25;
			cy_r26 = thread_arg->cy_r26;		cy_i26 = thread_arg->cy_i26;
			cy_r27 = thread_arg->cy_r27;		cy_i27 = thread_arg->cy_i27;

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
		#ifdef USE_AVX

			thread_arg->cy_r00 = cy_r00->d0;
			thread_arg->cy_r01 = cy_r00->d1;
			thread_arg->cy_r02 = cy_r00->d2;
			thread_arg->cy_r03 = cy_r00->d3;
			thread_arg->cy_r04 = cy_r04->d0;
			thread_arg->cy_r05 = cy_r04->d1;
			thread_arg->cy_r06 = cy_r04->d2;
			thread_arg->cy_r07 = cy_r04->d3;
			thread_arg->cy_r08 = cy_r08->d0;
			thread_arg->cy_r09 = cy_r08->d1;
			thread_arg->cy_r10 = cy_r08->d2;
			thread_arg->cy_r11 = cy_r08->d3;
			thread_arg->cy_r12 = cy_r12->d0;
			thread_arg->cy_r13 = cy_r12->d1;
			thread_arg->cy_r14 = cy_r12->d2;
			thread_arg->cy_r15 = cy_r12->d3;
			thread_arg->cy_r16 = cy_r16->d0;
			thread_arg->cy_r17 = cy_r16->d1;
			thread_arg->cy_r18 = cy_r16->d2;
			thread_arg->cy_r19 = cy_r16->d3;
			thread_arg->cy_r20 = cy_r20->d0;
			thread_arg->cy_r21 = cy_r20->d1;
			thread_arg->cy_r22 = cy_r20->d2;
			thread_arg->cy_r23 = cy_r20->d3;
			thread_arg->cy_r24 = cy_r24->d0;
			thread_arg->cy_r25 = cy_r24->d1;
			thread_arg->cy_r26 = cy_r24->d2;
			thread_arg->cy_r27 = cy_r24->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );

		#elif defined(USE_SSE2)

			thread_arg->cy_r00 = cy_r00->d0;
			thread_arg->cy_r01 = cy_r00->d1;
			thread_arg->cy_r02 = cy_r02->d0;
			thread_arg->cy_r03 = cy_r02->d1;
			thread_arg->cy_r04 = cy_r04->d0;
			thread_arg->cy_r05 = cy_r04->d1;
			thread_arg->cy_r06 = cy_r06->d0;
			thread_arg->cy_r07 = cy_r06->d1;
			thread_arg->cy_r08 = cy_r08->d0;
			thread_arg->cy_r09 = cy_r08->d1;
			thread_arg->cy_r10 = cy_r10->d0;
			thread_arg->cy_r11 = cy_r10->d1;
			thread_arg->cy_r12 = cy_r12->d0;
			thread_arg->cy_r13 = cy_r12->d1;
			thread_arg->cy_r14 = cy_r14->d0;
			thread_arg->cy_r15 = cy_r14->d1;
			thread_arg->cy_r16 = cy_r16->d0;
			thread_arg->cy_r17 = cy_r16->d1;
			thread_arg->cy_r18 = cy_r18->d0;
			thread_arg->cy_r19 = cy_r18->d1;
			thread_arg->cy_r20 = cy_r20->d0;
			thread_arg->cy_r21 = cy_r20->d1;
			thread_arg->cy_r22 = cy_r22->d0;
			thread_arg->cy_r23 = cy_r22->d1;
			thread_arg->cy_r24 = cy_r24->d0;
			thread_arg->cy_r25 = cy_r24->d1;
			thread_arg->cy_r26 = cy_r26->d0;
			thread_arg->cy_r27 = cy_r26->d1;
			maxerr = MAX(max_err->d0,max_err->d1);

		#else

			thread_arg->cy_r00 = cy_r00;
			thread_arg->cy_r01 = cy_r01;
			thread_arg->cy_r02 = cy_r02;
			thread_arg->cy_r03 = cy_r03;
			thread_arg->cy_r04 = cy_r04;
			thread_arg->cy_r05 = cy_r05;
			thread_arg->cy_r06 = cy_r06;
			thread_arg->cy_r07 = cy_r07;
			thread_arg->cy_r08 = cy_r08;
			thread_arg->cy_r09 = cy_r09;
			thread_arg->cy_r10 = cy_r10;
			thread_arg->cy_r11 = cy_r11;
			thread_arg->cy_r12 = cy_r12;
			thread_arg->cy_r13 = cy_r13;
			thread_arg->cy_r14 = cy_r14;
			thread_arg->cy_r15 = cy_r15;
			thread_arg->cy_r16 = cy_r16;
			thread_arg->cy_r17 = cy_r17;
			thread_arg->cy_r18 = cy_r18;
			thread_arg->cy_r19 = cy_r19;
			thread_arg->cy_r20 = cy_r20;
			thread_arg->cy_r21 = cy_r21;
			thread_arg->cy_r22 = cy_r22;
			thread_arg->cy_r23 = cy_r23;
			thread_arg->cy_r24 = cy_r24;
			thread_arg->cy_r25 = cy_r25;
			thread_arg->cy_r26 = cy_r26;
			thread_arg->cy_r27 = cy_r27;

		#endif	// AVX/SSE2
		}
		else
		{
		#ifdef USE_AVX

			thread_arg->cy_r00 = cy_r00->d0;	thread_arg->cy_i00 = cy_i00->d0;
			thread_arg->cy_r01 = cy_r00->d1;	thread_arg->cy_i01 = cy_i00->d1;
			thread_arg->cy_r02 = cy_r00->d2;	thread_arg->cy_i02 = cy_i00->d2;
			thread_arg->cy_r03 = cy_r00->d3;	thread_arg->cy_i03 = cy_i00->d3;
			thread_arg->cy_r04 = cy_r04->d0;	thread_arg->cy_i04 = cy_i04->d0;
			thread_arg->cy_r05 = cy_r04->d1;	thread_arg->cy_i05 = cy_i04->d1;
			thread_arg->cy_r06 = cy_r04->d2;	thread_arg->cy_i06 = cy_i04->d2;
			thread_arg->cy_r07 = cy_r04->d3;	thread_arg->cy_i07 = cy_i04->d3;
			thread_arg->cy_r08 = cy_r08->d0;	thread_arg->cy_i08 = cy_i08->d0;
			thread_arg->cy_r09 = cy_r08->d1;	thread_arg->cy_i09 = cy_i08->d1;
			thread_arg->cy_r10 = cy_r08->d2;	thread_arg->cy_i10 = cy_i08->d2;
			thread_arg->cy_r11 = cy_r08->d3;	thread_arg->cy_i11 = cy_i08->d3;
			thread_arg->cy_r12 = cy_r12->d0;	thread_arg->cy_i12 = cy_i12->d0;
			thread_arg->cy_r13 = cy_r12->d1;	thread_arg->cy_i13 = cy_i12->d1;
			thread_arg->cy_r14 = cy_r12->d2;	thread_arg->cy_i14 = cy_i12->d2;
			thread_arg->cy_r15 = cy_r12->d3;	thread_arg->cy_i15 = cy_i12->d3;
			thread_arg->cy_r16 = cy_r16->d0;	thread_arg->cy_i16 = cy_i16->d0;
			thread_arg->cy_r17 = cy_r16->d1;	thread_arg->cy_i17 = cy_i16->d1;
			thread_arg->cy_r18 = cy_r16->d2;	thread_arg->cy_i18 = cy_i16->d2;
			thread_arg->cy_r19 = cy_r16->d3;	thread_arg->cy_i19 = cy_i16->d3;
			thread_arg->cy_r20 = cy_r20->d0;	thread_arg->cy_i20 = cy_i20->d0;
			thread_arg->cy_r21 = cy_r20->d1;	thread_arg->cy_i21 = cy_i20->d1;
			thread_arg->cy_r22 = cy_r20->d2;	thread_arg->cy_i22 = cy_i20->d2;
			thread_arg->cy_r23 = cy_r20->d3;	thread_arg->cy_i23 = cy_i20->d3;
			thread_arg->cy_r24 = cy_r24->d0;	thread_arg->cy_i24 = cy_i24->d0;
			thread_arg->cy_r25 = cy_r24->d1;	thread_arg->cy_i25 = cy_i24->d1;
			thread_arg->cy_r26 = cy_r24->d2;	thread_arg->cy_i26 = cy_i24->d2;
			thread_arg->cy_r27 = cy_r24->d3;	thread_arg->cy_i27 = cy_i24->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );

		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			thread_arg->cy_r00 = cy_r00->d0;	thread_arg->cy_i00 = cy_r00->d1;
			thread_arg->cy_r01 = cy_r02->d0;	thread_arg->cy_i01 = cy_r02->d1;
			thread_arg->cy_r02 = cy_r04->d0;	thread_arg->cy_i02 = cy_r04->d1;
			thread_arg->cy_r03 = cy_r06->d0;	thread_arg->cy_i03 = cy_r06->d1;
			thread_arg->cy_r04 = cy_r08->d0;	thread_arg->cy_i04 = cy_r08->d1;
			thread_arg->cy_r05 = cy_r10->d0;	thread_arg->cy_i05 = cy_r10->d1;
			thread_arg->cy_r06 = cy_r12->d0;	thread_arg->cy_i06 = cy_r12->d1;
			thread_arg->cy_r07 = cy_r14->d0;	thread_arg->cy_i07 = cy_r14->d1;
			thread_arg->cy_r08 = cy_r16->d0;	thread_arg->cy_i08 = cy_r16->d1;
			thread_arg->cy_r09 = cy_r18->d0;	thread_arg->cy_i09 = cy_r18->d1;
			thread_arg->cy_r10 = cy_r20->d0;	thread_arg->cy_i10 = cy_r20->d1;
			thread_arg->cy_r11 = cy_r22->d0;	thread_arg->cy_i11 = cy_r22->d1;
			thread_arg->cy_r12 = cy_r24->d0;	thread_arg->cy_i12 = cy_r24->d1;
			thread_arg->cy_r13 = cy_r26->d0;	thread_arg->cy_i13 = cy_r26->d1;
			thread_arg->cy_r14 = cy_i00->d0;	thread_arg->cy_i14 = cy_i00->d1;
			thread_arg->cy_r15 = cy_i02->d0;	thread_arg->cy_i15 = cy_i02->d1;
			thread_arg->cy_r16 = cy_i04->d0;	thread_arg->cy_i16 = cy_i04->d1;
			thread_arg->cy_r17 = cy_i06->d0;	thread_arg->cy_i17 = cy_i06->d1;
			thread_arg->cy_r18 = cy_i08->d0;	thread_arg->cy_i18 = cy_i08->d1;
			thread_arg->cy_r19 = cy_i10->d0;	thread_arg->cy_i19 = cy_i10->d1;
			thread_arg->cy_r20 = cy_i12->d0;	thread_arg->cy_i20 = cy_i12->d1;
			thread_arg->cy_r21 = cy_i14->d0;	thread_arg->cy_i21 = cy_i14->d1;
			thread_arg->cy_r22 = cy_i16->d0;	thread_arg->cy_i22 = cy_i16->d1;
			thread_arg->cy_r23 = cy_i18->d0;	thread_arg->cy_i23 = cy_i18->d1;
			thread_arg->cy_r24 = cy_i20->d0;	thread_arg->cy_i24 = cy_i20->d1;
			thread_arg->cy_r25 = cy_i22->d0;	thread_arg->cy_i25 = cy_i22->d1;
			thread_arg->cy_r26 = cy_i24->d0;	thread_arg->cy_i26 = cy_i24->d1;
			thread_arg->cy_r27 = cy_i26->d0;	thread_arg->cy_i27 = cy_i26->d1;
			maxerr = MAX(max_err->d0,max_err->d1);

		#else

			thread_arg->cy_r00 = cy_r00;		thread_arg->cy_i00 = cy_i00;
			thread_arg->cy_r01 = cy_r01;		thread_arg->cy_i01 = cy_i01;
			thread_arg->cy_r02 = cy_r02;		thread_arg->cy_i02 = cy_i02;
			thread_arg->cy_r03 = cy_r03;		thread_arg->cy_i03 = cy_i03;
			thread_arg->cy_r04 = cy_r04;		thread_arg->cy_i04 = cy_i04;
			thread_arg->cy_r05 = cy_r05;		thread_arg->cy_i05 = cy_i05;
			thread_arg->cy_r06 = cy_r06;		thread_arg->cy_i06 = cy_i06;
			thread_arg->cy_r07 = cy_r07;		thread_arg->cy_i07 = cy_i07;
			thread_arg->cy_r08 = cy_r08;		thread_arg->cy_i08 = cy_i08;
			thread_arg->cy_r09 = cy_r09;		thread_arg->cy_i09 = cy_i09;
			thread_arg->cy_r10 = cy_r10;		thread_arg->cy_i10 = cy_i10;
			thread_arg->cy_r11 = cy_r11;		thread_arg->cy_i11 = cy_i11;
			thread_arg->cy_r12 = cy_r12;		thread_arg->cy_i12 = cy_i12;
			thread_arg->cy_r13 = cy_r13;		thread_arg->cy_i13 = cy_i13;
			thread_arg->cy_r14 = cy_r14;		thread_arg->cy_i14 = cy_i14;
			thread_arg->cy_r15 = cy_r15;		thread_arg->cy_i15 = cy_i15;
			thread_arg->cy_r16 = cy_r16;		thread_arg->cy_i16 = cy_i16;
			thread_arg->cy_r17 = cy_r17;		thread_arg->cy_i17 = cy_i17;
			thread_arg->cy_r18 = cy_r18;		thread_arg->cy_i18 = cy_i18;
			thread_arg->cy_r19 = cy_r19;		thread_arg->cy_i19 = cy_i19;
			thread_arg->cy_r20 = cy_r20;		thread_arg->cy_i20 = cy_i20;
			thread_arg->cy_r21 = cy_r21;		thread_arg->cy_i21 = cy_i21;
			thread_arg->cy_r22 = cy_r22;		thread_arg->cy_i22 = cy_i22;
			thread_arg->cy_r23 = cy_r23;		thread_arg->cy_i23 = cy_i23;
			thread_arg->cy_r24 = cy_r24;		thread_arg->cy_i24 = cy_i24;
			thread_arg->cy_r25 = cy_r25;		thread_arg->cy_i25 = cy_i25;
			thread_arg->cy_r26 = cy_r26;		thread_arg->cy_i26 = cy_i26;
			thread_arg->cy_r27 = cy_r27;		thread_arg->cy_i27 = cy_i27;

		#endif	// AVX/SSE2
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
