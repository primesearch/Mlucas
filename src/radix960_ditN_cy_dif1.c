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

#define RADIX 960	// Use #define rather than const int to ensure it's really a compile-time const in the C sense
#define ODD_RADIX 15	// ODD_RADIX = [radix >> trailz(radix)]

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

// See the radix28 version of this routine for details about the
// small-cyclic-array indexing scheme used in the fermat_carry_norm_errcheckB macros.

#ifdef USE_SSE2

  // RADIX = 0x3c0
  // For Mersenne-mod need (16 [SSE2] or 64 [AVX]) + (4 [HIACC] or 40 [LOACC]) added slots for half_arr lookup tables.
  // Max = (40 [SSE2]; 132 [AVX]),
  // For Fermat-mod in AVX mode we need RADIX*4 = 0x3c0*4 = 0xf00 [if HIACC] or 12 [if not] vec_dbl slots for the compact
  // negacyclic-roots chained-multiply scheme.
  // Add larger number in each case - i.e. max(132,3840) = 3840 if AVX, 40 if SSE2
  // to (half_arr_offset + RADIX) to get required value of radix960_creals_in_local_store:
  #ifdef USE_AVX512	// 0xf0 fewer carry slots than AVX:
	const int half_arr_offset960 = 0x1038;	// + RADIX = 0x13f8;  Used for thread local-storage-integrity checking
	const int radix960_creals_in_local_store = 0x2300;	// AVX-512: 0x13f8 + 0xf00, round up to next multiple of 0x10
  #elif defined(USE_AVX)
	const int half_arr_offset960 = 0x1128;	// + RADIX = 0x14e8;  Used for thread local-storage-integrity checking
	const int radix960_creals_in_local_store = 0x23f0;	// AVX: 0x14e8 + 0xf00, round up to next multiple of 0x10
  #else
	const int half_arr_offset960 = 0x1308;	// + RADIX = 0x16c8; Used for thread local-storage-integrity checking
	const int radix960_creals_in_local_store = 0x16f4;	// (half_arr_offset + RADIX) + 0x28 (=4), round up to next multiple of 0x10
  #endif

  #ifdef USE_AVX
	#include "radix960_avx_negadwt_consts.h"
  #endif
	#include "sse2_macro_gcc64.h"
	#include "radix15_sse_macro.h"

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

int radix960_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-960 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-960 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix960_ditN_cy_dif1";
	static int thr_id = 0;	// Master thread gets this special id
	const int pfetch_dist = PFETCH_DIST;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
	static double wts_mult[2], inv_mult[2];	// Const wts-multiplier and 2*(its multiplicative inverse)
	double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
	// Cleanup loop assumes carryins propagate at most 4 words up, but need at least 1 vec_cmplx
	// (2 vec_dbl)'s worth of doubles in wraparound step, hence AVX-512 needs mers-value bumped up:
  #ifdef USE_AVX512
	const int jhi_wrap_mers = 15;
	const int jhi_wrap_ferm = 15;
  #else
	const int jhi_wrap_mers =  7;
	const int jhi_wrap_ferm = 15;	// For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
  #endif
	int NDIVR,i,incr = 0,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,ntmp,outer,nbytes;
	// incr = Carry-chain wts-multipliers recurrence length, which must divide
	// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 60|120|240 for avx512,avx,sse, respectively:
	const int incr_long = 15,incr_med =10,incr_short = 5;
  // Have no specialized HIACC carry macro in USE_AVX512 and ARMv8 SIMD, use nonzero incr-value to differenetiate vs AVX/AVX2:
  #if defined(USE_AVX512) || defined(USE_ARM_V8_SIMD)
	const int incr_hiacc = 3;
  #else
	const int incr_hiacc = 0;
  #endif
	// Fermat-mod: For AVX+, define carry-subchain length in terms of 2^nfold subchains:
	const int *inc_arr = 0x0;
  #ifdef USE_AVX512
	// For nfold > 2,  RADIX/8 not divisible by 2^nfold, so use a more-general inner-loop scheme which can handle that:
	int nfold = USE_SHORT_CY_CHAIN + 2;
	const int nexec_long[] = {30,30,30,30}, nexec_med[] = {15,15,15,15,15,15,15,15}, nexec_short[] = {8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7}, nexec_hiacc[] = {4,4,4,3,4,4,4,3,4,4,4,3,4,4,4,3,4,4,4,3,4,4,4,3,4,4,4,3,4,4,4,3};
  #elif defined(USE_AVX)
	// For nfold > 4,  RADIX/4 not divisible by 2^nfold, so use a more-general inner-loop scheme which can handle that:
	int nfold = USE_SHORT_CY_CHAIN + 2;
	const int nexec_long[] = {60,60,60,60}, nexec_med[] = {30,30,30,30,30,30,30,30}, nexec_short[] = {15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15}, nexec_hiacc[] = {8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7};
  #endif
	// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
		if(USE_SHORT_CY_CHAIN == 0)
			incr = incr_long;
		else if(USE_SHORT_CY_CHAIN == 1)
			incr = incr_med;
		else if(USE_SHORT_CY_CHAIN == 2)
			incr = incr_short;
		else
			incr = incr_hiacc;
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
	static uint32 bw,sw,bjmodnini,
		p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
		p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,
		p100,p110,p120,p130,p140,p150,p160,p170,p180,p190,p1a0,p1b0,p1c0,p1d0,p1e0,p1f0,
		p200,p210,p220,p230,p240,p250,p260,p270,p280,p290,p2a0,p2b0,p2c0,p2d0,p2e0,p2f0,
		p300,p310,p320,p330,p340,p350,p360,p370,p380,p390,p3a0,p3b0, nsave = 0;
	static int poff[RADIX>>2];
#ifndef MULTITHREAD
	int kk,mask3, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf;
	static int plo[16], phi[RADIX>>4], p_out_hi[RADIX>>4];
	static int dif_o_offsets[64], dit_i_offsets[64], t_offsets[64];
	uint64 i64;
// DIF:
	const uint64 dif_perm16[60] = {
		0x0123456789abcdefull,0x1032547698badcfeull,0x23106754ab98efdcull,0x32017645ba89fecdull,
		0x54763201dcfeba89ull,0xcdefab9823106754ull,0x76450123fecd89abull,0xefdc98ba10325476ull,
		0xab98efdc67541032ull,0xba89fecd76450123ull,0x98badcfe54763201ull,0x45672310cdefab98ull,
		0xdcfeba8932017645ull,0x23106754ab98efdcull,0xfecd89ab01234567ull,0x1032547698badcfeull,
		0x67541032efdc98baull,0x76450123fecd89abull,0x54763201dcfeba89ull,0xcdefab9823106754ull,
		0x89abcdef45672310ull,0x98badcfe54763201ull,0xab98efdc67541032ull,0xba89fecd76450123ull,
		0xefdc98ba10325476ull,0xfecd89ab01234567ull,0xdcfeba8932017645ull,0x23106754ab98efdcull,
		0x45672310cdefab98ull,0x54763201dcfeba89ull,0x67541032efdc98baull,0x76450123fecd89abull,
		0x32017645ba89fecdull,0xab98efdc67541032ull,0x89abcdef45672310ull,0x98badcfe54763201ull,
		0xcdefab9823106754ull,0xdcfeba8932017645ull,0xefdc98ba10325476ull,0xfecd89ab01234567ull,
		0xba89fecd76450123ull,0x67541032efdc98baull,0x45672310cdefab98ull,0x54763201dcfeba89ull,
		0x1032547698badcfeull,0x89abcdef45672310ull,0x32017645ba89fecdull,0xab98efdc67541032ull,
		0x76450123fecd89abull,0xefdc98ba10325476ull,0xcdefab9823106754ull,0xdcfeba8932017645ull,
		0x98badcfe54763201ull,0x45672310cdefab98ull,0xba89fecd76450123ull,0x67541032efdc98baull,
		0x23106754ab98efdcull,0x32017645ba89fecdull,0x1032547698badcfeull,0x89abcdef45672310ull
	};
// DIT:
	// Low parts [p0-f] of input-index perms - 16 distinct ones, 4 lookups for each radix-64 DFT:
	const uint64 dit_perm16[16] = {
		0x01327654fedcba98ull,0xfedcba9876543210ull,0x5467102398abefcdull,0x98abefcd10236745ull,
		0xab89cdfe23014576ull,0x23014576cdfe89baull,0x76543210ba98dcefull,0x10236745efcdab89ull,
		0xcdfe89ba45760132ull,0xba98dcef32105467ull,0xefcdab8967452301ull,0x4576013289bafedcull,
		0x32105467dcef98abull,0x67452301ab89cdfeull,0x89bafedc01327654ull,0xdcef98ab54671023ull
	};
	// Pattern of access of above-encoded 16-perms, 4 for each DFT-64: Tempting to encoded via hex-hars,
	// but not worth it, just use byte array and get half the maximum possible compression ratio but simpler code:
	const uint8 dit_perm16_idx[60] = {
		0x0,0x1,0x1,0x1,0x2,0x2,0x2,0x3,0x4,0x5,0x4,0x5,0x6,0x6,0x6,
		0x6,0x3,0x3,0x3,0x7,0x5,0x8,0x5,0x8,0x9,0x9,0x9,0x9,0x7,0x7,
		0x7,0xa,0x8,0xb,0xb,0xb,0xc,0xc,0xc,0xc,0xa,0xd,0xa,0xd,0xb,
		0xe,0xe,0xe,0xf,0xf,0xf,0x2,0xd,0x4,0xd,0x4,0xe,0x0,0x0,0x0
	};
	// Pattern of access of phi[] elements for the above-encoded 16-perms, 4 for each DFT-64:
	const uint8 dit_perm16_phidx[60] = {
		0x00,0x01,0x03,0x02,0x09,0x08,0x0a,0x0b,0x06,0x07,0x04,0x05,0x17,0x16,0x15,
		0x14,0x11,0x10,0x12,0x13,0x0e,0x0f,0x0c,0x0d,0x1f,0x1e,0x1d,0x1c,0x19,0x18,
		0x1a,0x1b,0x20,0x21,0x23,0x22,0x27,0x26,0x25,0x24,0x2e,0x2f,0x2c,0x2d,0x28,
		0x29,0x2b,0x2a,0x39,0x38,0x3a,0x3b,0x36,0x37,0x34,0x35,0x30,0x31,0x33,0x32
	};
#endif
	static double radix_inv, n2inv;
	const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/
					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double scale,dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	double *add0, *addr, *addi;
	struct complex t[RADIX], *tptr;
	int *itmp,*itm2;	// Pointer into the bjmodn array
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	static int n_div_nwt;
	int col,co2,co3,m,m2;
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
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
  #endif
	double rt,it;
	// indices into weights arrays (mod NWT):
	static int ii[ODD_RADIX] = {
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1
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
	double *add1,*add2,*add3;
  #endif	// MULTITHREAD

	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
	const double crnd = 3.0*0x4000000*0x2000000;
	static vec_dbl *max_err, *sse2_rnd, *half_arr, *isrt2,*one,*two,
		*sse2_c3m1, *sse2_s, *sse2_cn1, *sse2_cn2, *sse2_ss3, *sse2_sn1, *sse2_sn2,	// Need radix-15 roots explicitly
		*r00,	// Head of RADIX*vec_cmplx-sized local store #1
		*s1p00,	// Head of RADIX*vec_cmplx-sized local store #2
		*x00,*x01,*x02,*x03,*x04,*x05,*x06,*x07,*x08,*x09,*x0a,*x0b,*x0c,*x0d,*x0e,	// Temps for 2x radix-15 DFT
		*y00,*y01,*y02,*y03,*y04,*y05,*y06,*y07,*y08,*y09,*y0a,*y0b,*y0c,*y0d,*y0e,
		*cy_r,*cy_i;	// Need RADIX slots for sse2 carries, RADIX/2 for avx
  #ifdef USE_AVX
	static vec_dbl *base_negacyclic_root;
  #endif

	vec_dbl *tmp,*tm0,*tm1,*tm2;	// Non-static utility ptrs
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer

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
	static task_control_t   task_control = {NULL, (void*)cy960_process_chunk, NULL, 0x0};

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
		cslots_in_local_store = radix960_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix960_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = sc_ptr;
		r00   = tmp;	tmp += 0x780;	// Head of RADIX*vec_cmplx-sized local store #1
		s1p00 = tmp;	tmp += 0x780;	// Head of RADIX*vec_cmplx-sized local store #2
		x00    = tmp + 0x00;
		x01    = tmp + 0x02;
		x02    = tmp + 0x04;
		x03    = tmp + 0x06;
		x04    = tmp + 0x08;
		x05    = tmp + 0x0a;
		x06    = tmp + 0x0c;
		x07    = tmp + 0x0e;
		x08    = tmp + 0x10;
		x09    = tmp + 0x12;
		x0a    = tmp + 0x14;
		x0b    = tmp + 0x16;
		x0c    = tmp + 0x18;
		x0d    = tmp + 0x1a;
		x0e    = tmp + 0x1c;
		tmp += 0x1e;
		y00    = tmp + 0x00;
		y01    = tmp + 0x02;
		y02    = tmp + 0x04;
		y03    = tmp + 0x06;
		y04    = tmp + 0x08;
		y05    = tmp + 0x0a;
		y06    = tmp + 0x0c;
		y07    = tmp + 0x0e;
		y08    = tmp + 0x10;
		y09    = tmp + 0x12;
		y0a    = tmp + 0x14;
		y0b    = tmp + 0x16;
		y0c    = tmp + 0x18;
		y0d    = tmp + 0x1a;
		y0e    = tmp + 0x1c;
		tmp += 0x1e;	// += 0x3c => sc_ptr + 0xf3c
		// DFT-15 roots needed explicitly, but throw in isrt2 as a pad-to-een and handy anchoring element:
		isrt2     = tmp + 0x00;
		sse2_c3m1 = tmp + 0x01;
		sse2_s    = tmp + 0x02;
		sse2_cn1  = tmp + 0x03;
		sse2_cn2  = tmp + 0x04;
		sse2_ss3  = tmp + 0x05;
		sse2_sn1  = tmp + 0x06;
		sse2_sn2  = tmp + 0x07;
		one       = tmp + 0x08;
		two       = tmp + 0x09;
		tmp += 0x0a;	// += 0xa => sc_ptr + 0xf46
	  #ifdef USE_AVX512
		cy_r = tmp;	cy_i = tmp+0x78;	tmp += 2*0x78;	// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
		cy_r = tmp;	cy_i = tmp+0x0f0;	tmp += 2*0x0f0;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0x1e0 + 2 => sc_ptr += 0x1128
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 96 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy_r = tmp;	cy_i = tmp+0x1e0;	tmp += 2*0x1e0;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0x3c0 + 2 => sc_ptr += 0x1308
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 32 x 16 bytes for Mersenne-mod, and [4*ODD_RADIX] x 16 for Fermat-mod */
	  #endif
		ASSERT(half_arr_offset960 == (uint32)(half_arr-sc_ptr), "half_arr_offset mismatches actual!");
		ASSERT((radix960_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r00) + (20 << L2_SZ_VD), "radix960_creals_in_local_store checksum failed!");

		/* These remain fixed: */
		VEC_DBL_INIT(isrt2,ISRT2);
		VEC_DBL_INIT(sse2_c3m1, c3m1);
		VEC_DBL_INIT(sse2_s   , s   );
		VEC_DBL_INIT(sse2_cn1 , cn1 );
		VEC_DBL_INIT(sse2_cn2 , cn2 );
		VEC_DBL_INIT(sse2_ss3 , ss3 );
		VEC_DBL_INIT(sse2_sn1 , sn1 );
		VEC_DBL_INIT(sse2_sn2 , sn2 );
		VEC_DBL_INIT(one  ,  1.0);
		VEC_DBL_INIT(two  ,  2.0);
		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		sse2_rnd->d0 = base[0]; sse2_rnd->d1 = baseinv[1]; sse2_rnd->d2 = wts_mult[1]; sse2_rnd->d3 = inv_mult[0];
	  #else
		VEC_DBL_INIT(sse2_rnd, crnd);
	  #endif

	// Init-mode calls to these functions which maintain an internal local-alloc static store:
		thr_id = -1;	// Use this special thread id for any macro-required thread-local-data inits...
		SSE2_RADIX_64_DIF(CY_THREADS,thr_id, 0,0,0,0,0,0);
		SSE2_RADIX_64_DIT(CY_THREADS,thr_id, 0,0,0,0,0);
		thr_id = 0;	// ...then revert to 0.

		// Propagate the above consts to the remaining threads:
		nbytes = (intptr_t)cy_r - (intptr_t)isrt2;	// #bytes in 1st of above block of consts
		tmp = isrt2;
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
		// Simple qfloat-based loop to crank out the roots which populate the radix960_avx_negadwt_consts table:
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
		tmp64 = radix960_avx_negadwt_consts[1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(  1*I*Pi/(2*RADIX)) = sin((RADIX-  1)*I*Pi/(2*RADIX)) */
		tmp64 = radix960_avx_negadwt_consts[2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(  2*I*Pi/(2*RADIX)) = sin((RADIX-  2)*I*Pi/(2*RADIX)) */
		tmp64 = radix960_avx_negadwt_consts[3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(  3*I*Pi/(2*RADIX)) = sin((RADIX-  3)*I*Pi/(2*RADIX)) */	tmp += 2;
		for(j = 4; j < RADIX; j += 4) {
			tmp64 = radix960_avx_negadwt_consts[j+0];	tmp->d0 = tm2->d0 = *(double *)&tmp64;	tm2 -= 2;
			tmp64 = radix960_avx_negadwt_consts[j+1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;
			tmp64 = radix960_avx_negadwt_consts[j+2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;
			tmp64 = radix960_avx_negadwt_consts[j+3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	tmp += 2;
		}
		tmp = base_negacyclic_root + RADIX*2;	// reset to point to start of above block
		nbytes = RADIX*SZ_VD/2;	// RADIX/4 AVX-register-sized complex data

	  #else	// HIACC = false:

		// lower-precision version, which yields slightly more roundoff error, but is simpler and more storage-compact.
	   #ifdef USE_AVX512

		// Init exp(j*I*Pi/2/RADIX), for j = 0-7:
		tmp = base_negacyclic_root + 16;	// First 16 slots reserved for Re/Im parts of the 8 base multipliers
		tmp->d0 = 1.0;
		tmp64 = radix960_avx_negadwt_consts[      1];	tmp->d1 = *(double *)&tmp64;	// cos(01*I*Pi/(2*RADIX))
		tmp64 = radix960_avx_negadwt_consts[      2];	tmp->d2 = *(double *)&tmp64;	// cos(02*I*Pi/(2*RADIX))
		tmp64 = radix960_avx_negadwt_consts[      3];	tmp->d3 = *(double *)&tmp64;	// cos(03*I*Pi/(2*RADIX))
		tmp64 = radix960_avx_negadwt_consts[      4];	tmp->d4 = *(double *)&tmp64;	// cos(04*I*Pi/(2*RADIX))
		tmp64 = radix960_avx_negadwt_consts[      5];	tmp->d5 = *(double *)&tmp64;	// cos(05*I*Pi/(2*RADIX))
		tmp64 = radix960_avx_negadwt_consts[      6];	tmp->d6 = *(double *)&tmp64;	// cos(06*I*Pi/(2*RADIX))
		tmp64 = radix960_avx_negadwt_consts[      7];	tmp->d7 = *(double *)&tmp64;	// cos(07*I*Pi/(2*RADIX))
		++tmp;
		tmp->d0 = 0.0;
		tmp64 = radix960_avx_negadwt_consts[RADIX-1];	tmp->d1 = *(double *)&tmp64;	// sin(01*I*Pi/(2*RADIX))
		tmp64 = radix960_avx_negadwt_consts[RADIX-2];	tmp->d2 = *(double *)&tmp64;	// sin(02*I*Pi/(2*RADIX))
		tmp64 = radix960_avx_negadwt_consts[RADIX-3];	tmp->d3 = *(double *)&tmp64;	// sin(03*I*Pi/(2*RADIX))
		tmp64 = radix960_avx_negadwt_consts[RADIX-4];	tmp->d4 = *(double *)&tmp64;	// sin(04*I*Pi/(2*RADIX))
		tmp64 = radix960_avx_negadwt_consts[RADIX-5];	tmp->d5 = *(double *)&tmp64;	// sin(05*I*Pi/(2*RADIX))
		tmp64 = radix960_avx_negadwt_consts[RADIX-6];	tmp->d6 = *(double *)&tmp64;	// sin(06*I*Pi/(2*RADIX))
		tmp64 = radix960_avx_negadwt_consts[RADIX-7];	tmp->d7 = *(double *)&tmp64;	// sin(07*I*Pi/(2*RADIX))
		++tmp;	// 0x480(base_negacyclic_root)
		tmp64 = radix960_avx_negadwt_consts[      8];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// cos(08*I*Pi/(2*RADIX))
		++tmp;	// 0x4c0(base_negacyclic_root)
		tmp64 = radix960_avx_negadwt_consts[RADIX-8];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// sin(08*I*Pi/(2*RADIX))
		tmp = base_negacyclic_root + 16;	// reset to point to start of above block

	   #elif defined(USE_AVX)

		// Init exp(j*I*Pi/2/RADIX), for j = 0-3:
		tmp = base_negacyclic_root + 8;	// First 8 slots reserved for Re/Im parts of the 4 base multipliers
		tmp->d0 = 1.0;
		tmp64 = radix960_avx_negadwt_consts[      1];	tmp->d1 = *(double *)&tmp64;	// cos(01*I*Pi/(2*RADIX))
		tmp64 = radix960_avx_negadwt_consts[      2];	tmp->d2 = *(double *)&tmp64;	// cos(02*I*Pi/(2*RADIX))
		tmp64 = radix960_avx_negadwt_consts[      3];	tmp->d3 = *(double *)&tmp64;	// cos(03*I*Pi/(2*RADIX))

		(++tmp)->d0 = 0.0;
		tmp64 = radix960_avx_negadwt_consts[RADIX-1];	tmp->d1 = *(double *)&tmp64;	// sin(01*I*Pi/(2*RADIX))
		tmp64 = radix960_avx_negadwt_consts[RADIX-2];	tmp->d2 = *(double *)&tmp64;	// sin(02*I*Pi/(2*RADIX))
		tmp64 = radix960_avx_negadwt_consts[RADIX-3];	tmp->d3 = *(double *)&tmp64;	// sin(03*I*Pi/(2*RADIX))
		++tmp;
		tmp64 = radix960_avx_negadwt_consts[      4];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// cos(04*I*Pi/(2*RADIX))
		++tmp;
		tmp64 = radix960_avx_negadwt_consts[RADIX-4];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// sin(04*I*Pi/(2*RADIX))
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
		ASSERT(tmp == half_arr, "tmp == half_arr check failed!");
	#ifdef USE_AVX512
		/* Each lookup-category in the 'mini-tables' used in AVX mode balloons from 16x32-bytes to 64x64-bytes,
			so switch to an opmask-based scheme which starts with e.g. a broadcast constant and onditional doubling.
			Here are the needed consts and opmasks:
			[1] Fwd-wt multipliers: Init = 0.50 x 8, anytime AVX-style lookup into 1st table below would have bit = 0, double the corr. datum
			[2] Inv-wt multipliers: Init = 0.25 x 8, anytime AVX-style lookup into 2nd table below would have bit = 0, double the corr. datum
			[3] Fwd-base mults: Init = base[0] x 8, anytime AVX-style lookup into 3rd table below would have bit = 1, double the corr. datum
			[4] Inv-base mults: Init = binv[1] x 8, anytime AVX-style lookup into 4th table below would have bit = 0, double the corr. datum
			[5] [LOACC] Init = wts_mult[1] x 8, anytime AVX-style lookup into 5th table below would have bit = 0, double the corr. datum
			[6] [LOACC] Init = inv_mult[0] x 8, anytime AVX-style lookup into 6th table below would have bit = 1, double the corr. datum
		*/
		nbytes = 2 << L2_SZ_VD;
		// Dec 2021: k10m / IMCI512 makes it really hard to on-the-fly vector-int zmm = 0.5 x 8 w/o supplying a new input-arg
		// the asm carry-macros containing (pointer to mem64 containing 0.5), so add inits of v8_double(0.5) and v8_double(0.25):
		VEC_DBL_INIT(tmp  , 0.50);
		VEC_DBL_INIT(tmp+1, 0.25);
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
	nbytes += 64;
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
	NDIVR <<= 4; p10 = NDIVR;	pf += ( (pf >> DAT_BITS) << PAD_BITS );
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
		pd0 = pc0 + NDIVR;		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pe0 = pd0 + NDIVR;		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		pf0 = pe0 + NDIVR;		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		p100= pf0 + NDIVR;		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
		p110= p100 + NDIVR;		p100+= ( (p100>> DAT_BITS) << PAD_BITS );
		p120= p110 + NDIVR;		p110+= ( (p110>> DAT_BITS) << PAD_BITS );
		p130= p120 + NDIVR;		p120+= ( (p120>> DAT_BITS) << PAD_BITS );
		p140= p130 + NDIVR;		p130+= ( (p130>> DAT_BITS) << PAD_BITS );
		p150= p140 + NDIVR;		p140+= ( (p140>> DAT_BITS) << PAD_BITS );
		p160= p150 + NDIVR;		p150+= ( (p150>> DAT_BITS) << PAD_BITS );
		p170= p160 + NDIVR;		p160+= ( (p160>> DAT_BITS) << PAD_BITS );
		p180= p170 + NDIVR;		p170+= ( (p170>> DAT_BITS) << PAD_BITS );
		p190= p180 + NDIVR;		p180+= ( (p180>> DAT_BITS) << PAD_BITS );
		p1a0= p190 + NDIVR;		p190+= ( (p190>> DAT_BITS) << PAD_BITS );
		p1b0= p1a0 + NDIVR;		p1a0+= ( (p1a0>> DAT_BITS) << PAD_BITS );
		p1c0= p1b0 + NDIVR;		p1b0+= ( (p1b0>> DAT_BITS) << PAD_BITS );
		p1d0= p1c0 + NDIVR;		p1c0+= ( (p1c0>> DAT_BITS) << PAD_BITS );
		p1e0= p1d0 + NDIVR;		p1d0+= ( (p1d0>> DAT_BITS) << PAD_BITS );
		p1f0= p1e0 + NDIVR;		p1e0+= ( (p1e0>> DAT_BITS) << PAD_BITS );
		p200= p1f0 + NDIVR;		p1f0+= ( (p1f0>> DAT_BITS) << PAD_BITS );
		p210= p200 + NDIVR;		p200+= ( (p200>> DAT_BITS) << PAD_BITS );
		p220= p210 + NDIVR;		p210+= ( (p210>> DAT_BITS) << PAD_BITS );
		p230= p220 + NDIVR;		p220+= ( (p220>> DAT_BITS) << PAD_BITS );
		p240= p230 + NDIVR;		p230+= ( (p230>> DAT_BITS) << PAD_BITS );
		p250= p240 + NDIVR;		p240+= ( (p240>> DAT_BITS) << PAD_BITS );
		p260= p250 + NDIVR;		p250+= ( (p250>> DAT_BITS) << PAD_BITS );
		p270= p260 + NDIVR;		p260+= ( (p260>> DAT_BITS) << PAD_BITS );
		p280= p270 + NDIVR;		p270+= ( (p270>> DAT_BITS) << PAD_BITS );
		p290= p280 + NDIVR;		p280+= ( (p280>> DAT_BITS) << PAD_BITS );
		p2a0= p290 + NDIVR;		p290+= ( (p290>> DAT_BITS) << PAD_BITS );
		p2b0= p2a0 + NDIVR;		p2a0+= ( (p2a0>> DAT_BITS) << PAD_BITS );
		p2c0= p2b0 + NDIVR;		p2b0+= ( (p2b0>> DAT_BITS) << PAD_BITS );
		p2d0= p2c0 + NDIVR;		p2c0+= ( (p2c0>> DAT_BITS) << PAD_BITS );
		p2e0= p2d0 + NDIVR;		p2d0+= ( (p2d0>> DAT_BITS) << PAD_BITS );
		p2f0= p2e0 + NDIVR;		p2e0+= ( (p2e0>> DAT_BITS) << PAD_BITS );
		p300= p2f0 + NDIVR;		p2f0+= ( (p2f0>> DAT_BITS) << PAD_BITS );
		p310= p300 + NDIVR;		p300+= ( (p300>> DAT_BITS) << PAD_BITS );
		p320= p310 + NDIVR;		p310+= ( (p310>> DAT_BITS) << PAD_BITS );
		p330= p320 + NDIVR;		p320+= ( (p320>> DAT_BITS) << PAD_BITS );
		p340= p330 + NDIVR;		p330+= ( (p330>> DAT_BITS) << PAD_BITS );
		p350= p340 + NDIVR;		p340+= ( (p340>> DAT_BITS) << PAD_BITS );
		p360= p350 + NDIVR;		p350+= ( (p350>> DAT_BITS) << PAD_BITS );
		p370= p360 + NDIVR;		p360+= ( (p360>> DAT_BITS) << PAD_BITS );
		p380= p370 + NDIVR;		p370+= ( (p370>> DAT_BITS) << PAD_BITS );
		p390= p380 + NDIVR;		p380+= ( (p380>> DAT_BITS) << PAD_BITS );
		p3a0= p390 + NDIVR;		p390+= ( (p390>> DAT_BITS) << PAD_BITS );
		p3b0= p3a0 + NDIVR;		p3a0+= ( (p3a0>> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			p3b0+= ( (p3b0>> DAT_BITS) << PAD_BITS );
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
		for(l = 0; l < 64; l++) {
			poff[l+ 64] = poff[l] + p100;
			poff[l+128] = poff[l] + p200;
			if(l < 48) {
				poff[l+192] = poff[l] + p300;
			}
		}

	#ifndef MULTITHREAD
		plo[ 0]= 0;plo[ 1]=p1;plo[ 2]=p2;plo[ 3]=p3;plo[ 4]=p4;plo[ 5]=p5;plo[ 6]=p6;plo[ 7]=p7;
		plo[ 8]=p8;plo[ 9]=p9;plo[10]=pa;plo[11]=pb;plo[12]=pc;plo[13]=pd;plo[14]=pe;plo[15]=pf;

		phi[0x00] =    0; phi[0x01] =  p10; phi[0x02] =  p20; phi[0x03] =  p30;
		phi[0x04] =  p40; phi[0x05] =  p50; phi[0x06] =  p60; phi[0x07] =  p70;
		phi[0x08] =  p80; phi[0x09] =  p90; phi[0x0a] =  pa0; phi[0x0b] =  pb0;
		phi[0x0c] =  pc0; phi[0x0d] =  pd0; phi[0x0e] =  pe0; phi[0x0f] =  pf0;
		phi[0x10] = p100; phi[0x11] = p110; phi[0x12] = p120; phi[0x13] = p130;
		phi[0x14] = p140; phi[0x15] = p150; phi[0x16] = p160; phi[0x17] = p170;
		phi[0x18] = p180; phi[0x19] = p190; phi[0x1a] = p1a0; phi[0x1b] = p1b0;
		phi[0x1c] = p1c0; phi[0x1d] = p1d0; phi[0x1e] = p1e0; phi[0x1f] = p1f0;
		phi[0x20] = p200; phi[0x21] = p210; phi[0x22] = p220; phi[0x23] = p230;
		phi[0x24] = p240; phi[0x25] = p250; phi[0x26] = p260; phi[0x27] = p270;
		phi[0x28] = p280; phi[0x29] = p290; phi[0x2a] = p2a0; phi[0x2b] = p2b0;
		phi[0x2c] = p2c0; phi[0x2d] = p2d0; phi[0x2e] = p2e0; phi[0x2f] = p2f0;
		phi[0x30] = p300; phi[0x31] = p310; phi[0x32] = p320; phi[0x33] = p330;
		phi[0x34] = p340; phi[0x35] = p350; phi[0x36] = p360; phi[0x37] = p370;
		phi[0x38] = p380; phi[0x39] = p390; phi[0x3a] = p3a0; phi[0x3b] = p3b0;

		p_out_hi[0x00] =    0; p_out_hi[0x01] =  p10; p_out_hi[0x02] =  p20; p_out_hi[0x03] =  p30;	//  p00
		p_out_hi[0x04] =  p10; p_out_hi[0x05] =    0; p_out_hi[0x06] =  p30; p_out_hi[0x07] =  p20;	//  p80
		p_out_hi[0x08] =  p20; p_out_hi[0x09] =  p30; p_out_hi[0x0a] =  p10; p_out_hi[0x0b] =    0;	//  p40
		p_out_hi[0x0c] =  p10; p_out_hi[0x0d] =    0; p_out_hi[0x0e] =  p30; p_out_hi[0x0f] =  p20;	// p380
		p_out_hi[0x10] =  p20; p_out_hi[0x11] =  p30; p_out_hi[0x12] =  p10; p_out_hi[0x13] =    0;	// p340
		p_out_hi[0x14] =    0; p_out_hi[0x15] =  p10; p_out_hi[0x16] =  p20; p_out_hi[0x17] =  p30;	// p300
		p_out_hi[0x18] =  p20; p_out_hi[0x19] =  p30; p_out_hi[0x1a] =  p10; p_out_hi[0x1b] =    0;	// p2c0
		p_out_hi[0x1c] =    0; p_out_hi[0x1d] =  p10; p_out_hi[0x1e] =  p20; p_out_hi[0x1f] =  p30;	// p280
		p_out_hi[0x20] =  p30; p_out_hi[0x21] =  p20; p_out_hi[0x22] =    0; p_out_hi[0x23] =  p10;	// p240
		p_out_hi[0x24] =    0; p_out_hi[0x25] =  p10; p_out_hi[0x26] =  p20; p_out_hi[0x27] =  p30;	// p200
		p_out_hi[0x28] =  p30; p_out_hi[0x29] =  p20; p_out_hi[0x2a] =    0; p_out_hi[0x2b] =  p10;	// p1c0
		p_out_hi[0x2c] =  p10; p_out_hi[0x2d] =    0; p_out_hi[0x2e] =  p30; p_out_hi[0x2f] =  p20;	// p180
		p_out_hi[0x30] =  p30; p_out_hi[0x31] =  p20; p_out_hi[0x32] =    0; p_out_hi[0x33] =  p10;	// p140
		p_out_hi[0x34] =  p10; p_out_hi[0x35] =    0; p_out_hi[0x36] =  p30; p_out_hi[0x37] =  p20;	// p100
		p_out_hi[0x38] =  p20; p_out_hi[0x39] =  p30; p_out_hi[0x3a] =  p10; p_out_hi[0x3b] =    0;	//  pc0

		for(l = 0; l < 64; l++) {
			t_offsets[l] = ((l<<4) - l);	// = l*15
			t_offsets[l] <<= 1;	// Need 2x incr due to the Radix-64 DFT casting t-array pointer to double.
								// In SIMD mode similarly need 2x due to vec_dbl type being used to access vec-complex data
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

/*...The radix-960 final DIT pass is here.	*/

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
		_i[0] = 0;		/* Pointer to the BASE and BASEINV arrays. If n divides p, lowest-order digit is always a smallword (_i[0] = 0).	*/
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
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		ASSERT(((tmp-1)->d0 == base[0] && (tmp-1)->d1 == baseinv[1] && (tmp-1)->d2 == wts_mult[1] && (tmp-1)->d3 == inv_mult[0]), "thread-local memcheck failed!");
		// ASSERT(((tmp+0)->d0 == 0.50 && (tmp+0)->d1 == 0.50 && (tmp+0)->d2 == 0.50 && (tmp+0)->d3 == 0.50 && (tmp+0)->d4 == 0.50 && (tmp+0)->d5 == 0.50 && (tmp+0)->d6 == 0.50 && (tmp+0)->d7 == 0.50, "thread-local memcheck failed!");
		// ASSERT(((tmp+1)->d0 == 0.25 && (tmp+1)->d1 == 0.25 && (tmp+1)->d2 == 0.25 && (tmp+1)->d3 == 0.25 && (tmp+1)->d4 == 0.25 && (tmp+1)->d5 == 0.25 && (tmp+1)->d6 == 0.25 && (tmp+1)->d7 == 0.25, "thread-local memcheck failed!");
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
		#include "radix960_main_carry_loop.h"

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
		ASSERT(0x0 == cy960_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
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

void radix960_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-960 complex DIF FFT pass on the data in the length-N real vector A.
!
!   The data are stored in a 1-D zero-offset array, with 2^PAD_BITS 8-byte padding elements inserted
!   between every block of 2^DAT_BITS contiguous data. The array padding is to prevent data being accessed
!   in strides that are large powers of two and thus to minimize cache thrashing
!   (in cache-based microprocessor architectures) or bank conflicts (in supercomputers.)
!
!   See the documentation in radix64_dif_pass for further details on storage, indexing and DFT decomposition.
*/
	int k,l,mask3, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf, j,j1,j2,jt,jp;
	static int NDIVR,first_entry=TRUE,
			p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
		p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,
		p100,p110,p120,p130,p140,p150,p160,p170,p180,p190,p1a0,p1b0,p1c0,p1d0,p1e0,p1f0,
		p200,p210,p220,p230,p240,p250,p260,p270,p280,p290,p2a0,p2b0,p2c0,p2d0,p2e0,p2f0,
		p300,p310,p320,p330,p340,p350,p360,p370,p380,p390,p3a0,p3b0;
	static int plo[16], phi[RADIX>>4], p_out_hi[RADIX>>4];
	static int dif_o_offsets[64], t_offsets[64];
	uint64 i64;
	// Low parts [p0-f] of output-index perms - need 64x4-bits for each radix-64 DFT:
	const uint64 dif_perm16[60] = {
		0x0123456789abcdefull,0x1032547698badcfeull,0x23106754ab98efdcull,0x32017645ba89fecdull,
		0x54763201dcfeba89ull,0xcdefab9823106754ull,0x76450123fecd89abull,0xefdc98ba10325476ull,
		0xab98efdc67541032ull,0xba89fecd76450123ull,0x98badcfe54763201ull,0x45672310cdefab98ull,
		0xdcfeba8932017645ull,0x23106754ab98efdcull,0xfecd89ab01234567ull,0x1032547698badcfeull,
		0x67541032efdc98baull,0x76450123fecd89abull,0x54763201dcfeba89ull,0xcdefab9823106754ull,
		0x89abcdef45672310ull,0x98badcfe54763201ull,0xab98efdc67541032ull,0xba89fecd76450123ull,
		0xefdc98ba10325476ull,0xfecd89ab01234567ull,0xdcfeba8932017645ull,0x23106754ab98efdcull,
		0x45672310cdefab98ull,0x54763201dcfeba89ull,0x67541032efdc98baull,0x76450123fecd89abull,
		0x32017645ba89fecdull,0xab98efdc67541032ull,0x89abcdef45672310ull,0x98badcfe54763201ull,
		0xcdefab9823106754ull,0xdcfeba8932017645ull,0xefdc98ba10325476ull,0xfecd89ab01234567ull,
		0xba89fecd76450123ull,0x67541032efdc98baull,0x45672310cdefab98ull,0x54763201dcfeba89ull,
		0x1032547698badcfeull,0x89abcdef45672310ull,0x32017645ba89fecdull,0xab98efdc67541032ull,
		0x76450123fecd89abull,0xefdc98ba10325476ull,0xcdefab9823106754ull,0xdcfeba8932017645ull,
		0x98badcfe54763201ull,0x45672310cdefab98ull,0xba89fecd76450123ull,0x67541032efdc98baull,
		0x23106754ab98efdcull,0x32017645ba89fecdull,0x1032547698badcfeull,0x89abcdef45672310ull
	};
	static double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/
					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in-place DFT macros:
	struct complex t[RADIX], *tptr;

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		ASSERT((double *)t == &(t[0x00].re), "Unexpected value for Tmp-array-start pointer!");
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
	NDIVR <<= 4; p10 = NDIVR;	pf += ( (pf >> DAT_BITS) << PAD_BITS );
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
		pd0 = pc0 + NDIVR;		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pe0 = pd0 + NDIVR;		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		pf0 = pe0 + NDIVR;		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		p100= pf0 + NDIVR;		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
		p110= p100 + NDIVR;		p100+= ( (p100>> DAT_BITS) << PAD_BITS );
		p120= p110 + NDIVR;		p110+= ( (p110>> DAT_BITS) << PAD_BITS );
		p130= p120 + NDIVR;		p120+= ( (p120>> DAT_BITS) << PAD_BITS );
		p140= p130 + NDIVR;		p130+= ( (p130>> DAT_BITS) << PAD_BITS );
		p150= p140 + NDIVR;		p140+= ( (p140>> DAT_BITS) << PAD_BITS );
		p160= p150 + NDIVR;		p150+= ( (p150>> DAT_BITS) << PAD_BITS );
		p170= p160 + NDIVR;		p160+= ( (p160>> DAT_BITS) << PAD_BITS );
		p180= p170 + NDIVR;		p170+= ( (p170>> DAT_BITS) << PAD_BITS );
		p190= p180 + NDIVR;		p180+= ( (p180>> DAT_BITS) << PAD_BITS );
		p1a0= p190 + NDIVR;		p190+= ( (p190>> DAT_BITS) << PAD_BITS );
		p1b0= p1a0 + NDIVR;		p1a0+= ( (p1a0>> DAT_BITS) << PAD_BITS );
		p1c0= p1b0 + NDIVR;		p1b0+= ( (p1b0>> DAT_BITS) << PAD_BITS );
		p1d0= p1c0 + NDIVR;		p1c0+= ( (p1c0>> DAT_BITS) << PAD_BITS );
		p1e0= p1d0 + NDIVR;		p1d0+= ( (p1d0>> DAT_BITS) << PAD_BITS );
		p1f0= p1e0 + NDIVR;		p1e0+= ( (p1e0>> DAT_BITS) << PAD_BITS );
		p200= p1f0 + NDIVR;		p1f0+= ( (p1f0>> DAT_BITS) << PAD_BITS );
		p210= p200 + NDIVR;		p200+= ( (p200>> DAT_BITS) << PAD_BITS );
		p220= p210 + NDIVR;		p210+= ( (p210>> DAT_BITS) << PAD_BITS );
		p230= p220 + NDIVR;		p220+= ( (p220>> DAT_BITS) << PAD_BITS );
		p240= p230 + NDIVR;		p230+= ( (p230>> DAT_BITS) << PAD_BITS );
		p250= p240 + NDIVR;		p240+= ( (p240>> DAT_BITS) << PAD_BITS );
		p260= p250 + NDIVR;		p250+= ( (p250>> DAT_BITS) << PAD_BITS );
		p270= p260 + NDIVR;		p260+= ( (p260>> DAT_BITS) << PAD_BITS );
		p280= p270 + NDIVR;		p270+= ( (p270>> DAT_BITS) << PAD_BITS );
		p290= p280 + NDIVR;		p280+= ( (p280>> DAT_BITS) << PAD_BITS );
		p2a0= p290 + NDIVR;		p290+= ( (p290>> DAT_BITS) << PAD_BITS );
		p2b0= p2a0 + NDIVR;		p2a0+= ( (p2a0>> DAT_BITS) << PAD_BITS );
		p2c0= p2b0 + NDIVR;		p2b0+= ( (p2b0>> DAT_BITS) << PAD_BITS );
		p2d0= p2c0 + NDIVR;		p2c0+= ( (p2c0>> DAT_BITS) << PAD_BITS );
		p2e0= p2d0 + NDIVR;		p2d0+= ( (p2d0>> DAT_BITS) << PAD_BITS );
		p2f0= p2e0 + NDIVR;		p2e0+= ( (p2e0>> DAT_BITS) << PAD_BITS );
		p300= p2f0 + NDIVR;		p2f0+= ( (p2f0>> DAT_BITS) << PAD_BITS );
		p310= p300 + NDIVR;		p300+= ( (p300>> DAT_BITS) << PAD_BITS );
		p320= p310 + NDIVR;		p310+= ( (p310>> DAT_BITS) << PAD_BITS );
		p330= p320 + NDIVR;		p320+= ( (p320>> DAT_BITS) << PAD_BITS );
		p340= p330 + NDIVR;		p330+= ( (p330>> DAT_BITS) << PAD_BITS );
		p350= p340 + NDIVR;		p340+= ( (p340>> DAT_BITS) << PAD_BITS );
		p360= p350 + NDIVR;		p350+= ( (p350>> DAT_BITS) << PAD_BITS );
		p370= p360 + NDIVR;		p360+= ( (p360>> DAT_BITS) << PAD_BITS );
		p380= p370 + NDIVR;		p370+= ( (p370>> DAT_BITS) << PAD_BITS );
		p390= p380 + NDIVR;		p380+= ( (p380>> DAT_BITS) << PAD_BITS );
		p3a0= p390 + NDIVR;		p390+= ( (p390>> DAT_BITS) << PAD_BITS );
		p3b0= p3a0 + NDIVR;		p3a0+= ( (p3a0>> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			p3b0+= ( (p3b0>> DAT_BITS) << PAD_BITS );

		plo[ 0]= 0;plo[ 1]=p1;plo[ 2]=p2;plo[ 3]=p3;plo[ 4]=p4;plo[ 5]=p5;plo[ 6]=p6;plo[ 7]=p7;
		plo[ 8]=p8;plo[ 9]=p9;plo[10]=pa;plo[11]=pb;plo[12]=pc;plo[13]=pd;plo[14]=pe;plo[15]=pf;

		phi[0x00] =    0; phi[0x01] =  p10; phi[0x02] =  p20; phi[0x03] =  p30;
		phi[0x04] =  p40; phi[0x05] =  p50; phi[0x06] =  p60; phi[0x07] =  p70;
		phi[0x08] =  p80; phi[0x09] =  p90; phi[0x0a] =  pa0; phi[0x0b] =  pb0;
		phi[0x0c] =  pc0; phi[0x0d] =  pd0; phi[0x0e] =  pe0; phi[0x0f] =  pf0;
		phi[0x10] = p100; phi[0x11] = p110; phi[0x12] = p120; phi[0x13] = p130;
		phi[0x14] = p140; phi[0x15] = p150; phi[0x16] = p160; phi[0x17] = p170;
		phi[0x18] = p180; phi[0x19] = p190; phi[0x1a] = p1a0; phi[0x1b] = p1b0;
		phi[0x1c] = p1c0; phi[0x1d] = p1d0; phi[0x1e] = p1e0; phi[0x1f] = p1f0;
		phi[0x20] = p200; phi[0x21] = p210; phi[0x22] = p220; phi[0x23] = p230;
		phi[0x24] = p240; phi[0x25] = p250; phi[0x26] = p260; phi[0x27] = p270;
		phi[0x28] = p280; phi[0x29] = p290; phi[0x2a] = p2a0; phi[0x2b] = p2b0;
		phi[0x2c] = p2c0; phi[0x2d] = p2d0; phi[0x2e] = p2e0; phi[0x2f] = p2f0;
		phi[0x30] = p300; phi[0x31] = p310; phi[0x32] = p320; phi[0x33] = p330;
		phi[0x34] = p340; phi[0x35] = p350; phi[0x36] = p360; phi[0x37] = p370;
		phi[0x38] = p380; phi[0x39] = p390; phi[0x3a] = p3a0; phi[0x3b] = p3b0;

		// When first running unpermuted-output version of the DFT in order to extract the needed output-perm,
		// simply use the above phi data ... p_out_hi data are (mod p40), rcol comment lists the omitted p40-multiple:
		p_out_hi[0x00] =    0; p_out_hi[0x01] =  p10; p_out_hi[0x02] =  p20; p_out_hi[0x03] =  p30;	//  p00
		p_out_hi[0x04] =  p10; p_out_hi[0x05] =    0; p_out_hi[0x06] =  p30; p_out_hi[0x07] =  p20;	//  p80
		p_out_hi[0x08] =  p20; p_out_hi[0x09] =  p30; p_out_hi[0x0a] =  p10; p_out_hi[0x0b] =    0;	//  p40
		p_out_hi[0x0c] =  p10; p_out_hi[0x0d] =    0; p_out_hi[0x0e] =  p30; p_out_hi[0x0f] =  p20;	// p380
		p_out_hi[0x10] =  p20; p_out_hi[0x11] =  p30; p_out_hi[0x12] =  p10; p_out_hi[0x13] =    0;	// p340
		p_out_hi[0x14] =    0; p_out_hi[0x15] =  p10; p_out_hi[0x16] =  p20; p_out_hi[0x17] =  p30;	// p300
		p_out_hi[0x18] =  p20; p_out_hi[0x19] =  p30; p_out_hi[0x1a] =  p10; p_out_hi[0x1b] =    0;	// p2c0
		p_out_hi[0x1c] =    0; p_out_hi[0x1d] =  p10; p_out_hi[0x1e] =  p20; p_out_hi[0x1f] =  p30;	// p280
		p_out_hi[0x20] =  p30; p_out_hi[0x21] =  p20; p_out_hi[0x22] =    0; p_out_hi[0x23] =  p10;	// p240
		p_out_hi[0x24] =    0; p_out_hi[0x25] =  p10; p_out_hi[0x26] =  p20; p_out_hi[0x27] =  p30;	// p200
		p_out_hi[0x28] =  p30; p_out_hi[0x29] =  p20; p_out_hi[0x2a] =    0; p_out_hi[0x2b] =  p10;	// p1c0
		p_out_hi[0x2c] =  p10; p_out_hi[0x2d] =    0; p_out_hi[0x2e] =  p30; p_out_hi[0x2f] =  p20;	// p180
		p_out_hi[0x30] =  p30; p_out_hi[0x31] =  p20; p_out_hi[0x32] =    0; p_out_hi[0x33] =  p10;	// p140
		p_out_hi[0x34] =  p10; p_out_hi[0x35] =    0; p_out_hi[0x36] =  p30; p_out_hi[0x37] =  p20;	// p100
		p_out_hi[0x38] =  p20; p_out_hi[0x39] =  p30; p_out_hi[0x3a] =  p10; p_out_hi[0x3b] =    0;	//  pc0
		// We can on-the-fly compute the p40*j parts by parametrizing the sequence of j's:
		//	j = 0,2,1,14,13,12,11,10,9,8,7,6,5,4,3, in terms of an anscending loop index i = 0,1.,2,...,14 . 2 parts to thid:
		// [1] first get the leading 0,2,1 part via (3-i)&(-(i>0)), where the mask = -0 = 0 if i==0, -1 = 0xfff...f if i > 0 .
		// [2] Get the remaining descending-terms via (17-i), where now need to 0-mask these for i < 3,
		// as we do the terms of [1] for i > 3 (or i > 3, since (3-i) = 0 for i=3.)
		// Define mask3 = -(i<3), which = 0xfff... for i<3 and -0 = 0 for i >= 3. We then multiply [1] by mask3, [2] by ~mask3
		// and sum to obtain the desired indexing formula.
		//	mask3 = -(i<3);
		//	j = ( ( (3-i)&(-(i>0)) ) & mask3 ) + ( (17-i) & ~mask3 );

		// Set array offsets for radix-64 DFT in/outputs:
		// set 1 (t_offsets) is w.r.to: t-array:					// set 2 (dif_o_offsets) w.r.to a-array:
		for(l = 0; l < 64; l++) {									// [*** compute these on the fly below ***]
			t_offsets[l] = ((l<<4) - l)<<1;	// = l*30 (Need to double the complex 15-incr
							// due to the Radix-64 DFT casting t-array pointer to double
		}
	}

/*...The radix-960 pass is here.	*/

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

	/*...gather the needed data (960 64-bit complex) and do 64 radix-15 transforms...
	Twiddleless version arranges 64 sets of radix-15 DFT inputs as follows:
	0 in upper left corner, decrement 64 horizontally and 15 vertically, indexing modulo 960
	(we can auto-generate these by compiling test_fft_radix.c with -DTTYPE=0 -DRADIX=960, running
	the resulting executable and snarfing the first set of index-outputs, "DIF/DIT input-scramble array").

	If we subtract each row's common (mod 16) low p-part from all terms of the row
	as we do in the implementation to reduce the number of index offsets needing to be stored,
	we decrement 16 horizontally and 16 vertically, all (mod 960 = 0x3c0):

	00,380,340,300,2c0,280,240,200,1c0,180,140,100,c0,80,40		  0,380,340,300,2c0,280,240,200,1c0,180,140,100, c0, 80, 40 + 0
	3b1,371,331,2f1,2b1,271,231,1f1,1b1,171,131,f1,b1,71,31		3b0,370,330,2f0,2b0,270,230,1f0,1b0,170,130, f0, b0, 70, 30 + 1
	3a2,362,322,2e2,2a2,262,222,1e2,1a2,162,122,e2,a2,62,22		3a0,360,320,2e0,2a0,260,220,1e0,1a0,160,120, e0, a0, 60, 20 + 2
	393,353,313,2d3,293,253,213,1d3,193,153,113,d3,93,53,13		390,350,310,2d0,290,250,210,1d0,190,150,110, d0, 90, 50, 10 + 3
	384,344,304,2c4,284,244,204,1c4,184,144,104,c4,84,44,04		380,340,300,2c0,280,240,200,1c0,180,140,100, c0, 80, 40,  0 + 4
	375,335,2f5,2b5,275,235,1f5,1b5,175,135,f5,b5,75,35,3b5		370,330,2f0,2b0,270,230,1f0,1b0,170,130, f0, b0, 70, 30,3b0 + 5
	366,326,2e6,2a6,266,226,1e6,1a6,166,126,e6,a6,66,26,3a6		360,320,2e0,2a0,260,220,1e0,1a0,160,120, e0, a0, 60, 20,3a0 + 6
	357,317,2d7,297,257,217,1d7,197,157,117,d7,97,57,17,397		350,310,2d0,290,250,210,1d0,190,150,110, d0, 90, 50, 10,390 + 7
	348,308,2c8,288,248,208,1c8,188,148,108,c8,88,48,08,388		340,300,2c0,280,240,200,1c0,180,140,100, c0, 80, 40,  0,380 + 8
	339,2f9,2b9,279,239,1f9,1b9,179,139,f9,b9,79,39,3b9,379	=	330,2f0,2b0,270,230,1f0,1b0,170,130, f0, b0, 70, 30,3b0,370 + 9
	32a,2ea,2aa,26a,22a,1ea,1aa,16a,12a,ea,aa,6a,2a,3aa,36a		320,2e0,2a0,260,220,1e0,1a0,160,120, e0, a0, 60, 20,3a0,360 + a
	31b,2db,29b,25b,21b,1db,19b,15b,11b,db,9b,5b,1b,39b,35b		310,2d0,290,250,210,1d0,190,150,110, d0, 90, 50, 10,390,350 + b
	30c,2cc,28c,24c,20c,1cc,18c,14c,10c,cc,8c,4c,0c,38c,34c		300,2c0,280,240,200,1c0,180,140,100, c0, 80, 40,  0,380,340 + c
	2fd,2bd,27d,23d,1fd,1bd,17d,13d,fd,bd,7d,3d,3bd,37d,33d		2f0,2b0,270,230,1f0,1b0,170,130, f0, b0, 70, 30,3b0,370,330 + d
	2ee,2ae,26e,22e,1ee,1ae,16e,12e,ee,ae,6e,2e,3ae,36e,32e		2e0,2a0,260,220,1e0,1a0,160,120, e0, a0, 60, 20,3a0,360,320 + e
	2df,29f,25f,21f,1df,19f,15f,11f,df,9f,5f,1f,39f,35f,31f	   *2d0,290,250,210,1d0,190,150,110, d0, 90, 50, 10,390,350,310 + f	<*** Note leading-parts repeat every time
	2d0,290,250,210,1d0,190,150,110,d0,90,50,10,390,350,310	   *2d0,290,250,210,1d0,190,150,110, d0, 90, 50, 10,390,350,310 + 0	<*** low part rests to 0!
	2c1,281,241,201,1c1,181,141,101,c1,81,41,01,381,341,301		2c0,280,240,200,1c0,180,140,100, c0, 80, 40,  0,380,340,300 + 1
	2b2,272,232,1f2,1b2,172,132,f2,b2,72,32,3b2,372,332,2f2		2b0,270,230,1f0,1b0,170,130, f0, b0, 70, 30,3b0,370,330,2f0 + 2
	2a3,263,223,1e3,1a3,163,123,e3,a3,63,23,3a3,363,323,2e3		2a0,260,220,1e0,1a0,160,120, e0, a0, 60, 20,3a0,360,320,2e0 + 3
	294,254,214,1d4,194,154,114,d4,94,54,14,394,354,314,2d4		290,250,210,1d0,190,150,110, d0, 90, 50, 10,390,350,310,2d0 + 4
	285,245,205,1c5,185,145,105,c5,85,45,05,385,345,305,2c5		280,240,200,1c0,180,140,100, c0, 80, 40,  0,380,340,300,2c0 + 5
	276,236,1f6,1b6,176,136,f6,b6,76,36,3b6,376,336,2f6,2b6		270,230,1f0,1b0,170,130, f0, b0, 70, 30,3b0,370,330,2f0,2b0 + 6
	267,227,1e7,1a7,167,127,e7,a7,67,27,3a7,367,327,2e7,2a7		260,220,1e0,1a0,160,120, e0, a0, 60, 20,3a0,360,320,2e0,2a0 + 7
	258,218,1d8,198,158,118,d8,98,58,18,398,358,318,2d8,298		250,210,1d0,190,150,110, d0, 90, 50, 10,390,350,310,2d0,290 + 8
	249,209,1c9,189,149,109,c9,89,49,09,389,349,309,2c9,289	=	240,200,1c0,180,140,100, c0, 80, 40,  0,380,340,300,2c0,280 + 9
	23a,1fa,1ba,17a,13a,fa,ba,7a,3a,3ba,37a,33a,2fa,2ba,27a		230,1f0,1b0,170,130, f0, b0, 70, 30,3b0,370,330,2f0,2b0,270 + a
	22b,1eb,1ab,16b,12b,eb,ab,6b,2b,3ab,36b,32b,2eb,2ab,26b		220,1e0,1a0,160,120, e0, a0, 60, 20,3a0,360,320,2e0,2a0,260 + b
	21c,1dc,19c,15c,11c,dc,9c,5c,1c,39c,35c,31c,2dc,29c,25c		210,1d0,190,150,110, d0, 90, 50, 10,390,350,310,2d0,290,250 + c
	20d,1cd,18d,14d,10d,cd,8d,4d,0d,38d,34d,30d,2cd,28d,24d		200,1c0,180,140,100, c0, 80, 40,  0,380,340,300,2c0,280,240 + d
	1fe,1be,17e,13e,fe,be,7e,3e,3be,37e,33e,2fe,2be,27e,23e		1f0,1b0,170,130, f0, b0, 70, 30,3b0,370,330,2f0,2b0,270,230 + e
	1ef,1af,16f,12f,ef,af,6f,2f,3af,36f,32f,2ef,2af,26f,22f	   *1e0,1a0,160,120, e0, a0, 60, 20,3a0,360,320,2e0,2a0,260,220 + f
	1e0,1a0,160,120,e0,a0,60,20,3a0,360,320,2e0,2a0,260,220	   *1e0,1a0,160,120, e0, a0, 60, 20,3a0,360,320,2e0,2a0,260,220 + 0
	1d1,191,151,111,d1,91,51,11,391,351,311,2d1,291,251,211		1d0,190,150,110, d0, 90, 50, 10,390,350,310,2d0,290,250,210 + 1
	1c2,182,142,102,c2,82,42,02,382,342,302,2c2,282,242,202		1c0,180,140,100, c0, 80, 40,  0,380,340,300,2c0,280,240,200 + 2
	1b3,173,133,f3,b3,73,33,3b3,373,333,2f3,2b3,273,233,1f3		1b0,170,130, f0, b0, 70, 30,3b0,370,330,2f0,2b0,270,230,1f0 + 3
	1a4,164,124,e4,a4,64,24,3a4,364,324,2e4,2a4,264,224,1e4		1a0,160,120, e0, a0, 60, 20,3a0,360,320,2e0,2a0,260,220,1e0 + 4
	195,155,115,d5,95,55,15,395,355,315,2d5,295,255,215,1d5		190,150,110, d0, 90, 50, 10,390,350,310,2d0,290,250,210,1d0 + 5
	186,146,106,c6,86,46,06,386,346,306,2c6,286,246,206,1c6		180,140,100, c0, 80, 40,  0,380,340,300,2c0,280,240,200,1c0 + 6
	177,137,f7,b7,77,37,3b7,377,337,2f7,2b7,277,237,1f7,1b7		170,130, f0, b0, 70, 30,3b0,370,330,2f0,2b0,270,230,1f0,1b0 + 7
	168,128,e8,a8,68,28,3a8,368,328,2e8,2a8,268,228,1e8,1a8		160,120, e0, a0, 60, 20,3a0,360,320,2e0,2a0,260,220,1e0,1a0 + 8
	159,119,d9,99,59,19,399,359,319,2d9,299,259,219,1d9,199	=	150,110, d0, 90, 50, 10,390,350,310,2d0,290,250,210,1d0,190 + 9
	14a,10a,ca,8a,4a,0a,38a,34a,30a,2ca,28a,24a,20a,1ca,18a		140,100, c0, 80, 40,  0,380,340,300,2c0,280,240,200,1c0,180 + a
	13b,fb,bb,7b,3b,3bb,37b,33b,2fb,2bb,27b,23b,1fb,1bb,17b		130, f0, b0, 70, 30,3b0,370,330,2f0,2b0,270,230,1f0,1b0,170 + b
	12c,ec,ac,6c,2c,3ac,36c,32c,2ec,2ac,26c,22c,1ec,1ac,16c		120, e0, a0, 60, 20,3a0,360,320,2e0,2a0,260,220,1e0,1a0,160 + c
	11d,dd,9d,5d,1d,39d,35d,31d,2dd,29d,25d,21d,1dd,19d,15d		110, d0, 90, 50, 10,390,350,310,2d0,290,250,210,1d0,190,150 + d
	10e,ce,8e,4e,0e,38e,34e,30e,2ce,28e,24e,20e,1ce,18e,14e		100, c0, 80, 40,  0,380,340,300,2c0,280,240,200,1c0,180,140 + e
	ff,bf,7f,3f,3bf,37f,33f,2ff,2bf,27f,23f,1ff,1bf,17f,13f	   * f0, b0, 70, 30,3b0,370,330,2f0,2b0,270,230,1f0,1b0,170,130 + f
	f0,b0,70,30,3b0,370,330,2f0,2b0,270,230,1f0,1b0,170,130	   * f0, b0, 70, 30,3b0,370,330,2f0,2b0,270,230,1f0,1b0,170,130 + 0
	e1,a1,61,21,3a1,361,321,2e1,2a1,261,221,1e1,1a1,161,121		 e0, a0, 60, 20,3a0,360,320,2e0,2a0,260,220,1e0,1a0,160,120 + 1
	d2,92,52,12,392,352,312,2d2,292,252,212,1d2,192,152,112		 d0, 90, 50, 10,390,350,310,2d0,290,250,210,1d0,190,150,110 + 2
	c3,83,43,03,383,343,303,2c3,283,243,203,1c3,183,143,103		 c0, 80, 40,  0,380,340,300,2c0,280,240,200,1c0,180,140,100 + 3
	b4,74,34,3b4,374,334,2f4,2b4,274,234,1f4,1b4,174,134,f4		 b0, 70, 30,3b0,370,330,2f0,2b0,270,230,1f0,1b0,170,130, f0 + 4
	a5,65,25,3a5,365,325,2e5,2a5,265,225,1e5,1a5,165,125,e5		 a0, 60, 20,3a0,360,320,2e0,2a0,260,220,1e0,1a0,160,120, e0 + 5
	96,56,16,396,356,316,2d6,296,256,216,1d6,196,156,116,d6		 90, 50, 10,390,350,310,2d0,290,250,210,1d0,190,150,110, d0 + 6
	87,47,07,387,347,307,2c7,287,247,207,1c7,187,147,107,c7		 80, 40,  0,380,340,300,2c0,280,240,200,1c0,180,140,100, c0 + 7
	78,38,3b8,378,338,2f8,2b8,278,238,1f8,1b8,178,138,f8,b8		 70, 30,3b0,370,330,2f0,2b0,270,230,1f0,1b0,170,130, f0, b0 + 8
	69,29,3a9,369,329,2e9,2a9,269,229,1e9,1a9,169,129,e9,a9	=	 60, 20,3a0,360,320,2e0,2a0,260,220,1e0,1a0,160,120, e0, a0 + 9
	5a,1a,39a,35a,31a,2da,29a,25a,21a,1da,19a,15a,11a,da,9a		 50, 10,390,350,310,2d0,290,250,210,1d0,190,150,110, d0, 90 + a
	4b,0b,38b,34b,30b,2cb,28b,24b,20b,1cb,18b,14b,10b,cb,8b		 40,  0,380,340,300,2c0,280,240,200,1c0,180,140,100, c0, 80 + b
	3c,3bc,37c,33c,2fc,2bc,27c,23c,1fc,1bc,17c,13c,fc,bc,7c		 30,3b0,370,330,2f0,2b0,270,230,1f0,1b0,170,130, f0, b0, 70 + c
	2d,3ad,36d,32d,2ed,2ad,26d,22d,1ed,1ad,16d,12d,ed,ad,6d		 20,3a0,360,320,2e0,2a0,260,220,1e0,1a0,160,120, e0, a0, 60 + d
	1e,39e,35e,31e,2de,29e,25e,21e,1de,19e,15e,11e,de,9e,5e		 10,390,350,310,2d0,290,250,210,1d0,190,150,110, d0, 90, 50 + e
	0f,38f,34f,30f,2cf,28f,24f,20f,1cf,18f,14f,10f,cf,8f,4f		  0,380,340,300,2c0,280,240,200,1c0,180,140,100, c0, 80, 40 + f

	In order to make the result amenable to loop-based execution, we need to encode the indices
	to the left & right of the + in easily-computable-index fashion. This scheme has 3 ingredients:

	[0] RHS is simply i = 0,..,f, repeated 4x.

	For the LHS we need 2 ingredients:

	[1] A formula yielding the leading part of the leftmost term of each row: 0,3b,3a,...,1,0, INCLUDING REPEATED-PAIRS .
		The repeats mean we can't use scheme like we did for Radix-240, i.e. [in decimal] idx = (60-i) & (-(i > 0)) .
		Instead init idx = 60, zero the result iff (i == 0) prior to computing the remaining high-part terms of the row,
		and at end of each loop pass decrement if i != 15 (mod 16):

		idx = 60;
		for(i = 0; i < 64; i++) {
			idx &= (-(i>0));	// Mask = -0 = 0 if i==0, -1 = 0xfff...f if i > 0
			...
			idx -= ((i & 0xf) != 0xf);
		}

	[2] An efficient decrement-4 (mod 60) scheme to yield the remaining leading parts of each row's elements:
		idx -= 4; idx += (-(idx < 0))&60;
	*/
		tptr = t; k = 0;	// k takes place of 'idx' in above notes
		for(l = 0; l < 64; l++) {
		// [0] here:
			jp = plo[l&0xf];	// p0,..,f, repeated 4x
			jt = j1 + jp; jp += j2;
		// [1],[2] here:
			// Now get the remaining row terms and the resulting p* offsets:
			k0 = k;	// Need a separate master index k because k0 will get overwritten with phi[k0]
			k1 = k0-4; k1 += (-(k1 < 0))&60;		k0 = phi[k0];
			k2 = k1-4; k2 += (-(k2 < 0))&60;		k1 = phi[k1];
			k3 = k2-4; k3 += (-(k3 < 0))&60;		k2 = phi[k2];
			k4 = k3-4; k4 += (-(k4 < 0))&60;		k3 = phi[k3];
			k5 = k4-4; k5 += (-(k5 < 0))&60;		k4 = phi[k4];
			k6 = k5-4; k6 += (-(k6 < 0))&60;		k5 = phi[k5];
			k7 = k6-4; k7 += (-(k7 < 0))&60;		k6 = phi[k6];
			k8 = k7-4; k8 += (-(k8 < 0))&60;		k7 = phi[k7];
			k9 = k8-4; k9 += (-(k9 < 0))&60;		k8 = phi[k8];
			ka = k9-4; ka += (-(ka < 0))&60;		k9 = phi[k9];
			kb = ka-4; kb += (-(kb < 0))&60;		ka = phi[ka];
			kc = kb-4; kc += (-(kc < 0))&60;		kb = phi[kb];
			kd = kc-4; kd += (-(kd < 0))&60;		kc = phi[kc];
			ke = kd-4; ke += (-(ke < 0))&60;		kd = phi[kd];
													ke = phi[ke];
			RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],
				tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im
			);	tptr += 0xf;
			k -= ((l & 0xf) != 0xf);
			k += (-(k < 0))&60;
		}

		/*...and now do 15 radix-64 transforms.
		The required output permutation is, grouping things into blocks of 64, one for each DFT output set:

		000,001,002,003,004,005,006,007,008,009,00a,00b,00c,00d,00e,00f = 0,1,2,3,4,5,6,7,8,9,a,b,c,d,e,f + p00
		011,010,013,012,015,014,017,016,019,018,01b,01a,01d,01c,01f,01e = 1,0,3,2,5,4,7,6,9,8,b,a,d,c,f,e + p10
		022,023,021,020,026,027,025,024,02a,02b,029,028,02e,02f,02d,02c = 2,3,1,0,6,7,5,4,a,b,9,8,e,f,d,c + p20
		033,032,030,031,037,036,034,035,03b,03a,038,039,03f,03e,03c,03d = 3,2,0,1,7,6,4,5,b,a,8,9,f,e,c,d + p30

		095,094,097,096,093,092,090,091,09d,09c,09f,09e,09b,09a,098,099 = 5,4,7,6,3,2,0,1,d,c,f,e,b,a,8,9 + p90
		08c,08d,08e,08f,08a,08b,089,088,082,083,081,080,086,087,085,084 = c,d,e,f,a,b,9,8,2,3,1,0,6,7,5,4 + p80
		0b7,0b6,0b4,0b5,0b0,0b1,0b2,0b3,0bf,0be,0bc,0bd,0b8,0b9,0ba,0bb = 7,6,4,5,0,1,2,3,f,e,c,d,8,9,a,b + pb0
		0ae,0af,0ad,0ac,0a9,0a8,0ab,0aa,0a1,0a0,0a3,0a2,0a5,0a4,0a7,0a6 = e,f,d,c,9,8,b,a,1,0,3,2,5,4,7,6 + pa0

		06a,06b,069,068,06e,06f,06d,06c,066,067,065,064,061,060,063,062 = a,b,9,8,e,f,d,c,6,7,5,4,1,0,3,2 + p60
		07b,07a,078,079,07f,07e,07c,07d,077,076,074,075,070,071,072,073 = b,a,8,9,f,e,c,d,7,6,4,5,0,1,2,3 + p70
		059,058,05b,05a,05d,05c,05f,05e,055,054,057,056,053,052,050,051 = 9,8,b,a,d,c,f,e,5,4,7,6,3,2,0,1 + p50
		044,045,046,047,042,043,041,040,04c,04d,04e,04f,04a,04b,049,048 = 4,5,6,7,2,3,1,0,c,d,e,f,a,b,9,8 + p40

		39d,39c,39f,39e,39b,39a,398,399,393,392,390,391,397,396,394,395 = d,c,f,e,b,a,8,9,3,2,0,1,7,6,4,5 + p390
		382,383,381,380,386,387,385,384,38a,38b,389,388,38e,38f,38d,38c = 2,3,1,0,6,7,5,4,a,b,9,8,e,f,d,c + p380
		3bf,3be,3bc,3bd,3b8,3b9,3ba,3bb,3b0,3b1,3b2,3b3,3b4,3b5,3b6,3b7 = f,e,c,d,8,9,a,b,0,1,2,3,4,5,6,7 + p3b0
		3a1,3a0,3a3,3a2,3a5,3a4,3a7,3a6,3a9,3a8,3ab,3aa,3ad,3ac,3af,3ae = 1,0,3,2,5,4,7,6,9,8,b,a,d,c,f,e + p3a0

		366,367,365,364,361,360,363,362,36e,36f,36d,36c,369,368,36b,36a = 6,7,5,4,1,0,3,2,e,f,d,c,9,8,b,a + p360
		377,376,374,375,370,371,372,373,37f,37e,37c,37d,378,379,37a,37b = 7,6,4,5,0,1,2,3,f,e,c,d,8,9,a,b + p370
		355,354,357,356,353,352,350,351,35d,35c,35f,35e,35b,35a,358,359 = 5,4,7,6,3,2,0,1,d,c,f,e,b,a,8,9 + p350
		34c,34d,34e,34f,34a,34b,349,348,342,343,341,340,346,347,345,344 = c,d,e,f,a,b,9,8,2,3,1,0,6,7,5,4 + p340

		308,309,30a,30b,30c,30d,30e,30f,304,305,306,307,302,303,301,300 = 8,9,a,b,c,d,e,f,4,5,6,7,2,3,1,0 + p300
		319,318,31b,31a,31d,31c,31f,31e,315,314,317,316,313,312,310,311 = 9,8,b,a,d,c,f,e,5,4,7,6,3,2,0,1 + p310
		32a,32b,329,328,32e,32f,32d,32c,326,327,325,324,321,320,323,322 = a,b,9,8,e,f,d,c,6,7,5,4,1,0,3,2 + p320
		33b,33a,338,339,33f,33e,33c,33d,337,336,334,335,330,331,332,333 = b,a,8,9,f,e,c,d,7,6,4,5,0,1,2,3 + p330

		2ee,2ef,2ed,2ec,2e9,2e8,2eb,2ea,2e1,2e0,2e3,2e2,2e5,2e4,2e7,2e6 = e,f,d,c,9,8,b,a,1,0,3,2,5,4,7,6 + p2e0
		2ff,2fe,2fc,2fd,2f8,2f9,2fa,2fb,2f0,2f1,2f2,2f3,2f4,2f5,2f6,2f7 = f,e,c,d,8,9,a,b,0,1,2,3,4,5,6,7 + p2f0
		2dd,2dc,2df,2de,2db,2da,2d8,2d9,2d3,2d2,2d0,2d1,2d7,2d6,2d4,2d5 = d,c,f,e,b,a,8,9,3,2,0,1,7,6,4,5 + p2d0
		2c2,2c3,2c1,2c0,2c6,2c7,2c5,2c4,2ca,2cb,2c9,2c8,2ce,2cf,2cd,2cc = 2,3,1,0,6,7,5,4,a,b,9,8,e,f,d,c + p2c0

		284,285,286,287,282,283,281,280,28c,28d,28e,28f,28a,28b,289,288 = 4,5,6,7,2,3,1,0,c,d,e,f,a,b,9,8 + p280
		295,294,297,296,293,292,290,291,29d,29c,29f,29e,29b,29a,298,299 = 5,4,7,6,3,2,0,1,d,c,f,e,b,a,8,9 + p290
		2a6,2a7,2a5,2a4,2a1,2a0,2a3,2a2,2ae,2af,2ad,2ac,2a9,2a8,2ab,2aa = 6,7,5,4,1,0,3,2,e,f,d,c,9,8,b,a + p2a0
		2b7,2b6,2b4,2b5,2b0,2b1,2b2,2b3,2bf,2be,2bc,2bd,2b8,2b9,2ba,2bb = 7,6,4,5,0,1,2,3,f,e,c,d,8,9,a,b + p2b0

		273,272,270,271,277,276,274,275,27b,27a,278,279,27f,27e,27c,27d = 3,2,0,1,7,6,4,5,b,a,8,9,f,e,c,d + p270
		26a,26b,269,268,26e,26f,26d,26c,266,267,265,264,261,260,263,262 = a,b,9,8,e,f,d,c,6,7,5,4,1,0,3,2 + p260
		248,249,24a,24b,24c,24d,24e,24f,244,245,246,247,242,243,241,240 = 8,9,a,b,c,d,e,f,4,5,6,7,2,3,1,0 + p240
		259,258,25b,25a,25d,25c,25f,25e,255,254,257,256,253,252,250,251 = 9,8,b,a,d,c,f,e,5,4,7,6,3,2,0,1 + p250

		20c,20d,20e,20f,20a,20b,209,208,202,203,201,200,206,207,205,204 = c,d,e,f,a,b,9,8,2,3,1,0,6,7,5,4 + p200
		21d,21c,21f,21e,21b,21a,218,219,213,212,210,211,217,216,214,215 = d,c,f,e,b,a,8,9,3,2,0,1,7,6,4,5 + p210
		22e,22f,22d,22c,229,228,22b,22a,221,220,223,222,225,224,227,226 = e,f,d,c,9,8,b,a,1,0,3,2,5,4,7,6 + p220
		23f,23e,23c,23d,238,239,23a,23b,230,231,232,233,234,235,236,237 = f,e,c,d,8,9,a,b,0,1,2,3,4,5,6,7 + p230

		1fb,1fa,1f8,1f9,1ff,1fe,1fc,1fd,1f7,1f6,1f4,1f5,1f0,1f1,1f2,1f3 = b,a,8,9,f,e,c,d,7,6,4,5,0,1,2,3 + p1f0
		1e6,1e7,1e5,1e4,1e1,1e0,1e3,1e2,1ee,1ef,1ed,1ec,1e9,1e8,1eb,1ea = 6,7,5,4,1,0,3,2,e,f,d,c,9,8,b,a + p1e0
		1c4,1c5,1c6,1c7,1c2,1c3,1c1,1c0,1cc,1cd,1ce,1cf,1ca,1cb,1c9,1c8 = 4,5,6,7,2,3,1,0,c,d,e,f,a,b,9,8 + p1c0
		1d5,1d4,1d7,1d6,1d3,1d2,1d0,1d1,1dd,1dc,1df,1de,1db,1da,1d8,1d9 = 5,4,7,6,3,2,0,1,d,c,f,e,b,a,8,9 + p1d0

		191,190,193,192,195,194,197,196,199,198,19b,19a,19d,19c,19f,19e = 1,0,3,2,5,4,7,6,9,8,b,a,d,c,f,e + p190
		188,189,18a,18b,18c,18d,18e,18f,184,185,186,187,182,183,181,180 = 8,9,a,b,c,d,e,f,4,5,6,7,2,3,1,0 + p180
		1b3,1b2,1b0,1b1,1b7,1b6,1b4,1b5,1bb,1ba,1b8,1b9,1bf,1be,1bc,1bd = 3,2,0,1,7,6,4,5,b,a,8,9,f,e,c,d + p1b0
		1aa,1ab,1a9,1a8,1ae,1af,1ad,1ac,1a6,1a7,1a5,1a4,1a1,1a0,1a3,1a2 = a,b,9,8,e,f,d,c,6,7,5,4,1,0,3,2 + p1a0

		177,176,174,175,170,171,172,173,17f,17e,17c,17d,178,179,17a,17b = 7,6,4,5,0,1,2,3,f,e,c,d,8,9,a,b + p170
		16e,16f,16d,16c,169,168,16b,16a,161,160,163,162,165,164,167,166 = e,f,d,c,9,8,b,a,1,0,3,2,5,4,7,6 + p160
		14c,14d,14e,14f,14a,14b,149,148,142,143,141,140,146,147,145,144 = c,d,e,f,a,b,9,8,2,3,1,0,6,7,5,4 + p140
		15d,15c,15f,15e,15b,15a,158,159,153,152,150,151,157,156,154,155 = d,c,f,e,b,a,8,9,3,2,0,1,7,6,4,5 + p150

		119,118,11b,11a,11d,11c,11f,11e,115,114,117,116,113,112,110,111 = 9,8,b,a,d,c,f,e,5,4,7,6,3,2,0,1 + p110
		104,105,106,107,102,103,101,100,10c,10d,10e,10f,10a,10b,109,108 = 4,5,6,7,2,3,1,0,c,d,e,f,a,b,9,8 + p100
		13b,13a,138,139,13f,13e,13c,13d,137,136,134,135,130,131,132,133 = b,a,8,9,f,e,c,d,7,6,4,5,0,1,2,3 + p130
		126,127,125,124,121,120,123,122,12e,12f,12d,12c,129,128,12b,12a = 6,7,5,4,1,0,3,2,e,f,d,c,9,8,b,a + p120

		0e2,0e3,0e1,0e0,0e6,0e7,0e5,0e4,0ea,0eb,0e9,0e8,0ee,0ef,0ed,0ec = 2,3,1,0,6,7,5,4,a,b,9,8,e,f,d,c + pe0
		0f3,0f2,0f0,0f1,0f7,0f6,0f4,0f5,0fb,0fa,0f8,0f9,0ff,0fe,0fc,0fd = 3,2,0,1,7,6,4,5,b,a,8,9,f,e,c,d + pf0
		0d1,0d0,0d3,0d2,0d5,0d4,0d7,0d6,0d9,0d8,0db,0da,0dd,0dc,0df,0de = 1,0,3,2,5,4,7,6,9,8,b,a,d,c,f,e + pd0
		0c8,0c9,0ca,0cb,0cc,0cd,0ce,0cf,0c4,0c5,0c6,0c7,0c2,0c3,0c1,0c0 = 8,9,a,b,c,d,e,f,4,5,6,7,2,3,1,0 + pc0

		As documented in e.g. radix240_dit_pass1, simply encode each length-16 subperm in the rcol as a hex-char string.
		This means we must extract each p-offset in little-endian fashion, e.g. low 4 bits have rightmost p-offset above.
		*/
		for(l = 0; l < 15; l++) {
			// Compute dif_o_offsets for current 64-block:
			for(k = 0; k < 4; k++) {
				jp = (l<<2) + k;
				i64 = dif_perm16[jp];
 				jp = p_out_hi[jp];	// Add the p[0123]0 term stored in the current p_out_hi quartet
				// And lastly add in the low parts - p-offset indices encoded in little-endian hex-char fashion:
				jt = (k<<4);
				k0 = plo[(i64 >> 60)&0xf];		dif_o_offsets[jt+0x0] = jp + k0;
				k1 = plo[(i64 >> 56)&0xf];		dif_o_offsets[jt+0x1] = jp + k1;
				k2 = plo[(i64 >> 52)&0xf];		dif_o_offsets[jt+0x2] = jp + k2;
				k3 = plo[(i64 >> 48)&0xf];		dif_o_offsets[jt+0x3] = jp + k3;
				k4 = plo[(i64 >> 44)&0xf];		dif_o_offsets[jt+0x4] = jp + k4;
				k5 = plo[(i64 >> 40)&0xf];		dif_o_offsets[jt+0x5] = jp + k5;
				k6 = plo[(i64 >> 36)&0xf];		dif_o_offsets[jt+0x6] = jp + k6;
				k7 = plo[(i64 >> 32)&0xf];		dif_o_offsets[jt+0x7] = jp + k7;
				k8 = plo[(i64 >> 28)&0xf];		dif_o_offsets[jt+0x8] = jp + k8;
				k9 = plo[(i64 >> 24)&0xf];		dif_o_offsets[jt+0x9] = jp + k9;
				ka = plo[(i64 >> 20)&0xf];		dif_o_offsets[jt+0xa] = jp + ka;
				kb = plo[(i64 >> 16)&0xf];		dif_o_offsets[jt+0xb] = jp + kb;
				kc = plo[(i64 >> 12)&0xf];		dif_o_offsets[jt+0xc] = jp + kc;
				kd = plo[(i64 >>  8)&0xf];		dif_o_offsets[jt+0xd] = jp + kd;
				ke = plo[(i64 >>  4)&0xf];		dif_o_offsets[jt+0xe] = jp + ke;
				kf = plo[(i64      )&0xf];		dif_o_offsets[jt+0xf] = jp + kf;
			}
			// Compute the p40 multiplier which acts as the base offset:
			mask3 = -(l<3);
			jp = ( ( (3-l)&(-(l>0)) ) & mask3 ) + ( (17-l) & ~mask3 );	// j = 0,2,1,14,13,12,11,10,9,8,7,6,5,4,3
			jt = j1 + phi[jp<<2];		//  Add main-array index  to j*p40 - dif_o_offsets are all relative to this
		/*
		In-order indexing we used initially, until extraction of operm:
			// The (l<<2) term effects a +p40 increment between radix-64 DFTs
			// (we cannot simply do += p40, since the array padding scheme means that e.g. p40*2 != p80):
			jt = j1 + phi[l<<2];	// = p0,p40,...,p380
			for(k = 0; k < 64; k++) { dif_o_offsets[k] = phi[(k>>4)] + plo[k&0xf]; }
		*/
			//	NOTE that RADIX_64_DIF outputs are IN-ORDER rather than BR:
			RADIX_64_DIF((double *)(t+l),t_offsets,1, (a+jt),dif_o_offsets,RE_IM_STRIDE);
		}
	}
}

/**************/

void radix960_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-960 complex inverse DIT FFT pass on the data in the length-N real vector A.
*/
	int k,l,mask3, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf, j,j1,j2,jt,jp;
	static int NDIVR,first_entry=TRUE,
			p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
		p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,
		p100,p110,p120,p130,p140,p150,p160,p170,p180,p190,p1a0,p1b0,p1c0,p1d0,p1e0,p1f0,
		p200,p210,p220,p230,p240,p250,p260,p270,p280,p290,p2a0,p2b0,p2c0,p2d0,p2e0,p2f0,
		p300,p310,p320,p330,p340,p350,p360,p370,p380,p390,p3a0,p3b0;
	static int plo[16], phi[RADIX>>4];
	static int t_offsets[64], dit_i_offsets[64];
	uint64 i64;
	// Low parts [p0-f] of input-index perms - 16 distinct ones, 4 lookups for each radix-64 DFT:
	const uint64 dit_perm16[16] = {
		0x01327654fedcba98ull,0xfedcba9876543210ull,0x5467102398abefcdull,0x98abefcd10236745ull,
		0xab89cdfe23014576ull,0x23014576cdfe89baull,0x76543210ba98dcefull,0x10236745efcdab89ull,
		0xcdfe89ba45760132ull,0xba98dcef32105467ull,0xefcdab8967452301ull,0x4576013289bafedcull,
		0x32105467dcef98abull,0x67452301ab89cdfeull,0x89bafedc01327654ull,0xdcef98ab54671023ull
	};
	// Pattern of access of above-encoded 16-perms, 4 for each DFT-64: Tempting to encoded via hex-hars,
	// but not worth it, just use byte array and get half the maximum possible compression ratio but simpler code:
	const uint8 dit_perm16_idx[60] = {
		0x0,0x1,0x1,0x1,0x2,0x2,0x2,0x3,0x4,0x5,0x4,0x5,0x6,0x6,0x6,
		0x6,0x3,0x3,0x3,0x7,0x5,0x8,0x5,0x8,0x9,0x9,0x9,0x9,0x7,0x7,
		0x7,0xa,0x8,0xb,0xb,0xb,0xc,0xc,0xc,0xc,0xa,0xd,0xa,0xd,0xb,
		0xe,0xe,0xe,0xf,0xf,0xf,0x2,0xd,0x4,0xd,0x4,0xe,0x0,0x0,0x0
	};
	// Pattern of access of phi[] elements for the above-encoded 16-perms, 4 for each DFT-64:
	const uint8 dit_perm16_phidx[60] = {
		0x00,0x01,0x03,0x02,0x09,0x08,0x0a,0x0b,0x06,0x07,0x04,0x05,0x17,0x16,0x15,
		0x14,0x11,0x10,0x12,0x13,0x0e,0x0f,0x0c,0x0d,0x1f,0x1e,0x1d,0x1c,0x19,0x18,
		0x1a,0x1b,0x20,0x21,0x23,0x22,0x27,0x26,0x25,0x24,0x2e,0x2f,0x2c,0x2d,0x28,
		0x29,0x2b,0x2a,0x39,0x38,0x3a,0x3b,0x36,0x37,0x34,0x35,0x30,0x31,0x33,0x32
	};
	static double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/
					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in-place DFT macros:
	struct complex t[RADIX], *tptr;

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		ASSERT((double *)t == &(t[0x00].re), "Unexpected value for Tmp-array-start pointer!");
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
	NDIVR <<= 4; p10 = NDIVR;	pf += ( (pf >> DAT_BITS) << PAD_BITS );
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
		pd0 = pc0 + NDIVR;		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pe0 = pd0 + NDIVR;		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		pf0 = pe0 + NDIVR;		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		p100= pf0 + NDIVR;		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
		p110= p100 + NDIVR;		p100+= ( (p100>> DAT_BITS) << PAD_BITS );
		p120= p110 + NDIVR;		p110+= ( (p110>> DAT_BITS) << PAD_BITS );
		p130= p120 + NDIVR;		p120+= ( (p120>> DAT_BITS) << PAD_BITS );
		p140= p130 + NDIVR;		p130+= ( (p130>> DAT_BITS) << PAD_BITS );
		p150= p140 + NDIVR;		p140+= ( (p140>> DAT_BITS) << PAD_BITS );
		p160= p150 + NDIVR;		p150+= ( (p150>> DAT_BITS) << PAD_BITS );
		p170= p160 + NDIVR;		p160+= ( (p160>> DAT_BITS) << PAD_BITS );
		p180= p170 + NDIVR;		p170+= ( (p170>> DAT_BITS) << PAD_BITS );
		p190= p180 + NDIVR;		p180+= ( (p180>> DAT_BITS) << PAD_BITS );
		p1a0= p190 + NDIVR;		p190+= ( (p190>> DAT_BITS) << PAD_BITS );
		p1b0= p1a0 + NDIVR;		p1a0+= ( (p1a0>> DAT_BITS) << PAD_BITS );
		p1c0= p1b0 + NDIVR;		p1b0+= ( (p1b0>> DAT_BITS) << PAD_BITS );
		p1d0= p1c0 + NDIVR;		p1c0+= ( (p1c0>> DAT_BITS) << PAD_BITS );
		p1e0= p1d0 + NDIVR;		p1d0+= ( (p1d0>> DAT_BITS) << PAD_BITS );
		p1f0= p1e0 + NDIVR;		p1e0+= ( (p1e0>> DAT_BITS) << PAD_BITS );
		p200= p1f0 + NDIVR;		p1f0+= ( (p1f0>> DAT_BITS) << PAD_BITS );
		p210= p200 + NDIVR;		p200+= ( (p200>> DAT_BITS) << PAD_BITS );
		p220= p210 + NDIVR;		p210+= ( (p210>> DAT_BITS) << PAD_BITS );
		p230= p220 + NDIVR;		p220+= ( (p220>> DAT_BITS) << PAD_BITS );
		p240= p230 + NDIVR;		p230+= ( (p230>> DAT_BITS) << PAD_BITS );
		p250= p240 + NDIVR;		p240+= ( (p240>> DAT_BITS) << PAD_BITS );
		p260= p250 + NDIVR;		p250+= ( (p250>> DAT_BITS) << PAD_BITS );
		p270= p260 + NDIVR;		p260+= ( (p260>> DAT_BITS) << PAD_BITS );
		p280= p270 + NDIVR;		p270+= ( (p270>> DAT_BITS) << PAD_BITS );
		p290= p280 + NDIVR;		p280+= ( (p280>> DAT_BITS) << PAD_BITS );
		p2a0= p290 + NDIVR;		p290+= ( (p290>> DAT_BITS) << PAD_BITS );
		p2b0= p2a0 + NDIVR;		p2a0+= ( (p2a0>> DAT_BITS) << PAD_BITS );
		p2c0= p2b0 + NDIVR;		p2b0+= ( (p2b0>> DAT_BITS) << PAD_BITS );
		p2d0= p2c0 + NDIVR;		p2c0+= ( (p2c0>> DAT_BITS) << PAD_BITS );
		p2e0= p2d0 + NDIVR;		p2d0+= ( (p2d0>> DAT_BITS) << PAD_BITS );
		p2f0= p2e0 + NDIVR;		p2e0+= ( (p2e0>> DAT_BITS) << PAD_BITS );
		p300= p2f0 + NDIVR;		p2f0+= ( (p2f0>> DAT_BITS) << PAD_BITS );
		p310= p300 + NDIVR;		p300+= ( (p300>> DAT_BITS) << PAD_BITS );
		p320= p310 + NDIVR;		p310+= ( (p310>> DAT_BITS) << PAD_BITS );
		p330= p320 + NDIVR;		p320+= ( (p320>> DAT_BITS) << PAD_BITS );
		p340= p330 + NDIVR;		p330+= ( (p330>> DAT_BITS) << PAD_BITS );
		p350= p340 + NDIVR;		p340+= ( (p340>> DAT_BITS) << PAD_BITS );
		p360= p350 + NDIVR;		p350+= ( (p350>> DAT_BITS) << PAD_BITS );
		p370= p360 + NDIVR;		p360+= ( (p360>> DAT_BITS) << PAD_BITS );
		p380= p370 + NDIVR;		p370+= ( (p370>> DAT_BITS) << PAD_BITS );
		p390= p380 + NDIVR;		p380+= ( (p380>> DAT_BITS) << PAD_BITS );
		p3a0= p390 + NDIVR;		p390+= ( (p390>> DAT_BITS) << PAD_BITS );
		p3b0= p3a0 + NDIVR;		p3a0+= ( (p3a0>> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			p3b0+= ( (p3b0>> DAT_BITS) << PAD_BITS );

		plo[ 0]= 0;plo[ 1]=p1;plo[ 2]=p2;plo[ 3]=p3;plo[ 4]=p4;plo[ 5]=p5;plo[ 6]=p6;plo[ 7]=p7;
		plo[ 8]=p8;plo[ 9]=p9;plo[10]=pa;plo[11]=pb;plo[12]=pc;plo[13]=pd;plo[14]=pe;plo[15]=pf;

		phi[0x00] =    0; phi[0x01] =  p10; phi[0x02] =  p20; phi[0x03] =  p30;
		phi[0x04] =  p40; phi[0x05] =  p50; phi[0x06] =  p60; phi[0x07] =  p70;
		phi[0x08] =  p80; phi[0x09] =  p90; phi[0x0a] =  pa0; phi[0x0b] =  pb0;
		phi[0x0c] =  pc0; phi[0x0d] =  pd0; phi[0x0e] =  pe0; phi[0x0f] =  pf0;
		phi[0x10] = p100; phi[0x11] = p110; phi[0x12] = p120; phi[0x13] = p130;
		phi[0x14] = p140; phi[0x15] = p150; phi[0x16] = p160; phi[0x17] = p170;
		phi[0x18] = p180; phi[0x19] = p190; phi[0x1a] = p1a0; phi[0x1b] = p1b0;
		phi[0x1c] = p1c0; phi[0x1d] = p1d0; phi[0x1e] = p1e0; phi[0x1f] = p1f0;
		phi[0x20] = p200; phi[0x21] = p210; phi[0x22] = p220; phi[0x23] = p230;
		phi[0x24] = p240; phi[0x25] = p250; phi[0x26] = p260; phi[0x27] = p270;
		phi[0x28] = p280; phi[0x29] = p290; phi[0x2a] = p2a0; phi[0x2b] = p2b0;
		phi[0x2c] = p2c0; phi[0x2d] = p2d0; phi[0x2e] = p2e0; phi[0x2f] = p2f0;
		phi[0x30] = p300; phi[0x31] = p310; phi[0x32] = p320; phi[0x33] = p330;
		phi[0x34] = p340; phi[0x35] = p350; phi[0x36] = p360; phi[0x37] = p370;
		phi[0x38] = p380; phi[0x39] = p390; phi[0x3a] = p3a0; phi[0x3b] = p3b0;

		for(l = 0; l < 64; l++) {
			t_offsets[l] = ((l<<4) - l)<<1;	// = l*30 (Need to double the complex 15-incr due to the Radix-64 DFT casting t-array pointer to double
		}
	}

/*...The radix-960 pass is here.	*/

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
	Gather the needed data (960 64-bit complex) and do 15 radix-64 transforms:

	Twiddleless version uses same linear-index-vector-form permutation as in DIF -
	Remember, inputs to DIT are bit-reversed, so use output of test_fft_radix() with
	TTYPE=0 to auto-generate needed input-index permutation:

	Combined DIT input-scramble array =
		 00, 01, 03, 02, 07, 06, 05, 04, 0f, 0e, 0d, 0c, 0b, 0a, 09, 08   0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p000   [a] + p000
		 1f, 1e, 1d, 1c, 1b, 1a, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10 = f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p010 = [b] + p010
		 3f, 3e, 3d, 3c, 3b, 3a, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30   f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p030   [b] + p030
		 2f, 2e, 2d, 2c, 2b, 2a, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20   f,e,d,c,b,a,9,8,7,6,5,4,3,2,1,0 + p020   [b] + p020
		 95, 94, 96, 97, 91, 90, 92, 93, 99, 98, 9a, 9b, 9e, 9f, 9c, 9d   5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p090   [c] + p090
		 85, 84, 86, 87, 81, 80, 82, 83, 89, 88, 8a, 8b, 8e, 8f, 8c, 8d = 5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p080 = [c] + p080
		 a5, a4, a6, a7, a1, a0, a2, a3, a9, a8, aa, ab, ae, af, ac, ad   5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p0a0   [c] + p0a0
		 b9, b8, ba, bb, be, bf, bc, bd, b1, b0, b2, b3, b6, b7, b4, b5   9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p0b0   [d] + p0b0
		 6a, 6b, 68, 69, 6c, 6d, 6f, 6e, 62, 63, 60, 61, 64, 65, 67, 66   a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p060   [e] + p060
		 72, 73, 70, 71, 74, 75, 77, 76, 7c, 7d, 7f, 7e, 78, 79, 7b, 7a = 2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p070 = [f] + p070
		 4a, 4b, 48, 49, 4c, 4d, 4f, 4e, 42, 43, 40, 41, 44, 45, 47, 46   a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p040   [e] + p040
		 52, 53, 50, 51, 54, 55, 57, 56, 5c, 5d, 5f, 5e, 58, 59, 5b, 5a   2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p050   [f] + p050
		177,176,175,174,173,172,171,170,17b,17a,179,178,17d,17c,17e,17f   7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + p170   [g] + p170
		167,166,165,164,163,162,161,160,16b,16a,169,168,16d,16c,16e,16f = 7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + p160 = [g] + p160
		157,156,155,154,153,152,151,150,15b,15a,159,158,15d,15c,15e,15f   7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + p150   [g] + p150
		147,146,145,144,143,142,141,140,14b,14a,149,148,14d,14c,14e,14f   7,6,5,4,3,2,1,0,b,a,9,8,d,c,e,f + p140   [g] + p140
		119,118,11a,11b,11e,11f,11c,11d,111,110,112,113,116,117,114,115   9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p110   [d] + p110
		109,108,10a,10b,10e,10f,10c,10d,101,100,102,103,106,107,104,105 = 9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p100 = [d] + p100
		129,128,12a,12b,12e,12f,12c,12d,121,120,122,123,126,127,124,125   9,8,a,b,e,f,c,d,1,0,2,3,6,7,4,5 + p120   [d] + p120
		131,130,132,133,136,137,134,135,13e,13f,13c,13d,13a,13b,138,139   1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + p130   [h] + p130
		 e2, e3, e0, e1, e4, e5, e7, e6, ec, ed, ef, ee, e8, e9, eb, ea   2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p0e0   [f] + p0e0
		 fc, fd, ff, fe, f8, f9, fb, fa, f4, f5, f7, f6, f0, f1, f3, f2 = c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p0f0 = [i] + p0f0
		 c2, c3, c0, c1, c4, c5, c7, c6, cc, cd, cf, ce, c8, c9, cb, ca   2,3,0,1,4,5,7,6,c,d,f,e,8,9,b,a + p0c0   [f] + p0c0
		 dc, dd, df, de, d8, d9, db, da, d4, d5, d7, d6, d0, d1, d3, d2   c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p0d0   [i] + p0d0
		1fb,1fa,1f9,1f8,1fd,1fc,1fe,1ff,1f3,1f2,1f1,1f0,1f5,1f4,1f6,1f7   b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p1f0   [j] + p1f0
		1eb,1ea,1e9,1e8,1ed,1ec,1ee,1ef,1e3,1e2,1e1,1e0,1e5,1e4,1e6,1e7 = b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p1e0 = [j] + p1e0
		1db,1da,1d9,1d8,1dd,1dc,1de,1df,1d3,1d2,1d1,1d0,1d5,1d4,1d6,1d7   b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p1d0   [j] + p1d0
		1cb,1ca,1c9,1c8,1cd,1cc,1ce,1cf,1c3,1c2,1c1,1c0,1c5,1c4,1c6,1c7   b,a,9,8,d,c,e,f,3,2,1,0,5,4,6,7 + p1c0   [j] + p1c0
		191,190,192,193,196,197,194,195,19e,19f,19c,19d,19a,19b,198,199   1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + p190   [h] + p190
		181,180,182,183,186,187,184,185,18e,18f,18c,18d,18a,18b,188,189 = 1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + p180 = [h] + p180
		1a1,1a0,1a2,1a3,1a6,1a7,1a4,1a5,1ae,1af,1ac,1ad,1aa,1ab,1a8,1a9   1,0,2,3,6,7,4,5,e,f,c,d,a,b,8,9 + p1a0   [h] + p1a0
		1be,1bf,1bc,1bd,1ba,1bb,1b8,1b9,1b6,1b7,1b4,1b5,1b2,1b3,1b0,1b1   e,f,c,d,a,b,8,9,6,7,4,5,2,3,0,1 + p1b0   [k] + p1b0
		20c,20d,20f,20e,208,209,20b,20a,204,205,207,206,200,201,203,202   c,d,f,e,8,9,b,a,4,5,7,6,0,1,3,2 + p200   [i] + p200
		214,215,217,216,210,211,213,212,218,219,21b,21a,21f,21e,21d,21c = 4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + p210 = [l] + p210
		234,235,237,236,230,231,233,232,238,239,23b,23a,23f,23e,23d,23c   4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + p230   [l] + p230
		224,225,227,226,220,221,223,222,228,229,22b,22a,22f,22e,22d,22c   4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + p220   [l] + p220
		273,272,271,270,275,274,276,277,27d,27c,27e,27f,279,278,27a,27b   3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p270   [m] + p270
		263,262,261,260,265,264,266,267,26d,26c,26e,26f,269,268,26a,26b = 3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p260 = [m] + p260
		253,252,251,250,255,254,256,257,25d,25c,25e,25f,259,258,25a,25b   3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p250   [m] + p250
		243,242,241,240,245,244,246,247,24d,24c,24e,24f,249,248,24a,24b   3,2,1,0,5,4,6,7,d,c,e,f,9,8,a,b + p240   [m] + p240
		2ee,2ef,2ec,2ed,2ea,2eb,2e8,2e9,2e6,2e7,2e4,2e5,2e2,2e3,2e0,2e1   e,f,c,d,a,b,8,9,6,7,4,5,2,3,0,1 + p2e0   [k] + p2e0
		2f6,2f7,2f4,2f5,2f2,2f3,2f0,2f1,2fa,2fb,2f8,2f9,2fc,2fd,2ff,2fe = 6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + p2f0 = [n] + p2f0
		2ce,2cf,2cc,2cd,2ca,2cb,2c8,2c9,2c6,2c7,2c4,2c5,2c2,2c3,2c0,2c1   e,f,c,d,a,b,8,9,6,7,4,5,2,3,0,1 + p2c0   [k] + p2c0
		2d6,2d7,2d4,2d5,2d2,2d3,2d0,2d1,2da,2db,2d8,2d9,2dc,2dd,2df,2de   6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + p2d0   [n] + p2d0
		284,285,287,286,280,281,283,282,288,289,28b,28a,28f,28e,28d,28c   4,5,7,6,0,1,3,2,8,9,b,a,f,e,d,c + p280   [l] + p280
		298,299,29b,29a,29f,29e,29d,29c,290,291,293,292,297,296,295,294 = 8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + p290 = [o] + p290
		2b8,2b9,2bb,2ba,2bf,2be,2bd,2bc,2b0,2b1,2b3,2b2,2b7,2b6,2b5,2b4   8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + p2b0   [o] + p2b0
		2a8,2a9,2ab,2aa,2af,2ae,2ad,2ac,2a0,2a1,2a3,2a2,2a7,2a6,2a5,2a4   8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + p2a0   [o] + p2a0
		39d,39c,39e,39f,399,398,39a,39b,395,394,396,397,391,390,392,393   d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + p390   [p] + p390
		38d,38c,38e,38f,389,388,38a,38b,385,384,386,387,381,380,382,383 = d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + p380 = [p] + p380
		3ad,3ac,3ae,3af,3a9,3a8,3aa,3ab,3a5,3a4,3a6,3a7,3a1,3a0,3a2,3a3   d,c,e,f,9,8,a,b,5,4,6,7,1,0,2,3 + p3a0   [p] + p3a0
		3b5,3b4,3b6,3b7,3b1,3b0,3b2,3b3,3b9,3b8,3ba,3bb,3be,3bf,3bc,3bd   5,4,6,7,1,0,2,3,9,8,a,b,e,f,c,d + p3b0   [c] + p3b0
		366,367,364,365,362,363,360,361,36a,36b,368,369,36c,36d,36f,36e   6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + p360   [n] + p360
		37a,37b,378,379,37c,37d,37f,37e,372,373,370,371,374,375,377,376 = a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p370 = [e] + p370
		346,347,344,345,342,343,340,341,34a,34b,348,349,34c,34d,34f,34e   6,7,4,5,2,3,0,1,a,b,8,9,c,d,f,e + p340   [n] + p340
		35a,35b,358,359,35c,35d,35f,35e,352,353,350,351,354,355,357,356   a,b,8,9,c,d,f,e,2,3,0,1,4,5,7,6 + p350   [e] + p350
		308,309,30b,30a,30f,30e,30d,30c,300,301,303,302,307,306,305,304   8,9,b,a,f,e,d,c,0,1,3,2,7,6,5,4 + p300   [o] + p300
		310,311,313,312,317,316,315,314,31f,31e,31d,31c,31b,31a,319,318 = 0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p310 = [a] + p310
		330,331,333,332,337,336,335,334,33f,33e,33d,33c,33b,33a,339,338   0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p330   [a] + p330
		320,321,323,322,327,326,325,324,32f,32e,32d,32c,32b,32a,329,328   0,1,3,2,7,6,5,4,f,e,d,c,b,a,9,8 + p320   [a] + p320

	As summarized in the rcol, there are 16 distinct low-part patterns a-p, which appear 4344444434344443 times, resp.
	We simply encode each distinct 16-perm as a hex-char string, with a length-60 aux array storing the index of the
	thus-encoded 16-perm needed for the current 16 data [Require 4 such lookups for each DFT-64].
	*/
		for(l = 0; l < 15; l++) {
			// Compute dit_i_offsets for current 64-block:
			for(k = 0; k < 4; k++) {
				jp = (l<<2) + k;
				i64 = dit_perm16[dit_perm16_idx[jp]];
 				jp = j1 + phi[dit_perm16_phidx[jp]];	// Main-array base index plus hi-part p-offset
				// And lastly add in the low parts - p-offset indices encoded in little-endian hex-char fashion:
				jt = (k<<4);
				k0 = plo[(i64 >> 60)&0xf];		dit_i_offsets[jt+0x0] = jp + k0;
				k1 = plo[(i64 >> 56)&0xf];		dit_i_offsets[jt+0x1] = jp + k1;
				k2 = plo[(i64 >> 52)&0xf];		dit_i_offsets[jt+0x2] = jp + k2;
				k3 = plo[(i64 >> 48)&0xf];		dit_i_offsets[jt+0x3] = jp + k3;
				k4 = plo[(i64 >> 44)&0xf];		dit_i_offsets[jt+0x4] = jp + k4;
				k5 = plo[(i64 >> 40)&0xf];		dit_i_offsets[jt+0x5] = jp + k5;
				k6 = plo[(i64 >> 36)&0xf];		dit_i_offsets[jt+0x6] = jp + k6;
				k7 = plo[(i64 >> 32)&0xf];		dit_i_offsets[jt+0x7] = jp + k7;
				k8 = plo[(i64 >> 28)&0xf];		dit_i_offsets[jt+0x8] = jp + k8;
				k9 = plo[(i64 >> 24)&0xf];		dit_i_offsets[jt+0x9] = jp + k9;
				ka = plo[(i64 >> 20)&0xf];		dit_i_offsets[jt+0xa] = jp + ka;
				kb = plo[(i64 >> 16)&0xf];		dit_i_offsets[jt+0xb] = jp + kb;
				kc = plo[(i64 >> 12)&0xf];		dit_i_offsets[jt+0xc] = jp + kc;
				kd = plo[(i64 >>  8)&0xf];		dit_i_offsets[jt+0xd] = jp + kd;
				ke = plo[(i64 >>  4)&0xf];		dit_i_offsets[jt+0xe] = jp + ke;
				kf = plo[(i64      )&0xf];		dit_i_offsets[jt+0xf] = jp + kf;
			}
			RADIX_64_DIT(a,dit_i_offsets,RE_IM_STRIDE, (double *)(t+l),t_offsets,1);
		}

		/*...and now do 64 radix-15 transforms.
		Using the initial out-index patterning:

			p40*[0x0,1,2,3,4,5,6,7,8,9,a,b,c,d,e] + plo[i], where i runs from 0-63
		=	phi[i/16] + p40*[0x0,1,2,3,4,5,6,7,8,9,a,b,c,d,e] + plo[i%16], for i =  0-63
		=	phi[i/16] + phi[4*j] + plo[i%16], for i =  0-63 and j = 0-14 for each i

		...the required output permutation is:

			000,0ff,1fe,2fd,03c,13b,23a,339,078,177,276,375,0b4,1b3,2b2,
			3b1,0f0,1ef,2ee,02d,12c,22b,32a,069,168,267,366,0a5,1a4,2a3,
			3a2,0e1,1e0,2df,01e,11d,21c,31b,05a,159,258,357,096,195,294,
			393,0d2,1d1,2d0,00f,10e,20d,30c,04b,14a,249,348,087,186,285,
			384,0c3,1c2,2c1,2c0,3bf,0fe,1fd,2fc,03b,13a,239,338,077,176,
			275,374,0b3,1b2,2b1,3b0,0ef,1ee,2ed,02c,12b,22a,329,068,167,
			266,365,0a4,1a3,2a2,3a1,0e0,1df,2de,01d,11c,21b,31a,059,158,
			257,356,095,194,293,392,0d1,1d0,2cf,00e,10d,20c,30b,04a,149,
			248,347,086,185,284,383,0c2,1c1,1c0,2bf,3be,0fd,1fc,2fb,03a,
			139,238,337,076,175,274,373,0b2,1b1,2b0,3af,0ee,1ed,2ec,02b,
			12a,229,328,067,166,265,364,0a3,1a2,2a1,3a0,0df,1de,2dd,01c,
			11b,21a,319,058,157,256,355,094,193,292,391,0d0,1cf,2ce,00d,
			10c,20b,30a,049,148,247,346,085,184,283,382,0c1,0c0,1bf,2be,
			3bd,0fc,1fb,2fa,039,138,237,336,075,174,273,372,0b1,1b0,2af,
			3ae,0ed,1ec,2eb,02a,129,228,327,066,165,264,363,0a2,1a1,2a0,
			39f,0de,1dd,2dc,01b,11a,219,318,057,156,255,354,093,192,291,
			390,0cf,1ce,2cd,00c,10b,20a,309,048,147,246,345,084,183,282,
			381,380,0bf,1be,2bd,3bc,0fb,1fa,2f9,038,137,236,335,074,173,
			272,371,0b0,1af,2ae,3ad,0ec,1eb,2ea,029,128,227,326,065,164,
			263,362,0a1,1a0,29f,39e,0dd,1dc,2db,01a,119,218,317,056,155,
			254,353,092,191,290,38f,0ce,1cd,2cc,00b,10a,209,308,047,146,
			245,344,083,182,281,280,37f,0be,1bd,2bc,3bb,0fa,1f9,2f8,037,
			136,235,334,073,172,271,370,0af,1ae,2ad,3ac,0eb,1ea,2e9,028,
			127,226,325,064,163,262,361,0a0,19f,29e,39d,0dc,1db,2da,019,
			118,217,316,055,154,253,352,091,190,28f,38e,0cd,1cc,2cb,00a,
			109,208,307,046,145,244,343,082,181,180,27f,37e,0bd,1bc,2bb,
			3ba,0f9,1f8,2f7,036,135,234,333,072,171,270,36f,0ae,1ad,2ac,
			3ab,0ea,1e9,2e8,027,126,225,324,063,162,261,360,09f,19e,29d,
			39c,0db,1da,2d9,018,117,216,315,054,153,252,351,090,18f,28e,
			38d,0cc,1cb,2ca,009,108,207,306,045,144,243,342,081,080,17f,
			27e,37d,0bc,1bb,2ba,3b9,0f8,1f7,2f6,035,134,233,332,071,170,
			26f,36e,0ad,1ac,2ab,3aa,0e9,1e8,2e7,026,125,224,323,062,161,
			260,35f,09e,19d,29c,39b,0da,1d9,2d8,017,116,215,314,053,152,
			251,350,08f,18e,28d,38c,0cb,1ca,2c9,008,107,206,305,044,143,
			242,341,340,07f,17e,27d,37c,0bb,1ba,2b9,3b8,0f7,1f6,2f5,034,
			133,232,331,070,16f,26e,36d,0ac,1ab,2aa,3a9,0e8,1e7,2e6,025,
			124,223,322,061,160,25f,35e,09d,19c,29b,39a,0d9,1d8,2d7,016,
			115,214,313,052,151,250,34f,08e,18d,28c,38b,0ca,1c9,2c8,007,
			106,205,304,043,142,241,240,33f,07e,17d,27c,37b,0ba,1b9,2b8,
			3b7,0f6,1f5,2f4,033,132,231,330,06f,16e,26d,36c,0ab,1aa,2a9,
			3a8,0e7,1e6,2e5,024,123,222,321,060,15f,25e,35d,09c,19b,29a,
			399,0d8,1d7,2d6,015,114,213,312,051,150,24f,34e,08d,18c,28b,
			38a,0c9,1c8,2c7,006,105,204,303,042,141,140,23f,33e,07d,17c,
			27b,37a,0b9,1b8,2b7,3b6,0f5,1f4,2f3,032,131,230,32f,06e,16d,
			26c,36b,0aa,1a9,2a8,3a7,0e6,1e5,2e4,023,122,221,320,05f,15e,
			25d,35c,09b,19a,299,398,0d7,1d6,2d5,014,113,212,311,050,14f,
			24e,34d,08c,18b,28a,389,0c8,1c7,2c6,005,104,203,302,041,040,
			13f,23e,33d,07c,17b,27a,379,0b8,1b7,2b6,3b5,0f4,1f3,2f2,031,
			130,22f,32e,06d,16c,26b,36a,0a9,1a8,2a7,3a6,0e5,1e4,2e3,022,
			121,220,31f,05e,15d,25c,35b,09a,199,298,397,0d6,1d5,2d4,013,
			112,211,310,04f,14e,24d,34c,08b,18a,289,388,0c7,1c6,2c5,004,
			103,202,301,300,03f,13e,23d,33c,07b,17a,279,378,0b7,1b6,2b5,
			3b4,0f3,1f2,2f1,030,12f,22e,32d,06c,16b,26a,369,0a8,1a7,2a6,
			3a5,0e4,1e3,2e2,021,120,21f,31e,05d,15c,25b,35a,099,198,297,
			396,0d5,1d4,2d3,012,111,210,30f,04e,14d,24c,34b,08a,189,288,
			387,0c6,1c5,2c4,003,102,201,200,2ff,03e,13d,23c,33b,07a,179,
			278,377,0b6,1b5,2b4,3b3,0f2,1f1,2f0,02f,12e,22d,32c,06b,16a,
			269,368,0a7,1a6,2a5,3a4,0e3,1e2,2e1,020,11f,21e,31d,05c,15b,
			25a,359,098,197,296,395,0d4,1d3,2d2,011,110,20f,30e,04d,14c,
			24b,34a,089,188,287,386,0c5,1c4,2c3,002,101,100,1ff,2fe,03d,
			13c,23b,33a,079,178,277,376,0b5,1b4,2b3,3b2,0f1,1f0,2ef,02e,
			12d,22c,32b,06a,169,268,367,0a6,1a5,2a4,3a3,0e2,1e1,2e0,01f,
			11e,21d,31c,05b,15a,259,358,097,196,295,394,0d3,1d2,2d1,010,
			10f,20e,30d,04c,14b,24a,349,088,187,286,385,0c4,1c3,2c2,001,

		but our initial patterning has unit strides going downward in each column, thus the real o-perm
		results from casting the above from a 64 x 15 matrix to a 15 x 64 [first convert to linear indexing,
		then chop into 64-element rows rather than 15-element ones]:

			000,0ff,1fe,2fd,03c,13b,23a,339,078,177,276,375,0b4,1b3,2b2,3b1,0f0,1ef,2ee,02d,12c,22b,32a,069,168,267,366,0a5,1a4,2a3,3a2,0e1,1e0,2df,01e,11d,21c,31b,05a,159,258,357,096,195,294,393,0d2,1d1,2d0,00f,10e,20d,30c,04b,14a,249,348,087,186,285,384,0c3,1c2,2c1,
			2c0,3bf,0fe,1fd,2fc,03b,13a,239,338,077,176,275,374,0b3,1b2,2b1,3b0,0ef,1ee,2ed,02c,12b,22a,329,068,167,266,365,0a4,1a3,2a2,3a1,0e0,1df,2de,01d,11c,21b,31a,059,158,257,356,095,194,293,392,0d1,1d0,2cf,00e,10d,20c,30b,04a,149,248,347,086,185,284,383,0c2,1c1,
			1c0,2bf,3be,0fd,1fc,2fb,03a,139,238,337,076,175,274,373,0b2,1b1,2b0,3af,0ee,1ed,2ec,02b,12a,229,328,067,166,265,364,0a3,1a2,2a1,3a0,0df,1de,2dd,01c,11b,21a,319,058,157,256,355,094,193,292,391,0d0,1cf,2ce,00d,10c,20b,30a,049,148,247,346,085,184,283,382,0c1,
			0c0,1bf,2be,3bd,0fc,1fb,2fa,039,138,237,336,075,174,273,372,0b1,1b0,2af,3ae,0ed,1ec,2eb,02a,129,228,327,066,165,264,363,0a2,1a1,2a0,39f,0de,1dd,2dc,01b,11a,219,318,057,156,255,354,093,192,291,390,0cf,1ce,2cd,00c,10b,20a,309,048,147,246,345,084,183,282,381,
			380,0bf,1be,2bd,3bc,0fb,1fa,2f9,038,137,236,335,074,173,272,371,0b0,1af,2ae,3ad,0ec,1eb,2ea,029,128,227,326,065,164,263,362,0a1,1a0,29f,39e,0dd,1dc,2db,01a,119,218,317,056,155,254,353,092,191,290,38f,0ce,1cd,2cc,00b,10a,209,308,047,146,245,344,083,182,281,
			280,37f,0be,1bd,2bc,3bb,0fa,1f9,2f8,037,136,235,334,073,172,271,370,0af,1ae,2ad,3ac,0eb,1ea,2e9,028,127,226,325,064,163,262,361,0a0,19f,29e,39d,0dc,1db,2da,019,118,217,316,055,154,253,352,091,190,28f,38e,0cd,1cc,2cb,00a,109,208,307,046,145,244,343,082,181,
			180,27f,37e,0bd,1bc,2bb,3ba,0f9,1f8,2f7,036,135,234,333,072,171,270,36f,0ae,1ad,2ac,3ab,0ea,1e9,2e8,027,126,225,324,063,162,261,360,09f,19e,29d,39c,0db,1da,2d9,018,117,216,315,054,153,252,351,090,18f,28e,38d,0cc,1cb,2ca,009,108,207,306,045,144,243,342,081,
			080,17f,27e,37d,0bc,1bb,2ba,3b9,0f8,1f7,2f6,035,134,233,332,071,170,26f,36e,0ad,1ac,2ab,3aa,0e9,1e8,2e7,026,125,224,323,062,161,260,35f,09e,19d,29c,39b,0da,1d9,2d8,017,116,215,314,053,152,251,350,08f,18e,28d,38c,0cb,1ca,2c9,008,107,206,305,044,143,242,341,
			340,07f,17e,27d,37c,0bb,1ba,2b9,3b8,0f7,1f6,2f5,034,133,232,331,070,16f,26e,36d,0ac,1ab,2aa,3a9,0e8,1e7,2e6,025,124,223,322,061,160,25f,35e,09d,19c,29b,39a,0d9,1d8,2d7,016,115,214,313,052,151,250,34f,08e,18d,28c,38b,0ca,1c9,2c8,007,106,205,304,043,142,241,
			240,33f,07e,17d,27c,37b,0ba,1b9,2b8,3b7,0f6,1f5,2f4,033,132,231,330,06f,16e,26d,36c,0ab,1aa,2a9,3a8,0e7,1e6,2e5,024,123,222,321,060,15f,25e,35d,09c,19b,29a,399,0d8,1d7,2d6,015,114,213,312,051,150,24f,34e,08d,18c,28b,38a,0c9,1c8,2c7,006,105,204,303,042,141,
			140,23f,33e,07d,17c,27b,37a,0b9,1b8,2b7,3b6,0f5,1f4,2f3,032,131,230,32f,06e,16d,26c,36b,0aa,1a9,2a8,3a7,0e6,1e5,2e4,023,122,221,320,05f,15e,25d,35c,09b,19a,299,398,0d7,1d6,2d5,014,113,212,311,050,14f,24e,34d,08c,18b,28a,389,0c8,1c7,2c6,005,104,203,302,041,
			040,13f,23e,33d,07c,17b,27a,379,0b8,1b7,2b6,3b5,0f4,1f3,2f2,031,130,22f,32e,06d,16c,26b,36a,0a9,1a8,2a7,3a6,0e5,1e4,2e3,022,121,220,31f,05e,15d,25c,35b,09a,199,298,397,0d6,1d5,2d4,013,112,211,310,04f,14e,24d,34c,08b,18a,289,388,0c7,1c6,2c5,004,103,202,301,
			300,03f,13e,23d,33c,07b,17a,279,378,0b7,1b6,2b5,3b4,0f3,1f2,2f1,030,12f,22e,32d,06c,16b,26a,369,0a8,1a7,2a6,3a5,0e4,1e3,2e2,021,120,21f,31e,05d,15c,25b,35a,099,198,297,396,0d5,1d4,2d3,012,111,210,30f,04e,14d,24c,34b,08a,189,288,387,0c6,1c5,2c4,003,102,201,
			200,2ff,03e,13d,23c,33b,07a,179,278,377,0b6,1b5,2b4,3b3,0f2,1f1,2f0,02f,12e,22d,32c,06b,16a,269,368,0a7,1a6,2a5,3a4,0e3,1e2,2e1,020,11f,21e,31d,05c,15b,25a,359,098,197,296,395,0d4,1d3,2d2,011,110,20f,30e,04d,14c,24b,34a,089,188,287,386,0c5,1c4,2c3,002,101,
			100,1ff,2fe,03d,13c,23b,33a,079,178,277,376,0b5,1b4,2b3,3b2,0f1,1f0,2ef,02e,12d,22c,32b,06a,169,268,367,0a6,1a5,2a4,3a3,0e2,1e1,2e0,01f,11e,21d,31c,05b,15a,259,358,097,196,295,394,0d3,1d2,2d1,010,10f,20e,30d,04c,14b,24a,349,088,187,286,385,0c4,1c3,2c2,001,

		...and transposing to yield a 64 x 15 perm-matrix with a simpler index pattern, to boot:

			000,2c0,1c0,0c0,380,280,180,080,340,240,140,040,300,200,100 + 0
			0f0,3b0,2b0,1b0,0b0,370,270,170,070,330,230,130,030,2f0,1f0 + f
			1f0,0f0,3b0,2b0,1b0,0b0,370,270,170,070,330,230,130,030,2f0 + e
			2f0,1f0,0f0,3b0,2b0,1b0,0b0,370,270,170,070,330,230,130,030 + d
			030,2f0,1f0,0f0,3b0,2b0,1b0,0b0,370,270,170,070,330,230,130 + c
			130,030,2f0,1f0,0f0,3b0,2b0,1b0,0b0,370,270,170,070,330,230 + b
			230,130,030,2f0,1f0,0f0,3b0,2b0,1b0,0b0,370,270,170,070,330 + a
			330,230,130,030,2f0,1f0,0f0,3b0,2b0,1b0,0b0,370,270,170,070 + 9
			070,330,230,130,030,2f0,1f0,0f0,3b0,2b0,1b0,0b0,370,270,170 + 8
			170,070,330,230,130,030,2f0,1f0,0f0,3b0,2b0,1b0,0b0,370,270 + 7
			270,170,070,330,230,130,030,2f0,1f0,0f0,3b0,2b0,1b0,0b0,370 + 6
			370,270,170,070,330,230,130,030,2f0,1f0,0f0,3b0,2b0,1b0,0b0 + 5
			0b0,370,270,170,070,330,230,130,030,2f0,1f0,0f0,3b0,2b0,1b0 + 4
			1b0,0b0,370,270,170,070,330,230,130,030,2f0,1f0,0f0,3b0,2b0 + 3
			2b0,1b0,0b0,370,270,170,070,330,230,130,030,2f0,1f0,0f0,3b0 + 2
			3b0,2b0,1b0,0b0,370,270,170,070,330,230,130,030,2f0,1f0,0f0 + 1
			0f0,3b0,2b0,1b0,0b0,370,270,170,070,330,230,130,030,2f0,1f0 + 0
		*	1e0,0e0,3a0,2a0,1a0,0a0,360,260,160,060,320,220,120,020,2e0 + f	<*** sub 0x10 from leading term whenever lo part passes thru 0!
			2e0,1e0,0e0,3a0,2a0,1a0,0a0,360,260,160,060,320,220,120,020 + e
			020,2e0,1e0,0e0,3a0,2a0,1a0,0a0,360,260,160,060,320,220,120 + d
			120,020,2e0,1e0,0e0,3a0,2a0,1a0,0a0,360,260,160,060,320,220 + c
			220,120,020,2e0,1e0,0e0,3a0,2a0,1a0,0a0,360,260,160,060,320 + b
			320,220,120,020,2e0,1e0,0e0,3a0,2a0,1a0,0a0,360,260,160,060 + a
			060,320,220,120,020,2e0,1e0,0e0,3a0,2a0,1a0,0a0,360,260,160 + 9
			160,060,320,220,120,020,2e0,1e0,0e0,3a0,2a0,1a0,0a0,360,260 + 8
			260,160,060,320,220,120,020,2e0,1e0,0e0,3a0,2a0,1a0,0a0,360 + 7
			360,260,160,060,320,220,120,020,2e0,1e0,0e0,3a0,2a0,1a0,0a0 + 6
			0a0,360,260,160,060,320,220,120,020,2e0,1e0,0e0,3a0,2a0,1a0 + 5
			1a0,0a0,360,260,160,060,320,220,120,020,2e0,1e0,0e0,3a0,2a0 + 4
			2a0,1a0,0a0,360,260,160,060,320,220,120,020,2e0,1e0,0e0,3a0 + 3
			3a0,2a0,1a0,0a0,360,260,160,060,320,220,120,020,2e0,1e0,0e0 + 2
			0e0,3a0,2a0,1a0,0a0,360,260,160,060,320,220,120,020,2e0,1e0 + 1
			1e0,0e0,3a0,2a0,1a0,0a0,360,260,160,060,320,220,120,020,2e0 + 0
		*	2d0,1d0,0d0,390,290,190,090,350,250,150,050,310,210,110,010 + f	<*** sub 0x10 from leading term whenever lo part passes thru 0!
			010,2d0,1d0,0d0,390,290,190,090,350,250,150,050,310,210,110 + e
			110,010,2d0,1d0,0d0,390,290,190,090,350,250,150,050,310,210 + d
			210,110,010,2d0,1d0,0d0,390,290,190,090,350,250,150,050,310 + c
			310,210,110,010,2d0,1d0,0d0,390,290,190,090,350,250,150,050 + b
			050,310,210,110,010,2d0,1d0,0d0,390,290,190,090,350,250,150 + a
			150,050,310,210,110,010,2d0,1d0,0d0,390,290,190,090,350,250 + 9
			250,150,050,310,210,110,010,2d0,1d0,0d0,390,290,190,090,350 + 8
			350,250,150,050,310,210,110,010,2d0,1d0,0d0,390,290,190,090 + 7
			090,350,250,150,050,310,210,110,010,2d0,1d0,0d0,390,290,190 + 6
			190,090,350,250,150,050,310,210,110,010,2d0,1d0,0d0,390,290 + 5
			290,190,090,350,250,150,050,310,210,110,010,2d0,1d0,0d0,390 + 4
			390,290,190,090,350,250,150,050,310,210,110,010,2d0,1d0,0d0 + 3
			0d0,390,290,190,090,350,250,150,050,310,210,110,010,2d0,1d0 + 2
			1d0,0d0,390,290,190,090,350,250,150,050,310,210,110,010,2d0 + 1
			2d0,1d0,0d0,390,290,190,090,350,250,150,050,310,210,110,010 + 0
		*	000,2c0,1c0,0c0,380,280,180,080,340,240,140,040,300,200,100 + f	<*** sub 0x10 from leading term whenever lo part passes thru 0!
			100,000,2c0,1c0,0c0,380,280,180,080,340,240,140,040,300,200 + e
			200,100,000,2c0,1c0,0c0,380,280,180,080,340,240,140,040,300 + d
			300,200,100,000,2c0,1c0,0c0,380,280,180,080,340,240,140,040 + c
			040,300,200,100,000,2c0,1c0,0c0,380,280,180,080,340,240,140 + b
			140,040,300,200,100,000,2c0,1c0,0c0,380,280,180,080,340,240 + a
			240,140,040,300,200,100,000,2c0,1c0,0c0,380,280,180,080,340 + 9
			340,240,140,040,300,200,100,000,2c0,1c0,0c0,380,280,180,080 + 8
			080,340,240,140,040,300,200,100,000,2c0,1c0,0c0,380,280,180 + 7
			180,080,340,240,140,040,300,200,100,000,2c0,1c0,0c0,380,280 + 6
			280,180,080,340,240,140,040,300,200,100,000,2c0,1c0,0c0,380 + 5
			380,280,180,080,340,240,140,040,300,200,100,000,2c0,1c0,0c0 + 4
			0c0,380,280,180,080,340,240,140,040,300,200,100,000,2c0,1c0 + 3
			1c0,0c0,380,280,180,080,340,240,140,040,300,200,100,000,2c0 + 2
			2c0,1c0,0c0,380,280,180,080,340,240,140,040,300,200,100,000 + 1

		[Compare to radix-240, which has
			00,e0,d0,c0,b0,a0,90,80,70,60,50,40,30,20,10 + 0
			00,e0,d0,c0,b0,a0,90,80,70,60,50,40,30,20,10 + f
			10,00,e0,d0,c0,b0,a0,90,80,70,60,50,40,30,20 + e
			20,10,00,e0,d0,c0,b0,a0,90,80,70,60,50,40,30 + d
			30,20,10,00,e0,d0,c0,b0,a0,90,80,70,60,50,40 + c
			40,30,20,10,00,e0,d0,c0,b0,a0,90,80,70,60,50 + b
			50,40,30,20,10,00,e0,d0,c0,b0,a0,90,80,70,60 + a
			60,50,40,30,20,10,00,e0,d0,c0,b0,a0,90,80,70 + 9
			70,60,50,40,30,20,10,00,e0,d0,c0,b0,a0,90,80 + 8
			80,70,60,50,40,30,20,10,00,e0,d0,c0,b0,a0,90 + 7
			90,80,70,60,50,40,30,20,10,00,e0,d0,c0,b0,a0 + 6
			a0,90,80,70,60,50,40,30,20,10,00,e0,d0,c0,b0 + 5
			b0,a0,90,80,70,60,50,40,30,20,10,00,e0,d0,c0 + 4
			c0,b0,a0,90,80,70,60,50,40,30,20,10,00,e0,d0 + 3
			d0,c0,b0,a0,90,80,70,60,50,40,30,20,10,00,e0 + 2
			e0,d0,c0,b0,a0,90,80,70,60,50,40,30,20,10,00 + 1 .
		]
		In order to make the above amenable to loop-based execution, we need to encode the indices
		to the left & right of the + in computable-index fashion. 3-part recipe;

			[0] RHS is simple: (64 - i)%16 for i = 0,..,63.

		For the LHS we need 2 ingredients:
			[1] A formula yielding the leftmost elt of each row. Init k = 0x3c0, decrementing by 0x2c0 each row
				(mod 960 = 0x3c0) and further -= 0x10 whenever lo part passes thru 0 gives the desired terms for i > 0;
				just need apply a bitmask which == 0 if i==0 and all-ones otherwise, thus k = k & (-(i > 0)) .
			[2] An efficient decrement-0x100 (mod 0x3c0) scheme to yield the remaining leading parts of each row's elements:
				k -= 0x100; k += (-(k < 0))&0x3c0;
		*/
		tptr = t;	k = 0x3c;	// Use k as idx into phi[], so drop the trailing hex 0
		for(l = 0; l < 64; l++) {
		#if 1
		// [0] here:
			jp = plo[(64 - l)&0xf];
		// [1] here:
			mask3 = (-(l > 0));
			k0 = k & mask3;
		// [2] here:
			// Remember: All indices /= 16 here:		Now get the resulting p* offsets:
			k1 = k0-0x10; k1 += (-(k1 < 0))&0x3c;		k0 = jp + phi[k0];
			k2 = k1-0x10; k2 += (-(k2 < 0))&0x3c;		k1 = jp + phi[k1];
			k3 = k2-0x10; k3 += (-(k3 < 0))&0x3c;		k2 = jp + phi[k2];
			k4 = k3-0x10; k4 += (-(k4 < 0))&0x3c;		k3 = jp + phi[k3];
			k5 = k4-0x10; k5 += (-(k5 < 0))&0x3c;		k4 = jp + phi[k4];
			k6 = k5-0x10; k6 += (-(k6 < 0))&0x3c;		k5 = jp + phi[k5];
			k7 = k6-0x10; k7 += (-(k7 < 0))&0x3c;		k6 = jp + phi[k6];
			k8 = k7-0x10; k8 += (-(k8 < 0))&0x3c;		k7 = jp + phi[k7];
			k9 = k8-0x10; k9 += (-(k9 < 0))&0x3c;		k8 = jp + phi[k8];
			ka = k9-0x10; ka += (-(ka < 0))&0x3c;		k9 = jp + phi[k9];
			kb = ka-0x10; kb += (-(kb < 0))&0x3c;		ka = jp + phi[ka];
			kc = kb-0x10; kc += (-(kc < 0))&0x3c;		kb = jp + phi[kb];
			kd = kc-0x10; kd += (-(kd < 0))&0x3c;		kc = jp + phi[kc];
			ke = kd-0x10; ke += (-(ke < 0))&0x3c;		kd = jp + phi[kd];
														ke = jp + phi[ke];
			// Set up for next loop execution:
			k -= 0x2c + ((l&0xf) == 0); //<*** further decr leading term whenever lo part passes thru 0
			k += (-(k < 0))&0x3c;
		#else
		// Initial ordered-outputs run computes phi[i/16] + phi[4*j] + plo[i%16], for i =  0-63 and j = 0-14 for each i
			jp = phi[l>>4] + plo[l&0xf];
			k0 = jp + phi[0x0<<2];
			k1 = jp + phi[0x1<<2];
			k2 = jp + phi[0x2<<2];
			k3 = jp + phi[0x3<<2];
			k4 = jp + phi[0x4<<2];
			k5 = jp + phi[0x5<<2];
			k6 = jp + phi[0x6<<2];
			k7 = jp + phi[0x7<<2];
			k8 = jp + phi[0x8<<2];
			k9 = jp + phi[0x9<<2];
			ka = jp + phi[0xa<<2];
			kb = jp + phi[0xb<<2];
			kc = jp + phi[0xc<<2];
			kd = jp + phi[0xd<<2];
			ke = jp + phi[0xe<<2];
		#endif
			jt = j1; jp = j2;
			RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,
				tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im,(tptr+0x5)->re,(tptr+0x5)->im,(tptr+0x6)->re,(tptr+0x6)->im,(tptr+0x7)->re,(tptr+0x7)->im,(tptr+0x8)->re,(tptr+0x8)->im,(tptr+0x9)->re,(tptr+0x9)->im,(tptr+0xa)->re,(tptr+0xa)->im,(tptr+0xb)->re,(tptr+0xb)->im,(tptr+0xc)->re,(tptr+0xc)->im,(tptr+0xd)->re,(tptr+0xd)->im,(tptr+0xe)->re,(tptr+0xe)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke]
			);	tptr += 0xf;
		}
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy960_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *add0, *addr,*addi;
		struct complex *tptr;
		const int pfetch_dist = PFETCH_DIST;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,pf,
			p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0,pb0,pc0,pd0,pe0,pf0,
			p100,p110,p120,p130,p140,p150,p160,p170,p180,p190,p1a0,p1b0,p1c0,p1d0,p1e0,p1f0,
			p200,p210,p220,p230,p240,p250,p260,p270,p280,p290,p2a0,p2b0,p2c0,p2d0,p2e0,p2f0,
			p300,p310,p320,p330,p340,p350,p360,p370,p380,p390,p3a0,p3b0;
		int poff[RADIX>>2];
	// DFT stuff:
		int kk,mask3, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf;
		int plo[16], phi[RADIX>>4], p_out_hi[RADIX>>4];
		int dif_o_offsets[64], dit_i_offsets[64], t_offsets[64];
		uint64 i64;
	// DIF:
		const uint64 dif_perm16[60] = {
			0x0123456789abcdefull,0x1032547698badcfeull,0x23106754ab98efdcull,0x32017645ba89fecdull,
			0x54763201dcfeba89ull,0xcdefab9823106754ull,0x76450123fecd89abull,0xefdc98ba10325476ull,
			0xab98efdc67541032ull,0xba89fecd76450123ull,0x98badcfe54763201ull,0x45672310cdefab98ull,
			0xdcfeba8932017645ull,0x23106754ab98efdcull,0xfecd89ab01234567ull,0x1032547698badcfeull,
			0x67541032efdc98baull,0x76450123fecd89abull,0x54763201dcfeba89ull,0xcdefab9823106754ull,
			0x89abcdef45672310ull,0x98badcfe54763201ull,0xab98efdc67541032ull,0xba89fecd76450123ull,
			0xefdc98ba10325476ull,0xfecd89ab01234567ull,0xdcfeba8932017645ull,0x23106754ab98efdcull,
			0x45672310cdefab98ull,0x54763201dcfeba89ull,0x67541032efdc98baull,0x76450123fecd89abull,
			0x32017645ba89fecdull,0xab98efdc67541032ull,0x89abcdef45672310ull,0x98badcfe54763201ull,
			0xcdefab9823106754ull,0xdcfeba8932017645ull,0xefdc98ba10325476ull,0xfecd89ab01234567ull,
			0xba89fecd76450123ull,0x67541032efdc98baull,0x45672310cdefab98ull,0x54763201dcfeba89ull,
			0x1032547698badcfeull,0x89abcdef45672310ull,0x32017645ba89fecdull,0xab98efdc67541032ull,
			0x76450123fecd89abull,0xefdc98ba10325476ull,0xcdefab9823106754ull,0xdcfeba8932017645ull,
			0x98badcfe54763201ull,0x45672310cdefab98ull,0xba89fecd76450123ull,0x67541032efdc98baull,
			0x23106754ab98efdcull,0x32017645ba89fecdull,0x1032547698badcfeull,0x89abcdef45672310ull
		};
	// DIT:
		// Low parts [p0-f] of input-index perms - 16 distinct ones, 4 lookups for each radix-64 DFT:
		const uint64 dit_perm16[16] = {
			0x01327654fedcba98ull,0xfedcba9876543210ull,0x5467102398abefcdull,0x98abefcd10236745ull,
			0xab89cdfe23014576ull,0x23014576cdfe89baull,0x76543210ba98dcefull,0x10236745efcdab89ull,
			0xcdfe89ba45760132ull,0xba98dcef32105467ull,0xefcdab8967452301ull,0x4576013289bafedcull,
			0x32105467dcef98abull,0x67452301ab89cdfeull,0x89bafedc01327654ull,0xdcef98ab54671023ull
		};
		// Pattern of access of above-encoded 16-perms, 4 for each DFT-64: Tempting to encoded via hex-hars,
		// but not worth it, just use byte array and get half the maximum possible compression ratio but simpler code:
		const uint8 dit_perm16_idx[60] = {
			0x0,0x1,0x1,0x1,0x2,0x2,0x2,0x3,0x4,0x5,0x4,0x5,0x6,0x6,0x6,
			0x6,0x3,0x3,0x3,0x7,0x5,0x8,0x5,0x8,0x9,0x9,0x9,0x9,0x7,0x7,
			0x7,0xa,0x8,0xb,0xb,0xb,0xc,0xc,0xc,0xc,0xa,0xd,0xa,0xd,0xb,
			0xe,0xe,0xe,0xf,0xf,0xf,0x2,0xd,0x4,0xd,0x4,0xe,0x0,0x0,0x0
		};
		// Pattern of access of phi[] elements for the above-encoded 16-perms, 4 for each DFT-64:
		const uint8 dit_perm16_phidx[60] = {
			0x00,0x01,0x03,0x02,0x09,0x08,0x0a,0x0b,0x06,0x07,0x04,0x05,0x17,0x16,0x15,
			0x14,0x11,0x10,0x12,0x13,0x0e,0x0f,0x0c,0x0d,0x1f,0x1e,0x1d,0x1c,0x19,0x18,
			0x1a,0x1b,0x20,0x21,0x23,0x22,0x27,0x26,0x25,0x24,0x2e,0x2f,0x2c,0x2d,0x28,
			0x29,0x2b,0x2a,0x39,0x38,0x3a,0x3b,0x36,0x37,0x34,0x35,0x30,0x31,0x33,0x32
		};

		int incr = 0,j,j1,j2,jt,jp,k,l,ntmp;
		// incr = Carry-chain wts-multipliers recurrence length, which must divide
		// RADIX/[n-wayness of carry macro], e.g. RADIX/[16|8|4] = 60|120|240 for avx512,avx,sse, respectively:
		const int incr_long = 15,incr_med =10,incr_short = 5;
	  // Have no specialized HIACC carry macro in USE_AVX512 and ARMv8 SIMD, use nonzero incr-value to differenetiate vs AVX/AVX2:
	  #if defined(USE_AVX512) || defined(USE_ARM_V8_SIMD)
		const int incr_hiacc = 3;
	  #else
		const int incr_hiacc = 0;
	  #endif
		// Fermat-mod: For AVX+, define carry-subchain length in terms of 2^nfold subchains:
		const int *inc_arr = 0x0;
	  #ifdef USE_AVX512
		// For nfold > 2,  RADIX/8 not divisible by 2^nfold, so use a more-general inner-loop scheme which can handle that:
		int nfold = USE_SHORT_CY_CHAIN + 2;
		const int nexec_long[] = {30,30,30,30}, nexec_med[] = {15,15,15,15,15,15,15,15}, nexec_short[] = {8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7}, nexec_hiacc[] = {4,4,4,3,4,4,4,3,4,4,4,3,4,4,4,3,4,4,4,3,4,4,4,3,4,4,4,3,4,4,4,3};
	  #elif defined(USE_AVX)
		// For nfold > 4,  RADIX/4 not divisible by 2^nfold, so use a more-general inner-loop scheme which can handle that:
		int nfold = USE_SHORT_CY_CHAIN + 2;
		const int nexec_long[] = {60,60,60,60}, nexec_med[] = {30,30,30,30,30,30,30,30}, nexec_short[] = {15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15}, nexec_hiacc[] = {8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7,8,7};
	  #endif
		// Allows cy-macro error data to be used to fiddle incr on the fly to a smaller, safer value if necessary
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
			if(USE_SHORT_CY_CHAIN == 0)
				incr = incr_long;
			else if(USE_SHORT_CY_CHAIN == 1)
				incr = incr_med;
			else if(USE_SHORT_CY_CHAIN == 2)
				incr = incr_short;
			else
				incr = incr_hiacc;
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
		double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
		double rt,it;

	#ifdef USE_SSE2

		const double crnd = 3.0*0x4000000*0x2000000;
		double *add1,*add2,*add3;
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *tmp,*tm0,*tm1,*tm2;	// utility ptrs
		int *itmp,*itm2;			// Pointer into the bjmodn array
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		vec_dbl *max_err, *sse2_rnd, *half_arr, *isrt2,*one,*two,
			*sse2_c3m1, *sse2_s, *sse2_cn1, *sse2_cn2, *sse2_ss3, *sse2_sn1, *sse2_sn2,	// Need radix-15 roots explicitly
			*r00,	// Head of RADIX*vec_cmplx-sized local store #1
			*s1p00,	// Head of RADIX*vec_cmplx-sized local store #2
			*x00,*x01,*x02,*x03,*x04,*x05,*x06,*x07,*x08,*x09,*x0a,*x0b,*x0c,*x0d,*x0e,	// Temps for 2x radix-15 DFT
			*y00,*y01,*y02,*y03,*y04,*y05,*y06,*y07,*y08,*y09,*y0a,*y0b,*y0c,*y0d,*y0e,
			*cy_r,*cy_i;	// Need RADIX slots for sse2 carries, RADIX/2 for avx
	  #ifdef USE_AVX
		vec_dbl *base_negacyclic_root;
	  #endif

		/* These are used in conjunction with the langth-ODD_RADIX arrays in the USE_SCALAR_CARRY #define below;
		In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
		int idx_offset, idx_incr;
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
	NDIVR <<= 4; p10 = NDIVR;	pf += ( (pf >> DAT_BITS) << PAD_BITS );
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
		pd0 = pc0 + NDIVR;		pc0 += ( (pc0 >> DAT_BITS) << PAD_BITS );
		pe0 = pd0 + NDIVR;		pd0 += ( (pd0 >> DAT_BITS) << PAD_BITS );
		pf0 = pe0 + NDIVR;		pe0 += ( (pe0 >> DAT_BITS) << PAD_BITS );
		p100= pf0 + NDIVR;		pf0 += ( (pf0 >> DAT_BITS) << PAD_BITS );
		p110= p100 + NDIVR;		p100+= ( (p100>> DAT_BITS) << PAD_BITS );
		p120= p110 + NDIVR;		p110+= ( (p110>> DAT_BITS) << PAD_BITS );
		p130= p120 + NDIVR;		p120+= ( (p120>> DAT_BITS) << PAD_BITS );
		p140= p130 + NDIVR;		p130+= ( (p130>> DAT_BITS) << PAD_BITS );
		p150= p140 + NDIVR;		p140+= ( (p140>> DAT_BITS) << PAD_BITS );
		p160= p150 + NDIVR;		p150+= ( (p150>> DAT_BITS) << PAD_BITS );
		p170= p160 + NDIVR;		p160+= ( (p160>> DAT_BITS) << PAD_BITS );
		p180= p170 + NDIVR;		p170+= ( (p170>> DAT_BITS) << PAD_BITS );
		p190= p180 + NDIVR;		p180+= ( (p180>> DAT_BITS) << PAD_BITS );
		p1a0= p190 + NDIVR;		p190+= ( (p190>> DAT_BITS) << PAD_BITS );
		p1b0= p1a0 + NDIVR;		p1a0+= ( (p1a0>> DAT_BITS) << PAD_BITS );
		p1c0= p1b0 + NDIVR;		p1b0+= ( (p1b0>> DAT_BITS) << PAD_BITS );
		p1d0= p1c0 + NDIVR;		p1c0+= ( (p1c0>> DAT_BITS) << PAD_BITS );
		p1e0= p1d0 + NDIVR;		p1d0+= ( (p1d0>> DAT_BITS) << PAD_BITS );
		p1f0= p1e0 + NDIVR;		p1e0+= ( (p1e0>> DAT_BITS) << PAD_BITS );
		p200= p1f0 + NDIVR;		p1f0+= ( (p1f0>> DAT_BITS) << PAD_BITS );
		p210= p200 + NDIVR;		p200+= ( (p200>> DAT_BITS) << PAD_BITS );
		p220= p210 + NDIVR;		p210+= ( (p210>> DAT_BITS) << PAD_BITS );
		p230= p220 + NDIVR;		p220+= ( (p220>> DAT_BITS) << PAD_BITS );
		p240= p230 + NDIVR;		p230+= ( (p230>> DAT_BITS) << PAD_BITS );
		p250= p240 + NDIVR;		p240+= ( (p240>> DAT_BITS) << PAD_BITS );
		p260= p250 + NDIVR;		p250+= ( (p250>> DAT_BITS) << PAD_BITS );
		p270= p260 + NDIVR;		p260+= ( (p260>> DAT_BITS) << PAD_BITS );
		p280= p270 + NDIVR;		p270+= ( (p270>> DAT_BITS) << PAD_BITS );
		p290= p280 + NDIVR;		p280+= ( (p280>> DAT_BITS) << PAD_BITS );
		p2a0= p290 + NDIVR;		p290+= ( (p290>> DAT_BITS) << PAD_BITS );
		p2b0= p2a0 + NDIVR;		p2a0+= ( (p2a0>> DAT_BITS) << PAD_BITS );
		p2c0= p2b0 + NDIVR;		p2b0+= ( (p2b0>> DAT_BITS) << PAD_BITS );
		p2d0= p2c0 + NDIVR;		p2c0+= ( (p2c0>> DAT_BITS) << PAD_BITS );
		p2e0= p2d0 + NDIVR;		p2d0+= ( (p2d0>> DAT_BITS) << PAD_BITS );
		p2f0= p2e0 + NDIVR;		p2e0+= ( (p2e0>> DAT_BITS) << PAD_BITS );
		p300= p2f0 + NDIVR;		p2f0+= ( (p2f0>> DAT_BITS) << PAD_BITS );
		p310= p300 + NDIVR;		p300+= ( (p300>> DAT_BITS) << PAD_BITS );
		p320= p310 + NDIVR;		p310+= ( (p310>> DAT_BITS) << PAD_BITS );
		p330= p320 + NDIVR;		p320+= ( (p320>> DAT_BITS) << PAD_BITS );
		p340= p330 + NDIVR;		p330+= ( (p330>> DAT_BITS) << PAD_BITS );
		p350= p340 + NDIVR;		p340+= ( (p340>> DAT_BITS) << PAD_BITS );
		p360= p350 + NDIVR;		p350+= ( (p350>> DAT_BITS) << PAD_BITS );
		p370= p360 + NDIVR;		p360+= ( (p360>> DAT_BITS) << PAD_BITS );
		p380= p370 + NDIVR;		p370+= ( (p370>> DAT_BITS) << PAD_BITS );
		p390= p380 + NDIVR;		p380+= ( (p380>> DAT_BITS) << PAD_BITS );
		p3a0= p390 + NDIVR;		p390+= ( (p390>> DAT_BITS) << PAD_BITS );
		p3b0= p3a0 + NDIVR;		p3a0+= ( (p3a0>> DAT_BITS) << PAD_BITS );
		NDIVR >>= 4;			p3b0+= ( (p3b0>> DAT_BITS) << PAD_BITS );
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
		for(l = 0; l < 64; l++) {
			poff[l+ 64] = poff[l] + p100;
			poff[l+128] = poff[l] + p200;
			if(l < 48) {
				poff[l+192] = poff[l] + p300;
			}
		}

		plo[ 0]= 0;plo[ 1]=p1;plo[ 2]=p2;plo[ 3]=p3;plo[ 4]=p4;plo[ 5]=p5;plo[ 6]=p6;plo[ 7]=p7;
		plo[ 8]=p8;plo[ 9]=p9;plo[10]=pa;plo[11]=pb;plo[12]=pc;plo[13]=pd;plo[14]=pe;plo[15]=pf;

		phi[0x00] =    0; phi[0x01] =  p10; phi[0x02] =  p20; phi[0x03] =  p30;
		phi[0x04] =  p40; phi[0x05] =  p50; phi[0x06] =  p60; phi[0x07] =  p70;
		phi[0x08] =  p80; phi[0x09] =  p90; phi[0x0a] =  pa0; phi[0x0b] =  pb0;
		phi[0x0c] =  pc0; phi[0x0d] =  pd0; phi[0x0e] =  pe0; phi[0x0f] =  pf0;
		phi[0x10] = p100; phi[0x11] = p110; phi[0x12] = p120; phi[0x13] = p130;
		phi[0x14] = p140; phi[0x15] = p150; phi[0x16] = p160; phi[0x17] = p170;
		phi[0x18] = p180; phi[0x19] = p190; phi[0x1a] = p1a0; phi[0x1b] = p1b0;
		phi[0x1c] = p1c0; phi[0x1d] = p1d0; phi[0x1e] = p1e0; phi[0x1f] = p1f0;
		phi[0x20] = p200; phi[0x21] = p210; phi[0x22] = p220; phi[0x23] = p230;
		phi[0x24] = p240; phi[0x25] = p250; phi[0x26] = p260; phi[0x27] = p270;
		phi[0x28] = p280; phi[0x29] = p290; phi[0x2a] = p2a0; phi[0x2b] = p2b0;
		phi[0x2c] = p2c0; phi[0x2d] = p2d0; phi[0x2e] = p2e0; phi[0x2f] = p2f0;
		phi[0x30] = p300; phi[0x31] = p310; phi[0x32] = p320; phi[0x33] = p330;
		phi[0x34] = p340; phi[0x35] = p350; phi[0x36] = p360; phi[0x37] = p370;
		phi[0x38] = p380; phi[0x39] = p390; phi[0x3a] = p3a0; phi[0x3b] = p3b0;

		p_out_hi[0x00] =    0; p_out_hi[0x01] =  p10; p_out_hi[0x02] =  p20; p_out_hi[0x03] =  p30;	//  p00
		p_out_hi[0x04] =  p10; p_out_hi[0x05] =    0; p_out_hi[0x06] =  p30; p_out_hi[0x07] =  p20;	//  p80
		p_out_hi[0x08] =  p20; p_out_hi[0x09] =  p30; p_out_hi[0x0a] =  p10; p_out_hi[0x0b] =    0;	//  p40
		p_out_hi[0x0c] =  p10; p_out_hi[0x0d] =    0; p_out_hi[0x0e] =  p30; p_out_hi[0x0f] =  p20;	// p380
		p_out_hi[0x10] =  p20; p_out_hi[0x11] =  p30; p_out_hi[0x12] =  p10; p_out_hi[0x13] =    0;	// p340
		p_out_hi[0x14] =    0; p_out_hi[0x15] =  p10; p_out_hi[0x16] =  p20; p_out_hi[0x17] =  p30;	// p300
		p_out_hi[0x18] =  p20; p_out_hi[0x19] =  p30; p_out_hi[0x1a] =  p10; p_out_hi[0x1b] =    0;	// p2c0
		p_out_hi[0x1c] =    0; p_out_hi[0x1d] =  p10; p_out_hi[0x1e] =  p20; p_out_hi[0x1f] =  p30;	// p280
		p_out_hi[0x20] =  p30; p_out_hi[0x21] =  p20; p_out_hi[0x22] =    0; p_out_hi[0x23] =  p10;	// p240
		p_out_hi[0x24] =    0; p_out_hi[0x25] =  p10; p_out_hi[0x26] =  p20; p_out_hi[0x27] =  p30;	// p200
		p_out_hi[0x28] =  p30; p_out_hi[0x29] =  p20; p_out_hi[0x2a] =    0; p_out_hi[0x2b] =  p10;	// p1c0
		p_out_hi[0x2c] =  p10; p_out_hi[0x2d] =    0; p_out_hi[0x2e] =  p30; p_out_hi[0x2f] =  p20;	// p180
		p_out_hi[0x30] =  p30; p_out_hi[0x31] =  p20; p_out_hi[0x32] =    0; p_out_hi[0x33] =  p10;	// p140
		p_out_hi[0x34] =  p10; p_out_hi[0x35] =    0; p_out_hi[0x36] =  p30; p_out_hi[0x37] =  p20;	// p100
		p_out_hi[0x38] =  p20; p_out_hi[0x39] =  p30; p_out_hi[0x3a] =  p10; p_out_hi[0x3b] =    0;	//  pc0

		for(l = 0; l < 64; l++) {
			t_offsets[l] = ((l<<4) - l);	// = l*15
			t_offsets[l] <<= 1;	// Need 2x incr due to the Radix-64 DFT casting t-array pointer to double.
								// In SIMD mode similarly need 2x due to vec_dbl type being used to access vec-complex data
		}

	#ifdef USE_SSE2
		uint32 nwt16 = nwt << L2_SZ_VD;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code

		tmp = thread_arg->r00;
		r00   = tmp;	tmp += 0x780;	// Head of RADIX*vec_cmplx-sized local store #1
		s1p00 = tmp;	tmp += 0x780;	// Head of RADIX*vec_cmplx-sized local store #2
		x00    = tmp + 0x00;
		x01    = tmp + 0x02;
		x02    = tmp + 0x04;
		x03    = tmp + 0x06;
		x04    = tmp + 0x08;
		x05    = tmp + 0x0a;
		x06    = tmp + 0x0c;
		x07    = tmp + 0x0e;
		x08    = tmp + 0x10;
		x09    = tmp + 0x12;
		x0a    = tmp + 0x14;
		x0b    = tmp + 0x16;
		x0c    = tmp + 0x18;
		x0d    = tmp + 0x1a;
		x0e    = tmp + 0x1c;
		tmp += 0x1e;
		y00    = tmp + 0x00;
		y01    = tmp + 0x02;
		y02    = tmp + 0x04;
		y03    = tmp + 0x06;
		y04    = tmp + 0x08;
		y05    = tmp + 0x0a;
		y06    = tmp + 0x0c;
		y07    = tmp + 0x0e;
		y08    = tmp + 0x10;
		y09    = tmp + 0x12;
		y0a    = tmp + 0x14;
		y0b    = tmp + 0x16;
		y0c    = tmp + 0x18;
		y0d    = tmp + 0x1a;
		y0e    = tmp + 0x1c;
		tmp += 0x1e;
		// DFT-15 roots needed explicitly, but throw in isrt2 as a pad-to-een and handy anchoring element:
		isrt2     = tmp + 0x00;
		sse2_c3m1 = tmp + 0x01;
		sse2_s    = tmp + 0x02;
		sse2_cn1  = tmp + 0x03;
		sse2_cn2  = tmp + 0x04;
		sse2_ss3  = tmp + 0x05;
		sse2_sn1  = tmp + 0x06;
		sse2_sn2  = tmp + 0x07;
		one       = tmp + 0x08;
		two       = tmp + 0x09;
		tmp += 0x0a;
	  #ifdef USE_AVX512
		cy_r = tmp;	cy_i = tmp+0x78;	tmp += 2*0x78;	// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod
	  #elif defined(USE_AVX)
		cy_r = tmp;	cy_i = tmp+0x0f0;	tmp += 2*0x0f0;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod in AVX mode
	  #else
		cy_r = tmp;	cy_i = tmp+0x1e0;	tmp += 2*0x1e0;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
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

		sign_mask = (uint64*)(r00 + radix960_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (#doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
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
		#include "radix960_main_carry_loop.h"

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
