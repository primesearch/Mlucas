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

#define RADIX 1008	// Use #define rather than const int to ensure it's really a compile-time const in the C sense
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

  // In hex, RADIX = 0x3f0, RADIX*4 = 0xfc0
  // For Mersenne-mod we need max(4*ODD_RADIX, (16 [SSE2] or 64 [AVX]) + 4) added slots for the half_arr lookup tables.
  // 4*ODD_RADIX = 252 = 0xfc here, trumps both elements of the (20,68), thus use extra slots in SSE2 mode.
  // In AVX mode.this gets supplemented by the extra storage we use for chained computation of the negacyclic-weights.
  // Add relevant number (half_arr_offset1008 + RADIX) to get required value of radix1008_creals_in_local_store:
  #ifdef USE_AVX512	// 0xfc fewer carry slots than AVX:
	const int half_arr_offset1008 = 0x10c4;	// + RADIX = 0x14b4; Used for thread local-storage-integrity checking
	const int radix1008_creals_in_local_store = 0x15b0;	// AVX+LOACC: (half_arr_offset1008 + RADIX) + 0xfc = 0x14b4 + 0xfc = 0x15b0 and round up to nearest multiple of 16
  #elif defined(USE_AVX)
	const int half_arr_offset1008 = 0x11c0;	// + RADIX = 0x15b0; Used for thread local-storage-integrity checking
   #if HIACC	// Need extra 4*RADIX = 0xfc0 vec_dbl slots in this mode
	const int radix1008_creals_in_local_store = 0x2670;	// AVX+HIACC: (half_arr_offset1008 + 5*RADIX) + 0xfc = 0x2570 + 0xfc = 0x266c and round up to nearest multiple of 16
   #else
	const int radix1008_creals_in_local_store = 0x16b0;	// AVX+LOACC: (half_arr_offset1008 + RADIX) + 0xfc = 0x15b0 + 0xfc = 0x16ac and round up to nearest multiple of 16
   #endif
  #else
	const int half_arr_offset1008 = 0x13b8;	// + RADIX = 0x17a8; Used for thread local-storage-integrity checking
const int radix1008_creals_in_local_store = 0x1ca0+0x1000;	// (half_arr_offset1008 + RADIX) + 0xfc = 0x1b98 + 0xfc = 0x1c94 and round up to nearest multiple of 16
  #endif

  #ifdef USE_AVX
	#include "radix1008_avx_negadwt_consts.h"
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

int radix1008_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-1008 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-1008 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix1008_ditN_cy_dif1";
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
	int kk, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf,po_kperm[16],*po_ptr = &(po_kperm[0]);
	static int plo[16], phi[ODD_RADIX], toff[ODD_RADIX];
	int ioff[ODD_RADIX];	// ioff elts get recomputed on-the-fly for each radix-63 call
	uint64 i64;
// DIF:
	// Need storage for circular-shifts perms of a basic 63-vector, with shift count in [0,60] that means 2*63-3 elts:
	const uint8 *iptr, dif_p10_cperms[128] = {	// Pad byte-array to next-higher 8-byte multiple
		0,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01,
		0,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x07,0x06,0x05,0x04
	},	dif_phi[64] = {	// Only 63 entries, but again pad out to nearest 8-multiple
		0,0x02,0x01,0x08,0x07,0x06,0x05,0x04,0x03,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09
	};
	// Low parts [p0-f] of output-index perms - need one uint64 (16 x 4-bits) for each radix-16 DFT:
	const uint64 dif16_oidx_lo[ODD_RADIX] = {
		0x0123456789abcdefull,0x54763201dcfeba89ull,0xab98efdc67541032ull,0x1032547698badcfeull,0xcdefab9823106754ull,
		0x76450123fecd89abull,0x89abcdef45672310ull,0x32017645ba89fecdull,0xefdc98ba10325476ull,0x76450123fecd89abull,
		0x98badcfe54763201ull,0x23106754ab98efdcull,0xefdc98ba10325476ull,0x45672310cdefab98ull,0xba89fecd76450123ull,
		0x0123456789abcdefull,0xdcfeba8932017645ull,0x67541032efdc98baull,0xba89fecd76450123ull,0x1032547698badcfeull,
		0xcdefab9823106754ull,0x67541032efdc98baull,0x89abcdef45672310ull,0x32017645ba89fecdull,0xfecd89ab01234567ull,
		0x54763201dcfeba89ull,0xab98efdc67541032ull,0x32017645ba89fecdull,0xefdc98ba10325476ull,0x45672310cdefab98ull,
		0xab98efdc67541032ull,0x0123456789abcdefull,0xdcfeba8932017645ull,0x76450123fecd89abull,0x98badcfe54763201ull,
		0x23106754ab98efdcull,0xdcfeba8932017645ull,0x67541032efdc98baull,0x89abcdef45672310ull,0x23106754ab98efdcull,
		0xfecd89ab01234567ull,0x54763201dcfeba89ull,0xba89fecd76450123ull,0x1032547698badcfeull,0xcdefab9823106754ull,
		0x54763201dcfeba89ull,0xab98efdc67541032ull,0x0123456789abcdefull,0xcdefab9823106754ull,0x76450123fecd89abull,
		0x98badcfe54763201ull,0x32017645ba89fecdull,0xefdc98ba10325476ull,0x45672310cdefab98ull,0x98badcfe54763201ull,
		0x23106754ab98efdcull,0xfecd89ab01234567ull,0x45672310cdefab98ull,0xba89fecd76450123ull,0x1032547698badcfeull,
		0xdcfeba8932017645ull,0x67541032efdc98baull,0x89abcdef45672310ull
	};
// DIT:
	// Need storage for circular-shifts perms of a basic 63-vector, with shift count in [0,15] that means 63+15 elts:
	const uint8 dit_p10_cperms[80] = {	// Pad byte-array to next-higher 8-byte multiple
		0,0x3b,0x37,0x33,0x2f,0x2b,0x27,0x23,0x1f,0x1b,0x17,0x13,0x0f,0x0b,0x07,0x03,0x3e,0x3a,0x36,0x32,0x2e,0x2a,0x26,0x22,0x1e,0x1a,0x16,0x12,0x0e,0x0a,0x06,0x02,0x3d,0x39,0x35,0x31,0x2d,0x29,0x25,0x21,0x1d,0x19,0x15,0x11,0x0d,0x09,0x05,0x01,0x3c,0x38,0x34,0x30,0x2c,0x28,0x24,0x20,0x1c,0x18,0x14,0x10,0x0c,0x08,0x04,
		0,0x3b,0x37,0x33,0x2f,0x2b,0x27,0x23,0x1f,0x1b,0x17,0x13,0x0f,0x0b,0x07
	},	dit_phi[64] = {	// Only 63 entries, but again pad out to nearest 8-multiple
		0,0x02,0x01,0x06,0x08,0x07,0x04,0x03,0x05,0x2f,0x2e,0x2d,0x35,0x34,0x33,0x30,0x32,0x31,0x1f,0x1e,0x20,0x1d,0x1c,0x1b,0x23,0x22,0x21,0x0f,0x11,0x10,0x0d,0x0c,0x0e,0x0b,0x0a,0x09,0x3e,0x3d,0x3c,0x39,0x3b,0x3a,0x37,0x36,0x38,0x26,0x25,0x24,0x2c,0x2b,0x2a,0x27,0x29,0x28,0x16,0x15,0x17,0x14,0x13,0x12,0x1a,0x19,0x18
	};
	// Low parts [p0-f] of input-index perms - need one uint64 (16 x 4-bits) for each radix-16 DFT:
	const uint64 dit16_iidx_lo[ODD_RADIX] = {
		0x01327654fedcba98ull,0x5467102398abefcdull,0xab89cdfe23014576ull,0x76543210ba98dcefull,0x10236745efcdab89ull,
		0xcdfe89ba45760132ull,0x32105467dcef98abull,0xefcdab8967452301ull,0x89bafedc01327654ull,0xfedcba9876543210ull,
		0x5467102398abefcdull,0xab89cdfe23014576ull,0xba98dcef32105467ull,0x10236745efcdab89ull,0xcdfe89ba45760132ull,
		0x32105467dcef98abull,0x67452301ab89cdfeull,0x89bafedc01327654ull,0xfedcba9876543210ull,0x5467102398abefcdull,
		0x23014576cdfe89baull,0xba98dcef32105467ull,0x10236745efcdab89ull,0xcdfe89ba45760132ull,0xdcef98ab54671023ull,
		0x67452301ab89cdfeull,0x89bafedc01327654ull,0xfedcba9876543210ull,0x98abefcd10236745ull,0x23014576cdfe89baull,
		0xba98dcef32105467ull,0x10236745efcdab89ull,0x4576013289bafedcull,0xdcef98ab54671023ull,0x67452301ab89cdfeull,
		0x89bafedc01327654ull,0x76543210ba98dcefull,0x98abefcd10236745ull,0x23014576cdfe89baull,0xba98dcef32105467ull,
		0xefcdab8967452301ull,0x4576013289bafedcull,0xdcef98ab54671023ull,0x67452301ab89cdfeull,0x01327654fedcba98ull,
		0x76543210ba98dcefull,0x98abefcd10236745ull,0x23014576cdfe89baull,0x32105467dcef98abull,0xefcdab8967452301ull,
		0x4576013289bafedcull,0xdcef98ab54671023ull,0xab89cdfe23014576ull,0x01327654fedcba98ull,0x76543210ba98dcefull,
		0x98abefcd10236745ull,0xcdfe89ba45760132ull,0x32105467dcef98abull,0xefcdab8967452301ull,0x4576013289bafedcull,
		0x5467102398abefcdull,0xab89cdfe23014576ull,0x01327654fedcba98ull
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
	double *add0,*add1,*add2,*add3,*add4,*add5,*add6,*add7,*add8,*add9,*adda,*addb,*addc,*addd,*adde,*addf;	/* Addresses into array sections */
  #endif

	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
	const double crnd = 3.0*0x4000000*0x2000000;
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	vec_dbl
		*tmp,*tm0,*tm1,*tm2;	// Non-static utility ptrs
	static vec_dbl *max_err, *sse2_rnd, *half_arr,
		*two,*one,*sqrt2,*isrt2,*cc0,*ss0,	// radix-16 DFT roots
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
	static task_control_t   task_control = {NULL, (void*)cy1008_process_chunk, NULL, 0x0};

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
		cslots_in_local_store = radix1008_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix1008_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = sc_ptr;	r00   = tmp;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x7e0;	s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
		tmp += 0x7e0;
		two       = tmp + 0x00;	// AVX+ versions of various DFT macros assume consts 2.0,1.0,isrt2 laid out thusly
		one       = tmp + 0x01;
		sqrt2     = tmp + 0x02;
		isrt2     = tmp + 0x03;
		cc0       = tmp + 0x04;
		ss0       = tmp + 0x05;
		tmp += 0x6;	// sc_ptr += 0xfc6
	  #ifdef USE_AVX512
		cy_r = tmp;	cy_i = tmp+0x07e;	tmp += 2*0x07e;	// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
	  #elif defined(USE_AVX)
		cy_r = tmp;	cy_i = tmp+0x0fc;	tmp += 2*0x0fc;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0x1f8 + 2 => sc_ptr += 0x11c0
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 68 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy_r = tmp;	cy_i = tmp+0x1f8;	tmp += 2*0x1f8;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0x3f0 + 2 => sc_ptr += 0x13b8
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*ODD_RADIX] x 16 for Fermat-mod */
	  #endif
		ASSERT(half_arr_offset1008 == (uint32)(half_arr-sc_ptr), "half_arr_offset mismatches actual!");
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
			j = (1<<(2*(L2_SZ_VD-2))) + 4;	// 16+4 for sse2, 64+4 for avx
		} else {
			j = ODD_RADIX<<2;				// 4*ODD_RADIX
		}
		ASSERT((radix1008_creals_in_local_store << L2_SZ_VD) >= ((intptr_t)half_arr - (intptr_t)r00) + (j << L2_SZ_VD), "radix1008_creals_in_local_store checksum failed!");

		// Roots for radix-16 DFTs:
		VEC_DBL_INIT(two  , 2.0  );	VEC_DBL_INIT(one  , 1.0  );
		VEC_DBL_INIT(sqrt2, SQRT2);	VEC_DBL_INIT(isrt2, ISRT2);
		VEC_DBL_INIT(cc0, c16  );	VEC_DBL_INIT(ss0, s16);	// Radix-16 DFT macros assume [sqrt2,isrt2,cc0,ss0] memory ordering
		VEC_DBL_INIT(sse2_rnd, crnd);		/* SSE2 math = 53-mantissa-bit IEEE double-float: */

	// Init-mode calls to these functions which maintain an internal local-alloc static store:
		thr_id = -1;	// Use this special thread id for any macro-required thread-local-data inits...
		SSE2_RADIX_63_DIF(CY_THREADS,thr_id, 0,0,0,0);
		SSE2_RADIX_63_DIT(CY_THREADS,thr_id, 0,0,0,0);
		thr_id = 0;	// ...then revert to 0.

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

		See the radix28 version of this routine for additional details.
		*/
		#if 0
		// Simple qfloat-based loop to crank out the roots which populate the radix1008_avx_negadwt_consts table:
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
		tmp64 = radix1008_avx_negadwt_consts[1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(  1*I*Pi/(2*RADIX)) = sin((RADIX-  1)*I*Pi/(2*RADIX)) */
		tmp64 = radix1008_avx_negadwt_consts[2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(  2*I*Pi/(2*RADIX)) = sin((RADIX-  2)*I*Pi/(2*RADIX)) */
		tmp64 = radix1008_avx_negadwt_consts[3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(  3*I*Pi/(2*RADIX)) = sin((RADIX-  3)*I*Pi/(2*RADIX)) */	tmp += 2;
		for(j = 4; j < RADIX; j += 4) {
			tmp64 = radix1008_avx_negadwt_consts[j+0];	tmp->d0 = tm2->d0 = *(double *)&tmp64;	tm2 -= 2;
			tmp64 = radix1008_avx_negadwt_consts[j+1];	tmp->d1 = tm2->d3 = *(double *)&tmp64;
			tmp64 = radix1008_avx_negadwt_consts[j+2];	tmp->d2 = tm2->d2 = *(double *)&tmp64;
			tmp64 = radix1008_avx_negadwt_consts[j+3];	tmp->d3 = tm2->d1 = *(double *)&tmp64;	tmp += 2;
		}
		tmp = base_negacyclic_root + RADIX*2;	// reset to point to start of above block
		nbytes = RADIX*SZ_VD/2;	// RADIX/4 AVX-register-sized complex data

	  #else	// HIACC = false:

		// lower-precision version, which yields slightly more roundoff error, but is simpler and more storage-compact.
	   #ifdef USE_AVX512

		// Init exp(j*I*Pi/2/RADIX), for j = 0-7:
		tmp = base_negacyclic_root + 16;	// First 16 slots reserved for Re/Im parts of the 8 base multipliers
		tmp->d0 = 1.0;
		tmp64 = radix1008_avx_negadwt_consts[      1];	tmp->d1 = *(double *)&tmp64;	// cos(01*I*Pi/(2*RADIX))
		tmp64 = radix1008_avx_negadwt_consts[      2];	tmp->d2 = *(double *)&tmp64;	// cos(02*I*Pi/(2*RADIX))
		tmp64 = radix1008_avx_negadwt_consts[      3];	tmp->d3 = *(double *)&tmp64;	// cos(03*I*Pi/(2*RADIX))
		tmp64 = radix1008_avx_negadwt_consts[      4];	tmp->d4 = *(double *)&tmp64;	// cos(04*I*Pi/(2*RADIX))
		tmp64 = radix1008_avx_negadwt_consts[      5];	tmp->d5 = *(double *)&tmp64;	// cos(05*I*Pi/(2*RADIX))
		tmp64 = radix1008_avx_negadwt_consts[      6];	tmp->d6 = *(double *)&tmp64;	// cos(06*I*Pi/(2*RADIX))
		tmp64 = radix1008_avx_negadwt_consts[      7];	tmp->d7 = *(double *)&tmp64;	// cos(07*I*Pi/(2*RADIX))
		++tmp;
		tmp->d0 = 0.0;
		tmp64 = radix1008_avx_negadwt_consts[RADIX-1];	tmp->d1 = *(double *)&tmp64;	// sin(01*I*Pi/(2*RADIX))
		tmp64 = radix1008_avx_negadwt_consts[RADIX-2];	tmp->d2 = *(double *)&tmp64;	// sin(02*I*Pi/(2*RADIX))
		tmp64 = radix1008_avx_negadwt_consts[RADIX-3];	tmp->d3 = *(double *)&tmp64;	// sin(03*I*Pi/(2*RADIX))
		tmp64 = radix1008_avx_negadwt_consts[RADIX-4];	tmp->d4 = *(double *)&tmp64;	// sin(04*I*Pi/(2*RADIX))
		tmp64 = radix1008_avx_negadwt_consts[RADIX-5];	tmp->d5 = *(double *)&tmp64;	// sin(05*I*Pi/(2*RADIX))
		tmp64 = radix1008_avx_negadwt_consts[RADIX-6];	tmp->d6 = *(double *)&tmp64;	// sin(06*I*Pi/(2*RADIX))
		tmp64 = radix1008_avx_negadwt_consts[RADIX-7];	tmp->d7 = *(double *)&tmp64;	// sin(07*I*Pi/(2*RADIX))
		++tmp;	// 0x480(base_negacyclic_root)
		tmp64 = radix1008_avx_negadwt_consts[      8];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// cos(08*I*Pi/(2*RADIX))
		++tmp;	// 0x4c0(base_negacyclic_root)
		tmp64 = radix1008_avx_negadwt_consts[RADIX-8];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// sin(08*I*Pi/(2*RADIX))
		tmp = base_negacyclic_root + 16;	// reset to point to start of above block

	   #elif defined(USE_AVX)

		// Init exp(j*I*Pi/2/RADIX), for j = 0-3:
		tmp = base_negacyclic_root + 8;	// First 8 slots reserved for Re/Im parts of the 4 base multipliers
		tmp->d0 = 1.0;
		tmp64 = radix1008_avx_negadwt_consts[      1];	tmp->d1 = *(double *)&tmp64;	// cos(01*I*Pi/(2*RADIX))
		tmp64 = radix1008_avx_negadwt_consts[      2];	tmp->d2 = *(double *)&tmp64;	// cos(02*I*Pi/(2*RADIX))
		tmp64 = radix1008_avx_negadwt_consts[      3];	tmp->d3 = *(double *)&tmp64;	// cos(03*I*Pi/(2*RADIX))

		(++tmp)->d0 = 0.0;
		tmp64 = radix1008_avx_negadwt_consts[RADIX-1];	tmp->d1 = *(double *)&tmp64;	// sin(01*I*Pi/(2*RADIX))
		tmp64 = radix1008_avx_negadwt_consts[RADIX-2];	tmp->d2 = *(double *)&tmp64;	// sin(02*I*Pi/(2*RADIX))
		tmp64 = radix1008_avx_negadwt_consts[RADIX-3];	tmp->d3 = *(double *)&tmp64;	// sin(03*I*Pi/(2*RADIX))
		++tmp;
		tmp64 = radix1008_avx_negadwt_consts[      4];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// cos(04*I*Pi/(2*RADIX))
		++tmp;
		tmp64 = radix1008_avx_negadwt_consts[RADIX-4];	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// sin(04*I*Pi/(2*RADIX))
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
		for(l = 0; l < 16; l++)
		{
			plo[l] = l*NDIVR;
			plo[l] += ( (plo[l] >> DAT_BITS) << PAD_BITS );
		}
		for(l = 0; l < ODD_RADIX; l++)
		{
			phi[l] = (l<<4)*NDIVR;
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

/*...The radix-1008 final DIT pass is here.	*/

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
		#include "radix1008_main_carry_loop.h"

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
		ASSERT(0x0 == cy1008_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
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

void radix1008_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-1008 = 63x16 complex DIF FFT pass on the data in the length-N real vector A.
*/
	int k,l, j,j1,j2,jp,jt, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf;
	static int NDIVR,first_entry=TRUE;
	struct complex t[RADIX], *tptr;
	// Need storage for circular-shifts perms of a basic 63-vector, with shift count in [0,60] that means 2*63-3 elts:
	const uint8 *iptr, dif_p10_cperms[128] = {	// Pad byte-array to next-higher 8-byte multiple
		0,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01,
		0,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x07,0x06,0x05,0x04
	},	dif_phi[64] = {	// Only 63 entries, but again pad out to nearest 8-multiple
		0,0x02,0x01,0x08,0x07,0x06,0x05,0x04,0x03,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09
	};
	static int plo[16], phi[ODD_RADIX], toff[ODD_RADIX];
	int ioff[ODD_RADIX];	// ioff elts get recomputed on-the-fly for each radix-63 call
	uint64 i64;
	// Low parts [p0-f] of output-index perms - need one uint64 (16 x 4-bits) for each radix-16 DFT:
	const uint64 dif16_oidx_lo[ODD_RADIX] = {
		0x0123456789abcdefull,0x54763201dcfeba89ull,0xab98efdc67541032ull,0x1032547698badcfeull,0xcdefab9823106754ull,
		0x76450123fecd89abull,0x89abcdef45672310ull,0x32017645ba89fecdull,0xefdc98ba10325476ull,0x76450123fecd89abull,
		0x98badcfe54763201ull,0x23106754ab98efdcull,0xefdc98ba10325476ull,0x45672310cdefab98ull,0xba89fecd76450123ull,
		0x0123456789abcdefull,0xdcfeba8932017645ull,0x67541032efdc98baull,0xba89fecd76450123ull,0x1032547698badcfeull,
		0xcdefab9823106754ull,0x67541032efdc98baull,0x89abcdef45672310ull,0x32017645ba89fecdull,0xfecd89ab01234567ull,
		0x54763201dcfeba89ull,0xab98efdc67541032ull,0x32017645ba89fecdull,0xefdc98ba10325476ull,0x45672310cdefab98ull,
		0xab98efdc67541032ull,0x0123456789abcdefull,0xdcfeba8932017645ull,0x76450123fecd89abull,0x98badcfe54763201ull,
		0x23106754ab98efdcull,0xdcfeba8932017645ull,0x67541032efdc98baull,0x89abcdef45672310ull,0x23106754ab98efdcull,
		0xfecd89ab01234567ull,0x54763201dcfeba89ull,0xba89fecd76450123ull,0x1032547698badcfeull,0xcdefab9823106754ull,
		0x54763201dcfeba89ull,0xab98efdc67541032ull,0x0123456789abcdefull,0xcdefab9823106754ull,0x76450123fecd89abull,
		0x98badcfe54763201ull,0x32017645ba89fecdull,0xefdc98ba10325476ull,0x45672310cdefab98ull,0x98badcfe54763201ull,
		0x23106754ab98efdcull,0xfecd89ab01234567ull,0x45672310cdefab98ull,0xba89fecd76450123ull,0x1032547698badcfeull,
		0xdcfeba8932017645ull,0x67541032efdc98baull,0x89abcdef45672310ull
	};

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
		for(l = 0; l < 16; l++)
		{
			plo[l] = l*NDIVR;
			plo[l] += ( (plo[l] >> DAT_BITS) << PAD_BITS );
		}
		for(l = 0; l < ODD_RADIX; l++)
		{
			phi[l] = (l<<4)*NDIVR;
			phi[l] += ( (phi[l] >> DAT_BITS) << PAD_BITS );
			toff[l] = l+l;
		}
	}

/*...The radix-1008 pass is here.	*/

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

	//...gather the needed data (1008 64-bit complex, i.e 2016 64-bit reals) and do 16 radix-63 transforms...
	/*
	Twiddleless version arranges 16 sets of radix-63 DIF inputs as follows: 0 in upper left corner,
	decrement 16 horizontally and 63 vertically, all resulting indexing modulo 1008 ( = 0x3f0):
	(we can auto-generate these by compiling test_fft_radix.c with -DTTYPE=0 -DRADIX=1008, running
	the resulting executable and snarfing the first set of index-outputs, "DIF/DIT input-scramble array").

	DIF/DIT input-scramble array = [in hex]
		000,3e0,3d0,3c0,3b0,3a0,390,380,370,360,350,340,330,320,310,300,2f0,2e0,2d0,2c0,2b0,2a0,290,280,270,260,250,240,230,220,210,200,1f0,1e0,1d0,1c0,1b0,1a0,190,180,170,160,150,140,130,120,110,100,0f0,0e0,0d0,0c0,0b0,0a0,090,080,070,060,050,040,030,020,010,
		3b1,3a1,391,381,371,361,351,341,331,321,311,301,2f1,2e1,2d1,2c1,2b1,2a1,291,281,271,261,251,241,231,221,211,201,1f1,1e1,1d1,1c1,1b1,1a1,191,181,171,161,151,141,131,121,111,101,0f1,0e1,0d1,0c1,0b1,0a1,091,081,071,061,051,041,031,021,011,001,3e1,3d1,3c1,
		372,362,352,342,332,322,312,302,2f2,2e2,2d2,2c2,2b2,2a2,292,282,272,262,252,242,232,222,212,202,1f2,1e2,1d2,1c2,1b2,1a2,192,182,172,162,152,142,132,122,112,102,0f2,0e2,0d2,0c2,0b2,0a2,092,082,072,062,052,042,032,022,012,002,3e2,3d2,3c2,3b2,3a2,392,382,
		333,323,313,303,2f3,2e3,2d3,2c3,2b3,2a3,293,283,273,263,253,243,233,223,213,203,1f3,1e3,1d3,1c3,1b3,1a3,193,183,173,163,153,143,133,123,113,103,0f3,0e3,0d3,0c3,0b3,0a3,093,083,073,063,053,043,033,023,013,003,3e3,3d3,3c3,3b3,3a3,393,383,373,363,353,343,
		2f4,2e4,2d4,2c4,2b4,2a4,294,284,274,264,254,244,234,224,214,204,1f4,1e4,1d4,1c4,1b4,1a4,194,184,174,164,154,144,134,124,114,104,0f4,0e4,0d4,0c4,0b4,0a4,094,084,074,064,054,044,034,024,014,004,3e4,3d4,3c4,3b4,3a4,394,384,374,364,354,344,334,324,314,304,
		2b5,2a5,295,285,275,265,255,245,235,225,215,205,1f5,1e5,1d5,1c5,1b5,1a5,195,185,175,165,155,145,135,125,115,105,0f5,0e5,0d5,0c5,0b5,0a5,095,085,075,065,055,045,035,025,015,005,3e5,3d5,3c5,3b5,3a5,395,385,375,365,355,345,335,325,315,305,2f5,2e5,2d5,2c5,
		276,266,256,246,236,226,216,206,1f6,1e6,1d6,1c6,1b6,1a6,196,186,176,166,156,146,136,126,116,106,0f6,0e6,0d6,0c6,0b6,0a6,096,086,076,066,056,046,036,026,016,006,3e6,3d6,3c6,3b6,3a6,396,386,376,366,356,346,336,326,316,306,2f6,2e6,2d6,2c6,2b6,2a6,296,286,
		237,227,217,207,1f7,1e7,1d7,1c7,1b7,1a7,197,187,177,167,157,147,137,127,117,107,0f7,0e7,0d7,0c7,0b7,0a7,097,087,077,067,057,047,037,027,017,007,3e7,3d7,3c7,3b7,3a7,397,387,377,367,357,347,337,327,317,307,2f7,2e7,2d7,2c7,2b7,2a7,297,287,277,267,257,247,
		1f8,1e8,1d8,1c8,1b8,1a8,198,188,178,168,158,148,138,128,118,108,0f8,0e8,0d8,0c8,0b8,0a8,098,088,078,068,058,048,038,028,018,008,3e8,3d8,3c8,3b8,3a8,398,388,378,368,358,348,338,328,318,308,2f8,2e8,2d8,2c8,2b8,2a8,298,288,278,268,258,248,238,228,218,208,
		1b9,1a9,199,189,179,169,159,149,139,129,119,109,0f9,0e9,0d9,0c9,0b9,0a9,099,089,079,069,059,049,039,029,019,009,3e9,3d9,3c9,3b9,3a9,399,389,379,369,359,349,339,329,319,309,2f9,2e9,2d9,2c9,2b9,2a9,299,289,279,269,259,249,239,229,219,209,1f9,1e9,1d9,1c9,
		17a,16a,15a,14a,13a,12a,11a,10a,0fa,0ea,0da,0ca,0ba,0aa,09a,08a,07a,06a,05a,04a,03a,02a,01a,00a,3ea,3da,3ca,3ba,3aa,39a,38a,37a,36a,35a,34a,33a,32a,31a,30a,2fa,2ea,2da,2ca,2ba,2aa,29a,28a,27a,26a,25a,24a,23a,22a,21a,20a,1fa,1ea,1da,1ca,1ba,1aa,19a,18a,
		13b,12b,11b,10b,0fb,0eb,0db,0cb,0bb,0ab,09b,08b,07b,06b,05b,04b,03b,02b,01b,00b,3eb,3db,3cb,3bb,3ab,39b,38b,37b,36b,35b,34b,33b,32b,31b,30b,2fb,2eb,2db,2cb,2bb,2ab,29b,28b,27b,26b,25b,24b,23b,22b,21b,20b,1fb,1eb,1db,1cb,1bb,1ab,19b,18b,17b,16b,15b,14b,
		0fc,0ec,0dc,0cc,0bc,0ac,09c,08c,07c,06c,05c,04c,03c,02c,01c,00c,3ec,3dc,3cc,3bc,3ac,39c,38c,37c,36c,35c,34c,33c,32c,31c,30c,2fc,2ec,2dc,2cc,2bc,2ac,29c,28c,27c,26c,25c,24c,23c,22c,21c,20c,1fc,1ec,1dc,1cc,1bc,1ac,19c,18c,17c,16c,15c,14c,13c,12c,11c,10c,
		0bd,0ad,09d,08d,07d,06d,05d,04d,03d,02d,01d,00d,3ed,3dd,3cd,3bd,3ad,39d,38d,37d,36d,35d,34d,33d,32d,31d,30d,2fd,2ed,2dd,2cd,2bd,2ad,29d,28d,27d,26d,25d,24d,23d,22d,21d,20d,1fd,1ed,1dd,1cd,1bd,1ad,19d,18d,17d,16d,15d,14d,13d,12d,11d,10d,0fd,0ed,0dd,0cd,
		07e,06e,05e,04e,03e,02e,01e,00e,3ee,3de,3ce,3be,3ae,39e,38e,37e,36e,35e,34e,33e,32e,31e,30e,2fe,2ee,2de,2ce,2be,2ae,29e,28e,27e,26e,25e,24e,23e,22e,21e,20e,1fe,1ee,1de,1ce,1be,1ae,19e,18e,17e,16e,15e,14e,13e,12e,11e,10e,0fe,0ee,0de,0ce,0be,0ae,09e,08e,
		03f,02f,01f,00f,3ef,3df,3cf,3bf,3af,39f,38f,37f,36f,35f,34f,33f,32f,31f,30f,2ff,2ef,2df,2cf,2bf,2af,29f,28f,27f,26f,25f,24f,23f,22f,21f,20f,1ff,1ef,1df,1cf,1bf,1af,19f,18f,17f,16f,15f,14f,13f,12f,11f,10f,0ff,0ef,0df,0cf,0bf,0af,09f,08f,07f,06f,05f,04f,

	If we subtract each row's common (mod 16) low p-part from all terms of the row
	as we do in the implementation to reduce the number of index offsets needing to be stored,
	we decrement 16 horizontally and 64 vertically, all (mod 1008), and where we elide the trailing 0 on all indices:

		00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01 + p0
		3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c + p1
		37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38 + p2
		33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34 + p3
		2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30 + p4
		2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c + p5
		27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28 + p6
		23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24 + p7
		1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20 + p8
		1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c + p9
		17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18 + pa
		13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14 + pb
		0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10 + pc
		0b,0a,09,08,07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c + pd
		07,06,05,04,03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08 + pe
		03,02,01,00,3e,3d,3c,3b,3a,39,38,37,36,35,34,33,32,31,30,2f,2e,2d,2c,2b,2a,29,28,27,26,25,24,23,22,21,20,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04 + pf

	In order to make the result amenable to loop-based execution, we need to encode the indices
	to the left & right of the + in easily-computable-index fashion. This scheme has 2 ingredients:
	[Use i as loop idx in comments, although use 'l' in actual code]

	[0] RHS is simply i = 0,..,f, i.e.plo[i]

	[1] Note that the (/= p10, where 10 = hex) row data are simply circular-shift perms of the basic (row 0) set.

	For row i (counting from 0 to 15), the leftward-shift count needing to be applied = 4*i.
	*/
		tptr = t;
		for(k = 0; k < 16; ++k)
		{
			iptr = dif_p10_cperms+(k<<2);
			for(l = 0; l < ODD_RADIX; l++)
			{
				ioff[l] = phi[iptr[l]];
			}
			RADIX_63_DIF(
				a+j1+plo[k], ioff, RE_IM_STRIDE,
				(double *)tptr, toff, 1
			);
			tptr += ODD_RADIX;
		}
	//...and now do 63 radix-16 transforms:
	/*
		O-perm in hex form - generated this with a view toward Impl 2-table-style streamlined indexing
		to match other large-radix DFTs in the compact-obj-code framework.
		The required output permutation is, grouping things into blocks of 16, one for each DFT output set:

		row:						o-Indices								...and in 2-table form:
		00	000,001,002,003,004,005,006,007,008,009,00a,00b,00c,00d,00e,00f   p[0123456789abcdef] + phi[00]
		01	025,024,027,026,023,022,020,021,02d,02c,02f,02e,02b,02a,028,029   p[54763201dcfeba89] + phi[02]
		02	01a,01b,019,018,01e,01f,01d,01c,016,017,015,014,011,010,013,012   p[ab98efdc67541032] + phi[01]
		03	081,080,083,082,085,084,087,086,089,088,08b,08a,08d,08c,08f,08e   p[1032547698badcfe] + phi[08]
		04	07c,07d,07e,07f,07a,07b,079,078,072,073,071,070,076,077,075,074   p[cdefab9823106754] + phi[07]
		05	067,066,064,065,060,061,062,063,06f,06e,06c,06d,068,069,06a,06b   p[76450123fecd89ab] + phi[06]
		06	058,059,05a,05b,05c,05d,05e,05f,054,055,056,057,052,053,051,050   p[89abcdef45672310] + phi[05]
		07	043,042,040,041,047,046,044,045,04b,04a,048,049,04f,04e,04c,04d   p[32017645ba89fecd] + phi[04]
		08	03e,03f,03d,03c,039,038,03b,03a,031,030,033,032,035,034,037,036   p[efdc98ba10325476] + phi[03]
		09	3e7,3e6,3e4,3e5,3e0,3e1,3e2,3e3,3ef,3ee,3ec,3ed,3e8,3e9,3ea,3eb   p[76450123fecd89ab] + phi[3e]
		0a	3d9,3d8,3db,3da,3dd,3dc,3df,3de,3d5,3d4,3d7,3d6,3d3,3d2,3d0,3d1   p[98badcfe54763201] + phi[3d]
		0b	3c2,3c3,3c1,3c0,3c6,3c7,3c5,3c4,3ca,3cb,3c9,3c8,3ce,3cf,3cd,3cc   p[23106754ab98efdc] + phi[3c]
		0c	3be,3bf,3bd,3bc,3b9,3b8,3bb,3ba,3b1,3b0,3b3,3b2,3b5,3b4,3b7,3b6   p[efdc98ba10325476] + phi[3b]
		0d	3a4,3a5,3a6,3a7,3a2,3a3,3a1,3a0,3ac,3ad,3ae,3af,3aa,3ab,3a9,3a8   p[45672310cdefab98] + phi[3a]
		0e	39b,39a,398,399,39f,39e,39c,39d,397,396,394,395,390,391,392,393   p[ba89fecd76450123] + phi[39]
		0f	380,381,382,383,384,385,386,387,388,389,38a,38b,38c,38d,38e,38f   p[0123456789abcdef] + phi[38]
		10	37d,37c,37f,37e,37b,37a,378,379,373,372,370,371,377,376,374,375   p[dcfeba8932017645] + phi[37]
		11	366,367,365,364,361,360,363,362,36e,36f,36d,36c,369,368,36b,36a   p[67541032efdc98ba] + phi[36]
		12	35b,35a,358,359,35f,35e,35c,35d,357,356,354,355,350,351,352,353   p[ba89fecd76450123] + phi[35]
		13	341,340,343,342,345,344,347,346,349,348,34b,34a,34d,34c,34f,34e   p[1032547698badcfe] + phi[34]
		14	33c,33d,33e,33f,33a,33b,339,338,332,333,331,330,336,337,335,334   p[cdefab9823106754] + phi[33]
		15	326,327,325,324,321,320,323,322,32e,32f,32d,32c,329,328,32b,32a   p[67541032efdc98ba] + phi[32]
		16	318,319,31a,31b,31c,31d,31e,31f,314,315,316,317,312,313,311,310   p[89abcdef45672310] + phi[31]
		17	303,302,300,301,307,306,304,305,30b,30a,308,309,30f,30e,30c,30d   p[32017645ba89fecd] + phi[30]
		18	2ff,2fe,2fc,2fd,2f8,2f9,2fa,2fb,2f0,2f1,2f2,2f3,2f4,2f5,2f6,2f7   p[fecd89ab01234567] + phi[2f]
		19	2e5,2e4,2e7,2e6,2e3,2e2,2e0,2e1,2ed,2ec,2ef,2ee,2eb,2ea,2e8,2e9   p[54763201dcfeba89] + phi[2e]
		1a	2da,2db,2d9,2d8,2de,2df,2dd,2dc,2d6,2d7,2d5,2d4,2d1,2d0,2d3,2d2   p[ab98efdc67541032] + phi[2d]
		1b	2c3,2c2,2c0,2c1,2c7,2c6,2c4,2c5,2cb,2ca,2c8,2c9,2cf,2ce,2cc,2cd   p[32017645ba89fecd] + phi[2c]
		1c	2be,2bf,2bd,2bc,2b9,2b8,2bb,2ba,2b1,2b0,2b3,2b2,2b5,2b4,2b7,2b6   p[efdc98ba10325476] + phi[2b]
		1d	2a4,2a5,2a6,2a7,2a2,2a3,2a1,2a0,2ac,2ad,2ae,2af,2aa,2ab,2a9,2a8   p[45672310cdefab98] + phi[2a]
		1e	29a,29b,299,298,29e,29f,29d,29c,296,297,295,294,291,290,293,292   p[ab98efdc67541032] + phi[29]
		1f	280,281,282,283,284,285,286,287,288,289,28a,28b,28c,28d,28e,28f = p[0123456789abcdef] + phi[28]
		20	27d,27c,27f,27e,27b,27a,278,279,273,272,270,271,277,276,274,275   p[dcfeba8932017645] + phi[27]
		21	267,266,264,265,260,261,262,263,26f,26e,26c,26d,268,269,26a,26b   p[76450123fecd89ab] + phi[26]
		22	259,258,25b,25a,25d,25c,25f,25e,255,254,257,256,253,252,250,251   p[98badcfe54763201] + phi[25]
		23	242,243,241,240,246,247,245,244,24a,24b,249,248,24e,24f,24d,24c   p[23106754ab98efdc] + phi[24]
		24	23d,23c,23f,23e,23b,23a,238,239,233,232,230,231,237,236,234,235   p[dcfeba8932017645] + phi[23]
		25	226,227,225,224,221,220,223,222,22e,22f,22d,22c,229,228,22b,22a   p[67541032efdc98ba] + phi[22]
		26	218,219,21a,21b,21c,21d,21e,21f,214,215,216,217,212,213,211,210   p[89abcdef45672310] + phi[21]
		27	202,203,201,200,206,207,205,204,20a,20b,209,208,20e,20f,20d,20c   p[23106754ab98efdc] + phi[20]
		28	1ff,1fe,1fc,1fd,1f8,1f9,1fa,1fb,1f0,1f1,1f2,1f3,1f4,1f5,1f6,1f7   p[fecd89ab01234567] + phi[1f]
		29	1e5,1e4,1e7,1e6,1e3,1e2,1e0,1e1,1ed,1ec,1ef,1ee,1eb,1ea,1e8,1e9   p[54763201dcfeba89] + phi[1e]
		2a	1db,1da,1d8,1d9,1df,1de,1dc,1dd,1d7,1d6,1d4,1d5,1d0,1d1,1d2,1d3   p[ba89fecd76450123] + phi[1d]
		2b	1c1,1c0,1c3,1c2,1c5,1c4,1c7,1c6,1c9,1c8,1cb,1ca,1cd,1cc,1cf,1ce   p[1032547698badcfe] + phi[1c]
		2c	1bc,1bd,1be,1bf,1ba,1bb,1b9,1b8,1b2,1b3,1b1,1b0,1b6,1b7,1b5,1b4   p[cdefab9823106754] + phi[1b]
		2d	1a5,1a4,1a7,1a6,1a3,1a2,1a0,1a1,1ad,1ac,1af,1ae,1ab,1aa,1a8,1a9   p[54763201dcfeba89] + phi[1a]
		2e	19a,19b,199,198,19e,19f,19d,19c,196,197,195,194,191,190,193,192   p[ab98efdc67541032] + phi[19]
		2f	180,181,182,183,184,185,186,187,188,189,18a,18b,18c,18d,18e,18f   p[0123456789abcdef] + phi[18]
		30	17c,17d,17e,17f,17a,17b,179,178,172,173,171,170,176,177,175,174   p[cdefab9823106754] + phi[17]
		31	167,166,164,165,160,161,162,163,16f,16e,16c,16d,168,169,16a,16b   p[76450123fecd89ab] + phi[16]
		32	159,158,15b,15a,15d,15c,15f,15e,155,154,157,156,153,152,150,151   p[98badcfe54763201] + phi[15]
		33	143,142,140,141,147,146,144,145,14b,14a,148,149,14f,14e,14c,14d   p[32017645ba89fecd] + phi[14]
		34	13e,13f,13d,13c,139,138,13b,13a,131,130,133,132,135,134,137,136   p[efdc98ba10325476] + phi[13]
		35	124,125,126,127,122,123,121,120,12c,12d,12e,12f,12a,12b,129,128   p[45672310cdefab98] + phi[12]
		36	119,118,11b,11a,11d,11c,11f,11e,115,114,117,116,113,112,110,111   p[98badcfe54763201] + phi[11]
		37	102,103,101,100,106,107,105,104,10a,10b,109,108,10e,10f,10d,10c   p[23106754ab98efdc] + phi[10]
		38	0ff,0fe,0fc,0fd,0f8,0f9,0fa,0fb,0f0,0f1,0f2,0f3,0f4,0f5,0f6,0f7   p[fecd89ab01234567] + phi[0f]
		39	0e4,0e5,0e6,0e7,0e2,0e3,0e1,0e0,0ec,0ed,0ee,0ef,0ea,0eb,0e9,0e8   p[45672310cdefab98] + phi[0e]
		3a	0db,0da,0d8,0d9,0df,0de,0dc,0dd,0d7,0d6,0d4,0d5,0d0,0d1,0d2,0d3   p[ba89fecd76450123] + phi[0d]
		3b	0c1,0c0,0c3,0c2,0c5,0c4,0c7,0c6,0c9,0c8,0cb,0ca,0cd,0cc,0cf,0ce   p[1032547698badcfe] + phi[0c]
		3c	0bd,0bc,0bf,0be,0bb,0ba,0b8,0b9,0b3,0b2,0b0,0b1,0b7,0b6,0b4,0b5   p[dcfeba8932017645] + phi[0b]
		3d	0a6,0a7,0a5,0a4,0a1,0a0,0a3,0a2,0ae,0af,0ad,0ac,0a9,0a8,0ab,0aa   p[67541032efdc98ba] + phi[0a]
		3e	098,099,09a,09b,09c,09d,09e,09f,094,095,096,097,092,093,091,090   p[89abcdef45672310] + phi[09]

		While the phi-indices may have a not-wildly-badly-exploitable pattern, simplest just to use a byte-array.
		*/
		tptr = t;
		for(l = 0; l < ODD_RADIX; ++l) {
		/* Initial-trial code for auto-gen of operm:
			k0 = plo[0x0];
			k1 = plo[0x1];
			k2 = plo[0x2];
			k3 = plo[0x3];
			k4 = plo[0x4];
			k5 = plo[0x5];
			k6 = plo[0x6];
			k7 = plo[0x7];
			k8 = plo[0x8];
			k9 = plo[0x9];
			ka = plo[0xa];
			kb = plo[0xb];
			kc = plo[0xc];
			kd = plo[0xd];
			ke = plo[0xe];
			kf = plo[0xf];
			jp = phi[l];
		*/
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
			jp = phi[dif_phi[l]];
			jt = j1 + jp; jp += j2;
			RADIX_16_DIF(
				tptr->re,tptr->im,(tptr+0x03f)->re,(tptr+0x03f)->im,(tptr+0x07e)->re,(tptr+0x07e)->im,(tptr+0x0bd)->re,(tptr+0x0bd)->im,(tptr+0x0fc)->re,(tptr+0x0fc)->im,(tptr+0x13b)->re,(tptr+0x13b)->im,(tptr+0x17a)->re,(tptr+0x17a)->im,(tptr+0x1b9)->re,(tptr+0x1b9)->im,(tptr+0x1f8)->re,(tptr+0x1f8)->im,(tptr+0x237)->re,(tptr+0x237)->im,(tptr+0x276)->re,(tptr+0x276)->im,(tptr+0x2b5)->re,(tptr+0x2b5)->im,(tptr+0x2f4)->re,(tptr+0x2f4)->im,(tptr+0x333)->re,(tptr+0x333)->im,(tptr+0x372)->re,(tptr+0x372)->im,(tptr+0x3b1)->re,(tptr+0x3b1)->im,
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],a[jt+kf],a[jp+kf],
				c16,s16);	tptr++;
		}
	}
}

/***************/

void radix1008_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-1008 = 63x16 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix1008_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
*/
	int k,l, j,j1,j2,jp,jt, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf;
	static int NDIVR,first_entry=TRUE;
	struct complex t[RADIX], *tptr;
	// Need storage for circular-shifts perms of a basic 63-vector, with shift count in [0,15] that means 63+15 elts:
	const uint8 *iptr, dit_p10_cperms[80] = {	// Pad byte-array to next-higher 8-byte multiple
		0,0x3b,0x37,0x33,0x2f,0x2b,0x27,0x23,0x1f,0x1b,0x17,0x13,0x0f,0x0b,0x07,0x03,0x3e,0x3a,0x36,0x32,0x2e,0x2a,0x26,0x22,0x1e,0x1a,0x16,0x12,0x0e,0x0a,0x06,0x02,0x3d,0x39,0x35,0x31,0x2d,0x29,0x25,0x21,0x1d,0x19,0x15,0x11,0x0d,0x09,0x05,0x01,0x3c,0x38,0x34,0x30,0x2c,0x28,0x24,0x20,0x1c,0x18,0x14,0x10,0x0c,0x08,0x04,
		0,0x3b,0x37,0x33,0x2f,0x2b,0x27,0x23,0x1f,0x1b,0x17,0x13,0x0f,0x0b,0x07
	},	dit_phi[64] = {	// Only 63 entries, but again pad out to nearest 8-multiple
		0,0x02,0x01,0x06,0x08,0x07,0x04,0x03,0x05,0x2f,0x2e,0x2d,0x35,0x34,0x33,0x30,0x32,0x31,0x1f,0x1e,0x20,0x1d,0x1c,0x1b,0x23,0x22,0x21,0x0f,0x11,0x10,0x0d,0x0c,0x0e,0x0b,0x0a,0x09,0x3e,0x3d,0x3c,0x39,0x3b,0x3a,0x37,0x36,0x38,0x26,0x25,0x24,0x2c,0x2b,0x2a,0x27,0x29,0x28,0x16,0x15,0x17,0x14,0x13,0x12,0x1a,0x19,0x18
	};
	static int plo[16], phi[ODD_RADIX], toff[ODD_RADIX];
	int ioff[ODD_RADIX];	// ioff elts get recomputed on-the-fly for each radix-63 call
	uint64 i64;
	// Low parts [p0-f] of input-index perms - need one uint64 (16 x 4-bits) for each radix-16 DFT:
	const uint64 dit16_iidx_lo[ODD_RADIX] = {
		0x01327654fedcba98ull,0x5467102398abefcdull,0xab89cdfe23014576ull,0x76543210ba98dcefull,0x10236745efcdab89ull,
		0xcdfe89ba45760132ull,0x32105467dcef98abull,0xefcdab8967452301ull,0x89bafedc01327654ull,0xfedcba9876543210ull,
		0x5467102398abefcdull,0xab89cdfe23014576ull,0xba98dcef32105467ull,0x10236745efcdab89ull,0xcdfe89ba45760132ull,
		0x32105467dcef98abull,0x67452301ab89cdfeull,0x89bafedc01327654ull,0xfedcba9876543210ull,0x5467102398abefcdull,
		0x23014576cdfe89baull,0xba98dcef32105467ull,0x10236745efcdab89ull,0xcdfe89ba45760132ull,0xdcef98ab54671023ull,
		0x67452301ab89cdfeull,0x89bafedc01327654ull,0xfedcba9876543210ull,0x98abefcd10236745ull,0x23014576cdfe89baull,
		0xba98dcef32105467ull,0x10236745efcdab89ull,0x4576013289bafedcull,0xdcef98ab54671023ull,0x67452301ab89cdfeull,
		0x89bafedc01327654ull,0x76543210ba98dcefull,0x98abefcd10236745ull,0x23014576cdfe89baull,0xba98dcef32105467ull,
		0xefcdab8967452301ull,0x4576013289bafedcull,0xdcef98ab54671023ull,0x67452301ab89cdfeull,0x01327654fedcba98ull,
		0x76543210ba98dcefull,0x98abefcd10236745ull,0x23014576cdfe89baull,0x32105467dcef98abull,0xefcdab8967452301ull,
		0x4576013289bafedcull,0xdcef98ab54671023ull,0xab89cdfe23014576ull,0x01327654fedcba98ull,0x76543210ba98dcefull,
		0x98abefcd10236745ull,0xcdfe89ba45760132ull,0x32105467dcef98abull,0xefcdab8967452301ull,0x4576013289bafedcull,
		0x5467102398abefcdull,0xab89cdfe23014576ull,0x01327654fedcba98ull
	};

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
		for(l = 0; l < 16; l++)
		{
			plo[l] = l*NDIVR;
			plo[l] += ( (plo[l] >> DAT_BITS) << PAD_BITS );
		}
		for(l = 0; l < ODD_RADIX; l++)
		{
			phi[l] = (l<<4)*NDIVR;
			phi[l] += ( (phi[l] >> DAT_BITS) << PAD_BITS );
			toff[l] = l+l;
		}
	}

/*...The radix-1008 pass is here.	*/

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
	Gather the needed data (1008 64-bit complex) and do 63 radix-16 transforms:

	Twiddleless version uses same linear-index-vector-form permutation as in DIF -
	Remember, inputs to DIT are bit-reversed, so use output of test_fft_radix() with
	TTYPE=0 to auto-generate needed input-index permutation:

	Combined DIT input-scramble array = [in hex]
		000,001,003,002,007,006,005,004,00f,00e,00d,00c,00b,00a,009,008   p[01327654fedcba98] + phi[00]
		025,024,026,027,021,020,022,023,029,028,02a,02b,02e,02f,02c,02d   p[5467102398abefcd] + phi[02]
		01a,01b,018,019,01c,01d,01f,01e,012,013,010,011,014,015,017,016   p[ab89cdfe23014576] + phi[01]
		067,066,065,064,063,062,061,060,06b,06a,069,068,06d,06c,06e,06f   p[76543210ba98dcef] + phi[06]
		081,080,082,083,086,087,084,085,08e,08f,08c,08d,08a,08b,088,089   p[10236745efcdab89] + phi[08]
		07c,07d,07f,07e,078,079,07b,07a,074,075,077,076,070,071,073,072   p[cdfe89ba45760132] + phi[07]
		043,042,041,040,045,044,046,047,04d,04c,04e,04f,049,048,04a,04b   p[32105467dcef98ab] + phi[04]
		03e,03f,03c,03d,03a,03b,038,039,036,037,034,035,032,033,030,031   p[efcdab8967452301] + phi[03]
		058,059,05b,05a,05f,05e,05d,05c,050,051,053,052,057,056,055,054   p[89bafedc01327654] + phi[05]
		2ff,2fe,2fd,2fc,2fb,2fa,2f9,2f8,2f7,2f6,2f5,2f4,2f3,2f2,2f1,2f0   p[fedcba9876543210] + phi[2f]
		2e5,2e4,2e6,2e7,2e1,2e0,2e2,2e3,2e9,2e8,2ea,2eb,2ee,2ef,2ec,2ed   p[5467102398abefcd] + phi[2e]
		2da,2db,2d8,2d9,2dc,2dd,2df,2de,2d2,2d3,2d0,2d1,2d4,2d5,2d7,2d6   p[ab89cdfe23014576] + phi[2d]
		35b,35a,359,358,35d,35c,35e,35f,353,352,351,350,355,354,356,357   p[ba98dcef32105467] + phi[35]
		341,340,342,343,346,347,344,345,34e,34f,34c,34d,34a,34b,348,349   p[10236745efcdab89] + phi[34]
		33c,33d,33f,33e,338,339,33b,33a,334,335,337,336,330,331,333,332   p[cdfe89ba45760132] + phi[33]
		303,302,301,300,305,304,306,307,30d,30c,30e,30f,309,308,30a,30b   p[32105467dcef98ab] + phi[30]
		326,327,324,325,322,323,320,321,32a,32b,328,329,32c,32d,32f,32e   p[67452301ab89cdfe] + phi[32]
		318,319,31b,31a,31f,31e,31d,31c,310,311,313,312,317,316,315,314   p[89bafedc01327654] + phi[31]
		1ff,1fe,1fd,1fc,1fb,1fa,1f9,1f8,1f7,1f6,1f5,1f4,1f3,1f2,1f1,1f0   p[fedcba9876543210] + phi[1f]
		1e5,1e4,1e6,1e7,1e1,1e0,1e2,1e3,1e9,1e8,1ea,1eb,1ee,1ef,1ec,1ed   p[5467102398abefcd] + phi[1e]
		202,203,200,201,204,205,207,206,20c,20d,20f,20e,208,209,20b,20a   p[23014576cdfe89ba] + phi[20]
		1db,1da,1d9,1d8,1dd,1dc,1de,1df,1d3,1d2,1d1,1d0,1d5,1d4,1d6,1d7   p[ba98dcef32105467] + phi[1d]
		1c1,1c0,1c2,1c3,1c6,1c7,1c4,1c5,1ce,1cf,1cc,1cd,1ca,1cb,1c8,1c9   p[10236745efcdab89] + phi[1c]
		1bc,1bd,1bf,1be,1b8,1b9,1bb,1ba,1b4,1b5,1b7,1b6,1b0,1b1,1b3,1b2   p[cdfe89ba45760132] + phi[1b]
		23d,23c,23e,23f,239,238,23a,23b,235,234,236,237,231,230,232,233   p[dcef98ab54671023] + phi[23]
		226,227,224,225,222,223,220,221,22a,22b,228,229,22c,22d,22f,22e   p[67452301ab89cdfe] + phi[22]
		218,219,21b,21a,21f,21e,21d,21c,210,211,213,212,217,216,215,214   p[89bafedc01327654] + phi[21]
		0ff,0fe,0fd,0fc,0fb,0fa,0f9,0f8,0f7,0f6,0f5,0f4,0f3,0f2,0f1,0f0   p[fedcba9876543210] + phi[0f]
		119,118,11a,11b,11e,11f,11c,11d,111,110,112,113,116,117,114,115   p[98abefcd10236745] + phi[11]
		102,103,100,101,104,105,107,106,10c,10d,10f,10e,108,109,10b,10a   p[23014576cdfe89ba] + phi[10]
		0db,0da,0d9,0d8,0dd,0dc,0de,0df,0d3,0d2,0d1,0d0,0d5,0d4,0d6,0d7   p[ba98dcef32105467] + phi[0d]
		0c1,0c0,0c2,0c3,0c6,0c7,0c4,0c5,0ce,0cf,0cc,0cd,0ca,0cb,0c8,0c9 = p[10236745efcdab89] + phi[0c]
		0e4,0e5,0e7,0e6,0e0,0e1,0e3,0e2,0e8,0e9,0eb,0ea,0ef,0ee,0ed,0ec   p[4576013289bafedc] + phi[0e]
		0bd,0bc,0be,0bf,0b9,0b8,0ba,0bb,0b5,0b4,0b6,0b7,0b1,0b0,0b2,0b3   p[dcef98ab54671023] + phi[0b]
		0a6,0a7,0a4,0a5,0a2,0a3,0a0,0a1,0aa,0ab,0a8,0a9,0ac,0ad,0af,0ae   p[67452301ab89cdfe] + phi[0a]
		098,099,09b,09a,09f,09e,09d,09c,090,091,093,092,097,096,095,094   p[89bafedc01327654] + phi[09]
		3e7,3e6,3e5,3e4,3e3,3e2,3e1,3e0,3eb,3ea,3e9,3e8,3ed,3ec,3ee,3ef   p[76543210ba98dcef] + phi[3e]
		3d9,3d8,3da,3db,3de,3df,3dc,3dd,3d1,3d0,3d2,3d3,3d6,3d7,3d4,3d5   p[98abefcd10236745] + phi[3d]
		3c2,3c3,3c0,3c1,3c4,3c5,3c7,3c6,3cc,3cd,3cf,3ce,3c8,3c9,3cb,3ca   p[23014576cdfe89ba] + phi[3c]
		39b,39a,399,398,39d,39c,39e,39f,393,392,391,390,395,394,396,397   p[ba98dcef32105467] + phi[39]
		3be,3bf,3bc,3bd,3ba,3bb,3b8,3b9,3b6,3b7,3b4,3b5,3b2,3b3,3b0,3b1   p[efcdab8967452301] + phi[3b]
		3a4,3a5,3a7,3a6,3a0,3a1,3a3,3a2,3a8,3a9,3ab,3aa,3af,3ae,3ad,3ac   p[4576013289bafedc] + phi[3a]
		37d,37c,37e,37f,379,378,37a,37b,375,374,376,377,371,370,372,373   p[dcef98ab54671023] + phi[37]
		366,367,364,365,362,363,360,361,36a,36b,368,369,36c,36d,36f,36e   p[67452301ab89cdfe] + phi[36]
		380,381,383,382,387,386,385,384,38f,38e,38d,38c,38b,38a,389,388   p[01327654fedcba98] + phi[38]
		267,266,265,264,263,262,261,260,26b,26a,269,268,26d,26c,26e,26f   p[76543210ba98dcef] + phi[26]
		259,258,25a,25b,25e,25f,25c,25d,251,250,252,253,256,257,254,255   p[98abefcd10236745] + phi[25]
		242,243,240,241,244,245,247,246,24c,24d,24f,24e,248,249,24b,24a   p[23014576cdfe89ba] + phi[24]
		2c3,2c2,2c1,2c0,2c5,2c4,2c6,2c7,2cd,2cc,2ce,2cf,2c9,2c8,2ca,2cb   p[32105467dcef98ab] + phi[2c]
		2be,2bf,2bc,2bd,2ba,2bb,2b8,2b9,2b6,2b7,2b4,2b5,2b2,2b3,2b0,2b1   p[efcdab8967452301] + phi[2b]
		2a4,2a5,2a7,2a6,2a0,2a1,2a3,2a2,2a8,2a9,2ab,2aa,2af,2ae,2ad,2ac   p[4576013289bafedc] + phi[2a]
		27d,27c,27e,27f,279,278,27a,27b,275,274,276,277,271,270,272,273   p[dcef98ab54671023] + phi[27]
		29a,29b,298,299,29c,29d,29f,29e,292,293,290,291,294,295,297,296   p[ab89cdfe23014576] + phi[29]
		280,281,283,282,287,286,285,284,28f,28e,28d,28c,28b,28a,289,288   p[01327654fedcba98] + phi[28]
		167,166,165,164,163,162,161,160,16b,16a,169,168,16d,16c,16e,16f   p[76543210ba98dcef] + phi[16]
		159,158,15a,15b,15e,15f,15c,15d,151,150,152,153,156,157,154,155   p[98abefcd10236745] + phi[15]
		17c,17d,17f,17e,178,179,17b,17a,174,175,177,176,170,171,173,172   p[cdfe89ba45760132] + phi[17]
		143,142,141,140,145,144,146,147,14d,14c,14e,14f,149,148,14a,14b   p[32105467dcef98ab] + phi[14]
		13e,13f,13c,13d,13a,13b,138,139,136,137,134,135,132,133,130,131   p[efcdab8967452301] + phi[13]
		124,125,127,126,120,121,123,122,128,129,12b,12a,12f,12e,12d,12c   p[4576013289bafedc] + phi[12]
		1a5,1a4,1a6,1a7,1a1,1a0,1a2,1a3,1a9,1a8,1aa,1ab,1ae,1af,1ac,1ad   p[5467102398abefcd] + phi[1a]
		19a,19b,198,199,19c,19d,19f,19e,192,193,190,191,194,195,197,196   p[ab89cdfe23014576] + phi[19]
		180,181,183,182,187,186,185,184,18f,18e,18d,18c,18b,18a,189,188   p[01327654fedcba98] + phi[18]
	*/
	/*...gather the needed data (1008 64-bit complex, i.e. 2016 64-bit reals) and do 63 radix-16 transforms...*/
		tptr = t;
		for(l = 0; l < ODD_RADIX; ++l) {
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
			jp = phi[dit_phi[l]];
			jt = j1 + jp; jp += j2;
			RADIX_16_DIT(
				a[jt+k0],a[jp+k0],a[jt+k1],a[jp+k1],a[jt+k2],a[jp+k2],a[jt+k3],a[jp+k3],a[jt+k4],a[jp+k4],a[jt+k5],a[jp+k5],a[jt+k6],a[jp+k6],a[jt+k7],a[jp+k7],a[jt+k8],a[jp+k8],a[jt+k9],a[jp+k9],a[jt+ka],a[jp+ka],a[jt+kb],a[jp+kb],a[jt+kc],a[jp+kc],a[jt+kd],a[jp+kd],a[jt+ke],a[jp+ke],a[jt+kf],a[jp+kf],
				tptr->re,tptr->im,(tptr+0x03f)->re,(tptr+0x03f)->im,(tptr+0x07e)->re,(tptr+0x07e)->im,(tptr+0x0bd)->re,(tptr+0x0bd)->im,(tptr+0x0fc)->re,(tptr+0x0fc)->im,(tptr+0x13b)->re,(tptr+0x13b)->im,(tptr+0x17a)->re,(tptr+0x17a)->im,(tptr+0x1b9)->re,(tptr+0x1b9)->im,(tptr+0x1f8)->re,(tptr+0x1f8)->im,(tptr+0x237)->re,(tptr+0x237)->im,(tptr+0x276)->re,(tptr+0x276)->im,(tptr+0x2b5)->re,(tptr+0x2b5)->im,(tptr+0x2f4)->re,(tptr+0x2f4)->im,(tptr+0x333)->re,(tptr+0x333)->im,(tptr+0x372)->re,(tptr+0x372)->im,(tptr+0x3b1)->re,(tptr+0x3b1)->im,
				c16,s16);	tptr++;
		}
		/*...and now do 16 radix-63 transforms. The required output permutation is as follows:

		000,3b0,370,330,2f0,2b0,270,230,1f0,1b0,170,130,0f0,0b0,070,030,3e0,3a0,360,320,2e0,2a0,260,220,1e0,1a0,160,120,0e0,0a0,060,020,3d0,390,350,310,2d0,290,250,210,1d0,190,150,110,0d0,090,050,010,3c0,380,340,300,2c0,280,240,200,1c0,180,140,100,0c0,080,040,
		03f,3ef,3af,36f,32f,2ef,2af,26f,22f,1ef,1af,16f,12f,0ef,0af,06f,02f,3df,39f,35f,31f,2df,29f,25f,21f,1df,19f,15f,11f,0df,09f,05f,01f,3cf,38f,34f,30f,2cf,28f,24f,20f,1cf,18f,14f,10f,0cf,08f,04f,00f,3bf,37f,33f,2ff,2bf,27f,23f,1ff,1bf,17f,13f,0ff,0bf,07f,
		07e,03e,3ee,3ae,36e,32e,2ee,2ae,26e,22e,1ee,1ae,16e,12e,0ee,0ae,06e,02e,3de,39e,35e,31e,2de,29e,25e,21e,1de,19e,15e,11e,0de,09e,05e,01e,3ce,38e,34e,30e,2ce,28e,24e,20e,1ce,18e,14e,10e,0ce,08e,04e,00e,3be,37e,33e,2fe,2be,27e,23e,1fe,1be,17e,13e,0fe,0be,
		0bd,07d,03d,3ed,3ad,36d,32d,2ed,2ad,26d,22d,1ed,1ad,16d,12d,0ed,0ad,06d,02d,3dd,39d,35d,31d,2dd,29d,25d,21d,1dd,19d,15d,11d,0dd,09d,05d,01d,3cd,38d,34d,30d,2cd,28d,24d,20d,1cd,18d,14d,10d,0cd,08d,04d,00d,3bd,37d,33d,2fd,2bd,27d,23d,1fd,1bd,17d,13d,0fd,
		0fc,0bc,07c,03c,3ec,3ac,36c,32c,2ec,2ac,26c,22c,1ec,1ac,16c,12c,0ec,0ac,06c,02c,3dc,39c,35c,31c,2dc,29c,25c,21c,1dc,19c,15c,11c,0dc,09c,05c,01c,3cc,38c,34c,30c,2cc,28c,24c,20c,1cc,18c,14c,10c,0cc,08c,04c,00c,3bc,37c,33c,2fc,2bc,27c,23c,1fc,1bc,17c,13c,
		13b,0fb,0bb,07b,03b,3eb,3ab,36b,32b,2eb,2ab,26b,22b,1eb,1ab,16b,12b,0eb,0ab,06b,02b,3db,39b,35b,31b,2db,29b,25b,21b,1db,19b,15b,11b,0db,09b,05b,01b,3cb,38b,34b,30b,2cb,28b,24b,20b,1cb,18b,14b,10b,0cb,08b,04b,00b,3bb,37b,33b,2fb,2bb,27b,23b,1fb,1bb,17b,
		17a,13a,0fa,0ba,07a,03a,3ea,3aa,36a,32a,2ea,2aa,26a,22a,1ea,1aa,16a,12a,0ea,0aa,06a,02a,3da,39a,35a,31a,2da,29a,25a,21a,1da,19a,15a,11a,0da,09a,05a,01a,3ca,38a,34a,30a,2ca,28a,24a,20a,1ca,18a,14a,10a,0ca,08a,04a,00a,3ba,37a,33a,2fa,2ba,27a,23a,1fa,1ba,
		1b9,179,139,0f9,0b9,079,039,3e9,3a9,369,329,2e9,2a9,269,229,1e9,1a9,169,129,0e9,0a9,069,029,3d9,399,359,319,2d9,299,259,219,1d9,199,159,119,0d9,099,059,019,3c9,389,349,309,2c9,289,249,209,1c9,189,149,109,0c9,089,049,009,3b9,379,339,2f9,2b9,279,239,1f9,
		1f8,1b8,178,138,0f8,0b8,078,038,3e8,3a8,368,328,2e8,2a8,268,228,1e8,1a8,168,128,0e8,0a8,068,028,3d8,398,358,318,2d8,298,258,218,1d8,198,158,118,0d8,098,058,018,3c8,388,348,308,2c8,288,248,208,1c8,188,148,108,0c8,088,048,008,3b8,378,338,2f8,2b8,278,238,
		237,1f7,1b7,177,137,0f7,0b7,077,037,3e7,3a7,367,327,2e7,2a7,267,227,1e7,1a7,167,127,0e7,0a7,067,027,3d7,397,357,317,2d7,297,257,217,1d7,197,157,117,0d7,097,057,017,3c7,387,347,307,2c7,287,247,207,1c7,187,147,107,0c7,087,047,007,3b7,377,337,2f7,2b7,277,
		276,236,1f6,1b6,176,136,0f6,0b6,076,036,3e6,3a6,366,326,2e6,2a6,266,226,1e6,1a6,166,126,0e6,0a6,066,026,3d6,396,356,316,2d6,296,256,216,1d6,196,156,116,0d6,096,056,016,3c6,386,346,306,2c6,286,246,206,1c6,186,146,106,0c6,086,046,006,3b6,376,336,2f6,2b6,
		2b5,275,235,1f5,1b5,175,135,0f5,0b5,075,035,3e5,3a5,365,325,2e5,2a5,265,225,1e5,1a5,165,125,0e5,0a5,065,025,3d5,395,355,315,2d5,295,255,215,1d5,195,155,115,0d5,095,055,015,3c5,385,345,305,2c5,285,245,205,1c5,185,145,105,0c5,085,045,005,3b5,375,335,2f5,
		2f4,2b4,274,234,1f4,1b4,174,134,0f4,0b4,074,034,3e4,3a4,364,324,2e4,2a4,264,224,1e4,1a4,164,124,0e4,0a4,064,024,3d4,394,354,314,2d4,294,254,214,1d4,194,154,114,0d4,094,054,014,3c4,384,344,304,2c4,284,244,204,1c4,184,144,104,0c4,084,044,004,3b4,374,334,
		333,2f3,2b3,273,233,1f3,1b3,173,133,0f3,0b3,073,033,3e3,3a3,363,323,2e3,2a3,263,223,1e3,1a3,163,123,0e3,0a3,063,023,3d3,393,353,313,2d3,293,253,213,1d3,193,153,113,0d3,093,053,013,3c3,383,343,303,2c3,283,243,203,1c3,183,143,103,0c3,083,043,003,3b3,373,
		372,332,2f2,2b2,272,232,1f2,1b2,172,132,0f2,0b2,072,032,3e2,3a2,362,322,2e2,2a2,262,222,1e2,1a2,162,122,0e2,0a2,062,022,3d2,392,352,312,2d2,292,252,212,1d2,192,152,112,0d2,092,052,012,3c2,382,342,302,2c2,282,242,202,1c2,182,142,102,0c2,082,042,002,3b2,
		3b1,371,331,2f1,2b1,271,231,1f1,1b1,171,131,0f1,0b1,071,031,3e1,3a1,361,321,2e1,2a1,261,221,1e1,1a1,161,121,0e1,0a1,061,021,3d1,391,351,311,2d1,291,251,211,1d1,191,151,111,0d1,091,051,011,3c1,381,341,301,2c1,281,241,201,1c1,181,141,101,0c1,081,041,001,

	If we subtract each row's common (mod 16) low p-part from all terms of the row
	as we do in the implementation to reduce the number of index offsets needing to be stored,
	we decrement 16 horizontally and 64 vertically, all (mod 1008), and where we elide the trailing 0 on all indices:

		00,3b,37,33,2f,2b,27,23,1f,1b,17,13,0f,0b,07,03,3e,3a,36,32,2e,2a,26,22,1e,1a,16,12,0e,0a,06,02,3d,39,35,31,2d,29,25,21,1d,19,15,11,0d,09,05,01,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04 + p0
		03,3e,3a,36,32,2e,2a,26,22,1e,1a,16,12,0e,0a,06,02,3d,39,35,31,2d,29,25,21,1d,19,15,11,0d,09,05,01,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,3b,37,33,2f,2b,27,23,1f,1b,17,13,0f,0b,07 + pf
		07,03,3e,3a,36,32,2e,2a,26,22,1e,1a,16,12,0e,0a,06,02,3d,39,35,31,2d,29,25,21,1d,19,15,11,0d,09,05,01,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,3b,37,33,2f,2b,27,23,1f,1b,17,13,0f,0b + pe
		0b,07,03,3e,3a,36,32,2e,2a,26,22,1e,1a,16,12,0e,0a,06,02,3d,39,35,31,2d,29,25,21,1d,19,15,11,0d,09,05,01,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,3b,37,33,2f,2b,27,23,1f,1b,17,13,0f + pd
		0f,0b,07,03,3e,3a,36,32,2e,2a,26,22,1e,1a,16,12,0e,0a,06,02,3d,39,35,31,2d,29,25,21,1d,19,15,11,0d,09,05,01,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,3b,37,33,2f,2b,27,23,1f,1b,17,13 + pc
		13,0f,0b,07,03,3e,3a,36,32,2e,2a,26,22,1e,1a,16,12,0e,0a,06,02,3d,39,35,31,2d,29,25,21,1d,19,15,11,0d,09,05,01,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,3b,37,33,2f,2b,27,23,1f,1b,17 + pb
		17,13,0f,0b,07,03,3e,3a,36,32,2e,2a,26,22,1e,1a,16,12,0e,0a,06,02,3d,39,35,31,2d,29,25,21,1d,19,15,11,0d,09,05,01,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,3b,37,33,2f,2b,27,23,1f,1b + pa
		1b,17,13,0f,0b,07,03,3e,3a,36,32,2e,2a,26,22,1e,1a,16,12,0e,0a,06,02,3d,39,35,31,2d,29,25,21,1d,19,15,11,0d,09,05,01,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,3b,37,33,2f,2b,27,23,1f + p9
		1f,1b,17,13,0f,0b,07,03,3e,3a,36,32,2e,2a,26,22,1e,1a,16,12,0e,0a,06,02,3d,39,35,31,2d,29,25,21,1d,19,15,11,0d,09,05,01,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,3b,37,33,2f,2b,27,23 + p8
		23,1f,1b,17,13,0f,0b,07,03,3e,3a,36,32,2e,2a,26,22,1e,1a,16,12,0e,0a,06,02,3d,39,35,31,2d,29,25,21,1d,19,15,11,0d,09,05,01,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,3b,37,33,2f,2b,27 + p7
		27,23,1f,1b,17,13,0f,0b,07,03,3e,3a,36,32,2e,2a,26,22,1e,1a,16,12,0e,0a,06,02,3d,39,35,31,2d,29,25,21,1d,19,15,11,0d,09,05,01,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,3b,37,33,2f,2b + p6
		2b,27,23,1f,1b,17,13,0f,0b,07,03,3e,3a,36,32,2e,2a,26,22,1e,1a,16,12,0e,0a,06,02,3d,39,35,31,2d,29,25,21,1d,19,15,11,0d,09,05,01,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,3b,37,33,2f + p5
		2f,2b,27,23,1f,1b,17,13,0f,0b,07,03,3e,3a,36,32,2e,2a,26,22,1e,1a,16,12,0e,0a,06,02,3d,39,35,31,2d,29,25,21,1d,19,15,11,0d,09,05,01,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,3b,37,33 + p4
		33,2f,2b,27,23,1f,1b,17,13,0f,0b,07,03,3e,3a,36,32,2e,2a,26,22,1e,1a,16,12,0e,0a,06,02,3d,39,35,31,2d,29,25,21,1d,19,15,11,0d,09,05,01,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,3b,37 + p3
		37,33,2f,2b,27,23,1f,1b,17,13,0f,0b,07,03,3e,3a,36,32,2e,2a,26,22,1e,1a,16,12,0e,0a,06,02,3d,39,35,31,2d,29,25,21,1d,19,15,11,0d,09,05,01,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00,3b + p2
		3b,37,33,2f,2b,27,23,1f,1b,17,13,0f,0b,07,03,3e,3a,36,32,2e,2a,26,22,1e,1a,16,12,0e,0a,06,02,3d,39,35,31,2d,29,25,21,1d,19,15,11,0d,09,05,01,3c,38,34,30,2c,28,24,20,1c,18,14,10,0c,08,04,00 + p1

	In order to make the result amenable to loop-based execution, we need to encode the indices
	to the left & right of the + in easily-computable-index fashion. This scheme has 2 ingredients:
	[Use i as loop idx in comments, although use 'k' in actual code]

	[0] RHS is simply i = 0,f,e,..,1, i.e. plo[(16-i) & (-(i > 0))]

	[1] Note that the (/= p10, where 10 = hex) row data are simply circular-shift perms of the basic (row 0) set;
	For row i (counting from 0 to 15), the leftward-shift count needing to be applied = 0,f,e,...,1, same as plo idx in [0].
		*/
		tptr = t;
		for(k = 0; k < 16; ++k)
		{
			k0 = (16-k) & (-(k > 0));	// Common index into plo[] and circ-shift count
			iptr = dit_p10_cperms+k0;
			for(l = 0; l < ODD_RADIX; l++)
			{
				ioff[l] = phi[iptr[l]];
			}
			RADIX_63_DIT(
				(double *)tptr, toff, 1,
				a+j1+plo[k0], ioff, RE_IM_STRIDE
			);
			tptr += ODD_RADIX;
		}
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy1008_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr,*addi;
		struct complex *tptr;
		const int pfetch_dist = PFETCH_DIST;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p1,p2,p3,p4;
		int poff[RADIX>>2];
	// DFT stuff:
		int kk, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf,po_kperm[16],*po_ptr = &(po_kperm[0]);
		int plo[16], phi[ODD_RADIX], toff[ODD_RADIX];
		int ioff[ODD_RADIX];	// ioff elts get recomputed on-the-fly for each radix-63 call
		uint64 i64;
	// DIF:
		// Need storage for circular-shifts perms of a basic 63-vector, with shift count in [0,60] that means 2*63-3 elts:
		const uint8 *iptr, dif_p10_cperms[128] = {	// Pad byte-array to next-higher 8-byte multiple
			0,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x07,0x06,0x05,0x04,0x03,0x02,0x01,
			0,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09,0x08,0x07,0x06,0x05,0x04
		},	dif_phi[64] = {	// Only 63 entries, but again pad out to nearest 8-multiple
			0,0x02,0x01,0x08,0x07,0x06,0x05,0x04,0x03,0x3e,0x3d,0x3c,0x3b,0x3a,0x39,0x38,0x37,0x36,0x35,0x34,0x33,0x32,0x31,0x30,0x2f,0x2e,0x2d,0x2c,0x2b,0x2a,0x29,0x28,0x27,0x26,0x25,0x24,0x23,0x22,0x21,0x20,0x1f,0x1e,0x1d,0x1c,0x1b,0x1a,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0x0f,0x0e,0x0d,0x0c,0x0b,0x0a,0x09
		};
		// Low parts [p0-f] of output-index perms - need one uint64 (16 x 4-bits) for each radix-16 DFT:
		const uint64 dif16_oidx_lo[ODD_RADIX] = {
			0x0123456789abcdefull,0x54763201dcfeba89ull,0xab98efdc67541032ull,0x1032547698badcfeull,0xcdefab9823106754ull,
			0x76450123fecd89abull,0x89abcdef45672310ull,0x32017645ba89fecdull,0xefdc98ba10325476ull,0x76450123fecd89abull,
			0x98badcfe54763201ull,0x23106754ab98efdcull,0xefdc98ba10325476ull,0x45672310cdefab98ull,0xba89fecd76450123ull,
			0x0123456789abcdefull,0xdcfeba8932017645ull,0x67541032efdc98baull,0xba89fecd76450123ull,0x1032547698badcfeull,
			0xcdefab9823106754ull,0x67541032efdc98baull,0x89abcdef45672310ull,0x32017645ba89fecdull,0xfecd89ab01234567ull,
			0x54763201dcfeba89ull,0xab98efdc67541032ull,0x32017645ba89fecdull,0xefdc98ba10325476ull,0x45672310cdefab98ull,
			0xab98efdc67541032ull,0x0123456789abcdefull,0xdcfeba8932017645ull,0x76450123fecd89abull,0x98badcfe54763201ull,
			0x23106754ab98efdcull,0xdcfeba8932017645ull,0x67541032efdc98baull,0x89abcdef45672310ull,0x23106754ab98efdcull,
			0xfecd89ab01234567ull,0x54763201dcfeba89ull,0xba89fecd76450123ull,0x1032547698badcfeull,0xcdefab9823106754ull,
			0x54763201dcfeba89ull,0xab98efdc67541032ull,0x0123456789abcdefull,0xcdefab9823106754ull,0x76450123fecd89abull,
			0x98badcfe54763201ull,0x32017645ba89fecdull,0xefdc98ba10325476ull,0x45672310cdefab98ull,0x98badcfe54763201ull,
			0x23106754ab98efdcull,0xfecd89ab01234567ull,0x45672310cdefab98ull,0xba89fecd76450123ull,0x1032547698badcfeull,
			0xdcfeba8932017645ull,0x67541032efdc98baull,0x89abcdef45672310ull
		};
	// DIT:
		// Need storage for circular-shifts perms of a basic 63-vector, with shift count in [0,15] that means 63+15 elts:
		const uint8 dit_p10_cperms[80] = {	// Pad byte-array to next-higher 8-byte multiple
			0,0x3b,0x37,0x33,0x2f,0x2b,0x27,0x23,0x1f,0x1b,0x17,0x13,0x0f,0x0b,0x07,0x03,0x3e,0x3a,0x36,0x32,0x2e,0x2a,0x26,0x22,0x1e,0x1a,0x16,0x12,0x0e,0x0a,0x06,0x02,0x3d,0x39,0x35,0x31,0x2d,0x29,0x25,0x21,0x1d,0x19,0x15,0x11,0x0d,0x09,0x05,0x01,0x3c,0x38,0x34,0x30,0x2c,0x28,0x24,0x20,0x1c,0x18,0x14,0x10,0x0c,0x08,0x04,
			0,0x3b,0x37,0x33,0x2f,0x2b,0x27,0x23,0x1f,0x1b,0x17,0x13,0x0f,0x0b,0x07
		},	dit_phi[64] = {	// Only 63 entries, but again pad out to nearest 8-multiple
			0,0x02,0x01,0x06,0x08,0x07,0x04,0x03,0x05,0x2f,0x2e,0x2d,0x35,0x34,0x33,0x30,0x32,0x31,0x1f,0x1e,0x20,0x1d,0x1c,0x1b,0x23,0x22,0x21,0x0f,0x11,0x10,0x0d,0x0c,0x0e,0x0b,0x0a,0x09,0x3e,0x3d,0x3c,0x39,0x3b,0x3a,0x37,0x36,0x38,0x26,0x25,0x24,0x2c,0x2b,0x2a,0x27,0x29,0x28,0x16,0x15,0x17,0x14,0x13,0x12,0x1a,0x19,0x18
		};
		// Low parts [p0-f] of input-index perms - need one uint64 (16 x 4-bits) for each radix-16 DFT:
		const uint64 dit16_iidx_lo[ODD_RADIX] = {
			0x01327654fedcba98ull,0x5467102398abefcdull,0xab89cdfe23014576ull,0x76543210ba98dcefull,0x10236745efcdab89ull,
			0xcdfe89ba45760132ull,0x32105467dcef98abull,0xefcdab8967452301ull,0x89bafedc01327654ull,0xfedcba9876543210ull,
			0x5467102398abefcdull,0xab89cdfe23014576ull,0xba98dcef32105467ull,0x10236745efcdab89ull,0xcdfe89ba45760132ull,
			0x32105467dcef98abull,0x67452301ab89cdfeull,0x89bafedc01327654ull,0xfedcba9876543210ull,0x5467102398abefcdull,
			0x23014576cdfe89baull,0xba98dcef32105467ull,0x10236745efcdab89ull,0xcdfe89ba45760132ull,0xdcef98ab54671023ull,
			0x67452301ab89cdfeull,0x89bafedc01327654ull,0xfedcba9876543210ull,0x98abefcd10236745ull,0x23014576cdfe89baull,
			0xba98dcef32105467ull,0x10236745efcdab89ull,0x4576013289bafedcull,0xdcef98ab54671023ull,0x67452301ab89cdfeull,
			0x89bafedc01327654ull,0x76543210ba98dcefull,0x98abefcd10236745ull,0x23014576cdfe89baull,0xba98dcef32105467ull,
			0xefcdab8967452301ull,0x4576013289bafedcull,0xdcef98ab54671023ull,0x67452301ab89cdfeull,0x01327654fedcba98ull,
			0x76543210ba98dcefull,0x98abefcd10236745ull,0x23014576cdfe89baull,0x32105467dcef98abull,0xefcdab8967452301ull,
			0x4576013289bafedcull,0xdcef98ab54671023ull,0xab89cdfe23014576ull,0x01327654fedcba98ull,0x76543210ba98dcefull,
			0x98abefcd10236745ull,0xcdfe89ba45760132ull,0x32105467dcef98abull,0xefcdab8967452301ull,0x4576013289bafedcull,
			0x5467102398abefcdull,0xab89cdfe23014576ull,0x01327654fedcba98ull
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
		double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7, *add8, *add9, *adda, *addb, *addc, *addd, *adde, *addf;
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
		vec_dbl *tmp,*tm0,*tm1,*tm2;	// utility ptrs
		int *itmp,*itm2;	// Pointer into the bjmodn array
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		vec_dbl *max_err, *sse2_rnd, *half_arr,
			*two,*one,*sqrt2,*isrt2,*cc0,*ss0,	// radix-16 DFT roots
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
		for(l = 0; l < 16; l++)
		{
			plo[l] = l*NDIVR;
			plo[l] += ( (plo[l] >> DAT_BITS) << PAD_BITS );
		}
		for(l = 0; l < ODD_RADIX; l++)
		{
			phi[l] = (l<<4)*NDIVR;
			phi[l] += ( (phi[l] >> DAT_BITS) << PAD_BITS );
			toff[l] = l+l;
		}

	#ifdef USE_SSE2
		uint32 nwt16 = nwt << L2_SZ_VD;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code

		r00 = tmp = thread_arg->r00;	// Head of RADIX*vec_cmplx-sized local store #1
		tmp += 0x7e0;	s1p00 = tmp;	// Head of RADIX*vec_cmplx-sized local store #2
		tmp += 0x7e0;
		two       = tmp + 0x00;	// AVX+ versions of various DFT macros assume consts 2.0,1.0,isrt2 laid out thusly
		one       = tmp + 0x01;
		sqrt2     = tmp + 0x02;
		isrt2     = tmp + 0x03;
		cc0       = tmp + 0x04;
		ss0       = tmp + 0x05;
		tmp += 0x6;	// sc_ptr += 0xfc6
	  #ifdef USE_AVX512
		cy_r = tmp;	cy_i = tmp+0x07e;	tmp += 2*0x07e;	// RADIX/8 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod in AVX mode
	  #elif defined(USE_AVX)
		cy_r = tmp;	cy_i = tmp+0x0fc;	tmp += 2*0x0fc;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod in AVX mode
	  #else
		cy_r = tmp;	cy_i = tmp+0x1f8;	tmp += 2*0x1f8;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
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

		sign_mask = (uint64*)(r00 + radix1008_creals_in_local_store);
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
		#include "radix1008_main_carry_loop.h"

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
