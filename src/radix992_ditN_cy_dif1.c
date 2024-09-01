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
#include "radix31.h"

#ifdef USE_SSE2	//*** NO SIMD support for this experimental radix
	#undef USE_SSE2
	#undef USE_AVX
	#undef USE_AVX2
	#undef USE_AVX512
#endif

#define RADIX 992	// Use #define rather than const int to ensure it's really a compile-time const in the C sense
#define ODD_RADIX 31	// ODD_RADIX = [radix >> trailz(radix)]

#define EPS 1e-10

#define USE_COMPACT_OBJ_CODE	1
/*
Apr 2014: Scalar-double tests of this radix indicate *really* bad ROE levels, worse than n = 960*2^k for same input!
Mersenne-mod: Test of default exponent @992K given by given_N_get_maxP barfs quickly:
....
M19388461 Roundoff warning on iteration       24, maxerr =   0.406250000000
M19388461 Roundoff warning on iteration       25, maxerr =   0.472656250000
 ERROR ERROR...Halting test of exponent 19388461

...so instead compare default exponent @960K using both that length and 992K:
100 iterations of M18776473 with FFT length 1015808 = 992 K
Res64: CA7D81B22AE24935. AvgMaxErr = 0.205201939. MaxErr = 0.250000000. <*** ROE almost same as 960k! ***
Res mod 2^36     =           9309407541
Res mod 2^35 - 1 =          24317941565
Res mod 2^36 - 1 =          67706175547

100 iterations of M18776473 with FFT length 983040 = 960 K
Res64: CA7D81B22AE24935. AvgMaxErr = 0.231762695. MaxErr = 0.312500000. Program: E3.0x
Res mod 2^36     =           9309407541
Res mod 2^35 - 1 =          24317941565
Res mod 2^36 - 1 =          67706175547

Fermat-mod:
100 iterations of F24 with FFT length 1015808 = 992 K
Res64: DB9F01963ED9DC8B. AvgMaxErr = 0.013013567. MaxErr = 0.017089844. <*** ROE almost same as 960k! ***
Res mod 2^36     =          26824268939
Res mod 2^35 - 1 =          27887793041
Res mod 2^36 - 1 =          13169874547

100 iterations of F24 with FFT length 983040 = 960 K
Res64: DB9F01963ED9DC8B. AvgMaxErr = 0.013779994. MaxErr = 0.017578125. Program: E3.0x
Res mod 2^36     =          26824268939
Res mod 2^35 - 1 =          27887793041
Res mod 2^36 - 1 =          13169874547

Compare time+ROE (1-thread, scalar-double) for n = 960 (R0 = 15*4), 992 (R0 = 31*32), 1008 (R0 == 63), 1024 (R0 = 64) K on F24:

N (K)	Radices used		ROE (avg, max)			T(sec for 100 iter)
 960	240,8,16,16		0.013624355. 0.015625000	12.6 (60,16,16,32 ==> T = 20.7, 960,16,32 ==> T=19.5, w/similar ROE)
 992	992,16,32		0.013013567. 0.017089844	19.6 <*** ROE worse than 960K makes this a nonstarter (But T not bad, comparable to 960,16,32)
1008	63,16,16,32		0.004876273. 0.005859375	19.5 (Was 24.5, but impl compact-obj-code cut .o from 220KB to 42KB and Core 2 timing by 20% as shown)
1024	256,8,16,16		0.002308873. 0.002441406	12.5 (64,16,16,32 ==> T = 19.0, 1024,16,32 ==> T=19.8, both w/similar ROE)
*/

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
	  #ifdef USE_AVX
		int kcycle[ODD_RADIX];
		int lcycle[ODD_RADIX];
	  #endif
	#endif

	// double data:
		double maxerr;
		double scale;
		double prp_mult;

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
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
		double cy_dat[2*RADIX+4] __attribute__ ((__aligned__(8)));	// Enforce min-alignment of 8 bytes in 32-bit builds.
	};

#endif

/***************/

int radix992_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-992 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-992 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix992_ditN_cy_dif1";
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
	int NDIVR,i,j,j1,jt,jstart,jhi,full_pass,khi,l,outer;
#ifndef MULTITHREAD
	int j2,jp,ntmp;
#endif
#ifdef USE_SSE2
	int nbytes;
	uint32 nwt16 = nwt << L2_SZ_VD;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code
#endif
	// Need these both in scalar mode and to ease the SSE2-array init...dimension = ODD_RADIX;
	// In order to ease the ptr-access for the || routine, lump these 4*ODD_RADIX doubles together with copies of
	// the 4 in the passed-in bs[2] and bsinv[2] arrays [and used in this 4-double form by the mersenne-mod carry macros]
	// into a single foo_array[4*(ODD_RADIX+1)], then convert what used to be disparate ODD_RADIX-sized arrays to pointers.
	static double foo_array[(ODD_RADIX+1)<<2], *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr;
#ifndef MULTITHREAD
	static double bs,bsinv;
#endif

	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
#ifdef MULTITHREAD
	double target_cy = 0;
#endif
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3, nsave = 0;
	static int poff[RADIX>>2];
#ifndef MULTITHREAD
	// Need storage for circular-shifts perms of a basic 31-vector, with shift count in [0,31] that means 2*31 elts:
	static int dif_p20_cperms[62], plo[32],phi[62],jj[32], *iptr;
	int idx,pidx,mask,lshift, is_even,is_odd, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf, o[32];	// o[] stores o-address offsets for current radix-32 DFT in the 31x-loop
	uint64 i64;
// DIF:
	// Low parts [p0-f] of output-index perms:
	const uint64 dif_perm16[16] = {
		0x0123456789abcdefull, 0xfecd89ab01234567ull, 0x76450123fecd89abull, 0xba89fecd76450123ull,
		0x32017645ba89fecdull, 0xdcfeba8932017645ull, 0x54763201dcfeba89ull, 0x98badcfe54763201ull,
		0x1032547698badcfeull, 0xefdc98ba10325476ull, 0x67541032efdc98baull, 0xab98efdc67541032ull,
		0x23106754ab98efdcull, 0xcdefab9823106754ull, 0x45672310cdefab98ull, 0x89abcdef45672310ull
	};
// DIT:
	// Low parts [p0-f] of output-index perms:
	const uint64 dit_perm16[16] = {
		0x01327654fedcba98ull,0xfedcba9876543210ull,0x76543210ba98dcefull,0xba98dcef32105467ull,
		0x32105467dcef98abull,0xdcef98ab54671023ull,0x5467102398abefcdull,0x98abefcd10236745ull,
		0x10236745efcdab89ull,0xefcdab8967452301ull,0x67452301ab89cdfeull,0xab89cdfe23014576ull,
		0x23014576cdfe89baull,0xcdfe89ba45760132ull,0x4576013289bafedcull,0x89bafedc01327654ull
	};
#endif
	static double radix_inv, n2inv;
	double scale,dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
  #ifndef MULTITHREAD
	double *addr, *addi;
	int *itmp;	// Pointer into the bjmodn array
  #endif
	struct complex t[RADIX]
  #ifndef MULTITHREAD
	  , *tptr
  #endif
	  ;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	static int n_div_nwt;
  #ifndef MULTITHREAD
	int col,co2,co3,m,m2;
  #endif
  #ifdef USE_AVX
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #elif !defined(MULTITHREAD)
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
  #endif
  #ifndef MULTITHREAD
	double rt,it, wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff
  #endif
	// indices into weights arrays (mod NWT):
	static int ii[ODD_RADIX] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	/* These are used in conjunction with the langth-ODD_RADIX arrays in the USE_SCALAR_CARRY code flow;
	In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
	static int wts_idx_incr = 0
		,icycle[ODD_RADIX];
  #ifndef MULTITHREAD
	static int ic;
  #endif
#ifdef USE_SSE2
	static int wts_idx_inc2 = 0;
	static int jcycle[ODD_RADIX],jc;
  #ifdef USE_AVX
	static int kcycle[ODD_RADIX];	// NB: kc already declared as part of k0-f set above
	static int lcycle[ODD_RADIX],lc;
  #endif

	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else
	double *add0,*add1,*add2,*add3;
  #endif	// MULTITHREAD

	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
	const double crnd = 3.0*0x4000000*0x2000000;
	static vec_dbl *max_err, *sse2_rnd, *half_arr, *isrt2,
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

#endif	// USE_SSE2

#ifdef MULTITHREAD

	static struct cy_thread_data_t *tdat = 0x0;
	// Threadpool-based dispatch stuff:
	//static int main_work_units = 0;
	static int pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)cy992_process_chunk, NULL, 0x0};

#elif 1//!defined(USE_SSE2)

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

	#ifdef MULTITHREAD

		/* #Chunks ||ized in carry step is ideally a power of 2, so use the smallest
		power of 2 that is >= the value of the global NTHREADS (but still <= MAX_THREADS):
		*/
		if(isPow2(NTHREADS))
			CY_THREADS = NTHREADS;
		else
		{
			i = leadz32(NTHREADS);
			CY_THREADS = (((uint32)NTHREADS << i) & 0x80000000) >> (i-1);
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
			tdat[ithread].si  = si;
			tdat[ithread].rn0 = rn0;
			tdat[ithread].rn1 = rn1;

		// This array pointer must be set based on vec_dbl-sized alignment at runtime for each thread:
			for(l = 0; l < 4; l++) {
				if( ((uintptr_t)&tdat[ithread].cy_dat[l] & SZ_VDM1) == 0 ) {
					tdat[ithread].cy_r = &tdat[ithread].cy_dat[l];
					tdat[ithread].cy_i = tdat[ithread].cy_r + RADIX;
				//	fprintf(stderr,"%d-byte-align cy_dat array at element[%d]\n",SZ_VD,l);
					break;
				}
			}
			ASSERT(l < 4, "Failed to align cy_dat array!");
		}
	#endif

	#ifdef USE_SSE2

		ASSERT(((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use vector-double type size (16 bytes for SSE2, 32 for AVX) to alloc a block of local storage
		// consisting of 128*2 vec_dbl and (8+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix992_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix992_creals_in_local_store);
		ASSERT(((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

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
		tmp += 0x08;	// += 0x8 => sc_ptr + 0xf44
	  #ifdef USE_AVX
		cy_r = tmp;	cy_i = tmp+0x0f0;	tmp += 0x1e0;	// RADIX/4 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0x1e0 + 2 => sc_ptr += 0x1126
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 68 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy_r = tmp;	cy_i = tmp+0x1e0;	tmp += 0x3c0;	// RADIX/2 vec_dbl slots for each of cy_r and cy_i carry sub-arrays
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// += 0x3c0 + 2 => sc_ptr += 0x1306
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*ODD_RADIX] x 16 for Fermat-mod */
	  #endif
		ASSERT(half_arr_offset992 == (uint32)(half_arr-sc_ptr), "half_arr_offset mismatches actual!");
		ASSERT((radix992_creals_in_local_store << L2_SZ_VD) >= ((long)half_arr - (long)r00) + (20 << L2_SZ_VD), "radix992_creals_in_local_store checksum failed!");

		/* These remain fixed: */
		VEC_DBL_INIT(isrt2,ISRT2);
		VEC_DBL_INIT(sse2_c3m1, c3m1);
		VEC_DBL_INIT(sse2_s   , s   );
		VEC_DBL_INIT(sse2_cn1 , cn1 );
		VEC_DBL_INIT(sse2_cn2 , cn2 );
		VEC_DBL_INIT(sse2_ss3 , ss3 );
		VEC_DBL_INIT(sse2_sn1 , sn1 );
		VEC_DBL_INIT(sse2_sn2 , sn2 );

		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
		VEC_DBL_INIT(sse2_rnd, crnd);

		// Propagate the above consts to the remaining threads:
		nbytes = (int)cy_r - (int)isrt2;	// #bytes in 1st of above block of consts
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
		/* Simple qfloat-based loop to crank out the roots: */
			struct qfloat qc,qs,qx,qy,qt,qn,qmul;
			qt = qfdiv(QPIHALF, i64_to_q((int64)RADIX));	// Pi/2/RADIX
			qc = qfcos(qt);	qs = qfsin(qt);
			qx = QONE;		qy = QZRO;
			for(j = 1; j <= RADIX; j++) {
				// Up-multiply the complex exponential:
				qn = qfmul(qx, qc); qt = qfmul(qy, qs); qmul = qfsub(qn, qt);	// Store qxnew in qmul for now.
				qn = qfmul(qx, qs); qt = qfmul(qy, qc); qy   = qfadd(qn, qt); qx = qmul;
				printf("j = %3u: cos = %#16" PRIX64 "\n",j,qfdbl_as_uint64(qx));
			}
			exit(0);
		#endif

		tmp = base_negacyclic_root + RADIX*2;	// First 2*RADIX slots reserved for RADIX/4 copies of the Re/Im parts of the 4 base multipliers
		tm2 = tmp + RADIX/2 - 1;
										tmp->d0 = 1.0;	(tmp+1)->d0 = 0.0;
		tmp64 = 0x3FEFFFD3151E5533ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(  1*I*Pi/480) = sin(239*I*Pi/480) */
		tmp64 = 0x3FEFFF4C54F76E1Cull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(  2*I*Pi/480) = sin(238*I*Pi/480) */
		tmp64 = 0x3FEFFE6BC105954Eull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(  3*I*Pi/480) = sin(237*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEFFD315BBF4275ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(  4*I*Pi/480) = sin(236*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEFFB9D2897136Eull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(  5*I*Pi/480) = sin(235*I*Pi/480) */
		tmp64 = 0x3FEFF9AF2BFBC297ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(  6*I*Pi/480) = sin(234*I*Pi/480) */
		tmp64 = 0x3FEFF7676B581A63ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(  7*I*Pi/480) = sin(233*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEFF4C5ED12E61Dull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(  8*I*Pi/480) = sin(232*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEFF1CAB88EDFF7ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(  9*I*Pi/480) = sin(231*I*Pi/480) */
		tmp64 = 0x3FEFEE75D62A9C46ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 10*I*Pi/480) = sin(230*I*Pi/480) */
		tmp64 = 0x3FEFEAC74F40720Cull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 11*I*Pi/480) = sin(229*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEFE6BF2E2660AFull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 12*I*Pi/480) = sin(228*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEFE25D7E2DF2F8ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 13*I*Pi/480) = sin(227*I*Pi/480) */
		tmp64 = 0x3FEFDDA24BA41F4Dull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 14*I*Pi/480) = sin(226*I*Pi/480) */
		tmp64 = 0x3FEFD88DA3D12526ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 15*I*Pi/480) = sin(225*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEFD31F94F867C6ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 16*I*Pi/480) = sin(224*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEFCD582E584632ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 17*I*Pi/480) = sin(223*I*Pi/480) */
		tmp64 = 0x3FEFC7378029F05Full;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 18*I*Pi/480) = sin(222*I*Pi/480) */
		tmp64 = 0x3FEFC0BD9BA139AEull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 19*I*Pi/480) = sin(221*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEFB9EA92EC689Bull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 20*I*Pi/480) = sin(220*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEFB2BE793403B9ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 21*I*Pi/480) = sin(219*I*Pi/480) */
		tmp64 = 0x3FEFAB39629A9BE3ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 22*I*Pi/480) = sin(218*I*Pi/480) */
		tmp64 = 0x3FEFA35B643C93B9ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 23*I*Pi/480) = sin(217*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEF9B24942FE45Cull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 24*I*Pi/480) = sin(216*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEF92950983DF6Bull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 25*I*Pi/480) = sin(215*I*Pi/480) */
		tmp64 = 0x3FEF89ACDC40EE4Bull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 26*I*Pi/480) = sin(214*I*Pi/480) */
		tmp64 = 0x3FEF806C25684EA8ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 27*I*Pi/480) = sin(213*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEF76D2FEF3CC4Bull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 28*I*Pi/480) = sin(212*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEF6CE183D57825ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 29*I*Pi/480) = sin(211*I*Pi/480) */
		tmp64 = 0x3FEF6297CFF75CB0ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 30*I*Pi/480) = sin(210*I*Pi/480) */
		tmp64 = 0x3FEF57F6003B2F91ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 31*I*Pi/480) = sin(209*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEF4CFC327A0080ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 32*I*Pi/480) = sin(208*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEF41AA8583E57Eull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 33*I*Pi/480) = sin(207*I*Pi/480) */
		tmp64 = 0x3FEF3601191FA459ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 34*I*Pi/480) = sin(206*I*Pi/480) */
		tmp64 = 0x3FEF2A000E0A5970ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 35*I*Pi/480) = sin(205*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEF1DA785F71BCEull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 36*I*Pi/480) = sin(204*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEF10F7A38E9E90ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 37*I*Pi/480) = sin(203*I*Pi/480) */
		tmp64 = 0x3FEF03F08A6ECF94ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 38*I*Pi/480) = sin(202*I*Pi/480) */
		tmp64 = 0x3FEEF6925F2A7380ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 39*I*Pi/480) = sin(201*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEEE8DD4748BF15ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 40*I*Pi/480) = sin(200*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEEDAD16944EDD0ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 41*I*Pi/480) = sin(199*I*Pi/480) */
		tmp64 = 0x3FEECC6EEC8DD5E9ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 42*I*Pi/480) = sin(198*I*Pi/480) */
		tmp64 = 0x3FEEBDB5F9857999ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 43*I*Pi/480) = sin(197*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEEAEA6B98095C0ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 44*I*Pi/480) = sin(196*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEE9F4156C62DDAull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 45*I*Pi/480) = sin(195*I*Pi/480) */
		tmp64 = 0x3FEE8F85FC8F1553ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 46*I*Pi/480) = sin(194*I*Pi/480) */
		tmp64 = 0x3FEE7F74D705762Bull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 47*I*Pi/480) = sin(193*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEE6F0E134454FFull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 48*I*Pi/480) = sin(192*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEE5E51DF571265ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 49*I*Pi/480) = sin(191*I*Pi/480) */
		tmp64 = 0x3FEE4D406A38E9ABull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 50*I*Pi/480) = sin(190*I*Pi/480) */
		tmp64 = 0x3FEE3BD9E3D46CEFull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 51*I*Pi/480) = sin(189*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEE2A1E7D02FE9Full;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 52*I*Pi/480) = sin(188*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEE180E678C4853ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 53*I*Pi/480) = sin(187*I*Pi/480) */
		tmp64 = 0x3FEE05A9D625AF0Full;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 54*I*Pi/480) = sin(186*I*Pi/480) */
		tmp64 = 0x3FEDF2F0FC71C4E5ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 55*I*Pi/480) = sin(185*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEDDFE40EFFB805ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 56*I*Pi/480) = sin(184*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEDCC83434ABF29ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 57*I*Pi/480) = sin(183*I*Pi/480) */
		tmp64 = 0x3FEDB8CECFB98376ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 58*I*Pi/480) = sin(182*I*Pi/480) */
		tmp64 = 0x3FEDA4C6EB9D87C2ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 59*I*Pi/480) = sin(181*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FED906BCF328D46ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 60*I*Pi/480) = sin(180*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FED7BBDB39DF5C3ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 61*I*Pi/480) = sin(179*I*Pi/480) */
		tmp64 = 0x3FED66BCD2EE2313ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 62*I*Pi/480) = sin(178*I*Pi/480) */
		tmp64 = 0x3FED51696819D42Bull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 63*I*Pi/480) = sin(177*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FED3BC3AEFF7F95ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 64*I*Pi/480) = sin(176*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FED25CBE464AB60ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 65*I*Pi/480) = sin(175*I*Pi/480) */
		tmp64 = 0x3FED0F8245F5427Full;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 66*I*Pi/480) = sin(174*I*Pi/480) */
		tmp64 = 0x3FECF8E71242E7ABull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 67*I*Pi/480) = sin(173*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FECE1FA88C445BBull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 68*I*Pi/480) = sin(172*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FECCABCE9D45D78ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 69*I*Pi/480) = sin(171*I*Pi/480) */
		tmp64 = 0x3FECB32E76B1D0F4ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 70*I*Pi/480) = sin(170*I*Pi/480) */
		tmp64 = 0x3FEC9B4F717E2C63ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 71*I*Pi/480) = sin(169*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEC83201D3D2C6Dull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 72*I*Pi/480) = sin(168*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEC6AA0BDD40210ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 73*I*Pi/480) = sin(167*I*Pi/480) */
		tmp64 = 0x3FEC51D198089406ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 74*I*Pi/480) = sin(166*I*Pi/480) */
		tmp64 = 0x3FEC38B2F180BDB1ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 75*I*Pi/480) = sin(165*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEC1F4510C18B95ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 76*I*Pi/480) = sin(164*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEC05883D2E7560ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 77*I*Pi/480) = sin(163*I*Pi/480) */
		tmp64 = 0x3FEBEB7CBF08957Dull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 78*I*Pi/480) = sin(162*I*Pi/480) */
		tmp64 = 0x3FEBD122DF6DDE43ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 79*I*Pi/480) = sin(161*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEBB67AE8584CAAull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 80*I*Pi/480) = sin(160*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEB9B85249D18A2ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 81*I*Pi/480) = sin(159*I*Pi/480) */
		tmp64 = 0x3FEB8041DFEBE2FBull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 82*I*Pi/480) = sin(158*I*Pi/480) */
		tmp64 = 0x3FEB64B166CDE0EEull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 83*I*Pi/480) = sin(157*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEB48D406A50540ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 84*I*Pi/480) = sin(156*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEB2CAA0DAB2702ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 85*I*Pi/480) = sin(155*I*Pi/480) */
		tmp64 = 0x3FEB1033CAF125F6ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 86*I*Pi/480) = sin(154*I*Pi/480) */
		tmp64 = 0x3FEAF3718E5E0C9Cull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 87*I*Pi/480) = sin(153*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEAD663A8AE2FDCull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 88*I*Pi/480) = sin(152*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEAB90A6B724C62ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 89*I*Pi/480) = sin(151*I*Pi/480) */
		tmp64 = 0x3FEA9B66290EA1A3ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 90*I*Pi/480) = sin(150*I*Pi/480) */
		tmp64 = 0x3FEA7D7734BA0A8Eull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 91*I*Pi/480) = sin(149*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FEA5F3DE27D13F2ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 92*I*Pi/480) = sin(148*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FEA40BA87311090ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 93*I*Pi/480) = sin(147*I*Pi/480) */
		tmp64 = 0x3FEA21ED787F2AEFull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 94*I*Pi/480) = sin(146*I*Pi/480) */
		tmp64 = 0x3FEA02D70CDF74DBull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 95*I*Pi/480) = sin(145*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE9E3779B97F4A8ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos( 96*I*Pi/480) = sin(144*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FE9C3CF7CBBB030ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos( 97*I*Pi/480) = sin(143*I*Pi/480) */
		tmp64 = 0x3FE9A3DF0929B594ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos( 98*I*Pi/480) = sin(142*I*Pi/480) */
		tmp64 = 0x3FE983A69A8C21B8ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos( 99*I*Pi/480) = sin(141*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE963268B572492ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(100*I*Pi/480) = sin(140*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FE9425F36C80335ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(101*I*Pi/480) = sin(139*I*Pi/480) */
		tmp64 = 0x3FE92150F8E417B1ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(102*I*Pi/480) = sin(138*I*Pi/480) */
		tmp64 = 0x3FE8FFFC2E77CEBAull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(103*I*Pi/480) = sin(137*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE8DE613515A328ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(104*I*Pi/480) = sin(136*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FE8BC806B151741ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(105*I*Pi/480) = sin(135*I*Pi/480) */
		tmp64 = 0x3FE89A5A2F91ABE4ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(106*I*Pi/480) = sin(134*I*Pi/480) */
		tmp64 = 0x3FE877EEE269D586ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(107*I*Pi/480) = sin(133*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE8553EE43DEF13ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(108*I*Pi/480) = sin(132*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FE8324A966F2AA5ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(109*I*Pi/480) = sin(131*I*Pi/480) */
		tmp64 = 0x3FE80F125B1E8028ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(110*I*Pi/480) = sin(130*I*Pi/480) */
		tmp64 = 0x3FE7EB96952B99DCull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(111*I*Pi/480) = sin(129*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE7C7D7A833BEC2ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(112*I*Pi/480) = sin(128*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FE7A3D5F890BAF9ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(113*I*Pi/480) = sin(127*I*Pi/480) */
		tmp64 = 0x3FE77F91EB57C602ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(114*I*Pi/480) = sin(126*I*Pi/480) */
		tmp64 = 0x3FE75B0BE65866FBull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(115*I*Pi/480) = sin(125*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE73644501B56CDull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(116*I*Pi/480) = sin(124*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FE7113B8FE16056ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(117*I*Pi/480) = sin(123*I*Pi/480) */
		tmp64 = 0x3FE6EBF20DA23E86ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(118*I*Pi/480) = sin(122*I*Pi/480) */
		tmp64 = 0x3FE6C668320B7884ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(119*I*Pi/480) = sin(121*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE6A09E667F3BCDull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(120*I*Pi/480) = sin(120*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FE67A951513345Cull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(121*I*Pi/480) = sin(119*I*Pi/480) */
		tmp64 = 0x3FE6544CA88F62DBull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(122*I*Pi/480) = sin(118*I*Pi/480) */
		tmp64 = 0x3FE62DC58C6CF0DBull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(123*I*Pi/480) = sin(117*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE607002CD5031Dull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(124*I*Pi/480) = sin(116*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FE5DFFCF69F89EDull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(125*I*Pi/480) = sin(115*I*Pi/480) */
		tmp64 = 0x3FE5B8BC57520F97ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(126*I*Pi/480) = sin(114*I*Pi/480) */
		tmp64 = 0x3FE5913EBD1E84E7ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(127*I*Pi/480) = sin(113*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE5698496E20BD8ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(128*I*Pi/480) = sin(112*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FE5418E5423C050ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(129*I*Pi/480) = sin(111*I*Pi/480) */
		tmp64 = 0x3FE5195C65137F0Cull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(130*I*Pi/480) = sin(110*I*Pi/480) */
		tmp64 = 0x3FE4F0EF3A88AAADull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(131*I*Pi/480) = sin(109*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE4C8474600EEEEull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(132*I*Pi/480) = sin(108*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FE49F64F99F0207ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(133*I*Pi/480) = sin(107*I*Pi/480) */
		tmp64 = 0x3FE47648C8296447ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(134*I*Pi/480) = sin(106*I*Pi/480) */
		tmp64 = 0x3FE44CF325091DD6ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(135*I*Pi/480) = sin(105*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE4236484487ABEull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(136*I*Pi/480) = sin(104*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FE3F99D5A91C51Full;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(137*I*Pi/480) = sin(103*I*Pi/480) */
		tmp64 = 0x3FE3CF9E1D2DFDB2ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(138*I*Pi/480) = sin(102*I*Pi/480) */
		tmp64 = 0x3FE3A56742039280ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(139*I*Pi/480) = sin(101*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE37AF93F9513EAull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(140*I*Pi/480) = sin(100*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FE350548CFFE7F2ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(141*I*Pi/480) = sin( 99*I*Pi/480) */
		tmp64 = 0x3FE32579A1FAFBDAull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(142*I*Pi/480) = sin( 98*I*Pi/480) */
		tmp64 = 0x3FE2FA68F6D5740Cull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(143*I*Pi/480) = sin( 97*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE2CF2304755A5Eull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(144*I*Pi/480) = sin( 96*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FE2A3A844564AA5ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(145*I*Pi/480) = sin( 95*I*Pi/480) */
		tmp64 = 0x3FE277F930881DAFull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(146*I*Pi/480) = sin( 94*I*Pi/480) */
		tmp64 = 0x3FE24C1643AD9295ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(147*I*Pi/480) = sin( 93*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE21FFFF8FAF674ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(148*I*Pi/480) = sin( 92*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FE1F3B6CC34CA8Bull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(149*I*Pi/480) = sin( 91*I*Pi/480) */
		tmp64 = 0x3FE1C73B39AE68C8ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(150*I*Pi/480) = sin( 90*I*Pi/480) */
		tmp64 = 0x3FE19A8DBE48A6C1ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(151*I*Pi/480) = sin( 89*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE16DAED770771Dull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(152*I*Pi/480) = sin( 88*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FE1409F031D897Eull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(153*I*Pi/480) = sin( 87*I*Pi/480) */
		tmp64 = 0x3FE1135EBFD0E8D7ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(154*I*Pi/480) = sin( 86*I*Pi/480) */
		tmp64 = 0x3FE0E5EE8C939850ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(155*I*Pi/480) = sin( 85*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE0B84EE8F52E9Dull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(156*I*Pi/480) = sin( 84*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FE08A80550A6FE5ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(157*I*Pi/480) = sin( 83*I*Pi/480) */
		tmp64 = 0x3FE05C83516BE635ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(158*I*Pi/480) = sin( 82*I*Pi/480) */
		tmp64 = 0x3FE02E585F347876ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(159*I*Pi/480) = sin( 81*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FE0000000000000ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(160*I*Pi/480) = sin( 80*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FDFA2F56BD3B979ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(161*I*Pi/480) = sin( 79*I*Pi/480) */
		tmp64 = 0x3FDF459207170FCEull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(162*I*Pi/480) = sin( 78*I*Pi/480) */
		tmp64 = 0x3FDEE7D6D7F64AD2ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(163*I*Pi/480) = sin( 77*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FDE89C4E59427B1ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(164*I*Pi/480) = sin( 76*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FDE2B5D3806F63Bull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(165*I*Pi/480) = sin( 75*I*Pi/480) */
		tmp64 = 0x3FDDCCA0D855B380ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(166*I*Pi/480) = sin( 74*I*Pi/480) */
		tmp64 = 0x3FDD6D90D07521CBull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(167*I*Pi/480) = sin( 73*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FDD0E2E2B44DE01ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(168*I*Pi/480) = sin( 72*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FDCAE79F48C726Cull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(169*I*Pi/480) = sin( 71*I*Pi/480) */
		tmp64 = 0x3FDC4E7538F866FCull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(170*I*Pi/480) = sin( 70*I*Pi/480) */
		tmp64 = 0x3FDBEE2106174F02ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(171*I*Pi/480) = sin( 69*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FDB8D7E6A56D476ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(172*I*Pi/480) = sin( 68*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FDB2C8E7500C0C6ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(173*I*Pi/480) = sin( 67*I*Pi/480) */
		tmp64 = 0x3FDACB523638033Bull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(174*I*Pi/480) = sin( 66*I*Pi/480) */
		tmp64 = 0x3FDA69CABEF5B501ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(175*I*Pi/480) = sin( 65*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FDA07F921061AD1ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(176*I*Pi/480) = sin( 64*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FD9A5DE6F05A44Bull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(177*I*Pi/480) = sin( 63*I*Pi/480) */
		tmp64 = 0x3FD9437BBC5DE90Aull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(178*I*Pi/480) = sin( 62*I*Pi/480) */
		tmp64 = 0x3FD8E0D21D42A377ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(179*I*Pi/480) = sin( 61*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FD87DE2A6AEA963ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(180*I*Pi/480) = sin( 60*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FD81AAE6E60E271ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(181*I*Pi/480) = sin( 59*I*Pi/480) */
		tmp64 = 0x3FD7B7368AD93C61ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(182*I*Pi/480) = sin( 58*I*Pi/480) */
		tmp64 = 0x3FD7537C13559D33ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(183*I*Pi/480) = sin( 57*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FD6EF801FCED33Cull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(184*I*Pi/480) = sin( 56*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FD68B43C8F5832Aull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(185*I*Pi/480) = sin( 55*I*Pi/480) */
		tmp64 = 0x3FD626C8282F1408ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(186*I*Pi/480) = sin( 54*I*Pi/480) */
		tmp64 = 0x3FD5C20E57929942ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(187*I*Pi/480) = sin( 53*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FD55D1771E5BAB9ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(188*I*Pi/480) = sin( 52*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FD4F7E492999AEEull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(189*I*Pi/480) = sin( 51*I*Pi/480) */
		tmp64 = 0x3FD49276D5C7BB48ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(190*I*Pi/480) = sin( 50*I*Pi/480) */
		tmp64 = 0x3FD42CCF582EDE82ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(191*I*Pi/480) = sin( 49*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FD3C6EF372FE950ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(192*I*Pi/480) = sin( 48*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FD360D790CAC12Eull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(193*I*Pi/480) = sin( 47*I*Pi/480) */
		tmp64 = 0x3FD2FA89839B2985ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(194*I*Pi/480) = sin( 46*I*Pi/480) */
		tmp64 = 0x3FD294062ED59F06ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(195*I*Pi/480) = sin( 45*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FD22D4EB2443163ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(196*I*Pi/480) = sin( 44*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FD1C6642E435B69ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(197*I*Pi/480) = sin( 43*I*Pi/480) */
		tmp64 = 0x3FD15F47C3BED971ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(198*I*Pi/480) = sin( 42*I*Pi/480) */
		tmp64 = 0x3FD0F7FA942E7E48ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(199*I*Pi/480) = sin( 41*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FD0907DC1930690ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(200*I*Pi/480) = sin( 40*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FD028D26E72EA99ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(201*I*Pi/480) = sin( 39*I*Pi/480) */
		tmp64 = 0x3FCF81F37BAE5D8Cull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(202*I*Pi/480) = sin( 38*I*Pi/480) */
		tmp64 = 0x3FCEB1E9A690650Eull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(203*I*Pi/480) = sin( 37*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FCDE189A594FBCCull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(204*I*Pi/480) = sin( 36*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FCD10D5C1B71B7Full;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(205*I*Pi/480) = sin( 35*I*Pi/480) */
		tmp64 = 0x3FCC3FD044DD3D45ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(206*I*Pi/480) = sin( 34*I*Pi/480) */
		tmp64 = 0x3FCB6E7B79D2ECC9ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(207*I*Pi/480) = sin( 33*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FCA9CD9AC4258F6ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(208*I*Pi/480) = sin( 32*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FC9CAED28ADE228ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(209*I*Pi/480) = sin( 31*I*Pi/480) */
		tmp64 = 0x3FC8F8B83C69A60Bull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(210*I*Pi/480) = sin( 30*I*Pi/480) */
		tmp64 = 0x3FC8263D35950926ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(211*I*Pi/480) = sin( 29*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FC7537E63143E2Eull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(212*I*Pi/480) = sin( 28*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FC6807E1489CB33ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(213*I*Pi/480) = sin( 27*I*Pi/480) */
		tmp64 = 0x3FC5AD3E9A500CADull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(214*I*Pi/480) = sin( 26*I*Pi/480) */
		tmp64 = 0x3FC4D9C24572B693ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(215*I*Pi/480) = sin( 25*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FC4060B67A85375ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(216*I*Pi/480) = sin( 24*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FC3321C534BC1BBull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(217*I*Pi/480) = sin( 23*I*Pi/480) */
		tmp64 = 0x3FC25DF75B55AF15ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(218*I*Pi/480) = sin( 22*I*Pi/480) */
		tmp64 = 0x3FC1899ED3561233ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(219*I*Pi/480) = sin( 21*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FC0B5150F6DA2D1ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(220*I*Pi/480) = sin( 20*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FBFC0B8C88EA05Aull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(221*I*Pi/480) = sin( 19*I*Pi/480) */
		tmp64 = 0x3FBE16EE4E236BF8ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(222*I*Pi/480) = sin( 18*I*Pi/480) */
		tmp64 = 0x3FBC6CCF5AF11FD1ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(223*I*Pi/480) = sin( 17*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FBAC2609B3C576Cull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(224*I*Pi/480) = sin( 16*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FB917A6BC29B42Cull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(225*I*Pi/480) = sin( 15*I*Pi/480) */
		tmp64 = 0x3FB76CA66BB0BC83ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(226*I*Pi/480) = sin( 14*I*Pi/480) */
		tmp64 = 0x3FB5C164588EB8DAull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(227*I*Pi/480) = sin( 13*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FB415E532398E49ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(228*I*Pi/480) = sin( 12*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FB26A2DA8D2974Aull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(229*I*Pi/480) = sin( 11*I*Pi/480) */
		tmp64 = 0x3FB0BE426D197A8Bull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(230*I*Pi/480) = sin( 10*I*Pi/480) */
		tmp64 = 0x3FAE245060BE0012ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(231*I*Pi/480) = sin(  9*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3FAACBC748EFC90Eull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(232*I*Pi/480) = sin(  8*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3FA772F2F75F573Cull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(233*I*Pi/480) = sin(  7*I*Pi/480) */
		tmp64 = 0x3FA419DCD176E0F7ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(234*I*Pi/480) = sin(  6*I*Pi/480) */
		tmp64 = 0x3FA0C08E3D596AEEull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(235*I*Pi/480) = sin(  5*I*Pi/480) */	tmp += 2;
		tmp64 = 0x3F9ACE214390CA91ull;	tmp->d0 = tm2->d0 = u64_to_f64(tmp64);	/* cos(236*I*Pi/480) = sin(  4*I*Pi/480) */	tm2 -= 2;
		tmp64 = 0x3F941ADACC128E22ull;	tmp->d1 = tm2->d3 = u64_to_f64(tmp64);	/* cos(237*I*Pi/480) = sin(  3*I*Pi/480) */
		tmp64 = 0x3F8ACEB7C72CA0A8ull;	tmp->d2 = tm2->d2 = u64_to_f64(tmp64);	/* cos(238*I*Pi/480) = sin(  2*I*Pi/480) */
		tmp64 = 0x3F7ACEDD6862D0D7ull;	tmp->d3 = tm2->d1 = u64_to_f64(tmp64);	/* cos(239*I*Pi/480) = sin(  1*I*Pi/480) */	tmp += 2;

		tmp = base_negacyclic_root + RADIX*2;	// reset to point to start of above block
		nbytes = RADIX*SZ_VD/2;	// RADIX/4 AVX-register-sized complex data

	  #else	// HIACC = false:

		// lower-precision version, which yields slightly more roundoff error, but is simpler and more storage-compact.
		// Init exp(j*I*Pi/2/RADIX), for j = 0-3:
		tmp = base_negacyclic_root + 8;	// First 8 slots reserved for Re/Im parts of the 4 base multipliers
		tmp->d0 = 1.0;
		tmp64 = 0x3FEFFFD3151E5533ull;	tmp->d1 = u64_to_f64(tmp64);	// cos(01*I*Pi/480)
		tmp64 = 0x3FEFFF4C54F76E1Cull;	tmp->d2 = u64_to_f64(tmp64);	// cos(02*I*Pi/480)
		tmp64 = 0x3FEFFE6BC105954Eull;	tmp->d3 = u64_to_f64(tmp64);	// cos(03*I*Pi/480)

		(++tmp)->d0 = 0.0;
		tmp64 = 0x3F941ADACC128E22ull;	tmp->d3 = u64_to_f64(tmp64);	// sin(03*I*Pi/480)
		tmp64 = 0x3F8ACEB7C72CA0A8ull;	tmp->d2 = u64_to_f64(tmp64);	// sin(02*I*Pi/480)
		tmp64 = 0x3F7ACEDD6862D0D7ull;	tmp->d1 = u64_to_f64(tmp64);	// sin(01*I*Pi/480)
		++tmp;
		tmp64 = 0x3FEFFD315BBF4275ull;	VEC_DBL_INIT(tmp, u64_to_f64(tmp64));	// cos(04*I*Pi/480)
		++tmp;
		tmp64 = 0x3F9ACE214390CA91ull;	VEC_DBL_INIT(tmp, u64_to_f64(tmp64));	// sin(04*I*Pi/480)
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

		nbytes = 64 << L2_SZ_VD;

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

		nbytes = 16 << L2_SZ_VD;
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

		for(l = 0; l < (RADIX>>2); l++) {
			poff[l] = (l<<2)*NDIVR;	// Corr. to every 4th plo[] term
			poff[l] += ( (poff[l] >> DAT_BITS) << PAD_BITS );
		}

	#ifndef MULTITHREAD
		// Init plo,phi arrays and circular-shifts perms of a basic 31-vector, with shift count in [0,31] that means 2 copies of a basic 31-elt set:
		for(l = 0; l < 31; l++)
		{
			plo[l] = l*NDIVR; phi[l] = plo[l]<<5;
			plo[l] += ( (plo[l] >> DAT_BITS) << PAD_BITS );
			phi[l] += ( (phi[l] >> DAT_BITS) << PAD_BITS );
			phi[31+l] = phi[l];	// Needed for the DIT cperm, which needs just a doubled version of phi[0-30]
			dif_p20_cperms[31-l] = dif_p20_cperms[62-l] = phi[l];	// p3c0,3a0,...,020; init in reverse order since init of needed phi elts run forward in same loop
		}	plo[l] = l*NDIVR;	plo[l] += ( (plo[l] >> DAT_BITS) << PAD_BITS );	// Need 1 more plo term than phi
		dif_p20_cperms[0] = dif_p20_cperms[31] = 0;	// Handle the leading 0 and its copy separately
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
				jcycle[i] = icycle[i] + wts_idx_incr;	jcycle[i] += ( (-(jcycle[i] < 0)) & nwt);
				tmp->d1 = wt_arr[jcycle[i]];
		  #ifdef USE_AVX
				kcycle[i] = jcycle[i] + wts_idx_incr;	kcycle[i] += ( (-(kcycle[i] < 0)) & nwt);
				tmp->d2 = wt_arr[kcycle[i]];
				lcycle[i] = kcycle[i] + wts_idx_incr;	lcycle[i] += ( (-(lcycle[i] < 0)) & nwt);
				tmp->d3 = wt_arr[lcycle[i]];
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
				icycle[i] <<= L2_SZ_VD;		jcycle[i] <<= L2_SZ_VD;
			#ifdef USE_AVX
				kcycle[i] <<= L2_SZ_VD;		lcycle[i] <<= L2_SZ_VD;
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
			tdat[ithread].r00      = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (vec_dbl *)((long)tdat[ithread].r00 + ((long)half_arr - (long)r00));
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
			  #ifdef USE_AVX
				tdat[0].kcycle[i] = kcycle[i];
				tdat[0].lcycle[i] = lcycle[i];
			  #endif
			//	printf("Thread 0 idx pair[%2d] = [%2d,%2d]\n",i,icycle[i],jcycle[i]);
			#endif
			}
			// For remaining threads, simulate the loop-evolution of the above indices.
			// Note that the non-thread-associated *cycle[] arry data will get changed fom their above-inited
			// values in the loop here, but that's OK because in || mode only the thread-associated values matter:
			for(ithread = 1; ithread < CY_THREADS; ++ithread) {
				jstart = 0;
				jhi = NDIVR/CY_THREADS;	// The earlier setting = NDIVR/CY_THREADS/2 was for simulating bjmodn evolution, must double that here
				// khi = 1 for Fermat-mod, thus no outer loop needed here
				for(j = jstart; j < jhi; j += stride)
				{
					for(i = 0; i < ODD_RADIX; i++) {
					#if 1//ndef USE_SSE2	// Scalar-double mode uses non-pointerized icycle values:
						icycle[i] += wts_idx_incr;		icycle[i] += ( (-(int)((uint32)icycle[i] >> 31)) & nwt);
					#else
						icycle[i] += wts_idx_inc2;		icycle[i] += ( (-(int)((uint32)icycle[i] >> 31)) & nwt16);
						jcycle[i] += wts_idx_inc2;		jcycle[i] += ( (-(int)((uint32)jcycle[i] >> 31)) & nwt16);
					  #ifdef USE_AVX
						kcycle[i] += wts_idx_inc2;		kcycle[i] += ( (-(int)((uint32)kcycle[i] >> 31)) & nwt16);
						lcycle[i] += wts_idx_inc2;		lcycle[i] += ( (-(int)((uint32)lcycle[i] >> 31)) & nwt16);
					  #endif
					#endif
					}
				}
				for(i = 0; i < ODD_RADIX; i++) {
					tdat[ithread].icycle[i] = icycle[i];
				#ifdef USE_SSE2
					tdat[ithread].wts_idx_inc2 = wts_idx_inc2;
					tdat[ithread].jcycle[i] = jcycle[i];
				  #ifdef USE_AVX
					tdat[ithread].kcycle[i] = kcycle[i];
					tdat[ithread].lcycle[i] = lcycle[i];
				  #endif
			//	printf("Thread %d idx pair[%2d] = [%2d,%2d]\n",ithread,i,icycle[i],jcycle[i]);
				#endif
				}
			}
			// Restore the original loop-start values of the cycle arrays, since we use these for init of inv-wts below:
			for(i = 0; i < ODD_RADIX; i++) {
				icycle[i] = tdat[0].icycle[i];
			#ifdef USE_SSE2
				jcycle[i] = tdat[0].jcycle[i];
			  #ifdef USE_AVX
				kcycle[i] = tdat[0].kcycle[i];
				lcycle[i] = tdat[0].lcycle[i];
			  #endif
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
		#ifdef MULTITHREAD
			// Compute wt = 2^(target_idx*sw % n)/n and its reciprocal:
			uint32 sw_idx_modn = ((uint64)target_idx*sw) % n;	// N is 32-bit, so only use 64-bit to hold intermediate product
			double target_wtfwd = pow(2.0, sw_idx_modn*0.5*n2inv);	// 0.5*n2inv = 0.5/(n/2) = 1.0/n
		#endif
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
		#ifdef MULTITHREAD
			target_cy  = target_wtfwd * (-(int)(2u << (itmp64 & 255)));
		#endif
		} else {
			target_idx = target_set = 0;
		#ifdef MULTITHREAD
			target_cy = -2.0;
		#endif
		}
	}

/*...The radix-992 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(i = 0; i < RADIX; i++) {
			_cy_r[i][ithread] = 0;
			_cy_i[i][ithread] = 0;
		}
	}
  #if 0	//ndef USE_SSE2	*** v20: Non-SIMD builds now also support shifted-residue
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
	// in order to so;ve the save/restore issue. We start from the (static, unmodified during loop) ii[]-index values:
	#ifndef MULTITHREAD
		for(i = 0; i < ODD_RADIX; i++) {
			/* Reinit *cycle indices their proper starting values - recall in SIMD mode these all are ( << 4): */
			icycle[i] = i;
		#ifdef USE_SSE2
			jcycle[i] = icycle[i] + wts_idx_incr;	jcycle[i] += ( (-(jcycle[i] < 0)) & nwt);
		  #ifdef USE_AVX
			kcycle[i] = jcycle[i] + wts_idx_incr;	kcycle[i] += ( (-(kcycle[i] < 0)) & nwt);
			lcycle[i] = kcycle[i] + wts_idx_incr;	lcycle[i] += ( (-(lcycle[i] < 0)) & nwt);
			kcycle[i] <<= L2_SZ_VD;		lcycle[i] <<= L2_SZ_VD;
		  #endif
			icycle[i] <<= L2_SZ_VD;		jcycle[i] <<= L2_SZ_VD;
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
			for(l = 0; l < RADIX; l++) {
				bjmodn[l] = _bjmodn[l][ithread];
			}
			/* init carries	*/
		#if 0//def USE_AVX	// AVX and AVX2 both use 256-bit registers
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
		#include "radix992_main_carry_loop.h"

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
		ASSERT(0x0 == cy992_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
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

/***************/

void radix992_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-992 = 31x32 complex DIF FFT pass on the data in the length-N real vector A.
*/
	int l, j,j1/* ,j2 */, jj[32], *iptr;
	static int NDIVR,first_entry=TRUE;
	struct complex t[RADIX], *tptr;
#if USE_COMPACT_OBJ_CODE
	int jp;
	// Need storage for circular-shifts perms of a basic 31-vector, with shift count in [0,31] that means 2*31 elts:
	static int dif_p20_cperms[62], plo[32],phi[31];
	int idx,pidx, is_even,is_odd, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf, o[32];	// o[] stores o-address offsets for current radix-32 DFT in the 31x-loop
	uint64 i64;
	// Low parts [p0-f] of output-index perms - need 64x4-bits for each radix-64 DFT:
	const uint64 dif_perm16[16] = {
		0x0123456789abcdefull, 0xfecd89ab01234567ull, 0x76450123fecd89abull, 0xba89fecd76450123ull,
		0x32017645ba89fecdull, 0xdcfeba8932017645ull, 0x54763201dcfeba89ull, 0x98badcfe54763201ull,
		0x1032547698badcfeull, 0xefdc98ba10325476ull, 0x67541032efdc98baull, 0xab98efdc67541032ull,
		0x23106754ab98efdcull, 0xcdefab9823106754ull, 0x45672310cdefab98ull, 0x89abcdef45672310ull
	};
#else
	static int p[RADIX],q[RADIX], o[RADIX] = {
	   0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
	 975,974,972,973,968,969,970,971,960,961,962,963,964,965,966,967,991,990,988,989,984,985,986,987,976,977,978,979,980,981,982,983,
	 951,950,948,949,944,945,946,947,959,958,956,957,952,953,954,955,943,942,940,941,936,937,938,939,928,929,930,931,932,933,934,935,
	 903,902,900,901,896,897,898,899,911,910,908,909,904,905,906,907,919,918,916,917,912,913,914,915,927,926,924,925,920,921,922,923,
	 891,890,888,889,895,894,892,893,887,886,884,885,880,881,882,883,871,870,868,869,864,865,866,867,879,878,876,877,872,873,874,875,
	 843,842,840,841,847,846,844,845,839,838,836,837,832,833,834,835,859,858,856,857,863,862,860,861,855,854,852,853,848,849,850,851,
	 819,818,816,817,823,822,820,821,827,826,824,825,831,830,828,829,811,810,808,809,815,814,812,813,807,806,804,805,800,801,802,803,
	 771,770,768,769,775,774,772,773,779,778,776,777,783,782,780,781,787,786,784,785,791,790,788,789,795,794,792,793,799,798,796,797,
	 765,764,767,766,763,762,760,761,755,754,752,753,759,758,756,757,739,738,736,737,743,742,740,741,747,746,744,745,751,750,748,749,
	 717,716,719,718,715,714,712,713,707,706,704,705,711,710,708,709,733,732,735,734,731,730,728,729,723,722,720,721,727,726,724,725,
	 693,692,695,694,691,690,688,689,701,700,703,702,699,698,696,697,685,684,687,686,683,682,680,681,675,674,672,673,679,678,676,677,
	 645,644,647,646,643,642,640,641,653,652,655,654,651,650,648,649,661,660,663,662,659,658,656,657,669,668,671,670,667,666,664,665,
	 633,632,635,634,637,636,639,638,629,628,631,630,627,626,624,625,613,612,615,614,611,610,608,609,621,620,623,622,619,618,616,617,
	 585,584,587,586,589,588,591,590,581,580,583,582,579,578,576,577,601,600,603,602,605,604,607,606,597,596,599,598,595,594,592,593,
	 561,560,563,562,565,564,567,566,569,568,571,570,573,572,575,574,553,552,555,554,557,556,559,558,549,548,551,550,547,546,544,545,
	 513,512,515,514,517,516,519,518,521,520,523,522,525,524,527,526,529,528,531,530,533,532,535,534,537,536,539,538,541,540,543,542,
	 510,511,509,508,505,504,507,506,497,496,499,498,501,500,503,502,481,480,483,482,485,484,487,486,489,488,491,490,493,492,495,494,
	 462,463,461,460,457,456,459,458,449,448,451,450,453,452,455,454,478,479,477,476,473,472,475,474,465,464,467,466,469,468,471,470,
	 438,439,437,436,433,432,435,434,446,447,445,444,441,440,443,442,430,431,429,428,425,424,427,426,417,416,419,418,421,420,423,422,
	 390,391,389,388,385,384,387,386,398,399,397,396,393,392,395,394,406,407,405,404,401,400,403,402,414,415,413,412,409,408,411,410,
	 378,379,377,376,382,383,381,380,374,375,373,372,369,368,371,370,358,359,357,356,353,352,355,354,366,367,365,364,361,360,363,362,
	 330,331,329,328,334,335,333,332,326,327,325,324,321,320,323,322,346,347,345,344,350,351,349,348,342,343,341,340,337,336,339,338,
	 306,307,305,304,310,311,309,308,314,315,313,312,318,319,317,316,298,299,297,296,302,303,301,300,294,295,293,292,289,288,291,290,
	 258,259,257,256,262,263,261,260,266,267,265,264,270,271,269,268,274,275,273,272,278,279,277,276,282,283,281,280,286,287,285,284,
	 252,253,254,255,250,251,249,248,242,243,241,240,246,247,245,244,226,227,225,224,230,231,229,228,234,235,233,232,238,239,237,236,
	 204,205,206,207,202,203,201,200,194,195,193,192,198,199,197,196,220,221,222,223,218,219,217,216,210,211,209,208,214,215,213,212,
	 180,181,182,183,178,179,177,176,188,189,190,191,186,187,185,184,172,173,174,175,170,171,169,168,162,163,161,160,166,167,165,164,
	 132,133,134,135,130,131,129,128,140,141,142,143,138,139,137,136,148,149,150,151,146,147,145,144,156,157,158,159,154,155,153,152,
	 120,121,122,123,124,125,126,127,116,117,118,119,114,115,113,112,100,101,102,103, 98, 99, 97, 96,108,109,110,111,106,107,105,104,
	  72, 73, 74, 75, 76, 77, 78, 79, 68, 69, 70, 71, 66, 67, 65, 64, 88, 89, 90, 91, 92, 93, 94, 95, 84, 85, 86, 87, 82, 83, 81, 80,
	  48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 40, 41, 42, 43, 44, 45, 46, 47, 36, 37, 38, 39, 34, 35, 33, 32};
#endif

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
	#if USE_COMPACT_OBJ_CODE
		// Init plo,phi arrays and circular-shifts perms of a basic 31-vector, with shift count in [0,31] that means 2 copies of a basic 31-elt set:
		for(l = 0; l < 31; l++)
		{
			plo[l] = l*NDIVR; phi[l] = plo[l]<<5;
			plo[l] += ( (plo[l] >> DAT_BITS) << PAD_BITS );
			phi[l] += ( (phi[l] >> DAT_BITS) << PAD_BITS );
			dif_p20_cperms[31-l] = dif_p20_cperms[62-l] = phi[l];	// p3c0,3a0,...,020; init in reverse order since init of needed phi elts run forward in same loop
		}	plo[l] = l*NDIVR;	plo[l] += ( (plo[l] >> DAT_BITS) << PAD_BITS );	// Need 1 more plo term than phi
		dif_p20_cperms[0] = dif_p20_cperms[31] = 0;	// Handle the leading 0 and its copy separately
	#else
		for(l = 0; l < RADIX; ++l)
		{
			p[l] = l*NDIVR;
			p[l] += ( (p[l] >> DAT_BITS) << PAD_BITS );
		}
		for(l = 0; l < RADIX; ++l)
		{
			o[l] = p[o[l]];
		}
		/* Now that have o-indices, mimic the radix-31-related indexing in the loop to predefine a suitable permutation array for those: */
		jj[0] = RADIX;
		for(l = 1; l < 31; ++l)
		{
			jj[l] = jj[l-1] - 32;
		}
		jj[0] = 0;
		for(int k = 0; k < 32; ++k)
		{
			/* Without the extra auxiliary q-array, the indexing here would be based on the p-array like so:
			RADIX_31_DIF(
				a+j1,p[jj[0]],p[jj[1]],p[jj[2]],p[jj[3]],p[jj[4]],p[jj[5]],p[jj[6]],p[jj[7]],p[jj[8]],p[jj[9]],p[jj[10]],p[jj[11]],p[jj[12]],p[jj[13]],p[jj[14]],p[jj[15]],p[jj[16]],p[jj[17]],p[jj[18]],p[jj[19]],p[jj[20]],p[jj[21]],p[jj[22]],p[jj[23]],p[jj[24]],p[jj[25]],p[jj[26]],p[jj[27]],p[jj[28]],p[jj[29]],p[jj[30]],
				t + k*31
			);
			Want to replace the 31 disparate p[jj[0-30]] indices w/block of 31 contiguous memory locations, so do like so:
			*/
			for(l = 0; l < 31; ++l)
			{
				q[k*31 + l] = p[jj[l]];
				jj[l] -= 31;
				jj[l] += ( (-(int)((uint32)jj[l] >> 31)) & RADIX);	/* (jj[l]-62)%1984; */
			}
		}
	#endif
	}

/*...The radix-992 pass is here.	*/

	for(j = 0; j < NDIVR; j += 2)
	{
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		//j2 = j1+RE_IM_STRIDE;

		/*...gather the needed data (992 64-bit complex, i.e 1984 64-bit reals) and do 32 radix-31 transforms...*/
	/*
	Twiddleless version arranges 32 sets of radix-31 DFT inputs as follows: 0 in upper left corner,
	decrement 32 horizontally and 31 vertically, all resulting indexing modulo 992:
	(we can auto-generate these by compiling test_fft_radix.c with -DTTYPE=0 -DRADIX=992, running
	the resulting executable and snarfing the first set of index-outputs, "DIF/DIT input-scramble array").

	DIF/DIT input-scramble array = [in hex]
		000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,
		3c1,3a1,381,361,341,321,301,2e1,2c1,2a1,281,261,241,221,201,1e1,1c1,1a1,181,161,141,121,101,0e1,0c1,0a1,081,061,041,021,001,
		3a2,382,362,342,322,302,2e2,2c2,2a2,282,262,242,222,202,1e2,1c2,1a2,182,162,142,122,102,0e2,0c2,0a2,082,062,042,022,002,3c2,
		383,363,343,323,303,2e3,2c3,2a3,283,263,243,223,203,1e3,1c3,1a3,183,163,143,123,103,0e3,0c3,0a3,083,063,043,023,003,3c3,3a3,
		364,344,324,304,2e4,2c4,2a4,284,264,244,224,204,1e4,1c4,1a4,184,164,144,124,104,0e4,0c4,0a4,084,064,044,024,004,3c4,3a4,384,
		345,325,305,2e5,2c5,2a5,285,265,245,225,205,1e5,1c5,1a5,185,165,145,125,105,0e5,0c5,0a5,085,065,045,025,005,3c5,3a5,385,365,
		326,306,2e6,2c6,2a6,286,266,246,226,206,1e6,1c6,1a6,186,166,146,126,106,0e6,0c6,0a6,086,066,046,026,006,3c6,3a6,386,366,346,
		307,2e7,2c7,2a7,287,267,247,227,207,1e7,1c7,1a7,187,167,147,127,107,0e7,0c7,0a7,087,067,047,027,007,3c7,3a7,387,367,347,327,
		2e8,2c8,2a8,288,268,248,228,208,1e8,1c8,1a8,188,168,148,128,108,0e8,0c8,0a8,088,068,048,028,008,3c8,3a8,388,368,348,328,308,
		2c9,2a9,289,269,249,229,209,1e9,1c9,1a9,189,169,149,129,109,0e9,0c9,0a9,089,069,049,029,009,3c9,3a9,389,369,349,329,309,2e9,
		2aa,28a,26a,24a,22a,20a,1ea,1ca,1aa,18a,16a,14a,12a,10a,0ea,0ca,0aa,08a,06a,04a,02a,00a,3ca,3aa,38a,36a,34a,32a,30a,2ea,2ca,
		28b,26b,24b,22b,20b,1eb,1cb,1ab,18b,16b,14b,12b,10b,0eb,0cb,0ab,08b,06b,04b,02b,00b,3cb,3ab,38b,36b,34b,32b,30b,2eb,2cb,2ab,
		26c,24c,22c,20c,1ec,1cc,1ac,18c,16c,14c,12c,10c,0ec,0cc,0ac,08c,06c,04c,02c,00c,3cc,3ac,38c,36c,34c,32c,30c,2ec,2cc,2ac,28c,
		24d,22d,20d,1ed,1cd,1ad,18d,16d,14d,12d,10d,0ed,0cd,0ad,08d,06d,04d,02d,00d,3cd,3ad,38d,36d,34d,32d,30d,2ed,2cd,2ad,28d,26d,
		22e,20e,1ee,1ce,1ae,18e,16e,14e,12e,10e,0ee,0ce,0ae,08e,06e,04e,02e,00e,3ce,3ae,38e,36e,34e,32e,30e,2ee,2ce,2ae,28e,26e,24e,
		20f,1ef,1cf,1af,18f,16f,14f,12f,10f,0ef,0cf,0af,08f,06f,04f,02f,00f,3cf,3af,38f,36f,34f,32f,30f,2ef,2cf,2af,28f,26f,24f,22f,
		1f0,1d0,1b0,190,170,150,130,110,0f0,0d0,0b0,090,070,050,030,010,3d0,3b0,390,370,350,330,310,2f0,2d0,2b0,290,270,250,230,210,
		1d1,1b1,191,171,151,131,111,0f1,0d1,0b1,091,071,051,031,011,3d1,3b1,391,371,351,331,311,2f1,2d1,2b1,291,271,251,231,211,1f1,
		1b2,192,172,152,132,112,0f2,0d2,0b2,092,072,052,032,012,3d2,3b2,392,372,352,332,312,2f2,2d2,2b2,292,272,252,232,212,1f2,1d2,
		193,173,153,133,113,0f3,0d3,0b3,093,073,053,033,013,3d3,3b3,393,373,353,333,313,2f3,2d3,2b3,293,273,253,233,213,1f3,1d3,1b3,
		174,154,134,114,0f4,0d4,0b4,094,074,054,034,014,3d4,3b4,394,374,354,334,314,2f4,2d4,2b4,294,274,254,234,214,1f4,1d4,1b4,194,
		155,135,115,0f5,0d5,0b5,095,075,055,035,015,3d5,3b5,395,375,355,335,315,2f5,2d5,2b5,295,275,255,235,215,1f5,1d5,1b5,195,175,
		136,116,0f6,0d6,0b6,096,076,056,036,016,3d6,3b6,396,376,356,336,316,2f6,2d6,2b6,296,276,256,236,216,1f6,1d6,1b6,196,176,156,
		117,0f7,0d7,0b7,097,077,057,037,017,3d7,3b7,397,377,357,337,317,2f7,2d7,2b7,297,277,257,237,217,1f7,1d7,1b7,197,177,157,137,
		0f8,0d8,0b8,098,078,058,038,018,3d8,3b8,398,378,358,338,318,2f8,2d8,2b8,298,278,258,238,218,1f8,1d8,1b8,198,178,158,138,118,
		0d9,0b9,099,079,059,039,019,3d9,3b9,399,379,359,339,319,2f9,2d9,2b9,299,279,259,239,219,1f9,1d9,1b9,199,179,159,139,119,0f9,
		0ba,09a,07a,05a,03a,01a,3da,3ba,39a,37a,35a,33a,31a,2fa,2da,2ba,29a,27a,25a,23a,21a,1fa,1da,1ba,19a,17a,15a,13a,11a,0fa,0da,
		09b,07b,05b,03b,01b,3db,3bb,39b,37b,35b,33b,31b,2fb,2db,2bb,29b,27b,25b,23b,21b,1fb,1db,1bb,19b,17b,15b,13b,11b,0fb,0db,0bb,
		07c,05c,03c,01c,3dc,3bc,39c,37c,35c,33c,31c,2fc,2dc,2bc,29c,27c,25c,23c,21c,1fc,1dc,1bc,19c,17c,15c,13c,11c,0fc,0dc,0bc,09c,
		05d,03d,01d,3dd,3bd,39d,37d,35d,33d,31d,2fd,2dd,2bd,29d,27d,25d,23d,21d,1fd,1dd,1bd,19d,17d,15d,13d,11d,0fd,0dd,0bd,09d,07d,
		03e,01e,3de,3be,39e,37e,35e,33e,31e,2fe,2de,2be,29e,27e,25e,23e,21e,1fe,1de,1be,19e,17e,15e,13e,11e,0fe,0de,0be,09e,07e,05e,
		01f,3df,3bf,39f,37f,35f,33f,31f,2ff,2df,2bf,29f,27f,25f,23f,21f,1ff,1df,1bf,19f,17f,15f,13f,11f,0ff,0df,0bf,09f,07f,05f,03f,

	If we subtract each row's common (mod 32) low p-part from all terms of the row
	as we do in the implementation to reduce the number of index offsets needing to be stored,
	we decrement 32 horizontally and 32 vertically, all (mod 992 = 0x3e0):

		000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020 + p00
		3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000 + p01
		3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0 + p02
		380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0 + p03
		360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380 + p04
		340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360 + p05
		320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340 + p06
		300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320 + p07
		2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300 + p08
		2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0 + p09
		2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0 + p0a
		280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0 + p0b
		260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280 + p0c
		240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260 + p0d
		220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240 + p0e
		200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220 + p0f
		1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200 + p10
		1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0 + p11
		1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0 + p12
		180,160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0 + p13
		160,140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180 + p14
		140,120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160 + p15
		120,100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140 + p16
		100,0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120 + p17
		0e0,0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100 + p18
		0c0,0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0 + p19
		0a0,080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0 + p1a
		080,060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0 + p1b
		060,040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080 + p1c
		040,020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060 + p1d
		020,000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040 + p1e
		000,3c0,3a0,380,360,340,320,300,2e0,2c0,2a0,280,260,240,220,200,1e0,1c0,1a0,180,160,140,120,100,0e0,0c0,0a0,080,060,040,020 + p1f

	In order to make the result amenable to loop-based execution, we need to encode the indices
	to the left & right of the + in easily-computable-index fashion. This scheme has 2 ingredients:
	[Use i as loop idx in comments, although use 'l' in actual code]

	[0] RHS is simply i = 0,..,1f, i.e.plo[i]

	[1] Note that the (mod p20, where 20 = hex) row data are simply circular-shift perms of the basic (row 0) set,
		which we simplify by dividing by p32, i.e. by reducing to indices into a small phi[i] = p20*i (i=0,..,1e) array:

		00,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01 + p00

	For row i, the leftward-shift count needing to be applied = i. Since i runs from 0 thru 1f (= decimal 31),
	by the time we get to the last of the 32 rows we reach a shift count = 31 which matches the number of elements
	in our basic array, thus the last row (mod p32) matches the first (0) row.
	*/
		tptr = t;
		for(l = 0; l < 32; ++l)
		{
		#if USE_COMPACT_OBJ_CODE
			iptr = dif_p20_cperms+l;
			RADIX_31_DIF(
				a+j1+plo[l], iptr,
				(double *)tptr
			);
		#else
			iptr = q + ((l<<5)-l);	// q + l*31
			RADIX_31_DIF(
				a+j1,iptr,
				(double *)tptr
			);
		#endif
			tptr += 31;
		}
	/*...and now do 31 radix-32 transforms:

		tptr = t;
		RADIX_32_DIF(tptr+ 0,31,62,...,961, a+j1,o[000],o[001],o[002],o[003],o[004],o[005],o[006],o[007],o[008],o[009],o[010],o[011],o[012],o[013],o[014],o[015],o[016],o[017],o[018],o[019],o[020],o[021],o[022],o[023],o[024],o[025],o[026],o[027],o[028],o[029],o[030],o[031]);
		RADIX_32_DIF(tptr+ 1,32,63,...,962, a+j1,o[031],...
		...                                        ...
		RADIX_32_DIF(tptr+30,61,92,...,991, a+j1,o[960],o[961],o[962],o[963],o[964],o[965],o[966],o[967],o[968],o[969],o[970],o[971],o[972],o[973],o[974],o[975],o[976],o[977],o[978],o[979],o[980],o[981],o[982],o[983],o[984],o[985],o[986],o[987],o[988],o[989],o[990],o[991]);

		The required output permutation is stored in the o-array.

		In the implementation of the above, the jj-offsets here are raw real-array offsets, so need to be doubled relative to complex-indexing;
		The iptr-offset is an index into the length-992 o-array, whose elements implicitly contain the needed doubling.
	*/
	/* Apr 2014: O-perm in hex form - generated this with a view toward Impl 2-table-style streamlined indexing
		to match other large-radix DFTs in the compact-obj-code framework.
		The required output permutation is, grouping things into blocks of 32, one for each DFT output set:

		000,001,002,003,004,005,006,007,008,009,00a,00b,00c,00d,00e,00f,	010,011,012,013,014,015,016,017,018,019,01a,01b,01c,01d,01e,01f,
		3cf,3ce,3cc,3cd,3c8,3c9,3ca,3cb,3c0,3c1,3c2,3c3,3c4,3c5,3c6,3c7,	3df,3de,3dc,3dd,3d8,3d9,3da,3db,3d0,3d1,3d2,3d3,3d4,3d5,3d6,3d7,
		3b7,3b6,3b4,3b5,3b0,3b1,3b2,3b3,3bf,3be,3bc,3bd,3b8,3b9,3ba,3bb,	3af,3ae,3ac,3ad,3a8,3a9,3aa,3ab,3a0,3a1,3a2,3a3,3a4,3a5,3a6,3a7,
		387,386,384,385,380,381,382,383,38f,38e,38c,38d,388,389,38a,38b,	397,396,394,395,390,391,392,393,39f,39e,39c,39d,398,399,39a,39b,
		37b,37a,378,379,37f,37e,37c,37d,377,376,374,375,370,371,372,373,	367,366,364,365,360,361,362,363,36f,36e,36c,36d,368,369,36a,36b,
		34b,34a,348,349,34f,34e,34c,34d,347,346,344,345,340,341,342,343,	35b,35a,358,359,35f,35e,35c,35d,357,356,354,355,350,351,352,353,
		333,332,330,331,337,336,334,335,33b,33a,338,339,33f,33e,33c,33d,	32b,32a,328,329,32f,32e,32c,32d,327,326,324,325,320,321,322,323,
		303,302,300,301,307,306,304,305,30b,30a,308,309,30f,30e,30c,30d,	313,312,310,311,317,316,314,315,31b,31a,318,319,31f,31e,31c,31d,
		2fd,2fc,2ff,2fe,2fb,2fa,2f8,2f9,2f3,2f2,2f0,2f1,2f7,2f6,2f4,2f5,	2e3,2e2,2e0,2e1,2e7,2e6,2e4,2e5,2eb,2ea,2e8,2e9,2ef,2ee,2ec,2ed,
		2cd,2cc,2cf,2ce,2cb,2ca,2c8,2c9,2c3,2c2,2c0,2c1,2c7,2c6,2c4,2c5,	2dd,2dc,2df,2de,2db,2da,2d8,2d9,2d3,2d2,2d0,2d1,2d7,2d6,2d4,2d5,
		2b5,2b4,2b7,2b6,2b3,2b2,2b0,2b1,2bd,2bc,2bf,2be,2bb,2ba,2b8,2b9,	2ad,2ac,2af,2ae,2ab,2aa,2a8,2a9,2a3,2a2,2a0,2a1,2a7,2a6,2a4,2a5,
		285,284,287,286,283,282,280,281,28d,28c,28f,28e,28b,28a,288,289,	295,294,297,296,293,292,290,291,29d,29c,29f,29e,29b,29a,298,299,
		279,278,27b,27a,27d,27c,27f,27e,275,274,277,276,273,272,270,271,	265,264,267,266,263,262,260,261,26d,26c,26f,26e,26b,26a,268,269,
		249,248,24b,24a,24d,24c,24f,24e,245,244,247,246,243,242,240,241,	259,258,25b,25a,25d,25c,25f,25e,255,254,257,256,253,252,250,251,
		231,230,233,232,235,234,237,236,239,238,23b,23a,23d,23c,23f,23e,	229,228,22b,22a,22d,22c,22f,22e,225,224,227,226,223,222,220,221,
		201,200,203,202,205,204,207,206,209,208,20b,20a,20d,20c,20f,20e,	211,210,213,212,215,214,217,216,219,218,21b,21a,21d,21c,21f,21e,
		1fe,1ff,1fd,1fc,1f9,1f8,1fb,1fa,1f1,1f0,1f3,1f2,1f5,1f4,1f7,1f6,	1e1,1e0,1e3,1e2,1e5,1e4,1e7,1e6,1e9,1e8,1eb,1ea,1ed,1ec,1ef,1ee,
		1ce,1cf,1cd,1cc,1c9,1c8,1cb,1ca,1c1,1c0,1c3,1c2,1c5,1c4,1c7,1c6,	1de,1df,1dd,1dc,1d9,1d8,1db,1da,1d1,1d0,1d3,1d2,1d5,1d4,1d7,1d6,
		1b6,1b7,1b5,1b4,1b1,1b0,1b3,1b2,1be,1bf,1bd,1bc,1b9,1b8,1bb,1ba,	1ae,1af,1ad,1ac,1a9,1a8,1ab,1aa,1a1,1a0,1a3,1a2,1a5,1a4,1a7,1a6,
		186,187,185,184,181,180,183,182,18e,18f,18d,18c,189,188,18b,18a,	196,197,195,194,191,190,193,192,19e,19f,19d,19c,199,198,19b,19a,
		17a,17b,179,178,17e,17f,17d,17c,176,177,175,174,171,170,173,172,	166,167,165,164,161,160,163,162,16e,16f,16d,16c,169,168,16b,16a,
		14a,14b,149,148,14e,14f,14d,14c,146,147,145,144,141,140,143,142,	15a,15b,159,158,15e,15f,15d,15c,156,157,155,154,151,150,153,152,
		132,133,131,130,136,137,135,134,13a,13b,139,138,13e,13f,13d,13c,	12a,12b,129,128,12e,12f,12d,12c,126,127,125,124,121,120,123,122,
		102,103,101,100,106,107,105,104,10a,10b,109,108,10e,10f,10d,10c,	112,113,111,110,116,117,115,114,11a,11b,119,118,11e,11f,11d,11c,
		0fc,0fd,0fe,0ff,0fa,0fb,0f9,0f8,0f2,0f3,0f1,0f0,0f6,0f7,0f5,0f4,	0e2,0e3,0e1,0e0,0e6,0e7,0e5,0e4,0ea,0eb,0e9,0e8,0ee,0ef,0ed,0ec,
		0cc,0cd,0ce,0cf,0ca,0cb,0c9,0c8,0c2,0c3,0c1,0c0,0c6,0c7,0c5,0c4,	0dc,0dd,0de,0df,0da,0db,0d9,0d8,0d2,0d3,0d1,0d0,0d6,0d7,0d5,0d4,
		0b4,0b5,0b6,0b7,0b2,0b3,0b1,0b0,0bc,0bd,0be,0bf,0ba,0bb,0b9,0b8,	0ac,0ad,0ae,0af,0aa,0ab,0a9,0a8,0a2,0a3,0a1,0a0,0a6,0a7,0a5,0a4,
		084,085,086,087,082,083,081,080,08c,08d,08e,08f,08a,08b,089,088,	094,095,096,097,092,093,091,090,09c,09d,09e,09f,09a,09b,099,098,
		078,079,07a,07b,07c,07d,07e,07f,074,075,076,077,072,073,071,070,	064,065,066,067,062,063,061,060,06c,06d,06e,06f,06a,06b,069,068,
		048,049,04a,04b,04c,04d,04e,04f,044,045,046,047,042,043,041,040,	058,059,05a,05b,05c,05d,05e,05f,054,055,056,057,052,053,051,050,
		030,031,032,033,034,035,036,037,038,039,03a,03b,03c,03d,03e,03f,	028,029,02a,02b,02c,02d,02e,02f,024,025,026,027,022,023,021,020,

		row:																			Indices of the 2 phi offsets per row:
																						+ := phi2 = phi1+p010, - := phi1 = phi2+p010
		00	p000 + p[0123456789abcdef]	p010 + p[0123456789abcdef]   p000 + [0]	p010 + [0]   00,01	+	0 = (row+1)/2 [compare to []-pattern index]
		01	p3c0 + p[fecd89ab01234567]	p3d0 + p[fecd89ab01234567]   p3c0 + [1]	p3d0 + [1]   3c,3d	+	1	<*** odd  row (and 0): []-idx = (row+1)/2 for both []
		02	p3b0 + p[76450123fecd89ab]	p3a0 + p[fecd89ab01234567]   p3b0 + [2]	p3a0 + [1]   3b,3a	-	1	<*** even row: 1st []-idx = 1+(row+1)/2, 2nd []-idx = (row+1)/2
		03	p380 + p[76450123fecd89ab]	p390 + p[76450123fecd89ab]   p380 + [2]	p390 + [2]   38,39	+	2
		04	p370 + p[ba89fecd76450123]	p360 + p[76450123fecd89ab]   p370 + [3]	p360 + [2]   37,36	-	2
		05	p340 + p[ba89fecd76450123]	p350 + p[ba89fecd76450123]   p340 + [3]	p350 + [3]   34,35	+	3
		06	p330 + p[32017645ba89fecd]	p320 + p[ba89fecd76450123]   p330 + [4]	p320 + [3]   33,32	-	3
		07	p300 + p[32017645ba89fecd]	p310 + p[32017645ba89fecd]   p300 + [4]	p310 + [4]   30,31	+	4
		08	p2f0 + p[dcfeba8932017645]	p2e0 + p[32017645ba89fecd]   p2f0 + [5]	p2e0 + [4]   2f,2e	-	4
		09	p2c0 + p[dcfeba8932017645]	p2d0 + p[dcfeba8932017645]   p2c0 + [5]	p2d0 + [5]   2c,2d	+	5
		0a	p2b0 + p[54763201dcfeba89]	p2a0 + p[dcfeba8932017645]   p2b0 + [6]	p2a0 + [5]   2b,2a	-	5
		0b	p280 + p[54763201dcfeba89]	p290 + p[54763201dcfeba89]   p280 + [6]	p290 + [6]   28,29	+	6
		0c	p270 + p[98badcfe54763201]	p260 + p[54763201dcfeba89]   p270 + [7]	p260 + [6]   27,26	-	6
		0d	p240 + p[98badcfe54763201]	p250 + p[98badcfe54763201]   p240 + [7]	p250 + [7]   24,25	+	7
		0e	p230 + p[1032547698badcfe]	p220 + p[98badcfe54763201]   p230 + [8]	p220 + [7]   23,22	-	7
		0f	p200 + p[1032547698badcfe]	p210 + p[1032547698badcfe] = p200 + [8]	p210 + [8] = 20,21	+	8
		10	p1f0 + p[efdc98ba10325476]	p1e0 + p[1032547698badcfe]   p1f0 + [9]	p1e0 + [8]   1f,1e	-	8
		11	p1c0 + p[efdc98ba10325476]	p1d0 + p[efdc98ba10325476]   p1c0 + [9]	p1d0 + [9]   1c,1d	+	9
		12	p1b0 + p[67541032efdc98ba]	p1a0 + p[efdc98ba10325476]   p1b0 + [A]	p1a0 + [9]   1b,1a	-	9
		13	p180 + p[67541032efdc98ba]	p190 + p[67541032efdc98ba]   p180 + [A]	p190 + [A]   18,19	+	A
		14	p170 + p[ab98efdc67541032]	p160 + p[67541032efdc98ba]   p170 + [B]	p160 + [A]   17,16	-	A
		15	p140 + p[ab98efdc67541032]	p150 + p[ab98efdc67541032]   p140 + [B]	p150 + [B]   14,15	+	B
		16	p130 + p[23106754ab98efdc]	p120 + p[ab98efdc67541032]   p130 + [C]	p120 + [B]   13,12	-	B
		17	p100 + p[23106754ab98efdc]	p110 + p[23106754ab98efdc]   p100 + [C]	p110 + [C]   10,11	+	C
		18	p0f0 + p[cdefab9823106754]	p0e0 + p[23106754ab98efdc]   p0f0 + [D]	p0e0 + [C]   0f,0e	-	C
		19	p0c0 + p[cdefab9823106754]	p0d0 + p[cdefab9823106754]   p0c0 + [D]	p0d0 + [D]   0c,0d	+	D
		1a	p0b0 + p[45672310cdefab98]	p0a0 + p[cdefab9823106754]   p0b0 + [E]	p0a0 + [D]   0b,0a	-	D
		1b	p080 + p[45672310cdefab98]	p090 + p[45672310cdefab98]   p080 + [E]	p090 + [E]   08,09	+	E
		1c	p070 + p[89abcdef45672310]	p060 + p[45672310cdefab98]   p070 + [F]	p060 + [E]   07,06	-	E
		1d	p040 + p[89abcdef45672310]	p050 + p[89abcdef45672310]   p040 + [F]	p050 + [F]   04,05	+	F
		1e	p030 + p[0123456789abcdef]	p020 + p[89abcdef45672310]   p030 + [0]	p020 + [F]   03,02	-	F	<*** Last row: need to do []-indexing (mod 16), thus 1st []-idx = 16 (mod 16) = 0

		Thus, aside from the 0-row pattern (which is generally 'special' for these kinds of DFTs), the 2 phi
		offsets per row follow a simple alternating pattern:

			idx = ((31-i)<<1) & (-(i>0));	// Index into phi[] for the 2 phi-offsets of this row = 00(3e & 0),3c,3a,...,4,2
			is_odd = (i&1) + (i==0); is_even = (~is_odd)&1;	// 0-row phi-pair relative offsets behave like odd-index row
			phi1 + phi[idx+is_even];  // = phi[idx]+p10 for even-idx rows, = phi[idx]     for odd]
			phi2 + phi[idx+is_odd ];  // = phi[idx]     for even-idx rows, = phi[idx]+p10 for odd]
		*/
		for(int k = 0; k < 32; ++k) {
			jj[k] = ((k<<5)-k)<<1;	// DFT macro takes *real*-double inputs, thus compute doubled offsets k*62
		}
		iptr = o;
		for(l = 0; l < 31; ++l) {
		#if USE_COMPACT_OBJ_CODE
			idx = (31-l) & (-(l>0));	// Base-offset Index into phi[] for the current row = 00,3c,3a,...,4,2
										// (Since our phi array here is in mults of 32, phi = 0x3c0 ==> idx = 3c/2 = 1e)
			is_odd = (l&1) + (l==0); is_even = (~is_odd)&1;	// 0-row phi-pair relative offsets behave like odd-index row
			// Compute index in the perm16 array needed for current row:
			pidx = (l+1)>>1;	//<*** odd  row (and 0): []-idx = (row+1)/2 for both []
								//<*** even row: 1st []-idx = 1+(row+1)/2, 2nd []-idx = (row+1)/2
			// First 16 offsets:
			jp = phi[idx] + plo[is_even<<4];  // = phi[idx]+p10 for even-idx rows, = phi[idx] for odd]
			i64 = dif_perm16[(pidx+is_even)&0xf];
			// And lastly add in the low parts - p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];		o[0x0] = jp + k0;
			k1 = plo[(i64 >> 56)&0xf];		o[0x1] = jp + k1;
			k2 = plo[(i64 >> 52)&0xf];		o[0x2] = jp + k2;
			k3 = plo[(i64 >> 48)&0xf];		o[0x3] = jp + k3;
			k4 = plo[(i64 >> 44)&0xf];		o[0x4] = jp + k4;
			k5 = plo[(i64 >> 40)&0xf];		o[0x5] = jp + k5;
			k6 = plo[(i64 >> 36)&0xf];		o[0x6] = jp + k6;
			k7 = plo[(i64 >> 32)&0xf];		o[0x7] = jp + k7;
			k8 = plo[(i64 >> 28)&0xf];		o[0x8] = jp + k8;
			k9 = plo[(i64 >> 24)&0xf];		o[0x9] = jp + k9;
			ka = plo[(i64 >> 20)&0xf];		o[0xa] = jp + ka;
			kb = plo[(i64 >> 16)&0xf];		o[0xb] = jp + kb;
			kc = plo[(i64 >> 12)&0xf];		o[0xc] = jp + kc;
			kd = plo[(i64 >>  8)&0xf];		o[0xd] = jp + kd;
			ke = plo[(i64 >>  4)&0xf];		o[0xe] = jp + ke;
			kf = plo[(i64      )&0xf];		o[0xf] = jp + kf;
			// Second 16 offsets:
			jp = phi[idx] + plo[is_odd<<4];  // = phi[idx] for even-idx rows, = phi[idx]+p10 for odd]
			i64 = dif_perm16[pidx];	// 2nd []-idx = (row+1)/2 for both even and odd rows
			// And lastly add in the low parts - p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];		o[0x10] = jp + k0;
			k1 = plo[(i64 >> 56)&0xf];		o[0x11] = jp + k1;
			k2 = plo[(i64 >> 52)&0xf];		o[0x12] = jp + k2;
			k3 = plo[(i64 >> 48)&0xf];		o[0x13] = jp + k3;
			k4 = plo[(i64 >> 44)&0xf];		o[0x14] = jp + k4;
			k5 = plo[(i64 >> 40)&0xf];		o[0x15] = jp + k5;
			k6 = plo[(i64 >> 36)&0xf];		o[0x16] = jp + k6;
			k7 = plo[(i64 >> 32)&0xf];		o[0x17] = jp + k7;
			k8 = plo[(i64 >> 28)&0xf];		o[0x18] = jp + k8;
			k9 = plo[(i64 >> 24)&0xf];		o[0x19] = jp + k9;
			ka = plo[(i64 >> 20)&0xf];		o[0x1a] = jp + ka;
			kb = plo[(i64 >> 16)&0xf];		o[0x1b] = jp + kb;
			kc = plo[(i64 >> 12)&0xf];		o[0x1c] = jp + kc;
			kd = plo[(i64 >>  8)&0xf];		o[0x1d] = jp + kd;
			ke = plo[(i64 >>  4)&0xf];		o[0x1e] = jp + ke;
			kf = plo[(i64      )&0xf];		o[0x1f] = jp + kf;
			RADIX_32_DIF(
				(double *)t,jj,   1,
				a+j1,o,RE_IM_STRIDE
			);
		#else
			RADIX_32_DIF(
				(double *)t,jj,   1,
				a+j1,iptr,RE_IM_STRIDE	// *** Apr 2014: Added data strdes needed by latest version of the radix-32 DFT macros
			);
			iptr += 32;
		#endif
			for(int k = 0; k < 32; ++k) {
				jj[k] += 2;
			}
		}
	}
	/* Totals: 32*radix31 + 31*radix32 = 32*(658 FADD, 234 FMUL) + 31*(376 FADD, 88 FMUL)
										= 32712 FADD, 10216 FMUL, or 16.5 FADD, 5.1 FMUL per real input.
	Compare this to radix 1024, which we can do via 32*radix32 + 32*[radix32 with twiddle-muls];
	we need 31 complex twiddle-muls for each of the 32 latter DFTs, each costing 2 FADD and 4 FMUL,
	thus add 62 FADD, 124 FMUL to each of the step-2 DFTs, thus our radix-1024 costs
	32*(376 FADD, 88 FMUL) + 32*(438 FADD, 212 FMUL)
										= 26048 FADD, 9600 FMUL, or 12.7 FADD, 4.7 FMUL per real input.
	Thus, even though the radix-31 DFT alone costs 10.6 FADD, 3.8 FMUL per real input which is nearly double
	the per-point FADD cost and triple the FMUL cost of a twiddleless complex radix-32 DFT [5.9 FADD, 1.4 FMUL per point],
	when we combine such prime-radix DFTs with power-of-2 ones in an optimized index-permuted twiddleless
	fashion, the resulting normalized opcounts quickly become much more competitive with those of a similar-sized
	power-of-2 transform: Our radix-992 transform has just 30% more floating-point adds and 10% more muls per point
	than does an similarly wel-optimized radix-1024 transform.

	But especially the add count is still higher than is desirable, so alternatively consider a length-1008 = 63*16 DFT.
	We assemble the radix-63 via a twiddleless combination of 9*radix7 + 7*radix9, with resulting opcount
	9*(72 FADD, 16 FMUL) + 7*(80 FADD, 40 FMUL) = 1208 FADD, 424 FMUL. This is combined twiddlelessly with a set
	of 63 radix-16 DFTs, for a total opcount of 16*(1208 FADD, 424 FMUL) + 63*(144 FADD, 24 FMUL)
	= 28400 FADD, 8296 FMUL, which is only 9% more floating-point adds and 14% fewer muls than our radix-1024,
	and moreover results in an appreciably lower level of roundoff error than does radix 992, which is especially crucial for F33.
	*/
}

/***************/

void radix992_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-992 = 31x32 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix992_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
*/
	int l, j,j1/* ,j2 */,jp, jj[32], *iptr;
	static int NDIVR,first_entry=TRUE;
	struct complex t[RADIX], *tptr;
#if USE_COMPACT_OBJ_CODE
	// Need storage for circular-shifts perms of a basic 31-vector, with shift count in [0,31] that means 2*31 elts:
	static int plo[32],phi[62];	// No need for separate dit_p20_cperms[] array; just need a doubled-sequence version of phi
	int idx,pidx,mask,lshift, is_even,is_odd, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf, o[32];	// o[] stores o-address offsets for current radix-32 DFT in the 31x-loop
	uint64 i64;
	// Low parts [p0-f] of output-index perms - need 64x4-bits for each radix-64 DFT:
	const uint64 dit_perm16[16] = {
		0x01327654fedcba98ull,0xfedcba9876543210ull,0x76543210ba98dcefull,0xba98dcef32105467ull,
		0x32105467dcef98abull,0xdcef98ab54671023ull,0x5467102398abefcdull,0x98abefcd10236745ull,
		0x10236745efcdab89ull,0xefcdab8967452301ull,0x67452301ab89cdfeull,0xab89cdfe23014576ull,
		0x23014576cdfe89baull,0xcdfe89ba45760132ull,0x4576013289bafedcull,0x89bafedc01327654ull
	};
#else
	static int p[RADIX], q[RADIX] = {0,1,3,2,7,6,5,4,15,14,13,12,11,10,9,8,31,30,29,28,27,26,25,24,
	23,22,21,20,19,18,17,16,975,974,973,972,971,970,969,968,967,966,965,964,963,962,961,960,983,982,981,
	980,979,978,977,976,987,986,985,984,989,988,990,991,951,950,949,948,947,946,945,944,955,954,953,952,
	957,956,958,959,935,934,933,932,931,930,929,928,939,938,937,936,941,940,942,943,903,902,901,900,899,
	898,897,896,907,906,905,904,909,908,910,911,923,922,921,920,925,924,926,927,915,914,913,912,917,916,
	918,919,891,890,889,888,893,892,894,895,883,882,881,880,885,884,886,887,875,874,873,872,877,876,878,
	879,867,866,865,864,869,868,870,871,843,842,841,840,845,844,846,847,835,834,833,832,837,836,838,839,
	851,850,849,848,853,852,854,855,861,860,862,863,857,856,858,859,819,818,817,816,821,820,822,823,829,
	828,830,831,825,824,826,827,803,802,801,800,805,804,806,807,813,812,814,815,809,808,810,811,771,770,
	769,768,773,772,774,775,781,780,782,783,777,776,778,779,797,796,798,799,793,792,794,795,789,788,790,
	791,785,784,786,787,765,764,766,767,761,760,762,763,757,756,758,759,753,752,754,755,749,748,750,751,
	745,744,746,747,741,740,742,743,737,736,738,739,717,716,718,719,713,712,714,715,709,708,710,711,705,
	704,706,707,725,724,726,727,721,720,722,723,729,728,730,731,734,735,732,733,693,692,694,695,689,688,
	690,691,697,696,698,699,702,703,700,701,677,676,678,679,673,672,674,675,681,680,682,683,686,687,684,
	685,645,644,646,647,641,640,642,643,649,648,650,651,654,655,652,653,665,664,666,667,670,671,668,669,
	657,656,658,659,662,663,660,661,633,632,634,635,638,639,636,637,625,624,626,627,630,631,628,629,617,
	616,618,619,622,623,620,621,609,608,610,611,614,615,612,613,585,584,586,587,590,591,588,589,577,576,
	578,579,582,583,580,581,593,592,594,595,598,599,596,597,606,607,604,605,602,603,600,601,561,560,562,
	563,566,567,564,565,574,575,572,573,570,571,568,569,545,544,546,547,550,551,548,549,558,559,556,557,
	554,555,552,553,513,512,514,515,518,519,516,517,526,527,524,525,522,523,520,521,542,543,540,541,538,
	539,536,537,534,535,532,533,530,531,528,529,510,511,508,509,506,507,504,505,502,503,500,501,498,499,
	496,497,494,495,492,493,490,491,488,489,486,487,484,485,482,483,480,481,462,463,460,461,458,459,456,
	457,454,455,452,453,450,451,448,449,470,471,468,469,466,467,464,465,474,475,472,473,476,477,479,478,
	438,439,436,437,434,435,432,433,442,443,440,441,444,445,447,446,422,423,420,421,418,419,416,417,426,
	427,424,425,428,429,431,430,390,391,388,389,386,387,384,385,394,395,392,393,396,397,399,398,410,411,
	408,409,412,413,415,414,402,403,400,401,404,405,407,406,378,379,376,377,380,381,383,382,370,371,368,
	369,372,373,375,374,362,363,360,361,364,365,367,366,354,355,352,353,356,357,359,358,330,331,328,329,
	332,333,335,334,322,323,320,321,324,325,327,326,338,339,336,337,340,341,343,342,348,349,351,350,344,
	345,347,346,306,307,304,305,308,309,311,310,316,317,319,318,312,313,315,314,290,291,288,289,292,293,
	295,294,300,301,303,302,296,297,299,298,258,259,256,257,260,261,263,262,268,269,271,270,264,265,267,
	266,284,285,287,286,280,281,283,282,276,277,279,278,272,273,275,274,252,253,255,254,248,249,251,250,
	244,245,247,246,240,241,243,242,236,237,239,238,232,233,235,234,228,229,231,230,224,225,227,226,204,
	205,207,206,200,201,203,202,196,197,199,198,192,193,195,194,212,213,215,214,208,209,211,210,216,217,
	219,218,223,222,221,220,180,181,183,182,176,177,179,178,184,185,187,186,191,190,189,188,164,165,167,
	166,160,161,163,162,168,169,171,170,175,174,173,172,132,133,135,134,128,129,131,130,136,137,139,138,
	143,142,141,140,152,153,155,154,159,158,157,156,144,145,147,146,151,150,149,148,120,121,123,122,127,
	126,125,124,112,113,115,114,119,118,117,116,104,105,107,106,111,110,109,108,96,97,99,98,103,102,101,
	100,72,73,75,74,79,78,77,76,64,65,67,66,71,70,69,68,80,81,83,82,87,86,85,84,95,94,93,92,91,90,89,88,
	48,49,51,50,55,54,53,52,63,62,61,60,59,58,57,56,32,33,35,34,39,38,37,36,47,46,45,44,43,42,41,40};
	static int o[RADIX] = {
	 0, 32, 64, 96,128,160,192,224,256,288,320,352,384,416,448,480,512,544,576,608,640,672,704,736,768,800,832,864,896,928,960,
	31, 63, 95,127,159,191,223,255,287,319,351,383,415,447,479,511,543,575,607,639,671,703,735,767,799,831,863,895,927,959,991,
	62, 94,126,158,190,222,254,286,318,350,382,414,446,478,510,542,574,606,638,670,702,734,766,798,830,862,894,926,958,990, 30,
	93,125,157,189,221,253,285,317,349,381,413,445,477,509,541,573,605,637,669,701,733,765,797,829,861,893,925,957,989, 29, 61,
	124,156,188,220,252,284,316,348,380,412,444,476,508,540,572,604,636,668,700,732,764,796,828,860,892,924,956,988, 28, 60, 92,
	155,187,219,251,283,315,347,379,411,443,475,507,539,571,603,635,667,699,731,763,795,827,859,891,923,955,987, 27, 59, 91,123,
	186,218,250,282,314,346,378,410,442,474,506,538,570,602,634,666,698,730,762,794,826,858,890,922,954,986, 26, 58, 90,122,154,
	217,249,281,313,345,377,409,441,473,505,537,569,601,633,665,697,729,761,793,825,857,889,921,953,985, 25, 57, 89,121,153,185,
	248,280,312,344,376,408,440,472,504,536,568,600,632,664,696,728,760,792,824,856,888,920,952,984, 24, 56, 88,120,152,184,216,
	279,311,343,375,407,439,471,503,535,567,599,631,663,695,727,759,791,823,855,887,919,951,983, 23, 55, 87,119,151,183,215,247,
	310,342,374,406,438,470,502,534,566,598,630,662,694,726,758,790,822,854,886,918,950,982, 22, 54, 86,118,150,182,214,246,278,
	341,373,405,437,469,501,533,565,597,629,661,693,725,757,789,821,853,885,917,949,981, 21, 53, 85,117,149,181,213,245,277,309,
	372,404,436,468,500,532,564,596,628,660,692,724,756,788,820,852,884,916,948,980, 20, 52, 84,116,148,180,212,244,276,308,340,
	403,435,467,499,531,563,595,627,659,691,723,755,787,819,851,883,915,947,979, 19, 51, 83,115,147,179,211,243,275,307,339,371,
	434,466,498,530,562,594,626,658,690,722,754,786,818,850,882,914,946,978, 18, 50, 82,114,146,178,210,242,274,306,338,370,402,
	465,497,529,561,593,625,657,689,721,753,785,817,849,881,913,945,977, 17, 49, 81,113,145,177,209,241,273,305,337,369,401,433,
	496,528,560,592,624,656,688,720,752,784,816,848,880,912,944,976, 16, 48, 80,112,144,176,208,240,272,304,336,368,400,432,464,
	527,559,591,623,655,687,719,751,783,815,847,879,911,943,975, 15, 47, 79,111,143,175,207,239,271,303,335,367,399,431,463,495,
	558,590,622,654,686,718,750,782,814,846,878,910,942,974, 14, 46, 78,110,142,174,206,238,270,302,334,366,398,430,462,494,526,
	589,621,653,685,717,749,781,813,845,877,909,941,973, 13, 45, 77,109,141,173,205,237,269,301,333,365,397,429,461,493,525,557,
	620,652,684,716,748,780,812,844,876,908,940,972, 12, 44, 76,108,140,172,204,236,268,300,332,364,396,428,460,492,524,556,588,
	651,683,715,747,779,811,843,875,907,939,971, 11, 43, 75,107,139,171,203,235,267,299,331,363,395,427,459,491,523,555,587,619,
	682,714,746,778,810,842,874,906,938,970, 10, 42, 74,106,138,170,202,234,266,298,330,362,394,426,458,490,522,554,586,618,650,
	713,745,777,809,841,873,905,937,969,  9, 41, 73,105,137,169,201,233,265,297,329,361,393,425,457,489,521,553,585,617,649,681,
	744,776,808,840,872,904,936,968,  8, 40, 72,104,136,168,200,232,264,296,328,360,392,424,456,488,520,552,584,616,648,680,712,
	775,807,839,871,903,935,967,  7, 39, 71,103,135,167,199,231,263,295,327,359,391,423,455,487,519,551,583,615,647,679,711,743,
	806,838,870,902,934,966,  6, 38, 70,102,134,166,198,230,262,294,326,358,390,422,454,486,518,550,582,614,646,678,710,742,774,
	837,869,901,933,965,  5, 37, 69,101,133,165,197,229,261,293,325,357,389,421,453,485,517,549,581,613,645,677,709,741,773,805,
	868,900,932,964,  4, 36, 68,100,132,164,196,228,260,292,324,356,388,420,452,484,516,548,580,612,644,676,708,740,772,804,836,
	899,931,963,  3, 35, 67, 99,131,163,195,227,259,291,323,355,387,419,451,483,515,547,579,611,643,675,707,739,771,803,835,867,
	930,962,  2, 34, 66, 98,130,162,194,226,258,290,322,354,386,418,450,482,514,546,578,610,642,674,706,738,770,802,834,866,898,
	961,  1, 33, 65, 97,129,161,193,225,257,289,321,353,385,417,449,481,513,545,577,609,641,673,705,737,769,801,833,865,897,929};
#endif

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
	#if USE_COMPACT_OBJ_CODE
		// Init plo,phi arrays - No need for separate dit_p20_cperms[] array; just need a doubled-sequence version of phi:
		for(l = 0; l < 31; l++)
		{
			plo[l] = l*NDIVR; phi[l] = plo[l]<<5;
			plo[l] += ( (plo[l] >> DAT_BITS) << PAD_BITS );
			phi[l] += ( (phi[l] >> DAT_BITS) << PAD_BITS );
			phi[31+l] = phi[l];
		}	plo[l] = l*NDIVR;	plo[l] += ( (plo[l] >> DAT_BITS) << PAD_BITS );	// Need 1 more plo term than phi
	#else
		for(j = 0; j < RADIX; ++j)
		{
			p[j] = j*NDIVR;
			p[j] = p[j] + ( (p[j] >> DAT_BITS) << PAD_BITS );
		}
		for(j = 0; j < RADIX; ++j)
		{
			q[j] = p[q[j]];
			o[j] = p[o[j]];
		}
	#endif
	}

/*...The radix-992 pass is here.	*/

	for(j = 0; j < NDIVR; j += 2)
	{
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		//j2 = j1+RE_IM_STRIDE;
	/*
	Gather the needed data (992 64-bit complex) and do 31 radix-32 transforms:

	Twiddleless version uses same linear-index-vector-form permutation as in DIF -
	Remember, inputs to DIT are bit-reversed, so use output of test_fft_radix() with
	TTYPE=0 to auto-generate needed input-index permutation:

	Combined DIT input-scramble array = [in hex]
		000,001,003,002,007,006,005,004,00f,00e,00d,00c,00b,00a,009,008,01f,01e,01d,01c,01b,01a,019,018,017,016,015,014,013,012,011,010,
		3cf,3ce,3cd,3cc,3cb,3ca,3c9,3c8,3c7,3c6,3c5,3c4,3c3,3c2,3c1,3c0,3d7,3d6,3d5,3d4,3d3,3d2,3d1,3d0,3db,3da,3d9,3d8,3dd,3dc,3de,3df,
		3b7,3b6,3b5,3b4,3b3,3b2,3b1,3b0,3bb,3ba,3b9,3b8,3bd,3bc,3be,3bf,3a7,3a6,3a5,3a4,3a3,3a2,3a1,3a0,3ab,3aa,3a9,3a8,3ad,3ac,3ae,3af,
		387,386,385,384,383,382,381,380,38b,38a,389,388,38d,38c,38e,38f,39b,39a,399,398,39d,39c,39e,39f,393,392,391,390,395,394,396,397,
		37b,37a,379,378,37d,37c,37e,37f,373,372,371,370,375,374,376,377,36b,36a,369,368,36d,36c,36e,36f,363,362,361,360,365,364,366,367,
		34b,34a,349,348,34d,34c,34e,34f,343,342,341,340,345,344,346,347,353,352,351,350,355,354,356,357,35d,35c,35e,35f,359,358,35a,35b,
		333,332,331,330,335,334,336,337,33d,33c,33e,33f,339,338,33a,33b,323,322,321,320,325,324,326,327,32d,32c,32e,32f,329,328,32a,32b,
		303,302,301,300,305,304,306,307,30d,30c,30e,30f,309,308,30a,30b,31d,31c,31e,31f,319,318,31a,31b,315,314,316,317,311,310,312,313,
		2fd,2fc,2fe,2ff,2f9,2f8,2fa,2fb,2f5,2f4,2f6,2f7,2f1,2f0,2f2,2f3,2ed,2ec,2ee,2ef,2e9,2e8,2ea,2eb,2e5,2e4,2e6,2e7,2e1,2e0,2e2,2e3,
		2cd,2cc,2ce,2cf,2c9,2c8,2ca,2cb,2c5,2c4,2c6,2c7,2c1,2c0,2c2,2c3,2d5,2d4,2d6,2d7,2d1,2d0,2d2,2d3,2d9,2d8,2da,2db,2de,2df,2dc,2dd,
		2b5,2b4,2b6,2b7,2b1,2b0,2b2,2b3,2b9,2b8,2ba,2bb,2be,2bf,2bc,2bd,2a5,2a4,2a6,2a7,2a1,2a0,2a2,2a3,2a9,2a8,2aa,2ab,2ae,2af,2ac,2ad,
		285,284,286,287,281,280,282,283,289,288,28a,28b,28e,28f,28c,28d,299,298,29a,29b,29e,29f,29c,29d,291,290,292,293,296,297,294,295,
		279,278,27a,27b,27e,27f,27c,27d,271,270,272,273,276,277,274,275,269,268,26a,26b,26e,26f,26c,26d,261,260,262,263,266,267,264,265,
		249,248,24a,24b,24e,24f,24c,24d,241,240,242,243,246,247,244,245,251,250,252,253,256,257,254,255,25e,25f,25c,25d,25a,25b,258,259,
		231,230,232,233,236,237,234,235,23e,23f,23c,23d,23a,23b,238,239,221,220,222,223,226,227,224,225,22e,22f,22c,22d,22a,22b,228,229,
		201,200,202,203,206,207,204,205,20e,20f,20c,20d,20a,20b,208,209,21e,21f,21c,21d,21a,21b,218,219,216,217,214,215,212,213,210,211,
		1fe,1ff,1fc,1fd,1fa,1fb,1f8,1f9,1f6,1f7,1f4,1f5,1f2,1f3,1f0,1f1,1ee,1ef,1ec,1ed,1ea,1eb,1e8,1e9,1e6,1e7,1e4,1e5,1e2,1e3,1e0,1e1,
		1ce,1cf,1cc,1cd,1ca,1cb,1c8,1c9,1c6,1c7,1c4,1c5,1c2,1c3,1c0,1c1,1d6,1d7,1d4,1d5,1d2,1d3,1d0,1d1,1da,1db,1d8,1d9,1dc,1dd,1df,1de,
		1b6,1b7,1b4,1b5,1b2,1b3,1b0,1b1,1ba,1bb,1b8,1b9,1bc,1bd,1bf,1be,1a6,1a7,1a4,1a5,1a2,1a3,1a0,1a1,1aa,1ab,1a8,1a9,1ac,1ad,1af,1ae,
		186,187,184,185,182,183,180,181,18a,18b,188,189,18c,18d,18f,18e,19a,19b,198,199,19c,19d,19f,19e,192,193,190,191,194,195,197,196,
		17a,17b,178,179,17c,17d,17f,17e,172,173,170,171,174,175,177,176,16a,16b,168,169,16c,16d,16f,16e,162,163,160,161,164,165,167,166,
		14a,14b,148,149,14c,14d,14f,14e,142,143,140,141,144,145,147,146,152,153,150,151,154,155,157,156,15c,15d,15f,15e,158,159,15b,15a,
		132,133,130,131,134,135,137,136,13c,13d,13f,13e,138,139,13b,13a,122,123,120,121,124,125,127,126,12c,12d,12f,12e,128,129,12b,12a,
		102,103,100,101,104,105,107,106,10c,10d,10f,10e,108,109,10b,10a,11c,11d,11f,11e,118,119,11b,11a,114,115,117,116,110,111,113,112,
		0fc,0fd,0ff,0fe,0f8,0f9,0fb,0fa,0f4,0f5,0f7,0f6,0f0,0f1,0f3,0f2,0ec,0ed,0ef,0ee,0e8,0e9,0eb,0ea,0e4,0e5,0e7,0e6,0e0,0e1,0e3,0e2,
		0cc,0cd,0cf,0ce,0c8,0c9,0cb,0ca,0c4,0c5,0c7,0c6,0c0,0c1,0c3,0c2,0d4,0d5,0d7,0d6,0d0,0d1,0d3,0d2,0d8,0d9,0db,0da,0df,0de,0dd,0dc,
		0b4,0b5,0b7,0b6,0b0,0b1,0b3,0b2,0b8,0b9,0bb,0ba,0bf,0be,0bd,0bc,0a4,0a5,0a7,0a6,0a0,0a1,0a3,0a2,0a8,0a9,0ab,0aa,0af,0ae,0ad,0ac,
		084,085,087,086,080,081,083,082,088,089,08b,08a,08f,08e,08d,08c,098,099,09b,09a,09f,09e,09d,09c,090,091,093,092,097,096,095,094,
		078,079,07b,07a,07f,07e,07d,07c,070,071,073,072,077,076,075,074,068,069,06b,06a,06f,06e,06d,06c,060,061,063,062,067,066,065,064,
		048,049,04b,04a,04f,04e,04d,04c,040,041,043,042,047,046,045,044,050,051,053,052,057,056,055,054,05f,05e,05d,05c,05b,05a,059,058,
		030,031,033,032,037,036,035,034,03f,03e,03d,03c,03b,03a,039,038,020,021,023,022,027,026,025,024,02f,02e,02d,02c,02b,02a,029,028,
	=
		00,01,03,02,07,06,05,04,0f,0e,0d,0c,0b,0a,09,08,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10   p000   phi[00]
		0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,00,17,16,15,14,13,12,11,10,1b,1a,19,18,1d,1c,1e,1f   p3c0   phi[1e]
		17,16,15,14,13,12,11,10,1b,1a,19,18,1d,1c,1e,1f,07,06,05,04,03,02,01,00,0b,0a,09,08,0d,0c,0e,0f   p3a0   phi[1d]
		07,06,05,04,03,02,01,00,0b,0a,09,08,0d,0c,0e,0f,1b,1a,19,18,1d,1c,1e,1f,13,12,11,10,15,14,16,17   p380   phi[1c]
		1b,1a,19,18,1d,1c,1e,1f,13,12,11,10,15,14,16,17,0b,0a,09,08,0d,0c,0e,0f,03,02,01,00,05,04,06,07   p360   phi[1b]
		0b,0a,09,08,0d,0c,0e,0f,03,02,01,00,05,04,06,07,13,12,11,10,15,14,16,17,1d,1c,1e,1f,19,18,1a,1b   p340   phi[1a]
		13,12,11,10,15,14,16,17,1d,1c,1e,1f,19,18,1a,1b,03,02,01,00,05,04,06,07,0d,0c,0e,0f,09,08,0a,0b   p320   phi[19]
		03,02,01,00,05,04,06,07,0d,0c,0e,0f,09,08,0a,0b,1d,1c,1e,1f,19,18,1a,1b,15,14,16,17,11,10,12,13   p300   phi[18]
		1d,1c,1e,1f,19,18,1a,1b,15,14,16,17,11,10,12,13,0d,0c,0e,0f,09,08,0a,0b,05,04,06,07,01,00,02,03   p2e0   phi[17]
		0d,0c,0e,0f,09,08,0a,0b,05,04,06,07,01,00,02,03,15,14,16,17,11,10,12,13,19,18,1a,1b,1e,1f,1c,1d   p2c0   phi[16]
		15,14,16,17,11,10,12,13,19,18,1a,1b,1e,1f,1c,1d,05,04,06,07,01,00,02,03,09,08,0a,0b,0e,0f,0c,0d   p2a0   phi[15]
		05,04,06,07,01,00,02,03,09,08,0a,0b,0e,0f,0c,0d,19,18,1a,1b,1e,1f,1c,1d,11,10,12,13,16,17,14,15   p280   phi[14]
		19,18,1a,1b,1e,1f,1c,1d,11,10,12,13,16,17,14,15,09,08,0a,0b,0e,0f,0c,0d,01,00,02,03,06,07,04,05   p260   phi[13]
		09,08,0a,0b,0e,0f,0c,0d,01,00,02,03,06,07,04,05,11,10,12,13,16,17,14,15,1e,1f,1c,1d,1a,1b,18,19   p240   phi[12]
		11,10,12,13,16,17,14,15,1e,1f,1c,1d,1a,1b,18,19,01,00,02,03,06,07,04,05,0e,0f,0c,0d,0a,0b,08,09   p220   phi[11]
		01,00,02,03,06,07,04,05,0e,0f,0c,0d,0a,0b,08,09,1e,1f,1c,1d,1a,1b,18,19,16,17,14,15,12,13,10,11 + p200 = phi[10]
		1e,1f,1c,1d,1a,1b,18,19,16,17,14,15,12,13,10,11,0e,0f,0c,0d,0a,0b,08,09,06,07,04,05,02,03,00,01   p1e0   phi[0f]
		0e,0f,0c,0d,0a,0b,08,09,06,07,04,05,02,03,00,01,16,17,14,15,12,13,10,11,1a,1b,18,19,1c,1d,1f,1e   p1c0   phi[0e]
		16,17,14,15,12,13,10,11,1a,1b,18,19,1c,1d,1f,1e,06,07,04,05,02,03,00,01,0a,0b,08,09,0c,0d,0f,0e   p1a0   phi[0d]
		06,07,04,05,02,03,00,01,0a,0b,08,09,0c,0d,0f,0e,1a,1b,18,19,1c,1d,1f,1e,12,13,10,11,14,15,17,16   p180   phi[0c]
		1a,1b,18,19,1c,1d,1f,1e,12,13,10,11,14,15,17,16,0a,0b,08,09,0c,0d,0f,0e,02,03,00,01,04,05,07,06   p160   phi[0b]
		0a,0b,08,09,0c,0d,0f,0e,02,03,00,01,04,05,07,06,12,13,10,11,14,15,17,16,1c,1d,1f,1e,18,19,1b,1a   p140   phi[0a]
		12,13,10,11,14,15,17,16,1c,1d,1f,1e,18,19,1b,1a,02,03,00,01,04,05,07,06,0c,0d,0f,0e,08,09,0b,0a   p120   phi[09]
		02,03,00,01,04,05,07,06,0c,0d,0f,0e,08,09,0b,0a,1c,1d,1f,1e,18,19,1b,1a,14,15,17,16,10,11,13,12   p100   phi[08]
		1c,1d,1f,1e,18,19,1b,1a,14,15,17,16,10,11,13,12,0c,0d,0f,0e,08,09,0b,0a,04,05,07,06,00,01,03,02   p0e0   phi[07]
		0c,0d,0f,0e,08,09,0b,0a,04,05,07,06,00,01,03,02,14,15,17,16,10,11,13,12,18,19,1b,1a,1f,1e,1d,1c   p0c0   phi[06]
		14,15,17,16,10,11,13,12,18,19,1b,1a,1f,1e,1d,1c,04,05,07,06,00,01,03,02,08,09,0b,0a,0f,0e,0d,0c   p0a0   phi[05]
		04,05,07,06,00,01,03,02,08,09,0b,0a,0f,0e,0d,0c,18,19,1b,1a,1f,1e,1d,1c,10,11,13,12,17,16,15,14   p080   phi[04]
		18,19,1b,1a,1f,1e,1d,1c,10,11,13,12,17,16,15,14,08,09,0b,0a,0f,0e,0d,0c,00,01,03,02,07,06,05,04   p060   phi[03]
		08,09,0b,0a,0f,0e,0d,0c,00,01,03,02,07,06,05,04,10,11,13,12,17,16,15,14,1f,1e,1d,1c,1b,1a,19,18   p040   phi[02]
		10,11,13,12,17,16,15,14,1f,1e,1d,1c,1b,1a,19,18,00,01,03,02,07,06,05,04,0f,0e,0d,0c,0b,0a,09,08   p020   phi[01]
	=
		phi[00]   [01327654fedcba98] + p00, [fedcba9876543210] + p10   phi[00]   [0] + p00, [1] + p10
		phi[1e]   [fedcba9876543210] + p00, [76543210ba98dcef] + p10   phi[1e]   [1] + p00, [2] + p10
		phi[1d]   [76543210ba98dcef] + p10, [76543210ba98dcef] + p00   phi[1d]   [2] + p10, [2] + p00
		phi[1c]   [76543210ba98dcef] + p00, [ba98dcef32105467] + p10   phi[1c]   [2] + p00, [3] + p10
		phi[1b]   [ba98dcef32105467] + p10, [ba98dcef32105467] + p00   phi[1b]   [3] + p10, [3] + p00
		phi[1a]   [ba98dcef32105467] + p00, [32105467dcef98ab] + p10   phi[1a]   [3] + p00, [4] + p10
		phi[19]   [32105467dcef98ab] + p10, [32105467dcef98ab] + p00   phi[19]   [4] + p10, [4] + p00
		phi[18]   [32105467dcef98ab] + p00, [dcef98ab54671023] + p10   phi[18]   [4] + p00, [5] + p10
		phi[17]   [dcef98ab54671023] + p10, [dcef98ab54671023] + p00   phi[17]   [5] + p10, [5] + p00
		phi[16]   [dcef98ab54671023] + p00, [5467102398abefcd] + p10   phi[16]   [5] + p00, [6] + p10
		phi[15]   [5467102398abefcd] + p10, [5467102398abefcd] + p00   phi[15]   [6] + p10, [6] + p00
		phi[14]   [5467102398abefcd] + p00, [98abefcd10236745] + p10   phi[14]   [6] + p00, [7] + p10
		phi[13]   [98abefcd10236745] + p10, [98abefcd10236745] + p00   phi[13]   [7] + p10, [7] + p00
		phi[12]   [98abefcd10236745] + p00, [10236745efcdab89] + p10   phi[12]   [7] + p00, [8] + p10
		phi[11]   [10236745efcdab89] + p10, [10236745efcdab89] + p00   phi[11]   [8] + p10, [8] + p00
		phi[10] + [10236745efcdab89] + p00, [efcdab8967452301] + p10 = phi[10] + [8] + p00, [9] + p10
		phi[0f]   [efcdab8967452301] + p10, [efcdab8967452301] + p00   phi[0f]   [9] + p10, [9] + p00
		phi[0e]   [efcdab8967452301] + p00, [67452301ab89cdfe] + p10   phi[0e]   [9] + p00, [A] + p10
		phi[0d]   [67452301ab89cdfe] + p10, [67452301ab89cdfe] + p00   phi[0d]   [A] + p10, [A] + p00
		phi[0c]   [67452301ab89cdfe] + p00, [ab89cdfe23014576] + p10   phi[0c]   [A] + p00, [B] + p10
		phi[0b]   [ab89cdfe23014576] + p10, [ab89cdfe23014576] + p00   phi[0b]   [B] + p10, [B] + p00
		phi[0a]   [ab89cdfe23014576] + p00, [23014576cdfe89ba] + p10   phi[0a]   [B] + p00, [C] + p10
		phi[09]   [23014576cdfe89ba] + p10, [23014576cdfe89ba] + p00   phi[09]   [C] + p10, [C] + p00
		phi[08]   [23014576cdfe89ba] + p00, [cdfe89ba45760132] + p10   phi[08]   [C] + p00, [D] + p10
		phi[07]   [cdfe89ba45760132] + p10, [cdfe89ba45760132] + p00   phi[07]   [D] + p10, [D] + p00
		phi[06]   [cdfe89ba45760132] + p00, [4576013289bafedc] + p10   phi[06]   [D] + p00, [E] + p10
		phi[05]   [4576013289bafedc] + p10, [4576013289bafedc] + p00   phi[05]   [E] + p10, [E] + p00
		phi[04]   [4576013289bafedc] + p00, [89bafedc01327654] + p10   phi[04]   [E] + p00, [F] + p10
		phi[03]   [89bafedc01327654] + p10, [89bafedc01327654] + p00   phi[03]   [F] + p10, [F] + p00
		phi[02]   [89bafedc01327654] + p00, [01327654fedcba98] + p10   phi[02]   [F] + p00, [0] + p10
		phi[01]   [01327654fedcba98] + p10, [01327654fedcba98] + p00   phi[01]   [0] + p10, [0] + p00

		The phi and p00/p10-offset patterns match that of the [31 x radix-32] phase of DIF, but the [0-F] perm16 patterns differ slightly:
		DIF:		DIT:
		[0],[0]		[0],[1] [leftmost] indices match for DIF+DIT = (row+1)/2 [compare to []-pattern index]
		[1],[1]		[1],[2]	<*** odd  row (and 0): 1st []-idx = (row+1)/2, 2nd []-idx = 1+(row+1)/2
		[2],[1]		[2],[2]	<*** even row: []-idx = 1+(row+1)/2 for both []
		[2],[2]		[2],[3]
		[3],[2]		[3],[3]
		[3],[3]		[3],[4]
		[4],[3]		[4],[4]
		[4],[4]		[4],[5]
		[5],[4]		[5],[5]
		[5],[5]		[5],[6]
		[6],[5]		[6],[6]
		[6],[6]		[6],[7]
		[7],[6]		[7],[7]
		[7],[7]		[7],[8]
		[8],[7]		[8],[8]
		[8],[8]		[8],[9]
		[9],[8]		[9],[9]
		[9],[9]		[9],[A]
		[A],[9]		[A],[A]
		[A],[A]		[A],[B]
		[B],[A]		[B],[B]
		[B],[B]		[B],[C]
		[C],[B]		[C],[C]
		[C],[C]		[C],[D]
		[D],[C]		[D],[D]
		[D],[D]		[D],[E]
		[E],[D]		[E],[E]
		[E],[E]		[E],[F]
		[F],[E]		[F],[F]
		[F],[F]		[F],[0]
		[0],[F]		[0],[0]	<*** Last row: need to do []-indexing (mod 16), thus 1st []-idx = 16 (mod 16) = 0
	*/
	/*...gather the needed data (992 64-bit complex, i.e. 1984 64-bit reals) and do 31 radix-32 transforms...*/
		/*
		The jj-offsets here are raw real-array offsets, so need to be doubled relative to complex-indexing;
		The iptr-offset is an index into the length-992 q-array, whose elements implicitly contain the needed doubling.
		*/
		for(int k = 0; k < 32; ++k) {
			jj[k] = ((k<<5)-k)<<1;	/* (k*62) = (k*31)<<1 */
		}
	#if USE_COMPACT_OBJ_CODE
		iptr = o;	// Name offs-array here 'o' rather than 'q' so can use same length-32 array for DIF and DIT index-offsets
	#else
		iptr = q;
	#endif
		for(l = 0; l < 31; ++l) {
		#if USE_COMPACT_OBJ_CODE
			idx = (31-l) & (-(l>0));	// Base-offset Index into phi[] for the current row = 00,3c,3a,...,4,2
										// (Since our phi array here is in mults of 32, phi = 0x3c0 ==> idx = 3c/2 = 1e)
			is_odd = (l&1) + (l==0); is_even = (~is_odd)&1;	// 0-row phi-pair relative offsets behave like odd-index row
			// Compute index in the perm16 array needed for current row:
			pidx = 1+((l+1)>>1);	//<*** odd  row (and 0): 1st []-idx = (row+1)/2, 2nd []-idx = 1+(row+1)/2
									//<*** even row: []-idx = 1+(row+1)/2 for both []
			// First 16 offsets:
			jp = phi[idx] + plo[is_even<<4];  // = phi[idx]+p10 for even-idx rows, = phi[idx] for odd]
			i64 = dit_perm16[(pidx-is_odd)&0xf];
			// And lastly add in the low parts - p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];		o[0x0] = jp + k0;
			k1 = plo[(i64 >> 56)&0xf];		o[0x1] = jp + k1;
			k2 = plo[(i64 >> 52)&0xf];		o[0x2] = jp + k2;
			k3 = plo[(i64 >> 48)&0xf];		o[0x3] = jp + k3;
			k4 = plo[(i64 >> 44)&0xf];		o[0x4] = jp + k4;
			k5 = plo[(i64 >> 40)&0xf];		o[0x5] = jp + k5;
			k6 = plo[(i64 >> 36)&0xf];		o[0x6] = jp + k6;
			k7 = plo[(i64 >> 32)&0xf];		o[0x7] = jp + k7;
			k8 = plo[(i64 >> 28)&0xf];		o[0x8] = jp + k8;
			k9 = plo[(i64 >> 24)&0xf];		o[0x9] = jp + k9;
			ka = plo[(i64 >> 20)&0xf];		o[0xa] = jp + ka;
			kb = plo[(i64 >> 16)&0xf];		o[0xb] = jp + kb;
			kc = plo[(i64 >> 12)&0xf];		o[0xc] = jp + kc;
			kd = plo[(i64 >>  8)&0xf];		o[0xd] = jp + kd;
			ke = plo[(i64 >>  4)&0xf];		o[0xe] = jp + ke;
			kf = plo[(i64      )&0xf];		o[0xf] = jp + kf;
			// Second 16 offsets:
			jp = phi[idx] + plo[is_odd<<4];  // = phi[idx] for even-idx rows, = phi[idx]+p10 for odd]
			i64 = dit_perm16[pidx&0xf];	// 2nd []-idx = (row+1)/2 for both even and odd rows
			// And lastly add in the low parts - p-offset indices encoded in little-endian hex-char fashion:
			k0 = plo[(i64 >> 60)&0xf];		o[0x10] = jp + k0;
			k1 = plo[(i64 >> 56)&0xf];		o[0x11] = jp + k1;
			k2 = plo[(i64 >> 52)&0xf];		o[0x12] = jp + k2;
			k3 = plo[(i64 >> 48)&0xf];		o[0x13] = jp + k3;
			k4 = plo[(i64 >> 44)&0xf];		o[0x14] = jp + k4;
			k5 = plo[(i64 >> 40)&0xf];		o[0x15] = jp + k5;
			k6 = plo[(i64 >> 36)&0xf];		o[0x16] = jp + k6;
			k7 = plo[(i64 >> 32)&0xf];		o[0x17] = jp + k7;
			k8 = plo[(i64 >> 28)&0xf];		o[0x18] = jp + k8;
			k9 = plo[(i64 >> 24)&0xf];		o[0x19] = jp + k9;
			ka = plo[(i64 >> 20)&0xf];		o[0x1a] = jp + ka;
			kb = plo[(i64 >> 16)&0xf];		o[0x1b] = jp + kb;
			kc = plo[(i64 >> 12)&0xf];		o[0x1c] = jp + kc;
			kd = plo[(i64 >>  8)&0xf];		o[0x1d] = jp + kd;
			ke = plo[(i64 >>  4)&0xf];		o[0x1e] = jp + ke;
			kf = plo[(i64      )&0xf];		o[0x1f] = jp + kf;
			RADIX_32_DIT(
				a+j1,iptr,RE_IM_STRIDE,
				(double *)t,jj,1
			);
		#else
			RADIX_32_DIT(
				a+j1,iptr,RE_IM_STRIDE,	// *** Apr 2014: Added data strdes needed by latest version of the radix-32 DFT macros
				(double *)t,jj,1
			);
			iptr += 32;
		#endif
			for(int k = 0; k < 32; ++k)
			{
				jj[k] += 2;
			}
		}
		/*...and now do 32 radix-31 transforms. The required output permutation is as follows:

			000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,
			01f,03f,05f,07f,09f,0bf,0df,0ff,11f,13f,15f,17f,19f,1bf,1df,1ff,21f,23f,25f,27f,29f,2bf,2df,2ff,31f,33f,35f,37f,39f,3bf,3df,
			03e,05e,07e,09e,0be,0de,0fe,11e,13e,15e,17e,19e,1be,1de,1fe,21e,23e,25e,27e,29e,2be,2de,2fe,31e,33e,35e,37e,39e,3be,3de,01e,
			05d,07d,09d,0bd,0dd,0fd,11d,13d,15d,17d,19d,1bd,1dd,1fd,21d,23d,25d,27d,29d,2bd,2dd,2fd,31d,33d,35d,37d,39d,3bd,3dd,01d,03d,
			07c,09c,0bc,0dc,0fc,11c,13c,15c,17c,19c,1bc,1dc,1fc,21c,23c,25c,27c,29c,2bc,2dc,2fc,31c,33c,35c,37c,39c,3bc,3dc,01c,03c,05c,
			09b,0bb,0db,0fb,11b,13b,15b,17b,19b,1bb,1db,1fb,21b,23b,25b,27b,29b,2bb,2db,2fb,31b,33b,35b,37b,39b,3bb,3db,01b,03b,05b,07b,
			0ba,0da,0fa,11a,13a,15a,17a,19a,1ba,1da,1fa,21a,23a,25a,27a,29a,2ba,2da,2fa,31a,33a,35a,37a,39a,3ba,3da,01a,03a,05a,07a,09a,
			0d9,0f9,119,139,159,179,199,1b9,1d9,1f9,219,239,259,279,299,2b9,2d9,2f9,319,339,359,379,399,3b9,3d9,019,039,059,079,099,0b9,
			0f8,118,138,158,178,198,1b8,1d8,1f8,218,238,258,278,298,2b8,2d8,2f8,318,338,358,378,398,3b8,3d8,018,038,058,078,098,0b8,0d8,
			117,137,157,177,197,1b7,1d7,1f7,217,237,257,277,297,2b7,2d7,2f7,317,337,357,377,397,3b7,3d7,017,037,057,077,097,0b7,0d7,0f7,
			136,156,176,196,1b6,1d6,1f6,216,236,256,276,296,2b6,2d6,2f6,316,336,356,376,396,3b6,3d6,016,036,056,076,096,0b6,0d6,0f6,116,
			155,175,195,1b5,1d5,1f5,215,235,255,275,295,2b5,2d5,2f5,315,335,355,375,395,3b5,3d5,015,035,055,075,095,0b5,0d5,0f5,115,135,
			174,194,1b4,1d4,1f4,214,234,254,274,294,2b4,2d4,2f4,314,334,354,374,394,3b4,3d4,014,034,054,074,094,0b4,0d4,0f4,114,134,154,
			193,1b3,1d3,1f3,213,233,253,273,293,2b3,2d3,2f3,313,333,353,373,393,3b3,3d3,013,033,053,073,093,0b3,0d3,0f3,113,133,153,173,
			1b2,1d2,1f2,212,232,252,272,292,2b2,2d2,2f2,312,332,352,372,392,3b2,3d2,012,032,052,072,092,0b2,0d2,0f2,112,132,152,172,192,
			1d1,1f1,211,231,251,271,291,2b1,2d1,2f1,311,331,351,371,391,3b1,3d1,011,031,051,071,091,0b1,0d1,0f1,111,131,151,171,191,1b1,
			1f0,210,230,250,270,290,2b0,2d0,2f0,310,330,350,370,390,3b0,3d0,010,030,050,070,090,0b0,0d0,0f0,110,130,150,170,190,1b0,1d0,
			20f,22f,24f,26f,28f,2af,2cf,2ef,30f,32f,34f,36f,38f,3af,3cf,00f,02f,04f,06f,08f,0af,0cf,0ef,10f,12f,14f,16f,18f,1af,1cf,1ef,
			22e,24e,26e,28e,2ae,2ce,2ee,30e,32e,34e,36e,38e,3ae,3ce,00e,02e,04e,06e,08e,0ae,0ce,0ee,10e,12e,14e,16e,18e,1ae,1ce,1ee,20e,
			24d,26d,28d,2ad,2cd,2ed,30d,32d,34d,36d,38d,3ad,3cd,00d,02d,04d,06d,08d,0ad,0cd,0ed,10d,12d,14d,16d,18d,1ad,1cd,1ed,20d,22d,
			26c,28c,2ac,2cc,2ec,30c,32c,34c,36c,38c,3ac,3cc,00c,02c,04c,06c,08c,0ac,0cc,0ec,10c,12c,14c,16c,18c,1ac,1cc,1ec,20c,22c,24c,
			28b,2ab,2cb,2eb,30b,32b,34b,36b,38b,3ab,3cb,00b,02b,04b,06b,08b,0ab,0cb,0eb,10b,12b,14b,16b,18b,1ab,1cb,1eb,20b,22b,24b,26b,
			2aa,2ca,2ea,30a,32a,34a,36a,38a,3aa,3ca,00a,02a,04a,06a,08a,0aa,0ca,0ea,10a,12a,14a,16a,18a,1aa,1ca,1ea,20a,22a,24a,26a,28a,
			2c9,2e9,309,329,349,369,389,3a9,3c9,009,029,049,069,089,0a9,0c9,0e9,109,129,149,169,189,1a9,1c9,1e9,209,229,249,269,289,2a9,
			2e8,308,328,348,368,388,3a8,3c8,008,028,048,068,088,0a8,0c8,0e8,108,128,148,168,188,1a8,1c8,1e8,208,228,248,268,288,2a8,2c8,
			307,327,347,367,387,3a7,3c7,007,027,047,067,087,0a7,0c7,0e7,107,127,147,167,187,1a7,1c7,1e7,207,227,247,267,287,2a7,2c7,2e7,
			326,346,366,386,3a6,3c6,006,026,046,066,086,0a6,0c6,0e6,106,126,146,166,186,1a6,1c6,1e6,206,226,246,266,286,2a6,2c6,2e6,306,
			345,365,385,3a5,3c5,005,025,045,065,085,0a5,0c5,0e5,105,125,145,165,185,1a5,1c5,1e5,205,225,245,265,285,2a5,2c5,2e5,305,325,
			364,384,3a4,3c4,004,024,044,064,084,0a4,0c4,0e4,104,124,144,164,184,1a4,1c4,1e4,204,224,244,264,284,2a4,2c4,2e4,304,324,344,
			383,3a3,3c3,003,023,043,063,083,0a3,0c3,0e3,103,123,143,163,183,1a3,1c3,1e3,203,223,243,263,283,2a3,2c3,2e3,303,323,343,363,
			3a2,3c2,002,022,042,062,082,0a2,0c2,0e2,102,122,142,162,182,1a2,1c2,1e2,202,222,242,262,282,2a2,2c2,2e2,302,322,342,362,382,
			3c1,001,021,041,061,081,0a1,0c1,0e1,101,121,141,161,181,1a1,1c1,1e1,201,221,241,261,281,2a1,2c1,2e1,301,321,341,361,381,3a1,
		=
			000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,   p00
			000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,   p1f
			020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,   p1e
			040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,   p1d
			060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,   p1c
			080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,   p1b
			0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,   p1a
			0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,   p19
			0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,   p18
			100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,   p17
			120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,   p16
			140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,   p15
			160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,   p14
			180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,   p13
			1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,   p12
			1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,   p11
			1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,   p10
			200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,   p0f
			220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,   p0e
			240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,   p0d
			260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,   p0c
			280,2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,   p0b
			2a0,2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,   p0a
			2c0,2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,   p09
			2e0,300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,   p08
			300,320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,   p07
			320,340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,   p06
			340,360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,   p05
			360,380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,   p04
			380,3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,   p03
			3a0,3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,   p02
			3c0,000,020,040,060,080,0a0,0c0,0e0,100,120,140,160,180,1a0,1c0,1e0,200,220,240,260,280,2a0,2c0,2e0,300,320,340,360,380,3a0,   p01

		The plo indices in the  rcol, 00,1f,1e,1d,1c,1b,1a,19,18,17,16,15,14,13,12,11,10,0f,0e,0d,0c,0b,0a,09,08,07,06,05,04,03,02,01,
		are computable based on the loop index as

			(32-l) & (-(i>0));

		Each row's high parts are just cshifts of the basic (row 0) perm31, which is strictly ascending, i.e.
		identical to the phi[] array; thus we create a cperm array simply by augmenting phi[] with a copy of its 31 elements.
		The per-row leftward shift is just 0,0,1,2,...,30, which we can easily compute on the fly via

			lshift = (i-1) & (-(i>0));
		*/
		tptr = t;
		for(int k = 0; k < 32; ++k)
		{
		#if USE_COMPACT_OBJ_CODE
			mask = (-(k>0));
			lshift = (k-1) & mask;
			iptr = phi + lshift;
			RADIX_31_DIT(
				(double *)tptr,
				a+j1+plo[(32-k) & mask], iptr
			);
		#else
			iptr = o + ((k<<5)-k);	// o + k*31
			RADIX_31_DIT(
				(double *)tptr,
				a+j1,iptr
			);
		#endif
			tptr += 31;
		}
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy992_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr,*addi;
		struct complex *tptr;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p1,p2,p3;
		int poff[RADIX>>2];
		// Need storage for circular-shifts perms of a basic 31-vector, with shift count in [0,31] that means 2*31 elts:
		int dif_p20_cperms[62], plo[32],phi[62],jj[32], *iptr;
		int idx,pidx,mask,lshift, is_even,is_odd, k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,ka,kb,kc,kd,ke,kf, o[32];	// o[] stores o-address offsets for current radix-32 DFT in the 31x-loop
		uint64 i64;
	// DIF:
		// Low parts [p0-f] of output-index perms:
		const uint64 dif_perm16[16] = {
			0x0123456789abcdefull, 0xfecd89ab01234567ull, 0x76450123fecd89abull, 0xba89fecd76450123ull,
			0x32017645ba89fecdull, 0xdcfeba8932017645ull, 0x54763201dcfeba89ull, 0x98badcfe54763201ull,
			0x1032547698badcfeull, 0xefdc98ba10325476ull, 0x67541032efdc98baull, 0xab98efdc67541032ull,
			0x23106754ab98efdcull, 0xcdefab9823106754ull, 0x45672310cdefab98ull, 0x89abcdef45672310ull
		};
	// DIT:
		// Low parts [p0-f] of output-index perms:
		const uint64 dit_perm16[16] = {
			0x01327654fedcba98ull,0xfedcba9876543210ull,0x76543210ba98dcefull,0xba98dcef32105467ull,
			0x32105467dcef98abull,0xdcef98ab54671023ull,0x5467102398abefcdull,0x98abefcd10236745ull,
			0x10236745efcdab89ull,0xefcdab8967452301ull,0x67452301ab89cdfeull,0xab89cdfe23014576ull,
			0x23014576cdfe89baull,0xcdfe89ba45760132ull,0x4576013289bafedcull,0x89bafedc01327654ull
		};

		int j,j1,j2,jt,jp,l;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
		double rt,it, wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff

		double *base, *baseinv, *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr, bs,bsinv, wts_idx_incr;
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

	// int data:
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

	// double data:
		double maxerr = thread_arg->maxerr;
		double scale = thread_arg->scale;
		double prp_mult = thread_arg->prp_mult;

	// pointer data:
		double *a = thread_arg->arrdat;
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
		int *si = thread_arg->si;
		struct complex *rn0 = thread_arg->rn0;
		struct complex *rn1 = thread_arg->rn1;

		/*   constant index offsets for array load/stores are here.	*/
		// Legacy code needs 3 lowest nonzero fixed-index p[] terms:
		p1 = 1*NDIVR;	 p1 += ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = 2*NDIVR;	 p2 += ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = 3*NDIVR;	 p3 += ( (p3 >> DAT_BITS) << PAD_BITS );

		for(l = 0; l < (RADIX>>2); l++) {
			poff[l] = (l<<2)*NDIVR;	// Corr. to every 4th plo[] term
			poff[l] += ( (poff[l] >> DAT_BITS) << PAD_BITS );
		}

		// Init plo,phi arrays and circular-shifts perms of a basic 31-vector, with shift count in [0,31] that means 2 copies of a basic 31-elt set:
		for(l = 0; l < 31; l++)
		{
			plo[l] = l*NDIVR; phi[l] = plo[l]<<5;
			plo[l] += ( (plo[l] >> DAT_BITS) << PAD_BITS );
			phi[l] += ( (phi[l] >> DAT_BITS) << PAD_BITS );
			phi[31+l] = phi[l];	// Needed for the DIT cperm, which needs just a doubled version of phi[0-30]
			dif_p20_cperms[31-l] = dif_p20_cperms[62-l] = phi[l];	// p3c0,3a0,...,020; init in reverse order since init of needed phi elts run forward in same loop
		}	plo[l] = l*NDIVR;	plo[l] += ( (plo[l] >> DAT_BITS) << PAD_BITS );	// Need 1 more plo term than phi
		dif_p20_cperms[0] = dif_p20_cperms[31] = 0;	// Handle the leading 0 and its copy separately

		// In scalar mode use these 2 ptrs to pass wts_idx_incr and the base/baseinv/etc array-ptrs:
		wts_idx_incr = *(int *)thread_arg->half_arr;
		base      = (double *)thread_arg->r00;
		baseinv   = base + 2;
		wt_arr    = base + 4;
		wtinv_arr = wt_arr    + ODD_RADIX;
		bs_arr    = wtinv_arr + ODD_RADIX;
		bsinv_arr = bs_arr    + ODD_RADIX;

		// Can't simply use thread-associated values of these *cycle index arrays here, since
		// thread values must be ***read-only*** so as to retain the proper first-init values
		// on each entry to this thread-task. Instead use the bjmodn data storage block - which
		// is otherwise unused in Fermat-Mod mode - for local storage of these cycle tables:
		int *icycle = bjmodn,ic;
		for(j = 0; j < ODD_RADIX; j++) {
			icycle[j] = thread_arg->icycle[j];
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
			// No_op in scalar case, since carry pattern matches that of thread data
		}

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix992_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
			// No_op in scalar case, since carry pattern matches that of thread data

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
