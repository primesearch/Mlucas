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

#ifdef USE_FGT61
	// Fwd decl:
	void mixed_carry(
		const int sw, const int bw, const int nm1, const int bits,
		const int j, const int col, const int co2, const int co3, const int n_minus_sil, const int n_minus_silp1, const int sinwt, const int sinwtm1,
		const double base[], const double baseinv[], const double wt1[], const double wtl, const double wtlp1, const double wtn, const double wtnm1,
		double*maxerr, double*fx, double*fy, uint64*ix, uint64*iy, double*cy, uint32 bjmodn, uint32 set, uint32 shift, uint32 l2_n2);
#endif

#define RADIX 16	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

#define EPS 1e-10

#ifndef USE_FGT61
	#define USE_SCALAR_DFT_MACRO	1
#else
	// v20: Make default - need outputs back in array to support shifted-residue carry injection:
	#define USE_SCALAR_DFT_MACRO	1
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

  // For Mersenne-mod need (16 [SSE2] or 64 [AVX]) + (24 [SSE2] or 68 [AVX])[LOACC] added slots for half_arr lookup tables.
  // For Fermat-mod we use RADIX*4 = 64 [note there is no Fermat-mod LOACC option for this power-of-2 DFT] more
  // slots in AVX mode for the compact negacyclic-roots chained-multiply scheme. Add larger of the 2 numbers -
  // 128 for AVX, 40 for SSE2 - to (half_arr_offset16 + RADIX) to get AVX value of radix16_creals_in_local_store:
  #ifdef USE_AVX512
	const int half_arr_offset16 = 42;	// + RADIX = 58; Used for thread local-storage-integrity checking
	const int radix16_creals_in_local_store = 196;	// May not need this many slots, but just use AVX-alloc for now
  #elif defined(USE_AVX)
	const int half_arr_offset16 = 46;	// + RADIX = 62; Used for thread local-storage-integrity checking
	const int radix16_creals_in_local_store = 196;	// (half_arr_offset16 + RADIX) + 128 and round up to nearest multiple of 4
  #else
	const int half_arr_offset16 = 54;	// + RADIX = 70; Used for thread local-storage-integrity checking
	const int radix16_creals_in_local_store = 112;	// (half_arr_offset16 + RADIX) + 40 and round up to nearest multiple of 4
  #endif

	#include "sse2_macro_gcc64.h"

	#ifdef COMPILER_TYPE_MSVC
		/*  */
	#else	/* GCC-style inline ASM: */
		#include "radix16_ditN_cy_dif1_asm.h"
	#endif

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
		int _pad_;

	// double data:
		double maxerr;
		double scale;
		double prp_mult;

	// pointer data:
		double *arrdat;			/* Main data array */
	#ifdef USE_FGT61
		uint64 *brrdat;
	#endif
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

		int bjmodn0;
		int bjmodn1;
		int bjmodn2;
		int bjmodn3;
		int bjmodn4;
		int bjmodn5;
		int bjmodn6;
		int bjmodn7;
		int bjmodn8;
		int bjmodn9;
		int bjmodnA;
		int bjmodnB;
		int bjmodnC;
		int bjmodnD;
		int bjmodnE;
		int bjmodnF;
		/* carries: */
		double cy_r0;
		double cy_r1;
		double cy_r2;
		double cy_r3;
		double cy_r4;
		double cy_r5;
		double cy_r6;
		double cy_r7;
		double cy_r8;
		double cy_r9;
		double cy_rA;
		double cy_rB;
		double cy_rC;
		double cy_rD;
		double cy_rE;
		double cy_rF;

		double cy_i0;
		double cy_i1;
		double cy_i2;
		double cy_i3;
		double cy_i4;
		double cy_i5;
		double cy_i6;
		double cy_i7;
		double cy_i8;
		double cy_i9;
		double cy_iA;
		double cy_iB;
		double cy_iC;
		double cy_iD;
		double cy_iE;
		double cy_iF;
	};

#endif

/***************/

#ifdef USE_FGT61
int radix16_ditN_cy_dif1		(double a[], uint64 b[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
#else
int radix16_ditN_cy_dif1		(double a[],             int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
#endif
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-16 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-16 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix16_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix16_ditN_cy_dif1";
	static int NDIVR;
	const int pfetch_dist = PFETCH_DIST;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
  #ifdef USE_AVX512
	const int jhi_wrap_mers = 15;
	const int jhi_wrap_ferm = 15;
  #else
	const int jhi_wrap_mers =  7;
	const int jhi_wrap_ferm = 15;	// For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
  #endif
	int i,j,j1,j2,jstart,jhi,full_pass,k,khi,l,outer,nbytes;
	int col,co2,co3;
  #ifdef USE_AVX512
//	double t0,t1,t2,t3;	This routine already has double t1-32 def'd below, so just prepend t0 to those
	static struct uint32x8 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #elif defined(USE_AVX)
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
  #endif
	static double wts_mult[2], inv_mult[2];	// Const wts-multiplier and 2*(its multiplicative inverse)
	double rt,it, wt_re,wt_im, wi_re,wi_im;	// Fermat-mod weights stuff, used in both scalar and AVX mode
	static uint32 bjmodnini, nsave = 0;
	// Jun 2018: Add support for residue shift. (Only LL-test needs intervention at carry-loop level).
	int target_idx = -1, target_set = 0,tidx_mod_stride;
	double target_cy = 0;
	static double ndivr_inv;
	uint64 itmp64;
	static uint64 psave = 0;
	static uint32 bw,sw,nm1,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15;
	static int poff[RADIX>>2],p0123[4];	// Store [RADIX/4] mults of p04 offset for loop control
	const double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;
	static double radix_inv,n2inv,scale;	/* Need scale to be static since in the MSVC-only pure-ASM vesrion of the carry step save the address in a static pointer below */
	double *addr;
	double dtmp, maxerr = 0.0;
	double t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;
	double temp,frac
		,a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r,a1p7r,a1p8r,a1p9r,a1pAr,a1pBr,a1pCr,a1pDr,a1pEr,a1pFr
		,a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i,a1p7i,a1p8i,a1p9i,a1pAi,a1pBi,a1pCi,a1pDi,a1pEi,a1pFi;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	static int n_div_nwt;

#ifdef USE_SSE2

	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_nm1;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else
	static int idx_offset, idx_incr;
	double *add0,*add1,*add2,*add3;	/* Addresses into array sections */
  #endif

	const double crnd = 3.0*0x4000000*0x2000000;
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	static int *bjmodn0,*bjmodn1,*bjmodn2,*bjmodn3,*bjmodn4,*bjmodn5,*bjmodn6,*bjmodn7,*bjmodn8,*bjmodn9,*bjmodnA,*bjmodnB,*bjmodnC,*bjmodnD,*bjmodnE,*bjmodnF;
	static double *wt0_ptr, *wt1_ptr, *scale_ptr = &scale;
	static vec_dbl *cc0, *ss0, *isrt2, *max_err, *sse2_rnd, *half_arr
		,*r00,*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09,*r0A,*r0B,*r0C,*r0D,*r0E,*r0F
		,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r1A,*r1B,*r1C,*r1D,*r1E,*r1F;
  #ifdef USE_AVX512
	static vec_dbl *cy_r0,*cy_r8, *cy_i0,*cy_i8;
  #elif defined(USE_AVX)
	static vec_dbl *cy_r0,*cy_r4,*cy_r8,*cy_rC, *cy_i0,*cy_i4,*cy_i8,*cy_iC;
  #else	// SSE2:
	static vec_dbl *cy_r0,*cy_r2,*cy_r4,*cy_r6,*cy_r8,*cy_rA,*cy_rC,*cy_rE, *cy_i0,*cy_i2,*cy_i4,*cy_i6,*cy_i8,*cy_iA,*cy_iC,*cy_iE;
  #endif

  #ifdef USE_AVX	// AVX and above:
	static vec_dbl *base_negacyclic_root;
	int k1,k2;
  #endif

	vec_dbl *tmp,*tm2;	// Non-static utility ptrs

#endif

#ifdef MULTITHREAD

	static struct cy_thread_data_t *tdat = 0x0;
	// Threadpool-based dispatch stuff:
	static int main_work_units = 0, pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)cy16_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	int jt,jp,k1,k2,m,m2,ntmp;
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
  #if PFETCH
	double *addr, *addp;
  #endif
	int bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodnA,bjmodnB,bjmodnC,bjmodnD,bjmodnE,bjmodnF;
	double cy_r0,cy_r1,cy_r2,cy_r3,cy_r4,cy_r5,cy_r6,cy_r7,cy_r8,cy_r9,cy_rA,cy_rB,cy_rC,cy_rD,cy_rE,cy_rF,
			cy_i0,cy_i1,cy_i2,cy_i3,cy_i4,cy_i5,cy_i6,cy_i7,cy_i8,cy_i9,cy_iA,cy_iB,cy_iC,cy_iD,cy_iE,cy_iF;
  #ifdef USE_FGT61
	const uint64 q  = 0x1FFFFFFFFFFFFFFFull, q2=q+q, q3=q2+q, q4=q2+q2, q5=q4+q, q8=q4+q4;	// q = 2^61 - 1, and needed small multiples
	// primitive 16th root of unity, scaled by *8:
	const uint64 cm = 1693317751237720973ull<<3, sm = 2283815672160731785ull<<3;
	uint64 rm,im;
	// varname m2 already used for float-carry macros, so replace mod-varname m2 with m$ in carry step
	uint64 m1,m$,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32
		,b1p0r,b1p1r,b1p2r,b1p3r,b1p4r,b1p5r,b1p6r,b1p7r,b1p8r,b1p9r,b1pAr,b1pBr,b1pCr,b1pDr,b1pEr,b1pFr
		,b1p0i,b1p1i,b1p2i,b1p3i,b1p4i,b1p5i,b1p6i,b1p7i,b1p8i,b1p9i,b1pAi,b1pBi,b1pCi,b1pDi,b1pEi,b1pFi;
	static uint8 mod_wt_exp[64];	// Pad the required 61 slots out to 64 for alignment reasons
	static uint32 swmod61,nmod61;
	uint32 bits,simodnmod61;
	static int l2_n2;	// n2 in modular form (i.e. base-2-logarithmically, that is in form of a shift count)
	const double inv61 = 1.0/61.0;
	const uint64 qhalf = 0x1000000000000000ull;	// (q+1)/2 = 2^60
  #endif

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 nt_save = 0xffffffff, CY_THREADS = 0,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_bjmodn0 = 0x0,*_bjmodn1 = 0x0,*_bjmodn2 = 0x0,*_bjmodn3 = 0x0,*_bjmodn4 = 0x0,*_bjmodn5 = 0x0,*_bjmodn6 = 0x0,*_bjmodn7 = 0x0,*_bjmodn8 = 0x0,*_bjmodn9 = 0x0,*_bjmodnA = 0x0,*_bjmodnB = 0x0,*_bjmodnC = 0x0,*_bjmodnD = 0x0,*_bjmodnE = 0x0,*_bjmodnF = 0x0;
	static int *_bjmodnini = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double
	*_cy_r0 = 0x0,*_cy_r1 = 0x0,*_cy_r2 = 0x0,*_cy_r3 = 0x0,*_cy_r4 = 0x0,*_cy_r5 = 0x0,*_cy_r6 = 0x0,*_cy_r7 = 0x0,*_cy_r8 = 0x0,*_cy_r9 = 0x0,*_cy_rA = 0x0,*_cy_rB = 0x0,*_cy_rC = 0x0,*_cy_rD = 0x0,*_cy_rE = 0x0,*_cy_rF = 0x0,
	*_cy_i0 = 0x0,*_cy_i1 = 0x0,*_cy_i2 = 0x0,*_cy_i3 = 0x0,*_cy_i4 = 0x0,*_cy_i5 = 0x0,*_cy_i6 = 0x0,*_cy_i7 = 0x0,*_cy_i8 = 0x0,*_cy_i9 = 0x0,*_cy_iA = 0x0,*_cy_iB = 0x0,*_cy_iC = 0x0,*_cy_iD = 0x0,*_cy_iE = 0x0,*_cy_iF = 0x0;

  #ifdef USE_IMCI512
	WARN(HERE, "radix16_ditN_cy_dif1: No k1om / IMCI-512 support; Skipping this leading radix.", "", 1); return(ERR_RADIX0_UNAVAILABLE);
  #endif

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
	NDIVR   = n/16;	ndivr_inv = (double)RADIX/n;
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

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		psave = p;	nsave = n;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)16));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/
	#ifdef USE_FGT61
		bits = p/n;		// #bits in a smallword
		swmod61 = sw % 61;	// Needed for the modular weights (see comments below).
		nmod61  = n  % 61;
		// Calculate the shift counts which substitute for the modular weight factors.
		/*
		Modular analog of WT(J) satisfies WM(J)^N == [2^(s*j mod N)] mod q = 2^[(s*j mod N) mod 61], since q = 2^61 - 1.
		Assume WM(J) = 2^A. Then WM(J)^N = 2^(A*N), so A satisfies the congruence A*N == (s*j mod N) mod 61, with A in [0,60]. These are then stored in a length-61 table,
		The appropriate element is accessed using the current value of (s*j mod N) mod 61, and sent (along with the
		digit to be multiplied by the weight factor) as a shift count to MUL_POW2_MODQ. Since WM is a power of 2, we need
		only store the shift count A, not the modular weight factor WM itself.
			Modular inverse weights are also easy: since WM*WMINV == 1 mod 2^61 - 1, if WM = 2^A, WMINV = 2^(61-A),
		i.e. to multiply by the inverse of 2^A, simply use 61-A as the shift count to be passed to MUL_POW2_MODQ. Sweet!
		*/
		for(simodnmod61 = 0; simodnmod61 < 61; simodnmod61++) {
			k = 0;		// K stores A*N mod 61.
			for(i = 0; i < 61; i++) {
				if(k == simodnmod61) {
					mod_wt_exp[simodnmod61] = i;
					break;
				}
				k += nmod61;
				if(k > 60) k -= 61;
			}
		}
		ASSERT(isPow2(N2), "N/2 not a power of 2!");
		l2_n2 = trailz32(N2);
// ******* For carry step, also need the 16 values of bimodnmod61 for i = j*(n/radix0), j = 0,...,15 ************
	#endif
		nm1   = n-1;

	  #define CARRY_8_WAY	// Make default due to higher accuracy

	  #ifdef USE_AVX	// AVX LOACC: can select between carry macros processing 4 and 8 independent carries
	   #ifdef CARRY_8_WAY
		i = 8;
	   #else
		i = 4;
	   #endif
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

				main_work_units = 0;
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
			tdat[ithread].iter = iter;
		// int data:
			tdat[ithread].tid = ithread;
			tdat[ithread].ndivr = NDIVR;

			tdat[ithread].sw  = sw;
			tdat[ithread].nwt = nwt;

		// pointer data:
			// Dec 2015: fast-GCD usage of this routine may involve multiple 'main' arrays
			// on successive calls, so set here at runtime rather than in init-only block:
		//	tdat[ithread].arrdat = a;			/* Main data array */
		#ifdef USE_FGT61
			tdat[ithread].brrdat = b;
		#endif
			tdat[ithread].wt0 = wt0;
			tdat[ithread].wt1 = wt1;
			tdat[ithread].wts_mult = wts_mult;
			tdat[ithread].inv_mult = inv_mult;
			tdat[ithread].si  = si;
			tdat[ithread].rn0 = rn0;
			tdat[ithread].rn1 = rn1;
		}
	#endif

	#ifdef USE_SSE2

		ASSERT(((intptr_t)wt0 & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(((intptr_t)wt1 & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use vector-double type size (16 bytes for SSE2, 32 for AVX) to alloc a block of local storage
		// consisting of 128 dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix16_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "ERROR: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(((intptr_t)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix16_creals_in_local_store);
		ASSERT(((intptr_t)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 32 16-byte slots of sc_arr for temporaries, next 3 for the nontrivial complex 16th roots,
	next 16 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	#ifdef USE_PTHREAD
		__r0 = sc_ptr;
	#endif
		r00 = sc_ptr + 0x00;		r01 = sc_ptr + 0x01;
		r02 = sc_ptr + 0x02;		r03 = sc_ptr + 0x03;
		r04 = sc_ptr + 0x04;		r05 = sc_ptr + 0x05;
		r06 = sc_ptr + 0x06;		r07 = sc_ptr + 0x07;
		r08 = sc_ptr + 0x08;		r09 = sc_ptr + 0x09;
		r0A = sc_ptr + 0x0a;		r0B = sc_ptr + 0x0b;
		r0C = sc_ptr + 0x0c;		r0D = sc_ptr + 0x0d;
		r0E = sc_ptr + 0x0e;		r0F = sc_ptr + 0x0f;
		r10 = sc_ptr + 0x10;		r11 = sc_ptr + 0x11;
		r12 = sc_ptr + 0x12;		r13 = sc_ptr + 0x13;
		r14 = sc_ptr + 0x14;		r15 = sc_ptr + 0x15;
		r16 = sc_ptr + 0x16;		r17 = sc_ptr + 0x17;
		r18 = sc_ptr + 0x18;		r19 = sc_ptr + 0x19;
		r1A = sc_ptr + 0x1a;		r1B = sc_ptr + 0x1b;
		r1C = sc_ptr + 0x1c;		r1D = sc_ptr + 0x1d;
		r1E = sc_ptr + 0x1e;		r1F = sc_ptr + 0x1f;
	  isrt2 = sc_ptr + 0x20;
		cc0 = sc_ptr + 0x21;
		ss0 = sc_ptr + 0x22;	// +34 vec_dbl slots
		tmp = ss0 + 0x2;	// +36; Only need incr = 1 but prefer an even offset
	#ifdef USE_AVX512
		 cy_r0 = tmp + 0x00;
		 cy_r8 = tmp + 0x01;
		 cy_i0 = tmp + 0x02;
		 cy_i8 = tmp + 0x03;
		max_err = tmp + 0x04;
		sse2_rnd= tmp + 0x05;
		half_arr= tmp + 0x06;	// 36 + 6 = 42 = half_arr_offset16
	#elif defined(USE_AVX)
		 cy_r0 = tmp + 0x00;
		 cy_r4 = tmp + 0x01;
		 cy_r8 = tmp + 0x02;
		 cy_rC = tmp + 0x03;
		 cy_i0 = tmp + 0x04;
		 cy_i4 = tmp + 0x05;
		 cy_i8 = tmp + 0x06;
		 cy_iC = tmp + 0x07;
		max_err = tmp + 0x08;
		sse2_rnd= tmp + 0x09;
		half_arr= tmp + 0x0a;	// 36 + 10 = 46 = half_arr_offset16; This table needs 20x16 bytes
	#else
		 cy_r0 = tmp + 0x00;		cy_r2 = tmp + 0x01;
		 cy_r4 = tmp + 0x02;		cy_r6 = tmp + 0x03;
		 cy_r8 = tmp + 0x04;		cy_rA = tmp + 0x05;
		 cy_rC = tmp + 0x06;		cy_rE = tmp + 0x07;
		 cy_i0 = tmp + 0x08;		cy_i2 = tmp + 0x09;
		 cy_i4 = tmp + 0x0a;		cy_i6 = tmp + 0x0b;
		 cy_i8 = tmp + 0x0c;		cy_iA = tmp + 0x0d;
		 cy_iC = tmp + 0x0e;		cy_iE = tmp + 0x0f;
		max_err = tmp + 0x10;
		sse2_rnd= tmp + 0x11;
		half_arr= tmp + 0x12;	// 36 + 18 = 54 = half_arr_offset16; This table needs 20x16 bytes
	#endif
		/* These remain fixed: */
		VEC_DBL_INIT(isrt2, ISRT2);
		VEC_DBL_INIT(cc0  , c	);
		VEC_DBL_INIT(ss0  , s	);
	  #ifdef USE_AVX512	// In AVX-512 mode, use VRNDSCALEPD for rounding and hijack this vector-data slot for the 4 base/baseinv-consts
		sse2_rnd->d0 = base[0]; sse2_rnd->d1 = baseinv[1]; sse2_rnd->d2 = wts_mult[1]; sse2_rnd->d3 = inv_mult[0];
	  #else
		VEC_DBL_INIT(sse2_rnd, crnd);
	  #endif

		// Propagate the above consts to the remaining threads:
		nbytes = (intptr_t)ss0 - (intptr_t)isrt2 + SZ_VD;	// #bytes in 1st of above block of consts
		tmp = isrt2;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

		nbytes = SZ_VD;	// #bytes in 2nd block of consts, which contains just sse2_rnd
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
		*/
		tmp = half_arr;

	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		/* In Fermat-mod mode, use first 2 128-bit slots for base and 1/base: */
		VEC_DBL_INIT(tmp, base   [0]);	++tmp;
		VEC_DBL_INIT(tmp, baseinv[0]);	++tmp;
		/* [+2] slot is for [scale,scale] */

		// Propagate the above consts to the remaining threads:
		nbytes = 2*SZ_VD;
		tmp = half_arr;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	#if 0
	/* Simple qfloat-based loop to crank out the roots needed by th AVX negacyclic-mults scheme below: */
		struct qfloat qc,qs,qx,qy,qt,qn,qmul;
		qt = qfdiv(QPIHALF, i64_to_q((int64)RADIX));	// Pi/2/RADIX
		qc = qfcos(qt);	qs = qfsin(qt);
		qx = QONE;		qy = QZRO;
		for(j = 1; j < RADIX; j++) {
			// Up-multiply the complex exponential:
			qn = qfmul(qx, qc); qt = qfmul(qy, qs); qmul = qfsub(qn, qt);	// Store qxnew in qmul for now.
			qn = qfmul(qx, qs); qt = qfmul(qy, qc); qy   = qfadd(qn, qt); qx = qmul;
			printf("j = %3u: cos[j*Pi/2] = %#16" PRIX64 ", sin[j*Pi/2] = %#16" PRIX64 "\n",j,qfdbl_as_uint64(qx),qfdbl_as_uint64(qy));
		}
		exit(0);
	#endif

	#ifdef USE_AVX

		base_negacyclic_root = half_arr + RADIX;
#if 0
===================
*** Ferm-mod: ***
start_ptr:			#slots		description
---------			------		-----------
half_arr			  radix		Low 2 slots used for base,binv, rest unused?? Mers-mod: ???
half_arr+  radix	2*radix		RADIX/4 copies of the Re/Im parts of the 4 base multipliers
half_arr+3*radix	4[avx],
					8[avx512]	32 scalar-double base multipliers

*** Mers-mod: ***
start_ptr:			#slots		description
---------			------		-----------
half_arr			radix		1.0/0.5 lookup-table [lut]
half_arr+  radix	radix		0.5/.25-lut
half_arr+2*radix	radix		base-lut
half_arr+3*radix	radix		binv-lut
half_arr+4*radix	radix		[LOACC-only] wts_mult-lut
half_arr+5*radix	radix		[LOACC-only] inv_mult-lut
===================
#endif
	  #ifdef USE_AVX512	// 8-way-double analog of AVX inits below:
		/*
		The pattern of the negacyclic-DWT-weights ("nDWTs") applied to the RADIX complex outputs of the final-DIT-pass is like so:
		The nDWTs multiplying each set of [RADIX] DIT DFT outputs are simply the product of a single complex "base multiplier" [rbase]
		(separately computed for each batch of [RADIX] DFT outputs), which multiplies the (0 - RADIX-1)st [4*RADIX]th roots of unity,
		i.e.
			 rbase * (j*I*Pi)/(2*RADIX), for j = 0, ..., RADIX-1 .

		As applied to the 8-way-SIMD data in AVX-512 mode, we will have 8 such base roots at a time in a pair of ZMM registers
		(real parts in one ZMM, imaginary parts in another. The above [4*RADIX]th roots will be precomputed and stored in 2*[RADIX/8]
		ZMM-sized local-data slots (RADIX/8 for the Re parts, RADIX/8 for the Im) like so, using the above j-index as a shorthand:

			Slot 0: Re[j =  0 + 0, 1, 2, ... , 7]
			Slot 1: Im[j =  0 + 0, 1, 2, ... , 7]

			Slot 2: Re[j =  8 + 0, 1, 2, ... , 7]
			Slot 3: Im[j =  8 + 0, 1, 2, ... , 7]

		Prior to performing the normalize-and-carry step on each set of RADIX AVX-complex ( = 8*RADIX double-complex) DIT DFT outputs,
		we compute the 8 base multipliers needed for that set of data:

			Re[rbase 0, 1, 2, ... , 7]
			Im[rbase 0, 1, 2, ... , 7]

		and then do a complex multiply of that octet of complex-double data with each of the above RADIX/8 precomputed SIMD-complex
		constants, storing the results in another set of local-mem slots and/or ZMM registers, as desired.
		*/
		tmp = base_negacyclic_root + 2*RADIX;	// First 32 = 2*RADIX slots reserved for RADIX/8 copies of the Re/Im parts of the 8 base multipliers
		tm2 = tmp + RADIX/4 - 1;	// tmp+3
		// Use tmp-pointer to init tmp+0,2; use tm2 to init tmp+3,1:
										tmp->d0 = 1.0;	(tmp+1)->d0 = 0.0;
		tmp64 = 0x3FEFD88DA3D12526ull;	tmp->d1 = tm2->d7 = *(double *)&tmp64;
		tmp64 = 0x3FEF6297CFF75CB0ull;	tmp->d2 = tm2->d6 = *(double *)&tmp64;
		tmp64 = 0x3FEE9F4156C62DDAull;	tmp->d3 = tm2->d5 = *(double *)&tmp64;
		tmp64 = 0x3FED906BCF328D46ull;	tmp->d4 = tm2->d4 = *(double *)&tmp64;
		tmp64 = 0x3FEC38B2F180BDB1ull;	tmp->d5 = tm2->d3 = *(double *)&tmp64;
		tmp64 = 0x3FEA9B66290EA1A3ull;	tmp->d6 = tm2->d2 = *(double *)&tmp64;
		tmp64 = 0x3FE8BC806B151741ull;	tmp->d7 = tm2->d1 = *(double *)&tmp64;	tmp += 2;	// tmp+2
		tmp64 = 0x3FE6A09E667F3BCDull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	tm2 -= 2;	// tmp+1
		tmp64 = 0x3FE44CF325091DD6ull;	tmp->d1 = tm2->d7 = *(double *)&tmp64;
		tmp64 = 0x3FE1C73B39AE68C8ull;	tmp->d2 = tm2->d6 = *(double *)&tmp64;
		tmp64 = 0x3FDE2B5D3806F63Bull;	tmp->d3 = tm2->d5 = *(double *)&tmp64;
		tmp64 = 0x3FD87DE2A6AEA963ull;	tmp->d4 = tm2->d4 = *(double *)&tmp64;
		tmp64 = 0x3FD294062ED59F06ull;	tmp->d5 = tm2->d3 = *(double *)&tmp64;
		tmp64 = 0x3FC8F8B83C69A60Bull;	tmp->d6 = tm2->d2 = *(double *)&tmp64;
		tmp64 = 0x3FB917A6BC29B42Cull;	tmp->d7 = tm2->d1 = *(double *)&tmp64;	tmp += 2;	// tmp+4 (unused)

		tmp = base_negacyclic_root + 2*RADIX;	// reset to point to start of above block

	  #else
		/*
		The pattern of the negacyclic-DWT-weights ("nDWTs") applied to the RADIX complex outputs of the final-DIT-pass is like so:
		The nDWTs multiplying each set of [RADIX] DIT DFT outputs are simply the product of a single complex-root "base multiplier" [rbase]
		(separately computed for each batch of DFT outputs), which "base root" multiplies the (0 - RADIX-1)st [4*RADIX]th roots of unity,
		i.e.
			 rbase * (j*I*Pi)/(2*RADIX), for j = 0, ..., RADIX-1 .

		See the radix28 version of this routine for additional details.
		*/
		tmp = base_negacyclic_root + RADIX*2;	// First 32 = 2*RADIX slots reserved for RADIX/4 copies of the Re/Im parts of the 4 base multipliers
		tm2 = tmp + RADIX/2 - 1;	// tmp+7
		// Use tmp-pointer to init tmp+0,2,4,6; use tm2 to init tmp+7,5,3,1:
										tmp->d0 = 1.0;	(tmp+1)->d0 = 0.0;
		tmp64 = 0x3FEFD88DA3D12526ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(01*I*Pi/32) = sin(31*I*Pi/32) */
		tmp64 = 0x3FEF6297CFF75CB0ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(02*I*Pi/32) = sin(30*I*Pi/32) */
		tmp64 = 0x3FEE9F4156C62DDAull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(03*I*Pi/32) = sin(29*I*Pi/32) */	tmp += 2;	// tmp+2
		tmp64 = 0x3FED906BCF328D46ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(04*I*Pi/32) = sin(28*I*Pi/32) */	tm2 -= 2;	// tmp+5
		tmp64 = 0x3FEC38B2F180BDB1ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(05*I*Pi/32) = sin(27*I*Pi/32) */
		tmp64 = 0x3FEA9B66290EA1A3ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(06*I*Pi/32) = sin(26*I*Pi/32) */
		tmp64 = 0x3FE8BC806B151741ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(07*I*Pi/32) = sin(25*I*Pi/32) */	tmp += 2;	// tmp+4
		tmp64 = 0x3FE6A09E667F3BCDull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(08*I*Pi/32) = sin(24*I*Pi/32) */	tm2 -= 2;	// tmp+3
		tmp64 = 0x3FE44CF325091DD6ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(09*I*Pi/32) = sin(23*I*Pi/32) */
		tmp64 = 0x3FE1C73B39AE68C8ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(10*I*Pi/32) = sin(22*I*Pi/32) */
		tmp64 = 0x3FDE2B5D3806F63Bull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(11*I*Pi/32) = sin(21*I*Pi/32) */	tmp += 2;	// tmp+6
		tmp64 = 0x3FD87DE2A6AEA963ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(12*I*Pi/32) = sin(20*I*Pi/32) */	tm2 -= 2;	// tmp+1
		tmp64 = 0x3FD294062ED59F06ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(13*I*Pi/32) = sin(19*I*Pi/32) */
		tmp64 = 0x3FC8F8B83C69A60Bull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(14*I*Pi/32) = sin(18*I*Pi/32) */
		tmp64 = 0x3FB917A6BC29B42Cull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(15*I*Pi/32) = sin(17*I*Pi/32) */	tmp += 2;	// tmp+8 (unused)

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
	n_minus_sil   = (struct uint32x8 *)sse_nm1 + 1;
	n_minus_silp1 = (struct uint32x8 *)sse_nm1 + 2;
	sinwt         = (struct uint32x8 *)sse_nm1 + 3;
	sinwtm1       = (struct uint32x8 *)sse_nm1 + 4;
	nbytes += 128;
#elif defined(USE_AVX)
	n_minus_sil   = (struct uint32x4 *)sse_nm1 + 1;
	n_minus_silp1 = (struct uint32x4 *)sse_nm1 + 2;
	sinwt         = (struct uint32x4 *)sse_nm1 + 3;
	sinwtm1       = (struct uint32x4 *)sse_nm1 + 4;
	nbytes += 64;
#endif

	// Propagate the above consts to the remaining threads:
	tmp = (vec_dbl *)sm_ptr;
	tm2 = tmp + cslots_in_local_store;
	for(ithread = 1; ithread < CY_THREADS; ++ithread) {
		memcpy(tm2, tmp, nbytes);
		tmp = tm2;		tm2 += cslots_in_local_store;
	}

	#ifdef USE_AVX
		bjmodn0 = (int*)(sinwtm1 + RE_IM_STRIDE);
	#else
		bjmodn0 = (int*)(sse_nm1 + RE_IM_STRIDE);
	#endif
		bjmodn1 = bjmodn0 + 1;
		bjmodn2 = bjmodn1 + 1;
		bjmodn3 = bjmodn2 + 1;
		bjmodn4 = bjmodn3 + 1;
		bjmodn5 = bjmodn4 + 1;
		bjmodn6 = bjmodn5 + 1;
		bjmodn7 = bjmodn6 + 1;
		bjmodn8 = bjmodn7 + 1;
		bjmodn9 = bjmodn8 + 1;
		bjmodnA = bjmodn9 + 1;
		bjmodnB = bjmodnA + 1;
		bjmodnC = bjmodnB + 1;
		bjmodnD = bjmodnC + 1;
		bjmodnE = bjmodnD + 1;
		bjmodnF = bjmodnE + 1;

	#endif	// USE_SSE2

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
		#ifdef USE_SSE2
			tdat[ithread].r00 = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (vec_dbl *)((intptr_t)tdat[ithread].r00 + ((intptr_t)half_arr - (intptr_t)r00));
		#else	// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
			tdat[ithread].r00      = (double *)base;
			tdat[ithread].half_arr = (double *)baseinv;
		#endif	// USE_SSE2
		}
	#endif

		/*   constant index offsets for load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		p1  = NDIVR;
		p2  = p1  + p1;
		p3  = p2  + p1;
		p4  = p3  + p1;
		p5  = p4  + p1;
		p6  = p5  + p1;
		p7  = p6  + p1;
		p8  = p7  + p1;
		p9  = p8  + p1;
		p10 = p9  + p1;
		p11 = p10 + p1;
		p12 = p11 + p1;
		p13 = p12 + p1;
		p14 = p13 + p1;
		p15 = p14 + p1;

		p1  = p1  + ( (p1  >> DAT_BITS) << PAD_BITS );
		p2  = p2  + ( (p2  >> DAT_BITS) << PAD_BITS );
		p3  = p3  + ( (p3  >> DAT_BITS) << PAD_BITS );
		p4  = p4  + ( (p4  >> DAT_BITS) << PAD_BITS );
		p5  = p5  + ( (p5  >> DAT_BITS) << PAD_BITS );
		p6  = p6  + ( (p6  >> DAT_BITS) << PAD_BITS );
		p7  = p7  + ( (p7  >> DAT_BITS) << PAD_BITS );
		p8  = p8  + ( (p8  >> DAT_BITS) << PAD_BITS );
		p9  = p9  + ( (p9  >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );

		p0123[0] = 0; p0123[1] = p1; p0123[2] = p2; p0123[3] = p3;
		poff[0x0] =   0; poff[0x1] = p4    ; poff[0x2] = p8; poff[0x3] = p12;

		if(_cy_r0)	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;

			free((void *)_bjmodn0); _bjmodn0 = 0x0;
			free((void *)_bjmodn1); _bjmodn1 = 0x0;
			free((void *)_bjmodn2); _bjmodn2 = 0x0;
			free((void *)_bjmodn3); _bjmodn3 = 0x0;
			free((void *)_bjmodn4); _bjmodn4 = 0x0;
			free((void *)_bjmodn5); _bjmodn5 = 0x0;
			free((void *)_bjmodn6); _bjmodn6 = 0x0;
			free((void *)_bjmodn7); _bjmodn7 = 0x0;
			free((void *)_bjmodn8); _bjmodn8 = 0x0;
			free((void *)_bjmodn9); _bjmodn9 = 0x0;
			free((void *)_bjmodnA); _bjmodnA = 0x0;
			free((void *)_bjmodnB); _bjmodnB = 0x0;
			free((void *)_bjmodnC); _bjmodnC = 0x0;
			free((void *)_bjmodnD); _bjmodnD = 0x0;
			free((void *)_bjmodnE); _bjmodnE = 0x0;
			free((void *)_bjmodnF); _bjmodnF = 0x0;

			free((void *)_cy_r0); _cy_r0 = 0x0;	free((void *)_cy_i0); _cy_i0 = 0x0;
			free((void *)_cy_r1); _cy_r1 = 0x0;	free((void *)_cy_i1); _cy_i1 = 0x0;
			free((void *)_cy_r2); _cy_r2 = 0x0;	free((void *)_cy_i2); _cy_i2 = 0x0;
			free((void *)_cy_r3); _cy_r3 = 0x0;	free((void *)_cy_i3); _cy_i3 = 0x0;
			free((void *)_cy_r4); _cy_r4 = 0x0;	free((void *)_cy_i4); _cy_i4 = 0x0;
			free((void *)_cy_r5); _cy_r5 = 0x0;	free((void *)_cy_i5); _cy_i5 = 0x0;
			free((void *)_cy_r6); _cy_r6 = 0x0;	free((void *)_cy_i6); _cy_i6 = 0x0;
			free((void *)_cy_r7); _cy_r7 = 0x0;	free((void *)_cy_i7); _cy_i7 = 0x0;
			free((void *)_cy_r8); _cy_r8 = 0x0;	free((void *)_cy_i8); _cy_i8 = 0x0;
			free((void *)_cy_r9); _cy_r9 = 0x0;	free((void *)_cy_i9); _cy_i9 = 0x0;
			free((void *)_cy_rA); _cy_rA = 0x0;	free((void *)_cy_iA); _cy_iA = 0x0;
			free((void *)_cy_rB); _cy_rB = 0x0;	free((void *)_cy_iB); _cy_iB = 0x0;
			free((void *)_cy_rC); _cy_rC = 0x0;	free((void *)_cy_iC); _cy_iC = 0x0;
			free((void *)_cy_rD); _cy_rD = 0x0;	free((void *)_cy_iD); _cy_iD = 0x0;
			free((void *)_cy_rE); _cy_rE = 0x0;	free((void *)_cy_iE); _cy_iE = 0x0;
			free((void *)_cy_rF); _cy_rF = 0x0;	free((void *)_cy_iF); _cy_iF = 0x0;

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
		_bjmodn0	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn0== 0x0);
		_bjmodn1	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn1== 0x0);
		_bjmodn2	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn2== 0x0);
		_bjmodn3	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn3== 0x0);
		_bjmodn4	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn4== 0x0);
		_bjmodn5	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn5== 0x0);
		_bjmodn6	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn6== 0x0);
		_bjmodn7	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn7== 0x0);
		_bjmodn8	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn8== 0x0);
		_bjmodn9	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn9== 0x0);
		_bjmodnA	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodnA== 0x0);
		_bjmodnB	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodnB== 0x0);
		_bjmodnC	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodnC== 0x0);
		_bjmodnD	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodnD== 0x0);
		_bjmodnE	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodnE== 0x0);
		_bjmodnF	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodnF== 0x0);
		_jstart  	= (int *)malloc(j);	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(j);	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(j);	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(j);	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(j);	ptr_prod += (uint32)(_co3     == 0x0);

		j = CY_THREADS*sizeof(double);
		_cy_r0	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r0== 0x0);
		_cy_r1	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r1== 0x0);
		_cy_r2	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r2== 0x0);
		_cy_r3	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r3== 0x0);
		_cy_r4	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r4== 0x0);
		_cy_r5	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r5== 0x0);
		_cy_r6	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r6== 0x0);
		_cy_r7	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r7== 0x0);
		_cy_r8	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r8== 0x0);
		_cy_r9	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r9== 0x0);
		_cy_rA	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_rA== 0x0);
		_cy_rB	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_rB== 0x0);
		_cy_rC	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_rC== 0x0);
		_cy_rD	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_rD== 0x0);
		_cy_rE	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_rE== 0x0);
		_cy_rF	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_rF== 0x0);

		_cy_i0	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i0== 0x0);
		_cy_i1	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i1== 0x0);
		_cy_i2	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i2== 0x0);
		_cy_i3	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i3== 0x0);
		_cy_i4	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i4== 0x0);
		_cy_i5	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i5== 0x0);
		_cy_i6	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i6== 0x0);
		_cy_i7	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i7== 0x0);
		_cy_i8	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i8== 0x0);
		_cy_i9	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i9== 0x0);
		_cy_iA	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_iA== 0x0);
		_cy_iB	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_iB== 0x0);
		_cy_iC	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_iC== 0x0);
		_cy_iD	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_iD== 0x0);
		_cy_iE	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_iE== 0x0);
		_cy_iF	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_iF== 0x0);

		ASSERT(ptr_prod == 0, "ERROR: unable to allocate one or more auxiliary arrays.");

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
			i.e. the one that n2/16-separated FFT outputs need:
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
		}

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

/*...The radix-16 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_cy_r0[ithread] = 0;	_cy_i0[ithread] = 0;
		_cy_r1[ithread] = 0;	_cy_i1[ithread] = 0;
		_cy_r2[ithread] = 0;	_cy_i2[ithread] = 0;
		_cy_r3[ithread] = 0;	_cy_i3[ithread] = 0;
		_cy_r4[ithread] = 0;	_cy_i4[ithread] = 0;
		_cy_r5[ithread] = 0;	_cy_i5[ithread] = 0;
		_cy_r6[ithread] = 0;	_cy_i6[ithread] = 0;
		_cy_r7[ithread] = 0;	_cy_i7[ithread] = 0;
		_cy_r8[ithread] = 0;	_cy_i8[ithread] = 0;
		_cy_r9[ithread] = 0;	_cy_i9[ithread] = 0;
		_cy_rA[ithread] = 0;	_cy_iA[ithread] = 0;
		_cy_rB[ithread] = 0;	_cy_iB[ithread] = 0;
		_cy_rC[ithread] = 0;	_cy_iC[ithread] = 0;
		_cy_rD[ithread] = 0;	_cy_iD[ithread] = 0;
		_cy_rE[ithread] = 0;	_cy_iE[ithread] = 0;
		_cy_rF[ithread] = 0;	_cy_iF[ithread] = 0;
	}
  #if 0	//ndef USE_SSE2	*** v20: Non-SIMD builds now also support shifted-residue
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy_r0[0] = -2;
	}
  #endif
	*fracmax = 0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	// init inverse-weight scale factor = 2/n for normal carry pass, 1 for wrapper pass

for(outer=0; outer <= 1; outer++)
{
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (_i[0] = 1).	*/

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

		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_bjmodn0[ithread] = _bjmodnini[ithread];	j = _bjmodnini[CY_THREADS];
			MOD_ADD32(_bjmodn0[ithread], j, n, _bjmodn1[ithread]);
			MOD_ADD32(_bjmodn1[ithread], j, n, _bjmodn2[ithread]);
			MOD_ADD32(_bjmodn2[ithread], j, n, _bjmodn3[ithread]);
			MOD_ADD32(_bjmodn3[ithread], j, n, _bjmodn4[ithread]);
			MOD_ADD32(_bjmodn4[ithread], j, n, _bjmodn5[ithread]);
			MOD_ADD32(_bjmodn5[ithread], j, n, _bjmodn6[ithread]);
			MOD_ADD32(_bjmodn6[ithread], j, n, _bjmodn7[ithread]);
			MOD_ADD32(_bjmodn7[ithread], j, n, _bjmodn8[ithread]);
			MOD_ADD32(_bjmodn8[ithread], j, n, _bjmodn9[ithread]);
			MOD_ADD32(_bjmodn9[ithread], j, n, _bjmodnA[ithread]);
			MOD_ADD32(_bjmodnA[ithread], j, n, _bjmodnB[ithread]);
			MOD_ADD32(_bjmodnB[ithread], j, n, _bjmodnC[ithread]);
			MOD_ADD32(_bjmodnC[ithread], j, n, _bjmodnD[ithread]);
			MOD_ADD32(_bjmodnD[ithread], j, n, _bjmodnE[ithread]);
			MOD_ADD32(_bjmodnE[ithread], j, n, _bjmodnF[ithread]);

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

#endif

	/* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
	if(!full_pass)
	{
		khi = 1;
	}

#ifdef USE_PTHREAD
	/* Populate the thread-specific data structs - use the invariant terms as memchecks: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
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
		// Dec 2015: fast-GCD usage of this routine may involve multiple 'main' arrays
		// on successive calls, so set here at runtime rather than in init-only block:
		tdat[ithread].arrdat = a;			/* Main data array */
	#ifdef USE_FGT61
		ASSERT(tdat[ithread].brrdat == b, "thread-local memcheck fail!");			/* Modular version of main data array */
	#endif
		ASSERT(tdat[ithread].wt0 == wt0, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].wt1 == wt1, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].si  == si, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].rn0 == rn0, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].rn1 == rn1, "thread-local memcheck fail!");
	#ifdef USE_SSE2
		ASSERT(tdat[ithread].r00 == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].r00;
		ASSERT(((tmp + 0x20)->d0 == ISRT2 && (tmp + 0x20)->d1 == ISRT2), "thread-local memcheck failed!");
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
			tdat[ithread].bjmodn0 = _bjmodn0[ithread];
			tdat[ithread].bjmodn1 = _bjmodn1[ithread];
			tdat[ithread].bjmodn2 = _bjmodn2[ithread];
			tdat[ithread].bjmodn3 = _bjmodn3[ithread];
			tdat[ithread].bjmodn4 = _bjmodn4[ithread];
			tdat[ithread].bjmodn5 = _bjmodn5[ithread];
			tdat[ithread].bjmodn6 = _bjmodn6[ithread];
			tdat[ithread].bjmodn7 = _bjmodn7[ithread];
			tdat[ithread].bjmodn8 = _bjmodn8[ithread];
			tdat[ithread].bjmodn9 = _bjmodn9[ithread];
			tdat[ithread].bjmodnA = _bjmodnA[ithread];
			tdat[ithread].bjmodnB = _bjmodnB[ithread];
			tdat[ithread].bjmodnC = _bjmodnC[ithread];
			tdat[ithread].bjmodnD = _bjmodnD[ithread];
			tdat[ithread].bjmodnE = _bjmodnE[ithread];
			tdat[ithread].bjmodnF = _bjmodnF[ithread];
			/* init carries	*/
			tdat[ithread].cy_r0 = _cy_r0[ithread];
			tdat[ithread].cy_r1 = _cy_r1[ithread];
			tdat[ithread].cy_r2 = _cy_r2[ithread];
			tdat[ithread].cy_r3 = _cy_r3[ithread];
			tdat[ithread].cy_r4 = _cy_r4[ithread];
			tdat[ithread].cy_r5 = _cy_r5[ithread];
			tdat[ithread].cy_r6 = _cy_r6[ithread];
			tdat[ithread].cy_r7 = _cy_r7[ithread];
			tdat[ithread].cy_r8 = _cy_r8[ithread];
			tdat[ithread].cy_r9 = _cy_r9[ithread];
			tdat[ithread].cy_rA = _cy_rA[ithread];
			tdat[ithread].cy_rB = _cy_rB[ithread];
			tdat[ithread].cy_rC = _cy_rC[ithread];
			tdat[ithread].cy_rD = _cy_rD[ithread];
			tdat[ithread].cy_rE = _cy_rE[ithread];
			tdat[ithread].cy_rF = _cy_rF[ithread];
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
		#ifdef USE_AVX512
			/* No-Op */
		#elif defined(USE_SSE2)
			// This is slightly different for power-of-2 DFTs: Here, scale is in the +2 slot, base & baseinv remain fixed in 0,+1 slots:
			dtmp = tmp->d0 * (tmp+1)->d0;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = tmp->d1 * (tmp+1)->d1;	ASSERT(fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			// scale gets set immediately prior to calling carry macro, hence no use checking it here.
		#endif
			/* init carries	*/
			tdat[ithread].cy_r0 = _cy_r0[ithread];	tdat[ithread].cy_i0 = _cy_i0[ithread];
			tdat[ithread].cy_r1 = _cy_r1[ithread];	tdat[ithread].cy_i1 = _cy_i1[ithread];
			tdat[ithread].cy_r2 = _cy_r2[ithread];	tdat[ithread].cy_i2 = _cy_i2[ithread];
			tdat[ithread].cy_r3 = _cy_r3[ithread];	tdat[ithread].cy_i3 = _cy_i3[ithread];
			tdat[ithread].cy_r4 = _cy_r4[ithread];	tdat[ithread].cy_i4 = _cy_i4[ithread];
			tdat[ithread].cy_r5 = _cy_r5[ithread];	tdat[ithread].cy_i5 = _cy_i5[ithread];
			tdat[ithread].cy_r6 = _cy_r6[ithread];	tdat[ithread].cy_i6 = _cy_i6[ithread];
			tdat[ithread].cy_r7 = _cy_r7[ithread];	tdat[ithread].cy_i7 = _cy_i7[ithread];
			tdat[ithread].cy_r8 = _cy_r8[ithread];	tdat[ithread].cy_i8 = _cy_i8[ithread];
			tdat[ithread].cy_r9 = _cy_r9[ithread];	tdat[ithread].cy_i9 = _cy_i9[ithread];
			tdat[ithread].cy_rA = _cy_rA[ithread];	tdat[ithread].cy_iA = _cy_iA[ithread];
			tdat[ithread].cy_rB = _cy_rB[ithread];	tdat[ithread].cy_iB = _cy_iB[ithread];
			tdat[ithread].cy_rC = _cy_rC[ithread];	tdat[ithread].cy_iC = _cy_iC[ithread];
			tdat[ithread].cy_rD = _cy_rD[ithread];	tdat[ithread].cy_iD = _cy_iD[ithread];
			tdat[ithread].cy_rE = _cy_rE[ithread];	tdat[ithread].cy_iE = _cy_iE[ithread];
			tdat[ithread].cy_rF = _cy_rF[ithread];	tdat[ithread].cy_iF = _cy_iF[ithread];
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

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_SSE2
			*bjmodn0 = _bjmodn0[ithread];
			*bjmodn1 = _bjmodn1[ithread];
			*bjmodn2 = _bjmodn2[ithread];
			*bjmodn3 = _bjmodn3[ithread];
			*bjmodn4 = _bjmodn4[ithread];
			*bjmodn5 = _bjmodn5[ithread];
			*bjmodn6 = _bjmodn6[ithread];
			*bjmodn7 = _bjmodn7[ithread];
			*bjmodn8 = _bjmodn8[ithread];
			*bjmodn9 = _bjmodn9[ithread];
			*bjmodnA = _bjmodnA[ithread];
			*bjmodnB = _bjmodnB[ithread];
			*bjmodnC = _bjmodnC[ithread];
			*bjmodnD = _bjmodnD[ithread];
			*bjmodnE = _bjmodnE[ithread];
			*bjmodnF = _bjmodnF[ithread];
		#else
			bjmodn0 = _bjmodn0[ithread];
			bjmodn1 = _bjmodn1[ithread];
			bjmodn2 = _bjmodn2[ithread];
			bjmodn3 = _bjmodn3[ithread];
			bjmodn4 = _bjmodn4[ithread];
			bjmodn5 = _bjmodn5[ithread];
			bjmodn6 = _bjmodn6[ithread];
			bjmodn7 = _bjmodn7[ithread];
			bjmodn8 = _bjmodn8[ithread];
			bjmodn9 = _bjmodn9[ithread];
			bjmodnA = _bjmodnA[ithread];
			bjmodnB = _bjmodnB[ithread];
			bjmodnC = _bjmodnC[ithread];
			bjmodnD = _bjmodnD[ithread];
			bjmodnE = _bjmodnE[ithread];
			bjmodnF = _bjmodnF[ithread];
		#endif
			/* init carries	*/
		#ifdef USE_AVX512
			cy_r0->d0 = _cy_r0[ithread];	cy_r0->d1 = _cy_r1[ithread];	cy_r0->d2 = _cy_r2[ithread];	cy_r0->d3 = _cy_r3[ithread];	cy_r0->d4 = _cy_r4[ithread];	cy_r0->d5 = _cy_r5[ithread];	cy_r0->d6 = _cy_r6[ithread];	cy_r0->d7 = _cy_r7[ithread];
			cy_r8->d0 = _cy_r8[ithread];	cy_r8->d1 = _cy_r9[ithread];	cy_r8->d2 = _cy_rA[ithread];	cy_r8->d3 = _cy_rB[ithread];	cy_r8->d4 = _cy_rC[ithread];	cy_r8->d5 = _cy_rD[ithread];	cy_r8->d6 = _cy_rE[ithread];	cy_r8->d7 = _cy_rF[ithread];
		#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers
			cy_r0->d0 = _cy_r0[ithread];	cy_r0->d1 = _cy_r1[ithread];	cy_r0->d2 = _cy_r2[ithread];	cy_r0->d3 = _cy_r3[ithread];
			cy_r4->d0 = _cy_r4[ithread];	cy_r4->d1 = _cy_r5[ithread];	cy_r4->d2 = _cy_r6[ithread];	cy_r4->d3 = _cy_r7[ithread];
			cy_r8->d0 = _cy_r8[ithread];	cy_r8->d1 = _cy_r9[ithread];	cy_r8->d2 = _cy_rA[ithread];	cy_r8->d3 = _cy_rB[ithread];
			cy_rC->d0 = _cy_rC[ithread];	cy_rC->d1 = _cy_rD[ithread];	cy_rC->d2 = _cy_rE[ithread];	cy_rC->d3 = _cy_rF[ithread];
		#elif defined(USE_SSE2)
			cy_r0->d0 = _cy_r0[ithread];	cy_r0->d1 = _cy_r1[ithread];
			cy_r2->d0 = _cy_r2[ithread];	cy_r2->d1 = _cy_r3[ithread];
			cy_r4->d0 = _cy_r4[ithread];	cy_r4->d1 = _cy_r5[ithread];
			cy_r6->d0 = _cy_r6[ithread];	cy_r6->d1 = _cy_r7[ithread];
			cy_r8->d0 = _cy_r8[ithread];	cy_r8->d1 = _cy_r9[ithread];
			cy_rA->d0 = _cy_rA[ithread];	cy_rA->d1 = _cy_rB[ithread];
			cy_rC->d0 = _cy_rC[ithread];	cy_rC->d1 = _cy_rD[ithread];
			cy_rE->d0 = _cy_rE[ithread];	cy_rE->d1 = _cy_rF[ithread];
		#else
			cy_r0 = _cy_r0[ithread];
			cy_r1 = _cy_r1[ithread];
			cy_r2 = _cy_r2[ithread];
			cy_r3 = _cy_r3[ithread];
			cy_r4 = _cy_r4[ithread];
			cy_r5 = _cy_r5[ithread];
			cy_r6 = _cy_r6[ithread];
			cy_r7 = _cy_r7[ithread];
			cy_r8 = _cy_r8[ithread];
			cy_r9 = _cy_r9[ithread];
			cy_rA = _cy_rA[ithread];
			cy_rB = _cy_rB[ithread];
			cy_rC = _cy_rC[ithread];
			cy_rD = _cy_rD[ithread];
			cy_rE = _cy_rE[ithread];
			cy_rF = _cy_rF[ithread];
		#endif
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
			/* init carries	*/
		#ifdef USE_AVX512
			cy_r0->d0 = _cy_r0[ithread];	cy_r0->d1 = _cy_r1[ithread];	cy_r0->d2 = _cy_r2[ithread];	cy_r0->d3 = _cy_r3[ithread];	cy_r0->d4 = _cy_r4[ithread];	cy_r0->d5 = _cy_r5[ithread];	cy_r0->d6 = _cy_r6[ithread];	cy_r0->d7 = _cy_r7[ithread];
			cy_r8->d0 = _cy_r8[ithread];	cy_r8->d1 = _cy_r9[ithread];	cy_r8->d2 = _cy_rA[ithread];	cy_r8->d3 = _cy_rB[ithread];	cy_r8->d4 = _cy_rC[ithread];	cy_r8->d5 = _cy_rD[ithread];	cy_r8->d6 = _cy_rE[ithread];	cy_r8->d7 = _cy_rF[ithread];

			cy_i0->d0 = _cy_i0[ithread];	cy_i0->d1 = _cy_i1[ithread];	cy_i0->d2 = _cy_i2[ithread];	cy_i0->d3 = _cy_i3[ithread];	cy_i0->d4 = _cy_i4[ithread];	cy_i0->d5 = _cy_i5[ithread];	cy_i0->d6 = _cy_i6[ithread];	cy_i0->d7 = _cy_i7[ithread];
			cy_i8->d0 = _cy_i8[ithread];	cy_i8->d1 = _cy_i9[ithread];	cy_i8->d2 = _cy_iA[ithread];	cy_i8->d3 = _cy_iB[ithread];	cy_i8->d4 = _cy_iC[ithread];	cy_i8->d5 = _cy_iD[ithread];	cy_i8->d6 = _cy_iE[ithread];	cy_i8->d7 = _cy_iF[ithread];
		#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers
			cy_r0->d0 = _cy_r0[ithread];	cy_r0->d1 = _cy_r1[ithread];	cy_r0->d2 = _cy_r2[ithread];	cy_r0->d3 = _cy_r3[ithread];
			cy_r4->d0 = _cy_r4[ithread];	cy_r4->d1 = _cy_r5[ithread];	cy_r4->d2 = _cy_r6[ithread];	cy_r4->d3 = _cy_r7[ithread];
			cy_r8->d0 = _cy_r8[ithread];	cy_r8->d1 = _cy_r9[ithread];	cy_r8->d2 = _cy_rA[ithread];	cy_r8->d3 = _cy_rB[ithread];
			cy_rC->d0 = _cy_rC[ithread];	cy_rC->d1 = _cy_rD[ithread];	cy_rC->d2 = _cy_rE[ithread];	cy_rC->d3 = _cy_rF[ithread];

			cy_i0->d0 = _cy_i0[ithread];	cy_i0->d1 = _cy_i1[ithread];	cy_i0->d2 = _cy_i2[ithread];	cy_i0->d3 = _cy_i3[ithread];
			cy_i4->d0 = _cy_i4[ithread];	cy_i4->d1 = _cy_i5[ithread];	cy_i4->d2 = _cy_i6[ithread];	cy_i4->d3 = _cy_i7[ithread];
			cy_i8->d0 = _cy_i8[ithread];	cy_i8->d1 = _cy_i9[ithread];	cy_i8->d2 = _cy_iA[ithread];	cy_i8->d3 = _cy_iB[ithread];
			cy_iC->d0 = _cy_iC[ithread];	cy_iC->d1 = _cy_iD[ithread];	cy_iC->d2 = _cy_iE[ithread];	cy_iC->d3 = _cy_iF[ithread];
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			cy_r0->d0 = _cy_r0[ithread];	cy_r0->d1 = _cy_i0[ithread];
			cy_r2->d0 = _cy_r1[ithread];	cy_r2->d1 = _cy_i1[ithread];
			cy_r4->d0 = _cy_r2[ithread];	cy_r4->d1 = _cy_i2[ithread];
			cy_r6->d0 = _cy_r3[ithread];	cy_r6->d1 = _cy_i3[ithread];
			cy_r8->d0 = _cy_r4[ithread];	cy_r8->d1 = _cy_i4[ithread];
			cy_rA->d0 = _cy_r5[ithread];	cy_rA->d1 = _cy_i5[ithread];
			cy_rC->d0 = _cy_r6[ithread];	cy_rC->d1 = _cy_i6[ithread];
			cy_rE->d0 = _cy_r7[ithread];	cy_rE->d1 = _cy_i7[ithread];
			cy_i0->d0 = _cy_r8[ithread];	cy_i0->d1 = _cy_i8[ithread];
			cy_i2->d0 = _cy_r9[ithread];	cy_i2->d1 = _cy_i9[ithread];
			cy_i4->d0 = _cy_rA[ithread];	cy_i4->d1 = _cy_iA[ithread];
			cy_i6->d0 = _cy_rB[ithread];	cy_i6->d1 = _cy_iB[ithread];
			cy_i8->d0 = _cy_rC[ithread];	cy_i8->d1 = _cy_iC[ithread];
			cy_iA->d0 = _cy_rD[ithread];	cy_iA->d1 = _cy_iD[ithread];
			cy_iC->d0 = _cy_rE[ithread];	cy_iC->d1 = _cy_iE[ithread];
			cy_iE->d0 = _cy_rF[ithread];	cy_iE->d1 = _cy_iF[ithread];
		#else
			cy_r0 = _cy_r0[ithread];
			cy_r1 = _cy_r1[ithread];
			cy_r2 = _cy_r2[ithread];
			cy_r3 = _cy_r3[ithread];
			cy_r4 = _cy_r4[ithread];
			cy_r5 = _cy_r5[ithread];
			cy_r6 = _cy_r6[ithread];
			cy_r7 = _cy_r7[ithread];
			cy_r8 = _cy_r8[ithread];
			cy_r9 = _cy_r9[ithread];
			cy_rA = _cy_rA[ithread];
			cy_rB = _cy_rB[ithread];
			cy_rC = _cy_rC[ithread];
			cy_rD = _cy_rD[ithread];
			cy_rE = _cy_rE[ithread];
			cy_rF = _cy_rF[ithread];

			cy_i0 = _cy_i0[ithread];
			cy_i1 = _cy_i1[ithread];
			cy_i2 = _cy_i2[ithread];
			cy_i3 = _cy_i3[ithread];
			cy_i4 = _cy_i4[ithread];
			cy_i5 = _cy_i5[ithread];
			cy_i6 = _cy_i6[ithread];
			cy_i7 = _cy_i7[ithread];
			cy_i8 = _cy_i8[ithread];
			cy_i9 = _cy_i9[ithread];
			cy_iA = _cy_iA[ithread];
			cy_iB = _cy_iB[ithread];
			cy_iC = _cy_iC[ithread];
			cy_iD = _cy_iD[ithread];
			cy_iE = _cy_iE[ithread];
			cy_iF = _cy_iF[ithread];
		#endif
		}

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix16_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX512
			_cy_r0[ithread] = cy_r0->d0;	_cy_r1[ithread] = cy_r0->d1;	_cy_r2[ithread] = cy_r0->d2;	_cy_r3[ithread] = cy_r0->d3;	_cy_r4[ithread] = cy_r0->d4;	_cy_r5[ithread] = cy_r0->d5;	_cy_r6[ithread] = cy_r0->d6;	_cy_r7[ithread] = cy_r0->d7;
			_cy_r8[ithread] = cy_r8->d0;	_cy_r9[ithread] = cy_r8->d1;	_cy_rA[ithread] = cy_r8->d2;	_cy_rB[ithread] = cy_r8->d3;	_cy_rC[ithread] = cy_r8->d4;	_cy_rD[ithread] = cy_r8->d5;	_cy_rE[ithread] = cy_r8->d6;	_cy_rF[ithread] = cy_r8->d7;
			if(full_pass) {
				t0 = MAX(max_err->d0,max_err->d1);
				t1 = MAX(max_err->d2,max_err->d3);
				t2 = MAX(max_err->d4,max_err->d5);
				t3 = MAX(max_err->d6,max_err->d7);
				maxerr = MAX( MAX(t0,t1), MAX(t2,t3) );
			}
		#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers
			_cy_r0[ithread] = cy_r0->d0;	_cy_r1[ithread] = cy_r0->d1;	_cy_r2[ithread] = cy_r0->d2;	_cy_r3[ithread] = cy_r0->d3;
			_cy_r4[ithread] = cy_r4->d0;	_cy_r5[ithread] = cy_r4->d1;	_cy_r6[ithread] = cy_r4->d2;	_cy_r7[ithread] = cy_r4->d3;
			_cy_r8[ithread] = cy_r8->d0;	_cy_r9[ithread] = cy_r8->d1;	_cy_rA[ithread] = cy_r8->d2;	_cy_rB[ithread] = cy_r8->d3;
			_cy_rC[ithread] = cy_rC->d0;	_cy_rD[ithread] = cy_rC->d1;	_cy_rE[ithread] = cy_rC->d2;	_cy_rF[ithread] = cy_rC->d3;
			if(full_pass) maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			_cy_r0[ithread] = cy_r0->d0;	_cy_r1[ithread] = cy_r0->d1;
			_cy_r2[ithread] = cy_r2->d0;	_cy_r3[ithread] = cy_r2->d1;
			_cy_r4[ithread] = cy_r4->d0;	_cy_r5[ithread] = cy_r4->d1;
			_cy_r6[ithread] = cy_r6->d0;	_cy_r7[ithread] = cy_r6->d1;
			_cy_r8[ithread] = cy_r8->d0;	_cy_r9[ithread] = cy_r8->d1;
			_cy_rA[ithread] = cy_rA->d0;	_cy_rB[ithread] = cy_rA->d1;
			_cy_rC[ithread] = cy_rC->d0;	_cy_rD[ithread] = cy_rC->d1;
			_cy_rE[ithread] = cy_rE->d0;	_cy_rF[ithread] = cy_rE->d1;
			if(full_pass) maxerr = MAX(max_err->d0,max_err->d1);
		#else
			_cy_r0[ithread] = cy_r0;
			_cy_r1[ithread] = cy_r1;
			_cy_r2[ithread] = cy_r2;
			_cy_r3[ithread] = cy_r3;
			_cy_r4[ithread] = cy_r4;
			_cy_r5[ithread] = cy_r5;
			_cy_r6[ithread] = cy_r6;
			_cy_r7[ithread] = cy_r7;
			_cy_r8[ithread] = cy_r8;
			_cy_r9[ithread] = cy_r9;
			_cy_rA[ithread] = cy_rA;
			_cy_rB[ithread] = cy_rB;
			_cy_rC[ithread] = cy_rC;
			_cy_rD[ithread] = cy_rD;
			_cy_rE[ithread] = cy_rE;
			_cy_rF[ithread] = cy_rF;
		#endif
		}
		else
		{
		#ifdef USE_AVX512
			_cy_r0[ithread] = cy_r0->d0;	_cy_r1[ithread] = cy_r0->d1;	_cy_r2[ithread] = cy_r0->d2;	_cy_r3[ithread] = cy_r0->d3;	_cy_r4[ithread] = cy_r0->d4;	_cy_r5[ithread] = cy_r0->d5;	_cy_r6[ithread] = cy_r0->d6;	_cy_r7[ithread] = cy_r0->d7;
			_cy_r8[ithread] = cy_r8->d0;	_cy_r9[ithread] = cy_r8->d1;	_cy_rA[ithread] = cy_r8->d2;	_cy_rB[ithread] = cy_r8->d3;	_cy_rC[ithread] = cy_r8->d4;	_cy_rD[ithread] = cy_r8->d5;	_cy_rE[ithread] = cy_r8->d6;	_cy_rF[ithread] = cy_r8->d7;

			_cy_i0[ithread] = cy_i0->d0;	_cy_i1[ithread] = cy_i0->d1;	_cy_i2[ithread] = cy_i0->d2;	_cy_i3[ithread] = cy_i0->d3;	_cy_i4[ithread] = cy_i0->d4;	_cy_i5[ithread] = cy_i0->d5;	_cy_i6[ithread] = cy_i0->d6;	_cy_i7[ithread] = cy_i0->d7;
			_cy_i8[ithread] = cy_i8->d0;	_cy_i9[ithread] = cy_i8->d1;	_cy_iA[ithread] = cy_i8->d2;	_cy_iB[ithread] = cy_i8->d3;	_cy_iC[ithread] = cy_i8->d4;	_cy_iD[ithread] = cy_i8->d5;	_cy_iE[ithread] = cy_i8->d6;	_cy_iF[ithread] = cy_i8->d7;
			if(full_pass) {
				t0 = MAX(max_err->d0,max_err->d1);
				t1 = MAX(max_err->d2,max_err->d3);
				t2 = MAX(max_err->d4,max_err->d5);
				t3 = MAX(max_err->d6,max_err->d7);
				maxerr = MAX( MAX(t0,t1), MAX(t2,t3) );
			}
		#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers
			_cy_r0[ithread] = cy_r0->d0;	_cy_r1[ithread] = cy_r0->d1;	_cy_r2[ithread] = cy_r0->d2;	_cy_r3[ithread] = cy_r0->d3;
			_cy_r4[ithread] = cy_r4->d0;	_cy_r5[ithread] = cy_r4->d1;	_cy_r6[ithread] = cy_r4->d2;	_cy_r7[ithread] = cy_r4->d3;
			_cy_r8[ithread] = cy_r8->d0;	_cy_r9[ithread] = cy_r8->d1;	_cy_rA[ithread] = cy_r8->d2;	_cy_rB[ithread] = cy_r8->d3;
			_cy_rC[ithread] = cy_rC->d0;	_cy_rD[ithread] = cy_rC->d1;	_cy_rE[ithread] = cy_rC->d2;	_cy_rF[ithread] = cy_rC->d3;

			_cy_i0[ithread] = cy_i0->d0;	_cy_i1[ithread] = cy_i0->d1;	_cy_i2[ithread] = cy_i0->d2;	_cy_i3[ithread] = cy_i0->d3;
			_cy_i4[ithread] = cy_i4->d0;	_cy_i5[ithread] = cy_i4->d1;	_cy_i6[ithread] = cy_i4->d2;	_cy_i7[ithread] = cy_i4->d3;
			_cy_i8[ithread] = cy_i8->d0;	_cy_i9[ithread] = cy_i8->d1;	_cy_iA[ithread] = cy_i8->d2;	_cy_iB[ithread] = cy_i8->d3;
			_cy_iC[ithread] = cy_iC->d0;	_cy_iD[ithread] = cy_iC->d1;	_cy_iE[ithread] = cy_iC->d2;	_cy_iF[ithread] = cy_iC->d3;
			if(full_pass) maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			_cy_r0[ithread] = cy_r0->d0;	_cy_i0[ithread] = cy_r0->d1;
			_cy_r1[ithread] = cy_r2->d0;	_cy_i1[ithread] = cy_r2->d1;
			_cy_r2[ithread] = cy_r4->d0;	_cy_i2[ithread] = cy_r4->d1;
			_cy_r3[ithread] = cy_r6->d0;	_cy_i3[ithread] = cy_r6->d1;
			_cy_r4[ithread] = cy_r8->d0;	_cy_i4[ithread] = cy_r8->d1;
			_cy_r5[ithread] = cy_rA->d0;	_cy_i5[ithread] = cy_rA->d1;
			_cy_r6[ithread] = cy_rC->d0;	_cy_i6[ithread] = cy_rC->d1;
			_cy_r7[ithread] = cy_rE->d0;	_cy_i7[ithread] = cy_rE->d1;
			_cy_r8[ithread] = cy_i0->d0;	_cy_i8[ithread] = cy_i0->d1;
			_cy_r9[ithread] = cy_i2->d0;	_cy_i9[ithread] = cy_i2->d1;
			_cy_rA[ithread] = cy_i4->d0;	_cy_iA[ithread] = cy_i4->d1;
			_cy_rB[ithread] = cy_i6->d0;	_cy_iB[ithread] = cy_i6->d1;
			_cy_rC[ithread] = cy_i8->d0;	_cy_iC[ithread] = cy_i8->d1;
			_cy_rD[ithread] = cy_iA->d0;	_cy_iD[ithread] = cy_iA->d1;
			_cy_rE[ithread] = cy_iC->d0;	_cy_iE[ithread] = cy_iC->d1;
			_cy_rF[ithread] = cy_iE->d0;	_cy_iF[ithread] = cy_iE->d1;
			if(full_pass) maxerr = MAX(max_err->d0,max_err->d1);
		#else
			_cy_r0[ithread] = cy_r0;
			_cy_r1[ithread] = cy_r1;
			_cy_r2[ithread] = cy_r2;
			_cy_r3[ithread] = cy_r3;
			_cy_r4[ithread] = cy_r4;
			_cy_r5[ithread] = cy_r5;
			_cy_r6[ithread] = cy_r6;
			_cy_r7[ithread] = cy_r7;
			_cy_r8[ithread] = cy_r8;
			_cy_r9[ithread] = cy_r9;
			_cy_rA[ithread] = cy_rA;
			_cy_rB[ithread] = cy_rB;
			_cy_rC[ithread] = cy_rC;
			_cy_rD[ithread] = cy_rD;
			_cy_rE[ithread] = cy_rE;
			_cy_rF[ithread] = cy_rF;

			_cy_i0[ithread] = cy_i0;
			_cy_i1[ithread] = cy_i1;
			_cy_i2[ithread] = cy_i2;
			_cy_i3[ithread] = cy_i3;
			_cy_i4[ithread] = cy_i4;
			_cy_i5[ithread] = cy_i5;
			_cy_i6[ithread] = cy_i6;
			_cy_i7[ithread] = cy_i7;
			_cy_i8[ithread] = cy_i8;
			_cy_i9[ithread] = cy_i9;
			_cy_iA[ithread] = cy_iA;
			_cy_iB[ithread] = cy_iB;
			_cy_iC[ithread] = cy_iC;
			_cy_iD[ithread] = cy_iD;
			_cy_iE[ithread] = cy_iE;
			_cy_iF[ithread] = cy_iF;
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
		ASSERT(0x0 == cy16_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
	ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
	ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(0 == mlucas_nanosleep(&ns_time), "nanosleep fail!");
	}
//	printf("radix16_ditN_cy_dif1 end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		if(maxerr < tdat[ithread].maxerr) {
			maxerr = tdat[ithread].maxerr;
		}

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			_cy_r0[ithread] = tdat[ithread].cy_r0;
			_cy_r1[ithread] = tdat[ithread].cy_r1;
			_cy_r2[ithread] = tdat[ithread].cy_r2;
			_cy_r3[ithread] = tdat[ithread].cy_r3;
			_cy_r4[ithread] = tdat[ithread].cy_r4;
			_cy_r5[ithread] = tdat[ithread].cy_r5;
			_cy_r6[ithread] = tdat[ithread].cy_r6;
			_cy_r7[ithread] = tdat[ithread].cy_r7;
			_cy_r8[ithread] = tdat[ithread].cy_r8;
			_cy_r9[ithread] = tdat[ithread].cy_r9;
			_cy_rA[ithread] = tdat[ithread].cy_rA;
			_cy_rB[ithread] = tdat[ithread].cy_rB;
			_cy_rC[ithread] = tdat[ithread].cy_rC;
			_cy_rD[ithread] = tdat[ithread].cy_rD;
			_cy_rE[ithread] = tdat[ithread].cy_rE;
			_cy_rF[ithread] = tdat[ithread].cy_rF;
		}
		else
		{
			_cy_r0[ithread] = tdat[ithread].cy_r0;	_cy_i0[ithread] = tdat[ithread].cy_i0;
			_cy_r1[ithread] = tdat[ithread].cy_r1;	_cy_i1[ithread] = tdat[ithread].cy_i1;
			_cy_r2[ithread] = tdat[ithread].cy_r2;	_cy_i2[ithread] = tdat[ithread].cy_i2;
			_cy_r3[ithread] = tdat[ithread].cy_r3;	_cy_i3[ithread] = tdat[ithread].cy_i3;
			_cy_r4[ithread] = tdat[ithread].cy_r4;	_cy_i4[ithread] = tdat[ithread].cy_i4;
			_cy_r5[ithread] = tdat[ithread].cy_r5;	_cy_i5[ithread] = tdat[ithread].cy_i5;
			_cy_r6[ithread] = tdat[ithread].cy_r6;	_cy_i6[ithread] = tdat[ithread].cy_i6;
			_cy_r7[ithread] = tdat[ithread].cy_r7;	_cy_i7[ithread] = tdat[ithread].cy_i7;
			_cy_r8[ithread] = tdat[ithread].cy_r8;	_cy_i8[ithread] = tdat[ithread].cy_i8;
			_cy_r9[ithread] = tdat[ithread].cy_r9;	_cy_i9[ithread] = tdat[ithread].cy_i9;
			_cy_rA[ithread] = tdat[ithread].cy_rA;	_cy_iA[ithread] = tdat[ithread].cy_iA;
			_cy_rB[ithread] = tdat[ithread].cy_rB;	_cy_iB[ithread] = tdat[ithread].cy_iB;
			_cy_rC[ithread] = tdat[ithread].cy_rC;	_cy_iC[ithread] = tdat[ithread].cy_iC;
			_cy_rD[ithread] = tdat[ithread].cy_rD;	_cy_iD[ithread] = tdat[ithread].cy_iD;
			_cy_rE[ithread] = tdat[ithread].cy_rE;	_cy_iE[ithread] = tdat[ithread].cy_iE;
			_cy_rF[ithread] = tdat[ithread].cy_rF;	_cy_iF[ithread] = tdat[ithread].cy_iF;
		}
	}
#endif

	if(full_pass) {
	//	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

	/*   Wraparound carry cleanup loop is here:

	The cleanup carries from the end of each length-N/16 block into the begining of the next
	can all be neatly processed as follows:

	(1) Invert the radix-16 forward DIF FFT of the first block of 16 complex elements in A and unweight;
	(2) Propagate cleanup carries among the real and imaginary parts of the 16 outputs of (1);
	(3) Reweight and perform a radix-16 forward DIF FFT on the result of (2);
	(4) If any of the exit carries from (2) are nonzero, advance to the next 16 elements and repeat (1-4).
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t1 = _cy_r0[CY_THREADS - 1];
		t3 = _cy_r1[CY_THREADS - 1];
		t5 = _cy_r2[CY_THREADS - 1];
		t7 = _cy_r3[CY_THREADS - 1];
		t9 = _cy_r4[CY_THREADS - 1];
		t11= _cy_r5[CY_THREADS - 1];
		t13= _cy_r6[CY_THREADS - 1];
		t15= _cy_r7[CY_THREADS - 1];
		t17= _cy_r8[CY_THREADS - 1];
		t19= _cy_r9[CY_THREADS - 1];
		t21= _cy_rA[CY_THREADS - 1];
		t23= _cy_rB[CY_THREADS - 1];
		t25= _cy_rC[CY_THREADS - 1];
		t27= _cy_rD[CY_THREADS - 1];
		t29= _cy_rE[CY_THREADS - 1];
		t31= _cy_rF[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(CY_THREADS > 1,"");	/* Make sure loop only gets executed if multiple threads */
			_cy_r0[ithread] = _cy_r0[ithread-1];
			_cy_r1[ithread] = _cy_r1[ithread-1];
			_cy_r2[ithread] = _cy_r2[ithread-1];
			_cy_r3[ithread] = _cy_r3[ithread-1];
			_cy_r4[ithread] = _cy_r4[ithread-1];
			_cy_r5[ithread] = _cy_r5[ithread-1];
			_cy_r6[ithread] = _cy_r6[ithread-1];
			_cy_r7[ithread] = _cy_r7[ithread-1];
			_cy_r8[ithread] = _cy_r8[ithread-1];
			_cy_r9[ithread] = _cy_r9[ithread-1];
			_cy_rA[ithread] = _cy_rA[ithread-1];
			_cy_rB[ithread] = _cy_rB[ithread-1];
			_cy_rC[ithread] = _cy_rC[ithread-1];
			_cy_rD[ithread] = _cy_rD[ithread-1];
			_cy_rE[ithread] = _cy_rE[ithread-1];
			_cy_rF[ithread] = _cy_rF[ithread-1];
		}

		_cy_r0[0] =+t31;	/* ...The wraparound carry is here: */
		_cy_r1[0] = t1 ;
		_cy_r2[0] = t3 ;
		_cy_r3[0] = t5 ;
		_cy_r4[0] = t7 ;
		_cy_r5[0] = t9 ;
		_cy_r6[0] = t11;
		_cy_r7[0] = t13;
		_cy_r8[0] = t15;
		_cy_r9[0] = t17;
		_cy_rA[0] = t19;
		_cy_rB[0] = t21;
		_cy_rC[0] = t23;
		_cy_rD[0] = t25;
		_cy_rE[0] = t27;
		_cy_rF[0] = t29;
	}
	else
	{
		t1 = _cy_r0[CY_THREADS - 1];	t2 = _cy_i0[CY_THREADS - 1];
		t3 = _cy_r1[CY_THREADS - 1];	t4 = _cy_i1[CY_THREADS - 1];
		t5 = _cy_r2[CY_THREADS - 1];	t6 = _cy_i2[CY_THREADS - 1];
		t7 = _cy_r3[CY_THREADS - 1];	t8 = _cy_i3[CY_THREADS - 1];
		t9 = _cy_r4[CY_THREADS - 1];	t10= _cy_i4[CY_THREADS - 1];
		t11= _cy_r5[CY_THREADS - 1];	t12= _cy_i5[CY_THREADS - 1];
		t13= _cy_r6[CY_THREADS - 1];	t14= _cy_i6[CY_THREADS - 1];
		t15= _cy_r7[CY_THREADS - 1];	t16= _cy_i7[CY_THREADS - 1];
		t17= _cy_r8[CY_THREADS - 1];	t18= _cy_i8[CY_THREADS - 1];
		t19= _cy_r9[CY_THREADS - 1];	t20= _cy_i9[CY_THREADS - 1];
		t21= _cy_rA[CY_THREADS - 1];	t22= _cy_iA[CY_THREADS - 1];
		t23= _cy_rB[CY_THREADS - 1];	t24= _cy_iB[CY_THREADS - 1];
		t25= _cy_rC[CY_THREADS - 1];	t26= _cy_iC[CY_THREADS - 1];
		t27= _cy_rD[CY_THREADS - 1];	t28= _cy_iD[CY_THREADS - 1];
		t29= _cy_rE[CY_THREADS - 1];	t30= _cy_iE[CY_THREADS - 1];
		t31= _cy_rF[CY_THREADS - 1];	t32= _cy_iF[CY_THREADS - 1];

		// Handle valid case of high Re or Im-word < 0 and corr. cyout = 1 here, by positivizing the word:
		if(MODULUS_TYPE == MODULUS_TYPE_GENFFTMUL && (t31 != 0.0 || t32 != 0.0))
		{
			// Must use NDIVR instead of p1 here since p1 may have pads which are not applied to element-2-slots-before
			j1 = NDIVR-2;	j1 += ( (j1 >> DAT_BITS) << PAD_BITS );
			j2 = j1+RE_IM_STRIDE;
			ASSERT(t31 <= 1.0 && t32 <= 1.0, "genFFTmul expects carryouts = 0 or 1 at top!");
			// Undo the initial dif pass just for the 16 complex terms in question:
			RADIX_16_DIT(a[j1],a[j2],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
						,a[j1],a[j2],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
						,c,s)
			a[j1    ] *= radix_inv;	a[j2    ] *= radix_inv;
			a[j1+p1 ] *= radix_inv;	a[j2+p1 ] *= radix_inv;
			a[j1+p2 ] *= radix_inv;	a[j2+p2 ] *= radix_inv;
			a[j1+p3 ] *= radix_inv;	a[j2+p3 ] *= radix_inv;
			a[j1+p4 ] *= radix_inv;	a[j2+p4 ] *= radix_inv;
			a[j1+p5 ] *= radix_inv;	a[j2+p5 ] *= radix_inv;
			a[j1+p6 ] *= radix_inv;	a[j2+p6 ] *= radix_inv;
			a[j1+p7 ] *= radix_inv;	a[j2+p7 ] *= radix_inv;
			a[j1+p8 ] *= radix_inv;	a[j2+p8 ] *= radix_inv;
			a[j1+p9 ] *= radix_inv;	a[j2+p9 ] *= radix_inv;
			a[j1+p10] *= radix_inv;	a[j2+p10] *= radix_inv;
			a[j1+p11] *= radix_inv;	a[j2+p11] *= radix_inv;
			a[j1+p12] *= radix_inv;	a[j2+p12] *= radix_inv;
			a[j1+p13] *= radix_inv;	a[j2+p13] *= radix_inv;
			a[j1+p14] *= radix_inv;	a[j2+p14] *= radix_inv;
			a[j1+p15] *= radix_inv;	a[j2+p15] *= radix_inv;
			printf("CYHI.re = %10.3f, im = %10.3f, High words: Re = %10.3f, Im = %10.3f\n",t31,t32,a[j1+p15],a[j2+p15]);
			// Verify that any cyout = 1 has the corresponding high word < 0,
			// then absorb cyout back into the high word and zero the carry:
			if(t31 == 1.0) {
				ASSERT(a[j1+p15] < 0.0, "genFFTmul: Legal Re-cyout = 1 must have the corresponding high word < 0!");
				a[j1+p15] += FFT_MUL_BASE;	t31 = 0.0;
			}
			if(t32 == 1.0) {
				ASSERT(a[j2+p15] < 0.0, "genFFTmul: Legal Im-cyout = 1 must have the corresponding high word < 0!");
				a[j2+p15] += FFT_MUL_BASE;	t32 = 0.0;
			}
			// Redo the initial dif pass just for the 16 complex terms in question:
			RADIX_16_DIF(a[j1],a[j2],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
						,a[j1],a[j2],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
						,c,s)
		}

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(CY_THREADS > 1,"");	/* Make sure loop only gets executed if multiple threads */
			_cy_r0[ithread] = _cy_r0[ithread-1];	_cy_i0[ithread] = _cy_i0[ithread-1];
			_cy_r1[ithread] = _cy_r1[ithread-1];	_cy_i1[ithread] = _cy_i1[ithread-1];
			_cy_r2[ithread] = _cy_r2[ithread-1];	_cy_i2[ithread] = _cy_i2[ithread-1];
			_cy_r3[ithread] = _cy_r3[ithread-1];	_cy_i3[ithread] = _cy_i3[ithread-1];
			_cy_r4[ithread] = _cy_r4[ithread-1];	_cy_i4[ithread] = _cy_i4[ithread-1];
			_cy_r5[ithread] = _cy_r5[ithread-1];	_cy_i5[ithread] = _cy_i5[ithread-1];
			_cy_r6[ithread] = _cy_r6[ithread-1];	_cy_i6[ithread] = _cy_i6[ithread-1];
			_cy_r7[ithread] = _cy_r7[ithread-1];	_cy_i7[ithread] = _cy_i7[ithread-1];
			_cy_r8[ithread] = _cy_r8[ithread-1];	_cy_i8[ithread] = _cy_i8[ithread-1];
			_cy_r9[ithread] = _cy_r9[ithread-1];	_cy_i9[ithread] = _cy_i9[ithread-1];
			_cy_rA[ithread] = _cy_rA[ithread-1];	_cy_iA[ithread] = _cy_iA[ithread-1];
			_cy_rB[ithread] = _cy_rB[ithread-1];	_cy_iB[ithread] = _cy_iB[ithread-1];
			_cy_rC[ithread] = _cy_rC[ithread-1];	_cy_iC[ithread] = _cy_iC[ithread-1];
			_cy_rD[ithread] = _cy_rD[ithread-1];	_cy_iD[ithread] = _cy_iD[ithread-1];
			_cy_rE[ithread] = _cy_rE[ithread-1];	_cy_iE[ithread] = _cy_iE[ithread-1];
			_cy_rF[ithread] = _cy_rF[ithread-1];	_cy_iF[ithread] = _cy_iF[ithread-1];
		}

		_cy_r0[0] =-t32;	_cy_i0[0] =+t31;	/* ...The 2 Mo"bius carries are here: */
		_cy_r1[0] = t1 ;	_cy_i1[0] = t2 ;
		_cy_r2[0] = t3 ;	_cy_i2[0] = t4 ;
		_cy_r3[0] = t5 ;	_cy_i3[0] = t6 ;
		_cy_r4[0] = t7 ;	_cy_i4[0] = t8 ;
		_cy_r5[0] = t9 ;	_cy_i5[0] = t10;
		_cy_r6[0] = t11;	_cy_i6[0] = t12;
		_cy_r7[0] = t13;	_cy_i7[0] = t14;
		_cy_r8[0] = t15;	_cy_i8[0] = t16;
		_cy_r9[0] = t17;	_cy_i9[0] = t18;
		_cy_rA[0] = t19;	_cy_iA[0] = t20;
		_cy_rB[0] = t21;	_cy_iB[0] = t22;
		_cy_rC[0] = t23;	_cy_iC[0] = t24;
		_cy_rD[0] = t25;	_cy_iD[0] = t26;
		_cy_rE[0] = t27;	_cy_iE[0] = t28;
		_cy_rF[0] = t29;	_cy_iF[0] = t30;
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
			j2 = j + ( (j >> DAT_BITS) << PAD_BITS );
	#ifdef USE_FGT61
		if(!j) {
			printf("J = 0, wraparound INputs:\n");
			printf("a1p0r,a1p0i, b1p0r,b1p0i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2    ],a[j2    +1], b[j2    ],b[j2    +1]);
			printf("a1p1r,a1p1i, b1p1r,b1p1i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p1 ],a[j2+p1 +1], b[j2+p1 ],b[j2+p1 +1]);
			printf("a1p2r,a1p2i, b1p2r,b1p2i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p2 ],a[j2+p2 +1], b[j2+p2 ],b[j2+p2 +1]);
			printf("a1p3r,a1p3i, b1p3r,b1p3i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p3 ],a[j2+p3 +1], b[j2+p3 ],b[j2+p3 +1]);
			printf("a1p4r,a1p4i, b1p4r,b1p4i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p4 ],a[j2+p4 +1], b[j2+p4 ],b[j2+p4 +1]);
			printf("a1p5r,a1p5i, b1p5r,b1p5i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p5 ],a[j2+p5 +1], b[j2+p5 ],b[j2+p5 +1]);
			printf("a1p6r,a1p6i, b1p6r,b1p6i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p6 ],a[j2+p6 +1], b[j2+p6 ],b[j2+p6 +1]);
			printf("a1p7r,a1p7i, b1p7r,b1p7i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p7 ],a[j2+p7 +1], b[j2+p7 ],b[j2+p7 +1]);
			printf("a1p8r,a1p8i, b1p8r,b1p8i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p8 ],a[j2+p8 +1], b[j2+p8 ],b[j2+p8 +1]);
			printf("a1p9r,a1p9i, b1p9r,b1p9i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p9 ],a[j2+p9 +1], b[j2+p9 ],b[j2+p9 +1]);
			printf("a1pAr,a1pAi, b1pAr,b1pAi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p10],a[j2+p10+1], b[j2+p10],b[j2+p10+1]);
			printf("a1pBr,a1pBi, b1pBr,b1pBi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p11],a[j2+p11+1], b[j2+p11],b[j2+p11+1]);
			printf("a1pCr,a1pCi, b1pCr,b1pCi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p12],a[j2+p12+1], b[j2+p12],b[j2+p12+1]);
			printf("a1pDr,a1pDi, b1pDr,b1pDi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p13],a[j2+p13+1], b[j2+p13],b[j2+p13+1]);
			printf("a1pEr,a1pEi, b1pEr,b1pEi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p14],a[j2+p14+1], b[j2+p14],b[j2+p14+1]);
			printf("a1pFr,a1pFi, b1pFr,b1pFi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p15],a[j2+p15+1], b[j2+p15],b[j2+p15+1]);
		}
	#endif
			a[j2    ] *= radix_inv;
			a[j2+p1 ] *= radix_inv;
			a[j2+p2 ] *= radix_inv;
			a[j2+p3 ] *= radix_inv;
			a[j2+p4 ] *= radix_inv;
			a[j2+p5 ] *= radix_inv;
			a[j2+p6 ] *= radix_inv;
			a[j2+p7 ] *= radix_inv;
			a[j2+p8 ] *= radix_inv;
			a[j2+p9 ] *= radix_inv;
			a[j2+p10] *= radix_inv;
			a[j2+p11] *= radix_inv;
			a[j2+p12] *= radix_inv;
			a[j2+p13] *= radix_inv;
			a[j2+p14] *= radix_inv;
			a[j2+p15] *= radix_inv;
		#ifdef USE_FGT61
			// Modmul *= 1/16 ==> mul_pow2_modq(*, 61-lg(16)) = mul_pow2_modq(*, 57):
			b[j2    ] = mul_pow2_modq( b[j2    ], 57);
			b[j2+p1 ] = mul_pow2_modq( b[j2+p1 ], 57);
			b[j2+p2 ] = mul_pow2_modq( b[j2+p2 ], 57);
			b[j2+p3 ] = mul_pow2_modq( b[j2+p3 ], 57);
			b[j2+p4 ] = mul_pow2_modq( b[j2+p4 ], 57);
			b[j2+p5 ] = mul_pow2_modq( b[j2+p5 ], 57);
			b[j2+p6 ] = mul_pow2_modq( b[j2+p6 ], 57);
			b[j2+p7 ] = mul_pow2_modq( b[j2+p7 ], 57);
			b[j2+p8 ] = mul_pow2_modq( b[j2+p8 ], 57);
			b[j2+p9 ] = mul_pow2_modq( b[j2+p9 ], 57);
			b[j2+p10] = mul_pow2_modq( b[j2+p10], 57);
			b[j2+p11] = mul_pow2_modq( b[j2+p11], 57);
			b[j2+p12] = mul_pow2_modq( b[j2+p12], 57);
			b[j2+p13] = mul_pow2_modq( b[j2+p13], 57);
			b[j2+p14] = mul_pow2_modq( b[j2+p14], 57);
			b[j2+p15] = mul_pow2_modq( b[j2+p15], 57);
		if(j==1) {
			printf("J = 0, wraparound OUTputs:\n");
			printf("a1p0r,a1p0i, b1p0r,b1p0i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2    -1],a[j2    ], b[j2    -1],b[j2    ]);
			printf("a1p1r,a1p1i, b1p1r,b1p1i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p1 -1],a[j2+p1 ], b[j2+p1 -1],b[j2+p1 ]);
			printf("a1p2r,a1p2i, b1p2r,b1p2i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p2 -1],a[j2+p2 ], b[j2+p2 -1],b[j2+p2 ]);
			printf("a1p3r,a1p3i, b1p3r,b1p3i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p3 -1],a[j2+p3 ], b[j2+p3 -1],b[j2+p3 ]);
			printf("a1p4r,a1p4i, b1p4r,b1p4i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p4 -1],a[j2+p4 ], b[j2+p4 -1],b[j2+p4 ]);
			printf("a1p5r,a1p5i, b1p5r,b1p5i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p5 -1],a[j2+p5 ], b[j2+p5 -1],b[j2+p5 ]);
			printf("a1p6r,a1p6i, b1p6r,b1p6i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p6 -1],a[j2+p6 ], b[j2+p6 -1],b[j2+p6 ]);
			printf("a1p7r,a1p7i, b1p7r,b1p7i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p7 -1],a[j2+p7 ], b[j2+p7 -1],b[j2+p7 ]);
			printf("a1p8r,a1p8i, b1p8r,b1p8i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p8 -1],a[j2+p8 ], b[j2+p8 -1],b[j2+p8 ]);
			printf("a1p9r,a1p9i, b1p9r,b1p9i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p9 -1],a[j2+p9 ], b[j2+p9 -1],b[j2+p9 ]);
			printf("a1pAr,a1pAi, b1pAr,b1pAi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p10-1],a[j2+p10], b[j2+p10-1],b[j2+p10]);
			printf("a1pBr,a1pBi, b1pBr,b1pBi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p11-1],a[j2+p11], b[j2+p11-1],b[j2+p11]);
			printf("a1pCr,a1pCi, b1pCr,b1pCi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p12-1],a[j2+p12], b[j2+p12-1],b[j2+p12]);
			printf("a1pDr,a1pDi, b1pDr,b1pDi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p13-1],a[j2+p13], b[j2+p13-1],b[j2+p13]);
			printf("a1pEr,a1pEi, b1pEr,b1pEi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p14-1],a[j2+p14], b[j2+p14-1],b[j2+p14]);
			printf("a1pFr,a1pFi, b1pFr,b1pFi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j2+p15-1],a[j2+p15], b[j2+p15-1],b[j2+p15]);
		}
		#endif
		}
	}
}	/* endfor(outer) */

	t1 = 0;
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		t1 += fabs(_cy_r0[0])+fabs(_cy_r1[0])+fabs(_cy_r2[0])+fabs(_cy_r3[0])+fabs(_cy_r4[0])+fabs(_cy_r5[0])+fabs(_cy_r6[0])+fabs(_cy_r7[0])+fabs(_cy_r8[0])+fabs(_cy_r9[0])+fabs(_cy_rA[0])+fabs(_cy_rB[0])+fabs(_cy_rC[0])+fabs(_cy_rD[0])+fabs(_cy_rE[0])+fabs(_cy_rF[0]);
		t1 += fabs(_cy_i0[0])+fabs(_cy_i1[0])+fabs(_cy_i2[0])+fabs(_cy_i3[0])+fabs(_cy_i4[0])+fabs(_cy_i5[0])+fabs(_cy_i6[0])+fabs(_cy_i7[0])+fabs(_cy_i8[0])+fabs(_cy_i9[0])+fabs(_cy_iA[0])+fabs(_cy_iB[0])+fabs(_cy_iC[0])+fabs(_cy_iD[0])+fabs(_cy_iE[0])+fabs(_cy_iF[0]);

		*fracmax = maxerr;
	}
//fprintf(stderr, "radix16_carry: A[0-3] = %20.5f %20.5f %20.5f %20.5f\n",a[0],a[1],a[2],a[3]);
	if(t1 != 0.0)
	{
		sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in %s - input wordsize may be too small.\n",iter,func);
		mlucas_fprint(cbuf,INTERACT);
		err = ERR_CARRY;
		return(err);
	}

	return(0);
}

/**************/

#ifdef USE_FGT61
void radix16_dif_pass1	(double a[], uint64 b[], int n)
#else
void radix16_dif_pass1	(double a[],             int n)
#endif
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-16 complex DIF FFT pass on the data in the length-N real vector A.
!
!   The data are stored in a 1-D zero-offset array, with 2^PAD_BITS 8-byte padding elements inserted
!   between every block of 2^DAT_BITS contiguous data. The array padding is to prevent data being accessed
!   in strides that are large powers of two and thus to minimize cache thrashing
!   (in cache-based microprocessor architectures) or bank conflicts (in supercomputers.)
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2;
	static int n16,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15, first_entry=TRUE;
#ifdef USE_FGT61
	const uint64 q  = 0x1FFFFFFFFFFFFFFFull, q2=q+q, q3=q2+q, q4=q2+q2, q5=q4+q;	// q = 2^61 - 1, and needed small multiples
	// primitive 16th root of unity, scaled by *8:
	const uint64 cm = 1693317751237720973ull<<3, sm = 2283815672160731785ull<<3;
#endif
	static double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;	/* exp[i*(twopi/16)] */
  #if !USE_SCALAR_DFT_MACRO
	double rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;
  #endif
  #ifdef USE_FGT61
	int jt,jp;
	uint64 rm,im;
	uint64 m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32;
  #endif

	if(!first_entry && (n >> 4) != n16)	/* New runlength?	*/
	{
	  first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
	  first_entry=FALSE;
	  n16=n/16;

	  p1 = n16 + ( (n16 >> DAT_BITS) << PAD_BITS );
	  p2 = p1 +p1;
	  p3 = p2 +p1;
	  p4 = p3 +p1;
	  p5 = p4 +p1;
	  p6 = p5 +p1;
	  p7 = p6 +p1;
	  p8 = p7 +p1;
	  p9 = p8 +p1;
	  p10= p9 +p1;
	  p11= p10+p1;
	  p12= p11+p1;
	  p13= p12+p1;
	  p14= p13+p1;
	  p15= p14+p1;
	}

/*...The radix-16 pass is here.	*/

	for(j = 0; j < n16; j += 2)
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

		RADIX_16_DIF(a[j1],a[j2],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
					,a[j1],a[j2],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
					,c,s)
		// Test FMA-based DIF macro, with twiddles set = unity:
//		RADIX_16_DIF_FMA(a[j1],a[j2],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
//						,a[j1],a[j2],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
//						,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,c,s)

	#else

	  #ifdef USE_FGT61

	/*...Block 1: */									// All inputs in 0,b - since b = q+7, i.e. exceeds q by < 1 part
		jt = j1;		jp = j2;						// in 2^58, treat as 0,q if it saves a nontrivial #qreduce calls
		t1 =a[jt	];	t2 =a[jp	];					m1 =b[jt	];	m2 =b[jp	];
		rt =a[jt+p8 ];	it =a[jp+p8 ];					rm =b[jt+p8 ];	im =b[jp+p8 ];
		t3 =t1 -rt;		t1 =t1 +rt;						m3 =m1 -rm;		m1 =m1 +rm;			// -:-b,b
		t4 =t2 -it;		t2 =t2 +it;						m4 =m2 -im;		m2 =m2 +im;			// +:0,2b

		t5 =a[jt+p4 ];	t6 =a[jp+p4 ];					m5 =b[jt+p4 ];	m6 =b[jp+p4 ];
		rt =a[jt+p12];	it =a[jp+p12];					rm =b[jt+p12];	im =b[jp+p12];
		t7 =t5 -rt;		t5 =t5 +rt;						m7 =m5 -rm;		m5 =m5 +rm;			// -:-b,b
		t8 =t6 -it;		t6 =t6 +it;						m8 =m6 -im;		m6 =m6 +im;			// +:0,2b
														rm =m5;			im =m6		;		// 0,2b
		rt =t5;	t5 =t1 -rt;			t1 =t1 +rt;			m5 =m1 -rm;		m6 =m2 -im	;		// m5,6: -2b,2b
		it =t6;	t6 =t2 -it;			t2 =t2 +it;			m1 =m1 +rm;		m2 =m2 +im	;		// m1,2: 0,4b
														rm =m7;			im =m8		;
		rt =t7;	t7 =t3 +t8;			t3 =t3 -t8;			m7 =m3 +im;		m8 =m4 -rm	;		// m3,4,7,8 all in -2b,2b
				t8 =t4 -rt;			t4 =t4 +rt;			m3 =m3 -im;		m4 =m4 +rm	;		// (Blocks 2-4 similar)

	/*...Block 2: */
		jt = j1 + p2;	jp = j2 + p2;
		t9 =a[jt    ];	t10=a[jp    ];					m9 =b[jt    ];	m10=b[jp    ];
		rt =a[jt+p8 ];	it =a[jp+p8 ];					rm =b[jt+p8 ];	im =b[jp+p8 ];
		t11=t9 -rt;		t9 =t9 +rt;						m11=m9 -rm;		m9 =m9 +rm;
		t12=t10-it;		t10=t10+it;						m12=m10-im;		m10=m10+im;

		t13=a[jt+p4 ];	t14=a[jp+p4 ];					m13=b[jt+p4 ];	m14=b[jp+p4 ];
		rt =a[jt+p12];	it =a[jp+p12];					rm =b[jt+p12];	im =b[jp+p12];
		t15=t13-rt;		t13=t13+rt;						m15=m13-rm;		m13=m13+rm;
		t16=t14-it;		t14=t14+it;						m16=m14-im;		m14=m14+im;
														rm =m13;		im =m14		;
		rt =t13;	t13=t9 -rt;		t9 =t9 +rt;			m13=m9 -rm;		m14=m10-im	;
		it =t14;	t14=t10-it;		t10=t10+it;			m9 =m9 +rm;		m10=m10+im	;
														rm =m15;		im =m16		;
		rt =t15;	t15=t11+t16;	t11=t11-t16;		m15=m11+im;		m16=m12-rm	;
					t16=t12-rt;		t12=t12+rt;			m11=m11-im;		m12=m12+rm	;

	/*...Block 3: */
		jt = j1 + p1;	jp = j2 + p1;
		t17=a[jt    ];	t18=a[jp    ];					m17=b[jt    ];	m18=b[jp    ];
		rt =a[jt+p8 ];	it =a[jp+p8 ];					rm =b[jt+p8 ];	im =b[jp+p8 ];
		t19=t17-rt;		t17=t17+rt;						m19=m17-rm;		m17=m17+rm;
		t20=t18-it;		t18=t18+it;						m20=m18-im;		m18=m18+im;

		t21=a[jt+p4 ];	t22=a[jp+p4 ];					m21=b[jt+p4 ];	m22=b[jp+p4 ];
		rt =a[jt+p12];	it =a[jp+p12];					rm =b[jt+p12];	im =b[jp+p12];
		t23=t21-rt;		t21=t21+rt;						m23=m21-rm;		m21=m21+rm;
		t24=t22-it;		t22=t22+it;						m24=m22-im;		m22=m22+im;
														rm =m21;					im =m22		;
		rt =t21;	t21=t17-rt;		t17=t17+rt;			m21=m17-rm;					m22=m18-im	;
		it =t22;	t22=t18-it;		t18=t18+it;			m17=m17+rm;					m18=m18+im	;
														rm =m23;					im =m24		;
		rt =t23;	t23=t19+t24;	t19=t19-t24;		m23=qreduce(m19+im+q4);		m24=qreduce(m20-rm+q4);	// all in -2b,2b
					t24=t20-rt;		t20=t20+rt;			m19=qreduce(m19-im+q4);		m20=qreduce(m20+rm+q4);	// prior to reduction
														// m19,20,23,24 are needed for CMUL, so reduce.
	/*...Block 4: */
		jt = j1 + p3;	jp = j2 + p3;
		t25=a[jt    ];	t26=a[jp    ];					m25=b[jt    ];	m26=b[jp    ];
		rt =a[jt+p8 ];	it =a[jp+p8 ];					rm =b[jt+p8 ];	im =b[jp+p8 ];
		t27=t25-rt;		t25=t25+rt;						m27=m25-rm;		m25=m25+rm;
		t28=t26-it;		t26=t26+it;						m28=m26-im;		m26=m26+im;

		t29=a[jt+p4 ];	t30=a[jp+p4 ];					m29=b[jt+p4 ];	m30=b[jp+p4 ];
		rt =a[jt+p12];	it =a[jp+p12];					rm =b[jt+p12];	im =b[jp+p12];
		t31=t29-rt;		t29=t29+rt;						m31=m29-rm;		m29=m29+rm;
		t32=t30-it;		t30=t30+it;						m32=m30-im;		m30=m30+im;
														rm =m29;					im =m30		;
		rt =t29;	t29=t25-rt;		t25=t25+rt;			m29=m25-rm;					m30=m26-im	;
		it =t30;	t30=t26-it;		t26=t26+it;			m25=m25+rm;					m26=m26+im	;
														rm =m31;					im =m32		;
		rt =t31;	t31=t27+t32;	t27=t27-t32;		m31=qreduce(m27+im+q4);		m32=qreduce(m28-rm+q4);	// all in -2b,2b
					t32=t28-rt;		t28=t28+rt;			m27=qreduce(m27-im+q4);		m28=qreduce(m28+rm+q4);	// prior to reduction
														// m27,28,31,32 are needed for CMUL, so reduce.
	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	*/
													/*===============
													Input bounds on modular terms:
														-2b,2b: m3-8, 11-16, 21,22, 29,30
														  0,4b: m1,2, 9 ,10, 17,18, 25,26
													These are all reduced (in 0,b) because they are inputs to CMUL:
																m19,20,23,24,27,28,31,32
													===============*/
	/*...Block 1: t1,9,17,25 */
		jt = j1;		jp = j2;
		/* Debug: check for overflow of + terms: */	ASSERT(m1+m9 >= m1 && m2+m10 >= m2,"Overflow of [0,8b] term!");
		rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;				rm =m9;	m9 =qreduce(m1 -rm+q4);	m1 =qreduce(m1 +rm   );	//  1, 2 in   0,8b -> 0,b
		it =t10;t10=t2 -it;	t2 =t2 +it;				im =m10;m10=qreduce(m2 -im+q4);	m2 =qreduce(m2 +im+q4);	//  9,10 in -4b,4b -> 0,b

		rt =t25;t25=t17-rt;	t17=t17+rt;				rm =m25;m25=qreduce(m17-rm+q4);	m17=qreduce(m17+rm   );	// 17,18 in   0,4b -> 0,b
		it =t26;t26=t18-it;	t18=t18+it;				im =m26;m26=qreduce(m18-im+q4);	m18=qreduce(m18+im+q2);	// 25,26 in -2b,2b -> 0,b
													// + terms in 0,2b; - tems in -b,b:
		a[jt    ]= t1+t17;	a[jp    ]= t2+t18;		b[jt    ]=qreduce( m1+m17   );	b[jp    ]=qreduce( m2+m18   );
		a[jt+p1 ]= t1-t17;	a[jp+p1 ]= t2-t18;		b[jt+p1 ]=qreduce( m1-m17+q2);	b[jp+p1 ]=qreduce( m2-m18+q2);
		// mpy by E^4=i is inlined here:
		a[jt+p2 ]=t9 -t26;	a[jp+p2 ]=t10+t25;		b[jt+p2 ]=qreduce(m9 -m26+q2);	b[jp+p2 ]=qreduce(m10+m25+q2);
		a[jt+p3 ]=t9 +t26;	a[jp+p3 ]=t10-t25;		b[jt+p3 ]=qreduce(m9 +m26+q2);	b[jp+p3 ]=qreduce(m10-m25+q2);

	/*...Block 3: t5,13,21,29 */
		jt = j1 + p4;		jp = j2 + p4;
		rt =t13;t13=t5 +t14;t5 =t5 -t14;			rm =m13;m13=qreduce(m5 +m14+q4);m5 =qreduce(m5 -m14+q4);	// m5,6,13,14 in
				t14=t6 -rt;	t6 =t6 +rt;						m14=qreduce(m6 -rm +q4);m6 =qreduce(m6 +rm +q4);	// -4b,4b -> 0,b
	// twiddle mpy by E^2:							// All 4 +- sums in -4b,4b; add q4 to ensure (more orless) MUL inputs >= 0:
		rt =(t21-t22)*ISRT2;t22=(t21+t22)*ISRT2;			rm = mul_i2(m21-m22+q4);m22= mul_i2(m21+m22+q4);	// All 4 outs
t21=rt;	rt =(t30+t29)*ISRT2;it =(t30-t29)*ISRT2;	m21=rm;	rm = mul_i2(m30+m29+q4);im = mul_i2(m30-m29+q4);	// in [0,b30]
		t29=t21+rt;			t21=t21-rt;						m29=m21+rm;			m21=m21-rm;		// m21,22 in [-b30,b30]
		t30=t22+it;			t22=t22-it;						m30=m22+im;			m22=m22-im;		// m29,30 in [0,2*b30]

		a[jt    ]= t5+t21;	a[jp    ]= t6+t22;		b[jt    ]=qreduce( m5+m21+q4);	b[jp    ]=qreduce( m6+m22+q4);	// + and - in [0,b] + [-b30,b30] = [-b30,b+b30]
		a[jt+p1 ]= t5-t21;	a[jp+p1 ]= t6-t22;		b[jt+p1 ]=qreduce( m5-m21+q4);	b[jp+p1 ]=qreduce( m6-m22+q4);
		// mpy by E^4=i is inlined here:
		a[jt+p2 ]=t13-t30;	a[jp+p2 ]=t14+t29;		b[jt+p2 ]=qreduce(m13-m30+q3);	b[jp+p2 ]=qreduce(m14+m29   );	// + in [0,b] + [0,2*b30] = [0,b+2*b30)]
		a[jt+p3 ]=t13+t30;	a[jp+p3 ]=t14-t29;		b[jt+p3 ]=qreduce(m13+m30   );	b[jp+p3 ]=qreduce(m14-m29+q3);	// - in [0,b] - [0,2*b30] = [-2*b30,b]

	/*...Block 2: t3,11,19,27 */
		jt = j1 + p8;		jp = j2 + p8;
	// twiddle mpy by E^2							// m3,4,11,12 all in -2b,2b -> m11+-m12 in -4b,4b -> rm,im in 0,b30:
		rt =(t11-t12)*ISRT2;it =(t11+t12)*ISRT2;	rm =mul_i2(m11-m12+q4);	im =mul_i2(m11+m12+q4);
		t11=t3 -rt;			t3 =t3 +rt;				m11=m3 -rm;				m12=m4 -im;	// + in [-2b,2b] + [0,b30] = -2b,2b+b30
		t12=t4 -it;			t4 =t4 +it;				m3 =m3 +rm;				m4 =m4 +im;	// - in [-2b,2b] - [0,b30] = -2b-b30,2b

		rt =t19*c - t20*s;	t20=t20*c + t19*s;		cmul_modq8(m19,m20, cm,sm, &m19,&m20);	// 0,4q
t19=rt;	rt =t27*s - t28*c;	it =t28*s + t27*c;		cmul_modq8(m27,m28, sm,cm,  &rm, &im);	// 0,4q
		t27=t19-rt;			t19=t19+rt;				m27=qreduce(m19-rm+q4);	m19=qreduce(m19+rm);	// +:   0,8q ==> 0,b
		t28=t20-it;			t20=t20+it;				m28=qreduce(m20-im+q4);	m20=qreduce(m20+im);	// -: -4q,4q ==> 0,b
													// m3,4 in -2b,2b+b30;  m19,20 in 0,b:
		a[jt    ]= t3+t19;	a[jp    ]= t4+t20;		b[jt    ]=qreduce(m3+m19+q3);	b[jp    ]=qreduce(m4+m20+q3);	// + in [-2b,2b+b30] + [0,b] = -2b,3b+b30
		a[jt+p1 ]= t3-t19;	a[jp+p1 ]= t4-t20;		b[jt+p1 ]=qreduce(m3-m19+q4);	b[jp+p1 ]=qreduce(m4-m20+q4);	// - in [-2b,2b+b30] - [0,b] = -3b,2b+b30
		// mpy by E^4=i is inlined here:			// m11,12 in -2b-b30,2b;  m27,28 in 0,b:
		a[jt+p2 ]=t11-t28;	a[jp+p2 ]=t12+t27;		b[jt+p2 ]=qreduce(m11-m28+q5);	b[jp+p2 ]=qreduce(m12+m27+q4);	// + in [-2b-b30,2b] + [0,b] = -2b-b30,3b
		a[jt+p3 ]=t11+t28;	a[jp+p3 ]=t12-t27;		b[jt+p3 ]=qreduce(m11+m28+q4);	b[jp+p3 ]=qreduce(m12-m27+q5);	// - in [-2b-b30,2b] - [0,b] = -3b-b30,2b

	/*...Block 4: t7,15,23,31 */
		jt = j1 + p12;		jp = j2 + p12;
		/* twiddle mpy by -E^6 is here... */		// m7,8,15,16 all in -2b,2b -> m16+-m15 in -4b,4b -> rm,im in 0,b30:
		rt =(t16+t15)*ISRT2;it =(t16-t15)*ISRT2;	rm =mul_i2(m16+m15+q4);	im =mul_i2(m16-m15+q4);
		t15=t7 +rt;			t7 =t7 -rt;				m15=m7 +rm;				m16=m8 +im;	// + in [-2b,2b] + [0,b30] = -2b,2b+b30
		t16=t8 +it;			t8 =t8 -it;				m7 =m7 -rm;				m8 =m8 -im;	// - in [-2b,2b] - [0,b30] = -2b-b30,2b
		// E^3 = s,c; E^9 = -c,-s = -E^1, do latter via twiddle-mul by E^1, then flip signs on ensuing +- [re,im]:
		rt =t23*s - t24*c;	t24=t24*s + t23*c;		cmul_modq8(m23,m24, sm,cm, &m23,&m24);	// 0,4q
t23=rt;	rt =t31*c - t32*s;	it =t32*c + t31*s;		cmul_modq8(m31,m32, cm,sm,  &rm, &im);	// 0,4q
		t31=t23+rt;			t23=t23-rt;				m31=qreduce(m23+rm);	m23=qreduce(m23-rm+q4);	// +:   0,8q ==> 0,b
		t32=t24+it;			t24=t24-it;				m32=qreduce(m24+im);	m24=qreduce(m24-im+q4);	// -: -4q,4q ==> 0,b
													// m7,8 in -2b-b30,2b;  m23,24 in 0,b:
		a[jt    ]= t7+t23;	a[jp    ]= t8+t24;		b[jt    ]=qreduce( m7+m23+q4);	b[jp    ]=qreduce( m8+m24+q4);	// + in [-2b-b30,2b] + [0,b] = -2b-b30,3b
		a[jt+p1 ]= t7-t23;	a[jp+p1 ]= t8-t24;		b[jt+p1 ]=qreduce( m7-m23+q5);	b[jp+p1 ]=qreduce( m8-m24+q5);	// - in [-2b-b30,2b] - [0,b] = -3b-b30,2b
		// mpy by E^4=i is inlined here:			// m15,16 in -2b,2b+b30;  m31,32 in 0,b:
		a[jt+p2 ]=t15-t32;	a[jp+p2 ]=t16+t31;		b[jt+p2 ]=qreduce(m15-m32+q4);	b[jp+p2 ]=qreduce(m16+m31+q3);	// + in [-2b,2b+b30] + [0,b] = -2b,3b+b30
		a[jt+p3 ]=t15+t32;	a[jp+p3 ]=t16-t31;		b[jt+p3 ]=qreduce(m15+m32+q3);	b[jp+p3 ]=qreduce(m16-m31+q4);	// - in [-2b,2b+b30] - [0,b] = -3b,2b+b30

			/**********************************************/
	  #else	// USE_FGT61 = False; Basic scalar-double mode:
			/**********************************************/

		/*...Block 1:	*/
		t1 =a[j1    ];	t2 =a[j2    ];
		rt =a[j1+p8 ];	it =a[j2+p8 ];
		t3 =t1 -rt;  	t1 =t1 +rt;
		t4 =t2 -it;		t2 =t2 +it;

		t5 =a[j1+p4 ];	t6 =a[j2+p4 ];
		rt =a[j1+p12];	it =a[j2+p12];
		t7 =t5 -rt;  	t5 =t5 +rt;
		t8 =t6 -it;  	t6 =t6 +it;

		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;

		rt =t7;	t7 =t3 +t8;	t3 =t3 -t8;
				t8 =t4 -rt;	t4 =t4 +rt;

		/*...Block 2:	*/
		t9 =a[j1+p2 ];	t10=a[j2+p2 ];
		rt =a[j1+p10];	it =a[j2+p10];
		t11=t9 -rt;  	t9 =t9 +rt;
		t12=t10-it;		t10=t10+it;

		t13=a[j1+p6 ];	t14=a[j2+p6 ];
		rt =a[j1+p14];	it =a[j2+p14];
		t15=t13-rt;  	t13=t13+rt;
		t16=t14-it;		t14=t14+it;

		rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
		it =t14;	t14=t10-it;	t10=t10+it;

		rt =t15;	t15=t11+t16;t11=t11-t16;
					t16=t12-rt;	t12=t12+rt;

		/*...Block 3:	*/
		t17=a[j1+p1 ];	t18=a[j2+p1 ];
		rt =a[j1+p9 ];	it =a[j2+p9 ];
		t19=t17-rt;  	t17=t17+rt;
		t20=t18-it;		t18=t18+it;

		t21=a[j1+p5 ];	t22=a[j2+p5 ];
		rt =a[j1+p13];	it =a[j2+p13];
		t23=t21-rt;  	t21=t21+rt;
		t24=t22-it;		t22=t22+it;

		rt =t21;	t21=t17-rt;	t17=t17+rt;
		it =t22;	t22=t18-it;	t18=t18+it;

		rt =t23;	t23=t19+t24;t19=t19-t24;
					t24=t20-rt;	t20=t20+rt;

		/*...Block 4:	*/
		t25=a[j1+p3 ];	t26=a[j2+p3 ];
		rt =a[j1+p11];	it =a[j2+p11];
		t27=t25-rt;  	t25=t25+rt;
		t28=t26-it;		t26=t26+it;

		t29=a[j1+p7 ];	t30=a[j2+p7 ];
		rt =a[j1+p15];	it =a[j2+p15];
		t31=t29-rt;  	t29=t29+rt;
		t32=t30-it;		t30=t30+it;

		rt =t29;	t29=t25-rt;	t25=t25+rt;
		it =t30;	t30=t26-it;	t26=t26+it;

		rt =t31;	t31=t27+t32;t27=t27-t32;
					t32=t28-rt;	t28=t28+rt;

		/*...and now do four more radix-4 transforms, including the internal twiddle factors:
			1, exp(i* 1*twopi/16) =       ( c, s), exp(i* 2*twopi/16) = ISRT2*( 1, 1), exp(i* 3*twopi/16) =       ( s, c) (for inputs to transform block 2)
			1, exp(i* 2*twopi/16) = ISRT2*( 1, 1), exp(i* 4*twopi/16) =       ( 0, 1), exp(i* 6*twopi/16) = ISRT2*(-1, 1) (for inputs to transform block 3)
			1, exp(i* 3*twopi/16) =       ( s, c), exp(i* 6*twopi/16) = ISRT2*(-1, 1), exp(i* 9*twopi/16) =       (-c,-s) (for inputs to transform block 4).
			(This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
			 I.e. do similar as above, except inputs a(j1  +p0:15:1) are replaced by t0:30:2,
													 a(j2+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
			 and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.	*/

		/*...Block 1: t1,9,17,25	*/
	    rt =t9 ;	t9 =t1 -rt;	t1 =t1 +rt;
	    it =t10;	t10=t2 -it;	t2 =t2 +it;

	    rt =t25;	t25=t17-rt;	t17=t17+rt;
	    it =t26;	t26=t18-it;	t18=t18+it;

	    a[j1    ]=t1 +t17;	a[j2    ]=t2 +t18;
	    a[j1+p1 ]=t1 -t17;	a[j2+p1 ]=t2 -t18;

		a[j1+p2 ]=t9 -t26;	a[j2+p2 ]=t10+t25;	/* mpy by E^4=i is inlined here...	*/
		a[j1+p3 ]=t9 +t26;	a[j2+p3 ]=t10-t25;

		/*...Block 3: t5,13,21,29	*/
		rt =t13;	t13=t5 +t14;	t5 =t5 -t14;		/* twiddle mpy by E^4 = I	*/
			t14=t6 -rt;	t6 =t6 +rt;

		rt =(t21-t22)*ISRT2;t22=(t21+t22)*ISRT2;	t21=rt;	/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t21,22,29,30	*/
		rt =(t30+t29)*ISRT2;it =(t30-t29)*ISRT2;		/* twiddle mpy by -E^6 is here...	*/
		t29=t21+rt;		t21=t21-rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here.	*/
		t30=t22+it;		t22=t22-it;

		a[j1+p4 ]=t5+t21;	a[j2+p4 ]=t6+t22;
		a[j1+p5 ]=t5-t21;	a[j2+p5 ]=t6-t22;

		a[j1+p6 ]=t13-t30;	a[j2+p6 ]=t14+t29;	/* mpy by E^4=i is inlined here...	*/
		a[j1+p7 ]=t13+t30;	a[j2+p7 ]=t14-t29;

		/*...Block 2: t3,11,19,27	*/
		rt =(t11-t12)*ISRT2;it =(t11+t12)*ISRT2;		/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t11,12	*/
		t11=t3 -rt;		t3 =t3 +rt;
		t12=t4 -it;		t4 =t4 +it;

		rt =t19*c - t20*s;	t20=t20*c + t19*s;	t19=rt;	/* twiddle mpy by E^1	*/
		rt =t27*s - t28*c;	it =t28*s + t27*c;		/* twiddle mpy by E^3	*/
		t27=t19-rt;		t19=t19+rt;
		t28=t20-it;		t20=t20+it;

		a[j1+p8 ]=t3+t19;	a[j2+p8 ]=t4+t20;
		a[j1+p9 ]=t3-t19;	a[j2+p9 ]=t4-t20;

		a[j1+p10]=t11-t28;	a[j2+p10]=t12+t27;	/* mpy by E^4=i is inlined here...	*/
		a[j1+p11]=t11+t28;	a[j2+p11]=t12-t27;

		/*...Block 4: t7,15,23,31	*/
		rt =(t16+t15)*ISRT2;it =(t16-t15)*ISRT2;		/* twiddle mpy by -E^6 is here...	*/
		t15=t7 +rt;		t7 =t7 -rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here.	*/
		t16=t8 +it;		t8 =t8 -it;

		rt =t23*s - t24*c;	t24=t24*s + t23*c;	t23=rt;	/* twiddle mpy by E^3	*/
		rt =t31*c - t32*s;	it =t32*c + t31*s;		/* twiddle mpy by E^1 = -E^9...	*/
		t31=t23+rt;		t23=t23-rt;			/* ...and get E^9 by flipping signs here.	*/
		t32=t24+it;		t24=t24-it;

		/* Note: t23+rt = t23*(s+1)	*/

		a[j1+p12]=t7+t23;	a[j2+p12]=t8+t24;
		a[j1+p13]=t7-t23;	a[j2+p13]=t8-t24;

		a[j1+p14]=t15-t32;	a[j2+p14]=t16+t31;	/* mpy by E^4=i is inlined here...	*/
		a[j1+p15]=t15+t32;	a[j2+p15]=t16-t31;

	  #endif	// USE_FGT61 ?

	#endif
	}
}

/**************/

#ifdef USE_FGT61
void radix16_dit_pass1	(double a[], uint64 b[], int n)
#else
void radix16_dit_pass1	(double a[],             int n)
#endif
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-16 complex inverse DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details.
*/
	int j,j1,j2;
	static int n16,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15, first_entry=TRUE;
#ifdef USE_FGT61
	const uint64 q  = 0x1FFFFFFFFFFFFFFFull, q2=q+q, q3=q2+q, q4=q2+q2, q5=q4+q, q8=q4+q4;	// q = 2^61 - 1, and needed small multiples
	// primitive 16th root of unity, scaled by *8:
	const uint64 cm = 1693317751237720973ull<<3, sm = 2283815672160731785ull<<3;
#endif
	static double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;	/* exp[i*(twopi/16)]	*/
  #if !USE_SCALAR_DFT_MACRO
	double rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;
  #endif
  #ifdef USE_FGT61
	int jt,jp;
	uint64 rm,im;
	uint64 m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32;
  #endif

	if(!first_entry && (n >> 4) != n16)	/* New runlength?	*/
	{
	  first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
	  first_entry=FALSE;
	  n16=n/16;

	  p1 = n16 + ( (n16 >> DAT_BITS) << PAD_BITS );
	  p2 = p1 +p1;
	  p3 = p2 +p1;
	  p4 = p3 +p1;
	  p5 = p4 +p1;
	  p6 = p5 +p1;
	  p7 = p6 +p1;
	  p8 = p7 +p1;
	  p9 = p8 +p1;
	  p10= p9 +p1;
	  p11= p10+p1;
	  p12= p11+p1;
	  p13= p12+p1;
	  p14= p13+p1;
	  p15= p14+p1;
	}

/*...The radix-16 pass is here.	*/

	for(j = 0; j < n16; j += 2)
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

		RADIX_16_DIT(a[j1],a[j2],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
					,a[j1],a[j2],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
					,c,s)

	#else

	  #ifdef USE_FGT61

	/* gather the needed data (8 64-bit complex, i.e. 16 64-bit reals) and do the first set of four length-4 transforms.
	   We process the sincos data in bit-reversed order.	*/

	// Since DIT applies twiddles to outputs, first set of 4 x radix-4 identical to radix16_dit_pass:
	/*...Block 1: */
		jt = j1;		jp = j2;
		t1 =a[jt ];		t2 =a[jp ];					m1 =b[jt ];		m2 =b[jp ];
		rt =a[jt+p1 ];	it =a[jp+p1 ];				rm =b[jt+p1 ];	im =b[jp+p1 ];
		t3 =t1 -rt;		t1 =t1 +rt;					m3 =m1 -rm;		m1 =m1 +rm;		// 1,2 in 0,2b
		t4 =t2 -it;		t2 =t2 +it;					m4 =m2 -im;		m2 =m2 +im;		// 3,4 in -b,b

		t5 =a[jt+p2 ];	t6 =a[jp+p2 ];				m5 =b[jt+p2 ];	m6 =b[jp+p2 ];
		rt =a[jt+p3 ];	it =a[jp+p3 ];				rm =b[jt+p3 ];	im =b[jp+p3 ];
		t7 =t5 -rt;		t5 =t5 +rt;					m7 =m5 -rm;		m5 =m5 +rm;		// 5,6 in 0,2b
		t8 =t6 -it;		t6 =t6 +it;					m8 =m6 -im;		m6 =m6 +im;		// 7,8 in -b,b

		rt =t5;	t5 =t1 -rt ;	t1 =t1 +rt;			rm =m5;	m5 =m1 -rm ;	m1 =m1 +rm;	// 1,2 in   0,4b
		it =t6;	t6 =t2 -it ;	t2 =t2 +it;			im =m6;	m6 =m2 -im ;	m2 =m2 +im;	// 5,6 in -2b,2b

		rt =t7;	t7 =t3 -t8 ;	t3 =t3 +t8;			rm =m7;	m7 =m3 -m8 ;	m3 =m3 +m8;	// 3,4,7,8 all in -2b,2b
				t8 =t4 +rt ;	t4 =t4 -rt;					m8 =m4 +rm ;	m4 =m4 -rm;	// rest of 4-DFTs same

	/*...Block 2: */
		jt = j1 + p4;		jp = j2 + p4;
		t9 =a[jt    ];	t10=a[jp    ];				m9 =b[jt    ];	m10=b[jp    ];
		rt =a[jt+p1 ];	it =a[jp+p1 ];				rm =b[jt+p1 ];	im =b[jp+p1 ];
		t11=t9 -rt;		t9 =t9 +rt;					m11=m9 -rm;		m9 =m9 +rm;		//  9,10 in 0,2b
		t12=t10-it;		t10=t10+it;					m12=m10-im;		m10=m10+im;		// 11,12 in -b,b

		t13=a[jt+p2 ];	t14=a[jp+p2 ];				m13=b[jt+p2 ];	m14=b[jp+p2 ];
		rt =a[jt+p3 ];	it =a[jp+p3 ];				rm =b[jt+p3 ];	im =b[jp+p3 ];
		t15=t13-rt;		t13=t13+rt;					m15=m13-rm;		m13=m13+rm;		// 13,14 in 0,2b
		t16=t14-it;		t14=t14+it;					m16=m14-im;		m14=m14+im;		// 15,16 in -b,b

		rt =t13;	t13=t9 -rt ;	t9 =t9 +rt;		rm =m13;	m13=m9 -rm ;	m9 =m9 +rm;
		it =t14;	t14=t10-it ;	t10=t10+it;		im =m14;	m14=m10-im ;	m10=m10+im;

		rt =t15;	t15=t11-t16;	t11=t11+t16;	rm =m15;	m15=m11-m16;	m11=m11+m16;
					t16=t12+rt ;	t12=t12-rt;					m16=m12+rm ;	m12=m12-rm;

	/*...Block 3: */
		jt = j1 + p8;		jp = j2 + p8;
		t17=a[jt    ];	t18=a[jp    ];				m17=b[jt    ];	m18=b[jp    ];
		rt =a[jt+p1 ];	it =a[jp+p1 ];				rm =b[jt+p1 ];	im =b[jp+p1 ];
		t19=t17-rt;		t17=t17+rt;					m19=m17-rm;		m17=m17+rm;		// 17,18 in 0,2b
		t20=t18-it;		t18=t18+it;					m20=m18-im;		m18=m18+im;		// 19,20 in -b,b

		t21=a[jt+p2 ];	t22=a[jp+p2 ];				m21=b[jt+p2 ];	m22=b[jp+p2 ];
		rt =a[jt+p3 ];	it =a[jp+p3 ];				rm =b[jt+p3 ];	im =b[jp+p3 ];
		t23=t21-rt;		t21=t21+rt;					m23=m21-rm;		m21=m21+rm;	// 21,22 in 0,2b
		t24=t22-it;		t22=t22+it;					m24=m22-im;		m22=m22+im;	// 23,24 in -b,b

		rt =t21;	t21=t17-rt ;	t17=t17+rt;		rm =m21;	m21=m17-rm ;	m17=m17+rm;	// 17,18 in   0,4b
		it =t22;	t22=t18-it ;	t18=t18+it;		im =m22;	m22=m18-im ;	m18=m18+im;	// 21,22 in -2b,2b
													// Prior to qreduce(...+q4), all 4 outs here in -2b,2b:
		rt =t23;	t23=t19-t24;	t19=t19+t24;	rm =m23;	m23=qreduce(m19-m24+q4);	m19=qreduce(m19+m24+q4);
					t24=t20+rt ;	t20=t20-rt;					m24=qreduce(m20+rm +q4);	m20=qreduce(m20-rm +q4);
													// m19,20,23,24 all needed for CMUL, so reduce.
	/*...Block 4: */
		jt = j1 + p12;		jp = j2 + p12;
		t25=a[jt    ];	t26=a[jp    ];				m25=b[jt    ];	m26=b[jp    ];
		rt =a[jt+p1 ];	it =a[jp+p1 ];				rm =b[jt+p1 ];	im =b[jp+p1 ];
		t27=t25-rt;		t25=t25+rt;					m27=m25-rm;		m25=m25+rm;	// 25,26 in 0,2b
		t28=t26-it;		t26=t26+it;					m28=m26-im;		m26=m26+im;	// 27,28 in -b,b

		t29=a[jt+p2 ];	t30=a[jp+p2 ];				m29=b[jt+p2 ];	m30=b[jp+p2 ];
		rt =a[jt+p3 ];	it =a[jp+p3 ];				rm =b[jt+p3 ];	im =b[jp+p3 ];
		t31=t29-rt;		t29=t29+rt;					m31=m29-rm;		m29=m29+rm;	// 29,30 in 0,2b
		t32=t30-it;		t30=t30+it;					m32=m30-im;		m30=m30+im;	// 31,32 in -b,b

		rt =t29;	t29=t25-rt ;	t25=t25+rt;		rm =m29;	m29=m25-rm ;	m25=m25+rm;
		it =t30;	t30=t26-it ;	t26=t26+it;		im =m30;	m30=m26-im ;	m26=m26+im;
													// Prior to qreduce(...+q4), all 4 outs here in -2b,2b:
		rt =t31;	t31=t27-t32;	t27=t27+t32;	rm =m31;	m31=qreduce(m27-m32+q4);	m27=qreduce(m27+m32+q4);
					t32=t28+rt ;	t28=t28-rt;					m32=qreduce(m28+rm +q4);	m28=qreduce(m28-rm +q4);
													// m27,28,31,32 all needed for CMUL, so reduce.
	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	*/
													/*===============
													Input bounds on modular terms:
														-2b,2b: m3-8, 11-16, 21,22, 29,30
														  0,4b: m1,2, 9 ,10, 17,18, 25,26
													These are all reduced (in 0,b) because they are inputs to CMUL:
																m19,20,23,24,27,28,31,32
													===============*/
		// And since radix16_dit_pass1 is twiddleless version, use it as template for phase 2:
	/*...Block 1: t1,9,17,25	*/
/*
printf("Block 1 float/int inputs:\n");
printf("1 ,2  float = [%10.5f,%10.5f]; int = [%" PRIu64 ",%" PRIu64 "]; neg = [%" PRIu64 ",%" PRIu64 "]\n",t1 ,t2 , m1 ,m2 , q-qreduce_full(m1 ),q-qreduce_full(m2 ));
printf("9 ,10 float = [%10.5f,%10.5f]; int = [%" PRIu64 ",%" PRIu64 "]; neg = [%" PRIu64 ",%" PRIu64 "]\n",t9 ,t10, m9 ,m10, q-qreduce_full(m9 ),q-qreduce_full(m10));
printf("17,18 float = [%10.5f,%10.5f]; int = [%" PRIu64 ",%" PRIu64 "]; neg = [%" PRIu64 ",%" PRIu64 "]\n",t17,t18, m17,m18, q-qreduce_full(m17),q-qreduce_full(m18));
printf("25,26 float = [%10.5f,%10.5f]; int = [%" PRIu64 ",%" PRIu64 "]; neg = [%" PRIu64 ",%" PRIu64 "]\n",t25,t26, m25,m26, q-qreduce_full(m25),q-qreduce_full(m26));
*/
		rt =t9 ;	t9 =t1 -rt;	t1 =t1 +rt;			rm =m9 ;	m9 =qreduce(m1 -rm+q4);	m1 =qreduce(m1 +rm);	// +:   0,8b -> 0,b
		it =t10;	t10=t2 -it;	t2 =t2 +it;			im =m10;	m10=qreduce(m2 -im+q4);	m2 =qreduce(m2 +im);	// -: -4b,4b -> 0,b

		rt =t25;	t25=t17-rt;	t17=t17+rt;			rm =m25;	m25=qreduce(m17-rm+q4);	m17=qreduce(m17+rm);	// same as above quartet
		it =t26;	t26=t18-it;	t18=t18+it;			im =m26;	m26=qreduce(m18-im+q4);	m18=qreduce(m18+im);

		a[j1    ]=t1 +t17;	a[j2    ]=t2 +t18;		b[j1    ]=qreduce(m1 +m17  );	b[j2    ]=qreduce(m2 +m18  );	// treat 0,b as == 0,q here
		a[j1+p8 ]=t1 -t17;	a[j2+p8 ]=t2 -t18;		b[j1+p8 ]=qreduce(m1 -m17+q);	b[j2+p8 ]=qreduce(m2 -m18+q);	// ~(1/2^58 odds of failure)

		a[j1+p4 ]=t9 +t26;	a[j2+p4 ]=t10-t25;		b[j1+p4 ]=qreduce(m9 +m26  );	b[j2+p4 ]=qreduce(m10-m25+q);
		a[j1+p12]=t9 -t26;	a[j2+p12]=t10+t25;		b[j1+p12]=qreduce(m9 -m26+q);	b[j2+p12]=qreduce(m10+m25  );

	/*...Block 3: t5,13,21,29	*/
/*
printf("Block 3 float/int inputs:\n");
printf("5 ,6  float = [%10.5f,%10.5f]; int = [%" PRIu64 ",%" PRIu64 "]; neg = [%" PRIu64 ",%" PRIu64 "]\n",t5 ,t6 , m5 ,m6 , q-qreduce_full(m5 ),q-qreduce_full(m6 ));
printf("13,14 float = [%10.5f,%10.5f]; int = [%" PRIu64 ",%" PRIu64 "]; neg = [%" PRIu64 ",%" PRIu64 "]\n",t13,t14, m13,m14, q-qreduce_full(m13),q-qreduce_full(m14));
printf("21,22 float = [%10.5f,%10.5f]; int = [%" PRIu64 ",%" PRIu64 "]; neg = [%" PRIu64 ",%" PRIu64 "]\n",t21,t22, m21,m22, q-qreduce_full(m21),q-qreduce_full(m22));
printf("29,30 float = [%10.5f,%10.5f]; int = [%" PRIu64 ",%" PRIu64 "]; neg = [%" PRIu64 ",%" PRIu64 "]\n",t29,t30, m29,m30, q-qreduce_full(m29),q-qreduce_full(m30));
*/
		rt =t13;	t13=t5 -t14;	t5 =t5 +t14;	rm =m13;m13=qreduce(m5-m14+q4);	m5 =qreduce(m5 +m14+q4);	// all 4 outs in -4b,4b;
					t14=t6 +rt;		t6 =t6 -rt;				m14=qreduce(m6+rm +q4);	m6 =qreduce(m6 -rm +q4);	// reduce all 4 to 0,b.

		rt =(t22+t21)*ISRT2;t22=(t22-t21)*ISRT2;	rm = mul_i2(m22+m21+q4);	m22= mul_i2(m22-m21+q4);	m21=rm;
t21=rt;	rt =(t29-t30)*ISRT2;it =(t29+t30)*ISRT2;	rm = mul_i2(m29-m30+q4);	im = mul_i2(m29+m30+q4);	// m21,22,rm,im in 0,b30
		t29=t21+rt;		t21=t21-rt;					m29 = m21+rm;		m21 = m21-rm;	// m21,22 in -b30,b30
		t30=t22+it;		t22=t22-it;					m30 = m22+im;		m22 = m22-im;	// m29,30 in 0,2*b30

		a[j1+p2 ]=t5+t21;	a[j2+p2 ]=t6+t22;		b[j1+p2 ]=qreduce(m5 +m21+q2);	b[j2+p2 ]=qreduce(m6 +m22+q2);	// + in -b30,b+b30
		a[j1+p10]=t5-t21;	a[j2+p10]=t6-t22;		b[j1+p10]=qreduce(m5 -m21+q2);	b[j2+p10]=qreduce(m6 -m22+q2);	// - same; reduce all to 0,b

		a[j1+p6 ]=t13+t30;	a[j2+p6 ]=t14-t29;		b[j1+p6 ]=qreduce(m13+m30   );	b[j2+p6 ]=qreduce(m14-m29+q3);	// + in 0,b+2*b30
		a[j1+p14]=t13-t30;	a[j2+p14]=t14+t29;		b[j1+p14]=qreduce(m13-m30+q3);	b[j2+p14]=qreduce(m14+m29   );	// - in -2*b30, b

	/*...Block 2: t3,11,19,27	*/
/*
printf("Block 2 float/int inputs:\n");
printf("3 ,4  float = [%10.5f,%10.5f]; int = [%" PRIu64 ",%" PRIu64 "]; neg = [%" PRIu64 ",%" PRIu64 "]\n",t3 ,t4 , m3 ,m4 , q-qreduce_full(m3 ),q-qreduce_full(m4 ));
printf("11,12 float = [%10.5f,%10.5f]; int = [%" PRIu64 ",%" PRIu64 "]; neg = [%" PRIu64 ",%" PRIu64 "]\n",t11,t12, m11,m12, q-qreduce_full(m11),q-qreduce_full(m12));
printf("19,20 float = [%10.5f,%10.5f]; int = [%" PRIu64 ",%" PRIu64 "]; neg = [%" PRIu64 ",%" PRIu64 "]\n",t19,t20, m19,m20, q-qreduce_full(m19),q-qreduce_full(m20));
printf("27,28 float = [%10.5f,%10.5f]; int = [%" PRIu64 ",%" PRIu64 "]; neg = [%" PRIu64 ",%" PRIu64 "]\n",t27,t28, m27,m28, q-qreduce_full(m27),q-qreduce_full(m28));
*/
		rt =(t12+t11)*ISRT2;it =(t12-t11)*ISRT2;	rm = mul_i2(m12+m11+q4);im = mul_i2(m12-m11+q4);	// 0,b30
		t11 = t3 -rt;		t3 = t3 +rt;			m11 = m3 -rm;		m3 = m3 +rm;	//  3, 4 in -2b,2b+b30
		t12 = t4 -it;		t4 = t4 +it;			m12 = m4 -im;		m4 = m4 +im;	// 11,12 in -2b-b30,2b
											// Remember, mod-root [cm,sm] premultiplied by 8, thus negate by subbing from q*8:
		rt =t19*c + t20*s;	t20=t20*c - t19*s;		cmul_modq8(m19,m20, cm,q8-sm, &m19,&m20);
t19=rt;	rt =t27*s + t28*c;	it =t28*s - t27*c;		cmul_modq8(m27,m28, sm,q8-cm, &rm ,&im );
		t27 = t19-rt;		t19 = t19+rt;			m27 =qreduce(m19-rm+q4);	m19 =qreduce(m19+rm);	// +:   0,8q ==> 0,b
		t28 = t20-it;		t20 = t20+it;			m28 =qreduce(m20-im+q4);	m20 =qreduce(m20+im);	// -: -4q,4q ==> 0,b
						// NB: +=3q for rm,im below too large, since odds of 4b+b30+3q overflowing 64-bits much larger than those of +2q not being enough to make result >= 0; analogously for m3,4)
		a[j1+p1 ]=t3+t19;	a[j2+p1 ]=t4+t20;		b[j1+p1 ]=qreduce(m3 +m19+q3);	b[j2+p1 ]=qreduce(m4 +m20+q3);	// + in [-2b,2b+b30] + [0,b] = [-2b,3b+b30], add 3q and reduce
		a[j1+p9 ]=t3-t19;	a[j2+p9 ]=t4-t20;		b[j1+p9 ]=qreduce(m3 -m19+q4);	b[j2+p9 ]=qreduce(m4 -m20+q4);	// - in [-2b,2b+b30] - [0,b] = [-3b,2b+b30], add 4q and reduce

		a[j1+p5 ]=t11+t28;	a[j2+p5 ]=t12-t27;		b[j1+p5 ]=qreduce(m11+m28+q4);	b[j2+p5 ]=qreduce(m12-m27+q4);	// +- in [-2b-b30,2b] +- [0,b] = -2b-b30,3b;
		a[j1+p13]=t11-t28;	a[j2+p13]=t12+t27;		b[j1+p13]=qreduce(m11-m28+q4);	b[j2+p13]=qreduce(m12+m27+q4);	// add 4q and reduce.
													/*===============
													Input bounds on modular terms:
														-2b,2b: m3-8, 11-16, 21,22, 29,30
														  0,4b: m1,2, 9 ,10, 17,18, 25,26
													These are all reduced (in 0,b) because they are inputs to CMUL:
																m19,20,23,24,27,28,31,32
													===============*/

	/*...Block 4: t7,15,23,31	*/
/*
printf("Block 4 float/int inputs:\n");
printf(" 7, 8 float = [%10.5f,%10.5f]; int = [%" PRIu64 ",%" PRIu64 "]; neg = [%" PRIu64 ",%" PRIu64 "]\n",t7 ,t8 , m7 ,m8 , q-qreduce_full(m7 ),q-qreduce_full(m8 ));
printf("15,16 float = [%10.5f,%10.5f]; int = [%" PRIu64 ",%" PRIu64 "]; neg = [%" PRIu64 ",%" PRIu64 "]\n",t15,t16, m15,m16, q-qreduce_full(m15),q-qreduce_full(m16));
printf("23,24 float = [%10.5f,%10.5f]; int = [%" PRIu64 ",%" PRIu64 "]; neg = [%" PRIu64 ",%" PRIu64 "]\n",t23,t24, m23,m24, q-qreduce_full(m23),q-qreduce_full(m24));
printf("31,32 float = [%10.5f,%10.5f]; int = [%" PRIu64 ",%" PRIu64 "]; neg = [%" PRIu64 ",%" PRIu64 "]\n",t31,t32, m31,m32, q-qreduce_full(m31),q-qreduce_full(m32));
exit(0);
*/
		rt =(t15-t16)*ISRT2;it =(t15+t16)*ISRT2;	rm = mul_i2(m15-m16+q4);im = mul_i2(m15+m16+q4);	// 0,b30
		t15 = t7 +rt;			t7 = t7 -rt;		m15 = m7 +rm;			m7 = m7 -rm;	//  7, 8 in -2b-b30,2b
		t16 = t8 +it;			t8 = t8 -it;		m16 = m8 +im;			m8 = m8 -im;	// 15,16 in -2b,2b+b30

		rt =t23*s + t24*c;	t24=t24*s - t23*c;		cmul_modq8(m23,m24, sm,q8-cm, &m23,&m24);
t23=rt;	rt =t31*c + t32*s;	it =t32*c - t31*s;		cmul_modq8(m31,m32, cm,q8-sm, &rm ,&im );
		t31 = t23+rt;			t23 = t23-rt;		m31 =qreduce(m23+rm);	m23 =qreduce(m23-rm+q4);// +:   0,8q ==> 0,b
		t32 = t24+it;			t24 = t24-it;		m32 =qreduce(m24+im);	m24 =qreduce(m24-im+q4);// -: -4q,4q ==> 0,b

		a[j1+p3 ]=t7 +t23;	a[j2+p3 ]=t8 +t24;		b[j1+p3 ]=qreduce(m7 +m23+q4);	b[j2+p3 ]=qreduce(m8 +m24+q4);	// + in [-2b-b30,2b] + [0,b] = -2b-b30,3b, add 4q and reduce
		a[j1+p11]=t7 -t23;	a[j2+p11]=t8 -t24;		b[j1+p11]=qreduce(m7 -m23+q5);	b[j2+p11]=qreduce(m8 -m24+q5);	// - in [-2b-b30,2b] - [0,b] = -3b-b30,2b, add 5q and reduce

		a[j1+p7 ]=t15+t32;	a[j2+p7 ]=t16-t31;		b[j1+p7 ]=qreduce(m15+m32+q3);	b[j2+p7 ]=qreduce(m16-m31+q4);	// + in [-2b,2b+b30] + [0,b] = -2b,3b+b30, add 3q and reduce
		a[j1+p15]=t15-t32;	a[j2+p15]=t16+t31;		b[j1+p15]=qreduce(m15-m32+q4);	b[j2+p15]=qreduce(m16+m31+q3);	// - in [-2b,2b+b30] - [0,b] = -3b,2b+b30, add 4q and reduce

			/**********************************************/
	  #else	// USE_FGT61 = False; Basic scalar-double mode:
			/**********************************************/

		/*...Block 1:	*/
		t1 =a[j1    ];	t2 =a[j2    ];
		rt =a[j1+p1 ];	it =a[j2+p1 ];
		t3 =t1 -rt;  	t1 =t1 +rt;
		t4 =t2 -it;		t2 =t2 +it;

		t5 =a[j1+p2 ];	t6 =a[j2+p2 ];
		rt =a[j1+p3 ];	it =a[j2+p3 ];
		t7 =t5 -rt;  	t5 =t5 +rt;
		t8 =t6 -it;  	t6 =t6 +it;

		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;

		rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
				t8 =t4 +rt;	t4 =t4 -rt;

		/*...Block 2:	*/
		t9 =a[j1+p4 ];	t10=a[j2+p4 ];
		rt =a[j1+p5 ];	it =a[j2+p5 ];
		t11=t9 -rt;  	t9 =t9 +rt;
		t12=t10-it;		t10=t10+it;

		t13=a[j1+p6 ];	t14=a[j2+p6 ];
		rt =a[j1+p7 ];	it =a[j2+p7 ];
		t15=t13-rt;  	t13=t13+rt;
		t16=t14-it;		t14=t14+it;

		rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
		it =t14;	t14=t10-it;	t10=t10+it;

		rt =t15;	t15=t11-t16;t11=t11+t16;
					t16=t12+rt;	t12=t12-rt;

		/*...Block 3:	*/
		t17=a[j1+p8 ];	t18=a[j2+p8 ];
		rt =a[j1+p9 ];	it =a[j2+p9 ];
		t19=t17-rt;  	t17=t17+rt;
		t20=t18-it;		t18=t18+it;

		t21=a[j1+p10];	t22=a[j2+p10];
		rt =a[j1+p11];	it =a[j2+p11];
		t23=t21-rt;  	t21=t21+rt;
		t24=t22-it;		t22=t22+it;

		rt =t21;	t21=t17-rt;	t17=t17+rt;
		it =t22;	t22=t18-it;	t18=t18+it;

		rt =t23;	t23=t19-t24;t19=t19+t24;
					t24=t20+rt;	t20=t20-rt;

		/*...Block 4:	*/
		t25=a[j1+p12];	t26=a[j2+p12];
		rt =a[j1+p13];	it =a[j2+p13];
		t27=t25-rt;  	t25=t25+rt;
		t28=t26-it;		t26=t26+it;

		t29=a[j1+p14];	t30=a[j2+p14];
		rt =a[j1+p15];	it =a[j2+p15];
		t31=t29-rt;  	t29=t29+rt;
		t32=t30-it;		t30=t30+it;

		rt =t29;	t29=t25-rt;	t25=t25+rt;
		it =t30;	t30=t26-it;	t26=t26+it;

		rt =t31;	t31=t27-t32;t27=t27+t32;
					t32=t28+rt;	t28=t28-rt;

		/*...and now do four more radix-4 transforms, including the internal twiddle factors:
			1, exp(-i* 1*twopi/16) =       ( c,-s), exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 3*twopi/16) =       ( s,-c) (for inputs to transform block 2)
			1, exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 4*twopi/16) =       ( 0,-1), exp(-i* 6*twopi/16) = ISRT2*(-1,-1) (for inputs to transform block 3)
			1, exp(-i* 3*twopi/16) =       ( s,-c), exp(-i* 6*twopi/16) = ISRT2*(-1,-1), exp(-i* 9*twopi/16) =       (-c, s) (for inputs to transform block 4).
			(This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
			 I.e. do similar as above, except inputs a(j1  +p0:15:1) are replaced by t0:30:2,
													 a(j2+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
			 and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.	*/

		/*...Block 1: t1,9,17,25	*/
		rt =t9 ;	t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		a[j1    ]=t1 +t17;	a[j2    ]=t2 +t18;
		a[j1+p8 ]=t1 -t17;	a[j2+p8 ]=t2 -t18;

		a[j1+p4 ]=t9 +t26;	a[j2+p4 ]=t10-t25;	/* mpy by E^-4 = -I is inlined here...	*/
		a[j1+p12]=t9 -t26;	a[j2+p12]=t10+t25;

		/*...Block 3: t5,13,21,29	*/
		rt =t13;	t13=t5 -t14;t5 =t5 +t14;		/* twiddle mpy by E^4 =-I	*/
					t14=t6 +rt;	t6 =t6 -rt;

		rt =(t22+t21)*ISRT2;t22=(t22-t21)*ISRT2;	t21=rt;	/* twiddle mpy by E^-2	*/
		rt =(t29-t30)*ISRT2;it =(t29+t30)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here...	*/
		t29=t21+rt;		t21=t21-rt;			/* ...and get E^-6 by flipping signs here.	*/
		t30=t22+it;		t22=t22-it;

		a[j1+p2 ]=t5+t21;	a[j2+p2 ]=t6+t22;
		a[j1+p10]=t5-t21;	a[j2+p10]=t6-t22;

		a[j1+p6 ]=t13+t30;	a[j2+p6 ]=t14-t29;	/* mpy by E^-4 = -I is inlined here...	*/
		a[j1+p14]=t13-t30;	a[j2+p14]=t14+t29;

		/*...Block 2: t3,11,19,27	*/
		rt =(t12+t11)*ISRT2;it =(t12-t11)*ISRT2;		/* twiddle mpy by E^-2	*/
		t11=t3 -rt;		t3 =t3 +rt;
		t12=t4 -it;		t4 =t4 +it;

		rt =t19*c + t20*s;	t20=t20*c - t19*s;	t19=rt;	/* twiddle mpy by E^-1	*/
		rt =t27*s + t28*c;	it =t28*s - t27*c;		/* twiddle mpy by E^-3	*/
		t27=t19-rt;		t19=t19+rt;
		t28=t20-it;		t20=t20+it;

		a[j1+p1 ]=t3+t19;	a[j2+p1 ]=t4+t20;
		a[j1+p9 ]=t3-t19;	a[j2+p9 ]=t4-t20;

		a[j1+p5 ]=t11+t28;	a[j2+p5 ]=t12-t27;	/* mpy by E^-4 = -I is inlined here...	*/
		a[j1+p13]=t11-t28;	a[j2+p13]=t12+t27;

		/*...Block 4: t7,15,23,31	*/
		rt =(t15-t16)*ISRT2;it =(t15+t16)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here...	*/
		t15=t7 +rt;		t7 =t7 -rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here.	*/
		t16=t8 +it;		t8 =t8 -it;

		rt =t23*s + t24*c;	t24=t24*s - t23*c;	t23=rt;	/* twiddle mpy by E^-3	*/
		rt =t31*c + t32*s;	it =t32*c - t31*s;		/* twiddle mpy by E^-1 = -E^-9...	*/
		t31=t23+rt;		t23=t23-rt;			/* ...and get E^9 by flipping signs here.	*/
		t32=t24+it;		t24=t24-it;

		a[j1+p3 ]=t7+t23;	a[j2+p3 ]=t8+t24;
		a[j1+p11]=t7-t23;	a[j2+p11]=t8-t24;

		a[j1+p7 ]=t15+t32;	a[j2+p7 ]=t16-t31;	/* mpy by E^-4 = -I is inlined here...	*/
		a[j1+p15]=t15-t32;	a[j2+p15]=t16+t31;

	  #endif	// USE_FGT61 ?

	#endif
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy16_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		const int pfetch_dist = PFETCH_DIST;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		int idx_offset,idx_incr;
		int j,j1,j2,k,l;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
		int poff[RADIX>>2],p0123[4];	// Store [RADIX/4] mults of p04 offset for loop control
	#ifdef USE_AVX512
		double t0,t1,t2,t3;
		struct uint32x8 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#elif defined(USE_AVX)
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#else
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	#endif
		double rt,it, wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and AVX mode
		double temp,frac
			,a1p0r,a1p1r,a1p2r,a1p3r,a1p4r,a1p5r,a1p6r,a1p7r,a1p8r,a1p9r,a1pAr,a1pBr,a1pCr,a1pDr,a1pEr,a1pFr
			,a1p0i,a1p1i,a1p2i,a1p3i,a1p4i,a1p5i,a1p6i,a1p7i,a1p8i,a1p9i,a1pAi,a1pBi,a1pCi,a1pDi,a1pEi,a1pFi;

	#ifdef USE_SSE2

		uint32 p1,p2,p3,p4,p8,p12;
		double *addr, *add0,*add1,*add2,*add3;
		const double crnd = 3.0*0x4000000*0x2000000;
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		vec_dbl *cc0, *ss0, *isrt2, *max_err, *sse2_rnd, *half_arr, *tmp,*tm2
			,*r00,*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09,*r0A,*r0B,*r0C,*r0D,*r0E,*r0F
			,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r1A,*r1B,*r1C,*r1D,*r1E,*r1F;
	  #ifdef USE_AVX512
		vec_dbl *cy_r0,*cy_r8, *cy_i0,*cy_i8;
	  #elif defined(USE_AVX)
		vec_dbl *cy_r0,*cy_r4,*cy_r8,*cy_rC, *cy_i0,*cy_i4,*cy_i8,*cy_iC;
	  #else	// SSE2:
		vec_dbl *cy_r0,*cy_r2,*cy_r4,*cy_r6,*cy_r8,*cy_rA,*cy_rC,*cy_rE, *cy_i0,*cy_i2,*cy_i4,*cy_i6,*cy_i8,*cy_iA,*cy_iC,*cy_iE;
	  #endif
	  #ifdef USE_AVX
		int k1,k2;
		vec_dbl *base_negacyclic_root;
	  #endif
		int *bjmodn0,*bjmodn1,*bjmodn2,*bjmodn3,*bjmodn4,*bjmodn5,*bjmodn6,*bjmodn7,*bjmodn8,*bjmodn9,*bjmodnA,*bjmodnB,*bjmodnC,*bjmodnD,*bjmodnE,*bjmodnF;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_nm1;
		double dtmp;

	#else

		const double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;
		double *addr, *base, *baseinv;
		uint32 p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15;
		// Vars needed in scalar mode only:
		const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int jt,jp,k1,k2,m,m2,ntmp;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	  #if !USE_SCALAR_DFT_MACRO
	   #if PFETCH
		double *addp;
	   #endif
		double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;
	  #endif
		int bjmodn0=0,bjmodn1=0,bjmodn2=0,bjmodn3=0,bjmodn4=0,bjmodn5=0,bjmodn6=0,bjmodn7=0,bjmodn8=0,bjmodn9=0,bjmodnA=0,bjmodnB=0,bjmodnC=0,bjmodnD=0,bjmodnE=0,bjmodnF=0;
		double cy_r0=0,cy_r1=0,cy_r2=0,cy_r3=0,cy_r4=0,cy_r5=0,cy_r6=0,cy_r7=0,cy_r8=0,cy_r9=0,cy_rA=0,cy_rB=0,cy_rC=0,cy_rD=0,cy_rE=0,cy_rF=0,
			cy_i0=0,cy_i1=0,cy_i2=0,cy_i3=0,cy_i4=0,cy_i5=0,cy_i6=0,cy_i7=0,cy_i8=0,cy_i9=0,cy_iA=0,cy_iB=0,cy_iC=0,cy_iD=0,cy_iE=0,cy_iF=0;
	#ifdef USE_FGT61
		const uint64 q  = 0x1FFFFFFFFFFFFFFFull, q2=q+q, q3=q2+q, q4=q2+q2, q5=q4+q, q8=q4+q4;	// q = 2^61 - 1, and needed small multiples
		// primitive 16th root of unity, scaled by *8:
		const uint64 cm = 1693317751237720973ull<<3, sm = 2283815672160731785ull<<3;
		uint64 rm,im;
		// varname m2 already used for float-carry macros, so replace mod-varname m2 with m$ in carry step
		uint64 m1,m$,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25,m26,m27,m28,m29,m30,m31,m32
			,b1p0r,b1p1r,b1p2r,b1p3r,b1p4r,b1p5r,b1p6r,b1p7r,b1p8r,b1p9r,b1pAr,b1pBr,b1pCr,b1pDr,b1pEr,b1pFr
			,b1p0i,b1p1i,b1p2i,b1p3i,b1p4i,b1p5i,b1p6i,b1p7i,b1p8i,b1p9i,b1pAi,b1pBi,b1pCi,b1pDi,b1pEi,b1pFi;
	#endif

	#endif

		struct cy_thread_data_t* thread_arg = targ;
	// int data:
		int iter = thread_arg->iter;
		int ithread = thread_arg->tid;	/* unique thread index (use for debug) */

		int NDIVR = thread_arg->ndivr;
		int n = NDIVR*RADIX, nm1 = n-1;
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
		int sw  = thread_arg->sw, bw = n - sw;
		int nwt = thread_arg->nwt;

	// double data:
		double maxerr = thread_arg->maxerr;
		double scale = thread_arg->scale;	int full_pass = scale < 0.5;
		double prp_mult = thread_arg->prp_mult;

	// pointer data:
		double *a = thread_arg->arrdat;
	#ifdef USE_FGT61
		uint64 *b = thread_arg->brrdat;
	#endif
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
		double *wts_mult = thread_arg->wts_mult;	// Const Intra-block wts-multiplier...
		double *inv_mult = thread_arg->inv_mult;	// ...and 2*(its multiplicative inverse).
		ASSERT(fabs(wts_mult[0]*inv_mult[0] - 1.0) < EPS, "wts_mults fail accuracy check!");
		ASSERT(fabs(wts_mult[1]*inv_mult[1] - 1.0) < EPS, "wts_mults fail accuracy check!");
		int *si = thread_arg->si;
		struct complex *rn0 = thread_arg->rn0;
		struct complex *rn1 = thread_arg->rn1;

		p1  = NDIVR;
		p2  = p1  + p1;
		p3  = p2  + p1;
		p4  = p3  + p1;
		p8  = p4  + p4;
		p12 = p8  + p4;
	#ifndef USE_SSE2
		p5  = p4  + p1;
		p6  = p5  + p1;
		p7  = p6  + p1;
		p9  = p7  + p2;
		p10 = p9  + p1;
		p11 = p10 + p1;
		p13 = p11 + p2;
		p14 = p13 + p1;
		p15 = p14 + p1;
	#endif	// USE_SSE2
		p1  = p1  + ( (p1  >> DAT_BITS) << PAD_BITS );
		p2  = p2  + ( (p2  >> DAT_BITS) << PAD_BITS );
		p3  = p3  + ( (p3  >> DAT_BITS) << PAD_BITS );
		p4  = p4  + ( (p4  >> DAT_BITS) << PAD_BITS );
		p8  = p8  + ( (p8  >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
	#ifndef USE_SSE2
		p5  = p5  + ( (p5  >> DAT_BITS) << PAD_BITS );
		p6  = p6  + ( (p6  >> DAT_BITS) << PAD_BITS );
		p7  = p7  + ( (p7  >> DAT_BITS) << PAD_BITS );
		p9  = p9  + ( (p9  >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
	#endif	// USE_SSE2
		p0123[0] = 0; p0123[1] = p1; p0123[2] = p2; p0123[3] = p3;
		poff[0x0] =   0; poff[0x1] = p4    ; poff[0x2] = p8; poff[0x3] = p12;

	#ifdef USE_SSE2
		r00 = thread_arg->r00;	r01 = r00 + 0x01;
		r02 = r00 + 0x02;		r03 = r00 + 0x03;
		r04 = r00 + 0x04;		r05 = r00 + 0x05;
		r06 = r00 + 0x06;		r07 = r00 + 0x07;
		r08 = r00 + 0x08;		r09 = r00 + 0x09;
		r0A = r00 + 0x0a;		r0B = r00 + 0x0b;
		r0C = r00 + 0x0c;		r0D = r00 + 0x0d;
		r0E = r00 + 0x0e;		r0F = r00 + 0x0f;
		r10 = r00 + 0x10;		r11 = r00 + 0x11;
		r12 = r00 + 0x12;		r13 = r00 + 0x13;
		r14 = r00 + 0x14;		r15 = r00 + 0x15;
		r16 = r00 + 0x16;		r17 = r00 + 0x17;
		r18 = r00 + 0x18;		r19 = r00 + 0x19;
		r1A = r00 + 0x1a;		r1B = r00 + 0x1b;
		r1C = r00 + 0x1c;		r1D = r00 + 0x1d;
		r1E = r00 + 0x1e;		r1F = r00 + 0x1f;
	  isrt2 = r00 + 0x20;
		cc0 = r00 + 0x21;
		ss0 = r00 + 0x22;
		tmp = ss0 + 0x2;	// Only need +1 but prefer an even offset
	#ifdef USE_AVX512
		 cy_r0 = tmp + 0x00;
		 cy_r8 = tmp + 0x01;
		 cy_i0 = tmp + 0x02;
		 cy_i8 = tmp + 0x03;
		max_err = tmp + 0x04;
		sse2_rnd= tmp + 0x05;
		half_arr= tmp + 0x06;	// 36 + 6 = 42 = half_arr_offset16
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod
	#elif defined(USE_AVX)
		 cy_r0 = tmp + 0x00;
		 cy_r4 = tmp + 0x01;
		 cy_r8 = tmp + 0x02;
		 cy_rC = tmp + 0x03;
		 cy_i0 = tmp + 0x04;
		 cy_i4 = tmp + 0x05;
		 cy_i8 = tmp + 0x06;
		 cy_iC = tmp + 0x07;
		max_err = tmp + 0x08;
		sse2_rnd= tmp + 0x09;
		half_arr= tmp + 0x0a;	/* This table needs 20x16 bytes */
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod
	  #else
		 cy_r0 = tmp + 0x00;
		 cy_r2 = tmp + 0x01;
		 cy_r4 = tmp + 0x02;
		 cy_r6 = tmp + 0x03;
		 cy_r8 = tmp + 0x04;
		 cy_rA = tmp + 0x05;
		 cy_rC = tmp + 0x06;
		 cy_rE = tmp + 0x07;
		 cy_i0 = tmp + 0x08;
		 cy_i2 = tmp + 0x09;
		 cy_i4 = tmp + 0x0a;
		 cy_i6 = tmp + 0x0b;
		 cy_i8 = tmp + 0x0c;
		 cy_iA = tmp + 0x0d;
		 cy_iC = tmp + 0x0e;
		 cy_iE = tmp + 0x0f;
		max_err = tmp + 0x10;
		sse2_rnd= tmp + 0x11;
		half_arr= tmp + 0x12;	/* This table needs 20x16 bytes */
	  #endif

		ASSERT((isrt2->d0 == ISRT2 && isrt2->d1 == ISRT2), "thread-local memcheck failed!");
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

		sign_mask = (uint64*)(r00 + radix16_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (#doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		sse_sw  = sse_bw    + RE_IM_STRIDE;
		sse_nm1 = sse_sw    + RE_IM_STRIDE;
	  #ifdef USE_AVX512
		n_minus_sil   = (struct uint32x8 *)sse_nm1 + 1;
		n_minus_silp1 = (struct uint32x8 *)sse_nm1 + 2;
		sinwt         = (struct uint32x8 *)sse_nm1 + 3;
		sinwtm1       = (struct uint32x8 *)sse_nm1 + 4;
		bjmodn0 = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #elif defined(USE_AVX)
		n_minus_sil   = (struct uint32x4 *)sse_nm1 + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_nm1 + 2;
		sinwt         = (struct uint32x4 *)sse_nm1 + 3;
		sinwtm1       = (struct uint32x4 *)sse_nm1 + 4;
		bjmodn0 = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn0 = (int*)(sse_nm1 + RE_IM_STRIDE);
	  #endif
		bjmodn1 = bjmodn0 + 0x01;
		bjmodn2 = bjmodn0 + 0x02;
		bjmodn3 = bjmodn0 + 0x03;
		bjmodn4 = bjmodn0 + 0x04;
		bjmodn5 = bjmodn0 + 0x05;
		bjmodn6 = bjmodn0 + 0x06;
		bjmodn7 = bjmodn0 + 0x07;
		bjmodn8 = bjmodn0 + 0x08;
		bjmodn9 = bjmodn0 + 0x09;
		bjmodnA = bjmodn0 + 0x0A;
		bjmodnB = bjmodn0 + 0x0B;
		bjmodnC = bjmodn0 + 0x0C;
		bjmodnD = bjmodn0 + 0x0D;
		bjmodnE = bjmodn0 + 0x0E;
		bjmodnF = bjmodn0 + 0x0F;

	#else

		// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
		base    = (double *)thread_arg->r00     ;
		baseinv = (double *)thread_arg->half_arr;

	#endif	// USE_SSE2 ?

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Init DWT-indices: */	/* init carries	*/
		#ifdef USE_AVX512
			*bjmodn0 = thread_arg->bjmodn0;		cy_r0->d0 = thread_arg->cy_r0;
			*bjmodn1 = thread_arg->bjmodn1;		cy_r0->d1 = thread_arg->cy_r1;
			*bjmodn2 = thread_arg->bjmodn2;		cy_r0->d2 = thread_arg->cy_r2;
			*bjmodn3 = thread_arg->bjmodn3;		cy_r0->d3 = thread_arg->cy_r3;
			*bjmodn4 = thread_arg->bjmodn4;		cy_r0->d4 = thread_arg->cy_r4;
			*bjmodn5 = thread_arg->bjmodn5;		cy_r0->d5 = thread_arg->cy_r5;
			*bjmodn6 = thread_arg->bjmodn6;		cy_r0->d6 = thread_arg->cy_r6;
			*bjmodn7 = thread_arg->bjmodn7;		cy_r0->d7 = thread_arg->cy_r7;
			*bjmodn8 = thread_arg->bjmodn8;		cy_r8->d0 = thread_arg->cy_r8;
			*bjmodn9 = thread_arg->bjmodn9;		cy_r8->d1 = thread_arg->cy_r9;
			*bjmodnA = thread_arg->bjmodnA;		cy_r8->d2 = thread_arg->cy_rA;
			*bjmodnB = thread_arg->bjmodnB;		cy_r8->d3 = thread_arg->cy_rB;
			*bjmodnC = thread_arg->bjmodnC;		cy_r8->d4 = thread_arg->cy_rC;
			*bjmodnD = thread_arg->bjmodnD;		cy_r8->d5 = thread_arg->cy_rD;
			*bjmodnE = thread_arg->bjmodnE;		cy_r8->d6 = thread_arg->cy_rE;
			*bjmodnF = thread_arg->bjmodnF;		cy_r8->d7 = thread_arg->cy_rF;
		#elif defined(USE_AVX)
			*bjmodn0 = thread_arg->bjmodn0;		cy_r0->d0 = thread_arg->cy_r0;
			*bjmodn1 = thread_arg->bjmodn1;		cy_r0->d1 = thread_arg->cy_r1;
			*bjmodn2 = thread_arg->bjmodn2;		cy_r0->d2 = thread_arg->cy_r2;
			*bjmodn3 = thread_arg->bjmodn3;		cy_r0->d3 = thread_arg->cy_r3;
			*bjmodn4 = thread_arg->bjmodn4;		cy_r4->d0 = thread_arg->cy_r4;
			*bjmodn5 = thread_arg->bjmodn5;		cy_r4->d1 = thread_arg->cy_r5;
			*bjmodn6 = thread_arg->bjmodn6;		cy_r4->d2 = thread_arg->cy_r6;
			*bjmodn7 = thread_arg->bjmodn7;		cy_r4->d3 = thread_arg->cy_r7;
			*bjmodn8 = thread_arg->bjmodn8;		cy_r8->d0 = thread_arg->cy_r8;
			*bjmodn9 = thread_arg->bjmodn9;		cy_r8->d1 = thread_arg->cy_r9;
			*bjmodnA = thread_arg->bjmodnA;		cy_r8->d2 = thread_arg->cy_rA;
			*bjmodnB = thread_arg->bjmodnB;		cy_r8->d3 = thread_arg->cy_rB;
			*bjmodnC = thread_arg->bjmodnC;		cy_rC->d0 = thread_arg->cy_rC;
			*bjmodnD = thread_arg->bjmodnD;		cy_rC->d1 = thread_arg->cy_rD;
			*bjmodnE = thread_arg->bjmodnE;		cy_rC->d2 = thread_arg->cy_rE;
			*bjmodnF = thread_arg->bjmodnF;		cy_rC->d3 = thread_arg->cy_rF;
		#elif defined(USE_SSE2)
			*bjmodn0 = thread_arg->bjmodn0;		cy_r0->d0 = thread_arg->cy_r0;
			*bjmodn1 = thread_arg->bjmodn1;		cy_r0->d1 = thread_arg->cy_r1;
			*bjmodn2 = thread_arg->bjmodn2;		cy_r2->d0 = thread_arg->cy_r2;
			*bjmodn3 = thread_arg->bjmodn3;		cy_r2->d1 = thread_arg->cy_r3;
			*bjmodn4 = thread_arg->bjmodn4;		cy_r4->d0 = thread_arg->cy_r4;
			*bjmodn5 = thread_arg->bjmodn5;		cy_r4->d1 = thread_arg->cy_r5;
			*bjmodn6 = thread_arg->bjmodn6;		cy_r6->d0 = thread_arg->cy_r6;
			*bjmodn7 = thread_arg->bjmodn7;		cy_r6->d1 = thread_arg->cy_r7;
			*bjmodn8 = thread_arg->bjmodn8;		cy_r8->d0 = thread_arg->cy_r8;
			*bjmodn9 = thread_arg->bjmodn9;		cy_r8->d1 = thread_arg->cy_r9;
			*bjmodnA = thread_arg->bjmodnA;		cy_rA->d0 = thread_arg->cy_rA;
			*bjmodnB = thread_arg->bjmodnB;		cy_rA->d1 = thread_arg->cy_rB;
			*bjmodnC = thread_arg->bjmodnC;		cy_rC->d0 = thread_arg->cy_rC;
			*bjmodnD = thread_arg->bjmodnD;		cy_rC->d1 = thread_arg->cy_rD;
			*bjmodnE = thread_arg->bjmodnE;		cy_rE->d0 = thread_arg->cy_rE;
			*bjmodnF = thread_arg->bjmodnF;		cy_rE->d1 = thread_arg->cy_rF;
		#else
			bjmodn0 = thread_arg->bjmodn0;		cy_r0 = thread_arg->cy_r0;
			bjmodn1 = thread_arg->bjmodn1;		cy_r1 = thread_arg->cy_r1;
			bjmodn2 = thread_arg->bjmodn2;		cy_r2 = thread_arg->cy_r2;
			bjmodn3 = thread_arg->bjmodn3;		cy_r3 = thread_arg->cy_r3;
			bjmodn4 = thread_arg->bjmodn4;		cy_r4 = thread_arg->cy_r4;
			bjmodn5 = thread_arg->bjmodn5;		cy_r5 = thread_arg->cy_r5;
			bjmodn6 = thread_arg->bjmodn6;		cy_r6 = thread_arg->cy_r6;
			bjmodn7 = thread_arg->bjmodn7;		cy_r7 = thread_arg->cy_r7;
			bjmodn8 = thread_arg->bjmodn8;		cy_r8 = thread_arg->cy_r8;
			bjmodn9 = thread_arg->bjmodn9;		cy_r9 = thread_arg->cy_r9;
			bjmodnA = thread_arg->bjmodnA;		cy_rA = thread_arg->cy_rA;
			bjmodnB = thread_arg->bjmodnB;		cy_rB = thread_arg->cy_rB;
			bjmodnC = thread_arg->bjmodnC;		cy_rC = thread_arg->cy_rC;
			bjmodnD = thread_arg->bjmodnD;		cy_rD = thread_arg->cy_rD;
			bjmodnE = thread_arg->bjmodnE;		cy_rE = thread_arg->cy_rE;
			bjmodnF = thread_arg->bjmodnF;		cy_rF = thread_arg->cy_rF;
		#endif
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
			/* init carries	*/
		#ifdef USE_AVX512
			cy_r0->d0 = thread_arg->cy_r0;		cy_i0->d0 = thread_arg->cy_i0;
			cy_r0->d1 = thread_arg->cy_r1;		cy_i0->d1 = thread_arg->cy_i1;
			cy_r0->d2 = thread_arg->cy_r2;		cy_i0->d2 = thread_arg->cy_i2;
			cy_r0->d3 = thread_arg->cy_r3;		cy_i0->d3 = thread_arg->cy_i3;
			cy_r0->d4 = thread_arg->cy_r4;		cy_i0->d4 = thread_arg->cy_i4;
			cy_r0->d5 = thread_arg->cy_r5;		cy_i0->d5 = thread_arg->cy_i5;
			cy_r0->d6 = thread_arg->cy_r6;		cy_i0->d6 = thread_arg->cy_i6;
			cy_r0->d7 = thread_arg->cy_r7;		cy_i0->d7 = thread_arg->cy_i7;
			cy_r8->d0 = thread_arg->cy_r8;		cy_i8->d0 = thread_arg->cy_i8;
			cy_r8->d1 = thread_arg->cy_r9;		cy_i8->d1 = thread_arg->cy_i9;
			cy_r8->d2 = thread_arg->cy_rA;		cy_i8->d2 = thread_arg->cy_iA;
			cy_r8->d3 = thread_arg->cy_rB;		cy_i8->d3 = thread_arg->cy_iB;
			cy_r8->d4 = thread_arg->cy_rC;		cy_i8->d4 = thread_arg->cy_iC;
			cy_r8->d5 = thread_arg->cy_rD;		cy_i8->d5 = thread_arg->cy_iD;
			cy_r8->d6 = thread_arg->cy_rE;		cy_i8->d6 = thread_arg->cy_iE;
			cy_r8->d7 = thread_arg->cy_rF;		cy_i8->d7 = thread_arg->cy_iF;
		#elif defined(USE_AVX)
			cy_r0->d0 = thread_arg->cy_r0;		cy_i0->d0 = thread_arg->cy_i0;
			cy_r0->d1 = thread_arg->cy_r1;		cy_i0->d1 = thread_arg->cy_i1;
			cy_r0->d2 = thread_arg->cy_r2;		cy_i0->d2 = thread_arg->cy_i2;
			cy_r0->d3 = thread_arg->cy_r3;		cy_i0->d3 = thread_arg->cy_i3;
			cy_r4->d0 = thread_arg->cy_r4;		cy_i4->d0 = thread_arg->cy_i4;
			cy_r4->d1 = thread_arg->cy_r5;		cy_i4->d1 = thread_arg->cy_i5;
			cy_r4->d2 = thread_arg->cy_r6;		cy_i4->d2 = thread_arg->cy_i6;
			cy_r4->d3 = thread_arg->cy_r7;		cy_i4->d3 = thread_arg->cy_i7;
			cy_r8->d0 = thread_arg->cy_r8;		cy_i8->d0 = thread_arg->cy_i8;
			cy_r8->d1 = thread_arg->cy_r9;		cy_i8->d1 = thread_arg->cy_i9;
			cy_r8->d2 = thread_arg->cy_rA;		cy_i8->d2 = thread_arg->cy_iA;
			cy_r8->d3 = thread_arg->cy_rB;		cy_i8->d3 = thread_arg->cy_iB;
			cy_rC->d0 = thread_arg->cy_rC;		cy_iC->d0 = thread_arg->cy_iC;
			cy_rC->d1 = thread_arg->cy_rD;		cy_iC->d1 = thread_arg->cy_iD;
			cy_rC->d2 = thread_arg->cy_rE;		cy_iC->d2 = thread_arg->cy_iE;
			cy_rC->d3 = thread_arg->cy_rF;		cy_iC->d3 = thread_arg->cy_iF;
		#elif defined(USE_SSE2)
			cy_r0->d0 = thread_arg->cy_r0;		cy_r0->d1 = thread_arg->cy_i0;
			cy_r2->d0 = thread_arg->cy_r1;		cy_r2->d1 = thread_arg->cy_i1;
			cy_r4->d0 = thread_arg->cy_r2;		cy_r4->d1 = thread_arg->cy_i2;
			cy_r6->d0 = thread_arg->cy_r3;		cy_r6->d1 = thread_arg->cy_i3;
			cy_r8->d0 = thread_arg->cy_r4;		cy_r8->d1 = thread_arg->cy_i4;
			cy_rA->d0 = thread_arg->cy_r5;		cy_rA->d1 = thread_arg->cy_i5;
			cy_rC->d0 = thread_arg->cy_r6;		cy_rC->d1 = thread_arg->cy_i6;
			cy_rE->d0 = thread_arg->cy_r7;		cy_rE->d1 = thread_arg->cy_i7;
			cy_i0->d0 = thread_arg->cy_r8;		cy_i0->d1 = thread_arg->cy_i8;
			cy_i2->d0 = thread_arg->cy_r9;		cy_i2->d1 = thread_arg->cy_i9;
			cy_i4->d0 = thread_arg->cy_rA;		cy_i4->d1 = thread_arg->cy_iA;
			cy_i6->d0 = thread_arg->cy_rB;		cy_i6->d1 = thread_arg->cy_iB;
			cy_i8->d0 = thread_arg->cy_rC;		cy_i8->d1 = thread_arg->cy_iC;
			cy_iA->d0 = thread_arg->cy_rD;		cy_iA->d1 = thread_arg->cy_iD;
			cy_iC->d0 = thread_arg->cy_rE;		cy_iC->d1 = thread_arg->cy_iE;
			cy_iE->d0 = thread_arg->cy_rF;		cy_iE->d1 = thread_arg->cy_iF;
		#else
			cy_r0 = thread_arg->cy_r0;			cy_i0 = thread_arg->cy_i0;
			cy_r1 = thread_arg->cy_r1;			cy_i1 = thread_arg->cy_i1;
			cy_r2 = thread_arg->cy_r2;			cy_i2 = thread_arg->cy_i2;
			cy_r3 = thread_arg->cy_r3;			cy_i3 = thread_arg->cy_i3;
			cy_r4 = thread_arg->cy_r4;			cy_i4 = thread_arg->cy_i4;
			cy_r5 = thread_arg->cy_r5;			cy_i5 = thread_arg->cy_i5;
			cy_r6 = thread_arg->cy_r6;			cy_i6 = thread_arg->cy_i6;
			cy_r7 = thread_arg->cy_r7;			cy_i7 = thread_arg->cy_i7;
			cy_r8 = thread_arg->cy_r8;			cy_i8 = thread_arg->cy_i8;
			cy_r9 = thread_arg->cy_r9;			cy_i9 = thread_arg->cy_i9;
			cy_rA = thread_arg->cy_rA;			cy_iA = thread_arg->cy_iA;
			cy_rB = thread_arg->cy_rB;			cy_iB = thread_arg->cy_iB;
			cy_rC = thread_arg->cy_rC;			cy_iC = thread_arg->cy_iC;
			cy_rD = thread_arg->cy_rD;			cy_iD = thread_arg->cy_iD;
			cy_rE = thread_arg->cy_rE;			cy_iE = thread_arg->cy_iE;
			cy_rF = thread_arg->cy_rF;			cy_iF = thread_arg->cy_iF;
		#endif
		}

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix16_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX512
			thread_arg->cy_r0 = cy_r0->d0;
			thread_arg->cy_r1 = cy_r0->d1;
			thread_arg->cy_r2 = cy_r0->d2;
			thread_arg->cy_r3 = cy_r0->d3;
			thread_arg->cy_r4 = cy_r0->d4;
			thread_arg->cy_r5 = cy_r0->d5;
			thread_arg->cy_r6 = cy_r0->d6;
			thread_arg->cy_r7 = cy_r0->d7;
			thread_arg->cy_r8 = cy_r8->d0;
			thread_arg->cy_r9 = cy_r8->d1;
			thread_arg->cy_rA = cy_r8->d2;
			thread_arg->cy_rB = cy_r8->d3;
			thread_arg->cy_rC = cy_r8->d4;
			thread_arg->cy_rD = cy_r8->d5;
			thread_arg->cy_rE = cy_r8->d6;
			thread_arg->cy_rF = cy_r8->d7;
			t0 = MAX(max_err->d0,max_err->d1);
			t1 = MAX(max_err->d2,max_err->d3);
			t2 = MAX(max_err->d4,max_err->d5);
			t3 = MAX(max_err->d6,max_err->d7);
			maxerr = MAX( MAX(t0,t1), MAX(t2,t3) );
		#elif defined(USE_AVX)
			thread_arg->cy_r0 = cy_r0->d0;
			thread_arg->cy_r1 = cy_r0->d1;
			thread_arg->cy_r2 = cy_r0->d2;
			thread_arg->cy_r3 = cy_r0->d3;
			thread_arg->cy_r4 = cy_r4->d0;
			thread_arg->cy_r5 = cy_r4->d1;
			thread_arg->cy_r6 = cy_r4->d2;
			thread_arg->cy_r7 = cy_r4->d3;
			thread_arg->cy_r8 = cy_r8->d0;
			thread_arg->cy_r9 = cy_r8->d1;
			thread_arg->cy_rA = cy_r8->d2;
			thread_arg->cy_rB = cy_r8->d3;
			thread_arg->cy_rC = cy_rC->d0;
			thread_arg->cy_rD = cy_rC->d1;
			thread_arg->cy_rE = cy_rC->d2;
			thread_arg->cy_rF = cy_rC->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			thread_arg->cy_r0 = cy_r0->d0;
			thread_arg->cy_r1 = cy_r0->d1;
			thread_arg->cy_r2 = cy_r2->d0;
			thread_arg->cy_r3 = cy_r2->d1;
			thread_arg->cy_r4 = cy_r4->d0;
			thread_arg->cy_r5 = cy_r4->d1;
			thread_arg->cy_r6 = cy_r6->d0;
			thread_arg->cy_r7 = cy_r6->d1;
			thread_arg->cy_r8 = cy_r8->d0;
			thread_arg->cy_r9 = cy_r8->d1;
			thread_arg->cy_rA = cy_rA->d0;
			thread_arg->cy_rB = cy_rA->d1;
			thread_arg->cy_rC = cy_rC->d0;
			thread_arg->cy_rD = cy_rC->d1;
			thread_arg->cy_rE = cy_rE->d0;
			thread_arg->cy_rF = cy_rE->d1;
			maxerr = MAX(max_err->d0,max_err->d1);
		#else
			thread_arg->cy_r0 = cy_r0;
			thread_arg->cy_r1 = cy_r1;
			thread_arg->cy_r2 = cy_r2;
			thread_arg->cy_r3 = cy_r3;
			thread_arg->cy_r4 = cy_r4;
			thread_arg->cy_r5 = cy_r5;
			thread_arg->cy_r6 = cy_r6;
			thread_arg->cy_r7 = cy_r7;
			thread_arg->cy_r8 = cy_r8;
			thread_arg->cy_r9 = cy_r9;
			thread_arg->cy_rA = cy_rA;
			thread_arg->cy_rB = cy_rB;
			thread_arg->cy_rC = cy_rC;
			thread_arg->cy_rD = cy_rD;
			thread_arg->cy_rE = cy_rE;
			thread_arg->cy_rF = cy_rF;
		#endif	// SSE2 or AVX?
		}
		else
		{
		#ifdef USE_AVX512
			thread_arg->cy_r0 = cy_r0->d0;		thread_arg->cy_i0 = cy_i0->d0;
			thread_arg->cy_r1 = cy_r0->d1;		thread_arg->cy_i1 = cy_i0->d1;
			thread_arg->cy_r2 = cy_r0->d2;		thread_arg->cy_i2 = cy_i0->d2;
			thread_arg->cy_r3 = cy_r0->d3;		thread_arg->cy_i3 = cy_i0->d3;
			thread_arg->cy_r4 = cy_r0->d4;		thread_arg->cy_i4 = cy_i0->d4;
			thread_arg->cy_r5 = cy_r0->d5;		thread_arg->cy_i5 = cy_i0->d5;
			thread_arg->cy_r6 = cy_r0->d6;		thread_arg->cy_i6 = cy_i0->d6;
			thread_arg->cy_r7 = cy_r0->d7;		thread_arg->cy_i7 = cy_i0->d7;
			thread_arg->cy_r8 = cy_r8->d0;		thread_arg->cy_i8 = cy_i8->d0;
			thread_arg->cy_r9 = cy_r8->d1;		thread_arg->cy_i9 = cy_i8->d1;
			thread_arg->cy_rA = cy_r8->d2;		thread_arg->cy_iA = cy_i8->d2;
			thread_arg->cy_rB = cy_r8->d3;		thread_arg->cy_iB = cy_i8->d3;
			thread_arg->cy_rC = cy_r8->d4;		thread_arg->cy_iC = cy_i8->d4;
			thread_arg->cy_rD = cy_r8->d5;		thread_arg->cy_iD = cy_i8->d5;
			thread_arg->cy_rE = cy_r8->d6;		thread_arg->cy_iE = cy_i8->d6;
			thread_arg->cy_rF = cy_r8->d7;		thread_arg->cy_iF = cy_i8->d7;
			t0 = MAX(max_err->d0,max_err->d1);
			t1 = MAX(max_err->d2,max_err->d3);
			t2 = MAX(max_err->d4,max_err->d5);
			t3 = MAX(max_err->d6,max_err->d7);
			maxerr = MAX( MAX(t0,t1), MAX(t2,t3) );
		#elif defined(USE_AVX)
			thread_arg->cy_r0 = cy_r0->d0;		thread_arg->cy_i0 = cy_i0->d0;
			thread_arg->cy_r1 = cy_r0->d1;		thread_arg->cy_i1 = cy_i0->d1;
			thread_arg->cy_r2 = cy_r0->d2;		thread_arg->cy_i2 = cy_i0->d2;
			thread_arg->cy_r3 = cy_r0->d3;		thread_arg->cy_i3 = cy_i0->d3;
			thread_arg->cy_r4 = cy_r4->d0;		thread_arg->cy_i4 = cy_i4->d0;
			thread_arg->cy_r5 = cy_r4->d1;		thread_arg->cy_i5 = cy_i4->d1;
			thread_arg->cy_r6 = cy_r4->d2;		thread_arg->cy_i6 = cy_i4->d2;
			thread_arg->cy_r7 = cy_r4->d3;		thread_arg->cy_i7 = cy_i4->d3;
			thread_arg->cy_r8 = cy_r8->d0;		thread_arg->cy_i8 = cy_i8->d0;
			thread_arg->cy_r9 = cy_r8->d1;		thread_arg->cy_i9 = cy_i8->d1;
			thread_arg->cy_rA = cy_r8->d2;		thread_arg->cy_iA = cy_i8->d2;
			thread_arg->cy_rB = cy_r8->d3;		thread_arg->cy_iB = cy_i8->d3;
			thread_arg->cy_rC = cy_rC->d0;		thread_arg->cy_iC = cy_iC->d0;
			thread_arg->cy_rD = cy_rC->d1;		thread_arg->cy_iD = cy_iC->d1;
			thread_arg->cy_rE = cy_rC->d2;		thread_arg->cy_iE = cy_iC->d2;
			thread_arg->cy_rF = cy_rC->d3;		thread_arg->cy_iF = cy_iC->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			thread_arg->cy_r0 = cy_r0->d0;		thread_arg->cy_i0 = cy_r0->d1;
			thread_arg->cy_r1 = cy_r2->d0;		thread_arg->cy_i1 = cy_r2->d1;
			thread_arg->cy_r2 = cy_r4->d0;		thread_arg->cy_i2 = cy_r4->d1;
			thread_arg->cy_r3 = cy_r6->d0;		thread_arg->cy_i3 = cy_r6->d1;
			thread_arg->cy_r4 = cy_r8->d0;		thread_arg->cy_i4 = cy_r8->d1;
			thread_arg->cy_r5 = cy_rA->d0;		thread_arg->cy_i5 = cy_rA->d1;
			thread_arg->cy_r6 = cy_rC->d0;		thread_arg->cy_i6 = cy_rC->d1;
			thread_arg->cy_r7 = cy_rE->d0;		thread_arg->cy_i7 = cy_rE->d1;
			thread_arg->cy_r8 = cy_i0->d0;		thread_arg->cy_i8 = cy_i0->d1;
			thread_arg->cy_r9 = cy_i2->d0;		thread_arg->cy_i9 = cy_i2->d1;
			thread_arg->cy_rA = cy_i4->d0;		thread_arg->cy_iA = cy_i4->d1;
			thread_arg->cy_rB = cy_i6->d0;		thread_arg->cy_iB = cy_i6->d1;
			thread_arg->cy_rC = cy_i8->d0;		thread_arg->cy_iC = cy_i8->d1;
			thread_arg->cy_rD = cy_iA->d0;		thread_arg->cy_iD = cy_iA->d1;
			thread_arg->cy_rE = cy_iC->d0;		thread_arg->cy_iE = cy_iC->d1;
			thread_arg->cy_rF = cy_iE->d0;		thread_arg->cy_iF = cy_iE->d1;
			maxerr = MAX(max_err->d0,max_err->d1);
		#else
			thread_arg->cy_r0 = cy_r0;			thread_arg->cy_i0 = cy_i0;
			thread_arg->cy_r1 = cy_r1;			thread_arg->cy_i1 = cy_i1;
			thread_arg->cy_r2 = cy_r2;			thread_arg->cy_i2 = cy_i2;
			thread_arg->cy_r3 = cy_r3;			thread_arg->cy_i3 = cy_i3;
			thread_arg->cy_r4 = cy_r4;			thread_arg->cy_i4 = cy_i4;
			thread_arg->cy_r5 = cy_r5;			thread_arg->cy_i5 = cy_i5;
			thread_arg->cy_r6 = cy_r6;			thread_arg->cy_i6 = cy_i6;
			thread_arg->cy_r7 = cy_r7;			thread_arg->cy_i7 = cy_i7;
			thread_arg->cy_r8 = cy_r8;			thread_arg->cy_i8 = cy_i8;
			thread_arg->cy_r9 = cy_r9;			thread_arg->cy_i9 = cy_i9;
			thread_arg->cy_rA = cy_rA;			thread_arg->cy_iA = cy_iA;
			thread_arg->cy_rB = cy_rB;			thread_arg->cy_iB = cy_iB;
			thread_arg->cy_rC = cy_rC;			thread_arg->cy_iC = cy_iC;
			thread_arg->cy_rD = cy_rD;			thread_arg->cy_iD = cy_iD;
			thread_arg->cy_rE = cy_rE;			thread_arg->cy_iE = cy_iE;
			thread_arg->cy_rF = cy_rF;			thread_arg->cy_iF = cy_iF;
		#endif	// SSE2 or AVX?
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

#ifdef USE_FGT61
	// Carry-macro is adapted from my original circa-2000 F90 mixed-carry subroutine:
//	#define mixed_carry(fx,fy,ix,iy,cy,bjmodn,set,shift)/*\*/
	void mixed_carry(
		const int sw, const int bw, const int nm1, const int bits,
		const int j, const int col, const int co2, const int co3, const int n_minus_sil, const int n_minus_silp1, const int sinwt, const int sinwtm1,
		const double base[], const double baseinv[], const double wt1[], const double wtl, const double wtlp1, const double wtn, const double wtnm1,
		double*maxerr, double*fx, double*fy, uint64*ix, uint64*iy, double*cy, uint32 bjmodn, uint32 set, uint32 shift, uint32 l2_n2)
	{/*@*/
	// Only needed for debuggable function-version:
	const uint64 q  = 0x1FFFFFFFFFFFFFFFull, q2=q+q, q3=q2+q, q4=q2+q2, q5=q4+q, q8=q4+q4;	// q = 2^61 - 1, and needed small multiples
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	int i,m,m2;
	double wt,wtinv,wtA,wtB,wtC,temp,frac;
	uint64 itmp;
		const uint64 two52 = 0x0010000000000000ull;	/* Use qfloat.h-defined TWOE52 here? *//*@*/
		const uint64 qhalf = 0x1000000000000000ull;	/* (q+1)/2 = 2^60 *//*@*/
		uint64 itmp64,ibase,rm,im,fvalmodq,fvlo,fvhi8,qval,ifval,mant,err;/*@*/
		int ishft;	/* Used for inverse-weighting shift count *//*@*/
		int bexp;	/* exponent of floating-double needs to be stored in a signed-int *//*@*/
		double fval = *fx;/*@*/
	/*@*/
		wtA=wt1[col+set];/*@*/
		wtB=wt1[co2-set];/*@*/
		wtC=wt1[co3-set];/*@*/
		wt   =wtl*wtA;/*@*/
		wtinv=wtn*wtB;/*@*/
		/* Branchlessly force bigword-ness for J = 0 by XORing the zero-loop-index flag (!J) with sign bit of (sw - bjmodn), *//*@*/
		i = (!j) ^ ((uint32)(sw - bjmodn) >> 31);	/* which == sw > 0 in this case. Otherwise (!J) == 0, thus XOR = no-op. *//*@*/
		m = ((uint32)(n_minus_sil-bjmodn) >> 31);/*@*/
		m2= 1 + ((uint32)(bjmodn - sinwt) >> 31);/*@*/
		wt   =wt   *one_half[m];/*@*/
		wtinv=wtinv*one_half[m2];/*@*/
		*fx  = *cy + *fx *wtinv;		/* floating (unmodded) part. *//*@*/
		temp = DNINT(*fx );				check_nint(temp, *fx );/*@*/
		frac = fabs(*fx -temp);/*@*/
		if(frac > *maxerr) *maxerr=frac;/*@*/
	/*============ Delay updates of fx and cy until mod-part is done: ============*//*@*/
		ibase = (uint64)base[i];/*@*/
		/* Shift count corr. to wtinv must be modified to include *= n2inv in mod-form, get that by subbing lg(N2): */\
		ishft = 61 - shift - l2_n2;	ishft += (-(ishft < 0)) & 61;/*&*/
		*ix = mul_pow2_modq(*ix,ishft );/*@*/
		qval = qreduce((uint64)*cy + q4 + *ix);	/* mod-q part, partially reduced... *//*@*/
		qval -= q; qval += (-(qval >= qhalf)) & q;	/* ...fully reduce and balance, i.e. put into [-qhalf, (q-1)/2]. *//*@*/
		ifval = *((uint64 *)&fval);				/* Transfer bit pattern of floating datum to integer register for bitwise manipulation... *//*@*/
		mant  = two52 + (ifval & (two52-1));	/* Extract explicit lower 52 bits of mantissa and restore hidden bit, ASSUME normalized input *//*@*/
		bexp  = ( (ifval>>52) & 2047 ) - 1023 - 52 + 3;	/* unbiased binary exponent  (+3 is a premultiply by 8 to make subsequent mod q reduction easier). *//*@*/
		/* whole-number magnitude part % 2^64 ... REPLACE WITH +- RND_TWOE64 after initial testing *//*@*/
		fvlo  = ishft64(mant,bexp);/*@*/
		fvlo >>= 3;/*@*/
		/* (whole-number magnitude part / 2^64)*8 : *//*@*/
		bexp -= 64;/*@*/
		fvhi8 = ishft64(mant,bexp);/*@*/
	/*@*/
		fvalmodq = fvhi8 + fvlo;			/* |fval| % q (unnormalized) *//*@*/
		if(fval < 0) {		/* restore signs *//*@*/
			fvalmodq = -fvalmodq;			/*  fval  % q *//*@*/
			fvlo     = -fvlo;/*@*/
		}/*@*/
		fvalmodq = qreduce( fvalmodq + q4 );		/* make positive and partially reduce... *//*@*/
		fvalmodq -= q; fvalmodq += (-(fvalmodq >= qhalf)) & q;	/* ...balance, i.e. put into [-qhalf, (q-1)/2]. *//*@*/
		/*
		Need both mod-q operands fully reduced, else could have, e.g. fval = q + 3, qval = 3, would give err = q, even though true error = 0.
		Need both mod-q operands balanced, else could have e.g. fval and qval small but of opposite sign,
		i.e. the negative one gets aliased to (x+q) and suddenly taking the difference indicates a big error:
		*/\
		err  = fvalmodq - qval;/*@*/
		/*
		Since fvalmodq is inexact, if qval ~ +-q/2, fvalmodq may get aliased to ~ -+q/2 owing to roundoff,
		and thus giving error ~ +-q. Following hack solves this error aliasing problem, but is slow.
		>>>NEED FAST VERSION<<<
		*/\
	/********* What does user-function kisign - which did not survive in my fragmentary legacy F90 code - do? *********
		err = merge( kisign( ABS((int64)err)-q,qval), err, ABS((int64)err) >= qhalf)
	******************************************************************************************/\
		itmp64   = fvlo - err;						/* Lower 64 bits of exact residue digit, signed, unnormalized. *//*@*/
		rm   = ABS((int64)itmp64) & (ibase - 1);	/* |Current digit| mod current base, i.e. in normalized positive-digit form. *//*@*/
		rm  -= ibase & (-( rm >= (ibase>>1) ));		/* |Current digit| in normalized balanced-digit form, i.e. in [-base/2, +base/2 - 1]. *//*@*/
		rm  -= ( -((int64)itmp64 < 0) & rm ) << 1;	/* restore sign by subbing rm - 2*rm if itmp64 < 0 *//*@*/
		im   = rm + ( q & (-( (int64)rm < 0 )) );	/* make a nonnegative version for the modular transform *//*@*/
	/*@*/
		*ix = mul_pow2_modq( im ,shift );		/* put im into uint64-version and forward weight *//*@*/
	/*@*/
		itmp64   = itmp64 - rm;					/* Lower 64 bits of carry, signed. *//*@*/
		im   = ABS((int64)itmp64) >> bits;/*@*/
		itmp  = fvhi8 << (61-bits);/*@*/
		im  -= ( -((int64)itmp64 < 0) & im ) << 1;	/* restore sign to both im and cy... *//*@*/
		itmp += im - ( ( -(fval < 0) & (uint64)itmp) << 1 );	/* ...and sum to get final carry. *//*@*/
		/* Now, finally can update fx and cy: */\
		*cy   = DNINT(temp*baseinv[i]);	check_nint(*cy, temp*baseinv[i]);/*@*/
		*fx  = (temp-*cy * base[i])*wt;/*@*/
	ASSERT(*fx == (double)rm * wt, "Bad mod-Xout!");	/* put rm into double-version and forward weight *//*@*/
	ASSERT(itmp == *cy, "Bad mod-carry!");
	/*========================*//*@*/
		bjmodn = (bjmodn + bw) & nm1;/*@*/
			wt   =wtlp1*wtA;/*@*/
			wtinv=wtnm1*wtC;/*@*/
			i =((uint32)(sw - bjmodn) >> 31);/*@*/
			m =((uint32)(n_minus_silp1-bjmodn) >> 31);/*@*/
			m2=1 + ((uint32)(bjmodn - sinwtm1) >> 31);/*@*/
			wt   =wt   *one_half[m];/*@*/
			wtinv=wtinv*one_half[m2];/*@*/
			*fy  = *cy+ *fy *wtinv;/*@*/
			temp = DNINT(*fy );				check_nint(temp, *fy );/*@*/
			frac = fabs(*fy -temp);/*@*/
			if(frac > *maxerr) *maxerr=frac;/*@*/
			*cy   = DNINT(temp*baseinv[i]);	check_nint(*cy, temp*baseinv[i]);/*@*/
			*fy  = (temp-*cy * base[i])*wt;/*@*/
		bjmodn = (bjmodn + bw) & nm1;/*@*/
	}
#endif
