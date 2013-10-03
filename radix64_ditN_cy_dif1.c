/*******************************************************************************
*                                                                              *
*   (C) 1997-2013 by Ernst W. Mayer.                                           *
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

#ifdef MULTITHREAD
	#ifndef USE_PTHREAD
		#error Pthreads is only thread model currently supported!
	#endif
#endif

// SIMD code only available for 64-bit GCC build - others simply use scalar DFT macros with SIMD-compatible data layout
#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)

	#define EPS 1e-10

	#include "sse2_macro.h"

  // For Mersenne-mod we need (16 [SSE2] or 64 [AVX]) + 4 added slots for the half_arr lookup tables.
  // For Fermat-mod we use RADIX*4 = 256 [note there is no LOACC option for this power-of-2 DFT] more
  // slots in AVX mode for the compact negacyclic-roots chained-multiply scheme. Add larger of the 2 numbers -
  // 256 for AVX, 20 for SSE2 - to (half_arr_offset64 + RADIX) to get AVX value of radix64_creals_in_local_store:
  #ifdef USE_AVX
	const int half_arr_offset64 = 336;	// + RADIX = 400; Used for thread local-storage-integrity checking
	const int radix64_creals_in_local_store = 656;	// (half_arr_offset64 + RADIX) + 256 and round up to nearest multiple of 4
  #else
	const int half_arr_offset64 = 368;	// + RADIX = 432; Used for thread local-storage-integrity checking
	const int radix64_creals_in_local_store = 452;	// (half_arr_offset64 + RADIX) + 20 and round up to nearest multiple of 4
  #endif

#elif defined(USE_SSE2)

	#error SIMD build only supported for GCC-compatible compiler under *nix/macOS in 64-bit mode!

#endif	/* USE_SSE2 */

#ifdef USE_PTHREAD

	// Use non-pooled simple spawn/rejoin thread-team model
	#include "threadpool.h"

	struct cy_thread_data_t{
	// int data - if needed, pad to yield an even number of these:
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
		int _pad_;

	// double data:
		double maxerr;
		double scale;

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
		int *si;
		struct complex *rn0;
		struct complex *rn1;
		vec_dbl *r00;
		vec_dbl *half_arr;

		int bjmodn00;	int bjmodn20;
		int bjmodn01;	int bjmodn21;
		int bjmodn02;	int bjmodn22;
		int bjmodn03;	int bjmodn23;
		int bjmodn04;	int bjmodn24;
		int bjmodn05;	int bjmodn25;
		int bjmodn06;	int bjmodn26;
		int bjmodn07;	int bjmodn27;
		int bjmodn08;	int bjmodn28;
		int bjmodn09;	int bjmodn29;
		int bjmodn0A;	int bjmodn2A;
		int bjmodn0B;	int bjmodn2B;
		int bjmodn0C;	int bjmodn2C;
		int bjmodn0D;	int bjmodn2D;
		int bjmodn0E;	int bjmodn2E;
		int bjmodn0F;	int bjmodn2F;
		int bjmodn10;	int bjmodn30;
		int bjmodn11;	int bjmodn31;
		int bjmodn12;	int bjmodn32;
		int bjmodn13;	int bjmodn33;
		int bjmodn14;	int bjmodn34;
		int bjmodn15;	int bjmodn35;
		int bjmodn16;	int bjmodn36;
		int bjmodn17;	int bjmodn37;
		int bjmodn18;	int bjmodn38;
		int bjmodn19;	int bjmodn39;
		int bjmodn1A;	int bjmodn3A;
		int bjmodn1B;	int bjmodn3B;
		int bjmodn1C;	int bjmodn3C;
		int bjmodn1D;	int bjmodn3D;
		int bjmodn1E;	int bjmodn3E;
		int bjmodn1F;	int bjmodn3F;
		/* carries: */
		double cy_r00;	double cy_r20;
		double cy_r01;	double cy_r21;
		double cy_r02;	double cy_r22;
		double cy_r03;	double cy_r23;
		double cy_r04;	double cy_r24;
		double cy_r05;	double cy_r25;
		double cy_r06;	double cy_r26;
		double cy_r07;	double cy_r27;
		double cy_r08;	double cy_r28;
		double cy_r09;	double cy_r29;
		double cy_r0A;	double cy_r2A;
		double cy_r0B;	double cy_r2B;
		double cy_r0C;	double cy_r2C;
		double cy_r0D;	double cy_r2D;
		double cy_r0E;	double cy_r2E;
		double cy_r0F;	double cy_r2F;
		double cy_r10;	double cy_r30;
		double cy_r11;	double cy_r31;
		double cy_r12;	double cy_r32;
		double cy_r13;	double cy_r33;
		double cy_r14;	double cy_r34;
		double cy_r15;	double cy_r35;
		double cy_r16;	double cy_r36;
		double cy_r17;	double cy_r37;
		double cy_r18;	double cy_r38;
		double cy_r19;	double cy_r39;
		double cy_r1A;	double cy_r3A;
		double cy_r1B;	double cy_r3B;
		double cy_r1C;	double cy_r3C;
		double cy_r1D;	double cy_r3D;
		double cy_r1E;	double cy_r3E;
		double cy_r1F;	double cy_r3F;

		double cy_i00;	double cy_i20;
		double cy_i01;	double cy_i21;
		double cy_i02;	double cy_i22;
		double cy_i03;	double cy_i23;
		double cy_i04;	double cy_i24;
		double cy_i05;	double cy_i25;
		double cy_i06;	double cy_i26;
		double cy_i07;	double cy_i27;
		double cy_i08;	double cy_i28;
		double cy_i09;	double cy_i29;
		double cy_i0A;	double cy_i2A;
		double cy_i0B;	double cy_i2B;
		double cy_i0C;	double cy_i2C;
		double cy_i0D;	double cy_i2D;
		double cy_i0E;	double cy_i2E;
		double cy_i0F;	double cy_i2F;
		double cy_i10;	double cy_i30;
		double cy_i11;	double cy_i31;
		double cy_i12;	double cy_i32;
		double cy_i13;	double cy_i33;
		double cy_i14;	double cy_i34;
		double cy_i15;	double cy_i35;
		double cy_i16;	double cy_i36;
		double cy_i17;	double cy_i37;
		double cy_i18;	double cy_i38;
		double cy_i19;	double cy_i39;
		double cy_i1A;	double cy_i3A;
		double cy_i1B;	double cy_i3B;
		double cy_i1C;	double cy_i3C;
		double cy_i1D;	double cy_i3D;
		double cy_i1E;	double cy_i3E;
		double cy_i1F;	double cy_i3F;
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
	const uint32 RADIX = 64;
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
#ifdef USE_SSE2
	const int sz_vd = sizeof(vec_dbl);
	// lg(sizeof(vec_dbl)):
  #ifdef USE_AVX
	const int l2_sz_vd = 5;
//	const int half_arr_offset =  90;	// Used for thread local-storage-integrity checking
  #else
	const int l2_sz_vd = 4;
//	const int half_arr_offset = 106;	// Used for thread local-storage-integrity checking
  #endif
#endif
	int NDIVR,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,outer,nbytes;
	int col,co2,co3;
  #if defined(USE_AVX) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
  #endif
	double rt,it,wt_re,wt_im;	/* Fermat-mod weights stuff, used in both scalar and AVX mode */
	static uint32 bjmodnini;
	static uint64 psave=0;
	static uint32 bw,sw,nm1,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30,p38;
	static double radix_inv, n2inv;
#if defined(USE_SSE2) || !defined(MULTITHREAD)
	const double c64_1= 0.99518472667219688624, s64_1 = 0.09801714032956060199,	/* exp(  I*twopi/64) */
				c64_2 = 0.98078528040323044913, s64_2 = 0.19509032201612826785,	/* exp(2*I*twopi/64) */
				c64_3 = 0.95694033573220886494, s64_3 = 0.29028467725446236764,	/* exp(3*I*twopi/64) */
				c64_4 = 0.92387953251128675613, s64_4 = 0.38268343236508977173,	/* exp(4*I*twopi/64) */
				c64_5 = 0.88192126434835502971, s64_5 = 0.47139673682599764856,	/* exp(5*I*twopi/64) */
				c64_6 = 0.83146961230254523708, s64_6 = 0.55557023301960222474,	/* exp(6*I*twopi/64) */
				c64_7 = 0.77301045336273696081, s64_7 = 0.63439328416364549822;	/* exp(7*I*twopi/64) */
#endif
	double scale,
		t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F,
		t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F,
		t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F,
		t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F,
		t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t4A,t4B,t4C,t4D,t4E,t4F,
		t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t5A,t5B,t5C,t5D,t5E,t5F,
		t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t6A,t6B,t6C,t6D,t6E,t6F,
		t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t7A,t7B,t7C,t7D,t7E,t7F;
	double dtmp, maxerr = 0.0;
	int err;
	static int first_entry=TRUE;

	int n_div_nwt;

// SIMD code only available for 64-bit GCC build - others simply use scalar DFT macros with SIMD-compatible data layout
#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)

	int idx_offset,idx_incr;
	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_nm1;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else

//	int i0,i1,m0,m1,m3;	/* m2 already def'd for regular carry sequence */
	double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7;	/* Addresses into array sections */

  #endif

	static int
	 *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn0A,*bjmodn0B,*bjmodn0C,*bjmodn0D,*bjmodn0E,*bjmodn0F
	,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,*bjmodn1A,*bjmodn1B,*bjmodn1C,*bjmodn1D,*bjmodn1E,*bjmodn1F
	,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27,*bjmodn28,*bjmodn29,*bjmodn2A,*bjmodn2B,*bjmodn2C,*bjmodn2D,*bjmodn2E,*bjmodn2F
	,*bjmodn30,*bjmodn31,*bjmodn32,*bjmodn33,*bjmodn34,*bjmodn35,*bjmodn36,*bjmodn37,*bjmodn38,*bjmodn39,*bjmodn3A,*bjmodn3B,*bjmodn3C,*bjmodn3D,*bjmodn3E,*bjmodn3F;
	const double crnd = 3.0*0x4000000*0x2000000;
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	static vec_dbl *max_err, *sse2_rnd, *half_arr, *cc0, *ss0,
		 *isrt2, *cc1, *ss1, *cc2, *ss2, *cc3, *ss3, *cc4, *ss4, *cc5, *ss5, *cc6, *ss6, *cc7, *ss7,
		*nisrt2,*ncc1,*nss1,*ncc2,*nss2,*ncc3,*nss3,*ncc4,*nss4,*ncc5,*nss5,*ncc6,*nss6,*ncc7,*nss7,	// each non-unity root now needs a negated counterpart
		*r00,*r02,*r04,*r06,*r08,*r0A,*r0C,*r0E,
		*r10,*r12,*r14,*r16,*r18,*r1A,*r1C,*r1E,
		*r20,*r22,*r24,*r26,*r28,*r2A,*r2C,*r2E,
		*r30,*r32,*r34,*r36,*r38,*r3A,*r3C,*r3E,
		*r40,*r42,*r44,*r46,*r48,*r4A,*r4C,*r4E,
		*r50,*r52,*r54,*r56,*r58,*r5A,*r5C,*r5E,
		*r60,*r62,*r64,*r66,*r68,*r6A,*r6C,*r6E,
		*r70,*r72,*r74,*r76,*r78,*r7A,*r7C,*r7E,
		*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p0ar,*s1p0br,*s1p0cr,*s1p0dr,*s1p0er,*s1p0fr,
		*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p1ar,*s1p1br,*s1p1cr,*s1p1dr,*s1p1er,*s1p1fr,
		*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r,*s1p28r,*s1p29r,*s1p2ar,*s1p2br,*s1p2cr,*s1p2dr,*s1p2er,*s1p2fr,
		*s1p30r,*s1p31r,*s1p32r,*s1p33r,*s1p34r,*s1p35r,*s1p36r,*s1p37r,*s1p38r,*s1p39r,*s1p3ar,*s1p3br,*s1p3cr,*s1p3dr,*s1p3er,*s1p3fr,
		*cy_r00,*cy_r04,*cy_r08,*cy_r0C,*cy_r10,*cy_r14,*cy_r18,*cy_r1C,*cy_r20,*cy_r24,*cy_r28,*cy_r2C,*cy_r30,*cy_r34,*cy_r38,*cy_r3C,
		*cy_i00,*cy_i04,*cy_i08,*cy_i0C,*cy_i10,*cy_i14,*cy_i18,*cy_i1C,*cy_i20,*cy_i24,*cy_i28,*cy_i2C,*cy_i30,*cy_i34,*cy_i38,*cy_i3C;
  #ifndef USE_AVX
	static vec_dbl
		*cy_r02,*cy_r06,*cy_r0A,*cy_r0E,*cy_r12,*cy_r16,*cy_r1A,*cy_r1E,*cy_r22,*cy_r26,*cy_r2A,*cy_r2E,*cy_r32,*cy_r36,*cy_r3A,*cy_r3E,
		*cy_i02,*cy_i06,*cy_i0A,*cy_i0E,*cy_i12,*cy_i16,*cy_i1A,*cy_i1E,*cy_i22,*cy_i26,*cy_i2A,*cy_i2E,*cy_i32,*cy_i36,*cy_i3A,*cy_i3E;
  #else
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
	static task_control_t   task_control = {NULL, (void*)cy64_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	int k1,k2,m,m2,ntmp;
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
  #if PFETCH
	double *addr;
	int prefetch_offset;
  #endif
	int  bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn0A,bjmodn0B,bjmodn0C,bjmodn0D,bjmodn0E,bjmodn0F
		,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn1A,bjmodn1B,bjmodn1C,bjmodn1D,bjmodn1E,bjmodn1F
		,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,bjmodn28,bjmodn29,bjmodn2A,bjmodn2B,bjmodn2C,bjmodn2D,bjmodn2E,bjmodn2F
		,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,bjmodn36,bjmodn37,bjmodn38,bjmodn39,bjmodn3A,bjmodn3B,bjmodn3C,bjmodn3D,bjmodn3E,bjmodn3F;
	double temp,frac
		,ar_p00,ar_p01,ar_p02,ar_p03,ar_p04,ar_p05,ar_p06,ar_p07,ar_p08,ar_p09,ar_p0a,ar_p0b,ar_p0c,ar_p0d,ar_p0e,ar_p0f
		,ar_p10,ar_p11,ar_p12,ar_p13,ar_p14,ar_p15,ar_p16,ar_p17,ar_p18,ar_p19,ar_p1a,ar_p1b,ar_p1c,ar_p1d,ar_p1e,ar_p1f
		,ar_p20,ar_p21,ar_p22,ar_p23,ar_p24,ar_p25,ar_p26,ar_p27,ar_p28,ar_p29,ar_p2a,ar_p2b,ar_p2c,ar_p2d,ar_p2e,ar_p2f
		,ar_p30,ar_p31,ar_p32,ar_p33,ar_p34,ar_p35,ar_p36,ar_p37,ar_p38,ar_p39,ar_p3a,ar_p3b,ar_p3c,ar_p3d,ar_p3e,ar_p3f
		,ai_p00,ai_p01,ai_p02,ai_p03,ai_p04,ai_p05,ai_p06,ai_p07,ai_p08,ai_p09,ai_p0a,ai_p0b,ai_p0c,ai_p0d,ai_p0e,ai_p0f
		,ai_p10,ai_p11,ai_p12,ai_p13,ai_p14,ai_p15,ai_p16,ai_p17,ai_p18,ai_p19,ai_p1a,ai_p1b,ai_p1c,ai_p1d,ai_p1e,ai_p1f
		,ai_p20,ai_p21,ai_p22,ai_p23,ai_p24,ai_p25,ai_p26,ai_p27,ai_p28,ai_p29,ai_p2a,ai_p2b,ai_p2c,ai_p2d,ai_p2e,ai_p2f
		,ai_p30,ai_p31,ai_p32,ai_p33,ai_p34,ai_p35,ai_p36,ai_p37,ai_p38,ai_p39,ai_p3a,ai_p3b,ai_p3c,ai_p3d,ai_p3e,ai_p3f
		,cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r0A,cy_r0B,cy_r0C,cy_r0D,cy_r0E,cy_r0F
		,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r1A,cy_r1B,cy_r1C,cy_r1D,cy_r1E,cy_r1F
		,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27,cy_r28,cy_r29,cy_r2A,cy_r2B,cy_r2C,cy_r2D,cy_r2E,cy_r2F
		,cy_r30,cy_r31,cy_r32,cy_r33,cy_r34,cy_r35,cy_r36,cy_r37,cy_r38,cy_r39,cy_r3A,cy_r3B,cy_r3C,cy_r3D,cy_r3E,cy_r3F
		,cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i0A,cy_i0B,cy_i0C,cy_i0D,cy_i0E,cy_i0F
		,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i1A,cy_i1B,cy_i1C,cy_i1D,cy_i1E,cy_i1F
		,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27,cy_i28,cy_i29,cy_i2A,cy_i2B,cy_i2C,cy_i2D,cy_i2E,cy_i2F
		,cy_i30,cy_i31,cy_i32,cy_i33,cy_i34,cy_i35,cy_i36,cy_i37,cy_i38,cy_i39,cy_i3A,cy_i3B,cy_i3C,cy_i3D,cy_i3E,cy_i3F;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int
	 *_bjmodn00 = 0x0,*_bjmodn01 = 0x0,*_bjmodn02 = 0x0,*_bjmodn03 = 0x0,*_bjmodn04 = 0x0,*_bjmodn05 = 0x0,*_bjmodn06 = 0x0,*_bjmodn07 = 0x0,*_bjmodn08 = 0x0,*_bjmodn09 = 0x0,*_bjmodn0A = 0x0,*_bjmodn0B = 0x0,*_bjmodn0C = 0x0,*_bjmodn0D = 0x0,*_bjmodn0E = 0x0,*_bjmodn0F = 0x0
	,*_bjmodn10 = 0x0,*_bjmodn11 = 0x0,*_bjmodn12 = 0x0,*_bjmodn13 = 0x0,*_bjmodn14 = 0x0,*_bjmodn15 = 0x0,*_bjmodn16 = 0x0,*_bjmodn17 = 0x0,*_bjmodn18 = 0x0,*_bjmodn19 = 0x0,*_bjmodn1A = 0x0,*_bjmodn1B = 0x0,*_bjmodn1C = 0x0,*_bjmodn1D = 0x0,*_bjmodn1E = 0x0,*_bjmodn1F = 0x0
	,*_bjmodn20 = 0x0,*_bjmodn21 = 0x0,*_bjmodn22 = 0x0,*_bjmodn23 = 0x0,*_bjmodn24 = 0x0,*_bjmodn25 = 0x0,*_bjmodn26 = 0x0,*_bjmodn27 = 0x0,*_bjmodn28 = 0x0,*_bjmodn29 = 0x0,*_bjmodn2A = 0x0,*_bjmodn2B = 0x0,*_bjmodn2C = 0x0,*_bjmodn2D = 0x0,*_bjmodn2E = 0x0,*_bjmodn2F = 0x0
	,*_bjmodn30 = 0x0,*_bjmodn31 = 0x0,*_bjmodn32 = 0x0,*_bjmodn33 = 0x0,*_bjmodn34 = 0x0,*_bjmodn35 = 0x0,*_bjmodn36 = 0x0,*_bjmodn37 = 0x0,*_bjmodn38 = 0x0,*_bjmodn39 = 0x0,*_bjmodn3A = 0x0,*_bjmodn3B = 0x0,*_bjmodn3C = 0x0,*_bjmodn3D = 0x0,*_bjmodn3E = 0x0,*_bjmodn3F = 0x0;
	static int *_bjmodnini = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_maxerr = 0x0
	,*_cy_r00 = 0x0,*_cy_r01 = 0x0,*_cy_r02 = 0x0,*_cy_r03 = 0x0,*_cy_r04 = 0x0,*_cy_r05 = 0x0,*_cy_r06 = 0x0,*_cy_r07 = 0x0,*_cy_r08 = 0x0,*_cy_r09 = 0x0,*_cy_r0A = 0x0,*_cy_r0B = 0x0,*_cy_r0C = 0x0,*_cy_r0D = 0x0,*_cy_r0E = 0x0,*_cy_r0F = 0x0
	,*_cy_r10 = 0x0,*_cy_r11 = 0x0,*_cy_r12 = 0x0,*_cy_r13 = 0x0,*_cy_r14 = 0x0,*_cy_r15 = 0x0,*_cy_r16 = 0x0,*_cy_r17 = 0x0,*_cy_r18 = 0x0,*_cy_r19 = 0x0,*_cy_r1A = 0x0,*_cy_r1B = 0x0,*_cy_r1C = 0x0,*_cy_r1D = 0x0,*_cy_r1E = 0x0,*_cy_r1F = 0x0
	,*_cy_r20 = 0x0,*_cy_r21 = 0x0,*_cy_r22 = 0x0,*_cy_r23 = 0x0,*_cy_r24 = 0x0,*_cy_r25 = 0x0,*_cy_r26 = 0x0,*_cy_r27 = 0x0,*_cy_r28 = 0x0,*_cy_r29 = 0x0,*_cy_r2A = 0x0,*_cy_r2B = 0x0,*_cy_r2C = 0x0,*_cy_r2D = 0x0,*_cy_r2E = 0x0,*_cy_r2F = 0x0
	,*_cy_r30 = 0x0,*_cy_r31 = 0x0,*_cy_r32 = 0x0,*_cy_r33 = 0x0,*_cy_r34 = 0x0,*_cy_r35 = 0x0,*_cy_r36 = 0x0,*_cy_r37 = 0x0,*_cy_r38 = 0x0,*_cy_r39 = 0x0,*_cy_r3A = 0x0,*_cy_r3B = 0x0,*_cy_r3C = 0x0,*_cy_r3D = 0x0,*_cy_r3E = 0x0,*_cy_r3F = 0x0
	,*_cy_i00 = 0x0,*_cy_i01 = 0x0,*_cy_i02 = 0x0,*_cy_i03 = 0x0,*_cy_i04 = 0x0,*_cy_i05 = 0x0,*_cy_i06 = 0x0,*_cy_i07 = 0x0,*_cy_i08 = 0x0,*_cy_i09 = 0x0,*_cy_i0A = 0x0,*_cy_i0B = 0x0,*_cy_i0C = 0x0,*_cy_i0D = 0x0,*_cy_i0E = 0x0,*_cy_i0F = 0x0
	,*_cy_i10 = 0x0,*_cy_i11 = 0x0,*_cy_i12 = 0x0,*_cy_i13 = 0x0,*_cy_i14 = 0x0,*_cy_i15 = 0x0,*_cy_i16 = 0x0,*_cy_i17 = 0x0,*_cy_i18 = 0x0,*_cy_i19 = 0x0,*_cy_i1A = 0x0,*_cy_i1B = 0x0,*_cy_i1C = 0x0,*_cy_i1D = 0x0,*_cy_i1E = 0x0,*_cy_i1F = 0x0
	,*_cy_i20 = 0x0,*_cy_i21 = 0x0,*_cy_i22 = 0x0,*_cy_i23 = 0x0,*_cy_i24 = 0x0,*_cy_i25 = 0x0,*_cy_i26 = 0x0,*_cy_i27 = 0x0,*_cy_i28 = 0x0,*_cy_i29 = 0x0,*_cy_i2A = 0x0,*_cy_i2B = 0x0,*_cy_i2C = 0x0,*_cy_i2D = 0x0,*_cy_i2E = 0x0,*_cy_i2F = 0x0
	,*_cy_i30 = 0x0,*_cy_i31 = 0x0,*_cy_i32 = 0x0,*_cy_i33 = 0x0,*_cy_i34 = 0x0,*_cy_i35 = 0x0,*_cy_i36 = 0x0,*_cy_i37 = 0x0,*_cy_i38 = 0x0,*_cy_i39 = 0x0,*_cy_i3A = 0x0,*_cy_i3B = 0x0,*_cy_i3C = 0x0,*_cy_i3D = 0x0,*_cy_i3E = 0x0,*_cy_i3F = 0x0;

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/RADIX;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/RADIX in radix64_ditN_cy_dif1.\n",iter);
		if(INTERACT)fprintf(stderr,"%s",cbuf);
		fp = fopen(   OFILE,"a");
		fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;
		err=ERR_CARRY;
		return(err);
	}

	if(p != psave)
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry: */

	if(first_entry)
	{
		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)RADIX));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

		nm1   = n-1;

	/* to-do: Add threading to sse2 code */
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
		ASSERT(HERE, CY_THREADS >= NTHREADS,"CY_THREADS < NTHREADS");
		ASSERT(HERE, isPow2(CY_THREADS)    ,"CY_THREADS not a power of 2!");
		if(CY_THREADS > 1)
		{
			ASSERT(HERE, NDIVR    %CY_THREADS == 0,"NDIVR    %CY_THREADS != 0");
			ASSERT(HERE, n_div_nwt%CY_THREADS == 0,"n_div_nwt%CY_THREADS != 0");
		}

	  #ifdef USE_PTHREAD

		j = (uint32)sizeof(struct cy_thread_data_t);
		if(0 != (j & 0xf)) {
			printf("sizeof(cy_thread_data_t) = %x\n",j);
			ASSERT(HERE, 0, "struct cy_thread_data_t not 16-byte size multiple!");
		}
		tdat = (struct cy_thread_data_t *)calloc(CY_THREADS, sizeof(struct cy_thread_data_t));

		// MacOS does weird things with threading (e.g. Idle" main thread burning 100% of 1 CPU)
		// so on that platform try to be clever and interleave main-thread and threadpool-work processing
		#ifdef OS_TYPE_MACOSX

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

			main_work_units = 0;
			pool_work_units = CY_THREADS;
			ASSERT(HERE, 0x0 != (tpool = threadpool_init(CY_THREADS, MAX_THREADS, CY_THREADS, &thread_control)), "threadpool_init failed!");

		#endif

	  #endif

		fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);

	#else
		CY_THREADS = 1;
	#endif

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
		// int data:
			tdat[ithread].tid = ithread;
			tdat[ithread].ndivr = NDIVR;

			tdat[ithread].sw  = sw;
			tdat[ithread].nwt = nwt;

		// pointer data:
			tdat[ithread].arrdat = a;			/* Main data array */
			tdat[ithread].wt0 = wt0;
			tdat[ithread].wt1 = wt1;
			tdat[ithread].si  = si;
			tdat[ithread].rn0 = rn0;
			tdat[ithread].rn1 = rn1;
		}
	#endif

	// SIMD code only available for 64-bit GCC build - others simply use scalar DFT macros with SIMD-compatible data layout
	#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)

		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use vector-double type size (16 bytes for SSE2, 32 for AVX) to alloc a block of local storage
		// consisting of 128 vec_dbl and ([8 if SSE2, 16 if AVX] + RADIX/2) uint64 element slots per thread
		cslots_in_local_store = radix64_creals_in_local_store + (20+RADIX/2)/2;	// Just add enough int64 space for both cases, plus some
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix64_creals_in_local_store);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 64 vec_ddl-sized slots of sc_arr for temporaries, next 7 for the nontrivial complex 16th roots,
	next 32 for the vector carries, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
		tmp = sc_ptr;
		r00 = tmp + 0x00;		r40 = tmp + 0x40;
		r02 = tmp + 0x02;		r42 = tmp + 0x42;
		r04 = tmp + 0x04;		r44 = tmp + 0x44;
		r06 = tmp + 0x06;		r46 = tmp + 0x46;
		r08 = tmp + 0x08;		r48 = tmp + 0x48;
		r0A = tmp + 0x0a;		r4A = tmp + 0x4a;
		r0C = tmp + 0x0c;		r4C = tmp + 0x4c;
		r0E = tmp + 0x0e;		r4E = tmp + 0x4e;
		r10 = tmp + 0x10;		r50 = tmp + 0x50;
		r12 = tmp + 0x12;		r52 = tmp + 0x52;
		r14 = tmp + 0x14;		r54 = tmp + 0x54;
		r16 = tmp + 0x16;		r56 = tmp + 0x56;
		r18 = tmp + 0x18;		r58 = tmp + 0x58;
		r1A = tmp + 0x1a;		r5A = tmp + 0x5a;
		r1C = tmp + 0x1c;		r5C = tmp + 0x5c;
		r1E = tmp + 0x1e;		r5E = tmp + 0x5e;
		r20 = tmp + 0x20;		r60 = tmp + 0x60;
		r22 = tmp + 0x22;		r62 = tmp + 0x62;
		r24 = tmp + 0x24;		r64 = tmp + 0x64;
		r26 = tmp + 0x26;		r66 = tmp + 0x66;
		r28 = tmp + 0x28;		r68 = tmp + 0x68;
		r2A = tmp + 0x2a;		r6A = tmp + 0x6a;
		r2C = tmp + 0x2c;		r6C = tmp + 0x6c;
		r2E = tmp + 0x2e;		r6E = tmp + 0x6e;
		r30 = tmp + 0x30;		r70 = tmp + 0x70;
		r32 = tmp + 0x32;		r72 = tmp + 0x72;
		r34 = tmp + 0x34;		r74 = tmp + 0x74;
		r36 = tmp + 0x36;		r76 = tmp + 0x76;
		r38 = tmp + 0x38;		r78 = tmp + 0x78;
		r3A = tmp + 0x3a;		r7A = tmp + 0x7a;
		r3C = tmp + 0x3c;		r7C = tmp + 0x7c;
		r3E = tmp + 0x3e;		r7E = tmp + 0x7e;
		tmp += 0x80;
		s1p00r = tmp + 0x00;	s1p20r = tmp + 0x40;
		s1p01r = tmp + 0x02;	s1p21r = tmp + 0x42;
		s1p02r = tmp + 0x04;	s1p22r = tmp + 0x44;
		s1p03r = tmp + 0x06;	s1p23r = tmp + 0x46;
		s1p04r = tmp + 0x08;	s1p24r = tmp + 0x48;
		s1p05r = tmp + 0x0a;	s1p25r = tmp + 0x4a;
		s1p06r = tmp + 0x0c;	s1p26r = tmp + 0x4c;
		s1p07r = tmp + 0x0e;	s1p27r = tmp + 0x4e;
		s1p08r = tmp + 0x10;	s1p28r = tmp + 0x50;
		s1p09r = tmp + 0x12;	s1p29r = tmp + 0x52;
		s1p0ar = tmp + 0x14;	s1p2ar = tmp + 0x54;
		s1p0br = tmp + 0x16;	s1p2br = tmp + 0x56;
		s1p0cr = tmp + 0x18;	s1p2cr = tmp + 0x58;
		s1p0dr = tmp + 0x1a;	s1p2dr = tmp + 0x5a;
		s1p0er = tmp + 0x1c;	s1p2er = tmp + 0x5c;
		s1p0fr = tmp + 0x1e;	s1p2fr = tmp + 0x5e;
		s1p10r = tmp + 0x20;	s1p30r = tmp + 0x60;
		s1p11r = tmp + 0x22;	s1p31r = tmp + 0x62;
		s1p12r = tmp + 0x24;	s1p32r = tmp + 0x64;
		s1p13r = tmp + 0x26;	s1p33r = tmp + 0x66;
		s1p14r = tmp + 0x28;	s1p34r = tmp + 0x68;
		s1p15r = tmp + 0x2a;	s1p35r = tmp + 0x6a;
		s1p16r = tmp + 0x2c;	s1p36r = tmp + 0x6c;
		s1p17r = tmp + 0x2e;	s1p37r = tmp + 0x6e;
		s1p18r = tmp + 0x30;	s1p38r = tmp + 0x70;
		s1p19r = tmp + 0x32;	s1p39r = tmp + 0x72;
		s1p1ar = tmp + 0x34;	s1p3ar = tmp + 0x74;
		s1p1br = tmp + 0x36;	s1p3br = tmp + 0x76;
		s1p1cr = tmp + 0x38;	s1p3cr = tmp + 0x78;
		s1p1dr = tmp + 0x3a;	s1p3dr = tmp + 0x7a;
		s1p1er = tmp + 0x3c;	s1p3er = tmp + 0x7c;
		s1p1fr = tmp + 0x3e;	s1p3fr = tmp + 0x7e;
		tmp += 0x80;
		// Each non-unity root now needs a negated counterpart:
		/* Stupidity: Since a truly general-purpose [in the sense that it can be used for our radix-64 internal-twiddles]
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
		tmp += 0x2e;
	  #ifdef USE_AVX
		cy_r00	= tmp + 0x00;
		cy_r04	= tmp + 0x01;
		cy_r08	= tmp + 0x02;
		cy_r0C	= tmp + 0x03;
		cy_r10	= tmp + 0x04;
		cy_r14	= tmp + 0x05;
		cy_r18	= tmp + 0x06;
		cy_r1C	= tmp + 0x07;
		cy_r20	= tmp + 0x08;
		cy_r24	= tmp + 0x09;
		cy_r28	= tmp + 0x0a;
		cy_r2C	= tmp + 0x0b;
		cy_r30	= tmp + 0x0c;
		cy_r34	= tmp + 0x0d;
		cy_r38	= tmp + 0x0e;
		cy_r3C	= tmp + 0x0f;
		cy_i00	= tmp + 0x10;
		cy_i04	= tmp + 0x11;
		cy_i08	= tmp + 0x12;
		cy_i0C	= tmp + 0x13;
		cy_i10	= tmp + 0x14;
		cy_i14	= tmp + 0x15;
		cy_i18	= tmp + 0x16;
		cy_i1C	= tmp + 0x17;
		cy_i20	= tmp + 0x18;
		cy_i24	= tmp + 0x19;
		cy_i28	= tmp + 0x1a;
		cy_i2C	= tmp + 0x1b;
		cy_i30	= tmp + 0x1c;
		cy_i34	= tmp + 0x1d;
		cy_i38	= tmp + 0x1e;
		cy_i3C	= tmp + 0x1f;
		tmp += 0x20;
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 336 vec_dbl
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 68 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy_r00	= tmp + 0x00;		cy_i00	= tmp + 0x20;
		cy_r02	= tmp + 0x01;		cy_i02	= tmp + 0x21;
		cy_r04	= tmp + 0x02;		cy_i04	= tmp + 0x22;
		cy_r06	= tmp + 0x03;		cy_i06	= tmp + 0x23;
		cy_r08	= tmp + 0x04;		cy_i08	= tmp + 0x24;
		cy_r0A	= tmp + 0x05;		cy_i0A	= tmp + 0x25;
		cy_r0C	= tmp + 0x06;		cy_i0C	= tmp + 0x26;
		cy_r0E	= tmp + 0x07;		cy_i0E	= tmp + 0x27;
		cy_r10	= tmp + 0x08;		cy_i10	= tmp + 0x28;
		cy_r12	= tmp + 0x09;		cy_i12	= tmp + 0x29;
		cy_r14	= tmp + 0x0a;		cy_i14	= tmp + 0x2a;
		cy_r16	= tmp + 0x0b;		cy_i16	= tmp + 0x2b;
		cy_r18	= tmp + 0x0c;		cy_i18	= tmp + 0x2c;
		cy_r1A	= tmp + 0x0d;		cy_i1A	= tmp + 0x2d;
		cy_r1C	= tmp + 0x0e;		cy_i1C	= tmp + 0x2e;
		cy_r1E	= tmp + 0x0f;		cy_i1E	= tmp + 0x2f;
		cy_r20	= tmp + 0x10;		cy_i20	= tmp + 0x30;
		cy_r22	= tmp + 0x11;		cy_i22	= tmp + 0x31;
		cy_r24	= tmp + 0x12;		cy_i24	= tmp + 0x32;
		cy_r26	= tmp + 0x13;		cy_i26	= tmp + 0x33;
		cy_r28	= tmp + 0x14;		cy_i28	= tmp + 0x34;
		cy_r2A	= tmp + 0x15;		cy_i2A	= tmp + 0x35;
		cy_r2C	= tmp + 0x16;		cy_i2C	= tmp + 0x36;
		cy_r2E	= tmp + 0x17;		cy_i2E	= tmp + 0x37;
		cy_r30	= tmp + 0x18;		cy_i30	= tmp + 0x38;
		cy_r32	= tmp + 0x19;		cy_i32	= tmp + 0x39;
		cy_r34	= tmp + 0x1a;		cy_i34	= tmp + 0x3a;
		cy_r36	= tmp + 0x1b;		cy_i36	= tmp + 0x3b;
		cy_r38	= tmp + 0x1c;		cy_i38	= tmp + 0x3c;
		cy_r3A	= tmp + 0x1d;		cy_i3A	= tmp + 0x3d;
		cy_r3C	= tmp + 0x1e;		cy_i3C	= tmp + 0x3e;
		cy_r3E	= tmp + 0x1f;		cy_i3E	= tmp + 0x3f;
		tmp += 0x40;
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 368 complex
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, 2 for Fermat-mod */
	  #endif
//		ASSERT(HERE, half_arr_offset == (uint32)(half_arr-sc_ptr), "half_arr_offset mismatches actual!");
		ASSERT(HERE, (radix64_creals_in_local_store << l2_sz_vd) >= ((long)half_arr - (long)r00) + (20 << l2_sz_vd), "radix64_creals_in_local_store checksum failed!");

		/* These remain fixed: */
		VEC_DBL_INIT(nisrt2,-ISRT2);
		VEC_DBL_INIT( isrt2, ISRT2);									// Copies of +ISRT2 needed for 30-asm-macro-operand-GCC-limit workaround:
		VEC_DBL_INIT( cc0,   1.0);		VEC_DBL_INIT( ss0,   0.0);		tmp =  cc0-1; ASSERT(HERE, tmp->d0 == ISRT2 && tmp->d1 == ISRT2, "tmp->d0,1 != ISRT2");
		VEC_DBL_INIT( cc1, c64_1);		VEC_DBL_INIT( ss1, s64_1);		tmp =  cc1-1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc2, c64_2);		VEC_DBL_INIT( ss2, s64_2);		tmp =  cc2-1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc3, c64_3);		VEC_DBL_INIT( ss3, s64_3);		tmp =  cc3-1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc4, c64_4);		VEC_DBL_INIT( ss4, s64_4);		tmp =  cc4-1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc5, c64_5);		VEC_DBL_INIT( ss5, s64_5);		tmp =  cc5-1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc6, c64_6);		VEC_DBL_INIT( ss6, s64_6);		tmp =  cc6-1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT( cc7, c64_7);		VEC_DBL_INIT( ss7, s64_7);		tmp =  cc7-1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc1,-c64_1);		VEC_DBL_INIT(nss1,-s64_1);		tmp = ncc1-1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc2,-c64_2);		VEC_DBL_INIT(nss2,-s64_2);		tmp = ncc2-1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc3,-c64_3);		VEC_DBL_INIT(nss3,-s64_3);		tmp = ncc3-1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc4,-c64_4);		VEC_DBL_INIT(nss4,-s64_4);		tmp = ncc4-1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc5,-c64_5);		VEC_DBL_INIT(nss5,-s64_5);		tmp = ncc5-1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc6,-c64_6);		VEC_DBL_INIT(nss6,-s64_6);		tmp = ncc6-1; VEC_DBL_INIT(tmp, ISRT2);
		VEC_DBL_INIT(ncc7,-c64_7);		VEC_DBL_INIT(nss7,-s64_7);		tmp = ncc7-1; VEC_DBL_INIT(tmp, ISRT2);

		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
		VEC_DBL_INIT(sse2_rnd, crnd);

		// Propagate the above consts to the remaining threads:
		nbytes = (int)nss7 - (int)nisrt2 + sz_vd;	// #bytes in 1st of above block of consts
		tmp = nisrt2;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

		nbytes = sz_vd;	// #bytes in 2nd block of consts, which contains just sse2_rnd
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
			nbytes = 2 << l2_sz_vd;
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
			/* Simple qfloat-based loop to crank out the roots: */
				struct qfloat qc,qs,qx,qy,qt,qn,qmul;
				qt = qfdiv(QPIHALF, i64_to_q((int64)RADIX));	// Pi/2/RADIX
				qc = qfcos(qt);	qs = qfsin(qt);
				qx = QONE;		qy = QZRO;
				for(j = 1; j <= RADIX; j++) {
					// Up-multiply the complex exponential:
					qn = qfmul(qx, qc); qt = qfmul(qy, qs); qmul = qfsub(qn, qt);	// Store qxnew in qmul for now.
					qn = qfmul(qx, qs); qt = qfmul(qy, qc); qy   = qfadd(qn, qt); qx = qmul;
					printf("j = %3u: cos = 0x%16llX\n",j,qfdbl_as_uint64(qx));
				}
				exit(0);
			#endif

			tmp = base_negacyclic_root + RADIX*2;	// First 128 = 2*RADIX slots reserved for RADIX/4 copies of the Re/Im parts of the 4 base multipliers
			tm2 = tmp + RADIX/2 - 1;
											tmp->d0 = 1.0;	(tmp+1)->d0 = 0.0;
			tmp64 = 0x3FEFFD886084CD0Dull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(0x01*I*Pi/128) */
			tmp64 = 0x3FEFF621E3796D7Eull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(0x02*I*Pi/128) */
			tmp64 = 0x3FEFE9CDAD01883Aull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(0x03*I*Pi/128) */	tmp += 2;
			tmp64 = 0x3FEFD88DA3D12526ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(0x04*I*Pi/128) */	tm2 -= 2;
			tmp64 = 0x3FEFC26470E19FD3ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(0x05*I*Pi/128) */
			tmp64 = 0x3FEFA7557F08A517ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(0x06*I*Pi/128) */
			tmp64 = 0x3FEF8764FA714BA9ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(0x07*I*Pi/128) */	tmp += 2;
			tmp64 = 0x3FEF6297CFF75CB0ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(0x08*I*Pi/128) */	tm2 -= 2;
			tmp64 = 0x3FEF38F3AC64E589ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(0x09*I*Pi/128) */
			tmp64 = 0x3FEF0A7EFB9230D7ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(0x0A*I*Pi/128) */
			tmp64 = 0x3FEED740E7684963ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(0x0B*I*Pi/128) */	tmp += 2;
			tmp64 = 0x3FEE9F4156C62DDAull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(0x0C*I*Pi/128) */	tm2 -= 2;
			tmp64 = 0x3FEE6288EC48E112ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(0x0D*I*Pi/128) */
			tmp64 = 0x3FEE212104F686E5ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(0x0E*I*Pi/128) */
			tmp64 = 0x3FEDDB13B6CCC23Cull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(0x0F*I*Pi/128) */	tmp += 2;
			tmp64 = 0x3FED906BCF328D46ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(0x10*I*Pi/128) */	tm2 -= 2;
			tmp64 = 0x3FED4134D14DC93Aull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(0x11*I*Pi/128) */
			tmp64 = 0x3FECED7AF43CC773ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(0x12*I*Pi/128) */
			tmp64 = 0x3FEC954B213411F5ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(0x13*I*Pi/128) */	tmp += 2;
			tmp64 = 0x3FEC38B2F180BDB1ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(0x14*I*Pi/128) */	tm2 -= 2;
			tmp64 = 0x3FEBD7C0AC6F952Aull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(0x15*I*Pi/128) */
			tmp64 = 0x3FEB728345196E3Eull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(0x16*I*Pi/128) */
			tmp64 = 0x3FEB090A58150200ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(0x17*I*Pi/128) */	tmp += 2;
			tmp64 = 0x3FEA9B66290EA1A3ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(0x18*I*Pi/128) */	tm2 -= 2;
			tmp64 = 0x3FEA29A7A0462782ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(0x19*I*Pi/128) */
			tmp64 = 0x3FE9B3E047F38741ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(0x1A*I*Pi/128) */
			tmp64 = 0x3FE93A22499263FBull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(0x1B*I*Pi/128) */	tmp += 2;
			tmp64 = 0x3FE8BC806B151741ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(0x1C*I*Pi/128) */	tm2 -= 2;
			tmp64 = 0x3FE83B0E0BFF976Eull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(0x1D*I*Pi/128) */
			tmp64 = 0x3FE7B5DF226AAFAFull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(0x1E*I*Pi/128) */
			tmp64 = 0x3FE72D0837EFFF96ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(0x1F*I*Pi/128) */	tmp += 2;
			tmp64 = 0x3FE6A09E667F3BCDull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(0x20*I*Pi/128) */	tm2 -= 2;
			tmp64 = 0x3FE610B7551D2CDFull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(0x21*I*Pi/128) */
			tmp64 = 0x3FE57D69348CECA0ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(0x22*I*Pi/128) */
			tmp64 = 0x3FE4E6CABBE3E5E9ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(0x23*I*Pi/128) */	tmp += 2;
			tmp64 = 0x3FE44CF325091DD6ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(0x24*I*Pi/128) */	tm2 -= 2;
			tmp64 = 0x3FE3AFFA292050B9ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(0x25*I*Pi/128) */
			tmp64 = 0x3FE30FF7FCE17035ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(0x26*I*Pi/128) */
			tmp64 = 0x3FE26D054CDD12DFull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(0x27*I*Pi/128) */	tmp += 2;
			tmp64 = 0x3FE1C73B39AE68C8ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(0x28*I*Pi/128) */	tm2 -= 2;
			tmp64 = 0x3FE11EB3541B4B23ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(0x29*I*Pi/128) */
			tmp64 = 0x3FE073879922FFEEull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(0x2A*I*Pi/128) */
			tmp64 = 0x3FDF8BA4DBF89ABAull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(0x2B*I*Pi/128) */	tmp += 2;
			tmp64 = 0x3FDE2B5D3806F63Bull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(0x2C*I*Pi/128) */	tm2 -= 2;
			tmp64 = 0x3FDCC66E9931C45Eull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(0x2D*I*Pi/128) */
			tmp64 = 0x3FDB5D1009E15CC0ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(0x2E*I*Pi/128) */
			tmp64 = 0x3FD9EF7943A8ED8Aull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(0x2F*I*Pi/128) */	tmp += 2;
			tmp64 = 0x3FD87DE2A6AEA963ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(0x30*I*Pi/128) */	tm2 -= 2;
			tmp64 = 0x3FD7088530FA459Full;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(0x31*I*Pi/128) */
			tmp64 = 0x3FD58F9A75AB1FDDull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(0x32*I*Pi/128) */
			tmp64 = 0x3FD4135C94176601ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(0x33*I*Pi/128) */	tmp += 2;
			tmp64 = 0x3FD294062ED59F06ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(0x34*I*Pi/128) */	tm2 -= 2;
			tmp64 = 0x3FD111D262B1F677ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(0x35*I*Pi/128) */
			tmp64 = 0x3FCF19F97B215F1Bull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(0x36*I*Pi/128) */
			tmp64 = 0x3FCC0B826A7E4F63ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(0x37*I*Pi/128) */	tmp += 2;
			tmp64 = 0x3FC8F8B83C69A60Bull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(0x38*I*Pi/128) */	tm2 -= 2;
			tmp64 = 0x3FC5E214448B3FC6ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(0x39*I*Pi/128) */
			tmp64 = 0x3FC2C8106E8E613Aull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(0x3A*I*Pi/128) */
			tmp64 = 0x3FBF564E56A9730Eull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(0x3B*I*Pi/128) */	tmp += 2;
			tmp64 = 0x3FB917A6BC29B42Cull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(0x3C*I*Pi/128) */	tm2 -= 2;
			tmp64 = 0x3FB2D52092CE19F6ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(0x3D*I*Pi/128) */
			tmp64 = 0x3FA91F65F10DD814ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(0x3E*I*Pi/128) */
			tmp64 = 0x3F992155F7A3667Eull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(0x3F*I*Pi/128) */	tmp += 2;

			tmp = base_negacyclic_root + RADIX*2;	// reset to point to start of above block
			nbytes = RADIX << (l2_sz_vd-1);	// RADIX*sz_vd/2; 7 AVX-register-sized complex data

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
			// For each of our four 4-bit [16-entry] lookup tables, the index of the .d* field selector
			// indicates the bit of the 4LUT, i.e.Low-order bits at left:
			/* Forward-weight multipliers: 1 for 0-bit, 0.5 for 1-bit: */			// Bitfield, bits ordered 0-3:
			tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;	// [0000]
			tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;	// [1000]
			tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;	// [0100]
			tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;	// [1100]
			tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;	// [0010]
			tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;	// [1010]
			tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;	// [0110]
			tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;	// [1110]
			tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;	// [0001]
			tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;	// [1001]
			tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;	// [0101]
			tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;	// [1101]
			tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;	// [0011]
			tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;	// [1011]
			tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;	// [0111]
			tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;	// [1111]
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

			nbytes = 64 << l2_sz_vd;

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

			nbytes = 16 << l2_sz_vd;

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

		nbytes = 4 << l2_sz_vd;

	  #ifdef USE_AVX
		n_minus_sil   = (struct uint32x4 *)sse_nm1 + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_nm1 + 2;
		sinwt         = (struct uint32x4 *)sse_nm1 + 3;
		sinwtm1       = (struct uint32x4 *)sse_nm1 + 4;
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
		bjmodn00 = (int*)(sse_nm1 + RE_IM_STRIDE);
	  #endif
										bjmodn20 = bjmodn00 + 0x20;
		bjmodn01 = bjmodn00 + 0x01;		bjmodn21 = bjmodn00 + 0x21;
		bjmodn02 = bjmodn00 + 0x02;		bjmodn22 = bjmodn00 + 0x22;
		bjmodn03 = bjmodn00 + 0x03;		bjmodn23 = bjmodn00 + 0x23;
		bjmodn04 = bjmodn00 + 0x04;		bjmodn24 = bjmodn00 + 0x24;
		bjmodn05 = bjmodn00 + 0x05;		bjmodn25 = bjmodn00 + 0x25;
		bjmodn06 = bjmodn00 + 0x06;		bjmodn26 = bjmodn00 + 0x26;
		bjmodn07 = bjmodn00 + 0x07;		bjmodn27 = bjmodn00 + 0x27;
		bjmodn08 = bjmodn00 + 0x08;		bjmodn28 = bjmodn00 + 0x28;
		bjmodn09 = bjmodn00 + 0x09;		bjmodn29 = bjmodn00 + 0x29;
		bjmodn0A = bjmodn00 + 0x0A;		bjmodn2A = bjmodn00 + 0x2A;
		bjmodn0B = bjmodn00 + 0x0B;		bjmodn2B = bjmodn00 + 0x2B;
		bjmodn0C = bjmodn00 + 0x0C;		bjmodn2C = bjmodn00 + 0x2C;
		bjmodn0D = bjmodn00 + 0x0D;		bjmodn2D = bjmodn00 + 0x2D;
		bjmodn0E = bjmodn00 + 0x0E;		bjmodn2E = bjmodn00 + 0x2E;
		bjmodn0F = bjmodn00 + 0x0F;		bjmodn2F = bjmodn00 + 0x2F;
		bjmodn10 = bjmodn00 + 0x10;		bjmodn30 = bjmodn00 + 0x30;
		bjmodn11 = bjmodn00 + 0x11;		bjmodn31 = bjmodn00 + 0x31;
		bjmodn12 = bjmodn00 + 0x12;		bjmodn32 = bjmodn00 + 0x32;
		bjmodn13 = bjmodn00 + 0x13;		bjmodn33 = bjmodn00 + 0x33;
		bjmodn14 = bjmodn00 + 0x14;		bjmodn34 = bjmodn00 + 0x34;
		bjmodn15 = bjmodn00 + 0x15;		bjmodn35 = bjmodn00 + 0x35;
		bjmodn16 = bjmodn00 + 0x16;		bjmodn36 = bjmodn00 + 0x36;
		bjmodn17 = bjmodn00 + 0x17;		bjmodn37 = bjmodn00 + 0x37;
		bjmodn18 = bjmodn00 + 0x18;		bjmodn38 = bjmodn00 + 0x38;
		bjmodn19 = bjmodn00 + 0x19;		bjmodn39 = bjmodn00 + 0x39;
		bjmodn1A = bjmodn00 + 0x1A;		bjmodn3A = bjmodn00 + 0x3A;
		bjmodn1B = bjmodn00 + 0x1B;		bjmodn3B = bjmodn00 + 0x3B;
		bjmodn1C = bjmodn00 + 0x1C;		bjmodn3C = bjmodn00 + 0x3C;
		bjmodn1D = bjmodn00 + 0x1D;		bjmodn3D = bjmodn00 + 0x3D;
		bjmodn1E = bjmodn00 + 0x1E;		bjmodn3E = bjmodn00 + 0x3E;
		bjmodn1F = bjmodn00 + 0x1F;		bjmodn3F = bjmodn00 + 0x3F;

	#endif	// USE_SSE2

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
		#ifdef USE_SSE2
			tdat[ithread].r00 = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (long)tdat[ithread].r00 + ((long)half_arr - (long)r00);
		#else	// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
			tdat[ithread].r00      = (vec_dbl *)base;
			tdat[ithread].half_arr = (vec_dbl *)baseinv;
		#endif	// USE_SSE2
		}
	#endif

		/*   constant index offsets for load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		pini += ( (pini >> DAT_BITS) << PAD_BITS );
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

		ASSERT(HERE, p01+p01 == p02, "p01+p01 != p02");
		ASSERT(HERE, p02+p02 == p04, "p02+p02 != p04");
		ASSERT(HERE, p04+p04 == p08, "p04+p04 != p08");
		ASSERT(HERE, p08+p08 == p10, "p08+p08 != p10");
		ASSERT(HERE, p10+p10 == p20, "p10+p10 != p20");
		ASSERT(HERE, p20+p18 == p38, "p20+p18 != p38");

		if(_cy_r00)	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;

			free((void *)_bjmodn00); _bjmodn00 = 0x0;	free((void *)_bjmodn20); _bjmodn20 = 0x0;
			free((void *)_bjmodn01); _bjmodn01 = 0x0;	free((void *)_bjmodn21); _bjmodn21 = 0x0;
			free((void *)_bjmodn02); _bjmodn02 = 0x0;	free((void *)_bjmodn22); _bjmodn22 = 0x0;
			free((void *)_bjmodn03); _bjmodn03 = 0x0;	free((void *)_bjmodn23); _bjmodn23 = 0x0;
			free((void *)_bjmodn04); _bjmodn04 = 0x0;	free((void *)_bjmodn24); _bjmodn24 = 0x0;
			free((void *)_bjmodn05); _bjmodn05 = 0x0;	free((void *)_bjmodn25); _bjmodn25 = 0x0;
			free((void *)_bjmodn06); _bjmodn06 = 0x0;	free((void *)_bjmodn26); _bjmodn26 = 0x0;
			free((void *)_bjmodn07); _bjmodn07 = 0x0;	free((void *)_bjmodn27); _bjmodn27 = 0x0;
			free((void *)_bjmodn08); _bjmodn08 = 0x0;	free((void *)_bjmodn28); _bjmodn28 = 0x0;
			free((void *)_bjmodn09); _bjmodn09 = 0x0;	free((void *)_bjmodn29); _bjmodn29 = 0x0;
			free((void *)_bjmodn0A); _bjmodn0A = 0x0;	free((void *)_bjmodn2A); _bjmodn2A = 0x0;
			free((void *)_bjmodn0B); _bjmodn0B = 0x0;	free((void *)_bjmodn2B); _bjmodn2B = 0x0;
			free((void *)_bjmodn0C); _bjmodn0C = 0x0;	free((void *)_bjmodn2C); _bjmodn2C = 0x0;
			free((void *)_bjmodn0D); _bjmodn0D = 0x0;	free((void *)_bjmodn2D); _bjmodn2D = 0x0;
			free((void *)_bjmodn0E); _bjmodn0E = 0x0;	free((void *)_bjmodn2E); _bjmodn2E = 0x0;
			free((void *)_bjmodn0F); _bjmodn0F = 0x0;	free((void *)_bjmodn2F); _bjmodn2F = 0x0;
			free((void *)_bjmodn10); _bjmodn10 = 0x0;	free((void *)_bjmodn30); _bjmodn30 = 0x0;
			free((void *)_bjmodn11); _bjmodn11 = 0x0;	free((void *)_bjmodn31); _bjmodn31 = 0x0;
			free((void *)_bjmodn12); _bjmodn12 = 0x0;	free((void *)_bjmodn32); _bjmodn32 = 0x0;
			free((void *)_bjmodn13); _bjmodn13 = 0x0;	free((void *)_bjmodn33); _bjmodn33 = 0x0;
			free((void *)_bjmodn14); _bjmodn14 = 0x0;	free((void *)_bjmodn34); _bjmodn34 = 0x0;
			free((void *)_bjmodn15); _bjmodn15 = 0x0;	free((void *)_bjmodn35); _bjmodn35 = 0x0;
			free((void *)_bjmodn16); _bjmodn16 = 0x0;	free((void *)_bjmodn36); _bjmodn36 = 0x0;
			free((void *)_bjmodn17); _bjmodn17 = 0x0;	free((void *)_bjmodn37); _bjmodn37 = 0x0;
			free((void *)_bjmodn18); _bjmodn18 = 0x0;	free((void *)_bjmodn38); _bjmodn38 = 0x0;
			free((void *)_bjmodn19); _bjmodn19 = 0x0;	free((void *)_bjmodn39); _bjmodn39 = 0x0;
			free((void *)_bjmodn1A); _bjmodn1A = 0x0;	free((void *)_bjmodn3A); _bjmodn3A = 0x0;
			free((void *)_bjmodn1B); _bjmodn1B = 0x0;	free((void *)_bjmodn3B); _bjmodn3B = 0x0;
			free((void *)_bjmodn1C); _bjmodn1C = 0x0;	free((void *)_bjmodn3C); _bjmodn3C = 0x0;
			free((void *)_bjmodn1D); _bjmodn1D = 0x0;	free((void *)_bjmodn3D); _bjmodn3D = 0x0;
			free((void *)_bjmodn1E); _bjmodn1E = 0x0;	free((void *)_bjmodn3E); _bjmodn3E = 0x0;
			free((void *)_bjmodn1F); _bjmodn1F = 0x0;	free((void *)_bjmodn3F); _bjmodn3F = 0x0;

			free((void *)_cy_r00); _cy_r00 = 0x0;	free((void *)_cy_i00); _cy_i00 = 0x0;
			free((void *)_cy_r01); _cy_r01 = 0x0;	free((void *)_cy_i01); _cy_i01 = 0x0;
			free((void *)_cy_r02); _cy_r02 = 0x0;	free((void *)_cy_i02); _cy_i02 = 0x0;
			free((void *)_cy_r03); _cy_r03 = 0x0;	free((void *)_cy_i03); _cy_i03 = 0x0;
			free((void *)_cy_r04); _cy_r04 = 0x0;	free((void *)_cy_i04); _cy_i04 = 0x0;
			free((void *)_cy_r05); _cy_r05 = 0x0;	free((void *)_cy_i05); _cy_i05 = 0x0;
			free((void *)_cy_r06); _cy_r06 = 0x0;	free((void *)_cy_i06); _cy_i06 = 0x0;
			free((void *)_cy_r07); _cy_r07 = 0x0;	free((void *)_cy_i07); _cy_i07 = 0x0;
			free((void *)_cy_r08); _cy_r08 = 0x0;	free((void *)_cy_i08); _cy_i08 = 0x0;
			free((void *)_cy_r09); _cy_r09 = 0x0;	free((void *)_cy_i09); _cy_i09 = 0x0;
			free((void *)_cy_r0A); _cy_r0A = 0x0;	free((void *)_cy_i0A); _cy_i0A = 0x0;
			free((void *)_cy_r0B); _cy_r0B = 0x0;	free((void *)_cy_i0B); _cy_i0B = 0x0;
			free((void *)_cy_r0C); _cy_r0C = 0x0;	free((void *)_cy_i0C); _cy_i0C = 0x0;
			free((void *)_cy_r0D); _cy_r0D = 0x0;	free((void *)_cy_i0D); _cy_i0D = 0x0;
			free((void *)_cy_r0E); _cy_r0E = 0x0;	free((void *)_cy_i0E); _cy_i0E = 0x0;
			free((void *)_cy_r0F); _cy_r0F = 0x0;	free((void *)_cy_i0F); _cy_i0F = 0x0;
			free((void *)_cy_r10); _cy_r10 = 0x0;	free((void *)_cy_i10); _cy_i10 = 0x0;
			free((void *)_cy_r11); _cy_r11 = 0x0;	free((void *)_cy_i11); _cy_i11 = 0x0;
			free((void *)_cy_r12); _cy_r12 = 0x0;	free((void *)_cy_i12); _cy_i12 = 0x0;
			free((void *)_cy_r13); _cy_r13 = 0x0;	free((void *)_cy_i13); _cy_i13 = 0x0;
			free((void *)_cy_r14); _cy_r14 = 0x0;	free((void *)_cy_i14); _cy_i14 = 0x0;
			free((void *)_cy_r15); _cy_r15 = 0x0;	free((void *)_cy_i15); _cy_i15 = 0x0;
			free((void *)_cy_r16); _cy_r16 = 0x0;	free((void *)_cy_i16); _cy_i16 = 0x0;
			free((void *)_cy_r17); _cy_r17 = 0x0;	free((void *)_cy_i17); _cy_i17 = 0x0;
			free((void *)_cy_r18); _cy_r18 = 0x0;	free((void *)_cy_i18); _cy_i18 = 0x0;
			free((void *)_cy_r19); _cy_r19 = 0x0;	free((void *)_cy_i19); _cy_i19 = 0x0;
			free((void *)_cy_r1A); _cy_r1A = 0x0;	free((void *)_cy_i1A); _cy_i1A = 0x0;
			free((void *)_cy_r1B); _cy_r1B = 0x0;	free((void *)_cy_i1B); _cy_i1B = 0x0;
			free((void *)_cy_r1C); _cy_r1C = 0x0;	free((void *)_cy_i1C); _cy_i1C = 0x0;
			free((void *)_cy_r1D); _cy_r1D = 0x0;	free((void *)_cy_i1D); _cy_i1D = 0x0;
			free((void *)_cy_r1E); _cy_r1E = 0x0;	free((void *)_cy_i1E); _cy_i1E = 0x0;
			free((void *)_cy_r1F); _cy_r1F = 0x0;	free((void *)_cy_i1F); _cy_i1F = 0x0;
			free((void *)_cy_r20); _cy_r20 = 0x0;	free((void *)_cy_i20); _cy_i20 = 0x0;
			free((void *)_cy_r21); _cy_r21 = 0x0;	free((void *)_cy_i21); _cy_i21 = 0x0;
			free((void *)_cy_r22); _cy_r22 = 0x0;	free((void *)_cy_i22); _cy_i22 = 0x0;
			free((void *)_cy_r23); _cy_r23 = 0x0;	free((void *)_cy_i23); _cy_i23 = 0x0;
			free((void *)_cy_r24); _cy_r24 = 0x0;	free((void *)_cy_i24); _cy_i24 = 0x0;
			free((void *)_cy_r25); _cy_r25 = 0x0;	free((void *)_cy_i25); _cy_i25 = 0x0;
			free((void *)_cy_r26); _cy_r26 = 0x0;	free((void *)_cy_i26); _cy_i26 = 0x0;
			free((void *)_cy_r27); _cy_r27 = 0x0;	free((void *)_cy_i27); _cy_i27 = 0x0;
			free((void *)_cy_r28); _cy_r28 = 0x0;	free((void *)_cy_i28); _cy_i28 = 0x0;
			free((void *)_cy_r29); _cy_r29 = 0x0;	free((void *)_cy_i29); _cy_i29 = 0x0;
			free((void *)_cy_r2A); _cy_r2A = 0x0;	free((void *)_cy_i2A); _cy_i2A = 0x0;
			free((void *)_cy_r2B); _cy_r2B = 0x0;	free((void *)_cy_i2B); _cy_i2B = 0x0;
			free((void *)_cy_r2C); _cy_r2C = 0x0;	free((void *)_cy_i2C); _cy_i2C = 0x0;
			free((void *)_cy_r2D); _cy_r2D = 0x0;	free((void *)_cy_i2D); _cy_i2D = 0x0;
			free((void *)_cy_r2E); _cy_r2E = 0x0;	free((void *)_cy_i2E); _cy_i2E = 0x0;
			free((void *)_cy_r2F); _cy_r2F = 0x0;	free((void *)_cy_i2F); _cy_i2F = 0x0;
			free((void *)_cy_r30); _cy_r30 = 0x0;	free((void *)_cy_i30); _cy_i30 = 0x0;
			free((void *)_cy_r31); _cy_r31 = 0x0;	free((void *)_cy_i31); _cy_i31 = 0x0;
			free((void *)_cy_r32); _cy_r32 = 0x0;	free((void *)_cy_i32); _cy_i32 = 0x0;
			free((void *)_cy_r33); _cy_r33 = 0x0;	free((void *)_cy_i33); _cy_i33 = 0x0;
			free((void *)_cy_r34); _cy_r34 = 0x0;	free((void *)_cy_i34); _cy_i34 = 0x0;
			free((void *)_cy_r35); _cy_r35 = 0x0;	free((void *)_cy_i35); _cy_i35 = 0x0;
			free((void *)_cy_r36); _cy_r36 = 0x0;	free((void *)_cy_i36); _cy_i36 = 0x0;
			free((void *)_cy_r37); _cy_r37 = 0x0;	free((void *)_cy_i37); _cy_i37 = 0x0;
			free((void *)_cy_r38); _cy_r38 = 0x0;	free((void *)_cy_i38); _cy_i38 = 0x0;
			free((void *)_cy_r39); _cy_r39 = 0x0;	free((void *)_cy_i39); _cy_i39 = 0x0;
			free((void *)_cy_r3A); _cy_r3A = 0x0;	free((void *)_cy_i3A); _cy_i3A = 0x0;
			free((void *)_cy_r3B); _cy_r3B = 0x0;	free((void *)_cy_i3B); _cy_i3B = 0x0;
			free((void *)_cy_r3C); _cy_r3C = 0x0;	free((void *)_cy_i3C); _cy_i3C = 0x0;
			free((void *)_cy_r3D); _cy_r3D = 0x0;	free((void *)_cy_i3D); _cy_i3D = 0x0;
			free((void *)_cy_r3E); _cy_r3E = 0x0;	free((void *)_cy_i3E); _cy_i3E = 0x0;
			free((void *)_cy_r3F); _cy_r3F = 0x0;	free((void *)_cy_i3F); _cy_i3F = 0x0;

			free((void *)_jstart ); _jstart  = 0x0;
			free((void *)_jhi    ); _jhi     = 0x0;
			free((void *)_maxerr); _maxerr = 0x0;
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
		_bjmodn0A	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn0A== 0x0);
		_bjmodn0B	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn0B== 0x0);
		_bjmodn0C	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn0C== 0x0);
		_bjmodn0D	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn0D== 0x0);
		_bjmodn0E	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn0E== 0x0);
		_bjmodn0F	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn0F== 0x0);
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
		_bjmodn1A	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn1A== 0x0);
		_bjmodn1B	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn1B== 0x0);
		_bjmodn1C	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn1C== 0x0);
		_bjmodn1D	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn1D== 0x0);
		_bjmodn1E	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn1E== 0x0);
		_bjmodn1F	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn1F== 0x0);
		_bjmodn20	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn20== 0x0);
		_bjmodn21	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn21== 0x0);
		_bjmodn22	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn22== 0x0);
		_bjmodn23	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn23== 0x0);
		_bjmodn24	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn24== 0x0);
		_bjmodn25	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn25== 0x0);
		_bjmodn26	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn26== 0x0);
		_bjmodn27	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn27== 0x0);
		_bjmodn28	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn28== 0x0);
		_bjmodn29	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn29== 0x0);
		_bjmodn2A	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn2A== 0x0);
		_bjmodn2B	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn2B== 0x0);
		_bjmodn2C	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn2C== 0x0);
		_bjmodn2D	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn2D== 0x0);
		_bjmodn2E	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn2E== 0x0);
		_bjmodn2F	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn2F== 0x0);
		_bjmodn30	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn30== 0x0);
		_bjmodn31	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn31== 0x0);
		_bjmodn32	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn32== 0x0);
		_bjmodn33	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn33== 0x0);
		_bjmodn34	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn34== 0x0);
		_bjmodn35	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn35== 0x0);
		_bjmodn36	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn36== 0x0);
		_bjmodn37	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn37== 0x0);
		_bjmodn38	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn38== 0x0);
		_bjmodn39	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn39== 0x0);
		_bjmodn3A	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn3A== 0x0);
		_bjmodn3B	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn3B== 0x0);
		_bjmodn3C	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn3C== 0x0);
		_bjmodn3D	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn3D== 0x0);
		_bjmodn3E	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn3E== 0x0);
		_bjmodn3F	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn3F== 0x0);
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
		_cy_r0A	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r0A== 0x0);
		_cy_r0B	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r0B== 0x0);
		_cy_r0C	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r0C== 0x0);
		_cy_r0D	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r0D== 0x0);
		_cy_r0E	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r0E== 0x0);
		_cy_r0F	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r0F== 0x0);
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
		_cy_r1A	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r1A== 0x0);
		_cy_r1B	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r1B== 0x0);
		_cy_r1C	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r1C== 0x0);
		_cy_r1D	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r1D== 0x0);
		_cy_r1E	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r1E== 0x0);
		_cy_r1F	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r1F== 0x0);
		_cy_r20	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r20== 0x0);
		_cy_r21	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r21== 0x0);
		_cy_r22	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r22== 0x0);
		_cy_r23	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r23== 0x0);
		_cy_r24	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r24== 0x0);
		_cy_r25	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r25== 0x0);
		_cy_r26	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r26== 0x0);
		_cy_r27	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r27== 0x0);
		_cy_r28	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r28== 0x0);
		_cy_r29	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r29== 0x0);
		_cy_r2A	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r2A== 0x0);
		_cy_r2B	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r2B== 0x0);
		_cy_r2C	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r2C== 0x0);
		_cy_r2D	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r2D== 0x0);
		_cy_r2E	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r2E== 0x0);
		_cy_r2F	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r2F== 0x0);
		_cy_r30	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r30== 0x0);
		_cy_r31	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r31== 0x0);
		_cy_r32	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r32== 0x0);
		_cy_r33	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r33== 0x0);
		_cy_r34	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r34== 0x0);
		_cy_r35	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r35== 0x0);
		_cy_r36	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r36== 0x0);
		_cy_r37	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r37== 0x0);
		_cy_r38	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r38== 0x0);
		_cy_r39	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r39== 0x0);
		_cy_r3A	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r3A== 0x0);
		_cy_r3B	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r3B== 0x0);
		_cy_r3C	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r3C== 0x0);
		_cy_r3D	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r3D== 0x0);
		_cy_r3E	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r3E== 0x0);
		_cy_r3F	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r3F== 0x0);

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
		_cy_i0A	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i0A== 0x0);
		_cy_i0B	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i0B== 0x0);
		_cy_i0C	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i0C== 0x0);
		_cy_i0D	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i0D== 0x0);
		_cy_i0E	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i0E== 0x0);
		_cy_i0F	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i0F== 0x0);
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
		_cy_i1A	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i1A== 0x0);
		_cy_i1B	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i1B== 0x0);
		_cy_i1C	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i1C== 0x0);
		_cy_i1D	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i1D== 0x0);
		_cy_i1E	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i1E== 0x0);
		_cy_i1F	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i1F== 0x0);
		_cy_i20	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i20== 0x0);
		_cy_i21	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i21== 0x0);
		_cy_i22	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i22== 0x0);
		_cy_i23	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i23== 0x0);
		_cy_i24	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i24== 0x0);
		_cy_i25	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i25== 0x0);
		_cy_i26	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i26== 0x0);
		_cy_i27	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i27== 0x0);
		_cy_i28	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i28== 0x0);
		_cy_i29	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i29== 0x0);
		_cy_i2A	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i2A== 0x0);
		_cy_i2B	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i2B== 0x0);
		_cy_i2C	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i2C== 0x0);
		_cy_i2D	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i2D== 0x0);
		_cy_i2E	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i2E== 0x0);
		_cy_i2F	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i2F== 0x0);
		_cy_i30	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i30== 0x0);
		_cy_i31	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i31== 0x0);
		_cy_i32	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i32== 0x0);
		_cy_i33	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i33== 0x0);
		_cy_i34	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i34== 0x0);
		_cy_i35	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i35== 0x0);
		_cy_i36	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i36== 0x0);
		_cy_i37	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i37== 0x0);
		_cy_i38	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i38== 0x0);
		_cy_i39	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i39== 0x0);
		_cy_i3A	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i3A== 0x0);
		_cy_i3B	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i3B== 0x0);
		_cy_i3C	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i3C== 0x0);
		_cy_i3D	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i3D== 0x0);
		_cy_i3E	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i3E== 0x0);
		_cy_i3F	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i3F== 0x0);

		_maxerr	= (double *)malloc(j);	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays!");

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
			i.e. the one that n2/32-separated FFT outputs need:
			*/
			_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
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
			ASSERT(HERE, _bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");
		}
	}	/* endif(first_entry) */

/*...The radix-64 final DIT pass is here.	*/

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
		_cy_r0A[ithread] = 0;	_cy_i0A[ithread] = 0;
		_cy_r0B[ithread] = 0;	_cy_i0B[ithread] = 0;
		_cy_r0C[ithread] = 0;	_cy_i0C[ithread] = 0;
		_cy_r0D[ithread] = 0;	_cy_i0D[ithread] = 0;
		_cy_r0E[ithread] = 0;	_cy_i0E[ithread] = 0;
		_cy_r0F[ithread] = 0;	_cy_i0F[ithread] = 0;
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
		_cy_r1A[ithread] = 0;	_cy_i1A[ithread] = 0;
		_cy_r1B[ithread] = 0;	_cy_i1B[ithread] = 0;
		_cy_r1C[ithread] = 0;	_cy_i1C[ithread] = 0;
		_cy_r1D[ithread] = 0;	_cy_i1D[ithread] = 0;
		_cy_r1E[ithread] = 0;	_cy_i1E[ithread] = 0;
		_cy_r1F[ithread] = 0;	_cy_i1F[ithread] = 0;
		_cy_r20[ithread] = 0;	_cy_i20[ithread] = 0;
		_cy_r21[ithread] = 0;	_cy_i21[ithread] = 0;
		_cy_r22[ithread] = 0;	_cy_i22[ithread] = 0;
		_cy_r23[ithread] = 0;	_cy_i23[ithread] = 0;
		_cy_r24[ithread] = 0;	_cy_i24[ithread] = 0;
		_cy_r25[ithread] = 0;	_cy_i25[ithread] = 0;
		_cy_r26[ithread] = 0;	_cy_i26[ithread] = 0;
		_cy_r27[ithread] = 0;	_cy_i27[ithread] = 0;
		_cy_r28[ithread] = 0;	_cy_i28[ithread] = 0;
		_cy_r29[ithread] = 0;	_cy_i29[ithread] = 0;
		_cy_r2A[ithread] = 0;	_cy_i2A[ithread] = 0;
		_cy_r2B[ithread] = 0;	_cy_i2B[ithread] = 0;
		_cy_r2C[ithread] = 0;	_cy_i2C[ithread] = 0;
		_cy_r2D[ithread] = 0;	_cy_i2D[ithread] = 0;
		_cy_r2E[ithread] = 0;	_cy_i2E[ithread] = 0;
		_cy_r2F[ithread] = 0;	_cy_i2F[ithread] = 0;
		_cy_r30[ithread] = 0;	_cy_i30[ithread] = 0;
		_cy_r31[ithread] = 0;	_cy_i31[ithread] = 0;
		_cy_r32[ithread] = 0;	_cy_i32[ithread] = 0;
		_cy_r33[ithread] = 0;	_cy_i33[ithread] = 0;
		_cy_r34[ithread] = 0;	_cy_i34[ithread] = 0;
		_cy_r35[ithread] = 0;	_cy_i35[ithread] = 0;
		_cy_r36[ithread] = 0;	_cy_i36[ithread] = 0;
		_cy_r37[ithread] = 0;	_cy_i37[ithread] = 0;
		_cy_r38[ithread] = 0;	_cy_i38[ithread] = 0;
		_cy_r39[ithread] = 0;	_cy_i39[ithread] = 0;
		_cy_r3A[ithread] = 0;	_cy_i3A[ithread] = 0;
		_cy_r3B[ithread] = 0;	_cy_i3B[ithread] = 0;
		_cy_r3C[ithread] = 0;	_cy_i3C[ithread] = 0;
		_cy_r3D[ithread] = 0;	_cy_i3D[ithread] = 0;
		_cy_r3E[ithread] = 0;	_cy_i3E[ithread] = 0;
		_cy_r3F[ithread] = 0;	_cy_i3F[ithread] = 0;

		_maxerr[ithread] = 0.0;
	}
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy_r00[0] = -2;
	}

	*fracmax=0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

for(outer=0; outer <= 1; outer++)
{
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		/*
		Moved this inside the outer-loop, so on cleanup pass can use it to reset _col,_co2,_co3 starting values,
		then simply overwrite it with 1 prior to starting the k-loop.
		*/
		_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (_i[0] = 1).	*/
		for(ithread = 1; ithread < CY_THREADS; ithread++)
		{
			_i[ithread] = ((uint32)(sw - _bjmodnini[ithread]) >> 31);
		}

		// Include 0-thread here ... bjmodn terms all 0 for that, but need jhi computed for all threads
		j = _bjmodnini[CY_THREADS];
		khi = n_div_nwt/CY_THREADS;
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
			MOD_ADD32(_bjmodn09[ithread], j, n, _bjmodn0A[ithread]);
			MOD_ADD32(_bjmodn0A[ithread], j, n, _bjmodn0B[ithread]);
			MOD_ADD32(_bjmodn0B[ithread], j, n, _bjmodn0C[ithread]);
			MOD_ADD32(_bjmodn0C[ithread], j, n, _bjmodn0D[ithread]);
			MOD_ADD32(_bjmodn0D[ithread], j, n, _bjmodn0E[ithread]);
			MOD_ADD32(_bjmodn0E[ithread], j, n, _bjmodn0F[ithread]);
			MOD_ADD32(_bjmodn0F[ithread], j, n, _bjmodn10[ithread]);
			MOD_ADD32(_bjmodn10[ithread], j, n, _bjmodn11[ithread]);
			MOD_ADD32(_bjmodn11[ithread], j, n, _bjmodn12[ithread]);
			MOD_ADD32(_bjmodn12[ithread], j, n, _bjmodn13[ithread]);
			MOD_ADD32(_bjmodn13[ithread], j, n, _bjmodn14[ithread]);
			MOD_ADD32(_bjmodn14[ithread], j, n, _bjmodn15[ithread]);
			MOD_ADD32(_bjmodn15[ithread], j, n, _bjmodn16[ithread]);
			MOD_ADD32(_bjmodn16[ithread], j, n, _bjmodn17[ithread]);
			MOD_ADD32(_bjmodn17[ithread], j, n, _bjmodn18[ithread]);
			MOD_ADD32(_bjmodn18[ithread], j, n, _bjmodn19[ithread]);
			MOD_ADD32(_bjmodn19[ithread], j, n, _bjmodn1A[ithread]);
			MOD_ADD32(_bjmodn1A[ithread], j, n, _bjmodn1B[ithread]);
			MOD_ADD32(_bjmodn1B[ithread], j, n, _bjmodn1C[ithread]);
			MOD_ADD32(_bjmodn1C[ithread], j, n, _bjmodn1D[ithread]);
			MOD_ADD32(_bjmodn1D[ithread], j, n, _bjmodn1E[ithread]);
			MOD_ADD32(_bjmodn1E[ithread], j, n, _bjmodn1F[ithread]);
			MOD_ADD32(_bjmodn1F[ithread], j, n, _bjmodn20[ithread]);
			MOD_ADD32(_bjmodn20[ithread], j, n, _bjmodn21[ithread]);
			MOD_ADD32(_bjmodn21[ithread], j, n, _bjmodn22[ithread]);
			MOD_ADD32(_bjmodn22[ithread], j, n, _bjmodn23[ithread]);
			MOD_ADD32(_bjmodn23[ithread], j, n, _bjmodn24[ithread]);
			MOD_ADD32(_bjmodn24[ithread], j, n, _bjmodn25[ithread]);
			MOD_ADD32(_bjmodn25[ithread], j, n, _bjmodn26[ithread]);
			MOD_ADD32(_bjmodn26[ithread], j, n, _bjmodn27[ithread]);
			MOD_ADD32(_bjmodn27[ithread], j, n, _bjmodn28[ithread]);
			MOD_ADD32(_bjmodn28[ithread], j, n, _bjmodn29[ithread]);
			MOD_ADD32(_bjmodn29[ithread], j, n, _bjmodn2A[ithread]);
			MOD_ADD32(_bjmodn2A[ithread], j, n, _bjmodn2B[ithread]);
			MOD_ADD32(_bjmodn2B[ithread], j, n, _bjmodn2C[ithread]);
			MOD_ADD32(_bjmodn2C[ithread], j, n, _bjmodn2D[ithread]);
			MOD_ADD32(_bjmodn2D[ithread], j, n, _bjmodn2E[ithread]);
			MOD_ADD32(_bjmodn2E[ithread], j, n, _bjmodn2F[ithread]);
			MOD_ADD32(_bjmodn2F[ithread], j, n, _bjmodn30[ithread]);
			MOD_ADD32(_bjmodn30[ithread], j, n, _bjmodn31[ithread]);
			MOD_ADD32(_bjmodn31[ithread], j, n, _bjmodn32[ithread]);
			MOD_ADD32(_bjmodn32[ithread], j, n, _bjmodn33[ithread]);
			MOD_ADD32(_bjmodn33[ithread], j, n, _bjmodn34[ithread]);
			MOD_ADD32(_bjmodn34[ithread], j, n, _bjmodn35[ithread]);
			MOD_ADD32(_bjmodn35[ithread], j, n, _bjmodn36[ithread]);
			MOD_ADD32(_bjmodn36[ithread], j, n, _bjmodn37[ithread]);
			MOD_ADD32(_bjmodn37[ithread], j, n, _bjmodn38[ithread]);
			MOD_ADD32(_bjmodn38[ithread], j, n, _bjmodn39[ithread]);
			MOD_ADD32(_bjmodn39[ithread], j, n, _bjmodn3A[ithread]);
			MOD_ADD32(_bjmodn3A[ithread], j, n, _bjmodn3B[ithread]);
			MOD_ADD32(_bjmodn3B[ithread], j, n, _bjmodn3C[ithread]);
			MOD_ADD32(_bjmodn3C[ithread], j, n, _bjmodn3D[ithread]);
			MOD_ADD32(_bjmodn3D[ithread], j, n, _bjmodn3E[ithread]);
			MOD_ADD32(_bjmodn3E[ithread], j, n, _bjmodn3F[ithread]);

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
	}

#if defined(USE_SSE2) && defined(USE_PTHREAD)

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
		tdat[ithread].maxerr = _maxerr[ithread];
		tdat[ithread].scale = scale;

	// pointer data:
		ASSERT(HERE, tdat[ithread].arrdat == a, "thread-local memcheck fail!");			/* Main data array */
		ASSERT(HERE, tdat[ithread].wt0 == wt0, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].wt1 == wt1, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].si  == si, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].rn0 == rn0, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].rn1 == rn1, "thread-local memcheck fail!");
	#ifdef USE_SSE2
		ASSERT(HERE, tdat[ithread].r00 == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].r00;		dtmp = -ISRT2;
		ASSERT(HERE, ((tmp + 0x100)->d0 == -ISRT2 && (tmp + 0x100)->d1 == -ISRT2), "thread-local memcheck failed!");
		ASSERT(HERE, ((tmp + 0x101)->d0 == +ISRT2 && (tmp + 0x101)->d1 == +ISRT2), "thread-local memcheck failed!");
		tmp = tdat[ithread].half_arr;
		ASSERT(HERE, ((tmp-1)->d0 == crnd && (tmp-1)->d1 == crnd), "thread-local memcheck failed!");
	#endif
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX
			// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
			dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#elif defined(USE_SSE2)
			dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#endif
			tdat[ithread].bjmodn00 = _bjmodn00[ithread];	tdat[ithread].bjmodn20 = _bjmodn20[ithread];
			tdat[ithread].bjmodn01 = _bjmodn01[ithread];	tdat[ithread].bjmodn21 = _bjmodn21[ithread];
			tdat[ithread].bjmodn02 = _bjmodn02[ithread];	tdat[ithread].bjmodn22 = _bjmodn22[ithread];
			tdat[ithread].bjmodn03 = _bjmodn03[ithread];	tdat[ithread].bjmodn23 = _bjmodn23[ithread];
			tdat[ithread].bjmodn04 = _bjmodn04[ithread];	tdat[ithread].bjmodn24 = _bjmodn24[ithread];
			tdat[ithread].bjmodn05 = _bjmodn05[ithread];	tdat[ithread].bjmodn25 = _bjmodn25[ithread];
			tdat[ithread].bjmodn06 = _bjmodn06[ithread];	tdat[ithread].bjmodn26 = _bjmodn26[ithread];
			tdat[ithread].bjmodn07 = _bjmodn07[ithread];	tdat[ithread].bjmodn27 = _bjmodn27[ithread];
			tdat[ithread].bjmodn08 = _bjmodn08[ithread];	tdat[ithread].bjmodn28 = _bjmodn28[ithread];
			tdat[ithread].bjmodn09 = _bjmodn09[ithread];	tdat[ithread].bjmodn29 = _bjmodn29[ithread];
			tdat[ithread].bjmodn0A = _bjmodn0A[ithread];	tdat[ithread].bjmodn2A = _bjmodn2A[ithread];
			tdat[ithread].bjmodn0B = _bjmodn0B[ithread];	tdat[ithread].bjmodn2B = _bjmodn2B[ithread];
			tdat[ithread].bjmodn0C = _bjmodn0C[ithread];	tdat[ithread].bjmodn2C = _bjmodn2C[ithread];
			tdat[ithread].bjmodn0D = _bjmodn0D[ithread];	tdat[ithread].bjmodn2D = _bjmodn2D[ithread];
			tdat[ithread].bjmodn0E = _bjmodn0E[ithread];	tdat[ithread].bjmodn2E = _bjmodn2E[ithread];
			tdat[ithread].bjmodn0F = _bjmodn0F[ithread];	tdat[ithread].bjmodn2F = _bjmodn2F[ithread];
			tdat[ithread].bjmodn10 = _bjmodn10[ithread];	tdat[ithread].bjmodn30 = _bjmodn30[ithread];
			tdat[ithread].bjmodn11 = _bjmodn11[ithread];	tdat[ithread].bjmodn31 = _bjmodn31[ithread];
			tdat[ithread].bjmodn12 = _bjmodn12[ithread];	tdat[ithread].bjmodn32 = _bjmodn32[ithread];
			tdat[ithread].bjmodn13 = _bjmodn13[ithread];	tdat[ithread].bjmodn33 = _bjmodn33[ithread];
			tdat[ithread].bjmodn14 = _bjmodn14[ithread];	tdat[ithread].bjmodn34 = _bjmodn34[ithread];
			tdat[ithread].bjmodn15 = _bjmodn15[ithread];	tdat[ithread].bjmodn35 = _bjmodn35[ithread];
			tdat[ithread].bjmodn16 = _bjmodn16[ithread];	tdat[ithread].bjmodn36 = _bjmodn36[ithread];
			tdat[ithread].bjmodn17 = _bjmodn17[ithread];	tdat[ithread].bjmodn37 = _bjmodn37[ithread];
			tdat[ithread].bjmodn18 = _bjmodn18[ithread];	tdat[ithread].bjmodn38 = _bjmodn38[ithread];
			tdat[ithread].bjmodn19 = _bjmodn19[ithread];	tdat[ithread].bjmodn39 = _bjmodn39[ithread];
			tdat[ithread].bjmodn1A = _bjmodn1A[ithread];	tdat[ithread].bjmodn3A = _bjmodn3A[ithread];
			tdat[ithread].bjmodn1B = _bjmodn1B[ithread];	tdat[ithread].bjmodn3B = _bjmodn3B[ithread];
			tdat[ithread].bjmodn1C = _bjmodn1C[ithread];	tdat[ithread].bjmodn3C = _bjmodn3C[ithread];
			tdat[ithread].bjmodn1D = _bjmodn1D[ithread];	tdat[ithread].bjmodn3D = _bjmodn3D[ithread];
			tdat[ithread].bjmodn1E = _bjmodn1E[ithread];	tdat[ithread].bjmodn3E = _bjmodn3E[ithread];
			tdat[ithread].bjmodn1F = _bjmodn1F[ithread];	tdat[ithread].bjmodn3F = _bjmodn3F[ithread];
			/* init carries	*/
			tdat[ithread].cy_r00 = _cy_r00[ithread];		tdat[ithread].cy_r20 = _cy_r20[ithread];
			tdat[ithread].cy_r01 = _cy_r01[ithread];		tdat[ithread].cy_r21 = _cy_r21[ithread];
			tdat[ithread].cy_r02 = _cy_r02[ithread];		tdat[ithread].cy_r22 = _cy_r22[ithread];
			tdat[ithread].cy_r03 = _cy_r03[ithread];		tdat[ithread].cy_r23 = _cy_r23[ithread];
			tdat[ithread].cy_r04 = _cy_r04[ithread];		tdat[ithread].cy_r24 = _cy_r24[ithread];
			tdat[ithread].cy_r05 = _cy_r05[ithread];		tdat[ithread].cy_r25 = _cy_r25[ithread];
			tdat[ithread].cy_r06 = _cy_r06[ithread];		tdat[ithread].cy_r26 = _cy_r26[ithread];
			tdat[ithread].cy_r07 = _cy_r07[ithread];		tdat[ithread].cy_r27 = _cy_r27[ithread];
			tdat[ithread].cy_r08 = _cy_r08[ithread];		tdat[ithread].cy_r28 = _cy_r28[ithread];
			tdat[ithread].cy_r09 = _cy_r09[ithread];		tdat[ithread].cy_r29 = _cy_r29[ithread];
			tdat[ithread].cy_r0A = _cy_r0A[ithread];		tdat[ithread].cy_r2A = _cy_r2A[ithread];
			tdat[ithread].cy_r0B = _cy_r0B[ithread];		tdat[ithread].cy_r2B = _cy_r2B[ithread];
			tdat[ithread].cy_r0C = _cy_r0C[ithread];		tdat[ithread].cy_r2C = _cy_r2C[ithread];
			tdat[ithread].cy_r0D = _cy_r0D[ithread];		tdat[ithread].cy_r2D = _cy_r2D[ithread];
			tdat[ithread].cy_r0E = _cy_r0E[ithread];		tdat[ithread].cy_r2E = _cy_r2E[ithread];
			tdat[ithread].cy_r0F = _cy_r0F[ithread];		tdat[ithread].cy_r2F = _cy_r2F[ithread];
			tdat[ithread].cy_r10 = _cy_r10[ithread];		tdat[ithread].cy_r30 = _cy_r30[ithread];
			tdat[ithread].cy_r11 = _cy_r11[ithread];		tdat[ithread].cy_r31 = _cy_r31[ithread];
			tdat[ithread].cy_r12 = _cy_r12[ithread];		tdat[ithread].cy_r32 = _cy_r32[ithread];
			tdat[ithread].cy_r13 = _cy_r13[ithread];		tdat[ithread].cy_r33 = _cy_r33[ithread];
			tdat[ithread].cy_r14 = _cy_r14[ithread];		tdat[ithread].cy_r34 = _cy_r34[ithread];
			tdat[ithread].cy_r15 = _cy_r15[ithread];		tdat[ithread].cy_r35 = _cy_r35[ithread];
			tdat[ithread].cy_r16 = _cy_r16[ithread];		tdat[ithread].cy_r36 = _cy_r36[ithread];
			tdat[ithread].cy_r17 = _cy_r17[ithread];		tdat[ithread].cy_r37 = _cy_r37[ithread];
			tdat[ithread].cy_r18 = _cy_r18[ithread];		tdat[ithread].cy_r38 = _cy_r38[ithread];
			tdat[ithread].cy_r19 = _cy_r19[ithread];		tdat[ithread].cy_r39 = _cy_r39[ithread];
			tdat[ithread].cy_r1A = _cy_r1A[ithread];		tdat[ithread].cy_r3A = _cy_r3A[ithread];
			tdat[ithread].cy_r1B = _cy_r1B[ithread];		tdat[ithread].cy_r3B = _cy_r3B[ithread];
			tdat[ithread].cy_r1C = _cy_r1C[ithread];		tdat[ithread].cy_r3C = _cy_r3C[ithread];
			tdat[ithread].cy_r1D = _cy_r1D[ithread];		tdat[ithread].cy_r3D = _cy_r3D[ithread];
			tdat[ithread].cy_r1E = _cy_r1E[ithread];		tdat[ithread].cy_r3E = _cy_r3E[ithread];
			tdat[ithread].cy_r1F = _cy_r1F[ithread];		tdat[ithread].cy_r3F = _cy_r3F[ithread];
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
		#ifdef USE_SSE2
			// This is slightly different for power-of-2 DFTs: Here, scale is in the +2 slot, base & baseinv remain fixed in 0,+1 slots:
			dtmp = tmp->d0 * (tmp+1)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = tmp->d1 * (tmp+1)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#endif
			// scale gets set immediately prior to calling carry macro, hence no use checking it here.
			/* init carries	*/
			tdat[ithread].cy_r00 = _cy_r00[ithread];		tdat[ithread].cy_r20 = _cy_r20[ithread];
			tdat[ithread].cy_r01 = _cy_r01[ithread];		tdat[ithread].cy_r21 = _cy_r21[ithread];
			tdat[ithread].cy_r02 = _cy_r02[ithread];		tdat[ithread].cy_r22 = _cy_r22[ithread];
			tdat[ithread].cy_r03 = _cy_r03[ithread];		tdat[ithread].cy_r23 = _cy_r23[ithread];
			tdat[ithread].cy_r04 = _cy_r04[ithread];		tdat[ithread].cy_r24 = _cy_r24[ithread];
			tdat[ithread].cy_r05 = _cy_r05[ithread];		tdat[ithread].cy_r25 = _cy_r25[ithread];
			tdat[ithread].cy_r06 = _cy_r06[ithread];		tdat[ithread].cy_r26 = _cy_r26[ithread];
			tdat[ithread].cy_r07 = _cy_r07[ithread];		tdat[ithread].cy_r27 = _cy_r27[ithread];
			tdat[ithread].cy_r08 = _cy_r08[ithread];		tdat[ithread].cy_r28 = _cy_r28[ithread];
			tdat[ithread].cy_r09 = _cy_r09[ithread];		tdat[ithread].cy_r29 = _cy_r29[ithread];
			tdat[ithread].cy_r0A = _cy_r0A[ithread];		tdat[ithread].cy_r2A = _cy_r2A[ithread];
			tdat[ithread].cy_r0B = _cy_r0B[ithread];		tdat[ithread].cy_r2B = _cy_r2B[ithread];
			tdat[ithread].cy_r0C = _cy_r0C[ithread];		tdat[ithread].cy_r2C = _cy_r2C[ithread];
			tdat[ithread].cy_r0D = _cy_r0D[ithread];		tdat[ithread].cy_r2D = _cy_r2D[ithread];
			tdat[ithread].cy_r0E = _cy_r0E[ithread];		tdat[ithread].cy_r2E = _cy_r2E[ithread];
			tdat[ithread].cy_r0F = _cy_r0F[ithread];		tdat[ithread].cy_r2F = _cy_r2F[ithread];
			tdat[ithread].cy_r10 = _cy_r10[ithread];		tdat[ithread].cy_r30 = _cy_r30[ithread];
			tdat[ithread].cy_r11 = _cy_r11[ithread];		tdat[ithread].cy_r31 = _cy_r31[ithread];
			tdat[ithread].cy_r12 = _cy_r12[ithread];		tdat[ithread].cy_r32 = _cy_r32[ithread];
			tdat[ithread].cy_r13 = _cy_r13[ithread];		tdat[ithread].cy_r33 = _cy_r33[ithread];
			tdat[ithread].cy_r14 = _cy_r14[ithread];		tdat[ithread].cy_r34 = _cy_r34[ithread];
			tdat[ithread].cy_r15 = _cy_r15[ithread];		tdat[ithread].cy_r35 = _cy_r35[ithread];
			tdat[ithread].cy_r16 = _cy_r16[ithread];		tdat[ithread].cy_r36 = _cy_r36[ithread];
			tdat[ithread].cy_r17 = _cy_r17[ithread];		tdat[ithread].cy_r37 = _cy_r37[ithread];
			tdat[ithread].cy_r18 = _cy_r18[ithread];		tdat[ithread].cy_r38 = _cy_r38[ithread];
			tdat[ithread].cy_r19 = _cy_r19[ithread];		tdat[ithread].cy_r39 = _cy_r39[ithread];
			tdat[ithread].cy_r1A = _cy_r1A[ithread];		tdat[ithread].cy_r3A = _cy_r3A[ithread];
			tdat[ithread].cy_r1B = _cy_r1B[ithread];		tdat[ithread].cy_r3B = _cy_r3B[ithread];
			tdat[ithread].cy_r1C = _cy_r1C[ithread];		tdat[ithread].cy_r3C = _cy_r3C[ithread];
			tdat[ithread].cy_r1D = _cy_r1D[ithread];		tdat[ithread].cy_r3D = _cy_r3D[ithread];
			tdat[ithread].cy_r1E = _cy_r1E[ithread];		tdat[ithread].cy_r3E = _cy_r3E[ithread];
			tdat[ithread].cy_r1F = _cy_r1F[ithread];		tdat[ithread].cy_r3F = _cy_r3F[ithread];

			tdat[ithread].cy_i00 = _cy_i00[ithread];		tdat[ithread].cy_i20 = _cy_i20[ithread];
			tdat[ithread].cy_i01 = _cy_i01[ithread];		tdat[ithread].cy_i21 = _cy_i21[ithread];
			tdat[ithread].cy_i02 = _cy_i02[ithread];		tdat[ithread].cy_i22 = _cy_i22[ithread];
			tdat[ithread].cy_i03 = _cy_i03[ithread];		tdat[ithread].cy_i23 = _cy_i23[ithread];
			tdat[ithread].cy_i04 = _cy_i04[ithread];		tdat[ithread].cy_i24 = _cy_i24[ithread];
			tdat[ithread].cy_i05 = _cy_i05[ithread];		tdat[ithread].cy_i25 = _cy_i25[ithread];
			tdat[ithread].cy_i06 = _cy_i06[ithread];		tdat[ithread].cy_i26 = _cy_i26[ithread];
			tdat[ithread].cy_i07 = _cy_i07[ithread];		tdat[ithread].cy_i27 = _cy_i27[ithread];
			tdat[ithread].cy_i08 = _cy_i08[ithread];		tdat[ithread].cy_i28 = _cy_i28[ithread];
			tdat[ithread].cy_i09 = _cy_i09[ithread];		tdat[ithread].cy_i29 = _cy_i29[ithread];
			tdat[ithread].cy_i0A = _cy_i0A[ithread];		tdat[ithread].cy_i2A = _cy_i2A[ithread];
			tdat[ithread].cy_i0B = _cy_i0B[ithread];		tdat[ithread].cy_i2B = _cy_i2B[ithread];
			tdat[ithread].cy_i0C = _cy_i0C[ithread];		tdat[ithread].cy_i2C = _cy_i2C[ithread];
			tdat[ithread].cy_i0D = _cy_i0D[ithread];		tdat[ithread].cy_i2D = _cy_i2D[ithread];
			tdat[ithread].cy_i0E = _cy_i0E[ithread];		tdat[ithread].cy_i2E = _cy_i2E[ithread];
			tdat[ithread].cy_i0F = _cy_i0F[ithread];		tdat[ithread].cy_i2F = _cy_i2F[ithread];
			tdat[ithread].cy_i10 = _cy_i10[ithread];		tdat[ithread].cy_i30 = _cy_i30[ithread];
			tdat[ithread].cy_i11 = _cy_i11[ithread];		tdat[ithread].cy_i31 = _cy_i31[ithread];
			tdat[ithread].cy_i12 = _cy_i12[ithread];		tdat[ithread].cy_i32 = _cy_i32[ithread];
			tdat[ithread].cy_i13 = _cy_i13[ithread];		tdat[ithread].cy_i33 = _cy_i33[ithread];
			tdat[ithread].cy_i14 = _cy_i14[ithread];		tdat[ithread].cy_i34 = _cy_i34[ithread];
			tdat[ithread].cy_i15 = _cy_i15[ithread];		tdat[ithread].cy_i35 = _cy_i35[ithread];
			tdat[ithread].cy_i16 = _cy_i16[ithread];		tdat[ithread].cy_i36 = _cy_i36[ithread];
			tdat[ithread].cy_i17 = _cy_i17[ithread];		tdat[ithread].cy_i37 = _cy_i37[ithread];
			tdat[ithread].cy_i18 = _cy_i18[ithread];		tdat[ithread].cy_i38 = _cy_i38[ithread];
			tdat[ithread].cy_i19 = _cy_i19[ithread];		tdat[ithread].cy_i39 = _cy_i39[ithread];
			tdat[ithread].cy_i1A = _cy_i1A[ithread];		tdat[ithread].cy_i3A = _cy_i3A[ithread];
			tdat[ithread].cy_i1B = _cy_i1B[ithread];		tdat[ithread].cy_i3B = _cy_i3B[ithread];
			tdat[ithread].cy_i1C = _cy_i1C[ithread];		tdat[ithread].cy_i3C = _cy_i3C[ithread];
			tdat[ithread].cy_i1D = _cy_i1D[ithread];		tdat[ithread].cy_i3D = _cy_i3D[ithread];
			tdat[ithread].cy_i1E = _cy_i1E[ithread];		tdat[ithread].cy_i3E = _cy_i3E[ithread];
			tdat[ithread].cy_i1F = _cy_i1F[ithread];		tdat[ithread].cy_i3F = _cy_i3F[ithread];
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
		/***** DEC/HP CC doesn't properly copy init value of maxerr = 0 into threads,
		so need to set once again explicitly for each: *****/
		maxerr = 0.0;
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
		#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
			*bjmodn00 = _bjmodn00[ithread];		*bjmodn20 = _bjmodn20[ithread];
			*bjmodn01 = _bjmodn01[ithread];		*bjmodn21 = _bjmodn21[ithread];
			*bjmodn02 = _bjmodn02[ithread];		*bjmodn22 = _bjmodn22[ithread];
			*bjmodn03 = _bjmodn03[ithread];		*bjmodn23 = _bjmodn23[ithread];
			*bjmodn04 = _bjmodn04[ithread];		*bjmodn24 = _bjmodn24[ithread];
			*bjmodn05 = _bjmodn05[ithread];		*bjmodn25 = _bjmodn25[ithread];
			*bjmodn06 = _bjmodn06[ithread];		*bjmodn26 = _bjmodn26[ithread];
			*bjmodn07 = _bjmodn07[ithread];		*bjmodn27 = _bjmodn27[ithread];
			*bjmodn08 = _bjmodn08[ithread];		*bjmodn28 = _bjmodn28[ithread];
			*bjmodn09 = _bjmodn09[ithread];		*bjmodn29 = _bjmodn29[ithread];
			*bjmodn0A = _bjmodn0A[ithread];		*bjmodn2A = _bjmodn2A[ithread];
			*bjmodn0B = _bjmodn0B[ithread];		*bjmodn2B = _bjmodn2B[ithread];
			*bjmodn0C = _bjmodn0C[ithread];		*bjmodn2C = _bjmodn2C[ithread];
			*bjmodn0D = _bjmodn0D[ithread];		*bjmodn2D = _bjmodn2D[ithread];
			*bjmodn0E = _bjmodn0E[ithread];		*bjmodn2E = _bjmodn2E[ithread];
			*bjmodn0F = _bjmodn0F[ithread];		*bjmodn2F = _bjmodn2F[ithread];
			*bjmodn10 = _bjmodn10[ithread];		*bjmodn30 = _bjmodn30[ithread];
			*bjmodn11 = _bjmodn11[ithread];		*bjmodn31 = _bjmodn31[ithread];
			*bjmodn12 = _bjmodn12[ithread];		*bjmodn32 = _bjmodn32[ithread];
			*bjmodn13 = _bjmodn13[ithread];		*bjmodn33 = _bjmodn33[ithread];
			*bjmodn14 = _bjmodn14[ithread];		*bjmodn34 = _bjmodn34[ithread];
			*bjmodn15 = _bjmodn15[ithread];		*bjmodn35 = _bjmodn35[ithread];
			*bjmodn16 = _bjmodn16[ithread];		*bjmodn36 = _bjmodn36[ithread];
			*bjmodn17 = _bjmodn17[ithread];		*bjmodn37 = _bjmodn37[ithread];
			*bjmodn18 = _bjmodn18[ithread];		*bjmodn38 = _bjmodn38[ithread];
			*bjmodn19 = _bjmodn19[ithread];		*bjmodn39 = _bjmodn39[ithread];
			*bjmodn1A = _bjmodn1A[ithread];		*bjmodn3A = _bjmodn3A[ithread];
			*bjmodn1B = _bjmodn1B[ithread];		*bjmodn3B = _bjmodn3B[ithread];
			*bjmodn1C = _bjmodn1C[ithread];		*bjmodn3C = _bjmodn3C[ithread];
			*bjmodn1D = _bjmodn1D[ithread];		*bjmodn3D = _bjmodn3D[ithread];
			*bjmodn1E = _bjmodn1E[ithread];		*bjmodn3E = _bjmodn3E[ithread];
			*bjmodn1F = _bjmodn1F[ithread];		*bjmodn3F = _bjmodn3F[ithread];
		#else
			bjmodn00 = _bjmodn00[ithread];		bjmodn20 = _bjmodn20[ithread];
			bjmodn01 = _bjmodn01[ithread];		bjmodn21 = _bjmodn21[ithread];
			bjmodn02 = _bjmodn02[ithread];		bjmodn22 = _bjmodn22[ithread];
			bjmodn03 = _bjmodn03[ithread];		bjmodn23 = _bjmodn23[ithread];
			bjmodn04 = _bjmodn04[ithread];		bjmodn24 = _bjmodn24[ithread];
			bjmodn05 = _bjmodn05[ithread];		bjmodn25 = _bjmodn25[ithread];
			bjmodn06 = _bjmodn06[ithread];		bjmodn26 = _bjmodn26[ithread];
			bjmodn07 = _bjmodn07[ithread];		bjmodn27 = _bjmodn27[ithread];
			bjmodn08 = _bjmodn08[ithread];		bjmodn28 = _bjmodn28[ithread];
			bjmodn09 = _bjmodn09[ithread];		bjmodn29 = _bjmodn29[ithread];
			bjmodn0A = _bjmodn0A[ithread];		bjmodn2A = _bjmodn2A[ithread];
			bjmodn0B = _bjmodn0B[ithread];		bjmodn2B = _bjmodn2B[ithread];
			bjmodn0C = _bjmodn0C[ithread];		bjmodn2C = _bjmodn2C[ithread];
			bjmodn0D = _bjmodn0D[ithread];		bjmodn2D = _bjmodn2D[ithread];
			bjmodn0E = _bjmodn0E[ithread];		bjmodn2E = _bjmodn2E[ithread];
			bjmodn0F = _bjmodn0F[ithread];		bjmodn2F = _bjmodn2F[ithread];
			bjmodn10 = _bjmodn10[ithread];		bjmodn30 = _bjmodn30[ithread];
			bjmodn11 = _bjmodn11[ithread];		bjmodn31 = _bjmodn31[ithread];
			bjmodn12 = _bjmodn12[ithread];		bjmodn32 = _bjmodn32[ithread];
			bjmodn13 = _bjmodn13[ithread];		bjmodn33 = _bjmodn33[ithread];
			bjmodn14 = _bjmodn14[ithread];		bjmodn34 = _bjmodn34[ithread];
			bjmodn15 = _bjmodn15[ithread];		bjmodn35 = _bjmodn35[ithread];
			bjmodn16 = _bjmodn16[ithread];		bjmodn36 = _bjmodn36[ithread];
			bjmodn17 = _bjmodn17[ithread];		bjmodn37 = _bjmodn37[ithread];
			bjmodn18 = _bjmodn18[ithread];		bjmodn38 = _bjmodn38[ithread];
			bjmodn19 = _bjmodn19[ithread];		bjmodn39 = _bjmodn39[ithread];
			bjmodn1A = _bjmodn1A[ithread];		bjmodn3A = _bjmodn3A[ithread];
			bjmodn1B = _bjmodn1B[ithread];		bjmodn3B = _bjmodn3B[ithread];
			bjmodn1C = _bjmodn1C[ithread];		bjmodn3C = _bjmodn3C[ithread];
			bjmodn1D = _bjmodn1D[ithread];		bjmodn3D = _bjmodn3D[ithread];
			bjmodn1E = _bjmodn1E[ithread];		bjmodn3E = _bjmodn3E[ithread];
			bjmodn1F = _bjmodn1F[ithread];		bjmodn3F = _bjmodn3F[ithread];
		#endif
			/* init carries	*/
		#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
		  #ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
			cy_r00->d0 = _cy_r00[ithread];	cy_r00->d1 = _cy_r01[ithread];	cy_r00->d2 = _cy_r02[ithread];	cy_r00->d3 = _cy_r03[ithread];
			cy_r04->d0 = _cy_r04[ithread];	cy_r04->d1 = _cy_r05[ithread];	cy_r04->d2 = _cy_r06[ithread];	cy_r04->d3 = _cy_r07[ithread];
			cy_r08->d0 = _cy_r08[ithread];	cy_r08->d1 = _cy_r09[ithread];	cy_r08->d2 = _cy_r0A[ithread];	cy_r08->d3 = _cy_r0B[ithread];
			cy_r0C->d0 = _cy_r0C[ithread];	cy_r0C->d1 = _cy_r0D[ithread];	cy_r0C->d2 = _cy_r0E[ithread];	cy_r0C->d3 = _cy_r0F[ithread];
			cy_r10->d0 = _cy_r10[ithread];	cy_r10->d1 = _cy_r11[ithread];	cy_r10->d2 = _cy_r12[ithread];	cy_r10->d3 = _cy_r13[ithread];
			cy_r14->d0 = _cy_r14[ithread];	cy_r14->d1 = _cy_r15[ithread];	cy_r14->d2 = _cy_r16[ithread];	cy_r14->d3 = _cy_r17[ithread];
			cy_r18->d0 = _cy_r18[ithread];	cy_r18->d1 = _cy_r19[ithread];	cy_r18->d2 = _cy_r1A[ithread];	cy_r18->d3 = _cy_r1B[ithread];
			cy_r1C->d0 = _cy_r1C[ithread];	cy_r1C->d1 = _cy_r1D[ithread];	cy_r1C->d2 = _cy_r1E[ithread];	cy_r1C->d3 = _cy_r1F[ithread];
			cy_r20->d0 = _cy_r20[ithread];	cy_r20->d1 = _cy_r21[ithread];	cy_r20->d2 = _cy_r22[ithread];	cy_r20->d3 = _cy_r23[ithread];
			cy_r24->d0 = _cy_r24[ithread];	cy_r24->d1 = _cy_r25[ithread];	cy_r24->d2 = _cy_r26[ithread];	cy_r24->d3 = _cy_r27[ithread];
			cy_r28->d0 = _cy_r28[ithread];	cy_r28->d1 = _cy_r29[ithread];	cy_r28->d2 = _cy_r2A[ithread];	cy_r28->d3 = _cy_r2B[ithread];
			cy_r2C->d0 = _cy_r2C[ithread];	cy_r2C->d1 = _cy_r2D[ithread];	cy_r2C->d2 = _cy_r2E[ithread];	cy_r2C->d3 = _cy_r2F[ithread];
			cy_r30->d0 = _cy_r30[ithread];	cy_r30->d1 = _cy_r31[ithread];	cy_r30->d2 = _cy_r32[ithread];	cy_r30->d3 = _cy_r33[ithread];
			cy_r34->d0 = _cy_r34[ithread];	cy_r34->d1 = _cy_r35[ithread];	cy_r34->d2 = _cy_r36[ithread];	cy_r34->d3 = _cy_r37[ithread];
			cy_r38->d0 = _cy_r38[ithread];	cy_r38->d1 = _cy_r39[ithread];	cy_r38->d2 = _cy_r3A[ithread];	cy_r38->d3 = _cy_r3B[ithread];
			cy_r3C->d0 = _cy_r3C[ithread];	cy_r3C->d1 = _cy_r3D[ithread];	cy_r3C->d2 = _cy_r3E[ithread];	cy_r3C->d3 = _cy_r3F[ithread];
		  #else	// USE_SSE2
			cy_r00->d0 = _cy_r00[ithread];	cy_r00->d1 = _cy_r01[ithread];
			cy_r02->d0 = _cy_r02[ithread];	cy_r02->d1 = _cy_r03[ithread];
			cy_r04->d0 = _cy_r04[ithread];	cy_r04->d1 = _cy_r05[ithread];
			cy_r06->d0 = _cy_r06[ithread];	cy_r06->d1 = _cy_r07[ithread];
			cy_r08->d0 = _cy_r08[ithread];	cy_r08->d1 = _cy_r09[ithread];
			cy_r0A->d0 = _cy_r0A[ithread];	cy_r0A->d1 = _cy_r0B[ithread];
			cy_r0C->d0 = _cy_r0C[ithread];	cy_r0C->d1 = _cy_r0D[ithread];
			cy_r0E->d0 = _cy_r0E[ithread];	cy_r0E->d1 = _cy_r0F[ithread];
			cy_r10->d0 = _cy_r10[ithread];	cy_r10->d1 = _cy_r11[ithread];
			cy_r12->d0 = _cy_r12[ithread];	cy_r12->d1 = _cy_r13[ithread];
			cy_r14->d0 = _cy_r14[ithread];	cy_r14->d1 = _cy_r15[ithread];
			cy_r16->d0 = _cy_r16[ithread];	cy_r16->d1 = _cy_r17[ithread];
			cy_r18->d0 = _cy_r18[ithread];	cy_r18->d1 = _cy_r19[ithread];
			cy_r1A->d0 = _cy_r1A[ithread];	cy_r1A->d1 = _cy_r1B[ithread];
			cy_r1C->d0 = _cy_r1C[ithread];	cy_r1C->d1 = _cy_r1D[ithread];
			cy_r1E->d0 = _cy_r1E[ithread];	cy_r1E->d1 = _cy_r1F[ithread];
			cy_r20->d0 = _cy_r20[ithread];	cy_r20->d1 = _cy_r21[ithread];
			cy_r22->d0 = _cy_r22[ithread];	cy_r22->d1 = _cy_r23[ithread];
			cy_r24->d0 = _cy_r24[ithread];	cy_r24->d1 = _cy_r25[ithread];
			cy_r26->d0 = _cy_r26[ithread];	cy_r26->d1 = _cy_r27[ithread];
			cy_r28->d0 = _cy_r28[ithread];	cy_r28->d1 = _cy_r29[ithread];
			cy_r2A->d0 = _cy_r2A[ithread];	cy_r2A->d1 = _cy_r2B[ithread];
			cy_r2C->d0 = _cy_r2C[ithread];	cy_r2C->d1 = _cy_r2D[ithread];
			cy_r2E->d0 = _cy_r2E[ithread];	cy_r2E->d1 = _cy_r2F[ithread];
			cy_r30->d0 = _cy_r30[ithread];	cy_r30->d1 = _cy_r31[ithread];
			cy_r32->d0 = _cy_r32[ithread];	cy_r32->d1 = _cy_r33[ithread];
			cy_r34->d0 = _cy_r34[ithread];	cy_r34->d1 = _cy_r35[ithread];
			cy_r36->d0 = _cy_r36[ithread];	cy_r36->d1 = _cy_r37[ithread];
			cy_r38->d0 = _cy_r38[ithread];	cy_r38->d1 = _cy_r39[ithread];
			cy_r3A->d0 = _cy_r3A[ithread];	cy_r3A->d1 = _cy_r3B[ithread];
			cy_r3C->d0 = _cy_r3C[ithread];	cy_r3C->d1 = _cy_r3D[ithread];
			cy_r3E->d0 = _cy_r3E[ithread];	cy_r3E->d1 = _cy_r3F[ithread];
		  #endif	// AVX/SSE2?
		#else
			cy_r00 = _cy_r00[ithread];	cy_r20 = _cy_r20[ithread];
			cy_r01 = _cy_r01[ithread];	cy_r21 = _cy_r21[ithread];
			cy_r02 = _cy_r02[ithread];	cy_r22 = _cy_r22[ithread];
			cy_r03 = _cy_r03[ithread];	cy_r23 = _cy_r23[ithread];
			cy_r04 = _cy_r04[ithread];	cy_r24 = _cy_r24[ithread];
			cy_r05 = _cy_r05[ithread];	cy_r25 = _cy_r25[ithread];
			cy_r06 = _cy_r06[ithread];	cy_r26 = _cy_r26[ithread];
			cy_r07 = _cy_r07[ithread];	cy_r27 = _cy_r27[ithread];
			cy_r08 = _cy_r08[ithread];	cy_r28 = _cy_r28[ithread];
			cy_r09 = _cy_r09[ithread];	cy_r29 = _cy_r29[ithread];
			cy_r0A = _cy_r0A[ithread];	cy_r2A = _cy_r2A[ithread];
			cy_r0B = _cy_r0B[ithread];	cy_r2B = _cy_r2B[ithread];
			cy_r0C = _cy_r0C[ithread];	cy_r2C = _cy_r2C[ithread];
			cy_r0D = _cy_r0D[ithread];	cy_r2D = _cy_r2D[ithread];
			cy_r0E = _cy_r0E[ithread];	cy_r2E = _cy_r2E[ithread];
			cy_r0F = _cy_r0F[ithread];	cy_r2F = _cy_r2F[ithread];
			cy_r10 = _cy_r10[ithread];	cy_r30 = _cy_r30[ithread];
			cy_r11 = _cy_r11[ithread];	cy_r31 = _cy_r31[ithread];
			cy_r12 = _cy_r12[ithread];	cy_r32 = _cy_r32[ithread];
			cy_r13 = _cy_r13[ithread];	cy_r33 = _cy_r33[ithread];
			cy_r14 = _cy_r14[ithread];	cy_r34 = _cy_r34[ithread];
			cy_r15 = _cy_r15[ithread];	cy_r35 = _cy_r35[ithread];
			cy_r16 = _cy_r16[ithread];	cy_r36 = _cy_r36[ithread];
			cy_r17 = _cy_r17[ithread];	cy_r37 = _cy_r37[ithread];
			cy_r18 = _cy_r18[ithread];	cy_r38 = _cy_r38[ithread];
			cy_r19 = _cy_r19[ithread];	cy_r39 = _cy_r39[ithread];
			cy_r1A = _cy_r1A[ithread];	cy_r3A = _cy_r3A[ithread];
			cy_r1B = _cy_r1B[ithread];	cy_r3B = _cy_r3B[ithread];
			cy_r1C = _cy_r1C[ithread];	cy_r3C = _cy_r3C[ithread];
			cy_r1D = _cy_r1D[ithread];	cy_r3D = _cy_r3D[ithread];
			cy_r1E = _cy_r1E[ithread];	cy_r3E = _cy_r3E[ithread];
			cy_r1F = _cy_r1F[ithread];	cy_r3F = _cy_r3F[ithread];
		#endif
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
			/* init carries	*/
		#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
		  #ifdef USE_AVX
			cy_r00->d0 = _cy_r00[ithread];	cy_i00->d0 = _cy_i00[ithread];
			cy_r00->d1 = _cy_r01[ithread];	cy_i00->d1 = _cy_i01[ithread];
			cy_r00->d2 = _cy_r02[ithread];	cy_i00->d2 = _cy_i02[ithread];
			cy_r00->d3 = _cy_r03[ithread];	cy_i00->d3 = _cy_i03[ithread];
			cy_r04->d0 = _cy_r04[ithread];	cy_i04->d0 = _cy_i04[ithread];
			cy_r04->d1 = _cy_r05[ithread];	cy_i04->d1 = _cy_i05[ithread];
			cy_r04->d2 = _cy_r06[ithread];	cy_i04->d2 = _cy_i06[ithread];
			cy_r04->d3 = _cy_r07[ithread];	cy_i04->d3 = _cy_i07[ithread];
			cy_r08->d0 = _cy_r08[ithread];	cy_i08->d0 = _cy_i08[ithread];
			cy_r08->d1 = _cy_r09[ithread];	cy_i08->d1 = _cy_i09[ithread];
			cy_r08->d2 = _cy_r0A[ithread];	cy_i08->d2 = _cy_i0A[ithread];
			cy_r08->d3 = _cy_r0B[ithread];	cy_i08->d3 = _cy_i0B[ithread];
			cy_r0C->d0 = _cy_r0C[ithread];	cy_i0C->d0 = _cy_i0C[ithread];
			cy_r0C->d1 = _cy_r0D[ithread];	cy_i0C->d1 = _cy_i0D[ithread];
			cy_r0C->d2 = _cy_r0E[ithread];	cy_i0C->d2 = _cy_i0E[ithread];
			cy_r0C->d3 = _cy_r0F[ithread];	cy_i0C->d3 = _cy_i0F[ithread];
			cy_r10->d0 = _cy_r10[ithread];	cy_i10->d0 = _cy_i10[ithread];
			cy_r10->d1 = _cy_r11[ithread];	cy_i10->d1 = _cy_i11[ithread];
			cy_r10->d2 = _cy_r12[ithread];	cy_i10->d2 = _cy_i12[ithread];
			cy_r10->d3 = _cy_r13[ithread];	cy_i10->d3 = _cy_i13[ithread];
			cy_r14->d0 = _cy_r14[ithread];	cy_i14->d0 = _cy_i14[ithread];
			cy_r14->d1 = _cy_r15[ithread];	cy_i14->d1 = _cy_i15[ithread];
			cy_r14->d2 = _cy_r16[ithread];	cy_i14->d2 = _cy_i16[ithread];
			cy_r14->d3 = _cy_r17[ithread];	cy_i14->d3 = _cy_i17[ithread];
			cy_r18->d0 = _cy_r18[ithread];	cy_i18->d0 = _cy_i18[ithread];
			cy_r18->d1 = _cy_r19[ithread];	cy_i18->d1 = _cy_i19[ithread];
			cy_r18->d2 = _cy_r1A[ithread];	cy_i18->d2 = _cy_i1A[ithread];
			cy_r18->d3 = _cy_r1B[ithread];	cy_i18->d3 = _cy_i1B[ithread];
			cy_r1C->d0 = _cy_r1C[ithread];	cy_i1C->d0 = _cy_i1C[ithread];
			cy_r1C->d1 = _cy_r1D[ithread];	cy_i1C->d1 = _cy_i1D[ithread];
			cy_r1C->d2 = _cy_r1E[ithread];	cy_i1C->d2 = _cy_i1E[ithread];
			cy_r1C->d3 = _cy_r1F[ithread];	cy_i1C->d3 = _cy_i1F[ithread];
			cy_r20->d0 = _cy_r20[ithread];	cy_i20->d0 = _cy_i20[ithread];
			cy_r20->d1 = _cy_r21[ithread];	cy_i20->d1 = _cy_i21[ithread];
			cy_r20->d2 = _cy_r22[ithread];	cy_i20->d2 = _cy_i22[ithread];
			cy_r20->d3 = _cy_r23[ithread];	cy_i20->d3 = _cy_i23[ithread];
			cy_r24->d0 = _cy_r24[ithread];	cy_i24->d0 = _cy_i24[ithread];
			cy_r24->d1 = _cy_r25[ithread];	cy_i24->d1 = _cy_i25[ithread];
			cy_r24->d2 = _cy_r26[ithread];	cy_i24->d2 = _cy_i26[ithread];
			cy_r24->d3 = _cy_r27[ithread];	cy_i24->d3 = _cy_i27[ithread];
			cy_r28->d0 = _cy_r28[ithread];	cy_i28->d0 = _cy_i28[ithread];
			cy_r28->d1 = _cy_r29[ithread];	cy_i28->d1 = _cy_i29[ithread];
			cy_r28->d2 = _cy_r2A[ithread];	cy_i28->d2 = _cy_i2A[ithread];
			cy_r28->d3 = _cy_r2B[ithread];	cy_i28->d3 = _cy_i2B[ithread];
			cy_r2C->d0 = _cy_r2C[ithread];	cy_i2C->d0 = _cy_i2C[ithread];
			cy_r2C->d1 = _cy_r2D[ithread];	cy_i2C->d1 = _cy_i2D[ithread];
			cy_r2C->d2 = _cy_r2E[ithread];	cy_i2C->d2 = _cy_i2E[ithread];
			cy_r2C->d3 = _cy_r2F[ithread];	cy_i2C->d3 = _cy_i2F[ithread];
			cy_r30->d0 = _cy_r30[ithread];	cy_i30->d0 = _cy_i30[ithread];
			cy_r30->d1 = _cy_r31[ithread];	cy_i30->d1 = _cy_i31[ithread];
			cy_r30->d2 = _cy_r32[ithread];	cy_i30->d2 = _cy_i32[ithread];
			cy_r30->d3 = _cy_r33[ithread];	cy_i30->d3 = _cy_i33[ithread];
			cy_r34->d0 = _cy_r34[ithread];	cy_i34->d0 = _cy_i34[ithread];
			cy_r34->d1 = _cy_r35[ithread];	cy_i34->d1 = _cy_i35[ithread];
			cy_r34->d2 = _cy_r36[ithread];	cy_i34->d2 = _cy_i36[ithread];
			cy_r34->d3 = _cy_r37[ithread];	cy_i34->d3 = _cy_i37[ithread];
			cy_r38->d0 = _cy_r38[ithread];	cy_i38->d0 = _cy_i38[ithread];
			cy_r38->d1 = _cy_r39[ithread];	cy_i38->d1 = _cy_i39[ithread];
			cy_r38->d2 = _cy_r3A[ithread];	cy_i38->d2 = _cy_i3A[ithread];
			cy_r38->d3 = _cy_r3B[ithread];	cy_i38->d3 = _cy_i3B[ithread];
			cy_r3C->d0 = _cy_r3C[ithread];	cy_i3C->d0 = _cy_i3C[ithread];
			cy_r3C->d1 = _cy_r3D[ithread];	cy_i3C->d1 = _cy_i3D[ithread];
			cy_r3C->d2 = _cy_r3E[ithread];	cy_i3C->d2 = _cy_i3E[ithread];
			cy_r3C->d3 = _cy_r3F[ithread];	cy_i3C->d3 = _cy_i3F[ithread];
		  #else	// USE_SSE2
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky - see the comments around the associated carry macros for explanation:
			cy_r00->d0 = _cy_r00[ithread];	cy_r00->d1 = _cy_i00[ithread];
			cy_r02->d0 = _cy_r01[ithread];	cy_r02->d1 = _cy_i01[ithread];
			cy_r04->d0 = _cy_r02[ithread];	cy_r04->d1 = _cy_i02[ithread];
			cy_r06->d0 = _cy_r03[ithread];	cy_r06->d1 = _cy_i03[ithread];
			cy_r08->d0 = _cy_r04[ithread];	cy_r08->d1 = _cy_i04[ithread];
			cy_r0A->d0 = _cy_r05[ithread];	cy_r0A->d1 = _cy_i05[ithread];
			cy_r0C->d0 = _cy_r06[ithread];	cy_r0C->d1 = _cy_i06[ithread];
			cy_r0E->d0 = _cy_r07[ithread];	cy_r0E->d1 = _cy_i07[ithread];
			cy_r10->d0 = _cy_r08[ithread];	cy_r10->d1 = _cy_i08[ithread];
			cy_r12->d0 = _cy_r09[ithread];	cy_r12->d1 = _cy_i09[ithread];
			cy_r14->d0 = _cy_r0A[ithread];	cy_r14->d1 = _cy_i0A[ithread];
			cy_r16->d0 = _cy_r0B[ithread];	cy_r16->d1 = _cy_i0B[ithread];
			cy_r18->d0 = _cy_r0C[ithread];	cy_r18->d1 = _cy_i0C[ithread];
			cy_r1A->d0 = _cy_r0D[ithread];	cy_r1A->d1 = _cy_i0D[ithread];
			cy_r1C->d0 = _cy_r0E[ithread];	cy_r1C->d1 = _cy_i0E[ithread];
			cy_r1E->d0 = _cy_r0F[ithread];	cy_r1E->d1 = _cy_i0F[ithread];
			cy_r20->d0 = _cy_r10[ithread];	cy_r20->d1 = _cy_i10[ithread];
			cy_r22->d0 = _cy_r11[ithread];	cy_r22->d1 = _cy_i11[ithread];
			cy_r24->d0 = _cy_r12[ithread];	cy_r24->d1 = _cy_i12[ithread];
			cy_r26->d0 = _cy_r13[ithread];	cy_r26->d1 = _cy_i13[ithread];
			cy_r28->d0 = _cy_r14[ithread];	cy_r28->d1 = _cy_i14[ithread];
			cy_r2A->d0 = _cy_r15[ithread];	cy_r2A->d1 = _cy_i15[ithread];
			cy_r2C->d0 = _cy_r16[ithread];	cy_r2C->d1 = _cy_i16[ithread];
			cy_r2E->d0 = _cy_r17[ithread];	cy_r2E->d1 = _cy_i17[ithread];
			cy_r30->d0 = _cy_r18[ithread];	cy_r30->d1 = _cy_i18[ithread];
			cy_r32->d0 = _cy_r19[ithread];	cy_r32->d1 = _cy_i19[ithread];
			cy_r34->d0 = _cy_r1A[ithread];	cy_r34->d1 = _cy_i1A[ithread];
			cy_r36->d0 = _cy_r1B[ithread];	cy_r36->d1 = _cy_i1B[ithread];
			cy_r38->d0 = _cy_r1C[ithread];	cy_r38->d1 = _cy_i1C[ithread];
			cy_r3A->d0 = _cy_r1D[ithread];	cy_r3A->d1 = _cy_i1D[ithread];
			cy_r3C->d0 = _cy_r1E[ithread];	cy_r3C->d1 = _cy_i1E[ithread];
			cy_r3E->d0 = _cy_r1F[ithread];	cy_r3E->d1 = _cy_i1F[ithread];

			cy_i00->d0 = _cy_r20[ithread];	cy_i00->d1 = _cy_i20[ithread];
			cy_i02->d0 = _cy_r21[ithread];	cy_i02->d1 = _cy_i21[ithread];
			cy_i04->d0 = _cy_r22[ithread];	cy_i04->d1 = _cy_i22[ithread];
			cy_i06->d0 = _cy_r23[ithread];	cy_i06->d1 = _cy_i23[ithread];
			cy_i08->d0 = _cy_r24[ithread];	cy_i08->d1 = _cy_i24[ithread];
			cy_i0A->d0 = _cy_r25[ithread];	cy_i0A->d1 = _cy_i25[ithread];
			cy_i0C->d0 = _cy_r26[ithread];	cy_i0C->d1 = _cy_i26[ithread];
			cy_i0E->d0 = _cy_r27[ithread];	cy_i0E->d1 = _cy_i27[ithread];
			cy_i10->d0 = _cy_r28[ithread];	cy_i10->d1 = _cy_i28[ithread];
			cy_i12->d0 = _cy_r29[ithread];	cy_i12->d1 = _cy_i29[ithread];
			cy_i14->d0 = _cy_r2A[ithread];	cy_i14->d1 = _cy_i2A[ithread];
			cy_i16->d0 = _cy_r2B[ithread];	cy_i16->d1 = _cy_i2B[ithread];
			cy_i18->d0 = _cy_r2C[ithread];	cy_i18->d1 = _cy_i2C[ithread];
			cy_i1A->d0 = _cy_r2D[ithread];	cy_i1A->d1 = _cy_i2D[ithread];
			cy_i1C->d0 = _cy_r2E[ithread];	cy_i1C->d1 = _cy_i2E[ithread];
			cy_i1E->d0 = _cy_r2F[ithread];	cy_i1E->d1 = _cy_i2F[ithread];
			cy_i20->d0 = _cy_r30[ithread];	cy_i20->d1 = _cy_i30[ithread];
			cy_i22->d0 = _cy_r31[ithread];	cy_i22->d1 = _cy_i31[ithread];
			cy_i24->d0 = _cy_r32[ithread];	cy_i24->d1 = _cy_i32[ithread];
			cy_i26->d0 = _cy_r33[ithread];	cy_i26->d1 = _cy_i33[ithread];
			cy_i28->d0 = _cy_r34[ithread];	cy_i28->d1 = _cy_i34[ithread];
			cy_i2A->d0 = _cy_r35[ithread];	cy_i2A->d1 = _cy_i35[ithread];
			cy_i2C->d0 = _cy_r36[ithread];	cy_i2C->d1 = _cy_i36[ithread];
			cy_i2E->d0 = _cy_r37[ithread];	cy_i2E->d1 = _cy_i37[ithread];
			cy_i30->d0 = _cy_r38[ithread];	cy_i30->d1 = _cy_i38[ithread];
			cy_i32->d0 = _cy_r39[ithread];	cy_i32->d1 = _cy_i39[ithread];
			cy_i34->d0 = _cy_r3A[ithread];	cy_i34->d1 = _cy_i3A[ithread];
			cy_i36->d0 = _cy_r3B[ithread];	cy_i36->d1 = _cy_i3B[ithread];
			cy_i38->d0 = _cy_r3C[ithread];	cy_i38->d1 = _cy_i3C[ithread];
			cy_i3A->d0 = _cy_r3D[ithread];	cy_i3A->d1 = _cy_i3D[ithread];
			cy_i3C->d0 = _cy_r3E[ithread];	cy_i3C->d1 = _cy_i3E[ithread];
			cy_i3E->d0 = _cy_r3F[ithread];	cy_i3E->d1 = _cy_i3F[ithread];
		  #endif	// AVX/SSE2?
		#else
			cy_r00 = _cy_r00[ithread];	cy_i00 = _cy_i00[ithread];
			cy_r01 = _cy_r01[ithread];	cy_i01 = _cy_i01[ithread];
			cy_r02 = _cy_r02[ithread];	cy_i02 = _cy_i02[ithread];
			cy_r03 = _cy_r03[ithread];	cy_i03 = _cy_i03[ithread];
			cy_r04 = _cy_r04[ithread];	cy_i04 = _cy_i04[ithread];
			cy_r05 = _cy_r05[ithread];	cy_i05 = _cy_i05[ithread];
			cy_r06 = _cy_r06[ithread];	cy_i06 = _cy_i06[ithread];
			cy_r07 = _cy_r07[ithread];	cy_i07 = _cy_i07[ithread];
			cy_r08 = _cy_r08[ithread];	cy_i08 = _cy_i08[ithread];
			cy_r09 = _cy_r09[ithread];	cy_i09 = _cy_i09[ithread];
			cy_r0A = _cy_r0A[ithread];	cy_i0A = _cy_i0A[ithread];
			cy_r0B = _cy_r0B[ithread];	cy_i0B = _cy_i0B[ithread];
			cy_r0C = _cy_r0C[ithread];	cy_i0C = _cy_i0C[ithread];
			cy_r0D = _cy_r0D[ithread];	cy_i0D = _cy_i0D[ithread];
			cy_r0E = _cy_r0E[ithread];	cy_i0E = _cy_i0E[ithread];
			cy_r0F = _cy_r0F[ithread];	cy_i0F = _cy_i0F[ithread];
			cy_r10 = _cy_r10[ithread];	cy_i10 = _cy_i10[ithread];
			cy_r11 = _cy_r11[ithread];	cy_i11 = _cy_i11[ithread];
			cy_r12 = _cy_r12[ithread];	cy_i12 = _cy_i12[ithread];
			cy_r13 = _cy_r13[ithread];	cy_i13 = _cy_i13[ithread];
			cy_r14 = _cy_r14[ithread];	cy_i14 = _cy_i14[ithread];
			cy_r15 = _cy_r15[ithread];	cy_i15 = _cy_i15[ithread];
			cy_r16 = _cy_r16[ithread];	cy_i16 = _cy_i16[ithread];
			cy_r17 = _cy_r17[ithread];	cy_i17 = _cy_i17[ithread];
			cy_r18 = _cy_r18[ithread];	cy_i18 = _cy_i18[ithread];
			cy_r19 = _cy_r19[ithread];	cy_i19 = _cy_i19[ithread];
			cy_r1A = _cy_r1A[ithread];	cy_i1A = _cy_i1A[ithread];
			cy_r1B = _cy_r1B[ithread];	cy_i1B = _cy_i1B[ithread];
			cy_r1C = _cy_r1C[ithread];	cy_i1C = _cy_i1C[ithread];
			cy_r1D = _cy_r1D[ithread];	cy_i1D = _cy_i1D[ithread];
			cy_r1E = _cy_r1E[ithread];	cy_i1E = _cy_i1E[ithread];
			cy_r1F = _cy_r1F[ithread];	cy_i1F = _cy_i1F[ithread];

			cy_r20 = _cy_r20[ithread];	cy_i20 = _cy_i20[ithread];
			cy_r21 = _cy_r21[ithread];	cy_i21 = _cy_i21[ithread];
			cy_r22 = _cy_r22[ithread];	cy_i22 = _cy_i22[ithread];
			cy_r23 = _cy_r23[ithread];	cy_i23 = _cy_i23[ithread];
			cy_r24 = _cy_r24[ithread];	cy_i24 = _cy_i24[ithread];
			cy_r25 = _cy_r25[ithread];	cy_i25 = _cy_i25[ithread];
			cy_r26 = _cy_r26[ithread];	cy_i26 = _cy_i26[ithread];
			cy_r27 = _cy_r27[ithread];	cy_i27 = _cy_i27[ithread];
			cy_r28 = _cy_r28[ithread];	cy_i28 = _cy_i28[ithread];
			cy_r29 = _cy_r29[ithread];	cy_i29 = _cy_i29[ithread];
			cy_r2A = _cy_r2A[ithread];	cy_i2A = _cy_i2A[ithread];
			cy_r2B = _cy_r2B[ithread];	cy_i2B = _cy_i2B[ithread];
			cy_r2C = _cy_r2C[ithread];	cy_i2C = _cy_i2C[ithread];
			cy_r2D = _cy_r2D[ithread];	cy_i2D = _cy_i2D[ithread];
			cy_r2E = _cy_r2E[ithread];	cy_i2E = _cy_i2E[ithread];
			cy_r2F = _cy_r2F[ithread];	cy_i2F = _cy_i2F[ithread];
			cy_r30 = _cy_r30[ithread];	cy_i30 = _cy_i30[ithread];
			cy_r31 = _cy_r31[ithread];	cy_i31 = _cy_i31[ithread];
			cy_r32 = _cy_r32[ithread];	cy_i32 = _cy_i32[ithread];
			cy_r33 = _cy_r33[ithread];	cy_i33 = _cy_i33[ithread];
			cy_r34 = _cy_r34[ithread];	cy_i34 = _cy_i34[ithread];
			cy_r35 = _cy_r35[ithread];	cy_i35 = _cy_i35[ithread];
			cy_r36 = _cy_r36[ithread];	cy_i36 = _cy_i36[ithread];
			cy_r37 = _cy_r37[ithread];	cy_i37 = _cy_i37[ithread];
			cy_r38 = _cy_r38[ithread];	cy_i38 = _cy_i38[ithread];
			cy_r39 = _cy_r39[ithread];	cy_i39 = _cy_i39[ithread];
			cy_r3A = _cy_r3A[ithread];	cy_i3A = _cy_i3A[ithread];
			cy_r3B = _cy_r3B[ithread];	cy_i3B = _cy_i3B[ithread];
			cy_r3C = _cy_r3C[ithread];	cy_i3C = _cy_i3C[ithread];
			cy_r3D = _cy_r3D[ithread];	cy_i3D = _cy_i3D[ithread];
			cy_r3E = _cy_r3E[ithread];	cy_i3E = _cy_i3E[ithread];
			cy_r3F = _cy_r3F[ithread];	cy_i3F = _cy_i3F[ithread];
		#endif
		}

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
			/* In SIMD mode, data are arranged in [re_0,...,re_n-1,im_0,...,im_n-1] groups, not the usual [re_0,im_0],...,[re_n-1,im_n-1] pairs.
			Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
			but to point to the proper real datum, we need to index-map e.g. [0,1,2,3] ==> [0,2,1,3] in 2-way SIMD mode.
			(But only ever need to explicitly do this in debug mode).
			*/
			for(j = jstart; j < jhi; j += stride)
			{
				j1 =  j;
				j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1 + RE_IM_STRIDE;

		/*...The radix-64 DIT pass is here:	*/

		// SIMD code only available for 64-bit GCC build - others simply use scalar DFT macros with SIMD-compatible data layout
		#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)

		/* Gather the needed data (64 64-bit complex, i.e. 128 64-bit reals) and do 8 twiddleless length-8 subtransforms: */
// Because of roots-sign diffs between SSE2 and scalar macros, 1/7 2/6 3/5 swapped in DIT_0TWIDDLE outputs!
			add0 = &a[j1]      ; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r00, isrt2)

			add0 = &a[j1] + p08; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r10, isrt2)

			add0 = &a[j1] + p10; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r20, isrt2)

			add0 = &a[j1] + p18; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r30, isrt2)

			add0 = &a[j1] + p20; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r40, isrt2)

			add0 = &a[j1] + p28; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r50, isrt2)

			add0 = &a[j1] + p30; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r60, isrt2)

			add0 = &a[j1] + p38; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r70, isrt2)

// Because of roots-sign diffs between SSE2 and scalar macros, 1/7 2/6 3/5 swapped in DIT_0TWIDDLE outputs!

		/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dit_pass1 for details: */
		// Note: 1st of the 15 sincos args in each call to SSE2_RADIX8_DIT_TWIDDLE_OOP is the basic isrt2 needed for
		// radix-8. This is a workaround of GCC's 30-arg limit for inline ASM macros, which proves a royal pain here.
		//...and another kludge for the 30-arg limit: put a copy of (vec_dbl)2.0 into the first of each s1p*r outputs:
			VEC_DBL_INIT(s1p01r,2.0);
			VEC_DBL_INIT(s1p02r,2.0);
			VEC_DBL_INIT(s1p03r,2.0);
			VEC_DBL_INIT(s1p04r,2.0);
			VEC_DBL_INIT(s1p05r,2.0);
			VEC_DBL_INIT(s1p06r,2.0);
			VEC_DBL_INIT(s1p07r,2.0);

			// Block 0: jt = j1;	jp = j2;
			SSE2_RADIX8_DIT_0TWIDDLE_OOP(	// This outputs o[07654321], so reverse o-index order of latter 7 outputs
				r00,r10,r20,r30,r40,r50,r60,r70,
				s1p00r,s1p38r,s1p30r,s1p28r,s1p20r,s1p18r,s1p10r,s1p08r, isrt2
			);
		/* 0-index block has all-unity twiddles: **
			VEC_DBL_INIT(s1p00r,2.0);
			SSE2_RADIX8_DIT_TWIDDLE_OOP(
				r00,r10,r20,r30,r40,r50,r60,r70,
				s1p00r,s1p20r,s1p10r,s1p30r,s1p08r,s1p28r,s1p18r,s1p38r,
				cc0,ss0,cc0,ss0,cc0,ss0,cc0,ss0,cc0,ss0,cc0,ss0,cc0,ss0
			);
		*/
			// Block 4: jt = j1 + p04;	jp = j2 + p04;
			SSE2_RADIX8_DIT_TWIDDLE_OOP(
				r08,r18,r28,r38,r48,r58,r68,r78,
				s1p04r,s1p24r,s1p14r,s1p34r,s1p0cr,s1p2cr,s1p1cr,s1p3cr,
				ss0,cc0, isrt2,isrt2, nisrt2,isrt2, cc4,ss4, nss4,cc4, ss4,cc4, ncc4,ss4
			);
			// Block 2: jt = j1 + p02;	jp = j2 + p02;
			SSE2_RADIX8_DIT_TWIDDLE_OOP(
				r0C,r1C,r2C,r3C,r4C,r5C,r6C,r7C,	// 2/6 swap ==> r*4 / r*C swap
				s1p02r,s1p22r,s1p12r,s1p32r,s1p0ar,s1p2ar,s1p1ar,s1p3ar,
				isrt2,isrt2,cc4,ss4,ss4,cc4,cc2,ss2,ss6,cc6,cc6,ss6,ss2,cc2
			);
			// Block 6: jt = j1 + p06;	jp = j2 + p06;
			SSE2_RADIX8_DIT_TWIDDLE_OOP(
				r04,r14,r24,r34,r44,r54,r64,r74,	// 2/6 swap ==> r*4 / r*C swap
				s1p06r,s1p26r,s1p16r,s1p36r,s1p0er,s1p2er,s1p1er,s1p3er,
				nisrt2,isrt2,ss4,cc4,ncc4,nss4,cc6,ss6,ncc2,ss2,nss2,cc2,nss6,ncc6
			);
			// Block 1: jt = j1 + p01;	jp = j2 + p01;
			SSE2_RADIX8_DIT_TWIDDLE_OOP(
				r0E,r1E,r2E,r3E,r4E,r5E,r6E,r7E,	// 1/7 swap ==> r*2 / r*E swap
				s1p01r,s1p21r,s1p11r,s1p31r,s1p09r,s1p29r,s1p19r,s1p39r,
				cc4,ss4,cc2,ss2,cc6,ss6,cc1,ss1,cc5,ss5,cc3,ss3,cc7,ss7
			);
			// Block 5: jt = j1 + p05;	jp = j2 + p05;
			SSE2_RADIX8_DIT_TWIDDLE_OOP(
				r06,r16,r26,r36,r46,r56,r66,r76,
				s1p05r,s1p25r,s1p15r,s1p35r,s1p0dr,s1p2dr,s1p1dr,s1p3dr,
				nss4,cc4,ss6,cc6,ncc2,ss2,cc5,ss5,ncc7,ss7,ss1,cc1,ncc3,nss3
			);
			// Block 3: jt = j1 + p03;	jp = j2 + p03;
			SSE2_RADIX8_DIT_TWIDDLE_OOP(
				r0A,r1A,r2A,r3A,r4A,r5A,r6A,r7A,	// 3/5 swap ==> r*6 / r*A swap
				s1p03r,s1p23r,s1p13r,s1p33r,s1p0br,s1p2br,s1p1br,s1p3br,
				ss4,cc4,cc6,ss6,nss2,cc2,cc3,ss3,ss1,cc1,ss7,cc7,nss5,cc5
			);
			// Block 7: jt = j1 + p07;	jp = j2 + p07;
			SSE2_RADIX8_DIT_TWIDDLE_OOP(
				r02,r12,r22,r32,r42,r52,r62,r72,	// 1/7 swap ==> r*2 / r*E swap
				s1p07r,s1p27r,s1p17r,s1p37r,s1p0fr,s1p2fr,s1p1fr,s1p3fr,
				ncc4,ss4,ss2,cc2,nss6,ncc6,cc7,ss7,ncc3,nss3,nss5,cc5,ss1,ncc1
			);

		#else	// USE_SSE2 = False, or non-64-bit-GCC:

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

		/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dit_pass1 for details: */

			// Block 0: jt = j1;	jp = j2;
			/* 0-index block has all-unity twiddles: Remember, the twiddleless DIT also bit-reverses its outputs, so a_p* terms appear in-order here: */
			RADIX_08_DIT_OOP(
				t00,t01,t10,t11,t20,t21,t30,t31,t40,t41,t50,t51,t60,t61,t70,t71,
				ar_p00,ai_p00,ar_p08,ai_p08,ar_p10,ai_p10,ar_p18,ai_p18,ar_p20,ai_p20,ar_p28,ai_p28,ar_p30,ai_p30,ar_p38,ai_p38
			);
			// Block 4: jt = j1 + p04;	jp = j2 + p04;
			RADIX_08_DIT_TWIDDLE_OOP(
				t08,t09,t18,t19,t28,t29,t38,t39,t48,t49,t58,t59,t68,t69,t78,t79,
				ar_p04,ai_p04,ar_p24,ai_p24,ar_p14,ai_p14,ar_p34,ai_p34,ar_p0c,ai_p0c,ar_p2c,ai_p2c,ar_p1c,ai_p1c,ar_p3c,ai_p3c,
				0,1,ISRT2,ISRT2,-ISRT2,ISRT2,c64_4,s64_4,-s64_4,c64_4,s64_4,c64_4,-c64_4,s64_4
			);
			// Block 2: jt = j1 + p02;	jp = j2 + p02;
			RADIX_08_DIT_TWIDDLE_OOP(
				t04,t05,t14,t15,t24,t25,t34,t35,t44,t45,t54,t55,t64,t65,t74,t75,
				ar_p02,ai_p02,ar_p22,ai_p22,ar_p12,ai_p12,ar_p32,ai_p32,ar_p0a,ai_p0a,ar_p2a,ai_p2a,ar_p1a,ai_p1a,ar_p3a,ai_p3a,
				ISRT2,ISRT2,c64_4,s64_4,s64_4,c64_4,c64_2,s64_2,s64_6,c64_6,c64_6,s64_6,s64_2,c64_2
			);
			// Block 6: jt = j1 + p06;	jp = j2 + p06;
			RADIX_08_DIT_TWIDDLE_OOP(
				t0C,t0D,t1C,t1D,t2C,t2D,t3C,t3D,t4C,t4D,t5C,t5D,t6C,t6D,t7C,t7D,
				ar_p06,ai_p06,ar_p26,ai_p26,ar_p16,ai_p16,ar_p36,ai_p36,ar_p0e,ai_p0e,ar_p2e,ai_p2e,ar_p1e,ai_p1e,ar_p3e,ai_p3e,
				-ISRT2,ISRT2,s64_4,c64_4,-c64_4,-s64_4,c64_6,s64_6,-c64_2,s64_2,-s64_2,c64_2,-s64_6,-c64_6
			);
			// Block 1: jt = j1 + p01;	jp = j2 + p01;
			RADIX_08_DIT_TWIDDLE_OOP(
				t02,t03,t12,t13,t22,t23,t32,t33,t42,t43,t52,t53,t62,t63,t72,t73,
				ar_p01,ai_p01,ar_p21,ai_p21,ar_p11,ai_p11,ar_p31,ai_p31,ar_p09,ai_p09,ar_p29,ai_p29,ar_p19,ai_p19,ar_p39,ai_p39,
				c64_4,s64_4,c64_2,s64_2,c64_6,s64_6,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7
			);
			// Block 5: jt = j1 + p05;	jp = j2 + p05;
			RADIX_08_DIT_TWIDDLE_OOP(
				t0A,t0B,t1A,t1B,t2A,t2B,t3A,t3B,t4A,t4B,t5A,t5B,t6A,t6B,t7A,t7B,
				ar_p05,ai_p05,ar_p25,ai_p25,ar_p15,ai_p15,ar_p35,ai_p35,ar_p0d,ai_p0d,ar_p2d,ai_p2d,ar_p1d,ai_p1d,ar_p3d,ai_p3d,
				-s64_4,c64_4,s64_6,c64_6,-c64_2,s64_2,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3
			);
			// Block 3: jt = j1 + p03;	jp = j2 + p03;
			RADIX_08_DIT_TWIDDLE_OOP(
				t06,t07,t16,t17,t26,t27,t36,t37,t46,t47,t56,t57,t66,t67,t76,t77,
				ar_p03,ai_p03,ar_p23,ai_p23,ar_p13,ai_p13,ar_p33,ai_p33,ar_p0b,ai_p0b,ar_p2b,ai_p2b,ar_p1b,ai_p1b,ar_p3b,ai_p3b,
				s64_4,c64_4,c64_6,s64_6,-s64_2,c64_2,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5
			);
			// Block 7: jt = j1 + p07;	jp = j2 + p07;
			RADIX_08_DIT_TWIDDLE_OOP(
				t0E,t0F,t1E,t1F,t2E,t2F,t3E,t3F,t4E,t4F,t5E,t5F,t6E,t6F,t7E,t7F,
				ar_p07,ai_p07,ar_p27,ai_p27,ar_p17,ai_p17,ar_p37,ai_p37,ar_p0f,ai_p0f,ar_p2f,ai_p2f,ar_p1f,ai_p1f,ar_p3f,ai_p3f,
				-c64_4,s64_4,s64_2,c64_2,-s64_6,-c64_6,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1
			);

		#endif	// USE_SSE2 ?

		/*...Now do the carries. Since the outputs would
		normally be getting dispatched to 32 separate blocks of the A-array, we need 32 separate carries.	*/

		/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/
			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
			#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
			  #ifdef USE_AVX

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

				l= j & (nwt-1);						tmp = half_arr + 64;	/* ptr to local storage for the doubled wtl,wtn terms: */
				n_minus_sil  ->d0 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d0 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d0 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d0 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+2) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d1 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d1 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d1 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d1 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+4) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d2 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d2 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d2 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d2 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+6) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d3 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d3 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d3 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d3 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				AVX_cmplx_carry_norm_pow2_errcheck0_X4(s1p00r,add1,add2,add3,cy_r00,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p04r,add1,add2,add3,cy_r04,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p08r,add1,add2,add3,cy_r08,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p0cr,add1,add2,add3,cy_r0C,bjmodn0C,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p10r,add1,add2,add3,cy_r10,bjmodn10,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p14r,add1,add2,add3,cy_r14,bjmodn14,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p18r,add1,add2,add3,cy_r18,bjmodn18,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p1cr,add1,add2,add3,cy_r1C,bjmodn1C,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p20r,add1,add2,add3,cy_r20,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p24r,add1,add2,add3,cy_r24,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p28r,add1,add2,add3,cy_r28,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p2cr,add1,add2,add3,cy_r2C,bjmodn2C,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p30r,add1,add2,add3,cy_r30,bjmodn30,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p34r,add1,add2,add3,cy_r34,bjmodn34,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p38r,add1,add2,add3,cy_r38,bjmodn38,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p3cr,add1,add2,add3,cy_r3C,bjmodn3C,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);

				co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

			  #else	// USE_SSE2

				l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 32 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				ctmp->re = ctmp->im = wtl;		++ctmp;
				ctmp->re = ctmp->im = wtn;		++ctmp;
				ctmp->re = ctmp->im = wtlp1;	++ctmp;
				ctmp->re = ctmp->im = wtnm1;

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_pow2_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p08r,add1,add2,add3,cy_r08,cy_r0A,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p0cr,add1,add2,add3,cy_r0C,cy_r0E,bjmodn0C,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p10r,add1,add2,add3,cy_r10,cy_r12,bjmodn10,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p14r,add1,add2,add3,cy_r14,cy_r16,bjmodn14,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p18r,add1,add2,add3,cy_r18,cy_r1A,bjmodn18,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p1cr,add1,add2,add3,cy_r1C,cy_r1E,bjmodn1C,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p28r,add1,add2,add3,cy_r28,cy_r2A,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p2cr,add1,add2,add3,cy_r2C,cy_r2E,bjmodn2C,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p30r,add1,add2,add3,cy_r30,cy_r32,bjmodn30,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p34r,add1,add2,add3,cy_r34,cy_r36,bjmodn34,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p38r,add1,add2,add3,cy_r38,cy_r3A,bjmodn38,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p3cr,add1,add2,add3,cy_r3C,cy_r3E,bjmodn3C,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			   #else
				SSE2_cmplx_carry_norm_pow2_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r0A,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p0cr,add1,add2,add3,cy_r0C,cy_r0E,bjmodn0C,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p10r,add1,add2,add3,cy_r10,cy_r12,bjmodn10,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p14r,add1,add2,add3,cy_r14,cy_r16,bjmodn14,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p18r,add1,add2,add3,cy_r18,cy_r1A,bjmodn18,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p1cr,add1,add2,add3,cy_r1C,cy_r1E,bjmodn1C,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p28r,add1,add2,add3,cy_r28,cy_r2A,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p2cr,add1,add2,add3,cy_r2C,cy_r2E,bjmodn2C,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p30r,add1,add2,add3,cy_r30,cy_r32,bjmodn30,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p34r,add1,add2,add3,cy_r34,cy_r36,bjmodn34,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p38r,add1,add2,add3,cy_r38,cy_r3A,bjmodn38,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p3cr,add1,add2,add3,cy_r3C,cy_r3E,bjmodn3C,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			   #endif

				l= (j+2) & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				ctmp->re = ctmp->im = wtl;		++ctmp;
				ctmp->re = ctmp->im = wtn;		++ctmp;
				ctmp->re = ctmp->im = wtlp1;	++ctmp;
				ctmp->re = ctmp->im = wtnm1;

				co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r0A,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p0cr,add1,add2,     cy_r0C,cy_r0E,bjmodn0C,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p10r,add1,add2,     cy_r10,cy_r12,bjmodn10,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p14r,add1,add2,     cy_r14,cy_r16,bjmodn14,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p18r,add1,add2,     cy_r18,cy_r1A,bjmodn18,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p1cr,add1,add2,     cy_r1C,cy_r1E,bjmodn1C,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p28r,add1,add2,     cy_r28,cy_r2A,bjmodn28,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p2cr,add1,add2,     cy_r2C,cy_r2E,bjmodn2C,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p30r,add1,add2,     cy_r30,cy_r32,bjmodn30,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p34r,add1,add2,     cy_r34,cy_r36,bjmodn34,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p38r,add1,add2,     cy_r38,cy_r3A,bjmodn38,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p3cr,add1,add2,     cy_r3C,cy_r3E,bjmodn3C,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			   #else
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p08r,add1,add2,     cy_r08,cy_r0A,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p0cr,add1,add2,     cy_r0C,cy_r0E,bjmodn0C,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p10r,add1,add2,     cy_r10,cy_r12,bjmodn10,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p14r,add1,add2,     cy_r14,cy_r16,bjmodn14,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p18r,add1,add2,     cy_r18,cy_r1A,bjmodn18,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p1cr,add1,add2,     cy_r1C,cy_r1E,bjmodn1C,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p28r,add1,add2,     cy_r28,cy_r2A,bjmodn28,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p2cr,add1,add2,     cy_r2C,cy_r2E,bjmodn2C,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p30r,add1,add2,     cy_r30,cy_r32,bjmodn30,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p34r,add1,add2,     cy_r34,cy_r36,bjmodn34,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p38r,add1,add2,     cy_r38,cy_r3A,bjmodn38,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p3cr,add1,add2,     cy_r3C,cy_r3E,bjmodn3C,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			   #endif

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

			  #endif	// AVX/SSE2?

			#else	// Scalar-double mode:

				l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 32 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				/*...set0 is slightly different from others:	*/
			   cmplx_carry_norm_pow2_errcheck0(ar_p00,ai_p00,cy_r00,bjmodn00);
				cmplx_carry_norm_pow2_errcheck(ar_p01,ai_p01,cy_r01,bjmodn01,0x01);
				cmplx_carry_norm_pow2_errcheck(ar_p02,ai_p02,cy_r02,bjmodn02,0x02);
				cmplx_carry_norm_pow2_errcheck(ar_p03,ai_p03,cy_r03,bjmodn03,0x03);
				cmplx_carry_norm_pow2_errcheck(ar_p04,ai_p04,cy_r04,bjmodn04,0x04);
				cmplx_carry_norm_pow2_errcheck(ar_p05,ai_p05,cy_r05,bjmodn05,0x05);
				cmplx_carry_norm_pow2_errcheck(ar_p06,ai_p06,cy_r06,bjmodn06,0x06);
				cmplx_carry_norm_pow2_errcheck(ar_p07,ai_p07,cy_r07,bjmodn07,0x07);
				cmplx_carry_norm_pow2_errcheck(ar_p08,ai_p08,cy_r08,bjmodn08,0x08);
				cmplx_carry_norm_pow2_errcheck(ar_p09,ai_p09,cy_r09,bjmodn09,0x09);
				cmplx_carry_norm_pow2_errcheck(ar_p0a,ai_p0a,cy_r0A,bjmodn0A,0x0A);
				cmplx_carry_norm_pow2_errcheck(ar_p0b,ai_p0b,cy_r0B,bjmodn0B,0x0B);
				cmplx_carry_norm_pow2_errcheck(ar_p0c,ai_p0c,cy_r0C,bjmodn0C,0x0C);
				cmplx_carry_norm_pow2_errcheck(ar_p0d,ai_p0d,cy_r0D,bjmodn0D,0x0D);
				cmplx_carry_norm_pow2_errcheck(ar_p0e,ai_p0e,cy_r0E,bjmodn0E,0x0E);
				cmplx_carry_norm_pow2_errcheck(ar_p0f,ai_p0f,cy_r0F,bjmodn0F,0x0F);
				cmplx_carry_norm_pow2_errcheck(ar_p10,ai_p10,cy_r10,bjmodn10,0x10);
				cmplx_carry_norm_pow2_errcheck(ar_p11,ai_p11,cy_r11,bjmodn11,0x11);
				cmplx_carry_norm_pow2_errcheck(ar_p12,ai_p12,cy_r12,bjmodn12,0x12);
				cmplx_carry_norm_pow2_errcheck(ar_p13,ai_p13,cy_r13,bjmodn13,0x13);
				cmplx_carry_norm_pow2_errcheck(ar_p14,ai_p14,cy_r14,bjmodn14,0x14);
				cmplx_carry_norm_pow2_errcheck(ar_p15,ai_p15,cy_r15,bjmodn15,0x15);
				cmplx_carry_norm_pow2_errcheck(ar_p16,ai_p16,cy_r16,bjmodn16,0x16);
				cmplx_carry_norm_pow2_errcheck(ar_p17,ai_p17,cy_r17,bjmodn17,0x17);
				cmplx_carry_norm_pow2_errcheck(ar_p18,ai_p18,cy_r18,bjmodn18,0x18);
				cmplx_carry_norm_pow2_errcheck(ar_p19,ai_p19,cy_r19,bjmodn19,0x19);
				cmplx_carry_norm_pow2_errcheck(ar_p1a,ai_p1a,cy_r1A,bjmodn1A,0x1A);
				cmplx_carry_norm_pow2_errcheck(ar_p1b,ai_p1b,cy_r1B,bjmodn1B,0x1B);
				cmplx_carry_norm_pow2_errcheck(ar_p1c,ai_p1c,cy_r1C,bjmodn1C,0x1C);
				cmplx_carry_norm_pow2_errcheck(ar_p1d,ai_p1d,cy_r1D,bjmodn1D,0x1D);
				cmplx_carry_norm_pow2_errcheck(ar_p1e,ai_p1e,cy_r1E,bjmodn1E,0x1E);
				cmplx_carry_norm_pow2_errcheck(ar_p1f,ai_p1f,cy_r1F,bjmodn1F,0x1F);
				cmplx_carry_norm_pow2_errcheck(ar_p20,ai_p20,cy_r20,bjmodn20,0x20);
				cmplx_carry_norm_pow2_errcheck(ar_p21,ai_p21,cy_r21,bjmodn21,0x21);
				cmplx_carry_norm_pow2_errcheck(ar_p22,ai_p22,cy_r22,bjmodn22,0x22);
				cmplx_carry_norm_pow2_errcheck(ar_p23,ai_p23,cy_r23,bjmodn23,0x23);
				cmplx_carry_norm_pow2_errcheck(ar_p24,ai_p24,cy_r24,bjmodn24,0x24);
				cmplx_carry_norm_pow2_errcheck(ar_p25,ai_p25,cy_r25,bjmodn25,0x25);
				cmplx_carry_norm_pow2_errcheck(ar_p26,ai_p26,cy_r26,bjmodn26,0x26);
				cmplx_carry_norm_pow2_errcheck(ar_p27,ai_p27,cy_r27,bjmodn27,0x27);
				cmplx_carry_norm_pow2_errcheck(ar_p28,ai_p28,cy_r28,bjmodn28,0x28);
				cmplx_carry_norm_pow2_errcheck(ar_p29,ai_p29,cy_r29,bjmodn29,0x29);
				cmplx_carry_norm_pow2_errcheck(ar_p2a,ai_p2a,cy_r2A,bjmodn2A,0x2A);
				cmplx_carry_norm_pow2_errcheck(ar_p2b,ai_p2b,cy_r2B,bjmodn2B,0x2B);
				cmplx_carry_norm_pow2_errcheck(ar_p2c,ai_p2c,cy_r2C,bjmodn2C,0x2C);
				cmplx_carry_norm_pow2_errcheck(ar_p2d,ai_p2d,cy_r2D,bjmodn2D,0x2D);
				cmplx_carry_norm_pow2_errcheck(ar_p2e,ai_p2e,cy_r2E,bjmodn2E,0x2E);
				cmplx_carry_norm_pow2_errcheck(ar_p2f,ai_p2f,cy_r2F,bjmodn2F,0x2F);
				cmplx_carry_norm_pow2_errcheck(ar_p30,ai_p30,cy_r30,bjmodn30,0x30);
				cmplx_carry_norm_pow2_errcheck(ar_p31,ai_p31,cy_r31,bjmodn31,0x31);
				cmplx_carry_norm_pow2_errcheck(ar_p32,ai_p32,cy_r32,bjmodn32,0x32);
				cmplx_carry_norm_pow2_errcheck(ar_p33,ai_p33,cy_r33,bjmodn33,0x33);
				cmplx_carry_norm_pow2_errcheck(ar_p34,ai_p34,cy_r34,bjmodn34,0x34);
				cmplx_carry_norm_pow2_errcheck(ar_p35,ai_p35,cy_r35,bjmodn35,0x35);
				cmplx_carry_norm_pow2_errcheck(ar_p36,ai_p36,cy_r36,bjmodn36,0x36);
				cmplx_carry_norm_pow2_errcheck(ar_p37,ai_p37,cy_r37,bjmodn37,0x37);
				cmplx_carry_norm_pow2_errcheck(ar_p38,ai_p38,cy_r38,bjmodn38,0x38);
				cmplx_carry_norm_pow2_errcheck(ar_p39,ai_p39,cy_r39,bjmodn39,0x39);
				cmplx_carry_norm_pow2_errcheck(ar_p3a,ai_p3a,cy_r3A,bjmodn3A,0x3A);
				cmplx_carry_norm_pow2_errcheck(ar_p3b,ai_p3b,cy_r3B,bjmodn3B,0x3B);
				cmplx_carry_norm_pow2_errcheck(ar_p3c,ai_p3c,cy_r3C,bjmodn3C,0x3C);
				cmplx_carry_norm_pow2_errcheck(ar_p3d,ai_p3d,cy_r3D,bjmodn3D,0x3D);
				cmplx_carry_norm_pow2_errcheck(ar_p3e,ai_p3e,cy_r3E,bjmodn3E,0x3E);
				cmplx_carry_norm_pow2_errcheck(ar_p3f,ai_p3f,cy_r3F,bjmodn3F,0x3F);

				i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
				co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
						 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

			#endif	// USE_AVX?
			}
			else	/* Fermat-mod carry in SIMD mode */
			{
			#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
			  #ifdef USE_AVX

				// For a description of the data movement for Fermat-mod carries in SSE2 mode, see radix16_ditN_cy_dif1.c.
				// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.

				tmp = half_arr+2;
				VEC_DBL_INIT(tmp, scale);
				/* Get the needed Nth root of -1: */
				add1 = (double *)&rn0[0];
				add2 = (double *)&rn1[0];

				idx_offset = j;
				idx_incr = NDIVR;

				tmp = base_negacyclic_root;	tm2 = tmp+1;

				// Hi-accuracy version needs 8 copies of each base root:
				l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				VEC_DBL_INIT(tmp+ 32,wt_re);	VEC_DBL_INIT(tm2+ 32,wt_im);
				VEC_DBL_INIT(tmp+ 40,wt_re);	VEC_DBL_INIT(tm2+ 40,wt_im);
				VEC_DBL_INIT(tmp+ 48,wt_re);	VEC_DBL_INIT(tm2+ 48,wt_im);
				VEC_DBL_INIT(tmp+ 56,wt_re);	VEC_DBL_INIT(tm2+ 56,wt_im);
				VEC_DBL_INIT(tmp+ 64,wt_re);	VEC_DBL_INIT(tm2+ 64,wt_im);
				VEC_DBL_INIT(tmp+ 72,wt_re);	VEC_DBL_INIT(tm2+ 72,wt_im);
				VEC_DBL_INIT(tmp+ 80,wt_re);	VEC_DBL_INIT(tm2+ 80,wt_im);
				VEC_DBL_INIT(tmp+ 88,wt_re);	VEC_DBL_INIT(tm2+ 88,wt_im);
				VEC_DBL_INIT(tmp+ 96,wt_re);	VEC_DBL_INIT(tm2+ 96,wt_im);
				VEC_DBL_INIT(tmp+104,wt_re);	VEC_DBL_INIT(tm2+104,wt_im);
				VEC_DBL_INIT(tmp+112,wt_re);	VEC_DBL_INIT(tm2+112,wt_im);
				VEC_DBL_INIT(tmp+120,wt_re);	VEC_DBL_INIT(tm2+120,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				VEC_DBL_INIT(tmp+ 32,wt_re);	VEC_DBL_INIT(tm2+ 32,wt_im);
				VEC_DBL_INIT(tmp+ 40,wt_re);	VEC_DBL_INIT(tm2+ 40,wt_im);
				VEC_DBL_INIT(tmp+ 48,wt_re);	VEC_DBL_INIT(tm2+ 48,wt_im);
				VEC_DBL_INIT(tmp+ 56,wt_re);	VEC_DBL_INIT(tm2+ 56,wt_im);
				VEC_DBL_INIT(tmp+ 64,wt_re);	VEC_DBL_INIT(tm2+ 64,wt_im);
				VEC_DBL_INIT(tmp+ 72,wt_re);	VEC_DBL_INIT(tm2+ 72,wt_im);
				VEC_DBL_INIT(tmp+ 80,wt_re);	VEC_DBL_INIT(tm2+ 80,wt_im);
				VEC_DBL_INIT(tmp+ 88,wt_re);	VEC_DBL_INIT(tm2+ 88,wt_im);
				VEC_DBL_INIT(tmp+ 96,wt_re);	VEC_DBL_INIT(tm2+ 96,wt_im);
				VEC_DBL_INIT(tmp+104,wt_re);	VEC_DBL_INIT(tm2+104,wt_im);
				VEC_DBL_INIT(tmp+112,wt_re);	VEC_DBL_INIT(tm2+112,wt_im);
				VEC_DBL_INIT(tmp+120,wt_re);	VEC_DBL_INIT(tm2+120,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				VEC_DBL_INIT(tmp+ 32,wt_re);	VEC_DBL_INIT(tm2+ 32,wt_im);
				VEC_DBL_INIT(tmp+ 40,wt_re);	VEC_DBL_INIT(tm2+ 40,wt_im);
				VEC_DBL_INIT(tmp+ 48,wt_re);	VEC_DBL_INIT(tm2+ 48,wt_im);
				VEC_DBL_INIT(tmp+ 56,wt_re);	VEC_DBL_INIT(tm2+ 56,wt_im);
				VEC_DBL_INIT(tmp+ 64,wt_re);	VEC_DBL_INIT(tm2+ 64,wt_im);
				VEC_DBL_INIT(tmp+ 72,wt_re);	VEC_DBL_INIT(tm2+ 72,wt_im);
				VEC_DBL_INIT(tmp+ 80,wt_re);	VEC_DBL_INIT(tm2+ 80,wt_im);
				VEC_DBL_INIT(tmp+ 88,wt_re);	VEC_DBL_INIT(tm2+ 88,wt_im);
				VEC_DBL_INIT(tmp+ 96,wt_re);	VEC_DBL_INIT(tm2+ 96,wt_im);
				VEC_DBL_INIT(tmp+104,wt_re);	VEC_DBL_INIT(tm2+104,wt_im);
				VEC_DBL_INIT(tmp+112,wt_re);	VEC_DBL_INIT(tm2+112,wt_im);
				VEC_DBL_INIT(tmp+120,wt_re);	VEC_DBL_INIT(tm2+120,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				VEC_DBL_INIT(tmp+ 32,wt_re);	VEC_DBL_INIT(tm2+ 32,wt_im);
				VEC_DBL_INIT(tmp+ 40,wt_re);	VEC_DBL_INIT(tm2+ 40,wt_im);
				VEC_DBL_INIT(tmp+ 48,wt_re);	VEC_DBL_INIT(tm2+ 48,wt_im);
				VEC_DBL_INIT(tmp+ 56,wt_re);	VEC_DBL_INIT(tm2+ 56,wt_im);
				VEC_DBL_INIT(tmp+ 64,wt_re);	VEC_DBL_INIT(tm2+ 64,wt_im);
				VEC_DBL_INIT(tmp+ 72,wt_re);	VEC_DBL_INIT(tm2+ 72,wt_im);
				VEC_DBL_INIT(tmp+ 80,wt_re);	VEC_DBL_INIT(tm2+ 80,wt_im);
				VEC_DBL_INIT(tmp+ 88,wt_re);	VEC_DBL_INIT(tm2+ 88,wt_im);
				VEC_DBL_INIT(tmp+ 96,wt_re);	VEC_DBL_INIT(tm2+ 96,wt_im);
				VEC_DBL_INIT(tmp+104,wt_re);	VEC_DBL_INIT(tm2+104,wt_im);
				VEC_DBL_INIT(tmp+112,wt_re);	VEC_DBL_INIT(tm2+112,wt_im);
				VEC_DBL_INIT(tmp+120,wt_re);	VEC_DBL_INIT(tm2+120,wt_im);

				// AVX-custom 4-way carry macro - each contains 4 of the RADIX stride-n/RADIX-separated carries
				// (processed independently in parallel), and steps through sequential-data indicies j,j+2,j+4,j+6:

				// The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl) = 2^12
				// which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
				// precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
				// uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls:
				tmp = base_negacyclic_root+  0;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p00r,tmp,0x1000,cy_r00,cy_i00,half_arr,sign_mask);
				tmp = base_negacyclic_root+  8;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p04r,tmp,0x0f40,cy_r04,cy_i04,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 16;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p08r,tmp,0x0e80,cy_r08,cy_i08,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 24;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p0cr,tmp,0x0dc0,cy_r0C,cy_i0C,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 32;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p10r,tmp,0x0d00,cy_r10,cy_i10,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 40;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p14r,tmp,0x0c40,cy_r14,cy_i14,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 48;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p18r,tmp,0x0b80,cy_r18,cy_i18,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 56;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p1cr,tmp,0x0ac0,cy_r1C,cy_i1C,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 64;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p20r,tmp,0x0a00,cy_r20,cy_i20,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 72;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p24r,tmp,0x0940,cy_r24,cy_i24,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 80;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p28r,tmp,0x0880,cy_r28,cy_i28,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 88;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p2cr,tmp,0x07c0,cy_r2C,cy_i2C,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 96;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p30r,tmp,0x0700,cy_r30,cy_i30,half_arr,sign_mask);
				tmp = base_negacyclic_root+104;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p34r,tmp,0x0640,cy_r34,cy_i34,half_arr,sign_mask);
				tmp = base_negacyclic_root+112;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p38r,tmp,0x0580,cy_r38,cy_i38,half_arr,sign_mask);
				tmp = base_negacyclic_root+120;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p3cr,tmp,0x04c0,cy_r3C,cy_i3C,half_arr,sign_mask);

			  #else	// USE_SSE2

				// For a description of the data movement for Fermat-mod carries in SSE2 mode, see radix16_ditN_cy_dif1.c.
				// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.

				tmp = half_arr+2;
				VEC_DBL_INIT(tmp, scale);
				/* Get the needed Nth root of -1: */
				add1 = (double *)&rn0[0];
				add2 = (double *)&rn1[0];

				idx_offset = j;
				idx_incr = NDIVR;

				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p00r,cy_r00,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p02r,cy_r04,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p04r,cy_r08,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p06r,cy_r0C,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p08r,cy_r10,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p0ar,cy_r14,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p0cr,cy_r18,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p0er,cy_r1C,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p10r,cy_r20,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p12r,cy_r24,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p14r,cy_r28,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p16r,cy_r2C,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p18r,cy_r30,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p1ar,cy_r34,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p1cr,cy_r38,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p1er,cy_r3C,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p20r,cy_i00,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p22r,cy_i04,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p24r,cy_i08,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p26r,cy_i0C,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p28r,cy_i10,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p2ar,cy_i14,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p2cr,cy_i18,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p2er,cy_i1C,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p30r,cy_i20,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p32r,cy_i24,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p34r,cy_i28,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p36r,cy_i2C,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p38r,cy_i30,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p3ar,cy_i34,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p3cr,cy_i38,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p3er,cy_i3C,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);

			  #endif	// AVX/SSE2?

			#else	// Scalar-double mode:

				fermat_carry_norm_pow2_errcheck(ar_p00,ai_p00,cy_r00,cy_i00,0x00*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p01,ai_p01,cy_r01,cy_i01,0x01*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p02,ai_p02,cy_r02,cy_i02,0x02*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p03,ai_p03,cy_r03,cy_i03,0x03*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p04,ai_p04,cy_r04,cy_i04,0x04*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p05,ai_p05,cy_r05,cy_i05,0x05*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p06,ai_p06,cy_r06,cy_i06,0x06*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p07,ai_p07,cy_r07,cy_i07,0x07*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p08,ai_p08,cy_r08,cy_i08,0x08*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p09,ai_p09,cy_r09,cy_i09,0x09*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p0a,ai_p0a,cy_r0A,cy_i0A,0x0A*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p0b,ai_p0b,cy_r0B,cy_i0B,0x0B*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p0c,ai_p0c,cy_r0C,cy_i0C,0x0C*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p0d,ai_p0d,cy_r0D,cy_i0D,0x0D*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p0e,ai_p0e,cy_r0E,cy_i0E,0x0E*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p0f,ai_p0f,cy_r0F,cy_i0F,0x0F*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p10,ai_p10,cy_r10,cy_i10,0x10*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p11,ai_p11,cy_r11,cy_i11,0x11*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p12,ai_p12,cy_r12,cy_i12,0x12*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p13,ai_p13,cy_r13,cy_i13,0x13*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p14,ai_p14,cy_r14,cy_i14,0x14*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p15,ai_p15,cy_r15,cy_i15,0x15*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p16,ai_p16,cy_r16,cy_i16,0x16*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p17,ai_p17,cy_r17,cy_i17,0x17*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p18,ai_p18,cy_r18,cy_i18,0x18*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p19,ai_p19,cy_r19,cy_i19,0x19*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p1a,ai_p1a,cy_r1A,cy_i1A,0x1A*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p1b,ai_p1b,cy_r1B,cy_i1B,0x1B*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p1c,ai_p1c,cy_r1C,cy_i1C,0x1C*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p1d,ai_p1d,cy_r1D,cy_i1D,0x1D*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p1e,ai_p1e,cy_r1E,cy_i1E,0x1E*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p1f,ai_p1f,cy_r1F,cy_i1F,0x1F*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p20,ai_p20,cy_r20,cy_i20,0x20*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p21,ai_p21,cy_r21,cy_i21,0x21*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p22,ai_p22,cy_r22,cy_i22,0x22*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p23,ai_p23,cy_r23,cy_i23,0x23*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p24,ai_p24,cy_r24,cy_i24,0x24*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p25,ai_p25,cy_r25,cy_i25,0x25*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p26,ai_p26,cy_r26,cy_i26,0x26*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p27,ai_p27,cy_r27,cy_i27,0x27*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p28,ai_p28,cy_r28,cy_i28,0x28*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p29,ai_p29,cy_r29,cy_i29,0x29*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p2a,ai_p2a,cy_r2A,cy_i2A,0x2A*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p2b,ai_p2b,cy_r2B,cy_i2B,0x2B*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p2c,ai_p2c,cy_r2C,cy_i2C,0x2C*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p2d,ai_p2d,cy_r2D,cy_i2D,0x2D*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p2e,ai_p2e,cy_r2E,cy_i2E,0x2E*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p2f,ai_p2f,cy_r2F,cy_i2F,0x2F*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p30,ai_p30,cy_r30,cy_i30,0x30*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p31,ai_p31,cy_r31,cy_i31,0x31*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p32,ai_p32,cy_r32,cy_i32,0x32*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p33,ai_p33,cy_r33,cy_i33,0x33*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p34,ai_p34,cy_r34,cy_i34,0x34*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p35,ai_p35,cy_r35,cy_i35,0x35*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p36,ai_p36,cy_r36,cy_i36,0x36*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p37,ai_p37,cy_r37,cy_i37,0x37*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p38,ai_p38,cy_r38,cy_i38,0x38*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p39,ai_p39,cy_r39,cy_i39,0x39*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p3a,ai_p3a,cy_r3A,cy_i3A,0x3A*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p3b,ai_p3b,cy_r3B,cy_i3B,0x3B*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p3c,ai_p3c,cy_r3C,cy_i3C,0x3C*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p3d,ai_p3d,cy_r3D,cy_i3D,0x3D*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p3e,ai_p3e,cy_r3E,cy_i3E,0x3E*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p3f,ai_p3f,cy_r3F,cy_i3F,0x3F*NDIVR,NRTM1,NRT_BITS);

			#endif	/* #ifdef USE_SSE2 */

			}	/* if(MODULUS_TYPE == ...) */

		/*...The radix-64 DIF pass is here:	*/

		#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
		  #ifdef USE_AVX

		/* Gather the needed data (64 64-bit complex, i.e. 128 64-bit reals) and do 8 twiddleless length-8 subtransforms: */
			//...Block 0: jt = j1;	jp = j2;
// Relative to RADIX_08_DIF_OOP, this SIMD macro produces outputs in BR order [04261537], so swap r-pointer pairs 2/8,6/C to handle that:
			SSE2_RADIX8_DIF_0TWIDDLE(
				s1p00r,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00,
				r00,r08,r04,r0C,r02,r0A,r06,r0E, isrt2
			);
			//...Block 1: jt = j1 + p04;	jp = j2 + p04;
			SSE2_RADIX8_DIF_0TWIDDLE(
				s1p04r,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00,
				r10,r18,r14,r1C,r12,r1A,r16,r1E, isrt2
			);
			//...Block 2: jt = j1 + p02;	jp = j2 + p02;
			SSE2_RADIX8_DIF_0TWIDDLE(
				s1p02r,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00,
				r20,r28,r24,r2C,r22,r2A,r26,r2E, isrt2
			);
			//...Block 3: jt = j1 + p06;	jp = j2 + p06;
			SSE2_RADIX8_DIF_0TWIDDLE(
				s1p06r,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00,
				r30,r38,r34,r3C,r32,r3A,r36,r3E, isrt2
			);
			//...Block 4: jt = j1 + p01;	jp = j2 + p01;
			SSE2_RADIX8_DIF_0TWIDDLE(
				s1p01r,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00,
				r40,r48,r44,r4C,r42,r4A,r46,r4E, isrt2
			);
			//...Block 5: jt = j1 + p05;	jp = j2 + p05;
			SSE2_RADIX8_DIF_0TWIDDLE(
				s1p05r,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00,
				r50,r58,r54,r5C,r52,r5A,r56,r5E, isrt2
			);
			//...Block 6: jt = j1 + p03;	jp = j2 + p03;
			SSE2_RADIX8_DIF_0TWIDDLE(
				s1p03r,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00,
				r60,r68,r64,r6C,r62,r6A,r66,r6E, isrt2
			);
			//...Block 7: jt = j1 + p07;	jp = j2 + p07;
			SSE2_RADIX8_DIF_0TWIDDLE(
				s1p07r,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00,
				r70,r78,r74,r7C,r72,r7A,r76,r7E, isrt2
			);

		/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dif_pass1 for details: */

			/* Block 0: */
			add0 = &a[j1]      ; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			/* 0-index block has all-unity twiddles: Remember, the twiddleless DIF bit-reverses both its in-and-outputs,
			so swap index-offset pairs 1/4 and 3/6 in t*-inputs and a-outputs: */
			SSE2_RADIX8_DIF_0TWIDDLE(
				r00, 0x800,0x400,0xc00,0x200,0xa00,0x600,0xe00,
				add0,add1,add2,add3,add4,add5,add6,add7, isrt2
			);
			/* Block 4: */
			add0 = &a[j1] + p08; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				r08,r48,r28,r68,r18,r58,r38,r78,
				add0,add1,add2,add3,add4,add5,add6,add7,
				ss0,cc0, isrt2,isrt2, nisrt2,isrt2, cc4,ss4, nss4,cc4, ss4,cc4, ncc4,ss4
			);
			/* Block 2: */
			add0 = &a[j1] + p10; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				r04,r44,r24,r64,r14,r54,r34,r74,
				add0,add1,add2,add3,add4,add5,add6,add7,
				isrt2,isrt2, cc4,ss4, ss4,cc4, cc2,ss2, ss6,cc6, cc6,ss6, ss2,cc2
			);
			/* Block 6: */
			add0 = &a[j1] + p18; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				r0C,r4C,r2C,r6C,r1C,r5C,r3C,r7C,
				add0,add1,add2,add3,add4,add5,add6,add7,
				nisrt2,isrt2, ss4,cc4, ncc4,nss4, cc6,ss6, ncc2,ss2 ,nss2,cc2, nss6,ncc6
			);
			/* Block 1: */
			add0 = &a[j1] + p20; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				r02,r42,r22,r62,r12,r52,r32,r72,
				add0,add1,add2,add3,add4,add5,add6,add7,
				cc4,ss4, cc2,ss2, cc6,ss6, cc1,ss1, cc5,ss5, cc3,ss3, cc7,ss7
			);
			/* Block 5: */
			add0 = &a[j1] + p28; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				r0A,r4A,r2A,r6A,r1A,r5A,r3A,r7A,
				add0,add1,add2,add3,add4,add5,add6,add7,
				nss4,cc4, ss6,cc6, ncc2,ss2, cc5,ss5, ncc7,ss7, ss1,cc1, ncc3,nss3
			);
			/* Block 3: */
			add0 = &a[j1] + p30; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				r06,r46,r26,r66,r16,r56,r36,r76,
				add0,add1,add2,add3,add4,add5,add6,add7,
				ss4,cc4, cc6,ss6, nss2,cc2, cc3,ss3, ss1,cc1, ss7,cc7, nss5,cc5
			);
			/* Block 7: */
			add0 = &a[j1] + p38; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				r0E,r4E,r2E,r6E,r1E,r5E,r3E,r7E,
				add0,add1,add2,add3,add4,add5,add6,add7,
				ncc4,ss4, ss2,cc2, nss6,ncc6, cc7,ss7, ncc3,nss3, nss5,cc5, ss1,ncc1
			);

		  #else	// USE_SSE2

		/* Gather the needed data (64 64-bit complex, i.e. 128 64-bit reals) and do 8 twiddleless length-8 subtransforms: */
			//...Block 0: jt = j1;	jp = j2;
// Relative to RADIX_08_DIF_OOP, this SIMD macro produces outputs in BR order [04261537], so swap r-pointer pairs 2/8,6/C to handle that:
			SSE2_RADIX8_DIF_0TWIDDLE(
				s1p00r,0x100,0x200,0x300,0x400,0x500,0x600,0x700,
				r00,r08,r04,r0C,r02,r0A,r06,r0E, isrt2
			);
			//...Block 1: jt = j1 + p04;	jp = j2 + p04;
			SSE2_RADIX8_DIF_0TWIDDLE(
				s1p04r,0x100,0x200,0x300,0x400,0x500,0x600,0x700,
				r10,r18,r14,r1C,r12,r1A,r16,r1E, isrt2
			);
			//...Block 2: jt = j1 + p02;	jp = j2 + p02;
			SSE2_RADIX8_DIF_0TWIDDLE(
				s1p02r,0x100,0x200,0x300,0x400,0x500,0x600,0x700,
				r20,r28,r24,r2C,r22,r2A,r26,r2E, isrt2
			);
			//...Block 3: jt = j1 + p06;	jp = j2 + p06;
			SSE2_RADIX8_DIF_0TWIDDLE(
				s1p06r,0x100,0x200,0x300,0x400,0x500,0x600,0x700,
				r30,r38,r34,r3C,r32,r3A,r36,r3E, isrt2
			);
			//...Block 4: jt = j1 + p01;	jp = j2 + p01;
			SSE2_RADIX8_DIF_0TWIDDLE(
				s1p01r,0x100,0x200,0x300,0x400,0x500,0x600,0x700,
				r40,r48,r44,r4C,r42,r4A,r46,r4E, isrt2
			);
			//...Block 5: jt = j1 + p05;	jp = j2 + p05;
			SSE2_RADIX8_DIF_0TWIDDLE(
				s1p05r,0x100,0x200,0x300,0x400,0x500,0x600,0x700,
				r50,r58,r54,r5C,r52,r5A,r56,r5E, isrt2
			);
			//...Block 6: jt = j1 + p03;	jp = j2 + p03;
			SSE2_RADIX8_DIF_0TWIDDLE(
				s1p03r,0x100,0x200,0x300,0x400,0x500,0x600,0x700,
				r60,r68,r64,r6C,r62,r6A,r66,r6E, isrt2
			);
			//...Block 7: jt = j1 + p07;	jp = j2 + p07;
			SSE2_RADIX8_DIF_0TWIDDLE(
				s1p07r,0x100,0x200,0x300,0x400,0x500,0x600,0x700,
				r70,r78,r74,r7C,r72,r7A,r76,r7E, isrt2
			);

		/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dif_pass1 for details: */

			/* Block 0: */
			add0 = &a[j1]      ; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			/* 0-index block has all-unity twiddles: Remember, the twiddleless DIF bit-reverses both its in-and-outputs,
			so swap index-offset pairs 1/4 and 3/6 in t*-inputs and a-outputs: */
			SSE2_RADIX8_DIF_0TWIDDLE(
				r00, 0x400,0x200,0x600,0x100,0x500,0x300,0x700,
				add0,add1,add2,add3,add4,add5,add6,add7, isrt2
			);
			/* Block 4: */
			add0 = &a[j1] + p08; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				r08,r48,r28,r68,r18,r58,r38,r78,
				add0,add1,add2,add3,add4,add5,add6,add7,
				ss0,cc0, isrt2,isrt2, nisrt2,isrt2, cc4,ss4, nss4,cc4, ss4,cc4, ncc4,ss4
			);
			/* Block 2: */
			add0 = &a[j1] + p10; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				r04,r44,r24,r64,r14,r54,r34,r74,
				add0,add1,add2,add3,add4,add5,add6,add7,
				isrt2,isrt2, cc4,ss4, ss4,cc4, cc2,ss2, ss6,cc6, cc6,ss6, ss2,cc2
			);
			/* Block 6: */
			add0 = &a[j1] + p18; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				r0C,r4C,r2C,r6C,r1C,r5C,r3C,r7C,
				add0,add1,add2,add3,add4,add5,add6,add7,
				nisrt2,isrt2, ss4,cc4, ncc4,nss4, cc6,ss6, ncc2,ss2 ,nss2,cc2, nss6,ncc6
			);
			/* Block 1: */
			add0 = &a[j1] + p20; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				r02,r42,r22,r62,r12,r52,r32,r72,
				add0,add1,add2,add3,add4,add5,add6,add7,
				cc4,ss4, cc2,ss2, cc6,ss6, cc1,ss1, cc5,ss5, cc3,ss3, cc7,ss7
			);
			/* Block 5: */
			add0 = &a[j1] + p28; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				r0A,r4A,r2A,r6A,r1A,r5A,r3A,r7A,
				add0,add1,add2,add3,add4,add5,add6,add7,
				nss4,cc4, ss6,cc6, ncc2,ss2, cc5,ss5, ncc7,ss7, ss1,cc1, ncc3,nss3
			);
			/* Block 3: */
			add0 = &a[j1] + p30; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				r06,r46,r26,r66,r16,r56,r36,r76,
				add0,add1,add2,add3,add4,add5,add6,add7,
				ss4,cc4, cc6,ss6, nss2,cc2, cc3,ss3, ss1,cc1, ss7,cc7, nss5,cc5
			);
			/* Block 7: */
			add0 = &a[j1] + p38; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
			SSE2_RADIX8_DIF_TWIDDLE_OOP(
				r0E,r4E,r2E,r6E,r1E,r5E,r3E,r7E,
				add0,add1,add2,add3,add4,add5,add6,add7,
				ncc4,ss4, ss2,cc2, nss6,ncc6, cc7,ss7, ncc3,nss3, nss5,cc5, ss1,ncc1
			);

		  #endif	// AVX/SSE2?

		#else	// USE_SSE2 = False:

		/* Gather the needed data (64 64-bit complex, i.e. 128 64-bit reals) and do 8 twiddleless length-8 subtransforms: */
			//...Block 0: jt = j1;	jp = j2;
			RADIX_08_DIF_OOP(
				ar_p00,ai_p00,ar_p08,ai_p08,ar_p10,ai_p10,ar_p18,ai_p18,ar_p20,ai_p20,ar_p28,ai_p28,ar_p30,ai_p30,ar_p38,ai_p38,
				t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F
			);
			//...Block 1: jt = j1 + p04;	jp = j2 + p04;
			RADIX_08_DIF_OOP(
				ar_p04,ai_p04,ar_p0c,ai_p0c,ar_p14,ai_p14,ar_p1c,ai_p1c,ar_p24,ai_p24,ar_p2c,ai_p2c,ar_p34,ai_p34,ar_p3c,ai_p3c,
				t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F
			);
			//...Block 2: jt = j1 + p02;	jp = j2 + p02;
			RADIX_08_DIF_OOP(
				ar_p02,ai_p02,ar_p0a,ai_p0a,ar_p12,ai_p12,ar_p1a,ai_p1a,ar_p22,ai_p22,ar_p2a,ai_p2a,ar_p32,ai_p32,ar_p3a,ai_p3a,
				t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F
			);
			//...Block 3: jt = j1 + p06;	jp = j2 + p06;
			RADIX_08_DIF_OOP(
				ar_p06,ai_p06,ar_p0e,ai_p0e,ar_p16,ai_p16,ar_p1e,ai_p1e,ar_p26,ai_p26,ar_p2e,ai_p2e,ar_p36,ai_p36,ar_p3e,ai_p3e,
				t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F
			);
			//...Block 4: jt = j1 + p01;	jp = j2 + p01;
			RADIX_08_DIF_OOP(
				ar_p01,ai_p01,ar_p09,ai_p09,ar_p11,ai_p11,ar_p19,ai_p19,ar_p21,ai_p21,ar_p29,ai_p29,ar_p31,ai_p31,ar_p39,ai_p39,
				t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t4A,t4B,t4C,t4D,t4E,t4F
			);
			//...Block 5: jt = j1 + p05;	jp = j2 + p05;
			RADIX_08_DIF_OOP(
				ar_p05,ai_p05,ar_p0d,ai_p0d,ar_p15,ai_p15,ar_p1d,ai_p1d,ar_p25,ai_p25,ar_p2d,ai_p2d,ar_p35,ai_p35,ar_p3d,ai_p3d,
				t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t5A,t5B,t5C,t5D,t5E,t5F
			);
			//...Block 6: jt = j1 + p03;	jp = j2 + p03;
			RADIX_08_DIF_OOP(
				ar_p03,ai_p03,ar_p0b,ai_p0b,ar_p13,ai_p13,ar_p1b,ai_p1b,ar_p23,ai_p23,ar_p2b,ai_p2b,ar_p33,ai_p33,ar_p3b,ai_p3b,
				t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t6A,t6B,t6C,t6D,t6E,t6F
			);
			//...Block 7: jt = j1 + p07;	jp = j2 + p07;
			RADIX_08_DIF_OOP(
				ar_p07,ai_p07,ar_p0f,ai_p0f,ar_p17,ai_p17,ar_p1f,ai_p1f,ar_p27,ai_p27,ar_p2f,ai_p2f,ar_p37,ai_p37,ar_p3f,ai_p3f,
				t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t7A,t7B,t7C,t7D,t7E,t7F
			);

		/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dif_pass1 for details: */

			/* Block 0: */
			jt = j1;	jp = j2;
			/* 0-index block has all-unity twiddles: Remember, the twiddleless DIF bit-reverses both its in-and-outputs,
			so swap index-offset pairs 1/4 and 3/6 in t*-inputs and a-outputs: */
			RADIX_08_DIF_OOP(
				t00,t01,t40,t41,t20,t21,t60,t61,t10,t11,t50,t51,t30,t31,t70,t71,
				a[jt],a[jp],a[jt+p04],a[jp+p04],a[jt+p02],a[jp+p02],a[jt+p06],a[jp+p06],a[jt+p01],a[jp+p01],a[jt+p05],a[jp+p05],a[jt+p03],a[jp+p03],a[jt+p07],a[jp+p07]
			);
			/* Block 4: */
			jt = j1 + p08;	jp = j2 + p08;
			RADIX_08_DIF_TWIDDLE_OOP(
				t08,t09,t18,t19,t28,t29,t38,t39,t48,t49,t58,t59,t68,t69,t78,t79,
				a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
				0,1,ISRT2,ISRT2,-ISRT2,ISRT2,c64_4,s64_4,-s64_4,c64_4,s64_4,c64_4,-c64_4,s64_4
			);
			/* Block 2: */
			jt = j1 + p10;	jp = j2 + p10;
			RADIX_08_DIF_TWIDDLE_OOP(
				t04,t05,t14,t15,t24,t25,t34,t35,t44,t45,t54,t55,t64,t65,t74,t75,
				a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
				ISRT2,ISRT2,c64_4,s64_4,s64_4,c64_4,c64_2,s64_2,s64_6,c64_6,c64_6,s64_6,s64_2,c64_2
			);
			/* Block 6: */
			jt = j1 + p18;	jp = j2 + p18;
			RADIX_08_DIF_TWIDDLE_OOP(
				t0C,t0D,t1C,t1D,t2C,t2D,t3C,t3D,t4C,t4D,t5C,t5D,t6C,t6D,t7C,t7D,
				a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
				-ISRT2,ISRT2,s64_4,c64_4,-c64_4,-s64_4,c64_6,s64_6,-c64_2,s64_2,-s64_2,c64_2,-s64_6,-c64_6
			);
			/* Block 1: */
			jt = j1 + p20;	jp = j2 + p20;
			RADIX_08_DIF_TWIDDLE_OOP(
				t02,t03,t12,t13,t22,t23,t32,t33,t42,t43,t52,t53,t62,t63,t72,t73,
				a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
				c64_4,s64_4,c64_2,s64_2,c64_6,s64_6,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7
			);
			/* Block 5: */
			jt = j1 + p28;	jp = j2 + p28;
			RADIX_08_DIF_TWIDDLE_OOP(
				t0A,t0B,t1A,t1B,t2A,t2B,t3A,t3B,t4A,t4B,t5A,t5B,t6A,t6B,t7A,t7B,
				a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
				-s64_4,c64_4,s64_6,c64_6,-c64_2,s64_2,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3
			);
			/* Block 3: */
			jt = j1 + p30;	jp = j2 + p30;
			RADIX_08_DIF_TWIDDLE_OOP(
				t06,t07,t16,t17,t26,t27,t36,t37,t46,t47,t56,t57,t66,t67,t76,t77,
				a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
				s64_4,c64_4,c64_6,s64_6,-s64_2,c64_2,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5
			);
			/* Block 7: */
			jt = j1 + p38;	jp = j2 + p38;
			RADIX_08_DIF_TWIDDLE_OOP(
				t0E,t0F,t1E,t1F,t2E,t2F,t3E,t3F,t4E,t4F,t5E,t5F,t6E,t6F,t7E,t7F,
				a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
				-c64_4,s64_4,s64_2,c64_2,-s64_6,-c64_6,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1
			);

		#endif	/* #ifdef USE_SSE2 */

			}	/* end for(j=_jstart[ithread]; j < _jhi[ithread]; j += 2) */

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
				jstart += nwt;
				jhi    += nwt;

				col += RADIX;
				co3 -= RADIX;
			}
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
		  #ifdef USE_AVX
			_cy_r00[ithread] = cy_r00->d0;	_cy_r01[ithread] = cy_r00->d1;	_cy_r02[ithread] = cy_r00->d2;	_cy_r03[ithread] = cy_r00->d3;
			_cy_r04[ithread] = cy_r04->d0;	_cy_r05[ithread] = cy_r04->d1;	_cy_r06[ithread] = cy_r04->d2;	_cy_r07[ithread] = cy_r04->d3;
			_cy_r08[ithread] = cy_r08->d0;	_cy_r09[ithread] = cy_r08->d1;	_cy_r0A[ithread] = cy_r08->d2;	_cy_r0B[ithread] = cy_r08->d3;
			_cy_r0C[ithread] = cy_r0C->d0;	_cy_r0D[ithread] = cy_r0C->d1;	_cy_r0E[ithread] = cy_r0C->d2;	_cy_r0F[ithread] = cy_r0C->d3;
			_cy_r10[ithread] = cy_r10->d0;	_cy_r11[ithread] = cy_r10->d1;	_cy_r12[ithread] = cy_r10->d2;	_cy_r13[ithread] = cy_r10->d3;
			_cy_r14[ithread] = cy_r14->d0;	_cy_r15[ithread] = cy_r14->d1;	_cy_r16[ithread] = cy_r14->d2;	_cy_r17[ithread] = cy_r14->d3;
			_cy_r18[ithread] = cy_r18->d0;	_cy_r19[ithread] = cy_r18->d1;	_cy_r1A[ithread] = cy_r18->d2;	_cy_r1B[ithread] = cy_r18->d3;
			_cy_r1C[ithread] = cy_r1C->d0;	_cy_r1D[ithread] = cy_r1C->d1;	_cy_r1E[ithread] = cy_r1C->d2;	_cy_r1F[ithread] = cy_r1C->d3;
			_cy_r20[ithread] = cy_r20->d0;	_cy_r21[ithread] = cy_r20->d1;	_cy_r22[ithread] = cy_r20->d2;	_cy_r23[ithread] = cy_r20->d3;
			_cy_r24[ithread] = cy_r24->d0;	_cy_r25[ithread] = cy_r24->d1;	_cy_r26[ithread] = cy_r24->d2;	_cy_r27[ithread] = cy_r24->d3;
			_cy_r28[ithread] = cy_r28->d0;	_cy_r29[ithread] = cy_r28->d1;	_cy_r2A[ithread] = cy_r28->d2;	_cy_r2B[ithread] = cy_r28->d3;
			_cy_r2C[ithread] = cy_r2C->d0;	_cy_r2D[ithread] = cy_r2C->d1;	_cy_r2E[ithread] = cy_r2C->d2;	_cy_r2F[ithread] = cy_r2C->d3;
			_cy_r30[ithread] = cy_r30->d0;	_cy_r31[ithread] = cy_r30->d1;	_cy_r32[ithread] = cy_r30->d2;	_cy_r33[ithread] = cy_r30->d3;
			_cy_r34[ithread] = cy_r34->d0;	_cy_r35[ithread] = cy_r34->d1;	_cy_r36[ithread] = cy_r34->d2;	_cy_r37[ithread] = cy_r34->d3;
			_cy_r38[ithread] = cy_r38->d0;	_cy_r39[ithread] = cy_r38->d1;	_cy_r3A[ithread] = cy_r38->d2;	_cy_r3B[ithread] = cy_r38->d3;
			_cy_r3C[ithread] = cy_r3C->d0;	_cy_r3D[ithread] = cy_r3C->d1;	_cy_r3E[ithread] = cy_r3C->d2;	_cy_r3F[ithread] = cy_r3C->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		  #else	// USE_SSE2
			_cy_r00[ithread] = cy_r00->d0;	_cy_r01[ithread] = cy_r00->d1;
			_cy_r02[ithread] = cy_r02->d0;	_cy_r03[ithread] = cy_r02->d1;
			_cy_r04[ithread] = cy_r04->d0;	_cy_r05[ithread] = cy_r04->d1;
			_cy_r06[ithread] = cy_r06->d0;	_cy_r07[ithread] = cy_r06->d1;
			_cy_r08[ithread] = cy_r08->d0;	_cy_r09[ithread] = cy_r08->d1;
			_cy_r0A[ithread] = cy_r0A->d0;	_cy_r0B[ithread] = cy_r0A->d1;
			_cy_r0C[ithread] = cy_r0C->d0;	_cy_r0D[ithread] = cy_r0C->d1;
			_cy_r0E[ithread] = cy_r0E->d0;	_cy_r0F[ithread] = cy_r0E->d1;
			_cy_r10[ithread] = cy_r10->d0;	_cy_r11[ithread] = cy_r10->d1;
			_cy_r12[ithread] = cy_r12->d0;	_cy_r13[ithread] = cy_r12->d1;
			_cy_r14[ithread] = cy_r14->d0;	_cy_r15[ithread] = cy_r14->d1;
			_cy_r16[ithread] = cy_r16->d0;	_cy_r17[ithread] = cy_r16->d1;
			_cy_r18[ithread] = cy_r18->d0;	_cy_r19[ithread] = cy_r18->d1;
			_cy_r1A[ithread] = cy_r1A->d0;	_cy_r1B[ithread] = cy_r1A->d1;
			_cy_r1C[ithread] = cy_r1C->d0;	_cy_r1D[ithread] = cy_r1C->d1;
			_cy_r1E[ithread] = cy_r1E->d0;	_cy_r1F[ithread] = cy_r1E->d1;
			_cy_r20[ithread] = cy_r20->d0;	_cy_r21[ithread] = cy_r20->d1;
			_cy_r22[ithread] = cy_r22->d0;	_cy_r23[ithread] = cy_r22->d1;
			_cy_r24[ithread] = cy_r24->d0;	_cy_r25[ithread] = cy_r24->d1;
			_cy_r26[ithread] = cy_r26->d0;	_cy_r27[ithread] = cy_r26->d1;
			_cy_r28[ithread] = cy_r28->d0;	_cy_r29[ithread] = cy_r28->d1;
			_cy_r2A[ithread] = cy_r2A->d0;	_cy_r2B[ithread] = cy_r2A->d1;
			_cy_r2C[ithread] = cy_r2C->d0;	_cy_r2D[ithread] = cy_r2C->d1;
			_cy_r2E[ithread] = cy_r2E->d0;	_cy_r2F[ithread] = cy_r2E->d1;
			_cy_r30[ithread] = cy_r30->d0;	_cy_r31[ithread] = cy_r30->d1;
			_cy_r32[ithread] = cy_r32->d0;	_cy_r33[ithread] = cy_r32->d1;
			_cy_r34[ithread] = cy_r34->d0;	_cy_r35[ithread] = cy_r34->d1;
			_cy_r36[ithread] = cy_r36->d0;	_cy_r37[ithread] = cy_r36->d1;
			_cy_r38[ithread] = cy_r38->d0;	_cy_r39[ithread] = cy_r38->d1;
			_cy_r3A[ithread] = cy_r3A->d0;	_cy_r3B[ithread] = cy_r3A->d1;
			_cy_r3C[ithread] = cy_r3C->d0;	_cy_r3D[ithread] = cy_r3C->d1;
			_cy_r3E[ithread] = cy_r3E->d0;	_cy_r3F[ithread] = cy_r3E->d1;
			maxerr = MAX(max_err->d0,max_err->d1);
		  #endif	// AVX/SSE2?
		#else
			_cy_r00[ithread] = cy_r00;		_cy_r20[ithread] = cy_r20;
			_cy_r01[ithread] = cy_r01;		_cy_r21[ithread] = cy_r21;
			_cy_r02[ithread] = cy_r02;		_cy_r22[ithread] = cy_r22;
			_cy_r03[ithread] = cy_r03;		_cy_r23[ithread] = cy_r23;
			_cy_r04[ithread] = cy_r04;		_cy_r24[ithread] = cy_r24;
			_cy_r05[ithread] = cy_r05;		_cy_r25[ithread] = cy_r25;
			_cy_r06[ithread] = cy_r06;		_cy_r26[ithread] = cy_r26;
			_cy_r07[ithread] = cy_r07;		_cy_r27[ithread] = cy_r27;
			_cy_r08[ithread] = cy_r08;		_cy_r28[ithread] = cy_r28;
			_cy_r09[ithread] = cy_r09;		_cy_r29[ithread] = cy_r29;
			_cy_r0A[ithread] = cy_r0A;		_cy_r2A[ithread] = cy_r2A;
			_cy_r0B[ithread] = cy_r0B;		_cy_r2B[ithread] = cy_r2B;
			_cy_r0C[ithread] = cy_r0C;		_cy_r2C[ithread] = cy_r2C;
			_cy_r0D[ithread] = cy_r0D;		_cy_r2D[ithread] = cy_r2D;
			_cy_r0E[ithread] = cy_r0E;		_cy_r2E[ithread] = cy_r2E;
			_cy_r0F[ithread] = cy_r0F;		_cy_r2F[ithread] = cy_r2F;
			_cy_r10[ithread] = cy_r10;		_cy_r30[ithread] = cy_r30;
			_cy_r11[ithread] = cy_r11;		_cy_r31[ithread] = cy_r31;
			_cy_r12[ithread] = cy_r12;		_cy_r32[ithread] = cy_r32;
			_cy_r13[ithread] = cy_r13;		_cy_r33[ithread] = cy_r33;
			_cy_r14[ithread] = cy_r14;		_cy_r34[ithread] = cy_r34;
			_cy_r15[ithread] = cy_r15;		_cy_r35[ithread] = cy_r35;
			_cy_r16[ithread] = cy_r16;		_cy_r36[ithread] = cy_r36;
			_cy_r17[ithread] = cy_r17;		_cy_r37[ithread] = cy_r37;
			_cy_r18[ithread] = cy_r18;		_cy_r38[ithread] = cy_r38;
			_cy_r19[ithread] = cy_r19;		_cy_r39[ithread] = cy_r39;
			_cy_r1A[ithread] = cy_r1A;		_cy_r3A[ithread] = cy_r3A;
			_cy_r1B[ithread] = cy_r1B;		_cy_r3B[ithread] = cy_r3B;
			_cy_r1C[ithread] = cy_r1C;		_cy_r3C[ithread] = cy_r3C;
			_cy_r1D[ithread] = cy_r1D;		_cy_r3D[ithread] = cy_r3D;
			_cy_r1E[ithread] = cy_r1E;		_cy_r3E[ithread] = cy_r3E;
			_cy_r1F[ithread] = cy_r1F;		_cy_r3F[ithread] = cy_r3F;
		#endif
		}
		else
		{
		#if defined(USE_SSE2) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
		  #ifdef USE_AVX
			_cy_r00[ithread] = cy_r00->d0;	_cy_r01[ithread] = cy_r00->d1;	_cy_r02[ithread] = cy_r00->d2;	_cy_r03[ithread] = cy_r00->d3;
			_cy_r04[ithread] = cy_r04->d0;	_cy_r05[ithread] = cy_r04->d1;	_cy_r06[ithread] = cy_r04->d2;	_cy_r07[ithread] = cy_r04->d3;
			_cy_r08[ithread] = cy_r08->d0;	_cy_r09[ithread] = cy_r08->d1;	_cy_r0A[ithread] = cy_r08->d2;	_cy_r0B[ithread] = cy_r08->d3;
			_cy_r0C[ithread] = cy_r0C->d0;	_cy_r0D[ithread] = cy_r0C->d1;	_cy_r0E[ithread] = cy_r0C->d2;	_cy_r0F[ithread] = cy_r0C->d3;
			_cy_r10[ithread] = cy_r10->d0;	_cy_r11[ithread] = cy_r10->d1;	_cy_r12[ithread] = cy_r10->d2;	_cy_r13[ithread] = cy_r10->d3;
			_cy_r14[ithread] = cy_r14->d0;	_cy_r15[ithread] = cy_r14->d1;	_cy_r16[ithread] = cy_r14->d2;	_cy_r17[ithread] = cy_r14->d3;
			_cy_r18[ithread] = cy_r18->d0;	_cy_r19[ithread] = cy_r18->d1;	_cy_r1A[ithread] = cy_r18->d2;	_cy_r1B[ithread] = cy_r18->d3;
			_cy_r1C[ithread] = cy_r1C->d0;	_cy_r1D[ithread] = cy_r1C->d1;	_cy_r1E[ithread] = cy_r1C->d2;	_cy_r1F[ithread] = cy_r1C->d3;
			_cy_r20[ithread] = cy_r20->d0;	_cy_r21[ithread] = cy_r20->d1;	_cy_r22[ithread] = cy_r20->d2;	_cy_r23[ithread] = cy_r20->d3;
			_cy_r24[ithread] = cy_r24->d0;	_cy_r25[ithread] = cy_r24->d1;	_cy_r26[ithread] = cy_r24->d2;	_cy_r27[ithread] = cy_r24->d3;
			_cy_r28[ithread] = cy_r28->d0;	_cy_r29[ithread] = cy_r28->d1;	_cy_r2A[ithread] = cy_r28->d2;	_cy_r2B[ithread] = cy_r28->d3;
			_cy_r2C[ithread] = cy_r2C->d0;	_cy_r2D[ithread] = cy_r2C->d1;	_cy_r2E[ithread] = cy_r2C->d2;	_cy_r2F[ithread] = cy_r2C->d3;
			_cy_r30[ithread] = cy_r30->d0;	_cy_r31[ithread] = cy_r30->d1;	_cy_r32[ithread] = cy_r30->d2;	_cy_r33[ithread] = cy_r30->d3;
			_cy_r34[ithread] = cy_r34->d0;	_cy_r35[ithread] = cy_r34->d1;	_cy_r36[ithread] = cy_r34->d2;	_cy_r37[ithread] = cy_r34->d3;
			_cy_r38[ithread] = cy_r38->d0;	_cy_r39[ithread] = cy_r38->d1;	_cy_r3A[ithread] = cy_r38->d2;	_cy_r3B[ithread] = cy_r38->d3;
			_cy_r3C[ithread] = cy_r3C->d0;	_cy_r3D[ithread] = cy_r3C->d1;	_cy_r3E[ithread] = cy_r3C->d2;	_cy_r3F[ithread] = cy_r3C->d3;

			_cy_i00[ithread] = cy_i00->d0;	_cy_i01[ithread] = cy_i00->d1;	_cy_i02[ithread] = cy_i00->d2;	_cy_i03[ithread] = cy_i00->d3;
			_cy_i04[ithread] = cy_i04->d0;	_cy_i05[ithread] = cy_i04->d1;	_cy_i06[ithread] = cy_i04->d2;	_cy_i07[ithread] = cy_i04->d3;
			_cy_i08[ithread] = cy_i08->d0;	_cy_i09[ithread] = cy_i08->d1;	_cy_i0A[ithread] = cy_i08->d2;	_cy_i0B[ithread] = cy_i08->d3;
			_cy_i0C[ithread] = cy_i0C->d0;	_cy_i0D[ithread] = cy_i0C->d1;	_cy_i0E[ithread] = cy_i0C->d2;	_cy_i0F[ithread] = cy_i0C->d3;
			_cy_i10[ithread] = cy_i10->d0;	_cy_i11[ithread] = cy_i10->d1;	_cy_i12[ithread] = cy_i10->d2;	_cy_i13[ithread] = cy_i10->d3;
			_cy_i14[ithread] = cy_i14->d0;	_cy_i15[ithread] = cy_i14->d1;	_cy_i16[ithread] = cy_i14->d2;	_cy_i17[ithread] = cy_i14->d3;
			_cy_i18[ithread] = cy_i18->d0;	_cy_i19[ithread] = cy_i18->d1;	_cy_i1A[ithread] = cy_i18->d2;	_cy_i1B[ithread] = cy_i18->d3;
			_cy_i1C[ithread] = cy_i1C->d0;	_cy_i1D[ithread] = cy_i1C->d1;	_cy_i1E[ithread] = cy_i1C->d2;	_cy_i1F[ithread] = cy_i1C->d3;
			_cy_i20[ithread] = cy_i20->d0;	_cy_i21[ithread] = cy_i20->d1;	_cy_i22[ithread] = cy_i20->d2;	_cy_i23[ithread] = cy_i20->d3;
			_cy_i24[ithread] = cy_i24->d0;	_cy_i25[ithread] = cy_i24->d1;	_cy_i26[ithread] = cy_i24->d2;	_cy_i27[ithread] = cy_i24->d3;
			_cy_i28[ithread] = cy_i28->d0;	_cy_i29[ithread] = cy_i28->d1;	_cy_i2A[ithread] = cy_i28->d2;	_cy_i2B[ithread] = cy_i28->d3;
			_cy_i2C[ithread] = cy_i2C->d0;	_cy_i2D[ithread] = cy_i2C->d1;	_cy_i2E[ithread] = cy_i2C->d2;	_cy_i2F[ithread] = cy_i2C->d3;
			_cy_i30[ithread] = cy_i30->d0;	_cy_i31[ithread] = cy_i30->d1;	_cy_i32[ithread] = cy_i30->d2;	_cy_i33[ithread] = cy_i30->d3;
			_cy_i34[ithread] = cy_i34->d0;	_cy_i35[ithread] = cy_i34->d1;	_cy_i36[ithread] = cy_i34->d2;	_cy_i37[ithread] = cy_i34->d3;
			_cy_i38[ithread] = cy_i38->d0;	_cy_i39[ithread] = cy_i38->d1;	_cy_i3A[ithread] = cy_i38->d2;	_cy_i3B[ithread] = cy_i38->d3;
			_cy_i3C[ithread] = cy_i3C->d0;	_cy_i3D[ithread] = cy_i3C->d1;	_cy_i3E[ithread] = cy_i3C->d2;	_cy_i3F[ithread] = cy_i3C->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		  #else // USE_SSE2
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			_cy_r00[ithread] = cy_r00->d0;	_cy_i00[ithread] = cy_r00->d1;
			_cy_r01[ithread] = cy_r02->d0;	_cy_i01[ithread] = cy_r02->d1;
			_cy_r02[ithread] = cy_r04->d0;	_cy_i02[ithread] = cy_r04->d1;
			_cy_r03[ithread] = cy_r06->d0;	_cy_i03[ithread] = cy_r06->d1;
			_cy_r04[ithread] = cy_r08->d0;	_cy_i04[ithread] = cy_r08->d1;
			_cy_r05[ithread] = cy_r0A->d0;	_cy_i05[ithread] = cy_r0A->d1;
			_cy_r06[ithread] = cy_r0C->d0;	_cy_i06[ithread] = cy_r0C->d1;
			_cy_r07[ithread] = cy_r0E->d0;	_cy_i07[ithread] = cy_r0E->d1;
			_cy_r08[ithread] = cy_r10->d0;	_cy_i08[ithread] = cy_r10->d1;
			_cy_r09[ithread] = cy_r12->d0;	_cy_i09[ithread] = cy_r12->d1;
			_cy_r0A[ithread] = cy_r14->d0;	_cy_i0A[ithread] = cy_r14->d1;
			_cy_r0B[ithread] = cy_r16->d0;	_cy_i0B[ithread] = cy_r16->d1;
			_cy_r0C[ithread] = cy_r18->d0;	_cy_i0C[ithread] = cy_r18->d1;
			_cy_r0D[ithread] = cy_r1A->d0;	_cy_i0D[ithread] = cy_r1A->d1;
			_cy_r0E[ithread] = cy_r1C->d0;	_cy_i0E[ithread] = cy_r1C->d1;
			_cy_r0F[ithread] = cy_r1E->d0;	_cy_i0F[ithread] = cy_r1E->d1;
			_cy_r10[ithread] = cy_r20->d0;	_cy_i10[ithread] = cy_r20->d1;
			_cy_r11[ithread] = cy_r22->d0;	_cy_i11[ithread] = cy_r22->d1;
			_cy_r12[ithread] = cy_r24->d0;	_cy_i12[ithread] = cy_r24->d1;
			_cy_r13[ithread] = cy_r26->d0;	_cy_i13[ithread] = cy_r26->d1;
			_cy_r14[ithread] = cy_r28->d0;	_cy_i14[ithread] = cy_r28->d1;
			_cy_r15[ithread] = cy_r2A->d0;	_cy_i15[ithread] = cy_r2A->d1;
			_cy_r16[ithread] = cy_r2C->d0;	_cy_i16[ithread] = cy_r2C->d1;
			_cy_r17[ithread] = cy_r2E->d0;	_cy_i17[ithread] = cy_r2E->d1;
			_cy_r18[ithread] = cy_r30->d0;	_cy_i18[ithread] = cy_r30->d1;
			_cy_r19[ithread] = cy_r32->d0;	_cy_i19[ithread] = cy_r32->d1;
			_cy_r1A[ithread] = cy_r34->d0;	_cy_i1A[ithread] = cy_r34->d1;
			_cy_r1B[ithread] = cy_r36->d0;	_cy_i1B[ithread] = cy_r36->d1;
			_cy_r1C[ithread] = cy_r38->d0;	_cy_i1C[ithread] = cy_r38->d1;
			_cy_r1D[ithread] = cy_r3A->d0;	_cy_i1D[ithread] = cy_r3A->d1;
			_cy_r1E[ithread] = cy_r3C->d0;	_cy_i1E[ithread] = cy_r3C->d1;
			_cy_r1F[ithread] = cy_r3E->d0;	_cy_i1F[ithread] = cy_r3E->d1;

			_cy_r20[ithread] = cy_i00->d0;	_cy_i20[ithread] = cy_i00->d1;
			_cy_r21[ithread] = cy_i02->d0;	_cy_i21[ithread] = cy_i02->d1;
			_cy_r22[ithread] = cy_i04->d0;	_cy_i22[ithread] = cy_i04->d1;
			_cy_r23[ithread] = cy_i06->d0;	_cy_i23[ithread] = cy_i06->d1;
			_cy_r24[ithread] = cy_i08->d0;	_cy_i24[ithread] = cy_i08->d1;
			_cy_r25[ithread] = cy_i0A->d0;	_cy_i25[ithread] = cy_i0A->d1;
			_cy_r26[ithread] = cy_i0C->d0;	_cy_i26[ithread] = cy_i0C->d1;
			_cy_r27[ithread] = cy_i0E->d0;	_cy_i27[ithread] = cy_i0E->d1;
			_cy_r28[ithread] = cy_i10->d0;	_cy_i28[ithread] = cy_i10->d1;
			_cy_r29[ithread] = cy_i12->d0;	_cy_i29[ithread] = cy_i12->d1;
			_cy_r2A[ithread] = cy_i14->d0;	_cy_i2A[ithread] = cy_i14->d1;
			_cy_r2B[ithread] = cy_i16->d0;	_cy_i2B[ithread] = cy_i16->d1;
			_cy_r2C[ithread] = cy_i18->d0;	_cy_i2C[ithread] = cy_i18->d1;
			_cy_r2D[ithread] = cy_i1A->d0;	_cy_i2D[ithread] = cy_i1A->d1;
			_cy_r2E[ithread] = cy_i1C->d0;	_cy_i2E[ithread] = cy_i1C->d1;
			_cy_r2F[ithread] = cy_i1E->d0;	_cy_i2F[ithread] = cy_i1E->d1;
			_cy_r30[ithread] = cy_i20->d0;	_cy_i30[ithread] = cy_i20->d1;
			_cy_r31[ithread] = cy_i22->d0;	_cy_i31[ithread] = cy_i22->d1;
			_cy_r32[ithread] = cy_i24->d0;	_cy_i32[ithread] = cy_i24->d1;
			_cy_r33[ithread] = cy_i26->d0;	_cy_i33[ithread] = cy_i26->d1;
			_cy_r34[ithread] = cy_i28->d0;	_cy_i34[ithread] = cy_i28->d1;
			_cy_r35[ithread] = cy_i2A->d0;	_cy_i35[ithread] = cy_i2A->d1;
			_cy_r36[ithread] = cy_i2C->d0;	_cy_i36[ithread] = cy_i2C->d1;
			_cy_r37[ithread] = cy_i2E->d0;	_cy_i37[ithread] = cy_i2E->d1;
			_cy_r38[ithread] = cy_i30->d0;	_cy_i38[ithread] = cy_i30->d1;
			_cy_r39[ithread] = cy_i32->d0;	_cy_i39[ithread] = cy_i32->d1;
			_cy_r3A[ithread] = cy_i34->d0;	_cy_i3A[ithread] = cy_i34->d1;
			_cy_r3B[ithread] = cy_i36->d0;	_cy_i3B[ithread] = cy_i36->d1;
			_cy_r3C[ithread] = cy_i38->d0;	_cy_i3C[ithread] = cy_i38->d1;
			_cy_r3D[ithread] = cy_i3A->d0;	_cy_i3D[ithread] = cy_i3A->d1;
			_cy_r3E[ithread] = cy_i3C->d0;	_cy_i3E[ithread] = cy_i3C->d1;
			_cy_r3F[ithread] = cy_i3E->d0;	_cy_i3F[ithread] = cy_i3E->d1;
			maxerr = MAX(max_err->d0,max_err->d1);
		  #endif	// AVX/SSE2?
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
			_cy_r0A[ithread] = cy_r0A;	_cy_i0A[ithread] = cy_i0A;
			_cy_r0B[ithread] = cy_r0B;	_cy_i0B[ithread] = cy_i0B;
			_cy_r0C[ithread] = cy_r0C;	_cy_i0C[ithread] = cy_i0C;
			_cy_r0D[ithread] = cy_r0D;	_cy_i0D[ithread] = cy_i0D;
			_cy_r0E[ithread] = cy_r0E;	_cy_i0E[ithread] = cy_i0E;
			_cy_r0F[ithread] = cy_r0F;	_cy_i0F[ithread] = cy_i0F;
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
			_cy_r1A[ithread] = cy_r1A;	_cy_i1A[ithread] = cy_i1A;
			_cy_r1B[ithread] = cy_r1B;	_cy_i1B[ithread] = cy_i1B;
			_cy_r1C[ithread] = cy_r1C;	_cy_i1C[ithread] = cy_i1C;
			_cy_r1D[ithread] = cy_r1D;	_cy_i1D[ithread] = cy_i1D;
			_cy_r1E[ithread] = cy_r1E;	_cy_i1E[ithread] = cy_i1E;
			_cy_r1F[ithread] = cy_r1F;	_cy_i1F[ithread] = cy_i1F;

			_cy_r20[ithread] = cy_r20;	_cy_i20[ithread] = cy_i20;
			_cy_r21[ithread] = cy_r21;	_cy_i21[ithread] = cy_i21;
			_cy_r22[ithread] = cy_r22;	_cy_i22[ithread] = cy_i22;
			_cy_r23[ithread] = cy_r23;	_cy_i23[ithread] = cy_i23;
			_cy_r24[ithread] = cy_r24;	_cy_i24[ithread] = cy_i24;
			_cy_r25[ithread] = cy_r25;	_cy_i25[ithread] = cy_i25;
			_cy_r26[ithread] = cy_r26;	_cy_i26[ithread] = cy_i26;
			_cy_r27[ithread] = cy_r27;	_cy_i27[ithread] = cy_i27;
			_cy_r28[ithread] = cy_r28;	_cy_i28[ithread] = cy_i28;
			_cy_r29[ithread] = cy_r29;	_cy_i29[ithread] = cy_i29;
			_cy_r2A[ithread] = cy_r2A;	_cy_i2A[ithread] = cy_i2A;
			_cy_r2B[ithread] = cy_r2B;	_cy_i2B[ithread] = cy_i2B;
			_cy_r2C[ithread] = cy_r2C;	_cy_i2C[ithread] = cy_i2C;
			_cy_r2D[ithread] = cy_r2D;	_cy_i2D[ithread] = cy_i2D;
			_cy_r2E[ithread] = cy_r2E;	_cy_i2E[ithread] = cy_i2E;
			_cy_r2F[ithread] = cy_r2F;	_cy_i2F[ithread] = cy_i2F;
			_cy_r30[ithread] = cy_r30;	_cy_i30[ithread] = cy_i30;
			_cy_r31[ithread] = cy_r31;	_cy_i31[ithread] = cy_i31;
			_cy_r32[ithread] = cy_r32;	_cy_i32[ithread] = cy_i32;
			_cy_r33[ithread] = cy_r33;	_cy_i33[ithread] = cy_i33;
			_cy_r34[ithread] = cy_r34;	_cy_i34[ithread] = cy_i34;
			_cy_r35[ithread] = cy_r35;	_cy_i35[ithread] = cy_i35;
			_cy_r36[ithread] = cy_r36;	_cy_i36[ithread] = cy_i36;
			_cy_r37[ithread] = cy_r37;	_cy_i37[ithread] = cy_i37;
			_cy_r38[ithread] = cy_r38;	_cy_i38[ithread] = cy_i38;
			_cy_r39[ithread] = cy_r39;	_cy_i39[ithread] = cy_i39;
			_cy_r3A[ithread] = cy_r3A;	_cy_i3A[ithread] = cy_i3A;
			_cy_r3B[ithread] = cy_r3B;	_cy_i3B[ithread] = cy_i3B;
			_cy_r3C[ithread] = cy_r3C;	_cy_i3C[ithread] = cy_i3C;
			_cy_r3D[ithread] = cy_r3D;	_cy_i3D[ithread] = cy_i3D;
			_cy_r3E[ithread] = cy_r3E;	_cy_i3E[ithread] = cy_i3E;
			_cy_r3F[ithread] = cy_r3F;	_cy_i3F[ithread] = cy_i3F;
		#endif
		}

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(_maxerr[ithread] < maxerr)
		{
			_maxerr[ithread] = maxerr;
		}

  #endif	// #ifdef USE_PTHREAD

	}	/******* END OF PARALLEL FOR-LOOP ********/

#ifdef USE_PTHREAD	// End of threadpool-based dispatch: Add a small wait-loop to ensure all threads complete

  #ifdef OS_TYPE_MACOSX

	/*** Main execution thread executes remaining chunks in serial fashion (but in || with the pool threads): ***/
	for(j = 0; j < main_work_units; ++j)
	{
	//	printf("adding main task %d\n",j + pool_work_units);
		ASSERT(HERE, 0x0 == cy64_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;
	ns_time.tv_sec  = 0.0001;// (time_t)seconds
	ns_time.tv_nsec = 0;	// (long)nanoseconds - At least allegedly, but under OS X it seems to be finer-grained than that

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(HERE, 0 == nanosleep(&ns_time, 0x0), "nanosleep fail!");
	}
//	printf("radix64_ditN_cy_dif1 end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = tdat[ithread].maxerr;
		if(maxerr < _maxerr[ithread]) {
			maxerr = _maxerr[ithread];
		}

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			_cy_r00[ithread] = tdat[ithread].cy_r00;	_cy_r20[ithread] = tdat[ithread].cy_r20;
			_cy_r01[ithread] = tdat[ithread].cy_r01;	_cy_r21[ithread] = tdat[ithread].cy_r21;
			_cy_r02[ithread] = tdat[ithread].cy_r02;	_cy_r22[ithread] = tdat[ithread].cy_r22;
			_cy_r03[ithread] = tdat[ithread].cy_r03;	_cy_r23[ithread] = tdat[ithread].cy_r23;
			_cy_r04[ithread] = tdat[ithread].cy_r04;	_cy_r24[ithread] = tdat[ithread].cy_r24;
			_cy_r05[ithread] = tdat[ithread].cy_r05;	_cy_r25[ithread] = tdat[ithread].cy_r25;
			_cy_r06[ithread] = tdat[ithread].cy_r06;	_cy_r26[ithread] = tdat[ithread].cy_r26;
			_cy_r07[ithread] = tdat[ithread].cy_r07;	_cy_r27[ithread] = tdat[ithread].cy_r27;
			_cy_r08[ithread] = tdat[ithread].cy_r08;	_cy_r28[ithread] = tdat[ithread].cy_r28;
			_cy_r09[ithread] = tdat[ithread].cy_r09;	_cy_r29[ithread] = tdat[ithread].cy_r29;
			_cy_r0A[ithread] = tdat[ithread].cy_r0A;	_cy_r2A[ithread] = tdat[ithread].cy_r2A;
			_cy_r0B[ithread] = tdat[ithread].cy_r0B;	_cy_r2B[ithread] = tdat[ithread].cy_r2B;
			_cy_r0C[ithread] = tdat[ithread].cy_r0C;	_cy_r2C[ithread] = tdat[ithread].cy_r2C;
			_cy_r0D[ithread] = tdat[ithread].cy_r0D;	_cy_r2D[ithread] = tdat[ithread].cy_r2D;
			_cy_r0E[ithread] = tdat[ithread].cy_r0E;	_cy_r2E[ithread] = tdat[ithread].cy_r2E;
			_cy_r0F[ithread] = tdat[ithread].cy_r0F;	_cy_r2F[ithread] = tdat[ithread].cy_r2F;
			_cy_r10[ithread] = tdat[ithread].cy_r10;	_cy_r30[ithread] = tdat[ithread].cy_r30;
			_cy_r11[ithread] = tdat[ithread].cy_r11;	_cy_r31[ithread] = tdat[ithread].cy_r31;
			_cy_r12[ithread] = tdat[ithread].cy_r12;	_cy_r32[ithread] = tdat[ithread].cy_r32;
			_cy_r13[ithread] = tdat[ithread].cy_r13;	_cy_r33[ithread] = tdat[ithread].cy_r33;
			_cy_r14[ithread] = tdat[ithread].cy_r14;	_cy_r34[ithread] = tdat[ithread].cy_r34;
			_cy_r15[ithread] = tdat[ithread].cy_r15;	_cy_r35[ithread] = tdat[ithread].cy_r35;
			_cy_r16[ithread] = tdat[ithread].cy_r16;	_cy_r36[ithread] = tdat[ithread].cy_r36;
			_cy_r17[ithread] = tdat[ithread].cy_r17;	_cy_r37[ithread] = tdat[ithread].cy_r37;
			_cy_r18[ithread] = tdat[ithread].cy_r18;	_cy_r38[ithread] = tdat[ithread].cy_r38;
			_cy_r19[ithread] = tdat[ithread].cy_r19;	_cy_r39[ithread] = tdat[ithread].cy_r39;
			_cy_r1A[ithread] = tdat[ithread].cy_r1A;	_cy_r3A[ithread] = tdat[ithread].cy_r3A;
			_cy_r1B[ithread] = tdat[ithread].cy_r1B;	_cy_r3B[ithread] = tdat[ithread].cy_r3B;
			_cy_r1C[ithread] = tdat[ithread].cy_r1C;	_cy_r3C[ithread] = tdat[ithread].cy_r3C;
			_cy_r1D[ithread] = tdat[ithread].cy_r1D;	_cy_r3D[ithread] = tdat[ithread].cy_r3D;
			_cy_r1E[ithread] = tdat[ithread].cy_r1E;	_cy_r3E[ithread] = tdat[ithread].cy_r3E;
			_cy_r1F[ithread] = tdat[ithread].cy_r1F;	_cy_r3F[ithread] = tdat[ithread].cy_r3F;
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
			_cy_r0A[ithread] = tdat[ithread].cy_r0A;	_cy_i0A[ithread] = tdat[ithread].cy_i0A;
			_cy_r0B[ithread] = tdat[ithread].cy_r0B;	_cy_i0B[ithread] = tdat[ithread].cy_i0B;
			_cy_r0C[ithread] = tdat[ithread].cy_r0C;	_cy_i0C[ithread] = tdat[ithread].cy_i0C;
			_cy_r0D[ithread] = tdat[ithread].cy_r0D;	_cy_i0D[ithread] = tdat[ithread].cy_i0D;
			_cy_r0E[ithread] = tdat[ithread].cy_r0E;	_cy_i0E[ithread] = tdat[ithread].cy_i0E;
			_cy_r0F[ithread] = tdat[ithread].cy_r0F;	_cy_i0F[ithread] = tdat[ithread].cy_i0F;
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
			_cy_r1A[ithread] = tdat[ithread].cy_r1A;	_cy_i1A[ithread] = tdat[ithread].cy_i1A;
			_cy_r1B[ithread] = tdat[ithread].cy_r1B;	_cy_i1B[ithread] = tdat[ithread].cy_i1B;
			_cy_r1C[ithread] = tdat[ithread].cy_r1C;	_cy_i1C[ithread] = tdat[ithread].cy_i1C;
			_cy_r1D[ithread] = tdat[ithread].cy_r1D;	_cy_i1D[ithread] = tdat[ithread].cy_i1D;
			_cy_r1E[ithread] = tdat[ithread].cy_r1E;	_cy_i1E[ithread] = tdat[ithread].cy_i1E;
			_cy_r1F[ithread] = tdat[ithread].cy_r1F;	_cy_i1F[ithread] = tdat[ithread].cy_i1F;

			_cy_r20[ithread] = tdat[ithread].cy_r20;	_cy_i20[ithread] = tdat[ithread].cy_i20;
			_cy_r21[ithread] = tdat[ithread].cy_r21;	_cy_i21[ithread] = tdat[ithread].cy_i21;
			_cy_r22[ithread] = tdat[ithread].cy_r22;	_cy_i22[ithread] = tdat[ithread].cy_i22;
			_cy_r23[ithread] = tdat[ithread].cy_r23;	_cy_i23[ithread] = tdat[ithread].cy_i23;
			_cy_r24[ithread] = tdat[ithread].cy_r24;	_cy_i24[ithread] = tdat[ithread].cy_i24;
			_cy_r25[ithread] = tdat[ithread].cy_r25;	_cy_i25[ithread] = tdat[ithread].cy_i25;
			_cy_r26[ithread] = tdat[ithread].cy_r26;	_cy_i26[ithread] = tdat[ithread].cy_i26;
			_cy_r27[ithread] = tdat[ithread].cy_r27;	_cy_i27[ithread] = tdat[ithread].cy_i27;
			_cy_r28[ithread] = tdat[ithread].cy_r28;	_cy_i28[ithread] = tdat[ithread].cy_i28;
			_cy_r29[ithread] = tdat[ithread].cy_r29;	_cy_i29[ithread] = tdat[ithread].cy_i29;
			_cy_r2A[ithread] = tdat[ithread].cy_r2A;	_cy_i2A[ithread] = tdat[ithread].cy_i2A;
			_cy_r2B[ithread] = tdat[ithread].cy_r2B;	_cy_i2B[ithread] = tdat[ithread].cy_i2B;
			_cy_r2C[ithread] = tdat[ithread].cy_r2C;	_cy_i2C[ithread] = tdat[ithread].cy_i2C;
			_cy_r2D[ithread] = tdat[ithread].cy_r2D;	_cy_i2D[ithread] = tdat[ithread].cy_i2D;
			_cy_r2E[ithread] = tdat[ithread].cy_r2E;	_cy_i2E[ithread] = tdat[ithread].cy_i2E;
			_cy_r2F[ithread] = tdat[ithread].cy_r2F;	_cy_i2F[ithread] = tdat[ithread].cy_i2F;
			_cy_r30[ithread] = tdat[ithread].cy_r30;	_cy_i30[ithread] = tdat[ithread].cy_i30;
			_cy_r31[ithread] = tdat[ithread].cy_r31;	_cy_i31[ithread] = tdat[ithread].cy_i31;
			_cy_r32[ithread] = tdat[ithread].cy_r32;	_cy_i32[ithread] = tdat[ithread].cy_i32;
			_cy_r33[ithread] = tdat[ithread].cy_r33;	_cy_i33[ithread] = tdat[ithread].cy_i33;
			_cy_r34[ithread] = tdat[ithread].cy_r34;	_cy_i34[ithread] = tdat[ithread].cy_i34;
			_cy_r35[ithread] = tdat[ithread].cy_r35;	_cy_i35[ithread] = tdat[ithread].cy_i35;
			_cy_r36[ithread] = tdat[ithread].cy_r36;	_cy_i36[ithread] = tdat[ithread].cy_i36;
			_cy_r37[ithread] = tdat[ithread].cy_r37;	_cy_i37[ithread] = tdat[ithread].cy_i37;
			_cy_r38[ithread] = tdat[ithread].cy_r38;	_cy_i38[ithread] = tdat[ithread].cy_i38;
			_cy_r39[ithread] = tdat[ithread].cy_r39;	_cy_i39[ithread] = tdat[ithread].cy_i39;
			_cy_r3A[ithread] = tdat[ithread].cy_r3A;	_cy_i3A[ithread] = tdat[ithread].cy_i3A;
			_cy_r3B[ithread] = tdat[ithread].cy_r3B;	_cy_i3B[ithread] = tdat[ithread].cy_i3B;
			_cy_r3C[ithread] = tdat[ithread].cy_r3C;	_cy_i3C[ithread] = tdat[ithread].cy_i3C;
			_cy_r3D[ithread] = tdat[ithread].cy_r3D;	_cy_i3D[ithread] = tdat[ithread].cy_i3D;
			_cy_r3E[ithread] = tdat[ithread].cy_r3E;	_cy_i3E[ithread] = tdat[ithread].cy_i3E;
			_cy_r3F[ithread] = tdat[ithread].cy_r3F;	_cy_i3F[ithread] = tdat[ithread].cy_i3F;
		}
	}

#endif

	if(full_pass) {
	//	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

	/*   Wraparound carry cleanup loop is here:

	The cleanup carries from the end of each length-N/32 block into the begining of the next
	can all be neatly processed as follows:

	(1) Invert the radix-64 forward DIF FFT of the first block of 32 complex elements in A and unweight;
	(2) Propagate cleanup carries among the real and imaginary parts of the 32 outputs of (1);
	(3) Reweight and perform a radix-64 forward DIF FFT on the result of (2);
	(4) If any of the exit carries from (2) are nonzero, advance to the next 32 elements and repeat (1-4).
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t00= _cy_r00[CY_THREADS - 1];	t40= _cy_r20[CY_THREADS - 1];
		t02= _cy_r01[CY_THREADS - 1];	t42= _cy_r21[CY_THREADS - 1];
		t04= _cy_r02[CY_THREADS - 1];	t44= _cy_r22[CY_THREADS - 1];
		t06= _cy_r03[CY_THREADS - 1];	t46= _cy_r23[CY_THREADS - 1];
		t08= _cy_r04[CY_THREADS - 1];	t48= _cy_r24[CY_THREADS - 1];
		t0A= _cy_r05[CY_THREADS - 1];	t4A= _cy_r25[CY_THREADS - 1];
		t0C= _cy_r06[CY_THREADS - 1];	t4C= _cy_r26[CY_THREADS - 1];
		t0E= _cy_r07[CY_THREADS - 1];	t4E= _cy_r27[CY_THREADS - 1];
		t10= _cy_r08[CY_THREADS - 1];	t50= _cy_r28[CY_THREADS - 1];
		t12= _cy_r09[CY_THREADS - 1];	t52= _cy_r29[CY_THREADS - 1];
		t14= _cy_r0A[CY_THREADS - 1];	t54= _cy_r2A[CY_THREADS - 1];
		t16= _cy_r0B[CY_THREADS - 1];	t56= _cy_r2B[CY_THREADS - 1];
		t18= _cy_r0C[CY_THREADS - 1];	t58= _cy_r2C[CY_THREADS - 1];
		t1A= _cy_r0D[CY_THREADS - 1];	t5A= _cy_r2D[CY_THREADS - 1];
		t1C= _cy_r0E[CY_THREADS - 1];	t5C= _cy_r2E[CY_THREADS - 1];
		t1E= _cy_r0F[CY_THREADS - 1];	t5E= _cy_r2F[CY_THREADS - 1];
		t20= _cy_r10[CY_THREADS - 1];	t60= _cy_r30[CY_THREADS - 1];
		t22= _cy_r11[CY_THREADS - 1];	t62= _cy_r31[CY_THREADS - 1];
		t24= _cy_r12[CY_THREADS - 1];	t64= _cy_r32[CY_THREADS - 1];
		t26= _cy_r13[CY_THREADS - 1];	t66= _cy_r33[CY_THREADS - 1];
		t28= _cy_r14[CY_THREADS - 1];	t68= _cy_r34[CY_THREADS - 1];
		t2A= _cy_r15[CY_THREADS - 1];	t6A= _cy_r35[CY_THREADS - 1];
		t2C= _cy_r16[CY_THREADS - 1];	t6C= _cy_r36[CY_THREADS - 1];
		t2E= _cy_r17[CY_THREADS - 1];	t6E= _cy_r37[CY_THREADS - 1];
		t30= _cy_r18[CY_THREADS - 1];	t70= _cy_r38[CY_THREADS - 1];
		t32= _cy_r19[CY_THREADS - 1];	t72= _cy_r39[CY_THREADS - 1];
		t34= _cy_r1A[CY_THREADS - 1];	t74= _cy_r3A[CY_THREADS - 1];
		t36= _cy_r1B[CY_THREADS - 1];	t76= _cy_r3B[CY_THREADS - 1];
		t38= _cy_r1C[CY_THREADS - 1];	t78= _cy_r3C[CY_THREADS - 1];
		t3A= _cy_r1D[CY_THREADS - 1];	t7A= _cy_r3D[CY_THREADS - 1];
		t3C= _cy_r1E[CY_THREADS - 1];	t7C= _cy_r3E[CY_THREADS - 1];
		t3E= _cy_r1F[CY_THREADS - 1];	t7E= _cy_r3F[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1,"cy_thread count error!");	/* Make sure loop only gets executed if multiple threads */
			_cy_r00[ithread] = _cy_r00[ithread-1];	_cy_r20[ithread] = _cy_r20[ithread-1];
			_cy_r01[ithread] = _cy_r01[ithread-1];	_cy_r21[ithread] = _cy_r21[ithread-1];
			_cy_r02[ithread] = _cy_r02[ithread-1];	_cy_r22[ithread] = _cy_r22[ithread-1];
			_cy_r03[ithread] = _cy_r03[ithread-1];	_cy_r23[ithread] = _cy_r23[ithread-1];
			_cy_r04[ithread] = _cy_r04[ithread-1];	_cy_r24[ithread] = _cy_r24[ithread-1];
			_cy_r05[ithread] = _cy_r05[ithread-1];	_cy_r25[ithread] = _cy_r25[ithread-1];
			_cy_r06[ithread] = _cy_r06[ithread-1];	_cy_r26[ithread] = _cy_r26[ithread-1];
			_cy_r07[ithread] = _cy_r07[ithread-1];	_cy_r27[ithread] = _cy_r27[ithread-1];
			_cy_r08[ithread] = _cy_r08[ithread-1];	_cy_r28[ithread] = _cy_r28[ithread-1];
			_cy_r09[ithread] = _cy_r09[ithread-1];	_cy_r29[ithread] = _cy_r29[ithread-1];
			_cy_r0A[ithread] = _cy_r0A[ithread-1];	_cy_r2A[ithread] = _cy_r2A[ithread-1];
			_cy_r0B[ithread] = _cy_r0B[ithread-1];	_cy_r2B[ithread] = _cy_r2B[ithread-1];
			_cy_r0C[ithread] = _cy_r0C[ithread-1];	_cy_r2C[ithread] = _cy_r2C[ithread-1];
			_cy_r0D[ithread] = _cy_r0D[ithread-1];	_cy_r2D[ithread] = _cy_r2D[ithread-1];
			_cy_r0E[ithread] = _cy_r0E[ithread-1];	_cy_r2E[ithread] = _cy_r2E[ithread-1];
			_cy_r0F[ithread] = _cy_r0F[ithread-1];	_cy_r2F[ithread] = _cy_r2F[ithread-1];
			_cy_r10[ithread] = _cy_r10[ithread-1];	_cy_r30[ithread] = _cy_r30[ithread-1];
			_cy_r11[ithread] = _cy_r11[ithread-1];	_cy_r31[ithread] = _cy_r31[ithread-1];
			_cy_r12[ithread] = _cy_r12[ithread-1];	_cy_r32[ithread] = _cy_r32[ithread-1];
			_cy_r13[ithread] = _cy_r13[ithread-1];	_cy_r33[ithread] = _cy_r33[ithread-1];
			_cy_r14[ithread] = _cy_r14[ithread-1];	_cy_r34[ithread] = _cy_r34[ithread-1];
			_cy_r15[ithread] = _cy_r15[ithread-1];	_cy_r35[ithread] = _cy_r35[ithread-1];
			_cy_r16[ithread] = _cy_r16[ithread-1];	_cy_r36[ithread] = _cy_r36[ithread-1];
			_cy_r17[ithread] = _cy_r17[ithread-1];	_cy_r37[ithread] = _cy_r37[ithread-1];
			_cy_r18[ithread] = _cy_r18[ithread-1];	_cy_r38[ithread] = _cy_r38[ithread-1];
			_cy_r19[ithread] = _cy_r19[ithread-1];	_cy_r39[ithread] = _cy_r39[ithread-1];
			_cy_r1A[ithread] = _cy_r1A[ithread-1];	_cy_r3A[ithread] = _cy_r3A[ithread-1];
			_cy_r1B[ithread] = _cy_r1B[ithread-1];	_cy_r3B[ithread] = _cy_r3B[ithread-1];
			_cy_r1C[ithread] = _cy_r1C[ithread-1];	_cy_r3C[ithread] = _cy_r3C[ithread-1];
			_cy_r1D[ithread] = _cy_r1D[ithread-1];	_cy_r3D[ithread] = _cy_r3D[ithread-1];
			_cy_r1E[ithread] = _cy_r1E[ithread-1];	_cy_r3E[ithread] = _cy_r3E[ithread-1];
			_cy_r1F[ithread] = _cy_r1F[ithread-1];	_cy_r3F[ithread] = _cy_r3F[ithread-1];
		}
		/* ...The wraparound carry is here: */
		_cy_r00[0] =+t7E;	_cy_r20[0] = t3E;
		_cy_r01[0] = t00;	_cy_r21[0] = t40;
		_cy_r02[0] = t02;	_cy_r22[0] = t42;
		_cy_r03[0] = t04;	_cy_r23[0] = t44;
		_cy_r04[0] = t06;	_cy_r24[0] = t46;
		_cy_r05[0] = t08;	_cy_r25[0] = t48;
		_cy_r06[0] = t0A;	_cy_r26[0] = t4A;
		_cy_r07[0] = t0C;	_cy_r27[0] = t4C;
		_cy_r08[0] = t0E;	_cy_r28[0] = t4E;
		_cy_r09[0] = t10;	_cy_r29[0] = t50;
		_cy_r0A[0] = t12;	_cy_r2A[0] = t52;
		_cy_r0B[0] = t14;	_cy_r2B[0] = t54;
		_cy_r0C[0] = t16;	_cy_r2C[0] = t56;
		_cy_r0D[0] = t18;	_cy_r2D[0] = t58;
		_cy_r0E[0] = t1A;	_cy_r2E[0] = t5A;
		_cy_r0F[0] = t1C;	_cy_r2F[0] = t5C;
		_cy_r10[0] = t1E;	_cy_r30[0] = t5E;
		_cy_r11[0] = t20;	_cy_r31[0] = t60;
		_cy_r12[0] = t22;	_cy_r32[0] = t62;
		_cy_r13[0] = t24;	_cy_r33[0] = t64;
		_cy_r14[0] = t26;	_cy_r34[0] = t66;
		_cy_r15[0] = t28;	_cy_r35[0] = t68;
		_cy_r16[0] = t2A;	_cy_r36[0] = t6A;
		_cy_r17[0] = t2C;	_cy_r37[0] = t6C;
		_cy_r18[0] = t2E;	_cy_r38[0] = t6E;
		_cy_r19[0] = t30;	_cy_r39[0] = t70;
		_cy_r1A[0] = t32;	_cy_r3A[0] = t72;
		_cy_r1B[0] = t34;	_cy_r3B[0] = t74;
		_cy_r1C[0] = t36;	_cy_r3C[0] = t76;
		_cy_r1D[0] = t38;	_cy_r3D[0] = t78;
		_cy_r1E[0] = t3A;	_cy_r3E[0] = t7A;
		_cy_r1F[0] = t3C;	_cy_r3F[0] = t7C;
	}
	else
	{
		t00= _cy_r00[CY_THREADS - 1];	t01= _cy_i00[CY_THREADS - 1];
		t02= _cy_r01[CY_THREADS - 1];	t03= _cy_i01[CY_THREADS - 1];
		t04= _cy_r02[CY_THREADS - 1];	t05= _cy_i02[CY_THREADS - 1];
		t06= _cy_r03[CY_THREADS - 1];	t07= _cy_i03[CY_THREADS - 1];
		t08= _cy_r04[CY_THREADS - 1];	t09= _cy_i04[CY_THREADS - 1];
		t0A= _cy_r05[CY_THREADS - 1];	t0B= _cy_i05[CY_THREADS - 1];
		t0C= _cy_r06[CY_THREADS - 1];	t0D= _cy_i06[CY_THREADS - 1];
		t0E= _cy_r07[CY_THREADS - 1];	t0F= _cy_i07[CY_THREADS - 1];
		t10= _cy_r08[CY_THREADS - 1];	t11= _cy_i08[CY_THREADS - 1];
		t12= _cy_r09[CY_THREADS - 1];	t13= _cy_i09[CY_THREADS - 1];
		t14= _cy_r0A[CY_THREADS - 1];	t15= _cy_i0A[CY_THREADS - 1];
		t16= _cy_r0B[CY_THREADS - 1];	t17= _cy_i0B[CY_THREADS - 1];
		t18= _cy_r0C[CY_THREADS - 1];	t19= _cy_i0C[CY_THREADS - 1];
		t1A= _cy_r0D[CY_THREADS - 1];	t1B= _cy_i0D[CY_THREADS - 1];
		t1C= _cy_r0E[CY_THREADS - 1];	t1D= _cy_i0E[CY_THREADS - 1];
		t1E= _cy_r0F[CY_THREADS - 1];	t1F= _cy_i0F[CY_THREADS - 1];
		t20= _cy_r10[CY_THREADS - 1];	t21= _cy_i10[CY_THREADS - 1];
		t22= _cy_r11[CY_THREADS - 1];	t23= _cy_i11[CY_THREADS - 1];
		t24= _cy_r12[CY_THREADS - 1];	t25= _cy_i12[CY_THREADS - 1];
		t26= _cy_r13[CY_THREADS - 1];	t27= _cy_i13[CY_THREADS - 1];
		t28= _cy_r14[CY_THREADS - 1];	t29= _cy_i14[CY_THREADS - 1];
		t2A= _cy_r15[CY_THREADS - 1];	t2B= _cy_i15[CY_THREADS - 1];
		t2C= _cy_r16[CY_THREADS - 1];	t2D= _cy_i16[CY_THREADS - 1];
		t2E= _cy_r17[CY_THREADS - 1];	t2F= _cy_i17[CY_THREADS - 1];
		t30= _cy_r18[CY_THREADS - 1];	t31= _cy_i18[CY_THREADS - 1];
		t32= _cy_r19[CY_THREADS - 1];	t33= _cy_i19[CY_THREADS - 1];
		t34= _cy_r1A[CY_THREADS - 1];	t35= _cy_i1A[CY_THREADS - 1];
		t36= _cy_r1B[CY_THREADS - 1];	t37= _cy_i1B[CY_THREADS - 1];
		t38= _cy_r1C[CY_THREADS - 1];	t39= _cy_i1C[CY_THREADS - 1];
		t3A= _cy_r1D[CY_THREADS - 1];	t3B= _cy_i1D[CY_THREADS - 1];
		t3C= _cy_r1E[CY_THREADS - 1];	t3D= _cy_i1E[CY_THREADS - 1];
		t3E= _cy_r1F[CY_THREADS - 1];	t3F= _cy_i1F[CY_THREADS - 1];
		t40= _cy_r20[CY_THREADS - 1];	t41= _cy_i20[CY_THREADS - 1];
		t42= _cy_r21[CY_THREADS - 1];	t43= _cy_i21[CY_THREADS - 1];
		t44= _cy_r22[CY_THREADS - 1];	t45= _cy_i22[CY_THREADS - 1];
		t46= _cy_r23[CY_THREADS - 1];	t47= _cy_i23[CY_THREADS - 1];
		t48= _cy_r24[CY_THREADS - 1];	t49= _cy_i24[CY_THREADS - 1];
		t4A= _cy_r25[CY_THREADS - 1];	t4B= _cy_i25[CY_THREADS - 1];
		t4C= _cy_r26[CY_THREADS - 1];	t4D= _cy_i26[CY_THREADS - 1];
		t4E= _cy_r27[CY_THREADS - 1];	t4F= _cy_i27[CY_THREADS - 1];
		t50= _cy_r28[CY_THREADS - 1];	t51= _cy_i28[CY_THREADS - 1];
		t52= _cy_r29[CY_THREADS - 1];	t53= _cy_i29[CY_THREADS - 1];
		t54= _cy_r2A[CY_THREADS - 1];	t55= _cy_i2A[CY_THREADS - 1];
		t56= _cy_r2B[CY_THREADS - 1];	t57= _cy_i2B[CY_THREADS - 1];
		t58= _cy_r2C[CY_THREADS - 1];	t59= _cy_i2C[CY_THREADS - 1];
		t5A= _cy_r2D[CY_THREADS - 1];	t5B= _cy_i2D[CY_THREADS - 1];
		t5C= _cy_r2E[CY_THREADS - 1];	t5D= _cy_i2E[CY_THREADS - 1];
		t5E= _cy_r2F[CY_THREADS - 1];	t5F= _cy_i2F[CY_THREADS - 1];
		t60= _cy_r30[CY_THREADS - 1];	t61= _cy_i30[CY_THREADS - 1];
		t62= _cy_r31[CY_THREADS - 1];	t63= _cy_i31[CY_THREADS - 1];
		t64= _cy_r32[CY_THREADS - 1];	t65= _cy_i32[CY_THREADS - 1];
		t66= _cy_r33[CY_THREADS - 1];	t67= _cy_i33[CY_THREADS - 1];
		t68= _cy_r34[CY_THREADS - 1];	t69= _cy_i34[CY_THREADS - 1];
		t6A= _cy_r35[CY_THREADS - 1];	t6B= _cy_i35[CY_THREADS - 1];
		t6C= _cy_r36[CY_THREADS - 1];	t6D= _cy_i36[CY_THREADS - 1];
		t6E= _cy_r37[CY_THREADS - 1];	t6F= _cy_i37[CY_THREADS - 1];
		t70= _cy_r38[CY_THREADS - 1];	t71= _cy_i38[CY_THREADS - 1];
		t72= _cy_r39[CY_THREADS - 1];	t73= _cy_i39[CY_THREADS - 1];
		t74= _cy_r3A[CY_THREADS - 1];	t75= _cy_i3A[CY_THREADS - 1];
		t76= _cy_r3B[CY_THREADS - 1];	t77= _cy_i3B[CY_THREADS - 1];
		t78= _cy_r3C[CY_THREADS - 1];	t79= _cy_i3C[CY_THREADS - 1];
		t7A= _cy_r3D[CY_THREADS - 1];	t7B= _cy_i3D[CY_THREADS - 1];
		t7C= _cy_r3E[CY_THREADS - 1];	t7D= _cy_i3E[CY_THREADS - 1];
		t7E= _cy_r3F[CY_THREADS - 1];	t7F= _cy_i3F[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1,"cy_thread count error!");	/* Make sure loop only gets executed if multiple threads */
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
			_cy_r0A[ithread] = _cy_r0A[ithread-1];		_cy_i0A[ithread] = _cy_i0A[ithread-1];
			_cy_r0B[ithread] = _cy_r0B[ithread-1];		_cy_i0B[ithread] = _cy_i0B[ithread-1];
			_cy_r0C[ithread] = _cy_r0C[ithread-1];		_cy_i0C[ithread] = _cy_i0C[ithread-1];
			_cy_r0D[ithread] = _cy_r0D[ithread-1];		_cy_i0D[ithread] = _cy_i0D[ithread-1];
			_cy_r0E[ithread] = _cy_r0E[ithread-1];		_cy_i0E[ithread] = _cy_i0E[ithread-1];
			_cy_r0F[ithread] = _cy_r0F[ithread-1];		_cy_i0F[ithread] = _cy_i0F[ithread-1];
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
			_cy_r1A[ithread] = _cy_r1A[ithread-1];		_cy_i1A[ithread] = _cy_i1A[ithread-1];
			_cy_r1B[ithread] = _cy_r1B[ithread-1];		_cy_i1B[ithread] = _cy_i1B[ithread-1];
			_cy_r1C[ithread] = _cy_r1C[ithread-1];		_cy_i1C[ithread] = _cy_i1C[ithread-1];
			_cy_r1D[ithread] = _cy_r1D[ithread-1];		_cy_i1D[ithread] = _cy_i1D[ithread-1];
			_cy_r1E[ithread] = _cy_r1E[ithread-1];		_cy_i1E[ithread] = _cy_i1E[ithread-1];
			_cy_r1F[ithread] = _cy_r1F[ithread-1];		_cy_i1F[ithread] = _cy_i1F[ithread-1];
			_cy_r20[ithread] = _cy_r20[ithread-1];		_cy_i20[ithread] = _cy_i20[ithread-1];
			_cy_r21[ithread] = _cy_r21[ithread-1];		_cy_i21[ithread] = _cy_i21[ithread-1];
			_cy_r22[ithread] = _cy_r22[ithread-1];		_cy_i22[ithread] = _cy_i22[ithread-1];
			_cy_r23[ithread] = _cy_r23[ithread-1];		_cy_i23[ithread] = _cy_i23[ithread-1];
			_cy_r24[ithread] = _cy_r24[ithread-1];		_cy_i24[ithread] = _cy_i24[ithread-1];
			_cy_r25[ithread] = _cy_r25[ithread-1];		_cy_i25[ithread] = _cy_i25[ithread-1];
			_cy_r26[ithread] = _cy_r26[ithread-1];		_cy_i26[ithread] = _cy_i26[ithread-1];
			_cy_r27[ithread] = _cy_r27[ithread-1];		_cy_i27[ithread] = _cy_i27[ithread-1];
			_cy_r28[ithread] = _cy_r28[ithread-1];		_cy_i28[ithread] = _cy_i28[ithread-1];
			_cy_r29[ithread] = _cy_r29[ithread-1];		_cy_i29[ithread] = _cy_i29[ithread-1];
			_cy_r2A[ithread] = _cy_r2A[ithread-1];		_cy_i2A[ithread] = _cy_i2A[ithread-1];
			_cy_r2B[ithread] = _cy_r2B[ithread-1];		_cy_i2B[ithread] = _cy_i2B[ithread-1];
			_cy_r2C[ithread] = _cy_r2C[ithread-1];		_cy_i2C[ithread] = _cy_i2C[ithread-1];
			_cy_r2D[ithread] = _cy_r2D[ithread-1];		_cy_i2D[ithread] = _cy_i2D[ithread-1];
			_cy_r2E[ithread] = _cy_r2E[ithread-1];		_cy_i2E[ithread] = _cy_i2E[ithread-1];
			_cy_r2F[ithread] = _cy_r2F[ithread-1];		_cy_i2F[ithread] = _cy_i2F[ithread-1];
			_cy_r30[ithread] = _cy_r30[ithread-1];		_cy_i30[ithread] = _cy_i30[ithread-1];
			_cy_r31[ithread] = _cy_r31[ithread-1];		_cy_i31[ithread] = _cy_i31[ithread-1];
			_cy_r32[ithread] = _cy_r32[ithread-1];		_cy_i32[ithread] = _cy_i32[ithread-1];
			_cy_r33[ithread] = _cy_r33[ithread-1];		_cy_i33[ithread] = _cy_i33[ithread-1];
			_cy_r34[ithread] = _cy_r34[ithread-1];		_cy_i34[ithread] = _cy_i34[ithread-1];
			_cy_r35[ithread] = _cy_r35[ithread-1];		_cy_i35[ithread] = _cy_i35[ithread-1];
			_cy_r36[ithread] = _cy_r36[ithread-1];		_cy_i36[ithread] = _cy_i36[ithread-1];
			_cy_r37[ithread] = _cy_r37[ithread-1];		_cy_i37[ithread] = _cy_i37[ithread-1];
			_cy_r38[ithread] = _cy_r38[ithread-1];		_cy_i38[ithread] = _cy_i38[ithread-1];
			_cy_r39[ithread] = _cy_r39[ithread-1];		_cy_i39[ithread] = _cy_i39[ithread-1];
			_cy_r3A[ithread] = _cy_r3A[ithread-1];		_cy_i3A[ithread] = _cy_i3A[ithread-1];
			_cy_r3B[ithread] = _cy_r3B[ithread-1];		_cy_i3B[ithread] = _cy_i3B[ithread-1];
			_cy_r3C[ithread] = _cy_r3C[ithread-1];		_cy_i3C[ithread] = _cy_i3C[ithread-1];
			_cy_r3D[ithread] = _cy_r3D[ithread-1];		_cy_i3D[ithread] = _cy_i3D[ithread-1];
			_cy_r3E[ithread] = _cy_r3E[ithread-1];		_cy_i3E[ithread] = _cy_i3E[ithread-1];
			_cy_r3F[ithread] = _cy_r3F[ithread-1];		_cy_i3F[ithread] = _cy_i3F[ithread-1];
		}
		/* ...The 2 Mo"bius carries are here: */
		_cy_r00[0] =-t7F;	_cy_i00[0] =+t7E;
		_cy_r01[0] = t00;	_cy_i01[0] = t01;
		_cy_r02[0] = t02;	_cy_i02[0] = t03;
		_cy_r03[0] = t04;	_cy_i03[0] = t05;
		_cy_r04[0] = t06;	_cy_i04[0] = t07;
		_cy_r05[0] = t08;	_cy_i05[0] = t09;
		_cy_r06[0] = t0A;	_cy_i06[0] = t0B;
		_cy_r07[0] = t0C;	_cy_i07[0] = t0D;
		_cy_r08[0] = t0E;	_cy_i08[0] = t0F;
		_cy_r09[0] = t10;	_cy_i09[0] = t11;
		_cy_r0A[0] = t12;	_cy_i0A[0] = t13;
		_cy_r0B[0] = t14;	_cy_i0B[0] = t15;
		_cy_r0C[0] = t16;	_cy_i0C[0] = t17;
		_cy_r0D[0] = t18;	_cy_i0D[0] = t19;
		_cy_r0E[0] = t1A;	_cy_i0E[0] = t1B;
		_cy_r0F[0] = t1C;	_cy_i0F[0] = t1D;
		_cy_r10[0] = t1E;	_cy_i10[0] = t1F;
		_cy_r11[0] = t20;	_cy_i11[0] = t21;
		_cy_r12[0] = t22;	_cy_i12[0] = t23;
		_cy_r13[0] = t24;	_cy_i13[0] = t25;
		_cy_r14[0] = t26;	_cy_i14[0] = t27;
		_cy_r15[0] = t28;	_cy_i15[0] = t29;
		_cy_r16[0] = t2A;	_cy_i16[0] = t2B;
		_cy_r17[0] = t2C;	_cy_i17[0] = t2D;
		_cy_r18[0] = t2E;	_cy_i18[0] = t2F;
		_cy_r19[0] = t30;	_cy_i19[0] = t31;
		_cy_r1A[0] = t32;	_cy_i1A[0] = t33;
		_cy_r1B[0] = t34;	_cy_i1B[0] = t35;
		_cy_r1C[0] = t36;	_cy_i1C[0] = t37;
		_cy_r1D[0] = t38;	_cy_i1D[0] = t39;
		_cy_r1E[0] = t3A;	_cy_i1E[0] = t3B;
		_cy_r1F[0] = t3C;	_cy_i1F[0] = t3D;
		_cy_r20[0] = t3E;	_cy_i20[0] = t3F;
		_cy_r21[0] = t40;	_cy_i21[0] = t41;
		_cy_r22[0] = t42;	_cy_i22[0] = t43;
		_cy_r23[0] = t44;	_cy_i23[0] = t45;
		_cy_r24[0] = t46;	_cy_i24[0] = t47;
		_cy_r25[0] = t48;	_cy_i25[0] = t49;
		_cy_r26[0] = t4A;	_cy_i26[0] = t4B;
		_cy_r27[0] = t4C;	_cy_i27[0] = t4D;
		_cy_r28[0] = t4E;	_cy_i28[0] = t4F;
		_cy_r29[0] = t50;	_cy_i29[0] = t51;
		_cy_r2A[0] = t52;	_cy_i2A[0] = t53;
		_cy_r2B[0] = t54;	_cy_i2B[0] = t55;
		_cy_r2C[0] = t56;	_cy_i2C[0] = t57;
		_cy_r2D[0] = t58;	_cy_i2D[0] = t59;
		_cy_r2E[0] = t5A;	_cy_i2E[0] = t5B;
		_cy_r2F[0] = t5C;	_cy_i2F[0] = t5D;
		_cy_r30[0] = t5E;	_cy_i30[0] = t5F;
		_cy_r31[0] = t60;	_cy_i31[0] = t61;
		_cy_r32[0] = t62;	_cy_i32[0] = t63;
		_cy_r33[0] = t64;	_cy_i33[0] = t65;
		_cy_r34[0] = t66;	_cy_i34[0] = t67;
		_cy_r35[0] = t68;	_cy_i35[0] = t69;
		_cy_r36[0] = t6A;	_cy_i36[0] = t6B;
		_cy_r37[0] = t6C;	_cy_i37[0] = t6D;
		_cy_r38[0] = t6E;	_cy_i38[0] = t6F;
		_cy_r39[0] = t70;	_cy_i39[0] = t71;
		_cy_r3A[0] = t72;	_cy_i3A[0] = t73;
		_cy_r3B[0] = t74;	_cy_i3B[0] = t75;
		_cy_r3C[0] = t76;	_cy_i3C[0] = t77;
		_cy_r3D[0] = t78;	_cy_i3D[0] = t79;
		_cy_r3E[0] = t7A;	_cy_i3E[0] = t7B;
		_cy_r3F[0] = t7C;	_cy_i3F[0] = t7D;
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
			jt = j;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			a[jt+p04] *= radix_inv;
			a[jt+p05] *= radix_inv;
			a[jt+p06] *= radix_inv;
			a[jt+p07] *= radix_inv;
			jt = j + p08;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			a[jt+p04] *= radix_inv;
			a[jt+p05] *= radix_inv;
			a[jt+p06] *= radix_inv;
			a[jt+p07] *= radix_inv;
			jt = j + p10;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			a[jt+p04] *= radix_inv;
			a[jt+p05] *= radix_inv;
			a[jt+p06] *= radix_inv;
			a[jt+p07] *= radix_inv;
			jt = j + p18;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			a[jt+p04] *= radix_inv;
			a[jt+p05] *= radix_inv;
			a[jt+p06] *= radix_inv;
			a[jt+p07] *= radix_inv;
			jt = j + p20;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			a[jt+p04] *= radix_inv;
			a[jt+p05] *= radix_inv;
			a[jt+p06] *= radix_inv;
			a[jt+p07] *= radix_inv;
			jt = j + p28;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			a[jt+p04] *= radix_inv;
			a[jt+p05] *= radix_inv;
			a[jt+p06] *= radix_inv;
			a[jt+p07] *= radix_inv;
			jt = j + p30;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			a[jt+p04] *= radix_inv;
			a[jt+p05] *= radix_inv;
			a[jt+p06] *= radix_inv;
			a[jt+p07] *= radix_inv;
			jt = j + p38;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			a[jt+p04] *= radix_inv;
			a[jt+p05] *= radix_inv;
			a[jt+p06] *= radix_inv;
			a[jt+p07] *= radix_inv;
		}
	}
}	/* endfor(outer) */

	t00 = 0;
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		t00 += fabs(_cy_r00[0])+fabs(_cy_r01[0])+fabs(_cy_r02[0])+fabs(_cy_r03[0])+fabs(_cy_r04[0])+fabs(_cy_r05[0])+fabs(_cy_r06[0])+fabs(_cy_r07[0])+fabs(_cy_r08[0])+fabs(_cy_r09[0])+fabs(_cy_r0A[0])+fabs(_cy_r0B[0])+fabs(_cy_r0C[0])+fabs(_cy_r0D[0])+fabs(_cy_r0E[0])+fabs(_cy_r0F[0])
			  +fabs(_cy_r10[0])+fabs(_cy_r11[0])+fabs(_cy_r12[0])+fabs(_cy_r13[0])+fabs(_cy_r14[0])+fabs(_cy_r15[0])+fabs(_cy_r16[0])+fabs(_cy_r17[0])+fabs(_cy_r18[0])+fabs(_cy_r19[0])+fabs(_cy_r1A[0])+fabs(_cy_r1B[0])+fabs(_cy_r1C[0])+fabs(_cy_r1D[0])+fabs(_cy_r1E[0])+fabs(_cy_r1F[0])
			  +fabs(_cy_r20[0])+fabs(_cy_r21[0])+fabs(_cy_r22[0])+fabs(_cy_r23[0])+fabs(_cy_r24[0])+fabs(_cy_r25[0])+fabs(_cy_r26[0])+fabs(_cy_r27[0])+fabs(_cy_r28[0])+fabs(_cy_r29[0])+fabs(_cy_r2A[0])+fabs(_cy_r2B[0])+fabs(_cy_r2C[0])+fabs(_cy_r2D[0])+fabs(_cy_r2E[0])+fabs(_cy_r2F[0])
			  +fabs(_cy_r30[0])+fabs(_cy_r31[0])+fabs(_cy_r32[0])+fabs(_cy_r33[0])+fabs(_cy_r34[0])+fabs(_cy_r35[0])+fabs(_cy_r36[0])+fabs(_cy_r37[0])+fabs(_cy_r38[0])+fabs(_cy_r39[0])+fabs(_cy_r3A[0])+fabs(_cy_r3B[0])+fabs(_cy_r3C[0])+fabs(_cy_r3D[0])+fabs(_cy_r3E[0])+fabs(_cy_r3F[0]);
		t00 += fabs(_cy_i00[0])+fabs(_cy_i01[0])+fabs(_cy_i02[0])+fabs(_cy_i03[0])+fabs(_cy_i04[0])+fabs(_cy_i05[0])+fabs(_cy_i06[0])+fabs(_cy_i07[0])+fabs(_cy_i08[0])+fabs(_cy_i09[0])+fabs(_cy_i0A[0])+fabs(_cy_i0B[0])+fabs(_cy_i0C[0])+fabs(_cy_i0D[0])+fabs(_cy_i0E[0])+fabs(_cy_i0F[0])
			  +fabs(_cy_i10[0])+fabs(_cy_i11[0])+fabs(_cy_i12[0])+fabs(_cy_i13[0])+fabs(_cy_i14[0])+fabs(_cy_i15[0])+fabs(_cy_i16[0])+fabs(_cy_i17[0])+fabs(_cy_i18[0])+fabs(_cy_i19[0])+fabs(_cy_i1A[0])+fabs(_cy_i1B[0])+fabs(_cy_i1C[0])+fabs(_cy_i1D[0])+fabs(_cy_i1E[0])+fabs(_cy_i1F[0])
			  +fabs(_cy_i20[0])+fabs(_cy_i21[0])+fabs(_cy_i22[0])+fabs(_cy_i23[0])+fabs(_cy_i24[0])+fabs(_cy_i25[0])+fabs(_cy_i26[0])+fabs(_cy_i27[0])+fabs(_cy_i28[0])+fabs(_cy_i29[0])+fabs(_cy_i2A[0])+fabs(_cy_i2B[0])+fabs(_cy_i2C[0])+fabs(_cy_i2D[0])+fabs(_cy_i2E[0])+fabs(_cy_i2F[0])
			  +fabs(_cy_i30[0])+fabs(_cy_i31[0])+fabs(_cy_i32[0])+fabs(_cy_i33[0])+fabs(_cy_i34[0])+fabs(_cy_i35[0])+fabs(_cy_i36[0])+fabs(_cy_i37[0])+fabs(_cy_i38[0])+fabs(_cy_i39[0])+fabs(_cy_i3A[0])+fabs(_cy_i3B[0])+fabs(_cy_i3C[0])+fabs(_cy_i3D[0])+fabs(_cy_i3E[0])+fabs(_cy_i3F[0]);

		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
	}
//*fracmax = 0;
	if(t00 != 0.0)
	{
		sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry = %20.10e - input wordsize may be too small.\n",iter,t00);
		if(INTERACT)fprintf(stderr,"%s",cbuf);
		fp = fopen(   OFILE,"a");
		fq = fopen(STATFILE,"a");
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
	int j,j1,j2,jt,jp;
	static int n64,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30,p38, first_entry=TRUE;
	static double
			c64_1 = 0.99518472667219688624, s64_1 = 0.09801714032956060199,	/* exp(  I*twopi/64) */
			c64_2 = 0.98078528040323044913, s64_2 = 0.19509032201612826785,	/* exp(2*I*twopi/64) */
			c64_3 = 0.95694033573220886494, s64_3 = 0.29028467725446236764,	/* exp(3*I*twopi/64) */
			c64_4 = 0.92387953251128675613, s64_4 = 0.38268343236508977173,	/* exp(4*I*twopi/64) */
			c64_5 = 0.88192126434835502971, s64_5 = 0.47139673682599764856,	/* exp(5*I*twopi/64) */
			c64_6 = 0.83146961230254523708, s64_6 = 0.55557023301960222474,	/* exp(6*I*twopi/64) */
			c64_7 = 0.77301045336273696081, s64_7 = 0.63439328416364549822;	/* exp(7*I*twopi/64) */
  double t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F
		,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F
		,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F
		,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t4A,t4B,t4C,t4D,t4E,t4F
		,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t5A,t5B,t5C,t5D,t5E,t5F
		,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t6A,t6B,t6C,t6D,t6E,t6F
		,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t7A,t7B,t7C,t7D,t7E,t7F;

	if(!first_entry && (n >> 6) != n64)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n64=n/64;

		p01 = n64;
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
	}

/*...The radix-64 pass is here.	*/

	for(j = 0; j < n64; j += 2)
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
		E^ 2 = exp(i* 2*twopi/64) =       ( c64_2, s64_2)
		E^ 3 = exp(i* 3*twopi/64) =       ( c64_3, s64_3)
		E^ 4 = exp(i* 4*twopi/64) =       ( c64_4, s64_4)
		E^ 5 = exp(i* 5*twopi/64) =       ( c64_5, s64_5)
		E^ 6 = exp(i* 6*twopi/64) =       ( c64_6, s64_6)
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
			1,0,1,0,1,0,1,0,1,0,1,0,1,0
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
			0,1,ISRT2,ISRT2,-ISRT2,ISRT2,c64_4,s64_4,-s64_4,c64_4,s64_4,c64_4,-c64_4,s64_4
		);
		/* Block 2: t*4,5, twiddles = {   1,  E^ 2,  E^ 4,  E^ 6,  E^ 8, *E^ 6, *E^ 4, *E^ 2}
							BR order: {   1,  E^ 8,  E^ 4, *E^ 4,  E^ 2, *E^ 6,  E^ 6, *E^ 2}
			In terms of [Re,Im] pairs:{[1,0],isrt2*[1,1],[c4,s4],[s4,c4],[c2,s2],[s6,c6],[c6,s6],[s2,c2]}
		*/
		jt = j1 + p10;	jp = j2 + p10;
		RADIX_08_DIF_TWIDDLE_OOP(
			t04,t05,t14,t15,t24,t25,t34,t35,t44,t45,t54,t55,t64,t65,t74,t75,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			ISRT2,ISRT2,c64_4,s64_4,s64_4,c64_4,c64_2,s64_2,s64_6,c64_6,c64_6,s64_6,s64_2,c64_2
		);
		/* Block 6: t*C,D, twiddles = {   1,  E^ 6, *E^ 4,*~E^ 2,*~E^ 8,-~E^ 2, -E^ 4,-*E^ 6}
							BR order: {   1,*~E^ 8, *E^ 4, -E^ 4,  E^ 6,-~E^ 2,*~E^ 2,-*E^ 6}
			In terms of [Re,Im] pairs:{[1,0],isrt2*[-1,1],[s4,c4],[-c4,-s4],[c6,s6],[-c2,s2],[-s2,c2],[-s6,-c6]}
		*/
		jt = j1 + p18;	jp = j2 + p18;
		RADIX_08_DIF_TWIDDLE_OOP(
			t0C,t0D,t1C,t1D,t2C,t2D,t3C,t3D,t4C,t4D,t5C,t5D,t6C,t6D,t7C,t7D,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			-ISRT2,ISRT2,s64_4,c64_4,-c64_4,-s64_4,c64_6,s64_6,-c64_2,s64_2,-s64_2,c64_2,-s64_6,-c64_6
		);
		/* Block 1: t*2,3, twiddles = {   1,  E^ 1,  E^ 2,  E^ 3,  E^ 4,  E^ 5,  E^ 6,  E^ 7}
							BR order: {   1,  E^ 4,  E^ 2,  E^ 6,  E^ 1,  E^ 5,  E^ 3,  E^ 7}
			In terms of [Re,Im] pairs:{[1,0],[c4,s4],[c2,s2],[c6,s6],[c1,s1],[c5,s5],[c3,s3],[c7,s7]}
		*/
		jt = j1 + p20;	jp = j2 + p20;
		RADIX_08_DIF_TWIDDLE_OOP(
			t02,t03,t12,t13,t22,t23,t32,t33,t42,t43,t52,t53,t62,t63,t72,t73,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			c64_4,s64_4,c64_2,s64_2,c64_6,s64_6,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7
		);
		/* Block 5: t*A,B, twiddles = {   1,  E^ 5, *E^ 6, *E^ 1,*~E^ 4,-~E^ 7,-~E^ 2, -E^ 3}
							BR order: {   1,*~E^ 4, *E^ 6,-~E^ 2,  E^ 5,-~E^ 7, *E^ 1, -E^ 3}
			In terms of [Re,Im] pairs:{[1,0],[-s4,c4],[s6,c6],[-c2,s2],[c5,s5],[-c7,s7],[s1,c1],[-c3,-s3]}
		*/
		jt = j1 + p28;	jp = j2 + p28;
		RADIX_08_DIF_TWIDDLE_OOP(
			t0A,t0B,t1A,t1B,t2A,t2B,t3A,t3B,t4A,t4B,t5A,t5B,t6A,t6B,t7A,t7B,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			-s64_4,c64_4,s64_6,c64_6,-c64_2,s64_2,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3
		);
		/* Block 3: t*6,7, twiddles = {   1,  E^ 3,  E^ 6, *E^ 7, *E^ 4, *E^ 1,*~E^ 2,*~E^ 5}
							BR order: {   1, *E^ 4,  E^ 6,*~E^ 2,  E^ 3, *E^ 1, *E^ 7,*~E^ 5}
			In terms of [Re,Im] pairs:{[1,0],[s4,c4],[c6,s6],[-s2,c2],[c3,s3],[s1,c1],[s7,c7],[-s5,c5]}
		*/
		jt = j1 + p30;	jp = j2 + p30;
		RADIX_08_DIF_TWIDDLE_OOP(
			t06,t07,t16,t17,t26,t27,t36,t37,t46,t47,t56,t57,t66,t67,t76,t77,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			s64_4,c64_4,c64_6,s64_6,-s64_2,c64_2,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5
		);
		/* Block 7: t*E,F, twiddles = {   1,  E^ 7, *E^ 2,*~E^ 5,-~E^ 4, -E^ 3,-*E^ 6,~*E^ 1} ,
							BR order: {   1,-~E^ 4, *E^ 2,-*E^ 6,  E^ 7, -E^ 3,*~E^ 5,~*E^ 1} ,
			In terms of [Re,Im] pairs:{[1,0],[-c4,s4],[s2,c2],[-s6,-c6],[c7,s7],[-c3,-s3],[-s5,c5],[s1,-c1]}
		*/
		jt = j1 + p38;	jp = j2 + p38;
		RADIX_08_DIF_TWIDDLE_OOP(
			t0E,t0F,t1E,t1F,t2E,t2F,t3E,t3F,t4E,t4F,t5E,t5F,t6E,t6F,t7E,t7F,
			a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
			-c64_4,s64_4,s64_2,c64_2,-s64_6,-c64_6,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1
		);

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
	int j,j1,j2,jt,jp;
	static int n64,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30,p38, first_entry=TRUE;
	static double
			c64_1 = 0.99518472667219688624, s64_1 = 0.09801714032956060199,	/* exp(  I*twopi/64) */
			c64_2 = 0.98078528040323044913, s64_2 = 0.19509032201612826785,	/* exp(2*I*twopi/64) */
			c64_3 = 0.95694033573220886494, s64_3 = 0.29028467725446236764,	/* exp(3*I*twopi/64) */
			c64_4 = 0.92387953251128675613, s64_4 = 0.38268343236508977173,	/* exp(4*I*twopi/64) */
			c64_5 = 0.88192126434835502971, s64_5 = 0.47139673682599764856,	/* exp(5*I*twopi/64) */
			c64_6 = 0.83146961230254523708, s64_6 = 0.55557023301960222474,	/* exp(6*I*twopi/64) */
			c64_7 = 0.77301045336273696081, s64_7 = 0.63439328416364549822;	/* exp(7*I*twopi/64) */
  double t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F
		,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F
		,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F
		,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t4A,t4B,t4C,t4D,t4E,t4F
		,t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t5A,t5B,t5C,t5D,t5E,t5F
		,t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t6A,t6B,t6C,t6D,t6E,t6F
		,t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t7A,t7B,t7C,t7D,t7E,t7F;

	if(!first_entry && (n >> 6) != n64)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n64=n/64;

		p01 = n64;
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
	}

/*...The radix-64 pass is here.	*/

	for(j = 0; j < n64; j += 2)
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
			1,0,1,0,1,0,1,0,1,0,1,0,1,0
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
			0,1,ISRT2,ISRT2,-ISRT2,ISRT2,c64_4,s64_4,-s64_4,c64_4,s64_4,c64_4,-c64_4,s64_4
		);
		/* Block 2: t*4,5, twiddles = ~{   1,  E^ 2,  E^ 4,  E^ 6,  E^ 8, *E^ 6, *E^ 4, *E^ 2}
							BR order: ~{   1,  E^ 8,  E^ 4, *E^ 4,  E^ 2, *E^ 6,  E^ 6, *E^ 2}
			In terms of [Re,Im] pairs:~{[1,0],isrt2*[1,1],[c4,s4],[s4,c4],[c2,s2],[s6,c6],[c6,s6],[s2,c2]}
		*/
		jt = j1 + p02;	jp = j2 + p02;
		RADIX_08_DIT_TWIDDLE_OOP(
			t04,t05,t14,t15,t24,t25,t34,t35,t44,t45,t54,t55,t64,t65,t74,t75,
			a[jt],a[jp],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt+p30],a[jp+p30],a[jt+p08],a[jp+p08],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p38],a[jp+p38],
			ISRT2,ISRT2,c64_4,s64_4,s64_4,c64_4,c64_2,s64_2,s64_6,c64_6,c64_6,s64_6,s64_2,c64_2
		);
		/* Block 6: t*C,D, twiddles = ~{   1,  E^ 6, *E^ 4,*~E^ 2,*~E^ 8,-~E^ 2, -E^ 4,-*E^ 6}
							BR order: ~{   1,*~E^ 8, *E^ 4, -E^ 4,  E^ 6,-~E^ 2,*~E^ 2,-*E^ 6}
			In terms of [Re,Im] pairs:~{[1,0],isrt2*[-1,1],[s4,c4],[-c4,-s4],[c6,s6],[-c2,s2],[-s2,c2],[-s6,-c6]}
		*/
		jt = j1 + p06;	jp = j2 + p06;
		RADIX_08_DIT_TWIDDLE_OOP(
			t0C,t0D,t1C,t1D,t2C,t2D,t3C,t3D,t4C,t4D,t5C,t5D,t6C,t6D,t7C,t7D,
			a[jt],a[jp],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt+p30],a[jp+p30],a[jt+p08],a[jp+p08],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p38],a[jp+p38],
			-ISRT2,ISRT2,s64_4,c64_4,-c64_4,-s64_4,c64_6,s64_6,-c64_2,s64_2,-s64_2,c64_2,-s64_6,-c64_6
		);
		/* Block 1: t*2,3, twiddles = ~{   1,  E^ 1,  E^ 2,  E^ 3,  E^ 4,  E^ 5,  E^ 6,  E^ 7}
							BR order: ~{   1,  E^ 4,  E^ 2,  E^ 6,  E^ 1,  E^ 5,  E^ 3,  E^ 7}
			In terms of [Re,Im] pairs:~{[1,0],[c4,s4],[c2,s2],[c6,s6],[c1,s1],[c5,s5],[c3,s3],[c7,s7]}
		*/
		jt = j1 + p01;	jp = j2 + p01;
		RADIX_08_DIT_TWIDDLE_OOP(
			t02,t03,t12,t13,t22,t23,t32,t33,t42,t43,t52,t53,t62,t63,t72,t73,
			a[jt],a[jp],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt+p30],a[jp+p30],a[jt+p08],a[jp+p08],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p38],a[jp+p38],
			c64_4,s64_4,c64_2,s64_2,c64_6,s64_6,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7
		);
		/* Block 5: t*A,B, twiddles = ~{   1,  E^ 5, *E^ 6, *E^ 1,*~E^ 4,-~E^ 7,-~E^ 2, -E^ 3}
							BR order: ~{   1,*~E^ 4, *E^ 6,-~E^ 2,  E^ 5,-~E^ 7, *E^ 1, -E^ 3}
			In terms of [Re,Im] pairs:~{[1,0],[-s4,c4],[s6,c6],[-c2,s2],[c5,s5],[-c7,s7],[s1,c1],[-c3,-s3]}
		*/
		jt = j1 + p05;	jp = j2 + p05;
		RADIX_08_DIT_TWIDDLE_OOP(
			t0A,t0B,t1A,t1B,t2A,t2B,t3A,t3B,t4A,t4B,t5A,t5B,t6A,t6B,t7A,t7B,
			a[jt],a[jp],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt+p30],a[jp+p30],a[jt+p08],a[jp+p08],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p38],a[jp+p38],
			-s64_4,c64_4,s64_6,c64_6,-c64_2,s64_2,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3
		);
		/* Block 3: t*6,7, twiddles = ~{   1,  E^ 3,  E^ 6, *E^ 7, *E^ 4, *E^ 1,*~E^ 2,*~E^ 5}
							BR order: ~{   1, *E^ 4,  E^ 6,*~E^ 2,  E^ 3, *E^ 1, *E^ 7,*~E^ 5}
			In terms of [Re,Im] pairs:~{[1,0],[s4,c4],[c6,s6],[-s2,c2],[c3,s3],[s1,c1],[s7,c7],[-s5,c5]}
		*/
		jt = j1 + p03;	jp = j2 + p03;
		RADIX_08_DIT_TWIDDLE_OOP(
			t06,t07,t16,t17,t26,t27,t36,t37,t46,t47,t56,t57,t66,t67,t76,t77,
			a[jt],a[jp],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt+p30],a[jp+p30],a[jt+p08],a[jp+p08],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p38],a[jp+p38],
			s64_4,c64_4,c64_6,s64_6,-s64_2,c64_2,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5
		);
		/* Block 7: t*E,F, twiddles = ~{   1,  E^ 7, *E^ 2,*~E^ 5,-~E^ 4, -E^ 3,-*E^ 6,~*E^ 1} ,
							BR order: ~{   1,-~E^ 4, *E^ 2,-*E^ 6,  E^ 7, -E^ 3,*~E^ 5,~*E^ 1} ,
			In terms of [Re,Im] pairs:~{[1,0],[-c4,s4],[s2,c2],[-s6,-c6],[c7,s7],[-c3,-s3],[-s5,c5],[s1,-c1]}
		*/
		jt = j1 + p07;	jp = j2 + p07;
		RADIX_08_DIT_TWIDDLE_OOP(
			t0E,t0F,t1E,t1F,t2E,t2F,t3E,t3F,t4E,t4F,t5E,t5F,t6E,t6F,t7E,t7F,
			a[jt],a[jp],a[jt+p20],a[jp+p20],a[jt+p10],a[jp+p10],a[jt+p30],a[jp+p30],a[jt+p08],a[jp+p08],a[jt+p28],a[jp+p28],a[jt+p18],a[jp+p18],a[jt+p38],a[jp+p38],
			-c64_4,s64_4,s64_2,c64_2,-s64_6,-c64_6,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1
		);

	}
}

/******************** Multithreaded function body: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void* 
	cy64_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		const uint32 RADIX = 64;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30,p38;
		int j,j1,j2,k,l;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
	#ifdef USE_AVX
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#else
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	#endif
		double rt,it,wt_re,wt_im;	/* Fermat-mod weights stuff, used in both scalar and AVX mode */
		int k1,k2;

	#ifdef USE_SSE2

		int idx_offset,idx_incr;
		const double crnd = 3.0*0x4000000*0x2000000;
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7;	/* Addresses into array sections */
		vec_dbl *max_err, *sse2_rnd, *half_arr, *cc0, *ss0,
		 *isrt2, *cc1, *ss1, *cc2, *ss2, *cc3, *ss3, *cc4, *ss4, *cc5, *ss5, *cc6, *ss6, *cc7, *ss7,
		*nisrt2,*ncc1,*nss1,*ncc2,*nss2,*ncc3,*nss3,*ncc4,*nss4,*ncc5,*nss5,*ncc6,*nss6,*ncc7,*nss7,	// each non-unity root now needs a negated counterpart
		*r00,*r02,*r04,*r06,*r08,*r0A,*r0C,*r0E,
		*r10,*r12,*r14,*r16,*r18,*r1A,*r1C,*r1E,
		*r20,*r22,*r24,*r26,*r28,*r2A,*r2C,*r2E,
		*r30,*r32,*r34,*r36,*r38,*r3A,*r3C,*r3E,
		*r40,*r42,*r44,*r46,*r48,*r4A,*r4C,*r4E,
		*r50,*r52,*r54,*r56,*r58,*r5A,*r5C,*r5E,
		*r60,*r62,*r64,*r66,*r68,*r6A,*r6C,*r6E,
		*r70,*r72,*r74,*r76,*r78,*r7A,*r7C,*r7E,
		*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p0ar,*s1p0br,*s1p0cr,*s1p0dr,*s1p0er,*s1p0fr,
		*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p1ar,*s1p1br,*s1p1cr,*s1p1dr,*s1p1er,*s1p1fr,
		*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r,*s1p28r,*s1p29r,*s1p2ar,*s1p2br,*s1p2cr,*s1p2dr,*s1p2er,*s1p2fr,
		*s1p30r,*s1p31r,*s1p32r,*s1p33r,*s1p34r,*s1p35r,*s1p36r,*s1p37r,*s1p38r,*s1p39r,*s1p3ar,*s1p3br,*s1p3cr,*s1p3dr,*s1p3er,*s1p3fr,
		*cy_r00,*cy_r04,*cy_r08,*cy_r0C,*cy_r10,*cy_r14,*cy_r18,*cy_r1C,*cy_r20,*cy_r24,*cy_r28,*cy_r2C,*cy_r30,*cy_r34,*cy_r38,*cy_r3C,
		*cy_i00,*cy_i04,*cy_i08,*cy_i0C,*cy_i10,*cy_i14,*cy_i18,*cy_i1C,*cy_i20,*cy_i24,*cy_i28,*cy_i2C,*cy_i30,*cy_i34,*cy_i38,*cy_i3C;
	  #ifndef USE_AVX
		vec_dbl
		*cy_r02,*cy_r06,*cy_r0A,*cy_r0E,*cy_r12,*cy_r16,*cy_r1A,*cy_r1E,*cy_r22,*cy_r26,*cy_r2A,*cy_r2E,*cy_r32,*cy_r36,*cy_r3A,*cy_r3E,
		*cy_i02,*cy_i06,*cy_i0A,*cy_i0E,*cy_i12,*cy_i16,*cy_i1A,*cy_i1E,*cy_i22,*cy_i26,*cy_i2A,*cy_i2E,*cy_i32,*cy_i36,*cy_i3A,*cy_i3E;
	  #else
		vec_dbl *base_negacyclic_root;
	  #endif

		vec_dbl *tmp,*tm2;	// Non-static utility ptrs

		int
		 *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn0A,*bjmodn0B,*bjmodn0C,*bjmodn0D,*bjmodn0E,*bjmodn0F
		,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,*bjmodn1A,*bjmodn1B,*bjmodn1C,*bjmodn1D,*bjmodn1E,*bjmodn1F
		,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27,*bjmodn28,*bjmodn29,*bjmodn2A,*bjmodn2B,*bjmodn2C,*bjmodn2D,*bjmodn2E,*bjmodn2F
		,*bjmodn30,*bjmodn31,*bjmodn32,*bjmodn33,*bjmodn34,*bjmodn35,*bjmodn36,*bjmodn37,*bjmodn38,*bjmodn39,*bjmodn3A,*bjmodn3B,*bjmodn3C,*bjmodn3D,*bjmodn3E,*bjmodn3F;

		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_nm1;

	#else

		const double c64_1= 0.99518472667219688624, s64_1 = 0.09801714032956060199,	/* exp(  I*twopi/64) */
					c64_2 = 0.98078528040323044913, s64_2 = 0.19509032201612826785,	/* exp(2*I*twopi/64) */
					c64_3 = 0.95694033573220886494, s64_3 = 0.29028467725446236764,	/* exp(3*I*twopi/64) */
					c64_4 = 0.92387953251128675613, s64_4 = 0.38268343236508977173,	/* exp(4*I*twopi/64) */
					c64_5 = 0.88192126434835502971, s64_5 = 0.47139673682599764856,	/* exp(5*I*twopi/64) */
					c64_6 = 0.83146961230254523708, s64_6 = 0.55557023301960222474,	/* exp(6*I*twopi/64) */
					c64_7 = 0.77301045336273696081, s64_7 = 0.63439328416364549822;	/* exp(7*I*twopi/64) */
		double *base, *baseinv;
		const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int jt,jp,m,m2;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		double temp,frac,
			ar_p00,ar_p01,ar_p02,ar_p03,ar_p04,ar_p05,ar_p06,ar_p07,ar_p08,ar_p09,ar_p0a,ar_p0b,ar_p0c,ar_p0d,ar_p0e,ar_p0f,
			ar_p10,ar_p11,ar_p12,ar_p13,ar_p14,ar_p15,ar_p16,ar_p17,ar_p18,ar_p19,ar_p1a,ar_p1b,ar_p1c,ar_p1d,ar_p1e,ar_p1f,
			ar_p20,ar_p21,ar_p22,ar_p23,ar_p24,ar_p25,ar_p26,ar_p27,ar_p28,ar_p29,ar_p2a,ar_p2b,ar_p2c,ar_p2d,ar_p2e,ar_p2f,
			ar_p30,ar_p31,ar_p32,ar_p33,ar_p34,ar_p35,ar_p36,ar_p37,ar_p38,ar_p39,ar_p3a,ar_p3b,ar_p3c,ar_p3d,ar_p3e,ar_p3f,
			ai_p00,ai_p01,ai_p02,ai_p03,ai_p04,ai_p05,ai_p06,ai_p07,ai_p08,ai_p09,ai_p0a,ai_p0b,ai_p0c,ai_p0d,ai_p0e,ai_p0f,
			ai_p10,ai_p11,ai_p12,ai_p13,ai_p14,ai_p15,ai_p16,ai_p17,ai_p18,ai_p19,ai_p1a,ai_p1b,ai_p1c,ai_p1d,ai_p1e,ai_p1f,
			ai_p20,ai_p21,ai_p22,ai_p23,ai_p24,ai_p25,ai_p26,ai_p27,ai_p28,ai_p29,ai_p2a,ai_p2b,ai_p2c,ai_p2d,ai_p2e,ai_p2f,
			ai_p30,ai_p31,ai_p32,ai_p33,ai_p34,ai_p35,ai_p36,ai_p37,ai_p38,ai_p39,ai_p3a,ai_p3b,ai_p3c,ai_p3d,ai_p3e,ai_p3f,
			cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r0A,cy_r0B,cy_r0C,cy_r0D,cy_r0E,cy_r0F,
			cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r1A,cy_r1B,cy_r1C,cy_r1D,cy_r1E,cy_r1F,
			cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27,cy_r28,cy_r29,cy_r2A,cy_r2B,cy_r2C,cy_r2D,cy_r2E,cy_r2F,
			cy_r30,cy_r31,cy_r32,cy_r33,cy_r34,cy_r35,cy_r36,cy_r37,cy_r38,cy_r39,cy_r3A,cy_r3B,cy_r3C,cy_r3D,cy_r3E,cy_r3F,
			cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i0A,cy_i0B,cy_i0C,cy_i0D,cy_i0E,cy_i0F,
			cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i1A,cy_i1B,cy_i1C,cy_i1D,cy_i1E,cy_i1F,
			cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27,cy_i28,cy_i29,cy_i2A,cy_i2B,cy_i2C,cy_i2D,cy_i2E,cy_i2F,
			cy_i30,cy_i31,cy_i32,cy_i33,cy_i34,cy_i35,cy_i36,cy_i37,cy_i38,cy_i39,cy_i3A,cy_i3B,cy_i3C,cy_i3D,cy_i3E,cy_i3F,
			t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F,
			t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F,
			t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F,
			t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F,
			t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t4A,t4B,t4C,t4D,t4E,t4F,
			t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t5A,t5B,t5C,t5D,t5E,t5F,
			t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t6A,t6B,t6C,t6D,t6E,t6F,
			t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t7A,t7B,t7C,t7D,t7E,t7F;
		int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn0A,bjmodn0B,bjmodn0C,bjmodn0D,bjmodn0E,bjmodn0F
			,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn1A,bjmodn1B,bjmodn1C,bjmodn1D,bjmodn1E,bjmodn1F
			,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,bjmodn28,bjmodn29,bjmodn2A,bjmodn2B,bjmodn2C,bjmodn2D,bjmodn2E,bjmodn2F
			,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,bjmodn36,bjmodn37,bjmodn38,bjmodn39,bjmodn3A,bjmodn3B,bjmodn3C,bjmodn3D,bjmodn3E,bjmodn3F;

	#endif

		struct cy_thread_data_t* thread_arg = targ;

	// int data:
		int NDIVR = thread_arg->ndivr;
		int n = NDIVR*RADIX, nm1 = n-1;
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
		double scale = thread_arg->scale;

	// pointer data:
		double *a = thread_arg->arrdat;
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
		int *si = thread_arg->si;
		struct complex *rn0 = thread_arg->rn0;
		struct complex *rn1 = thread_arg->rn1;

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

		ASSERT(HERE, p01+p01 == p02, "p01+p01 != p02");
		ASSERT(HERE, p02+p02 == p04, "p02+p02 != p04");
		ASSERT(HERE, p04+p04 == p08, "p04+p04 != p08");
		ASSERT(HERE, p08+p08 == p10, "p08+p08 != p10");
		ASSERT(HERE, p10+p10 == p20, "p10+p10 != p20");
		ASSERT(HERE, p20+p18 == p38, "p20+p18 != p38");

	#ifdef USE_SSE2
		r00 = thread_arg->r00;	// declared above
		tmp = r00;
		r00 = tmp + 0x00;		r40 = tmp + 0x40;
		r02 = tmp + 0x02;		r42 = tmp + 0x42;
		r04 = tmp + 0x04;		r44 = tmp + 0x44;
		r06 = tmp + 0x06;		r46 = tmp + 0x46;
		r08 = tmp + 0x08;		r48 = tmp + 0x48;
		r0A = tmp + 0x0a;		r4A = tmp + 0x4a;
		r0C = tmp + 0x0c;		r4C = tmp + 0x4c;
		r0E = tmp + 0x0e;		r4E = tmp + 0x4e;
		r10 = tmp + 0x10;		r50 = tmp + 0x50;
		r12 = tmp + 0x12;		r52 = tmp + 0x52;
		r14 = tmp + 0x14;		r54 = tmp + 0x54;
		r16 = tmp + 0x16;		r56 = tmp + 0x56;
		r18 = tmp + 0x18;		r58 = tmp + 0x58;
		r1A = tmp + 0x1a;		r5A = tmp + 0x5a;
		r1C = tmp + 0x1c;		r5C = tmp + 0x5c;
		r1E = tmp + 0x1e;		r5E = tmp + 0x5e;
		r20 = tmp + 0x20;		r60 = tmp + 0x60;
		r22 = tmp + 0x22;		r62 = tmp + 0x62;
		r24 = tmp + 0x24;		r64 = tmp + 0x64;
		r26 = tmp + 0x26;		r66 = tmp + 0x66;
		r28 = tmp + 0x28;		r68 = tmp + 0x68;
		r2A = tmp + 0x2a;		r6A = tmp + 0x6a;
		r2C = tmp + 0x2c;		r6C = tmp + 0x6c;
		r2E = tmp + 0x2e;		r6E = tmp + 0x6e;
		r30 = tmp + 0x30;		r70 = tmp + 0x70;
		r32 = tmp + 0x32;		r72 = tmp + 0x72;
		r34 = tmp + 0x34;		r74 = tmp + 0x74;
		r36 = tmp + 0x36;		r76 = tmp + 0x76;
		r38 = tmp + 0x38;		r78 = tmp + 0x78;
		r3A = tmp + 0x3a;		r7A = tmp + 0x7a;
		r3C = tmp + 0x3c;		r7C = tmp + 0x7c;
		r3E = tmp + 0x3e;		r7E = tmp + 0x7e;
		tmp += 0x80;
		s1p00r = tmp + 0x00;	s1p20r = tmp + 0x40;
		s1p01r = tmp + 0x02;	s1p21r = tmp + 0x42;
		s1p02r = tmp + 0x04;	s1p22r = tmp + 0x44;
		s1p03r = tmp + 0x06;	s1p23r = tmp + 0x46;
		s1p04r = tmp + 0x08;	s1p24r = tmp + 0x48;
		s1p05r = tmp + 0x0a;	s1p25r = tmp + 0x4a;
		s1p06r = tmp + 0x0c;	s1p26r = tmp + 0x4c;
		s1p07r = tmp + 0x0e;	s1p27r = tmp + 0x4e;
		s1p08r = tmp + 0x10;	s1p28r = tmp + 0x50;
		s1p09r = tmp + 0x12;	s1p29r = tmp + 0x52;
		s1p0ar = tmp + 0x14;	s1p2ar = tmp + 0x54;
		s1p0br = tmp + 0x16;	s1p2br = tmp + 0x56;
		s1p0cr = tmp + 0x18;	s1p2cr = tmp + 0x58;
		s1p0dr = tmp + 0x1a;	s1p2dr = tmp + 0x5a;
		s1p0er = tmp + 0x1c;	s1p2er = tmp + 0x5c;
		s1p0fr = tmp + 0x1e;	s1p2fr = tmp + 0x5e;
		s1p10r = tmp + 0x20;	s1p30r = tmp + 0x60;
		s1p11r = tmp + 0x22;	s1p31r = tmp + 0x62;
		s1p12r = tmp + 0x24;	s1p32r = tmp + 0x64;
		s1p13r = tmp + 0x26;	s1p33r = tmp + 0x66;
		s1p14r = tmp + 0x28;	s1p34r = tmp + 0x68;
		s1p15r = tmp + 0x2a;	s1p35r = tmp + 0x6a;
		s1p16r = tmp + 0x2c;	s1p36r = tmp + 0x6c;
		s1p17r = tmp + 0x2e;	s1p37r = tmp + 0x6e;
		s1p18r = tmp + 0x30;	s1p38r = tmp + 0x70;
		s1p19r = tmp + 0x32;	s1p39r = tmp + 0x72;
		s1p1ar = tmp + 0x34;	s1p3ar = tmp + 0x74;
		s1p1br = tmp + 0x36;	s1p3br = tmp + 0x76;
		s1p1cr = tmp + 0x38;	s1p3cr = tmp + 0x78;
		s1p1dr = tmp + 0x3a;	s1p3dr = tmp + 0x7a;
		s1p1er = tmp + 0x3c;	s1p3er = tmp + 0x7c;
		s1p1fr = tmp + 0x3e;	s1p3fr = tmp + 0x7e;
		tmp += 0x80;
		// Each non-unity root now needs a negated counterpart:
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
		tmp += 0x2e;
	  #ifdef USE_AVX
		cy_r00	= tmp + 0x00;
		cy_r04	= tmp + 0x01;
		cy_r08	= tmp + 0x02;
		cy_r0C	= tmp + 0x03;
		cy_r10	= tmp + 0x04;
		cy_r14	= tmp + 0x05;
		cy_r18	= tmp + 0x06;
		cy_r1C	= tmp + 0x07;
		cy_r20	= tmp + 0x08;
		cy_r24	= tmp + 0x09;
		cy_r28	= tmp + 0x0a;
		cy_r2C	= tmp + 0x0b;
		cy_r30	= tmp + 0x0c;
		cy_r34	= tmp + 0x0d;
		cy_r38	= tmp + 0x0e;
		cy_r3C	= tmp + 0x0f;
		cy_i00	= tmp + 0x10;
		cy_i04	= tmp + 0x11;
		cy_i08	= tmp + 0x12;
		cy_i0C	= tmp + 0x13;
		cy_i10	= tmp + 0x14;
		cy_i14	= tmp + 0x15;
		cy_i18	= tmp + 0x16;
		cy_i1C	= tmp + 0x17;
		cy_i20	= tmp + 0x18;
		cy_i24	= tmp + 0x19;
		cy_i28	= tmp + 0x1a;
		cy_i2C	= tmp + 0x1b;
		cy_i30	= tmp + 0x1c;
		cy_i34	= tmp + 0x1d;
		cy_i38	= tmp + 0x1e;
		cy_i3C	= tmp + 0x1f;
		tmp += 0x20;
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 336 vec_dbl
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 68 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy_r00	= tmp + 0x00;		cy_i00	= tmp + 0x20;
		cy_r02	= tmp + 0x01;		cy_i02	= tmp + 0x21;
		cy_r04	= tmp + 0x02;		cy_i04	= tmp + 0x22;
		cy_r06	= tmp + 0x03;		cy_i06	= tmp + 0x23;
		cy_r08	= tmp + 0x04;		cy_i08	= tmp + 0x24;
		cy_r0A	= tmp + 0x05;		cy_i0A	= tmp + 0x25;
		cy_r0C	= tmp + 0x06;		cy_i0C	= tmp + 0x26;
		cy_r0E	= tmp + 0x07;		cy_i0E	= tmp + 0x27;
		cy_r10	= tmp + 0x08;		cy_i10	= tmp + 0x28;
		cy_r12	= tmp + 0x09;		cy_i12	= tmp + 0x29;
		cy_r14	= tmp + 0x0a;		cy_i14	= tmp + 0x2a;
		cy_r16	= tmp + 0x0b;		cy_i16	= tmp + 0x2b;
		cy_r18	= tmp + 0x0c;		cy_i18	= tmp + 0x2c;
		cy_r1A	= tmp + 0x0d;		cy_i1A	= tmp + 0x2d;
		cy_r1C	= tmp + 0x0e;		cy_i1C	= tmp + 0x2e;
		cy_r1E	= tmp + 0x0f;		cy_i1E	= tmp + 0x2f;
		cy_r20	= tmp + 0x10;		cy_i20	= tmp + 0x30;
		cy_r22	= tmp + 0x11;		cy_i22	= tmp + 0x31;
		cy_r24	= tmp + 0x12;		cy_i24	= tmp + 0x32;
		cy_r26	= tmp + 0x13;		cy_i26	= tmp + 0x33;
		cy_r28	= tmp + 0x14;		cy_i28	= tmp + 0x34;
		cy_r2A	= tmp + 0x15;		cy_i2A	= tmp + 0x35;
		cy_r2C	= tmp + 0x16;		cy_i2C	= tmp + 0x36;
		cy_r2E	= tmp + 0x17;		cy_i2E	= tmp + 0x37;
		cy_r30	= tmp + 0x18;		cy_i30	= tmp + 0x38;
		cy_r32	= tmp + 0x19;		cy_i32	= tmp + 0x39;
		cy_r34	= tmp + 0x1a;		cy_i34	= tmp + 0x3a;
		cy_r36	= tmp + 0x1b;		cy_i36	= tmp + 0x3b;
		cy_r38	= tmp + 0x1c;		cy_i38	= tmp + 0x3c;
		cy_r3A	= tmp + 0x1d;		cy_i3A	= tmp + 0x3d;
		cy_r3C	= tmp + 0x1e;		cy_i3C	= tmp + 0x3e;
		cy_r3E	= tmp + 0x1f;		cy_i3E	= tmp + 0x3f;
		tmp += 0x40;
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 368 complex
		// This is where the value of half_arr_offset comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, 2 for Fermat-mod */
	  #endif

		ASSERT(HERE, (r00 == thread_arg->r00), "thread-local memcheck failed!");
		ASSERT(HERE, (half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
		ASSERT(HERE, (isrt2->d0 == ISRT2 && isrt2->d1 == ISRT2), "thread-local memcheck failed!");
		ASSERT(HERE, (sse2_rnd->d0 == crnd && sse2_rnd->d1 == crnd), "thread-local memcheck failed!");

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		tmp = half_arr;
	  #ifdef USE_AVX
		// Grab some elt of base-data [offset by +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif
	} else {
		dtmp = (half_arr)->d0 * (half_arr+1)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (half_arr)->d1 * (half_arr+1)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	}

		VEC_DBL_INIT(max_err, 0.0);

		sign_mask = (uint64*)(r00 + radix64_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (#doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		sse_sw  = sse_bw    + RE_IM_STRIDE;
		sse_nm1 = sse_sw    + RE_IM_STRIDE;
	  #ifdef USE_AVX
		n_minus_sil   = (struct uint32x4 *)sse_nm1 + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_nm1 + 2;
		sinwt         = (struct uint32x4 *)sse_nm1 + 3;
		sinwtm1       = (struct uint32x4 *)sse_nm1 + 4;

		bjmodn00 = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn00 = (int*)(sse_nm1 + RE_IM_STRIDE);
	  #endif
										bjmodn20 = bjmodn00 + 0x20;
		bjmodn01 = bjmodn00 + 0x01;		bjmodn21 = bjmodn00 + 0x21;
		bjmodn02 = bjmodn00 + 0x02;		bjmodn22 = bjmodn00 + 0x22;
		bjmodn03 = bjmodn00 + 0x03;		bjmodn23 = bjmodn00 + 0x23;
		bjmodn04 = bjmodn00 + 0x04;		bjmodn24 = bjmodn00 + 0x24;
		bjmodn05 = bjmodn00 + 0x05;		bjmodn25 = bjmodn00 + 0x25;
		bjmodn06 = bjmodn00 + 0x06;		bjmodn26 = bjmodn00 + 0x26;
		bjmodn07 = bjmodn00 + 0x07;		bjmodn27 = bjmodn00 + 0x27;
		bjmodn08 = bjmodn00 + 0x08;		bjmodn28 = bjmodn00 + 0x28;
		bjmodn09 = bjmodn00 + 0x09;		bjmodn29 = bjmodn00 + 0x29;
		bjmodn0A = bjmodn00 + 0x0A;		bjmodn2A = bjmodn00 + 0x2A;
		bjmodn0B = bjmodn00 + 0x0B;		bjmodn2B = bjmodn00 + 0x2B;
		bjmodn0C = bjmodn00 + 0x0C;		bjmodn2C = bjmodn00 + 0x2C;
		bjmodn0D = bjmodn00 + 0x0D;		bjmodn2D = bjmodn00 + 0x2D;
		bjmodn0E = bjmodn00 + 0x0E;		bjmodn2E = bjmodn00 + 0x2E;
		bjmodn0F = bjmodn00 + 0x0F;		bjmodn2F = bjmodn00 + 0x2F;
		bjmodn10 = bjmodn00 + 0x10;		bjmodn30 = bjmodn00 + 0x30;
		bjmodn11 = bjmodn00 + 0x11;		bjmodn31 = bjmodn00 + 0x31;
		bjmodn12 = bjmodn00 + 0x12;		bjmodn32 = bjmodn00 + 0x32;
		bjmodn13 = bjmodn00 + 0x13;		bjmodn33 = bjmodn00 + 0x33;
		bjmodn14 = bjmodn00 + 0x14;		bjmodn34 = bjmodn00 + 0x34;
		bjmodn15 = bjmodn00 + 0x15;		bjmodn35 = bjmodn00 + 0x35;
		bjmodn16 = bjmodn00 + 0x16;		bjmodn36 = bjmodn00 + 0x36;
		bjmodn17 = bjmodn00 + 0x17;		bjmodn37 = bjmodn00 + 0x37;
		bjmodn18 = bjmodn00 + 0x18;		bjmodn38 = bjmodn00 + 0x38;
		bjmodn19 = bjmodn00 + 0x19;		bjmodn39 = bjmodn00 + 0x39;
		bjmodn1A = bjmodn00 + 0x1A;		bjmodn3A = bjmodn00 + 0x3A;
		bjmodn1B = bjmodn00 + 0x1B;		bjmodn3B = bjmodn00 + 0x3B;
		bjmodn1C = bjmodn00 + 0x1C;		bjmodn3C = bjmodn00 + 0x3C;
		bjmodn1D = bjmodn00 + 0x1D;		bjmodn3D = bjmodn00 + 0x3D;
		bjmodn1E = bjmodn00 + 0x1E;		bjmodn3E = bjmodn00 + 0x3E;
		bjmodn1F = bjmodn00 + 0x1F;		bjmodn3F = bjmodn00 + 0x3F;
	#else

		// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
		base    = (double *)thread_arg->r00     ;
		baseinv = (double *)thread_arg->half_arr;

	#endif	// USE_SSE2 ?

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Init DWT-indices: */				/* init carries	*/
		#ifdef USE_AVX
			*bjmodn00 = thread_arg->bjmodn00;	cy_r00->d0 = thread_arg->cy_r00;
			*bjmodn01 = thread_arg->bjmodn01;	cy_r00->d1 = thread_arg->cy_r01;
			*bjmodn02 = thread_arg->bjmodn02;	cy_r00->d2 = thread_arg->cy_r02;
			*bjmodn03 = thread_arg->bjmodn03;	cy_r00->d3 = thread_arg->cy_r03;
			*bjmodn04 = thread_arg->bjmodn04;	cy_r04->d0 = thread_arg->cy_r04;
			*bjmodn05 = thread_arg->bjmodn05;	cy_r04->d1 = thread_arg->cy_r05;
			*bjmodn06 = thread_arg->bjmodn06;	cy_r04->d2 = thread_arg->cy_r06;
			*bjmodn07 = thread_arg->bjmodn07;	cy_r04->d3 = thread_arg->cy_r07;
			*bjmodn08 = thread_arg->bjmodn08;	cy_r08->d0 = thread_arg->cy_r08;
			*bjmodn09 = thread_arg->bjmodn09;	cy_r08->d1 = thread_arg->cy_r09;
			*bjmodn0A = thread_arg->bjmodn0A;	cy_r08->d2 = thread_arg->cy_r0A;
			*bjmodn0B = thread_arg->bjmodn0B;	cy_r08->d3 = thread_arg->cy_r0B;
			*bjmodn0C = thread_arg->bjmodn0C;	cy_r0C->d0 = thread_arg->cy_r0C;
			*bjmodn0D = thread_arg->bjmodn0D;	cy_r0C->d1 = thread_arg->cy_r0D;
			*bjmodn0E = thread_arg->bjmodn0E;	cy_r0C->d2 = thread_arg->cy_r0E;
			*bjmodn0F = thread_arg->bjmodn0F;	cy_r0C->d3 = thread_arg->cy_r0F;
			*bjmodn10 = thread_arg->bjmodn10;	cy_r10->d0 = thread_arg->cy_r10;
			*bjmodn11 = thread_arg->bjmodn11;	cy_r10->d1 = thread_arg->cy_r11;
			*bjmodn12 = thread_arg->bjmodn12;	cy_r10->d2 = thread_arg->cy_r12;
			*bjmodn13 = thread_arg->bjmodn13;	cy_r10->d3 = thread_arg->cy_r13;
			*bjmodn14 = thread_arg->bjmodn14;	cy_r14->d0 = thread_arg->cy_r14;
			*bjmodn15 = thread_arg->bjmodn15;	cy_r14->d1 = thread_arg->cy_r15;
			*bjmodn16 = thread_arg->bjmodn16;	cy_r14->d2 = thread_arg->cy_r16;
			*bjmodn17 = thread_arg->bjmodn17;	cy_r14->d3 = thread_arg->cy_r17;
			*bjmodn18 = thread_arg->bjmodn18;	cy_r18->d0 = thread_arg->cy_r18;
			*bjmodn19 = thread_arg->bjmodn19;	cy_r18->d1 = thread_arg->cy_r19;
			*bjmodn1A = thread_arg->bjmodn1A;	cy_r18->d2 = thread_arg->cy_r1A;
			*bjmodn1B = thread_arg->bjmodn1B;	cy_r18->d3 = thread_arg->cy_r1B;
			*bjmodn1C = thread_arg->bjmodn1C;	cy_r1C->d0 = thread_arg->cy_r1C;
			*bjmodn1D = thread_arg->bjmodn1D;	cy_r1C->d1 = thread_arg->cy_r1D;
			*bjmodn1E = thread_arg->bjmodn1E;	cy_r1C->d2 = thread_arg->cy_r1E;
			*bjmodn1F = thread_arg->bjmodn1F;	cy_r1C->d3 = thread_arg->cy_r1F;
			*bjmodn20 = thread_arg->bjmodn20;	cy_r20->d0 = thread_arg->cy_r20;
			*bjmodn21 = thread_arg->bjmodn21;	cy_r20->d1 = thread_arg->cy_r21;
			*bjmodn22 = thread_arg->bjmodn22;	cy_r20->d2 = thread_arg->cy_r22;
			*bjmodn23 = thread_arg->bjmodn23;	cy_r20->d3 = thread_arg->cy_r23;
			*bjmodn24 = thread_arg->bjmodn24;	cy_r24->d0 = thread_arg->cy_r24;
			*bjmodn25 = thread_arg->bjmodn25;	cy_r24->d1 = thread_arg->cy_r25;
			*bjmodn26 = thread_arg->bjmodn26;	cy_r24->d2 = thread_arg->cy_r26;
			*bjmodn27 = thread_arg->bjmodn27;	cy_r24->d3 = thread_arg->cy_r27;
			*bjmodn28 = thread_arg->bjmodn28;	cy_r28->d0 = thread_arg->cy_r28;
			*bjmodn29 = thread_arg->bjmodn29;	cy_r28->d1 = thread_arg->cy_r29;
			*bjmodn2A = thread_arg->bjmodn2A;	cy_r28->d2 = thread_arg->cy_r2A;
			*bjmodn2B = thread_arg->bjmodn2B;	cy_r28->d3 = thread_arg->cy_r2B;
			*bjmodn2C = thread_arg->bjmodn2C;	cy_r2C->d0 = thread_arg->cy_r2C;
			*bjmodn2D = thread_arg->bjmodn2D;	cy_r2C->d1 = thread_arg->cy_r2D;
			*bjmodn2E = thread_arg->bjmodn2E;	cy_r2C->d2 = thread_arg->cy_r2E;
			*bjmodn2F = thread_arg->bjmodn2F;	cy_r2C->d3 = thread_arg->cy_r2F;
			*bjmodn30 = thread_arg->bjmodn30;	cy_r30->d0 = thread_arg->cy_r30;
			*bjmodn31 = thread_arg->bjmodn31;	cy_r30->d1 = thread_arg->cy_r31;
			*bjmodn32 = thread_arg->bjmodn32;	cy_r30->d2 = thread_arg->cy_r32;
			*bjmodn33 = thread_arg->bjmodn33;	cy_r30->d3 = thread_arg->cy_r33;
			*bjmodn34 = thread_arg->bjmodn34;	cy_r34->d0 = thread_arg->cy_r34;
			*bjmodn35 = thread_arg->bjmodn35;	cy_r34->d1 = thread_arg->cy_r35;
			*bjmodn36 = thread_arg->bjmodn36;	cy_r34->d2 = thread_arg->cy_r36;
			*bjmodn37 = thread_arg->bjmodn37;	cy_r34->d3 = thread_arg->cy_r37;
			*bjmodn38 = thread_arg->bjmodn38;	cy_r38->d0 = thread_arg->cy_r38;
			*bjmodn39 = thread_arg->bjmodn39;	cy_r38->d1 = thread_arg->cy_r39;
			*bjmodn3A = thread_arg->bjmodn3A;	cy_r38->d2 = thread_arg->cy_r3A;
			*bjmodn3B = thread_arg->bjmodn3B;	cy_r38->d3 = thread_arg->cy_r3B;
			*bjmodn3C = thread_arg->bjmodn3C;	cy_r3C->d0 = thread_arg->cy_r3C;
			*bjmodn3D = thread_arg->bjmodn3D;	cy_r3C->d1 = thread_arg->cy_r3D;
			*bjmodn3E = thread_arg->bjmodn3E;	cy_r3C->d2 = thread_arg->cy_r3E;
			*bjmodn3F = thread_arg->bjmodn3F;	cy_r3C->d3 = thread_arg->cy_r3F;

		#elif defined(USE_SSE2)

			*bjmodn00 = thread_arg->bjmodn00;	cy_r00->d0 = thread_arg->cy_r00;
			*bjmodn01 = thread_arg->bjmodn01;	cy_r00->d1 = thread_arg->cy_r01;
			*bjmodn02 = thread_arg->bjmodn02;	cy_r02->d0 = thread_arg->cy_r02;
			*bjmodn03 = thread_arg->bjmodn03;	cy_r02->d1 = thread_arg->cy_r03;
			*bjmodn04 = thread_arg->bjmodn04;	cy_r04->d0 = thread_arg->cy_r04;
			*bjmodn05 = thread_arg->bjmodn05;	cy_r04->d1 = thread_arg->cy_r05;
			*bjmodn06 = thread_arg->bjmodn06;	cy_r06->d0 = thread_arg->cy_r06;
			*bjmodn07 = thread_arg->bjmodn07;	cy_r06->d1 = thread_arg->cy_r07;
			*bjmodn08 = thread_arg->bjmodn08;	cy_r08->d0 = thread_arg->cy_r08;
			*bjmodn09 = thread_arg->bjmodn09;	cy_r08->d1 = thread_arg->cy_r09;
			*bjmodn0A = thread_arg->bjmodn0A;	cy_r0A->d0 = thread_arg->cy_r0A;
			*bjmodn0B = thread_arg->bjmodn0B;	cy_r0A->d1 = thread_arg->cy_r0B;
			*bjmodn0C = thread_arg->bjmodn0C;	cy_r0C->d0 = thread_arg->cy_r0C;
			*bjmodn0D = thread_arg->bjmodn0D;	cy_r0C->d1 = thread_arg->cy_r0D;
			*bjmodn0E = thread_arg->bjmodn0E;	cy_r0E->d0 = thread_arg->cy_r0E;
			*bjmodn0F = thread_arg->bjmodn0F;	cy_r0E->d1 = thread_arg->cy_r0F;
			*bjmodn10 = thread_arg->bjmodn10;	cy_r10->d0 = thread_arg->cy_r10;
			*bjmodn11 = thread_arg->bjmodn11;	cy_r10->d1 = thread_arg->cy_r11;
			*bjmodn12 = thread_arg->bjmodn12;	cy_r12->d0 = thread_arg->cy_r12;
			*bjmodn13 = thread_arg->bjmodn13;	cy_r12->d1 = thread_arg->cy_r13;
			*bjmodn14 = thread_arg->bjmodn14;	cy_r14->d0 = thread_arg->cy_r14;
			*bjmodn15 = thread_arg->bjmodn15;	cy_r14->d1 = thread_arg->cy_r15;
			*bjmodn16 = thread_arg->bjmodn16;	cy_r16->d0 = thread_arg->cy_r16;
			*bjmodn17 = thread_arg->bjmodn17;	cy_r16->d1 = thread_arg->cy_r17;
			*bjmodn18 = thread_arg->bjmodn18;	cy_r18->d0 = thread_arg->cy_r18;
			*bjmodn19 = thread_arg->bjmodn19;	cy_r18->d1 = thread_arg->cy_r19;
			*bjmodn1A = thread_arg->bjmodn1A;	cy_r1A->d0 = thread_arg->cy_r1A;
			*bjmodn1B = thread_arg->bjmodn1B;	cy_r1A->d1 = thread_arg->cy_r1B;
			*bjmodn1C = thread_arg->bjmodn1C;	cy_r1C->d0 = thread_arg->cy_r1C;
			*bjmodn1D = thread_arg->bjmodn1D;	cy_r1C->d1 = thread_arg->cy_r1D;
			*bjmodn1E = thread_arg->bjmodn1E;	cy_r1E->d0 = thread_arg->cy_r1E;
			*bjmodn1F = thread_arg->bjmodn1F;	cy_r1E->d1 = thread_arg->cy_r1F;
			*bjmodn20 = thread_arg->bjmodn20;	cy_r20->d0 = thread_arg->cy_r20;
			*bjmodn21 = thread_arg->bjmodn21;	cy_r20->d1 = thread_arg->cy_r21;
			*bjmodn22 = thread_arg->bjmodn22;	cy_r22->d0 = thread_arg->cy_r22;
			*bjmodn23 = thread_arg->bjmodn23;	cy_r22->d1 = thread_arg->cy_r23;
			*bjmodn24 = thread_arg->bjmodn24;	cy_r24->d0 = thread_arg->cy_r24;
			*bjmodn25 = thread_arg->bjmodn25;	cy_r24->d1 = thread_arg->cy_r25;
			*bjmodn26 = thread_arg->bjmodn26;	cy_r26->d0 = thread_arg->cy_r26;
			*bjmodn27 = thread_arg->bjmodn27;	cy_r26->d1 = thread_arg->cy_r27;
			*bjmodn28 = thread_arg->bjmodn28;	cy_r28->d0 = thread_arg->cy_r28;
			*bjmodn29 = thread_arg->bjmodn29;	cy_r28->d1 = thread_arg->cy_r29;
			*bjmodn2A = thread_arg->bjmodn2A;	cy_r2A->d0 = thread_arg->cy_r2A;
			*bjmodn2B = thread_arg->bjmodn2B;	cy_r2A->d1 = thread_arg->cy_r2B;
			*bjmodn2C = thread_arg->bjmodn2C;	cy_r2C->d0 = thread_arg->cy_r2C;
			*bjmodn2D = thread_arg->bjmodn2D;	cy_r2C->d1 = thread_arg->cy_r2D;
			*bjmodn2E = thread_arg->bjmodn2E;	cy_r2E->d0 = thread_arg->cy_r2E;
			*bjmodn2F = thread_arg->bjmodn2F;	cy_r2E->d1 = thread_arg->cy_r2F;
			*bjmodn30 = thread_arg->bjmodn30;	cy_r30->d0 = thread_arg->cy_r30;
			*bjmodn31 = thread_arg->bjmodn31;	cy_r30->d1 = thread_arg->cy_r31;
			*bjmodn32 = thread_arg->bjmodn32;	cy_r32->d0 = thread_arg->cy_r32;
			*bjmodn33 = thread_arg->bjmodn33;	cy_r32->d1 = thread_arg->cy_r33;
			*bjmodn34 = thread_arg->bjmodn34;	cy_r34->d0 = thread_arg->cy_r34;
			*bjmodn35 = thread_arg->bjmodn35;	cy_r34->d1 = thread_arg->cy_r35;
			*bjmodn36 = thread_arg->bjmodn36;	cy_r36->d0 = thread_arg->cy_r36;
			*bjmodn37 = thread_arg->bjmodn37;	cy_r36->d1 = thread_arg->cy_r37;
			*bjmodn38 = thread_arg->bjmodn38;	cy_r38->d0 = thread_arg->cy_r38;
			*bjmodn39 = thread_arg->bjmodn39;	cy_r38->d1 = thread_arg->cy_r39;
			*bjmodn3A = thread_arg->bjmodn3A;	cy_r3A->d0 = thread_arg->cy_r3A;
			*bjmodn3B = thread_arg->bjmodn3B;	cy_r3A->d1 = thread_arg->cy_r3B;
			*bjmodn3C = thread_arg->bjmodn3C;	cy_r3C->d0 = thread_arg->cy_r3C;
			*bjmodn3D = thread_arg->bjmodn3D;	cy_r3C->d1 = thread_arg->cy_r3D;
			*bjmodn3E = thread_arg->bjmodn3E;	cy_r3E->d0 = thread_arg->cy_r3E;
			*bjmodn3F = thread_arg->bjmodn3F;	cy_r3E->d1 = thread_arg->cy_r3F;

		#else

			bjmodn00 = thread_arg->bjmodn00;	cy_r00 = thread_arg->cy_r00;
			bjmodn01 = thread_arg->bjmodn01;	cy_r01 = thread_arg->cy_r01;
			bjmodn02 = thread_arg->bjmodn02;	cy_r02 = thread_arg->cy_r02;
			bjmodn03 = thread_arg->bjmodn03;	cy_r03 = thread_arg->cy_r03;
			bjmodn04 = thread_arg->bjmodn04;	cy_r04 = thread_arg->cy_r04;
			bjmodn05 = thread_arg->bjmodn05;	cy_r05 = thread_arg->cy_r05;
			bjmodn06 = thread_arg->bjmodn06;	cy_r06 = thread_arg->cy_r06;
			bjmodn07 = thread_arg->bjmodn07;	cy_r07 = thread_arg->cy_r07;
			bjmodn08 = thread_arg->bjmodn08;	cy_r08 = thread_arg->cy_r08;
			bjmodn09 = thread_arg->bjmodn09;	cy_r09 = thread_arg->cy_r09;
			bjmodn0A = thread_arg->bjmodn0A;	cy_r0A = thread_arg->cy_r0A;
			bjmodn0B = thread_arg->bjmodn0B;	cy_r0B = thread_arg->cy_r0B;
			bjmodn0C = thread_arg->bjmodn0C;	cy_r0C = thread_arg->cy_r0C;
			bjmodn0D = thread_arg->bjmodn0D;	cy_r0D = thread_arg->cy_r0D;
			bjmodn0E = thread_arg->bjmodn0E;	cy_r0E = thread_arg->cy_r0E;
			bjmodn0F = thread_arg->bjmodn0F;	cy_r0F = thread_arg->cy_r0F;
			bjmodn10 = thread_arg->bjmodn10;	cy_r10 = thread_arg->cy_r10;
			bjmodn11 = thread_arg->bjmodn11;	cy_r11 = thread_arg->cy_r11;
			bjmodn12 = thread_arg->bjmodn12;	cy_r12 = thread_arg->cy_r12;
			bjmodn13 = thread_arg->bjmodn13;	cy_r13 = thread_arg->cy_r13;
			bjmodn14 = thread_arg->bjmodn14;	cy_r14 = thread_arg->cy_r14;
			bjmodn15 = thread_arg->bjmodn15;	cy_r15 = thread_arg->cy_r15;
			bjmodn16 = thread_arg->bjmodn16;	cy_r16 = thread_arg->cy_r16;
			bjmodn17 = thread_arg->bjmodn17;	cy_r17 = thread_arg->cy_r17;
			bjmodn18 = thread_arg->bjmodn18;	cy_r18 = thread_arg->cy_r18;
			bjmodn19 = thread_arg->bjmodn19;	cy_r19 = thread_arg->cy_r19;
			bjmodn1A = thread_arg->bjmodn1A;	cy_r1A = thread_arg->cy_r1A;
			bjmodn1B = thread_arg->bjmodn1B;	cy_r1B = thread_arg->cy_r1B;
			bjmodn1C = thread_arg->bjmodn1C;	cy_r1C = thread_arg->cy_r1C;
			bjmodn1D = thread_arg->bjmodn1D;	cy_r1D = thread_arg->cy_r1D;
			bjmodn1E = thread_arg->bjmodn1E;	cy_r1E = thread_arg->cy_r1E;
			bjmodn1F = thread_arg->bjmodn1F;	cy_r1F = thread_arg->cy_r1F;
			bjmodn20 = thread_arg->bjmodn20;	cy_r20 = thread_arg->cy_r20;
			bjmodn21 = thread_arg->bjmodn21;	cy_r21 = thread_arg->cy_r21;
			bjmodn22 = thread_arg->bjmodn22;	cy_r22 = thread_arg->cy_r22;
			bjmodn23 = thread_arg->bjmodn23;	cy_r23 = thread_arg->cy_r23;
			bjmodn24 = thread_arg->bjmodn24;	cy_r24 = thread_arg->cy_r24;
			bjmodn25 = thread_arg->bjmodn25;	cy_r25 = thread_arg->cy_r25;
			bjmodn26 = thread_arg->bjmodn26;	cy_r26 = thread_arg->cy_r26;
			bjmodn27 = thread_arg->bjmodn27;	cy_r27 = thread_arg->cy_r27;
			bjmodn28 = thread_arg->bjmodn28;	cy_r28 = thread_arg->cy_r28;
			bjmodn29 = thread_arg->bjmodn29;	cy_r29 = thread_arg->cy_r29;
			bjmodn2A = thread_arg->bjmodn2A;	cy_r2A = thread_arg->cy_r2A;
			bjmodn2B = thread_arg->bjmodn2B;	cy_r2B = thread_arg->cy_r2B;
			bjmodn2C = thread_arg->bjmodn2C;	cy_r2C = thread_arg->cy_r2C;
			bjmodn2D = thread_arg->bjmodn2D;	cy_r2D = thread_arg->cy_r2D;
			bjmodn2E = thread_arg->bjmodn2E;	cy_r2E = thread_arg->cy_r2E;
			bjmodn2F = thread_arg->bjmodn2F;	cy_r2F = thread_arg->cy_r2F;
			bjmodn30 = thread_arg->bjmodn30;	cy_r30 = thread_arg->cy_r30;
			bjmodn31 = thread_arg->bjmodn31;	cy_r31 = thread_arg->cy_r31;
			bjmodn32 = thread_arg->bjmodn32;	cy_r32 = thread_arg->cy_r32;
			bjmodn33 = thread_arg->bjmodn33;	cy_r33 = thread_arg->cy_r33;
			bjmodn34 = thread_arg->bjmodn34;	cy_r34 = thread_arg->cy_r34;
			bjmodn35 = thread_arg->bjmodn35;	cy_r35 = thread_arg->cy_r35;
			bjmodn36 = thread_arg->bjmodn36;	cy_r36 = thread_arg->cy_r36;
			bjmodn37 = thread_arg->bjmodn37;	cy_r37 = thread_arg->cy_r37;
			bjmodn38 = thread_arg->bjmodn38;	cy_r38 = thread_arg->cy_r38;
			bjmodn39 = thread_arg->bjmodn39;	cy_r39 = thread_arg->cy_r39;
			bjmodn3A = thread_arg->bjmodn3A;	cy_r3A = thread_arg->cy_r3A;
			bjmodn3B = thread_arg->bjmodn3B;	cy_r3B = thread_arg->cy_r3B;
			bjmodn3C = thread_arg->bjmodn3C;	cy_r3C = thread_arg->cy_r3C;
			bjmodn3D = thread_arg->bjmodn3D;	cy_r3D = thread_arg->cy_r3D;
			bjmodn3E = thread_arg->bjmodn3E;	cy_r3E = thread_arg->cy_r3E;
			bjmodn3F = thread_arg->bjmodn3F;	cy_r3F = thread_arg->cy_r3F;
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
			cy_r08->d2 = thread_arg->cy_r0A;	cy_i08->d2 = thread_arg->cy_i0A;
			cy_r08->d3 = thread_arg->cy_r0B;	cy_i08->d3 = thread_arg->cy_i0B;
			cy_r0C->d0 = thread_arg->cy_r0C;	cy_i0C->d0 = thread_arg->cy_i0C;
			cy_r0C->d1 = thread_arg->cy_r0D;	cy_i0C->d1 = thread_arg->cy_i0D;
			cy_r0C->d2 = thread_arg->cy_r0E;	cy_i0C->d2 = thread_arg->cy_i0E;
			cy_r0C->d3 = thread_arg->cy_r0F;	cy_i0C->d3 = thread_arg->cy_i0F;
			cy_r10->d0 = thread_arg->cy_r10;	cy_i10->d0 = thread_arg->cy_i10;
			cy_r10->d1 = thread_arg->cy_r11;	cy_i10->d1 = thread_arg->cy_i11;
			cy_r10->d2 = thread_arg->cy_r12;	cy_i10->d2 = thread_arg->cy_i12;
			cy_r10->d3 = thread_arg->cy_r13;	cy_i10->d3 = thread_arg->cy_i13;
			cy_r14->d0 = thread_arg->cy_r14;	cy_i14->d0 = thread_arg->cy_i14;
			cy_r14->d1 = thread_arg->cy_r15;	cy_i14->d1 = thread_arg->cy_i15;
			cy_r14->d2 = thread_arg->cy_r16;	cy_i14->d2 = thread_arg->cy_i16;
			cy_r14->d3 = thread_arg->cy_r17;	cy_i14->d3 = thread_arg->cy_i17;
			cy_r18->d0 = thread_arg->cy_r18;	cy_i18->d0 = thread_arg->cy_i18;
			cy_r18->d1 = thread_arg->cy_r19;	cy_i18->d1 = thread_arg->cy_i19;
			cy_r18->d2 = thread_arg->cy_r1A;	cy_i18->d2 = thread_arg->cy_i1A;
			cy_r18->d3 = thread_arg->cy_r1B;	cy_i18->d3 = thread_arg->cy_i1B;
			cy_r1C->d0 = thread_arg->cy_r1C;	cy_i1C->d0 = thread_arg->cy_i1C;
			cy_r1C->d1 = thread_arg->cy_r1D;	cy_i1C->d1 = thread_arg->cy_i1D;
			cy_r1C->d2 = thread_arg->cy_r1E;	cy_i1C->d2 = thread_arg->cy_i1E;
			cy_r1C->d3 = thread_arg->cy_r1F;	cy_i1C->d3 = thread_arg->cy_i1F;
			cy_r20->d0 = thread_arg->cy_r20;	cy_i20->d0 = thread_arg->cy_i20;
			cy_r20->d1 = thread_arg->cy_r21;	cy_i20->d1 = thread_arg->cy_i21;
			cy_r20->d2 = thread_arg->cy_r22;	cy_i20->d2 = thread_arg->cy_i22;
			cy_r20->d3 = thread_arg->cy_r23;	cy_i20->d3 = thread_arg->cy_i23;
			cy_r24->d0 = thread_arg->cy_r24;	cy_i24->d0 = thread_arg->cy_i24;
			cy_r24->d1 = thread_arg->cy_r25;	cy_i24->d1 = thread_arg->cy_i25;
			cy_r24->d2 = thread_arg->cy_r26;	cy_i24->d2 = thread_arg->cy_i26;
			cy_r24->d3 = thread_arg->cy_r27;	cy_i24->d3 = thread_arg->cy_i27;
			cy_r28->d0 = thread_arg->cy_r28;	cy_i28->d0 = thread_arg->cy_i28;
			cy_r28->d1 = thread_arg->cy_r29;	cy_i28->d1 = thread_arg->cy_i29;
			cy_r28->d2 = thread_arg->cy_r2A;	cy_i28->d2 = thread_arg->cy_i2A;
			cy_r28->d3 = thread_arg->cy_r2B;	cy_i28->d3 = thread_arg->cy_i2B;
			cy_r2C->d0 = thread_arg->cy_r2C;	cy_i2C->d0 = thread_arg->cy_i2C;
			cy_r2C->d1 = thread_arg->cy_r2D;	cy_i2C->d1 = thread_arg->cy_i2D;
			cy_r2C->d2 = thread_arg->cy_r2E;	cy_i2C->d2 = thread_arg->cy_i2E;
			cy_r2C->d3 = thread_arg->cy_r2F;	cy_i2C->d3 = thread_arg->cy_i2F;
			cy_r30->d0 = thread_arg->cy_r30;	cy_i30->d0 = thread_arg->cy_i30;
			cy_r30->d1 = thread_arg->cy_r31;	cy_i30->d1 = thread_arg->cy_i31;
			cy_r30->d2 = thread_arg->cy_r32;	cy_i30->d2 = thread_arg->cy_i32;
			cy_r30->d3 = thread_arg->cy_r33;	cy_i30->d3 = thread_arg->cy_i33;
			cy_r34->d0 = thread_arg->cy_r34;	cy_i34->d0 = thread_arg->cy_i34;
			cy_r34->d1 = thread_arg->cy_r35;	cy_i34->d1 = thread_arg->cy_i35;
			cy_r34->d2 = thread_arg->cy_r36;	cy_i34->d2 = thread_arg->cy_i36;
			cy_r34->d3 = thread_arg->cy_r37;	cy_i34->d3 = thread_arg->cy_i37;
			cy_r38->d0 = thread_arg->cy_r38;	cy_i38->d0 = thread_arg->cy_i38;
			cy_r38->d1 = thread_arg->cy_r39;	cy_i38->d1 = thread_arg->cy_i39;
			cy_r38->d2 = thread_arg->cy_r3A;	cy_i38->d2 = thread_arg->cy_i3A;
			cy_r38->d3 = thread_arg->cy_r3B;	cy_i38->d3 = thread_arg->cy_i3B;
			cy_r3C->d0 = thread_arg->cy_r3C;	cy_i3C->d0 = thread_arg->cy_i3C;
			cy_r3C->d1 = thread_arg->cy_r3D;	cy_i3C->d1 = thread_arg->cy_i3D;
			cy_r3C->d2 = thread_arg->cy_r3E;	cy_i3C->d2 = thread_arg->cy_i3E;
			cy_r3C->d3 = thread_arg->cy_r3F;	cy_i3C->d3 = thread_arg->cy_i3F;

		#elif defined(USE_SSE2)

			cy_r00->d0 = thread_arg->cy_r00;	cy_r00->d1 = thread_arg->cy_i00;
			cy_r02->d0 = thread_arg->cy_r01;	cy_r02->d1 = thread_arg->cy_i01;
			cy_r04->d0 = thread_arg->cy_r02;	cy_r04->d1 = thread_arg->cy_i02;
			cy_r06->d0 = thread_arg->cy_r03;	cy_r06->d1 = thread_arg->cy_i03;
			cy_r08->d0 = thread_arg->cy_r04;	cy_r08->d1 = thread_arg->cy_i04;
			cy_r0A->d0 = thread_arg->cy_r05;	cy_r0A->d1 = thread_arg->cy_i05;
			cy_r0C->d0 = thread_arg->cy_r06;	cy_r0C->d1 = thread_arg->cy_i06;
			cy_r0E->d0 = thread_arg->cy_r07;	cy_r0E->d1 = thread_arg->cy_i07;
			cy_r10->d0 = thread_arg->cy_r08;	cy_r10->d1 = thread_arg->cy_i08;
			cy_r12->d0 = thread_arg->cy_r09;	cy_r12->d1 = thread_arg->cy_i09;
			cy_r14->d0 = thread_arg->cy_r0A;	cy_r14->d1 = thread_arg->cy_i0A;
			cy_r16->d0 = thread_arg->cy_r0B;	cy_r16->d1 = thread_arg->cy_i0B;
			cy_r18->d0 = thread_arg->cy_r0C;	cy_r18->d1 = thread_arg->cy_i0C;
			cy_r1A->d0 = thread_arg->cy_r0D;	cy_r1A->d1 = thread_arg->cy_i0D;
			cy_r1C->d0 = thread_arg->cy_r0E;	cy_r1C->d1 = thread_arg->cy_i0E;
			cy_r1E->d0 = thread_arg->cy_r0F;	cy_r1E->d1 = thread_arg->cy_i0F;
			cy_r20->d0 = thread_arg->cy_r10;	cy_r20->d1 = thread_arg->cy_i10;
			cy_r22->d0 = thread_arg->cy_r11;	cy_r22->d1 = thread_arg->cy_i11;
			cy_r24->d0 = thread_arg->cy_r12;	cy_r24->d1 = thread_arg->cy_i12;
			cy_r26->d0 = thread_arg->cy_r13;	cy_r26->d1 = thread_arg->cy_i13;
			cy_r28->d0 = thread_arg->cy_r14;	cy_r28->d1 = thread_arg->cy_i14;
			cy_r2A->d0 = thread_arg->cy_r15;	cy_r2A->d1 = thread_arg->cy_i15;
			cy_r2C->d0 = thread_arg->cy_r16;	cy_r2C->d1 = thread_arg->cy_i16;
			cy_r2E->d0 = thread_arg->cy_r17;	cy_r2E->d1 = thread_arg->cy_i17;
			cy_r30->d0 = thread_arg->cy_r18;	cy_r30->d1 = thread_arg->cy_i18;
			cy_r32->d0 = thread_arg->cy_r19;	cy_r32->d1 = thread_arg->cy_i19;
			cy_r34->d0 = thread_arg->cy_r1A;	cy_r34->d1 = thread_arg->cy_i1A;
			cy_r36->d0 = thread_arg->cy_r1B;	cy_r36->d1 = thread_arg->cy_i1B;
			cy_r38->d0 = thread_arg->cy_r1C;	cy_r38->d1 = thread_arg->cy_i1C;
			cy_r3A->d0 = thread_arg->cy_r1D;	cy_r3A->d1 = thread_arg->cy_i1D;
			cy_r3C->d0 = thread_arg->cy_r1E;	cy_r3C->d1 = thread_arg->cy_i1E;
			cy_r3E->d0 = thread_arg->cy_r1F;	cy_r3E->d1 = thread_arg->cy_i1F;

			cy_i00->d0 = thread_arg->cy_r20;	cy_i00->d1 = thread_arg->cy_i20;
			cy_i02->d0 = thread_arg->cy_r21;	cy_i02->d1 = thread_arg->cy_i21;
			cy_i04->d0 = thread_arg->cy_r22;	cy_i04->d1 = thread_arg->cy_i22;
			cy_i06->d0 = thread_arg->cy_r23;	cy_i06->d1 = thread_arg->cy_i23;
			cy_i08->d0 = thread_arg->cy_r24;	cy_i08->d1 = thread_arg->cy_i24;
			cy_i0A->d0 = thread_arg->cy_r25;	cy_i0A->d1 = thread_arg->cy_i25;
			cy_i0C->d0 = thread_arg->cy_r26;	cy_i0C->d1 = thread_arg->cy_i26;
			cy_i0E->d0 = thread_arg->cy_r27;	cy_i0E->d1 = thread_arg->cy_i27;
			cy_i10->d0 = thread_arg->cy_r28;	cy_i10->d1 = thread_arg->cy_i28;
			cy_i12->d0 = thread_arg->cy_r29;	cy_i12->d1 = thread_arg->cy_i29;
			cy_i14->d0 = thread_arg->cy_r2A;	cy_i14->d1 = thread_arg->cy_i2A;
			cy_i16->d0 = thread_arg->cy_r2B;	cy_i16->d1 = thread_arg->cy_i2B;
			cy_i18->d0 = thread_arg->cy_r2C;	cy_i18->d1 = thread_arg->cy_i2C;
			cy_i1A->d0 = thread_arg->cy_r2D;	cy_i1A->d1 = thread_arg->cy_i2D;
			cy_i1C->d0 = thread_arg->cy_r2E;	cy_i1C->d1 = thread_arg->cy_i2E;
			cy_i1E->d0 = thread_arg->cy_r2F;	cy_i1E->d1 = thread_arg->cy_i2F;
			cy_i20->d0 = thread_arg->cy_r30;	cy_i20->d1 = thread_arg->cy_i30;
			cy_i22->d0 = thread_arg->cy_r31;	cy_i22->d1 = thread_arg->cy_i31;
			cy_i24->d0 = thread_arg->cy_r32;	cy_i24->d1 = thread_arg->cy_i32;
			cy_i26->d0 = thread_arg->cy_r33;	cy_i26->d1 = thread_arg->cy_i33;
			cy_i28->d0 = thread_arg->cy_r34;	cy_i28->d1 = thread_arg->cy_i34;
			cy_i2A->d0 = thread_arg->cy_r35;	cy_i2A->d1 = thread_arg->cy_i35;
			cy_i2C->d0 = thread_arg->cy_r36;	cy_i2C->d1 = thread_arg->cy_i36;
			cy_i2E->d0 = thread_arg->cy_r37;	cy_i2E->d1 = thread_arg->cy_i37;
			cy_i30->d0 = thread_arg->cy_r38;	cy_i30->d1 = thread_arg->cy_i38;
			cy_i32->d0 = thread_arg->cy_r39;	cy_i32->d1 = thread_arg->cy_i39;
			cy_i34->d0 = thread_arg->cy_r3A;	cy_i34->d1 = thread_arg->cy_i3A;
			cy_i36->d0 = thread_arg->cy_r3B;	cy_i36->d1 = thread_arg->cy_i3B;
			cy_i38->d0 = thread_arg->cy_r3C;	cy_i38->d1 = thread_arg->cy_i3C;
			cy_i3A->d0 = thread_arg->cy_r3D;	cy_i3A->d1 = thread_arg->cy_i3D;
			cy_i3C->d0 = thread_arg->cy_r3E;	cy_i3C->d1 = thread_arg->cy_i3E;
			cy_i3E->d0 = thread_arg->cy_r3F;	cy_i3E->d1 = thread_arg->cy_i3F;

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
			cy_r0A = thread_arg->cy_r0A;		cy_i0A = thread_arg->cy_i0A;
			cy_r0B = thread_arg->cy_r0B;		cy_i0B = thread_arg->cy_i0B;
			cy_r0C = thread_arg->cy_r0C;		cy_i0C = thread_arg->cy_i0C;
			cy_r0D = thread_arg->cy_r0D;		cy_i0D = thread_arg->cy_i0D;
			cy_r0E = thread_arg->cy_r0E;		cy_i0E = thread_arg->cy_i0E;
			cy_r0F = thread_arg->cy_r0F;		cy_i0F = thread_arg->cy_i0F;
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
			cy_r1A = thread_arg->cy_r1A;		cy_i1A = thread_arg->cy_i1A;
			cy_r1B = thread_arg->cy_r1B;		cy_i1B = thread_arg->cy_i1B;
			cy_r1C = thread_arg->cy_r1C;		cy_i1C = thread_arg->cy_i1C;
			cy_r1D = thread_arg->cy_r1D;		cy_i1D = thread_arg->cy_i1D;
			cy_r1E = thread_arg->cy_r1E;		cy_i1E = thread_arg->cy_i1E;
			cy_r1F = thread_arg->cy_r1F;		cy_i1F = thread_arg->cy_i1F;
			cy_r20 = thread_arg->cy_r20;		cy_i20 = thread_arg->cy_i20;
			cy_r21 = thread_arg->cy_r21;		cy_i21 = thread_arg->cy_i21;
			cy_r22 = thread_arg->cy_r22;		cy_i22 = thread_arg->cy_i22;
			cy_r23 = thread_arg->cy_r23;		cy_i23 = thread_arg->cy_i23;
			cy_r24 = thread_arg->cy_r24;		cy_i24 = thread_arg->cy_i24;
			cy_r25 = thread_arg->cy_r25;		cy_i25 = thread_arg->cy_i25;
			cy_r26 = thread_arg->cy_r26;		cy_i26 = thread_arg->cy_i26;
			cy_r27 = thread_arg->cy_r27;		cy_i27 = thread_arg->cy_i27;
			cy_r28 = thread_arg->cy_r28;		cy_i28 = thread_arg->cy_i28;
			cy_r29 = thread_arg->cy_r29;		cy_i29 = thread_arg->cy_i29;
			cy_r2A = thread_arg->cy_r2A;		cy_i2A = thread_arg->cy_i2A;
			cy_r2B = thread_arg->cy_r2B;		cy_i2B = thread_arg->cy_i2B;
			cy_r2C = thread_arg->cy_r2C;		cy_i2C = thread_arg->cy_i2C;
			cy_r2D = thread_arg->cy_r2D;		cy_i2D = thread_arg->cy_i2D;
			cy_r2E = thread_arg->cy_r2E;		cy_i2E = thread_arg->cy_i2E;
			cy_r2F = thread_arg->cy_r2F;		cy_i2F = thread_arg->cy_i2F;
			cy_r30 = thread_arg->cy_r30;		cy_i30 = thread_arg->cy_i30;
			cy_r31 = thread_arg->cy_r31;		cy_i31 = thread_arg->cy_i31;
			cy_r32 = thread_arg->cy_r32;		cy_i32 = thread_arg->cy_i32;
			cy_r33 = thread_arg->cy_r33;		cy_i33 = thread_arg->cy_i33;
			cy_r34 = thread_arg->cy_r34;		cy_i34 = thread_arg->cy_i34;
			cy_r35 = thread_arg->cy_r35;		cy_i35 = thread_arg->cy_i35;
			cy_r36 = thread_arg->cy_r36;		cy_i36 = thread_arg->cy_i36;
			cy_r37 = thread_arg->cy_r37;		cy_i37 = thread_arg->cy_i37;
			cy_r38 = thread_arg->cy_r38;		cy_i38 = thread_arg->cy_i38;
			cy_r39 = thread_arg->cy_r39;		cy_i39 = thread_arg->cy_i39;
			cy_r3A = thread_arg->cy_r3A;		cy_i3A = thread_arg->cy_i3A;
			cy_r3B = thread_arg->cy_r3B;		cy_i3B = thread_arg->cy_i3B;
			cy_r3C = thread_arg->cy_r3C;		cy_i3C = thread_arg->cy_i3C;
			cy_r3D = thread_arg->cy_r3D;		cy_i3D = thread_arg->cy_i3D;
			cy_r3E = thread_arg->cy_r3E;		cy_i3E = thread_arg->cy_i3E;
			cy_r3F = thread_arg->cy_r3F;		cy_i3F = thread_arg->cy_i3F;

		#endif
		}

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
			for(j = jstart; j < jhi; j += stride)	// Stride = 4 reals for SSE2, 8 for AVX
			{
				j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1+RE_IM_STRIDE;

			/*...The radix-64 DIT pass is here:	*/

			#ifdef USE_SSE2

			/* Gather the needed data (64 64-bit complex, i.e. 128 64-bit reals) and do 8 twiddleless length-8 subtransforms: */
	// Because of roots-sign diffs between SSE2 and scalar macros, 1/7 2/6 3/5 swapped in DIT_0TWIDDLE outputs!
				add0 = &a[j1]      ; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r00, isrt2)

				add0 = &a[j1] + p08; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r10, isrt2)

				add0 = &a[j1] + p10; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r20, isrt2)

				add0 = &a[j1] + p18; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r30, isrt2)

				add0 = &a[j1] + p20; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r40, isrt2)

				add0 = &a[j1] + p28; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r50, isrt2)

				add0 = &a[j1] + p30; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r60, isrt2)

				add0 = &a[j1] + p38; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r70, isrt2)

	// Because of roots-sign diffs between SSE2 and scalar macros, 1/7 2/6 3/5 swapped in DIT_0TWIDDLE outputs!

			/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dit_pass1 for details: */
			// Note: 1st of the 15 sincos args in each call to SSE2_RADIX8_DIT_TWIDDLE_OOP is the basic isrt2 needed for
			// radix-8. This is a workaround of GCC's 30-arg limit for inline ASM macros, which proves a royal pain here.
			//...and another kludge for the 30-arg limit: put a copy of (vec_dbl)2.0 into the first of each s1p*r outputs:
				VEC_DBL_INIT(s1p01r,2.0);
				VEC_DBL_INIT(s1p02r,2.0);
				VEC_DBL_INIT(s1p03r,2.0);
				VEC_DBL_INIT(s1p04r,2.0);
				VEC_DBL_INIT(s1p05r,2.0);
				VEC_DBL_INIT(s1p06r,2.0);
				VEC_DBL_INIT(s1p07r,2.0);

				// Block 0: jt = j1;	jp = j2;
				SSE2_RADIX8_DIT_0TWIDDLE_OOP(	// This outputs o[07654321], so reverse o-index order of latter 7 outputs
					r00,r10,r20,r30,r40,r50,r60,r70,
					s1p00r,s1p38r,s1p30r,s1p28r,s1p20r,s1p18r,s1p10r,s1p08r, isrt2
				);
			/* 0-index block has all-unity twiddles: **
				VEC_DBL_INIT(s1p00r,2.0);
				SSE2_RADIX8_DIT_TWIDDLE_OOP(
					r00,r10,r20,r30,r40,r50,r60,r70,
					s1p00r,s1p20r,s1p10r,s1p30r,s1p08r,s1p28r,s1p18r,s1p38r,
					cc0,ss0,cc0,ss0,cc0,ss0,cc0,ss0,cc0,ss0,cc0,ss0,cc0,ss0
				);
			*/
				// Block 4: jt = j1 + p04;	jp = j2 + p04;
				SSE2_RADIX8_DIT_TWIDDLE_OOP(
					r08,r18,r28,r38,r48,r58,r68,r78,
					s1p04r,s1p24r,s1p14r,s1p34r,s1p0cr,s1p2cr,s1p1cr,s1p3cr,
					ss0,cc0, isrt2,isrt2, nisrt2,isrt2, cc4,ss4, nss4,cc4, ss4,cc4, ncc4,ss4
				);
				// Block 2: jt = j1 + p02;	jp = j2 + p02;
				SSE2_RADIX8_DIT_TWIDDLE_OOP(
					r0C,r1C,r2C,r3C,r4C,r5C,r6C,r7C,	// 2/6 swap ==> r*4 / r*C swap
					s1p02r,s1p22r,s1p12r,s1p32r,s1p0ar,s1p2ar,s1p1ar,s1p3ar,
					isrt2,isrt2,cc4,ss4,ss4,cc4,cc2,ss2,ss6,cc6,cc6,ss6,ss2,cc2
				);
				// Block 6: jt = j1 + p06;	jp = j2 + p06;
				SSE2_RADIX8_DIT_TWIDDLE_OOP(
					r04,r14,r24,r34,r44,r54,r64,r74,	// 2/6 swap ==> r*4 / r*C swap
					s1p06r,s1p26r,s1p16r,s1p36r,s1p0er,s1p2er,s1p1er,s1p3er,
					nisrt2,isrt2,ss4,cc4,ncc4,nss4,cc6,ss6,ncc2,ss2,nss2,cc2,nss6,ncc6
				);
				// Block 1: jt = j1 + p01;	jp = j2 + p01;
				SSE2_RADIX8_DIT_TWIDDLE_OOP(
					r0E,r1E,r2E,r3E,r4E,r5E,r6E,r7E,	// 1/7 swap ==> r*2 / r*E swap
					s1p01r,s1p21r,s1p11r,s1p31r,s1p09r,s1p29r,s1p19r,s1p39r,
					cc4,ss4,cc2,ss2,cc6,ss6,cc1,ss1,cc5,ss5,cc3,ss3,cc7,ss7
				);
				// Block 5: jt = j1 + p05;	jp = j2 + p05;
				SSE2_RADIX8_DIT_TWIDDLE_OOP(
					r06,r16,r26,r36,r46,r56,r66,r76,
					s1p05r,s1p25r,s1p15r,s1p35r,s1p0dr,s1p2dr,s1p1dr,s1p3dr,
					nss4,cc4,ss6,cc6,ncc2,ss2,cc5,ss5,ncc7,ss7,ss1,cc1,ncc3,nss3
				);
				// Block 3: jt = j1 + p03;	jp = j2 + p03;
				SSE2_RADIX8_DIT_TWIDDLE_OOP(
					r0A,r1A,r2A,r3A,r4A,r5A,r6A,r7A,	// 3/5 swap ==> r*6 / r*A swap
					s1p03r,s1p23r,s1p13r,s1p33r,s1p0br,s1p2br,s1p1br,s1p3br,
					ss4,cc4,cc6,ss6,nss2,cc2,cc3,ss3,ss1,cc1,ss7,cc7,nss5,cc5
				);
				// Block 7: jt = j1 + p07;	jp = j2 + p07;
				SSE2_RADIX8_DIT_TWIDDLE_OOP(
					r02,r12,r22,r32,r42,r52,r62,r72,	// 1/7 swap ==> r*2 / r*E swap
					s1p07r,s1p27r,s1p17r,s1p37r,s1p0fr,s1p2fr,s1p1fr,s1p3fr,
					ncc4,ss4,ss2,cc2,nss6,ncc6,cc7,ss7,ncc3,nss3,nss5,cc5,ss1,ncc1
				);

			#else	// USE_SSE2 = False, or non-64-bit-GCC:

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

			/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dit_pass1 for details: */

				// Block 0: jt = j1;	jp = j2;
				/* 0-index block has all-unity twiddles: Remember, the twiddleless DIT also bit-reverses its outputs, so a_p* terms appear in-order here: */
				RADIX_08_DIT_OOP(
					t00,t01,t10,t11,t20,t21,t30,t31,t40,t41,t50,t51,t60,t61,t70,t71,
					ar_p00,ai_p00,ar_p08,ai_p08,ar_p10,ai_p10,ar_p18,ai_p18,ar_p20,ai_p20,ar_p28,ai_p28,ar_p30,ai_p30,ar_p38,ai_p38
				);
				// Block 4: jt = j1 + p04;	jp = j2 + p04;
				RADIX_08_DIT_TWIDDLE_OOP(
					t08,t09,t18,t19,t28,t29,t38,t39,t48,t49,t58,t59,t68,t69,t78,t79,
					ar_p04,ai_p04,ar_p24,ai_p24,ar_p14,ai_p14,ar_p34,ai_p34,ar_p0c,ai_p0c,ar_p2c,ai_p2c,ar_p1c,ai_p1c,ar_p3c,ai_p3c,
					0,1,ISRT2,ISRT2,-ISRT2,ISRT2,c64_4,s64_4,-s64_4,c64_4,s64_4,c64_4,-c64_4,s64_4
				);
				// Block 2: jt = j1 + p02;	jp = j2 + p02;
				RADIX_08_DIT_TWIDDLE_OOP(
					t04,t05,t14,t15,t24,t25,t34,t35,t44,t45,t54,t55,t64,t65,t74,t75,
					ar_p02,ai_p02,ar_p22,ai_p22,ar_p12,ai_p12,ar_p32,ai_p32,ar_p0a,ai_p0a,ar_p2a,ai_p2a,ar_p1a,ai_p1a,ar_p3a,ai_p3a,
					ISRT2,ISRT2,c64_4,s64_4,s64_4,c64_4,c64_2,s64_2,s64_6,c64_6,c64_6,s64_6,s64_2,c64_2
				);
				// Block 6: jt = j1 + p06;	jp = j2 + p06;
				RADIX_08_DIT_TWIDDLE_OOP(
					t0C,t0D,t1C,t1D,t2C,t2D,t3C,t3D,t4C,t4D,t5C,t5D,t6C,t6D,t7C,t7D,
					ar_p06,ai_p06,ar_p26,ai_p26,ar_p16,ai_p16,ar_p36,ai_p36,ar_p0e,ai_p0e,ar_p2e,ai_p2e,ar_p1e,ai_p1e,ar_p3e,ai_p3e,
					-ISRT2,ISRT2,s64_4,c64_4,-c64_4,-s64_4,c64_6,s64_6,-c64_2,s64_2,-s64_2,c64_2,-s64_6,-c64_6
				);
				// Block 1: jt = j1 + p01;	jp = j2 + p01;
				RADIX_08_DIT_TWIDDLE_OOP(
					t02,t03,t12,t13,t22,t23,t32,t33,t42,t43,t52,t53,t62,t63,t72,t73,
					ar_p01,ai_p01,ar_p21,ai_p21,ar_p11,ai_p11,ar_p31,ai_p31,ar_p09,ai_p09,ar_p29,ai_p29,ar_p19,ai_p19,ar_p39,ai_p39,
					c64_4,s64_4,c64_2,s64_2,c64_6,s64_6,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7
				);
				// Block 5: jt = j1 + p05;	jp = j2 + p05;
				RADIX_08_DIT_TWIDDLE_OOP(
					t0A,t0B,t1A,t1B,t2A,t2B,t3A,t3B,t4A,t4B,t5A,t5B,t6A,t6B,t7A,t7B,
					ar_p05,ai_p05,ar_p25,ai_p25,ar_p15,ai_p15,ar_p35,ai_p35,ar_p0d,ai_p0d,ar_p2d,ai_p2d,ar_p1d,ai_p1d,ar_p3d,ai_p3d,
					-s64_4,c64_4,s64_6,c64_6,-c64_2,s64_2,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3
				);
				// Block 3: jt = j1 + p03;	jp = j2 + p03;
				RADIX_08_DIT_TWIDDLE_OOP(
					t06,t07,t16,t17,t26,t27,t36,t37,t46,t47,t56,t57,t66,t67,t76,t77,
					ar_p03,ai_p03,ar_p23,ai_p23,ar_p13,ai_p13,ar_p33,ai_p33,ar_p0b,ai_p0b,ar_p2b,ai_p2b,ar_p1b,ai_p1b,ar_p3b,ai_p3b,
					s64_4,c64_4,c64_6,s64_6,-s64_2,c64_2,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5
				);
				// Block 7: jt = j1 + p07;	jp = j2 + p07;
				RADIX_08_DIT_TWIDDLE_OOP(
					t0E,t0F,t1E,t1F,t2E,t2F,t3E,t3F,t4E,t4F,t5E,t5F,t6E,t6F,t7E,t7F,
					ar_p07,ai_p07,ar_p27,ai_p27,ar_p17,ai_p17,ar_p37,ai_p37,ar_p0f,ai_p0f,ar_p2f,ai_p2f,ar_p1f,ai_p1f,ar_p3f,ai_p3f,
					-c64_4,s64_4,s64_2,c64_2,-s64_6,-c64_6,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1
				);

			#endif	// USE_SSE2 ?

			/*...Now do the carries. Since the outputs would
			normally be getting dispatched to 32 separate blocks of the A-array, we need 32 separate carries.	*/

			/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/
			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
			#ifdef USE_AVX

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

				l= j & (nwt-1);						tmp = half_arr + 64;	/* ptr to local storage for the doubled wtl,wtn terms: */
				n_minus_sil  ->d0 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d0 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d0 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d0 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+2) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d1 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d1 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d1 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d1 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+4) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d2 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d2 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d2 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d2 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+6) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d3 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d3 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d3 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d3 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				AVX_cmplx_carry_norm_pow2_errcheck0_X4(s1p00r,add1,add2,add3,cy_r00,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p04r,add1,add2,add3,cy_r04,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p08r,add1,add2,add3,cy_r08,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p0cr,add1,add2,add3,cy_r0C,bjmodn0C,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p10r,add1,add2,add3,cy_r10,bjmodn10,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p14r,add1,add2,add3,cy_r14,bjmodn14,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p18r,add1,add2,add3,cy_r18,bjmodn18,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p1cr,add1,add2,add3,cy_r1C,bjmodn1C,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p20r,add1,add2,add3,cy_r20,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p24r,add1,add2,add3,cy_r24,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p28r,add1,add2,add3,cy_r28,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p2cr,add1,add2,add3,cy_r2C,bjmodn2C,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p30r,add1,add2,add3,cy_r30,bjmodn30,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p34r,add1,add2,add3,cy_r34,bjmodn34,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p38r,add1,add2,add3,cy_r38,bjmodn38,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				AVX_cmplx_carry_norm_pow2_errcheck1_X4(s1p3cr,add1,add2,add3,cy_r3C,bjmodn3C,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);

				co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

			#elif defined(USE_SSE2)

				l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 32 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				ctmp->re = ctmp->im = wtl;		++ctmp;
				ctmp->re = ctmp->im = wtn;		++ctmp;
				ctmp->re = ctmp->im = wtlp1;	++ctmp;
				ctmp->re = ctmp->im = wtnm1;

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_pow2_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p08r,add1,add2,add3,cy_r08,cy_r0A,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p0cr,add1,add2,add3,cy_r0C,cy_r0E,bjmodn0C,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p10r,add1,add2,add3,cy_r10,cy_r12,bjmodn10,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p14r,add1,add2,add3,cy_r14,cy_r16,bjmodn14,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p18r,add1,add2,add3,cy_r18,cy_r1A,bjmodn18,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p1cr,add1,add2,add3,cy_r1C,cy_r1E,bjmodn1C,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p28r,add1,add2,add3,cy_r28,cy_r2A,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p2cr,add1,add2,add3,cy_r2C,cy_r2E,bjmodn2C,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p30r,add1,add2,add3,cy_r30,cy_r32,bjmodn30,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p34r,add1,add2,add3,cy_r34,cy_r36,bjmodn34,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p38r,add1,add2,add3,cy_r38,cy_r3A,bjmodn38,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck1_2B(s1p3cr,add1,add2,add3,cy_r3C,cy_r3E,bjmodn3C,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_pow2_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r0A,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p0cr,add1,add2,add3,cy_r0C,cy_r0E,bjmodn0C,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p10r,add1,add2,add3,cy_r10,cy_r12,bjmodn10,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p14r,add1,add2,add3,cy_r14,cy_r16,bjmodn14,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p18r,add1,add2,add3,cy_r18,cy_r1A,bjmodn18,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p1cr,add1,add2,add3,cy_r1C,cy_r1E,bjmodn1C,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p28r,add1,add2,add3,cy_r28,cy_r2A,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p2cr,add1,add2,add3,cy_r2C,cy_r2E,bjmodn2C,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p30r,add1,add2,add3,cy_r30,cy_r32,bjmodn30,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p34r,add1,add2,add3,cy_r34,cy_r36,bjmodn34,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p38r,add1,add2,add3,cy_r38,cy_r3A,bjmodn38,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck1_2B (s1p3cr,add1,add2,add3,cy_r3C,cy_r3E,bjmodn3C,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			  #endif

				l= (j+2) & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				ctmp->re = ctmp->im = wtl;		++ctmp;
				ctmp->re = ctmp->im = wtn;		++ctmp;
				ctmp->re = ctmp->im = wtlp1;	++ctmp;
				ctmp->re = ctmp->im = wtnm1;

				co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r0A,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p0cr,add1,add2,     cy_r0C,cy_r0E,bjmodn0C,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p10r,add1,add2,     cy_r10,cy_r12,bjmodn10,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p14r,add1,add2,     cy_r14,cy_r16,bjmodn14,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p18r,add1,add2,     cy_r18,cy_r1A,bjmodn18,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p1cr,add1,add2,     cy_r1C,cy_r1E,bjmodn1C,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p28r,add1,add2,     cy_r28,cy_r2A,bjmodn28,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p2cr,add1,add2,     cy_r2C,cy_r2E,bjmodn2C,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p30r,add1,add2,     cy_r30,cy_r32,bjmodn30,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p34r,add1,add2,     cy_r34,cy_r36,bjmodn34,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p38r,add1,add2,     cy_r38,cy_r3A,bjmodn38,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_errcheck2_2B(s1p3cr,add1,add2,     cy_r3C,cy_r3E,bjmodn3C,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p08r,add1,add2,     cy_r08,cy_r0A,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p0cr,add1,add2,     cy_r0C,cy_r0E,bjmodn0C,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p10r,add1,add2,     cy_r10,cy_r12,bjmodn10,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p14r,add1,add2,     cy_r14,cy_r16,bjmodn14,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p18r,add1,add2,     cy_r18,cy_r1A,bjmodn18,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p1cr,add1,add2,     cy_r1C,cy_r1E,bjmodn1C,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p28r,add1,add2,     cy_r28,cy_r2A,bjmodn28,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p2cr,add1,add2,     cy_r2C,cy_r2E,bjmodn2C,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p30r,add1,add2,     cy_r30,cy_r32,bjmodn30,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p34r,add1,add2,     cy_r34,cy_r36,bjmodn34,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p38r,add1,add2,     cy_r38,cy_r3A,bjmodn38,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
				SSE2_cmplx_carry_norm_pow2_nocheck2_2B (s1p3cr,add1,add2,     cy_r3C,cy_r3E,bjmodn3C,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
			  #endif

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

			#else	// Scalar-double mode:

				l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 32 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				/*...set0 is slightly different from others:	*/
			   cmplx_carry_norm_pow2_errcheck0(ar_p00,ai_p00,cy_r00,bjmodn00);
				cmplx_carry_norm_pow2_errcheck(ar_p01,ai_p01,cy_r01,bjmodn01,0x01);
				cmplx_carry_norm_pow2_errcheck(ar_p02,ai_p02,cy_r02,bjmodn02,0x02);
				cmplx_carry_norm_pow2_errcheck(ar_p03,ai_p03,cy_r03,bjmodn03,0x03);
				cmplx_carry_norm_pow2_errcheck(ar_p04,ai_p04,cy_r04,bjmodn04,0x04);
				cmplx_carry_norm_pow2_errcheck(ar_p05,ai_p05,cy_r05,bjmodn05,0x05);
				cmplx_carry_norm_pow2_errcheck(ar_p06,ai_p06,cy_r06,bjmodn06,0x06);
				cmplx_carry_norm_pow2_errcheck(ar_p07,ai_p07,cy_r07,bjmodn07,0x07);
				cmplx_carry_norm_pow2_errcheck(ar_p08,ai_p08,cy_r08,bjmodn08,0x08);
				cmplx_carry_norm_pow2_errcheck(ar_p09,ai_p09,cy_r09,bjmodn09,0x09);
				cmplx_carry_norm_pow2_errcheck(ar_p0a,ai_p0a,cy_r0A,bjmodn0A,0x0A);
				cmplx_carry_norm_pow2_errcheck(ar_p0b,ai_p0b,cy_r0B,bjmodn0B,0x0B);
				cmplx_carry_norm_pow2_errcheck(ar_p0c,ai_p0c,cy_r0C,bjmodn0C,0x0C);
				cmplx_carry_norm_pow2_errcheck(ar_p0d,ai_p0d,cy_r0D,bjmodn0D,0x0D);
				cmplx_carry_norm_pow2_errcheck(ar_p0e,ai_p0e,cy_r0E,bjmodn0E,0x0E);
				cmplx_carry_norm_pow2_errcheck(ar_p0f,ai_p0f,cy_r0F,bjmodn0F,0x0F);
				cmplx_carry_norm_pow2_errcheck(ar_p10,ai_p10,cy_r10,bjmodn10,0x10);
				cmplx_carry_norm_pow2_errcheck(ar_p11,ai_p11,cy_r11,bjmodn11,0x11);
				cmplx_carry_norm_pow2_errcheck(ar_p12,ai_p12,cy_r12,bjmodn12,0x12);
				cmplx_carry_norm_pow2_errcheck(ar_p13,ai_p13,cy_r13,bjmodn13,0x13);
				cmplx_carry_norm_pow2_errcheck(ar_p14,ai_p14,cy_r14,bjmodn14,0x14);
				cmplx_carry_norm_pow2_errcheck(ar_p15,ai_p15,cy_r15,bjmodn15,0x15);
				cmplx_carry_norm_pow2_errcheck(ar_p16,ai_p16,cy_r16,bjmodn16,0x16);
				cmplx_carry_norm_pow2_errcheck(ar_p17,ai_p17,cy_r17,bjmodn17,0x17);
				cmplx_carry_norm_pow2_errcheck(ar_p18,ai_p18,cy_r18,bjmodn18,0x18);
				cmplx_carry_norm_pow2_errcheck(ar_p19,ai_p19,cy_r19,bjmodn19,0x19);
				cmplx_carry_norm_pow2_errcheck(ar_p1a,ai_p1a,cy_r1A,bjmodn1A,0x1A);
				cmplx_carry_norm_pow2_errcheck(ar_p1b,ai_p1b,cy_r1B,bjmodn1B,0x1B);
				cmplx_carry_norm_pow2_errcheck(ar_p1c,ai_p1c,cy_r1C,bjmodn1C,0x1C);
				cmplx_carry_norm_pow2_errcheck(ar_p1d,ai_p1d,cy_r1D,bjmodn1D,0x1D);
				cmplx_carry_norm_pow2_errcheck(ar_p1e,ai_p1e,cy_r1E,bjmodn1E,0x1E);
				cmplx_carry_norm_pow2_errcheck(ar_p1f,ai_p1f,cy_r1F,bjmodn1F,0x1F);
				cmplx_carry_norm_pow2_errcheck(ar_p20,ai_p20,cy_r20,bjmodn20,0x20);
				cmplx_carry_norm_pow2_errcheck(ar_p21,ai_p21,cy_r21,bjmodn21,0x21);
				cmplx_carry_norm_pow2_errcheck(ar_p22,ai_p22,cy_r22,bjmodn22,0x22);
				cmplx_carry_norm_pow2_errcheck(ar_p23,ai_p23,cy_r23,bjmodn23,0x23);
				cmplx_carry_norm_pow2_errcheck(ar_p24,ai_p24,cy_r24,bjmodn24,0x24);
				cmplx_carry_norm_pow2_errcheck(ar_p25,ai_p25,cy_r25,bjmodn25,0x25);
				cmplx_carry_norm_pow2_errcheck(ar_p26,ai_p26,cy_r26,bjmodn26,0x26);
				cmplx_carry_norm_pow2_errcheck(ar_p27,ai_p27,cy_r27,bjmodn27,0x27);
				cmplx_carry_norm_pow2_errcheck(ar_p28,ai_p28,cy_r28,bjmodn28,0x28);
				cmplx_carry_norm_pow2_errcheck(ar_p29,ai_p29,cy_r29,bjmodn29,0x29);
				cmplx_carry_norm_pow2_errcheck(ar_p2a,ai_p2a,cy_r2A,bjmodn2A,0x2A);
				cmplx_carry_norm_pow2_errcheck(ar_p2b,ai_p2b,cy_r2B,bjmodn2B,0x2B);
				cmplx_carry_norm_pow2_errcheck(ar_p2c,ai_p2c,cy_r2C,bjmodn2C,0x2C);
				cmplx_carry_norm_pow2_errcheck(ar_p2d,ai_p2d,cy_r2D,bjmodn2D,0x2D);
				cmplx_carry_norm_pow2_errcheck(ar_p2e,ai_p2e,cy_r2E,bjmodn2E,0x2E);
				cmplx_carry_norm_pow2_errcheck(ar_p2f,ai_p2f,cy_r2F,bjmodn2F,0x2F);
				cmplx_carry_norm_pow2_errcheck(ar_p30,ai_p30,cy_r30,bjmodn30,0x30);
				cmplx_carry_norm_pow2_errcheck(ar_p31,ai_p31,cy_r31,bjmodn31,0x31);
				cmplx_carry_norm_pow2_errcheck(ar_p32,ai_p32,cy_r32,bjmodn32,0x32);
				cmplx_carry_norm_pow2_errcheck(ar_p33,ai_p33,cy_r33,bjmodn33,0x33);
				cmplx_carry_norm_pow2_errcheck(ar_p34,ai_p34,cy_r34,bjmodn34,0x34);
				cmplx_carry_norm_pow2_errcheck(ar_p35,ai_p35,cy_r35,bjmodn35,0x35);
				cmplx_carry_norm_pow2_errcheck(ar_p36,ai_p36,cy_r36,bjmodn36,0x36);
				cmplx_carry_norm_pow2_errcheck(ar_p37,ai_p37,cy_r37,bjmodn37,0x37);
				cmplx_carry_norm_pow2_errcheck(ar_p38,ai_p38,cy_r38,bjmodn38,0x38);
				cmplx_carry_norm_pow2_errcheck(ar_p39,ai_p39,cy_r39,bjmodn39,0x39);
				cmplx_carry_norm_pow2_errcheck(ar_p3a,ai_p3a,cy_r3A,bjmodn3A,0x3A);
				cmplx_carry_norm_pow2_errcheck(ar_p3b,ai_p3b,cy_r3B,bjmodn3B,0x3B);
				cmplx_carry_norm_pow2_errcheck(ar_p3c,ai_p3c,cy_r3C,bjmodn3C,0x3C);
				cmplx_carry_norm_pow2_errcheck(ar_p3d,ai_p3d,cy_r3D,bjmodn3D,0x3D);
				cmplx_carry_norm_pow2_errcheck(ar_p3e,ai_p3e,cy_r3E,bjmodn3E,0x3E);
				cmplx_carry_norm_pow2_errcheck(ar_p3f,ai_p3f,cy_r3F,bjmodn3F,0x3F);

				i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
				co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
						 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

			#endif	// USE_AVX?
			}
			else	/* Fermat-mod carry in SIMD mode */
			{
			#ifdef USE_AVX

				// For a description of the data movement for Fermat-mod carries in SSE2 mode, see radix16_ditN_cy_dif1.c.
				// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.

				tmp = half_arr+2;
				VEC_DBL_INIT(tmp, scale);
				/* Get the needed Nth root of -1: */
				add1 = (double *)&rn0[0];
				add2 = (double *)&rn1[0];

				idx_offset = j;
				idx_incr = NDIVR;

				tmp = base_negacyclic_root;	tm2 = tmp+1;

				// Hi-accuracy version needs 8 copies of each base root:
				l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				VEC_DBL_INIT(tmp+ 32,wt_re);	VEC_DBL_INIT(tm2+ 32,wt_im);
				VEC_DBL_INIT(tmp+ 40,wt_re);	VEC_DBL_INIT(tm2+ 40,wt_im);
				VEC_DBL_INIT(tmp+ 48,wt_re);	VEC_DBL_INIT(tm2+ 48,wt_im);
				VEC_DBL_INIT(tmp+ 56,wt_re);	VEC_DBL_INIT(tm2+ 56,wt_im);
				VEC_DBL_INIT(tmp+ 64,wt_re);	VEC_DBL_INIT(tm2+ 64,wt_im);
				VEC_DBL_INIT(tmp+ 72,wt_re);	VEC_DBL_INIT(tm2+ 72,wt_im);
				VEC_DBL_INIT(tmp+ 80,wt_re);	VEC_DBL_INIT(tm2+ 80,wt_im);
				VEC_DBL_INIT(tmp+ 88,wt_re);	VEC_DBL_INIT(tm2+ 88,wt_im);
				VEC_DBL_INIT(tmp+ 96,wt_re);	VEC_DBL_INIT(tm2+ 96,wt_im);
				VEC_DBL_INIT(tmp+104,wt_re);	VEC_DBL_INIT(tm2+104,wt_im);
				VEC_DBL_INIT(tmp+112,wt_re);	VEC_DBL_INIT(tm2+112,wt_im);
				VEC_DBL_INIT(tmp+120,wt_re);	VEC_DBL_INIT(tm2+120,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				VEC_DBL_INIT(tmp+ 32,wt_re);	VEC_DBL_INIT(tm2+ 32,wt_im);
				VEC_DBL_INIT(tmp+ 40,wt_re);	VEC_DBL_INIT(tm2+ 40,wt_im);
				VEC_DBL_INIT(tmp+ 48,wt_re);	VEC_DBL_INIT(tm2+ 48,wt_im);
				VEC_DBL_INIT(tmp+ 56,wt_re);	VEC_DBL_INIT(tm2+ 56,wt_im);
				VEC_DBL_INIT(tmp+ 64,wt_re);	VEC_DBL_INIT(tm2+ 64,wt_im);
				VEC_DBL_INIT(tmp+ 72,wt_re);	VEC_DBL_INIT(tm2+ 72,wt_im);
				VEC_DBL_INIT(tmp+ 80,wt_re);	VEC_DBL_INIT(tm2+ 80,wt_im);
				VEC_DBL_INIT(tmp+ 88,wt_re);	VEC_DBL_INIT(tm2+ 88,wt_im);
				VEC_DBL_INIT(tmp+ 96,wt_re);	VEC_DBL_INIT(tm2+ 96,wt_im);
				VEC_DBL_INIT(tmp+104,wt_re);	VEC_DBL_INIT(tm2+104,wt_im);
				VEC_DBL_INIT(tmp+112,wt_re);	VEC_DBL_INIT(tm2+112,wt_im);
				VEC_DBL_INIT(tmp+120,wt_re);	VEC_DBL_INIT(tm2+120,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				VEC_DBL_INIT(tmp+ 32,wt_re);	VEC_DBL_INIT(tm2+ 32,wt_im);
				VEC_DBL_INIT(tmp+ 40,wt_re);	VEC_DBL_INIT(tm2+ 40,wt_im);
				VEC_DBL_INIT(tmp+ 48,wt_re);	VEC_DBL_INIT(tm2+ 48,wt_im);
				VEC_DBL_INIT(tmp+ 56,wt_re);	VEC_DBL_INIT(tm2+ 56,wt_im);
				VEC_DBL_INIT(tmp+ 64,wt_re);	VEC_DBL_INIT(tm2+ 64,wt_im);
				VEC_DBL_INIT(tmp+ 72,wt_re);	VEC_DBL_INIT(tm2+ 72,wt_im);
				VEC_DBL_INIT(tmp+ 80,wt_re);	VEC_DBL_INIT(tm2+ 80,wt_im);
				VEC_DBL_INIT(tmp+ 88,wt_re);	VEC_DBL_INIT(tm2+ 88,wt_im);
				VEC_DBL_INIT(tmp+ 96,wt_re);	VEC_DBL_INIT(tm2+ 96,wt_im);
				VEC_DBL_INIT(tmp+104,wt_re);	VEC_DBL_INIT(tm2+104,wt_im);
				VEC_DBL_INIT(tmp+112,wt_re);	VEC_DBL_INIT(tm2+112,wt_im);
				VEC_DBL_INIT(tmp+120,wt_re);	VEC_DBL_INIT(tm2+120,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				VEC_DBL_INIT(tmp+ 32,wt_re);	VEC_DBL_INIT(tm2+ 32,wt_im);
				VEC_DBL_INIT(tmp+ 40,wt_re);	VEC_DBL_INIT(tm2+ 40,wt_im);
				VEC_DBL_INIT(tmp+ 48,wt_re);	VEC_DBL_INIT(tm2+ 48,wt_im);
				VEC_DBL_INIT(tmp+ 56,wt_re);	VEC_DBL_INIT(tm2+ 56,wt_im);
				VEC_DBL_INIT(tmp+ 64,wt_re);	VEC_DBL_INIT(tm2+ 64,wt_im);
				VEC_DBL_INIT(tmp+ 72,wt_re);	VEC_DBL_INIT(tm2+ 72,wt_im);
				VEC_DBL_INIT(tmp+ 80,wt_re);	VEC_DBL_INIT(tm2+ 80,wt_im);
				VEC_DBL_INIT(tmp+ 88,wt_re);	VEC_DBL_INIT(tm2+ 88,wt_im);
				VEC_DBL_INIT(tmp+ 96,wt_re);	VEC_DBL_INIT(tm2+ 96,wt_im);
				VEC_DBL_INIT(tmp+104,wt_re);	VEC_DBL_INIT(tm2+104,wt_im);
				VEC_DBL_INIT(tmp+112,wt_re);	VEC_DBL_INIT(tm2+112,wt_im);
				VEC_DBL_INIT(tmp+120,wt_re);	VEC_DBL_INIT(tm2+120,wt_im);

				// AVX-custom 4-way carry macro - each contains 4 of the RADIX stride-n/RADIX-separated carries
				// (processed independently in parallel), and steps through sequential-data indicies j,j+2,j+4,j+6:

				// The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl) = 2^12
				// which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
				// precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
				// uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls:
				tmp = base_negacyclic_root+  0;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p00r,tmp,0x1000,cy_r00,cy_i00,half_arr,sign_mask);
				tmp = base_negacyclic_root+  8;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p04r,tmp,0x0f40,cy_r04,cy_i04,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 16;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p08r,tmp,0x0e80,cy_r08,cy_i08,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 24;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p0cr,tmp,0x0dc0,cy_r0C,cy_i0C,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 32;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p10r,tmp,0x0d00,cy_r10,cy_i10,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 40;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p14r,tmp,0x0c40,cy_r14,cy_i14,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 48;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p18r,tmp,0x0b80,cy_r18,cy_i18,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 56;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p1cr,tmp,0x0ac0,cy_r1C,cy_i1C,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 64;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p20r,tmp,0x0a00,cy_r20,cy_i20,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 72;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p24r,tmp,0x0940,cy_r24,cy_i24,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 80;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p28r,tmp,0x0880,cy_r28,cy_i28,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 88;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p2cr,tmp,0x07c0,cy_r2C,cy_i2C,half_arr,sign_mask);
				tmp = base_negacyclic_root+ 96;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p30r,tmp,0x0700,cy_r30,cy_i30,half_arr,sign_mask);
				tmp = base_negacyclic_root+104;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p34r,tmp,0x0640,cy_r34,cy_i34,half_arr,sign_mask);
				tmp = base_negacyclic_root+112;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p38r,tmp,0x0580,cy_r38,cy_i38,half_arr,sign_mask);
				tmp = base_negacyclic_root+120;	SSE2_fermat_carry_norm_pow2_errcheck_X4(s1p3cr,tmp,0x04c0,cy_r3C,cy_i3C,half_arr,sign_mask);

			#elif defined(USE_SSE2)

				// For a description of the data movement for Fermat-mod carries in SSE2 mode, see radix16_ditN_cy_dif1.c.
				// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.

				tmp = half_arr+2;
				VEC_DBL_INIT(tmp, scale);
				/* Get the needed Nth root of -1: */
				add1 = (double *)&rn0[0];
				add2 = (double *)&rn1[0];

				idx_offset = j;
				idx_incr = NDIVR;

				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p00r,cy_r00,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p02r,cy_r04,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p04r,cy_r08,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p06r,cy_r0C,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p08r,cy_r10,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p0ar,cy_r14,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p0cr,cy_r18,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p0er,cy_r1C,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p10r,cy_r20,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p12r,cy_r24,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p14r,cy_r28,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p16r,cy_r2C,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p18r,cy_r30,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p1ar,cy_r34,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p1cr,cy_r38,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p1er,cy_r3C,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p20r,cy_i00,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p22r,cy_i04,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p24r,cy_i08,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p26r,cy_i0C,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p28r,cy_i10,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p2ar,cy_i14,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p2cr,cy_i18,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p2er,cy_i1C,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p30r,cy_i20,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p32r,cy_i24,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p34r,cy_i28,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p36r,cy_i2C,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p38r,cy_i30,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p3ar,cy_i34,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p3cr,cy_i38,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
				SSE2_fermat_carry_norm_pow2_errcheck_X2(s1p3er,cy_i3C,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);

			#else	// Scalar-double mode:

				fermat_carry_norm_pow2_errcheck(ar_p00,ai_p00,cy_r00,cy_i00,0x00*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p01,ai_p01,cy_r01,cy_i01,0x01*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p02,ai_p02,cy_r02,cy_i02,0x02*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p03,ai_p03,cy_r03,cy_i03,0x03*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p04,ai_p04,cy_r04,cy_i04,0x04*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p05,ai_p05,cy_r05,cy_i05,0x05*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p06,ai_p06,cy_r06,cy_i06,0x06*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p07,ai_p07,cy_r07,cy_i07,0x07*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p08,ai_p08,cy_r08,cy_i08,0x08*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p09,ai_p09,cy_r09,cy_i09,0x09*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p0a,ai_p0a,cy_r0A,cy_i0A,0x0A*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p0b,ai_p0b,cy_r0B,cy_i0B,0x0B*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p0c,ai_p0c,cy_r0C,cy_i0C,0x0C*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p0d,ai_p0d,cy_r0D,cy_i0D,0x0D*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p0e,ai_p0e,cy_r0E,cy_i0E,0x0E*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p0f,ai_p0f,cy_r0F,cy_i0F,0x0F*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p10,ai_p10,cy_r10,cy_i10,0x10*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p11,ai_p11,cy_r11,cy_i11,0x11*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p12,ai_p12,cy_r12,cy_i12,0x12*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p13,ai_p13,cy_r13,cy_i13,0x13*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p14,ai_p14,cy_r14,cy_i14,0x14*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p15,ai_p15,cy_r15,cy_i15,0x15*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p16,ai_p16,cy_r16,cy_i16,0x16*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p17,ai_p17,cy_r17,cy_i17,0x17*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p18,ai_p18,cy_r18,cy_i18,0x18*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p19,ai_p19,cy_r19,cy_i19,0x19*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p1a,ai_p1a,cy_r1A,cy_i1A,0x1A*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p1b,ai_p1b,cy_r1B,cy_i1B,0x1B*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p1c,ai_p1c,cy_r1C,cy_i1C,0x1C*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p1d,ai_p1d,cy_r1D,cy_i1D,0x1D*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p1e,ai_p1e,cy_r1E,cy_i1E,0x1E*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p1f,ai_p1f,cy_r1F,cy_i1F,0x1F*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p20,ai_p20,cy_r20,cy_i20,0x20*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p21,ai_p21,cy_r21,cy_i21,0x21*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p22,ai_p22,cy_r22,cy_i22,0x22*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p23,ai_p23,cy_r23,cy_i23,0x23*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p24,ai_p24,cy_r24,cy_i24,0x24*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p25,ai_p25,cy_r25,cy_i25,0x25*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p26,ai_p26,cy_r26,cy_i26,0x26*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p27,ai_p27,cy_r27,cy_i27,0x27*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p28,ai_p28,cy_r28,cy_i28,0x28*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p29,ai_p29,cy_r29,cy_i29,0x29*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p2a,ai_p2a,cy_r2A,cy_i2A,0x2A*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p2b,ai_p2b,cy_r2B,cy_i2B,0x2B*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p2c,ai_p2c,cy_r2C,cy_i2C,0x2C*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p2d,ai_p2d,cy_r2D,cy_i2D,0x2D*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p2e,ai_p2e,cy_r2E,cy_i2E,0x2E*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p2f,ai_p2f,cy_r2F,cy_i2F,0x2F*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p30,ai_p30,cy_r30,cy_i30,0x30*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p31,ai_p31,cy_r31,cy_i31,0x31*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p32,ai_p32,cy_r32,cy_i32,0x32*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p33,ai_p33,cy_r33,cy_i33,0x33*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p34,ai_p34,cy_r34,cy_i34,0x34*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p35,ai_p35,cy_r35,cy_i35,0x35*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p36,ai_p36,cy_r36,cy_i36,0x36*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p37,ai_p37,cy_r37,cy_i37,0x37*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p38,ai_p38,cy_r38,cy_i38,0x38*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p39,ai_p39,cy_r39,cy_i39,0x39*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p3a,ai_p3a,cy_r3A,cy_i3A,0x3A*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p3b,ai_p3b,cy_r3B,cy_i3B,0x3B*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p3c,ai_p3c,cy_r3C,cy_i3C,0x3C*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p3d,ai_p3d,cy_r3D,cy_i3D,0x3D*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p3e,ai_p3e,cy_r3E,cy_i3E,0x3E*NDIVR,NRTM1,NRT_BITS);
				fermat_carry_norm_pow2_errcheck(ar_p3f,ai_p3f,cy_r3F,cy_i3F,0x3F*NDIVR,NRTM1,NRT_BITS);

			#endif	/* #ifdef USE_SSE2 */

			}	/* if(MODULUS_TYPE == ...) */

			/*...The radix-64 DIF pass is here:	*/

			#ifdef USE_AVX

			/* Gather the needed data (64 64-bit complex, i.e. 128 64-bit reals) and do 8 twiddleless length-8 subtransforms: */
				//...Block 0: jt = j1;	jp = j2;
	// Relative to RADIX_08_DIF_OOP, this SIMD macro produces outputs in BR order [04261537], so swap r-pointer pairs 2/8,6/C to handle that:
				SSE2_RADIX8_DIF_0TWIDDLE(
					s1p00r,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00,
					r00,r08,r04,r0C,r02,r0A,r06,r0E, isrt2
				);
				//...Block 1: jt = j1 + p04;	jp = j2 + p04;
				SSE2_RADIX8_DIF_0TWIDDLE(
					s1p04r,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00,
					r10,r18,r14,r1C,r12,r1A,r16,r1E, isrt2
				);
				//...Block 2: jt = j1 + p02;	jp = j2 + p02;
				SSE2_RADIX8_DIF_0TWIDDLE(
					s1p02r,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00,
					r20,r28,r24,r2C,r22,r2A,r26,r2E, isrt2
				);
				//...Block 3: jt = j1 + p06;	jp = j2 + p06;
				SSE2_RADIX8_DIF_0TWIDDLE(
					s1p06r,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00,
					r30,r38,r34,r3C,r32,r3A,r36,r3E, isrt2
				);
				//...Block 4: jt = j1 + p01;	jp = j2 + p01;
				SSE2_RADIX8_DIF_0TWIDDLE(
					s1p01r,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00,
					r40,r48,r44,r4C,r42,r4A,r46,r4E, isrt2
				);
				//...Block 5: jt = j1 + p05;	jp = j2 + p05;
				SSE2_RADIX8_DIF_0TWIDDLE(
					s1p05r,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00,
					r50,r58,r54,r5C,r52,r5A,r56,r5E, isrt2
				);
				//...Block 6: jt = j1 + p03;	jp = j2 + p03;
				SSE2_RADIX8_DIF_0TWIDDLE(
					s1p03r,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00,
					r60,r68,r64,r6C,r62,r6A,r66,r6E, isrt2
				);
				//...Block 7: jt = j1 + p07;	jp = j2 + p07;
				SSE2_RADIX8_DIF_0TWIDDLE(
					s1p07r,0x200,0x400,0x600,0x800,0xa00,0xc00,0xe00,
					r70,r78,r74,r7C,r72,r7A,r76,r7E, isrt2
				);

			/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dif_pass1 for details: */

				/* Block 0: */
				add0 = &a[j1]      ; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				/* 0-index block has all-unity twiddles: Remember, the twiddleless DIF bit-reverses both its in-and-outputs,
				so swap index-offset pairs 1/4 and 3/6 in t*-inputs and a-outputs: */
				SSE2_RADIX8_DIF_0TWIDDLE(
					r00, 0x800,0x400,0xc00,0x200,0xa00,0x600,0xe00,
					add0,add1,add2,add3,add4,add5,add6,add7, isrt2
				);
				/* Block 4: */
				add0 = &a[j1] + p08; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIF_TWIDDLE_OOP(
					r08,r48,r28,r68,r18,r58,r38,r78,
					add0,add1,add2,add3,add4,add5,add6,add7,
					ss0,cc0, isrt2,isrt2, nisrt2,isrt2, cc4,ss4, nss4,cc4, ss4,cc4, ncc4,ss4
				);
				/* Block 2: */
				add0 = &a[j1] + p10; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIF_TWIDDLE_OOP(
					r04,r44,r24,r64,r14,r54,r34,r74,
					add0,add1,add2,add3,add4,add5,add6,add7,
					isrt2,isrt2, cc4,ss4, ss4,cc4, cc2,ss2, ss6,cc6, cc6,ss6, ss2,cc2
				);
				/* Block 6: */
				add0 = &a[j1] + p18; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIF_TWIDDLE_OOP(
					r0C,r4C,r2C,r6C,r1C,r5C,r3C,r7C,
					add0,add1,add2,add3,add4,add5,add6,add7,
					nisrt2,isrt2, ss4,cc4, ncc4,nss4, cc6,ss6, ncc2,ss2 ,nss2,cc2, nss6,ncc6
				);
				/* Block 1: */
				add0 = &a[j1] + p20; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIF_TWIDDLE_OOP(
					r02,r42,r22,r62,r12,r52,r32,r72,
					add0,add1,add2,add3,add4,add5,add6,add7,
					cc4,ss4, cc2,ss2, cc6,ss6, cc1,ss1, cc5,ss5, cc3,ss3, cc7,ss7
				);
				/* Block 5: */
				add0 = &a[j1] + p28; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIF_TWIDDLE_OOP(
					r0A,r4A,r2A,r6A,r1A,r5A,r3A,r7A,
					add0,add1,add2,add3,add4,add5,add6,add7,
					nss4,cc4, ss6,cc6, ncc2,ss2, cc5,ss5, ncc7,ss7, ss1,cc1, ncc3,nss3
				);
				/* Block 3: */
				add0 = &a[j1] + p30; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIF_TWIDDLE_OOP(
					r06,r46,r26,r66,r16,r56,r36,r76,
					add0,add1,add2,add3,add4,add5,add6,add7,
					ss4,cc4, cc6,ss6, nss2,cc2, cc3,ss3, ss1,cc1, ss7,cc7, nss5,cc5
				);
				/* Block 7: */
				add0 = &a[j1] + p38; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIF_TWIDDLE_OOP(
					r0E,r4E,r2E,r6E,r1E,r5E,r3E,r7E,
					add0,add1,add2,add3,add4,add5,add6,add7,
					ncc4,ss4, ss2,cc2, nss6,ncc6, cc7,ss7, ncc3,nss3, nss5,cc5, ss1,ncc1
				);

			#elif defined(USE_SSE2)

			/* Gather the needed data (64 64-bit complex, i.e. 128 64-bit reals) and do 8 twiddleless length-8 subtransforms: */
				//...Block 0: jt = j1;	jp = j2;
	// Relative to RADIX_08_DIF_OOP, this SIMD macro produces outputs in BR order [04261537], so swap r-pointer pairs 2/8,6/C to handle that:
				SSE2_RADIX8_DIF_0TWIDDLE(
					s1p00r,0x100,0x200,0x300,0x400,0x500,0x600,0x700,
					r00,r08,r04,r0C,r02,r0A,r06,r0E, isrt2
				);
				//...Block 1: jt = j1 + p04;	jp = j2 + p04;
				SSE2_RADIX8_DIF_0TWIDDLE(
					s1p04r,0x100,0x200,0x300,0x400,0x500,0x600,0x700,
					r10,r18,r14,r1C,r12,r1A,r16,r1E, isrt2
				);
				//...Block 2: jt = j1 + p02;	jp = j2 + p02;
				SSE2_RADIX8_DIF_0TWIDDLE(
					s1p02r,0x100,0x200,0x300,0x400,0x500,0x600,0x700,
					r20,r28,r24,r2C,r22,r2A,r26,r2E, isrt2
				);
				//...Block 3: jt = j1 + p06;	jp = j2 + p06;
				SSE2_RADIX8_DIF_0TWIDDLE(
					s1p06r,0x100,0x200,0x300,0x400,0x500,0x600,0x700,
					r30,r38,r34,r3C,r32,r3A,r36,r3E, isrt2
				);
				//...Block 4: jt = j1 + p01;	jp = j2 + p01;
				SSE2_RADIX8_DIF_0TWIDDLE(
					s1p01r,0x100,0x200,0x300,0x400,0x500,0x600,0x700,
					r40,r48,r44,r4C,r42,r4A,r46,r4E, isrt2
				);
				//...Block 5: jt = j1 + p05;	jp = j2 + p05;
				SSE2_RADIX8_DIF_0TWIDDLE(
					s1p05r,0x100,0x200,0x300,0x400,0x500,0x600,0x700,
					r50,r58,r54,r5C,r52,r5A,r56,r5E, isrt2
				);
				//...Block 6: jt = j1 + p03;	jp = j2 + p03;
				SSE2_RADIX8_DIF_0TWIDDLE(
					s1p03r,0x100,0x200,0x300,0x400,0x500,0x600,0x700,
					r60,r68,r64,r6C,r62,r6A,r66,r6E, isrt2
				);
				//...Block 7: jt = j1 + p07;	jp = j2 + p07;
				SSE2_RADIX8_DIF_0TWIDDLE(
					s1p07r,0x100,0x200,0x300,0x400,0x500,0x600,0x700,
					r70,r78,r74,r7C,r72,r7A,r76,r7E, isrt2
				);

			/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dif_pass1 for details: */

				/* Block 0: */
				add0 = &a[j1]      ; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				/* 0-index block has all-unity twiddles: Remember, the twiddleless DIF bit-reverses both its in-and-outputs,
				so swap index-offset pairs 1/4 and 3/6 in t*-inputs and a-outputs: */
				SSE2_RADIX8_DIF_0TWIDDLE(
					r00, 0x400,0x200,0x600,0x100,0x500,0x300,0x700,
					add0,add1,add2,add3,add4,add5,add6,add7, isrt2
				);
				/* Block 4: */
				add0 = &a[j1] + p08; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIF_TWIDDLE_OOP(
					r08,r48,r28,r68,r18,r58,r38,r78,
					add0,add1,add2,add3,add4,add5,add6,add7,
					ss0,cc0, isrt2,isrt2, nisrt2,isrt2, cc4,ss4, nss4,cc4, ss4,cc4, ncc4,ss4
				);
				/* Block 2: */
				add0 = &a[j1] + p10; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIF_TWIDDLE_OOP(
					r04,r44,r24,r64,r14,r54,r34,r74,
					add0,add1,add2,add3,add4,add5,add6,add7,
					isrt2,isrt2, cc4,ss4, ss4,cc4, cc2,ss2, ss6,cc6, cc6,ss6, ss2,cc2
				);
				/* Block 6: */
				add0 = &a[j1] + p18; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIF_TWIDDLE_OOP(
					r0C,r4C,r2C,r6C,r1C,r5C,r3C,r7C,
					add0,add1,add2,add3,add4,add5,add6,add7,
					nisrt2,isrt2, ss4,cc4, ncc4,nss4, cc6,ss6, ncc2,ss2 ,nss2,cc2, nss6,ncc6
				);
				/* Block 1: */
				add0 = &a[j1] + p20; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIF_TWIDDLE_OOP(
					r02,r42,r22,r62,r12,r52,r32,r72,
					add0,add1,add2,add3,add4,add5,add6,add7,
					cc4,ss4, cc2,ss2, cc6,ss6, cc1,ss1, cc5,ss5, cc3,ss3, cc7,ss7
				);
				/* Block 5: */
				add0 = &a[j1] + p28; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIF_TWIDDLE_OOP(
					r0A,r4A,r2A,r6A,r1A,r5A,r3A,r7A,
					add0,add1,add2,add3,add4,add5,add6,add7,
					nss4,cc4, ss6,cc6, ncc2,ss2, cc5,ss5, ncc7,ss7, ss1,cc1, ncc3,nss3
				);
				/* Block 3: */
				add0 = &a[j1] + p30; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIF_TWIDDLE_OOP(
					r06,r46,r26,r66,r16,r56,r36,r76,
					add0,add1,add2,add3,add4,add5,add6,add7,
					ss4,cc4, cc6,ss6, nss2,cc2, cc3,ss3, ss1,cc1, ss7,cc7, nss5,cc5
				);
				/* Block 7: */
				add0 = &a[j1] + p38; add1 = add0+p01; add2 = add0+p02; add3 = add0+p03; add4 = add0+p04; add5 = add0+p05; add6 = add0+p06; add7 = add0+p07;
				SSE2_RADIX8_DIF_TWIDDLE_OOP(
					r0E,r4E,r2E,r6E,r1E,r5E,r3E,r7E,
					add0,add1,add2,add3,add4,add5,add6,add7,
					ncc4,ss4, ss2,cc2, nss6,ncc6, cc7,ss7, ncc3,nss3, nss5,cc5, ss1,ncc1
				);

			#else	// USE_SSE2 = False:

			/* Gather the needed data (64 64-bit complex, i.e. 128 64-bit reals) and do 8 twiddleless length-8 subtransforms: */
				//...Block 0: jt = j1;	jp = j2;
				RADIX_08_DIF_OOP(
					ar_p00,ai_p00,ar_p08,ai_p08,ar_p10,ai_p10,ar_p18,ai_p18,ar_p20,ai_p20,ar_p28,ai_p28,ar_p30,ai_p30,ar_p38,ai_p38,
					t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F
				);
				//...Block 1: jt = j1 + p04;	jp = j2 + p04;
				RADIX_08_DIF_OOP(
					ar_p04,ai_p04,ar_p0c,ai_p0c,ar_p14,ai_p14,ar_p1c,ai_p1c,ar_p24,ai_p24,ar_p2c,ai_p2c,ar_p34,ai_p34,ar_p3c,ai_p3c,
					t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F
				);
				//...Block 2: jt = j1 + p02;	jp = j2 + p02;
				RADIX_08_DIF_OOP(
					ar_p02,ai_p02,ar_p0a,ai_p0a,ar_p12,ai_p12,ar_p1a,ai_p1a,ar_p22,ai_p22,ar_p2a,ai_p2a,ar_p32,ai_p32,ar_p3a,ai_p3a,
					t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F
				);
				//...Block 3: jt = j1 + p06;	jp = j2 + p06;
				RADIX_08_DIF_OOP(
					ar_p06,ai_p06,ar_p0e,ai_p0e,ar_p16,ai_p16,ar_p1e,ai_p1e,ar_p26,ai_p26,ar_p2e,ai_p2e,ar_p36,ai_p36,ar_p3e,ai_p3e,
					t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F
				);
				//...Block 4: jt = j1 + p01;	jp = j2 + p01;
				RADIX_08_DIF_OOP(
					ar_p01,ai_p01,ar_p09,ai_p09,ar_p11,ai_p11,ar_p19,ai_p19,ar_p21,ai_p21,ar_p29,ai_p29,ar_p31,ai_p31,ar_p39,ai_p39,
					t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t4A,t4B,t4C,t4D,t4E,t4F
				);
				//...Block 5: jt = j1 + p05;	jp = j2 + p05;
				RADIX_08_DIF_OOP(
					ar_p05,ai_p05,ar_p0d,ai_p0d,ar_p15,ai_p15,ar_p1d,ai_p1d,ar_p25,ai_p25,ar_p2d,ai_p2d,ar_p35,ai_p35,ar_p3d,ai_p3d,
					t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t5A,t5B,t5C,t5D,t5E,t5F
				);
				//...Block 6: jt = j1 + p03;	jp = j2 + p03;
				RADIX_08_DIF_OOP(
					ar_p03,ai_p03,ar_p0b,ai_p0b,ar_p13,ai_p13,ar_p1b,ai_p1b,ar_p23,ai_p23,ar_p2b,ai_p2b,ar_p33,ai_p33,ar_p3b,ai_p3b,
					t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t6A,t6B,t6C,t6D,t6E,t6F
				);
				//...Block 7: jt = j1 + p07;	jp = j2 + p07;
				RADIX_08_DIF_OOP(
					ar_p07,ai_p07,ar_p0f,ai_p0f,ar_p17,ai_p17,ar_p1f,ai_p1f,ar_p27,ai_p27,ar_p2f,ai_p2f,ar_p37,ai_p37,ar_p3f,ai_p3f,
					t70,t71,t72,t73,t74,t75,t76,t77,t78,t79,t7A,t7B,t7C,t7D,t7E,t7F
				);

			/*...and now do eight radix-8 subtransforms w/internal twiddles - cf. radix64_dif_pass1 for details: */

				/* Block 0: */
				jt = j1;	jp = j2;
				/* 0-index block has all-unity twiddles: Remember, the twiddleless DIF bit-reverses both its in-and-outputs,
				so swap index-offset pairs 1/4 and 3/6 in t*-inputs and a-outputs: */
				RADIX_08_DIF_OOP(
					t00,t01,t40,t41,t20,t21,t60,t61,t10,t11,t50,t51,t30,t31,t70,t71,
					a[jt],a[jp],a[jt+p04],a[jp+p04],a[jt+p02],a[jp+p02],a[jt+p06],a[jp+p06],a[jt+p01],a[jp+p01],a[jt+p05],a[jp+p05],a[jt+p03],a[jp+p03],a[jt+p07],a[jp+p07]
				);
				/* Block 4: */
				jt = j1 + p08;	jp = j2 + p08;
				RADIX_08_DIF_TWIDDLE_OOP(
					t08,t09,t18,t19,t28,t29,t38,t39,t48,t49,t58,t59,t68,t69,t78,t79,
					a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
					0,1,ISRT2,ISRT2,-ISRT2,ISRT2,c64_4,s64_4,-s64_4,c64_4,s64_4,c64_4,-c64_4,s64_4
				);
				/* Block 2: */
				jt = j1 + p10;	jp = j2 + p10;
				RADIX_08_DIF_TWIDDLE_OOP(
					t04,t05,t14,t15,t24,t25,t34,t35,t44,t45,t54,t55,t64,t65,t74,t75,
					a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
					ISRT2,ISRT2,c64_4,s64_4,s64_4,c64_4,c64_2,s64_2,s64_6,c64_6,c64_6,s64_6,s64_2,c64_2
				);
				/* Block 6: */
				jt = j1 + p18;	jp = j2 + p18;
				RADIX_08_DIF_TWIDDLE_OOP(
					t0C,t0D,t1C,t1D,t2C,t2D,t3C,t3D,t4C,t4D,t5C,t5D,t6C,t6D,t7C,t7D,
					a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
					-ISRT2,ISRT2,s64_4,c64_4,-c64_4,-s64_4,c64_6,s64_6,-c64_2,s64_2,-s64_2,c64_2,-s64_6,-c64_6
				);
				/* Block 1: */
				jt = j1 + p20;	jp = j2 + p20;
				RADIX_08_DIF_TWIDDLE_OOP(
					t02,t03,t12,t13,t22,t23,t32,t33,t42,t43,t52,t53,t62,t63,t72,t73,
					a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
					c64_4,s64_4,c64_2,s64_2,c64_6,s64_6,c64_1,s64_1,c64_5,s64_5,c64_3,s64_3,c64_7,s64_7
				);
				/* Block 5: */
				jt = j1 + p28;	jp = j2 + p28;
				RADIX_08_DIF_TWIDDLE_OOP(
					t0A,t0B,t1A,t1B,t2A,t2B,t3A,t3B,t4A,t4B,t5A,t5B,t6A,t6B,t7A,t7B,
					a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
					-s64_4,c64_4,s64_6,c64_6,-c64_2,s64_2,c64_5,s64_5,-c64_7,s64_7,s64_1,c64_1,-c64_3,-s64_3
				);
				/* Block 3: */
				jt = j1 + p30;	jp = j2 + p30;
				RADIX_08_DIF_TWIDDLE_OOP(
					t06,t07,t16,t17,t26,t27,t36,t37,t46,t47,t56,t57,t66,t67,t76,t77,
					a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
					s64_4,c64_4,c64_6,s64_6,-s64_2,c64_2,c64_3,s64_3,s64_1,c64_1,s64_7,c64_7,-s64_5,c64_5
				);
				/* Block 7: */
				jt = j1 + p38;	jp = j2 + p38;
				RADIX_08_DIF_TWIDDLE_OOP(
					t0E,t0F,t1E,t1F,t2E,t2F,t3E,t3F,t4E,t4F,t5E,t5F,t6E,t6F,t7E,t7F,
					a[jt],a[jp],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],
					-c64_4,s64_4,s64_2,c64_2,-s64_6,-c64_6,c64_7,s64_7,-c64_3,-s64_3,-s64_5,c64_5,s64_1,-c64_1
				);

			#endif	/* #ifdef USE_SSE2 */

			}	/* end for(j=_jstart; j < _jhi; j += 2) */

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
				jstart += nwt;
				jhi    += nwt;

				col += RADIX;
				co3 -= RADIX;
			}
		}	/* end for(k=1; k <= khi; k++) */

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
			thread_arg->cy_r0A = cy_r08->d2;
			thread_arg->cy_r0B = cy_r08->d3;
			thread_arg->cy_r0C = cy_r0C->d0;
			thread_arg->cy_r0D = cy_r0C->d1;
			thread_arg->cy_r0E = cy_r0C->d2;
			thread_arg->cy_r0F = cy_r0C->d3;
			thread_arg->cy_r10 = cy_r10->d0;
			thread_arg->cy_r11 = cy_r10->d1;
			thread_arg->cy_r12 = cy_r10->d2;
			thread_arg->cy_r13 = cy_r10->d3;
			thread_arg->cy_r14 = cy_r14->d0;
			thread_arg->cy_r15 = cy_r14->d1;
			thread_arg->cy_r16 = cy_r14->d2;
			thread_arg->cy_r17 = cy_r14->d3;
			thread_arg->cy_r18 = cy_r18->d0;
			thread_arg->cy_r19 = cy_r18->d1;
			thread_arg->cy_r1A = cy_r18->d2;
			thread_arg->cy_r1B = cy_r18->d3;
			thread_arg->cy_r1C = cy_r1C->d0;
			thread_arg->cy_r1D = cy_r1C->d1;
			thread_arg->cy_r1E = cy_r1C->d2;
			thread_arg->cy_r1F = cy_r1C->d3;
			thread_arg->cy_r20 = cy_r20->d0;
			thread_arg->cy_r21 = cy_r20->d1;
			thread_arg->cy_r22 = cy_r20->d2;
			thread_arg->cy_r23 = cy_r20->d3;
			thread_arg->cy_r24 = cy_r24->d0;
			thread_arg->cy_r25 = cy_r24->d1;
			thread_arg->cy_r26 = cy_r24->d2;
			thread_arg->cy_r27 = cy_r24->d3;
			thread_arg->cy_r28 = cy_r28->d0;
			thread_arg->cy_r29 = cy_r28->d1;
			thread_arg->cy_r2A = cy_r28->d2;
			thread_arg->cy_r2B = cy_r28->d3;
			thread_arg->cy_r2C = cy_r2C->d0;
			thread_arg->cy_r2D = cy_r2C->d1;
			thread_arg->cy_r2E = cy_r2C->d2;
			thread_arg->cy_r2F = cy_r2C->d3;
			thread_arg->cy_r30 = cy_r30->d0;
			thread_arg->cy_r31 = cy_r30->d1;
			thread_arg->cy_r32 = cy_r30->d2;
			thread_arg->cy_r33 = cy_r30->d3;
			thread_arg->cy_r34 = cy_r34->d0;
			thread_arg->cy_r35 = cy_r34->d1;
			thread_arg->cy_r36 = cy_r34->d2;
			thread_arg->cy_r37 = cy_r34->d3;
			thread_arg->cy_r38 = cy_r38->d0;
			thread_arg->cy_r39 = cy_r38->d1;
			thread_arg->cy_r3A = cy_r38->d2;
			thread_arg->cy_r3B = cy_r38->d3;
			thread_arg->cy_r3C = cy_r3C->d0;
			thread_arg->cy_r3D = cy_r3C->d1;
			thread_arg->cy_r3E = cy_r3C->d2;
			thread_arg->cy_r3F = cy_r3C->d3;
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
			thread_arg->cy_r0A = cy_r0A->d0;
			thread_arg->cy_r0B = cy_r0A->d1;
			thread_arg->cy_r0C = cy_r0C->d0;
			thread_arg->cy_r0D = cy_r0C->d1;
			thread_arg->cy_r0E = cy_r0E->d0;
			thread_arg->cy_r0F = cy_r0E->d1;
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
			thread_arg->cy_r1A = cy_r1A->d0;
			thread_arg->cy_r1B = cy_r1A->d1;
			thread_arg->cy_r1C = cy_r1C->d0;
			thread_arg->cy_r1D = cy_r1C->d1;
			thread_arg->cy_r1E = cy_r1E->d0;
			thread_arg->cy_r1F = cy_r1E->d1;
			thread_arg->cy_r20 = cy_r20->d0;
			thread_arg->cy_r21 = cy_r20->d1;
			thread_arg->cy_r22 = cy_r22->d0;
			thread_arg->cy_r23 = cy_r22->d1;
			thread_arg->cy_r24 = cy_r24->d0;
			thread_arg->cy_r25 = cy_r24->d1;
			thread_arg->cy_r26 = cy_r26->d0;
			thread_arg->cy_r27 = cy_r26->d1;
			thread_arg->cy_r28 = cy_r28->d0;
			thread_arg->cy_r29 = cy_r28->d1;
			thread_arg->cy_r2A = cy_r2A->d0;
			thread_arg->cy_r2B = cy_r2A->d1;
			thread_arg->cy_r2C = cy_r2C->d0;
			thread_arg->cy_r2D = cy_r2C->d1;
			thread_arg->cy_r2E = cy_r2E->d0;
			thread_arg->cy_r2F = cy_r2E->d1;
			thread_arg->cy_r30 = cy_r30->d0;
			thread_arg->cy_r31 = cy_r30->d1;
			thread_arg->cy_r32 = cy_r32->d0;
			thread_arg->cy_r33 = cy_r32->d1;
			thread_arg->cy_r34 = cy_r34->d0;
			thread_arg->cy_r35 = cy_r34->d1;
			thread_arg->cy_r36 = cy_r36->d0;
			thread_arg->cy_r37 = cy_r36->d1;
			thread_arg->cy_r38 = cy_r38->d0;
			thread_arg->cy_r39 = cy_r38->d1;
			thread_arg->cy_r3A = cy_r3A->d0;
			thread_arg->cy_r3B = cy_r3A->d1;
			thread_arg->cy_r3C = cy_r3C->d0;
			thread_arg->cy_r3D = cy_r3C->d1;
			thread_arg->cy_r3E = cy_r3E->d0;
			thread_arg->cy_r3F = cy_r3E->d1;
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
			thread_arg->cy_r0A = cy_r0A;
			thread_arg->cy_r0B = cy_r0B;
			thread_arg->cy_r0C = cy_r0C;
			thread_arg->cy_r0D = cy_r0D;
			thread_arg->cy_r0E = cy_r0E;
			thread_arg->cy_r0F = cy_r0F;
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
			thread_arg->cy_r1A = cy_r1A;
			thread_arg->cy_r1B = cy_r1B;
			thread_arg->cy_r1C = cy_r1C;
			thread_arg->cy_r1D = cy_r1D;
			thread_arg->cy_r1E = cy_r1E;
			thread_arg->cy_r1F = cy_r1F;
			thread_arg->cy_r20 = cy_r20;
			thread_arg->cy_r21 = cy_r21;
			thread_arg->cy_r22 = cy_r22;
			thread_arg->cy_r23 = cy_r23;
			thread_arg->cy_r24 = cy_r24;
			thread_arg->cy_r25 = cy_r25;
			thread_arg->cy_r26 = cy_r26;
			thread_arg->cy_r27 = cy_r27;
			thread_arg->cy_r28 = cy_r28;
			thread_arg->cy_r29 = cy_r29;
			thread_arg->cy_r2A = cy_r2A;
			thread_arg->cy_r2B = cy_r2B;
			thread_arg->cy_r2C = cy_r2C;
			thread_arg->cy_r2D = cy_r2D;
			thread_arg->cy_r2E = cy_r2E;
			thread_arg->cy_r2F = cy_r2F;
			thread_arg->cy_r30 = cy_r30;
			thread_arg->cy_r31 = cy_r31;
			thread_arg->cy_r32 = cy_r32;
			thread_arg->cy_r33 = cy_r33;
			thread_arg->cy_r34 = cy_r34;
			thread_arg->cy_r35 = cy_r35;
			thread_arg->cy_r36 = cy_r36;
			thread_arg->cy_r37 = cy_r37;
			thread_arg->cy_r38 = cy_r38;
			thread_arg->cy_r39 = cy_r39;
			thread_arg->cy_r3A = cy_r3A;
			thread_arg->cy_r3B = cy_r3B;
			thread_arg->cy_r3C = cy_r3C;
			thread_arg->cy_r3D = cy_r3D;
			thread_arg->cy_r3E = cy_r3E;
			thread_arg->cy_r3F = cy_r3F;

		#endif	// SSE2 or AVX?
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
			thread_arg->cy_r0A = cy_r08->d2;	thread_arg->cy_i0A = cy_i08->d2;
			thread_arg->cy_r0B = cy_r08->d3;	thread_arg->cy_i0B = cy_i08->d3;
			thread_arg->cy_r0C = cy_r0C->d0;	thread_arg->cy_i0C = cy_i0C->d0;
			thread_arg->cy_r0D = cy_r0C->d1;	thread_arg->cy_i0D = cy_i0C->d1;
			thread_arg->cy_r0E = cy_r0C->d2;	thread_arg->cy_i0E = cy_i0C->d2;
			thread_arg->cy_r0F = cy_r0C->d3;	thread_arg->cy_i0F = cy_i0C->d3;
			thread_arg->cy_r10 = cy_r10->d0;	thread_arg->cy_i10 = cy_i10->d0;
			thread_arg->cy_r11 = cy_r10->d1;	thread_arg->cy_i11 = cy_i10->d1;
			thread_arg->cy_r12 = cy_r10->d2;	thread_arg->cy_i12 = cy_i10->d2;
			thread_arg->cy_r13 = cy_r10->d3;	thread_arg->cy_i13 = cy_i10->d3;
			thread_arg->cy_r14 = cy_r14->d0;	thread_arg->cy_i14 = cy_i14->d0;
			thread_arg->cy_r15 = cy_r14->d1;	thread_arg->cy_i15 = cy_i14->d1;
			thread_arg->cy_r16 = cy_r14->d2;	thread_arg->cy_i16 = cy_i14->d2;
			thread_arg->cy_r17 = cy_r14->d3;	thread_arg->cy_i17 = cy_i14->d3;
			thread_arg->cy_r18 = cy_r18->d0;	thread_arg->cy_i18 = cy_i18->d0;
			thread_arg->cy_r19 = cy_r18->d1;	thread_arg->cy_i19 = cy_i18->d1;
			thread_arg->cy_r1A = cy_r18->d2;	thread_arg->cy_i1A = cy_i18->d2;
			thread_arg->cy_r1B = cy_r18->d3;	thread_arg->cy_i1B = cy_i18->d3;
			thread_arg->cy_r1C = cy_r1C->d0;	thread_arg->cy_i1C = cy_i1C->d0;
			thread_arg->cy_r1D = cy_r1C->d1;	thread_arg->cy_i1D = cy_i1C->d1;
			thread_arg->cy_r1E = cy_r1C->d2;	thread_arg->cy_i1E = cy_i1C->d2;
			thread_arg->cy_r1F = cy_r1C->d3;	thread_arg->cy_i1F = cy_i1C->d3;
			thread_arg->cy_r20 = cy_r20->d0;	thread_arg->cy_i20 = cy_i20->d0;
			thread_arg->cy_r21 = cy_r20->d1;	thread_arg->cy_i21 = cy_i20->d1;
			thread_arg->cy_r22 = cy_r20->d2;	thread_arg->cy_i22 = cy_i20->d2;
			thread_arg->cy_r23 = cy_r20->d3;	thread_arg->cy_i23 = cy_i20->d3;
			thread_arg->cy_r24 = cy_r24->d0;	thread_arg->cy_i24 = cy_i24->d0;
			thread_arg->cy_r25 = cy_r24->d1;	thread_arg->cy_i25 = cy_i24->d1;
			thread_arg->cy_r26 = cy_r24->d2;	thread_arg->cy_i26 = cy_i24->d2;
			thread_arg->cy_r27 = cy_r24->d3;	thread_arg->cy_i27 = cy_i24->d3;
			thread_arg->cy_r28 = cy_r28->d0;	thread_arg->cy_i28 = cy_i28->d0;
			thread_arg->cy_r29 = cy_r28->d1;	thread_arg->cy_i29 = cy_i28->d1;
			thread_arg->cy_r2A = cy_r28->d2;	thread_arg->cy_i2A = cy_i28->d2;
			thread_arg->cy_r2B = cy_r28->d3;	thread_arg->cy_i2B = cy_i28->d3;
			thread_arg->cy_r2C = cy_r2C->d0;	thread_arg->cy_i2C = cy_i2C->d0;
			thread_arg->cy_r2D = cy_r2C->d1;	thread_arg->cy_i2D = cy_i2C->d1;
			thread_arg->cy_r2E = cy_r2C->d2;	thread_arg->cy_i2E = cy_i2C->d2;
			thread_arg->cy_r2F = cy_r2C->d3;	thread_arg->cy_i2F = cy_i2C->d3;
			thread_arg->cy_r30 = cy_r30->d0;	thread_arg->cy_i30 = cy_i30->d0;
			thread_arg->cy_r31 = cy_r30->d1;	thread_arg->cy_i31 = cy_i30->d1;
			thread_arg->cy_r32 = cy_r30->d2;	thread_arg->cy_i32 = cy_i30->d2;
			thread_arg->cy_r33 = cy_r30->d3;	thread_arg->cy_i33 = cy_i30->d3;
			thread_arg->cy_r34 = cy_r34->d0;	thread_arg->cy_i34 = cy_i34->d0;
			thread_arg->cy_r35 = cy_r34->d1;	thread_arg->cy_i35 = cy_i34->d1;
			thread_arg->cy_r36 = cy_r34->d2;	thread_arg->cy_i36 = cy_i34->d2;
			thread_arg->cy_r37 = cy_r34->d3;	thread_arg->cy_i37 = cy_i34->d3;
			thread_arg->cy_r38 = cy_r38->d0;	thread_arg->cy_i38 = cy_i38->d0;
			thread_arg->cy_r39 = cy_r38->d1;	thread_arg->cy_i39 = cy_i38->d1;
			thread_arg->cy_r3A = cy_r38->d2;	thread_arg->cy_i3A = cy_i38->d2;
			thread_arg->cy_r3B = cy_r38->d3;	thread_arg->cy_i3B = cy_i38->d3;
			thread_arg->cy_r3C = cy_r3C->d0;	thread_arg->cy_i3C = cy_i3C->d0;
			thread_arg->cy_r3D = cy_r3C->d1;	thread_arg->cy_i3D = cy_i3C->d1;
			thread_arg->cy_r3E = cy_r3C->d2;	thread_arg->cy_i3E = cy_i3C->d2;
			thread_arg->cy_r3F = cy_r3C->d3;	thread_arg->cy_i3F = cy_i3C->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );

		#elif defined(USE_SSE2)

			thread_arg->cy_r00 = cy_r00->d0;	thread_arg->cy_i00 = cy_r00->d1;
			thread_arg->cy_r01 = cy_r02->d0;	thread_arg->cy_i01 = cy_r02->d1;
			thread_arg->cy_r02 = cy_r04->d0;	thread_arg->cy_i02 = cy_r04->d1;
			thread_arg->cy_r03 = cy_r06->d0;	thread_arg->cy_i03 = cy_r06->d1;
			thread_arg->cy_r04 = cy_r08->d0;	thread_arg->cy_i04 = cy_r08->d1;
			thread_arg->cy_r05 = cy_r0A->d0;	thread_arg->cy_i05 = cy_r0A->d1;
			thread_arg->cy_r06 = cy_r0C->d0;	thread_arg->cy_i06 = cy_r0C->d1;
			thread_arg->cy_r07 = cy_r0E->d0;	thread_arg->cy_i07 = cy_r0E->d1;
			thread_arg->cy_r08 = cy_r10->d0;	thread_arg->cy_i08 = cy_r10->d1;
			thread_arg->cy_r09 = cy_r12->d0;	thread_arg->cy_i09 = cy_r12->d1;
			thread_arg->cy_r0A = cy_r14->d0;	thread_arg->cy_i0A = cy_r14->d1;
			thread_arg->cy_r0B = cy_r16->d0;	thread_arg->cy_i0B = cy_r16->d1;
			thread_arg->cy_r0C = cy_r18->d0;	thread_arg->cy_i0C = cy_r18->d1;
			thread_arg->cy_r0D = cy_r1A->d0;	thread_arg->cy_i0D = cy_r1A->d1;
			thread_arg->cy_r0E = cy_r1C->d0;	thread_arg->cy_i0E = cy_r1C->d1;
			thread_arg->cy_r0F = cy_r1E->d0;	thread_arg->cy_i0F = cy_r1E->d1;
			thread_arg->cy_r10 = cy_r20->d0;	thread_arg->cy_i10 = cy_r20->d1;
			thread_arg->cy_r11 = cy_r22->d0;	thread_arg->cy_i11 = cy_r22->d1;
			thread_arg->cy_r12 = cy_r24->d0;	thread_arg->cy_i12 = cy_r24->d1;
			thread_arg->cy_r13 = cy_r26->d0;	thread_arg->cy_i13 = cy_r26->d1;
			thread_arg->cy_r14 = cy_r28->d0;	thread_arg->cy_i14 = cy_r28->d1;
			thread_arg->cy_r15 = cy_r2A->d0;	thread_arg->cy_i15 = cy_r2A->d1;
			thread_arg->cy_r16 = cy_r2C->d0;	thread_arg->cy_i16 = cy_r2C->d1;
			thread_arg->cy_r17 = cy_r2E->d0;	thread_arg->cy_i17 = cy_r2E->d1;
			thread_arg->cy_r18 = cy_r30->d0;	thread_arg->cy_i18 = cy_r30->d1;
			thread_arg->cy_r19 = cy_r32->d0;	thread_arg->cy_i19 = cy_r32->d1;
			thread_arg->cy_r1A = cy_r34->d0;	thread_arg->cy_i1A = cy_r34->d1;
			thread_arg->cy_r1B = cy_r36->d0;	thread_arg->cy_i1B = cy_r36->d1;
			thread_arg->cy_r1C = cy_r38->d0;	thread_arg->cy_i1C = cy_r38->d1;
			thread_arg->cy_r1D = cy_r3A->d0;	thread_arg->cy_i1D = cy_r3A->d1;
			thread_arg->cy_r1E = cy_r3C->d0;	thread_arg->cy_i1E = cy_r3C->d1;
			thread_arg->cy_r1F = cy_r3E->d0;	thread_arg->cy_i1F = cy_r3E->d1;

			thread_arg->cy_r20 = cy_i00->d0;	thread_arg->cy_i20 = cy_i00->d1;
			thread_arg->cy_r21 = cy_i02->d0;	thread_arg->cy_i21 = cy_i02->d1;
			thread_arg->cy_r22 = cy_i04->d0;	thread_arg->cy_i22 = cy_i04->d1;
			thread_arg->cy_r23 = cy_i06->d0;	thread_arg->cy_i23 = cy_i06->d1;
			thread_arg->cy_r24 = cy_i08->d0;	thread_arg->cy_i24 = cy_i08->d1;
			thread_arg->cy_r25 = cy_i0A->d0;	thread_arg->cy_i25 = cy_i0A->d1;
			thread_arg->cy_r26 = cy_i0C->d0;	thread_arg->cy_i26 = cy_i0C->d1;
			thread_arg->cy_r27 = cy_i0E->d0;	thread_arg->cy_i27 = cy_i0E->d1;
			thread_arg->cy_r28 = cy_i10->d0;	thread_arg->cy_i28 = cy_i10->d1;
			thread_arg->cy_r29 = cy_i12->d0;	thread_arg->cy_i29 = cy_i12->d1;
			thread_arg->cy_r2A = cy_i14->d0;	thread_arg->cy_i2A = cy_i14->d1;
			thread_arg->cy_r2B = cy_i16->d0;	thread_arg->cy_i2B = cy_i16->d1;
			thread_arg->cy_r2C = cy_i18->d0;	thread_arg->cy_i2C = cy_i18->d1;
			thread_arg->cy_r2D = cy_i1A->d0;	thread_arg->cy_i2D = cy_i1A->d1;
			thread_arg->cy_r2E = cy_i1C->d0;	thread_arg->cy_i2E = cy_i1C->d1;
			thread_arg->cy_r2F = cy_i1E->d0;	thread_arg->cy_i2F = cy_i1E->d1;
			thread_arg->cy_r30 = cy_i20->d0;	thread_arg->cy_i30 = cy_i20->d1;
			thread_arg->cy_r31 = cy_i22->d0;	thread_arg->cy_i31 = cy_i22->d1;
			thread_arg->cy_r32 = cy_i24->d0;	thread_arg->cy_i32 = cy_i24->d1;
			thread_arg->cy_r33 = cy_i26->d0;	thread_arg->cy_i33 = cy_i26->d1;
			thread_arg->cy_r34 = cy_i28->d0;	thread_arg->cy_i34 = cy_i28->d1;
			thread_arg->cy_r35 = cy_i2A->d0;	thread_arg->cy_i35 = cy_i2A->d1;
			thread_arg->cy_r36 = cy_i2C->d0;	thread_arg->cy_i36 = cy_i2C->d1;
			thread_arg->cy_r37 = cy_i2E->d0;	thread_arg->cy_i37 = cy_i2E->d1;
			thread_arg->cy_r38 = cy_i30->d0;	thread_arg->cy_i38 = cy_i30->d1;
			thread_arg->cy_r39 = cy_i32->d0;	thread_arg->cy_i39 = cy_i32->d1;
			thread_arg->cy_r3A = cy_i34->d0;	thread_arg->cy_i3A = cy_i34->d1;
			thread_arg->cy_r3B = cy_i36->d0;	thread_arg->cy_i3B = cy_i36->d1;
			thread_arg->cy_r3C = cy_i38->d0;	thread_arg->cy_i3C = cy_i38->d1;
			thread_arg->cy_r3D = cy_i3A->d0;	thread_arg->cy_i3D = cy_i3A->d1;
			thread_arg->cy_r3E = cy_i3C->d0;	thread_arg->cy_i3E = cy_i3C->d1;
			thread_arg->cy_r3F = cy_i3E->d0;	thread_arg->cy_i3F = cy_i3E->d1;
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
			thread_arg->cy_r0A = cy_r0A;		thread_arg->cy_i0A = cy_i0A;
			thread_arg->cy_r0B = cy_r0B;		thread_arg->cy_i0B = cy_i0B;
			thread_arg->cy_r0C = cy_r0C;		thread_arg->cy_i0C = cy_i0C;
			thread_arg->cy_r0D = cy_r0D;		thread_arg->cy_i0D = cy_i0D;
			thread_arg->cy_r0E = cy_r0E;		thread_arg->cy_i0E = cy_i0E;
			thread_arg->cy_r0F = cy_r0F;		thread_arg->cy_i0F = cy_i0F;
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
			thread_arg->cy_r1A = cy_r1A;		thread_arg->cy_i1A = cy_i1A;
			thread_arg->cy_r1B = cy_r1B;		thread_arg->cy_i1B = cy_i1B;
			thread_arg->cy_r1C = cy_r1C;		thread_arg->cy_i1C = cy_i1C;
			thread_arg->cy_r1D = cy_r1D;		thread_arg->cy_i1D = cy_i1D;
			thread_arg->cy_r1E = cy_r1E;		thread_arg->cy_i1E = cy_i1E;
			thread_arg->cy_r1F = cy_r1F;		thread_arg->cy_i1F = cy_i1F;
			thread_arg->cy_r20 = cy_r20;		thread_arg->cy_i20 = cy_i20;
			thread_arg->cy_r21 = cy_r21;		thread_arg->cy_i21 = cy_i21;
			thread_arg->cy_r22 = cy_r22;		thread_arg->cy_i22 = cy_i22;
			thread_arg->cy_r23 = cy_r23;		thread_arg->cy_i23 = cy_i23;
			thread_arg->cy_r24 = cy_r24;		thread_arg->cy_i24 = cy_i24;
			thread_arg->cy_r25 = cy_r25;		thread_arg->cy_i25 = cy_i25;
			thread_arg->cy_r26 = cy_r26;		thread_arg->cy_i26 = cy_i26;
			thread_arg->cy_r27 = cy_r27;		thread_arg->cy_i27 = cy_i27;
			thread_arg->cy_r28 = cy_r28;		thread_arg->cy_i28 = cy_i28;
			thread_arg->cy_r29 = cy_r29;		thread_arg->cy_i29 = cy_i29;
			thread_arg->cy_r2A = cy_r2A;		thread_arg->cy_i2A = cy_i2A;
			thread_arg->cy_r2B = cy_r2B;		thread_arg->cy_i2B = cy_i2B;
			thread_arg->cy_r2C = cy_r2C;		thread_arg->cy_i2C = cy_i2C;
			thread_arg->cy_r2D = cy_r2D;		thread_arg->cy_i2D = cy_i2D;
			thread_arg->cy_r2E = cy_r2E;		thread_arg->cy_i2E = cy_i2E;
			thread_arg->cy_r2F = cy_r2F;		thread_arg->cy_i2F = cy_i2F;
			thread_arg->cy_r30 = cy_r30;		thread_arg->cy_i30 = cy_i30;
			thread_arg->cy_r31 = cy_r31;		thread_arg->cy_i31 = cy_i31;
			thread_arg->cy_r32 = cy_r32;		thread_arg->cy_i32 = cy_i32;
			thread_arg->cy_r33 = cy_r33;		thread_arg->cy_i33 = cy_i33;
			thread_arg->cy_r34 = cy_r34;		thread_arg->cy_i34 = cy_i34;
			thread_arg->cy_r35 = cy_r35;		thread_arg->cy_i35 = cy_i35;
			thread_arg->cy_r36 = cy_r36;		thread_arg->cy_i36 = cy_i36;
			thread_arg->cy_r37 = cy_r37;		thread_arg->cy_i37 = cy_i37;
			thread_arg->cy_r38 = cy_r38;		thread_arg->cy_i38 = cy_i38;
			thread_arg->cy_r39 = cy_r39;		thread_arg->cy_i39 = cy_i39;
			thread_arg->cy_r3A = cy_r3A;		thread_arg->cy_i3A = cy_i3A;
			thread_arg->cy_r3B = cy_r3B;		thread_arg->cy_i3B = cy_i3B;
			thread_arg->cy_r3C = cy_r3C;		thread_arg->cy_i3C = cy_i3C;
			thread_arg->cy_r3D = cy_r3D;		thread_arg->cy_i3D = cy_i3D;
			thread_arg->cy_r3E = cy_r3E;		thread_arg->cy_i3E = cy_i3E;
			thread_arg->cy_r3F = cy_r3F;		thread_arg->cy_i3F = cy_i3F;

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

