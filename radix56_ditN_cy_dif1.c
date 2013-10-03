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

#define ERR_CHECK_ALL	/* #define this to do ROE checking of all convolution outputs, rather than just every 36th one */
#if defined(USE_AVX) && !defined(ERR_CHECK_ALL)
	#error ERR_CHECK_ALL *required* for AVX-mode builds!
#endif

/* Use for toggling higher-accuracy version of the twiddles computation */
//#define HIACC 0	<*** prefer to set via compile-time flag; default is FALSE [= LOACC]

#define EPS 1e-10

// See the radix28 version of this routine for details about the
// small-cyclic-array indexing scheme used in the fermat_carry_norm_errcheckB macros.

#ifdef USE_SSE2

	#include "sse2_macro.h"

  // For Mersenne-mod we need (16 [SSE2] or 64 [AVX]) + 4 added slots for the half_arr lookup tables.
  // For Fermat-mod in AVX mode we need RADIX*4 = 224 [if HIACC] or 12 [if not] slots for the compact
  // negacyclic-roots chained-multiply scheme.
  // Add larger number in each case - i.e. max(68,224) = 224 if AVX+HIACC, max(68,12) = 68 if AVX+LOACC, 20 if SSE2
  // to (half_arr_offset56 + RADIX) to get required value of radix56_creals_in_local_store:
  #ifdef USE_AVX
	const int half_arr_offset56 = 268;	// = 0x10c ... + RADIX = 324; Used for thread local-storage-integrity checking
   #if HIACC
	const int radix56_creals_in_local_store = 548;	// AVX+HIACC: 324 + 224 and round up to nearest multiple of 8
   #else
	const int radix56_creals_in_local_store = 396;	// AVX+LOACC: 324 +  68 and round up to nearest multiple of 8
   #endif
  #else
	const int half_arr_offset56 = 296;	// = 0x128 ... + RADIX = 352; Used for thread local-storage-integrity checking
	const int radix56_creals_in_local_store = 372;	// SSE2: 352 + 20 and round up to nearest multiple of 8
  #endif

#endif

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
		int *si;
		struct complex *rn0;
		struct complex *rn1;
		vec_dbl *r00r,*half_arr;

		int bjmodn00;	int bjmodn28;
		int bjmodn01;	int bjmodn29;
		int bjmodn02;	int bjmodn30;
		int bjmodn03;	int bjmodn31;
		int bjmodn04;	int bjmodn32;
		int bjmodn05;	int bjmodn33;
		int bjmodn06;	int bjmodn34;
		int bjmodn07;	int bjmodn35;
		int bjmodn08;	int bjmodn36;
		int bjmodn09;	int bjmodn37;
		int bjmodn10;	int bjmodn38;
		int bjmodn11;	int bjmodn39;
		int bjmodn12;	int bjmodn40;
		int bjmodn13;	int bjmodn41;
		int bjmodn14;	int bjmodn42;
		int bjmodn15;	int bjmodn43;
		int bjmodn16;	int bjmodn44;
		int bjmodn17;	int bjmodn45;
		int bjmodn18;	int bjmodn46;
		int bjmodn19;	int bjmodn47;
		int bjmodn20;	int bjmodn48;
		int bjmodn21;	int bjmodn49;
		int bjmodn22;	int bjmodn50;
		int bjmodn23;	int bjmodn51;
		int bjmodn24;	int bjmodn52;
		int bjmodn25;	int bjmodn53;
		int bjmodn26;	int bjmodn54;
		int bjmodn27;	int bjmodn55;

		// Pad to make size a multiple of 64 bytes:
		double dpad0;
		double dpad1;
		double dpad2;

		/* carries: */
		double cy_r00;	double cy_r28;
		double cy_r01;	double cy_r29;
		double cy_r02;	double cy_r30;
		double cy_r03;	double cy_r31;
		double cy_r04;	double cy_r32;
		double cy_r05;	double cy_r33;
		double cy_r06;	double cy_r34;
		double cy_r07;	double cy_r35;
		double cy_r08;	double cy_r36;
		double cy_r09;	double cy_r37;
		double cy_r10;	double cy_r38;
		double cy_r11;	double cy_r39;
		double cy_r12;	double cy_r40;
		double cy_r13;	double cy_r41;
		double cy_r14;	double cy_r42;
		double cy_r15;	double cy_r43;
		double cy_r16;	double cy_r44;
		double cy_r17;	double cy_r45;
		double cy_r18;	double cy_r46;
		double cy_r19;	double cy_r47;
		double cy_r20;	double cy_r48;
		double cy_r21;	double cy_r49;
		double cy_r22;	double cy_r50;
		double cy_r23;	double cy_r51;
		double cy_r24;	double cy_r52;
		double cy_r25;	double cy_r53;
		double cy_r26;	double cy_r54;
		double cy_r27;	double cy_r55;

		double cy_i00;	double cy_i28;
		double cy_i01;	double cy_i29;
		double cy_i02;	double cy_i30;
		double cy_i03;	double cy_i31;
		double cy_i04;	double cy_i32;
		double cy_i05;	double cy_i33;
		double cy_i06;	double cy_i34;
		double cy_i07;	double cy_i35;
		double cy_i08;	double cy_i36;
		double cy_i09;	double cy_i37;
		double cy_i10;	double cy_i38;
		double cy_i11;	double cy_i39;
		double cy_i12;	double cy_i40;
		double cy_i13;	double cy_i41;
		double cy_i14;	double cy_i42;
		double cy_i15;	double cy_i43;
		double cy_i16;	double cy_i44;
		double cy_i17;	double cy_i45;
		double cy_i18;	double cy_i46;
		double cy_i19;	double cy_i47;
		double cy_i20;	double cy_i48;
		double cy_i21;	double cy_i49;
		double cy_i22;	double cy_i50;
		double cy_i23;	double cy_i51;
		double cy_i24;	double cy_i52;
		double cy_i25;	double cy_i53;
		double cy_i26;	double cy_i54;
		double cy_i27;	double cy_i55;
	};

#endif

/**************/

int radix56_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-56 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-56 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix7/8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix56_ditN_cy_dif1";
	const int RADIX = 56, odd_radix = 7;	// odd_radix = [radix >> trailz(radix)]
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
	int NDIVR,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,ntmp,outer,nbytes;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30;
	static double radix_inv, n2inv;
#ifdef USE_SSE2
	uint32 nwt16 = nwt << l2_sz_vd;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code
#endif
	// Need these both in scalar mode and to ease the SSE2-array init...dimension = odd_radix;
	// In order to ease the ptr-access for the || routine, lump these 4*odd_radix doubles together with copies of
	// the 4 in the passed-in bs[2] and bsinv[2] arrays [and used in this 4-double form by the mersenne-mod carry macros]
	// into a single foo_array[4*(odd_radix+1)], then convert what used to be disparate odd_radix-sized arrays to pointers.
	static double foo_array[(odd_radix+1)<<2], *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr, bs,bsinv;
	// NOTE: If compiler squawks, Use literal value of odd_radix in dimensioning in order to allow us to explicitly make static:

#if defined(USE_SSE2) || !defined (LO_ADD)
	/* SSE2 version assumes LO_ADD = 0, i.e. the low-mul Nussbaumer-style DFT implementation: */
	const double	cx0 =-0.16666666666666666667,	/* (cc1+cc2+cc3)/3 */
				 	cx1 = 1.52445866976115265675, 	/*  cc1-cc3		*/
				 	cx2 = 0.67844793394610472196, 	/*  cc2-cc3		*/
				 	cx3 = 0.73430220123575245957,	/* (cc1+cc2-2*cc3)/3	*/
				/* Switch the sign of ss3 in these: */
				 	sx0 = 0.44095855184409843174,	/* (ss1+ss2-ss3)/3	*/
				 	sx1 = 1.21571522158558792920, 	/*  ss1+ss3		*/
				 	sx2 = 1.40881165129938172752, 	/*  ss2+ss3		*/
				 	sx3 = 0.87484229096165655224;	/* (ss1+ss2+2*ss3)/3	*/
#elif !defined(MULTITHREAD)
	/* Non-SSE2 version assumes LO_ADD = 1 and uses the corresponding versions of the sincos constants: */
	const double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
#endif

	double scale,
		t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F,
		t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F,
		t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F,
		t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F,
		t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t4A,t4B,t4C,t4D,t4E,t4F,
		t50,t51,t52,t53,t54,t55,t56,t57,t58,t59,t5A,t5B,t5C,t5D,t5E,t5F,
		t60,t61,t62,t63,t64,t65,t66,t67,t68,t69,t6A,t6B,t6C,t6D,t6E,t6F;
	double dtmp,maxerr = 0.0;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	static int n_div_nwt;
	int col,co2,co3,m,m2;
  #ifdef USE_AVX
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
  #endif
	double rt,it,wt_re,wt_im;	/* Fermat-mod weights stuff, used in both scalar and AVX mode */
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
	double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7;	/* Addresses into array sections */
  #endif

	const double crnd = 3.0*0x4000000*0x2000000;
	static vec_dbl *isrt2, *cc0, *ss0, *cc1, *ss1, *cc2, *ss2, *cc3, *ss3, *max_err, *sse2_rnd, *half_arr
	,*r00r,*r01r,*r02r,*r03r,*r04r,*r05r,*r06r,*r07r,*r08r,*r09r,*r10r,*r11r,*r12r,*r13r
	,*r14r,*r15r,*r16r,*r17r,*r18r,*r19r,*r20r,*r21r,*r22r,*r23r,*r24r,*r25r,*r26r,*r27r
	,*r28r,*r29r,*r30r,*r31r,*r32r,*r33r,*r34r,*r35r,*r36r,*r37r,*r38r,*r39r,*r40r,*r41r
	,*r42r,*r43r,*r44r,*r45r,*r46r,*r47r,*r48r,*r49r,*r50r,*r51r,*r52r,*r53r,*r54r,*r55r
	,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r
	,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r
	,*s1p28r,*s1p29r,*s1p30r,*s1p31r,*s1p32r,*s1p33r,*s1p34r,*s1p35r,*s1p36r,*s1p37r,*s1p38r,*s1p39r,*s1p40r,*s1p41r
	,*s1p42r,*s1p43r,*s1p44r,*s1p45r,*s1p46r,*s1p47r,*s1p48r,*s1p49r,*s1p50r,*s1p51r,*s1p52r,*s1p53r,*s1p54r,*s1p55r
  #ifndef USE_AVX
	;
  #else
	, *base_negacyclic_root;
  #endif

	vec_dbl *tmp,*tm2;	// Non-static utility ptrs
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	static int
	 *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13
	,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27
	,*bjmodn28,*bjmodn29,*bjmodn30,*bjmodn31,*bjmodn32,*bjmodn33,*bjmodn34,*bjmodn35,*bjmodn36,*bjmodn37,*bjmodn38,*bjmodn39,*bjmodn40,*bjmodn41
	,*bjmodn42,*bjmodn43,*bjmodn44,*bjmodn45,*bjmodn46,*bjmodn47,*bjmodn48,*bjmodn49,*bjmodn50,*bjmodn51,*bjmodn52,*bjmodn53,*bjmodn54,*bjmodn55;
	static vec_dbl
		*cy_r00,*cy_r04,*cy_r08,*cy_r12,*cy_r16,*cy_r20,*cy_r24,*cy_r28,*cy_r32,*cy_r36,*cy_r40,*cy_r44,*cy_r48,*cy_r52,
		*cy_i00,*cy_i04,*cy_i08,*cy_i12,*cy_i16,*cy_i20,*cy_i24,*cy_i28,*cy_i32,*cy_i36,*cy_i40,*cy_i44,*cy_i48,*cy_i52;
  #ifndef USE_AVX
	static vec_dbl
		*cy_r02,*cy_r06,*cy_r10,*cy_r14,*cy_r18,*cy_r22,*cy_r26,*cy_r30,*cy_r34,*cy_r38,*cy_r42,*cy_r46,*cy_r50,*cy_r54,
		*cy_i02,*cy_i06,*cy_i10,*cy_i14,*cy_i18,*cy_i22,*cy_i26,*cy_i30,*cy_i34,*cy_i38,*cy_i42,*cy_i46,*cy_i50,*cy_i54;
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
	static task_control_t   task_control = {NULL, (void*)cy56_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,
		bjmodn28,bjmodn29,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,bjmodn36,bjmodn37,bjmodn38,bjmodn39,bjmodn40,bjmodn41,bjmodn42,bjmodn43,bjmodn44,bjmodn45,bjmodn46,bjmodn47,bjmodn48,bjmodn49,bjmodn50,bjmodn51,bjmodn52,bjmodn53,bjmodn54,bjmodn55;
	double re,im,temp,frac
		,r00r,r01r,r02r,r03r,r04r,r05r,r06r,r07r,r08r,r09r,r10r,r11r,r12r,r13r,r14r,r15r,r16r,r17r,r18r,r19r,r20r,r21r,r22r,r23r,r24r,r25r,r26r,r27r
		,r28r,r29r,r30r,r31r,r32r,r33r,r34r,r35r,r36r,r37r,r38r,r39r,r40r,r41r,r42r,r43r,r44r,r45r,r46r,r47r,r48r,r49r,r50r,r51r,r52r,r53r,r54r,r55r
		,r00i,r01i,r02i,r03i,r04i,r05i,r06i,r07i,r08i,r09i,r10i,r11i,r12i,r13i,r14i,r15i,r16i,r17i,r18i,r19i,r20i,r21i,r22i,r23i,r24i,r25i,r26i,r27i
		,r28i,r29i,r30i,r31i,r32i,r33i,r34i,r35i,r36i,r37i,r38i,r39i,r40i,r41i,r42i,r43i,r44i,r45i,r46i,r47i,r48i,r49i,r50i,r51i,r52i,r53i,r54i,r55i
		,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r
		,a1p28r,a1p29r,a1p30r,a1p31r,a1p32r,a1p33r,a1p34r,a1p35r,a1p36r,a1p37r,a1p38r,a1p39r,a1p40r,a1p41r,a1p42r,a1p43r,a1p44r,a1p45r,a1p46r,a1p47r,a1p48r,a1p49r,a1p50r,a1p51r,a1p52r,a1p53r,a1p54r,a1p55r
		,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i
		,a1p28i,a1p29i,a1p30i,a1p31i,a1p32i,a1p33i,a1p34i,a1p35i,a1p36i,a1p37i,a1p38i,a1p39i,a1p40i,a1p41i,a1p42i,a1p43i,a1p44i,a1p45i,a1p46i,a1p47i,a1p48i,a1p49i,a1p50i,a1p51i,a1p52i,a1p53i,a1p54i,a1p55i
		,cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27
		,cy_r28,cy_r29,cy_r30,cy_r31,cy_r32,cy_r33,cy_r34,cy_r35,cy_r36,cy_r37,cy_r38,cy_r39,cy_r40,cy_r41,cy_r42,cy_r43,cy_r44,cy_r45,cy_r46,cy_r47,cy_r48,cy_r49,cy_r50,cy_r51,cy_r52,cy_r53,cy_r54,cy_r55
		,cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27
		,cy_i28,cy_i29,cy_i30,cy_i31,cy_i32,cy_i33,cy_i34,cy_i35,cy_i36,cy_i37,cy_i38,cy_i39,cy_i40,cy_i41,cy_i42,cy_i43,cy_i44,cy_i45,cy_i46,cy_i47,cy_i48,cy_i49,cy_i50,cy_i51,cy_i52,cy_i53,cy_i54,cy_i55;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_bjmodn00 = 0x0,*_bjmodn01 = 0x0,*_bjmodn02 = 0x0,*_bjmodn03 = 0x0,*_bjmodn04 = 0x0,*_bjmodn05 = 0x0,*_bjmodn06 = 0x0,*_bjmodn07 = 0x0,*_bjmodn08 = 0x0,*_bjmodn09 = 0x0,*_bjmodn10 = 0x0,*_bjmodn11 = 0x0,*_bjmodn12 = 0x0,*_bjmodn13 = 0x0,*_bjmodn14 = 0x0,*_bjmodn15 = 0x0,*_bjmodn16 = 0x0,*_bjmodn17 = 0x0,*_bjmodn18 = 0x0,*_bjmodn19 = 0x0,*_bjmodn20 = 0x0,*_bjmodn21 = 0x0,*_bjmodn22 = 0x0,*_bjmodn23 = 0x0,*_bjmodn24 = 0x0,*_bjmodn25 = 0x0,*_bjmodn26 = 0x0,*_bjmodn27 = 0x0
			  ,*_bjmodn28 = 0x0,*_bjmodn29 = 0x0,*_bjmodn30 = 0x0,*_bjmodn31 = 0x0,*_bjmodn32 = 0x0,*_bjmodn33 = 0x0,*_bjmodn34 = 0x0,*_bjmodn35 = 0x0,*_bjmodn36 = 0x0,*_bjmodn37 = 0x0,*_bjmodn38 = 0x0,*_bjmodn39 = 0x0,*_bjmodn40 = 0x0,*_bjmodn41 = 0x0,*_bjmodn42 = 0x0,*_bjmodn43 = 0x0,*_bjmodn44 = 0x0,*_bjmodn45 = 0x0,*_bjmodn46 = 0x0,*_bjmodn47 = 0x0,*_bjmodn48 = 0x0,*_bjmodn49 = 0x0,*_bjmodn50 = 0x0,*_bjmodn51 = 0x0,*_bjmodn52 = 0x0,*_bjmodn53 = 0x0,*_bjmodn54 = 0x0,*_bjmodn55 = 0x0;
	static int *_bjmodnini = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_maxerr = 0x0
	,*_cy_r00 = 0x0,*_cy_r01 = 0x0,*_cy_r02 = 0x0,*_cy_r03 = 0x0,*_cy_r04 = 0x0,*_cy_r05 = 0x0,*_cy_r06 = 0x0,*_cy_r07 = 0x0,*_cy_r08 = 0x0,*_cy_r09 = 0x0,*_cy_r10 = 0x0,*_cy_r11 = 0x0,*_cy_r12 = 0x0,*_cy_r13 = 0x0,*_cy_r14 = 0x0,*_cy_r15 = 0x0,*_cy_r16 = 0x0,*_cy_r17 = 0x0,*_cy_r18 = 0x0,*_cy_r19 = 0x0,*_cy_r20 = 0x0,*_cy_r21 = 0x0,*_cy_r22 = 0x0,*_cy_r23 = 0x0,*_cy_r24 = 0x0,*_cy_r25 = 0x0,*_cy_r26 = 0x0,*_cy_r27 = 0x0
	,*_cy_r28 = 0x0,*_cy_r29 = 0x0,*_cy_r30 = 0x0,*_cy_r31 = 0x0,*_cy_r32 = 0x0,*_cy_r33 = 0x0,*_cy_r34 = 0x0,*_cy_r35 = 0x0,*_cy_r36 = 0x0,*_cy_r37 = 0x0,*_cy_r38 = 0x0,*_cy_r39 = 0x0,*_cy_r40 = 0x0,*_cy_r41 = 0x0,*_cy_r42 = 0x0,*_cy_r43 = 0x0,*_cy_r44 = 0x0,*_cy_r45 = 0x0,*_cy_r46 = 0x0,*_cy_r47 = 0x0,*_cy_r48 = 0x0,*_cy_r49 = 0x0,*_cy_r50 = 0x0,*_cy_r51 = 0x0,*_cy_r52 = 0x0,*_cy_r53 = 0x0,*_cy_r54 = 0x0,*_cy_r55 = 0x0
	,*_cy_i00 = 0x0,*_cy_i01 = 0x0,*_cy_i02 = 0x0,*_cy_i03 = 0x0,*_cy_i04 = 0x0,*_cy_i05 = 0x0,*_cy_i06 = 0x0,*_cy_i07 = 0x0,*_cy_i08 = 0x0,*_cy_i09 = 0x0,*_cy_i10 = 0x0,*_cy_i11 = 0x0,*_cy_i12 = 0x0,*_cy_i13 = 0x0,*_cy_i14 = 0x0,*_cy_i15 = 0x0,*_cy_i16 = 0x0,*_cy_i17 = 0x0,*_cy_i18 = 0x0,*_cy_i19 = 0x0,*_cy_i20 = 0x0,*_cy_i21 = 0x0,*_cy_i22 = 0x0,*_cy_i23 = 0x0,*_cy_i24 = 0x0,*_cy_i25 = 0x0,*_cy_i26 = 0x0,*_cy_i27 = 0x0
	,*_cy_i28 = 0x0,*_cy_i29 = 0x0,*_cy_i30 = 0x0,*_cy_i31 = 0x0,*_cy_i32 = 0x0,*_cy_i33 = 0x0,*_cy_i34 = 0x0,*_cy_i35 = 0x0,*_cy_i36 = 0x0,*_cy_i37 = 0x0,*_cy_i38 = 0x0,*_cy_i39 = 0x0,*_cy_i40 = 0x0,*_cy_i41 = 0x0,*_cy_i42 = 0x0,*_cy_i43 = 0x0,*_cy_i44 = 0x0,*_cy_i45 = 0x0,*_cy_i46 = 0x0,*_cy_i47 = 0x0,*_cy_i48 = 0x0,*_cy_i49 = 0x0,*_cy_i50 = 0x0,*_cy_i51 = 0x0,*_cy_i52 = 0x0,*_cy_i53 = 0x0,*_cy_i54 = 0x0,*_cy_i55 = 0x0;

	foo_array[0] = base[0];
	foo_array[1] = base[1];
	foo_array[2] = baseinv[0];
	foo_array[3] = baseinv[1];
	wt_arr    = foo_array + 4;
	wtinv_arr = wt_arr    + odd_radix;
	bs_arr    = wtinv_arr + odd_radix;
	bsinv_arr = bs_arr    + odd_radix;

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=ii0=ii1=ii2=ii3=ii4=ii5=ii6=-1;
/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/RADIX;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/%d in %s.\n", iter,RADIX,func);
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
		ASSERT(HERE, LO_ADD,"LO_ADD");
		psave = p;
		first_entry=FALSE;
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
		ASSERT(HERE, CY_THREADS >= NTHREADS,"CY_THREADS < NTHREADS");
		ASSERT(HERE, isPow2(CY_THREADS)    ,"CY_THREADS not a power of 2!");
		if(CY_THREADS > 1)
		{
			ASSERT(HERE, NDIVR    %CY_THREADS == 0,"NDIVR    %CY_THREADS != 0");
			ASSERT(HERE, n_div_nwt%CY_THREADS == 0,"n_div_nwt%CY_THREADS != 0");
		}

	  #ifdef USE_PTHREAD

		j = (uint32)sizeof(struct cy_thread_data_t);
		tdat = (struct cy_thread_data_t *)calloc(CY_THREADS, j);

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

			pool_work_units = CY_THREADS;
			ASSERT(HERE, 0x0 != (tpool = threadpool_init(CY_THREADS, MAX_THREADS, CY_THREADS, &thread_control)), "threadpool_init failed!");

		#endif

		fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);

	  #endif

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

	#ifdef USE_SSE2

		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use vector-double type size (16 bytes for SSE2, 32 for AVX) to alloc a block of local storage
		// consisting of radix56_creals_in_local_store vec_dbl and (8+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix56_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix56_creals_in_local_store);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low slots of sc_arr for temporaries, next few for the nontrivial complex 16th roots,
	next few for the doubled carry pairs, next 2 for ROE and RND_CONST, next RADIX for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array. We also use a few more slots in AVX mode
	for the compact negacyclic-roots chained-multiply scheme, which drastically reduces accesses to the rn0,rn1 tables.
	*/
	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif
								tmp = sc_ptr + 0x38;
		r00r = sc_ptr + 0x00;	r28r = tmp + 0x00;
		r01r = sc_ptr + 0x02;	r29r = tmp + 0x02;
		r02r = sc_ptr + 0x04;	r30r = tmp + 0x04;
		r03r = sc_ptr + 0x06;	r31r = tmp + 0x06;
		r04r = sc_ptr + 0x08;	r32r = tmp + 0x08;
		r05r = sc_ptr + 0x0a;	r33r = tmp + 0x0a;
		r06r = sc_ptr + 0x0c;	r34r = tmp + 0x0c;
		r07r = sc_ptr + 0x0e;	r35r = tmp + 0x0e;
		r08r = sc_ptr + 0x10;	r36r = tmp + 0x10;
		r09r = sc_ptr + 0x12;	r37r = tmp + 0x12;
		r10r = sc_ptr + 0x14;	r38r = tmp + 0x14;
		r11r = sc_ptr + 0x16;	r39r = tmp + 0x16;
		r12r = sc_ptr + 0x18;	r40r = tmp + 0x18;
		r13r = sc_ptr + 0x1a;	r41r = tmp + 0x1a;
		r14r = sc_ptr + 0x1c;	r42r = tmp + 0x1c;
		r15r = sc_ptr + 0x1e;	r43r = tmp + 0x1e;
		r16r = sc_ptr + 0x20;	r44r = tmp + 0x20;
		r17r = sc_ptr + 0x22;	r45r = tmp + 0x22;
		r18r = sc_ptr + 0x24;	r46r = tmp + 0x24;
		r19r = sc_ptr + 0x26;	r47r = tmp + 0x26;
		r20r = sc_ptr + 0x28;	r48r = tmp + 0x28;
		r21r = sc_ptr + 0x2a;	r49r = tmp + 0x2a;
		r22r = sc_ptr + 0x2c;	r50r = tmp + 0x2c;
		r23r = sc_ptr + 0x2e;	r51r = tmp + 0x2e;
		r24r = sc_ptr + 0x30;	r52r = tmp + 0x30;
		r25r = sc_ptr + 0x32;	r53r = tmp + 0x32;
		r26r = sc_ptr + 0x34;	r54r = tmp + 0x34;
		r27r = sc_ptr + 0x36;	r55r = tmp + 0x36;
		tmp += 0x38;	// sc_ptr += 0x70
		tm2 = tmp + 0x38;
		s1p00r = tmp + 0x00;	s1p28r = tm2 + 0x00;
		s1p01r = tmp + 0x02;	s1p29r = tm2 + 0x02;
		s1p02r = tmp + 0x04;	s1p30r = tm2 + 0x04;
		s1p03r = tmp + 0x06;	s1p31r = tm2 + 0x06;
		s1p04r = tmp + 0x08;	s1p32r = tm2 + 0x08;
		s1p05r = tmp + 0x0a;	s1p33r = tm2 + 0x0a;
		s1p06r = tmp + 0x0c;	s1p34r = tm2 + 0x0c;
		s1p07r = tmp + 0x0e;	s1p35r = tm2 + 0x0e;
		s1p08r = tmp + 0x10;	s1p36r = tm2 + 0x10;
		s1p09r = tmp + 0x12;	s1p37r = tm2 + 0x12;
		s1p10r = tmp + 0x14;	s1p38r = tm2 + 0x14;
		s1p11r = tmp + 0x16;	s1p39r = tm2 + 0x16;
		s1p12r = tmp + 0x18;	s1p40r = tm2 + 0x18;
		s1p13r = tmp + 0x1a;	s1p41r = tm2 + 0x1a;
		s1p14r = tmp + 0x1c;	s1p42r = tm2 + 0x1c;
		s1p15r = tmp + 0x1e;	s1p43r = tm2 + 0x1e;
		s1p16r = tmp + 0x20;	s1p44r = tm2 + 0x20;
		s1p17r = tmp + 0x22;	s1p45r = tm2 + 0x22;
		s1p18r = tmp + 0x24;	s1p46r = tm2 + 0x24;
		s1p19r = tmp + 0x26;	s1p47r = tm2 + 0x26;
		s1p20r = tmp + 0x28;	s1p48r = tm2 + 0x28;
		s1p21r = tmp + 0x2a;	s1p49r = tm2 + 0x2a;
		s1p22r = tmp + 0x2c;	s1p50r = tm2 + 0x2c;
		s1p23r = tmp + 0x2e;	s1p51r = tm2 + 0x2e;
		s1p24r = tmp + 0x30;	s1p52r = tm2 + 0x30;
		s1p25r = tmp + 0x32;	s1p53r = tm2 + 0x32;
		s1p26r = tmp + 0x34;	s1p54r = tm2 + 0x34;
		s1p27r = tmp + 0x36;	s1p55r = tm2 + 0x36;
		tmp += 0x71;	// sc_ptr += 0xe1 [2x +0x38, plut jump of 1 for isrt2]
		isrt2   = tmp - 0x01;
		cc0		= tmp + 0x00;
		ss0		= tmp + 0x01;
		cc1		= tmp + 0x02;
		ss1		= tmp + 0x03;
		cc2		= tmp + 0x04;
		ss2		= tmp + 0x05;
		cc3  	= tmp + 0x06;
		ss3		= tmp + 0x07;	/* Pad with extra 4 slots for scratch storage needed by SSE2_RADIX_07_DFT macro here, plus 1 to make offset even */
		tmp += 0x08 + 0x5;	// sc_ptr += 0xee
	  #ifdef USE_AVX
		cy_r00 = tmp + 0x00;
		cy_r04 = tmp + 0x01;
		cy_r08 = tmp + 0x02;
		cy_r12 = tmp + 0x03;
		cy_r16 = tmp + 0x04;
		cy_r20 = tmp + 0x05;
		cy_r24 = tmp + 0x06;
		cy_r28 = tmp + 0x07;
		cy_r32 = tmp + 0x08;
		cy_r36 = tmp + 0x09;
		cy_r40 = tmp + 0x0a;
		cy_r44 = tmp + 0x0b;
		cy_r48 = tmp + 0x0c;
		cy_r52 = tmp + 0x0d;
		cy_i00 = tmp + 0x0e;
		cy_i04 = tmp + 0x0f;
		cy_i08 = tmp + 0x10;
		cy_i12 = tmp + 0x11;
		cy_i16 = tmp + 0x12;
		cy_i20 = tmp + 0x13;
		cy_i24 = tmp + 0x14;
		cy_i28 = tmp + 0x15;
		cy_i32 = tmp + 0x16;
		cy_i36 = tmp + 0x17;
		cy_i40 = tmp + 0x18;
		cy_i44 = tmp + 0x19;
		cy_i48 = tmp + 0x1a;
		cy_i52 = tmp + 0x1b;
		tmp += 0x1c;	// sc_ptr += 0x10a
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x10c; This is where the value of half_arr_offset56 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
	  #else
		cy_r00 = tmp + 0x00;	cy_r02 = tmp + 0x01;
		cy_r04 = tmp + 0x02;	cy_r06 = tmp + 0x03;
		cy_r08 = tmp + 0x04;	cy_r10 = tmp + 0x05;
		cy_r12 = tmp + 0x06;	cy_r14 = tmp + 0x07;
		cy_r16 = tmp + 0x08;	cy_r18 = tmp + 0x09;
		cy_r20 = tmp + 0x0a;	cy_r22 = tmp + 0x0b;
		cy_r24 = tmp + 0x0c;	cy_r26 = tmp + 0x0d;
		cy_r28 = tmp + 0x0e;	cy_r30 = tmp + 0x0f;
		cy_r32 = tmp + 0x10;	cy_r34 = tmp + 0x11;
		cy_r36 = tmp + 0x12;	cy_r38 = tmp + 0x13;
		cy_r40 = tmp + 0x14;	cy_r42 = tmp + 0x15;
		cy_r44 = tmp + 0x16;	cy_r46 = tmp + 0x17;
		cy_r48 = tmp + 0x18;	cy_r50 = tmp + 0x19;
		cy_r52 = tmp + 0x1a;	cy_r54 = tmp + 0x1b;
		tmp += 0x1c;
		cy_i00 = tmp + 0x00;	cy_i02 = tmp + 0x01;
		cy_i04 = tmp + 0x02;	cy_i06 = tmp + 0x03;
		cy_i08 = tmp + 0x04;	cy_i10 = tmp + 0x05;
		cy_i12 = tmp + 0x06;	cy_i14 = tmp + 0x07;
		cy_i16 = tmp + 0x08;	cy_i18 = tmp + 0x09;
		cy_i20 = tmp + 0x0a;	cy_i22 = tmp + 0x0b;
		cy_i24 = tmp + 0x0c;	cy_i26 = tmp + 0x0d;
		cy_i28 = tmp + 0x0e;	cy_i30 = tmp + 0x0f;
		cy_i32 = tmp + 0x10;	cy_i34 = tmp + 0x11;
		cy_i36 = tmp + 0x12;	cy_i38 = tmp + 0x13;
		cy_i40 = tmp + 0x14;	cy_i42 = tmp + 0x15;
		cy_i44 = tmp + 0x16;	cy_i46 = tmp + 0x17;
		cy_i48 = tmp + 0x18;	cy_i50 = tmp + 0x19;
		cy_i52 = tmp + 0x1a;	cy_i54 = tmp + 0x1b;
		tmp += 0x1c;	// sc_ptr += 0x126
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x128; This is where the value of half_arr_offset56 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*odd_radix] x 16 for Fermat-mod */
	  #endif

		/* These remain fixed: */
		VEC_DBL_INIT( isrt2, ISRT2);
		/* cc0 = (cc1+cc2+cc3)/3 - 1; subtract 1 from Nussbaumer's definition in order to ease in-place computation */
		VEC_DBL_INIT(cc0, cx0-1);	VEC_DBL_INIT(ss0, sx0);
		VEC_DBL_INIT(cc1, cx1  );	VEC_DBL_INIT(ss1, sx1);
		VEC_DBL_INIT(cc2, cx2  );	VEC_DBL_INIT(ss2, sx2);
		VEC_DBL_INIT(cc3, cx3  );	VEC_DBL_INIT(ss3, sx3);
		VEC_DBL_INIT(sse2_rnd, crnd);		/* SSE2 math = 53-mantissa-bit IEEE double-float: */

		// Propagate the above consts to the remaining threads:
		nbytes = (int)ss3 - (int)isrt2 + sz_vd;	// #bytes in 1st of above block of consts
		tmp = isrt2;
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
				printf("j = %3u: cos = 0x%16llX\n",j,qfdbl_as_uint64(qx));
			}
			exit(0);
		#endif

		tmp = base_negacyclic_root + RADIX*2;	// First 112 = 2*RADIX slots reserved for RADIX/4 copies of the Re/Im parts of the 4 base multipliers
		tm2 = tmp + RADIX/2 - 1;
										tmp->d0 = 1.0;	(tmp+1)->d0 = 0.0;
		tmp64 = 0x3FEFFCC70925C367ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(01*I*Pi/112) */
		tmp64 = 0x3FEFF31CCABE6E4Cull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(02*I*Pi/112) */
		tmp64 = 0x3FEFE303371EABCBull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(03*I*Pi/112) */	tmp += 2;
		tmp64 = 0x3FEFCC7D8C64135Full;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(04*I*Pi/112) */	tm2 -= 2;
		tmp64 = 0x3FEFAF9053CDF7E6ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(05*I*Pi/112) */
		tmp64 = 0x3FEF8C4160D38565ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(06*I*Pi/112) */
		tmp64 = 0x3FEF6297CFF75CB0ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(07*I*Pi/112) */	tmp += 2;
		tmp64 = 0x3FEF329C0558E969ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(08*I*Pi/112) */	tm2 -= 2;
		tmp64 = 0x3FEEFC57AB03BC37ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(09*I*Pi/112) */
		tmp64 = 0x3FEEBFD5AEFD405Aull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(10*I*Pi/112) */
		tmp64 = 0x3FEE7D22411130F5ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(11*I*Pi/112) */	tmp += 2;
		tmp64 = 0x3FEE344AD05D3F86ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(12*I*Pi/112) */	tm2 -= 2;
		tmp64 = 0x3FEDE55E089C6A31ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(13*I*Pi/112) */
		tmp64 = 0x3FED906BCF328D46ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(14*I*Pi/112) */
		tmp64 = 0x3FED35853FF8C869ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(15*I*Pi/112) */	tmp += 2;
		tmp64 = 0x3FECD4BCA9CB5C71ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(16*I*Pi/112) */	tm2 -= 2;
		tmp64 = 0x3FEC6E258AD9B3A1ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(17*I*Pi/112) */
		tmp64 = 0x3FEC01D48CB95263ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(18*I*Pi/112) */
		tmp64 = 0x3FEB8FDF803C7AE7ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(19*I*Pi/112) */	tmp += 2;
		tmp64 = 0x3FEB185D590D5A44ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(20*I*Pi/112) */	tm2 -= 2;
		tmp64 = 0x3FEA9B66290EA1A3ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(21*I*Pi/112) */
		tmp64 = 0x3FEA19131B8279C4ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(22*I*Pi/112) */
		tmp64 = 0x3FE9917E6FF8CAE3ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(23*I*Pi/112) */	tmp += 2;
		tmp64 = 0x3FE904C37505DE4Bull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(24*I*Pi/112) */	tm2 -= 2;
		tmp64 = 0x3FE872FE82C26A36ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(25*I*Pi/112) */
		tmp64 = 0x3FE7DC4CF5162385ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(26*I*Pi/112) */
		tmp64 = 0x3FE740CD25CDFBB2ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(27*I*Pi/112) */	tmp += 2;
		tmp64 = 0x3FE6A09E667F3BCDull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(28*I*Pi/112) */	tm2 -= 2;
		tmp64 = 0x3FE5FBE0FA38B7C6ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(29*I*Pi/112) */
		tmp64 = 0x3FE552B60F035F34ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(30*I*Pi/112) */
		tmp64 = 0x3FE4A53FB7337A9Bull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(31*I*Pi/112) */	tmp += 2;
		tmp64 = 0x3FE3F3A0E28BEDD1ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(32*I*Pi/112) */	tm2 -= 2;
		tmp64 = 0x3FE33DFD5734E146ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(33*I*Pi/112) */
		tmp64 = 0x3FE28479AA873CFEull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(34*I*Pi/112) */
		tmp64 = 0x3FE1C73B39AE68C8ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(35*I*Pi/112) */	tmp += 2;
		tmp64 = 0x3FE106682221CD8Aull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(36*I*Pi/112) */	tm2 -= 2;
		tmp64 = 0x3FE0422739F79BACull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(37*I*Pi/112) */
		tmp64 = 0x3FDEF5401024C4F4ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(38*I*Pi/112) */
		tmp64 = 0x3FDD5FF578561769ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(39*I*Pi/112) */	tmp += 2;
		tmp64 = 0x3FDBC4C04D71ABC1ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(40*I*Pi/112) */	tm2 -= 2;
		tmp64 = 0x3FDA23F36171DD2Dull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(41*I*Pi/112) */
		tmp64 = 0x3FD87DE2A6AEA963ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(42*I*Pi/112) */
		tmp64 = 0x3FD6D2E31EF560C0ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(43*I*Pi/112) */	tmp += 2;
		tmp64 = 0x3FD5234ACA69A9FEull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(44*I*Pi/112) */	tm2 -= 2;
		tmp64 = 0x3FD36F7096334C28ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(45*I*Pi/112) */
		tmp64 = 0x3FD1B7AC4AFC3C02ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(46*I*Pi/112) */
		tmp64 = 0x3FCFF8ACF684E6EDull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(47*I*Pi/112) */	tmp += 2;
		tmp64 = 0x3FCC7B90E3024582ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(48*I*Pi/112) */	tm2 -= 2;
		tmp64 = 0x3FC8F8B83C69A60Bull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(49*I*Pi/112) */
		tmp64 = 0x3FC570D80B7C3350ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(50*I*Pi/112) */
		tmp64 = 0x3FC1E4A65C4D04B5ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(51*I*Pi/112) */	tmp += 2;
		tmp64 = 0x3FBCA9B4332D6F61ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(52*I*Pi/112) */	tm2 -= 2;
		tmp64 = 0x3FB58455CFC8625Dull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(53*I*Pi/112) */
		tmp64 = 0x3FACB544024FC940ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(54*I*Pi/112) */
		tmp64 = 0x3F9CB82865EB9D38ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(55*I*Pi/112) */	tmp += 2;

		tmp = base_negacyclic_root + RADIX*2;	// reset to point to start of above block
		nbytes = RADIX*sz_vd/2;	// 7 AVX-register-sized complex data

	   #else	// HIACC = false:

		// lower-precision version, which yields slightly more roundoff error, but is simpler and more storage-compact.
		// Init exp(j*I*Pi/2/RADIX), for j = 0-3:
		tmp = base_negacyclic_root + 8;	// First 8 slots reserved for Re/Im parts of the 4 base multipliers
		tmp->d0 = 1.0;
		tmp64 = 0x3FEFFCC70925C367ull;	tmp->d1 = *(double *)&tmp64;	// cos(01*I*Pi/120)
		tmp64 = 0x3FEFF31CCABE6E4Cull;	tmp->d2 = *(double *)&tmp64;	// cos(02*I*Pi/120)
		tmp64 = 0x3FEFE303371EABCBull;	tmp->d3 = *(double *)&tmp64;	// cos(03*I*Pi/120)

		(++tmp)->d0 = 0.0;
		tmp64 = 0x3FB58455CFC8625Dull;	tmp->d3 = *(double *)&tmp64;	// sin(03*I*Pi/120) = cos(53*I*Pi/112)
		tmp64 = 0x3FACB544024FC940ull;	tmp->d2 = *(double *)&tmp64;	// sin(02*I*Pi/120) = cos(54*I*Pi/112)
		tmp64 = 0x3F9CB82865EB9D38ull;	tmp->d1 = *(double *)&tmp64;	// sin(01*I*Pi/120) = cos(55*I*Pi/112)
		++tmp;
		tmp64 = 0x3FEFCC7D8C64135Full;	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// cos(04*I*Pi/120)
		++tmp;
		tmp64 = 0x3FBCA9B4332D6F61ull;	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// sin(04*I*Pi/120) = cos(52*I*Pi/112)
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

		nbytes = 64 << l2_sz_vd;

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

		nbytes = 16 << l2_sz_vd;

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
									bjmodn28 = bjmodn00 + 28;
		bjmodn01 = bjmodn00 +  1;	bjmodn29 = bjmodn00 + 29;
		bjmodn02 = bjmodn00 +  2;	bjmodn30 = bjmodn00 + 30;
		bjmodn03 = bjmodn00 +  3;	bjmodn31 = bjmodn00 + 31;
		bjmodn04 = bjmodn00 +  4;	bjmodn32 = bjmodn00 + 32;
		bjmodn05 = bjmodn00 +  5;	bjmodn33 = bjmodn00 + 33;
		bjmodn06 = bjmodn00 +  6;	bjmodn34 = bjmodn00 + 34;
		bjmodn07 = bjmodn00 +  7;	bjmodn35 = bjmodn00 + 35;
		bjmodn08 = bjmodn00 +  8;	bjmodn36 = bjmodn00 + 36;
		bjmodn09 = bjmodn00 +  9;	bjmodn37 = bjmodn00 + 37;
		bjmodn10 = bjmodn00 + 10;	bjmodn38 = bjmodn00 + 38;
		bjmodn11 = bjmodn00 + 11;	bjmodn39 = bjmodn00 + 39;
		bjmodn12 = bjmodn00 + 12;	bjmodn40 = bjmodn00 + 40;
		bjmodn13 = bjmodn00 + 13;	bjmodn41 = bjmodn00 + 41;
		bjmodn14 = bjmodn00 + 14;	bjmodn42 = bjmodn00 + 42;
		bjmodn15 = bjmodn00 + 15;	bjmodn43 = bjmodn00 + 43;
		bjmodn16 = bjmodn00 + 16;	bjmodn44 = bjmodn00 + 44;
		bjmodn17 = bjmodn00 + 17;	bjmodn45 = bjmodn00 + 45;
		bjmodn18 = bjmodn00 + 18;	bjmodn46 = bjmodn00 + 46;
		bjmodn19 = bjmodn00 + 19;	bjmodn47 = bjmodn00 + 47;
		bjmodn20 = bjmodn00 + 20;	bjmodn48 = bjmodn00 + 48;
		bjmodn21 = bjmodn00 + 21;	bjmodn49 = bjmodn00 + 49;
		bjmodn22 = bjmodn00 + 22;	bjmodn50 = bjmodn00 + 50;
		bjmodn23 = bjmodn00 + 23;	bjmodn51 = bjmodn00 + 51;
		bjmodn24 = bjmodn00 + 24;	bjmodn52 = bjmodn00 + 52;
		bjmodn25 = bjmodn00 + 25;	bjmodn53 = bjmodn00 + 53;
		bjmodn26 = bjmodn00 + 26;	bjmodn54 = bjmodn00 + 54;
		bjmodn27 = bjmodn00 + 27;	bjmodn55 = bjmodn00 + 55;

	/*********** Defer the per-thread local-mem-block copy until after added wts-index precomputation below ************/
	#endif	/* USE_SSE2 */

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
		#ifdef USE_SSE2
			tdat[ithread].r00r = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (long)tdat[ithread].r00r + ((long)half_arr - (long)r00r);
		#else
			// In scalar mode use these 2 ptrs to pass wts_idx_incr and the base/baseinv/etc array-ptrs:
			tdat[ithread].r00r     = (vec_dbl *)foo_array;
			tdat[ithread].half_arr = (vec_dbl *)&wts_idx_incr;
		#endif	// USE_SSE2
		}
	#endif

		/*   constant index offsets for array load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		pini += ( (pini >> DAT_BITS) << PAD_BITS );
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p10 = p08 + p08;
		p18 = p10 + p08;
		p20 = p18 + p08;
		p28 = p20 + p08;
		p30 = p28 + p08;

		p01 += ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 += ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 += ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 += ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 += ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 += ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 += ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 += ( (p08 >> DAT_BITS) << PAD_BITS );
		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p18 += ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p28 += ( (p28 >> DAT_BITS) << PAD_BITS );
		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );

		if(_cy_r00)	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;

			free((void *)_bjmodn00); _bjmodn00 = 0x0;	free((void *)_bjmodn28); _bjmodn28 = 0x0;
			free((void *)_bjmodn01); _bjmodn01 = 0x0;	free((void *)_bjmodn29); _bjmodn29 = 0x0;
			free((void *)_bjmodn02); _bjmodn02 = 0x0;	free((void *)_bjmodn30); _bjmodn30 = 0x0;
			free((void *)_bjmodn03); _bjmodn03 = 0x0;	free((void *)_bjmodn31); _bjmodn31 = 0x0;
			free((void *)_bjmodn04); _bjmodn04 = 0x0;	free((void *)_bjmodn32); _bjmodn32 = 0x0;
			free((void *)_bjmodn05); _bjmodn05 = 0x0;	free((void *)_bjmodn33); _bjmodn33 = 0x0;
			free((void *)_bjmodn06); _bjmodn06 = 0x0;	free((void *)_bjmodn34); _bjmodn34 = 0x0;
			free((void *)_bjmodn07); _bjmodn07 = 0x0;	free((void *)_bjmodn35); _bjmodn35 = 0x0;
			free((void *)_bjmodn08); _bjmodn08 = 0x0;	free((void *)_bjmodn36); _bjmodn36 = 0x0;
			free((void *)_bjmodn09); _bjmodn09 = 0x0;	free((void *)_bjmodn37); _bjmodn37 = 0x0;
			free((void *)_bjmodn10); _bjmodn10 = 0x0;	free((void *)_bjmodn38); _bjmodn38 = 0x0;
			free((void *)_bjmodn11); _bjmodn11 = 0x0;	free((void *)_bjmodn39); _bjmodn39 = 0x0;
			free((void *)_bjmodn12); _bjmodn12 = 0x0;	free((void *)_bjmodn40); _bjmodn40 = 0x0;
			free((void *)_bjmodn13); _bjmodn13 = 0x0;	free((void *)_bjmodn41); _bjmodn41 = 0x0;
			free((void *)_bjmodn14); _bjmodn14 = 0x0;	free((void *)_bjmodn42); _bjmodn42 = 0x0;
			free((void *)_bjmodn15); _bjmodn15 = 0x0;	free((void *)_bjmodn43); _bjmodn43 = 0x0;
			free((void *)_bjmodn16); _bjmodn16 = 0x0;	free((void *)_bjmodn44); _bjmodn44 = 0x0;
			free((void *)_bjmodn17); _bjmodn17 = 0x0;	free((void *)_bjmodn45); _bjmodn45 = 0x0;
			free((void *)_bjmodn18); _bjmodn18 = 0x0;	free((void *)_bjmodn46); _bjmodn46 = 0x0;
			free((void *)_bjmodn19); _bjmodn19 = 0x0;	free((void *)_bjmodn47); _bjmodn47 = 0x0;
			free((void *)_bjmodn20); _bjmodn20 = 0x0;	free((void *)_bjmodn48); _bjmodn48 = 0x0;
			free((void *)_bjmodn21); _bjmodn21 = 0x0;	free((void *)_bjmodn49); _bjmodn49 = 0x0;
			free((void *)_bjmodn22); _bjmodn22 = 0x0;	free((void *)_bjmodn50); _bjmodn50 = 0x0;
			free((void *)_bjmodn23); _bjmodn23 = 0x0;	free((void *)_bjmodn51); _bjmodn51 = 0x0;
			free((void *)_bjmodn24); _bjmodn24 = 0x0;	free((void *)_bjmodn52); _bjmodn52 = 0x0;
			free((void *)_bjmodn25); _bjmodn25 = 0x0;	free((void *)_bjmodn53); _bjmodn53 = 0x0;
			free((void *)_bjmodn26); _bjmodn26 = 0x0;	free((void *)_bjmodn54); _bjmodn54 = 0x0;
			free((void *)_bjmodn27); _bjmodn27 = 0x0;	free((void *)_bjmodn55); _bjmodn55 = 0x0;

			free((void *)_cy_r00); _cy_r00 = 0x0;		free((void *)_cy_r28); _cy_r28 = 0x0;
			free((void *)_cy_r01); _cy_r01 = 0x0;		free((void *)_cy_r29); _cy_r29 = 0x0;
			free((void *)_cy_r02); _cy_r02 = 0x0;		free((void *)_cy_r30); _cy_r30 = 0x0;
			free((void *)_cy_r03); _cy_r03 = 0x0;		free((void *)_cy_r31); _cy_r31 = 0x0;
			free((void *)_cy_r04); _cy_r04 = 0x0;		free((void *)_cy_r32); _cy_r32 = 0x0;
			free((void *)_cy_r05); _cy_r05 = 0x0;		free((void *)_cy_r33); _cy_r33 = 0x0;
			free((void *)_cy_r06); _cy_r06 = 0x0;		free((void *)_cy_r34); _cy_r34 = 0x0;
			free((void *)_cy_r07); _cy_r07 = 0x0;		free((void *)_cy_r35); _cy_r35 = 0x0;
			free((void *)_cy_r08); _cy_r08 = 0x0;		free((void *)_cy_r36); _cy_r36 = 0x0;
			free((void *)_cy_r09); _cy_r09 = 0x0;		free((void *)_cy_r37); _cy_r37 = 0x0;
			free((void *)_cy_r10); _cy_r10 = 0x0;		free((void *)_cy_r38); _cy_r38 = 0x0;
			free((void *)_cy_r11); _cy_r11 = 0x0;		free((void *)_cy_r39); _cy_r39 = 0x0;
			free((void *)_cy_r12); _cy_r12 = 0x0;		free((void *)_cy_r40); _cy_r40 = 0x0;
			free((void *)_cy_r13); _cy_r13 = 0x0;		free((void *)_cy_r41); _cy_r41 = 0x0;
			free((void *)_cy_r14); _cy_r14 = 0x0;		free((void *)_cy_r42); _cy_r42 = 0x0;
			free((void *)_cy_r15); _cy_r15 = 0x0;		free((void *)_cy_r43); _cy_r43 = 0x0;
			free((void *)_cy_r16); _cy_r16 = 0x0;		free((void *)_cy_r44); _cy_r44 = 0x0;
			free((void *)_cy_r17); _cy_r17 = 0x0;		free((void *)_cy_r45); _cy_r45 = 0x0;
			free((void *)_cy_r18); _cy_r18 = 0x0;		free((void *)_cy_r46); _cy_r46 = 0x0;
			free((void *)_cy_r19); _cy_r19 = 0x0;		free((void *)_cy_r47); _cy_r47 = 0x0;
			free((void *)_cy_r20); _cy_r20 = 0x0;		free((void *)_cy_r48); _cy_r48 = 0x0;
			free((void *)_cy_r21); _cy_r21 = 0x0;		free((void *)_cy_r49); _cy_r49 = 0x0;
			free((void *)_cy_r22); _cy_r22 = 0x0;		free((void *)_cy_r50); _cy_r50 = 0x0;
			free((void *)_cy_r23); _cy_r23 = 0x0;		free((void *)_cy_r51); _cy_r51 = 0x0;
			free((void *)_cy_r24); _cy_r24 = 0x0;		free((void *)_cy_r52); _cy_r52 = 0x0;
			free((void *)_cy_r25); _cy_r25 = 0x0;		free((void *)_cy_r53); _cy_r53 = 0x0;
			free((void *)_cy_r26); _cy_r26 = 0x0;		free((void *)_cy_r54); _cy_r54 = 0x0;
			free((void *)_cy_r27); _cy_r27 = 0x0;		free((void *)_cy_r55); _cy_r55 = 0x0;

			free((void *)_cy_i00); _cy_i00 = 0x0;		free((void *)_cy_i28); _cy_i28 = 0x0;
			free((void *)_cy_i01); _cy_i01 = 0x0;		free((void *)_cy_i29); _cy_i29 = 0x0;
			free((void *)_cy_i02); _cy_i02 = 0x0;		free((void *)_cy_i30); _cy_i30 = 0x0;
			free((void *)_cy_i03); _cy_i03 = 0x0;		free((void *)_cy_i31); _cy_i31 = 0x0;
			free((void *)_cy_i04); _cy_i04 = 0x0;		free((void *)_cy_i32); _cy_i32 = 0x0;
			free((void *)_cy_i05); _cy_i05 = 0x0;		free((void *)_cy_i33); _cy_i33 = 0x0;
			free((void *)_cy_i06); _cy_i06 = 0x0;		free((void *)_cy_i34); _cy_i34 = 0x0;
			free((void *)_cy_i07); _cy_i07 = 0x0;		free((void *)_cy_i35); _cy_i35 = 0x0;
			free((void *)_cy_i08); _cy_i08 = 0x0;		free((void *)_cy_i36); _cy_i36 = 0x0;
			free((void *)_cy_i09); _cy_i09 = 0x0;		free((void *)_cy_i37); _cy_i37 = 0x0;
			free((void *)_cy_i10); _cy_i10 = 0x0;		free((void *)_cy_i38); _cy_i38 = 0x0;
			free((void *)_cy_i11); _cy_i11 = 0x0;		free((void *)_cy_i39); _cy_i39 = 0x0;
			free((void *)_cy_i12); _cy_i12 = 0x0;		free((void *)_cy_i40); _cy_i40 = 0x0;
			free((void *)_cy_i13); _cy_i13 = 0x0;		free((void *)_cy_i41); _cy_i41 = 0x0;
			free((void *)_cy_i14); _cy_i14 = 0x0;		free((void *)_cy_i42); _cy_i42 = 0x0;
			free((void *)_cy_i15); _cy_i15 = 0x0;		free((void *)_cy_i43); _cy_i43 = 0x0;
			free((void *)_cy_i16); _cy_i16 = 0x0;		free((void *)_cy_i44); _cy_i44 = 0x0;
			free((void *)_cy_i17); _cy_i17 = 0x0;		free((void *)_cy_i45); _cy_i45 = 0x0;
			free((void *)_cy_i18); _cy_i18 = 0x0;		free((void *)_cy_i46); _cy_i46 = 0x0;
			free((void *)_cy_i19); _cy_i19 = 0x0;		free((void *)_cy_i47); _cy_i47 = 0x0;
			free((void *)_cy_i20); _cy_i20 = 0x0;		free((void *)_cy_i48); _cy_i48 = 0x0;
			free((void *)_cy_i21); _cy_i21 = 0x0;		free((void *)_cy_i49); _cy_i49 = 0x0;
			free((void *)_cy_i22); _cy_i22 = 0x0;		free((void *)_cy_i50); _cy_i50 = 0x0;
			free((void *)_cy_i23); _cy_i23 = 0x0;		free((void *)_cy_i51); _cy_i51 = 0x0;
			free((void *)_cy_i24); _cy_i24 = 0x0;		free((void *)_cy_i52); _cy_i52 = 0x0;
			free((void *)_cy_i25); _cy_i25 = 0x0;		free((void *)_cy_i53); _cy_i53 = 0x0;
			free((void *)_cy_i26); _cy_i26 = 0x0;		free((void *)_cy_i54); _cy_i54 = 0x0;
			free((void *)_cy_i27); _cy_i27 = 0x0;		free((void *)_cy_i55); _cy_i55 = 0x0;

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
		_bjmodn28	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn28== 0x0);
		_bjmodn29	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn29== 0x0);
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
		_bjmodn40	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn40== 0x0);
		_bjmodn41	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn41== 0x0);
		_bjmodn42	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn42== 0x0);
		_bjmodn43	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn43== 0x0);
		_bjmodn44	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn44== 0x0);
		_bjmodn45	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn45== 0x0);
		_bjmodn46	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn46== 0x0);
		_bjmodn47	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn47== 0x0);
		_bjmodn48	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn48== 0x0);
		_bjmodn49	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn49== 0x0);
		_bjmodn50	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn50== 0x0);
		_bjmodn51	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn51== 0x0);
		_bjmodn52	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn52== 0x0);
		_bjmodn53	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn53== 0x0);
		_bjmodn54	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn54== 0x0);
		_bjmodn55	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn55== 0x0);
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
		_cy_r28	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r28== 0x0);
		_cy_r29	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r29== 0x0);
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
		_cy_r40	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r40== 0x0);
		_cy_r41	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r41== 0x0);
		_cy_r42	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r42== 0x0);
		_cy_r43	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r43== 0x0);
		_cy_r44	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r44== 0x0);
		_cy_r45	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r45== 0x0);
		_cy_r46	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r46== 0x0);
		_cy_r47	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r47== 0x0);
		_cy_r48	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r48== 0x0);
		_cy_r49	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r49== 0x0);
		_cy_r50	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r50== 0x0);
		_cy_r51	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r51== 0x0);
		_cy_r52	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r52== 0x0);
		_cy_r53	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r53== 0x0);
		_cy_r54	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r54== 0x0);
		_cy_r55	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r55== 0x0);

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
		_cy_i28	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i28== 0x0);
		_cy_i29	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i29== 0x0);
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
		_cy_i40	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i40== 0x0);
		_cy_i41	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i41== 0x0);
		_cy_i42	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i42== 0x0);
		_cy_i43	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i43== 0x0);
		_cy_i44	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i44== 0x0);
		_cy_i45	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i45== 0x0);
		_cy_i46	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i46== 0x0);
		_cy_i47	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i47== 0x0);
		_cy_i48	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i48== 0x0);
		_cy_i49	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i49== 0x0);
		_cy_i50	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i50== 0x0);
		_cy_i51	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i51== 0x0);
		_cy_i52	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i52== 0x0);
		_cy_i53	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i53== 0x0);
		_cy_i54	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i54== 0x0);
		_cy_i55	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i55== 0x0);

		_maxerr	= (double *)malloc(j);	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/radix-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
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

	}	/* endif(first_entry) */

/*...The radix-56 final DIT pass is here.	*/

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
		_cy_r28[ithread] = 0;	_cy_i28[ithread] = 0;
		_cy_r29[ithread] = 0;	_cy_i29[ithread] = 0;
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
		_cy_r40[ithread] = 0;	_cy_i40[ithread] = 0;
		_cy_r41[ithread] = 0;	_cy_i41[ithread] = 0;
		_cy_r42[ithread] = 0;	_cy_i42[ithread] = 0;
		_cy_r43[ithread] = 0;	_cy_i43[ithread] = 0;
		_cy_r44[ithread] = 0;	_cy_i44[ithread] = 0;
		_cy_r45[ithread] = 0;	_cy_i45[ithread] = 0;
		_cy_r46[ithread] = 0;	_cy_i46[ithread] = 0;
		_cy_r47[ithread] = 0;	_cy_i47[ithread] = 0;
		_cy_r48[ithread] = 0;	_cy_i48[ithread] = 0;
		_cy_r49[ithread] = 0;	_cy_i49[ithread] = 0;
		_cy_r50[ithread] = 0;	_cy_i50[ithread] = 0;
		_cy_r51[ithread] = 0;	_cy_i51[ithread] = 0;
		_cy_r52[ithread] = 0;	_cy_i52[ithread] = 0;
		_cy_r53[ithread] = 0;	_cy_i53[ithread] = 0;
		_cy_r54[ithread] = 0;	_cy_i54[ithread] = 0;
		_cy_r55[ithread] = 0;	_cy_i55[ithread] = 0;
	}
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy_r00[      0] = -2;
	}

	*fracmax=0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = 0.0;
	}

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
		MOD_ADD32(_bjmodn27[ithread], j, n, _bjmodn28[ithread]);
		MOD_ADD32(_bjmodn28[ithread], j, n, _bjmodn29[ithread]);
		MOD_ADD32(_bjmodn29[ithread], j, n, _bjmodn30[ithread]);
		MOD_ADD32(_bjmodn30[ithread], j, n, _bjmodn31[ithread]);
		MOD_ADD32(_bjmodn31[ithread], j, n, _bjmodn32[ithread]);
		MOD_ADD32(_bjmodn32[ithread], j, n, _bjmodn33[ithread]);
		MOD_ADD32(_bjmodn33[ithread], j, n, _bjmodn34[ithread]);
		MOD_ADD32(_bjmodn34[ithread], j, n, _bjmodn35[ithread]);
		MOD_ADD32(_bjmodn35[ithread], j, n, _bjmodn36[ithread]);
		MOD_ADD32(_bjmodn36[ithread], j, n, _bjmodn37[ithread]);
		MOD_ADD32(_bjmodn37[ithread], j, n, _bjmodn38[ithread]);
		MOD_ADD32(_bjmodn38[ithread], j, n, _bjmodn39[ithread]);
		MOD_ADD32(_bjmodn39[ithread], j, n, _bjmodn40[ithread]);
		MOD_ADD32(_bjmodn40[ithread], j, n, _bjmodn41[ithread]);
		MOD_ADD32(_bjmodn41[ithread], j, n, _bjmodn42[ithread]);
		MOD_ADD32(_bjmodn42[ithread], j, n, _bjmodn43[ithread]);
		MOD_ADD32(_bjmodn43[ithread], j, n, _bjmodn44[ithread]);
		MOD_ADD32(_bjmodn44[ithread], j, n, _bjmodn45[ithread]);
		MOD_ADD32(_bjmodn45[ithread], j, n, _bjmodn46[ithread]);
		MOD_ADD32(_bjmodn46[ithread], j, n, _bjmodn47[ithread]);
		MOD_ADD32(_bjmodn47[ithread], j, n, _bjmodn48[ithread]);
		MOD_ADD32(_bjmodn48[ithread], j, n, _bjmodn49[ithread]);
		MOD_ADD32(_bjmodn49[ithread], j, n, _bjmodn50[ithread]);
		MOD_ADD32(_bjmodn50[ithread], j, n, _bjmodn51[ithread]);
		MOD_ADD32(_bjmodn51[ithread], j, n, _bjmodn52[ithread]);
		MOD_ADD32(_bjmodn52[ithread], j, n, _bjmodn53[ithread]);
		MOD_ADD32(_bjmodn53[ithread], j, n, _bjmodn54[ithread]);
		MOD_ADD32(_bjmodn54[ithread], j, n, _bjmodn55[ithread]);

		// Every (odd_radix)th bjmodn initializer needs to be forced-to-bigword in fermat-mod DWT case:
		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		{
			/* Start this value off at N in Fermat-mod case, so (bjmodn >= sw) check in
			fermat_carry_norm_errcheck (cf. carry.h) yields a bigword (i == 1) for j= 0:
			*/
			_bjmodn00[ithread] = n;
			_bjmodn07[ithread] = n;
			_bjmodn14[ithread] = n;
			_bjmodn21[ithread] = n;
			_bjmodn28[ithread] = n;
			_bjmodn35[ithread] = n;
			_bjmodn42[ithread] = n;
			_bjmodn49[ithread] = n;
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

		tmp = half_arr + odd_radix*2;	/* Put the base-mini-arrays right after the weights */

	  #ifdef USE_AVX
		// Each transposed-data quartet in the AVX carry macro needs linearly incrementing bs_arr data (mod odd_radix);
		// Need all [odd_radix] possible such length-4 index subsequences, which will be accessed via their head element
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

/*	fprintf(stderr, "radix56_ditN_cy_dif1: wts_idx_incr = %d\n", wts_idx_incr);*/

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
		ASSERT(HERE, tdat[ithread].wts_idx_inc2 == wts_idx_inc2, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].r00r == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].r00r;
		ASSERT(HERE, ((tmp + half_arr_offset56-1)->d0 == crnd && (tmp + half_arr_offset56-1)->d1 == crnd), "thread-local memcheck failed!");
	#endif
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX
			// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
			dtmp = (tmp + half_arr_offset56+40)->d0 * (tmp + half_arr_offset56+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp + half_arr_offset56+40)->d1 * (tmp + half_arr_offset56+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#elif defined(USE_SSE2)
			dtmp = (tmp + half_arr_offset56+10)->d0 * (tmp + half_arr_offset56+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp + half_arr_offset56+10)->d1 * (tmp + half_arr_offset56+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
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
			tdat[ithread].bjmodn28 = _bjmodn28[ithread];
			tdat[ithread].bjmodn29 = _bjmodn29[ithread];
			tdat[ithread].bjmodn30 = _bjmodn30[ithread];
			tdat[ithread].bjmodn31 = _bjmodn31[ithread];
			tdat[ithread].bjmodn32 = _bjmodn32[ithread];
			tdat[ithread].bjmodn33 = _bjmodn33[ithread];
			tdat[ithread].bjmodn34 = _bjmodn34[ithread];
			tdat[ithread].bjmodn35 = _bjmodn35[ithread];
			tdat[ithread].bjmodn36 = _bjmodn36[ithread];
			tdat[ithread].bjmodn37 = _bjmodn37[ithread];
			tdat[ithread].bjmodn38 = _bjmodn38[ithread];
			tdat[ithread].bjmodn39 = _bjmodn39[ithread];
			tdat[ithread].bjmodn40 = _bjmodn40[ithread];
			tdat[ithread].bjmodn41 = _bjmodn41[ithread];
			tdat[ithread].bjmodn42 = _bjmodn42[ithread];
			tdat[ithread].bjmodn43 = _bjmodn43[ithread];
			tdat[ithread].bjmodn44 = _bjmodn44[ithread];
			tdat[ithread].bjmodn45 = _bjmodn45[ithread];
			tdat[ithread].bjmodn46 = _bjmodn46[ithread];
			tdat[ithread].bjmodn47 = _bjmodn47[ithread];
			tdat[ithread].bjmodn48 = _bjmodn48[ithread];
			tdat[ithread].bjmodn49 = _bjmodn49[ithread];
			tdat[ithread].bjmodn50 = _bjmodn50[ithread];
			tdat[ithread].bjmodn51 = _bjmodn51[ithread];
			tdat[ithread].bjmodn52 = _bjmodn52[ithread];
			tdat[ithread].bjmodn53 = _bjmodn53[ithread];
			tdat[ithread].bjmodn54 = _bjmodn54[ithread];
			tdat[ithread].bjmodn55 = _bjmodn55[ithread];
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
			tdat[ithread].cy_r28 = _cy_r28[ithread];
			tdat[ithread].cy_r29 = _cy_r29[ithread];
			tdat[ithread].cy_r30 = _cy_r30[ithread];
			tdat[ithread].cy_r31 = _cy_r31[ithread];
			tdat[ithread].cy_r32 = _cy_r32[ithread];
			tdat[ithread].cy_r33 = _cy_r33[ithread];
			tdat[ithread].cy_r34 = _cy_r34[ithread];
			tdat[ithread].cy_r35 = _cy_r35[ithread];
			tdat[ithread].cy_r36 = _cy_r36[ithread];
			tdat[ithread].cy_r37 = _cy_r37[ithread];
			tdat[ithread].cy_r38 = _cy_r38[ithread];
			tdat[ithread].cy_r39 = _cy_r39[ithread];
			tdat[ithread].cy_r40 = _cy_r40[ithread];
			tdat[ithread].cy_r41 = _cy_r41[ithread];
			tdat[ithread].cy_r42 = _cy_r42[ithread];
			tdat[ithread].cy_r43 = _cy_r43[ithread];
			tdat[ithread].cy_r44 = _cy_r44[ithread];
			tdat[ithread].cy_r45 = _cy_r45[ithread];
			tdat[ithread].cy_r46 = _cy_r46[ithread];
			tdat[ithread].cy_r47 = _cy_r47[ithread];
			tdat[ithread].cy_r48 = _cy_r48[ithread];
			tdat[ithread].cy_r49 = _cy_r49[ithread];
			tdat[ithread].cy_r50 = _cy_r50[ithread];
			tdat[ithread].cy_r51 = _cy_r51[ithread];
			tdat[ithread].cy_r52 = _cy_r52[ithread];
			tdat[ithread].cy_r53 = _cy_r53[ithread];
			tdat[ithread].cy_r54 = _cy_r54[ithread];
			tdat[ithread].cy_r55 = _cy_r55[ithread];
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
		#ifdef USE_SSE2
			dtmp = (tmp + half_arr_offset56)->d0 * (tmp + half_arr_offset56+odd_radix)->d0;	ASSERT(HERE, fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp + half_arr_offset56)->d1 * (tmp + half_arr_offset56+odd_radix)->d1;	ASSERT(HERE, fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
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
			tdat[ithread].cy_r28 = _cy_r28[ithread];	tdat[ithread].cy_i28 = _cy_i28[ithread];
			tdat[ithread].cy_r29 = _cy_r29[ithread];	tdat[ithread].cy_i29 = _cy_i29[ithread];
			tdat[ithread].cy_r30 = _cy_r30[ithread];	tdat[ithread].cy_i30 = _cy_i30[ithread];
			tdat[ithread].cy_r31 = _cy_r31[ithread];	tdat[ithread].cy_i31 = _cy_i31[ithread];
			tdat[ithread].cy_r32 = _cy_r32[ithread];	tdat[ithread].cy_i32 = _cy_i32[ithread];
			tdat[ithread].cy_r33 = _cy_r33[ithread];	tdat[ithread].cy_i33 = _cy_i33[ithread];
			tdat[ithread].cy_r34 = _cy_r34[ithread];	tdat[ithread].cy_i34 = _cy_i34[ithread];
			tdat[ithread].cy_r35 = _cy_r35[ithread];	tdat[ithread].cy_i35 = _cy_i35[ithread];
			tdat[ithread].cy_r36 = _cy_r36[ithread];	tdat[ithread].cy_i36 = _cy_i36[ithread];
			tdat[ithread].cy_r37 = _cy_r37[ithread];	tdat[ithread].cy_i37 = _cy_i37[ithread];
			tdat[ithread].cy_r38 = _cy_r38[ithread];	tdat[ithread].cy_i38 = _cy_i38[ithread];
			tdat[ithread].cy_r39 = _cy_r39[ithread];	tdat[ithread].cy_i39 = _cy_i39[ithread];
			tdat[ithread].cy_r40 = _cy_r40[ithread];	tdat[ithread].cy_i40 = _cy_i40[ithread];
			tdat[ithread].cy_r41 = _cy_r41[ithread];	tdat[ithread].cy_i41 = _cy_i41[ithread];
			tdat[ithread].cy_r42 = _cy_r42[ithread];	tdat[ithread].cy_i42 = _cy_i42[ithread];
			tdat[ithread].cy_r43 = _cy_r43[ithread];	tdat[ithread].cy_i43 = _cy_i43[ithread];
			tdat[ithread].cy_r44 = _cy_r44[ithread];	tdat[ithread].cy_i44 = _cy_i44[ithread];
			tdat[ithread].cy_r45 = _cy_r45[ithread];	tdat[ithread].cy_i45 = _cy_i45[ithread];
			tdat[ithread].cy_r46 = _cy_r46[ithread];	tdat[ithread].cy_i46 = _cy_i46[ithread];
			tdat[ithread].cy_r47 = _cy_r47[ithread];	tdat[ithread].cy_i47 = _cy_i47[ithread];
			tdat[ithread].cy_r48 = _cy_r48[ithread];	tdat[ithread].cy_i48 = _cy_i48[ithread];
			tdat[ithread].cy_r49 = _cy_r49[ithread];	tdat[ithread].cy_i49 = _cy_i49[ithread];
			tdat[ithread].cy_r50 = _cy_r50[ithread];	tdat[ithread].cy_i50 = _cy_i50[ithread];
			tdat[ithread].cy_r51 = _cy_r51[ithread];	tdat[ithread].cy_i51 = _cy_i51[ithread];
			tdat[ithread].cy_r52 = _cy_r52[ithread];	tdat[ithread].cy_i52 = _cy_i52[ithread];
			tdat[ithread].cy_r53 = _cy_r53[ithread];	tdat[ithread].cy_i53 = _cy_i53[ithread];
			tdat[ithread].cy_r54 = _cy_r54[ithread];	tdat[ithread].cy_i54 = _cy_i54[ithread];
			tdat[ithread].cy_r55 = _cy_r55[ithread];	tdat[ithread].cy_i55 = _cy_i55[ithread];
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

	#ifdef USE_SSE2
		*bjmodn00 = _bjmodn00[ithread];		*bjmodn28 = _bjmodn28[ithread];
		*bjmodn01 = _bjmodn01[ithread];		*bjmodn29 = _bjmodn29[ithread];
		*bjmodn02 = _bjmodn02[ithread];		*bjmodn30 = _bjmodn30[ithread];
		*bjmodn03 = _bjmodn03[ithread];		*bjmodn31 = _bjmodn31[ithread];
		*bjmodn04 = _bjmodn04[ithread];		*bjmodn32 = _bjmodn32[ithread];
		*bjmodn05 = _bjmodn05[ithread];		*bjmodn33 = _bjmodn33[ithread];
		*bjmodn06 = _bjmodn06[ithread];		*bjmodn34 = _bjmodn34[ithread];
		*bjmodn07 = _bjmodn07[ithread];		*bjmodn35 = _bjmodn35[ithread];
		*bjmodn08 = _bjmodn08[ithread];		*bjmodn36 = _bjmodn36[ithread];
		*bjmodn09 = _bjmodn09[ithread];		*bjmodn37 = _bjmodn37[ithread];
		*bjmodn10 = _bjmodn10[ithread];		*bjmodn38 = _bjmodn38[ithread];
		*bjmodn11 = _bjmodn11[ithread];		*bjmodn39 = _bjmodn39[ithread];
		*bjmodn12 = _bjmodn12[ithread];		*bjmodn40 = _bjmodn40[ithread];
		*bjmodn13 = _bjmodn13[ithread];		*bjmodn41 = _bjmodn41[ithread];
		*bjmodn14 = _bjmodn14[ithread];		*bjmodn42 = _bjmodn42[ithread];
		*bjmodn15 = _bjmodn15[ithread];		*bjmodn43 = _bjmodn43[ithread];
		*bjmodn16 = _bjmodn16[ithread];		*bjmodn44 = _bjmodn44[ithread];
		*bjmodn17 = _bjmodn17[ithread];		*bjmodn45 = _bjmodn45[ithread];
		*bjmodn18 = _bjmodn18[ithread];		*bjmodn46 = _bjmodn46[ithread];
		*bjmodn19 = _bjmodn19[ithread];		*bjmodn47 = _bjmodn47[ithread];
		*bjmodn20 = _bjmodn20[ithread];		*bjmodn48 = _bjmodn48[ithread];
		*bjmodn21 = _bjmodn21[ithread];		*bjmodn49 = _bjmodn49[ithread];
		*bjmodn22 = _bjmodn22[ithread];		*bjmodn50 = _bjmodn50[ithread];
		*bjmodn23 = _bjmodn23[ithread];		*bjmodn51 = _bjmodn51[ithread];
		*bjmodn24 = _bjmodn24[ithread];		*bjmodn52 = _bjmodn52[ithread];
		*bjmodn25 = _bjmodn25[ithread];		*bjmodn53 = _bjmodn53[ithread];
		*bjmodn26 = _bjmodn26[ithread];		*bjmodn54 = _bjmodn54[ithread];
		*bjmodn27 = _bjmodn27[ithread];		*bjmodn55 = _bjmodn55[ithread];
	#else
		bjmodn00 = _bjmodn00[ithread];		bjmodn28 = _bjmodn28[ithread];
		bjmodn01 = _bjmodn01[ithread];		bjmodn29 = _bjmodn29[ithread];
		bjmodn02 = _bjmodn02[ithread];		bjmodn30 = _bjmodn30[ithread];
		bjmodn03 = _bjmodn03[ithread];		bjmodn31 = _bjmodn31[ithread];
		bjmodn04 = _bjmodn04[ithread];		bjmodn32 = _bjmodn32[ithread];
		bjmodn05 = _bjmodn05[ithread];		bjmodn33 = _bjmodn33[ithread];
		bjmodn06 = _bjmodn06[ithread];		bjmodn34 = _bjmodn34[ithread];
		bjmodn07 = _bjmodn07[ithread];		bjmodn35 = _bjmodn35[ithread];
		bjmodn08 = _bjmodn08[ithread];		bjmodn36 = _bjmodn36[ithread];
		bjmodn09 = _bjmodn09[ithread];		bjmodn37 = _bjmodn37[ithread];
		bjmodn10 = _bjmodn10[ithread];		bjmodn38 = _bjmodn38[ithread];
		bjmodn11 = _bjmodn11[ithread];		bjmodn39 = _bjmodn39[ithread];
		bjmodn12 = _bjmodn12[ithread];		bjmodn40 = _bjmodn40[ithread];
		bjmodn13 = _bjmodn13[ithread];		bjmodn41 = _bjmodn41[ithread];
		bjmodn14 = _bjmodn14[ithread];		bjmodn42 = _bjmodn42[ithread];
		bjmodn15 = _bjmodn15[ithread];		bjmodn43 = _bjmodn43[ithread];
		bjmodn16 = _bjmodn16[ithread];		bjmodn44 = _bjmodn44[ithread];
		bjmodn17 = _bjmodn17[ithread];		bjmodn45 = _bjmodn45[ithread];
		bjmodn18 = _bjmodn18[ithread];		bjmodn46 = _bjmodn46[ithread];
		bjmodn19 = _bjmodn19[ithread];		bjmodn47 = _bjmodn47[ithread];
		bjmodn20 = _bjmodn20[ithread];		bjmodn48 = _bjmodn48[ithread];
		bjmodn21 = _bjmodn21[ithread];		bjmodn49 = _bjmodn49[ithread];
		bjmodn22 = _bjmodn22[ithread];		bjmodn50 = _bjmodn50[ithread];
		bjmodn23 = _bjmodn23[ithread];		bjmodn51 = _bjmodn51[ithread];
		bjmodn24 = _bjmodn24[ithread];		bjmodn52 = _bjmodn52[ithread];
		bjmodn25 = _bjmodn25[ithread];		bjmodn53 = _bjmodn53[ithread];
		bjmodn26 = _bjmodn26[ithread];		bjmodn54 = _bjmodn54[ithread];
		bjmodn27 = _bjmodn27[ithread];		bjmodn55 = _bjmodn55[ithread];
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
			cy_r28->d0 = _cy_r28[ithread];	cy_r28->d1 = _cy_r29[ithread];	cy_r28->d2 = _cy_r30[ithread];	cy_r28->d3 = _cy_r31[ithread];
			cy_r32->d0 = _cy_r32[ithread];	cy_r32->d1 = _cy_r33[ithread];	cy_r32->d2 = _cy_r34[ithread];	cy_r32->d3 = _cy_r35[ithread];
			cy_r36->d0 = _cy_r36[ithread];	cy_r36->d1 = _cy_r37[ithread];	cy_r36->d2 = _cy_r38[ithread];	cy_r36->d3 = _cy_r39[ithread];
			cy_r40->d0 = _cy_r40[ithread];	cy_r40->d1 = _cy_r41[ithread];	cy_r40->d2 = _cy_r42[ithread];	cy_r40->d3 = _cy_r43[ithread];
			cy_r44->d0 = _cy_r44[ithread];	cy_r44->d1 = _cy_r45[ithread];	cy_r44->d2 = _cy_r46[ithread];	cy_r44->d3 = _cy_r47[ithread];
			cy_r48->d0 = _cy_r48[ithread];	cy_r48->d1 = _cy_r49[ithread];	cy_r48->d2 = _cy_r50[ithread];	cy_r48->d3 = _cy_r51[ithread];
			cy_r52->d0 = _cy_r52[ithread];	cy_r52->d1 = _cy_r53[ithread];	cy_r52->d2 = _cy_r54[ithread];	cy_r52->d3 = _cy_r55[ithread];
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
			cy_r28->d0 = _cy_r28[ithread];	cy_r28->d1 = _cy_r29[ithread];
			cy_r30->d0 = _cy_r30[ithread];	cy_r30->d1 = _cy_r31[ithread];
			cy_r32->d0 = _cy_r32[ithread];	cy_r32->d1 = _cy_r33[ithread];
			cy_r34->d0 = _cy_r34[ithread];	cy_r34->d1 = _cy_r35[ithread];
			cy_r36->d0 = _cy_r36[ithread];	cy_r36->d1 = _cy_r37[ithread];
			cy_r38->d0 = _cy_r38[ithread];	cy_r38->d1 = _cy_r39[ithread];
			cy_r40->d0 = _cy_r40[ithread];	cy_r40->d1 = _cy_r41[ithread];
			cy_r42->d0 = _cy_r42[ithread];	cy_r42->d1 = _cy_r43[ithread];
			cy_r44->d0 = _cy_r44[ithread];	cy_r44->d1 = _cy_r45[ithread];
			cy_r46->d0 = _cy_r46[ithread];	cy_r46->d1 = _cy_r47[ithread];
			cy_r48->d0 = _cy_r48[ithread];	cy_r48->d1 = _cy_r49[ithread];
			cy_r50->d0 = _cy_r50[ithread];	cy_r50->d1 = _cy_r51[ithread];
			cy_r52->d0 = _cy_r52[ithread];	cy_r52->d1 = _cy_r53[ithread];
			cy_r54->d0 = _cy_r54[ithread];	cy_r54->d1 = _cy_r55[ithread];
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
			cy_r28 = _cy_r28[ithread];
			cy_r29 = _cy_r29[ithread];
			cy_r30 = _cy_r30[ithread];
			cy_r31 = _cy_r31[ithread];
			cy_r32 = _cy_r32[ithread];
			cy_r33 = _cy_r33[ithread];
			cy_r34 = _cy_r34[ithread];
			cy_r35 = _cy_r35[ithread];
			cy_r36 = _cy_r36[ithread];
			cy_r37 = _cy_r37[ithread];
			cy_r38 = _cy_r38[ithread];
			cy_r39 = _cy_r39[ithread];
			cy_r40 = _cy_r40[ithread];
			cy_r41 = _cy_r41[ithread];
			cy_r42 = _cy_r42[ithread];
			cy_r43 = _cy_r43[ithread];
			cy_r44 = _cy_r44[ithread];
			cy_r45 = _cy_r45[ithread];
			cy_r46 = _cy_r46[ithread];
			cy_r47 = _cy_r47[ithread];
			cy_r48 = _cy_r48[ithread];
			cy_r49 = _cy_r49[ithread];
			cy_r50 = _cy_r50[ithread];
			cy_r51 = _cy_r51[ithread];
			cy_r52 = _cy_r52[ithread];
			cy_r53 = _cy_r53[ithread];
			cy_r54 = _cy_r54[ithread];
			cy_r55 = _cy_r55[ithread];
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
			cy_r28->d0 = _cy_r28[ithread];	cy_r28->d1 = _cy_r29[ithread];	cy_r28->d2 = _cy_r30[ithread];	cy_r28->d3 = _cy_r31[ithread];
			cy_r32->d0 = _cy_r32[ithread];	cy_r32->d1 = _cy_r33[ithread];	cy_r32->d2 = _cy_r34[ithread];	cy_r32->d3 = _cy_r35[ithread];
			cy_r36->d0 = _cy_r36[ithread];	cy_r36->d1 = _cy_r37[ithread];	cy_r36->d2 = _cy_r38[ithread];	cy_r36->d3 = _cy_r39[ithread];
			cy_r40->d0 = _cy_r40[ithread];	cy_r40->d1 = _cy_r41[ithread];	cy_r40->d2 = _cy_r42[ithread];	cy_r40->d3 = _cy_r43[ithread];
			cy_r44->d0 = _cy_r44[ithread];	cy_r44->d1 = _cy_r45[ithread];	cy_r44->d2 = _cy_r46[ithread];	cy_r44->d3 = _cy_r47[ithread];
			cy_r48->d0 = _cy_r48[ithread];	cy_r48->d1 = _cy_r49[ithread];	cy_r48->d2 = _cy_r50[ithread];	cy_r48->d3 = _cy_r51[ithread];
			cy_r52->d0 = _cy_r52[ithread];	cy_r52->d1 = _cy_r53[ithread];	cy_r52->d2 = _cy_r54[ithread];	cy_r52->d3 = _cy_r55[ithread];

			cy_i00->d0 = _cy_i00[ithread];	cy_i00->d1 = _cy_i01[ithread];	cy_i00->d2 = _cy_i02[ithread];	cy_i00->d3 = _cy_i03[ithread];
			cy_i04->d0 = _cy_i04[ithread];	cy_i04->d1 = _cy_i05[ithread];	cy_i04->d2 = _cy_i06[ithread];	cy_i04->d3 = _cy_i07[ithread];
			cy_i08->d0 = _cy_i08[ithread];	cy_i08->d1 = _cy_i09[ithread];	cy_i08->d2 = _cy_i10[ithread];	cy_i08->d3 = _cy_i11[ithread];
			cy_i12->d0 = _cy_i12[ithread];	cy_i12->d1 = _cy_i13[ithread];	cy_i12->d2 = _cy_i14[ithread];	cy_i12->d3 = _cy_i15[ithread];
			cy_i16->d0 = _cy_i16[ithread];	cy_i16->d1 = _cy_i17[ithread];	cy_i16->d2 = _cy_i18[ithread];	cy_i16->d3 = _cy_i19[ithread];
			cy_i20->d0 = _cy_i20[ithread];	cy_i20->d1 = _cy_i21[ithread];	cy_i20->d2 = _cy_i22[ithread];	cy_i20->d3 = _cy_i23[ithread];
			cy_i24->d0 = _cy_i24[ithread];	cy_i24->d1 = _cy_i25[ithread];	cy_i24->d2 = _cy_i26[ithread];	cy_i24->d3 = _cy_i27[ithread];
			cy_i28->d0 = _cy_i28[ithread];	cy_i28->d1 = _cy_i29[ithread];	cy_i28->d2 = _cy_i30[ithread];	cy_i28->d3 = _cy_i31[ithread];
			cy_i32->d0 = _cy_i32[ithread];	cy_i32->d1 = _cy_i33[ithread];	cy_i32->d2 = _cy_i34[ithread];	cy_i32->d3 = _cy_i35[ithread];
			cy_i36->d0 = _cy_i36[ithread];	cy_i36->d1 = _cy_i37[ithread];	cy_i36->d2 = _cy_i38[ithread];	cy_i36->d3 = _cy_i39[ithread];
			cy_i40->d0 = _cy_i40[ithread];	cy_i40->d1 = _cy_i41[ithread];	cy_i40->d2 = _cy_i42[ithread];	cy_i40->d3 = _cy_i43[ithread];
			cy_i44->d0 = _cy_i44[ithread];	cy_i44->d1 = _cy_i45[ithread];	cy_i44->d2 = _cy_i46[ithread];	cy_i44->d3 = _cy_i47[ithread];
			cy_i48->d0 = _cy_i48[ithread];	cy_i48->d1 = _cy_i49[ithread];	cy_i48->d2 = _cy_i50[ithread];	cy_i48->d3 = _cy_i51[ithread];
			cy_i52->d0 = _cy_i52[ithread];	cy_i52->d1 = _cy_i53[ithread];	cy_i52->d2 = _cy_i54[ithread];	cy_i52->d3 = _cy_i55[ithread];
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
			cy_r28->d0 = _cy_r14[ithread];	cy_r28->d1 = _cy_i14[ithread];
			cy_r30->d0 = _cy_r15[ithread];	cy_r30->d1 = _cy_i15[ithread];
			cy_r32->d0 = _cy_r16[ithread];	cy_r32->d1 = _cy_i16[ithread];
			cy_r34->d0 = _cy_r17[ithread];	cy_r34->d1 = _cy_i17[ithread];
			cy_r36->d0 = _cy_r18[ithread];	cy_r36->d1 = _cy_i18[ithread];
			cy_r38->d0 = _cy_r19[ithread];	cy_r38->d1 = _cy_i19[ithread];
			cy_r40->d0 = _cy_r20[ithread];	cy_r40->d1 = _cy_i20[ithread];
			cy_r42->d0 = _cy_r21[ithread];	cy_r42->d1 = _cy_i21[ithread];
			cy_r44->d0 = _cy_r22[ithread];	cy_r44->d1 = _cy_i22[ithread];
			cy_r46->d0 = _cy_r23[ithread];	cy_r46->d1 = _cy_i23[ithread];
			cy_r48->d0 = _cy_r24[ithread];	cy_r48->d1 = _cy_i24[ithread];
			cy_r50->d0 = _cy_r25[ithread];	cy_r50->d1 = _cy_i25[ithread];
			cy_r52->d0 = _cy_r26[ithread];	cy_r52->d1 = _cy_i26[ithread];
			cy_r54->d0 = _cy_r27[ithread];	cy_r54->d1 = _cy_i27[ithread];
			cy_i00->d0 = _cy_r28[ithread];	cy_i00->d1 = _cy_i28[ithread];
			cy_i02->d0 = _cy_r29[ithread];	cy_i02->d1 = _cy_i29[ithread];
			cy_i04->d0 = _cy_r30[ithread];	cy_i04->d1 = _cy_i30[ithread];
			cy_i06->d0 = _cy_r31[ithread];	cy_i06->d1 = _cy_i31[ithread];
			cy_i08->d0 = _cy_r32[ithread];	cy_i08->d1 = _cy_i32[ithread];
			cy_i10->d0 = _cy_r33[ithread];	cy_i10->d1 = _cy_i33[ithread];
			cy_i12->d0 = _cy_r34[ithread];	cy_i12->d1 = _cy_i34[ithread];
			cy_i14->d0 = _cy_r35[ithread];	cy_i14->d1 = _cy_i35[ithread];
			cy_i16->d0 = _cy_r36[ithread];	cy_i16->d1 = _cy_i36[ithread];
			cy_i18->d0 = _cy_r37[ithread];	cy_i18->d1 = _cy_i37[ithread];
			cy_i20->d0 = _cy_r38[ithread];	cy_i20->d1 = _cy_i38[ithread];
			cy_i22->d0 = _cy_r39[ithread];	cy_i22->d1 = _cy_i39[ithread];
			cy_i24->d0 = _cy_r40[ithread];	cy_i24->d1 = _cy_i40[ithread];
			cy_i26->d0 = _cy_r41[ithread];	cy_i26->d1 = _cy_i41[ithread];
			cy_i28->d0 = _cy_r42[ithread];	cy_i28->d1 = _cy_i42[ithread];
			cy_i30->d0 = _cy_r43[ithread];	cy_i30->d1 = _cy_i43[ithread];
			cy_i32->d0 = _cy_r44[ithread];	cy_i32->d1 = _cy_i44[ithread];
			cy_i34->d0 = _cy_r45[ithread];	cy_i34->d1 = _cy_i45[ithread];
			cy_i36->d0 = _cy_r46[ithread];	cy_i36->d1 = _cy_i46[ithread];
			cy_i38->d0 = _cy_r47[ithread];	cy_i38->d1 = _cy_i47[ithread];
			cy_i40->d0 = _cy_r48[ithread];	cy_i40->d1 = _cy_i48[ithread];
			cy_i42->d0 = _cy_r49[ithread];	cy_i42->d1 = _cy_i49[ithread];
			cy_i44->d0 = _cy_r50[ithread];	cy_i44->d1 = _cy_i50[ithread];
			cy_i46->d0 = _cy_r51[ithread];	cy_i46->d1 = _cy_i51[ithread];
			cy_i48->d0 = _cy_r52[ithread];	cy_i48->d1 = _cy_i52[ithread];
			cy_i50->d0 = _cy_r53[ithread];	cy_i50->d1 = _cy_i53[ithread];
			cy_i52->d0 = _cy_r54[ithread];	cy_i52->d1 = _cy_i54[ithread];
			cy_i54->d0 = _cy_r55[ithread];	cy_i54->d1 = _cy_i55[ithread];
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
			cy_r28 = _cy_r28[ithread];		cy_i28 = _cy_i28[ithread];
			cy_r29 = _cy_r29[ithread];		cy_i29 = _cy_i29[ithread];
			cy_r30 = _cy_r30[ithread];		cy_i30 = _cy_i30[ithread];
			cy_r31 = _cy_r31[ithread];		cy_i31 = _cy_i31[ithread];
			cy_r32 = _cy_r32[ithread];		cy_i32 = _cy_i32[ithread];
			cy_r33 = _cy_r33[ithread];		cy_i33 = _cy_i33[ithread];
			cy_r34 = _cy_r34[ithread];		cy_i34 = _cy_i34[ithread];
			cy_r35 = _cy_r35[ithread];		cy_i35 = _cy_i35[ithread];
			cy_r36 = _cy_r36[ithread];		cy_i36 = _cy_i36[ithread];
			cy_r37 = _cy_r37[ithread];		cy_i37 = _cy_i37[ithread];
			cy_r38 = _cy_r38[ithread];		cy_i38 = _cy_i38[ithread];
			cy_r39 = _cy_r39[ithread];		cy_i39 = _cy_i39[ithread];
			cy_r40 = _cy_r40[ithread];		cy_i40 = _cy_i40[ithread];
			cy_r41 = _cy_r41[ithread];		cy_i41 = _cy_i41[ithread];
			cy_r42 = _cy_r42[ithread];		cy_i42 = _cy_i42[ithread];
			cy_r43 = _cy_r43[ithread];		cy_i43 = _cy_i43[ithread];
			cy_r44 = _cy_r44[ithread];		cy_i44 = _cy_i44[ithread];
			cy_r45 = _cy_r45[ithread];		cy_i45 = _cy_i45[ithread];
			cy_r46 = _cy_r46[ithread];		cy_i46 = _cy_i46[ithread];
			cy_r47 = _cy_r47[ithread];		cy_i47 = _cy_i47[ithread];
			cy_r48 = _cy_r48[ithread];		cy_i48 = _cy_i48[ithread];
			cy_r49 = _cy_r49[ithread];		cy_i49 = _cy_i49[ithread];
			cy_r50 = _cy_r50[ithread];		cy_i50 = _cy_i50[ithread];
			cy_r51 = _cy_r51[ithread];		cy_i51 = _cy_i51[ithread];
			cy_r52 = _cy_r52[ithread];		cy_i52 = _cy_i52[ithread];
			cy_r53 = _cy_r53[ithread];		cy_i53 = _cy_i53[ithread];
			cy_r54 = _cy_r54[ithread];		cy_i54 = _cy_i54[ithread];
			cy_r55 = _cy_r55[ithread];		cy_i55 = _cy_i55[ithread];
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

			/*...The radix-56 DIT pass is here:	*/

			#ifdef USE_SSE2

				/* Outputs of radix-4 are in consecutive memory locs, i.e. 0x20 bytes apart, e.g. the distance between r00r and r01r: */

				add0 = &a[j1    ];
				add1 = add0+p01;
				add2 = add0+p03;
				add3 = add0+p02;
				add4 = add0+p07;
				add5 = add0+p06;
				add6 = add0+p05;
				add7 = add0+p04;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r00r, isrt2)

				add3 = &a[j1+p30];
				add0 = add3+p03;
				add1 = add3+p02;
				add2 = add3+p01;
				add4 = add3+p05;
				add5 = add3+p04;
				add6 = add3+p06;
				add7 = add3+p07;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r08r, isrt2)

				add5 = &a[j1+p28];
				add0 = add5+p05;
				add1 = add5+p04;
				add2 = add5+p06;
				add3 = add5+p07;
				add4 = add5+p01;
				add6 = add5+p02;
				add7 = add5+p03;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r16r, isrt2)

				add1 = &a[j1+p20];
				add0 = add1+p01;
				add2 = add1+p02;
				add3 = add1+p03;
				add4 = add1+p06;
				add5 = add1+p07;
				add6 = add1+p04;
				add7 = add1+p05;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r24r, isrt2)

				add6 = &a[j1+p18];
				add0 = add6+p06;
				add1 = add6+p07;
				add2 = add6+p04;
				add3 = add6+p05;
				add4 = add6+p02;
				add5 = add6+p03;
				add7 = add6+p01;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r32r, isrt2)

				add2 = &a[j1+p10];
				add0 = add2+p02;
				add1 = add2+p03;
				add3 = add2+p01;
				add4 = add2+p04;
				add5 = add2+p05;
				add6 = add2+p07;
				add7 = add2+p06;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r40r, isrt2)

				add4 = &a[j1+p08];
				add0 = add4+p04;
				add1 = add4+p05;
				add2 = add4+p07;
				add3 = add4+p06;
				add5 = add4+p01;
				add6 = add4+p03;
				add7 = add4+p02;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r48r, isrt2)

			/*...and now do 8 radix-7 transforms, with the columns of r*[r,i] output pairs in the above 7x radix-8 set now acting as input rows: */
								/*            inputs        */  /* sincos ptr */ /*            outputs                   */
				SSE2_RADIX_07_DFT(r00r,r08r,r16r,r24r,r32r,r40r,r48r, cc0, s1p00r,s1p08r,s1p16r,s1p24r,s1p32r,s1p40r,s1p48r);
				SSE2_RADIX_07_DFT(r01r,r09r,r17r,r25r,r33r,r41r,r49r, cc0, s1p49r,s1p01r,s1p09r,s1p17r,s1p25r,s1p33r,s1p41r);
				SSE2_RADIX_07_DFT(r02r,r10r,r18r,r26r,r34r,r42r,r50r, cc0, s1p42r,s1p50r,s1p02r,s1p10r,s1p18r,s1p26r,s1p34r);
				SSE2_RADIX_07_DFT(r03r,r11r,r19r,r27r,r35r,r43r,r51r, cc0, s1p35r,s1p43r,s1p51r,s1p03r,s1p11r,s1p19r,s1p27r);
				SSE2_RADIX_07_DFT(r04r,r12r,r20r,r28r,r36r,r44r,r52r, cc0, s1p28r,s1p36r,s1p44r,s1p52r,s1p04r,s1p12r,s1p20r);
				SSE2_RADIX_07_DFT(r05r,r13r,r21r,r29r,r37r,r45r,r53r, cc0, s1p21r,s1p29r,s1p37r,s1p45r,s1p53r,s1p05r,s1p13r);
				SSE2_RADIX_07_DFT(r06r,r14r,r22r,r30r,r38r,r46r,r54r, cc0, s1p14r,s1p22r,s1p30r,s1p38r,s1p46r,s1p54r,s1p06r);
				SSE2_RADIX_07_DFT(r07r,r15r,r23r,r31r,r39r,r47r,r55r, cc0, s1p07r,s1p15r,s1p23r,s1p31r,s1p39r,s1p47r,s1p55r);

			#else	// USE_SSE2 = False:

			/*...gather the needed data (56 64-bit complex, i.e. 112 64-bit reals) and do 7 radix-8 transforms...*/
							/*                                    inputs                                  */ /*                         outputs                   */
				RADIX_08_DIT(a[j1    ],a[j2    ], a[j1+p01],a[j2+p01], a[j1+p03],a[j2+p03], a[j1+p02],a[j2+p02], a[j1+p07],a[j2+p07], a[j1+p06],a[j2+p06], a[j1+p05],a[j2+p05], a[j1+p04],a[j2+p04], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i, rt,it);	jt = j1+p30; jp = j2+p30;
				RADIX_08_DIT(a[jt+p03],a[jp+p03], a[jt+p02],a[jp+p02], a[jt+p01],a[jp+p01], a[jt    ],a[jp    ], a[jt+p05],a[jp+p05], a[jt+p04],a[jp+p04], a[jt+p06],a[jp+p06], a[jt+p07],a[jp+p07], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r08r,r08i,r09r,r09i,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i, rt,it);	jt = j1+p28; jp = j2+p28;
				RADIX_08_DIT(a[jt+p05],a[jp+p05], a[jt+p04],a[jp+p04], a[jt+p06],a[jp+p06], a[jt+p07],a[jp+p07], a[jt+p01],a[jp+p01], a[jt    ],a[jp    ], a[jt+p02],a[jp+p02], a[jt+p03],a[jp+p03], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i, rt,it);	jt = j1+p20; jp = j2+p20;
				RADIX_08_DIT(a[jt+p01],a[jp+p01], a[jt    ],a[jp    ], a[jt+p02],a[jp+p02], a[jt+p03],a[jp+p03], a[jt+p06],a[jp+p06], a[jt+p07],a[jp+p07], a[jt+p04],a[jp+p04], a[jt+p05],a[jp+p05], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i,r30r,r30i,r31r,r31i, rt,it);	jt = j1+p18; jp = j2+p18;
				RADIX_08_DIT(a[jt+p06],a[jp+p06], a[jt+p07],a[jp+p07], a[jt+p04],a[jp+p04], a[jt+p05],a[jp+p05], a[jt+p02],a[jp+p02], a[jt+p03],a[jp+p03], a[jt    ],a[jp    ], a[jt+p01],a[jp+p01], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i, rt,it);	jt = j1+p10; jp = j2+p10;
				RADIX_08_DIT(a[jt+p02],a[jp+p02], a[jt+p03],a[jp+p03], a[jt    ],a[jp    ], a[jt+p01],a[jp+p01], a[jt+p04],a[jp+p04], a[jt+p05],a[jp+p05], a[jt+p07],a[jp+p07], a[jt+p06],a[jp+p06], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i, rt,it);	jt = j1+p08; jp = j2+p08;
				RADIX_08_DIT(a[jt+p04],a[jp+p04], a[jt+p05],a[jp+p05], a[jt+p07],a[jp+p07], a[jt+p06],a[jp+p06], a[jt    ],a[jp    ], a[jt+p01],a[jp+p01], a[jt+p03],a[jp+p03], a[jt+p02],a[jp+p02], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r48r,r48i,r49r,r49i,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i, rt,it);

			/*...and now do 8 radix-7 transforms, with the columns of r*[r,i] output pairs in the above 7x radix-8 set now acting as input rows: */
							/*                            inputs                                */  /*                    intermediates                  */  /*                                                                     outputs                       */
				RADIX_07_DFT(r00r,r00i, r08r,r08i, r16r,r16i, r24r,r24i, r32r,r32i, r40r,r40i, r48r,r48i,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p00r,a1p00i, a1p08r,a1p08i, a1p16r,a1p16i, a1p24r,a1p24i, a1p32r,a1p32i, a1p40r,a1p40i, a1p48r,a1p48i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(r01r,r01i, r09r,r09i, r17r,r17i, r25r,r25i, r33r,r33i, r41r,r41i, r49r,r49i,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p49r,a1p49i, a1p01r,a1p01i, a1p09r,a1p09i, a1p17r,a1p17i, a1p25r,a1p25i, a1p33r,a1p33i, a1p41r,a1p41i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(r02r,r02i, r10r,r10i, r18r,r18i, r26r,r26i, r34r,r34i, r42r,r42i, r50r,r50i,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p42r,a1p42i, a1p50r,a1p50i, a1p02r,a1p02i, a1p10r,a1p10i, a1p18r,a1p18i, a1p26r,a1p26i, a1p34r,a1p34i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(r03r,r03i, r11r,r11i, r19r,r19i, r27r,r27i, r35r,r35i, r43r,r43i, r51r,r51i,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p35r,a1p35i, a1p43r,a1p43i, a1p51r,a1p51i, a1p03r,a1p03i, a1p11r,a1p11i, a1p19r,a1p19i, a1p27r,a1p27i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(r04r,r04i, r12r,r12i, r20r,r20i, r28r,r28i, r36r,r36i, r44r,r44i, r52r,r52i,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p28r,a1p28i, a1p36r,a1p36i, a1p44r,a1p44i, a1p52r,a1p52i, a1p04r,a1p04i, a1p12r,a1p12i, a1p20r,a1p20i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(r05r,r05i, r13r,r13i, r21r,r21i, r29r,r29i, r37r,r37i, r45r,r45i, r53r,r53i,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p21r,a1p21i, a1p29r,a1p29i, a1p37r,a1p37i, a1p45r,a1p45i, a1p53r,a1p53i, a1p05r,a1p05i, a1p13r,a1p13i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(r06r,r06i, r14r,r14i, r22r,r22i, r30r,r30i, r38r,r38i, r46r,r46i, r54r,r54i,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p14r,a1p14i, a1p22r,a1p22i, a1p30r,a1p30i, a1p38r,a1p38i, a1p46r,a1p46i, a1p54r,a1p54i, a1p06r,a1p06i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(r07r,r07i, r15r,r15i, r23r,r23i, r31r,r31i, r39r,r39i, r47r,r47i, r55r,r55i,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p07r,a1p07i, a1p15r,a1p15i, a1p23r,a1p23i, a1p31r,a1p31i, a1p39r,a1p39i, a1p47r,a1p47i, a1p55r,a1p55i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);

			#endif	/* USE_SSE2 */

			/*...Now do the carries. Since the outputs would
			normally be getting dispatched to 28 separate blocks of the A-array, we need 28 separate carries.	*/

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

				AVX_cmplx_carry_norm_errcheck0_X4(s1p00r,add1,add2,add3,cy_r00,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p04r,add1,add2,add3,cy_r04,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p08r,add1,add2,add3,cy_r08,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p12r,add1,add2,add3,cy_r12,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p16r,add1,add2,add3,cy_r16,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p20r,add1,add2,add3,cy_r20,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p24r,add1,add2,add3,cy_r24,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p28r,add1,add2,add3,cy_r28,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p32r,add1,add2,add3,cy_r32,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p36r,add1,add2,add3,cy_r36,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p40r,add1,add2,add3,cy_r40,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p44r,add1,add2,add3,cy_r44,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p48r,add1,add2,add3,cy_r48,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p52r,add1,add2,add3,cy_r52,bjmodn52,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);

				co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

			#elif defined(USE_SSE2)

				l= j & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
				ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
				ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
				ctmp->re = wtnm1;	ctmp->im = wtnm1;

				add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

			  #if defined(COMPILER_TYPE_MSVC)

			   #ifdef ERR_CHECK_ALL				/* Updating i prior to the 2nd-7th macro calls allows use of the same 0_2B macro for all */
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00);	// i =((uint32)(sw - *bjmodn04) >> 31);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04);	// i =((uint32)(sw - *bjmodn08) >> 31);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08);	// i =((uint32)(sw - *bjmodn12) >> 31);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12);	// i =((uint32)(sw - *bjmodn16) >> 31);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16);	// i =((uint32)(sw - *bjmodn20) >> 31);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20);	// i =((uint32)(sw - *bjmodn24) >> 31);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24);
			   #else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24);
			   #endif

			  #else	/* GCC-style inline ASM: */

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p28r,add1,add2,add3,cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p32r,add1,add2,add3,cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p36r,add1,add2,add3,cy_r36,cy_r38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p40r,add1,add2,add3,cy_r40,cy_r42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p44r,add1,add2,add3,cy_r44,cy_r46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p48r,add1,add2,add3,cy_r48,cy_r50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p52r,add1,add2,add3,cy_r52,cy_r54,bjmodn52,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			   #else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p28r,add1,add2,add3,cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p32r,add1,add2,add3,cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p36r,add1,add2,add3,cy_r36,cy_r38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p40r,add1,add2,add3,cy_r40,cy_r42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p44r,add1,add2,add3,cy_r44,cy_r46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p48r,add1,add2,add3,cy_r48,cy_r50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p52r,add1,add2,add3,cy_r52,cy_r54,bjmodn52,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			   #endif

				/* Bizarre - when I disabled the diagnostic prints above and below, the resulting GCC build immediately gave
					fatal roundoff errors starting on iteration #5 - so insert the bogus [never taken] if() here as a workaround.
					Equally bizarre, inserting the bogus if() *before* the 4 carry-macro calls above gave the correct result as well,
					but ran fully 10% slower. Good old GCC...
				Dec 2011: Suspect this was a side effect of my gcc asm macros not including cc/memory in the clobber list, because
				the code now runs correctly without this hack ... but the code runs sign. faster with iy left in. So still "bizarre" but in a new way.
				*/
				if(j < 0)
				{
					fprintf(stderr, "Iter %3d\n",iter);
				}

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
				ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
				ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
				ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
				ctmp->re = wtnm1;	ctmp->im = wtnm1;

			/*	i =((uint32)(sw - *bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];

			  #if defined(COMPILER_TYPE_MSVC)

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24);
			   #else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24);
			   #endif

			  #else	/* GCC-style inline ASM: */

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p28r,add1,add2,     cy_r28,cy_r30,bjmodn28,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p32r,add1,add2,     cy_r32,cy_r34,bjmodn32,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p36r,add1,add2,     cy_r36,cy_r38,bjmodn36,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p40r,add1,add2,     cy_r40,cy_r42,bjmodn40,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p44r,add1,add2,     cy_r44,cy_r46,bjmodn44,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p48r,add1,add2,     cy_r48,cy_r50,bjmodn48,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p52r,add1,add2,     cy_r52,cy_r54,bjmodn52,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			   #else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p28r,add1,add2,     cy_r28,cy_r30,bjmodn28,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p32r,add1,add2,     cy_r32,cy_r34,bjmodn32,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p36r,add1,add2,     cy_r36,cy_r38,bjmodn36,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p40r,add1,add2,     cy_r40,cy_r42,bjmodn40,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p44r,add1,add2,     cy_r44,cy_r46,bjmodn44,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p48r,add1,add2,     cy_r48,cy_r50,bjmodn48,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p52r,add1,add2,     cy_r52,cy_r54,bjmodn52,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			   #endif

			  #endif

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

			#else	// Scalar-double mode:

				l= j & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

			/*...set0 is slightly different from others:	*/
			   cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy_r00,bjmodn00   );
				cmplx_carry_norm_errcheck(a1p01r,a1p01i,cy_r01,bjmodn01,1 );
				cmplx_carry_norm_errcheck(a1p02r,a1p02i,cy_r02,bjmodn02,2 );
				cmplx_carry_norm_errcheck(a1p03r,a1p03i,cy_r03,bjmodn03,3 );
				cmplx_carry_norm_errcheck(a1p04r,a1p04i,cy_r04,bjmodn04,4 );
				cmplx_carry_norm_errcheck(a1p05r,a1p05i,cy_r05,bjmodn05,5 );
				cmplx_carry_norm_errcheck(a1p06r,a1p06i,cy_r06,bjmodn06,6 );
				cmplx_carry_norm_errcheck(a1p07r,a1p07i,cy_r07,bjmodn07,7 );
				cmplx_carry_norm_errcheck(a1p08r,a1p08i,cy_r08,bjmodn08,8 );
				cmplx_carry_norm_errcheck(a1p09r,a1p09i,cy_r09,bjmodn09,9 );
				cmplx_carry_norm_errcheck(a1p10r,a1p10i,cy_r10,bjmodn10,10);
				cmplx_carry_norm_errcheck(a1p11r,a1p11i,cy_r11,bjmodn11,11);
				cmplx_carry_norm_errcheck(a1p12r,a1p12i,cy_r12,bjmodn12,12);
				cmplx_carry_norm_errcheck(a1p13r,a1p13i,cy_r13,bjmodn13,13);
				cmplx_carry_norm_errcheck(a1p14r,a1p14i,cy_r14,bjmodn14,14);
				cmplx_carry_norm_errcheck(a1p15r,a1p15i,cy_r15,bjmodn15,15);
				cmplx_carry_norm_errcheck(a1p16r,a1p16i,cy_r16,bjmodn16,16);
				cmplx_carry_norm_errcheck(a1p17r,a1p17i,cy_r17,bjmodn17,17);
				cmplx_carry_norm_errcheck(a1p18r,a1p18i,cy_r18,bjmodn18,18);
				cmplx_carry_norm_errcheck(a1p19r,a1p19i,cy_r19,bjmodn19,19);
				cmplx_carry_norm_errcheck(a1p20r,a1p20i,cy_r20,bjmodn20,20);
				cmplx_carry_norm_errcheck(a1p21r,a1p21i,cy_r21,bjmodn21,21);
				cmplx_carry_norm_errcheck(a1p22r,a1p22i,cy_r22,bjmodn22,22);
				cmplx_carry_norm_errcheck(a1p23r,a1p23i,cy_r23,bjmodn23,23);
				cmplx_carry_norm_errcheck(a1p24r,a1p24i,cy_r24,bjmodn24,24);
				cmplx_carry_norm_errcheck(a1p25r,a1p25i,cy_r25,bjmodn25,25);
				cmplx_carry_norm_errcheck(a1p26r,a1p26i,cy_r26,bjmodn26,26);
				cmplx_carry_norm_errcheck(a1p27r,a1p27i,cy_r27,bjmodn27,27);
			    cmplx_carry_norm_errcheck(a1p28r,a1p28i,cy_r28,bjmodn28,28);
				cmplx_carry_norm_errcheck(a1p29r,a1p29i,cy_r29,bjmodn29,29);
				cmplx_carry_norm_errcheck(a1p30r,a1p30i,cy_r30,bjmodn30,30);
				cmplx_carry_norm_errcheck(a1p31r,a1p31i,cy_r31,bjmodn31,31);
				cmplx_carry_norm_errcheck(a1p32r,a1p32i,cy_r32,bjmodn32,32);
				cmplx_carry_norm_errcheck(a1p33r,a1p33i,cy_r33,bjmodn33,33);
				cmplx_carry_norm_errcheck(a1p34r,a1p34i,cy_r34,bjmodn34,34);
				cmplx_carry_norm_errcheck(a1p35r,a1p35i,cy_r35,bjmodn35,35);
				cmplx_carry_norm_errcheck(a1p36r,a1p36i,cy_r36,bjmodn36,36);
				cmplx_carry_norm_errcheck(a1p37r,a1p37i,cy_r37,bjmodn37,37);
				cmplx_carry_norm_errcheck(a1p38r,a1p38i,cy_r38,bjmodn38,38);
				cmplx_carry_norm_errcheck(a1p39r,a1p39i,cy_r39,bjmodn39,39);
				cmplx_carry_norm_errcheck(a1p40r,a1p40i,cy_r40,bjmodn40,40);
				cmplx_carry_norm_errcheck(a1p41r,a1p41i,cy_r41,bjmodn41,41);
				cmplx_carry_norm_errcheck(a1p42r,a1p42i,cy_r42,bjmodn42,42);
				cmplx_carry_norm_errcheck(a1p43r,a1p43i,cy_r43,bjmodn43,43);
				cmplx_carry_norm_errcheck(a1p44r,a1p44i,cy_r44,bjmodn44,44);
				cmplx_carry_norm_errcheck(a1p45r,a1p45i,cy_r45,bjmodn45,45);
				cmplx_carry_norm_errcheck(a1p46r,a1p46i,cy_r46,bjmodn46,46);
				cmplx_carry_norm_errcheck(a1p47r,a1p47i,cy_r47,bjmodn47,47);
				cmplx_carry_norm_errcheck(a1p48r,a1p48i,cy_r48,bjmodn48,48);
				cmplx_carry_norm_errcheck(a1p49r,a1p49i,cy_r49,bjmodn49,49);
				cmplx_carry_norm_errcheck(a1p50r,a1p50i,cy_r50,bjmodn50,50);
				cmplx_carry_norm_errcheck(a1p51r,a1p51i,cy_r51,bjmodn51,51);
				cmplx_carry_norm_errcheck(a1p52r,a1p52i,cy_r52,bjmodn52,52);
				cmplx_carry_norm_errcheck(a1p53r,a1p53i,cy_r53,bjmodn53,53);
				cmplx_carry_norm_errcheck(a1p54r,a1p54i,cy_r54,bjmodn54,54);
				cmplx_carry_norm_errcheck(a1p55r,a1p55i,cy_r55,bjmodn55,55);

				i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

			#endif	// USE_AVX?

			}
			else	/* MODULUS_TYPE_FERMAT */
			{

			#ifdef USE_AVX

				// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.

				/* Get the needed Nth root of -1: */
				add1 = (double *)&rn0[0];
				add2 = (double *)&rn1[0];

				idx_offset = j;
				idx_incr = NDIVR;

				tmp = base_negacyclic_root;	tm2 = tmp+1;

			  #if HIACC
				// Hi-accuracy version needs 7 copies of each base root:
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

			  #else	// HIACC = false:

				// Get the needed quartet of Nth roots of -1: This is the same code as in the scalar
				// fermat_carry_norm_errcheck() macro, with the single index j replaced by the quartet j,j+2,j+4,j+6:
				l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

				// The above need some inits to prepare for the AVX version of the Fermat-mod carry macro:
				SSE2_fermat_carry_init_loacc(base_negacyclic_root);

			  #endif

				// AVX-custom 4-way carry macro - each contains 4 of the RADIX stride-n/RADIX-separated carries
				// (processed independently in parallel), and steps through sequential-data indicies j,j+2,j+4,j+6:
			  #if HIACC
				// The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl)
				// which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
				// precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
				// uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls:
																																					// *cycle0 index increments by +4 (mod odd_radix) between macro calls
				tmp = base_negacyclic_root+  0;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p00r,tmp,0xe00,cy_r00,cy_i00,odd_radix,half_arr,sign_mask,icycle0,icycle1,icycle2,icycle3, jcycle0,kcycle0,lcycle0);
				tmp = base_negacyclic_root+  8;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p04r,tmp,0xd40,cy_r04,cy_i04,odd_radix,half_arr,sign_mask,icycle4,icycle5,icycle6,icycle0, jcycle4,kcycle4,lcycle4);
				tmp = base_negacyclic_root+ 16;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p08r,tmp,0xc80,cy_r08,cy_i08,odd_radix,half_arr,sign_mask,icycle1,icycle2,icycle3,icycle4, jcycle1,kcycle1,lcycle1);
				tmp = base_negacyclic_root+ 24;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p12r,tmp,0xbc0,cy_r12,cy_i12,odd_radix,half_arr,sign_mask,icycle5,icycle6,icycle0,icycle1, jcycle5,kcycle5,lcycle5);
				tmp = base_negacyclic_root+ 32;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p16r,tmp,0xb00,cy_r16,cy_i16,odd_radix,half_arr,sign_mask,icycle2,icycle3,icycle4,icycle5, jcycle2,kcycle2,lcycle2);
				tmp = base_negacyclic_root+ 40;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p20r,tmp,0xa40,cy_r20,cy_i20,odd_radix,half_arr,sign_mask,icycle6,icycle0,icycle1,icycle2, jcycle6,kcycle6,lcycle6);
				tmp = base_negacyclic_root+ 48;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p24r,tmp,0x980,cy_r24,cy_i24,odd_radix,half_arr,sign_mask,icycle3,icycle4,icycle5,icycle6, jcycle3,kcycle3,lcycle3);
				tmp = base_negacyclic_root+ 56;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p28r,tmp,0x8c0,cy_r28,cy_i28,odd_radix,half_arr,sign_mask,icycle0,icycle1,icycle2,icycle3, jcycle0,kcycle0,lcycle0);
				tmp = base_negacyclic_root+ 64;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p32r,tmp,0x800,cy_r32,cy_i32,odd_radix,half_arr,sign_mask,icycle4,icycle5,icycle6,icycle0, jcycle4,kcycle4,lcycle4);
				tmp = base_negacyclic_root+ 72;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p36r,tmp,0x740,cy_r36,cy_i36,odd_radix,half_arr,sign_mask,icycle1,icycle2,icycle3,icycle4, jcycle1,kcycle1,lcycle1);
				tmp = base_negacyclic_root+ 80;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p40r,tmp,0x680,cy_r40,cy_i40,odd_radix,half_arr,sign_mask,icycle5,icycle6,icycle0,icycle1, jcycle5,kcycle5,lcycle5);
				tmp = base_negacyclic_root+ 88;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p44r,tmp,0x5c0,cy_r44,cy_i44,odd_radix,half_arr,sign_mask,icycle2,icycle3,icycle4,icycle5, jcycle2,kcycle2,lcycle2);
				tmp = base_negacyclic_root+ 96;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p48r,tmp,0x500,cy_r48,cy_i48,odd_radix,half_arr,sign_mask,icycle6,icycle0,icycle1,icycle2, jcycle6,kcycle6,lcycle6);
				tmp = base_negacyclic_root+104;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p52r,tmp,0x440,cy_r52,cy_i52,odd_radix,half_arr,sign_mask,icycle3,icycle4,icycle5,icycle6, jcycle3,kcycle3,lcycle3);

			  #else	// HIACC = false:
																																					// *cycle0 index increments by +4 (mod odd_radix) between macro calls
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p00r,base_negacyclic_root,cy_r00,cy_i00,odd_radix,half_arr,sign_mask,                     icycle0,icycle1,icycle2,icycle3, jcycle0,kcycle0,lcycle0);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p04r,base_negacyclic_root,cy_r04,cy_i04,odd_radix,half_arr,sign_mask,                     icycle4,icycle5,icycle6,icycle0, jcycle4,kcycle4,lcycle4);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p08r,base_negacyclic_root,cy_r08,cy_i08,odd_radix,half_arr,sign_mask,                     icycle1,icycle2,icycle3,icycle4, jcycle1,kcycle1,lcycle1);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p12r,base_negacyclic_root,cy_r12,cy_i12,odd_radix,half_arr,sign_mask,                     icycle5,icycle6,icycle0,icycle1, jcycle5,kcycle5,lcycle5);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p16r,base_negacyclic_root,cy_r16,cy_i16,odd_radix,half_arr,sign_mask,                     icycle2,icycle3,icycle4,icycle5, jcycle2,kcycle2,lcycle2);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p20r,base_negacyclic_root,cy_r20,cy_i20,odd_radix,half_arr,sign_mask,                     icycle6,icycle0,icycle1,icycle2, jcycle6,kcycle6,lcycle6);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p24r,base_negacyclic_root,cy_r24,cy_i24,odd_radix,half_arr,sign_mask,                     icycle3,icycle4,icycle5,icycle6, jcycle3,kcycle3,lcycle3);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p28r,base_negacyclic_root,cy_r28,cy_i28,odd_radix,half_arr,sign_mask,                     icycle0,icycle1,icycle2,icycle3, jcycle0,kcycle0,lcycle0);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p32r,base_negacyclic_root,cy_r32,cy_i32,odd_radix,half_arr,sign_mask,                     icycle4,icycle5,icycle6,icycle0, jcycle4,kcycle4,lcycle4);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p36r,base_negacyclic_root,cy_r36,cy_i36,odd_radix,half_arr,sign_mask,                     icycle1,icycle2,icycle3,icycle4, jcycle1,kcycle1,lcycle1);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p40r,base_negacyclic_root,cy_r40,cy_i40,odd_radix,half_arr,sign_mask,                     icycle5,icycle6,icycle0,icycle1, jcycle5,kcycle5,lcycle5);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p44r,base_negacyclic_root,cy_r44,cy_i44,odd_radix,half_arr,sign_mask,                     icycle2,icycle3,icycle4,icycle5, jcycle2,kcycle2,lcycle2);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p48r,base_negacyclic_root,cy_r48,cy_i48,odd_radix,half_arr,sign_mask,                     icycle6,icycle0,icycle1,icycle2, jcycle6,kcycle6,lcycle6);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p52r,base_negacyclic_root,cy_r52,cy_i52,odd_radix,half_arr,sign_mask,                     icycle3,icycle4,icycle5,icycle6, jcycle3,kcycle3,lcycle3);

			  #endif	// HIACC?

			#elif defined(USE_SSE2)

				/* Get the needed Nth root of -1: */
				add1 = (double *)&rn0[0];
				add2 = (double *)&rn1[0];

				idx_offset = j;
				idx_incr = NDIVR;

			  #if defined(COMPILER_TYPE_MSVC)

			  /* The cy_[r|i]_idx[A|B] names here are not meaningful, each simply stores one [re,im] carry pair,
			  e.g. cy_r01 stores the carries our of [a0.re,a0.im], cy_r23 stores the carries our of [a1.re,a1.im], etc.
			  Here is the actual mapping between these SSE2-mode 2-vector carry pairs and the scalar carries:
																	2-vector				                                          Scalar
																	--------	 		                                           ------------- */
				SSE2_fermat_carry_norm_errcheck(s1p00r,cy_r00,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r00,cy_i00 */
				SSE2_fermat_carry_norm_errcheck(s1p01r,cy_r02,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r01,cy_i01 */
				SSE2_fermat_carry_norm_errcheck(s1p02r,cy_r04,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r02,cy_i02 */
				SSE2_fermat_carry_norm_errcheck(s1p03r,cy_r06,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r03,cy_i03 */
				SSE2_fermat_carry_norm_errcheck(s1p04r,cy_r08,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r04,cy_i04 */
				SSE2_fermat_carry_norm_errcheck(s1p05r,cy_r10,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r05,cy_i05 */
				SSE2_fermat_carry_norm_errcheck(s1p06r,cy_r12,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r06,cy_i06 */
				SSE2_fermat_carry_norm_errcheck(s1p07r,cy_r14,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r07,cy_i07 */
				SSE2_fermat_carry_norm_errcheck(s1p08r,cy_r16,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r08,cy_i08 */
				SSE2_fermat_carry_norm_errcheck(s1p09r,cy_r18,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r09,cy_i09 */
				SSE2_fermat_carry_norm_errcheck(s1p10r,cy_r20,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r10,cy_i10 */
				SSE2_fermat_carry_norm_errcheck(s1p11r,cy_r22,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r11,cy_i11 */
				SSE2_fermat_carry_norm_errcheck(s1p12r,cy_r24,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r12,cy_i12 */
				SSE2_fermat_carry_norm_errcheck(s1p13r,cy_r26,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r13,cy_i13 */
				SSE2_fermat_carry_norm_errcheck(s1p14r,cy_r28,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r14,cy_i14 */
				SSE2_fermat_carry_norm_errcheck(s1p15r,cy_r30,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r15,cy_i15 */
				SSE2_fermat_carry_norm_errcheck(s1p16r,cy_r32,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r16,cy_i16 */
				SSE2_fermat_carry_norm_errcheck(s1p17r,cy_r34,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r17,cy_i17 */
				SSE2_fermat_carry_norm_errcheck(s1p18r,cy_r36,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r18,cy_i18 */
				SSE2_fermat_carry_norm_errcheck(s1p19r,cy_r38,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r19,cy_i19 */
				SSE2_fermat_carry_norm_errcheck(s1p20r,cy_r40,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r20,cy_i20 */
				SSE2_fermat_carry_norm_errcheck(s1p21r,cy_r42,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r21,cy_i21 */
				SSE2_fermat_carry_norm_errcheck(s1p22r,cy_r44,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r22,cy_i22 */
				SSE2_fermat_carry_norm_errcheck(s1p23r,cy_r46,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r23,cy_i23 */
				SSE2_fermat_carry_norm_errcheck(s1p24r,cy_r48,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r24,cy_i24 */
				SSE2_fermat_carry_norm_errcheck(s1p25r,cy_r50,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r25,cy_i25 */
				SSE2_fermat_carry_norm_errcheck(s1p26r,cy_r52,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r26,cy_i26 */
				SSE2_fermat_carry_norm_errcheck(s1p27r,cy_r54,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r27,cy_i27 */

				SSE2_fermat_carry_norm_errcheck(s1pr28,cy_i00,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r28,cy_i28 */
				SSE2_fermat_carry_norm_errcheck(s1pr29,cy_i02,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r29,cy_i29 */
				SSE2_fermat_carry_norm_errcheck(s1pr30,cy_i04,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r30,cy_i30 */
				SSE2_fermat_carry_norm_errcheck(s1pr31,cy_i06,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r31,cy_i31 */
				SSE2_fermat_carry_norm_errcheck(s1pr32,cy_i08,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r32,cy_i32 */
				SSE2_fermat_carry_norm_errcheck(s1pr33,cy_i10,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r33,cy_i33 */
				SSE2_fermat_carry_norm_errcheck(s1pr34,cy_i12,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r34,cy_i34 */
				SSE2_fermat_carry_norm_errcheck(s1pr35,cy_i14,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r35,cy_i35 */
				SSE2_fermat_carry_norm_errcheck(s1pr36,cy_i16,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r36,cy_i36 */
				SSE2_fermat_carry_norm_errcheck(s1pr37,cy_i18,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r37,cy_i37 */
				SSE2_fermat_carry_norm_errcheck(s1pr38,cy_i20,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r38,cy_i38 */
				SSE2_fermat_carry_norm_errcheck(s1pr39,cy_i22,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r39,cy_i39 */
				SSE2_fermat_carry_norm_errcheck(s1pr40,cy_i24,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r40,cy_i40 */
				SSE2_fermat_carry_norm_errcheck(s1pr41,cy_i26,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r41,cy_i41 */
				SSE2_fermat_carry_norm_errcheck(s1pr42,cy_i28,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r42,cy_i42 */
				SSE2_fermat_carry_norm_errcheck(s1pr43,cy_i30,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r43,cy_i43 */
				SSE2_fermat_carry_norm_errcheck(s1pr44,cy_i32,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r44,cy_i44 */
				SSE2_fermat_carry_norm_errcheck(s1pr45,cy_i34,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r45,cy_i45 */
				SSE2_fermat_carry_norm_errcheck(s1pr46,cy_i36,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r46,cy_i46 */
				SSE2_fermat_carry_norm_errcheck(s1pr47,cy_i38,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r47,cy_i47 */
				SSE2_fermat_carry_norm_errcheck(s1pr48,cy_i40,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r48,cy_i48 */
				SSE2_fermat_carry_norm_errcheck(s1pr49,cy_i42,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r49,cy_i49 */
				SSE2_fermat_carry_norm_errcheck(s1pr50,cy_i44,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r50,cy_i50 */
				SSE2_fermat_carry_norm_errcheck(s1pr51,cy_i46,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r51,cy_i51 */
				SSE2_fermat_carry_norm_errcheck(s1pr52,cy_i48,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r52,cy_i52 */
				SSE2_fermat_carry_norm_errcheck(s1pr53,cy_i50,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r53,cy_i53 */
				SSE2_fermat_carry_norm_errcheck(s1pr54,cy_i52,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r54,cy_i54 */
				SSE2_fermat_carry_norm_errcheck(s1pr55,cy_i54,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r55,cy_i55 */

			  #elif (OS_BITS == 32)

				SSE2_fermat_carry_norm_errcheck(s1p00r,cy_r00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck(s1p01r,cy_r02,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck(s1p02r,cy_r04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck(s1p03r,cy_r06,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck(s1p04r,cy_r08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck(s1p05r,cy_r10,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck(s1p06r,cy_r12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
				SSE2_fermat_carry_norm_errcheck(s1p07r,cy_r14,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck(s1p08r,cy_r16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck(s1p09r,cy_r18,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck(s1p10r,cy_r20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck(s1p11r,cy_r22,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck(s1p12r,cy_r24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck(s1p13r,cy_r26,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
				SSE2_fermat_carry_norm_errcheck(s1p14r,cy_r28,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck(s1p15r,cy_r30,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck(s1p16r,cy_r32,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck(s1p17r,cy_r34,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck(s1p18r,cy_r36,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck(s1p19r,cy_r38,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck(s1p20r,cy_r40,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
				SSE2_fermat_carry_norm_errcheck(s1p21r,cy_r42,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck(s1p22r,cy_r44,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck(s1p23r,cy_r46,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck(s1p24r,cy_r48,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck(s1p25r,cy_r50,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck(s1p26r,cy_r52,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck(s1p27r,cy_r54,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);

				SSE2_fermat_carry_norm_errcheck(s1p28r,cy_i00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck(s1p29r,cy_i02,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck(s1p30r,cy_i04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck(s1p31r,cy_i06,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck(s1p32r,cy_i08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck(s1p33r,cy_i10,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck(s1p34r,cy_i12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
				SSE2_fermat_carry_norm_errcheck(s1p35r,cy_i14,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck(s1p36r,cy_i16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck(s1p37r,cy_i18,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck(s1p38r,cy_i20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck(s1p39r,cy_i22,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck(s1p40r,cy_i24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck(s1p41r,cy_i26,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
				SSE2_fermat_carry_norm_errcheck(s1p42r,cy_i28,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck(s1p43r,cy_i30,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck(s1p44r,cy_i32,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck(s1p45r,cy_i34,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck(s1p46r,cy_i36,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck(s1p47r,cy_i38,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck(s1p48r,cy_i40,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
				SSE2_fermat_carry_norm_errcheck(s1p49r,cy_i42,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck(s1p50r,cy_i44,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck(s1p51r,cy_i46,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck(s1p52r,cy_i48,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck(s1p53r,cy_i50,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck(s1p54r,cy_i52,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck(s1p55r,cy_i54,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);

			  #else	// 64-bit SSE2

				SSE2_fermat_carry_norm_errcheck_X2(s1p00r,cy_r00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck_X2(s1p02r,cy_r04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck_X2(s1p04r,cy_r08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck_X2(s1p06r,cy_r12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck_X2(s1p08r,cy_r16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck_X2(s1p10r,cy_r20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck_X2(s1p12r,cy_r24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5,icycle6,jcycle6);
				SSE2_fermat_carry_norm_errcheck_X2(s1p14r,cy_r28,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck_X2(s1p16r,cy_r32,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck_X2(s1p18r,cy_r36,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck_X2(s1p20r,cy_r40,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck_X2(s1p22r,cy_r44,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck_X2(s1p24r,cy_r48,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck_X2(s1p26r,cy_r52,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5,icycle6,jcycle6);
				SSE2_fermat_carry_norm_errcheck_X2(s1p28r,cy_i00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck_X2(s1p30r,cy_i04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck_X2(s1p32r,cy_i08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck_X2(s1p34r,cy_i12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck_X2(s1p36r,cy_i16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck_X2(s1p38r,cy_i20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck_X2(s1p40r,cy_i24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5,icycle6,jcycle6);
				SSE2_fermat_carry_norm_errcheck_X2(s1p42r,cy_i28,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0,icycle1,jcycle1);
				SSE2_fermat_carry_norm_errcheck_X2(s1p44r,cy_i32,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2,icycle3,jcycle3);
				SSE2_fermat_carry_norm_errcheck_X2(s1p46r,cy_i36,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4,icycle5,jcycle5);
				SSE2_fermat_carry_norm_errcheck_X2(s1p48r,cy_i40,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6,icycle0,jcycle0);
				SSE2_fermat_carry_norm_errcheck_X2(s1p50r,cy_i44,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1,icycle2,jcycle2);
				SSE2_fermat_carry_norm_errcheck_X2(s1p52r,cy_i48,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3,icycle4,jcycle4);
				SSE2_fermat_carry_norm_errcheck_X2(s1p54r,cy_i52,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5,icycle6,jcycle6);
			  #endif

			#else	// Scalar-double mode:

				ntmp = 0;
				fermat_carry_norm_errcheckB(a1p00r,a1p00i,cy_r00,cy_i00,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p01r,a1p01i,cy_r01,cy_i01,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p02r,a1p02i,cy_r02,cy_i02,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p03r,a1p03i,cy_r03,cy_i03,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p04r,a1p04i,cy_r04,cy_i04,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p05r,a1p05i,cy_r05,cy_i05,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p06r,a1p06i,cy_r06,cy_i06,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p07r,a1p07i,cy_r07,cy_i07,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p08r,a1p08i,cy_r08,cy_i08,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p09r,a1p09i,cy_r09,cy_i09,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p10r,a1p10i,cy_r10,cy_i10,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p11r,a1p11i,cy_r11,cy_i11,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p12r,a1p12i,cy_r12,cy_i12,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p13r,a1p13i,cy_r13,cy_i13,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p14r,a1p14i,cy_r14,cy_i14,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p15r,a1p15i,cy_r15,cy_i15,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p16r,a1p16i,cy_r16,cy_i16,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p17r,a1p17i,cy_r17,cy_i17,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p18r,a1p18i,cy_r18,cy_i18,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p19r,a1p19i,cy_r19,cy_i19,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p20r,a1p20i,cy_r20,cy_i20,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p21r,a1p21i,cy_r21,cy_i21,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p22r,a1p22i,cy_r22,cy_i22,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p23r,a1p23i,cy_r23,cy_i23,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p24r,a1p24i,cy_r24,cy_i24,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p25r,a1p25i,cy_r25,cy_i25,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p26r,a1p26i,cy_r26,cy_i26,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p27r,a1p27i,cy_r27,cy_i27,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p28r,a1p28i,cy_r28,cy_i28,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p29r,a1p29i,cy_r29,cy_i29,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p30r,a1p30i,cy_r30,cy_i30,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p31r,a1p31i,cy_r31,cy_i31,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p32r,a1p32i,cy_r32,cy_i32,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p33r,a1p33i,cy_r33,cy_i33,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p34r,a1p34i,cy_r34,cy_i34,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p35r,a1p35i,cy_r35,cy_i35,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p36r,a1p36i,cy_r36,cy_i36,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p37r,a1p37i,cy_r37,cy_i37,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p38r,a1p38i,cy_r38,cy_i38,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p39r,a1p39i,cy_r39,cy_i39,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p40r,a1p40i,cy_r40,cy_i40,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p41r,a1p41i,cy_r41,cy_i41,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p42r,a1p42i,cy_r42,cy_i42,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p43r,a1p43i,cy_r43,cy_i43,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p44r,a1p44i,cy_r44,cy_i44,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p45r,a1p45i,cy_r45,cy_i45,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p46r,a1p46i,cy_r46,cy_i46,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p47r,a1p47i,cy_r47,cy_i47,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p48r,a1p48i,cy_r48,cy_i48,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p49r,a1p49i,cy_r49,cy_i49,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p50r,a1p50i,cy_r50,cy_i50,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p51r,a1p51i,cy_r51,cy_i51,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p52r,a1p52i,cy_r52,cy_i52,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p53r,a1p53i,cy_r53,cy_i53,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p54r,a1p54i,cy_r54,cy_i54,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p55r,a1p55i,cy_r55,cy_i55,icycle6,ntmp,NRTM1,NRT_BITS);

				icycle0 += wts_idx_incr;	/* Inside the loop use this, as it is faster than general-mod '% nwt' */
				icycle1 += wts_idx_incr;
				icycle2 += wts_idx_incr;
				icycle3 += wts_idx_incr;
				icycle4 += wts_idx_incr;
				icycle5 += wts_idx_incr;
				icycle6 += wts_idx_incr;
				icycle0 += ( (-(int)((uint32)icycle0 >> 31)) & nwt);
				icycle1 += ( (-(int)((uint32)icycle1 >> 31)) & nwt);
				icycle2 += ( (-(int)((uint32)icycle2 >> 31)) & nwt);
				icycle3 += ( (-(int)((uint32)icycle3 >> 31)) & nwt);
				icycle4 += ( (-(int)((uint32)icycle4 >> 31)) & nwt);
				icycle5 += ( (-(int)((uint32)icycle5 >> 31)) & nwt);
				icycle6 += ( (-(int)((uint32)icycle6 >> 31)) & nwt);

			#endif	/* #ifdef USE_SSE2 */

			// Here we nest AVX inside SSE2 since i/jcycle updates are for both, k/l for AVX-only:
			#ifdef USE_SSE2

				icycle0 += wts_idx_inc2;		icycle0 += ( (-(icycle0 < 0)) & nwt16);
				icycle1 += wts_idx_inc2;		icycle1 += ( (-(icycle1 < 0)) & nwt16);
				icycle2 += wts_idx_inc2;		icycle2 += ( (-(icycle2 < 0)) & nwt16);
				icycle3 += wts_idx_inc2;		icycle3 += ( (-(icycle3 < 0)) & nwt16);
				icycle4 += wts_idx_inc2;		icycle4 += ( (-(icycle4 < 0)) & nwt16);
				icycle5 += wts_idx_inc2;		icycle5 += ( (-(icycle5 < 0)) & nwt16);
				icycle6 += wts_idx_inc2;		icycle6 += ( (-(icycle6 < 0)) & nwt16);

				jcycle0 += wts_idx_inc2;		jcycle0 += ( (-(jcycle0 < 0)) & nwt16);
				jcycle1 += wts_idx_inc2;		jcycle1 += ( (-(jcycle1 < 0)) & nwt16);
				jcycle2 += wts_idx_inc2;		jcycle2 += ( (-(jcycle2 < 0)) & nwt16);
				jcycle3 += wts_idx_inc2;		jcycle3 += ( (-(jcycle3 < 0)) & nwt16);
				jcycle4 += wts_idx_inc2;		jcycle4 += ( (-(jcycle4 < 0)) & nwt16);
				jcycle5 += wts_idx_inc2;		jcycle5 += ( (-(jcycle5 < 0)) & nwt16);
				jcycle6 += wts_idx_inc2;		jcycle6 += ( (-(jcycle6 < 0)) & nwt16);

			  #ifdef USE_AVX
				kcycle0 += wts_idx_inc2;		kcycle0 += ( (-(kcycle0 < 0)) & nwt16);
				kcycle1 += wts_idx_inc2;		kcycle1 += ( (-(kcycle1 < 0)) & nwt16);
				kcycle2 += wts_idx_inc2;		kcycle2 += ( (-(kcycle2 < 0)) & nwt16);
				kcycle3 += wts_idx_inc2;		kcycle3 += ( (-(kcycle3 < 0)) & nwt16);
				kcycle4 += wts_idx_inc2;		kcycle4 += ( (-(kcycle4 < 0)) & nwt16);
				kcycle5 += wts_idx_inc2;		kcycle5 += ( (-(kcycle5 < 0)) & nwt16);
				kcycle6 += wts_idx_inc2;		kcycle6 += ( (-(kcycle6 < 0)) & nwt16);

				lcycle0 += wts_idx_inc2;		lcycle0 += ( (-(lcycle0 < 0)) & nwt16);
				lcycle1 += wts_idx_inc2;		lcycle1 += ( (-(lcycle1 < 0)) & nwt16);
				lcycle2 += wts_idx_inc2;		lcycle2 += ( (-(lcycle2 < 0)) & nwt16);
				lcycle3 += wts_idx_inc2;		lcycle3 += ( (-(lcycle3 < 0)) & nwt16);
				lcycle4 += wts_idx_inc2;		lcycle4 += ( (-(lcycle4 < 0)) & nwt16);
				lcycle5 += wts_idx_inc2;		lcycle5 += ( (-(lcycle5 < 0)) & nwt16);
				lcycle6 += wts_idx_inc2;		lcycle6 += ( (-(lcycle6 < 0)) & nwt16);
			  #endif
			#endif

			}	/* if(MODULUS_TYPE == ...) */

		/*...The radix-56 DIF pass is here:	*/

			#ifdef USE_SSE2

			/*...gather the needed data (56 64-bit complex, i.e. 112 64-bit reals) and do 8 radix-7 transforms...*/
							 /*                   inputs                    */ /* sincos */ /*         outputs           */
				SSE2_RADIX_07_DFT(s1p00r,s1p48r,s1p40r,s1p32r,s1p24r,s1p16r,s1p08r, cc0, r00r,r08r,r16r,r24r,r32r,r40r,r48r);
				SSE2_RADIX_07_DFT(s1p49r,s1p41r,s1p33r,s1p25r,s1p17r,s1p09r,s1p01r, cc0, r01r,r09r,r17r,r25r,r33r,r41r,r49r);
				SSE2_RADIX_07_DFT(s1p42r,s1p34r,s1p26r,s1p18r,s1p10r,s1p02r,s1p50r, cc0, r02r,r10r,r18r,r26r,r34r,r42r,r50r);
				SSE2_RADIX_07_DFT(s1p35r,s1p27r,s1p19r,s1p11r,s1p03r,s1p51r,s1p43r, cc0, r03r,r11r,r19r,r27r,r35r,r43r,r51r);
				SSE2_RADIX_07_DFT(s1p28r,s1p20r,s1p12r,s1p04r,s1p52r,s1p44r,s1p36r, cc0, r04r,r12r,r20r,r28r,r36r,r44r,r52r);
				SSE2_RADIX_07_DFT(s1p21r,s1p13r,s1p05r,s1p53r,s1p45r,s1p37r,s1p29r, cc0, r05r,r13r,r21r,r29r,r37r,r45r,r53r);
				SSE2_RADIX_07_DFT(s1p14r,s1p06r,s1p54r,s1p46r,s1p38r,s1p30r,s1p22r, cc0, r06r,r14r,r22r,r30r,r38r,r46r,r54r);
				SSE2_RADIX_07_DFT(s1p07r,s1p55r,s1p47r,s1p39r,s1p31r,s1p23r,s1p15r, cc0, r07r,r15r,r23r,r31r,r39r,r47r,r55r);

			/*...and now do 7 radix-8 transforms: */
							 /*                                   inputs                                  */ /*                       intermediates                       */ /*                 outputs                   */
				add0 = &a[j1    ];
				add1 = add0+p01;
				add2 = add0+p02;
				add3 = add0+p03;
				add4 = add0+p04;
				add5 = add0+p05;
				add6 = add0+p06;
				add7 = add0+p07;
			  #ifdef USE_AVX
				SSE2_RADIX8_DIF_0TWIDDLE(r00r,0x40,0x80,0xc0,0x100,0x140,0x180,0x1c0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #else
				SSE2_RADIX8_DIF_0TWIDDLE(r00r,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #endif

				add2 = &a[j1+p30];
				add0 = add2+p03;
				add1 = add2+p02;
				add3 = add2+p01;
				add4 = add2+p07;
				add5 = add2+p06;
				add6 = add2+p04;
				add7 = add2+p05;
			  #ifdef USE_AVX
				SSE2_RADIX8_DIF_0TWIDDLE(r08r,0x40,0x80,0xc0,0x100,0x140,0x180,0x1c0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #else
				SSE2_RADIX8_DIF_0TWIDDLE(r08r,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #endif

				add6 = &a[j1+p28];
				add0 = add6+p05;
				add1 = add6+p04;
				add2 = add6+p07;
				add3 = add6+p06;
				add4 = add6+p03;
				add5 = add6+p02;
				add7 = add6+p01;
			  #ifdef USE_AVX
				SSE2_RADIX8_DIF_0TWIDDLE(r16r,0x40,0x80,0xc0,0x100,0x140,0x180,0x1c0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #else
				SSE2_RADIX8_DIF_0TWIDDLE(r16r,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #endif

				add1 = &a[j1+p20];
				add0 = add1+p01;
				add2 = add1+p03;
				add3 = add1+p02;
				add4 = add1+p05;
				add5 = add1+p04;
				add6 = add1+p07;
				add7 = add1+p06;
			  #ifdef USE_AVX
				SSE2_RADIX8_DIF_0TWIDDLE(r24r,0x40,0x80,0xc0,0x100,0x140,0x180,0x1c0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #else
				SSE2_RADIX8_DIF_0TWIDDLE(r24r,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #endif

				add5 = &a[j1+p18];
				add0 = add5+p06;
				add1 = add5+p07;
				add2 = add5+p05;
				add3 = add5+p04;
				add4 = add5+p01;
				add6 = add5+p03;
				add7 = add5+p02;
			  #ifdef USE_AVX
				SSE2_RADIX8_DIF_0TWIDDLE(r32r,0x40,0x80,0xc0,0x100,0x140,0x180,0x1c0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #else
				SSE2_RADIX8_DIF_0TWIDDLE(r32r,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #endif

				add3 = &a[j1+p10];
				add0 = add3+p02;
				add1 = add3+p03;
				add2 = add3+p01;
				add4 = add3+p06;
				add5 = add3+p07;
				add6 = add3+p05;
				add7 = add3+p04;
			  #ifdef USE_AVX
				SSE2_RADIX8_DIF_0TWIDDLE(r40r,0x40,0x80,0xc0,0x100,0x140,0x180,0x1c0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #else
				SSE2_RADIX8_DIF_0TWIDDLE(r40r,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #endif

				add7 = &a[j1+p08];
				add0 = add7+p04;
				add1 = add7+p05;
				add2 = add7+p06;
				add3 = add7+p07;
				add4 = add7+p02;
				add5 = add7+p03;
				add6 = add7+p01;
			  #ifdef USE_AVX
				SSE2_RADIX8_DIF_0TWIDDLE(r48r,0x40,0x80,0xc0,0x100,0x140,0x180,0x1c0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #else
				SSE2_RADIX8_DIF_0TWIDDLE(r48r,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #endif

			#else	// USE_SSE2 = False:

			/*...gather the needed data (56 64-bit complex, i.e. 112 64-bit reals) and do 8 radix-7 transforms...*/
							 /*                                               inputs                                              */  /*                   intermediates                   */  /*                                                  outputs        */  /*   sincos consts   */
				RADIX_07_DFT(a1p00r,a1p00i, a1p48r,a1p48i, a1p40r,a1p40i, a1p32r,a1p32i, a1p24r,a1p24i, a1p16r,a1p16i, a1p08r,a1p08i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r00r,r00i,r08r,r08i,r16r,r16i,r24r,r24i,r32r,r32i,r40r,r40i,r48r,r48i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(a1p49r,a1p49i, a1p41r,a1p41i, a1p33r,a1p33i, a1p25r,a1p25i, a1p17r,a1p17i, a1p09r,a1p09i, a1p01r,a1p01i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r01r,r01i,r09r,r09i,r17r,r17i,r25r,r25i,r33r,r33i,r41r,r41i,r49r,r49i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(a1p42r,a1p42i, a1p34r,a1p34i, a1p26r,a1p26i, a1p18r,a1p18i, a1p10r,a1p10i, a1p02r,a1p02i, a1p50r,a1p50i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r02r,r02i,r10r,r10i,r18r,r18i,r26r,r26i,r34r,r34i,r42r,r42i,r50r,r50i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(a1p35r,a1p35i, a1p27r,a1p27i, a1p19r,a1p19i, a1p11r,a1p11i, a1p03r,a1p03i, a1p51r,a1p51i, a1p43r,a1p43i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r03r,r03i,r11r,r11i,r19r,r19i,r27r,r27i,r35r,r35i,r43r,r43i,r51r,r51i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(a1p28r,a1p28i, a1p20r,a1p20i, a1p12r,a1p12i, a1p04r,a1p04i, a1p52r,a1p52i, a1p44r,a1p44i, a1p36r,a1p36i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r04r,r04i,r12r,r12i,r20r,r20i,r28r,r28i,r36r,r36i,r44r,r44i,r52r,r52i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(a1p21r,a1p21i, a1p13r,a1p13i, a1p05r,a1p05i, a1p53r,a1p53i, a1p45r,a1p45i, a1p37r,a1p37i, a1p29r,a1p29i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r05r,r05i,r13r,r13i,r21r,r21i,r29r,r29i,r37r,r37i,r45r,r45i,r53r,r53i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(a1p14r,a1p14i, a1p06r,a1p06i, a1p54r,a1p54i, a1p46r,a1p46i, a1p38r,a1p38i, a1p30r,a1p30i, a1p22r,a1p22i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r06r,r06i,r14r,r14i,r22r,r22i,r30r,r30i,r38r,r38i,r46r,r46i,r54r,r54i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(a1p07r,a1p07i, a1p55r,a1p55i, a1p47r,a1p47i, a1p39r,a1p39i, a1p31r,a1p31i, a1p23r,a1p23i, a1p15r,a1p15i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r07r,r07i,r15r,r15i,r23r,r23i,r31r,r31i,r39r,r39i,r47r,r47i,r55r,r55i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);

			/*...and now do 7 radix-8 transforms: */
							 /*                                   inputs                                  */ /*                       intermediates                       */ /*                 outputs                   */
				RADIX_08_DIF(r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05],a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],rt,it);	jt = j1+p30; jp = j2+p30;
				RADIX_08_DIF(r08r,r08i,r09r,r09i,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],rt,it);	jt = j1+p28; jp = j2+p28;
				RADIX_08_DIF(r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	jt = j1+p20; jp = j2+p20;
				RADIX_08_DIF(r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i,r30r,r30i,r31r,r31i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],rt,it);	jt = j1+p18; jp = j2+p18;
				RADIX_08_DIF(r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p10; jp = j2+p10;
				RADIX_08_DIF(r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],rt,it);	jt = j1+p08; jp = j2+p08;
				RADIX_08_DIF(r48r,r48i,r49r,r49i,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);

			#endif	/* !USE_SSE2 */

			}

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
			_cy_r00[ithread] = cy_r00->d0;	_cy_r01[ithread] = cy_r00->d1;	_cy_r02[ithread] = cy_r00->d2;	_cy_r03[ithread] = cy_r00->d3;
			_cy_r04[ithread] = cy_r04->d0;	_cy_r05[ithread] = cy_r04->d1;	_cy_r06[ithread] = cy_r04->d2;	_cy_r07[ithread] = cy_r04->d3;
			_cy_r08[ithread] = cy_r08->d0;	_cy_r09[ithread] = cy_r08->d1;	_cy_r10[ithread] = cy_r08->d2;	_cy_r11[ithread] = cy_r08->d3;
			_cy_r12[ithread] = cy_r12->d0;	_cy_r13[ithread] = cy_r12->d1;	_cy_r14[ithread] = cy_r12->d2;	_cy_r15[ithread] = cy_r12->d3;
			_cy_r16[ithread] = cy_r16->d0;	_cy_r17[ithread] = cy_r16->d1;	_cy_r18[ithread] = cy_r16->d2;	_cy_r19[ithread] = cy_r16->d3;
			_cy_r20[ithread] = cy_r20->d0;	_cy_r21[ithread] = cy_r20->d1;	_cy_r22[ithread] = cy_r20->d2;	_cy_r23[ithread] = cy_r20->d3;
			_cy_r24[ithread] = cy_r24->d0;	_cy_r25[ithread] = cy_r24->d1;	_cy_r26[ithread] = cy_r24->d2;	_cy_r27[ithread] = cy_r24->d3;
			_cy_r28[ithread] = cy_r28->d0;	_cy_r29[ithread] = cy_r28->d1;	_cy_r30[ithread] = cy_r28->d2;	_cy_r31[ithread] = cy_r28->d3;
			_cy_r32[ithread] = cy_r32->d0;	_cy_r33[ithread] = cy_r32->d1;	_cy_r34[ithread] = cy_r32->d2;	_cy_r35[ithread] = cy_r32->d3;
			_cy_r36[ithread] = cy_r36->d0;	_cy_r37[ithread] = cy_r36->d1;	_cy_r38[ithread] = cy_r36->d2;	_cy_r39[ithread] = cy_r36->d3;
			_cy_r40[ithread] = cy_r40->d0;	_cy_r41[ithread] = cy_r40->d1;	_cy_r42[ithread] = cy_r40->d2;	_cy_r43[ithread] = cy_r40->d3;
			_cy_r44[ithread] = cy_r44->d0;	_cy_r45[ithread] = cy_r44->d1;	_cy_r46[ithread] = cy_r44->d2;	_cy_r47[ithread] = cy_r44->d3;
			_cy_r48[ithread] = cy_r48->d0;	_cy_r49[ithread] = cy_r48->d1;	_cy_r50[ithread] = cy_r48->d2;	_cy_r51[ithread] = cy_r48->d3;
			_cy_r52[ithread] = cy_r52->d0;	_cy_r53[ithread] = cy_r52->d1;	_cy_r54[ithread] = cy_r52->d2;	_cy_r55[ithread] = cy_r52->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
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
			_cy_r28[ithread] = cy_r28->d0;	_cy_r29[ithread] = cy_r28->d1;
			_cy_r30[ithread] = cy_r30->d0;	_cy_r31[ithread] = cy_r30->d1;
			_cy_r32[ithread] = cy_r32->d0;	_cy_r33[ithread] = cy_r32->d1;
			_cy_r34[ithread] = cy_r34->d0;	_cy_r35[ithread] = cy_r34->d1;
			_cy_r36[ithread] = cy_r36->d0;	_cy_r37[ithread] = cy_r36->d1;
			_cy_r38[ithread] = cy_r38->d0;	_cy_r39[ithread] = cy_r38->d1;
			_cy_r40[ithread] = cy_r40->d0;	_cy_r41[ithread] = cy_r40->d1;
			_cy_r42[ithread] = cy_r42->d0;	_cy_r43[ithread] = cy_r42->d1;
			_cy_r44[ithread] = cy_r44->d0;	_cy_r45[ithread] = cy_r44->d1;
			_cy_r46[ithread] = cy_r46->d0;	_cy_r47[ithread] = cy_r46->d1;
			_cy_r48[ithread] = cy_r48->d0;	_cy_r49[ithread] = cy_r48->d1;
			_cy_r50[ithread] = cy_r50->d0;	_cy_r51[ithread] = cy_r50->d1;
			_cy_r52[ithread] = cy_r52->d0;	_cy_r53[ithread] = cy_r52->d1;
			_cy_r54[ithread] = cy_r54->d0;	_cy_r55[ithread] = cy_r54->d1;
			maxerr = MAX(max_err->d0,max_err->d1);
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
			_cy_r28[ithread] = cy_r28;
			_cy_r29[ithread] = cy_r29;
			_cy_r30[ithread] = cy_r30;
			_cy_r31[ithread] = cy_r31;
			_cy_r32[ithread] = cy_r32;
			_cy_r33[ithread] = cy_r33;
			_cy_r34[ithread] = cy_r34;
			_cy_r35[ithread] = cy_r35;
			_cy_r36[ithread] = cy_r36;
			_cy_r37[ithread] = cy_r37;
			_cy_r38[ithread] = cy_r38;
			_cy_r39[ithread] = cy_r39;
			_cy_r40[ithread] = cy_r40;
			_cy_r41[ithread] = cy_r41;
			_cy_r42[ithread] = cy_r42;
			_cy_r43[ithread] = cy_r43;
			_cy_r44[ithread] = cy_r44;
			_cy_r45[ithread] = cy_r45;
			_cy_r46[ithread] = cy_r46;
			_cy_r47[ithread] = cy_r47;
			_cy_r48[ithread] = cy_r48;
			_cy_r49[ithread] = cy_r49;
			_cy_r50[ithread] = cy_r50;
			_cy_r51[ithread] = cy_r51;
			_cy_r52[ithread] = cy_r52;
			_cy_r53[ithread] = cy_r53;
			_cy_r54[ithread] = cy_r54;
			_cy_r55[ithread] = cy_r55;
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
			_cy_r28[ithread] = cy_r28->d0;	_cy_r29[ithread] = cy_r28->d1;	_cy_r30[ithread] = cy_r28->d2;	_cy_r31[ithread] = cy_r28->d3;
			_cy_r32[ithread] = cy_r32->d0;	_cy_r33[ithread] = cy_r32->d1;	_cy_r34[ithread] = cy_r32->d2;	_cy_r35[ithread] = cy_r32->d3;
			_cy_r36[ithread] = cy_r36->d0;	_cy_r37[ithread] = cy_r36->d1;	_cy_r38[ithread] = cy_r36->d2;	_cy_r39[ithread] = cy_r36->d3;
			_cy_r40[ithread] = cy_r40->d0;	_cy_r41[ithread] = cy_r40->d1;	_cy_r42[ithread] = cy_r40->d2;	_cy_r43[ithread] = cy_r40->d3;
			_cy_r44[ithread] = cy_r44->d0;	_cy_r45[ithread] = cy_r44->d1;	_cy_r46[ithread] = cy_r44->d2;	_cy_r47[ithread] = cy_r44->d3;
			_cy_r48[ithread] = cy_r48->d0;	_cy_r49[ithread] = cy_r48->d1;	_cy_r50[ithread] = cy_r48->d2;	_cy_r51[ithread] = cy_r48->d3;
			_cy_r52[ithread] = cy_r52->d0;	_cy_r53[ithread] = cy_r52->d1;	_cy_r54[ithread] = cy_r52->d2;	_cy_r55[ithread] = cy_r52->d3;

			_cy_i00[ithread] = cy_i00->d0;	_cy_i01[ithread] = cy_i00->d1;	_cy_i02[ithread] = cy_i00->d2;	_cy_i03[ithread] = cy_i00->d3;
			_cy_i04[ithread] = cy_i04->d0;	_cy_i05[ithread] = cy_i04->d1;	_cy_i06[ithread] = cy_i04->d2;	_cy_i07[ithread] = cy_i04->d3;
			_cy_i08[ithread] = cy_i08->d0;	_cy_i09[ithread] = cy_i08->d1;	_cy_i10[ithread] = cy_i08->d2;	_cy_i11[ithread] = cy_i08->d3;
			_cy_i12[ithread] = cy_i12->d0;	_cy_i13[ithread] = cy_i12->d1;	_cy_i14[ithread] = cy_i12->d2;	_cy_i15[ithread] = cy_i12->d3;
			_cy_i16[ithread] = cy_i16->d0;	_cy_i17[ithread] = cy_i16->d1;	_cy_i18[ithread] = cy_i16->d2;	_cy_i19[ithread] = cy_i16->d3;
			_cy_i20[ithread] = cy_i20->d0;	_cy_i21[ithread] = cy_i20->d1;	_cy_i22[ithread] = cy_i20->d2;	_cy_i23[ithread] = cy_i20->d3;
			_cy_i24[ithread] = cy_i24->d0;	_cy_i25[ithread] = cy_i24->d1;	_cy_i26[ithread] = cy_i24->d2;	_cy_i27[ithread] = cy_i24->d3;
			_cy_i28[ithread] = cy_i28->d0;	_cy_i29[ithread] = cy_i28->d1;	_cy_i30[ithread] = cy_i28->d2;	_cy_i31[ithread] = cy_i28->d3;
			_cy_i32[ithread] = cy_i32->d0;	_cy_i33[ithread] = cy_i32->d1;	_cy_i34[ithread] = cy_i32->d2;	_cy_i35[ithread] = cy_i32->d3;
			_cy_i36[ithread] = cy_i36->d0;	_cy_i37[ithread] = cy_i36->d1;	_cy_i38[ithread] = cy_i36->d2;	_cy_i39[ithread] = cy_i36->d3;
			_cy_i40[ithread] = cy_i40->d0;	_cy_i41[ithread] = cy_i40->d1;	_cy_i42[ithread] = cy_i40->d2;	_cy_i43[ithread] = cy_i40->d3;
			_cy_i44[ithread] = cy_i44->d0;	_cy_i45[ithread] = cy_i44->d1;	_cy_i46[ithread] = cy_i44->d2;	_cy_i47[ithread] = cy_i44->d3;
			_cy_i48[ithread] = cy_i48->d0;	_cy_i49[ithread] = cy_i48->d1;	_cy_i50[ithread] = cy_i48->d2;	_cy_i51[ithread] = cy_i48->d3;
			_cy_i52[ithread] = cy_i52->d0;	_cy_i53[ithread] = cy_i52->d1;	_cy_i54[ithread] = cy_i52->d2;	_cy_i55[ithread] = cy_i52->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
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
			_cy_r14[ithread] = cy_r28->d0;	_cy_i14[ithread] = cy_r28->d1;
			_cy_r15[ithread] = cy_r30->d0;	_cy_i15[ithread] = cy_r30->d1;
			_cy_r16[ithread] = cy_r32->d0;	_cy_i16[ithread] = cy_r32->d1;
			_cy_r17[ithread] = cy_r34->d0;	_cy_i17[ithread] = cy_r34->d1;
			_cy_r18[ithread] = cy_r36->d0;	_cy_i18[ithread] = cy_r36->d1;
			_cy_r19[ithread] = cy_r38->d0;	_cy_i19[ithread] = cy_r38->d1;
			_cy_r20[ithread] = cy_r40->d0;	_cy_i20[ithread] = cy_r40->d1;
			_cy_r21[ithread] = cy_r42->d0;	_cy_i21[ithread] = cy_r42->d1;
			_cy_r22[ithread] = cy_r44->d0;	_cy_i22[ithread] = cy_r44->d1;
			_cy_r23[ithread] = cy_r46->d0;	_cy_i23[ithread] = cy_r46->d1;
			_cy_r24[ithread] = cy_r48->d0;	_cy_i24[ithread] = cy_r48->d1;
			_cy_r25[ithread] = cy_r50->d0;	_cy_i25[ithread] = cy_r50->d1;
			_cy_r26[ithread] = cy_r52->d0;	_cy_i26[ithread] = cy_r52->d1;
			_cy_r27[ithread] = cy_r54->d0;	_cy_i27[ithread] = cy_r54->d1;
			_cy_r28[ithread] = cy_i00->d0;	_cy_i28[ithread] = cy_i00->d1;
			_cy_r29[ithread] = cy_i02->d0;	_cy_i29[ithread] = cy_i02->d1;
			_cy_r30[ithread] = cy_i04->d0;	_cy_i30[ithread] = cy_i04->d1;
			_cy_r31[ithread] = cy_i06->d0;	_cy_i31[ithread] = cy_i06->d1;
			_cy_r32[ithread] = cy_i08->d0;	_cy_i32[ithread] = cy_i08->d1;
			_cy_r33[ithread] = cy_i10->d0;	_cy_i33[ithread] = cy_i10->d1;
			_cy_r34[ithread] = cy_i12->d0;	_cy_i34[ithread] = cy_i12->d1;
			_cy_r35[ithread] = cy_i14->d0;	_cy_i35[ithread] = cy_i14->d1;
			_cy_r36[ithread] = cy_i16->d0;	_cy_i36[ithread] = cy_i16->d1;
			_cy_r37[ithread] = cy_i18->d0;	_cy_i37[ithread] = cy_i18->d1;
			_cy_r38[ithread] = cy_i20->d0;	_cy_i38[ithread] = cy_i20->d1;
			_cy_r39[ithread] = cy_i22->d0;	_cy_i39[ithread] = cy_i22->d1;
			_cy_r40[ithread] = cy_i24->d0;	_cy_i40[ithread] = cy_i24->d1;
			_cy_r41[ithread] = cy_i26->d0;	_cy_i41[ithread] = cy_i26->d1;
			_cy_r42[ithread] = cy_i28->d0;	_cy_i42[ithread] = cy_i28->d1;
			_cy_r43[ithread] = cy_i30->d0;	_cy_i43[ithread] = cy_i30->d1;
			_cy_r44[ithread] = cy_i32->d0;	_cy_i44[ithread] = cy_i32->d1;
			_cy_r45[ithread] = cy_i34->d0;	_cy_i45[ithread] = cy_i34->d1;
			_cy_r46[ithread] = cy_i36->d0;	_cy_i46[ithread] = cy_i36->d1;
			_cy_r47[ithread] = cy_i38->d0;	_cy_i47[ithread] = cy_i38->d1;
			_cy_r48[ithread] = cy_i40->d0;	_cy_i48[ithread] = cy_i40->d1;
			_cy_r49[ithread] = cy_i42->d0;	_cy_i49[ithread] = cy_i42->d1;
			_cy_r50[ithread] = cy_i44->d0;	_cy_i50[ithread] = cy_i44->d1;
			_cy_r51[ithread] = cy_i46->d0;	_cy_i51[ithread] = cy_i46->d1;
			_cy_r52[ithread] = cy_i48->d0;	_cy_i52[ithread] = cy_i48->d1;
			_cy_r53[ithread] = cy_i50->d0;	_cy_i53[ithread] = cy_i50->d1;
			_cy_r54[ithread] = cy_i52->d0;	_cy_i54[ithread] = cy_i52->d1;
			_cy_r55[ithread] = cy_i54->d0;	_cy_i55[ithread] = cy_i54->d1;
			maxerr = MAX(max_err->d0,max_err->d1);
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
			_cy_r28[ithread] = cy_r28;	_cy_i28[ithread] = cy_i28;
			_cy_r29[ithread] = cy_r29;	_cy_i29[ithread] = cy_i29;
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
			_cy_r40[ithread] = cy_r40;	_cy_i40[ithread] = cy_i40;
			_cy_r41[ithread] = cy_r41;	_cy_i41[ithread] = cy_i41;
			_cy_r42[ithread] = cy_r42;	_cy_i42[ithread] = cy_i42;
			_cy_r43[ithread] = cy_r43;	_cy_i43[ithread] = cy_i43;
			_cy_r44[ithread] = cy_r44;	_cy_i44[ithread] = cy_i44;
			_cy_r45[ithread] = cy_r45;	_cy_i45[ithread] = cy_i45;
			_cy_r46[ithread] = cy_r46;	_cy_i46[ithread] = cy_i46;
			_cy_r47[ithread] = cy_r47;	_cy_i47[ithread] = cy_i47;
			_cy_r48[ithread] = cy_r48;	_cy_i48[ithread] = cy_i48;
			_cy_r49[ithread] = cy_r49;	_cy_i49[ithread] = cy_i49;
			_cy_r50[ithread] = cy_r50;	_cy_i50[ithread] = cy_i50;
			_cy_r51[ithread] = cy_r51;	_cy_i51[ithread] = cy_i51;
			_cy_r52[ithread] = cy_r52;	_cy_i52[ithread] = cy_i52;
			_cy_r53[ithread] = cy_r53;	_cy_i53[ithread] = cy_i53;
			_cy_r54[ithread] = cy_r54;	_cy_i54[ithread] = cy_i54;
			_cy_r55[ithread] = cy_r55;	_cy_i55[ithread] = cy_i55;
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
		ASSERT(HERE, 0x0 == cy56_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;
	ns_time.tv_sec  = 0.0001;// (time_t)seconds
	ns_time.tv_nsec = 0;	// (long)nanoseconds - At least allegedly, but under OS X it seems to be finer-grained than that

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(HERE, 0 == nanosleep(&ns_time, 0x0), "nanosleep fail!");
	}

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = tdat[ithread].maxerr;
		if(maxerr < _maxerr[ithread]) {
			maxerr = _maxerr[ithread];
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
			_cy_r28[ithread] = tdat[ithread].cy_r28;
			_cy_r29[ithread] = tdat[ithread].cy_r29;
			_cy_r30[ithread] = tdat[ithread].cy_r30;
			_cy_r31[ithread] = tdat[ithread].cy_r31;
			_cy_r32[ithread] = tdat[ithread].cy_r32;
			_cy_r33[ithread] = tdat[ithread].cy_r33;
			_cy_r34[ithread] = tdat[ithread].cy_r34;
			_cy_r35[ithread] = tdat[ithread].cy_r35;
			_cy_r36[ithread] = tdat[ithread].cy_r36;
			_cy_r37[ithread] = tdat[ithread].cy_r37;
			_cy_r38[ithread] = tdat[ithread].cy_r38;
			_cy_r39[ithread] = tdat[ithread].cy_r39;
			_cy_r40[ithread] = tdat[ithread].cy_r40;
			_cy_r41[ithread] = tdat[ithread].cy_r41;
			_cy_r42[ithread] = tdat[ithread].cy_r42;
			_cy_r43[ithread] = tdat[ithread].cy_r43;
			_cy_r44[ithread] = tdat[ithread].cy_r44;
			_cy_r45[ithread] = tdat[ithread].cy_r45;
			_cy_r46[ithread] = tdat[ithread].cy_r46;
			_cy_r47[ithread] = tdat[ithread].cy_r47;
			_cy_r48[ithread] = tdat[ithread].cy_r48;
			_cy_r49[ithread] = tdat[ithread].cy_r49;
			_cy_r50[ithread] = tdat[ithread].cy_r50;
			_cy_r51[ithread] = tdat[ithread].cy_r51;
			_cy_r52[ithread] = tdat[ithread].cy_r52;
			_cy_r53[ithread] = tdat[ithread].cy_r53;
			_cy_r54[ithread] = tdat[ithread].cy_r54;
			_cy_r55[ithread] = tdat[ithread].cy_r55;
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
			_cy_r28[ithread] = tdat[ithread].cy_r28;	_cy_i28[ithread] = tdat[ithread].cy_i28;
			_cy_r29[ithread] = tdat[ithread].cy_r29;	_cy_i29[ithread] = tdat[ithread].cy_i29;
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
			_cy_r40[ithread] = tdat[ithread].cy_r40;	_cy_i40[ithread] = tdat[ithread].cy_i40;
			_cy_r41[ithread] = tdat[ithread].cy_r41;	_cy_i41[ithread] = tdat[ithread].cy_i41;
			_cy_r42[ithread] = tdat[ithread].cy_r42;	_cy_i42[ithread] = tdat[ithread].cy_i42;
			_cy_r43[ithread] = tdat[ithread].cy_r43;	_cy_i43[ithread] = tdat[ithread].cy_i43;
			_cy_r44[ithread] = tdat[ithread].cy_r44;	_cy_i44[ithread] = tdat[ithread].cy_i44;
			_cy_r45[ithread] = tdat[ithread].cy_r45;	_cy_i45[ithread] = tdat[ithread].cy_i45;
			_cy_r46[ithread] = tdat[ithread].cy_r46;	_cy_i46[ithread] = tdat[ithread].cy_i46;
			_cy_r47[ithread] = tdat[ithread].cy_r47;	_cy_i47[ithread] = tdat[ithread].cy_i47;
			_cy_r48[ithread] = tdat[ithread].cy_r48;	_cy_i48[ithread] = tdat[ithread].cy_i48;
			_cy_r49[ithread] = tdat[ithread].cy_r49;	_cy_i49[ithread] = tdat[ithread].cy_i49;
			_cy_r50[ithread] = tdat[ithread].cy_r50;	_cy_i50[ithread] = tdat[ithread].cy_i50;
			_cy_r51[ithread] = tdat[ithread].cy_r51;	_cy_i51[ithread] = tdat[ithread].cy_i51;
			_cy_r52[ithread] = tdat[ithread].cy_r52;	_cy_i52[ithread] = tdat[ithread].cy_i52;
			_cy_r53[ithread] = tdat[ithread].cy_r53;	_cy_i53[ithread] = tdat[ithread].cy_i53;
			_cy_r54[ithread] = tdat[ithread].cy_r54;	_cy_i54[ithread] = tdat[ithread].cy_i54;
			_cy_r55[ithread] = tdat[ithread].cy_r55;	_cy_i55[ithread] = tdat[ithread].cy_i55;
		}

	}
#endif

	if(full_pass) {
	//	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
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
		t00	= _cy_r00[CY_THREADS - 1];
		t02	= _cy_r01[CY_THREADS - 1];
		t04	= _cy_r02[CY_THREADS - 1];
		t06	= _cy_r03[CY_THREADS - 1];
		t08	= _cy_r04[CY_THREADS - 1];
		t0A	= _cy_r05[CY_THREADS - 1];
		t0C	= _cy_r06[CY_THREADS - 1];
		t0E	= _cy_r07[CY_THREADS - 1];
		t10	= _cy_r08[CY_THREADS - 1];
		t12	= _cy_r09[CY_THREADS - 1];
		t14	= _cy_r10[CY_THREADS - 1];
		t16	= _cy_r11[CY_THREADS - 1];
		t18	= _cy_r12[CY_THREADS - 1];
		t1A	= _cy_r13[CY_THREADS - 1];
		t1C	= _cy_r14[CY_THREADS - 1];
		t1E	= _cy_r15[CY_THREADS - 1];
		t20	= _cy_r16[CY_THREADS - 1];
		t22	= _cy_r17[CY_THREADS - 1];
		t24	= _cy_r18[CY_THREADS - 1];
		t26	= _cy_r19[CY_THREADS - 1];
		t28	= _cy_r20[CY_THREADS - 1];
		t2A	= _cy_r21[CY_THREADS - 1];
		t2C	= _cy_r22[CY_THREADS - 1];
		t2E	= _cy_r23[CY_THREADS - 1];
		t30	= _cy_r24[CY_THREADS - 1];
		t32	= _cy_r25[CY_THREADS - 1];
		t34	= _cy_r26[CY_THREADS - 1];
		t36	= _cy_r27[CY_THREADS - 1];
		t38	= _cy_r28[CY_THREADS - 1];
		t3A	= _cy_r29[CY_THREADS - 1];
		t3C	= _cy_r30[CY_THREADS - 1];
		t3E	= _cy_r31[CY_THREADS - 1];
		t40	= _cy_r32[CY_THREADS - 1];
		t42	= _cy_r33[CY_THREADS - 1];
		t44	= _cy_r34[CY_THREADS - 1];
		t46	= _cy_r35[CY_THREADS - 1];
		t48	= _cy_r36[CY_THREADS - 1];
		t4A	= _cy_r37[CY_THREADS - 1];
		t4C	= _cy_r38[CY_THREADS - 1];
		t4E	= _cy_r39[CY_THREADS - 1];
		t50	= _cy_r40[CY_THREADS - 1];
		t52	= _cy_r41[CY_THREADS - 1];
		t54	= _cy_r42[CY_THREADS - 1];
		t56	= _cy_r43[CY_THREADS - 1];
		t58	= _cy_r44[CY_THREADS - 1];
		t5A	= _cy_r45[CY_THREADS - 1];
		t5C	= _cy_r46[CY_THREADS - 1];
		t5E	= _cy_r47[CY_THREADS - 1];
		t60	= _cy_r48[CY_THREADS - 1];
		t62	= _cy_r49[CY_THREADS - 1];
		t64	= _cy_r50[CY_THREADS - 1];
		t66	= _cy_r51[CY_THREADS - 1];
		t68	= _cy_r52[CY_THREADS - 1];
		t6A	= _cy_r53[CY_THREADS - 1];
		t6C	= _cy_r54[CY_THREADS - 1];
		t6E	= _cy_r55[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1, "CY_THREADS must be > 1!");	/* Make sure loop only gets executed if multiple threads */
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
			_cy_r28[ithread] = _cy_r28[ithread-1];
			_cy_r29[ithread] = _cy_r29[ithread-1];
			_cy_r30[ithread] = _cy_r30[ithread-1];
			_cy_r31[ithread] = _cy_r31[ithread-1];
			_cy_r32[ithread] = _cy_r32[ithread-1];
			_cy_r33[ithread] = _cy_r33[ithread-1];
			_cy_r34[ithread] = _cy_r34[ithread-1];
			_cy_r35[ithread] = _cy_r35[ithread-1];
			_cy_r36[ithread] = _cy_r36[ithread-1];
			_cy_r37[ithread] = _cy_r37[ithread-1];
			_cy_r38[ithread] = _cy_r38[ithread-1];
			_cy_r39[ithread] = _cy_r39[ithread-1];
			_cy_r40[ithread] = _cy_r40[ithread-1];
			_cy_r41[ithread] = _cy_r41[ithread-1];
			_cy_r42[ithread] = _cy_r42[ithread-1];
			_cy_r43[ithread] = _cy_r43[ithread-1];
			_cy_r44[ithread] = _cy_r44[ithread-1];
			_cy_r45[ithread] = _cy_r45[ithread-1];
			_cy_r46[ithread] = _cy_r46[ithread-1];
			_cy_r47[ithread] = _cy_r47[ithread-1];
			_cy_r48[ithread] = _cy_r48[ithread-1];
			_cy_r49[ithread] = _cy_r49[ithread-1];
			_cy_r50[ithread] = _cy_r50[ithread-1];
			_cy_r51[ithread] = _cy_r51[ithread-1];
			_cy_r52[ithread] = _cy_r52[ithread-1];
			_cy_r53[ithread] = _cy_r53[ithread-1];
			_cy_r54[ithread] = _cy_r54[ithread-1];
			_cy_r55[ithread] = _cy_r55[ithread-1];
		}

		_cy_r00[0] =+t6E;	/* ...The wraparound carry is here: */
		_cy_r01[0] = t00;
		_cy_r02[0] = t02;
		_cy_r03[0] = t04;
		_cy_r04[0] = t06;
		_cy_r05[0] = t08;
		_cy_r06[0] = t0A;
		_cy_r07[0] = t0C;
		_cy_r08[0] = t0E;
		_cy_r09[0] = t10;
		_cy_r10[0] = t12;
		_cy_r11[0] = t14;
		_cy_r12[0] = t16;
		_cy_r13[0] = t18;
		_cy_r14[0] = t1A;
		_cy_r15[0] = t1C;
		_cy_r16[0] = t1E;
		_cy_r17[0] = t20;
		_cy_r18[0] = t22;
		_cy_r19[0] = t24;
		_cy_r20[0] = t26;
		_cy_r21[0] = t28;
		_cy_r22[0] = t2A;
		_cy_r23[0] = t2C;
		_cy_r24[0] = t2E;
		_cy_r25[0] = t30;
		_cy_r26[0] = t32;
		_cy_r27[0] = t34;
		_cy_r28[0] = t36;
		_cy_r29[0] = t38;
		_cy_r30[0] = t3A;
		_cy_r31[0] = t3C;
		_cy_r32[0] = t3E;
		_cy_r33[0] = t40;
		_cy_r34[0] = t42;
		_cy_r35[0] = t44;
		_cy_r36[0] = t46;
		_cy_r37[0] = t48;
		_cy_r38[0] = t4A;
		_cy_r39[0] = t4C;
		_cy_r40[0] = t4E;
		_cy_r41[0] = t50;
		_cy_r42[0] = t52;
		_cy_r43[0] = t54;
		_cy_r44[0] = t56;
		_cy_r45[0] = t58;
		_cy_r46[0] = t5A;
		_cy_r47[0] = t5C;
		_cy_r48[0] = t5E;
		_cy_r49[0] = t60;
		_cy_r50[0] = t62;
		_cy_r51[0] = t64;
		_cy_r52[0] = t66;
		_cy_r53[0] = t68;
		_cy_r54[0] = t6A;
		_cy_r55[0] = t6C;
	}
	else
	{
		t00	= _cy_r00[CY_THREADS - 1];	t01 = _cy_i00[CY_THREADS - 1];
		t02	= _cy_r01[CY_THREADS - 1];	t03 = _cy_i01[CY_THREADS - 1];
		t04	= _cy_r02[CY_THREADS - 1];	t05 = _cy_i02[CY_THREADS - 1];
		t06	= _cy_r03[CY_THREADS - 1];	t07 = _cy_i03[CY_THREADS - 1];
		t08	= _cy_r04[CY_THREADS - 1];	t09 = _cy_i04[CY_THREADS - 1];
		t0A	= _cy_r05[CY_THREADS - 1];	t0B = _cy_i05[CY_THREADS - 1];
		t0C	= _cy_r06[CY_THREADS - 1];	t0D = _cy_i06[CY_THREADS - 1];
		t0E	= _cy_r07[CY_THREADS - 1];	t0F = _cy_i07[CY_THREADS - 1];
		t10	= _cy_r08[CY_THREADS - 1];	t11 = _cy_i08[CY_THREADS - 1];
		t12	= _cy_r09[CY_THREADS - 1];	t13 = _cy_i09[CY_THREADS - 1];
		t14	= _cy_r10[CY_THREADS - 1];	t15 = _cy_i10[CY_THREADS - 1];
		t16	= _cy_r11[CY_THREADS - 1];	t17 = _cy_i11[CY_THREADS - 1];
		t18	= _cy_r12[CY_THREADS - 1];	t19 = _cy_i12[CY_THREADS - 1];
		t1A	= _cy_r13[CY_THREADS - 1];	t1B = _cy_i13[CY_THREADS - 1];
		t1C	= _cy_r14[CY_THREADS - 1];	t1D = _cy_i14[CY_THREADS - 1];
		t1E	= _cy_r15[CY_THREADS - 1];	t1F = _cy_i15[CY_THREADS - 1];
		t20	= _cy_r16[CY_THREADS - 1];	t21 = _cy_i16[CY_THREADS - 1];
		t22	= _cy_r17[CY_THREADS - 1];	t23 = _cy_i17[CY_THREADS - 1];
		t24	= _cy_r18[CY_THREADS - 1];	t25 = _cy_i18[CY_THREADS - 1];
		t26	= _cy_r19[CY_THREADS - 1];	t27 = _cy_i19[CY_THREADS - 1];
		t28	= _cy_r20[CY_THREADS - 1];	t29 = _cy_i20[CY_THREADS - 1];
		t2A	= _cy_r21[CY_THREADS - 1];	t2B = _cy_i21[CY_THREADS - 1];
		t2C	= _cy_r22[CY_THREADS - 1];	t2D = _cy_i22[CY_THREADS - 1];
		t2E	= _cy_r23[CY_THREADS - 1];	t2F = _cy_i23[CY_THREADS - 1];
		t30	= _cy_r24[CY_THREADS - 1];	t31 = _cy_i24[CY_THREADS - 1];
		t32	= _cy_r25[CY_THREADS - 1];	t33 = _cy_i25[CY_THREADS - 1];
		t34	= _cy_r26[CY_THREADS - 1];	t35 = _cy_i26[CY_THREADS - 1];
		t36	= _cy_r27[CY_THREADS - 1];	t37 = _cy_i27[CY_THREADS - 1];
		t38	= _cy_r28[CY_THREADS - 1];	t39 = _cy_i28[CY_THREADS - 1];
		t3A	= _cy_r29[CY_THREADS - 1];	t3B = _cy_i29[CY_THREADS - 1];
		t3C	= _cy_r30[CY_THREADS - 1];	t3D = _cy_i30[CY_THREADS - 1];
		t3E	= _cy_r31[CY_THREADS - 1];	t3F = _cy_i31[CY_THREADS - 1];
		t40	= _cy_r32[CY_THREADS - 1];	t41 = _cy_i32[CY_THREADS - 1];
		t42	= _cy_r33[CY_THREADS - 1];	t43 = _cy_i33[CY_THREADS - 1];
		t44	= _cy_r34[CY_THREADS - 1];	t45 = _cy_i34[CY_THREADS - 1];
		t46	= _cy_r35[CY_THREADS - 1];	t47 = _cy_i35[CY_THREADS - 1];
		t48	= _cy_r36[CY_THREADS - 1];	t49 = _cy_i36[CY_THREADS - 1];
		t4A	= _cy_r37[CY_THREADS - 1];	t4B = _cy_i37[CY_THREADS - 1];
		t4C	= _cy_r38[CY_THREADS - 1];	t4D = _cy_i38[CY_THREADS - 1];
		t4E	= _cy_r39[CY_THREADS - 1];	t4F = _cy_i39[CY_THREADS - 1];
		t50	= _cy_r40[CY_THREADS - 1];	t51 = _cy_i40[CY_THREADS - 1];
		t52	= _cy_r41[CY_THREADS - 1];	t53 = _cy_i41[CY_THREADS - 1];
		t54	= _cy_r42[CY_THREADS - 1];	t55 = _cy_i42[CY_THREADS - 1];
		t56	= _cy_r43[CY_THREADS - 1];	t57 = _cy_i43[CY_THREADS - 1];
		t58	= _cy_r44[CY_THREADS - 1];	t59 = _cy_i44[CY_THREADS - 1];
		t5A	= _cy_r45[CY_THREADS - 1];	t5B = _cy_i45[CY_THREADS - 1];
		t5C	= _cy_r46[CY_THREADS - 1];	t5D = _cy_i46[CY_THREADS - 1];
		t5E	= _cy_r47[CY_THREADS - 1];	t5F = _cy_i47[CY_THREADS - 1];
		t60	= _cy_r48[CY_THREADS - 1];	t61 = _cy_i48[CY_THREADS - 1];
		t62	= _cy_r49[CY_THREADS - 1];	t63 = _cy_i49[CY_THREADS - 1];
		t64	= _cy_r50[CY_THREADS - 1];	t65 = _cy_i50[CY_THREADS - 1];
		t66	= _cy_r51[CY_THREADS - 1];	t67 = _cy_i51[CY_THREADS - 1];
		t68	= _cy_r52[CY_THREADS - 1];	t69 = _cy_i52[CY_THREADS - 1];
		t6A	= _cy_r53[CY_THREADS - 1];	t6B = _cy_i53[CY_THREADS - 1];
		t6C	= _cy_r54[CY_THREADS - 1];	t6D = _cy_i54[CY_THREADS - 1];
		t6E	= _cy_r55[CY_THREADS - 1];	t6F = _cy_i55[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1, "CY_THREADS must be > 1!");	/* Make sure loop only gets executed if multiple threads */
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
			_cy_r28[ithread] = _cy_r28[ithread-1];		_cy_i28[ithread] = _cy_i28[ithread-1];
			_cy_r29[ithread] = _cy_r29[ithread-1];		_cy_i29[ithread] = _cy_i29[ithread-1];
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
			_cy_r40[ithread] = _cy_r40[ithread-1];		_cy_i40[ithread] = _cy_i40[ithread-1];
			_cy_r41[ithread] = _cy_r41[ithread-1];		_cy_i41[ithread] = _cy_i41[ithread-1];
			_cy_r42[ithread] = _cy_r42[ithread-1];		_cy_i42[ithread] = _cy_i42[ithread-1];
			_cy_r43[ithread] = _cy_r43[ithread-1];		_cy_i43[ithread] = _cy_i43[ithread-1];
			_cy_r44[ithread] = _cy_r44[ithread-1];		_cy_i44[ithread] = _cy_i44[ithread-1];
			_cy_r45[ithread] = _cy_r45[ithread-1];		_cy_i45[ithread] = _cy_i45[ithread-1];
			_cy_r46[ithread] = _cy_r46[ithread-1];		_cy_i46[ithread] = _cy_i46[ithread-1];
			_cy_r47[ithread] = _cy_r47[ithread-1];		_cy_i47[ithread] = _cy_i47[ithread-1];
			_cy_r48[ithread] = _cy_r48[ithread-1];		_cy_i48[ithread] = _cy_i48[ithread-1];
			_cy_r49[ithread] = _cy_r49[ithread-1];		_cy_i49[ithread] = _cy_i49[ithread-1];
			_cy_r50[ithread] = _cy_r50[ithread-1];		_cy_i50[ithread] = _cy_i50[ithread-1];
			_cy_r51[ithread] = _cy_r51[ithread-1];		_cy_i51[ithread] = _cy_i51[ithread-1];
			_cy_r52[ithread] = _cy_r52[ithread-1];		_cy_i52[ithread] = _cy_i52[ithread-1];
			_cy_r53[ithread] = _cy_r53[ithread-1];		_cy_i53[ithread] = _cy_i53[ithread-1];
			_cy_r54[ithread] = _cy_r54[ithread-1];		_cy_i54[ithread] = _cy_i54[ithread-1];
			_cy_r55[ithread] = _cy_r55[ithread-1];		_cy_i55[ithread] = _cy_i55[ithread-1];
		}

		_cy_r00[0] =-t6F;	_cy_i00[0] =+t6E;	/* ...The 2 Mo"bius carries are here: */
		_cy_r01[0] = t00;	_cy_i01[0] = t01;
		_cy_r02[0] = t02;	_cy_i02[0] = t03;
		_cy_r03[0] = t04;	_cy_i03[0] = t05;
		_cy_r04[0] = t06;	_cy_i04[0] = t07;
		_cy_r05[0] = t08;	_cy_i05[0] = t09;
		_cy_r06[0] = t0A;	_cy_i06[0] = t0B;
		_cy_r07[0] = t0C;	_cy_i07[0] = t0D;
		_cy_r08[0] = t0E;	_cy_i08[0] = t0F;
		_cy_r09[0] = t10;	_cy_i09[0] = t11;
		_cy_r10[0] = t12;	_cy_i10[0] = t13;
		_cy_r11[0] = t14;	_cy_i11[0] = t15;
		_cy_r12[0] = t16;	_cy_i12[0] = t17;
		_cy_r13[0] = t18;	_cy_i13[0] = t19;
		_cy_r14[0] = t1A;	_cy_i14[0] = t1B;
		_cy_r15[0] = t1C;	_cy_i15[0] = t1D;
		_cy_r16[0] = t1E;	_cy_i16[0] = t1F;
		_cy_r17[0] = t20;	_cy_i17[0] = t21;
		_cy_r18[0] = t22;	_cy_i18[0] = t23;
		_cy_r19[0] = t24;	_cy_i19[0] = t25;
		_cy_r20[0] = t26;	_cy_i20[0] = t27;
		_cy_r21[0] = t28;	_cy_i21[0] = t29;
		_cy_r22[0] = t2A;	_cy_i22[0] = t2B;
		_cy_r23[0] = t2C;	_cy_i23[0] = t2D;
		_cy_r24[0] = t2E;	_cy_i24[0] = t2F;
		_cy_r25[0] = t30;	_cy_i25[0] = t31;
		_cy_r26[0] = t32;	_cy_i26[0] = t33;
		_cy_r27[0] = t34;	_cy_i27[0] = t35;
		_cy_r28[0] = t36;	_cy_i28[0] = t37;
		_cy_r29[0] = t38;	_cy_i29[0] = t39;
		_cy_r30[0] = t3A;	_cy_i30[0] = t3B;
		_cy_r31[0] = t3C;	_cy_i31[0] = t3D;
		_cy_r32[0] = t3E;	_cy_i32[0] = t3F;
		_cy_r33[0] = t40;	_cy_i33[0] = t41;
		_cy_r34[0] = t42;	_cy_i34[0] = t43;
		_cy_r35[0] = t44;	_cy_i35[0] = t45;
		_cy_r36[0] = t46;	_cy_i36[0] = t47;
		_cy_r37[0] = t48;	_cy_i37[0] = t49;
		_cy_r38[0] = t4A;	_cy_i38[0] = t4B;
		_cy_r39[0] = t4C;	_cy_i39[0] = t4D;
		_cy_r40[0] = t4E;	_cy_i40[0] = t4F;
		_cy_r41[0] = t50;	_cy_i41[0] = t51;
		_cy_r42[0] = t52;	_cy_i42[0] = t53;
		_cy_r43[0] = t54;	_cy_i43[0] = t55;
		_cy_r44[0] = t56;	_cy_i44[0] = t57;
		_cy_r45[0] = t58;	_cy_i45[0] = t59;
		_cy_r46[0] = t5A;	_cy_i46[0] = t5B;
		_cy_r47[0] = t5C;	_cy_i47[0] = t5D;
		_cy_r48[0] = t5E;	_cy_i48[0] = t5F;
		_cy_r49[0] = t60;	_cy_i49[0] = t61;
		_cy_r50[0] = t62;	_cy_i50[0] = t63;
		_cy_r51[0] = t64;	_cy_i51[0] = t65;
		_cy_r52[0] = t66;	_cy_i52[0] = t67;
		_cy_r53[0] = t68;	_cy_i53[0] = t69;
		_cy_r54[0] = t6A;	_cy_i54[0] = t6B;
		_cy_r55[0] = t6C;	_cy_i55[0] = t6D;
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
			jt = j;	// NB: pini *already* padded, so do not additionally pad index here!
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
		}
	}
}	/* endfor(outer) */

	t00 = 0;
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		t00 += fabs(_cy_r00[0])+fabs(_cy_i00[0])+fabs(_cy_r01[0])+fabs(_cy_i01[0])+fabs(_cy_r02[0])+fabs(_cy_i02[0])+fabs(_cy_r03[0])+fabs(_cy_i03[0])+fabs(_cy_r04[0])+fabs(_cy_i04[0])+fabs(_cy_r05[0])+fabs(_cy_i05[0])+fabs(_cy_r06[0])+fabs(_cy_i06[0])+fabs(_cy_r07[0])+fabs(_cy_i07[0])+fabs(_cy_r08[0])+fabs(_cy_i08[0])+fabs(_cy_r09[0])+fabs(_cy_i09[0]);
		t00 += fabs(_cy_r10[0])+fabs(_cy_i10[0])+fabs(_cy_r11[0])+fabs(_cy_i11[0])+fabs(_cy_r12[0])+fabs(_cy_i12[0])+fabs(_cy_r13[0])+fabs(_cy_i13[0])+fabs(_cy_r14[0])+fabs(_cy_i14[0])+fabs(_cy_r15[0])+fabs(_cy_i15[0])+fabs(_cy_r16[0])+fabs(_cy_i16[0])+fabs(_cy_r17[0])+fabs(_cy_i17[0])+fabs(_cy_r18[0])+fabs(_cy_i18[0])+fabs(_cy_r19[0])+fabs(_cy_i19[0]);
		t00 += fabs(_cy_r20[0])+fabs(_cy_i20[0])+fabs(_cy_r21[0])+fabs(_cy_i21[0])+fabs(_cy_r22[0])+fabs(_cy_i22[0])+fabs(_cy_r23[0])+fabs(_cy_i23[0])+fabs(_cy_r24[0])+fabs(_cy_i24[0])+fabs(_cy_r25[0])+fabs(_cy_i25[0])+fabs(_cy_r26[0])+fabs(_cy_i26[0])+fabs(_cy_r27[0])+fabs(_cy_i27[0])+fabs(_cy_r28[0])+fabs(_cy_i28[0])+fabs(_cy_r29[0])+fabs(_cy_i29[0]);
		t00 += fabs(_cy_r30[0])+fabs(_cy_i30[0])+fabs(_cy_r31[0])+fabs(_cy_i31[0])+fabs(_cy_r32[0])+fabs(_cy_i32[0])+fabs(_cy_r33[0])+fabs(_cy_i33[0])+fabs(_cy_r34[0])+fabs(_cy_i34[0])+fabs(_cy_r35[0])+fabs(_cy_i35[0])+fabs(_cy_r36[0])+fabs(_cy_i36[0])+fabs(_cy_r37[0])+fabs(_cy_i37[0])+fabs(_cy_r38[0])+fabs(_cy_i38[0])+fabs(_cy_r39[0])+fabs(_cy_i39[0]);
		t00 += fabs(_cy_r40[0])+fabs(_cy_i40[0])+fabs(_cy_r41[0])+fabs(_cy_i41[0])+fabs(_cy_r42[0])+fabs(_cy_i42[0])+fabs(_cy_r43[0])+fabs(_cy_i43[0])+fabs(_cy_r44[0])+fabs(_cy_i44[0])+fabs(_cy_r45[0])+fabs(_cy_i45[0])+fabs(_cy_r46[0])+fabs(_cy_i46[0])+fabs(_cy_r47[0])+fabs(_cy_i47[0])+fabs(_cy_r48[0])+fabs(_cy_i48[0])+fabs(_cy_r49[0])+fabs(_cy_i49[0]);
		t00 += fabs(_cy_r50[0])+fabs(_cy_i50[0])+fabs(_cy_r51[0])+fabs(_cy_i51[0])+fabs(_cy_r52[0])+fabs(_cy_i52[0])+fabs(_cy_r53[0])+fabs(_cy_i53[0])+fabs(_cy_r54[0])+fabs(_cy_i54[0])+fabs(_cy_r55[0])+fabs(_cy_i55[0]);

		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
	}

	if(t00 != 0.0)
	{
		sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in %s - input wordsize may be too small.\n",iter,func);
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

void radix56_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-56 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2,jt,jp;
	// p-indexing is hexadecimal here:
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30, first_entry=TRUE;
	static double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	double re,im,rt,it, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15
		,r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i
		,r08r,r08i,r09r,r09i,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i
		,r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i
		,r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i,r30r,r30i,r31r,r31i
		,r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i
		,r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i
		,r48r,r48i,r49r,r49i,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i;

	if(!first_entry && (n/56) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR = n/56;

/*   constant index offsets for array load/stores are here.	*/

		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p10 = p08 + p08;
		p18 = p10 + p08;
		p20 = p18 + p08;
		p28 = p20 + p08;
		p30 = p28 + p08;

		p01 += ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 += ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 += ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 += ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 += ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 += ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 += ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 += ( (p08 >> DAT_BITS) << PAD_BITS );
		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p18 += ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p28 += ( (p28 >> DAT_BITS) << PAD_BITS );
		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-56 pass is here.	*/

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
		j2 = j1+RE_IM_STRIDE;

	/*...gather the needed data (56 64-bit complex, i.e. 112 64-bit reals) and do 8 radix-7 transforms. */
	/*
	Twiddleless version arranges 8 sets of radix-7 DFT inputs as follows: 0 in upper left corner,
	decrement (mod 56) 8 horizontally and 7 vertically:

		RADIX_07_DFT(00,48,40,32,24,16,08)
		RADIX_07_DFT(49,41,33,25,17,09,01)
		RADIX_07_DFT(42,34,26,18,10,02,50)
		RADIX_07_DFT(35,27,19,11,03,51,43)
		RADIX_07_DFT(28,20,12,04,52,44,36)
		RADIX_07_DFT(21,13,05,53,45,37,29)
		RADIX_07_DFT(14,06,54,46,38,30,22)
		RADIX_07_DFT(07,55,47,39,31,23,15)
	*/
					 /*                                                                      inputs                                                           */ /*                   intermediates                   */ /*                                                  outputs                                    */ /*   sincos consts   */
		RADIX_07_DFT(a[j1    ],a[j2    ], a[j1+p30],a[j2+p30], a[j1+p28],a[j2+p28], a[j1+p20],a[j2+p20], a[j1+p18],a[j2+p18], a[j1+p10],a[j2+p10], a[j1+p08],a[j2+p08], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r00r,r00i,r08r,r08i,r16r,r16i,r24r,r24i,r32r,r32i,r40r,r40i,r48r,r48i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p01; jp = j2+p01;
		RADIX_07_DFT(a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r01r,r01i,r09r,r09i,r17r,r17i,r25r,r25i,r33r,r33i,r41r,r41i,r49r,r49i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p02; jp = j2+p02;
		RADIX_07_DFT(a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r02r,r02i,r10r,r10i,r18r,r18i,r26r,r26i,r34r,r34i,r42r,r42i,r50r,r50i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p03; jp = j2+p03;
		RADIX_07_DFT(a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r03r,r03i,r11r,r11i,r19r,r19i,r27r,r27i,r35r,r35i,r43r,r43i,r51r,r51i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p04; jp = j2+p04;
		RADIX_07_DFT(a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r04r,r04i,r12r,r12i,r20r,r20i,r28r,r28i,r36r,r36i,r44r,r44i,r52r,r52i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p05; jp = j2+p05;
		RADIX_07_DFT(a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r05r,r05i,r13r,r13i,r21r,r21i,r29r,r29i,r37r,r37i,r45r,r45i,r53r,r53i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p06; jp = j2+p06;
		RADIX_07_DFT(a[jt+p08],a[jp+p08], a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r06r,r06i,r14r,r14i,r22r,r22i,r30r,r30i,r38r,r38i,r46r,r46i,r54r,r54i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p07; jp = j2+p07;
		RADIX_07_DFT(a[jt    ],a[jp    ], a[jt+p30],a[jp+p30], a[jt+p28],a[jp+p28], a[jt+p20],a[jp+p20], a[jt+p18],a[jp+p18], a[jt+p10],a[jp+p10], a[jt+p08],a[jp+p08], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r07r,r07i,r15r,r15i,r23r,r23i,r31r,r31i,r39r,r39i,r47r,r47i,r55r,r55i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);

	/*...and now do 7 radix-8 transforms;
	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF]
	to properly permute the outputs of the ensuing 7 radix-8 DFTs to the required ordering, which in terms of our 7 DFTs is

		00,01,02,03,04,05,06,07		00,01,02,03,04,05,06,07 + p00[hex]
		51,50,48,49,55,54,52,53		03,02,00,01,07,06,04,05 + p30[hex]
		45,44,47,46,43,42,40,41	 =	05,04,07,06,03,02,00,01 + p28[hex]
		33,32,35,34,37,36,39,38		01,00,03,02,05,04,07,06 + p20[hex]
		30,31,29,28,25,24,27,26		06,07,05,04,01,00,03,02 + p18[hex]
		18,19,17,16,22,23,21,20		02,03,01,00,06,07,05,04 + p10[hex]
		12,13,14,15,10,11,09,08		04,05,06,07,02,03,01,00 + p08[hex]
	*/
					 /*                                                          inputs                                           */ /*                       intermediates                       */ /*                 outputs                   */
		RADIX_08_DIF(r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05],a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],rt,it);	jt = j1+p30; jp = j2+p30;
		RADIX_08_DIF(r08r,r08i,r09r,r09i,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],rt,it);	jt = j1+p28; jp = j2+p28;
		RADIX_08_DIF(r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	jt = j1+p20; jp = j2+p20;
		RADIX_08_DIF(r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i,r30r,r30i,r31r,r31i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],rt,it);	jt = j1+p18; jp = j2+p18;
		RADIX_08_DIF(r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p10; jp = j2+p10;
		RADIX_08_DIF(r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIF(r48r,r48i,r49r,r49i,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
	}
}

/***************/

void radix56_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-56 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2,jt,jp;
	// p-indexing is hexadecimal here:
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30, first_entry=TRUE;
	static double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	double re,im,rt,it, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15
		,r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i
		,r08r,r08i,r09r,r09i,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i
		,r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i
		,r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i,r30r,r30i,r31r,r31i
		,r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i
		,r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i
		,r48r,r48i,r49r,r49i,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i;

	if(!first_entry && (n/56) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR = n/56;

/*   constant index offsets for array load/stores are here.	*/

		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p10 = p08 + p08;
		p18 = p10 + p08;
		p20 = p18 + p08;
		p28 = p20 + p08;
		p30 = p28 + p08;

		p01 += ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 += ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 += ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 += ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 += ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 += ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 += ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 += ( (p08 >> DAT_BITS) << PAD_BITS );
		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p18 += ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p28 += ( (p28 >> DAT_BITS) << PAD_BITS );
		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-56 pass is here.	*/

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
		j2 = j1+RE_IM_STRIDE;
	/*
	Twiddleless version uses same linear-index-vector-form permutation as in DIF:

		00,01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,
		28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55
	map index-by-index [e.g. the index 01 gets replaced by 48, not the index at *position*...] to
		00,48,40,32,24,16,08,49,41,33,25,17,09,01,42,34,26,18,10,02,50,35,27,19,11,03,51,43,
		28,20,12,04,52,44,36,21,13,05,53,45,37,29,14,06,54,46,38,30,22,07,55,47,39,31,23,15.	(*)

	(***NOTE*** The following set of permutations can be auto-generated by running test_dft_radix()
	for the radix in question in TEST_TYPE = 1 [DIT] mode, skipping the actual DIT-pass step initially.)

	Remember, inputs to DIT are bit-reversed, so

		a[00,01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55] contain [cf. bit-reversal index array output by test_fft_radix()]
		x[00,28,14,42,07,35,21,49|01,29,15,43,08,36,22,50|02,30,16,44,09,37,23,51|03,31,17,45,10,38,24,52|04,32,18,46,11,39,25,53|05,33,19,47,12,40,26,54|06,34,20,48,13,41,27,55], which get swapped [using the permutation (*)] to
		x[00,28,42,14,49,21,35,07|48,20,34,06,41,13,27,55|40,12,26,54,33,05,19,47|32,04,18,46,25,53,11,39|24,52,10,38,17,45,03,31|16,44,02,30,09,37,51,23|08,36,50,22,01,29,43,15], which means the a-indices get swapped as
		a[00,01,03,02,07,06,05,04|51,50,49,48,53,52,54,55|45,44,46,47,41,40,42,43|33,32,34,35,38,39,36,37|30,31,28,29,26,27,24,25|18,19,16,17,20,21,23,22|12,13,15,14,08,09,11,10].

	These are the 7 octets going into the radix-8 DFTs:

		RADIX_08_DFT(00,01,03,02,07,06,05,04)
		RADIX_08_DFT(51,50,49,48,53,52,54,55)
		RADIX_08_DFT(45,44,46,47,41,40,42,43)
		RADIX_08_DFT(33,32,34,35,38,39,36,37)
		RADIX_08_DFT(30,31,28,29,26,27,24,25)
		RADIX_08_DFT(18,19,16,17,20,21,23,22)
		RADIX_08_DFT(12,13,15,14,08,09,11,10)

	Use the inverse supercalafragalistic Ancient Chinese Secret index-munging formula [iSACSIMPF]
	to properly permute the radix-5 DFT outputs into the required output ordering of the 8 radix-7 DFTs:

		00,08,16,24,32,40,48
		49,01,09,17,25,33,41
		42,50,02,10,18,26,34
		35,43,51,03,11,19,27
		28,36,44,52,04,12,20
		21,29,37,45,53,05,13
		14,22,30,38,46,54,06
		07,15,23,31,39,47,55
	*/
	/*...gather the needed data (56 64-bit complex, i.e. 112 64-bit reals) and do 7 radix-8 transforms...*/
					/*                                    inputs                                  */ /*                         outputs                   */
		RADIX_08_DIT(a[j1    ],a[j2    ], a[j1+p01],a[j2+p01], a[j1+p03],a[j2+p03], a[j1+p02],a[j2+p02], a[j1+p07],a[j2+p07], a[j1+p06],a[j2+p06], a[j1+p05],a[j2+p05], a[j1+p04],a[j2+p04], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i, rt,it);	jt = j1+p30; jp = j2+p30;
		RADIX_08_DIT(a[jt+p03],a[jp+p03], a[jt+p02],a[jp+p02], a[jt+p01],a[jp+p01], a[jt    ],a[jp    ], a[jt+p05],a[jp+p05], a[jt+p04],a[jp+p04], a[jt+p06],a[jp+p06], a[jt+p07],a[jp+p07], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r08r,r08i,r09r,r09i,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i, rt,it);	jt = j1+p28; jp = j2+p28;
		RADIX_08_DIT(a[jt+p05],a[jp+p05], a[jt+p04],a[jp+p04], a[jt+p06],a[jp+p06], a[jt+p07],a[jp+p07], a[jt+p01],a[jp+p01], a[jt    ],a[jp    ], a[jt+p02],a[jp+p02], a[jt+p03],a[jp+p03], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i, rt,it);	jt = j1+p20; jp = j2+p20;
		RADIX_08_DIT(a[jt+p01],a[jp+p01], a[jt    ],a[jp    ], a[jt+p02],a[jp+p02], a[jt+p03],a[jp+p03], a[jt+p06],a[jp+p06], a[jt+p07],a[jp+p07], a[jt+p04],a[jp+p04], a[jt+p05],a[jp+p05], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i,r30r,r30i,r31r,r31i, rt,it);	jt = j1+p18; jp = j2+p18;
		RADIX_08_DIT(a[jt+p06],a[jp+p06], a[jt+p07],a[jp+p07], a[jt+p04],a[jp+p04], a[jt+p05],a[jp+p05], a[jt+p02],a[jp+p02], a[jt+p03],a[jp+p03], a[jt    ],a[jp    ], a[jt+p01],a[jp+p01], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i, rt,it);	jt = j1+p10; jp = j2+p10;
		RADIX_08_DIT(a[jt+p02],a[jp+p02], a[jt+p03],a[jp+p03], a[jt    ],a[jp    ], a[jt+p01],a[jp+p01], a[jt+p04],a[jp+p04], a[jt+p05],a[jp+p05], a[jt+p07],a[jp+p07], a[jt+p06],a[jp+p06], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i, rt,it);	jt = j1+p08; jp = j2+p08;
		RADIX_08_DIT(a[jt+p04],a[jp+p04], a[jt+p05],a[jp+p05], a[jt+p07],a[jp+p07], a[jt+p06],a[jp+p06], a[jt    ],a[jp    ], a[jt+p01],a[jp+p01], a[jt+p03],a[jp+p03], a[jt+p02],a[jp+p02], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r48r,r48i,r49r,r49i,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i, rt,it);

	/*...and now do 8 radix-7 transforms, with the columns of a1p*[r,i] output pairs in the above 7x radix-8 set now acting as input rows: */
					/*                                              inputs                                          */ /*               intermediates              */ /*                                                                     outputs                                                           */
		RADIX_07_DFT(r00r,r00i, r08r,r08i, r16r,r16i, r24r,r24i, r32r,r32i, r40r,r40i, r48r,r48i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[j1    ],a[j2    ], a[j1+p08],a[j2+p08], a[j1+p10],a[j2+p10], a[j1+p18],a[j2+p18], a[j1+p20],a[j2+p20], a[j1+p28],a[j2+p28], a[j1+p30],a[j2+p30], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p01; jp = j2+p01;
		RADIX_07_DFT(r01r,r01i, r09r,r09i, r17r,r17i, r25r,r25i, r33r,r33i, r41r,r41i, r49r,r49i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p30],a[jp+p30], a[jt    ],a[jp    ], a[jt+p08],a[jp+p08], a[jt+p10],a[jp+p10], a[jt+p18],a[jp+p18], a[jt+p20],a[jp+p20], a[jt+p28],a[jp+p28], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p02; jp = j2+p02;
		RADIX_07_DFT(r02r,r02i, r10r,r10i, r18r,r18i, r26r,r26i, r34r,r34i, r42r,r42i, r50r,r50i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p28],a[jp+p28], a[jt+p30],a[jp+p30], a[jt    ],a[jp    ], a[jt+p08],a[jp+p08], a[jt+p10],a[jp+p10], a[jt+p18],a[jp+p18], a[jt+p20],a[jp+p20], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p03; jp = j2+p03;
		RADIX_07_DFT(r03r,r03i, r11r,r11i, r19r,r19i, r27r,r27i, r35r,r35i, r43r,r43i, r51r,r51i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p20],a[jp+p20], a[jt+p28],a[jp+p28], a[jt+p30],a[jp+p30], a[jt    ],a[jp    ], a[jt+p08],a[jp+p08], a[jt+p10],a[jp+p10], a[jt+p18],a[jp+p18], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p04; jp = j2+p04;
		RADIX_07_DFT(r04r,r04i, r12r,r12i, r20r,r20i, r28r,r28i, r36r,r36i, r44r,r44i, r52r,r52i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p18],a[jp+p18], a[jt+p20],a[jp+p20], a[jt+p28],a[jp+p28], a[jt+p30],a[jp+p30], a[jt    ],a[jp    ], a[jt+p08],a[jp+p08], a[jt+p10],a[jp+p10], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p05; jp = j2+p05;
		RADIX_07_DFT(r05r,r05i, r13r,r13i, r21r,r21i, r29r,r29i, r37r,r37i, r45r,r45i, r53r,r53i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p10],a[jp+p10], a[jt+p18],a[jp+p18], a[jt+p20],a[jp+p20], a[jt+p28],a[jp+p28], a[jt+p30],a[jp+p30], a[jt    ],a[jp    ], a[jt+p08],a[jp+p08], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p06; jp = j2+p06;
		RADIX_07_DFT(r06r,r06i, r14r,r14i, r22r,r22i, r30r,r30i, r38r,r38i, r46r,r46i, r54r,r54i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt+p08],a[jp+p08], a[jt+p10],a[jp+p10], a[jt+p18],a[jp+p18], a[jt+p20],a[jp+p20], a[jt+p28],a[jp+p28], a[jt+p30],a[jp+p30], a[jt    ],a[jp    ], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);	jt = j1+p07; jp = j2+p07;
		RADIX_07_DFT(r07r,r07i, r15r,r15i, r23r,r23i, r31r,r31i, r39r,r39i, r47r,r47i, r55r,r55i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[jt    ],a[jp    ], a[jt+p08],a[jp+p08], a[jt+p10],a[jp+p10], a[jt+p18],a[jp+p18], a[jt+p20],a[jp+p20], a[jt+p28],a[jp+p28], a[jt+p30],a[jp+p30], uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	}
}

/******************** Multithreaded function body: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void* 
	cy56_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
	const char func[] = "radix56_ditN_cy_dif1";
		const int RADIX = 56, odd_radix = 7;	// odd_radix = [radix >> trailz(radix)]
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		// lg(sizeof(vec_dbl)):
		uint32 p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,p20,p28,p30;
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

	  #ifdef USE_AVX
		const int l2_sz_vd = 5;
	  #else
		const int l2_sz_vd = 4;
	  #endif
		const double crnd = 3.0*0x4000000*0x2000000;
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7;	/* Addresses into array sections */
		vec_dbl *tmp,*tm2;	// Non-static utility ptrs
		vec_dbl *isrt2, *cc0, *ss0, *cc1, *ss1, *cc2, *ss2, *cc3, *ss3, *max_err, *sse2_rnd, *half_arr
		,*r00r,*r01r,*r02r,*r03r,*r04r,*r05r,*r06r,*r07r,*r08r,*r09r,*r10r,*r11r,*r12r,*r13r
		,*r14r,*r15r,*r16r,*r17r,*r18r,*r19r,*r20r,*r21r,*r22r,*r23r,*r24r,*r25r,*r26r,*r27r
		,*r28r,*r29r,*r30r,*r31r,*r32r,*r33r,*r34r,*r35r,*r36r,*r37r,*r38r,*r39r,*r40r,*r41r
		,*r42r,*r43r,*r44r,*r45r,*r46r,*r47r,*r48r,*r49r,*r50r,*r51r,*r52r,*r53r,*r54r,*r55r
		,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r
		,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r
		,*s1p28r,*s1p29r,*s1p30r,*s1p31r,*s1p32r,*s1p33r,*s1p34r,*s1p35r,*s1p36r,*s1p37r,*s1p38r,*s1p39r,*s1p40r,*s1p41r
		,*s1p42r,*s1p43r,*s1p44r,*s1p45r,*s1p46r,*s1p47r,*s1p48r,*s1p49r,*s1p50r,*s1p51r,*s1p52r,*s1p53r,*s1p54r,*s1p55r;
		int
		 *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13
		,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27
		,*bjmodn28,*bjmodn29,*bjmodn30,*bjmodn31,*bjmodn32,*bjmodn33,*bjmodn34,*bjmodn35,*bjmodn36,*bjmodn37,*bjmodn38,*bjmodn39,*bjmodn40,*bjmodn41
		,*bjmodn42,*bjmodn43,*bjmodn44,*bjmodn45,*bjmodn46,*bjmodn47,*bjmodn48,*bjmodn49,*bjmodn50,*bjmodn51,*bjmodn52,*bjmodn53,*bjmodn54,*bjmodn55;
		vec_dbl
		*cy_r00,*cy_r04,*cy_r08,*cy_r12,*cy_r16,*cy_r20,*cy_r24,*cy_r28,*cy_r32,*cy_r36,*cy_r40,*cy_r44,*cy_r48,*cy_r52,
		*cy_i00,*cy_i04,*cy_i08,*cy_i12,*cy_i16,*cy_i20,*cy_i24,*cy_i28,*cy_i32,*cy_i36,*cy_i40,*cy_i44,*cy_i48,*cy_i52;
	  #ifndef USE_AVX
		vec_dbl
		*cy_r02,*cy_r06,*cy_r10,*cy_r14,*cy_r18,*cy_r22,*cy_r26,*cy_r30,*cy_r34,*cy_r38,*cy_r42,*cy_r46,*cy_r50,*cy_r54,
		*cy_i02,*cy_i06,*cy_i10,*cy_i14,*cy_i18,*cy_i22,*cy_i26,*cy_i30,*cy_i34,*cy_i38,*cy_i42,*cy_i46,*cy_i50,*cy_i54;
	  #else
		vec_dbl *base_negacyclic_root;
	  #endif

		/* These are used in conjunction with the langth-7 arrays in the USE_SCALAR_CARRY #define below;
		In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
		int idx_offset, idx_incr;
		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else

		/* Non-SSE2 version assumes LO_ADD = 1 and uses the corresponding versions of the sincos constants: */
		const double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
						us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
						uc2 =-.22252093395631440426,	 /* cos(2u)	*/
						us2 = .97492791218182360702,	 /* sin(2u)	*/
						uc3 =-.90096886790241912622,	 /* cos(3u)	*/
						us3 = .43388373911755812050;	 /* sin(3u)	*/
		double *base, *baseinv, *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr, bs,bsinv;
		int wts_idx_incr;
		const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int jt,jp,m,m2,ntmp;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		double re,im,temp,frac,
			r00r,r01r,r02r,r03r,r04r,r05r,r06r,r07r,r08r,r09r,r10r,r11r,r12r,r13r,r14r,r15r,r16r,r17r,r18r,r19r,r20r,r21r,r22r,r23r,r24r,r25r,r26r,r27r,
			r28r,r29r,r30r,r31r,r32r,r33r,r34r,r35r,r36r,r37r,r38r,r39r,r40r,r41r,r42r,r43r,r44r,r45r,r46r,r47r,r48r,r49r,r50r,r51r,r52r,r53r,r54r,r55r,
			r00i,r01i,r02i,r03i,r04i,r05i,r06i,r07i,r08i,r09i,r10i,r11i,r12i,r13i,r14i,r15i,r16i,r17i,r18i,r19i,r20i,r21i,r22i,r23i,r24i,r25i,r26i,r27i,
			r28i,r29i,r30i,r31i,r32i,r33i,r34i,r35i,r36i,r37i,r38i,r39i,r40i,r41i,r42i,r43i,r44i,r45i,r46i,r47i,r48i,r49i,r50i,r51i,r52i,r53i,r54i,r55i,
			a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r,
			a1p28r,a1p29r,a1p30r,a1p31r,a1p32r,a1p33r,a1p34r,a1p35r,a1p36r,a1p37r,a1p38r,a1p39r,a1p40r,a1p41r,a1p42r,a1p43r,a1p44r,a1p45r,a1p46r,a1p47r,a1p48r,a1p49r,a1p50r,a1p51r,a1p52r,a1p53r,a1p54r,a1p55r,
			a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i,
			a1p28i,a1p29i,a1p30i,a1p31i,a1p32i,a1p33i,a1p34i,a1p35i,a1p36i,a1p37i,a1p38i,a1p39i,a1p40i,a1p41i,a1p42i,a1p43i,a1p44i,a1p45i,a1p46i,a1p47i,a1p48i,a1p49i,a1p50i,a1p51i,a1p52i,a1p53i,a1p54i,a1p55i,
			cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27,
			cy_r28,cy_r29,cy_r30,cy_r31,cy_r32,cy_r33,cy_r34,cy_r35,cy_r36,cy_r37,cy_r38,cy_r39,cy_r40,cy_r41,cy_r42,cy_r43,cy_r44,cy_r45,cy_r46,cy_r47,cy_r48,cy_r49,cy_r50,cy_r51,cy_r52,cy_r53,cy_r54,cy_r55,
			cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27,
			cy_i28,cy_i29,cy_i30,cy_i31,cy_i32,cy_i33,cy_i34,cy_i35,cy_i36,cy_i37,cy_i38,cy_i39,cy_i40,cy_i41,cy_i42,cy_i43,cy_i44,cy_i45,cy_i46,cy_i47,cy_i48,cy_i49,cy_i50,cy_i51,cy_i52,cy_i53,cy_i54,cy_i55,
			t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15;
		int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,
			bjmodn28,bjmodn29,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,bjmodn36,bjmodn37,bjmodn38,bjmodn39,bjmodn40,bjmodn41,bjmodn42,bjmodn43,bjmodn44,bjmodn45,bjmodn46,bjmodn47,bjmodn48,bjmodn49,bjmodn50,bjmodn51,bjmodn52,bjmodn53,bjmodn54,bjmodn55;

	#endif

		struct cy_thread_data_t* thread_arg = targ;
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
		double scale = thread_arg->scale;

	// pointer data:
		double *a = thread_arg->arrdat;
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
		int *si = thread_arg->si;
		struct complex *rn0 = thread_arg->rn0;
		struct complex *rn1 = thread_arg->rn1;

		/*   constant index offsets for array load/stores are here.	*/
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p10 = p08 + p08;
		p18 = p10 + p08;
		p20 = p18 + p08;
		p28 = p20 + p08;
		p30 = p28 + p08;

		p01 += ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 += ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 += ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 += ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 += ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 += ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 += ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 += ( (p08 >> DAT_BITS) << PAD_BITS );
		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p18 += ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p28 += ( (p28 >> DAT_BITS) << PAD_BITS );
		p30 += ( (p30 >> DAT_BITS) << PAD_BITS );

	#ifdef USE_SSE2
		uint32 nwt16 = nwt << l2_sz_vd;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code

		r00r = thread_arg->r00r;	// declared above
							tmp = r00r + 0x38;
		r00r = r00r + 0x00;	r28r = tmp + 0x00;
		r01r = r00r + 0x02;	r29r = tmp + 0x02;
		r02r = r00r + 0x04;	r30r = tmp + 0x04;
		r03r = r00r + 0x06;	r31r = tmp + 0x06;
		r04r = r00r + 0x08;	r32r = tmp + 0x08;
		r05r = r00r + 0x0a;	r33r = tmp + 0x0a;
		r06r = r00r + 0x0c;	r34r = tmp + 0x0c;
		r07r = r00r + 0x0e;	r35r = tmp + 0x0e;
		r08r = r00r + 0x10;	r36r = tmp + 0x10;
		r09r = r00r + 0x12;	r37r = tmp + 0x12;
		r10r = r00r + 0x14;	r38r = tmp + 0x14;
		r11r = r00r + 0x16;	r39r = tmp + 0x16;
		r12r = r00r + 0x18;	r40r = tmp + 0x18;
		r13r = r00r + 0x1a;	r41r = tmp + 0x1a;
		r14r = r00r + 0x1c;	r42r = tmp + 0x1c;
		r15r = r00r + 0x1e;	r43r = tmp + 0x1e;
		r16r = r00r + 0x20;	r44r = tmp + 0x20;
		r17r = r00r + 0x22;	r45r = tmp + 0x22;
		r18r = r00r + 0x24;	r46r = tmp + 0x24;
		r19r = r00r + 0x26;	r47r = tmp + 0x26;
		r20r = r00r + 0x28;	r48r = tmp + 0x28;
		r21r = r00r + 0x2a;	r49r = tmp + 0x2a;
		r22r = r00r + 0x2c;	r50r = tmp + 0x2c;
		r23r = r00r + 0x2e;	r51r = tmp + 0x2e;
		r24r = r00r + 0x30;	r52r = tmp + 0x30;
		r25r = r00r + 0x32;	r53r = tmp + 0x32;
		r26r = r00r + 0x34;	r54r = tmp + 0x34;
		r27r = r00r + 0x36;	r55r = tmp + 0x36;
		tmp += 0x38;	// sc_ptr += 0x70
		tm2 = tmp + 0x38;
		s1p00r = tmp + 0x00;	s1p28r = tm2 + 0x00;
		s1p01r = tmp + 0x02;	s1p29r = tm2 + 0x02;
		s1p02r = tmp + 0x04;	s1p30r = tm2 + 0x04;
		s1p03r = tmp + 0x06;	s1p31r = tm2 + 0x06;
		s1p04r = tmp + 0x08;	s1p32r = tm2 + 0x08;
		s1p05r = tmp + 0x0a;	s1p33r = tm2 + 0x0a;
		s1p06r = tmp + 0x0c;	s1p34r = tm2 + 0x0c;
		s1p07r = tmp + 0x0e;	s1p35r = tm2 + 0x0e;
		s1p08r = tmp + 0x10;	s1p36r = tm2 + 0x10;
		s1p09r = tmp + 0x12;	s1p37r = tm2 + 0x12;
		s1p10r = tmp + 0x14;	s1p38r = tm2 + 0x14;
		s1p11r = tmp + 0x16;	s1p39r = tm2 + 0x16;
		s1p12r = tmp + 0x18;	s1p40r = tm2 + 0x18;
		s1p13r = tmp + 0x1a;	s1p41r = tm2 + 0x1a;
		s1p14r = tmp + 0x1c;	s1p42r = tm2 + 0x1c;
		s1p15r = tmp + 0x1e;	s1p43r = tm2 + 0x1e;
		s1p16r = tmp + 0x20;	s1p44r = tm2 + 0x20;
		s1p17r = tmp + 0x22;	s1p45r = tm2 + 0x22;
		s1p18r = tmp + 0x24;	s1p46r = tm2 + 0x24;
		s1p19r = tmp + 0x26;	s1p47r = tm2 + 0x26;
		s1p20r = tmp + 0x28;	s1p48r = tm2 + 0x28;
		s1p21r = tmp + 0x2a;	s1p49r = tm2 + 0x2a;
		s1p22r = tmp + 0x2c;	s1p50r = tm2 + 0x2c;
		s1p23r = tmp + 0x2e;	s1p51r = tm2 + 0x2e;
		s1p24r = tmp + 0x30;	s1p52r = tm2 + 0x30;
		s1p25r = tmp + 0x32;	s1p53r = tm2 + 0x32;
		s1p26r = tmp + 0x34;	s1p54r = tm2 + 0x34;
		s1p27r = tmp + 0x36;	s1p55r = tm2 + 0x36;
		tmp += 0x71;	// sc_ptr += 0xe1 [2x +0x38, plut jump of 1 for isrt2]
		isrt2   = tmp - 0x01;
		cc0		= tmp + 0x00;
		ss0		= tmp + 0x01;
		cc1		= tmp + 0x02;
		ss1		= tmp + 0x03;
		cc2		= tmp + 0x04;
		ss2		= tmp + 0x05;
		cc3  	= tmp + 0x06;
		ss3		= tmp + 0x07;	/* Pad with extra 4 slots for scratch storage needed by SSE2_RADIX_07_DFT macro here, plus 1 to make offset even */
		tmp += 0x08 + 0x5;	// sc_ptr += 0xee
	  #ifdef USE_AVX
		cy_r00 = tmp + 0x00;
		cy_r04 = tmp + 0x01;
		cy_r08 = tmp + 0x02;
		cy_r12 = tmp + 0x03;
		cy_r16 = tmp + 0x04;
		cy_r20 = tmp + 0x05;
		cy_r24 = tmp + 0x06;
		cy_r28 = tmp + 0x07;
		cy_r32 = tmp + 0x08;
		cy_r36 = tmp + 0x09;
		cy_r40 = tmp + 0x0a;
		cy_r44 = tmp + 0x0b;
		cy_r48 = tmp + 0x0c;
		cy_r52 = tmp + 0x0d;
		cy_i00 = tmp + 0x0e;
		cy_i04 = tmp + 0x0f;
		cy_i08 = tmp + 0x10;
		cy_i12 = tmp + 0x11;
		cy_i16 = tmp + 0x12;
		cy_i20 = tmp + 0x13;
		cy_i24 = tmp + 0x14;
		cy_i28 = tmp + 0x15;
		cy_i32 = tmp + 0x16;
		cy_i36 = tmp + 0x17;
		cy_i40 = tmp + 0x18;
		cy_i44 = tmp + 0x19;
		cy_i48 = tmp + 0x1a;
		cy_i52 = tmp + 0x1b;
		tmp += 0x1c;	// sc_ptr += 0x10a
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x10c; This is where the value of half_arr_offset56 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */

		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod
	  #else
		cy_r00 = tmp + 0x00;	cy_r02 = tmp + 0x01;
		cy_r04 = tmp + 0x02;	cy_r06 = tmp + 0x03;
		cy_r08 = tmp + 0x04;	cy_r10 = tmp + 0x05;
		cy_r12 = tmp + 0x06;	cy_r14 = tmp + 0x07;
		cy_r16 = tmp + 0x08;	cy_r18 = tmp + 0x09;
		cy_r20 = tmp + 0x0a;	cy_r22 = tmp + 0x0b;
		cy_r24 = tmp + 0x0c;	cy_r26 = tmp + 0x0d;
		cy_r28 = tmp + 0x0e;	cy_r30 = tmp + 0x0f;
		cy_r32 = tmp + 0x10;	cy_r34 = tmp + 0x11;
		cy_r36 = tmp + 0x12;	cy_r38 = tmp + 0x13;
		cy_r40 = tmp + 0x14;	cy_r42 = tmp + 0x15;
		cy_r44 = tmp + 0x16;	cy_r46 = tmp + 0x17;
		cy_r48 = tmp + 0x18;	cy_r50 = tmp + 0x19;
		cy_r52 = tmp + 0x1a;	cy_r54 = tmp + 0x1b;
		tmp += 0x1c;
		cy_i00 = tmp + 0x00;	cy_i02 = tmp + 0x01;
		cy_i04 = tmp + 0x02;	cy_i06 = tmp + 0x03;
		cy_i08 = tmp + 0x04;	cy_i10 = tmp + 0x05;
		cy_i12 = tmp + 0x06;	cy_i14 = tmp + 0x07;
		cy_i16 = tmp + 0x08;	cy_i18 = tmp + 0x09;
		cy_i20 = tmp + 0x0a;	cy_i22 = tmp + 0x0b;
		cy_i24 = tmp + 0x0c;	cy_i26 = tmp + 0x0d;
		cy_i28 = tmp + 0x0e;	cy_i30 = tmp + 0x0f;
		cy_i32 = tmp + 0x10;	cy_i34 = tmp + 0x11;
		cy_i36 = tmp + 0x12;	cy_i38 = tmp + 0x13;
		cy_i40 = tmp + 0x14;	cy_i42 = tmp + 0x15;
		cy_i44 = tmp + 0x16;	cy_i46 = tmp + 0x17;
		cy_i48 = tmp + 0x18;	cy_i50 = tmp + 0x19;
		cy_i52 = tmp + 0x1a;	cy_i54 = tmp + 0x1b;
		tmp += 0x1c;	// sc_ptr += 0x126
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// sc_ptr += 0x128; This is where the value of half_arr_offset56 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*odd_radix] x 16 for Fermat-mod */
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
		dtmp = (tmp)->d0 * (tmp+odd_radix)->d0;	ASSERT(HERE, fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp)->d1 * (tmp+odd_radix)->d1;	ASSERT(HERE, fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
	}

		VEC_DBL_INIT(max_err, 0.0);

		sign_mask = (uint64*)(r00r + radix56_creals_in_local_store);
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
									bjmodn28 = bjmodn00 + 28;
		bjmodn01 = bjmodn00 +  1;	bjmodn29 = bjmodn00 + 29;
		bjmodn02 = bjmodn00 +  2;	bjmodn30 = bjmodn00 + 30;
		bjmodn03 = bjmodn00 +  3;	bjmodn31 = bjmodn00 + 31;
		bjmodn04 = bjmodn00 +  4;	bjmodn32 = bjmodn00 + 32;
		bjmodn05 = bjmodn00 +  5;	bjmodn33 = bjmodn00 + 33;
		bjmodn06 = bjmodn00 +  6;	bjmodn34 = bjmodn00 + 34;
		bjmodn07 = bjmodn00 +  7;	bjmodn35 = bjmodn00 + 35;
		bjmodn08 = bjmodn00 +  8;	bjmodn36 = bjmodn00 + 36;
		bjmodn09 = bjmodn00 +  9;	bjmodn37 = bjmodn00 + 37;
		bjmodn10 = bjmodn00 + 10;	bjmodn38 = bjmodn00 + 38;
		bjmodn11 = bjmodn00 + 11;	bjmodn39 = bjmodn00 + 39;
		bjmodn12 = bjmodn00 + 12;	bjmodn40 = bjmodn00 + 40;
		bjmodn13 = bjmodn00 + 13;	bjmodn41 = bjmodn00 + 41;
		bjmodn14 = bjmodn00 + 14;	bjmodn42 = bjmodn00 + 42;
		bjmodn15 = bjmodn00 + 15;	bjmodn43 = bjmodn00 + 43;
		bjmodn16 = bjmodn00 + 16;	bjmodn44 = bjmodn00 + 44;
		bjmodn17 = bjmodn00 + 17;	bjmodn45 = bjmodn00 + 45;
		bjmodn18 = bjmodn00 + 18;	bjmodn46 = bjmodn00 + 46;
		bjmodn19 = bjmodn00 + 19;	bjmodn47 = bjmodn00 + 47;
		bjmodn20 = bjmodn00 + 20;	bjmodn48 = bjmodn00 + 48;
		bjmodn21 = bjmodn00 + 21;	bjmodn49 = bjmodn00 + 49;
		bjmodn22 = bjmodn00 + 22;	bjmodn50 = bjmodn00 + 50;
		bjmodn23 = bjmodn00 + 23;	bjmodn51 = bjmodn00 + 51;
		bjmodn24 = bjmodn00 + 24;	bjmodn52 = bjmodn00 + 52;
		bjmodn25 = bjmodn00 + 25;	bjmodn53 = bjmodn00 + 53;
		bjmodn26 = bjmodn00 + 26;	bjmodn54 = bjmodn00 + 54;
		bjmodn27 = bjmodn00 + 27;	bjmodn55 = bjmodn00 + 55;
	#else

		// In scalar mode use these 2 ptrs to pass wts_idx_incr and the base/baseinv/etc array-ptrs:
		wts_idx_incr = *(int *)thread_arg->half_arr;
		base      = (double *)thread_arg->r00r;
		baseinv   = base + 2;
		wt_arr    = base + 4;
		wtinv_arr = wt_arr    + odd_radix;
		bs_arr    = wtinv_arr + odd_radix;
		bsinv_arr = bs_arr    + odd_radix;

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
			*bjmodn28 = thread_arg->bjmodn28;			cy_r28->d0 = thread_arg->cy_r28;
			*bjmodn29 = thread_arg->bjmodn29;			cy_r28->d1 = thread_arg->cy_r29;
			*bjmodn30 = thread_arg->bjmodn30;			cy_r28->d2 = thread_arg->cy_r30;
			*bjmodn31 = thread_arg->bjmodn31;			cy_r28->d3 = thread_arg->cy_r31;
			*bjmodn32 = thread_arg->bjmodn32;			cy_r32->d0 = thread_arg->cy_r32;
			*bjmodn33 = thread_arg->bjmodn33;			cy_r32->d1 = thread_arg->cy_r33;
			*bjmodn34 = thread_arg->bjmodn34;			cy_r32->d2 = thread_arg->cy_r34;
			*bjmodn35 = thread_arg->bjmodn35;			cy_r32->d3 = thread_arg->cy_r35;
			*bjmodn36 = thread_arg->bjmodn36;			cy_r36->d0 = thread_arg->cy_r36;
			*bjmodn37 = thread_arg->bjmodn37;			cy_r36->d1 = thread_arg->cy_r37;
			*bjmodn38 = thread_arg->bjmodn38;			cy_r36->d2 = thread_arg->cy_r38;
			*bjmodn39 = thread_arg->bjmodn39;			cy_r36->d3 = thread_arg->cy_r39;
			*bjmodn40 = thread_arg->bjmodn40;			cy_r40->d0 = thread_arg->cy_r40;
			*bjmodn41 = thread_arg->bjmodn41;			cy_r40->d1 = thread_arg->cy_r41;
			*bjmodn42 = thread_arg->bjmodn42;			cy_r40->d2 = thread_arg->cy_r42;
			*bjmodn43 = thread_arg->bjmodn43;			cy_r40->d3 = thread_arg->cy_r43;
			*bjmodn44 = thread_arg->bjmodn44;			cy_r44->d0 = thread_arg->cy_r44;
			*bjmodn45 = thread_arg->bjmodn45;			cy_r44->d1 = thread_arg->cy_r45;
			*bjmodn46 = thread_arg->bjmodn46;			cy_r44->d2 = thread_arg->cy_r46;
			*bjmodn47 = thread_arg->bjmodn47;			cy_r44->d3 = thread_arg->cy_r47;
			*bjmodn48 = thread_arg->bjmodn48;			cy_r48->d0 = thread_arg->cy_r48;
			*bjmodn49 = thread_arg->bjmodn49;			cy_r48->d1 = thread_arg->cy_r49;
			*bjmodn50 = thread_arg->bjmodn50;			cy_r48->d2 = thread_arg->cy_r50;
			*bjmodn51 = thread_arg->bjmodn51;			cy_r48->d3 = thread_arg->cy_r51;
			*bjmodn52 = thread_arg->bjmodn52;			cy_r52->d0 = thread_arg->cy_r52;
			*bjmodn53 = thread_arg->bjmodn53;			cy_r52->d1 = thread_arg->cy_r53;
			*bjmodn54 = thread_arg->bjmodn54;			cy_r52->d2 = thread_arg->cy_r54;
			*bjmodn55 = thread_arg->bjmodn55;			cy_r52->d3 = thread_arg->cy_r55;

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
			*bjmodn28 = thread_arg->bjmodn28;			cy_r28->d0 = thread_arg->cy_r28;
			*bjmodn29 = thread_arg->bjmodn29;			cy_r28->d1 = thread_arg->cy_r29;
			*bjmodn30 = thread_arg->bjmodn30;			cy_r30->d0 = thread_arg->cy_r30;
			*bjmodn31 = thread_arg->bjmodn31;			cy_r30->d1 = thread_arg->cy_r31;
			*bjmodn32 = thread_arg->bjmodn32;			cy_r32->d0 = thread_arg->cy_r32;
			*bjmodn33 = thread_arg->bjmodn33;			cy_r32->d1 = thread_arg->cy_r33;
			*bjmodn34 = thread_arg->bjmodn34;			cy_r34->d0 = thread_arg->cy_r34;
			*bjmodn35 = thread_arg->bjmodn35;			cy_r34->d1 = thread_arg->cy_r35;
			*bjmodn36 = thread_arg->bjmodn36;			cy_r36->d0 = thread_arg->cy_r36;
			*bjmodn37 = thread_arg->bjmodn37;			cy_r36->d1 = thread_arg->cy_r37;
			*bjmodn38 = thread_arg->bjmodn38;			cy_r38->d0 = thread_arg->cy_r38;
			*bjmodn39 = thread_arg->bjmodn39;			cy_r38->d1 = thread_arg->cy_r39;
			*bjmodn40 = thread_arg->bjmodn40;			cy_r40->d0 = thread_arg->cy_r40;
			*bjmodn41 = thread_arg->bjmodn41;			cy_r40->d1 = thread_arg->cy_r41;
			*bjmodn42 = thread_arg->bjmodn42;			cy_r42->d0 = thread_arg->cy_r42;
			*bjmodn43 = thread_arg->bjmodn43;			cy_r42->d1 = thread_arg->cy_r43;
			*bjmodn44 = thread_arg->bjmodn44;			cy_r44->d0 = thread_arg->cy_r44;
			*bjmodn45 = thread_arg->bjmodn45;			cy_r44->d1 = thread_arg->cy_r45;
			*bjmodn46 = thread_arg->bjmodn46;			cy_r46->d0 = thread_arg->cy_r46;
			*bjmodn47 = thread_arg->bjmodn47;			cy_r46->d1 = thread_arg->cy_r47;
			*bjmodn48 = thread_arg->bjmodn48;			cy_r48->d0 = thread_arg->cy_r48;
			*bjmodn49 = thread_arg->bjmodn49;			cy_r48->d1 = thread_arg->cy_r49;
			*bjmodn50 = thread_arg->bjmodn50;			cy_r50->d0 = thread_arg->cy_r50;
			*bjmodn51 = thread_arg->bjmodn51;			cy_r50->d1 = thread_arg->cy_r51;
			*bjmodn52 = thread_arg->bjmodn52;			cy_r52->d0 = thread_arg->cy_r52;
			*bjmodn53 = thread_arg->bjmodn53;			cy_r52->d1 = thread_arg->cy_r53;
			*bjmodn54 = thread_arg->bjmodn54;			cy_r54->d0 = thread_arg->cy_r54;
			*bjmodn55 = thread_arg->bjmodn55;			cy_r54->d1 = thread_arg->cy_r55;

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
			bjmodn28 = thread_arg->bjmodn28;			cy_r28 = thread_arg->cy_r28;
			bjmodn29 = thread_arg->bjmodn29;			cy_r29 = thread_arg->cy_r29;
			bjmodn30 = thread_arg->bjmodn30;			cy_r30 = thread_arg->cy_r30;
			bjmodn31 = thread_arg->bjmodn31;			cy_r31 = thread_arg->cy_r31;
			bjmodn32 = thread_arg->bjmodn32;			cy_r32 = thread_arg->cy_r32;
			bjmodn33 = thread_arg->bjmodn33;			cy_r33 = thread_arg->cy_r33;
			bjmodn34 = thread_arg->bjmodn34;			cy_r34 = thread_arg->cy_r34;
			bjmodn35 = thread_arg->bjmodn35;			cy_r35 = thread_arg->cy_r35;
			bjmodn36 = thread_arg->bjmodn36;			cy_r36 = thread_arg->cy_r36;
			bjmodn37 = thread_arg->bjmodn37;			cy_r37 = thread_arg->cy_r37;
			bjmodn38 = thread_arg->bjmodn38;			cy_r38 = thread_arg->cy_r38;
			bjmodn39 = thread_arg->bjmodn39;			cy_r39 = thread_arg->cy_r39;
			bjmodn40 = thread_arg->bjmodn40;			cy_r40 = thread_arg->cy_r40;
			bjmodn41 = thread_arg->bjmodn41;			cy_r41 = thread_arg->cy_r41;
			bjmodn42 = thread_arg->bjmodn42;			cy_r42 = thread_arg->cy_r42;
			bjmodn43 = thread_arg->bjmodn43;			cy_r43 = thread_arg->cy_r43;
			bjmodn44 = thread_arg->bjmodn44;			cy_r44 = thread_arg->cy_r44;
			bjmodn45 = thread_arg->bjmodn45;			cy_r45 = thread_arg->cy_r45;
			bjmodn46 = thread_arg->bjmodn46;			cy_r46 = thread_arg->cy_r46;
			bjmodn47 = thread_arg->bjmodn47;			cy_r47 = thread_arg->cy_r47;
			bjmodn48 = thread_arg->bjmodn48;			cy_r48 = thread_arg->cy_r48;
			bjmodn49 = thread_arg->bjmodn49;			cy_r49 = thread_arg->cy_r49;
			bjmodn50 = thread_arg->bjmodn50;			cy_r50 = thread_arg->cy_r50;
			bjmodn51 = thread_arg->bjmodn51;			cy_r51 = thread_arg->cy_r51;
			bjmodn52 = thread_arg->bjmodn52;			cy_r52 = thread_arg->cy_r52;
			bjmodn53 = thread_arg->bjmodn53;			cy_r53 = thread_arg->cy_r53;
			bjmodn54 = thread_arg->bjmodn54;			cy_r54 = thread_arg->cy_r54;
			bjmodn55 = thread_arg->bjmodn55;			cy_r55 = thread_arg->cy_r55;

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
			cy_r28->d0 = thread_arg->cy_r28;	cy_i28->d0 = thread_arg->cy_i28;
			cy_r28->d1 = thread_arg->cy_r29;	cy_i28->d1 = thread_arg->cy_i29;
			cy_r28->d2 = thread_arg->cy_r30;	cy_i28->d2 = thread_arg->cy_i30;
			cy_r28->d3 = thread_arg->cy_r31;	cy_i28->d3 = thread_arg->cy_i31;
			cy_r32->d0 = thread_arg->cy_r32;	cy_i32->d0 = thread_arg->cy_i32;
			cy_r32->d1 = thread_arg->cy_r33;	cy_i32->d1 = thread_arg->cy_i33;
			cy_r32->d2 = thread_arg->cy_r34;	cy_i32->d2 = thread_arg->cy_i34;
			cy_r32->d3 = thread_arg->cy_r35;	cy_i32->d3 = thread_arg->cy_i35;
			cy_r36->d0 = thread_arg->cy_r36;	cy_i36->d0 = thread_arg->cy_i36;
			cy_r36->d1 = thread_arg->cy_r37;	cy_i36->d1 = thread_arg->cy_i37;
			cy_r36->d2 = thread_arg->cy_r38;	cy_i36->d2 = thread_arg->cy_i38;
			cy_r36->d3 = thread_arg->cy_r39;	cy_i36->d3 = thread_arg->cy_i39;
			cy_r40->d0 = thread_arg->cy_r40;	cy_i40->d0 = thread_arg->cy_i40;
			cy_r40->d1 = thread_arg->cy_r41;	cy_i40->d1 = thread_arg->cy_i41;
			cy_r40->d2 = thread_arg->cy_r42;	cy_i40->d2 = thread_arg->cy_i42;
			cy_r40->d3 = thread_arg->cy_r43;	cy_i40->d3 = thread_arg->cy_i43;
			cy_r44->d0 = thread_arg->cy_r44;	cy_i44->d0 = thread_arg->cy_i44;
			cy_r44->d1 = thread_arg->cy_r45;	cy_i44->d1 = thread_arg->cy_i45;
			cy_r44->d2 = thread_arg->cy_r46;	cy_i44->d2 = thread_arg->cy_i46;
			cy_r44->d3 = thread_arg->cy_r47;	cy_i44->d3 = thread_arg->cy_i47;
			cy_r48->d0 = thread_arg->cy_r48;	cy_i48->d0 = thread_arg->cy_i48;
			cy_r48->d1 = thread_arg->cy_r49;	cy_i48->d1 = thread_arg->cy_i49;
			cy_r48->d2 = thread_arg->cy_r50;	cy_i48->d2 = thread_arg->cy_i50;
			cy_r48->d3 = thread_arg->cy_r51;	cy_i48->d3 = thread_arg->cy_i51;
			cy_r52->d0 = thread_arg->cy_r52;	cy_i52->d0 = thread_arg->cy_i52;
			cy_r52->d1 = thread_arg->cy_r53;	cy_i52->d1 = thread_arg->cy_i53;
			cy_r52->d2 = thread_arg->cy_r54;	cy_i52->d2 = thread_arg->cy_i54;
			cy_r52->d3 = thread_arg->cy_r55;	cy_i52->d3 = thread_arg->cy_i55;

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
			cy_r28->d0 = thread_arg->cy_r14;	cy_r28->d1 = thread_arg->cy_i14;
			cy_r30->d0 = thread_arg->cy_r15;	cy_r30->d1 = thread_arg->cy_i15;
			cy_r32->d0 = thread_arg->cy_r16;	cy_r32->d1 = thread_arg->cy_i16;
			cy_r34->d0 = thread_arg->cy_r17;	cy_r34->d1 = thread_arg->cy_i17;
			cy_r36->d0 = thread_arg->cy_r18;	cy_r36->d1 = thread_arg->cy_i18;
			cy_r38->d0 = thread_arg->cy_r19;	cy_r38->d1 = thread_arg->cy_i19;
			cy_r40->d0 = thread_arg->cy_r20;	cy_r40->d1 = thread_arg->cy_i20;
			cy_r42->d0 = thread_arg->cy_r21;	cy_r42->d1 = thread_arg->cy_i21;
			cy_r44->d0 = thread_arg->cy_r22;	cy_r44->d1 = thread_arg->cy_i22;
			cy_r46->d0 = thread_arg->cy_r23;	cy_r46->d1 = thread_arg->cy_i23;
			cy_r48->d0 = thread_arg->cy_r24;	cy_r48->d1 = thread_arg->cy_i24;
			cy_r50->d0 = thread_arg->cy_r25;	cy_r50->d1 = thread_arg->cy_i25;
			cy_r52->d0 = thread_arg->cy_r26;	cy_r52->d1 = thread_arg->cy_i26;
			cy_r54->d0 = thread_arg->cy_r27;	cy_r54->d1 = thread_arg->cy_i27;
			cy_i00->d0 = thread_arg->cy_r28;	cy_i00->d1 = thread_arg->cy_i28;
			cy_i02->d0 = thread_arg->cy_r29;	cy_i02->d1 = thread_arg->cy_i29;
			cy_i04->d0 = thread_arg->cy_r30;	cy_i04->d1 = thread_arg->cy_i30;
			cy_i06->d0 = thread_arg->cy_r31;	cy_i06->d1 = thread_arg->cy_i31;
			cy_i08->d0 = thread_arg->cy_r32;	cy_i08->d1 = thread_arg->cy_i32;
			cy_i10->d0 = thread_arg->cy_r33;	cy_i10->d1 = thread_arg->cy_i33;
			cy_i12->d0 = thread_arg->cy_r34;	cy_i12->d1 = thread_arg->cy_i34;
			cy_i14->d0 = thread_arg->cy_r35;	cy_i14->d1 = thread_arg->cy_i35;
			cy_i16->d0 = thread_arg->cy_r36;	cy_i16->d1 = thread_arg->cy_i36;
			cy_i18->d0 = thread_arg->cy_r37;	cy_i18->d1 = thread_arg->cy_i37;
			cy_i20->d0 = thread_arg->cy_r38;	cy_i20->d1 = thread_arg->cy_i38;
			cy_i22->d0 = thread_arg->cy_r39;	cy_i22->d1 = thread_arg->cy_i39;
			cy_i24->d0 = thread_arg->cy_r40;	cy_i24->d1 = thread_arg->cy_i40;
			cy_i26->d0 = thread_arg->cy_r41;	cy_i26->d1 = thread_arg->cy_i41;
			cy_i28->d0 = thread_arg->cy_r42;	cy_i28->d1 = thread_arg->cy_i42;
			cy_i30->d0 = thread_arg->cy_r43;	cy_i30->d1 = thread_arg->cy_i43;
			cy_i32->d0 = thread_arg->cy_r44;	cy_i32->d1 = thread_arg->cy_i44;
			cy_i34->d0 = thread_arg->cy_r45;	cy_i34->d1 = thread_arg->cy_i45;
			cy_i36->d0 = thread_arg->cy_r46;	cy_i36->d1 = thread_arg->cy_i46;
			cy_i38->d0 = thread_arg->cy_r47;	cy_i38->d1 = thread_arg->cy_i47;
			cy_i40->d0 = thread_arg->cy_r48;	cy_i40->d1 = thread_arg->cy_i48;
			cy_i42->d0 = thread_arg->cy_r49;	cy_i42->d1 = thread_arg->cy_i49;
			cy_i44->d0 = thread_arg->cy_r50;	cy_i44->d1 = thread_arg->cy_i50;
			cy_i46->d0 = thread_arg->cy_r51;	cy_i46->d1 = thread_arg->cy_i51;
			cy_i48->d0 = thread_arg->cy_r52;	cy_i48->d1 = thread_arg->cy_i52;
			cy_i50->d0 = thread_arg->cy_r53;	cy_i50->d1 = thread_arg->cy_i53;
			cy_i52->d0 = thread_arg->cy_r54;	cy_i52->d1 = thread_arg->cy_i54;
			cy_i54->d0 = thread_arg->cy_r55;	cy_i54->d1 = thread_arg->cy_i55;

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
			cy_r28 = thread_arg->cy_r28;		cy_i28 = thread_arg->cy_i28;
			cy_r29 = thread_arg->cy_r29;		cy_i29 = thread_arg->cy_i29;
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
			cy_r40 = thread_arg->cy_r40;		cy_i40 = thread_arg->cy_i40;
			cy_r41 = thread_arg->cy_r41;		cy_i41 = thread_arg->cy_i41;
			cy_r42 = thread_arg->cy_r42;		cy_i42 = thread_arg->cy_i42;
			cy_r43 = thread_arg->cy_r43;		cy_i43 = thread_arg->cy_i43;
			cy_r44 = thread_arg->cy_r44;		cy_i44 = thread_arg->cy_i44;
			cy_r45 = thread_arg->cy_r45;		cy_i45 = thread_arg->cy_i45;
			cy_r46 = thread_arg->cy_r46;		cy_i46 = thread_arg->cy_i46;
			cy_r47 = thread_arg->cy_r47;		cy_i47 = thread_arg->cy_i47;
			cy_r48 = thread_arg->cy_r48;		cy_i48 = thread_arg->cy_i48;
			cy_r49 = thread_arg->cy_r49;		cy_i49 = thread_arg->cy_i49;
			cy_r50 = thread_arg->cy_r50;		cy_i50 = thread_arg->cy_i50;
			cy_r51 = thread_arg->cy_r51;		cy_i51 = thread_arg->cy_i51;
			cy_r52 = thread_arg->cy_r52;		cy_i52 = thread_arg->cy_i52;
			cy_r53 = thread_arg->cy_r53;		cy_i53 = thread_arg->cy_i53;
			cy_r54 = thread_arg->cy_r54;		cy_i54 = thread_arg->cy_i54;
			cy_r55 = thread_arg->cy_r55;		cy_i55 = thread_arg->cy_i55;

		#endif
		}

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
			for(j = jstart; j < jhi; j += stride)	// Stride = 4 reals for SSE2, 8 for AVX
			{
				j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1+RE_IM_STRIDE;

			/*...The radix-56 DIT pass is here:	*/

			#ifdef USE_SSE2

				/* Outputs of radix-4 are in consecutive memory locs, i.e. 0x20 bytes apart, e.g. the distance between r00r and r01r: */

				add0 = &a[j1    ];
				add1 = add0+p01;
				add2 = add0+p03;
				add3 = add0+p02;
				add4 = add0+p07;
				add5 = add0+p06;
				add6 = add0+p05;
				add7 = add0+p04;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r00r, isrt2)

				add3 = &a[j1+p30];
				add0 = add3+p03;
				add1 = add3+p02;
				add2 = add3+p01;
				add4 = add3+p05;
				add5 = add3+p04;
				add6 = add3+p06;
				add7 = add3+p07;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r08r, isrt2)

				add5 = &a[j1+p28];
				add0 = add5+p05;
				add1 = add5+p04;
				add2 = add5+p06;
				add3 = add5+p07;
				add4 = add5+p01;
				add6 = add5+p02;
				add7 = add5+p03;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r16r, isrt2)

				add1 = &a[j1+p20];
				add0 = add1+p01;
				add2 = add1+p02;
				add3 = add1+p03;
				add4 = add1+p06;
				add5 = add1+p07;
				add6 = add1+p04;
				add7 = add1+p05;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r24r, isrt2)

				add6 = &a[j1+p18];
				add0 = add6+p06;
				add1 = add6+p07;
				add2 = add6+p04;
				add3 = add6+p05;
				add4 = add6+p02;
				add5 = add6+p03;
				add7 = add6+p01;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r32r, isrt2)

				add2 = &a[j1+p10];
				add0 = add2+p02;
				add1 = add2+p03;
				add3 = add2+p01;
				add4 = add2+p04;
				add5 = add2+p05;
				add6 = add2+p07;
				add7 = add2+p06;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r40r, isrt2)

				add4 = &a[j1+p08];
				add0 = add4+p04;
				add1 = add4+p05;
				add2 = add4+p07;
				add3 = add4+p06;
				add5 = add4+p01;
				add6 = add4+p03;
				add7 = add4+p02;
				SSE2_RADIX8_DIT_0TWIDDLE(add0,add1,add2,add3,add4,add5,add6,add7, r48r, isrt2)

			/*...and now do 8 radix-7 transforms, with the columns of r*[r,i] output pairs in the above 7x radix-8 set now acting as input rows: */
								/*            inputs        */  /* sincos ptr */ /*            outputs                   */
				SSE2_RADIX_07_DFT(r00r,r08r,r16r,r24r,r32r,r40r,r48r, cc0, s1p00r,s1p08r,s1p16r,s1p24r,s1p32r,s1p40r,s1p48r);
				SSE2_RADIX_07_DFT(r01r,r09r,r17r,r25r,r33r,r41r,r49r, cc0, s1p49r,s1p01r,s1p09r,s1p17r,s1p25r,s1p33r,s1p41r);
				SSE2_RADIX_07_DFT(r02r,r10r,r18r,r26r,r34r,r42r,r50r, cc0, s1p42r,s1p50r,s1p02r,s1p10r,s1p18r,s1p26r,s1p34r);
				SSE2_RADIX_07_DFT(r03r,r11r,r19r,r27r,r35r,r43r,r51r, cc0, s1p35r,s1p43r,s1p51r,s1p03r,s1p11r,s1p19r,s1p27r);
				SSE2_RADIX_07_DFT(r04r,r12r,r20r,r28r,r36r,r44r,r52r, cc0, s1p28r,s1p36r,s1p44r,s1p52r,s1p04r,s1p12r,s1p20r);
				SSE2_RADIX_07_DFT(r05r,r13r,r21r,r29r,r37r,r45r,r53r, cc0, s1p21r,s1p29r,s1p37r,s1p45r,s1p53r,s1p05r,s1p13r);
				SSE2_RADIX_07_DFT(r06r,r14r,r22r,r30r,r38r,r46r,r54r, cc0, s1p14r,s1p22r,s1p30r,s1p38r,s1p46r,s1p54r,s1p06r);
				SSE2_RADIX_07_DFT(r07r,r15r,r23r,r31r,r39r,r47r,r55r, cc0, s1p07r,s1p15r,s1p23r,s1p31r,s1p39r,s1p47r,s1p55r);

			#else	// USE_SSE2 = False:

			/*...gather the needed data (56 64-bit complex, i.e. 112 64-bit reals) and do 7 radix-8 transforms...*/
							/*                                    inputs                                  */ /*                         outputs                   */
				RADIX_08_DIT(a[j1    ],a[j2    ], a[j1+p01],a[j2+p01], a[j1+p03],a[j2+p03], a[j1+p02],a[j2+p02], a[j1+p07],a[j2+p07], a[j1+p06],a[j2+p06], a[j1+p05],a[j2+p05], a[j1+p04],a[j2+p04], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i, rt,it);	jt = j1+p30; jp = j2+p30;
				RADIX_08_DIT(a[jt+p03],a[jp+p03], a[jt+p02],a[jp+p02], a[jt+p01],a[jp+p01], a[jt    ],a[jp    ], a[jt+p05],a[jp+p05], a[jt+p04],a[jp+p04], a[jt+p06],a[jp+p06], a[jt+p07],a[jp+p07], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r08r,r08i,r09r,r09i,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i, rt,it);	jt = j1+p28; jp = j2+p28;
				RADIX_08_DIT(a[jt+p05],a[jp+p05], a[jt+p04],a[jp+p04], a[jt+p06],a[jp+p06], a[jt+p07],a[jp+p07], a[jt+p01],a[jp+p01], a[jt    ],a[jp    ], a[jt+p02],a[jp+p02], a[jt+p03],a[jp+p03], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i, rt,it);	jt = j1+p20; jp = j2+p20;
				RADIX_08_DIT(a[jt+p01],a[jp+p01], a[jt    ],a[jp    ], a[jt+p02],a[jp+p02], a[jt+p03],a[jp+p03], a[jt+p06],a[jp+p06], a[jt+p07],a[jp+p07], a[jt+p04],a[jp+p04], a[jt+p05],a[jp+p05], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i,r30r,r30i,r31r,r31i, rt,it);	jt = j1+p18; jp = j2+p18;
				RADIX_08_DIT(a[jt+p06],a[jp+p06], a[jt+p07],a[jp+p07], a[jt+p04],a[jp+p04], a[jt+p05],a[jp+p05], a[jt+p02],a[jp+p02], a[jt+p03],a[jp+p03], a[jt    ],a[jp    ], a[jt+p01],a[jp+p01], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i, rt,it);	jt = j1+p10; jp = j2+p10;
				RADIX_08_DIT(a[jt+p02],a[jp+p02], a[jt+p03],a[jp+p03], a[jt    ],a[jp    ], a[jt+p01],a[jp+p01], a[jt+p04],a[jp+p04], a[jt+p05],a[jp+p05], a[jt+p07],a[jp+p07], a[jt+p06],a[jp+p06], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i, rt,it);	jt = j1+p08; jp = j2+p08;
				RADIX_08_DIT(a[jt+p04],a[jp+p04], a[jt+p05],a[jp+p05], a[jt+p07],a[jp+p07], a[jt+p06],a[jp+p06], a[jt    ],a[jp    ], a[jt+p01],a[jp+p01], a[jt+p03],a[jp+p03], a[jt+p02],a[jp+p02], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15, r48r,r48i,r49r,r49i,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i, rt,it);

			/*...and now do 8 radix-7 transforms, with the columns of r*[r,i] output pairs in the above 7x radix-8 set now acting as input rows: */
							/*                            inputs                                */  /*                    intermediates                  */  /*                                                                     outputs                       */
				RADIX_07_DFT(r00r,r00i, r08r,r08i, r16r,r16i, r24r,r24i, r32r,r32i, r40r,r40i, r48r,r48i,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p00r,a1p00i, a1p08r,a1p08i, a1p16r,a1p16i, a1p24r,a1p24i, a1p32r,a1p32i, a1p40r,a1p40i, a1p48r,a1p48i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(r01r,r01i, r09r,r09i, r17r,r17i, r25r,r25i, r33r,r33i, r41r,r41i, r49r,r49i,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p49r,a1p49i, a1p01r,a1p01i, a1p09r,a1p09i, a1p17r,a1p17i, a1p25r,a1p25i, a1p33r,a1p33i, a1p41r,a1p41i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(r02r,r02i, r10r,r10i, r18r,r18i, r26r,r26i, r34r,r34i, r42r,r42i, r50r,r50i,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p42r,a1p42i, a1p50r,a1p50i, a1p02r,a1p02i, a1p10r,a1p10i, a1p18r,a1p18i, a1p26r,a1p26i, a1p34r,a1p34i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(r03r,r03i, r11r,r11i, r19r,r19i, r27r,r27i, r35r,r35i, r43r,r43i, r51r,r51i,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p35r,a1p35i, a1p43r,a1p43i, a1p51r,a1p51i, a1p03r,a1p03i, a1p11r,a1p11i, a1p19r,a1p19i, a1p27r,a1p27i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(r04r,r04i, r12r,r12i, r20r,r20i, r28r,r28i, r36r,r36i, r44r,r44i, r52r,r52i,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p28r,a1p28i, a1p36r,a1p36i, a1p44r,a1p44i, a1p52r,a1p52i, a1p04r,a1p04i, a1p12r,a1p12i, a1p20r,a1p20i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(r05r,r05i, r13r,r13i, r21r,r21i, r29r,r29i, r37r,r37i, r45r,r45i, r53r,r53i,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p21r,a1p21i, a1p29r,a1p29i, a1p37r,a1p37i, a1p45r,a1p45i, a1p53r,a1p53i, a1p05r,a1p05i, a1p13r,a1p13i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(r06r,r06i, r14r,r14i, r22r,r22i, r30r,r30i, r38r,r38i, r46r,r46i, r54r,r54i,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p14r,a1p14i, a1p22r,a1p22i, a1p30r,a1p30i, a1p38r,a1p38i, a1p46r,a1p46i, a1p54r,a1p54i, a1p06r,a1p06i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(r07r,r07i, r15r,r15i, r23r,r23i, r31r,r31i, r39r,r39i, r47r,r47i, r55r,r55i,  t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a1p07r,a1p07i, a1p15r,a1p15i, a1p23r,a1p23i, a1p31r,a1p31i, a1p39r,a1p39i, a1p47r,a1p47i, a1p55r,a1p55i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);

			#endif	/* USE_SSE2 */

			/*...Now do the carries. Since the outputs would
				normally be getting dispatched to 28 separate blocks of the A-array, we need 28 separate carries.	*/

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

					AVX_cmplx_carry_norm_errcheck0_X4(s1p00r,add1,add2,add3,cy_r00,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					AVX_cmplx_carry_norm_errcheck1_X4(s1p04r,add1,add2,add3,cy_r04,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					AVX_cmplx_carry_norm_errcheck1_X4(s1p08r,add1,add2,add3,cy_r08,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					AVX_cmplx_carry_norm_errcheck1_X4(s1p12r,add1,add2,add3,cy_r12,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					AVX_cmplx_carry_norm_errcheck1_X4(s1p16r,add1,add2,add3,cy_r16,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					AVX_cmplx_carry_norm_errcheck1_X4(s1p20r,add1,add2,add3,cy_r20,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					AVX_cmplx_carry_norm_errcheck1_X4(s1p24r,add1,add2,add3,cy_r24,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					AVX_cmplx_carry_norm_errcheck1_X4(s1p28r,add1,add2,add3,cy_r28,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					AVX_cmplx_carry_norm_errcheck1_X4(s1p32r,add1,add2,add3,cy_r32,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					AVX_cmplx_carry_norm_errcheck1_X4(s1p36r,add1,add2,add3,cy_r36,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					AVX_cmplx_carry_norm_errcheck1_X4(s1p40r,add1,add2,add3,cy_r40,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					AVX_cmplx_carry_norm_errcheck1_X4(s1p44r,add1,add2,add3,cy_r44,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					AVX_cmplx_carry_norm_errcheck1_X4(s1p48r,add1,add2,add3,cy_r48,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					AVX_cmplx_carry_norm_errcheck1_X4(s1p52r,add1,add2,add3,cy_r52,bjmodn52,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);

					co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
								// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

					i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

				#elif defined(USE_SSE2)

					l= j & (nwt-1);
					n_minus_sil   = n-si[l  ];
					n_minus_silp1 = n-si[l+1];
					sinwt   = si[nwt-l  ];
					sinwtm1 = si[nwt-l-1];

					wtl     =wt0[    l  ];
					wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
					wtlp1   =wt0[    l+1];
					wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

					ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
					ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
					ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
					ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
					ctmp->re = wtnm1;	ctmp->im = wtnm1;

					add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
					add2 = &wt1[co2-1];
					add3 = &wt1[co3-1];

				  #if defined(COMPILER_TYPE_MSVC)

				   #ifdef ERR_CHECK_ALL				/* Updating i prior to the 2nd-7th macro calls allows use of the same 0_2B macro for all */
					SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00);	// i =((uint32)(sw - *bjmodn04) >> 31);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04);	// i =((uint32)(sw - *bjmodn08) >> 31);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08);	// i =((uint32)(sw - *bjmodn12) >> 31);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12);	// i =((uint32)(sw - *bjmodn16) >> 31);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16);	// i =((uint32)(sw - *bjmodn20) >> 31);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20);	// i =((uint32)(sw - *bjmodn24) >> 31);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24);
				   #else
					SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24);
				   #endif

				  #else	/* GCC-style inline ASM: */

				   #ifdef ERR_CHECK_ALL
					SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p28r,add1,add2,add3,cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p32r,add1,add2,add3,cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p36r,add1,add2,add3,cy_r36,cy_r38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p40r,add1,add2,add3,cy_r40,cy_r42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p44r,add1,add2,add3,cy_r44,cy_r46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p48r,add1,add2,add3,cy_r48,cy_r50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck1_2B(s1p52r,add1,add2,add3,cy_r52,cy_r54,bjmodn52,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				   #else
					SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p28r,add1,add2,add3,cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p32r,add1,add2,add3,cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p36r,add1,add2,add3,cy_r36,cy_r38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p40r,add1,add2,add3,cy_r40,cy_r42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p44r,add1,add2,add3,cy_r44,cy_r46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p48r,add1,add2,add3,cy_r48,cy_r50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck1_2B (s1p52r,add1,add2,add3,cy_r52,cy_r54,bjmodn52,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				   #endif

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
					ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
					ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
					ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
					ctmp->re = wtnm1;	ctmp->im = wtnm1;

				/*	i =((uint32)(sw - *bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

					co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
								(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

					add1 = &wt1[col  ];
					add2 = &wt1[co2-1];

				  #if defined(COMPILER_TYPE_MSVC)

				   #ifdef ERR_CHECK_ALL
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24);
				   #else
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24);
				   #endif

				  #else	/* GCC-style inline ASM: */

				   #ifdef ERR_CHECK_ALL
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p28r,add1,add2,     cy_r28,cy_r30,bjmodn28,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p32r,add1,add2,     cy_r32,cy_r34,bjmodn32,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p36r,add1,add2,     cy_r36,cy_r38,bjmodn36,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p40r,add1,add2,     cy_r40,cy_r42,bjmodn40,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p44r,add1,add2,     cy_r44,cy_r46,bjmodn44,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p48r,add1,add2,     cy_r48,cy_r50,bjmodn48,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_errcheck2_2B(s1p52r,add1,add2,     cy_r52,cy_r54,bjmodn52,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				   #else
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p28r,add1,add2,     cy_r28,cy_r30,bjmodn28,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p32r,add1,add2,     cy_r32,cy_r34,bjmodn32,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p36r,add1,add2,     cy_r36,cy_r38,bjmodn36,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p40r,add1,add2,     cy_r40,cy_r42,bjmodn40,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p44r,add1,add2,     cy_r44,cy_r46,bjmodn44,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p48r,add1,add2,     cy_r48,cy_r50,bjmodn48,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
					SSE2_cmplx_carry_norm_nocheck2_2B (s1p52r,add1,add2,     cy_r52,cy_r54,bjmodn52,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				   #endif

				  #endif

					i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

				#else	// Scalar-double mode:

					l= j & (nwt-1);
					n_minus_sil   = n-si[l  ];
					n_minus_silp1 = n-si[l+1];
					sinwt   = si[nwt-l  ];
					sinwtm1 = si[nwt-l-1];

					wtl     =wt0[    l  ];
					wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
					wtlp1   =wt0[    l+1];
					wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				/*...set0 is slightly different from others:	*/
				   cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy_r00,bjmodn00   );
					cmplx_carry_norm_errcheck(a1p01r,a1p01i,cy_r01,bjmodn01,1 );
					cmplx_carry_norm_errcheck(a1p02r,a1p02i,cy_r02,bjmodn02,2 );
					cmplx_carry_norm_errcheck(a1p03r,a1p03i,cy_r03,bjmodn03,3 );
					cmplx_carry_norm_errcheck(a1p04r,a1p04i,cy_r04,bjmodn04,4 );
					cmplx_carry_norm_errcheck(a1p05r,a1p05i,cy_r05,bjmodn05,5 );
					cmplx_carry_norm_errcheck(a1p06r,a1p06i,cy_r06,bjmodn06,6 );
					cmplx_carry_norm_errcheck(a1p07r,a1p07i,cy_r07,bjmodn07,7 );
					cmplx_carry_norm_errcheck(a1p08r,a1p08i,cy_r08,bjmodn08,8 );
					cmplx_carry_norm_errcheck(a1p09r,a1p09i,cy_r09,bjmodn09,9 );
					cmplx_carry_norm_errcheck(a1p10r,a1p10i,cy_r10,bjmodn10,10);
					cmplx_carry_norm_errcheck(a1p11r,a1p11i,cy_r11,bjmodn11,11);
					cmplx_carry_norm_errcheck(a1p12r,a1p12i,cy_r12,bjmodn12,12);
					cmplx_carry_norm_errcheck(a1p13r,a1p13i,cy_r13,bjmodn13,13);
					cmplx_carry_norm_errcheck(a1p14r,a1p14i,cy_r14,bjmodn14,14);
					cmplx_carry_norm_errcheck(a1p15r,a1p15i,cy_r15,bjmodn15,15);
					cmplx_carry_norm_errcheck(a1p16r,a1p16i,cy_r16,bjmodn16,16);
					cmplx_carry_norm_errcheck(a1p17r,a1p17i,cy_r17,bjmodn17,17);
					cmplx_carry_norm_errcheck(a1p18r,a1p18i,cy_r18,bjmodn18,18);
					cmplx_carry_norm_errcheck(a1p19r,a1p19i,cy_r19,bjmodn19,19);
					cmplx_carry_norm_errcheck(a1p20r,a1p20i,cy_r20,bjmodn20,20);
					cmplx_carry_norm_errcheck(a1p21r,a1p21i,cy_r21,bjmodn21,21);
					cmplx_carry_norm_errcheck(a1p22r,a1p22i,cy_r22,bjmodn22,22);
					cmplx_carry_norm_errcheck(a1p23r,a1p23i,cy_r23,bjmodn23,23);
					cmplx_carry_norm_errcheck(a1p24r,a1p24i,cy_r24,bjmodn24,24);
					cmplx_carry_norm_errcheck(a1p25r,a1p25i,cy_r25,bjmodn25,25);
					cmplx_carry_norm_errcheck(a1p26r,a1p26i,cy_r26,bjmodn26,26);
					cmplx_carry_norm_errcheck(a1p27r,a1p27i,cy_r27,bjmodn27,27);
					cmplx_carry_norm_errcheck(a1p28r,a1p28i,cy_r28,bjmodn28,28);
					cmplx_carry_norm_errcheck(a1p29r,a1p29i,cy_r29,bjmodn29,29);
					cmplx_carry_norm_errcheck(a1p30r,a1p30i,cy_r30,bjmodn30,30);
					cmplx_carry_norm_errcheck(a1p31r,a1p31i,cy_r31,bjmodn31,31);
					cmplx_carry_norm_errcheck(a1p32r,a1p32i,cy_r32,bjmodn32,32);
					cmplx_carry_norm_errcheck(a1p33r,a1p33i,cy_r33,bjmodn33,33);
					cmplx_carry_norm_errcheck(a1p34r,a1p34i,cy_r34,bjmodn34,34);
					cmplx_carry_norm_errcheck(a1p35r,a1p35i,cy_r35,bjmodn35,35);
					cmplx_carry_norm_errcheck(a1p36r,a1p36i,cy_r36,bjmodn36,36);
					cmplx_carry_norm_errcheck(a1p37r,a1p37i,cy_r37,bjmodn37,37);
					cmplx_carry_norm_errcheck(a1p38r,a1p38i,cy_r38,bjmodn38,38);
					cmplx_carry_norm_errcheck(a1p39r,a1p39i,cy_r39,bjmodn39,39);
					cmplx_carry_norm_errcheck(a1p40r,a1p40i,cy_r40,bjmodn40,40);
					cmplx_carry_norm_errcheck(a1p41r,a1p41i,cy_r41,bjmodn41,41);
					cmplx_carry_norm_errcheck(a1p42r,a1p42i,cy_r42,bjmodn42,42);
					cmplx_carry_norm_errcheck(a1p43r,a1p43i,cy_r43,bjmodn43,43);
					cmplx_carry_norm_errcheck(a1p44r,a1p44i,cy_r44,bjmodn44,44);
					cmplx_carry_norm_errcheck(a1p45r,a1p45i,cy_r45,bjmodn45,45);
					cmplx_carry_norm_errcheck(a1p46r,a1p46i,cy_r46,bjmodn46,46);
					cmplx_carry_norm_errcheck(a1p47r,a1p47i,cy_r47,bjmodn47,47);
					cmplx_carry_norm_errcheck(a1p48r,a1p48i,cy_r48,bjmodn48,48);
					cmplx_carry_norm_errcheck(a1p49r,a1p49i,cy_r49,bjmodn49,49);
					cmplx_carry_norm_errcheck(a1p50r,a1p50i,cy_r50,bjmodn50,50);
					cmplx_carry_norm_errcheck(a1p51r,a1p51i,cy_r51,bjmodn51,51);
					cmplx_carry_norm_errcheck(a1p52r,a1p52i,cy_r52,bjmodn52,52);
					cmplx_carry_norm_errcheck(a1p53r,a1p53i,cy_r53,bjmodn53,53);
					cmplx_carry_norm_errcheck(a1p54r,a1p54i,cy_r54,bjmodn54,54);
					cmplx_carry_norm_errcheck(a1p55r,a1p55i,cy_r55,bjmodn55,55);

					i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/

					co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
								(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

				#endif	// USE_AVX?

				}
				else	/* MODULUS_TYPE_FERMAT */
				{

				#ifdef USE_AVX

					// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.

					/* Get the needed Nth root of -1: */
					add1 = (double *)&rn0[0];
					add2 = (double *)&rn1[0];

					idx_offset = j;
					idx_incr = NDIVR;

					tmp = base_negacyclic_root;	tm2 = tmp+1;

				  #if HIACC
					// Hi-accuracy version needs 7 copies of each base root:
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

				  #else	// HIACC = false:

					// Get the needed quartet of Nth roots of -1: This is the same code as in the scalar
					// fermat_carry_norm_errcheck() macro, with the single index j replaced by the quartet j,j+2,j+4,j+6:
					l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
					dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
					rt  =rn1[k2].re;			it   =rn1[k2].im;
					wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
					VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

					l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
					dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
					rt  =rn1[k2].re;			it   =rn1[k2].im;
					wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
					VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

					l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
					dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
					rt  =rn1[k2].re;			it   =rn1[k2].im;
					wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
					VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

					l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
					dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
					rt  =rn1[k2].re;			it   =rn1[k2].im;
					wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
					VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

					// The above need some inits to prepare for the AVX version of the Fermat-mod carry macro:
					SSE2_fermat_carry_init_loacc(base_negacyclic_root);

				  #endif

					// AVX-custom 4-way carry macro - each contains 4 of the RADIX stride-n/RADIX-separated carries
					// (processed independently in parallel), and steps through sequential-data indicies j,j+2,j+4,j+6:
				  #if HIACC
					// The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl)
					// which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
					// precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
					// uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls:
																																						// *cycle0 index increments by +4 (mod odd_radix) between macro calls
					tmp = base_negacyclic_root+  0;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p00r,tmp,0xe00,cy_r00,cy_i00,odd_radix,half_arr,sign_mask,icycle0,icycle1,icycle2,icycle3, jcycle0,kcycle0,lcycle0);
					tmp = base_negacyclic_root+  8;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p04r,tmp,0xd40,cy_r04,cy_i04,odd_radix,half_arr,sign_mask,icycle4,icycle5,icycle6,icycle0, jcycle4,kcycle4,lcycle4);
					tmp = base_negacyclic_root+ 16;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p08r,tmp,0xc80,cy_r08,cy_i08,odd_radix,half_arr,sign_mask,icycle1,icycle2,icycle3,icycle4, jcycle1,kcycle1,lcycle1);
					tmp = base_negacyclic_root+ 24;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p12r,tmp,0xbc0,cy_r12,cy_i12,odd_radix,half_arr,sign_mask,icycle5,icycle6,icycle0,icycle1, jcycle5,kcycle5,lcycle5);
					tmp = base_negacyclic_root+ 32;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p16r,tmp,0xb00,cy_r16,cy_i16,odd_radix,half_arr,sign_mask,icycle2,icycle3,icycle4,icycle5, jcycle2,kcycle2,lcycle2);
					tmp = base_negacyclic_root+ 40;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p20r,tmp,0xa40,cy_r20,cy_i20,odd_radix,half_arr,sign_mask,icycle6,icycle0,icycle1,icycle2, jcycle6,kcycle6,lcycle6);
					tmp = base_negacyclic_root+ 48;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p24r,tmp,0x980,cy_r24,cy_i24,odd_radix,half_arr,sign_mask,icycle3,icycle4,icycle5,icycle6, jcycle3,kcycle3,lcycle3);
					tmp = base_negacyclic_root+ 56;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p28r,tmp,0x8c0,cy_r28,cy_i28,odd_radix,half_arr,sign_mask,icycle0,icycle1,icycle2,icycle3, jcycle0,kcycle0,lcycle0);
					tmp = base_negacyclic_root+ 64;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p32r,tmp,0x800,cy_r32,cy_i32,odd_radix,half_arr,sign_mask,icycle4,icycle5,icycle6,icycle0, jcycle4,kcycle4,lcycle4);
					tmp = base_negacyclic_root+ 72;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p36r,tmp,0x740,cy_r36,cy_i36,odd_radix,half_arr,sign_mask,icycle1,icycle2,icycle3,icycle4, jcycle1,kcycle1,lcycle1);
					tmp = base_negacyclic_root+ 80;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p40r,tmp,0x680,cy_r40,cy_i40,odd_radix,half_arr,sign_mask,icycle5,icycle6,icycle0,icycle1, jcycle5,kcycle5,lcycle5);
					tmp = base_negacyclic_root+ 88;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p44r,tmp,0x5c0,cy_r44,cy_i44,odd_radix,half_arr,sign_mask,icycle2,icycle3,icycle4,icycle5, jcycle2,kcycle2,lcycle2);
					tmp = base_negacyclic_root+ 96;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p48r,tmp,0x500,cy_r48,cy_i48,odd_radix,half_arr,sign_mask,icycle6,icycle0,icycle1,icycle2, jcycle6,kcycle6,lcycle6);
					tmp = base_negacyclic_root+104;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p52r,tmp,0x440,cy_r52,cy_i52,odd_radix,half_arr,sign_mask,icycle3,icycle4,icycle5,icycle6, jcycle3,kcycle3,lcycle3);

				  #else	// HIACC = false:
																																						// *cycle0 index increments by +4 (mod odd_radix) between macro calls
					SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p00r,base_negacyclic_root,cy_r00,cy_i00,odd_radix,half_arr,sign_mask,                     icycle0,icycle1,icycle2,icycle3, jcycle0,kcycle0,lcycle0);
					SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p04r,base_negacyclic_root,cy_r04,cy_i04,odd_radix,half_arr,sign_mask,                     icycle4,icycle5,icycle6,icycle0, jcycle4,kcycle4,lcycle4);
					SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p08r,base_negacyclic_root,cy_r08,cy_i08,odd_radix,half_arr,sign_mask,                     icycle1,icycle2,icycle3,icycle4, jcycle1,kcycle1,lcycle1);
					SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p12r,base_negacyclic_root,cy_r12,cy_i12,odd_radix,half_arr,sign_mask,                     icycle5,icycle6,icycle0,icycle1, jcycle5,kcycle5,lcycle5);
					SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p16r,base_negacyclic_root,cy_r16,cy_i16,odd_radix,half_arr,sign_mask,                     icycle2,icycle3,icycle4,icycle5, jcycle2,kcycle2,lcycle2);
					SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p20r,base_negacyclic_root,cy_r20,cy_i20,odd_radix,half_arr,sign_mask,                     icycle6,icycle0,icycle1,icycle2, jcycle6,kcycle6,lcycle6);
					SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p24r,base_negacyclic_root,cy_r24,cy_i24,odd_radix,half_arr,sign_mask,                     icycle3,icycle4,icycle5,icycle6, jcycle3,kcycle3,lcycle3);
					SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p28r,base_negacyclic_root,cy_r28,cy_i28,odd_radix,half_arr,sign_mask,                     icycle0,icycle1,icycle2,icycle3, jcycle0,kcycle0,lcycle0);
					SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p32r,base_negacyclic_root,cy_r32,cy_i32,odd_radix,half_arr,sign_mask,                     icycle4,icycle5,icycle6,icycle0, jcycle4,kcycle4,lcycle4);
					SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p36r,base_negacyclic_root,cy_r36,cy_i36,odd_radix,half_arr,sign_mask,                     icycle1,icycle2,icycle3,icycle4, jcycle1,kcycle1,lcycle1);
					SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p40r,base_negacyclic_root,cy_r40,cy_i40,odd_radix,half_arr,sign_mask,                     icycle5,icycle6,icycle0,icycle1, jcycle5,kcycle5,lcycle5);
					SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p44r,base_negacyclic_root,cy_r44,cy_i44,odd_radix,half_arr,sign_mask,                     icycle2,icycle3,icycle4,icycle5, jcycle2,kcycle2,lcycle2);
					SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p48r,base_negacyclic_root,cy_r48,cy_i48,odd_radix,half_arr,sign_mask,                     icycle6,icycle0,icycle1,icycle2, jcycle6,kcycle6,lcycle6);
					SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p52r,base_negacyclic_root,cy_r52,cy_i52,odd_radix,half_arr,sign_mask,                     icycle3,icycle4,icycle5,icycle6, jcycle3,kcycle3,lcycle3);

				  #endif	// HIACC?

				#elif defined(USE_SSE2)

					/* Get the needed Nth root of -1: */
					add1 = (double *)&rn0[0];
					add2 = (double *)&rn1[0];

					idx_offset = j;
					idx_incr = NDIVR;

				  #if defined(COMPILER_TYPE_MSVC)

				  /* The cy_[r|i]_idx[A|B] names here are not meaningful, each simply stores one [re,im] carry pair,
				  e.g. cy_r01 stores the carries our of [a0.re,a0.im], cy_r23 stores the carries our of [a1.re,a1.im], etc.
				  Here is the actual mapping between these SSE2-mode 2-vector carry pairs and the scalar carries:
																		2-vector				                                          Scalar
																		--------	 		                                           ------------- */
					SSE2_fermat_carry_norm_errcheck(s1p00r,cy_r00,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r00,cy_i00 */
					SSE2_fermat_carry_norm_errcheck(s1p01r,cy_r02,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r01,cy_i01 */
					SSE2_fermat_carry_norm_errcheck(s1p02r,cy_r04,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r02,cy_i02 */
					SSE2_fermat_carry_norm_errcheck(s1p03r,cy_r06,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r03,cy_i03 */
					SSE2_fermat_carry_norm_errcheck(s1p04r,cy_r08,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r04,cy_i04 */
					SSE2_fermat_carry_norm_errcheck(s1p05r,cy_r10,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r05,cy_i05 */
					SSE2_fermat_carry_norm_errcheck(s1p06r,cy_r12,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r06,cy_i06 */
					SSE2_fermat_carry_norm_errcheck(s1p07r,cy_r14,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r07,cy_i07 */
					SSE2_fermat_carry_norm_errcheck(s1p08r,cy_r16,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r08,cy_i08 */
					SSE2_fermat_carry_norm_errcheck(s1p09r,cy_r18,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r09,cy_i09 */
					SSE2_fermat_carry_norm_errcheck(s1p10r,cy_r20,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r10,cy_i10 */
					SSE2_fermat_carry_norm_errcheck(s1p11r,cy_r22,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r11,cy_i11 */
					SSE2_fermat_carry_norm_errcheck(s1p12r,cy_r24,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r12,cy_i12 */
					SSE2_fermat_carry_norm_errcheck(s1p13r,cy_r26,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r13,cy_i13 */
					SSE2_fermat_carry_norm_errcheck(s1p14r,cy_r28,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r14,cy_i14 */
					SSE2_fermat_carry_norm_errcheck(s1p15r,cy_r30,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r15,cy_i15 */
					SSE2_fermat_carry_norm_errcheck(s1p16r,cy_r32,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r16,cy_i16 */
					SSE2_fermat_carry_norm_errcheck(s1p17r,cy_r34,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r17,cy_i17 */
					SSE2_fermat_carry_norm_errcheck(s1p18r,cy_r36,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r18,cy_i18 */
					SSE2_fermat_carry_norm_errcheck(s1p19r,cy_r38,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r19,cy_i19 */
					SSE2_fermat_carry_norm_errcheck(s1p20r,cy_r40,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r20,cy_i20 */
					SSE2_fermat_carry_norm_errcheck(s1p21r,cy_r42,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r21,cy_i21 */
					SSE2_fermat_carry_norm_errcheck(s1p22r,cy_r44,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r22,cy_i22 */
					SSE2_fermat_carry_norm_errcheck(s1p23r,cy_r46,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r23,cy_i23 */
					SSE2_fermat_carry_norm_errcheck(s1p24r,cy_r48,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r24,cy_i24 */
					SSE2_fermat_carry_norm_errcheck(s1p25r,cy_r50,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r25,cy_i25 */
					SSE2_fermat_carry_norm_errcheck(s1p26r,cy_r52,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r26,cy_i26 */
					SSE2_fermat_carry_norm_errcheck(s1p27r,cy_r54,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r27,cy_i27 */

					SSE2_fermat_carry_norm_errcheck(s1pr28,cy_i00,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r28,cy_i28 */
					SSE2_fermat_carry_norm_errcheck(s1pr29,cy_i02,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r29,cy_i29 */
					SSE2_fermat_carry_norm_errcheck(s1pr30,cy_i04,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r30,cy_i30 */
					SSE2_fermat_carry_norm_errcheck(s1pr31,cy_i06,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r31,cy_i31 */
					SSE2_fermat_carry_norm_errcheck(s1pr32,cy_i08,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r32,cy_i32 */
					SSE2_fermat_carry_norm_errcheck(s1pr33,cy_i10,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r33,cy_i33 */
					SSE2_fermat_carry_norm_errcheck(s1pr34,cy_i12,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r34,cy_i34 */
					SSE2_fermat_carry_norm_errcheck(s1pr35,cy_i14,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r35,cy_i35 */
					SSE2_fermat_carry_norm_errcheck(s1pr36,cy_i16,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r36,cy_i36 */
					SSE2_fermat_carry_norm_errcheck(s1pr37,cy_i18,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r37,cy_i37 */
					SSE2_fermat_carry_norm_errcheck(s1pr38,cy_i20,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r38,cy_i38 */
					SSE2_fermat_carry_norm_errcheck(s1pr39,cy_i22,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r39,cy_i39 */
					SSE2_fermat_carry_norm_errcheck(s1pr40,cy_i24,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r40,cy_i40 */
					SSE2_fermat_carry_norm_errcheck(s1pr41,cy_i26,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r41,cy_i41 */
					SSE2_fermat_carry_norm_errcheck(s1pr42,cy_i28,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r42,cy_i42 */
					SSE2_fermat_carry_norm_errcheck(s1pr43,cy_i30,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r43,cy_i43 */
					SSE2_fermat_carry_norm_errcheck(s1pr44,cy_i32,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r44,cy_i44 */
					SSE2_fermat_carry_norm_errcheck(s1pr45,cy_i34,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r45,cy_i45 */
					SSE2_fermat_carry_norm_errcheck(s1pr46,cy_i36,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r46,cy_i46 */
					SSE2_fermat_carry_norm_errcheck(s1pr47,cy_i38,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r47,cy_i47 */
					SSE2_fermat_carry_norm_errcheck(s1pr48,cy_i40,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r48,cy_i48 */
					SSE2_fermat_carry_norm_errcheck(s1pr49,cy_i42,idx_offset,idx_incr,odd_radix,icycle0,jcycle0,NRTM1,NRT_BITS);	/* cy_r49,cy_i49 */
					SSE2_fermat_carry_norm_errcheck(s1pr50,cy_i44,idx_offset,idx_incr,odd_radix,icycle1,jcycle1,NRTM1,NRT_BITS);	/* cy_r50,cy_i50 */
					SSE2_fermat_carry_norm_errcheck(s1pr51,cy_i46,idx_offset,idx_incr,odd_radix,icycle2,jcycle2,NRTM1,NRT_BITS);	/* cy_r51,cy_i51 */
					SSE2_fermat_carry_norm_errcheck(s1pr52,cy_i48,idx_offset,idx_incr,odd_radix,icycle3,jcycle3,NRTM1,NRT_BITS);	/* cy_r52,cy_i52 */
					SSE2_fermat_carry_norm_errcheck(s1pr53,cy_i50,idx_offset,idx_incr,odd_radix,icycle4,jcycle4,NRTM1,NRT_BITS);	/* cy_r53,cy_i53 */
					SSE2_fermat_carry_norm_errcheck(s1pr54,cy_i52,idx_offset,idx_incr,odd_radix,icycle5,jcycle5,NRTM1,NRT_BITS);	/* cy_r54,cy_i54 */
					SSE2_fermat_carry_norm_errcheck(s1pr55,cy_i54,idx_offset,idx_incr,odd_radix,icycle6,jcycle6,NRTM1,NRT_BITS);	/* cy_r55,cy_i55 */

				  #elif (OS_BITS == 32)

					SSE2_fermat_carry_norm_errcheck(s1p00r,cy_r00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck(s1p01r,cy_r02,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck(s1p02r,cy_r04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck(s1p03r,cy_r06,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck(s1p04r,cy_r08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck(s1p05r,cy_r10,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck(s1p06r,cy_r12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
					SSE2_fermat_carry_norm_errcheck(s1p07r,cy_r14,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck(s1p08r,cy_r16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck(s1p09r,cy_r18,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck(s1p10r,cy_r20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck(s1p11r,cy_r22,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck(s1p12r,cy_r24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck(s1p13r,cy_r26,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
					SSE2_fermat_carry_norm_errcheck(s1p14r,cy_r28,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck(s1p15r,cy_r30,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck(s1p16r,cy_r32,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck(s1p17r,cy_r34,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck(s1p18r,cy_r36,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck(s1p19r,cy_r38,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck(s1p20r,cy_r40,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
					SSE2_fermat_carry_norm_errcheck(s1p21r,cy_r42,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck(s1p22r,cy_r44,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck(s1p23r,cy_r46,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck(s1p24r,cy_r48,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck(s1p25r,cy_r50,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck(s1p26r,cy_r52,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck(s1p27r,cy_r54,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);

					SSE2_fermat_carry_norm_errcheck(s1p28r,cy_i00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck(s1p29r,cy_i02,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck(s1p30r,cy_i04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck(s1p31r,cy_i06,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck(s1p32r,cy_i08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck(s1p33r,cy_i10,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck(s1p34r,cy_i12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
					SSE2_fermat_carry_norm_errcheck(s1p35r,cy_i14,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck(s1p36r,cy_i16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck(s1p37r,cy_i18,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck(s1p38r,cy_i20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck(s1p39r,cy_i22,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck(s1p40r,cy_i24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck(s1p41r,cy_i26,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
					SSE2_fermat_carry_norm_errcheck(s1p42r,cy_i28,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck(s1p43r,cy_i30,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck(s1p44r,cy_i32,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck(s1p45r,cy_i34,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck(s1p46r,cy_i36,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck(s1p47r,cy_i38,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck(s1p48r,cy_i40,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);
					SSE2_fermat_carry_norm_errcheck(s1p49r,cy_i42,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck(s1p50r,cy_i44,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck(s1p51r,cy_i46,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck(s1p52r,cy_i48,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck(s1p53r,cy_i50,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck(s1p54r,cy_i52,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck(s1p55r,cy_i54,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6);

				  #else	// 64-bit SSE2

					SSE2_fermat_carry_norm_errcheck_X2(s1p00r,cy_r00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck_X2(s1p02r,cy_r04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck_X2(s1p04r,cy_r08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck_X2(s1p06r,cy_r12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck_X2(s1p08r,cy_r16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck_X2(s1p10r,cy_r20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck_X2(s1p12r,cy_r24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5,icycle6,jcycle6);
					SSE2_fermat_carry_norm_errcheck_X2(s1p14r,cy_r28,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck_X2(s1p16r,cy_r32,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck_X2(s1p18r,cy_r36,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck_X2(s1p20r,cy_r40,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck_X2(s1p22r,cy_r44,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck_X2(s1p24r,cy_r48,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck_X2(s1p26r,cy_r52,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5,icycle6,jcycle6);
					SSE2_fermat_carry_norm_errcheck_X2(s1p28r,cy_i00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck_X2(s1p30r,cy_i04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck_X2(s1p32r,cy_i08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck_X2(s1p34r,cy_i12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck_X2(s1p36r,cy_i16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck_X2(s1p38r,cy_i20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck_X2(s1p40r,cy_i24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5,icycle6,jcycle6);
					SSE2_fermat_carry_norm_errcheck_X2(s1p42r,cy_i28,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle0,jcycle0,icycle1,jcycle1);
					SSE2_fermat_carry_norm_errcheck_X2(s1p44r,cy_i32,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle2,jcycle2,icycle3,jcycle3);
					SSE2_fermat_carry_norm_errcheck_X2(s1p46r,cy_i36,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle4,jcycle4,icycle5,jcycle5);
					SSE2_fermat_carry_norm_errcheck_X2(s1p48r,cy_i40,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle6,jcycle6,icycle0,jcycle0);
					SSE2_fermat_carry_norm_errcheck_X2(s1p50r,cy_i44,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle1,jcycle1,icycle2,jcycle2);
					SSE2_fermat_carry_norm_errcheck_X2(s1p52r,cy_i48,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle3,jcycle3,icycle4,jcycle4);
					SSE2_fermat_carry_norm_errcheck_X2(s1p54r,cy_i52,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle5,jcycle5,icycle6,jcycle6);
				  #endif

				#else	// Scalar-double mode:

					ntmp = 0;
					fermat_carry_norm_errcheckB(a1p00r,a1p00i,cy_r00,cy_i00,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p01r,a1p01i,cy_r01,cy_i01,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p02r,a1p02i,cy_r02,cy_i02,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p03r,a1p03i,cy_r03,cy_i03,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p04r,a1p04i,cy_r04,cy_i04,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p05r,a1p05i,cy_r05,cy_i05,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p06r,a1p06i,cy_r06,cy_i06,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p07r,a1p07i,cy_r07,cy_i07,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p08r,a1p08i,cy_r08,cy_i08,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p09r,a1p09i,cy_r09,cy_i09,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p10r,a1p10i,cy_r10,cy_i10,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p11r,a1p11i,cy_r11,cy_i11,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p12r,a1p12i,cy_r12,cy_i12,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p13r,a1p13i,cy_r13,cy_i13,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p14r,a1p14i,cy_r14,cy_i14,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p15r,a1p15i,cy_r15,cy_i15,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p16r,a1p16i,cy_r16,cy_i16,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p17r,a1p17i,cy_r17,cy_i17,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p18r,a1p18i,cy_r18,cy_i18,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p19r,a1p19i,cy_r19,cy_i19,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p20r,a1p20i,cy_r20,cy_i20,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p21r,a1p21i,cy_r21,cy_i21,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p22r,a1p22i,cy_r22,cy_i22,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p23r,a1p23i,cy_r23,cy_i23,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p24r,a1p24i,cy_r24,cy_i24,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p25r,a1p25i,cy_r25,cy_i25,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p26r,a1p26i,cy_r26,cy_i26,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p27r,a1p27i,cy_r27,cy_i27,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p28r,a1p28i,cy_r28,cy_i28,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p29r,a1p29i,cy_r29,cy_i29,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p30r,a1p30i,cy_r30,cy_i30,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p31r,a1p31i,cy_r31,cy_i31,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p32r,a1p32i,cy_r32,cy_i32,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p33r,a1p33i,cy_r33,cy_i33,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p34r,a1p34i,cy_r34,cy_i34,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p35r,a1p35i,cy_r35,cy_i35,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p36r,a1p36i,cy_r36,cy_i36,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p37r,a1p37i,cy_r37,cy_i37,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p38r,a1p38i,cy_r38,cy_i38,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p39r,a1p39i,cy_r39,cy_i39,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p40r,a1p40i,cy_r40,cy_i40,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p41r,a1p41i,cy_r41,cy_i41,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p42r,a1p42i,cy_r42,cy_i42,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p43r,a1p43i,cy_r43,cy_i43,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p44r,a1p44i,cy_r44,cy_i44,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p45r,a1p45i,cy_r45,cy_i45,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p46r,a1p46i,cy_r46,cy_i46,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p47r,a1p47i,cy_r47,cy_i47,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p48r,a1p48i,cy_r48,cy_i48,icycle6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p49r,a1p49i,cy_r49,cy_i49,icycle0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p50r,a1p50i,cy_r50,cy_i50,icycle1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p51r,a1p51i,cy_r51,cy_i51,icycle2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p52r,a1p52i,cy_r52,cy_i52,icycle3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p53r,a1p53i,cy_r53,cy_i53,icycle4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p54r,a1p54i,cy_r54,cy_i54,icycle5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
					fermat_carry_norm_errcheckB(a1p55r,a1p55i,cy_r55,cy_i55,icycle6,ntmp,NRTM1,NRT_BITS);

					icycle0 += wts_idx_incr;	/* Inside the loop use this, as it is faster than general-mod '% nwt' */
					icycle1 += wts_idx_incr;
					icycle2 += wts_idx_incr;
					icycle3 += wts_idx_incr;
					icycle4 += wts_idx_incr;
					icycle5 += wts_idx_incr;
					icycle6 += wts_idx_incr;
					icycle0 += ( (-(int)((uint32)icycle0 >> 31)) & nwt);
					icycle1 += ( (-(int)((uint32)icycle1 >> 31)) & nwt);
					icycle2 += ( (-(int)((uint32)icycle2 >> 31)) & nwt);
					icycle3 += ( (-(int)((uint32)icycle3 >> 31)) & nwt);
					icycle4 += ( (-(int)((uint32)icycle4 >> 31)) & nwt);
					icycle5 += ( (-(int)((uint32)icycle5 >> 31)) & nwt);
					icycle6 += ( (-(int)((uint32)icycle6 >> 31)) & nwt);

				#endif	/* #ifdef USE_SSE2 */

				// Here we nest AVX inside SSE2 since i/jcycle updates are for both, k/l for AVX-only:
				#ifdef USE_SSE2

					icycle0 += wts_idx_inc2;		icycle0 += ( (-(icycle0 < 0)) & nwt16);
					icycle1 += wts_idx_inc2;		icycle1 += ( (-(icycle1 < 0)) & nwt16);
					icycle2 += wts_idx_inc2;		icycle2 += ( (-(icycle2 < 0)) & nwt16);
					icycle3 += wts_idx_inc2;		icycle3 += ( (-(icycle3 < 0)) & nwt16);
					icycle4 += wts_idx_inc2;		icycle4 += ( (-(icycle4 < 0)) & nwt16);
					icycle5 += wts_idx_inc2;		icycle5 += ( (-(icycle5 < 0)) & nwt16);
					icycle6 += wts_idx_inc2;		icycle6 += ( (-(icycle6 < 0)) & nwt16);

					jcycle0 += wts_idx_inc2;		jcycle0 += ( (-(jcycle0 < 0)) & nwt16);
					jcycle1 += wts_idx_inc2;		jcycle1 += ( (-(jcycle1 < 0)) & nwt16);
					jcycle2 += wts_idx_inc2;		jcycle2 += ( (-(jcycle2 < 0)) & nwt16);
					jcycle3 += wts_idx_inc2;		jcycle3 += ( (-(jcycle3 < 0)) & nwt16);
					jcycle4 += wts_idx_inc2;		jcycle4 += ( (-(jcycle4 < 0)) & nwt16);
					jcycle5 += wts_idx_inc2;		jcycle5 += ( (-(jcycle5 < 0)) & nwt16);
					jcycle6 += wts_idx_inc2;		jcycle6 += ( (-(jcycle6 < 0)) & nwt16);

				  #ifdef USE_AVX
					kcycle0 += wts_idx_inc2;		kcycle0 += ( (-(kcycle0 < 0)) & nwt16);
					kcycle1 += wts_idx_inc2;		kcycle1 += ( (-(kcycle1 < 0)) & nwt16);
					kcycle2 += wts_idx_inc2;		kcycle2 += ( (-(kcycle2 < 0)) & nwt16);
					kcycle3 += wts_idx_inc2;		kcycle3 += ( (-(kcycle3 < 0)) & nwt16);
					kcycle4 += wts_idx_inc2;		kcycle4 += ( (-(kcycle4 < 0)) & nwt16);
					kcycle5 += wts_idx_inc2;		kcycle5 += ( (-(kcycle5 < 0)) & nwt16);
					kcycle6 += wts_idx_inc2;		kcycle6 += ( (-(kcycle6 < 0)) & nwt16);

					lcycle0 += wts_idx_inc2;		lcycle0 += ( (-(lcycle0 < 0)) & nwt16);
					lcycle1 += wts_idx_inc2;		lcycle1 += ( (-(lcycle1 < 0)) & nwt16);
					lcycle2 += wts_idx_inc2;		lcycle2 += ( (-(lcycle2 < 0)) & nwt16);
					lcycle3 += wts_idx_inc2;		lcycle3 += ( (-(lcycle3 < 0)) & nwt16);
					lcycle4 += wts_idx_inc2;		lcycle4 += ( (-(lcycle4 < 0)) & nwt16);
					lcycle5 += wts_idx_inc2;		lcycle5 += ( (-(lcycle5 < 0)) & nwt16);
					lcycle6 += wts_idx_inc2;		lcycle6 += ( (-(lcycle6 < 0)) & nwt16);
				  #endif
				#endif

				}	/* if(MODULUS_TYPE == ...) */

		/*...The radix-56 DIF pass is here:	*/

			#ifdef USE_SSE2

			/*...gather the needed data (56 64-bit complex, i.e. 112 64-bit reals) and do 8 radix-7 transforms...*/
							 /*                   inputs                    */ /* sincos */ /*         outputs           */
				SSE2_RADIX_07_DFT(s1p00r,s1p48r,s1p40r,s1p32r,s1p24r,s1p16r,s1p08r, cc0, r00r,r08r,r16r,r24r,r32r,r40r,r48r);
				SSE2_RADIX_07_DFT(s1p49r,s1p41r,s1p33r,s1p25r,s1p17r,s1p09r,s1p01r, cc0, r01r,r09r,r17r,r25r,r33r,r41r,r49r);
				SSE2_RADIX_07_DFT(s1p42r,s1p34r,s1p26r,s1p18r,s1p10r,s1p02r,s1p50r, cc0, r02r,r10r,r18r,r26r,r34r,r42r,r50r);
				SSE2_RADIX_07_DFT(s1p35r,s1p27r,s1p19r,s1p11r,s1p03r,s1p51r,s1p43r, cc0, r03r,r11r,r19r,r27r,r35r,r43r,r51r);
				SSE2_RADIX_07_DFT(s1p28r,s1p20r,s1p12r,s1p04r,s1p52r,s1p44r,s1p36r, cc0, r04r,r12r,r20r,r28r,r36r,r44r,r52r);
				SSE2_RADIX_07_DFT(s1p21r,s1p13r,s1p05r,s1p53r,s1p45r,s1p37r,s1p29r, cc0, r05r,r13r,r21r,r29r,r37r,r45r,r53r);
				SSE2_RADIX_07_DFT(s1p14r,s1p06r,s1p54r,s1p46r,s1p38r,s1p30r,s1p22r, cc0, r06r,r14r,r22r,r30r,r38r,r46r,r54r);
				SSE2_RADIX_07_DFT(s1p07r,s1p55r,s1p47r,s1p39r,s1p31r,s1p23r,s1p15r, cc0, r07r,r15r,r23r,r31r,r39r,r47r,r55r);

			/*...and now do 7 radix-8 transforms: */
							 /*                                   inputs                                  */ /*                       intermediates                       */ /*                 outputs                   */
				add0 = &a[j1    ];
				add1 = add0+p01;
				add2 = add0+p02;
				add3 = add0+p03;
				add4 = add0+p04;
				add5 = add0+p05;
				add6 = add0+p06;
				add7 = add0+p07;
			  #ifdef USE_AVX
				SSE2_RADIX8_DIF_0TWIDDLE(r00r,0x40,0x80,0xc0,0x100,0x140,0x180,0x1c0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #else
				SSE2_RADIX8_DIF_0TWIDDLE(r00r,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #endif

				add2 = &a[j1+p30];
				add0 = add2+p03;
				add1 = add2+p02;
				add3 = add2+p01;
				add4 = add2+p07;
				add5 = add2+p06;
				add6 = add2+p04;
				add7 = add2+p05;
			  #ifdef USE_AVX
				SSE2_RADIX8_DIF_0TWIDDLE(r08r,0x40,0x80,0xc0,0x100,0x140,0x180,0x1c0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #else
				SSE2_RADIX8_DIF_0TWIDDLE(r08r,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #endif

				add6 = &a[j1+p28];
				add0 = add6+p05;
				add1 = add6+p04;
				add2 = add6+p07;
				add3 = add6+p06;
				add4 = add6+p03;
				add5 = add6+p02;
				add7 = add6+p01;
			  #ifdef USE_AVX
				SSE2_RADIX8_DIF_0TWIDDLE(r16r,0x40,0x80,0xc0,0x100,0x140,0x180,0x1c0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #else
				SSE2_RADIX8_DIF_0TWIDDLE(r16r,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #endif

				add1 = &a[j1+p20];
				add0 = add1+p01;
				add2 = add1+p03;
				add3 = add1+p02;
				add4 = add1+p05;
				add5 = add1+p04;
				add6 = add1+p07;
				add7 = add1+p06;
			  #ifdef USE_AVX
				SSE2_RADIX8_DIF_0TWIDDLE(r24r,0x40,0x80,0xc0,0x100,0x140,0x180,0x1c0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #else
				SSE2_RADIX8_DIF_0TWIDDLE(r24r,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #endif

				add5 = &a[j1+p18];
				add0 = add5+p06;
				add1 = add5+p07;
				add2 = add5+p05;
				add3 = add5+p04;
				add4 = add5+p01;
				add6 = add5+p03;
				add7 = add5+p02;
			  #ifdef USE_AVX
				SSE2_RADIX8_DIF_0TWIDDLE(r32r,0x40,0x80,0xc0,0x100,0x140,0x180,0x1c0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #else
				SSE2_RADIX8_DIF_0TWIDDLE(r32r,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #endif

				add3 = &a[j1+p10];
				add0 = add3+p02;
				add1 = add3+p03;
				add2 = add3+p01;
				add4 = add3+p06;
				add5 = add3+p07;
				add6 = add3+p05;
				add7 = add3+p04;
			  #ifdef USE_AVX
				SSE2_RADIX8_DIF_0TWIDDLE(r40r,0x40,0x80,0xc0,0x100,0x140,0x180,0x1c0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #else
				SSE2_RADIX8_DIF_0TWIDDLE(r40r,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #endif

				add7 = &a[j1+p08];
				add0 = add7+p04;
				add1 = add7+p05;
				add2 = add7+p06;
				add3 = add7+p07;
				add4 = add7+p02;
				add5 = add7+p03;
				add6 = add7+p01;
			  #ifdef USE_AVX
				SSE2_RADIX8_DIF_0TWIDDLE(r48r,0x40,0x80,0xc0,0x100,0x140,0x180,0x1c0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #else
				SSE2_RADIX8_DIF_0TWIDDLE(r48r,0x20,0x40,0x60,0x80,0xa0,0xc0,0xe0, add0,add1,add2,add3,add4,add5,add6,add7, isrt2);
			  #endif

			#else	// USE_SSE2 = False:

			/*...gather the needed data (56 64-bit complex, i.e. 112 64-bit reals) and do 8 radix-7 transforms...*/
							 /*                                               inputs                                              */  /*                   intermediates                   */  /*                                                  outputs        */  /*   sincos consts   */
				RADIX_07_DFT(a1p00r,a1p00i, a1p48r,a1p48i, a1p40r,a1p40i, a1p32r,a1p32i, a1p24r,a1p24i, a1p16r,a1p16i, a1p08r,a1p08i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r00r,r00i,r08r,r08i,r16r,r16i,r24r,r24i,r32r,r32i,r40r,r40i,r48r,r48i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(a1p49r,a1p49i, a1p41r,a1p41i, a1p33r,a1p33i, a1p25r,a1p25i, a1p17r,a1p17i, a1p09r,a1p09i, a1p01r,a1p01i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r01r,r01i,r09r,r09i,r17r,r17i,r25r,r25i,r33r,r33i,r41r,r41i,r49r,r49i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(a1p42r,a1p42i, a1p34r,a1p34i, a1p26r,a1p26i, a1p18r,a1p18i, a1p10r,a1p10i, a1p02r,a1p02i, a1p50r,a1p50i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r02r,r02i,r10r,r10i,r18r,r18i,r26r,r26i,r34r,r34i,r42r,r42i,r50r,r50i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(a1p35r,a1p35i, a1p27r,a1p27i, a1p19r,a1p19i, a1p11r,a1p11i, a1p03r,a1p03i, a1p51r,a1p51i, a1p43r,a1p43i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r03r,r03i,r11r,r11i,r19r,r19i,r27r,r27i,r35r,r35i,r43r,r43i,r51r,r51i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(a1p28r,a1p28i, a1p20r,a1p20i, a1p12r,a1p12i, a1p04r,a1p04i, a1p52r,a1p52i, a1p44r,a1p44i, a1p36r,a1p36i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r04r,r04i,r12r,r12i,r20r,r20i,r28r,r28i,r36r,r36i,r44r,r44i,r52r,r52i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(a1p21r,a1p21i, a1p13r,a1p13i, a1p05r,a1p05i, a1p53r,a1p53i, a1p45r,a1p45i, a1p37r,a1p37i, a1p29r,a1p29i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r05r,r05i,r13r,r13i,r21r,r21i,r29r,r29i,r37r,r37i,r45r,r45i,r53r,r53i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(a1p14r,a1p14i, a1p06r,a1p06i, a1p54r,a1p54i, a1p46r,a1p46i, a1p38r,a1p38i, a1p30r,a1p30i, a1p22r,a1p22i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r06r,r06i,r14r,r14i,r22r,r22i,r30r,r30i,r38r,r38i,r46r,r46i,r54r,r54i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
				RADIX_07_DFT(a1p07r,a1p07i, a1p55r,a1p55i, a1p47r,a1p47i, a1p39r,a1p39i, a1p31r,a1p31i, a1p23r,a1p23i, a1p15r,a1p15i, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, r07r,r07i,r15r,r15i,r23r,r23i,r31r,r31i,r39r,r39i,r47r,r47i,r55r,r55i, uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);

			/*...and now do 7 radix-8 transforms: */
							 /*                                   inputs                                  */ /*                       intermediates                       */ /*                 outputs                   */
				RADIX_08_DIF(r00r,r00i,r01r,r01i,r02r,r02i,r03r,r03i,r04r,r04i,r05r,r05i,r06r,r06i,r07r,r07i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05],a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],rt,it);	jt = j1+p30; jp = j2+p30;
				RADIX_08_DIF(r08r,r08i,r09r,r09i,r10r,r10i,r11r,r11i,r12r,r12i,r13r,r13i,r14r,r14i,r15r,r15i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],rt,it);	jt = j1+p28; jp = j2+p28;
				RADIX_08_DIF(r16r,r16i,r17r,r17i,r18r,r18i,r19r,r19i,r20r,r20i,r21r,r21i,r22r,r22i,r23r,r23i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	jt = j1+p20; jp = j2+p20;
				RADIX_08_DIF(r24r,r24i,r25r,r25i,r26r,r26i,r27r,r27i,r28r,r28i,r29r,r29i,r30r,r30i,r31r,r31i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],rt,it);	jt = j1+p18; jp = j2+p18;
				RADIX_08_DIF(r32r,r32i,r33r,r33i,r34r,r34i,r35r,r35i,r36r,r36i,r37r,r37i,r38r,r38i,r39r,r39i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p10; jp = j2+p10;
				RADIX_08_DIF(r40r,r40i,r41r,r41i,r42r,r42i,r43r,r43i,r44r,r44i,r45r,r45i,r46r,r46i,r47r,r47i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],rt,it);	jt = j1+p08; jp = j2+p08;
				RADIX_08_DIF(r48r,r48i,r49r,r49i,r50r,r50i,r51r,r51i,r52r,r52i,r53r,r53i,r54r,r54i,r55r,r55i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);

			#endif	// USE_SSE2 ?

			}	/* end for(j=_jstart; j < _jhi; j += stride) */

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
			thread_arg->cy_r28 = cy_r28->d0;
			thread_arg->cy_r29 = cy_r28->d1;
			thread_arg->cy_r30 = cy_r28->d2;
			thread_arg->cy_r31 = cy_r28->d3;
			thread_arg->cy_r32 = cy_r32->d0;
			thread_arg->cy_r33 = cy_r32->d1;
			thread_arg->cy_r34 = cy_r32->d2;
			thread_arg->cy_r35 = cy_r32->d3;
			thread_arg->cy_r36 = cy_r36->d0;
			thread_arg->cy_r37 = cy_r36->d1;
			thread_arg->cy_r38 = cy_r36->d2;
			thread_arg->cy_r39 = cy_r36->d3;
			thread_arg->cy_r40 = cy_r40->d0;
			thread_arg->cy_r41 = cy_r40->d1;
			thread_arg->cy_r42 = cy_r40->d2;
			thread_arg->cy_r43 = cy_r40->d3;
			thread_arg->cy_r44 = cy_r44->d0;
			thread_arg->cy_r45 = cy_r44->d1;
			thread_arg->cy_r46 = cy_r44->d2;
			thread_arg->cy_r47 = cy_r44->d3;
			thread_arg->cy_r48 = cy_r48->d0;
			thread_arg->cy_r49 = cy_r48->d1;
			thread_arg->cy_r50 = cy_r48->d2;
			thread_arg->cy_r51 = cy_r48->d3;
			thread_arg->cy_r52 = cy_r52->d0;
			thread_arg->cy_r53 = cy_r52->d1;
			thread_arg->cy_r54 = cy_r52->d2;
			thread_arg->cy_r55 = cy_r52->d3;
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
			thread_arg->cy_r28 = cy_r28->d0;
			thread_arg->cy_r29 = cy_r28->d1;
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
			thread_arg->cy_r40 = cy_r40->d0;
			thread_arg->cy_r41 = cy_r40->d1;
			thread_arg->cy_r42 = cy_r42->d0;
			thread_arg->cy_r43 = cy_r42->d1;
			thread_arg->cy_r44 = cy_r44->d0;
			thread_arg->cy_r45 = cy_r44->d1;
			thread_arg->cy_r46 = cy_r46->d0;
			thread_arg->cy_r47 = cy_r46->d1;
			thread_arg->cy_r48 = cy_r48->d0;
			thread_arg->cy_r49 = cy_r48->d1;
			thread_arg->cy_r50 = cy_r50->d0;
			thread_arg->cy_r51 = cy_r50->d1;
			thread_arg->cy_r52 = cy_r52->d0;
			thread_arg->cy_r53 = cy_r52->d1;
			thread_arg->cy_r54 = cy_r54->d0;
			thread_arg->cy_r55 = cy_r54->d1;
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
			thread_arg->cy_r28 = cy_r28;
			thread_arg->cy_r29 = cy_r29;
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
			thread_arg->cy_r40 = cy_r40;
			thread_arg->cy_r41 = cy_r41;
			thread_arg->cy_r42 = cy_r42;
			thread_arg->cy_r43 = cy_r43;
			thread_arg->cy_r44 = cy_r44;
			thread_arg->cy_r45 = cy_r45;
			thread_arg->cy_r46 = cy_r46;
			thread_arg->cy_r47 = cy_r47;
			thread_arg->cy_r48 = cy_r48;
			thread_arg->cy_r49 = cy_r49;
			thread_arg->cy_r50 = cy_r50;
			thread_arg->cy_r51 = cy_r51;
			thread_arg->cy_r52 = cy_r52;
			thread_arg->cy_r53 = cy_r53;
			thread_arg->cy_r54 = cy_r54;
			thread_arg->cy_r55 = cy_r55;

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
			thread_arg->cy_r28 = cy_r28->d0;	thread_arg->cy_i28 = cy_i28->d0;
			thread_arg->cy_r29 = cy_r28->d1;	thread_arg->cy_i29 = cy_i28->d1;
			thread_arg->cy_r30 = cy_r28->d2;	thread_arg->cy_i30 = cy_i28->d2;
			thread_arg->cy_r31 = cy_r28->d3;	thread_arg->cy_i31 = cy_i28->d3;
			thread_arg->cy_r32 = cy_r32->d0;	thread_arg->cy_i32 = cy_i32->d0;
			thread_arg->cy_r33 = cy_r32->d1;	thread_arg->cy_i33 = cy_i32->d1;
			thread_arg->cy_r34 = cy_r32->d2;	thread_arg->cy_i34 = cy_i32->d2;
			thread_arg->cy_r35 = cy_r32->d3;	thread_arg->cy_i35 = cy_i32->d3;
			thread_arg->cy_r36 = cy_r36->d0;	thread_arg->cy_i36 = cy_i36->d0;
			thread_arg->cy_r37 = cy_r36->d1;	thread_arg->cy_i37 = cy_i36->d1;
			thread_arg->cy_r38 = cy_r36->d2;	thread_arg->cy_i38 = cy_i36->d2;
			thread_arg->cy_r39 = cy_r36->d3;	thread_arg->cy_i39 = cy_i36->d3;
			thread_arg->cy_r40 = cy_r40->d0;	thread_arg->cy_i40 = cy_i40->d0;
			thread_arg->cy_r41 = cy_r40->d1;	thread_arg->cy_i41 = cy_i40->d1;
			thread_arg->cy_r42 = cy_r40->d2;	thread_arg->cy_i42 = cy_i40->d2;
			thread_arg->cy_r43 = cy_r40->d3;	thread_arg->cy_i43 = cy_i40->d3;
			thread_arg->cy_r44 = cy_r44->d0;	thread_arg->cy_i44 = cy_i44->d0;
			thread_arg->cy_r45 = cy_r44->d1;	thread_arg->cy_i45 = cy_i44->d1;
			thread_arg->cy_r46 = cy_r44->d2;	thread_arg->cy_i46 = cy_i44->d2;
			thread_arg->cy_r47 = cy_r44->d3;	thread_arg->cy_i47 = cy_i44->d3;
			thread_arg->cy_r48 = cy_r48->d0;	thread_arg->cy_i48 = cy_i48->d0;
			thread_arg->cy_r49 = cy_r48->d1;	thread_arg->cy_i49 = cy_i48->d1;
			thread_arg->cy_r50 = cy_r48->d2;	thread_arg->cy_i50 = cy_i48->d2;
			thread_arg->cy_r51 = cy_r48->d3;	thread_arg->cy_i51 = cy_i48->d3;
			thread_arg->cy_r52 = cy_r52->d0;	thread_arg->cy_i52 = cy_i52->d0;
			thread_arg->cy_r53 = cy_r52->d1;	thread_arg->cy_i53 = cy_i52->d1;
			thread_arg->cy_r54 = cy_r52->d2;	thread_arg->cy_i54 = cy_i52->d2;
			thread_arg->cy_r55 = cy_r52->d3;	thread_arg->cy_i55 = cy_i52->d3;
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
			thread_arg->cy_r14 = cy_r28->d0;	thread_arg->cy_i14 = cy_r28->d1;
			thread_arg->cy_r15 = cy_r30->d0;	thread_arg->cy_i15 = cy_r30->d1;
			thread_arg->cy_r16 = cy_r32->d0;	thread_arg->cy_i16 = cy_r32->d1;
			thread_arg->cy_r17 = cy_r34->d0;	thread_arg->cy_i17 = cy_r34->d1;
			thread_arg->cy_r18 = cy_r36->d0;	thread_arg->cy_i18 = cy_r36->d1;
			thread_arg->cy_r19 = cy_r38->d0;	thread_arg->cy_i19 = cy_r38->d1;
			thread_arg->cy_r20 = cy_r40->d0;	thread_arg->cy_i20 = cy_r40->d1;
			thread_arg->cy_r21 = cy_r42->d0;	thread_arg->cy_i21 = cy_r42->d1;
			thread_arg->cy_r22 = cy_r44->d0;	thread_arg->cy_i22 = cy_r44->d1;
			thread_arg->cy_r23 = cy_r46->d0;	thread_arg->cy_i23 = cy_r46->d1;
			thread_arg->cy_r24 = cy_r48->d0;	thread_arg->cy_i24 = cy_r48->d1;
			thread_arg->cy_r25 = cy_r50->d0;	thread_arg->cy_i25 = cy_r50->d1;
			thread_arg->cy_r26 = cy_r52->d0;	thread_arg->cy_i26 = cy_r52->d1;
			thread_arg->cy_r27 = cy_r54->d0;	thread_arg->cy_i27 = cy_r54->d1;
			thread_arg->cy_r28 = cy_i00->d0;	thread_arg->cy_i28 = cy_i00->d1;
			thread_arg->cy_r29 = cy_i02->d0;	thread_arg->cy_i29 = cy_i02->d1;
			thread_arg->cy_r30 = cy_i04->d0;	thread_arg->cy_i30 = cy_i04->d1;
			thread_arg->cy_r31 = cy_i06->d0;	thread_arg->cy_i31 = cy_i06->d1;
			thread_arg->cy_r32 = cy_i08->d0;	thread_arg->cy_i32 = cy_i08->d1;
			thread_arg->cy_r33 = cy_i10->d0;	thread_arg->cy_i33 = cy_i10->d1;
			thread_arg->cy_r34 = cy_i12->d0;	thread_arg->cy_i34 = cy_i12->d1;
			thread_arg->cy_r35 = cy_i14->d0;	thread_arg->cy_i35 = cy_i14->d1;
			thread_arg->cy_r36 = cy_i16->d0;	thread_arg->cy_i36 = cy_i16->d1;
			thread_arg->cy_r37 = cy_i18->d0;	thread_arg->cy_i37 = cy_i18->d1;
			thread_arg->cy_r38 = cy_i20->d0;	thread_arg->cy_i38 = cy_i20->d1;
			thread_arg->cy_r39 = cy_i22->d0;	thread_arg->cy_i39 = cy_i22->d1;
			thread_arg->cy_r40 = cy_i24->d0;	thread_arg->cy_i40 = cy_i24->d1;
			thread_arg->cy_r41 = cy_i26->d0;	thread_arg->cy_i41 = cy_i26->d1;
			thread_arg->cy_r42 = cy_i28->d0;	thread_arg->cy_i42 = cy_i28->d1;
			thread_arg->cy_r43 = cy_i30->d0;	thread_arg->cy_i43 = cy_i30->d1;
			thread_arg->cy_r44 = cy_i32->d0;	thread_arg->cy_i44 = cy_i32->d1;
			thread_arg->cy_r45 = cy_i34->d0;	thread_arg->cy_i45 = cy_i34->d1;
			thread_arg->cy_r46 = cy_i36->d0;	thread_arg->cy_i46 = cy_i36->d1;
			thread_arg->cy_r47 = cy_i38->d0;	thread_arg->cy_i47 = cy_i38->d1;
			thread_arg->cy_r48 = cy_i40->d0;	thread_arg->cy_i48 = cy_i40->d1;
			thread_arg->cy_r49 = cy_i42->d0;	thread_arg->cy_i49 = cy_i42->d1;
			thread_arg->cy_r50 = cy_i44->d0;	thread_arg->cy_i50 = cy_i44->d1;
			thread_arg->cy_r51 = cy_i46->d0;	thread_arg->cy_i51 = cy_i46->d1;
			thread_arg->cy_r52 = cy_i48->d0;	thread_arg->cy_i52 = cy_i48->d1;
			thread_arg->cy_r53 = cy_i50->d0;	thread_arg->cy_i53 = cy_i50->d1;
			thread_arg->cy_r54 = cy_i52->d0;	thread_arg->cy_i54 = cy_i52->d1;
			thread_arg->cy_r55 = cy_i54->d0;	thread_arg->cy_i55 = cy_i54->d1;
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
			thread_arg->cy_r28 = cy_r28;		thread_arg->cy_i28 = cy_i28;
			thread_arg->cy_r29 = cy_r29;		thread_arg->cy_i29 = cy_i29;
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
			thread_arg->cy_r40 = cy_r40;		thread_arg->cy_i40 = cy_i40;
			thread_arg->cy_r41 = cy_r41;		thread_arg->cy_i41 = cy_i41;
			thread_arg->cy_r42 = cy_r42;		thread_arg->cy_i42 = cy_i42;
			thread_arg->cy_r43 = cy_r43;		thread_arg->cy_i43 = cy_i43;
			thread_arg->cy_r44 = cy_r44;		thread_arg->cy_i44 = cy_i44;
			thread_arg->cy_r45 = cy_r45;		thread_arg->cy_i45 = cy_i45;
			thread_arg->cy_r46 = cy_r46;		thread_arg->cy_i46 = cy_i46;
			thread_arg->cy_r47 = cy_r47;		thread_arg->cy_i47 = cy_i47;
			thread_arg->cy_r48 = cy_r48;		thread_arg->cy_i48 = cy_i48;
			thread_arg->cy_r49 = cy_r49;		thread_arg->cy_i49 = cy_i49;
			thread_arg->cy_r50 = cy_r50;		thread_arg->cy_i50 = cy_i50;
			thread_arg->cy_r51 = cy_r51;		thread_arg->cy_i51 = cy_i51;
			thread_arg->cy_r52 = cy_r52;		thread_arg->cy_i52 = cy_i52;
			thread_arg->cy_r53 = cy_r53;		thread_arg->cy_i53 = cy_i53;
			thread_arg->cy_r54 = cy_r54;		thread_arg->cy_i54 = cy_i54;
			thread_arg->cy_r55 = cy_r55;		thread_arg->cy_i55 = cy_i55;

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

