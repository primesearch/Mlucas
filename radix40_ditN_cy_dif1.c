/*******************************************************************************
*                                                                              *
*   (C) 1997-2014 by Ernst W. Mayer.                                           *
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

#define RADIX 40	// Use #define rather than const int to ensure it's really a compile-time const in the C sense

#ifdef MULTITHREAD
	#ifndef USE_PTHREAD
		#error Pthreads is only thread model currently supported!
	#endif
#endif

#ifdef USE_SSE2

	#define ERR_CHECK_ALL	/* #define this to do ROE checking of all convolution outputs, rather than just every 36th one */
	#if defined(USE_AVX) && !defined(ERR_CHECK_ALL)
		#error ERR_CHECK_ALL *required* for AVX-mode builds!
	#endif

	#define EPS 1e-10

  // For Mersenne-mod we need (16 [SSE2] or 64 [AVX]) + 4 added slots for the half_arr lookup tables.
  // Add relevant number (half_arr_offset40 + RADIX) to get required value of radix40_creals_in_local_store:
  #ifdef USE_AVX
	const int half_arr_offset40 = 178;	// + RADIX = 218; Used for thread local-storage-integrity checking
	const int radix40_creals_in_local_store = 288;	// += 68 and round up to nearest multiple of 4
  #else
	const int half_arr_offset40 = 188;	// + RADIX = 228; Used for thread local-storage-integrity checking
	const int radix40_creals_in_local_store = 248;	// += 20 and round up to nearest multiple of 4
  #endif

	#include "sse2_macro.h"

	#ifdef COMPILER_TYPE_MSVC

	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 64

		  // #define to either (if left undefined) use small-macro form below, or (if defined) to inline the fused macros as single big blob of asm (64-bit only):
		  #define GCC_ASM_FULL_INLINE	1	// Enable xmm0-15-using version of radix-40 DFTs

		  #if GCC_ASM_FULL_INLINE
			#include "radix40_ditN_cy_dif1_gcc64.h"
		  #endif

		#endif

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

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
		int *si;
	#ifdef USE_SSE2
		vec_dbl *r00;
		vec_dbl *half_arr;
	#else
		double *r00;
		double *half_arr;
	#endif
		uint32 bjmodnini;
		int bjmodn0;
	// For large radix0 use thread-local arrays for DWT indices/carries - only caveat is these must be SIMD-aligned:
	#if GCC_EVER_GETS_ITS_ACT_TOGETHER_HERE
	/* Jan 2014: Bloody hell - turns out GCC uses __BIGGEST_ALIGNMENT__ = 16 on x86, which is too small to be useful for avx data!
		int bjmodn[768] __attribute__ ((aligned (32)));
		double cy[768] __attribute__ ((aligned (32)));
	*/
	#else
	// Thus, we are forced to resort to fugly hackage - add pad slots to a garbage-named struct-internal array along with
	// a pointer-to-be-inited-at-runtime, when we set ptr to the lowest-index array element having the desired alginment:
		double *cy;
		double cy_dat[RADIX+4] __attribute__ ((__aligned__(8)));	// Enforce min-alignment of 8 bytes in 32-bit builds.
	#endif
	};

#endif

/***************/

int radix40_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-40 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-40 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix40_ditN_cy_dif1";
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
#ifdef USE_SSE2
	const int sz_vd = sizeof(vec_dbl), sz_vd_m1 = sz_vd-1;
	// lg(sizeof(vec_dbl)):
  #ifdef USE_AVX
	const int l2_sz_vd = 5;
  #else
	const int l2_sz_vd = 4;
  #endif
#else
	const int sz_vd = sizeof(double), sz_vd_m1 = sz_vd-1;
#endif

	int NDIVR,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,ntmp,outer,nbytes;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p16,p24,p32;
	static int poff[RADIX>>2];	// Store mults of p04 offset for loop control
#if defined(USE_SSE2) || !defined(MULTITHREAD)
	const double cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
				 cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
				 s2  =  0.95105651629515357211,	/*  sin(u) */
				 ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
				 ss2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
#endif
	static double radix_inv, n2inv;
	double scale, dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	struct complex t[RADIX], *tptr;
	int *itmp;	// Pointer into the bjmodn array
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3;
  #ifdef USE_AVX
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
  #endif

#ifdef USE_SSE2

  #if !(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))
	#error SSE2 code not supported for this compiler!
  #endif

	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else
	double *add0, *add1, *add2, *add3, *add4, *add5, *add6, *add7;	/* Addresses into array sections */
  #endif

	static int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
	const double crnd = 3.0*0x4000000*0x2000000;
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	vec_dbl *tmp,*tm1,*tm2;	// Non-static utility ptrs
	static vec_dbl *isrt2, *xcc1, *xss1, *xcc2, *xss2, *xss3, *max_err, *sse2_rnd, *half_arr
		,*r00,*r02,*r04,*r06,*r08,*r10,*r12,*r14,*r16,*r18,*r20,*r22,*r24,*r26,*r28,*r30,*r32,*r34,*r36,*r38
		,*r40,*r42,*r44,*r46,*r48,*r50,*r52,*r54,*r56,*r58,*r60,*r62,*r64,*r66,*r68,*r70,*r72,*r74,*r76,*r78
		,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r
		,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r,*s1p28r,*s1p29r,*s1p30r,*s1p31r,*s1p32r,*s1p33r,*s1p34r,*s1p35r,*s1p36r,*s1p37r,*s1p38r,*s1p39r
		,*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx

#endif

#ifdef MULTITHREAD

	static struct cy_thread_data_t *tdat = 0x0;
	// Threadpool-based dispatch stuff:
	static int main_work_units = 0, pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)cy40_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	int m,m2;
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	double *addr;
	int bjmodn[RADIX];
	double rt,it,temp,frac,cy[RADIX]
		,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_bjmodnini = 0x0,*_bjmodn[RADIX];
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_maxerr = 0x0,*_cy[RADIX];
	if(!_maxerr) {
		_cy[0] = 0x0;	// First of these used as an "already inited consts?" sentinel, must init = 0x0 at same time do so for non-array static ptrs
	}

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, 0, "Fermat-mod only available for radices 7,8,9,15 and their multiples!");
	}

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
		if(CY_THREADS < NTHREADS)	{ WARN(HERE, "CY_THREADS < NTHREADS", "", 1); return(ERR_ASSERT); }
		if(!isPow2(CY_THREADS))		{ WARN(HERE, "CY_THREADS not a power of 2!", "", 1); return(ERR_ASSERT); }
		if(CY_THREADS > 1)
		{
			if(NDIVR    %CY_THREADS != 0) { WARN(HERE, "NDIVR    %CY_THREADS != 0", "", 1); return(ERR_ASSERT); }
			if(n_div_nwt%CY_THREADS != 0) { WARN(HERE, "n_div_nwt%CY_THREADS != 0", "", 1); return(ERR_ASSERT); }
		}

	  #ifdef USE_PTHREAD

		j = (uint32)sizeof(struct cy_thread_data_t);
		tdat = (struct cy_thread_data_t *)calloc(CY_THREADS, j);

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

		// This array pointer must be set based on vec_dbl-sized alignment at runtime for each thread:
			for(l = 0; l < 4; l++) {
				if( ((uint32)&tdat[ithread].cy_dat[l] & sz_vd_m1) == 0 ) {
					tdat[ithread].cy = &tdat[ithread].cy_dat[l];
				//	fprintf(stderr,"%d-byte-align cy_dat array at element[%d]\n",sz_vd,l);
					break;
				}
			}
			ASSERT(HERE, l < 4, "Failed to align cy_dat array!");
		}
	#endif

	#ifdef USE_SSE2

		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use double-complex type size (16 bytes) to alloc a block of local storage
		// consisting of 88 dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix40_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix40_creals_in_local_store);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 80 16-byte slots of sc_arr for temporaries, next 5 for the nontrivial complex 16th roots,
	next 80 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	#ifdef USE_PTHREAD
		__r0 = sc_ptr;
	#endif
									tmp	= sc_ptr + 0x50;
		r00	= sc_ptr + 0x00;		s1p00r = tmp + 0x00;
		r02	= sc_ptr + 0x02;		s1p01r = tmp + 0x02;
		r04	= sc_ptr + 0x04;		s1p02r = tmp + 0x04;
		r06	= sc_ptr + 0x06;		s1p03r = tmp + 0x06;
		r08	= sc_ptr + 0x08;		s1p04r = tmp + 0x08;
		r10	= sc_ptr + 0x0a;		s1p05r = tmp + 0x0a;
		r12	= sc_ptr + 0x0c;		s1p06r = tmp + 0x0c;
		r14	= sc_ptr + 0x0e;		s1p07r = tmp + 0x0e;
		r16	= sc_ptr + 0x10;		s1p08r = tmp + 0x10;
		r18	= sc_ptr + 0x12;		s1p09r = tmp + 0x12;
		r20	= sc_ptr + 0x14;		s1p10r = tmp + 0x14;
		r22	= sc_ptr + 0x16;		s1p11r = tmp + 0x16;
		r24	= sc_ptr + 0x18;		s1p12r = tmp + 0x18;
		r26	= sc_ptr + 0x1a;		s1p13r = tmp + 0x1a;
		r28	= sc_ptr + 0x1c;		s1p14r = tmp + 0x1c;
		r30	= sc_ptr + 0x1e;		s1p15r = tmp + 0x1e;
		r32	= sc_ptr + 0x20;		s1p16r = tmp + 0x20;
		r34	= sc_ptr + 0x22;		s1p17r = tmp + 0x22;
		r36	= sc_ptr + 0x24;		s1p18r = tmp + 0x24;
		r38	= sc_ptr + 0x26;		s1p19r = tmp + 0x26;
		r40	= sc_ptr + 0x28;		s1p20r = tmp + 0x28;
		r42	= sc_ptr + 0x2a;		s1p21r = tmp + 0x2a;
		r44	= sc_ptr + 0x2c;		s1p22r = tmp + 0x2c;
		r46	= sc_ptr + 0x2e;		s1p23r = tmp + 0x2e;
		r48	= sc_ptr + 0x30;		s1p24r = tmp + 0x30;
		r50	= sc_ptr + 0x32;		s1p25r = tmp + 0x32;
		r52	= sc_ptr + 0x34;		s1p26r = tmp + 0x34;
		r54	= sc_ptr + 0x36;		s1p27r = tmp + 0x36;
		r56	= sc_ptr + 0x38;		s1p28r = tmp + 0x38;
		r58	= sc_ptr + 0x3a;		s1p29r = tmp + 0x3a;
		r60	= sc_ptr + 0x3c;		s1p30r = tmp + 0x3c;
		r62	= sc_ptr + 0x3e;		s1p31r = tmp + 0x3e;
		r64	= sc_ptr + 0x40;		s1p32r = tmp + 0x40;
		r66	= sc_ptr + 0x42;		s1p33r = tmp + 0x42;
		r68	= sc_ptr + 0x44;		s1p34r = tmp + 0x44;
		r70	= sc_ptr + 0x46;		s1p35r = tmp + 0x46;
		r72	= sc_ptr + 0x48;		s1p36r = tmp + 0x48;
		r74	= sc_ptr + 0x4a;		s1p37r = tmp + 0x4a;
		r76	= sc_ptr + 0x4c;		s1p38r = tmp + 0x4c;
		r78	= sc_ptr + 0x4e;		s1p39r = tmp + 0x4e;
	// sc_ptr += 160
		isrt2 = tmp + 0x50;
		xcc1 = tmp + 0x51;
		xcc2 = tmp + 0x52;
		xss1 = tmp + 0x53;
		xss2 = tmp + 0x54;
		xss3 = tmp + 0x55;
	// sc_ptr += 166
		tmp += 0x56;
	#ifdef USE_AVX
		cy = tmp;		tmp += 0xa;
	// sc_ptr += 176
	#else
		cy = tmp;		tmp += 0x14;
	// sc_ptr += 186
	#endif
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
	// sc_ptr += 178 [AVX] or 188 [SSE2]; This is where the value of half_arr_offset40 comes from
		half_arr= tmp + 0x02;	/* This table needs 20x16 bytes */

		ASSERT(HERE, (radix40_creals_in_local_store << l2_sz_vd) >= ((long)half_arr - (long)r00) + (20 << l2_sz_vd), "radix40_creals_in_local_store checksum failed!");

		/* These remain fixed: */
		VEC_DBL_INIT(isrt2, ISRT2);
		VEC_DBL_INIT(xcc1 , cc1  );
		VEC_DBL_INIT(xcc2 , cc2  );
		VEC_DBL_INIT(xss1 , s2   );
		VEC_DBL_INIT(xss2 , ss1  );
		VEC_DBL_INIT(xss3 , ss2  );

		VEC_DBL_INIT(sse2_rnd,crnd);

		// Propagate the above consts to the remaining threads:
		nbytes = (int)xss3 - (int)isrt2 + sz_vd;	// #bytes in above sincos block of data
		tmp = isrt2;
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

		nbytes = 16 << l2_sz_vd;
	#endif

		// Propagate the above consts to the remaining threads:
		tmp = half_arr;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
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

	// For large radices, array-access to bjmodn means only init base-ptr here:
	#ifdef USE_AVX
		bjmodn = (int*)(sinwtm1 + RE_IM_STRIDE);
	#else
		bjmodn = (int*)(sse_n   + RE_IM_STRIDE);
	#endif

	#endif	// USE_SSE2

		/*   constant index offsets for array load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p16 = p08 + p08;
		p24 = p16 + p08;
		p32 = p24 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );

		poff[0x0] =   0; poff[0x1] = p04    ; poff[0x2] = p08; poff[0x3] = p04+p08;
		poff[0x4] = p16; poff[0x5] = p04+p16; poff[0x6] = p24; poff[0x7] = p04+p24;
		poff[0x8] = p32; poff[0x9] = p04+p32;

		if(_cy[0])	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;
			for(i = 0; i < RADIX; i++) {
				free((void *)_bjmodn[i]); _bjmodn[i] = 0x0;
				free((void *)    _cy[i]);     _cy[i] = 0x0;
			}
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
			_cy[i]	= (double *)malloc(j);	ptr_prod += (uint32)(_cy[i]== 0x0);
		}
		_maxerr	= (double *)malloc(j);	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/RADIX-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		_bjmodnini[0] = 0;
		_bjmodnini[1] = 0;

		jhi = NDIVR/CY_THREADS;

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

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			tdat[ithread].bjmodnini = _bjmodnini[CY_THREADS];
			tdat[ithread].bjmodn0 = _bjmodnini[ithread];
		#ifdef USE_SSE2
			tdat[ithread].r00 = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (long)tdat[ithread].r00 + ((long)half_arr - (long)r00);
		#else	// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
			tdat[ithread].r00      = (double *)base;
			tdat[ithread].half_arr = (double *)baseinv;
		#endif	// USE_SSE2
		}
	#endif

	}	/* endif(first_entry) */

/*...The radix-40 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(i = 0; i < RADIX; i++) {
			_cy[i][ithread] = 0;
		}
	}
	/* If an LL test, init the subtract-2: */
	if(TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy[0][0] = -2;
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
	_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. If n does not divide p, lowest-order digit is always a bigword (_i[0] = 1).	*/

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
	j = _bjmodnini[CY_THREADS];
	// Include 0-thread here ... bjmodn terms all 0 for that, but need jhi computed for all threads:
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_bjmodn[0][ithread] = _bjmodnini[ithread];
		for(i = 1; i < RADIX; i++) {
			MOD_ADD32(_bjmodn[i-1][ithread], j, n, _bjmodn[i][ithread]);
		}
		_jstart[ithread] = ithread*NDIVR/CY_THREADS;
		if(!full_pass)
			_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
		else
			_jhi[ithread] = _jstart[ithread] + nwt-1;

		_col[ithread] = ithread*(khi*RADIX);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
		_co2[ithread] = (n>>nwt_bits)-1+RADIX - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
		_co3[ithread] = _co2[ithread]-RADIX;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
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
		tdat[ithread].iter = iter;
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
	#ifdef USE_SSE2
		ASSERT(HERE, tdat[ithread].r00 == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].half_arr;
		ASSERT(HERE, ((tmp-1)->d0 == crnd && (tmp-1)->d1 == crnd), "thread-local memcheck failed!");
	  #ifdef USE_AVX
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif
	#endif
		/* init carries: */
		for(i = 0; i < RADIX; i++) {
			tdat[ithread].cy[i] = _cy[i][ithread];
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

		for(l = 0; l < RADIX; l++) {
			bjmodn[l] = _bjmodn[l][ithread];
		}
		/* init carries	*/
	#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
		tmp = cy;
		for(l = 0; l < RADIX; l += 4, ++tmp) {
			tmp->d0 = _cy[l  ][ithread];
			tmp->d1 = _cy[l+1][ithread];
			tmp->d2 = _cy[l+2][ithread];
			tmp->d3 = _cy[l+3][ithread];
		}
	#elif defined(USE_SSE2)
		tmp = cy;
		for(l = 0; l < RADIX; l += 2, ++tmp) {
			tmp->d0 = _cy[l  ][ithread];
			tmp->d1 = _cy[l+1][ithread];
		}
	#else
		for(l = 0; l < RADIX; l++) {
			cy[l] = _cy[l][ithread];
		}
	#endif

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix40_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
	#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
		tmp = cy;
		for(l = 0; l < RADIX; l += 4, ++tmp) {
			_cy[l  ][ithread] = tmp->d0;
			_cy[l+1][ithread] = tmp->d1;
			_cy[l+2][ithread] = tmp->d2;
			_cy[l+3][ithread] = tmp->d3;
		}
		maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
	#elif defined(USE_SSE2)
		tmp = cy;
		for(l = 0; l < RADIX; l += 2, ++tmp) {
			_cy[l  ][ithread] = tmp->d0;
			_cy[l+1][ithread] = tmp->d1;
		}
		maxerr = MAX(max_err->d0,max_err->d1);
	#else
		for(l = 0; l < RADIX; l++) {
			_cy[l][ithread] = cy[l];
		}
	#endif

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(_maxerr[ithread] < maxerr)
		{
			_maxerr[ithread] = maxerr;
		}

  #endif	// #ifdef USE_PTHREAD

	}	/******* END OF PARALLEL FOR-LOOP ********/

#ifdef USE_PTHREAD	// End of threadpool-based dispatch: Add a small wait-loop to ensure all threads complete

  #if 0//def OS_TYPE_MACOSX

	/*** Main execution thread executes remaining chunks in serial fashion (but in || with the pool threads): ***/
	for(j = 0; j < main_work_units; ++j)
	{
	//	printf("adding main task %d\n",j + pool_work_units);
		ASSERT(HERE, 0x0 == cy40_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
	ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
	ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(HERE, 0 == nanosleep(&ns_time, 0x0), "nanosleep fail!");
	}
//	printf("radix40_ditN_cy_dif1 end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = tdat[ithread].maxerr;
		if(maxerr < _maxerr[ithread]) {
			maxerr = _maxerr[ithread];
		}
		for(l = 0; l < RADIX; l++) {
			_cy[l][ithread] = tdat[ithread].cy[l];
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
	for(l = 0; l < RADIX; l++) {
		t[l].re = _cy[l][CY_THREADS - 1];
	}
	for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
	{
		for(l = 0; l < RADIX; l++) {
			_cy[l][ithread] = _cy[l][ithread-1];
		}
	}
	_cy[0][0] =+t[RADIX-1].re;	/* ...The wraparound carry is here: */
	for(l = 1; l < RADIX; l++) {
		_cy[l][0] = t[l-1].re;
	}

	full_pass = 0;
	scale = 1;
	j_jhi = 7;

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(j = ithread*pini; j <= ithread*pini + j_jhi; j++)
		{
			// Generate padded version of j, since prepadding pini is thread-count unsafe:
			j1 = j + ( (j >> DAT_BITS) << PAD_BITS );
			for(l = 0; l < RADIX>>2; l++) {
				jt = j1 + poff[l ];	// poff[] = p04,p08,...,p56
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
			dtmp += fabs(_cy[l][ithread]);
		}
		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
	}
	if(dtmp != 0.0)
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

void radix40_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-40 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in dft_macro.h for details on the radix-5,8 subtransforms.
*/
	int j,j1,j2,jp,jt;
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p16,p24,p32, first_entry=TRUE;
	const double cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
				cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
				s2  =  0.95105651629515357211,	/*  sin(u) */
				ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
				ss2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
	,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71;
	struct complex t[RADIX], *tptr;

	if(!first_entry && (n/40) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR=n/40;

/*   constant index offsets for array load/stores are here.	*/

		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p16 = p08 + p08;
		p24 = p16 + p08;
		p32 = p24 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-40 pass is here.	*/

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
	Twiddleless version arranges 8 sets of radix-5 DFT inputs as follows: 0 in upper left corner, decrement 8 horizontally and 5 vertically:

		RADIX_05_DFT(00,32,24,16,08)
		RADIX_05_DFT(35,27,19,11,03)
		RADIX_05_DFT(30,22,14,06,38)
		RADIX_05_DFT(25,17,09,01,33)
		RADIX_05_DFT(20,12,04,36,28)
		RADIX_05_DFT(15,07,39,31,23)
		RADIX_05_DFT(10,02,34,26,18)
		RADIX_05_DFT(05,37,29,21,13)

	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF] to properly permute the radix-8 DFT outputs.
	*/
	/*...gather the needed data (40 64-bit complex) and do 8 radix-5 transforms...*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, rt,it);	tptr += 5;
		jt = j1+p03; jp = j2+p03;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, rt,it);	tptr += 5;
		jt = j1+p06; jp = j2+p06;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, rt,it);	tptr += 5;
		jt = j1+p01; jp = j2+p01;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, rt,it);	tptr += 5;
		jt = j1+p04; jp = j2+p04;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, rt,it);	tptr += 5;
		jt = j1+p07; jp = j2+p07;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, rt,it);	tptr += 5;
		jt = j1+p02; jp = j2+p02;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt+p08],a[jp+p08],a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, rt,it);	tptr += 5;
		jt = j1+p05; jp = j2+p05;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2,a[jt    ],a[jp    ],a[jt+p32],a[jp+p32],a[jt+p24],a[jp+p24],a[jt+p16],a[jp+p16],a[jt+p08],a[jp+p08], tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, rt,it);

	/*...and now do 5 radix-8 transforms, ***SWAPPING THE ORDER OF THE FINAL 2*** to undo the last-2-outputs-swap in the RADIX_05_DFT macro:	*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],rt,it);	tptr++;
		jt = j1+p32; jp = j2+p32;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],rt,it);	tptr++;
		jt = j1+p24; jp = j2+p24;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	tptr++;
		jt = j1+p16; jp = j2+p16;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],rt,it);	tptr++;
		jt = j1+p08; jp = j2+p08;	RADIX_08_DIF(tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71,a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
	}
}

/***************/

void radix40_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in dft_macro.h for details on the radix-5,8 subtransforms.
*/
	int j,j1,j2,jp,jt;
	static int NDIVR,p01,p02,p03,p04,p05,p06,p07,p08,p16,p24,p32, first_entry=TRUE;
	const double cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
				cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
				s2  =  0.95105651629515357211,	/*  sin(u) */
				ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
				ss2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
	,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71;
	struct complex t[RADIX], *tptr;

	if(!first_entry && (n/40) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR=n/40;

/*   constant index offsets for array load/stores are here.	*/

		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p16 = p08 + p08;
		p24 = p16 + p08;
		p32 = p24 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-40 pass is here.	*/

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

		00,01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39

		00,32,24,16,08,35,27,19,11,03,30,22,14,06,38,25,17,09,01,33,20,12,04,36,28,15,07,39,31,23,10,02,34,26,18,05,37,29,21,13.	(*)

	Remember, inputs to DIT are bit-reversed, so

		a[00,01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39] contain [cf. bit-reversal index array output by test_fft_radix()]
		x[00,20,10,30,05,25,15,35,01,21,11,31,06,26,16,36,02,22,12,32,07,27,17,37,03,23,13,33,08,28,18,38,04,24,14,34,09,29,19,39], which get swapped [using the permutation (*)] to
		x[00,20,30,10,35,15,25,05,32,12,22,02,27,07,17,37,24,04,14,34,19,39,09,29,16,36,06,26,11,31,01,21,08,28,38,18,03,23,33,13], which means the a-indices get swapped as
		a[00,01,03,02,07,06,05,04|19,18,17,16,21,20,22,23|33,32,34,35,38,39,36,37|14,15,12,13,10,11,08,09|28,29,31,30,24,25,27,26]. These are the 5 octets going into the radix-8 DFTs.

		RADIX_08_DFT(00,01,03,02,07,06,05,04)
		RADIX_08_DFT(19,18,17,16,21,20,22,23)
		RADIX_08_DFT(33,32,34,35,38,39,36,37)
		RADIX_08_DFT(14,15,12,13,10,11,08,09)
		RADIX_08_DFT(28,29,31,30,24,25,27,26)

	Use the inverse supercalafragalistic Ancient Chinese Secret index-munging formula [iSACSIMPF] to properly permute the radix-5 DFT outputs.
	*/
	/*...gather the needed data (40 64-bit complex) and do 5 radix-8 transforms,	*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_08_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71, tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, rt,it);	tptr++;
		jt = j1+p16; jp = j2+p16;	RADIX_08_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p05],a[jp+p05],a[jt+p04],a[jp+p04],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71, tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, rt,it);	tptr++;
		jt = j1+p32; jp = j2+p32;	RADIX_08_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71, tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, rt,it);	tptr++;
		jt = j1+p08; jp = j2+p08;	RADIX_08_DIT(a[jt+p06],a[jp+p06],a[jt+p07],a[jp+p07],a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71, tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, rt,it);	tptr++;
		jt = j1+p24; jp = j2+p24;	RADIX_08_DIT(a[jt+p04],a[jp+p04],a[jt+p05],a[jp+p05],a[jt+p07],a[jp+p07],a[jt+p06],a[jp+p06],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71, tptr->re,tptr->im,(tptr+5)->re,(tptr+5)->im,(tptr+10)->re,(tptr+10)->im,(tptr+15)->re,(tptr+15)->im,(tptr+20)->re,(tptr+20)->im,(tptr+25)->re,(tptr+25)->im,(tptr+30)->re,(tptr+30)->im,(tptr+35)->re,(tptr+35)->im, rt,it)
	/*...and now do 8 radix-5 transforms.	*/
		tptr = t;
		jt = j1    ; jp = j2    ;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2, tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],rt,it);	tptr += 5;
		jt = j1+p01; jp = j2+p01;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2, tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],rt,it);	tptr += 5;
		jt = j1+p02; jp = j2+p02;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2, tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],rt,it);	tptr += 5;
		jt = j1+p03; jp = j2+p03;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2, tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],rt,it);	tptr += 5;
		jt = j1+p04; jp = j2+p04;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2, tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],rt,it);	tptr += 5;
		jt = j1+p05; jp = j2+p05;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2, tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],rt,it);	tptr += 5;
		jt = j1+p06; jp = j2+p06;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2, tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],a[jt+p08],a[jp+p08],rt,it);	tptr += 5;
		jt = j1+p07; jp = j2+p07;	RADIX_05_DFT(cc1,cc2,s2,ss1,ss2, tptr->re,tptr->im,(tptr+0x1)->re,(tptr+0x1)->im,(tptr+0x2)->re,(tptr+0x2)->im,(tptr+0x3)->re,(tptr+0x3)->im,(tptr+0x4)->re,(tptr+0x4)->im, a[jt+p08],a[jp+p08],a[jt+p24],a[jp+p24],a[jt    ],a[jp    ],a[jt+p16],a[jp+p16],a[jt+p32],a[jp+p32],rt,it)
	}
}

/******************** Multithreaded function body - NO STATIC VARS BELOW THIS POINT!: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void*
	cy40_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		struct cy_thread_data_t* thread_arg = targ;	// Move to top because scalar-mode carry pointers taken directly from it
		double *addr;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p01,p02,p03,p04,p05,p06,p07,p08,p16,p24,p32;
		int poff[RADIX>>2];
		int j,j1,j2,jt,jp,k,l,ntmp;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
	#ifdef USE_AVX
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#else
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	#endif

	#ifdef USE_SSE2

		const double crnd = 3.0*0x4000000*0x2000000;
		int *itmp;	// Pointer into the bjmodn array
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		double *add0, *add1, *add2, *add3;
		int *bjmodn;	// Alloc mem for this along with other 	SIMD stuff
	  #if defined(COMPILER_TYPE_MSVC) || !GCC_ASM_FULL_INLINE
		double *add4, *add5, *add6, *add7;
	  #endif
		vec_dbl *isrt2, *xcc1, *xss1, *xcc2, *xss2, *xss3, *max_err, *sse2_rnd, *half_arr
			,*r00,*r02,*r04,*r06,*r08,*r10,*r12,*r14,*r16,*r18,*r20,*r22,*r24,*r26,*r28,*r30,*r32,*r34,*r36,*r38
			,*r40,*r42,*r44,*r46,*r48,*r50,*r52,*r54,*r56,*r58,*r60,*r62,*r64,*r66,*r68,*r70,*r72,*r74,*r76,*r78
			,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r
			,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r,*s1p28r,*s1p29r,*s1p30r,*s1p31r,*s1p32r,*s1p33r,*s1p34r,*s1p35r,*s1p36r,*s1p37r,*s1p38r,*s1p39r
			,*cy;	// Need RADIX/2 slots for sse2 carries, RADIX/4 for avx
		vec_dbl *tmp,*tm1,*tm2;	// Non-static utility ptrs
		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else

		const double cc1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4; u = 2*pi/5 */
					 cc2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					 s2  =  0.95105651629515357211,	/*  sin(u) */
					 ss1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					 ss2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
		double *base, *baseinv;
		const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int m,m2;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		int bjmodn[RADIX];	// Thread only carries a base datum here, must alloc a local array for remaining values
		double *cy = thread_arg->cy, rt,it,temp,frac
			,x00,x01,x10,x11,x20,x21,x30,x31,x40,x41,x50,x51,x60,x61,x70,x71;
		struct complex t[RADIX], *tptr;
		int *itmp;	// Pointer into the bjmodn array

	#endif

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

	// double data:
		double maxerr = thread_arg->maxerr;
		double scale = thread_arg->scale;	int full_pass = scale < 0.5;

	// pointer data:
		double *a = thread_arg->arrdat;
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
		int *si = thread_arg->si;

		/*   constant index offsets for array load/stores are here.	*/
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p16 = p08 + p08;
		p24 = p16 + p08;
		p32 = p24 + p08;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );

		poff[0x0] =   0; poff[0x1] = p04    ; poff[0x2] = p08; poff[0x3] = p04+p08;
		poff[0x4] = p16; poff[0x5] = p04+p16; poff[0x6] = p24; poff[0x7] = p04+p24;
		poff[0x8] = p32; poff[0x9] = p04+p32;

	#ifdef USE_SSE2
		r00	= thread_arg->r00;		tmp = r00 + 0x50;
		r00	= r00 + 0x00;		s1p00r = tmp + 0x00;
		r02	= r00 + 0x02;		s1p01r = tmp + 0x02;
		r04	= r00 + 0x04;		s1p02r = tmp + 0x04;
		r06	= r00 + 0x06;		s1p03r = tmp + 0x06;
		r08	= r00 + 0x08;		s1p04r = tmp + 0x08;
		r10	= r00 + 0x0a;		s1p05r = tmp + 0x0a;
		r12	= r00 + 0x0c;		s1p06r = tmp + 0x0c;
		r14	= r00 + 0x0e;		s1p07r = tmp + 0x0e;
		r16	= r00 + 0x10;		s1p08r = tmp + 0x10;
		r18	= r00 + 0x12;		s1p09r = tmp + 0x12;
		r20	= r00 + 0x14;		s1p10r = tmp + 0x14;
		r22	= r00 + 0x16;		s1p11r = tmp + 0x16;
		r24	= r00 + 0x18;		s1p12r = tmp + 0x18;
		r26	= r00 + 0x1a;		s1p13r = tmp + 0x1a;
		r28	= r00 + 0x1c;		s1p14r = tmp + 0x1c;
		r30	= r00 + 0x1e;		s1p15r = tmp + 0x1e;
		r32	= r00 + 0x20;		s1p16r = tmp + 0x20;
		r34	= r00 + 0x22;		s1p17r = tmp + 0x22;
		r36	= r00 + 0x24;		s1p18r = tmp + 0x24;
		r38	= r00 + 0x26;		s1p19r = tmp + 0x26;
		r40	= r00 + 0x28;		s1p20r = tmp + 0x28;
		r42	= r00 + 0x2a;		s1p21r = tmp + 0x2a;
		r44	= r00 + 0x2c;		s1p22r = tmp + 0x2c;
		r46	= r00 + 0x2e;		s1p23r = tmp + 0x2e;
		r48	= r00 + 0x30;		s1p24r = tmp + 0x30;
		r50	= r00 + 0x32;		s1p25r = tmp + 0x32;
		r52	= r00 + 0x34;		s1p26r = tmp + 0x34;
		r54	= r00 + 0x36;		s1p27r = tmp + 0x36;
		r56	= r00 + 0x38;		s1p28r = tmp + 0x38;
		r58	= r00 + 0x3a;		s1p29r = tmp + 0x3a;
		r60	= r00 + 0x3c;		s1p30r = tmp + 0x3c;
		r62	= r00 + 0x3e;		s1p31r = tmp + 0x3e;
		r64	= r00 + 0x40;		s1p32r = tmp + 0x40;
		r66	= r00 + 0x42;		s1p33r = tmp + 0x42;
		r68	= r00 + 0x44;		s1p34r = tmp + 0x44;
		r70	= r00 + 0x46;		s1p35r = tmp + 0x46;
		r72	= r00 + 0x48;		s1p36r = tmp + 0x48;
		r74	= r00 + 0x4a;		s1p37r = tmp + 0x4a;
		r76	= r00 + 0x4c;		s1p38r = tmp + 0x4c;
		r78	= r00 + 0x4e;		s1p39r = tmp + 0x4e;

		isrt2 = tmp + 0x50;
		xcc1 = tmp + 0x51;
		xcc2 = tmp + 0x52;
		xss1 = tmp + 0x53;
		xss2 = tmp + 0x54;
		xss3 = tmp + 0x55;
		tmp += 0x56;
	  #ifdef USE_AVX
		cy = tmp;		tmp += 0xa;
	  #else
		cy = tmp;		tmp += 0x14;
	  #endif
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;	/* This table needs 20x16 bytes */

		ASSERT(HERE, (r00 == thread_arg->r00), "thread-local memcheck failed!");
		ASSERT(HERE, (half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
		ASSERT(HERE, (sse2_rnd->d0 == crnd && sse2_rnd->d1 == crnd), "thread-local memcheck failed!");
		tmp = half_arr;
	  #ifdef USE_AVX
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif

		VEC_DBL_INIT(max_err, 0.0);

		sign_mask = (uint64*)(r00 + radix40_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (  #doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		sse_sw  = sse_bw    + RE_IM_STRIDE;
		sse_n   = sse_sw    + RE_IM_STRIDE;
	  #ifdef USE_AVX
		n_minus_sil   = (struct uint32x4 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_n + 2;
		sinwt         = (struct uint32x4 *)sse_n + 3;
		sinwtm1       = (struct uint32x4 *)sse_n + 4;

		bjmodn = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn = (int*)(sse_n + RE_IM_STRIDE);
	  #endif

	#else

		// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
		base    = (double *)thread_arg->r00     ;
		baseinv = (double *)thread_arg->half_arr;

	#endif	// USE_SSE2 ?

		/* Init DWT-indices: */
		uint32 bjmodnini = thread_arg->bjmodnini;
		bjmodn[0] = thread_arg->bjmodn0;
		for(l = 1; l < RADIX; l++) {	// must use e.g. l for loop idx here as i is used for dwt indexing
			MOD_ADD32(bjmodn[l-1], bjmodnini, n, bjmodn[l]);
		}

		/* init carries	*/
		addr = thread_arg->cy;
	#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
		tmp = cy;
		for(l = 0; l < RADIX; l += 4, ++tmp) {
			tmp->d0 = *(addr+l  );
			tmp->d1 = *(addr+l+1);
			tmp->d2 = *(addr+l+2);
			tmp->d3 = *(addr+l+3);
		}
	#elif defined(USE_SSE2)
		tmp = cy;
		for(l = 0; l < RADIX; l += 2, ++tmp) {
			tmp->d0 = *(addr+l  );
			tmp->d1 = *(addr+l+1);
		}
	#elif 0	// No_op in scalar case, since carry pattern matches that of thread data
		for(l = 0; l < RADIX; l++) {
			cy[l] = *(addr+l);
		}
	#endif

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix40_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		addr = thread_arg->cy;
	#ifdef USE_AVX
		tmp = cy;
		for(l = 0; l < RADIX; l += 4, ++tmp) {
			*(addr+l  ) = tmp->d0;
			*(addr+l+1) = tmp->d1;
			*(addr+l+2) = tmp->d2;
			*(addr+l+3) = tmp->d3;
		}
		maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
	#elif defined(USE_SSE2)
		tmp = cy;
		for(l = 0; l < RADIX; l += 2, ++tmp) {
			*(addr+l  ) = tmp->d0;
			*(addr+l+1) = tmp->d1;
		}
		maxerr = MAX(max_err->d0,max_err->d1);
	#elif 0	// No_op in scalar case, since carry pattern matches that of thread data
		for(l = 0; l < RADIX; l++) {
			*(addr+l) = cy[l];
		}
	#endif

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(thread_arg->maxerr < maxerr)
		{
			thread_arg->maxerr = maxerr;
		}

		return 0x0;
	}
#endif

#undef RADIX
