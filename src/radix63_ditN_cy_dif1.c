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

#ifdef MULTITHREAD
	// No parallel impl of this experimental code yet! Please recompile without the USE_THREADS option.
	#undef MULTITHREAD
	#undef USE_PTHREAD
#endif	// MULTITHREAD

#define RADIX 63	// Use #define rather than const int to ensure it's really a compile-time const in the C sense
#define ODD_RADIX 63	// ODD_RADIX = [radix >> trailz(radix)]

/***************/

int radix63_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 pexp)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-63 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-63 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix63_ditN_cy_dif1";
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
	int NDIVR,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,ntmp,outer;

	// Need these both in scalar mode and to ease the SSE2-array init...dimension = ODD_RADIX;
	// In order to ease the ptr-access for the || routine, lump these 4*ODD_RADIX doubles together with copies of
	// the 4 in the passed-in bs[2] and bsinv[2] arrays [and used in this 4-double form by the mersenne-mod carry macros]
	// into a single foo_array[4*(ODD_RADIX+1)], then convert what used to be disparate ODD_RADIX-sized arrays to pointers.
	static double foo_array[(ODD_RADIX+1)<<2], *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr, bs,bsinv;

	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p[RADIX], nsave = 0;
#ifndef MULTITHREAD
	const double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	const double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
	int k0,k1,k2,k3,k4,k5,k6,k7,k8;
	const uint8 *iptr,
		dif_iperm[64] = {	// Only need 63, but pad to 64-bytes
			0x00,0x36,0x2d,0x24,0x1b,0x12,0x09,
			0x38,0x2f,0x26,0x1d,0x14,0x0b,0x02,
			0x31,0x28,0x1f,0x16,0x0d,0x04,0x3a,
			0x2a,0x21,0x18,0x0f,0x06,0x3c,0x33,
			0x23,0x1a,0x11,0x08,0x3e,0x35,0x2c,
			0x1c,0x13,0x0a,0x01,0x37,0x2e,0x25,
			0x15,0x0c,0x03,0x39,0x30,0x27,0x1e,
			0x0e,0x05,0x3b,0x32,0x29,0x20,0x17,
			0x07,0x3d,0x34,0x2b,0x22,0x19,0x10},
		dif_operm[64] = {	// ditto
			0x00,0x07,0x03,0x02,0x06,0x05,0x01,0x08,0x04,
			0x37,0x3e,0x3a,0x36,0x3d,0x39,0x38,0x3c,0x3b,
			0x32,0x2e,0x35,0x31,0x2d,0x34,0x30,0x2f,0x33,
			0x2a,0x29,0x25,0x2c,0x28,0x24,0x2b,0x27,0x26,
			0x1d,0x21,0x20,0x1c,0x23,0x1f,0x1b,0x22,0x1e,
			0x15,0x14,0x18,0x17,0x13,0x1a,0x16,0x12,0x19,
			0x10,0x0c,0x0b,0x0f,0x0e,0x0a,0x11,0x0d,0x09},
		dit_iperm[64] = {	// Only need 63, but pad to 64-bytes
			0x00,0x02,0x01,0x08,0x07,0x06,0x05,0x04,0x03,
			0x32,0x31,0x30,0x2f,0x2e,0x2d,0x34,0x33,0x35,
			0x1d,0x1c,0x1b,0x22,0x21,0x23,0x1f,0x1e,0x20,
			0x10,0x0f,0x11,0x0d,0x0c,0x0e,0x0a,0x09,0x0b,
			0x37,0x36,0x38,0x3c,0x3e,0x3d,0x39,0x3b,0x3a,
			0x2a,0x2c,0x2b,0x27,0x29,0x28,0x24,0x26,0x25,
			0x15,0x17,0x16,0x12,0x14,0x13,0x1a,0x19,0x18},
		dit_operm[64] = {	// ditto
			0x00,0x24,0x09,0x2d,0x12,0x36,0x1b,
			0x0e,0x32,0x17,0x3b,0x20,0x05,0x29,
			0x1c,0x01,0x25,0x0a,0x2e,0x13,0x37,
			0x2a,0x0f,0x33,0x18,0x3c,0x21,0x06,
			0x38,0x1d,0x02,0x26,0x0b,0x2f,0x14,
			0x07,0x2b,0x10,0x34,0x19,0x3d,0x22,
			0x15,0x39,0x1e,0x03,0x27,0x0c,0x30,
			0x23,0x08,0x2c,0x11,0x35,0x1a,0x3e,
			0x31,0x16,0x3a,0x1f,0x04,0x28,0x0d};
	double re,im,rt,it, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13;
#endif
	static double radix_inv, n2inv;
	double scale,dtmp, maxerr = 0.0;
	// Local storage: We must use an array here because scalars have no guarantees about relative address offsets
	// [and even if those are contiguous-as-hoped-for, they may run in reverse]; Make array type (struct complex)
	// to allow us to use the same offset-indexing as in the original radix-32 in-place DFT macros:
	double *addr, *addi;
	struct complex t[RADIX], *tptr;
	int *itmp;	// Pointer into the bjmodn array
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	static int n_div_nwt;
	int col,co2,co3,m,m2;
  #if 0//def USE_AVX
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
  #endif
	double wt_re,wt_im, wi_re,wi_im;	// Fermat-mod/LOACC weights stuff, used in both scalar and SIMD mode
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
	static int wts_idx_incr = 0, icycle[ODD_RADIX],ic;

#ifdef MULTITHREAD

	static struct cy_thread_data_t *tdat = 0x0;
	// Threadpool-based dispatch stuff:
	static int main_work_units = 0, pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)cy63_process_chunk, NULL, 0x0};

#else

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
	static double *_maxerr = 0x0, *_cy_r[RADIX],*_cy_i[RADIX];
	if(!_maxerr) {
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
	if((TEST_TYPE & 0xfffffffe) == TEST_TYPE_PRP) {	// Mask off low bit to lump together PRP and PRP-C tests
		i = (iter-1) % ITERS_BETWEEN_CHECKPOINTS;	// Bit we need to read...iter-counter is unit-offset w.r.to iter-interval, hence the -1
		if((BASE_MULTIPLIER_BITS[i>>6] >> (i&63)) & 1)
			prp_mult = PRP_BASE;
	}

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/RADIX;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		fprintf(stderr,"ERROR: iter = %10d; NWT_BITS does not divide N/RADIX in %s.\n",iter,func);
		err = ERR_SKIP_RADIX_SET;
		return(err);
	}

	if(pexp != psave || n != nsave) {	/* Exponent or array length change triggers re-init */
		first_entry=TRUE;
	}

/*...initialize things upon first entry: */

	if(first_entry)
	{
		psave = pexp;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)RADIX));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = pexp%n;	/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
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
				if( ((uint32)&tdat[ithread].cy_dat[l] & SZ_VDM1) == 0 ) {
					tdat[ithread].cy_r = &tdat[ithread].cy_dat[l];
					tdat[ithread].cy_i = tdat[ithread].cy_r + RADIX;
				//	fprintf(stderr,"%d-byte-align cy_dat array at element[%d]\n",SZ_VD,l);
					break;
				}
			}
			ASSERT(l < 4, "Failed to align cy_dat array!");
		}
	#endif

		/*   constant index offsets for array load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		// Legacy code needs 3 lowest nonzero fixed-index p[] terms:
		p1 = 1*NDIVR;	 p1 += ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = 2*NDIVR;	 p2 += ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = 3*NDIVR;	 p3 += ( (p3 >> DAT_BITS) << PAD_BITS );
		for(l = 0; l < RADIX; l++) {
			p[l] = l*NDIVR;
			p[l] += ( (p[l] >> DAT_BITS) << PAD_BITS );
		}

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
			_cy_r[i]	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r[i]== 0x0);
			_cy_i[i]	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i[i]== 0x0);
		}
		_maxerr	= (double *)malloc(j);	ptr_prod += (uint32)(_maxerr== 0x0);

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

			/* Subtract nwt from the increments to ease fast-mod */
			wts_idx_incr -= nwt;

			for(i = 0; i < ODD_RADIX; i++) {
				/* Need this both in scalar mode and to ease the SSE2-array init */
				j = _bjmodn[i][0] > sw;	bs_arr[i] = base[j];	bsinv_arr[i] = baseinv[j];
				wt_arr[i] = wt0[ii[i]];	// inverse wts must be reinited on each pass, since these have a *scale multiplier
				/* Give icycle indices their proper starting values: */
				icycle[i] = i;
			}
		}

	#ifdef USE_PTHREAD

		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			tdat[ithread].bjmodnini = _bjmodnini[CY_THREADS];
			tdat[ithread].bjmodn0 = _bjmodnini[ithread];
			// In scalar mode use these 2 ptrs to pass wts_idx_incr and the base/baseinv/etc array-ptrs:
			tdat[ithread].r00      = (vec_dbl *)foo_array;
			tdat[ithread].half_arr = (vec_dbl *)&wts_idx_incr;
		}

		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		{
			// These inits must occur just once, in loop-simulated full-pass mode,
			// in order to get the array-index-offset values of the icycle/jcycle indices right:
			for(i = 0; i < ODD_RADIX; i++) {
				tdat[0].icycle[i] = icycle[i];
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
						icycle[i] += wts_idx_incr;		icycle[i] += ( (-(int)((uint32)icycle[i] >> 31)) & nwt);
					}
				}
				for(i = 0; i < ODD_RADIX; i++) {
					tdat[ithread].icycle[i] = icycle[i];
				}
			}
			// Restore the original loop-start values of the cycle arrays, since we use these for init of inv-wts below:
			for(i = 0; i < ODD_RADIX; i++) {
				icycle[i] = tdat[0].icycle[i];
			}
		}

	#endif	// USE_PTHREAD

		first_entry=FALSE;
	}	/* endif(first_entry) */

/*...The radix-63 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(i = 0; i < RADIX; i++) {
			_cy_r[i][ithread] = 0;
			_cy_i[i][ithread] = 0;
		}
	}
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy_r[0][0] = -2;
	}

	*fracmax=0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	// init inverse-weight scale factor = 2/n for normal carry pass, 1 for wrapper pass

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
		}
	#endif
	}	// 	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)

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
		tdat[ithread].maxerr = _maxerr[ithread];
		tdat[ithread].scale = scale;
		tdat[ithread].prp_mult = prp_mult;

	// pointer data:
		tdat[ithread].arrdat = a;			/* Main data array */
		ASSERT(tdat[ithread].wt0 == wt0, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].wt1 == wt1, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].si  == si, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].rn0 == rn0, "thread-local memcheck fail!");
		ASSERT(tdat[ithread].rn1 == rn1, "thread-local memcheck fail!");
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* init carries	*/
			for(i = 0; i < RADIX; i++) {
				tdat[ithread].cy_r[i] = _cy_r[i][ithread];
			}
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
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
		/***** DEC/HP CC doesn't properly copy init value of maxerr = 0 into threads,
		so need to set once again explicitly for each: *****/
		maxerr = 0.0;
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
			for(l = 0; l < RADIX; l++) {
				cy_r[l] = _cy_r[l][ithread];
			}
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
			/* init carries	*/
			for(l = 0; l < RADIX; l++) {
				cy_r[l] = _cy_r[l][ithread];		cy_i[l] = _cy_i[l][ithread];
			}
		}

		/********************************************************************************/
		/* This main loop is same for un-and-multithreaded, so stick into a header file */
		/* (can't use a macro because of the #if-enclosed stuff).                       */
		/********************************************************************************/
		#include "radix63_main_carry_loop.h"

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			for(l = 0; l < RADIX; l++) {
				_cy_r[l][ithread] = cy_r[l];
			}
		}
		else
		{
			for(l = 0; l < RADIX; l++) {
				_cy_r[l][ithread] = cy_r[l];		_cy_i[l][ithread] = cy_i[l];
			}
		}

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
		ASSERT(0x0 == cy63_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
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
		_maxerr[ithread] = tdat[ithread].maxerr;
		if(maxerr < _maxerr[ithread]) {
			maxerr = _maxerr[ithread];
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
			// Radix == 3 (mod 4) here, so change 4-loop structure to handle that:
			for(l = 0; l < RADIX>>2; l++) {
				jt = j1 + p[l<<2];	// p0,4,8,...
				a[jt   ] *= radix_inv;
				a[jt+p1] *= radix_inv;
				a[jt+p2] *= radix_inv;
				a[jt+p3] *= radix_inv;
			}
			// Cleanuo loop for leftovers (mod 4):
			for(l = (l<<2); l < RADIX; l++) {
				a[j1+p[l]] *= radix_inv;
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
		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
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

void radix63_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-63 complex DIF FFT pass on the data in the length-N real vector A.
*/
	int j,j1,j2, l,k0,k1,k2,k3,k4,k5,k6,k7,k8;
	const uint8 *iptr,
		dif_iperm[64] = {	// Only need 63, but pad to 64-bytes
			0x00,0x36,0x2d,0x24,0x1b,0x12,0x09,
			0x38,0x2f,0x26,0x1d,0x14,0x0b,0x02,
			0x31,0x28,0x1f,0x16,0x0d,0x04,0x3a,
			0x2a,0x21,0x18,0x0f,0x06,0x3c,0x33,
			0x23,0x1a,0x11,0x08,0x3e,0x35,0x2c,
			0x1c,0x13,0x0a,0x01,0x37,0x2e,0x25,
			0x15,0x0c,0x03,0x39,0x30,0x27,0x1e,
			0x0e,0x05,0x3b,0x32,0x29,0x20,0x17,
			0x07,0x3d,0x34,0x2b,0x22,0x19,0x10},
		dif_operm[64] = {	// ditto
			0x00,0x07,0x03,0x02,0x06,0x05,0x01,0x08,0x04,
			0x37,0x3e,0x3a,0x36,0x3d,0x39,0x38,0x3c,0x3b,
			0x32,0x2e,0x35,0x31,0x2d,0x34,0x30,0x2f,0x33,
			0x2a,0x29,0x25,0x2c,0x28,0x24,0x2b,0x27,0x26,
			0x1d,0x21,0x20,0x1c,0x23,0x1f,0x1b,0x22,0x1e,
			0x15,0x14,0x18,0x17,0x13,0x1a,0x16,0x12,0x19,
			0x10,0x0c,0x0b,0x0f,0x0e,0x0a,0x11,0x0d,0x09};
	// p-indexing is hexadecimal here:
	static int NDIVR,p[RADIX], first_entry=TRUE;
	const double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	const double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
	double re,im,rt,it, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13;
	struct complex t[RADIX], *tptr;

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR = n/RADIX;
	/*   constant index offsets for array load/stores are here.	*/
		for(l = 0; l < RADIX; l++) {
			p[l] = l*NDIVR;
			p[l] += ( (p[l] >> DAT_BITS) << PAD_BITS );
		}
	}

/*...The radix-63 pass is here.	*/

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

	//...gather the needed data (63 64-bit complex, i.e. 126 64-bit reals) and do 9 radix-7 transforms:
	/*
	Twiddleless version arranges 9 sets of radix-7 DFT inputs as follows: 0 in upper left corner,
	decrement (mod 63) 9 horizontally and 7 vertically. Display result of DIF/DIT input-scramble array in hex:

		00,36,2d,24,1b,12,09
		38,2f,26,1d,14,0b,02
		31,28,1f,16,0d,04,3a
		2a,21,18,0f,06,3c,33
		23,1a,11,08,3e,35,2c
		1c,13,0a,01,37,2e,25
		15,0c,03,39,30,27,1e
		0e,05,3b,32,29,20,17
		07,3d,34,2b,22,19,10
	*/
	#if 1	// Prefer compact-obj-code scheme:
		tptr = t; iptr = dif_iperm;
		for(l = 0; l < 9; l++) {
			k0 = p[*iptr]; k1 = p[*(iptr+1)]; k2 = p[*(iptr+2)]; k3 = p[*(iptr+3)]; k4 = p[*(iptr+4)]; k5 = p[*(iptr+5)]; k6 = p[*(iptr+6)];
			RADIX_07_DFT(
				a[j1+k0],a[j2+k0],a[j1+k1],a[j2+k1],a[j1+k2],a[j2+k2],a[j1+k3],a[j2+k3],a[j1+k4],a[j2+k4],a[j1+k5],a[j2+k5],a[j1+k6],a[j2+k6],
				t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,
				tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im,
				uc1,us1,uc2,us2,uc3,us3, rt,it,re,im
			);	tptr++; iptr += 7;
		}
	#else
		tptr = t;
		RADIX_07_DFT(a[j1+p[0x00]],a[j2+p[0x00]],a[j1+p[0x36]],a[j2+p[0x36]],a[j1+p[0x2d]],a[j2+p[0x2d]],a[j1+p[0x24]],a[j2+p[0x24]],a[j1+p[0x1b]],a[j2+p[0x1b]],a[j1+p[0x12]],a[j2+p[0x12]],a[j1+p[0x09]],a[j2+p[0x09]], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		RADIX_07_DFT(a[j1+p[0x38]],a[j2+p[0x38]],a[j1+p[0x2f]],a[j2+p[0x2f]],a[j1+p[0x26]],a[j2+p[0x26]],a[j1+p[0x1d]],a[j2+p[0x1d]],a[j1+p[0x14]],a[j2+p[0x14]],a[j1+p[0x0b]],a[j2+p[0x0b]],a[j1+p[0x02]],a[j2+p[0x02]], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		RADIX_07_DFT(a[j1+p[0x31]],a[j2+p[0x31]],a[j1+p[0x28]],a[j2+p[0x28]],a[j1+p[0x1f]],a[j2+p[0x1f]],a[j1+p[0x16]],a[j2+p[0x16]],a[j1+p[0x0d]],a[j2+p[0x0d]],a[j1+p[0x04]],a[j2+p[0x04]],a[j1+p[0x3a]],a[j2+p[0x3a]], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		RADIX_07_DFT(a[j1+p[0x2a]],a[j2+p[0x2a]],a[j1+p[0x21]],a[j2+p[0x21]],a[j1+p[0x18]],a[j2+p[0x18]],a[j1+p[0x0f]],a[j2+p[0x0f]],a[j1+p[0x06]],a[j2+p[0x06]],a[j1+p[0x3c]],a[j2+p[0x3c]],a[j1+p[0x33]],a[j2+p[0x33]], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		RADIX_07_DFT(a[j1+p[0x23]],a[j2+p[0x23]],a[j1+p[0x1a]],a[j2+p[0x1a]],a[j1+p[0x11]],a[j2+p[0x11]],a[j1+p[0x08]],a[j2+p[0x08]],a[j1+p[0x3e]],a[j2+p[0x3e]],a[j1+p[0x35]],a[j2+p[0x35]],a[j1+p[0x2c]],a[j2+p[0x2c]], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		RADIX_07_DFT(a[j1+p[0x1c]],a[j2+p[0x1c]],a[j1+p[0x13]],a[j2+p[0x13]],a[j1+p[0x0a]],a[j2+p[0x0a]],a[j1+p[0x01]],a[j2+p[0x01]],a[j1+p[0x37]],a[j2+p[0x37]],a[j1+p[0x2e]],a[j2+p[0x2e]],a[j1+p[0x25]],a[j2+p[0x25]], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		RADIX_07_DFT(a[j1+p[0x15]],a[j2+p[0x15]],a[j1+p[0x0c]],a[j2+p[0x0c]],a[j1+p[0x03]],a[j2+p[0x03]],a[j1+p[0x39]],a[j2+p[0x39]],a[j1+p[0x30]],a[j2+p[0x30]],a[j1+p[0x27]],a[j2+p[0x27]],a[j1+p[0x1e]],a[j2+p[0x1e]], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		RADIX_07_DFT(a[j1+p[0x0e]],a[j2+p[0x0e]],a[j1+p[0x05]],a[j2+p[0x05]],a[j1+p[0x3b]],a[j2+p[0x3b]],a[j1+p[0x32]],a[j2+p[0x32]],a[j1+p[0x29]],a[j2+p[0x29]],a[j1+p[0x20]],a[j2+p[0x20]],a[j1+p[0x17]],a[j2+p[0x17]], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		RADIX_07_DFT(a[j1+p[0x07]],a[j2+p[0x07]],a[j1+p[0x3d]],a[j2+p[0x3d]],a[j1+p[0x34]],a[j2+p[0x34]],a[j1+p[0x2b]],a[j2+p[0x2b]],a[j1+p[0x22]],a[j2+p[0x22]],a[j1+p[0x19]],a[j2+p[0x19]],a[j1+p[0x10]],a[j2+p[0x10]], t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);
	#endif
	/*...and now do 7 radix-9 transforms. The required output permutation is

		00,07,03,02,06,05,01,08,04,
		37,3e,3a,36,3d,39,38,3c,3b,
		32,2e,35,31,2d,34,30,2f,33,
		2a,29,25,2c,28,24,2b,27,26,
		1d,21,20,1c,23,1f,1b,22,1e,
		15,14,18,17,13,1a,16,12,19,
		10,0c,0b,0f,0e,0a,11,0d,09.
	*/
	#if 1	// Prefer compact-obj-code scheme:
		tptr = t; iptr = dif_operm;
		for(l = 0; l < 7; l++) {
			k0 = p[*iptr]; k1 = p[*(iptr+1)]; k2 = p[*(iptr+2)]; k3 = p[*(iptr+3)]; k4 = p[*(iptr+4)]; k5 = p[*(iptr+5)]; k6 = p[*(iptr+6)]; k7 = p[*(iptr+7)]; k8 = p[*(iptr+8)];
			RADIX_09_DIF(
				tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im,
				a[j1+k0],a[j2+k0],a[j1+k1],a[j2+k1],a[j1+k2],a[j2+k2],a[j1+k3],a[j2+k3],a[j1+k4],a[j2+k4],a[j1+k5],a[j2+k5],a[j1+k6],a[j2+k6],a[j1+k7],a[j2+k7],a[j1+k8],a[j2+k8],
				rt,it,re
			);	tptr += 9; iptr += 9;
		}
	#else
		tptr = t;
		RADIX_09_DIF( tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im, a[j1+p[0x00]],a[j2+p[0x00]],a[j1+p[0x07]],a[j2+p[0x07]],a[j1+p[0x03]],a[j2+p[0x03]],a[j1+p[0x02]],a[j2+p[0x02]],a[j1+p[0x06]],a[j2+p[0x06]],a[j1+p[0x05]],a[j2+p[0x05]],a[j1+p[0x01]],a[j2+p[0x01]],a[j1+p[0x08]],a[j2+p[0x08]],a[j1+p[0x04]],a[j2+p[0x04]], rt,it,re);	tptr += 9;
		RADIX_09_DIF( tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im, a[j1+p[0x37]],a[j2+p[0x37]],a[j1+p[0x3e]],a[j2+p[0x3e]],a[j1+p[0x3a]],a[j2+p[0x3a]],a[j1+p[0x36]],a[j2+p[0x36]],a[j1+p[0x3d]],a[j2+p[0x3d]],a[j1+p[0x39]],a[j2+p[0x39]],a[j1+p[0x38]],a[j2+p[0x38]],a[j1+p[0x3c]],a[j2+p[0x3c]],a[j1+p[0x3b]],a[j2+p[0x3b]], rt,it,re);	tptr += 9;
		RADIX_09_DIF( tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im, a[j1+p[0x32]],a[j2+p[0x32]],a[j1+p[0x2e]],a[j2+p[0x2e]],a[j1+p[0x35]],a[j2+p[0x35]],a[j1+p[0x31]],a[j2+p[0x31]],a[j1+p[0x2d]],a[j2+p[0x2d]],a[j1+p[0x34]],a[j2+p[0x34]],a[j1+p[0x30]],a[j2+p[0x30]],a[j1+p[0x2f]],a[j2+p[0x2f]],a[j1+p[0x33]],a[j2+p[0x33]], rt,it,re);	tptr += 9;
		RADIX_09_DIF( tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im, a[j1+p[0x2a]],a[j2+p[0x2a]],a[j1+p[0x29]],a[j2+p[0x29]],a[j1+p[0x25]],a[j2+p[0x25]],a[j1+p[0x2c]],a[j2+p[0x2c]],a[j1+p[0x28]],a[j2+p[0x28]],a[j1+p[0x24]],a[j2+p[0x24]],a[j1+p[0x2b]],a[j2+p[0x2b]],a[j1+p[0x27]],a[j2+p[0x27]],a[j1+p[0x26]],a[j2+p[0x26]], rt,it,re);	tptr += 9;
		RADIX_09_DIF( tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im, a[j1+p[0x1d]],a[j2+p[0x1d]],a[j1+p[0x21]],a[j2+p[0x21]],a[j1+p[0x20]],a[j2+p[0x20]],a[j1+p[0x1c]],a[j2+p[0x1c]],a[j1+p[0x23]],a[j2+p[0x23]],a[j1+p[0x1f]],a[j2+p[0x1f]],a[j1+p[0x1b]],a[j2+p[0x1b]],a[j1+p[0x22]],a[j2+p[0x22]],a[j1+p[0x1e]],a[j2+p[0x1e]], rt,it,re);	tptr += 9;
		RADIX_09_DIF( tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im, a[j1+p[0x15]],a[j2+p[0x15]],a[j1+p[0x14]],a[j2+p[0x14]],a[j1+p[0x18]],a[j2+p[0x18]],a[j1+p[0x17]],a[j2+p[0x17]],a[j1+p[0x13]],a[j2+p[0x13]],a[j1+p[0x1a]],a[j2+p[0x1a]],a[j1+p[0x16]],a[j2+p[0x16]],a[j1+p[0x12]],a[j2+p[0x12]],a[j1+p[0x19]],a[j2+p[0x19]], rt,it,re);	tptr += 9;
		RADIX_09_DIF( tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im, a[j1+p[0x10]],a[j2+p[0x10]],a[j1+p[0x0c]],a[j2+p[0x0c]],a[j1+p[0x0b]],a[j2+p[0x0b]],a[j1+p[0x0f]],a[j2+p[0x0f]],a[j1+p[0x0e]],a[j2+p[0x0e]],a[j1+p[0x0a]],a[j2+p[0x0a]],a[j1+p[0x11]],a[j2+p[0x11]],a[j1+p[0x0d]],a[j2+p[0x0d]],a[j1+p[0x09]],a[j2+p[0x09]], rt,it,re);
	#endif
	}
}

/***************/

void radix63_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-63 complex DIT FFT pass on the data in the length-N real vector A.
*/
	int j,j1,j2, l,k0,k1,k2,k3,k4,k5,k6,k7,k8;
	const uint8 *iptr,
		dit_iperm[64] = {	// Only need 63, but pad to 64-bytes
			0x00,0x02,0x01,0x08,0x07,0x06,0x05,0x04,0x03,
			0x32,0x31,0x30,0x2f,0x2e,0x2d,0x34,0x33,0x35,
			0x1d,0x1c,0x1b,0x22,0x21,0x23,0x1f,0x1e,0x20,
			0x10,0x0f,0x11,0x0d,0x0c,0x0e,0x0a,0x09,0x0b,
			0x37,0x36,0x38,0x3c,0x3e,0x3d,0x39,0x3b,0x3a,
			0x2a,0x2c,0x2b,0x27,0x29,0x28,0x24,0x26,0x25,
			0x15,0x17,0x16,0x12,0x14,0x13,0x1a,0x19,0x18},
		dit_operm[64] = {	// ditto
			0x00,0x24,0x09,0x2d,0x12,0x36,0x1b,
			0x0e,0x32,0x17,0x3b,0x20,0x05,0x29,
			0x1c,0x01,0x25,0x0a,0x2e,0x13,0x37,
			0x2a,0x0f,0x33,0x18,0x3c,0x21,0x06,
			0x38,0x1d,0x02,0x26,0x0b,0x2f,0x14,
			0x07,0x2b,0x10,0x34,0x19,0x3d,0x22,
			0x15,0x39,0x1e,0x03,0x27,0x0c,0x30,
			0x23,0x08,0x2c,0x11,0x35,0x1a,0x3e,
			0x31,0x16,0x3a,0x1f,0x04,0x28,0x0d};
	// p-indexing is hexadecimal here:
	static int NDIVR,p[RADIX], first_entry=TRUE;
	const double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	const double	c   =  0.76604444311897803520,	/* cos(2*pi/9) */
					s   =  0.64278760968653932631,	/* sin(2*pi/9) */
					c2  =  0.17364817766693034887,	/* cos(2*u) */
					s2  =  0.98480775301220805936,	/* sin(2*u) */
					c3m1= -1.50000000000000000000,	/* cos(3*u)-1 */
					s3  =  0.86602540378443864677,	/* sin(3*u) */
					c4  = -0.93969262078590838404,	/* cos(4*u) */
					s4  =  0.34202014332566873307;	/* sin(4*u) */
	double re,im,rt,it, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13;
	struct complex t[RADIX], *tptr;

	if(!first_entry && (n/RADIX) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR = n/RADIX;
	/*   constant index offsets for array load/stores are here.	*/
		for(l = 0; l < RADIX; l++) {
			p[l] = l*NDIVR;
			p[l] += ( (p[l] >> DAT_BITS) << PAD_BITS );
		}
	}

/*...The radix-63 pass is here.	*/

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

	//...gather the needed data (63 64-bit complex, i.e. 126 64-bit reals) and do 7 radix-9 transforms:
	/*
	Twiddleless version arranges 7 sets of radix-9 DFT inputs as follows: 0 in upper left corner,
	decrement (mod 63) 9 horizontally and 7 vertically. Applying a further bit-reversal to that,
	Display result of Combined DIT input-scramble array in hex:

		00,02,01,08,07,06,05,04,03,
		32,31,30,2f,2e,2d,34,33,35,
		1d,1c,1b,22,21,23,1f,1e,20,
		10,0f,11,0d,0c,0e,0a,09,0b,
		37,36,38,3c,3e,3d,39,3b,3a,
		2a,2c,2b,27,29,28,24,26,25,
		15,17,16,12,14,13,1a,19,18.
	*/
	#if 1	// Prefer compact-obj-code scheme:
		tptr = t; iptr = dit_iperm;
		for(l = 0; l < 7; l++) {
			k0 = p[*iptr]; k1 = p[*(iptr+1)]; k2 = p[*(iptr+2)]; k3 = p[*(iptr+3)]; k4 = p[*(iptr+4)]; k5 = p[*(iptr+5)]; k6 = p[*(iptr+6)]; k7 = p[*(iptr+7)]; k8 = p[*(iptr+8)];
			RADIX_09_DIT(
				a[j1+k0],a[j2+k0],a[j1+k1],a[j2+k1],a[j1+k2],a[j2+k2],a[j1+k3],a[j2+k3],a[j1+k4],a[j2+k4],a[j1+k5],a[j2+k5],a[j1+k6],a[j2+k6],a[j1+k7],a[j2+k7],a[j1+k8],a[j2+k8],
				tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im,
				rt,it,re
			);	tptr += 9; iptr += 9;
		}
	#else
		tptr = t;
		RADIX_09_DIT(a[j1+p[0x00]],a[j2+p[0x00]],a[j1+p[0x02]],a[j2+p[0x02]],a[j1+p[0x01]],a[j2+p[0x01]],a[j1+p[0x08]],a[j2+p[0x08]],a[j1+p[0x07]],a[j2+p[0x07]],a[j1+p[0x06]],a[j2+p[0x06]],a[j1+p[0x05]],a[j2+p[0x05]],a[j1+p[0x04]],a[j2+p[0x04]],a[j1+p[0x03]],a[j2+p[0x03]], tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im, rt,it,re);	tptr += 9;
		RADIX_09_DIT(a[j1+p[0x32]],a[j2+p[0x32]],a[j1+p[0x31]],a[j2+p[0x31]],a[j1+p[0x30]],a[j2+p[0x30]],a[j1+p[0x2f]],a[j2+p[0x2f]],a[j1+p[0x2e]],a[j2+p[0x2e]],a[j1+p[0x2d]],a[j2+p[0x2d]],a[j1+p[0x34]],a[j2+p[0x34]],a[j1+p[0x33]],a[j2+p[0x33]],a[j1+p[0x35]],a[j2+p[0x35]], tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im, rt,it,re);	tptr += 9;
		RADIX_09_DIT(a[j1+p[0x1d]],a[j2+p[0x1d]],a[j1+p[0x1c]],a[j2+p[0x1c]],a[j1+p[0x1b]],a[j2+p[0x1b]],a[j1+p[0x22]],a[j2+p[0x22]],a[j1+p[0x21]],a[j2+p[0x21]],a[j1+p[0x23]],a[j2+p[0x23]],a[j1+p[0x1f]],a[j2+p[0x1f]],a[j1+p[0x1e]],a[j2+p[0x1e]],a[j1+p[0x20]],a[j2+p[0x20]], tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im, rt,it,re);	tptr += 9;
		RADIX_09_DIT(a[j1+p[0x10]],a[j2+p[0x10]],a[j1+p[0x0f]],a[j2+p[0x0f]],a[j1+p[0x11]],a[j2+p[0x11]],a[j1+p[0x0d]],a[j2+p[0x0d]],a[j1+p[0x0c]],a[j2+p[0x0c]],a[j1+p[0x0e]],a[j2+p[0x0e]],a[j1+p[0x0a]],a[j2+p[0x0a]],a[j1+p[0x09]],a[j2+p[0x09]],a[j1+p[0x0b]],a[j2+p[0x0b]], tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im, rt,it,re);	tptr += 9;
		RADIX_09_DIT(a[j1+p[0x37]],a[j2+p[0x37]],a[j1+p[0x36]],a[j2+p[0x36]],a[j1+p[0x38]],a[j2+p[0x38]],a[j1+p[0x3c]],a[j2+p[0x3c]],a[j1+p[0x3e]],a[j2+p[0x3e]],a[j1+p[0x3d]],a[j2+p[0x3d]],a[j1+p[0x39]],a[j2+p[0x39]],a[j1+p[0x3b]],a[j2+p[0x3b]],a[j1+p[0x3a]],a[j2+p[0x3a]], tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im, rt,it,re);	tptr += 9;
		RADIX_09_DIT(a[j1+p[0x2a]],a[j2+p[0x2a]],a[j1+p[0x2c]],a[j2+p[0x2c]],a[j1+p[0x2b]],a[j2+p[0x2b]],a[j1+p[0x27]],a[j2+p[0x27]],a[j1+p[0x29]],a[j2+p[0x29]],a[j1+p[0x28]],a[j2+p[0x28]],a[j1+p[0x24]],a[j2+p[0x24]],a[j1+p[0x26]],a[j2+p[0x26]],a[j1+p[0x25]],a[j2+p[0x25]], tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im, rt,it,re);	tptr += 9;
		RADIX_09_DIT(a[j1+p[0x15]],a[j2+p[0x15]],a[j1+p[0x17]],a[j2+p[0x17]],a[j1+p[0x16]],a[j2+p[0x16]],a[j1+p[0x12]],a[j2+p[0x12]],a[j1+p[0x14]],a[j2+p[0x14]],a[j1+p[0x13]],a[j2+p[0x13]],a[j1+p[0x1a]],a[j2+p[0x1a]],a[j1+p[0x19]],a[j2+p[0x19]],a[j1+p[0x18]],a[j2+p[0x18]], tptr->re,tptr->im,(tptr+1)->re,(tptr+1)->im,(tptr+2)->re,(tptr+2)->im,(tptr+3)->re,(tptr+3)->im,(tptr+4)->re,(tptr+4)->im,(tptr+5)->re,(tptr+5)->im,(tptr+6)->re,(tptr+6)->im,(tptr+7)->re,(tptr+7)->im,(tptr+8)->re,(tptr+8)->im, rt,it,re);
	#endif
	/*...and now do 9 radix-7 transforms. The required output permutation is

		00,24,09,2d,12,36,1b,
		0e,32,17,3b,20,05,29,
		1c,01,25,0a,2e,13,37,
		2a,0f,33,18,3c,21,06,
		38,1d,02,26,0b,2f,14,
		07,2b,10,34,19,3d,22,
		15,39,1e,03,27,0c,30,
		23,08,2c,11,35,1a,3e,
		31,16,3a,1f,04,28,0d
	*/
	#if 1	// Prefer compact-obj-code scheme:
		tptr = t; iptr = dit_operm;
		for(l = 0; l < 9; l++) {
			k0 = p[*iptr]; k1 = p[*(iptr+1)]; k2 = p[*(iptr+2)]; k3 = p[*(iptr+3)]; k4 = p[*(iptr+4)]; k5 = p[*(iptr+5)]; k6 = p[*(iptr+6)];
			RADIX_07_DFT(
				tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im,
				t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,
				a[j1+k0],a[j2+k0],a[j1+k1],a[j2+k1],a[j1+k2],a[j2+k2],a[j1+k3],a[j2+k3],a[j1+k4],a[j2+k4],a[j1+k5],a[j2+k5],a[j1+k6],a[j2+k6],
				uc1,us1,uc2,us2,uc3,us3, rt,it,re,im
			);	tptr++; iptr += 7;
		}
	#else
		tptr = t;
		RADIX_07_DFT(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[j1+p[0x00]],a[j2+p[0x00]],a[j1+p[0x24]],a[j2+p[0x24]],a[j1+p[0x09]],a[j2+p[0x09]],a[j1+p[0x2d]],a[j2+p[0x2d]],a[j1+p[0x12]],a[j2+p[0x12]],a[j1+p[0x36]],a[j2+p[0x36]],a[j1+p[0x1b]],a[j2+p[0x1b]], uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		RADIX_07_DFT(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[j1+p[0x0e]],a[j2+p[0x0e]],a[j1+p[0x32]],a[j2+p[0x32]],a[j1+p[0x17]],a[j2+p[0x17]],a[j1+p[0x3b]],a[j2+p[0x3b]],a[j1+p[0x20]],a[j2+p[0x20]],a[j1+p[0x05]],a[j2+p[0x05]],a[j1+p[0x29]],a[j2+p[0x29]], uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		RADIX_07_DFT(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[j1+p[0x1c]],a[j2+p[0x1c]],a[j1+p[0x01]],a[j2+p[0x01]],a[j1+p[0x25]],a[j2+p[0x25]],a[j1+p[0x0a]],a[j2+p[0x0a]],a[j1+p[0x2e]],a[j2+p[0x2e]],a[j1+p[0x13]],a[j2+p[0x13]],a[j1+p[0x37]],a[j2+p[0x37]], uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		RADIX_07_DFT(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[j1+p[0x2a]],a[j2+p[0x2a]],a[j1+p[0x0f]],a[j2+p[0x0f]],a[j1+p[0x33]],a[j2+p[0x33]],a[j1+p[0x18]],a[j2+p[0x18]],a[j1+p[0x3c]],a[j2+p[0x3c]],a[j1+p[0x21]],a[j2+p[0x21]],a[j1+p[0x06]],a[j2+p[0x06]], uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		RADIX_07_DFT(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[j1+p[0x38]],a[j2+p[0x38]],a[j1+p[0x1d]],a[j2+p[0x1d]],a[j1+p[0x02]],a[j2+p[0x02]],a[j1+p[0x26]],a[j2+p[0x26]],a[j1+p[0x0b]],a[j2+p[0x0b]],a[j1+p[0x2f]],a[j2+p[0x2f]],a[j1+p[0x14]],a[j2+p[0x14]], uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		RADIX_07_DFT(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[j1+p[0x07]],a[j2+p[0x07]],a[j1+p[0x2b]],a[j2+p[0x2b]],a[j1+p[0x10]],a[j2+p[0x10]],a[j1+p[0x34]],a[j2+p[0x34]],a[j1+p[0x19]],a[j2+p[0x19]],a[j1+p[0x3d]],a[j2+p[0x3d]],a[j1+p[0x22]],a[j2+p[0x22]], uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		RADIX_07_DFT(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[j1+p[0x15]],a[j2+p[0x15]],a[j1+p[0x39]],a[j2+p[0x39]],a[j1+p[0x1e]],a[j2+p[0x1e]],a[j1+p[0x03]],a[j2+p[0x03]],a[j1+p[0x27]],a[j2+p[0x27]],a[j1+p[0x0c]],a[j2+p[0x0c]],a[j1+p[0x30]],a[j2+p[0x30]], uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		RADIX_07_DFT(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[j1+p[0x23]],a[j2+p[0x23]],a[j1+p[0x08]],a[j2+p[0x08]],a[j1+p[0x2c]],a[j2+p[0x2c]],a[j1+p[0x11]],a[j2+p[0x11]],a[j1+p[0x35]],a[j2+p[0x35]],a[j1+p[0x1a]],a[j2+p[0x1a]],a[j1+p[0x3e]],a[j2+p[0x3e]], uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);	tptr++;
		RADIX_07_DFT(tptr->re,tptr->im,(tptr+9)->re,(tptr+9)->im,(tptr+18)->re,(tptr+18)->im,(tptr+27)->re,(tptr+27)->im,(tptr+36)->re,(tptr+36)->im,(tptr+45)->re,(tptr+45)->im,(tptr+54)->re,(tptr+54)->im, t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13, a[j1+p[0x31]],a[j2+p[0x31]],a[j1+p[0x16]],a[j2+p[0x16]],a[j1+p[0x3a]],a[j2+p[0x3a]],a[j1+p[0x1f]],a[j2+p[0x1f]],a[j1+p[0x04]],a[j2+p[0x04]],a[j1+p[0x28]],a[j2+p[0x28]],a[j1+p[0x0d]],a[j2+p[0x0d]], uc1,us1,uc2,us2,uc3,us3, rt,it,re,im);
	#endif
	}
}

#undef RADIX
#undef ODD_RADIX

