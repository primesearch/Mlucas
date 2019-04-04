/*******************************************************************************
*                                                                              *
*   (C) 1997-2018 by Ernst W. Mayer.                                           *
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

#undef RTIME
#undef CTIME

#ifdef MULTITHREAD
	#define	DBG_THREADS	0	/* Turn on to collect stats about how much work done by each thread */
	#define RTIME	/* In multithreaded mode, need to use real (wall-clock) time */

  #ifdef USE_PTHREAD

	#define USE_THREADPOOL	// Undef to use non-pooled simple spawn/rejoin thread-team model
	#ifdef USE_THREADPOOL
		#include "threadpool.h"
	#endif

	struct ferm_thread_data_t{
		int tid;
		int *retval;
		double *arrdat;			/* Main data array */
		int *arr_scratch;
		int n;					/* Chunksize */
		struct complex *rt0;	/* Roots table 1 */
		struct complex *rt1;	/* Roots table 2 */
		int *index;				/* Bit-reversal index array */
		int nradices_prim;
		int *radix_prim;
	#ifdef DBG_TIME	// define at compile time to enable internal timing diagnostics
		double dt_fwd;
		double dt_inv;
		double dt_sqr;
		double dt_tot;
	#endif
	};
  #endif

#else
	#define CTIME	/* In single-thread mode, prefer cycle-based time because of its finer granularity */
#endif

/*********************************************************************************/
/* Globals. Unless specified otherwise, these are declared in Mdata.h:           */
/*********************************************************************************/

uint32 SW_DIV_N;	/* Needed for the subset of radix* carry routines which support Fermat-mod. */

/***************/

/*...Subroutine to perform Fermat-mod squaring on the data in the length-N real vector A.

     Acronym: DIF = Decimation In Frequency, DIT = Decimation In Time.

     One-pass combined fwd-DWT-weighting/fwd-DIF-FFT/pointwise-squaring/inv-DIT-FFT/inv-DWT-weighting
     for use with a packed-vector complex transform to halve the runlength of the corresponding real-vector
     transform. transform. We use the standard weighting that allows an acyclic length-N convolution to be
     effected via a length-N cyclic, namely premultiplying the elements by weight factors which are the
     first N (2*N)th complex roots of unity (the array name 'wn' is a mnemonic for 'roots negacyclic'):

                   rn(j) = exp(i*j*pi/N), j=0,...,N-1,

     performing a length-N cyclic convolution on the weighted vector, then multiplying the
     outputs by the analogous inverse weights.

     We can exploit the symmetries of complex exponentials to reduce storage and improve
     performance in several ways here:

     Symmetry (1): Since each component of WT is simply a complex exponential, inverse weights
                   are just complex conjugates, and don't need to be stored separately.

     Symmetry (2): rn(j + N/2) = exp(i*(j + N/2)*pi/N) =  exp(i*pi/2) * exp(i*j*pi/N) = I * rn(j),
                   so the acyclic weighting leads naturally to a right-angle transform, i.e. by
                   grouping elements (j, j + N/2) of the original length-N real input vector
                   as complex pairs we can simply apply the weight factor rn(j) to each complex
                   datum that results from the pairing and then go ahead and do a length-(N/2)
                   complex transform on the paired, weighted data, without any additional
                   complex/real wrapper step as is normally required (e.g. for Mersenne-mod.)

     When the transform length is a power of 2 (which, using 64-bit floats, works well for numbers
     up to about the size of F35 or so), this is all that is required. For non-power-of-2 transform
     lengths we combine the above weighting with a Crandall/Fagin-style irrational-base discrete
     weighted transform (IBDWT) to effect the needed Fermat-mod autonegacyclic convolution.

     At the beginning of each loop execution, data are assumed to have already been forward-weighted
     and the initial-radix pass of the S-pass forward FFT done. The loop then does the following:

     (1) Performs the 2nd through (S-1)st passes of the complex DIF forward FFT;
     (2) Does the final-radix forward FFT pass, the complex-pointwise squaring step,
         and the initial-radix inverse FFT pass in a single pass through the data;
     (3) Performs the 2nd through (S-1)st passes of the complex DIT inverse FFT,
         with radices processed in reverse order from the forward FFT (this is
         not necessary for power-of-2 transform lengths, but ensures that DIF radix 1
         equals DIT radix S and vice versa, which is required for steps (2) and (4));
     (4) Does the final-radix inverse FFT pass, the inverse DWT weighting, the carry
         propagation step (with fractional roundoff error check), the forward DWT weighting,
         and the initial-radix forward FFT pass in a single pass through the data.

The scratch array (2nd input argument) is only needed for data table initializations, i.e. if first_entry = TRUE.
*/
int fermat_mod_square(double a[], int arr_scratch[], int n, int ilo, int ihi, uint64 p, int scrnFlag, double *tdiff)
{
	const char func[] = "fermat_mod_square";
	struct qfloat qmul, qwt, qt, qn;	/* qfloats used for DWT weights calculation. */
	struct qfloat qtheta, qr, qi, qc, qs;	/* qfloats used for FFT sincos  calculation. */
	double t1,t2,rt,it,re,im;
	/* Compare qfloat versions of precomputed-table data vs. stdlib double result: */
	const double err_threshold = 1e-12;
	long double theta,mt;
	double adiff, max_adiff = 0.0;	/* Use to store the max abs error between real*8 and real*16 computed values */
	 int64 i1,i2;
	uint64 idiff, max_idiff = 0;

	static int radix_set_save[10] = {1000,0,0,0,0,0,0,0,0,0};
	static int radix0, nchunks; 	/* Stores the first element, RADIX_VEC[0], to workaround an OpemMP loop-control problem */
#if DBG_THREADS
	int num_chunks[MAX_THREADS];		/* Collect stats about how much work done by each thread */
#endif
	static int nradices_prim,nradices_radix0,radix_prim[30];/* RADIX_PRIM stores sequence of complex FFT radices used, in terms of their prime factors.	*/
	static int *index = 0x0, *index_ptmp = 0x0;		/* Bit-reversal index array and array storing S*I mod N values for DWT weights.	*/

	int i,ii,ierr,iter,j,j1,j2,k,l,m,mm,k1,k2;
	static uint64 psave=0;
	static uint32 nsave=0, rad0save=0, new_runlength=0;
	static uint32 nwt,nwt_bits,bw,sw,bits_small;
	static double base[2],baseinv[2],radix_inv;
	/* roots of unity table pairs needed for FFT: */
	static struct complex *rt0 = 0x0, *rt0_ptmp = 0x0, *rt1 = 0x0, *rt1_ptmp = 0x0;
	/* roots of unity table pairs needed for cyclic->acyclic: */
	static struct complex *rn0 = 0x0, *rn0_ptmp = 0x0, *rn1 = 0x0, *rn1_ptmp = 0x0;
	/* IBDWT weights tables: only need if it's a non-power-of-2 FFT length.
	For Fermat-mod transform, wt1 is used to explicitly store inverse weights: */
	static double *wt0 = 0x0, *wt0_ptmp = 0x0,*wt1 = 0x0, *wt1_ptmp = 0x0;		/* DWT weights arrays - unlike mers_mod_square
																				(where it is used as part of a pair-of-small-tables-multiply
																				scheme and we store no inverse weights),
																				here wt1[] is used to store explicit inverse weights	*/
	static int pow2_fft;
	uint32 findex = 0;
	double fracmax,wt,wtinv;
	double max_fp = 0.0, frac_fp, atmp;
	static int first_entry=TRUE;

#ifdef CTIME
	clock_t clock1, clock2;
#else
/* Multithreaded needs wall-clock, not CPU time: */
//	time_t clock1, clock2;
	double clock1, clock2;	// Jun 2014: Switched to getRealTime() code
#endif
#ifdef DBG_TIME	// define at compile time to enable internal timing diagnostics
	const double ICPS = 1.0/CLOCKS_PER_SEC;
	double dt_fwd = 0.0, dt_inv = 0.0, dt_sqr = 0.0, dt_tot = 0.0;
#endif

	/* These came about as a result of multithreading, but now are needed whether built unthreaded or multithreaded */
	static int init_sse2 = FALSE;
	int thr_id;

#ifdef MULTITHREAD

	#ifdef USE_PTHREAD

		int isum;
		static int *thr_ret = 0x0;
		static pthread_t *thread = 0x0;
		static pthread_attr_t attr;
		static struct ferm_thread_data_t *tdat = 0x0;

	  #ifdef USE_THREADPOOL	// Threadpool-based dispatch:

		static int main_work_units = 0, pool_work_units = 0;
		static struct threadpool *tpool = 0x0;
		static int task_is_blocking = TRUE;
		static thread_control_t thread_control = {0,0,0};
		// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
		static task_control_t   task_control = {NULL, (void*)fermat_process_chunk, NULL, 0x0};

	  #else

		/* There are RADIX_VEC[0] independent chunks of work which will be done in parallel by teams of NTHREADS threads at a time
		('per shift', to borrow terminology from the factory floor): */
		int ioffset,nshift,rc;
		void *thr_status;

	  #endif

	#endif

#endif

	radix0 = RADIX_VEC[0];
	nchunks = radix0;
	ASSERT(HERE, TRANSFORM_TYPE == RIGHT_ANGLE, "fermat_mod_square: Incorrect TRANSFORM_TYPE!");

/*...initialize things upon first entry */

	/*...If a new exponent, runlength or radix set, deallocate any already-allocated
	allocatable arrays which are dependent on these values and set first_entry to true:
	*/
	if(n != nsave) new_runlength = TRUE;
	if(p != psave || new_runlength) first_entry=TRUE;

	for(i = 0; i < 10; i++) {
		if(RADIX_VEC[i] != radix_set_save[i]) {
			first_entry=TRUE;
			break;
		}
	}

	if(first_entry)
	{
		first_entry=FALSE;
		psave = p;
		nsave = n;
		N2 =n/2;		/* Complex vector length.	*/
		for(i = 0; i < NRADICES; i++) {
			if(RADIX_VEC[i] == 0) {
				sprintf(cbuf, "%s: RADIX_VEC[i = %d] zero, for i < [NRADICES = %d]!",func,i,NRADICES);
				ASSERT(HERE, 0, cbuf);
			}
			radix_set_save[i] = RADIX_VEC[i];
		}
		for(i = NRADICES; i < 10; i++) {
			if(RADIX_VEC[i] != 0) {
				sprintf(cbuf, "%s: RADIX_VEC[i = %d] nonzero, for i >= [NRADICES = %d]!",func,i,NRADICES);
				ASSERT(HERE, 0, cbuf);
			}
			radix_set_save[i] = 0;
		}

		/*...Check that the binary exponent corresponds to a proper Fermat index: */
		findex = trailz64(p);
		ASSERT(HERE, p >> findex == 1,"fermat_mod_square.c: p >> findex == 1");

		/* Check if the leading radix is one of the ones supported for Fermat-mod: */
		switch(radix0)
		{
		case    7: break;
		case    8: break;
		case   14: break;
		case   15: break;
		case   16: break;
		case   28: break;
		case   30: break;
		case   32: break;
		case   36: break;
		case   56: break;
		case   60: break;
		case   63: break;
		case   64: break;
		case  128: break;
		case  224: break;
		case  240: break;
		case  256: break;
		case  448: break;
		case  480: break;
		case  512: break;
		case  896: break;
		case  960: break;
		case  992: break;	// Include special test radix 31*32 for roundoff tests of F33
		case 1008: break;	// Include special test radix 63*16 for roundoff tests of F33
		case 1024: break;
		case 4032: break;	// Include special test radix 63*64 for roundoff tests of F33
		case 4096: break;	// Include special test radix 2^12 for roundoff tests of F33
		default :
			fprintf(stderr," ERROR : leading radix %d not currently supported for Fermat-mod\n",radix0);
			ierr=ERR_RADIX0_UNAVAILABLE;
			return(ierr);
		}
		/* Final radix must be 16 or 32: */
		// Nov 2017: Add support for emulated final-pass radices 256,1024:
		switch(RADIX_VEC[NRADICES-1])
		{
		case   16: break;
		case   32: break;
	/* Nov 2017: experimental-only:
		case  256: break;
		case 1024: break;
	*/
		default :
			fprintf(stderr," ERROR : final radix %d not currently supported for Fermat-mod\n",RADIX_VEC[NRADICES-1]);
			ierr=ERR_RADIX0_UNAVAILABLE;
			return(ierr);
		}

		/* My array padding scheme requires N/radix0 to be a power of 2, and to be >= 2^DAT_BITS, where the latter
		parameter is set in the Mdata.h file: */

		if(n%radix0 != 0)
		{
			sprintf(cbuf  ,"FATAL: radix0 does not divide N!\n");
			fprintf(stderr,"%s", cbuf);
			ASSERT(HERE, 0,cbuf);
		}

		/* Make sure n/radix0 is a power of 2: */
		i = n/radix0;
		if((i >> trailz32(i)) != 1)
		{
			sprintf(cbuf  ,"FATAL: n/radix0 not a power of 2!\n");
			fprintf(stderr,"%s", cbuf);
			ASSERT(HERE, 0,cbuf);
		}

		if(DAT_BITS < 31)
		{
			/* Now make sure n/radix0 is sufficiently large (unless n < 2^DAT_BITS, in which case it doesn't matter): */
			if(i < (1 << DAT_BITS))
			{
			//	sprintf(cbuf  ,"FATAL: n/radix0 must be >= %u!\n", (1 << DAT_BITS));	fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				// Mar 2018: Switch to 'soft' assertion error here, e.g. for timing tests at small FFT lengths:
				sprintf(cbuf  ,"n/radix0 must be >= %u! Skipping this radix combo.\n", (1 << DAT_BITS));	WARN(HERE, cbuf, "", 1); return(ERR_ASSERT);
			}

			/* We also have a lower limit on 2^DAT_BITS set by the wrapper_square routine: */
			// Nov 2017: This check applies only for non-synthesized final-pass radices, not the large synthesized ones, 256,512,1024.
			if((RADIX_VEC[NRADICES-1] <= 32) && (1 << DAT_BITS) < 2*RADIX_VEC[NRADICES-1])
			{
				sprintf(cbuf  ,"FATAL: Value of DAT_BITS means final FFT radix may not exceed = %u!\n", (1 << (DAT_BITS-1)));
				fprintf(stderr,"%s", cbuf);
				ASSERT(HERE, 0,cbuf);
			}
		}

		sprintf(cbuf,"Using complex FFT radices*");
		char_addr = strstr(cbuf,"*");
		for(i = 0; i < NRADICES; i++)
		{
			sprintf(char_addr,"%10d",RADIX_VEC[i]); char_addr += 10;
		}; sprintf(char_addr,"\n");

		if(INTERACT)
		{
			fprintf(stderr,"%s",cbuf);
		}
		else
		{
			fp = mlucas_fopen(STATFILE,"a");	fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;
			if (scrnFlag)	/* Echo output to stddev */
			{
				fprintf(stderr,"%s",cbuf);
			}
		}

	/*...******Forward FFT****** permuted sincos index array is here: first, calculate the needed dimension...	*/
		k =0;
		mm=radix0;			/* First radix requires no twiddle factors.	*/
		/* We do the final DIF FFT radix within the dyadic_square routine, so store
		that block of sincos data there, where they can be merged with the wrapper sincos data:
		*/
		for(i=1; i<NRADICES-1; i++)
		{
			k  += mm;
			mm *= RADIX_VEC[i];
		//	printf("\tR = %u, mm = %u\n",RADIX_VEC[i],mm);
		}
		if(RADIX_VEC[i] >= 256) k += mm;	// For emulated final-pass radices 256,512,1024, need to include loop-exit value of mm in the alloc

		if(mm*RADIX_VEC[NRADICES-1] != N2) {
			sprintf(cbuf  ,"FATAL: product of radices not equal to complex vector length\n");
			fprintf(stderr,"%s", cbuf);
			ASSERT(HERE, 0,cbuf);
		}

		index_ptmp = ALLOC_INT(index_ptmp, k);
		if(!index_ptmp)
		{
			sprintf(cbuf  ,"FATAL: unable to allocate array INDEX in %s.\n",func);
			fprintf(stderr,"%s", cbuf);
			ASSERT(HERE, 0,cbuf);
		}
		index = ALIGN_INT(index_ptmp);
	//	printf("Alloc index[%u]...\n",k);
		/*...Forward (DIF) FFT sincos data are in bit-reversed order. We define a separate last-pass twiddles
		array within the routine wrapper_square, since that allows us to merge those nicely with the wrapper sincos data.	*/
		k = 0; l =0; mm = 1;
		/*...First radix needs no twiddle factors, just need it for building the radix_prim array.	*/
		switch(radix0)
		{
		/*
		case 2 :
			nradices_radix0 = 1;
			radix_prim[l++] = 2; break;
		case 3 :
			nradices_radix0 = 1;
			radix_prim[l++] = 3; break;
		case 4 :
			nradices_radix0 = 2;
			radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		*/
		case 5 :
			nradices_radix0 = 1;
			radix_prim[l++] = 5; break;
		case 6 :
			nradices_radix0 = 2;
			radix_prim[l++] = 3; radix_prim[l++] = 2; break;
		case 7 :
			nradices_radix0 = 1;
			radix_prim[l++] = 7; break;
		case 8 :
			nradices_radix0 = 3;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 9 :
			nradices_radix0 = 2;
			radix_prim[l++] = 3; radix_prim[l++] = 3; break;
		case 10 :
			nradices_radix0 = 2;
			radix_prim[l++] = 5; radix_prim[l++] = 2; break;
		case 11 :
			nradices_radix0 = 1;
			radix_prim[l++] = 11; break;
		case 12 :
			nradices_radix0 = 3;
			radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 13 :
			nradices_radix0 = 1;
			radix_prim[l++] = 13; break;
		case 14 :
			nradices_radix0 = 2;
			radix_prim[l++] = 7; radix_prim[l++] = 2; break;
		case 15 :
			nradices_radix0 = 2;
			radix_prim[l++] = 5; radix_prim[l++] = 3; break;
		case 16 :
			nradices_radix0 = 4;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 18:
			nradices_radix0 = 3;
			radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; break;
		case 20:
			nradices_radix0 = 3;
			radix_prim[l++] = 5; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 22 :
			nradices_radix0 = 2;
			radix_prim[l++] =11; radix_prim[l++] = 2; break;
		case 24 :
			nradices_radix0 = 4;
			radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		/*
		case 25 :
			nradices_radix0 = 2;
			radix_prim[l++] = 5; radix_prim[l++] = 5; break;
		*/
		case 26 :
			nradices_radix0 = 2;
			radix_prim[l++] =13; radix_prim[l++] = 2; break;
		case 28:
			nradices_radix0 = 3;
			radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 30:
			nradices_radix0 = 3;
			radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; break;
		case 32 :
			nradices_radix0 = 5;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 36 :
			nradices_radix0 = 4;
			radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 56 :
			nradices_radix0 = 4;
			radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 60 :
			nradices_radix0 = 4;
			radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 63 :
			nradices_radix0 = 3;
			radix_prim[l++] = 7; radix_prim[l++] = 3; radix_prim[l++] = 3; break;
		case 64 :
			nradices_radix0 = 6;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		/*
		case 72 :
			nradices_radix0 = 5;
			radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 112:
			nradices_radix0 = 5;
			radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 120:
			nradices_radix0 = 5;
			radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		*/
		case 128 :
			nradices_radix0 = 7;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 224:
			nradices_radix0 = 6;
			radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 240:
			nradices_radix0 = 6;
			radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 256 :
			nradices_radix0 = 8;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 512 :
			nradices_radix0 = 9;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 960:
			nradices_radix0 = 8;
			radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 992:
			nradices_radix0 = 6;
			radix_prim[l++] =31; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 1008:
			nradices_radix0 = 7;
			radix_prim[l++] = 7; radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 1024:
			nradices_radix0 = 10;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 4032:
			nradices_radix0 = 9;
			radix_prim[l++] = 7; radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 4096:
			nradices_radix0 = 12;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		default :
			sprintf(cbuf  ,"FATAL: radix %d not available for Fermat-mod transform. Halting...\n",RADIX_VEC[i]);
			fprintf(stderr,"%s", cbuf);
			ASSERT(HERE, 0,cbuf);
		}
		nradices_prim = l;	for( ; l < 30; l++) { radix_prim[l] = 0; }	// Zero any higher elements which may have been previously set due
								// to use of a smoother n. (Not necessary since use nradices_prim to control array access, but nice to do.

		int iflag = (RADIX_VEC[i] >= 256);	// For emulated final-pass radices 256,1024, need to include loop-exit value of mm in the alloc
		for(i = 1; i < NRADICES-1+iflag; i++)
		{
			/*...Allocate and initialize an index array containing MM indices...	*/
			if(i < (NRADICES-1+iflag))
			{
				mm *= RADIX_VEC[i-1];	/* MM = product of all the preceding radices	*/
			//	printf("Init index[%u-%u]...\n",k,k+mm-1);
				for(m = 0; m < mm; m++)
				{
					index[k+m] = m;
				}
				/*...then bit-reverse INDEX with respect to the accumulated radices.
				The order of radices sent to bit_reverse_int is the reverse of that in which these radices are processed
				in the forward (decimation in frequency) FFT. This is moot for a power-of-2 FFT (or any FFT whose length
				is a prime power), but necessary for general vector lengths which are a product of 2 or more distinct primes.

				If the current (Ith) radix is composite with distinct prime factors (e.g. 15 = 3*5), we must specify these
				factors here in the opposite order from that which is used in the actual FFT-pass routine. For example,
				if the radix-15 pass implementation does 5 radix-3 DFTs, followed by 3 radix-5 DFTs, then we send (3,5)
				as the corresponding reverse-ordered prime radices to the bit-reversal routine, not (5,3).
				*/
				bit_reverse_int(&index[k],mm,l,&radix_prim[l-1],-1,(int *)arr_scratch);
			//	printf("index[%u-%u] = ",k,k+mm-1); { for(m = 0; m < mm; m++) { printf("%u.",index[k+m]); } printf("\n"); }
				k += mm;
			}
			/*...All radices beyond the initial-pass one are assumed to be powers of 2 in [8,32]:	*/
		  if(i < (NRADICES-1)) {	// Skip any emulated final-pass radices for this check:
			switch(RADIX_VEC[i])
			{
			case 8 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			case 16 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			case 32 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			default :
				sprintf(cbuf  ,"FATAL: intermediate radix %d not available. Halting...\n",RADIX_VEC[i]);
				fprintf(stderr,"%s", cbuf);
				ASSERT(HERE, 0,cbuf);
			}
		  }
		}
		// Final-pass radix must be one of the following values, where 256,1024 are the synthesized ones added in Nov 2017:
		switch(RADIX_VEC[i-iflag])
		{
		case 16 :
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 32 :
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 256:
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 1024:
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		default :
			sprintf(cbuf  ,"FATAL: final-pass radix %d not available. Halting...\n",RADIX_VEC[i]);
			fprintf(stderr,"%s", cbuf);
			ASSERT(HERE, 0,cbuf);
		}
		nradices_prim = l;	for( ; l < 30; l++) { radix_prim[l] = 0; }	// Zero any higher elements which may have been previously set due
								// to use of a smoother n. (Not necessary since use nradices_prim to control array access, but nice to do.
		/* Only need an IBDWT weights table if it's a non-power-of-2-length
		Fermat-mod transform, in which case the table has {odd part of N}
		distinct elements.
		*/
		nwt_bits = 0;	/* Always set = 0 for Fermat-mod */
		/* Vector length a power of 2? */
		nwt = (n >> trailz32(n));
		if(nwt == 1)
		{
			pow2_fft = TRUE;
		}
		else
		{
			pow2_fft = FALSE;
		}

		bw = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw = n - bw;	/* Number of smallwords.	*/

		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)radix0));

		bits_small = p/n;			/* number of bits in a smallword.	*/
		if(pow2_fft)
		{
			/* If power-of-2 runlength, no IBDWT gets done, make bases the same: */
			base   [0] = (double)(1 << bits_small);	base   [1] = base[0]	;
			baseinv[0] = 1.0/base[0];				baseinv[1] = baseinv[1]	;	/* don't need extended precision for this since both bases are powers of 2.	*/
		}
		else
		{
			base   [0] = (double)(1 << bits_small);	base   [1] = (double)(2*base[0]);
			baseinv[0] = (double)(1.0/base[0]    );	baseinv[1] = (double)(1.0/base[1]);	/* don't need extended precision for this since both bases are powers of 2.	*/

			/*...stuff for the reduced-length DWT weights arrays is here:	*/
			wt0_ptmp = ALLOC_DOUBLE(wt0_ptmp, nwt);	if(!wt0_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array WT0 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }; wt0 = ALIGN_DOUBLE(wt0_ptmp);
			wt1_ptmp = ALLOC_DOUBLE(wt1_ptmp, nwt);	if(!wt1_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array WT1 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }; wt1 = ALIGN_DOUBLE(wt1_ptmp);

		/******************************************************************/
		/* Crandall/Fagin weighting factors and number of bits per digit. */
		/******************************************************************/

			/* Double-check that sw*nwt (where nwt is the odd factor of N) is divisible by N: */
		//	printf("sw,nwt,n = %u,%u,%u; sw*nwt mod n = %u\n",sw,nwt,n, (uint64)sw*nwt % n);
			ASSERT(HERE, (uint64)sw*nwt % n == 0,"fermat_mod_square.c: sw*nwt % n == 0");
			SW_DIV_N = sw*nwt/n;

			qn   = i64_to_q((int64) nwt);
			qt   = qfinv(qn);			/* 1/nwt...	 */
			qt   = qfmul(qt, QLN2);		/* ...(1/nwt)*ln(2)...	*/
			qmul = qfexp(qt);			/* ...and get 2^(1/nwt) via exp[(1/nwt)*ln(2)]. */
			qwt  = QONE;				/* init weights multiplier chain. */

			t1 = qfdbl(qmul);
			/* Compare qfloat versions of precomputed-table data vs. stdlib double result: */
			t2 = pow(2.0, 1.0/nwt);
			adiff = ABS(t1-t2);
			if(adiff > max_adiff)
				max_adiff = adiff;
			if(adiff > err_threshold)
			{
				i1 = *(uint64 *)&t1;
				i2 = *(uint64 *)&t2;
				idiff = ABS(i1-i2);
				if(idiff > max_idiff)
					max_idiff = idiff;
				sprintf(cbuf,"INFO: QWT = %20.15f, DWT = %20.15f DIFFER BY %20.0f\n", t1, t2, (double)idiff);
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}

			for(i=0; i<nwt; i++)
			{
				/* Ith DWT weight factor = 2^[i/nwt], where the exponent is done using floating divide.	*/
				wt0[i] = qfdbl(qwt);
				t1 = wt0[i];
				t2 = pow(2.0, 1.0*i/nwt);
				adiff = ABS(t1-t2);
				if(adiff > max_adiff)
					max_adiff = adiff;
				if(adiff > err_threshold)
				{
					i1 = *(uint64 *)&t1;
					i2 = *(uint64 *)&t2;
					idiff = ABS(i1-i2);
					if(idiff > max_idiff)
						max_idiff = idiff;
					sprintf(cbuf,"INFO: I = %8d: QWT0 = %20.15f, DWT0 = %20.15f DIFFER BY %20.0f\n", i, t1, t2, (double)idiff);
					fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				}

				/* Inverse DWT weight factor:	*/
				wt1[i] = qfdbl(qfinv(qwt));
				t1 = wt1[i];
				t2 = 1.0/wt0[i];
				i1 = *(uint64 *)&t1;
				i2 = *(uint64 *)&t2;
				adiff = ABS(t1-t2);
				if(adiff > max_adiff)
					max_adiff = adiff;
				if(adiff > err_threshold)
				{
					i1 = *(uint64 *)&t1;
					i2 = *(uint64 *)&t2;
					idiff = ABS(i1-i2);
					if(idiff > max_idiff)
						max_idiff = idiff;
					sprintf(cbuf,"INFO: I = %8d: QWT1 = %20.15f, DWT1 = %20.15f DIFFER BY %20.0f\n", i, t1, t2, (double)idiff);
					fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				}
				qwt= qfmul(qwt, qmul);
			}
		}

		/**********************************************/
		/* Roots of unity table pairs needed for FFT: */
		/**********************************************/

		/* No need for a fancy NINT here: */
		NRT_BITS = (uint32)(log(sqrt(1.0*n))/log(2.0) + 0.5);
		NRT = 1 << NRT_BITS;
		if(n%NRT){ sprintf(cbuf,"FATAL: NRT does not divide N!\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		NRTM1 = NRT - 1;

		/*...The rt0 array stores the (0:NRT-1)th powers of the [N2]th root of unity
		(i.e. will be accessed using the lower (NRT) bits of the integer sincos index):
		*/
		rt0_ptmp = ALLOC_COMPLEX(rt0_ptmp, NRT);
		if(!rt0_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array RT0 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		rt0 = ALIGN_COMPLEX(rt0_ptmp);

		qt     = i64_to_q((int64)N2);
		qtheta = qfdiv(Q2PI, qt);	/* 2*pi/(N/2) */
		qr     = qfcos(qtheta);
		qi     = qfsin(qtheta);
		qc     = QONE; qs = QZRO;	/* init sincos multiplier chain. */

		t1 = qfdbl(qr);
		/* Compare qfloat versions of precomputed-table data vs. stdlib double result: */
		theta = qfdbl(Q2PI)/N2;
		t2 = cos(theta);
		adiff = ABS(t1-t2);
		if(adiff > max_adiff)
			max_adiff = adiff;
		if(adiff > err_threshold)
		{
			i1 = *(uint64 *)&t1;
			i2 = *(uint64 *)&t2;
			idiff = ABS(i1-i2);
			if(idiff > max_idiff)
				max_idiff = idiff;
			sprintf(cbuf,"INFO: QCOS1= %20.15f, DCOS = %20.15f DIFFER BY %20.0f\n", t1, t2, (double)idiff);
			fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}

		t1 = qfdbl(qi);
		/* Compare qfloat versions of precomputed-table data vs. stdlib double result: */
		t2 = sin(theta);
		adiff = ABS(t1-t2);
		if(adiff > max_adiff)
			max_adiff = adiff;
		if(adiff > err_threshold)
		{
			i1 = *(uint64 *)&t1;
			i2 = *(uint64 *)&t2;
			idiff = ABS(i1-i2);
			if(idiff > max_idiff)
				max_idiff = idiff;
			sprintf(cbuf,"INFO: QSIN1= %20.15f, DSIN = %20.15f DIFFER BY %20.0f\n", t1, t2, (double)idiff);
			fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}

		qt = QZRO;

		for(i=0; i<NRT; i++)
		{
			qc = qfcos(qt);
			t1 = qfdbl(qc);
			/* Compare qfloat versions of precomputed-table data vs. stdlib double result: */
			mt = i*theta;
			t2 = cos(mt);
			adiff = ABS(t1-t2);
			if(adiff > max_adiff)
				max_adiff = adiff;
			if(adiff > err_threshold)
			{
				i1 = *(uint64 *)&t1;
				i2 = *(uint64 *)&t2;
				idiff = ABS(i1-i2);
				if(idiff > max_idiff)
					max_idiff = idiff;
				sprintf(cbuf,"INFO: I = %8d: QCOS = %20.15f, DCOS = %20.15f DIFFER BY %20.0f\n", i, t1, t2, (double)idiff);
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}
			rt0[i].re = t1;

			qs = qfsin(qt);
			t1 = qfdbl(qs);
			/* Compare qfloat versions of precomputed-table data vs. stdlib double result: */
			t2 = sin(mt);
			adiff = ABS(t1-t2);
			if(adiff > max_adiff)
				max_adiff = adiff;
			if(adiff > err_threshold)
			{
				i1 = *(uint64 *)&t1;
				i2 = *(uint64 *)&t2;
				idiff = ABS(i1-i2);
				if(idiff > max_idiff)
					max_idiff = idiff;
				sprintf(cbuf,"INFO: I = %8d: QSIN = %20.15f, DSIN = %20.15f DIFFER BY %20.0f\n", i, t1, t2, (double)idiff);
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}
			rt0[i].im = t1;

			qt = qfadd(qt, qtheta);

			/* Get next terms of the recurrence:  qcnew = qcold*qr - qsold*qi,  qsnew = qcold*qi + qsold*qr:
			EWM - this needs further debug!
			qn = qfmul(qc, qr); qt = qfmul(qs, qi); qmul = qfsub(qn, qt);	// Store qcnew in qmul for now.
			qn = qfmul(qc, qi); qt = qfmul(qs, qr); qs   = qfadd(qn, qt); qc = qmul;
			*/
		}

		/*...The rt1 array stores the (0:(n/2)/NRT-1)th powers of the [(n/2)/NRT]th root of unity
		(and will be accessed using the upper bits, <NRT:31>, of the integer sincos index):
		*/
		rt1_ptmp = ALLOC_COMPLEX(rt1_ptmp, n/(2*NRT));
		if(!rt1_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array RT1 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		rt1 = ALIGN_COMPLEX(rt1_ptmp);

		qn     = i64_to_q((int64)NRT);
		qt     = i64_to_q((int64)N2);
		qt     = qfdiv(qn, qt);		/*      NRT/(N/2) */
		qtheta = qfmul(Q2PI, qt);	/* 2*pi*NRT/(N/2) */
		qr     = qfcos(qtheta);
		qi     = qfsin(qtheta);
		qc  = QONE; qs = QZRO;	/* init sincos multiplier chain. */

		t1 = qfdbl(qr);
		/* Compare qfloat versions of precomputed-table data vs. stdlib double result: */
		theta = qfdbl(Q2PI)*NRT/N2;
		t2 = cos(theta);
		adiff = ABS(t1-t2);
		if(adiff > max_adiff)
			max_adiff = adiff;
		if(adiff > err_threshold)
		{
			i1 = *(uint64 *)&t1;
			i2 = *(uint64 *)&t2;
			idiff = ABS(i1-i2);
			if(idiff > max_idiff)
				max_idiff = idiff;
			sprintf(cbuf,"INFO: QCOS2= %20.15f, DCOS = %20.15f DIFFER BY %20.0f\n", t1, t2, (double)idiff);
			fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}

		t1 = qfdbl(qi);
		/* Compare qfloat versions of precomputed-table data vs. stdlib double result: */
		t2 = sin(theta);
		adiff = ABS(t1-t2);
		if(adiff > max_adiff)
			max_adiff = adiff;
		if(adiff > err_threshold)
		{
			i1 = *(uint64 *)&t1;
			i2 = *(uint64 *)&t2;
			idiff = ABS(i1-i2);
			if(idiff > max_idiff)
				max_idiff = idiff;
			sprintf(cbuf,"INFO: QSIN2= %20.15f, DSIN = %20.15f DIFFER BY %20.0f\n", t1, t2, (double)idiff);
			fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}

		qt = QZRO;

		for(i=0; i<(N2/NRT); i++)
		{
			qc = qfcos(qt);
			t1 = qfdbl(qc);
			/* Compare qfloat versions of precomputed-table data vs. stdlib double result: */
			mt = i*theta;
			t2 = cos(mt);
			adiff = ABS(t1-t2);
			if(adiff > max_adiff)
				max_adiff = adiff;
			if(adiff > err_threshold)
			{
				i1 = *(uint64 *)&t1;
				i2 = *(uint64 *)&t2;
				idiff = ABS(i1-i2);
				if(idiff > max_idiff)
					max_idiff = idiff;
				sprintf(cbuf,"INFO: J = %8d: QCOS = %20.15f, DCOS = %20.15f DIFFER BY %20.0f\n", i, t1, t2, (double)idiff);
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}
			rt1[i].re = t1;

			qs = qfsin(qt);
			t1 = qfdbl(qs);
			/* Compare qfloat versions of precomputed-table data vs. stdlib double result: */
			t2 = sin(mt);
			adiff = ABS(t1-t2);
			if(adiff > max_adiff)
				max_adiff = adiff;
			if(adiff > err_threshold)
			{
				i1 = *(uint64 *)&t1;
				i2 = *(uint64 *)&t2;
				idiff = ABS(i1-i2);
				if(idiff > max_idiff)
					max_idiff = idiff;
				sprintf(cbuf,"INFO: J = %8d: QSIN = %20.15f, DSIN = %20.15f DIFFER BY %20.0f\n", i, t1, t2, (double)idiff);
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}
			rt1[i].im = t1;

			qt = qfadd(qt, qtheta);

			/* Get next terms of the recurrence:  qcnew = qcold*qr - qsold*qi,  qsnew = qcold*qi + qsold*qr. */
			/*qn = qfmul(qc, qr); qt = qfmul(qs, qi); qmul = qfsub(qn, qt);	* Store qcnew in qmul for now. */
			/*qn = qfmul(qc, qi); qt = qfmul(qs, qr); qs   = qfadd(qn, qt); qc = qmul;	*/
		}

		/****************************************************************/
		/* roots of unity table pairs needed for cyclic->acyclic        */
		/* are identical to roots of unity table pairs needed for FFT,  */
		/* except that we deal with Nth roots of -1 (i.e. (2*N)th roots */
		/* of unity) rather than (N/2)th roots of unity:                */
		/****************************************************************/

		/*...The rn0 array stores the (0:NRT-1)th powers of the [2*n]th root of unity
		(i.e. will be accessed using the lower (NRT) bits of the integer sincos index):
		*/
		rn0_ptmp = ALLOC_COMPLEX(rn0_ptmp, NRT);	if(!rn0_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array RN0 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); } rn0 = ALIGN_COMPLEX(rn0_ptmp);

		qt     = i64_to_q((int64)N2);
		qtheta = qfdiv(QPIHALF, qt);	/* (2*pi)/(2*N) = (pi/2)/(N/2) */
		qr     = qfcos(qtheta);
		qi     = qfsin(qtheta);
		qc     = QONE; qs = QZRO;	/* init sincos multiplier chain. */

		t1 = qfdbl(qr);
		theta = qfdbl(QPIHALF)/N2;
		t2 = cos(theta);
		adiff = ABS(t1-t2);
		if(adiff > max_adiff)
			max_adiff = adiff;
		if(adiff > err_threshold)
		{
			i1 = *(uint64 *)&t1;
			i2 = *(uint64 *)&t2;
			idiff = ABS(i1-i2);
			if(idiff > max_idiff)
				max_idiff = idiff;
			sprintf(cbuf,"INFO: QCOS1= %20.15f, DCOS = %20.15f DIFFER BY %20.0f\n", t1, t2, (double)idiff);
			fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}

		t1 = qfdbl(qi);
		t2 = sin(theta);
		adiff = ABS(t1-t2);
		if(adiff > max_adiff)
			max_adiff = adiff;
		if(adiff > err_threshold)
		{
			i1 = *(uint64 *)&t1;
			i2 = *(uint64 *)&t2;
			idiff = ABS(i1-i2);
			if(idiff > max_idiff)
				max_idiff = idiff;
			sprintf(cbuf,"INFO: QSIN1= %20.15f, DSIN = %20.15f DIFFER BY %20.0f\n", t1, t2, (double)idiff);
			fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}

		qt = QZRO;

		for(i=0; i<NRT; i++)
		{
			qc = qfcos(qt);
			t1 = qfdbl(qc);
			mt = i*theta;
			t2 = cos(mt);
			adiff = ABS(t1-t2);
			if(adiff > max_adiff)
				max_adiff = adiff;
			if(adiff > err_threshold)
			{
				i1 = *(uint64 *)&t1;
				i2 = *(uint64 *)&t2;
				idiff = ABS(i1-i2);
				if(idiff > max_idiff)
					max_idiff = idiff;
				sprintf(cbuf,"INFO: I = %8d: QCOS = %20.15f, DCOS = %20.15f DIFFER BY %20.0f\n", i, t1, t2, (double)idiff);
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}
			rn0[i].re = t1;

			qs = qfsin(qt);
			t1 = qfdbl(qs);
			t2 = sin(mt);
			adiff = ABS(t1-t2);
			if(adiff > max_adiff)
				max_adiff = adiff;
			if(adiff > err_threshold)
			{
				i1 = *(uint64 *)&t1;
				i2 = *(uint64 *)&t2;
				idiff = ABS(i1-i2);
				if(idiff > max_idiff)
					max_idiff = idiff;
				sprintf(cbuf,"INFO: I = %8d: QSIN = %20.15f, DSIN = %20.15f DIFFER BY %20.0f\n", i, t1, t2, (double)idiff);
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}
			rn0[i].im = t1;

			/*fprintf(stderr,"I = %d; RT0 = %20.10f %20.10f\n",i,rn0[i].re,rn0[i].im);	*/
			/*		errprint_sincos(&rn0[i].re,&rn0[i].im,(double)(mt));	* Workaround for the DEC Unix 4.0 real*16 sincos bug	*/

			qt = qfadd(qt, qtheta);
			/* Get next terms of the recurrence:  qcnew = qcold*qr - qsold*qi,  qsnew = qcold*qi + qsold*qr. */
			/*qn = qfmul(qc, qr); qt = qfmul(qs, qi); qmul = qfsub(qn, qt);	* Store qcnew in qmul for now. */
			/*qn = qfmul(qc, qi); qt = qfmul(qs, qr); qs   = qfadd(qn, qt); qc = qmul;	*/
		}

		/*...The rn1 array stores the (0:(n/2)/NRT-1)th powers of the [(n/2)/NRT]th root of unity
		(and will be accessed using the upper bits, <NRT:31>, of the integer sincos index):
		*/
		rn1_ptmp = ALLOC_COMPLEX(rn1_ptmp, N2/NRT);	if(!rn1_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array RN1 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); } rn1 = ALIGN_COMPLEX(rn1_ptmp);

		qn     = i64_to_q((int64)NRT);
		qt     = i64_to_q((int64)N2);
		qt     = qfdiv(qn, qt);			/*      NWT/(N/2) */
		qtheta = qfmul(QPIHALF, qt);	/* 2*pi*NWT/(2*N) = (pi/2)*NWT/(N/2)) */
		qr     = qfcos(qtheta);
		qi     = qfsin(qtheta);
		qc     = QONE; qs = QZRO;		/* init sincos multiplier chain. */

		t1 = qfdbl(qr);
		theta = qfdbl(QPIHALF)*NRT/N2;
		t2 = cos(theta);
		adiff = ABS(t1-t2);
		if(adiff > max_adiff)
			max_adiff = adiff;
		if(adiff > err_threshold)
		{
			i1 = *(uint64 *)&t1;
			i2 = *(uint64 *)&t2;
			idiff = ABS(i1-i2);
			if(idiff > max_idiff)
				max_idiff = idiff;
			sprintf(cbuf,"INFO: QCOS2= %20.15f, DCOS = %20.15f DIFFER BY %20.0f\n", t1, t2, (double)idiff);
			fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}

		t1 = qfdbl(qi);
		t2 = sin(theta);
		adiff = ABS(t1-t2);
		if(adiff > max_adiff)
			max_adiff = adiff;
		if(adiff > err_threshold)
		{
			i1 = *(uint64 *)&t1;
			i2 = *(uint64 *)&t2;
			idiff = ABS(i1-i2);
			if(idiff > max_idiff)
				max_idiff = idiff;
			sprintf(cbuf,"INFO: QSIN2= %20.15f, DSIN = %20.15f DIFFER BY %20.0f\n", t1, t2, (double)idiff);
			fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}

		qt = QZRO;

		for(i=0; i<(N2/NRT); i++)
		{
			qc = qfcos(qt);
			t1 = qfdbl(qc);
			mt = i*theta;
			t2 = cos(mt);
			adiff = ABS(t1-t2);
			if(adiff > max_adiff)
				max_adiff = adiff;
			if(adiff > err_threshold)
			{
				i1 = *(uint64 *)&t1;
				i2 = *(uint64 *)&t2;
				idiff = ABS(i1-i2);
				if(idiff > max_idiff)
					max_idiff = idiff;
				sprintf(cbuf,"INFO: J = %8d: QCOS = %20.15f, DCOS = %20.15f DIFFER BY %20.0f\n", i, t1, t2, (double)idiff);
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}
			rn1[i].re = t1;

			qs = qfsin(qt);
			t1 = qfdbl(qs);
			t2 = sin(mt);
			adiff = ABS(t1-t2);
			if(adiff > max_adiff)
				max_adiff = adiff;
			if(adiff > err_threshold)
			{
				i1 = *(uint64 *)&t1;
				i2 = *(uint64 *)&t2;
				idiff = ABS(i1-i2);
				if(idiff > max_idiff)
					max_idiff = idiff;
				sprintf(cbuf,"INFO: J = %8d: QSIN = %20.15f, DSIN = %20.15f DIFFER BY %20.0f\n", i, t1, t2, (double)idiff);
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}
			rn1[i].im = t1;

			/*fprintf(stderr,"I = %d; RT1 = %20.10f %20.10f\n",i,rn1[i].re,rn1[i].im);	*/
			/*		errprint_sincos(&rn1[i].re,&rn1[i].im,(double)(mt));	* Workaround for the DEC Unix 4.0 real*16 sincos bug	*/

			qt = qfadd(qt, qtheta);
			/* Get next terms of the recurrence:  qcnew = qcold*qr - qsold*qi,  qsnew = qcold*qi + qsold*qr. */
			/*qn = qfmul(qc, qr); qt = qfmul(qs, qi); qmul = qfsub(qn, qt);	* Store qcnew in qmul for now. */
			/*qn = qfmul(qc, qi); qt = qfmul(qs, qr); qs   = qfadd(qn, qt); qc = qmul;	*/
		}

		if(max_adiff > err_threshold)
		{
			fprintf(stderr, "%s:\n",func);
			fprintf(stderr, " Max abs error between real*8 and real*16 computed values = %20.15f\n",         max_adiff);
			fprintf(stderr, " Max bit error between real*8 and real*16 computed values = %20.0f \n", (double)max_idiff);

			ASSERT(HERE, (max_adiff < 100*err_threshold),"Max error between real*8 and real*16 unacceptably high - quitting.");
		}

	#ifdef MULTITHREAD

	  #ifdef USE_PTHREAD

		free((void *)thr_ret); thr_ret = 0x0;
		free((void *)thread ); thread  = 0x0;
		free((void *)tdat   ); tdat    = 0x0;

		thr_ret = (int *)calloc(radix0, sizeof(int));
		thread  = (pthread_t *)calloc(radix0, sizeof(pthread_t));
		tdat    = (struct ferm_thread_data_t *)calloc(radix0, sizeof(struct ferm_thread_data_t));

		/* Initialize and set thread detached attribute */
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		/* Populate the thread-specific data structs: */
		for(i = 0; i < radix0; ++i)
		{
			tdat[i].tid = i;
			tdat[i].retval = &thr_ret[i];
		//	tdat[i].arrdat = a;			/* Main data array */	Jan 2018: To support PRP-cofactor testing, must allow user to feed different input vectors to this modpow routine!
			tdat[i].arr_scratch = arr_scratch;
			tdat[i].n = n;					/* Chunksize */
			tdat[i].rt0 = rt0;	/* Roots table 1 */
			tdat[i].rt1 = rt1;	/* Roots table 2 */
			tdat[i].index = index;				/* Bit-reversal index array */
			tdat[i].nradices_prim = nradices_prim;
			tdat[i].radix_prim = radix_prim;
		#ifdef DBG_TIME
			tdat[i].dt_fwd = tdat[i].dt_inv = tdat[i].dt_sqr = tdat[i].dt_tot = 0.0;
		#endif
		}
	  #endif

	  #ifdef USE_THREADPOOL	// Threadpool-based dispatch:

		ASSERT(HERE, MAX_THREADS == get_num_cores(), "MAX_THREADS not set or incorrectly set!");

		if(radix0 % NTHREADS != 0) fprintf(stderr,"%s: radix0 not exactly divisible by NTHREADS - This will hurt performance.\n",func);

		// MacOS does weird things with threading (e.g. Idle" main thread burning 100% of 1 CPU)
		// so on that platform try to be clever and interleave main-thread and threadpool-work processing
		#if 0//def OS_TYPE_MACOSX

			if(NTHREADS > 1) {
				main_work_units = radix0/NTHREADS;
				pool_work_units = radix0 - main_work_units;
				ASSERT(HERE, 0x0 != (tpool = threadpool_init(NTHREADS-1, MAX_THREADS, pool_work_units, &thread_control)), "threadpool_init failed!");
				printf("%s: Init threadpool of %d threads\n",func,NTHREADS-1);
			} else {
				printf("%s: NTHREADS = 1: Using main execution thread, no threadpool needed.\n",func);
			}

		#else

			main_work_units = 0;
			pool_work_units = radix0;
			ASSERT(HERE, 0x0 != (tpool = threadpool_init(NTHREADS, MAX_THREADS, pool_work_units, &thread_control)), "threadpool_init failed!");
			printf("%s: Init threadpool of %d threads\n",func,NTHREADS);

		#endif

	  #endif

	#endif
	}

	/* 	This set of init-mode calls needs to go below above init-block because several
	of the inits need the primitive-radix data to have been inited.
	*/
	// New run, or need to up #threads in local-store inits:
	if(init_sse2 == FALSE || (init_sse2 < nchunks)) {
	//	init_sse2 = TRUE;
		init_sse2 = nchunks;	// Use *value* of init_sse2 to store #threads
		thr_id = -1;
		radix8_dif_pass (0x0, 0, 0x0, 0x0, 0x0, 0, 0, init_sse2, thr_id);
		radix16_dif_pass(0x0, 0, 0x0, 0x0, 0, 0, 0x0, init_sse2, thr_id);
		radix32_dif_pass(0x0, 0, 0x0, 0x0, 0, 0, 0x0, init_sse2, thr_id);

		/* The dyadic-square routines need a few more params passed in init-mode than do the standalone FFT-pass routines: */
		// Sep 2014: Add roots-table pointers to arglist since now need these to init modified 2-table-CMUL scheme:
		radix16_dyadic_square(0x0, arr_scratch, n, radix0, rt0, rt1, 0, nradices_prim, radix_prim, 0, init_sse2,thr_id);
		radix32_dyadic_square(0x0, arr_scratch, n, radix0, rt0, rt1, 0, nradices_prim, radix_prim, 0, init_sse2,thr_id);

		radix8_dit_pass (0x0, 0, 0x0, 0x0, 0, 0, 0x0, init_sse2, thr_id);
		radix16_dit_pass(0x0, 0, 0x0, 0x0, 0, 0, 0x0, init_sse2, thr_id);
		radix32_dit_pass(0x0, 0, 0x0, 0x0, 0, 0, 0x0, init_sse2, thr_id);
	}
	// New runlength or leading radix:
	else if(new_runlength || (rad0save != radix0)) {
		rad0save = radix0;		// Both radix 0 and nchunks get checked in the dyadic_square routines and only the needed inits done
		init_sse2 = nchunks;	// Use *value* of init_sse2 to store #threads
		thr_id = -1;
		/* The dyadic-square routines need a few more params passed in init-mode than do the standalone FFT-pass routines: */
		radix16_dyadic_square(0x0, arr_scratch, n, radix0, rt0, rt1, 0, nradices_prim, radix_prim, 0, init_sse2,thr_id);
		radix32_dyadic_square(0x0, arr_scratch, n, radix0, rt0, rt1, 0, nradices_prim, radix_prim, 0, init_sse2,thr_id);
	}
	new_runlength = FALSE;

	/* end of initialization sequence.	*/

/**********************************************************************/

	/*...Init clock counter:	*/
	ASSERT(HERE, tdiff != 0,"fermat_mod_square.c: tdiff != 0");

#ifdef CTIME
	clock1 = clock();
#else
//	clock1 = time(0x0);
	clock1 = getRealTime();
#endif

	*tdiff = 0.0;

	/*...At the start of each iteration cycle, need to forward-weight the array of integer residue digits.
	For the non-power-of-2 case, IBDWT weights must be applied to the non-acyclic-twisted transform inputs,
	i.e. first apply the IBDWT weights, then do the acyclic-twisting complex multiply:
	*/
	// Mar 2017: Can skip this step if it's the start of a production test (note that any initial-residue shift
	// in such cases is handled via single-array-word forward-DWT-weighting in the Mlucas.c shift_word() function),
	// but need it if add RNG-input-setting above for debug, hence also check a[1] for nonzero:
	if(ilo || a[1]) {
		if(!pow2_fft) {
			ii = 0;	/* index into wt0 array (mod NWT) is here: */
			/* Evens: */
			for(j = 0; j < n; j += 2)
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
				j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				wt = wt0[ii];
				a[j1] *= wt;
				ii += SW_DIV_N - nwt;
				ii += ( (-(int)((uint32)ii >> 31)) & nwt);\
			}
			/* Odds: */
			ASSERT(HERE, ii == 0,"fermat_mod_square.c: ii == 0");
			for(j = 0; j < n; j += 2)
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
				j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1 + RE_IM_STRIDE;
				wt = wt0[ii];
				a[j2] *= wt;
				ii += SW_DIV_N - nwt;
				ii += ( (-(int)((uint32)ii >> 31)) & nwt);\
			}
		}
		/* Acyclic twisting: */
		for(j = 0; j < n; j += 2)
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
			j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
			j2 = j1 + RE_IM_STRIDE;
			/* Get the needed Nth root of -1: */
			l = (j >> 1);	/* j/2 */
			k1=(l & NRTM1);
			k2=(l >> NRT_BITS);
			t1=rn0[k1].re;		t2=rn0[k1].im;
			rt=rn1[k2].re;		it=rn1[k2].im;
			re =t1*rt-t2*it;	im =t1*it+t2*rt;
			/* Do the cyclic -> acyclic weighting: */
			t1    = a[j1]*re - a[j2]*im;
			a[j2] = a[j2]*re + a[j1]*im;
			a[j1] = t1;
		}
	}

/*...and perform the initial pass of the forward transform.	*/

/*...NOTE: If the first radix to be processed is 2, 4 or 8, it is assumed that a power-of-2 FFT is being performed,
     hence no small-prime version of the corresponding pass1 routines is needed.	*/

	switch(radix0)
	{
	case 5 :
		 radix5_dif_pass1(a,n); break;
	case 6 :
		 radix6_dif_pass1(a,n); break;
	case 7 :
		 radix7_dif_pass1(a,n); break;
	case 8 :
		 radix8_dif_pass1(a,n); break;
	case 9 :
		 radix9_dif_pass1(a,n); break;
	case 10 :
		radix10_dif_pass1(a,n); break;
	case 11 :
		radix11_dif_pass1(a,n); break;
	case 12 :
		radix12_dif_pass1(a,n); break;
	case 13 :
		radix13_dif_pass1(a,n); break;
	case 14 :
		radix14_dif_pass1(a,n); break;
	case 15 :
		radix15_dif_pass1(a,n); break;
	case 16 :
		radix16_dif_pass1(a,n); break;
	case 18 :
		radix18_dif_pass1(a,n); break;
	case 20 :
		radix20_dif_pass1(a,n); break;
	case 22 :
		radix22_dif_pass1(a,n); break;
	case 24 :
		radix24_dif_pass1(a,n); break;
	case 26 :
		radix26_dif_pass1(a,n); break;
	case 28 :
		radix28_dif_pass1(a,n); break;
	case 30 :
		radix30_dif_pass1(a,n); break;
	case 32 :
		radix32_dif_pass1(a,n); break;
	case 36 :
		radix36_dif_pass1(a,n); break;
	case 56 :
		radix56_dif_pass1(a,n); break;
	case 60 :
		radix60_dif_pass1(a,n); break;
	case 63 :
		radix63_dif_pass1(a,n); break;
	case 64 :
		radix64_dif_pass1(a,n); break;
	/*
	case 72 :
		radix72_dif_pass1(a,n); break;
	case 112:
		radix112_dif_pass1(a,n); break;
	case 120:
		radix120_dif_pass1(a,n); break;
	*/
	case 128:
		radix128_dif_pass1(a,n); break;
	/*
	case 144:
		radix144_dif_pass1(a,n); break;
	*/
	case 224:
		radix224_dif_pass1(a,n); break;
	case 240 :
		radix240_dif_pass1(a,n); break;
	case 256 :
		radix256_dif_pass1(a,n); break;
	case 512 :
		radix512_dif_pass1(a,n); break;
	case 960 :
		radix960_dif_pass1(a,n); break;
	case 992 :
		radix992_dif_pass1(a,n); break;
	case 1008:
		radix1008_dif_pass1(a,n); break;
	case 1024:
		radix1024_dif_pass1(a,n); break;
	case 4032:
		radix4032_dif_pass1(a,n); break;
/*
	case 4096:
		radix4096_dif_pass1(a,n); break;
*/
	default :
		sprintf(cbuf,"FATAL: radix %d not available for dif_pass1. Halting...\n",radix0);
		fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
	}

/**********************************************************************/

	/* Main iteration loop is here. Do forward-FFT/pointwise-square/inverse-FFT, inverse weighting,
	carry propagation, fractional error checking and forward weighting in same loop:
	*/
	ierr = 0;	/* Any return-value error code (whether fatal or not) stored here */

	ASSERT(HERE, ihi > ilo,"ferm_mod_square.c: ihi <= ilo!");

#ifdef MULTITHREAD

  #ifdef USE_PTHREAD
	// Jan 2018: To support PRP-cofactor testing, must allow user to feed different input vectors
	// to this modpow routine, so move this struct-element init of Main data array pointer here:
	for(i = 0; i < radix0; ++i) { tdat[i].arrdat = a; }

  #elif(defined(USE_OMP))
	// OpenMP currently not supported - *** To-Do: Port mers_mod_square OpenMP stuff to here ****.
	#error OpenMP currently not supported - Please recheck thread-related defines in platform.h file.
	omp_set_num_threads(NTHREADS);
 	for(i=0; i < NTHREADS; i++) { num_chunks[i] = 0; }
  #else
	#error MULTITHREAD defined but USE_PTHREAD not - Please recheck thread-related defines in platform.h file.
  #endif

  #if DBG_THREADS
	fprintf(stderr,"%s: NTHREADS = %3d\n",func,NTHREADS);
  #endif

#endif

for(iter=ilo+1; iter <= ihi && MLUCAS_KEEP_RUNNING; iter++)
{

/*...perform the FFT-based squaring:
     Do last S-1 of S forward decimation-in-frequency transform passes.	*/

	/* Break the remaining portion of the FFT into radix0 blocks, each of which ideally
	should operate on a dataset which fits entirely into the L2 cache of the host machine.
	In a multithreaded implementation, process NTHREADS blocks in parallel fashion:
	If the number of available processor cores does not divide radix0, there will be one or more under-or-unutilized CPUs.
	*/

#ifdef MULTITHREAD

  #ifdef USE_PTHREAD

	/* create radix0 new threads each of which will execute 'fermat_process_chunk()' over some specified index subrange.
	In order to match the threads executing at any given time to the available CPUs, divide the thread execution
	into [NTHREADS] 'work shifts' ( <= #CPus), each with its threads starting and completing their work before the next shift
	comes online:
	*/
	if(NTHREADS > 0) {	/******* Change 0 --> to test thread-team/join overhead ******/
		isum = 0;

	// MacOS does weird things with threading (e.g. Idle" main thread burning 100% of 1 CPU)
	// so on that platform try to be clever and interleave main-thread and threadpool-work processing
	#if 0//def OS_TYPE_MACOSX

		if(NTHREADS > 1) {
			for(thr_id = 0; thr_id < pool_work_units; ++thr_id)
			{
			/*** As does each of the (NTHREADS - 1) pool threads: ***/
				task_control.data = (void*)(&tdat[thr_id]);
			//	printf("adding pool task %d\n",thr_id);
				threadpool_add_task(tpool, &task_control, task_is_blocking);
			}

			/*** Main execution thread executes remaining chunks in serial fashion (but in || with the pool threads): ***/
			for(j = 0; j < main_work_units; ++j)
			{
			//	printf("adding main task %d\n",j + pool_work_units);
				fermat_process_chunk( (void*)(&tdat[j + pool_work_units]) );
			}

			struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
			ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
			ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

			while(tpool->free_tasks_queue.num_tasks != pool_work_units) {
				// Finer-resolution, declared in <time.h>; cf. http://linux.die.net/man/2/nanosleep
				ASSERT(HERE, 0 == nanosleep(&ns_time, 0x0), "nanosleep fail!");
			}
		} else {
			for(thr_id = 0; thr_id < radix0; ++thr_id)
			{
				fermat_process_chunk( (void*)(&tdat[thr_id]) );
			}
		}
//printf("Iter = %d\n",iter);
		for(thr_id = 0; thr_id < radix0; ++thr_id)
		{
		#ifdef DBG_TIME
			dt_fwd += tdat[thr_id].dt_fwd;
			dt_inv += tdat[thr_id].dt_inv;
			dt_sqr += tdat[thr_id].dt_sqr;
			dt_tot += tdat[thr_id].dt_tot;
		#endif
		}

	#elif defined(USE_THREADPOOL)	// Threadpool-based dispatch for generic (non OS X) Linux

			for(thr_id = 0; thr_id < pool_work_units; ++thr_id)
			{
				task_control.data = (void*)(&tdat[thr_id]);
			// Get thread ID w.r.to the pool:
			//	j = ((struct thread_init *)(&task_control)->data)->thread_num;
			//	printf("adding pool task %d with pool ID [%d]\n",thr_id,j);
				threadpool_add_task(tpool, &task_control, task_is_blocking);
			//	printf("; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);
			}

		//	printf("start; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);
			struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
			ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
			ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

		//	while(tpool->tasks_queue.num_tasks != 0) {	//*** not safe, since can have #tasks == 0 with some tasks still in flight ***
			while(tpool->free_tasks_queue.num_tasks != pool_work_units) {
			//		sleep(1);	//*** too granular ***
				// Finer-resolution, declared in <time.h>; cf. http://linux.die.net/man/2/nanosleep
				ASSERT(HERE, 0 == nanosleep(&ns_time, 0x0), "nanosleep fail!");
			//	printf("sleep; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);
			}
		//	printf("end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

		for(thr_id = 0; thr_id < radix0; ++thr_id)
		{
		#ifdef DBG_TIME
			dt_fwd += tdat[thr_id].dt_fwd;
			dt_inv += tdat[thr_id].dt_inv;
			dt_sqr += tdat[thr_id].dt_sqr;
			dt_tot += tdat[thr_id].dt_tot;
		#endif
		}

	#elif 1	// Simple-and-stupid thread-dispatch:

		// Try a single shifts with one thread for each work chunk, let OS deal with dispatching 'em:
		for(thr_id = 0; thr_id < radix0; ++thr_id)
		{
			rc = pthread_create(&thread[thr_id], &attr, fermat_process_chunk, (void*)(&tdat[thr_id]));
			if (rc) {
				printf("ERROR; return code from pthread_create() is %d\n", rc);
				exit(-1);
			}
		}
		// As each thread finishes, add its result into an accumulator in non-blocking fashion (no mutexes needed):
		for(thr_id = 0; thr_id < radix0; ++thr_id)
		{
			rc = pthread_join(thread[thr_id], &thr_status);
			if (rc) {
				printf("ERROR; return code from pthread_join() is %d\n", rc);
				exit(-1);
			}
		//	printf("Main: completed join with thread %ld having a status of %ld\n",thr_id,(long)thr_status);
			isum += thr_ret[thr_id];
		}
		ASSERT(HERE, isum == 0, "Nonzero thread-team return checksum!");

	#else

		nshift = radix0/NTHREADS;	// Number of shifts with one thread for each CPU
		for(j = 0; j < nshift; ++j)
		{
			ioffset = j*NTHREADS;
			for(i = 0; i < NTHREADS; ++i)
			{
				thr_id = i+ioffset;
				rc = pthread_create(&thread[thr_id], &attr, fermat_process_chunk, (void*)(&tdat[thr_id]));
				if (rc) {
					printf("ERROR; return code from pthread_create() is %d\n", rc);
					exit(-1);
				}
			}
			// As each thread finishes, add its result into an accumulator in non-blocking fashion (no mutexes needed):
			for(i = 0; i < NTHREADS; ++i)
			{
				thr_id = i+ioffset;
				rc = pthread_join(thread[thr_id], &thr_status);
				if (rc) {
					printf("ERROR; return code from pthread_join() is %d\n", rc);
					exit(-1);
				}
			//	printf("Main: completed join with thread %ld having a status of %ld\n",thr_id,(long)thr_status);
				isum += thr_ret[thr_id];
			}
			ASSERT(HERE, isum == 0, "Nonzero thread-team return checksum!");
		}

	#endif
	} else {
		/* Single-threaded version: */
		for(i = 0; i < radix0; ++i)
		{
			fermat_process_chunk((void*)(&tdat[i]));
		#ifdef DBG_TIME
			dt_fwd += tdat[i].dt_fwd;
			dt_inv += tdat[i].dt_inv;
			dt_sqr += tdat[i].dt_sqr;
			dt_tot += tdat[i].dt_tot;
		#endif
		}
	}

  #endif

#else

	/* Unthreaded version: */
	for(ii = 0; ii < radix0; ++ii)
	{
		fermat_process_chunk(
			a,arr_scratch,n,rt0,rt1,index,ii,nradices_prim,radix_prim
		#ifdef DBG_TIME
			,&dt_fwd,&dt_inv,&dt_sqr,&dt_tot
		#endif
		);
	}

#endif	// threaded?

	// Update RES_SHIFT via mod-doubling-and-add-random-bit. If initial RES_SHIFT = 0 the 'random' bit array = 0, so RES_SHIFT remains 0:
	ASSERT(HERE, RES_SHIFT == 0, "Shifted residues currently unsupported for Fermat-mod!");	// Jul 2018: shifted-res works fine for Fermat-mod
								// until res-array gets close to full, then hit bizarre issue of shifted-run residue "going negative" for a few
								// iters before diverging completely). Example: -f 24 -fftlen 896 -radset 4, which I used as my first-debug case
								// here (after suitably modifying the radix-28 carry routine and the Fermat-mod carry macros to take a prp_mult
								// argument) residues ,atch thru iter 20, then iter 21,22 the shifted-run residue array is the negative of the
								// unshifted-run one, then starting at iter 23 things wholly diverge. Need to revisit later - maybe some aspect
								// of the right-angle transform in the presence of a shifted residue scheme?
	i = (iter % ITERS_BETWEEN_CHECKPOINTS) - 1;	// Bit we need to read...iter-counter is unit-offset w.r.to iter-interval, hence the -1
	MOD_ADD64(RES_SHIFT,RES_SHIFT,p,RES_SHIFT);
	RES_SHIFT += ((BASE_MULTIPLIER_BITS[i>>6] >> (i&63)) & 1);	// No mod needed on this add, since result of pvs line < p
//	printf(" shift = %d, ",RES_SHIFT);fflush(stdout);
/*...Do the final inverse FFT pass, carry propagation and initial forward FFT pass in one fell swoop, er, swell loop...	*/

	fracmax = 0.0;
//printf("Exit(0) from %s\n",func); exit(0);
	switch(radix0)
	{
		case  5 :
			ierr =  radix5_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case  6 :
			ierr =  radix6_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case  7 :
			ierr =  radix7_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case  8 :
			ierr =  radix8_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case  9 :
			ierr =  radix9_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 10 :
			ierr = radix10_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 11 :
			ierr = radix11_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 12 :
			ierr = radix12_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 13 :
			ierr = radix13_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 14 :
			ierr = radix14_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 15 :
			ierr = radix15_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 16 :
			ierr = radix16_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 18 :
			ierr = radix18_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 20 :
			ierr = radix20_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 22 :
			ierr = radix22_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 24 :
			ierr = radix24_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 26 :
			ierr = radix26_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 28 :
			ierr = radix28_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 30 :
			ierr = radix30_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 32 :
			ierr = radix32_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 36 :
			ierr = radix36_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 56 :
			ierr = radix56_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 60 :
			ierr = radix60_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 63 :
			ierr = radix63_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 64 :
			ierr = radix64_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 128 :
			ierr = radix128_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 224 :
			ierr = radix224_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 240 :
			ierr = radix240_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 256 :
			ierr = radix256_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 512 :
			ierr = radix512_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 960 :
			ierr = radix960_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 992 :
			ierr = radix992_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 1008:
			ierr = radix1008_ditN_cy_dif1    (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 1024:
			ierr = radix1024_ditN_cy_dif1    (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 4032:
			ierr = radix4032_ditN_cy_dif1    (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
	/*
		case 4096:
			ierr = radix4096_ditN_cy_dif1    (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
	*/
		default :
			sprintf(cbuf,"FATAL: radix %d not available for ditN_cy_dif1. Halting...\n",radix0); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
	}

	/* Nonzero remaining carries are instantly fatal: */
	if(ierr)
		return(ierr);
	/* Update Max. Max. Error: */
	if(fracmax > MME)
		MME  = fracmax;
	/* Accumulate Avg. Max. Error: */
	if(iter > AME_ITER_START)
		AME += fracmax;

/*...Now do the fractional error check. Any fractional part  > 0.40625 generates a warning...	*/
// Dec 2014: Bump threshold up from ( >= 0.4 ) to ( > 0.40625 ):
	if(fracmax > 0.40625)
	{
		sprintf(cbuf, "F%u Roundoff warning on iteration %8u, maxerr = %16.12f\n",findex,iter,fracmax);

	/*...Fractional parts close to 0.5 cause the program to quit.
		 We put the interactive-mode errlimit close to 0.5, to let people really push the limits if they want to...	*/
		if(INTERACT)
		{
			fprintf(stderr,"%s",cbuf);
			if(fracmax > 0.47 )
			{
				fprintf(stderr," FATAL ERROR...Halting test of F%u\n",findex);
				ierr = ERR_ROUNDOFF;
				return(ierr);
			}
		}
		else
		{
			fp = mlucas_fopen(   OFILE,"a");
			fq = mlucas_fopen(STATFILE,"a");
			fprintf(fp,"%s",cbuf);
			fprintf(fq,"%s",cbuf);
			if (scrnFlag)	/* Echo output to stddev */
			{
				fprintf(stderr,"%s",cbuf);
			}

			// In range test mode, any fractional part >= 0.4375 is cause for error exit:
			if(fracmax >= 0.4375 )
			{
			/* In order to cope with occasional spurious error - especially in multithreaded and GPU settings - need a mechanism
			to retry the current iteration interval starting from the last savefile and, based on the result, act accordingly.
			Here is the retry scheme in words - there are 4 distinct cases, divided into 2 main cases, the second of which has
			3 subcases:

			[1] [condition] (ROE_ITER == 0) : This is the initial attempt of the current interval.

				[action] Store the iteration and value (fracmax) of the initial-pass error in the (signed int)ROE_ITER and
				(double)ROE_VAL globals (both of which are inited to 0), then (via suitable logic in the calling routine)
				rerun the interval with the same FFT params:
				*/
				if(ROE_ITER == 0) {
					sprintf(cbuf," Retrying iteration interval to see if roundoff error is reproducible.\n");
					ROE_ITER = iter;
					ROE_VAL = fracmax;
				}
			/*
			[2] [condition] (ROE_ITER > 0) : This is a retry of the current interval. There are 3 subcases:
			*/
				else if(ROE_ITER > 0) {
				/*
				[2a] [condition] The interval-retry is successful, i.e. suffers no fatal ROE.
					In this case we will not hit this code snip, but rather will arrive at the return() with iter = ihi+1,
					hence the code handling for this condition is to be found there.

					[action] Prior to returning, print a "retry successful" informational and rezero ROE_ITER and ROE_VAL.
				*/
				/*
				[2b] [condition] The error is reproducible, i.e. recurs on the same iteration, with the same fractional value.

					[action] Negate ROE_ITER as a signal to the caller that the interval needs (and the rest of the run)
					needs to be with the next-larger available FFT length. Since ROE_VAL is no longer needed, rezero it,
					and expect the caller to do the same with ROE_ITER:
					*/
					if(ROE_ITER == iter && ROE_VAL == fracmax) {
						sprintf(cbuf," Roundoff error is reproducible ... switching to next-larger available FFT length and retrying.\n");
						ROE_ITER = -ROE_ITER;
						ROE_VAL = 0.0;
					}
				/*
				[2c] [condition; rare] The error is not reproducible, but encountered a different ROE in the retry of the interval.

					[action] Originally treated as fatal (usually a sign of data corruption), but switching to next-larger FFT length
					seems to work in most such cases, since it triggers an FFT-data-tables-reinit. Can run one or more iteration
					intervals at the longer length, then switch back:
					*/
					else {
						sprintf(cbuf," The error is not reproducible, encountered a different fatal ROE in interval-retry ... note this is an indicator of possible data corruption. Switching to next-larger FFT length, or next-smaller, if currently running at larger-than-default FFT length.\n");
						ROE_ITER = -ROE_ITER;
						ROE_VAL = fracmax;	// Use ROE_VAL-zero-or-not to distinguish from above case
					}
				} else if(ROE_ITER < 0) {
					sprintf(cbuf,"Unexpected condition (ROE_ITER < 0) in %s ... quitting. Please restart the program at your earliest convenience.\n",func);
					ierr = ERR_UNKNOWN_FATAL;
				}

				fprintf(fp,"%s",cbuf);
				fprintf(fq,"%s",cbuf);
				fclose(fp);	fp = 0x0;
				fclose(fq);	fq = 0x0;
				if (scrnFlag)	/* Echo output to stddev */
				{
					fprintf(stderr,"%s",cbuf);
				}
				ierr = ERR_ROUNDOFF;
				return(ierr);
			}
			fclose(fp);	fp = 0x0;
			fclose(fq);	fq = 0x0;
		}
	}

	//...Whew - that"s a lot of stuff that just happened. Luckily, computer chips don't understand the concept of "Miller time."

	// Accumulate the cycle count in a floating double on each pass to avoid problems
	// with integer overflow of the clock() result, if clock_t happens to be 32-bit int on the host platform:
#ifdef CTIME
	clock2 = clock();
	*tdiff += (double)(clock2 - clock1);
	clock1 = clock2;
#endif
	// Listen for interrupts:
	if (signal(SIGINT, sig_handler) == SIG_ERR)
		fprintf(stderr,"Can't catch SIGINT.\n");
	else if (signal(SIGTERM, sig_handler) == SIG_ERR)
		fprintf(stderr,"Can't catch SIGTERM.\n");
	else if (signal(SIGHUP, sig_handler) == SIG_ERR)
		fprintf(stderr,"Can't catch SIGHUP.\n");
}	/* End of main for(iter....) loop	*/

// On early-exit-due-to-interrupt, decrement iter since we didn't actually do the (iter)th iteration
if(!MLUCAS_KEEP_RUNNING) iter--;
if(iter < ihi) {
	ASSERT(HERE, !MLUCAS_KEEP_RUNNING, "Premature iteration-loop exit due to unexpected condition!");
	ierr = ERR_INTERRUPT;
	ROE_ITER = iter;	// Function return value used for error code, so save number of last-iteration-completed-before-interrupt here
}

#ifdef RTIME
//	clock2 = time(0x0);
//	*tdiff += difftime(clock2 , clock1);
	clock2 = getRealTime();
	*tdiff += clock2 - clock1;
#endif

#ifdef DBG_TIME
	printf("fermat_process_chunk times: total = %10.5f, fwd = %10.5f, inv = %10.5f, sqr = %10.5f\n", dt_tot*ICPS, dt_fwd*ICPS, dt_inv*ICPS, dt_sqr*ICPS);
#endif

#if DBG_THREADS
	fprintf(stderr,"%s: #Chunks processed by each thread: ",func);
	for(i=0; i < NTHREADS; i++)
	{
		fprintf(stderr,"%d[%d] ", i, num_chunks[i]);
	}
	fprintf(stderr,"\n");
#endif

/**********************************************************************/

/*...At the end of each iteration cycle, need to undo the initial DIF FFT pass...	*/

	switch(radix0)
	{
	case  5 :
		 radix5_dit_pass1(a,n); break;
	case  6 :
		 radix6_dit_pass1(a,n); break;
	case  7 :
		 radix7_dit_pass1(a,n); break;
	case  8 :
		 radix8_dit_pass1(a,n); break;
	case  9 :
		 radix9_dit_pass1(a,n); break;
	case 10 :
		radix10_dit_pass1(a,n); break;
	case 11 :
		radix11_dit_pass1(a,n); break;
	case 12 :
		radix12_dit_pass1(a,n); break;
	case 13 :
		radix13_dit_pass1(a,n); break;
	case 14 :
		radix14_dit_pass1(a,n); break;
	case 15 :
		radix15_dit_pass1(a,n); break;
	case 16 :
		radix16_dit_pass1(a,n); break;
	case 18 :
		radix18_dit_pass1(a,n); break;
	case 20 :
		radix20_dit_pass1(a,n); break;
	case 22 :
		radix22_dit_pass1(a,n); break;
	case 24 :
		radix24_dit_pass1(a,n); break;
	case 26 :
		radix26_dit_pass1(a,n); break;
	case 28 :
		radix28_dit_pass1(a,n); break;
	case 30 :
		radix30_dit_pass1(a,n); break;
	case 32 :
		radix32_dit_pass1(a,n); break;
	case 36 :
		radix36_dit_pass1(a,n); break;
	case 56 :
		radix56_dit_pass1(a,n); break;
	case 60 :
		radix60_dit_pass1(a,n); break;
	case 63 :
		radix63_dit_pass1(a,n); break;
	case 64 :
		radix64_dit_pass1(a,n); break;
/*
	case 72 :
		radix72_dit_pass1(a,n); break;
	case 112:
		radix112_dit_pass1(a,n); break;
	case 120:
		radix120_dit_pass1(a,n); break;
*/
	case 128:
		radix128_dit_pass1(a,n); break;
/*
	case 144:
		radix144_dit_pass1(a,n); break;
*/
	case 224 :
		radix224_dit_pass1(a,n); break;
	case 240 :
		radix240_dit_pass1(a,n); break;
	case 256 :
		radix256_dit_pass1(a,n); break;
	case 512 :
		radix512_dit_pass1(a,n); break;
	case 960 :
		radix960_dit_pass1(a,n); break;
	case 992 :
		radix992_dit_pass1(a,n); break;
	case 1008 :
		radix1008_dit_pass1(a,n); break;
	case 1024:
		radix1024_dit_pass1(a,n); break;
	case 4032 :
		radix4032_dit_pass1(a,n); break;
/*
	case 4096:
		radix4096_dit_pass1(a,n); break;
*/
	default :
		sprintf(cbuf,"FATAL: radix %d not available for dit_pass1. Halting...\n",radix0); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
	}

/*...and unweight the data array.	*/

	/*...At the end of each iteration cycle, need to inverse-weight the array of integer residue digits.
	For the non-power-of-2 case, IBDWT weights must be applied to the non-acyclic-twisted transform inputs,
	i.e. first do the acyclic-untwisting complex multiply, then apply the inverse IBDWT weights:
	*/
	max_fp = 0.0;
	if(!pow2_fft)
	{
		ii = 0;	/* index into wt1 array (mod NWT) is here: */
		/* Evens: */
		for(j = 0; j < n; j += 2)
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
			j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */

			wtinv = wt1[ii];
			a[j1] *= wtinv;
			ii += SW_DIV_N - nwt;
			ii += ( (-(int)((uint32)ii >> 31)) & nwt);\
		}
		/* Odds: */
		ASSERT(HERE, ii == 0,"fermat_mod_square.c: ii == 0");
		for(j = 0; j < n; j += 2)
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
			j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
			j2 = j1 + RE_IM_STRIDE;

			wtinv = wt1[ii];
			a[j2] *= wtinv;
			ii += SW_DIV_N - nwt;
			ii += ( (-(int)((uint32)ii >> 31)) & nwt);\
		}
	}

	/* Acyclic untwisting: */
	for(j = 0; j < n; j += 2)
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
		j1 = j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

		/* Get the needed Nth root of -1: */
		l = (j >> 1);	/* j/2 */
		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		t1=rn0[k1].re;		t2=rn0[k1].im;
		rt=rn1[k2].re;		it=rn1[k2].im;
		re =t1*rt-t2*it;	im =t1*it+t2*rt;

		/* Do the cyclic -> acyclic inverse weighting: */
		t1    = a[j1]*re + a[j2]*im;
		a[j2] = a[j2]*re - a[j1]*im;
		a[j1] = t1;

		atmp  = a[j1]*radix_inv;
		a[j1] = NINT(atmp);
		frac_fp = fabs(a[j1]-atmp);
		if(frac_fp > max_fp)
			max_fp = frac_fp;

		atmp  = a[j2]*radix_inv;
		a[j2] = NINT(atmp);
		frac_fp = fabs(a[j2]-atmp);
		if(frac_fp > max_fp)
			max_fp = frac_fp;
	}
	if(max_fp > 0.01)
	{
		fprintf(stderr,"%s: max_fp > 0.01! Value = %20.10f\n",func,max_fp);
		fprintf(stderr,"Check your build for inadvertent mixing of SIMD build modes!\n");
		ASSERT(HERE, max_fp < 0.01,"fermat_mod_square.c: max_fp < 0.01");
	}

	if(ierr == ERR_INTERRUPT) {	// In this need to bypass [2a] check below because ROE_ITER will be set to last-iteration-done
		return(ierr);
	}
	// Cf. [2a] above: The interval-retry is successful, i.e. suffers no fatal ROE.
	// [action] Prior to returning, print a "retry successful" informational and rezero ROE_ITER and ROE_VAL.
	if(ROE_ITER > 0) {
		ASSERT(HERE, (ierr == 0) && (iter = ihi+1), "[2a] sanity check failed!");
		ROE_ITER = 0;
		ROE_VAL = 0.0;
		fp = mlucas_fopen(   OFILE,"a");
		fq = mlucas_fopen(STATFILE,"a");
		sprintf(cbuf,"Retry of iteration interval with fatal roundoff error was successful.\n");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;
		if(scrnFlag)	/* Echo output to stddev */
		{
			fprintf(stderr,"%s",cbuf);
		}
	}

	return(ierr);
}

/***************/

#if(defined(MULTITHREAD) && defined(USE_PTHREAD))

void*
fermat_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
{
	struct ferm_thread_data_t* thread_arg = targ;
	int ii = thread_arg->tid, thr_id = ii;	// ii-value same as unique thread identifying number
//	int thr_id = ii%NTHREADS;	// 'thread ID' for calls to various subroutines is w.r.to current thread set, i.e. unique thread-ID modulo #threads-per-shift.
	double *a           = thread_arg->arrdat;
	int *arr_scratch    = thread_arg->arr_scratch;
	int n               = thread_arg->n;
	struct complex *rt0 = thread_arg->rt0;
	struct complex *rt1 = thread_arg->rt1;
	int*index           = thread_arg->index;
	int nradices_prim   = thread_arg->nradices_prim;
	int*radix_prim      = thread_arg->radix_prim;
#ifdef DBG_TIME
	double *dt_fwd = &(thread_arg->dt_fwd);
	double *dt_inv = &(thread_arg->dt_inv);
	double *dt_sqr = &(thread_arg->dt_sqr);
	double *dt_tot = &(thread_arg->dt_tot);
#endif

#else

void fermat_process_chunk(double a[], int arr_scratch[], int n, struct complex rt0[], struct complex rt1[], int index[], int ii, int nradices_prim, int radix_prim[]
#ifdef DBG_TIME
	, double *dt_fwd, double *dt_inv, double *dt_sqr, double *dt_tot
#endif
)
{
	int thr_id = 0;	/* In unthreaded mode this must always = 0 */

#endif	// #ifdef MULTITHREAD

	const char func[] = "fermat_process_chunk";
#ifdef DBG_TIME
	clock_t clock0, clock1, clock2, clock3;
	clock0 = clock();
#endif
	int radix0 = RADIX_VEC[0];
    int i,incr,istart,jstart,k,koffset,l,mm;
	int init_sse2 = FALSE;	// Init-calls to various radix-pass routines presumed done prior to entry into this routine

	l = ii;
	k    = 0;
	mm   = 1;
	// Starting location of current data-block-to-be-processed within A-array:
	incr = n/radix0;
	istart = l*incr;
	jstart = istart + ((istart >> DAT_BITS) << PAD_BITS );

	for(i=1; i <= NRADICES-2; i++)
	{
		/* Offset from base address of index array = L*NLOOPS = L*MM : */
		koffset = l*mm;
		switch(RADIX_VEC[i])
		{
		case  8 :
			 radix8_dif_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		case 16 :
			radix16_dif_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		case 32 :
			radix32_dif_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		default :
			sprintf(cbuf,"FATAL: radix %d not available for dif_pass. Halting...\n",RADIX_VEC[i]); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}
		k    += mm*radix0;
		mm   *= RADIX_VEC[i];
		incr /= RADIX_VEC[i];
	}	/* end i-loop. */

#ifdef DBG_TIME
	clock1 = clock();
	*dt_fwd += (double)(clock1 - clock0);
#endif
	/*...Final DIF pass, dyadic squaring and initial DIT pass are all done via a fused 1-pass procedure: */
	koffset = l*mm;
	/* The roots-block-re-use param mm not needed for innermost pass, since there each set of inputs gets its own set of roots: */
	switch(RADIX_VEC[NRADICES-1])
	{
		case 16 :
			radix16_dyadic_square(&a[jstart],arr_scratch,n,radix0,rt0,rt1,ii,nradices_prim,radix_prim,incr,init_sse2,thr_id); break;
		case 32 :
			radix32_dyadic_square(&a[jstart],arr_scratch,n,radix0,rt0,rt1,ii,nradices_prim,radix_prim,incr,init_sse2,thr_id); break;
	#if 0
		/*
		In synthesized-large-final-pass-radix scheme with penultimate radix R added to the supported final-pass
		radix 16 or 32, we divide incr /= R as usual to properly set the data-access stride, but we DO NOT increment
		mm *= R as normal - instead increment mm based on desired L1-cache-target-size and call DFT-pass routine as
		many times as needed to process same total #data as in the above radix-16 and radix-32 passes, and then wrap
		the R-DIF/16-or-32-DYAD-SQUARE/R-DIT sequence of the synthesized approach in a loop, inside of which we incr
		jstart and adjust koffset as needed between these mini-pass function calls:

		NOTE: have already done a final set of updates: k += mm*radix0, mm *= RADIX_VEC[i], incr /= RADIX_VEC[i]
		prior to exit from the above multi-DIF-pass loop, we need to UNDO those first, noting that [i] needs to be
		decremented -= 1 first to undo final incrementing which triggers loop exit:
		*/
		case 256 :
			koffset = l*mm;
			int jdx = istart,jpad,loop;	// Init jdx = istart (unpadded input index) here since will be successively adding +=incr and padding the result within the loop below (i.e. to avoid double-padding the starting offset)
		  for(loop = 0; loop < mm; loop++, jdx += incr, koffset++) {
			jpad = jdx + ( (jdx >> DAT_BITS) << PAD_BITS );
			radix16_dif_pass     (&a[jpad],n,rt0,rt1,&index[k+koffset],1,incr,init_sse2,thr_id);
			radix16_dyadic_square(&a[jpad],arr_scratch,n,radix0,rt0,rt1,ii,nradices_prim,radix_prim,incr>>4,init_sse2,thr_id);
			radix16_dit_pass     (&a[jpad],n,rt0,rt1,&index[k+koffset],1,incr,init_sse2,thr_id);
		  }
			break;
		case 1024 :
			koffset = l*mm;
			jdx = istart;
		  for(loop = 0; loop < mm; loop++, jdx += incr, koffset++) {
			jpad = jdx + ( (jdx >> DAT_BITS) << PAD_BITS );
			radix32_dif_pass     (&a[jpad],n,rt0,rt1,&index[k+koffset],1,incr,init_sse2,thr_id);
			radix32_dyadic_square(&a[jpad],arr_scratch,n,radix0,rt0,rt1,ii,nradices_prim,radix_prim,incr>>5,init_sse2,thr_id);
			radix32_dit_pass     (&a[jpad],n,rt0,rt1,&index[k+koffset],1,incr,init_sse2,thr_id);
		  }
			break;
	#endif
		default :
			sprintf(cbuf,"FATAL: radix %d not available for wrapper/square. Halting...\n",RADIX_VEC[NRADICES-1]); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
	}

#ifdef DBG_TIME
	clock2 = clock();
	*dt_sqr += (double)(clock2 - clock1);
#endif

	/*...Rest of inverse decimation-in-time (DIT) transform. Note that during IFFT we process the radices in reverse
     order. The first array sent to each pass routine is assumed to contain the bit-reversed floating data.	*/

	/* Quick-n-dirty way of generating the correct starting values of k, mm and incr -
	simply use the skeleton of the forward (DIF) loop, sans the i = NRADICES-2 pass
	(since k, mm and incr are post-incremented):
	*/
	k    = 0;
	mm   = 1;
	incr = n/radix0;

	for(i=1; i <= NRADICES-2; i++)
	{
		k    += mm*radix0;
		mm   *= RADIX_VEC[i];
		incr /= RADIX_VEC[i];
	}

	/* Now do the DIT loop, running the radices (and hence the values of k, mm and incr) in reverse: */

	for(i=NRADICES-2; i >= 1; i--)
	{
		incr *= RADIX_VEC[i];
		mm   /= RADIX_VEC[i];
		k    -= mm*radix0;
		koffset = l*mm;
		switch(RADIX_VEC[i])
		{
		case  8 :
			 radix8_dit_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		case 16 :
			radix16_dit_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		case 32 :
			radix32_dit_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		default :
			sprintf(cbuf,"FATAL: radix %d not available for dit_pass. Halting...\n",RADIX_VEC[i]); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}
	}	/* end i-loop */

#ifdef DBG_TIME
	clock3 = clock();
	*dt_inv += (double)(clock3 - clock2);
	*dt_tot += (double)(clock3 - clock0);
#endif

#ifdef MULTITHREAD

	*(thread_arg->retval) = 0;	// 0 indicates successful return of current thread
//	printf("Return from Thread %d ... ", ii);
  #ifdef USE_THREADPOOL
	return 0x0;
  #else
	pthread_exit(NULL);
  #endif

#else
	return;
#endif
}

