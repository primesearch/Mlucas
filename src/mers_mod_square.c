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
	#else
		#error Requires threadpool!
	#endif

	struct mers_thread_data_t{
		int tid;
		int*retval;
		double*arrdat;			// Main data array
	#ifdef USE_FGT61
		uint64*brrdat;			// Modular data in here
	#endif
		int*arr_scratch;
		int n;					// Chunksize
		struct complex*rt0;		// Roots table 1
		struct complex*rt1;		// Roots table 2
	#ifdef USE_FGT61
		uint128 *mt0;
		uint128 *mt1;
	#endif
		int*index;				// Bit-reversal index array
		int*block_index;		// 2-data-blocks-per-thread indexing needed by the complex/real FFT wrapper step used by mers-mod
		int nradices_prim;
		int*radix_prim;
		int*ws_i;
		int*ws_j1;
		int*ws_j2;
		int*ws_j2_start;
		int*ws_k;
		int*ws_m;
		int*ws_blocklen;
		int*ws_blocklen_sum;
	};
  #endif

#else
	#define CTIME	/* In single-thread mode, prefer cycle-based time because of its finer granularity */
#endif

/* Extra vars for storing time spent in wrapper/dyadic-square and carry steps: */
#ifdef CTIME
	double dt_fwd, dt_inv, dt_sqr, dt_supp;
	clock_t clock_supp;
#endif

/***************/

#ifdef USE_FGT61
int	mers_mod_square(double a[], uint64 b[], int arr_scratch[], int n, int ilo, int ihi, uint64 p, int scrnFlag, double *tdiff)
#else
int	mers_mod_square(double a[],             int arr_scratch[], int n, int ilo, int ihi, uint64 p, int scrnFlag, double *tdiff)
#endif
{
	const char func[] = "mers_mod_square";
/*...Subroutine to perform Mersenne-mod squaring using Crandall and Fagin's discrete weighted transform (DWT)
     on the data in the length-N real vector A.

     Acronym: DIF = Decimation In Frequency, DIT = Decimation In Time.

     One-pass combined fwd-DWT-weighting/fwd-DIF-FFT/pointwise-squaring/inv-DIT-FFT/inv-DWT-weighting routine
     for use with a packed-vector complex transform to halve the runlength of the corresponding real-vector transform.

     At the beginning of each loop execution, data are assumed to have already been forward-weighted
     and the initial-radix pass of the S-pass forward FFT done. The loop then does the following:

     (1) Performs the 2nd through (S-1)st passes of the complex DIF forward FFT;
     (2) Does the final-radix forward FFT pass, the real/complex-wrapper/pointwise squaring step,
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
	struct qfloat qmul, qwt, qt, qn;	/* qfloats used for DWT weights calculation. */
	struct qfloat qtheta, qr, qi, qc, qs;	/* qfloats used for FFT sincos  calculation. */
	double t1,t2;
	/* Compare qfloat versions of precomputed-table data vs. stdlib double result: */
	const double err_threshold = 1e-12;
	long double theta,mt;
	double adiff, max_adiff = 0.0;	/* Use to store the max abs error between real*8 and real*16 computed values */
	 int64 i1,i2;
	uint64 idiff, max_idiff = 0;
	const double mult[2] = {1.0,-1.0};
	static double nh_inv,nq_inv;	// Needed for "which complex quadrant?" computation
	static int nh,nq;			// #rt1 elts in each quadrant
	int qodd;
	double *re_im_ptr;
	static int radix_set_save[10] = {1000,0,0,0,0,0,0,0,0,0};
	static int radix0, nchunks; 	// Store frequently-used RADIX_VEC[0] and number-of-independently-doable work units
#if DBG_THREADS
	int num_chunks[16];		/* Collect stats about how much work done by each thread */
#endif
	static int nradices_prim,nradices_radix0,radix_prim[30];/* RADIX_PRIM stores sequence of complex FFT radices used, in terms of their prime factors.	*/
	static int *index = 0x0, *si = 0x0, *index_ptmp = 0x0, *si_ptmp = 0x0;		/* Bit-reversal index array and array storing S*I mod N values for DWT weights.	*/
	static int *block_index;				/* array storing the RADIX_VEC[0] data-block indices for pass-2 of the FFT.	*/
	/* arrays storing the index values needed for the paired-block wrapper/square scheme: */
	static int *ws_i,*ws_j1,*ws_j2,*ws_j2_start,*ws_k,*ws_m,*ws_blocklen,*ws_blocklen_sum;
	int bimodn,simodn;					/* Mnemonic: BIMODN stands for "B times I mod N", nothing to do with bimodal.	*/
	int i,ii,ierr,iter,j,jhi,k,l,m,mm,k2,m2,l1,l2,l2_start,blocklen,blocklen_sum,outer;
	static uint64 psave=0;
	static uint32 nsave=0, new_runlength=0;
	static uint32 nwt,nwt_bits,bw,sw,bits_small;
	const  double one_half[3] = {1.0, 0.5, 0.25};		/* Needed for small-weights-tables scheme */
	static double base[2],baseinv[2],radix_inv;
	static struct complex *rt0 = 0x0, *rt1 = 0x0, *rt0_ptmp = 0x0, *rt1_ptmp = 0x0;		/* reduced-size roots of unity arrays	*/
#ifdef USE_FGT61
	const uint64 q  = 0x1FFFFFFFFFFFFFFFull;	// 2^61 - 1
	static uint128 *mt0 = 0x0, *mt1 = 0x0, *mt0_ptmp = 0x0, *mt1_ptmp = 0x0;		// Mod-M61 analogs of same
	static uint8 mod_wt_exp[64];	// Pad the required 61 slots out to 64 for alignment reasons
	uint64 rm,im, rtmp,itmp, order;
	static uint32 swmod61,nmod61;
	uint32 simodnmod61;
	const double inv61 = 1.0/61.0;
	const uint64 qhalf = 0x1000000000000000ull;	// (q+1)/2 = 2^60
#endif
	static double *wt0 = 0x0, *wt1 = 0x0, *tmp = 0x0, *wt0_ptmp = 0x0, *wt1_ptmp = 0x0, *tmp_ptmp = 0x0;		/* reduced-size DWT weights arrays	*/
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

	/* These came about as a result of multithreading, but now are needed whether built unthreaded or multithreaded */
	static int init_sse2 = FALSE;
	int saved_init_sse2, thr_id = -1;	// No multithread support yet.

#ifdef MULTITHREAD

	#ifdef USE_PTHREAD

		int isum;
		static int *thr_ret = 0x0;
		static pthread_t *thread = 0x0;
		static pthread_attr_t attr;
		static struct mers_thread_data_t *tdat = 0x0;

	  #ifdef USE_THREADPOOL	// Threadpool-based dispatch:

		static int main_work_units = 0, pool_work_units = 0;
		static struct threadpool *tpool = 0x0;
		static int task_is_blocking = TRUE;
		static thread_control_t thread_control = {0,0,0};
		// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
		static task_control_t   task_control = {NULL, (void*)mers_process_chunk, NULL, 0x0};

	  #endif

	#endif

#endif

	radix0 = RADIX_VEC[0];
	nchunks = radix0>>1;
	ASSERT(HERE, TRANSFORM_TYPE == REAL_WRAPPER, "mers_mod_square: Incorrect TRANSFORM_TYPE!");

/*...initialize things upon first entry */

	/*...If a new exponent, runlength or radix set, deallocate any already-allocated
	allocatable arrays and set first_entry to true:	*/

	if(n != nsave) new_runlength = TRUE;
	if(p != psave || new_runlength) first_entry=TRUE;

	for(i = 0; i < 10; i++)
	{
		if(RADIX_VEC[i] != radix_set_save[i])
		{
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

		for(i = 0; i < NRADICES; i++)
		{
			if(RADIX_VEC[i] == 0)
			{
				sprintf(cbuf, "%s: RADIX_VEC[i = %d] zero, for i < [NRADICES = %d]!",func,i,NRADICES);
				ASSERT(HERE, 0, cbuf);
			}
			radix_set_save[i] = RADIX_VEC[i];
		}
		for(i = NRADICES; i < 10; i++)
		{
			if(RADIX_VEC[i] != 0)
			{
				sprintf(cbuf, "%s: RADIX_VEC[i = %d] nonzero, for i >= [NRADICES = %d]!",func,i,NRADICES);
				ASSERT(HERE, 0, cbuf);
			}
			radix_set_save[i] = 0;
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
			if((1 << DAT_BITS) < 2*RADIX_VEC[NRADICES-1])
			{
				sprintf(cbuf  ,"FATAL: final FFT radix may not exceed = %u!\n", (1 << (DAT_BITS-1)));
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

		/* We do the final DIF FFT radix within the wrapper_square routine, so store
		that block of sincos data there, where they can be merged with the wrapper sincos data:
		*/
		for(i=1; i<NRADICES-1; i++)
		{
			k =k+mm;
			mm=mm*RADIX_VEC[i];
		}

		if(mm*RADIX_VEC[NRADICES-1] != N2)
		{
			sprintf(cbuf  ,"FATAL: product of radices not equal to complex vector length\n");
			fprintf(stderr,"%s", cbuf);
			ASSERT(HERE, 0,cbuf);
		}

/*		index = (int *)calloc(k,sizeof(int));	*/
		index_ptmp = ALLOC_INT(index_ptmp, k);
		if(!index_ptmp)
		{
			sprintf(cbuf  ,"FATAL: unable to allocate array INDEX in %s.\n",func);
			fprintf(stderr,"%s", cbuf);
			ASSERT(HERE, 0,cbuf);
		}
		index = ALIGN_INT(index_ptmp);

		/*...Forward (DIF) FFT sincos data are in bit-reversed order. We define a separate last-pass twiddles
		array within the routine wrapper_square, since that allows us to merge those nicely with the wrapper sincos data.	*/

		k =0;
		l =0;
		mm=1;

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
		case 40:
			nradices_radix0 = 4;
			radix_prim[l++] = 5; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 44 :
			nradices_radix0 = 3;
			radix_prim[l++] =11; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 48 :
			nradices_radix0 = 5;
			radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 52 :
			nradices_radix0 = 3;
			radix_prim[l++] =13; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
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
		case 80 :
			nradices_radix0 = 5;
			radix_prim[l++] = 5; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 88 :
			nradices_radix0 = 4;
			radix_prim[l++] =11; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 96 :
			nradices_radix0 = 6;
			radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 104:
			nradices_radix0 = 4;
			radix_prim[l++] =13; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
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
		case 144:
			nradices_radix0 = 6;
			radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 160:
			nradices_radix0 = 6;
			radix_prim[l++] = 5; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 176:
			nradices_radix0 = 5;
			radix_prim[l++] =11; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 192:
			nradices_radix0 = 7;
			radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 208:
			nradices_radix0 = 5;
			radix_prim[l++] =13; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 224:
			nradices_radix0 = 6;
			radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 240:
			nradices_radix0 = 6;
			radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 256 :
			nradices_radix0 = 8;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 288:
			nradices_radix0 = 7;
			radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 320:
			nradices_radix0 = 7;
			radix_prim[l++] = 5; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 512 :
			nradices_radix0 = 9;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 768:
			nradices_radix0 = 9;
			radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
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
			sprintf(cbuf  ,"FATAL: radix %d not available. Halting...\n",radix0);
			fprintf(stderr,"%s", cbuf);
			ASSERT(HERE, 0,cbuf);
		}

		for(i=1; i < NRADICES; i++)
		{
			/*...Allocate and initialize an index array containing MM indices...	*/

			if(i<(NRADICES-1))
			{
				mm=mm*RADIX_VEC[i-1];	/* MM = product of all the preceding radices	*/

				for(m=0; m < mm; m++)
				{
					index[k+m]=m;
				}

				/*...then bit-reverse INDEX with respect to the accumulated radices.
				The order of radices sent to bit_reverse_int is the reverse of that in which these radices are processed
				in the forward (decimation in frequency) FFT. This is moot for a power-of-2 FFT (or any FFT whose length
				is a prime power), but necessary for general vector lengths which are a product of 2 or more distinct primes.

				If the current (Ith) radix is composite with distinct prime factors (e.g. 15 = 3*5), we must specify these
				factors here in the opposite order from that which is used in the actual FFT-pass routine. For example,
				if the radix-15 pass implementation does 5 radix-3 DFTs, followed by 3 radix-5 DFTs, then we send (3,5)
				as the corresponding reverse-ordered prime radices to the bit-reversal routine, not (5,3).	*/

				bit_reverse_int(&index[k],mm,l,&radix_prim[l-1],-1,(int *)arr_scratch);

				k += mm;
			}

			/*...All radices beyond the initial-pass one are assumed to be powers of 2 in [8,32]:	*/

			switch(RADIX_VEC[i]) {
		/*
			case 2 :
				radix_prim[l++] = 2; break;
			case 4 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		*/
			case 8 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			case 16 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			case 32 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		/*
			case 64 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			case 128 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2;
		*/
			default :
				sprintf(cbuf  ,"FATAL: intermediate radix %d not available. Halting...\n",RADIX_VEC[i]);
				fprintf(stderr,"%s", cbuf);
				ASSERT(HERE, 0,cbuf);
			}

			/* Final radix must be 16 or 32: */
			if(i == NRADICES-1 && RADIX_VEC[i] < 16)
			{
				sprintf(cbuf  ,"FATAL: final radix %d not available. Halting...\n",RADIX_VEC[i]);
				fprintf(stderr,"%s", cbuf);
				ASSERT(HERE, 0,cbuf);
			}
		}
		nradices_prim = l;	for( ; l < 30; l++) { radix_prim[l] = 0; }	// Zero any higher elements which may have been previously set due
								// to use of a smoother n. (Not necessary since use nradices_prim to control array access, but nice to do.
		bw = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw = n - bw;	/* Number of smallwords.	*/
	#ifdef USE_FGT61
		swmod61 = sw % 61;	// Needed for the modular weights (see comments below).
		nmod61  = n  % 61;
		// Calculate the shift counts which substitute for the modular weight factors.
		/*
		Modular analog of WT(J) satisfies WM(J)^N == [2^(s*j mod N)] mod q = 2^[(s*j mod N) mod 61], since q = 2^61 - 1.
		Assume WM(J) = 2^A. Then WM(J)^N = 2^(A*N), so A satisfies the congruence A*N == (s*j mod N) mod 61.
		Example:
			N=512, s*j mod N = 137. (s*j mod N) mod 61 = 15, so A*512 == 15 mod 61, hence A = 54 and WM(J) = 2^54.

		We can afford a brute-force search for the A's since A in [0,60]. These are then stored in a length-61 table,
		and the appropriate element accessed using the current value of (s*j mod N) mod 61, and sent (along with the
		digit to be multiplied by the weight factor) as a shift count to MUL_POW2_MODQ. Since WM is a power of 2, we need
		only store the shift count A, not the modular weight factor WM itself. Modular inverse weights are also easy: since
		WM*WMINV == 1 mod 2^61 - 1, if WM = 2^A, WMINV = 2^(61-A), i.e. to multiply by the inverse of 2^A, simply use 61-A
		as the shift count to be passed to MUL_POW2_MODQ. Sweet indeed!
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
	#endif

		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)radix0));

		bits_small = p/n;			/* number of bits in a smallword.	*/
		base   [0] = (double)(1 << bits_small);	base   [1] = (double)(2*base[0]);
		baseinv[0] = (double)(1.0/base[0]    );	baseinv[1] = (double)(1.0/base[1]);	/* don't need extended precision for this since both bases are powers of 2.	*/

		/*...stuff for the reduced-length DWT weights arrays is here:	*/

		/* No need for a fancy NINT here: */
nwt_bits = (uint32)(log(sqrt(1.0*n))/log(2.0) + 0.5) - 2;	// Jan 2014: -1; reduces nwt to allow more threads to be used at a given N
		nwt = 1 << nwt_bits;	/* To save on storage, we calculate the first NWT weights directly and then re-use
								them N/NWT times, each time multiplying the basic weights by a single scalar multiplier (and times 0.5
								or 1.0). Thus, the total number of weights data is NWT + (N/NWT). To minimize this, we find the positive
								minimum of the function f(x) = x + N/x, which occurs when f' = 1 - N/(x^2) = 0, or x = NWT = sqrt(N),
								which gives the total number of weights data as 2*sqrt(N). However, to make the algorithm efficient,
								we want NWT a power of 2, so we take NWT = 2^[nint(log2(sqrt(N)))]. In the worst case this needs two
								arrays, one with sqrt(2*N) elements, the other with sqrt(2*N)/2, for a total of 3*sqrt(2)*sqrt(N)/2
								elements, only a few percent more than the minimum possible number.
								The NWT basic weights are stored in the WT0 array; the N/NWT scalar multipliers are in WT1.	*/

		if(n%nwt)
		{
			sprintf(cbuf  ,"FATAL: NWT does not divide N!\n");
			fprintf(stderr,"%s", cbuf);
			ASSERT(HERE, 0,cbuf);
		}

		/*...The roots arrays need only be half the dimension of the weights arrays (since we need n/2 complex roots
		vs. n real weights), but have the same total storage since each entry is complex:	*/
		/*
		wt0 = (double *)calloc(nwt+1         ,sizeof(double));
		wt1 = (double *)calloc(n/nwt+radix0  ,sizeof(double));
		tmp = (double *)calloc(n/nwt+1       ,sizeof(double));
		si  = (   int *)calloc(nwt+1         ,sizeof(   int));
		*/
		wt0_ptmp = ALLOC_DOUBLE(wt0_ptmp, nwt+1         );	if(!wt0_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array WT0 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }; wt0 = ALIGN_DOUBLE(wt0_ptmp);
		wt1_ptmp = ALLOC_DOUBLE(wt1_ptmp, n/nwt+radix0  );	if(!wt1_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array WT1 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }; wt1 = ALIGN_DOUBLE(wt1_ptmp);
		tmp_ptmp = ALLOC_DOUBLE(tmp_ptmp, n/nwt+1       );	if(!tmp_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array TMP in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }; tmp = ALIGN_DOUBLE(tmp_ptmp);
		si_ptmp  = ALLOC_INT   ( si_ptmp, nwt+1         );	if(!si_ptmp ){ sprintf(cbuf,"FATAL: unable to allocate array SI  in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }; si  = ALIGN_INT   (si_ptmp );

		/******************************************************************/
		/* Crandall/Fagin weighting factors and number of bits per digit. */
		/******************************************************************/

		/* The TMP array gets rearranged
		in a cache-friendly way to obtain the WT1 array, but keep the TMP array around to simplify the weighting
		steps at the beginning and end of each iteration cycle, where speed is not crucial.	*/

		/* For qfloat implementation, speed things by defining a constant multiplier qmul = 2^(s/N), repeatedly
		multiplying the previous weight factor by that, and dividing by 2 whenever simodn >= N.	*/

		qt   = i64_to_q((int64)sw);
		qn   = i64_to_q((int64) n);
		qt   = qfdiv(qt, qn);		/* s/N...	 */
		qt   = qfmul(qt, QLN2);		/* ...(s/N)*ln(2)...	*/
		qmul = qfexp(qt);			/* ...and get 2^(s/N) via exp[(s/N)*ln(2)]. */
		qwt  = QONE;				/* init weights multiplier chain. */

		t1 = qfdbl(qmul);
		/* Compare qfloat versions of precomputed-table data vs. stdlib double result: */
		t2 = pow(2.0, 1.0*sw/n);
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
			sprintf(cbuf,"INFO: QWT1= %20.15f, DWT = %20.15f DIFFER BY %20.0f\n", t1, t2, (double)idiff);
			fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}

		simodn=0;
		for(i=0; i<nwt; i++)
		{
			t1 = qfdbl(qwt);
			/* Compare qfloat versions of precomputed-table data vs. stdlib double result: */
			t2 = pow(2.0, 1.0*simodn/n);
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
				sprintf(cbuf,"INFO: I = %8d: QWT = %20.15f, DWT = %20.15f DIFFER BY %20.0f\n", i, t1, t2, (double)idiff);
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}

			wt0  [i] = t1;	/* Ith DWT weight factor = 2^[(s*i mod N)/N], where the exponent is done using floating divide.	*/
			si   [i] = simodn;
			/*fprintf(stderr,"I = %d; WT0 = %20.10f; SI = %d\n",i,wt0[i],si[i]);	*/
			simodn = simodn + sw;

			qwt= qfmul(qwt, qmul);

			if(simodn >= n)
			{
				simodn = simodn - n;
				qwt= qfmul_pow2(qwt, -1);
			}
		}

		wt0[nwt] = 2*wt0[0];	/* This is so the case L = 0 comes out right.	*/
		si [nwt] = n;		/* When L = 0, want to ensure (B*J mod N) - SI(NWT) < 0, so set SI(NWT) := N	*/

		k=(int)(1.0*sw*nwt - 1.0*((int)(1.0*sw*nwt/n))*n);	/* The Jth scalar weights multiplier is 2^[(J*S*NWT mod N)/N].	*/
								/* Use real*8 to do the MOD, since S*NWT may be > 31 bits.	*/
		qt   = i64_to_q((int64) k);
		qn   = i64_to_q((int64) n);
		qt   = qfdiv(qt, qn);		/* k/N...	 */
		qt   = qfmul(qt, QLN2);	/* ...(k/N)*ln(2)...	*/
		qmul = qfexp(qt);		/* ...and get 2^(k/N) via exp[(k/N)*ln(2)]. */
		qwt  = QONE;	/* init weights multiplier chain. */

		t1 = qfdbl(qmul);
		/* Compare qfloat versions of precomputed-table data vs. stdlib double result: */
		t2 = pow(2.0, 1.0*k/n);
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
			sprintf(cbuf,"INFO: QWT2= %20.15f, DWT = %20.15f DIFFER BY %20.0f\n", t1, t2, (double)idiff);
			fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}

		j=0;	/* Store I*K mod NN here. We don't directly calculate I*K, since that can overflow a 32-bit integer at large runlengths.	*/
		for(i=0; i<n/nwt; i++)
		{
			t1 = qfdbl(qwt);
			/* Compare qfloat versions of precomputed-table data vs. stdlib double result: */
			t2 = pow(2.0, 1.0*j/n);
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
				sprintf(cbuf,"INFO: J = %8d: QWT = %20.15f, DWT = %20.15f DIFFER BY %20.0f\n", i, t1, t2, (double)idiff);
				fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}
			tmp[i] = t1;
		/*fprintf(stderr,"I = %d; TMP = %20.10f\n",i,tmp[i]);	*/
			j=j+k;

			qwt= qfmul(qwt, qmul);

			if(j >= n)
			{
				j = j - n;
				qwt= qfmul_pow2(qwt, -1);
			}
		}
		tmp[n/nwt] = 2*tmp[0];	/* This is so the case L = 0 comes out right.	*/

		/*...In the actual radix*_ditN_cy_dif1 carry propagation routine, the elements of the second weights table
		as constructed above are accessed in strides of length n/(radix0*nwt), so it makes sense to prearrange
		them so as to replace these long strides with unit strides, and thus to be accessing contiguous data instead:
		*/
	//	printf("%s: Grouping WTS1 elements into contiguous blocks of length %d\n",func,n/(nwt*radix0));
		for(i = 0; i < n/(nwt*radix0); i++)
		{
			for(j = i*radix0, k = i; j < (i+1)*radix0; j++, k += n/(nwt*radix0))
			{
				wt1[j] = tmp[k];	/* Gather (radix0) stride-[n/(radix0*nwt)]-separated data into a contiguous block.		*/
			}
		}
		for(j = n/nwt, k = n/(nwt*radix0); j < (n/nwt+radix0); j++, k += n/(nwt*radix0))
		{
			wt1[j] = tmp[k];	/* This is so the case L = 0 comes out right.	*/
		}

		/**********************************************/
		/* Roots of unity table pairs needed for FFT: */
		/**********************************************/

		/*...The roots arrays need only be half the dimension of the weights arrays (since we need n/2 complex roots
		vs. n real weights), which for the 2-table scheme is reflected in the halved size of table 2:
		*/
	#if 0	// default:
		NRT = nwt;	NRT_BITS = nwt_bits;
		NRTM1 = NRT - 1;
	#else
		i = 1;	// Causes rt0-table to be 2^i times larger than default, rt1-table to be 2^i times smaller
		NRT = nwt<<i;	NRT_BITS = nwt_bits + i;
		NRTM1 = NRT - 1;
	#endif
		/*...The rt0 array stores the (0:NRT-1)th powers of the [N2]th root of unity
		(i.e. will be accessed using the lower lg(NRT) bits of the integer sincos index):
		*/
	#ifdef USE_FGT61
		mt0_ptmp = ALLOC_UINT128(mt0_ptmp, NRT);
		if(!mt0_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array MT0 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		mt0 = ALIGN_UINT128(mt0_ptmp);

		order = N2;
		prim_root_q(order, &rm,&im);	// primitive root of order n/2
		mt0[i].d0 = 1ull;	mt0[i].d1 = 0ull;
		for(i = 1; i < NRT; i++)
		{
			cmul_modq(mt0[i-1].d0,mt0[i-1].d1, rm,im, &rtmp,&itmp);
			mt0[i].d0 = qreduce_full(rtmp);	mt0[i].d1 = qreduce_full(itmp);
		}
	#endif
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

		for(i = 0; i < NRT; i++)
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
	//	printf("%s: Complex-roots arrays have %u, %u elements.\n",func,NRT,n/(2*NRT));
		/*...The rt1 array stores the (0:(n/2)/NRT-1)th powers of the [(n/2)/NRT]th root of unity
		(and will be accessed using the upper bits, <lg(NRT):31>, of the integer sincos index):
		*/
	#ifdef USE_FGT61
		order = N2/NRT;
		mt1_ptmp = ALLOC_UINT128(mt1_ptmp, order);
		if(!mt1_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array MT1 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		mt1 = ALIGN_UINT128(mt1_ptmp);

		prim_root_q(order, &rm,&im);	// primitive root of order n/2
		mt1[i].d0 = 1ull;	mt1[i].d1 = 0ull;
		for(i = 1; i < order; i++)
		{
			cmul_modq(mt1[i-1].d0,mt1[i-1].d1, rm,im, &rtmp,&itmp);
			mt1[i].d0 = qreduce_full(rtmp);	mt1[i].d1 = qreduce_full(itmp);
		}
	#endif
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
//if((i & 63) ==0)printf("rt1[%3u] = %20.15f, %20.15f\n",i,rt1[i].re,rt1[i].im);
			qt = qfadd(qt, qtheta);
			/* Get next terms of the recurrence:  qcnew = qcold*qr - qsold*qi,  qsnew = qcold*qi + qsold*qr. */
			/*qn = qfmul(qc, qr); qt = qfmul(qs, qi); qmul = qfsub(qn, qt);	* Store qcnew in qmul for now. */
			/*qn = qfmul(qc, qi); qt = qfmul(qs, qr); qs   = qfadd(qn, qt); qc = qmul;	*/
		}

/*
Oct 2014:
Exploitable symmetries which can be used to cut size of rt1, based on the quadrant. Let n := nrt/4,
let q: = 0-3, corr. to int(theta/(pi/2), gamma := theta - q*(pi/2), and iq := i - q*n = "index within quadrant":

I-range	qud		root in terms of iq:		Alternative (#else below)
-------	---	 ----------------------------	-----------------------
[ 0, n)	q=0	 rt1[iq].re, rt1[iq].im			 rt1[iq].re, rt1[iq].im
[ n,2n)	q=1	-rt1[nq-iq].re, rt1[nq-iq].im	-rt1[iq].im, rt1[iq].re
[2n,3n)	q=2	-rt1[iq].re,-rt1[iq].im			-rt1[iq].re,-rt1[iq].im
[3n,4n)	q=3	 rt1[nq-iq].re,-rt1[nq-iq].im	 rt1[iq].im,-rt1[iq].re

Cutting the rt1 size by a factor of 4 also means tweaking our relative-table-size formula:
For the original (rt1-size-unreduced) scheme we want the 2 table to be equal-sized in the case of N an even power of 2
So e.g. rt0,1 have B,B elts each, total = 2B.
Now get rt0,1 have B,B/4 elts each, total = (5/4)*B. If instead started with (rt1-unreduced) sizes B/2,2B,
    get rt0,1 with B/2,B/2 elts each, total = B, half the original 2B, and the minimum possible,
 since starting with (rt1-unreduced) sizes B/4,4B ==> B/4,B, total = (5/4)*B, i.e. the total again starts to rise.

We can even go one better by using the fact that all roots can be mapped to data in the first octant: Let n := nrt/8,
let q: = 0-7, corr. to int(theta/(pi/4), gamma := theta - q*(pi/4), and iq := i - q*n = "index within octant":

I-range	oct		root in terms of iq:
-------	---	 ----------------------------
[ 0, n)	q=0	 rt1[   iq].re, rt1[   iq].im
[ n,2n)	q=1	 rt1[nq-iq].im, rt1[nq-iq].re
[2n,3n)	q=2	-rt1[   iq].im, rt1[   iq].re
[3n,4n)	q=3	-rt1[nq-iq].re, rt1[nq-iq].im
[4n,5n)	q=4	-rt1[   iq].re,-rt1[   iq].im
[5n,6n)	q=5	-rt1[nq-iq].im,-rt1[nq-iq].re
[6n,7n)	q=6	 rt1[   iq].im,-rt1[   iq].re
[7n,8n)	q=7	 rt1[nq-iq].re,-rt1[nq-iq].im

Now if start with (rt1-unreduced) sizes B/2,2B,
    get rt0,1 with B/2,B/4 elts each, total = (3/4)*B, vs the original 2B, and the minimum possible.
    (Get same total if start with (rt1-unreduced) sizes B/4,4B ==> B/4,B/2, total = (3/4)*B.)
*/
/*
for(i=0; i < NRT; i++) {
	printf("I = %3d: rt0[i].re,im = %20.15f, %20.15f\n",i,rt0[i].re,rt0[i].im);
}
*/
#if 0
	#define SYMM	2	// "foldness" of the symmetry scheme used: 2 = half-plane, 4 = quadrans, 8 = half-quads.
	#if SYMM == 2
		nh = n/(NRT<<2);	// #rt1 elts in each quadrant
		nh_inv = 1.0/(double)nh;
		printf("half-plane #elts = %d\n",nh);
	#elif SYMM == 4
		nq = n/(NRT<<3);	// #rt1 elts in each quadrant
		nq_inv = 1.0/(double)nq;
		printf("quadrant #elts = %d\n",nq);
	#else
		Value of SYMM unsupported!
	#endif

	printf("rt1 #elts = %d\n",n/(2*NRT));
	for(i=0; i < n/(2*NRT); i++) {
	#if SYMM == 2
		qodd = i >= nh;
		t1 = mult[qodd];	// -1 if root in lower half-plane
		j = i - ((-qodd) & nh);	// i % nh
		double c = rt1[j].re;
		double s = rt1[j].im;
		c *= t1;
		s *= t1;
	if(i > (nh-3) && i < (nh+3))
		printf("I = %3d, J = %3d [h = %1d]: rt1[i].re,im = %17.15f, %17.15f; V2 = %17.15f, %17.15f\n",i,j,qodd,rt1[i].re,rt1[i].im, c,s);
	#elif SYMM == 4
		uint32 iq = (int)((double)i*nq_inv);	// Efficient way of computing i*(NRT<<3)/n
		j = i - iq*nq;		// i % nq
	  #if 0
		qodd = -(iq&1);			// Negate (iq odd?) result to turn into bitmask
		j += qodd & (nq-j-j);	// If quadrant index odd (1 or 3), need nq-j instead of j
		double c = rt1[j].re;
		double s = rt1[j].im;
	  #else
		qodd = iq&1;			// quadrant index odd (1 or 3)?
		re_im_ptr = rt1 + j;	// Cast to double*
		double c = *(re_im_ptr +    qodd );
		double s = *(re_im_ptr + (1-qodd));
	  #endif
		c *= mult[(iq-1) < 2];	// re part negated for quadrants 1,2 (this is why iq needs to be unsigned)
		s *= mult[ iq    > 1];	// im part negated for quadrants 2,3
		printf("I = %3d [q = %1d]: rt1[i].re,im = %17.15f, %17.15f; V2[m1,m2 = %d,%d] = %17.15f, %17.15f\n",i,iq,rt1[i].re,rt1[i].im, (iq-1) < 2,(iq > 1),c,s);
	#endif
	}
	exit(0);
#endif
		/**********************************************/
		/*************** MERSENNE-ONLY: ***************/
		/**********************************************/

		/* Break the remaining portion of the FFT into radix0 blocks, each of which ideally
			 should operate on a dataset which fits entirely into the L2 cache of the host machine. */

		/* 8/23/2004: Need to allocate an extra element here to account for the padding element that gets inserted when radix0 is odd: */
		block_index = (int *)calloc((radix0+1),sizeof(int));
		if(!block_index){ sprintf(cbuf,"FATAL: unable to allocate array BLOCK_INDEX in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		/*
		Examples:

			radix0 = 3: want 2 sets of l1,2 pairs = {0,-},{1,2}
			radix0 = 4: want 2 sets of l1,2 pairs = {0,1},{2,3}
			radix0 = 5: want 2 sets of l1,2 pairs = {0,-},{1,4,2,3}
			radix0 = 6: want 2 sets of l1,2 pairs = {0,1},{2,5,3,4}
			radix0 = 7: want 2 sets of l1,2 pairs = {0,-},{1,6,2,5,3,4}
			radix0 = 8: want 3 sets of l1,2 pairs = {0,1},{2,3},{4,7,5,6}
			radix0 = 9: want 3 sets of l1,2 pairs = {0,-},{1,2},{3,8,4,7,5,6}
			radix0 =10: want 2 sets of l1,2 pairs = {0,1},{2,9,3,8,4,7,5,6}
			radix0 =11: want 2 sets of l1,2 pairs = {0,-},{1,10,2,9,3,8,4,7,5,6}
			radix0 =12: want 3 sets of l1,2 pairs = {0,1},{2,3},{4,11,5,10,6,9,7,8}
			radix0 =13: want 2 sets of l1,2 pairs = {0,-},{1,12,2,11,3,10,4,9,5,8,6,7}
			radix0 =14: want 2 sets of l1,2 pairs = {0,1},{2,13,3,12,4,11,5,10,6,9,7,8}				blocklen = 2,12		l2_start = 1,13		increments by 12
			radix0 =15: want 3 sets of l1,2 pairs = {0,-},{1,2},{3,14,4,13,5,12,6,11,7,10,8,9}		blocklen = 1,2,12	l2_start = 0,2,14	increments by 2,12
			radix0 =16: want 4 sets of l1,2 pairs = {0,1},{2,3},{4,7,5,6},{8,15,9,14,10,13,11,12}	blocklen = 2,2,4,8	l2_start = 1,3,7,15	increments by 2,4,8

		Rules:
			1) For even radix0 use [nradices_radix0] sets of blocks; for odd radices use [nradices_radix0+1] sets,
				the first of which consists of a single block.
			2) Outer loop over prime factors of radix0, in reverse order of which these subradices appear in forward FFT:
				- blocklen (number of l-values of) first set = 2 - radix0%2
				- blocklen of sets 2...nradices_radix0 = (number of blocks done in previous sets) * (current subradix - 1) .
			3) Throughout block processing, l1 increases monotonically from 0. Within each block, for each value of l1, the
				corresponding l2_start starts at 1 - radix0%2and increases by blocklen_newblock for each new block.

		In a multithreaded implementation, each thread will be responsible for 2 of the above blocks of data (unless radix0
		odd, in which case one of the threads only gets one block to chew on.) E.g. if NTHREADS = 4 and radix0 = 14, the
		chunk-processing loop occurs in 2 parallel passes, with each of the 4 threads processing the following blockpairs:

			Pass 1:
			Thread 0: ( 0, 1)
			Thread 1: ( 2,13)
			Thread 2: ( 3,12)
			Thread 3: ( 4,11)

			Pass 2:
			Thread 0: ( 5,10)
			Thread 1: ( 6, 9)
			Thread 2: ( 7, 8)
			Thread 3:   ---		(i.e. thread is idle)

		 For radix0 = 15, things look as follows:

			Pass 1:
			Thread 0: ( 0, -)	(i.e. thread is only 50% utilized)
			Thread 1: ( 1, 2)
			Thread 2: ( 3,14)
			Thread 3: ( 4,13)

			Pass 2:
			Thread 0: ( 5,12)
			Thread 1: ( 6,11)
			Thread 2: ( 7,10)
			Thread 3: ( 8, 9)

		Thus, if (2*NTHREADS) does not divide radix0, there will always be one or more under-or-unutilized threads.
		*/
		blocklen = 2 - (radix0 & 1);	/* blocklen = 2 for even radix0, 1 for odd. */
		blocklen_sum=0;
		l1=0;
		l2=blocklen - 1; l2_start=l2;		/* l2_start = 1 for even radix0, 0 for odd. */

		/* Init the block_index array: */
		ii = 0;	/* Need a linear index to provide access into the block_index array: */

		for(outer = nradices_radix0 - l2; outer >= 0; outer-- )	/* This mimics the increasing-blocksize loop in the wrapper_square routines. */
		{
			/*fprintf(stderr,"Block %d : blocklen = %d  blocklen_sum = %d\n",nradices_radix0 - outer,blocklen,blocklen_sum);*/

			for(m = 0; m <= (blocklen-1)>>1; m++)	/* Since we now process TWO blocks per loop execution, only execute the loop half as many times as before. */
			{
				/* Now execute body of loop once with l = l1, once with l = l2 (only do l1 if blocklen == 1).
				Since blocklen is even except when radix0 odd and first pass, loop from 0 to 1 - blocklen%2.
				*/
				l = l1;
				/*for(j = 0; j <= 1 - (blocklen & 1); j++)*/
				/* Do two loop executions. For the case where blocklen is odd, insert a dummy padding element
				(indicated by its having a negative value) into the second (j = 1) slot of the block_index array:
				*/
				for(j = 0; j < 2; j++)
				{
					if(!(l >= 0 && l < radix0)) { sprintf(cbuf,"ERROR 10 in %s.c\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }

					if((blocklen & 1) && j == 1)
					{
						block_index[ii] = -1;
						/*fprintf(stderr,"%3d ---\n",ii);*/
					}
					else
					{
						block_index[ii] = l;
						/*fprintf(stderr,"%3d %3d\n",ii,l);*/
					}

					ii++;	/* every time we execute this innermost loop (which corresponds to
						 one block of FFT data being processed), increment the linear array index */
					l += l2-l1;
				}	/* end j-loop */
				l1++;
				l2--;
			}

			if(outer == 0)break;

			blocklen_sum += blocklen;
			blocklen = (radix_prim[outer-1] - 1)*blocklen_sum;

			/*...Next j2_start is previous one plus the length of the current block */

			l1 = l2_start + 1;
			l2_start += blocklen;
			l2 = l2_start;			/* Reset j2 for start of the next block. */
		}		/* End of Main loop */

		/* arrays storing the index values needed for the parallel-block wrapper/square scheme: */
		if( !(ws_i            = (int *)calloc(radix0,sizeof(int))) ) { sprintf(cbuf,"FATAL: unable to allocate array WS_I            in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		if( !(ws_j1           = (int *)calloc(radix0,sizeof(int))) ) { sprintf(cbuf,"FATAL: unable to allocate array WS_J1           in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		if( !(ws_j2           = (int *)calloc(radix0,sizeof(int))) ) { sprintf(cbuf,"FATAL: unable to allocate array WS_J2           in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		if( !(ws_j2_start     = (int *)calloc(radix0,sizeof(int))) ) { sprintf(cbuf,"FATAL: unable to allocate array WS_J2_START     in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		if( !(ws_k            = (int *)calloc(radix0,sizeof(int))) ) { sprintf(cbuf,"FATAL: unable to allocate array WS_K            in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		if( !(ws_m            = (int *)calloc(radix0,sizeof(int))) ) { sprintf(cbuf,"FATAL: unable to allocate array WS_M            in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		if( !(ws_blocklen     = (int *)calloc(radix0,sizeof(int))) ) { sprintf(cbuf,"FATAL: unable to allocate array WS_BLOCKLEN     in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		if( !(ws_blocklen_sum = (int *)calloc(radix0,sizeof(int))) ) { sprintf(cbuf,"FATAL: unable to allocate array WS_BLOCKLEN_SUM in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }

		for(ii = 0; ii < radix0; ii += 2)
		{
			/*...Final DIF pass, wrapper/squaring and initial DIT pass are all done in-place.
				 This combines data from both the l1 and l2-block, except in the case ii = 0
				 for even radix0, for which the l1 = 0 and l2 = 1 blocks are processed separately within
				 wrapper_square, i.e. we must call this routine a second time to process data in the l2-block.
			*/
			if(ii == 0 && !(radix0 & 1))
				jhi = 2;
			else
				jhi = 1;

			for(j = 0; j < jhi; j++)
			{
				l = ii + j;	/* This will actually leave some 'holes' in the init arrays, but we
						don't care, as these are small and we only use the init'ed elements. */

				switch(RADIX_VEC[NRADICES-1])
				{
					case 16 :
						radix16_wrapper_ini(n, radix0, l, nradices_prim, radix_prim, ws_i, ws_j1, ws_j2, ws_j2_start, ws_k, ws_m, ws_blocklen, ws_blocklen_sum);
					//	radix16_wrapper_square(a,arr_scratch,n,radix0,rt0,rt1,nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l],TRUE);
						break;
					case 32 :
						radix32_wrapper_ini(n, radix0, l, nradices_prim, radix_prim, ws_i, ws_j1, ws_j2, ws_j2_start, ws_k, ws_m, ws_blocklen, ws_blocklen_sum);
					//	radix32_wrapper_square(a,arr_scratch,n,radix0,rt0,rt1,nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l],TRUE);
						break;
					/*
					case 64 :
						radix64_wrapper_ini(n, radix0, l, nradices_prim, radix_prim, ws_i, ws_j1, ws_j2, ws_j2_start, ws_k, ws_m, ws_blocklen, ws_blocklen_sum);	break;
					//	radix64_wrapper_square(a,arr_scratch,n,radix0,rt0,rt1,nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l],TRUE);
						break;
					*/
					default :
					sprintf(cbuf,"FATAL: radix %d not available for wrapper_square. Halting...\n",RADIX_VEC[NRADICES-1]);
					fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
				}
			}
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

		thr_ret = (int *)calloc(nchunks, sizeof(int));
		thread  = (pthread_t *)calloc(nchunks, sizeof(pthread_t));
		tdat    = (struct mers_thread_data_t *)calloc(nchunks, sizeof(struct mers_thread_data_t));

		/* Initialize and set thread detached attribute */
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		/* Populate the thread-specific data structs: */
		for(i = 0; i < nchunks; ++i)
		{
			tdat[i].tid = i;
			tdat[i].retval = &thr_ret[i];
			tdat[i].arrdat = a;			// Main data array
		#ifdef USE_FGT61
			tdat[i].brrdat = b;			// Modular data in here
		#endif
			tdat[i].arr_scratch = arr_scratch;
			tdat[i].n = n;					// Chunksize
			tdat[i].rt0 = rt0;	// Roots table 1
			tdat[i].rt1 = rt1;	// Roots table 2
		#ifdef USE_FGT61
			tdat[i].mt0 = mt0;
			tdat[i].mt1 = mt1;
		#endif
			tdat[i].index = index;				// Bit-reversal index array
			tdat[i].block_index = block_index;	// 2-data-blocks-per-thread indexing needed by the complex/real FFT wrapper step used by mers-mod
			tdat[i].nradices_prim = nradices_prim;
			tdat[i].radix_prim = radix_prim;
			tdat[i].ws_i = ws_i;
			tdat[i].ws_j1 = ws_j1;
			tdat[i].ws_j2 = ws_j2;
			tdat[i].ws_j2_start = ws_j2_start;
			tdat[i].ws_k = ws_k;
			tdat[i].ws_m = ws_m;
			tdat[i].ws_blocklen = ws_blocklen;
			tdat[i].ws_blocklen_sum = ws_blocklen_sum;
		}
	  #endif

	  #ifdef USE_THREADPOOL	// Threadpool-based dispatch:

		// MAX_THREADS is the max. no. of threads we expect to be able to make use of, at 1 thread per core.
		ASSERT(HERE, MAX_THREADS == get_num_cores(), "MAX_THREADS not set or incorrectly set!");

		if(nchunks % NTHREADS != 0) fprintf(stderr,"%s: radix0/2 not exactly divisible by NTHREADS - This will hurt performance.\n",func);

		// MacOS does weird things with threading (e.g. Idle" main thread burning 100% of 1 CPU)
		// so on that platform try to be clever and interleave main-thread and threadpool-work processing
		#if 0//def OS_TYPE_MACOSX

			// Under OS X we want one core dispatching serial non-pool threads while the others crunch pool threads,
			// in a manner which balances the load optimally:
			if(NTHREADS > 1) {
				main_work_units = nchunks/MAX_THREADS;
				pool_work_units = nchunks - main_work_units;
				ASSERT(HERE, 0x0 != (tpool = threadpool_init(MAX_THREADS-1, MAX_THREADS, pool_work_units, &thread_control)), "threadpool_init failed!");
				printf("%s: Init threadpool of %d threads\n",func,MAX_THREADS-1);
			} else {
				printf("%s: NTHREADS = 1: Using main execution thread, no threadpool needed.\n",func);
			}

		#else

			main_work_units = 0;
			pool_work_units = nchunks;
			ASSERT(HERE, 0x0 != (tpool = threadpool_init(NTHREADS, MAX_THREADS, pool_work_units, &thread_control)), "threadpool_init failed!");
			printf("%s: Init threadpool of %d threads\n",func,NTHREADS);

		#endif

	  #endif

	#endif
	}

	/* 	This set of init-mode calls needs to go below above init-block because several
	of the inits need the primitive-radix data to have been inited.
	*/
	// Apr 2014: Thanks to Stephen Searle [SMUS] for the init_sse2-related bugfix:
	saved_init_sse2 = init_sse2;	// SMJS init_sse2 gets changed in first if, so need to store its current value for second if statement
	if(new_runlength && init_sse2) {	// Pvsly inited SSE2 local storage, but now have new runlength
		init_sse2 = nchunks;	// Use *value* of init_sse2 to store #threads
		thr_id = -1;
		/* The dyadic-square routines need a few more params passed in init-mode than do the standalone FFT-pass routines: */
	#ifdef USE_FGT61
		radix16_wrapper_square(0x0,0x0, arr_scratch, n, radix0, 0x0,0x0, 0x0,0x0, nradices_prim, radix_prim, 0,0,0,0,0,0,0,0, init_sse2, thr_id);
	#else
		radix16_wrapper_square(0x0,     arr_scratch, n, radix0, 0x0,0x0,          nradices_prim, radix_prim, 0,0,0,0,0,0,0,0, init_sse2, thr_id);
	#endif
		radix32_wrapper_square(0x0, arr_scratch, n, radix0, 0x0, 0x0, nradices_prim, radix_prim, 0,0,0,0,0,0,0,0, init_sse2, thr_id);
	}
	if(init_sse2 == FALSE || (saved_init_sse2 < nchunks)) {		// New run, or need to up #threads in local-store inits
	//	init_sse2 = TRUE;
		init_sse2 = nchunks;	// Use *value* of init_sse2 to store #threads
		thr_id = -1;
		radix8_dif_pass (0x0, 0, 0x0, 0x0, 0x0, 0, 0, init_sse2, thr_id);
	#ifdef USE_FGT61
		radix16_dif_pass(0x0,0x0, 0, 0x0,0x0, 0x0,0x0, 0, 0, 0x0, init_sse2, thr_id);
	#else
		radix16_dif_pass(0x0,     0, 0x0,0x0,          0, 0, 0x0, init_sse2, thr_id);
	#endif
		radix32_dif_pass(0x0, 0, 0x0, 0x0, 0, 0, 0x0, init_sse2, thr_id);
		/* The dyadic-square routines need a few more params passed in init-mode than do the standalone FFT-pass routines: */
		// Dec 2017: Add rt0,rt1-pointers to the wrapper_square init calls to support USE_PRECOMPUTED_TWIDDLES build option:
	#ifdef USE_FGT61
		radix16_wrapper_square(0x0,0x0, arr_scratch, n, radix0, rt0,rt1, 0x0,0x0, nradices_prim, radix_prim, 0,0,0,0,0,0,0,0, init_sse2, thr_id);
	#else
		radix16_wrapper_square(0x0,     arr_scratch, n, radix0, rt0,rt1,          nradices_prim, radix_prim, 0,0,0,0,0,0,0,0, init_sse2, thr_id);
	#endif
		radix32_wrapper_square(0x0, arr_scratch, n, radix0, rt0,rt1, nradices_prim, radix_prim, 0,0,0,0,0,0,0,0, init_sse2, thr_id);

		radix8_dit_pass (0x0, 0, 0x0, 0x0, 0, 0, 0x0, init_sse2, thr_id);
	#ifdef USE_FGT61
		radix16_dit_pass(0x0,0x0, 0, 0x0,0x0, 0x0,0x0, 0, 0, 0x0, init_sse2, thr_id);
	#else
		radix16_dit_pass(0x0,     0, 0x0,0x0,          0, 0, 0x0, init_sse2, thr_id);
	#endif
		radix32_dit_pass(0x0, 0, 0x0, 0x0, 0, 0, 0x0, init_sse2, thr_id);
	}
	new_runlength = FALSE;

	/* end of initialization sequence.	*/

/**********************************************************************/

/*...Notes on generation of the (radix-1) complex sincos data needed for each block of (radix) complex array data.

     Begin with two small tables: RT0 contains the first NRT = O(sqrt(N)) Nth complex roots of unity;
                                  RT1 contains all of the (N/NRT)th complex roots of unity.
     NRT is chosen to be a power of 2, to make table lookups fast. Then, to generate any desired Jth power of the primitive
     Nth root of unity (0 <= J < N), one can simply do a complex multiply of elements from each of these small tables:

  		exp(i*2*pi*J/N) = RT0(int[J/NRT]) * RT1(J mod NRT) .	(cmul)

    Since NRT is a power of 2, the integer divide-and-truncate int[J/NRT] needs just a rightward binary shift,
     and the remainder (J mod NRT) needs just a mask, i.e. a bitwise AND(J,NRT-1).  Both of these are fast,
     and the smallness of the tables means that they can completely fit into the machine's data cache, so that
     the array accesses are fast even if the index patterns are not sequential.

     But we can do even better: if we are willing to accept a (very) slight increase in the level of rounding error,
     we can drastically reduce the number of even the relatively cheap roots array accesses, and reduce the average
     number of floating-point operations needed to get each needed complex root, from 6 for the above scheme (4 FMUL, 2 FADD)
     to something closer to 4. We also improve the balance of floating multiplies and adds at the same time. Here's how:

     In each case we allow at most binary products of roots gotten via table multiply (in order to keep RO error in check),
     and we seek combinations that are fast to calculate, i.e. have a relatively small number of floating-point operations,
     and a good balance of FADD and FMUL.

     Summary of sincos generation for the various radices supported by Mlucas:

     RADIX-5: need powers 1-4:
     Use small-tables multiply to get powers 1,3:
     ------
     1x		standard cmul					4 mul, 2 add
     2=3-1	(c1,-s1)*(c3,s3) = (c1*c3+s1*s3, c1*s3-s1*c3)	4 mul, 2 add; store the 4 scalar products...
     3x		standard cmul					4 mul, 2 add
     4=3+1	(c1,+s1)*(c3,s3) = (c1*c3-s1*s3, c1*s3+s1*c3)	0 mul, 2 add ...and then need just 4 adds to get both 2,4
     							totals: 12 mul, 8 add	(5.0 ops per complex datum) + 2 table look-ups

     RADIX-6: need powers 1-5:
     Use small-tables multiply to get powers 1,4:		opcount:
     ------	----------------------------			------------
     1x		standard cmul					4 mul, 2 add + table look-up
     2=1+1	standard csqr					2 mul, 3 add
     3=4-1	get together with 5 using just...		4 mul, 4 add
     4x		standard cmul					4 mul, 2 add + table look-up
     5=4+1						totals: 14 mul, 11 add	(5.0 ops per complex datum) + 2 table look-ups

     RADIX-7: need powers 1-6:
     Use small-tables multiply to get powers 2,3:		opcount:
     ------	----------------------------			------------
     1=3-2	get together with 5 using just...		4 mul, 4 add
     2x		standard cmul					4 mul, 2 add + table look-up
     3x		standard cmul					4 mul, 2 add + table look-up
     4=2+2	standard csqr					2 mul, 3 add
     5=3+2
     6=3+3	standard csqr					2 mul, 3 add
     							totals: 16 mul, 14 add	(5.0 ops per complex datum) + 2 table look-ups

     RADIX-8: need powers 1-7:
     Use small-tables multiply to get powers 1,2,5:		opcount:
     ------	----------------------------			------------
     1x		standard cmul					4 mul, 2 add + table look-up
     2x		standard cmul					4 mul, 2 add + table look-up
     3=5-2	get together with 7 using just...		4 mul, 4 add
     4=5-1	get together with 6 using just...		4 mul, 4 add
     5x		standard cmul					4 mul, 2 add + table look-up
     6=5+1
     7=5+2						totals: 20 mul, 14 add	(4.86 ops per complex datum) + 3 table look-ups

     On the other hand, if slightly less (but still decent) accuracy is needed, we could use one less table look-up:
     1x		standard cmul					4 mul, 2 add + table look-up
     2=1+1	standard csqr					2 mul, 3 add
     3=5-2	get together with 7 using just...		4 mul, 4 add
     4=2+2	standard csqr					2 mul, 3 add
     5x		standard cmul					4 mul, 2 add + table look-up
     6=3+3	standard csqr					2 mul, 3 add
     7=5+2						totals: 18 mul, 17 add	(5.00 ops per complex datum) + 2 table look-ups
     ...but the very slight time savings this yields are usually not worth the decrease in accuracy.

     RADIX-9: need powers 1-8:
     Use small-tables multiply to get powers 1,2,6:		opcount:
     ------	----------------------------			------------
     1x		standard cmul					4 mul, 2 add + table look-up
     2x		standard cmul					4 mul, 2 add + table look-up
     3=1+2	standard cmul					4 mul, 2 add
     4=6-2	get together with 8 using just...		4 mul, 4 add
     5=6-1	get together with 7 using just...		4 mul, 4 add
     6x		standard cmul					4 mul, 2 add + table look-up
     7=6+1
     8=6+2						totals: 24 mul, 16 add	(5.0 ops per complex datum) + 3 table look-ups

     RADIX-10: need powers 1-9:
     Use small-tables multiply to get powers 1,2,7:		opcount:
     ------	----------------------------			------------
     1x		standard cmul					4 mul, 2 add + table look-up
     2x		standard cmul					4 mul, 2 add + table look-up
     3=2+1	standard cmul					4 mul, 2 add
     4=2+2	standard csqr					2 mul, 3 add
     5=7-2	get together with 9 using just...		4 mul, 4 add
     6=7-1	get together with 8 using just...		4 mul, 4 add
     7x		standard cmul					4 mul, 2 add + table look-up
     8=7+1
     9=7+2						totals: 26 mul, 19 add	(5.0 ops per complex datum) + 3 table look-ups

     RADIX-11: need powers 1-10:
     Use small-tables multiply to get powers 1,2,7:		opcount:
     ------	----------------------------			------------
     1x		standard cmul					4 mul, 2 add + table look-up
     2=3-1	get together with 4 using just...		4 mul, 4 add
     3x		standard cmul					4 mul, 2 add + table look-up
     4=3+1
     5=6-1	get together with 7 using just...		4 mul, 4 add
     6x		standard cmul					4 mul, 2 add + table look-up
     7=6+1
     8=9-1	get together with 10 using just...		4 mul, 4 add
     9x		standard cmul					4 mul, 2 add + table look-up
     10=9+1						totals: 28 mul, 20 add	(4.8 ops per complex datum) + 4 table look-ups

     RADIX-12: need powers 1-11:
     Use small-tables multiply to get powers 1,2,3,8:		opcount:
     --------	----------------------------			------------
     1x		standard cmul					4 mul, 2 add + table look-up
     2x		standard cmul					4 mul, 2 add + table look-up
     3x		standard cmul					4 mul, 2 add + table look-up
     4 =2+2	standard csqr					2 mul, 3 add
     5 =8-3	get together with 11 using just...		4 mul, 4 add
     6 =8-2	get together with 10 using just...		4 mul, 4 add
     7 =8-1	get together with 9  using just...		4 mul, 4 add
     8x		standard cmul					4 mul, 2 add + table look-up
     9 =8+1
     10=8+2
     11=8+3						totals: 30 mul, 23 add	(4.82 ops per complex datum) + 4 table look-ups

     RADIX-14: need powers 1-13:
     Use small-tables multiply to get powers 1,2,8,11:
     ---------
     1x		standard cmul					4 mul, 2 add + table look-up
     2x		standard cmul					4 mul, 2 add + table look-up
     3 =5-2	get together with 7 using just...		4 mul, 4 add
     4 =5-1	get together with 6 using just...		4 mul, 4 add
     5x		standard cmul					4 mul, 2 add + table look-up
     6 =5+1
     7 =5+2
     8=4+4	standard csqr					2 mul, 3 add
     9 =11-2	get together with 13 using just...		4 mul, 4 add
     10=11-1	get together with 12 using just...		4 mul, 4 add
     11x	standard cmul					4 mul, 2 add + table look-up
     12=11+1
     13=11+2						totals: 34 mul, 27 add	(4.69 ops per complex datum) + 4 table look-ups

     RADIX-15: need powers 1-14:
     Use small-tables multiply to get powers 1,2,3,7,12:
     ---------
     1x		standard cmul					4 mul, 2 add + table look-up
     2x		standard cmul					4 mul, 2 add + table look-up
     3x		standard cmul					4 mul, 2 add + table look-up
     4 =7-3	get together with 10 using just...		4 mul, 4 add
     5 =7-2	get together with 9 using just...		4 mul, 4 add
     6 =7-1	get together with 8 using just...		4 mul, 4 add
     7x		standard cmul					4 mul, 2 add + table look-up
     8 =7+1
     9 =7+2
     10=7+3
     11=12-1	get together with 13 using just...		4 mul, 4 add
     12x	standard cmul					4 mul, 2 add + table look-up
     13=12+1
     14=7+7	standard csqr					2 mul, 3 add
  							totals: 38 mul, 29 add	(4.79 ops per complex datum) + 5 table look-ups

     RADIX-16: need powers 1-15:
     Use small-tables multiply to get powers 1,2,8,13:
     ---------
     1x		standard cmul					4 mul, 2 add + table look-up
     2x		standard cmul					4 mul, 2 add + table look-up
     3 =2+1	standard cmul					4 mul, 2 add
     4 =2+2	standard csqr					2 mul, 3 add
     5 =13-8	standard cmul					4 mul, 2 add
     6 =8-2	get together with 10 using just...		4 mul, 4 add
     7 =8-1	get together with 9 using just...		4 mul, 4 add
     8x		standard cmul					4 mul, 2 add + table look-up
     9 =8+1
     10=8+2
     11=13-2	get together with 15 using just...		4 mul, 4 add
     12=13-1	get together with 14 using just...		4 mul, 4 add
     13x	standard cmul					4 mul, 2 add + table look-up
     14=13+1
     15=13+2						totals: 42 mul, 31 add	(4.87 ops per complex datum) + 4 table look-ups

     73 ops seems a bit steep, though - perhaps we can save quite a few by improving the symmetry of the recurrence,
     i.e. getting rid of those expensive unpaired cmuls: e.g. 8,13 partially wasted above since don't use 13+8, just 13-8.

     Use small-tables multiply to get powers 1,2,4,8,13:
     -----------
     1x		standard cmul					4 mul, 2 add + table look-up
     2x		standard cmul					4 mul, 2 add + table look-up
     3 =4-1	get together with 5 using just...		4 mul, 4 add
     4x		standard cmul					4 mul, 2 add + table look-up
     5 =4+1
     6 =8-2	get together with 10 using just...		4 mul, 4 add
     7 =8-1	get together with 9 using just...		4 mul, 4 add
     8x		standard cmul					4 mul, 2 add + table look-up
     9 =8+1
     10=8+2
     11=13-2	get together with 15 using just...		4 mul, 4 add
     12=13-1	get together with 14 using just...		4 mul, 4 add
     13x	standard cmul					4 mul, 2 add + table look-up
     14=13+1
     15=13+2						totals: 40 mul, 30 add	(4.67 ops per complex datum) + 5 table look-ups

     Thus, for one extra CMUL (4x) save ourselves 3 ops, and get slightly better accuracy.

     RADIX-25: need powers 1-24:
     Use small-tables multiply to get powers 1,2,3,7,14,21:
     -----------
     1x		standard cmul					4 mul, 2 add + table look-up
     2x		standard cmul					4 mul, 2 add + table look-up
     3x		standard cmul					4 mul, 2 add + table look-up
     4 =7-3	get together with 10 using just...		4 mul, 4 add
     5 =7-2	get together with 9  using just...		4 mul, 4 add
     6 =7-1	get together with 8  using just...		4 mul, 4 add
     7x		standard cmul					4 mul, 2 add + table look-up
     8 =7+1
     9 =7+2
     10=7+3
     11=14-3	get together with 17 using just...		4 mul, 4 add
     12=14-2	get together with 16 using just...		4 mul, 4 add
     13=14-1	get together with 15 using just...		4 mul, 4 add
     14x	standard cmul					4 mul, 2 add + table look-up
     15=14+1
     16=14+2
     17=14+3
     18=21-3	get together with 24 using just...		4 mul, 4 add
     19=21-2	get together with 23 using just...		4 mul, 4 add
     20=21-1	get together with 22 using just...		4 mul, 4 add
     21x	standard cmul					4 mul, 2 add + table look-up
     22=21+1
     23=21+2
     24=21+3						totals: 60 mul, 48 add	(4.50 ops per complex datum) + 6 table look-ups

     RADIX-32: need powers 1-31:
     Use small-tables multiply to get powers 1,2,3,7,14,21,28:
     -----------
     1x		standard cmul					4 mul, 2 add + table look-up
     2x		standard cmul					4 mul, 2 add + table look-up
     3x		standard cmul					4 mul, 2 add + table look-up
     4 =7-3	get together with 10 using just...		4 mul, 4 add
     5 =7-2	get together with 9  using just...		4 mul, 4 add
     6 =7-1	get together with 8  using just...		4 mul, 4 add
     7x		standard cmul					4 mul, 2 add + table look-up
     8 =7+1
     9 =7+2
     10=7+3
     11=14-3	get together with 17 using just...		4 mul, 4 add
     12=14-2	get together with 16 using just...		4 mul, 4 add
     13=14-1	get together with 15 using just...		4 mul, 4 add
     14x	standard cmul					4 mul, 2 add + table look-up
     15=14+1
     16=14+2
     17=14+3
     18=21-3	get together with 24 using just...		4 mul, 4 add
     19=21-2	get together with 23 using just...		4 mul, 4 add
     20=21-1	get together with 22 using just...		4 mul, 4 add
     21x	standard cmul					4 mul, 2 add + table look-up
     22=21+1
     23=21+2
     24=21+3
     25=28-3	get together with 31 using just...		4 mul, 4 add
     26=28-2	get together with 30 using just...		4 mul, 4 add
     27=28-1	get together with 29 using just...		4 mul, 4 add
     28x	standard cmul					4 mul, 2 add + table look-up
     29=28+1
     30=28+2
     31=28+3						totals: 76 mul, 62 add	(4.45 ops per complex datum) + 7 table look-ups

     RADIX-64: need powers 1-63:
     Use small-tables multiply to get powers 1,2,3,6,11,18,25,32,39,46,53,60;
     get 4 via 2^2, {5,7} via 6-+1; {8,14,9,13,10,12} via 11-+{3,2,1}, rest similar.
*/

/**********************************************************************/

	/*...Init clock counter:	*/
	ASSERT(HERE, tdiff != 0,"mers_mod_square.c: NULL tdiff ptr!");

#ifdef CTIME
	clock1 = clock();
#else
//	clock1 = time(0x0);
	clock1 = getRealTime();
#endif

	*tdiff = 0.0;
#ifdef CTIME
	dt_fwd = dt_inv = dt_sqr = 0.0;
#endif

/******************* AVX debug stuff: *******************/
	int ipad;
	double avg_abs_val = 0;
#if 0
	// Use RNG to populate data array:
	rng_isaac_init(TRUE);
	double pow2_dmult = 1024.0*128.0;	// Restrict inputs to 18 bits, which in balanced-digit representation
										// means restricting multiplier of random-inputs-in-[-1,+1] below to 2^17
	for(i = 0; i < n; i += 16) {
		ipad = i + ( (i >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		// All the inits are w.r.to an un-SIMD-rearranged ...,re,im,re,im,... pattern:
	#ifdef USE_AVX512
		a[ipad+br16[ 0]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re0
		a[ipad+br16[ 1]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im0
		a[ipad+br16[ 2]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re1
		a[ipad+br16[ 3]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im1
		a[ipad+br16[ 4]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re2
		a[ipad+br16[ 5]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im2
		a[ipad+br16[ 6]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re3
		a[ipad+br16[ 7]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im3
		a[ipad+br16[ 8]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re4
		a[ipad+br16[ 9]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im4
		a[ipad+br16[10]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re5
		a[ipad+br16[11]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im5
		a[ipad+br16[12]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re6
		a[ipad+br16[13]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im6
		a[ipad+br16[14]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re7
		a[ipad+br16[15]] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im7
	#elif defined(USE_AVX)
		a[ipad+br8[0]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re0
		a[ipad+br8[1]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im0
		a[ipad+br8[2]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re1
		a[ipad+br8[3]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im1
		a[ipad+br8[4]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re2
		a[ipad+br8[5]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im2
		a[ipad+br8[6]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re3
		a[ipad+br8[7]  ] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im3
		a[ipad+br8[0]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re4
		a[ipad+br8[1]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im4
		a[ipad+br8[2]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re5
		a[ipad+br8[3]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im5
		a[ipad+br8[4]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re6
		a[ipad+br8[5]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im6
		a[ipad+br8[6]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// re7
		a[ipad+br8[7]+8] = DNINT( rng_isaac_rand_double_norm_pm1() * pow2_dmult );	// im7
	#else
		#error Debug only enabled for AVX and above!
	#endif
  #if 0
	  if(i < 1024) {
	#ifdef USE_AVX512
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+ 0, a[ipad+br16[ 0]],ipad+ 0, a[ipad+ 0]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+ 1, a[ipad+br16[ 1]],ipad+ 1, a[ipad+ 1]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+ 2, a[ipad+br16[ 2]],ipad+ 2, a[ipad+ 2]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+ 3, a[ipad+br16[ 3]],ipad+ 3, a[ipad+ 3]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+ 4, a[ipad+br16[ 4]],ipad+ 4, a[ipad+ 4]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+ 5, a[ipad+br16[ 5]],ipad+ 5, a[ipad+ 5]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+ 6, a[ipad+br16[ 6]],ipad+ 6, a[ipad+ 6]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+ 7, a[ipad+br16[ 7]],ipad+ 7, a[ipad+ 7]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+ 8, a[ipad+br16[ 8]],ipad+ 8, a[ipad+ 8]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+ 9, a[ipad+br16[ 9]],ipad+ 9, a[ipad+ 9]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+10, a[ipad+br16[10]],ipad+10, a[ipad+10]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+11, a[ipad+br16[11]],ipad+11, a[ipad+11]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+12, a[ipad+br16[12]],ipad+12, a[ipad+12]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+13, a[ipad+br16[13]],ipad+13, a[ipad+13]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+14, a[ipad+br16[14]],ipad+14, a[ipad+14]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+15, a[ipad+br16[15]],ipad+15, a[ipad+15]);
	#elif defined(USE_AVX)
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+0  ,a[ipad+br8[0]  ],ipad+0  ,a[ipad+0  ]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+1  ,a[ipad+br8[1]  ],ipad+1  ,a[ipad+1  ]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+2  ,a[ipad+br8[2]  ],ipad+2  ,a[ipad+2  ]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+3  ,a[ipad+br8[3]  ],ipad+3  ,a[ipad+3  ]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+4  ,a[ipad+br8[4]  ],ipad+4  ,a[ipad+4  ]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+5  ,a[ipad+br8[5]  ],ipad+5  ,a[ipad+5  ]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+6  ,a[ipad+br8[6]  ],ipad+6  ,a[ipad+6  ]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+7  ,a[ipad+br8[7]  ],ipad+7  ,a[ipad+7  ]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+0+8,a[ipad+br8[0]+8],ipad+0+8,a[ipad+0+8]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+1+8,a[ipad+br8[1]+8],ipad+1+8,a[ipad+1+8]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+2+8,a[ipad+br8[2]+8],ipad+2+8,a[ipad+2+8]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+3+8,a[ipad+br8[3]+8],ipad+3+8,a[ipad+3+8]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+4+8,a[ipad+br8[4]+8],ipad+4+8,a[ipad+4+8]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+5+8,a[ipad+br8[5]+8],ipad+5+8,a[ipad+5+8]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+6+8,a[ipad+br8[6]+8],ipad+6+8,a[ipad+6+8]);
		printf("A_in[%2d] = %20.5f; SIMD: A_in[%2d] = %20.5f\n",ipad+7+8,a[ipad+br8[7]+8],ipad+7+8,a[ipad+7+8]);
	#endif
	  }
  #endif
	}
  #if 0
	for(i = 0; i < n; i += 16) {
		ipad = i + ( (i >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		avg_abs_val += fabs(a[ipad])+fabs(a[ipad+1])+fabs(a[ipad+2])+fabs(a[ipad+3])+fabs(a[ipad+4])+fabs(a[ipad+5])+fabs(a[ipad+6])+fabs(a[ipad+7])+fabs(a[ipad+8])+fabs(a[ipad+9])+fabs(a[ipad+10])+fabs(a[ipad+11])+fabs(a[ipad+12])+fabs(a[ipad+13])+fabs(a[ipad+14])+fabs(a[ipad+15]);
	}
	printf("Avg abs-val of RNG inputs = %20.10f\n",avg_abs_val);
  #endif
#endif
/********************************************************/
	/*...At the start of each iteration cycle, need to forward-weight the array of integer residue digits.
	*/
	// Mar 2017: Can skip this step if it's the start of a production test (note that any initial-residue shift
	// in such cases is handled via single-array-word forward-DWT-weighting in the Mlucas.c shift_word() function),
	// but need it if add RNG-input-setting above for debug, hence also check a[1] for nonzero:
	if(ilo || a[1]) {
		simodn = bimodn = 0;	// Init both = 0,but note for all but the 0-element they will satisfy simodn = (n - bimodn).
	#ifdef USE_FGT61
		simodnmod61 = 0;
	uint32 bimodnmod61 = 0;
	#endif
		ii     = 1;		/* Pointer to the BASE and BASEINV arrays. If n does not divide p, lowest-order digit is always a bigword (ii = 1).	*/
		for(i=0; i < n; i++)
		{
		// Apr 2014: Thanks to Stephen Searle [SMUS] for the missing-AVX-index-munge bugfix:
		#ifdef USE_AVX512
			j = (i & mask03) + br16[i&15];
		#elif defined(USE_AVX)
			j = (i & mask02) + br8[i&7];
		#elif defined(USE_SSE2)
			j = (i & mask01) + br4[i&3];
		#else
			j = i;
		#endif
			j += ( (j>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */

			l = i & (nwt-1);
			k =    i  >> nwt_bits;
			k2= (n-i) >> nwt_bits;	/* Inv-wt stuff not needed here, but gives a cheap debug check (plus, bizarrely, GCC build ~3% faster with it) */
			m =       (uint32)(simodn-si[    l])  >> 31;	/* Cast diff to unsigned int, to ensure vacated bits filled with 0 on right-shifting. */
			m2= 1 + (((uint32)(bimodn-si[nwt-l])) >> 31);
			wt    = wt0[    l]*tmp[k ]*one_half[m ];
			wtinv = wt0[nwt-l]*tmp[k2]*one_half[m2]*radix_inv;
			simodn += sw;	if(simodn >= n) simodn -= n;
			bimodn += bw;	if(bimodn >= n) bimodn -= n;
		//	if(simodn != n - bimodn) printf("I = %d: simodn[%u] != n - bimodn[%u]\n",i,simodn,n - bimodn);
		//	ASSERT(HERE, simodn == n - bimodn, "simodn != n - bimodn");	<*** cannot require this because (for i = n-1) have simodn = 0, bimodn = n,
			ASSERT(HERE, DNINT(a[j]) == a[j],"mers_mod_square.c: Input a[j] noninteger!");
		#ifdef USE_FGT61
			b[j] = a[j];						// First copy the floating-double datum into a uint64 slot...
			b[j] += (-(b[j] >> 63) & q);		// ...add q if it's < 0...
			b[j] = mul_pow2_modq(b[j], mod_wt_exp[simodnmod61]);	// ...and apply modular forward weighting.
		//simodnmod61=mod(simodn,61)								// Is there a faster way to do the mod 61?
			simodnmod61 = simodn*inv61;	simodnmod61 = simodn - 61*simodnmod61;	// Indeed there is.
		bimodnmod61 = bimodn*inv61;	bimodnmod61 = bimodn - 61*bimodnmod61;
		if(i < 20)printf("I = %2u: s,bimodnmod61 = %u,%u\n",i,bimodnmod61,simodnmod61);
		#endif
			fracmax = fabs( wt*wtinv*radix0 - 1.0 );
			ASSERT(HERE, fracmax < 1e-10, "wt*wtinv check failed!");
			a[j] *= wt;
			ii =((uint32)(sw - bimodn) >> 31);
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
	#ifdef USE_FGT61
		radix16_dif_pass1(a,b,n); break;
	#else
		radix16_dif_pass1(a,  n); break;
	#endif
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
	case 40 :
		radix40_dif_pass1(a,n); break;
	case 44 :
		radix44_dif_pass1(a,n); break;
	case 48 :
		radix48_dif_pass1(a,n); break;
	case 52 :
		radix52_dif_pass1(a,n); break;
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
	case 80 :
		radix80_dif_pass1(a,n); break;
	case 88 :
		radix88_dif_pass1(a,n); break;
	case 96 :
		radix96_dif_pass1(a,n); break;
	case 104:
		radix104_dif_pass1(a,n); break;
	case 112:
		radix112_dif_pass1(a,n); break;
	case 120:
		radix120_dif_pass1(a,n); break;
	*/
	case 128:
		radix128_dif_pass1(a,n); break;
	case 144:
		radix144_dif_pass1(a,n); break;
	case 160:
		radix160_dif_pass1(a,n); break;
	case 176:
		radix176_dif_pass1(a,n); break;
	case 192:
		radix192_dif_pass1(a,n); break;
	case 208:
		radix208_dif_pass1(a,n); break;
	case 224:
		radix224_dif_pass1(a,n); break;
	case 240 :
		radix240_dif_pass1(a,n); break;
	case 256 :
		radix256_dif_pass1(a,n); break;
	case 288:
		radix288_dif_pass1(a,n); break;
	case 320:
		radix320_dif_pass1(a,n); break;
	case 512 :
		radix512_dif_pass1(a,n); break;
	case 768 :
		radix768_dif_pass1(a,n); break;
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

	ASSERT(HERE, ihi > ilo,"mers_mod_square.c: ihi <= ilo!");

#ifdef MULTITHREAD

  #if(defined(USE_PTHREAD))
	/* Pthread stuff defined above, with other variables */
  #elif(defined(USE_OMP))
	// OpenMP currently not supported - *** To-Do: Port mers_mod_square OpenMP stuff to here ****.
	#error OpenMP currently not supported - Please recheck thread-related defines in platform.h file.
	omp_set_num_threads(NTHREADS);
	for(i=0; i < NTHREADS; i++)
	{
		num_chunks[i] = 0;
	}
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

	/* Process (radix0/2) pairs of same-sized data blocks.
	In a multithreaded implementation, process NTHREADS block pairs in parallel fashion.

	If NTHREADS does not divide (radix0/2), there will be one or more under-or-unutilized threads.
	*/

#ifdef MULTITHREAD

  #ifdef USE_PTHREAD

	/* create nchunks = radix0/2 new threads each of which will execute 'mers_process_chunk()' over some specified index
	subrange. In order to match the threads executing at any given time to the available CPUs, divide the thread execution into
	[NTHREADS] 'work shifts' ( <= #CPus), each with its threads starting and completing their work before the next shift begins:
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
				mers_process_chunk( (void*)(&tdat[j + pool_work_units]) );
			}

			struct timespec ns_time;	// We want a sleep interval of 0.1 mSec here...
			ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
			ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

			while(tpool->free_tasks_queue.num_tasks != pool_work_units) {
				// Finer-resolution, declared in <time.h>; cf. http://linux.die.net/man/2/nanosleep
				ASSERT(HERE, 0 == nanosleep(&ns_time, 0x0), "nanosleep fail!");
			}
		} else {
			for(thr_id = 0; thr_id < nchunks; ++thr_id)
			{
				mers_process_chunk( (void*)(&tdat[thr_id]) );
			}
		}

	#elif defined(USE_THREADPOOL)	// Threadpool-based dispatch for generic (non OS X) Linux

			for(thr_id = 0; thr_id < pool_work_units; ++thr_id)
			{
				task_control.data = (void*)(&tdat[thr_id]);
			//	printf("adding pool task %d\n",thr_id);
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

	#endif

	} else {
		/* Single-threaded version: */
		for(ii = 0; ii < nchunks; ++ii)
		{
			mers_process_chunk((void*)(&tdat[ii]));
		}
	}

  #endif

#else

	/* Unthreaded version: */
	for(ii = 0; ii < radix0; ii += 2)
	{
		mers_process_chunk(a,arr_scratch,n,rt0,rt1,index,block_index,ii,nradices_prim,radix_prim,ws_i,ws_j1,ws_j2,ws_j2_start,ws_k,ws_m,ws_blocklen,ws_blocklen_sum);
	}

#endif

	// Update RES_SHIFT via mod-doubling:
	MOD_ADD64(RES_SHIFT,RES_SHIFT,p,RES_SHIFT);

/*...Do the final inverse FFT pass, carry propagation and initial forward FFT pass in one fell swoop, er, swell loop...	*/

	fracmax = 0.0;

	switch(radix0)
	{
		case  5 :
			ierr =  radix5_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case  6 :
			ierr =  radix6_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case  7 :
			ierr =  radix7_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case  8 :
			ierr =  radix8_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case  9 :
			ierr =  radix9_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 10 :
			ierr = radix10_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 11 :
			ierr = radix11_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 12 :
			ierr = radix12_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 13 :
			ierr = radix13_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 14 :
			ierr = radix14_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 15 :
			ierr = radix15_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 16 :
		#ifdef USE_FGT61
			ierr = radix16_ditN_cy_dif1    (a,b,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		#else
			ierr = radix16_ditN_cy_dif1    (a,  n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		#endif
		case 18 :
			ierr = radix18_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 20 :
			ierr = radix20_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 22 :
			ierr = radix22_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 24 :
			ierr = radix24_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 26 :
			ierr = radix26_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 28 :
			ierr = radix28_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 30 :
			ierr = radix30_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 32 :
			ierr = radix32_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 36 :
			ierr = radix36_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 40 :
			ierr = radix40_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 44 :
			ierr = radix44_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 48 :
			ierr = radix48_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 52 :
			ierr = radix52_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 56 :
			ierr = radix56_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 60 :
			ierr = radix60_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 63 :
			ierr = radix63_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 64 :
			ierr = radix64_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 128 :
			ierr = radix128_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 144:
			ierr = radix144_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 160:
			ierr = radix160_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 176:
			ierr = radix176_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 192:
			ierr = radix192_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 208:
			ierr = radix208_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 224 :
			ierr = radix224_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 240 :
			ierr = radix240_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 256 :
			ierr = radix256_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 288:
			ierr = radix288_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 320:
			ierr = radix320_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 512 :
			ierr = radix512_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 768 :
			ierr = radix768_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 960 :
			ierr = radix960_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 992 :
			ierr = radix992_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 1008:
			ierr = radix1008_ditN_cy_dif1    (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 1024:
			ierr = radix1024_ditN_cy_dif1    (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
		case 4032:
			ierr = radix4032_ditN_cy_dif1    (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
	/*
		case 4096:
			ierr = radix4096_ditN_cy_dif1    (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
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
		sprintf(cbuf, "M%u Roundoff warning on iteration %8u, maxerr = %16.12f\n",(uint32)p,iter,fracmax);

	/*...Fractional parts close to 0.5 cause the program to quit.
		We put the interactive-mode errlimit close to 0.5, to let people really push the limits if they want to,
		but also store data re. less-dire ROEs, e.g. for use by timing-test 'accept this radix set?' logic: */
		if(INTERACT)
		{
			fprintf(stderr,"%s",cbuf);
			if(fracmax > 0.47 ) {
				fprintf(stderr," FATAL ERROR...Halting test of exponent %u\n",(uint32)p);
				ierr = ERR_ROUNDOFF;
				return(ierr);
			} else {
				// ***To-do:*** Accumulate number-of-worrisome-ROEs (rather than a specific iteration number) in a new global
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
			if(fracmax >= 0.4375 ) {
				// Roundoff-retry scheme detailed in comments for above fermat_mod_square() function:
				if(ROE_ITER == 0) {
					sprintf(cbuf," Retrying iteration interval to see if roundoff error is reproducible.\n");
					ROE_ITER = iter;
					ROE_VAL = fracmax;
					ierr = ERR_ROUNDOFF;
				} else if(ROE_ITER > 0) {
					if(ROE_ITER == iter && ROE_VAL == fracmax) {
						sprintf(cbuf," Roundoff error is reproducible ... switching to next-larger available FFT length and retrying.\n");
						ROE_ITER = -ROE_ITER;
						ROE_VAL = 0.0;
						ierr = ERR_ROUNDOFF;
					} else {
						sprintf(cbuf," The error is not reproducible, encountered a different fatal ROE in interval-retry ... note this is an indicator of possible data corruption. Switching to next-larger FFT length, or next-smaller, if currently running at larger-than-default FFT length.\n");
						ROE_ITER = -ROE_ITER;
						ROE_VAL = fracmax;	// Use ROE_VAL-zero-or-not to distinguish from above case
						ierr = ERR_ROUNDOFF;
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
#if 0
	for(i = 0; i < n; i += 16) {
		ipad = i + ( (i >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		avg_abs_val += fabs(a[ipad])+fabs(a[ipad+1])+fabs(a[ipad+2])+fabs(a[ipad+3])+fabs(a[ipad+4])+fabs(a[ipad+5])+fabs(a[ipad+6])+fabs(a[ipad+7])+fabs(a[ipad+8])+fabs(a[ipad+9])+fabs(a[ipad+10])+fabs(a[ipad+11])+fabs(a[ipad+12])+fabs(a[ipad+13])+fabs(a[ipad+14])+fabs(a[ipad+15]);
	}
	printf("Final-iteration: Avg abs-val of residue-array digits = %20.10f\n",avg_abs_val);
#endif

#ifdef RTIME
//	clock2 = time(0x0);
//	*tdiff += difftime(clock2 , clock1);
	clock2 = getRealTime();
	*tdiff += clock2 - clock1;
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
	#ifdef USE_FGT61
		radix16_dit_pass1(a,b,n); break;
	#else
		radix16_dit_pass1(a,  n); break;
	#endif
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
	case 40 :
		radix40_dit_pass1(a,n); break;
	case 44 :
		radix44_dit_pass1(a,n); break;
	case 48 :
		radix48_dit_pass1(a,n); break;
	case 52 :
		radix52_dit_pass1(a,n); break;
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
	case 80 :
		radix80_dit_pass1(a,n); break;
	case 88 :
		radix88_dit_pass1(a,n); break;
	case 96 :
		radix96_dit_pass1(a,n); break;
	case 104:
		radix104_dit_pass1(a,n); break;
	case 112:
		radix112_dit_pass1(a,n); break;
	case 120:
		radix120_dit_pass1(a,n); break;
*/
	case 128:
		radix128_dit_pass1(a,n); break;
	case 144:
		radix144_dit_pass1(a,n); break;
	case 160:
		radix160_dit_pass1(a,n); break;
	case 176:
		radix176_dit_pass1(a,n); break;
	case 192:
		radix192_dit_pass1(a,n); break;
	case 208:
		radix208_dit_pass1(a,n); break;
	case 224:
		radix224_dit_pass1(a,n); break;
	case 240 :
		radix240_dit_pass1(a,n); break;
	case 256 :
		radix256_dit_pass1(a,n); break;
	case 288:
		radix288_dit_pass1(a,n); break;
	case 320:
		radix320_dit_pass1(a,n); break;
	case 512 :
		radix512_dit_pass1(a,n); break;
	case 768 :
		radix768_dit_pass1(a,n); break;
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

	bimodn = 0;
#ifdef USE_FGT61
	simodnmod61 = 0;
#endif
	ii     = 1;		/* Pointer to the BASE and BASEINV arrays. If n does not divide p, lowest-order digit is always a bigword (ii = 1).	*/
	for(i=0; i < n; i++)
	{
	#ifdef USE_AVX512
		j = (i & mask03) + br16[i&15];
	#elif defined(USE_AVX)
		j = (i & mask02) + br8[i&7];
	#elif defined(USE_SSE2)
		j = (i & mask01) + br4[i&3];
	#else
		j = i;
	#endif
		j = j + ( (j>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */

		l = i & (nwt-1);
		k2= (n-i) >> nwt_bits;

		/* Cast result of subtract to unsigned int, to ensure vacated bits filled with 0 on right-shifting.
		Of course C has these screwed-up precedence rules which mean we've got to enclose the shift in (), too.
		*/
		m2= 1 + (((uint32)(bimodn-si[nwt-l])) >> 31);
		wtinv=wt0[nwt-l]*tmp[k2]*one_half[m2]*radix_inv;
		atmp = a[j]*wtinv;
		a[j] = DNINT(atmp);	// Jul 2015: This was an 'NINT' - not sure why I was using that slow macro, as opposed to DNINT
		frac_fp = fabs(a[j]-atmp);
		if(frac_fp > max_fp)
			max_fp = frac_fp;
		if(fabs(2*a[j]) > base[ii]) {
		#ifdef USE_AVX
			printf("Output a[i = %d, j = %d, jpad = %d] = %20.5f out of range: base[%d] = %15.0f\n",i,(i & mask02) + br8[i&7],j,a[j],ii,base[ii]);
		#elif defined(USE_SSE2)
			printf("Output a[i = %d, j = %d, jpad = %d] = %20.5f out of range: base[%d] = %15.0f\n",i,(i & mask01) + br4[i&3],j,a[j],ii,base[ii]);
		#else
			printf("Output a[i = %d, jpad = %d] = %20.5f out of range: base[%d] = %15.0f\n",i,j,a[j],ii,base[ii]);
		#endif
			ASSERT(HERE, 0, "Output out of range!");
		}
		bimodn += bw;	if(bimodn >= n) bimodn -= n;
		simodn = n - bimodn;
	#ifdef USE_FGT61
		b[j] += (-(b[j] >> 63) & q);		// Add q if modular digit < 0...
		b[j] = mul_pow2_modq(b[j], 61-mod_wt_exp[simodnmod61]);	// ...and apply modular inverse weighting.
		b[j] = mul_pow2_modq(b[j], 61-radix0);				// ...multiply by the modular inverse of radix0...
		if(b[j] >= qhalf) b[j] -= 0x1FFFFFFFFFFFFFFFull;		// ...and balance, i.e. put into [-qhalf, (q-1)/2].
		simodnmod61 = simodn*inv61;	simodnmod61 = simodn - 61*simodnmod61;
	// DEBUG: Compare vs floating result:
		ASSERT(HERE, a[j] == b[j], "ERROR: Floating and Modular outputs mismatch!");
	#endif
		ii =((uint32)(sw - bimodn) >> 31);
	}
	if(max_fp > 0.01)
	{
		fprintf(stderr,"%s: max_fp > 0.01! Value = %20.10f\n",func,max_fp);
		fprintf(stderr,"Check your build for inadvertent mixing of SIMD build modes!\n");
		ASSERT(HERE, max_fp < 0.01,"mers_mod_square.c: max_fp < 0.01");
	}

#if 0//CTIME
	dt_supp = dt_fwd;
	for(j=0; j<3; ++j)
	{
		if(j==1)dt_supp = dt_inv;
		if(j==2)dt_supp = dt_sqr;

		sprintf(cbuf, "Time spent inside loop[%d] =%s\n",j,get_time_str(dt_supp),
		fprintf(stderr,"%s",cbuf);
	}
#endif

	if(ierr == ERR_INTERRUPT) {	// In this need to bypass [2a] check below because ROE_ITER will be set to last-iteration-done
		return(ierr);
	}
	// Cf. [2a] in fermat_mod_square() function: The interval-retry is successful, i.e. suffers no fatal ROE.
	// [action] Prior to returning, print a "retry successful" informational and rezero ROE_ITER and ROE_VAL.
	if(!INTERACT && ROE_ITER > 0) {	// In interactive (timing-test) mode, use ROE_ITER to accumulate #iters-with-dangerous-ROEs
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
mers_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
{
	struct mers_thread_data_t* thread_arg = targ;
	int thr_id = thread_arg->tid, ii = thr_id<<1;	// Since mers-mod processes 2 data blocks per thread, ii-value = 2x unique thread identifying number
	double *a           = thread_arg->arrdat;		// Main data array
#ifdef USE_FGT61
	uint64 *b           = thread_arg->brrdat;		// Modular data in here
#endif
	int *arr_scratch    = thread_arg->arr_scratch;
	int n               = thread_arg->n;
	struct complex *rt0 = thread_arg->rt0;
	struct complex *rt1 = thread_arg->rt1;
#ifdef USE_FGT61
	uint128 *mt0        = thread_arg->mt0;
	uint128 *mt1        = thread_arg->mt1;
#endif
	int*index           = thread_arg->index;
	int*block_index     = thread_arg->block_index;
	int nradices_prim   = thread_arg->nradices_prim;
	int*radix_prim      = thread_arg->radix_prim;
	int*ws_i            = thread_arg->ws_i;
	int*ws_j1           = thread_arg->ws_j1;
	int*ws_j2           = thread_arg->ws_j2;
	int*ws_j2_start     = thread_arg->ws_j2_start;
	int*ws_k            = thread_arg->ws_k;
	int*ws_m            = thread_arg->ws_m;
	int*ws_blocklen     = thread_arg->ws_blocklen;
	int*ws_blocklen_sum = thread_arg->ws_blocklen_sum;

#else

#ifdef USE_FGT61
void mers_process_chunk(
	double a[], int arr_scratch[], int n, struct complex rt0[], struct complex rt1[], uint128 mt0[], uint128 mt1[],
	int index[], int block_index[], int ii, int nradices_prim, int radix_prim[],
	int ws_i[], int ws_j1[], int ws_j2[], int ws_j2_start[], int ws_k[], int ws_m[], int ws_blocklen[], int ws_blocklen_sum[]
)
#else
void mers_process_chunk(
	double a[], int arr_scratch[], int n, struct complex rt0[], struct complex rt1[],
	int index[], int block_index[], int ii, int nradices_prim, int radix_prim[],
	int ws_i[], int ws_j1[], int ws_j2[], int ws_j2_start[], int ws_k[], int ws_m[], int ws_blocklen[], int ws_blocklen_sum[]
)
#endif
{
	int thr_id = 0;	/* In unthreaded mode this must always = 0 */

#endif	// #ifdef MULTITHREAD

	const char func[] = "mers_process_chunk";
	int radix0 = RADIX_VEC[0];
	int i,incr,istart,j,jhi,jstart,k,koffset,l,mm;
	int init_sse2 = FALSE;	// Init-calls to various radix-pass routines presumed done prior to entry into this routine

	/* If radix0 odd and i = 0, process just one block of data, otherwise do two: */
	if(ii == 0 && (radix0 & 1))
		jhi = 1;
	else
		jhi = 2;

	for(j = 0; j < jhi; j++)
	{
		/* Get block index of the chunk of contiguous data to be processed: */
		l = block_index[ii + j];
		k    = 0;
		mm   = 1;
		incr = n/radix0;
		istart = l*incr;	/* Starting location of current data-block-to-be-processed within A-array. */
		jstart = istart + ((istart >> DAT_BITS) << PAD_BITS );

#ifdef CTIME
	clock_supp = clock();
#endif

		for(i=1; i <= NRADICES-2; i++)
		{
			/* Offset from base address of index array = L*NLOOPS = L*MM : */
			koffset = l*mm;

			switch(RADIX_VEC[i])
			{
			case  8 :
				 radix8_dif_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
			case 16 :
			#ifdef USE_FGT61
				radix16_dif_pass(&a[jstart],&b[jstart],n,rt0,rt1,mt0,mt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
			#else
				radix16_dif_pass(&a[jstart],           n,rt0,rt1        ,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
			#endif
			case 32 :
				radix32_dif_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
			default :
				sprintf(cbuf,"FATAL: radix %d not available for dif_pass. Halting...\n",RADIX_VEC[i]); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}

			k    += mm*radix0;
			mm   *= RADIX_VEC[i];
			incr /= RADIX_VEC[i];
		}	/* end i-loop. */

#ifdef CTIME
	dt_fwd += (double)(clock() - clock_supp);
#endif

	}	/* end j-loop */

	/*...Final DIF pass, wrapper/squaring and initial DIT pass are all done in-place.
	This combines data from both the l1 and l2-block, except in the case ii = 0
	for even radix0, for which the l1 = 0 and l2 = 1 blocks are processed separately within
	wrapper_square, i.e. we must call this routine a second time to process data in the l2-block.
	*/
#ifdef CTIME
	clock_supp = clock();
#endif

	for(j = 0; j < jhi; j++)
	{
		l = ii + j;

		switch(RADIX_VEC[NRADICES-1])
		{
		case 16 :
		#ifdef USE_FGT61
			radix16_wrapper_square(a,b,arr_scratch,n,radix0,rt0,rt1,mt0,mt1,nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l],init_sse2,thr_id); break;
		#else
			radix16_wrapper_square(a,  arr_scratch,n,radix0,rt0,rt1,        nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l],init_sse2,thr_id); break;
		#endif
		case 32 :
			radix32_wrapper_square(a,arr_scratch,n,radix0,rt0,rt1,nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l],init_sse2,thr_id); break;
	/*	case 64 :
			radix64_wrapper_square(a,arr_scratch,n,radix0,rt0,rt1,nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l],init_sse2,thr_id); break;
		*/
		default :
			sprintf(cbuf,"FATAL: radix %d not available for wrapper/square. Halting...\n",RADIX_VEC[NRADICES-1]); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}
	}

#ifdef CTIME
	dt_sqr += (double)(clock() - clock_supp);
#endif

	/*...Rest of inverse decimation-in-time (DIT) transform. Note that during IFFT we process the radices in reverse
	order. The first array sent to each pass routine is assumed to contain the bit-reversed floating data.	*/

	for(j = 0; j < jhi; j++)
	{
		/* Get block index of the chunk of contiguous data to be processed: */
		l = block_index[ii + j];
		ASSERT(HERE, l >= 0,"mers_mod_square.c: l >= 0");

		/* Quick-n-dirty way of generating the correct starting values of k, mm and incr -
		simply use the skeleton of the forward (DIF) loop, sans the i = NRADICES-2 pass
		(since k, mm and incr are post-incremented):
		*/
		k    = 0;
		mm   = 1;
		incr = n/radix0;

		/* calculate main-array index offset here, before incr gets changed: */
		istart = l*incr;
		jstart = istart + ((istart >> DAT_BITS) << PAD_BITS );

		for(i=1; i <= NRADICES-2; i++)
		{
			k    += mm*radix0;
			mm   *= RADIX_VEC[i];
			incr /= RADIX_VEC[i];
		}

		/* Now do the DIT loop, running the radices (and hence the values of k, mm and incr) in reverse: */

#ifdef CTIME
	clock_supp = clock();
#endif

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
			#ifdef USE_FGT61
				radix16_dit_pass(&a[jstart],&b[jstart],n,rt0,rt1,mt0,mt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
			#else
				radix16_dit_pass(&a[jstart],           n,rt0,rt1,        &index[k+koffset],mm,incr,init_sse2,thr_id); break;
			#endif
			case 32 :
				radix32_dit_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
			default :
				sprintf(cbuf,"FATAL: radix %d not available for dit_pass. Halting...\n",RADIX_VEC[i]); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}
		}	/* end i-loop */

#ifdef CTIME
	dt_inv += (double)(clock() - clock_supp);
#endif

	}	/* end j-loop */

#ifdef MULTITHREAD

	*(thread_arg->retval) = 0;	// 0 indicates successful return of current thread
//	printf("Return from Thread %d ... ", ii);
  #ifdef USE_THREADPOOL
	return 0x0;
  #else
	pthread_exit(NULL);
  #endif

#endif
}

