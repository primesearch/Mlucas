/*******************************************************************************
*                                                                              *
*   (C) 1997-2009 by Ernst W. Mayer.                                           *
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

#undef FFT_DEBUG
#define FFT_DEBUG	0

#undef RTIME
#undef CTIME

#ifdef MULTITHREAD
	#define	DBG_THREADS	0	/* Turn on to collect stats about how much work done by each thread */
	#define RTIME	/* In multithreaded mode, need to use real (wall-clock) time */
#else
	#define CTIME	/* In single-thread mode, prefer cycle-based time because of its finer granularity */
#endif

/* Extra vars for storing time spent in wrapper/dyadic-square and carry steps: */
#ifdef CTIME
	double dt_fwd, dt_inv, dt_sqr, dt_supp;
	clock_t clock_supp;
#endif

/***************/

int mers_mod_square(double a[], int arr_scratch[], int n, int ilo, int ihi, uint64 p, uint32 *err_iter, int scrnFlag, double *tdiff)
{

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

	static int radix_set_save[10] = {1000,0,0,0,0,0,0,0,0,0};
	static int radix_vec0; 	/* Stores the first element, RADIX_VEC[0], to workaround an OpemMP loop-control problem */
#if DBG_THREADS
	int num_chunks[16];		/* Collect stats about how much work done by each thread */
#endif
	static int nradices_prim,nradices_radix0,radix_prim[30];/* RADIX_PRIM stores sequence of complex FFT radices used, in terms of their prime factors.	*/
	static int *index = 0x0, *si = 0x0, *index_ptmp = 0x0, *si_ptmp = 0x0;		/* Bit-reversal index array and array storing S*I mod N values for DWT weights.	*/
	static int *block_index;				/* array storing the RADIX_VEC[0] data-block indices for pass-2 of the FFT.	*/

	/* arrays storing the index values needed for the parallel-block wrapper/square scheme: */
	static int *ws_i,*ws_j1,*ws_j2,*ws_j2_start,*ws_k,*ws_m,*ws_blocklen,*ws_blocklen_sum;
	int bimodn,simodn;					/* Mnemonic: BIMODN stands for "B times I mod N", nothing to do with bimodal.	*/
	int i,ii,ierr,iter,j,jhi,k,l,m,mm,k2,m2,l1,l2,l2_start,blocklen,blocklen_sum,outer;
	static uint64 psave=0;
	static uint32 nsave=0;
	static uint32 nwt,nwt_bits,bw,sw,bits_small;
	const  double one_half[3] = {1.0, 0.5, 0.25};		/* Needed for small-weights-tables scheme */
	static double base[2],baseinv[2],radix_inv;
	static struct complex *rt0 = 0x0, *rt1 = 0x0, *rt0_ptmp = 0x0, *rt1_ptmp = 0x0;		/* reduced-size roots of unity arrays	*/
	static double *wt0 = 0x0, *wt1 = 0x0, *tmp = 0x0, *wt0_ptmp = 0x0, *wt1_ptmp = 0x0, *tmp_ptmp = 0x0;		/* reduced-size DWT weights arrays	*/
	double fracmax,wt,wtinv;
	static int first_entry=TRUE;

#ifdef CTIME
	clock_t clock1, clock2;
#else
/* Multithreaded needs wall-clock, not CPU time: */
	time_t clock1, clock2;
#endif

	ASSERT(HERE, TRANSFORM_TYPE == REAL_WRAPPER, "mers_mod_square: Incorrect TRANSFORM_TYPE!");

/*...initialize things upon first entry */

	/*...If a new exponent, runlength or radix set, deallocate any already-allocated
	allocatable arrays and set first_entry to true:	*/

	if(p != psave || n != nsave) first_entry=TRUE;

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
		for(i = 0; i < NRADICES; i++)
		{
			if(RADIX_VEC[i] == 0)
			{
				sprintf(cbuf, "mers_mod_square: RADIX_VEC[i = %d] zero, for i < [NRADICES = %d]!",i,NRADICES);
				ASSERT(HERE, 0, cbuf);
			}
			radix_set_save[i] = RADIX_VEC[i];
		}
		for(i = NRADICES; i < 10; i++)
		{
			if(RADIX_VEC[i] != 0)
			{
				sprintf(cbuf, "mers_mod_square: RADIX_VEC[i = %d] nonzero, for i >= [NRADICES = %d]!",i,NRADICES);
				ASSERT(HERE, 0, cbuf);
			}
			radix_set_save[i] = 0;
		}

		if(0)//rt0)	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)index_ptmp     ); index_ptmp      = 0x0; index = 0x0;
			free((void *)si_ptmp        ); si_ptmp         = 0x0; si    = 0x0;
			free((void *)block_index    ); block_index     = 0x0;
			free((void *)wt0_ptmp       ); wt0_ptmp        = 0x0; wt0   = 0x0;
			free((void *)wt1_ptmp       ); wt1_ptmp        = 0x0; wt1   = 0x0;
			free((void *)tmp_ptmp       ); tmp_ptmp        = 0x0; tmp   = 0x0;
			free((void *)rt0_ptmp       ); rt0_ptmp        = 0x0; rt0   = 0x0;
			free((void *)rt1_ptmp       ); rt1_ptmp        = 0x0; rt1   = 0x0;
			free((void *)ws_i           ); ws_i            = 0x0;
			free((void *)ws_j1          ); ws_j1           = 0x0;
			free((void *)ws_j2          ); ws_j2           = 0x0;
			free((void *)ws_j2_start    ); ws_j2_start     = 0x0;
			free((void *)ws_k           ); ws_k            = 0x0;
			free((void *)ws_m           ); ws_m            = 0x0;
			free((void *)ws_blocklen    ); ws_blocklen     = 0x0;
			free((void *)ws_blocklen_sum); ws_blocklen_sum = 0x0;
		}

		N2 =n/2;		/* Complex vector length.	*/

	/* no longer needed due to above direct setting of RADIX_VEC: */
	#if 0
		/* This call sets NRADICES and the first (NRADICES) elements of RADIX_VEC: */
		int retval = get_fft_radices(n>>10, RADIX_SET, &NRADICES, RADIX_VEC, 10);

		if(retval == ERR_FFTLENGTH_ILLEGAL)
		{
			sprintf(cbuf,"ERROR: mers_mod_square: length %d = %d K not available.\n",n,n>>10);
			fp = fopen(STATFILE,"a");	fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;
			fprintf(stderr,"%s", cbuf);
			ASSERT(HERE, 0,cbuf);
		}
		else if(retval == ERR_RADIXSET_UNAVAILABLE)
		{
			/* Since the FFT length is supported, radix set 0 should be available: */
			if(get_fft_radices(n>>10, 0, &NRADICES, RADIX_VEC, 10)
			{
				sprintf(cbuf, "mers_mod_square: get_fft_radices fails with default RADIX_SET = 0 at FFT length %u K\n", n);
				fprintf(stderr,"%s",cbuf);	ASSERT(HERE, 0, cbuf);
			}

			sprintf(cbuf,"WARN: radix set %10d not available - using default.\n",RADIX_SET);
			RADIX_SET = 0;

			if(INTERACT)
			{
				fprintf(stderr,"%s",cbuf);
			}
			else
			{
				fp = fopen(STATFILE,"a");	fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;
				if(scrnFlag)                               /* Echo output to stddev */
				{
					fprintf(stderr,"%s",cbuf);
				}
			}
		}
		else if(retval != 0)
		{
			sprintf(cbuf  ,"ERROR: unknown return value %d from get_fft_radix; N = %d, kblocks = %u, radset = %u.\n", retval, n, n>>10, RADIX_SET);
			fprintf(stderr,"%s", cbuf);
			ASSERT(HERE, 0,cbuf);
		}
	#endif

		/* My array padding scheme requires N/RADIX_VEC[0] to be a power of 2, and to be >= 2^DAT_BITS, where the latter
		parameter is set in the Mdata.h file: */

		if(n%RADIX_VEC[0] != 0)
		{
			sprintf(cbuf  ,"FATAL: RADIX_VEC[0] does not divide N!\n");
			fprintf(stderr,"%s", cbuf);
			ASSERT(HERE, 0,cbuf);
		}

		/* Make sure n/RADIX_VEC[0] is a power of 2: */
		i = n/RADIX_VEC[0];
		if((i >> trailz32(i)) != 1)
		{
			sprintf(cbuf  ,"FATAL: n/RADIX_VEC[0] not a power of 2!\n");
			fprintf(stderr,"%s", cbuf);
			ASSERT(HERE, 0,cbuf);
		}

		if(DAT_BITS < 31)
		{
			/* Now make sure n/RADIX_VEC[0] is sufficiently large (unless n < 2^DAT_BITS, in which case it doesn't matter): */
			if(i < (1 << DAT_BITS))
			{
				sprintf(cbuf  ,"FATAL: n/RADIX_VEC[0] must be >= %u!\n", (1 << DAT_BITS));
				fprintf(stderr,"%s", cbuf);
				ASSERT(HERE, 0,cbuf);
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
			fp = fopen(STATFILE,"a");	fprintf(fp,"%s",cbuf);	fclose(fp);	fp = 0x0;
			if (scrnFlag)	/* Echo output to stddev */
			{
				fprintf(stderr,"%s",cbuf);
			}
		}

	/*...******Forward FFT****** permuted sincos index array is here: first, calculate the needed dimension...	*/
		k =0;
		mm=RADIX_VEC[0];			/* First radix requires no twiddle factors.	*/

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
			sprintf(cbuf  ,"FATAL: unable to allocate array INDEX in mers_mod_square.\n");
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

		switch(RADIX_VEC[0])
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
		/*
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
		case 64 :
			nradices_radix0 = 6;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
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
		case 128 :
			nradices_radix0 = 7;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		*/
		default :
			sprintf(cbuf  ,"FATAL: radix %d not available. Halting...\n",RADIX_VEC[i]);
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

			switch(RADIX_VEC[i])
			{
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
		nradices_prim = l;

		bw     = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw     = n - bw;	/* Number of smallwords.	*/

		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)RADIX_VEC[0]));

		bits_small = p/n;			/* number of bits in a smallword.	*/
		base   [0] = (double)(1 << bits_small);	base   [1] = (double)(2*base[0]);
		baseinv[0] = (double)(1.0/base[0]    );	baseinv[1] = (double)(1.0/base[1]);	/* don't need extended precision for this since both bases are powers of 2.	*/

		/*...stuff for the reduced-length DWT weights arrays is here:	*/

		/* No need for a fancy NINT here: */
		nwt_bits = (uint32)(log(sqrt(1.0*n))/log(2.0) + 0.5);
		nwt    = 1 << nwt_bits;	/* To save on storage, we calculate the first NWT weights directly and then re-use
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
		wt1 = (double *)calloc(n/nwt+RADIX_VEC[0],sizeof(double));
		tmp = (double *)calloc(n/nwt+1       ,sizeof(double));
		si  = (   int *)calloc(nwt+1         ,sizeof(   int));
		*/
		wt0_ptmp = ALLOC_DOUBLE(wt0_ptmp, nwt+1         );	if(!wt0_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array WT0 in mers_mod_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }; wt0 = ALIGN_DOUBLE(wt0_ptmp);
		wt1_ptmp = ALLOC_DOUBLE(wt1_ptmp, n/nwt+RADIX_VEC[0]);if(!wt1_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array WT1 in mers_mod_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }; wt1 = ALIGN_DOUBLE(wt1_ptmp);
		tmp_ptmp = ALLOC_DOUBLE(tmp_ptmp, n/nwt+1       );	if(!tmp_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array TMP in mers_mod_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }; tmp = ALIGN_DOUBLE(tmp_ptmp);
		si_ptmp  = ALLOC_INT   ( si_ptmp, nwt+1         );	if(!si_ptmp ){ sprintf(cbuf,"FATAL: unable to allocate array SI  in mers_mod_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }; si  = ALIGN_INT   (si_ptmp );

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
		as constructed above are accessed in strides of length n/(radix1*nwt), so it makes sense to prearrange
		them so as to replace these long strides with unit strides, and thus to be accessing contiguous data instead.	*/

		for(i=0; i<n/(nwt*RADIX_VEC[0]); i++)
		{
			for(j=i*RADIX_VEC[0], k=i; j<(i+1)*RADIX_VEC[0]; j++, k += n/(nwt*RADIX_VEC[0]))
			{
				wt1[j] = tmp[k];	/* Gather (radix1) stride-[n/(radix1*nwt)]-separated data into a contiguous block.		*/
			}
		}
		for(j=n/nwt, k=n/(nwt*RADIX_VEC[0]); j < (n/nwt+RADIX_VEC[0]); j++, k += n/(nwt*RADIX_VEC[0]))
		{
			wt1[j] = tmp[k];	/* This is so the case L = 0 comes out right.	*/
		}

		/**********************************************/
		/* Roots of unity table pairs needed for FFT: */
		/**********************************************/

		/*...The roots arrays need only be half the dimension of the weights arrays (since we need n/2 complex roots
		vs. n real weights), but have the same total storage since each entry is complex:	*/

		NRT=nwt;	NRT_BITS=nwt_bits;
		NRTM1 = NRT - 1;

		/*
		rt0 = (struct complex *)calloc(nwt      ,sizeof(struct complex));
		rt1 = (struct complex *)calloc(n/(2*nwt),sizeof(struct complex));
		*/

		/*...The rt0 array stores the (0:NRT-1)th powers of the [N2]th root of unity
		(i.e. will be accessed using the lower (NRT) bits of the integer sincos index):
		*/
		rt0_ptmp = ALLOC_COMPLEX(rt0_ptmp, NRT);
		if(!rt0_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array RT0 in mers_mod_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
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
		if(!rt1_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array RT1 in mers_mod_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
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

		/**********************************************/
		/*************** MERSENNE-ONLY: ***************/
		/**********************************************/

		/* Break the remaining portion of the FFT into RADIX_VEC[0] blocks, each of which ideally
			 should operate on a dataset which fits entirely into the L2 cache of the host machine. */

		/* 8/23/2004: Need to allocate an extra element here to account for the padding element that gets inserted when RADIX_VEC[0] is odd: */

		block_index     = (int *)malloc((RADIX_VEC[0]+1)*sizeof(int));
		if(!block_index){ sprintf(cbuf,"FATAL: unable to allocate array BLOCK_INDEX in mers_mod_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
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
			radix0 =14: want 2 sets of l1,2 pairs = {0,1},{2,13,3,12,4,11,5,10,6,9,7,8}		blocklen = 2,12		l2_start = 1,13		increments by 12
			radix0 =15: want 3 sets of l1,2 pairs = {0,-},{1,2},{3,14,4,13,5,12,6,11,7,10,8,9}	blocklen = 1,2,12	l2_start = 0,2,14,	increments by 2,12
			radix0 =16: want 4 sets of l1,2 pairs = {0,1},{2,3},{4,7,5,6},{8,15,9,14,10,13,11,12}	blocklen = 2,2,4,8	l2_start = 1,3,7,15,	increments by 2,4,8

		Rules:
			1) For even radix0 use [nradices_radix0] sets of blocks; for odd radices use [nradices_radix0+1] sets,
				the first of which consists of a single block.
			2) Outer loop over prime factors of radix0, in reverse order of which these subradices appear in forward FFT:
				- blocklen (number of l-values of) first set = 2 - radix0%2
				- blocklen of sets 2...nradices_radix0 = (number of blocks done in previous sets) * (current subradix - 1) .
			3) Throughout block processing, l1 increases monotonically from 0. Within each block, for each value of l1, the corresponding
				l2_start starts at 1 - radix0%2and increases by blocklen_newblock for each new block.

		In a multithreaded implementation, each thread will be responsible for 2 of the above blocks of data (unless radix0 is odd,
		in which case one of the threads only gets one block to chew on.) For example, if NTHREADS = 4 and radix0 = 14, the chunk-processing
		loop would occur in 2 parallel passes, with each of the 4 threads having the following blocks to process:

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

		blocklen = 2 - (RADIX_VEC[0] & 1);	/* blocklen = 2 for even RADIX_VEC[0], 1 for odd. */
		blocklen_sum=0;
		l1=0;
		l2=blocklen - 1; l2_start=l2;		/* l2_start = 1 for even RADIX_VEC[0], 0 for odd. */

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
					if(!(l >= 0 && l < RADIX_VEC[0])) { sprintf(cbuf,"ERROR 10 in mers_mod_square.c\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }

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
		ws_i            = (int *)malloc(RADIX_VEC[0]*sizeof(int));	if(!ws_i           ){ sprintf(cbuf,"FATAL: unable to allocate array WS_I            in mers_mod_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		ws_j1           = (int *)malloc(RADIX_VEC[0]*sizeof(int));	if(!ws_j1          ){ sprintf(cbuf,"FATAL: unable to allocate array WS_J1           in mers_mod_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		ws_j2           = (int *)malloc(RADIX_VEC[0]*sizeof(int));	if(!ws_j2          ){ sprintf(cbuf,"FATAL: unable to allocate array WS_J2           in mers_mod_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		ws_j2_start     = (int *)malloc(RADIX_VEC[0]*sizeof(int));	if(!ws_j2_start    ){ sprintf(cbuf,"FATAL: unable to allocate array WS_J2_START     in mers_mod_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		ws_k            = (int *)malloc(RADIX_VEC[0]*sizeof(int));	if(!ws_k           ){ sprintf(cbuf,"FATAL: unable to allocate array WS_K            in mers_mod_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		ws_m            = (int *)malloc(RADIX_VEC[0]*sizeof(int));	if(!ws_m           ){ sprintf(cbuf,"FATAL: unable to allocate array WS_M            in mers_mod_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		ws_blocklen     = (int *)malloc(RADIX_VEC[0]*sizeof(int));	if(!ws_blocklen    ){ sprintf(cbuf,"FATAL: unable to allocate array WS_BLOCKLEN     in mers_mod_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		ws_blocklen_sum = (int *)malloc(RADIX_VEC[0]*sizeof(int));	if(!ws_blocklen_sum){ sprintf(cbuf,"FATAL: unable to allocate array WS_BLOCKLEN_SUM in mers_mod_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }

		for(ii = 0; ii < RADIX_VEC[0]; ii += 2)
		{
			/*...Final DIF pass, wrapper/squaring and initial DIT pass are all done in-place.
				 This combines data from both the l1 and l2-block, except in the case ii = 0
				 for even RADIX_VEC[0], for which the l1 = 0 and l2 = 1 blocks are processed separately within
				 wrapper_square, i.e. we must call this routine a second time to process data in the l2-block.
			*/
			if(ii == 0 && !(RADIX_VEC[0] & 1))
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
						radix16_wrapper_ini(n, RADIX_VEC[0], l, nradices_prim, radix_prim, ws_i, ws_j1, ws_j2, ws_j2_start, ws_k, ws_m, ws_blocklen, ws_blocklen_sum);
						radix16_wrapper_square(a,arr_scratch,n,RADIX_VEC[0],rt0,rt1,nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l],TRUE);
						break;
					case 32 :
						radix32_wrapper_ini(n, RADIX_VEC[0], l, nradices_prim, radix_prim, ws_i, ws_j1, ws_j2, ws_j2_start, ws_k, ws_m, ws_blocklen, ws_blocklen_sum);
						radix32_wrapper_square(a,arr_scratch,n,RADIX_VEC[0],rt0,rt1,nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l],TRUE);
						break;
					/*
					case 64 :
						radix64_wrapper_ini(n, RADIX_VEC[0], l, nradices_prim, radix_prim, ws_i, ws_j1, ws_j2, ws_j2_start, ws_k, ws_m, ws_blocklen, ws_blocklen_sum);	break;
						radix64_wrapper_square(a,arr_scratch,n,RADIX_VEC[0],rt0,rt1,nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l],TRUE);
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
			fprintf(stderr, "mers_mod_square:\n");
			fprintf(stderr, " Max abs error between real*8 and real*16 computed values = %20.15f\n",         max_adiff);
			fprintf(stderr, " Max bit error between real*8 and real*16 computed values = %20.0f \n", (double)max_idiff);

			ASSERT(HERE, (max_adiff < 100*err_threshold),"Max error between real*8 and real*16 unacceptably high - quitting.");
		}
	}
	/* end of initialization sequence.	*/


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
	clock1 = time(0x0);
#endif

	*tdiff = 0.0;
#ifdef CTIME
	dt_fwd = dt_inv = dt_sqr = 0.0;
#endif

/*...At the start of each iteration cycle, need to forward weight the array of integer residue digits...	*/

/* This is a workaround for yet another DEC/Compaq TruUnix 4.0 bug - in this case, the
4.0 C compiler causes one of the elements of the A-array to be set to some small multiple
of machine epsilon, rather tha to zero. */
/*
if(ilo == 0)
{
fprintf(stderr,"test loop...\n");
	for(i=1; i < n; i++)	* Expect a[0] to be nonzero, so skip it. *
	{
		j = i + ((i >> DAT_BITS) << PAD_BITS );
		if(a[j] != 0.0)
		{
			 fprintf(stderr,"warning: nonzero a[%8d] = %25.16e\n",j,a[j]);
			 a[j] = 0.0;
		}
	}
fprintf(stderr,"done.\n");
}
*/
#if FFT_DEBUG
	sprintf(cbuf, "RT0 vector:\n");
	write_fft_debug_data((double *)rt0,0,2*nwt);

	sprintf(cbuf, "RT1 vector:\n");
	write_fft_debug_data((double *)rt1,0,n/nwt);

	sprintf(cbuf, "WT0 vector:\n");
	write_fft_debug_data((double *)wt0,0,nwt+1);

	sprintf(cbuf, "WT1 vector:\n");
	write_fft_debug_data((double *)wt1,0,n/nwt+RADIX_VEC[0]);

	rng_isaac_init(TRUE);	/* init RNG */
	for(i=0; i < n; i++)
	{
		j = i + ((i >> DAT_BITS) << PAD_BITS );
	#ifdef USE_SSE2
		j = (j & mask01) + br4[j&3];
	#endif
		a[j] = (double)((int)(100.0*rng_isaac_rand_double_norm_pm1()));
	}

	sprintf(cbuf, "Initial data vector:\n");
	write_fft_debug_data(a,0,n);
#endif

	simodn = 0;
	for(i=0; i < n; i++)
	{
		j = i + ((i >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#ifdef USE_SSE2
		ASSERT(HERE, ((j-i)&0x1) == 0,"mers_mod_square.c: Array padding non-even!");	/* Remember the stupid &#@@#^(% C precedence rules here ... bitwise &,^,| lower than arithemtic-compare. */
		j = (j & mask01) + br4[j&3];	/* As long as the array padding is always by an EVEN number of elements [which we just checked], the order here is unimportant - do parity-twiddling after padding in this case to save an extra index-store step */
	#endif

		l = i & (nwt-1);
		k = i >> nwt_bits;
		m = (uint32)(simodn-si[l]) >> 31;	/* Cast result of subtract to unsigned int, to ensure vacated bits filled with 0 on right-shifting.	*/
		wt= wt0[l]*tmp[k]*one_half[m];
		ASSERT(HERE, (double)((int32)a[j]) == a[j],"mers_mod_square.c: Input a[j] noninteger!");
		a[j] *= wt;
		simodn += sw;
		if(simodn >= n) simodn -= n;
	}

/*...and perform the initial pass of the forward transform.	*/

/*...NOTE: If the first radix to be processed is 2, 4 or 8, it is assumed that a power-of-2 FFT is being performed,
     hence no small-prime version of the corresponding pass1 routines is needed.	*/

	switch(RADIX_VEC[0])
	{
/*	case 3 :			*/
/*		 radix3_dif_pass1(a,n)	*/
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
	case 40 :
		radix40_dif_pass1(a,n); break;
	case 44 :
		radix44_dif_pass1(a,n); break;
	/*
		case 48 :
			radix48_dif_pass1(a,n); break;
		case 52 :
			radix52_dif_pass1(a,n); break;
		case 56 :
			radix56_dif_pass1(a,n); break;
		case 60 :
			radix60_dif_pass1(a,n); break;
		case 64 :
			radix64_dif_pass1(a,n); break;
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
		case 128:
			radix128_dif_pass1(a,n); break;
	*/
	default :
		sprintf(cbuf,"FATAL: radix %d not available for dif_pass1. Halting...\n",RADIX_VEC[0]);
		fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
	}

#if FFT_DEBUG
	sprintf(cbuf, "After radix[%2d]_dif_pass1:\n",RADIX_VEC[0]);
	write_fft_debug_data(a,0,n);
#endif

/**********************************************************************/

/* Main iteration loop is here. Do forward-FFT/pointwise-square/inverse-FFT, inverse weighting,
carry propagation, fractional error checking and forward weighting in same loop:
*/
ierr = 0;	/* Any return-value error code (whether fatal or not) stored here */

ASSERT(HERE, ihi > ilo,"mers_mod_square.c: ihi <= ilo!");

#ifdef MULTITHREAD

  omp_set_num_threads(NTHREADS);

  #if DBG_THREADS
	fprintf(stderr,"mers_mod_square: NTHREADS = %3d\n", NTHREADS);
	for(i=0; i < NTHREADS; i++)
	{
		num_chunks[i] = 0;
	}
  #endif
#endif

for(iter=ilo+1; iter <= ihi; iter++)
{

#if DEBUG
	fprintf(stderr,"Iter = %d\n",iter);
#endif
/*...perform the FFT-based squaring:
     Do last S-1 of S forward decimation-in-frequency transform passes.	*/

    /* Process (radix0/2) pairs of same-sized data blocks.
    In a multithreaded implementation, process NTHREADS block pairs in parallel fashion.

	If NTHREADS does not divide (radix0/2), there will be one or more under-or-unutilized threads.
    */

	radix_vec0 = RADIX_VEC[0]; 	/* Stores RADIX_VEC[0] in a scalar to work around an OpemMP loop-control problem */

#ifdef MULTITHREAD

  #pragma omp parallel for default(shared) schedule(static)

/* Options:
	#pragma omp for schedule(static)
	#pragma omp for schedule(static) wait/nowait
	#pragma omp for schedule(dynamic)
*/
#endif

    for(ii = 0; ii < radix_vec0; ii += 2)
    {
	#if DBG_THREADS
	/*	if(iter == ilo+1) fprintf(stderr,"Thread %3d : ii = %d\n", omp_get_thread_num(), ii); */
		++num_chunks[omp_get_thread_num()];
	#endif

		mers_process_chunk(a,arr_scratch,n,rt0,rt1,index,block_index,ii,nradices_prim,radix_prim,ws_i,ws_j1,ws_j2,ws_j2_start,ws_k,ws_m,ws_blocklen,ws_blocklen_sum);
    }	/* End of ii-loop. */

#ifdef MULTITHREAD
	/* end of #pragma omp parallel for{} */
#endif

/*...Do the final inverse FFT pass, carry propagation and initial forward FFT pass in one fell swoop, er, swell loop...	*/

	fracmax = 0.0;

/* Only define 2nd version of carry routine[s] with ROE checking disabled in non-SSE2 mode, as SSE2 ROE checking is cheap: */
#ifndef USE_SSE2
	if(iter <= *err_iter)	/* Determine whether to do RO error checking in carry step, depending on iteration number.	*/
	{
#endif
		switch(RADIX_VEC[0])
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
				ierr = radix16_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
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
/*	case 52 :
				ierr = radix52_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
			case 60 :
				ierr = radix60_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;*/
			default :
				sprintf(cbuf,"FATAL: radix %d not available for ditN_cy_dif1. Halting...\n",RADIX_VEC[0]); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}
#ifndef USE_SSE2
	}
	else
	{
		switch(RADIX_VEC[0])
		{
			case  5 :
				ierr =  radix5_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,         p); break;
			case  6 :
				ierr =  radix6_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,         p); break;
			case  7 :
				ierr =  radix7_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,         p); break;
			case  8 :
				ierr =  radix8_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,         p); break;
			case  9 :
				ierr =  radix9_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,         p); break;
			case 10 :
				ierr = radix10_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,         p); break;
			case 11 :
				ierr = radix11_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,         p); break;
			case 12 :
				ierr = radix12_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,         p); break;
			case 13 :
				ierr = radix13_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,         p); break;
			case 14 :
				ierr = radix14_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,         p); break;
			case 15 :
				ierr = radix15_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,         p); break;
			case 16 :
				ierr = radix16_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,         p); break;
			case 18 :
				ierr = radix18_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,         p); break;
			case 20 :
				ierr = radix20_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,         p); break;
			case 22 :
				ierr = radix22_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,         p); break;
			case 24 :
				ierr = radix24_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,         p); break;
			case 26 :
				ierr = radix26_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,         p); break;
			case 28 :
				ierr = radix28_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,         p); break;
			case 30 :
				ierr = radix30_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,         p); break;
			case 32 :
				ierr = radix32_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
			case 36 :
				ierr = radix36_ditN_cy_dif1_nochk(a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,         p); break;
			case 40 :
				ierr = radix40_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
			case 44 :
				ierr = radix44_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
/*	case 52 :
				ierr = radix52_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
			case 60 :
				ierr = radix60_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;*/
			default :
				sprintf(cbuf,"FATAL: radix %d not available for ditN_cy_dif1_nochk. Halting...\n",RADIX_VEC[0]); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}
	}
#endif	/* #ifndef USE_SSE2 */

#if FFT_DEBUG
	sprintf(cbuf, "After radix[%2d]_ditN_cy_dif1:\n",RADIX_VEC[0]);
	write_fft_debug_data(a,0,n);
#endif

	/* Nonzero remaining carries are instantly fatal: */
	if(ierr)
		return(ierr);

	/* Update     Max. Max. Error: */
	if(fracmax > MME)
		MME  = fracmax;
	/* Accumulate Avg. Max. Error: */
	if(iter > AME_ITER_START)
		AME += fracmax;

/*...Now do the fractional error check. Any fractional part in [0.4,0.6] generates a warning...	*/

	if(fracmax >= 0.4)
	{
		sprintf(cbuf, "M%u Roundoff warning on iteration %8u, maxerr = %16.12f\n",(uint32)p,iter,fracmax);

		/*...Fractional parts close to 0.5 cause the program to quit.
		We put the interactive-mode errlimit close to 0.5, to let people really push the limits if they want to:
		*/
		if(INTERACT)
		{
			fprintf(stderr,"%s",cbuf);
			if(fracmax >= 0.40625) *err_iter = p-1;	/*...If RO > 0.40625 warning issued at any point of the initial error-checked
													segment, require error checking on each iteration, even if iter > err_iter.	*/
			if(fracmax > 0.47 )
			{
				fprintf(stderr," FATAL ERROR...Halting test of exponent %u\n",(uint32)p);
				ierr = ERR_ROUNDOFF;
				return(ierr);
			}
		}
		else
		{
			fp = fopen(   OFILE,"a");
			fq = fopen(STATFILE,"a");
			fprintf(fp,"%s",cbuf);
			fprintf(fq,"%s",cbuf);
			if (scrnFlag)			/* Echo output to stddev */
			{
				fprintf(stderr,"%s",cbuf);
			}

			if(fracmax >= 0.40625) *err_iter = p-1;

/*...In range test mode, any fractional part > 0.4375 is cause for error exit.	*/
			if(fracmax > 0.4375 ) {
				fprintf(fp," FATAL ERROR...Halting test of exponent %u\n",(uint32)p);
				fprintf(fq," FATAL ERROR...Halting test of exponent %u\n",(uint32)p);
				fclose(fp);	fp = 0x0;
				fclose(fq);	fq = 0x0;
				if (scrnFlag)		/* Echo output to stddev */
				{
					fprintf(stderr," FATAL ERROR...Halting test of exponent %u\n",(uint32)p);
				}
				ierr = ERR_ROUNDOFF;
				return(ierr);
			}
			fclose(fp);	fp = 0x0;
			fclose(fq);	fq = 0x0;
		}
	}

	/*...Whew - that"s a lot of stuff that just happened. Luckily,
	computer chips don't understand the concept of "Miller time."	*/

	/* Accumulate the cycle count in a floating double on each pass to avoid problems
	with integer overflow of the clock() result, if clock_t happens to be 32-bit int on the host platform:
	*/
#ifdef CTIME
	clock2 = clock();
	*tdiff += (double)(clock2 - clock1);
	clock1 = clock2;
#endif

#if FFT_DEBUG
	exit(0);
#endif

}	/* End of main loop	*/

#ifdef RTIME
	clock2 = time(0x0);
	*tdiff += difftime(clock2 , clock1);
#endif

#if DBG_THREADS
	fprintf(stderr,"mers_mod_square: #Chunks processed by each thread: ");
	for(i=0; i < NTHREADS; i++)
	{
		fprintf(stderr,"%d[%d] ", i, num_chunks[i]);
	}
	fprintf(stderr,"\n");
#endif

/**********************************************************************/

/*...At the end of each iteration cycle, need to undo the initial DIF FFT pass...	*/

	switch(RADIX_VEC[0])
	{
	/*
		case  3 :
			 radix3_dit_pass1(a,n); break;
	*/
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
		case 40 :
			radix40_dit_pass1(a,n); break;
		case 44 :
			radix44_dit_pass1(a,n); break;
	/*
		case 48 :
			radix48_dit_pass1(a,n); break;
		case 52 :
			radix52_dit_pass1(a,n); break;
		case 56 :
			radix56_dit_pass1(a,n); break;
		case 60 :
			radix60_dit_pass1(a,n); break;
		case 64 :
			radix64_dit_pass1(a,n); break;
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
		case 128:
			radix128_dit_pass1(a,n); break;
	*/
		default :
			sprintf(cbuf,"FATAL: radix %d not available for dit_pass1. Halting...\n",RADIX_VEC[0]); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
	}

/*...and unweight the data array.	*/

	bimodn = 0;
	ii     = 1;		/* Pointer to the BASE and BASEINV arrays. If n does not divide p, lowest-order digit is always a bigword (ii = 1).	*/
	for(i=0; i < n; i++)
	{
		j = i + ((i >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#ifdef USE_SSE2
		j = (j & mask01) + br4[j&3];
	#endif

		l = i & (nwt-1);
		k2= (n-i) >> nwt_bits;

		/* Cast result of subtract to unsigned int, to ensure vacated bits filled with 0 on right-shifting.
		Of course C has these screwed-up precedence rules which mean we've got to enclose the shift in (), too.
		*/
		m2= 1 + (((uint32)(bimodn-si[nwt-l])) >> 31);
		wtinv=wt0[nwt-l]*tmp[k2]*one_half[m2]*radix_inv;
//if(i<4) { fprintf(stderr, "radix16_carry: A[%1u] * wtinv: %20.5f * %15.13f = %20.5f",j,a[j],wtinv,a[j]*wtinv); }
		a[j] = NINT(a[j]*wtinv);
//if(i<4) { fprintf(stderr, "; RND = %20.5f\n",a[j]); }
		wt = a[j];	/* Use wt for temporary storage here */
		ASSERT(HERE, fabs(wt+wt) <= base[ii], "Output out of range!");

		bimodn += bw;
		if(bimodn - n >= 0)bimodn -= n;
		ii =((uint32)(sw - bimodn) >> 31);
	}

#if 0//CTIME
	dt_supp = dt_fwd;
	for(j=0; j<3; ++j)
	{
		if(j==1)dt_supp = dt_inv;
		if(j==2)dt_supp = dt_sqr;

		dt_supp /= CLOCKS_PER_SEC;
		sprintf(cbuf, "Time spent inside loop[%d] =%2d%1d:%1d%1d:%1d%1d.%1d%1d%1d\n", j
		,(int)dt_supp/36000,((int)dt_supp%36000)/3600
		,((int)dt_supp%3600)/600,((int)dt_supp%600)/60
		,((int)dt_supp%60)/10,(int)dt_supp%10
		,(int)(10*(dt_supp-(int)dt_supp)),(int)(100*(dt_supp-(int)dt_supp))%10,(int)(1000*(dt_supp-(int)dt_supp))%10);
		fprintf(stderr,"%s",cbuf);
	}
#endif

#ifdef CTIME
	*tdiff /= CLOCKS_PER_SEC;
#endif

	return(ierr);
}

/***************/

void mers_process_chunk(double a[], int arr_scratch[], int n, struct complex rt0[], struct complex rt1[], int index[], int block_index[], int ii, int nradices_prim, int radix_prim[], int ws_i[], int ws_j1[], int ws_j2[], int ws_j2_start[], int ws_k[], int ws_m[], int ws_blocklen[], int ws_blocklen_sum[])
{
    int i,incr,istart,j,jhi,jstart,k,koffset,l,mm;

#if FFT_DEBUG
	int l0,ilo,ihi;
	int l1,klo,khi;
#endif

	/* If radix0 odd and i = 0, process just one block of data, otherwise do two: */
	if(ii == 0 && (RADIX_VEC[0] & 1))
		jhi = 1;
	else
		jhi = 2;

	for(j = 0; j < jhi; j++)
	{
		/* Get block index of the chunk of contiguous data to be processed: */
		l = block_index[ii + j];
		k    = 0;
		mm   = 1;
		incr = n/RADIX_VEC[0];

		istart = l*incr;	/* Starting location of current data-block-to-be-processed within A-array. */
		jstart = istart + ((istart >> DAT_BITS) << PAD_BITS );

	#if FFT_DEBUG
		if(j==0)
		{
			l0 = l;
			sprintf(cbuf, "On entry to mers_process_chunk: block %d:\n",l0);
			ilo = istart;	ihi = istart+incr;
			write_fft_debug_data(a,ilo,ihi);
		}
		else
		{
			l1 = l;
			sprintf(cbuf, "On entry to mers_process_chunk: block %d:\n",l1);
			klo = istart;	khi = istart+incr;
			write_fft_debug_data(a,klo,khi);
		}
	#endif

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
				 radix8_dif_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr); break;
			case 16 :
				radix16_dif_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr); break;
			case 32 :
				radix32_dif_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr); break;
			default :
				sprintf(cbuf,"FATAL: radix %d not available for dif_pass. Halting...\n",RADIX_VEC[i]); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}

			k    += mm*RADIX_VEC[0];
			mm   *= RADIX_VEC[i];
			incr /= RADIX_VEC[i];
		}	/* end i-loop. */

#ifdef CTIME
	dt_fwd += (double)(clock() - clock_supp);
#endif

	}	/* end j-loop */

#if FFT_DEBUG
	sprintf(cbuf, "After radix_dif_pass: block %d:\n",l0);
	write_fft_debug_data(a,ilo,ihi);
	sprintf(cbuf, "After radix_dif_pass: block %d:\n",l1);
	write_fft_debug_data(a,klo,khi);
#endif

	/*...Final DIF pass, wrapper/squaring and initial DIT pass are all done in-place.
	This combines data from both the l1 and l2-block, except in the case ii = 0
	for even RADIX_VEC[0], for which the l1 = 0 and l2 = 1 blocks are processed separately within
	wrapper_square, i.e. we must call this routine a second time to process data in the l2-block.
	*/
	if(ii == 0 && !(RADIX_VEC[0] & 1))
		jhi = 2;
	else
		jhi = 1;

#ifdef CTIME
	clock_supp = clock();
#endif

	for(j = 0; j < jhi; j++)
	{
		l = ii + j;

		switch(RADIX_VEC[NRADICES-1])
		{
		case 16 :
			radix16_wrapper_square(a,arr_scratch,n,RADIX_VEC[0],rt0,rt1,nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l],FALSE); break;
		case 32 :
			radix32_wrapper_square(a,arr_scratch,n,RADIX_VEC[0],rt0,rt1,nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l],FALSE); break;
		/*
		case 64 :
			radix64_wrapper_square(a,arr_scratch,n,RADIX_VEC[0],rt0,rt1,nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l],FALSE); break;
		*/
		default :
			sprintf(cbuf,"FATAL: radix %d not available for wrapper/square. Halting...\n",RADIX_VEC[NRADICES-1]); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
		}
	}

#ifdef CTIME
	dt_sqr += (double)(clock() - clock_supp);
#endif

#if FFT_DEBUG
	sprintf(cbuf, "After radix[%1d = %2d]_wrapper_square: block %d:\n",i,RADIX_VEC[NRADICES-1],l0);
	write_fft_debug_data(a,ilo,ihi);
	sprintf(cbuf, "After radix[%1d = %2d]_wrapper_square: block %d:\n",i,RADIX_VEC[NRADICES-1],l1);
	write_fft_debug_data(a,klo,khi);
#endif

	/*...Rest of inverse decimation-in-time (DIT) transform. Note that during IFFT we process the radices in reverse
	order. The first array sent to each pass routine is assumed to contain the bit-reversed floating data.	*/

	/* If radix0 odd and i = 0, process just one block of data, otherwise do two: */
	if(ii == 0 && (RADIX_VEC[0] & 1))
		jhi = 1;
	else
		jhi = 2;

	for(j = 0; j < jhi; j++)
	{
		/* Get block index of the chunk of contiguous data to be processed: */
		l = block_index[ii + j];
		ASSERT(HERE, l >= 0,"mers_mod_square.c: l >= 0");

		/* Quick-n-dirty way of generating the correct starting values of k, mm and incr -
		simply use the skeleton of the forward (DIF) loop, sans the i = NRADICES-2 pass (since k, mm and incr are post-incremented):
		*/
		k    = 0;
		mm   = 1;
		incr = n/RADIX_VEC[0];

		/* calculate main-array index offset here, before incr gets changed: */
		istart = l*incr;
		jstart = istart + ((istart >> DAT_BITS) << PAD_BITS );

		for(i=1; i <= NRADICES-2; i++)
		{
			k    += mm*RADIX_VEC[0];
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
			k    -= mm*RADIX_VEC[0];

			koffset = l*mm;

			switch(RADIX_VEC[i])
			{
			case  8 :
				 radix8_dit_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr); break;
			case 16 :
				radix16_dit_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr); break;
			case 32 :
				radix32_dit_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr); break;
			default :
				sprintf(cbuf,"FATAL: radix %d not available for dit_pass. Halting...\n",RADIX_VEC[i]); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf);
			}
		}	/* end i-loop */

#ifdef CTIME
	dt_inv += (double)(clock() - clock_supp);
#endif

	}	/* end j-loop */

#if FFT_DEBUG
	sprintf(cbuf, "After radix_dit_pass: block %d:\n",l0);
	write_fft_debug_data(a,ilo,ihi);
	sprintf(cbuf, "After radix_dit_pass: block %d:\n",l1);
	write_fft_debug_data(a,klo,khi);
#endif

}

