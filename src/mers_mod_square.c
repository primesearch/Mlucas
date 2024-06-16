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

#undef RTIME
#undef CTIME

#ifdef MULTITHREAD
	#define	DBG_THREADS	0	/* Turn on to collect stats about how much work done by each thread */
	#define RTIME	/* In multithreaded mode, need to use real (wall-clock) time */

	#include "threadpool.h"

	struct mers_thread_data_t{
		int tid;
		int*retval;
		double*arrdat;			// Main data array
		int*arr_scratch;
		int n;					// Chunksize
		struct complex*rt0;		// Roots table 1
		struct complex*rt1;		// Roots table 2
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
		uint64 fwd_fft;
		double *c;		// Optional auxiliary submul array
	};
#else
	#define CTIME	/* In single-thread mode, prefer cycle-based time because of its finer granularity */
#endif

/* Extra vars for storing time spent in wrapper/dyadic-square and carry steps: */
#ifdef CTIME
	double dt_fwd, dt_inv, dt_sqr, dt_supp;
	clock_t clock_supp;
#endif

/***************/

/* v19: To support Gerbicz check (and later, p-1) and its need for 2-input FFT-mul, added fwd_fft_only flag, which consists
of 2 separate sub-bitfields:

	Low 2 bits:	a mode_flag whose low and high bit controls whether [forward-DWT-weighting and initial-fFFT-pass] and
				undo-[forward-DWT-weighting and initial-fFFT-pass] are performed on entry and exit, respectively:
		mode_flag		description/effects:
		---------	---------------------------------------
		bit 0		Value [0|1] indicates that input [a]-array [!has | has] initial pass of fwd-FFT done upon entry.
					0 means [a]-data need to undergo initial-DWT-weighting and initial-fwd-FFT-pass, as is
					the case at the start of a typical multi-iteration interval. 1 means [a]-data are in the
					same form they have as a result of a call to the applicable fused-carry-pass routine.
		bit 1		Value [0|1] indicates that output [a]-array [has | !has] initial pass of fwd-FFT undone prior to exit.
					0 means [a]-data undergo initial-fwd-FFT-pass-rollback and initial-DWT-unweighting prior
					to return, as is the case at the end of a typical multi-iteration interval.
					1 means [a]-data are left in the same form they have as a result of the preceding call to the applicable
					fused-carry-pass routine, i.e. [a]-data left in a [forward-DWT-weighting and initial-fFFT-pass]-done state.

	High 62 bits (store these in a uint64 fwd_fft):
		o If = 0, do (ihi - ilo) FFT-based autosquarings, with input and output-handling as encoded by the mode_flag;
		o If = 1, this means "Do forward FFT only and store result in a[]". In this case we expect ilo == ihi, and in this case
			only the low bit of above-described mode_flag is meaningful, i.e. do input-handling as encoded by bit0 of mode_flag;
		o If = 2, this means "Undo first pass of forward FFT and DWT-weighting only and store result in a[]". Expect ilo == ihi,
			and in this case the above-described mode_flag is ignored. This special case is needed to support interrupt handling;
		o Otherwise: fwd_fft contains a pointer to an array b[] in already-forward-FFTed form:
				(double *)(fwd_fft_only - mode_flag) = FFT(b).
			In this case we compute FFT(a), then dyadic-multiply FFT(a) * FFT(b) (storing the result in a[]) in place
			of the usual pair_square step, then iFFT the result. The 2 bits of the mode_flag contain details about the
			state of the input [a]-array (the one needing some form of fwd-FFTing) and of the state of the resulting
			iFFTed data array, respectively.
		o [v20] If = 3, this means "dyadic-multiply of 2 inputs a and b, with both already forward-FFTed:
				(double *)a = FFT(a), (double *)(fwd_fft_only - mode_flag) = FFT(b).
			In this case we dyadic-multiply FFT(a) * FFT(b) and iFFT the product, storing the result in a[].
			Bit 0 of the mode_flag is expected = 0; bit 1 contains details about the state of the resulting
			iFFTed data array.

	In a typical PRP-test/Gerbicz-check scenario, a long sequence of PRP-test FFT mod-squarings - say iteration 0 to 10000 -
	will be interrupted every 1000 iterations to update the Gerbicz checkproduct. Here is the state of the two residues -
	PRP-test one [a] and the Gerbicz-checkproduct [b] - at the key multiple-of-1000 iterations in said interval:

		iteration
		---------	---------------------------------------
			0		Both [a],[b] in pure-integer form; do fwt-weight and initial-fwd-FFT-pass of [a], followed by 1000 FFT mod-squarings.
					mode_flag = 0 for all this since we are FFT-auto-squaring, not FFT-based 2-input modmul.
			1000	Set up for first G-checkproduct update of the current 10000-iteration interval, by calling this routine with [b] as
					the main-array argument and fwd_fft_only = 100_2: Do fwt-weight and initial-fwd-FFT-pass of [b], then fwd-FFT and exit.
					We then put a copy of [a] into a third array [c] and call this function with [c] as
					the main-array argument and fwd_fft_only = ([b] + 01_2), i.e. with mode_flag indicating
					that [c] has had the fwd-weighting and initial pass of fwd-FFT done on entry: Skip those preliminary steps,
					proceed directly to FFT-based 2-input modmul FFT(c) * FFT(b), storing the result in [c], i.e. in the
					input a[]-array argument), iFFT the result and - since the high bit of mode_flag = 0 - leave the output-array
					data in forward-weighted and initial pass-of-fwd-FFT-done form, i.e. in the form they have as a result of the
					call to the applicable fused-carry-pass routine.
					Then proceed to next 1000 autosquare-[a] iterations.
		[2-9]000	Set up for first G-checkproduct update of the current 10000-iteration interval, by calling this routine with [b] as
					the main-array argument and fwd_fft_only = 101_2: Skip fwt-weight and initial-fwd-FFT-pass of [b], proceed directly
					to remaining passes of fwd-FFT and exit. Remainder of checkproduct update is just as for the 1000-iter case:
					Again put a copy of [a] into a third array [c], &c.
			10000	Just as for the [2-9]000-iteration case, but end by calling this function with [c] as
					the main-array argument and fwd_fft_only = ([b] + 11_2), i.e. with mode_flag indicating
					that [c] has had the fwd-weighting and initial pass of fwd-FFT done on entry, but at end of the iFFT we now
					undo the initial pass-of-fwd-FFT and the forward-weighting, i.e. our output is the updated checkproduct array
					in pure-integer form.
					Then proceed to the interim-checkpointing step, consisting of savefile-writing, *without* an accompanying Gerbicz check.
			...
			k*10^6	Every millionth iteration, accompany the G-checkproduct update with an actual Gerbicz-check secondary powering. If the
					check passes, update the usual checkpoint-savefiles and also write the current residue+checkproduct to a special
					p[exponent].M savefile.
					if the check fails, skip the savefile-writing and instead roll back and restart the PRP test from the p[exponent].M file.
*/
// v20.1: Add auxiliary submul-array arg c[] to support FFT-submul a[]*(b[] - c[]) used by p-1 stage 2, with b[],c[] fwd-FFTed on entry:
int	mers_mod_square(double a[],             int arr_scratch[], int n, int ilo, int ihi, uint64 fwd_fft_only, uint64 p, int scrnFlag, double *tdiff, int update_shift, double c[])
{
	const char func[] = "mers_mod_square";
	const char*arr_sml[] = {"long","medium","short","hiacc"};	// Carry-chain length descriptors
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
	int i,ii,ierr = 0,iter,j,jhi,k,l,m,mm,k2,m2,l1,l2,l2_start,blocklen,blocklen_sum,outer;
	static uint64 psave=0;
	static uint32 nsave=0, new_runlength=0;
	static uint32 nwt,nwt_bits,bw,sw,bits_small;
	const double one_half[3] = {1.0, 0.5, 0.25};		/* Needed for small-weights-tables scheme */
	static double base[2],baseinv[2],radix_inv;
	static struct complex *rt0 = 0x0, *rt1 = 0x0, *rt0_ptmp = 0x0, *rt1_ptmp = 0x0;		/* reduced-size roots of unity arrays	*/
	static double *wt0 = 0x0, *wt1 = 0x0, *tmp = 0x0, *wt0_ptmp = 0x0, *wt1_ptmp = 0x0, *tmp_ptmp = 0x0;		/* reduced-size DWT weights arrays	*/
	double fracmax,wt,wtinv;
	double max_fp = 0.0, frac_fp, atmp;
	static int first_entry = TRUE;
	// Function pointers for DIF|DIT pass1; get set in init-block based on value of radix0:
	static void (*func_dif1)(double [], int) = 0x0;
	static void (*func_dit1)(double [], int) = 0x0;

#ifdef CTIME
	clock_t clock1, clock2;
#else
/* Multithreaded needs wall-clock, not CPU time: */
//	time_t clock1, clock2;
	double clock1, clock2;	// Jun 2014: Switched to getRealTime() code
#endif
	uint32 mode_flag = fwd_fft_only & 3;
	uint64 fwd_fft = fwd_fft_only - (uint64)mode_flag;	// fwd_fft = bits-0:1-cleared version of fwd_fft_only
	// fwd_fft_only == 0x4 yields fwd_fft = 1, "Do forward FFT only and store result in a[]"
	// fwd_fft_only == 0x8 yields fwd_fft = 2, "Undo first pass of forward FFT and DWT-weighting only and store result in a[]". Expect ilo == ihi"
	// v20: fwd_fft_only == 0xC: "fwd_fft &= ~0xC yields pointer to FFT(b)", and rely on the dyadic-mul routine to read and then clear bits 2:3.
	if(fwd_fft > 3 && fwd_fft < 0xC) {
		fwd_fft >>= 2;
	// v20: got rid of 1st constraint, so we can use a single mode_flag value in p-1 stage 2 for both vecs we want to fwd-FFT-only
	//      but input in fwd-FFT-pass-1-already-done mode and ones where we do both FFTs, input in said form and left so on return:
	//	if(fwd_fft == 1ull)
	//		ASSERT(mode_flag < 2, "Only low bit of mode_flag field may be used in this case!");
	}

	/* These came about as a result of multithreading, but now are needed whether built unthreaded or multithreaded */
	static int init_sse2 = FALSE;
	int saved_init_sse2, thr_id = -1;	// No multithread support yet.

#ifdef MULTITHREAD

	static int *thr_ret = 0x0;
	static pthread_t *thread = 0x0;
	static pthread_attr_t attr;
	static struct mers_thread_data_t *tdat = 0x0;

	// Threadpool-based dispatch:
	static int main_work_units = 0, pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)mers_process_chunk, NULL, 0x0};

#endif

	radix0 = RADIX_VEC[0];
	nchunks = radix0>>1;
	ASSERT(TRANSFORM_TYPE == REAL_WRAPPER, "mers_mod_square: Incorrect TRANSFORM_TYPE!");

/*...initialize things upon first entry */

	/*...If a new exponent, runlength or radix set, deallocate any already-allocated
	allocatable arrays and set first_entry to true:	*/

	if(n != nsave)
		new_runlength = TRUE;
	if(p != psave || new_runlength)
		first_entry=TRUE;

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
		if(!arr_scratch) {
			sprintf(cbuf, "Init portion of %s requires non-null scratch array!",func);
			ASSERT(0, cbuf);
		}
		first_entry=FALSE;
		psave = p;
		nsave = n;
		N2 =n/2;		/* Complex vector length.	*/

		for(i = 0; i < NRADICES; i++)
		{
			if(RADIX_VEC[i] == 0)
			{
				sprintf(cbuf, "%s: RADIX_VEC[i = %d] zero, for i < [NRADICES = %d]!",func,i,NRADICES);
				ASSERT(0, cbuf);
			}
			radix_set_save[i] = RADIX_VEC[i];
		}
		for(i = NRADICES; i < 10; i++)
		{
			if(RADIX_VEC[i] != 0)
			{
				sprintf(cbuf, "%s: RADIX_VEC[i = %d] nonzero, for i >= [NRADICES = %d]!",func,i,NRADICES);
				ASSERT(0, cbuf);
			}
			radix_set_save[i] = 0;
		}

		// Set function pointers for DIF|DIT pass1:
		dif1_dit1_func_name( radix0, &func_dif1, &func_dit1 );

		/* My array padding scheme requires N/radix0 to be a power of 2, and to be >= 2^DAT_BITS, where the latter
		parameter is set in the Mdata.h file: */
		if(n%radix0 != 0) {
			sprintf(cbuf  ,"ERROR: radix0 does not divide N!\n"); fprintf(stderr,"%s", cbuf); ASSERT(0,cbuf);
		}
		/* Make sure n/radix0 is a power of 2: */
		i = n/radix0;
		if((i >> trailz32(i)) != 1) {
			sprintf(cbuf  ,"ERROR: n/radix0 not a power of 2!\n"); fprintf(stderr,"%s", cbuf); ASSERT(0,cbuf);
		}

		if(DAT_BITS < 31)
		{
			/* Now make sure n/radix0 is sufficiently large (unless n < 2^DAT_BITS, in which case it doesn't matter): */
			if(i < (1 << DAT_BITS))
			{
			//	sprintf(cbuf  ,"ERROR: n/radix0 must be >= %u!\n", (1 << DAT_BITS));	fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
				// Mar 2018: Switch to 'soft' assertion error here, e.g. for timing tests at small FFT lengths:
				sprintf(cbuf  ,"n/radix0 must be >= %u! Skipping this radix combo.\n", (1 << DAT_BITS));	WARN(HERE, cbuf, "", 1); return(ERR_ASSERT);
			}

			/* We also have a lower limit on 2^DAT_BITS set by the wrapper_square routine: */
			if((1 << DAT_BITS) < 2*RADIX_VEC[NRADICES-1])
			{
				sprintf(cbuf  ,"ERROR: final FFT radix may not exceed = %u!\n", (1 << (DAT_BITS-1)));
				fprintf(stderr,"%s", cbuf);
				ASSERT(0,cbuf);
			}
		}

		sprintf(cbuf,"Using complex FFT radices*");
		char_addr = strstr(cbuf,"*");
		for(i = 0; i < NRADICES; i++)
		{
			sprintf(char_addr,"%10d",RADIX_VEC[i]); char_addr += 10;
		}; sprintf(char_addr,"\n");

		if(INTERACT)
			fprintf(stderr,"%s",cbuf);
		else
			mlucas_fprint(cbuf,0);	// In production-run mode, write to logfile

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
			sprintf(cbuf  ,"ERROR: product of radices not equal to complex vector length\n");
			fprintf(stderr,"%s", cbuf);
			ASSERT(0,cbuf);
		}

/*		index = (int *)calloc(k,sizeof(int));	*/
		index_ptmp = ALLOC_INT(index_ptmp, k);
		if(!index_ptmp)
		{
			sprintf(cbuf  ,"ERROR: unable to allocate array INDEX in %s.\n",func);
			fprintf(stderr,"%s", cbuf);
			ASSERT(0,cbuf);
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
		case 352:
			nradices_radix0 = 6;
			radix_prim[l++] =11; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 384:
			nradices_radix0 = 8;
			radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
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
			sprintf(cbuf  ,"ERROR: radix %d not available. Halting...\n",radix0);
			fprintf(stderr,"%s", cbuf);
			ASSERT(0,cbuf);
		}

		for(i = 1; i < NRADICES; i++)
		{
			/*...Allocate and initialize an index array containing MM indices...	*/
			if(i < (NRADICES-1))
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
			switch(RADIX_VEC[i]) {
			case 8 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			case 16 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			case 32 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			default :
				sprintf(cbuf  ,"ERROR: intermediate radix %d not available. Halting...\n",RADIX_VEC[i]);
				fprintf(stderr,"%s", cbuf);
				ASSERT(0,cbuf);
			}

			/* Final radix must be 16 or 32: */
			if(i == NRADICES-1 && RADIX_VEC[i] < 16)
			{
				sprintf(cbuf  ,"ERROR: final radix %d not available. Halting...\n",RADIX_VEC[i]);
				fprintf(stderr,"%s", cbuf);
				ASSERT(0,cbuf);
			}
		}
		nradices_prim = l;	for( ; l < 30; l++) { radix_prim[l] = 0; }	// Zero any higher elements which may have been previously set due
								// to use of a smoother n. (Not necessary since use nradices_prim to control array access, but nice to do.
		bw = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw = n - bw;	/* Number of smallwords.	*/
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
			sprintf(cbuf  ,"ERROR: NWT does not divide N!\n");
			fprintf(stderr,"%s", cbuf);
			ASSERT(0,cbuf);
		}

		/*...The roots arrays need only be half the dimension of the weights arrays (since we need n/2 complex roots
		vs. n real weights), but have the same total storage since each entry is complex:	*/
		/*
		wt0 = (double *)calloc(nwt+1         ,sizeof(double));
		wt1 = (double *)calloc(n/nwt+radix0  ,sizeof(double));
		tmp = (double *)calloc(n/nwt+1       ,sizeof(double));
		si  = (   int *)calloc(nwt+1         ,sizeof(   int));
		*/
		wt0_ptmp = ALLOC_DOUBLE(wt0_ptmp, nwt+1         );	if(!wt0_ptmp){ sprintf(cbuf,"ERROR: unable to allocate array WT0 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }; wt0 = ALIGN_DOUBLE(wt0_ptmp);
		wt1_ptmp = ALLOC_DOUBLE(wt1_ptmp, n/nwt+radix0  );	if(!wt1_ptmp){ sprintf(cbuf,"ERROR: unable to allocate array WT1 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }; wt1 = ALIGN_DOUBLE(wt1_ptmp);
		tmp_ptmp = ALLOC_DOUBLE(tmp_ptmp, n/nwt+1       );	if(!tmp_ptmp){ sprintf(cbuf,"ERROR: unable to allocate array TMP in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }; tmp = ALIGN_DOUBLE(tmp_ptmp);
		si_ptmp  = ALLOC_INT   ( si_ptmp, nwt+1         );	if(!si_ptmp ){ sprintf(cbuf,"ERROR: unable to allocate array SI  in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }; si  = ALIGN_INT   (si_ptmp );

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
			fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
				fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
			fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
				fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
		rt0_ptmp = ALLOC_COMPLEX(rt0_ptmp, NRT);
		if(!rt0_ptmp){ sprintf(cbuf,"ERROR: unable to allocate array RT0 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
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
			fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
			fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
				fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
				fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
		rt1_ptmp = ALLOC_COMPLEX(rt1_ptmp, n/(2*NRT));
		if(!rt1_ptmp){ sprintf(cbuf,"ERROR: unable to allocate array RT1 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
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
			fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
			fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
				fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
				fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
			}
			rt1[i].im = t1;
//if((i & 63) ==0)printf("rt1[%3u] = %20.15f, %20.15f\n",i,rt1[i].re,rt1[i].im);
			qt = qfadd(qt, qtheta);
			/* Get next terms of the recurrence:  qcnew = qcold*qr - qsold*qi,  qsnew = qcold*qi + qsold*qr. */
			/*qn = qfmul(qc, qr); qt = qfmul(qs, qi); qmul = qfsub(qn, qt);	* Store qcnew in qmul for now. */
			/*qn = qfmul(qc, qi); qt = qfmul(qs, qr); qs   = qfadd(qn, qt); qc = qmul;	*/
		}
	//	printf("%s: Complex-roots array 1 has %u elements, theta < %18.15f.\n",func,n/(2*NRT),(double)(i*theta));
	//	exit(0);
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
		if(!block_index){ sprintf(cbuf,"ERROR: unable to allocate array BLOCK_INDEX in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
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
					if(!(l >= 0 && l < radix0)) { sprintf(cbuf,"ERROR 10 in %s.c\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }

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
		if( !(ws_i            = (int *)calloc(radix0,sizeof(int))) ) { sprintf(cbuf,"ERROR: unable to allocate array WS_I            in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		if( !(ws_j1           = (int *)calloc(radix0,sizeof(int))) ) { sprintf(cbuf,"ERROR: unable to allocate array WS_J1           in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		if( !(ws_j2           = (int *)calloc(radix0,sizeof(int))) ) { sprintf(cbuf,"ERROR: unable to allocate array WS_J2           in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		if( !(ws_j2_start     = (int *)calloc(radix0,sizeof(int))) ) { sprintf(cbuf,"ERROR: unable to allocate array WS_J2_START     in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		if( !(ws_k            = (int *)calloc(radix0,sizeof(int))) ) { sprintf(cbuf,"ERROR: unable to allocate array WS_K            in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		if( !(ws_m            = (int *)calloc(radix0,sizeof(int))) ) { sprintf(cbuf,"ERROR: unable to allocate array WS_M            in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		if( !(ws_blocklen     = (int *)calloc(radix0,sizeof(int))) ) { sprintf(cbuf,"ERROR: unable to allocate array WS_BLOCKLEN     in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		if( !(ws_blocklen_sum = (int *)calloc(radix0,sizeof(int))) ) { sprintf(cbuf,"ERROR: unable to allocate array WS_BLOCKLEN_SUM in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }

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
					sprintf(cbuf,"ERROR: radix %d not available for wrapper_square. Halting...\n",RADIX_VEC[NRADICES-1]);
					fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
				}
			}
		}

		if(max_adiff > err_threshold)
		{
			fprintf(stderr, "%s:\n",func);
			fprintf(stderr, " Max abs error between real*8 and real*16 computed values = %20.15f\n",         max_adiff);
			fprintf(stderr, " Max bit error between real*8 and real*16 computed values = %20.0f \n", (double)max_idiff);
			ASSERT((max_adiff < 100*err_threshold),"Max error between real*8 and real*16 unacceptably high - quitting.");
		}

	#ifdef MULTITHREAD

		free((void *)thr_ret); thr_ret = 0x0;
		free((void *)thread ); thread  = 0x0;
		free((void *)tdat   ); tdat    = 0x0;

		thr_ret = (int *)calloc(nchunks, sizeof(int));
		thread  = (pthread_t *)calloc(nchunks, sizeof(pthread_t));
		tdat    = (struct mers_thread_data_t *)calloc(nchunks, sizeof(struct mers_thread_data_t));

		/* Initialize and set thread detached attribute */
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		/* Populate the parts of the thread-specific data structs which remain fixed from one call to the next: */
		for(i = 0; i < nchunks; ++i)
		{
			tdat[i].tid = i;
			tdat[i].retval = &thr_ret[i];
		//	tdat[i].arrdat = a;			// Main data array *** v19: Need to make updatable to support mix of mod_square and mod_mul calls ***
			tdat[i].arr_scratch = arr_scratch;
			tdat[i].n = n;					// Chunksize
			tdat[i].rt0 = rt0;	// Roots table 1
			tdat[i].rt1 = rt1;	// Roots table 2
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

	// Threadpool-based dispatch:

		// MAX_THREADS is the max. no. of threads we expect to be able to make use of, at 1 thread per core.
		ASSERT(MAX_THREADS == get_num_cores(), "MAX_THREADS not set or incorrectly set!");

		if(nchunks % NTHREADS != 0) fprintf(stderr,"%s: radix0/2 not exactly divisible by NTHREADS - This will hurt performance.\n",func);

		main_work_units = 0;
		pool_work_units = nchunks;
		ASSERT(0x0 != (tpool = threadpool_init(NTHREADS, MAX_THREADS, pool_work_units, &thread_control)), "threadpool_init failed!");
		printf("%s: Init threadpool of %d threads\n",func,NTHREADS);

	#endif	// MULTITHREAD?
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
		radix16_wrapper_square(0x0, arr_scratch, n, radix0, 0x0,0x0, nradices_prim, radix_prim, 0,0,0,0,0,0,0,0, init_sse2, thr_id, FALSE, 0x0);
		radix32_wrapper_square(0x0, arr_scratch, n, radix0, 0x0,0x0, nradices_prim, radix_prim, 0,0,0,0,0,0,0,0, init_sse2, thr_id, FALSE, 0x0);
	}
	if(init_sse2 == FALSE || (saved_init_sse2 < nchunks)) {		// New run, or need to up #threads in local-store inits
	//	init_sse2 = TRUE;
		init_sse2 = nchunks;	// Use *value* of init_sse2 to store #threads
		thr_id = -1;
		radix8_dif_pass (0x0, 0, 0x0, 0x0, 0x0, 0, 0, init_sse2, thr_id);
		radix16_dif_pass(0x0,     0, 0x0,0x0,          0, 0, 0x0, init_sse2, thr_id);
		radix32_dif_pass(0x0, 0, 0x0, 0x0, 0, 0, 0x0, init_sse2, thr_id);
		/* The dyadic-square routines need a few more params passed in init-mode than do the standalone FFT-pass routines: */
		// Dec 2017: Add rt0,rt1-pointers to the wrapper_square init calls to support USE_PRECOMPUTED_TWIDDLES build option:
		radix16_wrapper_square(0x0, arr_scratch, n, radix0, rt0,rt1, nradices_prim, radix_prim, 0,0,0,0,0,0,0,0, init_sse2, thr_id, FALSE, 0x0);
		radix32_wrapper_square(0x0, arr_scratch, n, radix0, rt0,rt1, nradices_prim, radix_prim, 0,0,0,0,0,0,0,0, init_sse2, thr_id, FALSE, 0x0);

		radix8_dit_pass (0x0, 0, 0x0, 0x0, 0, 0, 0x0, init_sse2, thr_id);
		radix16_dit_pass(0x0,     0, 0x0,0x0,          0, 0, 0x0, init_sse2, thr_id);
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

	// v19: Add support for mod_mul with one input being in precomputed fwd-FFTed form:
#ifdef MULTITHREAD
	for(i = 0; i < nchunks; ++i) { tdat[i].arrdat = a; tdat[i].fwd_fft = fwd_fft; tdat[i].c = c; }
//	printf("Thread 0: arrdat = %#" PRIX64 ", fwd_fft = %#" PRIX64 "\n",tdat[0].arrdat,tdat[0].fwd_fft);
#endif

	/*...Init clock counter:	*/
	ASSERT(tdiff != 0,"mers_mod_square.c: NULL tdiff ptr!");

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

	/*...At the start of each iteration cycle, need to forward-weight the array of integer residue digits.
	*/
	// Sep 2019: Add support for fwd_fft_only|mode_flag as described in top-of-function comments
	if(fwd_fft == 2ull) {
	//	fprintf(stderr,"[ilo,ihi] = [%u,%u], Array = %#" PRIX64 ": jumping directly to undo_initial_ffft_pass.\n",ilo,ihi,(uint64)a);
		goto undo_initial_ffft_pass;
	}
	if((mode_flag & 1) == 0)
	{
//	if(ihi<1000 && !fwd_fft)fprintf(stderr,"[ilo,ihi] = [%u,%u], Array = %#" PRIX64 ", Fwd-WT: mode_flag = %#X, a[1] = %18.10f\n",ilo,ihi,(uint64)a,mode_flag,a[1]);
		// Mar 2017: Can skip this step if it's the start of a production test (note that any initial-residue shift
		// in such cases is handled via single-array-word forward-DWT-weighting in the Mlucas.c shift_word() function),
		// but need it if add RNG-input-setting above for debug, hence also check a[1] for nonzero:
		if(ilo || a[1]) {
			simodn = bimodn = 0;	// Init both = 0,but note for all but the 0-element they will satisfy simodn = (n - bimodn).
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
			//	ASSERT(simodn == n - bimodn, "simodn != n - bimodn");	<*** cannot require this because (for i = n-1) have simodn = 0, bimodn = n,
				ASSERT(DNINT(a[j]) == a[j],"mers_mod_square.c: Input a[j] noninteger!");
				fracmax = fabs( wt*wtinv*radix0 - 1.0 );
				ASSERT(fracmax < 1e-10, "wt*wtinv check failed!");
				a[j] *= wt;
				ii =((uint32)(sw - bimodn) >> 31);
			}
		}

	/*...and perform the initial pass of the forward transform.
	NOTE: If the first radix to be processed is 2, 4 or 8, it is assumed that a power-of-2 FFT is being performed,
	hence no small-prime version of the corresponding pass1 routines is needed: */
		func_dif1(a,n);
	}	// if(low bit of mode_flag unset)

	// p-1 stage 2 restart-from-savefile needs extra option, "do just forward-weighting and fwd-FFT-pass1;
	// use the sentinel fwd_fft == (uint64)-1 to handle. That means fwd_fft_only == 0xFFFFFFFFFFFFFFFC = (uint64)-4:
	if(fwd_fft == -4ull)
		return 0;

/**********************************************************************/

	/* Main iteration loop is here. Do forward-FFT/pointwise-square/inverse-FFT, inverse weighting,
	carry propagation, fractional error checking and forward weighting in same loop:
	*/
	ierr = 0;	/* Any return-value error code (whether fatal or not) stored here */

	ASSERT(ihi > ilo,"mers_mod_square.c: ihi <= ilo!");

  #if DBG_THREADS
	fprintf(stderr,"%s: NTHREADS = %3d\n",func,NTHREADS);
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

	/* create nchunks = radix0/2 new threads each of which will execute 'mers_process_chunk()' over some specified index
	subrange. In order to match the threads executing at any given time to the available CPUs, divide the thread execution into
	[NTHREADS] 'work shifts' ( <= #CPus), each with its threads starting and completing their work before the next shift begins:
	*/
	// Threadpool-based dispatch
	for(thr_id = 0; thr_id < pool_work_units; ++thr_id)
	{
		task_control.data = (void*)(&tdat[thr_id]);
	//	printf("adding pool task %d\n",thr_id);
		threadpool_add_task(tpool, &task_control, task_is_blocking);
	//	printf("; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);
	}

//	printf("start; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);
	int	ns_retval;
	struct timespec ns_time,ns_err;	// We want a sleep interval of 0.1 mSec here...
	ns_time.tv_sec  =      0;	// (time_t)seconds - Don't use this because under OS X it's of type __darwin_time_t, which is long rather than double as under most linux distros
	ns_time.tv_nsec = 100000;	// (long)nanoseconds - Get our desired 0.1 mSec as 10^5 nSec here

//	while(tpool->tasks_queue.num_tasks != 0) {	//*** not safe, since can have #tasks == 0 with some tasks still in flight ***
	while(tpool->free_tasks_queue.num_tasks != pool_work_units) {
	//		sleep(1);	//*** too granular ***
		// Finer-resolution, declared in <time.h>; cf. http://linux.die.net/man/2/nanosleep
		ASSERT(0 == mlucas_nanosleep(&ns_time), "nanosleep re-call-on-signal fail!");
	//	printf("sleep; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);
	}
//	printf("end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

#else

	/* Unthreaded version: */
	for(ii = 0; ii < radix0; ii += 2)
	{
		mers_process_chunk(a,arr_scratch,n,rt0,rt1,index,block_index,ii,nradices_prim,radix_prim,ws_i,ws_j1,ws_j2,ws_j2_start,ws_k,ws_m,ws_blocklen,ws_blocklen_sum, fwd_fft, c);
	}

#endif

	if(fwd_fft == 1) {
	//	fprintf(stderr,"[ilo,ihi] = [%u,%u]: fwd_fft = %" PRIu64 ", mode_flag = %u: exiting after fwd-FFT.\n",ilo,ihi,fwd_fft,mode_flag);
		return 0;	// Skip carry step [and preceding inverse-FFT] in this case
	}
	// Update RES_SHIFT via mod-doubling, *** BUT ONLY IF IT'S AN AUTOSQUARE ***:
	if(update_shift) {
		MOD_ADD64(RES_SHIFT,RES_SHIFT,p,RES_SHIFT);
	}

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
			ierr = radix16_ditN_cy_dif1    (a,  n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
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
		case 352:
			ierr = radix352_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
		case 384:
			ierr = radix384_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,si,        base,baseinv,iter,&fracmax,p); break;
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
			sprintf(cbuf,"ERROR: radix %d not available for ditN_cy_dif1. Halting...\n",radix0); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
	}

	// v19: Nonzero exit carries used to be fatal, added retry-from-last-savefile handling for these
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
				fprintf(stderr," ERROR ERROR...Halting test of exponent %u\n",(uint32)p);
				ierr = ERR_ROUNDOFF;
				return(ierr);
			} else if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {
				USE_SHORT_CY_CHAIN++;	MME = 0.4;	// On switch to more-accurate-DWT-weights mode, reset MME to 0.4 to 'erase' history
				sprintf(cbuf," Reducing DWT-multipliers chain length from [%s] to [%s] in carry step to see if this prevents further excessive fractional parts.\n",arr_sml[USE_SHORT_CY_CHAIN-1],arr_sml[USE_SHORT_CY_CHAIN]);
				fprintf(stderr,"%s",cbuf);
			} else {
				// ***To-do:*** Accumulate number-of-worrisome-ROEs (rather than a specific iteration number) in a new global
			}
		}
		else
		{
			mlucas_fprint(cbuf,scrnFlag);	// Echo output to stddev

			// For <= 0.4375, before switching to next-larger FFT length, try reducing DWT-multipliers chain length:
			if(fracmax == 0.4375 && USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {
				USE_SHORT_CY_CHAIN++;
				sprintf(cbuf," Reducing DWT-multipliers chain length from [%s] to [%s] in carry step to see if this prevents further excessive fractional parts.\n",arr_sml[USE_SHORT_CY_CHAIN-1],arr_sml[USE_SHORT_CY_CHAIN]);
				mlucas_fprint(cbuf,scrnFlag);	// Echo output to stderr
			// v20: For == 0.4375, in PRP-test mode, if already at shortest chain length, simply allow further 0.4375 errors,
			// since Gerbicz check will catch rare cases of a wrong-way digit rounding, i.e. 0.4375 really being a 0.5625 in disguise:
			} else if(fracmax == 0.4375 && DO_GCHECK) {
				/* No-op */
			} else if(fracmax >= 0.4375 ) {	// already at shortest chain length
				// In range test mode, any fractional part > 0.4375 is cause for error exit, triggering switch to next-larger FFT length:
				ROE_ITER = iter;
				ROE_VAL = fracmax;
				ierr = ERR_ROUNDOFF;
			}

			if(ierr)
				return(ierr);
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
#ifndef NO_USE_SIGNALS
	// Listen for interrupts:
	if (signal(SIGINT, sig_handler) == SIG_ERR)
		fprintf(stderr,"Can't catch SIGINT.\n");
	else if (signal(SIGTERM, sig_handler) == SIG_ERR)
		fprintf(stderr,"Can't catch SIGTERM.\n");
	#ifndef __MINGW32__
	else if (signal(SIGHUP, sig_handler) == SIG_ERR)
		fprintf(stderr,"Can't catch SIGHUP.\n");
	else if (signal(SIGALRM, sig_handler) == SIG_ERR)
		fprintf(stderr,"Can't catch SIGALRM.\n");
	else if (signal(SIGUSR1, sig_handler) == SIG_ERR)
		fprintf(stderr,"Can't catch SIGUSR1.\n");
	else if (signal(SIGUSR2, sig_handler) == SIG_ERR)
		fprintf(stderr,"Can't catch SIGUSR2.\n");
	#endif
#endif
}	/* End of main for(iter....) loop	*/

// On early-exit-due-to-interrupt, decrement iter since we didn't actually do the (iter)th iteration
if(!MLUCAS_KEEP_RUNNING) {
	iter--;
//	fprintf(stderr,"%s: fwd_fft_only = %#016X, fwd_fft = %X; Caught interrupt at iter = %u; --iter = %u\n",func,fwd_fft_only,fwd_fft,iter+1,iter);
}
if(iter < ihi) {
	ASSERT(!MLUCAS_KEEP_RUNNING, "Premature iteration-loop exit due to unexpected condition!");
	ierr = ERR_INTERRUPT;
	ROE_ITER = iter;	// Function return value used for error code, so save number of last-iteration-completed-before-interrupt here
//	fprintf(stderr,"Caught signal at iter = %u; mode_flag = %#X\n",iter,mode_flag);
	mode_flag &= 0xfffffffd;	// v20: In case of interrupt-exit override any mode_flag "skip undo of initial DIF pass" setting
//	fprintf(stderr,"After ^2-toggle, mode_flag = %#X, (mode_flag >> 1) = %#X\n",mode_flag,mode_flag>>1);
}

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
	// Sep 2019: Add support for fwd_fft_only|mode_flag as described in top-of-function comments
undo_initial_ffft_pass:

	if((mode_flag >> 1) == 0)
	{
	//	fprintf(stderr,"[ilo,ihi] = [%u,%u], Array = %#" PRIX64 ", Inv-WT: mode_flag = %#X\n",ilo,ihi,(uint64)a,mode_flag);
		func_dit1(a,n);

	/*...and unweight the data array.	*/
		k = 0;	// Counts # of FP outputs which are out of range in the sense that [digit] > base/2
		bimodn = 0;
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
			k += (fabs(2*a[j]) > base[ii]);
			bimodn += bw;	if(bimodn >= n) bimodn -= n;
			simodn = n - bimodn;
			ii =((uint32)(sw - bimodn) >> 31);
		}
		if(k) {
			fprintf(stderr,"%s: Warning: %u of %u out-of-range output digits detected.\n",func,k,n);
		}
		if(max_fp > 0.01)
		{
			fprintf(stderr,"%s: max_fp > 0.01! Value = %20.10f\n",func,max_fp);
			fprintf(stderr,"Check your build for inadvertent mixing of SIMD build modes!\n");
			sprintf(cbuf,"%s: max_fp < 0.01 encountered during DWT-unweighting step!\n",func);	WARN(HERE, cbuf, "", 1); return(ERR_ASSERT);
		}
	}	// if(high bit of mode_flag unset)

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
	// *** v19: For PRP-test Must make sure we are at end of checkpoint-file iteration interval, not one of the Gerbicz-update subintervals ***
	if(!INTERACT && ROE_ITER > 0 && ihi%ITERS_BETWEEN_CHECKPOINTS == 0) {	// In interactive (timing-test) mode, use ROE_ITER to accumulate #iters-with-dangerous-ROEs
		ASSERT((ierr == 0) && (iter = ihi+1), "[2a] sanity check failed!");
		ROE_ITER = 0;
		ROE_VAL = 0.0;
		sprintf(cbuf,"Retry of iteration interval with fatal roundoff error was successful.\n");
		mlucas_fprint(cbuf,scrnFlag);	// Echo output to stddev
	}
	return(ierr);
}

/***************/

#ifdef MULTITHREAD

void*
mers_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
{
	struct mers_thread_data_t* thread_arg = targ;
	int thr_id = thread_arg->tid, ii = thr_id<<1;	// Since mers-mod processes 2 data blocks per thread, ii-value = 2x unique thread identifying number
	double *a           = thread_arg->arrdat;		// Main data array
	int *arr_scratch    = thread_arg->arr_scratch;
	int n               = thread_arg->n;
	struct complex *rt0 = thread_arg->rt0;
	struct complex *rt1 = thread_arg->rt1;
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
	uint64 fwd_fft      = thread_arg->fwd_fft;
	double*c            = thread_arg->c;

#else

void mers_process_chunk(
	double a[], int arr_scratch[], int n, struct complex rt0[], struct complex rt1[],
	int index[], int block_index[], int ii, int nradices_prim, int radix_prim[],
	int ws_i[], int ws_j1[], int ws_j2[], int ws_j2_start[], int ws_k[], int ws_m[], int ws_blocklen[], int ws_blocklen_sum[],
	uint64 fwd_fft, double c[]
)
{
	int thr_id = 0;	/* In unthreaded mode this must always = 0 */

#endif	// #ifdef MULTITHREAD

	const char func[] = "mers_process_chunk";
	int radix0 = RADIX_VEC[0];
	int i,incr,istart,j,jhi,jstart,k,koffset,l,mm;
	int init_sse2 = FALSE;	// Init-calls to various radix-pass routines presumed done prior to entry into this routine
	/*** Unlike fermat_mod_square, no need for separate cptr = c + [offset] here, since c-array offsets computed inside radix*_wrapper_square routines ***/

	/* If radix0 odd and i = 0, process just one block of data, otherwise do two: */
	if(ii == 0 && (radix0 & 1))
		jhi = 1;
	else
		jhi = 2;

	/* v20: [Bits 2:3 of fwd_fft] = 3 means "dyadic-multiply of 2 inputs a and b, with both already forward-FFTed:
		(double *)a = FFT(a), (double *)(fwd_fft_only - mode_flag) = FFT(b).
	In this case fwd_fft &= ~0xC yields pointer to FFT(b) we skip over fwd-FFT directly to
	dyadic-multiply FFT(a) * FFT(b) and iFFT the product, storing the result in a[].
	*/
  if((fwd_fft & 0xC) != 0) {
	ASSERT(((fwd_fft & 0xF) == 0xC) && ((fwd_fft>>4) != 0x0), "Bits 2:3 of fwd_fft == 3: Expect Bits 0:1 == 0 and nonzero b[] = hi60! *");
  }	else {
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
				radix16_dif_pass(&a[jstart],           n,rt0,rt1        ,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
			case 32 :
				radix32_dif_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
			default :
				sprintf(cbuf,"ERROR: radix %d not available for dif_pass. Halting...\n",RADIX_VEC[i]); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
			}

			k    += mm*radix0;
			mm   *= RADIX_VEC[i];
			incr /= RADIX_VEC[i];
		}	/* end i-loop. */
#ifdef CTIME
	dt_fwd += (double)(clock() - clock_supp);
#endif
	}	/* end j-loop */
  }	// v20: endif((fwd_fft & 0xC) != 0)

	/*...Final DIF pass, wrapper/squaring and initial DIT pass are all done in-place.
	wrapper_square, i.e. we must call this routine a second time to process data in the l2-block.
	This combines data from both the l1 and l2-block, except in the case ii = 0
	for even radix0, for which the l1 = 0 and l2 = 1 blocks are processed separately within
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
			radix16_wrapper_square(a,arr_scratch,n,radix0,rt0,rt1,nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l],init_sse2,thr_id, fwd_fft, c); break;
		case 32 :
			radix32_wrapper_square(a,arr_scratch,n,radix0,rt0,rt1,nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l],init_sse2,thr_id, fwd_fft, c); break;
	/*	case 64 :
			radix64_wrapper_square(a,arr_scratch,n,radix0,rt0,rt1,nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l],init_sse2,thr_id, fwd_fft, c); break;
		*/
		default :
			sprintf(cbuf,"ERROR: radix %d not available for wrapper/square. Halting...\n",RADIX_VEC[NRADICES-1]); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
		}
	}

#ifdef CTIME
	dt_sqr += (double)(clock() - clock_supp);
#endif

	if(fwd_fft == 1) {
	#ifdef MULTITHREAD
		*(thread_arg->retval) = 0;	// 0 indicates successful return of current thread
		return 0x0;
	#else
		return;
	#endif
	}

	/*...Rest of inverse decimation-in-time (DIT) transform. Note that during IFFT we process the radices in reverse
	order. The first array sent to each pass routine is assumed to contain the bit-reversed floating data.	*/

	for(j = 0; j < jhi; j++)
	{
		/* Get block index of the chunk of contiguous data to be processed: */
		l = block_index[ii + j];
		ASSERT(l >= 0,"mers_mod_square.c: l >= 0");

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
				radix16_dit_pass(&a[jstart],           n,rt0,rt1,        &index[k+koffset],mm,incr,init_sse2,thr_id); break;
			case 32 :
				radix32_dit_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
			default :
				sprintf(cbuf,"ERROR: radix %d not available for dit_pass. Halting...\n",RADIX_VEC[i]); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
			}
		}	/* end i-loop */

#ifdef CTIME
	dt_inv += (double)(clock() - clock_supp);
#endif

	}	/* end j-loop */

#ifdef MULTITHREAD
	*(thread_arg->retval) = 0;	// 0 indicates successful return of current thread
//	printf("Return from Thread %d ... ", ii);
	return 0x0;
#endif
}

