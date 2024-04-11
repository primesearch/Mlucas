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
		uint64 fwd_fft;
		double *c;		// Optional auxiliary submul array
	};

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
/* v20: To support Gerbicz check (and later, p-1) and its need for 2-input FFT-mul, added fwd_fft_only flag, which consists
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
int fermat_mod_square(double a[], int arr_scratch[], int n, int ilo, int ihi, uint64 fwd_fft_only, uint64 p, int scrnFlag, double *tdiff, int update_shift, double c[])
{
	const char func[] = "fermat_mod_square";
	const char*arr_sml[] = {"long","medium","short","hiacc"};	// Carry-chain length descriptors
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
	static int radix0, nchunks; 	// Store frequently-used RADIX_VEC[0] and number-of-independently-doable work units
#if DBG_THREADS
	int num_chunks[MAX_THREADS];		/* Collect stats about how much work done by each thread */
#endif
	static int nradices_prim,nradices_radix0,radix_prim[30];/* RADIX_PRIM stores sequence of complex FFT radices used, in terms of their prime factors.	*/
	static int *index = 0x0, *index_ptmp = 0x0;		/* Bit-reversal index array and array storing S*I mod N values for DWT weights.	*/

	int bimodn,i,ii,ierr = 0,iter,j,j1,j2,k,l,m,mm,k1,k2;
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
	static uint32 findex = 0;
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
	int thr_id;

#ifdef MULTITHREAD

	static int *thr_ret = 0x0;
	static pthread_t *thread = 0x0;
	static pthread_attr_t attr;
	static struct ferm_thread_data_t *tdat = 0x0;

	// Threadpool-based dispatch:
	static int main_work_units = 0, pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)fermat_process_chunk, NULL, 0x0};

#endif

#ifdef USE_IMCI512	// 1st-gen Xeon Phi - Use modified 8x8 doubles-transpose algo [1a] from util.c:test_simd_transpose_8x8()
	ASSERT(0,"Fermat-mod unsupported in k1om / IMCI-512 build mode!");
	exit(1);
#endif
	radix0 = RADIX_VEC[0];
	nchunks = radix0;
	ASSERT(TRANSFORM_TYPE == RIGHT_ANGLE, "fermat_mod_square: Incorrect TRANSFORM_TYPE!");

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
		if(!arr_scratch) {
			sprintf(cbuf, "Init portion of %s requires non-null scratch array!",func);
			ASSERT(0, cbuf);
		}
		for(i = 0; i < NRADICES; i++) {
			if(RADIX_VEC[i] == 0) {
				sprintf(cbuf, "%s: RADIX_VEC[i = %d] zero, for i < [NRADICES = %d]!",func,i,NRADICES);
				ASSERT(0, cbuf);
			}
			radix_set_save[i] = RADIX_VEC[i];
		}
		for(i = NRADICES; i < 10; i++) {
			if(RADIX_VEC[i] != 0) {
				sprintf(cbuf, "%s: RADIX_VEC[i = %d] nonzero, for i >= [NRADICES = %d]!",func,i,NRADICES);
				ASSERT(0, cbuf);
			}
			radix_set_save[i] = 0;
		}

		/*...Check that the binary exponent corresponds to a proper Fermat index: */
		findex = trailz64(p);
		ASSERT(p >> findex == 1,"fermat_mod_square.c: p >> findex == 1");

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

		if(DAT_BITS < 31) {
			/* Now make sure n/radix0 is sufficiently large (unless n < 2^DAT_BITS, in which case it doesn't matter): */
			if(i < (1 << DAT_BITS)) {
			//	sprintf(cbuf  ,"ERROR: n/radix0 must be >= %u!\n", (1 << DAT_BITS));	fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
				// Mar 2018: Switch to 'soft' assertion error here, e.g. for timing tests at small FFT lengths:
				sprintf(cbuf  ,"n/radix0 must be >= %u! Skipping this radix combo.\n", (1 << DAT_BITS));	WARN(HERE, cbuf, "", 1); return(ERR_ASSERT);
			}
			/* We also have a lower limit on 2^DAT_BITS set by the wrapper_square routine: */
			if((1 << DAT_BITS) < 2*RADIX_VEC[NRADICES-1]) {
				sprintf(cbuf  ,"ERROR: Value of DAT_BITS means final FFT radix may not exceed = %u!\n", (1 << (DAT_BITS-1)));
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
		k  = 0;
		mm = radix0;			/* First radix requires no twiddle factors.	*/
		/* We do the final DIF FFT radix within the dyadic_square routine, so store
		that block of sincos data there, where they can be merged with the wrapper sincos data:
		*/
		for(i = 1; i < NRADICES-1; i++)
		{
			k  += mm;
			mm *= RADIX_VEC[i];
		//	printf("\tR = %u, mm = %u\n",RADIX_VEC[i],mm);
		}

		if(mm*RADIX_VEC[NRADICES-1] != N2) {
			sprintf(cbuf  ,"ERROR: product of radices not equal to complex vector length\n");
			fprintf(stderr,"%s", cbuf);
			ASSERT(0,cbuf);
		}
	/*	index_ptmp = ALLOC_INT(index_ptmp, k);	<*** Jan 2020: Started getting this error here, NFC as to why:
			malloc: *** error for object 0x100802608: incorrect checksum for freed object - object was probably modified after being freed.
					*** set a breakpoint in malloc_error_break to debug
		if(!index_ptmp)
		{
			sprintf(cbuf  ,"ERROR: unable to allocate array INDEX in %s.\n",func);
			fprintf(stderr,"%s", cbuf);
			ASSERT(0,cbuf);
		}
		index = ALIGN_INT(index_ptmp);
	*/
		if(index) {
			free((void *)index); index = 0x0;
		}
		index = (int *)calloc(k,sizeof(int));
	//	printf("Alloc index[%u]...\n",k);
		/*...Forward (DIF) FFT sincos data are in bit-reversed order. We define a separate last-pass twiddles
		array within the routine wrapper_square, since that allows us to merge those nicely with the wrapper sincos data.	*/
		k = 0; l = 0; mm = 1;
		/*...First radix needs no twiddle factors, just need it for building the radix_prim array.	*/
		switch(radix0)
		{
		/*
		case 2:
			nradices_radix0 = 1;
			radix_prim[l++] = 2; break;
		case 3:
			nradices_radix0 = 1;
			radix_prim[l++] = 3; break;
		case 4:
			nradices_radix0 = 2;
			radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		*/
		case 5:
			nradices_radix0 = 1;
			radix_prim[l++] = 5; break;
		case 6:
			nradices_radix0 = 2;
			radix_prim[l++] = 3; radix_prim[l++] = 2; break;
		case 7:
			nradices_radix0 = 1;
			radix_prim[l++] = 7; break;
		case 8:
			nradices_radix0 = 3;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 9:
			nradices_radix0 = 2;
			radix_prim[l++] = 3; radix_prim[l++] = 3; break;
		case 10:
			nradices_radix0 = 2;
			radix_prim[l++] = 5; radix_prim[l++] = 2; break;
		case 11:
			nradices_radix0 = 1;
			radix_prim[l++] = 11; break;
		case 12:
			nradices_radix0 = 3;
			radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 13:
			nradices_radix0 = 1;
			radix_prim[l++] = 13; break;
		case 14:
			nradices_radix0 = 2;
			radix_prim[l++] = 7; radix_prim[l++] = 2; break;
		case 15:
			nradices_radix0 = 2;
			radix_prim[l++] = 5; radix_prim[l++] = 3; break;
		case 16:
			nradices_radix0 = 4;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 18:
			nradices_radix0 = 3;
			radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; break;
		case 20:
			nradices_radix0 = 3;
			radix_prim[l++] = 5; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 22:
			nradices_radix0 = 2;
			radix_prim[l++] =11; radix_prim[l++] = 2; break;
		case 24:
			nradices_radix0 = 4;
			radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		/*
		case 25:
			nradices_radix0 = 2;
			radix_prim[l++] = 5; radix_prim[l++] = 5; break;
		*/
		case 26:
			nradices_radix0 = 2;
			radix_prim[l++] =13; radix_prim[l++] = 2; break;
		case 28:
			nradices_radix0 = 3;
			radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 30:
			nradices_radix0 = 3;
			radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; break;
		case 32:
			nradices_radix0 = 5;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 36:
			nradices_radix0 = 4;
			radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 56:
			nradices_radix0 = 4;
			radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 60:
			nradices_radix0 = 4;
			radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 63:
			nradices_radix0 = 3;
			radix_prim[l++] = 7; radix_prim[l++] = 3; radix_prim[l++] = 3; break;
		case 64:
			nradices_radix0 = 6;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		/*
		case 72:
			nradices_radix0 = 5;
			radix_prim[l++] = 3; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 112:
			nradices_radix0 = 5;
			radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 120:
			nradices_radix0 = 5;
			radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		*/
		case 128:
			nradices_radix0 = 7;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 224:
			nradices_radix0 = 6;
			radix_prim[l++] = 7; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 240:
			nradices_radix0 = 6;
			radix_prim[l++] = 5; radix_prim[l++] = 3; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 256:
			nradices_radix0 = 8;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 512:
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
		default:
			sprintf(cbuf  ,"ERROR: radix %d not available for Fermat-mod transform. Halting...\n",RADIX_VEC[i]);
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
			switch(RADIX_VEC[i])
			{
			case 8:
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			case 16:
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			case 32:
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			default:
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

		/* Only need an IBDWT weights table if it's a non-power-of-2-length
		Fermat-mod transform, in which case the table has {odd part of N}
		distinct elements.
		*/
		nwt_bits = 0;	/* Always set = 0 for Fermat-mod */
		/* Vector length a power of 2? */
		nwt = (n >> trailz32(n));
		if(nwt == 1)
			pow2_fft = TRUE;
		else
			pow2_fft = FALSE;

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
			wt0_ptmp = ALLOC_DOUBLE(wt0_ptmp, nwt);	if(!wt0_ptmp){ sprintf(cbuf,"ERROR: unable to allocate array WT0 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }; wt0 = ALIGN_DOUBLE(wt0_ptmp);
			wt1_ptmp = ALLOC_DOUBLE(wt1_ptmp, nwt);	if(!wt1_ptmp){ sprintf(cbuf,"ERROR: unable to allocate array WT1 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }; wt1 = ALIGN_DOUBLE(wt1_ptmp);

		/******************************************************************/
		/* Crandall/Fagin weighting factors and number of bits per digit. */
		/******************************************************************/

			/* Double-check that sw*nwt (where nwt is the odd factor of N) is divisible by N: */
		//	printf("sw,nwt,n = %u,%u,%u; sw*nwt mod n = %u\n",sw,nwt,n, (uint64)sw*nwt % n);
			ASSERT((uint64)sw*nwt % n == 0,"fermat_mod_square.c: sw*nwt % n == 0");
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
				fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
			}

			for(i = 0; i < nwt; i++)
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
					fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
					fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
		if(n%NRT){ sprintf(cbuf,"ERROR: NRT does not divide N!\n"); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		NRTM1 = NRT - 1;

		/*...The rt0 array stores the (0:NRT-1)th powers of the [N2]th root of unity
		(i.e. will be accessed using the lower (NRT) bits of the integer sincos index):
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

		/*...The rt1 array stores the (0:(n/2)/NRT-1)th powers of the [(n/2)/NRT]th root of unity
		(and will be accessed using the upper bits, <NRT:31>, of the integer sincos index):
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
		for(i = 0; i < (N2/NRT); i++)
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
		rn0_ptmp = ALLOC_COMPLEX(rn0_ptmp, NRT);	if(!rn0_ptmp){ sprintf(cbuf,"ERROR: unable to allocate array RN0 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); } rn0 = ALIGN_COMPLEX(rn0_ptmp);

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
			fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
			fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
		}

		qt = QZRO;

		for(i = 0; i < NRT; i++)
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
				fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
				fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
		rn1_ptmp = ALLOC_COMPLEX(rn1_ptmp, N2/NRT);	if(!rn1_ptmp){ sprintf(cbuf,"ERROR: unable to allocate array RN1 in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); } rn1 = ALIGN_COMPLEX(rn1_ptmp);

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
			fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
			fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
		}

		qt = QZRO;

		for(i = 0; i < (N2/NRT); i++)
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
				fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
				fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
			ASSERT((max_adiff < 100*err_threshold),"Max error between real*8 and real*16 unacceptably high - quitting.");
		}

	#ifdef MULTITHREAD

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
		//	tdat[i].arrdat = a;			// Main data array *** v19: Need to make updatable to support mix of mod_square and mod_mul calls ***
			tdat[i].arr_scratch = arr_scratch;
			tdat[i].n = n;					/* Chunksize */
			tdat[i].rt0 = rt0;	/* Roots table 1 */
			tdat[i].rt1 = rt1;	/* Roots table 2 */
			tdat[i].index = index;				/* Bit-reversal index array */
			tdat[i].nradices_prim = nradices_prim;
			tdat[i].radix_prim = radix_prim;
		}

		// Threadpool-based dispatch:
		ASSERT(MAX_THREADS == get_num_cores(), "MAX_THREADS not set or incorrectly set!");

		if(radix0 % NTHREADS != 0) fprintf(stderr,"%s: radix0 not exactly divisible by NTHREADS - This will hurt performance.\n",func);

		main_work_units = 0;
		pool_work_units = radix0;
		ASSERT(0x0 != (tpool = threadpool_init(NTHREADS, MAX_THREADS, pool_work_units, &thread_control)), "threadpool_init failed!");
		printf("%s: Init threadpool of %d threads\n",func,NTHREADS);

	#endif	// MULTITHREAD?
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
		radix16_dyadic_square(0x0, arr_scratch, n, radix0, rt0, rt1, 0, nradices_prim, radix_prim, 0, init_sse2,thr_id, FALSE, 0x0);
		radix32_dyadic_square(0x0, arr_scratch, n, radix0, rt0, rt1, 0, nradices_prim, radix_prim, 0, init_sse2,thr_id, FALSE, 0x0);

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
		radix16_dyadic_square(0x0, arr_scratch, n, radix0, rt0, rt1, 0, nradices_prim, radix_prim, 0, init_sse2,thr_id, FALSE, 0x0);
		radix32_dyadic_square(0x0, arr_scratch, n, radix0, rt0, rt1, 0, nradices_prim, radix_prim, 0, init_sse2,thr_id, FALSE, 0x0);
	}
	new_runlength = FALSE;

	/* end of initialization sequence.	*/

/**********************************************************************/

	// v20: Add support for mod_mul with one input being in precomputed fwd-FFTed form:
#ifdef MULTITHREAD
	for(i = 0; i < nchunks; ++i) { tdat[i].arrdat = a; tdat[i].fwd_fft = fwd_fft; tdat[i].c = c; }
//	printf("Thread 0: arrdat = %#" PRIX64 ", fwd_fft = %#" PRIX64 "\n",tdat[0].arrdat,tdat[0].fwd_fft);
#endif

	/*...Init clock counter:	*/
	ASSERT(tdiff != 0,"fermat_mod_square.c: tdiff != 0");

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
	// Sep 2019: Add support for fwd_fft_only|mode_flag as described in top-of-function comments
	if(fwd_fft == 2ull)
		goto undo_initial_ffft_pass;
	if((mode_flag & 1) == 0)
	{
	//	fprintf(stderr,"Array = %#" PRIX64 ", Iter = %u, Fwd-WT: mode_flag = %#X, ilo = %u, a[1] = %18.10f\n",(uint64)a,ilo+1,mode_flag,ilo,a[1]);
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
					ASSERT(DNINT(a[j1]) == a[j1],"fermat_mod_square.c: Input a[j] noninteger!");
					wt = wt0[ii];
					a[j1] *= wt;
					ii += SW_DIV_N - nwt;
					ii += ( (-(int)((uint32)ii >> 31)) & nwt);
				}
				/* Odds: */
				ASSERT(ii == 0,"fermat_mod_square.c: ii == 0");
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
					ASSERT(DNINT(a[j2]) == a[j2],"fermat_mod_square.c: Input a[j] noninteger!");
					wt = wt0[ii];
					a[j2] *= wt;
					ii += SW_DIV_N - nwt;
					ii += ( (-(int)((uint32)ii >> 31)) & nwt);
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

	ASSERT(ihi > ilo,"ferm_mod_square.c: ihi <= ilo!");

  #if DBG_THREADS
	fprintf(stderr,"%s: NTHREADS = %3d\n",func,NTHREADS);
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

	/* create radix0 new threads each of which will execute 'fermat_process_chunk()' over some specified index subrange.
	In order to match the threads executing at any given time to the available CPUs, divide the thread execution
	into [NTHREADS] 'work shifts' ( <= #CPus), each with its threads starting and completing their work before the next shift
	comes online:
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
	for(ii = 0; ii < radix0; ++ii)
	{
		fermat_process_chunk(a,arr_scratch,n,rt0,rt1,index,ii,nradices_prim,radix_prim, fwd_fft, c);
	}

#endif	// threaded?

	if(fwd_fft == 1)
		return 0;	// Skip carry step [and preceding inverse-FFT] in this case

	// Update RES_SHIFT via mod-doubling, *** BUT ONLY IF IT'S AN AUTOSQUARE ***:
	if(update_shift) {
		/* Update RES_SHIFT via mod-doubling-and-add-random-bit. If initial RES_SHIFT = 0 the 'random' bit array = 0, so RES_SHIFT remains 0:
		Since - unlike LL - this PRP test involves just repeated mod-squaring, WE DON'T CARE if the sign of the residue on any
		particular intermediate iteration is flipped, since the subsequent mod-squaring will flip it back to +. The only thing
		we care about is whether the sign of the FINAL RESIDUE has flipped w.r.to that of its immediate predecessor, the
		penultimate residue, since such a final-squaring sign-flip has no ensuing mod-squaring to autocorrect it. Thus,
		if there is a final-squaring sign flip we need to unflip the sign of the resulting residue manually.
		*/
	#if 1
		i = (iter-1) % ITERS_BETWEEN_CHECKPOINTS;	// Bit we need to read...iter-counter is unit-offset w.r.to iter-interval, hence the -1
		RES_SIGN = (RES_SHIFT >= (p>>1));	// For Fermat-mod arithmetic, each modular reduction here implies a sign flip of the shifted
											// residue; but only care about whether such a sign flip occurred on the current mod-squaring.
											// Thus, RES_FLIP would probably be a better name for this global.
		MOD_ADD64(RES_SHIFT,RES_SHIFT,p,RES_SHIFT);
		RES_SHIFT += ((BASE_MULTIPLIER_BITS[i>>6] >> (i&63)) & 1);	// No mod needed on this add, since result of pvs line even and < p, which is itself even in the Fermat-mod case (p = 2^m)
		const char flip[2] = {' ','*'};
//	printf("Iter %d: shift = [%c]%" PRIu64 "\n",iter,flip[RES_SIGN],RES_SHIFT);
	#endif
	}
/*...Do the final inverse FFT pass, carry propagation and initial forward FFT pass in one fell swoop, er, swell loop...	*/

	fracmax = 0.0;
//printf("Exit(0) from %s\n",func); exit(0);
	switch(radix0)
	{
		case  5:
			ierr =  radix5_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case  6:
			ierr =  radix6_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case  7:
			ierr =  radix7_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case  8:
			ierr =  radix8_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case  9:
			ierr =  radix9_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 10:
			ierr = radix10_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 11:
			ierr = radix11_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 12:
			ierr = radix12_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 13:
			ierr = radix13_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 14:
			ierr = radix14_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 15:
			ierr = radix15_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 16:
			ierr = radix16_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 18:
			ierr = radix18_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 20:
			ierr = radix20_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 22:
			ierr = radix22_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 24:
			ierr = radix24_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 26:
			ierr = radix26_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,        base,baseinv,iter,&fracmax,p); break;
		case 28:
			ierr = radix28_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 30:
			ierr = radix30_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 32:
			ierr = radix32_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 36:
			ierr = radix36_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 56:
			ierr = radix56_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 60:
			ierr = radix60_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 63:
			ierr = radix63_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 64:
			ierr = radix64_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 128:
			ierr = radix128_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 224:
			ierr = radix224_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 240:
			ierr = radix240_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 256:
			ierr = radix256_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 512:
			ierr = radix512_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 960:
			ierr = radix960_ditN_cy_dif1     (a,n,nwt,nwt_bits,wt0,wt1,0x0,rn0,rn1,base,baseinv,iter,&fracmax,p); break;
		case 992:
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
		default:
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
		sprintf(cbuf, "F%u Roundoff warning on iteration %8u, maxerr = %16.12f\n",findex,iter,fracmax);

	/*...Fractional parts close to 0.5 cause the program to quit.
		 We put the interactive-mode errlimit close to 0.5, to let people really push the limits if they want to...	*/
		if(INTERACT)
		{
			fprintf(stderr,"%s",cbuf);
			if(fracmax > 0.47 ) {
				fprintf(stderr," ERROR ERROR...Halting test of F%u\n",findex);
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
			} else if(fracmax == 0.4375 && DO_GCHECK) {	// v21: Now do G-check for Fermat-mod, too
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
if(!MLUCAS_KEEP_RUNNING) iter--;
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
#if 0
	#warning Debug only!
	fprintf(stderr, "Skipping postprocessing ... %s: %u iter, clocks =%s [%8.4f msec/iter]. AvgMaxErr = %10.9f. MaxErr = %10.9f.\n"
	,PSTRING,(ihi-ilo),get_time_str(*tdiff),1000*get_time(*tdiff)/(ihi-ilo), AME/(ihi-ilo), MME);	exit(0);
#endif
/*...At the end of each iteration cycle, need to undo the initial DIF FFT pass...	*/

	// v20: Add support for fwd_fft_only|mode_flag as described in top-of-function comments
undo_initial_ffft_pass:
//	printf("Iter %u: ierr = %u, fwd_fft = %" PRIu64 ", mode_flag = %u\n",iter,ierr,fwd_fft,mode_flag);
	if((mode_flag >> 1) == 0)
	{
	//	fprintf(stderr,"Array = %#" PRIX64 ", Iter = %u, Inv-WT: mode_flag = %#X\n",(uint64)a,iter,mode_flag);
		func_dit1(a,n);

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
			ASSERT(ii == 0,"fermat_mod_square.c: ii == 0");
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
		k = 0;	// Counts # of FP outputs which are out of range in the sense that [digit] > base/2
		bimodn = n;
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

			ii = (bimodn > sw);					/*       i = 1 if a bigword,   0 if a smallword */
			bimodn -= sw;						/* result >= 0 if a bigword, < 0 if a smallword */
			bimodn += ( ((int)ii-1) & n);		/*       add 0 if a bigword,   N if a smallword */

			atmp  = a[j1]*radix_inv;
			a[j1] = DNINT(atmp);
			k += (fabs(2*a[j1]) > base[ii]);
			frac_fp = fabs(a[j1]-atmp);
			if(frac_fp > max_fp)
				max_fp = frac_fp;

			atmp  = a[j2]*radix_inv;
			a[j2] = DNINT(atmp);
			k += (fabs(2*a[j1]) > base[ii]);
			frac_fp = fabs(a[j2]-atmp);
			if(frac_fp > max_fp)
				max_fp = frac_fp;
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

	if(ierr == ERR_INTERRUPT) {	// In this need to bypass [2a] check below because ROE_ITER will be set to last-iteration-done
		return(ierr);
	}
	// Cf. [2a] above: The interval-retry is successful, i.e. suffers no fatal ROE.
	// [action] Prior to returning, print a "retry successful" informational and rezero ROE_ITER and ROE_VAL.
	// *** v20: For PRP-test Must make sure we are at end of checkpoint-file iteration interval, not one of the Gerbicz-update subintervals ***
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
	uint64 fwd_fft      = thread_arg->fwd_fft;
	double*c            = thread_arg->c;

#else

void fermat_process_chunk(
	double a[], int arr_scratch[], int n, struct complex rt0[], struct complex rt1[],
	int index[], int ii, int nradices_prim, int radix_prim[],
	uint64 fwd_fft, double c[])
{
	int thr_id = 0;	/* In unthreaded mode this must always = 0 */

#endif	// #ifdef MULTITHREAD

	const char func[] = "fermat_process_chunk";
	int radix0 = RADIX_VEC[0];
	int i,incr,istart,jstart,k,koffset,l,mm;
	int init_sse2 = FALSE;	// Init-calls to various radix-pass routines presumed done prior to entry into this routine
	uint64 bptr = 0x0;	// Pointer to B-array, if one is supplied in guise of the uint64 fwd_fft arg
	double*cptr = 0x0;

	l = ii;
	k    = 0;
	mm   = 1;
	// Starting location of current data-block-to-be-processed within A-array:
	incr = n/radix0;
	istart = l*incr;
	jstart = istart + ((istart >> DAT_BITS) << PAD_BITS );
	if(fwd_fft > 1)	// v20: Add support for 2-input modmul, as did for Mersenne-mod case already in v19
		bptr = (uint64)((double*)fwd_fft + jstart);
	else
		bptr = fwd_fft;
	if(c) cptr = c + jstart;
	/* v20: [Bits 2:3 of fwd_fft] = 3 means "dyadic-multiply of 2 inputs a and b, with both already forward-FFTed:
		(double *)a = FFT(a), (double *)(fwd_fft_only - mode_flag) = FFT(b).
	In this case fwd_fft &= ~0xC yields pointer to FFT(b) we skip over fwd-FFT directly to
	dyadic-multiply FFT(a) * FFT(b) and iFFT the product, storing the result in a[].
	*/
  if((fwd_fft & 0xC) != 0) {
	ASSERT(((fwd_fft & 0xF) == 0xC) && ((fwd_fft>>4) != 0x0), "Bits 2:3 of fwd_fft == 3: Expect Bits 0:1 == 0 and nonzero b[] = hi60! *");
	incr = RADIX_VEC[NRADICES-1]<<1;
  }	else {
	for(i=1; i <= NRADICES-2; i++)
	{
		/* Offset from base address of index array = L*NLOOPS = L*MM: */
		koffset = l*mm;
		switch(RADIX_VEC[i])
		{
		case  8:
			 radix8_dif_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		case 16:
			radix16_dif_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		case 32:
			radix32_dif_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		default:
			sprintf(cbuf,"ERROR: radix %d not available for dif_pass. Halting...\n",RADIX_VEC[i]); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
		}
		k    += mm*radix0;
		mm   *= RADIX_VEC[i];
		incr /= RADIX_VEC[i];
	}	/* end i-loop. */
  }	// v20: endif((fwd_fft & 0xC) != 0)

#ifdef DBG_TIME
	clock1 = clock();
	*dt_fwd += (double)(clock1 - clock0);
#endif
	/*...Final DIF pass, dyadic squaring and initial DIT pass are all done via a fused 1-pass procedure: */
	koffset = l*mm;
	/* The roots-block-re-use param mm not needed for innermost pass, since there each set of inputs gets its own set of roots: */
	switch(RADIX_VEC[NRADICES-1])
	{
		case 16:
			radix16_dyadic_square(&a[jstart],arr_scratch,n,radix0,rt0,rt1,ii,nradices_prim,radix_prim,incr,init_sse2,thr_id, bptr, cptr); break;
		case 32:
			radix32_dyadic_square(&a[jstart],arr_scratch,n,radix0,rt0,rt1,ii,nradices_prim,radix_prim,incr,init_sse2,thr_id, bptr, cptr); break;
		default:
			sprintf(cbuf,"ERROR: radix %d not available for wrapper/square. Halting...\n",RADIX_VEC[NRADICES-1]); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
	}
#ifdef DBG_TIME
	clock2 = clock();
	*dt_sqr += (double)(clock2 - clock1);
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
		case  8:
			 radix8_dit_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		case 16:
			radix16_dit_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		case 32:
			radix32_dit_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		default:
			sprintf(cbuf,"ERROR: radix %d not available for dit_pass. Halting...\n",RADIX_VEC[i]); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
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
	return 0x0;
#endif
}

