/*******************************************************************************
*                                                                              *
*   (C) 1997-2012 by Ernst W. Mayer.                                           *
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
#include "genFFT_mul.h"

/*** DEBUG NOTES: ***/

/* Set FFT_DEBUG = 1 (here or globally in in masterdefs.h) and call with FORWARD_FFT_ONLY = FALSE
	to invoke a simple FFT/IFFT-returns-original-inputs test
	(sans weighting and dyadic squaring) using pseudorandom inputs:
*/
#define FFT_DEBUG	0	// Setting this nonzero breaks the top-level GCD test of this routine, is just for
						// basic FFT test, in the sense that if the function returns sans ASSERT, basic FFT is OK.
#if FFT_DEBUG
	char *char_addr;
#endif

/* Set this = 1 if you suspect something is qoing awry with the DWT weights and/or FFT sincos data. */
#define QFDBG 0

/***************/

/*
If INIT_ARRAYS = TRUE, must have valid length-n x-array pointer, which is used for
scratch space for the init and whose contents are thus assumed to be destroyed!
*/
void pairFFT_mul(double x[], double y[], int n, int INIT_ARRAYS, int FORWARD_FFT_ONLY)
{
/*
	Depending on the precise way the function is called, uses a length-n
	paired-real-vector FFT (i.e. length-n/2 complex FFT of a pair of
	length-n/2 even/odd-interleaved real input vectors) to perform either a pair of
	in-place multiprecision (MP) integer linear combinations of the form

		( u' )   ( a  -b ) ( u )
		(    ) = (       )*(   )		(*)
		( v' )   ( c  -d ) ( v )

	or a pair of multiprecision real-vector-times-real-vector multiplies of form

		( u' )   ( a   0 ) ( u )
		(    ) = (       )*(   )		(**)
		( v' )   ( 0   b ) ( v )

	This function can be called in either init or active-data mode,
	based on the boolean value of the input INIT_ARRAYS:

	If INIT_ARRAYS = TRUE, the function simply uses the input (n) to init the FFT sincos
	data tables and returns; the x-array is assumed to be valid but uninitialized, and is
	used for scratch space for the init procedures. All other input arguments are ignored.
	If FFT_DEBUG is set, the y-input is assumed to contain nontrivial test data for the FFT
	routines, and instead of immediately returning after doing the table initializations,
	the function copies the y-data into the scratch x-array, does a fwd+inv FFT on them,
	and compares the result to the original y-inputs.

	Such an init call must always precede any calls in which actual u'/v' computation occurs,
	the latter being effected by calling with INIT_ARRAYS = FALSE.

	In the above formula (*), a, b, c, d, u, v are nonnegative multiword integers (input in
	extracted-to-double form), that is, converted to base-W balanced-digit floating-point
	form with some wordsize W allowing for direct length-n FFT-based multiply. (W = 2^16
	is a convenient and safe wordsize for most applications). It is assumed that a given
	set of matrix coefficients a/b/c/d in (*) will generally be applied to multiple sets
	of u/v data. Thus it is desirable to precompute and store the forward-transforms of
	a/b/c/d just once and then reuse same as many times as needed. For that purpose
	(and to minimize code duplication) we provide 2 active-data modes (as opposed to
	precomputed-data initialization mode) of calling this function, based on the boolean
	value of the input FORWARD_FFT_ONLY:

	FORWARD_FFT_ONLY = TRUE: The a/b and c/d multiplier data are assumed to be input in
	packed even/odd form in the x and y-arrays, respectively - if both array pointers are
	valid we will set up to compute (*), otherwise if only ab_mul[] is valid we set up for (**).
	The forward FFTs of these input vectors are computed in-place. Such an init call must precede
	any call to compute (*) or (**), i.e. a call with FORWARD_FFT_ONLY = FALSE
	and with the u/v data input in the x-array.

	FORWARD_FFT_ONLY = FALSE: The u/v multiplicands are assumed to be input in
	packed even/odd form in the x-array (y is *required* to be null). FFT(u/v) is
	computed in in-place form, and the appropriate linear combination (*) or (**) (depending
	on whether the static pointer cd_mul is set or not) is obtained by combining FFT(u/v)
	with FFT(a/b, c/d), followed by in-place IFFT in packed even/odd floating-point form.

*/
	static double *ab_mul = 0x0, *cd_mul = 0x0;
	double *a = 0x0;
	double *ivec[2] = {0x0, 0x0};
	int kblocks = n>>10, n_inputs = 0, input;

	struct qfloat qtheta,qr,qi,qn,qt,qc,qs;	/* qfloats used for FFT sincos  calculation. */
	double t1;
#if QFDBG
	double t2;
	long double theta,mt;
	uint64 i1,i2,idiff;
#endif
	static int radix_set, radix_set_save[10] = {1000,0,0,0,0,0,0,0,0,0};
	static int nradices_prim,nradices_radix0,radix_prim[30];/* radix_prim stores sequence of complex FFT radices used, in terms of their prime factors.	*/
	static int *index = 0x0, *index_ptmp = 0x0;		/* Bit-reversal index array and array storing S*I mod N values for DWT weights.	*/

	int i,ii,j,k,l,m,mm,retval;
	int skip_square;	/* Set = 1 to skip dyadic square (e.g. for FFT/IFFT debug) */
	static double n2inv;
	/* roots of unity table pairs needed for FFT: */
	static struct complex *rt0 = 0x0, *rt0_ptmp = 0x0, *rt1 = 0x0, *rt1_ptmp = 0x0;
	double ftmp, re,im, cy_re,cy_im, frac, fracmax;
	char char_str[STR_MAX_LEN];
	static int first_entry=TRUE;		/* Master variable controlling whether the init sequences in subroutines get done (.true.) or not.	*/
#if FFT_DEBUG
	const char dbg_fname[] = "FFT_DEBUG.txt";
#endif

	/* These came about as a result of multithreading, but now are needed whether built unthreaded or multithreaded */
	static int init_sse2 = FALSE;
	int thr_id = -1;	// No multithread support yet.

	ASSERT(HERE, ((uint32)FFT_MUL_BASE >> 16) == 1, "FFT_MUL_BASE != 2^16");

	/***
	Having a separate init block for the big index array allows us to init this prior
	to anything else, using the A-array for scratch space in the call to bit_reverse_int:
	***/
	if(INIT_ARRAYS)
	{
	#if FFT_DEBUG
		/* In FFT-debug mode, need both input arrays to be non-null: */
		ASSERT(HERE, x && y, "if FFT_DEBUG = TRUE, Both x and y-input arrays must be non-null!");
		/* One input to be processed (other used for temporary storage): */
		ivec[n_inputs++] = x;
	#else
		/* In non-debug mode, x-input array needed for temporary storage: */
		ASSERT(HERE, x != 0x0, "if INIT_ARRAYS = TRUE, x-input array must be non-null!");
	#endif
		/* Reset this on an INIT_ARRAYS call to ensure that the
		radix_set != radix_set_save code below gets executed in that case: */
		radix_set_save[0] = 1000;

		/* Clear the a/b/c/d multiplier-data pointers: */
		ab_mul = cd_mul = 0x0;

		/*...******Forward FFT****** sincos index array is here: first, calculate the needed dimension...	*/
		free((void *)index_ptmp); index_ptmp = 0x0; index = 0x0;

		N2 =n/2;		/* Complex vector length.	*/
		n2inv = 1.0/(N2);

		/* Only power-of-2 FFT lengths supported for now: */
		ASSERT(HERE, (n>>trailz32(n)) == 1,"Only power-of-2 FFT lengths supported!");

		/* Use get_fft_radices' zero-index radix set (guaranteed to be available if the FFT length is supported) for now */
		radix_set = 0;
		/* This call sets NRADICES and the first (NRADICES) elements of RADIX_VEC: */
		retval = get_fft_radices(n>>10, radix_set, &NRADICES, RADIX_VEC, 10);

		if(retval == ERR_FFTLENGTH_ILLEGAL)
		{
			sprintf(char_str, "ERROR: length %d = %d K not available.\n", n, n>>10);
			ASSERT(HERE, 0, char_str);
		}
		else if(retval == ERR_RADIXSET_UNAVAILABLE)
		{
			sprintf(char_str, "ERROR: radix set %10d not available.\n",radix_set);
			ASSERT(HERE, 0, char_str);
		}
		else if(retval != 0)
		{
			sprintf(char_str, "ERROR: unknown return value %d from get_fft_radix; N = %d, kblocks = %u, radset = %u.\n", retval, n, kblocks, radix_set);
			ASSERT(HERE, 0, char_str);
		}

		/* My array padding scheme requires N/RADIX_VEC[0] to be a power of 2, and to be >= 2^DAT_BITS, where the latter
		parameter is set in the Mdata.h file: */
		if(n%RADIX_VEC[0] != 0)
		{
			ASSERT(HERE, 0, "FATAL: RADIX_VEC[0] does not divide N!\n");
		}

		/* Make sure n/RADIX_VEC[0] is a power of 2: */
		i = n/RADIX_VEC[0];
		if((i >> trailz32(i)) != 1)
		{
			ASSERT(HERE, 0, "FATAL: n/RADIX_VEC[0] not a power of 2!\n");
		}

		/*...Set the array padding parameters - only use array padding elements for runlengths > 32K. */
		if(kblocks > 32)
		{
			DAT_BITS = DAT_BITS_DEF;
			PAD_BITS = PAD_BITS_DEF;
			/*...If array padding turned on, check that the blocklength divides the unpadded runlength...	*/
			if((DAT_BITS < 31) && ((n >> DAT_BITS) << DAT_BITS) != n)
			{
				ASSERT(HERE, 0,"ERROR: blocklength does not divide runlength!");
			}
		}
		else
		{
			DAT_BITS =31;	/* This causes the padding to go away */
			PAD_BITS = 0;
		}

		if(DAT_BITS < 31)
		{
			/* Now make sure n/RADIX_VEC[0] is sufficiently large (unless n < 2^DAT_BITS, in which case it doesn't matter): */
			if(i < (1 << DAT_BITS))
			{
				sprintf(char_str, "FATAL: n/RADIX_VEC[0] must be >= %u!\n", (1 << DAT_BITS));
				ASSERT(HERE, 0, char_str);
			}
		}

		/* Moved the if(RADIX_SET != radix_set_save) block of code to the if(first_entry)
		init block since otherwise that block won't get executed if we have a new radix set:
		*/

		/* Now init the permuted-index array: */
		k =0;
		mm=RADIX_VEC[0];			/* First radix requires no twiddle factors.	*/

		/* We do the final DIF FFT radix within the wrapper_square routine, so store	*/
		for(i=1; i<NRADICES-1; i++)
		{
			k =k+mm;			/* that block of sincos data there, where they can be merged with the wrapper sincos data.	*/
			mm=mm*RADIX_VEC[i];
		}

		if(mm*RADIX_VEC[NRADICES-1] != N2) {
			ASSERT(HERE, 0, "product of radices not equal to complex vector length\n");
		}

/*		index = (int *)calloc(k,sizeof(int));	*/
		index_ptmp = ALLOC_INT(index_ptmp, k);	if(!index_ptmp){ ASSERT(HERE, 0, "unable to allocate array INDEX in mers_mod_square.\n"); }
		index      = ALIGN_INT(index_ptmp);

		/*...Forward (DIF) FFT sincos data are in bit-reversed order. We define a separate last-pass twiddles
		array within the routine wrapper_square, since that allows us to merge those nicely with the wrapper sincos data.	*/

		k =0;
		l =0;
		mm=1;

		/*...First radix needs no twiddle factors, just need it for building the radix_prim array.	*/

		switch(RADIX_VEC[0]){
		case 8 :
			nradices_radix0 = 3;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 16 :
			nradices_radix0 = 4;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		case 32 :
			nradices_radix0 = 5;
			radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
		default :
			sprintf(char_str, "radix[0] = %d not available.\n",RADIX_VEC[i]);
			ASSERT(HERE, 0, char_str);
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

				bit_reverse_int(&index[k],mm,l,&radix_prim[l-1],-1,(int *)x);

				k += mm;
			}

			/*...All radices beyond the initial-pass one are assumed to be powers of 2:	*/

			switch(RADIX_VEC[i]){
			case 2 :
				radix_prim[l++] = 2; break;
			case 4 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			case 8 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			case 16 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			case 32 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			case 64 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; break;
			case 128 :
				radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2; radix_prim[l++] = 2;
			default :
				sprintf(char_str, "radix %d not available. Halting...\n",RADIX_VEC[i]);
				ASSERT(HERE, 0, char_str);
			}
		}
		nradices_prim = l;

		/* Once we have the FFT radices, call the analogous
		bit-reverse-index-array-init version of the dyadic square function:
		*/
		switch(RADIX_VEC[NRADICES-1])
		{
		  case 16 :
				radix16_pairFFT_mul(x, 0x0, 0x0, n,RADIX_VEC[0],0X0,0X0,0,nradices_prim,radix_prim,0,0, TRUE, FALSE, FALSE); break;
/*
		  case 32 :
				radix32_pairFFT_mul(x, 0x0, 0x0, n,RADIX_VEC[0],0X0,0X0,0,nradices_prim,radix_prim,0,0, TRUE, FALSE, FALSE); break;
*/
		  default :
				sprintf(char_str, "FATAL: radix %d not available for _pairFFT dyadic-mul step.\n",RADIX_VEC[NRADICES-1]);
				ASSERT(HERE, 0, char_str);
		}

	#if !FFT_DEBUG
		return;
	#endif
	}
	else	/* if(INIT_ARRAYS = FALSE) */
	{
	#if FFT_DEBUG
		/* FFT-debug mode not allowed if INIT_ARRAYS = FALSE: */
		ASSERT(HERE, 0, "FFT-debug mode not allowed if INIT_ARRAYS = FALSE!");
	#endif
		/* If FORWARD_FFT_ONLY = TRUE, both input arrays should be non-null - otherwise only X should be valid: */
		ASSERT(HERE, x != 0x0, "X-input null!");
		n_inputs = 1;
		if(FORWARD_FFT_ONLY)
		{
			/* Nullify the a/b/c/d ptrs, in case we're inputting fresh a/b/c/d data to be FFTed: */
			ab_mul = 0x0;
			cd_mul = 0x0;
			/* One or two inputs to be processed? */
			ivec[0] = x;
			ivec[1] = y;
			n_inputs += (y != 0x0);
		}
		else
		{
			ASSERT(HERE,!y, "Non-null Y-input requires FORWARD_FFT_ONLY!");
			/* One input to be processed: */
			ivec[0] = x;
		}
	}	/* end if(INIT_ARRAYS)	*/

	/*...If a new runlength or radix set, deallocate any already-allocated
	allocatable arrays and set first_entry to true:	*/
	for(i = 0; i < 10; i++)
	{
		if(RADIX_VEC[i] != radix_set_save[i])
		{
			first_entry=TRUE;
			break;
		}
	}

/*...on first entry, deallocate any already-allocated allocatable arrays and init DWT weights and FFT sincos tables:	*/
	if(first_entry)
	{
		first_entry=FALSE;

		for(i = 0; i < NRADICES; i++)
		{
			if(RADIX_VEC[i] == 0)
			{
				sprintf(cbuf, "RADIX_VEC[i = %d] zero, for i < [NRADICES = %d]!",i,NRADICES);
				ASSERT(HERE, 0, cbuf);
			}
			radix_set_save[i] = RADIX_VEC[i];
		}
		for(i = NRADICES; i < 10; i++)
		{
			if(RADIX_VEC[i] != 0)
			{
				sprintf(cbuf, "RADIX_VEC[i = %d] nonzero, for i >= [NRADICES = %d]!",i,NRADICES);
				ASSERT(HERE, 0, cbuf);
			}
			radix_set_save[i] = 0;
		}
		/* My array padding scheme requires N/RADIX_VEC[0] to be a power of 2, and to be >= 2^DAT_BITS, where the latter
		parameter is set in the Mdata.h file: */

		if(n%RADIX_VEC[0] != 0)
		{
			sprintf(cbuf  ,"RADIX_VEC[0] does not divide N!\n");
			fprintf(stderr,"%s", cbuf);
			ASSERT(HERE, 0,cbuf);
		}

		/* Make sure n/RADIX_VEC[0] is a power of 2: */
		i = n/RADIX_VEC[0];
		if((i >> trailz32(i)) != 1)
		{
			sprintf(cbuf  ,"n/RADIX_VEC[0] not a power of 2!\n");
			fprintf(stderr,"%s", cbuf);
			ASSERT(HERE, 0,cbuf);
		}

		if(DAT_BITS < 31)
		{
			/* Now make sure n/RADIX_VEC[0] is sufficiently large (unless n < 2^DAT_BITS, in which case it doesn't matter): */
			if(i < (1 << DAT_BITS))
			{
				sprintf(cbuf  ,"vn/RADIX_VEC[0] must be >= %u!\n", (1 << DAT_BITS));
				fprintf(stderr,"%s", cbuf);
				ASSERT(HERE, 0,cbuf);
			}

			/* We also have a lower limit on 2^DAT_BITS set by the wrapper_square routine: */
			if((1 << DAT_BITS) < 2*RADIX_VEC[NRADICES-1])
			{
				sprintf(cbuf  ,"final FFT radix may not exceed = %u!\n", (1 << (DAT_BITS-1)));
				fprintf(stderr,"%s", cbuf);
				ASSERT(HERE, 0,cbuf);
			}
		}
	#if FFT_DEBUG
		sprintf(cbuf,"Using complex FFT radices*");
		char_addr = strstr(cbuf,"*");
		for(i = 0; i < NRADICES; i++)
		{
			sprintf(char_addr,"%10d",RADIX_VEC[i]); char_addr += 10;
		}; sprintf(char_addr,"\n");
		fprintf(stderr,"%s",cbuf);
	#endif

		if(rt0)	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)rt0_ptmp); rt0_ptmp = 0x0; rt0 = 0x0;
			free((void *)rt1_ptmp); rt1_ptmp = 0x0; rt1 = 0x0;
		}

		/**********************************************/
		/* roots of unity table pairs needed for FFT: */
		/**********************************************/

		/* No need for a fancy NINT here: */
		NRT_BITS = (uint32)(log(sqrt(1.0*n))/log(2.0) + 0.5);
		NRT      = 1 << NRT_BITS;
		if(n%NRT) {
			sprintf(cbuf,"FATAL: NRT does not divide N!\n");
			ASSERT(HERE, 0,cbuf);
		}
		NRTM1 = NRT - 1;

		/*...The rt0 array stores the (0:NRT-1)th powers of the [N2]th root of unity
		(i.e. will be accessed using the lower (NRT) bits of the integer sincos index):
		*/
		rt0_ptmp = ALLOC_COMPLEX(rt0_ptmp, NRT);
		if(!rt0_ptmp) {
			sprintf(cbuf,"FATAL: unable to allocate array RT0 in mers_mod_square.\n");
			ASSERT(HERE, 0,cbuf);
		}
		rt0 = ALIGN_COMPLEX(rt0_ptmp);

		qt     = i64_to_q((int64)N2);
		qtheta = qfdiv(Q2PI, qt);	/* 2*pi/(N/2) */
		qr     = qfcos(qtheta);
		qi     = qfsin(qtheta);
		qc     = QONE; qs = QZRO;	/* init sincos multiplier chain. */

		t1 = qfdbl(qr);
	#if QFDBG
		theta = qfdbl(Q2PI)/N2;
		t2 = cos(theta);
		i1 = *(uint64 *)&t1;
		i2 = *(uint64 *)&t2;
		if(i1 != i2)
		{
			idiff = ABS(i1-i2);
			printf("INFO: QCOS1= %16llX, DCOS = %16llX DIFFER BY %20llu\n", i1, i2, idiff);
		}
	#endif

		t1 = qfdbl(qi);
	#if QFDBG
		t2 = sin(theta);
		i1 = *(uint64 *)&t1;
		i2 = *(uint64 *)&t2;
		if(i1 != i2)
		{
			idiff = ABS(i1-i2);
			printf("INFO: QSIN1= %16llX, DSIN = %16llX DIFFER BY %20llu\n", i1, i2, idiff);
		}
	#endif

		qt = QZRO;

		for(i=0; i<NRT; i++)
		{
			qc = qfcos(qt);
			t1 = qfdbl(qc);
		#if QFDBG
			mt = i*theta;
			t2 = cos(mt);
			i1 = *(uint64 *)&t1;
			i2 = *(uint64 *)&t2;
			if(i1 != i2)
			{
				idiff = ABS(i1-i2);
				if(idiff >= (uint64)8)
				{
					printf("INFO: I = %8d: QCOS = %16llX, DCOS = %16llX DIFFER BY %20llu\n", i, i1, i2, idiff);
					/*	ASSERT(HERE, 0,"0");	*/
				}
			}
		#endif
			rt0[i].re = t1;

			qs = qfsin(qt);
			t1 = qfdbl(qs);
		#if QFDBG
			t2 = sin(mt);
			i1 = *(uint64 *)&t1;
			i2 = *(uint64 *)&t2;
			if(i1 != i2)
			{
				idiff = ABS(i1-i2);
				if(idiff >= (uint64)8)
				{
					printf("INFO: I = %8d: QSIN = %16llX, DSIN = %16llX DIFFER BY %20llu\n", i, i1, i2, idiff);
					/*	ASSERT(HERE, 0,"0");	*/
				}
			}
		#endif
			rt0[i].im = t1;
		/*printf("I = %d; RT0 = %20.10f %20.10f\n",i,rt0[i].re,rt0[i].im);	*/
		/*errprint_sincos(&rt0[i].re,&rt0[i].im,(double)(mt));	* Workaround for the DEC Unix 4.0 real*16 sincos bug	*/

			qt = qfadd(qt, qtheta);
			/* Get next terms of the recurrence:  qcnew = qcold*qr - qsold*qi,  qsnew = qcold*qi + qsold*qr. */
			/*qn = qfmul(qc, qr); qt = qfmul(qs, qi); qmul = qfsub(qn, qt);	* Store qcnew in qmul for now. */
			/*qn = qfmul(qc, qi); qt = qfmul(qs, qr); qs   = qfadd(qn, qt); qc = qmul;	*/
		}

		/*...The rt1 array stores the (0:(n/2)/NRT-1)th powers of the [(n/2)/NRT]th root of unity
		(and will be accessed using the upper bits, <NRT:31>, of the integer sincos index):
		*/
		rt1_ptmp = ALLOC_COMPLEX(rt1_ptmp, n/(2*NRT));
		if(!rt1_ptmp) {
			sprintf(cbuf,"FATAL: unable to allocate array RT1 in mers_mod_square.\n");
			ASSERT(HERE, 0,cbuf);
		}
		rt1 = ALIGN_COMPLEX(rt1_ptmp);

		qn     = i64_to_q((int64)NRT);
		qt     = i64_to_q((int64)N2);
		qt     = qfdiv(qn, qt);		/*      NWT/(N/2) */
		qtheta = qfmul(Q2PI, qt);	/* 2*pi*NWT/(N/2) */
		qr     = qfcos(qtheta);
		qi     = qfsin(qtheta);
		qc     = QONE; qs = QZRO;	/* init sincos multiplier chain. */

		t1 = qfdbl(qr);
	#if QFDBG
		theta = qfdbl(Q2PI)*NRT/N2;
		t2 = cos(theta);
		i1 = *(uint64 *)&t1;
		i2 = *(uint64 *)&t2;
		if(i1 != i2)
		{
			idiff = ABS(i1-i2);
			printf("INFO: QCOS2= %16llX, DCOS = %16llX DIFFER BY %20llu\n", i1, i2, idiff);
		}
	#endif

		t1 = qfdbl(qi);
	#if QFDBG
		t2 = sin(theta);
		i1 = *(uint64 *)&t1;
		i2 = *(uint64 *)&t2;
		if(i1 != i2)
		{
			idiff = ABS(i1-i2);
			printf("INFO: QSIN2= %16llX, DSIN = %16llX DIFFER BY %20llu\n", i1, i2, idiff);
		}
	#endif

		qt = QZRO;

		for(i=0; i<(N2/NRT); i++)
		{
			qc = qfcos(qt);
			t1 = qfdbl(qc);
		#if QFDBG
			mt = i*theta;
			t2 = cos(mt);
			i1 = *(uint64 *)&t1;
			i2 = *(uint64 *)&t2;
			if(i1 != i2)
			{
				idiff = ABS(i1-i2);
				if(idiff >= (uint64)8)
				{
					printf("INFO: J = %8d: QCOS = %16llX, DCOS = %16llX DIFFER BY %20llu\n", i, i1, i2, idiff);
					/*	ASSERT(HERE, 0,"0");	*/
				}
			}
		#endif
			rt1[i].re = t1;

			qs = qfsin(qt);
			t1 = qfdbl(qs);
		#if QFDBG
			t2 = sin(mt);
			i1 = *(uint64 *)&t1;
			i2 = *(uint64 *)&t2;
			if(i1 != i2)
			{
				idiff = ABS(i1-i2);
				if(idiff >= (uint64)8)
				{
					printf("INFO: J = %8d: QSIN = %16llX, DSIN = %16llX DIFFER BY %20llu\n", i, i1, i2, idiff);
					/*	ASSERT(HERE, 0,"0");	*/
				}
			}
		#endif
			rt1[i].im = t1;
		/*printf("I = %d; RT1 = %20.10f %20.10f\n",i,rt1[i].re,rt1[i].im);	*/
		/*errprint_sincos(&rt1[i].re,&rt1[i].im,(double)(mt));	* Workaround for the DEC Unix 4.0 real*16 sincos bug	*/

			qt = qfadd(qt, qtheta);
			/* Get next terms of the recurrence:  qcnew = qcold*qr - qsold*qi,  qsnew = qcold*qi + qsold*qr. */
			/*qn = qfmul(qc, qr); qt = qfmul(qs, qi); qmul = qfsub(qn, qt);	* Store qcnew in qmul for now. */
			/*qn = qfmul(qc, qi); qt = qfmul(qs, qr); qs   = qfadd(qn, qt); qc = qmul;	*/
		}

	}	/* end of initialization sequence.	*/

/********************************************************************************************/
/* Do forward-FFT/pointwise-square/inverse-FFT/carry-propagation/fractional-error-checking:	*/
/********************************************************************************************/

	fracmax = 0.0;

/******* DEBUG: do FFT/IFFT on the y-array data, using x-array for temporary storage: *******/
#if FFT_DEBUG
	for(i = 0; i < n; i++)
	{
		j = i + ((i >> DAT_BITS) << PAD_BITS );
		x[j] = y[j];
	}
	skip_square = TRUE;
#else
	skip_square = FALSE;
#endif

	for(input = 0; input < n_inputs; input++)
	{
		/* Data to be processed on current pass will be either x or y-input: */
		a = ivec[input];
	
	#if FFT_DEBUG
		ASSERT(HERE,a==x, "FFT_DEBUG but a != x!");
		ASSERT(HERE, dbg_file == 0x0, "dbg_file != 0x0 prior to fopen");
		dbg_file = fopen(dbg_fname, "a");
		ASSERT(HERE, dbg_file != 0x0, "Unable to open dbg_file!");
		fprintf(dbg_file, "input X-array values:\n");
		for(i = 0; i < n; i+=2)
		{
			j = i + ((i >> DAT_BITS) << PAD_BITS );
			fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j, a[j], a[j+1]);
		}
		fclose(dbg_file); dbg_file = 0x0;
	#endif
	
		/* Perform the initial pass of the forward transform:	*/
	
		switch(RADIX_VEC[0])
		{
		case 8 :
			 radix8_dif_pass1(a,n); break;
		case 16 :
			radix16_dif_pass1(a,n); break;
		case 32 :
			radix32_dif_pass1(a,n); break;
		default :
			sprintf(cbuf,"FATAL: radix %d not available for dif_pass1. Halting...\n",RADIX_VEC[0]);
			ASSERT(HERE, 0,cbuf);
		}
	
	#if FFT_DEBUG
		ASSERT(HERE,a==x, "FFT_DEBUG but a != x!");
		ASSERT(HERE, dbg_file == 0x0, "dbg_file != 0x0 prior to fopen");
		dbg_file = fopen(dbg_fname, "a");
		ASSERT(HERE, dbg_file != 0x0, "Unable to open dbg_file!");
		fprintf(dbg_file, "A-array values after radix%u_dif_pass1, n = %u :\n",RADIX_VEC[0], n);
		for(i = 0; i < n; i+=2)
		{
			j = i + ((i >> DAT_BITS) << PAD_BITS );
			fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j, a[j], a[j+1]);
		}
		fclose(dbg_file); dbg_file = 0x0;
	#endif
	
		/* Break the remaining portion of the FFT into RADIX_VEC[0] blocks: */
	
		for(ii = 0; ii < RADIX_VEC[0]; ++ii)
		{
		#if FFT_DEBUG
			ASSERT(HERE,a==x, "FFT_DEBUG but a != x!");
			dbg_file = fopen(dbg_fname, "a");
			fprintf(dbg_file, "pairFFT_mul_process_chunk with ii = %u:\n",ii);
			fclose(dbg_file); dbg_file = 0x0;
		#endif

			pairFFT_mul_process_chunk(a,ab_mul,cd_mul,n,rt0,rt1,index,ii,nradices_prim,radix_prim, FORWARD_FFT_ONLY, skip_square);
		}
	
		/* In forward-FFT-only mode (assume pairFFT_mul_process_chunk has done none of the IFFT passes),
		save the ptr to the just-FFTed vector in appropriate one of thr a/b or c/d static pointers and
		either cycle or return: */
		if(FORWARD_FFT_ONLY)
		{
			if(input == 0)
				ab_mul = a;	// ab_mul = ivec[0]
			else
				cd_mul = a;	// cd_mul = ivec[1]
			continue;
		}
	
	/*...Do the final inverse FFT pass, carry propagation and initial forward FFT pass in one fell swoop, er, swell loop...	*/
	
		switch(RADIX_VEC[0])
		{
		  case  8 :
			 radix8_dit_pass1(a,n);	break;
		  case 16 :
			radix16_dit_pass1(a,n);	break;
		  case 32 :
			radix32_dit_pass1(a,n);	break;
		  default :
			sprintf(char_str, "radix %d not available for final IFFT pass!\n",RADIX_VEC[0]);
			ASSERT(HERE, 0, char_str);
		}
	
	#if FFT_DEBUG
		ASSERT(HERE,a==x, "FFT_DEBUG but a != x!");
		dbg_file = fopen(dbg_fname, "a");
		fprintf(dbg_file, "A-array values after final radix%u_dit_pass :\n",RADIX_VEC[0]);
		for(i = 0; i < n; i+=2)
		{
			j = i + ((i >> DAT_BITS) << PAD_BITS );
			fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j, a[j], a[j+1]);
		}
		fclose(dbg_file); dbg_file = 0x0;
	#endif

	}	/* endfor(input = 0; input < n_inputs; input++) */

	if(FORWARD_FFT_ONLY)
		return;

	/************* If FFT_DEBUG = TRUE, do a simple FFT/IFFT-returns-original-inputs test: **************/
	// Carry step proceeds in separate even(real)/odd(imag) fashion:
	cy_re = cy_im = 0;
	for(i = 0; i < n; i+=2)
	{
		j = i + ((i >> DAT_BITS) << PAD_BITS );
	// Real part:
		re = a[j]*n2inv + cy_re;
		ftmp = DNINT(re);
		frac = fabs(re-ftmp);
		if(frac > fracmax)
			fracmax = frac;
	#if FFT_DEBUG
		if(ftmp != y[j])
		{
			sprintf(char_str, "%20.5f != %20.5f in FFT_DEBUG post-test loop!\n",ftmp, y[j]);
			ASSERT(HERE, 0, char_str);
		}
		a[j] = ftmp;
	#else
		cy_re = DNINT(ftmp*FFT_MUL_BASE_INV);
		a[j] = ftmp - cy_re*FFT_MUL_BASE;
	#endif
		j++;
	// Imag part:
		im = a[j]*n2inv + cy_im;
		ftmp = DNINT(im);
		frac = fabs(im-ftmp);
		if(frac > fracmax)
			fracmax = frac;
	#if FFT_DEBUG
		if(ftmp != y[j])
		{
			sprintf(char_str, "%20.5f != %20.5f in FFT_DEBUG post-test loop!\n",ftmp, y[j]);
			ASSERT(HERE, 0, char_str);
		}
		a[j] = ftmp;
	#else
		cy_im = DNINT(ftmp*FFT_MUL_BASE_INV);
		a[j] = ftmp - cy_im*FFT_MUL_BASE;
	#endif
	}
	ASSERT(HERE, fabs(cy_re) + fabs(cy_im) == 0, " Fatal: Nonzero exit carry!");

#if FFT_DEBUG
printf("FFT_DEBUG: FFT/IFFT-returns-original-inputs test PASS\n",fracmax);
#endif
printf("fracmax = %15.10f\n",fracmax);

	/*...Now do the fractional error check. */
	if(fracmax >= 0.4)
	{
		sprintf(char_str, "Roundoff warning, maxerr = %16.12f\n", fracmax);
		WARN(HERE, char_str, "", 0);

		/*...Fractional parts close to 0.5 cause the program to quit.
		We put the interactive-mode errlimit close to 0.5, to let people really push the limits if they want to...	*/
		if(fracmax > 0.47 )
		{
			ASSERT(HERE, 0, " Fatal Roundoff error...Halting execution.");
		}
	}

	return;
}

/***************/

void pairFFT_mul_process_chunk(double a[], double ab_mul[], double cd_mul[], int n, struct complex rt0[], struct complex rt1[], int index[], int ii, int nradices_prim, int radix_prim[], int FORWARD_FFT_ONLY, int skip_square)
{
	int thr_id = 0;	/* In unthreaded mode this must always = 0 */
	int init_sse2 = FALSE;	/* Init-calls to the various radix-pass routines presumed done prior to entry into this routine */

	int i,incr,istart,jstart,k,koffset,l,mm;
	double *ab_mul_ptr,*cd_mul_ptr;
	char char_str[STR_MAX_LEN];
#if (FFT_DEBUG && 0)
	const char dbg_fname[] = "FFT_DEBUG.txt";
#endif

	l = ii;
	k    = 0;
	mm   = 1;
	incr = n/RADIX_VEC[0];

	istart = l*incr;	/* Starting location of current data-block-to-be-processed within A-array. */
	jstart = istart + ((istart >> DAT_BITS) << PAD_BITS );

	if(!ab_mul)
	{
		ASSERT(HERE, !cd_mul, "pairFFT_mul_process_chunk: AB_MUL = 0x0 but CD_MUL non-null!");
	#if !FFT_DEBUG
		ASSERT(HERE, FORWARD_FFT_ONLY == TRUE, "pairFFT_mul_process_chunk: AB_MUL = 0x0 but FORWARD_FFT_ONLY not set!");
	#endif
		ab_mul_ptr = 0x0;
		cd_mul_ptr = 0x0;
	}
	else if(!cd_mul)
	{
		ab_mul_ptr = &ab_mul[jstart];
		cd_mul_ptr = 0x0;
	}
	else
	{
		ab_mul_ptr = &ab_mul[jstart];
		cd_mul_ptr = &cd_mul[jstart];
	}

	for(i=1; i <= NRADICES-2; i++)
	{
		/* Offset from base address of index array = L*NLOOPS = L*MM : */
		koffset = l*mm;

#if (FFT_DEBUG && 0)
	ASSERT(HERE, dbg_file == 0x0, "dbg_file != 0x0 prior to fopen");
	dbg_file = fopen(dbg_fname, "a");
	ASSERT(HERE, dbg_file != 0x0, "Unable to open dbg_file!");
	fprintf(dbg_file, "pairFFT_mul_process_chunk: i = %u, calling radix%u_dif_pass with jstart = %u, koffset = %u, mm = %u, incr = %u:\n",i,RADIX_VEC[0],jstart,koffset,mm,incr);
	fclose(dbg_file); dbg_file = 0x0;
#endif
		switch(RADIX_VEC[i])
		{
		case  8 :
			 radix8_dif_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		case 16 :
			radix16_dif_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		case 32 :
			radix32_dif_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		default :
			sprintf(cbuf,"pairFFT_mul_process_chunk: FATAL: radix %d not available for dif_pass. Halting...\n",RADIX_VEC[i]);
			ASSERT(HERE, 0,cbuf);
		}

		k    += mm*RADIX_VEC[0];
		mm   *= RADIX_VEC[i];
		incr /= RADIX_VEC[i];

	}	/* end i-loop. */

/*...Final DIF pass, dyadic squaring and initial DIT pass are all done in-place. */
	koffset = l*mm;

	switch(RADIX_VEC[NRADICES-1])
	{
		/* Call radix*_pairFFT_mul with INIT_ARRAYS = 0, FORWARD_FFT_ONLY and skip_square as passed: */
	  case 16 :
			radix16_pairFFT_mul(&a[jstart],ab_mul_ptr,cd_mul_ptr,n,RADIX_VEC[0],rt0,rt1,ii,nradices_prim,radix_prim,mm,incr, FALSE, FORWARD_FFT_ONLY, skip_square); break;
/*
	  case 32 :
			radix32_pairFFT_mul(&a[jstart],ab_mul_ptr,cd_mul_ptr,n,RADIX_VEC[0],rt0,rt1,ii,nradices_prim,radix_prim,mm,incr, FALSE, FORWARD_FFT_ONLY, skip_square); break;
*/
	  default :
			sprintf(char_str, "pairFFT_mul_process_chunk: FATAL: radix %d not available for dyadic mul step.\n",RADIX_VEC[NRADICES-1]);
			ASSERT(HERE, 0, char_str);
	}

	/* In forward-FFT-only mode, do none of the IFFT passes: */
	if(FORWARD_FFT_ONLY)
		return;

/*...Rest of inverse decimation-in-time (DIT) transform. Note that during IFFT we process the radices in reverse
	order. The first array sent to each pass routine is assumed to contain the bit-reversed floating data.	*/

	/* Quick-n-dirty way of generating the correct starting values of k, mm and incr -
	simply use the skeleton of the forward (DIF) loop, sans the i = NRADICES-2 pass
	(since k, mm and incr are post-incremented):
	*/
	k    = 0;
	mm   = 1;
	incr = n/RADIX_VEC[0];

	for(i=1; i <= NRADICES-2; i++)
	{
		k    += mm*RADIX_VEC[0];
		mm   *= RADIX_VEC[i];
		incr /= RADIX_VEC[i];
	}

	/* Now do the DIT loop, running the radices (and hence the values of k, mm and incr) in reverse: */

	for(i=NRADICES-2; i >= 1; i--)
	{
		incr *= RADIX_VEC[i];
		mm   /= RADIX_VEC[i];
		k    -= mm*RADIX_VEC[0];

		koffset = l*mm;

#if (FFT_DEBUG && 0)
	dbg_file = fopen(dbg_fname, "a");
	fprintf(dbg_file, "pairFFT_mul_process_chunk: i = %u, calling radix%u_dit_pass with jstart = %u, koffset = %u, mm = %u, incr = %u:\n",i,RADIX_VEC[0],jstart,koffset,mm,incr);
	fclose(dbg_file); dbg_file = 0x0;
#endif
		switch(RADIX_VEC[i])
		{
		case  8 :
			 radix8_dit_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		case 16 :
			radix16_dit_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		case 32 :
			radix32_dit_pass(&a[jstart],n,rt0,rt1,&index[k+koffset],mm,incr,init_sse2,thr_id); break;
		default :
			sprintf(cbuf,"pairFFT_mul_process_chunk: FATAL: radix %d not available for dit_pass. Halting...\n",RADIX_VEC[i]);
			ASSERT(HERE, 0,cbuf);
		}
	}	/* end i-loop */
}

