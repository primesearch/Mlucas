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

/***************/

/*
If INIT_ARRAYS = TRUE, must have valid length-n x-array pointer, which is used for
scratch space for the init and whose contents are thus assumed to be destroyed!
*/
void pairFFT_mul(double x[], double y[], double z[], int n, int INIT_ARRAYS, int FORWARD_FFT_ONLY)
{
	const char func[] = "pairFFT_mul";
/*
	Depending on the precise way the function is called, either computes the forward FFT
	of a length-n paired-real-vector FFT (i.e. length-n/2 complex FFT of a pair of
	length-n/2 even/odd-interleaved real input vectors) or uses a pointer to 1 or 2 such
	precomputed fFFT vectors to perform either a pair of in-place multiprecision (MP)
	integer linear combinations of the form

		( u' )   ( a  -b ) ( u )
		(    ) = (       )*(   )		(*)
		( v' )   ( c  -d ) ( v )

	or a pair of multiprecision real-vector-times-real-vector multiplies of form

		( u' )   ( a   0 ) ( u )
		(    ) = (       )*(   )		(**)
		( v' )   ( 0   b ) ( v )

	All FFTs are done in-place. This function can be called in either init or active-data mode,
	based on the boolean value of the input INIT_ARRAYS:

	If INIT_ARRAYS = TRUE, the function simply uses the input (n) to init the FFT sincos
	data tables and returns; the x-array is assumed to be valid but uninitialized, and is
	used for scratch space for the init procedures. All other input arguments are ignored.
		Such an init call MUST PRECEDE any calls in which actual FFT computation occurs,
	the latter being effected by calling with INIT_ARRAYS = FALSE.

	INIT_ARRAYS = FALSE:

	In the above formula (*), a, b, c, d, u, v are multiword integers input in
	extracted-to-double form, that is, converted to base-W balanced-digit floating-point
	form with some wordsize W allowing for direct length-n FFT-based multiply. (W = 2^16
	is a convenient and safe wordsize for most applications). It is assumed that a given
	set of matrix coefficients a/b/c/d in (*) will generally be applied to multiple sets
	of u/v data. Thus it is desirable to precompute and store the forward-transforms of
	a/b/c/d just once and then reuse same as many times as needed. For that purpose
	(and to minimize code duplication) we provide 2 active-data modes (as opposed to
	precomputed-data initialization mode) of calling this function, based on the boolean
	value of the input FORWARD_FFT_ONLY:

	FORWARD_FFT_ONLY = TRUE: ***** Note that by 'TRUE' here we specifically mean = 1 *****
		The a/b and c/d multiplier data are assumed to be input in
	packed even/odd form in the x and y-arrays, respectively - if both array pointers are
	valid we will compute the fFFTs of both input vectors, otherwise if only x != 0 we
	compute the fFFT of it only, in-place. Such precompute-fFFT calls must precede
	any call to compute (*) or (**), i.e. a call with FORWARD_FFT_ONLY = FALSE
	and with the u/v data input in the x-array.

	FORWARD_FFT_ONLY = 2: This is just like FORWARD_FFT_ONLY = FALSE, but with x-input also
					assumed to be already fFFTed on entry.

	FORWARD_FFT_ONLY = FALSE: The u/v multiplicands are assumed to be input in
	packed even/odd form in the x-array (y is *required* to be null), and at least the y-array
	input pointer must be non-null and point to precomputed fFFT data. The forward FFT(u/v) is
	computed in-place, and the appropriate linear combination (*) or (**) is computed
	depending on the y and z-pointers (y must be be non-null in this scenario):

	[1]	y != 0, z = 0: The presumed-to-already-be-fwd-FFTed data in y take the place of a/b in (*)
					and the result is obtained by combining FFT(u/v) with FFT(a/b),
					followed by in-place iFFT in packed even/odd floating-point form;

	[2] y,z both != 0: The presumed-to-already-be-fwd-FFTed data in y,z take the place of a/b and
					c,d in (**) and the result is obtained by combining FFT(u/v) with FFT(a/b, c/d),
					followed by in-place iFFT in packed even/odd floating-point form.

*/
	double *ab_mul = 0x0, *cd_mul = 0x0;
	double *a = 0x0;
	double *ivec[2] = {0x0, 0x0};
	int kblocks = n>>10, n_inputs = 0, input;

	struct qfloat qtheta,qr,qi,qn,qt,qc,qs;	/* qfloats used for FFT sincos  calculation. */
	static int radix_set, nradsets, radix_set_save[10] = {1000,0,0,0,0,0,0,0,0,0};
	static int radix_vec0, nchunks, modtype_save; 	// Store frequently-used RADIX_VEC[0] and modulus type on entry
	static int nradices_prim,nradices_radix0,radix_prim[30];/* radix_prim stores sequence of complex FFT radices used, in terms of their prime factors.	*/
	static int *index = 0x0, *index_ptmp = 0x0;		/* Bit-reversal index array and array storing S*I mod N values for DWT weights.	*/
	static int *block_index;				/* array storing the RADIX_VEC[0] data-block indices for pass-2 of the FFT.	*/
	/* arrays storing the index values needed for the paired-block wrapper/square scheme: */
	static int *ws_i,*ws_j1,*ws_j2,*ws_j2_start,*ws_k,*ws_m,*ws_blocklen,*ws_blocklen_sum;
	int i,ii,ierr,j,jhi,k,l,m,mm,k2,m2,l1,l2,l2_start,blocklen,blocklen_sum,outer,retval;
	static double n2inv,radix_inv;
	/* roots of unity table pairs needed for FFT: */
	static struct complex *rt0 = 0x0, *rt0_ptmp = 0x0, *rt1 = 0x0, *rt1_ptmp = 0x0;
	double ftmp, re,im, cy_re,cy_im, frac, fracmax;
	char char_str[STR_MAX_LEN];
	static int first_entry=TRUE;		/* Master variable controlling whether the init sequences in subroutines get done (.true.) or not.	*/
	char *char_addr;

	/* These came about as a result of multithreading, but now are needed whether built unthreaded or multithreaded */
	static int init_sse2 = FALSE;
	int thr_id = -1;	// No multithread support yet.

	ASSERT(((uint32)FFT_MUL_BASE >> 16) == 1, "FFT_MUL_BASE != 2^16");

	/***
	Having a separate init block for the big index array allows us to init this prior
	to anything else, using the A-array for scratch space in the call to bit_reverse_int:
	***/
	if(INIT_ARRAYS)
	{
		/* In init mode, x-input array used for temporary storage: */
		ASSERT(x != 0x0, "if INIT_ARRAYS = TRUE, x-input array must be non-null!");

		/* Reset this on an INIT_ARRAYS call to ensure that the
		radix_set != radix_set_save code below gets executed in that case: */
		radix_set_save[0] = 1000;
	#warning Add support for multiple FFT lengths, each with associated saved data in order to obviate need for redundant re-initing!

		/*...******Forward FFT****** sincos index array is here: first, calculate the needed dimension...	*/
		free((void *)index_ptmp); index_ptmp = 0x0; index = 0x0;

		N2 =n/2;		/* Complex vector length.	*/
		n2inv = 1.0/(N2);

		/* Only power-of-2 FFT lengths supported for now: */
		ASSERT((n>>trailz32(n)) == 1,"Only power-of-2 FFT lengths supported!");

		// Use get_fft_radices' zero-index radix set (guaranteed to be available if the FFT length is supported)
		// to find how many different radsets available at this length, then loop over them (including the 0-one)
		// in order to find the first one satisfying the radix-set constraints of this routine [carry-radix = 16
		// or 32, dyad-mul-surrounding radix = 16]:
		radix_set = -1;	retval = get_fft_radices(n>>10, radix_set, &nradsets, RADIX_VEC, 10);
		for(radix_set = 0; radix_set < nradsets; radix_set++) {
			/* This call sets NRADICES and the first (NRADICES) elements of RADIX_VEC: */
			retval = get_fft_radices(n>>10, radix_set, &NRADICES, RADIX_VEC, 10);
			if(retval == ERR_FFTLENGTH_ILLEGAL) {
				sprintf(char_str, "ERROR: length %d = %d K not available.\n", n, n>>10);
				ASSERT(0, char_str);
			} else if(retval == ERR_RADIXSET_UNAVAILABLE) {
				sprintf(char_str, "ERROR: radix set %10d not available.\n",radix_set);
				ASSERT(0, char_str);
			} else if(retval != 0) {
				sprintf(char_str, "ERROR: unknown return value %d from get_fft_radix; N = %d, kblocks = %u, radset = %u.\n", retval, n, kblocks, radix_set);
				ASSERT(0, char_str);
			}
			// Make sure n/radix_vec0 >= 1024:
			if(n/RADIX_VEC[0] < 1024)
				continue;
			if( (RADIX_VEC[NRADICES-1] == 16) && (RADIX_VEC[0] == 8 || RADIX_VEC[0] == 16 || RADIX_VEC[0] == 32) )
				break;
		}
		ASSERT(radix_set < nradsets, "Unable to find suitable radix set!");
		radix_vec0 = RADIX_VEC[0];
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)radix_vec0));
		nchunks = radix_vec0>>1;

		/* My array padding scheme requires N/radix_vec0 to be a power of 2, and to be >= 2^DAT_BITS, where the latter
		parameter is set in the Mdata.h file: */
		if(n%radix_vec0 != 0) {
			ASSERT(0, "ERROR: RADIX_VEC[0] does not divide N!\n");
		}

		/* Make sure n/radix_vec0 is a power of 2: */
		i = n/radix_vec0;
		if((i >> trailz32(i)) != 1) {
			ASSERT(0, "ERROR: n/RADIX_VEC[0] not a power of 2!\n");
		}

		/*...Set the array padding parameters - only use array padding elements for runlengths > 32K. */
		if(DAT_BITS < 31) {
			/*...If array padding turned on, check that the blocklength divides the unpadded runlength...	*/
			ASSERT(((n >> DAT_BITS) << DAT_BITS) == n,"ERROR: blocklength does not divide runlength!");

			/* Now make sure n/RADIX_VEC[0] is sufficiently large (unless n < 2^DAT_BITS, in which case it doesn't matter): */
			if(i < (1 << DAT_BITS)) {
				sprintf(char_str, "ERROR: n/RADIX_VEC[0] must be >= %u!\n", (1 << DAT_BITS));
				ASSERT(0, char_str);
			}
		}

		/* Moved the if(RADIX_SET != radix_set_save) block of code to the if(first_entry)
		init block since otherwise that block won't get executed if we have a new radix set:
		*/

		/* Now init the permuted-index array: */
		k = 0;
		mm = radix_vec0;			/* First radix requires no twiddle factors.	*/

		/* We do the final DIF FFT radix within the pairFFT_mul routine, so store
		that block of sincos data there, where they can be merged with the wrapper sincos data:
		*/
		for(i = 1; i < NRADICES-1; i++) {
			k = k+mm;
			mm = mm*RADIX_VEC[i];
		}

		if(mm*RADIX_VEC[NRADICES-1] != N2) {
			ASSERT(0, "product of radices not equal to complex vector length\n");
		}

/*		index = (int *)calloc(k,sizeof(int));	*/
		index_ptmp = ALLOC_INT(index_ptmp, k);	if(!index_ptmp){ ASSERT(0, "unable to allocate array INDEX in pairFFT_mul.\n"); }
		index      = ALIGN_INT(index_ptmp);

		/*...Forward (DIF) FFT sincos data are in bit-reversed order. We define a separate last-pass twiddles
		array within the routine pairFFT_mul, since that allows us to merge those nicely with the wrapper sincos data.	*/

		k = 0;	l = 0;	mm = 1;

		/*...First radix needs no twiddle factors, just need it for building the radix_prim array.	*/
		switch(radix_vec0){
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
			ASSERT(0, char_str);
		}

		for(i=1; i < NRADICES; i++)
		{
			/*...Allocate and initialize an index array containing MM indices...	*/
			if(i<(NRADICES-1)) {
				mm=mm*RADIX_VEC[i-1];	/* MM = product of all the preceding radices	*/
				for(m=0; m < mm; m++) {
					index[k+m]=m;
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
				bit_reverse_int(&index[k],mm,l,&radix_prim[l-1],-1,(int *)x);
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
				sprintf(char_str, "radix %d not available. Halting...\n",RADIX_VEC[i]);
				ASSERT(0, char_str);
			}
		}
		nradices_prim = l;

		/* Once we have the FFT radices, call the analogous
		bit-reverse-index-array-init version of the dyadic square function:
		*/
		switch(RADIX_VEC[NRADICES-1])
		{
		  case 16 :
			radix16_pairFFT_mul(x, 0x0, 0x0, n,radix_vec0,0X0,0X0, nradices_prim,radix_prim, 0,0,0,0,0,0,0,0, TRUE, FALSE); break;
/*
		  case 32 :
			radix32_pairFFT_mul(x, 0x0, 0x0, n,radix_vec0,0X0,0X0,0,nradices_prim,radix_prim,0,0, TRUE, FALSE, FALSE); break;
*/
		  default :
			sprintf(char_str, "ERROR: radix %d not available for _pairFFT dyadic-mul step.\n",RADIX_VEC[NRADICES-1]);
			ASSERT(0, char_str);
		}

		return;
	}
	else	/* if(INIT_ARRAYS = FALSE) */
	{
		/* If FORWARD_FFT_ONLY = TRUE, at least the X-ptr should be valid: */
		n_inputs = 1;
		if((uint32)FORWARD_FFT_ONLY > 2) {
			ASSERT(0, "FORWARD_FFT_ONLY not a any-nonzero-denotes-TRUE param: legal TRUE-values are 1 and 2!");
		} else if(FORWARD_FFT_ONLY == 1) {
			ASSERT(x != 0x0 && z == 0x0, "FORWARD_FFT_ONLY requires X-input nonzero and Z-input null!");
			/* One or two inputs to be processed? */
			ivec[0] = x;
			ivec[1] = y;
			n_inputs += (y != 0x0);
		} else {	// FORWARD_FFT_ONLY = 0 and 2 behave similarly
			ASSERT(x != 0x0 && y != 0x0, "FORWARD_FFT_ONLY = FALSE requires Non-null X,Y-inputs!");
			/* One input to be processed: */
			ivec[0] = x;
			ab_mul = y; cd_mul = z;
		}
	}	/* end if(INIT_ARRAYS)	*/

	// Since will typically be calling this routine to do GCD in between some specific-modulus work,
	// save the transform type of the latter and temporarily switch to the generic-FFT-mul type for the carry-step call:
	modtype_save = MODULUS_TYPE;
	MODULUS_TYPE = MODULUS_TYPE_GENFFTMUL;

	/*...If a new runlength or radix set, deallocate any already-allocated
	allocatable arrays and set first_entry to true:	*/
	for(i = 0; i < 10; i++) {
		if(RADIX_VEC[i] != radix_set_save[i]) {
			first_entry=TRUE;
			break;
		}
	}

/*...on first entry, deallocate any already-allocated allocatable arrays and init DWT weights and FFT sincos tables:	*/
	if(first_entry) {
		first_entry=FALSE;

		for(i = 0; i < NRADICES; i++) {
			if(RADIX_VEC[i] == 0) {
				sprintf(cbuf, "RADIX_VEC[i = %d] zero, for i < [NRADICES = %d]!",i,NRADICES);
				ASSERT(0, cbuf);
			}
			radix_set_save[i] = RADIX_VEC[i];
		}
		for(i = NRADICES; i < 10; i++) {
			if(RADIX_VEC[i] != 0) {
				sprintf(cbuf, "RADIX_VEC[i = %d] nonzero, for i >= [NRADICES = %d]!",i,NRADICES);
				ASSERT(0, cbuf);
			}
			radix_set_save[i] = 0;
		}

		/* My array padding scheme requires N/RADIX_VEC[0] to be a power of 2, and to be >= 2^DAT_BITS, where the latter
		parameter is set in the Mdata.h file: */
		if(n%radix_vec0 != 0) {
			sprintf(cbuf  ,"RADIX_VEC[0] does not divide N!\n");
			fprintf(stderr,"%s", cbuf);
			ASSERT(0,cbuf);
		}

		/* Make sure n/RADIX_VEC[0] is a power of 2: */
		i = n/radix_vec0;
		if((i >> trailz32(i)) != 1) {
			sprintf(cbuf  ,"n/RADIX_VEC[0] not a power of 2!\n");
			fprintf(stderr,"%s", cbuf);
			ASSERT(0,cbuf);
		}

		if(DAT_BITS < 31) {
			/* Now make sure n/RADIX_VEC[0] is sufficiently large (unless n < 2^DAT_BITS, in which case it doesn't matter): */
			if(i < (1 << DAT_BITS)) {
				sprintf(cbuf  ,"vn/RADIX_VEC[0] must be >= %u!\n", (1 << DAT_BITS));
				fprintf(stderr,"%s", cbuf);
				ASSERT(0,cbuf);
			}

			/* We also have a lower limit on 2^DAT_BITS set by the pairFFT_mul routine: */
			if((1 << DAT_BITS) < 2*RADIX_VEC[NRADICES-1]) {
				sprintf(cbuf  ,"final FFT radix may not exceed = %u!\n", (1 << (DAT_BITS-1)));
				fprintf(stderr,"%s", cbuf);
				ASSERT(0,cbuf);
			}
		}

		sprintf(cbuf,"%s: Using complex FFT radices*",func);
		char_addr = strstr(cbuf,"*");
		for(i = 0; i < NRADICES; i++) {
			sprintf(char_addr,"%10d",RADIX_VEC[i]); char_addr += 10;
		}; sprintf(char_addr,"\n");
		fprintf(stderr,"%s",cbuf);

		if(rt0) {	/* If it's a new exponent of a range test, need to deallocate these. */
			free((void *)rt0_ptmp); rt0_ptmp = 0x0; rt0 = 0x0;
			free((void *)rt1_ptmp); rt1_ptmp = 0x0; rt1 = 0x0;
			free((void *)ws_i           ); ws_i            = 0x0;
			free((void *)ws_j1          ); ws_j1           = 0x0;
			free((void *)ws_j2          ); ws_j2           = 0x0;
			free((void *)ws_j2_start    ); ws_j2_start     = 0x0;
			free((void *)ws_k           ); ws_k            = 0x0;
			free((void *)ws_m           ); ws_m            = 0x0;
			free((void *)ws_blocklen    ); ws_blocklen     = 0x0;
			free((void *)ws_blocklen_sum); ws_blocklen_sum = 0x0;
		}

		/**********************************************/
		/* Roots of unity table pairs needed for FFT: */
		/**********************************************/

		/* No need for a fancy NINT here: */
		NRT_BITS = (uint32)(log(sqrt(1.0*n))/log(2.0) + 0.5);	NRT = 1 << NRT_BITS;	NRTM1 = NRT - 1;
		if(n%NRT) {
			sprintf(cbuf,"ERROR: NRT does not divide N!\n");
			ASSERT(0,cbuf);
		}

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
		qc = QONE; qs = qt = QZRO;	/* init sincos multiplier chain. */
		for(i = 0; i < NRT; i++) {
			qc = qfcos(qt);	rt0[i].re = qfdbl(qc);
			qs = qfsin(qt);	rt0[i].im = qfdbl(qs);
			qt = qfadd(qt, qtheta);
			/* Get next terms of the recurrence:  qcnew = qcold*qr - qsold*qi,  qsnew = qcold*qi + qsold*qr. */
			/*qn = qfmul(qc, qr); qt = qfmul(qs, qi); qmul = qfsub(qn, qt);	* Store qcnew in qmul for now. */
			/*qn = qfmul(qc, qi); qt = qfmul(qs, qr); qs   = qfadd(qn, qt); qc = qmul;	*/
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
		qc     = QONE; qs = qt = QZRO;	/* init sincos multiplier chain. */
		for(i=0; i<(N2/NRT); i++) {
			qc = qfcos(qt);	rt1[i].re = qfdbl(qc);
			qs = qfsin(qt);	rt1[i].im = qfdbl(qs);
			qt = qfadd(qt, qtheta);
			/* Get next terms of the recurrence:  qcnew = qcold*qr - qsold*qi,  qsnew = qcold*qi + qsold*qr. */
			/*qn = qfmul(qc, qr); qt = qfmul(qs, qi); qmul = qfsub(qn, qt);	* Store qcnew in qmul for now. */
			/*qn = qfmul(qc, qi); qt = qfmul(qs, qr); qs   = qfadd(qn, qt); qc = qmul;	*/
		}

		/* Break the remaining portion of the FFT into radix_vec0 blocks, each of which ideally
			 should operate on a dataset which fits entirely into the L2 cache of the host machine. */

		/* 8/23/2004: Need to allocate an extra element here to account for the padding element that gets inserted when radix_vec0 is odd: */

		block_index = (int *)calloc((radix_vec0+1),sizeof(int));
		if(!block_index){ sprintf(cbuf,"ERROR: unable to allocate array BLOCK_INDEX in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		/*
		Examples - We only allow powers of 2 here, for the more general case cf. mers_mod_square.c:

			radix0 = 4: want 2 sets of l1,2 pairs = {0,1},{2,3}
			radix0 = 8: want 3 sets of l1,2 pairs = {0,1},{2,3},{4,7,5,6}
			radix0 =16: want 4 sets of l1,2 pairs = {0,1},{2,3},{4,7,5,6},{8,15,9,14,10,13,11,12}, loop control
											has blocklen = 2,2,4,8	l2_start = 1,3,7,15	increments by 2,4,8

		Rules:
			1) For even radix0 use [nradices_radix0] sets of blocks; for odd radices use [nradices_radix0+1] sets,
				the first of which consists of a single block.
			2) Outer loop over prime factors of radix0, in reverse order of which these subradices appear in forward FFT:
				- blocklen (number of l-values of) first set = 2 - radix0%2
				- blocklen of sets 2...nradices_radix0 = (number of blocks done in previous sets) * (current subradix - 1) .
			3) Throughout block processing, l1 increases monotonically from 0. Within each block, for each value of l1, the
				corresponding l2_start starts at 1 - radix0%2and increases by blocklen_newblock for each new block.

		In a multithreaded implementation, each thread does 2 of the above blocks. E.g. if NTHREADS = 4 and radix0 = 16, the
		chunk-processing loop occurs in 2 parallel passes, with each of the 4 threads processing the following blockpairs:

			Pass 1:
			Thread 0: ( 0, 1)
			Thread 1: ( 2, 3)
			Thread 2: ( 4, 7)
			Thread 3: ( 5, 6)

			Pass 2:
			Thread 0: ( 8,15)
			Thread 1: ( 9,14)
			Thread 2: (10,13)
			Thread 3: (11,12)

		If NTHREADS does not divide radix0/2, there will always be one or more under-or-unutilized threads.
		*/
		blocklen = 2;			// blocklen = 2 for even radix_vec0
		blocklen_sum = 0;
		l1 = 0;
		l2 = 1; l2_start = l2;	// l2_start = 1 for even radix_vec0

		// Init the block_index array, which allows a linear index to provide access into the nonmonotone block_index array:
		ii = 0;
		for(outer = nradices_radix0 - l2; outer >= 0; outer-- )	/* This mimics the increasing-blocksize loop in the pairFFT_mul routines. */
		{
		//	fprintf(stderr,"Block %d : blocklen = %d  blocklen_sum = %d\n",nradices_radix0 - outer,blocklen,blocklen_sum);

			// Process 2 blocks per loop execution, thus only execute the loop half as many times as otherwise:
			for(m = 0; m <= (blocklen-1)>>1; m++)
			{
				// Now execute body of loop once with l = l1, once with l = l2:
				l = l1;
				// Do two loop executions:
				for(j = 0; j < 2; j++)
				{
					if(!(l >= 0 && l < radix_vec0)) { sprintf(cbuf,"ERROR 10 in %s.c\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
					block_index[ii] = l;	//fprintf(stderr,"%3d %3d\n",ii,l);
					ii++;	// every time we execute this innermost loop (which corresponds to one
							// block of FFT data being processed), increment the linear array index
					l += l2-l1;
				}	/* end j-loop */
				l1++;
				l2--;
			}
			if(outer == 0) break;

			blocklen_sum += blocklen;
			blocklen = (radix_prim[outer-1] - 1)*blocklen_sum;

			// Next j2_start is previous one plus the length of the current block:
			l1 = l2_start + 1;
			l2_start += blocklen;
			l2 = l2_start;	// Reset j2 for start of the next block.
		}		/* End of Main loop */

		/* arrays storing the index values needed for the parallel-block wrapper/square scheme: */
		if( !(ws_i            = (int *)calloc(radix_vec0,sizeof(int))) ) { sprintf(cbuf,"ERROR: unable to allocate array WS_I            in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		if( !(ws_j1           = (int *)calloc(radix_vec0,sizeof(int))) ) { sprintf(cbuf,"ERROR: unable to allocate array WS_J1           in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		if( !(ws_j2           = (int *)calloc(radix_vec0,sizeof(int))) ) { sprintf(cbuf,"ERROR: unable to allocate array WS_J2           in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		if( !(ws_j2_start     = (int *)calloc(radix_vec0,sizeof(int))) ) { sprintf(cbuf,"ERROR: unable to allocate array WS_J2_START     in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		if( !(ws_k            = (int *)calloc(radix_vec0,sizeof(int))) ) { sprintf(cbuf,"ERROR: unable to allocate array WS_K            in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		if( !(ws_m            = (int *)calloc(radix_vec0,sizeof(int))) ) { sprintf(cbuf,"ERROR: unable to allocate array WS_M            in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		if( !(ws_blocklen     = (int *)calloc(radix_vec0,sizeof(int))) ) { sprintf(cbuf,"ERROR: unable to allocate array WS_BLOCKLEN     in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }
		if( !(ws_blocklen_sum = (int *)calloc(radix_vec0,sizeof(int))) ) { sprintf(cbuf,"ERROR: unable to allocate array WS_BLOCKLEN_SUM in %s.\n",func); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf); }

		/*...Final DIF pass, wrapper/squaring and initial DIT pass are all done in-place.
			 This combines data from both the l1 and l2-block, except in the case ii = 0
			 for even radix_vec0, for which the l1 = 0 and l2 = 1 blocks are processed separately within
			 pairFFT_mul, i.e. we must call this routine a second time to process data in the l2-block.
		*/
		for(ii = 0; ii < radix_vec0; ii += 2)
		{
			if(ii == 0 && !(radix_vec0 & 1))	// Don't need the 2nd clause here since all radices even
				jhi = 2;
			else
				jhi = 1;
			// Loop-control-inits are same as for single-real-vector wrapper_square:
			for(j = 0; j < jhi; j++)
			{
				l = ii + j;	// This will actually leave some 'holes' in the init arrays, but we
							// don't care, as these are small and we only use the init'ed elements.
				switch(RADIX_VEC[NRADICES-1])
				{
					case 16 :
						radix16_wrapper_ini(n, radix_vec0, l, nradices_prim, radix_prim, ws_i, ws_j1, ws_j2, ws_j2_start, ws_k, ws_m, ws_blocklen, ws_blocklen_sum);
						break;
					/*
					case 32 :
						radix32_wrapper_ini(n, radix_vec0, l, nradices_prim, radix_prim, ws_i, ws_j1, ws_j2, ws_j2_start, ws_k, ws_m, ws_blocklen, ws_blocklen_sum);
						break;
					case 64 :
						radix64_wrapper_ini(n, radix_vec0, l, nradices_prim, radix_prim, ws_i, ws_j1, ws_j2, ws_j2_start, ws_k, ws_m, ws_blocklen, ws_blocklen_sum);	break;
						break;
					*/
					default :
						sprintf(cbuf,"ERROR: Final radix %d not available for %s. Halting...\n",RADIX_VEC[NRADICES-1],func);
						fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
				}
			}
		}

	}	/* end of initialization sequence.	*/

/********************************************************************************************/
/* Do forward-FFT/pointwise-square/inverse-FFT/carry-propagation/fractional-error-checking:	*/
/********************************************************************************************/
	fracmax = 0.0;

	for(input = 0; input < n_inputs; input++)
	{
	//	printf("Input vec %d of %d\n",input+1,n_inputs);
		/* Data to be processed on current pass will be either x or y-input: */
		a = ivec[input];

	  if(FORWARD_FFT_ONLY != 2)
	  {
		/* Perform the initial pass of the forward transform:	*/
		switch(radix_vec0)
		{
		case 8 :
			 radix8_dif_pass1(a,n); break;
		case 16 :
			radix16_dif_pass1(a,n); break;
		case 32 :
			radix32_dif_pass1(a,n); break;
		default :
			sprintf(cbuf,"ERROR: radix %d not available for dif_pass1. Halting...\n",radix_vec0);
			ASSERT(0,cbuf);
		}
	  }
		/* Break the remaining portion of the FFT into radix0 blocks, and in each pass of the resulting loop
		process (radix0/2) pairs of same-sized data blocks.

		In a multithreaded implementation of the loop, process NTHREADS block pairs in parallel fashion.
		If NTHREADS does not divide (radix0/2), there will be one or more under-or-unutilized threads.
		*/
		for(ii = 0; ii < radix_vec0; ii += 2) {
			pairFFT_mul_process_chunk(
				a,ab_mul,cd_mul,n,rt0,rt1,
				index,block_index,ii,nradices_prim,radix_prim,
				ws_i,ws_j1,ws_j2,ws_j2_start,ws_k,ws_m,ws_blocklen,ws_blocklen_sum,
				FORWARD_FFT_ONLY
			);
		}

		if(FORWARD_FFT_ONLY == 1) {
			// Skip DIT passes (including final one below) and proceed to next chunk:
			continue;
		} else {
			/*...Do the final inverse FFT pass, carry propagation and initial forward FFT pass in one fell swoop, er, swell loop...	*/
			ierr = 0; fracmax = 0.0;
			switch(radix_vec0) {
				// For non-modded FFTmul, no need for:    *** ******** *** *** **         **** ******* ****          *
				//	ierr = radix16_ditN_cy_dif1      (a,n,nwt,nwt_bits,wt0,wt1,si,0x0,0x0,base,baseinv,iter,&fracmax,p); break;
				case  8 :
					ierr =  radix8_ditN_cy_dif1      (a,n,  0,       0,0x0,0x0,0x0,0x0,0x0,0x0,     0x0,   0,&fracmax,0); break;
				case 16 :
					ierr = radix16_ditN_cy_dif1      (a,n,  0,       0,0x0,0x0,0x0,0x0,0x0,0x0,     0x0,   0,&fracmax,0); break;
				case 32 :
					ierr = radix32_ditN_cy_dif1      (a,n,  0,       0,0x0,0x0,0x0,0x0,0x0,0x0,     0x0,   0,&fracmax,0); break;
				default :
					sprintf(cbuf,"ERROR: radix %d not available for ditN_cy_dif1. Halting...\n",radix_vec0); fprintf(stderr,"%s", cbuf);	ASSERT(0,cbuf);
			}
			/* Nonzero remaining carries are instantly fatal: */
			ASSERT(ierr == 0, "pairFFT_mul: Fatal: carry routine return error!");

		/*...Now do the fractional error check. Any fractional part  > 0.40625 generates a warning...	*/
		// Dec 2014: Bump threshold up from ( >= 0.4 ) to ( > 0.40625 ):
			if(fracmax > 0.40625)
			{
				sprintf(cbuf, "pairFFT_mul: Roundoff warning for length %u, maxerr = %16.12f\n",n,fracmax);
				fprintf(stderr,"%s",cbuf);
			}
		}	// FORWARD_FFT_ONLY == 1 ?
	}	/* endfor(input = 0; input < n_inputs; input++) */

	// In fwd-only mode return here:
	if(FORWARD_FFT_ONLY == 1)
		return;

	/*...At the end of each iteration cycle, need to undo the initial DIF FFT pass...	*/
	switch(radix_vec0)
	{
	  case  8 :
		 radix8_dit_pass1(a,n);	break;
	  case 16 :
		radix16_dit_pass1(a,n);	break;
	  case 32 :
		radix32_dit_pass1(a,n);	break;
	  default :
		sprintf(char_str, "radix %d not available for final IFFT pass!\n",radix_vec0);
		ASSERT(0, char_str);
	}

	/*...And re-NINT the 'undo pass' data, which may differ from pure-int by some tiny amount: */
	double atmp, frac_fp, max_fp = 0.0;
	for(i=0; i < n; i++)
	{
	#ifdef USE_AVX
		j = (i & mask02) + br8[i&7];
	#elif defined(USE_SSE2)
		j = (i & mask01) + br4[i&3];
	#else
		j = i;
	#endif
		j = j + ( (j>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		atmp = a[j]*radix_inv;
		a[j] = DNINT(atmp);
		frac_fp = fabs(a[j]-atmp);
		if(frac_fp > max_fp)
			max_fp = frac_fp;
	}
	if(max_fp > 0.01)
	{
		fprintf(stderr,"%s: max_fp > 0.01! Value = %20.10f\n",func,max_fp);
		fprintf(stderr,"Check your build for inadvertent mixing of SSE2 and non-SSE2-enabled files!\n");
		ASSERT(max_fp < 0.01,"max_fp < 0.01");
	}

	// Restore input value of MODULUS_TYPE:
	MODULUS_TYPE = modtype_save;
	return;
}

/***************/

void pairFFT_mul_process_chunk(
	double a[], double ab_mul[], double cd_mul[], int n, struct complex rt0[], struct complex rt1[],
	int index[], int block_index[], int ii, int nradices_prim, int radix_prim[],
	int ws_i[], int ws_j1[], int ws_j2[], int ws_j2_start[], int ws_k[], int ws_m[], int ws_blocklen[], int ws_blocklen_sum[],
	int FORWARD_FFT_ONLY
)
{
	int thr_id = 0;	/* In unthreaded mode this must always = 0 */
	int init_sse2 = FALSE;	/* Init-calls to the various radix-pass routines presumed done prior to entry into this routine */
	int radix_vec0 = RADIX_VEC[0];
	int i,incr,istart,j,jhi,jstart,k,koffset,l,mm;
	char char_str[STR_MAX_LEN];
	jhi = 2;	// Since radix0 even (in fact power of 2) here, always do 2 blocks per loop pass

if(FORWARD_FFT_ONLY != 2)	// Cf. comments in pairFFT_mul about this
{
	for(j = 0; j < jhi; j++)
	{
		/* Get block index of the chunk of contiguous data to be processed: */
		l = block_index[ii + j];
		k    = 0;
		mm   = 1;
		incr = n/radix_vec0;

		istart = l*incr;	/* Starting location of current data-block-to-be-processed within A-array. */
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
				sprintf(cbuf,"pairFFT_mul_process_chunk: ERROR: radix %d not available for dif_pass. Halting...\n",RADIX_VEC[i]);
				ASSERT(0,cbuf);
			}

			k    += mm*radix_vec0;
			mm   *= RADIX_VEC[i];
			incr /= RADIX_VEC[i];
		}	/* end i-loop. */
	}	/* end j-loop */
}	// if(FORWARD_FFT_ONLY != 2)

	/*...Final DIF pass, wrapper/squaring and initial DIT pass are all done in-place.
	This combines data from both the l1 and l2-block, except in the case ii = 0
	for even radix_vec0, for which the l1 = 0 and l2 = 1 blocks are processed separately within
	pairFFT_mul, i.e. we must call this routine a second time to process data in the l2-block.
	*/
	for(j = 0; j < jhi; j++)
	{
		l = ii + j;

		switch(RADIX_VEC[NRADICES-1])
		{
			/* Call radix*_pairFFT_mul with INIT_ARRAYS = 0 and FORWARD_FFT_ONLY as passed: */
		  case 16 :
			radix16_pairFFT_mul(a,ab_mul,cd_mul,n,radix_vec0,rt0,rt1, nradices_prim,radix_prim, ws_i[l], ws_j1[l], ws_j2[l], ws_j2_start[l], ws_k[l], ws_m[l], ws_blocklen[l], ws_blocklen_sum[l], FALSE, FORWARD_FFT_ONLY); break;
	/*	  case 32 :
			radix32_pairFFT_mul(&a[jstart],ab_mul,cd_mul,n,radix_vec0,rt0,rt1,ii,nradices_prim,radix_prim,mm,incr, FALSE, FORWARD_FFT_ONLY); break;
	*/
		  default :
			sprintf(char_str, "pairFFT_mul_process_chunk: ERROR: radix %d not available for dyadic mul step.\n",RADIX_VEC[NRADICES-1]);
			ASSERT(0, char_str);
		}
	}
	/* In forward-FFT-only mode, do none of the IFFT passes: */
	if(FORWARD_FFT_ONLY == 1)
		return;

	/*...Rest of inverse decimation-in-time (DIT) transform. Note that during IFFT we process the radices in reverse
	order. The first array sent to each pass routine is assumed to contain the bit-reversed floating data.	*/

	for(j = 0; j < jhi; j++)
	{
		/* Get block index of the chunk of contiguous data to be processed: */
		l = block_index[ii + j];
		ASSERT(l >= 0,"pair_FFTmul_process_chunk: l >= 0");

		/* Quick-n-dirty way of generating the correct starting values of k, mm and incr -
		simply use the skeleton of the forward (DIF) loop, sans the i = NRADICES-2 pass
		(since k, mm and incr are post-incremented):
		*/
		k    = 0;
		mm   = 1;
		incr = n/radix_vec0;

		/* calculate main-array index offset here, before incr gets changed: */
		istart = l*incr;
		jstart = istart + ((istart >> DAT_BITS) << PAD_BITS );

		for(i=1; i <= NRADICES-2; i++)
		{
			k    += mm*radix_vec0;
			mm   *= RADIX_VEC[i];
			incr /= RADIX_VEC[i];
		}

		/* Now do the DIT loop, running the radices (and hence the values of k, mm and incr) in reverse: */

		for(i=NRADICES-2; i >= 1; i--)
		{
			incr *= RADIX_VEC[i];
			mm   /= RADIX_VEC[i];
			k    -= mm*radix_vec0;

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
				sprintf(cbuf,"pairFFT_mul_process_chunk: ERROR: radix %d not available for dit_pass. Halting...\n",RADIX_VEC[i]);
				ASSERT(0,cbuf);
			}
		}	/* end i-loop */
	}	/* end j-loop */
}

