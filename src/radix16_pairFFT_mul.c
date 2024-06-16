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
#if 0 //*************** To-Do: ******************
Sep 2016: See if dyadic-MUL impls here can be speeded by the following (late-1997 notes by me,
rediscovered being used as a bookmark in a diff-geom text):

First et us excerpt from the pair_square code comments:

	Given complex scalars H[j] = (x1,y1) and H[N-j] = (x2,y2) along with complex exponential
	E := exp(pi*I*j/N) = (c,s), calculates I[j] = H[j]^2 + (1 + E^2)*(H[j] - H~[N-j])^2/4
	and its complex conjugate I~, returning the former in H[j] and the latter in H[N-j].

Now my notes-to-self say this re. the analogous pair_mul of a pair of fwd-FFTed inputs H and I:

	Given complex scalar-pairs H[j],H[N-j] and I[j],I[N-j] along with complex exponential E,
	calculates M[j] = H[j]*I[j] - (1 + E^2)*(H[j] - H~[N-j])*(I[j] - I~[N-j])/4
	and its complex conjugate M~, returning the former in H[j] and the latter in H[N-j].

#endif

#include "Mlucas.h"

#define FFT_DEBUG	0

/* Use this macro in case we ever want to use this routine for auto-squaring of a pair of
independent vectors, one in the Re-parts of our main data array, the other in the Im-parts:

	A0^[1] = (C[1] + C~[2])/ 2   = [(re1+re2) + I*(im1-im2)]/2 := [rt, it],
	A1^[1] = (C[1] - C~[2])/(2I) = [(im1+im2) - I*(re1-re2)]/2 := [re,-im],

	A0^[2] = (C[2] + C~[1])/ 2   = [(re1+re2) - I*(im1-im2)]/2 := [rt,-it],
	A1^[2] = (C[2] - C~[1])/(2I) = [(im1+im2) + I*(re1-re2)]/2 := [re, im].

Outputs have Re = A0^2, Im = A1^2:

	Out1 = A0[1]^2 + I*A1[1]^2 = [rt, it]^2 + I*[re,-im]^2
	                           = [rt^2-it^2, 2*rt*it] + I*[re^2-im^2, -2*re*im] := [re0, im0]+I*[re1,-im1]
	                           = [rt^2-it^2, 2*rt*it] + [2*re*im, re^2-im^2]
	                           = [re0+im1, re1+im0]
	Out2 = A0[2]^2 + I*A1[2]^2 = [rt,-it]^2 + I*[re, im]^2
	                           = [rt^2-it^2,-2*rt*it] + I*[re^2-im^2, +2*re*im]  = [re0,-im0]+I*[re1,+im1]
	                           = [rt^2-it^2,-2*rt*it] + [-2*re*im, re^2-im^2]
	                           = [re0-im1, re1-im0]
*/
#define PAIR_SQR(_x1, _y1, _x2, _y2)\
{\
	double _rt,_it, _re,_im, _re0,_im0, _re1,_im1;\
	_rt = 0.5*(_x1 + _x2);	_it = 0.5*(_y1 - _y2);\
	_re = 0.5*(_y1 + _y2);	_im = 0.5*(_x1 - _x2);\
	_re0 = (_rt+_it)*(_rt-_it);	_rt *= _it; _im0 = _rt+_rt;	/* [re0,im0] = [rt, it]^2 = [rt^2-it^2, 2*rt*it] */\
	_re1 = (_re+_im)*(_re-_im);	_re *= _im; _im1 = _re+_re;	/* [re1,im1] = [re, im]^2 = [re^2-im^2, 2*re*im] */\
	_x1 = _re0+_im1;	_y1 = _re1+_im0;	/* outA = [re0, im0]+I*[re1,-im1] = [re0+im1]+I*[re1+im0] */\
	_x2 = _re0-_im1;	_y2 = _re1-_im0;	/* outB = [re0,-im0]+I*[re1, im1] = [re0-im1]+I*[re1-im0] */\
}

/*
In the case where we have 2 separate complex arrays, each containing the fFFT of a pair of independent
vectors as above, we have fot the 1st complex array (same as above but tweaked notation on RHS):
*** to-do: fiddle notation here to match u/v, A/B form of macro vars ***
	A0^[1] = (C[1] + C~[2])/ 2   = [(C1r+C2r) + I*(C1i-C2i)]/2 := [A0r, A0i],
	A1^[1] = (C[1] - C~[2])/(2I) = [(C1i+C2i) - I*(C1r-C2r)]/2 := [A1r,-A1i],

	A0^[2] = (C[2] + C~[1])/ 2   = [(C1r+C2r) - I*(C1i-C2i)]/2 := [A0r,-A0i],
	A1^[2] = (C[2] - C~[1])/(2I) = [(C1i+C2i) + I*(C1r-C2r)]/2 := [A1r, A1i], analogously for the B-terms.

Outputs have Re = A0*B0, Im = A1*B1:

	Out1 = A0[1]*B0[1] + I*A1[1]*B1[1] = [A0r, A0i]*[B0r, B0i] + I*[A1r,-A1i]*[B1r,-B1i]
					   = [A0r*B0r-A0i*B0i, A0r*B0i+A0i*B0r] + I*[A1r*B1r-A1i*B1i,-A1r*B1i-A1i*B1r] := [re0, im0]+I*[re1,-im1]
					   = [A0r*B0r-A0i*B0i, A0r*B0i+A0i*B0r] + [A1r*B1i+A1i*B1r, A1r*B1r-A1i*B1i] = [re0, im0]+[im1, re1]
					   = [re0+im1, re1+im0]

	Out2 = A0[2]*B0[2] + I*A1[2]*B1[2] = [A0r,-A0i]*[B0r,-B0i] + I*[A1r, A1i]*[B1r, B1i]
					   = [A0r*B0r-A0i*B0i,-A0r*B0i-A0i*B0r] + I*[A1r*B1r-A1i*B1i, A1r*B1i+A1i*B1r]  = [re0,-im0]+I*[re1,+im1]
					   = [A0r*B0r-A0i*B0i,-A0r*B0i-A0i*B0r] + [-A1r*B1i-A1i*B1r, A1r*B1r-A1i*B1i] = [re0,-im0]+[-im1,re1]
					   = [re0-im1, re1-im0];

same final form as for the squaring, just the re0,re1,im0,im1 terms require true complex-muls instead of
a complex squaring - cost [4 mul, 2 add] versus [2 mul, 3 add] for the squaring.
*/
#define PAIR_MUL(_u1,_v1, _u2,_v2, _Ar,_Ai, _Br,_Bi)\
{\
	double _ar,_ai, _br,_bi, _cr,_ci, _dr,_di, _re0,_im0, _re1,_im1;\
	/* Absorb both sets of div-by-2s into the 2nd set (precomputed-fFFT) inputs: */\
	_ar = (_u1 + _u2);	_ai = (_v1 - _v2);\
	_br = (_v1 + _v2);	_bi = (_u1 - _u2);\
	/* These assumed to come from the fixed (precomputed-fFFT) inputs: */\
	_cr = 0.25*(_Ar + _Br);	_ci = 0.25*(_Ai - _Bi);/* WILL EVENTUALLY DO THESE LIN-COMBOS AS PART OF THE fFFT */\
	_dr = 0.25*(_Ai + _Bi);	_di = 0.25*(_Ar - _Br);\
	_re0 = _ar*_cr - _ai*_ci;	_im0 = _ar*_ci + _ai*_cr;\
	_re1 = _br*_dr - _bi*_di;	_im1 = _br*_di + _bi*_dr;\
	_u1 = _re0+_im1;	_v1 = _re1+_im0;	/* outA = [re0, im0]+I*[re1,-im1] = [re0+im1]+I*[re1+im0] */\
	_u2 = _re0-_im1;	_v2 = _re1-_im0;	/* outB = [re0,-im0]+I*[re1, im1] = [re0-im1]+I*[re1-im0] */\
}

// Analogous to above macro to do [Re,Im] = [A*u, B*v] dyad-mul, this computes [Re,Im] = [A*u-B*v, C*u-D*v]:
#define ABCD_MUL(_u1,_v1, _u2,_v2, _Ar,_Ai, _Br,_Bi, _Cr,_Ci, _Dr,_Di)\
{\
	double _ar,_ai, _br,_bi, _cr,_ci, _dr,_di, _re0,_im0, _re1,_im1;\
	/* Absorb both sets of div-by-2s into the 2nd set (precomputed-fFFT) inputs: */\
	_ar = (_u1 + _u2);	_ai = (_v1 - _v2);\
	_br = (_v1 + _v2);	_bi = (_u1 - _u2);\
	/* These assumed to come from the fixed (precomputed-fFFT) inputs: */\
	_cr = 0.25*(_Ar + _Br);	_ci = 0.25*(_Ai - _Bi);\
	_dr = 0.25*(_Ai + _Bi);	_di = 0.25*(_Ar - _Br);\
	_re0 = _ar*_cr - _ai*_ci;	_im0 = _ar*_ci + _ai*_cr;\
	_re1 = _br*_dr - _bi*_di;	_im1 = _br*_di + _bi*_dr;\
	_u1 = (_re0+_im1) - (_re1+_im0);	/* [Re,Im] parts of [A*u-B*v] go into [u1,v1] */\
	_v1 = (_re0-_im1) - (_re1-_im0);\
	/* Now do C*u-D*v: */\
	_cr = 0.25*(_Cr + _Dr);	_ci = 0.25*(_Ci - _Di);\
	_dr = 0.25*(_Ci + _Di);	_di = 0.25*(_Cr - _Dr);\
	_re0 = _ar*_cr - _ai*_ci;	_im0 = _ar*_ci + _ai*_cr;\
	_re1 = _br*_dr - _bi*_di;	_im1 = _br*_di + _bi*_dr;\
	_u2 = (_re0+_im1) - (_re1+_im0);	/* [Re,Im] parts of [C*u-D*v] go into [u1,v1] */\
	_v2 = (_re0-_im1) - (_re1-_im0);\
}

/**********************/

void radix16_pairFFT_mul(
	double uv[], double ab_mul[], double cd_mul[], int n, int radix0, struct complex rt0[], struct complex rt1[],
	int nradices_prim, int radix_prim[], int ws_i, int ws_j1, int ws_j2, int ws_j2_start, int ws_k, int ws_m, int ws_blocklen, int ws_blocklen_sum,
	int INIT_ARRAYS, int FORWARD_FFT_ONLY
)
{
	const char func[] = "radix16_pairFFT_mul";
/*
!   NOTE: In the following commentary, N refers to the COMPLEX vector length (N2 in the code),
!   which is half the real vector length.
!
!...Acronyms:   DIT = Decimation In Time
!               DIF = Decimation In Frequency
!               FFT = Fast Fourier Transform, i.e. a discrete FT over the complex numbers
!
!...Final complex-radix-16 pass of a length-N complex DIF forward transform of 2 length-N *real* vectors
	(one stored in the even-index slots of A[], the other in the odds), pointwise squaring and initial
	pass of the inverse transform, all performed in one fell swoop (or better, one swell loop :) .

	***NOTES***: Since this routine is designed to be used primarily for FFT-based multiply of multiple-
	length real arrays within the same process (e.g. we might build a pair of length-N product vectors via
	successive parallel multiplies of shorter-length subproducts), instead of initializing a single set
	of bit-reversed sincos-data index arrays and storing them via static pointers within the function,
	we assume the caller has declared pointers to as many sets of such arrays (e.g. one pair of index arrays
	for length-N, another for length-N/2, etc.) as needed and passes the corresponding pointers in the
	argument list. The caller then inits said arrays via a series of INIT_ARRAYS = TRUE calls to this
	function, each of which causes the function to allocate the needed memory for a pair of such index
	tables (suitable for the specified FFT length) init the tables, and pass the pointers resulting from
	the allocs back - no other computation is done when INIT_ARRAYS = TRUE. Note that in such cases the
	init routines use the input a[] array for scratch space, i.e. it is assumed the array has been allocated
	but does not yet contain valid (or needed) data.

	When INIT_ARRAYS = FALSE, the aforementioned alloc/init step is assumed to have been done previously,
	i.e. the passed index-table pointers are assumed to point to valid index-table data (and the rt0[]
	and rt1[] arrays to valid sincos-table data) for the FFT length in question, which need not be the
	same from one call to the next - as long as the user has initialized the needed set of tables (via
	a call in INIT_ARRAYS = TRUE mode) for the FFT length in question and passes the proper table pointers
	when calling the function for each desire FFT length, there is no problem using multiple FFT lengths.

	Since we also expect to have many instances where one of the multiplicand arrays will be re-used many
	times, support two functional modes via the TRUE/FALSE status of the input argument FORWARD_FFT_ONLY:

		[1] FORWARD_FFT_ONLY = TRUE: ***** Note that by 'TRUE' here we specifically mean = 1 *****
		Compute final radix-16 pass of forward-FFT of A-array only - do not do dyadic mul, nor IFFT.

		[2a] FORWARD_FFT_ONLY = 2: This is just like FORWARD_FFT_ONLY = FALSE, but with uv-inputs also
					assumed to be already fFFTed on entry, i.e. no need for final radix-16 DIF fFFT pass.

		[2b] FORWARD_FFT_ONLY = FALSE:
		In this case uv[] is assumed to contain data that have just been partially forward-transformed (of which
		this routine does the final radix-16 DIF fFFT pass, the dyad-mul and the initial radix-16 DIT iFFT pass),
		and ab_mul[] and cd_mul[] a pair of precomputed forward-transformed datasets (e.g. corresponding to
		a quartet of a/b/c/d input multiplicands stored in packed even/odd form in this pair of vectors,
		which will be reused multiple times.) In this case, the routine	performs the final radix-16 pass
		of the forward-FFT of the uv-vector, dyad-muls the result with the corresponding a/b/c/d-array data
		to obtain FFT(a*u-b*v, c*u-d*v), and performs the initial radix-16 pass of the inverse FFT of these
		linear combinations, returning the result in uv[]. The ab_mul and cd_mul data are unaffected.

		One variant of [2]: If the cd_mul-array pointer is null, we instead compute the final radix-16 pass
		of the forward-FFT of the uv-vector, combine the result with the corresponding a/b-array data
		to obtain FFT(a*u, b*v), and perform the initial radix-16 pass of the inverse FFT of this [Re,Im] pair
		of independent pure-real dyadic products, and performs the initial radix-16 pass of the inverse FFT of these
		linear combinations, returning the result in uv[]. The ab_mul and cd_mul data are unaffected.

	Note that the value of FORWARD_FFT_ONLY is only relevant if INIT_ARRAYS = FALSE.

*/
	double *a = uv;	/* Handy alias pointer to the uv[] vector */
#ifdef USE_SSE2
	const int stride = (int)RE_IM_STRIDE << 4;	// main-array loop stride = 32 for sse2, 64 for avx
#else
	const int stride = 32;	// In this particular routine, scalar mode has same stride as SSE2
#endif
	static int nsave = 0;
	static int *index = 0x0, *index_ptmp = 0x0;	/* N2/16-length Bit-reversal index array. */
	int *itmp = 0x0;

	int rdum,idum, j1pad,j2pad,kp,l,iroot,k1,k2;
	int i,j1,j2,j2_start,k,m,blocklen,blocklen_sum;
	/*int ndivrad0m1;*/
	const double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;	/* exp[i*(twopi/16)] */
#if SYMM == 2	// Use complex-plane symmetries to reduce fraction of rt1 array actually needed
	const double mult[2] = {1.0,-1.0};
	static double nh_inv;	// Needed for "which complex quadrant?" computation
	static int nh;			// #rt1 elts in each quadrant
	int qodd;
	struct complex *cptr;
#elif defined(SYMM)
	#error Only SYMM = 2 is supported!
#endif
	double rt,it, re,im;
	double re0,im0,re1,im1;
	double cA1,cA2,cA3,cA4,cA5,cA6,cA7,cA8,cA9,cA10,cA11,cA12,cA13,cA14,cA15,sA1,sA2,sA3,sA4,sA5,sA6,sA7,sA8,sA9,sA10,sA11,sA12,sA13,sA14,sA15;
	double cB1,cB2,cB3,cB4,cB5,cB6,cB7,cB8,cB9,cB10,cB11,cB12,cB13,cB14,cB15,sB1,sB2,sB3,sB4,sB5,sB6,sB7,sB8,sB9,sB10,sB11,sB12,sB13,sB14,sB15;
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32;
	double aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r,aj1p6r,aj1p7r,aj1p8r,aj1p9r,aj1p10r,aj1p11r,aj1p12r,aj1p13r,aj1p14r,aj1p15r
			,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i,aj1p6i,aj1p7i,aj1p8i,aj1p9i,aj1p10i,aj1p11i,aj1p12i,aj1p13i,aj1p14i,aj1p15i;
	double aj2p0r,aj2p1r,aj2p2r,aj2p3r,aj2p4r,aj2p5r,aj2p6r,aj2p7r,aj2p8r,aj2p9r,aj2p10r,aj2p11r,aj2p12r,aj2p13r,aj2p14r,aj2p15r
			,aj2p0i,aj2p1i,aj2p2i,aj2p3i,aj2p4i,aj2p5i,aj2p6i,aj2p7i,aj2p8i,aj2p9i,aj2p10i,aj2p11i,aj2p12i,aj2p13i,aj2p14i,aj2p15i;
#if PFETCH
	double *addr;
#endif

	/***
	Having a separate init block for the big index array allows us to init this prior
	to anything else, using the A-array for scratch space in the call to bit_reverse_int:
	***/
	if(INIT_ARRAYS)
	{
		nsave = n;
		ASSERT(N2 == n/2, "N2 bad!");

	#if SYMM == 2	// Use complex-plane symmetries to reduce fraction of rt1 array actually needed
		nh = n/(NRT<<2);	// #rt1 elts in each quadrant
		nh_inv = 1.0/(double)nh;
	#endif

		/*
		!...Final forward (DIF) FFT pass sincos data start out in bit-reversed order.
		!   Initialize a scratch array containing N2/16 indices - again use big
		!   as-yet-unused A-array for this, but making sure the section of A used
		!   for the itmp space and that sent to the bit_reverse_int for scratch space
		!   don't overlap:
		*/
		itmp = (int *)&a[N2/16];	/* Conservatively assume an int might be as long as 8 bytes here */
		for(i=0; i < N2/16; i++)
		{
			itmp[i]=i;
		}

		/*...then bit-reverse INDEX with respect to N/16. Radices are sent to BR routine
		in opposite order from which they are processed in forward FFT.
		*/
		bit_reverse_int(itmp, N2/16, nradices_prim-4, &radix_prim[nradices_prim-5], -1, (int *)a);

	/*
	!...Allocate and initialize an index array containing N/16 indices to store the sequentially-rearranged
	!   FFT sincos data indices. (Only need N/16 of these since we need one base sincos index for every 16 complex data).
	!   We don't need a separate sincos array for the rea/complex wrapper phase, since this uses the same sincos datum
	!   as is used for the the first of each of the two blocks of 16 complex FFT data.
	*/
		if(index_ptmp != 0x0) {	// Have previously-malloc'ed local storage
			free((void *)index_ptmp);	index_ptmp=0x0;
		}
		index_ptmp = ALLOC_INT(index_ptmp, N2/16);
		ASSERT(index_ptmp != 0,"ERROR: unable to allocate array INDEX!");
		index = ALIGN_INT(index_ptmp);
	/*
	!...Now rearrange FFT sincos indices using the main loop structure as a template.
	!   The first length-2 block is a little different, since there we process the 0-15 and 16-31 array
	!   elements separately.
	*/
		index[0]=itmp [0];
		index[1]=itmp [1];

		k1 =16;	/* init uncompressed wrapper sincos array index */
		k  =2;	/* init   compressed wrapper sincos array index */

		blocklen=16;
		blocklen_sum=16;
		j2_start=96;

		for(i = nradices_prim-6; i >= 0; i--)   /* Radices get processed in reverse order here as in forward FFT. */
		{
		  kp = k1 + blocklen;

		  for(m = 0; m < (blocklen-1)>>1; m += 8)	/* Since we now process TWO 16-element sets per loop execution, only execute the loop half as many times as before. */
		  {
	/*...grab the next 2 FFT sincos indices: one for the lower (j1, in the actual loop) data block, one for the upper (j2) data block... */

			index[k  ]=itmp[(k1  )>>3];
			index[k+1]=itmp[(kp-8)>>3];

			k1 = k1+ 8;
			kp = kp- 8;

			k  = k + 2;
		  }

		  k1 = k1 + (blocklen >> 1);

		  if(j2_start == n-32)break;

		  blocklen_sum = blocklen_sum + blocklen;
		  ASSERT(i != 0,"ERROR 10!");
		  blocklen = (radix_prim[i-1]-1)*blocklen_sum;

		  j2_start = j2_start+(blocklen<<2);
		}

		j1 = 0;

		/* Restore zeros here, to prevent any barfing due to interpretation of the above integer values as floats,
		in case at some future point we find it useful to be able to re-use a part of the main a-array for scratch: */
		for(i=0; i < N2/16; i++)
		{
			itmp[i]=0;
		}

		return;
	}	/* end of inits. */

	/*...If a new runlength, should not get to this point: */
	if(n != nsave) {
		sprintf(cbuf,"ERROR: %s: INIT_ARRAYS not invoked for new runlength!",func);
		ASSERT(0,cbuf);
	}

	/* If precomputing a forward FFT of a set of inputs, make sure
	they're in the uv-vector and the abcd-multiplier vectors are null: */
	if(FORWARD_FFT_ONLY == 1 && (ab_mul != 0x0 || cd_mul != 0x0)) {
		sprintf(cbuf,"%s: FORWARD_FFT_ONLY = TRUE but non-null abcd-multiplier vectors!",func);
		ASSERT(0,cbuf);
	}

/* Init the loop-control variables: */

	i            = ws_i           ;
	j1           = ws_j1          ;
	j2           = ws_j2          ;
	j2_start     = ws_j2_start    ;
	k            = ws_k           ;
	m            = ws_m           ;
	blocklen     = ws_blocklen    ;
	blocklen_sum = ws_blocklen_sum;

	/* If j1 == 0 we need to init the loop counters; otherwise, just jump
	   right in and pick up where we left off on the previous pair of blocks:
	*/
	if(j1 > 0) {
		goto jump_in;
	}

/*
!...All but the first two radix-16 blocks are done on Mr. Henry Ford's assembly line. From there onward the blocklength
!   will always be an integer multiple of 16, i.e. we can process each block using pairs of nonoverlapping blocks of 16
!   complex data each, which is compatible to fusion with radix-16 pass routines.
*/

for(i = nradices_prim-5; i >= 0; i-- )	/* Main loop: lower bound = nradices_prim-radix_now. */
{						/* Remember, radices get processed in reverse order here as in forward FFT. */
	for(m = 0; m < (blocklen-1)>>1; m += 8) /* Since we now process TWO 16-element sets per loop execution, only execute the loop half as many times as before. */
	{
		// This tells us when we've reached the end of the current data block:
		// Apr 2014: Must store intermediate product j1*radix0 in a 64-bit int to prevent overflow!
		if(j1 && ((uint64)j1*radix0)%n == 0)
		{
			return;
		}

jump_in:	/* Entry point for all blocks but the first. */

	  j1pad = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* floating padded-array 1st element index is here */
	  j2pad = j2 + ( (j2 >> DAT_BITS) << PAD_BITS );	/* floating padded-array 2nd element index is here */

	/*************************************************************/
	/*                  1st set of sincos:                       */
	/*************************************************************/
		iroot = index[k++];
		l = iroot;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;		 /* 2*iroot */
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
	  #if SYMM == 2	// Use complex-plane symmetries to reduce fraction of rt1 array actually needed
		qodd = k2 >= nh;	k2 -= ((-qodd) & nh);	rt = mult[qodd];
		re1 = rt*rt1[k2].re;im1 = rt*rt1[k2].im;	// Mult = -+1 if root is/is-not in lower half-plane
	  #else
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
	  #endif
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cA1 =rt;	sA1 =it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 1);	/* 4*iroot */
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
	  #if SYMM == 2
		if(k2 >= nh) {	// Bizarrely, on x86 this branched version is faster than above branchless
			cptr = rt1 + k2-nh;
			re1 = -cptr->re;	im1 = -cptr->im;
		} else {
			cptr = rt1 + k2;
			re1 =  cptr->re;	im1 =  cptr->im;
		}
	  #else
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
	  #endif
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cA2 =rt;	sA2 =it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 2);	/* 8*iroot */
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
	  #if SYMM == 2
		if(k2 >= nh) {
			cptr = rt1 + k2-nh;
			re1 = -cptr->re;	im1 = -cptr->im;
		} else {
			cptr = rt1 + k2;
			re1 =  cptr->re;	im1 =  cptr->im;
		}
	  #else
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
	  #endif
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cA4 =rt;	sA4 =it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 2) + iroot;	/* 13*iroot */
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
	  #if SYMM == 2
		if(k2 >= nh) {
			cptr = rt1 + k2-nh;
			re1 = -cptr->re;	im1 = -cptr->im;
		} else {
			cptr = rt1 + k2;
			re1 =  cptr->re;	im1 =  cptr->im;
		}
	  #else
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
	  #endif
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cA8 =rt;	sA8 =it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
	  #if SYMM == 2
		if(k2 >= nh) {
			cptr = rt1 + k2-nh;
			re1 = -cptr->re;	im1 = -cptr->im;
		} else {
			cptr = rt1 + k2;
			re1 =  cptr->re;	im1 =  cptr->im;
		}
	  #else
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
	  #endif
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cA13=rt;	sA13=it;
		/* c3,5 */
		t1=cA1 *cA4 ;	t2=cA1 *sA4 ;	rt=sA1 *cA4 ;	it=sA1 *sA4;
		cA3 =t1 +it;	sA3 =t2 -rt;	cA5 =t1 -it;	sA5 =t2 +rt;

		/* c6,7,9,10 */
		t1=cA1 *cA8 ;	t2=cA1 *sA8 ;	rt=sA1 *cA8 ;	it=sA1 *sA8;
		cA7 =t1 +it;	sA7 =t2 -rt;	cA9 =t1 -it;	sA9 =t2 +rt;

		t1=cA2 *cA8 ;	t2=cA2 *sA8 ;	rt=sA2 *cA8 ;	it=sA2 *sA8;
		cA6 =t1 +it;	sA6 =t2 -rt;	cA10=t1 -it;	sA10=t2 +rt;

		/* c11,12,14,15 */
		t1=cA1 *cA13;	t2=cA1 *sA13;	rt=sA1 *cA13;	it=sA1 *sA13;
		cA12=t1 +it;	sA12=t2 -rt;	cA14=t1 -it;	sA14=t2 +rt;

		t1=cA2 *cA13;	t2=cA2 *sA13;	rt=sA2 *cA13;	it=sA2 *sA13;
		cA11=t1 +it;	sA11=t2 -rt;	cA15=t1 -it;	sA15=t2 +rt;

	/*************************************************************/
	/*                  2nd set of sincos:                       */
	/*************************************************************/
		iroot = index[k++];
		l = iroot;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;		 /* 2*iroot */
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
	  #if SYMM == 2
		if(k2 >= nh) {
			cptr = rt1 + k2-nh;
			re1 = -cptr->re;	im1 = -cptr->im;
		} else {
			cptr = rt1 + k2;
			re1 =  cptr->re;	im1 =  cptr->im;
		}
	  #else
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
	  #endif
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cB1 =rt;	sB1 =it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 1);	/* 4*iroot */
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
	  #if SYMM == 2
		if(k2 >= nh) {
			cptr = rt1 + k2-nh;
			re1 = -cptr->re;	im1 = -cptr->im;
		} else {
			cptr = rt1 + k2;
			re1 =  cptr->re;	im1 =  cptr->im;
		}
	  #else
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
	  #endif
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cB2 =rt;	sB2 =it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 2);	/* 8*iroot */
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
	  #if SYMM == 2
		if(k2 >= nh) {
			cptr = rt1 + k2-nh;
			re1 = -cptr->re;	im1 = -cptr->im;
		} else {
			cptr = rt1 + k2;
			re1 =  cptr->re;	im1 =  cptr->im;
		}
	  #else
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
	  #endif
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cB4 =rt;	sB4 =it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 2) + iroot;	/* 13*iroot */
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
	  #if SYMM == 2
		if(k2 >= nh) {
			cptr = rt1 + k2-nh;
			re1 = -cptr->re;	im1 = -cptr->im;
		} else {
			cptr = rt1 + k2;
			re1 =  cptr->re;	im1 =  cptr->im;
		}
	  #else
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
	  #endif
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cB8 =rt;	sB8 =it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		re0 = rt0[k1].re;	im0 = rt0[k1].im;
	  #if SYMM == 2
		if(k2 >= nh) {
			cptr = rt1 + k2-nh;
			re1 = -cptr->re;	im1 = -cptr->im;
		} else {
			cptr = rt1 + k2;
			re1 =  cptr->re;	im1 =  cptr->im;
		}
	  #else
		re1 = rt1[k2].re;	im1 = rt1[k2].im;
	  #endif
		rt = re0*re1 - im0*im1;	it = re0*im1 + im0*re1;
		cB13=rt;	sB13=it;
		/* c3,5 */
		t1=cB1 *cB4 ;	t2=cB1 *sB4 ;	rt=sB1 *cB4 ;	it=sB1 *sB4;
		cB3 =t1 +it;	sB3 =t2 -rt;	cB5 =t1 -it;	sB5 =t2 +rt;

		/* c6,7,9,10 */
		t1=cB1 *cB8 ;	t2=cB1 *sB8 ;	rt=sB1 *cB8 ;	it=sB1 *sB8;
		cB7 =t1 +it;	sB7 =t2 -rt;	cB9 =t1 -it;	sB9 =t2 +rt;

		t1=cB2 *cB8 ;	t2=cB2 *sB8 ;	rt=sB2 *cB8 ;	it=sB2 *sB8;
		cB6 =t1 +it;	sB6 =t2 -rt;	cB10=t1 -it;	sB10=t2 +rt;

		/* c11,12,14,15 */
		t1=cB1 *cB13;	t2=cB1 *sB13;	rt=sB1 *cB13;	it=sB1 *sB13;
		cB12=t1 +it;	sB12=t2 -rt;	cB14=t1 -it;	sB14=t2 +rt;

		t1=cB2 *cB13;	t2=cB2 *sB13;	rt=sB2 *cB13;	it=sB2 *sB13;
		cB11=t1 +it;	sB11=t2 -rt;	cB15=t1 -it;	sB15=t2 +rt;

	if(FORWARD_FFT_ONLY == 2)	// Cf. comments in pairFFT_mul about this
	{
		aj1p0r  = a[j1pad+ 0];	aj1p0i  = a[j1pad+ 1];
		aj1p1r  = a[j1pad+ 2];	aj1p1i  = a[j1pad+ 3];
		aj1p2r  = a[j1pad+ 4];	aj1p2i  = a[j1pad+ 5];
		aj1p3r  = a[j1pad+ 6];	aj1p3i  = a[j1pad+ 7];
		aj1p4r  = a[j1pad+ 8];	aj1p4i  = a[j1pad+ 9];
		aj1p5r  = a[j1pad+10];	aj1p5i  = a[j1pad+11];
		aj1p6r  = a[j1pad+12];	aj1p6i  = a[j1pad+13];
		aj1p7r  = a[j1pad+14];	aj1p7i  = a[j1pad+15];
		aj1p8r  = a[j1pad+16];	aj1p8i  = a[j1pad+17];
		aj1p9r  = a[j1pad+18];	aj1p9i  = a[j1pad+19];
		aj1p10r = a[j1pad+20];	aj1p10i = a[j1pad+21];
		aj1p11r = a[j1pad+22];	aj1p11i = a[j1pad+23];
		aj1p12r = a[j1pad+24];	aj1p12i = a[j1pad+25];
		aj1p13r = a[j1pad+26];	aj1p13i = a[j1pad+27];
		aj1p14r = a[j1pad+28];	aj1p14i = a[j1pad+29];
		aj1p15r = a[j1pad+30];	aj1p15i = a[j1pad+31];
	
		aj2p0r  = a[j2pad+ 0];	aj2p0i  = a[j2pad+ 1];
		aj2p1r  = a[j2pad+ 2];	aj2p1i  = a[j2pad+ 3];
		aj2p2r  = a[j2pad+ 4];	aj2p2i  = a[j2pad+ 5];
		aj2p3r  = a[j2pad+ 6];	aj2p3i  = a[j2pad+ 7];
		aj2p4r  = a[j2pad+ 8];	aj2p4i  = a[j2pad+ 9];
		aj2p5r  = a[j2pad+10];	aj2p5i  = a[j2pad+11];
		aj2p6r  = a[j2pad+12];	aj2p6i  = a[j2pad+13];
		aj2p7r  = a[j2pad+14];	aj2p7i  = a[j2pad+15];
		aj2p8r  = a[j2pad+16];	aj2p8i  = a[j2pad+17];
		aj2p9r  = a[j2pad+18];	aj2p9i  = a[j2pad+19];
		aj2p10r = a[j2pad+20];	aj2p10i = a[j2pad+21];
		aj2p11r = a[j2pad+22];	aj2p11i = a[j2pad+23];
		aj2p12r = a[j2pad+24];	aj2p12i = a[j2pad+25];
		aj2p13r = a[j2pad+26];	aj2p13i = a[j2pad+27];
		aj2p14r = a[j2pad+28];	aj2p14i = a[j2pad+29];
		aj2p15r = a[j2pad+30];	aj2p15i = a[j2pad+31];
		// Now proceed directly to dyadic-square step:
	} else {
	/*************************************************************/
	/*                  1st set of inputs:                       */
	/*************************************************************/
	/*
	Data layout comparison:
	A-index:0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31
	SSE:	r0	r1	i0	i1	r2	r3	i2	i3	r4	r5	i4	i5	r6	r7	i6	i7	r8	r9	i8	i9	r10	r11	i10	i11	r12	r13	i12	i13	r14	r15	i14	i15
	AVX:	r0	r1	r2	r3	i0	i1	i2	i3	r4	r5	r6	r7	i4	i5	i6	i7	r8	r9	r10	r11	i8	i9	i10	i11	r12	r13	r14	r15	i12	i13	i14	i15
	*/
		rdum = j1pad;
		idum = j1pad+RE_IM_STRIDE;

	/*...Block 1: */
		t1 =a[rdum   ];							t2 =a[idum   ];						// z0
		rt =a[rdum+16]*cA8 -a[idum+16]*sA8 ;	it =a[idum+16]*cA8 +a[rdum+16]*sA8;	// z8
		t3 =t1 -rt;	t1 =t1 +rt;
		t4 =t2 -it;	t2 =t2 +it;

		t5 =a[rdum+8 ]*cA4 -a[idum+8 ]*sA4;		t6 =a[idum+8 ]*cA4 +a[rdum+8 ]*sA4;	// z4
		rt =a[rdum+24]*cA12-a[idum+24]*sA12;	it =a[idum+24]*cA12+a[rdum+24]*sA12;// z12
		t7 =t5 -rt;	t5 =t5 +rt;
		t8 =t6 -it;	t6 =t6 +it;

		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;

		rt =t7;	t7 =t3 +t8;	t3 =t3 -t8;
				t8 =t4 -rt;	t4 =t4 +rt;

	/*...Block 2: */
		t9 =a[rdum+4 ]*cA2 -a[idum+4 ]*sA2 ;	t10=a[idum+4 ]*cA2 +a[rdum+4 ]*sA2 ;// z2
		rt =a[rdum+20]*cA10-a[idum+20]*sA10;	it =a[idum+20]*cA10+a[rdum+20]*sA10;// z10
		t11=t9 -rt;	t9 =t9 +rt;
		t12=t10-it;	t10=t10+it;

		t13=a[rdum+12]*cA6 -a[idum+12]*sA6 ;	t14=a[idum+12]*cA6 +a[rdum+12]*sA6 ;// z6
		rt =a[rdum+28]*cA14-a[idum+28]*sA14;	it =a[idum+28]*cA14+a[rdum+28]*sA14;// z14
		t15=t13-rt;	t13=t13+rt;
		t16=t14-it;	t14=t14+it;

		rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
		it =t14;	t14=t10-it;	t10=t10+it;

		rt =t15;	t15=t11+t16;t11=t11-t16;
					t16=t12-rt;	t12=t12+rt;

	/*...Block 3: */
		t17=a[rdum+2 ]*cA1 -a[idum+2 ]*sA1 ;	t18=a[idum+2 ]*cA1 +a[rdum+2 ]*sA1 ;// z1
		rt =a[rdum+18]*cA9 -a[idum+18]*sA9 ;	it =a[idum+18]*cA9 +a[rdum+18]*sA9 ;// z9
		t19=t17-rt;	t17=t17+rt;
		t20=t18-it;	t18=t18+it;

		t21=a[rdum+10]*cA5 -a[idum+10]*sA5 ;	t22=a[idum+10]*cA5 +a[rdum+10]*sA5 ;// z5
		rt =a[rdum+26]*cA13-a[idum+26]*sA13;	it =a[idum+26]*cA13+a[rdum+26]*sA13;// z13
		t23=t21-rt;	t21=t21+rt;
		t24=t22-it;	t22=t22+it;

		rt =t21;	t21=t17-rt;	t17=t17+rt;
		it =t22;	t22=t18-it;	t18=t18+it;

		rt =t23;	t23=t19+t24;t19=t19-t24;
					t24=t20-rt;	t20=t20+rt;

	/*...Block 4: */
		t25=a[rdum+6 ]*cA3 -a[idum+6 ]*sA3 ;	t26=a[idum+6 ]*cA3 +a[rdum+6 ]*sA3 ;// z3
		rt =a[rdum+22]*cA11-a[idum+22]*sA11;	it =a[idum+22]*cA11+a[rdum+22]*sA11;// z11
		t27=t25-rt;	t25=t25+rt;
		t28=t26-it;	t26=t26+it;

		t29=a[rdum+14]*cA7 -a[idum+14]*sA7 ;	t30=a[idum+14]*cA7 +a[rdum+14]*sA7 ;// z7
		rt =a[rdum+30]*cA15-a[idum+30]*sA15;	it =a[idum+30]*cA15+a[rdum+30]*sA15;// z15
		t31=t29-rt;	t29=t29+rt;
		t32=t30-it;	t30=t30+it;

		rt =t29;	t29=t25-rt;	t25=t25+rt;
		it =t30;	t30=t26-it;	t26=t26+it;

		rt =t31;	t31=t27+t32;t27=t27-t32;
					t32=t28-rt;	t28=t28+rt;
	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	!	1, exp(i* 1*twopi/16) =       ( c, s), exp(i* 2*twopi/16) = ISRT2*( 1, 1), exp(i* 3*twopi/16) =       ( s, c) (for inputs to transform block 2)
	!	1, exp(i* 2*twopi/16) = ISRT2*( 1, 1), exp(i* 4*twopi/16) =       ( 0, 1), exp(i* 6*twopi/16) = ISRT2*(-1, 1) (for inputs to transform block 3)
	!	1, exp(i* 3*twopi/16) =       ( s, c), exp(i* 6*twopi/16) = ISRT2*(-1, 1), exp(i* 9*twopi/16) =       (-c,-s) (for inputs to transform block 4).
	!  (This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
	!   I.e. do similar as above, except inputs a[j1  +p0:15:1] are replaced by t0:30:2,
	!					   a[j1+1+p0:15:1] are replaced by t1:31:2, and v.v. for outputs,
	!   and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.
	*/
	/*...Block 1: t1,9,17,25 */
		rt =t9;		t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		aj1p0r =t1+t17;	aj1p0i =t2+t18;
		aj1p1r =t1-t17;	aj1p1i =t2-t18;

		aj1p2r =t9 -t26;	aj1p2i =t10+t25;
		aj1p3r =t9 +t26;	aj1p3i =t10-t25;

	/*...Block 3: t5,13,21,29 */
		rt =t13;	t13=t5 +t14;t5 =t5 -t14;
					t14=t6 -rt;	t6 =t6 +rt;

		rt =(t21-t22)*ISRT2;t22=(t21+t22)*ISRT2;	t21=rt;
		rt =(t30+t29)*ISRT2;it =(t30-t29)*ISRT2;
		t29=t21+rt;	t21=t21-rt;
		t30=t22+it;	t22=t22-it;

		aj1p4r =t5+t21;	aj1p4i =t6+t22;
		aj1p5r =t5-t21;	aj1p5i =t6-t22;

		aj1p6r =t13-t30;	aj1p6i =t14+t29;
		aj1p7r =t13+t30;	aj1p7i =t14-t29;

	/*...Block 2: t3,11,19,27 */
		rt =(t11-t12)*ISRT2;it =(t11+t12)*ISRT2;
		t11=t3 -rt;	t3 =t3 +rt;
		t12=t4 -it;	t4 =t4 +it;

		rt =t19*c - t20*s;	t20=t20*c + t19*s;	t19=rt;
		rt =t27*s - t28*c;	it =t28*s + t27*c;
		t27=t19-rt;	t19=t19+rt;
		t28=t20-it;	t20=t20+it;

		aj1p8r =t3+t19;	aj1p8i =t4+t20;
		aj1p9r =t3-t19;	aj1p9i =t4-t20;

		aj1p10r=t11-t28;	aj1p10i=t12+t27;
		aj1p11r=t11+t28;	aj1p11i=t12-t27;

	/*...Block 4: t7,15,23,31 */
		rt =(t16+t15)*ISRT2;it =(t16-t15)*ISRT2;
		t15=t7 +rt;	t7 =t7 -rt;
		t16=t8 +it;	t8 =t8 -it;

		rt =t23*s - t24*c;	t24=t24*s + t23*c;	t23=rt;
		rt =t31*c - t32*s;	it =t32*c + t31*s;
		t31=t23+rt;	t23=t23-rt;
		t32=t24+it;	t24=t24-it;

		aj1p12r=t7+t23;	aj1p12i=t8+t24;
		aj1p13r=t7-t23;	aj1p13i=t8-t24;

		aj1p14r=t15-t32;	aj1p14i=t16+t31;
		aj1p15r=t15+t32;	aj1p15i=t16-t31;

	/*************************************************************/
	/*                  2nd set of inputs:                       */
	/*************************************************************/

		rdum = j2pad;
		idum = j2pad+RE_IM_STRIDE;

	/*...Block 1: */
		t1 =a[rdum   ];							t2 =a[idum   ];
		rt =a[rdum+16]*cB8 -a[idum+16]*sB8 ;	it =a[idum+16]*cB8 +a[rdum+16]*sB8;
		t3 =t1 -rt;	t1 =t1 +rt;
		t4 =t2 -it;	t2 =t2 +it;

		t5 =a[rdum+8 ]*cB4 -a[idum+8 ]*sB4;		t6 =a[idum+8 ]*cB4 +a[rdum+8 ]*sB4;
		rt =a[rdum+24]*cB12-a[idum+24]*sB12;	it =a[idum+24]*cB12+a[rdum+24]*sB12;
		t7 =t5 -rt;	t5 =t5 +rt;
		t8 =t6 -it;	t6 =t6 +it;

		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;

		rt =t7;	t7 =t3 +t8;	t3 =t3 -t8;
				t8 =t4 -rt;	t4 =t4 +rt;

	/*...Block 2: */
		t9 =a[rdum+4 ]*cB2 -a[idum+4 ]*sB2 ;	t10=a[idum+4 ]*cB2 +a[rdum+4 ]*sB2 ;
		rt =a[rdum+20]*cB10-a[idum+20]*sB10;	it =a[idum+20]*cB10+a[rdum+20]*sB10;
		t11=t9 -rt;	t9 =t9 +rt;
		t12=t10-it;	t10=t10+it;

		t13=a[rdum+12]*cB6 -a[idum+12]*sB6 ;	t14=a[idum+12]*cB6 +a[rdum+12]*sB6 ;
		rt =a[rdum+28]*cB14-a[idum+28]*sB14;	it =a[idum+28]*cB14+a[rdum+28]*sB14;
		t15=t13-rt;	t13=t13+rt;
		t16=t14-it;	t14=t14+it;

		rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
		it =t14;	t14=t10-it;	t10=t10+it;

		rt =t15;	t15=t11+t16;t11=t11-t16;
					t16=t12-rt;	t12=t12+rt;

	/*...Block 3: */
		t17=a[rdum+2 ]*cB1 -a[idum+2 ]*sB1 ;	t18=a[idum+2 ]*cB1 +a[rdum+2 ]*sB1 ;
		rt =a[rdum+18]*cB9 -a[idum+18]*sB9 ;	it =a[idum+18]*cB9 +a[rdum+18]*sB9 ;
		t19=t17-rt;	t17=t17+rt;
		t20=t18-it;	t18=t18+it;

		t21=a[rdum+10]*cB5 -a[idum+10]*sB5 ;	t22=a[idum+10]*cB5 +a[rdum+10]*sB5 ;
		rt =a[rdum+26]*cB13-a[idum+26]*sB13;	it =a[idum+26]*cB13+a[rdum+26]*sB13;
		t23=t21-rt;	t21=t21+rt;
		t24=t22-it;	t22=t22+it;

		rt =t21;	t21=t17-rt;	t17=t17+rt;
		it =t22;	t22=t18-it;	t18=t18+it;

		rt =t23;	t23=t19+t24;t19=t19-t24;
					t24=t20-rt;	t20=t20+rt;

	/*...Block 4: */
		t25=a[rdum+6 ]*cB3 -a[idum+6 ]*sB3 ;	t26=a[idum+6 ]*cB3 +a[rdum+6 ]*sB3 ;
		rt =a[rdum+22]*cB11-a[idum+22]*sB11;	it =a[idum+22]*cB11+a[rdum+22]*sB11;
		t27=t25-rt;	t25=t25+rt;
		t28=t26-it;	t26=t26+it;

		t29=a[rdum+14]*cB7 -a[idum+14]*sB7 ;	t30=a[idum+14]*cB7 +a[rdum+14]*sB7 ;
		rt =a[rdum+30]*cB15-a[idum+30]*sB15;	it =a[idum+30]*cB15+a[rdum+30]*sB15;
		t31=t29-rt;	t29=t29+rt;
		t32=t30-it;	t30=t30+it;

		rt =t29;	t29=t25-rt;	t25=t25+rt;
		it =t30;	t30=t26-it;	t26=t26+it;

		rt =t31;	t31=t27+t32;t27=t27-t32;
					t32=t28-rt;	t28=t28+rt;
	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	*/
	/*...Block 1: t1,9,17,25 */
		rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		aj2p0r =t1+t17;	aj2p0i =t2+t18;
		aj2p1r =t1-t17;	aj2p1i =t2-t18;

		aj2p2r =t9 -t26;	aj2p2i =t10+t25;
		aj2p3r =t9 +t26;	aj2p3i =t10-t25;

	/*...Block 3: t5,13,21,29 */
		rt =t13;	t13=t5 +t14;	t5 =t5 -t14;
			t14=t6 -rt;	t6 =t6 +rt;

		rt =(t21-t22)*ISRT2;t22=(t21+t22)*ISRT2;	t21=rt;
		rt =(t30+t29)*ISRT2;it =(t30-t29)*ISRT2;
		t29=t21+rt;	t21=t21-rt;
		t30=t22+it;	t22=t22-it;

		aj2p4r =t5+t21;	aj2p4i =t6+t22;
		aj2p5r =t5-t21;	aj2p5i =t6-t22;

		aj2p6r =t13-t30;	aj2p6i =t14+t29;
		aj2p7r =t13+t30;	aj2p7i =t14-t29;

	/*...Block 2: t3,11,19,27 */
		rt =(t11-t12)*ISRT2;it =(t11+t12)*ISRT2;
		t11=t3 -rt;	t3 =t3 +rt;
		t12=t4 -it;	t4 =t4 +it;

		rt =t19*c - t20*s;	t20=t20*c + t19*s;	t19=rt;
		rt =t27*s - t28*c;	it =t28*s + t27*c;
		t27=t19-rt;	t19=t19+rt;
		t28=t20-it;	t20=t20+it;

		aj2p8r =t3+t19;	aj2p8i =t4+t20;
		aj2p9r =t3-t19;	aj2p9i =t4-t20;

		aj2p10r=t11-t28;	aj2p10i=t12+t27;
		aj2p11r=t11+t28;	aj2p11i=t12-t27;

	/*...Block 4: t7,15,23,31 */
		rt =(t16+t15)*ISRT2;it =(t16-t15)*ISRT2;
		t15=t7 +rt;	t7 =t7 -rt;
		t16=t8 +it;	t8 =t8 -it;

		rt =t23*s - t24*c;	t24=t24*s + t23*c;	t23=rt;
		rt =t31*c - t32*s;	it =t32*c + t31*s;
		t31=t23+rt;	t23=t23-rt;
		t32=t24+it;	t24=t24-it;

		aj2p12r=t7+t23;	aj2p12i=t8+t24;
		aj2p13r=t7-t23;	aj2p13i=t8-t24;

		aj2p14r=t15-t32;	aj2p14i=t16+t31;
		aj2p15r=t15+t32;	aj2p15i=t16-t31;
	}	// if(FORWARD_FFT_ONLY == 2)

	if(FORWARD_FFT_ONLY == 1) {
		a[j1pad+ 0] = aj1p0r ;	a[j1pad+ 1] = aj1p0i ;
		a[j1pad+ 2] = aj1p1r ;	a[j1pad+ 3] = aj1p1i ;
		a[j1pad+ 4] = aj1p2r ;	a[j1pad+ 5] = aj1p2i ;
		a[j1pad+ 6] = aj1p3r ;	a[j1pad+ 7] = aj1p3i ;
		a[j1pad+ 8] = aj1p4r ;	a[j1pad+ 9] = aj1p4i ;
		a[j1pad+10] = aj1p5r ;	a[j1pad+11] = aj1p5i ;
		a[j1pad+12] = aj1p6r ;	a[j1pad+13] = aj1p6i ;
		a[j1pad+14] = aj1p7r ;	a[j1pad+15] = aj1p7i ;
		a[j1pad+16] = aj1p8r ;	a[j1pad+17] = aj1p8i ;
		a[j1pad+18] = aj1p9r ;	a[j1pad+19] = aj1p9i ;
		a[j1pad+20] = aj1p10r;	a[j1pad+21] = aj1p10i;
		a[j1pad+22] = aj1p11r;	a[j1pad+23] = aj1p11i;
		a[j1pad+24] = aj1p12r;	a[j1pad+25] = aj1p12i;
		a[j1pad+26] = aj1p13r;	a[j1pad+27] = aj1p13i;
		a[j1pad+28] = aj1p14r;	a[j1pad+29] = aj1p14i;
		a[j1pad+30] = aj1p15r;	a[j1pad+31] = aj1p15i;

		a[j2pad+ 0] = aj2p0r ;	a[j2pad+ 1] = aj2p0i ;
		a[j2pad+ 2] = aj2p1r ;	a[j2pad+ 3] = aj2p1i ;
		a[j2pad+ 4] = aj2p2r ;	a[j2pad+ 5] = aj2p2i ;
		a[j2pad+ 6] = aj2p3r ;	a[j2pad+ 7] = aj2p3i ;
		a[j2pad+ 8] = aj2p4r ;	a[j2pad+ 9] = aj2p4i ;
		a[j2pad+10] = aj2p5r ;	a[j2pad+11] = aj2p5i ;
		a[j2pad+12] = aj2p6r ;	a[j2pad+13] = aj2p6i ;
		a[j2pad+14] = aj2p7r ;	a[j2pad+15] = aj2p7i ;
		a[j2pad+16] = aj2p8r ;	a[j2pad+17] = aj2p8i ;
		a[j2pad+18] = aj2p9r ;	a[j2pad+19] = aj2p9i ;
		a[j2pad+20] = aj2p10r;	a[j2pad+21] = aj2p10i;
		a[j2pad+22] = aj2p11r;	a[j2pad+23] = aj2p11i;
		a[j2pad+24] = aj2p12r;	a[j2pad+25] = aj2p12i;
		a[j2pad+26] = aj2p13r;	a[j2pad+27] = aj2p13i;
		a[j2pad+28] = aj2p14r;	a[j2pad+29] = aj2p14i;
		a[j2pad+30] = aj2p15r;	a[j2pad+31] = aj2p15i;
/*
	printf("Writing fFFT data j1 = %u, j2 = %u:\n",j1pad,j2pad);
	printf("[%2u]: a1 = %20.5f %20.5f, a2 = %20.5f %20.5f\n", 0,a[j1pad+ 0],a[j1pad+ 1],a[j2pad+ 0],a[j2pad+ 1]);
	printf("[%2u]: a1 = %20.5f %20.5f, a2 = %20.5f %20.5f\n", 2,a[j1pad+ 2],a[j1pad+ 3],a[j2pad+ 2],a[j2pad+ 3]);
	printf("[%2u]: a1 = %20.5f %20.5f, a2 = %20.5f %20.5f\n", 4,a[j1pad+ 4],a[j1pad+ 5],a[j2pad+ 4],a[j2pad+ 5]);
	printf("[%2u]: a1 = %20.5f %20.5f, a2 = %20.5f %20.5f\n", 6,a[j1pad+ 6],a[j1pad+ 7],a[j2pad+ 6],a[j2pad+ 7]);
	printf("[%2u]: a1 = %20.5f %20.5f, a2 = %20.5f %20.5f\n", 8,a[j1pad+ 8],a[j1pad+ 9],a[j2pad+ 8],a[j2pad+ 9]);
	printf("[%2u]: a1 = %20.5f %20.5f, a2 = %20.5f %20.5f\n",10,a[j1pad+10],a[j1pad+11],a[j2pad+10],a[j2pad+11]);
	printf("[%2u]: a1 = %20.5f %20.5f, a2 = %20.5f %20.5f\n",12,a[j1pad+12],a[j1pad+13],a[j2pad+12],a[j2pad+13]);
	printf("[%2u]: a1 = %20.5f %20.5f, a2 = %20.5f %20.5f\n",14,a[j1pad+14],a[j1pad+15],a[j2pad+14],a[j2pad+15]);
	printf("[%2u]: a1 = %20.5f %20.5f, a2 = %20.5f %20.5f\n",16,a[j1pad+16],a[j1pad+17],a[j2pad+16],a[j2pad+17]);
	printf("[%2u]: a1 = %20.5f %20.5f, a2 = %20.5f %20.5f\n",18,a[j1pad+18],a[j1pad+19],a[j2pad+18],a[j2pad+19]);
	printf("[%2u]: a1 = %20.5f %20.5f, a2 = %20.5f %20.5f\n",20,a[j1pad+20],a[j1pad+21],a[j2pad+20],a[j2pad+21]);
	printf("[%2u]: a1 = %20.5f %20.5f, a2 = %20.5f %20.5f\n",22,a[j1pad+22],a[j1pad+23],a[j2pad+22],a[j2pad+23]);
	printf("[%2u]: a1 = %20.5f %20.5f, a2 = %20.5f %20.5f\n",24,a[j1pad+24],a[j1pad+25],a[j2pad+24],a[j2pad+25]);
	printf("[%2u]: a1 = %20.5f %20.5f, a2 = %20.5f %20.5f\n",26,a[j1pad+26],a[j1pad+27],a[j2pad+26],a[j2pad+27]);
	printf("[%2u]: a1 = %20.5f %20.5f, a2 = %20.5f %20.5f\n",28,a[j1pad+28],a[j1pad+29],a[j2pad+28],a[j2pad+29]);
	printf("[%2u]: a1 = %20.5f %20.5f, a2 = %20.5f %20.5f\n",30,a[j1pad+30],a[j1pad+31],a[j2pad+30],a[j2pad+31]);
*/
		goto loop;
	}
	else	// [Re,Im] parts of fFFT result contain (j,n-j)-scrambled [u,v]-data: Proceed to the dyadic-mul step, followed by iFFT:
	{
		if(!cd_mul) {	// cd_mul = 0x0: dyadic-muls to compute [Re,Im] = [A*u, B*v]:
			/* Assumes we have forward FFTs of a second pair of input vectors
			precomputed and stored in (Re,Im)ab_mul[], and performs dyadic muls
			of the forward FFT outputs with the corresponding a/b-vector
			data so as to obtain FFT(a*u, b*v).
			Since we assume we're dealing with pairs of real vectors (separately stored in the
			even/odd array slots on input to the complex FFT), we first need to unpack the separate
			outouts of the 2 real-data FFTs we have just done, then separately dyadic-mul each,
			then repack the results in preparation for a complex iFFT.

			For a length-n real-vector a[] whose elements are stored as the real parts of a length-n complex
			vector, the outputs of the length-n complex fFFT vector A[] satisfy the conjugate-symmetry
			(~ denotes complex-conjugate):

				A[n-j] = A~[j], j = 1, ..., n-1. (j = 0 is special case where no conjugate symmetry applies, nor is needed)

			Similarly, for a length-n real-vector b[] whose elements are stored as the imaginary parts of a length-n complex
			vector, the outputs - expressed in terms of the REAL-VECTOR fFFT outputs B, by explicitly pulling out the I* -
			satisfy:

				I*B[n-j] = I*B~[j], j = 1, ..., n-1. (j = 0 is again a special case).

			For length-n real-vectors A and B packed in interleaved [re = A, im = B] fashion into a length-n complex
			vector, linearity of the DFT means that the outputs of the length-n complex fFFT of the packed-input vector -
			let us call the fFFT result C - satisfy

				C[  j] = A [j] + I*B [j]	(linearity), and
				C[n-j] = A~[j] + I*B~[j]	(linearity) and conjugate symmetry), for j = 1, ..., n-1.

			Thus to unpack the separate components of A = DFT(a) and B = DFT(b) we take the conjugate of the 2nd symmetry

				C [  j] = A[j] + I*B[j],
				C~[n-j] = A[j] - I*B[j], for j = 1, ..., n-1,

			Whence

				A[j] =    (C[j] + C~[n-j])/2,
				B[j] = -I*(C[j] - C~[n-j])/2, for j = 1, ..., n-1,

			and the zero-elements are simply A[0] = Re(C[0]), B[0] = Im(C[0]) .
			We then do a pointwise (dyadic) complex multiply of each such unpacked output pair with the corresponding
			pair of data D[j] and E[j] which are similarly unpacked-on-the-fly from the precomputed-fFFT vector F[]:

				D[j] =    (F[j] + F~[n-j])/2,
				E[j] = -I*(F[j] - F~[n-j])/2, for j = 1, ..., n-1 (with zero-elements D[0] = Re(F[0]), E[0] = Im(F[0]).

			The dyadic-muls are just separate complex muls which give A*D = G, B*E = H, with {Re,Im} parts

				G[j] = { (A[j].re*D[j].re - A[j].im*D[j].im), (A[j].re*D[j].im + A[j].im*D[j].re) },
				H[j] = { (B[j].re*E[j].re - B[j].im*E[j].im), (B[j].re*E[j].im + B[j].im*E[j].re) }.

			We then re-pack these as G[j] + I*H[j] in preparation for the inverse FFT. The only additional nuance
			is that out fFFT is DIF, thus our forward-transform result is bit-reversed, and the [j,n-j] index-pair
			operations must take this into account. The details are the same as described in radi16_wrapper_square.c,
			just with the simpler dyadic-mul step which results from the separate real input vectors combined as A+I*B.
			[The examples below express the transforms as matrix-multiply DFTs, thus all outputs remain ordered.]

			EXAMPLE 1: Convolution-based multiply of a pair of pure-real (single) input vector:
			Multiply 12 and 23 the FFT way, representing them as base-10 vectors, i.e. we form our input
			vectors simply by separating out the digits of the above base-10 representation of the numbers.
			Our zero-padded input vectors are A = [2,1,0,0] and B = [3,2,0,0], where the least-significant
			digit of each number is leftmost. The fFFT of a length-4 vector [a,b,c,d] is

			/+1 +1 +1 +1\ /a\    a+b  +    c+d
			|+1 +I -1 -I|*|b| = (a-c) + I*(b-d)
			|+1 -1 +1 -1| |c|    a-b  +    c-d
			\+1 -I -1 +I/ \d/   (a-c) - I*(b-d)

			where I = sqrt(-1) is the usual imaginary constant. Doing this for our two input vectors, our
			forward-transformed input vectors are A^ = [3, 2+I, 1, 2-I] and B^ = [5, 3+2I, 1, 3-2I]. The
			Fourier-transformed discrete convolution of A and B is then simply A^*B^, where the '*' means
			component-by-component (the fancy word is 'dyadic') multiply, which gives [15, 4+7I, 1, 4-7I].
			To get the digits of the product we're after, we need to inverse-Fourier-transform this vector.
			For length-4 vectors, the iFFT looks just like the fFFT, but with the signs on the i-terms in
			the above matrix switched and a factor of one-fourth multiplying the whole thing. (We always
			have the factor of 1/N for a length-N inverse transform, since we require that doing an fFFT
			of a vector followed by an iFFT simply give us back our original vector).

			That is, the length-4 iFFT of our vector [15, 4+7I, 1, 4-7I] has entries

			[(19+7I) +   (5- 7I)]/4 = 24/4 = 6,
			[(14   ) - I*(  14I)]/4 = 28/4 = 7,
			[(11-7I) +   (-3+7I)]/4 =  8/4 = 2,
			[(14   ) + I*(  14I)]/4 =  0/4 = 0.

			Since all the output digits happen to be less than 10 and nonnegative, need no carry step, and
			since the outputs are least-significant first, result (written in normal decimal order) is 276.

			-------------------------------------------------------------------------------------------------

			EXAMPLE 2: Again a convo-multiply of a pair of pure-real (single) input vectors, for 31*77 = 2387.
			A = [1,3,0,0] and B = [7,7,0,0], fFFT gives A^ = [4, 1+3I, -2, 1-3I], B^ = [14, 7+7I, 0, 7-7I],
			A^*B^ = [56, -14+28I, 0, -14-28I], iFFT gives [56-28, 56-I*(56I), 56+28, 56+I*(56I)]/4 = [7,28,21,0],
			carry step gives [7,8,1,0] + [0,0,2,2] = [7,8,3,2], as expected.

			-------------------------------------------------------------------------------------------------

			EXAMPLE 3: Now do same pair of real-multiplies, but pairwise pack the input-vectors into complex
			vectors. 12 and 31 get packed into C = [2+I,1+3I,0,0]; 23 and 77 into D = [3+7I,2+7I,0,0].
			The fFFT gives
				C^ = [3, 2+ I, 1, 2- I] + I*[ 4, 1+3I, -2, 1-3I] = [3+ 4I, -1+2I, 1-2I,  5   ],
				D^ = [5, 3+2I, 1, 3-2I] + I*[14, 7+7I,  0, 7-7I] = [5+14I, -4+9I, 1   , 10+5I].
			Now before doing the dyadic-mul we need to unpack the separate fFFTs of our packed-input-vectors:

			DC: A0^[0] = Re(C[0]), A1^[0] = Im(C[0]) .
				A0^[j] = (C[j] + C~[n-j])/2,	A1^[j] = (C[j] - C~[n-j])/2, for j = 1, ..., n-1.
			Term-by-term this gives
				A0^[0] = Re(C[0]) =         3  ,	A1^[0] = Im(C[0]) =                        4   ,
				A0^[1] = (C[1] + C~[3])/2 = 2+I,	A1^[1] = (C[1] - C~[3])/(2I) = -I*(-3+I) = 1+3I,
				A0^[2] = Re(C[2]) =         1  ,	A1^[2] = Im(C[2]) =                       -2   ,
				A0^[3] = (C[3] + C~[1])/2 = 2-I,	A1^[3] = (C[3] - C~[1])/(2I) = -I*( 3+I) = 1-3I, all of which check.
			Subbing A -> B and C -> D in the above we get
				B0^[0] = Re(D[0]) =         5   ,	B1^[0] = Im(D[0]) =                        14   ,
				B0^[1] = (D[1] + D~[3])/2 = 3+2I,	B1^[1] = (D[1] - D~[3])/(2I) = -I*(-7+7I) = 7+7I,
				B0^[2] = Re(D[2]) =         1   ,	B1^[2] = Im(D[2]) =                         0   ,
				B0^[3] = (D[3] + D~[1])/2 = 3-2I,	B1^[3] = (D[3] - D~[1])/(2I) = -I*( 7+7I) = 7-7I, all check.

			Our separate dyadic-mul outputs are then again simply recombined as [A0^ * B0^] + I*[A1^ * B1^] to give the
			input vector for our iFFT:
			[NOTE: No need to 'precombine' terms to prevent iFFT outputs from ending similarly linear-combined
			as the fFFT ones, nor any need for post-transform data unscrambling - very nice!]

				 = [3*5, (2+I)*(3+2I), 1*1, (2-I)*(3-2I)] + I*[4*14, (1+3I)*(7+7I), -2*0, (1-3I)*(7-7I)]
				 = [15, 4+7I, 1, 4-7I] + I*[56, -14+28I, 0, -14-28I]
				 = [15+56I, -24-7I, 1, 32-21I] .
			iFFT ==> [24+28I, 28+112I, 8+84I, 0], and simply separately carry the Re,Im-parts to get our 2 mul results.

			If instead we we simply FFT-squaring 2 real vectors we would instead compute [A0^ * A0^] + I*[A1^ * A1^]
				 = [3^2, (2+I)^2, 1^2, (2-I)^2] + I*[4^2, (1+3I)^2, -2^2, (1-3I)^2]
				 = [9, 3+4I, 1, 3-4I] + I*[16, -8+6I, 4, -8-6I]
				 = [9, 3+4I, 1, 3-4I] + [16I, -6-8I, 4I, 6-8I]
				 = [9+16I, -3-4I, 1+4I, 9-12I]
			iFFT ==> [16+4I, 16+24I, 4+36I, 0]/4 = [4+I, 4+6I, 1+9I, 0], again separately carry Re,Im-parts to get our
			2 mul results, 12^2 = 144 and 31^2 = 961.
			*/
			// First pair of blocks is special, in that there is no cross-mixing of data between the blocks:
			if(j1==0) {
			/* Block 1: */
				aj1p0r *= ab_mul[j1pad+ 0];	aj1p0i *= ab_mul[j1pad+ 1];	//...j = 0 is done separately...
				aj1p1r *= ab_mul[j1pad+ 2];	aj1p1i *= ab_mul[j1pad+ 3];	//...as is j = N/2, which under BR ends up in j=1 complex data slot [real-data slots 2,3]
				/* Remaining index-pairings are as in radi16_wrapper_square.c, but with above simplified dyadic-mul, e.g.:
					A0^[2] = (C[2] + C~[3])/ 2   = [(re2+re3) + I*(im2-im3)]/2 = [rt, it],
					A1^[2] = (C[2] - C~[3])/(2I) = [(im2+im3) - I*(re2-re3)]/2 = [re,-im],
	
					A0^[3] = (C[3] + C~[2])/ 2   = [(re2+re3) - I*(im2-im3)]/2 = [rt,-it],
					A1^[3] = (C[3] - C~[2])/(2I) = [(im2+im3) + I*(re2-re3)]/2 = [re, im], analogously for the B-terms.
				*/
				PAIR_MUL(aj1p2r ,aj1p2i , aj1p3r ,aj1p3i , ab_mul[j1pad+ 4],ab_mul[j1pad+ 5], ab_mul[j1pad+ 6],ab_mul[j1pad+ 7]);
				PAIR_MUL(aj1p4r ,aj1p4i , aj1p7r ,aj1p7i , ab_mul[j1pad+ 8],ab_mul[j1pad+ 9], ab_mul[j1pad+14],ab_mul[j1pad+15]);
				PAIR_MUL(aj1p6r ,aj1p6i , aj1p5r ,aj1p5i , ab_mul[j1pad+12],ab_mul[j1pad+13], ab_mul[j1pad+10],ab_mul[j1pad+11]);
				PAIR_MUL(aj1p8r ,aj1p8i , aj1p15r,aj1p15i, ab_mul[j1pad+16],ab_mul[j1pad+17], ab_mul[j1pad+30],ab_mul[j1pad+31]);
				PAIR_MUL(aj1p10r,aj1p10i, aj1p13r,aj1p13i, ab_mul[j1pad+20],ab_mul[j1pad+21], ab_mul[j1pad+26],ab_mul[j1pad+27]);
				PAIR_MUL(aj1p12r,aj1p12i, aj1p11r,aj1p11i, ab_mul[j1pad+24],ab_mul[j1pad+25], ab_mul[j1pad+22],ab_mul[j1pad+23]);
				PAIR_MUL(aj1p14r,aj1p14i, aj1p9r ,aj1p9i , ab_mul[j1pad+28],ab_mul[j1pad+29], ab_mul[j1pad+18],ab_mul[j1pad+19]);
			/* Block 2: */				
				PAIR_MUL(aj2p0r ,aj2p0i , aj2p15r,aj2p15i, ab_mul[j2pad+ 0],ab_mul[j2pad+ 1], ab_mul[j2pad+30],ab_mul[j2pad+31]);
				PAIR_MUL(aj2p2r ,aj2p2i , aj2p13r,aj2p13i, ab_mul[j2pad+ 4],ab_mul[j2pad+ 5], ab_mul[j2pad+26],ab_mul[j2pad+27]);
				PAIR_MUL(aj2p4r ,aj2p4i , aj2p11r,aj2p11i, ab_mul[j2pad+ 8],ab_mul[j2pad+ 9], ab_mul[j2pad+22],ab_mul[j2pad+23]);
				PAIR_MUL(aj2p6r ,aj2p6i , aj2p9r ,aj2p9i , ab_mul[j2pad+12],ab_mul[j2pad+13], ab_mul[j2pad+18],ab_mul[j2pad+19]);
				PAIR_MUL(aj2p8r ,aj2p8i , aj2p7r ,aj2p7i , ab_mul[j2pad+16],ab_mul[j2pad+17], ab_mul[j2pad+14],ab_mul[j2pad+15]);
				PAIR_MUL(aj2p10r,aj2p10i, aj2p5r ,aj2p5i , ab_mul[j2pad+20],ab_mul[j2pad+21], ab_mul[j2pad+10],ab_mul[j2pad+11]);
				PAIR_MUL(aj2p12r,aj2p12i, aj2p3r ,aj2p3i , ab_mul[j2pad+24],ab_mul[j2pad+25], ab_mul[j2pad+ 6],ab_mul[j2pad+ 7]);
				PAIR_MUL(aj2p14r,aj2p14i, aj2p1r ,aj2p1i , ab_mul[j2pad+28],ab_mul[j2pad+29], ab_mul[j2pad+ 2],ab_mul[j2pad+ 3]);
			} else {	// (j1 == 0) = false:
				// All but the first 2 blocks mix each others' data:
				PAIR_MUL(aj1p0r ,aj1p0i , aj2p15r,aj2p15i, ab_mul[j1pad+ 0],ab_mul[j1pad+ 1], ab_mul[j2pad+30],ab_mul[j2pad+31]);
				PAIR_MUL(aj1p2r ,aj1p2i , aj2p13r,aj2p13i, ab_mul[j1pad+ 4],ab_mul[j1pad+ 5], ab_mul[j2pad+26],ab_mul[j2pad+27]);
				PAIR_MUL(aj1p4r ,aj1p4i , aj2p11r,aj2p11i, ab_mul[j1pad+ 8],ab_mul[j1pad+ 9], ab_mul[j2pad+22],ab_mul[j2pad+23]);
				PAIR_MUL(aj1p6r ,aj1p6i , aj2p9r ,aj2p9i , ab_mul[j1pad+12],ab_mul[j1pad+13], ab_mul[j2pad+18],ab_mul[j2pad+19]);
				PAIR_MUL(aj1p8r ,aj1p8i , aj2p7r ,aj2p7i , ab_mul[j1pad+16],ab_mul[j1pad+17], ab_mul[j2pad+14],ab_mul[j2pad+15]);
				PAIR_MUL(aj1p10r,aj1p10i, aj2p5r ,aj2p5i , ab_mul[j1pad+20],ab_mul[j1pad+21], ab_mul[j2pad+10],ab_mul[j2pad+11]);
				PAIR_MUL(aj1p12r,aj1p12i, aj2p3r ,aj2p3i , ab_mul[j1pad+24],ab_mul[j1pad+25], ab_mul[j2pad+ 6],ab_mul[j2pad+ 7]);
				PAIR_MUL(aj1p14r,aj1p14i, aj2p1r ,aj2p1i , ab_mul[j1pad+28],ab_mul[j1pad+29], ab_mul[j2pad+ 2],ab_mul[j2pad+ 3]);

				PAIR_MUL(aj2p0r ,aj2p0i , aj1p15r,aj1p15i, ab_mul[j2pad+ 0],ab_mul[j2pad+ 1], ab_mul[j1pad+30],ab_mul[j1pad+31]);
				PAIR_MUL(aj2p2r ,aj2p2i , aj1p13r,aj1p13i, ab_mul[j2pad+ 4],ab_mul[j2pad+ 5], ab_mul[j1pad+26],ab_mul[j1pad+27]);
				PAIR_MUL(aj2p4r ,aj2p4i , aj1p11r,aj1p11i, ab_mul[j2pad+ 8],ab_mul[j2pad+ 9], ab_mul[j1pad+22],ab_mul[j1pad+23]);
				PAIR_MUL(aj2p6r ,aj2p6i , aj1p9r ,aj1p9i , ab_mul[j2pad+12],ab_mul[j2pad+13], ab_mul[j1pad+18],ab_mul[j1pad+19]);
				PAIR_MUL(aj2p8r ,aj2p8i , aj1p7r ,aj1p7i , ab_mul[j2pad+16],ab_mul[j2pad+17], ab_mul[j1pad+14],ab_mul[j1pad+15]);
				PAIR_MUL(aj2p10r,aj2p10i, aj1p5r ,aj1p5i , ab_mul[j2pad+20],ab_mul[j2pad+21], ab_mul[j1pad+10],ab_mul[j1pad+11]);
				PAIR_MUL(aj2p12r,aj2p12i, aj1p3r ,aj1p3i , ab_mul[j2pad+24],ab_mul[j2pad+25], ab_mul[j1pad+ 6],ab_mul[j1pad+ 7]);
				PAIR_MUL(aj2p14r,aj2p14i, aj1p1r ,aj1p1i , ab_mul[j2pad+28],ab_mul[j2pad+29], ab_mul[j1pad+ 2],ab_mul[j1pad+ 3]);
			}
		} else {	// cd_mul != 0x0: dyadic-muls to compute [Re,Im] = [a*u-b*v, c*u-d*v]:
			// Dec 2015: Despite all my efforts, simply not yet able to wring out remaining bug(s) in indexing scheme
			// here. If and when I do finally get things working, also need to fuse the 2 x PAIR_MUL occurrences on
			// each line into a working single ABCD_MUL macro, which avoids the work-duplication of the 2 x PAIR_MUL:
			ASSERT(0, "Linear-combo algorithm not yet working!");
			/*
			Dyadic muls of the forward FFT outputs with the corresponding a/b and c/d-vector data so as to
			obtain FFT(a*u-b*v, c*u-d*v). u,v in ajp*r,i; a,b in ab_mul[even,odd]; c,d in cd_mul[even,odd]:

			1. Complex fFFT of a pair of real invecs [a,b] packed into a complex array as [a+I*b] -> ab_mul[] data here;
			2. Complex fFFT of a pair of real invecs [c,d] packed into a complex array as [c+I*d] -> cd_mul[] data here;
			3. Complex fFFT of a pair of real invecs [x,y] packed into a complex array as [x+I*y] -> aj1p*r,i data here;
			4. The PAIR_MUL macro takes a quartet of complex fFFT outputs - 2 from [3] and 2 from either [1] or [2] and
			   effects a dyadic-mul sequence which prepares the results for iFFT of packed-real-vector-pair, said
			   results overwrite the data pair of complex outputs from [3], i.e. the 4 float-double aj1p*r,i data.

			Now more in-depth, the PAIR_MUL macro takes 4 complex data (xr,xi, yr,yi, ar,ai, br,bi) and dyad-muls them,
			under the assumption that the r,i-terms are Re,Im-part outputs of the aforementione real-input-vec-pair fFFTs
			and the r,i-parts of the 2 complex PAIR_MUL outputs will be interleaved in the Re,Im slots of the input vec
			to another real-input-vec-pair FFT, this time an inverse one.
			
			When using the same machinery to handle the linear-combination case needed for the fast FFT-GCD-muls, we are
			mixing Re and Im outputs of the PAIR_MUL macro so as to generate the inputs to the ensuing iFFT intended to
			produce the desired linear combos [a*u-b*v, c*u-d*v] in the Re,Im parts of the output. BUT, the FFT-of-pair-
			of-real-input-vectors scheme is such that the 2nd vector, whose elements get put into he Im-parts of our
			complex-FFT-vector, PRODUCES THE VECTOR-#2 OUTPUTS IN INDEX-REVERSED ORDER. When we keep the Re,Im-parts
			separate as in the [a*u, b*v] case above we need no explicit handling of this because the fFFT's reversal
			simply gets undone by the ensuing iFFT, but in the linear-combo [a*u-b*v, c*u-d*v] case we must combine
			terms of form (a*u)[j] - (b*v)[n-j]. Plus our fFFT outputs here are bit-reversed - but this seemingly-
			nightmare scenario is mitigated by the fact the the index-pair patterns of terms fed to the PAIR_MUL macro
			already embody the scrambling needed to handle [j,n-j] + BR index scrambling, as noted in commentary below.
			*/
			if(j1==0) {
				/* All but the first 2 blocks mix each others' data: */
			// Reverse-indexing of Im-parts of fFFT outputs means we must it <-> im and im0 <-> im1 in our re* - im* following the PAIR_MULs:
			/* Block 1: */
				re = aj1p0r*ab_mul[j1pad+ 0] - aj1p0i*ab_mul[j1pad+ 1];	aj1p0i = aj1p0r*cd_mul[j1pad+ 0] - aj1p0i*cd_mul[j1pad+ 1];	aj1p0r = re;
				re = aj1p1r*ab_mul[j1pad+ 2] - aj1p1i*ab_mul[j1pad+ 3];	aj1p1i = aj1p1r*cd_mul[j1pad+ 2] - aj1p1i*cd_mul[j1pad+ 3];	aj1p1r = re;
				rt = re0 = aj1p2r ; it = im0 = aj1p2i ; re = re1 = aj1p3r ; im = im1 = aj1p3i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j1pad+ 4],ab_mul[j1pad+ 5], ab_mul[j1pad+ 6],ab_mul[j1pad+ 7]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j1pad+ 4],cd_mul[j1pad+ 5], cd_mul[j1pad+ 6],cd_mul[j1pad+ 7]);	aj1p2r = rt -im ; aj1p2i = re -it ; aj1p3r = re1-im0; aj1p3i = re0-im1;
				rt = re0 = aj1p4r ; it = im0 = aj1p4i ; re = re1 = aj1p7r ; im = im1 = aj1p7i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j1pad+ 8],ab_mul[j1pad+ 9], ab_mul[j1pad+14],ab_mul[j1pad+15]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j1pad+ 8],cd_mul[j1pad+ 9], cd_mul[j1pad+14],cd_mul[j1pad+15]);	aj1p4r = rt -im ; aj1p4i = re -it ; aj1p7r = re1-im0; aj1p7i = re0-im1;
				rt = re0 = aj1p6r ; it = im0 = aj1p6i ; re = re1 = aj1p5r ; im = im1 = aj1p5i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j1pad+12],ab_mul[j1pad+13], ab_mul[j1pad+10],ab_mul[j1pad+11]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j1pad+12],cd_mul[j1pad+13], cd_mul[j1pad+10],cd_mul[j1pad+11]);	aj1p6r = rt -im ; aj1p6i = re -it ; aj1p5r = re1-im0; aj1p5i = re0-im1;
				rt = re0 = aj1p8r ; it = im0 = aj1p8i ; re = re1 = aj1p15r; im = im1 = aj1p15i;	PAIR_MUL(rt ,it , re ,im , ab_mul[j1pad+16],ab_mul[j1pad+17], ab_mul[j1pad+30],ab_mul[j1pad+31]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j1pad+16],cd_mul[j1pad+17], cd_mul[j1pad+30],cd_mul[j1pad+31]);	aj1p8r = rt -im ; aj1p8i = re -it ; aj1p15r= re1-im0; aj1p15i= re0-im1;
				rt = re0 = aj1p10r; it = im0 = aj1p10i; re = re1 = aj1p13r; im = im1 = aj1p13i;	PAIR_MUL(rt ,it , re ,im , ab_mul[j1pad+20],ab_mul[j1pad+21], ab_mul[j1pad+26],ab_mul[j1pad+27]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j1pad+20],cd_mul[j1pad+21], cd_mul[j1pad+26],cd_mul[j1pad+27]);	aj1p10r= rt -im ; aj1p10i= re -it ; aj1p13r= re1-im0; aj1p13i= re0-im1;
				rt = re0 = aj1p12r; it = im0 = aj1p12i; re = re1 = aj1p11r; im = im1 = aj1p11i;	PAIR_MUL(rt ,it , re ,im , ab_mul[j1pad+24],ab_mul[j1pad+25], ab_mul[j1pad+22],ab_mul[j1pad+23]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j1pad+24],cd_mul[j1pad+25], cd_mul[j1pad+22],cd_mul[j1pad+23]);	aj1p12r= rt -im ; aj1p12i= re -it ; aj1p11r= re1-im0; aj1p11i= re0-im1;
				rt = re0 = aj1p14r; it = im0 = aj1p14i; re = re1 = aj1p9r ; im = im1 = aj1p9i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j1pad+28],ab_mul[j1pad+29], ab_mul[j1pad+18],ab_mul[j1pad+19]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j1pad+28],cd_mul[j1pad+29], cd_mul[j1pad+18],cd_mul[j1pad+19]);	aj1p14r= rt -im ; aj1p14i= re -it ; aj1p9r = re1-im0; aj1p9i = re0-im1;
			/* Block 2: */
				rt = re0 = aj2p0r ; it = im0 = aj2p0i ; re = re1 = aj2p15r; im = im1 = aj2p15i;	PAIR_MUL(rt ,it , re ,im , ab_mul[j2pad+ 0],ab_mul[j2pad+ 1], ab_mul[j2pad+30],ab_mul[j2pad+31]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j2pad+ 0],cd_mul[j2pad+ 1], cd_mul[j2pad+30],cd_mul[j2pad+31]);	aj2p0r = rt -im ; aj2p0i = re -it ; aj2p15r= re1-im0; aj2p15i= re0-im1;
				rt = re0 = aj2p2r ; it = im0 = aj2p2i ; re = re1 = aj2p13r; im = im1 = aj2p13i;	PAIR_MUL(rt ,it , re ,im , ab_mul[j2pad+ 4],ab_mul[j2pad+ 5], ab_mul[j2pad+26],ab_mul[j2pad+27]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j2pad+ 4],cd_mul[j2pad+ 5], cd_mul[j2pad+26],cd_mul[j2pad+27]);	aj2p2r = rt -im ; aj2p2i = re -it ; aj2p13r= re1-im0; aj2p13i= re0-im1;
				rt = re0 = aj2p4r ; it = im0 = aj2p4i ; re = re1 = aj2p11r; im = im1 = aj2p11i;	PAIR_MUL(rt ,it , re ,im , ab_mul[j2pad+ 8],ab_mul[j2pad+ 9], ab_mul[j2pad+22],ab_mul[j2pad+23]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j2pad+ 8],cd_mul[j2pad+ 9], cd_mul[j2pad+22],cd_mul[j2pad+23]);	aj2p4r = rt -im ; aj2p4i = re -it ; aj2p11r= re1-im0; aj2p11i= re0-im1;
				rt = re0 = aj2p6r ; it = im0 = aj2p6i ; re = re1 = aj2p9r ; im = im1 = aj2p9i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j2pad+12],ab_mul[j2pad+13], ab_mul[j2pad+18],ab_mul[j2pad+19]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j2pad+12],cd_mul[j2pad+13], cd_mul[j2pad+18],cd_mul[j2pad+19]);	aj2p6r = rt -im ; aj2p6i = re -it ; aj2p9r = re1-im0; aj2p9i = re0-im1;
				rt = re0 = aj2p8r ; it = im0 = aj2p8i ; re = re1 = aj2p7r ; im = im1 = aj2p7i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j2pad+16],ab_mul[j2pad+17], ab_mul[j2pad+14],ab_mul[j2pad+15]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j2pad+16],cd_mul[j2pad+17], cd_mul[j2pad+14],cd_mul[j2pad+15]);	aj2p8r = rt -im ; aj2p8i = re -it ; aj2p7r = re1-im0; aj2p7i = re0-im1;
				rt = re0 = aj2p10r; it = im0 = aj2p10i; re = re1 = aj2p5r ; im = im1 = aj2p5i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j2pad+20],ab_mul[j2pad+21], ab_mul[j2pad+10],ab_mul[j2pad+11]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j2pad+20],cd_mul[j2pad+21], cd_mul[j2pad+10],cd_mul[j2pad+11]);	aj2p10r= rt -im ; aj2p10i= re -it ; aj2p5r = re1-im0; aj2p5i = re0-im1;
				rt = re0 = aj2p12r; it = im0 = aj2p12i; re = re1 = aj2p3r ; im = im1 = aj2p3i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j2pad+24],ab_mul[j2pad+25], ab_mul[j2pad+ 6],ab_mul[j2pad+ 7]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j2pad+24],cd_mul[j2pad+25], cd_mul[j2pad+ 6],cd_mul[j2pad+ 7]);	aj2p12r= rt -im ; aj2p12i= re -it ; aj2p3r = re1-im0; aj2p3i = re0-im1;
				rt = re0 = aj2p14r; it = im0 = aj2p14i; re = re1 = aj2p1r ; im = im1 = aj2p1i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j2pad+28],ab_mul[j2pad+29], ab_mul[j2pad+ 2],ab_mul[j2pad+ 3]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j2pad+28],cd_mul[j2pad+29], cd_mul[j2pad+ 2],cd_mul[j2pad+ 3]);	aj2p14r= rt -im ; aj2p14i= re -it ; aj2p1r = re1-im0; aj2p1i = re0-im1;
			} else {	/* (j1 == 0) = false: */
				/* All but the first 2 blocks mix each others' data: */
				rt = re0 = aj1p0r ; it = im0 = aj1p0i ; re = re1 = aj2p15r; im = im1 = aj2p15i;	PAIR_MUL(rt ,it , re ,im , ab_mul[j1pad+ 0],ab_mul[j1pad+ 1], ab_mul[j2pad+30],ab_mul[j2pad+31]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j1pad+ 0],cd_mul[j1pad+ 1], cd_mul[j2pad+30],cd_mul[j2pad+31]);	aj1p0r = rt -im ; aj1p0i = re -it ; aj2p15r= re1-im0; aj2p15i= re0-im1;
				rt = re0 = aj1p2r ; it = im0 = aj1p2i ; re = re1 = aj2p13r; im = im1 = aj2p13i;	PAIR_MUL(rt ,it , re ,im , ab_mul[j1pad+ 4],ab_mul[j1pad+ 5], ab_mul[j2pad+26],ab_mul[j2pad+27]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j1pad+ 4],cd_mul[j1pad+ 5], cd_mul[j2pad+26],cd_mul[j2pad+27]);	aj1p2r = rt -im ; aj1p2i = re -it ; aj2p13r= re1-im0; aj2p13i= re0-im1;
				rt = re0 = aj1p4r ; it = im0 = aj1p4i ; re = re1 = aj2p11r; im = im1 = aj2p11i;	PAIR_MUL(rt ,it , re ,im , ab_mul[j1pad+ 8],ab_mul[j1pad+ 9], ab_mul[j2pad+22],ab_mul[j2pad+23]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j1pad+ 8],cd_mul[j1pad+ 9], cd_mul[j2pad+22],cd_mul[j2pad+23]);	aj1p4r = rt -im ; aj1p4i = re -it ; aj2p11r= re1-im0; aj2p11i= re0-im1;
				rt = re0 = aj1p6r ; it = im0 = aj1p6i ; re = re1 = aj2p9r ; im = im1 = aj2p9i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j1pad+12],ab_mul[j1pad+13], ab_mul[j2pad+18],ab_mul[j2pad+19]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j1pad+12],cd_mul[j1pad+13], cd_mul[j2pad+18],cd_mul[j2pad+19]);	aj1p6r = rt -im ; aj1p6i = re -it ; aj2p9r = re1-im0; aj2p9i = re0-im1;
				rt = re0 = aj1p8r ; it = im0 = aj1p8i ; re = re1 = aj2p7r ; im = im1 = aj2p7i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j1pad+16],ab_mul[j1pad+17], ab_mul[j2pad+14],ab_mul[j2pad+15]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j1pad+16],cd_mul[j1pad+17], cd_mul[j2pad+14],cd_mul[j2pad+15]);	aj1p8r = rt -im ; aj1p8i = re -it ; aj2p7r = re1-im0; aj2p7i = re0-im1;
				rt = re0 = aj1p10r; it = im0 = aj1p10i; re = re1 = aj2p5r ; im = im1 = aj2p5i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j1pad+20],ab_mul[j1pad+21], ab_mul[j2pad+10],ab_mul[j2pad+11]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j1pad+20],cd_mul[j1pad+21], cd_mul[j2pad+10],cd_mul[j2pad+11]);	aj1p10r= rt -im ; aj1p10i= re -it ; aj2p5r = re1-im0; aj2p5i = re0-im1;
				rt = re0 = aj1p12r; it = im0 = aj1p12i; re = re1 = aj2p3r ; im = im1 = aj2p3i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j1pad+24],ab_mul[j1pad+25], ab_mul[j2pad+ 6],ab_mul[j2pad+ 7]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j1pad+24],cd_mul[j1pad+25], cd_mul[j2pad+ 6],cd_mul[j2pad+ 7]);	aj1p12r= rt -im ; aj1p12i= re -it ; aj2p3r = re1-im0; aj2p3i = re0-im1;
				rt = re0 = aj1p14r; it = im0 = aj1p14i; re = re1 = aj2p1r ; im = im1 = aj2p1i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j1pad+28],ab_mul[j1pad+29], ab_mul[j2pad+ 2],ab_mul[j2pad+ 3]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j1pad+28],cd_mul[j1pad+29], cd_mul[j2pad+ 2],cd_mul[j2pad+ 3]);	aj1p14r= rt -im ; aj1p14i= re -it ; aj2p1r = re1-im0; aj2p1i = re0-im1;

				rt = re0 = aj2p0r ; it = im0 = aj2p0i ; re = re1 = aj1p15r; im = im1 = aj1p15i;	PAIR_MUL(rt ,it , re ,im , ab_mul[j2pad+ 0],ab_mul[j2pad+ 1], ab_mul[j1pad+30],ab_mul[j1pad+31]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j2pad+ 0],cd_mul[j2pad+ 1], cd_mul[j1pad+30],cd_mul[j1pad+31]);	aj2p0r = rt -im ; aj2p0i = re -it ; aj1p15r= re1-im0; aj1p15i= re0-im1;
				rt = re0 = aj2p2r ; it = im0 = aj2p2i ; re = re1 = aj1p13r; im = im1 = aj1p13i;	PAIR_MUL(rt ,it , re ,im , ab_mul[j2pad+ 4],ab_mul[j2pad+ 5], ab_mul[j1pad+26],ab_mul[j1pad+27]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j2pad+ 4],cd_mul[j2pad+ 5], cd_mul[j1pad+26],cd_mul[j1pad+27]);	aj2p2r = rt -im ; aj2p2i = re -it ; aj1p13r= re1-im0; aj1p13i= re0-im1;
				rt = re0 = aj2p4r ; it = im0 = aj2p4i ; re = re1 = aj1p11r; im = im1 = aj1p11i;	PAIR_MUL(rt ,it , re ,im , ab_mul[j2pad+ 8],ab_mul[j2pad+ 9], ab_mul[j1pad+22],ab_mul[j1pad+23]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j2pad+ 8],cd_mul[j2pad+ 9], cd_mul[j1pad+22],cd_mul[j1pad+23]);	aj2p4r = rt -im ; aj2p4i = re -it ; aj1p11r= re1-im0; aj1p11i= re0-im1;
				rt = re0 = aj2p6r ; it = im0 = aj2p6i ; re = re1 = aj1p9r ; im = im1 = aj1p9i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j2pad+12],ab_mul[j2pad+13], ab_mul[j1pad+18],ab_mul[j1pad+19]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j2pad+12],cd_mul[j2pad+13], cd_mul[j1pad+18],cd_mul[j1pad+19]);	aj2p6r = rt -im ; aj2p6i = re -it ; aj1p9r = re1-im0; aj1p9i = re0-im1;
				rt = re0 = aj2p8r ; it = im0 = aj2p8i ; re = re1 = aj1p7r ; im = im1 = aj1p7i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j2pad+16],ab_mul[j2pad+17], ab_mul[j1pad+14],ab_mul[j1pad+15]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j2pad+16],cd_mul[j2pad+17], cd_mul[j1pad+14],cd_mul[j1pad+15]);	aj2p8r = rt -im ; aj2p8i = re -it ; aj1p7r = re1-im0; aj1p7i = re0-im1;
				rt = re0 = aj2p10r; it = im0 = aj2p10i; re = re1 = aj1p5r ; im = im1 = aj1p5i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j2pad+20],ab_mul[j2pad+21], ab_mul[j1pad+10],ab_mul[j1pad+11]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j2pad+20],cd_mul[j2pad+21], cd_mul[j1pad+10],cd_mul[j1pad+11]);	aj2p10r= rt -im ; aj2p10i= re -it ; aj1p5r = re1-im0; aj1p5i = re0-im1;
				rt = re0 = aj2p12r; it = im0 = aj2p12i; re = re1 = aj1p3r ; im = im1 = aj1p3i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j2pad+24],ab_mul[j2pad+25], ab_mul[j1pad+ 6],ab_mul[j1pad+ 7]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j2pad+24],cd_mul[j2pad+25], cd_mul[j1pad+ 6],cd_mul[j1pad+ 7]);	aj2p12r= rt -im ; aj2p12i= re -it ; aj1p3r = re1-im0; aj1p3i = re0-im1;
				rt = re0 = aj2p14r; it = im0 = aj2p14i; re = re1 = aj1p1r ; im = im1 = aj1p1i ;	PAIR_MUL(rt ,it , re ,im , ab_mul[j2pad+28],ab_mul[j2pad+29], ab_mul[j1pad+ 2],ab_mul[j1pad+ 3]);	PAIR_MUL(re0,im0, re1,im1, cd_mul[j2pad+28],cd_mul[j2pad+29], cd_mul[j1pad+ 2],cd_mul[j1pad+ 3]);	aj2p14r= rt -im ; aj2p14i= re -it ; aj1p1r = re1-im0; aj1p1i = re0-im1;
			}
		}
	}

/*...And do an inverse DIT radix-16 pass on the squared-data blocks. */

	/*************************************************************/
	/*                  1st set of inputs:                       */
	/*************************************************************/
		rdum = j1pad;
		idum = j1pad+RE_IM_STRIDE;
	#if PFETCH
		addr = &a[rdum+32];
	#endif
	/*   gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do the first set of four length-4 IDIT transforms... */

	/*...Block 1: */
	#if PFETCH
		prefetch_p_doubles(addr);
		addr += 4;
	#endif
		t3 =aj1p0r -aj1p1r;	t4 =aj1p0i -aj1p1i;
		t1 =aj1p0r +aj1p1r;	t2 =aj1p0i +aj1p1i;

		t7 =aj1p2r -aj1p3r;	t8 =aj1p2i -aj1p3i;
		t5 =aj1p2r +aj1p3r;	t6 =aj1p2i +aj1p3i;

		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;

		rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
				t8 =t4 +rt;	t4 =t4 -rt;

	/*...Block 2: */
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
			prefetch_p_doubles(addr);
		#endif
		addr += 4;
	#endif
		t11=aj1p4r -aj1p5r;	t12=aj1p4i-aj1p5i;
		t9 =aj1p4r +aj1p5r;	t10=aj1p4i+aj1p5i;

		t15=aj1p6r-aj1p7r;	t16=aj1p6i-aj1p7i;
		t13=aj1p6r+aj1p7r;	t14=aj1p6i+aj1p7i;

		rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
		it =t14;	t14=t10-it;	t10=t10+it;

		rt =t15;	t15=t11-t16;t11=t11+t16;
					t16=t12+rt;	t12=t12-rt;

	/*...Block 3: */
	#if PFETCH
		prefetch_p_doubles(addr);
		addr += 4;
	#endif
		t19=aj1p8r-aj1p9r;	t20=aj1p8i-aj1p9i;
		t17=aj1p8r+aj1p9r;	t18=aj1p8i+aj1p9i;

		t23=aj1p10r-aj1p11r;	t24=aj1p10i-aj1p11i;
		t21=aj1p10r+aj1p11r;	t22=aj1p10i+aj1p11i;

		rt =t21;	t21=t17-rt;	t17=t17+rt;
		it =t22;	t22=t18-it;	t18=t18+it;

		rt =t23;	t23=t19-t24;t19=t19+t24;
					t24=t20+rt;	t20=t20-rt;

	/*...Block 4: */
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
			prefetch_p_doubles(addr);
		#endif
		addr += 4;
	#endif
		t27=aj1p12r-aj1p13r;	t28=aj1p12i-aj1p13i;
		t25=aj1p12r+aj1p13r;	t26=aj1p12i+aj1p13i;

		t31=aj1p14r-aj1p15r;	t32=aj1p14i-aj1p15i;
		t29=aj1p14r+aj1p15r;	t30=aj1p14i+aj1p15i;

		rt =t29;	t29=t25-rt;	t25=t25+rt;
		it =t30;	t30=t26-it;	t26=t26+it;

		rt =t31;	t31=t27-t32;t27=t27+t32;
					t32=t28+rt;	t28=t28-rt;
	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	!	1, exp(-i* 1*twopi/16) =       ( c,-s), exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 3*twopi/16) =       ( s,-c) (for inputs to transform block 2)
	!	1, exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 4*twopi/16) =       ( 0,-1), exp(-i* 6*twopi/16) = ISRT2*(-1,-1) (for inputs to transform block 3)
	!	1, exp(-i* 3*twopi/16) =       ( s,-c), exp(-i* 6*twopi/16) = ISRT2*(-1,-1), exp(-i* 9*twopi/16) =       (-c, s) (for inputs to transform block 4).
	!  (This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
	!   I.e. do similar as above, except inputs a[j1  +p0:15:1] are replaced by t0:30:2,
	!					   a[j1+1+p0:15:1] are replaced by t1:31:2, and v.v. for outputs,
	!   and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.
	*/
	/*...Block 1: t1,9,17,25 */
	#if PFETCH
		prefetch_p_doubles(addr);
		addr += 4;
	#endif
		rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		a[rdum   ]=t1+t17;				a[idum   ]=t2+t18;
		t1	=t1-t17;					t2	=t2-t18;
		a[rdum+16]=t1 *cA8 +t2 *sA8 ;	a[idum+16]=t2 *cA8 -t1 *sA8;

		rt	=t9 +t26;					it	=t10-t25;
		t9	=t9 -t26;					t10	=t10+t25;
		a[rdum+8 ]=rt *cA4 +it *sA4 ;	a[idum+8 ]=it *cA4 -rt *sA4;
		a[rdum+24]=t9 *cA12+t10*sA12;	a[idum+24]=t10*cA12-t9 *sA12;

/*...Block 3: t5,13,21,29 */
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
			prefetch_p_doubles(addr);
		#endif
		addr += 4;
	#endif
		rt =t13;	t13=t5 -t14;	t5 =t5 +t14;
					t14=t6 +rt;		t6 =t6 -rt;

		rt =(t22+t21)*ISRT2;t22=(t22-t21)*ISRT2;	t21=rt;
		rt =(t29-t30)*ISRT2;it =(t29+t30)*ISRT2;
		t29=t21+rt;	t21=t21-rt;
		t30=t22+it;	t22=t22-it;

		rt	=t5 +t21;					it	=t6 +t22;
		t5	=t5 -t21;					t6	=t6 -t22;
		a[rdum+4 ]=rt *cA2 +it *sA2 ;	a[idum+4 ]=it *cA2 -rt *sA2 ;
		a[rdum+20]=t5 *cA10+t6 *sA10;	a[idum+20]=t6 *cA10-t5 *sA10;

		rt	=t13+t30;					it	=t14-t29;
		t13	=t13-t30;					t14	=t14+t29;
		a[rdum+12]=rt *cA6 +it *sA6 ;	a[idum+12]=it *cA6 -rt *sA6 ;
		a[rdum+28]=t13*cA14+t14*sA14;	a[idum+28]=t14*cA14-t13*sA14;

/*...Block 2: t3,11,19,27 */
	#if PFETCH
		prefetch_p_doubles(addr);
		addr += 4;
	#endif
		rt =(t12+t11)*ISRT2;it =(t12-t11)*ISRT2;
		t11=t3 -rt;	t3 =t3 +rt;
		t12=t4 -it;	t4 =t4 +it;

		rt =t19*c + t20*s;	t20=t20*c - t19*s;	t19=rt;
		rt =t27*s + t28*c;	it =t28*s - t27*c;
		t27=t19-rt;	t19=t19+rt;
		t28=t20-it;	t20=t20+it;

		rt	=t3 +t19;					it	=t4 +t20;
		t3	=t3 -t19;					t4	=t4 -t20;
		a[rdum+2 ]=rt *cA1 +it *sA1 ;	a[idum+2 ]=it *cA1 -rt *sA1 ;
		a[rdum+18]=t3 *cA9 +t4 *sA9 ;	a[idum+18]=t4 *cA9 -t3 *sA9 ;

		rt	=t11+t28;					it	=t12-t27;
		t11	=t11-t28;					t12	=t12+t27;
		a[rdum+10]=rt *cA5 +it *sA5 ;	a[idum+10]=it *cA5 -rt *sA5 ;
		a[rdum+26]=t11*cA13+t12*sA13;	a[idum+26]=t12*cA13-t11*sA13;

/*...Block 4: t7,15,23,31 */

	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
			prefetch_p_doubles(addr);
		#endif
		addr += 4;
	#endif
		rt =(t15-t16)*ISRT2;it =(t15+t16)*ISRT2;
		t15=t7 +rt;	t7 =t7 -rt;
		t16=t8 +it;	t8 =t8 -it;

		rt =t23*s + t24*c;	t24=t24*s - t23*c;	t23=rt;
		rt =t31*c + t32*s;	it =t32*c - t31*s;
		t31=t23+rt;	t23=t23-rt;
		t32=t24+it;	t24=t24-it;

		rt	=t7 +t23;					it	=t8 +t24;
		t7	=t7 -t23;					t8	=t8 -t24;
		a[rdum+6 ]=rt *cA3 +it *sA3 ;	a[idum+6 ]=it *cA3 -rt *sA3 ;
		a[rdum+22]=t7 *cA11+t8 *sA11;	a[idum+22]=t8 *cA11-t7 *sA11;

		rt	=t15+t32;					it	=t16-t31;
		t15	=t15-t32;					t16	=t16+t31;
		a[rdum+14]=rt *cA7 +it *sA7 ;	a[idum+14]=it *cA7 -rt *sA7 ;
		a[rdum+30]=t15*cA15+t16*sA15;	a[idum+30]=t16*cA15-t15*sA15;

	/*************************************************************/
	/*                  2nd set of inputs:                       */
	/*************************************************************/

		rdum = j2pad;
		idum = j2pad+RE_IM_STRIDE;
	#if PFETCH
		addr = &a[rdum-32];
	#endif
	/*...Block 1: */
	#if PFETCH
		prefetch_p_doubles(addr);
		addr += 4;
	#endif
		t3 =aj2p0r -aj2p1r;	t4 =aj2p0i -aj2p1i;
		t1 =aj2p0r +aj2p1r;	t2 =aj2p0i +aj2p1i;

		t7 =aj2p2r -aj2p3r;	t8 =aj2p2i -aj2p3i;
		t5 =aj2p2r +aj2p3r;	t6 =aj2p2i +aj2p3i;

		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;

		rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
			t8 =t4 +rt;	t4 =t4 -rt;

	/*...Block 2: */
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
			prefetch_p_doubles(addr);
		#endif
		addr += 4;
	#endif
		t11=aj2p4r -aj2p5r;	t12=aj2p4i-aj2p5i;
		t9 =aj2p4r +aj2p5r;	t10=aj2p4i+aj2p5i;

		t15=aj2p6r-aj2p7r;	t16=aj2p6i-aj2p7i;
		t13=aj2p6r+aj2p7r;	t14=aj2p6i+aj2p7i;

		rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
		it =t14;	t14=t10-it;	t10=t10+it;

		rt =t15;	t15=t11-t16;	t11=t11+t16;
			t16=t12+rt;	t12=t12-rt;

	/*...Block 3: */
	#if PFETCH
		prefetch_p_doubles(addr);
		addr += 4;
	#endif
		t19=aj2p8r-aj2p9r;	t20=aj2p8i-aj2p9i;
		t17=aj2p8r+aj2p9r;	t18=aj2p8i+aj2p9i;

		t23=aj2p10r-aj2p11r;	t24=aj2p10i-aj2p11i;
		t21=aj2p10r+aj2p11r;	t22=aj2p10i+aj2p11i;

		rt =t21;	t21=t17-rt;	t17=t17+rt;
		it =t22;	t22=t18-it;	t18=t18+it;

		rt =t23;	t23=t19-t24;	t19=t19+t24;
			t24=t20+rt;	t20=t20-rt;

	/*...Block 4: */
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
			prefetch_p_doubles(addr);
		#endif
		addr += 4;
	#endif
		t27=aj2p12r-aj2p13r;	t28=aj2p12i-aj2p13i;
		t25=aj2p12r+aj2p13r;	t26=aj2p12i+aj2p13i;

		t31=aj2p14r-aj2p15r;	t32=aj2p14i-aj2p15i;
		t29=aj2p14r+aj2p15r;	t30=aj2p14i+aj2p15i;

		rt =t29;	t29=t25-rt;	t25=t25+rt;
		it =t30;	t30=t26-it;	t26=t26+it;

		rt =t31;	t31=t27-t32;	t27=t27+t32;
			t32=t28+rt;	t28=t28-rt;

	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	*/
	/*...Block 1: t1,9,17,25 */
	#if PFETCH
		prefetch_p_doubles(addr);
		addr += 4;
	#endif
		rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		a[rdum   ]=t1+t17;				a[idum   ]=t2+t18;
		t1	=t1-t17;					t2	=t2-t18;
		a[rdum+16]=t1 *cB8 +t2 *sB8 ;	a[idum+16]=t2 *cB8 -t1 *sB8;

		rt	=t9 +t26;					it	=t10-t25;
		t9	=t9 -t26;					t10	=t10+t25;
		a[rdum+8 ]=rt *cB4 +it *sB4;	a[idum+8 ]=it *cB4 -rt *sB4;
		a[rdum+24]=t9 *cB12+t10*sB12;	a[idum+24]=t10*cB12-t9 *sB12;

	/*...Block 3: t5,13,21,29 */
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
			prefetch_p_doubles(addr);
		#endif
		addr += 4;
	#endif
		rt =t13;	t13=t5 -t14;	t5 =t5 +t14;
			t14=t6 +rt;	t6 =t6 -rt;

		rt =(t22+t21)*ISRT2;t22=(t22-t21)*ISRT2;	t21=rt;
		rt =(t29-t30)*ISRT2;it =(t29+t30)*ISRT2;
		t29=t21+rt;	t21=t21-rt;
		t30=t22+it;	t22=t22-it;

		rt	=t5 +t21;					it	=t6 +t22;
		t5	=t5 -t21;					t6	=t6 -t22;
		a[rdum+4 ]=rt *cB2 +it *sB2 ;	a[idum+4 ]=it *cB2 -rt *sB2 ;
		a[rdum+20]=t5 *cB10+t6 *sB10;	a[idum+20]=t6 *cB10-t5 *sB10;

		rt	=t13+t30;					it	=t14-t29;
		t13	=t13-t30;					t14	=t14+t29;
		a[rdum+12]=rt *cB6 +it *sB6 ;	a[idum+12]=it *cB6 -rt *sB6 ;
		a[rdum+28]=t13*cB14+t14*sB14;	a[idum+28]=t14*cB14-t13*sB14;

	/*...Block 2: t3,11,19,27 */
	#if PFETCH
		prefetch_p_doubles(addr);
		addr += 4;
	#endif
		rt =(t12+t11)*ISRT2;it =(t12-t11)*ISRT2;
		t11=t3 -rt;	t3 =t3 +rt;
		t12=t4 -it;	t4 =t4 +it;

		rt =t19*c + t20*s;	t20=t20*c - t19*s;	t19=rt;
		rt =t27*s + t28*c;	it =t28*s - t27*c;
		t27=t19-rt;	t19=t19+rt;
		t28=t20-it;	t20=t20+it;

		rt	=t3 +t19;					it	=t4 +t20;
		t3	=t3 -t19;					t4	=t4 -t20;
		a[rdum+2 ]=rt *cB1 +it *sB1 ;	a[idum+2 ]=it *cB1 -rt *sB1 ;
		a[rdum+18]=t3 *cB9 +t4 *sB9 ;	a[idum+18]=t4 *cB9 -t3 *sB9 ;

		rt	=t11+t28;					it	=t12-t27;
		t11	=t11-t28;					t12	=t12+t27;
		a[rdum+10]=rt *cB5 +it *sB5 ;	a[idum+10]=it *cB5 -rt *sB5 ;
		a[rdum+26]=t11*cB13+t12*sB13;	a[idum+26]=t12*cB13-t11*sB13;

	/*...Block 4: t7,15,23,31 */
	#if PFETCH
		#if(CACHE_LINE_DOUBLES == 4)
			prefetch_p_doubles(addr);
		#endif
		addr += 4;
	#endif
		rt =(t15-t16)*ISRT2;it =(t15+t16)*ISRT2;
		t15=t7 +rt;	t7 =t7 -rt;
		t16=t8 +it;	t8 =t8 -it;

		rt =t23*s + t24*c;	t24=t24*s - t23*c;	t23=rt;
		rt =t31*c + t32*s;	it =t32*c - t31*s;
		t31=t23+rt;	t23=t23-rt;
		t32=t24+it;	t24=t24-it;

		rt	=t7 +t23;					it	=t8 +t24;
		t7	=t7 -t23;					t8	=t8 -t24;
		a[rdum+6 ]=rt *cB3 +it *sB3 ;	a[idum+6 ]=it *cB3 -rt *sB3 ;
		a[rdum+22]=t7 *cB11+t8 *sB11;	a[idum+22]=t8 *cB11-t7 *sB11;

		rt	=t15+t32;					it	=t16-t31;
		t15	=t15-t32;					t16	=t16+t31;
		a[rdum+14]=rt *cB7 +it *sB7 ;	a[idum+14]=it *cB7 -rt *sB7 ;
		a[rdum+30]=t15*cB15+t16*sB15;	a[idum+30]=t16*cB15-t15*sB15;

/*...Update the data (j1 and j2) array indices. */
loop:
	  if(j1 <= 64) {
		// Use scalar code (with index offsets properly fiddled) for j1 == 0 case in SIMD mode
		j1 = j1+32;
		j2 = j2-32;
	  } else {
		// Scalar and SSE2 mode both use same increment of 32; avx uses 64:
		j1 = j1+stride;
		j2 = j2-stride;
	  }

	}	/* endfor(m-loop) */

/*
!...Since the foregoing loop only gets executed half as many times as in the simple version, to properly position
!   ourselves in the data array for the start of the next block, need to bump up j1 by as much as would occur in a
!   second execution of the above loop. The exception is the first loop execution, where j1 needs to be doubled (32 x 2).
*/

	j1 = j1+(blocklen << 1);
	if(j2_start == n-32) {
		j1 = 0;
		return;
	}

	// Reset half-complex-blocklength for next pass. If K >> 1 has a zero trailing bit, we multiply
	// blocklength *= (K >> 1) in preparation for the final block:
	blocklen_sum = blocklen_sum + blocklen;
	blocklen = (radix_prim[i-1]-1)*blocklen_sum;

	/*...Next j2_start is previous one plus the (real) length of the current block = 4*(half-complex-blocklength) */
	j2_start = j2_start+(blocklen<<2);
	j2 = j2_start;			    /* Reset j2 for start of the next block. */

}	 /* End of Main (i) loop */

}

