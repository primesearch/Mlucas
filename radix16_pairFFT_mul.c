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

void radix16_pairFFT_mul(double uv[], double ab_mul[], double cd_mul[], int n, int radix0, struct complex rt0[], struct complex rt1[], int ii, int nradices_prim, int radix_prim[], int nloops, int incr, int INIT_ARRAYS, int FORWARD_FFT_ONLY, int skip_square)
{

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

		1) FORWARD_FFT_ONLY = TRUE:
		Compute final radix-16 pass of forward-FFT of A-array only - do not do dyadic mul, do not do IFFT.

		2) FORWARD_FFT_ONLY = FALSE:
		In this case uv[] is assumed to contain data that have just been newly forward-transformed,
		and ab_mul[] and cd_mul[] a pair of precomputed forward-transformed datasets (e.g. corresponding to
		a quartet of a/b/c/d input multiplicands stored in packed even/odd form in this pair of vectors,
		which will be reused multiple times.) In this case, the routine	performs the final radix-16 pass
		of the forward-FFT of the uv-vector, combines the result with the corresponding a/b/c/d-array data
		to obtain FFT(a*u-b*v, c*u-d*v), and performs the initial radix-16 pass of the inverse FFT of these
		linear combinations, returning the result in uv[]. The ab_mul and cd_mul data are unaffected.

		If the cd_mul-array pointer is null, we instead compute the final radix-16 pass
		of the forward-FFT of the uv-vector, combine the result with the corresponding a/b-array data
		to obtain FFT(a*u, b*v), and perform the initial radix-16 pass of the inverse FFT of this pair
		of dyadic products.

	Note that the value of FORWARD_FFT_ONLY is only relevant if INIT_ARRAYS = FALSE.

*/
	double *a = uv;	/* Handy alias pointer to the uv[] vector */

	static int nsave = 0;
/*	static int *index = 0x0;	OBSOLETE: full-N2/16-length Bit-reversal index array. */
	static int *index0 = 0x0, *index1 = 0x0, *index_ptmp0 = 0x0, *index_ptmp1 = 0x0;
	int index0_idx, index0_mod, index1_idx, index1_mod;
	int nradices_prim_radix0;

	int i,j1,j1pad,l,iroot,k1,k2,ndivrad0;
	const double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;	/* exp[i*(twopi/16)] */
	double rt,it
	,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15
	,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32
	,aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r,aj1p6r,aj1p7r,aj1p8r,aj1p9r,aj1p10r,aj1p11r,aj1p12r,aj1p13r,aj1p14r,aj1p15r
	,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i,aj1p6i,aj1p7i,aj1p8i,aj1p9i,aj1p10i,aj1p11i,aj1p12i,aj1p13i,aj1p14i,aj1p15i;
#if PFETCH
	double *add0;
#endif
#if (FFT_DEBUG && 0)
	FILE *dbg_file;
	const char dbg_fname[] = "FFT_DEBUG.txt";
#endif

	index0_mod =        radix0;
	index1_mod = (n>>5)/radix0;	/* complex length requires an additional divide by 2 */

	/***
	Having a separate init block for the big index array allows us to init this prior
	to anything else, using the A-array for scratch space in the call to bit_reverse_int:
	***/
	if(INIT_ARRAYS)
	{
		nsave = n;

		free((void *)index_ptmp0);	index_ptmp0=0x0;	index0=0x0;
		free((void *)index_ptmp1);	index_ptmp1=0x0;	index1=0x0;

		/*
		!...Final forward (DIF) FFT pass sincos data start out in bit-reversed order.
		!   Allocate and initialize an index array containing N/16 indices...

		index_ptmp = ALLOC_INT(N2/16);
		if(!index_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array INDEX in radix16_dyadic_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		index = ALIGN_INT(index_ptmp);
		if(!index){ sprintf(cbuf,"FATAL: unable to allocate array ITMP in radix16_dyadic_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		for(i=0; i < N2/16; i++)
		{
			index[i]=i;
		}
		*/
		index_ptmp0 = ALLOC_INT(index_ptmp0, index0_mod);	index0 = ALIGN_INT(index_ptmp0);
		index_ptmp1 = ALLOC_INT(index_ptmp1, index1_mod);	index1 = ALIGN_INT(index_ptmp1);

		for(i=0; i < index0_mod; i++){index0[i]=       i;}
		for(i=0; i < index1_mod; i++){index1[i]=radix0*i;}

		/*...then bit-reverse INDEX with respect to N/16. Radices are sent to BR routine
		in opposite order from which they are processed in forward FFT.

		bit_reverse_int(index, N2/16, nradices_prim-4, &radix_prim[nradices_prim-5], -1, (int *)a);
		*/

		i = 1;
		for(nradices_prim_radix0 = 1; nradices_prim_radix0 < nradices_prim; nradices_prim_radix0++)
		{
			i *= radix_prim[nradices_prim_radix0-1];
			if(i == radix0)
				break;
		}
		ASSERT(HERE, nradices_prim_radix0 < nradices_prim,"radix16_dyadic_square.c: nradices_prim_radix0 < nradices_prim");

		bit_reverse_int(index0, index0_mod,                 nradices_prim_radix0, &radix_prim[nradices_prim_radix0-1], -1, (int *)a);
		bit_reverse_int(index1, index1_mod, nradices_prim-4-nradices_prim_radix0, &radix_prim[nradices_prim       -5], -1, (int *)a);

		return;
	}

	/*...If a new runlength, should not get to this point: */
	if(n != nsave)
	{
		printf("FATAL: radix16_pairFFT_mul: INIT_ARRAYS not invoked for new runlength!\n");
		ASSERT(HERE, 0,"0");
	}

	/* If precomputing a forward FFT of a set of inputs, make sure
	they're in the uv-vector and the abcd-multiplier vectors are null: */
	if(FORWARD_FFT_ONLY)
	{
		ASSERT(HERE, (ab_mul == 0x0 && cd_mul == 0x0),"radix16_pairFFT_mul: FORWARD_FFT_ONLY = TRUE but non-null abcd-multiplier vectors!");
	}

	ASSERT(HERE, incr == 32,"radix16_pairFFT_mul: incr != 32");
	ndivrad0 = n/radix0;
	/*
	k = ii*(ndivrad0 >> 5);
	*/
	index0_idx = ii;
	index1_idx = 0;

#if (FFT_DEBUG && 0)
	dbg_file = fopen(dbg_fname, "a");
	fprintf(dbg_file, "radix16_pairFFT_mul: jhi = %u, incr = 32:\n",ndivrad0);
#endif

	for(j1 = 0; j1 < ndivrad0; j1 += 32)
	{
	  j1pad = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );
		/*
		iroot = index[k];
		ASSERT(HERE, iroot == index0[index0_idx] + index1[index1_idx],"radix16_pairFFT_mul: iroot == index0[index0_idx] + index1[index1_idx]");
		k = k + 1;	// increment sincos array index
		*/
		iroot = index0[index0_idx] + index1[index1_idx];

		if(++index1_idx >= index1_mod)
		{
			index1_idx -= index1_mod;
		/*	ASSERT(HERE, index0_idx < index0_mod,"radix16_pairFFT.c: index0_idx < index0_mod");	*/
			++index0_idx;
		}

		l = iroot;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;		 /* 2*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		rt=rt1[k2].re;	it=rt1[k2].im;
		c1 =t1*rt-t2*it;	s1 =t1*it+t2*rt;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 1);	/* 4*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		rt=rt1[k2].re;	it=rt1[k2].im;
		c2 =t1*rt-t2*it;	s2 =t1*it+t2*rt;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 2);	/* 8*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		rt=rt1[k2].re;	it=rt1[k2].im;
		c4 =t1*rt-t2*it;	s4 =t1*it+t2*rt;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 2) + iroot;	/* 13*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		rt=rt1[k2].re;	it=rt1[k2].im;
		c8 =t1*rt-t2*it;	s8 =t1*it+t2*rt;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		t1=rt0[k1].re;	t2=rt0[k1].im;
		rt=rt1[k2].re;	it=rt1[k2].im;
		c13=t1*rt-t2*it;	s13=t1*it+t2*rt;

		/* c3,5 */
		t1=c1 *c4 ;	t2=c1 *s4 ;	rt=s1 *c4 ;	it=s1 *s4;
		c3 =t1 +it;	s3 =t2 -rt;	c5 =t1 -it;	s5 =t2 +rt;

		/* c6,7,9,10 */
		t1=c1 *c8 ;	t2=c1 *s8 ;	rt=s1 *c8 ;	it=s1 *s8;
		c7 =t1 +it;	s7 =t2 -rt;	c9 =t1 -it;	s9 =t2 +rt;

		t1=c2 *c8 ;	t2=c2 *s8 ;	rt=s2 *c8 ;	it=s2 *s8;
		c6 =t1 +it;	s6 =t2 -rt;	c10=t1 -it;	s10=t2 +rt;

		/* c11,12,14,15 */
		t1=c1 *c13;	t2=c1 *s13;	rt=s1 *c13;	it=s1 *s13;
		c12=t1 +it;	s12=t2 -rt;	c14=t1 -it;	s14=t2 +rt;

		t1=c2 *c13;	t2=c2 *s13;	rt=s2 *c13;	it=s2 *s13;
		c11=t1 +it;	s11=t2 -rt;	c15=t1 -it;	s15=t2 +rt;

/*   gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do the first set of four length-4 transforms... */

/*...Block 1: */
		t1 =a[j1pad   ];						t2 =a[j1pad+1 ];
		rt =a[j1pad+16]*c8 -a[j1pad+17]*s8 ;	it =a[j1pad+17]*c8 +a[j1pad+16]*s8;
		t3 =t1 -rt;	t1 =t1 +rt;
		t4 =t2 -it;	t2 =t2 +it;

		t5 =a[j1pad+8 ]*c4 -a[j1pad+9 ]*s4;		t6 =a[j1pad+9 ]*c4 +a[j1pad+8 ]*s4;
		rt =a[j1pad+24]*c12-a[j1pad+25]*s12;	it =a[j1pad+25]*c12+a[j1pad+24]*s12;
		t7 =t5 -rt;	t5 =t5 +rt;
		t8 =t6 -it;	t6 =t6 +it;

		rt =t5;	t5 =t1 -rt;			t1 =t1 +rt;
		it =t6;	t6 =t2 -it;			t2 =t2 +it;

		rt =t7;	t7 =t3 +t8;			t3 =t3 -t8;
				t8 =t4 -rt;			t4 =t4 +rt;

/*...Block 2: */
		t9 =a[j1pad+4 ]*c2 -a[j1pad+5 ]*s2 ;	t10=a[j1pad+5 ]*c2 +a[j1pad+4 ]*s2 ;
		rt =a[j1pad+20]*c10-a[j1pad+21]*s10;	it =a[j1pad+21]*c10+a[j1pad+20]*s10;
		t11=t9 -rt;	t9 =t9 +rt;
		t12=t10-it;	t10=t10+it;

		t13=a[j1pad+12]*c6 -a[j1pad+13]*s6 ;	t14=a[j1pad+13]*c6 +a[j1pad+12]*s6 ;
		rt =a[j1pad+28]*c14-a[j1pad+29]*s14;	it =a[j1pad+29]*c14+a[j1pad+28]*s14;
		t15=t13-rt;	t13=t13+rt;
		t16=t14-it;	t14=t14+it;

		rt =t13;	t13=t9 -rt;		t9 =t9 +rt;
		it =t14;	t14=t10-it;		t10=t10+it;

		rt =t15;	t15=t11+t16;	t11=t11-t16;
		t16=t12-rt;	t12=t12+rt;

/*...Block 3: */
		t17=a[j1pad+2 ]*c1 -a[j1pad+3 ]*s1 ;	t18=a[j1pad+3 ]*c1 +a[j1pad+2 ]*s1 ;
		rt =a[j1pad+18]*c9 -a[j1pad+19]*s9 ;	it =a[j1pad+19]*c9 +a[j1pad+18]*s9 ;
		t19=t17-rt;	t17=t17+rt;
		t20=t18-it;	t18=t18+it;

		t21=a[j1pad+10]*c5 -a[j1pad+11]*s5 ;	t22=a[j1pad+11]*c5 +a[j1pad+10]*s5 ;
		rt =a[j1pad+26]*c13-a[j1pad+27]*s13;	it =a[j1pad+27]*c13+a[j1pad+26]*s13;
		t23=t21-rt;	t21=t21+rt;
		t24=t22-it;	t22=t22+it;

		rt =t21;	t21=t17-rt;		t17=t17+rt;
		it =t22;	t22=t18-it;		t18=t18+it;

		rt =t23;	t23=t19+t24;	t19=t19-t24;
					t24=t20-rt;		t20=t20+rt;

/*...Block 4: */
		t25=a[j1pad+6 ]*c3 -a[j1pad+7 ]*s3 ;	t26=a[j1pad+7 ]*c3 +a[j1pad+6 ]*s3 ;
		rt =a[j1pad+22]*c11-a[j1pad+23]*s11;	it =a[j1pad+23]*c11+a[j1pad+22]*s11;
		t27=t25-rt;	t25=t25+rt;
		t28=t26-it;	t26=t26+it;

		t29=a[j1pad+14]*c7 -a[j1pad+15]*s7 ;	t30=a[j1pad+15]*c7 +a[j1pad+14]*s7 ;
		rt =a[j1pad+30]*c15-a[j1pad+31]*s15;	it =a[j1pad+31]*c15+a[j1pad+30]*s15;
		t31=t29-rt;	t29=t29+rt;
		t32=t30-it;	t30=t30+it;

		rt =t29;	t29=t25-rt;		t25=t25+rt;
		it =t30;	t30=t26-it;		t26=t26+it;

		rt =t31;	t31=t27+t32;	t27=t27-t32;
					t32=t28-rt;		t28=t28+rt;
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
		rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		aj1p0r =t1 +t17;	aj1p0i =t2 +t18;
		aj1p1r =t1 -t17;	aj1p1i =t2 -t18;

		aj1p2r =t9 -t26;	aj1p2i =t10+t25;
		aj1p3r =t9 +t26;	aj1p3i =t10-t25;

/*...Block 3: t5,13,21,29 */
		rt =t13;	t13=t5 +t14;	t5 =t5 -t14;
		t14=t6 -rt;	t6 =t6 +rt;

		rt =(t21-t22)*ISRT2;	t22=(t21+t22)*ISRT2;	t21=rt;
		rt =(t30+t29)*ISRT2;	it =(t30-t29)*ISRT2;
		t29=t21+rt;	t21=t21-rt;
		t30=t22+it;	t22=t22-it;

		aj1p4r =t5 +t21;	aj1p4i =t6 +t22;
		aj1p5r =t5 -t21;	aj1p5i =t6 -t22;

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

		aj1p8r =t3 +t19;	aj1p8i =t4 +t20;
		aj1p9r =t3 -t19;	aj1p9i =t4 -t20;

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

		aj1p12r=t7 +t23;	aj1p12i=t8 +t24;
		aj1p13r=t7 -t23;	aj1p13i=t8 -t24;

		aj1p14r=t15-t32;	aj1p14i=t16+t31;
		aj1p15r=t15+t32;	aj1p15i=t16-t31;

	if(FORWARD_FFT_ONLY)
	{
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

		continue;
	}
	else
	{
		if(skip_square)
		{
			/* no-op */
		}
		else if(!cd_mul)
		{
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
			(^ denotes complex-conjugate):
			
				A[n-j] = A^[j], j = 1, ..., n-1. (j = 0 is special case where no conjugate symmetry applies, nor is needed)
			
			Similarly, for a length-n real-vector b[] whose elements are stored as the imaginary parts of a length-n complex
			vector, the outputs satisfy:
			
				B[n-j] = -B^[j], j = 1, ..., n-1. (j = 0 is again a special case).
			
			For length-n real-vectors A and B packed in interleaved [re = A, im = B] fashion into a length-n complex
			vector, linearity of the DFT means that the outputs of the length-n complex fFFT of the packed-input vector -
			call the fFFT result C - satisfy
			
				C[j] = A[j] + B[j]	(linearity), and

				C[n-j] = A^[j] - B^[j]	(linearity) and conjugate symmetry), for j = 1, ..., n-1.
			
			Thus to unpack the separate components of A = DFT(a) and B = DFT(b) we take the conjugate of the 2nd symmetry
				
				A[j] + B[j]	= C[j],

				A[j] - B[j]	= C^[n-j], for j = 1, ..., n-1,
			
			Whence

				A[j] = (C[j] + C^[n-j])/2,
			
				B[j] = (C[j] - C^[n-j])/2, for j = 1, ..., n-1,
				
			and the zero-elements are simply A[0] = Re(C[0]), B[0[ = Im(C[0]) .
			
			*/
		#if 0	/************* This needs a redo: *************************/
			aj1p0r  *= ab_mul[j1pad+ 0];	aj1p0i  *= ab_mul[j1pad+ 1];
			aj1p1r  *= ab_mul[j1pad+ 2];	aj1p1i  *= ab_mul[j1pad+ 3];
			aj1p2r  *= ab_mul[j1pad+ 4];	aj1p2i  *= ab_mul[j1pad+ 5];
			aj1p3r  *= ab_mul[j1pad+ 6];	aj1p3i  *= ab_mul[j1pad+ 7];
			aj1p4r  *= ab_mul[j1pad+ 8];	aj1p4i  *= ab_mul[j1pad+ 9];
			aj1p5r  *= ab_mul[j1pad+10];	aj1p5i  *= ab_mul[j1pad+11];
			aj1p6r  *= ab_mul[j1pad+12];	aj1p6i  *= ab_mul[j1pad+13];
			aj1p7r  *= ab_mul[j1pad+14];	aj1p7i  *= ab_mul[j1pad+15];
			aj1p8r  *= ab_mul[j1pad+16];	aj1p8i  *= ab_mul[j1pad+17];
			aj1p9r  *= ab_mul[j1pad+18];	aj1p9i  *= ab_mul[j1pad+19];
			aj1p10r *= ab_mul[j1pad+20];	aj1p10i *= ab_mul[j1pad+21];
			aj1p11r *= ab_mul[j1pad+22];	aj1p11i *= ab_mul[j1pad+23];
			aj1p12r *= ab_mul[j1pad+24];	aj1p12i *= ab_mul[j1pad+25];
			aj1p13r *= ab_mul[j1pad+26];	aj1p13i *= ab_mul[j1pad+27];
			aj1p14r *= ab_mul[j1pad+28];	aj1p14i *= ab_mul[j1pad+29];
			aj1p15r *= ab_mul[j1pad+30];	aj1p15i *= ab_mul[j1pad+31];
		#else	/**************** Try auto-sqr code from Fermat modmul ****************/
			rt = aj1p0r *aj1p0i ;	aj1p0r  = (aj1p0r  + aj1p0i )*(aj1p0r  - aj1p0i );	aj1p0i  = (rt + rt);
			rt = aj1p1r *aj1p1i ;	aj1p1r  = (aj1p1r  + aj1p1i )*(aj1p1r  - aj1p1i );	aj1p1i  = (rt + rt);
			rt = aj1p2r *aj1p2i ;	aj1p2r  = (aj1p2r  + aj1p2i )*(aj1p2r  - aj1p2i );	aj1p2i  = (rt + rt);
			rt = aj1p3r *aj1p3i ;	aj1p3r  = (aj1p3r  + aj1p3i )*(aj1p3r  - aj1p3i );	aj1p3i  = (rt + rt);
			rt = aj1p4r *aj1p4i ;	aj1p4r  = (aj1p4r  + aj1p4i )*(aj1p4r  - aj1p4i );	aj1p4i  = (rt + rt);
			rt = aj1p5r *aj1p5i ;	aj1p5r  = (aj1p5r  + aj1p5i )*(aj1p5r  - aj1p5i );	aj1p5i  = (rt + rt);
			rt = aj1p6r *aj1p6i ;	aj1p6r  = (aj1p6r  + aj1p6i )*(aj1p6r  - aj1p6i );	aj1p6i  = (rt + rt);
			rt = aj1p7r *aj1p7i ;	aj1p7r  = (aj1p7r  + aj1p7i )*(aj1p7r  - aj1p7i );	aj1p7i  = (rt + rt);
			rt = aj1p8r *aj1p8i ;	aj1p8r  = (aj1p8r  + aj1p8i )*(aj1p8r  - aj1p8i );	aj1p8i  = (rt + rt);
			rt = aj1p9r *aj1p9i ;	aj1p9r  = (aj1p9r  + aj1p9i )*(aj1p9r  - aj1p9i );	aj1p9i  = (rt + rt);
			rt = aj1p10r*aj1p10i;	aj1p10r = (aj1p10r + aj1p10i)*(aj1p10r - aj1p10i);	aj1p10i = (rt + rt);
			rt = aj1p11r*aj1p11i;	aj1p11r = (aj1p11r + aj1p11i)*(aj1p11r - aj1p11i);	aj1p11i = (rt + rt);
			rt = aj1p12r*aj1p12i;	aj1p12r = (aj1p12r + aj1p12i)*(aj1p12r - aj1p12i);	aj1p12i = (rt + rt);
			rt = aj1p13r*aj1p13i;	aj1p13r = (aj1p13r + aj1p13i)*(aj1p13r - aj1p13i);	aj1p13i = (rt + rt);
			rt = aj1p14r*aj1p14i;	aj1p14r = (aj1p14r + aj1p14i)*(aj1p14r - aj1p14i);	aj1p14i = (rt + rt);
			rt = aj1p15r*aj1p15i;	aj1p15r = (aj1p15r + aj1p15i)*(aj1p15r - aj1p15i);	aj1p15i = (rt + rt);
		#endif
		}
		else
		{
			/*
			Dyadic muls of the forward FFT outputs with the corresponding a/b and c/d-vector
			data so as to obtain FFT(a*u-b*v, c*u-d*v).
			Since we assume we're dealing with pairs of FFT vectors (separately stored in the
			even/odd array slots), MULs of Re/Im data are segregated.
			*/
			/* Store u/v in re/im:			u <== a*u - b*v is here:								v <== c*u - d*v is here: */
			rt = aj1p0r ; it = aj1p0i ; aj1p0r  = rt*ab_mul[j1pad+ 0] - it*ab_mul[j1pad+ 1];	aj1p0i  = rt*cd_mul[j1pad+ 0] - it*cd_mul[j1pad+ 1];
			rt = aj1p1r ; it = aj1p1i ; aj1p1r  = rt*ab_mul[j1pad+ 2] - it*ab_mul[j1pad+ 3];	aj1p1i  = rt*cd_mul[j1pad+ 2] - it*cd_mul[j1pad+ 3];
			rt = aj1p2r ; it = aj1p2i ; aj1p2r  = rt*ab_mul[j1pad+ 4] - it*ab_mul[j1pad+ 5];	aj1p2i  = rt*cd_mul[j1pad+ 4] - it*cd_mul[j1pad+ 5];
			rt = aj1p3r ; it = aj1p3i ; aj1p3r  = rt*ab_mul[j1pad+ 6] - it*ab_mul[j1pad+ 7];	aj1p3i  = rt*cd_mul[j1pad+ 6] - it*cd_mul[j1pad+ 7];
			rt = aj1p4r ; it = aj1p4i ; aj1p4r  = rt*ab_mul[j1pad+ 8] - it*ab_mul[j1pad+ 9];	aj1p4i  = rt*cd_mul[j1pad+ 8] - it*cd_mul[j1pad+ 9];
			rt = aj1p5r ; it = aj1p5i ; aj1p5r  = rt*ab_mul[j1pad+10] - it*ab_mul[j1pad+11];	aj1p5i  = rt*cd_mul[j1pad+10] - it*cd_mul[j1pad+11];
			rt = aj1p6r ; it = aj1p6i ; aj1p6r  = rt*ab_mul[j1pad+12] - it*ab_mul[j1pad+13];	aj1p6i  = rt*cd_mul[j1pad+12] - it*cd_mul[j1pad+13];
			rt = aj1p7r ; it = aj1p7i ; aj1p7r  = rt*ab_mul[j1pad+14] - it*ab_mul[j1pad+15];	aj1p7i  = rt*cd_mul[j1pad+14] - it*cd_mul[j1pad+15];
			rt = aj1p8r ; it = aj1p8i ; aj1p8r  = rt*ab_mul[j1pad+16] - it*ab_mul[j1pad+17];	aj1p8i  = rt*cd_mul[j1pad+16] - it*cd_mul[j1pad+17];
			rt = aj1p9r ; it = aj1p9i ; aj1p9r  = rt*ab_mul[j1pad+18] - it*ab_mul[j1pad+19];	aj1p9i  = rt*cd_mul[j1pad+18] - it*cd_mul[j1pad+19];
			rt = aj1p10r; it = aj1p10i; aj1p10r = rt*ab_mul[j1pad+20] - it*ab_mul[j1pad+21];	aj1p10i = rt*cd_mul[j1pad+20] - it*cd_mul[j1pad+21];
			rt = aj1p11r; it = aj1p11i; aj1p11r = rt*ab_mul[j1pad+22] - it*ab_mul[j1pad+23];	aj1p11i = rt*cd_mul[j1pad+22] - it*cd_mul[j1pad+23];
			rt = aj1p12r; it = aj1p12i; aj1p12r = rt*ab_mul[j1pad+24] - it*ab_mul[j1pad+25];	aj1p12i = rt*cd_mul[j1pad+24] - it*cd_mul[j1pad+25];
			rt = aj1p13r; it = aj1p13i; aj1p13r = rt*ab_mul[j1pad+26] - it*ab_mul[j1pad+27];	aj1p13i = rt*cd_mul[j1pad+26] - it*cd_mul[j1pad+27];
			rt = aj1p14r; it = aj1p14i; aj1p14r = rt*ab_mul[j1pad+28] - it*ab_mul[j1pad+29];	aj1p14i = rt*cd_mul[j1pad+28] - it*cd_mul[j1pad+29];
			rt = aj1p15r; it = aj1p15i; aj1p15r = rt*ab_mul[j1pad+30] - it*ab_mul[j1pad+31];	aj1p15i = rt*cd_mul[j1pad+30] - it*cd_mul[j1pad+31];
		}
	}

/*...And do an inverse DIT radix-16 pass on the squared-data blocks. */

/* 1st set of inputs: */
#if PFETCH
add0 = &a[j1pad+32];
#endif
/*   gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do the first set of four length-4 IDIT transforms... */

/*...Block 1: */
#if PFETCH
prefetch_p_doubles(add0);
add0 += 4;
#endif
		t3 =aj1p0r -aj1p1r ;	t4 =aj1p0i -aj1p1i ;
		t1 =aj1p0r +aj1p1r ;	t2 =aj1p0i +aj1p1i ;

		t7 =aj1p2r -aj1p3r ;	t8 =aj1p2i -aj1p3i ;
		t5 =aj1p2r +aj1p3r ;	t6 =aj1p2i +aj1p3i ;

		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;

		rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
		t8 =t4 +rt;	t4 =t4 -rt;

/*...Block 2: */
#if PFETCH
	#if(CACHE_LINE_DOUBLES == 4)
	prefetch_p_doubles(add0);
	#endif
add0 += 4;
#endif
		t11=aj1p4r -aj1p5r ;	t12=aj1p4i -aj1p5i ;
		t9 =aj1p4r +aj1p5r ;	t10=aj1p4i +aj1p5i ;

		t15=aj1p6r -aj1p7r ;	t16=aj1p6i -aj1p7i ;
		t13=aj1p6r +aj1p7r ;	t14=aj1p6i +aj1p7i ;

		rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
		it =t14;	t14=t10-it;	t10=t10+it;

		rt =t15;	t15=t11-t16;	t11=t11+t16;
		t16=t12+rt;	t12=t12-rt;

/*...Block 3: */
#if PFETCH
prefetch_p_doubles(add0);
add0 += 4;
#endif
		t19=aj1p8r -aj1p9r ;	t20=aj1p8i -aj1p9i ;
		t17=aj1p8r +aj1p9r ;	t18=aj1p8i +aj1p9i ;

		t23=aj1p10r-aj1p11r;	t24=aj1p10i-aj1p11i;
		t21=aj1p10r+aj1p11r;	t22=aj1p10i+aj1p11i;

		rt =t21;	t21=t17-rt;	t17=t17+rt;
		it =t22;	t22=t18-it;	t18=t18+it;

		rt =t23;	t23=t19-t24;	t19=t19+t24;
		t24=t20+rt;	t20=t20-rt;

/*...Block 4: */
#if PFETCH
	#if(CACHE_LINE_DOUBLES == 4)
	prefetch_p_doubles(add0);
	#endif
add0 += 4;
#endif
		t27=aj1p12r-aj1p13r;	t28=aj1p12i-aj1p13i;
		t25=aj1p12r+aj1p13r;	t26=aj1p12i+aj1p13i;

		t31=aj1p14r-aj1p15r;	t32=aj1p14i-aj1p15i;
		t29=aj1p14r+aj1p15r;	t30=aj1p14i+aj1p15i;

		rt =t29;	t29=t25-rt;	t25=t25+rt;
		it =t30;	t30=t26-it;	t26=t26+it;

		rt =t31;	t31=t27-t32;	t27=t27+t32;
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
prefetch_p_doubles(add0);
add0 += 4;
#endif
		rt =t9;		t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		a[j1pad   ]=t1+t17;				a[j1pad+1 ]=t2+t18;
		t1	=t1-t17;					t2	=t2-t18;
		a[j1pad+16]=t1 *c8 +t2 *s8 ;	a[j1pad+17]=t2 *c8 -t1 *s8;

		rt	=t9 +t26;					it	=t10-t25;
		t9	=t9 -t26;					t10	=t10+t25;
		a[j1pad+8 ]=rt *c4 +it *s4 ;	a[j1pad+9 ]=it *c4 -rt *s4;
		a[j1pad+24]=t9 *c12+t10*s12;	a[j1pad+25]=t10*c12-t9 *s12;

/*...Block 3: t5,13,21,29 */
#if PFETCH
	#if(CACHE_LINE_DOUBLES == 4)
	prefetch_p_doubles(add0);
	#endif
add0 += 4;
#endif
		rt =t13;	t13=t5 -t14;	t5 =t5 +t14;
					t14=t6 +rt;		t6 =t6 -rt;

		rt =(t22+t21)*ISRT2;	t22=(t22-t21)*ISRT2;	t21=rt;
		rt =(t29-t30)*ISRT2;	it =(t29+t30)*ISRT2;
		t29=t21+rt;	t21=t21-rt;
		t30=t22+it;	t22=t22-it;

		rt	=t5 +t21;					it	=t6 +t22;
		t5	=t5 -t21;					t6	=t6 -t22;
		a[j1pad+4 ]=rt *c2 +it *s2 ;	a[j1pad+5 ]=it *c2 -rt *s2 ;
		a[j1pad+20]=t5 *c10+t6 *s10;	a[j1pad+21]=t6 *c10-t5 *s10;

		rt	=t13+t30;					it	=t14-t29;
		t13	=t13-t30;					t14	=t14+t29;
		a[j1pad+12]=rt *c6 +it *s6 ;	a[j1pad+13]=it *c6 -rt *s6 ;
		a[j1pad+28]=t13*c14+t14*s14;	a[j1pad+29]=t14*c14-t13*s14;

/*...Block 2: t3,11,19,27 */
#if PFETCH
prefetch_p_doubles(add0);
add0 += 4;
#endif
		rt =(t12+t11)*ISRT2;	it =(t12-t11)*ISRT2;
		t11=t3 -rt;	t3 =t3 +rt;
		t12=t4 -it;	t4 =t4 +it;

		rt =t19*c + t20*s;	t20=t20*c - t19*s;	t19=rt;
		rt =t27*s + t28*c;	it =t28*s - t27*c;
		t27=t19-rt;	t19=t19+rt;
		t28=t20-it;	t20=t20+it;

		rt	=t3 +t19;					it	=t4 +t20;
		t3	=t3 -t19;					t4	=t4 -t20;
		a[j1pad+2 ]=rt *c1 +it *s1 ;	a[j1pad+3 ]=it *c1 -rt *s1 ;
		a[j1pad+18]=t3 *c9 +t4 *s9 ;	a[j1pad+19]=t4 *c9 -t3 *s9 ;

		rt	=t11+t28;					it	=t12-t27;
		t11	=t11-t28;					t12	=t12+t27;
		a[j1pad+10]=rt *c5 +it *s5 ;	a[j1pad+11]=it *c5 -rt *s5 ;
		a[j1pad+26]=t11*c13+t12*s13;	a[j1pad+27]=t12*c13-t11*s13;

/*...Block 4: t7,15,23,31 */
#if PFETCH
	#if(CACHE_LINE_DOUBLES == 4)
	prefetch_p_doubles(add0);
	#endif
add0 += 4;
#endif
		rt =(t15-t16)*ISRT2;	it =(t15+t16)*ISRT2;
		rt =(t15-t16)*ISRT2;	it =(t15+t16)*ISRT2;
		t15=t7 +rt;	t7 =t7 -rt;
		t16=t8 +it;	t8 =t8 -it;

		rt =t23*s + t24*c;	t24=t24*s - t23*c;	t23=rt;
		rt =t31*c + t32*s;	it =t32*c - t31*s;
		t31=t23+rt;	t23=t23-rt;
		t32=t24+it;	t24=t24-it;

		rt	=t7 +t23;					it	=t8 +t24;
		t7	=t7 -t23;					t8	=t8 -t24;
		a[j1pad+6 ]=rt *c3 +it *s3 ;	a[j1pad+7 ]=it *c3 -rt *s3 ;
		a[j1pad+22]=t7 *c11+t8 *s11;	a[j1pad+23]=t8 *c11-t7 *s11;

		rt	=t15+t32;					it	=t16-t31;
		t15	=t15-t32;					t16	=t16+t31;
		a[j1pad+14]=rt *c7 +it *s7 ;	a[j1pad+15]=it *c7 -rt *s7 ;
		a[j1pad+30]=t15*c15+t16*s15;	a[j1pad+31]=t16*c15-t15*s15;

#if (FFT_DEBUG && 0)
	fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j1pad+ 0, a[j1pad+ 0], a[j1pad+ 1]);
	fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j1pad+ 2, a[j1pad+ 2], a[j1pad+ 3]);
	fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j1pad+ 4, a[j1pad+ 4], a[j1pad+ 5]);
	fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j1pad+ 6, a[j1pad+ 6], a[j1pad+ 7]);
	fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j1pad+ 8, a[j1pad+ 8], a[j1pad+ 9]);
	fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j1pad+10, a[j1pad+10], a[j1pad+11]);
	fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j1pad+12, a[j1pad+12], a[j1pad+13]);
	fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j1pad+14, a[j1pad+14], a[j1pad+15]);
	fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j1pad+16, a[j1pad+16], a[j1pad+17]);
	fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j1pad+18, a[j1pad+18], a[j1pad+19]);
	fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j1pad+20, a[j1pad+20], a[j1pad+21]);
	fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j1pad+22, a[j1pad+22], a[j1pad+23]);
	fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j1pad+24, a[j1pad+24], a[j1pad+25]);
	fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j1pad+26, a[j1pad+26], a[j1pad+27]);
	fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j1pad+28, a[j1pad+28], a[j1pad+29]);
	fprintf(dbg_file, "%8u : %20.5f  %20.5f\n", j1pad+30, a[j1pad+30], a[j1pad+31]);
	fprintf(dbg_file, "\n");
#endif
	}

#if (FFT_DEBUG && 0)
fclose(dbg_file); dbg_file = 0x0;
#endif
}

