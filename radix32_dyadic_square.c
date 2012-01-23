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

/***************/

#include "Mlucas.h"

void radix32_dyadic_square(double a[], int arr_scratch[], int n, int radix0, struct complex rt0[], struct complex rt1[], int ii, int nradices_prim, int radix_prim[], int nloops, int incr)
{

/*
!   NOTE: In the following commentary, N refers to the COMPLEX vector length (N2 in the code),
!   which is half the real vector length.
!
!...Acronyms:   DIT = Decimation In Time
!               DIF = Decimation In Frequency
!               FFT = Fast Fourier Transform, i.e. a discrete FT over the complex numbers
!
!...Final pass of the DIF forward transform, pointwise squaring and initial pass of
!   the inverse transform, all performed in one fell swoop (or better, one swell loop :)
*/
	static int nsave = 0;
	static int rad0save = 0;
/*	static int *index = 0x0;	OBSOLETE: full-N2/16-length Bit-reversal index array. */
	static int *index0 = 0x0, *index1 = 0x0, *index_ptmp0 = 0x0, *index_ptmp1 = 0x0;
	static int index0_idx, index0_mod, index1_idx, index1_mod;
	int nradices_prim_radix0;

	int i,j,j1,j2,k,l,iroot,k1,k2,ndivrad0;
	const double c = 0.92387953251128675613, s     = 0.38268343236508977173	/* exp[  i*(twopi/16)]	*/
			,c32_1 = 0.98078528040323044912, s32_1 = 0.19509032201612826784	/* exp(  i*twopi/32), the radix-32 fundamental sincos datum	*/
			,c32_3 = 0.83146961230254523708, s32_3 = 0.55557023301960222473;/* exp(3*i*twopi/32)	*/

	double rt,it
			,c01,c02,c03,c04,c05,c06,c07,c08,c09,c0A,c0B,c0C,c0D,c0E,c0F
	,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c1A,c1B,c1C,c1D,c1E,c1F
			,s01,s02,s03,s04,s05,s06,s07,s08,s09,s0A,s0B,s0C,s0D,s0E,s0F
	,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s1A,s1B,s1C,s1D,s1E,s1F
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F
	,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F
	,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F
	,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p0Ar,a1p0Br,a1p0Cr,a1p0Dr,a1p0Er,a1p0Fr
	,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p1Ar,a1p1Br,a1p1Cr,a1p1Dr,a1p1Er,a1p1Fr
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p0Ai,a1p0Bi,a1p0Ci,a1p0Di,a1p0Ei,a1p0Fi
	,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p1Ai,a1p1Bi,a1p1Ci,a1p1Di,a1p1Ei,a1p1Fi;
#if PFETCH
	double *add0;
#endif

	static int first_entry=TRUE;

/*...initialize things upon first entry */
/*...If a new runlength or first-pass radix, set first_entry to true:	*/

	if(n != nsave || radix0 != rad0save) first_entry=TRUE;

	if(first_entry)
	{
		first_entry=FALSE;
		nsave = n;
		rad0save = radix0;

		free((void *)index_ptmp0);	index_ptmp0=0x0;	index0=0x0;
		free((void *)index_ptmp1);	index_ptmp1=0x0;	index1=0x0;

		/*
		!...Final forward (DIF) FFT pass sincos data start out in bit-reversed order.
		!   Allocate and initialize an index array containing N/32 indices...

		index_ptmp = ALLOC_INT(N2/32);
		if(!index_ptmp){ sprintf(cbuf,"FATAL: unable to allocate array INDEX in radix32_dyadic_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		index = ALIGN_INT(index_ptmp);
		if(!index){ sprintf(cbuf,"FATAL: unable to allocate array ITMP in radix32_dyadic_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }

		for(i=0; i < N2/32; i++)
		{
			index[i]=i;
		}
		*/
		index0_mod =        radix0;
		index1_mod = (n>>6)/radix0;	/* complex length requires an additional divide by 2 */

		index_ptmp0 = ALLOC_INT(index_ptmp0, index0_mod);
		if(!index_ptmp0){ sprintf(cbuf,"FATAL: unable to allocate array INDEX_PTMP0 in radix32_dyadic_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		index0 = ALIGN_INT(index_ptmp0);

		index_ptmp1 = ALLOC_INT(index_ptmp1, index1_mod);
		if(!index_ptmp1){ sprintf(cbuf,"FATAL: unable to allocate array INDEX_PTMP1 in radix32_dyadic_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		index1 = ALIGN_INT(index_ptmp1);

		for(i=0; i < index0_mod; i++){index0[i]=       i;}
		for(i=0; i < index1_mod; i++){index1[i]=radix0*i;}

		/*...then bit-reverse INDEX with respect to N/32. Radices are sent to BR routine
		in opposite order from which they are processed in forward FFT.
		bit_reverse_int(index, N2/32, nradices_prim-5, &radix_prim[nradices_prim-6], -1, (int *)arr_scratch);
		*/

		i = 1;
		for(nradices_prim_radix0 = 1; nradices_prim_radix0 < nradices_prim; nradices_prim_radix0++)
		{
			i *= radix_prim[nradices_prim_radix0-1];
			if(i == radix0)
				break;
		}
		ASSERT(HERE, nradices_prim_radix0 < nradices_prim,"radix32_dyadic_square.c: nradices_prim_radix0 < nradices_prim");

		bit_reverse_int(index0, index0_mod,                 nradices_prim_radix0, &radix_prim[nradices_prim_radix0-1], -1, (int *)arr_scratch);
		bit_reverse_int(index1, index1_mod, nradices_prim-5-nradices_prim_radix0, &radix_prim[nradices_prim       -6], -1, (int *)arr_scratch);

		return;
	}

	/*...If a new runlength, should not get to this point: */
	if(n != nsave)
	{
		printf("FATAL: radix32_dyadic_square: INIT_ARRAYS not invoked for new runlength!\n");
		ASSERT(HERE, 0,"0");
	}

	ASSERT(HERE, incr == 64,"radix32_dyadic_square.c: incr == 64");
	ndivrad0 = n/radix0;
	/*
	k = ii*(ndivrad0 >> 6);
	*/
	index0_idx = ii;
	index1_idx = 0;

	for(j = 0; j < ndivrad0; j += 64)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
	    j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
	    j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;
	/*
		iroot = index[k];
		ASSERT(HERE, iroot == index0[index0_idx] + index1[index1_idx],"radix32_dyadic_square.c: iroot == index0[index0_idx] + index1[index1_idx]");
		k = k + 1;	// increment sincos array index
	*/
		iroot = index0[index0_idx] + index1[index1_idx];

		if(++index1_idx >= index1_mod)
		{
			index1_idx -= index1_mod;
		/*	ASSERT(HERE, index0_idx < index0_mod,"radix32_dyadic_square.c: index0_idx < index0_mod");	*/
			++index0_idx;
		}

		l = iroot;

			k1=(l & NRTM1);
			k2=(l >> NRT_BITS);
			l += iroot;			/* 2*iroot	*/
			t00=rt0[k1].re;	t01=rt0[k1].im;
			rt =rt1[k2].re;	it =rt1[k2].im;
			c01=t00*rt -t01*it;	s01=t00*it +t01*rt;

			k1=(l & NRTM1);
			k2=(l >> NRT_BITS);
			l += iroot;				/* 3*iroot	*/
			t00=rt0[k1].re;	t01=rt0[k1].im;
			rt =rt1[k2].re;	it =rt1[k2].im;
			c02=t00*rt -t01*it;	s02=t00*it +t01*rt;

			k1=(l & NRTM1);
			k2=(l >> NRT_BITS);
			l += 4*iroot;			/* 7*iroot	*/
			iroot = l;				/* 7*iroot	*/
			t00=rt0[k1].re;	t01=rt0[k1].im;
			rt =rt1[k2].re;	it =rt1[k2].im;
			c03=t00*rt -t01*it;	s03=t00*it +t01*rt;

			k1=(l & NRTM1);
			k2=(l >> NRT_BITS);
			l += iroot;				/* 14*iroot	*/
			t00=rt0[k1].re;	t01=rt0[k1].im;
			rt =rt1[k2].re;	it =rt1[k2].im;
			c07=t00*rt -t01*it;	s07=t00*it +t01*rt;

			k1=(l & NRTM1);
			k2=(l >> NRT_BITS);
			l += iroot;				/* 21*iroot	*/
			t00=rt0[k1].re;	t01=rt0[k1].im;
			rt =rt1[k2].re;	it =rt1[k2].im;
			c0E=t00*rt -t01*it;	s0E=t00*it +t01*rt;

			k1=(l & NRTM1);
			k2=(l >> NRT_BITS);
			l += iroot;				/* 28*iroot	*/
			t00=rt0[k1].re;	t01=rt0[k1].im;
			rt =rt1[k2].re;	it =rt1[k2].im;
			c15=t00*rt -t01*it;	s15=t00*it +t01*rt;

			k1=(l & NRTM1);
			k2=(l >> NRT_BITS);
			t00=rt0[k1].re;	t01=rt0[k1].im;
			rt =rt1[k2].re;	it =rt1[k2].im;
			c1C=t00*rt -t01*it;	s1C=t00*it +t01*rt;

/*(c,s)4-10 (decimal indices), 4-0A (hexadecimal indices):	*/
			t00=c01*c07; t01=c01*s07; t02=s01*c07; t03=s01*s07;
			c06=t00+t03; s06=t01-t02; c08=t00-t03; s08=t01+t02;

			t00=c02*c07; t01=c02*s07; t02=s02*c07; t03=s02*s07;
			c05=t00+t03; s05=t01-t02; c09=t00-t03; s09=t01+t02;

			t00=c03*c07; t01=c03*s07; t02=s03*c07; t03=s03*s07;
			c04=t00+t03; s04=t01-t02; c0A=t00-t03; s0A=t01+t02;

/*(c,s)11-17 (decimal indices), 0B-11 (hexadecimal indices):;	*/
			t00=c01*c0E; t01=c01*s0E; t02=s01*c0E; t03=s01*s0E;
			c0D=t00+t03; s0D=t01-t02; c0F=t00-t03; s0F=t01+t02;

			t00=c02*c0E; t01=c02*s0E; t02=s02*c0E; t03=s02*s0E;
			c0C=t00+t03; s0C=t01-t02; c10=t00-t03; s10=t01+t02;

			t00=c03*c0E; t01=c03*s0E; t02=s03*c0E; t03=s03*s0E;
			c0B=t00+t03; s0B=t01-t02; c11=t00-t03; s11=t01+t02;

/*(c,s)18-24 (decimal indices), 12-18 (hexadecimal indices):	*/
			t00=c01*c15; t01=c01*s15; t02=s01*c15; t03=s01*s15;
			c14=t00+t03; s14=t01-t02; c16=t00-t03; s16=t01+t02;

			t00=c02*c15; t01=c02*s15; t02=s02*c15; t03=s02*s15;
			c13=t00+t03; s13=t01-t02; c17=t00-t03; s17=t01+t02;

			t00=c03*c15; t01=c03*s15; t02=s03*c15; t03=s03*s15;
			c12=t00+t03; s12=t01-t02; c18=t00-t03; s18=t01+t02;

/*(c,s)25-31 (decimal indices), 19-1F (hexadecimal indices):	*/
			t00=c01*c1C; t01=c01*s1C; t02=s01*c1C; t03=s01*s1C;
			c1B=t00+t03; s1B=t01-t02; c1D=t00-t03; s1D=t01+t02;

			t00=c02*c1C; t01=c02*s1C; t02=s02*c1C; t03=s02*s1C;
			c1A=t00+t03; s1A=t01-t02; c1E=t00-t03; s1E=t01+t02;

			t00=c03*c1C; t01=c03*s1C; t02=s03*c1C; t03=s03*s1C;
			c19=t00+t03; s19=t01-t02; c1F=t00-t03; s1F=t01+t02;

/*       gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...
         We process the sincos data in bit-reversed order.	*/

/*...Block 1:	*/
			t00=a[j1   ];						t01=a[j2     ];
			rt =a[j1+32]*c10 - a[j2+32]*s10;	it =a[j2+32]*c10 + a[j1+32]*s10;
			t02=t00-rt;		t03=t01-it;
			t00=t00+rt;		t01=t01+it;

			t04=a[j1+16]*c08 - a[j2+16]*s08;	t05=a[j2+16]*c08 + a[j1+16]*s08;
			rt =a[j1+48]*c18 - a[j2+48]*s18;	it =a[j2+48]*c18 + a[j1+48]*s18;
			t06=t04-rt;		t07=t05-it;
			t04=t04+rt;		t05=t05+it;

			rt =t04;		it =t05;
			t04=t00-rt;		t05=t01-it;
			t00=t00+rt;		t01=t01+it;

			rt =t06;		it =t07;
			t06=t02+it;		t07=t03-rt;
			t02=t02-it;		t03=t03+rt;

			t08=a[j1+8 ]*c04 - a[j2+8 ]*s04;	t09=a[j2+8 ]*c04 + a[j1+8 ]*s04;
			rt =a[j1+40]*c14 - a[j2+40]*s14;	it =a[j2+40]*c14 + a[j1+40]*s14;
			t0A=t08-rt;		t0B=t09-it;
			t08=t08+rt;		t09=t09+it;

			t0C=a[j1+24]*c0C - a[j2+24]*s0C;	t0D=a[j2+24]*c0C + a[j1+24]*s0C;
			rt =a[j1+56]*c1C - a[j2+56]*s1C;	it =a[j2+56]*c1C + a[j1+56]*s1C;
			t0E=t0C-rt;		t0F=t0D-it;
			t0C=t0C+rt;		t0D=t0D+it;

			rt =t0C;		it =t0D;
			t0C=t08-rt;		t0D=t09-it;
			t08=t08+rt;		t09=t09+it;

			rt =t0E;		it =t0F;
			t0E=t0A+it;		t0F=t0B-rt;
			t0A=t0A-it;		t0B=t0B+rt;

			rt =t08;		it =t09;
			t08=t00-rt;		t09=t01-it;
			t00=t00+rt;		t01=t01+it;

			rt =t0C;		it =t0D;
			t0C=t04+it;		t0D=t05-rt;
			t04=t04-it;		t05=t05+rt;

			rt =(t0A-t0B)*ISRT2;it =(t0A+t0B)*ISRT2;
			t0A=t02-rt;		t0B=t03-it;
			t02=t02+rt;		t03=t03+it;

			rt =(t0E+t0F)*ISRT2;it =(t0F-t0E)*ISRT2;
			t0E=t06+rt;		t0F=t07+it;
			t06=t06-rt;		t07=t07-it;

/*...Block 2:	*/
			t10=a[j1+4 ]*c02 - a[j2+4 ]*s02;	t11=a[j2+4 ]*c02 + a[j1+4 ]*s02;
			rt =a[j1+36]*c12 - a[j2+36]*s12;	it =a[j2+36]*c12 + a[j1+36]*s12;
			t12=t10-rt;		t13=t11-it;
			t10=t10+rt;		t11=t11+it;

			t14=a[j1+20]*c0A - a[j2+20]*s0A;	t15=a[j2+20]*c0A + a[j1+20]*s0A;
			rt =a[j1+52]*c1A - a[j2+52]*s1A;	it =a[j2+52]*c1A + a[j1+52]*s1A;
			t16=t14-rt;		t17=t15-it;
			t14=t14+rt;		t15=t15+it;

			rt =t14;		it =t15;
			t14=t10-rt;		t15=t11-it;
			t10=t10+rt;		t11=t11+it;

			rt =t16;		it =t17;
			t16=t12+it;		t17=t13-rt;
			t12=t12-it;		t13=t13+rt;

			t18=a[j1+12]*c06 - a[j2+12]*s06;	t19=a[j2+12]*c06 + a[j1+12]*s06;
			rt =a[j1+44]*c16 - a[j2+44]*s16;	it =a[j2+44]*c16 + a[j1+44]*s16;
			t1A=t18-rt;		t1B=t19-it;
			t18=t18+rt;		t19=t19+it;

			t1C=a[j1+28]*c0E - a[j2+28]*s0E;	t1D=a[j2+28]*c0E + a[j1+28]*s0E;
			rt =a[j1+60]*c1E - a[j2+60]*s1E;	it =a[j2+60]*c1E + a[j1+60]*s1E;
			t1E=t1C-rt;		t1F=t1D-it;
			t1C=t1C+rt;		t1D=t1D+it;

			rt =t1C;		it =t1D;
			t1C=t18-rt;		t1D=t19-it;
			t18=t18+rt;		t19=t19+it;

			rt =t1E;		it =t1F;
			t1E=t1A+it;		t1F=t1B-rt;
			t1A=t1A-it;		t1B=t1B+rt;

			rt =t18;		it =t19;
			t18=t10-rt;		t19=t11-it;
			t10=t10+rt;		t11=t11+it;

			rt =t1C;		it =t1D;
			t1C=t14+it;		t1D=t15-rt;
			t14=t14-it;		t15=t15+rt;

			rt =(t1A-t1B)*ISRT2;it =(t1A+t1B)*ISRT2;
			t1A=t12-rt;		t1B=t13-it;
			t12=t12+rt;		t13=t13+it;

			rt =(t1E+t1F)*ISRT2;it =(t1F-t1E)*ISRT2;
			t1E=t16+rt;		t1F=t17+it;
			t16=t16-rt;		t17=t17-it;

/*...Block 3:	*/
			t20=a[j1+2 ]*c01 - a[j2+2 ]*s01;	t21=a[j2+2 ]*c01 + a[j1+2 ]*s01;
			rt =a[j1+34]*c11 - a[j2+34]*s11;	it =a[j2+34]*c11 + a[j1+34]*s11;
			t22=t20-rt;		t23=t21-it;
			t20=t20+rt;		t21=t21+it;

			t24=a[j1+18]*c09 - a[j2+18]*s09;	t25=a[j2+18]*c09 + a[j1+18]*s09;
			rt =a[j1+50]*c19 - a[j2+50]*s19;	it =a[j2+50]*c19 + a[j1+50]*s19;
			t26=t24-rt;		t27=t25-it;
			t24=t24+rt;		t25=t25+it;

			rt =t24;		it =t25;
			t24=t20-rt;		t25=t21-it;
			t20=t20+rt;		t21=t21+it;

			rt =t26;		it =t27;
			t26=t22+it;		t27=t23-rt;
			t22=t22-it;		t23=t23+rt;

			t28=a[j1+10]*c05 - a[j2+10]*s05;	t29=a[j2+10]*c05 + a[j1+10]*s05;
			rt =a[j1+42]*c15 - a[j2+42]*s15;	it =a[j2+42]*c15 + a[j1+42]*s15;
			t2A=t28-rt;		t2B=t29-it;
			t28=t28+rt;		t29=t29+it;

			t2C=a[j1+26]*c0D - a[j2+26]*s0D;	t2D=a[j2+26]*c0D + a[j1+26]*s0D;
			rt =a[j1+58]*c1D - a[j2+58]*s1D;	it =a[j2+58]*c1D + a[j1+58]*s1D;
			t2E=t2C-rt;		t2F=t2D-it;
			t2C=t2C+rt;		t2D=t2D+it;

			rt =t2C;		it =t2D;
			t2C=t28-rt;		t2D=t29-it;
			t28=t28+rt;		t29=t29+it;

			rt =t2E;		it =t2F;
			t2E=t2A+it;		t2F=t2B-rt;
			t2A=t2A-it;		t2B=t2B+rt;

			rt =t28;		it =t29;
			t28=t20-rt;		t29=t21-it;
			t20=t20+rt;		t21=t21+it;

			rt =t2C;		it =t2D;
			t2C=t24+it;		t2D=t25-rt;
			t24=t24-it;		t25=t25+rt;

			rt =(t2A-t2B)*ISRT2;it =(t2A+t2B)*ISRT2;
			t2A=t22-rt;		t2B=t23-it;
			t22=t22+rt;		t23=t23+it;

			rt =(t2E+t2F)*ISRT2;it =(t2F-t2E)*ISRT2;
			t2E=t26+rt;		t2F=t27+it;
			t26=t26-rt;		t27=t27-it;

/*...Block 4:	*/
			t30=a[j1+6 ]*c03 - a[j2+6 ]*s03;	t31=a[j2+6 ]*c03 + a[j1+6 ]*s03;
			rt =a[j1+38]*c13 - a[j2+38]*s13;	it =a[j2+38]*c13 + a[j1+38]*s13;
			t32=t30-rt;		t33=t31-it;
			t30=t30+rt;		t31=t31+it;

			t34=a[j1+22]*c0B - a[j2+22]*s0B;	t35=a[j2+22]*c0B + a[j1+22]*s0B;
			rt =a[j1+54]*c1B - a[j2+54]*s1B;	it =a[j2+54]*c1B + a[j1+54]*s1B;
			t36=t34-rt;		t37=t35-it;
			t34=t34+rt;		t35=t35+it;

			rt =t34;		it =t35;
			t34=t30-rt;		t35=t31-it;
			t30=t30+rt;		t31=t31+it;

			rt =t36;		it =t37;
			t36=t32+it;		t37=t33-rt;
			t32=t32-it;		t33=t33+rt;

			t38=a[j1+14]*c07 - a[j2+14]*s07;	t39=a[j2+14]*c07 + a[j1+14]*s07;
			rt =a[j1+46]*c17 - a[j2+46]*s17;	it =a[j2+46]*c17 + a[j1+46]*s17;
			t3A=t38-rt;		t3B=t39-it;
			t38=t38+rt;		t39=t39+it;

			t3C=a[j1+30]*c0F - a[j2+30]*s0F;	t3D=a[j2+30]*c0F + a[j1+30]*s0F;
			rt =a[j1+62]*c1F - a[j2+62]*s1F;	it =a[j2+62]*c1F + a[j1+62]*s1F;
			t3E=t3C-rt;		t3F=t3D-it;
			t3C=t3C+rt;		t3D=t3D+it;

			rt =t3C;		it =t3D;
			t3C=t38-rt;		t3D=t39-it;
			t38=t38+rt;		t39=t39+it;

			rt =t3E;		it =t3F;
			t3E=t3A+it;		t3F=t3B-rt;
			t3A=t3A-it;		t3B=t3B+rt;

			rt =t38;		it =t39;
			t38=t30-rt;		t39=t31-it;
			t30=t30+rt;		t31=t31+it;

			rt =t3C;		it =t3D;
			t3C=t34+it;		t3D=t35-rt;
			t34=t34-it;		t35=t35+rt;

			rt =(t3A-t3B)*ISRT2;it =(t3A+t3B)*ISRT2;
			t3A=t32-rt;		t3B=t33-it;
			t32=t32+rt;		t33=t33+it;

			rt =(t3E+t3F)*ISRT2;it =(t3F-t3E)*ISRT2;
			t3E=t36+rt;		t3F=t37+it;
			t36=t36-rt;		t37=t37-it;

/*...and now do eight radix-4 transforms, including the internal twiddle factors:
	1, exp(i* 1*twopi/32) =       ( c32_1, s32_1), exp(i* 2*twopi/32) =       ( c    , s    ), exp(i* 3*twopi/32) =       ( c32_3, s32_3) (for inputs to transform block 2),
	1, exp(i* 2*twopi/32) =       ( c    , s    ), exp(i* 4*twopi/32) = ISRT2*( 1    , 1    ), exp(i* 3*twopi/32) =       ( s    , c    ) (for inputs to transform block 3),
	1, exp(i* 3*twopi/32) =       ( c32_3, s32_3), exp(i* 6*twopi/32) =       ( s    , c    ), exp(i* 9*twopi/32) =       (-s32_1, c32_1) (for inputs to transform block 4),
	1, exp(i* 4*twopi/32) = ISRT2*( 1    , 1    ), exp(i* 8*twopi/32) =       ( 0    , 1    ), exp(i*12*twopi/32) = ISRT2*(-1    , 1    ) (for inputs to transform block 5),
	1, exp(i* 5*twopi/32) =       ( s32_3, c32_3), exp(i*10*twopi/32) =       (-s    , c    ), exp(i*15*twopi/32) =       (-c32_1, s32_1) (for inputs to transform block 6),
	1, exp(i* 6*twopi/32) =       ( s    , c    ), exp(i*12*twopi/32) = ISRT2*(-1    , 1    ), exp(i*18*twopi/32) =       (-c    ,-s    ) (for inputs to transform block 7),
	1, exp(i* 7*twopi/32) =       ( s32_1, c32_1), exp(i*14*twopi/32) =       (-c    , s    ), exp(i*21*twopi/32) =       (-s32_3,-c32_3) (for inputs to transform block 8),
   and only the last 3 inputs to each of the radix-4 transforms 2 through 8 are multiplied by non-unity twiddles.
   Within each block we process the 8 needed twiddles in bit-reversed order.	*/

/*...Block 1: t00,t10,t20,t30	*/
			rt =t10;	t10=t00-rt;	t00=t00+rt;
			it =t11;	t11=t01-it;	t01=t01+it;

			rt =t30;	t30=t20-rt;	t20=t20+rt;
			it =t31;	t31=t21-it;	t21=t21+it;

			a1p00r = t00+t20;		a1p00i = t01+t21;
			a1p01r = t00-t20;		a1p01i = t01-t21;

			a1p02r = t10-t31;		a1p02i = t11+t30;
			a1p03r = t10+t31;		a1p03i = t11-t30;

/*...Block 5: t08,t18,t28,t38	*/
			rt =t18;	t18=t08+t19;	t08=t08-t19;
			t19=t09-rt;	t09=t09+rt;

			rt =(t28-t29)*ISRT2;	t29=(t28+t29)*ISRT2;		t28=rt;
			rt =(t39+t38)*ISRT2;	it =(t39-t38)*ISRT2;
			t38=t28+rt;			t28=t28-rt;
			t39=t29+it;			t29=t29-it;

			a1p04r = t08+t28;		a1p04i = t09+t29;
			a1p05r = t08-t28;		a1p05i = t09-t29;

			a1p06r = t18-t39;		a1p06i = t19+t38;
			a1p07r = t18+t39;		a1p07i = t19-t38;

/*...Block 3: t04,t14,t24,t34	*/
			rt =(t14-t15)*ISRT2;	it =(t14+t15)*ISRT2;
			t14=t04-rt;			t04=t04+rt;
			t15=t05-it;			t05=t05+it;

			rt =t24*c - t25*s;		t25=t25*c + t24*s;		t24=rt;
			rt =t34*s - t35*c;		it =t35*s + t34*c;
			t34=t24-rt;			t24=t24+rt;
			t35=t25-it;			t25=t25+it;

			a1p08r = t04+t24;		a1p08i = t05+t25;
			a1p09r = t04-t24;		a1p09i = t05-t25;

			a1p0Ar = t14-t35;		a1p0Ai = t15+t34;
			a1p0Br = t14+t35;		a1p0Bi = t15-t34;

/*...Block 7: t0C,t1C,t2C,t3C	*/
			rt =(t1D+t1C)*ISRT2;	it =(t1D-t1C)*ISRT2;
			t1C=t0C+rt;			t0C=t0C-rt;
			t1D=t0D+it;			t0D=t0D-it;

			rt =t2C*s - t2D*c;		t2D=t2D*s + t2C*c;		t2C=rt;
			rt =t3C*c - t3D*s;		it =t3D*c + t3C*s;
			t3C=t2C+rt;			t2C=t2C-rt;
			t3D=t2D+it;			t2D=t2D-it;

			a1p0Cr = t0C+t2C;		a1p0Ci = t0D+t2D;
			a1p0Dr = t0C-t2C;		a1p0Di = t0D-t2D;

			a1p0Er = t1C-t3D;		a1p0Ei = t1D+t3C;
			a1p0Fr = t1C+t3D;		a1p0Fi = t1D-t3C;

/*...Block 2: t02,t12,t22,t32	*/
			rt =t12*c - t13*s;		it =t13*c + t12*s;
			t12=t02-rt;			t02=t02+rt;
			t13=t03-it;			t03=t03+it;

			rt =t22*c32_1 - t23*s32_1;	t23=t23*c32_1 + t22*s32_1;	t22=rt;
			rt =t32*c32_3 - t33*s32_3;	it =t33*c32_3 + t32*s32_3;
			t32=t22-rt;			t22=t22+rt;
			t33=t23-it;			t23=t23+it;

			a1p10r = t02+t22;		a1p10i = t03+t23;
			a1p11r = t02-t22;		a1p11i = t03-t23;

			a1p12r = t12-t33;		a1p12i = t13+t32;
			a1p13r = t12+t33;		a1p13i = t13-t32;

/*...Block 6: t0A,t1A,t2A,t3A	*/
			rt =t1A*s + t1B*c;		it =t1B*s - t1A*c;
			t1A=t0A+rt;			t0A =t0A-rt;
			t1B=t0B+it;			t0B =t0B-it;

			rt =t2A*s32_3 - t2B*c32_3;	t2B=t2B*s32_3 + t2A*c32_3;	t2A=rt;
			rt =t3A*c32_1 + t3B*s32_1;	it =t3B*c32_1 - t3A*s32_1;
			t3A=t2A+rt;			t2A=t2A-rt;
			t3B=t2B+it;			t2B=t2B-it;

			a1p14r = t0A+t2A;		a1p14i = t0B+t2B;
			a1p15r = t0A-t2A;		a1p15i = t0B-t2B;

			a1p16r = t1A-t3B;		a1p16i = t1B+t3A;
			a1p17r = t1A+t3B;		a1p17i = t1B-t3A;

/*...Block 4: t06,t16,t26,t36	*/
			rt =t16*s - t17*c;		it =t17*s + t16*c;
			t16=t06-rt;			t06 =t06+rt;
			t17=t07-it;			t07 =t07+it;

			rt =t26*c32_3 - t27*s32_3;	t27=t27*c32_3 + t26*s32_3;	t26=rt;
			rt =t36*s32_1 + t37*c32_1;	it =t37*s32_1 - t36*c32_1;
			t36=t26+rt;			t26=t26-rt;
			t37=t27+it;			t27=t27-it;

			a1p18r = t06+t26;		a1p18i = t07+t27;
			a1p19r = t06-t26;		a1p19i = t07-t27;

			a1p1Ar = t16-t37;		a1p1Ai = t17+t36;
			a1p1Br = t16+t37;		a1p1Bi = t17-t36;

/*...Block 8: t0E,t1E,t2E,t3E	*/
			rt =t1E*c + t1F*s;		it =t1F*c - t1E*s;
			t1E=t0E+rt;			t0E =t0E-rt;
			t1F=t0F+it;			t0F =t0F-it;

			rt =t2E*s32_1 - t2F*c32_1;	t2F=t2F*s32_1 + t2E*c32_1;	t2E=rt;
			rt =t3E*s32_3 - t3F*c32_3;	it =t3F*s32_3 + t3E*c32_3;
			t3E=t2E+rt;			t2E=t2E-rt;
			t3F=t2F+it;			t2F=t2F-it;

			a1p1Cr = t0E+t2E;		a1p1Ci = t0F+t2F;
			a1p1Dr = t0E-t2E;		a1p1Di = t0F-t2F;

			a1p1Er = t1E-t3F;		a1p1Ei = t1F+t3E;
			a1p1Fr = t1E+t3F;		a1p1Fi = t1F-t3E;

/*...Dyadic square of the forward FFT outputs: */

	rt = a1p00r*a1p00i;	a1p00r = (a1p00r + a1p00i)*(a1p00r - a1p00i);	a1p00i = (rt + rt);
	rt = a1p01r*a1p01i;	a1p01r = (a1p01r + a1p01i)*(a1p01r - a1p01i);	a1p01i = (rt + rt);
	rt = a1p02r*a1p02i;	a1p02r = (a1p02r + a1p02i)*(a1p02r - a1p02i);	a1p02i = (rt + rt);
	rt = a1p03r*a1p03i;	a1p03r = (a1p03r + a1p03i)*(a1p03r - a1p03i);	a1p03i = (rt + rt);
	rt = a1p04r*a1p04i;	a1p04r = (a1p04r + a1p04i)*(a1p04r - a1p04i);	a1p04i = (rt + rt);
	rt = a1p05r*a1p05i;	a1p05r = (a1p05r + a1p05i)*(a1p05r - a1p05i);	a1p05i = (rt + rt);
	rt = a1p06r*a1p06i;	a1p06r = (a1p06r + a1p06i)*(a1p06r - a1p06i);	a1p06i = (rt + rt);
	rt = a1p07r*a1p07i;	a1p07r = (a1p07r + a1p07i)*(a1p07r - a1p07i);	a1p07i = (rt + rt);
	rt = a1p08r*a1p08i;	a1p08r = (a1p08r + a1p08i)*(a1p08r - a1p08i);	a1p08i = (rt + rt);
	rt = a1p09r*a1p09i;	a1p09r = (a1p09r + a1p09i)*(a1p09r - a1p09i);	a1p09i = (rt + rt);
	rt = a1p0Ar*a1p0Ai;	a1p0Ar = (a1p0Ar + a1p0Ai)*(a1p0Ar - a1p0Ai);	a1p0Ai = (rt + rt);
	rt = a1p0Br*a1p0Bi;	a1p0Br = (a1p0Br + a1p0Bi)*(a1p0Br - a1p0Bi);	a1p0Bi = (rt + rt);
	rt = a1p0Cr*a1p0Ci;	a1p0Cr = (a1p0Cr + a1p0Ci)*(a1p0Cr - a1p0Ci);	a1p0Ci = (rt + rt);
	rt = a1p0Dr*a1p0Di;	a1p0Dr = (a1p0Dr + a1p0Di)*(a1p0Dr - a1p0Di);	a1p0Di = (rt + rt);
	rt = a1p0Er*a1p0Ei;	a1p0Er = (a1p0Er + a1p0Ei)*(a1p0Er - a1p0Ei);	a1p0Ei = (rt + rt);
	rt = a1p0Fr*a1p0Fi;	a1p0Fr = (a1p0Fr + a1p0Fi)*(a1p0Fr - a1p0Fi);	a1p0Fi = (rt + rt);

	rt = a1p10r*a1p10i;	a1p10r = (a1p10r + a1p10i)*(a1p10r - a1p10i);	a1p10i = (rt + rt);
	rt = a1p11r*a1p11i;	a1p11r = (a1p11r + a1p11i)*(a1p11r - a1p11i);	a1p11i = (rt + rt);
	rt = a1p12r*a1p12i;	a1p12r = (a1p12r + a1p12i)*(a1p12r - a1p12i);	a1p12i = (rt + rt);
	rt = a1p13r*a1p13i;	a1p13r = (a1p13r + a1p13i)*(a1p13r - a1p13i);	a1p13i = (rt + rt);
	rt = a1p14r*a1p14i;	a1p14r = (a1p14r + a1p14i)*(a1p14r - a1p14i);	a1p14i = (rt + rt);
	rt = a1p15r*a1p15i;	a1p15r = (a1p15r + a1p15i)*(a1p15r - a1p15i);	a1p15i = (rt + rt);
	rt = a1p16r*a1p16i;	a1p16r = (a1p16r + a1p16i)*(a1p16r - a1p16i);	a1p16i = (rt + rt);
	rt = a1p17r*a1p17i;	a1p17r = (a1p17r + a1p17i)*(a1p17r - a1p17i);	a1p17i = (rt + rt);
	rt = a1p18r*a1p18i;	a1p18r = (a1p18r + a1p18i)*(a1p18r - a1p18i);	a1p18i = (rt + rt);
	rt = a1p19r*a1p19i;	a1p19r = (a1p19r + a1p19i)*(a1p19r - a1p19i);	a1p19i = (rt + rt);
	rt = a1p1Ar*a1p1Ai;	a1p1Ar = (a1p1Ar + a1p1Ai)*(a1p1Ar - a1p1Ai);	a1p1Ai = (rt + rt);
	rt = a1p1Br*a1p1Bi;	a1p1Br = (a1p1Br + a1p1Bi)*(a1p1Br - a1p1Bi);	a1p1Bi = (rt + rt);
	rt = a1p1Cr*a1p1Ci;	a1p1Cr = (a1p1Cr + a1p1Ci)*(a1p1Cr - a1p1Ci);	a1p1Ci = (rt + rt);
	rt = a1p1Dr*a1p1Di;	a1p1Dr = (a1p1Dr + a1p1Di)*(a1p1Dr - a1p1Di);	a1p1Di = (rt + rt);
	rt = a1p1Er*a1p1Ei;	a1p1Er = (a1p1Er + a1p1Ei)*(a1p1Er - a1p1Ei);	a1p1Ei = (rt + rt);
	rt = a1p1Fr*a1p1Fi;	a1p1Fr = (a1p1Fr + a1p1Fi)*(a1p1Fr - a1p1Fi);	a1p1Fi = (rt + rt);

/*...And do an inverse DIT radix-32 pass on the squared-data blocks.	*/

/* 1st set of inputs: */
#if PFETCH
add0 = &a[j1+64];
#endif
/*   gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 IDIT transforms... */

/*...Block 1: */
#if PFETCH
prefetch_p_doubles(add0);
add0 += 4;
#endif
			t02=a1p00r-a1p01r;	t03=a1p00i-a1p01i;
			t00=a1p00r+a1p01r;	t01=a1p00i+a1p01i;

			t06=a1p02r-a1p03r;	t07=a1p02i-a1p03i;
			t04=a1p02r+a1p03r;	t05=a1p02i+a1p03i;

			rt =t04;		it =t05;
			t04=t00-rt;		t05=t01-it;
			t00=t00+rt;		t01=t01+it;

			rt =t06;		it =t07;
			t06=t02-it;		t07=t03+rt;
			t02=t02+it;		t03=t03-rt;

			t0A=a1p04r-a1p05r;	t0B=a1p04i-a1p05i;
			t08=a1p04r+a1p05r;	t09=a1p04i+a1p05i;

			t0E=a1p06r-a1p07r;	t0F=a1p06i-a1p07i;
			t0C=a1p06r+a1p07r;	t0D=a1p06i+a1p07i;

			rt =t0C;		it =t0D;
			t0C=t08-rt;		t0D=t09-it;
			t08=t08+rt;		t09=t09+it;

			rt =t0E;		it =t0F;
			t0E=t0A-it;		t0F=t0B+rt;
			t0A=t0A+it;		t0B=t0B-rt;
#if PFETCH
	#if(CACHE_LINE_DOUBLES == 4)
	prefetch_p_doubles(add0);
	#endif
add0 += 4;
#endif
			rt =t08;		it =t09;
			t08=t00-rt;		t09=t01-it;
			t00=t00+rt;		t01=t01+it;

			rt =t0C;		it =t0D;
			t0C=t04-it;		t0D=t05+rt;
			t04=t04+it;		t05=t05-rt;

			rt =(t0A+t0B)*ISRT2;it =(t0A-t0B)*ISRT2;
			t0A=t02-rt;		t0B=t03+it;
			t02=t02+rt;		t03=t03-it;

			rt =(t0E-t0F)*ISRT2;it =(t0F+t0E)*ISRT2;
			t0E=t06+rt;		t0F=t07+it;
			t06=t06-rt;		t07=t07-it;

/*...Block 2:	*/
#if PFETCH
prefetch_p_doubles(add0);
add0 += 4;
#endif
			t12=a1p08r-a1p09r;	t13=a1p08i-a1p09i;
			t10=a1p08r+a1p09r;	t11=a1p08i+a1p09i;

			t16=a1p0Ar-a1p0Br;	t17=a1p0Ai-a1p0Bi;
			t14=a1p0Ar+a1p0Br;	t15=a1p0Ai+a1p0Bi;

			rt =t14;		it =t15;
			t14=t10-rt;		t15=t11-it;
			t10=t10+rt;		t11=t11+it;

			rt =t16;		it =t17;
			t16=t12-it;		t17=t13+rt;
			t12=t12+it;		t13=t13-rt;

			t1A=a1p0Cr-a1p0Dr;	t1B=a1p0Ci-a1p0Di;
			t18=a1p0Cr+a1p0Dr;	t19=a1p0Ci+a1p0Di;

			t1E=a1p0Er-a1p0Fr;	t1F=a1p0Ei-a1p0Fi;
			t1C=a1p0Er+a1p0Fr;	t1D=a1p0Ei+a1p0Fi;

			rt =t1C;		it =t1D;
			t1C=t18-rt;		t1D=t19-it;
			t18=t18+rt;		t19=t19+it;

			rt =t1E;		it =t1F;
			t1E=t1A-it;		t1F=t1B+rt;
			t1A=t1A+it;		t1B=t1B-rt;
#if PFETCH
	#if(CACHE_LINE_DOUBLES == 4)
	prefetch_p_doubles(add0);
	#endif
add0 += 4;
#endif
			rt =t18;		it =t19;
			t18=t10-rt;		t19=t11-it;
			t10=t10+rt;		t11=t11+it;

			rt =t1C;		it =t1D;
			t1C=t14-it;		t1D=t15+rt;
			t14=t14+it;		t15=t15-rt;

			rt =(t1A+t1B)*ISRT2;it =(t1A-t1B)*ISRT2;
			t1A=t12-rt;		t1B=t13+it;
			t12=t12+rt;		t13=t13-it;

			rt =(t1E-t1F)*ISRT2;it =(t1F+t1E)*ISRT2;
			t1E=t16+rt;		t1F=t17+it;
			t16=t16-rt;		t17=t17-it;

/*...Block 3:	*/
#if PFETCH
prefetch_p_doubles(add0);
add0 += 4;
#endif
			t22=a1p10r-a1p11r;	t23=a1p10i-a1p11i;
			t20=a1p10r+a1p11r;	t21=a1p10i+a1p11i;

			t26=a1p12r-a1p13r;	t27=a1p12i-a1p13i;
			t24=a1p12r+a1p13r;	t25=a1p12i+a1p13i;

			rt =t24;		it =t25;
			t24=t20-rt;		t25=t21-it;
			t20=t20+rt;		t21=t21+it;

			rt =t26;		it =t27;
			t26=t22-it;		t27=t23+rt;
			t22=t22+it;		t23=t23-rt;

			t2A=a1p14r-a1p15r;	t2B=a1p14i-a1p15i;
			t28=a1p14r+a1p15r;	t29=a1p14i+a1p15i;

			t2E=a1p16r-a1p17r;	t2F=a1p16i-a1p17i;
			t2C=a1p16r+a1p17r;	t2D=a1p16i+a1p17i;

			rt =t2C;		it =t2D;
			t2C=t28-rt;		t2D=t29-it;
			t28=t28+rt;		t29=t29+it;

			rt =t2E;		it =t2F;
			t2E=t2A-it;		t2F=t2B+rt;
			t2A=t2A+it;		t2B=t2B-rt;
#if PFETCH
	#if(CACHE_LINE_DOUBLES == 4)
	prefetch_p_doubles(add0);
	#endif
add0 += 4;
#endif
			rt =t28;		it =t29;
			t28=t20-rt;		t29=t21-it;
			t20=t20+rt;		t21=t21+it;

			rt =t2C;		it =t2D;
			t2C=t24-it;		t2D=t25+rt;
			t24=t24+it;		t25=t25-rt;

			rt =(t2A+t2B)*ISRT2;it =(t2A-t2B)*ISRT2;
			t2A=t22-rt;		t2B=t23+it;
			t22=t22+rt;		t23=t23-it;

			rt =(t2E-t2F)*ISRT2;it =(t2F+t2E)*ISRT2;
			t2E=t26+rt;		t2F=t27+it;
			t26=t26-rt;		t27=t27-it;

/*...Block 4:	*/
#if PFETCH
prefetch_p_doubles(add0);
add0 += 4;
#endif
			t32=a1p18r-a1p19r;	t33=a1p18i-a1p19i;
			t30=a1p18r+a1p19r;	t31=a1p18i+a1p19i;

			t36=a1p1Ar-a1p1Br;	t37=a1p1Ai-a1p1Bi;
			t34=a1p1Ar+a1p1Br;	t35=a1p1Ai+a1p1Bi;

			rt =t34;		it =t35;
			t34=t30-rt;		t35=t31-it;
			t30=t30+rt;		t31=t31+it;

			rt =t36;		it =t37;
			t36=t32-it;		t37=t33+rt;
			t32=t32+it;		t33=t33-rt;

			t3A=a1p1Cr-a1p1Dr;	t3B=a1p1Ci-a1p1Di;
			t38=a1p1Cr+a1p1Dr;	t39=a1p1Ci+a1p1Di;

			t3E=a1p1Er-a1p1Fr;	t3F=a1p1Ei-a1p1Fi;
			t3C=a1p1Er+a1p1Fr;	t3D=a1p1Ei+a1p1Fi;

			rt =t3C;		it =t3D;
			t3C=t38-rt;		t3D=t39-it;
			t38=t38+rt;		t39=t39+it;

			rt =t3E;		it =t3F;
			t3E=t3A-it;		t3F=t3B+rt;
			t3A=t3A+it;		t3B=t3B-rt;
#if PFETCH
	#if(CACHE_LINE_DOUBLES == 4)
	prefetch_p_doubles(add0);
	#endif
add0 += 4;
#endif
			rt =t38;		it =t39;
			t38=t30-rt;		t39=t31-it;
			t30=t30+rt;		t31=t31+it;

			rt =t3C;		it =t3D;
			t3C=t34-it;		t3D=t35+rt;
			t34=t34+it;		t35=t35-rt;

			rt =(t3A+t3B)*ISRT2;it =(t3A-t3B)*ISRT2;
			t3A=t32-rt;		t3B=t33+it;
			t32=t32+rt;		t33=t33-it;

			rt =(t3E-t3F)*ISRT2;it =(t3F+t3E)*ISRT2;
			t3E=t36+rt;		t3F=t37+it;
			t36=t36-rt;		t37=t37-it;

/*...and now do eight radix-4 transforms, including the internal twiddle factors:
	1, exp(-i* 1*twopi/32) =       ( c32_1,-s32_1), exp(-i* 2*twopi/32) =       ( c    ,-s    ), exp(-i* 3*twopi/32) =       ( c32_3,-s32_3) (for inputs to transform block 2),
	1, exp(-i* 2*twopi/32) =       ( c    ,-s    ), exp(-i* 4*twopi/32) = ISRT2*( 1    ,-1    ), exp(-i* 3*twopi/32) =       ( s    ,-c    ) (for inputs to transform block 3),
	1, exp(-i* 3*twopi/32) =       ( c32_3,-s32_3), exp(-i* 6*twopi/32) =       ( s    ,-c    ), exp(-i* 9*twopi/32) =       (-s32_1,-c32_1) (for inputs to transform block 4),
	1, exp(-i* 4*twopi/32) = ISRT2*( 1    ,-1    ), exp(-i* 8*twopi/32) =       ( 0    ,-1    ), exp(-i*12*twopi/32) = ISRT2*(-1    ,-1    ) (for inputs to transform block 5),
	1, exp(-i* 5*twopi/32) =       ( s32_3,-c32_3), exp(-i*10*twopi/32) =       (-s    ,-c    ), exp(-i*15*twopi/32) =       (-c32_1,-s32_1) (for inputs to transform block 6),
	1, exp(-i* 6*twopi/32) =       ( s    ,-c    ), exp(-i*12*twopi/32) = ISRT2*(-1    ,-1    ), exp(-i*18*twopi/32) =       (-c    , s    ) (for inputs to transform block 7),
	1, exp(-i* 7*twopi/32) =       ( s32_1,-c32_1), exp(-i*14*twopi/32) =       (-c    ,-s    ), exp(-i*21*twopi/32) =       (-s32_3, c32_3) (for inputs to transform block 8),
   and only the last 3 inputs to each of the radix-4 transforms 2 through 8 are multiplied by non-unity twiddles.	*/

/*...Block 1: t00,t10,t20,t30	*/
#if PFETCH
prefetch_p_doubles(add0);
add0 += 4;
#endif
			rt =t10;	t10=t00-rt;	t00=t00+rt;
			it =t11;	t11=t01-it;	t01=t01+it;

			rt =t30;	t30=t20-rt;	t20=t20+rt;
			it =t31;	t31=t21-it;	t21=t21+it;

			a[j1   ]=t00+t20;			a[j2   ]=t01+t21;
			t00       =t00-t20;				t01       =t01-t21;
			a[j1+32]=t00*c10+t01*s10;	a[j2+32]=t01*c10-t00*s10;

			rt        =t10+t31;				it        =t11-t30;
			t10       =t10-t31;				t11       =t11+t30;
			a[j1+16]=rt *c08+it *s08;	a[j2+16]=it *c08-rt *s08;
			a[j1+48]=t10*c18+t11*s18;	a[j2+48]=t11*c18-t10*s18;

/*...Block 5: t08,t18,t28,t38	*/
#if PFETCH
	#if(CACHE_LINE_DOUBLES == 4)
	prefetch_p_doubles(add0);
	#endif
add0 += 4;
#endif
			rt =t18;	t18=t08-t19;		t08=t08+t19;
					t19=t09+rt;			t09=t09-rt;

			rt =(t29+t28)*ISRT2;			t29=(t29-t28)*ISRT2;		t28=rt;
			rt =(t38-t39)*ISRT2;			it =(t38+t39)*ISRT2;
			t38=t28+rt;						t28=t28-rt;
			t39=t29+it;						t29=t29-it;

			rt        =t08+t28;				it        =t09+t29;
			t08       =t08-t28;				t09       =t09-t29;
			a[j1+8 ]=rt *c04+it *s04;	a[j2+8 ]=it *c04-rt *s04;
			a[j1+40]=t08*c14+t09*s14;	a[j2+40]=t09*c14-t08*s14;

			rt        =t18+t39;				it        =t19-t38;
			t18       =t18-t39;				t19       =t19+t38;
			a[j1+24]=rt *c0C+it *s0C;	a[j2+24]=it *c0C-rt *s0C;
			a[j1+56]=t18*c1C+t19*s1C;	a[j2+56]=t19*c1C-t18*s1C;

/*...Block 3: t04,t14,t24,t34	*/
#if PFETCH
prefetch_p_doubles(add0);
add0 += 4;
#endif
			rt =(t15+t14)*ISRT2;			it =(t15-t14)*ISRT2;
			t14=t04-rt;						t04=t04+rt;
			t15=t05-it;						t05=t05+it;

			rt =t24*c + t25*s;				t25=t25*c - t24*s;		t24=rt;
			rt =t34*s + t35*c;				it =t35*s - t34*c;
			t34=t24-rt;						t24=t24+rt;
			t35=t25-it;						t25=t25+it;

			rt        =t04+t24;				it        =t05+t25;
			t04       =t04-t24;				t05       =t05-t25;
			a[j1+4 ]=rt *c02+it *s02;	a[j2+4 ]=it *c02-rt *s02;
			a[j1+36]=t04*c12+t05*s12;	a[j2+36]=t05*c12-t04*s12;

			rt        =t14+t35;				it        =t15-t34;
			t14       =t14-t35;				t15       =t15+t34;
			a[j1+20]=rt *c0A+it *s0A;	a[j2+20]=it *c0A-rt *s0A;
			a[j1+52]=t14*c1A+t15*s1A;	a[j2+52]=t15*c1A-t14*s1A;

/*...Block 7: t0C,t1C,t2C,t3C	*/
#if PFETCH
	#if(CACHE_LINE_DOUBLES == 4)
	prefetch_p_doubles(add0);
	#endif
add0 += 4;
#endif
			rt =(t1C-t1D)*ISRT2;			it =(t1C+t1D)*ISRT2;
			t1C=t0C+rt;						t0C=t0C-rt;
			t1D=t0D+it;						t0D=t0D-it;

			rt =t2C*s + t2D*c;				t2D=t2D*s - t2C*c;		t2C=rt;
			rt =t3C*c + t3D*s;				it =t3D*c - t3C*s;
			t3C=t2C+rt;						t2C=t2C-rt;
			t3D=t2D+it;						t2D=t2D-it;

			rt        =t0C+t2C;				it        =t0D+t2D;
			t0C       =t0C-t2C;				t0D       =t0D-t2D;
			a[j1+12]=rt *c06+it *s06;	a[j2+12]=it *c06-rt *s06;
			a[j1+44]=t0C*c16+t0D*s16;	a[j2+44]=t0D*c16-t0C*s16;

			rt        =t1C+t3D;				it        =t1D-t3C;
			t1C       =t1C-t3D;				t1D       =t1D+t3C;
			a[j1+28]=rt *c0E+it *s0E;	a[j2+28]=it *c0E-rt *s0E;
			a[j1+60]=t1C*c1E+t1D*s1E;	a[j2+60]=t1D*c1E-t1C*s1E;

/*...Block 2: t02,t12,t22,t32	*/
#if PFETCH
prefetch_p_doubles(add0);
add0 += 4;
#endif
			rt =t12*c + t13*s;				it =t13*c - t12*s;
			t12=t02-rt;						t02=t02+rt;
			t13=t03-it;						t03=t03+it;

			rt =t22*c32_1 + t23*s32_1;		t23=t23*c32_1 - t22*s32_1;	t22=rt;
			rt =t32*c32_3 + t33*s32_3;		it =t33*c32_3 - t32*s32_3;
			t32=t22-rt;						t22=t22+rt;
			t33=t23-it;						t23=t23+it;

			rt        =t02+t22;				it        =t03+t23;
			t02       =t02-t22;				t03       =t03-t23;
			a[j1+2 ]=rt *c01+it *s01;	a[j2+2 ]=it *c01-rt *s01;
			a[j1+34]=t02*c11+t03*s11;	a[j2+34]=t03*c11-t02*s11;

			rt        =t12+t33;				it        =t13-t32;
			t12       =t12-t33;				t13       =t13+t32;
			a[j1+18]=rt *c09+it *s09;	a[j2+18]=it *c09-rt *s09;
			a[j1+50]=t12*c19+t13*s19;	a[j2+50]=t13*c19-t12*s19;

/*...Block 6: t0A,t1A,t2A,t3A	*/
#if PFETCH
	#if(CACHE_LINE_DOUBLES == 4)
	prefetch_p_doubles(add0);
	#endif
add0 += 4;
#endif
			rt =t1A*s - t1B*c;				it =t1B*s + t1A*c;
			t1A=t0A+rt;						t0A =t0A-rt;
			t1B=t0B+it;						t0B =t0B-it;

			rt =t2A*s32_3 + t2B*c32_3;		t2B=t2B*s32_3 - t2A*c32_3;	t2A=rt;
			rt =t3A*c32_1 - t3B*s32_1;		it =t3B*c32_1 + t3A*s32_1;
			t3A=t2A+rt;						t2A=t2A-rt;
			t3B=t2B+it;						t2B=t2B-it;

			rt        =t0A+t2A;				it        =t0B+t2B;
			t0A       =t0A-t2A;				t0B       =t0B-t2B;
			a[j1+10]=rt *c05+it *s05;	a[j2+10]=it *c05-rt *s05;
			a[j1+42]=t0A*c15+t0B*s15;	a[j2+42]=t0B*c15-t0A*s15;

			rt        =t1A+t3B;				it        =t1B-t3A;
			t1A       =t1A-t3B;				t1B       =t1B+t3A;
			a[j1+26]=rt *c0D+it *s0D;	a[j2+26]=it *c0D-rt *s0D;
			a[j1+58]=t1A*c1D+t1B*s1D;	a[j2+58]=t1B*c1D-t1A*s1D;

/*...Block 4: t06,t16,t26,t36	*/
#if PFETCH
prefetch_p_doubles(add0);
add0 += 4;
#endif
			rt =t16*s + t17*c;				it =t17*s - t16*c;
			t16=t06-rt;						t06 =t06+rt;
			t17=t07-it;						t07 =t07+it;

			rt =t26*c32_3 + t27*s32_3;		t27=t27*c32_3 - t26*s32_3;	t26=rt;
			rt =t36*s32_1 - t37*c32_1;		it =t37*s32_1 + t36*c32_1;
			t36=t26+rt;						t26=t26-rt;
			t37=t27+it;						t27=t27-it;

			rt        =t06+t26;				it        =t07+t27;
			t06       =t06-t26;				t07       =t07-t27;
			a[j1+6 ]=rt *c03+it *s03;	a[j2+6 ]=it *c03-rt *s03;
			a[j1+38]=t06*c13+t07*s13;	a[j2+38]=t07*c13-t06*s13;

			rt        =t16+t37;				it        =t17-t36;
			t16       =t16-t37;				t17       =t17+t36;
			a[j1+22]=rt *c0B+it *s0B;	a[j2+22]=it *c0B-rt *s0B;
			a[j1+54]=t16*c1B+t17*s1B;	a[j2+54]=t17*c1B-t16*s1B;

/*...Block 8: t0E,t1E,t2E,t3E	*/
#if PFETCH
	#if(CACHE_LINE_DOUBLES == 4)
	prefetch_p_doubles(add0);
	#endif
add0 += 4;
#endif
			rt =t1E*c - t1F*s;				it =t1F*c + t1E*s;
			t1E=t0E+rt;						t0E =t0E-rt;
			t1F=t0F+it;						t0F =t0F-it;

			rt =t2E*s32_1 + t2F*c32_1;		t2F=t2F*s32_1 - t2E*c32_1;	t2E=rt;
			rt =t3E*s32_3 + t3F*c32_3;		it =t3F*s32_3 - t3E*c32_3;
			t3E=t2E+rt;						t2E=t2E-rt;
			t3F=t2F+it;						t2F=t2F-it;

			rt        =t0E+t2E;				it        =t0F+t2F;
			t0E       =t0E-t2E;				t0F       =t0F-t2F;
			a[j1+14]=rt *c07+it *s07;	a[j2+14]=it *c07-rt *s07;
			a[j1+46]=t0E*c17+t0F*s17;	a[j2+46]=t0F*c17-t0E*s17;

			rt        =t1E+t3F;				it        =t1F-t3E;
			t1E       =t1E-t3F;				t1F       =t1F+t3E;
			a[j1+30]=rt *c0F+it *s0F;	a[j2+30]=it *c0F-rt *s0F;
			a[j1+62]=t1E*c1F+t1F*s1F;	a[j2+62]=t1F*c1F-t1E*s1F;

	}
}

