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

	#undef DEBUG_SSE2
//	#define DEBUG_SSE2

#if defined(USE_SSE2)

	#include "sse2_macro.h"

	#if defined(COMPILER_TYPE_MSVC)

		/* DIT radix-4 subconvolution, sans twiddles - this is the in-place version needed by the wrapper_square routines.
		Assumes the base address __add0 enters in eax:
		*/
		#define SSE2_RADIX4_DIT_IN_PLACE_B()\
		{\
			__asm	movaps	xmm0,[eax+0x100]	/* a[jt+p2] */\
			__asm	movaps	xmm1,[eax+0x110]	/* a[jp+p2] */\
			__asm	movaps	xmm2,[eax      ]	/* a[jt   ] */\
			__asm	movaps	xmm3,[eax+0x010]	/* a[jp   ] */\
			\
			__asm	subpd	xmm2,xmm0	/* t3 */				__asm	movaps	xmm4,[eax+0x180]	/* a[jt+p3] */\
			__asm	subpd	xmm3,xmm1	/* t4 */				__asm	movaps	xmm5,[eax+0x190]	/* a[jp+p3] */\
			__asm	addpd	xmm0,xmm0	/* 2*y */				__asm	movaps	xmm6,[eax+0x080]	/* a[jt+p1] */\
			__asm	addpd	xmm1,xmm1	/* 2*y */				__asm	movaps	xmm7,[eax+0x090]	/* a[jp+p1] */\
			__asm	addpd	xmm0,xmm2	/* t1 */				\
			__asm	addpd	xmm1,xmm3	/* t2 */				__asm	subpd	xmm6,xmm4	/* t7 */\
																__asm	subpd	xmm7,xmm5	/* t8 */\
			__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */		__asm	addpd	xmm4,xmm4	/* 2*y */\
			__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */		__asm	addpd	xmm5,xmm5	/* 2*y */\
			__asm	movaps	[eax+0x090],xmm3	/* <- ~t4 */	__asm	addpd	xmm4,xmm6	/* t4 */\
			__asm	movaps	[eax+0x180],xmm2	/* <- ~t7 */	__asm	addpd	xmm5,xmm7	/* t5 */\
			\
			__asm	addpd	xmm7,xmm7	/*          2*t8 */		__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */\
			__asm	addpd	xmm6,xmm6	/*          2*t7 */		__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */\
			__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */		__asm	movaps	[eax+0x100],xmm0	/* <- ~t5 */\
			__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */		__asm	movaps	[eax+0x110],xmm1	/* <- ~t6 */\
			__asm	movaps	[eax+0x080],xmm7	/* <- ~t3 */	__asm	addpd	xmm4,xmm4	/*          2*t5 */\
			__asm	movaps	[eax+0x190],xmm6	/* <- ~t8 */	__asm	addpd	xmm5,xmm5	/*          2*t6 */\
																__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */\
																__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */\
																__asm	movaps	[eax      ],xmm4	/* <- ~t1 */\
																__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */\
		}
		/* Assumes a[jt,jp]+[p0,p2,p1,p3] enter in __r0,__r1,__r2,__r3,__r4,__r5,__r6,__r7, and *two in esi: */
		#define SSE2_RADIX4_DIT_IN_PLACE_C(__r0,__r1,__r2,__r3,__r4,__r5,__r6,__r7)\
		{\
			__asm	subpd	__r0,__r2	/* t3 */				\
			__asm	subpd	__r1,__r3	/* t4 */				\
			__asm	addpd	__r2,__r2	/* 2*y */				\
			__asm	addpd	__r3,__r3	/* 2*y */				\
			__asm	addpd	__r2,__r0	/* t1 */				\
			__asm	addpd	__r3,__r1	/* t2 */				__asm	subpd	__r4,__r6	/* t7 */\
																__asm	subpd	__r5,__r7	/* t8 */\
			__asm	subpd	__r1,__r4	/* ~t4 <- t4 -t7 */		__asm	mulpd	__r6,[esi]	/* 2*y */\
			__asm	subpd	__r0,__r5	/* ~t7 <- t3 -t8 */		__asm	mulpd	__r7,[esi]	/* 2*y */\
			__asm	movaps	[eax+0x090],__r1	/* <- ~t4 */	__asm	addpd	__r6,__r4	/* t4 */\
			__asm	movaps	[eax+0x180],__r0	/* <- ~t7 */	__asm	addpd	__r7,__r5	/* t5 */\
			\
			__asm	addpd	__r5,__r5	/*          2*t8 */		__asm	subpd	__r2,__r6	/* ~t5 <- t1 -t5 */\
			__asm	addpd	__r4,__r4	/*          2*t7 */		__asm	subpd	__r3,__r7	/* ~t6 <- t2 -t6 */\
			__asm	addpd	__r5,__r0	/* ~t3 <- t3 +t8 */		__asm	movaps	[eax+0x100],__r2	/* <- ~t5 */\
			__asm	addpd	__r4,__r1	/* ~t8 <- t4 +t7 */		__asm	movaps	[eax+0x110],__r3	/* <- ~t6 */\
			__asm	movaps	[eax+0x080],__r5	/* <- ~t3 */	__asm	addpd	__r6,__r6	/*          2*t5 */\
			__asm	movaps	[eax+0x190],__r4	/* <- ~t8 */	__asm	addpd	__r7,__r7	/*          2*t6 */\
																__asm	addpd	__r6,__r2	/* ~t1 <- t1 +t5 */\
																__asm	addpd	__r7,__r3	/* ~t2 <- t2 +t6 */\
																__asm	movaps	[eax      ],__r6	/* <- ~t1 */\
																__asm	movaps	[eax+0x010],__r7	/* <- ~t2 */\
		}

	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "radix16_wrapper_square_gcc32.h"

		#else

			#include "radix16_wrapper_square_gcc64.h"

		#endif

	#endif

#endif

/***************/

void radix16_dyadic_square(double a[], int arr_scratch[], int n, int radix0, struct complex rt0[], struct complex rt1[], int ii, int nradices_prim, int radix_prim[], int nloops, int incr)
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

The scratch array (2nd input argument) is only needed for data table initializations, i.e. if first_entry = TRUE.
*/
	static int nsave = 0;
	static int rad0save = 0;
/*	static int *index = 0x0;	OBSOLETE: full-N2/16-length Bit-reversal index array. */
	static int *index0 = 0x0, *index1 = 0x0, *index_ptmp0 = 0x0, *index_ptmp1 = 0x0;
	static int index0_idx, index0_mod, index1_idx, index1_mod;
	int nradices_prim_radix0;

	int i,j,j1,j2,k,l,iroot,k1,k2,ndivrad0;
	const double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;	/* exp[i*(twopi/16)] */
	double re0,im0,re1,im1;
	double rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32
	,aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r,aj1p6r,aj1p7r,aj1p8r,aj1p9r,aj1p10r,aj1p11r,aj1p12r,aj1p13r,aj1p14r,aj1p15r
	,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i,aj1p6i,aj1p7i,aj1p8i,aj1p9i,aj1p10i,aj1p11i,aj1p12i,aj1p13i,aj1p14i,aj1p15i;

#ifdef USE_SSE2

	double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
	static struct complex *sc_arr = 0x0, *sc_ptr;
	static struct complex *cc0, *ss0, *isrt2, *two;
	static struct complex *c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7,*c8,*c9,*c10,*c11,*c12,*c13,*c14,*c15,*s0,*s1,*s2,*s3,*s4,*s5,*s6,*s7,*s8,*s9,*s10,*s11,*s12,*s13,*s14,*s15;
	static struct complex *r1,*r2,*r3,*r4,*r5,*r6,*r7,*r8,*r9,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r30,*r31,*r32;

#else

  #if PFETCH
	double *add0;
  #endif
	double c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15;

#endif

	static int first_entry=TRUE;

#ifdef DEBUG_SSE2
	int iloop;
#endif

/*...initialize things upon first entry */
/*...If a new runlength or first-pass radix, set first_entry to true:	*/

	if(n != nsave || radix0 != rad0save) first_entry=TRUE;

	if(first_entry)
	{
		first_entry=FALSE;
		nsave = n;
		rad0save = radix0;

#ifdef USE_SSE2
		sc_arr = ALLOC_COMPLEX(sc_arr, 72);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 32 16-byte slots of sc_arr for temporaries, next 3 for the nontrivial complex 16th roots,
	last 30 for the doubled sincos twiddles, plus at least 3 more slots to allow for 64-byte alignment of the array.
	*/
		r1  = sc_ptr + 0x00;	  isrt2 = sc_ptr + 0x20;
		r2  = sc_ptr + 0x01;		cc0 = sc_ptr + 0x21;
		r3  = sc_ptr + 0x02;		ss0 = sc_ptr + 0x22;
		r4  = sc_ptr + 0x03;		c0  = sc_ptr + 0x23;
		r5  = sc_ptr + 0x04;		s0  = sc_ptr + 0x24;
		r6  = sc_ptr + 0x05;		c8  = sc_ptr + 0x25;
		r7  = sc_ptr + 0x06;		s8  = sc_ptr + 0x26;
		r8  = sc_ptr + 0x07;		c4  = sc_ptr + 0x27;
		r9  = sc_ptr + 0x08;		s4  = sc_ptr + 0x28;
		r10 = sc_ptr + 0x09;		c12 = sc_ptr + 0x29;
		r11 = sc_ptr + 0x0a;		s12 = sc_ptr + 0x2a;
		r12 = sc_ptr + 0x0b;		c2  = sc_ptr + 0x2b;
		r13 = sc_ptr + 0x0c;		s2  = sc_ptr + 0x2c;
		r14 = sc_ptr + 0x0d;		c10 = sc_ptr + 0x2d;
		r15 = sc_ptr + 0x0e;		s10 = sc_ptr + 0x2e;
		r16 = sc_ptr + 0x0f;		c6  = sc_ptr + 0x2f;
		r17 = sc_ptr + 0x10;		s6  = sc_ptr + 0x30;
		r18 = sc_ptr + 0x11;		c14 = sc_ptr + 0x31;
		r19 = sc_ptr + 0x12;		s14 = sc_ptr + 0x32;
		r20 = sc_ptr + 0x13;		c1  = sc_ptr + 0x33;
		r21 = sc_ptr + 0x14;		s1  = sc_ptr + 0x34;
		r22 = sc_ptr + 0x15;		c9  = sc_ptr + 0x35;
		r23 = sc_ptr + 0x16;		s9  = sc_ptr + 0x36;
		r24 = sc_ptr + 0x17;		c5  = sc_ptr + 0x37;
		r25 = sc_ptr + 0x18;		s5  = sc_ptr + 0x38;
		r26 = sc_ptr + 0x19;		c13 = sc_ptr + 0x39;
		r27 = sc_ptr + 0x1a;		s13 = sc_ptr + 0x3a;
		r28 = sc_ptr + 0x1b;		c3  = sc_ptr + 0x3b;
		r29 = sc_ptr + 0x1c;		s3  = sc_ptr + 0x3c;
		r30 = sc_ptr + 0x1d;		c11 = sc_ptr + 0x3d;
		r31 = sc_ptr + 0x1e;		s11 = sc_ptr + 0x3e;
		r32 = sc_ptr + 0x1f;		c7  = sc_ptr + 0x3f;
									s7  = sc_ptr + 0x40;
									c15 = sc_ptr + 0x41;
									s15 = sc_ptr + 0x42;
									two = sc_ptr + 0x43;

		/* These remain fixed: */
		isrt2->re = ISRT2;	isrt2->im = ISRT2;
		two  ->re = 2.0;	two  ->im = 2.0;
		cc0  ->re = c	;	cc0  ->im = c	;
		ss0  ->re = s	;	ss0  ->im = s	;

		c0->re=1.;	s0->re=0.;
		c0->im=1.;	s0->im=0.;
#endif

		free((void *)index_ptmp0);	index_ptmp0=0x0;	index0=0x0;
		free((void *)index_ptmp1);	index_ptmp1=0x0;	index1=0x0;

		/*
		!...Final forward (DIF) FFT pass sincos data start out in bit-reversed order.
		!   Allocate and initialize an index array containing N/16 indices...

		index_ptmp = ALLOC_INT(N2/16);
		index = ALIGN_INT(index_ptmp);
		if(!index){ sprintf(cbuf,"FATAL: unable to allocate array ITMP in radix16_dyadic_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		for(i=0; i < N2/16; i++)
		{
			index[i]=i;
		}
		*/
		index0_mod =        radix0;
		index1_mod = (n>>5)/radix0;	/* complex length requires an additional divide by 2 */

		index_ptmp0 = ALLOC_INT(index_ptmp0, index0_mod);
		if(!index_ptmp0){ sprintf(cbuf,"FATAL: unable to allocate array INDEX_PTMP0 in radix16_dyadic_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		index0 = ALIGN_INT(index_ptmp0);

		index_ptmp1 = ALLOC_INT(index_ptmp1, index1_mod);
		if(!index_ptmp1){ sprintf(cbuf,"FATAL: unable to allocate array INDEX_PTMP1 in radix16_dyadic_square.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		index1 = ALIGN_INT(index_ptmp1);

		for(i=0; i < index0_mod; i++){index0[i]=       i;}
		for(i=0; i < index1_mod; i++){index1[i]=radix0*i;}

		/*...then bit-reverse INDEX with respect to N/16. Radices are sent to BR routine
		in opposite order from which they are processed in forward FFT.

		bit_reverse_int(index, N2/16, nradices_prim-4, &radix_prim[nradices_prim-5], -1,(int *)arr_scratch);
		*/

		i = 1;
		for(nradices_prim_radix0 = 1; nradices_prim_radix0 < nradices_prim; nradices_prim_radix0++)
		{
			i *= radix_prim[nradices_prim_radix0-1];
			if(i == radix0)
				break;
		}
		ASSERT(HERE, nradices_prim_radix0 < nradices_prim,"radix16_dyadic_square.c: nradices_prim_radix0 < nradices_prim");

		bit_reverse_int(index0, index0_mod,                 nradices_prim_radix0, &radix_prim[nradices_prim_radix0-1], -1,(int *)arr_scratch);
		bit_reverse_int(index1, index1_mod, nradices_prim-4-nradices_prim_radix0, &radix_prim[nradices_prim       -5], -1,(int *)arr_scratch);

		/*****DEBUG: ******/
		/*
		fp = fopen("index.txt","a");
		fprintf(fp,"FFT length = %d\n",n);
		fprintf(fp,"Using complex FFT radices"); for(i = 0; i < NRADICES; i++){ fprintf(fp,"%10d",RADIX_VEC[i]); }; fprintf(fp,"\n");
		fprintf(fp,"radix16_dyadic_square: index array elements are:\n");
		for(i=0; i < N2/16; i++)
		{
			fprintf(fp,"%10d",index[i]);
		}
		fprintf(fp,"\n");
		fprintf(fp,"index0 array:\n");
		for(i=0; i <       radix0; i++){fprintf(fp,"%10d",index0[i]);}fprintf(fp,"\n");
		for(i=0; i < N2/16/radix0; i++){fprintf(fp,"%10d",index1[i]);}fprintf(fp,"\n");
		fclose(fp);	fp = 0x0;
		*/
		/******************/

	}

	/*...If a new runlength, should not get to this point: */
	ASSERT(HERE, n == nsave,"radix16_dyadic_square: n != nsave");
	ASSERT(HERE, incr == 32,"radix16_dyadic_square.c: incr != 32");
	ndivrad0 = n/radix0;
	/*
	k = ii*(ndivrad0 >> 5);
	*/
	index0_idx = ii;
	index1_idx = 0;

	#ifdef USE_SSE2
	for(j = 0; j < ndivrad0; j += 64)
	{
		j1 = (j & mask01) + br4[j&3];
	    j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
	for(j = 0; j < ndivrad0; j += 32)
	{
	    j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;
	/*
		iroot = index[k];
		ASSERT(HERE, iroot == index0[index0_idx] + index1[index1_idx],"radix16_dyadic_square.c: iroot == index0[index0_idx] + index1[index1_idx]");
		k = k + 1;	// increment sincos array index
	*/

		iroot = index0[index0_idx] + index1[index1_idx];

		if(++index1_idx >= index1_mod)
		{
			index1_idx -= index1_mod;
		/*	ASSERT(HERE, index0_idx < index0_mod,"radix16_dyadic_square.c: index0_idx < index0_mod");	*/
			++index0_idx;
		}

		l = iroot;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;		 /* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c1 ->re=rt;	s1 ->re=it;
	#else
		c1 =rt;	s1 =it;
	#endif

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 1);	/* 4*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c2 ->re=rt;	s2 ->re=it;
	#else
		c2 =rt;	s2 =it;
	#endif

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 2);	/* 8*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c4 ->re=rt;	s4 ->re=it;
	#else
		c4 =rt;	s4 =it;
	#endif

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 2) + iroot;	/* 13*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c8 ->re=rt;	s8 ->re=it;
	#else
		c8 =rt;	s8 =it;
	#endif

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c13->re=rt;	s13->re=it;
	#else
		c13 =rt;	s13 =it;
	#endif

	/* In SSE2 mode, also need next set of sincos to put into "Imaginary" slots: */
	#ifdef USE_SSE2
		iroot = index0[index0_idx] + index1[index1_idx];

		if(++index1_idx >= index1_mod)
		{
			index1_idx -= index1_mod;
		/*	ASSERT(HERE, index0_idx < index0_mod,"radix16_dyadic_square.c: index0_idx < index0_mod");	*/
			++index0_idx;
		}

		l = iroot;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += iroot;		 /* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c1 ->im=rt;	s1 ->im=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 1);	/* 4*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c2 ->im=rt;	s2 ->im=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 2);	/* 8*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c4 ->im=rt;	s4 ->im=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		l += (iroot << 2) + iroot;	/* 13*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c8 ->im=rt;	s8 ->im=it;

		k1=(l & NRTM1);
		k2=(l >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
		c13->im=rt;	s13->im=it;

		SSE2_CMUL_EXPO(c1,c4 ,c3 ,c5 )
		SSE2_CMUL_EXPO(c1,c8 ,c7 ,c9 )
		SSE2_CMUL_EXPO(c2,c8 ,c6 ,c10)
		SSE2_CMUL_EXPO(c1,c13,c12,c14)
		SSE2_CMUL_EXPO(c2,c13,c11,c15)
	#else
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
	#endif

	#ifdef DEBUG_SSE2
		rng_isaac_init(TRUE);
		for(iloop = 0; iloop < 64; iloop += 4)
		{
	  #ifdef USE_SSE2
			a[j1+iloop  ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[j1+iloop+2] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[j1+iloop+1] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[j1+iloop+3] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			fprintf(stderr, "radix16_dyadic_square: A_in[%2d] = %20.5f, %20.5f\n",iloop,a[j1+iloop+0],a[j1+iloop+2]);
			fprintf(stderr, "radix16_dyadic_square: A_in[%2d] = %20.5f, %20.5f\n",iloop,a[j1+iloop+1],a[j1+iloop+3]);
	  #else
			a[j1+iloop  ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[j1+iloop+1] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[j1+iloop+2] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[j1+iloop+3] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			fprintf(stderr, "radix16_dyadic_square: A_in[%2d] = %20.5f, %20.5f\n",iloop,a[j1+iloop+0],a[j1+iloop+1]);
			fprintf(stderr, "radix16_dyadic_square: A_in[%2d] = %20.5f, %20.5f\n",iloop,a[j1+iloop+2],a[j1+iloop+3]);
	  #endif
		}
	#endif

/*   gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do the first set of four length-4 transforms... */

#ifdef USE_SSE2

	/* In sse2 mode, the input data are arranged in memory like so, where we view things in 16-byte chunks:

		&a[j1]	a0.re,a1.re  &a[j1+32]	b0.re,b1.re
		+0x010	a0.im,a1.im		+0x010	b0.im,b1.im
		+0x020	a2.re,a3.re		+0x020	b2.re,b3.re
		+0x030	a2.im,a3.im		+0x030	b2.im,b3.im
		+0x040	a4.re,a5.re		+0x040	b4.re,b5.re
		+0x050	a4.im,a5.im		+0x050	b4.im,b5.im
		+0x060	a6.re,a7.re		+0x060	b6.re,b7.re
		+0x070	a6.im,a7.im		+0x070	b6.im,b7.im
		+0x080	a8.re,a9.re		+0x080	b8.re,b9.re
		+0x090	a8.im,a9.im		+0x090	b8.im,b9.im
		+0x0a0	aA.re,aB.re		+0x0a0	bA.re,bB.re
		+0x0b0	aA.im,aB.im		+0x0b0	bA.im,bB.im
		+0x0c0	aC.re,aD.re		+0x0c0	bC.re,bD.re
		+0x0d0	aC.im,aD.im		+0x0d0	bC.im,bD.im
		+0x0e0	aE.re,aF.re		+0x0e0	bE.re,bF.re
		+0x0f0	aE.im,aF.im		+0x0f0	bE.im,bF.im

	We need to interleave these pairwise so as to swap the high word of each A-pair with the low word of the corresponding B-pair, e.g

				low		high	low		high
				[a0.re,a1.re]	[b0.re,b1.re]
				   |      \       /      |
				   |        \   /        |
				   |          x          |
				   |        /   \        |
				   V      /       \      V
				[a0.re,b0.re]	[a1.re,b1.re]
			=	[ A.lo, B.lo]	[ A.hi, B.hi]
					(1)				(2)

	The instructions needed for this permutation [assuming A-term in memA, B-term in memB] is

		(0)		movaps		xmm0,memA
		(1)		movaps		xmm1,memA
		(1)		unpcklpd	xmm0,xmmB		[xmm0.lo,xmm0.hi] <- [A.lo,B.lo]
		(2)		unpckhpd	xmm1,xmmB		[xmm1.lo,xmm1.hi] <- [A.hi,B.hi] ,

	Alternatively, we can use the SHUFPD instruction.
	SHUFPD xmm,mem,imm produces [xmm.lo,xmm.hi] <- [xmm.(lo or hi),xmm.(lo or hi), depending on the value of imm:

		IMM:
		0		shufpd	xmm,mem		[xmm.lo,xmm.hi] <- [xmm.lo,mem.lo]
		1		shufpd	xmm,mem		[xmm.lo,xmm.hi] <- [xmm.hi,mem.lo]
		2		shufpd	xmm,mem		[xmm.lo,xmm.hi] <- [xmm.lo,mem.hi]
		3		shufpd	xmm,mem		[xmm.lo,xmm.hi] <- [xmm.hi,mem.hi]

	So what we need is

		(0)		movaps		xmm0,memA
		(1)		movaps		xmm1,memA
		(2)		shufpd		xmm0,memB,0		[xmm0.lo,xmm0.hi] <- [A.lo,B.lo]
		(3)		shufpd		xmm1,memB,3		[xmm1.lo,xmm1.hi] <- [A.hi,B.hi]

	It's not clear whether there is a preference for one or the other instruction sequence based on resulting performance.
	*/

		add0 = &a[j1];
		add1 = &a[j1+32];

	#ifdef COMPILER_TYPE_MSVC

	/*...Block 1: */
		__asm	mov	eax, add0
		__asm	mov	ebx, add1
		__asm	mov	ecx, r1

		__asm	mov	edx, c4

		/* For interleaved [j1,j2] version, replace e.g.

			__asm	movaps	xmm0,[eax+0x40]	// a[jt+p4]
			__asm	movaps	xmm1,[eax+0x50]	// a[jp+p4]
			__asm	movaps	xmm2,[eax+0x40]	// xmm2 <- cpy a[jt+p4]
			__asm	movaps	xmm3,[eax+0x50]	// xmm3 <- cpy a[jp+p4]

		by the following:
		*/
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x40]	/* a[j1+p4], this is the scratch xmm register  */
		__asm	movaps		xmm0,xmm6		/* a[j1+p4] copy, his is the active  xmm register */
		__asm	movaps		xmm2,[ebx+0x40]	/* a[j2+p4] */
		__asm	movaps		xmm3,xmm2		/* a[j2+p4] copy */
		__asm	unpckhpd	xmm6,xmm2
		__asm	unpcklpd	xmm0,xmm3
		__asm	movaps	[ecx+0x140],xmm6	/* Store hi real in t21 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x50]
		__asm	movaps		xmm1,xmm7
		__asm	movaps		xmm4,[ebx+0x50]
		__asm	movaps		xmm5,xmm4
		__asm	unpckhpd	xmm7,xmm4
		__asm	unpcklpd	xmm1,xmm5
		__asm	movaps	[ecx+0x150],xmm7	/* Store hi imag in t22 */

		__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy a[jt+p4] */
		__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy a[jp+p4] */
		/*************************************************************************************/
		/******** From here on, things are identical to the code in radix16_dif_pass: ********/
		/*************************************************************************************/
		__asm	mulpd	xmm0,[edx     ]	/* a[jt+p4]*c4 */
		__asm	mulpd	xmm1,[edx     ]	/* a[jp+p4]*c4 */
		__asm	mulpd	xmm2,[edx+0x10]	/* a[jt+p4]*s4 */
		__asm	mulpd	xmm3,[edx+0x10]	/* a[jp+p4]*s4 */
		__asm	addpd	xmm1,xmm2	/* xmm1 <- t6 */
		__asm	subpd	xmm0,xmm3	/* xmm0 <- t5 */
		__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy t6 */
		__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy t5 */

		__asm	mov	edx, c12
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0xc0]	/* a[j1+p12], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0xc0]	/* a[j1+p12], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0xc0]	/* a[j2+p12] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0xc0]	/* a[jt+p12] */
		__asm	movaps	[ecx+0x160],xmm6	/* Store hi real in t23 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0xd0]
		__asm	movaps		xmm5,[eax+0xd0]
		__asm	unpckhpd	xmm7,[ebx+0xd0]
		__asm	unpcklpd	xmm5,[ebx+0xd0]	/* a[jp+p12] */
		__asm	movaps	[ecx+0x170],xmm7	/* Store hi imag in t24 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p12] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p12] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p12]*c12 */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p12]*c12 */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p12]*s12 */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p12]*s12 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt		xmm6,7 free */

		__asm	addpd	xmm0,xmm4	/* ~t5 <- t5 +rt */
		__asm	addpd	xmm1,xmm5	/* ~t6 <- t6 +it */
		__asm	subpd	xmm2,xmm4	/* ~t7 <- t5 -rt */
		__asm	subpd	xmm3,xmm5	/* ~t8 <- t6 -it	xmm4,5 free */

		/* Now do the p0,8 combo: */
		__asm	mov	edx, c8
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x80]	/* a[j1+p8 ], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0x80]	/* a[j1+p8 ], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x80]	/* a[j2+p8 ] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0x80]	/* a[jt+p8 ] */
		__asm	movaps	[ecx+0x120],xmm6	/* Store hi real in t19 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x90]
		__asm	movaps		xmm5,[eax+0x90]
		__asm	unpckhpd	xmm7,[ebx+0x90]
		__asm	unpcklpd	xmm5,[ebx+0x90]	/* a[jp+p8 ] */
		__asm	movaps	[ecx+0x130],xmm7	/* Store hi imag in t20 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p8] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p8] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p8]*c8 */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p8]*c8 */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p8]*s8 */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p8]*s8 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt 	xmm6,7 free - stick t1,t2 in those */

		/* Real parts: */
		__asm	movaps		xmm6,[eax     ]	/* a[j1    ], this is the scratch xmm register  */
		__asm	movaps		xmm7,[eax     ]	/* a[j1    ], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx     ]	/* a[j2    ] gets read twice */
		__asm	unpcklpd	xmm7,[ebx     ]	/* a[jt] = t1*/
		__asm	movaps	[ecx+0x100],xmm6	/* Store hi real in t17 */
		__asm	movaps	[ecx      ],xmm7	/* Store active  in t1  */
		/* Imag parts: */
		__asm	movaps		xmm6,[eax+0x10]
		__asm	movaps		xmm7,[eax+0x10]
		__asm	unpckhpd	xmm6,[ebx+0x10]
		__asm	unpcklpd	xmm7,[ebx+0x10]	/* a[jp] = t2*/
		__asm	movaps	[ecx+0x110],xmm6	/* Store hi imag in t18... */
		__asm	movaps	xmm6,[ecx      ]	/* ...and reload t1. */

		__asm	subpd	xmm6,xmm4	/* ~t3 <- t1 -rt */
		__asm	subpd	xmm7,xmm5	/* ~t4 <- t2 -it */
		__asm	addpd	xmm4,xmm4	/*          2*rt */
		__asm	addpd	xmm5,xmm5	/*          2*it */
		__asm	addpd	xmm4,xmm6	/* ~t1 <- t1 +rt */
		__asm	addpd	xmm5,xmm7	/* ~t2 <- t2 +it	xmm4,5 free */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */
		__asm	mov	eax, r1
		/*
		~t5 =t1 -t5;		~t1 =t1 +t5;
		~t6 =t2 -t6;		~t2 =t2 +t6;
		*/
		__asm	subpd	xmm4,xmm0	/*~t5 =t1 -t5 */
		__asm	subpd	xmm5,xmm1	/*~t6 =t2 -t6 */
		__asm	movaps	[eax+0x040],xmm4	/* a[jt+p8 ] <- ~t5 */
		__asm	movaps	[eax+0x050],xmm5	/* a[jp+p8 ] <- ~t6 */
		__asm	addpd	xmm0,xmm0	/* 2*t5 */
		__asm	addpd	xmm1,xmm1	/* 2*t6 */
		__asm	addpd	xmm0,xmm4	/*~t1 =t1 +t5 */
		__asm	addpd	xmm1,xmm5	/*~t2 =t2 +t6 */
		__asm	movaps	[eax      ],xmm0	/* a[jt    ] <- ~t1 */
		__asm	movaps	[eax+0x010],xmm1	/* a[jp    ] <- ~t2 */

		/*
		~t7 =t3 +t8;		~t3 =t3 -t8;
		~t8 =t4 -t7;		~t4 =t4 +t7;
		*/
		__asm	subpd	xmm6,xmm3	/*~t3 =t3 -t8 */
		__asm	subpd	xmm7,xmm2	/*~t8 =t4 -t7 */
		__asm	movaps	[eax+0x020],xmm6	/* a[jt+p4 ] <- ~t3 */
		__asm	movaps	[eax+0x070],xmm7	/* a[jp+p12] <- ~t8 */
		__asm	addpd	xmm3,xmm3	/* 2*t8 */
		__asm	addpd	xmm2,xmm2	/* 2*t7 */
		__asm	addpd	xmm3,xmm6	/*~t7 =t3 +t8 */
		__asm	addpd	xmm2,xmm7	/*~t4 =t4 +t7 */
		__asm	movaps	[eax+0x060],xmm3	/* a[jt+p12] <- ~t7 */
		__asm	movaps	[eax+0x030],xmm2	/* a[jp+p4 ] <- ~t4 */

	/*...Block 2:		Cost: 46 MOVapd, 16 UNPCKHPD, 28 ADD/SUBpd, 16 MULpd */
		__asm	mov	eax, add0
		__asm	mov	ebx, add1
		__asm	mov	ecx, r9

	/* Do the p2,10 combo: */
		__asm	mov	edx, c2
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x20]	/* a[j1+p2 ], this is the scratch xmm register */
		__asm	movaps		xmm0,[eax+0x20]	/* a[j1+p2 ], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x20]	/* a[j2+p2 ] gets read twice */
		__asm	unpcklpd	xmm0,[ebx+0x20]	/* a[jt+p2 ] */
		__asm	movaps	[ecx+0x100],xmm6	/* Store hi real in t9 +16 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x30]
		__asm	movaps		xmm1,[eax+0x30]
		__asm	unpckhpd	xmm7,[ebx+0x30]
		__asm	unpcklpd	xmm1,[ebx+0x30]	/* a[jp+p2 ] */
		__asm	movaps	[ecx+0x110],xmm7	/* Store hi imag in t10+16 */

		__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy a[jt+p2] */
		__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy a[jp+p2] */

		__asm	mulpd	xmm0,[edx     ]	/* a[jt+p2]*c2 */
		__asm	mulpd	xmm1,[edx     ]	/* a[jp+p2]*c2 */
		__asm	mulpd	xmm2,[edx+0x10]	/* a[jt+p2]*s2 */
		__asm	mulpd	xmm3,[edx+0x10]	/* a[jp+p2]*s2 */
		__asm	addpd	xmm1,xmm2	/* xmm1 <- t10*/
		__asm	subpd	xmm0,xmm3	/* xmm0 <- t9 */
		__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy t10*/
		__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy t9 */

		__asm	mov	edx, c10
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0xa0]	/* a[j1+p10], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0xa0]	/* a[j1+p10], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0xa0]	/* a[j2+p10] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0xa0]	/* a[jt+p10] */
		__asm	movaps	[ecx+0x120],xmm6	/* Store hi real in t11+16 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0xb0]
		__asm	movaps		xmm5,[eax+0xb0]
		__asm	unpckhpd	xmm7,[ebx+0xb0]
		__asm	unpcklpd	xmm5,[ebx+0xb0]	/* a[jp+p10] */
		__asm	movaps	[ecx+0x130],xmm7	/* Store hi imag in t12+16 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p10] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p10] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p10]*c10 */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p10]*c10 */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p10]*s10 */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p10]*s10 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt		xmm6,7 free */

		__asm	addpd	xmm0,xmm4	/* ~t13<- t13+rt */
		__asm	addpd	xmm1,xmm5	/* ~t14<- t14+it */
		__asm	subpd	xmm2,xmm4	/* ~t15<- t13-rt */
		__asm	subpd	xmm3,xmm5	/* ~t16<- t14-it	xmm4,5 free */

	/* Do the p6,14 combo - do p14 first so register assignments come out in same relative order as for p2,10 */
		__asm	mov	edx, c14
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0xe0]	/* a[j1+p14], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0xe0]	/* a[j1+p14], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0xe0]	/* a[j2+p14] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0xe0]	/* a[jt+p14] */
		__asm	movaps	[ecx+0x160],xmm6	/* Store hi real in t15+16 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0xf0]
		__asm	movaps		xmm5,[eax+0xf0]
		__asm	unpckhpd	xmm7,[ebx+0xf0]
		__asm	unpcklpd	xmm5,[ebx+0xf0]	/* a[jp+p14] */
		__asm	movaps	[ecx+0x170],xmm7	/* Store hi imag in t16+16 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p14] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p14] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p14]*c14 */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p14]*c14 */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p14]*s14 */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p14]*s14 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt		xmm6,7 free */
		__asm	movaps	[ecx+0x070],xmm5	/* Store it in t16*/
		__asm	movaps	[ecx+0x060],xmm4	/* Store rt in t15*/

		__asm	mov	edx, c6
		/* Real parts: */
		__asm	movaps		xmm6,[eax+0x60]	/* a[j1+p6 ], this is the scratch xmm register  */
		__asm	movaps		xmm4,[eax+0x60]	/* a[j1+p6 ], this is the active  xmm register */
		__asm	unpckhpd	xmm6,[ebx+0x60]	/* a[j2+p6 ] gets read twice */
		__asm	unpcklpd	xmm4,[ebx+0x60]	/* a[jt+p6 ] */
		__asm	movaps	[ecx+0x140],xmm6	/* Store hi real in t13+16 */
		/* Imag parts: */
		__asm	movaps		xmm7,[eax+0x70]
		__asm	movaps		xmm5,[eax+0x70]
		__asm	unpckhpd	xmm7,[ebx+0x70]
		__asm	unpcklpd	xmm5,[ebx+0x70]	/* a[jp+p6 ] */
		__asm	movaps	[ecx+0x150],xmm7	/* Store hi imag in t14+16 */

		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p6] */
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p6] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p6]*c6 */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p6]*c6 */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p6]*s6 */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p6]*s6 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- t14*/
		__asm	subpd	xmm4,xmm7	/* xmm4 <- t13*/
		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy t14*/
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy t13*/

		__asm	subpd	xmm4,[ecx+0x060]	/* ~t15<- t13-rt */
		__asm	subpd	xmm5,[ecx+0x070]	/* ~t16<- t14-it */
		__asm	addpd	xmm6,[ecx+0x060]	/* ~t13<- t13+rt */
		__asm	addpd	xmm7,[ecx+0x070]	/* ~t14<- t14+it */

		/* Finish radix-4 butterfly and store results into temporary-array slots: */
		/*
		~t13=t9 -t5;		~t9 =t9 +t5;
		~t14=t10-t6;		~t10=t10+t6;
		*/
		__asm	subpd	xmm0,xmm6	/*~t13*/
		__asm	subpd	xmm1,xmm7	/*~t14*/
		__asm	movaps	[ecx+0x040],xmm0	/* a[jt+p8 ] <- ~t13*/
		__asm	movaps	[ecx+0x050],xmm1	/* a[jp+p8 ] <- ~t14*/
		__asm	addpd	xmm6,xmm6	/* 2*t13*/
		__asm	addpd	xmm7,xmm7	/* 2*t14*/
		__asm	addpd	xmm6,xmm0	/*~t9 */
		__asm	addpd	xmm7,xmm1	/*~t10*/
		__asm	movaps	[ecx      ],xmm6	/* a[jt    ] <- ~t9 */
		__asm	movaps	[ecx+0x010],xmm7	/* a[jp    ] <- ~t10*/

		/*
		~t15=t11+t8;		~t11=t11-t8;
		~t16=t12-t7;		~t12=t12+t7;
		*/
		__asm	subpd	xmm2,xmm5	/*~t11*/
		__asm	subpd	xmm3,xmm4	/*~t16*/
		__asm	movaps	[ecx+0x020],xmm2	/* a[jt+p4 ] <- ~t11*/
		__asm	movaps	[ecx+0x070],xmm3	/* a[jp+p12] <- ~t16*/
		__asm	addpd	xmm5,xmm5	/* 2*t16*/
		__asm	addpd	xmm4,xmm4	/* 2*t15*/
		__asm	addpd	xmm5,xmm2	/*~t15*/
		__asm	addpd	xmm4,xmm3	/*~t12*/
		__asm	movaps	[ecx+0x060],xmm5	/* a[jt+p12] <- ~t15*/
		__asm	movaps	[ecx+0x030],xmm4	/* a[jp+p4 ] <- ~t12*/

	/***************************************************************************************************************************
	Each of the above 2 radix-4 blocks is unique, but the ensuing 2 blocks
	[operating on the odd-indexed elements from the unpck*pd commands which were stored to temporaries can use a common macro:
	***************************************************************************************************************************/
	/*...Block 3: */
		SSE2_RADIX4_DIF_4TWIDDLE(r17,r21,r19,r23, r17,c1)

	/*...Block 4: */
		SSE2_RADIX4_DIF_4TWIDDLE(r25,r29,r27,r31, r25,c3)

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
		__asm	mov	esi, two
		__asm	mov	eax, r1

		__asm	movaps	xmm0,[eax      ]	/* t1  */
		__asm	movaps	xmm1,[eax+0x010]	/* t2  */
		__asm	movaps	xmm2,[eax+0x080]	/* t9  */
		__asm	movaps	xmm3,[eax+0x090]	/* t10 */

		__asm	subpd	xmm0,[eax+0x080]	/* t9 =t1 -rt */
		__asm	subpd	xmm1,[eax+0x090]	/* t10=t2 -it */
		__asm	addpd	xmm2,[eax      ]	/* t1 =t1 +rt */
		__asm	addpd	xmm3,[eax+0x010]	/* t2 =t2 +it */

		__asm	movaps	xmm4,[eax+0x100]	/* t17 */
		__asm	movaps	xmm5,[eax+0x110]	/* t18 */
		__asm	movaps	xmm6,[eax+0x180]	/* t25 */
		__asm	movaps	xmm7,[eax+0x190]	/* t26 */

		__asm	subpd	xmm4,[eax+0x180]	/* t25=t17-rt */
		__asm	subpd	xmm5,[eax+0x190]	/* t26=t18-it */
		__asm	addpd	xmm6,[eax+0x100]	/* t17=t17+rt */
		__asm	addpd	xmm7,[eax+0x110]	/* t18=t18+it */

		__asm	subpd	xmm2,xmm6		/* t1  <- t1 -t17 */
		__asm	subpd	xmm3,xmm7		/* t2  <- t2 -t18 */
		__asm	movaps	[eax+0x100],xmm2	/* a[jt+p1 ], store in t17 */
		__asm	addpd	xmm6,xmm6		/*          2*t17 */
		__asm	addpd	xmm7,xmm7		/*          2*t18 */
		__asm	addpd	xmm6,xmm2		/* t17 <- t1 +t17 */
		__asm	addpd	xmm7,xmm3		/* t18 <- t2 +t18 */
		__asm	movaps	[eax      ],xmm6	/* a[jt+p0 ], store in t0  */
	__asm	movaps	xmm6,xmm2	/* cpy x */
	__asm	addpd	xmm2,xmm3	/* x+y */
	__asm	subpd	xmm6,xmm3	/* x-y */
	__asm	addpd	xmm3,xmm3	/* 2*y */
	__asm	mulpd	xmm2,xmm6	/* x^2-y^2 */
	__asm	mulpd	xmm3,[eax+0x100]/* 2xy */
	__asm	movaps	xmm6,[eax      ]	/* a[jt+p0 ], reload */
		__asm	movaps	[eax+0x100],xmm2	/* a[jt+p1 ], store in t17 */
		__asm	movaps	[eax+0x110],xmm3	/* a[jp+p1 ], store in t18 */
	__asm	movaps	xmm2,xmm6	/* cpy x */
	__asm	movaps	xmm3,xmm6	/* cpy x */
	__asm	addpd	xmm6,xmm7	/* x+y */
	__asm	subpd	xmm2,xmm7	/* x-y */
	__asm	addpd	xmm7,xmm7	/* 2*y */
	__asm	mulpd	xmm6,xmm2	/* x^2-y^2 */
	__asm	mulpd	xmm7,xmm3	/* 2xy */
//		__asm	movaps	[eax      ],xmm6	/* a[jt+p0 ], store in t0 */
//		__asm	movaps	[eax+0x010],xmm7	/* a[jp+p0 ], store in t1 */

		__asm	subpd	xmm0,xmm5		/* t9  <- t9 -t26 */
		__asm	subpd	xmm1,xmm4		/* t10 <- t10-t25 */
		__asm	addpd	xmm5,xmm5		/*          2*t26 */
		__asm	addpd	xmm4,xmm4		/*          2*t25 */
		__asm	addpd	xmm5,xmm0		/* t26 <- t9 +t26 */
		__asm	addpd	xmm4,xmm1		/* t25 <- t10+t25 */
	__asm	movaps	xmm2,xmm5	/* cpy x */
	__asm	movaps	xmm3,xmm5	/* cpy x */
	__asm	addpd	xmm5,xmm1	/* x+y */
	__asm	subpd	xmm2,xmm1	/* x-y */
	__asm	addpd	xmm1,xmm1	/* 2*y */
	__asm	mulpd	xmm5,xmm2	/* x^2-y^2 */
	__asm	mulpd	xmm1,xmm3	/* 2xy */
//		__asm	movaps	[eax+0x180],xmm5	/* a[jt+p3 ], store in t25 */
//		__asm	movaps	[eax+0x190],xmm1	/* a[jp+p3 ], store in t26 */
	__asm	movaps	xmm2,xmm0	/* cpy x */
	__asm	movaps	xmm3,xmm0	/* cpy x */
	__asm	addpd	xmm0,xmm4	/* x+y */
	__asm	subpd	xmm2,xmm4	/* x-y */
	__asm	addpd	xmm4,xmm4	/* 2*y */
	__asm	mulpd	xmm0,xmm2	/* x^2-y^2 */
	__asm	mulpd	xmm4,xmm3	/* 2xy */
//		__asm	movaps	[eax+0x080],xmm0	/* a[jt+p2 ], store in t9  */
//		__asm	movaps	[eax+0x090],xmm4	/* a[jp+p2 ], store in t10 */

	__asm	movaps	xmm2,[eax+0x100]	/* a[jt+p1 ], reload */
	__asm	movaps	xmm3,[eax+0x110]	/* a[jp+p1 ], reload */
		SSE2_RADIX4_DIT_IN_PLACE_C(xmm6,xmm7,xmm2,xmm3,xmm0,xmm4,xmm5,xmm1)

	/*...Block 3: t5,13,21,29 */
		__asm	mov	eax, r5

		__asm	movaps	xmm0,[eax      ]	/* t5  */
		__asm	movaps	xmm1,[eax+0x010]	/* t6  */
		__asm	movaps	xmm2,[eax+0x080]	/* t13 */
		__asm	movaps	xmm3,[eax+0x090]	/* t14 */

		__asm	subpd	xmm0,[eax+0x090]	/* t5 =t5 -t14*/
		__asm	subpd	xmm1,[eax+0x080]	/* t14=t6 -t13*/
		__asm	addpd	xmm2,[eax+0x010]	/* t6 =t13+t6 */
		__asm	addpd	xmm3,[eax      ]	/* t13=t14+t5 */
		__asm	mov	ebx, isrt2

		__asm	movaps	xmm4,[eax+0x100]	/* t21 */
		__asm	movaps	xmm5,[eax+0x110]	/* t22 */
		__asm	movaps	xmm6,[eax+0x180]	/* t29 */
		__asm	movaps	xmm7,[eax+0x190]	/* t30 */

		__asm	subpd	xmm4,[eax+0x110]	/* t21-t22 */
		__asm	addpd	xmm5,[eax+0x100]	/* t22+t21 */
		__asm	mulpd	xmm4,[ebx]	/* t21 = (t21-t22)*ISRT2 */
		__asm	mulpd	xmm5,[ebx]	/* t22 = (t22+t21)*ISRT2 */

		__asm	addpd	xmm6,[eax+0x190]	/* t29+t30 */
		__asm	subpd	xmm7,[eax+0x180]	/* t30-t29 */
		__asm	mulpd	xmm6,[ebx]	/*  rt = (t29+t30)*ISRT2 */
		__asm	mulpd	xmm7,[ebx]	/*  it = (t30-t29)*ISRT2 */

		__asm	subpd	xmm4,xmm6		/* t21=t21-rt */
		__asm	subpd	xmm5,xmm7		/* t22=t22-it */
		__asm	addpd	xmm6,xmm6		/*      2* rt */
		__asm	addpd	xmm7,xmm7		/*      2* it */
		__asm	addpd	xmm6,xmm4		/* t29=t21+rt */
		__asm	addpd	xmm7,xmm5		/* t30=t22+it */

		__asm	subpd	xmm0,xmm4		/* t5 -t21 */
		__asm	subpd	xmm2,xmm5		/* t6 -t22 */
		__asm	addpd	xmm4,xmm4		/*   2*t21 */
		__asm	addpd	xmm5,xmm5		/*   2*t22 */
		__asm	movaps	[eax+0x100],xmm0	/* a[jt+p1 ] */
		__asm	addpd	xmm4,xmm0		/* t5 +t21 */
		__asm	addpd	xmm5,xmm2		/* t6 +t22 */
		__asm	movaps	[eax      ],xmm4	/* a[jt+p0 ] */
	__asm	movaps	xmm4,xmm0	/* cpy x */
	__asm	addpd	xmm0,xmm2	/* x+y */
	__asm	subpd	xmm4,xmm2	/* x-y */
	__asm	addpd	xmm2,xmm2	/* 2*y */
	__asm	mulpd	xmm0,xmm4	/* x^2-y^2 */
	__asm	mulpd	xmm2,[eax+0x100]/* 2xy */
	__asm	movaps	xmm4,[eax      ]	/* a[jt+p0 ], reload */
		__asm	movaps	[eax+0x100],xmm0	/* a[jt+p1 ] */
		__asm	movaps	[eax+0x110],xmm2	/* a[jp+p1 ] */
	__asm	movaps	xmm0,xmm4	/* cpy x */
	__asm	movaps	xmm2,xmm4	/* cpy x */
	__asm	addpd	xmm4,xmm5	/* x+y */
	__asm	subpd	xmm0,xmm5	/* x-y */
	__asm	addpd	xmm5,xmm5	/* 2*y */
	__asm	mulpd	xmm4,xmm0	/* x^2-y^2 */
	__asm	mulpd	xmm5,xmm2	/* 2xy */
//		__asm	movaps	[eax      ],xmm4	/* a[jt+p0 ] */
//		__asm	movaps	[eax+0x010],xmm5	/* a[jp+p0 ] */

		__asm	subpd	xmm3,xmm7		/* t13-t30 */
		__asm	subpd	xmm1,xmm6		/* t14-t29 */
		__asm	addpd	xmm7,xmm7		/*   2*t30 */
		__asm	addpd	xmm6,xmm6		/*   2*t29 */
		__asm	addpd	xmm7,xmm3		/* t13+t30 */
		__asm	addpd	xmm6,xmm1		/* t14+t29 */
	__asm	movaps	xmm0,xmm7	/* cpy x */
	__asm	movaps	xmm2,xmm7	/* cpy x */
	__asm	addpd	xmm7,xmm1	/* x+y */
	__asm	subpd	xmm0,xmm1	/* x-y */
	__asm	addpd	xmm1,xmm1	/* 2*y */
	__asm	mulpd	xmm7,xmm0	/* x^2-y^2 */
	__asm	mulpd	xmm1,xmm2	/* 2xy */
//		__asm	movaps	[eax+0x180],xmm7	/* a[jt+p3 ] */
//		__asm	movaps	[eax+0x190],xmm1	/* a[jp+p3 ] */
	__asm	movaps	xmm0,xmm3	/* cpy x */
	__asm	movaps	xmm2,xmm3	/* cpy x */
	__asm	addpd	xmm3,xmm6	/* x+y */
	__asm	subpd	xmm0,xmm6	/* x-y */
	__asm	addpd	xmm6,xmm6	/* 2*y */
	__asm	mulpd	xmm3,xmm0	/* x^2-y^2 */
	__asm	mulpd	xmm6,xmm2	/* 2xy */
//		__asm	movaps	[eax+0x080],xmm3	/* a[jt+p2 ] */
//		__asm	movaps	[eax+0x090],xmm6	/* a[jp+p2 ] */

	__asm	movaps	xmm0,[eax+0x100]	/* a[jt+p1 ], reload */
	__asm	movaps	xmm2,[eax+0x110]	/* a[jp+p1 ], reload */
		SSE2_RADIX4_DIT_IN_PLACE_C(xmm4,xmm5,xmm0,xmm2,xmm3,xmm6,xmm7,xmm1)

	/*...Block 2: t3,11,19,27 */
		__asm	mov	eax, r3
		__asm	mov	ebx, cc0

		__asm	movaps	xmm4,[eax+0x100]	/* t19 */		__asm	movaps	xmm6,[eax+0x180]	/* t27 */
		__asm	movaps	xmm5,[eax+0x110]	/* t20 */		__asm	movaps	xmm7,[eax+0x190]	/* t28 */
		__asm	movaps	xmm0,[eax+0x100]	/* copy t19 */	__asm	movaps	xmm2,[eax+0x180]	/* copy t27 */
		__asm	movaps	xmm1,[eax+0x110]	/* copy t20 */	__asm	movaps	xmm3,[eax+0x190]	/* copy t28 */

		__asm	mulpd	xmm4,[ebx     ]	/* t19*c */			__asm	mulpd	xmm6,[ebx+0x10]	/* t27*s */
		__asm	mulpd	xmm1,[ebx+0x10]	/* t20*s */			__asm	mulpd	xmm3,[ebx     ]	/* t28*c */
		__asm	mulpd	xmm5,[ebx     ]	/* t20*c */			__asm	mulpd	xmm7,[ebx+0x10]	/* t28*s */
		__asm	mulpd	xmm0,[ebx+0x10]	/* t19*s */			__asm	mulpd	xmm2,[ebx     ]	/* t27*c */
		__asm	subpd	xmm4,xmm1	/* ~t19 */				__asm	subpd	xmm6,xmm3	/* rt */
		__asm	addpd	xmm5,xmm0	/* ~t20 */				__asm	addpd	xmm7,xmm2	/* it */

		__asm	subpd	xmm4,xmm6		/*~t27=t19-rt */
		__asm	subpd	xmm5,xmm7		/*~t28=t20-it */
		__asm	addpd	xmm6,xmm6		/*      2* rt */
		__asm	addpd	xmm7,xmm7		/*      2* it */
		__asm	addpd	xmm6,xmm4		/*~t19=t19+rt */
		__asm	addpd	xmm7,xmm5		/*~t20=t20+it */

		__asm	mov	ebx, isrt2
		__asm	movaps	xmm2,[eax+0x080]	/* t11 */
		__asm	movaps	xmm3,[eax+0x090]	/* t12 */
		__asm	subpd	xmm2,[eax+0x090]	/* t11-t12 */
		__asm	addpd	xmm3,[eax+0x080]	/* t12+t11 */
		__asm	mulpd	xmm2,[ebx]	/* rt = (t11-t12)*ISRT2 */
		__asm	mulpd	xmm3,[ebx]	/* it = (t12+t11)*ISRT2 */

		__asm	movaps	xmm0,[eax      ]	/* t3  */
		__asm	movaps	xmm1,[eax+0x010]	/* t4  */

		__asm	subpd	xmm0,xmm2			/*~t11=t3 -rt */
		__asm	subpd	xmm1,xmm3			/*~t12=t4 -it */
		__asm	addpd	xmm2,[eax      ]	/*~t3 =rt +t3 */
		__asm	addpd	xmm3,[eax+0x010]	/*~t4 =it +t4 */

		__asm	subpd	xmm2,xmm6		/* t3 -t19 */
		__asm	subpd	xmm3,xmm7		/* t4 -t20 */
		__asm	addpd	xmm6,xmm6		/*   2*t19 */
		__asm	addpd	xmm7,xmm7		/*   2*t20 */
		__asm	movaps	[eax+0x100],xmm2	/* a[jt+p1 ] */
		__asm	addpd	xmm6,xmm2		/* t3 +t19 */
		__asm	addpd	xmm7,xmm3		/* t4 +t20 */
		__asm	movaps	[eax      ],xmm6	/* a[jt+p0 ] */
	__asm	movaps	xmm6,xmm2	/* cpy x */
	__asm	addpd	xmm2,xmm3	/* x+y */
	__asm	subpd	xmm6,xmm3	/* x-y */
	__asm	addpd	xmm3,xmm3	/* 2*y */
	__asm	mulpd	xmm2,xmm6	/* x^2-y^2 */
	__asm	mulpd	xmm3,[eax+0x100]/* 2xy */
	__asm	movaps	xmm6,[eax      ]	/* a[jt+p0 ], reload */
		__asm	movaps	[eax+0x100],xmm2	/* a[jt+p1 ] */
		__asm	movaps	[eax+0x110],xmm3	/* a[jp+p1 ] */
	__asm	movaps	xmm2,xmm6	/* cpy x */
	__asm	movaps	xmm3,xmm6	/* cpy x */
	__asm	addpd	xmm6,xmm7	/* x+y */
	__asm	subpd	xmm2,xmm7	/* x-y */
	__asm	addpd	xmm7,xmm7	/* 2*y */
	__asm	mulpd	xmm6,xmm2	/* x^2-y^2 */
	__asm	mulpd	xmm7,xmm3	/* 2xy */
//		__asm	movaps	[eax      ],xmm6	/* a[jt+p0 ] */
//		__asm	movaps	[eax+0x010],xmm7	/* a[jp+p0 ] */

		__asm	subpd	xmm0,xmm5		/* t11-t28 */
		__asm	subpd	xmm1,xmm4		/* t12-t27 */
		__asm	addpd	xmm5,xmm5		/*          2*t28 */
		__asm	addpd	xmm4,xmm4		/*          2*t27 */
		__asm	addpd	xmm5,xmm0		/* t11+t28 */
		__asm	addpd	xmm4,xmm1		/* t12+t27 */
	__asm	movaps	xmm2,xmm5	/* cpy x */
	__asm	movaps	xmm3,xmm5	/* cpy x */
	__asm	addpd	xmm5,xmm1	/* x+y */
	__asm	subpd	xmm2,xmm1	/* x-y */
	__asm	addpd	xmm1,xmm1	/* 2*y */
	__asm	mulpd	xmm5,xmm2	/* x^2-y^2 */
	__asm	mulpd	xmm1,xmm3	/* 2xy */
//		__asm	movaps	[eax+0x180],xmm5	/* a[jt+p3 ] */
//		__asm	movaps	[eax+0x190],xmm1	/* a[jp+p3 ] */
	__asm	movaps	xmm2,xmm0	/* cpy x */
	__asm	movaps	xmm3,xmm0	/* cpy x */
	__asm	addpd	xmm0,xmm4	/* x+y */
	__asm	subpd	xmm2,xmm4	/* x-y */
	__asm	addpd	xmm4,xmm4	/* 2*y */
	__asm	mulpd	xmm0,xmm2	/* x^2-y^2 */
	__asm	mulpd	xmm4,xmm3	/* 2xy */
//		__asm	movaps	[eax+0x080],xmm0	/* a[jt+p2 ] */
//		__asm	movaps	[eax+0x090],xmm4	/* a[jp+p2 ] */

	__asm	movaps	xmm2,[eax+0x100]	/* a[jt+p1 ], reload */
	__asm	movaps	xmm3,[eax+0x110]	/* a[jp+p1 ], reload */
		SSE2_RADIX4_DIT_IN_PLACE_C(xmm6,xmm7,xmm2,xmm3,xmm0,xmm4,xmm5,xmm1)

	/*...Block 4: t7,15,23,31 */
		__asm	mov	eax, r7
		__asm	mov	ebx, cc0

		__asm	movaps	xmm4,[eax+0x100]	/* t23 */			__asm	movaps	xmm6,[eax+0x180]	/* t31 */
		__asm	movaps	xmm5,[eax+0x110]	/* t24 */			__asm	movaps	xmm7,[eax+0x190]	/* t32 */
		__asm	movaps	xmm0,[eax+0x100]	/* copy t23 */		__asm	movaps	xmm2,[eax+0x180]	/* copy t31 */
		__asm	movaps	xmm1,[eax+0x110]	/* copy t24 */		__asm	movaps	xmm3,[eax+0x190]	/* copy t32 */

		__asm	mulpd	xmm4,[ebx+0x10]	/* t23*s */			__asm	mulpd	xmm6,[ebx     ]	/* t31*c */
		__asm	mulpd	xmm1,[ebx     ]	/* t24*c */			__asm	mulpd	xmm3,[ebx+0x10]	/* t32*s */
		__asm	mulpd	xmm5,[ebx+0x10]	/* t24*s */			__asm	mulpd	xmm7,[ebx     ]	/* t32*c */
		__asm	mulpd	xmm0,[ebx     ]	/* t23*c */			__asm	mulpd	xmm2,[ebx+0x10]	/* t31*s */
		__asm	subpd	xmm4,xmm1	/* ~t23 */				__asm	subpd	xmm6,xmm3	/* rt */
		__asm	addpd	xmm5,xmm0	/* ~t24 */				__asm	addpd	xmm7,xmm2	/* it */

		__asm	subpd	xmm4,xmm6		/*~t23=t23-rt */
		__asm	subpd	xmm5,xmm7		/*~t24=t24-it */
		__asm	addpd	xmm6,xmm6		/*      2* rt */
		__asm	addpd	xmm7,xmm7		/*      2* it */
		__asm	addpd	xmm6,xmm4		/*~t31=t23+rt */
		__asm	addpd	xmm7,xmm5		/*~t32=t24+it */

		__asm	mov	ebx, isrt2
		__asm	movaps	xmm2,[eax+0x080]	/* t15 */
		__asm	movaps	xmm3,[eax+0x090]	/* t16 */
		__asm	addpd	xmm2,[eax+0x090]	/* t15+t16 */
		__asm	subpd	xmm3,[eax+0x080]	/* t16-t15 */
		__asm	mulpd	xmm2,[ebx]	/* rt = (t15+t16)*ISRT2 */
		__asm	mulpd	xmm3,[ebx]	/* it = (t16-t15)*ISRT2 */

		__asm	movaps	xmm0,[eax      ]	/* t7  */
		__asm	movaps	xmm1,[eax+0x010]	/* t8  */

		__asm	subpd	xmm0,xmm2			/*~t7 =t7 -rt */
		__asm	subpd	xmm1,xmm3			/*~t8 =t8 -it */
		__asm	addpd	xmm2,[eax      ]	/*~t15=rt +t7 */
		__asm	addpd	xmm3,[eax+0x010]	/*~t16=it +t8 */

		__asm	subpd	xmm0,xmm4		/* t7 -t23 */
		__asm	subpd	xmm1,xmm5		/* t8 -t24 */
		__asm	addpd	xmm4,xmm4		/*   2*t23 */
		__asm	addpd	xmm5,xmm5		/*   2*t24 */
		__asm	movaps	[eax+0x100],xmm0	/* a[jt+p1 ] */
		__asm	addpd	xmm4,xmm0		/* t7 +t23 */
		__asm	addpd	xmm5,xmm1		/* t8 +t24 */
		__asm	movaps	[eax      ],xmm4	/* a[jt+p0 ] */
	__asm	movaps	xmm4,xmm0	/* cpy x */
	__asm	addpd	xmm0,xmm1	/* x+y */
	__asm	subpd	xmm4,xmm1	/* x-y */
	__asm	addpd	xmm1,xmm1	/* 2*y */
	__asm	mulpd	xmm0,xmm4	/* x^2-y^2 */
	__asm	mulpd	xmm1,[eax+0x100]/* 2xy */
	__asm	movaps	xmm4,[eax      ]	/* a[jt+p0 ], reload */
		__asm	movaps	[eax+0x100],xmm0	/* a[jt+p1 ] */
		__asm	movaps	[eax+0x110],xmm1	/* a[jp+p1 ] */
	__asm	movaps	xmm0,xmm4	/* cpy x */
	__asm	movaps	xmm1,xmm4	/* cpy x */
	__asm	addpd	xmm4,xmm5	/* x+y */
	__asm	subpd	xmm0,xmm5	/* x-y */
	__asm	addpd	xmm5,xmm5	/* 2*y */
	__asm	mulpd	xmm4,xmm0	/* x^2-y^2 */
	__asm	mulpd	xmm5,xmm1	/* 2xy */
//		__asm	movaps	[eax      ],xmm4	/* a[jt+p0 ] */
//		__asm	movaps	[eax+0x010],xmm5	/* a[jp+p0 ] */

		__asm	subpd	xmm2,xmm7		/* t15-t32 */
		__asm	subpd	xmm3,xmm6		/* t16-t31 */
		__asm	addpd	xmm7,xmm7		/*   2*t32 */
		__asm	addpd	xmm6,xmm6		/*   2*t31 */
		__asm	addpd	xmm7,xmm2		/* t15+t32 */
		__asm	addpd	xmm6,xmm3		/* t16+t31 */
	__asm	movaps	xmm0,xmm7	/* cpy x */
	__asm	movaps	xmm1,xmm7	/* cpy x */
	__asm	addpd	xmm7,xmm3	/* x+y */
	__asm	subpd	xmm0,xmm3	/* x-y */
	__asm	addpd	xmm3,xmm3	/* 2*y */
	__asm	mulpd	xmm7,xmm0	/* x^2-y^2 */
	__asm	mulpd	xmm3,xmm1	/* 2xy */
//		__asm	movaps	[eax+0x180],xmm7	/* a[jt+p3 ] */
//		__asm	movaps	[eax+0x190],xmm3	/* a[jp+p3 ] */
	__asm	movaps	xmm0,xmm2	/* cpy x */
	__asm	movaps	xmm1,xmm2	/* cpy x */
	__asm	addpd	xmm2,xmm6	/* x+y */
	__asm	subpd	xmm0,xmm6	/* x-y */
	__asm	addpd	xmm6,xmm6	/* 2*y */
	__asm	mulpd	xmm2,xmm0	/* x^2-y^2 */
	__asm	mulpd	xmm6,xmm1	/* 2xy */
//		__asm	movaps	[eax+0x080],xmm2	/* a[jt+p2 ] */
//		__asm	movaps	[eax+0x090],xmm6	/* a[jp+p2 ] */

	__asm	movaps	xmm0,[eax+0x100]	/* a[jt+p1 ], reload */
	__asm	movaps	xmm1,[eax+0x110]	/* a[jp+p1 ], reload */
		SSE2_RADIX4_DIT_IN_PLACE_C(xmm4,xmm5,xmm0,xmm1,xmm2,xmm6,xmm7,xmm3)

	#else	/* GCC-style inline ASM: */

		SSE2_RADIX16_WRAPPER_DIF(add0,add1,r1,r9,r17,r25,isrt2,cc0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15)

	#endif

	#if 0//DEBUG_SSE2
		fprintf(stderr, "Outputs B:\n");
		fprintf(stderr, "radix16_dyadic_square: R1:= %20.5f, %20.5f\n",r1 ->re,r2 ->re);
		fprintf(stderr, "radix16_dyadic_square: R17= %20.5f, %20.5f\n",r17->re,r18->re);
		fprintf(stderr, "radix16_dyadic_square: R9:= %20.5f, %20.5f\n",r9 ->re,r10->re);
		fprintf(stderr, "radix16_dyadic_square: R25= %20.5f, %20.5f\n",r25->re,r26->re);
		fprintf(stderr, "radix16_dyadic_square: R5:= %20.5f, %20.5f\n",r5 ->re,r6 ->re);
		fprintf(stderr, "radix16_dyadic_square: R21= %20.5f, %20.5f\n",r21->re,r22->re);
		fprintf(stderr, "radix16_dyadic_square: R13= %20.5f, %20.5f\n",r13->re,r14->re);
		fprintf(stderr, "radix16_dyadic_square: R29= %20.5f, %20.5f\n",r29->re,r30->re);
		fprintf(stderr, "radix16_dyadic_square: R3:= %20.5f, %20.5f\n",r3 ->re,r4 ->re);
		fprintf(stderr, "radix16_dyadic_square: R19= %20.5f, %20.5f\n",r19->re,r20->re);
		fprintf(stderr, "radix16_dyadic_square: R11= %20.5f, %20.5f\n",r11->re,r12->re);
		fprintf(stderr, "radix16_dyadic_square: R27= %20.5f, %20.5f\n",r27->re,r28->re);
		fprintf(stderr, "radix16_dyadic_square: R7:= %20.5f, %20.5f\n",r7 ->re,r8 ->re);
		fprintf(stderr, "radix16_dyadic_square: R23= %20.5f, %20.5f\n",r23->re,r24->re);
		fprintf(stderr, "radix16_dyadic_square: R15= %20.5f, %20.5f\n",r15->re,r16->re);
		fprintf(stderr, "radix16_dyadic_square: R31= %20.5f, %20.5f\n",r31->re,r32->re);
		fprintf(stderr, "\n");
		exit(0);
	#endif

	#ifdef COMPILER_TYPE_MSVC
	  #if 0	// Roll squarings in with radix-4 DFTs
		__asm	mov	eax, r1
		__asm	movaps	xmm0,[eax      ]	/* x0*/		__asm	movaps	xmm0,[eax+0x100]	/* x8*/
		__asm	movaps	xmm1,[eax+0x010]	/* y0*/		__asm	movaps	xmm1,[eax+0x110]	/* y8*/
		__asm	movaps	xmm2,xmm0		/* cpy x */		__asm	movaps	xmm2,xmm0		/* cpy x */
		__asm	addpd	xmm0,xmm1		/* x+y */		__asm	addpd	xmm0,xmm1		/* x+y */
		__asm	subpd	xmm2,xmm1		/* x-y */		__asm	subpd	xmm2,xmm1		/* x-y */
		__asm	addpd	xmm1,xmm1		/* 2*y */		__asm	addpd	xmm1,xmm1		/* 2*y */
		__asm	mulpd	xmm0,xmm2	/* x^2-y^2 */		__asm	mulpd	xmm0,xmm2	/* x^2-y^2 */
		__asm	mulpd	xmm1,[eax      ]/* 2xy */		__asm	mulpd	xmm1,[eax+0x100]/* 2xy */
		__asm	movaps	[eax      ],xmm0	/* Re */	__asm	movaps	[eax+0x100],xmm0	/* Re */
		__asm	movaps	[eax+0x010],xmm1	/* Im */	__asm	movaps	[eax+0x110],xmm1	/* Im */

		__asm	movaps	xmm0,[eax+0x080]	/* x4*/		__asm	movaps	xmm0,[eax+0x180]	/* xC*/
		__asm	movaps	xmm1,[eax+0x090]	/* y4*/		__asm	movaps	xmm1,[eax+0x190]	/* yC*/
		__asm	movaps	xmm2,xmm0		/* cpy x */		__asm	movaps	xmm2,xmm0		/* cpy x */
		__asm	addpd	xmm0,xmm1		/* x+y */		__asm	addpd	xmm0,xmm1		/* x+y */
		__asm	subpd	xmm2,xmm1		/* x-y */		__asm	subpd	xmm2,xmm1		/* x-y */
		__asm	addpd	xmm1,xmm1		/* 2*y */		__asm	addpd	xmm1,xmm1		/* 2*y */
		__asm	mulpd	xmm0,xmm2	/* x^2-y^2 */		__asm	mulpd	xmm0,xmm2	/* x^2-y^2 */
		__asm	mulpd	xmm1,[eax+0x080]/* 2xy */		__asm	mulpd	xmm1,[eax+0x180]/* 2xy */
		__asm	movaps	[eax+0x080],xmm0	/* Re */	__asm	movaps	[eax+0x180],xmm0	/* Re */
		__asm	movaps	[eax+0x090],xmm1	/* Im */	__asm	movaps	[eax+0x190],xmm1	/* Im */

		__asm	movaps	xmm0,[eax+0x040]	/* x2*/		__asm	movaps	xmm3,[eax+0x140]
		__asm	movaps	xmm1,[eax+0x050]	/* y2*/		__asm	movaps	xmm4,[eax+0x150]
		__asm	movaps	xmm2,xmm0		/* cpy x */		__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1		/* x+y */		__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1		/* x-y */		__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1		/* 2*y */		__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2	/* x^2-y^2 */		__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x040]/* 2xy */		__asm	mulpd	xmm4,[eax+0x140]
		__asm	movaps	[eax+0x040],xmm0	/* Re */	__asm	movaps	[eax+0x140],xmm3
		__asm	movaps	[eax+0x050],xmm1	/* Im */	__asm	movaps	[eax+0x150],xmm4

		__asm	movaps	xmm0,[eax+0x0c0]	/* x6*/		__asm	movaps	xmm3,[eax+0x1c0]
		__asm	movaps	xmm1,[eax+0x0d0]	/* y6*/		__asm	movaps	xmm4,[eax+0x1d0]
		__asm	movaps	xmm2,xmm0		/* cpy x */		__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1		/* x+y */		__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1		/* x-y */		__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1		/* 2*y */		__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2	/* x^2-y^2 */		__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x0c0]/* 2xy */		__asm	mulpd	xmm4,[eax+0x1c0]
		__asm	movaps	[eax+0x0c0],xmm0	/* Re */	__asm	movaps	[eax+0x1c0],xmm3
		__asm	movaps	[eax+0x0d0],xmm1	/* Im */	__asm	movaps	[eax+0x1d0],xmm4

		__asm	movaps	xmm0,[eax+0x020]				__asm	movaps	xmm3,[eax+0x120]
		__asm	movaps	xmm1,[eax+0x030]				__asm	movaps	xmm4,[eax+0x130]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x020]				__asm	mulpd	xmm4,[eax+0x120]
		__asm	movaps	[eax+0x020],xmm0				__asm	movaps	[eax+0x120],xmm3
		__asm	movaps	[eax+0x030],xmm1				__asm	movaps	[eax+0x130],xmm4

		__asm	movaps	xmm0,[eax+0x0a0]				__asm	movaps	xmm3,[eax+0x1a0]
		__asm	movaps	xmm1,[eax+0x0b0]				__asm	movaps	xmm4,[eax+0x1b0]
		__asm	movaps	xmm2,xmm0						__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1						__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1						__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1						__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2						__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x0a0]				__asm	mulpd	xmm4,[eax+0x1a0]
		__asm	movaps	[eax+0x0a0],xmm0				__asm	movaps	[eax+0x1a0],xmm3
		__asm	movaps	[eax+0x0b0],xmm1				__asm	movaps	[eax+0x1b0],xmm4

		__asm	movaps	xmm0,[eax+0x060]	/* x */		__asm	movaps	xmm3,[eax+0x160]
		__asm	movaps	xmm1,[eax+0x070]	/* y */		__asm	movaps	xmm4,[eax+0x170]
		__asm	movaps	xmm2,xmm0		/* cpy x */		__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1		/* x+y */		__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1		/* x-y */		__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1		/* 2*y */		__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2	/* x^2-y^2 */		__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x060]/* 2xy */		__asm	mulpd	xmm4,[eax+0x160]
		__asm	movaps	[eax+0x060],xmm0	/* Re */	__asm	movaps	[eax+0x160],xmm3
		__asm	movaps	[eax+0x070],xmm1	/* Im */	__asm	movaps	[eax+0x170],xmm4

		__asm	movaps	xmm0,[eax+0x0e0]	/* xE*/		__asm	movaps	xmm3,[eax+0x1e0]
		__asm	movaps	xmm1,[eax+0x0f0]	/* yE*/		__asm	movaps	xmm4,[eax+0x1f0]
		__asm	movaps	xmm2,xmm0		/* cpy x */		__asm	movaps	xmm5,xmm3
		__asm	addpd	xmm0,xmm1		/* x+y */		__asm	addpd	xmm3,xmm4
		__asm	subpd	xmm2,xmm1		/* x-y */		__asm	subpd	xmm5,xmm4
		__asm	addpd	xmm1,xmm1		/* 2*y */		__asm	addpd	xmm4,xmm4
		__asm	mulpd	xmm0,xmm2	/* x^2-y^2 */		__asm	mulpd	xmm3,xmm5
		__asm	mulpd	xmm1,[eax+0x0e0]/* 2xy */		__asm	mulpd	xmm4,[eax+0x1e0]
		__asm	movaps	[eax+0x0e0],xmm0	/* Re */	__asm	movaps	[eax+0x1e0],xmm3
		__asm	movaps	[eax+0x0f0],xmm1	/* Im */	__asm	movaps	[eax+0x1f0],xmm4
	  #endif
	#else	/* GCC-style inline ASM: */

	  #if OS_BITS == 32

		__asm__ volatile (\
			"movl	%[__r1],%%eax		\n\t"\
			"/* z0^2: */				\n\t"\
			"movaps	     (%%eax),%%xmm0	\n\t"\
			"movaps	0x010(%%eax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	     (%%eax),%%xmm1	\n\t"\
			"movaps	%%xmm0,     (%%eax)	\n\t"\
			"movaps	%%xmm1,0x010(%%eax)	\n\t"\
			"/* z1^2: */				\n\t"\
			"movaps	0x020(%%eax),%%xmm0	\n\t"\
			"movaps	0x030(%%eax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x020(%%eax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x020(%%eax)	\n\t"\
			"movaps	%%xmm1,0x030(%%eax)	\n\t"\
			"/* z2^2: */				\n\t"\
			"movaps	0x040(%%eax),%%xmm0	\n\t"\
			"movaps	0x050(%%eax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x050(%%eax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x040(%%eax)	\n\t"\
			"movaps	%%xmm1,0x050(%%eax)	\n\t"\
			"/* z3^2: */				\n\t"\
			"movaps	0x060(%%eax),%%xmm0	\n\t"\
			"movaps	0x070(%%eax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x060(%%eax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x060(%%eax)	\n\t"\
			"movaps	%%xmm1,0x070(%%eax)	\n\t"\
			"/* z4^2: */				\n\t"\
			"movaps	0x080(%%eax),%%xmm0	\n\t"\
			"movaps	0x090(%%eax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x080(%%eax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x080(%%eax)	\n\t"\
			"movaps	%%xmm1,0x090(%%eax)	\n\t"\
			"/* z5^2: */				\n\t"\
			"movaps	0x0a0(%%eax),%%xmm0	\n\t"\
			"movaps	0x0b0(%%eax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x0a0(%%eax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x0a0(%%eax)	\n\t"\
			"movaps	%%xmm1,0x0b0(%%eax)	\n\t"\
			"/* z6^2: */				\n\t"\
			"movaps	0x0c0(%%eax),%%xmm0	\n\t"\
			"movaps	0x0d0(%%eax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x0c0(%%eax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x0c0(%%eax)	\n\t"\
			"movaps	%%xmm1,0x0d0(%%eax)	\n\t"\
			"/* z7^2: */				\n\t"\
			"movaps	0x0e0(%%eax),%%xmm0	\n\t"\
			"movaps	0x0f0(%%eax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x0e0(%%eax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x0e0(%%eax)	\n\t"\
			"movaps	%%xmm1,0x0f0(%%eax)	\n\t"\
			"/* z8^2: */				\n\t"\
			"movaps	0x100(%%eax),%%xmm0	\n\t"\
			"movaps	0x110(%%eax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x100(%%eax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x100(%%eax)	\n\t"\
			"movaps	%%xmm1,0x110(%%eax)	\n\t"\
			"/* z9^2: */				\n\t"\
			"movaps	0x120(%%eax),%%xmm0	\n\t"\
			"movaps	0x130(%%eax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x120(%%eax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x120(%%eax)	\n\t"\
			"movaps	%%xmm1,0x130(%%eax)	\n\t"\
			"/* zA^2: */				\n\t"\
			"movaps	0x140(%%eax),%%xmm0	\n\t"\
			"movaps	0x150(%%eax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x150(%%eax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x140(%%eax)	\n\t"\
			"movaps	%%xmm1,0x150(%%eax)	\n\t"\
			"/* zB^2: */				\n\t"\
			"movaps	0x160(%%eax),%%xmm0	\n\t"\
			"movaps	0x170(%%eax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x160(%%eax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x160(%%eax)	\n\t"\
			"movaps	%%xmm1,0x170(%%eax)	\n\t"\
			"/* zC^2: */				\n\t"\
			"movaps	0x180(%%eax),%%xmm0	\n\t"\
			"movaps	0x190(%%eax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x180(%%eax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x180(%%eax)	\n\t"\
			"movaps	%%xmm1,0x190(%%eax)	\n\t"\
			"/* zD^2: */				\n\t"\
			"movaps	0x1a0(%%eax),%%xmm0	\n\t"\
			"movaps	0x1b0(%%eax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x1a0(%%eax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x1a0(%%eax)	\n\t"\
			"movaps	%%xmm1,0x1b0(%%eax)	\n\t"\
			"/* zE^2: */				\n\t"\
			"movaps	0x1c0(%%eax),%%xmm0	\n\t"\
			"movaps	0x1d0(%%eax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x1c0(%%eax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x1c0(%%eax)	\n\t"\
			"movaps	%%xmm1,0x1d0(%%eax)	\n\t"\
			"/* zF^2: */				\n\t"\
			"movaps	0x1e0(%%eax),%%xmm0	\n\t"\
			"movaps	0x1f0(%%eax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x1e0(%%eax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x1e0(%%eax)	\n\t"\
			"movaps	%%xmm1,0x1f0(%%eax)	\n\t"\
			:					// outputs: none
			: [__r1] "m" (r1)	// All inputs from memory addresses here
			: "eax"				// Clobbered registers
		);

	  #else

		__asm__ volatile (\
			"movq	%[__r1],%%rax		\n\t"\
			"/* z0^2: */				\n\t"\
			"movaps	     (%%rax),%%xmm0	\n\t"\
			"movaps	0x010(%%rax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	     (%%rax),%%xmm1	\n\t"\
			"movaps	%%xmm0,     (%%rax)	\n\t"\
			"movaps	%%xmm1,0x010(%%rax)	\n\t"\
			"/* z1^2: */				\n\t"\
			"movaps	0x020(%%rax),%%xmm0	\n\t"\
			"movaps	0x030(%%rax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x020(%%rax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x020(%%rax)	\n\t"\
			"movaps	%%xmm1,0x030(%%rax)	\n\t"\
			"/* z2^2: */				\n\t"\
			"movaps	0x040(%%rax),%%xmm0	\n\t"\
			"movaps	0x050(%%rax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x050(%%rax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x040(%%rax)	\n\t"\
			"movaps	%%xmm1,0x050(%%rax)	\n\t"\
			"/* z3^2: */				\n\t"\
			"movaps	0x060(%%rax),%%xmm0	\n\t"\
			"movaps	0x070(%%rax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x060(%%rax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x060(%%rax)	\n\t"\
			"movaps	%%xmm1,0x070(%%rax)	\n\t"\
			"/* z4^2: */				\n\t"\
			"movaps	0x080(%%rax),%%xmm0	\n\t"\
			"movaps	0x090(%%rax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x080(%%rax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x080(%%rax)	\n\t"\
			"movaps	%%xmm1,0x090(%%rax)	\n\t"\
			"/* z5^2: */				\n\t"\
			"movaps	0x0a0(%%rax),%%xmm0	\n\t"\
			"movaps	0x0b0(%%rax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x0a0(%%rax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x0a0(%%rax)	\n\t"\
			"movaps	%%xmm1,0x0b0(%%rax)	\n\t"\
			"/* z6^2: */				\n\t"\
			"movaps	0x0c0(%%rax),%%xmm0	\n\t"\
			"movaps	0x0d0(%%rax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x0c0(%%rax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x0c0(%%rax)	\n\t"\
			"movaps	%%xmm1,0x0d0(%%rax)	\n\t"\
			"/* z7^2: */				\n\t"\
			"movaps	0x0e0(%%rax),%%xmm0	\n\t"\
			"movaps	0x0f0(%%rax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x0e0(%%rax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x0e0(%%rax)	\n\t"\
			"movaps	%%xmm1,0x0f0(%%rax)	\n\t"\
			"/* z8^2: */				\n\t"\
			"movaps	0x100(%%rax),%%xmm0	\n\t"\
			"movaps	0x110(%%rax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x100(%%rax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x100(%%rax)	\n\t"\
			"movaps	%%xmm1,0x110(%%rax)	\n\t"\
			"/* z9^2: */				\n\t"\
			"movaps	0x120(%%rax),%%xmm0	\n\t"\
			"movaps	0x130(%%rax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x120(%%rax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x120(%%rax)	\n\t"\
			"movaps	%%xmm1,0x130(%%rax)	\n\t"\
			"/* zA^2: */				\n\t"\
			"movaps	0x140(%%rax),%%xmm0	\n\t"\
			"movaps	0x150(%%rax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x150(%%rax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x140(%%rax)	\n\t"\
			"movaps	%%xmm1,0x150(%%rax)	\n\t"\
			"/* zB^2: */				\n\t"\
			"movaps	0x160(%%rax),%%xmm0	\n\t"\
			"movaps	0x170(%%rax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x160(%%rax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x160(%%rax)	\n\t"\
			"movaps	%%xmm1,0x170(%%rax)	\n\t"\
			"/* zC^2: */				\n\t"\
			"movaps	0x180(%%rax),%%xmm0	\n\t"\
			"movaps	0x190(%%rax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x180(%%rax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x180(%%rax)	\n\t"\
			"movaps	%%xmm1,0x190(%%rax)	\n\t"\
			"/* zD^2: */				\n\t"\
			"movaps	0x1a0(%%rax),%%xmm0	\n\t"\
			"movaps	0x1b0(%%rax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x1a0(%%rax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x1a0(%%rax)	\n\t"\
			"movaps	%%xmm1,0x1b0(%%rax)	\n\t"\
			"/* zE^2: */				\n\t"\
			"movaps	0x1c0(%%rax),%%xmm0	\n\t"\
			"movaps	0x1d0(%%rax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x1c0(%%rax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x1c0(%%rax)	\n\t"\
			"movaps	%%xmm1,0x1d0(%%rax)	\n\t"\
			"/* zF^2: */				\n\t"\
			"movaps	0x1e0(%%rax),%%xmm0	\n\t"\
			"movaps	0x1f0(%%rax),%%xmm1	\n\t"\
			"movaps	      %%xmm0,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm0	\n\t"\
			"subpd	      %%xmm1,%%xmm2	\n\t"\
			"addpd	      %%xmm1,%%xmm1	\n\t"\
			"mulpd	      %%xmm2,%%xmm0	\n\t"\
			"mulpd	0x1e0(%%rax),%%xmm1	\n\t"\
			"movaps	%%xmm0,0x1e0(%%rax)	\n\t"\
			"movaps	%%xmm1,0x1f0(%%rax)	\n\t"\
			:					// outputs: none
			: [__r1] "m" (r1)	// All inputs from memory addresses here
			: "rax","xmm0","xmm1","xmm2"	// Clobbered registers
		);

	  #endif

	#endif

	#if 0//DEBUG_SSE2
		fprintf(stderr, "Outputs B^2:\n");
		fprintf(stderr, "radix16_dyadic_square: R1:= %20.5f, %20.5f\n",r1 ->re,r2 ->re);
		fprintf(stderr, "radix16_dyadic_square: R17= %20.5f, %20.5f\n",r17->re,r18->re);
		fprintf(stderr, "radix16_dyadic_square: R9:= %20.5f, %20.5f\n",r9 ->re,r10->re);
		fprintf(stderr, "radix16_dyadic_square: R25= %20.5f, %20.5f\n",r25->re,r26->re);
		fprintf(stderr, "radix16_dyadic_square: R5:= %20.5f, %20.5f\n",r5 ->re,r6 ->re);
		fprintf(stderr, "radix16_dyadic_square: R21= %20.5f, %20.5f\n",r21->re,r22->re);
		fprintf(stderr, "radix16_dyadic_square: R13= %20.5f, %20.5f\n",r13->re,r14->re);
		fprintf(stderr, "radix16_dyadic_square: R29= %20.5f, %20.5f\n",r29->re,r30->re);
		fprintf(stderr, "radix16_dyadic_square: R3:= %20.5f, %20.5f\n",r3 ->re,r4 ->re);
		fprintf(stderr, "radix16_dyadic_square: R19= %20.5f, %20.5f\n",r19->re,r20->re);
		fprintf(stderr, "radix16_dyadic_square: R11= %20.5f, %20.5f\n",r11->re,r12->re);
		fprintf(stderr, "radix16_dyadic_square: R27= %20.5f, %20.5f\n",r27->re,r28->re);
		fprintf(stderr, "radix16_dyadic_square: R7:= %20.5f, %20.5f\n",r7 ->re,r8 ->re);
		fprintf(stderr, "radix16_dyadic_square: R23= %20.5f, %20.5f\n",r23->re,r24->re);
		fprintf(stderr, "radix16_dyadic_square: R15= %20.5f, %20.5f\n",r15->re,r16->re);
		fprintf(stderr, "radix16_dyadic_square: R31= %20.5f, %20.5f\n",r31->re,r32->re);
		fprintf(stderr, "\n");
		exit(0);
	#endif

	/*...And do an inverse DIT radix-16 pass on the squared-data blocks. */

	#if defined(COMPILER_TYPE_MSVC)

	/*************************************************************/
	/*                  1st set of inputs:                       */
	/*************************************************************/
	  #if 0
	/*...Block 1: */
		/* eax,ebx,ecx,edx = r1,r17,r9,r25: */
		__asm	mov eax, r1
		__asm	mov ebx, eax
		__asm	mov ecx, eax
		__asm	mov edx, eax
		__asm	add ebx, 0x100
		__asm	add ecx, 0x080
		__asm	add edx, 0x180
		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		SSE2_RADIX4_DIT_IN_PLACE_B()

	/*...Block 2: */
		/* eax,ebx,ecx,edx = r5,r21,r13,r29: */
		__asm	add eax, 0x040
		__asm	add ebx, 0x040
		__asm	add ecx, 0x040
		__asm	add edx, 0x040
		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		SSE2_RADIX4_DIT_IN_PLACE_B()

	/*...Block 3: */
		/* eax,ebx,ecx,edx = r3,r19,r11,r27: */
		__asm	sub eax, 0x020
		__asm	sub ebx, 0x020
		__asm	sub ecx, 0x020
		__asm	sub edx, 0x020
		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		SSE2_RADIX4_DIT_IN_PLACE_B()

	/*...Block 4: */
		/* eax,ebx,ecx,edx = r7,r23,r15,r31: */
		__asm	add eax, 0x040
		__asm	add ebx, 0x040
		__asm	add ecx, 0x040
		__asm	add edx, 0x040
		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		SSE2_RADIX4_DIT_IN_PLACE_B()
	  #endif
	/****************************************************************************************************
	!...and now do four more radix-4 transforms, including the internal and external twiddle factors.   !
	!   Write even-index 16-byte output pairs to a[j1], odd-index to a[j2], unpack same as on inputs.   !
	!   We do the last 2 radix-4 blocks first, to make the unpack-interleaving of outputs smoother.     !
	****************************************************************************************************/

	/* Main-array addresses still in add0,1, no need to re-init:
		add0 = &a[j1pad];
		add1 = &a[j2pad];
	*/
		__asm	mov	eax, r9
		__asm	mov	ebx, isrt2
		__asm	mov	ecx, cc0

		__asm	movaps	xmm4,[eax+0x020]	/* t19/r11 */			__asm	movaps	xmm0,[eax+0x060]	/* t27/r15 */
		__asm	movaps	xmm5,[eax+0x030]	/* t20/r12 */			__asm	movaps	xmm1,[eax+0x070]	/* t28/r16 */
		__asm	movaps	xmm6,[eax+0x020]	/* xmm2 <- cpy t19 */	__asm	movaps	xmm2,[eax+0x060]	/* xmm6 <- cpy t27 */
		__asm	movaps	xmm7,[eax+0x030]	/* xmm3 <- cpy t20 */	__asm	movaps	xmm3,[eax+0x070]	/* xmm7 <- cpy t28 */

		__asm	mulpd	xmm4,[ecx     ]	/* t19*c */					__asm	mulpd	xmm0,[ecx+0x10]	/* t27*s */
		__asm	mulpd	xmm5,[ecx     ]	/* t20*c */					__asm	mulpd	xmm1,[ecx+0x10]	/* t28*s */
		__asm	mulpd	xmm6,[ecx+0x10]	/* t19*s */					__asm	mulpd	xmm2,[ecx     ]	/* t27*c */
		__asm	mulpd	xmm7,[ecx+0x10]	/* t20*s */					__asm	mulpd	xmm3,[ecx     ]	/* t28*c */
		__asm	subpd	xmm5,xmm6	/* xmm1 <-~t20*/				__asm	subpd	xmm1,xmm2	/* xmm5 <- it */
		__asm	addpd	xmm4,xmm7	/* xmm0 <-~t19*/				__asm	addpd	xmm0,xmm3	/* xmm4 <- rt */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy~t20*/
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy~t19*/

		__asm	addpd	xmm4,xmm0	/* ~t19 <- t19+rt */
		__asm	addpd	xmm5,xmm1	/* ~t20 <- t20+it */
		__asm	subpd	xmm6,xmm0	/* ~t27 <- t19-rt */
		__asm	subpd	xmm7,xmm1	/* ~t28 <- t20-it */

		__asm	movaps	xmm2,[eax+0x040]	/* t11/r13 */
		__asm	movaps	xmm3,[eax+0x050]	/* t12/r14 */
		__asm	movaps	xmm0,[eax      ]	/* t3 /r9  */
		__asm	movaps	xmm1,[eax+0x010]	/* t4 /r10 */
		__asm	addpd	xmm2,[eax+0x050]	/*~t11=t11+t12*/
		__asm	subpd	xmm3,[eax+0x040]	/*~t12=t12-t11*/
		__asm	mulpd	xmm2,[ebx]	/* rt */
		__asm	mulpd	xmm3,[ebx]	/* it */

		__asm	subpd	xmm0,xmm2	/*~t11 <- t3 - rt */
		__asm	subpd	xmm1,xmm3	/*~t12 <- t4 - it */
		__asm	addpd	xmm2,xmm2	/*          2* rt */
		__asm	addpd	xmm3,xmm3	/*          2* it */
		__asm	addpd	xmm2,xmm0	/*~t3  <- t3 + rt */
		__asm	addpd	xmm3,xmm1	/*~t4  <- t4 + it */

		__asm	mov	ebx, add1	/* restore main-array index [now into j2 block], but keep &r9 in eax for temporary storage */
		__asm	mov	ecx, c1
		__asm	mov	edx, c9
	/*
		rt       =t3 +t19;			it       =t4 +t20;
		t19      =t3 -t19;			t20      =t4 -t20;
		a[jt    ]=rt *c1 +it *s1;	a[jp    ]=it *c1 -rt *s1;
		a[jt+p8 ]=t19*c9 +t20*s9;	a[jp+p8 ]=t19*c9 -t20*s9;
	*/
		__asm	subpd	xmm2,xmm4		/*~t19 <- t3 -t19 */
		__asm	subpd	xmm3,xmm5		/*~t20 <- t4 -t20 */
		__asm	addpd	xmm4,xmm4		/*          2*t19 */
		__asm	addpd	xmm5,xmm5		/*          2*t20 */
		__asm	addpd	xmm4,xmm2		/* rt  <- t3 +t19 */
		__asm	addpd	xmm5,xmm3		/* it  <- t4 +t20 */
		__asm	movaps	[eax      ],xmm2	/* tmp store ~t3 */
		__asm	movaps	[eax+0x010],xmm3	/* tmp store ~t4 */
		__asm	movaps	xmm2,xmm4		/* rt copy */
		__asm	movaps	xmm3,xmm5		/* it copy */
		__asm	mulpd	xmm4,[ecx     ]	/* rt *c1 */
		__asm	mulpd	xmm5,[ecx     ]	/* it *c1 */
		__asm	mulpd	xmm2,[ecx+0x10]	/* rt *s1 */
		__asm	mulpd	xmm3,[ecx+0x10]	/* it *s1 */
		__asm	subpd	xmm5,xmm2	/* xmm5 <- im */
		__asm	addpd	xmm4,xmm3	/* xmm4 <- re, xmm2,3 free */

		__asm	movaps	[ebx+0x10],xmm5	/* a[jp    ] */
		__asm	movaps	[ebx     ],xmm4	/* a[jt    ] */

		__asm	movaps	xmm4,[eax      ]	/* load ~t3 */
		__asm	movaps	xmm5,[eax+0x010]	/* load ~t4 */
		__asm	movaps	xmm2,xmm4		/* re copy */
		__asm	movaps	xmm3,xmm5		/* im copy */
		__asm	mulpd	xmm4,[edx     ]	/* re *c9 */
		__asm	mulpd	xmm5,[edx     ]	/* im *c9 */
		__asm	mulpd	xmm2,[edx+0x10]	/* re *s9 */
		__asm	mulpd	xmm3,[edx+0x10]	/* im *s9 */
		__asm	subpd	xmm5,xmm2	/* xmm5 <- im */
		__asm	addpd	xmm4,xmm3	/* xmm4 <- re */

		__asm	movaps	[ebx+0x90],xmm5	/* a[jp+p8 ] */
		__asm	movaps	[ebx+0x80],xmm4	/* a[jt+p8 ] */

		__asm	mov	ecx, c5
		__asm	mov	edx, c13
	/*
		rt       =t11+t28;			it       =t12-t27;	// mpy by E^-4 = -I is inlined here...
		t28      =t11-t28;			t27      =t12+t27;
		a[jt+p4 ]=rt *c5 +it *s5;	a[jp+p4 ]=it *c5 -rt *s5;
		a[jt+p12]=t28*c13+t27*s13;	a[jp+p12]=t28*c13-t27*s13;
	*/
		__asm	subpd	xmm0,xmm7		/*~t28 <- t11-t28 */
		__asm	subpd	xmm1,xmm6		/* it  <- t12-t27 */
		__asm	addpd	xmm7,xmm7		/*          2*t28 */
		__asm	addpd	xmm6,xmm6		/*          2*t27 */
		__asm	addpd	xmm7,xmm0		/* rt  <- t11+t28 */
		__asm	addpd	xmm6,xmm1		/*~t27 <- t12+t27 */
		__asm	movaps	xmm4,xmm7		/* rt copy */
		__asm	movaps	xmm5,xmm1		/* it copy */
		__asm	mulpd	xmm7,[ecx     ]	/* rt*c5 */
		__asm	mulpd	xmm1,[ecx     ]	/* it*c5 */
		__asm	mulpd	xmm4,[ecx+0x10]	/* rt*s5 */
		__asm	mulpd	xmm5,[ecx+0x10]	/* it*s5 */
		__asm	subpd	xmm1,xmm4	/* xmm1 <- im */
		__asm	addpd	xmm7,xmm5	/* xmm7 <- re */

		__asm	movaps	[ebx+0x50],xmm1	/* a[jp+p4 ] */
		__asm	movaps	[ebx+0x40],xmm7	/* a[jt+p4 ] */

		__asm	movaps	xmm4,xmm0		/*t28 copy */
		__asm	movaps	xmm5,xmm6		/*t27 copy */
		__asm	mulpd	xmm0,[edx     ]	/*t28*c13 */
		__asm	mulpd	xmm6,[edx     ]	/*t27*c13 */
		__asm	mulpd	xmm4,[edx+0x10]	/*t28*s13 */
		__asm	mulpd	xmm5,[edx+0x10]	/*t27*s13 */
		__asm	subpd	xmm6,xmm4	/* xmm6 <- im */
		__asm	addpd	xmm0,xmm5	/* xmm7 <- re */

		__asm	movaps	[ebx+0xd0],xmm6	/* a[jp+p12] */
		__asm	movaps	[ebx+0xc0],xmm0	/* a[jt+p12] */

	/*...Block 4: t7,15,23,31 -> r25,29,27,31: */
		__asm	mov	eax, r25
		__asm	mov	ebx, isrt2
		__asm	mov	ecx, cc0

		__asm	movaps	xmm4,[eax+0x020]	/* t23/r27 */			__asm	movaps	xmm0,[eax+0x060]	/* t31/r31 */
		__asm	movaps	xmm5,[eax+0x030]	/* t24/r28 */			__asm	movaps	xmm1,[eax+0x070]	/* t32/r32 */
		__asm	movaps	xmm6,[eax+0x020]	/* xmm2 <- cpy t23 */	__asm	movaps	xmm2,[eax+0x060]	/* xmm6 <- cpy t31 */
		__asm	movaps	xmm7,[eax+0x030]	/* xmm3 <- cpy t24 */	__asm	movaps	xmm3,[eax+0x070]	/* xmm7 <- cpy t32 */

		__asm	mulpd	xmm4,[ecx+0x10]	/* t23*s */					__asm	mulpd	xmm0,[ecx     ]	/* t31*c */
		__asm	mulpd	xmm5,[ecx+0x10]	/* t24*s */					__asm	mulpd	xmm1,[ecx     ]	/* t32*c */
		__asm	mulpd	xmm6,[ecx     ]	/* t23*c */					__asm	mulpd	xmm2,[ecx+0x10]	/* t31*s */
		__asm	mulpd	xmm7,[ecx     ]	/* t24*c */					__asm	mulpd	xmm3,[ecx+0x10]	/* t32*s */
		__asm	subpd	xmm5,xmm6	/* xmm1 <-~t24*/				__asm	subpd	xmm1,xmm2	/* xmm5 <- it */
		__asm	addpd	xmm4,xmm7	/* xmm0 <-~t23*/				__asm	addpd	xmm0,xmm3	/* xmm4 <- rt */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy~t24*/
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy~t23*/

		__asm	addpd	xmm4,xmm0	/* ~t31 <- t23+rt */
		__asm	addpd	xmm5,xmm1	/* ~t32 <- t24+it */
		__asm	subpd	xmm6,xmm0	/* ~t23 <- t23-rt */
		__asm	subpd	xmm7,xmm1	/* ~t24 <- t24-it */

		__asm	movaps	xmm2,[eax+0x040]	/* t15/r29 */
		__asm	movaps	xmm3,[eax+0x050]	/* t16/r30 */
		__asm	movaps	xmm0,[eax      ]	/* t7 /r25 */
		__asm	movaps	xmm1,[eax+0x010]	/* t8 /r26 */
		__asm	subpd	xmm2,[eax+0x050]	/*~t15=t15-t16*/
		__asm	addpd	xmm3,[eax+0x040]	/*~t16=t16+t15*/
		__asm	mulpd	xmm2,[ebx]	/* rt */
		__asm	mulpd	xmm3,[ebx]	/* it */

		__asm	subpd	xmm0,xmm2	/*~t7  <- t7 - rt */
		__asm	subpd	xmm1,xmm3	/*~t8  <- t8 - it */
		__asm	addpd	xmm2,xmm2	/*          2* rt */
		__asm	addpd	xmm3,xmm3	/*          2* it */
		__asm	addpd	xmm2,xmm0	/*~t15 <- t7 + rt */
		__asm	addpd	xmm3,xmm1	/*~t16 <- t8 + it */

		__asm	mov	ebx, add1	/* restore main-array index [now into j2 block], but keep &r25 in eax for temporary storage */
		__asm	mov	ecx, c3
		__asm	mov	edx, c11
	/*
		rt       =t7 +t23;			it       =t8 +t24;
		t23      =t7 -t23;			t24      =t8 -t24;
		a[jt    ]=rt *c3 +it *s3;	a[jp    ]=it *c3 -rt *s3;
		a[jt+p8 ]=t23*c11+t24*s11;	a[jp+p8 ]=t23*c11-t24*s11;
	*/
		__asm	subpd	xmm0,xmm6		/*~t23 <- t7 -t23 */
		__asm	subpd	xmm1,xmm7		/*~t24 <- t8 -t24 */
		__asm	addpd	xmm6,xmm6		/*          2*t23 */
		__asm	addpd	xmm7,xmm7		/*          2*t24 */
		__asm	addpd	xmm6,xmm0		/* rt  <- t7 +t23 */
		__asm	addpd	xmm7,xmm1		/* it  <- t8 +t24 */
		__asm	movaps	[eax      ],xmm0	/* tmp store ~t23*/
		__asm	movaps	[eax+0x010],xmm1	/* tmp store ~t24*/
		__asm	movaps	xmm0,xmm6		/* rt copy */
		__asm	movaps	xmm1,xmm7		/* it copy */
		__asm	mulpd	xmm6,[ecx     ]	/* rt *c3 */
		__asm	mulpd	xmm7,[ecx     ]	/* it *c3 */
		__asm	mulpd	xmm0,[ecx+0x10]	/* rt *s3 */
		__asm	mulpd	xmm1,[ecx+0x10]	/* it *s3 */
		__asm	subpd	xmm7,xmm0	/* xmm7 <- im */
		__asm	addpd	xmm6,xmm1	/* xmm6 <- re, xmm0,1 free */

		__asm	movaps	[ebx+0x30],xmm7	/* a[jp    ] */
		__asm	movaps	[ebx+0x20],xmm6	/* a[jt    ] */

		__asm	movaps	xmm6,[eax      ]	/* load ~t23*/
		__asm	movaps	xmm7,[eax+0x010]	/* load ~t24*/
		__asm	movaps	xmm0,xmm6		/* re copy */
		__asm	movaps	xmm1,xmm7		/* im copy */
		__asm	mulpd	xmm6,[edx     ]	/* re *c11*/
		__asm	mulpd	xmm7,[edx     ]	/* im *c11*/
		__asm	mulpd	xmm0,[edx+0x10]	/* re *s11*/
		__asm	mulpd	xmm1,[edx+0x10]	/* im *s11*/
		__asm	subpd	xmm7,xmm0	/* xmm7 <- im */
		__asm	addpd	xmm6,xmm1	/* xmm6 <- re */

		__asm	movaps	[ebx+0xb0],xmm7	/* a[jp+p8 ] */
		__asm	movaps	[ebx+0xa0],xmm6	/* a[jt+p8 ] */

		__asm	mov	ecx, c7
		__asm	mov	edx, c15
	/*
		rt		 =t15+t32;			it		 =t16-t31;	// mpy by E^-4 = -I is inlined here...
		t32		 =t15-t32;			t31		 =t16+t31;
		a[jt+p4 ]=rt *c7 +it *s7;	a[jp+p4 ]=it *c7 -rt *s7;
		a[jt+p12]=t32*c15+t31*s15;	a[jp+p12]=t32*c15-t31*s15;
	*/
		__asm	subpd	xmm2,xmm5		/*~t32 <- t15-t32 */
		__asm	subpd	xmm3,xmm4		/* it  <- t16-t31 */
		__asm	addpd	xmm5,xmm5		/*          2*t32 */
		__asm	addpd	xmm4,xmm4		/*          2*t31 */
		__asm	addpd	xmm5,xmm2		/* rt  <- t15+t32 */
		__asm	addpd	xmm4,xmm3		/*~t31 <- t16+t31 */
		__asm	movaps	xmm0,xmm5		/* rt copy */
		__asm	movaps	xmm1,xmm3		/* it copy */
		__asm	mulpd	xmm5,[ecx     ]	/* rt*c7 */
		__asm	mulpd	xmm3,[ecx     ]	/* it*c7 */
		__asm	mulpd	xmm0,[ecx+0x10]	/* rt*s7 */
		__asm	mulpd	xmm1,[ecx+0x10]	/* it*s7 */
		__asm	subpd	xmm3,xmm0	/* xmm3 <- im */
		__asm	addpd	xmm5,xmm1	/* xmm5 <- re, xmm6,7 free */

		__asm	movaps	[ebx+0x70],xmm3	/* a[jp+p4 ] */
		__asm	movaps	[ebx+0x60],xmm5	/* a[jt+p4 ] */

		__asm	movaps	xmm0,xmm2		/*t32 copy */
		__asm	movaps	xmm1,xmm4		/*t31 copy */
		__asm	mulpd	xmm2,[edx     ]	/*t32*c15 */
		__asm	mulpd	xmm4,[edx     ]	/*t31*c15 */
		__asm	mulpd	xmm0,[edx+0x10]	/*t32*s15 */
		__asm	mulpd	xmm1,[edx+0x10]	/*t31*s15 */
		__asm	subpd	xmm4,xmm0	/* xmm4 <- im */
		__asm	addpd	xmm2,xmm1	/* xmm2 <- re */

		__asm	movaps	[ebx+0xf0],xmm4	/* a[jp+p12] */
		__asm	movaps	[ebx+0xe0],xmm2	/* a[jt+p12] */

	/*...Block 1: t1,9,17,25 -> r1,5,3,7: */
		__asm	mov	eax, r1

		__asm	movaps	xmm0,[eax      ]	/* t1 /r1 */
		__asm	movaps	xmm1,[eax+0x010]	/* t2 /r2 */
		__asm	movaps	xmm2,[eax+0x040]	/* t9 /r5 */
		__asm	movaps	xmm3,[eax+0x050]	/* t10/r6 */

		__asm	subpd	xmm0,[eax+0x040]	/*~t9 =t1 -t9 */
		__asm	subpd	xmm1,[eax+0x050]	/*~t10=t2 -t10*/
		__asm	addpd	xmm2,[eax      ]	/*~t1 =t9 +t1 */
		__asm	addpd	xmm3,[eax+0x010]	/*~t2 =t10+t2 */

		__asm	movaps	xmm4,[eax+0x020]	/* t17/r3 */
		__asm	movaps	xmm5,[eax+0x030]	/* t18/r4 */
		__asm	movaps	xmm6,[eax+0x060]	/* t25/r7 */
		__asm	movaps	xmm7,[eax+0x070]	/* t26/r8 */

		__asm	subpd	xmm4,[eax+0x060]	/*~t25=t17-t25*/
		__asm	subpd	xmm5,[eax+0x070]	/*~t26=t18-t26*/
		__asm	addpd	xmm6,[eax+0x020]	/*~t17=t25+t17*/
		__asm	addpd	xmm7,[eax+0x030]	/*~t18=t26+t18*/

		__asm	mov	eax, add0	/* restore main-array index */
		__asm	mov	edx, c8
	/*
		a[jt    ]=t1+t17;			a[jp    ]=t2+t18;
		t17      =t1-t17;			t18      =t2-t18;
		a[jt+p8 ]=t17*c8 +t18*s8;	a[jp+p8 ]=t18*c8 -t17*s8;
	*/
		__asm	addpd	xmm2,xmm6		/* t1 +t17 */
		__asm	addpd	xmm3,xmm7		/* t2 +t18 */

		__asm	movaps	[eax     ],xmm2	/* tmp store non-interleaved a[jt+p0 ] */
		__asm	movaps	[eax+0x10],xmm3	/* tmp store non-interleaved a[jp+p0 ] */
		__asm	addpd	xmm6,xmm6		/*   2*t17 */
		__asm	addpd	xmm7,xmm7		/*   2*t18 */
		__asm	subpd	xmm2,xmm6		/*~t17 <- t1 -t17 */
		__asm	subpd	xmm3,xmm7		/*~t18 <- t2 -t18 */
		__asm	movaps	xmm6,xmm2		/*~t17 copy */
		__asm	movaps	xmm7,xmm3		/*~t18 copy */
		__asm	mulpd	xmm2,[edx     ]	/*~t17*c8 */
		__asm	mulpd	xmm3,[edx     ]	/*~t18*c8 */
		__asm	mulpd	xmm6,[edx+0x10]	/*~t17*s8 */
		__asm	mulpd	xmm7,[edx+0x10]	/*~t18*s8 */
		__asm	subpd	xmm3,xmm6	/* xmm3 <- it */
		__asm	addpd	xmm2,xmm7	/* xmm2 <- rt 	xmm6,7 free */

		__asm	mov	ecx, add1	/* Need address of j2-block for butterfly-swap-on-write */
		__asm	movaps		xmm7,xmm3	/* cpy a[jp    ] */
		__asm	movaps		xmm6,xmm2	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm7,[ecx+0x90]
		__asm	unpcklpd	xmm3,[ecx+0x90]
		__asm	movaps	[ecx+0x90],xmm7	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm6,[ecx+0x80]
		__asm	unpcklpd	xmm2,[ecx+0x80]
		__asm	movaps	[ecx+0x80],xmm6	/* Store hi real in aj2 */

		__asm	movaps	[eax+0x90],xmm3	/* a[jp+p8 ] */
		__asm	movaps	[eax+0x80],xmm2	/* a[jt+p8 ] */

		__asm	movaps	xmm3,[eax+0x10]	/* reload a[jp+p0 ] */
		__asm	movaps	xmm2,[eax     ]	/* reload a[jt+p0 ] */
		__asm	movaps		xmm7,xmm3	/* cpy a[jp    ] */
		__asm	movaps		xmm6,xmm2	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm7,[ecx+0x10]
		__asm	unpcklpd	xmm3,[ecx+0x10]
		__asm	movaps	[ecx+0x10],xmm7	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm6,[ecx     ]
		__asm	unpcklpd	xmm2,[ecx     ]
		__asm	movaps	[ecx     ],xmm6	/* Store hi real in aj2 */

		__asm	movaps	[eax+0x10],xmm3	/* a[jp+p0 ] */
		__asm	movaps	[eax     ],xmm2	/* a[jt+p0 ] */

		__asm	mov	ecx, c4
		__asm	mov	edx, c12
	/* mpy by E^-4 = -I is inlined here...
		rt       =t9 +t26;			it       =t10-t25;
		t26      =t9 -t26;			t25      =t10+t25;
		a[jt+p4 ]=rt *c4 +it *s4;	a[jp+p4 ]=it *c4 -rt *s4;
		a[jt+p12]=t26*c12+t25*s12;	a[jp+p12]=t25*c12-t26*s12;
	*/
		__asm	addpd	xmm0,xmm5		/* rt  <- t9 +t26 */
		__asm	subpd	xmm1,xmm4		/* it  <- t10-t25 */
		__asm	movaps	xmm2,xmm0		/* rt  copy */
		__asm	movaps	xmm3,xmm1		/* it  copy */
		__asm	addpd	xmm5,xmm5		/*          2*t26 */
		__asm	addpd	xmm4,xmm4		/*          2*t25 */
		__asm	movaps	xmm6,xmm0		/* rt  copy */
		__asm	movaps	xmm7,xmm1		/* it  copy */
		__asm	mulpd	xmm2,[ecx     ]	/* rt *c4 */
		__asm	mulpd	xmm3,[ecx     ]	/* it *c4 */
		__asm	mulpd	xmm6,[ecx+0x10]	/* rt *s4 */
		__asm	mulpd	xmm7,[ecx+0x10]	/* it *s4 */
		__asm	subpd	xmm3,xmm6	/* xmm3 <- im */
		__asm	addpd	xmm2,xmm7	/* xmm2 <- re 	xmm6,7 free */

		__asm	mov	ecx, add1	/* Need address of j2-block for butterfly-swap-on-write */
		__asm	movaps		xmm7,xmm3	/* cpy a[jp    ] */
		__asm	movaps		xmm6,xmm2	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm7,[ecx+0x50]
		__asm	unpcklpd	xmm3,[ecx+0x50]
		__asm	movaps	[ecx+0x50],xmm7	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm6,[ecx+0x40]
		__asm	unpcklpd	xmm2,[ecx+0x40]
		__asm	movaps	[ecx+0x40],xmm6	/* Store hi real in aj2 */

		__asm	movaps	[eax+0x50],xmm3	/* a[jp+p4 ] */
		__asm	movaps	[eax+0x40],xmm2	/* a[jt+p4 ] */

		__asm	subpd	xmm0,xmm5		/*~t26 <- t9 -t26 */
		__asm	addpd	xmm1,xmm4		/*~t25 <- t10+t25 */
		__asm	movaps	xmm6,xmm0		/*~t26 copy */
		__asm	movaps	xmm7,xmm1		/*~t25 copy */
		__asm	mulpd	xmm0,[edx     ]	/*~t26*c12*/
		__asm	mulpd	xmm1,[edx     ]	/*~t25*c12*/
		__asm	mulpd	xmm6,[edx+0x10]	/*~t26*s12*/
		__asm	mulpd	xmm7,[edx+0x10]	/*~t25*s12*/
		__asm	subpd	xmm1,xmm6	/* xmm3 <- im */
		__asm	addpd	xmm0,xmm7	/* xmm2 <- re 	xmm6,7 free */

		__asm	mov	ecx, add1	/* Need address of j2-block for butterfly-swap-on-write */
		__asm	movaps		xmm7,xmm1	/* cpy a[jp    ] */
		__asm	movaps		xmm6,xmm0	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm7,[ecx+0xd0]
		__asm	unpcklpd	xmm1,[ecx+0xd0]
		__asm	movaps	[ecx+0xd0],xmm7	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm6,[ecx+0xc0]
		__asm	unpcklpd	xmm0,[ecx+0xc0]
		__asm	movaps	[ecx+0xc0],xmm6	/* Store hi real in aj2 */

		__asm	movaps	[eax+0xd0],xmm1	/* a[jp+p12] */
		__asm	movaps	[eax+0xc0],xmm0	/* a[jt+p12] */

	/*...Block 3: t5,13,21,29 -> r17,21,19,23: */
		__asm	mov	eax, r17
		__asm	mov	ebx, isrt2
		__asm	movaps	xmm2,[ebx]	/* isrt2 */

		__asm	movaps	xmm4,[eax+0x020]	/* t21/r19 */
		__asm	movaps	xmm5,[eax+0x030]	/* t22/r20 */
		__asm	movaps	xmm0,[eax+0x060]	/* t29/r23 */
		__asm	movaps	xmm1,[eax+0x070]	/* t30/r24 */

		__asm	addpd	xmm4,[eax+0x030]	/*~t21=t21+t22*/
		__asm	subpd	xmm5,[eax+0x020]	/*~t22=t22-t21*/
		__asm	subpd	xmm0,[eax+0x070]	/* rt =t29-t30*/
		__asm	addpd	xmm1,[eax+0x060]	/* it =t30+t29*/
		__asm	mulpd	xmm4,xmm2
		__asm	mulpd	xmm5,xmm2
		__asm	mulpd	xmm0,xmm2
		__asm	mulpd	xmm1,xmm2
		__asm	movaps	xmm6,xmm4			/* t21 copy */
		__asm	movaps	xmm7,xmm5			/* t22 copy */

		__asm	subpd	xmm4,xmm0			/*~t21=t21-rt */
		__asm	subpd	xmm5,xmm1			/*~t22=t22-it */
		__asm	addpd	xmm6,xmm0			/*~t29=t21+rt */
		__asm	addpd	xmm7,xmm1			/*~t30=t22+it */

		__asm	movaps	xmm0,[eax      ]	/* t5 /r17 */
		__asm	movaps	xmm1,[eax+0x010]	/* t6 /r18 */
		__asm	movaps	xmm2,[eax+0x040]	/* t13/r21 */
		__asm	movaps	xmm3,[eax+0x050]	/* t14/r22 */

		__asm	subpd	xmm0,[eax+0x050]	/*~t13=t5 -t14*/
		__asm	subpd	xmm1,[eax+0x040]	/*~t6 =t6 -t13*/
		__asm	addpd	xmm3,[eax      ]	/*~t5 =t14+t5 */
		__asm	addpd	xmm2,[eax+0x010]	/*~t14=t13+t6 */

		__asm	mov	ebx, add0	/* restore main-array index, but keep &r17 in eax for temporary storage */
		__asm	mov	ecx, c2
		__asm	mov	edx, c10
	/*
		rt   =t5 +t21;				it   =t6 +t22;
		t21  =t5 -t21;				t22  =t6 -t22;
		a[jt    ]=rt *c2 +it *s2;	a[jp    ]=it *c2 -rt *s2;
		a[jt+p8 ]=t21*c10+t22*s10;	a[jp+p8 ]=t21*c10-t22*s10;
	*/
		__asm	subpd	xmm3,xmm4		/*~t21 <- t5 -t21 */
		__asm	subpd	xmm1,xmm5		/*~t22 <- t6 -t22 */
		__asm	addpd	xmm4,xmm4		/*          2*t21 */
		__asm	addpd	xmm5,xmm5		/*          2*t22 */
		__asm	addpd	xmm4,xmm3		/* rt  <- t5 +t21 */
		__asm	addpd	xmm5,xmm1		/* it  <- t6 +t22 */
		__asm	movaps	[eax      ],xmm3	/* tmp store ~t21*/
		__asm	movaps	[eax+0x010],xmm1	/* tmp store ~t22*/
		__asm	movaps	xmm3,xmm4		/* rt copy */
		__asm	movaps	xmm1,xmm5		/* it copy */
		__asm	mulpd	xmm4,[ecx     ]	/* rt*c2 */
		__asm	mulpd	xmm5,[ecx     ]	/* it*c2 */
		__asm	mulpd	xmm3,[ecx+0x10]	/* rt*s2 */
		__asm	mulpd	xmm1,[ecx+0x10]	/* it*s2 */
		__asm	subpd	xmm5,xmm3	/* xmm5 <- im */
		__asm	addpd	xmm4,xmm1	/* xmm4 <- re 	xmm1,3 free */

		__asm	mov	ecx, add1	/* Need address of j2-block for butterfly-swap-on-write */
		__asm	movaps		xmm3,xmm5	/* cpy a[jp    ] */
		__asm	movaps		xmm1,xmm4	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm3,[ecx+0x30]
		__asm	unpcklpd	xmm5,[ecx+0x30]
		__asm	movaps	[ecx+0x30],xmm3	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm1,[ecx+0x20]
		__asm	unpcklpd	xmm4,[ecx+0x20]
		__asm	movaps	[ecx+0x20],xmm1	/* Store hi real in aj2 */

		__asm	movaps	[ebx+0x30],xmm5	/* a[jp    ] */
		__asm	movaps	[ebx+0x20],xmm4	/* a[jt    ] */

		__asm	movaps	xmm4,[eax      ]	/* load ~t5 */
		__asm	movaps	xmm5,[eax+0x010]	/* load ~t6 */
		__asm	movaps	xmm3,xmm4		/* re copy */
		__asm	movaps	xmm1,xmm5		/* im copy */
		__asm	mulpd	xmm4,[edx     ]	/* re*c10 */
		__asm	mulpd	xmm5,[edx     ]	/* im*c10 */
		__asm	mulpd	xmm3,[edx+0x10]	/* re*s10 */
		__asm	mulpd	xmm1,[edx+0x10]	/* im*s10 */
		__asm	subpd	xmm5,xmm3	/* xmm5 <- im */
		__asm	addpd	xmm4,xmm1	/* xmm4 <- re 	xmm1,3 free */

		__asm	movaps		xmm3,xmm5	/* cpy a[jp    ] */
		__asm	movaps		xmm1,xmm4	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm3,[ecx+0xb0]
		__asm	unpcklpd	xmm5,[ecx+0xb0]
		__asm	movaps	[ecx+0xb0],xmm3	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm1,[ecx+0xa0]
		__asm	unpcklpd	xmm4,[ecx+0xa0]
		__asm	movaps	[ecx+0xa0],xmm1	/* Store hi real in aj2 */

		__asm	movaps	[ebx+0xb0],xmm5	/* a[jp+p8 ] */
		__asm	movaps	[ebx+0xa0],xmm4	/* a[jt+p8 ] */

		__asm	mov	ecx, c6
		__asm	mov	edx, c14
	/*
		rt		 =t13+t30;			it		 =t14-t29;	// mpy by E^-4 = -I is inlined here...
		t30		 =t13-t30;			t29		 =t14+t29;
		a[jt+p4 ]=rt *c6 +it *s6;	a[jp+p4 ]=it *c6 -rt *s6;
		a[jt+p12]=t30*c14+t29*s14;	a[jp+p12]=t30*c14-t29*s14;
	*/
		__asm	subpd	xmm0,xmm7		/*~t30 <- t13-t30 */
		__asm	subpd	xmm2,xmm6		/* it  <- t14-t29 */
		__asm	addpd	xmm7,xmm7		/*          2*t30 */
		__asm	addpd	xmm6,xmm6		/*          2*t29 */
		__asm	addpd	xmm7,xmm0		/* rt  <- t13+t30 */
		__asm	addpd	xmm6,xmm2		/*~t29 <- t14+t29 */
		__asm	movaps	xmm4,xmm7		/* rt copy */
		__asm	movaps	xmm5,xmm2		/* it copy */
		__asm	mulpd	xmm7,[ecx     ]	/* rt*c6 */
		__asm	mulpd	xmm2,[ecx     ]	/* it*c6 */
		__asm	mulpd	xmm4,[ecx+0x10]	/* rt*s6 */
		__asm	mulpd	xmm5,[ecx+0x10]	/* it*s6 */
		__asm	subpd	xmm2,xmm4	/* xmm2 <- im */
		__asm	addpd	xmm7,xmm5	/* xmm7 <- re 	xmm4,5 free */

		__asm	mov	ecx, add1	/* Need address of j2-block for butterfly-swap-on-write */
		__asm	movaps		xmm5,xmm2	/* cpy a[jp    ] */
		__asm	movaps		xmm4,xmm7	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm5,[ecx+0x70]
		__asm	unpcklpd	xmm2,[ecx+0x70]
		__asm	movaps	[ecx+0x70],xmm5	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm4,[ecx+0x60]
		__asm	unpcklpd	xmm7,[ecx+0x60]
		__asm	movaps	[ecx+0x60],xmm4	/* Store hi real in aj2 */

		__asm	movaps	[ebx+0x70],xmm2	/* a[jp+p4 ] */
		__asm	movaps	[ebx+0x60],xmm7	/* a[jt+p4 ] */

		__asm	movaps	xmm4,xmm0		/* t30 copy */
		__asm	movaps	xmm5,xmm6		/* t29 copy */
		__asm	mulpd	xmm0,[edx     ]	/* t30*c14 */
		__asm	mulpd	xmm6,[edx     ]	/* t29*c14 */
		__asm	mulpd	xmm4,[edx+0x10]	/* t30*s14 */
		__asm	mulpd	xmm5,[edx+0x10]	/* t29*s14 */
		__asm	subpd	xmm6,xmm4	/* xmm6 <- im */
		__asm	addpd	xmm0,xmm5	/* xmm7 <- re 	xmm4,5 free */

		__asm	movaps		xmm5,xmm6	/* cpy a[jp    ] */
		__asm	movaps		xmm4,xmm0	/* cpy a[jt    ] */
		__asm	unpckhpd	xmm5,[ecx+0xf0]
		__asm	unpcklpd	xmm6,[ecx+0xf0]
		__asm	movaps	[ecx+0xf0],xmm5	/* Store hi imag in aj2 */
		__asm	unpckhpd	xmm4,[ecx+0xe0]
		__asm	unpcklpd	xmm0,[ecx+0xe0]
		__asm	movaps	[ecx+0xe0],xmm4	/* Store hi real in aj2 */

		__asm	movaps	[ebx+0xf0],xmm6	/* a[jp+p12] */
		__asm	movaps	[ebx+0xe0],xmm0	/* a[jt+p12] */

	#else	/* GCC-style inline ASM: */

		SSE2_RADIX16_WRAPPER_DIT(add0,add1,r1,r9,r17,r25,isrt2,cc0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15)

	#endif

#else	/* if(!USE_SSE2) */

	/*...Block 1: */
		t1 =a[j1   ];					t2 =a[j2   ];
		rt =a[j1+16]*c8 -a[j2+16]*s8 ;	it =a[j2+16]*c8 +a[j1+16]*s8;
		t3 =t1 -rt;	t1 =t1 +rt;
		t4 =t2 -it;	t2 =t2 +it;

		t5 =a[j1+8 ]*c4 -a[j2+8 ]*s4;	t6 =a[j2+8 ]*c4 +a[j1+8 ]*s4;
		rt =a[j1+24]*c12-a[j2+24]*s12;	it =a[j2+24]*c12+a[j1+24]*s12;
		t7 =t5 -rt;	t5 =t5 +rt;
		t8 =t6 -it;	t6 =t6 +it;

		rt =t5;	t5 =t1 -rt;		t1 =t1 +rt;
		it =t6;	t6 =t2 -it;		t2 =t2 +it;

		rt =t7;	t7 =t3 +t8;		t3 =t3 -t8;
				t8 =t4 -rt;		t4 =t4 +rt;

	/*...Block 2: */
		t9 =a[j1+4 ]*c2 -a[j2+4 ]*s2 ;	t10=a[j2+4 ]*c2 +a[j1+4 ]*s2 ;
		rt =a[j1+20]*c10-a[j2+20]*s10;	it =a[j2+20]*c10+a[j1+20]*s10;
		t11=t9 -rt;	t9 =t9 +rt;
		t12=t10-it;	t10=t10+it;

		t13=a[j1+12]*c6 -a[j2+12]*s6 ;	t14=a[j2+12]*c6 +a[j1+12]*s6 ;
		rt =a[j1+28]*c14-a[j2+28]*s14;	it =a[j2+28]*c14+a[j1+28]*s14;
		t15=t13-rt;	t13=t13+rt;
		t16=t14-it;	t14=t14+it;

		rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
		it =t14;	t14=t10-it;	t10=t10+it;

		rt =t15;	t15=t11+t16;t11=t11-t16;
					t16=t12-rt;	t12=t12+rt;

	/*...Block 3: */
		t17=a[j1+2 ]*c1 -a[j2+2 ]*s1 ;	t18=a[j2+2 ]*c1 +a[j1+2 ]*s1 ;
		rt =a[j1+18]*c9 -a[j2+18]*s9 ;	it =a[j2+18]*c9 +a[j1+18]*s9 ;
		t19=t17-rt;	t17=t17+rt;
		t20=t18-it;	t18=t18+it;

		t21=a[j1+10]*c5 -a[j2+10]*s5 ;	t22=a[j2+10]*c5 +a[j1+10]*s5 ;
		rt =a[j1+26]*c13-a[j2+26]*s13;	it =a[j2+26]*c13+a[j1+26]*s13;
		t23=t21-rt;	t21=t21+rt;
		t24=t22-it;	t22=t22+it;

		rt =t21;	t21=t17-rt;	t17=t17+rt;
		it =t22;	t22=t18-it;	t18=t18+it;

		rt =t23;	t23=t19+t24;t19=t19-t24;
					t24=t20-rt;	t20=t20+rt;

	/*...Block 4: */
		t25=a[j1+6 ]*c3 -a[j2+6 ]*s3 ;	t26=a[j2+6 ]*c3 +a[j1+6 ]*s3 ;
		rt =a[j1+22]*c11-a[j2+22]*s11;	it =a[j2+22]*c11+a[j1+22]*s11;
		t27=t25-rt;	t25=t25+rt;
		t28=t26-it;	t26=t26+it;

		t29=a[j1+14]*c7 -a[j2+14]*s7 ;	t30=a[j2+14]*c7 +a[j1+14]*s7 ;
		rt =a[j1+30]*c15-a[j2+30]*s15;	it =a[j2+30]*c15+a[j1+30]*s15;
		t31=t29-rt;	t29=t29+rt;
		t32=t30-it;	t30=t30+it;

		rt =t29;	t29=t25-rt;	t25=t25+rt;
		it =t30;	t30=t26-it;	t26=t26+it;

		rt =t31;	t31=t27+t32;t27=t27-t32;
					t32=t28-rt;	t28=t28+rt;

	#if 0//DEBUG_SSE2
		fprintf(stderr, "Outputs A:\n");
		fprintf(stderr, "radix16_dyadic_square: t1:= %20.5f, %20.5f\n",t1 ,t2 );
		fprintf(stderr, "radix16_dyadic_square: t3:= %20.5f, %20.5f\n",t3 ,t4 );
		fprintf(stderr, "radix16_dyadic_square: t5:= %20.5f, %20.5f\n",t5 ,t6 );
		fprintf(stderr, "radix16_dyadic_square: t7:= %20.5f, %20.5f\n",t7 ,t8 );
		fprintf(stderr, "radix16_dyadic_square: t9:= %20.5f, %20.5f\n",t9 ,t10);
		fprintf(stderr, "radix16_dyadic_square: t11= %20.5f, %20.5f\n",t11,t12);
		fprintf(stderr, "radix16_dyadic_square: t13= %20.5f, %20.5f\n",t13,t14);
		fprintf(stderr, "radix16_dyadic_square: t15= %20.5f, %20.5f\n",t15,t16);
		fprintf(stderr, "radix16_dyadic_square: t17= %20.5f, %20.5f\n",t17,t18);
		fprintf(stderr, "radix16_dyadic_square: t19= %20.5f, %20.5f\n",t19,t20);
		fprintf(stderr, "radix16_dyadic_square: t21= %20.5f, %20.5f\n",t21,t22);
		fprintf(stderr, "radix16_dyadic_square: t23= %20.5f, %20.5f\n",t23,t24);
		fprintf(stderr, "radix16_dyadic_square: t25= %20.5f, %20.5f\n",t25,t26);
		fprintf(stderr, "radix16_dyadic_square: t27= %20.5f, %20.5f\n",t27,t28);
		fprintf(stderr, "radix16_dyadic_square: t29= %20.5f, %20.5f\n",t29,t30);
		fprintf(stderr, "radix16_dyadic_square: t31= %20.5f, %20.5f\n",t31,t32);
		fprintf(stderr, "\n");
		exit(0);
	#endif

	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	!	1, exp(i* 1*twopi/16) =       ( c, s), exp(i* 2*twopi/16) = ISRT2*( 1, 1), exp(i* 3*twopi/16) =       ( s, c) (for inputs to transform block 2)
	!	1, exp(i* 2*twopi/16) = ISRT2*( 1, 1), exp(i* 4*twopi/16) =       ( 0, 1), exp(i* 6*twopi/16) = ISRT2*(-1, 1) (for inputs to transform block 3)
	!	1, exp(i* 3*twopi/16) =       ( s, c), exp(i* 6*twopi/16) = ISRT2*(-1, 1), exp(i* 9*twopi/16) =       (-c,-s) (for inputs to transform block 4).
	!  (This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
	!   I.e. do similar as above, except inputs a[j1  +p0:15:1] are replaced by t0:30:2,
	!						 a[j2+p0:15:1] are replaced by t1:31:2, and v.v. for outputs,
	!   and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.
	*/
	/*...Block 1: t1,9,17,25 */
		rt =t9;		t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		aj1p0r =t1 +t17;	aj1p0i =t2 +t18;
		aj1p1r =t1 -t17;	aj1p1i =t2 -t18;

		aj1p2r =t9 -t26;	aj1p2i =t10+t25;
		aj1p3r =t9 +t26;	aj1p3i =t10-t25;
	/*...Block 3: t5,13,21,29 */
		rt =t13;	t13=t5 +t14;t5 =t5 -t14;
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

	#if 0//DEBUG_SSE2
		fprintf(stderr, "Outputs A:\n");
		fprintf(stderr, "radix16_dyadic_square: aj1p0 = %20.5f, %20.5f\n",aj1p0r ,aj1p0i );
		fprintf(stderr, "radix16_dyadic_square: aj1p1 = %20.5f, %20.5f\n",aj1p1r ,aj1p1i );
		fprintf(stderr, "radix16_dyadic_square: aj1p2 = %20.5f, %20.5f\n",aj1p2r ,aj1p2i );
		fprintf(stderr, "radix16_dyadic_square: aj1p3 = %20.5f, %20.5f\n",aj1p3r ,aj1p3i );
		fprintf(stderr, "radix16_dyadic_square: aj1p4 = %20.5f, %20.5f\n",aj1p4r ,aj1p4i );
		fprintf(stderr, "radix16_dyadic_square: aj1p5 = %20.5f, %20.5f\n",aj1p5r ,aj1p5i );
		fprintf(stderr, "radix16_dyadic_square: aj1p6 = %20.5f, %20.5f\n",aj1p6r ,aj1p6i );
		fprintf(stderr, "radix16_dyadic_square: aj1p7 = %20.5f, %20.5f\n",aj1p7r ,aj1p7i );
		fprintf(stderr, "radix16_dyadic_square: aj1p8 = %20.5f, %20.5f\n",aj1p8r ,aj1p8i );
		fprintf(stderr, "radix16_dyadic_square: aj1p9 = %20.5f, %20.5f\n",aj1p9r ,aj1p9i );
		fprintf(stderr, "radix16_dyadic_square: aj1p10= %20.5f, %20.5f\n",aj1p10r,aj1p10i);
		fprintf(stderr, "radix16_dyadic_square: aj1p11= %20.5f, %20.5f\n",aj1p11r,aj1p11i);
		fprintf(stderr, "radix16_dyadic_square: aj1p12= %20.5f, %20.5f\n",aj1p12r,aj1p12i);
		fprintf(stderr, "radix16_dyadic_square: aj1p13= %20.5f, %20.5f\n",aj1p13r,aj1p13i);
		fprintf(stderr, "radix16_dyadic_square: aj1p14= %20.5f, %20.5f\n",aj1p14r,aj1p14i);
		fprintf(stderr, "radix16_dyadic_square: aj1p15= %20.5f, %20.5f\n",aj1p15r,aj1p15i);
		exit(0);
	#endif

	/*...Dyadic square of the forward FFT outputs: */
	#if !FFT_DEBUG
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
	/*...And do an inverse DIT radix-16 pass on the squared-data blocks. */

	#if 0// DEBUG_SSE2
		fprintf(stderr, "Outputs A^2:\n");
		fprintf(stderr, "radix16_dyadic_square: aj1p0 = %20.5f, %20.5f\n",aj1p0r ,aj1p0i );
		fprintf(stderr, "radix16_dyadic_square: aj1p1 = %20.5f, %20.5f\n",aj1p1r ,aj1p1i );
		fprintf(stderr, "radix16_dyadic_square: aj1p2 = %20.5f, %20.5f\n",aj1p2r ,aj1p2i );
		fprintf(stderr, "radix16_dyadic_square: aj1p3 = %20.5f, %20.5f\n",aj1p3r ,aj1p3i );
		fprintf(stderr, "radix16_dyadic_square: aj1p4 = %20.5f, %20.5f\n",aj1p4r ,aj1p4i );
		fprintf(stderr, "radix16_dyadic_square: aj1p5 = %20.5f, %20.5f\n",aj1p5r ,aj1p5i );
		fprintf(stderr, "radix16_dyadic_square: aj1p6 = %20.5f, %20.5f\n",aj1p6r ,aj1p6i );
		fprintf(stderr, "radix16_dyadic_square: aj1p7 = %20.5f, %20.5f\n",aj1p7r ,aj1p7i );
		fprintf(stderr, "radix16_dyadic_square: aj1p8 = %20.5f, %20.5f\n",aj1p8r ,aj1p8i );
		fprintf(stderr, "radix16_dyadic_square: aj1p9 = %20.5f, %20.5f\n",aj1p9r ,aj1p9i );
		fprintf(stderr, "radix16_dyadic_square: aj1p10= %20.5f, %20.5f\n",aj1p10r,aj1p10i);
		fprintf(stderr, "radix16_dyadic_square: aj1p11= %20.5f, %20.5f\n",aj1p11r,aj1p11i);
		fprintf(stderr, "radix16_dyadic_square: aj1p12= %20.5f, %20.5f\n",aj1p12r,aj1p12i);
		fprintf(stderr, "radix16_dyadic_square: aj1p13= %20.5f, %20.5f\n",aj1p13r,aj1p13i);
		fprintf(stderr, "radix16_dyadic_square: aj1p14= %20.5f, %20.5f\n",aj1p14r,aj1p14i);
		fprintf(stderr, "radix16_dyadic_square: aj1p15= %20.5f, %20.5f\n",aj1p15r,aj1p15i);
		exit(0);
	#endif

	/* 1st set of inputs: */
	#if PFETCH
	add0 = &a[j1+32];
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
	!						 a[j2+p0:15:1] are replaced by t1:31:2, and v.v. for outputs,
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

		a[j1   ]=t1+t17;				a[j2   ]=t2+t18;
		t1	=t1-t17;					t2	=t2-t18;
		a[j1+16]=t1 *c8 +t2 *s8 ;	a[j2+16]=t2 *c8 -t1 *s8;

		rt	=t9 +t26;					it	=t10-t25;
		t9	=t9 -t26;					t10	=t10+t25;
		a[j1+8 ]=rt *c4 +it *s4 ;	a[j2+8 ]=it *c4 -rt *s4;
		a[j1+24]=t9 *c12+t10*s12;	a[j2+24]=t10*c12-t9 *s12;

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
		a[j1+4 ]=rt *c2 +it *s2 ;	a[j2+4 ]=it *c2 -rt *s2 ;
		a[j1+20]=t5 *c10+t6 *s10;	a[j2+20]=t6 *c10-t5 *s10;

		rt	=t13+t30;					it	=t14-t29;
		t13	=t13-t30;					t14	=t14+t29;
		a[j1+12]=rt *c6 +it *s6 ;	a[j2+12]=it *c6 -rt *s6 ;
		a[j1+28]=t13*c14+t14*s14;	a[j2+28]=t14*c14-t13*s14;

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
		a[j1+2 ]=rt *c1 +it *s1 ;	a[j2+2 ]=it *c1 -rt *s1 ;
		a[j1+18]=t3 *c9 +t4 *s9 ;	a[j2+18]=t4 *c9 -t3 *s9 ;

		rt	=t11+t28;					it	=t12-t27;
		t11	=t11-t28;					t12	=t12+t27;
		a[j1+10]=rt *c5 +it *s5 ;	a[j2+10]=it *c5 -rt *s5 ;
		a[j1+26]=t11*c13+t12*s13;	a[j2+26]=t12*c13-t11*s13;

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
		a[j1+6 ]=rt *c3 +it *s3 ;	a[j2+6 ]=it *c3 -rt *s3 ;
		a[j1+22]=t7 *c11+t8 *s11;	a[j2+22]=t8 *c11-t7 *s11;

		rt	=t15+t32;					it	=t16-t31;
		t15	=t15-t32;					t16	=t16+t31;
		a[j1+14]=rt *c7 +it *s7 ;	a[j2+14]=it *c7 -rt *s7 ;
		a[j1+30]=t15*c15+t16*s15;	a[j2+30]=t16*c15-t15*s15;
#endif	/* USE_SSE2 */

	#ifdef DEBUG_SSE2
		for(iloop = 0; iloop < 64; iloop += 4)
		{
	  #ifdef USE_SSE2
			fprintf(stderr, "radix16_dyadic_square: A_out[%2d] = %20.5f, %20.5f\n",iloop,a[j1+iloop+0],a[j1+iloop+2]);
			fprintf(stderr, "radix16_dyadic_square: A_out[%2d] = %20.5f, %20.5f\n",iloop,a[j1+iloop+1],a[j1+iloop+3]);
	  #else
			fprintf(stderr, "radix16_dyadic_square: A_out[%2d] = %20.5f, %20.5f\n",iloop,a[j1+iloop+0],a[j1+iloop+1]);
			fprintf(stderr, "radix16_dyadic_square: A_out[%2d] = %20.5f, %20.5f\n",iloop,a[j1+iloop+2],a[j1+iloop+3]);
	  #endif
		}
		exit(0);
	#endif

	}
}

