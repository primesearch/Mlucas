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

/* Use for testing higher-accuracy version of the twiddles computation */
#define HIACC 1

#ifdef USE_SSE2

	#undef DEBUG_SSE2
//	#define DEBUG_SSE2

/*	Recipe for MSVC --> GCC inline ASM conversion:

Before you begin tranmslation:
	- Max. number of input variables GCC allows = 30 ... if you're using more than that,
	trying reducing the count e.g. by using var2 = var1 + memoffset in the ASM.
	DO THIS USING THE MSVC CODE, i.e. only *after* you've successfully reduced
	the inline ASM macro arg count should you proceed with syntax translation.
	That allows you to work through small chunks of inline ASM at a time, doing
	quick-build-and-debug to check the changes, i.e. greatly eases debug.

	0. Remove all but most-crucial comments to ease conversion, as follows:
		[blockmode] space all "keeper" comments to extreme left
		multistatement __asm lines --> one __asm pre line, realign __asm to left-justify, delete __asm\t
		For non-keeper comments: /* --> @@
		[regexp] @@ --> \t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t@@
		(delete all @@... stuff)\
		\t\n --> \n (repeat until no more trailing tabs)
		[/regexp]
		Repeat /* --> @@, [regexp] @@ --> \t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t@@ steps for keeper comments, to move out of harm's way.
	1. [...] --> (...)
	2. ALU ops [e.g. mov, add, shl] --> spaces to tabs, then append "l" [if e*x] or "q" [if r*x] to instruction name
	3. Numeric literals in above kinds of instructions: Prepend "$" [",0x" --> ",$0x"]
	4. Address offsets of form (...+0x100) --> 0x100(...), (...-0x100) --> -0x100(...)
	5. External variable names get wrapped in %[]
	6. Line up commas in vertically stacked columns, then reverse operand order columnwise [for both 2 and 3-operand instructions].
	7. Prepend "%%" to all register names
	8. Only e*x/r*x registers appear in clobber list, not special regs like mmx and xmm.

Additional Notes:
	- Need to strip off any leading white space from named vars inside [], e.g. for "movl %[  c4],%%ecx \n\t" get "undefined named operand '  c4'" error;
	- Offsets with explicit + sign, e.g. "+0x10(%%eax)", not allowed
*/

	#ifdef COMPILER_TYPE_MSVC

		#include "sse2_macro.h"

	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "radix16_dif_dit_pass_gcc32.h"

		#else

			#include "radix16_dif_dit_pass_gcc64.h"

		#endif

	#endif

#endif

/***************/

void radix16_dif_pass(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr)
{

/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform a single radix-16 complex DIF FFT pass on the data in the length-N real vector A.
!
!   The data are stored in a 1-D zero-offset array, with 2^PAD_BITS 8-byte padding elements inserted
!   between every block of 2^DAT_BITS contiguous data. The array padding is to prevent data being accessed
!   in strides that are large powers of two and thus to minimize cache thrashing
!   (in cache-based microprocessor architectures) or bank conflicts (in supercomputers.)
*/
	static double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;	/* exp[i*(twopi/16)]*/
	int i,j,j1,j2,jt,jp,jlo,jhi,m,iroot_prim,iroot,k1,k2;
	static int p1,p2,p3,p4,p8,p12;
	double rt,it;
	double re0,im0,re1,im1;
	double *addr, *addp;
	int prefetch_offset;

#ifdef USE_SSE2

	static int	first_entry = TRUE;
	static double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
	static struct complex *sc_arr = 0x0, *sc_ptr;
	static struct complex *cc0, *ss0, *isrt2, *two;
	static struct complex *c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7,*c8,*c9,*c10,*c11,*c12,*c13,*c14,*c15,*s0,*s1,*s2,*s3,*s4,*s5,*s6,*s7,*s8,*s9,*s10,*s11,*s12,*s13,*s14,*s15;
	static struct complex *r1,*r2,*r3,*r4,*r5,*r6,*r7,*r8,*r9,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r30,*r31,*r32;

#else

	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32
	,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15;

#endif

#ifdef USE_SSE2

	if(first_entry)
	{
		first_entry = FALSE;

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

	}	/* end of inits */

#endif

/*
!   Here's how our padded-array indexing works: given an unpadded-array index J,
!   we convert that to a padded-array index J1 using the formula
!
!	J1 = floor[J/(NBLOCK)] x [NBLOCK+PAD] + mod[J, NBLOCK],
!
!   where NBLOCK is the number of 8-byte array data in each contiguous-data block = 2^DAT_BITS,
!   and PAD is the number of 8-byte padding data between data block = 2^PAD_BITS.
!   Since N is an integer, we use the fact that mod[J, NBLOCK] = J - [NBLOCK] x floor[J/(NBLOCK)] to derive
!
!	J1 = J + PAD x floor[J/(NBLOCK)],
!
!   which expresses the fact that J1-J is PAD times the number of padding blocks we skip over.
!
!   Since NBLOCK is a power of 2, i.e. DAT_BITS = log_2(NBLOCK) is an integer, then floor[J/(NBLOCK)] = J >> DAT_BITS,
!   and J1 can be calculated using a shift, an integer multiply, and an add:
!
!	J1 = J + PAD x (J >> DAT_BITS).
!
!   If both NBLOCK and PAD are powers of 2 (as here), this can be done using two shifts and an add:
!
!	J1 = J + (J >> DAT_BITS) << PAD_BITS.
!
!   If PAD is 4 or 8, this can be done using a shift and a scaled add on Alpha.
!
!   Note that we only need to do this kind of calculation for the first elements of each block of 16 complex data.
!
!...The radix-16 pass is here. The data are processed in NLOOPS blocks of INCR real elements each
!   stride between elements processed together.
*/
	p1 = incr >> 4;
	p2 = p1 +p1;
	p3 = p2 +p1;
	p4 = p3 +p1;
	p8 = p4 +p4;
	p12= p8 +p4;

	p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
	p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
	p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );

	iroot_prim=(incr >> 5);		/* (incr/2)/radix_now */

	for(m=0; m < nloops; m++)	/* NLOOPS may range from 1 (if first pass radix = 16) to P*N/32 (last pass radix = 16).	 */
	{				/* NLOOPS satisfies the identity NLOOPS * INCR = P*N, which says that we process the entire */
					/* array each time this subroutine is executed (since P*N = vector length, sans padding.)   */

/*	here are the needed sincos data - these are processed below in bit-reversed order. */
	  iroot = index[m]*iroot_prim;
	  i = iroot;

/* In the DIF pass we may be able to afford a more-accurate twiddles computation,
   since we only compute each set of twiddles once and re-use it for many data blocks.
*/
#if HIACC
	#ifdef USE_SSE2
		c0->re=1.;	s0->re=0.;
		c0->im=1.;	s0->im=0.;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c1 ->re=rt;	s1 ->re=it;
		c1 ->im=rt;	s1 ->im=it;
	#else
		c1 =rt;		s1 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 3*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c2 ->re=rt;	s2 ->re=it;
		c2 ->im=rt;	s2 ->im=it;
	#else
		c2 =rt;		s2 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 4*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c3 ->re=rt;	s3 ->re=it;
		c3 ->im=rt;	s3 ->im=it;
	#else
		c3 =rt;		s3 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 5*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c4 ->re=rt;	s4 ->re=it;
		c4 ->im=rt;	s4 ->im=it;
	#else
		c4 =rt;		s4 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 6*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c5 ->re=rt;	s5 ->re=it;
		c5 ->im=rt;	s5 ->im=it;
	#else
		c5 =rt;		s5 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 7*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c6 ->re=rt;	s6 ->re=it;
		c6 ->im=rt;	s6 ->im=it;
	#else
		c6 =rt;		s6 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 8*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c7 ->re=rt;	s7 ->re=it;
		c7 ->im=rt;	s7 ->im=it;
	#else
		c7 =rt;		s7 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 9*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c8 ->re=rt;	s8 ->re=it;
		c8 ->im=rt;	s8 ->im=it;
	#else
		c8 =rt;		s8 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*10*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c9 ->re=rt;	s9 ->re=it;
		c9 ->im=rt;	s9 ->im=it;
	#else
		c9 =rt;		s9 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*11*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c10->re=rt;	s10->re=it;
		c10->im=rt;	s10->im=it;
	#else
		c10=rt;		s10=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*12*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c11->re=rt;	s11->re=it;
		c11->im=rt;	s11->im=it;
	#else
		c11=rt;		s11=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*13*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c12->re=rt;	s12->re=it;
		c12->im=rt;	s12->im=it;
	#else
		c12=rt;		s12=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*14*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c13->re=rt;	s13->re=it;
		c13->im=rt;	s13->im=it;
	#else
		c13=rt;		s13=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*15*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c14->re=rt;	s14->re=it;
		c14->im=rt;	s14->im=it;
	#else
		c14=rt;		s14=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c15->re=rt;	s15->re=it;
		c15->im=rt;	s15->im=it;
	#else
		c15=rt;		s15=it;
	#endif

#else
		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 2*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c1=t1*t3-t2*t4;	s1=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += i;				/* 4*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c2=t1*t3-t2*t4;	s2=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += i;				/* 8*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c4=t1*t3-t2*t4;	s4=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += (iroot << 2)+iroot;		/* 13*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c8=t1*t3-t2*t4;	s8=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c13=t1*t3-t2*t4;	s13=t1*t4+t2*t3;

		/* c3,5 */
		t1=c1*c4; t2=c1*s4; t3=s1*c4; t4=s1*s4;
		c3=t1+t4; s3=t2-t3; c5=t1-t4; s5=t2+t3;

		/* c6,7,9,10 */
		t1=c1*c8; t2=c1*s8; t3=s1*c8; t4=s1*s8;
		c7=t1+t4; s7=t2-t3; c9=t1-t4; s9=t2+t3;

		t1=c2*c8; t2=c2*s8; t3=s2*c8; t4=s2*s8;
		c6=t1+t4; s6=t2-t3; c10=t1-t4; s10=t2+t3;

		/* c11,12,14,15 */
		t1=c1*c13; t2=c1*s13; t3=s1*c13; t4=s1*s13;
		c12=t1+t4; s12=t2-t3; c14=t1-t4; s14=t2+t3;

		t1=c2*c13; t2=c2*s13; t3=s2*c13; t4=s2*s13;
		c11=t1+t4; s11=t2-t3; c15=t1-t4; s15=t2+t3;
#endif

	/* Define the inner-loop parameters in terms of the outer-loop ones to make OpenMP's job easier: */
	  jlo = m*incr;
	  jhi = jlo+(incr >> 4);

	#ifdef USE_SSE2
	  for(j=jlo; j < jhi; j += 4)
	  {
		/* In SSE2 mode, data are arranged in [re0,re1,im0,im1] quartets, not the usual [re0,im0],[re1,im1] pairs.
		Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
		but to point to the proper real datum, we need to bit-reverse bits <0:1> of j, i.e. [0,1,2,3] ==> [0,2,1,3].
		*/
		j1 = (j & mask01) + br4[j&3];
	#elif defined(USE_SSE2)	/* This allows us to use #if 0 above and disable sse2-based *computation*, while still using sse2-style data layout */
	  for(j=jlo; j < jhi; j += 2)
	  {
		j1 = (j & mask01) + br4[j&3];
	#else
	  for(j=jlo; j < jhi; j += 2)
	  {
		j1 =  j;
	#endif
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;

#ifdef USE_SSE2

	#ifdef DEBUG_SSE2
		rng_isaac_init(TRUE);
		jt = j1;		jp = j2;
		for(i = 0; i < 16; i++)
		{
			a[jt] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			fprintf(stderr, "radix16_dif_pass: A_in[%2d] = %20.5f, %20.5f\n",i,a[jt],a[jp]);
			jt += p1;	jp += p1;
		}
	#endif

	#ifdef COMPILER_TYPE_MSVC

	  #if 1	// if(1) - test out pure-asm version

		add0 = &a[j1];
		__asm	mov	eax, add0
		__asm	mov	ebx, p4		// Can't get these via simple load-one-and-shift-as-needed due to array padding scheme
		__asm	mov	ecx, p8
		__asm	mov	edx, p12
		__asm	shl	ebx,  3
		__asm	shl	ecx,  3
		__asm	shl	edx,  3
		__asm	add	ebx, eax
		__asm	add	ecx, eax
		__asm	add	edx, eax
		SSE2_RADIX4_DIF_4TWIDDLE_B(r1 ,c0)
		__asm	mov	edi, p2		// Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p1+p1 not guaranteed.
		__asm	shl	edi,  3		// 8-bytes for array-of-doubles
		__asm	add	eax, edi	// &a[j1+p2];
		__asm	add	ebx, edi
		__asm	add	ecx, edi
		__asm	add	edx, edi
		SSE2_RADIX4_DIF_4TWIDDLE_B(r9 ,c2)
		__asm	sub	eax, edi	// &a[j1];
		__asm	sub	ebx, edi
		__asm	sub	ecx, edi
		__asm	sub	edx, edi
		__asm	mov	edi, p1
		__asm	shl	edi,  3		// 8-bytes for array-of-doubles
		__asm	add	eax, edi	// &a[j1+p1];
		__asm	add	ebx, edi
		__asm	add	ecx, edi
		__asm	add	edx, edi
		SSE2_RADIX4_DIF_4TWIDDLE_B(r17,c1)
		__asm	sub	eax, edi	// &a[j1];
		__asm	sub	ebx, edi
		__asm	sub	ecx, edi
		__asm	sub	edx, edi
		__asm	mov	edi, p3
		__asm	shl	edi,  3		// 8-bytes for array-of-doubles
		__asm	add	eax, edi	// &a[j1+p3];
		__asm	add	ebx, edi
		__asm	add	ecx, edi
		__asm	add	edx, edi
		SSE2_RADIX4_DIF_4TWIDDLE_B(r25,c3)

	  #else

	/*...Block 1: */
		add0 = &a[j1];
		add1 = add0+p4;
		add2 = add0+p8;
		add3 = add0+p12;

		/* Do the p4,12 combo first. 	Cost: 24 MOVapd, 28 ADD/SUBpd, 12 MULpd */
		__asm	mov	eax, add1
		__asm	mov	ebx, c4
		__asm	mov	ecx, add3
		__asm	mov	edx, c12

		__asm	movaps	xmm0,[eax     ]	/* a[jt+p4] */				__asm	movaps	xmm4,[ecx     ]	/* a[jt+p12] */
		__asm	movaps	xmm1,[eax+0x10]	/* a[jp+p4] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p12] */
		__asm	movaps	xmm2,[eax     ]	/* xmm2 <- cpy a[jt+p4] */	__asm	movaps	xmm6,[ecx     ]	/* xmm6 <- cpy a[jt+p12] */
		__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[jp+p4] */	__asm	movaps	xmm7,[ecx+0x10]	/* xmm7 <- cpy a[jp+p12] */

		__asm	mulpd	xmm0,[ebx     ]	/* a[jt+p4]*c4 */			__asm	mulpd	xmm4,[edx     ]	/* a[jt+p12]*c12 */
		__asm	mulpd	xmm1,[ebx     ]	/* a[jp+p4]*c4 */			__asm	mulpd	xmm5,[edx     ]	/* a[jp+p12]*c12 */
		__asm	mulpd	xmm2,[ebx+0x10]	/* a[jt+p4]*s4 */			__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p12]*s12 */
		__asm	mulpd	xmm3,[ebx+0x10]	/* a[jp+p4]*s4 */			__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p12]*s12 */
		__asm	addpd	xmm1,xmm2	/* xmm1 <- t6 */				__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm0,xmm3	/* xmm0 <- t5 */				__asm	subpd	xmm4,xmm7	/* xmm4 <- rt 	xmm6,7 free */
		__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy t6 */
		__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy t5 */

		__asm	addpd	xmm0,xmm4	/* ~t5 <- t5 +rt */
		__asm	addpd	xmm1,xmm5	/* ~t6 <- t6 +it */
		__asm	subpd	xmm2,xmm4	/* ~t7 <- t5 -rt */
		__asm	subpd	xmm3,xmm5	/* ~t8 <- t6 -it	xmm4,5 free */

		/* Now do the p0,8 combo: */
		__asm	mov	eax, add0
		__asm	mov	ecx, add2
		__asm	mov	edx, c8

		__asm	movaps	xmm4,[ecx     ]	/* a[jt+p8] */
		__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p8] */
		__asm	movaps	xmm6,[ecx     ]	/* xmm6 <- cpy a[jt+p8] */
		__asm	movaps	xmm7,[ecx+0x10]	/* xmm7 <- cpy a[jp+p8] */

		__asm	mulpd	xmm4,[edx     ]	/* a[jt+p8]*c8 */
		__asm	mulpd	xmm5,[edx     ]	/* a[jp+p8]*c8 */
		__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p8]*s8 */
		__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p8]*s8 */
		__asm	addpd	xmm5,xmm6	/* xmm5 <- it */
		__asm	subpd	xmm4,xmm7	/* xmm4 <- rt 	xmm6,7 free - stick t1,t2 in those */

		__asm	movaps	xmm6,[eax     ]	/* t1 = a[jt   ] */
		__asm	movaps	xmm7,[eax+0x10]	/* t2 = a[jp   ] */

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

	/*...Block 2: */
		add0 = &a[j1+p2];
		add1 = add0+p4;
		add2 = add0+p8;
		add3 = add0+p12;
		/* 	Cost: 30 MOVapd, 28 ADD/SUBpd, 16 MULpd */
		SSE2_RADIX4_DIF_4TWIDDLE(add0,add1,add2,add3,r9 ,c2)

	/*...Block 3: */
		add0 = &a[j1+p1];
		add1 = add0+p4;
		add2 = add0+p8;
		add3 = add0+p12;
		/* 	Cost: 30 MOVapd, 28 ADD/SUBpd, 16 MULpd */
		SSE2_RADIX4_DIF_4TWIDDLE(add0,add1,add2,add3,r17,c1 )

	/*...Block 4: */
		add0 = &a[j1+p3];
		add1 = add0+p4;
		add2 = add0+p8;
		add3 = add0+p12;
		/* 	Cost: 30 MOVapd, 28 ADD/SUBpd, 16 MULpd */
		SSE2_RADIX4_DIF_4TWIDDLE(add0,add1,add2,add3,r25,c3)

	  #endif	/* if(1) */

	/**************************************************************************************
	!...and now do four more radix-4 transforms, including the internal twiddle factors:  !
	**************************************************************************************/
	/*...Block 1: t1,9,17,25	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
	  #if 1
		__asm	mov eax, add0	// &a[j1]
		__asm	mov ebx, p1
		__asm	mov ecx, p2
		__asm	mov edx, p3
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	shl	ecx, 3
		__asm	shl	edx, 3
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
	  #else
		add0 = &a[j1];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		__asm	mov	eax, add0	/* restore main-array indices */
		__asm	mov	ebx, add1
		__asm	mov	ecx, add2
		__asm	mov	edx, add3
	  #endif
		__asm	mov	edi, r1
		__asm	movaps	xmm0,[edi      ]	/* t1  */			__asm	movaps	xmm4,[edi+0x100]	/* t17 */
		__asm	movaps	xmm1,[edi+0x010]	/* t2  */			__asm	movaps	xmm5,[edi+0x110]	/* t18 */
		__asm	movaps	xmm2,[edi+0x080]	/* t9  */			__asm	movaps	xmm6,[edi+0x180]	/* t25 */
		__asm	movaps	xmm3,[edi+0x090]	/* t10 */			__asm	movaps	xmm7,[edi+0x190]	/* t26 */

		__asm	subpd	xmm0,xmm2		/* ~t9 = t1-t9 */		__asm	subpd	xmm4,xmm6		/* ~t25=t17-t25*/
		__asm	subpd	xmm1,xmm3		/* ~t10= t2-t10*/		__asm	subpd	xmm5,xmm7		/* ~t26=t18-t26*/
		__asm	addpd	xmm2,xmm2		/*        2*t9 */		__asm	addpd	xmm6,xmm6		/*        2*t25*/
		__asm	addpd	xmm3,xmm3		/*        2*t10*/		__asm	addpd	xmm7,xmm7		/*        2*t26*/
		__asm	addpd	xmm2,xmm0		/* ~t1 = t1+t9 */		__asm	addpd	xmm6,xmm4		/* ~t17=t17+t25*/
		__asm	addpd	xmm3,xmm1		/* ~t2 = t2+t10*/		__asm	addpd	xmm7,xmm5		/* ~t18=t18+t26*/

		__asm	subpd	xmm2,xmm6		/* t1  <- t1 -t17 */
		__asm	subpd	xmm3,xmm7		/* t2  <- t2 -t18 */
		__asm	addpd	xmm6,xmm6		/*          2*t17 */
		__asm	addpd	xmm7,xmm7		/*          2*t18 */
		__asm	movaps	[ebx     ],xmm2	/* a[jt+p1 ] */
		__asm	movaps	[ebx+0x10],xmm3	/* a[jp+p1 ] */
		__asm	addpd	xmm6,xmm2		/* t17 <- t1 +t17 */
		__asm	addpd	xmm7,xmm3		/* t18 <- t2 +t18 */
		__asm	movaps	[eax     ],xmm6	/* a[jt+p0 ] */
		__asm	movaps	[eax+0x10],xmm7	/* a[jp+p0 ] */

		__asm	subpd	xmm0,xmm5		/* t9  <- t9 -t26 */
		__asm	subpd	xmm1,xmm4		/* t10 <- t10-t25 */
		__asm	addpd	xmm5,xmm5		/*          2*t26 */
		__asm	addpd	xmm4,xmm4		/*          2*t25 */
		__asm	movaps	[ecx     ],xmm0	/* a[jt+p2 ] */
		__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
		__asm	addpd	xmm5,xmm0		/* t26 <- t9 +t26 */
		__asm	addpd	xmm4,xmm1		/* t25 <- t10+t25 */
		__asm	movaps	[edx     ],xmm5	/* a[jt+p3 ] */
		__asm	movaps	[ecx+0x10],xmm4	/* a[jp+p2 ] */

	/*...Block 3: t5,13,21,29	Cost: 16 MOVapd, 26 ADD/SUBpd,  4 MULpd */
	  #if 1
		__asm	mov	edi, p4
		__asm	shl	edi, 3
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
	  #else
		add0 = &a[j1+p4];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		__asm	mov	eax, add0	/* restore main-array indices */
		__asm	mov	ebx, add1
		__asm	mov	ecx, add2
		__asm	mov	edx, add3
	  #endif
		__asm	mov	edi, r5
		__asm	mov	esi, isrt2
		__asm	movaps	xmm3,[esi ]	/* isrt2 */
																__asm	movaps	xmm4,[edi+0x100]	/* t21 */
																__asm	movaps	xmm5,[edi+0x110]	/* t22 */
																__asm	movaps	xmm6,[edi+0x180]	/* t29 */
																__asm	movaps	xmm7,[edi+0x190]	/* t30 */
																__asm	mulpd	xmm4,xmm3	/* t21 *ISRT2 */
		__asm	movaps	xmm0,[edi      ]	/* t5  */			__asm	mulpd	xmm5,xmm3	/* t22 *ISRT2 */
		__asm	movaps	xmm1,[edi+0x010]	/* t6  */			__asm	mulpd	xmm6,xmm3	/* t29 *ISRT2 */
		__asm	movaps	xmm2,[edi+0x080]	/* t13 */			__asm	mulpd	xmm7,xmm3	/* t30 *ISRT2 */
		__asm	movaps	xmm3,[edi+0x090]	/* t14; this must execute after the last mul-by-ISRT2 above */

		__asm	subpd	xmm0,xmm3		/* ~t5 = t5 -t14*/		__asm	subpd	xmm4,xmm5		/* ~t21=t21-t22*/
		__asm	subpd	xmm1,xmm2		/* ~t14= t6 -t13*/		__asm	subpd	xmm7,xmm6		/*  it =t30-t29*/
		__asm	addpd	xmm3,xmm3		/*         2*t14*/		__asm	addpd	xmm5,xmm5		/*        2*t22*/
		__asm	addpd	xmm2,xmm2		/*         2*t13*/		__asm	addpd	xmm6,xmm6		/*        2*t29*/
		__asm	addpd	xmm3,xmm0		/* ~t13= t14+t5 */		__asm	addpd	xmm5,xmm4		/* ~t22=t22+t21*/
		__asm	addpd	xmm2,xmm1		/* ~t6 = t13+t6 */		__asm	addpd	xmm6,xmm7		/*  rt =t29+t30*/

		__asm	subpd	xmm4,xmm6		/* t21=t21-rt */
		__asm	subpd	xmm5,xmm7		/* t22=t22-it */
		__asm	addpd	xmm6,xmm6		/*      2* rt */
		__asm	addpd	xmm7,xmm7		/*      2* it */
		__asm	addpd	xmm6,xmm4		/* t29=t21+rt */
		__asm	addpd	xmm7,xmm5		/* t30=t22+it */

		__asm	subpd	xmm0,xmm4		/* t5 -t21 */			__asm	subpd	xmm3,xmm7		/* t13-t30 */
		__asm	subpd	xmm2,xmm5		/* t6 -t22 */			__asm	subpd	xmm1,xmm6		/* t14-t29 */
		__asm	addpd	xmm4,xmm4		/*   2*t21 */			__asm	addpd	xmm7,xmm7		/*   2*t30 */
		__asm	addpd	xmm5,xmm5		/*   2*t22 */			__asm	addpd	xmm6,xmm6		/*   2*t29 */
		__asm	movaps	[ebx     ],xmm0	/* a[jt+p1 ] */			__asm	movaps	[ecx     ],xmm3	/* a[jt+p2 ] */
		__asm	movaps	[ebx+0x10],xmm2	/* a[jp+p1 ] */			__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
		__asm	addpd	xmm4,xmm0		/* t5 +t21 */			__asm	addpd	xmm7,xmm3		/* t13+t30 */
		__asm	addpd	xmm5,xmm2		/* t6 +t22 */			__asm	addpd	xmm6,xmm1		/* t14+t29 */
		__asm	movaps	[eax     ],xmm4	/* a[jt+p0 ] */			__asm	movaps	[edx     ],xmm7	/* a[jt+p3 ] */
		__asm	movaps	[eax+0x10],xmm5	/* a[jp+p0 ] */			__asm	movaps	[ecx+0x10],xmm6	/* a[jp+p2 ] */

	/*...Block 2: t3,11,19,27	Cost: 18 MOVapd, 28 ADD/SUBpd, 10 MULpd */
	  #if 1
		__asm	mov	edi, p4
		__asm	shl	edi, 3
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
	  #else
		add0 = &a[j1+p8];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		__asm	mov	eax, add0	/* restore main-array indices */
		__asm	mov	ebx, add1
		__asm	mov	ecx, add2
		__asm	mov	edx, add3
	  #endif
		__asm	mov	edi, r3
		__asm	mov	esi, cc0
		__asm	movaps	xmm4,[edi+0x100]	/* t19 */
		__asm	movaps	xmm5,[edi+0x110]	/* t20 */
		__asm	movaps	xmm3,[esi      ]	/* c */
		__asm	movaps	xmm2,[esi+0x010]	/* s */
		__asm	movaps	xmm6,xmm4		/* copy t19 */
		__asm	movaps	xmm7,xmm5		/* copy t20 */

		__asm	mulpd	xmm4,xmm3		/* t19*c */
		__asm	mulpd	xmm5,xmm3		/* t20*c */
		__asm	mulpd	xmm6,xmm2		/* t19*s */				__asm	movaps	xmm0,[edi+0x180]	/* t27 */
		__asm	mulpd	xmm7,xmm2		/* t20*s */				__asm	movaps	xmm1,[edi+0x190]	/* t28 */
		__asm	addpd	xmm5,xmm6	/* ~t20 */					__asm	movaps	xmm6,xmm0		/* copy t27 */
		__asm	subpd	xmm4,xmm7	/* ~t19 */					__asm	movaps	xmm7,xmm1		/* copy t28 */

																__asm	mulpd	xmm6,xmm2		/* t27*s */
																__asm	mulpd	xmm7,xmm2		/* t28*s */
																__asm	mulpd	xmm0,xmm3		/* t27*c */
																__asm	mulpd	xmm1,xmm3		/* t28*c */
																__asm	addpd	xmm7,xmm0	/* it */
																__asm	subpd	xmm6,xmm1	/* rt */

		__asm	movaps	xmm2,xmm4		/* copy t19 */
		__asm	movaps	xmm3,xmm5		/* copy t20 */
		__asm	subpd	xmm4,xmm6		/*~t27=t19-rt */
		__asm	subpd	xmm5,xmm7		/*~t28=t20-it */
		__asm	addpd	xmm6,xmm2		/*~t19=t19+rt */
		__asm	addpd	xmm7,xmm3		/*~t20=t20+it */

		__asm	mov	esi, isrt2
		__asm	movaps	xmm2,[edi+0x080]	/* t11 */
		__asm	movaps	xmm3,[edi+0x090]	/* t12 */
		__asm	movaps	xmm1,[esi]	/* isrt2 */
		__asm	movaps	xmm0,xmm2	/* cpy t11 */
		__asm	subpd	xmm2,xmm3	/*~t11=t11-t12 */
		__asm	addpd	xmm3,xmm0	/*~t12=t12+t11 */
		__asm	mulpd	xmm2,xmm1	/* rt = (t11-t12)*ISRT2 */
		__asm	mulpd	xmm3,xmm1	/* it = (t12+t11)*ISRT2 */

		__asm	movaps	xmm0,[edi      ]	/* t3  */
		__asm	movaps	xmm1,[edi+0x010]	/* t4  */

		__asm	subpd	xmm0,xmm2			/*~t11=t3 -rt */
		__asm	subpd	xmm1,xmm3			/*~t12=t4 -it */
		__asm	addpd	xmm2,xmm2			/*      2* rt */
		__asm	addpd	xmm3,xmm3			/*      2* it */
		__asm	addpd	xmm2,xmm0			/*~t3 =rt +t3 */
		__asm	addpd	xmm3,xmm1			/*~t4 =it +t4 */

		__asm	subpd	xmm2,xmm6		/* t3 -t19 */			__asm	subpd	xmm0,xmm5		/* t11-t28 */
		__asm	subpd	xmm3,xmm7		/* t4 -t20 */			__asm	subpd	xmm1,xmm4		/* t12-t27 */
		__asm	addpd	xmm6,xmm6		/*   2*t19 */			__asm	addpd	xmm5,xmm5		/*          2*t28 */
		__asm	addpd	xmm7,xmm7		/*   2*t20 */			__asm	addpd	xmm4,xmm4		/*          2*t27 */
		__asm	movaps	[ebx     ],xmm2	/* a[jt+p1 ] */			__asm	movaps	[ecx     ],xmm0	/* a[jt+p2 ] */
		__asm	movaps	[ebx+0x10],xmm3	/* a[jp+p1 ] */			__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
		__asm	addpd	xmm6,xmm2		/* t3 +t19 */			__asm	addpd	xmm5,xmm0		/* t11+t28 */
		__asm	addpd	xmm7,xmm3		/* t4 +t20 */			__asm	addpd	xmm4,xmm1		/* t12+t27 */
		__asm	movaps	[eax     ],xmm6	/* a[jt+p0 ] */			__asm	movaps	[edx     ],xmm5	/* a[jt+p3 ] */
		__asm	movaps	[eax+0x10],xmm7	/* a[jp+p0 ] */			__asm	movaps	[ecx+0x10],xmm4	/* a[jp+p2 ] */

	/*...Block 4: t7,15,23,31	Cost: 18 MOVapd, 28 ADD/SUBpd, 10 MULpd */
	  #if 1
		__asm	mov	edi, p4
		__asm	shl	edi, 3
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
	  #else
		add0 = &a[j1+p12];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		__asm	mov	eax, add0	/* restore main-array indices */
		__asm	mov	ebx, add1
		__asm	mov	ecx, add2
		__asm	mov	edx, add3
	  #endif
		__asm	mov	edi, r7
		__asm	mov	esi, cc0
		__asm	movaps	xmm4,[edi+0x100]	/* t23 */
		__asm	movaps	xmm5,[edi+0x110]	/* t24 */
		__asm	movaps	xmm2,[esi      ]	/* c */
		__asm	movaps	xmm3,[esi+0x010]	/* s */
		__asm	movaps	xmm6,xmm4		/* copy t23 */
		__asm	movaps	xmm7,xmm5		/* copy t24 */

		__asm	mulpd	xmm4,xmm3		/* t23*s */
		__asm	mulpd	xmm5,xmm3		/* t24*s */
		__asm	mulpd	xmm6,xmm2		/* t23*c */				__asm	movaps	xmm0,[edi+0x180]	/* t31 */
		__asm	mulpd	xmm7,xmm2		/* t24*c */				__asm	movaps	xmm1,[edi+0x190]	/* t32 */
		__asm	addpd	xmm5,xmm6	/* ~t24 */					__asm	movaps	xmm6,xmm0		/* copy t31 */
		__asm	subpd	xmm4,xmm7	/* ~t23 */					__asm	movaps	xmm7,xmm1		/* copy t32 */

																__asm	mulpd	xmm6,xmm2		/* t31*c */
																__asm	mulpd	xmm7,xmm2		/* t32*c */
																__asm	mulpd	xmm0,xmm3		/* t31*s */
																__asm	mulpd	xmm1,xmm3		/* t32*s */
																__asm	addpd	xmm7,xmm0	/* it */
																__asm	subpd	xmm6,xmm1	/* rt */

		__asm	movaps	xmm2,xmm4		/* copy t23 */
		__asm	movaps	xmm3,xmm5		/* copy t24 */
		__asm	subpd	xmm4,xmm6		/*~t23=t23-rt */
		__asm	subpd	xmm5,xmm7		/*~t24=t24-it */
		__asm	addpd	xmm6,xmm2		/*~t31=t23+rt */
		__asm	addpd	xmm7,xmm3		/*~t32=t24+it */

		__asm	mov	esi, isrt2
		__asm	movaps	xmm2,[edi+0x080]	/* t15 */
		__asm	movaps	xmm3,[edi+0x090]	/* t16 */
		__asm	movaps	xmm1,[esi]	/* isrt2 */
		__asm	movaps	xmm0,xmm2	/* cpy t15 */
		__asm	addpd	xmm2,xmm3	/*~t15=t15+t16 */
		__asm	subpd	xmm3,xmm0	/*~t16=t16-t15 */
		__asm	mulpd	xmm2,xmm1	/* rt = (t15+t16)*ISRT2 */
		__asm	mulpd	xmm3,xmm1	/* it = (t16-t15)*ISRT2 */

		__asm	movaps	xmm0,[edi      ]	/* t7  */
		__asm	movaps	xmm1,[edi+0x010]	/* t8  */

		__asm	subpd	xmm0,xmm2			/*~t7 =t7 -rt */
		__asm	subpd	xmm1,xmm3			/*~t8 =t8 -it */
		__asm	addpd	xmm2,xmm2			/*      2* rt */
		__asm	addpd	xmm3,xmm3			/*      2* it */
		__asm	addpd	xmm2,xmm0			/*~t15=rt +t7 */
		__asm	addpd	xmm3,xmm1			/*~t16=it +t8 */

		__asm	subpd	xmm0,xmm4		/* t7 -t23 */			__asm	subpd	xmm2,xmm7		/* t15-t32 */
		__asm	subpd	xmm1,xmm5		/* t8 -t24 */			__asm	subpd	xmm3,xmm6		/* t16-t31 */
		__asm	addpd	xmm4,xmm4		/*   2*t23 */			__asm	addpd	xmm7,xmm7		/*   2*t32 */
		__asm	addpd	xmm5,xmm5		/*   2*t24 */			__asm	addpd	xmm6,xmm6		/*   2*t31 */
		__asm	movaps	[ebx     ],xmm0	/* a[jt+p1 ] */			__asm	movaps	[ecx     ],xmm2	/* a[jt+p2 ] */
		__asm	movaps	[ebx+0x10],xmm1	/* a[jp+p1 ] */			__asm	movaps	[edx+0x10],xmm3	/* a[jp+p3 ] */
		__asm	addpd	xmm4,xmm0		/* t7 +t23 */			__asm	addpd	xmm7,xmm2		/* t15+t32 */
		__asm	addpd	xmm5,xmm1		/* t8 +t24 */			__asm	addpd	xmm6,xmm3		/* t16+t31 */
		__asm	movaps	[eax     ],xmm4	/* a[jt+p0 ] */			__asm	movaps	[edx     ],xmm7	/* a[jt+p3 ] */
		__asm	movaps	[eax+0x10],xmm5	/* a[jp+p0 ] */			__asm	movaps	[ecx+0x10],xmm6	/* a[jp+p2 ] */

		/***************************************************/
		/* Total Cost: 182 MOVapd, 214 ADD/SUBpd, 84 MULpd */
		/***************************************************/

	#else	/* GCC-style inline ASM: */

		add0 = &a[j1];
		SSE2_RADIX16_DIF_TWIDDLE(add0,p1,p2,p3,p4,p8,p12,r1,r3,r5,r7,r9,r17,r25,isrt2,cc0,c0,c1,c2,c3);

	#endif

	#ifdef DEBUG_SSE2
			jt = j1;		jp = j2;
			fprintf(stderr, "radix16_dit: Aout[0] = %20.5f, %20.5f\n",a[jt   ],a[jp   ]);
			fprintf(stderr, "radix16_dit: Aout[1] = %20.5f, %20.5f\n",a[jt+p1],a[jp+p1]);
			fprintf(stderr, "radix16_dit: Aout[2] = %20.5f, %20.5f\n",a[jt+p2],a[jp+p2]);
			fprintf(stderr, "radix16_dit: Aout[3] = %20.5f, %20.5f\n",a[jt+p3],a[jp+p3]);	jt = j1 + p4;	jp = j2 + p4;
			fprintf(stderr, "radix16_dit: Aout[4] = %20.5f, %20.5f\n",a[jt   ],a[jp   ]);
			fprintf(stderr, "radix16_dit: Aout[5] = %20.5f, %20.5f\n",a[jt+p1],a[jp+p1]);
			fprintf(stderr, "radix16_dit: Aout[6] = %20.5f, %20.5f\n",a[jt+p2],a[jp+p2]);
			fprintf(stderr, "radix16_dit: Aout[7] = %20.5f, %20.5f\n",a[jt+p3],a[jp+p3]);	jt = j1 + p8;	jp = j2 + p8;
			fprintf(stderr, "radix16_dit: Aout[8] = %20.5f, %20.5f\n",a[jt   ],a[jp   ]);
			fprintf(stderr, "radix16_dit: Aout[9] = %20.5f, %20.5f\n",a[jt+p1],a[jp+p1]);
			fprintf(stderr, "radix16_dit: Aout[A] = %20.5f, %20.5f\n",a[jt+p2],a[jp+p2]);
			fprintf(stderr, "radix16_dit: Aout[B] = %20.5f, %20.5f\n",a[jt+p3],a[jp+p3]);	jt = j1 + p12;	jp = j2 + p12;
			fprintf(stderr, "radix16_dit: Aout[C] = %20.5f, %20.5f\n",a[jt   ],a[jp   ]);
			fprintf(stderr, "radix16_dit: Aout[D] = %20.5f, %20.5f\n",a[jt+p1],a[jp+p1]);
			fprintf(stderr, "radix16_dit: Aout[E] = %20.5f, %20.5f\n",a[jt+p2],a[jp+p2]);
			fprintf(stderr, "radix16_dit: Aout[F] = %20.5f, %20.5f\n",a[jt+p3],a[jp+p3]);
			exit(0);
	#endif

#else	/* USE_SSE2 */

	/* gather the needed data (8 64-bit complex, i.e. 16 64-bit reals) and do the first set of four length-4 transforms.
	   We process the sincos data in bit-reversed order.	*/
	#ifdef PFETCH_AGGRESSIVE
		addp = &a[j1];
	#elif PFETCH
		prefetch_offset = ((j >> 1) & 3)*p4 + 4;	/* Cycle among p0, p4, p8 and p12. */
	#endif
	/*...Block 1: */

		jt = j1;
		jp = j2;

	#ifdef PFETCH_AGGRESSIVE
		addr = addp;
		prefetch_p_doubles(addr);
	#elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	#endif
		t1 =a[jt	];						t2 =a[jp	];
		rt =a[jt+p8 ]*c8 -a[jp+p8 ]*s8 ;	it =a[jp+p8 ]*c8 +a[jt+p8 ]*s8;
		t3 =t1 -rt;							t1 =t1 +rt;
		t4 =t2 -it;							t2 =t2 +it;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		t5 =a[jt+p4 ]*c4 -a[jp+p4 ]*s4 ;	t6 =a[jp+p4 ]*c4 +a[jt+p4 ]*s4;
		rt =a[jt+p12]*c12-a[jp+p12]*s12;	it =a[jp+p12]*c12+a[jt+p12]*s12;
		t7 =t5 -rt;							t5 =t5 +rt;
		t8 =t6 -it;							t6 =t6 +it;

		rt =t5;	t5 =t1 -rt;					t1 =t1 +rt;
		it =t6;	t6 =t2 -it;					t2 =t2 +it;

		rt =t7;	t7 =t3 +t8;					t3 =t3 -t8;
				t8 =t4 -rt;					t4 =t4 +rt;

	/*...Block 2: */

		jt = j1 + p2;
		jp = j2 + p2;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	#endif
		t9 =a[jt    ]*c2 -a[jp    ]*s2 ;	t10=a[jp    ]*c2 +a[jt    ]*s2;
		rt =a[jt+p8 ]*c10-a[jp+p8 ]*s10;	it =a[jp+p8 ]*c10+a[jt+p8 ]*s10;
		t11=t9 -rt;							t9 =t9 +rt;
		t12=t10-it;							t10=t10+it;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		t13=a[jt+p4 ]*c6 -a[jp+p4 ]*s6 ;	t14=a[jp+p4 ]*c6 +a[jt+p4 ]*s6;
		rt =a[jt+p12]*c14-a[jp+p12]*s14;	it =a[jp+p12]*c14+a[jt+p12]*s14;
		t15=t13-rt;							t13=t13+rt;
		t16=t14-it;							t14=t14+it;

		rt =t13;	t13=t9 -rt;				t9 =t9 +rt;
		it =t14;	t14=t10-it;				t10=t10+it;

		rt =t15;	t15=t11+t16;			t11=t11-t16;
					t16=t12-rt;				t12=t12+rt;

	/*...Block 3: */

		jt = j1 + p1;
		jp = j2 + p1;

	#ifdef PFETCH_AGGRESSIVE
		addr = addp + p4;
		prefetch_p_doubles(addr);
	#elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	#endif
		t17=a[jt    ]*c1 -a[jp    ]*s1 ;	t18=a[jp    ]*c1 +a[jt    ]*s1;
		rt =a[jt+p8 ]*c9 -a[jp+p8 ]*s9 ;	it =a[jp+p8 ]*c9 +a[jt+p8 ]*s9;
		t19=t17-rt;							t17=t17+rt;
		t20=t18-it;							t18=t18+it;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		t21=a[jt+p4 ]*c5 -a[jp+p4 ]*s5 ;	t22=a[jp+p4 ]*c5 +a[jt+p4 ]*s5;
		rt =a[jt+p12]*c13-a[jp+p12]*s13;	it =a[jp+p12]*c13+a[jt+p12]*s13;
		t23=t21-rt;							t21=t21+rt;
		t24=t22-it;							t22=t22+it;

		rt =t21;	t21=t17-rt;				t17=t17+rt;
		it =t22;	t22=t18-it;				t18=t18+it;

		rt =t23;	t23=t19+t24;			t19=t19-t24;
					t24=t20-rt;				t20=t20+rt;

	/*...Block 4: */

		jt = j1 + p3;
		jp = j2 + p3;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	#endif
		t25=a[jt    ]*c3 -a[jp    ]*s3 ;	t26=a[jp    ]*c3 +a[jt    ]*s3;
		rt =a[jt+p8 ]*c11-a[jp+p8 ]*s11;	it =a[jp+p8 ]*c11+a[jt+p8 ]*s11;
		t27=t25-rt;							t25=t25+rt;
		t28=t26-it;							t26=t26+it;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		t29=a[jt+p4 ]*c7 -a[jp+p4 ]*s7 ;	t30=a[jp+p4 ]*c7 +a[jt+p4 ]*s7;
		rt =a[jt+p12]*c15-a[jp+p12]*s15;	it =a[jp+p12]*c15+a[jt+p12]*s15;
		t31=t29-rt;							t29=t29+rt;
		t32=t30-it;							t30=t30+it;

		rt =t29;	t29=t25-rt;				t25=t25+rt;
		it =t30;	t30=t26-it;				t26=t26+it;

		rt =t31;	t31=t27+t32;			t27=t27-t32;
					t32=t28-rt;				t28=t28+rt;
	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	!	1, exp(i* 1*twopi/16) =		( c, s), exp(i* 2*twopi/16) = isqrt2*( 1, 1), exp(i* 3*twopi/16) =		( s, c) (for inputs to transform block 2)
	!	1, exp(i* 2*twopi/16) = isqrt2*( 1, 1), exp(i* 4*twopi/16) =		( 0, 1), exp(i* 6*twopi/16) = isqrt2*(-1, 1) (for inputs to transform block 3)
	!	1, exp(i* 3*twopi/16) =		( s, c), exp(i* 6*twopi/16) = isqrt2*(-1, 1), exp(i* 9*twopi/16) =		(-c,-s) (for inputs to transform block 4).
	!  (This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
	!   I.e. do similar as above, except inputs a[j1,2+p0:15:1) are replaced by t0:30:2,
	!										   a[j1,2+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
	!   and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.
	*/
	/*...Block 1: t1,9,17,25 */

		jt = j1;
		jp = j2;

	#ifdef PFETCH_AGGRESSIVE
		addr = addp + p8 ;
		prefetch_p_doubles(addr);
	#endif
		rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;t10=t2 -it;	t2 =t2 +it;

		rt =t25;t25=t17-rt;	t17=t17+rt;
		it =t26;t26=t18-it;	t18=t18+it;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		a[jt    ]=t1+t17;	a[jp    ]=t2+t18;
		a[jt+p1 ]=t1-t17;	a[jp+p1 ]=t2-t18;

		a[jt+p2 ]=t9 -t26;	a[jp+p2 ]=t10+t25;	/* mpy by E^4=i is inlined here... */
		a[jt+p3 ]=t9 +t26;	a[jp+p3 ]=t10-t25;

	/*...Block 3: t5,13,21,29 */

		jt = j1 + p4;
		jp = j2 + p4;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		rt =t13;t13=t5 +t14;t5 =t5 -t14;		/* twiddle mpy by E^4 = I */
				t14=t6 -rt;	t6 =t6 +rt;

		rt =(t21-t22)*ISRT2;t22=(t21+t22)*ISRT2;	t21=rt;	/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t21,22,29,30 */
		rt =(t30+t29)*ISRT2;it =(t30-t29)*ISRT2;		/* twiddle mpy by -E^6 is here... */
		t29=t21+rt;			t21=t21-rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here. */
		t30=t22+it;			t22=t22-it;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		a[jt    ]=t5+t21;	a[jp    ]=t6+t22;
		a[jt+p1 ]=t5-t21;	a[jp+p1 ]=t6-t22;

		a[jt+p2 ]=t13-t30;	a[jp+p2 ]=t14+t29;	/* mpy by E^4=i is inlined here... */
		a[jt+p3 ]=t13+t30;	a[jp+p3 ]=t14-t29;

	/*...Block 2: t3,11,19,27 */

		jt = j1 + p8;
		jp = j2 + p8;

	#ifdef PFETCH_AGGRESSIVE
		addr = addp + p12;
		prefetch_p_doubles(addr);
	#endif
		rt =(t11-t12)*ISRT2;it =(t11+t12)*ISRT2;		/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t11,12 */
		t11=t3 -rt;			t3 =t3 +rt;
		t12=t4 -it;			t4 =t4 +it;

		rt =t19*c - t20*s;	t20=t20*c + t19*s;	t19=rt;	/* twiddle mpy by E^1 */
		rt =t27*s - t28*c;	it =t28*s + t27*c;		/* twiddle mpy by E^3 */
		t27=t19-rt;			t19=t19+rt;
		t28=t20-it;			t20=t20+it;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		a[jt    ]=t3+t19;	a[jp    ]=t4+t20;
		a[jt+p1 ]=t3-t19;	a[jp+p1 ]=t4-t20;

		a[jt+p2 ]=t11-t28;	a[jp+p2 ]=t12+t27;	/* mpy by E^4=i is inlined here... */
		a[jt+p3 ]=t11+t28;	a[jp+p3 ]=t12-t27;

	/*...Block 4: t7,15,23,31 */

		jt = j1 + p12;
		jp = j2 + p12;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		rt =(t16+t15)*ISRT2;it =(t16-t15)*ISRT2;		/* twiddle mpy by -E^6 is here... */
		t15=t7 +rt;			t7 =t7 -rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here. */
		t16=t8 +it;			t8 =t8 -it;

		rt =t23*s - t24*c;	t24=t24*s + t23*c;	t23=rt;	/* twiddle mpy by E^3 */
		rt =t31*c - t32*s;	it =t32*c + t31*s;		/* twiddle mpy by E^1 = -E^9... */
		t31=t23+rt;			t23=t23-rt;			/* ...and get E^9 by flipping signs here. */
		t32=t24+it;			t24=t24-it;			/* Note: t23+rt = t23*(s+1) */

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		a[jt    ]=t7+t23;	a[jp    ]=t8+t24;
		a[jt+p1 ]=t7-t23;	a[jp+p1 ]=t8-t24;

		a[jt+p2 ]=t15-t32;	a[jp+p2 ]=t16+t31;	/* mpy by E^4=i is inlined here... */
		a[jt+p3 ]=t15+t32;	a[jp+p3 ]=t16-t31;

#endif	/* USE_SSE2 */

	  }	/* endfor(j=jlo; j < jhi; j += 4) */

	  /*jlo=jlo+incr; jhi=jhi+incr; See my note about OpenMP above. */

	}	/* endfor(m=0; m < nloops; m++) */

}

/***************/

void radix16_dit_pass(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr)
{
/*
!...Acronym: DIT = Decimation In Time
!...Post-twiddles implementation of radix16_dit_pass.
*/
	static double c = 0.9238795325112867561281831, s = 0.3826834323650897717284599;	/* exp[i*(twopi/16)]*/
	int i,j,j1,j2,jt,jp,jlo,jhi,m,iroot_prim,iroot,k1,k2;
	static int p1,p2,p3,p4,p8,p12;
	double rt,it;
	double re0,im0,re1,im1;
	double *addr, *addp;
	int prefetch_offset;

#ifdef USE_SSE2

	static int	first_entry = TRUE;
	static double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
	static struct complex *sc_arr = 0x0, *sc_ptr;
	static struct complex *cc0, *ss0, *isrt2, *two;
	static struct complex *c0,*c1,*c2,*c3,*c4,*c5,*c6,*c7,*c8,*c9,*c10,*c11,*c12,*c13,*c14,*c15,*s0,*s1,*s2,*s3,*s4,*s5,*s6,*s7,*s8,*s9,*s10,*s11,*s12,*s13,*s14,*s15;
	static struct complex *r1,*r2,*r3,*r4,*r5,*r6,*r7,*r8,*r9,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r30,*r31,*r32;

#else

	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32
	,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15;

#endif

#ifdef USE_SSE2

	if(first_entry)
	{
		first_entry = FALSE;

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
		r6  = sc_ptr + 0x05;		c4  = sc_ptr + 0x25;
		r7  = sc_ptr + 0x06;		s4  = sc_ptr + 0x26;
		r8  = sc_ptr + 0x07;		c8  = sc_ptr + 0x27;
		r9  = sc_ptr + 0x08;		s8  = sc_ptr + 0x28;
		r10 = sc_ptr + 0x09;		c12 = sc_ptr + 0x29;
		r11 = sc_ptr + 0x0a;		s12 = sc_ptr + 0x2a;
		r12 = sc_ptr + 0x0b;		c2  = sc_ptr + 0x2b;
		r13 = sc_ptr + 0x0c;		s2  = sc_ptr + 0x2c;
		r14 = sc_ptr + 0x0d;		c6  = sc_ptr + 0x2d;
		r15 = sc_ptr + 0x0e;		s6  = sc_ptr + 0x2e;
		r16 = sc_ptr + 0x0f;		c10 = sc_ptr + 0x2f;
		r17 = sc_ptr + 0x10;		s10 = sc_ptr + 0x30;
		r18 = sc_ptr + 0x11;		c14 = sc_ptr + 0x31;
		r19 = sc_ptr + 0x12;		s14 = sc_ptr + 0x32;
		r20 = sc_ptr + 0x13;		c1  = sc_ptr + 0x33;
		r21 = sc_ptr + 0x14;		s1  = sc_ptr + 0x34;
		r22 = sc_ptr + 0x15;		c5  = sc_ptr + 0x35;
		r23 = sc_ptr + 0x16;		s5  = sc_ptr + 0x36;
		r24 = sc_ptr + 0x17;		c9  = sc_ptr + 0x37;
		r25 = sc_ptr + 0x18;		s9  = sc_ptr + 0x38;
		r26 = sc_ptr + 0x19;		c13 = sc_ptr + 0x39;
		r27 = sc_ptr + 0x1a;		s13 = sc_ptr + 0x3a;
		r28 = sc_ptr + 0x1b;		c3  = sc_ptr + 0x3b;
		r29 = sc_ptr + 0x1c;		s3  = sc_ptr + 0x3c;
		r30 = sc_ptr + 0x1d;		c7  = sc_ptr + 0x3d;
		r31 = sc_ptr + 0x1e;		s7  = sc_ptr + 0x3e;
		r32 = sc_ptr + 0x1f;		c11 = sc_ptr + 0x3f;
									s11 = sc_ptr + 0x40;
									c15 = sc_ptr + 0x41;
									s15 = sc_ptr + 0x42;
									two = sc_ptr + 0x43;
		/* These remain fixed: */
		isrt2->re = ISRT2;	isrt2->im = ISRT2;
		two  ->re = 2.0;	two  ->im = 2.0;
		cc0  ->re = c	;	cc0  ->im = c	;
		ss0  ->re = s	;	ss0  ->im = s	;

	}	/* end of inits */

#endif

	p1 = incr >> 4;
	p2 = p1 +p1;
	p3 = p2 +p1;
	p4 = p3 +p1;
	p8 = p4 +p4;
	p12= p8 +p4;

	p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
	p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
	p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
	p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
	p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
	p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );
	ASSERT(HERE, p8  == p4+p4, "radix16_dit_pass: p8  != p4+p4!");
	ASSERT(HERE, p12 == p4+p8, "radix16_dit_pass: p12 != p4+p8!");

	iroot_prim=(incr >> 5);		/* (incr/2)/radix_now */

	for(m=0; m < nloops; m++)	/* NLOOPS may range from 1 (if first pass radix = 16) to P*N/32 (last pass radix = 16).	 */
	{				/* NLOOPS satisfies the identity NLOOPS * INCR = P*N, which says that we process the entire */
					/* array each time this subroutine is executed (since P*N = vector length, sans padding.)   */

/*	here are the needed sincos data - these are processed below in bit-reversed order. */
	  iroot = index[m]*iroot_prim;
	  i = iroot;

/* In the post-twiddle-mul version of the DIT pass we may be able to afford a more-accurate twiddles computation,
   since we only compute each set of twiddles once and re-use it for many data blocks.
*/
#if HIACC
	#ifdef USE_SSE2
		c0->re=1.;	s0->re=0.;
		c0->im=1.;	s0->im=0.;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 2*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c1 ->re=rt;	s1 ->re=it;
		c1 ->im=rt;	s1 ->im=it;
	#else
		c1 =rt;		s1 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 3*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c2 ->re=rt;	s2 ->re=it;
		c2 ->im=rt;	s2 ->im=it;
	#else
		c2 =rt;		s2 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 4*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c3 ->re=rt;	s3 ->re=it;
		c3 ->im=rt;	s3 ->im=it;
	#else
		c3 =rt;		s3 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 5*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c4 ->re=rt;	s4 ->re=it;
		c4 ->im=rt;	s4 ->im=it;
	#else
		c4 =rt;		s4 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 6*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c5 ->re=rt;	s5 ->re=it;
		c5 ->im=rt;	s5 ->im=it;
	#else
		c5 =rt;		s5 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 7*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c6 ->re=rt;	s6 ->re=it;
		c6 ->im=rt;	s6 ->im=it;
	#else
		c6 =rt;		s6 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 8*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c7 ->re=rt;	s7 ->re=it;
		c7 ->im=rt;	s7 ->im=it;
	#else
		c7 =rt;		s7 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 9*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c8 ->re=rt;	s8 ->re=it;
		c8 ->im=rt;	s8 ->im=it;
	#else
		c8 =rt;		s8 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*10*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c9 ->re=rt;	s9 ->re=it;
		c9 ->im=rt;	s9 ->im=it;
	#else
		c9 =rt;		s9 =it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*11*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c10->re=rt;	s10->re=it;
		c10->im=rt;	s10->im=it;
	#else
		c10=rt;		s10=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*12*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c11->re=rt;	s11->re=it;
		c11->im=rt;	s11->im=it;
	#else
		c11=rt;		s11=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*13*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c12->re=rt;	s12->re=it;
		c12->im=rt;	s12->im=it;
	#else
		c12=rt;		s12=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*14*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c13->re=rt;	s13->re=it;
		c13->im=rt;	s13->im=it;
	#else
		c13=rt;		s13=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/*15*iroot */
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c14->re=rt;	s14->re=it;
		c14->im=rt;	s14->im=it;
	#else
		c14=rt;		s14=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c15->re=rt;	s15->re=it;
		c15->im=rt;	s15->im=it;
	#else
		c15=rt;		s15=it;
	#endif

#else
		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;			/* 2*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c1=t1*t3-t2*t4;	s1=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += i;				/* 4*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c2=t1*t3-t2*t4;	s2=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += i;				/* 8*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c4=t1*t3-t2*t4;	s4=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += (iroot << 2)+iroot;		/* 13*iroot */
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c8=t1*t3-t2*t4;	s8=t1*t4+t2*t3;

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		t1=rt0[k1].re;	t2=rt0[k1].im;
		t3=rt1[k2].re;	t4=rt1[k2].im;
		c13=t1*t3-t2*t4;	s13=t1*t4+t2*t3;

		/* c3,5 */
		t1=c1*c4; t2=c1*s4; t3=s1*c4; t4=s1*s4;
		c3=t1+t4; s3=t2-t3; c5=t1-t4; s5=t2+t3;

		/* c6,7,9,10 */
		t1=c1*c8; t2=c1*s8; t3=s1*c8; t4=s1*s8;
		c7=t1+t4; s7=t2-t3; c9=t1-t4; s9=t2+t3;

		t1=c2*c8; t2=c2*s8; t3=s2*c8; t4=s2*s8;
		c6=t1+t4; s6=t2-t3; c10=t1-t4; s10=t2+t3;

		/* c11,12,14,15 */
		t1=c1*c13; t2=c1*s13; t3=s1*c13; t4=s1*s13;
		c12=t1+t4; s12=t2-t3; c14=t1-t4; s14=t2+t3;

		t1=c2*c13; t2=c2*s13; t3=s2*c13; t4=s2*s13;
		c11=t1+t4; s11=t2-t3; c15=t1-t4; s15=t2+t3;
#endif

	/* Define the inner-loop parameters in terms of the outer-loop ones to make OpenMP's job easier: */
	  jlo = m*incr;
	  jhi = jlo+(incr >> 4);

	#ifdef USE_SSE2
	  for(j=jlo; j < jhi; j += 4)
	  {
		/* In SSE2 mode, data are arranged in [re0,re1,im0,im1] quartets, not the usual [re0,im0],[re1,im1] pairs.
		Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
		but to point to the proper real datum, we need to bit-reverse bits <0:1> of j, i.e. [0,1,2,3] ==> [0,2,1,3].
		*/
		j1 = (j & mask01) + br4[j&3];
	#elif defined(USE_SSE2)	/* This allows us to use #if 0 above and disable sse2-based *computation*, while still using sse2-style data layout */
	  for(j=jlo; j < jhi; j += 2)
	  {
		j1 = (j & mask01) + br4[j&3];
	#else
	  for(j=jlo; j < jhi; j += 2)
	  {
		j1 =  j;
	#endif
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;

#ifdef USE_SSE2

	#ifdef DEBUG_SSE2
		rng_isaac_init(TRUE);
		jt = j1;		jp = j2;
		for(i = 0; i < 16; i++)
		{
			a[jt] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			fprintf(stderr, "radix16_dit_pass: A_in[%2d] = %20.5f, %20.5f\n",i,a[jt],a[jp]);
			jt += p1;	jp += p1;
		}
	#endif

	#ifdef COMPILER_TYPE_MSVC

	/*...Block 1: */
	#if 1
		add0 = &a[j1];
		__asm	mov eax, add0
		__asm	mov ebx, p1
		__asm	mov ecx, p2
		__asm	mov edx, p3
		__asm	mov	edi, p4		/* edi will store copy of p4 throughout */
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	shl	ecx, 3
		__asm	shl	edx, 3
		__asm	shl	edi, 3
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
		SSE2_RADIX4_DIT_0TWIDDLE_B(r1)
	#else
		add0 = &a[j1];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r1)
	#endif

	/*...Block 2: */
	#if 1
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
		SSE2_RADIX4_DIT_0TWIDDLE_B(r9)
	#else
		add0 = &a[j1+p4];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r9)
	#endif

	/*...Block 3: */
	#if 1
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
		SSE2_RADIX4_DIT_0TWIDDLE_B(r17)
	#else
		add0 = &a[j1+p8];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r17)
	#endif

	/*...Block 4: */
	#if 1
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
		SSE2_RADIX4_DIT_0TWIDDLE_B(r25)
	#else
		add0 = &a[j1+p12];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r25)
	#endif

	/****************************************************************************************************
	!...and now do four more radix-4 transforms, including the internal and external twiddle factors:   !
	****************************************************************************************************/

	/*...Block 1: Cost: 24 MOVapd, 26 ADD/SUBpd, 12 MULpd */
		__asm	mov eax, add0
		__asm	mov ebx, p4
		__asm	mov ecx, r1
		__asm	mov edx, r9
		__asm	mov edi, p8		/* edi will store copy of p8 throughout */
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	shl	edi, 3
		__asm	add ebx, eax	/* add1 = add0+p4 */

		__asm	movaps	xmm2,[edx      ]	/* t10 */				__asm	movaps	xmm4,[edx+0x100]	/* t30 */
		__asm	movaps	xmm3,[edx+0x010]	/* t11 */				__asm	movaps	xmm5,[edx+0x110]	/* t31 */
		__asm	movaps	xmm0,[ecx      ]	/* t00 */				__asm	movaps	xmm6,[ecx+0x100]	/* t20 */
		__asm	movaps	xmm1,[ecx+0x010]	/* t01 */				__asm	movaps	xmm7,[ecx+0x110]	/* t21 */

		__asm	subpd	xmm0,xmm2			/*~t10=t00-t10*/		__asm	subpd	xmm6,xmm4			/*~t30=t20-t30*/
		__asm	subpd	xmm1,xmm3			/*~t11=t01-t11*/		__asm	subpd	xmm7,xmm5			/*~t31=t21-t31*/
		__asm	addpd	xmm2,xmm2			/*       2*t10*/		__asm	addpd	xmm4,xmm4			/*       2*t30*/
		__asm	addpd	xmm3,xmm3			/*       2*t11*/		__asm	addpd	xmm5,xmm5			/*       2*t31*/
		__asm	addpd	xmm2,xmm0			/*~t00=t00+t10*/		__asm	addpd	xmm4,xmm6			/*~t20=t20+t30*/
		__asm	addpd	xmm3,xmm1			/*~t01=t01+t11*/		__asm	addpd	xmm5,xmm7			/*~t21=t21+t31*/

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c0)

	/*...Block 2: Cost: 34 MOVapd, 36 ADD/SUBpd, 26 MULpd */
		__asm	mov eax, add0
		__asm	mov esi, p1
		__asm	mov ebx, p4
		__asm	shl	esi, 3
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	add eax, esi	/* add0 = &a[j1+p1] */
		__asm	add ebx, eax	/* add1 = add0+p4 */
		__asm	mov ecx, r3
		__asm	mov edx, r11
		__asm	add	ecx, 0x100
		__asm	add	edx, 0x100
		__asm	mov	esi, cc0

		__asm	movaps	xmm4,[ecx      ]	/* t24 */
		__asm	movaps	xmm5,[ecx+0x010]	/* t25 */
		__asm	movaps	xmm2,[esi      ]	/* c */
		__asm	movaps	xmm3,[esi+0x010]	/* s */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy t24 */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy t25 */

		__asm	mulpd	xmm4,xmm2		/* t24*c */
		__asm	mulpd	xmm5,xmm2		/* t25*c */
		__asm	mulpd	xmm6,xmm3		/* t24*s */		__asm	movaps	xmm0,[edx      ]	/* t34 */
		__asm	mulpd	xmm7,xmm3		/* t25*s */		__asm	movaps	xmm1,[edx+0x010]	/* t35 */
		__asm	subpd	xmm5,xmm6	/* xmm1 <-~t25 */	__asm	movaps	xmm6,xmm0	/* xmm6 <- cpy t34 */
		__asm	addpd	xmm4,xmm7	/* xmm0 <-~t24 */	__asm	movaps	xmm7,xmm1	/* xmm7 <- cpy t35 */

														__asm	mulpd	xmm0,xmm3		/* t34*s */
														__asm	mulpd	xmm1,xmm3		/* t35*s */
														__asm	mulpd	xmm6,xmm2		/* t34*c */
														__asm	mulpd	xmm7,xmm2		/* t35*c */
														__asm	subpd	xmm1,xmm6	/* xmm5 <- it */
														__asm	addpd	xmm0,xmm7	/* xmm4 <- rt */

		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy~t25*/
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy~t24*/

		__asm	addpd	xmm4,xmm0	/* ~t24 <- t24+rt */
		__asm	addpd	xmm5,xmm1	/* ~t25 <- t25+it */
		__asm	subpd	xmm6,xmm0	/* ~t34 <- t24-rt */
		__asm	subpd	xmm7,xmm1	/* ~t35 <- t25-it */

		__asm	sub	ecx, 0x100
		__asm	sub	edx, 0x100
		__asm	mov	esi, isrt2
		__asm	movaps	xmm2,[edx      ]	/* t14 */
		__asm	movaps	xmm3,[edx+0x010]	/* t15 */
		__asm	movaps	xmm1,[esi]	/* isrt2 */
		__asm	movaps	xmm0,xmm3	/* cpy t15 */
		__asm	subpd	xmm3,xmm2	/*~t15=t15-t14 */
		__asm	addpd	xmm2,xmm0	/*~t14=t14+t15 */
		__asm	mulpd	xmm2,xmm1	/* rt */
		__asm	mulpd	xmm3,xmm1	/* it */

		__asm	movaps	xmm0,[ecx      ]	/* t04 */
		__asm	movaps	xmm1,[ecx+0x010]	/* t05 */
		__asm	subpd	xmm0,xmm2	/*~t14 <- t04- rt */
		__asm	subpd	xmm1,xmm3	/*~t15 <- t05- it */
		__asm	addpd	xmm2,xmm2	/*          2* rt */
		__asm	addpd	xmm3,xmm3	/*          2* it */
		__asm	addpd	xmm2,xmm0	/*~t04 <- t04+ rt */
		__asm	addpd	xmm3,xmm1	/*~t05 <- t05+ it */

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c1)

	/*...Block 3: Cost: 31 MOVapd, 32 ADD/SUBpd, 20 MULpd */
		__asm	mov eax, add0
		__asm	mov esi, p2
		__asm	mov ebx, p4
		__asm	shl	esi, 3
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	add eax, esi	/* add0 = &a[j1+p2] */
		__asm	add ebx, eax	/* add1 = add0+p4 */
		__asm	mov ecx, r5
		__asm	mov edx, r13
		__asm	add	ecx, 0x100
		__asm	add	edx, 0x100
		__asm	mov	esi, isrt2
		__asm	movaps	xmm2,[esi]	/* isrt2 */

		__asm	movaps	xmm4,[ecx      ]	/* t28 */
		__asm	movaps	xmm5,[ecx+0x010]	/* t29 */
		__asm	movaps	xmm6,[edx      ]	/* t38 */
		__asm	movaps	xmm7,[edx+0x010]	/* t39 */
		__asm	sub	ecx, 0x100
		__asm	sub	edx, 0x100
		__asm	mulpd	xmm4,xmm2
		__asm	mulpd	xmm5,xmm2
		__asm	mulpd	xmm6,xmm2
		__asm	mulpd	xmm7,xmm2

		__asm	subpd	xmm5,xmm4			/*~t29=t29-t28*/		__asm	movaps	xmm0,[ecx      ]	/* t08 */
		__asm	subpd	xmm6,xmm7			/* rt =t38-t39*/		__asm	movaps	xmm2,[edx+0x010]	/* t19 */
		__asm	addpd	xmm4,xmm4			/*       2*t28*/		__asm	movaps	xmm3,[ecx+0x010]	/* t09 */
		__asm	addpd	xmm7,xmm7			/*       2*t39*/		__asm	movaps	xmm1,[edx      ]	/* t18 */
		__asm	addpd	xmm4,xmm5			/*~t28=t28+t29*/
		__asm	addpd	xmm7,xmm6			/* it =t39+t38*/

		__asm	subpd	xmm4,xmm6			/*~t28=t28-rt */		__asm	subpd	xmm0,xmm2			/*~t18=t08-t19*/
		__asm	subpd	xmm5,xmm7			/*~t29=t29-it */		__asm	subpd	xmm3,xmm1			/*~t09=t09-t18*/
		__asm	addpd	xmm6,xmm6			/*       2*rt */		__asm	addpd	xmm2,xmm2			/*       2*t08*/
		__asm	addpd	xmm7,xmm7			/*       2*it */		__asm	addpd	xmm1,xmm1			/*       2*t09*/
		__asm	addpd	xmm6,xmm4			/*~t38=t28+rt */		__asm	addpd	xmm2,xmm0			/*~t08=t19+t08*/
		__asm	addpd	xmm7,xmm5			/*~t39=t29+it */		__asm	addpd	xmm1,xmm3			/*~t19=t18+t09*/

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c2)

	/*...Block 4: Cost: 34 MOVapd, 36 ADD/SUBpd, 26 MULpd */
	#if 0
		add0 = &a[j1+p3];
		add1 = add0+p4;
		add2 = add0+p8;
		add3 = add0+p12;

		__asm	mov	eax, r7
		__asm	mov	edx, eax
		__asm	add	eax, 0x100
		__asm	add	edx, 0x180
		__asm	mov	ecx, cc0

		__asm	movaps	xmm4,[eax      ]	/* t2C */
		__asm	movaps	xmm5,[eax+0x010]	/* t2D */
		__asm	movaps	xmm3,[ecx      ]	/* c */
		__asm	movaps	xmm2,[ecx+0x010]	/* s */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy t2C */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy t2D */

		__asm	mulpd	xmm4,xmm2		/* t2C*s */
		__asm	mulpd	xmm5,xmm2		/* t2D*s */
		__asm	mulpd	xmm6,xmm3		/* t2C*c */		__asm	movaps	xmm0,[edx      ]	/* t3C */
		__asm	mulpd	xmm7,xmm3		/* t2D*c */		__asm	movaps	xmm1,[edx+0x010]	/* t3D */
		__asm	subpd	xmm5,xmm6	/* xmm5 <-~t2D */	__asm	movaps	xmm6,xmm0	/* xmm6 <- cpy t3C */
		__asm	addpd	xmm4,xmm7	/* xmm4 <-~t2C */	__asm	movaps	xmm7,xmm1	/* xmm7 <- cpy t3D */

														__asm	mulpd	xmm0,xmm3		/* t3C*c */
														__asm	mulpd	xmm1,xmm3		/* t3D*c */
														__asm	mulpd	xmm6,xmm2		/* t3C*s */
														__asm	mulpd	xmm7,xmm2		/* t3D*s */
														__asm	subpd	xmm1,xmm6	/* xmm1 <- it */
														__asm	addpd	xmm0,xmm7	/* xmm0 <- rt */

		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy~t2D*/
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy~t2C*/

		__asm	addpd	xmm6,xmm0	/* ~t3C <- t2C+rt */
		__asm	addpd	xmm7,xmm1	/* ~t3D <- t2D+it */
		__asm	subpd	xmm4,xmm0	/* ~t2C <- t2C-rt */
		__asm	subpd	xmm5,xmm1	/* ~t2D <- t2D-it */

		__asm	sub	eax, 0x100
		__asm	sub	edx, 0x100
		__asm	mov	ecx, isrt2
		__asm	movaps	xmm0,[edx      ]	/* t1C */
		__asm	movaps	xmm1,[edx+0x010]	/* t1D */
		__asm	movaps	xmm3,[ecx]	/* isrt2 */
		__asm	movaps	xmm2,xmm0	/* cpy t1C */
		__asm	subpd	xmm0,xmm1	/*~t1C=t1C-t1D */
		__asm	addpd	xmm1,xmm2	/*~t1D=t1D+t1C */
		__asm	mulpd	xmm0,xmm3	/* it */
		__asm	mulpd	xmm1,xmm3	/* rt */

		__asm	movaps	xmm2,[eax      ]	/* t0C */
		__asm	movaps	xmm3,[eax+0x010]	/* t0D */
		__asm	subpd	xmm2,xmm0	/*~t0C <- t0C- rt */
		__asm	subpd	xmm3,xmm1	/*~t0D <- t0D- it */
		__asm	addpd	xmm0,xmm0	/*          2* rt */
		__asm	addpd	xmm1,xmm1	/*          2* it */
		__asm	addpd	xmm0,xmm2	/*~t1C <- t0C+ rt */
		__asm	addpd	xmm1,xmm3	/*~t1D <- t0D+ it */

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF(add0, add1, add2, add3, c3)
	#else
		__asm	mov eax, add0
		__asm	mov esi, p3
		__asm	mov ebx, p4
		__asm	shl	esi, 3
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	add eax, esi	/* add0 = &a[j1+p2] */
		__asm	add ebx, eax	/* add1 = add0+p4 */
		__asm	mov ecx, r7
		__asm	mov edx, r15
		__asm	add	ecx, 0x100
		__asm	add	edx, 0x100
		__asm	mov	esi, cc0

		__asm	movaps	xmm4,[ecx      ]	/* t2C */
		__asm	movaps	xmm5,[ecx+0x010]	/* t2D */
		__asm	movaps	xmm3,[esi      ]	/* c */
		__asm	movaps	xmm2,[esi+0x010]	/* s */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy t2C */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy t2D */

		__asm	mulpd	xmm4,xmm2		/* t2C*s */
		__asm	mulpd	xmm5,xmm2		/* t2D*s */
		__asm	mulpd	xmm6,xmm3		/* t2C*c */		__asm	movaps	xmm0,[edx      ]	/* t3C */
		__asm	mulpd	xmm7,xmm3		/* t2D*c */		__asm	movaps	xmm1,[edx+0x010]	/* t3D */
		__asm	subpd	xmm5,xmm6	/* xmm5 <-~t2D */	__asm	movaps	xmm6,xmm0	/* xmm6 <- cpy t3C */
		__asm	addpd	xmm4,xmm7	/* xmm4 <-~t2C */	__asm	movaps	xmm7,xmm1	/* xmm7 <- cpy t3D */

														__asm	mulpd	xmm0,xmm3		/* t3C*c */
														__asm	mulpd	xmm1,xmm3		/* t3D*c */
														__asm	mulpd	xmm6,xmm2		/* t3C*s */
														__asm	mulpd	xmm7,xmm2		/* t3D*s */
														__asm	subpd	xmm1,xmm6	/* xmm1 <- it */
														__asm	addpd	xmm0,xmm7	/* xmm0 <- rt */

		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy~t2D*/
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy~t2C*/

		__asm	addpd	xmm6,xmm0	/* ~t3C <- t2C+rt */
		__asm	addpd	xmm7,xmm1	/* ~t3D <- t2D+it */
		__asm	subpd	xmm4,xmm0	/* ~t2C <- t2C-rt */
		__asm	subpd	xmm5,xmm1	/* ~t2D <- t2D-it */

		__asm	sub	ecx, 0x100
		__asm	sub	edx, 0x100
		__asm	mov	esi, isrt2
		__asm	movaps	xmm0,[edx      ]	/* t1C */
		__asm	movaps	xmm1,[edx+0x010]	/* t1D */
		__asm	movaps	xmm3,[esi]	/* isrt2 */
		__asm	movaps	xmm2,xmm0	/* cpy t1C */
		__asm	subpd	xmm0,xmm1	/*~t1C=t1C-t1D */
		__asm	addpd	xmm1,xmm2	/*~t1D=t1D+t1C */
		__asm	mulpd	xmm0,xmm3	/* it */
		__asm	mulpd	xmm1,xmm3	/* rt */

		__asm	movaps	xmm2,[ecx      ]	/* t0C */
		__asm	movaps	xmm3,[ecx+0x010]	/* t0D */
		__asm	subpd	xmm2,xmm0	/*~t0C <- t0C- rt */
		__asm	subpd	xmm3,xmm1	/*~t0D <- t0D- it */
		__asm	addpd	xmm0,xmm0	/*          2* rt */
		__asm	addpd	xmm1,xmm1	/*          2* it */
		__asm	addpd	xmm0,xmm2	/*~t1C <- t0C+ rt */
		__asm	addpd	xmm1,xmm3	/*~t1D <- t0D+ it */

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c3)
	#endif

		/***************************************************/
		/* Total Cost: 187 MOVapd, 210 ADD/SUBpd, 84 MULpd */
		/***************************************************/

	#else	/* GCC-style inline ASM: */

		add0 = &a[j1];
		SSE2_RADIX16_DIT_TWIDDLE(add0,p1,p2,p3,p4,p8,r1,r3,r5,r7,r9,r11,r13,r15,r17,r25,isrt2,cc0,c0,c1,c2,c3);

	#endif

	#ifdef DEBUG_SSE2
		jt = j1;		jp = j2;
		for(i = 0; i < 16; i++)
		{
			fprintf(stderr, "radix16_dit_pass: Aout[%2d] = %20.5f, %20.5f\n",i,a[jt],a[jp]);
			jt += p1;	jp += p1;
		}
		exit(0);
	#endif

#else	/* USE_SSE2 */

	/* gather the needed data (8 64-bit complex, i.e. 16 64-bit reals) and do the first set of four length-4 transforms.
	   We process the sincos data in bit-reversed order.	*/
	#ifdef PFETCH_AGGRESSIVE
		addp = &a[j1];
	#elif PFETCH
		prefetch_offset = ((j >> 1) & 3)*p1 + 4;	/* Cycle among p0, p1, p2 and p3. */
	#endif
	/*...Block 1: */

		jt = j1;
		jp = j2;

	#ifdef PFETCH_AGGRESSIVE
		addr = addp;
		prefetch_p_doubles(addr);
	#elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	#endif
		t1 =a[jt ];	t2 =a[jp ];
		rt =a[jt+p1 ];	it =a[jp+p1 ];
		t3 =t1 -rt;		t1 =t1 +rt;
		t4 =t2 -it;		t2 =t2 +it;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		rt =a[jt+p2 ];	t6 =a[jp+p2 ];	t5 =rt;
		rt =a[jt+p3 ];	it =a[jp+p3 ];
		t7 =t5 -rt;		t5 =t5 +rt;
		t8 =t6 -it;		t6 =t6 +it;

		rt =t5;	t5 =t1 -rt ;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it ;	t2 =t2 +it;

		rt =t7;	t7 =t3 -t8 ;	t3 =t3 +t8;
				t8 =t4 +rt ;	t4 =t4 -rt;

	/*...Block 2: */

		jt = j1 + p4;
		jp = j2 + p4;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	#endif
		rt =a[jt    ];	t10=a[jp    ];	t9 =rt;
		rt =a[jt+p1 ];	it =a[jp+p1 ];
		t11=t9 -rt;		t9 =t9 +rt;
		t12=t10-it;		t10=t10+it;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		rt =a[jt+p2 ];	t14=a[jp+p2 ];	t13=rt;
		rt =a[jt+p3 ];	it =a[jp+p3 ];
		t15=t13-rt;		t13=t13+rt;
		t16=t14-it;		t14=t14+it;

		rt =t13;	t13=t9 -rt ;	t9 =t9 +rt;
		it =t14;	t14=t10-it ;	t10=t10+it;

		rt =t15;	t15=t11-t16;	t11=t11+t16;
					t16=t12+rt ;	t12=t12-rt;

	/*...Block 3: */

		jt = j1 + p8;
		jp = j2 + p8;

	#ifdef PFETCH_AGGRESSIVE
		addr = addp + p4 ;
		prefetch_p_doubles(addr);
	#elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	#endif
		rt =a[jt    ];	t18=a[jp    ];	t17=rt;
		rt =a[jt+p1 ];	it =a[jp+p1 ];
		t19=t17-rt;		t17=t17+rt;
		t20=t18-it;		t18=t18+it;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		rt =a[jt+p2 ];	t22=a[jp+p2 ];	t21=rt;
		rt =a[jt+p3 ];	it =a[jp+p3 ];
		t23=t21-rt;		t21=t21+rt;
		t24=t22-it;		t22=t22+it;

		rt =t21;	t21=t17-rt ;	t17=t17+rt;
		it =t22;	t22=t18-it ;	t18=t18+it;

		rt =t23;	t23=t19-t24;	t19=t19+t24;
					t24=t20+rt ;	t20=t20-rt;

	/*...Block 4: */

		jt = j1 + p12;
		jp = j2 + p12;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#elif PFETCH
		addp = &a[jt];
		addr = addp+prefetch_offset;
		prefetch_p_doubles(addr);
	#endif
		rt =a[jt    ];	t26=a[jp    ];	t25=rt;
		rt =a[jt+p1 ];	it =a[jp+p1 ];
		t27=t25-rt;		t25=t25+rt;
		t28=t26-it;		t26=t26+it;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		rt =a[jt+p2 ];	t30=a[jp+p2 ];	t29=rt;
		rt =a[jt+p3 ];	it =a[jp+p3 ];
		t31=t29-rt;		t29=t29+rt;
		t32=t30-it;		t30=t30+it;

		rt =t29;	t29=t25-rt ;	t25=t25+rt;
		it =t30;	t30=t26-it ;	t26=t26+it;

		rt =t31;	t31=t27-t32;	t27=t27+t32;
					t32=t28+rt ;	t28=t28-rt;
	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	!	1, exp(-i* 1*twopi/16) =	   ( c,-s), exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 3*twopi/16) =	   ( s,-c) (for inputs to transform block 2)
	!	1, exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 4*twopi/16) =	   ( 0,-1), exp(-i* 6*twopi/16) = ISRT2*(-1,-1) (for inputs to transform block 3)
	!	1, exp(-i* 3*twopi/16) =	   ( s,-c), exp(-i* 6*twopi/16) = ISRT2*(-1,-1), exp(-i* 9*twopi/16) =	   (-c, s) (for inputs to transform block 4).
	!  (This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
	!   I.e. do similar as above, except inputs a[jt  +p0:15:1) are replaced by t0:30:2,
	!										   a[jp+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
	!   and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.
	*/
	/*...Block 1: t1,9,17,25 */

		jt = j1;
		jp = j2;

	#ifdef PFETCH_AGGRESSIVE
		addr = addp + p8 ;
		prefetch_p_doubles(addr);
	#endif
		rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		a[jt    ]=t1+t17;			a[jp    ]=t2+t18;
		t1	     =t1-t17;			t2		 =t2-t18;
		a[jt+p8 ]=t1 *c8 +t2 *s8;	a[jp+p8 ]=t2 *c8 -t1 *s8;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		rt	   =t9 +t26;		it		 =t10-t25;	/* mpy by E^-4 = -I is inlined here... */
		t9	   =t9 -t26;		t10		=t10+t25;
		a[jt+p4 ]=rt *c4 +it *s4;	a[jp+p4 ]=it *c4 -rt *s4;
		a[jt+p12]=t9 *c12+t10*s12;	a[jp+p12]=t10*c12-t9 *s12;

	/*...Block 3: t5,13,21,29 */

		jt = j1 + p2;
		jp = j2 + p2;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		rt =t13;	t13=t5 -t14;	t5 =t5 +t14;		/* twiddle mpy by E^4 =-I */
			t14=t6 +rt;	t6 =t6 -rt;

		rt =(t22+t21)*ISRT2;t22=(t22-t21)*ISRT2;	t21=rt;	/* twiddle mpy by E^-2 */
		rt =(t29-t30)*ISRT2;it =(t29+t30)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here... */
		t29=t21+rt;		t21=t21-rt;			/* ...and get E^-6 by flipping signs here. */
		t30=t22+it;		t22=t22-it;

		rt	   =t5 +t21;		it		 =t6 +t22;
		t5	   =t5 -t21;		t6		 =t6 -t22;
		a[jt    ]=rt *c2 +it *s2;	a[jp    ]=it *c2 -rt *s2;
		a[jt+p8 ]=t5 *c10+t6 *s10;	a[jp+p8 ]=t6 *c10-t5 *s10;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		rt	   =t13+t30;		it		 =t14-t29;	/* mpy by E^-4 = -I is inlined here... */
		t13	  =t13-t30;		t14		=t14+t29;
		a[jt+p4 ]=rt *c6 +it *s6;	a[jp+p4 ]=it *c6 -rt *s6;
		a[jt+p12]=t13*c14+t14*s14;	a[jp+p12]=t14*c14-t13*s14;

	/*...Block 2: t3,11,19,27 */

		jt = j1 + p1;
		jp = j2 + p1;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		rt =(t12+t11)*ISRT2;it =(t12-t11)*ISRT2;		/* twiddle mpy by E^-2 */
		t11=t3 -rt;		t3 =t3 +rt;
		t12=t4 -it;		t4 =t4 +it;

		rt =t19*c + t20*s;	t20=t20*c - t19*s;	t19=rt;	/* twiddle mpy by E^-1 */
		rt =t27*s + t28*c;	it =t28*s - t27*c;		/* twiddle mpy by E^-3 */
		t27=t19-rt;		t19=t19+rt;
		t28=t20-it;		t20=t20+it;

		rt	   =t3 +t19;		it		 =t4 +t20;
		t3	   =t3 -t19;		t4		 =t4 -t20;
		a[jt    ]=rt *c1 +it *s1;	a[jp    ]=it *c1 -rt *s1;
		a[jt+p8 ]=t3 *c9 +t4 *s9;	a[jp+p8 ]=t4 *c9 -t3 *s9;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		rt	   =t11+t28;		it		 =t12-t27;	/* mpy by E^-4 = -I is inlined here... */
		t11	  =t11-t28;		t12		=t12+t27;
		a[jt+p4 ]=rt *c5 +it *s5;	a[jp+p4 ]=it *c5 -rt *s5;
		a[jt+p12]=t11*c13+t12*s13;	a[jp+p12]=t12*c13-t11*s13;

	/*...Block 4: t7,15,23,31 */

		jt = j1 + p3;
		jp = j2 + p3;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		rt =(t15-t16)*ISRT2;it =(t15+t16)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here... */
		t15=t7 +rt;		t7 =t7 -rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here. */
		t16=t8 +it;		t8 =t8 -it;

		rt =t23*s + t24*c;	t24=t24*s - t23*c;	t23=rt;	/* twiddle mpy by E^-3 */
		rt =t31*c + t32*s;	it =t32*c - t31*s;		/* twiddle mpy by E^-1 = -E^-9... */
		t31=t23+rt;		t23=t23-rt;			/* ...and get E^9 by flipping signs here. */
		t32=t24+it;		t24=t24-it;

		rt	   =t7 +t23;		it		 =t8 +t24;
		t7	   =t7 -t23;		t8		 =t8 -t24;
		a[jt    ]=rt *c3 +it *s3;	a[jp    ]=it *c3 -rt *s3;
		a[jt+p8 ]=t7 *c11+t8 *s11;	a[jp+p8 ]=t8 *c11-t7 *s11;

	#ifdef PFETCH_AGGRESSIVE
		addr += p1;
		prefetch_p_doubles(addr);
	#endif
		rt	   =t15+t32;		it		 =t16-t31;	/* mpy by E^-4 = -I is inlined here... */
		t15	  =t15-t32;		t16		=t16+t31;
		a[jt+p4 ]=rt *c7 +it *s7;	a[jp+p4 ]=it *c7 -rt *s7;
		a[jt+p12]=t15*c15+t16*s15;	a[jp+p12]=t16*c15-t15*s15;
#endif	/* USE_SSE2 */

	  }	/* endfor(j=jlo; j < jhi; j += 4) */

	  /*jlo=jlo+incr; jhi=jhi+incr; See my note about OpenMP above. */

	}	/* endfor(m=0; m < nloops; m++) */

}

