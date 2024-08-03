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

/*******************************************************************************
   We now include this header file if it was not included before.
*******************************************************************************/
#ifndef radix16_dif_dit_pass_h_included
#define radix16_dif_dit_pass_h_included

/*** GCC_BUG: Jul 2013:
Note that the various prefetch-implementing versions of the radix-16 DFTs fail to build under -O0 (i.e. my standard
debug-build settings), giving the following error:

	error: impossible constraint in ‘asm’

This GCC bug is detailed at.http://gcc.gnu.org/bugzilla/show_bug.cgi?id=23200

The workaround is to use -O1 or higher, whether one is building a debuggable binary or not.
***/
#ifdef USE_AVX2	// These non-inline-asm macros shared by avx2 and avx512 builds:

	// See the AVX2 comments in radix16_dif_dit_pass.c for details on the data layout here:
	#define FMA_TWIDDLE_FIDDLE(\
	          __c8,__s8,__c4,__s4,__cC,__sC,\
	__c2,__s2,__cA,__sA,__c6,__s6,__cE,__sE,\
	__c1,__s1,__c9,__s9,__c5,__s5,__cD,__sD,\
	__c3,__s3,__cB,__sB,__c7,__s7,__cF,__sF,\
	__c,__tan,\
	__twid_ptr)\
	{\
		double *add0,*add1,*add2;\
		add0 = (double *)__twid_ptr;	/* add0 points to 16 cos-data-to-be-inverted; Need a double-ptr on lhs here */\
		ASSERT(add0 != 0x0, "Null add0 pointer!");\
		add1 = add0 + 16;	/* add1 points to block of memory temporarily used to store the corresponding sine data */\
		add2 = add0 + 32;	/* add2 points to block of memory temporarily used to store the 11 [0-padded to 12]
							cosine data which need to be divided by other cosines (i.e. multiplied by inverses) */\
	/* The add2-addressed cosine ratios are arranged in 3 YMM-register/memory-sized slots like so:
	  once we have filled 4 YYMs with inverses 1/[c3,c1-15] and used those to get the 16 tangents (1st set = 1/c3
	  and discarded) we will do as described in the right column to set up for the cosine-ratios computation:

		double __c31 = __c3/__c1;
		double __c51 = __c5/__c1;
		double __c62 = __c6/__c2;
		[0 pad]						shuffle YMM with 1/[c3,c1,c2,c3] to get 1/[c1,c1,c2,c3], then *= [c3,c5,c6,0]

		double __c73 = __c7/__c3;
		double __c91 = __c9/__c1;
		double __cA2 = __cA/__c2;
		double __cB3 = __cB/__c3;	initialize YMM with 1/[c3,c1,c2,c3], then *= [c7,c9,cA,cB]

		double __cC4 = __cC/__c4;
		double __cD5 = __cD/__c5;
		double __cE6 = __cE/__c6;
		double __cF7 = __cF/__c7;	Multiply YMM with 1/[c4-7] *= [cC-F]
	*/\
		*add0++ = 0.0;	/* Since tan0 defined as const, use this pair of double slots to hold 1/c3 (via c3,1 on input, then invert c3 and multiply */\
		*add1++ = 1.0;	/* them together), which extra 1/c3 copy saves some really awkward permuting, at least in terms of the idiotic x86 ISA. */\
\
		*add0++ = __c1;	/* c1, for inversion  */\
		*add1++ = __s1;	/* s1  slot will hold __r1 = s1 /c1  */\
\
		*add0++ = __c2;	/* c2, for inversion  */\
		*add1++ = __s2;	/* s2  slot will hold __r2 = s2 /c2  */\
\
		*(add0-3) = __c3;	/* c3, for inversion ... */\
		*add0++   = __c3;	/* place extra copy in 0-slot as described above - put on separate line to avoid ambiguity of *(add0-3) = *add0++ = __c3. */\
		*add1++ = __s3;	/* s3  slot will hold __r3 = s3 /c3  */\
		*add2++ = __c3;	/* c3, will get multiplied by 1/c1 to yield __c31 */\
\
		*add0++ = __c4;	/* c4, for inversion  */\
		*add1++ = __s4;	/* s4  slot will hold __r4 = s4 /c4  */\
\
		*add0++ = __c5;	/* c5, for inversion  */\
		*add1++ = __s5;	/* s5  slot will hold __r5 = s5 /c5  */\
		*add2++ = __c5;	/* c5, will get multiplied by 1/c1 to yield __c51 */\
\
		*add0++ = __c6;	/* c6, for inversion  */\
		*add1++ = __s6;	/* s6  slot will hold __r6 = s6 /c6  */\
		*add2++ = __c6;	/* c6, will get multiplied by 1/c2 to yield __c62 */\
		*add2++ = 0.0;	/* 0-pad will get multiplied by 1/c3 term, remains 0-pad. */\
\
		*add0++ = __c7;	/* c7, for inversion  */\
		*add1++ = __s7;	/* s7  slot will hold __r7 = s7 /c7  */\
		*add2++ = __c7;	/* c7, will get multiplied by 1/c3 to yield __c73 */\
\
		*add0++ = __c8;	/* c8, for inversion  */\
		*add1++ = __s8;	/* s8  slot will hold __r8 = s8 /c8  */\
\
		*add0++ = __c9;	/* c9, for inversion  */\
		*add1++ = __s9;	/* s9  slot will hold __r9 = s9 /c9  */\
		*add2++ = __c9;	/* c9, will get multiplied by 1/c1 to yield __c91 */\
\
		*add0++ = __cA;	/* c10, for inversion  */\
		*add1++ = __sA;	/* s10 slot will hold __rA = s10/c10 */\
		*add2++ = __cA;	/* cA, will get multiplied by 1/c2 to yield __cA2 */\
\
		*add0++ = __cB;	/* c11, for inversion  */\
		*add1++ = __sB;	/* s11 slot will hold __rB = s11/c11 */\
		*add2++ = __cB;	/* cB, will get multiplied by 1/c3 to yield __cB3 */\
\
		*add0++ = __cC;	/* c12, for inversion  */\
		*add1++ = __sC;	/* s12 slot will hold __rC = s12/c12 */\
		*add2++ = __cC;	/* cC, will get multiplied by 1/c4 to yield __cC4 */\
\
		*add0++ = __cD;	/* c13, for inversion  */\
		*add1++ = __sD;	/* s13 slot will hold __rD = s13/c13 */\
		*add2++ = __cD;	/* cD, will get multiplied by 1/c5 to yield __cD5 */\
\
		*add0++ = __cE;	/* c14, for inversion  */\
		*add1++ = __sE;	/* s14 slot will hold __rE = s14/c14 */\
		*add2++ = __cE;	/* cE, will get multiplied by 1/c6 to yield __cE6 */\
\
		*add0++ = __cF;	/* c15, for inversion  */\
		*add1++ = __sF;	/* s15 slot will hold __rF = s15/c15 */\
		*add2++ = __cF;	/* cF, will get multiplied by 1/c7 to yield __cF7 */\
\
		/* This places us at add0 == c8 and add1 = c12. */\
		ASSERT(add0 == (double *)__twid_ptr+16 && add1 == (double *)__twid_ptr+32 && add2 == (double *)__twid_ptr+44, "add0,1,2 checksum failed in AVX2 sincos inits!");\
	/*
	At this point, the 11 ymm-sized [32-byte] chunks starting at &__twid_ptr contain the following scalar-double data:

		0:	c3,c1-3
		1:	c4-7
		2:	c8-11
		3:	c12-c15
		4:	1.0,s1-3
		5:	s4-7
		6:	s8-11
		7:	s12-s15
		8:	c3,5,6,[0-pad]
		9:	c7,9-B
		A:	cC-F
	*/\
\
	/* Now send the cosine terms to the inversion routine, which also does the combine-and-populate-SIMD-slots step. */\
		RADIX16_COMPUTE_FMA_SINCOS_DIF_2(__twid_ptr,one);\
		add0 = (double *)__twid_ptr;\
\
	/* Scalar data starting at add0 = __twid_ptr laid out as below:

		a[0x00,0x01,0x02,0x03]: add0 + 0x[  0,  8, 10, 18]: [c3 ,c1 ,c2 ,c3 ] Cosines:
		a[0x04,0x05,0x06,0x07]: add0 + 0x[ 20, 28, 30, 38]: [c4 ,c5 ,c6 ,c7 ]
		a[0x08,0x09,0x0a,0x0b]: add0 + 0x[ 40, 48, 50, 58]: [c8 ,c9 ,cA ,cB ]
		a[0x0c,0x0d,0x0e,0x0f]: add0 + 0x[ 60, 68, 70, 78]: [cC ,cD ,cE ,cF ]
		a[0x10,0x11,0x12,0x13]: add0 + 0x[ 80, 88, 90, 98]: [-- ,r1 ,r2 ,r3 ] Tangents:
		a[0x14,0x15,0x16,0x17]: add0 + 0x[ a0, a8, b0, b8]: [r4 ,r5 ,r6 ,r7 ]
		a[0x18,0x19,0x1a,0x1b]: add0 + 0x[ c0, c8, d0, d8]: [r8 ,r9 ,rA ,rB ]
		a[0x1c,0x1d,0x1e,0x1f]: add0 + 0x[ e0, e8, f0, f8]: [rC ,rD ,rE ,rF ]
		a[0x20,0x21,0x22,0x23]: add0 + 0x[100,108,110,118]: [c31,c51,c62,  0] Cosine ratios:
		a[0x24,0x25,0x26,0x27]: add0 + 0x[120,128,130,138]: [c73,c91,cA2,cB3]
		a[0x28,0x29,0x2a,0x2b]: add0 + 0x[140,148,150,158]: [cC4,cD5,cE6,cF7]

	Ensuing C code massages the above into a scalar-data analog of the 4-copy layout.
	*/\
	/* put all overwrites-of-no-longer-needed data first, to minimize conflicts later.
	Data which will not be used in the FMA-based radix-16 DIF are at indices 0x[0,3,5-7,9-15,16]:
	Arrange the rest so RHS (read-elt) indices are ascending, then manually move
	up those which overwrite indices appearing further down in the RHS:
	*/\
		add0[0x00] = add0[0x01];	/* c1, copy to *= __c */\
		add0[0x03] = add0[0x02];	/* c2, copy to *= ISRT2 */\
		add0[0x0a] = add0[0x02];	/* copy c2 to final loc before overwriting */\
		add0[0x02] = add0[0x01];	/* c1, copy to *= ISRT2 */\
		add0[0x01] = __tan;\
		add0[0x06] = add0[0x04];	/* c4 */\
		add0[0x04] = add0[0x08];	/* c8 */\
		add0[0x05] = add0[0x18];	/* r8 */\
		add0[0x07] = add0[0x14];	/* r4 */\
		add0[0x08] = add0[0x28];	/* cC4 */\
		add0[0x09] = add0[0x1c];	/* rC */\
		add0[0x0b] = add0[0x12];	/* r2 */\
		add0[0x0c] = add0[0x26];	/* cA2 */\
		add0[0x0d] = add0[0x1a];	/* rA */\
		add0[0x0e] = add0[0x22];	/* c62 */\
		add0[0x0f] = add0[0x16];	/* r6 */\
		add0[0x10] = add0[0x2a];	/* cE6 */\
		add0[0x16] = add0[0x21];	/* c51 */\
		add0[0x21] = add0[0x1f];	/* rF */\
		add0[0x1f] = add0[0x17];	/* r7 */\
		add0[0x17] = add0[0x15];	/* r5 */\
		add0[0x15] = add0[0x19];	/* r9 */\
		add0[0x19] = add0[0x1d];	/* rD */\
		add0[0x1d] = add0[0x1b];	/* rB */\
		add0[0x1b] = add0[0x13];	/* r3 */\
		add0[0x13] = add0[0x11];	/* r1 */\
		add0[0x11] = add0[0x1e];	/* rE */\
		add0[0x12] = add0[0x00];	/* c1 now in [0] */\
		add0[0x14] = add0[0x25];	/* c91 */\
		add0[0x18] = add0[0x29];	/* cD5 */\
		add0[0x1a] = add0[0x20];	/* c31 */\
		add0[0x1c] = add0[0x27];	/* cB3 */\
		add0[0x1e] = add0[0x24];	/* c73 */\
		add0[0x20] = add0[0x2b];	/* cF7 */\
		/* Now mpy elts in slots 0,2,3 by __c, ISRT2, ISRT2, respectively: */\
		add0[0x00] *= __c;\
		add0[0x02] *= ISRT2;\
		add0[0x03] *= ISRT2;\
		/* And stick a 1.0 at the end of the above block-of-doubles: */\
		add0[0x22] = 1.0;\
\
	/* Yielding the following layout-of-scalar-doubles, with data above index 0x22 unused in the  DIF DFT:

		a[0x00,0x01,0x02,0x03]: add0 + 0x[  0,  8, 10, 18]: [c1*c,s/c,c1*ISRT2,c2*ISRT2]
		a[0x04,0x05,0x06,0x07]: add0 + 0x[ 20, 28, 30, 38]: [c8 ,r8 ,c4 ,r4 ]
		a[0x08,0x09,0x0a,0x0b]: add0 + 0x[ 40, 48, 50, 58]: [cC4,rC ,c2 ,r2 ]
		a[0x0c,0x0d,0x0e,0x0f]: add0 + 0x[ 60, 68, 70, 78]: [cA2,rA ,c62,r6 ]
		a[0x10,0x11,0x12,0x13]: add0 + 0x[ 80, 88, 90, 98]: [cE6,rE ,c1 ,r1 ]
		a[0x14,0x15,0x16,0x17]: add0 + 0x[ a0, a8, b0, b8]: [c91,r9 ,c51,r5 ]
		a[0x18,0x19,0x1a,0x1b]: add0 + 0x[ c0, c8, d0, d8]: [cD5,rD ,c31,r3 ]
		a[0x1c,0x1d,0x1e,0x1f]: add0 + 0x[ e0, e8, f0, f8]: [cB3,rB ,c73,r7 ]
		a[0x20,0x21,0x22,0x23]: add0 + 0x[100,108,110,118]: [cF7,rF ,1.0,  0]
		a[0x24,0x25,0x26,0x27]: add0 + 0x[120,128,130,138]: [c73,c91,cA2,cB3]
		a[0x28,0x29,0x2a,0x2b]: add0 + 0x[140,148,150,158]: [cC4,cD5,cE6,cF7]
	*/\
	}

#endif	// avx2?

#ifdef USE_ARM_V8_SIMD

	// Macro name here is from the analogous SSE2 macro; for ARMv8 'V2' is in fact the only 16-DFT macro.
	// Vector-opcount: 50 LDP, 32 STP, 128 ADD, 96 MUL
	#define SSE2_RADIX16_DIF_TWIDDLE_V2(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr0,Xtwo,Xcc0,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	"ldr	x15,%[__add0]		\n\t"/* base address for 4 prefetches-from-main-data-array spread through this macro */\
	"ldr	w14,%[__pfetch_dist]	\n\t"/* data-fetch-ahead byte-offset */\
	"add	x15,x15,x14			\n\t"\
	"prfm	PLDL1KEEP,[x15]		\n\t"/* [base-address + data-fetch-ahead byte-offset] + p0 */\
	/*
		k1 = p2*8;	// Set k1 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k2 = 0x80;	// Set k2 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		// Pass1-DFTs all use same twiddle-triplet, so lit-byte address offset between ptrs to 1st,2nd DFT's twiddles = 0
		i0 = (vec_dbl *)add0; i1 = (vec_dbl *)(add0+p4); i2 = (vec_dbl *)(add0+p8); i3 = (vec_dbl *)(add0+p12);
		o0 = r1; o1 = r1+2; o2 = r1+4; o3 = r1+6;
		c_tmp = cc0;	// c8,4,C
		SSE2_RADIX_04_DIF_3TWIDDLE_X1(i0,i1,i2,i3,k1, two,c_tmp, o0,o1,o2,o3,k2)
	*/\
		"ldr	x4,%[__cc0]			\n\t	ldr	x14,%[__r0]				\n\t"\
		"ldr	x0,%[__add0]		\n\t	ldr	w6,%[__p1]				\n\t"\
		"ldr	w1,%[__p4]			\n\t	add	x1 , x0,x1,lsl #3		\n\t"/* add0 + p4 */\
		"ldr	w2,%[__p8]			\n\t	add	x2 , x0,x2,lsl #3		\n\t"/* add0 + p8 */\
		"ldr	w3,%[__p12]			\n\t	add	x3 , x0,x3,lsl #3		\n\t"/* add0 + p12 */\
		"ldr	w5,%[__p2]			\n\t	add	x10, x0,x5,lsl #3		\n\t"/* add0 + p2 */\
		"ldp	q4,q5,[x2]			\n\t	add	x11, x1,x5,lsl #3		\n\t"/* add0 + p6 */\
/* c0 */"ldp	q8,q9,[x4] 			\n\t	add	x12, x2,x5,lsl #3		\n\t"/* add0 + p10 */\
		"ldp	q0,q1,[x0]			\n\t	add	x13, x3,x5,lsl #3		\n\t"/* add0 + p14 */\
		"fmul	v6.2d,v4.2d,v8.2d	\n\t	ldp	q16,q17,[x12]			\n\t"\
		"fmul	v7.2d,v5.2d,v8.2d	\n\t	ldp	q12,q13,[x10]			\n\t"\
		"fmls	v6.2d,v5.2d,v9.2d	\n\t"\
		"fmla	v7.2d,v4.2d,v9.2d	\n\t"\
		"fsub	v2.2d ,v0.2d,v6.2d	\n\t fmul	v18.2d,v16.2d,v8.2d		\n\t"/* twiddle-mul: */\
		"fsub	v3.2d ,v1.2d,v7.2d	\n\t fmul	v19.2d,v17.2d,v8.2d		\n\t"\
		"fadd	v10.2d,v0.2d,v6.2d	\n\t fmls	v18.2d,v17.2d,v9.2d		\n\t"\
		"fadd	v11.2d,v1.2d,v7.2d	\n\t fmla	v19.2d,v16.2d,v9.2d		\n\t"\
/*c0+4*/"ldp	q8,q9,[x4,#0x40]	\n\t fsub	v14.2d,v12.2d,v18.2d	\n\t"/* 2 x 2 complex butterfly: */\
		"ldp	q6,q7,[x3]			\n\t fsub	v15.2d,v13.2d,v19.2d	\n\t"\
		"fmul	v0.2d,v6.2d,v8.2d	\n\t fadd	v22.2d,v12.2d,v18.2d	\n\t"\
		"fmul	v1.2d,v7.2d,v8.2d	\n\t fadd	v23.2d,v13.2d,v19.2d	\n\t"\
		"fmls	v0.2d,v7.2d,v9.2d	\n\t	mov	v20.16b,v8.16b	\n\t"/* Need copies of v8,9 twiddle-data because */\
		"fmla	v1.2d,v6.2d,v9.2d	\n\t	mov	v21.16b,v9.16b	\n\t"/* overwrite those in lcol immediately below */\
/*c0+2*/"ldp	q8,q9,[x4,#0x20]	\n\t	ldp	q18,q19,[x13]			\n\t"\
		"ldp	q6,q7,[x1]			\n\t fmul	v12.2d,v18.2d,v20.2d	\n\t"/* twiddle-mul: */\
		"fmul	v4.2d,v6.2d,v8.2d	\n\t fmul	v13.2d,v19.2d,v20.2d	\n\t"\
		"fmul	v5.2d,v7.2d,v8.2d	\n\t fmls	v12.2d,v19.2d,v21.2d	\n\t"\
		"fmls	v4.2d,v7.2d,v9.2d	\n\t fmla	v13.2d,v18.2d,v21.2d	\n\t"\
		"fmla	v5.2d,v6.2d,v9.2d	\n\t	ldp	q18,q19,[x11]			\n\t"\
		"fadd	v6.2d,v4.2d,v0.2d	\n\t fmul	v16.2d,v18.2d,v8.2d		\n\t"/* twiddle-mul: */\
		"fadd	v7.2d,v5.2d,v1.2d	\n\t fmul	v17.2d,v19.2d,v8.2d		\n\t"\
		"fsub	v4.2d,v4.2d,v0.2d	\n\t fmls	v16.2d,v19.2d,v9.2d		\n\t"\
		"fsub	v5.2d,v5.2d,v1.2d	\n\t fmla	v17.2d,v18.2d,v9.2d		\n\t"\
		"fsub	v8.2d,v10.2d,v6.2d	\n\t fadd	v18.2d,v16.2d,v12.2d	\n\t"/* 2 x 2 complex butterfly: */\
		"fsub	v9.2d,v11.2d,v7.2d	\n\t fadd	v19.2d,v17.2d,v13.2d	\n\t"\
		"fsub	v1.2d,v3.2d,v4.2d	\n\t fsub	v16.2d,v16.2d,v12.2d	\n\t"\
		"fsub	v0.2d,v2.2d,v5.2d	\n\t fsub	v17.2d,v17.2d,v13.2d	\n\t"\
		"fadd	v6.2d,v6.2d,v10.2d	\n\t fsub	v20.2d,v22.2d,v18.2d	\n\t"\
		"fadd	v7.2d,v7.2d,v11.2d	\n\t fsub	v21.2d,v23.2d,v19.2d	\n\t"\
		"fadd	v4.2d,v4.2d,v3.2d	\n\t fsub	v13.2d,v15.2d,v16.2d	\n\t"\
		"fadd	v5.2d,v5.2d,v2.2d	\n\t fsub	v12.2d,v14.2d,v17.2d	\n\t"\
		"stp	q6,q7,[x14      ]	\n\t fadd	v18.2d,v18.2d,v22.2d	\n\t"\
		"stp	q0,q4,[x14,#0x20]	\n\t fadd	v19.2d,v19.2d,v23.2d	\n\t"\
		"stp	q8,q9,[x14,#0x40]	\n\t fadd	v16.2d,v16.2d,v15.2d	\n\t"\
		"stp	q5,q1,[x14,#0x60]	\n\t fadd	v17.2d,v17.2d,v14.2d	\n\t"\
		"add	x0, x0,x6,lsl #3	\n\t	stp	q18,q19,[x14,#0x80]		\n\t"/* out 0 */\
		"add	x1, x1,x6,lsl #3	\n\t	stp	q12,q16,[x14,#0xa0]		\n\t"/* out 1 */\
		"add	x2, x2,x6,lsl #3	\n\t	stp	q20,q21,[x14,#0xc0]		\n\t"/* out 2 */\
		"add	x3, x3,x6,lsl #3	\n\t	stp	q17,q13,[x14,#0xe0]		\n\t"/* out 3 */\
	/*
		i0 = (vec_dbl *)(add0+p1); i1 = (vec_dbl *)(add0+p4+p1); i2 = (vec_dbl *)(add0+p8+p1); i3 = (vec_dbl *)(add0+p12+p1);
		o0 += 16; o1 += 16; o2 += 16; o3 += 16;
		SSE2_RADIX_04_DIF_3TWIDDLE_X1(i0,i1,i2,i3,k1, two,c_tmp, o0,o1,o2,o3,k2)
	*/\
	"prfm	PLDL1KEEP,[x15,x2,LSL #3]	\n\t"/* ... + p8 */\
		"ldp	q4,q5,[x2]			\n\t	add	x10,x10,x6,lsl #3		\n\t"/* add0 + p1 ,3 */\
/* c0 */"ldp	q8,q9,[x4] 			\n\t	add	x11,x11,x6,lsl #3		\n\t"/* add0 + p5 ,7 */\
		"ldp	q0,q1,[x0]			\n\t	add	x12,x12,x6,lsl #3		\n\t"/* add0 + p9 ,11 */\
		"fmul	v6.2d,v4.2d,v8.2d	\n\t	add	x13,x13,x6,lsl #3		\n\t"/* add0 + p13,15 */\
		"fmul	v7.2d,v5.2d,v8.2d	\n\t	add	x14,x14,#0x100			\n\t"\
		"fmls	v6.2d,v5.2d,v9.2d	\n\t	ldp	q16,q17,[x12]			\n\t"\
		"fmla	v7.2d,v4.2d,v9.2d	\n\t	ldp	q12,q13,[x10]			\n\t"\
		"fsub	v2.2d ,v0.2d,v6.2d	\n\t fmul	v18.2d,v16.2d,v8.2d		\n\t"/* twiddle-mul: */\
		"fsub	v3.2d ,v1.2d,v7.2d	\n\t fmul	v19.2d,v17.2d,v8.2d		\n\t"\
		"fadd	v10.2d,v0.2d,v6.2d	\n\t fmls	v18.2d,v17.2d,v9.2d		\n\t"\
		"fadd	v11.2d,v1.2d,v7.2d	\n\t fmla	v19.2d,v16.2d,v9.2d		\n\t"\
/*c0+4*/"ldp	q8,q9,[x4,#0x40]	\n\t fsub	v14.2d,v12.2d,v18.2d	\n\t"/* 2 x 2 complex butterfly: */\
		"ldp	q6,q7,[x3]			\n\t fsub	v15.2d,v13.2d,v19.2d	\n\t"\
		"fmul	v0.2d,v6.2d,v8.2d	\n\t fadd	v22.2d,v12.2d,v18.2d	\n\t"\
		"fmul	v1.2d,v7.2d,v8.2d	\n\t fadd	v23.2d,v13.2d,v19.2d	\n\t"\
		"fmls	v0.2d,v7.2d,v9.2d	\n\t	mov	v20.16b,v8.16b	\n\t"/* Need copies of v8,9 twiddle-data because */\
		"fmla	v1.2d,v6.2d,v9.2d	\n\t	mov	v21.16b,v9.16b	\n\t"/* overwrite those in lcol immediately below */\
/*c0+2*/"ldp	q8,q9,[x4,#0x20]	\n\t	ldp	q18,q19,[x13]			\n\t"\
		"ldp	q6,q7,[x1]			\n\t fmul	v12.2d,v18.2d,v20.2d	\n\t"/* twiddle-mul: */\
		"fmul	v4.2d,v6.2d,v8.2d	\n\t fmul	v13.2d,v19.2d,v20.2d	\n\t"\
		"fmul	v5.2d,v7.2d,v8.2d	\n\t fmls	v12.2d,v19.2d,v21.2d	\n\t"\
		"fmls	v4.2d,v7.2d,v9.2d	\n\t fmla	v13.2d,v18.2d,v21.2d	\n\t"\
		"fmla	v5.2d,v6.2d,v9.2d	\n\t	ldp	q18,q19,[x11]			\n\t"\
		"fadd	v6.2d,v4.2d,v0.2d	\n\t fmul	v16.2d,v18.2d,v8.2d		\n\t"/* twiddle-mul: */\
		"fadd	v7.2d,v5.2d,v1.2d	\n\t fmul	v17.2d,v19.2d,v8.2d		\n\t"\
		"fsub	v4.2d,v4.2d,v0.2d	\n\t fmls	v16.2d,v19.2d,v9.2d		\n\t"\
		"fsub	v5.2d,v5.2d,v1.2d	\n\t fmla	v17.2d,v18.2d,v9.2d		\n\t"\
		"fsub	v8.2d,v10.2d,v6.2d	\n\t fadd	v18.2d,v16.2d,v12.2d	\n\t"/* 2 x 2 complex butterfly: */\
		"fsub	v9.2d,v11.2d,v7.2d	\n\t fadd	v19.2d,v17.2d,v13.2d	\n\t"\
		"fsub	v1.2d,v3.2d,v4.2d	\n\t fsub	v16.2d,v16.2d,v12.2d	\n\t"\
		"fsub	v0.2d,v2.2d,v5.2d	\n\t fsub	v17.2d,v17.2d,v13.2d	\n\t"\
		"fadd	v6.2d,v6.2d,v10.2d	\n\t fsub	v20.2d,v22.2d,v18.2d	\n\t"\
		"fadd	v7.2d,v7.2d,v11.2d	\n\t fsub	v21.2d,v23.2d,v19.2d	\n\t"\
		"fadd	v4.2d,v4.2d,v3.2d	\n\t fsub	v13.2d,v15.2d,v16.2d	\n\t"\
		"fadd	v5.2d,v5.2d,v2.2d	\n\t fsub	v12.2d,v14.2d,v17.2d	\n\t"\
		"stp	q6,q7,[x14      ]	\n\t fadd	v18.2d,v18.2d,v22.2d	\n\t"\
		"stp	q0,q4,[x14,#0x20]	\n\t fadd	v19.2d,v19.2d,v23.2d	\n\t"\
		"stp	q8,q9,[x14,#0x40]	\n\t fadd	v16.2d,v16.2d,v15.2d	\n\t"\
		"stp	q5,q1,[x14,#0x60]	\n\t fadd	v17.2d,v17.2d,v14.2d	\n\t"\
		"mov x2, x0				\n\t"/* o2 = add0 + p1 */"stp q18,q19,[x14,#0x80]	\n\t"/* out 0 */\
		"mov x3,x10				\n\t"/* o3 = add0 + p3 */"stp q12,q16,[x14,#0xa0]	\n\t"/* out 1 */\
		"sub x0, x0,x6,lsl #3	\n\t"/* o0 = add0 + p0 */"stp q20,q21,[x14,#0xc0]	\n\t"/* out 2 */\
		"add x1, x0,x5,lsl #3	\n\t"/* o1 = add0 + p2 */"stp q17,q13,[x14,#0xe0]	\n\t"/* out 3 */\
	/*
		// Pass 2:
		k1 = 0x40;	// Set k1 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k2 = p4*8;	// Set k2 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		// 0x60 = literal-bytewise address offset between the pointers to the first and second DFT's twiddles = 0:
		i0 = r1; i1 = r1+16; i2 = r1+8; i3 = r1+24;
		o0 = (vec_dbl *)add0; o1 = (vec_dbl *)(add0+p2); o2 = (vec_dbl *)(add0+p1); o3 = (vec_dbl *)(add0+p3);
		c_tmp += 6;	// c2,A,6
		SSE2_RADIX_04_DIF_3TWIDDLE_X2(i0,i1,i2,i3,k1, two,c_tmp,0x60, o0,o1,o2,o3,k2)
	*/\
	"prfm	PLDL1KEEP,[x15,x1,LSL #3]	\n\t"/* ... + p4 */\
		"ldr	x14,%[__r0]			\n\t	add	x4, x4,#0x60			\n\t"/* cc0 + 0x60 */\
		"ldr	w5,%[__p4]			\n\t	add	x10, x0,x5,lsl #3		\n\t"/* add0 + p4 */\
		"ldp	q4,q5,[x14,#0x080]	\n\t	add	x11, x1,x5,lsl #3		\n\t"/* add0 + p6 */\
		"ldp	q8,q9,[x4]			\n\t	add	x12, x2,x5,lsl #3		\n\t"/* add0 + p5 */\
		"ldp	q0,q1,[x14       ]	\n\t	add	x13, x3,x5,lsl #3		\n\t"/* add0 + p7 */\
		"fmul	v6.2d,v4.2d,v8.2d	\n\t	ldr	w5,%[__p8]		\n\t"/* For next 4-DFT */\
		"fmul	v7.2d,v5.2d,v8.2d	\n\t	ldp	q16,q17,[x14,#0x0c0]	\n\t"\
		"fmls	v6.2d,v5.2d,v9.2d	\n\t	ldp	q20,q21,[x4,#0x60]		\n\t"/* cc0 */\
		"fmla	v7.2d,v4.2d,v9.2d	\n\t	ldp	q12,q13,[x14,#0x040]	\n\t"\
		"fsub	v2.2d ,v0.2d,v6.2d	\n\t fmul	v18.2d,v16.2d,v20.2d	\n\t"/* twiddle-mul: */\
		"fsub	v3.2d ,v1.2d,v7.2d	\n\t fmul	v19.2d,v17.2d,v20.2d	\n\t"\
		"fadd	v10.2d,v0.2d,v6.2d	\n\t fmls	v18.2d,v17.2d,v21.2d	\n\t"\
		"fadd	v11.2d,v1.2d,v7.2d	\n\t fmla	v19.2d,v16.2d,v21.2d	\n\t"\
		"ldp	q8,q9,[x4,#0x40]	\n\t fsub	v14.2d,v12.2d,v18.2d	\n\t"/* 2 x 2 complex butterfly: */\
		"ldp	q6,q7,[x14,#0x180]	\n\t fsub	v15.2d,v13.2d,v19.2d	\n\t"\
		"fmul	v0.2d,v6.2d,v8.2d	\n\t fadd	v22.2d,v12.2d,v18.2d	\n\t"\
		"fmul	v1.2d,v7.2d,v8.2d	\n\t fadd	v23.2d,v13.2d,v19.2d	\n\t"\
		"fmls	v0.2d,v7.2d,v9.2d	\n\t	ldp	q20,q21,[x4,#0xa0]		\n\t"/* cc0+4 */\
		"fmla	v1.2d,v6.2d,v9.2d	\n\t	ldp	q18,q19,[x14,#0x1c0]	\n\t"\
		"ldp	q8,q9,[x4,#0x20]	\n\t fmul	v12.2d,v18.2d,v20.2d	\n\t"/* twiddle-mul: */\
		"ldp	q6,q7,[x14,#0x100]	\n\t fmul	v13.2d,v19.2d,v20.2d	\n\t"\
		"fmul	v4.2d,v6.2d,v8.2d	\n\t fmls	v12.2d,v19.2d,v21.2d	\n\t"\
		"fmul	v5.2d,v7.2d,v8.2d	\n\t fmla	v13.2d,v18.2d,v21.2d	\n\t"\
		"fmls	v4.2d,v7.2d,v9.2d	\n\t	ldp	q20,q21,[x4,#0x80]		\n\t"/* cc0+2 */\
		"fmla	v5.2d,v6.2d,v9.2d	\n\t	ldp	q18,q19,[x14,#0x140]	\n\t"\
		"fadd	v6.2d,v4.2d,v0.2d	\n\t fmul	v16.2d,v18.2d,v20.2d	\n\t"/* twiddle-mul: */\
		"fadd	v7.2d,v5.2d,v1.2d	\n\t fmul	v17.2d,v19.2d,v20.2d	\n\t"\
		"fsub	v4.2d,v4.2d,v0.2d	\n\t fmls	v16.2d,v19.2d,v21.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v1.2d	\n\t fmla	v17.2d,v18.2d,v21.2d	\n\t"\
		"fsub	v8.2d,v10.2d,v6.2d	\n\t fadd	v18.2d,v16.2d,v12.2d	\n\t"/* 2 x 2 complex butterfly: */\
		"fsub	v9.2d,v11.2d,v7.2d	\n\t fadd	v19.2d,v17.2d,v13.2d	\n\t"\
		"fsub	v1.2d,v3.2d,v4.2d	\n\t fsub	v16.2d,v16.2d,v12.2d	\n\t"\
		"fsub	v0.2d,v2.2d,v5.2d	\n\t fsub	v17.2d,v17.2d,v13.2d	\n\t"\
		"fadd	v6.2d,v6.2d,v10.2d	\n\t fsub	v20.2d,v22.2d,v18.2d	\n\t"\
		"fadd	v7.2d,v7.2d,v11.2d	\n\t fsub	v21.2d,v23.2d,v19.2d	\n\t"\
		"fadd	v4.2d,v4.2d,v3.2d	\n\t fsub	v13.2d,v15.2d,v16.2d	\n\t"\
		"fadd	v5.2d,v5.2d,v2.2d	\n\t fsub	v12.2d,v14.2d,v17.2d	\n\t"\
		"stp	q6,q7,[x0]			\n\t fadd	v18.2d,v18.2d,v22.2d	\n\t"\
		"stp	q0,q4,[x1]			\n\t fadd	v19.2d,v19.2d,v23.2d	\n\t"\
		"stp	q8,q9,[x2]			\n\t fadd	v16.2d,v16.2d,v15.2d	\n\t"\
		"stp	q5,q1,[x3]			\n\t fadd	v17.2d,v17.2d,v14.2d	\n\t"\
		"add x0, x0,x5,lsl #3	\n\t"/* add0 + p8  */"stp q18,q19,[x10]	\n\t"/* out 0 */\
		"add x1, x1,x5,lsl #3	\n\t"/* add0 + p10 */"stp q12,q16,[x11]	\n\t"/* out 1 */\
		"add x2, x2,x5,lsl #3	\n\t"/* add0 + p9  */"stp q20,q21,[x12]	\n\t"/* out 2 */\
		"add x3, x3,x5,lsl #3	\n\t"/* add0 + p11 */"stp q17,q13,[x13]	\n\t"/* out 3 */\
	/*
		i0 += 2; i1 += 2; i2 += 2; i3 += 2;
		o0 = (vec_dbl *)(add0+p8); o1 = (vec_dbl *)(add0+p8+p2); o2 = (vec_dbl *)(add0+p8+p1); o3 = (vec_dbl *)(add0+p8+p3);
		c_tmp += 12;	// c5,D,3
		SSE2_RADIX_04_DIF_3TWIDDLE_X2(i0,i1,i2,i3,k1, two,c_tmp,0x60, o0,o1,o2,o3,k2)
	*/\
	"prfm	PLDL1KEEP,[x15,x3,LSL #3]	\n\t"/* ... + pC */\
		"add	x4, x4,#0xc0		\n\t	add	x10,x10,x5,lsl #3		\n\t"/* add0 + p12 */\
		"ldp	q4,q5,[x14,#0x0a0]	\n\t	add	x11,x11,x5,lsl #3		\n\t"/* add0 + p14 */\
		"ldp	q8,q9,[x4]			\n\t	add	x12,x12,x5,lsl #3		\n\t"/* add0 + p13 */\
		"ldp	q0,q1,[x14,#0x020]	\n\t	add	x13,x13,x5,lsl #3		\n\t"/* add0 + p15 */\
		"fmul	v6.2d,v4.2d,v8.2d	\n\t"\
		"fmul	v7.2d,v5.2d,v8.2d	\n\t	ldp	q16,q17,[x14,#0x0e0]	\n\t"\
		"fmls	v6.2d,v5.2d,v9.2d	\n\t	ldp	q20,q21,[x4,#0x60]		\n\t"/* cc0 */\
		"fmla	v7.2d,v4.2d,v9.2d	\n\t	ldp	q12,q13,[x14,#0x060]	\n\t"\
		"fsub	v2.2d ,v0.2d,v6.2d	\n\t fmul	v18.2d,v16.2d,v20.2d	\n\t"/* twiddle-mul: */\
		"fsub	v3.2d ,v1.2d,v7.2d	\n\t fmul	v19.2d,v17.2d,v20.2d	\n\t"\
		"fadd	v10.2d,v0.2d,v6.2d	\n\t fmls	v18.2d,v17.2d,v21.2d	\n\t"\
		"fadd	v11.2d,v1.2d,v7.2d	\n\t fmla	v19.2d,v16.2d,v21.2d	\n\t"\
		"ldp	q8,q9,[x4,#0x40]	\n\t fsub	v14.2d,v12.2d,v18.2d	\n\t"/* 2 x 2 complex butterfly: */\
		"ldp	q6,q7,[x14,#0x1a0]	\n\t fsub	v15.2d,v13.2d,v19.2d	\n\t"\
		"fmul	v0.2d,v6.2d,v8.2d	\n\t fadd	v22.2d,v12.2d,v18.2d	\n\t"\
		"fmul	v1.2d,v7.2d,v8.2d	\n\t fadd	v23.2d,v13.2d,v19.2d	\n\t"\
		"fmls	v0.2d,v7.2d,v9.2d	\n\t	ldp	q20,q21,[x4,#0xa0]		\n\t"/* cc0+4 */\
		"fmla	v1.2d,v6.2d,v9.2d	\n\t	ldp	q18,q19,[x14,#0x1e0]	\n\t"\
		"ldp	q8,q9,[x4,#0x20]	\n\t fmul	v12.2d,v18.2d,v20.2d	\n\t"/* twiddle-mul: */\
		"ldp	q6,q7,[x14,#0x120]	\n\t fmul	v13.2d,v19.2d,v20.2d	\n\t"\
		"fmul	v4.2d,v6.2d,v8.2d	\n\t fmls	v12.2d,v19.2d,v21.2d	\n\t"\
		"fmul	v5.2d,v7.2d,v8.2d	\n\t fmla	v13.2d,v18.2d,v21.2d	\n\t"\
		"fmls	v4.2d,v7.2d,v9.2d	\n\t	ldp	q20,q21,[x4,#0x80]		\n\t"/* cc0+2 */\
		"fmla	v5.2d,v6.2d,v9.2d	\n\t	ldp	q18,q19,[x14,#0x160]	\n\t"\
		"fadd	v6.2d,v4.2d,v0.2d	\n\t fmul	v16.2d,v18.2d,v20.2d	\n\t"/* twiddle-mul: */\
		"fadd	v7.2d,v5.2d,v1.2d	\n\t fmul	v17.2d,v19.2d,v20.2d	\n\t"\
		"fsub	v4.2d,v4.2d,v0.2d	\n\t fmls	v16.2d,v19.2d,v21.2d	\n\t"\
		"fsub	v5.2d,v5.2d,v1.2d	\n\t fmla	v17.2d,v18.2d,v21.2d	\n\t"\
		"fsub	v8.2d,v10.2d,v6.2d	\n\t fadd	v18.2d,v16.2d,v12.2d	\n\t"/* 2 x 2 complex butterfly: */\
		"fsub	v9.2d,v11.2d,v7.2d	\n\t fadd	v19.2d,v17.2d,v13.2d	\n\t"\
		"fsub	v1.2d,v3.2d,v4.2d	\n\t fsub	v16.2d,v16.2d,v12.2d	\n\t"\
		"fsub	v0.2d,v2.2d,v5.2d	\n\t fsub	v17.2d,v17.2d,v13.2d	\n\t"\
		"fadd	v6.2d,v6.2d,v10.2d	\n\t fsub	v20.2d,v22.2d,v18.2d	\n\t"\
		"fadd	v7.2d,v7.2d,v11.2d	\n\t fsub	v21.2d,v23.2d,v19.2d	\n\t"\
		"fadd	v4.2d,v4.2d,v3.2d	\n\t fsub	v13.2d,v15.2d,v16.2d	\n\t"\
		"fadd	v5.2d,v5.2d,v2.2d	\n\t fsub	v12.2d,v14.2d,v17.2d	\n\t"\
		"stp	q6,q7,[x0]			\n\t fadd	v18.2d,v18.2d,v22.2d	\n\t"\
		"stp	q0,q4,[x1]			\n\t fadd	v19.2d,v19.2d,v23.2d	\n\t"\
		"stp	q8,q9,[x2]			\n\t fadd	v16.2d,v16.2d,v15.2d	\n\t"\
		"stp	q5,q1,[x3]			\n\t fadd	v17.2d,v17.2d,v14.2d	\n\t"\
		"										stp	q18,q19,[x10]		\n\t"/* out 0 */\
		"										stp	q12,q16,[x11]		\n\t"/* out 1 */\
		"										stp	q20,q21,[x12]		\n\t"/* out 2 */\
		"										stp	q17,q13,[x13]		\n\t"/* out 3 */\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r0] "m" (Xr0)\
		 ,[__two] "m" (Xtwo)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "m" (Xpfetch_dist)\
		: "cc","memory","x0","x1","x2","x3","x4","x5","x6","v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11",\
		"x10","x11","x12","x13","x14","x15","v12","v13","v14","v15","v16","v17","v18","v19","v20","v21","v22","v23"	/* Clobbered registers */\
	);\
	}

	// Vector-opcount: 56 LDP, 32 STP, 128 ADD, 96 MUL; 6 more loads than DIF
	#define SSE2_RADIX16_DIT_TWIDDLE_V2(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr0,Xtwo,Xcc0,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/*
		// Pass 1:
		j = p4*8;	// Set j to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k = 0x80;	// Set k to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		// 0x60 = literal-bytewise address offset between the pointers to the first and second DFT's twiddles:
		i0 = (vec_dbl *)a; i1 = (vec_dbl *)(a+p1); i2 = (vec_dbl *)(a+p2); i3 = (vec_dbl *)(a+p3);
		o0 = r00; o1 = r00+2; o2 = r00+4; o3 = r00+6;
		c_tmp = cc0+6;	// c2,1,3
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,0x60, o0,o1,o2,o3,k)
	*/\
		"ldr	x4,%[__cc0]			\n\t	ldr	x14,%[__r0]				\n\t	add	x4, x4,0x60			\n\t"/* cc0 + 6 */\
		"ldr	x0,%[__add0]		\n\t	ldr	w5,%[__p4]				\n\t	ldr	w6,%[__p8]			\n\t"\
		"ldr	w1,%[__p1]			\n\t	add	x1 , x0,x1,lsl #3		\n\t"/* add0 + p1 */\
		"ldr	w2,%[__p2]			\n\t	add	x2 , x0,x2,lsl #3		\n\t"/* add0 + p2 */\
		"ldr	w3,%[__p3]			\n\t	add	x3 , x0,x3,lsl #3		\n\t"/* add0 + p3 */\
		"ldp q0,q1,[x0]	\n\t"/* Ar,i0 */"	add	x10, x0,x5,lsl #3		\n\t"/* add0 + p4 */\
		"ldp q2,q3,[x1]	\n\t"/* Ar,i1 */"	add	x11, x1,x5,lsl #3		\n\t"/* add0 + p5 */\
		"ldp q4,q5,[x2]	\n\t"/* Ar,i2 */"	add	x12, x2,x5,lsl #3		\n\t"/* add0 + p6 */\
		"ldp q6,q7,[x3]	\n\t"/* Ar,i3 */"	add	x13, x3,x5,lsl #3		\n\t"/* add0 + p7 */\
		"fsub	v8.2d, v0.2d,v2.2d	\n\t	ldp	q10,q11,[x10]		\n\t"/* Ar,i4 */\
		"fsub	v9.2d, v1.2d,v3.2d	\n\t	ldp	q12,q13,[x11]		\n\t"/* Ar,i5 */\
		"fsub	v20.2d,v4.2d,v6.2d	\n\t	ldp	q14,q15,[x12]		\n\t"/* Ar,i6 */\
		"fsub	v21.2d,v5.2d,v7.2d	\n\t	ldp	q16,q17,[x13]		\n\t"/* Ar,i7 */\
		"fadd	v2.2d,v0.2d,v2.2d	\n\t	fsub	v18.2d,v10.2d,v12.2d	\n\t"\
		"fadd	v3.2d,v1.2d,v3.2d	\n\t	fsub	v19.2d,v11.2d,v13.2d	\n\t"\
		"fadd	v6.2d,v4.2d,v6.2d	\n\t	fsub	v22.2d,v14.2d,v16.2d	\n\t"\
		"fadd	v7.2d,v5.2d,v7.2d	\n\t	fsub	v23.2d,v15.2d,v17.2d	\n\t"\
		"fsub	v4.2d,v2.2d,v6.2d	\n\t	fadd	v12.2d,v10.2d,v12.2d	\n\t"\
		"fsub	v5.2d,v3.2d,v7.2d	\n\t	fadd	v13.2d,v11.2d,v13.2d	\n\t"\
		"fsub	v0.2d,v8.2d,v21.2d	\n\t	fadd	v16.2d,v14.2d,v16.2d	\n\t"\
		"fsub	v1.2d,v9.2d,v20.2d	\n\t	fadd	v17.2d,v15.2d,v17.2d	\n\t"\
		"fadd	v6.2d,v2.2d,v6.2d	\n\t	fsub	v14.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v7.2d,v3.2d,v7.2d	\n\t	fsub	v15.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v3.2d,v8.2d,v21.2d	\n\t	fsub	v10.2d,v18.2d,v23.2d	\n\t"\
		"fadd	v2.2d,v9.2d,v20.2d	\n\t	fsub	v11.2d,v19.2d,v22.2d	\n\t"\
		"stp	q6,q7,[x14      ]	\n\t	fadd	v16.2d,v12.2d,v16.2d	\n\t"\
		"ldp	q8,q9,[x4]			\n\t	fadd	v17.2d,v13.2d,v17.2d	\n\t"\
		"fmul	v6.2d,v4.2d,v8.2d	\n\t	fadd	v13.2d,v18.2d,v23.2d	\n\t"\
		"fmul	v7.2d,v5.2d,v8.2d	\n\t	fadd	v12.2d,v19.2d,v22.2d	\n\t"\
		"fmla	v6.2d,v5.2d,v9.2d	\n\t	stp	q16,q17,[x14,#0x80]	\n\t"\
		"fmls	v7.2d,v4.2d,v9.2d	\n\t	ldp	q18,q19,[x4,#0x60]	\n\t"\
		"stp	q6,q7,[x14,#0x40]	\n\t	fmul	v16.2d,v14.2d,v18.2d	\n\t"\
		"ldp	q8,q9,[x4,#0x20]	\n\t	fmul	v17.2d,v15.2d,v18.2d	\n\t"\
		"fmul	v6.2d,v3.2d,v8.2d	\n\t	fmla	v16.2d,v15.2d,v19.2d	\n\t"\
		"fmul	v7.2d,v1.2d,v8.2d	\n\t	fmls	v17.2d,v14.2d,v19.2d	\n\t"\
		"fmla	v6.2d,v1.2d,v9.2d	\n\t	stp	q16,q17,[x14,#0xc0]	\n\t"\
		"fmls	v7.2d,v3.2d,v9.2d	\n\t	ldp	q18,q19,[x4,#0x80]	\n\t"\
		"stp	q6,q7,[x14,#0x20]	\n\t	fmul	v16.2d,v13.2d,v18.2d	\n\t"\
		"ldp	q8,q9,[x4,#0x40]	\n\t	fmul	v17.2d,v11.2d,v18.2d	\n\t"\
		"fmul	v6.2d,v0.2d,v8.2d	\n\t	fmla	v16.2d,v11.2d,v19.2d	\n\t"\
		"fmul	v7.2d,v2.2d,v8.2d	\n\t	fmls	v17.2d,v13.2d,v19.2d	\n\t"\
		"fmla	v6.2d,v2.2d,v9.2d	\n\t	stp	q16,q17,[x14,#0xa0]	\n\t"\
		"fmls	v7.2d,v0.2d,v9.2d	\n\t	ldp	q18,q19,[x4,#0xa0]	\n\t"\
		"stp	q6,q7,[x14,#0x60]	\n\t	fmul	v16.2d,v10.2d,v18.2d	\n\t"\
		"add	x0, x0,x6,lsl #3	\n\t	fmul	v17.2d,v12.2d,v18.2d	\n\t"\
		"add	x1, x1,x6,lsl #3	\n\t	fmla	v16.2d,v12.2d,v19.2d	\n\t"\
		"add	x2, x2,x6,lsl #3	\n\t	fmls	v17.2d,v10.2d,v19.2d	\n\t"\
		"add	x3, x3,x6,lsl #3	\n\t	stp	q16,q17,[x14,#0xe0]	\n\t"\
	/*
		i0 = (vec_dbl *)(a+p8); i1 = (vec_dbl *)(a+p9); i2 = (vec_dbl *)(a+pA); i3 = (vec_dbl *)(a+pB);
		o0 += 16; o1 += 16; o2 += 16; o3 += 16;
		c_tmp += 12;		// cA,9,B
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,0x60, o0,o1,o2,o3,k)
	*/\
		"add	x14, x14,#0x100		\n\t	add	x4, x4,0xc0				\n\t"/* r0 += 16, cc0+6 += 12 */\
		"ldp q0,q1,[x0]	\n\t"/* Ar,i0 */"	add	x10,x10,x6,lsl #3		\n\t"/* x0,10 = add0 + p8 ,12 */\
		"ldp q2,q3,[x1]	\n\t"/* Ar,i1 */"	add	x11,x11,x6,lsl #3		\n\t"/* x1,11 = add0 + p9 ,13 */\
		"ldp q4,q5,[x2]	\n\t"/* Ar,i2 */"	add	x12,x12,x6,lsl #3		\n\t"/* x2,12 = add0 + p10,14 */\
		"ldp q6,q7,[x3]	\n\t"/* Ar,i3 */"	add	x13,x13,x6,lsl #3		\n\t"/* x3,13 = add0 + p11,15 */\
		"fsub	v8.2d, v0.2d,v2.2d	\n\t	ldp	q10,q11,[x10]		\n\t"/* Ar,i4 */\
		"fsub	v9.2d, v1.2d,v3.2d	\n\t	ldp	q12,q13,[x11]		\n\t"/* Ar,i5 */\
		"fsub	v20.2d,v4.2d,v6.2d	\n\t	ldp	q14,q15,[x12]		\n\t"/* Ar,i6 */\
		"fsub	v21.2d,v5.2d,v7.2d	\n\t	ldp	q16,q17,[x13]		\n\t"/* Ar,i7 */\
		"fadd	v2.2d,v0.2d,v2.2d	\n\t	fsub	v18.2d,v10.2d,v12.2d	\n\t"\
		"fadd	v3.2d,v1.2d,v3.2d	\n\t	fsub	v19.2d,v11.2d,v13.2d	\n\t"\
		"fadd	v6.2d,v4.2d,v6.2d	\n\t	fsub	v22.2d,v14.2d,v16.2d	\n\t"\
		"fadd	v7.2d,v5.2d,v7.2d	\n\t	fsub	v23.2d,v15.2d,v17.2d	\n\t"\
		"fsub	v4.2d,v2.2d,v6.2d	\n\t	fadd	v12.2d,v10.2d,v12.2d	\n\t"\
		"fsub	v5.2d,v3.2d,v7.2d	\n\t	fadd	v13.2d,v11.2d,v13.2d	\n\t"\
		"fsub	v0.2d,v8.2d,v21.2d	\n\t	fadd	v16.2d,v14.2d,v16.2d	\n\t"\
		"fsub	v1.2d,v9.2d,v20.2d	\n\t	fadd	v17.2d,v15.2d,v17.2d	\n\t"\
		"fadd	v6.2d,v2.2d,v6.2d	\n\t	fsub	v14.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v7.2d,v3.2d,v7.2d	\n\t	fsub	v15.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v3.2d,v8.2d,v21.2d	\n\t	fsub	v10.2d,v18.2d,v23.2d	\n\t"\
		"fadd	v2.2d,v9.2d,v20.2d	\n\t	fsub	v11.2d,v19.2d,v22.2d	\n\t"\
		"stp	q6,q7,[x14      ]	\n\t	fadd	v16.2d,v12.2d,v16.2d	\n\t"\
		"ldp	q8,q9,[x4]			\n\t	fadd	v17.2d,v13.2d,v17.2d	\n\t"\
		"fmul	v6.2d,v4.2d,v8.2d	\n\t	fadd	v13.2d,v18.2d,v23.2d	\n\t"\
		"fmul	v7.2d,v5.2d,v8.2d	\n\t	fadd	v12.2d,v19.2d,v22.2d	\n\t"\
		"fmla	v6.2d,v5.2d,v9.2d	\n\t	stp	q16,q17,[x14,#0x80]	\n\t"\
		"fmls	v7.2d,v4.2d,v9.2d	\n\t	ldp	q18,q19,[x4,#0x60]	\n\t"\
		"stp	q6,q7,[x14,#0x40]	\n\t	fmul	v16.2d,v14.2d,v18.2d	\n\t"\
		"ldp	q8,q9,[x4,#0x20]	\n\t	fmul	v17.2d,v15.2d,v18.2d	\n\t"\
		"fmul	v6.2d,v3.2d,v8.2d	\n\t	fmla	v16.2d,v15.2d,v19.2d	\n\t"\
		"fmul	v7.2d,v1.2d,v8.2d	\n\t	fmls	v17.2d,v14.2d,v19.2d	\n\t"\
		"fmla	v6.2d,v1.2d,v9.2d	\n\t	stp	q16,q17,[x14,#0xc0]	\n\t"\
		"fmls	v7.2d,v3.2d,v9.2d	\n\t	ldp	q18,q19,[x4,#0x80]	\n\t"\
		"stp	q6,q7,[x14,#0x20]	\n\t	fmul	v16.2d,v13.2d,v18.2d	\n\t"\
		"ldp	q8,q9,[x4,#0x40]	\n\t	fmul	v17.2d,v11.2d,v18.2d	\n\t"\
		"fmul	v6.2d,v0.2d,v8.2d	\n\t	fmla	v16.2d,v11.2d,v19.2d	\n\t"\
		"fmul	v7.2d,v2.2d,v8.2d	\n\t	fmls	v17.2d,v13.2d,v19.2d	\n\t"\
		"fmla	v6.2d,v2.2d,v9.2d	\n\t	stp	q16,q17,[x14,#0xa0]	\n\t"\
		"fmls	v7.2d,v0.2d,v9.2d	\n\t	ldp	q18,q19,[x4,#0xa0]	\n\t"\
		"stp	q6,q7,[x14,#0x60]	\n\t			fmul	v16.2d,v10.2d,v18.2d	\n\t"\
		"mov x2, x0			\n\t"/* x2=add0+p8  */"	fmul	v17.2d,v12.2d,v18.2d	\n\t"\
		"mov x3,x10			\n\t"/* x3=add0+p12 */"	fmla	v16.2d,v12.2d,v19.2d	\n\t"\
		"sub x0,x0,x6,lsl #3\n\t"/* x0=add0     */"	fmls	v17.2d,v10.2d,v19.2d	\n\t"\
		"add x1,x0,x5,lsl #3\n\t"/* x1=add0+p4  */"	stp	q16,q17,[x14,#0xe0]	\n\t"\
	/*
		// Pass 2:
		j = 0x40;	// Set j to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k = p2*8;	// Set k to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		c_tmp = cc0;	// c8,4,C
		// Pass1-DFTs all use same twiddle-triplet, so lit-byte address offset between ptrs to 1st,2nd DFT's twiddles = 0
		i0 = r00; i1 = r00+8; i2 = r00+16; i3 = r00+24;
		o0 = (vec_dbl *)a; o1 = (vec_dbl *)(a+p4); o2 = (vec_dbl *)(a+p8); o3 = (vec_dbl *)(a+pC);
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,   0, o0,o1,o2,o3,k)
	*/\
		"ldr	x4,%[__cc0]			\n\t	ldr	x14,%[__r0]					\n\t"\
		"ldr	w5,%[__p1]			\n\t	ldr	w6 ,%[__p2]					\n\t"\
	"ldp q0,q1,[x14       ]\n\t"/* Ar,i0 */"add	x10, x0,x6,lsl #3			\n\t"/* add0 + p2  */\
	"ldp q2,q3,[x14,#0x080]\n\t"/* Ar,i1 */"add	x11, x1,x6,lsl #3			\n\t"/* add0 + p6  */\
	"ldp q4,q5,[x14,#0x100]\n\t"/* Ar,i2 */"add	x12, x2,x6,lsl #3			\n\t"/* add0 + p10 */\
	"ldp q6,q7,[x14,#0x180]\n\t"/* Ar,i3 */"add	x13, x3,x6,lsl #3			\n\t"/* add0 + p14 */\
		"fsub	v8.2d, v0.2d,v2.2d	\n\t	ldp	q10,q11,[x14,#0x040]		\n\t"/* Ar,i4 */\
		"fsub	v9.2d, v1.2d,v3.2d	\n\t	ldp	q12,q13,[x14,#0x0c0]		\n\t"/* Ar,i5 */\
		"fsub	v20.2d,v4.2d,v6.2d	\n\t	ldp	q14,q15,[x14,#0x140]		\n\t"/* Ar,i6 */\
		"fsub	v21.2d,v5.2d,v7.2d	\n\t	ldp	q16,q17,[x14,#0x1c0]		\n\t"/* Ar,i7 */\
		"fadd	v2.2d,v0.2d,v2.2d	\n\t	fsub	v18.2d,v10.2d,v12.2d	\n\t"\
		"fadd	v3.2d,v1.2d,v3.2d	\n\t	fsub	v19.2d,v11.2d,v13.2d	\n\t"\
		"fadd	v6.2d,v4.2d,v6.2d	\n\t	fsub	v22.2d,v14.2d,v16.2d	\n\t"\
		"fadd	v7.2d,v5.2d,v7.2d	\n\t	fsub	v23.2d,v15.2d,v17.2d	\n\t"\
		"fsub	v4.2d,v2.2d,v6.2d	\n\t	fadd	v12.2d,v10.2d,v12.2d	\n\t"\
		"fsub	v5.2d,v3.2d,v7.2d	\n\t	fadd	v13.2d,v11.2d,v13.2d	\n\t"\
		"fsub	v0.2d,v8.2d,v21.2d	\n\t	fadd	v16.2d,v14.2d,v16.2d	\n\t"\
		"fsub	v1.2d,v9.2d,v20.2d	\n\t	fadd	v17.2d,v15.2d,v17.2d	\n\t"\
		"fadd	v6.2d,v2.2d,v6.2d	\n\t	fsub	v14.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v7.2d,v3.2d,v7.2d	\n\t	fsub	v15.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v3.2d,v8.2d,v21.2d	\n\t	fsub	v10.2d,v18.2d,v23.2d	\n\t"\
		"fadd	v2.2d,v9.2d,v20.2d	\n\t	fsub	v11.2d,v19.2d,v22.2d	\n\t"\
		"stp	q6,q7,[x0]			\n\t	fadd	v16.2d,v12.2d,v16.2d	\n\t"\
		"ldp	q8,q9,[x4      ]	\n\t	fadd	v17.2d,v13.2d,v17.2d	\n\t"\
		"fmul	v6.2d,v4.2d,v8.2d	\n\t	fadd	v13.2d,v18.2d,v23.2d	\n\t"\
		"fmul	v7.2d,v5.2d,v8.2d	\n\t	fadd	v12.2d,v19.2d,v22.2d	\n\t"\
		"fmla	v6.2d,v5.2d,v9.2d	\n\t	stp	q16,q17,[x10]		\n\t"\
		"fmls	v7.2d,v4.2d,v9.2d	\n\t	ldp	q18,q19,[x4      ]	\n\t"\
		"stp	q6,q7,[x2]			\n\t	fmul	v16.2d,v14.2d,v18.2d	\n\t"\
		"ldp	q8,q9,[x4,#0x20]	\n\t	fmul	v17.2d,v15.2d,v18.2d	\n\t"\
		"fmul	v6.2d,v3.2d,v8.2d	\n\t	fmla	v16.2d,v15.2d,v19.2d	\n\t"\
		"fmul	v7.2d,v1.2d,v8.2d	\n\t	fmls	v17.2d,v14.2d,v19.2d	\n\t"\
		"fmla	v6.2d,v1.2d,v9.2d	\n\t	stp	q16,q17,[x12]		\n\t"\
		"fmls	v7.2d,v3.2d,v9.2d	\n\t	ldp	q18,q19,[x4,#0x20]	\n\t"\
		"stp	q6,q7,[x1]			\n\t	fmul	v16.2d,v13.2d,v18.2d	\n\t"\
		"ldp	q8,q9,[x4,#0x40]	\n\t	fmul	v17.2d,v11.2d,v18.2d	\n\t"\
		"fmul	v6.2d,v0.2d,v8.2d	\n\t	fmla	v16.2d,v11.2d,v19.2d	\n\t"\
		"fmul	v7.2d,v2.2d,v8.2d	\n\t	fmls	v17.2d,v13.2d,v19.2d	\n\t"\
		"fmla	v6.2d,v2.2d,v9.2d	\n\t	stp	q16,q17,[x11]		\n\t"\
		"fmls	v7.2d,v0.2d,v9.2d	\n\t	ldp	q18,q19,[x4,#0x40]	\n\t"\
		"stp	q6,q7,[x3]			\n\t	fmul	v16.2d,v10.2d,v18.2d	\n\t"\
		"add	x0, x0,x5,lsl #3	\n\t	fmul	v17.2d,v12.2d,v18.2d	\n\t"\
		"add	x1, x1,x5,lsl #3	\n\t	fmla	v16.2d,v12.2d,v19.2d	\n\t"\
		"add	x2, x2,x5,lsl #3	\n\t	fmls	v17.2d,v10.2d,v19.2d	\n\t"\
		"add	x3, x3,x5,lsl #3	\n\t	stp	q16,q17,[x13]		\n\t"\
	/*
		i0 += 2; i1 += 2; i2 += 2; i3 += 2;
		o0 = (vec_dbl *)(a+p1); o1 = (vec_dbl *)(a+p5); o2 = (vec_dbl *)(a+p9); o3 = (vec_dbl *)(a+pD);
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,   0, o0,o1,o2,o3,k)
	*/\
		"add	x14,x14,#0x20		\n\t"/* r0 + 2 */\
	"ldp q0,q1,[x14       ]\n\t"/* Ar,i0 */"add	x10,x10,x5,lsl #3			\n\t"/* x0,10 = add0 + p1 ,3  */\
	"ldp q2,q3,[x14,#0x080]\n\t"/* Ar,i1 */"add	x11,x11,x5,lsl #3			\n\t"/* x1,11 = add0 + p5 ,7  */\
	"ldp q4,q5,[x14,#0x100]\n\t"/* Ar,i2 */"add	x12,x12,x5,lsl #3			\n\t"/* x2,12 = add0 + p9 ,11 */\
	"ldp q6,q7,[x14,#0x180]\n\t"/* Ar,i3 */"add	x13,x13,x5,lsl #3			\n\t"/* x3,13 = add0 + p13,15 */\
		"fsub	v8.2d, v0.2d,v2.2d	\n\t	ldp	q10,q11,[x14,#0x040]		\n\t"/* Ar,i4 */\
		"fsub	v9.2d, v1.2d,v3.2d	\n\t	ldp	q12,q13,[x14,#0x0c0]		\n\t"/* Ar,i5 */\
		"fsub	v20.2d,v4.2d,v6.2d	\n\t	ldp	q14,q15,[x14,#0x140]		\n\t"/* Ar,i6 */\
		"fsub	v21.2d,v5.2d,v7.2d	\n\t	ldp	q16,q17,[x14,#0x1c0]		\n\t"/* Ar,i7 */\
		"fadd	v2.2d,v0.2d,v2.2d	\n\t	fsub	v18.2d,v10.2d,v12.2d	\n\t"\
		"fadd	v3.2d,v1.2d,v3.2d	\n\t	fsub	v19.2d,v11.2d,v13.2d	\n\t"\
		"fadd	v6.2d,v4.2d,v6.2d	\n\t	fsub	v22.2d,v14.2d,v16.2d	\n\t"\
		"fadd	v7.2d,v5.2d,v7.2d	\n\t	fsub	v23.2d,v15.2d,v17.2d	\n\t"\
		"fsub	v4.2d,v2.2d,v6.2d	\n\t	fadd	v12.2d,v10.2d,v12.2d	\n\t"\
		"fsub	v5.2d,v3.2d,v7.2d	\n\t	fadd	v13.2d,v11.2d,v13.2d	\n\t"\
		"fsub	v0.2d,v8.2d,v21.2d	\n\t	fadd	v16.2d,v14.2d,v16.2d	\n\t"\
		"fsub	v1.2d,v9.2d,v20.2d	\n\t	fadd	v17.2d,v15.2d,v17.2d	\n\t"\
		"fadd	v6.2d,v2.2d,v6.2d	\n\t	fsub	v14.2d,v12.2d,v16.2d	\n\t"\
		"fadd	v7.2d,v3.2d,v7.2d	\n\t	fsub	v15.2d,v13.2d,v17.2d	\n\t"\
		"fadd	v3.2d,v8.2d,v21.2d	\n\t	fsub	v10.2d,v18.2d,v23.2d	\n\t"\
		"fadd	v2.2d,v9.2d,v20.2d	\n\t	fsub	v11.2d,v19.2d,v22.2d	\n\t"\
		"stp	q6,q7,[x0]			\n\t	fadd	v16.2d,v12.2d,v16.2d	\n\t"\
		"ldp	q8,q9,[x4      ]	\n\t	fadd	v17.2d,v13.2d,v17.2d	\n\t"\
		"fmul	v6.2d,v4.2d,v8.2d	\n\t	fadd	v13.2d,v18.2d,v23.2d	\n\t"\
		"fmul	v7.2d,v5.2d,v8.2d	\n\t	fadd	v12.2d,v19.2d,v22.2d	\n\t"\
		"fmla	v6.2d,v5.2d,v9.2d	\n\t	stp	q16,q17,[x10]		\n\t"\
		"fmls	v7.2d,v4.2d,v9.2d	\n\t	ldp	q18,q19,[x4      ]	\n\t"\
		"stp	q6,q7,[x2]			\n\t	fmul	v16.2d,v14.2d,v18.2d	\n\t"\
		"ldp	q8,q9,[x4,#0x20]	\n\t	fmul	v17.2d,v15.2d,v18.2d	\n\t"\
		"fmul	v6.2d,v3.2d,v8.2d	\n\t	fmla	v16.2d,v15.2d,v19.2d	\n\t"\
		"fmul	v7.2d,v1.2d,v8.2d	\n\t	fmls	v17.2d,v14.2d,v19.2d	\n\t"\
		"fmla	v6.2d,v1.2d,v9.2d	\n\t	stp	q16,q17,[x12]		\n\t"\
		"fmls	v7.2d,v3.2d,v9.2d	\n\t	ldp	q18,q19,[x4,#0x20]	\n\t"\
		"stp	q6,q7,[x1]			\n\t	fmul	v16.2d,v13.2d,v18.2d	\n\t"\
		"ldp	q8,q9,[x4,#0x40]	\n\t	fmul	v17.2d,v11.2d,v18.2d	\n\t"\
		"fmul	v6.2d,v0.2d,v8.2d	\n\t	fmla	v16.2d,v11.2d,v19.2d	\n\t"\
		"fmul	v7.2d,v2.2d,v8.2d	\n\t	fmls	v17.2d,v13.2d,v19.2d	\n\t"\
		"fmla	v6.2d,v2.2d,v9.2d	\n\t	stp	q16,q17,[x11]		\n\t"\
		"fmls	v7.2d,v0.2d,v9.2d	\n\t	ldp	q18,q19,[x4,#0x40]	\n\t"\
		"stp	q6,q7,[x3]			\n\t	fmul	v16.2d,v10.2d,v18.2d	\n\t"\
		"									fmul	v17.2d,v12.2d,v18.2d	\n\t"\
		"									fmla	v16.2d,v12.2d,v19.2d	\n\t"\
		"									fmls	v17.2d,v10.2d,v19.2d	\n\t"\
		"									stp	q16,q17,[x13]		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r0] "m" (Xr0)\
		 ,[__two] "m" (Xtwo)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "m" (Xpfetch_dist)\
		: "cc","memory","x0","x1","x2","x3","x4","x5","x6","v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11",\
		"x10","x11","x12","x13","x14","v12","v13","v14","v15","v16","v17","v18","v19","v20","v21","v22","v23"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_AVX512)	// See AVX2 commentary for mem-layout details

	/* In avx512 can use VRCP14PD to compute 8 DP-approx-recips good to 14 bits, no need for
	prior double -> single conversion, use pair of 2nd-order Newton updates rather than 2nd/3rd-order combo.
	Further, the resulting mults-preprocessing should be fast enough to seriously consider doing the dyadic-
	square-step-wrapping radix-16 DFTs in this fashion, as well.

	Compare arithmetic cost of the 16-double iterative inversion used to get the tangents:
		AVX2   Cost: 4 vcvtpd2ps, 4 vrcpps  , 4 vcvtps2pd, 12 fma, 4 sub, 8 mul
		AVX512 Cost:              2 vrcp14pd,               4 fma,        4 mul
	*/
	#define RADIX16_COMPUTE_FMA_SINCOS_DIF(Xadd0,Xone,Xcos,Xtan)\
	{\
	__asm__ volatile (\
	/*************** Keep opening section in same form as AVX2, not much to be gained from 8-doubles-wide **************/\
		"movq	%[__add0],%%rax			\n\t"/* base-ptr to 32-byte-aligned 128-byte-large memchunk 16 doubles needing iterative inversion */\
		"movq	%[__one] ,%%rbx			\n\t"/* 1.0 in 4-fold-double form */\
		"vmovaps	    (%%rbx),%%ymm14	\n\t"/* 1.0 */\
		"vaddpd	%%ymm14,%%ymm14,%%ymm14	\n\t"/* 2.0 */\
		"vmovapd	    (%%rax),%%ymm4	\n\t"/* load inputs-to-be-inverted (d's) into 4 AVX registers */\
		"vmovapd	0x20(%%rax),%%ymm5	\n\t"\
		"vmovapd	0x40(%%rax),%%ymm6	\n\t"\
		"vmovapd	0x60(%%rax),%%ymm7	\n\t"\
		"vrcp14pd	%%zmm4,%%zmm0	\n\t"	/* ainv := approx 1/d to >= 14 bits of precision */\
		"vrcp14pd	%%zmm5,%%zmm1	\n\t"	/* NOTE: AVX-512F requires full-width zmm-regs here, but only use lower 256 bits */\
		"vrcp14pd	%%zmm6,%%zmm2	\n\t"\
		"vrcp14pd	%%zmm7,%%zmm3	\n\t"\
		/* 1st NR iteration gives ~28 bits of precision: */\
		"vmovaps	%%ymm0,%%ymm8 	\n\t"	/* make a copy of ainv */\
		"vmovaps	%%ymm1,%%ymm9 	\n\t"\
		"vmovaps	%%ymm2,%%ymm10	\n\t"\
		"vmovaps	%%ymm3,%%ymm11	\n\t"\
		"vfnmadd132pd	%%ymm4,%%ymm14,%%ymm0 	\n\t"	/* 2 - d*ainv, overwrites ainv */\
		"vfnmadd132pd	%%ymm5,%%ymm14,%%ymm1 	\n\t"\
		"vfnmadd132pd	%%ymm6,%%ymm14,%%ymm2 	\n\t"\
		"vfnmadd132pd	%%ymm7,%%ymm14,%%ymm3 	\n\t"\
		"vmulpd		%%ymm0,%%ymm8 ,%%ymm0	\n\t"	/* ainv*(2 - d*ainv) = 1/d accurate to ~28 bits */\
		"vmulpd		%%ymm1,%%ymm9 ,%%ymm1	\n\t"\
		"vmulpd		%%ymm2,%%ymm10,%%ymm2	\n\t"\
		"vmulpd		%%ymm3,%%ymm11,%%ymm3	\n\t"\
		/* 2nd NR iteration gives maximal 53 bits of precision: */\
		"vmovaps	%%ymm0,%%ymm8 	\n\t"	/* make a copy of ainv */\
		"vmovaps	%%ymm1,%%ymm9 	\n\t"\
		"vmovaps	%%ymm2,%%ymm10	\n\t"\
		"vmovaps	%%ymm3,%%ymm11	\n\t"\
		"vfnmadd132pd	%%ymm4,%%ymm14,%%ymm0 	\n\t"	/* 2 - d*ainv, overwrites ainv */\
		"vfnmadd132pd	%%ymm5,%%ymm14,%%ymm1 	\n\t"\
		"vfnmadd132pd	%%ymm6,%%ymm14,%%ymm2 	\n\t"\
		"vfnmadd132pd	%%ymm7,%%ymm14,%%ymm3 	\n\t"\
		"vmulpd		%%ymm0,%%ymm8 ,%%ymm0	\n\t"	/* ainv*(2 - d*ainv) = 1/d accurate to ~53 bits */\
		"vmulpd		%%ymm1,%%ymm9 ,%%ymm1	\n\t"\
		"vmulpd		%%ymm2,%%ymm10,%%ymm2	\n\t"\
		"vmulpd		%%ymm3,%%ymm11,%%ymm3	\n\t"\
		/* Multiply by the sine terms to get quotients in ymm4-7 (still need inverses in ymm0-3 for cosine-ratios): */\
		"vmulpd	0x80(%%rax),%%ymm0,%%ymm4	\n\t"\
		"vmulpd	0xa0(%%rax),%%ymm1,%%ymm5	\n\t"\
		"vmulpd	0xc0(%%rax),%%ymm2,%%ymm6	\n\t"\
		"vmulpd	0xe0(%%rax),%%ymm3,%%ymm7	\n\t"\
		/* Outputs overwrite sines: */\
		"vmovaps	%%ymm4,0x80(%%rax)	\n\t"\
		"vmovaps	%%ymm5,0xa0(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0xc0(%%rax)	\n\t"\
		"vmovaps	%%ymm7,0xe0(%%rax)	\n\t"\
		/* Set up for cosine ratios: */\
		"vpermilpd	$11,%%ymm0,%%ymm2 	\n\t"/* permute ymm0 = 1/c[3123] to get ymm2 = 1/c[1123], then *= [c3,c5,c6,0] */\
		"vmulpd	0x100(%%rax),%%ymm2,%%ymm2	\n\t"/* ymm2 = [c3,c5,c6, 0] * 1/c[1123] = [c31,c51,c62,  0] */\
		"vmulpd	0x120(%%rax),%%ymm0,%%ymm0	\n\t"/* ymm0 = [c7,c9,cA,cB] * 1/c[3123] = [c73,c91,cA2,cB3] */\
		"vmulpd	0x140(%%rax),%%ymm1,%%ymm1	\n\t"/* ymm1 = [cC,cD,cE,cF] * 1/c[4567] = [cC4,cD5,cE6,cF7] */\
		/* Outputs into slots imm. above above tangents: */\
		"vmovaps	%%ymm2,0x100(%%rax)	\n\t"\
		"vmovaps	%%ymm0,0x120(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x140(%%rax)	\n\t"\
		/* Now propagate the scalar doubles in the 11 ymm-sized memlocs referenced above to their final 4x-copied locs. Scalar data laid out as */\
		/* add0 + 0x[  0,  8, 10, 18]: [c3 ,c1 ,c2 ,c3 ] Cosines: */\
		/* add0 + 0x[ 20, 28, 30, 38]: [c4 ,c5 ,c6 ,c7 ] */\
		/* add0 + 0x[ 40, 48, 50, 58]: [c8 ,c9 ,cA ,cB ] */\
		/* add0 + 0x[ 60, 68, 70, 78]: [cC ,cD ,cE ,cF ] */\
		/* add0 + 0x[ 80, 88, 90, 98]: [-- ,r1 ,r2 ,r3 ] Tangents: */\
		/* add0 + 0x[ a0, a8, b0, b8]: [r4 ,r5 ,r6 ,r7 ] */\
		/* add0 + 0x[ c0, c8, d0, d8]: [r8 ,r9 ,rA ,rB ] */\
		/* add0 + 0x[ e0, e8, f0, f8]: [rC ,rD ,rE ,rF ] */\
		/* add0 + 0x[100,108,110,118]: [c31,c51,c62,  0] Cosine ratios: */\
		/* add0 + 0x[120,128,130,138]: [c73,c91,cA2,cB3] */\
		/* add0 + 0x[140,148,150,158]: [cC4,cD5,cE6,cF7] */\
		/* We start at the upper end of the memory chunk allocated for sincos data and proceed downward, so the 11 memory */\
		/* slots inited with scalar data above are only overwritten after their scalar-data contents have been copied: */\
		"vbroadcastsd	0x0f8(%%rax),%%zmm15	/* sF*/	\n\t"/* Allow several cycles for each vbroadcast to complete before writing result */\
		"vbroadcastsd	0x158(%%rax),%%zmm14	/* cF*/	\n\t	movq	%[__cos] ,%%rbx			\n\t"\
		"vbroadcastsd	0x0b8(%%rax),%%zmm13	/* s7*/	\n\t	movq	%[__tan] ,%%rcx			\n\t"\
		"vbroadcastsd	0x120(%%rax),%%zmm12	/* c7*/	\n\t"\
		"vbroadcastsd	0x0d8(%%rax),%%zmm11	/* sB*/	\n\t	vmovaps	%%zmm15,0x840(%%rax)	\n\t"/* s15 = cc0 + 0x21;	__rF = s15/c15		*/\
		"vbroadcastsd	0x138(%%rax),%%zmm10	/* cB*/	\n\t	vmovaps	%%zmm14,0x800(%%rax)	\n\t"/* c15 = cc0 + 0x20;	__cF7 = __cF/__c7	*/\
		"vbroadcastsd	0x098(%%rax),%%zmm9 	/* s3*/	\n\t	vmovaps	%%zmm13,0x7c0(%%rax)	\n\t"/* s7  = cc0 + 0x1f;	__r7 = s7 /c7 		*/\
		"vbroadcastsd	0x100(%%rax),%%zmm8 	/* c3*/	\n\t	vmovaps	%%zmm12,0x780(%%rax)	\n\t"/* c7  = cc0 + 0x1e;	__c73 = __c7/__c3	*/\
		"vbroadcastsd	0x0e8(%%rax),%%zmm7 	/* sD*/	\n\t	vmovaps	%%zmm11,0x740(%%rax)	\n\t"/* s11 = cc0 + 0x1d;	__rB = s11/c11		*/\
		"vbroadcastsd	0x148(%%rax),%%zmm6 	/* cD*/	\n\t	vmovaps	%%zmm10,0x700(%%rax)	\n\t"/* c11 = cc0 + 0x1c;	__cB3 = __cB/__c3	*/\
		"vbroadcastsd	0x0a8(%%rax),%%zmm5 	/* s5*/	\n\t	vmovaps	%%zmm9 ,0x6c0(%%rax)	\n\t"/* s3  = cc0 + 0x1b;	__r3 = s3 /c3 		*/\
		"vbroadcastsd	0x108(%%rax),%%zmm4 	/* c5*/	\n\t	vmovaps	%%zmm8 ,0x680(%%rax)	\n\t"/* c3  = cc0 + 0x1a;	__c31 = __c3/__c1	*/\
		"vbroadcastsd	0x0c8(%%rax),%%zmm3 	/* s9*/	\n\t	vmovaps	%%zmm7 ,0x640(%%rax)	\n\t"/* s13 = cc0 + 0x19;	__rD = s13/c13		*/\
		"vbroadcastsd	0x128(%%rax),%%zmm2 	/* c9*/	\n\t	vmovaps	%%zmm6 ,0x600(%%rax)	\n\t"/* c13 = cc0 + 0x18;	__cD5 = __cD/__c5	*/\
		"vbroadcastsd	0x088(%%rax),%%zmm1 	/* s1*/	\n\t	vmovaps	%%zmm5 ,0x5c0(%%rax)	\n\t"/* s5  = cc0 + 0x17;	__r5 = s5 /c5 		*/\
		"vbroadcastsd	0x0f0(%%rax),%%zmm0 	/* sE*/	\n\t	vmovaps	%%zmm4 ,0x580(%%rax)	\n\t"/* c5  = cc0 + 0x16;	__c51 = __c5/__c1	*/\
		"vbroadcastsd	0x150(%%rax),%%zmm15	/* cE*/	\n\t	vmovaps	%%zmm3 ,0x540(%%rax)	\n\t"/* s9  = cc0 + 0x15;	__r9 = s9 /c9 		*/\
		"vbroadcastsd	0x0b0(%%rax),%%zmm14	/* s6*/	\n\t	vmovaps	%%zmm2 ,0x500(%%rax)	\n\t"/* c9  = cc0 + 0x14;	__c91 = __c9/__c1	*/\
		"vbroadcastsd	0x110(%%rax),%%zmm13	/* c6*/	\n\t	vmovaps	%%zmm1 ,0x4c0(%%rax)	\n\t"/* s1  = cc0 + 0x13;	__r1 = s1 /c1 		*/\
		"vbroadcastsd	0x0d0(%%rax),%%zmm12	/* sA*/	\n\t	vmovaps	%%zmm0 ,0x440(%%rax)	\n\t"/* s14 = cc0 + 0x11;	__rE = s14/c14		*/\
		"vbroadcastsd	0x130(%%rax),%%zmm11	/* cA*/	\n\t	vmovaps	%%zmm15,0x400(%%rax)	\n\t"/* c14 = cc0 + 0x10;	__cE6 = __cE/__c6	*/\
		"vbroadcastsd	0x090(%%rax),%%zmm10	/* s2*/	\n\t	vmovaps	%%zmm14,0x3c0(%%rax)	\n\t"/* s6  = cc0 + 0x0f;	__r6 = s6 /c6 		*/\
		"vbroadcastsd	0x0e0(%%rax),%%zmm9 	/* sC*/	\n\t	vmovaps	%%zmm13,0x380(%%rax)	\n\t"/* c6  = cc0 + 0x0e;	__c62 = __c6/__c2	*/\
		"vbroadcastsd	0x140(%%rax),%%zmm8 	/* cC*/	\n\t	vmovaps	%%zmm12,0x340(%%rax)	\n\t"/* s10 = cc0 + 0x0d;	__rA = s10/c10		*/\
		"vbroadcastsd	0x0a0(%%rax),%%zmm7 	/* s4*/	\n\t	vmovaps	%%zmm11,0x300(%%rax)	\n\t"/* c10 = cc0 + 0x0c;	__cA2 = __cA/__c2	*/\
		"vbroadcastsd	0x020(%%rax),%%zmm6 	/* c4*/	\n\t	vmovaps	%%zmm10,0x2c0(%%rax)	\n\t"/* s2  = cc0 + 0x0b;	__r2 = s2 /c2 		*/\
		"vbroadcastsd	0x0c0(%%rax),%%zmm5 	/* s8*/	\n\t	vmovaps	%%zmm9 ,0x240(%%rax)	\n\t"/* s12 = cc0 + 0x09;	__rC = s12/c12		*/\
		"vbroadcastsd	0x040(%%rax),%%zmm4 	/* c8*/	\n\t	vmovaps	%%zmm8 ,0x200(%%rax)	\n\t"/* c12 = cc0 + 0x08;	__cC4 = __cC/__c4	*/\
		/* Lowest-addressed few data need special handling, and are written in ascending-address order: */\
		"vbroadcastsd	     (%%rbx),%%zmm1	/*__c */	\n\t	vmovaps	%%zmm7 ,0x1c0(%%rax)	\n\t"/* s4  = cc0 + 0x07;	__r4 = s4 /c4 		*/\
		"vbroadcastsd	0x008(%%rax),%%zmm0	/* c1 */	\n\t	vmovaps	%%zmm6 ,0x180(%%rax)	\n\t"/* c4  = cc0 + 0x06;	__c4 [unchanged]	*/\
		"vbroadcastsd	0x010(%%rax),%%zmm3	/* c2 */	\n\t	vmovaps	%%zmm5 ,0x140(%%rax)	\n\t"/* s8  = cc0 + 0x05;	__r8 = s8 /c8 		*/\
		"vmulpd		   %%zmm1,%%zmm0,%%zmm1	/* c1*__c */\n\t	vmovaps	%%zmm4 ,0x100(%%rax)	\n\t"/* c8  = cc0 + 0x04;	__c8 [unchanged]	*/\
		"vbroadcastsd	     (%%rcx),%%zmm2	/*__s/__c */\n\t	vmovaps	%%zmm0 ,0x480(%%rax)	\n\t"/* c1  = cc0 + 0x12;	__c1 [unchanged]	*/\
		"vmulpd	-0x040(%%rax),%%zmm0,%%zmm0	/*c1*ISRT2*/\n\t	vmovaps	%%zmm3 ,0x280(%%rax)	\n\t"/* c2  = cc0 + 0x0a;	__c2 [unchanged]	*/\
		"vmulpd	-0x040(%%rax),%%zmm3,%%zmm3	/*c2*ISRT2*/\n\t	vmovaps	%%zmm1 ,     (%%rax)	\n\t"/* cc0 = cc0 + 0x00;	__c1_c = c1*__c		*/\
		"														vmovaps	%%zmm2 ,0x040(%%rax)	\n\t"/* ss0 = cc0 + 0x01;	__sc = __s/__c		*/\
		"														vmovaps	%%zmm0 ,0x080(%%rax)	\n\t"/* c0  = cc0 + 0x02;	__c1i2 = c1*ISRT2	*/\
		"														vmovaps	%%zmm3 ,0x0c0(%%rax)	\n\t"/* s0  = cc0 + 0x03;	__c2i2 = c2*ISRT2	*/\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__one] "m" (Xone)\
		 ,[__cos] "m" (Xcos)\
		 ,[__tan] "m" (Xtan)\
		: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	// 2nd version of trig-data and DFT-pass routines uses 1-copy trig data, read in and broadcast to full YMM register on-the-fly:
	#define RADIX16_COMPUTE_FMA_SINCOS_DIF_2(Xadd0,Xone)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax			\n\t"/* base-ptr to 32-byte-aligned 128-byte-large memchunk 16 doubles needing iterative inversion */\
		"movq	%[__one] ,%%rbx			\n\t"/* 1.0 in 4-fold-double form */\
		"vmovaps	    (%%rbx),%%ymm14	\n\t"/* 1.0 */\
		"vaddpd	%%ymm14,%%ymm14,%%ymm14	\n\t"/* 2.0 */\
		"vmovapd	    (%%rax),%%ymm4	\n\t"/* load inputs-to-be-inverted (d's) into 4 AVX registers */\
		"vmovapd	0x20(%%rax),%%ymm5	\n\t"\
		"vmovapd	0x40(%%rax),%%ymm6	\n\t"\
		"vmovapd	0x60(%%rax),%%ymm7	\n\t"\
		"vrcp14pd	%%zmm4,%%zmm0	\n\t"	/* ainv := approx 1/d to >= 14 bits of precision */\
		"vrcp14pd	%%zmm5,%%zmm1	\n\t"	/* NOTE: AVX-512F requires full-width zmm-regs here, but only use lower 256 bits */\
		"vrcp14pd	%%zmm6,%%zmm2	\n\t"\
		"vrcp14pd	%%zmm7,%%zmm3	\n\t"\
		/* 1st NR iteration gives ~28 bits of precision: */\
		"vmovaps	%%ymm0,%%ymm8 	\n\t"	/* make a copy of ainv */\
		"vmovaps	%%ymm1,%%ymm9 	\n\t"\
		"vmovaps	%%ymm2,%%ymm10	\n\t"\
		"vmovaps	%%ymm3,%%ymm11	\n\t"\
		"vfnmadd132pd	%%ymm4,%%ymm14,%%ymm0 	\n\t"	/* 2 - d*ainv, overwrites ainv */\
		"vfnmadd132pd	%%ymm5,%%ymm14,%%ymm1 	\n\t"\
		"vfnmadd132pd	%%ymm6,%%ymm14,%%ymm2 	\n\t"\
		"vfnmadd132pd	%%ymm7,%%ymm14,%%ymm3 	\n\t"\
		"vmulpd		%%ymm0,%%ymm8 ,%%ymm0	\n\t"	/* ainv*(2 - d*ainv) = 1/d accurate to ~28 bits */\
		"vmulpd		%%ymm1,%%ymm9 ,%%ymm1	\n\t"\
		"vmulpd		%%ymm2,%%ymm10,%%ymm2	\n\t"\
		"vmulpd		%%ymm3,%%ymm11,%%ymm3	\n\t"\
		/* 2nd NR iteration gives maximal 53 bits of precision: */\
		"vmovaps	%%ymm0,%%ymm8 	\n\t"	/* make a copy of ainv */\
		"vmovaps	%%ymm1,%%ymm9 	\n\t"\
		"vmovaps	%%ymm2,%%ymm10	\n\t"\
		"vmovaps	%%ymm3,%%ymm11	\n\t"\
		"vfnmadd132pd	%%ymm4,%%ymm14,%%ymm0 	\n\t"	/* 2 - d*ainv, overwrites ainv */\
		"vfnmadd132pd	%%ymm5,%%ymm14,%%ymm1 	\n\t"\
		"vfnmadd132pd	%%ymm6,%%ymm14,%%ymm2 	\n\t"\
		"vfnmadd132pd	%%ymm7,%%ymm14,%%ymm3 	\n\t"\
		"vmulpd		%%ymm0,%%ymm8 ,%%ymm0	\n\t"	/* ainv*(2 - d*ainv) = 1/d accurate to ~53 bits */\
		"vmulpd		%%ymm1,%%ymm9 ,%%ymm1	\n\t"\
		"vmulpd		%%ymm2,%%ymm10,%%ymm2	\n\t"\
		"vmulpd		%%ymm3,%%ymm11,%%ymm3	\n\t"\
		/* Multiply by the sine terms to get quotients in ymm4-7 (still need inverses in ymm0-3 for cosine-ratios): */\
		"vmulpd	0x80(%%rax),%%ymm0,%%ymm4	\n\t"\
		"vmulpd	0xa0(%%rax),%%ymm1,%%ymm5	\n\t"\
		"vmulpd	0xc0(%%rax),%%ymm2,%%ymm6	\n\t"\
		"vmulpd	0xe0(%%rax),%%ymm3,%%ymm7	\n\t"\
		/* Outputs overwrite sines: */\
		"vmovaps	%%ymm4,0x80(%%rax)	\n\t"\
		"vmovaps	%%ymm5,0xa0(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0xc0(%%rax)	\n\t"\
		"vmovaps	%%ymm7,0xe0(%%rax)	\n\t"\
		/* Set up for cosine ratios: */\
		"vpermilpd	$11,%%ymm0,%%ymm2 	\n\t"/* permute ymm0 = 1/c[3123] to get ymm2 = 1/c[1123], then *= [c3,c5,c6,0] */\
		"vmulpd	0x100(%%rax),%%ymm2,%%ymm2	\n\t"/* ymm2 = [c3,c5,c6, 0] * 1/c[1123] = [c31,c51,c62,  0] */\
		"vmulpd	0x120(%%rax),%%ymm0,%%ymm0	\n\t"/* ymm0 = [c7,c9,cA,cB] * 1/c[3123] = [c73,c91,cA2,cB3] */\
		"vmulpd	0x140(%%rax),%%ymm1,%%ymm1	\n\t"/* ymm1 = [cC,cD,cE,cF] * 1/c[4567] = [cC4,cD5,cE6,cF7] */\
		/* Outputs into slots imm. above above tangents: */\
		"vmovaps	%%ymm2,0x100(%%rax)	\n\t"\
		"vmovaps	%%ymm0,0x120(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x140(%%rax)	\n\t"\
		/* Scalar data starting at add0 = cc0 laid out as below. Ensuing C code assumed to massage into a scalar-data analog of the 4-copy layout. */\
		/* add0 + 0x[  0,  8, 10, 18]: [c3 ,c1 ,c2 ,c3 ] Cosines: */\
		/* add0 + 0x[ 20, 28, 30, 38]: [c4 ,c5 ,c6 ,c7 ] */\
		/* add0 + 0x[ 40, 48, 50, 58]: [c8 ,c9 ,cA ,cB ] */\
		/* add0 + 0x[ 60, 68, 70, 78]: [cC ,cD ,cE ,cF ] */\
		/* add0 + 0x[ 80, 88, 90, 98]: [-- ,r1 ,r2 ,r3 ] Tangents: */\
		/* add0 + 0x[ a0, a8, b0, b8]: [r4 ,r5 ,r6 ,r7 ] */\
		/* add0 + 0x[ c0, c8, d0, d8]: [r8 ,r9 ,rA ,rB ] */\
		/* add0 + 0x[ e0, e8, f0, f8]: [rC ,rD ,rE ,rF ] */\
		/* add0 + 0x[100,108,110,118]: [c31,c51,c62,  0] Cosine ratios: */\
		/* add0 + 0x[120,128,130,138]: [c73,c91,cA2,cB3] */\
		/* add0 + 0x[140,148,150,158]: [cC4,cD5,cE6,cF7] */\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__one] "m" (Xone)\
		: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm14"	/* Clobbered registers */\
	);\
	}

	// Here are the hexadecimal byte offsets w.r.to the __c = cc0 base-root address of the various derived sincos terms:
	// Datum	Offset
	// ------	------
	// __c1_c	0x000
	// __sc		0x040
	// __c1i2	0x080
	// __c2i2	0x0c0
	// __c8		0x100	__r8		0x140
	// __c4		0x180	__r4		0x1C0
	// __cC4	0x200	__rC		0x240
	// __c2		0x280	__r2		0x2C0
	// __cA2	0x300	__rA		0x340
	// __c62	0x380	__r6		0x3C0
	// __cE6	0x400	__rE		0x440
	// __c1		0x480	__r1		0x4C0
	// __c91	0x500	__r9		0x540
	// __c51	0x580	__r5		0x5C0
	// __cD5	0x600	__rD		0x640
	// __c31	0x680	__r3		0x6C0
	// __cB3	0x700	__rB		0x740
	// __c73	0x780	__r7		0x7C0
	// __cF7	0x800	__rF		0x840

	// [NB: column-folded 32-register version which processes pairs of 4-DFT blocks side-by-side of this macro was a bust,
	// but the address-comp streamlining I did in prep for it cut cycles by 7%! Did similar for the DIT version, as well.]

	/**** Cf. AVX2 version of this macro for the inline comments ****/
	#define SSE2_RADIX16_DIF_TWIDDLE_1(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr1,Xisrt2,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
		/*...Block 1:	*/\
		"movq	%[__add0],%%rax				\n\t"\
	"movslq	%[__pfetch_addr1],%%r13	\n\t"	/* Prefetch base-index offset */\
	"leaq	(%%rax,%%r13,8),%%r13	\n\t"	/* base-addr + fetch-ahead idx */\
		"movslq	%[__p4],%%rbx				\n\t"\
		"movslq	%[__p8],%%rcx				\n\t"\
		"movslq	%[__p12],%%rdx				\n\t"\
		"movq	%[__isrt2],%%rsi 			\n\t"\
	"xorq	%%r14,%%r14	\n\t"/* Zero r14 and include in LEA to eliminate "scale factor of 8 without an index register" assembler warnings */\
	"leaq	(%%r14,%%rbx,8),%%r14	\n\t"/* Save copy of p4 pointer offset for prefetch */\
	"prefetcht1	%c[__pfetch_dist](%%r13)\n\t"/* [base-addr + fetch-ahead] + p0 */\
		"leaq	(%%rax,%%rbx,8),%%rbx		\n\t"\
		"leaq	(%%rax,%%rcx,8),%%rcx		\n\t"\
		"leaq	(%%rax,%%rdx,8),%%rdx		\n\t"\
		"movq	%[__r1],%%rdi	\n\t"/* ptr to local-mem */\
		"addq	$0x040,%%rsi 	\n\t"/* cc0 */\
		"vmovaps		     (%%rcx),%%zmm4 		\n\t	vmovaps			0x040(%%rcx),%%zmm5 		\n\t"\
		"vmovaps		     (%%rax),%%zmm0 		\n\t	vmovaps			0x040(%%rax),%%zmm1 		\n\t"\
		"vmovaps		%%zmm4,%%zmm6				\n\t	vmovaps			0x140(%%rsi),%%zmm13		\n\t"\
		"vmovaps		0x1c0(%%rsi),%%zmm14		\n\t	vmovaps			0x240(%%rsi),%%zmm15		\n\t"\
		"vfnmadd231pd	%%zmm5 ,%%zmm13,%%zmm4 		\n\t	 vfmadd231pd	%%zmm6 ,%%zmm13,%%zmm5 		\n\t"\
		"vmovaps		     (%%rbx),%%zmm8			\n\t	vmovaps			0x040(%%rbx),%%zmm9 		\n\t"\
		"vfnmadd231pd	0x040(%%rbx),%%zmm14,%%zmm8 \n\t	 vfmadd231pd	     (%%rbx),%%zmm14,%%zmm9 \n\t"\
		"vmovaps		     (%%rdx),%%zmm6			\n\t	vmovaps			0x040(%%rdx),%%zmm7 		\n\t"\
		"vfnmadd231pd	0x040(%%rdx),%%zmm15,%%zmm6 \n\t	 vfmadd231pd	     (%%rdx),%%zmm15,%%zmm7 \n\t"\
		"vmovaps		0x200(%%rsi),%%zmm13		\n\t	vmovaps			0x100(%%rsi),%%zmm14		\n\t"\
		"vmovaps		%%zmm8 ,%%zmm10				\n\t	vmovaps			%%zmm0,%%zmm2 				\n\t"\
		" vfmadd231pd	%%zmm6 ,%%zmm13,%%zmm8 		\n\t	 vfmadd231pd	%%zmm4 ,%%zmm14,%%zmm0 		\n\t"\
		"vmovaps		%%zmm9 ,%%zmm11				\n\t	vmovaps			%%zmm1 ,%%zmm3 				\n\t"\
		" vfmadd231pd	%%zmm7 ,%%zmm13,%%zmm9 		\n\t	 vfmadd231pd	%%zmm5 ,%%zmm14,%%zmm1 		\n\t"\
		"vfnmadd231pd	%%zmm6 ,%%zmm13,%%zmm10		\n\t	vfnmadd231pd	%%zmm4 ,%%zmm14,%%zmm2 		\n\t"\
		"vfnmadd231pd	%%zmm7 ,%%zmm13,%%zmm11		\n\t	vfnmadd231pd	%%zmm5 ,%%zmm14,%%zmm3 		\n\t"\
		"vmovaps		0x180(%%rsi),%%zmm15		\n\t												\n\t"\
		"vmovaps		%%zmm0 ,%%zmm4 				\n\t	vmovaps			%%zmm1 ,%%zmm5 				\n\t"\
		"vfnmadd231pd	%%zmm8 ,%%zmm15,%%zmm4 		\n\t	 vfmadd231pd	%%zmm8 ,%%zmm15,%%zmm0 		\n\t"\
		"vfnmadd231pd	%%zmm9 ,%%zmm15,%%zmm5 		\n\t	 vfmadd231pd	%%zmm9 ,%%zmm15,%%zmm1 		\n\t"\
		"vmovaps		%%zmm2 ,%%zmm6 				\n\t	vmovaps			%%zmm3 ,%%zmm7 				\n\t"\
		" vfmadd231pd	%%zmm11,%%zmm15,%%zmm6 		\n\t	vfnmadd231pd	%%zmm11,%%zmm15,%%zmm2 		\n\t"\
		"vfnmadd231pd	%%zmm10,%%zmm15,%%zmm7 		\n\t	 vfmadd231pd	%%zmm10,%%zmm15,%%zmm3 		\n\t"\
		"vmovaps		%%zmm4 ,0x100(%%rdi)		\n\t	vmovaps			%%zmm0 ,     (%%rdi)		\n\t"\
		"vmovaps		%%zmm5 ,0x140(%%rdi)		\n\t	vmovaps			%%zmm1 ,0x040(%%rdi)		\n\t"\
		"vmovaps		%%zmm6 ,0x180(%%rdi)		\n\t	vmovaps			%%zmm2 ,0x080(%%rdi)		\n\t"\
		"vmovaps		%%zmm7 ,0x1c0(%%rdi)		\n\t	vmovaps			%%zmm3 ,0x0c0(%%rdi)		\n\t"\
		"\n\t"\
		/*...Block 2: Register indices for the 8 t-data = [t-index - 8]: */\
		"movslq	%[__p2],%%r9		\n\t"\
		"vmovaps		(%%rax,%%r9,8),%%zmm0		\n\t	vmovaps		0x040(%%rax,%%r9,8),%%zmm1		\n\t"\
		"vmovaps		(%%rcx,%%r9,8),%%zmm4		\n\t	vmovaps		0x040(%%rcx,%%r9,8),%%zmm5		\n\t"\
		"vmovaps		(%%rbx,%%r9,8),%%zmm8		\n\t	vmovaps		0x040(%%rbx,%%r9,8),%%zmm9		\n\t"\
		"vmovaps		(%%rdx,%%r9,8),%%zmm6		\n\t	vmovaps		0x040(%%rdx,%%r9,8),%%zmm7		\n\t"\
		"vmovaps		0x2c0(%%rsi),%%zmm12		\n\t	vmovaps		%%zmm0,%%zmm2					\n\t"\
		"vmovaps		0x340(%%rsi),%%zmm13		\n\t	vmovaps		%%zmm4,%%zmm3					\n\t"\
		"vmovaps		0x3c0(%%rsi),%%zmm14		\n\t	vmovaps		%%zmm8,%%zmm10					\n\t"\
		"vmovaps		0x440(%%rsi),%%zmm15		\n\t	vmovaps		%%zmm6,%%zmm11					\n\t"\
		"vfnmadd231pd	%%zmm1,%%zmm12,%%zmm0		\n\t	vfmadd231pd	%%zmm2 ,%%zmm12,%%zmm1			\n\t"\
		"vfnmadd231pd	%%zmm5,%%zmm13,%%zmm4		\n\t	vfmadd231pd	%%zmm3 ,%%zmm13,%%zmm5			\n\t"\
		"vfnmadd231pd	%%zmm9,%%zmm14,%%zmm8		\n\t	vfmadd231pd	%%zmm10,%%zmm14,%%zmm9			\n\t"\
		"vfnmadd231pd	%%zmm7,%%zmm15,%%zmm6		\n\t	vfmadd231pd	%%zmm11,%%zmm15,%%zmm7			\n\t"\
		"vmovaps		0x400(%%rsi),%%zmm13		\n\t	vmovaps			0x300(%%rsi),%%zmm14		\n\t"\
		"vmovaps		%%zmm8 ,%%zmm10				\n\t	vmovaps			%%zmm0,%%zmm2 				\n\t"\
		" vfmadd231pd	%%zmm6 ,%%zmm13,%%zmm8 		\n\t	 vfmadd231pd	%%zmm4 ,%%zmm14,%%zmm0 		\n\t"\
		"vmovaps		%%zmm9 ,%%zmm11				\n\t	vmovaps			%%zmm1 ,%%zmm3 				\n\t"\
		" vfmadd231pd	%%zmm7 ,%%zmm13,%%zmm9 		\n\t	 vfmadd231pd	%%zmm5 ,%%zmm14,%%zmm1 		\n\t"\
		"vfnmadd231pd	%%zmm6 ,%%zmm13,%%zmm10		\n\t	vfnmadd231pd	%%zmm4 ,%%zmm14,%%zmm2 		\n\t"\
		"vfnmadd231pd	%%zmm7 ,%%zmm13,%%zmm11		\n\t	vfnmadd231pd	%%zmm5 ,%%zmm14,%%zmm3 		\n\t"\
		"vmovaps		0x380(%%rsi),%%zmm15		\n\t												\n\t"\
		"vmovaps		%%zmm0 ,%%zmm4 				\n\t	vmovaps			%%zmm1 ,%%zmm5 				\n\t"\
		"vfnmadd231pd	%%zmm8 ,%%zmm15,%%zmm4 		\n\t	 vfmadd231pd	%%zmm8 ,%%zmm15,%%zmm0 		\n\t"\
		"vfnmadd231pd	%%zmm9 ,%%zmm15,%%zmm5 		\n\t	 vfmadd231pd	%%zmm9 ,%%zmm15,%%zmm1 		\n\t"\
		"vmovaps		%%zmm2 ,%%zmm6 				\n\t	vmovaps			%%zmm3 ,%%zmm7 				\n\t"\
		" vfmadd231pd	%%zmm11,%%zmm15,%%zmm6 		\n\t	vfnmadd231pd	%%zmm11,%%zmm15,%%zmm2 		\n\t"\
		"vfnmadd231pd	%%zmm10,%%zmm15,%%zmm7 		\n\t	 vfmadd231pd	%%zmm10,%%zmm15,%%zmm3 		\n\t"\
		"vmovaps		%%zmm4 ,0x300(%%rdi)		\n\t	vmovaps			%%zmm0 ,0x200(%%rdi)		\n\t"\
		"vmovaps		%%zmm5 ,0x340(%%rdi)		\n\t	vmovaps			%%zmm1 ,0x240(%%rdi)		\n\t"\
		"vmovaps		%%zmm6 ,0x380(%%rdi)		\n\t	vmovaps			%%zmm2 ,0x280(%%rdi)		\n\t"\
		"vmovaps		%%zmm7 ,0x3c0(%%rdi)		\n\t	vmovaps			%%zmm3 ,0x2c0(%%rdi)		\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%r14)\n\t"/* ...+p4 */\
		"\n\t"\
		/*...Block 3: Register indices for the 8 t-date = [t-index - 16]: */\
		"movslq	%[__p1],%%r8		\n\t"\
		"vmovaps		(%%rax,%%r8,8),%%zmm0		\n\t	vmovaps		0x040(%%rax,%%r8,8),%%zmm1		\n\t"\
		"vmovaps		(%%rcx,%%r8,8),%%zmm4		\n\t	vmovaps		0x040(%%rcx,%%r8,8),%%zmm5		\n\t"\
		"vmovaps		(%%rbx,%%r8,8),%%zmm8		\n\t	vmovaps		0x040(%%rbx,%%r8,8),%%zmm9		\n\t"\
		"vmovaps		(%%rdx,%%r8,8),%%zmm6		\n\t	vmovaps		0x040(%%rdx,%%r8,8),%%zmm7		\n\t"\
		"vmovaps		0x4c0(%%rsi),%%zmm12		\n\t	vmovaps		%%zmm0,%%zmm2					\n\t"\
		"vmovaps		0x540(%%rsi),%%zmm13		\n\t	vmovaps		%%zmm4,%%zmm3					\n\t"\
		"vmovaps		0x5c0(%%rsi),%%zmm14		\n\t	vmovaps		%%zmm8,%%zmm10					\n\t"\
		"vmovaps		0x640(%%rsi),%%zmm15		\n\t	vmovaps		%%zmm6,%%zmm11					\n\t"\
		"vfnmadd231pd	%%zmm1,%%zmm12,%%zmm0		\n\t	vfmadd231pd	%%zmm2 ,%%zmm12,%%zmm1			\n\t"\
		"vfnmadd231pd	%%zmm5,%%zmm13,%%zmm4		\n\t	vfmadd231pd	%%zmm3 ,%%zmm13,%%zmm5			\n\t"\
		"vfnmadd231pd	%%zmm9,%%zmm14,%%zmm8		\n\t	vfmadd231pd	%%zmm10,%%zmm14,%%zmm9			\n\t"\
		"vfnmadd231pd	%%zmm7,%%zmm15,%%zmm6		\n\t	vfmadd231pd	%%zmm11,%%zmm15,%%zmm7			\n\t"\
		"vmovaps		0x600(%%rsi),%%zmm13		\n\t	vmovaps			0x500(%%rsi),%%zmm14		\n\t"\
		"vmovaps		%%zmm8 ,%%zmm10				\n\t	vmovaps			%%zmm0,%%zmm2 				\n\t"\
		" vfmadd231pd	%%zmm6 ,%%zmm13,%%zmm8 		\n\t	 vfmadd231pd	%%zmm4 ,%%zmm14,%%zmm0 		\n\t"\
		"vmovaps		%%zmm9 ,%%zmm11				\n\t	vmovaps			%%zmm1 ,%%zmm3 				\n\t"\
		" vfmadd231pd	%%zmm7 ,%%zmm13,%%zmm9 		\n\t	 vfmadd231pd	%%zmm5 ,%%zmm14,%%zmm1 		\n\t"\
		"vfnmadd231pd	%%zmm6 ,%%zmm13,%%zmm10		\n\t	vfnmadd231pd	%%zmm4 ,%%zmm14,%%zmm2 		\n\t"\
		"vfnmadd231pd	%%zmm7 ,%%zmm13,%%zmm11		\n\t	vfnmadd231pd	%%zmm5 ,%%zmm14,%%zmm3 		\n\t"\
		"vmovaps		0x580(%%rsi),%%zmm15		\n\t												\n\t"\
		"vmovaps		%%zmm0 ,%%zmm4 				\n\t	vmovaps			%%zmm1 ,%%zmm5 				\n\t"\
		"vfnmadd231pd	%%zmm8 ,%%zmm15,%%zmm4 		\n\t	 vfmadd231pd	%%zmm8 ,%%zmm15,%%zmm0 		\n\t"\
		"vfnmadd231pd	%%zmm9 ,%%zmm15,%%zmm5 		\n\t	 vfmadd231pd	%%zmm9 ,%%zmm15,%%zmm1 		\n\t"\
		"vmovaps		%%zmm2 ,%%zmm6 				\n\t	vmovaps			%%zmm3 ,%%zmm7 				\n\t"\
		" vfmadd231pd	%%zmm11,%%zmm15,%%zmm6 		\n\t	vfnmadd231pd	%%zmm11,%%zmm15,%%zmm2 		\n\t"\
		"vfnmadd231pd	%%zmm10,%%zmm15,%%zmm7 		\n\t	 vfmadd231pd	%%zmm10,%%zmm15,%%zmm3 		\n\t"\
		"vmovaps		%%zmm4 ,0x500(%%rdi)		\n\t	vmovaps			%%zmm0 ,0x400(%%rdi)		\n\t"\
		"vmovaps		%%zmm5 ,0x540(%%rdi)		\n\t	vmovaps			%%zmm1 ,0x440(%%rdi)		\n\t"\
		"vmovaps		%%zmm6 ,0x580(%%rdi)		\n\t	vmovaps			%%zmm2 ,0x480(%%rdi)		\n\t"\
		"vmovaps		%%zmm7 ,0x5c0(%%rdi)		\n\t	vmovaps			%%zmm3 ,0x4c0(%%rdi)		\n\t"\
		"\n\t"\
		/*...Block 4: Register indices for the 8 t-date = [t-index - 24]; __p2 << 3 still in %%r9: */\
		"addq	%%r8,%%r9			\n\t"/* p3 */\
		"vmovaps		(%%rax,%%r9,8),%%zmm0		\n\t	vmovaps		0x040(%%rax,%%r9,8),%%zmm1		\n\t"\
		"vmovaps		(%%rcx,%%r9,8),%%zmm4		\n\t	vmovaps		0x040(%%rcx,%%r9,8),%%zmm5		\n\t"\
		"vmovaps		(%%rbx,%%r9,8),%%zmm8		\n\t	vmovaps		0x040(%%rbx,%%r9,8),%%zmm9		\n\t"\
		"vmovaps		(%%rdx,%%r9,8),%%zmm6		\n\t	vmovaps		0x040(%%rdx,%%r9,8),%%zmm7		\n\t"\
		"vmovaps		0x6c0(%%rsi),%%zmm12		\n\t	vmovaps		%%zmm0,%%zmm2					\n\t"\
		"vmovaps		0x740(%%rsi),%%zmm13		\n\t	vmovaps		%%zmm4,%%zmm3					\n\t"\
		"vmovaps		0x7c0(%%rsi),%%zmm14		\n\t	vmovaps		%%zmm8,%%zmm10					\n\t"\
		"vmovaps		0x840(%%rsi),%%zmm15		\n\t	vmovaps		%%zmm6,%%zmm11					\n\t"\
		"vfnmadd231pd	%%zmm1,%%zmm12,%%zmm0		\n\t	vfmadd231pd	%%zmm2 ,%%zmm12,%%zmm1			\n\t"\
		"vfnmadd231pd	%%zmm5,%%zmm13,%%zmm4		\n\t	vfmadd231pd	%%zmm3 ,%%zmm13,%%zmm5			\n\t"\
		"vfnmadd231pd	%%zmm9,%%zmm14,%%zmm8		\n\t	vfmadd231pd	%%zmm10,%%zmm14,%%zmm9			\n\t"\
		"vfnmadd231pd	%%zmm7,%%zmm15,%%zmm6		\n\t	vfmadd231pd	%%zmm11,%%zmm15,%%zmm7			\n\t"\
		"vmovaps		0x800(%%rsi),%%zmm13		\n\t	vmovaps			0x700(%%rsi),%%zmm14		\n\t"\
		"vmovaps		%%zmm8 ,%%zmm10				\n\t	vmovaps			%%zmm0,%%zmm2 				\n\t"\
		" vfmadd231pd	%%zmm6 ,%%zmm13,%%zmm8 		\n\t	 vfmadd231pd	%%zmm4 ,%%zmm14,%%zmm0 		\n\t"\
		"vmovaps		%%zmm9 ,%%zmm11				\n\t	vmovaps			%%zmm1 ,%%zmm3 				\n\t"\
		" vfmadd231pd	%%zmm7 ,%%zmm13,%%zmm9 		\n\t	 vfmadd231pd	%%zmm5 ,%%zmm14,%%zmm1 		\n\t"\
		"vfnmadd231pd	%%zmm6 ,%%zmm13,%%zmm10		\n\t	vfnmadd231pd	%%zmm4 ,%%zmm14,%%zmm2 		\n\t"\
		"vfnmadd231pd	%%zmm7 ,%%zmm13,%%zmm11		\n\t	vfnmadd231pd	%%zmm5 ,%%zmm14,%%zmm3 		\n\t"\
		"vmovaps		0x780(%%rsi),%%zmm15		\n\t												\n\t"\
		"vmovaps		%%zmm0 ,%%zmm4 				\n\t	vmovaps			%%zmm1 ,%%zmm5 				\n\t"\
		"vfnmadd231pd	%%zmm8 ,%%zmm15,%%zmm4 		\n\t	 vfmadd231pd	%%zmm8 ,%%zmm15,%%zmm0 		\n\t"\
		"vfnmadd231pd	%%zmm9 ,%%zmm15,%%zmm5 		\n\t	 vfmadd231pd	%%zmm9 ,%%zmm15,%%zmm1 		\n\t"\
		"vmovaps		%%zmm2 ,%%zmm6 				\n\t	vmovaps			%%zmm3 ,%%zmm7 				\n\t"\
		" vfmadd231pd	%%zmm11,%%zmm15,%%zmm6 		\n\t	vfnmadd231pd	%%zmm11,%%zmm15,%%zmm2 		\n\t"\
		"vfnmadd231pd	%%zmm10,%%zmm15,%%zmm7 		\n\t	 vfmadd231pd	%%zmm10,%%zmm15,%%zmm3 		\n\t"\
		"vmovaps		%%zmm4 ,0x700(%%rdi)		\n\t	vmovaps			%%zmm0 ,0x600(%%rdi)		\n\t"\
		"vmovaps		%%zmm5 ,0x740(%%rdi)		\n\t	vmovaps			%%zmm1 ,0x640(%%rdi)		\n\t"\
		"vmovaps		%%zmm6 ,0x780(%%rdi)		\n\t	vmovaps			%%zmm2 ,0x680(%%rdi)		\n\t"\
		"vmovaps		%%zmm7 ,0x7c0(%%rdi)		\n\t	vmovaps			%%zmm3 ,0x6c0(%%rdi)		\n\t"\
		"\n\t"\
	/*************************************************************************************/\
	/*  And now do four more radix-4 transforms, including the internal twiddle factors: */\
	/*************************************************************************************/\
		/* Block 1: r[abcd]x = add0+p[0,4,8,12], r8 = p1, r9 = p3: */\
		"leaq	(%%rax,%%r9,8),%%rdx	\n\t"/* add0+p3 */\
		"subq	%%r8,%%r9				\n\t"/* p2 */\
		"leaq	(%%rax,%%r8,8),%%rbx	\n\t"/* add0+p1 */\
		"leaq	(%%rax,%%r9,8),%%rcx	\n\t"/* add0+p2 */\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%r14,2)\n\t"/* ...+p8 */\
	"leaq	(%%r14,%%r14,2),%%r14	\n\t"	/* p4 + (p4*2) = p12, ptr-offset form */\
		"\n\t"\
		/* Read t0,8,16,24 from local store ... Do 4 Im-part FMAs first, because their outputs needed 1st */\
		"vmovaps		0x280(%%rsi),%%zmm14		\n\t	vmovaps			0x680(%%rsi),%%zmm15		\n\t"\
		"vmovaps		     (%%rdi),%%zmm0 		\n\t	vmovaps			0x200(%%rdi),%%zmm2 		\n\t"\
		"vmovaps		0x400(%%rdi),%%zmm4 		\n\t	vmovaps			0x600(%%rdi),%%zmm6 		\n\t"\
		"vmovaps			 %%zmm0 ,%%zmm8 		\n\t	vmovaps				 %%zmm4 ,%%zmm10		\n\t"\
		" vfmadd231pd	%%zmm2 ,%%zmm14,%%zmm0 		\n\t	 vfmadd231pd	%%zmm6 ,%%zmm15,%%zmm4 		\n\t"\
		"vmovaps		0x040(%%rdi),%%zmm1 		\n\t	vmovaps			0x240(%%rdi),%%zmm3 		\n\t"\
		"vmovaps		0x440(%%rdi),%%zmm5 		\n\t	vmovaps			0x640(%%rdi),%%zmm7 		\n\t"\
		"vmovaps			 %%zmm1 ,%%zmm9 		\n\t	vmovaps				 %%zmm5 ,%%zmm11		\n\t"\
		" vfmadd231pd	%%zmm3 ,%%zmm14,%%zmm1 		\n\t	 vfmadd231pd	%%zmm7 ,%%zmm15,%%zmm5 		\n\t"\
		"vfnmadd231pd	%%zmm2 ,%%zmm14,%%zmm8 		\n\t	vfnmadd231pd	%%zmm6 ,%%zmm15,%%zmm10		\n\t"\
		"vfnmadd231pd	%%zmm3 ,%%zmm14,%%zmm9 		\n\t	vfnmadd231pd	%%zmm7 ,%%zmm15,%%zmm11		\n\t"\
		"vmovaps		0x480(%%rsi),%%zmm15		\n\t												\n\t"\
		"vmovaps			 %%zmm0 ,%%zmm12		\n\t	vmovaps				 %%zmm1 ,%%zmm13		\n\t"\
		" vfmadd231pd	%%zmm4 ,%%zmm15,%%zmm0 		\n\t	 vfmadd231pd	%%zmm5 ,%%zmm15,%%zmm1 		\n\t"\
		"vfnmadd231pd	%%zmm4 ,%%zmm15,%%zmm12		\n\t	vfnmadd231pd	%%zmm5 ,%%zmm15,%%zmm13		\n\t"\
		"vmovaps			 %%zmm8 ,%%zmm2 		\n\t	vmovaps				 %%zmm9 ,%%zmm3 		\n\t"\
		"vfnmadd231pd	%%zmm11,%%zmm15,%%zmm2 		\n\t	 vfmadd231pd	%%zmm10,%%zmm15,%%zmm3 		\n\t"\
		" vfmadd231pd	%%zmm11,%%zmm15,%%zmm8 		\n\t	vfnmadd231pd	%%zmm10,%%zmm15,%%zmm9 		\n\t"\
		/* Write outputs back to main array: */\
		"vmovaps		%%zmm0 ,     (%%rax)		\n\t	vmovaps			%%zmm1 ,0x040(%%rax)		\n\t"\
		"vmovaps		%%zmm12,     (%%rbx)		\n\t	vmovaps			%%zmm13,0x040(%%rbx)		\n\t"\
		"vmovaps		%%zmm2 ,     (%%rcx)		\n\t	vmovaps			%%zmm3 ,0x040(%%rcx)		\n\t"\
		"vmovaps		%%zmm8 ,     (%%rdx)		\n\t	vmovaps			%%zmm9 ,0x040(%%rdx)		\n\t"\
		"\n\t"\
		/*...Block 3: t4,12,20,28 */\
		"movslq	%[__p4],%%r8		\n\t"\
		"vmovaps	0x880(%%rsi),%%zmm15	\n\t"/* cc0 + 0x44 = __two; Actually holds 1.0 in AVX2 mode */\
		"vmovaps		0x500(%%rdi),%%zmm4 		\n\t	vmovaps			0x540(%%rdi),%%zmm5 		\n\t"\
		"vmovaps		0x700(%%rdi),%%zmm6 		\n\t	vmovaps			0x740(%%rdi),%%zmm7 		\n\t"\
		"vmovaps			 %%zmm4 ,%%zmm10 		\n\t	vmovaps				 %%zmm7 ,%%zmm11		\n\t"\
		"vfnmadd231pd	%%zmm5 ,%%zmm15,%%zmm10		\n\t	 vfmadd231pd	%%zmm5 ,%%zmm15,%%zmm4 		\n\t"\
		" vfmadd231pd	%%zmm6 ,%%zmm15,%%zmm11		\n\t	vfnmadd231pd	%%zmm6 ,%%zmm15,%%zmm7 		\n\t"\
		"vmovaps		0x280(%%rsi),%%zmm14		\n\t	vmovaps			0x680(%%rsi),%%zmm15		\n\t"\
		"vmovaps		0x100(%%rdi),%%zmm0 		\n\t	vmovaps			0x140(%%rdi),%%zmm1 		\n\t"\
		"vmovaps		0x300(%%rdi),%%zmm2 		\n\t	vmovaps			0x340(%%rdi),%%zmm3 		\n\t"\
		"vmovaps			 %%zmm0 ,%%zmm8 		\n\t	vmovaps				 %%zmm1 ,%%zmm9 		\n\t"\
		"vmovaps			 %%zmm10,%%zmm5 		\n\t	vmovaps				 %%zmm4 ,%%zmm12		\n\t"\
		" vfmadd231pd	%%zmm3 ,%%zmm14,%%zmm8 		\n\t	 vfmadd231pd	%%zmm11,%%zmm15,%%zmm5 		\n\t"\
		"vfnmadd231pd	%%zmm2 ,%%zmm14,%%zmm9 		\n\t	 vfmadd231pd	%%zmm7 ,%%zmm15,%%zmm4 		\n\t"\
		"vfnmadd231pd	%%zmm3 ,%%zmm14,%%zmm0 		\n\t	vfnmadd231pd	%%zmm11,%%zmm15,%%zmm10		\n\t"\
		" vfmadd231pd	%%zmm2 ,%%zmm14,%%zmm1 		\n\t	vfnmadd231pd	%%zmm7 ,%%zmm15,%%zmm12		\n\t"\
		"vmovaps		0x080(%%rsi),%%zmm15		\n\t												\n\t"\
		"vmovaps			 %%zmm8 ,%%zmm2 		\n\t	vmovaps				 %%zmm9 ,%%zmm3 		\n\t"\
		"vfnmadd231pd	%%zmm4 ,%%zmm15,%%zmm2 		\n\t	 vfmadd231pd	%%zmm5 ,%%zmm15,%%zmm3 		\n\t"\
		" vfmadd231pd	%%zmm4 ,%%zmm15,%%zmm8 		\n\t	vfnmadd231pd	%%zmm5 ,%%zmm15,%%zmm9 		\n\t"\
		"vmovaps			 %%zmm0 ,%%zmm11		\n\t	vmovaps				 %%zmm1 ,%%zmm7 		\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm15,%%zmm0 		\n\t	 vfmadd231pd	%%zmm12,%%zmm15,%%zmm1 		\n\t"\
		"vfnmadd231pd	%%zmm10,%%zmm15,%%zmm11		\n\t	vfnmadd231pd	%%zmm12,%%zmm15,%%zmm7 		\n\t"\
		/* Write outputs back to main array: */\
		"vmovaps		%%zmm2 ,     (%%rcx,%%r8,8)	\n\t	vmovaps			%%zmm3 ,0x040(%%rcx,%%r8,8)	\n\t"\
		"vmovaps		%%zmm8 ,     (%%rdx,%%r8,8)	\n\t	vmovaps			%%zmm9 ,0x040(%%rdx,%%r8,8)	\n\t"\
		"vmovaps		%%zmm0 ,     (%%rax,%%r8,8)	\n\t	vmovaps			%%zmm1 ,0x040(%%rax,%%r8,8)	\n\t"\
		"vmovaps		%%zmm11,     (%%rbx,%%r8,8)	\n\t	vmovaps			%%zmm7 ,0x040(%%rbx,%%r8,8)	\n\t"\
		"\n\t"\
		/*...Block 2: t2,10,18,26 */\
		"leaq	(%%r8,%%r8),%%r9	\n\t"/* p8 */\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%r14)\n\t"/* ...+p12 */\
		"vmovaps		0x280(%%rdi),%%zmm2 		\n\t	vmovaps			0x2c0(%%rdi),%%zmm3 		\n\t"\
		"vmovaps		0x480(%%rdi),%%zmm4 		\n\t	vmovaps			0x4c0(%%rdi),%%zmm5 		\n\t"\
		"vmovaps		0x680(%%rdi),%%zmm6 		\n\t	vmovaps			0x6c0(%%rdi),%%zmm7 		\n\t"\
		"vsubpd		%%zmm3 ,%%zmm2 ,%%zmm12 		\n\t"/* _e = t10-t11; */\
		"vmovaps		0x040(%%rsi),%%zmm15		\n\t"/* load __sc into reg */\
		"vaddpd		%%zmm2 ,%%zmm3 ,%%zmm13			\n\t"/* _f = t10+t11; */\
		"vmovaps			 %%zmm4 ,%%zmm10 		\n\t	vmovaps				 %%zmm7 ,%%zmm8 		\n\t"\
		"vfnmadd231pd	%%zmm5 ,%%zmm15,%%zmm10		\n\t	 vfmsub231pd	%%zmm6 ,%%zmm15,%%zmm8 		\n\t"\
		"vmovaps			 %%zmm5 ,%%zmm11 		\n\t	vmovaps				 %%zmm6 ,%%zmm9 		\n\t"\
		" vfmadd231pd	%%zmm4 ,%%zmm15,%%zmm11		\n\t	 vfmadd231pd	%%zmm7 ,%%zmm15,%%zmm9 		\n\t"\
		"vmovaps		0x080(%%rdi),%%zmm0 		\n\t	vmovaps			0x0c0(%%rdi),%%zmm1 		\n\t"\
		"vmovaps		0x680(%%rsi),%%zmm14		\n\t	vmovaps			0x0c0(%%rsi),%%zmm15		\n\t"\
		"vmovaps			 %%zmm10,%%zmm4 		\n\t	vmovaps				 %%zmm0 ,%%zmm2 		\n\t"\
		" vfmadd231pd	%%zmm8 ,%%zmm14,%%zmm4 		\n\t	 vfmadd231pd	%%zmm12,%%zmm15,%%zmm0 		\n\t"\
		"vmovaps			 %%zmm11,%%zmm5 		\n\t	vmovaps				 %%zmm1 ,%%zmm3 		\n\t"\
		" vfmadd231pd	%%zmm9 ,%%zmm14,%%zmm5 		\n\t	 vfmadd231pd	%%zmm13,%%zmm15,%%zmm1 		\n\t"\
		"vfnmadd231pd	%%zmm8 ,%%zmm14,%%zmm10		\n\t	vfnmadd231pd	%%zmm12,%%zmm15,%%zmm2 		\n\t"\
		"vfnmadd231pd	%%zmm9 ,%%zmm14,%%zmm11		\n\t	vfnmadd231pd	%%zmm13,%%zmm15,%%zmm3 		\n\t"\
		"vmovaps		(%%rsi),%%zmm15				\n\t												\n\t"\
		"vmovaps			 %%zmm0 ,%%zmm8 		\n\t	vmovaps				 %%zmm1 ,%%zmm9 		\n\t"\
		" vfmadd231pd	%%zmm4 ,%%zmm15,%%zmm0 		\n\t	 vfmadd231pd	%%zmm5 ,%%zmm15,%%zmm1 		\n\t"\
		"vfnmadd231pd	%%zmm4 ,%%zmm15,%%zmm8 		\n\t	vfnmadd231pd	%%zmm5 ,%%zmm15,%%zmm9 		\n\t"\
		"vmovaps			 %%zmm2 ,%%zmm12		\n\t	vmovaps				 %%zmm3 ,%%zmm13		\n\t"\
		"vfnmadd231pd	%%zmm11,%%zmm15,%%zmm2 		\n\t	 vfmadd231pd	%%zmm10,%%zmm15,%%zmm3 		\n\t"\
		" vfmadd231pd	%%zmm11,%%zmm15,%%zmm12		\n\t	vfnmadd231pd	%%zmm10,%%zmm15,%%zmm13		\n\t"\
		/* Write outputs back to main array: */\
		"vmovaps		%%zmm0 ,     (%%rax,%%r9,8)	\n\t	vmovaps			%%zmm1 ,0x040(%%rax,%%r9,8)	\n\t"\
		"vmovaps		%%zmm8 ,     (%%rbx,%%r9,8)	\n\t	vmovaps			%%zmm9 ,0x040(%%rbx,%%r9,8)	\n\t"\
		"vmovaps		%%zmm2 ,     (%%rcx,%%r9,8)	\n\t	vmovaps			%%zmm3 ,0x040(%%rcx,%%r9,8)	\n\t"\
		"vmovaps		%%zmm12,     (%%rdx,%%r9,8)	\n\t	vmovaps			%%zmm13,0x040(%%rdx,%%r9,8)	\n\t"\
		"\n\t"\
		/*...Block 4: t6,14,22,30 */\
		"addq	%%r8,%%r9		\n\t"/* pC */\
		"vmovaps		0x380(%%rdi),%%zmm2 		\n\t	vmovaps			0x3c0(%%rdi),%%zmm3 		\n\t"\
		"vmovaps		0x580(%%rdi),%%zmm4 		\n\t	vmovaps			0x5c0(%%rdi),%%zmm5 		\n\t"\
		"vmovaps		0x780(%%rdi),%%zmm6 		\n\t	vmovaps			0x7c0(%%rdi),%%zmm7 		\n\t"\
		"vaddpd		%%zmm2 ,%%zmm3 ,%%zmm10 		\n\t"/* _c = t14+t5; */\
		"vmovaps		0x040(%%rsi),%%zmm15		\n\t"/* load __sc into reg */\
		"vsubpd		%%zmm2 ,%%zmm3 ,%%zmm11			\n\t"/* _d = t15-t14; */\
		"vmovaps			 %%zmm5 ,%%zmm12 		\n\t	vmovaps				 %%zmm4 ,%%zmm13		\n\t"\
		" vfmsub231pd	%%zmm4 ,%%zmm15,%%zmm12		\n\t	 vfmadd231pd	%%zmm5 ,%%zmm15,%%zmm13		\n\t"\
		"vmovaps			 %%zmm6 ,%%zmm8  		\n\t	vmovaps				 %%zmm7 ,%%zmm9 		\n\t"\
		"vfnmadd231pd	%%zmm7 ,%%zmm15,%%zmm8 		\n\t	 vfmadd231pd	%%zmm6 ,%%zmm15,%%zmm9 		\n\t"\
		"vmovaps		0x180(%%rdi),%%zmm0 		\n\t	vmovaps			0x1c0(%%rdi),%%zmm1 		\n\t"\
		"vmovaps		0x680(%%rsi),%%zmm14		\n\t	vmovaps			0x0c0(%%rsi),%%zmm15		\n\t"\
		"vmovaps			 %%zmm1 ,%%zmm3 		\n\t	vmovaps				 %%zmm0 ,%%zmm2 		\n\t"\
		"vfnmadd231pd	%%zmm11,%%zmm15,%%zmm1 		\n\t	vfnmadd231pd	%%zmm10,%%zmm15,%%zmm0 		\n\t"\
		"vmovaps			 %%zmm12,%%zmm4 		\n\t	vmovaps				 %%zmm13,%%zmm5 		\n\t"\
		"vfnmadd231pd	%%zmm8 ,%%zmm14,%%zmm4 		\n\t	vfnmadd231pd	%%zmm9 ,%%zmm14,%%zmm5 		\n\t"\
		" vfmadd231pd	%%zmm8 ,%%zmm14,%%zmm12		\n\t	 vfmadd231pd	%%zmm9 ,%%zmm14,%%zmm13		\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm15,%%zmm2 		\n\t	 vfmadd231pd	%%zmm11,%%zmm15,%%zmm3 		\n\t"\
		"vmovaps		(%%rsi),%%zmm15				\n\t												\n\t"\
		"vmovaps			 %%zmm0 ,%%zmm8 		\n\t	vmovaps				 %%zmm1 ,%%zmm9 		\n\t"\
		" vfmadd231pd	%%zmm4 ,%%zmm15,%%zmm0 		\n\t	 vfmadd231pd	%%zmm5 ,%%zmm15,%%zmm1 		\n\t"\
		"vfnmadd231pd	%%zmm4 ,%%zmm15,%%zmm8 		\n\t	vfnmadd231pd	%%zmm5 ,%%zmm15,%%zmm9 		\n\t"\
		"vmovaps			 %%zmm2 ,%%zmm10		\n\t	vmovaps				 %%zmm3 ,%%zmm11		\n\t"\
		"vfnmadd231pd	%%zmm13,%%zmm15,%%zmm2 		\n\t	 vfmadd231pd	%%zmm12,%%zmm15,%%zmm3 		\n\t"\
		" vfmadd231pd	%%zmm13,%%zmm15,%%zmm10		\n\t	vfnmadd231pd	%%zmm12,%%zmm15,%%zmm11		\n\t"\
		/* Write outputs back to main array: */\
		"vmovaps		%%zmm0 ,     (%%rax,%%r9,8)	\n\t	vmovaps			%%zmm1 ,0x040(%%rax,%%r9,8)	\n\t"\
		"vmovaps		%%zmm8 ,     (%%rbx,%%r9,8)	\n\t	vmovaps			%%zmm9 ,0x040(%%rbx,%%r9,8)	\n\t"\
		"vmovaps		%%zmm2 ,     (%%rcx,%%r9,8)	\n\t	vmovaps			%%zmm3 ,0x040(%%rcx,%%r9,8)	\n\t"\
		"vmovaps		%%zmm10,     (%%rdx,%%r9,8)	\n\t	vmovaps			%%zmm11,0x040(%%rdx,%%r9,8)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r1] "m" (Xr1)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	// The non-name-suffixed DIF,DIT macros in AVX2 and AVX512 mode are based on the AVX-mode _V2 macros:
	// 96 ADD, 80 FMA, 48 pure-MUL, 316 MOVAPS
	/*** Dec 2021: Delete the prefetcht1 instructions in order to allow use on IMCI512 (first-gen Xeon Phi) ***/
	#define SSE2_RADIX16_DIF_TWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr0,Xtwo,Xcc0,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/*
		k1 = p2*8;	// Set k1 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k2 = 0x80;	// Set k2 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		// Pass1-DFTs all use same twiddle-triplet, so lit-byte address offset between ptrs to 1st,2nd DFT's twiddles = 0
		i0 = (vec_dbl *)add0; i1 = (vec_dbl *)(add0+p4); i2 = (vec_dbl *)(add0+p8); i3 = (vec_dbl *)(add0+p12);
		o0 = r1; o1 = r1+2; o2 = r1+4; o3 = r1+6;
		c_tmp = cc0;	// c8,4,C
		SSE2_RADIX_04_DIF_3TWIDDLE_X1(i0,i1,i2,i3,k1, two,c_tmp, o0,o1,o2,o3,k2)
	*/\
		"movq	%[__cc0],%%rsi 			\n\t	movslq	%[__p2],%%rdi		\n\t"\
		"movq	%[__add0],%%rax			\n\t	movq	%[__r0],%%r10		\n\t"\
		"movslq	%[__p4],%%rbx			\n\t	leaq	(%%rax,%%rbx,8),%%rbx	\n\t"\
		"movslq	%[__p8],%%rcx			\n\t	leaq	(%%rax,%%rcx,8),%%rcx	\n\t"\
		"movslq	%[__p12],%%rdx			\n\t	leaq	(%%rax,%%rdx,8),%%rdx	\n\t"\
	/*"prefetcht1	%c[__pfetch_dist](%%rax)\n\t"// [base-address + data-fetch-ahead index] + p0 */\
		/* Do the p0,p2 combo: */\
		"vmovaps	     (%%rcx),%%zmm4	\n\t	shlq	$3,%%rdi			\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5	\n\t"\
		"vmovaps	     (%%rsi),%%zmm10\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm11\n\t"\
		"vmovaps	%%zmm4,%%zmm6		\n\t	vmovaps	     (%%rcx,%%rdi),%%zmm12	\n\t"\
		"vmovaps	%%zmm5,%%zmm7		\n\t	vmovaps	0x040(%%rcx,%%rdi),%%zmm13	\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t	vmovaps	     (%%rax,%%rdi),%%zmm8 	\n\t	vmovaps	%%zmm12,%%zmm14		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t	vmovaps	0x040(%%rax,%%rdi),%%zmm9 	\n\t	vmovaps	%%zmm13,%%zmm15		\n\t"\
		"vmulpd	%%zmm10,%%zmm4,%%zmm4	\n\t	vmulpd	%%zmm10,%%zmm12,%%zmm12		\n\t	vmovaps	%%zmm0,%%zmm2		\n\t"\
		"vmulpd	%%zmm10,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm10,%%zmm13,%%zmm13		\n\t	vmovaps	%%zmm1,%%zmm3		\n\t"\
	"vfnmadd231pd %%zmm11,%%zmm7,%%zmm4	\n\t vfnmadd231pd %%zmm11,%%zmm15,%%zmm12	\n\t	vmovaps	%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd231pd  %%zmm11,%%zmm6,%%zmm5	\n\t vfmadd231pd  %%zmm11,%%zmm14,%%zmm13	\n\t	vmovaps	%%zmm9 ,%%zmm11		\n\t"\
		"vaddpd	%%zmm4 ,%%zmm0,%%zmm0	\n\t	vaddpd	%%zmm12,%%zmm8 ,%%zmm8 		\n\t"\
		"vaddpd	%%zmm5 ,%%zmm1,%%zmm1	\n\t	vaddpd	%%zmm13,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd	%%zmm4 ,%%zmm2,%%zmm2	\n\t	vsubpd	%%zmm12,%%zmm10,%%zmm10		\n\t"\
		"vsubpd	%%zmm5 ,%%zmm3,%%zmm3	\n\t	vsubpd	%%zmm13,%%zmm11,%%zmm11		\n\t"\
		"vmovaps	%%zmm0,     (%%r10)	\n\t	vmovaps	%%zmm8 ,0x200(%%r10)		\n\t"/* Spill 1: free up zmm0,1 */\
		"vmovaps	%%zmm1,0x040(%%r10)	\n\t	vmovaps	%%zmm9 ,0x240(%%r10)		\n\t"/* Do the p1,3 combo: */\
		"vmovaps	0x100(%%rsi),%%zmm0	\n\t"\
		"vmovaps	0x140(%%rsi),%%zmm1	\n\t"\
		"vmovaps	     (%%rdx),%%zmm6	\n\t	vmovaps	     (%%rdx,%%rdi),%%zmm14	\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm7	\n\t	vmovaps	0x040(%%rdx,%%rdi),%%zmm15	\n\t"\
		"vmovaps	%%zmm6,%%zmm4		\n\t	vmovaps	%%zmm14,%%zmm12				\n\t"\
		"vmovaps	%%zmm7,%%zmm5		\n\t	vmovaps	%%zmm15,%%zmm13				\n\t"\
		"vmulpd	%%zmm0,%%zmm4,%%zmm4	\n\t	vmulpd	%%zmm0 ,%%zmm12,%%zmm12		\n\t"\
		"vmulpd	%%zmm0,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm0 ,%%zmm13,%%zmm13		\n\t"\
	"vfnmadd231pd %%zmm1,%%zmm7,%%zmm4	\n\t vfnmadd231pd %%zmm1,%%zmm15,%%zmm12	\n\t"\
	"vfmadd231pd  %%zmm1,%%zmm6,%%zmm5	\n\t vfmadd231pd  %%zmm1,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	%%zmm5,0x140(%%r10)	\n\t	vmovaps	%%zmm13,0x340(%%r10)		\n\t"/* Spill 2 */\
		"vmovaps	%%zmm4,0x100(%%r10)	\n\t	vmovaps	%%zmm12,0x300(%%r10)		\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm8	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm9	\n\t"\
		"vmovaps	     (%%rbx),%%zmm6	\n\t	vmovaps	     (%%rbx,%%rdi),%%zmm14	\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm7	\n\t	vmovaps	0x040(%%rbx,%%rdi),%%zmm15	\n\t"\
		"vmovaps	%%zmm6,%%zmm4		\n\t	vmovaps	%%zmm14,%%zmm12				\n\t"\
		"vmovaps	%%zmm7,%%zmm5		\n\t	vmovaps	%%zmm15,%%zmm13				\n\t"\
		"vmulpd	%%zmm8,%%zmm4,%%zmm4	\n\t	vmulpd	%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vmulpd	%%zmm8,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm8 ,%%zmm13,%%zmm13		\n\t"\
	"vfnmadd231pd %%zmm9,%%zmm7,%%zmm4	\n\t vfnmadd231pd %%zmm9,%%zmm15,%%zmm12	\n\t"\
	"vfmadd231pd  %%zmm9,%%zmm6,%%zmm5	\n\t vfmadd231pd  %%zmm9,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	0x100(%%r10),%%zmm0	\n\t	vmovaps	0x300(%%r10),%%zmm8 		\n\t"/* Restore 2 */\
		"vmovaps	0x140(%%r10),%%zmm1	\n\t	vmovaps	0x340(%%r10),%%zmm9 		\n\t"\
		"vmovaps	%%zmm5,%%zmm7		\n\t	vmovaps	%%zmm13,%%zmm15				\n\t	 movq	%[__two],%%rsi	\n\t"\
		"vmovaps	%%zmm4,%%zmm6		\n\t	vmovaps	%%zmm12,%%zmm14				\n\t"\
		"vsubpd	%%zmm0,%%zmm4,%%zmm4	\n\t	vsubpd	%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vsubpd	%%zmm1,%%zmm5,%%zmm5	\n\t	vsubpd	%%zmm9 ,%%zmm13,%%zmm13		\n\t"\
		"vaddpd	%%zmm0,%%zmm6,%%zmm6	\n\t	vaddpd	%%zmm8 ,%%zmm14,%%zmm14		\n\t"\
		"vaddpd	%%zmm1,%%zmm7,%%zmm7	\n\t	vaddpd	%%zmm9 ,%%zmm15,%%zmm15		\n\t"\
		/* Finish radix-4 butterfly and store results: */\
		"vmovaps	     (%%r10),%%zmm0	\n\t	vmovaps	0x200(%%r10),%%zmm8 		\n\t"/* Restore 1 */\
		"vmovaps	0x040(%%r10),%%zmm1	\n\t	vmovaps	0x240(%%r10),%%zmm9 		\n\t"\
		"vsubpd	%%zmm6,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm14,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd	%%zmm5,%%zmm2,%%zmm2	\n\t	vsubpd	%%zmm13,%%zmm10,%%zmm10		\n\t"\
		"vsubpd	%%zmm7,%%zmm1,%%zmm1	\n\t	vsubpd	%%zmm15,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vsubpd	%%zmm12,%%zmm11,%%zmm11		\n\t"\
		"vmovaps	%%zmm0,0x100(%%r10)	\n\t	vmovaps	%%zmm8 ,0x300(%%r10)		\n\t"	/* 2.0, shared by both columns: */\
		"vmovaps	%%zmm2,0x080(%%r10)	\n\t	vmovaps	%%zmm10,0x280(%%r10)		\n\t	vmovaps	(%%rsi),%%zmm0	\n\t"\
		"vmovaps	%%zmm1,0x140(%%r10)	\n\t	vmovaps	%%zmm9 ,0x340(%%r10)		\n\t"\
		"vmovaps	%%zmm3,0x1c0(%%r10)	\n\t	vmovaps	%%zmm11,0x3c0(%%r10)		\n\t"\
	"vfmadd213pd 0x100(%%r10),%%zmm0,%%zmm6 \n\t	vfmadd213pd	%%zmm8 ,%%zmm0,%%zmm14	\n\t"\
	"vfmadd213pd      %%zmm2,%%zmm0,%%zmm5 \n\t	vfmadd213pd	%%zmm10,%%zmm0,%%zmm13	\n\t"\
	"vfmadd213pd      %%zmm1,%%zmm0,%%zmm7 \n\t	vfmadd213pd	%%zmm9 ,%%zmm0,%%zmm15	\n\t"\
	"vfmadd213pd      %%zmm3,%%zmm0,%%zmm4 \n\t	vfmadd213pd	%%zmm11,%%zmm0,%%zmm12	\n\t"\
		"vmovaps	%%zmm6,     (%%r10)	\n\t	vmovaps	%%zmm14,0x200(%%r10)		\n\t"/* don't need reload-from-mem of zmm8/0x300(%%r10) as we do in lcol, but 1 cycle faster with it. [!?] */\
		"vmovaps	%%zmm5,0x180(%%r10)	\n\t	vmovaps	%%zmm13,0x380(%%r10)		\n\t"\
		"vmovaps	%%zmm7,0x040(%%r10)	\n\t	vmovaps	%%zmm15,0x240(%%r10)		\n\t"\
		"vmovaps	%%zmm4,0x0c0(%%r10)	\n\t	vmovaps	%%zmm12,0x2c0(%%r10)		\n\t"\
	/*"prefetcht1	%c[__pfetch_dist](%%rcx)\n\t"// ... + p8 */\
	/*
		i0 = (vec_dbl *)(add0+p1); i1 = (vec_dbl *)(add0+p4+p1); i2 = (vec_dbl *)(add0+p8+p1); i3 = (vec_dbl *)(add0+p12+p1);
		o0 += 16; o1 += 16; o2 += 16; o3 += 16;
		SSE2_RADIX_04_DIF_3TWIDDLE_X1(i0,i1,i2,i3,k1, two,c_tmp, o0,o1,o2,o3,k2)
	*/\
		"movslq	%[__p1],%%r9			\n\t	addq	$0x400,%%r10		\n\t"\
		"leaq	(%%rax,%%r9,8),%%rax	\n\t	movq	%[__cc0],%%rsi 		\n\t"/* repoint rsi from two -> cc0 */\
		"leaq	(%%rbx,%%r9,8),%%rbx	\n\t"\
		"leaq	(%%rcx,%%r9,8),%%rcx	\n\t"\
		"leaq	(%%rdx,%%r9,8),%%rdx	\n\t"\
		/* Do the p0,p2 combo: */\
		"vmovaps	     (%%rcx),%%zmm4	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5	\n\t"\
		"vmovaps	     (%%rsi),%%zmm10\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm11\n\t"\
		"vmovaps	%%zmm4,%%zmm6		\n\t	vmovaps	     (%%rcx,%%rdi),%%zmm12	\n\t"\
		"vmovaps	%%zmm5,%%zmm7		\n\t	vmovaps	0x040(%%rcx,%%rdi),%%zmm13	\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t	vmovaps	     (%%rax,%%rdi),%%zmm8 	\n\t	vmovaps	%%zmm12,%%zmm14		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t	vmovaps	0x040(%%rax,%%rdi),%%zmm9 	\n\t	vmovaps	%%zmm13,%%zmm15		\n\t"\
		"vmulpd	%%zmm10,%%zmm4,%%zmm4	\n\t	vmulpd	%%zmm10,%%zmm12,%%zmm12		\n\t	vmovaps	%%zmm0,%%zmm2		\n\t"\
		"vmulpd	%%zmm10,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm10,%%zmm13,%%zmm13		\n\t	vmovaps	%%zmm1,%%zmm3		\n\t"\
	"vfnmadd231pd %%zmm11,%%zmm7,%%zmm4	\n\t vfnmadd231pd %%zmm11,%%zmm15,%%zmm12	\n\t	vmovaps	%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd231pd  %%zmm11,%%zmm6,%%zmm5	\n\t vfmadd231pd  %%zmm11,%%zmm14,%%zmm13	\n\t	vmovaps	%%zmm9 ,%%zmm11		\n\t"\
		"vaddpd	%%zmm4 ,%%zmm0,%%zmm0	\n\t	vaddpd	%%zmm12,%%zmm8 ,%%zmm8 		\n\t"\
		"vaddpd	%%zmm5 ,%%zmm1,%%zmm1	\n\t	vaddpd	%%zmm13,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd	%%zmm4 ,%%zmm2,%%zmm2	\n\t	vsubpd	%%zmm12,%%zmm10,%%zmm10		\n\t"\
		"vsubpd	%%zmm5 ,%%zmm3,%%zmm3	\n\t	vsubpd	%%zmm13,%%zmm11,%%zmm11		\n\t"\
		"vmovaps	%%zmm0,     (%%r10)	\n\t	vmovaps	%%zmm8 ,0x200(%%r10)		\n\t"/* Spill 1: free up zmm0,1 */\
		"vmovaps	%%zmm1,0x040(%%r10)	\n\t	vmovaps	%%zmm9 ,0x240(%%r10)		\n\t"/* Do the p1,3 combo: */\
		"vmovaps	0x100(%%rsi),%%zmm0	\n\t"\
		"vmovaps	0x140(%%rsi),%%zmm1	\n\t"\
		"vmovaps	     (%%rdx),%%zmm6	\n\t	vmovaps	     (%%rdx,%%rdi),%%zmm14	\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm7	\n\t	vmovaps	0x040(%%rdx,%%rdi),%%zmm15	\n\t"\
		"vmovaps	%%zmm6,%%zmm4		\n\t	vmovaps	%%zmm14,%%zmm12				\n\t"\
		"vmovaps	%%zmm7,%%zmm5		\n\t	vmovaps	%%zmm15,%%zmm13				\n\t"\
		"vmulpd	%%zmm0,%%zmm4,%%zmm4	\n\t	vmulpd	%%zmm0 ,%%zmm12,%%zmm12		\n\t"\
		"vmulpd	%%zmm0,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm0 ,%%zmm13,%%zmm13		\n\t"\
	"vfnmadd231pd %%zmm1,%%zmm7,%%zmm4	\n\t vfnmadd231pd %%zmm1,%%zmm15,%%zmm12	\n\t"\
	"vfmadd231pd  %%zmm1,%%zmm6,%%zmm5	\n\t vfmadd231pd  %%zmm1,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	%%zmm5,0x140(%%r10)	\n\t	vmovaps	%%zmm13,0x340(%%r10)		\n\t"/* Spill 2 */\
		"vmovaps	%%zmm4,0x100(%%r10)	\n\t	vmovaps	%%zmm12,0x300(%%r10)		\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm8	\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm9	\n\t"\
		"vmovaps	     (%%rbx),%%zmm6	\n\t	vmovaps	     (%%rbx,%%rdi),%%zmm14	\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm7	\n\t	vmovaps	0x040(%%rbx,%%rdi),%%zmm15	\n\t"\
		"vmovaps	%%zmm6,%%zmm4		\n\t	vmovaps	%%zmm14,%%zmm12				\n\t"\
		"vmovaps	%%zmm7,%%zmm5		\n\t	vmovaps	%%zmm15,%%zmm13				\n\t"\
		"vmulpd	%%zmm8,%%zmm4,%%zmm4	\n\t	vmulpd	%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vmulpd	%%zmm8,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm8 ,%%zmm13,%%zmm13		\n\t"\
	"vfnmadd231pd %%zmm9,%%zmm7,%%zmm4	\n\t vfnmadd231pd %%zmm9,%%zmm15,%%zmm12	\n\t"\
	"vfmadd231pd  %%zmm9,%%zmm6,%%zmm5	\n\t vfmadd231pd  %%zmm9,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	0x100(%%r10),%%zmm0	\n\t	vmovaps	0x300(%%r10),%%zmm8 		\n\t"/* Restore 2 */\
		"vmovaps	0x140(%%r10),%%zmm1	\n\t	vmovaps	0x340(%%r10),%%zmm9 		\n\t"\
		"vmovaps	%%zmm5,%%zmm7		\n\t	vmovaps	%%zmm13,%%zmm15				\n\t	 movq	%[__two],%%rsi	\n\t"\
		"vmovaps	%%zmm4,%%zmm6		\n\t	vmovaps	%%zmm12,%%zmm14				\n\t"\
		"vsubpd	%%zmm0,%%zmm4,%%zmm4	\n\t	vsubpd	%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vsubpd	%%zmm1,%%zmm5,%%zmm5	\n\t	vsubpd	%%zmm9 ,%%zmm13,%%zmm13		\n\t"\
		"vaddpd	%%zmm0,%%zmm6,%%zmm6	\n\t	vaddpd	%%zmm8 ,%%zmm14,%%zmm14		\n\t"\
		"vaddpd	%%zmm1,%%zmm7,%%zmm7	\n\t	vaddpd	%%zmm9 ,%%zmm15,%%zmm15		\n\t"\
		/* Finish radix-4 butterfly and store results: */\
		"vmovaps	     (%%r10),%%zmm0	\n\t	vmovaps	0x200(%%r10),%%zmm8 		\n\t"/* Restore 1 */\
		"vmovaps	0x040(%%r10),%%zmm1	\n\t	vmovaps	0x240(%%r10),%%zmm9 		\n\t"\
		"vsubpd	%%zmm6,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm14,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd	%%zmm5,%%zmm2,%%zmm2	\n\t	vsubpd	%%zmm13,%%zmm10,%%zmm10		\n\t"\
		"vsubpd	%%zmm7,%%zmm1,%%zmm1	\n\t	vsubpd	%%zmm15,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vsubpd	%%zmm12,%%zmm11,%%zmm11		\n\t"\
		"vmovaps	%%zmm0,0x100(%%r10)	\n\t	vmovaps	%%zmm8 ,0x300(%%r10)		\n\t"	/* 2.0, shared by both columns: */\
		"vmovaps	%%zmm2,0x080(%%r10)	\n\t	vmovaps	%%zmm10,0x280(%%r10)		\n\t	vmovaps	(%%rsi),%%zmm0	\n\t"\
		"vmovaps	%%zmm1,0x140(%%r10)	\n\t	vmovaps	%%zmm9 ,0x340(%%r10)		\n\t"\
		"vmovaps	%%zmm3,0x1c0(%%r10)	\n\t	vmovaps	%%zmm11,0x3c0(%%r10)		\n\t"\
	"vfmadd213pd 0x100(%%r10),%%zmm0,%%zmm6 \n\t	vfmadd213pd	%%zmm8 ,%%zmm0,%%zmm14	\n\t"\
	"vfmadd213pd      %%zmm2,%%zmm0,%%zmm5 \n\t	vfmadd213pd	%%zmm10,%%zmm0,%%zmm13	\n\t"\
	"vfmadd213pd      %%zmm1,%%zmm0,%%zmm7 \n\t	vfmadd213pd	%%zmm9 ,%%zmm0,%%zmm15	\n\t"\
	"vfmadd213pd      %%zmm3,%%zmm0,%%zmm4 \n\t	vfmadd213pd	%%zmm11,%%zmm0,%%zmm12	\n\t"\
		"vmovaps	%%zmm6,     (%%r10)	\n\t	vmovaps	%%zmm14,0x200(%%r10)		\n\t"/* don't need reload-from-mem of zmm8/0x300(%%r10) as we do in lcol, but 1 cycle faster with it. [!?] */\
		"vmovaps	%%zmm5,0x180(%%r10)	\n\t	vmovaps	%%zmm13,0x380(%%r10)		\n\t"\
		"vmovaps	%%zmm7,0x040(%%r10)	\n\t	vmovaps	%%zmm15,0x240(%%r10)		\n\t"\
		"vmovaps	%%zmm4,0x0c0(%%r10)	\n\t	vmovaps	%%zmm12,0x2c0(%%r10)		\n\t"\
	/*
		// Pass 2:
		k1 = 0x100;	// Set k1 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k2 = p4*8;	// Set k2 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		// 0x180 = literal-bytewise address offset between the pointers to the first and second DFT's twiddles = 0:
		i0 = r1; i1 = r1+16; i2 = r1+8; i3 = r1+24;
		o0 = (vec_dbl *)add0; o1 = (vec_dbl *)(add0+p2); o2 = (vec_dbl *)(add0+p1); o3 = (vec_dbl *)(add0+p3);
		c_tmp += 6;	// c2,A,6
		SSE2_RADIX_04_DIF_3TWIDDLE_X2(i0,i1,i2,i3,k1, two,c_tmp,0x180, o0,o1,o2,o3,k2)
	*/\
		/* Take advantage of fact that i1 = r0 + 16 already in %%r10 and o2 = add0+p1 in %%rax : */\
		/* i-offset = 0x100, so inline with r[a-d]x base-address offsets rather than adding via (r[a-d]x,rdi) */\
		"movq	%%r10,%%rbx				\n\t	movq	%%rax,%%r12			\n\t"/* rbx:i1 = r0+16; r12:o2 = add0+p1 */\
		"movq	%[__cc0],%%rsi 			\n\t	addq	$0x180,%%rsi 		\n\t"/* repoint rsi from two -> cc0, then increment for Pass 2 */\
		"leaq	 0x180(%%rsi),%%r8		\n\t	movslq	%[__p2],%%r11		\n\t"\
		"leaq	-0x400(%%rbx),%%rax		\n\t	movslq	%[__p4],%%r9		\n\t"\
		"leaq	-0x200(%%rbx),%%rcx		\n\t	movq	%[__add0]	,%%r10		\n\t"/* o0 */\
		"leaq	 0x200(%%rbx),%%rdx		\n\t	leaq	(%%r12,%%r11,8),%%r13	\n\t"/* o3 = (add0+p1)+p2 */\
		"vmovaps	     (%%rax),%%zmm0	\n\t	leaq	(%%r10,%%r11,8),%%r11	\n\t"/* o1 = (add0+p2) */\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t	shlq	$3,%%r9			\n\t"/* p4 in byte-offset form */\
		/* Do the p0,p2 combo: */\
		"vmovaps	     (%%rcx),%%zmm4	\n\t	vmovaps	0x100(%%rcx),%%zmm12		\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5	\n\t	vmovaps	0x140(%%rcx),%%zmm13		\n\t"\
		"vmovaps	     (%%rsi),%%zmm2	\n\t	vmovaps	     (%%r8),%%zmm10			\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm3	\n\t	vmovaps	0x040(%%r8),%%zmm11			\n\t"\
	/*"prefetcht1	%c[__pfetch_dist](%%r10,%%r9)\n\t"// ... + p4 */\
		"vmovaps	%%zmm4,%%zmm6		\n\t	vmovaps	%%zmm12,%%zmm14				\n\t"\
		"vmovaps	%%zmm5,%%zmm7		\n\t	vmovaps	%%zmm13,%%zmm15				\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t	vmovaps	0x100(%%rax),%%zmm8 		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t	vmovaps	0x140(%%rax),%%zmm9 		\n\t"\
		"vmulpd	%%zmm2,%%zmm4,%%zmm4	\n\t	vmulpd	%%zmm10,%%zmm12,%%zmm12		\n\t"\
		"vmulpd	%%zmm2,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm10,%%zmm13,%%zmm13		\n\t	vmovaps	%%zmm0,%%zmm2		\n\t"\
	"vfnmadd231pd %%zmm3,%%zmm7,%%zmm4	\n\t vfnmadd231pd %%zmm11,%%zmm15,%%zmm12	\n\t	vmovaps	%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd231pd  %%zmm3,%%zmm6,%%zmm5	\n\t vfmadd231pd  %%zmm11,%%zmm14,%%zmm13	\n\t	vmovaps	%%zmm1,%%zmm3		\n\t"\
		"vaddpd	%%zmm4,%%zmm0,%%zmm0	\n\t	vaddpd	%%zmm12,%%zmm8 ,%%zmm8 		\n\t	vmovaps	%%zmm9 ,%%zmm11		\n\t"\
		"vaddpd	%%zmm5,%%zmm1,%%zmm1	\n\t	vaddpd	%%zmm13,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd	%%zmm4,%%zmm2,%%zmm2	\n\t	vsubpd	%%zmm12,%%zmm10,%%zmm10		\n\t"\
		"vsubpd	%%zmm5,%%zmm3,%%zmm3	\n\t	vsubpd	%%zmm13,%%zmm11,%%zmm11		\n\t"\
		"vmovaps	%%zmm0,     (%%r10)	\n\t	vmovaps	%%zmm8 ,     (%%r10,%%r9)	\n\t"/* Spill 1: free up zmm0,1 */\
		"vmovaps	%%zmm1,0x040(%%r10)	\n\t	vmovaps	%%zmm9 ,0x040(%%r10,%%r9)	\n\t"/* Do the p1,3 combo: */\
		"vmovaps	0x100(%%rsi),%%zmm0	\n\t	vmovaps	0x100(%%r8),%%zmm8	 		\n\t"\
		"vmovaps	0x140(%%rsi),%%zmm1	\n\t	vmovaps	0x140(%%r8),%%zmm9 			\n\t"\
		"vmovaps	     (%%rdx),%%zmm6	\n\t	vmovaps	0x100(%%rdx),%%zmm14		\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm7	\n\t	vmovaps	0x140(%%rdx),%%zmm15		\n\t"\
		"vmovaps	%%zmm6,%%zmm4		\n\t	vmovaps	%%zmm14,%%zmm12				\n\t"\
		"vmovaps	%%zmm7,%%zmm5		\n\t	vmovaps	%%zmm15,%%zmm13				\n\t"\
		"vmulpd	%%zmm0,%%zmm4,%%zmm4	\n\t	vmulpd	%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vmulpd	%%zmm0,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm8 ,%%zmm13,%%zmm13		\n\t"\
	"vfnmadd231pd %%zmm1,%%zmm7,%%zmm4	\n\t vfnmadd231pd %%zmm9 ,%%zmm15,%%zmm12	\n\t"\
	"vfmadd231pd  %%zmm1,%%zmm6,%%zmm5	\n\t vfmadd231pd  %%zmm9 ,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	%%zmm5,0x040(%%r12)	\n\t	vmovaps	%%zmm13,0x040(%%r12,%%r9)	\n\t"/* Spill 2 */\
		"vmovaps	%%zmm4,     (%%r12)	\n\t	vmovaps	%%zmm12,     (%%r12,%%r9)	\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm0	\n\t	vmovaps	0x080(%%r8),%%zmm8	 		\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm1	\n\t	vmovaps	0x0c0(%%r8),%%zmm9 			\n\t"\
		"vmovaps	     (%%rbx),%%zmm6	\n\t	vmovaps	0x100(%%rbx),%%zmm14		\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm7	\n\t	vmovaps	0x140(%%rbx),%%zmm15		\n\t"\
		"vmovaps	%%zmm6,%%zmm4		\n\t	vmovaps	%%zmm14,%%zmm12				\n\t"\
		"vmovaps	%%zmm7,%%zmm5		\n\t	vmovaps	%%zmm15,%%zmm13				\n\t"\
		"vmulpd	%%zmm0,%%zmm4,%%zmm4	\n\t	vmulpd	%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vmulpd	%%zmm0,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm8 ,%%zmm13,%%zmm13		\n\t"\
	"vfnmadd231pd %%zmm1,%%zmm7,%%zmm4	\n\t vfnmadd231pd %%zmm9 ,%%zmm15,%%zmm12	\n\t"\
	"vfmadd231pd  %%zmm1,%%zmm6,%%zmm5	\n\t vfmadd231pd  %%zmm9 ,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	     (%%r12),%%zmm0	\n\t	vmovaps	     (%%r12,%%r9),%%zmm8	\n\t"/* Restore 2 */\
		"vmovaps	0x040(%%r12),%%zmm1	\n\t	vmovaps	0x040(%%r12,%%r9),%%zmm9	\n\t"\
		"vmovaps	%%zmm5,%%zmm7		\n\t	vmovaps	%%zmm13,%%zmm15				\n\t	 movq	%[__two],%%rsi	\n\t"\
		"vmovaps	%%zmm4,%%zmm6		\n\t	vmovaps	%%zmm12,%%zmm14				\n\t"\
		"vsubpd	%%zmm0,%%zmm4,%%zmm4	\n\t	vsubpd	%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vsubpd	%%zmm1,%%zmm5,%%zmm5	\n\t	vsubpd	%%zmm9 ,%%zmm13,%%zmm13		\n\t"\
		"vaddpd	%%zmm0,%%zmm6,%%zmm6	\n\t	vaddpd	%%zmm8 ,%%zmm14,%%zmm14		\n\t"\
		"vaddpd	%%zmm1,%%zmm7,%%zmm7	\n\t	vaddpd	%%zmm9 ,%%zmm15,%%zmm15		\n\t"\
		/* Finish radix-4 butterfly and store results: */\
		"vmovaps	     (%%r10),%%zmm0	\n\t	vmovaps	     (%%r10,%%r9),%%zmm8	\n\t"/* Restore 1 */\
		"vmovaps	0x040(%%r10),%%zmm1	\n\t	vmovaps	0x040(%%r10,%%r9),%%zmm9	\n\t"\
		"vsubpd	%%zmm6,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm14,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd	%%zmm5,%%zmm2,%%zmm2	\n\t	vsubpd	%%zmm13,%%zmm10,%%zmm10		\n\t"\
		"vsubpd	%%zmm7,%%zmm1,%%zmm1	\n\t	vsubpd	%%zmm15,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vsubpd	%%zmm12,%%zmm11,%%zmm11		\n\t"\
		"vmovaps	%%zmm0,     (%%r12)	\n\t	vmovaps	%%zmm8 ,     (%%r12,%%r9)	\n\t"	/* 2.0, shared by both columns ... moved +- until found best cycle count: */\
		"vmovaps	%%zmm2,     (%%r11)	\n\t	vmovaps	%%zmm10,     (%%r11,%%r9)	\n\t	vmovaps	(%%rsi),%%zmm0	\n\t"\
		"vmovaps	%%zmm1,0x040(%%r12)	\n\t	vmovaps	%%zmm9 ,0x040(%%r12,%%r9)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%r13)	\n\t	vmovaps	%%zmm11,0x040(%%r13,%%r9)	\n\t"\
	"vfmadd213pd (%%r12),%%zmm0,%%zmm6	\n\t	vfmadd213pd	%%zmm8 ,%%zmm0,%%zmm14	\n\t"\
	"vfmadd213pd %%zmm2 ,%%zmm0,%%zmm5	\n\t	vfmadd213pd	%%zmm10,%%zmm0,%%zmm13	\n\t"\
	"vfmadd213pd %%zmm1 ,%%zmm0,%%zmm7	\n\t	vfmadd213pd	%%zmm9 ,%%zmm0,%%zmm15	\n\t"\
	"vfmadd213pd %%zmm3 ,%%zmm0,%%zmm4	\n\t	vfmadd213pd	%%zmm11,%%zmm0,%%zmm12	\n\t"\
		"vmovaps	%%zmm6,     (%%r10)	\n\t	vmovaps	%%zmm14,     (%%r10,%%r9)	\n\t"/* don't need reload-from-mem of zmm8/0x300(%%r12,%%r9) as we do in lcol, but 1 cycle faster with it. [!?] */\
		"vmovaps	%%zmm5,     (%%r13)	\n\t	vmovaps	%%zmm13,     (%%r13,%%r9)	\n\t"\
		"vmovaps	%%zmm7,0x040(%%r10)	\n\t	vmovaps	%%zmm15,0x040(%%r10,%%r9)	\n\t"\
		"vmovaps	%%zmm4,0x040(%%r11)	\n\t	vmovaps	%%zmm12,0x040(%%r11,%%r9)	\n\t"\
	/*
		i0 += 2; i1 += 2; i2 += 2; i3 += 2;
		o0 = (vec_dbl *)(add0+p8); o1 = (vec_dbl *)(add0+p8+p2); o2 = (vec_dbl *)(add0+p8+p1); o3 = (vec_dbl *)(add0+p8+p3);
		c_tmp += 12;	// c5,D,3
		SSE2_RADIX_04_DIF_3TWIDDLE_X2(i0,i1,i2,i3,k1, two,c_tmp,0x180, o0,o1,o2,o3,k2)
	*/\
		/* Take advantage of fact that i1 = r0 + 16 already in %%r10 and o2 = add0+p1 in %%rax : */\
		"movq	%[__cc0],%%rsi 			\n\t	addq	$0x480,%%rsi 		\n\t"/* repoint rsi from two -> cc0, then increment for Pass, set 2 */\
		"addq	$0x080,%%rax			\n\t	movslq	%[__p8],%%r8		\n\t"\
		"addq	$0x080,%%rbx			\n\t"	/* r9 still has p4<<3 */\
		"addq	$0x080,%%rcx			\n\t	leaq	(%%r10,%%r8,8),%%r10	\n\t"/* o0 =  add0     +p8 */\
		"addq	$0x080,%%rdx			\n\t	leaq	(%%r11,%%r8,8),%%r11	\n\t"/* o1 = (add0+p2)+p8 */\
		"vmovaps	     (%%rax),%%zmm0	\n\t	leaq	(%%r12,%%r8,8),%%r12	\n\t"/* o2 = (add0+p1)+p8 */\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t	leaq	(%%r13,%%r8,8),%%r13	\n\t"/* o3 = (add0+p3)+p8 */\
		"leaq	0x180(%%rsi),%%r8		\n\t"\
		/* Do the p0,p2 combo: */\
		"vmovaps	     (%%rcx),%%zmm4	\n\t	vmovaps	0x100(%%rcx),%%zmm12		\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5	\n\t	vmovaps	0x140(%%rcx),%%zmm13		\n\t"\
		"vmovaps	     (%%rsi),%%zmm2	\n\t	vmovaps	     (%%r8),%%zmm10			\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm3	\n\t	vmovaps	0x040(%%r8),%%zmm11			\n\t"\
	/*"prefetcht1	%c[__pfetch_dist](%%r10,%%r9)\n\t"// ... + pC */\
		"vmovaps	%%zmm4,%%zmm6		\n\t	vmovaps	%%zmm12,%%zmm14				\n\t"\
		"vmovaps	%%zmm5,%%zmm7		\n\t	vmovaps	%%zmm13,%%zmm15				\n\t"\
		"vmovaps	     (%%rax),%%zmm0	\n\t	vmovaps	0x100(%%rax),%%zmm8 		\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t	vmovaps	0x140(%%rax),%%zmm9 		\n\t"\
		"vmulpd	%%zmm2,%%zmm4,%%zmm4	\n\t	vmulpd	%%zmm10,%%zmm12,%%zmm12		\n\t"\
		"vmulpd	%%zmm2,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm10,%%zmm13,%%zmm13		\n\t	vmovaps	%%zmm0,%%zmm2		\n\t"\
	"vfnmadd231pd %%zmm3,%%zmm7,%%zmm4	\n\t vfnmadd231pd %%zmm11,%%zmm15,%%zmm12	\n\t	vmovaps	%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd231pd  %%zmm3,%%zmm6,%%zmm5	\n\t vfmadd231pd  %%zmm11,%%zmm14,%%zmm13	\n\t	vmovaps	%%zmm1,%%zmm3		\n\t"\
		"vaddpd	%%zmm4,%%zmm0,%%zmm0	\n\t	vaddpd	%%zmm12,%%zmm8 ,%%zmm8 		\n\t	vmovaps	%%zmm9 ,%%zmm11		\n\t"\
		"vaddpd	%%zmm5,%%zmm1,%%zmm1	\n\t	vaddpd	%%zmm13,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd	%%zmm4,%%zmm2,%%zmm2	\n\t	vsubpd	%%zmm12,%%zmm10,%%zmm10		\n\t"\
		"vsubpd	%%zmm5,%%zmm3,%%zmm3	\n\t	vsubpd	%%zmm13,%%zmm11,%%zmm11		\n\t"\
		"vmovaps	%%zmm0,     (%%r10)	\n\t	vmovaps	%%zmm8 ,     (%%r10,%%r9)	\n\t"/* Spill 1: free up zmm0,1 */\
		"vmovaps	%%zmm1,0x040(%%r10)	\n\t	vmovaps	%%zmm9 ,0x040(%%r10,%%r9)	\n\t"/* Do the p1,3 combo: */\
		"vmovaps	0x100(%%rsi),%%zmm0	\n\t	vmovaps	0x100(%%r8),%%zmm8	 		\n\t"\
		"vmovaps	0x140(%%rsi),%%zmm1	\n\t	vmovaps	0x140(%%r8),%%zmm9 			\n\t"\
		"vmovaps	     (%%rdx),%%zmm6	\n\t	vmovaps	0x100(%%rdx),%%zmm14		\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm7	\n\t	vmovaps	0x140(%%rdx),%%zmm15		\n\t"\
		"vmovaps	%%zmm6,%%zmm4		\n\t	vmovaps	%%zmm14,%%zmm12				\n\t"\
		"vmovaps	%%zmm7,%%zmm5		\n\t	vmovaps	%%zmm15,%%zmm13				\n\t"\
		"vmulpd	%%zmm0,%%zmm4,%%zmm4	\n\t	vmulpd	%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vmulpd	%%zmm0,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm8 ,%%zmm13,%%zmm13		\n\t"\
	"vfnmadd231pd %%zmm1,%%zmm7,%%zmm4	\n\t vfnmadd231pd %%zmm9 ,%%zmm15,%%zmm12	\n\t"\
	"vfmadd231pd  %%zmm1,%%zmm6,%%zmm5	\n\t vfmadd231pd  %%zmm9 ,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	%%zmm5,0x040(%%r12)	\n\t	vmovaps	%%zmm13,0x040(%%r12,%%r9)	\n\t"/* Spill 2 */\
		"vmovaps	%%zmm4,     (%%r12)	\n\t	vmovaps	%%zmm12,     (%%r12,%%r9)	\n\t"\
		"vmovaps	0x080(%%rsi),%%zmm0	\n\t	vmovaps	0x080(%%r8),%%zmm8	 		\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm1	\n\t	vmovaps	0x0c0(%%r8),%%zmm9 			\n\t"\
		"vmovaps	     (%%rbx),%%zmm6	\n\t	vmovaps	0x100(%%rbx),%%zmm14		\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm7	\n\t	vmovaps	0x140(%%rbx),%%zmm15		\n\t"\
		"vmovaps	%%zmm6,%%zmm4		\n\t	vmovaps	%%zmm14,%%zmm12				\n\t"\
		"vmovaps	%%zmm7,%%zmm5		\n\t	vmovaps	%%zmm15,%%zmm13				\n\t"\
		"vmulpd	%%zmm0,%%zmm4,%%zmm4	\n\t	vmulpd	%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vmulpd	%%zmm0,%%zmm5,%%zmm5	\n\t	vmulpd	%%zmm8 ,%%zmm13,%%zmm13		\n\t"\
	"vfnmadd231pd %%zmm1,%%zmm7,%%zmm4	\n\t vfnmadd231pd %%zmm9 ,%%zmm15,%%zmm12	\n\t"\
	"vfmadd231pd  %%zmm1,%%zmm6,%%zmm5	\n\t vfmadd231pd  %%zmm9 ,%%zmm14,%%zmm13	\n\t"\
		"vmovaps	     (%%r12),%%zmm0	\n\t	vmovaps	     (%%r12,%%r9),%%zmm8	\n\t"/* Restore 2 */\
		"vmovaps	0x040(%%r12),%%zmm1	\n\t	vmovaps	0x040(%%r12,%%r9),%%zmm9	\n\t"\
		"vmovaps	%%zmm5,%%zmm7		\n\t	vmovaps	%%zmm13,%%zmm15				\n\t	 movq	%[__two],%%rsi	\n\t"\
		"vmovaps	%%zmm4,%%zmm6		\n\t	vmovaps	%%zmm12,%%zmm14				\n\t"\
		"vsubpd	%%zmm0,%%zmm4,%%zmm4	\n\t	vsubpd	%%zmm8 ,%%zmm12,%%zmm12		\n\t"\
		"vsubpd	%%zmm1,%%zmm5,%%zmm5	\n\t	vsubpd	%%zmm9 ,%%zmm13,%%zmm13		\n\t"\
		"vaddpd	%%zmm0,%%zmm6,%%zmm6	\n\t	vaddpd	%%zmm8 ,%%zmm14,%%zmm14		\n\t"\
		"vaddpd	%%zmm1,%%zmm7,%%zmm7	\n\t	vaddpd	%%zmm9 ,%%zmm15,%%zmm15		\n\t"\
		/* Finish radix-4 butterfly and store results: */\
		"vmovaps	     (%%r10),%%zmm0	\n\t	vmovaps	     (%%r10,%%r9),%%zmm8	\n\t"/* Restore 1 */\
		"vmovaps	0x040(%%r10),%%zmm1	\n\t	vmovaps	0x040(%%r10,%%r9),%%zmm9	\n\t"\
		"vsubpd	%%zmm6,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm14,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd	%%zmm5,%%zmm2,%%zmm2	\n\t	vsubpd	%%zmm13,%%zmm10,%%zmm10		\n\t"\
		"vsubpd	%%zmm7,%%zmm1,%%zmm1	\n\t	vsubpd	%%zmm15,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vsubpd	%%zmm12,%%zmm11,%%zmm11		\n\t"\
		"vmovaps	%%zmm0,     (%%r12)	\n\t	vmovaps	%%zmm8 ,     (%%r12,%%r9)	\n\t"	/* 2.0, shared by both columns ... moved +- until found best cycle count: */\
		"vmovaps	%%zmm2,     (%%r11)	\n\t	vmovaps	%%zmm10,     (%%r11,%%r9)	\n\t	vmovaps	(%%rsi),%%zmm0	\n\t"\
		"vmovaps	%%zmm1,0x040(%%r12)	\n\t	vmovaps	%%zmm9 ,0x040(%%r12,%%r9)	\n\t"\
		"vmovaps	%%zmm3,0x040(%%r13)	\n\t	vmovaps	%%zmm11,0x040(%%r13,%%r9)	\n\t"\
	"vfmadd213pd (%%r12),%%zmm0,%%zmm6	\n\t	vfmadd213pd	%%zmm8 ,%%zmm0,%%zmm14	\n\t"\
	"vfmadd213pd %%zmm2 ,%%zmm0,%%zmm5	\n\t	vfmadd213pd	%%zmm10,%%zmm0,%%zmm13	\n\t"\
	"vfmadd213pd %%zmm1 ,%%zmm0,%%zmm7	\n\t	vfmadd213pd	%%zmm9 ,%%zmm0,%%zmm15	\n\t"\
	"vfmadd213pd %%zmm3 ,%%zmm0,%%zmm4	\n\t	vfmadd213pd	%%zmm11,%%zmm0,%%zmm12	\n\t"\
		"vmovaps	%%zmm6,     (%%r10)	\n\t	vmovaps	%%zmm14,     (%%r10,%%r9)	\n\t"/* don't need reload-from-mem of zmm8/0x300(%%r12,%%r9) as we do in lcol, but 1 cycle faster with it. [!?] */\
		"vmovaps	%%zmm5,     (%%r13)	\n\t	vmovaps	%%zmm13,     (%%r13,%%r9)	\n\t"\
		"vmovaps	%%zmm7,0x040(%%r10)	\n\t	vmovaps	%%zmm15,0x040(%%r10,%%r9)	\n\t"\
		"vmovaps	%%zmm4,0x040(%%r11)	\n\t	vmovaps	%%zmm12,0x040(%%r11,%%r9)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r0] "m" (Xr0)\
		 ,[__two] "m" (Xtwo)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	/***********************************************************************/
	/********************** Radix-16 DIT macros ****************************/
	/***********************************************************************/
	#define RADIX16_COMPUTE_FMA_SINCOS_DIT(Xadd0,Xone)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax			\n\t"/* base-ptr to 32-byte-aligned 128-byte-large memchunk 16 doubles needing iterative inversion */\
		"movq	%[__one] ,%%rbx			\n\t"/* 1.0 in 4-fold-double form */\
		"vmovaps	    (%%rbx),%%ymm14	\n\t"/* 1.0 */\
		"vaddpd	%%ymm14,%%ymm14,%%ymm14	\n\t"/* 2.0 */\
		"vmovapd	    (%%rax),%%ymm4	\n\t"/* load inputs-to-be-inverted (d's) into 4 AVX registers */\
		"vmovapd	0x20(%%rax),%%ymm5	\n\t"\
		"vmovapd	0x40(%%rax),%%ymm6	\n\t"\
		"vmovapd	0x60(%%rax),%%ymm7	\n\t"\
		"vrcp14pd	%%zmm4,%%zmm0	\n\t"	/* ainv := approx 1/d to >= 14 bits of precision */\
		"vrcp14pd	%%zmm5,%%zmm1	\n\t"	/* NOTE: AVX-512F requires full-width zmm-regs here, but only use lower 256 bits */\
		"vrcp14pd	%%zmm6,%%zmm2	\n\t"\
		"vrcp14pd	%%zmm7,%%zmm3	\n\t"\
		/* 1st NR iteration gives ~28 bits of precision: */\
		"vmovaps	%%ymm0,%%ymm8 	\n\t"	/* make a copy of ainv */\
		"vmovaps	%%ymm1,%%ymm9 	\n\t"\
		"vmovaps	%%ymm2,%%ymm10	\n\t"\
		"vmovaps	%%ymm3,%%ymm11	\n\t"\
		"vfnmadd132pd	%%ymm4,%%ymm14,%%ymm0 	\n\t"	/* 2 - d*ainv, overwrites ainv */\
		"vfnmadd132pd	%%ymm5,%%ymm14,%%ymm1 	\n\t"\
		"vfnmadd132pd	%%ymm6,%%ymm14,%%ymm2 	\n\t"\
		"vfnmadd132pd	%%ymm7,%%ymm14,%%ymm3 	\n\t"\
		"vmulpd		%%ymm0,%%ymm8 ,%%ymm0	\n\t"	/* ainv*(2 - d*ainv) = 1/d accurate to ~28 bits */\
		"vmulpd		%%ymm1,%%ymm9 ,%%ymm1	\n\t"\
		"vmulpd		%%ymm2,%%ymm10,%%ymm2	\n\t"\
		"vmulpd		%%ymm3,%%ymm11,%%ymm3	\n\t"\
		/* 2nd NR iteration gives the maximal ~53 bits of precision: */\
		"vmovaps	%%ymm0,%%ymm8 	\n\t"	/* make a copy of ainv */\
		"vmovaps	%%ymm1,%%ymm9 	\n\t"\
		"vmovaps	%%ymm2,%%ymm10	\n\t"\
		"vmovaps	%%ymm3,%%ymm11	\n\t"\
		"vfnmadd132pd	%%ymm4,%%ymm14,%%ymm0 	\n\t"	/* 2 - d*ainv, overwrites ainv */\
		"vfnmadd132pd	%%ymm5,%%ymm14,%%ymm1 	\n\t"\
		"vfnmadd132pd	%%ymm6,%%ymm14,%%ymm2 	\n\t"\
		"vfnmadd132pd	%%ymm7,%%ymm14,%%ymm3 	\n\t"\
		"vmulpd		%%ymm0,%%ymm8 ,%%ymm0	\n\t"	/* ainv*(2 - d*ainv) = 1/d accurate to ~53 bits */\
		"vmulpd		%%ymm1,%%ymm9 ,%%ymm1	\n\t"\
		"vmulpd		%%ymm2,%%ymm10,%%ymm2	\n\t"\
		"vmulpd		%%ymm3,%%ymm11,%%ymm3	\n\t"\
		/* Multiply by the sine terms to get quotients in ymm4-7: */\
		"vmulpd	0x80(%%rax),%%ymm0,%%ymm4	\n\t"\
		"vmulpd	0xa0(%%rax),%%ymm1,%%ymm5	\n\t"\
		"vmulpd	0xc0(%%rax),%%ymm2,%%ymm6	\n\t"\
		"vmulpd	0xe0(%%rax),%%ymm3,%%ymm7	\n\t"\
		/* Outputs overwrite sines: */\
		"vmovaps	%%ymm4,0x80(%%rax)	\n\t"\
		"vmovaps	%%ymm5,0xa0(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0xc0(%%rax)	\n\t"\
		"vmovaps	%%ymm7,0xe0(%%rax)	\n\t"\
		/* Now propagate the scalar doubles in the 8 ymm-sized memlocs referenced above to their final 4x-copied locs. Scalar data laid out as */\
		/* add0 + 0x[  0,  8, 10, 18]: [c0 ,c1 ,c2 ,c3 ] Cosines: */\
		/* add0 + 0x[ 20, 28, 30, 38]: [c4 ,c5 ,c6 ,c7 ] */\
		/* add0 + 0x[ 40, 48, 50, 58]: [c8 ,c9 ,cA ,cB ] */\
		/* add0 + 0x[ 60, 68, 70, 78]: [cC ,cD ,cE ,cF ] */\
		/* add0 + 0x[ 80, 88, 90, 98]: [r0 ,r1 ,r2 ,r3 ] Tangents: */\
		/* add0 + 0x[ a0, a8, b0, b8]: [r4 ,r5 ,r6 ,r7 ] */\
		/* add0 + 0x[ c0, c8, d0, d8]: [r8 ,r9 ,rA ,rB ] */\
		/* add0 + 0x[ e0, e8, f0, f8]: [rC ,rD ,rE ,rF ] */\
		/* We start at the upper end of the memory chunk allocated for sincos data and proceed downward, so the 8 memory */\
		/* slots inited with scalar data above are only overwritten after their scalar-data contents have been copied: */\
		"vbroadcastsd	0x0f8(%%rax),%%zmm15	/* rF */\n\t"/* Allow several cycles for each vbroadcast to complete before writing result */\
		"vbroadcastsd	0x078(%%rax),%%zmm14	/* cF */\n\t"\
		"vbroadcastsd	0x0d8(%%rax),%%zmm13	/* rB */\n\t"\
		"vbroadcastsd	0x058(%%rax),%%zmm12	/* cB */\n\t"\
		"vbroadcastsd	0x0b8(%%rax),%%zmm11	/* r7 */\n\t	vmovaps	%%zmm15,0x840(%%rax)	\n\t"/* s15 = cc0 + 0x21;	__rF = s15/c15		*/\
		"vbroadcastsd	0x038(%%rax),%%zmm10	/* c7 */\n\t	vmovaps	%%zmm14,0x800(%%rax)	\n\t"/* c15 = cc0 + 0x20;	__cF				*/\
		"vbroadcastsd	0x098(%%rax),%%zmm9 	/* r3 */\n\t	vmovaps	%%zmm13,0x7c0(%%rax)	\n\t"/* s11 = cc0 + 0x1f;	__rB = s11/c11		*/\
		"vbroadcastsd	0x018(%%rax),%%zmm8 	/* c3 */\n\t	vmovaps	%%zmm12,0x780(%%rax)	\n\t"/* c11 = cc0 + 0x1e;	__cB				*/\
		"vbroadcastsd	0x0e8(%%rax),%%zmm7 	/* rD */\n\t	vmovaps	%%zmm11,0x740(%%rax)	\n\t"/* s7  = cc0 + 0x1d;	__r7 = s7 /c7 		*/\
		"vbroadcastsd	0x068(%%rax),%%zmm6 	/* cD */\n\t	vmovaps	%%zmm10,0x700(%%rax)	\n\t"/* c7  = cc0 + 0x1c;	__c7				*/\
		"vbroadcastsd	0x0c8(%%rax),%%zmm5 	/* r9 */\n\t	vmovaps	%%zmm9 ,0x6c0(%%rax)	\n\t"/* s3  = cc0 + 0x1b;	__r3 = s3 /c3 		*/\
		"vbroadcastsd	0x048(%%rax),%%zmm4 	/* c9 */\n\t	vmovaps	%%zmm8 ,0x680(%%rax)	\n\t"/* c3  = cc0 + 0x1a;	__c3				*/\
		"vbroadcastsd	0x0a8(%%rax),%%zmm3 	/* r5 */\n\t	vmovaps	%%zmm7 ,0x640(%%rax)	\n\t"/* s13 = cc0 + 0x19;	__rD = s13/c13		*/\
		"vbroadcastsd	0x028(%%rax),%%zmm2 	/* c5 */\n\t	vmovaps	%%zmm6 ,0x600(%%rax)	\n\t"/* c13 = cc0 + 0x18;	__cD				*/\
		"vbroadcastsd	0x088(%%rax),%%zmm1 	/* r1 */\n\t	vmovaps	%%zmm5 ,0x5c0(%%rax)	\n\t"/* s9  = cc0 + 0x17;	__r9 = s9 /c9 		*/\
		"vbroadcastsd	0x008(%%rax),%%zmm0 	/* c1 */\n\t	vmovaps	%%zmm4 ,0x580(%%rax)	\n\t"/* c9  = cc0 + 0x16;	__c9				*/\
		"vbroadcastsd	0x0f0(%%rax),%%zmm15	/* rE */\n\t	vmovaps	%%zmm3 ,0x540(%%rax)	\n\t"/* s5  = cc0 + 0x15;	__r5 = s5 /c5 		*/\
		"vbroadcastsd	0x070(%%rax),%%zmm14	/* cE */\n\t	vmovaps	%%zmm2 ,0x500(%%rax)	\n\t"/* c5  = cc0 + 0x14;	__c5				*/\
		"vbroadcastsd	0x0d0(%%rax),%%zmm13	/* rA */\n\t	vmovaps	%%zmm1 ,0x4c0(%%rax)	\n\t"/* s1  = cc0 + 0x13;	__r1 = s1 /c1 		*/\
		"vbroadcastsd	0x050(%%rax),%%zmm12	/* cA */\n\t	vmovaps	%%zmm0 ,0x480(%%rax)	\n\t"/* s1  = cc0 + 0x12;	__c1				*/\
		"vbroadcastsd	0x0b0(%%rax),%%zmm11	/* r6 */\n\t	vmovaps	%%zmm15,0x440(%%rax)	\n\t"/* s14 = cc0 + 0x11;	__rE = s14/c14		*/\
		"vbroadcastsd	0x030(%%rax),%%zmm10	/* c6 */\n\t	vmovaps	%%zmm14,0x400(%%rax)	\n\t"/* c14 = cc0 + 0x10;	__cE				*/\
		"vbroadcastsd	0x090(%%rax),%%zmm9 	/* r2 */\n\t	vmovaps	%%zmm13,0x3c0(%%rax)	\n\t"/* s10 = cc0 + 0x0f;	__rA = s10/c10		*/\
		"vbroadcastsd	0x010(%%rax),%%zmm8 	/* c2 */\n\t	vmovaps	%%zmm12,0x380(%%rax)	\n\t"/* c10 = cc0 + 0x0e;	__cA				*/\
		"vbroadcastsd	0x0e0(%%rax),%%zmm7 	/* rC */\n\t	vmovaps	%%zmm11,0x340(%%rax)	\n\t"/* s6  = cc0 + 0x0d;	__r6 = s6 /c6 		*/\
		"vbroadcastsd	0x060(%%rax),%%zmm6 	/* cC */\n\t	vmovaps	%%zmm10,0x300(%%rax)	\n\t"/* c6  = cc0 + 0x0c;	__c6				*/\
		"vbroadcastsd	0x0c0(%%rax),%%zmm5 	/* r8 */\n\t	vmovaps	%%zmm9 ,0x2c0(%%rax)	\n\t"/* s2  = cc0 + 0x0b;	__r2 = s2 /c2 		*/\
		"vbroadcastsd	0x040(%%rax),%%zmm4 	/* c8 */\n\t	vmovaps	%%zmm8 ,0x280(%%rax)	\n\t"/* c2  = cc0 + 0x0a;	__c2				*/\
		"vbroadcastsd	0x0a0(%%rax),%%zmm3 	/* r4 */\n\t	vmovaps	%%zmm7 ,0x240(%%rax)	\n\t"/* s12 = cc0 + 0x09;	__rC = s12/c12		*/\
		"vbroadcastsd	0x020(%%rax),%%zmm2 	/* c4 */\n\t	vmovaps	%%zmm6 ,0x200(%%rax)	\n\t"/* c12 = cc0 + 0x08;	__cC				*/\
		"vbroadcastsd	0x080(%%rax),%%zmm1		/*__t */\n\t	vmovaps	%%zmm5 ,0x1c0(%%rax)	\n\t"/* s8  = cc0 + 0x07;	__r8 = s8 /c8 		*/\
		"vbroadcastsd	     (%%rax),%%zmm0		/*__c */\n\t	vmovaps	%%zmm4 ,0x180(%%rax)	\n\t"/* c8  = cc0 + 0x06;	__c8 [unchanged]	*/\
		"														vmovaps	%%zmm3 ,0x140(%%rax)	\n\t"/* s4  = cc0 + 0x05;	__r4 = s4 /c4 		*/\
		"														vmovaps	%%zmm2 ,0x100(%%rax)	\n\t"/* c4  = cc0 + 0x04;	__c4 [unchanged]	*/\
		"														vmovaps	%%zmm1 ,0x040(%%rax)	\n\t"/* ss0 = cc0 + 0x01;	__sc = __s/__c		*/\
		"														vmovaps	%%zmm0 ,     (%%rax)	\n\t"/* cc0 = cc0 + 0x00;	__c					*/\
		"vmovaps	    (%%rbx),%%zmm14	\n\t"/* 1.0 */\
		"vpxorq	%%zmm0,%%zmm0,%%zmm0	\n\t"/* 0.0 */\
		"vmovaps	%%zmm14,0x080(%%rax)	\n\t"/* c0 = cc0 + 0x02;	don't use in DFT but init = 1.0 */\
		"vmovaps	%%zmm0 ,0x0c0(%%rax)	\n\t"/* s0 = cc0 + 0x03;	don't use in DFT but init = 0.0 */\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__one] "m" (Xone)\
		: "cc","memory","rax","rbx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	/* Note no cosine ratios needed in DIT version. */\
	#define SSE2_RADIX16_DIT_TWIDDLE_1(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr1,Xcc0,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/* Gather the needed data and do first set of four length-4 transforms. Each complex radix-4 block needs 11 registers: */\
	/* 8 t-temps mapped to zmm0-7, plus rt,it (zmm8,9), plus const 1.0 (zmm10), though we add 2 more tmps re,im (zmm11,12) */\
	/* to allow independent Re/Im subsections to overlap w/o introducing false dependencies. */\
	/*...Block 1: */\
		"movslq	%[__p4],%%r10	\n\t"/* This serves as main-array stride-between radix-4 blocks */\
		"movq	%[__add0],%%rax				\n\t"\
	"movslq	%[__pfetch_addr1],%%r13	\n\t"	/* Prefetch base-index offset cycles among add0 + p0,1,2,3 on successive macro calls */\
	"leaq	(%%rax,%%r13,8),%%r13	\n\t"	/* [base-address + data-fetch-ahead index] */\
		"movslq	%[__p1],%%rbx	\n\t"\
		"movslq	%[__p2],%%rcx	\n\t"\
		"movslq	%[__p3],%%rdx	\n\t"\
		"movq	%[__cc0],%%rsi 	\n\t"\
		"vmovaps 0x880(%%rsi),%%zmm10	\n\t"/* cc0 + 0x44 = __two; Actually holds 1.0 in AVX2 mode */\
		"leaq	(%%rax,%%rbx,8),%%rbx	\n\t"/* add0+p1 */\
		"leaq	(%%rax,%%rcx,8),%%rcx	\n\t"/* add0+p2 */\
		"leaq	(%%rax,%%rdx,8),%%rdx	\n\t"/* add0+p3 */\
		"movq	%[__r1],%%rdi	\n\t"/* ptr to local-mem [t1] */\
	"xorq	%%r14,%%r14	\n\t"/* Zero r14 and include in LEA to eliminate "scale factor of 8 without an index register" assembler warnings */\
	"leaq	(%%r14,%%r10,8),%%r14	\n\t"/* Save a copy of p4 in ptr-offset form. Will prefetch from [base-address + data-fetch-ahead index] + [0,p4,p8,p12] on each macro call. */\
	"prefetcht1	%c[__pfetch_dist](%%r13)\n\t"/* [base-address + data-fetch-ahead index] + p0 */\
		"vmovaps		     (%%rax),%%zmm0 	\n\t	vmovaps		     (%%rbx),%%zmm8 	\n\t"/*	t1,rt =__A0,1r; */\
		"vmovaps		     (%%rcx),%%zmm4 	\n\t	vmovaps		     (%%rdx),%%zmm9 	\n\t"/*	t5,it =__A2,3r; */\
		"vsubpd			%%zmm8 ,%%zmm0 ,%%zmm2 	\n\t"/*	FNMA231(1.0,rt ,t3 ); */\
		" vfmadd231pd	%%zmm10,%%zmm8 ,%%zmm0 	\n\t"/*	 FMA231(1.0,rt ,t1 ); */\
		"vsubpd			%%zmm9 ,%%zmm4 ,%%zmm7 	\n\t"/*	FNMA231(1.0,it ,t8 ); */\
		" vfmadd231pd	%%zmm10,%%zmm9 ,%%zmm4 	\n\t"/*	 FMA231(1.0,it ,t5 ); */\
		"vmovaps		0x040(%%rax),%%zmm1 	\n\t	vmovaps		0x040(%%rbx),%%zmm11	\n\t"/*	t2,re =__A0,1i; */\
		"vmovaps		0x040(%%rcx),%%zmm5 	\n\t	vmovaps		0x040(%%rdx),%%zmm12	\n\t"/*	t6,im =__A2,3i; */\
		"vsubpd			%%zmm11,%%zmm1 ,%%zmm3 	\n\t"/*	FNMA231(1.0,re ,t4 ); */\
		" vfmadd231pd	%%zmm10,%%zmm11,%%zmm1 	\n\t"/*	 FMA231(1.0,re ,t2 ); */\
		"vsubpd			%%zmm12,%%zmm5 ,%%zmm6 	\n\t"/*	FNMA231(1.0,im ,t7 ); */\
		" vfmadd231pd	%%zmm10,%%zmm12,%%zmm5 	\n\t"/*	 FMA231(1.0,im ,t6 ); */\
		"vmovaps		%%zmm4,%%zmm8			\n\t	vmovaps		%%zmm2,%%zmm9	\n\t"/*	rt = t5; it = t3; */\
		/*** REMEMBER - For FMA132, rightmost 2 operands in ASM are reverse-ordered w.r.to the prototyping code in the comments! ***/\
		"vfnmadd132pd	%%zmm10,%%zmm0 ,%%zmm4 	\n\t"/*	FNMA132(1.0,t5 ,t1 ); */\
		" vfmadd231pd	%%zmm10,%%zmm8 ,%%zmm0 	\n\t"/*	 FMA231(1.0,rt ,t1 ); */\
		" vfmadd231pd	%%zmm10,%%zmm6 ,%%zmm2 	\n\t"/*	 FMA231(1.0,t7 ,t3 ); */\
		"vfnmadd132pd	%%zmm10,%%zmm9 ,%%zmm6 	\n\t"/*	FNMA132(1.0,t7 ,it ); */\
		"vmovaps		%%zmm5,%%zmm11			\n\t	vmovaps		%%zmm3,%%zmm12		\n\t"/*	re = t6; im = t4; */\
		"vfnmadd132pd	%%zmm10,%%zmm1 ,%%zmm5 	\n\t	vmovaps	%%zmm4 ,0x100(%%rdi)	\n\t"/* FNMA132(1.0,t6 ,t2 ); write t5 */\
		" vfmadd231pd	%%zmm10,%%zmm11,%%zmm1 	\n\t	vmovaps	%%zmm0 ,     (%%rdi)	\n\t"/*  FMA231(1.0,re ,t2 ); write t1 */\
		"vfnmadd231pd	%%zmm10,%%zmm7 ,%%zmm3 	\n\t	vmovaps	%%zmm2 ,0x080(%%rdi)	\n\t"/* FNMA231(1.0,t8 ,t4 ); write t3 */\
		" vfmadd132pd	%%zmm10,%%zmm12,%%zmm7 	\n\t	vmovaps	%%zmm6 ,0x180(%%rdi)	\n\t"/*  FMA132(1.0,t8 ,im ); write t7 */\
		"vmovaps		%%zmm5 ,0x140(%%rdi)  	\n\t"/* write t6 */\
		"vmovaps		%%zmm1 ,0x040(%%rdi)  	\n\t"/* write t2 */\
		"vmovaps		%%zmm3 ,0x0c0(%%rdi)  	\n\t"/* write t4 */\
		"vmovaps		%%zmm7 ,0x1c0(%%rdi)  	\n\t"/* write t8 */\
		/**/\
	/*...Block 2: __A4-7, t9-16 */\
		"addq	$0x200,%%rdi	\n\t"/* ptr to t9 */\
		"vmovaps	     (%%rax,%%r10,8),%%zmm0 \n\t	vmovaps	     (%%rbx,%%r10,8),%%zmm8 \n\t"\
		"vmovaps	     (%%rcx,%%r10,8),%%zmm4 \n\t	vmovaps	     (%%rdx,%%r10,8),%%zmm9 \n\t"\
		"vsubpd			%%zmm8 ,%%zmm0 ,%%zmm2 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm8 ,%%zmm0 	\n\t"\
		"vsubpd			%%zmm9 ,%%zmm4 ,%%zmm7 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm9 ,%%zmm4 	\n\t"\
		"vmovaps	0x040(%%rax,%%r10,8),%%zmm1 \n\t	vmovaps	0x040(%%rbx,%%r10,8),%%zmm11\n\t"\
		"vmovaps	0x040(%%rcx,%%r10,8),%%zmm5 \n\t	vmovaps	0x040(%%rdx,%%r10,8),%%zmm12\n\t"\
		"vsubpd			%%zmm11,%%zmm1 ,%%zmm3 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm11,%%zmm1 	\n\t"\
		"vsubpd			%%zmm12,%%zmm5 ,%%zmm6 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm12,%%zmm5 	\n\t"\
		"vmovaps		%%zmm4,%%zmm8			\n\t	vmovaps		%%zmm2,%%zmm9	\n\t"\
		"vfnmadd132pd	%%zmm10,%%zmm0 ,%%zmm4 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm8 ,%%zmm0 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm6 ,%%zmm2 	\n\t"\
		"vfnmadd132pd	%%zmm10,%%zmm9 ,%%zmm6 	\n\t"\
		"vmovaps		%%zmm5,%%zmm11			\n\t	vmovaps		%%zmm3,%%zmm12	\n\t"\
		"vfnmadd132pd	%%zmm10,%%zmm1 ,%%zmm5 	\n\t	vmovaps		%%zmm4 ,0x100(%%rdi) 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm11,%%zmm1 	\n\t	vmovaps		%%zmm0 ,     (%%rdi) 	\n\t"\
		"vfnmadd231pd	%%zmm10,%%zmm7 ,%%zmm3 	\n\t	vmovaps		%%zmm2 ,0x080(%%rdi) 	\n\t"\
		" vfmadd132pd	%%zmm10,%%zmm12,%%zmm7 	\n\t	vmovaps		%%zmm6 ,0x180(%%rdi) 	\n\t"\
		"vmovaps		%%zmm5 ,0x140(%%rdi)  	\n\t"\
		"vmovaps		%%zmm1 ,0x040(%%rdi)  	\n\t"\
		"vmovaps		%%zmm3 ,0x0c0(%%rdi)  	\n\t"\
		"vmovaps		%%zmm7 ,0x1c0(%%rdi)  	\n\t"\
		/**/\
	/*...Block 3: __A8-B, t17-24 */\
		"movslq	%[__p8],%%r11	\n\t"\
		"addq	$0x200,%%rdi	\n\t"/* ptr to t17 */\
		"vmovaps		     (%%rax,%%r11,8),%%zmm0 	\n\t	vmovaps		     (%%rbx,%%r11,8),%%zmm8 	\n\t"\
		"vmovaps		     (%%rcx,%%r11,8),%%zmm4 	\n\t	vmovaps		     (%%rdx,%%r11,8),%%zmm9 	\n\t"\
		"vsubpd			%%zmm8 ,%%zmm0 ,%%zmm2 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm8 ,%%zmm0 	\n\t"\
		"vsubpd			%%zmm9 ,%%zmm4 ,%%zmm7 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm9 ,%%zmm4 	\n\t"\
		"vmovaps	0x040(%%rax,%%r11,8),%%zmm1 \n\t	vmovaps	0x040(%%rbx,%%r11,8),%%zmm11\n\t"\
		"vmovaps	0x040(%%rcx,%%r11,8),%%zmm5 \n\t	vmovaps	0x040(%%rdx,%%r11,8),%%zmm12\n\t"\
		"vsubpd			%%zmm11,%%zmm1 ,%%zmm3 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm11,%%zmm1 	\n\t"\
		"vsubpd			%%zmm12,%%zmm5 ,%%zmm6 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm12,%%zmm5 	\n\t"\
		"vmovaps		%%zmm4,%%zmm8			\n\t	vmovaps		%%zmm2,%%zmm9	\n\t"\
		"vfnmadd132pd	%%zmm10,%%zmm0 ,%%zmm4 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm8 ,%%zmm0 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm6 ,%%zmm2 	\n\t"\
		"vfnmadd132pd	%%zmm10,%%zmm9 ,%%zmm6 	\n\t"\
		"vmovaps		%%zmm5,%%zmm11			\n\t	vmovaps		%%zmm3,%%zmm12	\n\t"\
		"vfnmadd132pd	%%zmm10,%%zmm1 ,%%zmm5 	\n\t	vmovaps		%%zmm4 ,0x100(%%rdi) \n\t"\
		" vfmadd231pd	%%zmm10,%%zmm11,%%zmm1 	\n\t	vmovaps		%%zmm0 ,     (%%rdi) \n\t"\
		"vfnmadd231pd	%%zmm10,%%zmm7 ,%%zmm3 	\n\t	vmovaps		%%zmm2 ,0x080(%%rdi) \n\t"\
		" vfmadd132pd	%%zmm10,%%zmm12,%%zmm7 	\n\t	vmovaps		%%zmm6 ,0x180(%%rdi) \n\t"\
		"vmovaps		%%zmm5 ,0x140(%%rdi)  	\n\t"\
		"vmovaps		%%zmm1 ,0x040(%%rdi)  	\n\t"\
		"vmovaps		%%zmm3 ,0x0c0(%%rdi)  	\n\t"\
		"vmovaps		%%zmm7 ,0x1c0(%%rdi)  	\n\t"\
		/**/\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%r14)\n\t"/* ...+p4 */\
	/*...Block 4: __AC-F, t25-32 */\
		"movslq	%[__p12],%%r12	\n\t"\
		"addq	$0x200,%%rdi	\n\t"/* ptr to t25 */\
		"vmovaps		     (%%rax,%%r12,8),%%zmm0 	\n\t	vmovaps		     (%%rbx,%%r12,8),%%zmm8 	\n\t"\
		"vmovaps		     (%%rcx,%%r12,8),%%zmm4 	\n\t	vmovaps		     (%%rdx,%%r12,8),%%zmm9 	\n\t"\
		"vsubpd			%%zmm8 ,%%zmm0 ,%%zmm2 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm8 ,%%zmm0 	\n\t"\
		"vsubpd			%%zmm9 ,%%zmm4 ,%%zmm7 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm9 ,%%zmm4 	\n\t"\
		"vmovaps	0x040(%%rax,%%r12,8),%%zmm1 \n\t	vmovaps	0x040(%%rbx,%%r12,8),%%zmm11\n\t"\
		"vmovaps	0x040(%%rcx,%%r12,8),%%zmm5 \n\t	vmovaps	0x040(%%rdx,%%r12,8),%%zmm12\n\t"\
		"vsubpd			%%zmm11,%%zmm1 ,%%zmm3 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm11,%%zmm1 	\n\t"\
		"vsubpd			%%zmm12,%%zmm5 ,%%zmm6 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm12,%%zmm5 	\n\t"\
		"vmovaps		%%zmm4,%%zmm8			\n\t	vmovaps		%%zmm2,%%zmm9	\n\t"\
		"vfnmadd132pd	%%zmm10,%%zmm0 ,%%zmm4 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm8 ,%%zmm0 	\n\t"\
		" vfmadd231pd	%%zmm10,%%zmm6 ,%%zmm2 	\n\t"\
		"vfnmadd132pd	%%zmm10,%%zmm9 ,%%zmm6 	\n\t"\
		"vmovaps		%%zmm5,%%zmm11			\n\t	vmovaps		%%zmm3,%%zmm12	\n\t"\
		"vfnmadd132pd	%%zmm10,%%zmm1 ,%%zmm5 	\n\t	vmovaps		%%zmm4 ,0x100(%%rdi) \n\t"\
		" vfmadd231pd	%%zmm10,%%zmm11,%%zmm1 	\n\t	vmovaps		%%zmm0 ,     (%%rdi) \n\t"\
		"vfnmadd231pd	%%zmm10,%%zmm7 ,%%zmm3 	\n\t	vmovaps		%%zmm2 ,0x080(%%rdi) \n\t"\
		" vfmadd132pd	%%zmm10,%%zmm12,%%zmm7 	\n\t	vmovaps		%%zmm6 ,0x180(%%rdi) \n\t"\
		"vmovaps		%%zmm5 ,0x140(%%rdi)  	\n\t"\
		"vmovaps		%%zmm1 ,0x040(%%rdi)  	\n\t"\
		"vmovaps		%%zmm3 ,0x0c0(%%rdi)  	\n\t"\
		"vmovaps		%%zmm7 ,0x1c0(%%rdi)  	\n\t"\
		"subq	$0x600,%%rdi	\n\t"/* Revert rdi to point at t1 */\
	/**************************************************************************************/\
	/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
	/**************************************************************************************/\
		/**/\
	/*...Block 1: t1/2,9/10,17/18,25/26 in zmm0-7, resp.; add0 in rax, p4,8,12 in r10,11,12: */\
		"leaq	(%%rax,%%r10,8),%%rbx	\n\t"/* add0+p4 */\
		"leaq	(%%rax,%%r11,8),%%rcx	\n\t"/* add0+p8 */\
		"leaq	(%%rax,%%r12,8),%%rdx	\n\t"/* add0+pC */\
		"vmovaps		     (%%rdi),%%zmm0 	\n\t	vmovaps		0x200(%%rdi),%%zmm2 	\n\t"/*	t1,9 ; */\
		"vmovaps		0x040(%%rdi),%%zmm1 	\n\t	vmovaps		0x240(%%rdi),%%zmm3 	\n\t"/*	t2,10; */\
		"vmovaps		%%zmm0,%%zmm8			\n\t	vmovaps		%%zmm1,%%zmm9	\n\t"/*	rt = t1; it = t2; */\
		" vfmadd231pd	%%zmm10,%%zmm2 ,%%zmm0 	\n\t"/*	 FMA231(1.0,t9 ,t1 ); */\
		"vfnmadd132pd	%%zmm10,%%zmm8 ,%%zmm2 	\n\t"/*	FNMA132(1.0,t9 ,rt ); */\
		" vfmadd231pd	%%zmm10,%%zmm3 ,%%zmm1 	\n\t"/*	 FMA231(1.0,t10,t2 ); */\
		"vfnmadd132pd	%%zmm10,%%zmm9 ,%%zmm3 	\n\t"/*	FNMA132(1.0,t10,it ); */\
		"vmovaps		0x400(%%rdi),%%zmm4 	\n\t	vmovaps		0x600(%%rdi),%%zmm6 	\n\t"/*	t17,25; */\
		"vmovaps		0x440(%%rdi),%%zmm5 	\n\t	vmovaps		0x640(%%rdi),%%zmm7 	\n\t"/*	t18,26; */\
		"vmovaps		%%zmm4,%%zmm11			\n\t	vmovaps		%%zmm5,%%zmm12	\n\t"/*	re = t17; im = t18; */\
		" vfmadd231pd	%%zmm10,%%zmm6 ,%%zmm4 	\n\t"/*	 FMA231(1.0,t25,t17); */\
		"vfnmadd132pd	%%zmm10,%%zmm11,%%zmm6 	\n\t"/*	FNMA132(1.0,t25,re ); */\
		" vfmadd231pd	%%zmm10,%%zmm7 ,%%zmm5 	\n\t"/*	 FMA231(1.0,t26,t18); */\
		"vfnmadd132pd	%%zmm10,%%zmm12,%%zmm7 	\n\t"/*	FNMA132(1.0,t26,im ); */\
		"vmovaps		%%zmm0,%%zmm8			\n\t	vmovaps		%%zmm1,%%zmm9	\n\t"/*	rt = t1; it = t2; */\
		" vfmadd231pd	%%zmm10,%%zmm4 ,%%zmm8 	\n\t"/*	 FMA231(1.0,t17,rt ); */\
		" vfmadd231pd	%%zmm10,%%zmm5 ,%%zmm9 	\n\t"/*	 FMA231(1.0,t18,it ); */\
		"vfnmadd231pd	%%zmm10,%%zmm4 ,%%zmm0 	\n\t"/*	FNMA231(1.0,t17,t1 ); */\
		"vfnmadd231pd	%%zmm10,%%zmm5 ,%%zmm1 	\n\t"/*	FNMA231(1.0,t18,t2 ); */\
		"vmovaps		%%zmm2,%%zmm11			\n\t	vmovaps		%%zmm3,%%zmm12	\n\t"/*	re = t9 ; im = t10; */\
		"vmovaps 0x140(%%rsi),%%zmm13 \n\t vmovaps 0x1c0(%%rsi),%%zmm14 \n\t vmovaps 0x240(%%rsi),%%zmm15 \n\t"/* __t4,8,C; */\
		"vfnmadd231pd	%%zmm10,%%zmm7 ,%%zmm2 	\n\t"/*	FNMA231(1.0,t26,t9 ); */\
		" vfmadd132pd	%%zmm10,%%zmm11,%%zmm7 	\n\t"/*	 FMA132(1.0,t26,re ); */\
		" vfmadd231pd	%%zmm10,%%zmm6 ,%%zmm3 	\n\t"/*	 FMA231(1.0,t25,t10); */\
		"vfnmadd132pd	%%zmm10,%%zmm12,%%zmm6 	\n\t"/*	FNMA132(1.0,t25,im ); */\
		"vmovaps		%%zmm8 ,     (%%rax)  	\n\t"/* __B0r = rt; */\
		"vmovaps		%%zmm9 ,0x040(%%rax)  	\n\t"/* __B0i = it; */\
		"vmovaps %%zmm0,%%zmm8	\n\t	vmovaps %%zmm7,%%zmm11 \n\t vmovaps %%zmm2,%%zmm12	\n\t"/*	rt =t1; re =t26; im =t9; */\
		"vmovaps 0x100(%%rsi),%%zmm4  \n\t vmovaps 0x180(%%rsi),%%zmm5  \n\t vmovaps 0x200(%%rsi),%%zmm9  \n\t"/* __c4,8,C; use unused zmms for these */\
		" vfmadd231pd	%%zmm14,%%zmm1 ,%%zmm0 	\n\t"/*	 FMA231(__t8,t2 ,t1 ); */\
		"vfnmadd231pd	%%zmm14,%%zmm8 ,%%zmm1 	\n\t"/*	FNMA231(__t8,rt ,t2 ); */\
		" vfmadd231pd	%%zmm13,%%zmm6 ,%%zmm7 	\n\t"/*	 FMA231(__t4,t25,t26); */\
		"vfnmadd231pd	%%zmm13,%%zmm11,%%zmm6 	\n\t"/*	FNMA231(__t4,re ,t25); */\
		" vfmadd231pd	%%zmm15,%%zmm3 ,%%zmm2 	\n\t"/*	 FMA231(__tC,t10,t9 ); */\
		"vfnmadd231pd	%%zmm15,%%zmm12,%%zmm3 	\n\t"/*	FNMA231(__tC,im ,t10); */\
		"vmulpd			%%zmm5 ,%%zmm0 ,%%zmm0 	\n\t"/*	t1  *= __c8; */\
		"vmulpd			%%zmm5 ,%%zmm1 ,%%zmm1 	\n\t"/*	t2  *= __c8; */\
		"vmulpd			%%zmm4 ,%%zmm7 ,%%zmm7 	\n\t"/*	t26 *= __c4; */\
		"vmulpd			%%zmm4 ,%%zmm6 ,%%zmm6 	\n\t"/*	t25 *= __c4; */\
		"vmulpd			%%zmm9 ,%%zmm2 ,%%zmm2 	\n\t"/*	t9  *= __cC; */\
		"vmulpd			%%zmm9 ,%%zmm3 ,%%zmm3 	\n\t"/*	t10 *= __cC; */\
		"vmovaps		%%zmm0 ,     (%%rcx)  	\n\t"/* __B8r = t1 ; */\
		"vmovaps		%%zmm1 ,0x040(%%rcx)  	\n\t"/* __B8i = t2 ; */\
		"vmovaps		%%zmm7 ,     (%%rbx)  	\n\t"/* __B4r = t26; */\
		"vmovaps		%%zmm6 ,0x040(%%rbx)  	\n\t"/* __B4i = t25; */\
		"vmovaps		%%zmm2 ,     (%%rdx)  	\n\t"/* __BCr = t9 ; */\
		"vmovaps		%%zmm3 ,0x040(%%rdx)  	\n\t"/* __BCi = t10; */\
		/**/\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%r14,2)\n\t"/* ...+p8 */\
	"leaq	(%%r14,%%r14,2),%%r14	\n\t"	/* p4 + (p4*2) = p12, ptr-offset form */\
	/* Swap Blocks 2,3 here to allow simple-incrementation of array ptrs (decr. not safe w.r.to array padding). */\
		"movslq	%[__p1],%%r10	\n\t"\
	/*...Block 2: t3/4,11/12,19/20,27/28 in zmm0-7: */\
		"addq	$0x080,%%rdi		\n\t"/* t3 */\
		"vmovaps		0x200(%%rdi),%%zmm2 	\n\t	vmovaps		0x240(%%rdi),%%zmm3 	\n\t"/*	t11,12; */\
		"vmovaps		%%zmm2,%%zmm8			\n\t	vmovaps		0x040(%%rsi),%%zmm13	\n\t"/* rt = t11; load __sc */\
		" vfmadd231pd	%%zmm10,%%zmm3 ,%%zmm2 	\n\t"/*	 FMA231(1.0,t12,t11); */\
		"vfnmadd231pd	%%zmm10,%%zmm8 ,%%zmm3 	\n\t"/*	FNMA231(1.0,rt ,t12); */\
		"vmovaps		0x400(%%rdi),%%zmm4 	\n\t	vmovaps		0x440(%%rdi),%%zmm5 	\n\t"/*	t19,20; */\
		"vmovaps		0x600(%%rdi),%%zmm6 	\n\t	vmovaps		0x640(%%rdi),%%zmm7 	\n\t"/*	t27,28; */\
		"vmovaps		%%zmm4,%%zmm11			\n\t	vmovaps		%%zmm6,%%zmm12	\n\t"/*	re = t19; im = t27; */\
		" vfmadd231pd	%%zmm13,%%zmm5 ,%%zmm4 	\n\t"/*	 FMA231(__sc,t20,t19); */\
		"vfnmadd231pd	%%zmm13,%%zmm11,%%zmm5 	\n\t"/*	FNMA231(__sc,re ,t20); */\
		" vfmadd132pd	%%zmm13,%%zmm7 ,%%zmm6 	\n\t"/*	 FMA132(__sc,t27,t28); */\
		" vfmsub132pd	%%zmm13,%%zmm12,%%zmm7 	\n\t"/*	 FMS132(__sc,t28,im ); */\
		"vmovaps		     (%%rdi),%%zmm0 	\n\t	vmovaps		0x040(%%rdi),%%zmm1 	\n\t"/*	t3,4 ; */\
		"vmovaps		-0x040(%%rsi),%%zmm14	\n\t	vmovaps		     (%%rsi),%%zmm15	\n\t"/* load ISRT2, __c */\
		"vmovaps		%%zmm0,%%zmm8			\n\t	vmovaps		%%zmm1,%%zmm9	\n\t"/*	rt = t3; it = t4; */\
		" vfmadd231pd	%%zmm14,%%zmm2 ,%%zmm0 	\n\t"/*	 FMA231(ISRT2,t11,t3 ); */\
		"vfnmadd132pd	%%zmm14,%%zmm8 ,%%zmm2 	\n\t"/*	FNMA132(ISRT2,t11,rt ); */\
		" vfmadd231pd	%%zmm14,%%zmm3 ,%%zmm1 	\n\t"/*	 FMA231(ISRT2,t12,t4 ); */\
		"vfnmadd132pd	%%zmm14,%%zmm9 ,%%zmm3 	\n\t"/*	FNMA132(ISRT2,t12,it ); */\
		"vmulpd			%%zmm15,%%zmm6 ,%%zmm6 	\n\t"/*	t27 *= __c; */\
		"vmulpd			%%zmm15,%%zmm7 ,%%zmm7 	\n\t"/*	t28 *= __c; */\
		"vmovaps		%%zmm6,%%zmm11			\n\t	vmovaps		%%zmm7,%%zmm12	\n\t"/*	re = t27; im = t28; */\
		" vfmsub231pd	%%zmm15,%%zmm4 ,%%zmm6 	\n\t"/*	 FMS231(__c,t19,t27); */\
		" vfmadd132pd	%%zmm15,%%zmm11,%%zmm4 	\n\t"/*	 FMA132(__c,t19,re ); */\
		" vfmsub231pd	%%zmm15,%%zmm5 ,%%zmm7 	\n\t"/*	 FMS231(__c,t20,t28); */\
		" vfmadd132pd	%%zmm15,%%zmm12,%%zmm5 	\n\t"/*	 FMA132(__c,t20,im ); */\
		"vmovaps		%%zmm0,%%zmm8			\n\t	vmovaps		%%zmm1,%%zmm9	\n\t"/*	rt = t3; it = t4; */\
		" vfmsub132pd	%%zmm10,%%zmm4 ,%%zmm0 	\n\t"/*	 FMS132(1.0,t3 ,t19); */\
		" vfmsub132pd	%%zmm10,%%zmm5 ,%%zmm1 	\n\t"/*	 FMS132(1.0,t4 ,t20); */\
		" vfmadd231pd	%%zmm10,%%zmm8 ,%%zmm4 	\n\t"/*	 FMA231(1.0,rt ,t19); */\
		" vfmadd231pd	%%zmm10,%%zmm9 ,%%zmm5 	\n\t"/*	 FMA231(1.0,it ,t20); */\
		"vmovaps		0x4c0(%%rsi),%%zmm14	\n\t	vmovaps		0x5c0(%%rsi),%%zmm15	\n\t"/* load __t1,9 */\
		"vmovaps		%%zmm2,%%zmm11			\n\t	vmovaps		%%zmm3,%%zmm12	\n\t"/*	re = t11; im = t12; */\
		" vfmsub132pd	%%zmm10,%%zmm7 ,%%zmm2 	\n\t"/*	 FMS132(1.0,t11,t28); */\
		" vfmadd132pd	%%zmm10,%%zmm6 ,%%zmm3 	\n\t"/*	 FMA132(1.0,t12,t27); */\
		" vfmadd231pd	%%zmm10,%%zmm11,%%zmm7 	\n\t"/*	 FMA231(1.0,re ,t28); */\
		" vfmsub231pd	%%zmm10,%%zmm12,%%zmm6 	\n\t"/*	 FMS231(1.0,im ,t27); */\
		"vmovaps		%%zmm4,%%zmm8			\n\t	vmovaps		%%zmm0,%%zmm9	\n\t"/*	rt = t19; it = t3; */\
		"vmovaps		0x540(%%rsi),%%zmm10	\n\t	vmovaps		0x640(%%rsi),%%zmm13	\n\t"/* load __t5,D; use 1.0 slot for __t5 */\
		" vfmadd231pd	%%zmm14,%%zmm5 ,%%zmm4 	\n\t"/*	 FMA231(__t1,t20,t19); */\
		"vfnmadd231pd	%%zmm14,%%zmm8 ,%%zmm5 	\n\t"/*	FNMA231(__t1,rt ,t20); */\
		" vfmadd231pd	%%zmm15,%%zmm1 ,%%zmm0 	\n\t"/*	 FMA231(__t9,t4 ,t3 ); */\
		"vfnmadd231pd	%%zmm15,%%zmm9 ,%%zmm1 	\n\t"/*	FNMA231(__t9,it ,t4 ); */\
		"vmovaps		0x480(%%rsi),%%zmm14	\n\t	vmovaps		0x580(%%rsi),%%zmm15	\n\t"/* load __c1,9 */\
		"vmovaps		%%zmm7,%%zmm11			\n\t	vmovaps		%%zmm2,%%zmm12	\n\t"/*	re = t28; im = t11; */\
		" vfmadd231pd	%%zmm10,%%zmm6 ,%%zmm7 	\n\t"/*	 FMA231(__t5,t27,t28); */\
		"vfnmadd231pd	%%zmm10,%%zmm11,%%zmm6 	\n\t"/*	FNMA231(__t5,re ,t27); */\
		" vfmadd231pd	%%zmm13,%%zmm3 ,%%zmm2 	\n\t"/*	 FMA231(__tD,t12,t11); */\
		"vfnmadd231pd	%%zmm13,%%zmm12,%%zmm3 	\n\t"/*	FNMA231(__tD,im ,t12); */\
		"vmovaps		0x500(%%rsi),%%zmm10	\n\t	vmovaps		0x600(%%rsi),%%zmm13	\n\t"/* load __c5,D; use 1.0 slot for __t5 */\
		"vmulpd			%%zmm14,%%zmm4 ,%%zmm4 	\n\t"/*	t19 *= __c1; */\
		"vmulpd			%%zmm14,%%zmm5 ,%%zmm5 	\n\t"/*	t20 *= __c1; */\
		"vmulpd			%%zmm15,%%zmm0 ,%%zmm0 	\n\t"/*	t3  *= __c9; */\
		"vmulpd			%%zmm15,%%zmm1 ,%%zmm1 	\n\t"/*	t4  *= __c9; */\
		"vmulpd			%%zmm10,%%zmm7 ,%%zmm7 	\n\t"/*	t28 *= __c5; */\
		"vmulpd			%%zmm10,%%zmm6 ,%%zmm6 	\n\t"/*	t27 *= __c5; */\
		"vmulpd			%%zmm13,%%zmm2 ,%%zmm2 	\n\t"/*	t11 *= __cD; */\
		"vmulpd			%%zmm13,%%zmm3 ,%%zmm3 	\n\t"/*	t12 *= __cD; */\
		"vmovaps	%%zmm4 ,     (%%rax,%%r10,8)  	\n\t"/* __B1r = t19; */\
		"vmovaps	%%zmm5 ,0x040(%%rax,%%r10,8)  	\n\t"/* __B1i = t20; */\
		"vmovaps	%%zmm0 ,     (%%rcx,%%r10,8)  	\n\t"/* __B9r = t3 ; */\
		"vmovaps	%%zmm1 ,0x040(%%rcx,%%r10,8)  	\n\t"/* __B9i = t4 ;	Since cannot assume p2 = (p1+p1) due to array padding, revert:*/\
		"vmovaps	%%zmm7 ,     (%%rbx,%%r10,8)  	\n\t"/* __B5r = t28 */\
		"vmovaps	%%zmm6 ,0x040(%%rbx,%%r10,8)  	\n\t"/* __B5i = t27 */\
		"vmovaps	%%zmm2 ,     (%%rdx,%%r10,8)  	\n\t"/* __BDr = t11 */\
		"vmovaps	%%zmm3 ,0x040(%%rdx,%%r10,8)  	\n\t"/* __BDi = t12 */\
		/**/\
	/*...Block 3: t5/6,13/14,21/22,29/30 in zmm0-7: */\
		"vmovaps 0x880(%%rsi),%%zmm10	\n\t"/* retore 1.0 to zmm10 */\
		"addq	$0x080,%%rdi		\n\t"/* t5 */\
		"movslq	%[__p2],%%r11	\n\t"\
		"vmovaps		     (%%rdi),%%zmm0 	\n\t	vmovaps		0x200(%%rdi),%%zmm2 	\n\t"/*	t5,13; */\
		"vmovaps		0x040(%%rdi),%%zmm1 	\n\t	vmovaps		0x240(%%rdi),%%zmm3 	\n\t"/*	t6,14; */\
		"vmovaps		%%zmm0,%%zmm8			\n\t	vmovaps		%%zmm1,%%zmm9	\n\t"/*	rt = t5; it = t6; */\
		" vfmadd231pd	%%zmm10,%%zmm3 ,%%zmm0 	\n\t"/*	 FMA231(1.0,t14,t5 ); */\
		"vfnmadd132pd	%%zmm10,%%zmm8 ,%%zmm3 	\n\t"/*	FNMA132(1.0,t14,rt ); */\
		"vfnmadd231pd	%%zmm10,%%zmm2 ,%%zmm1 	\n\t"/*	FNMA231(1.0,t13,t6 ); */\
		" vfmadd132pd	%%zmm10,%%zmm9 ,%%zmm2 	\n\t"/*	 FMA132(1.0,t13,it ); */\
		"vmovaps		0x400(%%rdi),%%zmm4 	\n\t	vmovaps		0x440(%%rdi),%%zmm5 	\n\t"/*	t21,22; */\
		"vmovaps		0x600(%%rdi),%%zmm6 	\n\t	vmovaps		0x640(%%rdi),%%zmm7 	\n\t"/*	t29,30; */\
		"vmovaps		%%zmm5,%%zmm11			\n\t	vmovaps		%%zmm7,%%zmm12	\n\t"/*	re = t22; im = t30; */\
		"vfnmadd231pd	%%zmm10,%%zmm4 ,%%zmm5 	\n\t"/*	FNMA231(1.0,t21,t22); */\
		" vfmadd132pd	%%zmm10,%%zmm11,%%zmm4 	\n\t"/*	 FMA132(1.0,t21,re ); */\
		" vfmadd231pd	%%zmm10,%%zmm6 ,%%zmm7 	\n\t"/*	 FMA231(1.0,t29,t30); */\
		" vfmsub132pd	%%zmm10,%%zmm12,%%zmm6 	\n\t"/*	 FMS132(1.0,t29,im ); */\
		"vmovaps		%%zmm4,%%zmm8			\n\t	vmovaps		%%zmm5,%%zmm9	\n\t"/*	rt = t21; it = t22; */\
		"vfnmadd231pd	%%zmm10,%%zmm6 ,%%zmm4 	\n\t"/*	FNMA231(1.0,t29,t21); */\
		" vfmadd132pd	%%zmm10,%%zmm8 ,%%zmm6 	\n\t"/*	 FMA132(1.0,t29,rt ); */\
		"vfnmadd231pd	%%zmm10,%%zmm7 ,%%zmm5 	\n\t"/*	FNMA231(1.0,t30,t22); */\
		" vfmadd132pd	%%zmm10,%%zmm9 ,%%zmm7 	\n\t"/*	 FMA132(1.0,t30,it ); */\
		"vmovaps		-0x040(%%rsi),%%zmm13	\n\t"/* load ISRT2 */\
		"vmovaps		%%zmm0,%%zmm8			\n\t	vmovaps		%%zmm1,%%zmm9	\n\t"/*	rt = t5; it = t6; */\
		"vfnmadd231pd	%%zmm13,%%zmm4 ,%%zmm0 	\n\t"/*	FNMA231(ISRT2,t21,t5 ); */\
		" vfmadd132pd	%%zmm13,%%zmm8 ,%%zmm4 	\n\t"/*	 FMA132(ISRT2,t21,rt ); */\
		"vfnmadd231pd	%%zmm13,%%zmm5 ,%%zmm1 	\n\t"/*	FNMA231(ISRT2,t22,t6 ); */\
		" vfmadd132pd	%%zmm13,%%zmm9 ,%%zmm5 	\n\t"/*	 FMA132(ISRT2,t22,it ); */\
		"vmovaps		0x2c0(%%rsi),%%zmm14	\n\t	vmovaps		0x3c0(%%rsi),%%zmm15	\n\t"/* load __t2,A */\
		"vmovaps		%%zmm2,%%zmm11			\n\t	vmovaps		%%zmm3,%%zmm12	\n\t"/*	re =t13; im =t14; */\
		" vfmadd231pd	%%zmm13,%%zmm6 ,%%zmm2 	\n\t"/*	 FMA231(ISRT2,t29,t13); */\
		"vfnmadd132pd	%%zmm13,%%zmm11,%%zmm6 	\n\t"/*	FNMA132(ISRT2,t29,re ); */\
		"vfnmadd231pd	%%zmm13,%%zmm7 ,%%zmm3 	\n\t"/*	FNMA231(ISRT2,t30,t14); */\
		" vfmadd132pd	%%zmm13,%%zmm12,%%zmm7 	\n\t"/*	 FMA132(ISRT2,t30,im ); */\
		"vmovaps		%%zmm4,%%zmm8			\n\t	vmovaps		%%zmm0,%%zmm9	\n\t"/*	rt =t21; it = t5; */\
		"vmovaps		0x340(%%rsi),%%zmm10	\n\t	vmovaps		0x440(%%rsi),%%zmm13	\n\t"/* load __t6,E; use 1.0 slot for __t5 */\
		" vfmadd231pd	%%zmm14,%%zmm5 ,%%zmm4 	\n\t"/*	 FMA231(__t2,t22,t21); */\
		"vfnmadd231pd	%%zmm14,%%zmm8 ,%%zmm5 	\n\t"/*	FNMA231(__t2,rt ,t22); */\
		" vfmadd231pd	%%zmm15,%%zmm1 ,%%zmm0 	\n\t"/*	 FMA231(__tA,t6 ,t5 ); */\
		"vfnmadd231pd	%%zmm15,%%zmm9 ,%%zmm1 	\n\t"/*	FNMA231(__tA,it ,t6 ); */\
		"vmovaps		0x280(%%rsi),%%zmm14	\n\t	vmovaps		0x380(%%rsi),%%zmm15	\n\t"/* load __c2,A */\
		"vmovaps		%%zmm7,%%zmm11			\n\t	vmovaps		%%zmm3,%%zmm12	\n\t"/*	re =t30; im =t14; */\
		" vfmadd231pd	%%zmm10,%%zmm6 ,%%zmm7 	\n\t"/*	 FMA231(__t6,t29,t30); */\
		"vfnmadd231pd	%%zmm10,%%zmm11,%%zmm6 	\n\t"/*	FNMA231(__t6,re ,t29); */\
		" vfmadd231pd	%%zmm13,%%zmm2 ,%%zmm3 	\n\t"/*	 FMA231(__tE,t13,t14); */\
		"vfnmadd231pd	%%zmm13,%%zmm12,%%zmm2 	\n\t"/*	FNMA231(__tE,im ,t13); */\
		"vmovaps		0x300(%%rsi),%%zmm10	\n\t	vmovaps		0x400(%%rsi),%%zmm13	\n\t"/* load __c6,E; use 1.0 slot for __t5 */\
		"vmulpd			%%zmm14,%%zmm4 ,%%zmm4 	\n\t"/*	t21 *= __c2; */\
		"vmulpd			%%zmm14,%%zmm5 ,%%zmm5 	\n\t"/*	t22 *= __c2; */\
		"vmulpd			%%zmm15,%%zmm0 ,%%zmm0 	\n\t"/*	t5  *= __cA; */\
		"vmulpd			%%zmm15,%%zmm1 ,%%zmm1 	\n\t"/*	t6  *= __cA; */\
		"vmulpd			%%zmm10,%%zmm7 ,%%zmm7 	\n\t"/*	t30 *= __c6; */\
		"vmulpd			%%zmm10,%%zmm6 ,%%zmm6 	\n\t"/*	t29 *= __c6; */\
		"vmulpd			%%zmm13,%%zmm3 ,%%zmm3 	\n\t"/*	t14 *= __cE; */\
		"vmulpd			%%zmm13,%%zmm2 ,%%zmm2 	\n\t"/*	t13 *= __cE; */\
		"vmovaps	%%zmm4 ,     (%%rax,%%r11,8)  	\n\t"/* __B2r = t21; */\
		"vmovaps	%%zmm5 ,0x040(%%rax,%%r11,8)  	\n\t"/* __B2i = t22; */\
		"vmovaps	%%zmm0 ,     (%%rcx,%%r11,8)  	\n\t"/* __BAr = t5 ; */\
		"vmovaps	%%zmm1 ,0x040(%%rcx,%%r11,8)  	\n\t"/* __BAi = t6 ; */\
		"vmovaps	%%zmm7 ,     (%%rbx,%%r11,8)  	\n\t"/* __B6r = t30; */\
		"vmovaps	%%zmm6 ,0x040(%%rbx,%%r11,8)  	\n\t"/* __B6i = t29; */\
		"vmovaps	%%zmm3 ,     (%%rdx,%%r11,8)  	\n\t"/* __BEr = t14; */\
		"vmovaps	%%zmm2 ,0x040(%%rdx,%%r11,8)  	\n\t"/* __BEi = t13; */\
		/**/\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%r14)\n\t"/* ...+p12 */\
	/*...Block 4: t7/8,15/16,23/24,31/32 in zmm0-7: */\
		"movslq	%[__p3],%%r12	\n\t"\
		"vmovaps 0x880(%%rsi),%%zmm10	\n\t"/* retore 1.0 to zmm10 */\
		"addq	$0x080,%%rdi		\n\t"/* t7 */\
		"vmovaps		0x200(%%rdi),%%zmm2 	\n\t	vmovaps		0x240(%%rdi),%%zmm3 	\n\t"/*	t15,16; */\
		"vmovaps		%%zmm3,%%zmm8			\n\t	vmovaps		0x040(%%rsi),%%zmm13	\n\t"/* rt = t16; load __sc */\
		" vfmadd231pd	%%zmm10,%%zmm2 ,%%zmm3 	\n\t"/*	 FMA231(1.0,t15,t16); */\
		" vfmsub132pd	%%zmm10,%%zmm8 ,%%zmm2 	\n\t"/*	 FMS132(1.0,t15,rt ); */\
		"vmovaps		0x400(%%rdi),%%zmm4 	\n\t	vmovaps		0x440(%%rdi),%%zmm5 	\n\t"/*	t23,24; */\
		"vmovaps		0x600(%%rdi),%%zmm6 	\n\t	vmovaps		0x640(%%rdi),%%zmm7 	\n\t"/*	t31,32; */\
		"vmovaps		%%zmm4,%%zmm11			\n\t	vmovaps		%%zmm6,%%zmm12	\n\t"/*	re = t23; im = t31; */\
		" vfmadd231pd	%%zmm13,%%zmm7 ,%%zmm6 	\n\t"/*	 FMA231(__sc,t32,t31); */\
		"vfnmadd231pd	%%zmm13,%%zmm12,%%zmm7 	\n\t"/*	FNMA231(__sc,im ,t32); */\
		" vfmadd132pd	%%zmm13,%%zmm5 ,%%zmm4 	\n\t"/*	 FMA132(__sc,t23,t24); */\
		" vfmsub132pd	%%zmm13,%%zmm11,%%zmm5 	\n\t"/*	 FMS132(__sc,t24,re ); */\
		"vmovaps		-0x040(%%rsi),%%zmm14	\n\t	vmovaps		     (%%rsi),%%zmm15	\n\t"/* load ISRT2, __c */\
		"vmovaps		     (%%rdi),%%zmm0 	\n\t	vmovaps		0x040(%%rdi),%%zmm1 	\n\t"/*	t7,8 ; */\
		"vmovaps		%%zmm0,%%zmm8			\n\t	vmovaps		%%zmm1,%%zmm9	\n\t"/*	rt = t7; it = t8; */\
		"vfnmadd231pd	%%zmm14,%%zmm2 ,%%zmm0 	\n\t"/*	FNMA231(ISRT2,t15,t7 ); */\
		" vfmadd132pd	%%zmm14,%%zmm8 ,%%zmm2 	\n\t"/*	 FMA132(ISRT2,t15,rt ); */\
		"vfnmadd231pd	%%zmm14,%%zmm3 ,%%zmm1 	\n\t"/*	FNMA231(ISRT2,t16,t8 ); */\
		" vfmadd132pd	%%zmm14,%%zmm9 ,%%zmm3 	\n\t"/*	 FMA132(ISRT2,t16,it ); */\
		"vmulpd			%%zmm15,%%zmm6 ,%%zmm6 	\n\t"/*	t31 *= __c; */\
		"vmulpd			%%zmm15,%%zmm7 ,%%zmm7 	\n\t"/*	t32 *= __c; */\
		"vmovaps		%%zmm6,%%zmm11			\n\t	vmovaps		%%zmm7,%%zmm12	\n\t"/*	re = t31; im = t32; */\
		" vfmadd231pd	%%zmm15,%%zmm4 ,%%zmm6 	\n\t"/*	 FMA231(__c,t23,t31); */\
		" vfmsub132pd	%%zmm15,%%zmm11,%%zmm4 	\n\t"/*	 FMS132(__c,t23,re ); */\
		" vfmadd231pd	%%zmm15,%%zmm5 ,%%zmm7 	\n\t"/*	 FMA231(__c,t24,t32); */\
		" vfmsub132pd	%%zmm15,%%zmm12,%%zmm5 	\n\t"/*	 FMS132(__c,t24,im ); */\
		"vmovaps		%%zmm2,%%zmm8			\n\t	vmovaps		%%zmm3,%%zmm9	\n\t"/*	rt = t15; it = t16; */\
		" vfmsub132pd	%%zmm10,%%zmm7 ,%%zmm2 	\n\t"/*	 FMS132(1.0,t15,t32); */\
		" vfmadd132pd	%%zmm10,%%zmm6 ,%%zmm3 	\n\t"/*	 FMA132(1.0,t16,t31); */\
		" vfmadd231pd	%%zmm10,%%zmm8 ,%%zmm7 	\n\t"/*	 FMA231(1.0,rt ,t32); */\
		" vfmsub231pd	%%zmm10,%%zmm9 ,%%zmm6 	\n\t"/*	 FMS231(1.0,it ,t31); */\
		"vmovaps		0x6c0(%%rsi),%%zmm14	\n\t	vmovaps		0x7c0(%%rsi),%%zmm15	\n\t"/* load __t3,B */\
		"vmovaps		%%zmm0,%%zmm11			\n\t	vmovaps		%%zmm1,%%zmm12	\n\t"/*	re = t7; im = t8; */\
		" vfmsub132pd	%%zmm10,%%zmm4 ,%%zmm0 	\n\t"/*	 FMS132(1.0,t7 ,t23 ); */\
		" vfmsub132pd	%%zmm10,%%zmm5 ,%%zmm1 	\n\t"/*	 FMS132(1.0,t8 ,t24 ); */\
		" vfmadd231pd	%%zmm10,%%zmm11,%%zmm4 	\n\t"/*	 FMA231(1.0,re ,t23 ); */\
		" vfmadd231pd	%%zmm10,%%zmm12,%%zmm5 	\n\t"/*	 FMA231(1.0,im ,t24 ); */\
		"vmovaps		%%zmm4,%%zmm8			\n\t	vmovaps		%%zmm0,%%zmm9	\n\t"/*	rt = t23; it = t7; */\
		"vmovaps		0x740(%%rsi),%%zmm10	\n\t	vmovaps		0x840(%%rsi),%%zmm13	\n\t"/* load __t7,F; use 1.0 slot for __t7 */\
		" vfmadd231pd	%%zmm14,%%zmm5 ,%%zmm4 	\n\t"/*	 FMA231(__t3,t24,t23); */\
		"vfnmadd231pd	%%zmm14,%%zmm8 ,%%zmm5 	\n\t"/*	FNMA231(__t3,rt ,t24); */\
		" vfmadd231pd	%%zmm15,%%zmm1 ,%%zmm0 	\n\t"/*	 FMA231(__tB,t8 ,t7 ); */\
		"vfnmadd231pd	%%zmm15,%%zmm9 ,%%zmm1 	\n\t"/*	FNMA231(__tB,it ,t8 ); */\
		"vmovaps		0x680(%%rsi),%%zmm14	\n\t	vmovaps		0x780(%%rsi),%%zmm15	\n\t"/* load __c3,B */\
		"vmovaps		%%zmm7,%%zmm11			\n\t	vmovaps		%%zmm2,%%zmm12	\n\t"/*	re = t32; im = t15; */\
		" vfmadd231pd	%%zmm10,%%zmm6 ,%%zmm7 	\n\t"/*	 FMA231(__t7,t31,t32); */\
		"vfnmadd231pd	%%zmm10,%%zmm11,%%zmm6 	\n\t"/*	FNMA231(__t7,re ,t31); */\
		" vfmadd231pd	%%zmm13,%%zmm3 ,%%zmm2 	\n\t"/*	 FMA231(__tF,t16,t15); */\
		"vfnmadd231pd	%%zmm13,%%zmm12,%%zmm3 	\n\t"/*	FNMA231(__tF,im ,t16); */\
		"vmovaps		0x700(%%rsi),%%zmm10	\n\t	vmovaps		0x800(%%rsi),%%zmm13	\n\t"/* load __c7,F; use 1.0 slot for __t5 */\
		"vmulpd			%%zmm14,%%zmm4 ,%%zmm4 	\n\t"/*	t23 *= __c3; */\
		"vmulpd			%%zmm14,%%zmm5 ,%%zmm5 	\n\t"/*	t24 *= __c3; */\
		"vmulpd			%%zmm15,%%zmm0 ,%%zmm0 	\n\t"/*	t7  *= __cB; */\
		"vmulpd			%%zmm15,%%zmm1 ,%%zmm1 	\n\t"/*	t8  *= __cB; */\
		"vmulpd			%%zmm10,%%zmm7 ,%%zmm7 	\n\t"/*	t32 *= __c7; */\
		"vmulpd			%%zmm10,%%zmm6 ,%%zmm6 	\n\t"/*	t31 *= __c7; */\
		"vmulpd			%%zmm13,%%zmm2 ,%%zmm2 	\n\t"/*	t15 *= __cF; */\
		"vmulpd			%%zmm13,%%zmm3 ,%%zmm3 	\n\t"/*	t16 *= __cF; */\
		"vmovaps	%%zmm4 ,     (%%rax,%%r12,8)  	\n\t"/* __B3r = t23; */\
		"vmovaps	%%zmm5 ,0x040(%%rax,%%r12,8)  	\n\t"/* __B3i = t24; */\
		"vmovaps	%%zmm0 ,     (%%rcx,%%r12,8)  	\n\t"/* __BBr = t7 ; */\
		"vmovaps	%%zmm1 ,0x040(%%rcx,%%r12,8)  	\n\t"/* __BBi = t8 ; */\
		"vmovaps	%%zmm7 ,     (%%rbx,%%r12,8)  	\n\t"/* __B7r = t32; */\
		"vmovaps	%%zmm6 ,0x040(%%rbx,%%r12,8)  	\n\t"/* __B7i = t31; */\
		"vmovaps	%%zmm2 ,     (%%rdx,%%r12,8)  	\n\t"/* __BFr = t15; */\
		"vmovaps	%%zmm3 ,0x040(%%rdx,%%r12,8)  	\n\t"/* __BFi = t16; */\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r1] "m" (Xr1)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r9","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	// The non-name-suffixed DIF,DIT macros in AVX2 mode are based on the AVX-mode _V2 macros:
	// 64 ADD, 112 FMA, 48 pure-MUL, 216 MOVAPS
	#define SSE2_RADIX16_DIT_TWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr0,Xtwo,Xcc0,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/*
		// Pass 1:
		j = p4*8;	// Set j to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k = 0x100;	// Set k to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		// 0x0c0 = literal-bytewise address offset between the pointers to the first and second DFT's twiddles:
		i0 = (vec_dbl *)a; i1 = (vec_dbl *)(a+p1); i2 = (vec_dbl *)(a+p2); i3 = (vec_dbl *)(a+p3);
		o0 = r00; o1 = r00+2; o2 = r00+4; o3 = r00+6;
		c_tmp = cc0+6;	// c2,1,3
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,0x0c0, o0,o1,o2,o3,k)
	*/\
		"movq	%[__two],%%rsi 				\n\t	movslq	%[__p4],%%rdi		\n\t"\
		"movq	%[__add0],%%rax				\n\t	movq	%[__r0],%%r10		\n\t"\
		"movslq	%[__p1],%%rbx				\n\t	leaq	(%%rax,%%rbx,8),%%rbx	\n\t"\
		"movslq	%[__p2],%%rcx				\n\t	leaq	(%%rax,%%rcx,8),%%rcx	\n\t"\
		"movslq	%[__p3],%%rdx				\n\t	leaq	(%%rax,%%rdx,8),%%rdx	\n\t"\
/* two: */"vmovaps (%%rsi),%%zmm15 \n\t shlq $3,%%rdi \n\t	leaq 0x040(%%rdx,%%rdi),%%r11 \n\t"\
		"vmovaps	     (%%rax),%%zmm0		\n\t	vmovaps	     (%%rax,%%rdi),%%zmm8 	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1		\n\t	vmovaps	0x040(%%rax,%%rdi),%%zmm9 	\n\t"\
		"vmovaps	     (%%rbx),%%zmm2		\n\t	vmovaps	     (%%rbx,%%rdi),%%zmm10	\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm3		\n\t	vmovaps	0x040(%%rbx,%%rdi),%%zmm11	\n\t"\
		"vmovaps	     (%%rcx),%%zmm4		\n\t	vmovaps	     (%%rcx,%%rdi),%%zmm12	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5		\n\t	vmovaps	0x040(%%rcx,%%rdi),%%zmm13	\n\t"\
		"vmovaps	     (%%rdx),%%zmm6		\n\t	vmovaps	     (%%rdx,%%rdi),%%zmm14	\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm7		\n\t"/*	vmovaps	0x040(%%rdx,%%rdi),%%zmm15	\n\t"*/\
		"vsubpd	%%zmm2 ,%%zmm0,%%zmm0		\n\t	vsubpd	%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd	%%zmm3 ,%%zmm1,%%zmm1		\n\t	vsubpd	%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd	%%zmm6 ,%%zmm4,%%zmm4		\n\t	vsubpd	%%zmm14,%%zmm12,%%zmm12		\n\t"\
		"vsubpd	%%zmm7 ,%%zmm5,%%zmm5		\n\t	vsubpd	(%%r11),%%zmm13,%%zmm13		\n\t"\
	"vfmadd132pd %%zmm15,%%zmm0,%%zmm2		\n\t	vfmadd132pd	%%zmm15,%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd132pd %%zmm15,%%zmm1,%%zmm3		\n\t	vfmadd132pd	%%zmm15,%%zmm9 ,%%zmm11		\n\t"\
	"vfmadd132pd %%zmm15,%%zmm4,%%zmm6		\n\t	vfmadd132pd	%%zmm15,%%zmm12,%%zmm14		\n\t"\
	"vfmadd132pd %%zmm15,%%zmm5,%%zmm7		\n\t	vfmadd132pd	(%%r11),%%zmm13,%%zmm15		\n\t"\
		"vmovaps	%%zmm12,    (%%r10)		\n\t	vmovaps (%%rsi),%%zmm12			\n\t"/* spill zmm12 to free up reg for 2.0 */\
		"vsubpd	%%zmm6 ,%%zmm2,%%zmm2		\n\t	vsubpd	%%zmm14,%%zmm10,%%zmm10		\n\t"\
		"vsubpd	%%zmm7 ,%%zmm3,%%zmm3		\n\t	vsubpd	%%zmm15,%%zmm11,%%zmm11		\n\t"\
		"vsubpd	%%zmm5 ,%%zmm0,%%zmm0		\n\t	vsubpd	%%zmm13,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd	%%zmm4 ,%%zmm1,%%zmm1		\n\t	vsubpd	(%%r10),%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd %%zmm12,%%zmm2,%%zmm6		\n\t	vfmadd132pd	%%zmm12,%%zmm10,%%zmm14		\n\t"\
	"vfmadd132pd %%zmm12,%%zmm3,%%zmm7		\n\t	vfmadd132pd	%%zmm12,%%zmm11,%%zmm15		\n\t"\
	"vfmadd132pd %%zmm12,%%zmm0,%%zmm5		\n\t	vfmadd132pd	%%zmm12,%%zmm8 ,%%zmm13		\n\t"\
	"vfmadd132pd %%zmm12,%%zmm1,%%zmm4		\n\t	vfmadd132pd	(%%r10),%%zmm9 ,%%zmm12		\n\t"\
		"vmovaps	%%zmm6,    (%%r10)		\n\t	vmovaps	%%zmm7,0x040(%%r10)		\n\t"/* Write B0 to free up 2 regs */\
		/* B2 = t0*~w2 = t0*[c2-I.s1] */\
		"movq	%[__cc0],%%rsi				\n\t	vmovaps %%zmm14,0x200(%%r10)			\n\t"\
/* cc0+6: */"addq	$0x180,%%rsi			\n\t	vmovaps %%zmm15,0x240(%%r10)			\n\t"\
		"vmovaps	%%zmm2,%%zmm6			\n\t	vmovaps	%%zmm10,%%zmm14					\n\t"\
		"vmovaps	%%zmm3,%%zmm7			\n\t	vmovaps	%%zmm11,%%zmm15					\n\t"\
		"vmulpd	     (%%rsi),%%zmm6,%%zmm6	\n\t	vmulpd	0x180(%%rsi),%%zmm14,%%zmm14	\n\t"\
		"vmulpd	     (%%rsi),%%zmm7,%%zmm7	\n\t	vmulpd	0x180(%%rsi),%%zmm15,%%zmm15	\n\t"\
	" vfmadd231pd 0x40(%%rsi),%%zmm3,%%zmm6	\n\t  vfmadd231pd 0x1c0(%%rsi),%%zmm11,%%zmm14	\n\t"\
	"vfnmadd231pd 0x40(%%rsi),%%zmm2,%%zmm7	\n\t vfnmadd231pd 0x1c0(%%rsi),%%zmm10,%%zmm15	\n\t"\
		"vmovaps	%%zmm6,0x100(%%r10)		\n\t	vmovaps %%zmm14,0x300(%%r10)			\n\t"\
		"vmovaps	%%zmm7,0x140(%%r10)		\n\t	vmovaps %%zmm15,0x340(%%r10)			\n\t"\
		/* B1 = t1*~w1 = t1*[c1-I.s1] */\
		"vmovaps	0x080(%%rsi),%%zmm2		\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm3		\n\t"\
		"vmovaps	%%zmm5,%%zmm6			\n\t	vmovaps	%%zmm13,%%zmm14					\n\t"\
		"vmovaps	%%zmm1,%%zmm7			\n\t	vmovaps	%%zmm9 ,%%zmm15					\n\t"\
		"vmulpd	%%zmm2,%%zmm6,%%zmm6		\n\t	vmulpd	0x200(%%rsi),%%zmm14,%%zmm14	\n\t"\
		"vmulpd	%%zmm2,%%zmm7,%%zmm7		\n\t	vmulpd	0x200(%%rsi),%%zmm15,%%zmm15	\n\t"\
	" vfmadd231pd	%%zmm3,%%zmm1,%%zmm6	\n\t  vfmadd231pd 0x240(%%rsi),%%zmm9 ,%%zmm14	\n\t"\
	"vfnmadd231pd	%%zmm3,%%zmm5,%%zmm7	\n\t vfnmadd231pd 0x240(%%rsi),%%zmm13,%%zmm15	\n\t"\
		"vmovaps	0x100(%%rsi),%%zmm2		\n\t"\
		"vmovaps	0x140(%%rsi),%%zmm3		\n\t"\
		"vmovaps	%%zmm6,0x080(%%r10)		\n\t	vmovaps %%zmm14,0x280(%%r10)			\n\t"\
		"vmovaps	%%zmm7,0x0c0(%%r10)		\n\t	vmovaps %%zmm15,0x2c0(%%r10)			\n\t"\
		/* B3 = t3*~w3 = t3*[c3-I.s3] */\
		"vmovaps	%%zmm0,%%zmm6			\n\t	vmovaps	%%zmm8 ,%%zmm14					\n\t"\
		"vmovaps	%%zmm4,%%zmm7			\n\t	vmovaps	%%zmm12,%%zmm15					\n\t"\
		"vmulpd	%%zmm2,%%zmm6,%%zmm6		\n\t	vmulpd	0x280(%%rsi),%%zmm14,%%zmm14	\n\t"\
		"vmulpd	%%zmm2,%%zmm7,%%zmm7		\n\t	vmulpd	0x280(%%rsi),%%zmm15,%%zmm15	\n\t"\
	" vfmadd231pd	%%zmm3,%%zmm4,%%zmm6	\n\t  vfmadd231pd 0x2c0(%%rsi),%%zmm12,%%zmm14	\n\t"\
	"vfnmadd231pd	%%zmm3,%%zmm0,%%zmm7	\n\t vfnmadd231pd 0x2c0(%%rsi),%%zmm8 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm6,0x180(%%r10)		\n\t	vmovaps	%%zmm14,0x380(%%r10)			\n\t"\
		"vmovaps	%%zmm7,0x1c0(%%r10)		\n\t	vmovaps	%%zmm15,0x3c0(%%r10)			\n\t"\
	/*
		i0 = (vec_dbl *)(a+p8); i1 = (vec_dbl *)(a+p9); i2 = (vec_dbl *)(a+pA); i3 = (vec_dbl *)(a+pB);
		o0 += 16; o1 += 16; o2 += 16; o3 += 16;
		c_tmp += 12;		// cA,9,B
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,0x180, o0,o1,o2,o3,k)
	*/\
		"movq	%[__two],%%rsi 		\n\t	movslq	%[__p8],%%r11		\n\t"\
		"leaq (%%rax,%%r11,8),%%rax	\n\t	addq	 $0x400,%%r10		\n\t"/* r0 + 16 */\
		"leaq (%%rbx,%%r11,8),%%rbx	\n\t"\
		"leaq (%%rcx,%%r11,8),%%rcx	\n\t"\
		"leaq (%%rdx,%%r11,8),%%rdx	\n\t"\
	"vmovaps	(%%rsi),%%zmm15		\n\t"/* two */"	leaq	0x040(%%rdx,%%rdi),%%r11	\n\t"\
		"vmovaps	     (%%rax),%%zmm0		\n\t	vmovaps	     (%%rax,%%rdi),%%zmm8 	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1		\n\t	vmovaps	0x040(%%rax,%%rdi),%%zmm9 	\n\t"\
		"vmovaps	     (%%rbx),%%zmm2		\n\t	vmovaps	     (%%rbx,%%rdi),%%zmm10	\n\t"\
		"vmovaps	0x040(%%rbx),%%zmm3		\n\t	vmovaps	0x040(%%rbx,%%rdi),%%zmm11	\n\t"\
		"vmovaps	     (%%rcx),%%zmm4		\n\t	vmovaps	     (%%rcx,%%rdi),%%zmm12	\n\t"\
		"vmovaps	0x040(%%rcx),%%zmm5		\n\t	vmovaps	0x040(%%rcx,%%rdi),%%zmm13	\n\t"\
		"vmovaps	     (%%rdx),%%zmm6		\n\t	vmovaps	     (%%rdx,%%rdi),%%zmm14	\n\t"\
		"vmovaps	0x040(%%rdx),%%zmm7		\n\t"/*	vmovaps	0x040(%%rdx,%%rdi),%%zmm15	\n\t"*/\
		"vsubpd	%%zmm2 ,%%zmm0,%%zmm0		\n\t	vsubpd	%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd	%%zmm3 ,%%zmm1,%%zmm1		\n\t	vsubpd	%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd	%%zmm6 ,%%zmm4,%%zmm4		\n\t	vsubpd	%%zmm14,%%zmm12,%%zmm12		\n\t"\
		"vsubpd	%%zmm7 ,%%zmm5,%%zmm5		\n\t	vsubpd	(%%r11),%%zmm13,%%zmm13		\n\t"\
	"vfmadd132pd %%zmm15,%%zmm0,%%zmm2		\n\t	vfmadd132pd	%%zmm15,%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd132pd %%zmm15,%%zmm1,%%zmm3		\n\t	vfmadd132pd	%%zmm15,%%zmm9 ,%%zmm11		\n\t"\
	"vfmadd132pd %%zmm15,%%zmm4,%%zmm6		\n\t	vfmadd132pd	%%zmm15,%%zmm12,%%zmm14		\n\t"\
	"vfmadd132pd %%zmm15,%%zmm5,%%zmm7		\n\t	vfmadd132pd	(%%r11),%%zmm13,%%zmm15		\n\t"\
		"vmovaps	%%zmm12,    (%%r10)		\n\t	vmovaps (%%rsi),%%zmm12			\n\t"/* spill zmm12 to free up reg for 2.0 */\
		"vsubpd	%%zmm6 ,%%zmm2,%%zmm2		\n\t	vsubpd	%%zmm14,%%zmm10,%%zmm10		\n\t"\
		"vsubpd	%%zmm7 ,%%zmm3,%%zmm3		\n\t	vsubpd	%%zmm15,%%zmm11,%%zmm11		\n\t"\
		"vsubpd	%%zmm5 ,%%zmm0,%%zmm0		\n\t	vsubpd	%%zmm13,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd	%%zmm4 ,%%zmm1,%%zmm1		\n\t	vsubpd	(%%r10),%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd %%zmm12,%%zmm2,%%zmm6		\n\t	vfmadd132pd	%%zmm12,%%zmm10,%%zmm14		\n\t"\
	"vfmadd132pd %%zmm12,%%zmm3,%%zmm7		\n\t	vfmadd132pd	%%zmm12,%%zmm11,%%zmm15		\n\t"\
	"vfmadd132pd %%zmm12,%%zmm0,%%zmm5		\n\t	vfmadd132pd	%%zmm12,%%zmm8 ,%%zmm13		\n\t"\
	"vfmadd132pd %%zmm12,%%zmm1,%%zmm4		\n\t	vfmadd132pd	(%%r10),%%zmm9 ,%%zmm12		\n\t"\
		"vmovaps	%%zmm6,    (%%r10)		\n\t	vmovaps	%%zmm7,0x040(%%r10)		\n\t"/* Write B0 to free up 2 regs */\
		/* B2 = t0*~w2 = t0*[c2-I.s1] */\
		"movq	%[__cc0],%%rsi				\n\t	vmovaps %%zmm14,0x200(%%r10)			\n\t"\
/* cc0+18: */"addq	$0x480,%%rsi			\n\t	vmovaps %%zmm15,0x240(%%r10)			\n\t"\
		"vmovaps	%%zmm2,%%zmm6			\n\t	vmovaps	%%zmm10,%%zmm14					\n\t"\
		"vmovaps	%%zmm3,%%zmm7			\n\t	vmovaps	%%zmm11,%%zmm15					\n\t"\
		"vmulpd	     (%%rsi),%%zmm6,%%zmm6	\n\t	vmulpd	0x180(%%rsi),%%zmm14,%%zmm14	\n\t"\
		"vmulpd	     (%%rsi),%%zmm7,%%zmm7	\n\t	vmulpd	0x180(%%rsi),%%zmm15,%%zmm15	\n\t"\
	" vfmadd231pd 0x40(%%rsi),%%zmm3,%%zmm6	\n\t  vfmadd231pd 0x1c0(%%rsi),%%zmm11,%%zmm14	\n\t"\
	"vfnmadd231pd 0x40(%%rsi),%%zmm2,%%zmm7	\n\t vfnmadd231pd 0x1c0(%%rsi),%%zmm10,%%zmm15	\n\t"\
		"vmovaps	%%zmm6,0x100(%%r10)		\n\t	vmovaps %%zmm14,0x300(%%r10)			\n\t"\
		"vmovaps	%%zmm7,0x140(%%r10)		\n\t	vmovaps %%zmm15,0x340(%%r10)			\n\t"\
		/* B1 = t1*~w1 = t1*[c1-I.s1] */\
		"vmovaps	0x080(%%rsi),%%zmm2		\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm3		\n\t"\
		"vmovaps	%%zmm5,%%zmm6			\n\t	vmovaps	%%zmm13,%%zmm14					\n\t"\
		"vmovaps	%%zmm1,%%zmm7			\n\t	vmovaps	%%zmm9 ,%%zmm15					\n\t"\
		"vmulpd	%%zmm2,%%zmm6,%%zmm6		\n\t	vmulpd	0x200(%%rsi),%%zmm14,%%zmm14	\n\t"\
		"vmulpd	%%zmm2,%%zmm7,%%zmm7		\n\t	vmulpd	0x200(%%rsi),%%zmm15,%%zmm15	\n\t"\
	" vfmadd231pd	%%zmm3,%%zmm1,%%zmm6	\n\t  vfmadd231pd 0x240(%%rsi),%%zmm9 ,%%zmm14	\n\t"\
	"vfnmadd231pd	%%zmm3,%%zmm5,%%zmm7	\n\t vfnmadd231pd 0x240(%%rsi),%%zmm13,%%zmm15	\n\t"\
		"vmovaps	0x100(%%rsi),%%zmm2		\n\t"\
		"vmovaps	0x140(%%rsi),%%zmm3		\n\t"\
		"vmovaps	%%zmm6,0x080(%%r10)		\n\t	vmovaps %%zmm14,0x280(%%r10)			\n\t"\
		"vmovaps	%%zmm7,0x0c0(%%r10)		\n\t	vmovaps %%zmm15,0x2c0(%%r10)			\n\t"\
		/* B3 = t3*~w3 = t3*[c3-I.s3] */\
		"vmovaps	%%zmm0,%%zmm6			\n\t	vmovaps	%%zmm8 ,%%zmm14					\n\t"\
		"vmovaps	%%zmm4,%%zmm7			\n\t	vmovaps	%%zmm12,%%zmm15					\n\t"\
		"vmulpd	%%zmm2,%%zmm6,%%zmm6		\n\t	vmulpd	0x280(%%rsi),%%zmm14,%%zmm14	\n\t"\
		"vmulpd	%%zmm2,%%zmm7,%%zmm7		\n\t	vmulpd	0x280(%%rsi),%%zmm15,%%zmm15	\n\t"\
	" vfmadd231pd	%%zmm3,%%zmm4,%%zmm6	\n\t  vfmadd231pd 0x2c0(%%rsi),%%zmm12,%%zmm14	\n\t"\
	"vfnmadd231pd	%%zmm3,%%zmm0,%%zmm7	\n\t vfnmadd231pd 0x2c0(%%rsi),%%zmm8 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm6,0x180(%%r10)		\n\t	vmovaps	%%zmm14,0x380(%%r10)			\n\t"\
		"vmovaps	%%zmm7,0x1c0(%%r10)		\n\t	vmovaps	%%zmm15,0x3c0(%%r10)			\n\t"\
	/*
		// Pass 2:
		j = 0x100;	// Set j to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k = p2*8;	// Set k to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		c_tmp = cc0;	// c8,4,C
		// Pass1-DFTs all use same twiddle-triplet, so lit-byte address offset between ptrs to 1st,2nd DFT's twiddles = 0
		i0 = r00; i1 = r00+8; i2 = r00+16; i3 = r00+24;
		o0 = (vec_dbl *)a; o1 = (vec_dbl *)(a+p4); o2 = (vec_dbl *)(a+p8); o3 = (vec_dbl *)(a+pC);
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,   0, o0,o1,o2,o3,k)
	*/\
		"movq	%[__two],%%rsi 				\n\t	movslq	%[__p2],%%r9		\n\t"\
		"movq	%[__add0],%%r10				\n\t	movq	%[__r0],%%rax		\n\t"\
		"movslq	%[__p4]  ,%%r11				\n\t	leaq	(%%r10,%%r11,8),%%r11	\n\t"\
		"movslq	%[__p8]  ,%%r12				\n\t	leaq	(%%r10,%%r12,8),%%r12	\n\t"\
		"movslq	%[__p12] ,%%r13				\n\t	leaq	(%%r10,%%r13,8),%%r13	\n\t"\
		"vmovaps	(%%rsi),%%zmm15	\n\t"/* two */"	leaq	0x100(%%rax),%%rbx	\n\t"\
		"vmovaps	     (%%rax),%%zmm0		\n\t	vmovaps	     (%%rbx),%%zmm8 	\n\t"/* r0 */\
		"vmovaps	0x040(%%rax),%%zmm1		\n\t	vmovaps	0x040(%%rbx),%%zmm9 	\n\t"\
		"vmovaps	0x200(%%rax),%%zmm2		\n\t	vmovaps	0x200(%%rbx),%%zmm10	\n\t"/* r0 + 8 */\
		"vmovaps	0x240(%%rax),%%zmm3		\n\t	vmovaps	0x240(%%rbx),%%zmm11	\n\t"\
		"vmovaps	0x400(%%rax),%%zmm4		\n\t	vmovaps	0x400(%%rbx),%%zmm12	\n\t"/* r0 + 16 */\
		"vmovaps	0x440(%%rax),%%zmm5		\n\t	vmovaps	0x440(%%rbx),%%zmm13	\n\t"\
		"vmovaps	0x600(%%rax),%%zmm6		\n\t	vmovaps	0x600(%%rbx),%%zmm14	\n\t"/* r0 + 24 */\
		"vmovaps	0x640(%%rax),%%zmm7		\n\t"/*	vmovaps	0x640(%%rbx),%%zmm15	\n\t*/"	leaq	0x640(%%rbx),%%r8	\n\t"\
		"vsubpd	%%zmm2 ,%%zmm0,%%zmm0		\n\t	vsubpd	%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd	%%zmm3 ,%%zmm1,%%zmm1		\n\t	vsubpd	%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd	%%zmm6 ,%%zmm4,%%zmm4		\n\t	vsubpd	%%zmm14,%%zmm12,%%zmm12		\n\t"\
		"vsubpd	%%zmm7 ,%%zmm5,%%zmm5		\n\t	vsubpd	(%%r8) ,%%zmm13,%%zmm13		\n\t"\
	"vfmadd132pd %%zmm15,%%zmm0,%%zmm2		\n\t	vfmadd132pd	%%zmm15,%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd132pd %%zmm15,%%zmm1,%%zmm3		\n\t	vfmadd132pd	%%zmm15,%%zmm9 ,%%zmm11		\n\t"\
	"vfmadd132pd %%zmm15,%%zmm4,%%zmm6		\n\t	vfmadd132pd	%%zmm15,%%zmm12,%%zmm14		\n\t"\
	"vfmadd132pd %%zmm15,%%zmm5,%%zmm7		\n\t	vfmadd132pd	(%%r8),%%zmm13,%%zmm15		\n\t"\
		"vmovaps	%%zmm12,    (%%rax)		\n\t	vmovaps (%%rsi),%%zmm12			\n\t"/* spill zmm12 to free up reg for 2.0 */\
		"vsubpd	%%zmm6 ,%%zmm2,%%zmm2		\n\t	vsubpd	%%zmm14,%%zmm10,%%zmm10		\n\t"\
		"vsubpd	%%zmm7 ,%%zmm3,%%zmm3		\n\t	vsubpd	%%zmm15,%%zmm11,%%zmm11		\n\t"\
		"vsubpd	%%zmm5 ,%%zmm0,%%zmm0		\n\t	vsubpd	%%zmm13,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd	%%zmm4 ,%%zmm1,%%zmm1		\n\t	vsubpd	(%%rax),%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd %%zmm12,%%zmm2,%%zmm6		\n\t	vfmadd132pd	%%zmm12,%%zmm10,%%zmm14		\n\t"\
	"vfmadd132pd %%zmm12,%%zmm3,%%zmm7		\n\t	vfmadd132pd	%%zmm12,%%zmm11,%%zmm15		\n\t"\
	"vfmadd132pd %%zmm12,%%zmm0,%%zmm5		\n\t	vfmadd132pd	%%zmm12,%%zmm8 ,%%zmm13		\n\t"\
	"vfmadd132pd %%zmm12,%%zmm1,%%zmm4		\n\t	vfmadd132pd	(%%rax),%%zmm9 ,%%zmm12		\n\t"\
		"vmovaps	%%zmm6,    (%%r10)		\n\t	vmovaps	%%zmm7,0x040(%%r10)	\n\t"/* Write B0 to free up 2 regs */\
		/* B2 = t0*~w2 = t0*[c2-I.s1] */"			shlq	$3,%%r9				\n\t"/* p2*8: */\
		"vmovaps	%%zmm0,    (%%rax)		\n\t	vmovaps	%%zmm4,0x040(%%rax)	\n\t"/* Write mm0,4 to free up 2 regs */\
		"movq	%[__cc0],%%rsi			 	\n\t	vmovaps %%zmm14,     (%%r10,%%r9)		\n\t"\
		"											vmovaps %%zmm15,0x040(%%r10,%%r9)		\n\t"\
		"vmovaps	     (%%rsi),%%zmm0		\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm4		\n\t"\
		"vmovaps	%%zmm2,%%zmm6			\n\t	vmovaps	%%zmm10,%%zmm14					\n\t"\
		"vmovaps	%%zmm3,%%zmm7			\n\t	vmovaps	%%zmm11,%%zmm15					\n\t"\
		"vmulpd	  %%zmm0,%%zmm6,%%zmm6		\n\t	vmulpd	  %%zmm0,%%zmm14,%%zmm14	\n\t"\
		"vmulpd	  %%zmm0,%%zmm7,%%zmm7		\n\t	vmulpd	  %%zmm0,%%zmm15,%%zmm15	\n\t"\
	" vfmadd231pd %%zmm4,%%zmm3,%%zmm6		\n\t  vfmadd231pd %%zmm4,%%zmm11,%%zmm14	\n\t"\
	"vfnmadd231pd %%zmm4,%%zmm2,%%zmm7		\n\t vfnmadd231pd %%zmm4,%%zmm10,%%zmm15	\n\t"\
		"vmovaps	%%zmm6,     (%%r12)		\n\t	vmovaps %%zmm14,     (%%r12,%%r9)		\n\t"\
		"vmovaps	%%zmm7,0x040(%%r12)		\n\t	vmovaps %%zmm15,0x040(%%r12,%%r9)		\n\t"\
		"vmovaps	    (%%rax),%%zmm0		\n\t	vmovaps	0x040(%%rax),%%zmm4	\n\t"/* Restore spill of mm0,4 */\
		/* B1 = t1*~w1 = t1*[c1-I.s1] */\
		"vmovaps	0x080(%%rsi),%%zmm2		\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm3		\n\t"\
		"vmovaps	     %%zmm5 ,%%zmm6		\n\t	vmovaps	%%zmm13,%%zmm14				\n\t"\
		"vmovaps	     %%zmm1 ,%%zmm7		\n\t	vmovaps	%%zmm9 ,%%zmm15				\n\t"\
		"vmulpd	  %%zmm2,%%zmm6,%%zmm6		\n\t	vmulpd	  %%zmm2,%%zmm14,%%zmm14	\n\t"\
		"vmulpd	  %%zmm2,%%zmm7,%%zmm7		\n\t	vmulpd	  %%zmm2,%%zmm15,%%zmm15	\n\t"\
	" vfmadd231pd %%zmm3,%%zmm1,%%zmm6		\n\t  vfmadd231pd %%zmm3,%%zmm9 ,%%zmm14	\n\t"\
	"vfnmadd231pd %%zmm3,%%zmm5,%%zmm7		\n\t vfnmadd231pd %%zmm3,%%zmm13,%%zmm15	\n\t"\
		"vmovaps	%%zmm6,     (%%r11)		\n\t	vmovaps %%zmm14,     (%%r11,%%r9)	\n\t"\
		"vmovaps	%%zmm7,0x040(%%r11)		\n\t	vmovaps %%zmm15,0x040(%%r11,%%r9)	\n\t"\
		/* B3 = t3*~w3 = t3*[c3-I.s3] */\
		"vmovaps	0x100(%%rsi),%%zmm2		\n\t"\
		"vmovaps	0x140(%%rsi),%%zmm3		\n\t"\
		"vmovaps	%%zmm0,%%zmm6			\n\t	vmovaps	%%zmm8 ,%%zmm14				\n\t"\
		"vmovaps	%%zmm4,%%zmm7			\n\t	vmovaps	%%zmm12,%%zmm15				\n\t"\
		"vmulpd	  %%zmm2,%%zmm6,%%zmm6		\n\t	vmulpd	  %%zmm2,%%zmm14,%%zmm14	\n\t"\
		"vmulpd	  %%zmm2,%%zmm7,%%zmm7		\n\t	vmulpd	  %%zmm2,%%zmm15,%%zmm15	\n\t"\
	" vfmadd231pd %%zmm3,%%zmm4,%%zmm6		\n\t  vfmadd231pd %%zmm3,%%zmm12,%%zmm14	\n\t"\
	"vfnmadd231pd %%zmm3,%%zmm0,%%zmm7		\n\t vfnmadd231pd %%zmm3,%%zmm8 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm6,     (%%r13)		\n\t	vmovaps	%%zmm14,     (%%r13,%%r9)	\n\t"\
		"vmovaps	%%zmm7,0x040(%%r13)		\n\t	vmovaps	%%zmm15,0x040(%%r13,%%r9)	\n\t"\
	/*
		i0 += 2; i1 += 2; i2 += 2; i3 += 2;
		o0 = (vec_dbl *)(a+p1); o1 = (vec_dbl *)(a+p5); o2 = (vec_dbl *)(a+p9); o3 = (vec_dbl *)(a+pD);
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,   0, o0,o1,o2,o3,k)
	*/\
		"movq	%[__two],%%rsi 		\n\t	movslq	%[__p1],%%rdi		\n\t"\
		"leaq (%%r10,%%rdi,8),%%r10	\n\t	addq	$0x080,%%rax		\n\t"\
		"leaq (%%r11,%%rdi,8),%%r11	\n\t"\
		"leaq (%%r12,%%rdi,8),%%r12	\n\t"\
		"leaq (%%r13,%%rdi,8),%%r13	\n\t"\
		"vmovaps	(%%rsi),%%zmm15	\n\t"/* two */"	leaq	0x100(%%rax),%%rbx	\n\t"\
		"vmovaps	     (%%rax),%%zmm0		\n\t	vmovaps	     (%%rbx),%%zmm8 	\n\t"/* r0 */\
		"vmovaps	0x040(%%rax),%%zmm1		\n\t	vmovaps	0x040(%%rbx),%%zmm9 	\n\t"\
		"vmovaps	0x200(%%rax),%%zmm2		\n\t	vmovaps	0x200(%%rbx),%%zmm10	\n\t"/* r0 + 8 */\
		"vmovaps	0x240(%%rax),%%zmm3		\n\t	vmovaps	0x240(%%rbx),%%zmm11	\n\t"\
		"vmovaps	0x400(%%rax),%%zmm4		\n\t	vmovaps	0x400(%%rbx),%%zmm12	\n\t"/* r0 + 16 */\
		"vmovaps	0x440(%%rax),%%zmm5		\n\t	vmovaps	0x440(%%rbx),%%zmm13	\n\t"\
		"vmovaps	0x600(%%rax),%%zmm6		\n\t	vmovaps	0x600(%%rbx),%%zmm14	\n\t"/* r0 + 24 */\
		"vmovaps	0x640(%%rax),%%zmm7		\n\t"/*	vmovaps	0x640(%%rbx),%%zmm15	\n\t*/"	leaq	0x640(%%rbx),%%r8	\n\t"\
		"vsubpd	%%zmm2 ,%%zmm0,%%zmm0		\n\t	vsubpd	%%zmm10,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd	%%zmm3 ,%%zmm1,%%zmm1		\n\t	vsubpd	%%zmm11,%%zmm9 ,%%zmm9 		\n\t"\
		"vsubpd	%%zmm6 ,%%zmm4,%%zmm4		\n\t	vsubpd	%%zmm14,%%zmm12,%%zmm12		\n\t"\
		"vsubpd	%%zmm7 ,%%zmm5,%%zmm5		\n\t	vsubpd	(%%r8) ,%%zmm13,%%zmm13		\n\t"\
	"vfmadd132pd %%zmm15,%%zmm0,%%zmm2		\n\t	vfmadd132pd	%%zmm15,%%zmm8 ,%%zmm10		\n\t"\
	"vfmadd132pd %%zmm15,%%zmm1,%%zmm3		\n\t	vfmadd132pd	%%zmm15,%%zmm9 ,%%zmm11		\n\t"\
	"vfmadd132pd %%zmm15,%%zmm4,%%zmm6		\n\t	vfmadd132pd	%%zmm15,%%zmm12,%%zmm14		\n\t"\
	"vfmadd132pd %%zmm15,%%zmm5,%%zmm7		\n\t	vfmadd132pd	(%%r8),%%zmm13,%%zmm15		\n\t"\
		"vmovaps	%%zmm12,    (%%rax)		\n\t	vmovaps (%%rsi),%%zmm12			\n\t"/* spill zmm12 to free up reg for 2.0 */\
		"vsubpd	%%zmm6 ,%%zmm2,%%zmm2		\n\t	vsubpd	%%zmm14,%%zmm10,%%zmm10		\n\t"\
		"vsubpd	%%zmm7 ,%%zmm3,%%zmm3		\n\t	vsubpd	%%zmm15,%%zmm11,%%zmm11		\n\t"\
		"vsubpd	%%zmm5 ,%%zmm0,%%zmm0		\n\t	vsubpd	%%zmm13,%%zmm8 ,%%zmm8 		\n\t"\
		"vsubpd	%%zmm4 ,%%zmm1,%%zmm1		\n\t	vsubpd	(%%rax),%%zmm9 ,%%zmm9 		\n\t"\
	"vfmadd132pd %%zmm12,%%zmm2,%%zmm6		\n\t	vfmadd132pd	%%zmm12,%%zmm10,%%zmm14		\n\t"\
	"vfmadd132pd %%zmm12,%%zmm3,%%zmm7		\n\t	vfmadd132pd	%%zmm12,%%zmm11,%%zmm15		\n\t"\
	"vfmadd132pd %%zmm12,%%zmm0,%%zmm5		\n\t	vfmadd132pd	%%zmm12,%%zmm8 ,%%zmm13		\n\t"\
	"vfmadd132pd %%zmm12,%%zmm1,%%zmm4		\n\t	vfmadd132pd	(%%rax),%%zmm9 ,%%zmm12		\n\t"\
		"vmovaps	%%zmm6,    (%%r10)		\n\t	vmovaps	%%zmm7,0x040(%%r10)	\n\t"/* Write B0 to free up 2 regs */\
		/* B2 = t0*~w2 = t0*[c2-I.s1] */\
		"vmovaps	%%zmm0,    (%%rax)		\n\t	vmovaps	%%zmm4,0x040(%%rax)	\n\t"/* Write mm0,4 to free up 2 regs */\
		"movq	%[__cc0],%%rsi			 	\n\t	vmovaps %%zmm14,     (%%r10,%%r9)		\n\t"\
		"											vmovaps %%zmm15,0x040(%%r10,%%r9)		\n\t"\
		"vmovaps	     (%%rsi),%%zmm0		\n\t"\
		"vmovaps	0x040(%%rsi),%%zmm4		\n\t"\
		"vmovaps	%%zmm2,%%zmm6			\n\t	vmovaps	%%zmm10,%%zmm14					\n\t"\
		"vmovaps	%%zmm3,%%zmm7			\n\t	vmovaps	%%zmm11,%%zmm15					\n\t"\
		"vmulpd	  %%zmm0,%%zmm6,%%zmm6		\n\t	vmulpd	  %%zmm0,%%zmm14,%%zmm14	\n\t"\
		"vmulpd	  %%zmm0,%%zmm7,%%zmm7		\n\t	vmulpd	  %%zmm0,%%zmm15,%%zmm15	\n\t"\
	" vfmadd231pd %%zmm4,%%zmm3,%%zmm6		\n\t  vfmadd231pd %%zmm4,%%zmm11,%%zmm14	\n\t"\
	"vfnmadd231pd %%zmm4,%%zmm2,%%zmm7		\n\t vfnmadd231pd %%zmm4,%%zmm10,%%zmm15	\n\t"\
		"vmovaps	%%zmm6,     (%%r12)		\n\t	vmovaps %%zmm14,     (%%r12,%%r9)		\n\t"\
		"vmovaps	%%zmm7,0x040(%%r12)		\n\t	vmovaps %%zmm15,0x040(%%r12,%%r9)		\n\t"\
		"vmovaps	    (%%rax),%%zmm0		\n\t	vmovaps	0x040(%%rax),%%zmm4	\n\t"/* Restore spill of mm0,4 */\
		/* B1 = t1*~w1 = t1*[c1-I.s1] */\
		"vmovaps	0x080(%%rsi),%%zmm2		\n\t"\
		"vmovaps	0x0c0(%%rsi),%%zmm3		\n\t"\
		"vmovaps	     %%zmm5 ,%%zmm6		\n\t	vmovaps	%%zmm13,%%zmm14				\n\t"\
		"vmovaps	     %%zmm1 ,%%zmm7		\n\t	vmovaps	%%zmm9 ,%%zmm15				\n\t"\
		"vmulpd	  %%zmm2,%%zmm6,%%zmm6		\n\t	vmulpd	  %%zmm2,%%zmm14,%%zmm14	\n\t"\
		"vmulpd	  %%zmm2,%%zmm7,%%zmm7		\n\t	vmulpd	  %%zmm2,%%zmm15,%%zmm15	\n\t"\
	" vfmadd231pd %%zmm3,%%zmm1,%%zmm6		\n\t  vfmadd231pd %%zmm3,%%zmm9 ,%%zmm14	\n\t"\
	"vfnmadd231pd %%zmm3,%%zmm5,%%zmm7		\n\t vfnmadd231pd %%zmm3,%%zmm13,%%zmm15	\n\t"\
		"vmovaps	%%zmm6,     (%%r11)		\n\t	vmovaps %%zmm14,     (%%r11,%%r9)	\n\t"\
		"vmovaps	%%zmm7,0x040(%%r11)		\n\t	vmovaps %%zmm15,0x040(%%r11,%%r9)	\n\t"\
		/* B3 = t3*~w3 = t3*[c3-I.s3] */\
		"vmovaps	0x100(%%rsi),%%zmm2		\n\t"\
		"vmovaps	0x140(%%rsi),%%zmm3		\n\t"\
		"vmovaps	%%zmm0,%%zmm6			\n\t	vmovaps	%%zmm8 ,%%zmm14				\n\t"\
		"vmovaps	%%zmm4,%%zmm7			\n\t	vmovaps	%%zmm12,%%zmm15				\n\t"\
		"vmulpd	  %%zmm2,%%zmm6,%%zmm6		\n\t	vmulpd	  %%zmm2,%%zmm14,%%zmm14	\n\t"\
		"vmulpd	  %%zmm2,%%zmm7,%%zmm7		\n\t	vmulpd	  %%zmm2,%%zmm15,%%zmm15	\n\t"\
	" vfmadd231pd %%zmm3,%%zmm4,%%zmm6		\n\t  vfmadd231pd %%zmm3,%%zmm12,%%zmm14	\n\t"\
	"vfnmadd231pd %%zmm3,%%zmm0,%%zmm7		\n\t vfnmadd231pd %%zmm3,%%zmm8 ,%%zmm15	\n\t"\
		"vmovaps	%%zmm6,     (%%r13)		\n\t	vmovaps	%%zmm14,     (%%r13,%%r9)	\n\t"\
		"vmovaps	%%zmm7,0x040(%%r13)		\n\t	vmovaps	%%zmm15,0x040(%%r13,%%r9)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r0] "m" (Xr0)\
		 ,[__two] "m" (Xtwo)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_AVX2)	// AVX and AVX2 both use 256-bit registers

	// The non-name-suffixed DIF,DIT macros in AVX2 mode are based on the AVX-mode _V2 macros:
	// 96 ADD, 80 FMA, 48 pure-MUL, 316 MOVAPS
	#define SSE2_RADIX16_DIF_TWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr0,Xtwo,Xcc0,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/*
		k1 = p2*8;	// Set k1 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k2 = 0x80;	// Set k2 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		// Pass1-DFTs all use same twiddle-triplet, so lit-byte address offset between ptrs to 1st,2nd DFT's twiddles = 0
		i0 = (vec_dbl *)add0; i1 = (vec_dbl *)(add0+p4); i2 = (vec_dbl *)(add0+p8); i3 = (vec_dbl *)(add0+p12);
		o0 = r1; o1 = r1+2; o2 = r1+4; o3 = r1+6;
		c_tmp = cc0;	// c8,4,C
		SSE2_RADIX_04_DIF_3TWIDDLE_X1(i0,i1,i2,i3,k1, two,c_tmp, o0,o1,o2,o3,k2)
	*/\
		"movq	%[__cc0],%%rsi 			\n\t	movslq	%[__p2],%%rdi		\n\t"\
		"movq	%[__add0],%%rax			\n\t	movq	%[__r0],%%r10		\n\t"\
		"movslq	%[__p4],%%rbx			\n\t	leaq	(%%rax,%%rbx,8),%%rbx	\n\t"\
		"movslq	%[__p8],%%rcx			\n\t	leaq	(%%rax,%%rcx,8),%%rcx	\n\t"\
		"movslq	%[__p12],%%rdx			\n\t	leaq	(%%rax,%%rdx,8),%%rdx	\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%rax)\n\t"/* [base-address + data-fetch-ahead index] + p0 */\
		/* Do the p0,p2 combo: */\
		"vmovaps	     (%%rcx),%%ymm4	\n\t	shlq	$3,%%rdi			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rsi),%%ymm10\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm11\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t	vmovaps	     (%%rcx,%%rdi),%%ymm12	\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	0x020(%%rcx,%%rdi),%%ymm13	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t	vmovaps	     (%%rax,%%rdi),%%ymm8 	\n\t	vmovaps	%%ymm12,%%ymm14		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t	vmovaps	0x020(%%rax,%%rdi),%%ymm9 	\n\t	vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmulpd	%%ymm10,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12		\n\t	vmovaps	%%ymm0,%%ymm2		\n\t"\
		"vmulpd	%%ymm10,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13		\n\t	vmovaps	%%ymm1,%%ymm3		\n\t"\
	"vfnmadd231pd %%ymm11,%%ymm7,%%ymm4	\n\t vfnmadd231pd %%ymm11,%%ymm15,%%ymm12	\n\t	vmovaps	%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd231pd  %%ymm11,%%ymm6,%%ymm5	\n\t vfmadd231pd  %%ymm11,%%ymm14,%%ymm13	\n\t	vmovaps	%%ymm9 ,%%ymm11		\n\t"\
		"vaddpd	%%ymm4 ,%%ymm0,%%ymm0	\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm1,%%ymm1	\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm4 ,%%ymm2,%%ymm2	\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm3,%%ymm3	\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,     (%%r10)	\n\t	vmovaps	%%ymm8 ,0x100(%%r10)		\n\t"/* Spill 1: free up ymm0,1 */\
		"vmovaps	%%ymm1,0x020(%%r10)	\n\t	vmovaps	%%ymm9 ,0x120(%%r10)		\n\t"/* Do the p1,3 combo: */\
		"vmovaps	0x080(%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm1	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t	vmovaps	     (%%rdx,%%rdi),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t	vmovaps	0x020(%%rdx,%%rdi),%%ymm15	\n\t"\
		"vmovaps	%%ymm6,%%ymm4		\n\t	vmovaps	%%ymm14,%%ymm12				\n\t"\
		"vmovaps	%%ymm7,%%ymm5		\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm0 ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm0 ,%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd %%ymm1,%%ymm7,%%ymm4	\n\t vfnmadd231pd %%ymm1,%%ymm15,%%ymm12	\n\t"\
	"vfmadd231pd  %%ymm1,%%ymm6,%%ymm5	\n\t vfmadd231pd  %%ymm1,%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm5,0x0a0(%%r10)	\n\t	vmovaps	%%ymm13,0x1a0(%%r10)		\n\t"/* Spill 2 */\
		"vmovaps	%%ymm4,0x080(%%r10)	\n\t	vmovaps	%%ymm12,0x180(%%r10)		\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm8	\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm9	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6	\n\t	vmovaps	     (%%rbx,%%rdi),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t	vmovaps	0x020(%%rbx,%%rdi),%%ymm15	\n\t"\
		"vmovaps	%%ymm6,%%ymm4		\n\t	vmovaps	%%ymm14,%%ymm12				\n\t"\
		"vmovaps	%%ymm7,%%ymm5		\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"\
		"vmulpd	%%ymm8,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	%%ymm8,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd %%ymm9,%%ymm7,%%ymm4	\n\t vfnmadd231pd %%ymm9,%%ymm15,%%ymm12	\n\t"\
	"vfmadd231pd  %%ymm9,%%ymm6,%%ymm5	\n\t vfmadd231pd  %%ymm9,%%ymm14,%%ymm13	\n\t"\
		"vmovaps	0x080(%%r10),%%ymm0	\n\t	vmovaps	0x180(%%r10),%%ymm8 		\n\t"/* Restore 2 */\
		"vmovaps	0x0a0(%%r10),%%ymm1	\n\t	vmovaps	0x1a0(%%r10),%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	%%ymm13,%%ymm15				\n\t	 movq	%[__two],%%rsi	\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t	vmovaps	%%ymm12,%%ymm14				\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5	\n\t	vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		/* Finish radix-4 butterfly and store results: */\
		"vmovaps	     (%%r10),%%ymm0	\n\t	vmovaps	0x100(%%r10),%%ymm8 		\n\t"/* Restore 1 */\
		"vmovaps	0x020(%%r10),%%ymm1	\n\t	vmovaps	0x120(%%r10),%%ymm9 		\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1	\n\t	vsubpd	%%ymm15,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,0x080(%%r10)	\n\t	vmovaps	%%ymm8 ,0x180(%%r10)		\n\t"	/* 2.0, shared by both columns: */\
		"vmovaps	%%ymm2,0x040(%%r10)	\n\t	vmovaps	%%ymm10,0x140(%%r10)		\n\t	vmovaps	(%%rsi),%%ymm0	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%r10)	\n\t	vmovaps	%%ymm9 ,0x1a0(%%r10)		\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%r10)	\n\t	vmovaps	%%ymm11,0x1e0(%%r10)		\n\t"\
	"vfmadd213pd 0x80(%%r10),%%ymm0,%%ymm6 \n\t	vfmadd213pd	%%ymm8 ,%%ymm0,%%ymm14	\n\t"\
	"vfmadd213pd      %%ymm2,%%ymm0,%%ymm5 \n\t	vfmadd213pd	%%ymm10,%%ymm0,%%ymm13	\n\t"\
	"vfmadd213pd      %%ymm1,%%ymm0,%%ymm7 \n\t	vfmadd213pd	%%ymm9 ,%%ymm0,%%ymm15	\n\t"\
	"vfmadd213pd      %%ymm3,%%ymm0,%%ymm4 \n\t	vfmadd213pd	%%ymm11,%%ymm0,%%ymm12	\n\t"\
		"vmovaps	%%ymm6,     (%%r10)	\n\t	vmovaps	%%ymm14,0x100(%%r10)		\n\t"/* don't need reload-from-mem of ymm8/0x180(%%r10) as we do in lcol, but 1 cycle faster with it. [!?] */\
		"vmovaps	%%ymm5,0x0c0(%%r10)	\n\t	vmovaps	%%ymm13,0x1c0(%%r10)		\n\t"\
		"vmovaps	%%ymm7,0x020(%%r10)	\n\t	vmovaps	%%ymm15,0x120(%%r10)		\n\t"\
		"vmovaps	%%ymm4,0x060(%%r10)	\n\t	vmovaps	%%ymm12,0x160(%%r10)		\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%rcx)\n\t"/* ... + p8 */\
	/*
		i0 = (vec_dbl *)(add0+p1); i1 = (vec_dbl *)(add0+p4+p1); i2 = (vec_dbl *)(add0+p8+p1); i3 = (vec_dbl *)(add0+p12+p1);
		o0 += 16; o1 += 16; o2 += 16; o3 += 16;
		SSE2_RADIX_04_DIF_3TWIDDLE_X1(i0,i1,i2,i3,k1, two,c_tmp, o0,o1,o2,o3,k2)
	*/\
		"movslq	%[__p1],%%r9			\n\t	addq	$0x200,%%r10		\n\t"\
		"leaq	(%%rax,%%r9,8),%%rax	\n\t	movq	%[__cc0],%%rsi 		\n\t"/* repoint rsi from two -> cc0 */\
		"leaq	(%%rbx,%%r9,8),%%rbx	\n\t"\
		"leaq	(%%rcx,%%r9,8),%%rcx	\n\t"\
		"leaq	(%%rdx,%%r9,8),%%rdx	\n\t"\
		/* Do the p0,p2 combo: */\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rsi),%%ymm10\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm11\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t	vmovaps	     (%%rcx,%%rdi),%%ymm12	\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	0x020(%%rcx,%%rdi),%%ymm13	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t	vmovaps	     (%%rax,%%rdi),%%ymm8 	\n\t	vmovaps	%%ymm12,%%ymm14		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t	vmovaps	0x020(%%rax,%%rdi),%%ymm9 	\n\t	vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmulpd	%%ymm10,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12		\n\t	vmovaps	%%ymm0,%%ymm2		\n\t"\
		"vmulpd	%%ymm10,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13		\n\t	vmovaps	%%ymm1,%%ymm3		\n\t"\
	"vfnmadd231pd %%ymm11,%%ymm7,%%ymm4	\n\t vfnmadd231pd %%ymm11,%%ymm15,%%ymm12	\n\t	vmovaps	%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd231pd  %%ymm11,%%ymm6,%%ymm5	\n\t vfmadd231pd  %%ymm11,%%ymm14,%%ymm13	\n\t	vmovaps	%%ymm9 ,%%ymm11		\n\t"\
		"vaddpd	%%ymm4 ,%%ymm0,%%ymm0	\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm1,%%ymm1	\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm4 ,%%ymm2,%%ymm2	\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm3,%%ymm3	\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,     (%%r10)	\n\t	vmovaps	%%ymm8 ,0x100(%%r10)		\n\t"/* Spill 1: free up ymm0,1 */\
		"vmovaps	%%ymm1,0x020(%%r10)	\n\t	vmovaps	%%ymm9 ,0x120(%%r10)		\n\t"/* Do the p1,3 combo: */\
		"vmovaps	0x080(%%rsi),%%ymm0	\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm1	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t	vmovaps	     (%%rdx,%%rdi),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t	vmovaps	0x020(%%rdx,%%rdi),%%ymm15	\n\t"\
		"vmovaps	%%ymm6,%%ymm4		\n\t	vmovaps	%%ymm14,%%ymm12				\n\t"\
		"vmovaps	%%ymm7,%%ymm5		\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm0 ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm0 ,%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd %%ymm1,%%ymm7,%%ymm4	\n\t vfnmadd231pd %%ymm1,%%ymm15,%%ymm12	\n\t"\
	"vfmadd231pd  %%ymm1,%%ymm6,%%ymm5	\n\t vfmadd231pd  %%ymm1,%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm5,0x0a0(%%r10)	\n\t	vmovaps	%%ymm13,0x1a0(%%r10)		\n\t"/* Spill 2 */\
		"vmovaps	%%ymm4,0x080(%%r10)	\n\t	vmovaps	%%ymm12,0x180(%%r10)		\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm8	\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm9	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6	\n\t	vmovaps	     (%%rbx,%%rdi),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t	vmovaps	0x020(%%rbx,%%rdi),%%ymm15	\n\t"\
		"vmovaps	%%ymm6,%%ymm4		\n\t	vmovaps	%%ymm14,%%ymm12				\n\t"\
		"vmovaps	%%ymm7,%%ymm5		\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"\
		"vmulpd	%%ymm8,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	%%ymm8,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd %%ymm9,%%ymm7,%%ymm4	\n\t vfnmadd231pd %%ymm9,%%ymm15,%%ymm12	\n\t"\
	"vfmadd231pd  %%ymm9,%%ymm6,%%ymm5	\n\t vfmadd231pd  %%ymm9,%%ymm14,%%ymm13	\n\t"\
		"vmovaps	0x080(%%r10),%%ymm0	\n\t	vmovaps	0x180(%%r10),%%ymm8 		\n\t"/* Restore 2 */\
		"vmovaps	0x0a0(%%r10),%%ymm1	\n\t	vmovaps	0x1a0(%%r10),%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	%%ymm13,%%ymm15				\n\t	 movq	%[__two],%%rsi	\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t	vmovaps	%%ymm12,%%ymm14				\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5	\n\t	vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		/* Finish radix-4 butterfly and store results: */\
		"vmovaps	     (%%r10),%%ymm0	\n\t	vmovaps	0x100(%%r10),%%ymm8 		\n\t"/* Restore 1 */\
		"vmovaps	0x020(%%r10),%%ymm1	\n\t	vmovaps	0x120(%%r10),%%ymm9 		\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1	\n\t	vsubpd	%%ymm15,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,0x080(%%r10)	\n\t	vmovaps	%%ymm8 ,0x180(%%r10)		\n\t"	/* 2.0, shared by both columns: */\
		"vmovaps	%%ymm2,0x040(%%r10)	\n\t	vmovaps	%%ymm10,0x140(%%r10)		\n\t	vmovaps	(%%rsi),%%ymm0	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%r10)	\n\t	vmovaps	%%ymm9 ,0x1a0(%%r10)		\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%r10)	\n\t	vmovaps	%%ymm11,0x1e0(%%r10)		\n\t"\
	"vfmadd213pd 0x80(%%r10),%%ymm0,%%ymm6 \n\t	vfmadd213pd	%%ymm8 ,%%ymm0,%%ymm14	\n\t"\
	"vfmadd213pd      %%ymm2,%%ymm0,%%ymm5 \n\t	vfmadd213pd	%%ymm10,%%ymm0,%%ymm13	\n\t"\
	"vfmadd213pd      %%ymm1,%%ymm0,%%ymm7 \n\t	vfmadd213pd	%%ymm9 ,%%ymm0,%%ymm15	\n\t"\
	"vfmadd213pd      %%ymm3,%%ymm0,%%ymm4 \n\t	vfmadd213pd	%%ymm11,%%ymm0,%%ymm12	\n\t"\
		"vmovaps	%%ymm6,     (%%r10)	\n\t	vmovaps	%%ymm14,0x100(%%r10)		\n\t"/* don't need reload-from-mem of ymm8/0x180(%%r10) as we do in lcol, but 1 cycle faster with it. [!?] */\
		"vmovaps	%%ymm5,0x0c0(%%r10)	\n\t	vmovaps	%%ymm13,0x1c0(%%r10)		\n\t"\
		"vmovaps	%%ymm7,0x020(%%r10)	\n\t	vmovaps	%%ymm15,0x120(%%r10)		\n\t"\
		"vmovaps	%%ymm4,0x060(%%r10)	\n\t	vmovaps	%%ymm12,0x160(%%r10)		\n\t"\
	/*
		// Pass 2:
		k1 = 0x080;	// Set k1 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k2 = p4*8;	// Set k2 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		// 0x0c0 = literal-bytewise address offset between the pointers to the first and second DFT's twiddles = 0:
		i0 = r1; i1 = r1+16; i2 = r1+8; i3 = r1+24;
		o0 = (vec_dbl *)add0; o1 = (vec_dbl *)(add0+p2); o2 = (vec_dbl *)(add0+p1); o3 = (vec_dbl *)(add0+p3);
		c_tmp += 6;	// c2,A,6
		SSE2_RADIX_04_DIF_3TWIDDLE_X2(i0,i1,i2,i3,k1, two,c_tmp,0x0c0, o0,o1,o2,o3,k2)
	*/\
		/* Take advantage of fact that i1 = r0 + 16 already in %%r10 and o2 = add0+p1 in %%rax : */\
		/* i-offset = 0x080, so inline with r[a-d]x base-address offsets rather than adding via (r[a-d]x,rdi) */\
		"movq	%%r10,%%rbx				\n\t	movq	%%rax,%%r12			\n\t"/* rbx:i1 = r0+16; r12:o2 = add0+p1 */\
		"movq	%[__cc0],%%rsi 			\n\t	addq	$0x0c0,%%rsi 		\n\t"/* repoint rsi from two -> cc0, then increment for Pass 2 */\
		"leaq	 0x0c0(%%rsi),%%r8		\n\t	movslq	%[__p2],%%r11		\n\t"\
		"leaq	-0x200(%%rbx),%%rax		\n\t	movslq	%[__p4],%%r9		\n\t"\
		"leaq	-0x100(%%rbx),%%rcx		\n\t	movq	%[__add0]	,%%r10		\n\t"/* o0 */\
		"leaq	 0x100(%%rbx),%%rdx		\n\t	leaq	(%%r12,%%r11,8),%%r13	\n\t"/* o3 = (add0+p1)+p2 */\
		"vmovaps	     (%%rax),%%ymm0	\n\t	leaq	(%%r10,%%r11,8),%%r11	\n\t"/* o1 = (add0+p2) */\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t	shlq	$3,%%r9			\n\t"/* p4 in byte-offset form */\
		/* Do the p0,p2 combo: */\
		"vmovaps	     (%%rcx),%%ymm4	\n\t	vmovaps	0x080(%%rcx),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t	vmovaps	0x0a0(%%rcx),%%ymm13		\n\t"\
		"vmovaps	     (%%rsi),%%ymm2	\n\t	vmovaps	     (%%r8),%%ymm10			\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3	\n\t	vmovaps	0x020(%%r8),%%ymm11			\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r9)\n\t"/* ... + p4 */\
		"vmovaps	%%ymm4,%%ymm6		\n\t	vmovaps	%%ymm12,%%ymm14				\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	%%ymm13,%%ymm15				\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t	vmovaps	0x080(%%rax),%%ymm8 		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t	vmovaps	0x0a0(%%rax),%%ymm9 		\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13		\n\t	vmovaps	%%ymm0,%%ymm2		\n\t"\
	"vfnmadd231pd %%ymm3,%%ymm7,%%ymm4	\n\t vfnmadd231pd %%ymm11,%%ymm15,%%ymm12	\n\t	vmovaps	%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd231pd  %%ymm3,%%ymm6,%%ymm5	\n\t vfmadd231pd  %%ymm11,%%ymm14,%%ymm13	\n\t	vmovaps	%%ymm1,%%ymm3		\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0	\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t	vmovaps	%%ymm9 ,%%ymm11		\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1	\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3	\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,     (%%r10)	\n\t	vmovaps	%%ymm8 ,     (%%r10,%%r9)	\n\t"/* Spill 1: free up ymm0,1 */\
		"vmovaps	%%ymm1,0x020(%%r10)	\n\t	vmovaps	%%ymm9 ,0x020(%%r10,%%r9)	\n\t"/* Do the p1,3 combo: */\
		"vmovaps	0x080(%%rsi),%%ymm0	\n\t	vmovaps	0x080(%%r8),%%ymm8	 		\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm1	\n\t	vmovaps	0x0a0(%%r8),%%ymm9 			\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t	vmovaps	0x080(%%rdx),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t	vmovaps	0x0a0(%%rdx),%%ymm15		\n\t"\
		"vmovaps	%%ymm6,%%ymm4		\n\t	vmovaps	%%ymm14,%%ymm12				\n\t"\
		"vmovaps	%%ymm7,%%ymm5		\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd %%ymm1,%%ymm7,%%ymm4	\n\t vfnmadd231pd %%ymm9 ,%%ymm15,%%ymm12	\n\t"\
	"vfmadd231pd  %%ymm1,%%ymm6,%%ymm5	\n\t vfmadd231pd  %%ymm9 ,%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm5,0x020(%%r12)	\n\t	vmovaps	%%ymm13,0x020(%%r12,%%r9)	\n\t"/* Spill 2 */\
		"vmovaps	%%ymm4,     (%%r12)	\n\t	vmovaps	%%ymm12,     (%%r12,%%r9)	\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm0	\n\t	vmovaps	0x040(%%r8),%%ymm8	 		\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm1	\n\t	vmovaps	0x060(%%r8),%%ymm9 			\n\t"\
		"vmovaps	     (%%rbx),%%ymm6	\n\t	vmovaps	0x080(%%rbx),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t	vmovaps	0x0a0(%%rbx),%%ymm15		\n\t"\
		"vmovaps	%%ymm6,%%ymm4		\n\t	vmovaps	%%ymm14,%%ymm12				\n\t"\
		"vmovaps	%%ymm7,%%ymm5		\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd %%ymm1,%%ymm7,%%ymm4	\n\t vfnmadd231pd %%ymm9 ,%%ymm15,%%ymm12	\n\t"\
	"vfmadd231pd  %%ymm1,%%ymm6,%%ymm5	\n\t vfmadd231pd  %%ymm9 ,%%ymm14,%%ymm13	\n\t"\
		"vmovaps	     (%%r12),%%ymm0	\n\t	vmovaps	     (%%r12,%%r9),%%ymm8	\n\t"/* Restore 2 */\
		"vmovaps	0x020(%%r12),%%ymm1	\n\t	vmovaps	0x020(%%r12,%%r9),%%ymm9	\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	%%ymm13,%%ymm15				\n\t	 movq	%[__two],%%rsi	\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t	vmovaps	%%ymm12,%%ymm14				\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5	\n\t	vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		/* Finish radix-4 butterfly and store results: */\
		"vmovaps	     (%%r10),%%ymm0	\n\t	vmovaps	     (%%r10,%%r9),%%ymm8	\n\t"/* Restore 1 */\
		"vmovaps	0x020(%%r10),%%ymm1	\n\t	vmovaps	0x020(%%r10,%%r9),%%ymm9	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1	\n\t	vsubpd	%%ymm15,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,     (%%r12)	\n\t	vmovaps	%%ymm8 ,     (%%r12,%%r9)	\n\t"	/* 2.0, shared by both columns ... moved +- until found best cycle count: */\
		"vmovaps	%%ymm2,     (%%r11)	\n\t	vmovaps	%%ymm10,     (%%r11,%%r9)	\n\t	vmovaps	(%%rsi),%%ymm0	\n\t"\
		"vmovaps	%%ymm1,0x020(%%r12)	\n\t	vmovaps	%%ymm9 ,0x020(%%r12,%%r9)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%r13)	\n\t	vmovaps	%%ymm11,0x020(%%r13,%%r9)	\n\t"\
	"vfmadd213pd (%%r12),%%ymm0,%%ymm6	\n\t	vfmadd213pd	%%ymm8 ,%%ymm0,%%ymm14	\n\t"\
	"vfmadd213pd %%ymm2 ,%%ymm0,%%ymm5	\n\t	vfmadd213pd	%%ymm10,%%ymm0,%%ymm13	\n\t"\
	"vfmadd213pd %%ymm1 ,%%ymm0,%%ymm7	\n\t	vfmadd213pd	%%ymm9 ,%%ymm0,%%ymm15	\n\t"\
	"vfmadd213pd %%ymm3 ,%%ymm0,%%ymm4	\n\t	vfmadd213pd	%%ymm11,%%ymm0,%%ymm12	\n\t"\
		"vmovaps	%%ymm6,     (%%r10)	\n\t	vmovaps	%%ymm14,     (%%r10,%%r9)	\n\t"/* don't need reload-from-mem of ymm8/0x180(%%r12,%%r9) as we do in lcol, but 1 cycle faster with it. [!?] */\
		"vmovaps	%%ymm5,     (%%r13)	\n\t	vmovaps	%%ymm13,     (%%r13,%%r9)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%r10)	\n\t	vmovaps	%%ymm15,0x020(%%r10,%%r9)	\n\t"\
		"vmovaps	%%ymm4,0x020(%%r11)	\n\t	vmovaps	%%ymm12,0x020(%%r11,%%r9)	\n\t"\
	/*
		i0 += 2; i1 += 2; i2 += 2; i3 += 2;
		o0 = (vec_dbl *)(add0+p8); o1 = (vec_dbl *)(add0+p8+p2); o2 = (vec_dbl *)(add0+p8+p1); o3 = (vec_dbl *)(add0+p8+p3);
		c_tmp += 12;	// c5,D,3
		SSE2_RADIX_04_DIF_3TWIDDLE_X2(i0,i1,i2,i3,k1, two,c_tmp,0x0c0, o0,o1,o2,o3,k2)
	*/\
		/* Take advantage of fact that i1 = r0 + 16 already in %%r10 and o2 = add0+p1 in %%rax : */\
		"movq	%[__cc0],%%rsi 			\n\t	addq	$0x240,%%rsi 		\n\t"/* repoint rsi from two -> cc0, then increment for Pass, set 2 */\
		"addq	$0x040,%%rax			\n\t	movslq	%[__p8],%%r8		\n\t"\
		"addq	$0x040,%%rbx			\n\t"	/* r9 still has p4<<3 */\
		"addq	$0x040,%%rcx			\n\t	leaq	(%%r10,%%r8,8),%%r10	\n\t"/* o0 =  add0     +p8 */\
		"addq	$0x040,%%rdx			\n\t	leaq	(%%r11,%%r8,8),%%r11	\n\t"/* o1 = (add0+p2)+p8 */\
		"vmovaps	     (%%rax),%%ymm0	\n\t	leaq	(%%r12,%%r8,8),%%r12	\n\t"/* o2 = (add0+p1)+p8 */\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t	leaq	(%%r13,%%r8,8),%%r13	\n\t"/* o3 = (add0+p3)+p8 */\
		"leaq	0x0c0(%%rsi),%%r8		\n\t"\
		/* Do the p0,p2 combo: */\
		"vmovaps	     (%%rcx),%%ymm4	\n\t	vmovaps	0x080(%%rcx),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t	vmovaps	0x0a0(%%rcx),%%ymm13		\n\t"\
		"vmovaps	     (%%rsi),%%ymm2	\n\t	vmovaps	     (%%r8),%%ymm10			\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3	\n\t	vmovaps	0x020(%%r8),%%ymm11			\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r9)\n\t"/* ... + pC */\
		"vmovaps	%%ymm4,%%ymm6		\n\t	vmovaps	%%ymm12,%%ymm14				\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	%%ymm13,%%ymm15				\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t	vmovaps	0x080(%%rax),%%ymm8 		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t	vmovaps	0x0a0(%%rax),%%ymm9 		\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13		\n\t	vmovaps	%%ymm0,%%ymm2		\n\t"\
	"vfnmadd231pd %%ymm3,%%ymm7,%%ymm4	\n\t vfnmadd231pd %%ymm11,%%ymm15,%%ymm12	\n\t	vmovaps	%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd231pd  %%ymm3,%%ymm6,%%ymm5	\n\t vfmadd231pd  %%ymm11,%%ymm14,%%ymm13	\n\t	vmovaps	%%ymm1,%%ymm3		\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0	\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t	vmovaps	%%ymm9 ,%%ymm11		\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1	\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3	\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,     (%%r10)	\n\t	vmovaps	%%ymm8 ,     (%%r10,%%r9)	\n\t"/* Spill 1: free up ymm0,1 */\
		"vmovaps	%%ymm1,0x020(%%r10)	\n\t	vmovaps	%%ymm9 ,0x020(%%r10,%%r9)	\n\t"/* Do the p1,3 combo: */\
		"vmovaps	0x080(%%rsi),%%ymm0	\n\t	vmovaps	0x080(%%r8),%%ymm8	 		\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm1	\n\t	vmovaps	0x0a0(%%r8),%%ymm9 			\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t	vmovaps	0x080(%%rdx),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t	vmovaps	0x0a0(%%rdx),%%ymm15		\n\t"\
		"vmovaps	%%ymm6,%%ymm4		\n\t	vmovaps	%%ymm14,%%ymm12				\n\t"\
		"vmovaps	%%ymm7,%%ymm5		\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd %%ymm1,%%ymm7,%%ymm4	\n\t vfnmadd231pd %%ymm9 ,%%ymm15,%%ymm12	\n\t"\
	"vfmadd231pd  %%ymm1,%%ymm6,%%ymm5	\n\t vfmadd231pd  %%ymm9 ,%%ymm14,%%ymm13	\n\t"\
		"vmovaps	%%ymm5,0x020(%%r12)	\n\t	vmovaps	%%ymm13,0x020(%%r12,%%r9)	\n\t"/* Spill 2 */\
		"vmovaps	%%ymm4,     (%%r12)	\n\t	vmovaps	%%ymm12,     (%%r12,%%r9)	\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm0	\n\t	vmovaps	0x040(%%r8),%%ymm8	 		\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm1	\n\t	vmovaps	0x060(%%r8),%%ymm9 			\n\t"\
		"vmovaps	     (%%rbx),%%ymm6	\n\t	vmovaps	0x080(%%rbx),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t	vmovaps	0x0a0(%%rbx),%%ymm15		\n\t"\
		"vmovaps	%%ymm6,%%ymm4		\n\t	vmovaps	%%ymm14,%%ymm12				\n\t"\
		"vmovaps	%%ymm7,%%ymm5		\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
	"vfnmadd231pd %%ymm1,%%ymm7,%%ymm4	\n\t vfnmadd231pd %%ymm9 ,%%ymm15,%%ymm12	\n\t"\
	"vfmadd231pd  %%ymm1,%%ymm6,%%ymm5	\n\t vfmadd231pd  %%ymm9 ,%%ymm14,%%ymm13	\n\t"\
		"vmovaps	     (%%r12),%%ymm0	\n\t	vmovaps	     (%%r12,%%r9),%%ymm8	\n\t"/* Restore 2 */\
		"vmovaps	0x020(%%r12),%%ymm1	\n\t	vmovaps	0x020(%%r12,%%r9),%%ymm9	\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	%%ymm13,%%ymm15				\n\t	 movq	%[__two],%%rsi	\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t	vmovaps	%%ymm12,%%ymm14				\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5	\n\t	vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		/* Finish radix-4 butterfly and store results: */\
		"vmovaps	     (%%r10),%%ymm0	\n\t	vmovaps	     (%%r10,%%r9),%%ymm8	\n\t"/* Restore 1 */\
		"vmovaps	0x020(%%r10),%%ymm1	\n\t	vmovaps	0x020(%%r10,%%r9),%%ymm9	\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1	\n\t	vsubpd	%%ymm15,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm0,     (%%r12)	\n\t	vmovaps	%%ymm8 ,     (%%r12,%%r9)	\n\t"	/* 2.0, shared by both columns ... moved +- until found best cycle count: */\
		"vmovaps	%%ymm2,     (%%r11)	\n\t	vmovaps	%%ymm10,     (%%r11,%%r9)	\n\t	vmovaps	(%%rsi),%%ymm0	\n\t"\
		"vmovaps	%%ymm1,0x020(%%r12)	\n\t	vmovaps	%%ymm9 ,0x020(%%r12,%%r9)	\n\t"\
		"vmovaps	%%ymm3,0x020(%%r13)	\n\t	vmovaps	%%ymm11,0x020(%%r13,%%r9)	\n\t"\
	"vfmadd213pd (%%r12),%%ymm0,%%ymm6	\n\t	vfmadd213pd	%%ymm8 ,%%ymm0,%%ymm14	\n\t"\
	"vfmadd213pd %%ymm2 ,%%ymm0,%%ymm5	\n\t	vfmadd213pd	%%ymm10,%%ymm0,%%ymm13	\n\t"\
	"vfmadd213pd %%ymm1 ,%%ymm0,%%ymm7	\n\t	vfmadd213pd	%%ymm9 ,%%ymm0,%%ymm15	\n\t"\
	"vfmadd213pd %%ymm3 ,%%ymm0,%%ymm4	\n\t	vfmadd213pd	%%ymm11,%%ymm0,%%ymm12	\n\t"\
		"vmovaps	%%ymm6,     (%%r10)	\n\t	vmovaps	%%ymm14,     (%%r10,%%r9)	\n\t"/* don't need reload-from-mem of ymm8/0x180(%%r12,%%r9) as we do in lcol, but 1 cycle faster with it. [!?] */\
		"vmovaps	%%ymm5,     (%%r13)	\n\t	vmovaps	%%ymm13,     (%%r13,%%r9)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%r10)	\n\t	vmovaps	%%ymm15,0x020(%%r10,%%r9)	\n\t"\
		"vmovaps	%%ymm4,0x020(%%r11)	\n\t	vmovaps	%%ymm12,0x020(%%r11,%%r9)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r0] "m" (Xr0)\
		 ,[__two] "m" (Xtwo)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	// The non-name-suffixed DIF,DIT macros in AVX2 mode are based on the AVX-mode _V2 macros:
	// 64 ADD, 112 FMA, 48 pure-MUL, 216 MOVAPS
	#define SSE2_RADIX16_DIT_TWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr0,Xtwo,Xcc0,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/*
		// Pass 1:
		j = p4*8;	// Set j to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k = 0x100;	// Set k to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		// 0x0c0 = literal-bytewise address offset between the pointers to the first and second DFT's twiddles:
		i0 = (vec_dbl *)a; i1 = (vec_dbl *)(a+p1); i2 = (vec_dbl *)(a+p2); i3 = (vec_dbl *)(a+p3);
		o0 = r00; o1 = r00+2; o2 = r00+4; o3 = r00+6;
		c_tmp = cc0+6;	// c2,1,3
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,0x0c0, o0,o1,o2,o3,k)
	*/\
		"movq	%[__two],%%rsi 				\n\t	movslq	%[__p4],%%rdi		\n\t"\
		"movq	%[__add0],%%rax				\n\t	movq	%[__r0],%%r10		\n\t"\
		"movslq	%[__p1],%%rbx				\n\t	leaq	(%%rax,%%rbx,8),%%rbx	\n\t"\
		"movslq	%[__p2],%%rcx				\n\t	leaq	(%%rax,%%rcx,8),%%rcx	\n\t"\
		"movslq	%[__p3],%%rdx				\n\t	leaq	(%%rax,%%rdx,8),%%rdx	\n\t"\
/* two: */"vmovaps (%%rsi),%%ymm15 \n\t shlq $3,%%rdi \n\t	leaq 0x020(%%rdx,%%rdi),%%r11 \n\t"\
		"vmovaps	     (%%rax),%%ymm0		\n\t	vmovaps	     (%%rax,%%rdi),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t	vmovaps	0x020(%%rax,%%rdi),%%ymm9 	\n\t"\
		"vmovaps	     (%%rbx),%%ymm2		\n\t	vmovaps	     (%%rbx,%%rdi),%%ymm10	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3		\n\t	vmovaps	0x020(%%rbx,%%rdi),%%ymm11	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4		\n\t	vmovaps	     (%%rcx,%%rdi),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5		\n\t	vmovaps	0x020(%%rcx,%%rdi),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6		\n\t	vmovaps	     (%%rdx,%%rdi),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7		\n\t"/*	vmovaps	0x020(%%rdx,%%rdi),%%ymm15	\n\t"*/\
		"vsubpd	%%ymm2 ,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm3 ,%%ymm1,%%ymm1		\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm6 ,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm7 ,%%ymm5,%%ymm5		\n\t	vsubpd	(%%r11),%%ymm13,%%ymm13		\n\t"\
	"vfmadd132pd %%ymm15,%%ymm0,%%ymm2		\n\t	vfmadd132pd	%%ymm15,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd %%ymm15,%%ymm1,%%ymm3		\n\t	vfmadd132pd	%%ymm15,%%ymm9 ,%%ymm11		\n\t"\
	"vfmadd132pd %%ymm15,%%ymm4,%%ymm6		\n\t	vfmadd132pd	%%ymm15,%%ymm12,%%ymm14		\n\t"\
	"vfmadd132pd %%ymm15,%%ymm5,%%ymm7		\n\t	vfmadd132pd	(%%r11),%%ymm13,%%ymm15		\n\t"\
		"vmovaps	%%ymm12,    (%%r10)		\n\t	vmovaps (%%rsi),%%ymm12			\n\t"/* spill ymm12 to free up reg for 2.0 */\
		"vsubpd	%%ymm6 ,%%ymm2,%%ymm2		\n\t	vsubpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7 ,%%ymm3,%%ymm3		\n\t	vsubpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm13,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm4 ,%%ymm1,%%ymm1		\n\t	vsubpd	(%%r10),%%ymm9 ,%%ymm9 		\n\t"\
	"vfmadd132pd %%ymm12,%%ymm2,%%ymm6		\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm14		\n\t"\
	"vfmadd132pd %%ymm12,%%ymm3,%%ymm7		\n\t	vfmadd132pd	%%ymm12,%%ymm11,%%ymm15		\n\t"\
	"vfmadd132pd %%ymm12,%%ymm0,%%ymm5		\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm13		\n\t"\
	"vfmadd132pd %%ymm12,%%ymm1,%%ymm4		\n\t	vfmadd132pd	(%%r10),%%ymm9 ,%%ymm12		\n\t"\
		"vmovaps	%%ymm6,    (%%r10)		\n\t	vmovaps	%%ymm7,0x020(%%r10)		\n\t"/* Write B0 to free up 2 regs */\
		/* B2 = t0*~w2 = t0*[c2-I.s1] */\
		"movq	%[__cc0],%%rsi				\n\t	vmovaps %%ymm14,0x100(%%r10)			\n\t"\
/* cc0+6: */"addq	$0x0c0,%%rsi			\n\t	vmovaps %%ymm15,0x120(%%r10)			\n\t"\
		"vmovaps	%%ymm2,%%ymm6			\n\t	vmovaps	%%ymm10,%%ymm14					\n\t"\
		"vmovaps	%%ymm3,%%ymm7			\n\t	vmovaps	%%ymm11,%%ymm15					\n\t"\
		"vmulpd	     (%%rsi),%%ymm6,%%ymm6	\n\t	vmulpd	0x0c0(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	     (%%rsi),%%ymm7,%%ymm7	\n\t	vmulpd	0x0c0(%%rsi),%%ymm15,%%ymm15	\n\t"\
	" vfmadd231pd 0x20(%%rsi),%%ymm3,%%ymm6	\n\t  vfmadd231pd 0xe0(%%rsi),%%ymm11,%%ymm14	\n\t"\
	"vfnmadd231pd 0x20(%%rsi),%%ymm2,%%ymm7	\n\t vfnmadd231pd 0xe0(%%rsi),%%ymm10,%%ymm15	\n\t"\
		"vmovaps	%%ymm6,0x080(%%r10)		\n\t	vmovaps %%ymm14,0x180(%%r10)			\n\t"\
		"vmovaps	%%ymm7,0x0a0(%%r10)		\n\t	vmovaps %%ymm15,0x1a0(%%r10)			\n\t"\
		/* B1 = t1*~w1 = t1*[c1-I.s1] */\
		"vmovaps	0x040(%%rsi),%%ymm2		\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm3		\n\t"\
		"vmovaps	%%ymm5,%%ymm6			\n\t	vmovaps	%%ymm13,%%ymm14					\n\t"\
		"vmovaps	%%ymm1,%%ymm7			\n\t	vmovaps	%%ymm9 ,%%ymm15					\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6		\n\t	vmulpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7		\n\t	vmulpd	0x100(%%rsi),%%ymm15,%%ymm15	\n\t"\
	" vfmadd231pd	%%ymm3,%%ymm1,%%ymm6	\n\t  vfmadd231pd 0x120(%%rsi),%%ymm9 ,%%ymm14	\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm5,%%ymm7	\n\t vfnmadd231pd 0x120(%%rsi),%%ymm13,%%ymm15	\n\t"\
		"vmovaps	0x080(%%rsi),%%ymm2		\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm3		\n\t"\
		"vmovaps	%%ymm6,0x040(%%r10)		\n\t	vmovaps %%ymm14,0x140(%%r10)			\n\t"\
		"vmovaps	%%ymm7,0x060(%%r10)		\n\t	vmovaps %%ymm15,0x160(%%r10)			\n\t"\
		/* B3 = t3*~w3 = t3*[c3-I.s3] */\
		"vmovaps	%%ymm0,%%ymm6			\n\t	vmovaps	%%ymm8 ,%%ymm14					\n\t"\
		"vmovaps	%%ymm4,%%ymm7			\n\t	vmovaps	%%ymm12,%%ymm15					\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6		\n\t	vmulpd	0x140(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7		\n\t	vmulpd	0x140(%%rsi),%%ymm15,%%ymm15	\n\t"\
	" vfmadd231pd	%%ymm3,%%ymm4,%%ymm6	\n\t  vfmadd231pd 0x160(%%rsi),%%ymm12,%%ymm14	\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm0,%%ymm7	\n\t vfnmadd231pd 0x160(%%rsi),%%ymm8 ,%%ymm15	\n\t"\
		"vmovaps	%%ymm6,0x0c0(%%r10)		\n\t	vmovaps	%%ymm14,0x1c0(%%r10)			\n\t"\
		"vmovaps	%%ymm7,0x0e0(%%r10)		\n\t	vmovaps	%%ymm15,0x1e0(%%r10)			\n\t"\
	/*
		i0 = (vec_dbl *)(a+p8); i1 = (vec_dbl *)(a+p9); i2 = (vec_dbl *)(a+pA); i3 = (vec_dbl *)(a+pB);
		o0 += 16; o1 += 16; o2 += 16; o3 += 16;
		c_tmp += 12;		// cA,9,B
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,0x0c0, o0,o1,o2,o3,k)
	*/\
		"movq	%[__two],%%rsi 		\n\t	movslq	%[__p8],%%r11		\n\t"\
		"leaq (%%rax,%%r11,8),%%rax	\n\t	addq	 $0x200,%%r10		\n\t"/* r0 + 16 */\
		"leaq (%%rbx,%%r11,8),%%rbx	\n\t"\
		"leaq (%%rcx,%%r11,8),%%rcx	\n\t"\
		"leaq (%%rdx,%%r11,8),%%rdx	\n\t"\
	"vmovaps	(%%rsi),%%ymm15		\n\t"/* two */"	leaq	0x020(%%rdx,%%rdi),%%r11	\n\t"\
		"vmovaps	     (%%rax),%%ymm0		\n\t	vmovaps	     (%%rax,%%rdi),%%ymm8 	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t	vmovaps	0x020(%%rax,%%rdi),%%ymm9 	\n\t"\
		"vmovaps	     (%%rbx),%%ymm2		\n\t	vmovaps	     (%%rbx,%%rdi),%%ymm10	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3		\n\t	vmovaps	0x020(%%rbx,%%rdi),%%ymm11	\n\t"\
		"vmovaps	     (%%rcx),%%ymm4		\n\t	vmovaps	     (%%rcx,%%rdi),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5		\n\t	vmovaps	0x020(%%rcx,%%rdi),%%ymm13	\n\t"\
		"vmovaps	     (%%rdx),%%ymm6		\n\t	vmovaps	     (%%rdx,%%rdi),%%ymm14	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7		\n\t"/*	vmovaps	0x020(%%rdx,%%rdi),%%ymm15	\n\t"*/\
		"vsubpd	%%ymm2 ,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm3 ,%%ymm1,%%ymm1		\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm6 ,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm7 ,%%ymm5,%%ymm5		\n\t	vsubpd	(%%r11),%%ymm13,%%ymm13		\n\t"\
	"vfmadd132pd %%ymm15,%%ymm0,%%ymm2		\n\t	vfmadd132pd	%%ymm15,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd %%ymm15,%%ymm1,%%ymm3		\n\t	vfmadd132pd	%%ymm15,%%ymm9 ,%%ymm11		\n\t"\
	"vfmadd132pd %%ymm15,%%ymm4,%%ymm6		\n\t	vfmadd132pd	%%ymm15,%%ymm12,%%ymm14		\n\t"\
	"vfmadd132pd %%ymm15,%%ymm5,%%ymm7		\n\t	vfmadd132pd	(%%r11),%%ymm13,%%ymm15		\n\t"\
		"vmovaps	%%ymm12,    (%%r10)		\n\t	vmovaps (%%rsi),%%ymm12			\n\t"/* spill ymm12 to free up reg for 2.0 */\
		"vsubpd	%%ymm6 ,%%ymm2,%%ymm2		\n\t	vsubpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7 ,%%ymm3,%%ymm3		\n\t	vsubpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm13,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm4 ,%%ymm1,%%ymm1		\n\t	vsubpd	(%%r10),%%ymm9 ,%%ymm9 		\n\t"\
	"vfmadd132pd %%ymm12,%%ymm2,%%ymm6		\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm14		\n\t"\
	"vfmadd132pd %%ymm12,%%ymm3,%%ymm7		\n\t	vfmadd132pd	%%ymm12,%%ymm11,%%ymm15		\n\t"\
	"vfmadd132pd %%ymm12,%%ymm0,%%ymm5		\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm13		\n\t"\
	"vfmadd132pd %%ymm12,%%ymm1,%%ymm4		\n\t	vfmadd132pd	(%%r10),%%ymm9 ,%%ymm12		\n\t"\
		"vmovaps	%%ymm6,    (%%r10)		\n\t	vmovaps	%%ymm7,0x020(%%r10)		\n\t"/* Write B0 to free up 2 regs */\
		/* B2 = t0*~w2 = t0*[c2-I.s1] */\
		"movq	%[__cc0],%%rsi				\n\t	vmovaps %%ymm14,0x100(%%r10)			\n\t"\
/* cc0+18: */"addq	$0x240,%%rsi			\n\t	vmovaps %%ymm15,0x120(%%r10)			\n\t"\
		"vmovaps	%%ymm2,%%ymm6			\n\t	vmovaps	%%ymm10,%%ymm14					\n\t"\
		"vmovaps	%%ymm3,%%ymm7			\n\t	vmovaps	%%ymm11,%%ymm15					\n\t"\
		"vmulpd	     (%%rsi),%%ymm6,%%ymm6	\n\t	vmulpd	0x0c0(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	     (%%rsi),%%ymm7,%%ymm7	\n\t	vmulpd	0x0c0(%%rsi),%%ymm15,%%ymm15	\n\t"\
	" vfmadd231pd 0x20(%%rsi),%%ymm3,%%ymm6	\n\t  vfmadd231pd 0xe0(%%rsi),%%ymm11,%%ymm14	\n\t"\
	"vfnmadd231pd 0x20(%%rsi),%%ymm2,%%ymm7	\n\t vfnmadd231pd 0xe0(%%rsi),%%ymm10,%%ymm15	\n\t"\
		"vmovaps	%%ymm6,0x080(%%r10)		\n\t	vmovaps %%ymm14,0x180(%%r10)			\n\t"\
		"vmovaps	%%ymm7,0x0a0(%%r10)		\n\t	vmovaps %%ymm15,0x1a0(%%r10)			\n\t"\
		/* B1 = t1*~w1 = t1*[c1-I.s1] */\
		"vmovaps	0x040(%%rsi),%%ymm2		\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm3		\n\t"\
		"vmovaps	%%ymm5,%%ymm6			\n\t	vmovaps	%%ymm13,%%ymm14					\n\t"\
		"vmovaps	%%ymm1,%%ymm7			\n\t	vmovaps	%%ymm9 ,%%ymm15					\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6		\n\t	vmulpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7		\n\t	vmulpd	0x100(%%rsi),%%ymm15,%%ymm15	\n\t"\
	" vfmadd231pd	%%ymm3,%%ymm1,%%ymm6	\n\t  vfmadd231pd 0x120(%%rsi),%%ymm9 ,%%ymm14	\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm5,%%ymm7	\n\t vfnmadd231pd 0x120(%%rsi),%%ymm13,%%ymm15	\n\t"\
		"vmovaps	0x080(%%rsi),%%ymm2		\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm3		\n\t"\
		"vmovaps	%%ymm6,0x040(%%r10)		\n\t	vmovaps %%ymm14,0x140(%%r10)			\n\t"\
		"vmovaps	%%ymm7,0x060(%%r10)		\n\t	vmovaps %%ymm15,0x160(%%r10)			\n\t"\
		/* B3 = t3*~w3 = t3*[c3-I.s3] */\
		"vmovaps	%%ymm0,%%ymm6			\n\t	vmovaps	%%ymm8 ,%%ymm14					\n\t"\
		"vmovaps	%%ymm4,%%ymm7			\n\t	vmovaps	%%ymm12,%%ymm15					\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6		\n\t	vmulpd	0x140(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7		\n\t	vmulpd	0x140(%%rsi),%%ymm15,%%ymm15	\n\t"\
	" vfmadd231pd	%%ymm3,%%ymm4,%%ymm6	\n\t  vfmadd231pd 0x160(%%rsi),%%ymm12,%%ymm14	\n\t"\
	"vfnmadd231pd	%%ymm3,%%ymm0,%%ymm7	\n\t vfnmadd231pd 0x160(%%rsi),%%ymm8 ,%%ymm15	\n\t"\
		"vmovaps	%%ymm6,0x0c0(%%r10)		\n\t	vmovaps	%%ymm14,0x1c0(%%r10)			\n\t"\
		"vmovaps	%%ymm7,0x0e0(%%r10)		\n\t	vmovaps	%%ymm15,0x1e0(%%r10)			\n\t"\
	/*
		// Pass 2:
		j = 0x080;	// Set j to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k = p2*8;	// Set k to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		c_tmp = cc0;	// c8,4,C
		// Pass1-DFTs all use same twiddle-triplet, so lit-byte address offset between ptrs to 1st,2nd DFT's twiddles = 0
		i0 = r00; i1 = r00+8; i2 = r00+16; i3 = r00+24;
		o0 = (vec_dbl *)a; o1 = (vec_dbl *)(a+p4); o2 = (vec_dbl *)(a+p8); o3 = (vec_dbl *)(a+pC);
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,   0, o0,o1,o2,o3,k)
	*/\
		"movq	%[__two],%%rsi 				\n\t	movslq	%[__p2],%%r9		\n\t"\
		"movq	%[__add0],%%r10				\n\t	movq	%[__r0],%%rax		\n\t"\
		"movslq	%[__p4]  ,%%r11				\n\t	leaq	(%%r10,%%r11,8),%%r11	\n\t"\
		"movslq	%[__p8]  ,%%r12				\n\t	leaq	(%%r10,%%r12,8),%%r12	\n\t"\
		"movslq	%[__p12] ,%%r13				\n\t	leaq	(%%r10,%%r13,8),%%r13	\n\t"\
		"vmovaps	(%%rsi),%%ymm15	\n\t"/* two */"	leaq	0x080(%%rax),%%rbx	\n\t"\
		"vmovaps	     (%%rax),%%ymm0		\n\t	vmovaps	     (%%rbx),%%ymm8 	\n\t"/* r0 */\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t	vmovaps	0x020(%%rbx),%%ymm9 	\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2		\n\t	vmovaps	0x100(%%rbx),%%ymm10	\n\t"/* r0 + 8 */\
		"vmovaps	0x120(%%rax),%%ymm3		\n\t	vmovaps	0x120(%%rbx),%%ymm11	\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4		\n\t	vmovaps	0x200(%%rbx),%%ymm12	\n\t"/* r0 + 16 */\
		"vmovaps	0x220(%%rax),%%ymm5		\n\t	vmovaps	0x220(%%rbx),%%ymm13	\n\t"\
		"vmovaps	0x300(%%rax),%%ymm6		\n\t	vmovaps	0x300(%%rbx),%%ymm14	\n\t"/* r0 + 24 */\
		"vmovaps	0x320(%%rax),%%ymm7		\n\t"/*	vmovaps	0x320(%%rbx),%%ymm15	\n\t*/"	leaq	0x320(%%rbx),%%r8	\n\t"\
		"vsubpd	%%ymm2 ,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm3 ,%%ymm1,%%ymm1		\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm6 ,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm7 ,%%ymm5,%%ymm5		\n\t	vsubpd	(%%r8) ,%%ymm13,%%ymm13		\n\t"\
	"vfmadd132pd %%ymm15,%%ymm0,%%ymm2		\n\t	vfmadd132pd	%%ymm15,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd %%ymm15,%%ymm1,%%ymm3		\n\t	vfmadd132pd	%%ymm15,%%ymm9 ,%%ymm11		\n\t"\
	"vfmadd132pd %%ymm15,%%ymm4,%%ymm6		\n\t	vfmadd132pd	%%ymm15,%%ymm12,%%ymm14		\n\t"\
	"vfmadd132pd %%ymm15,%%ymm5,%%ymm7		\n\t	vfmadd132pd	(%%r8),%%ymm13,%%ymm15		\n\t"\
		"vmovaps	%%ymm12,    (%%rax)		\n\t	vmovaps (%%rsi),%%ymm12			\n\t"/* spill ymm12 to free up reg for 2.0 */\
		"vsubpd	%%ymm6 ,%%ymm2,%%ymm2		\n\t	vsubpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7 ,%%ymm3,%%ymm3		\n\t	vsubpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm13,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm4 ,%%ymm1,%%ymm1		\n\t	vsubpd	(%%rax),%%ymm9 ,%%ymm9 		\n\t"\
	"vfmadd132pd %%ymm12,%%ymm2,%%ymm6		\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm14		\n\t"\
	"vfmadd132pd %%ymm12,%%ymm3,%%ymm7		\n\t	vfmadd132pd	%%ymm12,%%ymm11,%%ymm15		\n\t"\
	"vfmadd132pd %%ymm12,%%ymm0,%%ymm5		\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm13		\n\t"\
	"vfmadd132pd %%ymm12,%%ymm1,%%ymm4		\n\t	vfmadd132pd	(%%rax),%%ymm9 ,%%ymm12		\n\t"\
		"vmovaps	%%ymm6,    (%%r10)		\n\t	vmovaps	%%ymm7,0x020(%%r10)	\n\t"/* Write B0 to free up 2 regs */\
		/* B2 = t0*~w2 = t0*[c2-I.s1] */"			shlq	$3,%%r9				\n\t"/* p2*8: */\
		"vmovaps	%%ymm0,    (%%rax)		\n\t	vmovaps	%%ymm4,0x020(%%rax)	\n\t"/* Write mm0,4 to free up 2 regs */\
		"movq	%[__cc0],%%rsi			 	\n\t	vmovaps %%ymm14,     (%%r10,%%r9)		\n\t"\
		"											vmovaps %%ymm15,0x020(%%r10,%%r9)		\n\t"\
		"vmovaps	     (%%rsi),%%ymm0		\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm4		\n\t"\
		"vmovaps	%%ymm2,%%ymm6			\n\t	vmovaps	%%ymm10,%%ymm14					\n\t"\
		"vmovaps	%%ymm3,%%ymm7			\n\t	vmovaps	%%ymm11,%%ymm15					\n\t"\
		"vmulpd	  %%ymm0,%%ymm6,%%ymm6		\n\t	vmulpd	  %%ymm0,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	  %%ymm0,%%ymm7,%%ymm7		\n\t	vmulpd	  %%ymm0,%%ymm15,%%ymm15	\n\t"\
	" vfmadd231pd %%ymm4,%%ymm3,%%ymm6		\n\t  vfmadd231pd %%ymm4,%%ymm11,%%ymm14	\n\t"\
	"vfnmadd231pd %%ymm4,%%ymm2,%%ymm7		\n\t vfnmadd231pd %%ymm4,%%ymm10,%%ymm15	\n\t"\
		"vmovaps	%%ymm6,     (%%r12)		\n\t	vmovaps %%ymm14,     (%%r12,%%r9)		\n\t"\
		"vmovaps	%%ymm7,0x020(%%r12)		\n\t	vmovaps %%ymm15,0x020(%%r12,%%r9)		\n\t"\
		"vmovaps	    (%%rax),%%ymm0		\n\t	vmovaps	0x020(%%rax),%%ymm4	\n\t"/* Restore spill of mm0,4 */\
		/* B1 = t1*~w1 = t1*[c1-I.s1] */\
		"vmovaps	0x040(%%rsi),%%ymm2		\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm3		\n\t"\
		"vmovaps	     %%ymm5 ,%%ymm6		\n\t	vmovaps	%%ymm13,%%ymm14				\n\t"\
		"vmovaps	     %%ymm1 ,%%ymm7		\n\t	vmovaps	%%ymm9 ,%%ymm15				\n\t"\
		"vmulpd	  %%ymm2,%%ymm6,%%ymm6		\n\t	vmulpd	  %%ymm2,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	  %%ymm2,%%ymm7,%%ymm7		\n\t	vmulpd	  %%ymm2,%%ymm15,%%ymm15	\n\t"\
	" vfmadd231pd %%ymm3,%%ymm1,%%ymm6		\n\t  vfmadd231pd %%ymm3,%%ymm9 ,%%ymm14	\n\t"\
	"vfnmadd231pd %%ymm3,%%ymm5,%%ymm7		\n\t vfnmadd231pd %%ymm3,%%ymm13,%%ymm15	\n\t"\
		"vmovaps	%%ymm6,     (%%r11)		\n\t	vmovaps %%ymm14,     (%%r11,%%r9)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%r11)		\n\t	vmovaps %%ymm15,0x020(%%r11,%%r9)	\n\t"\
		/* B3 = t3*~w3 = t3*[c3-I.s3] */\
		"vmovaps	0x080(%%rsi),%%ymm2		\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm3		\n\t"\
		"vmovaps	%%ymm0,%%ymm6			\n\t	vmovaps	%%ymm8 ,%%ymm14				\n\t"\
		"vmovaps	%%ymm4,%%ymm7			\n\t	vmovaps	%%ymm12,%%ymm15				\n\t"\
		"vmulpd	  %%ymm2,%%ymm6,%%ymm6		\n\t	vmulpd	  %%ymm2,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	  %%ymm2,%%ymm7,%%ymm7		\n\t	vmulpd	  %%ymm2,%%ymm15,%%ymm15	\n\t"\
	" vfmadd231pd %%ymm3,%%ymm4,%%ymm6		\n\t  vfmadd231pd %%ymm3,%%ymm12,%%ymm14	\n\t"\
	"vfnmadd231pd %%ymm3,%%ymm0,%%ymm7		\n\t vfnmadd231pd %%ymm3,%%ymm8 ,%%ymm15	\n\t"\
		"vmovaps	%%ymm6,     (%%r13)		\n\t	vmovaps	%%ymm14,     (%%r13,%%r9)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%r13)		\n\t	vmovaps	%%ymm15,0x020(%%r13,%%r9)	\n\t"\
	/*
		i0 += 2; i1 += 2; i2 += 2; i3 += 2;
		o0 = (vec_dbl *)(a+p1); o1 = (vec_dbl *)(a+p5); o2 = (vec_dbl *)(a+p9); o3 = (vec_dbl *)(a+pD);
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,   0, o0,o1,o2,o3,k)
	*/\
		"movq	%[__two],%%rsi 		\n\t	movslq	%[__p1],%%rdi		\n\t"\
		"leaq (%%r10,%%rdi,8),%%r10	\n\t	addq	$0x040,%%rax		\n\t"\
		"leaq (%%r11,%%rdi,8),%%r11	\n\t"\
		"leaq (%%r12,%%rdi,8),%%r12	\n\t"\
		"leaq (%%r13,%%rdi,8),%%r13	\n\t"\
		"vmovaps	(%%rsi),%%ymm15	\n\t"/* two */"	leaq	0x080(%%rax),%%rbx	\n\t"\
		"vmovaps	     (%%rax),%%ymm0		\n\t	vmovaps	     (%%rbx),%%ymm8 	\n\t"/* r0 */\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t	vmovaps	0x020(%%rbx),%%ymm9 	\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2		\n\t	vmovaps	0x100(%%rbx),%%ymm10	\n\t"/* r0 + 8 */\
		"vmovaps	0x120(%%rax),%%ymm3		\n\t	vmovaps	0x120(%%rbx),%%ymm11	\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4		\n\t	vmovaps	0x200(%%rbx),%%ymm12	\n\t"/* r0 + 16 */\
		"vmovaps	0x220(%%rax),%%ymm5		\n\t	vmovaps	0x220(%%rbx),%%ymm13	\n\t"\
		"vmovaps	0x300(%%rax),%%ymm6		\n\t	vmovaps	0x300(%%rbx),%%ymm14	\n\t"/* r0 + 24 */\
		"vmovaps	0x320(%%rax),%%ymm7		\n\t"/*	vmovaps	0x320(%%rbx),%%ymm15	\n\t*/"	leaq	0x320(%%rbx),%%r8	\n\t"\
		"vsubpd	%%ymm2 ,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm3 ,%%ymm1,%%ymm1		\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm6 ,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm7 ,%%ymm5,%%ymm5		\n\t	vsubpd	(%%r8) ,%%ymm13,%%ymm13		\n\t"\
	"vfmadd132pd %%ymm15,%%ymm0,%%ymm2		\n\t	vfmadd132pd	%%ymm15,%%ymm8 ,%%ymm10		\n\t"\
	"vfmadd132pd %%ymm15,%%ymm1,%%ymm3		\n\t	vfmadd132pd	%%ymm15,%%ymm9 ,%%ymm11		\n\t"\
	"vfmadd132pd %%ymm15,%%ymm4,%%ymm6		\n\t	vfmadd132pd	%%ymm15,%%ymm12,%%ymm14		\n\t"\
	"vfmadd132pd %%ymm15,%%ymm5,%%ymm7		\n\t	vfmadd132pd	(%%r8),%%ymm13,%%ymm15		\n\t"\
		"vmovaps	%%ymm12,    (%%rax)		\n\t	vmovaps (%%rsi),%%ymm12			\n\t"/* spill ymm12 to free up reg for 2.0 */\
		"vsubpd	%%ymm6 ,%%ymm2,%%ymm2		\n\t	vsubpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7 ,%%ymm3,%%ymm3		\n\t	vsubpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm13,%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm4 ,%%ymm1,%%ymm1		\n\t	vsubpd	(%%rax),%%ymm9 ,%%ymm9 		\n\t"\
	"vfmadd132pd %%ymm12,%%ymm2,%%ymm6		\n\t	vfmadd132pd	%%ymm12,%%ymm10,%%ymm14		\n\t"\
	"vfmadd132pd %%ymm12,%%ymm3,%%ymm7		\n\t	vfmadd132pd	%%ymm12,%%ymm11,%%ymm15		\n\t"\
	"vfmadd132pd %%ymm12,%%ymm0,%%ymm5		\n\t	vfmadd132pd	%%ymm12,%%ymm8 ,%%ymm13		\n\t"\
	"vfmadd132pd %%ymm12,%%ymm1,%%ymm4		\n\t	vfmadd132pd	(%%rax),%%ymm9 ,%%ymm12		\n\t"\
		"vmovaps	%%ymm6,    (%%r10)		\n\t	vmovaps	%%ymm7,0x020(%%r10)	\n\t"/* Write B0 to free up 2 regs */\
		/* B2 = t0*~w2 = t0*[c2-I.s1] */\
		"vmovaps	%%ymm0,    (%%rax)		\n\t	vmovaps	%%ymm4,0x020(%%rax)	\n\t"/* Write mm0,4 to free up 2 regs */\
		"movq	%[__cc0],%%rsi			 	\n\t	vmovaps %%ymm14,     (%%r10,%%r9)		\n\t"\
		"											vmovaps %%ymm15,0x020(%%r10,%%r9)		\n\t"\
		"vmovaps	     (%%rsi),%%ymm0		\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm4		\n\t"\
		"vmovaps	%%ymm2,%%ymm6			\n\t	vmovaps	%%ymm10,%%ymm14					\n\t"\
		"vmovaps	%%ymm3,%%ymm7			\n\t	vmovaps	%%ymm11,%%ymm15					\n\t"\
		"vmulpd	  %%ymm0,%%ymm6,%%ymm6		\n\t	vmulpd	  %%ymm0,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	  %%ymm0,%%ymm7,%%ymm7		\n\t	vmulpd	  %%ymm0,%%ymm15,%%ymm15	\n\t"\
	" vfmadd231pd %%ymm4,%%ymm3,%%ymm6		\n\t  vfmadd231pd %%ymm4,%%ymm11,%%ymm14	\n\t"\
	"vfnmadd231pd %%ymm4,%%ymm2,%%ymm7		\n\t vfnmadd231pd %%ymm4,%%ymm10,%%ymm15	\n\t"\
		"vmovaps	%%ymm6,     (%%r12)		\n\t	vmovaps %%ymm14,     (%%r12,%%r9)		\n\t"\
		"vmovaps	%%ymm7,0x020(%%r12)		\n\t	vmovaps %%ymm15,0x020(%%r12,%%r9)		\n\t"\
		"vmovaps	    (%%rax),%%ymm0		\n\t	vmovaps	0x020(%%rax),%%ymm4	\n\t"/* Restore spill of mm0,4 */\
		/* B1 = t1*~w1 = t1*[c1-I.s1] */\
		"vmovaps	0x040(%%rsi),%%ymm2		\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm3		\n\t"\
		"vmovaps	     %%ymm5 ,%%ymm6		\n\t	vmovaps	%%ymm13,%%ymm14				\n\t"\
		"vmovaps	     %%ymm1 ,%%ymm7		\n\t	vmovaps	%%ymm9 ,%%ymm15				\n\t"\
		"vmulpd	  %%ymm2,%%ymm6,%%ymm6		\n\t	vmulpd	  %%ymm2,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	  %%ymm2,%%ymm7,%%ymm7		\n\t	vmulpd	  %%ymm2,%%ymm15,%%ymm15	\n\t"\
	" vfmadd231pd %%ymm3,%%ymm1,%%ymm6		\n\t  vfmadd231pd %%ymm3,%%ymm9 ,%%ymm14	\n\t"\
	"vfnmadd231pd %%ymm3,%%ymm5,%%ymm7		\n\t vfnmadd231pd %%ymm3,%%ymm13,%%ymm15	\n\t"\
		"vmovaps	%%ymm6,     (%%r11)		\n\t	vmovaps %%ymm14,     (%%r11,%%r9)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%r11)		\n\t	vmovaps %%ymm15,0x020(%%r11,%%r9)	\n\t"\
		/* B3 = t3*~w3 = t3*[c3-I.s3] */\
		"vmovaps	0x080(%%rsi),%%ymm2		\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm3		\n\t"\
		"vmovaps	%%ymm0,%%ymm6			\n\t	vmovaps	%%ymm8 ,%%ymm14				\n\t"\
		"vmovaps	%%ymm4,%%ymm7			\n\t	vmovaps	%%ymm12,%%ymm15				\n\t"\
		"vmulpd	  %%ymm2,%%ymm6,%%ymm6		\n\t	vmulpd	  %%ymm2,%%ymm14,%%ymm14	\n\t"\
		"vmulpd	  %%ymm2,%%ymm7,%%ymm7		\n\t	vmulpd	  %%ymm2,%%ymm15,%%ymm15	\n\t"\
	" vfmadd231pd %%ymm3,%%ymm4,%%ymm6		\n\t  vfmadd231pd %%ymm3,%%ymm12,%%ymm14	\n\t"\
	"vfnmadd231pd %%ymm3,%%ymm0,%%ymm7		\n\t vfnmadd231pd %%ymm3,%%ymm8 ,%%ymm15	\n\t"\
		"vmovaps	%%ymm6,     (%%r13)		\n\t	vmovaps	%%ymm14,     (%%r13,%%r9)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%r13)		\n\t	vmovaps	%%ymm15,0x020(%%r13,%%r9)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r0] "m" (Xr0)\
		 ,[__two] "m" (Xtwo)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	// Initial version of trig-data and DFT-pass routines uses 4-copied trig data, as usual.
	// My timing tests indicate that using the more-intricate and much-more-instruction-heavy
	// Newtonian iterative inversion is appreciably faster than using vdivpd: ~20% faster than
	// using vdivpd to compute the same 4 reciprocals and shuffling the 1st of those prior to
	// multiplying to get the 3 needed cosine ratios, and a full 2x faster than shuffling the
	// first cosine-vector on input and doing 7 vdivpd [using 2 of the input cosine vectors twice,
	// and the shuffled-input-cosine once] directly of the 4 sine and 3 cosine-divisor vectors.
	//
	// For intermediate DFT passes, where we re-use each set of the trig data many times, this speed
	// difference will likely be negligible, but for the wrapper-square step where we use each set
	// of trigs just once, every cycle counts.
	//
	#define RADIX16_COMPUTE_FMA_SINCOS_DIF(Xadd0,Xone,Xcos,Xtan)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax			\n\t"/* base-ptr to 32-byte-aligned 128-byte-large memchunk 16 doubles needing iterative inversion */\
		"movq	%[__one] ,%%rbx			\n\t"/* 1.0 in 4-fold-double form */\
		"vmovaps	    (%%rbx),%%ymm14	\n\t"/* 1.0 */\
		"vmovaps	    (%%rax),%%ymm4	\n\t"/* load inputs-to-be-inverted (d's) into 4 AVX registers */\
		"vmovaps	0x20(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x40(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x60(%%rax),%%ymm7	\n\t"\
		"vaddpd	%%ymm14,%%ymm14,%%ymm14	\n\t"/* 2.0 */\
		"vcvtpd2ps	%%ymm4,%%xmm0	\n\t"/* convert d's to SP ... Note in AVX mode output register *must* be a 128-bit xmm! */\
		"vcvtpd2ps	%%ymm5,%%xmm1	\n\t"\
		"vcvtpd2ps	%%ymm6,%%xmm2	\n\t"\
		"vcvtpd2ps	%%ymm7,%%xmm3	\n\t"\
		"vrcpps		%%xmm0,%%xmm0	\n\t"	/* ainv := approx 1/d to 11-12 bits of precision */\
		"vrcpps		%%xmm1,%%xmm1	\n\t"\
		"vrcpps		%%xmm2,%%xmm2	\n\t"\
		"vrcpps		%%xmm3,%%xmm3	\n\t"\
		"vcvtps2pd	%%xmm0,%%ymm0	\n\t"	/* convert ~1/d back to DP  ... Note in AVX mode input register *must* be a 128-bit xmm!*/\
		"vcvtps2pd	%%xmm1,%%ymm1	\n\t"\
		"vcvtps2pd	%%xmm2,%%ymm2	\n\t"\
		"vcvtps2pd	%%xmm3,%%ymm3	\n\t"\
		/* 1st NR iteration gives ~23 bits of precision: */\
		"vmovaps	%%ymm0,%%ymm8 	\n\t"	/* make a copy of ainv */\
		"vmovaps	%%ymm1,%%ymm9 	\n\t"\
		"vmovaps	%%ymm2,%%ymm10	\n\t"\
		"vmovaps	%%ymm3,%%ymm11	\n\t"\
		"vfnmadd132pd	%%ymm4,%%ymm14,%%ymm0 	\n\t"	/* 2 - d*ainv, overwrites ainv */\
		"vfnmadd132pd	%%ymm5,%%ymm14,%%ymm1 	\n\t"\
		"vfnmadd132pd	%%ymm6,%%ymm14,%%ymm2 	\n\t"\
		"vfnmadd132pd	%%ymm7,%%ymm14,%%ymm3 	\n\t"\
		"vmulpd		%%ymm0,%%ymm8 ,%%ymm0	\n\t"	/* ainv*(2 - d*ainv) = 1/d accurate to ~23 bits */\
		"vmulpd		%%ymm1,%%ymm9 ,%%ymm1	\n\t"\
		"vmulpd		%%ymm2,%%ymm10,%%ymm2	\n\t"\
		"vmulpd		%%ymm3,%%ymm11,%%ymm3	\n\t"\
		/* 3rd-order update of 23-bit result needs just 2 FMA, 1 SUB, 1 MUL: */\
		"vaddpd	    (%%rbx),%%ymm14,%%ymm14	\n\t"/* 3.0, needed by 3rd-order update step */\
		"vmovaps	%%ymm0,%%ymm8 	\n\t"	/* make a copy of ainv */\
		"vmovaps	%%ymm1,%%ymm9 	\n\t"\
		"vmovaps	%%ymm2,%%ymm10	\n\t"\
		"vmovaps	%%ymm3,%%ymm11	\n\t"\
		"vfnmadd132pd	%%ymm0,%%ymm14,%%ymm4 	\n\t"/* 1st FMA overwrites d data (inputs) with (3 - d*ainv) */\
		"vfnmadd132pd	%%ymm1,%%ymm14,%%ymm5 	\n\t"\
		"vfnmadd132pd	%%ymm2,%%ymm14,%%ymm6 	\n\t"\
		"vfnmadd132pd	%%ymm3,%%ymm14,%%ymm7 	\n\t"\
		"vsubpd		%%ymm14,%%ymm4,%%ymm0 	\n\t"/* Subtract 3 from (3 - d*ainv) to get -y = -d*ainv terms in ymm0-3 */\
		"vsubpd		%%ymm14,%%ymm5,%%ymm1 	\n\t"\
		"vsubpd		%%ymm14,%%ymm6,%%ymm2 	\n\t"\
		"vsubpd		%%ymm14,%%ymm7,%%ymm3 	\n\t"\
		"vfmadd132pd	%%ymm4,%%ymm14,%%ymm0 	\n\t"/* Positive-product FMA gives (3 - y*(3 - d*ainv)) in ymm0-3*/\
		"vfmadd132pd	%%ymm5,%%ymm14,%%ymm1 	\n\t"\
		"vfmadd132pd	%%ymm6,%%ymm14,%%ymm2 	\n\t"\
		"vfmadd132pd	%%ymm7,%%ymm14,%%ymm3 	\n\t"\
		"vmulpd		%%ymm0,%%ymm8 ,%%ymm0	\n\t"	/* ainv*(3 - y*(3 - d*ainv)) = 1/d accurate to ~53 bits */\
		"vmulpd		%%ymm1,%%ymm9 ,%%ymm1	\n\t"\
		"vmulpd		%%ymm2,%%ymm10,%%ymm2	\n\t"\
		"vmulpd		%%ymm3,%%ymm11,%%ymm3	\n\t"\
		/* Multiply by the sine terms to get quotients in ymm4-7 (still need inverses in ymm0-3 for cosine-ratios): */\
		"vmulpd	0x80(%%rax),%%ymm0,%%ymm4	\n\t"\
		"vmulpd	0xa0(%%rax),%%ymm1,%%ymm5	\n\t"\
		"vmulpd	0xc0(%%rax),%%ymm2,%%ymm6	\n\t"\
		"vmulpd	0xe0(%%rax),%%ymm3,%%ymm7	\n\t"\
		/* Outputs overwrite sines: */\
		"vmovaps	%%ymm4,0x80(%%rax)	\n\t"\
		"vmovaps	%%ymm5,0xa0(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0xc0(%%rax)	\n\t"\
		"vmovaps	%%ymm7,0xe0(%%rax)	\n\t"\
		/* Set up for cosine ratios: */\
		"vpermilpd	$11,%%ymm0,%%ymm2 	\n\t"/* permute ymm0 = 1/c[3123] to get ymm2 = 1/c[1123], then *= [c3,c5,c6,0] */\
		"vmulpd	0x100(%%rax),%%ymm2,%%ymm2	\n\t"/* ymm2 = [c3,c5,c6, 0] * 1/c[1123] = [c31,c51,c62,  0] */\
		"vmulpd	0x120(%%rax),%%ymm0,%%ymm0	\n\t"/* ymm0 = [c7,c9,cA,cB] * 1/c[3123] = [c73,c91,cA2,cB3] */\
		"vmulpd	0x140(%%rax),%%ymm1,%%ymm1	\n\t"/* ymm1 = [cC,cD,cE,cF] * 1/c[4567] = [cC4,cD5,cE6,cF7] */\
		/* Outputs into slots imm. above above tangents: */\
		"vmovaps	%%ymm2,0x100(%%rax)	\n\t"\
		"vmovaps	%%ymm0,0x120(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x140(%%rax)	\n\t"\
		/* Now propagate the scalar doubles in the 11 ymm-sized memlocs referenced above to their final 4x-copied locs. Scalar data laid out as */\
		/* add0 + 0x[  0,  8, 10, 18]: [c3 ,c1 ,c2 ,c3 ] Cosines: */\
		/* add0 + 0x[ 20, 28, 30, 38]: [c4 ,c5 ,c6 ,c7 ] */\
		/* add0 + 0x[ 40, 48, 50, 58]: [c8 ,c9 ,cA ,cB ] */\
		/* add0 + 0x[ 60, 68, 70, 78]: [cC ,cD ,cE ,cF ] */\
		/* add0 + 0x[ 80, 88, 90, 98]: [-- ,r1 ,r2 ,r3 ] Tangents: */\
		/* add0 + 0x[ a0, a8, b0, b8]: [r4 ,r5 ,r6 ,r7 ] */\
		/* add0 + 0x[ c0, c8, d0, d8]: [r8 ,r9 ,rA ,rB ] */\
		/* add0 + 0x[ e0, e8, f0, f8]: [rC ,rD ,rE ,rF ] */\
		/* add0 + 0x[100,108,110,118]: [c31,c51,c62,  0] Cosine ratios: */\
		/* add0 + 0x[120,128,130,138]: [c73,c91,cA2,cB3] */\
		/* add0 + 0x[140,148,150,158]: [cC4,cD5,cE6,cF7] */\
		/* We start at the upper end of the memory chunk allocated for sincos data and proceed downward, so the 11 memory */\
		/* slots inited with scalar data above are only overwritten after their scalar-data contents have been copied: */\
		"vbroadcastsd	0x0f8(%%rax),%%ymm15	/* sF*/	\n\t"/* Allow several cycles for each vbroadcast to complete before writing result */\
		"vbroadcastsd	0x158(%%rax),%%ymm14	/* cF*/	\n\t	movq	%[__cos] ,%%rbx			\n\t"\
		"vbroadcastsd	0x0b8(%%rax),%%ymm13	/* s7*/	\n\t	movq	%[__tan] ,%%rcx			\n\t"\
		"vbroadcastsd	0x120(%%rax),%%ymm12	/* c7*/	\n\t"\
		"vbroadcastsd	0x0d8(%%rax),%%ymm11	/* sB*/	\n\t	vmovaps	%%ymm15,0x420(%%rax)	\n\t"/* s15 = cc0 + 0x21;	__rF = s15/c15		*/\
		"vbroadcastsd	0x138(%%rax),%%ymm10	/* cB*/	\n\t	vmovaps	%%ymm14,0x400(%%rax)	\n\t"/* c15 = cc0 + 0x20;	__cF7 = __cF/__c7	*/\
		"vbroadcastsd	0x098(%%rax),%%ymm9 	/* s3*/	\n\t	vmovaps	%%ymm13,0x3e0(%%rax)	\n\t"/* s7  = cc0 + 0x1f;	__r7 = s7 /c7 		*/\
		"vbroadcastsd	0x100(%%rax),%%ymm8 	/* c3*/	\n\t	vmovaps	%%ymm12,0x3c0(%%rax)	\n\t"/* c7  = cc0 + 0x1e;	__c73 = __c7/__c3	*/\
		"vbroadcastsd	0x0e8(%%rax),%%ymm7 	/* sD*/	\n\t	vmovaps	%%ymm11,0x3a0(%%rax)	\n\t"/* s11 = cc0 + 0x1d;	__rB = s11/c11		*/\
		"vbroadcastsd	0x148(%%rax),%%ymm6 	/* cD*/	\n\t	vmovaps	%%ymm10,0x380(%%rax)	\n\t"/* c11 = cc0 + 0x1c;	__cB3 = __cB/__c3	*/\
		"vbroadcastsd	0x0a8(%%rax),%%ymm5 	/* s5*/	\n\t	vmovaps	%%ymm9 ,0x360(%%rax)	\n\t"/* s3  = cc0 + 0x1b;	__r3 = s3 /c3 		*/\
		"vbroadcastsd	0x108(%%rax),%%ymm4 	/* c5*/	\n\t	vmovaps	%%ymm8 ,0x340(%%rax)	\n\t"/* c3  = cc0 + 0x1a;	__c31 = __c3/__c1	*/\
		"vbroadcastsd	0x0c8(%%rax),%%ymm3 	/* s9*/	\n\t	vmovaps	%%ymm7 ,0x320(%%rax)	\n\t"/* s13 = cc0 + 0x19;	__rD = s13/c13		*/\
		"vbroadcastsd	0x128(%%rax),%%ymm2 	/* c9*/	\n\t	vmovaps	%%ymm6 ,0x300(%%rax)	\n\t"/* c13 = cc0 + 0x18;	__cD5 = __cD/__c5	*/\
		"vbroadcastsd	0x088(%%rax),%%ymm1 	/* s1*/	\n\t	vmovaps	%%ymm5 ,0x2e0(%%rax)	\n\t"/* s5  = cc0 + 0x17;	__r5 = s5 /c5 		*/\
		"vbroadcastsd	0x0f0(%%rax),%%ymm0 	/* sE*/	\n\t	vmovaps	%%ymm4 ,0x2c0(%%rax)	\n\t"/* c5  = cc0 + 0x16;	__c51 = __c5/__c1	*/\
		"vbroadcastsd	0x150(%%rax),%%ymm15	/* cE*/	\n\t	vmovaps	%%ymm3 ,0x2a0(%%rax)	\n\t"/* s9  = cc0 + 0x15;	__r9 = s9 /c9 		*/\
		"vbroadcastsd	0x0b0(%%rax),%%ymm14	/* s6*/	\n\t	vmovaps	%%ymm2 ,0x280(%%rax)	\n\t"/* c9  = cc0 + 0x14;	__c91 = __c9/__c1	*/\
		"vbroadcastsd	0x110(%%rax),%%ymm13	/* c6*/	\n\t	vmovaps	%%ymm1 ,0x260(%%rax)	\n\t"/* s1  = cc0 + 0x13;	__r1 = s1 /c1 		*/\
		"vbroadcastsd	0x0d0(%%rax),%%ymm12	/* sA*/	\n\t	vmovaps	%%ymm0 ,0x220(%%rax)	\n\t"/* s14 = cc0 + 0x11;	__rE = s14/c14		*/\
		"vbroadcastsd	0x130(%%rax),%%ymm11	/* cA*/	\n\t	vmovaps	%%ymm15,0x200(%%rax)	\n\t"/* c14 = cc0 + 0x10;	__cE6 = __cE/__c6	*/\
		"vbroadcastsd	0x090(%%rax),%%ymm10	/* s2*/	\n\t	vmovaps	%%ymm14,0x1e0(%%rax)	\n\t"/* s6  = cc0 + 0x0f;	__r6 = s6 /c6 		*/\
		"vbroadcastsd	0x0e0(%%rax),%%ymm9 	/* sC*/	\n\t	vmovaps	%%ymm13,0x1c0(%%rax)	\n\t"/* c6  = cc0 + 0x0e;	__c62 = __c6/__c2	*/\
		"vbroadcastsd	0x140(%%rax),%%ymm8 	/* cC*/	\n\t	vmovaps	%%ymm12,0x1a0(%%rax)	\n\t"/* s10 = cc0 + 0x0d;	__rA = s10/c10		*/\
		"vbroadcastsd	0x0a0(%%rax),%%ymm7 	/* s4*/	\n\t	vmovaps	%%ymm11,0x180(%%rax)	\n\t"/* c10 = cc0 + 0x0c;	__cA2 = __cA/__c2	*/\
		"vbroadcastsd	0x020(%%rax),%%ymm6 	/* c4*/	\n\t	vmovaps	%%ymm10,0x160(%%rax)	\n\t"/* s2  = cc0 + 0x0b;	__r2 = s2 /c2 		*/\
		"vbroadcastsd	0x0c0(%%rax),%%ymm5 	/* s8*/	\n\t	vmovaps	%%ymm9 ,0x120(%%rax)	\n\t"/* s12 = cc0 + 0x09;	__rC = s12/c12		*/\
		"vbroadcastsd	0x040(%%rax),%%ymm4 	/* c8*/	\n\t	vmovaps	%%ymm8 ,0x100(%%rax)	\n\t"/* c12 = cc0 + 0x08;	__cC4 = __cC/__c4	*/\
		/* Lowest-addressed few data need special handling, and are written in ascending-address order: */\
		"vbroadcastsd	     (%%rbx),%%ymm1	/*__c */	\n\t	vmovaps	%%ymm7 ,0x0e0(%%rax)	\n\t"/* s4  = cc0 + 0x07;	__r4 = s4 /c4 		*/\
		"vbroadcastsd	0x008(%%rax),%%ymm0	/* c1 */	\n\t	vmovaps	%%ymm6 ,0x0c0(%%rax)	\n\t"/* c4  = cc0 + 0x06;	__c4 [unchanged]	*/\
		"vbroadcastsd	0x010(%%rax),%%ymm3	/* c2 */	\n\t	vmovaps	%%ymm5 ,0x0a0(%%rax)	\n\t"/* s8  = cc0 + 0x05;	__r8 = s8 /c8 		*/\
		"vmulpd		   %%ymm1,%%ymm0,%%ymm1	/* c1*__c */\n\t	vmovaps	%%ymm4 ,0x080(%%rax)	\n\t"/* c8  = cc0 + 0x04;	__c8 [unchanged]	*/\
		"vbroadcastsd	     (%%rcx),%%ymm2	/*__s/__c */\n\t	vmovaps	%%ymm0 ,0x240(%%rax)	\n\t"/* c1  = cc0 + 0x12;	__c1 [unchanged]	*/\
		"vmulpd	-0x020(%%rax),%%ymm0,%%ymm0	/*c1*ISRT2*/\n\t	vmovaps	%%ymm3 ,0x140(%%rax)	\n\t"/* c2  = cc0 + 0x0a;	__c2 [unchanged]	*/\
		"vmulpd	-0x020(%%rax),%%ymm3,%%ymm3	/*c2*ISRT2*/\n\t	vmovaps	%%ymm1 ,     (%%rax)	\n\t"/* cc0 = cc0 + 0x00;	__c1_c = c1*__c		*/\
		"														vmovaps	%%ymm2 ,0x020(%%rax)	\n\t"/* ss0 = cc0 + 0x01;	__sc = __s/__c		*/\
		"														vmovaps	%%ymm0 ,0x040(%%rax)	\n\t"/* c0  = cc0 + 0x02;	__c1i2 = c1*ISRT2	*/\
		"														vmovaps	%%ymm3 ,0x060(%%rax)	\n\t"/* s0  = cc0 + 0x03;	__c2i2 = c2*ISRT2	*/\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__one] "m" (Xone)\
		 ,[__cos] "m" (Xcos)\
		 ,[__tan] "m" (Xtan)\
		: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	// Here are the hexadecimal byte offsets w.r.to the __c = cc0 base-root address of the various derived sincos terms:
	// Datum	Offset
	// ------	------
	// __c1_c	0x000
	// __sc		0x020
	// __c1i2	0x040
	// __c2i2	0x060
	// __c8		0x080	__r8		0x0A0
	// __c4		0x0C0	__r4		0x0E0
	// __cC4	0x100	__rC		0x120
	// __c2		0x140	__r2		0x160
	// __cA2	0x180	__rA		0x1A0
	// __c62	0x1C0	__r6		0x1E0
	// __cE6	0x200	__rE		0x220
	// __c1		0x240	__r1		0x260
	// __c91	0x280	__r9		0x2A0
	// __c51	0x2C0	__r5		0x2E0
	// __cD5	0x300	__rD		0x320
	// __c31	0x340	__r3		0x360
	// __cB3	0x380	__rB		0x3A0
	// __c73	0x3C0	__r7		0x3E0
	// __cF7	0x400	__rF		0x420

	// Remember that for AVX2-style 3-operand FMA in AT&T syntax, the result overwrites the rightmost input!
	// 4 ADD, 170 FMA, 0 pure-MUL, 269 MOVAPS
	#define SSE2_RADIX16_DIF_TWIDDLE_1(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr1,Xisrt2,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/*...Block 1:	*/\
		"movq	%[__add0],%%rax				\n\t"\
	"movslq	%[__pfetch_addr1],%%r13	\n\t"	/* Prefetch base-index offset cycles among add0 + p0,1,2,3 on successive macro calls */\
	"leaq	(%%rax,%%r13,8),%%r13	\n\t"	/* [base-address + data-fetch-ahead index] */\
		"movslq	%[__p4],%%rbx				\n\t"\
		"movslq	%[__p8],%%rcx				\n\t"\
	"shlq	$3,%%rbx	\n\t"\
		"movslq	%[__p12],%%rdx				\n\t"\
		"movq	%[__isrt2],%%rsi 			\n\t"\
	"movq	%%rbx,%%r14	\n\t"	/* Save a copy of p4 pointer offset. Will prefetch from [base-address + data-fetch-ahead index] + [0,p4,p8,p12] on each macro call. */\
	"prefetcht1	%c[__pfetch_dist](%%r13)\n\t"/* [base-address + data-fetch-ahead index] + p0 */\
		"addq	 %%rax,%%rbx		\n\t"\
		"leaq	(%%rax,%%rcx,8),%%rcx		\n\t"\
		"leaq	(%%rax,%%rdx,8),%%rdx		\n\t"\
		"movq	%[__r1],%%rdi	\n\t"/* ptr to local-mem */\
		"addq	$0x20,%%rsi 	\n\t"/* cc0 */\
		"vmovaps		     (%%rcx),%%ymm4 		\n\t	vmovaps			0x020(%%rcx),%%ymm5 		\n\t"/*	t04 =__A8r;					t05 =__A8i; */\
		"vmovaps		     (%%rax),%%ymm0 		\n\t	vmovaps			0x020(%%rax),%%ymm1 		\n\t"/*	t00 =__A0r;					t01 =__A0i; */\
		"vmovaps		%%ymm4,%%ymm6				\n\t	vmovaps			0x0a0(%%rsi),%%ymm13		\n\t"/*	t06 = t04; load __r8 into ymm13 */\
		"vmovaps		0x0e0(%%rsi),%%ymm14		\n\t	vmovaps			0x120(%%rsi),%%ymm15		\n\t"/* load __r4,rC into ymm14,ymm15 */\
		"vfnmadd231pd	%%ymm5 ,%%ymm13,%%ymm4 		\n\t	 vfmadd231pd	%%ymm6 ,%%ymm13,%%ymm5 		\n\t"/*	FNMA231(  t05,__r8,t04);	 FMA231(  t06,__r8,t05); */\
		"vmovaps		     (%%rbx),%%ymm8			\n\t	vmovaps			0x020(%%rbx),%%ymm9 		\n\t"/*	_a =__A4r;					_b =__A4i; */\
		"vfnmadd231pd	0x020(%%rbx),%%ymm14,%%ymm8 \n\t	 vfmadd231pd	     (%%rbx),%%ymm14,%%ymm9 \n\t"/*	FNMA231(__A4i,__r4,_a );	 FMA231(__A4r,__r4,_b ); */\
		"vmovaps		     (%%rdx),%%ymm6			\n\t	vmovaps			0x020(%%rdx),%%ymm7 		\n\t"/*	t06 =__ACr;					t07 =__ACi; */\
		"vfnmadd231pd	0x020(%%rdx),%%ymm15,%%ymm6 \n\t	 vfmadd231pd	     (%%rdx),%%ymm15,%%ymm7 \n\t"/*	FNMA231(__ACi,__rC,t06);	 FMA231(__ACr,__rC,t07); */\
		"vmovaps		0x100(%%rsi),%%ymm13		\n\t	vmovaps			0x080(%%rsi),%%ymm14		\n\t"/* load __cC4,c8 into pair of regs */\
		"vmovaps		%%ymm8 ,%%ymm10				\n\t	vmovaps			%%ymm0,%%ymm2 				\n\t"/*	_c = _a;	t02 = t00; */\
		" vfmadd231pd	%%ymm6 ,%%ymm13,%%ymm8 		\n\t	 vfmadd231pd	%%ymm4 ,%%ymm14,%%ymm0 		\n\t"/*	 FMA231(t06,__cC4,_a);		 FMA231(t04,__c8,t00); */\
		"vmovaps		%%ymm9 ,%%ymm11				\n\t	vmovaps			%%ymm1 ,%%ymm3 				\n\t"/*	_d = _b;	t03 = t01; */\
		" vfmadd231pd	%%ymm7 ,%%ymm13,%%ymm9 		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm14,%%ymm1 		\n\t"/*	 FMA231(t07,__cC4,_b);		 FMA231(t05,__c8,t01); */\
		"vfnmadd231pd	%%ymm6 ,%%ymm13,%%ymm10		\n\t	vfnmadd231pd	%%ymm4 ,%%ymm14,%%ymm2 		\n\t"/*	FNMA231(t06,__cC4,_c);		FNMA231(t04,__c8,t02); */\
		"vfnmadd231pd	%%ymm7 ,%%ymm13,%%ymm11		\n\t	vfnmadd231pd	%%ymm5 ,%%ymm14,%%ymm3 		\n\t"/*	FNMA231(t07,__cC4,_d);		FNMA231(t05,__c8,t03); */\
		"vmovaps		0x0c0(%%rsi),%%ymm15		\n\t												\n\t"/* load __c4 into reg */\
		"vmovaps		%%ymm0 ,%%ymm4 				\n\t	vmovaps			%%ymm1 ,%%ymm5 				\n\t"/*	t04 =t00; t05 =t01; */\
		"vfnmadd231pd	%%ymm8 ,%%ymm15,%%ymm4 		\n\t	 vfmadd231pd	%%ymm8 ,%%ymm15,%%ymm0 		\n\t"/*	FNMA231(_a ,__c4 ,t04);		 FMA231(_a ,__c4 ,t00); */\
		"vfnmadd231pd	%%ymm9 ,%%ymm15,%%ymm5 		\n\t	 vfmadd231pd	%%ymm9 ,%%ymm15,%%ymm1 		\n\t"/*	FNMA231(_b ,__c4 ,t05);		 FMA231(_b ,__c4 ,t01); */\
		"vmovaps		%%ymm2 ,%%ymm6 				\n\t	vmovaps			%%ymm3 ,%%ymm7 				\n\t"/*	t06 =t02;	t07 =t03; */\
		" vfmadd231pd	%%ymm11,%%ymm15,%%ymm6 		\n\t	vfnmadd231pd	%%ymm11,%%ymm15,%%ymm2 		\n\t"/*	 FMA231(_d ,__c4 ,t06);		FNMA231(_d ,__c4 ,t02); */\
		"vfnmadd231pd	%%ymm10,%%ymm15,%%ymm7 		\n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm3 		\n\t"/*	FNMA231(_c ,__c4 ,t07);		 FMA231(_c ,__c4 ,t03); */\
		/* Write outputs into local store, and preload 4 inputs for the next block: */\
		"movslq	%[__p2],%%r9		\n\t"\
		"shlq	$3,%%r9				\n\t"\
		"vmovaps	%%ymm4,0x080(%%rdi)	\n\t	vmovaps	%%ymm0 ,     (%%rdi)	\n\t	vmovaps (%%r9,%%rax),%%ymm0  \n\t"\
		"vmovaps	%%ymm5,0x0a0(%%rdi)	\n\t	vmovaps	%%ymm1 ,0x020(%%rdi)	\n\t	vmovaps (%%r9,%%rcx),%%ymm4  \n\t"\
		"vmovaps	%%ymm6,0x0c0(%%rdi)	\n\t	vmovaps	%%ymm2 ,0x040(%%rdi)	\n\t	vmovaps (%%r9,%%rbx),%%ymm8  \n\t"\
		"vmovaps	%%ymm7,0x0e0(%%rdi)	\n\t	vmovaps	%%ymm3 ,0x060(%%rdi)	\n\t	vmovaps (%%r9,%%rdx),%%ymm6  \n\t"\
		"\n\t"\
	/*...Block 2: Register indices for the 8 t-date = [t-index - 8]: */\
		"addq	$0x100,%%rsi 		\n\t"/* cc += 8 */\
		"addq	$0x100,%%rdi 		\n\t"/* r1 += 8 */\
		"vmovaps	0x020(%%r9,%%rax),%%ymm1	\n\t	vmovaps 0x060(%%rsi),%%ymm12	\n\t"/* t08 =__A2r; t09 =__A2i; load __r2 into ymm12 */\
		"vmovaps	0x020(%%r9,%%rcx),%%ymm5	\n\t	vmovaps 0x0a0(%%rsi),%%ymm13	\n\t"/* t12 =__AAr; t13 =__AAi; load __rA into ymm13 */\
		"vmovaps	0x020(%%r9,%%rbx),%%ymm9	\n\t	vmovaps 0x0e0(%%rsi),%%ymm14	\n\t"/* _a  =__A6r; _b  =__A6i; load __r6 into ymm14 */\
		"vmovaps	0x020(%%r9,%%rdx),%%ymm7	\n\t	vmovaps 0x120(%%rsi),%%ymm15	\n\t"/* t14 =__AEr; t15 =__AEi; load __rE into ymm15 */\
		/* 2,3,10,11 unused - use to hold copies of mem-operand data in the ensuing 8 FMAs: */\
		"vmovaps	%%ymm0,%%ymm2 	\n\t	vmovaps		%%ymm1,%%ymm3 	\n\t"/*	m0,1 copies */\
		"vmovaps	%%ymm4,%%ymm10	\n\t	vmovaps		%%ymm5,%%ymm11	\n\t"/*	m4,5 copies */\
		"vfnmadd231pd	%%ymm3 ,%%ymm12,%%ymm0 \n\t	 vfmadd231pd	%%ymm2 ,%%ymm12,%%ymm1  	\n\t"/* FNMA231(__A2i,__r2,t08); FMA231(__A2r,__r2,t09); */\
		"vfnmadd231pd	%%ymm11,%%ymm13,%%ymm4 \n\t	 vfmadd231pd	%%ymm10,%%ymm13,%%ymm5  	\n\t"/* FNMA231(__AAi,__rA,t12); FMA231(__AAr,__rA,t13); */\
		"vmovaps	%%ymm8,%%ymm2 	\n\t	vmovaps		%%ymm9,%%ymm3 	\n\t"/*	m8,9 copies */\
		"vmovaps	%%ymm6,%%ymm10	\n\t	vmovaps		%%ymm7,%%ymm11	\n\t"/*	m6,7 copies */\
		"vfnmadd231pd	%%ymm3 ,%%ymm14,%%ymm8 \n\t	 vfmadd231pd	%%ymm2 ,%%ymm14,%%ymm9  	\n\t"/* FNMA231(__A6i,__r6,_a ); FMA231(__A6r,__r6,_b ); */\
		"vfnmadd231pd	%%ymm11,%%ymm15,%%ymm6 \n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm7  	\n\t"/* FNMA231(__AEi,__rE,t14); FMA231(__AEr,__rE,t15); */\
		"vmovaps		0x100(%%rsi),%%ymm13		\n\t	vmovaps			0x080(%%rsi),%%ymm14		\n\t"/* load __cE6,cA2 into pair of regs */\
		"vmovaps		%%ymm8 ,%%ymm10				\n\t	vmovaps			%%ymm0,%%ymm2 				\n\t"/*	_c = _a;	t10 = t08; */\
		" vfmadd231pd	%%ymm6 ,%%ymm13,%%ymm8 		\n\t	 vfmadd231pd	%%ymm4 ,%%ymm14,%%ymm0 		\n\t"/*	 FMA231(t14,__cE6,_a);		 FMA231(t12,__cA2,t08); */\
		"vmovaps		%%ymm9 ,%%ymm11				\n\t	vmovaps			%%ymm1 ,%%ymm3 				\n\t"/*	_d = _b;	t11 = t09; */\
		" vfmadd231pd	%%ymm7 ,%%ymm13,%%ymm9 		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm14,%%ymm1 		\n\t"/*	 FMA231(t15,__cE6,_b);		 FMA231(t13,__cA2,t09); */\
		"vfnmadd231pd	%%ymm6 ,%%ymm13,%%ymm10		\n\t	vfnmadd231pd	%%ymm4 ,%%ymm14,%%ymm2 		\n\t"/*	FNMA231(t14,__cE6,_c);		FNMA231(t12,__cA2,t10); */\
		"vfnmadd231pd	%%ymm7 ,%%ymm13,%%ymm11		\n\t	vfnmadd231pd	%%ymm5 ,%%ymm14,%%ymm3 		\n\t"/*	FNMA231(t15,__cE6,_d);		FNMA231(t13,__cA2,t11); */\
		"vmovaps		0x0c0(%%rsi),%%ymm15		\n\t												\n\t"/* load __c62 into reg */\
		"vmovaps		%%ymm0 ,%%ymm4 				\n\t	vmovaps			%%ymm1 ,%%ymm5 				\n\t"/*	t12 =t08 ;	t13 =t09; */\
		"vfnmadd231pd	%%ymm8 ,%%ymm15,%%ymm4 		\n\t	 vfmadd231pd	%%ymm8 ,%%ymm15,%%ymm0 		\n\t"/*	FNMA231(_a,__c62,t12);		 FMA231( _a,__c62,t08); */\
		"vfnmadd231pd	%%ymm9 ,%%ymm15,%%ymm5 		\n\t	 vfmadd231pd	%%ymm9 ,%%ymm15,%%ymm1 		\n\t"/*	FNMA231(_b,__c62,t13);		 FMA231( _b,__c62,t09); */\
		"vmovaps		%%ymm2 ,%%ymm6 				\n\t	vmovaps			%%ymm3 ,%%ymm7 				\n\t"/*	t14 =t10;	t15 =t11; */\
		" vfmadd231pd	%%ymm11,%%ymm15,%%ymm6 		\n\t	vfnmadd231pd	%%ymm11,%%ymm15,%%ymm2 		\n\t"/*	 FMA231(_d,__c62,t14);		FNMA231( _d,__c62,t10); */\
		"vfnmadd231pd	%%ymm10,%%ymm15,%%ymm7 		\n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm3 		\n\t"/*	FNMA231(_c,__c62,t15);		 FMA231( _c,__c62,t11); */\
		/* Write outputs into local store, and preload 4 inputs for the next block: */\
		"movslq	%[__p1],%%r8		\n\t"\
		"shlq	$3,%%r8				\n\t"\
		"vmovaps %%ymm4,0x080(%%rdi) \n\t vmovaps %%ymm0 ,     (%%rdi) \n\t vmovaps (%%r8,%%rax),%%ymm0  \n\t"\
		"vmovaps %%ymm5,0x0a0(%%rdi) \n\t vmovaps %%ymm1 ,0x020(%%rdi) \n\t vmovaps (%%r8,%%rcx),%%ymm4  \n\t"\
		"vmovaps %%ymm6,0x0c0(%%rdi) \n\t vmovaps %%ymm2 ,0x040(%%rdi) \n\t vmovaps (%%r8,%%rbx),%%ymm8  \n\t"\
		"vmovaps %%ymm7,0x0e0(%%rdi) \n\t vmovaps %%ymm3 ,0x060(%%rdi) \n\t vmovaps (%%r8,%%rdx),%%ymm6  \n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%r14)\n\t"/* ...+p4 */\
		"\n\t"\
	/*...Block 3: Register indices for the 8 t-date = [t-index - 16]: */\
		"addq	$0x100,%%rsi 		\n\t"/* cc += 8 */\
		"addq	$0x100,%%rdi 		\n\t"/* r1 += 8 */\
		"vmovaps	0x020(%%r8,%%rax),%%ymm1	\n\t	vmovaps 0x060(%%rsi),%%ymm12	\n\t"/* t16 =__A1r;	t17 =__A1i; load __r1 into ymm12  */\
		"vmovaps	0x020(%%r8,%%rcx),%%ymm5	\n\t	vmovaps 0x0a0(%%rsi),%%ymm13	\n\t"/* t20 =__A9r;	t21 =__A9i; load __r9 into ymm13  */\
		"vmovaps	0x020(%%r8,%%rbx),%%ymm9	\n\t	vmovaps 0x0e0(%%rsi),%%ymm14	\n\t"/* _a=  __A5r;	_b  =__A5i; load __r5 into ymm14  */\
		"vmovaps	0x020(%%r8,%%rdx),%%ymm7	\n\t	vmovaps 0x120(%%rsi),%%ymm15	\n\t"/* t22 =__ADr;	t23 =__ADi; load __rD into ymm15  */\
		"vmovaps	%%ymm0,%%ymm2 	\n\t	vmovaps		%%ymm1,%%ymm3 	\n\t"/*	m0,1 copies */\
		"vmovaps	%%ymm4,%%ymm10	\n\t	vmovaps		%%ymm5,%%ymm11	\n\t"/*	m4,5 copies */\
		"vfnmadd231pd	%%ymm3 ,%%ymm12,%%ymm0 \n\t	 vfmadd231pd	%%ymm2 ,%%ymm12,%%ymm1  	\n\t"/* FNMA231(__A1i,__r1,t16);	 FMA231(__A1r,__r1,t17); */\
		"vfnmadd231pd	%%ymm11,%%ymm13,%%ymm4 \n\t	 vfmadd231pd	%%ymm10,%%ymm13,%%ymm5  	\n\t"/* FNMA231(__A9i,__r9,t20);	 FMA231(__A9r,__r9,t21); */\
		"vmovaps	%%ymm8,%%ymm2 	\n\t	vmovaps		%%ymm9,%%ymm3 	\n\t"/*	m8,9 copies */\
		"vmovaps	%%ymm6,%%ymm10	\n\t	vmovaps		%%ymm7,%%ymm11	\n\t"/*	m6,7 copies */\
		"vfnmadd231pd	%%ymm3 ,%%ymm14,%%ymm8 \n\t	 vfmadd231pd	%%ymm2 ,%%ymm14,%%ymm9  	\n\t"/* FNMA231(__A5i,__r5,_a );	 FMA231(__A5r,__r5,_b ); */\
		"vfnmadd231pd	%%ymm11,%%ymm15,%%ymm6 \n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm7  	\n\t"/* FNMA231(__ADi,__rD,t22);	 FMA231(__ADr,__rD,t23); */\
		"vmovaps		0x100(%%rsi),%%ymm13		\n\t	vmovaps			0x080(%%rsi),%%ymm14		\n\t"/* load __cD5,c91 into pair of regs */\
		"vmovaps		%%ymm8 ,%%ymm10				\n\t	vmovaps			%%ymm0,%%ymm2 				\n\t"/*	_c= _a;	t18= t16; */\
		" vfmadd231pd	%%ymm6 ,%%ymm13,%%ymm8 		\n\t	 vfmadd231pd	%%ymm4 ,%%ymm14,%%ymm0 		\n\t"/*	 FMA231(t22,__cD5,_a);		 FMA231(t20,__c91,t16); */\
		"vmovaps		%%ymm9 ,%%ymm11				\n\t	vmovaps			%%ymm1 ,%%ymm3 				\n\t"/*	_d= _b;	t19= t17; */\
		" vfmadd231pd	%%ymm7 ,%%ymm13,%%ymm9 		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm14,%%ymm1 		\n\t"/*	 FMA231(t23,__cD5,_b);		 FMA231(t21,__c91,t17); */\
		"vfnmadd231pd	%%ymm6 ,%%ymm13,%%ymm10		\n\t	vfnmadd231pd	%%ymm4 ,%%ymm14,%%ymm2 		\n\t"/*	FNMA231(t22,__cD5,_c);		FNMA231(t20,__c91,t18); */\
		"vfnmadd231pd	%%ymm7 ,%%ymm13,%%ymm11		\n\t	vfnmadd231pd	%%ymm5 ,%%ymm14,%%ymm3 		\n\t"/*	FNMA231(t23,__cD5,_d);		FNMA231(t21,__c91,t19); */\
		"vmovaps		0x0c0(%%rsi),%%ymm15		\n\t												\n\t"/* load __c51 into reg */\
		"vmovaps		%%ymm0 ,%%ymm4 				\n\t	vmovaps			%%ymm1 ,%%ymm5 				\n\t"/*	t20 =t16;	t21 =t17; */\
		"vfnmadd231pd	%%ymm8 ,%%ymm15,%%ymm4 		\n\t	 vfmadd231pd	%%ymm8 ,%%ymm15,%%ymm0 		\n\t"/*	FNMA231(_a,__c51,t20);		 FMA231(_a,__c51,t16); */\
		"vfnmadd231pd	%%ymm9 ,%%ymm15,%%ymm5 		\n\t	 vfmadd231pd	%%ymm9 ,%%ymm15,%%ymm1 		\n\t"/*	FNMA231(_b,__c51,t21);		 FMA231(_b,__c51,t17); */\
		"vmovaps		%%ymm2 ,%%ymm6 				\n\t	vmovaps			%%ymm3 ,%%ymm7 				\n\t"/*	t22 =t18;	t23 =t19; */\
		" vfmadd231pd	%%ymm11,%%ymm15,%%ymm6 		\n\t	vfnmadd231pd	%%ymm11,%%ymm15,%%ymm2 		\n\t"/*	 FMA231(_d,__c51,t22);		FNMA231(_d,__c51,t18); */\
		"vfnmadd231pd	%%ymm10,%%ymm15,%%ymm7 		\n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm3 		\n\t"/*	FNMA231(_c,__c51,t23);		 FMA231(_c,__c51,t19); */\
		/* Write outputs into local store, and preload 4 inputs for the next block: */\
		"addq	%%r8,%%r9				\n\t"/* p3<<3, overwrites p2<<3 */\
		"vmovaps %%ymm4,0x080(%%rdi) \n\t vmovaps %%ymm0 ,     (%%rdi) \n\t vmovaps (%%r9,%%rax),%%ymm0  \n\t"\
		"vmovaps %%ymm5,0x0a0(%%rdi) \n\t vmovaps %%ymm1 ,0x020(%%rdi) \n\t vmovaps (%%r9,%%rcx),%%ymm4  \n\t"\
		"vmovaps %%ymm6,0x0c0(%%rdi) \n\t vmovaps %%ymm2 ,0x040(%%rdi) \n\t vmovaps (%%r9,%%rbx),%%ymm8  \n\t"\
		"vmovaps %%ymm7,0x0e0(%%rdi) \n\t vmovaps %%ymm3 ,0x060(%%rdi) \n\t vmovaps (%%r9,%%rdx),%%ymm6  \n\t"\
		"\n\t"\
	/*...Block 4: Register indices for the 8 t-date = [t-index - 24]: */\
		"addq	$0x100,%%rsi 		\n\t"/* cc += 8 */\
		"addq	$0x100,%%rdi 		\n\t"/* r1 += 8 */\
		"vmovaps	0x020(%%r9,%%rax),%%ymm1	\n\t	vmovaps 0x060(%%rsi),%%ymm12	\n\t"/* t24 =__A3r;	t25 =__A3i; load __r3 into ymm12 */\
		"vmovaps	0x020(%%r9,%%rcx),%%ymm5	\n\t	vmovaps 0x0a0(%%rsi),%%ymm13	\n\t"/* t28 =__ABr;	t29 =__ABi; load __rB into ymm13 */\
		"vmovaps	0x020(%%r9,%%rbx),%%ymm9	\n\t	vmovaps 0x0e0(%%rsi),%%ymm14	\n\t"/* _a  =__A7r;	_b  =__A7i; load __r7 into ymm14 */\
		"vmovaps	0x020(%%r9,%%rdx),%%ymm7	\n\t	vmovaps 0x120(%%rsi),%%ymm15	\n\t"/* t30 =__AFr;	t31 =__AFi; load __rF into ymm15 */\
		"vmovaps	%%ymm0,%%ymm2 	\n\t	vmovaps		%%ymm1,%%ymm3 	\n\t"/*	m0,1 copies */\
		"vmovaps	%%ymm4,%%ymm10	\n\t	vmovaps		%%ymm5,%%ymm11	\n\t"/*	m4,5 copies */\
		"vfnmadd231pd	%%ymm3 ,%%ymm12,%%ymm0 \n\t	 vfmadd231pd	%%ymm2 ,%%ymm12,%%ymm1  	\n\t"/* FNMA231(__A3i,__r3,t24);	 FMA231(__A3r,__r3,t25); */\
		"vfnmadd231pd	%%ymm11,%%ymm13,%%ymm4 \n\t	 vfmadd231pd	%%ymm10,%%ymm13,%%ymm5  	\n\t"/* FNMA231(__ABi,__rB,t28 );	 FMA231(__ABr,__rB,t29 ); */\
		"vmovaps	%%ymm8,%%ymm2 	\n\t	vmovaps		%%ymm9,%%ymm3 	\n\t"/*	m8,9 copies */\
		"vmovaps	%%ymm6,%%ymm10	\n\t	vmovaps		%%ymm7,%%ymm11	\n\t"/*	m6,7 copies */\
		"vfnmadd231pd	%%ymm3 ,%%ymm14,%%ymm8 \n\t	 vfmadd231pd	%%ymm2 ,%%ymm14,%%ymm9  	\n\t"/* FNMA231(__A7i,__r7,_a);		 FMA231(__A7r,__r7,_b); */\
		"vfnmadd231pd	%%ymm11,%%ymm15,%%ymm6 \n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm7  	\n\t"/* FNMA231(__AFi,__rF,t30 );	 FMA231(__AFr,__rF,t31 ); */\
		"vmovaps		0x100(%%rsi),%%ymm13		\n\t	vmovaps			0x080(%%rsi),%%ymm14		\n\t"/* load __cF7,cB3 into pair of regs */\
		"vmovaps		%%ymm8 ,%%ymm10				\n\t	vmovaps			%%ymm0,%%ymm2 				\n\t"/*	_c= _a;	t26= t24; */\
		" vfmadd231pd	%%ymm6 ,%%ymm13,%%ymm8 		\n\t	 vfmadd231pd	%%ymm4 ,%%ymm14,%%ymm0 		\n\t"/*	 FMA231(t30,__cF7,_a);		 FMA231(t28,__cB3,t24); */\
		"vmovaps		%%ymm9 ,%%ymm11				\n\t	vmovaps			%%ymm1 ,%%ymm3 				\n\t"/*	_d= _b;	t27= t25; */\
		" vfmadd231pd	%%ymm7 ,%%ymm13,%%ymm9 		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm14,%%ymm1 		\n\t"/*	 FMA231(t31,__cF7,_b);		 FMA231(t29,__cB3,t25); */\
		"vfnmadd231pd	%%ymm6 ,%%ymm13,%%ymm10		\n\t	vfnmadd231pd	%%ymm4 ,%%ymm14,%%ymm2 		\n\t"/*	FNMA231(t30,__cF7,_c);		FNMA231(t28,__cB3,t26); */\
		"vfnmadd231pd	%%ymm7 ,%%ymm13,%%ymm11		\n\t	vfnmadd231pd	%%ymm5 ,%%ymm14,%%ymm3 		\n\t"/*	FNMA231(t31,__cF7,_d);		FNMA231(t29,__cB3,t27); */\
		"vmovaps		0x0c0(%%rsi),%%ymm15		\n\t												\n\t"/* load __c73 into reg */\
		"vmovaps		%%ymm0 ,%%ymm4 				\n\t	vmovaps			%%ymm1 ,%%ymm5 				\n\t"/*	t28 =t24;	t29 =t25; */\
		"vfnmadd231pd	%%ymm8 ,%%ymm15,%%ymm4 		\n\t	 vfmadd231pd	%%ymm8 ,%%ymm15,%%ymm0 		\n\t"/*	FNMA231(_a,__c73,t28);		 FMA231(_a,__c73,t24); */\
		"vfnmadd231pd	%%ymm9 ,%%ymm15,%%ymm5 		\n\t	 vfmadd231pd	%%ymm9 ,%%ymm15,%%ymm1 		\n\t"/*	FNMA231(_b,__c73,t29);		 FMA231(_b,__c73,t25); */\
		"vmovaps		%%ymm2 ,%%ymm6 				\n\t	vmovaps			%%ymm3 ,%%ymm7 				\n\t"/*	t30 =t26;	t31 =t27; */\
		" vfmadd231pd	%%ymm11,%%ymm15,%%ymm6 		\n\t	vfnmadd231pd	%%ymm11,%%ymm15,%%ymm2 		\n\t"/*	 FMA231(_d,__c73,t30);		FNMA231(_d,__c73,t26); */\
		"vfnmadd231pd	%%ymm10,%%ymm15,%%ymm7 		\n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm3 		\n\t"/*	FNMA231(_c,__c73,t31);		 FMA231(_c,__c73,t27); */\
		"vmovaps		%%ymm4 ,0x080(%%rdi)		\n\t	vmovaps			%%ymm0 ,     (%%rdi)		\n\t"/* Write outputs into local store */\
		"vmovaps		%%ymm5 ,0x0a0(%%rdi)		\n\t	vmovaps			%%ymm1 ,0x020(%%rdi)		\n\t"\
		"vmovaps		%%ymm6 ,0x0c0(%%rdi)		\n\t	vmovaps			%%ymm2 ,0x040(%%rdi)		\n\t"\
		"vmovaps		%%ymm7 ,0x0e0(%%rdi)		\n\t	vmovaps			%%ymm3 ,0x060(%%rdi)		\n\t"\
		"\n\t"\
	/*************************************************************************************/\
	/*  And now do four more radix-4 transforms, including the internal twiddle factors: */\
	/*************************************************************************************/\
	/* Block 1: */\
		"movq	%[__add0],%%rax		\n\t"\
		"movq	%%r9,%%rdx			\n\t"/* copy p3<<3 to rdx... */\
		"subq	%%r8,%%r9			\n\t"/* and restore p2<<3 to r9. */\
		"leaq	(%%rax,%%r8),%%rbx	\n\t"/* add0+p1 */\
		"leaq	(%%rax,%%r9),%%rcx	\n\t"/* add0+p2 */\
		"addq	%%rax,%%rdx			\n\t"/* add0+p3 */\
		"subq	$0x300,%%rsi 		\n\t"/* revert cc-ptr to base value */\
		"subq	$0x300,%%rdi 		\n\t"/* revert r1-ptr to base value */\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%r14,2)\n\t"/* ...+p8 */\
	"leaq	(%%r14,%%r14,2),%%r14	\n\t"	/* p4 + (p4*2) = p12, ptr-offset form */\
		"\n\t"\
		/*...Read t0,8,16,24 from local store ... Do the 4 Im-part FMAs first, because their results needed 1st below */\
		"vmovaps		0x140(%%rsi),%%ymm14		\n\t	vmovaps			0x340(%%rsi),%%ymm15		\n\t"/* load __c2,c31 into pair of regs */\
		"vmovaps		     (%%rdi),%%ymm0 		\n\t	vmovaps			0x100(%%rdi),%%ymm2 		\n\t"/*    t00;    t08; */\
		"vmovaps		0x200(%%rdi),%%ymm4 		\n\t	vmovaps			0x300(%%rdi),%%ymm6 		\n\t"/*    t16;    t24; */\
		"vmovaps			 %%ymm0 ,%%ymm8 		\n\t	vmovaps				 %%ymm4 ,%%ymm10		\n\t"/* _a=t00; _c=t16; */\
		" vfmadd231pd	%%ymm2 ,%%ymm14,%%ymm0 		\n\t	 vfmadd231pd	%%ymm6 ,%%ymm15,%%ymm4 		\n\t"/*	 FMA231(t08,__c2 ,t00);		 FMA231(t24,__c31,t16); */\
		"vmovaps		0x020(%%rdi),%%ymm1 		\n\t	vmovaps			0x120(%%rdi),%%ymm3 		\n\t"/*    t01;    t09; */\
		"vmovaps		0x220(%%rdi),%%ymm5 		\n\t	vmovaps			0x320(%%rdi),%%ymm7 		\n\t"/*    t17;    t25; */\
		"vmovaps			 %%ymm1 ,%%ymm9 		\n\t	vmovaps				 %%ymm5 ,%%ymm11		\n\t"/* _b=t01; _d=t17; */\
		" vfmadd231pd	%%ymm3 ,%%ymm14,%%ymm1 		\n\t	 vfmadd231pd	%%ymm7 ,%%ymm15,%%ymm5 		\n\t"/*	 FMA231(t09,__c2 ,t01);		 FMA231(t25,__c31,t17); */\
		"vfnmadd231pd	%%ymm2 ,%%ymm14,%%ymm8 		\n\t	vfnmadd231pd	%%ymm6 ,%%ymm15,%%ymm10		\n\t"/*	FNMA231(t08,__c2 ,_a );		FNMA231(t24,__c31,_c ); */\
		"vfnmadd231pd	%%ymm3 ,%%ymm14,%%ymm9 		\n\t	vfnmadd231pd	%%ymm7 ,%%ymm15,%%ymm11		\n\t"/*	FNMA231(t09,__c2 ,_b );		FNMA231(t25,__c31,_d ); */\
		"vmovaps		0x240(%%rsi),%%ymm15		\n\t												\n\t"/* load __c1 into reg */\
		"vmovaps			 %%ymm0 ,%%ymm12		\n\t	vmovaps				 %%ymm1 ,%%ymm13		\n\t"/* _e = t00; _f = t01; */\
		" vfmadd231pd	%%ymm4 ,%%ymm15,%%ymm0 		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm15,%%ymm1 		\n\t"/*	 FMA231(t16,__c1 ,t00);		 FMA231(t17,__c1 ,t01); */\
		"vfnmadd231pd	%%ymm4 ,%%ymm15,%%ymm12		\n\t	vfnmadd231pd	%%ymm5 ,%%ymm15,%%ymm13		\n\t"/*	FNMA231(t16,__c1 ,_e );		FNMA231(t17,__c1 ,_f ); */\
		"vmovaps			 %%ymm8 ,%%ymm2 		\n\t	vmovaps				 %%ymm9 ,%%ymm3 		\n\t"/* t08 = _a ; t09 = _b; */\
		"vfnmadd231pd	%%ymm11,%%ymm15,%%ymm2 		\n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm3 		\n\t"/*	FNMA231(_d ,__c1 ,t08);		 FMA231(_c ,__c1 ,t09); */\
		" vfmadd231pd	%%ymm11,%%ymm15,%%ymm8 		\n\t	vfnmadd231pd	%%ymm10,%%ymm15,%%ymm9 		\n\t"/*	 FMA231(_d ,__c1 ,_a );		FNMA231(_c ,__c1 ,_b ); */\
		/* Write outputs back to main array: */\
		"vmovaps		%%ymm0 ,     (%%rax)		\n\t	vmovaps			%%ymm1 ,0x020(%%rax)		\n\t"/* __B0r= t00;		__B0i= t01; */\
		"vmovaps		%%ymm12,     (%%rbx)		\n\t	vmovaps			%%ymm13,0x020(%%rbx)		\n\t"/* __B1r= _e ;		__B1i= _f ; */\
		"vmovaps		%%ymm2 ,     (%%rcx)		\n\t	vmovaps			%%ymm3 ,0x020(%%rcx)		\n\t"/* __B2r= t08;		__B2i= t09; */\
		"vmovaps		%%ymm8 ,     (%%rdx)		\n\t	vmovaps			%%ymm9 ,0x020(%%rdx)		\n\t"/* __B3r= _a ;		__B3i= _b ; */\
		"\n\t"\
	/*...Block 3: t4,12,20,28 */\
		"movslq	%[__p4],%%r8	\n\t"\
		"shlq	$3,%%r8			\n\t"\
		"addq	%%r8,%%rax		\n\t"/* add0+p4 */\
		"addq	%%r8,%%rbx		\n\t"/* add0+p5 */\
		"addq	%%r8,%%rcx		\n\t"/* add0+p6 */\
		"addq	%%r8,%%rdx		\n\t"/* add0+p7 */\
		"addq	$0x80,%%rdi		\n\t"/* r5 */\
		"vmovaps	0x440(%%rsi),%%ymm15	\n\t"/* cc0 + 0x22 = __two; Actually holds 1.0 in AVX2 mode */\
		"vmovaps		0x200(%%rdi),%%ymm4 		\n\t	vmovaps			0x220(%%rdi),%%ymm5 		\n\t"/*    t20;    t21; */\
		"vmovaps		0x300(%%rdi),%%ymm6 		\n\t	vmovaps			0x320(%%rdi),%%ymm7 		\n\t"/*    t28;    t29; */\
		"vmovaps			 %%ymm4 ,%%ymm10 		\n\t	vmovaps				 %%ymm7 ,%%ymm11		\n\t"/* _c=t20; _d=t29; */\
		"vfnmadd231pd	%%ymm5 ,%%ymm15,%%ymm10		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm15,%%ymm4 		\n\t"/*	FNMA231(t21,1.0,_c );		 FMA231(t21,1.0,t20); */\
		" vfmadd231pd	%%ymm6 ,%%ymm15,%%ymm11		\n\t	vfnmadd231pd	%%ymm6 ,%%ymm15,%%ymm7 		\n\t"/*	 FMA231(t28,1.0,_d );		FNMA231(t28,1.0,t29); */\
		"vmovaps		0x140(%%rsi),%%ymm14		\n\t	vmovaps			0x340(%%rsi),%%ymm15		\n\t"/* load __c2,31 into pair of regs */\
		"vmovaps		     (%%rdi),%%ymm0 		\n\t	vmovaps			0x020(%%rdi),%%ymm1 		\n\t"/*    t04;    t05; */\
		"vmovaps		0x100(%%rdi),%%ymm2 		\n\t	vmovaps			0x120(%%rdi),%%ymm3 		\n\t"/*    t12;    t13; */\
		"vmovaps			 %%ymm0 ,%%ymm8 		\n\t	vmovaps				 %%ymm1 ,%%ymm9 		\n\t"/* _a = t04; _b = t05; */\
		"vmovaps			 %%ymm10,%%ymm5 		\n\t	vmovaps				 %%ymm4 ,%%ymm12		\n\t"/* t21 = _c; _e = t20; */\
		" vfmadd231pd	%%ymm3 ,%%ymm14,%%ymm8 		\n\t	 vfmadd231pd	%%ymm11,%%ymm15,%%ymm5 		\n\t"/*	 FMA231(t13,__c2 ,_a );		 FMA231(_d ,__c31,t21); */\
		"vfnmadd231pd	%%ymm2 ,%%ymm14,%%ymm9 		\n\t	 vfmadd231pd	%%ymm7 ,%%ymm15,%%ymm4 		\n\t"/*	FNMA231(t12,__c2 ,_b );		 FMA231(t29,__c31,t20); */\
		"vfnmadd231pd	%%ymm3 ,%%ymm14,%%ymm0 		\n\t	vfnmadd231pd	%%ymm11,%%ymm15,%%ymm10		\n\t"/*	FNMA231(t13,__c2 ,t04);		FNMA231(_d ,__c31,_c ); */\
		" vfmadd231pd	%%ymm2 ,%%ymm14,%%ymm1 		\n\t	vfnmadd231pd	%%ymm7 ,%%ymm15,%%ymm12		\n\t"/*	 FMA231(t12,__c2 ,t05);		FNMA231(t29,__c31,_e ); */\
		"vmovaps		0x040(%%rsi),%%ymm15		\n\t												\n\t"/* load __c1i2 into reg */\
		"vmovaps			 %%ymm8 ,%%ymm2 		\n\t	vmovaps				 %%ymm9 ,%%ymm3 		\n\t"/* t12 = _a; t13 = _b; */\
		"vfnmadd231pd	%%ymm4 ,%%ymm15,%%ymm2 		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm15,%%ymm3 		\n\t"/*	FNMA231(t20,__c1i2,t12);	 FMA231(t21,__c1i2,t13); */\
		" vfmadd231pd	%%ymm4 ,%%ymm15,%%ymm8 		\n\t	vfnmadd231pd	%%ymm5 ,%%ymm15,%%ymm9 		\n\t"/*	 FMA231(t20,__c1i2,_a );	FNMA231(t21,__c1i2,_b ); */\
		"vmovaps			 %%ymm0 ,%%ymm11		\n\t	vmovaps				 %%ymm1 ,%%ymm7 		\n\t"/* _d = t04; t29 = t05; */\
		" vfmadd231pd	%%ymm10,%%ymm15,%%ymm0 		\n\t	 vfmadd231pd	%%ymm12,%%ymm15,%%ymm1 		\n\t"/*	 FMA231(_c ,__c1i2,t04);	 FMA231(_e ,__c1i2,t05); */\
		"vfnmadd231pd	%%ymm10,%%ymm15,%%ymm11		\n\t	vfnmadd231pd	%%ymm12,%%ymm15,%%ymm7 		\n\t"/*	FNMA231(_c ,__c1i2,_d );	FNMA231(_e ,__c1i2,t29); */\
		/* Write outputs back to main array: */\
		"vmovaps		%%ymm2 ,     (%%rcx)		\n\t	vmovaps			%%ymm3 ,0x020(%%rcx)		\n\t"/* __B6r= t12;		__B6i= t13; */\
		"vmovaps		%%ymm8 ,     (%%rdx)		\n\t	vmovaps			%%ymm9 ,0x020(%%rdx)		\n\t"/* __B7r= _a ;		__B7i= _b ; */\
		"vmovaps		%%ymm0 ,     (%%rax)		\n\t	vmovaps			%%ymm1 ,0x020(%%rax)		\n\t"/* __B4r= t04;		__B4i= t05; */\
		"vmovaps		%%ymm11,     (%%rbx)		\n\t	vmovaps			%%ymm7 ,0x020(%%rbx)		\n\t"/* __B5r= _d ;		__B5i= t29; */\
		"\n\t"\
	/*...Block 2: t2,10,18,26 */\
		"subq	$0x40,%%rdi		\n\t"/* r3 */\
		"vmovaps		0x100(%%rdi),%%ymm2 		\n\t	vmovaps			0x120(%%rdi),%%ymm3 		\n\t"/*    t10;    t11; */\
		"addq	%%r8,%%rax		\n\t"/* add0+p8 */\
		"addq	%%r8,%%rbx		\n\t"/* add0+p9 */\
		"addq	%%r8,%%rcx		\n\t"/* add0+pA */\
		"addq	%%r8,%%rdx		\n\t"/* add0+pB */\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%r14)\n\t"/* ...+p12 */\
		"vmovaps		0x200(%%rdi),%%ymm4 		\n\t	vmovaps			0x220(%%rdi),%%ymm5 		\n\t"/*    t18;    t19; */\
		"vmovaps		0x300(%%rdi),%%ymm6 		\n\t	vmovaps			0x320(%%rdi),%%ymm7 		\n\t"/*    t26;    t27; */\
		"vsubpd		%%ymm3 ,%%ymm2 ,%%ymm12 		\n\t"/* _e = t10-t11; */\
		"vmovaps		0x020(%%rsi),%%ymm15		\n\t"/* load __sc into reg */\
		"vaddpd		%%ymm2 ,%%ymm3 ,%%ymm13			\n\t"/* _f = t10+t11; */\
		"vmovaps			 %%ymm4 ,%%ymm10 		\n\t	vmovaps				 %%ymm7 ,%%ymm8 		\n\t"/* _c = t18; _a = t27; */\
		"vfnmadd231pd	%%ymm5 ,%%ymm15,%%ymm10		\n\t	 vfmsub231pd	%%ymm6 ,%%ymm15,%%ymm8 		\n\t"/*	FNMA231(t19,__sc,_c );		 FMS231(t26,__sc,_a ); */\
		"vmovaps			 %%ymm5 ,%%ymm11 		\n\t	vmovaps				 %%ymm6 ,%%ymm9 		\n\t"/* _d = t19; _b = t26; */\
		" vfmadd231pd	%%ymm4 ,%%ymm15,%%ymm11		\n\t	 vfmadd231pd	%%ymm7 ,%%ymm15,%%ymm9 		\n\t"/*	 FMA231(t18,__sc,_d );		 FMA231(t27,__sc,_b ); */\
		"vmovaps		     (%%rdi),%%ymm0 		\n\t	vmovaps			0x020(%%rdi),%%ymm1 		\n\t"/*    t02;    t03; */\
		"vmovaps		0x340(%%rsi),%%ymm14		\n\t	vmovaps			0x060(%%rsi),%%ymm15		\n\t"/* load __c31,c2i2 into pair of regs */\
		"vmovaps			 %%ymm10,%%ymm4 		\n\t	vmovaps				 %%ymm0 ,%%ymm2 		\n\t"/* t18 = _c;	t10 = t02; */\
		" vfmadd231pd	%%ymm8 ,%%ymm14,%%ymm4 		\n\t	 vfmadd231pd	%%ymm12,%%ymm15,%%ymm0 		\n\t"/*	 FMA231(_a ,__c31,t18);		 FMA231(_e ,__c2i2,t02); */\
		"vmovaps			 %%ymm11,%%ymm5 		\n\t	vmovaps				 %%ymm1 ,%%ymm3 		\n\t"/* t19 = _d;	t11 = t03; */\
		" vfmadd231pd	%%ymm9 ,%%ymm14,%%ymm5 		\n\t	 vfmadd231pd	%%ymm13,%%ymm15,%%ymm1 		\n\t"/*	 FMA231(_b ,__c31,t19);		 FMA231(_f ,__c2i2,t03); */\
		"vfnmadd231pd	%%ymm8 ,%%ymm14,%%ymm10		\n\t	vfnmadd231pd	%%ymm12,%%ymm15,%%ymm2 		\n\t"/*	FNMA231(_a ,__c31,_c );		FNMA231(_e ,__c2i2,t10); */\
		"vfnmadd231pd	%%ymm9 ,%%ymm14,%%ymm11		\n\t	vfnmadd231pd	%%ymm13,%%ymm15,%%ymm3 		\n\t"/*	FNMA231(_b ,__c31,_d );		FNMA231(_f ,__c2i2,t11); */\
		"vmovaps		(%%rsi),%%ymm15				\n\t												\n\t"/* load __c1_c into reg */\
		"vmovaps			 %%ymm0 ,%%ymm8 		\n\t	vmovaps				 %%ymm1 ,%%ymm9 		\n\t"/* _a = t02; _b = t03; */\
		" vfmadd231pd	%%ymm4 ,%%ymm15,%%ymm0 		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm15,%%ymm1 		\n\t"/*	 FMA231(t18,__c1_c,t02);	 FMA231(t19,__c1_c,t03); */\
		"vfnmadd231pd	%%ymm4 ,%%ymm15,%%ymm8 		\n\t	vfnmadd231pd	%%ymm5 ,%%ymm15,%%ymm9 		\n\t"/*	FNMA231(t18,__c1_c,_a );	FNMA231(t19,__c1_c,_b ); */\
		"vmovaps			 %%ymm2 ,%%ymm12		\n\t	vmovaps				 %%ymm3 ,%%ymm13		\n\t"/* _e = t10; _f = t11; */\
		"vfnmadd231pd	%%ymm11,%%ymm15,%%ymm2 		\n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm3 		\n\t"/*	FNMA231(_d ,__c1_c,t10);	 FMA231(_c ,__c1_c,t11); */\
		" vfmadd231pd	%%ymm11,%%ymm15,%%ymm12		\n\t	vfnmadd231pd	%%ymm10,%%ymm15,%%ymm13		\n\t"/*	 FMA231(_d ,__c1_c,_e );	FNMA231(_c ,__c1_c,_f ); */\
		/* Write outputs back to main array: */\
		"vmovaps		%%ymm0 ,     (%%rax)		\n\t	vmovaps			%%ymm1 ,0x020(%%rax)		\n\t"/* __B8r= t02;		__B8i= t03; */\
		"vmovaps		%%ymm8 ,     (%%rbx)		\n\t	vmovaps			%%ymm9 ,0x020(%%rbx)		\n\t"/* __B9r= _a ;		__B9i= _b ; */\
		"vmovaps		%%ymm2 ,     (%%rcx)		\n\t	vmovaps			%%ymm3 ,0x020(%%rcx)		\n\t"/* __BAr= t10;		__BAi= t11; */\
		"vmovaps		%%ymm12,     (%%rdx)		\n\t	vmovaps			%%ymm13,0x020(%%rdx)		\n\t"/* __BBr= _e ;		__BBi= _f ; */\
		"\n\t"\
	/*...Block 4: t6,14,22,30 */\
		"addq	$0x80,%%rdi		\n\t"/* r7 */\
		"vmovaps		0x100(%%rdi),%%ymm2 		\n\t	vmovaps			0x120(%%rdi),%%ymm3 		\n\t"/*    t14;    t15; */\
		"addq	%%r8,%%rax		\n\t"/* add0+pC */\
		"addq	%%r8,%%rbx		\n\t"/* add0+pD */\
		"addq	%%r8,%%rcx		\n\t"/* add0+pE */\
		"addq	%%r8,%%rdx		\n\t"/* add0+pF */\
		"vmovaps		0x200(%%rdi),%%ymm4 		\n\t	vmovaps			0x220(%%rdi),%%ymm5 		\n\t"/*    t22;    t23; */\
		"vmovaps		0x300(%%rdi),%%ymm6 		\n\t	vmovaps			0x320(%%rdi),%%ymm7 		\n\t"/*    t30;    t31; */\
		"vaddpd		%%ymm2 ,%%ymm3 ,%%ymm10 		\n\t"/* _c = t14+t5; */\
		"vmovaps		0x020(%%rsi),%%ymm15		\n\t"/* load __sc into reg */\
		"vsubpd		%%ymm2 ,%%ymm3 ,%%ymm11			\n\t"/* _d = t15-t14; */\
		"vmovaps			 %%ymm5 ,%%ymm12 		\n\t	vmovaps				 %%ymm4 ,%%ymm13		\n\t"/* _e = t23; _f = t22;*/\
		" vfmsub231pd	%%ymm4 ,%%ymm15,%%ymm12		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm15,%%ymm13		\n\t"/*	 FMS231(t22,__sc,_e );		 FMA231(t23,__sc,_f );*/\
		"vmovaps			 %%ymm6 ,%%ymm8  		\n\t	vmovaps				 %%ymm7 ,%%ymm9 		\n\t"/* _a = t30; _b = t31; */\
		"vfnmadd231pd	%%ymm7 ,%%ymm15,%%ymm8 		\n\t	 vfmadd231pd	%%ymm6 ,%%ymm15,%%ymm9 		\n\t"/*	FNMA231(t31,__sc,_a );		 FMA231(t30,__sc,_b );*/\
		"vmovaps		     (%%rdi),%%ymm0 		\n\t	vmovaps			0x020(%%rdi),%%ymm1 		\n\t"/*    t06;    t07; */\
		"vmovaps		0x340(%%rsi),%%ymm14		\n\t	vmovaps			0x060(%%rsi),%%ymm15		\n\t"/* load __c31,c2i2 into pair of regs */\
		"vmovaps			 %%ymm1 ,%%ymm3 		\n\t	vmovaps				 %%ymm0 ,%%ymm2 		\n\t"/* t15= t07;	t14= t06; */\
		"vfnmadd231pd	%%ymm11,%%ymm15,%%ymm1 		\n\t	vfnmadd231pd	%%ymm10,%%ymm15,%%ymm0 		\n\t"/*	FNMA231(_d ,__c2i2,t07);	FNMA231(_c ,__c2i2,t06); */\
		"vmovaps			 %%ymm12,%%ymm4 		\n\t	vmovaps				 %%ymm13,%%ymm5 		\n\t"/* t22= _e; t23= _f; */\
		"vfnmadd231pd	%%ymm8 ,%%ymm14,%%ymm4 		\n\t	vfnmadd231pd	%%ymm9 ,%%ymm14,%%ymm5 		\n\t"/*	FNMA231(_a ,__c31 ,t22);	FNMA231(_b ,__c31 ,t23); */\
		" vfmadd231pd	%%ymm8 ,%%ymm14,%%ymm12		\n\t	 vfmadd231pd	%%ymm9 ,%%ymm14,%%ymm13		\n\t"/*	 FMA231(_a ,__c31 ,_e );	 FMA231(_b ,__c31 ,_f ); */\
		" vfmadd231pd	%%ymm10,%%ymm15,%%ymm2 		\n\t	 vfmadd231pd	%%ymm11,%%ymm15,%%ymm3 		\n\t"/*	 FMA231(_c ,__c2i2,t14);	 FMA231(_d ,__c2i2,t15); */\
		"vmovaps		(%%rsi),%%ymm15				\n\t												\n\t"/* load __c1_c into reg */\
		"vmovaps			 %%ymm0 ,%%ymm8 		\n\t	vmovaps				 %%ymm1 ,%%ymm9 		\n\t"/* _a = t06; _b = t07; */\
		" vfmadd231pd	%%ymm4 ,%%ymm15,%%ymm0 		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm15,%%ymm1 		\n\t"/*	 FMA231(t22,__c1_c,t06);	 FMA231(t23,__c1_c,t07); */\
		"vfnmadd231pd	%%ymm4 ,%%ymm15,%%ymm8 		\n\t	vfnmadd231pd	%%ymm5 ,%%ymm15,%%ymm9 		\n\t"/*	FNMA231(t22,__c1_c,_a );	FNMA231(t23,__c1_c,_b ); */\
		"vmovaps			 %%ymm2 ,%%ymm10		\n\t	vmovaps				 %%ymm3 ,%%ymm11		\n\t"/* _c = t14; _d = t15; */\
		"vfnmadd231pd	%%ymm13,%%ymm15,%%ymm2 		\n\t	 vfmadd231pd	%%ymm12,%%ymm15,%%ymm3 		\n\t"/*	FNMA231(_f ,__c1_c,t14);	 FMA231(_e ,__c1_c,t15); */\
		" vfmadd231pd	%%ymm13,%%ymm15,%%ymm10		\n\t	vfnmadd231pd	%%ymm12,%%ymm15,%%ymm11		\n\t"/*	 FMA231(_f ,__c1_c,_c );	FNMA231(_e ,__c1_c,_d ); */\
		/* Write outputs back to main array: */\
		"vmovaps		%%ymm0 ,     (%%rax)		\n\t	vmovaps			%%ymm1 ,0x020(%%rax)		\n\t"/* __BCr= t06;		__BCi= t07; */\
		"vmovaps		%%ymm8 ,     (%%rbx)		\n\t	vmovaps			%%ymm9 ,0x020(%%rbx)		\n\t"/* __BDr= _a ;		__BDi= _b ; */\
		"vmovaps		%%ymm2 ,     (%%rcx)		\n\t	vmovaps			%%ymm3 ,0x020(%%rcx)		\n\t"/* __BEr= t14;		__BEi= t15; */\
		"vmovaps		%%ymm10,     (%%rdx)		\n\t	vmovaps			%%ymm11,0x020(%%rdx)		\n\t"/* __BFr= _c ;		__BFi= _d ; */\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r1] "m" (Xr1)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	// 2nd version of trig-data and DFT-pass routines uses 1-copy trig data, read in and broadcast to full YMM register on-the-fly:
	#define RADIX16_COMPUTE_FMA_SINCOS_DIF_2(Xadd0,Xone)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax			\n\t"/* base-ptr to 32-byte-aligned 128-byte-large memchunk 16 doubles needing iterative inversion */\
		"movq	%[__one] ,%%rbx			\n\t"/* 1.0 in 4-fold-double form */\
		"vmovaps	    (%%rbx),%%ymm14	\n\t"/* 1.0 */\
		"vmovaps	    (%%rax),%%ymm4	\n\t"/* load inputs-to-be-inverted (d's) into 4 AVX registers */\
		"vmovaps	0x20(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x40(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x60(%%rax),%%ymm7	\n\t"\
		"vaddpd	%%ymm14,%%ymm14,%%ymm14	\n\t"/* 2.0 */\
		"vcvtpd2ps	%%ymm4,%%xmm0	\n\t"/* convert d's to SP ... Note in AVX mode output register *must* be a 128-bit xmm! */\
		"vcvtpd2ps	%%ymm5,%%xmm1	\n\t"\
		"vcvtpd2ps	%%ymm6,%%xmm2	\n\t"\
		"vcvtpd2ps	%%ymm7,%%xmm3	\n\t"\
		"vrcpps		%%xmm0,%%xmm0	\n\t"	/* ainv := approx 1/d to 11-12 bits of precision */\
		"vrcpps		%%xmm1,%%xmm1	\n\t"\
		"vrcpps		%%xmm2,%%xmm2	\n\t"\
		"vrcpps		%%xmm3,%%xmm3	\n\t"\
		"vcvtps2pd	%%xmm0,%%ymm0	\n\t"	/* convert ~1/d back to DP  ... Note in AVX mode input register *must* be a 128-bit xmm!*/\
		"vcvtps2pd	%%xmm1,%%ymm1	\n\t"\
		"vcvtps2pd	%%xmm2,%%ymm2	\n\t"\
		"vcvtps2pd	%%xmm3,%%ymm3	\n\t"\
		/* 1st NR iteration gives ~23 bits of precision: */\
		"vmovaps	%%ymm0,%%ymm8 	\n\t"	/* make a copy of ainv */\
		"vmovaps	%%ymm1,%%ymm9 	\n\t"\
		"vmovaps	%%ymm2,%%ymm10	\n\t"\
		"vmovaps	%%ymm3,%%ymm11	\n\t"\
		"vfnmadd132pd	%%ymm4,%%ymm14,%%ymm0 	\n\t"	/* 2 - d*ainv, overwrites ainv */\
		"vfnmadd132pd	%%ymm5,%%ymm14,%%ymm1 	\n\t"\
		"vfnmadd132pd	%%ymm6,%%ymm14,%%ymm2 	\n\t"\
		"vfnmadd132pd	%%ymm7,%%ymm14,%%ymm3 	\n\t"\
		"vmulpd		%%ymm0,%%ymm8 ,%%ymm0	\n\t"	/* ainv*(2 - d*ainv) = 1/d accurate to ~23 bits */\
		"vmulpd		%%ymm1,%%ymm9 ,%%ymm1	\n\t"\
		"vmulpd		%%ymm2,%%ymm10,%%ymm2	\n\t"\
		"vmulpd		%%ymm3,%%ymm11,%%ymm3	\n\t"\
		/* 3rd-order update of 23-bit result needs just 2 FMA, 1 SUB, 1 MUL: */\
		"vaddpd	    (%%rbx),%%ymm14,%%ymm14	\n\t"/* 3.0, needed by 3rd-order update step */\
		"vmovaps	%%ymm0,%%ymm8 	\n\t"	/* make a copy of ainv */\
		"vmovaps	%%ymm1,%%ymm9 	\n\t"\
		"vmovaps	%%ymm2,%%ymm10	\n\t"\
		"vmovaps	%%ymm3,%%ymm11	\n\t"\
		"vfnmadd132pd	%%ymm0,%%ymm14,%%ymm4 	\n\t"/* 1st FMA overwrites d data (inputs) with (3 - d*ainv) */\
		"vfnmadd132pd	%%ymm1,%%ymm14,%%ymm5 	\n\t"\
		"vfnmadd132pd	%%ymm2,%%ymm14,%%ymm6 	\n\t"\
		"vfnmadd132pd	%%ymm3,%%ymm14,%%ymm7 	\n\t"\
		"vsubpd		%%ymm14,%%ymm4,%%ymm0 	\n\t"/* Subtract 3 from (3 - d*ainv) to get -y = -d*ainv terms in ymm0-3 */\
		"vsubpd		%%ymm14,%%ymm5,%%ymm1 	\n\t"\
		"vsubpd		%%ymm14,%%ymm6,%%ymm2 	\n\t"\
		"vsubpd		%%ymm14,%%ymm7,%%ymm3 	\n\t"\
		"vfmadd132pd	%%ymm4,%%ymm14,%%ymm0 	\n\t"/* Positive-product FMA gives (3 - y*(3 - d*ainv)) in ymm0-3*/\
		"vfmadd132pd	%%ymm5,%%ymm14,%%ymm1 	\n\t"\
		"vfmadd132pd	%%ymm6,%%ymm14,%%ymm2 	\n\t"\
		"vfmadd132pd	%%ymm7,%%ymm14,%%ymm3 	\n\t"\
		"vmulpd		%%ymm0,%%ymm8 ,%%ymm0	\n\t"	/* ainv*(3 - y*(3 - d*ainv)) = 1/d accurate to ~53 bits */\
		"vmulpd		%%ymm1,%%ymm9 ,%%ymm1	\n\t"\
		"vmulpd		%%ymm2,%%ymm10,%%ymm2	\n\t"\
		"vmulpd		%%ymm3,%%ymm11,%%ymm3	\n\t"\
		/* Multiply by the sine terms to get quotients in ymm4-7 (still need inverses in ymm0-3 for cosine-ratios): */\
		"vmulpd	0x80(%%rax),%%ymm0,%%ymm4	\n\t"\
		"vmulpd	0xa0(%%rax),%%ymm1,%%ymm5	\n\t"\
		"vmulpd	0xc0(%%rax),%%ymm2,%%ymm6	\n\t"\
		"vmulpd	0xe0(%%rax),%%ymm3,%%ymm7	\n\t"\
		/* Outputs overwrite sines: */\
		"vmovaps	%%ymm4,0x80(%%rax)	\n\t"\
		"vmovaps	%%ymm5,0xa0(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0xc0(%%rax)	\n\t"\
		"vmovaps	%%ymm7,0xe0(%%rax)	\n\t"\
		/* Set up for cosine ratios: */\
		"vpermilpd	$11,%%ymm0,%%ymm2 	\n\t"/* permute ymm0 = 1/c[3123] to get ymm2 = 1/c[1123], then *= [c3,c5,c6,0] */\
		"vmulpd	0x100(%%rax),%%ymm2,%%ymm2	\n\t"/* ymm2 = [c3,c5,c6, 0] * 1/c[1123] = [c31,c51,c62,  0] */\
		"vmulpd	0x120(%%rax),%%ymm0,%%ymm0	\n\t"/* ymm0 = [c7,c9,cA,cB] * 1/c[3123] = [c73,c91,cA2,cB3] */\
		"vmulpd	0x140(%%rax),%%ymm1,%%ymm1	\n\t"/* ymm1 = [cC,cD,cE,cF] * 1/c[4567] = [cC4,cD5,cE6,cF7] */\
		/* Outputs into slots imm. above above tangents: */\
		"vmovaps	%%ymm2,0x100(%%rax)	\n\t"\
		"vmovaps	%%ymm0,0x120(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x140(%%rax)	\n\t"\
		/* Scalar data starting at add0 = cc0 laid out as below. Ensuing C code assumed to massage into a scalar-data analog of the 4-copy layout. */\
		/* add0 + 0x[  0,  8, 10, 18]: [c3 ,c1 ,c2 ,c3 ] Cosines: */\
		/* add0 + 0x[ 20, 28, 30, 38]: [c4 ,c5 ,c6 ,c7 ] */\
		/* add0 + 0x[ 40, 48, 50, 58]: [c8 ,c9 ,cA ,cB ] */\
		/* add0 + 0x[ 60, 68, 70, 78]: [cC ,cD ,cE ,cF ] */\
		/* add0 + 0x[ 80, 88, 90, 98]: [-- ,r1 ,r2 ,r3 ] Tangents: */\
		/* add0 + 0x[ a0, a8, b0, b8]: [r4 ,r5 ,r6 ,r7 ] */\
		/* add0 + 0x[ c0, c8, d0, d8]: [r8 ,r9 ,rA ,rB ] */\
		/* add0 + 0x[ e0, e8, f0, f8]: [rC ,rD ,rE ,rF ] */\
		/* add0 + 0x[100,108,110,118]: [c31,c51,c62,  0] Cosine ratios: */\
		/* add0 + 0x[120,128,130,138]: [c73,c91,cA2,cB3] */\
		/* add0 + 0x[140,148,150,158]: [cC4,cD5,cE6,cF7] */\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__one] "m" (Xone)\
		: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm14"	/* Clobbered registers */\
	);\
	}

	/* Array-of-doubles-index and byte offsets w.r.to the __c = cc0 base-root address
	of the various derived sincos terms:

		a[0x00,0x01,0x02,0x03]: add0 + 0x[  0,  8, 10, 18]: [1.0,c1 ,c2 ,c1*c]
		a[0x04,0x05,0x06,0x07]: add0 + 0x[ 20, 28, 30, 38]: [c4 ,c1I,c2I,---]
		a[0x08,0x09,0x0a,0x0b]: add0 + 0x[ 40, 48, 50, 58]: [c8 ,---,---,---]
		a[0x0c,0x0d,0x0e,0x0f]: add0 + 0x[ 60, 68, 70, 78]: [---,---,---,---]
		a[0x10,0x11,0x12,0x13]: add0 + 0x[ 80, 88, 90, 98]: [s/c,r1 ,r2 ,r3 ] Tangents:
		a[0x14,0x15,0x16,0x17]: add0 + 0x[ a0, a8, b0, b8]: [r4 ,r5 ,r6 ,r7 ]
		a[0x18,0x19,0x1a,0x1b]: add0 + 0x[ c0, c8, d0, d8]: [r8 ,r9 ,rA ,rB ]
		a[0x1c,0x1d,0x1e,0x1f]: add0 + 0x[ e0, e8, f0, f8]: [rC ,rD ,rE ,rF ]
		a[0x20,0x21,0x22,0x23]: add0 + 0x[100,108,110,118]: [c31,c51,c62,---] Cosine ratios:
		a[0x24,0x25,0x26,0x27]: add0 + 0x[120,128,130,138]: [c73,c91,cA2,cB3]
		a[0x28,0x29,0x2a,0x2b]: add0 + 0x[140,148,150,158]: [cC4,cD5,cE6,cF7]
	*/
	// Remember that for AVX2-style 3-operand FMA in AT&T syntax, the result overwrites the rightmost input!
	// 4 ADD, 170 FMA, 0 pure-MUL], 203 MOVAPS [but lots of mem-operand-using FMAs]
	#define SSE2_RADIX16_DIF_TWIDDLE_2(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr1,Xcc0,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
		/*...Block 1:	*/\
		"movq	%[__add0],%%rax				\n\t"\
	"movq	%%rax,%%r10	\n\t"	/* Prefetch base-index offset cycles among add0 + p0,1,2,3 on successive macro calls */\
		"movslq	%[__p4],%%rbx				\n\t"\
		"movslq	%[__p8],%%rcx				\n\t"\
	"shlq	$3,%%rbx	\n\t"\
		"movslq	%[__p12],%%rdx				\n\t"\
	"mov	%%rbx,%%r14		\n\t"/* Save copies of p2,p4 pointer offsets. */\
	"movslq	%[__p2],%%r13	\n\t"/* Will prefetch from [base-address + data-fetch-ahead index] */\
	"shlq	$3,%%r13		\n\t"/* + p[2,4,6,8,10,12,14] in first half of macro, p1 + p[2,4,6,8,10,12,14] in 2nd half. */\
	"prefetcht1	%c[__pfetch_dist](%%r10)\n\t"/* [base-address + data-fetch-ahead index] + p0 */\
		"movq	%[__cc0],%%rsi 			\n\t"\
		"vbroadcastsd	0x0c0(%%rsi),%%ymm13 \n\t vbroadcastsd	0x0a0(%%rsi),%%ymm14 \n\t vbroadcastsd	0x0e0(%%rsi),%%ymm15 \n\t"/* load __r8,r4,rC into ymm13-15 */\
		"addq	 %%rax,%%rbx		\n\t"\
		"leaq	(%%rax,%%rcx,8),%%rcx		\n\t"\
		"leaq	(%%rax,%%rdx,8),%%rdx		\n\t"\
		"movq	%[__r1],%%rdi	\n\t"/* ptr to local-mem */\
		"vmovaps		     (%%rcx),%%ymm4 		\n\t	vmovaps			0x020(%%rcx),%%ymm5 		\n\t"/*	t04 =__A8r;					t05 =__A8i; */\
		"vmovaps		     (%%rax),%%ymm0 		\n\t	vmovaps			0x020(%%rax),%%ymm1 		\n\t"/*	t00 =__A0r;					t01 =__A0i; */\
		"vmovaps		%%ymm4,%%ymm6				\n\t"/*	t06 = t04; */\
		"vfnmadd231pd	%%ymm5 ,%%ymm13,%%ymm4 		\n\t	 vfmadd231pd	%%ymm6 ,%%ymm13,%%ymm5 		\n\t"/*	FNMA231(  t05,__r8,t04);	 FMA231(  t06,__r8,t05); */\
		"vmovaps		     (%%rbx),%%ymm8			\n\t	vmovaps			0x020(%%rbx),%%ymm9 		\n\t"/*	_a =__A4r;					_b =__A4i; */\
		"vfnmadd231pd	0x020(%%rbx),%%ymm14,%%ymm8 \n\t	 vfmadd231pd	     (%%rbx),%%ymm14,%%ymm9 \n\t"/*	FNMA231(__A4i,__r4,_a );	 FMA231(__A4r,__r4,_b ); */\
		"vbroadcastsd	0x140(%%rsi),%%ymm13		\n\t	vbroadcastsd	0x040(%%rsi),%%ymm14		\n\t"/* load __cC4,c8 into pair of regs */\
		"vmovaps		     (%%rdx),%%ymm6			\n\t	vmovaps			0x020(%%rdx),%%ymm7 		\n\t"/*	t06 =__ACr;					t07 =__ACi; */\
		"vfnmadd231pd	0x020(%%rdx),%%ymm15,%%ymm6 \n\t	 vfmadd231pd	     (%%rdx),%%ymm15,%%ymm7 \n\t"/*	FNMA231(__ACi,__rC,t06);	 FMA231(__ACr,__rC,t07); */\
		"vbroadcastsd	0x020(%%rsi),%%ymm15		\n\t"/* load __c4 */\
		"vmovaps		%%ymm8 ,%%ymm10				\n\t	vmovaps			%%ymm0,%%ymm2 				\n\t"/*	_c = _a;	t02 = t00; */\
		" vfmadd231pd	%%ymm6 ,%%ymm13,%%ymm8 		\n\t	 vfmadd231pd	%%ymm4 ,%%ymm14,%%ymm0 		\n\t"/*	 FMA231(t06,__cC4,_a);		 FMA231(t04,__c8,t00); */\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r13)\n\t"/* ...+p2 */\
		"vmovaps		%%ymm9 ,%%ymm11				\n\t	vmovaps			%%ymm1 ,%%ymm3 				\n\t"/*	_d = _b;	t03 = t01; */\
		" vfmadd231pd	%%ymm7 ,%%ymm13,%%ymm9 		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm14,%%ymm1 		\n\t"/*	 FMA231(t07,__cC4,_b);		 FMA231(t05,__c8,t01); */\
		"vfnmadd231pd	%%ymm6 ,%%ymm13,%%ymm10		\n\t	vfnmadd231pd	%%ymm4 ,%%ymm14,%%ymm2 		\n\t"/*	FNMA231(t06,__cC4,_c);		FNMA231(t04,__c8,t02); */\
		"vfnmadd231pd	%%ymm7 ,%%ymm13,%%ymm11		\n\t	vfnmadd231pd	%%ymm5 ,%%ymm14,%%ymm3 		\n\t"/*	FNMA231(t07,__cC4,_d);		FNMA231(t05,__c8,t03); */\
		"vmovaps		%%ymm0 ,%%ymm4 				\n\t	vmovaps			%%ymm1 ,%%ymm5 				\n\t"/*	t04 =t00; t05 =t01; */\
		"vfnmadd231pd	%%ymm8 ,%%ymm15,%%ymm4 		\n\t	 vfmadd231pd	%%ymm8 ,%%ymm15,%%ymm0 		\n\t"/*	FNMA231(_a ,__c4 ,t04);		 FMA231(_a ,__c4 ,t00); */\
		"vfnmadd231pd	%%ymm9 ,%%ymm15,%%ymm5 		\n\t	 vfmadd231pd	%%ymm9 ,%%ymm15,%%ymm1 		\n\t"/*	FNMA231(_b ,__c4 ,t05);		 FMA231(_b ,__c4 ,t01); */\
		"vmovaps		%%ymm2 ,%%ymm6 				\n\t	vmovaps			%%ymm3 ,%%ymm7 				\n\t"/*	t06 =t02;	t07 =t03; */\
		" vfmadd231pd	%%ymm11,%%ymm15,%%ymm6 		\n\t	vfnmadd231pd	%%ymm11,%%ymm15,%%ymm2 		\n\t"/*	 FMA231(_d ,__c4 ,t06);		FNMA231(_d ,__c4 ,t02); */\
		"vfnmadd231pd	%%ymm10,%%ymm15,%%ymm7 		\n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm3 		\n\t"/*	FNMA231(_c ,__c4 ,t07);		 FMA231(_c ,__c4 ,t03); */\
		"vmovaps		%%ymm4 ,0x080(%%rdi)		\n\t	vmovaps			%%ymm0 ,     (%%rdi)		\n\t"/* Write outputs into local store */\
		"vmovaps		%%ymm5 ,0x0a0(%%rdi)		\n\t	vmovaps			%%ymm1 ,0x020(%%rdi)		\n\t"\
		"vmovaps		%%ymm6 ,0x0c0(%%rdi)		\n\t	vmovaps			%%ymm2 ,0x040(%%rdi)		\n\t"\
		"vmovaps		%%ymm7 ,0x0e0(%%rdi)		\n\t	vmovaps			%%ymm3 ,0x060(%%rdi)		\n\t"\
		"\n\t"\
		/*...Block 2: Register indices for the 8 t-date = [t-index - 8]: */\
		"movslq	%[__p2],%%r9		\n\t"\
		"shlq	$3,%%r9				\n\t"\
		"addq	%%r9,%%rax			\n\t"\
		"addq	%%r9,%%rbx			\n\t"\
		"addq	%%r9,%%rcx			\n\t"\
		"addq	%%r9,%%rdx			\n\t"\
		"addq	$0x100,%%rdi 		/* r1 += 8 */\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r14)\n\t"/* ...+p4 */\
	"addq	%%r14,%%r13 		\n\t"/* p2 -> p6 */\
		"vbroadcastsd	0x090(%%rsi),%%ymm12 	\n\t"/* load __r2 into ymm12 */\
		"vbroadcastsd	0x0d0(%%rsi),%%ymm13 	\n\t"/* load __rA into ymm13 */\
		"vbroadcastsd	0x0b0(%%rsi),%%ymm14 	\n\t"/* load __r6 into ymm14 */\
		"vbroadcastsd	0x0f0(%%rsi),%%ymm15 	\n\t"/* load __rE into ymm15 */\
		"vmovaps (%%rax),%%ymm0  \n\t vmovaps 0x020(%%rax),%%ymm1  \n\t"/* t08 =__A2r; t09 =__A2i; */\
		"vmovaps (%%rcx),%%ymm4  \n\t vmovaps 0x020(%%rcx),%%ymm5  \n\t"/* t12 =__AAr; t13 =__AAi; */\
		"vmovaps (%%rbx),%%ymm8  \n\t vmovaps 0x020(%%rbx),%%ymm9  \n\t"/* _a  =__A6r; _b  =__A6i; */\
		"vmovaps (%%rdx),%%ymm6  \n\t vmovaps 0x020(%%rdx),%%ymm7  \n\t"/* t14 =__AEr; t15 =__AEi; */\
		"vfnmadd231pd	0x020(%%rax),%%ymm12,%%ymm0 \n\t	 vfmadd231pd	(%%rax),%%ymm12,%%ymm1  	\n\t"/* FNMA231(__A2i,__r2,t08); FMA231(__A2r,__r2,t09); */\
		"vfnmadd231pd	0x020(%%rcx),%%ymm13,%%ymm4 \n\t	 vfmadd231pd	(%%rcx),%%ymm13,%%ymm5  	\n\t"/* FNMA231(__AAi,__rA,t12); FMA231(__AAr,__rA,t13); */\
		"vbroadcastsd	0x150(%%rsi),%%ymm13		\n\t"/* load __cE6 */\
		"vfnmadd231pd	0x020(%%rbx),%%ymm14,%%ymm8 \n\t	 vfmadd231pd	(%%rbx),%%ymm14,%%ymm9  	\n\t"/* FNMA231(__A6i,__r6,_a ); FMA231(__A6r,__r6,_b ); */\
		"vbroadcastsd	0x130(%%rsi),%%ymm14		\n\t"/* load __cA2 */\
		"vfnmadd231pd	0x020(%%rdx),%%ymm15,%%ymm6 \n\t	 vfmadd231pd	(%%rdx),%%ymm15,%%ymm7  	\n\t"/* FNMA231(__AEi,__rE,t14); FMA231(__AEr,__rE,t15); */\
		"vbroadcastsd	0x110(%%rsi),%%ymm15		\n\t"/* load __c62 */\
		"vmovaps		%%ymm8 ,%%ymm10				\n\t	vmovaps			%%ymm0,%%ymm2 				\n\t"/*	_c = _a;	t10 = t08; */\
		" vfmadd231pd	%%ymm6 ,%%ymm13,%%ymm8 		\n\t	 vfmadd231pd	%%ymm4 ,%%ymm14,%%ymm0 		\n\t"/*	 FMA231(t14,__cE6,_a);		 FMA231(t12,__cA2,t08); */\
		"vmovaps		%%ymm9 ,%%ymm11				\n\t	vmovaps			%%ymm1 ,%%ymm3 				\n\t"/*	_d = _b;	t11 = t09; */\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r13)\n\t"/* ...+p6 */\
		" vfmadd231pd	%%ymm7 ,%%ymm13,%%ymm9 		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm14,%%ymm1 		\n\t"/*	 FMA231(t15,__cE6,_b);		 FMA231(t13,__cA2,t09); */\
		"vfnmadd231pd	%%ymm6 ,%%ymm13,%%ymm10		\n\t	vfnmadd231pd	%%ymm4 ,%%ymm14,%%ymm2 		\n\t"/*	FNMA231(t14,__cE6,_c);		FNMA231(t12,__cA2,t10); */\
		"vfnmadd231pd	%%ymm7 ,%%ymm13,%%ymm11		\n\t	vfnmadd231pd	%%ymm5 ,%%ymm14,%%ymm3 		\n\t"/*	FNMA231(t15,__cE6,_d);		FNMA231(t13,__cA2,t11); */\
		"vmovaps		%%ymm0 ,%%ymm4 				\n\t	vmovaps			%%ymm1 ,%%ymm5 				\n\t"/*	t12 =t08 ;	t13 =t09; */\
		"vfnmadd231pd	%%ymm8 ,%%ymm15,%%ymm4 		\n\t	 vfmadd231pd	%%ymm8 ,%%ymm15,%%ymm0 		\n\t"/*	FNMA231(_a,__c62,t12);		 FMA231( _a,__c62,t08); */\
		"vfnmadd231pd	%%ymm9 ,%%ymm15,%%ymm5 		\n\t	 vfmadd231pd	%%ymm9 ,%%ymm15,%%ymm1 		\n\t"/*	FNMA231(_b,__c62,t13);		 FMA231( _b,__c62,t09); */\
		"vmovaps		%%ymm2 ,%%ymm6 				\n\t	vmovaps			%%ymm3 ,%%ymm7 				\n\t"/*	t14 =t10;	t15 =t11; */\
		" vfmadd231pd	%%ymm11,%%ymm15,%%ymm6 		\n\t	vfnmadd231pd	%%ymm11,%%ymm15,%%ymm2 		\n\t"/*	 FMA231(_d,__c62,t14);		FNMA231( _d,__c62,t10); */\
		"vfnmadd231pd	%%ymm10,%%ymm15,%%ymm7 		\n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm3 		\n\t"/*	FNMA231(_c,__c62,t15);		 FMA231( _c,__c62,t11); */\
		"vmovaps		%%ymm4 ,0x080(%%rdi)		\n\t	vmovaps			%%ymm0 ,     (%%rdi)		\n\t"/* Write outputs into local store */\
		"vmovaps		%%ymm5 ,0x0a0(%%rdi)		\n\t	vmovaps			%%ymm1 ,0x020(%%rdi)		\n\t"\
		"vmovaps		%%ymm6 ,0x0c0(%%rdi)		\n\t	vmovaps			%%ymm2 ,0x040(%%rdi)		\n\t"\
		"vmovaps		%%ymm7 ,0x0e0(%%rdi)		\n\t	vmovaps			%%ymm3 ,0x060(%%rdi)		\n\t"\
		"/* Due to array-padding scheme cannot assume add0 + p2 - p1 == add0 + p1, must subtract p2-offsets, then add p1: */\n\t"\
		"subq	%%r9,%%rax			\n\t"\
		"subq	%%r9,%%rbx			\n\t"\
		"subq	%%r9,%%rcx			\n\t"\
		"subq	%%r9,%%rdx			\n\t"\
		"\n\t"\
		/*...Block 3: Register indices for the 8 t-date = [t-index - 16]: */\
		"movslq	%[__p1],%%r8		\n\t"\
		"shlq	$3,%%r8				\n\t"\
		"addq	%%r8,%%rax			\n\t"\
		"addq	%%r8,%%rbx			\n\t"\
		"addq	%%r8,%%rcx			\n\t"\
		"addq	%%r8,%%rdx			\n\t"\
		"addq	$0x100,%%rdi 		/* r1 += 8 */\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r14,2)\n\t"/* ...+p8 */\
	"addq	%%r14,%%r13 		\n\t"/* p6 -> p10 */\
		"vbroadcastsd	0x088(%%rsi),%%ymm12 	\n\t"/* load __r1 into ymm12 */\
		"vbroadcastsd	0x0c8(%%rsi),%%ymm13 	\n\t"/* load __r9 into ymm13 */\
		"vbroadcastsd	0x0a8(%%rsi),%%ymm14 	\n\t"/* load __r5 into ymm14 */\
		"vbroadcastsd	0x0e8(%%rsi),%%ymm15 	\n\t"/* load __rD into ymm15 */\
		"vmovaps (%%rax),%%ymm0  \n\t vmovaps 0x020(%%rax),%%ymm1  \n\t"/* t16 =__A1r;	t17 =__A1i; */\
		"vmovaps (%%rcx),%%ymm4  \n\t vmovaps 0x020(%%rcx),%%ymm5  \n\t"/* t20 =__A9r;	t21 =__A9i; */\
		"vmovaps (%%rbx),%%ymm8  \n\t vmovaps 0x020(%%rbx),%%ymm9  \n\t"/* _a=  __A5r;	_b  =__A5i; */\
		"vmovaps (%%rdx),%%ymm6  \n\t vmovaps 0x020(%%rdx),%%ymm7  \n\t"/* t22 =__ADr;	t23 =__ADi; */\
		"vfnmadd231pd	0x020(%%rax),%%ymm12,%%ymm0 \n\t	 vfmadd231pd	(%%rax),%%ymm12,%%ymm1  	\n\t"/* FNMA231(__A1i,__r1,t16);	 FMA231(__A1r,__r1,t17); */\
		"vfnmadd231pd	0x020(%%rcx),%%ymm13,%%ymm4 \n\t	 vfmadd231pd	(%%rcx),%%ymm13,%%ymm5  	\n\t"/* FNMA231(__A9i,__r9,t20);	 FMA231(__A9r,__r9,t21); */\
		"vbroadcastsd	0x148(%%rsi),%%ymm13		\n\t"/* load __cD5 */\
		"vfnmadd231pd	0x020(%%rbx),%%ymm14,%%ymm8 \n\t	 vfmadd231pd	(%%rbx),%%ymm14,%%ymm9  	\n\t"/* FNMA231(__A5i,__r5,_a );	 FMA231(__A5r,__r5,_b ); */\
		"vbroadcastsd	0x128(%%rsi),%%ymm14		\n\t"/* load __c91 */\
		"vfnmadd231pd	0x020(%%rdx),%%ymm15,%%ymm6 \n\t	 vfmadd231pd	(%%rdx),%%ymm15,%%ymm7  	\n\t"/* FNMA231(__ADi,__rD,t22);	 FMA231(__ADr,__rD,t23); */\
		"vbroadcastsd	0x108(%%rsi),%%ymm15		\n\t"/* load __c51 */\
		"vmovaps		%%ymm8 ,%%ymm10				\n\t	vmovaps			%%ymm0,%%ymm2 				\n\t"/*	_c= _a;	t18= t16; */\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r13)\n\t"/* ...+p10 */\
		" vfmadd231pd	%%ymm6 ,%%ymm13,%%ymm8 		\n\t	 vfmadd231pd	%%ymm4 ,%%ymm14,%%ymm0 		\n\t"/*	 FMA231(t22,__cD5,_a);		 FMA231(t20,__c91,t16); */\
		"vmovaps		%%ymm9 ,%%ymm11				\n\t	vmovaps			%%ymm1 ,%%ymm3 				\n\t"/*	_d= _b;	t19= t17; */\
		" vfmadd231pd	%%ymm7 ,%%ymm13,%%ymm9 		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm14,%%ymm1 		\n\t"/*	 FMA231(t23,__cD5,_b);		 FMA231(t21,__c91,t17); */\
		"vfnmadd231pd	%%ymm6 ,%%ymm13,%%ymm10		\n\t	vfnmadd231pd	%%ymm4 ,%%ymm14,%%ymm2 		\n\t"/*	FNMA231(t22,__cD5,_c);		FNMA231(t20,__c91,t18); */\
		"vfnmadd231pd	%%ymm7 ,%%ymm13,%%ymm11		\n\t	vfnmadd231pd	%%ymm5 ,%%ymm14,%%ymm3 		\n\t"/*	FNMA231(t23,__cD5,_d);		FNMA231(t21,__c91,t19); */\
		"vmovaps		%%ymm0 ,%%ymm4 				\n\t	vmovaps			%%ymm1 ,%%ymm5 				\n\t"/*	t20 =t16;	t21 =t17; */\
		"vfnmadd231pd	%%ymm8 ,%%ymm15,%%ymm4 		\n\t	 vfmadd231pd	%%ymm8 ,%%ymm15,%%ymm0 		\n\t"/*	FNMA231(_a,__c51,t20);		 FMA231(_a,__c51,t16); */\
		"vfnmadd231pd	%%ymm9 ,%%ymm15,%%ymm5 		\n\t	 vfmadd231pd	%%ymm9 ,%%ymm15,%%ymm1 		\n\t"/*	FNMA231(_b,__c51,t21);		 FMA231(_b,__c51,t17); */\
		"vmovaps		%%ymm2 ,%%ymm6 				\n\t	vmovaps			%%ymm3 ,%%ymm7 				\n\t"/*	t22 =t18;	t23 =t19; */\
		" vfmadd231pd	%%ymm11,%%ymm15,%%ymm6 		\n\t	vfnmadd231pd	%%ymm11,%%ymm15,%%ymm2 		\n\t"/*	 FMA231(_d,__c51,t22);		FNMA231(_d,__c51,t18); */\
		"vfnmadd231pd	%%ymm10,%%ymm15,%%ymm7 		\n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm3 		\n\t"/*	FNMA231(_c,__c51,t23);		 FMA231(_c,__c51,t19); */\
		"vmovaps		%%ymm4 ,0x080(%%rdi)		\n\t	vmovaps			%%ymm0 ,     (%%rdi)		\n\t"/* Write outputs into local store */\
		"vmovaps		%%ymm5 ,0x0a0(%%rdi)		\n\t	vmovaps			%%ymm1 ,0x020(%%rdi)		\n\t"\
		"vmovaps		%%ymm6 ,0x0c0(%%rdi)		\n\t	vmovaps			%%ymm2 ,0x040(%%rdi)		\n\t"\
		"vmovaps		%%ymm7 ,0x0e0(%%rdi)		\n\t	vmovaps			%%ymm3 ,0x060(%%rdi)		\n\t"\
		"\n\t"\
		/*...Block 4: Register indices for the 8 t-date = [t-index - 24]: */\
		"addq	%%r9,%%rax			\n\t"/* __p2 << 3 still in %%r9: */\
		"addq	%%r9,%%rbx			\n\t"\
		"addq	%%r9,%%rcx			\n\t"\
		"addq	%%r9,%%rdx			\n\t"\
		"addq	$0x100,%%rdi 		/* r1 += 8 */\n\t"\
	"addq	%%r14,%%r13 		\n\t"/* p10 -> p14 */\
	"leaq	(%%r14,%%r14,2),%%r14 \n\t"/* p4  -> p12 */\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r14)\n\t"/* ...+p12 */\
		"vbroadcastsd	0x098(%%rsi),%%ymm12 	\n\t"/* load __r3 into ymm12 */\
		"vbroadcastsd	0x0d8(%%rsi),%%ymm13 	\n\t"/* load __rB into ymm13 */\
		"vbroadcastsd	0x0b8(%%rsi),%%ymm14 	\n\t"/* load __r7 into ymm14 */\
		"vbroadcastsd	0x0f8(%%rsi),%%ymm15 	\n\t"/* load __rF into ymm15 */\
		"vmovaps (%%rax),%%ymm0  \n\t vmovaps 0x020(%%rax),%%ymm1  \n\t"/* t24 =__A3r;	t25 =__A3i; */\
		"vmovaps (%%rcx),%%ymm4  \n\t vmovaps 0x020(%%rcx),%%ymm5  \n\t"/* t28 =__ABr;	t29 =__ABi; */\
		"vmovaps (%%rbx),%%ymm8  \n\t vmovaps 0x020(%%rbx),%%ymm9  \n\t"/* _a  =__A7r;	_b  =__A7i; */\
		"vmovaps (%%rdx),%%ymm6  \n\t vmovaps 0x020(%%rdx),%%ymm7  \n\t"/* t30 =__AFr;	t31 =__AFi; */\
		"vfnmadd231pd	0x020(%%rax),%%ymm12,%%ymm0 \n\t	 vfmadd231pd	(%%rax),%%ymm12,%%ymm1  	\n\t"/* FNMA231(__A3i,__r3,t24);	 FMA231(__A3r,__r3,t25); */\
		"vfnmadd231pd	0x020(%%rcx),%%ymm13,%%ymm4 \n\t	 vfmadd231pd	(%%rcx),%%ymm13,%%ymm5  	\n\t"/* FNMA231(__ABi,__rB,t28 );	 FMA231(__ABr,__rB,t29 ); */\
		"vbroadcastsd	0x158(%%rsi),%%ymm13		\n\t"/* load __cF7 */\
		"vfnmadd231pd	0x020(%%rbx),%%ymm14,%%ymm8 \n\t	 vfmadd231pd	(%%rbx),%%ymm14,%%ymm9  	\n\t"/* FNMA231(__A7i,__r7,_a);		 FMA231(__A7r,__r7,_b); */\
		"vbroadcastsd	0x138(%%rsi),%%ymm14		\n\t"/* load __cB3 */\
		"vfnmadd231pd	0x020(%%rdx),%%ymm15,%%ymm6 \n\t	 vfmadd231pd	(%%rdx),%%ymm15,%%ymm7  	\n\t"/* FNMA231(__AFi,__rF,t30 );	 FMA231(__AFr,__rF,t31 ); */\
		"vbroadcastsd	0x120(%%rsi),%%ymm15		\n\t"/* load __c73 */\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r13)\n\t"/* ...+p14 */\
		"vmovaps		%%ymm8 ,%%ymm10				\n\t	vmovaps			%%ymm0,%%ymm2 				\n\t"/*	_c= _a;	t26= t24; */\
		" vfmadd231pd	%%ymm6 ,%%ymm13,%%ymm8 		\n\t	 vfmadd231pd	%%ymm4 ,%%ymm14,%%ymm0 		\n\t"/*	 FMA231(t30,__cF7,_a);		 FMA231(t28,__cB3,t24); */\
		"vmovaps		%%ymm9 ,%%ymm11				\n\t	vmovaps			%%ymm1 ,%%ymm3 				\n\t"/*	_d= _b;	t27= t25; */\
		" vfmadd231pd	%%ymm7 ,%%ymm13,%%ymm9 		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm14,%%ymm1 		\n\t"/*	 FMA231(t31,__cF7,_b);		 FMA231(t29,__cB3,t25); */\
		"vfnmadd231pd	%%ymm6 ,%%ymm13,%%ymm10		\n\t	vfnmadd231pd	%%ymm4 ,%%ymm14,%%ymm2 		\n\t"/*	FNMA231(t30,__cF7,_c);		FNMA231(t28,__cB3,t26); */\
		"vfnmadd231pd	%%ymm7 ,%%ymm13,%%ymm11		\n\t	vfnmadd231pd	%%ymm5 ,%%ymm14,%%ymm3 		\n\t"/*	FNMA231(t31,__cF7,_d);		FNMA231(t29,__cB3,t27); */\
		"vmovaps		%%ymm0 ,%%ymm4 				\n\t	vmovaps			%%ymm1 ,%%ymm5 				\n\t"/*	t28 =t24;	t29 =t25; */\
		"vfnmadd231pd	%%ymm8 ,%%ymm15,%%ymm4 		\n\t	 vfmadd231pd	%%ymm8 ,%%ymm15,%%ymm0 		\n\t"/*	FNMA231(_a,__c73,t28);		 FMA231(_a,__c73,t24); */\
		"vfnmadd231pd	%%ymm9 ,%%ymm15,%%ymm5 		\n\t	 vfmadd231pd	%%ymm9 ,%%ymm15,%%ymm1 		\n\t"/*	FNMA231(_b,__c73,t29);		 FMA231(_b,__c73,t25); */\
		"vmovaps		%%ymm2 ,%%ymm6 				\n\t	vmovaps			%%ymm3 ,%%ymm7 				\n\t"/*	t30 =t26;	t31 =t27; */\
		" vfmadd231pd	%%ymm11,%%ymm15,%%ymm6 		\n\t	vfnmadd231pd	%%ymm11,%%ymm15,%%ymm2 		\n\t"/*	 FMA231(_d,__c73,t30);		FNMA231(_d,__c73,t26); */\
		"vfnmadd231pd	%%ymm10,%%ymm15,%%ymm7 		\n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm3 		\n\t"/*	FNMA231(_c,__c73,t31);		 FMA231(_c,__c73,t27); */\
		"vmovaps		%%ymm4 ,0x080(%%rdi)		\n\t	vmovaps			%%ymm0 ,     (%%rdi)		\n\t"/* Write outputs into local store */\
		"vmovaps		%%ymm5 ,0x0a0(%%rdi)		\n\t	vmovaps			%%ymm1 ,0x020(%%rdi)		\n\t"\
		"vmovaps		%%ymm6 ,0x0c0(%%rdi)		\n\t	vmovaps			%%ymm2 ,0x040(%%rdi)		\n\t"\
		"vmovaps		%%ymm7 ,0x0e0(%%rdi)		\n\t	vmovaps			%%ymm3 ,0x060(%%rdi)		\n\t"\
		"\n\t"\
	/*************************************************************************************/\
	/*  And now do four more radix-4 transforms, including the internal twiddle factors: */\
	/*************************************************************************************/\
	"movslq	%[__p2],%%r13	\n\t"/* Restore p2,4 ptr-offsets in prep for 2nd-hald prefetches */\
	"movslq	%[__p4],%%r14	\n\t"\
	"shlq	$3,%%r13		\n\t"\
	"shlq	$3,%%r14		\n\t"\
		/* Block 1: */\
		"movq	%[__add0],%%rax		\n\t"\
		"leaq	(%%r8,%%r9),%%rdx	/* [pointer arithmetic] p2 += p1 to give p3-ptr in rdx */\n\t"\
		"leaq	(%%rax,%%r8),%%rbx	/* add0+p1 */\n\t"\
		"leaq	(%%rax,%%r9),%%rcx	/* add0+p2 */\n\t"\
		"addq	%%rax,%%rdx			/* add0+p3 */\n\t"\
		"subq	$0x300,%%rdi 		/* revert r1-ptr to base value */\n\t"\
		"\n\t"\
	"addq	%%r8,%%r10	\n\t"	/*  __p1 << 3 still in %%r8: [base-address + data-fetch-ahead index + p1] */\
	"prefetcht1	%c[__pfetch_dist](%%r10)\n\t"/* ...+p1 */\
		/*...Read t0,8,16,24 from local store ... Do the 4 Im-part FMAs first, because their results needed 1st below */\
		"vbroadcastsd	0x010(%%rsi),%%ymm14		\n\t	vbroadcastsd	0x100(%%rsi),%%ymm15		\n\t"/* load __c2,c31 into pair of regs */\
		"vmovaps		     (%%rdi),%%ymm0 		\n\t	vmovaps			0x100(%%rdi),%%ymm2 		\n\t"/*    t00;    t08; */\
		"vmovaps		0x200(%%rdi),%%ymm4 		\n\t	vmovaps			0x300(%%rdi),%%ymm6 		\n\t"/*    t16;    t24; */\
		"vmovaps			 %%ymm0 ,%%ymm8 		\n\t	vmovaps				 %%ymm4 ,%%ymm10		\n\t"/* _a=t00; _c=t16; */\
		" vfmadd231pd	%%ymm2 ,%%ymm14,%%ymm0 		\n\t	 vfmadd231pd	%%ymm6 ,%%ymm15,%%ymm4 		\n\t"/*	 FMA231(t08,__c2 ,t00);		 FMA231(t24,__c31,t16); */\
		"vmovaps		0x020(%%rdi),%%ymm1 		\n\t	vmovaps			0x120(%%rdi),%%ymm3 		\n\t"/*    t01;    t09; */\
		"vmovaps		0x220(%%rdi),%%ymm5 		\n\t	vmovaps			0x320(%%rdi),%%ymm7 		\n\t"/*    t17;    t25; */\
		"vmovaps			 %%ymm1 ,%%ymm9 		\n\t	vmovaps				 %%ymm5 ,%%ymm11		\n\t"/* _b=t01; _d=t17; */\
		" vfmadd231pd	%%ymm3 ,%%ymm14,%%ymm1 		\n\t	 vfmadd231pd	%%ymm7 ,%%ymm15,%%ymm5 		\n\t"/*	 FMA231(t09,__c2 ,t01);		 FMA231(t25,__c31,t17); */\
		"vfnmadd231pd	%%ymm2 ,%%ymm14,%%ymm8 		\n\t	vfnmadd231pd	%%ymm6 ,%%ymm15,%%ymm10		\n\t"/*	FNMA231(t08,__c2 ,_a );		FNMA231(t24,__c31,_c ); */\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r13)\n\t"/* ...+p3 */\
		"vbroadcastsd	0x008(%%rsi),%%ymm6 		\n\t"/* load __c1 */\
		"vfnmadd231pd	%%ymm3 ,%%ymm14,%%ymm9 		\n\t	vfnmadd231pd	%%ymm7 ,%%ymm15,%%ymm11		\n\t"/*	FNMA231(t09,__c2 ,_b );		FNMA231(t25,__c31,_d ); */\
		"vmovaps			 %%ymm0 ,%%ymm12		\n\t	vmovaps				 %%ymm1 ,%%ymm13		\n\t"/* _e = t00; _f = t01; */\
		" vfmadd231pd	%%ymm4 ,%%ymm6 ,%%ymm0 		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm6 ,%%ymm1 		\n\t"/*	 FMA231(t16,__c1 ,t00);		 FMA231(t17,__c1 ,t01); */\
		"vfnmadd231pd	%%ymm4 ,%%ymm6 ,%%ymm12		\n\t	vfnmadd231pd	%%ymm5 ,%%ymm6 ,%%ymm13		\n\t"/*	FNMA231(t16,__c1 ,_e );		FNMA231(t17,__c1 ,_f ); */\
		"vmovaps			 %%ymm8 ,%%ymm2 		\n\t	vmovaps				 %%ymm9 ,%%ymm3 		\n\t"/* t08 = _a ; t09 = _b; */\
		"vfnmadd231pd	%%ymm11,%%ymm6 ,%%ymm2 		\n\t	 vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm3 		\n\t"/*	FNMA231(_d ,__c1 ,t08);		 FMA231(_c ,__c1 ,t09); */\
		" vfmadd231pd	%%ymm11,%%ymm6 ,%%ymm8 		\n\t	vfnmadd231pd	%%ymm10,%%ymm6 ,%%ymm9 		\n\t"/*	 FMA231(_d ,__c1 ,_a );		FNMA231(_c ,__c1 ,_b ); */\
		/* Write outputs back to main array: */\
		"vmovaps		%%ymm0 ,     (%%rax)		\n\t	vmovaps			%%ymm1 ,0x020(%%rax)		\n\t"/* __B0r= t00;		__B0i= t01; */\
		"vmovaps		%%ymm12,     (%%rbx)		\n\t	vmovaps			%%ymm13,0x020(%%rbx)		\n\t"/* __B1r= _e ;		__B1i= _f ; */\
		"vmovaps		%%ymm2 ,     (%%rcx)		\n\t	vmovaps			%%ymm3 ,0x020(%%rcx)		\n\t"/* __B2r= t08;		__B2i= t09; */\
		"vmovaps		%%ymm8 ,     (%%rdx)		\n\t	vmovaps			%%ymm9 ,0x020(%%rdx)		\n\t"/* __B3r= _a ;		__B3i= _b ; */\
		"\n\t"\
		/*...Block 3: t4,12,20,28 */\
		"vbroadcastsd	(%%rsi),%%ymm13	\n\t"/* cc0 holds 1.0 in AVX2/DFT_V2 mode */\
		"movslq	%[__p4],%%r8		\n\t"\
		"shlq	$3,%%r8				\n\t"\
		"addq	%%r8,%%rax		/* add0+p4 */\n\t"\
		"addq	%%r8,%%rbx		/* add0+p5 */\n\t"\
		"addq	%%r8,%%rcx		/* add0+p6 */\n\t"\
		"addq	%%r8,%%rdx		/* add0+p7 */\n\t"\
		"addq	$0x80,%%rdi		/* r5 */\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r14)\n\t"/* ...+p5 */\
	"addq	%%r14,%%r13 		\n\t"/* p2 -> p6 */\
		"vbroadcastsd	0x010(%%rsi),%%ymm14		\n\t	vbroadcastsd	0x100(%%rsi),%%ymm15		\n\t"/* load __c2,31 into pair of regs */\
		"vmovaps		0x200(%%rdi),%%ymm4 		\n\t	vmovaps			0x220(%%rdi),%%ymm5 		\n\t"/*    t20;    t21; */\
		"vmovaps		0x300(%%rdi),%%ymm6 		\n\t	vmovaps			0x320(%%rdi),%%ymm7 		\n\t"/*    t28;    t29; */\
		"vmovaps			 %%ymm4 ,%%ymm10 		\n\t	vmovaps				 %%ymm7 ,%%ymm11		\n\t"/* _c=t20; _d=t29; */\
		"vfnmadd231pd	%%ymm5 ,%%ymm13,%%ymm10		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm13,%%ymm4 		\n\t"/*	FNMA231(t21,1.0,_c );		 FMA231(t21,1.0,t20); */\
		" vfmadd231pd	%%ymm6 ,%%ymm13,%%ymm11		\n\t	vfnmadd231pd	%%ymm6 ,%%ymm13,%%ymm7 		\n\t"/*	 FMA231(t28,1.0,_d );		FNMA231(t28,1.0,t29); */\
		"vbroadcastsd	0x028(%%rsi),%%ymm13		\n\t"/* load __c1i2 */\
		"vmovaps		     (%%rdi),%%ymm0 		\n\t	vmovaps			0x020(%%rdi),%%ymm1 		\n\t"/*    t04;    t05; */\
		"vmovaps		0x100(%%rdi),%%ymm2 		\n\t	vmovaps			0x120(%%rdi),%%ymm3 		\n\t"/*    t12;    t13; */\
		"vmovaps			 %%ymm0 ,%%ymm8 		\n\t	vmovaps				 %%ymm1 ,%%ymm9 		\n\t"/* _a = t04; _b = t05; */\
		"vmovaps			 %%ymm10,%%ymm5 		\n\t	vmovaps				 %%ymm4 ,%%ymm12		\n\t"/* t21 = _c; _e = t20; */\
		" vfmadd231pd	%%ymm3 ,%%ymm14,%%ymm8 		\n\t	 vfmadd231pd	%%ymm11,%%ymm15,%%ymm5 		\n\t"/*	 FMA231(t13,__c2 ,_a );		 FMA231(_d ,__c31,t21); */\
		"vfnmadd231pd	%%ymm2 ,%%ymm14,%%ymm9 		\n\t	 vfmadd231pd	%%ymm7 ,%%ymm15,%%ymm4 		\n\t"/*	FNMA231(t12,__c2 ,_b );		 FMA231(t29,__c31,t20); */\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r13)\n\t"/* ...+p7 */\
		"vfnmadd231pd	%%ymm3 ,%%ymm14,%%ymm0 		\n\t	vfnmadd231pd	%%ymm11,%%ymm15,%%ymm10		\n\t"/*	FNMA231(t13,__c2 ,t04);		FNMA231(_d ,__c31,_c ); */\
		" vfmadd231pd	%%ymm2 ,%%ymm14,%%ymm1 		\n\t	vfnmadd231pd	%%ymm7 ,%%ymm15,%%ymm12		\n\t"/*	 FMA231(t12,__c2 ,t05);		FNMA231(t29,__c31,_e ); */\
		"vmovaps			 %%ymm8 ,%%ymm2 		\n\t	vmovaps				 %%ymm9 ,%%ymm3 		\n\t"/* t12 = _a; t13 = _b; */\
		"vfnmadd231pd	%%ymm4 ,%%ymm13,%%ymm2 		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm13,%%ymm3 		\n\t"/*	FNMA231(t20,__c1i2,t12);	 FMA231(t21,__c1i2,t13); */\
		" vfmadd231pd	%%ymm4 ,%%ymm13,%%ymm8 		\n\t	vfnmadd231pd	%%ymm5 ,%%ymm13,%%ymm9 		\n\t"/*	 FMA231(t20,__c1i2,_a );	FNMA231(t21,__c1i2,_b ); */\
		"vmovaps			 %%ymm0 ,%%ymm11		\n\t	vmovaps				 %%ymm1 ,%%ymm7 		\n\t"/* _d = t04; t29 = t05; */\
		" vfmadd231pd	%%ymm10,%%ymm13,%%ymm0 		\n\t	 vfmadd231pd	%%ymm12,%%ymm13,%%ymm1 		\n\t"/*	 FMA231(_c ,__c1i2,t04);	 FMA231(_e ,__c1i2,t05); */\
		"vfnmadd231pd	%%ymm10,%%ymm13,%%ymm11		\n\t	vfnmadd231pd	%%ymm12,%%ymm13,%%ymm7 		\n\t"/*	FNMA231(_c ,__c1i2,_d );	FNMA231(_e ,__c1i2,t29); */\
		/* Write outputs back to main array: */\
		"vmovaps		%%ymm2 ,     (%%rcx)		\n\t	vmovaps			%%ymm3 ,0x020(%%rcx)		\n\t"/* __B6r= t12;		__B6i= t13; */\
		"vmovaps		%%ymm8 ,     (%%rdx)		\n\t	vmovaps			%%ymm9 ,0x020(%%rdx)		\n\t"/* __B7r= _a ;		__B7i= _b ; */\
		"vmovaps		%%ymm0 ,     (%%rax)		\n\t	vmovaps			%%ymm1 ,0x020(%%rax)		\n\t"/* __B4r= t04;		__B4i= t05; */\
		"vmovaps		%%ymm11,     (%%rbx)		\n\t	vmovaps			%%ymm7 ,0x020(%%rbx)		\n\t"/* __B5r= _d ;		__B5i= t29; */\
		"\n\t"\
		/*...Block 2: t2,10,18,26 */\
		"subq	$0x40,%%rdi		/* r3 */\n\t"\
		"vmovaps		0x100(%%rdi),%%ymm2 		\n\t	vmovaps			0x120(%%rdi),%%ymm3 		\n\t"/*    t10;    t11; */\
		"addq	%%r8,%%rax		/* add0+p8 */\n\t"\
		"addq	%%r8,%%rbx		/* add0+p9 */\n\t"\
		"addq	%%r8,%%rcx		/* add0+pA */\n\t"\
		"addq	%%r8,%%rdx		/* add0+pB */\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r14,2)\n\t"/* ...+p9 */\
	"addq	%%r14,%%r13 		\n\t"/* p6 -> p10 */\
		"vbroadcastsd	0x080(%%rsi),%%ymm15		\n\t"/* load __sc  */\
		"vbroadcastsd	0x100(%%rsi),%%ymm14		\n\t"/* load __c31 */\
		"vmovaps		0x200(%%rdi),%%ymm4 		\n\t	vmovaps			0x220(%%rdi),%%ymm5 		\n\t"/*    t18;    t19; */\
		"vmovaps		0x300(%%rdi),%%ymm6 		\n\t	vmovaps			0x320(%%rdi),%%ymm7 		\n\t"/*    t26;    t27; */\
		"vsubpd		%%ymm3 ,%%ymm2 ,%%ymm12 		\n\t"/* _e = t10-t11; */\
		"vaddpd		%%ymm2 ,%%ymm3 ,%%ymm13			\n\t"/* _f = t10+t11; */\
		"vmovaps			 %%ymm4 ,%%ymm10 		\n\t	vmovaps				 %%ymm7 ,%%ymm8 		\n\t"/* _c = t18; _a = t27; */\
		"vfnmadd231pd	%%ymm5 ,%%ymm15,%%ymm10		\n\t	 vfmsub231pd	%%ymm6 ,%%ymm15,%%ymm8 		\n\t"/*	FNMA231(t19,__sc,_c );		 FMS231(t26,__sc,_a ); */\
		"vmovaps			 %%ymm5 ,%%ymm11 		\n\t	vmovaps				 %%ymm6 ,%%ymm9 		\n\t"/* _d = t19; _b = t26; */\
		"vbroadcastsd	0x030(%%rsi),%%ymm6 		\n\t"/* load __c2i2 */\
		" vfmadd231pd	%%ymm4 ,%%ymm15,%%ymm11		\n\t	 vfmadd231pd	%%ymm7 ,%%ymm15,%%ymm9 		\n\t"/*	 FMA231(t18,__sc,_d );		 FMA231(t27,__sc,_b ); */\
		"vbroadcastsd	0x018(%%rsi),%%ymm15		\n\t"/* load __c1_c */\
		"vmovaps		     (%%rdi),%%ymm0 		\n\t	vmovaps			0x020(%%rdi),%%ymm1 		\n\t"/*    t02;    t03; */\
		"vmovaps			 %%ymm10,%%ymm4 		\n\t	vmovaps				 %%ymm0 ,%%ymm2 		\n\t"/* t18 = _c;	t10 = t02; */\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r13)\n\t"/* ...+p11 */\
		" vfmadd231pd	%%ymm8 ,%%ymm14,%%ymm4 		\n\t	 vfmadd231pd	%%ymm12,%%ymm6 ,%%ymm0 		\n\t"/*	 FMA231(_a ,__c31,t18);		 FMA231(_e ,__c2i2,t02); */\
		"vmovaps			 %%ymm11,%%ymm5 		\n\t	vmovaps				 %%ymm1 ,%%ymm3 		\n\t"/* t19 = _d;	t11 = t03; */\
		" vfmadd231pd	%%ymm9 ,%%ymm14,%%ymm5 		\n\t	 vfmadd231pd	%%ymm13,%%ymm6 ,%%ymm1 		\n\t"/*	 FMA231(_b ,__c31,t19);		 FMA231(_f ,__c2i2,t03); */\
		"vfnmadd231pd	%%ymm8 ,%%ymm14,%%ymm10		\n\t	vfnmadd231pd	%%ymm12,%%ymm6 ,%%ymm2 		\n\t"/*	FNMA231(_a ,__c31,_c );		FNMA231(_e ,__c2i2,t10); */\
		"vfnmadd231pd	%%ymm9 ,%%ymm14,%%ymm11		\n\t	vfnmadd231pd	%%ymm13,%%ymm6 ,%%ymm3 		\n\t"/*	FNMA231(_b ,__c31,_d );		FNMA231(_f ,__c2i2,t11); */\
		"vmovaps			 %%ymm0 ,%%ymm8 		\n\t	vmovaps				 %%ymm1 ,%%ymm9 		\n\t"/* _a = t02; _b = t03; */\
		" vfmadd231pd	%%ymm4 ,%%ymm15,%%ymm0 		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm15,%%ymm1 		\n\t"/*	 FMA231(t18,__c1_c,t02);	 FMA231(t19,__c1_c,t03); */\
		"vfnmadd231pd	%%ymm4 ,%%ymm15,%%ymm8 		\n\t	vfnmadd231pd	%%ymm5 ,%%ymm15,%%ymm9 		\n\t"/*	FNMA231(t18,__c1_c,_a );	FNMA231(t19,__c1_c,_b ); */\
		"vmovaps			 %%ymm2 ,%%ymm12		\n\t	vmovaps				 %%ymm3 ,%%ymm13		\n\t"/* _e = t10; _f = t11; */\
		"vfnmadd231pd	%%ymm11,%%ymm15,%%ymm2 		\n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm3 		\n\t"/*	FNMA231(_d ,__c1_c,t10);	 FMA231(_c ,__c1_c,t11); */\
		" vfmadd231pd	%%ymm11,%%ymm15,%%ymm12		\n\t	vfnmadd231pd	%%ymm10,%%ymm15,%%ymm13		\n\t"/*	 FMA231(_d ,__c1_c,_e );	FNMA231(_c ,__c1_c,_f ); */\
		/* Write outputs back to main array: */\
		"vmovaps		%%ymm0 ,     (%%rax)		\n\t	vmovaps			%%ymm1 ,0x020(%%rax)		\n\t"/* __B8r= t02;		__B8i= t03; */\
		"vmovaps		%%ymm8 ,     (%%rbx)		\n\t	vmovaps			%%ymm9 ,0x020(%%rbx)		\n\t"/* __B9r= _a ;		__B9i= _b ; */\
		"vmovaps		%%ymm2 ,     (%%rcx)		\n\t	vmovaps			%%ymm3 ,0x020(%%rcx)		\n\t"/* __BAr= t10;		__BAi= t11; */\
		"vmovaps		%%ymm12,     (%%rdx)		\n\t	vmovaps			%%ymm13,0x020(%%rdx)		\n\t"/* __BBr= _e ;		__BBi= _f ; */\
		"\n\t"\
		/*...Block 4: t6,14,22,30 */\
		"addq	$0x80,%%rdi		/* r7 */\n\t"\
		"vmovaps		0x100(%%rdi),%%ymm2 		\n\t	vmovaps			0x120(%%rdi),%%ymm3 		\n\t"/*    t14;    t15; */\
		"addq	%%r8,%%rax		/* add0+pC */\n\t"\
		"addq	%%r8,%%rbx		/* add0+pD */\n\t"\
		"addq	%%r8,%%rcx		/* add0+pE */\n\t"\
		"addq	%%r8,%%rdx		/* add0+pF */\n\t"\
	"addq	%%r14,%%r13 		\n\t"/* p10 -> p14 */\
	"leaq	(%%r14,%%r14,2),%%r14 \n\t"/* p4  -> p12 */\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r14)\n\t"/* ...+p13 */\
		"vbroadcastsd	0x080(%%rsi),%%ymm15		\n\t"/* load __sc  */\
		"vmovaps		0x200(%%rdi),%%ymm4 		\n\t	vmovaps			0x220(%%rdi),%%ymm5 		\n\t"/*    t22;    t23; */\
		"vmovaps		0x300(%%rdi),%%ymm6 		\n\t	vmovaps			0x320(%%rdi),%%ymm7 		\n\t"/*    t30;    t31; */\
		"vaddpd		%%ymm2 ,%%ymm3 ,%%ymm10 		\n\t"/* _c = t14+t15; */\
		"vsubpd		%%ymm2 ,%%ymm3 ,%%ymm11			\n\t"/* _d = t15-t14; */\
		"vbroadcastsd	0x100(%%rsi),%%ymm14		\n\t"/* load __c31 */\
		"vmovaps			 %%ymm5 ,%%ymm12 		\n\t	vmovaps				 %%ymm4 ,%%ymm13		\n\t"/* _e = t23; _f = t22;*/\
		" vfmsub231pd	%%ymm4 ,%%ymm15,%%ymm12		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm15,%%ymm13		\n\t"/*	 FMS231(t22,__sc,_e );		 FMA231(t23,__sc,_f );*/\
		"vmovaps			 %%ymm6 ,%%ymm8  		\n\t	vmovaps				 %%ymm7 ,%%ymm9 		\n\t"/* _a = t30; _b = t31; */\
		"vfnmadd231pd	%%ymm7 ,%%ymm15,%%ymm8 		\n\t	 vfmadd231pd	%%ymm6 ,%%ymm15,%%ymm9 		\n\t"/*	FNMA231(t31,__sc,_a );		 FMA231(t30,__sc,_b );*/\
		"vbroadcastsd	0x030(%%rsi),%%ymm6 		\n\t"/* load __c2i2 */\
		"vbroadcastsd	0x018(%%rsi),%%ymm15		\n\t"/* load __c1_c */\
		"vmovaps		     (%%rdi),%%ymm0 		\n\t	vmovaps			0x020(%%rdi),%%ymm1 		\n\t"/*    t06;    t07; */\
		"vmovaps			 %%ymm1 ,%%ymm3 		\n\t	vmovaps				 %%ymm0 ,%%ymm2 		\n\t"/* t15= t07;	t14= t06; */\
		"vfnmadd231pd	%%ymm11,%%ymm6 ,%%ymm1 		\n\t	vfnmadd231pd	%%ymm10,%%ymm6 ,%%ymm0 		\n\t"/*	FNMA231(_d ,__c2i2,t07);	FNMA231(_c ,__c2i2,t06); */\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r13)\n\t"/* ...+p15 */\
		"vmovaps			 %%ymm12,%%ymm4 		\n\t	vmovaps				 %%ymm13,%%ymm5 		\n\t"/* t22= _e; t23= _f; */\
		"vfnmadd231pd	%%ymm8 ,%%ymm14,%%ymm4 		\n\t	vfnmadd231pd	%%ymm9 ,%%ymm14,%%ymm5 		\n\t"/*	FNMA231(_a ,__c31 ,t22);	FNMA231(_b ,__c31 ,t23); */\
		" vfmadd231pd	%%ymm8 ,%%ymm14,%%ymm12		\n\t	 vfmadd231pd	%%ymm9 ,%%ymm14,%%ymm13		\n\t"/*	 FMA231(_a ,__c31 ,_e );	 FMA231(_b ,__c31 ,_f ); */\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm2 		\n\t	 vfmadd231pd	%%ymm11,%%ymm6 ,%%ymm3 		\n\t"/*	 FMA231(_c ,__c2i2,t14);	 FMA231(_d ,__c2i2,t15); */\
		"vmovaps			 %%ymm0 ,%%ymm8 		\n\t	vmovaps				 %%ymm1 ,%%ymm9 		\n\t"/* _a = t06; _b = t07; */\
		" vfmadd231pd	%%ymm4 ,%%ymm15,%%ymm0 		\n\t	 vfmadd231pd	%%ymm5 ,%%ymm15,%%ymm1 		\n\t"/*	 FMA231(t22,__c1_c,t06);	 FMA231(t23,__c1_c,t07); */\
		"vfnmadd231pd	%%ymm4 ,%%ymm15,%%ymm8 		\n\t	vfnmadd231pd	%%ymm5 ,%%ymm15,%%ymm9 		\n\t"/*	FNMA231(t22,__c1_c,_a );	FNMA231(t23,__c1_c,_b ); */\
		"vmovaps			 %%ymm2 ,%%ymm10		\n\t	vmovaps				 %%ymm3 ,%%ymm11		\n\t"/* _c = t14; _d = t15; */\
		"vfnmadd231pd	%%ymm13,%%ymm15,%%ymm2 		\n\t	 vfmadd231pd	%%ymm12,%%ymm15,%%ymm3 		\n\t"/*	FNMA231(_f ,__c1_c,t14);	 FMA231(_e ,__c1_c,t15); */\
		" vfmadd231pd	%%ymm13,%%ymm15,%%ymm10		\n\t	vfnmadd231pd	%%ymm12,%%ymm15,%%ymm11		\n\t"/*	 FMA231(_f ,__c1_c,_c );	FNMA231(_e ,__c1_c,_d ); */\
		/* Write outputs back to main array: */\
		"vmovaps		%%ymm0 ,     (%%rax)		\n\t	vmovaps			%%ymm1 ,0x020(%%rax)		\n\t"/* __BCr= t06;		__BCi= t07; */\
		"vmovaps		%%ymm8 ,     (%%rbx)		\n\t	vmovaps			%%ymm9 ,0x020(%%rbx)		\n\t"/* __BDr= _a ;		__BDi= _b ; */\
		"vmovaps		%%ymm2 ,     (%%rcx)		\n\t	vmovaps			%%ymm3 ,0x020(%%rcx)		\n\t"/* __BEr= t14;		__BEi= t15; */\
		"vmovaps		%%ymm10,     (%%rdx)		\n\t	vmovaps			%%ymm11,0x020(%%rdx)		\n\t"/* __BFr= _c ;		__BFi= _d ; */\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r1] "m" (Xr1)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	#define RADIX16_COMPUTE_FMA_SINCOS_DIT(Xadd0,Xone)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax			\n\t"/* base-ptr to 32-byte-aligned 128-byte-large memchunk 16 doubles needing iterative inversion */\
		"movq	%[__one] ,%%rbx			\n\t"/* 1.0 in 4-fold-double form */\
		"vmovaps	    (%%rbx),%%ymm14	\n\t"/* 1.0 */\
		"vmovaps	    (%%rax),%%ymm4	\n\t"/* load inputs-to-be-inverted (d's) into 4 AVX registers */\
		"vmovaps	0x20(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x40(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x60(%%rax),%%ymm7	\n\t"\
		"vaddpd	%%ymm14,%%ymm14,%%ymm14	\n\t"/* 2.0 */\
		"vcvtpd2ps	%%ymm4,%%xmm0	\n\t"/* convert d's to SP ... Note in AVX mode output register *must* be a 128-bit xmm! */\
		"vcvtpd2ps	%%ymm5,%%xmm1	\n\t"\
		"vcvtpd2ps	%%ymm6,%%xmm2	\n\t"\
		"vcvtpd2ps	%%ymm7,%%xmm3	\n\t"\
		"vrcpps		%%xmm0,%%xmm0	\n\t"	/* ainv := approx 1/d to 11-12 bits of precision */\
		"vrcpps		%%xmm1,%%xmm1	\n\t"\
		"vrcpps		%%xmm2,%%xmm2	\n\t"\
		"vrcpps		%%xmm3,%%xmm3	\n\t"\
		"vcvtps2pd	%%xmm0,%%ymm0	\n\t"	/* convert ~1/d back to DP  ... Note in AVX mode input register *must* be a 128-bit xmm!*/\
		"vcvtps2pd	%%xmm1,%%ymm1	\n\t"\
		"vcvtps2pd	%%xmm2,%%ymm2	\n\t"\
		"vcvtps2pd	%%xmm3,%%ymm3	\n\t"\
		/* 1st NR iteration gives ~23 bits of precision: */\
		"vmovaps	%%ymm0,%%ymm8 	\n\t"	/* make a copy of ainv */\
		"vmovaps	%%ymm1,%%ymm9 	\n\t"\
		"vmovaps	%%ymm2,%%ymm10	\n\t"\
		"vmovaps	%%ymm3,%%ymm11	\n\t"\
		"vfnmadd132pd	%%ymm4,%%ymm14,%%ymm0 	\n\t"	/* 2 - d*ainv, overwrites ainv */\
		"vfnmadd132pd	%%ymm5,%%ymm14,%%ymm1 	\n\t"\
		"vfnmadd132pd	%%ymm6,%%ymm14,%%ymm2 	\n\t"\
		"vfnmadd132pd	%%ymm7,%%ymm14,%%ymm3 	\n\t"\
		"vmulpd		%%ymm0,%%ymm8 ,%%ymm0	\n\t"	/* ainv*(2 - d*ainv) = 1/d accurate to ~23 bits */\
		"vmulpd		%%ymm1,%%ymm9 ,%%ymm1	\n\t"\
		"vmulpd		%%ymm2,%%ymm10,%%ymm2	\n\t"\
		"vmulpd		%%ymm3,%%ymm11,%%ymm3	\n\t"\
		/* 3rd-order update of 23-bit result needs just 2 FMA, 1 SUB, 1 MUL: */\
		"vaddpd	    (%%rbx),%%ymm14,%%ymm14	\n\t"/* 3.0, needed by 3rd-order update step */\
		"vmovaps	%%ymm0,%%ymm8 	\n\t"	/* make a copy of ainv */\
		"vmovaps	%%ymm1,%%ymm9 	\n\t"\
		"vmovaps	%%ymm2,%%ymm10	\n\t"\
		"vmovaps	%%ymm3,%%ymm11	\n\t"\
		"vfnmadd132pd	%%ymm0,%%ymm14,%%ymm4 	\n\t"/* 1st FMA overwrites d data (inputs) with (3 - d*ainv) */\
		"vfnmadd132pd	%%ymm1,%%ymm14,%%ymm5 	\n\t"\
		"vfnmadd132pd	%%ymm2,%%ymm14,%%ymm6 	\n\t"\
		"vfnmadd132pd	%%ymm3,%%ymm14,%%ymm7 	\n\t"\
		"vsubpd		%%ymm14,%%ymm4,%%ymm0 	\n\t"/* Subtract 3 from (3 - d*ainv) to get -y = -d*ainv terms in ymm0-3 */\
		"vsubpd		%%ymm14,%%ymm5,%%ymm1 	\n\t"\
		"vsubpd		%%ymm14,%%ymm6,%%ymm2 	\n\t"\
		"vsubpd		%%ymm14,%%ymm7,%%ymm3 	\n\t"\
		"vfmadd132pd	%%ymm4,%%ymm14,%%ymm0 	\n\t"/* Positive-product FMA gives (3 - y*(3 - d*ainv)) in ymm0-3*/\
		"vfmadd132pd	%%ymm5,%%ymm14,%%ymm1 	\n\t"\
		"vfmadd132pd	%%ymm6,%%ymm14,%%ymm2 	\n\t"\
		"vfmadd132pd	%%ymm7,%%ymm14,%%ymm3 	\n\t"\
		"vmulpd		%%ymm0,%%ymm8 ,%%ymm0	\n\t"	/* ainv*(3 - y*(3 - d*ainv)) = 1/d accurate to ~53 bits */\
		"vmulpd		%%ymm1,%%ymm9 ,%%ymm1	\n\t"\
		"vmulpd		%%ymm2,%%ymm10,%%ymm2	\n\t"\
		"vmulpd		%%ymm3,%%ymm11,%%ymm3	\n\t"\
		/* Multiply by the sine terms to get quotients in ymm4-7: */\
		"vmulpd	0x80(%%rax),%%ymm0,%%ymm4	\n\t"\
		"vmulpd	0xa0(%%rax),%%ymm1,%%ymm5	\n\t"\
		"vmulpd	0xc0(%%rax),%%ymm2,%%ymm6	\n\t"\
		"vmulpd	0xe0(%%rax),%%ymm3,%%ymm7	\n\t"\
		/* Outputs overwrite sines: */\
		"vmovaps	%%ymm4,0x80(%%rax)	\n\t"\
		"vmovaps	%%ymm5,0xa0(%%rax)	\n\t"\
		"vmovaps	%%ymm6,0xc0(%%rax)	\n\t"\
		"vmovaps	%%ymm7,0xe0(%%rax)	\n\t"\
		/* Now propagate the scalar doubles in the 8 ymm-sized memlocs referenced above to their final 4x-copied locs. Scalar data laid out as */\
		/* add0 + 0x[  0,  8, 10, 18]: [c0 ,c1 ,c2 ,c3 ] Cosines: */\
		/* add0 + 0x[ 20, 28, 30, 38]: [c4 ,c5 ,c6 ,c7 ] */\
		/* add0 + 0x[ 40, 48, 50, 58]: [c8 ,c9 ,cA ,cB ] */\
		/* add0 + 0x[ 60, 68, 70, 78]: [cC ,cD ,cE ,cF ] */\
		/* add0 + 0x[ 80, 88, 90, 98]: [r0 ,r1 ,r2 ,r3 ] Tangents: */\
		/* add0 + 0x[ a0, a8, b0, b8]: [r4 ,r5 ,r6 ,r7 ] */\
		/* add0 + 0x[ c0, c8, d0, d8]: [r8 ,r9 ,rA ,rB ] */\
		/* add0 + 0x[ e0, e8, f0, f8]: [rC ,rD ,rE ,rF ] */\
		/* We start at the upper end of the memory chunk allocated for sincos data and proceed downward, so the 8 memory */\
		/* slots inited with scalar data above are only overwritten after their scalar-data contents have been copied: */\
		"vbroadcastsd	0x0f8(%%rax),%%ymm15	/* rF */\n\t"/* Allow several cycles for each vbroadcast to complete before writing result */\
		"vbroadcastsd	0x078(%%rax),%%ymm14	/* cF */\n\t"\
		"vbroadcastsd	0x0d8(%%rax),%%ymm13	/* rB */\n\t"\
		"vbroadcastsd	0x058(%%rax),%%ymm12	/* cB */\n\t"\
		"vbroadcastsd	0x0b8(%%rax),%%ymm11	/* r7 */\n\t	vmovaps	%%ymm15,0x420(%%rax)	\n\t"/* s15 = cc0 + 0x21;	__rF = s15/c15		*/\
		"vbroadcastsd	0x038(%%rax),%%ymm10	/* c7 */\n\t	vmovaps	%%ymm14,0x400(%%rax)	\n\t"/* c15 = cc0 + 0x20;	__cF				*/\
		"vbroadcastsd	0x098(%%rax),%%ymm9 	/* r3 */\n\t	vmovaps	%%ymm13,0x3e0(%%rax)	\n\t"/* s11 = cc0 + 0x1f;	__rB = s11/c11		*/\
		"vbroadcastsd	0x018(%%rax),%%ymm8 	/* c3 */\n\t	vmovaps	%%ymm12,0x3c0(%%rax)	\n\t"/* c11 = cc0 + 0x1e;	__cB				*/\
		"vbroadcastsd	0x0e8(%%rax),%%ymm7 	/* rD */\n\t	vmovaps	%%ymm11,0x3a0(%%rax)	\n\t"/* s7  = cc0 + 0x1d;	__r7 = s7 /c7 		*/\
		"vbroadcastsd	0x068(%%rax),%%ymm6 	/* cD */\n\t	vmovaps	%%ymm10,0x380(%%rax)	\n\t"/* c7  = cc0 + 0x1c;	__c7				*/\
		"vbroadcastsd	0x0c8(%%rax),%%ymm5 	/* r9 */\n\t	vmovaps	%%ymm9 ,0x360(%%rax)	\n\t"/* s3  = cc0 + 0x1b;	__r3 = s3 /c3 		*/\
		"vbroadcastsd	0x048(%%rax),%%ymm4 	/* c9 */\n\t	vmovaps	%%ymm8 ,0x340(%%rax)	\n\t"/* c3  = cc0 + 0x1a;	__c3				*/\
		"vbroadcastsd	0x0a8(%%rax),%%ymm3 	/* r5 */\n\t	vmovaps	%%ymm7 ,0x320(%%rax)	\n\t"/* s13 = cc0 + 0x19;	__rD = s13/c13		*/\
		"vbroadcastsd	0x028(%%rax),%%ymm2 	/* c5 */\n\t	vmovaps	%%ymm6 ,0x300(%%rax)	\n\t"/* c13 = cc0 + 0x18;	__cD				*/\
		"vbroadcastsd	0x088(%%rax),%%ymm1 	/* r1 */\n\t	vmovaps	%%ymm5 ,0x2e0(%%rax)	\n\t"/* s9  = cc0 + 0x17;	__r9 = s9 /c9 		*/\
		"vbroadcastsd	0x008(%%rax),%%ymm0 	/* c1 */\n\t	vmovaps	%%ymm4 ,0x2c0(%%rax)	\n\t"/* c9  = cc0 + 0x16;	__c9				*/\
		"vbroadcastsd	0x0f0(%%rax),%%ymm15	/* rE */\n\t	vmovaps	%%ymm3 ,0x2a0(%%rax)	\n\t"/* s5  = cc0 + 0x15;	__r5 = s5 /c5 		*/\
		"vbroadcastsd	0x070(%%rax),%%ymm14	/* cE */\n\t	vmovaps	%%ymm2 ,0x280(%%rax)	\n\t"/* c5  = cc0 + 0x14;	__c5				*/\
		"vbroadcastsd	0x0d0(%%rax),%%ymm13	/* rA */\n\t	vmovaps	%%ymm1 ,0x260(%%rax)	\n\t"/* s1  = cc0 + 0x13;	__r1 = s1 /c1 		*/\
		"vbroadcastsd	0x050(%%rax),%%ymm12	/* cA */\n\t	vmovaps	%%ymm0 ,0x240(%%rax)	\n\t"/* s1  = cc0 + 0x12;	__c1				*/\
		"vbroadcastsd	0x0b0(%%rax),%%ymm11	/* r6 */\n\t	vmovaps	%%ymm15,0x220(%%rax)	\n\t"/* s14 = cc0 + 0x11;	__rE = s14/c14		*/\
		"vbroadcastsd	0x030(%%rax),%%ymm10	/* c6 */\n\t	vmovaps	%%ymm14,0x200(%%rax)	\n\t"/* c14 = cc0 + 0x10;	__cE				*/\
		"vbroadcastsd	0x090(%%rax),%%ymm9 	/* r2 */\n\t	vmovaps	%%ymm13,0x1e0(%%rax)	\n\t"/* s10 = cc0 + 0x0f;	__rA = s10/c10		*/\
		"vbroadcastsd	0x010(%%rax),%%ymm8 	/* c2 */\n\t	vmovaps	%%ymm12,0x1c0(%%rax)	\n\t"/* c10 = cc0 + 0x0e;	__cA				*/\
		"vbroadcastsd	0x0e0(%%rax),%%ymm7 	/* rC */\n\t	vmovaps	%%ymm11,0x1a0(%%rax)	\n\t"/* s6  = cc0 + 0x0d;	__r6 = s6 /c6 		*/\
		"vbroadcastsd	0x060(%%rax),%%ymm6 	/* cC */\n\t	vmovaps	%%ymm10,0x180(%%rax)	\n\t"/* c6  = cc0 + 0x0c;	__c6				*/\
		"vbroadcastsd	0x0c0(%%rax),%%ymm5 	/* r8 */\n\t	vmovaps	%%ymm9 ,0x160(%%rax)	\n\t"/* s2  = cc0 + 0x0b;	__r2 = s2 /c2 		*/\
		"vbroadcastsd	0x040(%%rax),%%ymm4 	/* c8 */\n\t	vmovaps	%%ymm8 ,0x140(%%rax)	\n\t"/* c2  = cc0 + 0x0a;	__c2				*/\
		"vbroadcastsd	0x0a0(%%rax),%%ymm3 	/* r4 */\n\t	vmovaps	%%ymm7 ,0x120(%%rax)	\n\t"/* s12 = cc0 + 0x09;	__rC = s12/c12		*/\
		"vbroadcastsd	0x020(%%rax),%%ymm2 	/* c4 */\n\t	vmovaps	%%ymm6 ,0x100(%%rax)	\n\t"/* c12 = cc0 + 0x08;	__cC				*/\
		"vbroadcastsd	0x080(%%rax),%%ymm1		/*__t */\n\t	vmovaps	%%ymm5 ,0x0e0(%%rax)	\n\t"/* s8  = cc0 + 0x07;	__r8 = s8 /c8 		*/\
		"vbroadcastsd	     (%%rax),%%ymm0		/*__c */\n\t	vmovaps	%%ymm4 ,0x0c0(%%rax)	\n\t"/* c8  = cc0 + 0x06;	__c8 [unchanged]	*/\
		"														vmovaps	%%ymm3 ,0x0a0(%%rax)	\n\t"/* s4  = cc0 + 0x05;	__r4 = s4 /c4 		*/\
		"														vmovaps	%%ymm2 ,0x080(%%rax)	\n\t"/* c4  = cc0 + 0x04;	__c4 [unchanged]	*/\
		"														vmovaps	%%ymm1 ,0x020(%%rax)	\n\t"/* ss0 = cc0 + 0x01;	__sc = __s/__c		*/\
		"														vmovaps	%%ymm0 ,     (%%rax)	\n\t"/* cc0 = cc0 + 0x00;	__c					*/\
		"vmovaps	    (%%rbx),%%ymm14	\n\t"/* 1.0 */\
		"vxorpd	%%ymm0,%%ymm0,%%ymm0	\n\t"/* 0.0 */\
		"vmovaps	%%ymm14,0x040(%%rax)	\n\t"/* c0 = cc0 + 0x02;	don't use in DFT but init = 1.0 */\
		"vmovaps	%%ymm0 ,0x060(%%rax)	\n\t"/* s0 = cc0 + 0x03;	don't use in DFT but init = 0.0 */\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__one] "m" (Xone)\
		: "cc","memory","rax","rbx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	/* Note no cosine ratios needed in DIT version. */\
	// 0 ADD, 174 FMA [most involving a unity multiplicand], 34 pure-MUL, 255 MOVAPS
	#define SSE2_RADIX16_DIT_TWIDDLE_1(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr1,Xcc0,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/* Gather the needed data and do first set of four length-4 transforms. Each complex radix-4 block needs 11 registers: */\
	/* 8 t-temps mapped to ymm0-7, plus rt,it (ymm8,9), plus const 1.0 (ymm10), though we add 2 more tmps re,im (ymm11,12) */\
	/* to allow independent Re/Im subsections to overlap w/o introducing false dependencies. */\
	/*...Block 1: */\
		"movslq	%[__p4],%%r10	\n\t"/* This serves as main-array stride-between radix-4 blocks */\
		"movq	%[__add0],%%rax				\n\t"\
	"movslq	%[__pfetch_addr1],%%r13	\n\t"	/* Prefetch base-index offset cycles among add0 + p0,1,2,3 on successive macro calls */\
	"leaq	(%%rax,%%r13,8),%%r13	\n\t"	/* [base-address + data-fetch-ahead index] */\
		"movslq	%[__p1],%%rbx	\n\t"\
		"movslq	%[__p2],%%rcx	\n\t"\
		"movslq	%[__p3],%%rdx	\n\t"\
		"movq	%[__cc0],%%rsi 	\n\t"\
		"shlq	$3,%%r10		\n\t"/* p4 in bytewise ptr-arithmetic form */\
		"vmovaps 0x440(%%rsi),%%ymm10	\n\t"/* cc0 + 0x22 = __two; Actually holds 1.0 in AVX2 mode */\
		"leaq	(%%rax,%%rbx,8),%%rbx	\n\t"/* add0+p1 */\
		"leaq	(%%rax,%%rcx,8),%%rcx	\n\t"/* add0+p2 */\
		"leaq	(%%rax,%%rdx,8),%%rdx	\n\t"/* add0+p3 */\
		"movq	%[__r1],%%rdi	\n\t"/* ptr to local-mem [t1] */\
	"movq	%%r10,%%r14	\n\t"	/* Save a copy of p4 pointer offset. Will prefetch from [base-address + data-fetch-ahead index] + [0,p4,p8,p12] on each macro call. */\
	"prefetcht1	%c[__pfetch_dist](%%r13)\n\t"/* [base-address + data-fetch-ahead index] + p0 */\
		"vmovaps		     (%%rax),%%ymm0 	\n\t	vmovaps		     (%%rbx),%%ymm8 	\n\t"/*	t1,rt =__A0,1r; */\
		"vmovaps		     (%%rcx),%%ymm4 	\n\t	vmovaps		     (%%rdx),%%ymm9 	\n\t"/*	t5,it =__A2,3r; */\
		"vmovaps		%%ymm0,%%ymm2			\n\t	vmovaps		%%ymm4,%%ymm7	\n\t"/*	t3 = t1; t8 = t5; */\
		"vfnmadd231pd	%%ymm10,%%ymm8 ,%%ymm2 	\n\t"/*	FNMA231(1.0,rt ,t3 ); */\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm0 	\n\t"/*	 FMA231(1.0,rt ,t1 ); */\
		"vfnmadd231pd	%%ymm10,%%ymm9 ,%%ymm7 	\n\t"/*	FNMA231(1.0,it ,t8 ); */\
		" vfmadd231pd	%%ymm10,%%ymm9 ,%%ymm4 	\n\t"/*	 FMA231(1.0,it ,t5 ); */\
		"vmovaps		0x020(%%rax),%%ymm1 	\n\t	vmovaps		0x020(%%rbx),%%ymm11	\n\t"/*	t2,re =__A0,1i; */\
		"vmovaps		0x020(%%rcx),%%ymm5 	\n\t	vmovaps		0x020(%%rdx),%%ymm12	\n\t"/*	t6,im =__A2,3i; */\
		"vmovaps		%%ymm1,%%ymm3			\n\t	vmovaps		%%ymm5,%%ymm6	\n\t"/*	t4 = t2; t7 = t6; */\
		"vfnmadd231pd	%%ymm10,%%ymm11,%%ymm3 	\n\t"/*	FNMA231(1.0,re ,t4 ); */\
		" vfmadd231pd	%%ymm10,%%ymm11,%%ymm1 	\n\t"/*	 FMA231(1.0,re ,t2 ); */\
		"vfnmadd231pd	%%ymm10,%%ymm12,%%ymm6 	\n\t"/*	FNMA231(1.0,im ,t7 ); */\
		" vfmadd231pd	%%ymm10,%%ymm12,%%ymm5 	\n\t"/*	 FMA231(1.0,im ,t6 ); */\
		"vmovaps		%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm2,%%ymm9	\n\t"/*	rt = t5; it = t3; */\
		/*** REMEMBER - For FMA132, rightmost 2 operands in ASM are reverse-ordered w.r.to the prototyping code in the comments! ***/\
		"vfnmadd132pd	%%ymm10,%%ymm0 ,%%ymm4 	\n\t"/*	FNMA132(1.0,t5 ,t1 ); */\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm0 	\n\t"/*	 FMA231(1.0,rt ,t1 ); */\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm2 	\n\t"/*	 FMA231(1.0,t7 ,t3 ); */\
		"vfnmadd132pd	%%ymm10,%%ymm9 ,%%ymm6 	\n\t"/*	FNMA132(1.0,t7 ,it ); */\
		"vmovaps		%%ymm5,%%ymm11			\n\t	vmovaps		%%ymm3,%%ymm12		\n\t"/*	re = t6; im = t4; */\
		"vfnmadd132pd	%%ymm10,%%ymm1 ,%%ymm5 	\n\t	vmovaps	%%ymm4 ,0x080(%%rdi)	\n\t"/* FNMA132(1.0,t6 ,t2 ); write t5 */\
		" vfmadd231pd	%%ymm10,%%ymm11,%%ymm1 	\n\t	vmovaps	%%ymm0 ,     (%%rdi)	\n\t"/*  FMA231(1.0,re ,t2 ); write t1 */\
		"vfnmadd231pd	%%ymm10,%%ymm7 ,%%ymm3 	\n\t	vmovaps	%%ymm2 ,0x040(%%rdi)	\n\t"/* FNMA231(1.0,t8 ,t4 ); write t3 */\
		" vfmadd132pd	%%ymm10,%%ymm12,%%ymm7 	\n\t	vmovaps	%%ymm6 ,0x0c0(%%rdi)	\n\t"/*  FMA132(1.0,t8 ,im ); write t7 */\
		"vmovaps		%%ymm5 ,0x0a0(%%rdi)  	\n\t	addq	%%r10,%%rax	\n\t"/* write t6; add0+p4 */\
		"vmovaps		%%ymm1 ,0x020(%%rdi)  	\n\t	addq	%%r10,%%rbx	\n\t"/* write t2; add0+p5 */\
		"vmovaps		%%ymm3 ,0x060(%%rdi)  	\n\t	addq	%%r10,%%rcx	\n\t"/* write t4; add0+p6 */\
		"vmovaps		%%ymm7 ,0x0e0(%%rdi)  	\n\t	addq	%%r10,%%rdx	\n\t"/* write t8; add0+p7 */\
		/**/\
	/*...Block 2: __A4-7, t9-16 */\
		"addq	$0x100,%%rdi	\n\t"/* ptr to t9 */\
		"vmovaps		     (%%rax),%%ymm0 	\n\t	vmovaps		     (%%rbx),%%ymm8 	\n\t"\
		"vmovaps		     (%%rcx),%%ymm4 	\n\t	vmovaps		     (%%rdx),%%ymm9 	\n\t"\
		"vmovaps		%%ymm0,%%ymm2			\n\t	vmovaps		%%ymm4,%%ymm7	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm8 ,%%ymm2 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm0 	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm9 ,%%ymm7 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm9 ,%%ymm4 	\n\t"\
		"vmovaps		0x020(%%rax),%%ymm1 	\n\t	vmovaps		0x020(%%rbx),%%ymm11	\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm5 	\n\t	vmovaps		0x020(%%rdx),%%ymm12	\n\t"\
		"vmovaps		%%ymm1,%%ymm3			\n\t	vmovaps		%%ymm5,%%ymm6	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm11,%%ymm3 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm11,%%ymm1 	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm12,%%ymm6 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm12,%%ymm5 	\n\t"\
		"vmovaps		%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm2,%%ymm9	\n\t"\
		"vfnmadd132pd	%%ymm10,%%ymm0 ,%%ymm4 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm0 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm2 	\n\t"\
		"vfnmadd132pd	%%ymm10,%%ymm9 ,%%ymm6 	\n\t"\
		"vmovaps		%%ymm5,%%ymm11			\n\t	vmovaps		%%ymm3,%%ymm12	\n\t"\
		"vfnmadd132pd	%%ymm10,%%ymm1 ,%%ymm5 	\n\t	vmovaps		%%ymm4 ,0x080(%%rdi) \n\t"\
		" vfmadd231pd	%%ymm10,%%ymm11,%%ymm1 	\n\t	vmovaps		%%ymm0 ,     (%%rdi) \n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm7 ,%%ymm3 	\n\t	vmovaps		%%ymm2 ,0x040(%%rdi) \n\t"\
		" vfmadd132pd	%%ymm10,%%ymm12,%%ymm7 	\n\t	vmovaps		%%ymm6 ,0x0c0(%%rdi) \n\t"\
		"vmovaps		%%ymm5 ,0x0a0(%%rdi)  	\n\t	addq	%%r10,%%rax	\n\t"/*add0+p8 */\
		"vmovaps		%%ymm1 ,0x020(%%rdi)  	\n\t	addq	%%r10,%%rbx	\n\t"/*add0+p9 */\
		"vmovaps		%%ymm3 ,0x060(%%rdi)  	\n\t	addq	%%r10,%%rcx	\n\t"/*add0+pA */\
		"vmovaps		%%ymm7 ,0x0e0(%%rdi)  	\n\t	addq	%%r10,%%rdx	\n\t"/*add0+pB */\
		/**/\
	/*...Block 3: __A8-B, t17-24 */\
		"addq	$0x100,%%rdi	\n\t"/* ptr to t17 */\
		"vmovaps		     (%%rax),%%ymm0 	\n\t	vmovaps		     (%%rbx),%%ymm8 	\n\t"\
		"vmovaps		     (%%rcx),%%ymm4 	\n\t	vmovaps		     (%%rdx),%%ymm9 	\n\t"\
		"vmovaps		%%ymm0,%%ymm2			\n\t	vmovaps		%%ymm4,%%ymm7	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm8 ,%%ymm2 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm0 	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm9 ,%%ymm7 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm9 ,%%ymm4 	\n\t"\
		"vmovaps		0x020(%%rax),%%ymm1 	\n\t	vmovaps		0x020(%%rbx),%%ymm11	\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm5 	\n\t	vmovaps		0x020(%%rdx),%%ymm12	\n\t"\
		"vmovaps		%%ymm1,%%ymm3			\n\t	vmovaps		%%ymm5,%%ymm6	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm11,%%ymm3 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm11,%%ymm1 	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm12,%%ymm6 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm12,%%ymm5 	\n\t"\
		"vmovaps		%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm2,%%ymm9	\n\t"\
		"vfnmadd132pd	%%ymm10,%%ymm0 ,%%ymm4 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm0 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm2 	\n\t"\
		"vfnmadd132pd	%%ymm10,%%ymm9 ,%%ymm6 	\n\t"\
		"vmovaps		%%ymm5,%%ymm11			\n\t	vmovaps		%%ymm3,%%ymm12	\n\t"\
		"vfnmadd132pd	%%ymm10,%%ymm1 ,%%ymm5 	\n\t	vmovaps		%%ymm4 ,0x080(%%rdi) \n\t"\
		" vfmadd231pd	%%ymm10,%%ymm11,%%ymm1 	\n\t	vmovaps		%%ymm0 ,     (%%rdi) \n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm7 ,%%ymm3 	\n\t	vmovaps		%%ymm2 ,0x040(%%rdi) \n\t"\
		" vfmadd132pd	%%ymm10,%%ymm12,%%ymm7 	\n\t	vmovaps		%%ymm6 ,0x0c0(%%rdi) \n\t"\
		"vmovaps		%%ymm5 ,0x0a0(%%rdi)  	\n\t	addq	%%r10,%%rax	\n\t"/*add0+pC */\
		"vmovaps		%%ymm1 ,0x020(%%rdi)  	\n\t	addq	%%r10,%%rbx	\n\t"/*add0+pD */\
		"vmovaps		%%ymm3 ,0x060(%%rdi)  	\n\t	addq	%%r10,%%rcx	\n\t"/*add0+pE */\
		"vmovaps		%%ymm7 ,0x0e0(%%rdi)  	\n\t	addq	%%r10,%%rdx	\n\t"/*add0+pF */\
		/**/\
	/*...Block 4: __AC-F, t25-32 */\
		"addq	$0x100,%%rdi	\n\t"/* ptr to t25 */\
		"vmovaps		     (%%rax),%%ymm0 	\n\t	vmovaps		     (%%rbx),%%ymm8 	\n\t"\
		"vmovaps		     (%%rcx),%%ymm4 	\n\t	vmovaps		     (%%rdx),%%ymm9 	\n\t"\
		"vmovaps		%%ymm0,%%ymm2			\n\t	vmovaps		%%ymm4,%%ymm7	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm8 ,%%ymm2 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm0 	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm9 ,%%ymm7 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm9 ,%%ymm4 	\n\t"\
		"vmovaps		0x020(%%rax),%%ymm1 	\n\t	vmovaps		0x020(%%rbx),%%ymm11	\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm5 	\n\t	vmovaps		0x020(%%rdx),%%ymm12	\n\t"\
		"vmovaps		%%ymm1,%%ymm3			\n\t	vmovaps		%%ymm5,%%ymm6	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm11,%%ymm3 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm11,%%ymm1 	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm12,%%ymm6 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm12,%%ymm5 	\n\t"\
		"vmovaps		%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm2,%%ymm9	\n\t"\
		"vfnmadd132pd	%%ymm10,%%ymm0 ,%%ymm4 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm0 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm2 	\n\t"\
		"vfnmadd132pd	%%ymm10,%%ymm9 ,%%ymm6 	\n\t"\
		"vmovaps		%%ymm5,%%ymm11			\n\t	vmovaps		%%ymm3,%%ymm12	\n\t"\
		"vfnmadd132pd	%%ymm10,%%ymm1 ,%%ymm5 	\n\t	vmovaps		%%ymm4 ,0x080(%%rdi) \n\t"\
		" vfmadd231pd	%%ymm10,%%ymm11,%%ymm1 	\n\t	vmovaps		%%ymm0 ,     (%%rdi) \n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm7 ,%%ymm3 	\n\t	vmovaps		%%ymm2 ,0x040(%%rdi) \n\t"\
		" vfmadd132pd	%%ymm10,%%ymm12,%%ymm7 	\n\t	vmovaps		%%ymm6 ,0x0c0(%%rdi) \n\t"\
		"vmovaps		%%ymm5 ,0x0a0(%%rdi)  	\n\t"\
		"vmovaps		%%ymm1 ,0x020(%%rdi)  	\n\t"\
		"vmovaps		%%ymm3 ,0x060(%%rdi)  	\n\t"\
		"vmovaps		%%ymm7 ,0x0e0(%%rdi)  	\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%r14)\n\t"/* ...+p4 */\
		"subq	$0x300,%%rdi	\n\t"/* Revert rdi to point at t1 */\
		/**/\
	/**************************************************************************************/\
	/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
	/**************************************************************************************/\
		"movq	%%r10,%%rbx	/* p4 already in ptr-offset form in r10 */\n\t"\
		"movslq	%[__p1],%%r10	\n\t"/* This serves as main-array stride-between radix-4 blocks */\
		"shlq	$3,%%r10		\n\t"\
		/**/\
	/*...Block 1: t1/2,9/10,17/18,25/26 in ymm0-7, resp. */\
		"movq	%[__add0],%%rax			\n\t"\
		"movslq	%[__p12],%%rdx			\n\t"\
		"leaq	(%%rax,%%rbx,2),%%rcx	\n\t"/* add0+p8 */\
		"addq	%%rax,%%rbx				\n\t"/* add0+p4 */\
		"leaq	(%%rax,%%rdx,8),%%rdx	\n\t"/* add0+pC */\
		"vmovaps		     (%%rdi),%%ymm0 	\n\t	vmovaps		0x100(%%rdi),%%ymm2 	\n\t"/*	t1,9 ; */\
		"vmovaps		0x020(%%rdi),%%ymm1 	\n\t	vmovaps		0x120(%%rdi),%%ymm3 	\n\t"/*	t2,10; */\
		"vmovaps		%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9	\n\t"/*	rt = t1; it = t2; */\
		" vfmadd231pd	%%ymm10,%%ymm2 ,%%ymm0 	\n\t"/*	 FMA231(1.0,t9 ,t1 ); */\
		"vfnmadd132pd	%%ymm10,%%ymm8 ,%%ymm2 	\n\t"/*	FNMA132(1.0,t9 ,rt ); */\
		" vfmadd231pd	%%ymm10,%%ymm3 ,%%ymm1 	\n\t"/*	 FMA231(1.0,t10,t2 ); */\
		"vfnmadd132pd	%%ymm10,%%ymm9 ,%%ymm3 	\n\t"/*	FNMA132(1.0,t10,it ); */\
		"vmovaps		0x200(%%rdi),%%ymm4 	\n\t	vmovaps		0x300(%%rdi),%%ymm6 	\n\t"/*	t17,25; */\
		"vmovaps		0x220(%%rdi),%%ymm5 	\n\t	vmovaps		0x320(%%rdi),%%ymm7 	\n\t"/*	t18,26; */\
		"vmovaps		%%ymm4,%%ymm11			\n\t	vmovaps		%%ymm5,%%ymm12	\n\t"/*	re = t17; im = t18; */\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm4 	\n\t"/*	 FMA231(1.0,t25,t17); */\
		"vfnmadd132pd	%%ymm10,%%ymm11,%%ymm6 	\n\t"/*	FNMA132(1.0,t25,re ); */\
		" vfmadd231pd	%%ymm10,%%ymm7 ,%%ymm5 	\n\t"/*	 FMA231(1.0,t26,t18); */\
		"vfnmadd132pd	%%ymm10,%%ymm12,%%ymm7 	\n\t"/*	FNMA132(1.0,t26,im ); */\
		"vmovaps		%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9	\n\t"/*	rt = t1; it = t2; */\
		" vfmadd231pd	%%ymm10,%%ymm4 ,%%ymm8 	\n\t"/*	 FMA231(1.0,t17,rt ); */\
		" vfmadd231pd	%%ymm10,%%ymm5 ,%%ymm9 	\n\t"/*	 FMA231(1.0,t18,it ); */\
		"vfnmadd231pd	%%ymm10,%%ymm4 ,%%ymm0 	\n\t"/*	FNMA231(1.0,t17,t1 ); */\
		"vfnmadd231pd	%%ymm10,%%ymm5 ,%%ymm1 	\n\t"/*	FNMA231(1.0,t18,t2 ); */\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%r14,2)\n\t"/* ...+p8 */\
		"vmovaps		%%ymm2,%%ymm11			\n\t	vmovaps		%%ymm3,%%ymm12	\n\t"/*	re = t9 ; im = t10; */\
		"vmovaps 0x0a0(%%rsi),%%ymm13 \n\t vmovaps 0x0e0(%%rsi),%%ymm14 \n\t vmovaps 0x120(%%rsi),%%ymm15 \n\t"/* __t4,8,C; */\
		"vfnmadd231pd	%%ymm10,%%ymm7 ,%%ymm2 	\n\t"/*	FNMA231(1.0,t26,t9 ); */\
		" vfmadd132pd	%%ymm10,%%ymm11,%%ymm7 	\n\t"/*	 FMA132(1.0,t26,re ); */\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm3 	\n\t"/*	 FMA231(1.0,t25,t10); */\
		"vfnmadd132pd	%%ymm10,%%ymm12,%%ymm6 	\n\t"/*	FNMA132(1.0,t25,im ); */\
		"vmovaps		%%ymm8 ,     (%%rax)  	\n\t"/* __B0r = rt; */\
		"vmovaps		%%ymm9 ,0x020(%%rax)  	\n\t"/* __B0i = it; */\
		"vmovaps %%ymm0,%%ymm8	\n\t	vmovaps %%ymm7,%%ymm11 \n\t vmovaps %%ymm2,%%ymm12	\n\t"/*	rt =t1; re =t26; im =t9; */\
		"vmovaps 0x080(%%rsi),%%ymm4  \n\t vmovaps 0x0c0(%%rsi),%%ymm5  \n\t vmovaps 0x100(%%rsi),%%ymm9  \n\t"/* __c4,8,C; use unused ymms for these */\
		" vfmadd231pd	%%ymm14,%%ymm1 ,%%ymm0 	\n\t"/*	 FMA231(__t8,t2 ,t1 ); */\
		"vfnmadd231pd	%%ymm14,%%ymm8 ,%%ymm1 	\n\t"/*	FNMA231(__t8,rt ,t2 ); */\
		" vfmadd231pd	%%ymm13,%%ymm6 ,%%ymm7 	\n\t"/*	 FMA231(__t4,t25,t26); */\
		"vfnmadd231pd	%%ymm13,%%ymm11,%%ymm6 	\n\t"/*	FNMA231(__t4,re ,t25); */\
		" vfmadd231pd	%%ymm15,%%ymm3 ,%%ymm2 	\n\t"/*	 FMA231(__tC,t10,t9 ); */\
		"vfnmadd231pd	%%ymm15,%%ymm12,%%ymm3 	\n\t"/*	FNMA231(__tC,im ,t10); */\
		"vmulpd			%%ymm5 ,%%ymm0 ,%%ymm0 	\n\t"/*	t1  *= __c8; */\
		"vmulpd			%%ymm5 ,%%ymm1 ,%%ymm1 	\n\t"/*	t2  *= __c8; */\
		"vmulpd			%%ymm4 ,%%ymm7 ,%%ymm7 	\n\t"/*	t26 *= __c4; */\
		"vmulpd			%%ymm4 ,%%ymm6 ,%%ymm6 	\n\t"/*	t25 *= __c4; */\
		"vmulpd			%%ymm9 ,%%ymm2 ,%%ymm2 	\n\t"/*	t9  *= __cC; */\
		"vmulpd			%%ymm9 ,%%ymm3 ,%%ymm3 	\n\t"/*	t10 *= __cC; */\
		"vmovaps		%%ymm0 ,     (%%rcx)  	\n\t"/* __B8r = t1 ; */\
		"vmovaps		%%ymm1 ,0x020(%%rcx)  	\n\t"/* __B8i = t2 ; */\
		"vmovaps		%%ymm7 ,     (%%rbx)  	\n\t"/* __B4r = t26; */\
		"vmovaps		%%ymm6 ,0x020(%%rbx)  	\n\t"/* __B4i = t25; */\
		"vmovaps		%%ymm2 ,     (%%rdx)  	\n\t"/* __BCr = t9 ; */\
		"vmovaps		%%ymm3 ,0x020(%%rdx)  	\n\t"/* __BCi = t10; */\
		/**/\
	"leaq	(%%r14,%%r14,2),%%r14	\n\t"	/* p4 + (p4*2) = p12, ptr-offset form */\
/* Swap Blocks 2,3 here to allow simple-incrementation of array ptrs (decr. not safe w.r.to array padding). */\
	/*...Block 2: t3/4,11/12,19/20,27/28 in ymm0-7: */\
		"addq	$0x40,%%rdi		\n\t"/* t3 */\
		"addq	%%r10,%%rax		\n\t"/* add0+p1 */\
		"addq	%%r10,%%rbx		\n\t"/* add0+p5 */\
		"addq	%%r10,%%rcx		\n\t"/* add0+p9 */\
		"addq	%%r10,%%rdx		\n\t"/* add0+pD */\
		"vmovaps		0x100(%%rdi),%%ymm2 	\n\t	vmovaps		0x120(%%rdi),%%ymm3 	\n\t"/*	t11,12; */\
		"vmovaps		%%ymm2,%%ymm8			\n\t	vmovaps		0x020(%%rsi),%%ymm13	\n\t"/* rt = t11; load __sc */\
		" vfmadd231pd	%%ymm10,%%ymm3 ,%%ymm2 	\n\t"/*	 FMA231(1.0,t12,t11); */\
		"vfnmadd231pd	%%ymm10,%%ymm8 ,%%ymm3 	\n\t"/*	FNMA231(1.0,rt ,t12); */\
		"vmovaps		0x200(%%rdi),%%ymm4 	\n\t	vmovaps		0x220(%%rdi),%%ymm5 	\n\t"/*	t19,20; */\
		"vmovaps		0x300(%%rdi),%%ymm6 	\n\t	vmovaps		0x320(%%rdi),%%ymm7 	\n\t"/*	t27,28; */\
		"vmovaps		%%ymm4,%%ymm11			\n\t	vmovaps		%%ymm6,%%ymm12	\n\t"/*	re = t19; im = t27; */\
		" vfmadd231pd	%%ymm13,%%ymm5 ,%%ymm4 	\n\t"/*	 FMA231(__sc,t20,t19); */\
		"vfnmadd231pd	%%ymm13,%%ymm11,%%ymm5 	\n\t"/*	FNMA231(__sc,re ,t20); */\
		" vfmadd132pd	%%ymm13,%%ymm7 ,%%ymm6 	\n\t"/*	 FMA132(__sc,t27,t28); */\
		" vfmsub132pd	%%ymm13,%%ymm12,%%ymm7 	\n\t"/*	 FMS132(__sc,t28,im ); */\
		"vmovaps		     (%%rdi),%%ymm0 	\n\t	vmovaps		0x020(%%rdi),%%ymm1 	\n\t"/*	t3,4 ; */\
		"vmovaps		-0x20(%%rsi),%%ymm14	\n\t	vmovaps		     (%%rsi),%%ymm15	\n\t"/* load ISRT2, __c */\
		"vmovaps		%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9	\n\t"/*	rt = t3; it = t4; */\
		" vfmadd231pd	%%ymm14,%%ymm2 ,%%ymm0 	\n\t"/*	 FMA231(ISRT2,t11,t3 ); */\
		"vfnmadd132pd	%%ymm14,%%ymm8 ,%%ymm2 	\n\t"/*	FNMA132(ISRT2,t11,rt ); */\
		" vfmadd231pd	%%ymm14,%%ymm3 ,%%ymm1 	\n\t"/*	 FMA231(ISRT2,t12,t4 ); */\
		"vfnmadd132pd	%%ymm14,%%ymm9 ,%%ymm3 	\n\t"/*	FNMA132(ISRT2,t12,it ); */\
		"vmulpd			%%ymm15,%%ymm6 ,%%ymm6 	\n\t"/*	t27 *= __c; */\
		"vmulpd			%%ymm15,%%ymm7 ,%%ymm7 	\n\t"/*	t28 *= __c; */\
		"vmovaps		%%ymm6,%%ymm11			\n\t	vmovaps		%%ymm7,%%ymm12	\n\t"/*	re = t27; im = t28; */\
		" vfmsub231pd	%%ymm15,%%ymm4 ,%%ymm6 	\n\t"/*	 FMS231(__c,t19,t27); */\
		" vfmadd132pd	%%ymm15,%%ymm11,%%ymm4 	\n\t"/*	 FMA132(__c,t19,re ); */\
		" vfmsub231pd	%%ymm15,%%ymm5 ,%%ymm7 	\n\t"/*	 FMS231(__c,t20,t28); */\
		" vfmadd132pd	%%ymm15,%%ymm12,%%ymm5 	\n\t"/*	 FMA132(__c,t20,im ); */\
		"vmovaps		%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9	\n\t"/*	rt = t3; it = t4; */\
		" vfmsub132pd	%%ymm10,%%ymm4 ,%%ymm0 	\n\t"/*	 FMS132(1.0,t3 ,t19); */\
		" vfmsub132pd	%%ymm10,%%ymm5 ,%%ymm1 	\n\t"/*	 FMS132(1.0,t4 ,t20); */\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm4 	\n\t"/*	 FMA231(1.0,rt ,t19); */\
		" vfmadd231pd	%%ymm10,%%ymm9 ,%%ymm5 	\n\t"/*	 FMA231(1.0,it ,t20); */\
		"vmovaps		0x260(%%rsi),%%ymm14	\n\t	vmovaps		0x2e0(%%rsi),%%ymm15	\n\t"/* load __t1,9 */\
		"vmovaps		%%ymm2,%%ymm11			\n\t	vmovaps		%%ymm3,%%ymm12	\n\t"/*	re = t11; im = t12; */\
		" vfmsub132pd	%%ymm10,%%ymm7 ,%%ymm2 	\n\t"/*	 FMS132(1.0,t11,t28); */\
		" vfmadd132pd	%%ymm10,%%ymm6 ,%%ymm3 	\n\t"/*	 FMA132(1.0,t12,t27); */\
		" vfmadd231pd	%%ymm10,%%ymm11,%%ymm7 	\n\t"/*	 FMA231(1.0,re ,t28); */\
		" vfmsub231pd	%%ymm10,%%ymm12,%%ymm6 	\n\t"/*	 FMS231(1.0,im ,t27); */\
		"vmovaps		%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm0,%%ymm9	\n\t"/*	rt = t19; it = t3; */\
		"vmovaps		0x2a0(%%rsi),%%ymm10	\n\t	vmovaps		0x320(%%rsi),%%ymm13	\n\t"/* load __t5,D; use 1.0 slot for __t5 */\
		" vfmadd231pd	%%ymm14,%%ymm5 ,%%ymm4 	\n\t"/*	 FMA231(__t1,t20,t19); */\
		"vfnmadd231pd	%%ymm14,%%ymm8 ,%%ymm5 	\n\t"/*	FNMA231(__t1,rt ,t20); */\
		" vfmadd231pd	%%ymm15,%%ymm1 ,%%ymm0 	\n\t"/*	 FMA231(__t9,t4 ,t3 ); */\
		"vfnmadd231pd	%%ymm15,%%ymm9 ,%%ymm1 	\n\t"/*	FNMA231(__t9,it ,t4 ); */\
		"vmovaps		0x240(%%rsi),%%ymm14	\n\t	vmovaps		0x2c0(%%rsi),%%ymm15	\n\t"/* load __c1,9 */\
		"vmovaps		%%ymm7,%%ymm11			\n\t	vmovaps		%%ymm2,%%ymm12	\n\t"/*	re = t28; im = t11; */\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm7 	\n\t"/*	 FMA231(__t5,t27,t28); */\
		"vfnmadd231pd	%%ymm10,%%ymm11,%%ymm6 	\n\t"/*	FNMA231(__t5,re ,t27); */\
		" vfmadd231pd	%%ymm13,%%ymm3 ,%%ymm2 	\n\t"/*	 FMA231(__tD,t12,t11); */\
		"vfnmadd231pd	%%ymm13,%%ymm12,%%ymm3 	\n\t"/*	FNMA231(__tD,im ,t12); */\
		"vmovaps		0x280(%%rsi),%%ymm10	\n\t	vmovaps		0x300(%%rsi),%%ymm13	\n\t"/* load __c5,D; use 1.0 slot for __t5 */\
		"vmulpd			%%ymm14,%%ymm4 ,%%ymm4 	\n\t"/*	t19 *= __c1; */\
		"vmulpd			%%ymm14,%%ymm5 ,%%ymm5 	\n\t"/*	t20 *= __c1; */\
		"vmulpd			%%ymm15,%%ymm0 ,%%ymm0 	\n\t"/*	t3  *= __c9; */\
		"vmulpd			%%ymm15,%%ymm1 ,%%ymm1 	\n\t"/*	t4  *= __c9; */\
		"vmulpd			%%ymm10,%%ymm7 ,%%ymm7 	\n\t"/*	t28 *= __c5; */\
		"vmulpd			%%ymm10,%%ymm6 ,%%ymm6 	\n\t"/*	t27 *= __c5; */\
		"vmulpd			%%ymm13,%%ymm2 ,%%ymm2 	\n\t"/*	t11 *= __cD; */\
		"vmulpd			%%ymm13,%%ymm3 ,%%ymm3 	\n\t"/*	t12 *= __cD; */\
		"vmovaps		%%ymm4 ,     (%%rax)  	\n\t"/* __B1r = t19; */\
		"vmovaps		%%ymm5 ,0x020(%%rax)  	\n\t"/* __B1i = t20; */\
		"vmovaps		%%ymm0 ,     (%%rcx)  	\n\t"/* __B9r = t3 ; */\
		"vmovaps		%%ymm1 ,0x020(%%rcx)  	\n\t"/* __B9i = t4 ;	Since cannot assume p2 = (p1+p1) due to array padding, revert:*/\
		"vmovaps		%%ymm7 ,     (%%rbx)  	\n\t	subq	%%r10,%%rax	\n\t"/* __B5r = t28; add0+p1 -> add0+p0 */\
		"vmovaps		%%ymm6 ,0x020(%%rbx)  	\n\t	subq	%%r10,%%rbx	\n\t"/* __B5i = t27; add0+p5 -> add0+p4 */\
		"vmovaps		%%ymm2 ,     (%%rdx)  	\n\t	subq	%%r10,%%rcx	\n\t"/* __BDr = t11; add0+p9 -> add0+p8 */\
		"vmovaps		%%ymm3 ,0x020(%%rdx)  	\n\t	subq	%%r10,%%rdx	\n\t"/* __BDi = t12; add0+pD -> add0+pC */\
		/**/\
	/*...Block 3: t5/6,13/14,21/22,29/30 in ymm0-7: */\
		"vmovaps 0x440(%%rsi),%%ymm10	\n\t"/* retore 1.0 to ymm10 */\
		"addq	$0x40,%%rdi		\n\t"/* t5 */\
		"movslq	%[__p2],%%r9	\n\t"\
		"shlq	$3,%%r9			\n\t"\
		"addq	%%r9,%%rax		\n\t"/* add0+p2 */\
		"addq	%%r9,%%rbx		\n\t"/* add0+p6 */\
		"addq	%%r9,%%rcx		\n\t"/* add0+pA */\
		"addq	%%r9,%%rdx		\n\t"/* add0+pE */\
		"vmovaps		     (%%rdi),%%ymm0 	\n\t	vmovaps		0x100(%%rdi),%%ymm2 	\n\t"/*	t5,13; */\
		"vmovaps		0x020(%%rdi),%%ymm1 	\n\t	vmovaps		0x120(%%rdi),%%ymm3 	\n\t"/*	t6,14; */\
		"vmovaps		%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9	\n\t"/*	rt = t5; it = t6; */\
		" vfmadd231pd	%%ymm10,%%ymm3 ,%%ymm0 	\n\t"/*	 FMA231(1.0,t14,t5 ); */\
		"vfnmadd132pd	%%ymm10,%%ymm8 ,%%ymm3 	\n\t"/*	FNMA132(1.0,t14,rt ); */\
		"vfnmadd231pd	%%ymm10,%%ymm2 ,%%ymm1 	\n\t"/*	FNMA231(1.0,t13,t6 ); */\
		" vfmadd132pd	%%ymm10,%%ymm9 ,%%ymm2 	\n\t"/*	 FMA132(1.0,t13,it ); */\
		"vmovaps		0x200(%%rdi),%%ymm4 	\n\t	vmovaps		0x220(%%rdi),%%ymm5 	\n\t"/*	t21,22; */\
		"vmovaps		0x300(%%rdi),%%ymm6 	\n\t	vmovaps		0x320(%%rdi),%%ymm7 	\n\t"/*	t29,30; */\
		"vmovaps		%%ymm5,%%ymm11			\n\t	vmovaps		%%ymm7,%%ymm12	\n\t"/*	re = t22; im = t30; */\
		"vfnmadd231pd	%%ymm10,%%ymm4 ,%%ymm5 	\n\t"/*	FNMA231(1.0,t21,t22); */\
		" vfmadd132pd	%%ymm10,%%ymm11,%%ymm4 	\n\t"/*	 FMA132(1.0,t21,re ); */\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm7 	\n\t"/*	 FMA231(1.0,t29,t30); */\
		" vfmsub132pd	%%ymm10,%%ymm12,%%ymm6 	\n\t"/*	 FMS132(1.0,t29,im ); */\
		"vmovaps		%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm5,%%ymm9	\n\t"/*	rt = t21; it = t22; */\
		"vfnmadd231pd	%%ymm10,%%ymm6 ,%%ymm4 	\n\t"/*	FNMA231(1.0,t29,t21); */\
		" vfmadd132pd	%%ymm10,%%ymm8 ,%%ymm6 	\n\t"/*	 FMA132(1.0,t29,rt ); */\
		"vfnmadd231pd	%%ymm10,%%ymm7 ,%%ymm5 	\n\t"/*	FNMA231(1.0,t30,t22); */\
		" vfmadd132pd	%%ymm10,%%ymm9 ,%%ymm7 	\n\t"/*	 FMA132(1.0,t30,it ); */\
		"vmovaps		-0x20(%%rsi),%%ymm13	\n\t"/* load ISRT2 */\
		"vmovaps		%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9	\n\t"/*	rt = t5; it = t6; */\
		"vfnmadd231pd	%%ymm13,%%ymm4 ,%%ymm0 	\n\t"/*	FNMA231(ISRT2,t21,t5 ); */\
		" vfmadd132pd	%%ymm13,%%ymm8 ,%%ymm4 	\n\t"/*	 FMA132(ISRT2,t21,rt ); */\
		"vfnmadd231pd	%%ymm13,%%ymm5 ,%%ymm1 	\n\t"/*	FNMA231(ISRT2,t22,t6 ); */\
		" vfmadd132pd	%%ymm13,%%ymm9 ,%%ymm5 	\n\t"/*	 FMA132(ISRT2,t22,it ); */\
		"vmovaps		0x160(%%rsi),%%ymm14	\n\t	vmovaps		0x1e0(%%rsi),%%ymm15	\n\t"/* load __t2,A */\
		"vmovaps		%%ymm2,%%ymm11			\n\t	vmovaps		%%ymm3,%%ymm12	\n\t"/*	re =t13; im =t14; */\
		" vfmadd231pd	%%ymm13,%%ymm6 ,%%ymm2 	\n\t"/*	 FMA231(ISRT2,t29,t13); */\
		"vfnmadd132pd	%%ymm13,%%ymm11,%%ymm6 	\n\t"/*	FNMA132(ISRT2,t29,re ); */\
		"vfnmadd231pd	%%ymm13,%%ymm7 ,%%ymm3 	\n\t"/*	FNMA231(ISRT2,t30,t14); */\
		" vfmadd132pd	%%ymm13,%%ymm12,%%ymm7 	\n\t"/*	 FMA132(ISRT2,t30,im ); */\
		"vmovaps		%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm0,%%ymm9	\n\t"/*	rt =t21; it = t5; */\
		"vmovaps		0x1a0(%%rsi),%%ymm10	\n\t	vmovaps		0x220(%%rsi),%%ymm13	\n\t"/* load __t6,E; use 1.0 slot for __t5 */\
		" vfmadd231pd	%%ymm14,%%ymm5 ,%%ymm4 	\n\t"/*	 FMA231(__t2,t22,t21); */\
		"vfnmadd231pd	%%ymm14,%%ymm8 ,%%ymm5 	\n\t"/*	FNMA231(__t2,rt ,t22); */\
		" vfmadd231pd	%%ymm15,%%ymm1 ,%%ymm0 	\n\t"/*	 FMA231(__tA,t6 ,t5 ); */\
		"vfnmadd231pd	%%ymm15,%%ymm9 ,%%ymm1 	\n\t"/*	FNMA231(__tA,it ,t6 ); */\
		"vmovaps		0x140(%%rsi),%%ymm14	\n\t	vmovaps		0x1c0(%%rsi),%%ymm15	\n\t"/* load __c2,A */\
		"vmovaps		%%ymm7,%%ymm11			\n\t	vmovaps		%%ymm3,%%ymm12	\n\t"/*	re =t30; im =t14; */\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm7 	\n\t"/*	 FMA231(__t6,t29,t30); */\
		"vfnmadd231pd	%%ymm10,%%ymm11,%%ymm6 	\n\t"/*	FNMA231(__t6,re ,t29); */\
		" vfmadd231pd	%%ymm13,%%ymm2 ,%%ymm3 	\n\t"/*	 FMA231(__tE,t13,t14); */\
		"vfnmadd231pd	%%ymm13,%%ymm12,%%ymm2 	\n\t"/*	FNMA231(__tE,im ,t13); */\
		"vmovaps		0x180(%%rsi),%%ymm10	\n\t	vmovaps		0x200(%%rsi),%%ymm13	\n\t"/* load __c6,E; use 1.0 slot for __t5 */\
		"vmulpd			%%ymm14,%%ymm4 ,%%ymm4 	\n\t"/*	t21 *= __c2; */\
		"vmulpd			%%ymm14,%%ymm5 ,%%ymm5 	\n\t"/*	t22 *= __c2; */\
		"vmulpd			%%ymm15,%%ymm0 ,%%ymm0 	\n\t"/*	t5  *= __cA; */\
		"vmulpd			%%ymm15,%%ymm1 ,%%ymm1 	\n\t"/*	t6  *= __cA; */\
		"vmulpd			%%ymm10,%%ymm7 ,%%ymm7 	\n\t"/*	t30 *= __c6; */\
		"vmulpd			%%ymm10,%%ymm6 ,%%ymm6 	\n\t"/*	t29 *= __c6; */\
		"vmulpd			%%ymm13,%%ymm3 ,%%ymm3 	\n\t"/*	t14 *= __cE; */\
		"vmulpd			%%ymm13,%%ymm2 ,%%ymm2 	\n\t"/*	t13 *= __cE; */\
		"vmovaps		%%ymm4 ,     (%%rax)  	\n\t"/* __B2r = t21; */\
		"vmovaps		%%ymm5 ,0x020(%%rax)  	\n\t"/* __B2i = t22; */\
		"vmovaps		%%ymm0 ,     (%%rcx)  	\n\t"/* __BAr = t5 ; */\
		"vmovaps		%%ymm1 ,0x020(%%rcx)  	\n\t"/* __BAi = t6 ; */\
		"vmovaps		%%ymm7 ,     (%%rbx)  	\n\t"/* __B6r = t30; */\
		"vmovaps		%%ymm6 ,0x020(%%rbx)  	\n\t"/* __B6i = t29; */\
		"vmovaps		%%ymm3 ,     (%%rdx)  	\n\t"/* __BEr = t14; */\
		"vmovaps		%%ymm2 ,0x020(%%rdx)  	\n\t"/* __BEi = t13; */\
		/**/\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%r14)\n\t"/* ...+p12 */\
	/*...Block 4: t7/8,15/16,23/24,31/32 in ymm0-7: */\
		"vmovaps 0x440(%%rsi),%%ymm10	\n\t"/* retore 1.0 to ymm10 */\
		"addq	$0x40,%%rdi		\n\t"/* t7 */\
		"addq	%%r10,%%rax		\n\t"/* add0+p3 */\
		"addq	%%r10,%%rbx		\n\t"/* add0+p7 */\
		"addq	%%r10,%%rcx		\n\t"/* add0+pB */\
		"addq	%%r10,%%rdx		\n\t"/* add0+pF */\
		"vmovaps		0x100(%%rdi),%%ymm2 	\n\t	vmovaps		0x120(%%rdi),%%ymm3 	\n\t"/*	t15,16; */\
		"vmovaps		%%ymm3,%%ymm8			\n\t	vmovaps		0x020(%%rsi),%%ymm13	\n\t"/* rt = t16; load __sc */\
		" vfmadd231pd	%%ymm10,%%ymm2 ,%%ymm3 	\n\t"/*	 FMA231(1.0,t15,t16); */\
		" vfmsub132pd	%%ymm10,%%ymm8 ,%%ymm2 	\n\t"/*	 FMS132(1.0,t15,rt ); */\
		"vmovaps		0x200(%%rdi),%%ymm4 	\n\t	vmovaps		0x220(%%rdi),%%ymm5 	\n\t"/*	t23,24; */\
		"vmovaps		0x300(%%rdi),%%ymm6 	\n\t	vmovaps		0x320(%%rdi),%%ymm7 	\n\t"/*	t31,32; */\
		"vmovaps		%%ymm4,%%ymm11			\n\t	vmovaps		%%ymm6,%%ymm12	\n\t"/*	re = t23; im = t31; */\
		" vfmadd231pd	%%ymm13,%%ymm7 ,%%ymm6 	\n\t"/*	 FMA231(__sc,t32,t31); */\
		"vfnmadd231pd	%%ymm13,%%ymm12,%%ymm7 	\n\t"/*	FNMA231(__sc,im ,t32); */\
		" vfmadd132pd	%%ymm13,%%ymm5 ,%%ymm4 	\n\t"/*	 FMA132(__sc,t23,t24); */\
		" vfmsub132pd	%%ymm13,%%ymm11,%%ymm5 	\n\t"/*	 FMS132(__sc,t24,re ); */\
		"vmovaps		-0x20(%%rsi),%%ymm14	\n\t	vmovaps		     (%%rsi),%%ymm15	\n\t"/* load ISRT2, __c */\
		"vmovaps		     (%%rdi),%%ymm0 	\n\t	vmovaps		0x020(%%rdi),%%ymm1 	\n\t"/*	t7,8 ; */\
		"vmovaps		%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9	\n\t"/*	rt = t7; it = t8; */\
		"vfnmadd231pd	%%ymm14,%%ymm2 ,%%ymm0 	\n\t"/*	FNMA231(ISRT2,t15,t7 ); */\
		" vfmadd132pd	%%ymm14,%%ymm8 ,%%ymm2 	\n\t"/*	 FMA132(ISRT2,t15,rt ); */\
		"vfnmadd231pd	%%ymm14,%%ymm3 ,%%ymm1 	\n\t"/*	FNMA231(ISRT2,t16,t8 ); */\
		" vfmadd132pd	%%ymm14,%%ymm9 ,%%ymm3 	\n\t"/*	 FMA132(ISRT2,t16,it ); */\
		"vmulpd			%%ymm15,%%ymm6 ,%%ymm6 	\n\t"/*	t31 *= __c; */\
		"vmulpd			%%ymm15,%%ymm7 ,%%ymm7 	\n\t"/*	t32 *= __c; */\
		"vmovaps		%%ymm6,%%ymm11			\n\t	vmovaps		%%ymm7,%%ymm12	\n\t"/*	re = t31; im = t32; */\
		" vfmadd231pd	%%ymm15,%%ymm4 ,%%ymm6 	\n\t"/*	 FMA231(__c,t23,t31); */\
		" vfmsub132pd	%%ymm15,%%ymm11,%%ymm4 	\n\t"/*	 FMS132(__c,t23,re ); */\
		" vfmadd231pd	%%ymm15,%%ymm5 ,%%ymm7 	\n\t"/*	 FMA231(__c,t24,t32); */\
		" vfmsub132pd	%%ymm15,%%ymm12,%%ymm5 	\n\t"/*	 FMS132(__c,t24,im ); */\
		"vmovaps		%%ymm2,%%ymm8			\n\t	vmovaps		%%ymm3,%%ymm9	\n\t"/*	rt = t15; it = t16; */\
		" vfmsub132pd	%%ymm10,%%ymm7 ,%%ymm2 	\n\t"/*	 FMS132(1.0,t15,t32); */\
		" vfmadd132pd	%%ymm10,%%ymm6 ,%%ymm3 	\n\t"/*	 FMA132(1.0,t16,t31); */\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm7 	\n\t"/*	 FMA231(1.0,rt ,t32); */\
		" vfmsub231pd	%%ymm10,%%ymm9 ,%%ymm6 	\n\t"/*	 FMS231(1.0,it ,t31); */\
		"vmovaps		0x360(%%rsi),%%ymm14	\n\t	vmovaps		0x3e0(%%rsi),%%ymm15	\n\t"/* load __t3,B */\
		"vmovaps		%%ymm0,%%ymm11			\n\t	vmovaps		%%ymm1,%%ymm12	\n\t"/*	re = t7; im = t8; */\
		" vfmsub132pd	%%ymm10,%%ymm4 ,%%ymm0 	\n\t"/*	 FMS132(1.0,t7 ,t23 ); */\
		" vfmsub132pd	%%ymm10,%%ymm5 ,%%ymm1 	\n\t"/*	 FMS132(1.0,t8 ,t24 ); */\
		" vfmadd231pd	%%ymm10,%%ymm11,%%ymm4 	\n\t"/*	 FMA231(1.0,re ,t23 ); */\
		" vfmadd231pd	%%ymm10,%%ymm12,%%ymm5 	\n\t"/*	 FMA231(1.0,im ,t24 ); */\
		"vmovaps		%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm0,%%ymm9	\n\t"/*	rt = t23; it = t7; */\
		"vmovaps		0x3a0(%%rsi),%%ymm10	\n\t	vmovaps		0x420(%%rsi),%%ymm13	\n\t"/* load __t7,F; use 1.0 slot for __t7 */\
		" vfmadd231pd	%%ymm14,%%ymm5 ,%%ymm4 	\n\t"/*	 FMA231(__t3,t24,t23); */\
		"vfnmadd231pd	%%ymm14,%%ymm8 ,%%ymm5 	\n\t"/*	FNMA231(__t3,rt ,t24); */\
		" vfmadd231pd	%%ymm15,%%ymm1 ,%%ymm0 	\n\t"/*	 FMA231(__tB,t8 ,t7 ); */\
		"vfnmadd231pd	%%ymm15,%%ymm9 ,%%ymm1 	\n\t"/*	FNMA231(__tB,it ,t8 ); */\
		"vmovaps		0x340(%%rsi),%%ymm14	\n\t	vmovaps		0x3c0(%%rsi),%%ymm15	\n\t"/* load __c3,B */\
		"vmovaps		%%ymm7,%%ymm11			\n\t	vmovaps		%%ymm2,%%ymm12	\n\t"/*	re = t32; im = t15; */\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm7 	\n\t"/*	 FMA231(__t7,t31,t32); */\
		"vfnmadd231pd	%%ymm10,%%ymm11,%%ymm6 	\n\t"/*	FNMA231(__t7,re ,t31); */\
		" vfmadd231pd	%%ymm13,%%ymm3 ,%%ymm2 	\n\t"/*	 FMA231(__tF,t16,t15); */\
		"vfnmadd231pd	%%ymm13,%%ymm12,%%ymm3 	\n\t"/*	FNMA231(__tF,im ,t16); */\
		"vmovaps		0x380(%%rsi),%%ymm10	\n\t	vmovaps		0x400(%%rsi),%%ymm13	\n\t"/* load __c7,F; use 1.0 slot for __t5 */\
		"vmulpd			%%ymm14,%%ymm4 ,%%ymm4 	\n\t"/*	t23 *= __c3; */\
		"vmulpd			%%ymm14,%%ymm5 ,%%ymm5 	\n\t"/*	t24 *= __c3; */\
		"vmulpd			%%ymm15,%%ymm0 ,%%ymm0 	\n\t"/*	t7  *= __cB; */\
		"vmulpd			%%ymm15,%%ymm1 ,%%ymm1 	\n\t"/*	t8  *= __cB; */\
		"vmulpd			%%ymm10,%%ymm7 ,%%ymm7 	\n\t"/*	t32 *= __c7; */\
		"vmulpd			%%ymm10,%%ymm6 ,%%ymm6 	\n\t"/*	t31 *= __c7; */\
		"vmulpd			%%ymm13,%%ymm2 ,%%ymm2 	\n\t"/*	t15 *= __cF; */\
		"vmulpd			%%ymm13,%%ymm3 ,%%ymm3 	\n\t"/*	t16 *= __cF; */\
		"vmovaps		%%ymm4 ,     (%%rax)  	\n\t"/* __B3r = t23; */\
		"vmovaps		%%ymm5 ,0x020(%%rax)  	\n\t"/* __B3i = t24; */\
		"vmovaps		%%ymm0 ,     (%%rcx)  	\n\t"/* __BBr = t7 ; */\
		"vmovaps		%%ymm1 ,0x020(%%rcx)  	\n\t"/* __BBi = t8 ; */\
		"vmovaps		%%ymm7 ,     (%%rbx)  	\n\t"/* __B7r = t32; */\
		"vmovaps		%%ymm6 ,0x020(%%rbx)  	\n\t"/* __B7i = t31; */\
		"vmovaps		%%ymm2 ,     (%%rdx)  	\n\t"/* __BFr = t15; */\
		"vmovaps		%%ymm3 ,0x020(%%rdx)  	\n\t"/* __BFi = t16; */\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r1] "m" (Xr1)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r9","r10","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	// 2nd version of trig-data and DFT-pass routines uses 1-copy trig data, read in and broadcast to full YMM register on-the-fly.
	// Compare the simplicity - at the cost of ~10-cycle latency of vdivpd - of this routine to its V1 vector-data counterpart,
	// which uses Newtonian-iterative inversion:
	#define RADIX16_COMPUTE_FMA_SINCOS_DIT_2(Xadd0)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax			\n\t"/* base-ptr to 32-byte-aligned vector-cosine data holding 16 double-precision divisors */\
		"vmovaps	0x80(%%rax),%%ymm0	\n\t"/* load sines [dividends] into 4 AVX registers */\
		"vmovaps	0xa0(%%rax),%%ymm1	\n\t"\
		"vmovaps	0xc0(%%rax),%%ymm2	\n\t"\
		"vmovaps	0xe0(%%rax),%%ymm3	\n\t"\
		/* Divide by cosine terms to get quotients in ymm0-3: */\
		"vdivpd	    (%%rax),%%ymm0,%%ymm0	\n\t"\
		"vdivpd	0x20(%%rax),%%ymm1,%%ymm1	\n\t"\
		"vdivpd	0x40(%%rax),%%ymm2,%%ymm2	\n\t"\
		"vdivpd	0x60(%%rax),%%ymm3,%%ymm3	\n\t"\
		/* Outputs overwrite sines: */\
		"vmovaps	%%ymm0,0x80(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0xa0(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0xc0(%%rax)	\n\t"\
		"vmovaps	%%ymm3,0xe0(%%rax)	\n\t"\
		/* Scalar data starting at add0 = cc0 laid out as below. Ensuing C code assumed to replace c0 by cc0 and r0 by tan0 = cc0/ss0. */\
		/* add0 + 0x[  0,  8, 10, 18]: [c0 ,c1 ,c2 ,c3 ] Cosines: */\
		/* add0 + 0x[ 20, 28, 30, 38]: [c4 ,c5 ,c6 ,c7 ] */\
		/* add0 + 0x[ 40, 48, 50, 58]: [c8 ,c9 ,cA ,cB ] */\
		/* add0 + 0x[ 60, 68, 70, 78]: [cC ,cD ,cE ,cF ] */\
		/* add0 + 0x[ 80, 88, 90, 98]: [r0 ,r1 ,r2 ,r3 ] Tangents: */\
		/* add0 + 0x[ a0, a8, b0, b8]: [r4 ,r5 ,r6 ,r7 ] */\
		/* add0 + 0x[ c0, c8, d0, d8]: [r8 ,r9 ,rA ,rB ] */\
		/* add0 + 0x[ e0, e8, f0, f8]: [rC ,rD ,rE ,rF ] */\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3"	/* Clobbered registers */\
	);\
	}

	// Array-of-doubles index offsets w.r.to the __c = cc0 base-root address of the various derived sincos terms:
    // Datum  Offset          Datum  Offset
    //        idx   byte             idx   byte
    // ----   ----------      ----   ----------
    // ISRT2 -0x01 -0x08
    // __c    0x00  0x00      __sc   0x10  0x80
    // __c1   0x01  0x08      __t1   0x11  0x88
    // __c2   0x02  0x10      __t2   0x12  0x90
    // __c3   0x03  0x18      __t3   0x13  0x98
    // __c4   0x04  0x20      __t4   0x14  0xa0
    // __c5   0x05  0x28      __t5   0x15  0xa8
    // __c6   0x06  0x30      __t6   0x16  0xb0
    // __c7   0x07  0x38      __t7   0x17  0xb8
    // __c8   0x08  0x40      __t8   0x18  0xc0
    // __c9   0x09  0x48      __t9   0x19  0xc8
    // __cA   0x0a  0x50      __tA   0x1a  0xd0
    // __cB   0x0b  0x58      __tB   0x1b  0xd8
    // __cC   0x0c  0x60      __tC   0x1c  0xe0
    // __cD   0x0d  0x68      __tD   0x1d  0xe8
    // __cE   0x0e  0x70      __tE   0x1e  0xf0
    // __cF   0x0f  0x78      __tF   0x1f  0xf8
	//                        1.0    0x20  0x100

	// Version 2 of the DIT macro:
	//
	// [1] Implements 'Option B' of RADIX_16_DIT_FMA, which delays the *__c multiply of
	//     the t23/24,t31/32 outputs of the initial radix-2 butterflies in the final radix-4 subtransform
	//     in an attempt to make for better scheduling;
	// [2] Uses scalar-double/vbroadcastsd scheme for trig data, inited by above RADIX16_COMPUTE_FMA_SINCOS_DIT_2 macro.
	//
	/* Sep 2014: Since my DIT impl is post-twiddles, that is twiddles-applied-to-outputs, original all-FMA version of DIT
	is dominated by "trivial-MUL" FMAs in which one of the 2 multiplicands is unity. Haswell can in theory issue FMAs at
	2x the rate of ADD/SUB (i.e. 2 per cycle vs 1) so this should be better in terms of throughput, but due to longer
	latency of FMA (5 cycles vs 3) it's worth also trying a variant in which we replace the trivial-MUL FMAs with their
	equivalent ADD/SUB. Tried such a modified version, replacing 110 trivial-MUL FMAs with equivalent ADD/SUB - slower.
	*/
	// 0 ADD, 174 FMA, 34 pure-MUL, 219 MOVAPS
	#define SSE2_RADIX16_DIT_TWIDDLE_2(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr1,Xcc0,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/* Gather the needed data and do first set of four length-4 transforms. Each complex radix-4 block needs 11 registers: */\
	/* 8 t-temps mapped to ymm0-7, plus rt,it (ymm8,9), plus const 1.0 (ymm10), though we add 2 more tmps re,im (ymm11,12) */\
	/* to allow independent Re/Im subsections to overlap w/o introducing false dependencies. */\
	/*...Block 1: */\
		"movslq	%[__p4],%%r10	\n\t"/* This serves as main-array stride-between radix-4 blocks */\
		"movq	%[__add0],%%rax				\n\t"\
	"movslq	%[__pfetch_addr1],%%r13	\n\t"	/* Prefetch base-index offset cycles among add0 + p0,1,2,3 on successive macro calls */\
	"leaq	(%%rax,%%r13,8),%%r13	\n\t"	/* [base-address + data-fetch-ahead index] */\
		"movslq	%[__p1],%%rbx	\n\t"\
		"movslq	%[__p2],%%rcx	\n\t"\
		"movslq	%[__p3],%%rdx	\n\t"\
		"movq	%[__cc0],%%rsi 	\n\t"\
		"shlq	$3,%%r10		\n\t"/* p4 in bytewise ptr-arithmetic form */\
		"vmovaps 0x440(%%rsi),%%ymm10	\n\t"/* cc0 + 0x22 = __two; ACTUALLY HOLDS 1.0 IN AVX2 MODE */\
		"leaq	(%%rax,%%rbx,8),%%rbx	\n\t"/* add0+p1 */\
		"leaq	(%%rax,%%rcx,8),%%rcx	\n\t"/* add0+p2 */\
		"leaq	(%%rax,%%rdx,8),%%rdx	\n\t"/* add0+p3 */\
	"movq	%%r10,%%r14	\n\t"	/* Save a copy of p4 pointer offset. Will prefetch from [base-address + data-fetch-ahead index] + [0,p4,p8,p12] on each macro call. */\
	"prefetcht1	%c[__pfetch_dist](%%r13)\n\t"/* [base-address + data-fetch-ahead index] + p0 */\
		"movq	%[__r1],%%rdi	\n\t"/* ptr to local-mem [t1] */\
		"vmovaps		     (%%rax),%%ymm0 	\n\t	vmovaps		     (%%rbx),%%ymm8 	\n\t"/*	t1,rt =__A0,1r; */\
		"vmovaps		0x020(%%rax),%%ymm1 	\n\t	vmovaps		0x020(%%rbx),%%ymm11	\n\t"/*	t2,re =__A0,1i; */\
		"vmovaps		     (%%rcx),%%ymm4 	\n\t	vmovaps		     (%%rdx),%%ymm9 	\n\t"/*	t5,it =__A2,3r; */\
		"vmovaps		0x020(%%rcx),%%ymm5 	\n\t	vmovaps		0x020(%%rdx),%%ymm12	\n\t"/*	t6,im =__A2,3i; */\
		"vmovaps		%%ymm0,%%ymm2			\n\t	vmovaps		%%ymm4,%%ymm7	\n\t"/*	t3 = t1; t8 = t5; */\
		"vmovaps		%%ymm1,%%ymm3			\n\t	vmovaps		%%ymm5,%%ymm6	\n\t"/*	t4 = t2; t7 = t6; */\
		"vfnmadd231pd	%%ymm10,%%ymm8 ,%%ymm2 	\n\t	vfnmadd231pd	%%ymm10,%%ymm11,%%ymm3 	\n\t"/*	FNMA231(1.0,rt ,t3 ); FNMA231(1.0,re ,t4 ); */\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm0 	\n\t	 vfmadd231pd	%%ymm10,%%ymm11,%%ymm1 	\n\t"/*	 FMA231(1.0,rt ,t1 );  FMA231(1.0,re ,t2 ); */\
		"vfnmadd231pd	%%ymm10,%%ymm9 ,%%ymm7 	\n\t	vfnmadd231pd	%%ymm10,%%ymm12,%%ymm6 	\n\t"/*	FNMA231(1.0,it ,t8 ); FNMA231(1.0,im ,t7 ); */\
		" vfmadd231pd	%%ymm10,%%ymm9 ,%%ymm4 	\n\t	 vfmadd231pd	%%ymm10,%%ymm12,%%ymm5 	\n\t"/*	 FMA231(1.0,it ,t5 );  FMA231(1.0,im ,t6 ); */\
		"vmovaps		%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm2,%%ymm9	\n\t"/*	rt = t5; it = t3; */\
		/*** REMEMBER - For FMA132, rightmost 2 operands in ASM are reverse-ordered w.r.to the prototyping code in the comments! ***/\
		"vmovaps		%%ymm5,%%ymm11			\n\t	vmovaps		%%ymm3,%%ymm12		\n\t"/*	re = t6; im = t4; */\
		"vfnmadd132pd	%%ymm10,%%ymm0 ,%%ymm4 	\n\t	vfnmadd132pd	%%ymm10,%%ymm1 ,%%ymm5 	\n\t"/*	FNMA132(1.0,t5 ,t1 ); FNMA132(1.0,t6 ,t2 ); */\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm0 	\n\t	 vfmadd231pd	%%ymm10,%%ymm11,%%ymm1 	\n\t"/*	 FMA231(1.0,rt ,t1 );  FMA231(1.0,re ,t2 ); */\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm2 	\n\t	vfnmadd231pd	%%ymm10,%%ymm7 ,%%ymm3 	\n\t"/*	 FMA231(1.0,t7 ,t3 ); FNMA231(1.0,t8 ,t4 ); */\
		"vfnmadd132pd	%%ymm10,%%ymm9 ,%%ymm6 	\n\t	 vfmadd132pd	%%ymm10,%%ymm12,%%ymm7 	\n\t"/*	FNMA132(1.0,t7 ,it );  FMA132(1.0,t8 ,im ); */\
		"vmovaps %%ymm4 ,0x080(%%rdi) \n\t vmovaps %%ymm5 ,0x0a0(%%rdi)   \n\t addq %%r10,%%rax \n\t"/* write t5; write t6; add0+p4 */\
		"vmovaps %%ymm0 ,     (%%rdi) \n\t vmovaps %%ymm1 ,0x020(%%rdi)   \n\t addq %%r10,%%rbx \n\t"/* write t1; write t2; add0+p5 */\
		"vmovaps %%ymm2 ,0x040(%%rdi) \n\t vmovaps %%ymm3 ,0x060(%%rdi)   \n\t addq %%r10,%%rcx \n\t"/* write t3; write t4; add0+p6 */\
		"vmovaps %%ymm6 ,0x0c0(%%rdi) \n\t vmovaps %%ymm7 ,0x0e0(%%rdi)   \n\t addq %%r10,%%rdx \n\t"/* write t7; write t8; add0+p7 */\
		/**/\
	/*...Block 2: __A4-7, t9-16 */\
		"addq	$0x100,%%rdi	\n\t"/* ptr to t9 */\
		"vmovaps		     (%%rax),%%ymm0 	\n\t	vmovaps		     (%%rbx),%%ymm8 	\n\t"\
		"vmovaps		0x020(%%rax),%%ymm1 	\n\t	vmovaps		0x020(%%rbx),%%ymm11	\n\t"\
		"vmovaps		     (%%rcx),%%ymm4 	\n\t	vmovaps		     (%%rdx),%%ymm9 	\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm5 	\n\t	vmovaps		0x020(%%rdx),%%ymm12	\n\t"\
		"vmovaps		%%ymm0,%%ymm2			\n\t	vmovaps		%%ymm4,%%ymm7	\n\t"\
		"vmovaps		%%ymm1,%%ymm3			\n\t	vmovaps		%%ymm5,%%ymm6	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm8 ,%%ymm2 	\n\t	vfnmadd231pd	%%ymm10,%%ymm11,%%ymm3 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm0 	\n\t	 vfmadd231pd	%%ymm10,%%ymm11,%%ymm1 	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm9 ,%%ymm7 	\n\t	vfnmadd231pd	%%ymm10,%%ymm12,%%ymm6 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm9 ,%%ymm4 	\n\t	 vfmadd231pd	%%ymm10,%%ymm12,%%ymm5 	\n\t"\
		"vmovaps		%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm2,%%ymm9	\n\t"\
		"vmovaps		%%ymm5,%%ymm11			\n\t	vmovaps		%%ymm3,%%ymm12	\n\t"\
		"vfnmadd132pd	%%ymm10,%%ymm0 ,%%ymm4 	\n\t	vfnmadd132pd	%%ymm10,%%ymm1 ,%%ymm5 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm0 	\n\t	 vfmadd231pd	%%ymm10,%%ymm11,%%ymm1 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm2 	\n\t	vfnmadd231pd	%%ymm10,%%ymm7 ,%%ymm3 	\n\t"\
		"vfnmadd132pd	%%ymm10,%%ymm9 ,%%ymm6 	\n\t	 vfmadd132pd	%%ymm10,%%ymm12,%%ymm7 	\n\t"\
		"vmovaps %%ymm4 ,0x080(%%rdi) \n\t vmovaps %%ymm5 ,0x0a0(%%rdi)   \n\t addq %%r10,%%rax \n\t"/*add0+p8 */\
		"vmovaps %%ymm0 ,     (%%rdi) \n\t vmovaps %%ymm1 ,0x020(%%rdi)   \n\t addq %%r10,%%rbx \n\t"/*add0+p9 */\
		"vmovaps %%ymm2 ,0x040(%%rdi) \n\t vmovaps %%ymm3 ,0x060(%%rdi)   \n\t addq %%r10,%%rcx \n\t"/*add0+pA */\
		"vmovaps %%ymm6 ,0x0c0(%%rdi) \n\t vmovaps %%ymm7 ,0x0e0(%%rdi)   \n\t addq %%r10,%%rdx \n\t"/*add0+pB */\
		/**/\
	/*...Block 3: __A8-B, t17-24 */\
		"addq	$0x100,%%rdi	\n\t"/* ptr to t17 */\
		"vmovaps		     (%%rax),%%ymm0 	\n\t	vmovaps		     (%%rbx),%%ymm8 	\n\t"\
		"vmovaps		0x020(%%rax),%%ymm1 	\n\t	vmovaps		0x020(%%rbx),%%ymm11	\n\t"\
		"vmovaps		     (%%rcx),%%ymm4 	\n\t	vmovaps		     (%%rdx),%%ymm9 	\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm5 	\n\t	vmovaps		0x020(%%rdx),%%ymm12	\n\t"\
		"vmovaps		%%ymm0,%%ymm2			\n\t	vmovaps		%%ymm4,%%ymm7	\n\t"\
		"vmovaps		%%ymm1,%%ymm3			\n\t	vmovaps		%%ymm5,%%ymm6	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm8 ,%%ymm2 	\n\t	vfnmadd231pd	%%ymm10,%%ymm11,%%ymm3 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm0 	\n\t	 vfmadd231pd	%%ymm10,%%ymm11,%%ymm1 	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm9 ,%%ymm7 	\n\t	vfnmadd231pd	%%ymm10,%%ymm12,%%ymm6 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm9 ,%%ymm4 	\n\t	 vfmadd231pd	%%ymm10,%%ymm12,%%ymm5 	\n\t"\
		"vmovaps		%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm2,%%ymm9	\n\t"\
		"vmovaps		%%ymm5,%%ymm11			\n\t	vmovaps		%%ymm3,%%ymm12	\n\t"\
		"vfnmadd132pd	%%ymm10,%%ymm0 ,%%ymm4 	\n\t	vfnmadd132pd	%%ymm10,%%ymm1 ,%%ymm5 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm0 	\n\t	 vfmadd231pd	%%ymm10,%%ymm11,%%ymm1 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm2 	\n\t	vfnmadd231pd	%%ymm10,%%ymm7 ,%%ymm3 	\n\t"\
		"vfnmadd132pd	%%ymm10,%%ymm9 ,%%ymm6 	\n\t	 vfmadd132pd	%%ymm10,%%ymm12,%%ymm7 	\n\t"\
		"vmovaps %%ymm4 ,0x080(%%rdi) \n\t vmovaps %%ymm5 ,0x0a0(%%rdi)   \n\t addq %%r10,%%rax \n\t"/*add0+pC */\
		"vmovaps %%ymm0 ,     (%%rdi) \n\t vmovaps %%ymm1 ,0x020(%%rdi)   \n\t addq %%r10,%%rbx \n\t"/*add0+pD */\
		"vmovaps %%ymm2 ,0x040(%%rdi) \n\t vmovaps %%ymm3 ,0x060(%%rdi)   \n\t addq %%r10,%%rcx \n\t"/*add0+pE */\
		"vmovaps %%ymm6 ,0x0c0(%%rdi) \n\t vmovaps %%ymm7 ,0x0e0(%%rdi)   \n\t addq %%r10,%%rdx \n\t"/*add0+pF */\
		/**/\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%r14)\n\t"/* ...+p4 */\
	/*...Block 4: __AC-F, t25-32 */\
		"addq	$0x100,%%rdi	\n\t"/* ptr to t25 */\
		"vmovaps		     (%%rax),%%ymm0 	\n\t	vmovaps		     (%%rbx),%%ymm8 	\n\t"\
		"vmovaps		0x020(%%rax),%%ymm1 	\n\t	vmovaps		0x020(%%rbx),%%ymm11	\n\t"\
		"vmovaps		     (%%rcx),%%ymm4 	\n\t	vmovaps		     (%%rdx),%%ymm9 	\n\t"\
		"vmovaps		0x020(%%rcx),%%ymm5 	\n\t	vmovaps		0x020(%%rdx),%%ymm12	\n\t"\
		"vmovaps		%%ymm0,%%ymm2			\n\t	vmovaps		%%ymm4,%%ymm7	\n\t"\
		"vmovaps		%%ymm1,%%ymm3			\n\t	vmovaps		%%ymm5,%%ymm6	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm8 ,%%ymm2 	\n\t	vfnmadd231pd	%%ymm10,%%ymm11,%%ymm3 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm0 	\n\t	 vfmadd231pd	%%ymm10,%%ymm11,%%ymm1 	\n\t"\
		"vfnmadd231pd	%%ymm10,%%ymm9 ,%%ymm7 	\n\t	vfnmadd231pd	%%ymm10,%%ymm12,%%ymm6 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm9 ,%%ymm4 	\n\t	 vfmadd231pd	%%ymm10,%%ymm12,%%ymm5 	\n\t"\
		"vmovaps		%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm2,%%ymm9	\n\t"\
		"vmovaps		%%ymm5,%%ymm11			\n\t	vmovaps		%%ymm3,%%ymm12	\n\t"\
		"vfnmadd132pd	%%ymm10,%%ymm0 ,%%ymm4 	\n\t	vfnmadd132pd	%%ymm10,%%ymm1 ,%%ymm5 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm0 	\n\t	 vfmadd231pd	%%ymm10,%%ymm11,%%ymm1 	\n\t"\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm2 	\n\t	vfnmadd231pd	%%ymm10,%%ymm7 ,%%ymm3 	\n\t"\
		"vfnmadd132pd	%%ymm10,%%ymm9 ,%%ymm6 	\n\t	 vfmadd132pd	%%ymm10,%%ymm12,%%ymm7 	\n\t"\
		"vmovaps %%ymm4 ,0x080(%%rdi) \n\t vmovaps %%ymm5 ,0x0a0(%%rdi)   \n\t"\
		"vmovaps %%ymm0 ,     (%%rdi) \n\t vmovaps %%ymm1 ,0x020(%%rdi)   \n\t"\
		"vmovaps %%ymm2 ,0x040(%%rdi) \n\t vmovaps %%ymm3 ,0x060(%%rdi)   \n\t"\
		"vmovaps %%ymm6 ,0x0c0(%%rdi) \n\t vmovaps %%ymm7 ,0x0e0(%%rdi)   \n\t"\
		"subq	$0x300,%%rdi	\n\t"/* Revert rdi to point at t1 */\
		/**/\
/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\
		"movq	%%r10,%%rbx	/* p4 already in ptr-offset form in r10 */\n\t"\
		"movslq	%[__p1],%%r10	\n\t"/* This serves as main-array stride-between radix-4 blocks */\
		"shlq	$3,%%r10		\n\t"\
		/**/\
	/*...Block 1: t1/2,9/10,17/18,25/26 in ymm0-7, resp. */\
		"movq	%[__add0],%%rax			\n\t"\
		"movslq	%[__p12],%%rdx			\n\t"\
		"leaq	(%%rax,%%rbx,2),%%rcx	\n\t"/* add0+p8 */\
		"addq	%%rax,%%rbx				\n\t"/* add0+p4 */\
		"leaq	(%%rax,%%rdx,8),%%rdx	\n\t"/* add0+pC */\
		"vmovaps		     (%%rdi),%%ymm0 	\n\t	vmovaps		0x100(%%rdi),%%ymm2 	\n\t"/*	t1,9 ; */\
		"vmovaps		0x020(%%rdi),%%ymm1 	\n\t	vmovaps		0x120(%%rdi),%%ymm3 	\n\t"/*	t2,10; */\
		"vmovaps		%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9	\n\t"/*	rt = t1; it = t2; */\
		" vfmadd231pd	%%ymm10,%%ymm2 ,%%ymm0 	\n\t"/*	 FMA231(1.0,t9 ,t1 ); */\
		"vfnmadd132pd	%%ymm10,%%ymm8 ,%%ymm2 	\n\t"/*	FNMA132(1.0,t9 ,rt ); */\
		" vfmadd231pd	%%ymm10,%%ymm3 ,%%ymm1 	\n\t"/*	 FMA231(1.0,t10,t2 ); */\
		"vfnmadd132pd	%%ymm10,%%ymm9 ,%%ymm3 	\n\t"/*	FNMA132(1.0,t10,it ); */\
		"vmovaps		0x200(%%rdi),%%ymm4 	\n\t	vmovaps		0x300(%%rdi),%%ymm6 	\n\t"/*	t17,25; */\
		"vmovaps		0x220(%%rdi),%%ymm5 	\n\t	vmovaps		0x320(%%rdi),%%ymm7 	\n\t"/*	t18,26; */\
		"vmovaps		%%ymm4,%%ymm11			\n\t	vmovaps		%%ymm5,%%ymm12	\n\t"/*	re = t17; im = t18; */\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm4 	\n\t"/*	 FMA231(1.0,t25,t17); */\
		"vfnmadd132pd	%%ymm10,%%ymm11,%%ymm6 	\n\t"/*	FNMA132(1.0,t25,re ); */\
		" vfmadd231pd	%%ymm10,%%ymm7 ,%%ymm5 	\n\t"/*	 FMA231(1.0,t26,t18); */\
		"vfnmadd132pd	%%ymm10,%%ymm12,%%ymm7 	\n\t"/*	FNMA132(1.0,t26,im ); */\
		"vmovaps		%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9	\n\t"/*	rt = t1; it = t2; */\
		" vfmadd231pd	%%ymm10,%%ymm4 ,%%ymm8 	\n\t"/*	 FMA231(1.0,t17,rt ); */\
		" vfmadd231pd	%%ymm10,%%ymm5 ,%%ymm9 	\n\t"/*	 FMA231(1.0,t18,it ); */\
		"vfnmadd231pd	%%ymm10,%%ymm4 ,%%ymm0 	\n\t"/*	FNMA231(1.0,t17,t1 ); */\
		"vfnmadd231pd	%%ymm10,%%ymm5 ,%%ymm1 	\n\t"/*	FNMA231(1.0,t18,t2 ); */\
		"vmovaps		%%ymm2,%%ymm11			\n\t	vmovaps		%%ymm3,%%ymm12	\n\t"/*	re = t9 ; im = t10; */\
		"vbroadcastsd 0xa0(%%rsi),%%ymm13 \n\t vbroadcastsd 0xc0(%%rsi),%%ymm14 \n\t vbroadcastsd 0xe0(%%rsi),%%ymm15 \n\t"/* __t4,8,C; */\
		"vfnmadd231pd	%%ymm10,%%ymm7 ,%%ymm2 	\n\t"/*	FNMA231(1.0,t26,t9 ); */\
		" vfmadd132pd	%%ymm10,%%ymm11,%%ymm7 	\n\t"/*	 FMA132(1.0,t26,re ); */\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm3 	\n\t"/*	 FMA231(1.0,t25,t10); */\
		"vfnmadd132pd	%%ymm10,%%ymm12,%%ymm6 	\n\t"/*	FNMA132(1.0,t25,im ); */\
		"vmovaps		%%ymm8 ,     (%%rax)  	\n\t"/* __B0r = rt; */\
		"vmovaps		%%ymm9 ,0x020(%%rax)  	\n\t"/* __B0i = it; */\
		"vmovaps %%ymm0,%%ymm8	\n\t	vmovaps %%ymm7,%%ymm11 \n\t vmovaps %%ymm2,%%ymm12	\n\t"/*	rt =t1; re =t26; im =t9; */\
		"vbroadcastsd 0x20(%%rsi),%%ymm4  \n\t vbroadcastsd 0x040(%%rsi),%%ymm5  \n\t vbroadcastsd 0x60(%%rsi),%%ymm9  \n\t"/* __c4,8,C; use unused ymms for these */\
		" vfmadd231pd	%%ymm14,%%ymm1 ,%%ymm0 	\n\t"/*	 FMA231(__t8,t2 ,t1 ); */\
		"vfnmadd231pd	%%ymm14,%%ymm8 ,%%ymm1 	\n\t"/*	FNMA231(__t8,rt ,t2 ); */\
		" vfmadd231pd	%%ymm13,%%ymm6 ,%%ymm7 	\n\t"/*	 FMA231(__t4,t25,t26); */\
		"vfnmadd231pd	%%ymm13,%%ymm11,%%ymm6 	\n\t"/*	FNMA231(__t4,re ,t25); */\
		" vfmadd231pd	%%ymm15,%%ymm3 ,%%ymm2 	\n\t"/*	 FMA231(__tC,t10,t9 ); */\
		"vfnmadd231pd	%%ymm15,%%ymm12,%%ymm3 	\n\t"/*	FNMA231(__tC,im ,t10); */\
		"vmulpd			%%ymm5 ,%%ymm0 ,%%ymm0 	\n\t"/*	t1  *= __c8; */\
		"vmulpd			%%ymm5 ,%%ymm1 ,%%ymm1 	\n\t"/*	t2  *= __c8; */\
		"vmulpd			%%ymm4 ,%%ymm7 ,%%ymm7 	\n\t"/*	t26 *= __c4; */\
		"vmulpd			%%ymm4 ,%%ymm6 ,%%ymm6 	\n\t"/*	t25 *= __c4; */\
		"vmulpd			%%ymm9 ,%%ymm2 ,%%ymm2 	\n\t"/*	t9  *= __cC; */\
		"vmulpd			%%ymm9 ,%%ymm3 ,%%ymm3 	\n\t"/*	t10 *= __cC; */\
		"vmovaps		%%ymm0 ,     (%%rcx)  	\n\t"/* __B8r = t1 ; */\
		"vmovaps		%%ymm1 ,0x020(%%rcx)  	\n\t"/* __B8i = t2 ; */\
		"vmovaps		%%ymm7 ,     (%%rbx)  	\n\t"/* __B4r = t26; */\
		"vmovaps		%%ymm6 ,0x020(%%rbx)  	\n\t"/* __B4i = t25; */\
		"vmovaps		%%ymm2 ,     (%%rdx)  	\n\t"/* __BCr = t9 ; */\
		"vmovaps		%%ymm3 ,0x020(%%rdx)  	\n\t"/* __BCi = t10; */\
		/**/\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%r14,2)\n\t"/* ...+p8 */\
	"leaq	(%%r14,%%r14,2),%%r14	\n\t"	/* p4 + (p4*2) = p12, ptr-offset form */\
/* Swap Blocks 2,3 here to allow simple-incrementation of array ptrs (decr. not safe w.r.to array padding). */\
	/*...Block 2: t3/4,11/12,19/20,27/28 in ymm0-7: */\
		"addq	$0x40,%%rdi		\n\t"/* t3 */\
		"addq	%%r10,%%rax		\n\t"/* add0+p1 */\
		"addq	%%r10,%%rbx		\n\t"/* add0+p5 */\
		"addq	%%r10,%%rcx		\n\t"/* add0+p9 */\
		"addq	%%r10,%%rdx		\n\t"/* add0+pD */\
		"vmovaps		0x100(%%rdi),%%ymm2 	\n\t	vmovaps		0x120(%%rdi),%%ymm3 	\n\t"/*	t11,12; */\
		"vmovaps		%%ymm2,%%ymm8			\n\t	vbroadcastsd 0x80(%%rsi),%%ymm13	\n\t"/* rt = t11; load __sc */\
		" vfmadd231pd	%%ymm10,%%ymm3 ,%%ymm2 	\n\t"/*	 FMA231(1.0,t12,t11); */\
		"vfnmadd231pd	%%ymm10,%%ymm8 ,%%ymm3 	\n\t"/*	FNMA231(1.0,rt ,t12); */\
		"vmovaps		0x200(%%rdi),%%ymm4 	\n\t	vmovaps		0x220(%%rdi),%%ymm5 	\n\t"/*	t19,20; */\
		"vmovaps		0x300(%%rdi),%%ymm6 	\n\t	vmovaps		0x320(%%rdi),%%ymm7 	\n\t"/*	t27,28; */\
		"vmovaps		%%ymm4,%%ymm11			\n\t	vmovaps		%%ymm6,%%ymm12	\n\t"/*	re = t19; im = t27; */\
		" vfmadd231pd	%%ymm13,%%ymm5 ,%%ymm4 	\n\t"/*	 FMA231(__sc,t20,t19); */\
		"vfnmadd231pd	%%ymm13,%%ymm11,%%ymm5 	\n\t"/*	FNMA231(__sc,re ,t20); */\
		" vfmadd132pd	%%ymm13,%%ymm7 ,%%ymm6 	\n\t"/*	 FMA132(__sc,t27,t28); */\
		" vfmsub132pd	%%ymm13,%%ymm12,%%ymm7 	\n\t"/*	 FMS132(__sc,t28,im ); */\
		"vmovaps		     (%%rdi),%%ymm0 	\n\t	vmovaps		0x020(%%rdi),%%ymm1 	\n\t"/*	t3,4 ; */\
		"vbroadcastsd	-0x08(%%rsi),%%ymm14	\n\t	vbroadcastsd     (%%rsi),%%ymm15	\n\t"/* load ISRT2, __c */\
		"vmovaps		%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9	\n\t"/*	rt = t3; it = t4; */\
		" vfmadd231pd	%%ymm14,%%ymm2 ,%%ymm0 	\n\t"/*	 FMA231(ISRT2,t11,t3 ); */\
		"vfnmadd132pd	%%ymm14,%%ymm8 ,%%ymm2 	\n\t"/*	FNMA132(ISRT2,t11,rt ); */\
		" vfmadd231pd	%%ymm14,%%ymm3 ,%%ymm1 	\n\t"/*	 FMA231(ISRT2,t12,t4 ); */\
		"vfnmadd132pd	%%ymm14,%%ymm9 ,%%ymm3 	\n\t"/*	FNMA132(ISRT2,t12,it ); */\
		"vmulpd			%%ymm15,%%ymm6 ,%%ymm6 	\n\t"/*	t27 *= __c; */\
		"vmulpd			%%ymm15,%%ymm7 ,%%ymm7 	\n\t"/*	t28 *= __c; */\
		"vmovaps		%%ymm6,%%ymm11			\n\t	vmovaps		%%ymm7,%%ymm12	\n\t"/*	re = t27; im = t28; */\
		" vfmsub231pd	%%ymm15,%%ymm4 ,%%ymm6 	\n\t"/*	 FMS231(__c,t19,t27); */\
		" vfmadd132pd	%%ymm15,%%ymm11,%%ymm4 	\n\t"/*	 FMA132(__c,t19,re ); */\
		" vfmsub231pd	%%ymm15,%%ymm5 ,%%ymm7 	\n\t"/*	 FMS231(__c,t20,t28); */\
		" vfmadd132pd	%%ymm15,%%ymm12,%%ymm5 	\n\t"/*	 FMA132(__c,t20,im ); */\
		"vmovaps		%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9	\n\t"/*	rt = t3; it = t4; */\
		" vfmsub132pd	%%ymm10,%%ymm4 ,%%ymm0 	\n\t"/*	 FMS132(1.0,t3 ,t19); */\
		" vfmsub132pd	%%ymm10,%%ymm5 ,%%ymm1 	\n\t"/*	 FMS132(1.0,t4 ,t20); */\
		" vfmadd231pd	%%ymm10,%%ymm8 ,%%ymm4 	\n\t"/*	 FMA231(1.0,rt ,t19); */\
		" vfmadd231pd	%%ymm10,%%ymm9 ,%%ymm5 	\n\t"/*	 FMA231(1.0,it ,t20); */\
		"vbroadcastsd		0x88(%%rsi),%%ymm14	\n\t	vbroadcastsd		0xc8(%%rsi),%%ymm15	\n\t"/* load __t1,9 */\
		"vmovaps		%%ymm2,%%ymm11			\n\t	vmovaps		%%ymm3,%%ymm12	\n\t"/*	re = t11; im = t12; */\
		" vfmsub132pd	%%ymm10,%%ymm7 ,%%ymm2 	\n\t"/*	 FMS132(1.0,t11,t28); */\
		" vfmadd132pd	%%ymm10,%%ymm6 ,%%ymm3 	\n\t"/*	 FMA132(1.0,t12,t27); */\
		" vfmadd231pd	%%ymm10,%%ymm11,%%ymm7 	\n\t"/*	 FMA231(1.0,re ,t28); */\
		" vfmsub231pd	%%ymm10,%%ymm12,%%ymm6 	\n\t"/*	 FMS231(1.0,im ,t27); */\
		"vmovaps		%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm0,%%ymm9	\n\t"/*	rt = t19; it = t3; */\
		"vbroadcastsd		0xa8(%%rsi),%%ymm10	\n\t	vbroadcastsd		0xe8(%%rsi),%%ymm13	\n\t"/* load __t5,D; use 1.0 slot for __t5 */\
		" vfmadd231pd	%%ymm14,%%ymm5 ,%%ymm4 	\n\t"/*	 FMA231(__t1,t20,t19); */\
		"vfnmadd231pd	%%ymm14,%%ymm8 ,%%ymm5 	\n\t"/*	FNMA231(__t1,rt ,t20); */\
		" vfmadd231pd	%%ymm15,%%ymm1 ,%%ymm0 	\n\t"/*	 FMA231(__t9,t4 ,t3 ); */\
		"vfnmadd231pd	%%ymm15,%%ymm9 ,%%ymm1 	\n\t"/*	FNMA231(__t9,it ,t4 ); */\
		"vbroadcastsd		0x08(%%rsi),%%ymm14	\n\t	vbroadcastsd		0x48(%%rsi),%%ymm15	\n\t"/* load __c1,9 */\
		"vmovaps		%%ymm7,%%ymm11			\n\t	vmovaps		%%ymm2,%%ymm12	\n\t"/*	re = t28; im = t11; */\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm7 	\n\t"/*	 FMA231(__t5,t27,t28); */\
		"vfnmadd231pd	%%ymm10,%%ymm11,%%ymm6 	\n\t"/*	FNMA231(__t5,re ,t27); */\
		" vfmadd231pd	%%ymm13,%%ymm3 ,%%ymm2 	\n\t"/*	 FMA231(__tD,t12,t11); */\
		"vfnmadd231pd	%%ymm13,%%ymm12,%%ymm3 	\n\t"/*	FNMA231(__tD,im ,t12); */\
		"vbroadcastsd		0x28(%%rsi),%%ymm10	\n\t	vbroadcastsd		0x68(%%rsi),%%ymm13	\n\t"/* load __c5,D; use 1.0 slot for __t5 */\
		"vmulpd			%%ymm14,%%ymm4 ,%%ymm4 	\n\t"/*	t19 *= __c1; */\
		"vmulpd			%%ymm14,%%ymm5 ,%%ymm5 	\n\t"/*	t20 *= __c1; */\
		"vmulpd			%%ymm15,%%ymm0 ,%%ymm0 	\n\t"/*	t3  *= __c9; */\
		"vmulpd			%%ymm15,%%ymm1 ,%%ymm1 	\n\t"/*	t4  *= __c9; */\
		"vmulpd			%%ymm10,%%ymm7 ,%%ymm7 	\n\t"/*	t28 *= __c5; */\
		"vmulpd			%%ymm10,%%ymm6 ,%%ymm6 	\n\t"/*	t27 *= __c5; */\
		"vmulpd			%%ymm13,%%ymm2 ,%%ymm2 	\n\t"/*	t11 *= __cD; */\
		"vmulpd			%%ymm13,%%ymm3 ,%%ymm3 	\n\t"/*	t12 *= __cD; */\
		"vmovaps		%%ymm4 ,     (%%rax)  	\n\t"/* __B1r = t19; */\
		"vmovaps		%%ymm5 ,0x020(%%rax)  	\n\t"/* __B1i = t20; */\
		"vmovaps		%%ymm0 ,     (%%rcx)  	\n\t"/* __B9r = t3 ; */\
		"vmovaps		%%ymm1 ,0x020(%%rcx)  	\n\t"/* __B9i = t4 ;	Since cannot assume p2 = (p1+p1) due to array padding, revert:*/\
		"vmovaps		%%ymm7 ,     (%%rbx)  	\n\t	subq	%%r10,%%rax	\n\t"/* __B5r = t28; add0+p1 -> add0+p0 */\
		"vmovaps		%%ymm6 ,0x020(%%rbx)  	\n\t	subq	%%r10,%%rbx	\n\t"/* __B5i = t27; add0+p5 -> add0+p4 */\
		"vmovaps		%%ymm2 ,     (%%rdx)  	\n\t	subq	%%r10,%%rcx	\n\t"/* __BDr = t11; add0+p9 -> add0+p8 */\
		"vmovaps		%%ymm3 ,0x020(%%rdx)  	\n\t	subq	%%r10,%%rdx	\n\t"/* __BDi = t12; add0+pD -> add0+pC */\
		/**/\
	/*...Block 3: t5/6,13/14,21/22,29/30 in ymm0-7: */\
		"vmovaps 0x440(%%rsi),%%ymm10	\n\t"/* retore 1.0 to ymm10 */\
		"addq	$0x40,%%rdi		\n\t"/* t5 */\
		"movslq	%[__p2],%%r9	\n\t"\
		"shlq	$3,%%r9			\n\t"\
		"addq	%%r9,%%rax		\n\t"/* add0+p2 */\
		"addq	%%r9,%%rbx		\n\t"/* add0+p6 */\
		"addq	%%r9,%%rcx		\n\t"/* add0+pA */\
		"addq	%%r9,%%rdx		\n\t"/* add0+pE */\
		"vmovaps		     (%%rdi),%%ymm0 	\n\t	vmovaps		0x100(%%rdi),%%ymm2 	\n\t"/*	t5,13; */\
		"vmovaps		0x020(%%rdi),%%ymm1 	\n\t	vmovaps		0x120(%%rdi),%%ymm3 	\n\t"/*	t6,14; */\
		"vmovaps		%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9	\n\t"/*	rt = t5; it = t6; */\
		" vfmadd231pd	%%ymm10,%%ymm3 ,%%ymm0 	\n\t"/*	 FMA231(1.0,t14,t5 ); */\
		"vfnmadd132pd	%%ymm10,%%ymm8 ,%%ymm3 	\n\t"/*	FNMA132(1.0,t14,rt ); */\
		"vfnmadd231pd	%%ymm10,%%ymm2 ,%%ymm1 	\n\t"/*	FNMA231(1.0,t13,t6 ); */\
		" vfmadd132pd	%%ymm10,%%ymm9 ,%%ymm2 	\n\t"/*	 FMA132(1.0,t13,it ); */\
		"vmovaps		0x200(%%rdi),%%ymm4 	\n\t	vmovaps		0x220(%%rdi),%%ymm5 	\n\t"/*	t21,22; */\
		"vmovaps		0x300(%%rdi),%%ymm6 	\n\t	vmovaps		0x320(%%rdi),%%ymm7 	\n\t"/*	t29,30; */\
		"vmovaps		%%ymm5,%%ymm11			\n\t	vmovaps		%%ymm7,%%ymm12	\n\t"/*	re = t22; im = t30; */\
		"vfnmadd231pd	%%ymm10,%%ymm4 ,%%ymm5 	\n\t"/*	FNMA231(1.0,t21,t22); */\
		" vfmadd132pd	%%ymm10,%%ymm11,%%ymm4 	\n\t"/*	 FMA132(1.0,t21,re ); */\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm7 	\n\t"/*	 FMA231(1.0,t29,t30); */\
		" vfmsub132pd	%%ymm10,%%ymm12,%%ymm6 	\n\t"/*	 FMS132(1.0,t29,im ); */\
		"vmovaps		%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm5,%%ymm9	\n\t"/*	rt = t21; it = t22; */\
		"vfnmadd231pd	%%ymm10,%%ymm6 ,%%ymm4 	\n\t"/*	FNMA231(1.0,t29,t21); */\
		" vfmadd132pd	%%ymm10,%%ymm8 ,%%ymm6 	\n\t"/*	 FMA132(1.0,t29,rt ); */\
		"vfnmadd231pd	%%ymm10,%%ymm7 ,%%ymm5 	\n\t"/*	FNMA231(1.0,t30,t22); */\
		" vfmadd132pd	%%ymm10,%%ymm9 ,%%ymm7 	\n\t"/*	 FMA132(1.0,t30,it ); */\
		"vbroadcastsd      -0x08(%%rsi),%%ymm13	\n\t"/* load ISRT2 */\
		"vmovaps		%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9	\n\t"/*	rt = t5; it = t6; */\
		"vfnmadd231pd	%%ymm13,%%ymm4 ,%%ymm0 	\n\t"/*	FNMA231(ISRT2,t21,t5 ); */\
		" vfmadd132pd	%%ymm13,%%ymm8 ,%%ymm4 	\n\t"/*	 FMA132(ISRT2,t21,rt ); */\
		"vfnmadd231pd	%%ymm13,%%ymm5 ,%%ymm1 	\n\t"/*	FNMA231(ISRT2,t22,t6 ); */\
		" vfmadd132pd	%%ymm13,%%ymm9 ,%%ymm5 	\n\t"/*	 FMA132(ISRT2,t22,it ); */\
		"vbroadcastsd		0x90(%%rsi),%%ymm14	\n\t	vbroadcastsd		0xd0(%%rsi),%%ymm15	\n\t"/* load __t2,A */\
		"vmovaps		%%ymm2,%%ymm11			\n\t	vmovaps		%%ymm3,%%ymm12	\n\t"/*	re =t13; im =t14; */\
		" vfmadd231pd	%%ymm13,%%ymm6 ,%%ymm2 	\n\t"/*	 FMA231(ISRT2,t29,t13); */\
		"vfnmadd132pd	%%ymm13,%%ymm11,%%ymm6 	\n\t"/*	FNMA132(ISRT2,t29,re ); */\
		"vfnmadd231pd	%%ymm13,%%ymm7 ,%%ymm3 	\n\t"/*	FNMA231(ISRT2,t30,t14); */\
		" vfmadd132pd	%%ymm13,%%ymm12,%%ymm7 	\n\t"/*	 FMA132(ISRT2,t30,im ); */\
		"vmovaps		%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm0,%%ymm9	\n\t"/*	rt =t21; it = t5; */\
		"vbroadcastsd		0xb0(%%rsi),%%ymm10	\n\t	vbroadcastsd		0xf0(%%rsi),%%ymm13	\n\t"/* load __t6,E; use 1.0 slot for __t5 */\
		" vfmadd231pd	%%ymm14,%%ymm5 ,%%ymm4 	\n\t"/*	 FMA231(__t2,t22,t21); */\
		"vfnmadd231pd	%%ymm14,%%ymm8 ,%%ymm5 	\n\t"/*	FNMA231(__t2,rt ,t22); */\
		" vfmadd231pd	%%ymm15,%%ymm1 ,%%ymm0 	\n\t"/*	 FMA231(__tA,t6 ,t5 ); */\
		"vfnmadd231pd	%%ymm15,%%ymm9 ,%%ymm1 	\n\t"/*	FNMA231(__tA,it ,t6 ); */\
		"vbroadcastsd		0x10(%%rsi),%%ymm14	\n\t	vbroadcastsd		0x50(%%rsi),%%ymm15	\n\t"/* load __c2,A */\
		"vmovaps		%%ymm7,%%ymm11			\n\t	vmovaps		%%ymm3,%%ymm12	\n\t"/*	re =t30; im =t14; */\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm7 	\n\t"/*	 FMA231(__t6,t29,t30); */\
		"vfnmadd231pd	%%ymm10,%%ymm11,%%ymm6 	\n\t"/*	FNMA231(__t6,re ,t29); */\
		" vfmadd231pd	%%ymm13,%%ymm2 ,%%ymm3 	\n\t"/*	 FMA231(__tE,t13,t14); */\
		"vfnmadd231pd	%%ymm13,%%ymm12,%%ymm2 	\n\t"/*	FNMA231(__tE,im ,t13); */\
		"vbroadcastsd		0x30(%%rsi),%%ymm10	\n\t	vbroadcastsd		0x70(%%rsi),%%ymm13	\n\t"/* load __c6,E; use 1.0 slot for __t5 */\
		"vmulpd			%%ymm14,%%ymm4 ,%%ymm4 	\n\t"/*	t21 *= __c2; */\
		"vmulpd			%%ymm14,%%ymm5 ,%%ymm5 	\n\t"/*	t22 *= __c2; */\
		"vmulpd			%%ymm15,%%ymm0 ,%%ymm0 	\n\t"/*	t5  *= __cA; */\
		"vmulpd			%%ymm15,%%ymm1 ,%%ymm1 	\n\t"/*	t6  *= __cA; */\
		"vmulpd			%%ymm10,%%ymm7 ,%%ymm7 	\n\t"/*	t30 *= __c6; */\
		"vmulpd			%%ymm10,%%ymm6 ,%%ymm6 	\n\t"/*	t29 *= __c6; */\
		"vmulpd			%%ymm13,%%ymm3 ,%%ymm3 	\n\t"/*	t14 *= __cE; */\
		"vmulpd			%%ymm13,%%ymm2 ,%%ymm2 	\n\t"/*	t13 *= __cE; */\
		"vmovaps		%%ymm4 ,     (%%rax)  	\n\t"/* __B2r = t21; */\
		"vmovaps		%%ymm5 ,0x020(%%rax)  	\n\t"/* __B2i = t22; */\
		"vmovaps		%%ymm0 ,     (%%rcx)  	\n\t"/* __BAr = t5 ; */\
		"vmovaps		%%ymm1 ,0x020(%%rcx)  	\n\t"/* __BAi = t6 ; */\
		"vmovaps		%%ymm7 ,     (%%rbx)  	\n\t"/* __B6r = t30; */\
		"vmovaps		%%ymm6 ,0x020(%%rbx)  	\n\t"/* __B6i = t29; */\
		"vmovaps		%%ymm3 ,     (%%rdx)  	\n\t"/* __BEr = t14; */\
		"vmovaps		%%ymm2 ,0x020(%%rdx)  	\n\t"/* __BEi = t13; */\
		/**/\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%r14)\n\t"/* ...+p12 */\
	/*...Block 4: t7/8,15/16,23/24,31/32 in ymm0-7: */\
		"vmovaps 0x440(%%rsi),%%ymm10	\n\t"/* retore 1.0 to ymm10 */\
		"addq	$0x40,%%rdi		\n\t"/* t7 */\
		"addq	%%r10,%%rax		\n\t"/* add0+p3 */\
		"addq	%%r10,%%rbx		\n\t"/* add0+p7 */\
		"addq	%%r10,%%rcx		\n\t"/* add0+pB */\
		"addq	%%r10,%%rdx		\n\t"/* add0+pF */\
		"vmovaps		0x100(%%rdi),%%ymm2 	\n\t	vmovaps		0x120(%%rdi),%%ymm3 	\n\t"/*	t15,16; */\
		"vmovaps		%%ymm3,%%ymm8			\n\t	vbroadcastsd 0x80(%%rsi),%%ymm13	\n\t"/* rt = t16; load __sc */\
		" vfmadd231pd	%%ymm10,%%ymm2 ,%%ymm3 	\n\t"/*	 FMA231(1.0,t15,t16); */\
		" vfmsub132pd	%%ymm10,%%ymm8 ,%%ymm2 	\n\t"/*	 FMS132(1.0,t15,rt ); */\
		"vmovaps		0x200(%%rdi),%%ymm4 	\n\t	vmovaps		0x220(%%rdi),%%ymm5 	\n\t"/*	t23,24; */\
		"vmovaps		0x300(%%rdi),%%ymm6 	\n\t	vmovaps		0x320(%%rdi),%%ymm7 	\n\t"/*	t31,32; */\
		"vmovaps		%%ymm4,%%ymm11			\n\t	vmovaps		%%ymm6,%%ymm12	\n\t"/*	re = t23; im = t31; */\
		" vfmadd231pd	%%ymm13,%%ymm7 ,%%ymm6 	\n\t"/*	 FMA231(__sc,t32,t31); */\
		"vfnmadd231pd	%%ymm13,%%ymm12,%%ymm7 	\n\t"/*	FNMA231(__sc,im ,t32); */\
		" vfmadd132pd	%%ymm13,%%ymm5 ,%%ymm4 	\n\t"/*	 FMA132(__sc,t23,t24); */\
		" vfmsub132pd	%%ymm13,%%ymm11,%%ymm5 	\n\t"/*	 FMS132(__sc,t24,re ); */\
		"vbroadcastsd	-0x08(%%rsi),%%ymm14	\n\t	vbroadcastsd     (%%rsi),%%ymm15	\n\t"/* load ISRT2, __c */\
		"vmovaps		     (%%rdi),%%ymm0 	\n\t	vmovaps		0x020(%%rdi),%%ymm1 	\n\t"/*	t7,8 ; */\
		"vmovaps		%%ymm0,%%ymm8			\n\t	vmovaps		%%ymm1,%%ymm9	\n\t"/*	rt = t7; it = t8; */\
		"vfnmadd231pd	%%ymm14,%%ymm2 ,%%ymm0 	\n\t"/*	FNMA231(ISRT2,t15,t7 ); */\
		" vfmadd132pd	%%ymm14,%%ymm8 ,%%ymm2 	\n\t"/*	 FMA132(ISRT2,t15,rt ); */\
		"vfnmadd231pd	%%ymm14,%%ymm3 ,%%ymm1 	\n\t"/*	FNMA231(ISRT2,t16,t8 ); */\
		" vfmadd132pd	%%ymm14,%%ymm9 ,%%ymm3 	\n\t"/*	 FMA132(ISRT2,t16,it ); */\
		"vmovaps		%%ymm6,%%ymm11			\n\t	vmovaps		%%ymm7,%%ymm12	\n\t"/*	re = t31; im = t32; */\
		" vfmadd231pd	%%ymm10,%%ymm4 ,%%ymm6 	\n\t"/*	 FMA231(1.0,t23,t31); */\
		" vfmsub132pd	%%ymm10,%%ymm11,%%ymm4 	\n\t"/*	 FMS132(1.0,t23,re ); */\
		" vfmadd231pd	%%ymm10,%%ymm5 ,%%ymm7 	\n\t"/*	 FMA231(1.0,t24,t32); */\
		" vfmsub132pd	%%ymm10,%%ymm12,%%ymm5 	\n\t"/*	 FMS132(1.0,t24,im ); */\
		"vmovaps		%%ymm15,%%ymm10	\n\t"/*	copy __c into 1.0-holding register; */\
		"vmovaps		%%ymm2,%%ymm8			\n\t	vmovaps		%%ymm3,%%ymm9	\n\t"/*	rt = t15; it = t16; */\
		"vfnmadd231pd	%%ymm10,%%ymm7 ,%%ymm2 	\n\t"/*	FNMA231(__c,t32,t15); */\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm3 	\n\t"/*	 FMA231(__c,t31,t16); */\
		" vfmadd132pd	%%ymm10,%%ymm8 ,%%ymm7 	\n\t"/*	 FMA132(__c,t32,rt ); */\
		"vfnmadd132pd	%%ymm10,%%ymm9 ,%%ymm6 	\n\t"/*	FNMA132(__c,t31,it ); */\
		"vbroadcastsd		0x98(%%rsi),%%ymm14	\n\t	vbroadcastsd		0xd8(%%rsi),%%ymm15	\n\t"/* load __t3,B */\
		"vmovaps		%%ymm0,%%ymm11			\n\t	vmovaps		%%ymm1,%%ymm12	\n\t"/*	re = t7; im = t8; */\
		"vfnmadd231pd	%%ymm10,%%ymm4 ,%%ymm0 	\n\t"/*	FNMA231(__c,t23,t7 ); */\
		"vfnmadd231pd	%%ymm10,%%ymm5 ,%%ymm1 	\n\t"/*	FNMA231(__c,t24,t8 ); */\
		" vfmadd132pd	%%ymm10,%%ymm11,%%ymm4 	\n\t"/*	 FMA132(__c,t23,re ); */\
		" vfmadd132pd	%%ymm10,%%ymm12,%%ymm5 	\n\t"/*	 FMA132(__c,t24,im ); */\
		"vmovaps		%%ymm4,%%ymm8			\n\t	vmovaps		%%ymm0,%%ymm9	\n\t"/*	rt = t23; it = t7; */\
		"vbroadcastsd		0xb8(%%rsi),%%ymm10	\n\t	vbroadcastsd		0xf8(%%rsi),%%ymm13	\n\t"/* load __t7,F; use 1.0 slot for __t7 */\
		" vfmadd231pd	%%ymm14,%%ymm5 ,%%ymm4 	\n\t"/*	 FMA231(__t3,t24,t23); */\
		"vfnmadd231pd	%%ymm14,%%ymm8 ,%%ymm5 	\n\t"/*	FNMA231(__t3,rt ,t24); */\
		" vfmadd231pd	%%ymm15,%%ymm1 ,%%ymm0 	\n\t"/*	 FMA231(__tB,t8 ,t7 ); */\
		"vfnmadd231pd	%%ymm15,%%ymm9 ,%%ymm1 	\n\t"/*	FNMA231(__tB,it ,t8 ); */\
		"vbroadcastsd		0x18(%%rsi),%%ymm14	\n\t	vbroadcastsd		0x58(%%rsi),%%ymm15	\n\t"/* load __c3,B */\
		"vmovaps		%%ymm7,%%ymm11			\n\t	vmovaps		%%ymm2,%%ymm12	\n\t"/*	re = t32; im = t15; */\
		" vfmadd231pd	%%ymm10,%%ymm6 ,%%ymm7 	\n\t"/*	 FMA231(__t7,t31,t32); */\
		"vfnmadd231pd	%%ymm10,%%ymm11,%%ymm6 	\n\t"/*	FNMA231(__t7,re ,t31); */\
		" vfmadd231pd	%%ymm13,%%ymm3 ,%%ymm2 	\n\t"/*	 FMA231(__tF,t16,t15); */\
		"vfnmadd231pd	%%ymm13,%%ymm12,%%ymm3 	\n\t"/*	FNMA231(__tF,im ,t16); */\
		"vbroadcastsd		0x38(%%rsi),%%ymm10	\n\t	vbroadcastsd		0x78(%%rsi),%%ymm13	\n\t"/* load __c7,F; use 1.0 slot for __t5 */\
		"vmulpd			%%ymm14,%%ymm4 ,%%ymm4 	\n\t"/*	t23 *= __c3; */\
		"vmulpd			%%ymm14,%%ymm5 ,%%ymm5 	\n\t"/*	t24 *= __c3; */\
		"vmulpd			%%ymm15,%%ymm0 ,%%ymm0 	\n\t"/*	t7  *= __cB; */\
		"vmulpd			%%ymm15,%%ymm1 ,%%ymm1 	\n\t"/*	t8  *= __cB; */\
		"vmulpd			%%ymm10,%%ymm7 ,%%ymm7 	\n\t"/*	t32 *= __c7; */\
		"vmulpd			%%ymm10,%%ymm6 ,%%ymm6 	\n\t"/*	t31 *= __c7; */\
		"vmulpd			%%ymm13,%%ymm2 ,%%ymm2 	\n\t"/*	t15 *= __cF; */\
		"vmulpd			%%ymm13,%%ymm3 ,%%ymm3 	\n\t"/*	t16 *= __cF; */\
		"vmovaps		%%ymm4 ,     (%%rax)  	\n\t"/* __B3r = t23; */\
		"vmovaps		%%ymm5 ,0x020(%%rax)  	\n\t"/* __B3i = t24; */\
		"vmovaps		%%ymm0 ,     (%%rcx)  	\n\t"/* __BBr = t7 ; */\
		"vmovaps		%%ymm1 ,0x020(%%rcx)  	\n\t"/* __BBi = t8 ; */\
		"vmovaps		%%ymm7 ,     (%%rbx)  	\n\t"/* __B7r = t32; */\
		"vmovaps		%%ymm6 ,0x020(%%rbx)  	\n\t"/* __B7i = t31; */\
		"vmovaps		%%ymm2 ,     (%%rdx)  	\n\t"/* __BFr = t15; */\
		"vmovaps		%%ymm3 ,0x020(%%rdx)  	\n\t"/* __BFi = t16; */\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r1] "m" (Xr1)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r9","r10","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	// Default is just an FMA-greedy-optimized version of the non-FMA AVX-based DIT:
	// 82 ADD, 102 FMA, 48 pure-MUL, 216 MOVAPS
	#define SSE2_RADIX16_DIT_TWIDDLE_0(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xr1,Xisrt2,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
		"movq	%[__isrt2],%%r8					\n\t"\
		"movq	%[__add0],%%rax					\n\t"\
	"movslq	%[__pfetch_addr1],%%r13	\n\t"	/* Prefetch base-index offset cycles among add0 + p0,1,2,3 on successive macro calls */\
	"leaq	(%%rax,%%r13,8),%%r13	\n\t"	/* [base-address + data-fetch-ahead index] */\
		"movslq	%[__p1],%%rbx					\n\t"\
		"movslq	%[__p2],%%rcx					\n\t"\
		"movslq	%[__p3],%%rdx					\n\t"\
		"movslq	%[__p4],%%rdi					\n\t		movslq	%[__p8],%%r10		\n\t"\
		"shlq	$3,%%rbx						\n\t		shlq	$3,%%r10			\n\t"\
		"shlq	$3,%%rcx						\n\t		movq	%%r10,%%r11			\n\t"\
		"shlq	$3,%%rdx						\n\t		movq	%%r10,%%r12			\n\t"\
		"shlq	$3,%%rdi						\n\t		movq	%%r10,%%r14			\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r13)\n\t"/* [base-address + data-fetch-ahead index] + p0 */\
		"addq	$0x460,%%r8	/* two */			\n\t		addq	%%rax,%%r10			\n\t"\
		"addq	%%rax,%%rbx						\n\t		addq	%%rbx,%%r11			\n\t"\
		"addq	%%rax,%%rcx						\n\t		addq	%%rcx,%%r12			\n\t"\
		"addq	%%rax,%%rdx						\n\t		addq	%%rdx,%%r14			\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r1) */	\n\t		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r17) */\n\t"\
		"movq	%[__r1],%%rsi					\n\t"/* r17 at rsi + 0x200 */\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	     (%%r11),%%ymm8 	\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	     (%%r14),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		vmovaps	0x020(%%r11),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		vmovaps	0x020(%%r14),%%ymm13	\n\t"\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	-0x400(%%r8),%%ymm10	\n\t"/* load 1.0 and */\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vmovaps	%%ymm10,%%ymm14			\n\t"/* propagate copies */\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	%%ymm10,%%ymm11			\n\t"/* to all 4 output */\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vmovaps	%%ymm10,%%ymm15			\n\t"/* registers of FMAs: */\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2			\n\t		vfnmadd213pd	     (%%r10),%%ymm8 ,%%ymm10\n\t"/* FNMA_213: MUL rightmost 2 ops, */\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t		vfnmadd213pd	     (%%r12),%%ymm12,%%ymm14\n\t"/* [i.e. 1.0*ymm8/10/14/11/15], SUB */\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3			\n\t		vfnmadd213pd	0x020(%%r10),%%ymm9 ,%%ymm11\n\t"/* product from left memory operand, */\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t		vfnmadd213pd	0x020(%%r12),%%ymm13,%%ymm15\n\t"/* result overwrites rightmost operand. */\
		" vfmadd132pd	(%%r8),%%ymm2,%%ymm0	\n\t"\
		" vfmadd132pd	(%%r8),%%ymm6,%%ymm4	\n\t"\
		" vfmadd132pd	(%%r8),%%ymm3,%%ymm1	\n\t		 vfmadd132pd	(%%r8),%%ymm10,%%ymm8 		\n\t"/* FMA_132: MUL the bracketing 2 ops */\
		" vfmadd132pd	(%%r8),%%ymm7,%%ymm5	\n\t		 vfmadd132pd	(%%r8),%%ymm14,%%ymm12		\n\t"/* [i.e. 2.0*ymm8/12/9/13], ADD */\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t		 vfmadd132pd	(%%r8),%%ymm11,%%ymm9 		\n\t"/* product to middle register operand, */\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t		 vfmadd132pd	(%%r8),%%ymm15,%%ymm13		\n\t"/* result overwrites rightmost operand. */\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)			\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm2,0x0c0(%%rsi)			\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)			\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)			\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		" vfmadd132pd	(%%r8),%%ymm0,%%ymm4	\n\t		vmovaps	%%ymm8 ,0x280(%%rsi)\n\t"\
		" vfmadd132pd	(%%r8),%%ymm2,%%ymm7	\n\t		vmovaps	%%ymm10,0x2c0(%%rsi)\n\t"\
		" vfmadd132pd	(%%r8),%%ymm1,%%ymm5	\n\t		vmovaps	%%ymm9 ,0x2a0(%%rsi)\n\t"\
		" vfmadd132pd	(%%r8),%%ymm3,%%ymm6	\n\t		vmovaps	%%ymm11,0x260(%%rsi)\n\t"\
		"vmovaps	%%ymm4,(%%rsi) \n\t vmovaps (%%r8),%%ymm4 \n\t	 vfmadd132pd %%ymm4,%%ymm8 ,%%ymm12	\n\t"\
		"vmovaps	%%ymm7,0x040(%%rsi)			\n\t				 vfmadd132pd %%ymm4,%%ymm10,%%ymm15	\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)			\n\t				 vfmadd132pd %%ymm4,%%ymm9 ,%%ymm13	\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)			\n\t				 vfmadd132pd %%ymm4,%%ymm11,%%ymm14	\n\t"\
		"addq	%%rdi,%%rax						\n\t		vmovaps	%%ymm12,0x200(%%rsi)\n\t"\
		"addq	%%rdi,%%rbx						\n\t		vmovaps	%%ymm15,0x240(%%rsi)\n\t"\
		"addq	%%rdi,%%rcx						\n\t		vmovaps	%%ymm13,0x220(%%rsi)\n\t"\
		"addq	%%rdi,%%rdx						\n\t		vmovaps	%%ymm14,0x2e0(%%rsi)\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r9) */	\n\t		addq	%%rdi,%%r10			\n\t"\
		"addq	$0x100,%%rsi		/* r9 */	\n\t		addq	%%rdi,%%r11			\n\t"\
		"vmovaps	     (%%rax),%%ymm2			\n\t		addq	%%rdi,%%r12			\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		addq	%%rdi,%%r14			\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r25) */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t"		/* r25 at rsi + 0x200 */\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 	\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		vmovaps	     (%%r14),%%ymm12	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		vmovaps	0x020(%%r14),%%ymm13	\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%rdi)\n\t"/* ...+p4 */\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2			\n\t		vmovaps	-0x400(%%r8),%%ymm10	\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t		vmovaps	%%ymm10,%%ymm14			\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3			\n\t		vmovaps	%%ymm10,%%ymm11			\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t		vmovaps	%%ymm10,%%ymm15			\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0			\n\t		vfnmadd213pd	     (%%r10),%%ymm8 ,%%ymm10\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4			\n\t		vfnmadd213pd	     (%%r12),%%ymm12,%%ymm14\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t		vfnmadd213pd	0x020(%%r10),%%ymm9 ,%%ymm11\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5			\n\t		vfnmadd213pd	0x020(%%r12),%%ymm13,%%ymm15\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0			\n\t		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4			\n\t		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1			\n\t		 vfmadd132pd	(%%r8),%%ymm10,%%ymm8 		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5			\n\t		 vfmadd132pd	(%%r8),%%ymm14,%%ymm12		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t		 vfmadd132pd	(%%r8),%%ymm11,%%ymm9 		\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t		 vfmadd132pd	(%%r8),%%ymm15,%%ymm13		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)			\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm2,0x0c0(%%rsi)			\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)			\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)			\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4			\n\t		vmovaps	%%ymm8 ,0x280(%%rsi)\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7			\n\t		vmovaps	%%ymm10,0x2c0(%%rsi)\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5			\n\t		vmovaps	%%ymm9 ,0x2a0(%%rsi)\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6			\n\t		vmovaps	%%ymm11,0x260(%%rsi)\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4			\n\t		 vfmadd132pd	(%%r8),%%ymm8 ,%%ymm12		\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7			\n\t		 vfmadd132pd	(%%r8),%%ymm10,%%ymm15		\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5			\n\t		 vfmadd132pd	(%%r8),%%ymm9 ,%%ymm13		\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6			\n\t		 vfmadd132pd	(%%r8),%%ymm11,%%ymm14		\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)			\n\t		vmovaps	%%ymm12,0x200(%%rsi)\n\t"\
		"vmovaps	%%ymm7,0x040(%%rsi)			\n\t		vmovaps	%%ymm15,0x240(%%rsi)\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)			\n\t		vmovaps	%%ymm13,0x220(%%rsi)\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)			\n\t		vmovaps	%%ymm14,0x2e0(%%rsi)\n\t"\
	/***************************************************************************************************/\
	/*..and now do four more radix-4 transforms, including the internal and external twiddle factors:  */\
	/***************************************************************************************************/\
	/* %%rsi contains address of r9 here */\
	"/*...Block 1: 24 MOVapd, 26 ADD/SUBpd, 12 MULpd */	/*...Block 3: 31 MOVapd, 32 ADD/SUBpd, 20 MULpd */\n\t"\
		"movq	%[__add0],%%rax							\n\t	movslq	%[__p2],%%r10			\n\t"\
		"movslq	%[__p4],%%rbx							\n\t	leaq	-0x080(%%rsi),%%r12	/* r5  */\n\t"\
		"movslq	%[__p8],%%rdi							\n\t"	/* r13 in %%r12 + 0x100*/\
		"shlq	$3,%%rbx								\n\t	movq	%%rbx,%%r11		/* p4 */\n\t"\
		"shlq	$3,%%rdi								\n\t	shlq	$3,%%r10		/* p2 */\n\t"\
		"movq	%[__r1],%%rcx							\n\t	vmovaps	-0x460(%%r8),%%ymm10	/* isrt2 */\n\t"\
		"movq	%%rsi,%%rdx		/* r9 */				\n\t	vmovaps	0x200(%%r12),%%ymm12	\n\t"\
		"			addq	%%rax,%%r10	/* a[j+p2 ] */	\n\t	vmovaps	0x220(%%r12),%%ymm13	\n\t"\
		"addq	%%rax,%%rbx	/* a[j+p4 ] */				\n\t	vmovaps	0x300(%%r12),%%ymm14	\n\t"\
		"			addq	%%r10,%%r11	/* a[j+p6 ] */	\n\t	vmovaps	0x320(%%r12),%%ymm15	\n\t"\
		"vmovaps	     (%%rdx),%%ymm2					\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	0x200(%%rdx),%%ymm4					\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3					\n\t	vmulpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	0x220(%%rdx),%%ymm5					\n\t	vmulpd	%%ymm10,%%ymm15,%%ymm15			\n\t"\
		"vmovaps	     (%%rcx),%%ymm0					\n\t	vmovaps	     (%%r12),%%ymm8 		\n\t"\
		"vmovaps	0x200(%%rcx),%%ymm6					\n\t	vmovaps	0x120(%%r12),%%ymm10	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t	vmovaps	0x020(%%r12),%%ymm11	\n\t"\
		"vmovaps	0x220(%%rcx),%%ymm7					\n\t	vmovaps	0x100(%%r12),%%ymm9 		\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm12,%%ymm13,%%ymm13			\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6					\n\t	vsubpd	%%ymm15,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7					\n\t	vsubpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		" vfmadd132pd	(%%r8),%%ymm0,%%ymm2			\n\t	 vfmadd132pd	(%%r8),%%ymm13,%%ymm12	\n\t"\
		" vfmadd132pd	(%%r8),%%ymm6,%%ymm4			\n\t	 vfmadd132pd	(%%r8),%%ymm14,%%ymm15	\n\t"\
		" vfmadd132pd	(%%r8),%%ymm1,%%ymm3			\n\t	\n\t"\
		" vfmadd132pd	(%%r8),%%ymm7,%%ymm5			\n\t	\n\t"\
	"addq %%rdi,%%r13	\n\t"/* ...+p8 */\
	"prefetcht1	%c[__pfetch_dist](%%r13)\n\t"/* ...+p8 */\
	"addq %%rbx,%%r13	\n\t"/* ...+p12 */\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c0) */	\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"addq	$0x360,%%rsi	/* c0, from r9 */		\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13			\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t	 vfmadd132pd	 (%%r8),%%ymm8 ,%%ymm10	\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0					\n\t	 vfmadd132pd	 (%%r8),%%ymm11,%%ymm9 	\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t	 vfmadd132pd	 (%%r8),%%ymm12,%%ymm14	\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1					\n\t	 vfmadd132pd	 (%%r8),%%ymm13,%%ymm15	\n\t"\
		" vfmadd132pd	(%%r8),%%ymm2,%%ymm4			\n\t	/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c2) */\n\t"\
		" vfmadd132pd	(%%r8),%%ymm0,%%ymm7			\n\t	addq	$0x260,%%r14	/* c2, from r25 */\n\t"\
		" vfmadd132pd	(%%r8),%%ymm3,%%ymm5			\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11			\n\t"\
		" vfmadd132pd	(%%r8),%%ymm1,%%ymm6			\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)					\n\t	vsubpd	%%ymm14,%%ymm9 ,%%ymm9 			\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)					\n\t	vsubpd	%%ymm15,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)	\n\t vmovaps (%%r8),%%ymm3 \n\t	 vfmadd132pd %%ymm3,%%ymm10,%%ymm12	\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)					\n\t			 vfmadd132pd %%ymm3,%%ymm8 ,%%ymm15	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t			 vfmadd132pd %%ymm3,%%ymm11,%%ymm13	\n\t"\
		"vmovaps	%%ymm7,%%ymm0						\n\t			 vfmadd132pd %%ymm3,%%ymm9 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t	vmovaps	%%ymm10,     (%%r12)	\n\t"\
		"vmovaps	%%ymm1,%%ymm6						\n\t	vmovaps	%%ymm8 ,0x120(%%r12)	\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4				\n\t	vmovaps	%%ymm11,0x020(%%r12)	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm7,%%ymm7				\n\t	vmovaps	%%ymm14,0x100(%%r12)	\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5				\n\t	vmovaps	%%ymm12,%%ymm10			\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm1,%%ymm1				\n\t	vmovaps	%%ymm15,%%ymm8 			\n\t"\
		" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4		\n\t	vmovaps	%%ymm13,%%ymm11			\n\t"\
		" vfmadd231pd	0x060(%%rsi),%%ymm6,%%ymm7		\n\t	vmovaps	%%ymm9 ,%%ymm14			\n\t"\
		"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5		\n\t	vmulpd	     (%%r14),%%ymm12,%%ymm12		\n\t"\
		"vfnmadd231pd	0x060(%%rsi),%%ymm0,%%ymm1		\n\t	vmulpd	0x040(%%r14),%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t	vmulpd	     (%%r14),%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)					\n\t	vmulpd	0x040(%%r14),%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t	 vfmadd231pd	0x020(%%r14),%%ymm11,%%ymm12\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)					\n\t	 vfmadd231pd	0x060(%%r14),%%ymm14,%%ymm15\n\t"\
		"addq	%%rdi,%%rax	/* a[j+p8 ] */				\n\t	vfnmadd231pd	0x020(%%r14),%%ymm10,%%ymm13\n\t"\
		"addq	%%rdi,%%rbx	/* a[j+p12] */				\n\t	vfnmadd231pd	0x060(%%r14),%%ymm8 ,%%ymm9 \n\t"\
		"addq	$0x080,%%rsi	/* c4 */				\n\t	vmovaps	%%ymm12,     (%%r10)		\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	vmovaps	%%ymm15,     (%%r11)		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t	vmovaps	%%ymm13,0x020(%%r10)		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	%%ymm9 ,0x020(%%r11)		\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	addq	%%rdi,%%r10	/* a[j+p10] */	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t	addq	%%rdi,%%r11	/* a[j+p14] */	\n\t"\
		"vmovaps	%%ymm0,%%ymm1						\n\t	addq	$0x080,%%r14	/* c2 */	\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmovaps	%%ymm6,%%ymm7						\n\t	vmovaps	0x120(%%r12),%%ymm8 		\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4				\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm0,%%ymm0				\n\t	vmovaps	0x100(%%r12),%%ymm14	\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5				\n\t	vmovaps	%%ymm12,%%ymm10			\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm6,%%ymm6				\n\t	vmovaps	%%ymm8 ,%%ymm9 			\n\t"\
		" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4		\n\t	vmovaps	%%ymm13,%%ymm11			\n\t"\
		" vfmadd231pd	0x060(%%rsi),%%ymm7,%%ymm0		\n\t	vmovaps	%%ymm14,%%ymm15			\n\t"\
		"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5		\n\t	vmulpd	     (%%r14),%%ymm12,%%ymm12		\n\t"\
		"vfnmadd231pd	0x060(%%rsi),%%ymm1,%%ymm6		\n\t	vmulpd	0x040(%%r14),%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t	vmulpd	     (%%r14),%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)					\n\t	vmulpd	0x040(%%r14),%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t	 vfmadd231pd	0x020(%%r14),%%ymm11,%%ymm12	\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)					\n\t	 vfmadd231pd	0x060(%%r14),%%ymm15,%%ymm8 	\n\t"\
		"/*...Block 2: 34 MOVapd, 36 ADD/SUB, 26 MUL */	\n\t	vfnmadd231pd	0x020(%%r14),%%ymm10,%%ymm13	\n\t"\
		"subq	%%rdi,%%rax	/* a[j+p0 ] */				\n\t	vfnmadd231pd	0x060(%%r14),%%ymm9 ,%%ymm14	\n\t"\
		"movslq	%[__p1],%%r14							\n\t	vmovaps	%%ymm12,     (%%r10)		\n\t"\
		"subq	%%rdi,%%rbx	/* a[j+p4 ] */				\n\t	vmovaps	%%ymm8 ,     (%%r11)		\n\t"\
		"shlq	$3,%%r14								\n\t	vmovaps	%%ymm13,0x020(%%r10)		\n\t"\
		"addq	%%r14,%%rax	/* a[j+p1 ] */				\n\t	vmovaps	%%ymm14,0x020(%%r11)		\n\t"\
		"addq	%%r14,%%rbx	/* a[j+p5 ] */				\n\t/*...Block 4: 34 MOVapd, 36 ADD/SUBpd, 26 MULpd */\n\t"\
		"addq	$0x040,%%rcx		/* r3 , from r1 */	\n\t	subq	%%rdi,%%r10	/* a[j+p2 ] */	\n\t"\
		"leaq	0x100(%%rcx),%%rdx	/* r11, from r3 */	\n\t	subq	%%rdi,%%r11	/* a[j+p6 ] */	\n\t"\
		"subq	$0x0c0,%%rsi		/* cc0, from c4 */	\n\t	addq	%%r14,%%r10	/* a[j+p3 ] */	\n\t"\
		"vmovaps	0x200(%%rcx),%%ymm4					\n\t	addq	%%r14,%%r11	/* a[j+p7 ] */	\n\t"\
		"vmovaps	0x220(%%rcx),%%ymm5					\n\t	addq	$0x240,%%r12	/* r7 , from r5  */\n\t"\
		"vmovaps	     (%%rsi),%%ymm2					\n\t"	/* r15 in %%r12 + 0x100 */\
		"vmovaps	0x020(%%rsi),%%ymm3					\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmovaps	     (%%rsi),%%ymm11	\n\t"\
		"vmovaps	0x200(%%rdx),%%ymm0					\n\t	vmovaps	0x020(%%rsi),%%ymm10	\n\t"\
		"vmovaps	0x220(%%rdx),%%ymm1					\n\t	vmovaps	%%ymm12,%%ymm14			\n\t"\
		"vfnmadd231pd	%%ymm3,%%ymm6,%%ymm5			\n\t	vmovaps	%%ymm13,%%ymm15			\n\t"\
		" vfmadd231pd	%%ymm3,%%ymm7,%%ymm4			\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t	vmovaps	0x100(%%r12),%%ymm8 		\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0					\n\t	vmovaps	0x120(%%r12),%%ymm9 		\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1					\n\t	vfnmadd231pd	%%ymm11,%%ymm14,%%ymm13	\n\t"\
		"vfnmadd231pd	%%ymm2,%%ymm6,%%ymm1			\n\t	 vfmadd231pd	%%ymm11,%%ymm15,%%ymm12	\n\t"\
		" vfmadd231pd	%%ymm2,%%ymm7,%%ymm0			\n\t	vmovaps	%%ymm8 ,%%ymm14					\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	%%ymm9 ,%%ymm15			\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmulpd	%%ymm11,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4					\n\t	vmulpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5					\n\t	vfnmadd231pd	%%ymm10,%%ymm14,%%ymm9 	\n\t"\
		"vsubpd	%%ymm0,%%ymm6,%%ymm6					\n\t	 vfmadd231pd	%%ymm10,%%ymm15,%%ymm8 	\n\t"\
		"vsubpd	%%ymm1,%%ymm7,%%ymm7					\n\t	vmovaps	%%ymm13,%%ymm15			\n\t"\
		"vmovaps	     (%%rdx),%%ymm2					\n\t	vmovaps	%%ymm12,%%ymm14			\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3					\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	-0x460(%%r8),%%ymm1					\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"vmovaps	%%ymm3,%%ymm0						\n\t	vsubpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vsubpd	%%ymm2,%%ymm3,%%ymm3					\n\t	vsubpd	%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2					\n\t	subq	$0x200,%%r12			\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2					\n\t"\
		"vmulpd	%%ymm1,%%ymm3,%%ymm3					\n\t	vmovaps	0x100(%%r12),%%ymm8 		\n\t"\
		"vmovaps	     (%%rcx),%%ymm0					\n\t	vmovaps	0x120(%%r12),%%ymm9 		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t	vmovaps	-0x460(%%r8),%%ymm11	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t	vmovaps	%%ymm8 ,%%ymm10			\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm9 ,%%ymm8 ,%%ymm8 			\n\t"\
		" vfmadd132pd	(%%r8),%%ymm0,%%ymm2			\n\t	vaddpd	%%ymm10,%%ymm9 ,%%ymm9 			\n\t"\
		" vfmadd132pd	(%%r8),%%ymm1,%%ymm3			\n\t	vmulpd	%%ymm11,%%ymm8 ,%%ymm8 			\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r13)\n\t"/* ...+p12 */\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c1) */	\n\t	vmulpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"addq	$0x240,%%rsi	/* c1, from cc0 */		\n\t	vmovaps	     (%%r12),%%ymm10	\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t	vmovaps	0x020(%%r12),%%ymm11	\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t	vsubpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1					\n\t	 vfmadd132pd	(%%r8),%%ymm10,%%ymm8 	\n\t"\
		" vfmadd132pd	(%%r8),%%ymm2,%%ymm4			\n\t	 vfmadd132pd	(%%r8),%%ymm11,%%ymm9 	\n\t"\
		" vfmadd132pd	(%%r8),%%ymm0,%%ymm7			\n\t/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c3) */\n\t"\
		" vfmadd132pd	(%%r8),%%ymm3,%%ymm5			\n\t	leaq	0x100(%%rsi),%%r14	/* c3, from c1 */\n\t"\
		" vfmadd132pd	(%%r8),%%ymm1,%%ymm6			\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)					\n\t	vsubpd	%%ymm15,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)					\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11			\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)					\n\t	vsubpd	%%ymm14,%%ymm9 ,%%ymm9 			\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)					\n\t	 vfmadd132pd	(%%r8),%%ymm10,%%ymm12	\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t	 vfmadd132pd	(%%r8),%%ymm8 ,%%ymm15	\n\t"\
		"vmovaps	%%ymm7,%%ymm0						\n\t	 vfmadd132pd	(%%r8),%%ymm11,%%ymm13	\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t	 vfmadd132pd	(%%r8),%%ymm9 ,%%ymm14	\n\t"\
		"vmovaps	%%ymm1,%%ymm6						\n\t	vmovaps	%%ymm10,     (%%r12)	\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4				\n\t	vmovaps	%%ymm8 ,0x120(%%r12)	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm7,%%ymm7				\n\t	vmovaps	%%ymm11,0x020(%%r12)	\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5				\n\t	vmovaps	%%ymm14,0x100(%%r12)	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm1,%%ymm1				\n\t	vmovaps	%%ymm12,%%ymm10			\n\t"\
		"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5		\n\t	vmovaps	%%ymm15,%%ymm8 			\n\t"\
		"vfnmadd231pd	0x060(%%rsi),%%ymm0,%%ymm1		\n\t	vmovaps	%%ymm13,%%ymm11			\n\t"\
		" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4		\n\t	vmovaps	%%ymm9 ,%%ymm14			\n\t"\
		" vfmadd231pd	0x060(%%rsi),%%ymm6,%%ymm7		\n\t	vmulpd	     (%%r14),%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t	vmulpd	0x040(%%r14),%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)					\n\t	vmulpd	     (%%r14),%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t	vmulpd	0x040(%%r14),%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)					\n\t	vfnmadd231pd	0x020(%%r14),%%ymm10,%%ymm13\n\t"\
		"addq	%%rdi,%%rax								\n\t	vfnmadd231pd	0x060(%%r14),%%ymm8 ,%%ymm9 \n\t"\
		"addq	%%rdi,%%rbx								\n\t	 vfmadd231pd	0x020(%%r14),%%ymm11,%%ymm12\n\t"\
		"addq	$0x080,%%rsi	/* c2 */				\n\t	 vfmadd231pd	0x060(%%r14),%%ymm14,%%ymm15\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	vmovaps	%%ymm13,0x020(%%r10)		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t	vmovaps	%%ymm9 ,0x020(%%r11)		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	%%ymm12,     (%%r10)		\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	vmovaps	%%ymm15,     (%%r11)		\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t	addq	%%rdi,%%r10				\n\t"\
		"vmovaps	%%ymm0,%%ymm1						\n\t	addq	%%rdi,%%r11				\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t	addq	$0x080,%%r14	/* c2 */	\n\t"\
		"vmovaps	%%ymm6,%%ymm7						\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4				\n\t	vmovaps	0x120(%%r12),%%ymm8 		\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm0,%%ymm0				\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5				\n\t	vmovaps	0x100(%%r12),%%ymm14	\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm6,%%ymm6				\n\t	vmovaps	%%ymm12,%%ymm10			\n\t"\
		" vfmadd231pd	0x020(%%rsi),%%ymm3,%%ymm4		\n\t	vmovaps	%%ymm8 ,%%ymm9 			\n\t"\
		" vfmadd231pd	0x060(%%rsi),%%ymm7,%%ymm0		\n\t	vmovaps	%%ymm13,%%ymm11			\n\t"\
		"vfnmadd231pd	0x020(%%rsi),%%ymm2,%%ymm5		\n\t	vmovaps	%%ymm14,%%ymm15			\n\t"\
		"vfnmadd231pd	0x060(%%rsi),%%ymm1,%%ymm6		\n\t	vmulpd	     (%%r14),%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t	vmulpd	0x040(%%r14),%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)					\n\t	vmulpd	     (%%r14),%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t	vmulpd	0x040(%%r14),%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)					\n\t	 vfmadd231pd	0x020(%%r14),%%ymm11,%%ymm12\n\t"\
		"												\n\t	 vfmadd231pd	0x060(%%r14),%%ymm15,%%ymm8 \n\t"\
		"												\n\t	vfnmadd231pd	0x020(%%r14),%%ymm10,%%ymm13\n\t"\
		"												\n\t	vfnmadd231pd	0x060(%%r14),%%ymm9 ,%%ymm14\n\t"\
		"												\n\t	vmovaps	%%ymm12,     (%%r10)		\n\t"\
		"												\n\t	vmovaps	%%ymm8 ,     (%%r11)		\n\t"\
		"												\n\t	vmovaps	%%ymm13,0x020(%%r10)		\n\t"\
		"												\n\t	vmovaps	%%ymm14,0x020(%%r11)		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__r1] "m" (Xr1)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers

	/* Do 4 prefetches from main-array address offsets ...+p0,4,8,c in this macro, and want to spread them out
	roughly evenly, which works out to 1 every 100 lines or so. Here is the i/o register usage summary for the
	four distinct SSE2_RADIX_04_DIF_3TWIDDLE_X1,2 - corresponding sub-blocks [1-4] and the actual address-offset
	we end up using in both a 16-prefetch (as in SSE2_RADIX16_DIF_TWIDDLE(), but that seems slower in _V2 here)
	and 4-prefetch paradigm:
									prefetch offsets:
								 possible:				actual:
												16-fetch	4-fetch
	[1]	[a-d]x: 0,4,8,c	rsi: 2	0,2,4,6,8,a,c,e	0,4,8,c		0,8 [one at start of block, other at end]
	[2]	[a-d]x: 1,5,9,d	rsi: 2	1,3,5,7,9,b,d,f	1,5,9,d
	[3]	r10-13: 0,2,1,3	r9 : 4	0-7				2,3,6,7		4
	[4]	r10-13: 8,a,9,b	r9 : 4	8-f				a,b,e,f		C
	*/
	#define SSE2_RADIX16_DIF_TWIDDLE_V2(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr0,Xtwo,Xcc0,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/*
		k1 = p2*8;	// Set k1 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k2 = 0x80;	// Set k2 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		// Pass1-DFTs all use same twiddle-triplet, so lit-byte address offset between ptrs to 1st,2nd DFT's twiddles = 0
		i0 = (vec_dbl *)add0; i1 = (vec_dbl *)(add0+p4); i2 = (vec_dbl *)(add0+p8); i3 = (vec_dbl *)(add0+p12);
		o0 = r1; o1 = r1+2; o2 = r1+4; o3 = r1+6;
		c_tmp = cc0;	// c8,4,C
		SSE2_RADIX_04_DIF_3TWIDDLE_X1(i0,i1,i2,i3,k1, two,c_tmp, o0,o1,o2,o3,k2)
	*/\
		"movq	%[__cc0],%%rsi 			\n\t	movslq	%[__p2],%%rdi		\n\t"\
		"movq	%[__add0],%%rax			\n\t	movq	%[__r0],%%r10		\n\t"\
		"movslq	%[__p4],%%rbx			\n\t	leaq	(%%rax,%%rbx,8),%%rbx	\n\t"\
		"movslq	%[__p8],%%rcx			\n\t	leaq	(%%rax,%%rcx,8),%%rcx	\n\t"\
		"movslq	%[__p12],%%rdx			\n\t	leaq	(%%rax,%%rdx,8),%%rdx	\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%rax)\n\t"/* [base-address + data-fetch-ahead index] + p0 */\
		/* Do the p0,p2 combo: */\
		"vmovaps	     (%%rcx),%%ymm4	\n\t	shlq	$3,%%rdi			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rsi),%%ymm10\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm11\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t	vmovaps	     (%%rcx,%%rdi),%%ymm12	\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	0x020(%%rcx,%%rdi),%%ymm13	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmulpd	%%ymm10,%%ymm4,%%ymm4	\n\t	vmovaps	%%ymm12,%%ymm14				\n\t"\
		"vmulpd	%%ymm10,%%ymm5,%%ymm5	\n\t	vmovaps	%%ymm13,%%ymm15				\n\t"\
		"vmulpd	%%ymm11,%%ymm6,%%ymm6	\n\t	vmovaps	     (%%rax,%%rdi),%%ymm8 	\n\t	vmovaps	%%ymm0,%%ymm2		\n\t"\
		"vmulpd	%%ymm11,%%ymm7,%%ymm7	\n\t	vmovaps	0x020(%%rax,%%rdi),%%ymm9 	\n\t	vmovaps	%%ymm1,%%ymm3		\n\t"\
		"vaddpd	%%ymm6 ,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm7 ,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm4 ,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm11,%%ymm14,%%ymm14		\n\t	vmovaps	%%ymm8 ,%%ymm10		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm11,%%ymm15,%%ymm15		\n\t	vmovaps	%%ymm9 ,%%ymm11		\n\t"\
		"vsubpd	%%ymm4 ,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm3,%%ymm3	\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm0,     (%%r10)	\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"/* Spill 1: free up ymm0,1 */\
		"vmovaps	%%ymm1,0x020(%%r10)	\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"/* Do the p1,3 combo: */\
		"vmovaps	0x080(%%rsi),%%ymm0	\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm1	\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t	vmovaps	%%ymm8 ,0x100(%%r10)		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t	vmovaps	%%ymm9 ,0x120(%%r10)		\n\t"\
		"vmovaps	%%ymm6,%%ymm4		\n\t"\
		"vmovaps	%%ymm7,%%ymm5		\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vmovaps	     (%%rdx,%%rdi),%%ymm14	\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t	vmovaps	0x020(%%rdx,%%rdi),%%ymm15	\n\t"\
		"vmulpd	%%ymm1,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm14,%%ymm12				\n\t"\
		"vmulpd	%%ymm1,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm0 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm0 ,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm5,0x0a0(%%r10)	\n\t	vmulpd	%%ymm1 ,%%ymm14,%%ymm14		\n\t"/* Spill 2 */\
		"vmovaps	%%ymm4,0x080(%%r10)	\n\t	vmulpd	%%ymm1 ,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm8	\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm9	\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	     (%%rbx),%%ymm6	\n\t	vmovaps	%%ymm13,0x1a0(%%r10)		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t	vmovaps	%%ymm12,0x180(%%r10)		\n\t"\
		"vmovaps	%%ymm6,%%ymm4		\n\t"\
		"vmovaps	%%ymm7,%%ymm5		\n\t"\
		"vmulpd	%%ymm8,%%ymm4,%%ymm4	\n\t	vmovaps	     (%%rbx,%%rdi),%%ymm14	\n\t"\
		"vmulpd	%%ymm8,%%ymm5,%%ymm5	\n\t	vmovaps	0x020(%%rbx,%%rdi),%%ymm15	\n\t"\
		"vmulpd	%%ymm9,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm14,%%ymm12				\n\t"\
		"vmulpd	%%ymm9,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"\
		"vmovaps	0x080(%%r10),%%ymm0	\n\t	vmulpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"/* Restore 2 */\
		"vmovaps	0x0a0(%%r10),%%ymm1	\n\t	vmulpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t	 movq	%[__two],%%rsi	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	0x180(%%r10),%%ymm8 		\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t	vmovaps	0x1a0(%%r10),%%ymm9 		\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5	\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm13,%%ymm15				\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm12,%%ymm14				\n\t"\
		/* Finish radix-4 butterfly and store results: */\
		"vmovaps	     (%%r10),%%ymm0	\n\t	vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"/* Restore 1 */\
		"vmovaps	0x020(%%r10),%%ymm1	\n\t	vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0	\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1	\n\t	vmovaps	0x100(%%r10),%%ymm8 		\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmovaps	0x120(%%r10),%%ymm9 		\n\t"\
		"vmovaps	%%ymm0,0x080(%%r10)	\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"	/* 2.0, shared by both columns ... moved +- until found best cycle count: */\
		"vmovaps	%%ymm2,0x040(%%r10)	\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10		\n\t	vmovaps	(%%rsi),%%ymm0	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%r10)	\n\t	vsubpd	%%ymm15,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%r10)	\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	%%ymm0,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm8 ,0x180(%%r10)		\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t	vmovaps	%%ymm10,0x140(%%r10)		\n\t"\
		"vmulpd	%%ymm0,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm9 ,0x1a0(%%r10)		\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vmovaps	%%ymm11,0x1e0(%%r10)		\n\t"\
	"vaddpd	0x080(%%r10),%%ymm6,%%ymm6	\n\t	vmulpd	%%ymm0 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	 %%ymm2 ,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm0 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	 %%ymm1 ,%%ymm7,%%ymm7	\n\t	vmulpd	%%ymm0 ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	 %%ymm3 ,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm0 ,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm6,     (%%r10)	\n\t	vaddpd 0x180(%%r10),%%ymm14,%%ymm14	\n\t"/* don't need reload-from-mem of ymm8/0x180(%%r10) as we do in lcol, but 1 cycle faster with it. [!?] */\
		"vmovaps	%%ymm5,0x0c0(%%r10)	\n\t	vaddpd  %%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm7,0x020(%%r10)	\n\t	vaddpd  %%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,0x060(%%r10)	\n\t	vaddpd  %%ymm11,%%ymm12,%%ymm12		\n\t"\
		"										vmovaps	%%ymm14,0x100(%%r10)		\n\t"\
		"										vmovaps	%%ymm13,0x1c0(%%r10)		\n\t"\
		"										vmovaps	%%ymm15,0x120(%%r10)		\n\t"\
		"										vmovaps	%%ymm12,0x160(%%r10)		\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%rcx)\n\t"/* ... + p8 */\
	/*
		i0 = (vec_dbl *)(add0+p1); i1 = (vec_dbl *)(add0+p4+p1); i2 = (vec_dbl *)(add0+p8+p1); i3 = (vec_dbl *)(add0+p12+p1);
		o0 += 16; o1 += 16; o2 += 16; o3 += 16;
		SSE2_RADIX_04_DIF_3TWIDDLE_X1(i0,i1,i2,i3,k1, two,c_tmp, o0,o1,o2,o3,k2)
	*/\
		"movslq	%[__p1],%%r9			\n\t	addq	$0x200,%%r10		\n\t"\
		"leaq	(%%rax,%%r9,8),%%rax	\n\t	movq	%[__cc0],%%rsi 		\n\t"/* repoint rsi from two -> cc0 */\
		"leaq	(%%rbx,%%r9,8),%%rbx	\n\t"\
		"leaq	(%%rcx,%%r9,8),%%rcx	\n\t"\
		"leaq	(%%rdx,%%r9,8),%%rdx	\n\t"\
		/* Do the p0,p2 combo: */\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rsi),%%ymm10\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm11\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t	vmovaps	     (%%rcx,%%rdi),%%ymm12	\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	0x020(%%rcx,%%rdi),%%ymm13	\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t"\
		"vmulpd	%%ymm10,%%ymm4,%%ymm4	\n\t	vmovaps	%%ymm12,%%ymm14				\n\t"\
		"vmulpd	%%ymm10,%%ymm5,%%ymm5	\n\t	vmovaps	%%ymm13,%%ymm15				\n\t"\
		"vmulpd	%%ymm11,%%ymm6,%%ymm6	\n\t	vmovaps	     (%%rax,%%rdi),%%ymm8 	\n\t	vmovaps	%%ymm0,%%ymm2		\n\t"\
		"vmulpd	%%ymm11,%%ymm7,%%ymm7	\n\t	vmovaps	0x020(%%rax,%%rdi),%%ymm9 	\n\t	vmovaps	%%ymm1,%%ymm3		\n\t"\
		"vaddpd	%%ymm6 ,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm7 ,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm4 ,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm11,%%ymm14,%%ymm14		\n\t	vmovaps	%%ymm8 ,%%ymm10		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm11,%%ymm15,%%ymm15		\n\t	vmovaps	%%ymm9 ,%%ymm11		\n\t"\
		"vsubpd	%%ymm4 ,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm3,%%ymm3	\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm0,     (%%r10)	\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"/* Spill 1: free up ymm0,1 */\
		"vmovaps	%%ymm1,0x020(%%r10)	\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"/* Do the p1,3 combo: */\
		"vmovaps	0x080(%%rsi),%%ymm0	\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm1	\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t	vmovaps	%%ymm8 ,0x100(%%r10)		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t	vmovaps	%%ymm9 ,0x120(%%r10)		\n\t"\
		"vmovaps	%%ymm6,%%ymm4		\n\t"\
		"vmovaps	%%ymm7,%%ymm5		\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vmovaps	     (%%rdx,%%rdi),%%ymm14	\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t	vmovaps	0x020(%%rdx,%%rdi),%%ymm15	\n\t"\
		"vmulpd	%%ymm1,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm14,%%ymm12				\n\t"\
		"vmulpd	%%ymm1,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm0 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm0 ,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm5,0x0a0(%%r10)	\n\t	vmulpd	%%ymm1 ,%%ymm14,%%ymm14		\n\t"/* Spill 2 */\
		"vmovaps	%%ymm4,0x080(%%r10)	\n\t	vmulpd	%%ymm1 ,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm8	\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm9	\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	     (%%rbx),%%ymm6	\n\t	vmovaps	%%ymm13,0x1a0(%%r10)		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t	vmovaps	%%ymm12,0x180(%%r10)		\n\t"\
		"vmovaps	%%ymm6,%%ymm4		\n\t"\
		"vmovaps	%%ymm7,%%ymm5		\n\t"\
		"vmulpd	%%ymm8,%%ymm4,%%ymm4	\n\t	vmovaps	     (%%rbx,%%rdi),%%ymm14	\n\t"\
		"vmulpd	%%ymm8,%%ymm5,%%ymm5	\n\t	vmovaps	0x020(%%rbx,%%rdi),%%ymm15	\n\t"\
		"vmulpd	%%ymm9,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm14,%%ymm12				\n\t"\
		"vmulpd	%%ymm9,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"\
		"vmovaps	0x080(%%r10),%%ymm0	\n\t	vmulpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"/* Restore 2 */\
		"vmovaps	0x0a0(%%r10),%%ymm1	\n\t	vmulpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t	 movq	%[__two],%%rsi	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	0x180(%%r10),%%ymm8 		\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t	vmovaps	0x1a0(%%r10),%%ymm9 		\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5	\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm13,%%ymm15				\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm12,%%ymm14				\n\t"\
		/* Finish radix-4 butterfly and store results: */\
		"vmovaps	     (%%r10),%%ymm0	\n\t	vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"/* Restore 1 */\
		"vmovaps	0x020(%%r10),%%ymm1	\n\t	vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0	\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1	\n\t	vmovaps	0x100(%%r10),%%ymm8 		\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmovaps	0x120(%%r10),%%ymm9 		\n\t"\
		"vmovaps	%%ymm0,0x080(%%r10)	\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"	/* 2.0, shared by both columns ... moved +- until found best cycle count: */\
		"vmovaps	%%ymm2,0x040(%%r10)	\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10		\n\t	vmovaps	(%%rsi),%%ymm0	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%r10)	\n\t	vsubpd	%%ymm15,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%r10)	\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	%%ymm0,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm8 ,0x180(%%r10)		\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t	vmovaps	%%ymm10,0x140(%%r10)		\n\t"\
		"vmulpd	%%ymm0,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm9 ,0x1a0(%%r10)		\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vmovaps	%%ymm11,0x1e0(%%r10)		\n\t"\
	"vaddpd	0x080(%%r10),%%ymm6,%%ymm6	\n\t	vmulpd	%%ymm0 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	 %%ymm2 ,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm0 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	 %%ymm1 ,%%ymm7,%%ymm7	\n\t	vmulpd	%%ymm0 ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	 %%ymm3 ,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm0 ,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm6,     (%%r10)	\n\t	vaddpd 0x180(%%r10),%%ymm14,%%ymm14	\n\t"/* don't need reload-from-mem of ymm8/0x180(%%r10) as we do in lcol, but 1 cycle faster with it. [!?] */\
		"vmovaps	%%ymm5,0x0c0(%%r10)	\n\t	vaddpd  %%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm7,0x020(%%r10)	\n\t	vaddpd  %%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,0x060(%%r10)	\n\t	vaddpd  %%ymm11,%%ymm12,%%ymm12		\n\t"\
		"										vmovaps	%%ymm14,0x100(%%r10)		\n\t"\
		"										vmovaps	%%ymm13,0x1c0(%%r10)		\n\t"\
		"										vmovaps	%%ymm15,0x120(%%r10)		\n\t"\
		"										vmovaps	%%ymm12,0x160(%%r10)		\n\t"\
	/*
		// Pass 2:
		k1 = 0x080;	// Set k1 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k2 = p4*8;	// Set k2 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		// 0x0c0 = literal-bytewise address offset between the pointers to the first and second DFT's twiddles = 0:
		i0 = r1; i1 = r1+16; i2 = r1+8; i3 = r1+24;
		o0 = (vec_dbl *)add0; o1 = (vec_dbl *)(add0+p2); o2 = (vec_dbl *)(add0+p1); o3 = (vec_dbl *)(add0+p3);
		c_tmp += 6;	// c2,A,6
		SSE2_RADIX_04_DIF_3TWIDDLE_X2(i0,i1,i2,i3,k1, two,c_tmp,0x0c0, o0,o1,o2,o3,k2)
	*/\
		/* Take advantage of fact that i1 = r0 + 16 already in %%r10 and o2 = add0+p1 in %%rax : */\
		/* i-offset = 0x080, so inline with r[a-d]x base-address offsets rather than adding via (r[a-d]x,rdi) */\
		"movq	%%r10,%%rbx				\n\t	movq	%%rax,%%r12			\n\t"/* rbx:i1 = r0+16; r12:o2 = add0+p1 */\
		"movq	%[__cc0],%%rsi 			\n\t	addq	$0x0c0,%%rsi 		\n\t"/* repoint rsi from two -> cc0, then increment for Pass 2 */\
		"leaq	 0x0c0(%%rsi),%%r8		\n\t	movslq	%[__p2],%%r11		\n\t"\
		"leaq	-0x200(%%rbx),%%rax		\n\t	movslq	%[__p4],%%r9		\n\t"\
		"leaq	-0x100(%%rbx),%%rcx		\n\t	movq	%[__add0]	,%%r10		\n\t"/* o0 */\
		"leaq	 0x100(%%rbx),%%rdx		\n\t	leaq	(%%r12,%%r11,8),%%r13	\n\t"/* o3 = (add0+p1)+p2 */\
		"vmovaps	     (%%rax),%%ymm0	\n\t	leaq	(%%r10,%%r11,8),%%r11	\n\t"/* o1 = (add0+p2) */\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t	shlq	$3,%%r9			\n\t"/* p4 in byte-offset form */\
		/* Do the p0,p2 combo: */\
		"vmovaps	     (%%rcx),%%ymm4	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rsi),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3	\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r9)\n\t"/* ... + p4 */\
		"vmovaps	%%ymm4,%%ymm6		\n\t	vmovaps	0x080(%%rcx),%%ymm12		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	0x0a0(%%rcx),%%ymm13		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t	vmovaps	     (%%r8),%%ymm10			\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t	vmovaps	0x020(%%r8),%%ymm11			\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4	\n\t	vmovaps	%%ymm12,%%ymm14				\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5	\n\t	vmovaps	%%ymm13,%%ymm15				\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	\n\t	vmovaps	0x080(%%rax),%%ymm8 		\n\t	vmovaps	%%ymm0,%%ymm2		\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	\n\t	vmovaps	0x0a0(%%rax),%%ymm9 		\n\t	vmovaps	%%ymm1,%%ymm3		\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm11,%%ymm14,%%ymm14		\n\t	vmovaps	%%ymm8 ,%%ymm10		\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm11,%%ymm15,%%ymm15		\n\t	vmovaps	%%ymm9 ,%%ymm11		\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3	\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm0,     (%%r10)	\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"/* Spill 1: free up ymm0,1 */\
		"vmovaps	%%ymm1,0x020(%%r10)	\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"/* Do the p1,3 combo: */\
		"vmovaps	0x080(%%rsi),%%ymm0	\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm1	\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t	vmovaps	%%ymm8 ,     (%%r10,%%r9)	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t	vmovaps	%%ymm9 ,0x020(%%r10,%%r9)	\n\t"\
		"vmovaps	%%ymm6,%%ymm4		\n\t	vmovaps	0x080(%%r8),%%ymm8	 		\n\t"\
		"vmovaps	%%ymm7,%%ymm5		\n\t	vmovaps	0x0a0(%%r8),%%ymm9 			\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vmovaps	0x080(%%rdx),%%ymm14		\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t	vmovaps	0x0a0(%%rdx),%%ymm15		\n\t"\
		"vmulpd	%%ymm1,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm14,%%ymm12				\n\t"\
		"vmulpd	%%ymm1,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm5,0x020(%%r12)	\n\t	vmulpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t"/* Spill 2 */\
		"vmovaps	%%ymm4,     (%%r12)	\n\t	vmulpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm0	\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm1	\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	     (%%rbx),%%ymm6	\n\t	vmovaps	%%ymm13,0x020(%%r12,%%r9)	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t	vmovaps	%%ymm12,     (%%r12,%%r9)	\n\t"\
		"vmovaps	%%ymm6,%%ymm4		\n\t	vmovaps	0x040(%%r8),%%ymm8	 		\n\t"\
		"vmovaps	%%ymm7,%%ymm5		\n\t	vmovaps	0x060(%%r8),%%ymm9 			\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vmovaps	0x080(%%rbx),%%ymm14		\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t	vmovaps	0x0a0(%%rbx),%%ymm15		\n\t"\
		"vmulpd	%%ymm1,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm14,%%ymm12				\n\t"\
		"vmulpd	%%ymm1,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"\
		"vmovaps	     (%%r12),%%ymm0	\n\t	vmulpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"/* Restore 2 */\
		"vmovaps	0x020(%%r12),%%ymm1	\n\t	vmulpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t	 movq	%[__two],%%rsi	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	     (%%r12,%%r9),%%ymm8	\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t	vmovaps	0x020(%%r12,%%r9),%%ymm9	\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5	\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm13,%%ymm15				\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm12,%%ymm14				\n\t"\
		/* Finish radix-4 butterfly and store results: */\
		"vmovaps	     (%%r10),%%ymm0	\n\t	vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"/* Restore 1 */\
		"vmovaps	0x020(%%r10),%%ymm1	\n\t	vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0	\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1	\n\t	vmovaps	     (%%r10,%%r9),%%ymm8	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmovaps	0x020(%%r10,%%r9),%%ymm9	\n\t"\
		"vmovaps	%%ymm0,     (%%r12)	\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"	/* 2.0, shared by both columns ... moved +- until found best cycle count: */\
		"vmovaps	%%ymm2,     (%%r11)	\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10		\n\t	vmovaps	(%%rsi),%%ymm0	\n\t"\
		"vmovaps	%%ymm1,0x020(%%r12)	\n\t	vsubpd	%%ymm15,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm3,0x020(%%r13)	\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	%%ymm0,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm8 ,     (%%r12,%%r9)	\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t	vmovaps	%%ymm10,     (%%r11,%%r9)	\n\t"\
		"vmulpd	%%ymm0,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm9 ,0x020(%%r12,%%r9)	\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vmovaps	%%ymm11,0x020(%%r13,%%r9)	\n\t"\
		"vaddpd	(%%r12),%%ymm6,%%ymm6	\n\t	vmulpd	%%ymm0 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm2 ,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm0 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm1 ,%%ymm7,%%ymm7	\n\t	vmulpd	%%ymm0 ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm3 ,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm0 ,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm6,     (%%r10)	\n\t	vaddpd (%%r12,%%r9),%%ymm14,%%ymm14	\n\t"/* don't need reload-from-mem of ymm8/0x180(%%r12,%%r9) as we do in lcol, but 1 cycle faster with it. [!?] */\
		"vmovaps	%%ymm5,     (%%r13)	\n\t	vaddpd      %%ymm10,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm7,0x020(%%r10)	\n\t	vaddpd      %%ymm9 ,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm4,0x020(%%r11)	\n\t	vaddpd      %%ymm11,%%ymm12,%%ymm12	\n\t"\
		"										vmovaps	%%ymm14,     (%%r10,%%r9)	\n\t"\
		"										vmovaps	%%ymm13,     (%%r13,%%r9)	\n\t"\
		"										vmovaps	%%ymm15,0x020(%%r10,%%r9)	\n\t"\
		"										vmovaps	%%ymm12,0x020(%%r11,%%r9)	\n\t"\
	/*
		i0 += 2; i1 += 2; i2 += 2; i3 += 2;
		o0 = (vec_dbl *)(add0+p8); o1 = (vec_dbl *)(add0+p8+p2); o2 = (vec_dbl *)(add0+p8+p1); o3 = (vec_dbl *)(add0+p8+p3);
		c_tmp += 12;	// c5,D,3
		SSE2_RADIX_04_DIF_3TWIDDLE_X2(i0,i1,i2,i3,k1, two,c_tmp,0x0c0, o0,o1,o2,o3,k2)
	*/\
		/* Take advantage of fact that i1 = r0 + 16 already in %%r10 and o2 = add0+p1 in %%rax : */\
		"movq	%[__cc0],%%rsi 			\n\t	addq	$0x240,%%rsi 		\n\t"/* repoint rsi from two -> cc0, then increment for Pass, set 2 */\
		"addq	$0x040,%%rax			\n\t	movslq	%[__p8],%%r8		\n\t"\
		"addq	$0x040,%%rbx			\n\t"	/* r9 still has p4<<3 */\
		"addq	$0x040,%%rcx			\n\t	leaq	(%%r10,%%r8,8),%%r10	\n\t"/* o0 =  add0     +p8 */\
		"addq	$0x040,%%rdx			\n\t	leaq	(%%r11,%%r8,8),%%r11	\n\t"/* o1 = (add0+p2)+p8 */\
		"vmovaps	     (%%rax),%%ymm0	\n\t	leaq	(%%r12,%%r8,8),%%r12	\n\t"/* o2 = (add0+p1)+p8 */\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t	leaq	(%%r13,%%r8,8),%%r13	\n\t"/* o3 = (add0+p3)+p8 */\
		/* Do the p0,p2 combo: */\
		"vmovaps	     (%%rcx),%%ymm4	\n\t	leaq	 0x0c0(%%rsi),%%r8	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
		"vmovaps	     (%%rsi),%%ymm2	\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3	\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r9)\n\t"/* ... + pC */\
		"vmovaps	%%ymm4,%%ymm6		\n\t	vmovaps	0x080(%%rcx),%%ymm12		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	0x0a0(%%rcx),%%ymm13		\n\t"\
		"vmovaps	     (%%rax),%%ymm0	\n\t	vmovaps	     (%%r8),%%ymm10			\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t	vmovaps	0x020(%%r8),%%ymm11			\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4	\n\t	vmovaps	%%ymm12,%%ymm14				\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5	\n\t	vmovaps	%%ymm13,%%ymm15				\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6	\n\t	vmovaps	0x080(%%rax),%%ymm8 		\n\t	vmovaps	%%ymm0,%%ymm2		\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7	\n\t	vmovaps	0x0a0(%%rax),%%ymm9 		\n\t	vmovaps	%%ymm1,%%ymm3		\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm11,%%ymm14,%%ymm14		\n\t	vmovaps	%%ymm8 ,%%ymm10		\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm11,%%ymm15,%%ymm15		\n\t	vmovaps	%%ymm9 ,%%ymm11		\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3	\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm0,     (%%r10)	\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"/* Spill 1: free up ymm0,1 */\
		"vmovaps	%%ymm1,0x020(%%r10)	\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"/* Do the p1,3 combo: */\
		"vmovaps	0x080(%%rsi),%%ymm0	\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm1	\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	     (%%rdx),%%ymm6	\n\t	vmovaps	%%ymm8 ,     (%%r10,%%r9)	\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7	\n\t	vmovaps	%%ymm9 ,0x020(%%r10,%%r9)	\n\t"\
		"vmovaps	%%ymm6,%%ymm4		\n\t	vmovaps	0x080(%%r8),%%ymm8	 		\n\t"\
		"vmovaps	%%ymm7,%%ymm5		\n\t	vmovaps	0x0a0(%%r8),%%ymm9 			\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vmovaps	0x080(%%rdx),%%ymm14		\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t	vmovaps	0x0a0(%%rdx),%%ymm15		\n\t"\
		"vmulpd	%%ymm1,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm14,%%ymm12				\n\t"\
		"vmulpd	%%ymm1,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm5,0x020(%%r12)	\n\t	vmulpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t"/* Spill 2 */\
		"vmovaps	%%ymm4,     (%%r12)	\n\t	vmulpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm0	\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm1	\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	     (%%rbx),%%ymm6	\n\t	vmovaps	%%ymm13,0x020(%%r12,%%r9)	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7	\n\t	vmovaps	%%ymm12,     (%%r12,%%r9)	\n\t"\
		"vmovaps	%%ymm6,%%ymm4		\n\t	vmovaps	0x040(%%r8),%%ymm8	 		\n\t"\
		"vmovaps	%%ymm7,%%ymm5		\n\t	vmovaps	0x060(%%r8),%%ymm9 			\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vmovaps	0x080(%%rbx),%%ymm14		\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t	vmovaps	0x0a0(%%rbx),%%ymm15		\n\t"\
		"vmulpd	%%ymm1,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm14,%%ymm12				\n\t"\
		"vmulpd	%%ymm1,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm15,%%ymm13				\n\t"\
		"vmovaps	     (%%r12),%%ymm0	\n\t	vmulpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"/* Restore 2 */\
		"vmovaps	0x020(%%r12),%%ymm1	\n\t	vmulpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t	 movq	%[__two],%%rsi	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm5,%%ymm7		\n\t	vmovaps	     (%%r12,%%r9),%%ymm8	\n\t"\
		"vmovaps	%%ymm4,%%ymm6		\n\t	vmovaps	0x020(%%r12,%%r9),%%ymm9	\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5	\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm13,%%ymm15				\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm12,%%ymm14				\n\t"\
		/* Finish radix-4 butterfly and store results: */\
		"vmovaps	     (%%r10),%%ymm0	\n\t	vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"/* Restore 1 */\
		"vmovaps	0x020(%%r10),%%ymm1	\n\t	vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0	\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1	\n\t	vmovaps	     (%%r10,%%r9),%%ymm8	\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmovaps	0x020(%%r10,%%r9),%%ymm9	\n\t"\
		"vmovaps	%%ymm0,     (%%r12)	\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"	/* 2.0, shared by both columns ... moved +- until found best cycle count: */\
		"vmovaps	%%ymm2,     (%%r11)	\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10		\n\t	vmovaps	(%%rsi),%%ymm0	\n\t"\
		"vmovaps	%%ymm1,0x020(%%r12)	\n\t	vsubpd	%%ymm15,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm3,0x020(%%r13)	\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	%%ymm0,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm8 ,     (%%r12,%%r9)	\n\t"\
		"vmulpd	%%ymm0,%%ymm5,%%ymm5	\n\t	vmovaps	%%ymm10,     (%%r11,%%r9)	\n\t"\
		"vmulpd	%%ymm0,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm9 ,0x020(%%r12,%%r9)	\n\t"\
		"vmulpd	%%ymm0,%%ymm4,%%ymm4	\n\t	vmovaps	%%ymm11,0x020(%%r13,%%r9)	\n\t"\
		"vaddpd	(%%r12),%%ymm6,%%ymm6	\n\t	vmulpd	%%ymm0 ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm2 ,%%ymm5,%%ymm5	\n\t	vmulpd	%%ymm0 ,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm1 ,%%ymm7,%%ymm7	\n\t	vmulpd	%%ymm0 ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm3 ,%%ymm4,%%ymm4	\n\t	vmulpd	%%ymm0 ,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm6,     (%%r10)	\n\t	vaddpd (%%r12,%%r9),%%ymm14,%%ymm14	\n\t"/* don't need reload-from-mem of ymm8/0x180(%%r12,%%r9) as we do in lcol, but 1 cycle faster with it. [!?] */\
		"vmovaps	%%ymm5,     (%%r13)	\n\t	vaddpd      %%ymm10,%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm7,0x020(%%r10)	\n\t	vaddpd      %%ymm9 ,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm4,0x020(%%r11)	\n\t	vaddpd      %%ymm11,%%ymm12,%%ymm12	\n\t"\
		"										vmovaps	%%ymm14,     (%%r10,%%r9)	\n\t"\
		"										vmovaps	%%ymm13,     (%%r13,%%r9)	\n\t"\
		"										vmovaps	%%ymm15,0x020(%%r10,%%r9)	\n\t"\
		"										vmovaps	%%ymm12,0x020(%%r11,%%r9)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r0] "m" (Xr0)\
		 ,[__two] "m" (Xtwo)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_DIT_TWIDDLE_V2(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr0,Xtwo,Xcc0,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/*
		// Pass 1:
		j = p4*8;	// Set j to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k = 0x100;	// Set k to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		// 0x0c0 = literal-bytewise address offset between the pointers to the first and second DFT's twiddles:
		i0 = (vec_dbl *)a; i1 = (vec_dbl *)(a+p1); i2 = (vec_dbl *)(a+p2); i3 = (vec_dbl *)(a+p3);
		o0 = r00; o1 = r00+2; o2 = r00+4; o3 = r00+6;
		c_tmp = cc0+6;	// c2,1,3
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,0x0c0, o0,o1,o2,o3,k)
	*/\
		"movq	%[__two],%%rsi 				\n\t	movslq	%[__p4],%%rdi		\n\t"\
		"movq	%[__add0],%%rax				\n\t	movq	%[__r0],%%r10		\n\t"\
		"movslq	%[__p1],%%rbx				\n\t	leaq	(%%rax,%%rbx,8),%%rbx	\n\t"\
		"movslq	%[__p2],%%rcx				\n\t	leaq	(%%rax,%%rcx,8),%%rcx	\n\t"\
		"movslq	%[__p3],%%rdx				\n\t	leaq	(%%rax,%%rdx,8),%%rdx	\n\t"\
	"vmovaps	(%%rsi),%%ymm15		\n\t"/* two */"	shlq	$3,%%rdi			\n\t"\
		"vmovaps	     (%%rax),%%ymm0		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t"\
		"vmovaps	     (%%rbx),%%ymm2		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3		\n\t"\
		"vmovaps	     (%%rcx),%%ymm4		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5		\n\t"\
		"vmovaps	     (%%rdx),%%ymm6		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7		\n\t	leaq	0x020(%%rdx,%%rdi),%%r11	\n\t"\
		"vsubpd	%%ymm2 ,%%ymm0,%%ymm0		\n\t	vmovaps	     (%%rax,%%rdi),%%ymm8 	\n\t"\
		"vsubpd	%%ymm3 ,%%ymm1,%%ymm1		\n\t	vmovaps	0x020(%%rax,%%rdi),%%ymm9 	\n\t"\
		"vsubpd	%%ymm6 ,%%ymm4,%%ymm4		\n\t	vmovaps	     (%%rbx,%%rdi),%%ymm10	\n\t"\
		"vsubpd	%%ymm7 ,%%ymm5,%%ymm5		\n\t	vmovaps	0x020(%%rbx,%%rdi),%%ymm11	\n\t"\
		"vmulpd	%%ymm15,%%ymm2,%%ymm2		\n\t	vmovaps	     (%%rcx,%%rdi),%%ymm12	\n\t"\
		"vmulpd	%%ymm15,%%ymm3,%%ymm3		\n\t	vmovaps	0x020(%%rcx,%%rdi),%%ymm13	\n\t"\
		"vmulpd	%%ymm15,%%ymm6,%%ymm6		\n\t	vmovaps	     (%%rdx,%%rdi),%%ymm14	\n\t"\
		"vmulpd	%%ymm15,%%ymm7,%%ymm7		\n\t"/*	vmovaps	0x020(%%rdx,%%rdi),%%ymm15	\n\t"*/\
		"vaddpd	%%ymm0 ,%%ymm2,%%ymm2		\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm1 ,%%ymm3,%%ymm3		\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm4 ,%%ymm6,%%ymm6		\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm7,%%ymm7		\n\t	vsubpd	(%%r11),%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm6 ,%%ymm2,%%ymm2		\n\t	vmulpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7 ,%%ymm3,%%ymm3		\n\t	vmulpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm0,%%ymm0		\n\t	vmulpd	%%ymm15,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm4 ,%%ymm1,%%ymm1		\n\t	vmulpd	(%%r11),%%ymm15,%%ymm15		\n\t"\
		"vmulpd	(%%rsi),%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	(%%rsi),%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	(%%rsi),%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	(%%rsi),%%ymm4,%%ymm4		\n\t	vaddpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm2 ,%%ymm6,%%ymm6		\n\t	vsubpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm3 ,%%ymm7,%%ymm7		\n\t	vsubpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm0 ,%%ymm5,%%ymm5		\n\t	vsubpd	%%ymm13,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm1 ,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm6,    (%%r10)		\n\t	vmovaps	%%ymm7,0x020(%%r10)		\n\t"/* Write B0 to free up 2 regs */\
		"vmovaps	(%%rsi),%%ymm6			\n\t	vmulpd	%%ymm6 ,%%ymm14,%%ymm14		\n\t"\
		"											vmulpd	%%ymm6 ,%%ymm15,%%ymm15		\n\t"\
		"											vmulpd	%%ymm6 ,%%ymm13,%%ymm13		\n\t"\
		/* B2 = t0*~w2 = t0*[c2-I.s1] */"			vmulpd	%%ymm6 ,%%ymm12,%%ymm12		\n\t	movq	%[__cc0],%%rsi 	\n\t"\
		"vmovaps	%%ymm2,%%ymm6			\n\t	vaddpd	%%ymm10,%%ymm14,%%ymm14		\n\t	addq	$0x0c0,%%rsi 	\n\t"/* cc0 + 6 */\
		"vmovaps	%%ymm3,%%ymm7			\n\t	vaddpd	%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	     (%%rsi),%%ymm6,%%ymm6	\n\t	vaddpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t	vaddpd	%%ymm9 ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	     (%%rsi),%%ymm7,%%ymm7	\n\t	vmovaps %%ymm14,0x100(%%r10)	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2	\n\t	vmovaps %%ymm15,0x120(%%r10)	\n\t"\
		"vaddpd	     %%ymm3 ,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm10,%%ymm14			\n\t"\
		"vsubpd	     %%ymm2 ,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm11,%%ymm15			\n\t"\
		"vmovaps	%%ymm6,0x080(%%r10)		\n\t	vmulpd	0x0c0(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm7,0x0a0(%%r10)		\n\t	vmulpd	0x0e0(%%rsi),%%ymm11,%%ymm11	\n\t"\
		/* B1 = t1*~w1 = t1*[c1-I.s1] */"			vmulpd	0x0c0(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm2		\n\t	vmulpd	0x0e0(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm3		\n\t	vaddpd	     %%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm6			\n\t	vsubpd	     %%ymm10,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm1,%%ymm7			\n\t	vmovaps %%ymm14,0x180(%%r10)	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6		\n\t	vmovaps %%ymm15,0x1a0(%%r10)	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vmovaps	%%ymm13,%%ymm14			\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7		\n\t	vmovaps	%%ymm9 ,%%ymm15			\n\t"\
		"vmulpd	%%ymm3,%%ymm5,%%ymm5		\n\t	vmulpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps	0x080(%%rsi),%%ymm2		\n\t	vmulpd	0x120(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm3		\n\t	vmulpd	0x100(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6		\n\t	vmulpd	0x120(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7		\n\t	vaddpd	     %%ymm9 ,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm6,0x040(%%r10)		\n\t	vsubpd	     %%ymm13,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm7,0x060(%%r10)		\n\t	vmovaps %%ymm14,0x140(%%r10)	\n\t"\
		/* B3 = t3*~w3 = t3*[c3-I.s3] */"			vmovaps %%ymm15,0x160(%%r10)	\n\t"\
		"vmovaps	%%ymm0,%%ymm6			\n\t	vmovaps	%%ymm8 ,%%ymm14			\n\t"\
		"vmovaps	%%ymm4,%%ymm7			\n\t	vmovaps	%%ymm12,%%ymm15			\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6		\n\t	vmulpd	0x140(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm3,%%ymm4,%%ymm4		\n\t	vmulpd	0x160(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7		\n\t	vmulpd	0x140(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0		\n\t	vmulpd	0x160(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	     %%ymm12,%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm0,%%ymm7,%%ymm7		\n\t	vsubpd	     %%ymm8 ,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm6,0x0c0(%%r10)		\n\t	vmovaps	%%ymm14,0x1c0(%%r10)	\n\t"\
		"vmovaps	%%ymm7,0x0e0(%%r10)		\n\t	vmovaps	%%ymm15,0x1e0(%%r10)	\n\t"\
	/*
		i0 = (vec_dbl *)(a+p8); i1 = (vec_dbl *)(a+p9); i2 = (vec_dbl *)(a+pA); i3 = (vec_dbl *)(a+pB);
		o0 += 16; o1 += 16; o2 += 16; o3 += 16;
		c_tmp += 12;		// cA,9,B
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,0x0c0, o0,o1,o2,o3,k)
	*/\
		"movq	%[__two],%%rsi 		\n\t	movslq	%[__p8],%%r11		\n\t"\
		"leaq (%%rax,%%r11,8),%%rax	\n\t	addq	 $0x200,%%r10		\n\t"/* r0 + 16 */\
		"leaq (%%rbx,%%r11,8),%%rbx	\n\t"\
		"leaq (%%rcx,%%r11,8),%%rcx	\n\t"\
		"leaq (%%rdx,%%r11,8),%%rdx	\n\t"\
	"vmovaps	(%%rsi),%%ymm15	\n\t"/* two */"	"\
		"vmovaps	     (%%rax),%%ymm0		\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t"\
		"vmovaps	     (%%rbx),%%ymm2		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm3		\n\t"\
		"vmovaps	     (%%rcx),%%ymm4		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5		\n\t"\
		"vmovaps	     (%%rdx),%%ymm6		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7		\n\t	leaq	0x020(%%rdx,%%rdi),%%r11	\n\t"\
		"vsubpd	%%ymm2 ,%%ymm0,%%ymm0		\n\t	vmovaps	     (%%rax,%%rdi),%%ymm8 	\n\t"\
		"vsubpd	%%ymm3 ,%%ymm1,%%ymm1		\n\t	vmovaps	0x020(%%rax,%%rdi),%%ymm9 	\n\t"\
		"vsubpd	%%ymm6 ,%%ymm4,%%ymm4		\n\t	vmovaps	     (%%rbx,%%rdi),%%ymm10	\n\t"\
		"vsubpd	%%ymm7 ,%%ymm5,%%ymm5		\n\t	vmovaps	0x020(%%rbx,%%rdi),%%ymm11	\n\t"\
		"vmulpd	%%ymm15,%%ymm2,%%ymm2		\n\t	vmovaps	     (%%rcx,%%rdi),%%ymm12	\n\t"\
		"vmulpd	%%ymm15,%%ymm3,%%ymm3		\n\t	vmovaps	0x020(%%rcx,%%rdi),%%ymm13	\n\t"\
		"vmulpd	%%ymm15,%%ymm6,%%ymm6		\n\t	vmovaps	     (%%rdx,%%rdi),%%ymm14	\n\t"\
		"vmulpd	%%ymm15,%%ymm7,%%ymm7		\n\t"/*	vmovaps	0x020(%%rdx,%%rdi),%%ymm15	\n\t"*/\
		"vaddpd	%%ymm0 ,%%ymm2,%%ymm2		\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm1 ,%%ymm3,%%ymm3		\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm4 ,%%ymm6,%%ymm6		\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm7,%%ymm7		\n\t	vsubpd	(%%r11),%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm6 ,%%ymm2,%%ymm2		\n\t	vmulpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7 ,%%ymm3,%%ymm3		\n\t	vmulpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm0,%%ymm0		\n\t	vmulpd	%%ymm15,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm4 ,%%ymm1,%%ymm1		\n\t	vmulpd	(%%r11),%%ymm15,%%ymm15		\n\t"\
		"vmulpd	(%%rsi),%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	(%%rsi),%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	(%%rsi),%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	(%%rsi),%%ymm4,%%ymm4		\n\t	vaddpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm2 ,%%ymm6,%%ymm6		\n\t	vsubpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm3 ,%%ymm7,%%ymm7		\n\t	vsubpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm0 ,%%ymm5,%%ymm5		\n\t	vsubpd	%%ymm13,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm1 ,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm6,    (%%r10)		\n\t	vmovaps	%%ymm7,0x020(%%r10)		\n\t"/* Write B0 to free up 2 regs */\
		"vmovaps	(%%rsi),%%ymm6			\n\t	vmulpd	%%ymm6 ,%%ymm14,%%ymm14		\n\t"\
		"											vmulpd	%%ymm6 ,%%ymm15,%%ymm15		\n\t"\
		"											vmulpd	%%ymm6 ,%%ymm13,%%ymm13		\n\t"\
		/* B2 = t0*~w2 = t0*[c2-I.s1] */"			vmulpd	%%ymm6 ,%%ymm12,%%ymm12		\n\t	movq	%[__cc0],%%rsi 	\n\t"\
		"vmovaps	%%ymm2,%%ymm6			\n\t	vaddpd	%%ymm10,%%ymm14,%%ymm14		\n\t	addq	$0x240,%%rsi 	\n\t"/* cc0 + 18 */\
		"vmovaps	%%ymm3,%%ymm7			\n\t	vaddpd	%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	     (%%rsi),%%ymm6,%%ymm6	\n\t	vaddpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t	vaddpd	%%ymm9 ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	     (%%rsi),%%ymm7,%%ymm7	\n\t	vmovaps %%ymm14,0x100(%%r10)	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2	\n\t	vmovaps %%ymm15,0x120(%%r10)	\n\t"\
		"vaddpd	     %%ymm3 ,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm10,%%ymm14			\n\t"\
		"vsubpd	     %%ymm2 ,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm11,%%ymm15			\n\t"\
		"vmovaps	%%ymm6,0x080(%%r10)		\n\t	vmulpd	0x0c0(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm7,0x0a0(%%r10)		\n\t	vmulpd	0x0e0(%%rsi),%%ymm11,%%ymm11	\n\t"\
		/* B1 = t1*~w1 = t1*[c1-I.s1] */"			vmulpd	0x0c0(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm2		\n\t	vmulpd	0x0e0(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm3		\n\t	vaddpd	     %%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm5,%%ymm6			\n\t	vsubpd	     %%ymm10,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm1,%%ymm7			\n\t	vmovaps %%ymm14,0x180(%%r10)	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6		\n\t	vmovaps %%ymm15,0x1a0(%%r10)	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vmovaps	%%ymm13,%%ymm14			\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7		\n\t	vmovaps	%%ymm9 ,%%ymm15			\n\t"\
		"vmulpd	%%ymm3,%%ymm5,%%ymm5		\n\t	vmulpd	0x100(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps	0x080(%%rsi),%%ymm2		\n\t	vmulpd	0x120(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm3		\n\t	vmulpd	0x100(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6		\n\t	vmulpd	0x120(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7		\n\t	vaddpd	     %%ymm9 ,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm6,0x040(%%r10)		\n\t	vsubpd	     %%ymm13,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm7,0x060(%%r10)		\n\t	vmovaps %%ymm14,0x140(%%r10)	\n\t"\
		/* B3 = t3*~w3 = t3*[c3-I.s3] */"			vmovaps %%ymm15,0x160(%%r10)	\n\t"\
		"vmovaps	%%ymm0,%%ymm6			\n\t	vmovaps	%%ymm8 ,%%ymm14			\n\t"\
		"vmovaps	%%ymm4,%%ymm7			\n\t	vmovaps	%%ymm12,%%ymm15			\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6		\n\t	vmulpd	0x140(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm3,%%ymm4,%%ymm4		\n\t	vmulpd	0x160(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7		\n\t	vmulpd	0x140(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0		\n\t	vmulpd	0x160(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	     %%ymm12,%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm0,%%ymm7,%%ymm7		\n\t	vsubpd	     %%ymm8 ,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm6,0x0c0(%%r10)		\n\t	vmovaps	%%ymm14,0x1c0(%%r10)	\n\t"\
		"vmovaps	%%ymm7,0x0e0(%%r10)		\n\t	vmovaps	%%ymm15,0x1e0(%%r10)	\n\t"\
	/*
		// Pass 2:
		j = 0x080;	// Set j to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k = p2*8;	// Set k to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		c_tmp = cc0;	// c8,4,C
		// Pass1-DFTs all use same twiddle-triplet, so lit-byte address offset between ptrs to 1st,2nd DFT's twiddles = 0
		i0 = r00; i1 = r00+8; i2 = r00+16; i3 = r00+24;
		o0 = (vec_dbl *)a; o1 = (vec_dbl *)(a+p4); o2 = (vec_dbl *)(a+p8); o3 = (vec_dbl *)(a+pC);
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,   0, o0,o1,o2,o3,k)
	*/\
		"movq	%[__two],%%rsi 				\n\t	movslq	%[__p2],%%r9		\n\t"\
		"movq	%[__add0],%%r10				\n\t	movq	%[__r0],%%rax		\n\t"\
		"movslq	%[__p4]  ,%%r11				\n\t	leaq	(%%r10,%%r11,8),%%r11	\n\t"\
		"movslq	%[__p8]  ,%%r12				\n\t	leaq	(%%r10,%%r12,8),%%r12	\n\t"\
		"movslq	%[__p12] ,%%r13				\n\t	leaq	(%%r10,%%r13,8),%%r13	\n\t"\
		"vmovaps	(%%rsi),%%ymm15			\n\t"/* two */\
		"vmovaps	     (%%rax),%%ymm0		\n\t"/* r0 */"	leaq	0x080(%%rax),%%rbx	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2		\n\t"/* r0 + 8 */\
		"vmovaps	0x120(%%rax),%%ymm3		\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4		\n\t"/* r0 + 16 */\
		"vmovaps	0x220(%%rax),%%ymm5		\n\t"\
		"vmovaps	0x300(%%rax),%%ymm6		\n\t"/* r0 + 24 */\
		"vmovaps	0x320(%%rax),%%ymm7		\n\t	leaq	0x320(%%rbx),%%r8		\n\t"\
		"vsubpd	%%ymm2 ,%%ymm0,%%ymm0		\n\t	vmovaps	     (%%rbx),%%ymm8 	\n\t"\
		"vsubpd	%%ymm3 ,%%ymm1,%%ymm1		\n\t	vmovaps	0x020(%%rbx),%%ymm9 	\n\t"\
		"vsubpd	%%ymm6 ,%%ymm4,%%ymm4		\n\t	vmovaps	0x100(%%rbx),%%ymm10	\n\t"\
		"vsubpd	%%ymm7 ,%%ymm5,%%ymm5		\n\t	vmovaps	0x120(%%rbx),%%ymm11	\n\t"\
		"vmulpd	%%ymm15,%%ymm2,%%ymm2		\n\t	vmovaps	0x200(%%rbx),%%ymm12	\n\t"\
		"vmulpd	%%ymm15,%%ymm3,%%ymm3		\n\t	vmovaps	0x220(%%rbx),%%ymm13	\n\t"\
		"vmulpd	%%ymm15,%%ymm6,%%ymm6		\n\t	vmovaps	0x300(%%rbx),%%ymm14	\n\t"\
		"vmulpd	%%ymm15,%%ymm7,%%ymm7		\n\t"/*	vmovaps	0x320(%%rbx),%%ymm15	\n\t"*/\
		"vaddpd	%%ymm0 ,%%ymm2,%%ymm2		\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm1 ,%%ymm3,%%ymm3		\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm4 ,%%ymm6,%%ymm6		\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm7,%%ymm7		\n\t	vsubpd	(%%r8) ,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm6 ,%%ymm2,%%ymm2		\n\t	vmulpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7 ,%%ymm3,%%ymm3		\n\t	vmulpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm0,%%ymm0		\n\t	vmulpd	%%ymm15,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm4 ,%%ymm1,%%ymm1		\n\t	vmulpd	(%%r8) ,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	(%%rsi),%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	(%%rsi),%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	(%%rsi),%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	(%%rsi),%%ymm4,%%ymm4		\n\t	vaddpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm2 ,%%ymm6,%%ymm6		\n\t	vsubpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm3 ,%%ymm7,%%ymm7		\n\t	vsubpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm0 ,%%ymm5,%%ymm5		\n\t	vsubpd	%%ymm13,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm1 ,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm6,    (%%r10)	\n\t	vmovaps	%%ymm7,0x020(%%r10)	\n\t"/* Write B0 to free up 2 regs */\
		"vmovaps	(%%rsi),%%ymm6			\n\t	vmulpd	%%ymm6 ,%%ymm14,%%ymm14		\n\t"\
		"											vmulpd	%%ymm6 ,%%ymm15,%%ymm15		\n\t"\
		"											vmulpd	%%ymm6 ,%%ymm13,%%ymm13		\n\t"\
		/* B2 = t0*~w2 = t0*[c2-I.s1] */"			vmulpd	%%ymm6 ,%%ymm12,%%ymm12		\n\t	movq	%[__cc0],%%rsi 	\n\t"\
		"vmovaps	%%ymm2,%%ymm6			\n\t	vaddpd	%%ymm10,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm3,%%ymm7			\n\t	vaddpd	%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	     (%%rsi),%%ymm6,%%ymm6	\n\t	vaddpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t	vaddpd	%%ymm9 ,%%ymm12,%%ymm12		\n\t	shlq	$3,%%r9	\n\t"/* p2*8 */\
		"vmulpd	     (%%rsi),%%ymm7,%%ymm7	\n\t	vmovaps %%ymm14,     (%%r10,%%r9)	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2	\n\t	vmovaps %%ymm15,0x020(%%r10,%%r9)	\n\t"\
		"vaddpd	     %%ymm3 ,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm10,%%ymm14				\n\t"\
		"vsubpd	     %%ymm2 ,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm11,%%ymm15				\n\t"\
		"vmovaps	%%ymm6,     (%%r12)		\n\t	vmulpd	     (%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm7,0x020(%%r12)		\n\t	vmulpd	0x020(%%rsi),%%ymm11,%%ymm11	\n\t"\
		/* B1 = t1*~w1 = t1*[c1-I.s1] */"			vmulpd	     (%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm2		\n\t	vmulpd	0x020(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm3		\n\t	vaddpd	     %%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	     %%ymm5 ,%%ymm6		\n\t	vsubpd	     %%ymm10,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	     %%ymm1 ,%%ymm7		\n\t	vmovaps %%ymm14,     (%%r12,%%r9)	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6		\n\t	vmovaps %%ymm15,0x020(%%r12,%%r9)	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vmovaps	%%ymm13,%%ymm14				\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7		\n\t	vmovaps	%%ymm9 ,%%ymm15				\n\t"\
		"vmulpd	%%ymm3,%%ymm5,%%ymm5		\n\t	vmulpd	0x040(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps	0x080(%%rsi),%%ymm2		\n\t	vmulpd	0x060(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm3		\n\t	vmulpd	0x040(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6		\n\t	vmulpd	0x060(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7		\n\t	vaddpd	     %%ymm9 ,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm6,     (%%r11)		\n\t	vsubpd	     %%ymm13,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm7,0x020(%%r11)		\n\t	vmovaps %%ymm14,     (%%r11,%%r9)	\n\t"\
		/* B3 = t3*~w3 = t3*[c3-I.s3] */"			vmovaps %%ymm15,0x020(%%r11,%%r9)	\n\t"\
		"vmovaps	%%ymm0,%%ymm6			\n\t	vmovaps	%%ymm8 ,%%ymm14				\n\t"\
		"vmovaps	%%ymm4,%%ymm7			\n\t	vmovaps	%%ymm12,%%ymm15				\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6		\n\t	vmulpd	0x080(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm3,%%ymm4,%%ymm4		\n\t	vmulpd	0x0a0(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7		\n\t	vmulpd	0x080(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0		\n\t	vmulpd	0x0a0(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	     %%ymm12,%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm0,%%ymm7,%%ymm7		\n\t	vsubpd	     %%ymm8 ,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm6,     (%%r13)		\n\t	vmovaps	%%ymm14,     (%%r13,%%r9)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%r13)		\n\t	vmovaps	%%ymm15,0x020(%%r13,%%r9)	\n\t"\
	/*
		i0 += 2; i1 += 2; i2 += 2; i3 += 2;
		o0 = (vec_dbl *)(a+p1); o1 = (vec_dbl *)(a+p5); o2 = (vec_dbl *)(a+p9); o3 = (vec_dbl *)(a+pD);
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,   0, o0,o1,o2,o3,k)
	*/\
		"movq	%[__two],%%rsi 		\n\t	movslq	%[__p1],%%rdi		\n\t"\
		"leaq (%%r10,%%rdi,8),%%r10	\n\t	addq	$0x040,%%rax		\n\t"\
		"leaq (%%r11,%%rdi,8),%%r11	\n\t"\
		"leaq (%%r12,%%rdi,8),%%r12	\n\t"\
		"leaq (%%r13,%%rdi,8),%%r13	\n\t"\
		"vmovaps	(%%rsi),%%ymm15			\n\t"/* two */\
		"vmovaps	     (%%rax),%%ymm0		\n\t"/* r0 */"	leaq	0x080(%%rax),%%rbx	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1		\n\t"\
		"vmovaps	0x100(%%rax),%%ymm2		\n\t"/* r0 + 8 */\
		"vmovaps	0x120(%%rax),%%ymm3		\n\t"\
		"vmovaps	0x200(%%rax),%%ymm4		\n\t"/* r0 + 16 */\
		"vmovaps	0x220(%%rax),%%ymm5		\n\t"\
		"vmovaps	0x300(%%rax),%%ymm6		\n\t"/* r0 + 24 */\
		"vmovaps	0x320(%%rax),%%ymm7		\n\t	leaq	0x320(%%rbx),%%r8		\n\t"\
		"vsubpd	%%ymm2 ,%%ymm0,%%ymm0		\n\t	vmovaps	     (%%rbx),%%ymm8 	\n\t"\
		"vsubpd	%%ymm3 ,%%ymm1,%%ymm1		\n\t	vmovaps	0x020(%%rbx),%%ymm9 	\n\t"\
		"vsubpd	%%ymm6 ,%%ymm4,%%ymm4		\n\t	vmovaps	0x100(%%rbx),%%ymm10	\n\t"\
		"vsubpd	%%ymm7 ,%%ymm5,%%ymm5		\n\t	vmovaps	0x120(%%rbx),%%ymm11	\n\t"\
		"vmulpd	%%ymm15,%%ymm2,%%ymm2		\n\t	vmovaps	0x200(%%rbx),%%ymm12	\n\t"\
		"vmulpd	%%ymm15,%%ymm3,%%ymm3		\n\t	vmovaps	0x220(%%rbx),%%ymm13	\n\t"\
		"vmulpd	%%ymm15,%%ymm6,%%ymm6		\n\t	vmovaps	0x300(%%rbx),%%ymm14	\n\t"\
		"vmulpd	%%ymm15,%%ymm7,%%ymm7		\n\t"/*	vmovaps	0x320(%%rbx),%%ymm15	\n\t"*/\
		"vaddpd	%%ymm0 ,%%ymm2,%%ymm2		\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm1 ,%%ymm3,%%ymm3		\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm4 ,%%ymm6,%%ymm6		\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm5 ,%%ymm7,%%ymm7		\n\t	vsubpd	(%%r8) ,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm6 ,%%ymm2,%%ymm2		\n\t	vmulpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7 ,%%ymm3,%%ymm3		\n\t	vmulpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm5 ,%%ymm0,%%ymm0		\n\t	vmulpd	%%ymm15,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm4 ,%%ymm1,%%ymm1		\n\t	vmulpd	(%%r8) ,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	(%%rsi),%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	(%%rsi),%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	(%%rsi),%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	(%%rsi),%%ymm4,%%ymm4		\n\t	vaddpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm2 ,%%ymm6,%%ymm6		\n\t	vsubpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm3 ,%%ymm7,%%ymm7		\n\t	vsubpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm0 ,%%ymm5,%%ymm5		\n\t	vsubpd	%%ymm13,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm1 ,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm12,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm6,    (%%r10)	\n\t	vmovaps	%%ymm7,0x020(%%r10)	\n\t"/* Write B0 to free up 2 regs */\
		"vmovaps	(%%rsi),%%ymm6			\n\t	vmulpd	%%ymm6 ,%%ymm14,%%ymm14		\n\t"\
		"											vmulpd	%%ymm6 ,%%ymm15,%%ymm15		\n\t"\
		"											vmulpd	%%ymm6 ,%%ymm13,%%ymm13		\n\t"\
		/* B2 = t0*~w2 = t0*[c2-I.s1] */"			vmulpd	%%ymm6 ,%%ymm12,%%ymm12		\n\t	movq	%[__cc0],%%rsi 	\n\t"\
		"vmovaps	%%ymm2,%%ymm6			\n\t	vaddpd	%%ymm10,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm3,%%ymm7			\n\t	vaddpd	%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"vmulpd	     (%%rsi),%%ymm6,%%ymm6	\n\t	vaddpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3	\n\t	vaddpd	%%ymm9 ,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	     (%%rsi),%%ymm7,%%ymm7	\n\t	vmovaps %%ymm14,     (%%r10,%%r9)	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2	\n\t	vmovaps %%ymm15,0x020(%%r10,%%r9)	\n\t"\
		"vaddpd	     %%ymm3 ,%%ymm6,%%ymm6	\n\t	vmovaps	%%ymm10,%%ymm14				\n\t"\
		"vsubpd	     %%ymm2 ,%%ymm7,%%ymm7	\n\t	vmovaps	%%ymm11,%%ymm15				\n\t"\
		"vmovaps	%%ymm6,     (%%r12)		\n\t	vmulpd	     (%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm7,0x020(%%r12)		\n\t	vmulpd	0x020(%%rsi),%%ymm11,%%ymm11	\n\t"\
		/* B1 = t1*~w1 = t1*[c1-I.s1] */"			vmulpd	     (%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm2		\n\t	vmulpd	0x020(%%rsi),%%ymm10,%%ymm10	\n\t"\
		"vmovaps	0x060(%%rsi),%%ymm3		\n\t	vaddpd	     %%ymm11,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	     %%ymm5 ,%%ymm6		\n\t	vsubpd	     %%ymm10,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	     %%ymm1 ,%%ymm7		\n\t	vmovaps %%ymm14,     (%%r12,%%r9)	\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6		\n\t	vmovaps %%ymm15,0x020(%%r12,%%r9)	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vmovaps	%%ymm13,%%ymm14				\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7		\n\t	vmovaps	%%ymm9 ,%%ymm15				\n\t"\
		"vmulpd	%%ymm3,%%ymm5,%%ymm5		\n\t	vmulpd	0x040(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps	0x080(%%rsi),%%ymm2		\n\t	vmulpd	0x060(%%rsi),%%ymm9 ,%%ymm9 	\n\t"\
		"vmovaps	0x0a0(%%rsi),%%ymm3		\n\t	vmulpd	0x040(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6		\n\t	vmulpd	0x060(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7		\n\t	vaddpd	     %%ymm9 ,%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm6,     (%%r11)		\n\t	vsubpd	     %%ymm13,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm7,0x020(%%r11)		\n\t	vmovaps %%ymm14,     (%%r11,%%r9)	\n\t"\
		/* B3 = t3*~w3 = t3*[c3-I.s3] */"			vmovaps %%ymm15,0x020(%%r11,%%r9)	\n\t"\
		"vmovaps	%%ymm0,%%ymm6			\n\t	vmovaps	%%ymm8 ,%%ymm14				\n\t"\
		"vmovaps	%%ymm4,%%ymm7			\n\t	vmovaps	%%ymm12,%%ymm15				\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6		\n\t	vmulpd	0x080(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmulpd	%%ymm3,%%ymm4,%%ymm4		\n\t	vmulpd	0x0a0(%%rsi),%%ymm12,%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7		\n\t	vmulpd	0x080(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0		\n\t	vmulpd	0x0a0(%%rsi),%%ymm8 ,%%ymm8 	\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	     %%ymm12,%%ymm14,%%ymm14	\n\t"\
		"vsubpd	%%ymm0,%%ymm7,%%ymm7		\n\t	vsubpd	     %%ymm8 ,%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm6,     (%%r13)		\n\t	vmovaps	%%ymm14,     (%%r13,%%r9)	\n\t"\
		"vmovaps	%%ymm7,0x020(%%r13)		\n\t	vmovaps	%%ymm15,0x020(%%r13,%%r9)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r0] "m" (Xr0)\
		 ,[__two] "m" (Xtwo)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_DIF_TWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr1,Xisrt2,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
		"movq	%[__r1],%%r9							\n\t"\
		"movq	%[__add0],%%rax				\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%rax)\n\t"/* [base-address + data-fetch-ahead index] + p0 */\
		"movslq	%[__p4],%%rbx							\n\t"\
		"movslq	%[__p8],%%rcx							\n\t	movq	%[__isrt2],%%r8			\n\t"\
		"movslq	%[__p12],%%rdx							\n\t	movslq	%[__p1],%%r10			\n\t"\
		"shlq	$3,%%rbx								\n\t	movslq	%[__p2],%%rdi			\n\t"\
		"shlq	$3,%%rcx								\n\t	shlq	$3,%%r10				\n\t"\
		"shlq	$3,%%rdx								\n\t	shlq	$3,%%rdi				\n\t"\
	"movq	%%rbx,%%r14	\n\t"	/* Save a copy of p4 pointer offset. Will prefetch from [base-address + data-fetch-ahead index] + [0,p4,p8,p12] on each macro call. */\
		"addq	%%rax,%%rbx								\n\t	addq	$0x460,%%r8	/* two */	\n\t"\
		"addq	%%rax,%%rcx								\n\t	movq	%%r10,%%r11				\n\t"\
		"addq	%%rax,%%rdx								\n\t	movq	%%r10,%%r12				\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r1 ,c0) */		\n\t	movq	%%r10,%%r13				\n\t"\
		"/* Do	the p0,p8 combo: */						\n\t	addq	%%rax,%%r10				\n\t"\
		"leaq	0x460(%%r9),%%rsi	/* c0, from r1 */	\n\t	addq	%%rbx,%%r11				\n\t"\
		"vmovaps	     (%%rax),%%ymm0					\n\t	addq	%%rcx,%%r12				\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	addq	%%rdx,%%r13				\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1					\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE_B(r17,c1) */\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	/* Do	the p0,p8 combo: */		\n\t"\
		"vmovaps	     (%%rsi),%%ymm6		/* c0 */	\n\t	vmovaps	     (%%r10),%%ymm8 		\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm7					\n\t	vmovaps	     (%%r12),%%ymm12		\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	0x020(%%r10),%%ymm9 		\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r10)\n\t"/* ...+p1 */\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmovaps	0x020(%%r12),%%ymm13		\n\t"\
		"vmulpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vmovaps	0x200(%%rsi),%%ymm14	/* c1 */\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3					\n\t	vmovaps	0x220(%%rsi),%%ymm15	\n\t"\
		"vmulpd	%%ymm6,%%ymm1,%%ymm1					\n\t	vmovaps	%%ymm8 ,%%ymm10			\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2					\n\t	vmovaps	%%ymm9 ,%%ymm11			\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmulpd	%%ymm14,%%ymm9 ,%%ymm9 			\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmulpd	%%ymm15,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm3,%%ymm0,%%ymm0					\n\t	vmulpd	%%ymm14,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm2,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm15,%%ymm11,%%ymm11			\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm4,%%ymm4		/* c8 */\n\t	vmovaps	%%ymm12,%%ymm14			\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm7,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15			\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm5,%%ymm5				\n\t	vaddpd	%%ymm10,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm6,%%ymm6				\n\t	vsubpd	%%ymm11,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmulpd	0x240(%%rsi),%%ymm12,%%ymm12	/* c9 */\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmulpd	0x260(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vmulpd	0x240(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vmulpd	0x260(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0					\n\t	vmovaps	%%ymm8 ,%%ymm10			\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1					\n\t	vmovaps	%%ymm9 ,%%ymm11			\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12			\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%rbx)\n\t"/* ...+p4 */\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13			\n\t"\
		"/* Do	the p4,12 combo: */						\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps	     (%%rdx),%%ymm4					\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5					\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	      %%ymm4,%%ymm6					\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11			\n\t"\
		"vmovaps	      %%ymm5,%%ymm7					\n\t	/* Do	the p4,12 combo: */		\n\t"\
		"vmulpd	0x0c0(%%rsi),%%ymm4,%%ymm4		/* c12*/\n\t	vmovaps	     (%%r13),%%ymm12		\n\t"\
		"vmulpd	0x0e0(%%rsi),%%ymm7,%%ymm7				\n\t	vmovaps	0x020(%%r13),%%ymm13		\n\t"\
		"vmulpd	0x0c0(%%rsi),%%ymm5,%%ymm5				\n\t	vmovaps	      %%ymm12,%%ymm14		\n\t"\
		"vmulpd	0x0e0(%%rsi),%%ymm6,%%ymm6				\n\t	vmovaps	      %%ymm13,%%ymm15		\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vmulpd	0x2c0(%%rsi),%%ymm12,%%ymm12	/* c13*/\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vmulpd	0x2e0(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm7					\n\t	vmulpd	0x2c0(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmovaps	     (%%rbx),%%ymm6					\n\t	vmulpd	0x2e0(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm4,     (%%r9)					\n\t	vsubpd	     %%ymm15,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	%%ymm5,0x020(%%r9)					\n\t	vaddpd	     %%ymm14,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	     (%%rbx),%%ymm4					\n\t	vmovaps	     (%%r11),%%ymm14		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5					\n\t	vmovaps	0x020(%%r11),%%ymm15		\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r11)\n\t"/* ...+p5 */\
		"vmulpd	0x080(%%rsi),%%ymm4,%%ymm4		/* c4 */\n\t	vmovaps	%%ymm12,0x200(%%r9)		\n\t"\
		"vmulpd	0x0a0(%%rsi),%%ymm7,%%ymm7				\n\t	vmovaps	%%ymm13,0x220(%%r9)		/* r17 */\n\t"\
		"vmulpd	0x080(%%rsi),%%ymm5,%%ymm5				\n\t	vmovaps	      %%ymm14,%%ymm12		\n\t"\
		"vmulpd	0x0a0(%%rsi),%%ymm6,%%ymm6				\n\t	vmovaps	      %%ymm15,%%ymm13		\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vmulpd	0x280(%%rsi),%%ymm12,%%ymm12	/* c5 */\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vmulpd	0x2a0(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmulpd	0x280(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmulpd	0x2a0(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vsubpd	     (%%r9),%%ymm4,%%ymm4				\n\t	vsubpd	     %%ymm15,%%ymm12,%%ymm12			\n\t"\
		"vsubpd	0x020(%%r9),%%ymm5,%%ymm5				\n\t	vaddpd	     %%ymm14,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	     (%%r9),%%ymm6,%%ymm6				\n\t	vmovaps	%%ymm12,%%ymm14			\n\t"\
		"vaddpd	0x020(%%r9),%%ymm7,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15			\n\t"\
		"/* Finish radix-4 butterfly, store: */			\n\t	vsubpd	0x200(%%r9),%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3					\n\t	vsubpd	0x220(%%r9),%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2					\n\t	vaddpd	0x200(%%r9),%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vaddpd	0x220(%%r9),%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1					\n\t	/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%r9)					\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11			\n\t"\
		"vmovaps	%%ymm2,0x040(%%r9)					\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	%%ymm0,0x080(%%r9)					\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%r9)					\n\t	vsubpd	%%ymm15,%%ymm9 ,%%ymm9 			\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%rcx)\n\t"/* ...+p8 */\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6					\n\t	vmovaps	%%ymm8 ,0x280(%%r9)		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5					\n\t	vmovaps	%%ymm10,0x240(%%r9)		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7					\n\t	vmovaps	%%ymm9 ,0x2a0(%%r9)		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4					\n\t	vmovaps	%%ymm11,0x2e0(%%r9)		\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t	vmulpd	(%%r8) ,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmulpd	(%%r8) ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7					\n\t	vmulpd	(%%r8) ,%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4					\n\t	vmulpd	(%%r8) ,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	%%ymm6,     (%%r9)					\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	%%ymm5,0x0c0(%%r9)					\n\t	vaddpd	%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	%%ymm7,0x020(%%r9)					\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"vmovaps	%%ymm4,0x060(%%r9)					\n\t	vaddpd	%%ymm11,%%ymm12,%%ymm12			\n\t"\
		"addq	%%rdi,%%rax								\n\t	vmovaps	%%ymm14,0x200(%%r9)		\n\t"\
		"addq	%%rdi,%%rbx								\n\t	vmovaps	%%ymm13,0x2c0(%%r9)		\n\t"\
		"addq	%%rdi,%%rcx								\n\t	vmovaps	%%ymm15,0x220(%%r9)		\n\t"\
		"addq	%%rdi,%%rdx								\n\t	vmovaps	%%ymm12,0x260(%%r9)		\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r12)\n\t"/* ...+p9 */\
		"/* SSE2_RADIX4_DIF_4TWIDDLE_B(r9 ,c2) */		\n\t	addq	%%rdi,%%r10				\n\t"\
		"/* Do	the p0,p8 combo: */						\n\t	addq	%%rdi,%%r11				\n\t"\
		"addq	$0x100,%%rsi 		/* c2 */			\n\t	addq	%%rdi,%%r12				\n\t"\
		"vmovaps	     (%%rax),%%ymm0					\n\t	addq	%%rdi,%%r13				\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	/* SSE2_RADIX4_DIF_4TWIDDLE_B(r25,c3) */\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1					\n\t	/* Do	the p0,p8 combo: */		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	     (%%r10),%%ymm8 		\n\t"\
		"vmovaps	     (%%rsi),%%ymm6					\n\t	vmovaps	     (%%r12),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm7					\n\t	vmovaps	0x020(%%r10),%%ymm9 		\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmovaps	0x020(%%r12),%%ymm13		\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmovaps	0x200(%%rsi),%%ymm14	/* c3 */\n\t"\
		"vmulpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vmovaps	0x220(%%rsi),%%ymm15	\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3					\n\t	vmovaps	%%ymm8 ,%%ymm10			\n\t"\
		"vmulpd	%%ymm6,%%ymm1,%%ymm1					\n\t	vmovaps	%%ymm9 ,%%ymm11			\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2					\n\t	vmulpd	 %%ymm14,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmulpd	 %%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmulpd	 %%ymm14,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm3,%%ymm0,%%ymm0					\n\t	vmulpd	 %%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm2,%%ymm1,%%ymm1					\n\t	vmovaps	%%ymm12,%%ymm14			\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm4,%%ymm4	/* c10*/	\n\t	vmovaps	%%ymm13,%%ymm15			\n\t"\
"prefetcht1	%c[__pfetch_dist](%%rax)\n\t"/*..+p2 */\
		"vmulpd	0x060(%%rsi),%%ymm7,%%ymm7				\n\t	vsubpd	 %%ymm11,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm5,%%ymm5				\n\t	vaddpd	 %%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm6,%%ymm6				\n\t	vmulpd	0x240(%%rsi),%%ymm12,%%ymm12	/* c11*/\n\t"\
		"vmovaps	%%ymm0,%%ymm2						\n\t	vmulpd	0x260(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vmovaps	%%ymm1,%%ymm3						\n\t	vmulpd	0x240(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vmulpd	0x260(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vmovaps	%%ymm8 ,%%ymm10			\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0					\n\t	vmovaps	%%ymm9 ,%%ymm11			\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12			\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13			\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t	vaddpd	%%ymm12,%%ymm8 ,%%ymm8 			\n\t"\
		"/* Do	the p4,12 combo: */						\n\t	vaddpd	%%ymm13,%%ymm9 ,%%ymm9 			\n\t"\
		"addq	$0x100,%%r9 			/* r9 */		\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	     (%%rdx),%%ymm4					\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11			\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5					\n\t	/* Do	the p4,12 combo: */		\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	vmovaps	     (%%r13),%%ymm12		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm7					\n\t	vmovaps	0x020(%%r13),%%ymm13		\n\t"\
		"vmulpd	0x0c0(%%rsi),%%ymm4,%%ymm4	/* c14*/	\n\t	vmovaps	     (%%r13),%%ymm14		\n\t"\
		"vmulpd	0x0e0(%%rsi),%%ymm7,%%ymm7				\n\t	vmovaps	0x020(%%r13),%%ymm15		\n\t"\
		"vmulpd	0x0c0(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0x2c0(%%rsi),%%ymm12,%%ymm12	/* c15*/\n\t"\
"prefetcht1	%c[__pfetch_dist](%%r10)\n\t"/*..+p3 */\
		"vmulpd	0x0e0(%%rsi),%%ymm6,%%ymm6				\n\t	vmulpd	0x2e0(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vmulpd	0x2c0(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vmulpd	0x2e0(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm4,     (%%r9)					\n\t	vsubpd	     %%ymm15,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	%%ymm5,0x020(%%r9)					\n\t	vaddpd	     %%ymm14,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	     (%%rbx),%%ymm4					\n\t	vmovaps	%%ymm12,0x200(%%r9)		\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm5					\n\t	vmovaps	%%ymm13,0x220(%%r9)		/* r25*/\n\t"\
		"vmovaps	      %%ymm4,%%ymm6					\n\t	vmovaps	     (%%r11),%%ymm12		\n\t"\
		"vmovaps	      %%ymm5,%%ymm7					\n\t	vmovaps	0x020(%%r11),%%ymm13		\n\t"\
		"vmulpd	0x080(%%rsi),%%ymm4,%%ymm4	/* c6 */	\n\t	vmovaps	%%ymm12,%%ymm14		\n\t"\
		"vmulpd	0x0a0(%%rsi),%%ymm7,%%ymm7				\n\t	vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmulpd	0x080(%%rsi),%%ymm5,%%ymm5				\n\t	vmulpd	0x280(%%rsi),%%ymm12,%%ymm12	/* c7 */\n\t"\
		"vmulpd	0x0a0(%%rsi),%%ymm6,%%ymm6				\n\t	vmulpd	0x2a0(%%rsi),%%ymm15,%%ymm15	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vmulpd	0x280(%%rsi),%%ymm13,%%ymm13	\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vmulpd	0x2a0(%%rsi),%%ymm14,%%ymm14	\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vsubpd	     %%ymm15,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vaddpd	     %%ymm14,%%ymm13,%%ymm13			\n\t"\
		"vsubpd	     (%%r9),%%ymm4,%%ymm4				\n\t	vmovaps	%%ymm12,%%ymm14			\n\t"\
		"vsubpd	0x020(%%r9),%%ymm5,%%ymm5				\n\t	vmovaps	%%ymm13,%%ymm15			\n\t"\
		"vaddpd	     (%%r9),%%ymm6,%%ymm6				\n\t	vsubpd	0x200(%%r9),%%ymm12,%%ymm12		\n\t"\
		"vaddpd	0x020(%%r9),%%ymm7,%%ymm7				\n\t	vsubpd	0x220(%%r9),%%ymm13,%%ymm13		\n\t"\
"prefetcht1	%c[__pfetch_dist](%%rbx)\n\t"/*..+p6 */\
		"/* Finish radix-4 butterfly, store: */			\n\t	vaddpd	0x200(%%r9),%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm4,%%ymm3,%%ymm3					\n\t	vaddpd	0x220(%%r9),%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2					\n\t	/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"vsubpd	%%ymm6,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm12,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm7,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10			\n\t"\
		"vmovaps	%%ymm3,0x0e0(%%r9)					\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps	%%ymm2,0x040(%%r9)					\n\t	vsubpd	%%ymm15,%%ymm9 ,%%ymm9 			\n\t"\
		"vmovaps	%%ymm0,0x080(%%r9)					\n\t	vmovaps	%%ymm11,0x2e0(%%r9)		\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%r9)					\n\t	vmovaps	%%ymm10,0x240(%%r9)		\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4					\n\t	vmovaps	%%ymm8 ,0x280(%%r9)		\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5					\n\t	vmovaps	%%ymm9 ,0x2a0(%%r9)		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6					\n\t	vmulpd	(%%r8) ,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7					\n\t	vmulpd	(%%r8) ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm3,%%ymm4,%%ymm4					\n\t	vmulpd	(%%r8) ,%%ymm15,%%ymm15			\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmulpd	(%%r8) ,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm1,%%ymm7,%%ymm7					\n\t	vaddpd	%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	%%ymm4,0x060(%%r9)					\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"vmovaps	%%ymm5,0x0c0(%%r9)					\n\t	vaddpd	%%ymm11,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	%%ymm6,     (%%r9)					\n\t	vmovaps	%%ymm14,0x200(%%r9)		\n\t"\
		"vmovaps	%%ymm7,0x020(%%r9)					\n\t	vmovaps	%%ymm13,0x2c0(%%r9)		\n\t"\
		"												\n\t	vmovaps	%%ymm15,0x220(%%r9)		\n\t"\
		"												\n\t	vmovaps	%%ymm12,0x260(%%r9)		\n\t"\
		"/*************************************************************************************/\n\t"\
		"/*  And now do four more radix-4 transforms, including the internal twiddle factors: */\n\t"\
		"/*************************************************************************************/\n\t"\
		"/* Block 1 - p2 in rdi, add0+p2,3 already in rax,r10: */	\n\t"\
		"movq	%%r10,%%rbx		/* cpy add0+p3*/		\n\t"\
		"movq	%%rax,%%rcx		/* add0+p2 */			\n\t"\
		"movq	%%r10,%%rdx		/* add0+p3*/			/* Block 3: */				\n\t"\
		"subq	%%rdi,%%rax		/* add0    */			\n\t	movslq	%[__p4],%%r10			\n\t"\
		"subq	%%rdi,%%rbx		/* add0+p1 */			\n\t	movslq	%[__p8],%%rdi			\n\t"\
		"movq	%[__isrt2],%%rsi						\n\t	shlq	$3,%%r10				\n\t"\
		"subq	$0x100,%%r9		/* r1 */				\n\t	shlq	$3,%%rdi				\n\t"\
		"vmovaps	     (%%r9),%%ymm0					\n\t	movq	%%r10,%%r11				\n\t"\
		"vmovaps	0x200(%%r9),%%ymm4					\n\t	movq	%%r10,%%r12				\n\t"\
		"vmovaps	0x020(%%r9),%%ymm1					\n\t	movq	%%r10,%%r13				\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r11)\n\t"/* ...+p7 */\
		"vmovaps	0x220(%%r9),%%ymm5					\n\t	addq	%%rax,%%r10				\n\t"\
		"vmovaps	0x100(%%r9),%%ymm2					\n\t	addq	%%rbx,%%r11				\n\t"\
		"vmovaps	0x300(%%r9),%%ymm6					\n\t	addq	%%rcx,%%r12				\n\t"\
		"vmovaps	0x120(%%r9),%%ymm3					\n\t	addq	%%rdx,%%r13				\n\t"\
		"vmovaps	0x320(%%r9),%%ymm7					\n\t	vmovaps	(%%rsi),%%ymm11	/* isrt2 */	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t	vmovaps	0x280(%%r9),%%ymm12	/* r5 */\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4					\n\t	vmovaps	0x2a0(%%r9),%%ymm13	\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t	vmovaps	0x380(%%r9),%%ymm14	\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5					\n\t	vmovaps	0x3a0(%%r9),%%ymm15	\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm11,%%ymm12,%%ymm12		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6					\n\t	vmulpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3					\n\t	vmulpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7					\n\t	vmulpd	%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2					\n\t	vmovaps	0x080(%%r9),%%ymm8 	\n\t"\
		"vaddpd	%%ymm4,%%ymm6,%%ymm6					\n\t	vmovaps	0x1a0(%%r9),%%ymm11	\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3					\n\t	vmovaps	0x0a0(%%r9),%%ymm9 	\n\t"\
		"vaddpd	%%ymm5,%%ymm7,%%ymm7					\n\t	vmovaps	0x180(%%r9),%%ymm10	\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2					\n\t	vsubpd	%%ymm11,%%ymm8 ,%%ymm8 		\n\t"\
"prefetcht1	%c[__pfetch_dist](%%r12,%%r14)\n\t"/* ..+p10 */\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3					\n\t	vsubpd	%%ymm13,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t	vsubpd	%%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t	vsubpd	%%ymm14,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)					\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)					\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm2,	%%ymm6,	%%ymm6					\n\t	vmulpd	(%%r8) ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm3,	%%ymm7,	%%ymm7					\n\t	vmulpd	(%%r8) ,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm6,     (%%rax)					\n\t	vaddpd	%%ymm8 ,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)					\n\t	vaddpd	%%ymm12,%%ymm13,%%ymm13		\n\t"\
		"vsubpd	%%ymm5,	%%ymm0,	%%ymm0					\n\t	vaddpd	%%ymm9 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm4,	%%ymm1,	%%ymm1					\n\t	vaddpd	%%ymm15,%%ymm14,%%ymm14		\n\t"\
		"vmulpd	(%%r8),	%%ymm5,	%%ymm5					\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	(%%r8),	%%ymm4,	%%ymm4					\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)					\n\t	vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)					\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm0,	%%ymm5,	%%ymm5					\n\t	vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm1,	%%ymm4,	%%ymm4					\n\t	vaddpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)					\n\t	vsubpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)					\n\t	vsubpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"												\n\t	vsubpd	%%ymm13,%%ymm10,%%ymm10		\n\t"\
		"												\n\t	vsubpd	%%ymm14,%%ymm9 ,%%ymm9 		\n\t"\
"prefetcht1	%c[__pfetch_dist](%%r13,%%r14)\n\t"/*...+p11 */\
		"												\n\t	vmulpd	(%%r8) ,%%ymm12,%%ymm12		\n\t"\
		"												\n\t	vmulpd	(%%r8) ,%%ymm15,%%ymm15		\n\t"\
		"												\n\t	vmulpd	(%%r8) ,%%ymm13,%%ymm13		\n\t"\
		"												\n\t	vmulpd	(%%r8) ,%%ymm14,%%ymm14		\n\t"\
		"												\n\t	vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"												\n\t	vmovaps	%%ymm11,     (%%r12)	\n\t"\
		"												\n\t	vmovaps	%%ymm10,0x020(%%r11)	\n\t"\
		"												\n\t	vmovaps	%%ymm9 ,0x020(%%r13)	\n\t"\
		"/* Block 2: */									\n\t	vaddpd	%%ymm8 ,%%ymm12,%%ymm12	\n\t"\
		"addq	%%rdi,%%rax								\n\t	vaddpd	%%ymm11,%%ymm15,%%ymm15	\n\t"\
		"addq	%%rdi,%%rbx								\n\t	vaddpd	%%ymm10,%%ymm13,%%ymm13	\n\t"\
		"addq	%%rdi,%%rcx								\n\t	vaddpd	%%ymm9 ,%%ymm14,%%ymm14	\n\t"\
		"addq	%%rdi,%%rdx								\n\t	vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"addq	$0x040,%%r9	/* r3 */					\n\t	vmovaps	%%ymm15,     (%%r13)	\n\t"\
		"vmovaps	0x200(%%r9),%%ymm4					\n\t	vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"vmovaps	0x220(%%r9),%%ymm5					\n\t	vmovaps	%%ymm14,0x020(%%r12)	\n\t"\
		"/* Share cc0/ss0 between 2 halves: */			\n\t	/* Block 4: */				\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm11	/* cc0 */		\n\t	addq	%%rdi,%%r10			\n\t"\
		"vmovaps	0x040(%%rsi),%%ymm10	/* ss0 */		\n\t	addq	%%rdi,%%r11			\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r10)\n\t"/* ...+p12 */\
		"vmovaps	%%ymm4,%%ymm6						\n\t	addq	%%rdi,%%r12			\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	addq	%%rdi,%%r13			\n\t"\
		"vmulpd	%%ymm11,%%ymm4,%%ymm4					\n\t	vmovaps	0x280(%%r9),%%ymm12	/* r7 */\n\t"\
		"vmulpd	%%ymm10,%%ymm7,%%ymm7					\n\t	vmovaps	0x2a0(%%r9),%%ymm13	\n\t"\
		"vmulpd	%%ymm11,%%ymm5,%%ymm5					\n\t	vmovaps	%%ymm12,%%ymm14		\n\t"\
		"vmulpd	%%ymm10,%%ymm6,%%ymm6					\n\t	vmovaps	%%ymm13,%%ymm15		\n\t"\
		"vmovaps	0x300(%%r9),%%ymm0					\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	0x320(%%r9),%%ymm1					\n\t	vmulpd	%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13		\n\t"\
		"vaddpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vmulpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t	vmovaps	0x380(%%r9),%%ymm8 	\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t	vmovaps	0x3a0(%%r9),%%ymm9 	\n\t"\
		"vmulpd	%%ymm10,%%ymm6,%%ymm6					\n\t	vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmulpd	%%ymm11,%%ymm1,%%ymm1					\n\t	vaddpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vmulpd	%%ymm10,%%ymm7,%%ymm7					\n\t	vmovaps	%%ymm8 ,%%ymm14		\n\t"\
		"vmulpd	%%ymm11,%%ymm0,%%ymm0					\n\t	vmovaps	%%ymm9 ,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t	vmulpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t	vmulpd	%%ymm10,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm1 ,%%ymm6,%%ymm6					\n\t	vmulpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm0 ,%%ymm7,%%ymm7					\n\t	vmulpd	%%ymm11,%%ymm15,%%ymm15		\n\t"\
"prefetcht1	%c[__pfetch_dist](%%r11)\n\t"/* ...+p13 */\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4					\n\t	vsubpd	%%ymm9 ,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5					\n\t	vaddpd	%%ymm8 ,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm2,%%ymm6,%%ymm6					\n\t	vmovaps	%%ymm12,%%ymm10		\n\t"\
		"vaddpd	%%ymm3,%%ymm7,%%ymm7					\n\t	vmovaps	%%ymm13,%%ymm11		\n\t"\
		"vmovaps	0x100(%%r9),%%ymm2					\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	0x120(%%r9),%%ymm3					\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	(%%rsi),%%ymm1	/* isrt2 */			\n\t	vaddpd	%%ymm10,%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm2,%%ymm0						\n\t	vaddpd	%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm3,%%ymm2,%%ymm2					\n\t	vmovaps	0x180(%%r9),%%ymm10	\n\t"\
		"vaddpd	%%ymm0,%%ymm3,%%ymm3					\n\t	vmovaps	0x1a0(%%r9),%%ymm11	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2					\n\t	vmovaps	(%%rsi),%%ymm9 	/* isrt2 */\n\t"\
		"vmulpd	%%ymm1,%%ymm3,%%ymm3					\n\t	vmovaps	%%ymm10,%%ymm8 		\n\t"\
		"vmovaps	     (%%r9),%%ymm0					\n\t	vaddpd	%%ymm11,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	0x020(%%r9),%%ymm1					\n\t	vsubpd	%%ymm8 ,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t	vmulpd	%%ymm9 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2					\n\t	vmovaps	0x080(%%r9),%%ymm8 	\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3					\n\t	vmovaps	0x0a0(%%r9),%%ymm9 	\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r12)\n\t"/* ...+p14 */\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2					\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3					\n\t	vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm5,%%ymm0,%%ymm0					\n\t	vmulpd	(%%r8) ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1					\n\t	vmulpd	(%%r8) ,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm6,%%ymm2,%%ymm2					\n\t	vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm7,%%ymm3,%%ymm3					\n\t	vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t	vsubpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t	vsubpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t	vsubpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t	vsubpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm2,     (%%rbx)					\n\t	vmulpd	(%%r8) ,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm0,     (%%rcx)					\n\t	vmulpd	(%%r8) ,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm3,0x020(%%rbx)					\n\t	vmulpd	(%%r8) ,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm1,0x020(%%rdx)					\n\t	vmulpd	(%%r8) ,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm2,	%%ymm6,	%%ymm6					\n\t	vmovaps	%%ymm8 ,     (%%r11)	\n\t"\
		"vaddpd	%%ymm0,	%%ymm5,	%%ymm5					\n\t	vmovaps	%%ymm10,     (%%r12)	\n\t"\
		"vaddpd	%%ymm3,	%%ymm7,	%%ymm7					\n\t	vmovaps	%%ymm9 ,0x020(%%r11)	\n\t"\
		"vaddpd	%%ymm1,	%%ymm4,	%%ymm4					\n\t	vmovaps	%%ymm11,0x020(%%r13)	\n\t"\
		"vmovaps	%%ymm6,     (%%rax)					\n\t	vaddpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm5,     (%%rdx)					\n\t	vaddpd	%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm7,0x020(%%rax)					\n\t	vaddpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm4,0x020(%%rcx)					\n\t	vaddpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r13)\n\t"/* ...+p15 */\
		"														vmovaps	%%ymm12,     (%%r10)	\n\t"\
		"														vmovaps	%%ymm15,     (%%r13)	\n\t"\
		"														vmovaps	%%ymm13,0x020(%%r10)	\n\t"\
		"														vmovaps	%%ymm14,0x020(%%r12)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r1] "m" (Xr1)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_DIT_TWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xr1,Xisrt2,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
		"movq	%[__isrt2],%%r8					\n\t"\
		"movq	%[__add0],%%rax				\n\t"\
		"movslq	%[__p1],%%rbx					\n\t"\
		"movslq	%[__p2],%%rcx					\n\t"\
		"movslq	%[__p3],%%rdx					\n\t"\
		"movslq	%[__p4],%%rdi					\n\t"\
		"shlq	$3,%%rbx						\n\t"\
		"shlq	$3,%%rcx						\n\t"\
		"shlq	$3,%%rdx						\n\t"\
		"shlq	$3,%%rdi						\n\t"\
		"addq	$0x460,%%r8	/* two */			\n\t"\
		"addq	%%rax,%%rbx						\n\t"\
		"addq	%%rax,%%rcx						\n\t"\
		"addq	%%rax,%%rdx						\n\t		movslq	%[__p8],%%r10		\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r1) */	\n\t		shlq	$3,%%r10			\n\t"\
		"movq	%[__r1],%%rsi					\n\t		movq	%%r10,%%r11			\n\t"\
		"vmovaps	     (%%rax),%%ymm2			\n\t		movq	%%r10,%%r12			\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		movq	%%r10,%%r13			\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		addq	%%rax,%%r10			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		addq	%%rbx,%%r11			\n\t"\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		addq	%%rcx,%%r12			\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		addq	%%rdx,%%r13			\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r17) */\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		leaq	0x200(%%rsi),%%r14	/* r17 */\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10	\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14	\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11	\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t		vmovaps	0x020(%%r12),%%ymm15	\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4			\n\t		vmovaps	     (%%r13),%%ymm12	\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5			\n\t		vmovaps	0x020(%%r13),%%ymm13	\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0			\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4			\n\t		vsubpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1			\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5			\n\t		vsubpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t		vmulpd	 (%%r8),%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t		vmulpd	 (%%r8),%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t		vmulpd	 (%%r8),%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t		vmulpd	 (%%r8),%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)			\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm2,0x0c0(%%rsi)			\n\t		vaddpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)			\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)			\n\t		vaddpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4			\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7			\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5			\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6			\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4			\n\t		vmovaps	%%ymm8 ,0x080(%%r14)\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7			\n\t		vmovaps	%%ymm10,0x0c0(%%r14)\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5			\n\t		vmovaps	%%ymm9 ,0x0a0(%%r14)\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6			\n\t		vmovaps	%%ymm11,0x060(%%r14)\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)			\n\t		vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm7,0x040(%%rsi)			\n\t		vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)			\n\t		vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)			\n\t		vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"addq	%%rdi,%%rax						\n\t		vaddpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"addq	%%rdi,%%rbx						\n\t		vaddpd	%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"addq	%%rdi,%%rcx						\n\t		vaddpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"addq	%%rdi,%%rdx						\n\t		vaddpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r9) */	\n\t		vmovaps	%%ymm12,     (%%r14)\n\t"\
		"addq	$0x100,%%rsi		/* r9 */	\n\t		vmovaps	%%ymm15,0x040(%%r14)\n\t"\
		"vmovaps	     (%%rax),%%ymm2			\n\t		vmovaps	%%ymm13,0x020(%%r14)\n\t"\
		"vmovaps	     (%%rcx),%%ymm6			\n\t		vmovaps	%%ymm14,0x0e0(%%r14)\n\t"\
		"vmovaps	0x020(%%rax),%%ymm3			\n\t		addq	%%rdi,%%r10			\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm7			\n\t		addq	%%rdi,%%r11			\n\t"\
		"vmovaps	     (%%rbx),%%ymm0			\n\t		addq	%%rdi,%%r12			\n\t"\
		"vmovaps	     (%%rdx),%%ymm4			\n\t		addq	%%rdi,%%r13			\n\t"\
		"vmovaps	0x020(%%rbx),%%ymm1			\n\t		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r25) */\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm5			\n\t		leaq	0x200(%%rsi),%%r14	/* r25 */\n\t"\
		"vsubpd	%%ymm0,%%ymm2,%%ymm2			\n\t		vmovaps	     (%%r10),%%ymm10	\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6			\n\t		vmovaps	     (%%r12),%%ymm14	\n\t"\
		"vsubpd	%%ymm1,%%ymm3,%%ymm3			\n\t		vmovaps	0x020(%%r10),%%ymm11	\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7			\n\t		vmovaps	0x020(%%r12),%%ymm15	\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0			\n\t		vmovaps	     (%%r11),%%ymm8 	\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4			\n\t		vmovaps	     (%%r13),%%ymm12	\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1			\n\t		vmovaps	0x020(%%r11),%%ymm9 	\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5			\n\t		vmovaps	0x020(%%r13),%%ymm13	\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0			\n\t		vsubpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4			\n\t		vsubpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1			\n\t		vsubpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5			\n\t		vsubpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm4,%%ymm0,%%ymm0			\n\t		vmulpd	 (%%r8),%%ymm8 ,%%ymm8 		\n\t"\
		"vsubpd	%%ymm7,%%ymm2,%%ymm2			\n\t		vmulpd	 (%%r8),%%ymm12,%%ymm12		\n\t"\
		"vsubpd	%%ymm5,%%ymm1,%%ymm1			\n\t		vmulpd	 (%%r8),%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	%%ymm6,%%ymm3,%%ymm3			\n\t		vmulpd	 (%%r8),%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rsi)			\n\t		vaddpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm2,0x0c0(%%rsi)			\n\t		vaddpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rsi)			\n\t		vaddpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm3,0x060(%%rsi)			\n\t		vaddpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4			\n\t		vsubpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7			\n\t		vsubpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5			\n\t		vsubpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6			\n\t		vsubpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4			\n\t		vmovaps	%%ymm8 ,0x080(%%r14)\n\t"\
		"vaddpd	%%ymm2,%%ymm7,%%ymm7			\n\t		vmovaps	%%ymm10,0x0c0(%%r14)\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5			\n\t		vmovaps	%%ymm9 ,0x0a0(%%r14)\n\t"\
		"vaddpd	%%ymm3,%%ymm6,%%ymm6			\n\t		vmovaps	%%ymm11,0x060(%%r14)\n\t"\
		"vmovaps	%%ymm4,     (%%rsi)			\n\t		vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm7,0x040(%%rsi)			\n\t		vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rsi)			\n\t		vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm6,0x0e0(%%rsi)			\n\t		vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"													vaddpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"													vaddpd	%%ymm10,%%ymm15,%%ymm15		\n\t"\
		"													vaddpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"													vaddpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
		"													vmovaps	%%ymm12,     (%%r14)\n\t"\
		"													vmovaps	%%ymm15,0x040(%%r14)\n\t"\
		"													vmovaps	%%ymm13,0x020(%%r14)\n\t"\
		"													vmovaps	%%ymm14,0x0e0(%%r14)\n\t"\
	/***************************************************************************************************/\
	/*..and now do four more radix-4 transforms, including the internal and external twiddle factors:  */\
	/***************************************************************************************************/\
	/* May 2016: With genuine-BR-ordered roots (to match DIF data layout, have
	SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c0,2,1,3) base-address offsets = r25+0xb0,isr2+0x130,0xb0,0x1b0,
	and base-plus-* byte offsets in the 2 cmul blocks changed from 0x[00,40],[80,c0] to 0x[00,80],[40,c0].
	We also get rid of the sincos-address-incrementing previously used in between the 2 CMUL-pairs done
	in the 2ND_HALF_B portion of each posttwiddle-4-DFT, so w.r.to the address offsets in said 2 pairs this means
	address mungings:
		Set 1: 0246 -> 028a
		Set 2: 0246 -> 46ce .
	*/\
	/*...Block 1: 24 MOVapd, 26 ADD/SUBpd, 12 MULpd */		/*...Block 3: 31 MOVapd, 32 ADD/SUBpd, 20 MULpd */\
		"movq	%[__add0],%%rax							\n\t	movslq	%[__p2],%%r10			\n\t"\
		"movslq	%[__p4],%%rbx							\n\t	leaq	-0x080(%%rsi),%%r12	/* r5  */\n\t"\
		"movslq	%[__p8],%%rdi							\n\t	leaq	 0x080(%%rsi),%%r13	/* r13 */\n\t"\
		"shlq	$3,%%rbx								\n\t	movq	%%rbx,%%r11		/* p4 */\n\t"\
		"shlq	$3,%%rdi								\n\t	shlq	$3,%%r10		/* p2 */\n\t"\
		"movq	%[__r1],%%rcx							\n\t	vmovaps	-0x460(%%r8),%%ymm10	/* isrt2 */\n\t"\
		"movq	%%rsi,%%rdx		/* r9 */				\n\t	vmovaps	0x200(%%r12),%%ymm12	\n\t"\
		"addq	%%rax,%%r10	/* a[j+p2 ] */				\n\t	vmovaps	0x220(%%r12),%%ymm13	\n\t"\
		"addq	%%rax,%%rbx	/* a[j+p4 ] */				\n\t	vmovaps	0x200(%%r13),%%ymm14	\n\t"\
		"addq	%%r10,%%r11	/* a[j+p6 ] */				\n\t	vmovaps	0x220(%%r13),%%ymm15	\n\t"\
		"vmovaps	     (%%rdx),%%ymm2					\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	0x200(%%rdx),%%ymm4					\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3					\n\t	vmulpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	0x220(%%rdx),%%ymm5					\n\t	vmulpd	%%ymm10,%%ymm15,%%ymm15			\n\t"\
		"vmovaps	     (%%rcx),%%ymm0					\n\t	vmovaps	     (%%r12),%%ymm8 		\n\t"\
		"vmovaps	0x200(%%rcx),%%ymm6					\n\t	vmovaps	0x020(%%r13),%%ymm10	\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t	vmovaps	0x020(%%r12),%%ymm11	\n\t"\
		"vmovaps	0x220(%%rcx),%%ymm7					\n\t	vmovaps	     (%%r13),%%ymm9 		\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t	vsubpd	%%ymm12,%%ymm13,%%ymm13			\n\t"\
		"vsubpd	%%ymm4,%%ymm6,%%ymm6					\n\t	vsubpd	%%ymm15,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12			\n\t"\
		"vsubpd	%%ymm5,%%ymm7,%%ymm7					\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8),%%ymm2,%%ymm2					\n\t	vaddpd	%%ymm13,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t	vaddpd	%%ymm14,%%ymm15,%%ymm15			\n\t"\
		"vmulpd	(%%r8),%%ymm3,%%ymm3					\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t	vsubpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2					\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm6,%%ymm4,%%ymm4					\n\t	vsubpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3					\n\t	vmulpd	 (%%r8),%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm7,%%ymm5,%%ymm5					\n\t	vmulpd	 (%%r8),%%ymm10,%%ymm10			\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c0) */	\n\t	vmulpd	 (%%r8),%%ymm15,%%ymm15			\n\t"\
		"addq	$0x360,%%rsi	/* c0, from r9 */		\n\t	vmulpd	 (%%r8),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t	vaddpd	%%ymm12,%%ymm14,%%ymm14			\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0					\n\t	vaddpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t	vaddpd	%%ymm13,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1					\n\t	vaddpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t	/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c2) */\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t	addq	$0x260,%%r14	/* c2, from r25 */\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10			\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t	vsubpd	%%ymm15,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4					\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11			\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7					\n\t	vsubpd	%%ymm14,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5					\n\t	vmulpd	 (%%r8),%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6					\n\t	vmulpd	 (%%r8),%%ymm15,%%ymm15			\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)					\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)					\n\t	vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)					\n\t	vaddpd	%%ymm10,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)					\n\t	vaddpd	%%ymm8 ,%%ymm15,%%ymm15			\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t	vaddpd	%%ymm11,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	%%ymm7,%%ymm0						\n\t	vaddpd	%%ymm9 ,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t	vmovaps	%%ymm10,     (%%r12)	\n\t"\
		"vmovaps	%%ymm1,%%ymm6						\n\t	vmovaps	%%ymm8 ,0x020(%%r13)	\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5				\n\t	vmovaps	%%ymm11,0x020(%%r12)	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2				\n\t	vmovaps	%%ymm14,     (%%r13)	\n\t"\
		"vmulpd	0x0a0(%%rsi),%%ymm0,%%ymm0				\n\t	vmovaps	%%ymm12,%%ymm10			\n\t"\
		"vmulpd	0x080(%%rsi),%%ymm1,%%ymm1				\n\t	vmovaps	%%ymm15,%%ymm8 			\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3				\n\t	vmovaps	%%ymm13,%%ymm11			\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4				\n\t	vmovaps	%%ymm9 ,%%ymm14			\n\t"\
		"vmulpd	0x080(%%rsi),%%ymm7,%%ymm7				\n\t	vmulpd	0x020(%%r14),%%ymm10,%%ymm10		\n\t"\
		"vmulpd	0x0a0(%%rsi),%%ymm6,%%ymm6				\n\t	vmulpd	     (%%r14),%%ymm13,%%ymm13		\n\t"\
		"vsubpd	     %%ymm2,%%ymm5,%%ymm5				\n\t	vmulpd	0x080(%%r14),%%ymm9 ,%%ymm9 		\n\t"\
		"vsubpd	     %%ymm0,%%ymm1,%%ymm1				\n\t	vmulpd	0x0a0(%%r14),%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	     %%ymm3,%%ymm4,%%ymm4				\n\t	vmulpd	0x020(%%r14),%%ymm11,%%ymm11		\n\t"\
		"vaddpd	     %%ymm6,%%ymm7,%%ymm7				\n\t	vmulpd	     (%%r14),%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t	vmulpd	0x080(%%r14),%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)					\n\t	vmulpd	0x0a0(%%r14),%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t	vsubpd	%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)					\n\t	vsubpd	%%ymm8 ,%%ymm9 ,%%ymm9 			\n\t"\
		"addq	%%rdi,%%rax	/* a[j+p8 ] */				\n\t	vaddpd	%%ymm11,%%ymm12,%%ymm12			\n\t"\
		"addq	%%rdi,%%rbx	/* a[j+p12] */				\n\t	vaddpd	%%ymm14,%%ymm15,%%ymm15			\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	vmovaps	%%ymm13,0x020(%%r10)		\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t	vmovaps	%%ymm9 ,0x020(%%r11)		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	%%ymm12,     (%%r10)		\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	vmovaps	%%ymm15,     (%%r11)		\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t	addq	%%rdi,%%r10	/* a[j+p10] */	\n\t"\
		"vmovaps	%%ymm0,%%ymm1						\n\t	addq	%%rdi,%%r11	/* a[j+p14] */	\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmovaps	%%ymm6,%%ymm7						\n\t	vmovaps	0x020(%%r13),%%ymm8 		\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm5,%%ymm5				\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm2,%%ymm2				\n\t	vmovaps	     (%%r13),%%ymm14	\n\t"\
		"vmulpd	0x0c0(%%rsi),%%ymm6,%%ymm6				\n\t	vmovaps	%%ymm12,%%ymm10			\n\t"\
		"vmulpd	0x0e0(%%rsi),%%ymm1,%%ymm1				\n\t	vmovaps	%%ymm8 ,%%ymm9 			\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm3,%%ymm3				\n\t	vmovaps	%%ymm13,%%ymm11			\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm4,%%ymm4				\n\t	vmovaps	%%ymm14,%%ymm15			\n\t"\
		"vmulpd	0x0c0(%%rsi),%%ymm0,%%ymm0				\n\t	vmulpd	0x060(%%r14),%%ymm10,%%ymm10		\n\t"\
		"vmulpd	0x0e0(%%rsi),%%ymm7,%%ymm7				\n\t	vmulpd	0x040(%%r14),%%ymm13,%%ymm13		\n\t"\
		"vsubpd	     %%ymm2,%%ymm5,%%ymm5				\n\t	vmulpd	0x0c0(%%r14),%%ymm14,%%ymm14		\n\t"\
		"vsubpd	     %%ymm1,%%ymm6,%%ymm6				\n\t	vmulpd	0x0e0(%%r14),%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	     %%ymm3,%%ymm4,%%ymm4				\n\t	vmulpd	0x060(%%r14),%%ymm11,%%ymm11		\n\t"\
		"vaddpd	     %%ymm7,%%ymm0,%%ymm0				\n\t	vmulpd	0x040(%%r14),%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t	vmulpd	0x0c0(%%r14),%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)					\n\t	vmulpd	0x0e0(%%r14),%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t	vsubpd	%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)					\n\t	vsubpd	%%ymm9 ,%%ymm14,%%ymm14			\n\t"\
		/*...Block 2: 34 MOVapd, 36 ADD/SUB, 26 MUL */		"	vaddpd	%%ymm11,%%ymm12,%%ymm12			\n\t"\
		"movslq	%[__p1],%%r14							\n\t	vaddpd	%%ymm15,%%ymm8 ,%%ymm8 			\n\t"\
		"subq	%%rdi,%%rax	/* a[j+p0 ] */				\n\t	vmovaps	%%ymm13,0x020(%%r10)		\n\t"\
		"subq	%%rdi,%%rbx	/* a[j+p4 ] */				\n\t	vmovaps	%%ymm14,0x020(%%r11)		\n\t"\
		"shlq	$3,%%r14								\n\t	vmovaps	%%ymm12,     (%%r10)		\n\t"\
		"addq	%%r14,%%rax	/* a[j+p1 ] */				\n\t	vmovaps	%%ymm8 ,     (%%r11)		\n\t"\
		"addq	%%r14,%%rbx	/* a[j+p5 ] */				\n\t"	/*...Block 4: 34 MOVapd, 36 ADD/SUBpd, 26 MULpd */\
		"addq	$0x040,%%rcx		/* r3 , from r1 */	\n\t	subq	%%rdi,%%r10	/* a[j+p2 ] */	\n\t"\
		"leaq	0x100(%%rcx),%%rdx	/* r11, from r3 */	\n\t	subq	%%rdi,%%r11	/* a[j+p6 ] */	\n\t"\
		"subq	$0x040,%%rsi		/* cc0, from c0 */	\n\t	addq	%%r14,%%r10	/* a[j+p3 ] */	\n\t"\
		"vmovaps	0x200(%%rcx),%%ymm4					\n\t	addq	%%r14,%%r11	/* a[j+p7 ] */	\n\t"\
		"vmovaps	0x220(%%rcx),%%ymm5					\n\t	addq	$0x040,%%r12	/* r7 , from r5  */\n\t"\
		"vmovaps	     (%%rsi),%%ymm2					\n\t	addq	$0x040,%%r13	/* r15, from r13 */\n\t"\
		"vmovaps	0x020(%%rsi),%%ymm3					\n\t	addq	$0x200,%%r12			\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	addq	$0x200,%%r13			\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t	vmovaps	     (%%rsi),%%ymm11	\n\t"\
		"vmulpd	%%ymm3,%%ymm6,%%ymm6					\n\t	vmovaps	0x020(%%rsi),%%ymm10	\n\t"\
		"vmulpd	%%ymm3,%%ymm7,%%ymm7					\n\t	vmovaps	%%ymm12,%%ymm14			\n\t"\
		"vmovaps	0x200(%%rdx),%%ymm0					\n\t	vmovaps	%%ymm13,%%ymm15			\n\t"\
		"vmovaps	0x220(%%rdx),%%ymm1					\n\t	vmulpd	%%ymm10,%%ymm12,%%ymm12			\n\t"\
		"vaddpd	%%ymm7,%%ymm4,%%ymm4					\n\t	vmulpd	%%ymm11,%%ymm15,%%ymm15			\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5					\n\t	vmulpd	%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	%%ymm0,%%ymm6						\n\t	vmulpd	%%ymm11,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	%%ymm1,%%ymm7						\n\t	vmovaps	     (%%r13),%%ymm8 		\n\t"\
		"vmulpd	%%ymm3,%%ymm0,%%ymm0					\n\t	vmovaps	0x020(%%r13),%%ymm9 		\n\t"\
		"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t	vaddpd	%%ymm15,%%ymm12,%%ymm12			\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm14,%%ymm13,%%ymm13			\n\t"\
		"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t	vmovaps	%%ymm8 ,%%ymm14			\n\t"\
		"vaddpd	%%ymm7,%%ymm0,%%ymm0					\n\t	vmovaps	%%ymm9 ,%%ymm15			\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1					\n\t	vmulpd	%%ymm11,%%ymm8 ,%%ymm8 			\n\t"\
		"vmovaps	%%ymm5,%%ymm7						\n\t	vmulpd	%%ymm10,%%ymm15,%%ymm15			\n\t"\
		"vmovaps	%%ymm4,%%ymm6						\n\t	vmulpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4					\n\t	vmulpd	%%ymm10,%%ymm14,%%ymm14			\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5					\n\t	vaddpd	%%ymm15,%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm0,%%ymm6,%%ymm6					\n\t	vsubpd	%%ymm14,%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm1,%%ymm7,%%ymm7					\n\t	vmovaps	%%ymm13,%%ymm15			\n\t"\
		"vmovaps	     (%%rdx),%%ymm2					\n\t	vmovaps	%%ymm12,%%ymm14			\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm3					\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	-0x460(%%r8),%%ymm1					\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15			\n\t"\
		"vmovaps	%%ymm3,%%ymm0						\n\t	vsubpd	%%ymm8 ,%%ymm12,%%ymm12			\n\t"\
		"vsubpd	%%ymm2,%%ymm3,%%ymm3					\n\t	vsubpd	%%ymm9 ,%%ymm13,%%ymm13			\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2					\n\t	subq	$0x200,%%r12			\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2					\n\t	subq	$0x200,%%r13			\n\t"\
		"vmulpd	%%ymm1,%%ymm3,%%ymm3					\n\t	vmovaps	     (%%r13),%%ymm8 		\n\t"\
		"vmovaps	     (%%rcx),%%ymm0					\n\t	vmovaps	0x020(%%r13),%%ymm9 		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm1					\n\t	vmovaps	-0x460(%%r8),%%ymm11	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0					\n\t	vmovaps	%%ymm8 ,%%ymm10			\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1					\n\t	vsubpd	%%ymm9 ,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd	(%%r8),%%ymm2,%%ymm2					\n\t	vaddpd	%%ymm10,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd	(%%r8),%%ymm3,%%ymm3					\n\t	vmulpd	%%ymm11,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2					\n\t	vmulpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3					\n\t	vmovaps	     (%%r12),%%ymm10	\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c1) */	\n\t	vmovaps	0x020(%%r12),%%ymm11	\n\t"\
		"addq	$0x240,%%rsi	/* c1, from cc0 */		\n\t	vsubpd	%%ymm8 ,%%ymm10,%%ymm10			\n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2					\n\t	vsubpd	%%ymm9 ,%%ymm11,%%ymm11			\n\t"\
		"vsubpd	%%ymm7,%%ymm0,%%ymm0					\n\t	vmulpd	 (%%r8),%%ymm8 ,%%ymm8 			\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3					\n\t	vmulpd	 (%%r8),%%ymm9 ,%%ymm9 			\n\t"\
		"vsubpd	%%ymm6,%%ymm1,%%ymm1					\n\t	vaddpd	%%ymm10,%%ymm8 ,%%ymm8 			\n\t"\
		"vmulpd	(%%r8),%%ymm4,%%ymm4					\n\t	vaddpd	%%ymm11,%%ymm9 ,%%ymm9 			\n\t"\
		"vmulpd	(%%r8),%%ymm7,%%ymm7					\n\t/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c3) */\n\t"\
		"vmulpd	(%%r8),%%ymm5,%%ymm5					\n\t	leaq	0x100(%%rsi),%%r14	/* c3, from c1 */\n\t"\
		"vmulpd	(%%r8),%%ymm6,%%ymm6					\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10			\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4					\n\t	vsubpd	%%ymm15,%%ymm8 ,%%ymm8 			\n\t"\
		"vaddpd	%%ymm0,%%ymm7,%%ymm7					\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11			\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5					\n\t	vsubpd	%%ymm14,%%ymm9 ,%%ymm9 			\n\t"\
		"vaddpd	%%ymm1,%%ymm6,%%ymm6					\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	%%ymm2,     (%%rcx)					\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15			\n\t"\
		"vmovaps	%%ymm0,0x020(%%rdx)					\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	%%ymm3,0x020(%%rcx)					\n\t	vaddpd	%%ymm14,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	%%ymm6,     (%%rdx)					\n\t	vaddpd	%%ymm10,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t	vaddpd	%%ymm8 ,%%ymm15,%%ymm15			\n\t"\
		"vmovaps	%%ymm7,%%ymm0						\n\t	vaddpd	%%ymm11,%%ymm13,%%ymm13			\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t	vaddpd	%%ymm9 ,%%ymm14,%%ymm14			\n\t"\
		"vmovaps	%%ymm1,%%ymm6						\n\t	vmovaps	%%ymm10,     (%%r12)	\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm2,%%ymm2				\n\t	vmovaps	%%ymm8 ,0x020(%%r13)	\n\t"\
		"vmulpd	     (%%rsi),%%ymm5,%%ymm5				\n\t	vmovaps	%%ymm11,0x020(%%r12)	\n\t"\
		"vmulpd	0x080(%%rsi),%%ymm1,%%ymm1				\n\t	vmovaps	%%ymm14,     (%%r13)	\n\t"\
		"vmulpd	0x0a0(%%rsi),%%ymm0,%%ymm0				\n\t	vmovaps	%%ymm12,%%ymm10			\n\t"\
		"vmulpd	0x020(%%rsi),%%ymm3,%%ymm3				\n\t	vmovaps	%%ymm15,%%ymm8 			\n\t"\
		"vmulpd	     (%%rsi),%%ymm4,%%ymm4				\n\t	vmovaps	%%ymm13,%%ymm11			\n\t"\
		"vmulpd	0x080(%%rsi),%%ymm7,%%ymm7				\n\t	vmovaps	%%ymm9 ,%%ymm14			\n\t"\
		"vmulpd	0x0a0(%%rsi),%%ymm6,%%ymm6				\n\t	vmulpd	0x020(%%r14),%%ymm10,%%ymm10		\n\t"\
		"vsubpd	     %%ymm2,%%ymm5,%%ymm5				\n\t	vmulpd	     (%%r14),%%ymm13,%%ymm13		\n\t"\
		"vsubpd	     %%ymm0,%%ymm1,%%ymm1				\n\t	vmulpd	0x080(%%r14),%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	     %%ymm3,%%ymm4,%%ymm4				\n\t	vmulpd	0x0a0(%%r14),%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	     %%ymm6,%%ymm7,%%ymm7				\n\t	vmulpd	0x020(%%r14),%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t	vmulpd	     (%%r14),%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm1,0x020(%%rbx)					\n\t	vmulpd	0x080(%%r14),%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t	vmulpd	0x0a0(%%r14),%%ymm14,%%ymm14		\n\t"\
		"vmovaps	%%ymm7,     (%%rbx)					\n\t	vsubpd	%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"addq	%%rdi,%%rax								\n\t	vsubpd	%%ymm8 ,%%ymm9 ,%%ymm9 			\n\t"\
		"addq	%%rdi,%%rbx								\n\t	vaddpd	%%ymm11,%%ymm12,%%ymm12			\n\t"\
		"vmovaps	     (%%rcx),%%ymm4					\n\t	vaddpd	%%ymm14,%%ymm15,%%ymm15			\n\t"\
		"vmovaps	0x020(%%rdx),%%ymm0					\n\t	vmovaps	%%ymm13,0x020(%%r10)		\n\t"\
		"vmovaps	0x020(%%rcx),%%ymm5					\n\t	vmovaps	%%ymm9 ,0x020(%%r11)		\n\t"\
		"vmovaps	     (%%rdx),%%ymm6					\n\t	vmovaps	%%ymm12,     (%%r10)		\n\t"\
		"vmovaps	%%ymm4,%%ymm2						\n\t	vmovaps	%%ymm15,     (%%r11)		\n\t"\
		"vmovaps	%%ymm0,%%ymm1						\n\t	addq	%%rdi,%%r10				\n\t"\
		"vmovaps	%%ymm5,%%ymm3						\n\t	addq	%%rdi,%%r11				\n\t"\
		"vmovaps	%%ymm6,%%ymm7						\n\t	vmovaps	     (%%r12),%%ymm12	\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm2,%%ymm2				\n\t	vmovaps	0x020(%%r13),%%ymm8 		\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm5,%%ymm5				\n\t	vmovaps	0x020(%%r12),%%ymm13	\n\t"\
		"vmulpd	0x0c0(%%rsi),%%ymm6,%%ymm6				\n\t	vmovaps	     (%%r13),%%ymm14	\n\t"\
		"vmulpd	0x0e0(%%rsi),%%ymm1,%%ymm1				\n\t	vmovaps	%%ymm12,%%ymm10			\n\t"\
		"vmulpd	0x060(%%rsi),%%ymm3,%%ymm3				\n\t	vmovaps	%%ymm8 ,%%ymm9 			\n\t"\
		"vmulpd	0x040(%%rsi),%%ymm4,%%ymm4				\n\t	vmovaps	%%ymm13,%%ymm11			\n\t"\
		"vmulpd	0x0c0(%%rsi),%%ymm0,%%ymm0				\n\t	vmovaps	%%ymm14,%%ymm15			\n\t"\
		"vmulpd	0x0e0(%%rsi),%%ymm7,%%ymm7				\n\t	vmulpd	0x060(%%r14),%%ymm10,%%ymm10		\n\t"\
		"vsubpd	     %%ymm2,%%ymm5,%%ymm5				\n\t	vmulpd	0x040(%%r14),%%ymm13,%%ymm13		\n\t"\
		"vsubpd	     %%ymm1,%%ymm6,%%ymm6				\n\t	vmulpd	0x0c0(%%r14),%%ymm14,%%ymm14		\n\t"\
		"vaddpd	     %%ymm3,%%ymm4,%%ymm4				\n\t	vmulpd	0x0e0(%%r14),%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	     %%ymm7,%%ymm0,%%ymm0				\n\t	vmulpd	0x060(%%r14),%%ymm11,%%ymm11		\n\t"\
		"vmovaps	%%ymm5,0x020(%%rax)					\n\t	vmulpd	0x040(%%r14),%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm6,0x020(%%rbx)					\n\t	vmulpd	0x0c0(%%r14),%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm4,     (%%rax)					\n\t	vmulpd	0x0e0(%%r14),%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm0,     (%%rbx)					\n\t	vsubpd	%%ymm10,%%ymm13,%%ymm13			\n\t"\
		"												\n\t	vsubpd	%%ymm9 ,%%ymm14,%%ymm14			\n\t"\
		"												\n\t	vaddpd	%%ymm11,%%ymm12,%%ymm12			\n\t"\
		"												\n\t	vaddpd	%%ymm15,%%ymm8 ,%%ymm8 			\n\t"\
		"												\n\t	vmovaps	%%ymm13,0x020(%%r10)		\n\t"\
		"												\n\t	vmovaps	%%ymm14,0x020(%%r11)		\n\t"\
		"												\n\t	vmovaps	%%ymm12,     (%%r10)		\n\t"\
		"												\n\t	vmovaps	%%ymm8 ,     (%%r11)		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__r1] "m" (Xr1)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_SSE2) && (OS_BITS == 64)

	/* Do 4 prefetches from main-array address offsets ...+p0,4,8,c in this macro, and want to spread them out
	roughly evenly, which works out to 1 every 100 lines or so. Here is the i/o register usage summary for the
	four distinct SSE2_RADIX_04_DIF_3TWIDDLE_X1,2 - corresponding sub-blocks [1-4] and the actual address-offset
	we end up using in both a 16-prefetch (as in SSE2_RADIX16_DIF_TWIDDLE(), but that seems slower in _V2 here)
	and 4-prefetch paradigm:
									prefetch offsets:
								 possible:				actual:
												16-fetch	4-fetch
	[1]	[a-d]x: 0,4,8,c	rsi: 2	0,2,4,6,8,a,c,e	0,4,8,c		0,8 [one at start of block, other at end]
	[2]	[a-d]x: 1,5,9,d	rsi: 2	1,3,5,7,9,b,d,f	1,5,9,d
	[3]	r10-13: 0,2,1,3	r9 : 4	0-7				2,3,6,7		4
	[4]	r10-13: 8,a,9,b	r9 : 4	8-f				a,b,e,f		C
	*/
	#define SSE2_RADIX16_DIF_TWIDDLE_V2(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr0,Xtwo,Xcc0,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/*
		k1 = p2*8;	// Set k1 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k2 = 0x80;	// Set k2 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		// Pass1-DFTs all use same twiddle-triplet, so lit-byte address offset between ptrs to 1st,2nd DFT's twiddles = 0
		i0 = (vec_dbl *)add0; i1 = (vec_dbl *)(add0+p4); i2 = (vec_dbl *)(add0+p8); i3 = (vec_dbl *)(add0+p12);
		o0 = r1; o1 = r1+2; o2 = r1+4; o3 = r1+6;
		c_tmp = cc0;	// c8,4,C
		SSE2_RADIX_04_DIF_3TWIDDLE_X1(i0,i1,i2,i3,k1, two,c_tmp, o0,o1,o2,o3,k2)
	*/\
		"movq	%[__cc0],%%rsi 		\n\t	movslq	%[__p2],%%rdi		\n\t"\
		"movq	%[__add0],%%rax		\n\t	movq	%[__r0],%%r10		\n\t"\
		"movslq	%[__p4],%%rbx		\n\t	leaq	(%%rax,%%rbx,8),%%rbx	\n\t"\
		"movslq	%[__p8],%%rcx		\n\t	leaq	(%%rax,%%rcx,8),%%rcx	\n\t"\
		"movslq	%[__p12],%%rdx		\n\t	leaq	(%%rax,%%rdx,8),%%rdx	\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%rax)\n\t"/* [base-address + data-fetch-ahead index] + p0 */\
		/* Do	the p0,p2 combo: */\
		"movaps	    (%%rcx),%%xmm4	\n\t	shlq	$3,%%rdi			\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rsi),%%xmm10	\n\t"\
		"movaps	0x10(%%rsi),%%xmm11	\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t	movaps	    (%%rcx,%%rdi),%%xmm12	\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t	movaps	0x10(%%rcx,%%rdi),%%xmm13	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"mulpd	%%xmm10,%%xmm4		\n\t	movaps	%%xmm12,%%xmm14		\n\t"\
		"mulpd	%%xmm10,%%xmm5		\n\t	movaps	%%xmm13,%%xmm15		\n\t"\
		"mulpd	%%xmm11,%%xmm6		\n\t	movaps	    (%%rax,%%rdi),%%xmm8 	\n\t	movaps	%%xmm0,%%xmm2		\n\t"\
		"mulpd	%%xmm11,%%xmm7		\n\t	movaps	0x10(%%rax,%%rdi),%%xmm9 	\n\t	movaps	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t	mulpd	%%xmm10,%%xmm12		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t	mulpd	%%xmm10,%%xmm13		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t	mulpd	%%xmm11,%%xmm14		\n\t	movaps	%%xmm8 ,%%xmm10		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t	mulpd	%%xmm11,%%xmm15		\n\t	movaps	%%xmm9 ,%%xmm11		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t	addpd	%%xmm14,%%xmm13		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t	subpd	%%xmm15,%%xmm12		\n\t"\
		"movaps	%%xmm0,    (%%r10)	\n\t	addpd	%%xmm12,%%xmm8 		\n\t"/* Spill 1: free up xmm0,1 */\
		"movaps	%%xmm1,0x10(%%r10)	\n\t	addpd	%%xmm13,%%xmm9 		\n\t"/* Do	the p1,3 combo: */\
		"movaps	0x40(%%rsi),%%xmm0	\n\t	subpd	%%xmm12,%%xmm10		\n\t"\
		"movaps	0x50(%%rsi),%%xmm1	\n\t	subpd	%%xmm13,%%xmm11		\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t	movaps	%%xmm8 ,0x80(%%r10)	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t	movaps	%%xmm9 ,0x90(%%r10)	\n\t"\
		"movaps	%%xmm6,%%xmm4		\n\t"\
		"movaps	%%xmm7,%%xmm5		\n\t"\
		"mulpd	%%xmm0,%%xmm4		\n\t	movaps	    (%%rdx,%%rdi),%%xmm14	\n\t"\
		"mulpd	%%xmm0,%%xmm5		\n\t	movaps	0x10(%%rdx,%%rdi),%%xmm15	\n\t"\
		"mulpd	%%xmm1,%%xmm6		\n\t	movaps	%%xmm14,%%xmm12		\n\t"\
		"mulpd	%%xmm1,%%xmm7		\n\t	movaps	%%xmm15,%%xmm13		\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t	mulpd	%%xmm0 ,%%xmm12		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t	mulpd	%%xmm0 ,%%xmm13		\n\t"\
		"movaps	%%xmm5,0x50(%%r10)	\n\t	mulpd	%%xmm1 ,%%xmm14		\n\t"/* Spill 2 */\
		"movaps	%%xmm4,0x40(%%r10)	\n\t	mulpd	%%xmm1 ,%%xmm15		\n\t"\
		"movaps	0x20(%%rsi),%%xmm8	\n\t	addpd	%%xmm14,%%xmm13		\n\t"\
		"movaps	0x30(%%rsi),%%xmm9	\n\t	subpd	%%xmm15,%%xmm12		\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t	movaps	%%xmm13,0xd0(%%r10)	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t	movaps	%%xmm12,0xc0(%%r10)	\n\t"\
		"movaps	%%xmm6,%%xmm4		\n\t"\
		"movaps	%%xmm7,%%xmm5		\n\t"\
		"mulpd	%%xmm8,%%xmm4		\n\t	movaps	    (%%rbx,%%rdi),%%xmm14	\n\t"\
		"mulpd	%%xmm8,%%xmm5		\n\t	movaps	0x10(%%rbx,%%rdi),%%xmm15	\n\t"\
		"mulpd	%%xmm9,%%xmm6		\n\t	movaps	%%xmm14,%%xmm12		\n\t"\
		"mulpd	%%xmm9,%%xmm7		\n\t	movaps	%%xmm15,%%xmm13		\n\t"\
		"movaps	0x40(%%r10),%%xmm0	\n\t	mulpd	%%xmm8 ,%%xmm12		\n\t"/* Restore 2 */\
		"movaps	0x50(%%r10),%%xmm1	\n\t	mulpd	%%xmm8 ,%%xmm13		\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t	mulpd	%%xmm9 ,%%xmm14		\n\t	 movq	%[__two],%%rsi	\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t	mulpd	%%xmm9 ,%%xmm15		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t	movaps	0xc0(%%r10),%%xmm8 	\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t	movaps	0xd0(%%r10),%%xmm9 	\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t	addpd	%%xmm14,%%xmm13		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t	subpd	%%xmm15,%%xmm12		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t	movaps	%%xmm13,%%xmm15		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t	movaps	%%xmm12,%%xmm14		\n\t"\
		/* Finish radix-4 butterfly and store results: */\
		"movaps	    (%%r10),%%xmm0	\n\t	subpd	%%xmm8 ,%%xmm12		\n\t"/* Restore 1 */\
		"movaps	0x10(%%r10),%%xmm1	\n\t	subpd	%%xmm9 ,%%xmm13		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t	addpd	%%xmm8 ,%%xmm14		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t	addpd	%%xmm9 ,%%xmm15		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t	movaps	0x80(%%r10),%%xmm8 	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t	movaps	0x90(%%r10),%%xmm9 	\n\t"\
		"movaps	%%xmm0,0x40(%%r10)	\n\t	subpd	%%xmm14,%%xmm8 		\n\t"	/* 2.0, shared by both columns ... moved +- until found best cycle count: */\
		"movaps	%%xmm2,0x20(%%r10)	\n\t	subpd	%%xmm13,%%xmm10		\n\t	movaps	(%%rsi),%%xmm0	\n\t"\
		"movaps	%%xmm1,0x50(%%r10)	\n\t	subpd	%%xmm15,%%xmm9 		\n\t"\
		"movaps	%%xmm3,0x70(%%r10)	\n\t	subpd	%%xmm12,%%xmm11		\n\t"\
		"mulpd	%%xmm0,%%xmm6		\n\t	movaps	%%xmm8 ,0xc0(%%r10)	\n\t"\
		"mulpd	%%xmm0,%%xmm5		\n\t	movaps	%%xmm10,0xa0(%%r10)	\n\t"\
		"mulpd	%%xmm0,%%xmm7		\n\t	movaps	%%xmm9 ,0xd0(%%r10)	\n\t"\
		"mulpd	%%xmm0,%%xmm4		\n\t	movaps	%%xmm11,0xf0(%%r10)	\n\t"\
		"addpd	0x40(%%r10),%%xmm6	\n\t	mulpd	%%xmm0 ,%%xmm14		\n\t"\
		"addpd		%%xmm2 ,%%xmm5	\n\t	mulpd	%%xmm0 ,%%xmm13		\n\t"\
		"addpd		%%xmm1 ,%%xmm7	\n\t	mulpd	%%xmm0 ,%%xmm15		\n\t"\
		"addpd		%%xmm3 ,%%xmm4	\n\t	mulpd	%%xmm0 ,%%xmm12		\n\t"\
		"movaps	%%xmm6,    (%%r10)	\n\t	addpd	0xc0(%%r10),%%xmm14	\n\t"/* don't need reload-from-mem of xmm8/0xc0(%%r10) as we do in lcol, but 1 cycle faster with it. [!?] */\
		"movaps	%%xmm5,0x60(%%r10)	\n\t	addpd		%%xmm10,%%xmm13	\n\t"\
		"movaps	%%xmm7,0x10(%%r10)	\n\t	addpd		%%xmm9 ,%%xmm15	\n\t"\
		"movaps	%%xmm4,0x30(%%r10)	\n\t	addpd		%%xmm11,%%xmm12	\n\t"\
		"									movaps	%%xmm14,0x80(%%r10)	\n\t"\
		"									movaps	%%xmm13,0xe0(%%r10)	\n\t"\
		"									movaps	%%xmm15,0x90(%%r10)	\n\t"\
		"									movaps	%%xmm12,0xb0(%%r10)	\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%rcx)\n\t"/* ... + p8 */\
	/*
		i0 = (vec_dbl *)(add0+p1); i1 = (vec_dbl *)(add0+p4+p1); i2 = (vec_dbl *)(add0+p8+p1); i3 = (vec_dbl *)(add0+p12+p1);
		o0 += 16; o1 += 16; o2 += 16; o3 += 16;
		SSE2_RADIX_04_DIF_3TWIDDLE_X1(i0,i1,i2,i3,k1, two,c_tmp, o0,o1,o2,o3,k2)
	*/\
		"movslq	%[__p1],%%r9			\n\t	addq	$0x100,%%r10		\n\t"\
		"leaq	(%%rax,%%r9,8),%%rax	\n\t	movq	%[__cc0],%%rsi 		\n\t"/* repoint rsi from two -> cc0 */\
		"leaq	(%%rbx,%%r9,8),%%rbx	\n\t"\
		"leaq	(%%rcx,%%r9,8),%%rcx	\n\t"\
		"leaq	(%%rdx,%%r9,8),%%rdx	\n\t"\
		/* Do	the p0,p2 combo: */\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rsi),%%xmm10	\n\t"\
		"movaps	0x10(%%rsi),%%xmm11	\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t	movaps	    (%%rcx,%%rdi),%%xmm12	\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t	movaps	0x10(%%rcx,%%rdi),%%xmm13	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"mulpd	%%xmm10,%%xmm4		\n\t	movaps	%%xmm12,%%xmm14		\n\t"\
		"mulpd	%%xmm10,%%xmm5		\n\t	movaps	%%xmm13,%%xmm15		\n\t"\
		"mulpd	%%xmm11,%%xmm6		\n\t	movaps	    (%%rax,%%rdi),%%xmm8 	\n\t	movaps	%%xmm0,%%xmm2		\n\t"\
		"mulpd	%%xmm11,%%xmm7		\n\t	movaps	0x10(%%rax,%%rdi),%%xmm9 	\n\t	movaps	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t	mulpd	%%xmm10,%%xmm12		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t	mulpd	%%xmm10,%%xmm13		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t	mulpd	%%xmm11,%%xmm14		\n\t	movaps	%%xmm8 ,%%xmm10		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t	mulpd	%%xmm11,%%xmm15		\n\t	movaps	%%xmm9 ,%%xmm11		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t	addpd	%%xmm14,%%xmm13		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t	subpd	%%xmm15,%%xmm12		\n\t"\
		"movaps	%%xmm0,    (%%r10)	\n\t	addpd	%%xmm12,%%xmm8 		\n\t"/* Spill 1: free up xmm0,1 */\
		"movaps	%%xmm1,0x10(%%r10)	\n\t	addpd	%%xmm13,%%xmm9 		\n\t"/* Do	the p1,3 combo: */\
		"movaps	0x40(%%rsi),%%xmm0	\n\t	subpd	%%xmm12,%%xmm10		\n\t"\
		"movaps	0x50(%%rsi),%%xmm1	\n\t	subpd	%%xmm13,%%xmm11		\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t	movaps	%%xmm8 ,0x80(%%r10)	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t	movaps	%%xmm9 ,0x90(%%r10)	\n\t"\
		"movaps	%%xmm6,%%xmm4		\n\t"\
		"movaps	%%xmm7,%%xmm5		\n\t"\
		"mulpd	%%xmm0,%%xmm4		\n\t	movaps	    (%%rdx,%%rdi),%%xmm14	\n\t"\
		"mulpd	%%xmm0,%%xmm5		\n\t	movaps	0x10(%%rdx,%%rdi),%%xmm15	\n\t"\
		"mulpd	%%xmm1,%%xmm6		\n\t	movaps	%%xmm14,%%xmm12		\n\t"\
		"mulpd	%%xmm1,%%xmm7		\n\t	movaps	%%xmm15,%%xmm13		\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t	mulpd	%%xmm0 ,%%xmm12		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t	mulpd	%%xmm0 ,%%xmm13		\n\t"\
		"movaps	%%xmm5,0x50(%%r10)	\n\t	mulpd	%%xmm1 ,%%xmm14		\n\t"/* Spill 2 */\
		"movaps	%%xmm4,0x40(%%r10)	\n\t	mulpd	%%xmm1 ,%%xmm15		\n\t"\
		"movaps	0x20(%%rsi),%%xmm8	\n\t	addpd	%%xmm14,%%xmm13		\n\t"\
		"movaps	0x30(%%rsi),%%xmm9	\n\t	subpd	%%xmm15,%%xmm12		\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t	movaps	%%xmm13,0xd0(%%r10)	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t	movaps	%%xmm12,0xc0(%%r10)	\n\t"\
		"movaps	%%xmm6,%%xmm4		\n\t"\
		"movaps	%%xmm7,%%xmm5		\n\t"\
		"mulpd	%%xmm8,%%xmm4		\n\t	movaps	    (%%rbx,%%rdi),%%xmm14	\n\t"\
		"mulpd	%%xmm8,%%xmm5		\n\t	movaps	0x10(%%rbx,%%rdi),%%xmm15	\n\t"\
		"mulpd	%%xmm9,%%xmm6		\n\t	movaps	%%xmm14,%%xmm12		\n\t"\
		"mulpd	%%xmm9,%%xmm7		\n\t	movaps	%%xmm15,%%xmm13		\n\t"\
		"movaps	0x40(%%r10),%%xmm0	\n\t	mulpd	%%xmm8 ,%%xmm12		\n\t"/* Restore 2 */\
		"movaps	0x50(%%r10),%%xmm1	\n\t	mulpd	%%xmm8 ,%%xmm13		\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t	mulpd	%%xmm9 ,%%xmm14		\n\t	 movq	%[__two],%%rsi	\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t	mulpd	%%xmm9 ,%%xmm15		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t	movaps	0xc0(%%r10),%%xmm8 	\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t	movaps	0xd0(%%r10),%%xmm9 	\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t	addpd	%%xmm14,%%xmm13		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t	subpd	%%xmm15,%%xmm12		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t	movaps	%%xmm13,%%xmm15		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t	movaps	%%xmm12,%%xmm14		\n\t"\
		/* Finish radix-4 butterfly and store results: */\
		"movaps	    (%%r10),%%xmm0	\n\t	subpd	%%xmm8 ,%%xmm12		\n\t"/* Restore 1 */\
		"movaps	0x10(%%r10),%%xmm1	\n\t	subpd	%%xmm9 ,%%xmm13		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t	addpd	%%xmm8 ,%%xmm14		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t	addpd	%%xmm9 ,%%xmm15		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t	movaps	0x80(%%r10),%%xmm8 	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t	movaps	0x90(%%r10),%%xmm9 	\n\t"\
		"movaps	%%xmm0,0x40(%%r10)	\n\t	subpd	%%xmm14,%%xmm8 		\n\t"	/* 2.0, shared by both columns ... moved +- until found best cycle count: */\
		"movaps	%%xmm2,0x20(%%r10)	\n\t	subpd	%%xmm13,%%xmm10		\n\t	movaps	(%%rsi),%%xmm0	\n\t"\
		"movaps	%%xmm1,0x50(%%r10)	\n\t	subpd	%%xmm15,%%xmm9 		\n\t"\
		"movaps	%%xmm3,0x70(%%r10)	\n\t	subpd	%%xmm12,%%xmm11		\n\t"\
		"mulpd	%%xmm0,%%xmm6		\n\t	movaps	%%xmm8 ,0xc0(%%r10)	\n\t"\
		"mulpd	%%xmm0,%%xmm5		\n\t	movaps	%%xmm10,0xa0(%%r10)	\n\t"\
		"mulpd	%%xmm0,%%xmm7		\n\t	movaps	%%xmm9 ,0xd0(%%r10)	\n\t"\
		"mulpd	%%xmm0,%%xmm4		\n\t	movaps	%%xmm11,0xf0(%%r10)	\n\t"\
		"addpd	0x40(%%r10),%%xmm6	\n\t	mulpd	%%xmm0 ,%%xmm14		\n\t"\
		"addpd		%%xmm2 ,%%xmm5	\n\t	mulpd	%%xmm0 ,%%xmm13		\n\t"\
		"addpd		%%xmm1 ,%%xmm7	\n\t	mulpd	%%xmm0 ,%%xmm15		\n\t"\
		"addpd		%%xmm3 ,%%xmm4	\n\t	mulpd	%%xmm0 ,%%xmm12		\n\t"\
		"movaps	%%xmm6,    (%%r10)	\n\t	addpd	0xc0(%%r10),%%xmm14	\n\t"/* don't need reload-from-mem of xmm8/0xc0(%%r10) as we do in lcol, but 1 cycle faster with it. [!?] */\
		"movaps	%%xmm5,0x60(%%r10)	\n\t	addpd		%%xmm10,%%xmm13	\n\t"\
		"movaps	%%xmm7,0x10(%%r10)	\n\t	addpd		%%xmm9 ,%%xmm15	\n\t"\
		"movaps	%%xmm4,0x30(%%r10)	\n\t	addpd		%%xmm11,%%xmm12	\n\t"\
		"									movaps	%%xmm14,0x80(%%r10)	\n\t"\
		"									movaps	%%xmm13,0xe0(%%r10)	\n\t"\
		"									movaps	%%xmm15,0x90(%%r10)	\n\t"\
		"									movaps	%%xmm12,0xb0(%%r10)	\n\t"\
	/*
		// Pass 2:
		k1 = 0x40;	// Set k1 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k2 = p4*8;	// Set k2 to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		// 0x60 = literal-bytewise address offset between the pointers to the first and second DFT's twiddles = 0:
		i0 = r1; i1 = r1+16; i2 = r1+8; i3 = r1+24;
		o0 = (vec_dbl *)add0; o1 = (vec_dbl *)(add0+p2); o2 = (vec_dbl *)(add0+p1); o3 = (vec_dbl *)(add0+p3);
		c_tmp += 6;	// c2,A,6
		SSE2_RADIX_04_DIF_3TWIDDLE_X2(i0,i1,i2,i3,k1, two,c_tmp,0x60, o0,o1,o2,o3,k2)
	*/\
		/* Take advantage of fact that i1 = r0 + 16 already in %%r10 and o2 = add0+p1 in %%rax : */\
		/* i-offset = 0x40, so inline with r[a-d]x base-address offsets rather than adding via (r[a-d]x,rdi) */\
		"movq	%%r10,%%rbx			\n\t	movq	%%rax,%%r12		\n\t"/* rbx:i1 = r0+16; r12:o2 = add0+p1 */\
		"movq	%[__cc0],%%rsi 		\n\t	addq	$0x60,%%rsi 		\n\t"/* repoint rsi from two -> cc0, then increment for Pass 2 */\
		"leaq	 0x60(%%rsi),%%r8	\n\t	movslq	%[__p2],%%r11		\n\t"\
		"leaq	-0x100(%%rbx),%%rax	\n\t	movslq	%[__p4],%%r9		\n\t"\
		"leaq	-0x080(%%rbx),%%rcx	\n\t	movq	%[__add0]	,%%r10		\n\t"/* o0 */\
		"leaq	 0x080(%%rbx),%%rdx	\n\t	leaq	(%%r12,%%r11,8),%%r13	\n\t"/* o3 = (add0+p1)+p2 */\
		"movaps	    (%%rax),%%xmm0	\n\t	leaq	(%%r10,%%r11,8),%%r11	\n\t"/* o1 = (add0+p2) */\
		"movaps	0x10(%%rax),%%xmm1	\n\t	shlq	$3,%%r9			\n\t"/* p4 in byte-offset form */\
		/* Do	the p0,p2 combo: */\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rsi),%%xmm2	\n\t"\
		"movaps	0x10(%%rsi),%%xmm3	\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r9)\n\t"/* ... + p4 */\
		"movaps	%%xmm4,%%xmm6		\n\t	movaps	0x40(%%rcx),%%xmm12	\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t	movaps	0x50(%%rcx),%%xmm13	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t	movaps	    (%%r8),%%xmm10	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t	movaps	0x10(%%r8),%%xmm11	\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t	movaps	%%xmm12,%%xmm14		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t	movaps	%%xmm13,%%xmm15		\n\t"\
		"mulpd	%%xmm3,%%xmm6		\n\t	movaps	0x40(%%rax),%%xmm8 	\n\t	movaps	%%xmm0,%%xmm2		\n\t"\
		"mulpd	%%xmm3,%%xmm7		\n\t	movaps	0x50(%%rax),%%xmm9 	\n\t	movaps	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t	mulpd	%%xmm10,%%xmm12		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t	mulpd	%%xmm10,%%xmm13		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t	mulpd	%%xmm11,%%xmm14		\n\t	movaps	%%xmm8 ,%%xmm10		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t	mulpd	%%xmm11,%%xmm15		\n\t	movaps	%%xmm9 ,%%xmm11		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t	addpd	%%xmm14,%%xmm13		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t	subpd	%%xmm15,%%xmm12		\n\t"\
		"movaps	%%xmm0,    (%%r10)	\n\t	addpd	%%xmm12,%%xmm8 		\n\t"/* Spill 1: free up xmm0,1 */\
		"movaps	%%xmm1,0x10(%%r10)	\n\t	addpd	%%xmm13,%%xmm9 		\n\t"/* Do	the p1,3 combo: */\
		"movaps	0x40(%%rsi),%%xmm0	\n\t	subpd	%%xmm12,%%xmm10		\n\t"\
		"movaps	0x50(%%rsi),%%xmm1	\n\t	subpd	%%xmm13,%%xmm11		\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t	movaps	%%xmm8 ,    (%%r10,%%r9)	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t	movaps	%%xmm9 ,0x10(%%r10,%%r9)	\n\t"\
		"movaps	%%xmm6,%%xmm4		\n\t	movaps	0x40(%%r8),%%xmm8 	\n\t"\
		"movaps	%%xmm7,%%xmm5		\n\t	movaps	0x50(%%r8),%%xmm9 	\n\t"\
		"mulpd	%%xmm0,%%xmm4		\n\t	movaps	0x40(%%rdx),%%xmm14	\n\t"\
		"mulpd	%%xmm0,%%xmm5		\n\t	movaps	0x50(%%rdx),%%xmm15	\n\t"\
		"mulpd	%%xmm1,%%xmm6		\n\t	movaps	%%xmm14,%%xmm12		\n\t"\
		"mulpd	%%xmm1,%%xmm7		\n\t	movaps	%%xmm15,%%xmm13		\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t	mulpd	%%xmm8 ,%%xmm12		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t	mulpd	%%xmm8 ,%%xmm13		\n\t"\
		"movaps	%%xmm5,0x10(%%r12)	\n\t	mulpd	%%xmm9 ,%%xmm14		\n\t"/* Spill 2 */\
		"movaps	%%xmm4,    (%%r12)	\n\t	mulpd	%%xmm9 ,%%xmm15		\n\t"\
		"movaps	0x20(%%rsi),%%xmm0	\n\t	addpd	%%xmm14,%%xmm13		\n\t"\
		"movaps	0x30(%%rsi),%%xmm1	\n\t	subpd	%%xmm15,%%xmm12		\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t	movaps	%%xmm13,0x10(%%r12,%%r9)	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t	movaps	%%xmm12,    (%%r12,%%r9)	\n\t"\
		"movaps	%%xmm6,%%xmm4		\n\t	movaps	0x20(%%r8),%%xmm8 	\n\t"\
		"movaps	%%xmm7,%%xmm5		\n\t	movaps	0x30(%%r8),%%xmm9 	\n\t"\
		"mulpd	%%xmm0,%%xmm4		\n\t	movaps	0x40(%%rbx),%%xmm14	\n\t"\
		"mulpd	%%xmm0,%%xmm5		\n\t	movaps	0x50(%%rbx),%%xmm15	\n\t"\
		"mulpd	%%xmm1,%%xmm6		\n\t	movaps	%%xmm14,%%xmm12		\n\t"\
		"mulpd	%%xmm1,%%xmm7		\n\t	movaps	%%xmm15,%%xmm13		\n\t"\
		"movaps	    (%%r12),%%xmm0	\n\t	mulpd	%%xmm8 ,%%xmm12		\n\t"/* Restore 2 */\
		"movaps	0x10(%%r12),%%xmm1	\n\t	mulpd	%%xmm8 ,%%xmm13		\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t	mulpd	%%xmm9 ,%%xmm14		\n\t	 movq	%[__two],%%rsi	\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t	mulpd	%%xmm9 ,%%xmm15		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t	movaps	    (%%r12,%%r9),%%xmm8 	\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t	movaps	0x10(%%r12,%%r9),%%xmm9 	\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t	addpd	%%xmm14,%%xmm13		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t	subpd	%%xmm15,%%xmm12		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t	movaps	%%xmm13,%%xmm15		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t	movaps	%%xmm12,%%xmm14		\n\t"\
		/* Finish radix-4 butterfly and store results: */\
		"movaps	    (%%r10),%%xmm0	\n\t	subpd	%%xmm8 ,%%xmm12		\n\t"/* Restore 1 */\
		"movaps	0x10(%%r10),%%xmm1	\n\t	subpd	%%xmm9 ,%%xmm13		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t	addpd	%%xmm8 ,%%xmm14		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t	addpd	%%xmm9 ,%%xmm15		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t	movaps	    (%%r10,%%r9),%%xmm8 	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t	movaps	0x10(%%r10,%%r9),%%xmm9 	\n\t"\
		"movaps	%%xmm0,    (%%r12)	\n\t	subpd	%%xmm14,%%xmm8 		\n\t"	/* 2.0, shared by both columns ... moved +- until found best cycle count: */\
		"movaps	%%xmm2,    (%%r11)	\n\t	subpd	%%xmm13,%%xmm10		\n\t	movaps	(%%rsi),%%xmm0	\n\t"\
		"movaps	%%xmm1,0x10(%%r12)	\n\t	subpd	%%xmm15,%%xmm9 		\n\t"\
		"movaps	%%xmm3,0x10(%%r13)	\n\t	subpd	%%xmm12,%%xmm11		\n\t"\
		"mulpd	%%xmm0,%%xmm6		\n\t	movaps	%%xmm8 ,    (%%r12,%%r9)	\n\t"\
		"mulpd	%%xmm0,%%xmm5		\n\t	movaps	%%xmm10,    (%%r11,%%r9)	\n\t"\
		"mulpd	%%xmm0,%%xmm7		\n\t	movaps	%%xmm9 ,0x10(%%r12,%%r9)	\n\t"\
		"mulpd	%%xmm0,%%xmm4		\n\t	movaps	%%xmm11,0x10(%%r13,%%r9)	\n\t"\
		"addpd	    (%%r12),%%xmm6	\n\t	mulpd	%%xmm0 ,%%xmm14		\n\t"\
		"addpd		%%xmm2 ,%%xmm5	\n\t	mulpd	%%xmm0 ,%%xmm13		\n\t"\
		"addpd		%%xmm1 ,%%xmm7	\n\t	mulpd	%%xmm0 ,%%xmm15		\n\t"\
		"addpd		%%xmm3 ,%%xmm4	\n\t	mulpd	%%xmm0 ,%%xmm12		\n\t"\
		"movaps	%%xmm6,    (%%r10)	\n\t	addpd	    (%%r12,%%r9),%%xmm14	\n\t"/* don't need reload-from-mem of xmm8/0xc0(%%r12,%%r9) as we do in lcol, but 1 cycle faster with it. [!?] */\
		"movaps	%%xmm5,    (%%r13)	\n\t	addpd		%%xmm10,%%xmm13	\n\t"\
		"movaps	%%xmm7,0x10(%%r10)	\n\t	addpd		%%xmm9 ,%%xmm15	\n\t"\
		"movaps	%%xmm4,0x10(%%r11)	\n\t	addpd		%%xmm11,%%xmm12	\n\t"\
		"									movaps	%%xmm14,    (%%r10,%%r9)	\n\t"\
		"									movaps	%%xmm13,    (%%r13,%%r9)	\n\t"\
		"									movaps	%%xmm15,0x10(%%r10,%%r9)	\n\t"\
		"									movaps	%%xmm12,0x10(%%r11,%%r9)	\n\t"\
	/*
		i0 += 2; i1 += 2; i2 += 2; i3 += 2;
		o0 = (vec_dbl *)(add0+p8); o1 = (vec_dbl *)(add0+p8+p2); o2 = (vec_dbl *)(add0+p8+p1); o3 = (vec_dbl *)(add0+p8+p3);
		c_tmp += 12;	// c5,D,3
		SSE2_RADIX_04_DIF_3TWIDDLE_X2(i0,i1,i2,i3,k1, two,c_tmp,0x60, o0,o1,o2,o3,k2)
	*/\
		/* Take advantage of fact that i1 = r0 + 16 already in %%r10 and o2 = add0+p1 in %%rax : */\
		"movq	%[__cc0],%%rsi 		\n\t	addq	$0x120,%%rsi 		\n\t"/* repoint rsi from two -> cc0, then increment for Pass, set 2 */\
		"addq	$0x20,%%rax			\n\t	movslq	%[__p8],%%r8		\n\t"\
		"addq	$0x20,%%rbx			\n\t"	/* r9 still has p4<<3 */\
		"addq	$0x20,%%rcx			\n\t	leaq	(%%r10,%%r8,8),%%r10	\n\t"/* o0 =  add0    +p8 */\
		"addq	$0x20,%%rdx			\n\t	leaq	(%%r11,%%r8,8),%%r11	\n\t"/* o1 = (add0+p2)+p8 */\
		"movaps	    (%%rax),%%xmm0	\n\t	leaq	(%%r12,%%r8,8),%%r12	\n\t"/* o2 = (add0+p1)+p8 */\
		"movaps	0x10(%%rax),%%xmm1	\n\t	leaq	(%%r13,%%r8,8),%%r13	\n\t"/* o3 = (add0+p3)+p8 */\
		/* Do	the p0,p2 combo: */\
		"movaps	    (%%rcx),%%xmm4	\n\t	leaq	 0x60(%%rsi),%%r8	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rsi),%%xmm2	\n\t"\
		"movaps	0x10(%%rsi),%%xmm3	\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r10,%%r9)\n\t"/* ... + pC */\
		"movaps	%%xmm4,%%xmm6		\n\t	movaps	0x40(%%rcx),%%xmm12	\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t	movaps	0x50(%%rcx),%%xmm13	\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t	movaps	    (%%r8),%%xmm10	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t	movaps	0x10(%%r8),%%xmm11	\n\t"\
		"mulpd	%%xmm2,%%xmm4		\n\t	movaps	%%xmm12,%%xmm14		\n\t"\
		"mulpd	%%xmm2,%%xmm5		\n\t	movaps	%%xmm13,%%xmm15		\n\t"\
		"mulpd	%%xmm3,%%xmm6		\n\t	movaps	0x40(%%rax),%%xmm8 	\n\t	movaps	%%xmm0,%%xmm2		\n\t"\
		"mulpd	%%xmm3,%%xmm7		\n\t	movaps	0x50(%%rax),%%xmm9 	\n\t	movaps	%%xmm1,%%xmm3		\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t	mulpd	%%xmm10,%%xmm12		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t	mulpd	%%xmm10,%%xmm13		\n\t"\
		"addpd	%%xmm4,%%xmm0		\n\t	mulpd	%%xmm11,%%xmm14		\n\t	movaps	%%xmm8 ,%%xmm10		\n\t"\
		"addpd	%%xmm5,%%xmm1		\n\t	mulpd	%%xmm11,%%xmm15		\n\t	movaps	%%xmm9 ,%%xmm11		\n\t"\
		"subpd	%%xmm4,%%xmm2		\n\t	addpd	%%xmm14,%%xmm13		\n\t"\
		"subpd	%%xmm5,%%xmm3		\n\t	subpd	%%xmm15,%%xmm12		\n\t"\
		"movaps	%%xmm0,    (%%r10)	\n\t	addpd	%%xmm12,%%xmm8 		\n\t"/* Spill 1: free up xmm0,1 */\
		"movaps	%%xmm1,0x10(%%r10)	\n\t	addpd	%%xmm13,%%xmm9 		\n\t"/* Do	the p1,3 combo: */\
		"movaps	0x40(%%rsi),%%xmm0	\n\t	subpd	%%xmm12,%%xmm10		\n\t"\
		"movaps	0x50(%%rsi),%%xmm1	\n\t	subpd	%%xmm13,%%xmm11		\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t	movaps	%%xmm8 ,    (%%r10,%%r9)	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t	movaps	%%xmm9 ,0x10(%%r10,%%r9)	\n\t"\
		"movaps	%%xmm6,%%xmm4		\n\t	movaps	0x40(%%r8),%%xmm8 	\n\t"\
		"movaps	%%xmm7,%%xmm5		\n\t	movaps	0x50(%%r8),%%xmm9 	\n\t"\
		"mulpd	%%xmm0,%%xmm4		\n\t	movaps	0x40(%%rdx),%%xmm14	\n\t"\
		"mulpd	%%xmm0,%%xmm5		\n\t	movaps	0x50(%%rdx),%%xmm15	\n\t"\
		"mulpd	%%xmm1,%%xmm6		\n\t	movaps	%%xmm14,%%xmm12		\n\t"\
		"mulpd	%%xmm1,%%xmm7		\n\t	movaps	%%xmm15,%%xmm13		\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t	mulpd	%%xmm8 ,%%xmm12		\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t	mulpd	%%xmm8 ,%%xmm13		\n\t"\
		"movaps	%%xmm5,0x10(%%r12)	\n\t	mulpd	%%xmm9 ,%%xmm14		\n\t"/* Spill 2 */\
		"movaps	%%xmm4,    (%%r12)	\n\t	mulpd	%%xmm9 ,%%xmm15		\n\t"\
		"movaps	0x20(%%rsi),%%xmm0	\n\t	addpd	%%xmm14,%%xmm13		\n\t"\
		"movaps	0x30(%%rsi),%%xmm1	\n\t	subpd	%%xmm15,%%xmm12		\n\t"\
		"movaps	    (%%rbx),%%xmm6	\n\t	movaps	%%xmm13,0x10(%%r12,%%r9)	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7	\n\t	movaps	%%xmm12,    (%%r12,%%r9)	\n\t"\
		"movaps	%%xmm6,%%xmm4		\n\t	movaps	0x20(%%r8),%%xmm8 	\n\t"\
		"movaps	%%xmm7,%%xmm5		\n\t	movaps	0x30(%%r8),%%xmm9 	\n\t"\
		"mulpd	%%xmm0,%%xmm4		\n\t	movaps	0x40(%%rbx),%%xmm14	\n\t"\
		"mulpd	%%xmm0,%%xmm5		\n\t	movaps	0x50(%%rbx),%%xmm15	\n\t"\
		"mulpd	%%xmm1,%%xmm6		\n\t	movaps	%%xmm14,%%xmm12		\n\t"\
		"mulpd	%%xmm1,%%xmm7		\n\t	movaps	%%xmm15,%%xmm13		\n\t"\
		"movaps	    (%%r12),%%xmm0	\n\t	mulpd	%%xmm8 ,%%xmm12		\n\t"/* Restore 2 */\
		"movaps	0x10(%%r12),%%xmm1	\n\t	mulpd	%%xmm8 ,%%xmm13		\n\t"\
		"addpd	%%xmm6,%%xmm5		\n\t	mulpd	%%xmm9 ,%%xmm14		\n\t	 movq	%[__two],%%rsi	\n\t"\
		"subpd	%%xmm7,%%xmm4		\n\t	mulpd	%%xmm9 ,%%xmm15		\n\t"\
		"movaps	%%xmm5,%%xmm7		\n\t	movaps	    (%%r12,%%r9),%%xmm8 	\n\t"\
		"movaps	%%xmm4,%%xmm6		\n\t	movaps	0x10(%%r12,%%r9),%%xmm9 	\n\t"\
		"subpd	%%xmm0,%%xmm4		\n\t	addpd	%%xmm14,%%xmm13		\n\t"\
		"subpd	%%xmm1,%%xmm5		\n\t	subpd	%%xmm15,%%xmm12		\n\t"\
		"addpd	%%xmm0,%%xmm6		\n\t	movaps	%%xmm13,%%xmm15		\n\t"\
		"addpd	%%xmm1,%%xmm7		\n\t	movaps	%%xmm12,%%xmm14		\n\t"\
		/* Finish radix-4 butterfly and store results: */\
		"movaps	    (%%r10),%%xmm0	\n\t	subpd	%%xmm8 ,%%xmm12		\n\t"/* Restore 1 */\
		"movaps	0x10(%%r10),%%xmm1	\n\t	subpd	%%xmm9 ,%%xmm13		\n\t"\
		"subpd	%%xmm6,%%xmm0		\n\t	addpd	%%xmm8 ,%%xmm14		\n\t"\
		"subpd	%%xmm5,%%xmm2		\n\t	addpd	%%xmm9 ,%%xmm15		\n\t"\
		"subpd	%%xmm7,%%xmm1		\n\t	movaps	    (%%r10,%%r9),%%xmm8 	\n\t"\
		"subpd	%%xmm4,%%xmm3		\n\t	movaps	0x10(%%r10,%%r9),%%xmm9 	\n\t"\
		"movaps	%%xmm0,    (%%r12)	\n\t	subpd	%%xmm14,%%xmm8 		\n\t"	/* 2.0, shared by both columns ... moved +- until found best cycle count: */\
		"movaps	%%xmm2,    (%%r11)	\n\t	subpd	%%xmm13,%%xmm10		\n\t	movaps	(%%rsi),%%xmm0	\n\t"\
		"movaps	%%xmm1,0x10(%%r12)	\n\t	subpd	%%xmm15,%%xmm9 		\n\t"\
		"movaps	%%xmm3,0x10(%%r13)	\n\t	subpd	%%xmm12,%%xmm11		\n\t"\
		"mulpd	%%xmm0,%%xmm6		\n\t	movaps	%%xmm8 ,    (%%r12,%%r9)	\n\t"\
		"mulpd	%%xmm0,%%xmm5		\n\t	movaps	%%xmm10,    (%%r11,%%r9)	\n\t"\
		"mulpd	%%xmm0,%%xmm7		\n\t	movaps	%%xmm9 ,0x10(%%r12,%%r9)	\n\t"\
		"mulpd	%%xmm0,%%xmm4		\n\t	movaps	%%xmm11,0x10(%%r13,%%r9)	\n\t"\
		"addpd	    (%%r12),%%xmm6	\n\t	mulpd	%%xmm0 ,%%xmm14		\n\t"\
		"addpd		%%xmm2 ,%%xmm5	\n\t	mulpd	%%xmm0 ,%%xmm13		\n\t"\
		"addpd		%%xmm1 ,%%xmm7	\n\t	mulpd	%%xmm0 ,%%xmm15		\n\t"\
		"addpd		%%xmm3 ,%%xmm4	\n\t	mulpd	%%xmm0 ,%%xmm12		\n\t"\
		"movaps	%%xmm6,    (%%r10)	\n\t	addpd	    (%%r12,%%r9),%%xmm14	\n\t"/* don't need reload-from-mem of xmm8/0xc0(%%r12,%%r9) as we do in lcol, but 1 cycle faster with it. [!?] */\
		"movaps	%%xmm5,    (%%r13)	\n\t	addpd		%%xmm10,%%xmm13	\n\t"\
		"movaps	%%xmm7,0x10(%%r10)	\n\t	addpd		%%xmm9 ,%%xmm15	\n\t"\
		"movaps	%%xmm4,0x10(%%r11)	\n\t	addpd		%%xmm11,%%xmm12	\n\t"\
		"									movaps	%%xmm14,    (%%r10,%%r9)	\n\t"\
		"									movaps	%%xmm13,    (%%r13,%%r9)	\n\t"\
		"									movaps	%%xmm15,0x10(%%r10,%%r9)	\n\t"\
		"									movaps	%%xmm12,0x10(%%r11,%%r9)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r0] "m" (Xr0)\
		 ,[__two] "m" (Xtwo)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_DIT_TWIDDLE_V2(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr0,Xtwo,Xcc0,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
	/*
		// Pass 1:
		j = p4*8;	// Set j to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k = 0x80;	// Set k to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		// 0x60 = literal-bytewise address offset between the pointers to the first and second DFT's twiddles:
		i0 = (vec_dbl *)a; i1 = (vec_dbl *)(a+p1); i2 = (vec_dbl *)(a+p2); i3 = (vec_dbl *)(a+p3);
		o0 = r00; o1 = r00+2; o2 = r00+4; o3 = r00+6;
		c_tmp = cc0+6;	// c2,1,3
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,0x60, o0,o1,o2,o3,k)
	*/\
		"movq	%[__two],%%rsi 		\n\t	movslq	%[__p4],%%rdi		\n\t"\
		"movq	%[__add0],%%rax		\n\t	movq	%[__r0],%%r10		\n\t"\
		"movslq	%[__p1],%%rbx		\n\t	leaq	(%%rax,%%rbx,8),%%rbx	\n\t"\
		"movslq	%[__p2],%%rcx		\n\t	leaq	(%%rax,%%rcx,8),%%rcx	\n\t"\
		"movslq	%[__p3],%%rdx		\n\t	leaq	(%%rax,%%rdx,8),%%rdx	\n\t"\
	"movaps	(%%rsi),%%xmm15	\n\t"/* two */"	shlq	$3,%%rdi			\n\t"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t	leaq	0x10(%%rdx,%%rdi),%%r11		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t	movaps	    (%%rax,%%rdi),%%xmm8 	\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t	movaps	0x10(%%rax,%%rdi),%%xmm9 	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t	movaps	    (%%rbx,%%rdi),%%xmm10	\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t	movaps	0x10(%%rbx,%%rdi),%%xmm11	\n\t"\
		"mulpd	%%xmm15,%%xmm2		\n\t	movaps	    (%%rcx,%%rdi),%%xmm12	\n\t"\
		"mulpd	%%xmm15,%%xmm3		\n\t	movaps	0x10(%%rcx,%%rdi),%%xmm13	\n\t"\
		"mulpd	%%xmm15,%%xmm6		\n\t	movaps	    (%%rdx,%%rdi),%%xmm14	\n\t"\
		"mulpd	%%xmm15,%%xmm7		\n\t"/*	movaps	0x10(%%rdx,%%rdi),%%xmm15	\n\t"*/\
		"addpd	%%xmm0,%%xmm2		\n\t	subpd	%%xmm10,%%xmm8 		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t	subpd	%%xmm11,%%xmm9 		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t	subpd	%%xmm14,%%xmm12		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t	subpd	(%%r11),%%xmm13		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t	mulpd	%%xmm15,%%xmm10		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t	mulpd	%%xmm15,%%xmm11		\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t	mulpd	%%xmm15,%%xmm14		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t	mulpd	(%%r11),%%xmm15		\n\t"\
		"mulpd	(%%rsi),%%xmm6		\n\t	addpd	%%xmm8 ,%%xmm10		\n\t"\
		"mulpd	(%%rsi),%%xmm7		\n\t	addpd	%%xmm9 ,%%xmm11		\n\t"\
		"mulpd	(%%rsi),%%xmm5		\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
		"mulpd	(%%rsi),%%xmm4		\n\t	addpd	%%xmm13,%%xmm15		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t	subpd	%%xmm14,%%xmm10		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t	subpd	%%xmm15,%%xmm11		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t	subpd	%%xmm13,%%xmm8 		\n\t"\
		"addpd	%%xmm1,%%xmm4		\n\t	subpd	%%xmm12,%%xmm9 		\n\t"\
		"movaps	%%xmm6,    (%%r10)	\n\t	movaps	%%xmm7,0x10(%%r10)	\n\t"/* Write B0 to free up 2 regs */\
		"movaps	(%%rsi),%%xmm6		\n\t	mulpd	%%xmm6 ,%%xmm14		\n\t"\
		"									mulpd	%%xmm6 ,%%xmm15		\n\t"\
		"									mulpd	%%xmm6 ,%%xmm13		\n\t"\
		/* B2 = t0*~w2 = t0*[c2-I.s1] */"	mulpd	%%xmm6 ,%%xmm12		\n\t	movq	%[__cc0],%%rsi 	\n\t"\
		"movaps	%%xmm2,%%xmm6		\n\t	addpd	%%xmm10,%%xmm14		\n\t	addq	$0x60,%%rsi 	\n\t"/* cc0 + 6 */\
		"movaps	%%xmm3,%%xmm7		\n\t	addpd	%%xmm11,%%xmm15		\n\t"\
		"mulpd	    (%%rsi),%%xmm6	\n\t	addpd	%%xmm8 ,%%xmm13		\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t	addpd	%%xmm9 ,%%xmm12		\n\t"\
		"mulpd	    (%%rsi),%%xmm7	\n\t	movaps %%xmm14,0x80(%%r10)	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t	movaps %%xmm15,0x90(%%r10)	\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t	movaps	%%xmm10,%%xmm14		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t	movaps	%%xmm11,%%xmm15		\n\t"\
		"movaps	%%xmm6,0x40(%%r10)	\n\t	mulpd	0x60(%%rsi),%%xmm14	\n\t"\
		"movaps	%%xmm7,0x50(%%r10)	\n\t	mulpd	0x70(%%rsi),%%xmm11	\n\t"\
		/* B1 = t1*~w1 = t1*[c1-I.s1] */"	mulpd	0x60(%%rsi),%%xmm15	\n\t"\
		"movaps	0x20(%%rsi),%%xmm2	\n\t	mulpd	0x70(%%rsi),%%xmm10	\n\t"\
		"movaps	0x30(%%rsi),%%xmm3	\n\t	addpd	%%xmm11,%%xmm14		\n\t"\
		"movaps	%%xmm5,%%xmm6		\n\t	subpd	%%xmm10,%%xmm15		\n\t"\
		"movaps	%%xmm1,%%xmm7		\n\t	movaps %%xmm14,0xc0(%%r10)	\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t	movaps %%xmm15,0xd0(%%r10)	\n\t"\
		"mulpd	%%xmm3,%%xmm1		\n\t	movaps	%%xmm13,%%xmm14		\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t	movaps	%%xmm9 ,%%xmm15		\n\t"\
		"mulpd	%%xmm3,%%xmm5		\n\t	mulpd	0x80(%%rsi),%%xmm14	\n\t"\
		"movaps	0x40(%%rsi),%%xmm2	\n\t	mulpd	0x90(%%rsi),%%xmm9 	\n\t"\
		"movaps	0x50(%%rsi),%%xmm3	\n\t	mulpd	0x80(%%rsi),%%xmm15	\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t	mulpd	0x90(%%rsi),%%xmm13	\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t	addpd	%%xmm9 ,%%xmm14		\n\t"\
		"movaps	%%xmm6,0x20(%%r10)	\n\t	subpd	%%xmm13,%%xmm15		\n\t"\
		"movaps	%%xmm7,0x30(%%r10)	\n\t	movaps %%xmm14,0xa0(%%r10)	\n\t"\
		/* B3 = t3*~w3 = t3*[c3-I.s3] */"	movaps %%xmm15,0xb0(%%r10)	\n\t"\
		"movaps	%%xmm0,%%xmm6		\n\t	movaps	%%xmm8 ,%%xmm14		\n\t"\
		"movaps	%%xmm4,%%xmm7		\n\t	movaps	%%xmm12,%%xmm15		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t	mulpd	0xa0(%%rsi),%%xmm14	\n\t"\
		"mulpd	%%xmm3,%%xmm4		\n\t	mulpd	0xb0(%%rsi),%%xmm12	\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t	mulpd	0xa0(%%rsi),%%xmm15	\n\t"\
		"mulpd	%%xmm3,%%xmm0		\n\t	mulpd	0xb0(%%rsi),%%xmm8 	\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
		"subpd	%%xmm0,%%xmm7		\n\t	subpd	%%xmm8 ,%%xmm15		\n\t"\
		"movaps	%%xmm6,0x60(%%r10)	\n\t	movaps	%%xmm14,0xe0(%%r10)	\n\t"\
		"movaps	%%xmm7,0x70(%%r10)	\n\t	movaps	%%xmm15,0xf0(%%r10)	\n\t"\
	/*
		i0 = (vec_dbl *)(a+p8); i1 = (vec_dbl *)(a+p9); i2 = (vec_dbl *)(a+pA); i3 = (vec_dbl *)(a+pB);
		o0 += 16; o1 += 16; o2 += 16; o3 += 16;
		c_tmp += 12;		// cA,9,B
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,0x60, o0,o1,o2,o3,k)
	*/\
		"movq	%[__two],%%rsi 		\n\t	movslq	%[__p8],%%r11		\n\t"\
		"leaq (%%rax,%%r11,8),%%rax	\n\t	addq	 $0x100,%%r10		\n\t"/* r0 + 16 */\
		"leaq (%%rbx,%%r11,8),%%rbx	\n\t"\
		"leaq (%%rcx,%%r11,8),%%rcx	\n\t"\
		"leaq (%%rdx,%%r11,8),%%rdx	\n\t"\
	"movaps	(%%rsi),%%xmm15	\n\t"/* two */"	"\
		"movaps	    (%%rax),%%xmm0	\n\t"\
		"movaps	0x10(%%rax),%%xmm1	\n\t"\
		"movaps	    (%%rbx),%%xmm2	\n\t"\
		"movaps	0x10(%%rbx),%%xmm3	\n\t"\
		"movaps	    (%%rcx),%%xmm4	\n\t"\
		"movaps	0x10(%%rcx),%%xmm5	\n\t"\
		"movaps	    (%%rdx),%%xmm6	\n\t"\
		"movaps	0x10(%%rdx),%%xmm7	\n\t	leaq	0x10(%%rdx,%%rdi),%%r11		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t	movaps	    (%%rax,%%rdi),%%xmm8 	\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t	movaps	0x10(%%rax,%%rdi),%%xmm9 	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t	movaps	    (%%rbx,%%rdi),%%xmm10	\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t	movaps	0x10(%%rbx,%%rdi),%%xmm11	\n\t"\
		"mulpd	%%xmm15,%%xmm2		\n\t	movaps	    (%%rcx,%%rdi),%%xmm12	\n\t"\
		"mulpd	%%xmm15,%%xmm3		\n\t	movaps	0x10(%%rcx,%%rdi),%%xmm13	\n\t"\
		"mulpd	%%xmm15,%%xmm6		\n\t	movaps	    (%%rdx,%%rdi),%%xmm14	\n\t"\
		"mulpd	%%xmm15,%%xmm7		\n\t"/*	movaps	0x10(%%rdx,%%rdi),%%xmm15	\n\t"*/\
		"addpd	%%xmm0,%%xmm2		\n\t	subpd	%%xmm10,%%xmm8 		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t	subpd	%%xmm11,%%xmm9 		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t	subpd	%%xmm14,%%xmm12		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t	subpd	(%%r11),%%xmm13		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t	mulpd	%%xmm15,%%xmm10		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t	mulpd	%%xmm15,%%xmm11		\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t	mulpd	%%xmm15,%%xmm14		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t	mulpd	(%%r11),%%xmm15		\n\t"\
		"mulpd	(%%rsi),%%xmm6		\n\t	addpd	%%xmm8 ,%%xmm10		\n\t"\
		"mulpd	(%%rsi),%%xmm7		\n\t	addpd	%%xmm9 ,%%xmm11		\n\t"\
		"mulpd	(%%rsi),%%xmm5		\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
		"mulpd	(%%rsi),%%xmm4		\n\t	addpd	%%xmm13,%%xmm15		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t	subpd	%%xmm14,%%xmm10		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t	subpd	%%xmm15,%%xmm11		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t	subpd	%%xmm13,%%xmm8 		\n\t"\
		"addpd	%%xmm1,%%xmm4		\n\t	subpd	%%xmm12,%%xmm9 		\n\t"\
		"movaps	%%xmm6,    (%%r10)	\n\t	movaps	%%xmm7,0x10(%%r10)	\n\t"/* Write B0 to free up 2 regs */\
		"movaps	(%%rsi),%%xmm6		\n\t	mulpd	%%xmm6 ,%%xmm14		\n\t"\
		"									mulpd	%%xmm6 ,%%xmm15		\n\t"\
		"									mulpd	%%xmm6 ,%%xmm13		\n\t"\
		/* B2 = t0*~w2 = t0*[c2-I.s1] */"	mulpd	%%xmm6 ,%%xmm12		\n\t	movq	%[__cc0],%%rsi 	\n\t"\
		"movaps	%%xmm2,%%xmm6		\n\t	addpd	%%xmm10,%%xmm14		\n\t	addq	$0x120,%%rsi 	\n\t"/* cc0 + 18 */\
		"movaps	%%xmm3,%%xmm7		\n\t	addpd	%%xmm11,%%xmm15		\n\t"\
		"mulpd	    (%%rsi),%%xmm6	\n\t	addpd	%%xmm8 ,%%xmm13		\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t	addpd	%%xmm9 ,%%xmm12		\n\t"\
		"mulpd	    (%%rsi),%%xmm7	\n\t	movaps %%xmm14,0x80(%%r10)	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t	movaps %%xmm15,0x90(%%r10)	\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t	movaps	%%xmm10,%%xmm14		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t	movaps	%%xmm11,%%xmm15		\n\t"\
		"movaps	%%xmm6,0x40(%%r10)	\n\t	mulpd	0x60(%%rsi),%%xmm14	\n\t"\
		"movaps	%%xmm7,0x50(%%r10)	\n\t	mulpd	0x70(%%rsi),%%xmm11	\n\t"\
		/* B1 = t1*~w1 = t1*[c1-I.s1] */"	mulpd	0x60(%%rsi),%%xmm15	\n\t"\
		"movaps	0x20(%%rsi),%%xmm2	\n\t	mulpd	0x70(%%rsi),%%xmm10	\n\t"\
		"movaps	0x30(%%rsi),%%xmm3	\n\t	addpd	%%xmm11,%%xmm14		\n\t"\
		"movaps	%%xmm5,%%xmm6		\n\t	subpd	%%xmm10,%%xmm15		\n\t"\
		"movaps	%%xmm1,%%xmm7		\n\t	movaps %%xmm14,0xc0(%%r10)	\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t	movaps %%xmm15,0xd0(%%r10)	\n\t"\
		"mulpd	%%xmm3,%%xmm1		\n\t	movaps	%%xmm13,%%xmm14		\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t	movaps	%%xmm9 ,%%xmm15		\n\t"\
		"mulpd	%%xmm3,%%xmm5		\n\t	mulpd	0x80(%%rsi),%%xmm14	\n\t"\
		"movaps	0x40(%%rsi),%%xmm2	\n\t	mulpd	0x90(%%rsi),%%xmm9 	\n\t"\
		"movaps	0x50(%%rsi),%%xmm3	\n\t	mulpd	0x80(%%rsi),%%xmm15	\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t	mulpd	0x90(%%rsi),%%xmm13	\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t	addpd	%%xmm9 ,%%xmm14		\n\t"\
		"movaps	%%xmm6,0x20(%%r10)	\n\t	subpd	%%xmm13,%%xmm15		\n\t"\
		"movaps	%%xmm7,0x30(%%r10)	\n\t	movaps %%xmm14,0xa0(%%r10)	\n\t"\
		/* B3 = t3*~w3 = t3*[c3-I.s3] */"	movaps %%xmm15,0xb0(%%r10)	\n\t"\
		"movaps	%%xmm0,%%xmm6		\n\t	movaps	%%xmm8 ,%%xmm14		\n\t"\
		"movaps	%%xmm4,%%xmm7		\n\t	movaps	%%xmm12,%%xmm15		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t	mulpd	0xa0(%%rsi),%%xmm14	\n\t"\
		"mulpd	%%xmm3,%%xmm4		\n\t	mulpd	0xb0(%%rsi),%%xmm12	\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t	mulpd	0xa0(%%rsi),%%xmm15	\n\t"\
		"mulpd	%%xmm3,%%xmm0		\n\t	mulpd	0xb0(%%rsi),%%xmm8 	\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
		"subpd	%%xmm0,%%xmm7		\n\t	subpd	%%xmm8 ,%%xmm15		\n\t"\
		"movaps	%%xmm6,0x60(%%r10)	\n\t	movaps	%%xmm14,0xe0(%%r10)	\n\t"\
		"movaps	%%xmm7,0x70(%%r10)	\n\t	movaps	%%xmm15,0xf0(%%r10)	\n\t"\
	/*
		// Pass 2:
		j = 0x40;	// Set j to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's inputs
		k = p2*8;	// Set k to the bytewise address offset between the pointers to the 1st,2nd 4-DFT's outputs
		c_tmp = cc0;	// c8,4,C
		// Pass1-DFTs all use same twiddle-triplet, so lit-byte address offset between ptrs to 1st,2nd DFT's twiddles = 0
		i0 = r00; i1 = r00+8; i2 = r00+16; i3 = r00+24;
		o0 = (vec_dbl *)a; o1 = (vec_dbl *)(a+p4); o2 = (vec_dbl *)(a+p8); o3 = (vec_dbl *)(a+pC);
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,   0, o0,o1,o2,o3,k)
	*/\
		"movq	%[__two],%%rsi 		\n\t	movslq	%[__p2],%%r9		\n\t"\
		"movq	%[__add0],%%r10		\n\t	movq	%[__r0],%%rax		\n\t"\
		"movslq	%[__p4]  ,%%r11		\n\t	leaq	(%%r10,%%r11,8),%%r11	\n\t"\
		"movslq	%[__p8]  ,%%r12		\n\t	leaq	(%%r10,%%r12,8),%%r12	\n\t"\
		"movslq	%[__p12] ,%%r13		\n\t	leaq	(%%r10,%%r13,8),%%r13	\n\t"\
		"movaps	(%%rsi),%%xmm15	\n\t"/* two */\
		"movaps	     (%%rax),%%xmm0	\n\t"/* r0 */"	leaq	0x40(%%rax),%%rbx	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	0x080(%%rax),%%xmm2	\n\t"/* r0 + 8 */\
		"movaps	0x090(%%rax),%%xmm3	\n\t"\
		"movaps	0x100(%%rax),%%xmm4	\n\t"/* r0 + 16 */\
		"movaps	0x110(%%rax),%%xmm5	\n\t"\
		"movaps	0x180(%%rax),%%xmm6	\n\t"/* r0 + 24 */\
		"movaps	0x190(%%rax),%%xmm7	\n\t	leaq	0x190(%%rbx),%%r8		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t	movaps	     (%%rbx),%%xmm8 	\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t	movaps	0x010(%%rbx),%%xmm9 	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t	movaps	0x080(%%rbx),%%xmm10	\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t	movaps	0x090(%%rbx),%%xmm11	\n\t"\
		"mulpd	%%xmm15,%%xmm2		\n\t	movaps	0x100(%%rbx),%%xmm12	\n\t"\
		"mulpd	%%xmm15,%%xmm3		\n\t	movaps	0x110(%%rbx),%%xmm13	\n\t"\
		"mulpd	%%xmm15,%%xmm6		\n\t	movaps	0x180(%%rbx),%%xmm14	\n\t"\
		"mulpd	%%xmm15,%%xmm7		\n\t"/*	movaps	0x190(%%rbx),%%xmm15	\n\t"*/\
		"addpd	%%xmm0,%%xmm2		\n\t	subpd	%%xmm10,%%xmm8 		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t	subpd	%%xmm11,%%xmm9 		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t	subpd	%%xmm14,%%xmm12		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t	subpd	(%%r8),%%xmm13		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t	mulpd	%%xmm15,%%xmm10		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t	mulpd	%%xmm15,%%xmm11		\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t	mulpd	%%xmm15,%%xmm14		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t	mulpd	(%%r8),%%xmm15		\n\t"\
		"mulpd	(%%rsi),%%xmm6		\n\t	addpd	%%xmm8 ,%%xmm10		\n\t"\
		"mulpd	(%%rsi),%%xmm7		\n\t	addpd	%%xmm9 ,%%xmm11		\n\t"\
		"mulpd	(%%rsi),%%xmm5		\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
		"mulpd	(%%rsi),%%xmm4		\n\t	addpd	%%xmm13,%%xmm15		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t	subpd	%%xmm14,%%xmm10		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t	subpd	%%xmm15,%%xmm11		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t	subpd	%%xmm13,%%xmm8 		\n\t"\
		"addpd	%%xmm1,%%xmm4		\n\t	subpd	%%xmm12,%%xmm9 		\n\t"\
		"movaps	%%xmm6,    (%%r10)	\n\t	movaps	%%xmm7,0x10(%%r10)	\n\t"/* Write B0 to free up 2 regs */\
		"movaps	(%%rsi),%%xmm6		\n\t	mulpd	%%xmm6 ,%%xmm14		\n\t"\
		"									mulpd	%%xmm6 ,%%xmm15		\n\t"\
		"									mulpd	%%xmm6 ,%%xmm13		\n\t"\
		/* B2 = t0*~w2 = t0*[c2-I.s1] */"	mulpd	%%xmm6 ,%%xmm12		\n\t	movq	%[__cc0],%%rsi 	\n\t"\
		"movaps	%%xmm2,%%xmm6		\n\t	addpd	%%xmm10,%%xmm14		\n\t"\
		"movaps	%%xmm3,%%xmm7		\n\t	addpd	%%xmm11,%%xmm15		\n\t"\
		"mulpd	    (%%rsi),%%xmm6	\n\t	addpd	%%xmm8 ,%%xmm13		\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t	addpd	%%xmm9 ,%%xmm12		\n\t	shlq	$3,%%r9	\n\t"/* p2*8 */\
		"mulpd	    (%%rsi),%%xmm7	\n\t	movaps %%xmm14,    (%%r10,%%r9)	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t	movaps %%xmm15,0x10(%%r10,%%r9)	\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t	movaps	%%xmm10,%%xmm14		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t	movaps	%%xmm11,%%xmm15		\n\t"\
		"movaps	%%xmm6,    (%%r12)	\n\t	mulpd	    (%%rsi),%%xmm14	\n\t"\
		"movaps	%%xmm7,0x10(%%r12)	\n\t	mulpd	0x10(%%rsi),%%xmm11	\n\t"\
		/* B1 = t1*~w1 = t1*[c1-I.s1] */"	mulpd	    (%%rsi),%%xmm15	\n\t"\
		"movaps	0x20(%%rsi),%%xmm2	\n\t	mulpd	0x10(%%rsi),%%xmm10	\n\t"\
		"movaps	0x30(%%rsi),%%xmm3	\n\t	addpd	%%xmm11,%%xmm14		\n\t"\
		"movaps	%%xmm5,%%xmm6		\n\t	subpd	%%xmm10,%%xmm15		\n\t"\
		"movaps	%%xmm1,%%xmm7		\n\t	movaps %%xmm14,    (%%r12,%%r9)	\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t	movaps %%xmm15,0x10(%%r12,%%r9)	\n\t"\
		"mulpd	%%xmm3,%%xmm1		\n\t	movaps	%%xmm13,%%xmm14		\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t	movaps	%%xmm9 ,%%xmm15		\n\t"\
		"mulpd	%%xmm3,%%xmm5		\n\t	mulpd	0x20(%%rsi),%%xmm14	\n\t"\
		"movaps	0x40(%%rsi),%%xmm2	\n\t	mulpd	0x30(%%rsi),%%xmm9 	\n\t"\
		"movaps	0x50(%%rsi),%%xmm3	\n\t	mulpd	0x20(%%rsi),%%xmm15	\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t	mulpd	0x30(%%rsi),%%xmm13	\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t	addpd	%%xmm9 ,%%xmm14		\n\t"\
		"movaps	%%xmm6,    (%%r11)	\n\t	subpd	%%xmm13,%%xmm15		\n\t"\
		"movaps	%%xmm7,0x10(%%r11)	\n\t	movaps %%xmm14,    (%%r11,%%r9)	\n\t"\
		/* B3 = t3*~w3 = t3*[c3-I.s3] */"	movaps %%xmm15,0x10(%%r11,%%r9)	\n\t"\
		"movaps	%%xmm0,%%xmm6		\n\t	movaps	%%xmm8 ,%%xmm14		\n\t"\
		"movaps	%%xmm4,%%xmm7		\n\t	movaps	%%xmm12,%%xmm15		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t	mulpd	0x40(%%rsi),%%xmm14	\n\t"\
		"mulpd	%%xmm3,%%xmm4		\n\t	mulpd	0x50(%%rsi),%%xmm12	\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t	mulpd	0x40(%%rsi),%%xmm15	\n\t"\
		"mulpd	%%xmm3,%%xmm0		\n\t	mulpd	0x50(%%rsi),%%xmm8 	\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
		"subpd	%%xmm0,%%xmm7		\n\t	subpd	%%xmm8 ,%%xmm15		\n\t"\
		"movaps	%%xmm6,    (%%r13)	\n\t	movaps	%%xmm14,    (%%r13,%%r9)	\n\t"\
		"movaps	%%xmm7,0x10(%%r13)	\n\t	movaps	%%xmm15,0x10(%%r13,%%r9)	\n\t"\
	/*
		i0 += 2; i1 += 2; i2 += 2; i3 += 2;
		o0 = (vec_dbl *)(a+p1); o1 = (vec_dbl *)(a+p5); o2 = (vec_dbl *)(a+p9); o3 = (vec_dbl *)(a+pD);
		SSE2_RADIX_04_DIT_3TWIDDLE_X2(i0,i1,i2,i3,j, two,c_tmp,   0, o0,o1,o2,o3,k)
	*/\
		"movq	%[__two],%%rsi 		\n\t	movslq	%[__p1],%%rdi		\n\t"\
		"leaq (%%r10,%%rdi,8),%%r10	\n\t	addq	$0x20,%%rax		\n\t"\
		"leaq (%%r11,%%rdi,8),%%r11	\n\t"\
		"leaq (%%r12,%%rdi,8),%%r12	\n\t"\
		"leaq (%%r13,%%rdi,8),%%r13	\n\t"\
		"movaps	(%%rsi),%%xmm15	\n\t"/* two */\
		"movaps	     (%%rax),%%xmm0	\n\t"/* r0 */"	leaq	0x40(%%rax),%%rbx	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t"\
		"movaps	0x080(%%rax),%%xmm2	\n\t"/* r0 + 8 */\
		"movaps	0x090(%%rax),%%xmm3	\n\t"\
		"movaps	0x100(%%rax),%%xmm4	\n\t"/* r0 + 16 */\
		"movaps	0x110(%%rax),%%xmm5	\n\t"\
		"movaps	0x180(%%rax),%%xmm6	\n\t"/* r0 + 24 */\
		"movaps	0x190(%%rax),%%xmm7	\n\t	leaq	0x190(%%rbx),%%r8		\n\t"\
		"subpd	%%xmm2,%%xmm0		\n\t	movaps	     (%%rbx),%%xmm8 	\n\t"\
		"subpd	%%xmm3,%%xmm1		\n\t	movaps	0x010(%%rbx),%%xmm9 	\n\t"\
		"subpd	%%xmm6,%%xmm4		\n\t	movaps	0x080(%%rbx),%%xmm10	\n\t"\
		"subpd	%%xmm7,%%xmm5		\n\t	movaps	0x090(%%rbx),%%xmm11	\n\t"\
		"mulpd	%%xmm15,%%xmm2		\n\t	movaps	0x100(%%rbx),%%xmm12	\n\t"\
		"mulpd	%%xmm15,%%xmm3		\n\t	movaps	0x110(%%rbx),%%xmm13	\n\t"\
		"mulpd	%%xmm15,%%xmm6		\n\t	movaps	0x180(%%rbx),%%xmm14	\n\t"\
		"mulpd	%%xmm15,%%xmm7		\n\t"/*	movaps	0x190(%%rbx),%%xmm15	\n\t"*/\
		"addpd	%%xmm0,%%xmm2		\n\t	subpd	%%xmm10,%%xmm8 		\n\t"\
		"addpd	%%xmm1,%%xmm3		\n\t	subpd	%%xmm11,%%xmm9 		\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t	subpd	%%xmm14,%%xmm12		\n\t"\
		"addpd	%%xmm5,%%xmm7		\n\t	subpd	(%%r8),%%xmm13		\n\t"\
		"subpd	%%xmm6,%%xmm2		\n\t	mulpd	%%xmm15,%%xmm10		\n\t"\
		"subpd	%%xmm7,%%xmm3		\n\t	mulpd	%%xmm15,%%xmm11		\n\t"\
		"subpd	%%xmm5,%%xmm0		\n\t	mulpd	%%xmm15,%%xmm14		\n\t"\
		"subpd	%%xmm4,%%xmm1		\n\t	mulpd	(%%r8),%%xmm15		\n\t"\
		"mulpd	(%%rsi),%%xmm6		\n\t	addpd	%%xmm8 ,%%xmm10		\n\t"\
		"mulpd	(%%rsi),%%xmm7		\n\t	addpd	%%xmm9 ,%%xmm11		\n\t"\
		"mulpd	(%%rsi),%%xmm5		\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
		"mulpd	(%%rsi),%%xmm4		\n\t	addpd	%%xmm13,%%xmm15		\n\t"\
		"addpd	%%xmm2,%%xmm6		\n\t	subpd	%%xmm14,%%xmm10		\n\t"\
		"addpd	%%xmm3,%%xmm7		\n\t	subpd	%%xmm15,%%xmm11		\n\t"\
		"addpd	%%xmm0,%%xmm5		\n\t	subpd	%%xmm13,%%xmm8 		\n\t"\
		"addpd	%%xmm1,%%xmm4		\n\t	subpd	%%xmm12,%%xmm9 		\n\t"\
		"movaps	%%xmm6,    (%%r10)	\n\t	movaps	%%xmm7,0x10(%%r10)	\n\t"/* Write B0 to free up 2 regs */\
		"movaps	(%%rsi),%%xmm6		\n\t	mulpd	%%xmm6 ,%%xmm14		\n\t"\
		"									mulpd	%%xmm6 ,%%xmm15		\n\t"\
		"									mulpd	%%xmm6 ,%%xmm13		\n\t"\
		/* B2 = t0*~w2 = t0*[c2-I.s1] */"	mulpd	%%xmm6 ,%%xmm12		\n\t	movq	%[__cc0],%%rsi 	\n\t"\
		"movaps	%%xmm2,%%xmm6		\n\t	addpd	%%xmm10,%%xmm14		\n\t"\
		"movaps	%%xmm3,%%xmm7		\n\t	addpd	%%xmm11,%%xmm15		\n\t"\
		"mulpd	    (%%rsi),%%xmm6	\n\t	addpd	%%xmm8 ,%%xmm13		\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3	\n\t	addpd	%%xmm9 ,%%xmm12		\n\t"\
		"mulpd	    (%%rsi),%%xmm7	\n\t	movaps %%xmm14,    (%%r10,%%r9)	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2	\n\t	movaps %%xmm15,0x10(%%r10,%%r9)	\n\t"\
		"addpd	%%xmm3,%%xmm6		\n\t	movaps	%%xmm10,%%xmm14		\n\t"\
		"subpd	%%xmm2,%%xmm7		\n\t	movaps	%%xmm11,%%xmm15		\n\t"\
		"movaps	%%xmm6,    (%%r12)	\n\t	mulpd	    (%%rsi),%%xmm14	\n\t"\
		"movaps	%%xmm7,0x10(%%r12)	\n\t	mulpd	0x10(%%rsi),%%xmm11	\n\t"\
		/* B1 = t1*~w1 = t1*[c1-I.s1] */"	mulpd	    (%%rsi),%%xmm15	\n\t"\
		"movaps	0x20(%%rsi),%%xmm2	\n\t	mulpd	0x10(%%rsi),%%xmm10	\n\t"\
		"movaps	0x30(%%rsi),%%xmm3	\n\t	addpd	%%xmm11,%%xmm14		\n\t"\
		"movaps	%%xmm5,%%xmm6		\n\t	subpd	%%xmm10,%%xmm15		\n\t"\
		"movaps	%%xmm1,%%xmm7		\n\t	movaps %%xmm14,    (%%r12,%%r9)	\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t	movaps %%xmm15,0x10(%%r12,%%r9)	\n\t"\
		"mulpd	%%xmm3,%%xmm1		\n\t	movaps	%%xmm13,%%xmm14		\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t	movaps	%%xmm9 ,%%xmm15		\n\t"\
		"mulpd	%%xmm3,%%xmm5		\n\t	mulpd	0x20(%%rsi),%%xmm14	\n\t"\
		"movaps	0x40(%%rsi),%%xmm2	\n\t	mulpd	0x30(%%rsi),%%xmm9 	\n\t"\
		"movaps	0x50(%%rsi),%%xmm3	\n\t	mulpd	0x20(%%rsi),%%xmm15	\n\t"\
		"addpd	%%xmm1,%%xmm6		\n\t	mulpd	0x30(%%rsi),%%xmm13	\n\t"\
		"subpd	%%xmm5,%%xmm7		\n\t	addpd	%%xmm9 ,%%xmm14		\n\t"\
		"movaps	%%xmm6,    (%%r11)	\n\t	subpd	%%xmm13,%%xmm15		\n\t"\
		"movaps	%%xmm7,0x10(%%r11)	\n\t	movaps %%xmm14,    (%%r11,%%r9)	\n\t"\
		/* B3 = t3*~w3 = t3*[c3-I.s3] */"	movaps %%xmm15,0x10(%%r11,%%r9)	\n\t"\
		"movaps	%%xmm0,%%xmm6		\n\t	movaps	%%xmm8 ,%%xmm14		\n\t"\
		"movaps	%%xmm4,%%xmm7		\n\t	movaps	%%xmm12,%%xmm15		\n\t"\
		"mulpd	%%xmm2,%%xmm6		\n\t	mulpd	0x40(%%rsi),%%xmm14	\n\t"\
		"mulpd	%%xmm3,%%xmm4		\n\t	mulpd	0x50(%%rsi),%%xmm12	\n\t"\
		"mulpd	%%xmm2,%%xmm7		\n\t	mulpd	0x40(%%rsi),%%xmm15	\n\t"\
		"mulpd	%%xmm3,%%xmm0		\n\t	mulpd	0x50(%%rsi),%%xmm8 	\n\t"\
		"addpd	%%xmm4,%%xmm6		\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
		"subpd	%%xmm0,%%xmm7		\n\t	subpd	%%xmm8 ,%%xmm15		\n\t"\
		"movaps	%%xmm6,    (%%r13)	\n\t	movaps	%%xmm14,    (%%r13,%%r9)	\n\t"\
		"movaps	%%xmm7,0x10(%%r13)	\n\t	movaps	%%xmm15,0x10(%%r13,%%r9)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r0] "m" (Xr0)\
		 ,[__two] "m" (Xtwo)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_DIF_TWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xp12,Xr1,Xisrt2,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
		"movq	%[__r1],%%r9						\n\t"\
		"movq	%[__add0],%%rax						\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%rax)\n\t"/* [base-address + data-fetch-ahead index] + p0 */\
		"movslq	%[__p4],%%rbx						\n\t"\
		"movslq	%[__p8],%%rcx						\n\t	movq	%[__isrt2],%%r8			\n\t"\
		"movslq	%[__p12],%%rdx						\n\t	movslq	%[__p1],%%r10			\n\t"\
		"shlq	$3,%%rbx							\n\t	movslq	%[__p2],%%rdi			\n\t"\
		"shlq	$3,%%rcx							\n\t	shlq	$3,%%r10				\n\t"\
		"shlq	$3,%%rdx							\n\t	shlq	$3,%%rdi				\n\t"\
	"movq	%%rbx,%%r14	\n\t"	/* Save a copy of p4 pointer offset. Will prefetch from [base-address + data-fetch-ahead index] + [0,p4,p8,p12] on each macro call. */\
		"addq	%%rax,%%rbx			/* add0 + p4  */\n\t	addq	$0x230,%%r8	/* two */	\n\t"\
		"addq	%%rax,%%rcx			/* add0 + p8  */\n\t	movq	%%r10,%%r11				\n\t"\
		"addq	%%rax,%%rdx			/* add0 + p12 */\n\t	movq	%%r10,%%r12				\n\t"\
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r1 ,c0) */"	\n\t	movq	%%r10,%%r13				\n\t"\
		"/* Do	the p0,p8 combo: */					\n\t	addq	%%rax,%%r10	/* add0 + p1  */\n\t"\
		"leaq	0x230(%%r9),%%rsi	/* c0, from r1 */\n\t	addq	%%rbx,%%r11	/* add0 + p5  */\n\t"\
		"movaps	    (%%rax),%%xmm0					\n\t	addq	%%rcx,%%r12	/* add0 + p9  */\n\t"\
		"movaps	    (%%rcx),%%xmm4					\n\t	addq	%%rdx,%%r13	/* add0 + p13 */\n\t"\
		"movaps	0x10(%%rax),%%xmm1					\n\t"	/* SSE2_RADIX4_DIF_4TWIDDLE_B(r17,c1) */\
		"movaps	0x10(%%rcx),%%xmm5					\n\t	/* Do	the p0,p8 combo: */		\n\t"\
		"movaps	    (%%rsi),%%xmm6		/* c0 */	\n\t	movaps	    (%%r10),%%xmm8 		\n\t"\
		"movaps	0x10(%%rsi),%%xmm7					\n\t	movaps	    (%%r12),%%xmm12		\n\t"\
		"movaps	%%xmm0,%%xmm2						\n\t	movaps	0x10(%%r10),%%xmm9 		\n\t"\
		"movaps	%%xmm1,%%xmm3						\n\t	movaps	0x10(%%r12),%%xmm13		\n\t"\
		"mulpd	%%xmm6,%%xmm0						\n\t	movaps	0x100(%%rsi),%%xmm14	/* c1 */\n\t"\
		"mulpd	%%xmm6,%%xmm1						\n\t	movaps	0x110(%%rsi),%%xmm15	\n\t"\
		"mulpd	%%xmm7,%%xmm2						\n\t	movaps	%%xmm8 ,%%xmm10			\n\t"\
		"mulpd	%%xmm7,%%xmm3						\n\t	movaps	%%xmm9 ,%%xmm11			\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	mulpd	%%xmm14,%%xmm9 			\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	mulpd	%%xmm15,%%xmm10			\n\t"\
		"subpd	%%xmm3,%%xmm0						\n\t	mulpd	%%xmm14,%%xmm8 			\n\t"\
		"addpd	%%xmm2,%%xmm1						\n\t	mulpd	%%xmm15,%%xmm11			\n\t"\
		"mulpd	0x20(%%rsi),%%xmm4		/* c8 */	\n\t	movaps	%%xmm12,%%xmm14			\n\t"\
		"mulpd	0x30(%%rsi),%%xmm7					\n\t	movaps	%%xmm13,%%xmm15			\n\t"\
		"mulpd	0x20(%%rsi),%%xmm5					\n\t	addpd	%%xmm10,%%xmm9 			\n\t"\
		"mulpd	0x30(%%rsi),%%xmm6					\n\t	subpd	%%xmm11,%%xmm8 			\n\t"\
		"movaps	%%xmm0,%%xmm2						\n\t	mulpd	0x120(%%rsi),%%xmm12	/* c9 */\n\t"\
		"movaps	%%xmm1,%%xmm3						\n\t	mulpd	0x130(%%rsi),%%xmm15	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	mulpd	0x120(%%rsi),%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	mulpd	0x130(%%rsi),%%xmm14	\n\t"\
		"addpd	%%xmm4,%%xmm0						\n\t	movaps	%%xmm8 ,%%xmm10			\n\t"\
		"addpd	%%xmm5,%%xmm1						\n\t	movaps	%%xmm9 ,%%xmm11			\n\t"\
		"subpd	%%xmm4,%%xmm2						\n\t	subpd	%%xmm15,%%xmm12			\n\t"\
		"subpd	%%xmm5,%%xmm3						\n\t	addpd	%%xmm14,%%xmm13			\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%rbx)\n\t"/* ...+p4 */\
		"/* Do	the p4,12 combo: */					\n\t	addpd	%%xmm12,%%xmm8 			\n\t"\
		"movaps	    (%%rdx),%%xmm4					\n\t	addpd	%%xmm13,%%xmm9 			\n\t"\
		"movaps	0x10(%%rdx),%%xmm5					\n\t	subpd	%%xmm12,%%xmm10			\n\t"\
		"movaps	    (%%rdx),%%xmm6					\n\t	subpd	%%xmm13,%%xmm11			\n\t"\
		"movaps	0x10(%%rdx),%%xmm7					\n\t	/* Do	the p4,12 combo: */		\n\t"\
		"mulpd	0x60(%%rsi),%%xmm4		/* c12*/	\n\t	movaps	    (%%r13),%%xmm12		\n\t"\
		"mulpd	0x70(%%rsi),%%xmm7					\n\t	movaps	0x10(%%r13),%%xmm13		\n\t"\
		"mulpd	0x60(%%rsi),%%xmm5					\n\t	movaps	    (%%r13),%%xmm14		\n\t"\
		"mulpd	0x70(%%rsi),%%xmm6					\n\t	movaps	0x10(%%r13),%%xmm15		\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	mulpd	0x160(%%rsi),%%xmm12	/* c13*/\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	mulpd	0x170(%%rsi),%%xmm15	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7					\n\t	mulpd	0x160(%%rsi),%%xmm13	\n\t"\
		"movaps	    (%%rbx),%%xmm6					\n\t	mulpd	0x170(%%rsi),%%xmm14	\n\t"\
		"movaps	%%xmm4,     (%%r9)					\n\t	subpd	%%xmm15,%%xmm12			\n\t"\
		"movaps	%%xmm5,0x010(%%r9)					\n\t	addpd	%%xmm14,%%xmm13			\n\t"\
		"movaps	    (%%rbx),%%xmm4					\n\t	movaps	0x10(%%r11),%%xmm15		\n\t"\
		"movaps	0x10(%%rbx),%%xmm5					\n\t	movaps	    (%%r11),%%xmm14		\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r11)\n\t"/* ...+p5 */\
		"mulpd	0x40(%%rsi),%%xmm4		/* c4 */	\n\t	movaps	%%xmm12,0x100(%%r9)		\n\t"\
		"mulpd	0x50(%%rsi),%%xmm7					\n\t	movaps	%%xmm13,0x110(%%r9)		/* r17 */\n\t"\
		"mulpd	0x40(%%rsi),%%xmm5					\n\t	movaps	    (%%r11),%%xmm12		\n\t"\
		"mulpd	0x50(%%rsi),%%xmm6					\n\t	movaps	0x10(%%r11),%%xmm13		\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	mulpd	0x140(%%rsi),%%xmm12	/* c5 */\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	mulpd	0x150(%%rsi),%%xmm15	\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	mulpd	0x140(%%rsi),%%xmm13	\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	mulpd	0x150(%%rsi),%%xmm14	\n\t"\
		"subpd	     (%%r9),%%xmm4					\n\t	subpd	%%xmm15,%%xmm12			\n\t"\
		"subpd	0x010(%%r9),%%xmm5					\n\t	addpd	%%xmm14,%%xmm13			\n\t"\
		"addpd	     (%%r9),%%xmm6					\n\t	movaps	%%xmm12,%%xmm14			\n\t"\
		"addpd	0x010(%%r9),%%xmm7					\n\t	movaps	%%xmm13,%%xmm15			\n\t"\
		"/* Finish radix-4 butterfly, store: */		\n\t	subpd	0x100(%%r9),%%xmm12		\n\t"\
		"subpd	%%xmm4,%%xmm3						\n\t	subpd	0x110(%%r9),%%xmm13		\n\t"\
		"subpd	%%xmm5,%%xmm2						\n\t	addpd	0x100(%%r9),%%xmm14		\n\t"\
		"subpd	%%xmm6,%%xmm0						\n\t	addpd	0x110(%%r9),%%xmm15		\n\t"\
		"subpd	%%xmm7,%%xmm1						\n\t	/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"movaps	%%xmm3,0x070(%%r9)					\n\t	subpd	%%xmm12,%%xmm11			\n\t"\
		"movaps	%%xmm2,0x020(%%r9)					\n\t	subpd	%%xmm13,%%xmm10			\n\t"\
		"movaps	%%xmm0,0x040(%%r9)					\n\t	subpd	%%xmm14,%%xmm8 			\n\t"\
		"movaps	%%xmm1,0x050(%%r9)					\n\t	subpd	%%xmm15,%%xmm9 			\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%rcx)\n\t"/* ...+p8 */\
		"addpd	%%xmm6,%%xmm6						\n\t	movaps	%%xmm11,0x170(%%r9)		\n\t"\
		"addpd	%%xmm5,%%xmm5						\n\t	movaps	%%xmm10,0x120(%%r9)		\n\t"\
		"addpd	%%xmm7,%%xmm7						\n\t	movaps	%%xmm8 ,0x140(%%r9)		\n\t"\
		"addpd	%%xmm4,%%xmm4						\n\t	movaps	%%xmm9 ,0x150(%%r9)		\n\t"\
		"addpd	%%xmm0,%%xmm6						\n\t	mulpd	(%%r8),%%xmm14			\n\t"\
		"addpd	%%xmm2,%%xmm5						\n\t	mulpd	(%%r8),%%xmm13			\n\t"\
		"addpd	%%xmm1,%%xmm7						\n\t	mulpd	(%%r8),%%xmm15			\n\t"\
		"addpd	%%xmm3,%%xmm4						\n\t	mulpd	(%%r8),%%xmm12			\n\t"\
		"movaps	%%xmm6,     (%%r9)					\n\t	addpd	%%xmm8 ,%%xmm14			\n\t"\
		"movaps	%%xmm5,0x060(%%r9)					\n\t	addpd	%%xmm10,%%xmm13			\n\t"\
		"movaps	%%xmm7,0x010(%%r9)					\n\t	addpd	%%xmm9 ,%%xmm15			\n\t"\
		"movaps	%%xmm4,0x030(%%r9)					\n\t	addpd	%%xmm11,%%xmm12			\n\t"\
		"addq	%%rdi,%%rax		/* add0 + p2  */	\n\t	movaps	%%xmm14,0x100(%%r9)		\n\t"\
		"addq	%%rdi,%%rbx		/* add0 + p6  */	\n\t	movaps	%%xmm13,0x160(%%r9)		\n\t"\
		"addq	%%rdi,%%rcx		/* add0 + p10 */	\n\t	movaps	%%xmm15,0x110(%%r9)		\n\t"\
		"addq	%%rdi,%%rdx		/* add0 + p14 */	\n\t	movaps	%%xmm12,0x130(%%r9)		\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r12)\n\t"/* ...+p9 */\
		/* SSE2_RADIX4_DIF_4TWIDDLE_B(r9 ,c2) */"	\n\t	addq	%%rdi,%%r10	/* add0 + p3  */\n\t"\
		"/* Do	the p0,p8 combo: */					\n\t	addq	%%rdi,%%r11	/* add0 + p7  */\n\t"\
		"addq	$0x80,%%rsi 		/* c2 */		\n\t	addq	%%rdi,%%r12	/* add0 + p11 */\n\t"\
		"movaps	    (%%rax),%%xmm0					\n\t	addq	%%rdi,%%r13	/* add0 + p15 */\n\t"\
		"movaps	    (%%rcx),%%xmm4					\n\t"	/* SSE2_RADIX4_DIF_4TWIDDLE_B(r25,c3) */\
		"movaps	0x10(%%rax),%%xmm1					\n\t	/* Do	the p0,p8 combo: */		\n\t"\
		"movaps	0x10(%%rcx),%%xmm5					\n\t	movaps	    (%%r10),%%xmm8 		\n\t"\
		"movaps	    (%%rsi),%%xmm6					\n\t	movaps	    (%%r12),%%xmm12		\n\t"\
		"movaps	0x10(%%rsi),%%xmm7					\n\t	movaps	0x10(%%r10),%%xmm9 		\n\t"\
		"movaps	%%xmm0,%%xmm2						\n\t	movaps	0x10(%%r12),%%xmm13		\n\t"\
		"movaps	%%xmm1,%%xmm3						\n\t	movaps	0x100(%%rsi),%%xmm14	/* c3 */\n\t"\
		"mulpd	%%xmm6,%%xmm0						\n\t	movaps	0x110(%%rsi),%%xmm15	\n\t"\
		"mulpd	%%xmm7,%%xmm3						\n\t	movaps	%%xmm8 ,%%xmm10			\n\t"\
		"mulpd	%%xmm6,%%xmm1						\n\t	movaps	%%xmm9 ,%%xmm11			\n\t"\
		"mulpd	%%xmm7,%%xmm2						\n\t	mulpd	 %%xmm14,%%xmm8 		\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	mulpd	 %%xmm15,%%xmm11		\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	mulpd	 %%xmm14,%%xmm9 		\n\t"\
		"subpd	%%xmm3,%%xmm0						\n\t	mulpd	 %%xmm15,%%xmm10		\n\t"\
		"addpd	%%xmm2,%%xmm1						\n\t	movaps	%%xmm12,%%xmm14			\n\t"\
		"mulpd	0x20(%%rsi),%%xmm4	/* c10*/		\n\t	movaps	%%xmm13,%%xmm15			\n\t"\
"prefetcht1	%c[__pfetch_dist](%%rax)\n\t"/*..+p2 */\
		"mulpd	0x30(%%rsi),%%xmm7					\n\t	subpd	 %%xmm11,%%xmm8 		\n\t"\
		"mulpd	0x20(%%rsi),%%xmm5					\n\t	addpd	 %%xmm10,%%xmm9 		\n\t"\
		"mulpd	0x30(%%rsi),%%xmm6					\n\t	mulpd	0x120(%%rsi),%%xmm12	/* c11*/\n\t"\
		"movaps	%%xmm0,%%xmm2						\n\t	mulpd	0x130(%%rsi),%%xmm15	\n\t"\
		"movaps	%%xmm1,%%xmm3						\n\t	mulpd	0x120(%%rsi),%%xmm13	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	mulpd	0x130(%%rsi),%%xmm14	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	movaps	%%xmm8 ,%%xmm10			\n\t"\
		"addpd	%%xmm4,%%xmm0						\n\t	movaps	%%xmm9 ,%%xmm11			\n\t"\
		"addpd	%%xmm5,%%xmm1						\n\t	subpd	%%xmm15,%%xmm12			\n\t"\
		"subpd	%%xmm4,%%xmm2						\n\t	addpd	%%xmm14,%%xmm13			\n\t"\
		"subpd	%%xmm5,%%xmm3						\n\t	addpd	%%xmm12,%%xmm8 			\n\t"\
		"/* Do	the p4,12 combo: */					\n\t	addpd	%%xmm13,%%xmm9 			\n\t"\
		"addq	$0x80,%%r9 			/* r9 */		\n\t	subpd	%%xmm12,%%xmm10			\n\t"\
		"movaps	    (%%rdx),%%xmm4					\n\t	subpd	%%xmm13,%%xmm11			\n\t"\
		"movaps	0x10(%%rdx),%%xmm5					\n\t	/* Do	the p4,12 combo: */		\n\t"\
		"movaps	    (%%rdx),%%xmm6					\n\t	movaps	    (%%r13),%%xmm12		\n\t"\
		"movaps	0x10(%%rdx),%%xmm7					\n\t	movaps	0x10(%%r13),%%xmm13		\n\t"\
		"mulpd	0x60(%%rsi),%%xmm4	/* c14*/		\n\t	movaps	    (%%r13),%%xmm14		\n\t"\
		"mulpd	0x70(%%rsi),%%xmm7					\n\t	movaps	0x10(%%r13),%%xmm15		\n\t"\
"prefetcht1	%c[__pfetch_dist](%%r10)\n\t"/*..+p3 */\
		"mulpd	0x60(%%rsi),%%xmm5					\n\t	mulpd	0x160(%%rsi),%%xmm12	/* c15*/\n\t"\
		"mulpd	0x70(%%rsi),%%xmm6					\n\t	mulpd	0x170(%%rsi),%%xmm15	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	mulpd	0x160(%%rsi),%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	mulpd	0x170(%%rsi),%%xmm14	\n\t"\
		"movaps	0x10(%%rbx),%%xmm7					\n\t	subpd	%%xmm15,%%xmm12			\n\t"\
		"movaps	    (%%rbx),%%xmm6					\n\t	addpd	%%xmm14,%%xmm13			\n\t"\
		"movaps	%%xmm4,     (%%r9)					\n\t	movaps	%%xmm12,0x100(%%r9)		\n\t"\
		"movaps	%%xmm5,0x010(%%r9)					\n\t	movaps	%%xmm13,0x110(%%r9)		/* r25*/\n\t"\
		"movaps	    (%%rbx),%%xmm4					\n\t	movaps	    (%%r11),%%xmm12		\n\t"\
		"movaps	0x10(%%rbx),%%xmm5					\n\t	movaps	0x10(%%r11),%%xmm13		\n\t"\
		"mulpd	0x40(%%rsi),%%xmm4	/* c6 */		\n\t	movaps	    (%%r11),%%xmm14		\n\t"\
		"mulpd	0x50(%%rsi),%%xmm7					\n\t	movaps	0x10(%%r11),%%xmm15		\n\t"\
		"mulpd	0x40(%%rsi),%%xmm5					\n\t	mulpd	0x140(%%rsi),%%xmm12	/* c7 */\n\t"\
		"mulpd	0x50(%%rsi),%%xmm6					\n\t	mulpd	0x150(%%rsi),%%xmm15	\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	mulpd	0x140(%%rsi),%%xmm13	\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	mulpd	0x150(%%rsi),%%xmm14	\n\t"\
		"movaps	%%xmm4,%%xmm6						\n\t	subpd	%%xmm15,%%xmm12			\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	addpd	%%xmm14,%%xmm13			\n\t"\
		"subpd	     (%%r9),%%xmm4					\n\t	movaps	%%xmm12,%%xmm14			\n\t"\
		"subpd	0x010(%%r9),%%xmm5					\n\t	movaps	%%xmm13,%%xmm15			\n\t"\
		"addpd	     (%%r9),%%xmm6					\n\t	subpd	0x100(%%r9),%%xmm12		\n\t"\
		"addpd	0x010(%%r9),%%xmm7					\n\t	subpd	0x110(%%r9),%%xmm13		\n\t"\
"prefetcht1	%c[__pfetch_dist](%%rbx)\n\t"/*..+p6 */\
		"/* Finish radix-4 butterfly, store: */		\n\t	addpd	0x100(%%r9),%%xmm14		\n\t"\
		"subpd	%%xmm4,%%xmm3						\n\t	addpd	0x110(%%r9),%%xmm15		\n\t"\
		"subpd	%%xmm5,%%xmm2						\n\t	/* Finish radix-4 butterfly and store results into temporary-array slots: */\n\t"\
		"subpd	%%xmm6,%%xmm0						\n\t	subpd	%%xmm12,%%xmm11			\n\t"\
		"subpd	%%xmm7,%%xmm1						\n\t	subpd	%%xmm13,%%xmm10			\n\t"\
		"movaps	%%xmm3,0x070(%%r9)					\n\t	subpd	%%xmm14,%%xmm8 			\n\t"\
		"movaps	%%xmm2,0x020(%%r9)					\n\t	subpd	%%xmm15,%%xmm9 			\n\t"\
		"movaps	%%xmm0,0x040(%%r9)					\n\t	movaps	%%xmm11,0x170(%%r9)		\n\t"\
		"movaps	%%xmm1,0x050(%%r9)					\n\t	movaps	%%xmm10,0x120(%%r9)		\n\t"\
		"addpd	%%xmm6,%%xmm6						\n\t	movaps	%%xmm8 ,0x140(%%r9)		\n\t"\
		"addpd	%%xmm5,%%xmm5						\n\t	movaps	%%xmm9 ,0x150(%%r9)		\n\t"\
		"addpd	%%xmm7,%%xmm7						\n\t	mulpd	(%%r8),%%xmm14			\n\t"\
		"addpd	%%xmm4,%%xmm4						\n\t	mulpd	(%%r8),%%xmm13			\n\t"\
		"addpd	%%xmm0,%%xmm6						\n\t	mulpd	(%%r8),%%xmm15			\n\t"\
		"addpd	%%xmm2,%%xmm5						\n\t	mulpd	(%%r8),%%xmm12			\n\t"\
		"addpd	%%xmm1,%%xmm7						\n\t	addpd	%%xmm8 ,%%xmm14			\n\t"\
		"addpd	%%xmm3,%%xmm4						\n\t	addpd	%%xmm10,%%xmm13			\n\t"\
		"movaps	%%xmm6,     (%%r9)					\n\t	addpd	%%xmm9 ,%%xmm15			\n\t"\
		"movaps	%%xmm5,0x060(%%r9)					\n\t	addpd	%%xmm11,%%xmm12			\n\t"\
		"movaps	%%xmm7,0x010(%%r9)					\n\t	movaps	%%xmm14,0x100(%%r9)		\n\t"\
		"movaps	%%xmm4,0x030(%%r9)					\n\t	movaps	%%xmm13,0x160(%%r9)		\n\t"\
		"											\n\t	movaps	%%xmm15,0x110(%%r9)		\n\t"\
		"											\n\t	movaps	%%xmm12,0x130(%%r9)		\n\t"\
		/*************************************************************************************/\
		/*  And now do four more radix-4 transforms, including the internal twiddle factors: */\
		/*************************************************************************************/\
		"/* Block 1 - p2 in rdi, add0+p2,3 already in rax,r10: */	\n\t"\
		"movq	%%r10,%%rbx		/* cpy add0+p3*/	\n\t"\
		"movq	%%rax,%%rcx		/* add0+p2  */		\n\t"\
		"movq	%%r10,%%rdx		/* add0+p3  */		/* Block 3: */				\n\t"\
		"subq	%%rdi,%%rax		/* add0+p0  */		\n\t	movslq	%[__p4],%%r10			\n\t"\
		"subq	%%rdi,%%rbx		/* add0+p1  */		\n\t	movslq	%[__p8],%%rdi			\n\t"\
		"movq	%[__isrt2],%%rsi					\n\t	shlq	$3,%%r10				\n\t"\
		"subq	$0x80,%%r9		/* r1 */			\n\t	shlq	$3,%%rdi				\n\t"\
		"movaps	     (%%r9),%%xmm0					\n\t	movq	%%r10,%%r11				\n\t"\
		"movaps	0x100(%%r9),%%xmm4					\n\t	movq	%%r10,%%r12				\n\t"\
		"movaps	0x010(%%r9),%%xmm1					\n\t	movq	%%r10,%%r13				\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r11)\n\t"/* ...+p7 */\
		"movaps	0x110(%%r9),%%xmm5					\n\t	addq	%%rax,%%r10	/* add0 + p4  */\n\t"\
		"movaps	0x080(%%r9),%%xmm2					\n\t	addq	%%rbx,%%r11	/* add0 + p5  */\n\t"\
		"movaps	0x180(%%r9),%%xmm6					\n\t	addq	%%rcx,%%r12	/* add0 + p6  */\n\t"\
		"movaps	0x090(%%r9),%%xmm3					\n\t	addq	%%rdx,%%r13	/* add0 + p7  */\n\t"\
		"movaps	0x190(%%r9),%%xmm7					\n\t	movaps	(%%rsi),%%xmm11	/* isrt2 */	\n\t"\
		"subpd	%%xmm2,%%xmm0						\n\t	movaps	0x140(%%r9),%%xmm12	/* r5 */\n\t"\
		"subpd	%%xmm6,%%xmm4						\n\t	movaps	0x150(%%r9),%%xmm13	\n\t"\
		"subpd	%%xmm3,%%xmm1						\n\t	movaps	0x1c0(%%r9),%%xmm14	\n\t"\
		"subpd	%%xmm7,%%xmm5						\n\t	movaps	0x1d0(%%r9),%%xmm15	\n\t"\
		"addpd	%%xmm2,%%xmm2						\n\t	mulpd	%%xmm11,%%xmm12		\n\t"\
		"addpd	%%xmm6,%%xmm6						\n\t	mulpd	%%xmm11,%%xmm13		\n\t"\
		"addpd	%%xmm3,%%xmm3						\n\t	mulpd	%%xmm11,%%xmm14		\n\t"\
		"addpd	%%xmm7,%%xmm7						\n\t	mulpd	%%xmm11,%%xmm15		\n\t"\
		"addpd	%%xmm0,%%xmm2						\n\t	movaps	0x040(%%r9),%%xmm8 	\n\t"\
		"addpd	%%xmm4,%%xmm6						\n\t	movaps	0x0d0(%%r9),%%xmm11	\n\t"\
		"addpd	%%xmm1,%%xmm3						\n\t	movaps	0x050(%%r9),%%xmm9 	\n\t"\
		"addpd	%%xmm5,%%xmm7						\n\t	movaps	0x0c0(%%r9),%%xmm10	\n\t"\
		"subpd	%%xmm6,%%xmm2						\n\t	subpd	%%xmm13,%%xmm12		\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r12,%%r14)\n\t"/*..+p10 */\
		"subpd	%%xmm7,%%xmm3						\n\t	subpd	%%xmm14,%%xmm15		\n\t"\
		"mulpd	(%%r8),%%xmm6						\n\t	subpd	%%xmm11,%%xmm8 		\n\t"\
		"mulpd	(%%r8),%%xmm7						\n\t	subpd	%%xmm10,%%xmm9 		\n\t"\
		"movaps	%%xmm2,    (%%rbx)					\n\t	addpd	%%xmm11,%%xmm11		\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)					\n\t	addpd	%%xmm13,%%xmm13		\n\t"\
		"addpd	%%xmm2,	%%xmm6						\n\t	mulpd	(%%r8),%%xmm10		\n\t"\
		"addpd	%%xmm3,	%%xmm7						\n\t	mulpd	(%%r8),%%xmm14		\n\t"\
		"movaps	%%xmm6,    (%%rax)					\n\t	addpd	%%xmm8 ,%%xmm11		\n\t"\
		"movaps	%%xmm7,0x10(%%rax)					\n\t	addpd	%%xmm12,%%xmm13		\n\t"\
		"subpd	%%xmm5,	%%xmm0						\n\t	addpd	%%xmm9 ,%%xmm10		\n\t"\
		"subpd	%%xmm4,	%%xmm1						\n\t	addpd	%%xmm15,%%xmm14		\n\t"\
		"mulpd	(%%r8),	%%xmm5						\n\t	subpd	%%xmm15,%%xmm13		\n\t"\
		"mulpd	(%%r8),	%%xmm4						\n\t	addpd	%%xmm15,%%xmm15		\n\t"\
		"movaps	%%xmm0,    (%%rcx)					\n\t	subpd	%%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)					\n\t	addpd	%%xmm14,%%xmm14		\n\t"\
		"addpd	%%xmm0,	%%xmm5						\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
		"addpd	%%xmm1,	%%xmm4						\n\t	addpd	%%xmm13,%%xmm15		\n\t"\
		"movaps	%%xmm5,    (%%rdx)					\n\t	subpd	%%xmm12,%%xmm8 		\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)					\n\t	subpd	%%xmm15,%%xmm11		\n\t"\
		"											\n\t	subpd	%%xmm13,%%xmm10		\n\t"\
		"											\n\t	subpd	%%xmm14,%%xmm9 		\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r13,%%r14)\n\t"/* ...+p11 */\
		"											\n\t	addpd	%%xmm12,%%xmm12		\n\t"\
		"											\n\t	addpd	%%xmm15,%%xmm15		\n\t"\
		"											\n\t	addpd	%%xmm13,%%xmm13		\n\t"\
		"											\n\t	addpd	%%xmm14,%%xmm14		\n\t"\
		"											\n\t	movaps	%%xmm8 ,    (%%r11)	\n\t"\
		"											\n\t	movaps	%%xmm11,    (%%r12)	\n\t"\
		"											\n\t	movaps	%%xmm10,0x10(%%r11)	\n\t"\
		"											\n\t	movaps	%%xmm9 ,0x10(%%r13)	\n\t"\
		"/* Block 2: */								\n\t	addpd	%%xmm8 ,	%%xmm12	\n\t"\
		"addq	%%rdi,%%rax			/* add0 + p8  */\n\t	addpd	%%xmm11,	%%xmm15	\n\t"\
		"addq	%%rdi,%%rbx			/* add0 + p9  */\n\t	addpd	%%xmm10,	%%xmm13	\n\t"\
		"addq	%%rdi,%%rcx			/* add0 + p10 */\n\t	addpd	%%xmm9 ,	%%xmm14	\n\t"\
		"addq	%%rdi,%%rdx			/* add0 + p11 */\n\t	movaps	%%xmm12,    (%%r10)	\n\t"\
		"addq	$0x20,%%r9	/* r3 */				\n\t	movaps	%%xmm15,    (%%r13)	\n\t"\
		"movaps	0x100(%%r9),%%xmm4					\n\t	movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"movaps	0x110(%%r9),%%xmm5					\n\t	movaps	%%xmm14,0x10(%%r12)	\n\t"\
		"/* Share cc0/ss0 between 2 halves: */		\n\t	/* Block 4: */				\n\t"\
		"movaps	0x10(%%rsi),%%xmm11	/* cc0 */		\n\t	addq	%%rdi,%%r10	/* add0 + p12 */\n\t"\
		"movaps	0x20(%%rsi),%%xmm10	/* ss0 */		\n\t	addq	%%rdi,%%r11	/* add0 + p13 */\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r10)\n\t"/* ...+p12 */\
		"movaps	%%xmm4,%%xmm6						\n\t	addq	%%rdi,%%r12	/* add0 + p14 */\n\t"\
		"movaps	%%xmm5,%%xmm7						\n\t	addq	%%rdi,%%r13	/* add0 + p15 */\n\t"\
		"mulpd	%%xmm11,%%xmm4						\n\t	movaps	0x140(%%r9),%%xmm12	/* r7 */\n\t"\
		"mulpd	%%xmm10,%%xmm7						\n\t	movaps	0x150(%%r9),%%xmm13	\n\t"\
		"mulpd	%%xmm11,%%xmm5						\n\t	movaps	%%xmm12,%%xmm14		\n\t"\
		"mulpd	%%xmm10,%%xmm6						\n\t	movaps	%%xmm13,%%xmm15		\n\t"\
		"movaps	0x180(%%r9),%%xmm0					\n\t	mulpd	%%xmm10,%%xmm12		\n\t"\
		"movaps	0x190(%%r9),%%xmm1					\n\t	mulpd	%%xmm11,%%xmm15		\n\t"\
		"subpd	%%xmm7,%%xmm4						\n\t	mulpd	%%xmm10,%%xmm13		\n\t"\
		"addpd	%%xmm6,%%xmm5						\n\t	mulpd	%%xmm11,%%xmm14		\n\t"\
		"movaps	%%xmm1,%%xmm7						\n\t	movaps	0x1c0(%%r9),%%xmm8 	\n\t"\
		"movaps	%%xmm0,%%xmm6						\n\t	movaps	0x1d0(%%r9),%%xmm9 	\n\t"\
		"mulpd	%%xmm11,%%xmm0						\n\t	subpd	%%xmm15,%%xmm12		\n\t"\
		"mulpd	%%xmm10,%%xmm7						\n\t	addpd	%%xmm14,%%xmm13		\n\t"\
		"mulpd	%%xmm10,%%xmm6						\n\t	movaps	%%xmm8 ,%%xmm14		\n\t"\
		"mulpd	%%xmm11,%%xmm1						\n\t	movaps	%%xmm9 ,%%xmm15		\n\t"\
		"addpd	%%xmm0,%%xmm7						\n\t	mulpd	%%xmm11,%%xmm15		\n\t"\
		"subpd	%%xmm1,%%xmm6						\n\t	mulpd	%%xmm10,%%xmm8 		\n\t"\
		"movaps	%%xmm4,%%xmm2						\n\t	mulpd	%%xmm11,%%xmm14		\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r11)\n\t"/*..+p13*/\
		"movaps	%%xmm5,%%xmm3						\n\t	mulpd	%%xmm10,%%xmm9 		\n\t"\
		"subpd	%%xmm7,%%xmm5						\n\t	movaps	%%xmm12,%%xmm10		\n\t"\
		"subpd	%%xmm6,%%xmm4						\n\t	movaps	%%xmm13,%%xmm11		\n\t"\
		"addpd	%%xmm3,%%xmm7						\n\t	addpd	%%xmm8 ,%%xmm15		\n\t"\
		"addpd	%%xmm2,%%xmm6						\n\t	subpd	%%xmm9 ,%%xmm14		\n\t"\
		"movaps	0x090(%%r9),%%xmm3					\n\t	subpd	%%xmm15,%%xmm13		\n\t"\
		"movaps	0x080(%%r9),%%xmm2					\n\t	subpd	%%xmm14,%%xmm12		\n\t"\
		"movaps	(%%rsi),%%xmm1	/* isrt2 */			\n\t	addpd	%%xmm11,%%xmm15		\n\t"\
		"movaps	%%xmm2,%%xmm0						\n\t	addpd	%%xmm10,%%xmm14		\n\t"\
		"subpd	%%xmm3,%%xmm2						\n\t	movaps	0x0c0(%%r9),%%xmm10	\n\t"\
		"addpd	%%xmm0,%%xmm3						\n\t	movaps	0x0d0(%%r9),%%xmm11	\n\t"\
		"mulpd	%%xmm1,%%xmm2						\n\t	movaps	(%%rsi),%%xmm9 	/* isrt2 */\n\t"\
		"mulpd	%%xmm1,%%xmm3						\n\t	movaps	%%xmm10,%%xmm8 		\n\t"\
		"movaps	     (%%r9),%%xmm0					\n\t	addpd	%%xmm11,%%xmm10		\n\t"\
		"movaps	0x010(%%r9),%%xmm1					\n\t	subpd	%%xmm8 ,%%xmm11		\n\t"\
		"subpd	%%xmm2,%%xmm0						\n\t	mulpd	%%xmm9 ,%%xmm10		\n\t"\
		"subpd	%%xmm3,%%xmm1						\n\t	mulpd	%%xmm9 ,%%xmm11		\n\t"\
		"addpd	%%xmm2,%%xmm2						\n\t	movaps	0x040(%%r9),%%xmm8 	\n\t"\
		"addpd	%%xmm3,%%xmm3						\n\t	movaps	0x050(%%r9),%%xmm9 	\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r12)\n\t"/* ...+p14 */\
		"addpd	%%xmm0,%%xmm2						\n\t	subpd	%%xmm10,%%xmm8 		\n\t"\
		"addpd	%%xmm1,%%xmm3						\n\t	addpd	%%xmm10,%%xmm10		\n\t"\
		"subpd	%%xmm5,%%xmm0						\n\t	subpd	%%xmm11,%%xmm9 		\n\t"\
		"subpd	%%xmm4,%%xmm1						\n\t	addpd	%%xmm11,%%xmm11		\n\t"\
		"subpd	%%xmm6,%%xmm2						\n\t	addpd	%%xmm8 ,%%xmm10		\n\t"\
		"subpd	%%xmm7,%%xmm3						\n\t	addpd	%%xmm9 ,%%xmm11		\n\t"\
		"mulpd	(%%r8),%%xmm6						\n\t	subpd	%%xmm12,%%xmm8 		\n\t"\
		"mulpd	(%%r8),%%xmm5						\n\t	subpd	%%xmm13,%%xmm9 		\n\t"\
		"mulpd	(%%r8),%%xmm7						\n\t	subpd	%%xmm15,%%xmm10		\n\t"\
		"mulpd	(%%r8),%%xmm4						\n\t	subpd	%%xmm14,%%xmm11		\n\t"\
		"movaps	%%xmm2,    (%%rbx)					\n\t	addpd	%%xmm12,%%xmm12		\n\t"\
		"movaps	%%xmm0,    (%%rcx)					\n\t	addpd	%%xmm15,%%xmm15		\n\t"\
		"movaps	%%xmm3,0x10(%%rbx)					\n\t	addpd	%%xmm13,%%xmm13		\n\t"\
		"movaps	%%xmm1,0x10(%%rdx)					\n\t	addpd	%%xmm14,%%xmm14		\n\t"\
		"addpd	%%xmm2,	%%xmm6						\n\t	movaps	%%xmm8 ,    (%%r11)	\n\t"\
		"addpd	%%xmm0,	%%xmm5						\n\t	movaps	%%xmm10,    (%%r12)	\n\t"\
		"addpd	%%xmm3,	%%xmm7						\n\t	movaps	%%xmm9 ,0x10(%%r11)	\n\t"\
		"addpd	%%xmm1,	%%xmm4						\n\t	movaps	%%xmm11,0x10(%%r13)	\n\t"\
		"movaps	%%xmm6,    (%%rax)					\n\t	addpd	%%xmm8 ,%%xmm12		\n\t"\
		"movaps	%%xmm5,    (%%rdx)					\n\t	addpd	%%xmm10,%%xmm15		\n\t"\
		"movaps	%%xmm7,0x10(%%rax)					\n\t	addpd	%%xmm9 ,%%xmm13		\n\t"\
		"movaps	%%xmm4,0x10(%%rcx)					\n\t	addpd	%%xmm11,%%xmm14		\n\t"\
	"prefetcht1	%c[__pfetch_dist](%%r13)\n\t"/* ...+p15 */\
		"													movaps	%%xmm12,    (%%r10)	\n\t"\
		"													movaps	%%xmm15,    (%%r13)	\n\t"\
		"													movaps	%%xmm13,0x10(%%r10)	\n\t"\
		"													movaps	%%xmm14,0x10(%%r12)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__p12] "m" (Xp12)\
		 ,[__r1] "m" (Xr1)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r9","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX16_DIT_TWIDDLE(Xadd0,Xp1,Xp2,Xp3,Xp4,Xp8,Xr1,Xisrt2,Xpfetch_addr1,Xpfetch_dist)\
	{\
	__asm__ volatile (\
		"movq	%[__isrt2],%%r8		\n\t"\
		"movq	%[__add0],%%rax		\n\t"\
		"movslq	%[__p1],%%rbx		\n\t"\
		"movslq	%[__p2],%%rcx		\n\t"\
		"movslq	%[__p3],%%rdx		\n\t"\
		"movslq	%[__p4],%%rdi		\n\t"\
		"shlq	$3,%%rbx			\n\t"\
		"shlq	$3,%%rcx			\n\t"\
		"shlq	$3,%%rdx			\n\t"\
		"shlq	$3,%%rdi			\n\t"\
		"addq	$0x230,%%r8	/* two */\n\t"\
		"addq	%%rax,%%rbx			\n\t"\
		"addq	%%rax,%%rcx			\n\t"\
		"addq	%%rax,%%rdx			\n\t			movslq	%[__p8],%%r10		\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r1) */\n\t	shlq	$3,%%r10			\n\t"\
		"movq	%[__r1],%%rsi			\n\t		movq	%%r10,%%r11			\n\t"\
		"movaps	    (%%rax),%%xmm2		\n\t		movq	%%r10,%%r12			\n\t"\
		"movaps	    (%%rcx),%%xmm6		\n\t		movq	%%r10,%%r13			\n\t"\
		"movaps	0x10(%%rax),%%xmm3		\n\t		addq	%%rax,%%r10			\n\t"\
		"movaps	0x10(%%rcx),%%xmm7		\n\t		addq	%%rbx,%%r11			\n\t"\
		"movaps	    (%%rbx),%%xmm0		\n\t		addq	%%rcx,%%r12			\n\t"\
		"movaps	    (%%rdx),%%xmm4		\n\t		addq	%%rdx,%%r13			\n\t"\
		"movaps	0x10(%%rbx),%%xmm1		\n\t		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r17) */\n\t"\
		"movaps	0x10(%%rdx),%%xmm5		\n\t		leaq	0x100(%%rsi),%%r14	/* r17 */\n\t"\
		"subpd	%%xmm0,%%xmm2			\n\t		movaps	    (%%r10),%%xmm10	\n\t"\
		"subpd	%%xmm4,%%xmm6			\n\t		movaps	    (%%r12),%%xmm14	\n\t"\
		"subpd	%%xmm1,%%xmm3			\n\t		movaps	0x10(%%r10),%%xmm11	\n\t"\
		"subpd	%%xmm5,%%xmm7			\n\t		movaps	0x10(%%r12),%%xmm15	\n\t"\
		"addpd	%%xmm0,%%xmm0			\n\t		movaps	    (%%r11),%%xmm8	\n\t"\
		"addpd	%%xmm4,%%xmm4			\n\t		movaps	    (%%r13),%%xmm12	\n\t"\
		"addpd	%%xmm1,%%xmm1			\n\t		movaps	0x10(%%r11),%%xmm9	\n\t"\
		"addpd	%%xmm5,%%xmm5			\n\t		movaps	0x10(%%r13),%%xmm13	\n\t"\
		"addpd	%%xmm2,%%xmm0			\n\t		subpd	%%xmm8 ,%%xmm10		\n\t"\
		"addpd	%%xmm6,%%xmm4			\n\t		subpd	%%xmm12,%%xmm14		\n\t"\
		"addpd	%%xmm3,%%xmm1			\n\t		subpd	%%xmm9 ,%%xmm11		\n\t"\
		"addpd	%%xmm7,%%xmm5			\n\t		subpd	%%xmm13,%%xmm15		\n\t"\
		"subpd	%%xmm4,%%xmm0			\n\t		mulpd	(%%r8),%%xmm8		\n\t"\
		"subpd	%%xmm7,%%xmm2			\n\t		mulpd	(%%r8),%%xmm12		\n\t"\
		"subpd	%%xmm5,%%xmm1			\n\t		mulpd	(%%r8),%%xmm9		\n\t"\
		"subpd	%%xmm6,%%xmm3			\n\t		mulpd	(%%r8),%%xmm13		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)		\n\t		addpd	%%xmm10,%%xmm8		\n\t"\
		"movaps	%%xmm2,0x060(%%rsi)		\n\t		addpd	%%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)		\n\t		addpd	%%xmm11,%%xmm9		\n\t"\
		"movaps	%%xmm3,0x030(%%rsi)		\n\t		addpd	%%xmm15,%%xmm13		\n\t"\
		"mulpd	(%%r8),%%xmm4			\n\t		subpd	%%xmm12,%%xmm8		\n\t"\
		"mulpd	(%%r8),%%xmm7			\n\t		subpd	%%xmm15,%%xmm10		\n\t"\
		"mulpd	(%%r8),%%xmm5			\n\t		subpd	%%xmm13,%%xmm9		\n\t"\
		"mulpd	(%%r8),%%xmm6			\n\t		subpd	%%xmm14,%%xmm11		\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t		movaps	%%xmm8 ,0x040(%%r14)\n\t"\
		"addpd	%%xmm2,%%xmm7			\n\t		movaps	%%xmm10,0x060(%%r14)\n\t"\
		"addpd	%%xmm1,%%xmm5			\n\t		movaps	%%xmm9 ,0x050(%%r14)\n\t"\
		"addpd	%%xmm3,%%xmm6			\n\t		movaps	%%xmm11,0x030(%%r14)\n\t"\
		"movaps	%%xmm4,     (%%rsi)		\n\t		addpd	%%xmm12,%%xmm12		\n\t"\
		"movaps	%%xmm7,0x020(%%rsi)		\n\t		addpd	%%xmm15,%%xmm15		\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)		\n\t		addpd	%%xmm13,%%xmm13		\n\t"\
		"movaps	%%xmm6,0x070(%%rsi)		\n\t		addpd	%%xmm14,%%xmm14		\n\t"\
		"addq	%%rdi,%%rax				\n\t		addpd	%%xmm8,%%xmm12		\n\t"\
		"addq	%%rdi,%%rbx				\n\t		addpd	%%xmm10,%%xmm15		\n\t"\
		"addq	%%rdi,%%rcx				\n\t		addpd	%%xmm9,%%xmm13		\n\t"\
		"addq	%%rdi,%%rdx				\n\t		addpd	%%xmm11,%%xmm14		\n\t"\
		"/* SSE2_RADIX4_DIT_0TWIDDLE_B(r9) */\n\t	movaps	%%xmm12,     (%%r14)\n\t"\
		"addq	$0x80,%%rsi		/* r9 */\n\t		movaps	%%xmm15,0x020(%%r14)\n\t"\
		"movaps	    (%%rax),%%xmm2		\n\t		movaps	%%xmm13,0x010(%%r14)\n\t"\
		"movaps	    (%%rcx),%%xmm6		\n\t		movaps	%%xmm14,0x070(%%r14)\n\t"\
		"movaps	0x10(%%rax),%%xmm3		\n\t		addq	%%rdi,%%r10			\n\t"\
		"movaps	0x10(%%rcx),%%xmm7		\n\t		addq	%%rdi,%%r11			\n\t"\
		"movaps	    (%%rbx),%%xmm0		\n\t		addq	%%rdi,%%r12			\n\t"\
		"movaps	    (%%rdx),%%xmm4		\n\t		addq	%%rdi,%%r13			\n\t"\
		"movaps	0x10(%%rbx),%%xmm1		\n\t		/* SSE2_RADIX4_DIT_0TWIDDLE_B(r25) */\n\t"\
		"movaps	0x10(%%rdx),%%xmm5		\n\t		leaq	0x100(%%rsi),%%r14	/* r25 */\n\t"\
		"subpd	%%xmm0,%%xmm2			\n\t		movaps	    (%%r10),%%xmm10	\n\t"\
		"subpd	%%xmm4,%%xmm6			\n\t		movaps	    (%%r12),%%xmm14	\n\t"\
		"subpd	%%xmm1,%%xmm3			\n\t		movaps	0x10(%%r10),%%xmm11	\n\t"\
		"subpd	%%xmm5,%%xmm7			\n\t		movaps	0x10(%%r12),%%xmm15	\n\t"\
		"addpd	%%xmm0,%%xmm0			\n\t		movaps	    (%%r11),%%xmm8	\n\t"\
		"addpd	%%xmm4,%%xmm4			\n\t		movaps	    (%%r13),%%xmm12	\n\t"\
		"addpd	%%xmm1,%%xmm1			\n\t		movaps	0x10(%%r11),%%xmm9	\n\t"\
		"addpd	%%xmm5,%%xmm5			\n\t		movaps	0x10(%%r13),%%xmm13	\n\t"\
		"addpd	%%xmm2,%%xmm0			\n\t		subpd	%%xmm8 ,%%xmm10		\n\t"\
		"addpd	%%xmm6,%%xmm4			\n\t		subpd	%%xmm12,%%xmm14		\n\t"\
		"addpd	%%xmm3,%%xmm1			\n\t		subpd	%%xmm9 ,%%xmm11		\n\t"\
		"addpd	%%xmm7,%%xmm5			\n\t		subpd	%%xmm13,%%xmm15		\n\t"\
		"subpd	%%xmm4,%%xmm0			\n\t		mulpd	(%%r8),%%xmm8		\n\t"\
		"subpd	%%xmm7,%%xmm2			\n\t		mulpd	(%%r8),%%xmm12		\n\t"\
		"subpd	%%xmm5,%%xmm1			\n\t		mulpd	(%%r8),%%xmm9		\n\t"\
		"subpd	%%xmm6,%%xmm3			\n\t		mulpd	(%%r8),%%xmm13		\n\t"\
		"movaps	%%xmm0,0x040(%%rsi)		\n\t		addpd	%%xmm10,%%xmm8		\n\t"\
		"movaps	%%xmm2,0x060(%%rsi)		\n\t		addpd	%%xmm14,%%xmm12		\n\t"\
		"movaps	%%xmm1,0x050(%%rsi)		\n\t		addpd	%%xmm11,%%xmm9		\n\t"\
		"movaps	%%xmm3,0x030(%%rsi)		\n\t		addpd	%%xmm15,%%xmm13		\n\t"\
		"mulpd	(%%r8),%%xmm4			\n\t		subpd	%%xmm12,%%xmm8		\n\t"\
		"mulpd	(%%r8),%%xmm7			\n\t		subpd	%%xmm15,%%xmm10		\n\t"\
		"mulpd	(%%r8),%%xmm5			\n\t		subpd	%%xmm13,%%xmm9		\n\t"\
		"mulpd	(%%r8),%%xmm6			\n\t		subpd	%%xmm14,%%xmm11		\n\t"\
		"addpd	%%xmm0,%%xmm4			\n\t		movaps	%%xmm8 ,0x040(%%r14)\n\t"\
		"addpd	%%xmm2,%%xmm7			\n\t		movaps	%%xmm10,0x060(%%r14)\n\t"\
		"addpd	%%xmm1,%%xmm5			\n\t		movaps	%%xmm9 ,0x050(%%r14)\n\t"\
		"addpd	%%xmm3,%%xmm6			\n\t		movaps	%%xmm11,0x030(%%r14)\n\t"\
		"movaps	%%xmm4,     (%%rsi)		\n\t		addpd	%%xmm12,%%xmm12		\n\t"\
		"movaps	%%xmm7,0x020(%%rsi)		\n\t		addpd	%%xmm15,%%xmm15		\n\t"\
		"movaps	%%xmm5,0x010(%%rsi)		\n\t		addpd	%%xmm13,%%xmm13		\n\t"\
		"movaps	%%xmm6,0x070(%%rsi)		\n\t		addpd	%%xmm14,%%xmm14		\n\t"\
		"											addpd	%%xmm8 ,%%xmm12		\n\t"\
		"											addpd	%%xmm10,%%xmm15		\n\t"\
		"											addpd	%%xmm9 ,%%xmm13		\n\t"\
		"											addpd	%%xmm11,%%xmm14		\n\t"\
		"											movaps	%%xmm12,     (%%r14)\n\t"\
		"											movaps	%%xmm15,0x020(%%r14)\n\t"\
		"											movaps	%%xmm13,0x010(%%r14)\n\t"\
		"											movaps	%%xmm14,0x070(%%r14)\n\t"\
	/***************************************************************************************************/\
	/*..and now do four more radix-4 transforms, including the internal and external twiddle factors:  */\
	/***************************************************************************************************/\
	/*...Block 1: 24 MOVapd, 26 ADD/SUBpd, 12 MULpd */	/*...Block 3: 31 MOVapd, 32 ADD/SUBpd, 20 MULpd */\
		"movq	%[__add0],%%rax							\n\t	movslq	%[__p2],%%r10			\n\t"\
		"movslq	%[__p4],%%rbx							\n\t	leaq	-0x40(%%rsi),%%r12	/* r5  */\n\t"\
		"movslq	%[__p8],%%rdi							\n\t	leaq	 0x40(%%rsi),%%r13	/* r13 */\n\t"\
		"shlq	$3,%%rbx								\n\t	movq	%%rbx,%%r11		/* p4 */\n\t"\
		"shlq	$3,%%rdi								\n\t	shlq	$3,%%r10		/* p2 */\n\t"\
		"movq	%[__r1],%%rcx							\n\t	movaps	-0x230(%%r8),%%xmm10	\n\t"/* isrt2 */\
		"movq	%%rsi,%%rdx		/* r9 */				\n\t	movaps	0x100(%%r12),%%xmm12	\n\t"\
		"addq	%%rax,%%r10	/* a[j+p2 ] */				\n\t	movaps	0x110(%%r12),%%xmm13	\n\t"\
		"addq	%%rax,%%rbx	/* a[j+p4 ] */				\n\t	movaps	0x100(%%r13),%%xmm14	\n\t"\
		"addq	%%r10,%%r11	/* a[j+p6 ] */				\n\t	movaps	0x110(%%r13),%%xmm15	\n\t"\
		"movaps	     (%%rdx),%%xmm2						\n\t	mulpd	%%xmm10,%%xmm12			\n\t"\
		"movaps	0x100(%%rdx),%%xmm4						\n\t	mulpd	%%xmm10,%%xmm13			\n\t"\
		"movaps	0x010(%%rdx),%%xmm3						\n\t	mulpd	%%xmm10,%%xmm14			\n\t"\
		"movaps	0x110(%%rdx),%%xmm5						\n\t	mulpd	%%xmm10,%%xmm15			\n\t"\
		"movaps	     (%%rcx),%%xmm0						\n\t	movaps	     (%%r12),%%xmm8		\n\t"\
		"movaps	0x100(%%rcx),%%xmm6						\n\t	movaps	0x010(%%r13),%%xmm10	\n\t"\
		"movaps	0x010(%%rcx),%%xmm1						\n\t	movaps	0x010(%%r12),%%xmm11	\n\t"\
		"movaps	0x110(%%rcx),%%xmm7						\n\t	movaps	     (%%r13),%%xmm9		\n\t"\
		"subpd	%%xmm2,%%xmm0							\n\t	subpd	%%xmm12,%%xmm13			\n\t"\
		"subpd	%%xmm4,%%xmm6							\n\t	subpd	%%xmm15,%%xmm14			\n\t"\
		"subpd	%%xmm3,%%xmm1							\n\t	mulpd	(%%r8),%%xmm12			\n\t"\
		"subpd	%%xmm5,%%xmm7							\n\t	mulpd	(%%r8),%%xmm15			\n\t"\
		"mulpd	(%%r8),%%xmm2							\n\t	addpd	%%xmm13,%%xmm12			\n\t"\
		"mulpd	(%%r8),%%xmm4							\n\t	addpd	%%xmm14,%%xmm15			\n\t"\
		"mulpd	(%%r8),%%xmm3							\n\t	subpd	%%xmm14,%%xmm12			\n\t"\
		"mulpd	(%%r8),%%xmm5							\n\t	subpd	%%xmm10,%%xmm8			\n\t"\
		"addpd	%%xmm0,%%xmm2							\n\t	subpd	%%xmm15,%%xmm13			\n\t"\
		"addpd	%%xmm6,%%xmm4							\n\t	subpd	%%xmm9 ,%%xmm11			\n\t"\
		"addpd	%%xmm1,%%xmm3							\n\t	mulpd	(%%r8),%%xmm14			\n\t"\
		"addpd	%%xmm7,%%xmm5							\n\t	mulpd	(%%r8),%%xmm10			\n\t"\
	/* May 2016: With genuine-BR-ordered roots (to match DIF data layout, have
	SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c0,2,1,3) base-address offsets = r25+0xb0,cc0+0x120,0xa0,0x1a0,
	and base-plus-* byte offsets in the 2 cmul blocks changed from 0x[00,20],[40,60] to 0x[00,40],[20,60]: */\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c0) */	\n\t	mulpd	(%%r8),%%xmm15			\n\t"\
		"addq	$0x1b0,%%rsi	/* c0, from r9 */		\n\t	mulpd	(%%r8),%%xmm9			\n\t"\
		"subpd	%%xmm4,%%xmm2							\n\t	addpd	%%xmm12,%%xmm14			\n\t"\
		"subpd	%%xmm7,%%xmm0							\n\t	addpd	%%xmm8 ,%%xmm10			\n\t"\
		"subpd	%%xmm5,%%xmm3							\n\t	addpd	%%xmm13,%%xmm15			\n\t"\
		"subpd	%%xmm6,%%xmm1							\n\t	addpd	%%xmm11,%%xmm9			\n\t"\
		"mulpd	(%%r8),%%xmm4							\n\t	/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c2) */\n\t"\
		"mulpd	(%%r8),%%xmm7							\n\t	leaq	0x80(%%rsi),%%r14	/* c2, from c0 */\n\t"\
		"mulpd	(%%r8),%%xmm5							\n\t	subpd	%%xmm12,%%xmm10			\n\t"\
		"mulpd	(%%r8),%%xmm6							\n\t	subpd	%%xmm15,%%xmm8			\n\t"\
		"addpd	%%xmm2,%%xmm4							\n\t	subpd	%%xmm13,%%xmm11			\n\t"\
		"addpd	%%xmm0,%%xmm7							\n\t	subpd	%%xmm14,%%xmm9			\n\t"\
		"addpd	%%xmm3,%%xmm5							\n\t	mulpd	(%%r8),%%xmm12			\n\t"\
		"addpd	%%xmm1,%%xmm6							\n\t	mulpd	(%%r8),%%xmm15			\n\t"\
		"movaps	%%xmm2,     (%%rcx)						\n\t	addpd	%%xmm13,%%xmm13			\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)						\n\t	addpd	%%xmm14,%%xmm14			\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)						\n\t	addpd	%%xmm10,%%xmm12			\n\t"\
		"movaps	%%xmm6,     (%%rdx)						\n\t	addpd	%%xmm8 ,%%xmm15			\n\t"\
		"movaps	%%xmm4,%%xmm2							\n\t	addpd	%%xmm11,%%xmm13			\n\t"\
		"movaps	%%xmm7,%%xmm0							\n\t	addpd	%%xmm9 ,%%xmm14			\n\t"\
		"movaps	%%xmm5,%%xmm3							\n\t	movaps	%%xmm10,     (%%r12)	\n\t"\
		"movaps	%%xmm1,%%xmm6							\n\t	movaps	%%xmm8 ,0x010(%%r13)	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2		/* c0,4 */		\n\t	movaps	%%xmm11,0x010(%%r12)	\n\t"\
		"mulpd	    (%%rsi),%%xmm5						\n\t	movaps	%%xmm14,     (%%r13)	\n\t"\
		"mulpd	0x40(%%rsi),%%xmm1						\n\t	movaps	%%xmm12,%%xmm10			\n\t"\
		"mulpd	0x50(%%rsi),%%xmm0						\n\t	movaps	%%xmm15,%%xmm8			\n\t"\
		"mulpd	    (%%rsi),%%xmm4						\n\t	movaps	%%xmm13,%%xmm11			\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3						\n\t	movaps	%%xmm9 ,%%xmm14			\n\t"\
		"mulpd	0x50(%%rsi),%%xmm6						\n\t	mulpd	    (%%r14),%%xmm13		\n\t"/* c2,6 */\
		"mulpd	0x40(%%rsi),%%xmm7						\n\t	mulpd	0x10(%%r14),%%xmm10		\n\t"\
		"subpd	%%xmm2,%%xmm5							\n\t	mulpd	0x40(%%r14),%%xmm9		\n\t"\
		"subpd	%%xmm0,%%xmm1							\n\t	mulpd	0x50(%%r14),%%xmm8		\n\t"\
		"addpd	%%xmm3,%%xmm4							\n\t	mulpd	0x10(%%r14),%%xmm11		\n\t"\
		"addpd	%%xmm6,%%xmm7							\n\t	mulpd	    (%%r14),%%xmm12		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)						\n\t	mulpd	0x40(%%r14),%%xmm15		\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)						\n\t	mulpd	0x50(%%r14),%%xmm14		\n\t"\
		"movaps	%%xmm4,    (%%rax)						\n\t	subpd	%%xmm10,%%xmm13			\n\t"\
		"movaps	%%xmm7,    (%%rbx)						\n\t	subpd	%%xmm8 ,%%xmm9			\n\t"\
		"addq	%%rdi,%%rax	/* a[j+p8 ] */				\n\t	addpd	%%xmm11,%%xmm12			\n\t"\
		"addq	%%rdi,%%rbx	/* a[j+p12] */				\n\t	addpd	%%xmm14,%%xmm15			\n\t"\
		"movaps	     (%%rcx),%%xmm4						\n\t	movaps	%%xmm13,0x10(%%r10)		\n\t"\
		"movaps	0x010(%%rdx),%%xmm0						\n\t	movaps	%%xmm9 ,0x10(%%r11)		\n\t"\
		"movaps	0x010(%%rcx),%%xmm5						\n\t	movaps	%%xmm12,    (%%r10)		\n\t"\
		"movaps	     (%%rdx),%%xmm6						\n\t	movaps	%%xmm15,    (%%r11)		\n\t"\
		"movaps	%%xmm4,%%xmm2							\n\t	addq	%%rdi,%%r10	/* a[j+p10] */	\n\t"\
		"movaps	%%xmm0,%%xmm1							\n\t	addq	%%rdi,%%r11	/* a[j+p14] */	\n\t"\
		"movaps	%%xmm5,%%xmm3							\n\t	movaps	     (%%r12),%%xmm12	\n\t"\
		"movaps	%%xmm6,%%xmm7							\n\t	movaps	0x010(%%r13),%%xmm8		\n\t"\
		"mulpd	0x20(%%rsi),%%xmm5		/* c8,C */		\n\t	movaps	0x010(%%r12),%%xmm13	\n\t"\
		"mulpd	0x30(%%rsi),%%xmm2						\n\t	movaps	     (%%r13),%%xmm14	\n\t"\
		"mulpd	0x70(%%rsi),%%xmm1						\n\t	movaps	%%xmm12,%%xmm10			\n\t"\
		"mulpd	0x60(%%rsi),%%xmm6						\n\t	movaps	%%xmm8 ,%%xmm9			\n\t"\
		"mulpd	0x20(%%rsi),%%xmm4						\n\t	movaps	%%xmm13,%%xmm11			\n\t"\
		"mulpd	0x30(%%rsi),%%xmm3						\n\t	movaps	%%xmm14,%%xmm15			\n\t"\
		"mulpd	0x60(%%rsi),%%xmm0						\n\t	mulpd	0x20(%%r14),%%xmm13		\n\t"/* cA,E */\
		"mulpd	0x70(%%rsi),%%xmm7						\n\t	mulpd	0x30(%%r14),%%xmm10		\n\t"\
		"subpd	%%xmm2,%%xmm5							\n\t	mulpd	0x60(%%r14),%%xmm14		\n\t"\
		"subpd	%%xmm1,%%xmm6							\n\t	mulpd	0x70(%%r14),%%xmm9		\n\t"\
		"addpd	%%xmm3,%%xmm4							\n\t	mulpd	0x20(%%r14),%%xmm12		\n\t"\
		"addpd	%%xmm7,%%xmm0							\n\t	mulpd	0x30(%%r14),%%xmm11		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)						\n\t	mulpd	0x60(%%r14),%%xmm8		\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)						\n\t	mulpd	0x70(%%r14),%%xmm15		\n\t"\
		"movaps	%%xmm4,    (%%rax)						\n\t	subpd	%%xmm10,%%xmm13			\n\t"\
		"movaps	%%xmm0,    (%%rbx)						\n\t	subpd	%%xmm9 ,%%xmm14			\n\t"\
		/*...Block 2: 34 MOVapd, 36 ADD/SUBpd, 26 MULpd */"		addpd	%%xmm11,%%xmm12			\n\t"\
		"movslq	%[__p1],%%r14							\n\t	addpd	%%xmm15,%%xmm8			\n\t"\
		"subq	%%rdi,%%rax	/* a[j+p0 ] */				\n\t	movaps	%%xmm13,0x10(%%r10)		\n\t"\
		"subq	%%rdi,%%rbx	/* a[j+p4 ] */				\n\t	movaps	%%xmm14,0x10(%%r11)		\n\t"\
		"shlq	$3,%%r14								\n\t	movaps	%%xmm12,    (%%r10)		\n\t"\
		"addq	%%r14,%%rax	/* a[j+p1 ] */				\n\t	movaps	%%xmm8 ,    (%%r11)		\n\t"\
		"addq	%%r14,%%rbx	/* a[j+p5 ] */				\n\t"/*...Block 4: 34 MOVapd, 36 ADD/SUBpd, 26 MULpd */\
		"addq	$0x20,%%rcx			/* r3 , from r1 */	\n\t	subq	%%rdi,%%r10	/* a[j+p2 ] */	\n\t"\
		"leaq	0x80(%%rcx),%%rdx	/* r11, from r3 */	\n\t	subq	%%rdi,%%r11	/* a[j+p6 ] */	\n\t"\
		"subq	$0x20,%%rsi		/* cc0, from c0 */		\n\t	addq	%%r14,%%r10	/* a[j+p3 ] */	\n\t"\
		"movaps	0x100(%%rcx),%%xmm4						\n\t	addq	%%r14,%%r11	/* a[j+p7 ] */	\n\t"\
		"movaps	0x110(%%rcx),%%xmm5						\n\t	addq	$0x20,%%r12	/* r7 , from r5  */\n\t"\
		"movaps	     (%%rsi),%%xmm2						\n\t	addq	$0x20,%%r13	/* r15, from r13 */\n\t"\
		"movaps	0x010(%%rsi),%%xmm3						\n\t	addq	$0x100,%%r12			\n\t"\
		"movaps	%%xmm4,%%xmm6							\n\t	addq	$0x100,%%r13			\n\t"\
		"movaps	%%xmm5,%%xmm7							\n\t	movaps	     (%%r12),%%xmm12	\n\t"\
		"mulpd	%%xmm3,%%xmm6							\n\t	movaps	0x010(%%r12),%%xmm13	\n\t"\
		"mulpd	%%xmm2,%%xmm5							\n\t	movaps	     (%%rsi),%%xmm11	\n\t"\
		"mulpd	%%xmm2,%%xmm4							\n\t	movaps	0x010(%%rsi),%%xmm10	\n\t"\
		"mulpd	%%xmm3,%%xmm7							\n\t	movaps	%%xmm12,%%xmm14			\n\t"\
		"movaps	0x100(%%rdx),%%xmm0						\n\t	movaps	%%xmm13,%%xmm15			\n\t"\
		"movaps	0x110(%%rdx),%%xmm1						\n\t	mulpd	%%xmm11,%%xmm14			\n\t"\
		"subpd	%%xmm6,%%xmm5							\n\t	mulpd	%%xmm10,%%xmm13			\n\t"\
		"movaps	%%xmm0,%%xmm6							\n\t	mulpd	%%xmm10,%%xmm12			\n\t"\
		"addpd	%%xmm7,%%xmm4							\n\t	mulpd	%%xmm11,%%xmm15			\n\t"\
		"movaps	%%xmm1,%%xmm7							\n\t	movaps	     (%%r13),%%xmm8		\n\t"\
		"mulpd	%%xmm2,%%xmm6							\n\t	movaps	0x010(%%r13),%%xmm9		\n\t"\
		"mulpd	%%xmm3,%%xmm1							\n\t	subpd	%%xmm14,%%xmm13			\n\t"\
		"mulpd	%%xmm3,%%xmm0							\n\t	movaps	%%xmm8 ,%%xmm14			\n\t"\
		"mulpd	%%xmm2,%%xmm7							\n\t	addpd	%%xmm15,%%xmm12			\n\t"\
		"subpd	%%xmm6,%%xmm1							\n\t	movaps	%%xmm9 ,%%xmm15			\n\t"\
		"addpd	%%xmm7,%%xmm0							\n\t	mulpd	%%xmm10,%%xmm14			\n\t"\
		"movaps	%%xmm5,%%xmm7							\n\t	mulpd	%%xmm11,%%xmm9 			\n\t"\
		"movaps	%%xmm4,%%xmm6							\n\t	mulpd	%%xmm11,%%xmm8 			\n\t"\
		"addpd	%%xmm0,%%xmm4							\n\t	mulpd	%%xmm10,%%xmm15			\n\t"\
		"addpd	%%xmm1,%%xmm5							\n\t	subpd	%%xmm14,%%xmm9			\n\t"\
		"subpd	%%xmm0,%%xmm6							\n\t	addpd	%%xmm15,%%xmm8			\n\t"\
		"subpd	%%xmm1,%%xmm7							\n\t	movaps	%%xmm13,%%xmm15			\n\t"\
		"movaps	     (%%rdx),%%xmm2						\n\t	movaps	%%xmm12,%%xmm14			\n\t"\
		"movaps	0x010(%%rdx),%%xmm3						\n\t	addpd	%%xmm8 ,%%xmm14			\n\t"\
		"movaps	-0x230(%%r8),%%xmm1						\n\t	addpd	%%xmm9 ,%%xmm15			\n\t"\
		"movaps	%%xmm3,%%xmm0							\n\t	subpd	%%xmm8 ,%%xmm12			\n\t"\
		"subpd	%%xmm2,%%xmm3							\n\t	subpd	%%xmm9 ,%%xmm13			\n\t"\
		"addpd	%%xmm0,%%xmm2							\n\t	subq	$0x100,%%r12			\n\t"\
		"mulpd	%%xmm1,%%xmm2							\n\t	subq	$0x100,%%r13			\n\t"\
		"mulpd	%%xmm1,%%xmm3							\n\t	movaps	     (%%r13),%%xmm8		\n\t"\
		"movaps	     (%%rcx),%%xmm0						\n\t	movaps	0x010(%%r13),%%xmm9		\n\t"\
		"movaps	0x010(%%rcx),%%xmm1						\n\t	movaps	-0x230(%%r8),%%xmm11	\n\t"\
		"subpd	%%xmm2,%%xmm0							\n\t	movaps	%%xmm8 ,%%xmm10			\n\t"\
		"subpd	%%xmm3,%%xmm1							\n\t	subpd	%%xmm9 ,%%xmm8			\n\t"\
		"addpd	%%xmm2,%%xmm2							\n\t	addpd	%%xmm10,%%xmm9			\n\t"\
		"addpd	%%xmm3,%%xmm3							\n\t	mulpd	%%xmm11,%%xmm8			\n\t"\
		"addpd	%%xmm0,%%xmm2							\n\t	mulpd	%%xmm11,%%xmm9			\n\t"\
		"addpd	%%xmm1,%%xmm3							\n\t	movaps	     (%%r12),%%xmm10	\n\t"\
		"/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c1) */	\n\t	movaps	0x010(%%r12),%%xmm11	\n\t"\
		"addq	$0x120,%%rsi	/* c1, from cc0 */		\n\t	subpd	%%xmm8 ,%%xmm10			\n\t"\
		"subpd	%%xmm4,%%xmm2							\n\t	subpd	%%xmm9 ,%%xmm11			\n\t"\
		"subpd	%%xmm7,%%xmm0							\n\t	addpd	%%xmm8 ,%%xmm8			\n\t"\
		"subpd	%%xmm5,%%xmm3							\n\t	addpd	%%xmm9 ,%%xmm9			\n\t"\
		"subpd	%%xmm6,%%xmm1							\n\t	addpd	%%xmm10,%%xmm8			\n\t"\
		"mulpd	(%%r8),%%xmm4							\n\t	addpd	%%xmm11,%%xmm9			\n\t"\
		"mulpd	(%%r8),%%xmm7							\n\t/* SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c3) */\n\t"\
		"mulpd	(%%r8),%%xmm5							\n\t	leaq	0x80(%%rsi),%%r14	/* c3, from c1 */\n\t"\
		"mulpd	(%%r8),%%xmm6							\n\t	subpd	%%xmm12,%%xmm10			\n\t"\
		"addpd	%%xmm2,%%xmm4							\n\t	subpd	%%xmm15,%%xmm8			\n\t"\
		"addpd	%%xmm0,%%xmm7							\n\t	subpd	%%xmm13,%%xmm11			\n\t"\
		"addpd	%%xmm3,%%xmm5							\n\t	subpd	%%xmm14,%%xmm9			\n\t"\
		"addpd	%%xmm1,%%xmm6							\n\t	addpd	%%xmm12,%%xmm12			\n\t"\
		"movaps	%%xmm2,     (%%rcx)						\n\t	addpd	%%xmm15,%%xmm15			\n\t"\
		"movaps	%%xmm0,0x010(%%rdx)						\n\t	addpd	%%xmm13,%%xmm13			\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)						\n\t	addpd	%%xmm14,%%xmm14			\n\t"\
		"movaps	%%xmm6,     (%%rdx)						\n\t	addpd	%%xmm10,%%xmm12			\n\t"\
		"movaps	%%xmm4,%%xmm2							\n\t	addpd	%%xmm8 ,%%xmm15			\n\t"\
		"movaps	%%xmm7,%%xmm0							\n\t	addpd	%%xmm11,%%xmm13			\n\t"\
		"movaps	%%xmm5,%%xmm3							\n\t	addpd	%%xmm9 ,%%xmm14			\n\t"\
		"movaps	%%xmm1,%%xmm6							\n\t	movaps	%%xmm10,     (%%r12)	\n\t"\
		"mulpd	    (%%rsi),%%xmm5		/* c1,5 */		\n\t	movaps	%%xmm8 ,0x010(%%r13)	\n\t"\
		"mulpd	0x10(%%rsi),%%xmm2						\n\t	movaps	%%xmm11,0x010(%%r12)	\n\t"\
		"mulpd	0x40(%%rsi),%%xmm1						\n\t	movaps	%%xmm14,     (%%r13)	\n\t"\
		"mulpd	0x50(%%rsi),%%xmm0						\n\t	movaps	%%xmm12,%%xmm10			\n\t"\
		"mulpd	    (%%rsi),%%xmm4						\n\t	movaps	%%xmm15,%%xmm8			\n\t"\
		"mulpd	0x10(%%rsi),%%xmm3						\n\t	movaps	%%xmm13,%%xmm11			\n\t"\
		"mulpd	0x50(%%rsi),%%xmm6						\n\t	movaps	%%xmm9 ,%%xmm14			\n\t"\
		"mulpd	0x40(%%rsi),%%xmm7						\n\t	mulpd	    (%%r14),%%xmm13		\n\t"/* c3,7 */\
		"subpd	%%xmm2,%%xmm5							\n\t	mulpd	0x10(%%r14),%%xmm10		\n\t"\
		"subpd	%%xmm0,%%xmm1							\n\t	mulpd	0x40(%%r14),%%xmm9		\n\t"\
		"addpd	%%xmm3,%%xmm4							\n\t	mulpd	0x50(%%r14),%%xmm8		\n\t"\
		"addpd	%%xmm6,%%xmm7							\n\t	mulpd	    (%%r14),%%xmm12		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)						\n\t	mulpd	0x10(%%r14),%%xmm11		\n\t"\
		"movaps	%%xmm1,0x10(%%rbx)						\n\t	mulpd	0x40(%%r14),%%xmm15		\n\t"\
		"movaps	%%xmm4,    (%%rax)						\n\t	mulpd	0x50(%%r14),%%xmm14		\n\t"\
		"movaps	%%xmm7,    (%%rbx)						\n\t	subpd	%%xmm10,%%xmm13			\n\t"\
		"addq	%%rdi,%%rax								\n\t	subpd	%%xmm8 ,%%xmm9			\n\t"\
		"addq	%%rdi,%%rbx								\n\t	addpd	%%xmm11,%%xmm12			\n\t"\
		"movaps	     (%%rcx),%%xmm4						\n\t	addpd	%%xmm14,%%xmm15			\n\t"\
		"movaps	0x010(%%rdx),%%xmm0						\n\t	movaps	%%xmm13,0x10(%%r10)		\n\t"\
		"movaps	0x010(%%rcx),%%xmm5						\n\t	movaps	%%xmm9 ,0x10(%%r11)		\n\t"\
		"movaps	     (%%rdx),%%xmm6						\n\t	movaps	%%xmm12,    (%%r10)		\n\t"\
		"movaps	%%xmm4,%%xmm2							\n\t	movaps	%%xmm15,    (%%r11)		\n\t"\
		"movaps	%%xmm0,%%xmm1							\n\t	addq	%%rdi,%%r10				\n\t"\
		"movaps	%%xmm5,%%xmm3							\n\t	addq	%%rdi,%%r11				\n\t"\
		"movaps	%%xmm6,%%xmm7							\n\t	movaps	     (%%r12),%%xmm12	\n\t"\
		"mulpd	0x20(%%rsi),%%xmm5		/* c9,D */		\n\t	movaps	0x010(%%r13),%%xmm8		\n\t"\
		"mulpd	0x30(%%rsi),%%xmm2						\n\t	movaps	0x010(%%r12),%%xmm13	\n\t"\
		"mulpd	0x60(%%rsi),%%xmm6						\n\t	movaps	     (%%r13),%%xmm14	\n\t"\
		"mulpd	0x70(%%rsi),%%xmm1						\n\t	movaps	%%xmm12,%%xmm10			\n\t"\
		"mulpd	0x20(%%rsi),%%xmm4						\n\t	movaps	%%xmm8 ,%%xmm9			\n\t"\
		"mulpd	0x30(%%rsi),%%xmm3						\n\t	movaps	%%xmm13,%%xmm11			\n\t"\
		"mulpd	0x60(%%rsi),%%xmm0						\n\t	movaps	%%xmm14,%%xmm15			\n\t"\
		"mulpd	0x70(%%rsi),%%xmm7						\n\t	mulpd	0x20(%%r14),%%xmm13		\n\t"/* cB,F */\
		"subpd	%%xmm2,%%xmm5							\n\t	mulpd	0x30(%%r14),%%xmm10		\n\t"\
		"subpd	%%xmm1,%%xmm6							\n\t	mulpd	0x60(%%r14),%%xmm14		\n\t"\
		"addpd	%%xmm3,%%xmm4							\n\t	mulpd	0x70(%%r14),%%xmm9		\n\t"\
		"addpd	%%xmm7,%%xmm0							\n\t	mulpd	0x20(%%r14),%%xmm12		\n\t"\
		"movaps	%%xmm5,0x10(%%rax)						\n\t	mulpd	0x30(%%r14),%%xmm11		\n\t"\
		"movaps	%%xmm6,0x10(%%rbx)						\n\t	mulpd	0x60(%%r14),%%xmm8		\n\t"\
		"movaps	%%xmm4,    (%%rax)						\n\t	mulpd	0x70(%%r14),%%xmm15		\n\t"\
		"movaps	%%xmm0,    (%%rbx)						\n\t	subpd	%%xmm10,%%xmm13			\n\t"\
		"												\n\t	subpd	%%xmm9 ,%%xmm14			\n\t"\
		"												\n\t	addpd	%%xmm11,%%xmm12			\n\t"\
		"												\n\t	addpd	%%xmm15,%%xmm8			\n\t"\
		"												\n\t	movaps	%%xmm13,0x10(%%r10)		\n\t"\
		"												\n\t	movaps	%%xmm14,0x10(%%r11)		\n\t"\
		"												\n\t	movaps	%%xmm12,    (%%r10)		\n\t"\
		"												\n\t	movaps	%%xmm8 ,    (%%r11)		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__p1] "m" (Xp1)\
		 ,[__p2] "m" (Xp2)\
		 ,[__p3] "m" (Xp3)\
		 ,[__p4] "m" (Xp4)\
		 ,[__p8] "m" (Xp8)\
		 ,[__r1] "m" (Xr1)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__pfetch_addr1] "m" (Xpfetch_addr1)\
		 ,[__pfetch_dist]  "e" (Xpfetch_dist)\
		: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r8","r10","r11","r12","r13","r14","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
	);\
	}

#elif defined(USE_SSE2) && (OS_BITS == 32)

	#error 32-bit OSes no longer supported for SIMD builds!

#else

	#error Unhandled combination of preprocessor flags!

#endif	// x86 simd version ?

#endif	/* radix16_dif_dit_pass_h_included */

