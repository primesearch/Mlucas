/*******************************************************************************
*                                                                              *
*   (C) 1997-2014 by Ernst W. Mayer.                                           *
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
#ifndef radix15_sse_macro_h_included
#define radix15_sse_macro_h_included

	#include "sse2_macro.h"

	// *NOTE:* Currently, GCC_ASM_FULL_INLINE only toggles inlining of the 15 radix-4 subDFTs done as part of each radix-60 DFT:
	#if defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
	  #ifdef USE_AVX
		#define GCC_ASM_FULL_INLINE  0	// 0 to use small-macros to assemble radix-60 DFTs, 1 to inline fuse macros as a few big blobs of asm (64-bit only)
		#define USE_64BIT_ASM_STYLE  1
	  #else
		#define GCC_ASM_FULL_INLINE  0	// 0 to use small-macros to assemble radix-60 DFTs, 1 to inline fuse macros as a few big blobs of asm (64-bit only)
		#define USE_64BIT_ASM_STYLE  1
	  #endif
	#else
		#undef GCC_ASM_FULL_INLINE
	#endif

	#ifndef USE_64BIT_ASM_STYLE	// Default is to use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just xmm0-7.

		/* General indexing for twiddleless radix-15 done as 3*radix-5 followed by 5*radix-3 is as for the scalar macro above:
		RADIX_15_DIF(00,01,02,03,04,05,06,07,08,09,0A,0B,0C,0D,0E)
		->
			RADIX_05_DFT(i0,iC,i9,i6,i3, t0,t1,t2,t3,t4)
			RADIX_05_DFT(iA,i7,i4,i1,iD, t5,t6,t7,t8,t9)
			RADIX_05_DFT(i5,i2,iE,iB,i8, tA,tB,tC,tD,tE)

			RADIX_03_DFT(t0,t5,tA, o0,o1,o2,)
			RADIX_03_DFT(t1,t6,tB, oD,oE,oB,)
			RADIX_03_DFT(t2,t7,tC, o9,oA,oB,)
			RADIX_03_DFT(t3,t8,tD, o8,o6,o7,)
			RADIX_03_DFT(t4,t9,tE, o4,o5,o3,)

		In our impl below, the __i are input pointers, which may overlap the __o outputs;
		..cc0 and cc1 are ptrs to the radix-3 and radix-5 SSE2 sincos constants (c3m1 and cn1);
		__t0-E are ptr to scratch local storage (i.e. the address block pointed to by r00-r3e).
		*/
	// Aug 2014: Need arbitrary-pointer-offsets to support I/O permutations needed by
	// larger-radix DFTs of length 15 * 2^n
	//	#define USE_LITERAL_BYTE_OFFSETS	// undef this to revert to the 1st set of macros below:

	  #ifndef USE_LITERAL_BYTE_OFFSETS

		#define SSE2_RADIX_15_DIF(\
			__cc0, __cc1,\
			__i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA,__iB,__iC,__iD,__iE,\
			__t0,__t1,__t2,__t3,__t4,__t5,__t6,__t7,__t8,__t9,__tA,__tB,__tC,__tD,__tE,\
			__o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA,__oB,__oC,__oD,__oE)\
		{\
			SSE2_RADIX_05_DFT_0TWIDDLE(__i0,__iC,__i9,__i6,__i3, __cc1, __t0,__t1,__t2,__t3,__t4);\
			SSE2_RADIX_05_DFT_0TWIDDLE(__iA,__i7,__i4,__i1,__iD, __cc1, __t5,__t6,__t7,__t8,__t9);\
			SSE2_RADIX_05_DFT_0TWIDDLE(__i5,__i2,__iE,__iB,__i8, __cc1, __tA,__tB,__tC,__tD,__tE);\
		\
			SSE2_RADIX_03_DFT(__t0,__t5,__tA, __cc0, __o0,__o1,__o2);\
			SSE2_RADIX_03_DFT(__t1,__t6,__tB, __cc0, __oD,__oE,__oC);\
			SSE2_RADIX_03_DFT(__t2,__t7,__tC, __cc0, __o9,__oA,__oB);\
			SSE2_RADIX_03_DFT(__t3,__t8,__tD, __cc0, __o8,__o6,__o7);\
			SSE2_RADIX_03_DFT(__t4,__t9,__tE, __cc0, __o4,__o5,__o3);\
		/*\
			fprintf(stderr, "radix15 dif: out[0] = %24.5f, %24.5f\n",__o0->re,(__o0+1)->re);\
			fprintf(stderr, "radix15 dif: out[1] = %24.5f, %24.5f\n",__o1->re,(__o1+1)->re);\
			fprintf(stderr, "radix15 dif: out[2] = %24.5f, %24.5f\n",__o2->re,(__o2+1)->re);\
			fprintf(stderr, "radix15 dif: out[3] = %24.5f, %24.5f\n",__o3->re,(__o3+1)->re);\
			fprintf(stderr, "radix15 dif: out[4] = %24.5f, %24.5f\n",__o4->re,(__o4+1)->re);\
			fprintf(stderr, "radix15 dif: out[5] = %24.5f, %24.5f\n",__o5->re,(__o5+1)->re);\
			fprintf(stderr, "radix15 dif: out[6] = %24.5f, %24.5f\n",__o6->re,(__o6+1)->re);\
			fprintf(stderr, "radix15 dif: out[7] = %24.5f, %24.5f\n",__o7->re,(__o7+1)->re);\
			fprintf(stderr, "radix15 dif: out[8] = %24.5f, %24.5f\n",__o8->re,(__o8+1)->re);\
			fprintf(stderr, "radix15 dif: out[9] = %24.5f, %24.5f\n",__o9->re,(__o9+1)->re);\
			fprintf(stderr, "radix15 dif: out[A] = %24.5f, %24.5f\n",__oA->re,(__oA+1)->re);\
			fprintf(stderr, "radix15 dif: out[B] = %24.5f, %24.5f\n",__oB->re,(__oB+1)->re);\
			fprintf(stderr, "radix15 dif: out[C] = %24.5f, %24.5f\n",__oC->re,(__oC+1)->re);\
			fprintf(stderr, "radix15 dif: out[D] = %24.5f, %24.5f\n",__oD->re,(__oD+1)->re);\
			fprintf(stderr, "radix15 dif: out[E] = %24.5f, %24.5f\n\n",__oE->re,(__oE+1)->re);\
		*/\
		}

		#define SSE2_RADIX_15_DIT(\
			__cc0, __cc1,\
			__i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA,__iB,__iC,__iD,__iE,\
			__t0,__t1,__t2,__t3,__t4,__t5,__t6,__t7,__t8,__t9,__tA,__tB,__tC,__tD,__tE,\
			__o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA,__oB,__oC,__oD,__oE)\
		{\
		/* Swap the 2nd pair of each output triplet to effect iDFT: */\
			SSE2_RADIX_03_DFT(__i0,__i2,__i1, __cc0, __t0,__t2,__t1);\
			SSE2_RADIX_03_DFT(__i8,__i7,__i6, __cc0, __t3,__t5,__t4);\
			SSE2_RADIX_03_DFT(__iD,__iC,__iE, __cc0, __t6,__t8,__t7);\
			SSE2_RADIX_03_DFT(__i4,__i3,__i5, __cc0, __t9,__tB,__tA);\
			SSE2_RADIX_03_DFT(__i9,__iB,__iA, __cc0, __tC,__tE,__tD);\
		\
		/* Output perm here is 0123456789abcde --> 05a6b1c2738d9e4: */\
			SSE2_RADIX_05_DFT_0TWIDDLE(__t0,__t3,__t6,__t9,__tC, __cc1, __o0,__o6,__oC,__o3,__o9);\
			SSE2_RADIX_05_DFT_0TWIDDLE(__t1,__t4,__t7,__tA,__tD, __cc1, __o5,__oB,__o2,__o8,__oE);\
			SSE2_RADIX_05_DFT_0TWIDDLE(__t2,__t5,__t8,__tB,__tE, __cc1, __oA,__o1,__o7,__oD,__o4);\
		/*\
			fprintf(stderr, "radix15 dit: out[0] = %24.5f, %24.5f\n",__o0->re,(__o0+1)->re);\
			fprintf(stderr, "radix15 dit: out[1] = %24.5f, %24.5f\n",__o1->re,(__o1+1)->re);\
			fprintf(stderr, "radix15 dit: out[2] = %24.5f, %24.5f\n",__o2->re,(__o2+1)->re);\
			fprintf(stderr, "radix15 dit: out[3] = %24.5f, %24.5f\n",__o3->re,(__o3+1)->re);\
			fprintf(stderr, "radix15 dit: out[4] = %24.5f, %24.5f\n",__o4->re,(__o4+1)->re);\
			fprintf(stderr, "radix15 dit: out[5] = %24.5f, %24.5f\n",__o5->re,(__o5+1)->re);\
			fprintf(stderr, "radix15 dit: out[6] = %24.5f, %24.5f\n",__o6->re,(__o6+1)->re);\
			fprintf(stderr, "radix15 dit: out[7] = %24.5f, %24.5f\n",__o7->re,(__o7+1)->re);\
			fprintf(stderr, "radix15 dit: out[8] = %24.5f, %24.5f\n",__o8->re,(__o8+1)->re);\
			fprintf(stderr, "radix15 dit: out[9] = %24.5f, %24.5f\n",__o9->re,(__o9+1)->re);\
			fprintf(stderr, "radix15 dit: out[A] = %24.5f, %24.5f\n",__oA->re,(__oA+1)->re);\
			fprintf(stderr, "radix15 dit: out[B] = %24.5f, %24.5f\n",__oB->re,(__oB+1)->re);\
			fprintf(stderr, "radix15 dit: out[C] = %24.5f, %24.5f\n",__oC->re,(__oC+1)->re);\
			fprintf(stderr, "radix15 dit: out[D] = %24.5f, %24.5f\n",__oD->re,(__oD+1)->re);\
			fprintf(stderr, "radix15 dit: out[E] = %24.5f, %24.5f\n\n",__oE->re,(__oE+1)->re);\
		*/\
		}

	  #else	// USE_LITERAL_BYTE_OFFSETS = True: Versions using base-address-plus-literal-byte-offsets:

	   #if (OS_BITS == 32)

		#define SSE2_RADIX_03_DFT_X1(Xcc, XI,Xi0,Xi1,Xi2, XO,Xo0,Xo1,Xo2)\
		{\
		__asm__ volatile (\
		"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
			"movl	%[__I],%%eax			\n\t"\
			"leal	0x10(%%eax),%%ebx		\n\t"\
			"movaps	%c[__i1](%%eax),%%xmm4	\n\t"\
			"movaps	%c[__i1](%%ebx),%%xmm5	\n\t"\
			"subpd	%c[__i2](%%eax),%%xmm4	\n\t"\
			"subpd	%c[__i2](%%ebx),%%xmm5	\n\t"\
			"movaps	%c[__i2](%%eax),%%xmm2	\n\t"\
			"movaps	%c[__i2](%%ebx),%%xmm3	\n\t"\
			"addpd	%%xmm2,%%xmm2			\n\t"\
			"addpd	%%xmm3,%%xmm3			\n\t"\
			"movl	%[__cc],%%esi			\n\t"\
			"addpd	%%xmm4,%%xmm2			\n\t"\
			"addpd	%%xmm5,%%xmm3			\n\t"\
			"movaps	%c[__i0](%%eax),%%xmm0	\n\t"\
			"movaps	%c[__i0](%%ebx),%%xmm1	\n\t"\
			"movl	%[__O],%%ecx			\n\t"\
			"leal	0x10(%%ecx),%%edx		\n\t"\
			"addpd	%%xmm2,%%xmm0			\n\t"\
			"addpd	%%xmm3,%%xmm1			\n\t"\
			"movaps	    (%%esi),%%xmm6		\n\t"\
			"movaps	0x10(%%esi),%%xmm7		\n\t"\
			"movaps	%%xmm0,%c[__o0](%%ecx)	\n\t"\
			"movaps	%%xmm1,%c[__o0](%%edx)	\n\t"\
			"mulpd	%%xmm6,%%xmm2			\n\t"\
			"mulpd	%%xmm6,%%xmm3			\n\t"\
			"mulpd	%%xmm7,%%xmm4			\n\t"\
			"mulpd	%%xmm7,%%xmm5			\n\t"\
			"addpd	%%xmm0,%%xmm2			\n\t"\
			"addpd	%%xmm1,%%xmm3			\n\t"\
			"movaps	%%xmm2,%%xmm0			\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"subpd	%%xmm5,%%xmm2			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"addpd	%%xmm5,%%xmm0			\n\t"\
			"subpd	%%xmm4,%%xmm1			\n\t"\
			"movaps	%%xmm2,%c[__o1](%%ecx)	\n\t"\
			"movaps	%%xmm3,%c[__o1](%%edx)	\n\t"\
			"movaps	%%xmm0,%c[__o2](%%ecx)	\n\t"\
			"movaps	%%xmm1,%c[__o2](%%edx)	\n\t"\
		"popl %%ebx	\n\t"\
			:					/* outputs: none */\
			: [__cc] "m" (Xcc)	/* m-inputs from memory addresses, e-inputs literal byte offsets here */\
			 ,[__I] "m" (XI)\
			 ,[__i0] "e" (Xi0)\
			 ,[__i1] "e" (Xi1)\
			 ,[__i2] "e" (Xi2)\
			 ,[__O] "m" (XO)\
			 ,[__o0] "e" (Xo0)\
			 ,[__o1] "e" (Xo1)\
			 ,[__o2] "e" (Xo2)\
			: "cc","memory","eax",/*"ebx",*/"ecx","edx","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
		}

		#define SSE2_RADIX_05_DFT_0TWIDDLE_X1(Xcc, XI,Xi0,Xi1,Xi2,Xi3,Xi4, XO,Xo0,Xo1,Xo2,Xo3,Xo4)\
		{\
		__asm__ volatile (\
		"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
			"movl	%[__I],%%eax			\n\t"\
			"leal	0x10(%%eax),%%ebx		\n\t"\
			"movaps	%c[__i1](%%eax),%%xmm0	\n\t"\
			"movaps	%c[__i1](%%ebx),%%xmm1	\n\t"\
			"movaps	%c[__i4](%%eax),%%xmm6	\n\t"\
			"movaps	%c[__i4](%%ebx),%%xmm7	\n\t"\
			"subpd	%%xmm6,%%xmm0			\n\t"\
			"subpd	%%xmm7,%%xmm1			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm0,%%xmm6			\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
		/* Diff-indent the spill/fill stuff to highlight the diffs from the xmm8/9-using 64-bit version */\
		"movaps	%%xmm0,%c[__i1](%%eax)	\n\t"/* spill A */\
		"movaps	%%xmm1,%c[__i1](%%ebx)	\n\t"\
			"movaps	%c[__i2](%%eax),%%xmm2	\n\t"\
			"movaps	%c[__i2](%%ebx),%%xmm3	\n\t"\
			"movaps	%c[__i3](%%eax),%%xmm4	\n\t"\
			"movaps	%c[__i3](%%ebx),%%xmm5	\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"subpd	%%xmm5,%%xmm3			\n\t"\
			"addpd	%%xmm4,%%xmm4			\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm2,%%xmm4			\n\t"\
			"addpd	%%xmm3,%%xmm5			\n\t"\
			"movl	%[__O],%%ecx			\n\t"\
			"leal	0x10(%%ecx),%%edx		\n\t"\
			"subpd	%%xmm4,%%xmm6			\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
		"movaps	%c[__i0](%%eax),%%xmm0	\n\t"\
		"movaps	%c[__i0](%%ebx),%%xmm1	\n\t"\
			"addpd	%%xmm4,%%xmm4			\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm4			\n\t"\
			"addpd	%%xmm7,%%xmm5			\n\t"\
			"movl	%[__cc],%%esi			\n\t"\
		"addpd	%%xmm0,%%xmm4	\n\t"\
		"addpd	%%xmm1,%%xmm5	\n\t"\
			"movaps	%%xmm4,%c[__o0](%%ecx)	\n\t"\
			"movaps	%%xmm5,%c[__o0](%%edx)	\n\t"\
			"mulpd	0x10(%%esi),%%xmm6		\n\t"\
			"mulpd	0x10(%%esi),%%xmm7		\n\t"\
			"subpd	%%xmm4,%%xmm0	\n\t"\
			"subpd	%%xmm5,%%xmm1	\n\t"\
			"mulpd	    (%%esi),%%xmm0		\n\t"\
			"mulpd	    (%%esi),%%xmm1		\n\t"\
			"subpd	%%xmm0,%%xmm4	\n\t"\
			"subpd	%%xmm1,%%xmm5	\n\t"\
			"subpd	%%xmm6,%%xmm4			\n\t"\
			"subpd	%%xmm7,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm4,%%xmm6			\n\t"\
			"addpd	%%xmm5,%%xmm7			\n\t"\
		"movaps	%%xmm4,%c[__i0](%%eax)	\n\t"/* spill B */\
		"movaps	%%xmm5,%c[__i0](%%ebx)	\n\t"\
		"movaps	%c[__i1](%%eax),%%xmm0	\n\t"/* restore spill A */\
		"movaps	%c[__i1](%%ebx),%%xmm1	\n\t"\
		"movaps	%%xmm0,%%xmm4			\n\t"\
		"movaps	%%xmm1,%%xmm5			\n\t"\
			"subpd	%%xmm2,%%xmm0			\n\t"\
			"subpd	%%xmm3,%%xmm1			\n\t"\
			"mulpd	0x20(%%esi),%%xmm0		\n\t"\
			"mulpd	0x20(%%esi),%%xmm1		\n\t"\
			"mulpd	0x30(%%esi),%%xmm2		\n\t"\
			"mulpd	0x30(%%esi),%%xmm3		\n\t"\
			"mulpd	0x40(%%esi),%%xmm4		\n\t"\
			"mulpd	0x40(%%esi),%%xmm5		\n\t"\
			"addpd	%%xmm0,%%xmm2			\n\t"\
			"addpd	%%xmm1,%%xmm3			\n\t"\
			"subpd	%%xmm4,%%xmm0			\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t"\
		"movaps	%c[__i0](%%eax),%%xmm4	\n\t"/* restore spill B */\
		"movaps	%c[__i0](%%ebx),%%xmm5	\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"subpd	%%xmm2,%%xmm7			\n\t"\
			"addpd	%%xmm3,%%xmm3			\n\t"\
			"addpd	%%xmm2,%%xmm2			\n\t"\
			"movaps	%%xmm6,%c[__o1](%%ecx)	\n\t"\
			"movaps	%%xmm7,%c[__o4](%%edx)	\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t"\
			"addpd	%%xmm7,%%xmm2			\n\t"\
			"movaps	%%xmm3,%c[__o4](%%ecx)	\n\t"\
			"movaps	%%xmm2,%c[__o1](%%edx)	\n\t"\
			"subpd	%%xmm1,%%xmm4			\n\t"\
			"subpd	%%xmm0,%%xmm5			\n\t"\
			"addpd	%%xmm1,%%xmm1			\n\t"\
			"addpd	%%xmm0,%%xmm0			\n\t"\
			"movaps	%%xmm4,%c[__o2](%%ecx)	\n\t"\
			"movaps	%%xmm5,%c[__o3](%%edx)	\n\t"\
			"addpd	%%xmm4,%%xmm1			\n\t"\
			"addpd	%%xmm5,%%xmm0			\n\t"\
			"movaps	%%xmm1,%c[__o3](%%ecx)	\n\t"\
			"movaps	%%xmm0,%c[__o2](%%edx)	\n\t"\
		"popl %%ebx	\n\t"\
			:					/* outputs: none */\
			: [__cc] "m" (Xcc)	/* m-inputs from memory addresses, e-inputs literal byte offsets here */\
			 ,[__I] "m" (XI)\
			 ,[__i0] "e" (Xi0)\
			 ,[__i1] "e" (Xi1)\
			 ,[__i2] "e" (Xi2)\
			 ,[__i3] "e" (Xi3)\
			 ,[__i4] "e" (Xi4)\
			 ,[__O] "m" (XO)\
			 ,[__o0] "e" (Xo0)\
			 ,[__o1] "e" (Xo1)\
			 ,[__o2] "e" (Xo2)\
			 ,[__o3] "e" (Xo3)\
			 ,[__o4] "e" (Xo4)\
			: "cc","memory","eax",/*"ebx",*/"ecx","edx","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
		}

	   #else

		#define SSE2_RADIX_03_DFT_X1(Xcc, XI,Xi0,Xi1,Xi2, XO,Xo0,Xo1,Xo2)\
		{\
		__asm__ volatile (\
			"movq	%[__I],%%rax			\n\t"\
			"leaq	0x10(%%rax),%%rbx		\n\t"\
			"movaps	%c[__i1](%%rax),%%xmm4	\n\t"\
			"movaps	%c[__i1](%%rbx),%%xmm5	\n\t"\
			"subpd	%c[__i2](%%rax),%%xmm4	\n\t"\
			"subpd	%c[__i2](%%rbx),%%xmm5	\n\t"\
			"movaps	%c[__i2](%%rax),%%xmm2	\n\t"\
			"movaps	%c[__i2](%%rbx),%%xmm3	\n\t"\
			"addpd	%%xmm2,%%xmm2			\n\t"\
			"addpd	%%xmm3,%%xmm3			\n\t"\
			"movq	%[__cc],%%rsi			\n\t"\
			"addpd	%%xmm4,%%xmm2			\n\t"\
			"addpd	%%xmm5,%%xmm3			\n\t"\
			"movaps	%c[__i0](%%rax),%%xmm0	\n\t"\
			"movaps	%c[__i0](%%rbx),%%xmm1	\n\t"\
			"movq	%[__O],%%rcx			\n\t"\
			"leaq	0x10(%%rcx),%%rdx		\n\t"\
			"addpd	%%xmm2,%%xmm0			\n\t"\
			"addpd	%%xmm3,%%xmm1			\n\t"\
			"movaps	    (%%rsi),%%xmm6		\n\t"\
			"movaps	0x10(%%rsi),%%xmm7		\n\t"\
			"movaps	%%xmm0,%c[__o0](%%rcx)	\n\t"\
			"movaps	%%xmm1,%c[__o0](%%rdx)	\n\t"\
			"mulpd	%%xmm6,%%xmm2			\n\t"\
			"mulpd	%%xmm6,%%xmm3			\n\t"\
			"mulpd	%%xmm7,%%xmm4			\n\t"\
			"mulpd	%%xmm7,%%xmm5			\n\t"\
			"addpd	%%xmm0,%%xmm2			\n\t"\
			"addpd	%%xmm1,%%xmm3			\n\t"\
			"movaps	%%xmm2,%%xmm0			\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"subpd	%%xmm5,%%xmm2			\n\t"\
			"addpd	%%xmm4,%%xmm3			\n\t"\
			"addpd	%%xmm5,%%xmm0			\n\t"\
			"subpd	%%xmm4,%%xmm1			\n\t"\
			"movaps	%%xmm2,%c[__o1](%%rcx)	\n\t"\
			"movaps	%%xmm3,%c[__o1](%%rdx)	\n\t"\
			"movaps	%%xmm0,%c[__o2](%%rcx)	\n\t"\
			"movaps	%%xmm1,%c[__o2](%%rdx)	\n\t"\
			:					/* outputs: none */\
			: [__cc] "m" (Xcc)	/* m-inputs from memory addresses, e-inputs literal byte offsets here */\
			 ,[__I] "m" (XI)\
			 ,[__i0] "e" (Xi0)\
			 ,[__i1] "e" (Xi1)\
			 ,[__i2] "e" (Xi2)\
			 ,[__O] "m" (XO)\
			 ,[__o0] "e" (Xo0)\
			 ,[__o1] "e" (Xo1)\
			 ,[__o2] "e" (Xo2)\
			: "cc","memory","rax","rbx","rcx","rdx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
		}

		#define SSE2_RADIX_05_DFT_0TWIDDLE_X1(Xcc, XI,Xi0,Xi1,Xi2,Xi3,Xi4, XO,Xo0,Xo1,Xo2,Xo3,Xo4)\
		{\
		__asm__ volatile (\
			"movq	%[__I],%%rax			\n\t"\
			"leaq	0x10(%%rax),%%rbx		\n\t"\
			"movaps	%c[__i1](%%rax),%%xmm0	\n\t"\
			"movaps	%c[__i1](%%rbx),%%xmm1	\n\t"\
			"movaps	%c[__i4](%%rax),%%xmm6	\n\t"\
			"movaps	%c[__i4](%%rbx),%%xmm7	\n\t"\
			"subpd	%%xmm6,%%xmm0			\n\t"\
			"subpd	%%xmm7,%%xmm1			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm0,%%xmm6			\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"movaps	%c[__i2](%%rax),%%xmm2	\n\t"\
			"movaps	%c[__i2](%%rbx),%%xmm3	\n\t"\
			"movaps	%c[__i3](%%rax),%%xmm4	\n\t"\
			"movaps	%c[__i3](%%rbx),%%xmm5	\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"subpd	%%xmm5,%%xmm3			\n\t"\
			"addpd	%%xmm4,%%xmm4			\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm2,%%xmm4			\n\t"\
			"addpd	%%xmm3,%%xmm5			\n\t"\
			"movq	%[__O],%%rcx			\n\t"\
			"leaq	0x10(%%rcx),%%rdx		\n\t"\
			"subpd	%%xmm4,%%xmm6			\n\t"\
			"subpd	%%xmm5,%%xmm7			\n\t"\
			"movaps	%c[__i0](%%rax),%%xmm8	\n\t"\
			"movaps	%c[__i0](%%rbx),%%xmm9	\n\t"\
			"addpd	%%xmm4,%%xmm4			\n\t"\
			"addpd	%%xmm5,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm4			\n\t"\
			"addpd	%%xmm7,%%xmm5			\n\t"\
			"movq	%[__cc],%%rsi			\n\t"\
			"addpd	%%xmm8,%%xmm4	\n\t"\
			"addpd	%%xmm9,%%xmm5	\n\t"\
			"movaps	%%xmm4,%c[__o0](%%rcx)	\n\t"\
			"movaps	%%xmm5,%c[__o0](%%rdx)	\n\t"\
			"mulpd	0x10(%%rsi),%%xmm6		\n\t"\
			"mulpd	0x10(%%rsi),%%xmm7		\n\t"\
			"subpd	%%xmm4,%%xmm8	\n\t"\
			"subpd	%%xmm5,%%xmm9	\n\t"\
			"mulpd	    (%%rsi),%%xmm8		\n\t"\
			"mulpd	    (%%rsi),%%xmm9		\n\t"\
			"subpd	%%xmm8,%%xmm4	\n\t"\
			"subpd	%%xmm9,%%xmm5	\n\t"\
			"subpd	%%xmm6,%%xmm4			\n\t"\
			"subpd	%%xmm7,%%xmm5			\n\t"\
			"addpd	%%xmm6,%%xmm6			\n\t"\
			"addpd	%%xmm7,%%xmm7			\n\t"\
			"addpd	%%xmm4,%%xmm6			\n\t"\
			"addpd	%%xmm5,%%xmm7			\n\t"\
			"movaps	%%xmm4,%%xmm8	\n\t"\
			"movaps	%%xmm5,%%xmm9	\n\t"\
			"movaps	%%xmm0,%%xmm4			\n\t"\
			"movaps	%%xmm1,%%xmm5			\n\t"\
			"subpd	%%xmm2,%%xmm0			\n\t"\
			"subpd	%%xmm3,%%xmm1			\n\t"\
			"mulpd	0x20(%%rsi),%%xmm0		\n\t"\
			"mulpd	0x20(%%rsi),%%xmm1		\n\t"\
			"mulpd	0x30(%%rsi),%%xmm2		\n\t"\
			"mulpd	0x30(%%rsi),%%xmm3		\n\t"\
			"mulpd	0x40(%%rsi),%%xmm4		\n\t"\
			"mulpd	0x40(%%rsi),%%xmm5		\n\t"\
			"addpd	%%xmm0,%%xmm2			\n\t"\
			"addpd	%%xmm1,%%xmm3			\n\t"\
			"subpd	%%xmm4,%%xmm0			\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t"\
			"movaps	%%xmm8,%%xmm4	\n\t"\
			"movaps	%%xmm9,%%xmm5	\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"subpd	%%xmm2,%%xmm7			\n\t"\
			"addpd	%%xmm3,%%xmm3			\n\t"\
			"addpd	%%xmm2,%%xmm2			\n\t"\
			"movaps	%%xmm6,%c[__o1](%%rcx)	\n\t"\
			"movaps	%%xmm7,%c[__o4](%%rdx)	\n\t"\
			"addpd	%%xmm6,%%xmm3			\n\t"\
			"addpd	%%xmm7,%%xmm2			\n\t"\
			"movaps	%%xmm3,%c[__o4](%%rcx)	\n\t"\
			"movaps	%%xmm2,%c[__o1](%%rdx)	\n\t"\
			"subpd	%%xmm1,%%xmm4			\n\t"\
			"subpd	%%xmm0,%%xmm5			\n\t"\
			"addpd	%%xmm1,%%xmm1			\n\t"\
			"addpd	%%xmm0,%%xmm0			\n\t"\
			"movaps	%%xmm4,%c[__o2](%%rcx)	\n\t"\
			"movaps	%%xmm5,%c[__o3](%%rdx)	\n\t"\
			"addpd	%%xmm4,%%xmm1			\n\t"\
			"addpd	%%xmm5,%%xmm0			\n\t"\
			"movaps	%%xmm1,%c[__o3](%%rcx)	\n\t"\
			"movaps	%%xmm0,%c[__o2](%%rdx)	\n\t"\
			:					/* outputs: none */\
			: [__cc] "m" (Xcc)	/* m-inputs from memory addresses, e-inputs literal byte offsets here */\
			 ,[__I] "m" (XI)\
			 ,[__i0] "e" (Xi0)\
			 ,[__i1] "e" (Xi1)\
			 ,[__i2] "e" (Xi2)\
			 ,[__i3] "e" (Xi3)\
			 ,[__i4] "e" (Xi4)\
			 ,[__O] "m" (XO)\
			 ,[__o0] "e" (Xo0)\
			 ,[__o1] "e" (Xo1)\
			 ,[__o2] "e" (Xo2)\
			 ,[__o3] "e" (Xo3)\
			 ,[__o4] "e" (Xo4)\
			: "cc","memory","rax","rbx","rcx","rdx","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"		/* Clobbered registers */\
		);\
		}

	   #endif	// (OS_BITS == 32)

		#define SSE2_RADIX_15_DIF(\
			__cc0, __cc1,\
			__I,__i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA,__iB,__iC,__iD,__iE,\
			__T,__t0,__t1,__t2,__t3,__t4,__t5,__t6,__t7,__t8,__t9,__tA,__tB,__tC,__tD,__tE,\
			__O,__o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA,__oB,__oC,__oD,__oE)\
		{\
			SSE2_RADIX_05_DFT_0TWIDDLE_X1(__cc1, __I,__i0,__iC,__i9,__i6,__i3, __T,__t0,__t1,__t2,__t3,__t4);\
			SSE2_RADIX_05_DFT_0TWIDDLE_X1(__cc1, __I,__iA,__i7,__i4,__i1,__iD, __T,__t5,__t6,__t7,__t8,__t9);\
			SSE2_RADIX_05_DFT_0TWIDDLE_X1(__cc1, __I,__i5,__i2,__iE,__iB,__i8, __T,__tA,__tB,__tC,__tD,__tE);\
		\
			SSE2_RADIX_03_DFT_X1(__cc0, __T,__t0,__t5,__tA, __O,__o0,__o1,__o2);\
			SSE2_RADIX_03_DFT_X1(__cc0, __T,__t1,__t6,__tB, __O,__oD,__oE,__oC);\
			SSE2_RADIX_03_DFT_X1(__cc0, __T,__t2,__t7,__tC, __O,__o9,__oA,__oB);\
			SSE2_RADIX_03_DFT_X1(__cc0, __T,__t3,__t8,__tD, __O,__o8,__o6,__o7);\
			SSE2_RADIX_03_DFT_X1(__cc0, __T,__t4,__t9,__tE, __O,__o4,__o5,__o3);\
		/*\
			fprintf(stderr, "radix15 rad3: out[0] = %24.5f, %24.5f\n",(__O+(__o0>>4))->re,(__O+(__o0>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[1] = %24.5f, %24.5f\n",(__O+(__o1>>4))->re,(__O+(__o1>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[2] = %24.5f, %24.5f\n",(__O+(__o2>>4))->re,(__O+(__o2>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[3] = %24.5f, %24.5f\n",(__O+(__o3>>4))->re,(__O+(__o3>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[4] = %24.5f, %24.5f\n",(__O+(__o4>>4))->re,(__O+(__o4>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[5] = %24.5f, %24.5f\n",(__O+(__o5>>4))->re,(__O+(__o5>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[6] = %24.5f, %24.5f\n",(__O+(__o6>>4))->re,(__O+(__o6>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[7] = %24.5f, %24.5f\n",(__O+(__o7>>4))->re,(__O+(__o7>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[8] = %24.5f, %24.5f\n",(__O+(__o8>>4))->re,(__O+(__o8>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[9] = %24.5f, %24.5f\n",(__O+(__o9>>4))->re,(__O+(__o9>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[A] = %24.5f, %24.5f\n",(__O+(__oA>>4))->re,(__O+(__oA>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[B] = %24.5f, %24.5f\n",(__O+(__oB>>4))->re,(__O+(__oB>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[C] = %24.5f, %24.5f\n",(__O+(__oC>>4))->re,(__O+(__oC>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[D] = %24.5f, %24.5f\n",(__O+(__oD>>4))->re,(__O+(__oD>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[E] = %24.5f, %24.5f\n\n",(__O+(__oE>>4))->re,(__O+(__oE>>4)+1)->re);\
		*/\
		}

		#define SSE2_RADIX_15_DIT(\
			__cc0, __cc1,\
			__I,__i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA,__iB,__iC,__iD,__iE,\
			__T,__t0,__t1,__t2,__t3,__t4,__t5,__t6,__t7,__t8,__t9,__tA,__tB,__tC,__tD,__tE,\
			__O,__o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA,__oB,__oC,__oD,__oE)\
		{\
		/* Swap the 2nd pair of each output triplet to effect iDFT: */\
			SSE2_RADIX_03_DFT_X1(__cc0, __I,__i0,__i2,__i1, __T,__t0,__t2,__t1);\
			SSE2_RADIX_03_DFT_X1(__cc0, __I,__i8,__i7,__i6, __T,__t3,__t5,__t4);\
			SSE2_RADIX_03_DFT_X1(__cc0, __I,__iD,__iC,__iE, __T,__t6,__t8,__t7);\
			SSE2_RADIX_03_DFT_X1(__cc0, __I,__i4,__i3,__i5, __T,__t9,__tB,__tA);\
			SSE2_RADIX_03_DFT_X1(__cc0, __I,__i9,__iB,__iA, __T,__tC,__tE,__tD);\
		/*\
			fprintf(stderr, "radix15 rad3: out[0] = %24.5f, %24.5f\n",(__T+(__t0>>4))->re,(__T+(__t0>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[1] = %24.5f, %24.5f\n",(__T+(__t1>>4))->re,(__T+(__t1>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[2] = %24.5f, %24.5f\n",(__T+(__t2>>4))->re,(__T+(__t2>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[3] = %24.5f, %24.5f\n",(__T+(__t3>>4))->re,(__T+(__t3>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[4] = %24.5f, %24.5f\n",(__T+(__t4>>4))->re,(__T+(__t4>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[5] = %24.5f, %24.5f\n",(__T+(__t5>>4))->re,(__T+(__t5>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[6] = %24.5f, %24.5f\n",(__T+(__t6>>4))->re,(__T+(__t6>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[7] = %24.5f, %24.5f\n",(__T+(__t7>>4))->re,(__T+(__t7>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[8] = %24.5f, %24.5f\n",(__T+(__t8>>4))->re,(__T+(__t8>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[9] = %24.5f, %24.5f\n",(__T+(__t9>>4))->re,(__T+(__t9>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[A] = %24.5f, %24.5f\n",(__T+(__tA>>4))->re,(__T+(__tA>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[B] = %24.5f, %24.5f\n",(__T+(__tB>>4))->re,(__T+(__tB>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[C] = %24.5f, %24.5f\n",(__T+(__tC>>4))->re,(__T+(__tC>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[D] = %24.5f, %24.5f\n",(__T+(__tD>>4))->re,(__T+(__tD>>4)+1)->re);\
			fprintf(stderr, "radix15 rad3: out[E] = %24.5f, %24.5f\n\n",(__T+(__tE>>4))->re,(__T+(__tE>>4)+1)->re);\
		exit(0);\
		*/\
		/* Output perm here is 0123456789abcde --> 05a6b1c2738d9e4: */\
			SSE2_RADIX_05_DFT_0TWIDDLE_X1(__cc1, __T,__t0,__t3,__t6,__t9,__tC, __O,__o0,__o6,__oC,__o3,__o9);\
			SSE2_RADIX_05_DFT_0TWIDDLE_X1(__cc1, __T,__t1,__t4,__t7,__tA,__tD, __O,__o5,__oB,__o2,__o8,__oE);\
			SSE2_RADIX_05_DFT_0TWIDDLE_X1(__cc1, __T,__t2,__t5,__t8,__tB,__tE, __O,__oA,__o1,__o7,__oD,__o4);\
		}

	  #endif

	#else	// USE_64BIT_ASM_STYLE = True; 16-register version:

		#define SSE2_RADIX_15_DIF_X2(\
			__cc0, __cc1, __two,\
			__i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA,__iB,__iC,__iD,__iE,\
			__s0,__s1,__s2,__s3,__s4,__s5,__s6,__s7,__s8,__s9,__sA,__sB,__sC,__sD,__sE,\
			__o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA,__oB,__oC,__oD,__oE,\
					__j0,__j1,__j2,__j3,__j4,__j5,__j6,__j7,__j8,__j9,__jA,__jB,__jC,__jD,__jE,\
					__t0,__t1,__t2,__t3,__t4,__t5,__t6,__t7,__t8,__t9,__tA,__tB,__tC,__tD,__tE,\
					__u0,__u1,__u2,__u3,__u4,__u5,__u6,__u7,__u8,__u9,__uA,__uB,__uC,__uD,__uE)\
		{\
			SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1,__two, __i0,__iC,__i9,__i6,__i3, __s0,__s1,__s2,__s3,__s4,	__j0,__jC,__j9,__j6,__j3, __t0,__t1,__t2,__t3,__t4);\
			SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1,__two, __iA,__i7,__i4,__i1,__iD, __s5,__s6,__s7,__s8,__s9,	__jA,__j7,__j4,__j1,__jD, __t5,__t6,__t7,__t8,__t9);\
			SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1,__two, __i5,__i2,__iE,__iB,__i8, __sA,__sB,__sC,__sD,__sE,	__j5,__j2,__jE,__jB,__j8, __tA,__tB,__tC,__tD,__tE);\
		\
			SSE2_RADIX_03_DFT_X2(__cc0, __s0,__s5,__sA, __o0,__o1,__o2,		__t0,__t5,__tA, __u0,__u1,__u2);\
			SSE2_RADIX_03_DFT_X2(__cc0, __s1,__s6,__sB, __oD,__oE,__oC,		__t1,__t6,__tB, __uD,__uE,__uC);\
			SSE2_RADIX_03_DFT_X2(__cc0, __s2,__s7,__sC, __o9,__oA,__oB,		__t2,__t7,__tC, __u9,__uA,__uB);\
			SSE2_RADIX_03_DFT_X2(__cc0, __s3,__s8,__sD, __o8,__o6,__o7,		__t3,__t8,__tD, __u8,__u6,__u7);\
			SSE2_RADIX_03_DFT_X2(__cc0, __s4,__s9,__sE, __o4,__o5,__o3,		__t4,__t9,__tE, __u4,__u5,__u3);\
		}

		#define SSE2_RADIX_15_DIT_X2(\
			__cc0, __cc1,__two,\
			__i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA,__iB,__iC,__iD,__iE,\
			__s0,__s1,__s2,__s3,__s4,__s5,__s6,__s7,__s8,__s9,__sA,__sB,__sC,__sD,__sE,\
			__o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA,__oB,__oC,__oD,__oE,\
					__j0,__j1,__j2,__j3,__j4,__j5,__j6,__j7,__j8,__j9,__jA,__jB,__jC,__jD,__jE,\
					__t0,__t1,__t2,__t3,__t4,__t5,__t6,__t7,__t8,__t9,__tA,__tB,__tC,__tD,__tE,\
					__u0,__u1,__u2,__u3,__u4,__u5,__u6,__u7,__u8,__u9,__uA,__uB,__uC,__uD,__uE)\
		{\
		/* Swap the 2nd pair of each output triplet to effect iDFT: */\
			SSE2_RADIX_03_DFT_X2(__cc0, __i0,__i2,__i1, __s0,__s2,__s1,		__j0,__j2,__j1, __t0,__t2,__t1);\
			SSE2_RADIX_03_DFT_X2(__cc0, __i8,__i7,__i6, __s3,__s5,__s4,		__j8,__j7,__j6, __t3,__t5,__t4);\
			SSE2_RADIX_03_DFT_X2(__cc0, __iD,__iC,__iE, __s6,__s8,__s7,		__jD,__jC,__jE, __t6,__t8,__t7);\
			SSE2_RADIX_03_DFT_X2(__cc0, __i4,__i3,__i5, __s9,__sB,__sA,		__j4,__j3,__j5, __t9,__tB,__tA);\
			SSE2_RADIX_03_DFT_X2(__cc0, __i9,__iB,__iA, __sC,__sE,__sD,		__j9,__jB,__jA, __tC,__tE,__tD);\
		\
		/* Output perm here is 0123456789abcde --> 05a6b1c2738d9e4: */\
			SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1,__two, __s0,__s3,__s6,__s9,__sC, __o0,__o6,__oC,__o3,__o9,	__t0,__t3,__t6,__t9,__tC, __u0,__u6,__uC,__u3,__u9);\
			SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1,__two, __s1,__s4,__s7,__sA,__sD, __o5,__oB,__o2,__o8,__oE,	__t1,__t4,__t7,__tA,__tD, __u5,__uB,__u2,__u8,__uE);\
			SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1,__two, __s2,__s5,__s8,__sB,__sE, __oA,__o1,__o7,__oD,__o4,	__t2,__t5,__t8,__tB,__tE, __uA,__u1,__u7,__uD,__u4);\
		}

	#endif	// USE_64BIT_ASM_STYLE ?

#endif	/* radix15_sse_macro_h_included */

