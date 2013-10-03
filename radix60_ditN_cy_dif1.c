/*******************************************************************************
*                                                                              *
*   (C) 1997-2013 by Ernst W. Mayer.                                           *
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

#ifdef MULTITHREAD
	#ifndef USE_PTHREAD
		#error Pthreads is only thread model currently supported!
	#endif
#endif

#define ERR_CHECK_ALL	/* #define this to do ROE checking of all convolution outputs, rather than just every 36th one */
#if defined(USE_AVX) && !defined(ERR_CHECK_ALL)
	#error ERR_CHECK_ALL *required* for AVX-mode builds!
#endif

/* Use for toggling higher-accuracy version of the twiddles computation */
//#define HIACC 0	<*** prefer to set via compile-time flag; default is FALSE [= LOACC]

#define EPS 1e-10

	/* Scalar-data test macros for comparing with the SSE2 radix-15 DFTs: */
	#define RADIX_15_DIF_B(\
		__s,__c3m1,\
		__cn1, __cn2, __ss3, __sn1, __sn2,\
		__Ar00,__Ai00,__Ar01,__Ai01,__Ar02,__Ai02,__Ar03,__Ai03,__Ar04,__Ai04,__Ar05,__Ai05,__Ar06,__Ai06,__Ar07,__Ai07,__Ar08,__Ai08,__Ar09,__Ai09,__Ar10,__Ai10,__Ar11,__Ai11,__Ar12,__Ai12,__Ar13,__Ai13,__Ar14,__Ai14,\
		__tr00,__ti00,__tr01,__ti01,__tr02,__ti02,__tr03,__ti03,__tr04,__ti04,__tr05,__ti05,__tr06,__ti06,__tr07,__ti07,__tr08,__ti08,__tr09,__ti09,__tr10,__ti10,__tr11,__ti11,__tr12,__ti12,__tr13,__ti13,__tr14,__ti14,\
		__Br00,__Bi00,__Br01,__Bi01,__Br02,__Bi02,__Br03,__Bi03,__Br04,__Bi04,__Br05,__Bi05,__Br06,__Bi06,__Br07,__Bi07,__Br08,__Bi08,__Br09,__Bi09,__Br10,__Bi10,__Br11,__Bi11,__Br12,__Bi12,__Br13,__Bi13,__Br14,__Bi14,\
		__rt,__it)\
	{\
		double x0,x1,x2,x3,x4,x5;\
		RADIX_05_DFT(__cn1,__cn2,__ss3,__sn1,__sn2,__Ar00,__Ai00,__Ar12,__Ai12,__Ar09,__Ai09,__Ar06,__Ai06,__Ar03,__Ai03,__tr00,__ti00,__tr01,__ti01,__tr02,__ti02,__tr03,__ti03,__tr04,__ti04,__rt,__it);\
		RADIX_05_DFT(__cn1,__cn2,__ss3,__sn1,__sn2,__Ar10,__Ai10,__Ar07,__Ai07,__Ar04,__Ai04,__Ar01,__Ai01,__Ar13,__Ai13,__tr05,__ti05,__tr06,__ti06,__tr07,__ti07,__tr08,__ti08,__tr09,__ti09,__rt,__it);\
		RADIX_05_DFT(__cn1,__cn2,__ss3,__sn1,__sn2,__Ar05,__Ai05,__Ar02,__Ai02,__Ar14,__Ai14,__Ar11,__Ai11,__Ar08,__Ai08,__tr10,__ti10,__tr11,__ti11,__tr12,__ti12,__tr13,__ti13,__tr14,__ti14,__rt,__it);\
	\
		RADIX_03_DFT(__s,__c3m1,__tr00,__ti00,__tr05,__ti05,__tr10,__ti10,x0,x1,x2,x3,x4,x5,__Br00,__Bi00,__Br01,__Bi01,__Br02,__Bi02);\
		RADIX_03_DFT(__s,__c3m1,__tr01,__ti01,__tr06,__ti06,__tr11,__ti11,x0,x1,x2,x3,x4,x5,__Br13,__Bi13,__Br14,__Bi14,__Br12,__Bi12);\
		RADIX_03_DFT(__s,__c3m1,__tr02,__ti02,__tr07,__ti07,__tr12,__ti12,x0,x1,x2,x3,x4,x5,__Br09,__Bi09,__Br10,__Bi10,__Br11,__Bi11);\
		RADIX_03_DFT(__s,__c3m1,__tr03,__ti03,__tr08,__ti08,__tr13,__ti13,x0,x1,x2,x3,x4,x5,__Br08,__Bi08,__Br06,__Bi06,__Br07,__Bi07);\
		RADIX_03_DFT(__s,__c3m1,__tr04,__ti04,__tr09,__ti09,__tr14,__ti14,x0,x1,x2,x3,x4,x5,__Br04,__Bi04,__Br05,__Bi05,__Br03,__Bi03);\
	}

	#define RADIX_15_DIT_B(\
		__s,__c3m1,\
		__cn1, __cn2, __ss3, __sn1, __sn2,\
		__Ar00,__Ai00,__Ar01,__Ai01,__Ar02,__Ai02,__Ar03,__Ai03,__Ar04,__Ai04,__Ar05,__Ai05,__Ar06,__Ai06,__Ar07,__Ai07,__Ar08,__Ai08,__Ar09,__Ai09,__Ar10,__Ai10,__Ar11,__Ai11,__Ar12,__Ai12,__Ar13,__Ai13,__Ar14,__Ai14,\
		__tr00,__ti00,__tr01,__ti01,__tr02,__ti02,__tr03,__ti03,__tr04,__ti04,__tr05,__ti05,__tr06,__ti06,__tr07,__ti07,__tr08,__ti08,__tr09,__ti09,__tr10,__ti10,__tr11,__ti11,__tr12,__ti12,__tr13,__ti13,__tr14,__ti14,\
		__Br00,__Bi00,__Br01,__Bi01,__Br02,__Bi02,__Br03,__Bi03,__Br04,__Bi04,__Br05,__Bi05,__Br06,__Bi06,__Br07,__Bi07,__Br08,__Bi08,__Br09,__Bi09,__Br10,__Bi10,__Br11,__Bi11,__Br12,__Bi12,__Br13,__Bi13,__Br14,__Bi14,\
		__rt,__it)\
	{\
		double x0,x1,x2,x3,x4,x5;\
		/* Swap the 2nd pair of each output triplet to effect iDFT: */\
		RADIX_03_DFT(__s,__c3m1,__Ar00,__Ai00,__Ar02,__Ai02,__Ar01,__Ai01,x0,x1,x2,x3,x4,x5,__tr00,__ti00,__tr02,__ti02,__tr01,__ti01);\
		RADIX_03_DFT(__s,__c3m1,__Ar08,__Ai08,__Ar07,__Ai07,__Ar06,__Ai06,x0,x1,x2,x3,x4,x5,__tr03,__ti03,__tr05,__ti05,__tr04,__ti04);\
		RADIX_03_DFT(__s,__c3m1,__Ar13,__Ai13,__Ar12,__Ai12,__Ar14,__Ai14,x0,x1,x2,x3,x4,x5,__tr06,__ti06,__tr08,__ti08,__tr07,__ti07);\
		RADIX_03_DFT(__s,__c3m1,__Ar04,__Ai04,__Ar03,__Ai03,__Ar05,__Ai05,x0,x1,x2,x3,x4,x5,__tr09,__ti09,__tr11,__ti11,__tr10,__ti10);\
		RADIX_03_DFT(__s,__c3m1,__Ar09,__Ai09,__Ar11,__Ai11,__Ar10,__Ai10,x0,x1,x2,x3,x4,x5,__tr12,__ti12,__tr14,__ti14,__tr13,__ti13);\
		/* Output perm here is 0123456789abcde --> 05a6b1c2738d9e4: */\
		RADIX_05_DFT(__cn1,__cn2,__ss3,__sn1,__sn2,__tr00,__ti00,__tr03,__ti03,__tr06,__ti06,__tr09,__ti09,__tr12,__ti12,__Br00,__Bi00,__Br06,__Bi06,__Br12,__Bi12,__Br03,__Bi03,__Br09,__Bi09,__rt,__it);\
		RADIX_05_DFT(__cn1,__cn2,__ss3,__sn1,__sn2,__tr01,__ti01,__tr04,__ti04,__tr07,__ti07,__tr10,__ti10,__tr13,__ti13,__Br05,__Bi05,__Br11,__Bi11,__Br02,__Bi02,__Br08,__Bi08,__Br14,__Bi14,__rt,__it);\
		RADIX_05_DFT(__cn1,__cn2,__ss3,__sn1,__sn2,__tr02,__ti02,__tr05,__ti05,__tr08,__ti08,__tr11,__ti11,__tr14,__ti14,__Br10,__Bi10,__Br01,__Bi01,__Br07,__Bi07,__Br13,__Bi13,__Br04,__Bi04,__rt,__it);\
	}

#ifdef USE_SSE2

  // For Mersenne-mod we need (16 [SSE2] or 64 [AVX]) + 4 added slots for the half_arr lookup tables.
  // For Fermat-mod in AVX mode we need RADIX*4 = 240 [if HIACC] or 12 [if not] slots for the compact
  // negacyclic-roots chained-multiply scheme.
  // Add larger number in each case - i.e. max(68,240) = 240 if AVX+HIACC, max(68,12) = 68 if AVX+LOACC, 20 if SSE2
  // to (half_arr_offset60 + RADIX) to get required value of radix60_creals_in_local_store:
  #ifdef USE_AVX
	const int half_arr_offset60 = 340;	// + RADIX = 400;  Used for thread local-storage-integrity checking
   #if HIACC
	const int radix60_creals_in_local_store = 640;	// AVX+HIACC: 400 + 240
   #else
	const int radix60_creals_in_local_store = 468;	// AVX+LOACC: 400 + 68
   #endif
  #else
	const int half_arr_offset60 = 370;	// + RADIX = 430; Used for thread local-storage-integrity checking
	const int radix60_creals_in_local_store = 452;	// (half_arr_offset60 + RADIX) + 20 and round up to nearest multiple of 4
  #endif

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

	#if (OS_BITS == 32) || !defined(USE_64BIT_ASM_STYLE)	// Default is to use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just xmm0-7.

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
		#define USE_LITERAL_BYTE_OFFSETS	// undef this to revert to the 1st set of macros below:
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

		// SSE2_RADIX_03_DFT_X2 also used by other radices, to be found in sse2_macro*.h

	  #if USE_AVX

		/* This is slightly faster - for unknown reasons - than the SSE2_RADIX_05_DFT_0TWID_X2 defined in sse_macro_gcc64.h: */
		#define SSE2_RADIX_05_DFT_0TWIDDLE_X2(Xcc1, Xi0,Xi1,Xi2,Xi3,Xi4, Xo0,Xo1,Xo2,Xo3,Xo4, Xj0,Xj1,Xj2,Xj3,Xj4, Xu0,Xu1,Xu2,Xu3,Xu4)\
		{\
		__asm__ volatile (\
			"movq	%[__i0],%%rsi				\n\t	movq	%[__j0],%%r10		\n\t"\
			"movq	%[__i1],%%rax				\n\t	movq	%[__j1],%%r11		\n\t"\
			"movq	%[__i2],%%rbx				\n\t	movq	%[__j2],%%r12		\n\t"\
			"movq	%[__i3],%%rcx				\n\t	movq	%[__j3],%%r13		\n\t"\
			"movq	%[__i4],%%rdx				\n\t	movq	%[__j4],%%r14		\n\t"\
			"movq	%[__o0],%%rdi				\n\t	movq	%[__u0],%%r15		\n\t"\
			"vmovaps	    (%%rax),%%ymm0		\n\t	vmovaps	    (%%r11),%%ymm8 	\n\t"\
			"vmovaps	0x20(%%rax),%%ymm1		\n\t	vmovaps	0x20(%%r11),%%ymm9 	\n\t"\
			"vmovaps	    (%%rbx),%%ymm2		\n\t	vmovaps	    (%%r12),%%ymm10	\n\t"\
			"vmovaps	0x20(%%rbx),%%ymm3		\n\t	vmovaps	0x20(%%r12),%%ymm11	\n\t"\
			"vmovaps	    (%%rcx),%%ymm4		\n\t	vmovaps	    (%%r13),%%ymm12	\n\t"\
			"vmovaps	0x20(%%rcx),%%ymm5		\n\t	vmovaps	0x20(%%r13),%%ymm13	\n\t"\
			"vmovaps	    (%%rdx),%%ymm6		\n\t	vmovaps	    (%%r14),%%ymm14	\n\t"\
			"vmovaps	0x20(%%rdx),%%ymm7		\n\t	vmovaps	0x20(%%r14),%%ymm15	\n\t"\
			"vsubpd	%%ymm6,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm14,%%ymm8 ,%%ymm8 		\n\t"\
			"vsubpd	%%ymm7,%%ymm1,%%ymm1		\n\t	vsubpd	%%ymm15,%%ymm9 ,%%ymm9 		\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
			"vaddpd	%%ymm0,%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm8 ,%%ymm14,%%ymm14		\n\t"\
			"vaddpd	%%ymm1,%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm9 ,%%ymm15,%%ymm15		\n\t"\
			"vsubpd	%%ymm4,%%ymm2,%%ymm2		\n\t	vsubpd	%%ymm12,%%ymm10,%%ymm10		\n\t"\
			"vsubpd	%%ymm5,%%ymm3,%%ymm3		\n\t	vsubpd	%%ymm13,%%ymm11,%%ymm11		\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
			"vaddpd	%%ymm2,%%ymm4,%%ymm4		\n\t	vaddpd	%%ymm10,%%ymm12,%%ymm12		\n\t"\
			"vaddpd	%%ymm3,%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm11,%%ymm13,%%ymm13		\n\t"\
			"movq	%[__cc1],%%rax				\n\t"\
			"vsubpd	%%ymm4,%%ymm6,%%ymm6		\n\t	vsubpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
			"vsubpd	%%ymm5,%%ymm7,%%ymm7		\n\t	vsubpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4		\n\t	vaddpd	%%ymm12,%%ymm12,%%ymm12		\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm13,%%ymm13,%%ymm13		\n\t"\
			"vaddpd	%%ymm6,%%ymm4,%%ymm4		\n\t	vaddpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
			"vaddpd	%%ymm7,%%ymm5,%%ymm5		\n\t	vaddpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
			"vaddpd	    (%%rsi),%%ymm4,%%ymm4	\n\t	vaddpd	    (%%r10),%%ymm12,%%ymm12	\n\t"\
			"vaddpd	0x20(%%rsi),%%ymm5,%%ymm5	\n\t	vaddpd	0x20(%%r10),%%ymm13,%%ymm13	\n\t"\
			"vmovaps	%%ymm4,    (%%rdi)		\n\t	vmovaps	%%ymm12,    (%%r15)		\n\t"\
			"vmovaps	%%ymm5,0x20(%%rdi)		\n\t	vmovaps	%%ymm13,0x20(%%r15)		\n\t"\
			"vmulpd	0x20(%%rax),%%ymm6,%%ymm6	\n\t	vmulpd	0x20(%%rax),%%ymm14,%%ymm14	\n\t"\
			"vmulpd	0x20(%%rax),%%ymm7,%%ymm7	\n\t	vmulpd	0x20(%%rax),%%ymm15,%%ymm15	\n\t"\
			"vsubpd	    (%%rsi),%%ymm4,%%ymm4	\n\t	vsubpd	    (%%r10),%%ymm12,%%ymm12	\n\t"\
			"vsubpd	0x20(%%rsi),%%ymm5,%%ymm5	\n\t	vsubpd	0x20(%%r10),%%ymm13,%%ymm13	\n\t"\
			"vmulpd	    (%%rax),%%ymm4,%%ymm4	\n\t	vmulpd	    (%%rax),%%ymm12,%%ymm12	\n\t"\
			"vmulpd	    (%%rax),%%ymm5,%%ymm5	\n\t	vmulpd	    (%%rax),%%ymm13,%%ymm13	\n\t"\
			"vaddpd	    (%%rdi),%%ymm4,%%ymm4	\n\t	vaddpd	    (%%r15),%%ymm12,%%ymm12	\n\t"\
			"vaddpd	0x20(%%rdi),%%ymm5,%%ymm5	\n\t	vaddpd	0x20(%%r15),%%ymm13,%%ymm13	\n\t"\
			"vsubpd	%%ymm6,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm14,%%ymm12,%%ymm12		\n\t"\
			"vsubpd	%%ymm7,%%ymm5,%%ymm5		\n\t	vsubpd	%%ymm15,%%ymm13,%%ymm13		\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
			"vaddpd	%%ymm4,%%ymm6,%%ymm6		\n\t	vaddpd	%%ymm12,%%ymm14,%%ymm14		\n\t"\
			"vaddpd	%%ymm5,%%ymm7,%%ymm7		\n\t	vaddpd	%%ymm13,%%ymm15,%%ymm15		\n\t"\
			"vmovaps	%%ymm4,    (%%rsi)		\n\t	vmovaps	%%ymm12,    (%%r10)		\n\t"\
			"vmovaps	%%ymm5,0x20(%%rsi)		\n\t	vmovaps	%%ymm13,0x20(%%r10)		\n\t"\
			"vmovaps	%%ymm0,%%ymm4			\n\t	vmovaps	%%ymm8 ,%%ymm12			\n\t"\
			"vmovaps	%%ymm1,%%ymm5			\n\t	vmovaps	%%ymm9 ,%%ymm13			\n\t"\
			"vsubpd	%%ymm2,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm10,%%ymm8,%%ymm8 		\n\t"\
			"vsubpd	%%ymm3,%%ymm1,%%ymm1		\n\t	vsubpd	%%ymm11,%%ymm9,%%ymm9 		\n\t"\
			"vmulpd	0x40(%%rax),%%ymm0,%%ymm0	\n\t	vmulpd	0x40(%%rax),%%ymm8 ,%%ymm8 	\n\t"\
			"vmulpd	0x40(%%rax),%%ymm1,%%ymm1	\n\t	vmulpd	0x40(%%rax),%%ymm9 ,%%ymm9 	\n\t"\
			"vmulpd	0x60(%%rax),%%ymm2,%%ymm2	\n\t	vmulpd	0x60(%%rax),%%ymm10,%%ymm10	\n\t"\
			"vmulpd	0x60(%%rax),%%ymm3,%%ymm3	\n\t	vmulpd	0x60(%%rax),%%ymm11,%%ymm11	\n\t"\
			"vmulpd	0x80(%%rax),%%ymm4,%%ymm4	\n\t	vmulpd	0x80(%%rax),%%ymm12,%%ymm12	\n\t"\
			"vmulpd	0x80(%%rax),%%ymm5,%%ymm5	\n\t	vmulpd	0x80(%%rax),%%ymm13,%%ymm13	\n\t"\
			"vaddpd	%%ymm0,%%ymm2,%%ymm2		\n\t	vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
			"vaddpd	%%ymm1,%%ymm3,%%ymm3		\n\t	vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0		\n\t	vsubpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1		\n\t	vsubpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
			"vmovaps	    (%%rsi),%%ymm4		\n\t	vmovaps	    (%%r10),%%ymm12		\n\t"\
			"vmovaps	0x20(%%rsi),%%ymm5		\n\t	vmovaps	0x20(%%r10),%%ymm13		\n\t"\
			"movq	%[__o1],%%rax				\n\t	movq	%[__u1],%%r11			\n\t"\
			"movq	%[__o4],%%rdx				\n\t	movq	%[__u4],%%r14			\n\t"\
			"vsubpd	%%ymm3,%%ymm6,%%ymm6		\n\t	vsubpd	%%ymm11,%%ymm14,%%ymm14		\n\t"\
			"vsubpd	%%ymm2,%%ymm7,%%ymm7		\n\t	vsubpd	%%ymm10,%%ymm15,%%ymm15		\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3		\n\t	vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2		\n\t	vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
			"vmovaps	%%ymm6,    (%%rax)		\n\t	vmovaps	%%ymm14,    (%%r11)		\n\t"\
			"vmovaps	%%ymm7,0x20(%%rdx)		\n\t	vmovaps	%%ymm15,0x20(%%r14)		\n\t"\
			"vaddpd	%%ymm6,%%ymm3,%%ymm3		\n\t	vaddpd	%%ymm14,%%ymm11,%%ymm11		\n\t"\
			"vaddpd	%%ymm7,%%ymm2,%%ymm2		\n\t	vaddpd	%%ymm15,%%ymm10,%%ymm10		\n\t"\
			"vmovaps	%%ymm3,    (%%rdx)		\n\t	vmovaps	%%ymm11,    (%%r14)		\n\t"\
			"vmovaps	%%ymm2,0x20(%%rax)		\n\t	vmovaps	%%ymm10,0x20(%%r11)		\n\t"\
			"movq	%[__o2],%%rbx				\n\t	movq	%[__u2],%%r12			\n\t"\
			"movq	%[__o3],%%rcx				\n\t	movq	%[__u3],%%r13			\n\t"\
			"vsubpd	%%ymm1,%%ymm4,%%ymm4		\n\t	vsubpd	%%ymm9 ,%%ymm12,%%ymm12		\n\t"\
			"vsubpd	%%ymm0,%%ymm5,%%ymm5		\n\t	vsubpd	%%ymm8 ,%%ymm13,%%ymm13		\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 		\n\t"\
			"vaddpd	%%ymm0,%%ymm0,%%ymm0		\n\t	vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 		\n\t"\
			"vmovaps	%%ymm4,    (%%rbx)		\n\t	vmovaps	%%ymm12,    (%%r12)		\n\t"\
			"vmovaps	%%ymm5,0x20(%%rcx)		\n\t	vmovaps	%%ymm13,0x20(%%r13)		\n\t"\
			"vaddpd	%%ymm4,%%ymm1,%%ymm1		\n\t	vaddpd	%%ymm12,%%ymm9,%%ymm9 		\n\t"\
			"vaddpd	%%ymm5,%%ymm0,%%ymm0		\n\t	vaddpd	%%ymm13,%%ymm8,%%ymm8 		\n\t"\
			"vmovaps	%%ymm1,    (%%rcx)		\n\t	vmovaps	%%ymm9 ,    (%%r13)		\n\t"\
			"vmovaps	%%ymm0,0x20(%%rbx)		\n\t	vmovaps	%%ymm8 ,0x20(%%r12)		\n\t"\
			:					/* outputs: none */\
			: [__cc1] "m" (Xcc1)	/* All inputs from memory addresses here */\
			 ,[__i0] "m" (Xi0)\
			 ,[__i1] "m" (Xi1)\
			 ,[__i2] "m" (Xi2)\
			 ,[__i3] "m" (Xi3)\
			 ,[__i4] "m" (Xi4)\
			 ,[__o0] "m" (Xo0)\
			 ,[__o1] "m" (Xo1)\
			 ,[__o2] "m" (Xo2)\
			 ,[__o3] "m" (Xo3)\
			 ,[__o4] "m" (Xo4)\
			 ,[__j0] "m" (Xj0)\
			 ,[__j1] "m" (Xj1)\
			 ,[__j2] "m" (Xj2)\
			 ,[__j3] "m" (Xj3)\
			 ,[__j4] "m" (Xj4)\
			 ,[__u0] "m" (Xu0)\
			 ,[__u1] "m" (Xu1)\
			 ,[__u2] "m" (Xu2)\
			 ,[__u3] "m" (Xu3)\
			 ,[__u4] "m" (Xu4)\
			: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
		);\
		}

	  #else	// USE_SSE2

		/* This is slightly faster - for unknown reasons - than the SSE2_RADIX_05_DFT_0TWID_X2 defined in sse_macro_gcc64.h: */
		#define SSE2_RADIX_05_DFT_0TWIDDLE_X2(Xcc1, Xi0,Xi1,Xi2,Xi3,Xi4, Xo0,Xo1,Xo2,Xo3,Xo4, Xj0,Xj1,Xj2,Xj3,Xj4, Xu0,Xu1,Xu2,Xu3,Xu4)\
		{\
		__asm__ volatile (\
			"movq	%[__i0],%%rsi		\n\t	movq	%[__j0],%%r10		\n\t"\
			"movq	%[__i1],%%rax		\n\t	movq	%[__j1],%%r11		\n\t"\
			"movq	%[__i2],%%rbx		\n\t	movq	%[__j2],%%r12		\n\t"\
			"movq	%[__i3],%%rcx		\n\t	movq	%[__j3],%%r13		\n\t"\
			"movq	%[__i4],%%rdx		\n\t	movq	%[__j4],%%r14		\n\t"\
			"movq	%[__o0],%%rdi		\n\t	movq	%[__u0],%%r15		\n\t"\
			"movaps	    (%%rax),%%xmm0	\n\t	movaps	    (%%r11),%%xmm8 	\n\t"\
			"movaps	0x10(%%rax),%%xmm1	\n\t	movaps	0x10(%%r11),%%xmm9 	\n\t"\
			"movaps	    (%%rbx),%%xmm2	\n\t	movaps	    (%%r12),%%xmm10	\n\t"\
			"movaps	0x10(%%rbx),%%xmm3	\n\t	movaps	0x10(%%r12),%%xmm11	\n\t"\
			"movaps	    (%%rcx),%%xmm4	\n\t	movaps	    (%%r13),%%xmm12	\n\t"\
			"movaps	0x10(%%rcx),%%xmm5	\n\t	movaps	0x10(%%r13),%%xmm13	\n\t"\
			"movaps	    (%%rdx),%%xmm6	\n\t	movaps	    (%%r14),%%xmm14	\n\t"\
			"movaps	0x10(%%rdx),%%xmm7	\n\t	movaps	0x10(%%r14),%%xmm15	\n\t"\
			"subpd	%%xmm6,%%xmm0		\n\t	subpd	%%xmm14,%%xmm8 		\n\t"\
			"subpd	%%xmm7,%%xmm1		\n\t	subpd	%%xmm15,%%xmm9 		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t	addpd	%%xmm14,%%xmm14		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t	addpd	%%xmm15,%%xmm15		\n\t"\
			"addpd	%%xmm0,%%xmm6		\n\t	addpd	%%xmm8 ,%%xmm14		\n\t"\
			"addpd	%%xmm1,%%xmm7		\n\t	addpd	%%xmm9 ,%%xmm15		\n\t"\
			"subpd	%%xmm4,%%xmm2		\n\t	subpd	%%xmm12,%%xmm10		\n\t"\
			"subpd	%%xmm5,%%xmm3		\n\t	subpd	%%xmm13,%%xmm11		\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t	addpd	%%xmm12,%%xmm12		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t	addpd	%%xmm13,%%xmm13		\n\t"\
			"addpd	%%xmm2,%%xmm4		\n\t	addpd	%%xmm10,%%xmm12		\n\t"\
			"addpd	%%xmm3,%%xmm5		\n\t	addpd	%%xmm11,%%xmm13		\n\t"\
			"movq	%[__cc1],%%rax		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t	subpd	%%xmm12,%%xmm14		\n\t"\
			"subpd	%%xmm5,%%xmm7		\n\t	subpd	%%xmm13,%%xmm15		\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t	addpd	%%xmm12,%%xmm12		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t	addpd	%%xmm13,%%xmm13		\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t	addpd	%%xmm14,%%xmm12		\n\t"\
			"addpd	%%xmm7,%%xmm5		\n\t	addpd	%%xmm15,%%xmm13		\n\t"\
			"addpd	    (%%rsi),%%xmm4	\n\t	addpd	    (%%r10),%%xmm12	\n\t"\
			"addpd	0x10(%%rsi),%%xmm5	\n\t	addpd	0x10(%%r10),%%xmm13	\n\t"\
			"movaps	%%xmm4,    (%%rdi)	\n\t	movaps	%%xmm12,    (%%r15)	\n\t"\
			"movaps	%%xmm5,0x10(%%rdi)	\n\t	movaps	%%xmm13,0x10(%%r15)	\n\t"\
			"mulpd	0x10(%%rax),%%xmm6	\n\t	mulpd	0x10(%%rax),%%xmm14	\n\t"\
			"mulpd	0x10(%%rax),%%xmm7	\n\t	mulpd	0x10(%%rax),%%xmm15	\n\t"\
			"subpd	    (%%rsi),%%xmm4	\n\t	subpd	    (%%r10),%%xmm12	\n\t"\
			"subpd	0x10(%%rsi),%%xmm5	\n\t	subpd	0x10(%%r10),%%xmm13	\n\t"\
			"mulpd	    (%%rax),%%xmm4	\n\t	mulpd	    (%%rax),%%xmm12	\n\t"\
			"mulpd	    (%%rax),%%xmm5	\n\t	mulpd	    (%%rax),%%xmm13	\n\t"\
			"addpd	    (%%rdi),%%xmm4	\n\t	addpd	    (%%r15),%%xmm12	\n\t"\
			"addpd	0x10(%%rdi),%%xmm5	\n\t	addpd	0x10(%%r15),%%xmm13	\n\t"\
			"subpd	%%xmm6,%%xmm4		\n\t	subpd	%%xmm14,%%xmm12		\n\t"\
			"subpd	%%xmm7,%%xmm5		\n\t	subpd	%%xmm15,%%xmm13		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t	addpd	%%xmm14,%%xmm14		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t	addpd	%%xmm15,%%xmm15		\n\t"\
			"addpd	%%xmm4,%%xmm6		\n\t	addpd	%%xmm12,%%xmm14		\n\t"\
			"addpd	%%xmm5,%%xmm7		\n\t	addpd	%%xmm13,%%xmm15		\n\t"\
			"movaps	%%xmm4,    (%%rsi)	\n\t	movaps	%%xmm12,    (%%r10)	\n\t"\
			"movaps	%%xmm5,0x10(%%rsi)	\n\t	movaps	%%xmm13,0x10(%%r10)	\n\t"\
			"movaps	%%xmm0,%%xmm4		\n\t	movaps	%%xmm8 ,%%xmm12		\n\t"\
			"movaps	%%xmm1,%%xmm5		\n\t	movaps	%%xmm9 ,%%xmm13		\n\t"\
			"subpd	%%xmm2,%%xmm0		\n\t	subpd	%%xmm10,%%xmm8 		\n\t"\
			"subpd	%%xmm3,%%xmm1		\n\t	subpd	%%xmm11,%%xmm9 		\n\t"\
			"mulpd	0x20(%%rax),%%xmm0	\n\t	mulpd	0x20(%%rax),%%xmm8 	\n\t"\
			"mulpd	0x20(%%rax),%%xmm1	\n\t	mulpd	0x20(%%rax),%%xmm9 	\n\t"\
			"mulpd	0x30(%%rax),%%xmm2	\n\t	mulpd	0x30(%%rax),%%xmm10	\n\t"\
			"mulpd	0x30(%%rax),%%xmm3	\n\t	mulpd	0x30(%%rax),%%xmm11	\n\t"\
			"mulpd	0x40(%%rax),%%xmm4	\n\t	mulpd	0x40(%%rax),%%xmm12	\n\t"\
			"mulpd	0x40(%%rax),%%xmm5	\n\t	mulpd	0x40(%%rax),%%xmm13	\n\t"\
			"addpd	%%xmm0,%%xmm2		\n\t	addpd	%%xmm8 ,%%xmm10		\n\t"\
			"addpd	%%xmm1,%%xmm3		\n\t	addpd	%%xmm9 ,%%xmm11		\n\t"\
			"subpd	%%xmm4,%%xmm0		\n\t	subpd	%%xmm12,%%xmm8 		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t	subpd	%%xmm13,%%xmm9 		\n\t"\
			"movaps	    (%%rsi),%%xmm4	\n\t	movaps	    (%%r10),%%xmm12	\n\t"\
			"movaps	0x10(%%rsi),%%xmm5	\n\t	movaps	0x10(%%r10),%%xmm13	\n\t"\
			"movq	%[__o1],%%rax		\n\t	movq	%[__u1],%%r11		\n\t"\
			"movq	%[__o4],%%rdx		\n\t	movq	%[__u4],%%r14		\n\t"\
			"subpd	%%xmm3,%%xmm6		\n\t	subpd	%%xmm11,%%xmm14		\n\t"\
			"subpd	%%xmm2,%%xmm7		\n\t	subpd	%%xmm10,%%xmm15		\n\t"\
			"addpd	%%xmm3,%%xmm3		\n\t	addpd	%%xmm11,%%xmm11		\n\t"\
			"addpd	%%xmm2,%%xmm2		\n\t	addpd	%%xmm10,%%xmm10		\n\t"\
			"movaps	%%xmm6,    (%%rax)	\n\t	movaps	%%xmm14,    (%%r11)	\n\t"\
			"movaps	%%xmm7,0x10(%%rdx)	\n\t	movaps	%%xmm15,0x10(%%r14)	\n\t"\
			"addpd	%%xmm6,%%xmm3		\n\t	addpd	%%xmm14,%%xmm11		\n\t"\
			"addpd	%%xmm7,%%xmm2		\n\t	addpd	%%xmm15,%%xmm10		\n\t"\
			"movaps	%%xmm3,    (%%rdx)	\n\t	movaps	%%xmm11,    (%%r14)	\n\t"\
			"movaps	%%xmm2,0x10(%%rax)	\n\t	movaps	%%xmm10,0x10(%%r11)	\n\t"\
			"movq	%[__o2],%%rbx		\n\t	movq	%[__u2],%%r12		\n\t"\
			"movq	%[__o3],%%rcx		\n\t	movq	%[__u3],%%r13		\n\t"\
			"subpd	%%xmm1,%%xmm4		\n\t	subpd	%%xmm9 ,%%xmm12		\n\t"\
			"subpd	%%xmm0,%%xmm5		\n\t	subpd	%%xmm8 ,%%xmm13		\n\t"\
			"addpd	%%xmm1,%%xmm1		\n\t	addpd	%%xmm9 ,%%xmm9 		\n\t"\
			"addpd	%%xmm0,%%xmm0		\n\t	addpd	%%xmm8 ,%%xmm8 		\n\t"\
			"movaps	%%xmm4,    (%%rbx)	\n\t	movaps	%%xmm12,    (%%r12)	\n\t"\
			"movaps	%%xmm5,0x10(%%rcx)	\n\t	movaps	%%xmm13,0x10(%%r13)	\n\t"\
			"addpd	%%xmm4,%%xmm1		\n\t	addpd	%%xmm12,%%xmm9 		\n\t"\
			"addpd	%%xmm5,%%xmm0		\n\t	addpd	%%xmm13,%%xmm8 		\n\t"\
			"movaps	%%xmm1,    (%%rcx)	\n\t	movaps	%%xmm9 ,    (%%r13)	\n\t"\
			"movaps	%%xmm0,0x10(%%rbx)	\n\t	movaps	%%xmm8 ,0x10(%%r12)	\n\t"\
			:					/* outputs: none */\
			: [__cc1] "m" (Xcc1)	/* All inputs from memory addresses here */\
			 ,[__i0] "m" (Xi0)\
			 ,[__i1] "m" (Xi1)\
			 ,[__i2] "m" (Xi2)\
			 ,[__i3] "m" (Xi3)\
			 ,[__i4] "m" (Xi4)\
			 ,[__o0] "m" (Xo0)\
			 ,[__o1] "m" (Xo1)\
			 ,[__o2] "m" (Xo2)\
			 ,[__o3] "m" (Xo3)\
			 ,[__o4] "m" (Xo4)\
			 ,[__j0] "m" (Xj0)\
			 ,[__j1] "m" (Xj1)\
			 ,[__j2] "m" (Xj2)\
			 ,[__j3] "m" (Xj3)\
			 ,[__j4] "m" (Xj4)\
			 ,[__u0] "m" (Xu0)\
			 ,[__u1] "m" (Xu1)\
			 ,[__u2] "m" (Xu2)\
			 ,[__u3] "m" (Xu3)\
			 ,[__u4] "m" (Xu4)\
			: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","r14","r15","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
		);\
		}

	  #endif	// AVX or SSE2?

		#define SSE2_RADIX_15_DIF_X2(\
			__cc0, __cc1,\
			__i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA,__iB,__iC,__iD,__iE,\
			__s0,__s1,__s2,__s3,__s4,__s5,__s6,__s7,__s8,__s9,__sA,__sB,__sC,__sD,__sE,\
			__o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA,__oB,__oC,__oD,__oE,\
					__j0,__j1,__j2,__j3,__j4,__j5,__j6,__j7,__j8,__j9,__jA,__jB,__jC,__jD,__jE,\
					__t0,__t1,__t2,__t3,__t4,__t5,__t6,__t7,__t8,__t9,__tA,__tB,__tC,__tD,__tE,\
					__u0,__u1,__u2,__u3,__u4,__u5,__u6,__u7,__u8,__u9,__uA,__uB,__uC,__uD,__uE)\
		{\
			SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1, __i0,__iC,__i9,__i6,__i3, __s0,__s1,__s2,__s3,__s4,	__j0,__jC,__j9,__j6,__j3, __t0,__t1,__t2,__t3,__t4);\
			SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1, __iA,__i7,__i4,__i1,__iD, __s5,__s6,__s7,__s8,__s9,	__jA,__j7,__j4,__j1,__jD, __t5,__t6,__t7,__t8,__t9);\
			SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1, __i5,__i2,__iE,__iB,__i8, __sA,__sB,__sC,__sD,__sE,	__j5,__j2,__jE,__jB,__j8, __tA,__tB,__tC,__tD,__tE);\
		\
			SSE2_RADIX_03_DFT_X2(__cc0, __s0,__s5,__sA, __o0,__o1,__o2,		__t0,__t5,__tA, __u0,__u1,__u2);\
			SSE2_RADIX_03_DFT_X2(__cc0, __s1,__s6,__sB, __oD,__oE,__oC,		__t1,__t6,__tB, __uD,__uE,__uC);\
			SSE2_RADIX_03_DFT_X2(__cc0, __s2,__s7,__sC, __o9,__oA,__oB,		__t2,__t7,__tC, __u9,__uA,__uB);\
			SSE2_RADIX_03_DFT_X2(__cc0, __s3,__s8,__sD, __o8,__o6,__o7,		__t3,__t8,__tD, __u8,__u6,__u7);\
			SSE2_RADIX_03_DFT_X2(__cc0, __s4,__s9,__sE, __o4,__o5,__o3,		__t4,__t9,__tE, __u4,__u5,__u3);\
		}

		#define SSE2_RADIX_15_DIT_X2(\
			__cc0, __cc1,\
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
			SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1, __s0,__s3,__s6,__s9,__sC, __o0,__o6,__oC,__o3,__o9,	__t0,__t3,__t6,__t9,__tC, __u0,__u6,__uC,__u3,__u9);\
			SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1, __s1,__s4,__s7,__sA,__sD, __o5,__oB,__o2,__o8,__oE,	__t1,__t4,__t7,__tA,__tD, __u5,__uB,__u2,__u8,__uE);\
			SSE2_RADIX_05_DFT_0TWIDDLE_X2(__cc1, __s2,__s5,__s8,__sB,__sE, __oA,__o1,__o7,__oD,__o4,	__t2,__t5,__t8,__tB,__tE, __uA,__u1,__u7,__uD,__u4);\
		}

	#endif	// USE_64BIT_ASM_STYLE ?

	#ifdef COMPILER_TYPE_GCC

		#if OS_BITS == 32

		//	#include "radix60_ditN_cy_dif1_gcc32.h"

		#else

			#include "radix60_ditN_cy_dif1_gcc64.h"

		#endif

	#endif

#endif	/* USE_SSE2 */

#ifdef USE_PTHREAD

	// Use non-pooled simple spawn/rejoin thread-team model
	#include "threadpool.h"

	struct cy_thread_data_t{
	// int data - if needed, pad to yield an even number of these:
		int tid;
		int ndivr;
		uint32 nrtm1;
		uint32 nrt_bits;

		int khi;
		int i;
		int jstart;
		int jhi;
		int col;
		int co2;
		int co3;
		int sw;
		int nwt;
		int wts_idx_inc2;
		int _pad;
		int icycle00;
		int icycle01;
		int icycle02;
		int icycle03;
		int icycle04;
		int icycle05;
		int icycle06;
		int icycle07;
		int icycle08;
		int icycle09;
		int icycle10;
		int icycle11;
		int icycle12;
		int icycle13;
		int icycle14;

	#ifdef USE_SSE2
		int jcycle00;
		int jcycle01;
		int jcycle02;
		int jcycle03;
		int jcycle04;
		int jcycle05;
		int jcycle06;
		int jcycle07;
		int jcycle08;
		int jcycle09;
		int jcycle10;
		int jcycle11;
		int jcycle12;
		int jcycle13;
		int jcycle14;

	  #ifdef USE_AVX
		int kcycle00;
		int kcycle01;
		int kcycle02;
		int kcycle03;
		int kcycle04;
		int kcycle05;
		int kcycle06;
		int kcycle07;
		int kcycle08;
		int kcycle09;
		int kcycle10;
		int kcycle11;
		int kcycle12;
		int kcycle13;
		int kcycle14;

		int lcycle00;
		int lcycle01;
		int lcycle02;
		int lcycle03;
		int lcycle04;
		int lcycle05;
		int lcycle06;
		int lcycle07;
		int lcycle08;
		int lcycle09;
		int lcycle10;
		int lcycle11;
		int lcycle12;
		int lcycle13;
		int lcycle14;
	  #endif
	#endif

	// double data:
		double maxerr;
		double scale;

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
		int *si;
		struct complex *rn0;
		struct complex *rn1;
		vec_dbl *s1p00r;
		vec_dbl *r00;
		vec_dbl *half_arr;
		uint64*sm_ptr;

		int bjmodn00;
		int bjmodn01;
		int bjmodn02;
		int bjmodn03;
		int bjmodn04;
		int bjmodn05;
		int bjmodn06;
		int bjmodn07;
		int bjmodn08;
		int bjmodn09;
		int bjmodn10;
		int bjmodn11;
		int bjmodn12;
		int bjmodn13;
		int bjmodn14;
		int bjmodn15;
		int bjmodn16;
		int bjmodn17;
		int bjmodn18;
		int bjmodn19;
		int bjmodn20;
		int bjmodn21;
		int bjmodn22;
		int bjmodn23;
		int bjmodn24;
		int bjmodn25;
		int bjmodn26;
		int bjmodn27;
		int bjmodn28;
		int bjmodn29;
		int bjmodn30;
		int bjmodn31;
		int bjmodn32;
		int bjmodn33;
		int bjmodn34;
		int bjmodn35;
		int bjmodn36;
		int bjmodn37;
		int bjmodn38;
		int bjmodn39;
		int bjmodn40;
		int bjmodn41;
		int bjmodn42;
		int bjmodn43;
		int bjmodn44;
		int bjmodn45;
		int bjmodn46;
		int bjmodn47;
		int bjmodn48;
		int bjmodn49;
		int bjmodn50;
		int bjmodn51;
		int bjmodn52;
		int bjmodn53;
		int bjmodn54;
		int bjmodn55;
		int bjmodn56;
		int bjmodn57;
		int bjmodn58;
		int bjmodn59;
		/* carries: */
		double cy_r00;
		double cy_r01;
		double cy_r02;
		double cy_r03;
		double cy_r04;
		double cy_r05;
		double cy_r06;
		double cy_r07;
		double cy_r08;
		double cy_r09;
		double cy_r10;
		double cy_r11;
		double cy_r12;
		double cy_r13;
		double cy_r14;
		double cy_r15;
		double cy_r16;
		double cy_r17;
		double cy_r18;
		double cy_r19;
		double cy_r20;
		double cy_r21;
		double cy_r22;
		double cy_r23;
		double cy_r24;
		double cy_r25;
		double cy_r26;
		double cy_r27;
		double cy_r28;
		double cy_r29;
		double cy_r30;
		double cy_r31;
		double cy_r32;
		double cy_r33;
		double cy_r34;
		double cy_r35;
		double cy_r36;
		double cy_r37;
		double cy_r38;
		double cy_r39;
		double cy_r40;
		double cy_r41;
		double cy_r42;
		double cy_r43;
		double cy_r44;
		double cy_r45;
		double cy_r46;
		double cy_r47;
		double cy_r48;
		double cy_r49;
		double cy_r50;
		double cy_r51;
		double cy_r52;
		double cy_r53;
		double cy_r54;
		double cy_r55;
		double cy_r56;
		double cy_r57;
		double cy_r58;
		double cy_r59;

		double cy_i00;
		double cy_i01;
		double cy_i02;
		double cy_i03;
		double cy_i04;
		double cy_i05;
		double cy_i06;
		double cy_i07;
		double cy_i08;
		double cy_i09;
		double cy_i10;
		double cy_i11;
		double cy_i12;
		double cy_i13;
		double cy_i14;
		double cy_i15;
		double cy_i16;
		double cy_i17;
		double cy_i18;
		double cy_i19;
		double cy_i20;
		double cy_i21;
		double cy_i22;
		double cy_i23;
		double cy_i24;
		double cy_i25;
		double cy_i26;
		double cy_i27;
		double cy_i28;
		double cy_i29;
		double cy_i30;
		double cy_i31;
		double cy_i32;
		double cy_i33;
		double cy_i34;
		double cy_i35;
		double cy_i36;
		double cy_i37;
		double cy_i38;
		double cy_i39;
		double cy_i40;
		double cy_i41;
		double cy_i42;
		double cy_i43;
		double cy_i44;
		double cy_i45;
		double cy_i46;
		double cy_i47;
		double cy_i48;
		double cy_i49;
		double cy_i50;
		double cy_i51;
		double cy_i52;
		double cy_i53;
		double cy_i54;
		double cy_i55;
		double cy_i56;
		double cy_i57;
		double cy_i58;
		double cy_i59;
	};

#endif

/**************/

int radix60_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-[radix] complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-[radix] complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const char func[] = "radix60_ditN_cy_dif1";
	const uint32 RADIX = 60, odd_radix = 15;	// odd_radix = [radix >> trailz(radix)]
	const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
#ifdef USE_SSE2
	const int sz_vd = sizeof(vec_dbl);
	// lg(sizeof(vec_dbl)):
  #ifdef USE_AVX
	const int l2_sz_vd = 5;
  #else
	const int l2_sz_vd = 4;
  #endif
#endif
	int NDIVR,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,ntmp,outer,nbytes;
#ifdef USE_SSE2
	uint32 nwt16 = nwt << l2_sz_vd;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code
#endif
	// Need these both in scalar mode and to ease the SSE2-array init...dimension = odd_radix;
	// In order to ease the ptr-access for the || routine, lump these 4*odd_radix doubles together with copies of
	// the 4 in the passed-in bs[2] and bsinv[2] arrays [and used in this 4-double form by the mersenne-mod carry macros]
	// into a single foo_array[4*(odd_radix+1)], then convert what used to be disparate odd_radix-sized arrays to pointers.
	static double foo_array[(odd_radix+1)<<2], *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr, bs,bsinv;
	// NOTE: If compiler squawks, Use literal value of odd_radix in dimensioning in order to allow us to explicitly make static:

	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40,p44,p48,p52,p56;
	static double radix_inv, n2inv;
	const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/
					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double scale,
	t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,
	t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,
	t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,
	t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei;
	double dtmp, maxerr = 0.0;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	static int n_div_nwt;
	int col,co2,co3,k1,k2,m,m2;
  #ifdef USE_AVX
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
  #endif
	double rt,it,wt_re,wt_im;	/* Fermat-mod weights stuff, used in both scalar and AVX mode */
	int ii00,ii01,ii02,ii03,ii04,ii05,ii06,ii07,ii08,ii09,ii10,ii11,ii12,ii13,ii14;	/* indices into weights arrays (mod NWT) */
	/* These are used in conjunction with the langth-7 arrays in the USE_SCALAR_CARRY code flow;
	In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
	static int idx_offset, idx_incr, wts_idx_incr = 0, wts_idx_inc2 = 0
		,icycle00,icycle01,icycle02,icycle03,icycle04,icycle05,icycle06,icycle07,icycle08,icycle09,icycle10,icycle11,icycle12,icycle13,icycle14;

#ifdef USE_SSE2

	static int jcycle00,jcycle01,jcycle02,jcycle03,jcycle04,jcycle05,jcycle06,jcycle07,jcycle08,jcycle09,jcycle10,jcycle11,jcycle12,jcycle13,jcycle14;
  #ifdef USE_AVX
	static int kcycle00,kcycle01,kcycle02,kcycle03,kcycle04,kcycle05,kcycle06,kcycle07,kcycle08,kcycle09,kcycle10,kcycle11,kcycle12,kcycle13,kcycle14;
	static int lcycle00,lcycle01,lcycle02,lcycle03,lcycle04,lcycle05,lcycle06,lcycle07,lcycle08,lcycle09,lcycle10,lcycle11,lcycle12,lcycle13,lcycle14;
  #endif

	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_arr = 0x0, *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else
	double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
  #endif	// MULTITHREAD

	const double crnd = 3.0*0x4000000*0x2000000;
	static int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,
		*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,
		*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27,*bjmodn28,*bjmodn29,
		*bjmodn30,*bjmodn31,*bjmodn32,*bjmodn33,*bjmodn34,*bjmodn35,*bjmodn36,*bjmodn37,*bjmodn38,*bjmodn39,
		*bjmodn40,*bjmodn41,*bjmodn42,*bjmodn43,*bjmodn44,*bjmodn45,*bjmodn46,*bjmodn47,*bjmodn48,*bjmodn49,
		*bjmodn50,*bjmodn51,*bjmodn52,*bjmodn53,*bjmodn54,*bjmodn55,*bjmodn56,*bjmodn57,*bjmodn58,*bjmodn59;
	static vec_dbl *max_err, *sse2_rnd, *half_arr, *sse2_c3m1, *sse2_s, *sse2_cn1, *sse2_cn2, *sse2_ss3, *sse2_sn1, *sse2_sn2,
		*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p0ar,*s1p0br,*s1p0cr,*s1p0dr,*s1p0er,
		*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p1ar,*s1p1br,*s1p1cr,*s1p1dr,*s1p1er,
		*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r,*s1p28r,*s1p29r,*s1p2ar,*s1p2br,*s1p2cr,*s1p2dr,*s1p2er,
		*s1p30r,*s1p31r,*s1p32r,*s1p33r,*s1p34r,*s1p35r,*s1p36r,*s1p37r,*s1p38r,*s1p39r,*s1p3ar,*s1p3br,*s1p3cr,*s1p3dr,*s1p3er,
		*r00,*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09,*r0a,*r0b,*r0c,*r0d,*r0e,
		*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r1a,*r1b,*r1c,*r1d,*r1e,
		*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r2a,*r2b,*r2c,*r2d,*r2e,
		*r30,*r31,*r32,*r33,*r34,*r35,*r36,*r37,*r38,*r39,*r3a,*r3b,*r3c,*r3d,*r3e,
		*x00,*x01,*x02,*x03,*x04,*x05,*x06,*x07,*x08,*x09,*x0a,*x0b,*x0c,*x0d,*x0e,
		*y00,*y01,*y02,*y03,*y04,*y05,*y06,*y07,*y08,*y09,*y0a,*y0b,*y0c,*y0d,*y0e,
		*cy_r00,*cy_r04,*cy_r08,*cy_r12,*cy_r16,*cy_r20,*cy_r24,*cy_r28,*cy_r32,*cy_r36,*cy_r40,*cy_r44,*cy_r48,*cy_r52,*cy_r56,
		*cy_i00,*cy_i04,*cy_i08,*cy_i12,*cy_i16,*cy_i20,*cy_i24,*cy_i28,*cy_i32,*cy_i36,*cy_i40,*cy_i44,*cy_i48,*cy_i52,*cy_i56;
  #ifndef USE_AVX
	static vec_dbl
		*cy_r02,*cy_r06,*cy_r10,*cy_r14,*cy_r18,*cy_r22,*cy_r26,*cy_r30,*cy_r34,*cy_r38,*cy_r42,*cy_r46,*cy_r50,*cy_r54,*cy_r58,
		*cy_i02,*cy_i06,*cy_i10,*cy_i14,*cy_i18,*cy_i22,*cy_i26,*cy_i30,*cy_i34,*cy_i38,*cy_i42,*cy_i46,*cy_i50,*cy_i54,*cy_i58;
  #else
	static vec_dbl *base_negacyclic_root;
  #endif

	vec_dbl *tmp,*tm2;	// Non-static utility ptrs
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer

#endif	// USE_SSE2

#ifdef MULTITHREAD

	static struct cy_thread_data_t *tdat = 0x0;
	// Threadpool-based dispatch stuff:
	static int main_work_units = 0, pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)cy60_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
  #if PFETCH
	double *addr, *addp;
  #endif
	int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,
		bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,
		bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,bjmodn28,bjmodn29,
		bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,bjmodn36,bjmodn37,bjmodn38,bjmodn39,
		bjmodn40,bjmodn41,bjmodn42,bjmodn43,bjmodn44,bjmodn45,bjmodn46,bjmodn47,bjmodn48,bjmodn49,
		bjmodn50,bjmodn51,bjmodn52,bjmodn53,bjmodn54,bjmodn55,bjmodn56,bjmodn57,bjmodn58,bjmodn59;
	double re,im,temp,frac,
		a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,a1p08r,a1p08i,a1p09r,a1p09i,a1p0ar,a1p0ai,a1p0br,a1p0bi,a1p0cr,a1p0ci,a1p0dr,a1p0di,a1p0er,a1p0ei,
		a1p10r,a1p10i,a1p11r,a1p11i,a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a1p1ar,a1p1ai,a1p1br,a1p1bi,a1p1cr,a1p1ci,a1p1dr,a1p1di,a1p1er,a1p1ei,
		a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,a1p24r,a1p24i,a1p25r,a1p25i,a1p26r,a1p26i,a1p27r,a1p27i,a1p28r,a1p28i,a1p29r,a1p29i,a1p2ar,a1p2ai,a1p2br,a1p2bi,a1p2cr,a1p2ci,a1p2dr,a1p2di,a1p2er,a1p2ei,
		a1p30r,a1p30i,a1p31r,a1p31i,a1p32r,a1p32i,a1p33r,a1p33i,a1p34r,a1p34i,a1p35r,a1p35i,a1p36r,a1p36i,a1p37r,a1p37i,a1p38r,a1p38i,a1p39r,a1p39i,a1p3ar,a1p3ai,a1p3br,a1p3bi,a1p3cr,a1p3ci,a1p3dr,a1p3di,a1p3er,a1p3ei,
		cy_r00,cy_i00,cy_r01,cy_i01,cy_r02,cy_i02,cy_r03,cy_i03,cy_r04,cy_i04,cy_r05,cy_i05,cy_r06,cy_i06,cy_r07,cy_i07,cy_r08,cy_i08,cy_r09,cy_i09,
		cy_r10,cy_i10,cy_r11,cy_i11,cy_r12,cy_i12,cy_r13,cy_i13,cy_r14,cy_i14,cy_r15,cy_i15,cy_r16,cy_i16,cy_r17,cy_i17,cy_r18,cy_i18,cy_r19,cy_i19,
		cy_r20,cy_i20,cy_r21,cy_i21,cy_r22,cy_i22,cy_r23,cy_i23,cy_r24,cy_i24,cy_r25,cy_i25,cy_r26,cy_i26,cy_r27,cy_i27,cy_r28,cy_i28,cy_r29,cy_i29,
		cy_r30,cy_i30,cy_r31,cy_i31,cy_r32,cy_i32,cy_r33,cy_i33,cy_r34,cy_i34,cy_r35,cy_i35,cy_r36,cy_i36,cy_r37,cy_i37,cy_r38,cy_i38,cy_r39,cy_i39,
		cy_r40,cy_i40,cy_r41,cy_i41,cy_r42,cy_i42,cy_r43,cy_i43,cy_r44,cy_i44,cy_r45,cy_i45,cy_r46,cy_i46,cy_r47,cy_i47,cy_r48,cy_i48,cy_r49,cy_i49,
		cy_r50,cy_i50,cy_r51,cy_i51,cy_r52,cy_i52,cy_r53,cy_i53,cy_r54,cy_i54,cy_r55,cy_i55,cy_r56,cy_i56,cy_r57,cy_i57,cy_r58,cy_i58,cy_r59,cy_i59,
		x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static int *_bjmodnini = 0x0,
		*_bjmodn00 = 0x0,*_bjmodn01 = 0x0,*_bjmodn02 = 0x0,*_bjmodn03 = 0x0,*_bjmodn04 = 0x0,*_bjmodn05 = 0x0,*_bjmodn06 = 0x0,*_bjmodn07 = 0x0,*_bjmodn08 = 0x0,*_bjmodn09 = 0x0,
		*_bjmodn10 = 0x0,*_bjmodn11 = 0x0,*_bjmodn12 = 0x0,*_bjmodn13 = 0x0,*_bjmodn14 = 0x0,*_bjmodn15 = 0x0,*_bjmodn16 = 0x0,*_bjmodn17 = 0x0,*_bjmodn18 = 0x0,*_bjmodn19 = 0x0,
		*_bjmodn20 = 0x0,*_bjmodn21 = 0x0,*_bjmodn22 = 0x0,*_bjmodn23 = 0x0,*_bjmodn24 = 0x0,*_bjmodn25 = 0x0,*_bjmodn26 = 0x0,*_bjmodn27 = 0x0,*_bjmodn28 = 0x0,*_bjmodn29 = 0x0,
		*_bjmodn30 = 0x0,*_bjmodn31 = 0x0,*_bjmodn32 = 0x0,*_bjmodn33 = 0x0,*_bjmodn34 = 0x0,*_bjmodn35 = 0x0,*_bjmodn36 = 0x0,*_bjmodn37 = 0x0,*_bjmodn38 = 0x0,*_bjmodn39 = 0x0,
		*_bjmodn40 = 0x0,*_bjmodn41 = 0x0,*_bjmodn42 = 0x0,*_bjmodn43 = 0x0,*_bjmodn44 = 0x0,*_bjmodn45 = 0x0,*_bjmodn46 = 0x0,*_bjmodn47 = 0x0,*_bjmodn48 = 0x0,*_bjmodn49 = 0x0,
		*_bjmodn50 = 0x0,*_bjmodn51 = 0x0,*_bjmodn52 = 0x0,*_bjmodn53 = 0x0,*_bjmodn54 = 0x0,*_bjmodn55 = 0x0,*_bjmodn56 = 0x0,*_bjmodn57 = 0x0,*_bjmodn58 = 0x0,*_bjmodn59 = 0x0;
	static double *_maxerr = 0x0,
		*_cy_r00 = 0x0,*_cy_i00 = 0x0,*_cy_r01 = 0x0,*_cy_i01 = 0x0,*_cy_r02 = 0x0,*_cy_i02 = 0x0,*_cy_r03 = 0x0,*_cy_i03 = 0x0,*_cy_r04 = 0x0,*_cy_i04 = 0x0,*_cy_r05 = 0x0,*_cy_i05 = 0x0,*_cy_r06 = 0x0,*_cy_i06 = 0x0,*_cy_r07 = 0x0,*_cy_i07 = 0x0,*_cy_r08 = 0x0,*_cy_i08 = 0x0,*_cy_r09 = 0x0,*_cy_i09 = 0x0,
		*_cy_r10 = 0x0,*_cy_i10 = 0x0,*_cy_r11 = 0x0,*_cy_i11 = 0x0,*_cy_r12 = 0x0,*_cy_i12 = 0x0,*_cy_r13 = 0x0,*_cy_i13 = 0x0,*_cy_r14 = 0x0,*_cy_i14 = 0x0,*_cy_r15 = 0x0,*_cy_i15 = 0x0,*_cy_r16 = 0x0,*_cy_i16 = 0x0,*_cy_r17 = 0x0,*_cy_i17 = 0x0,*_cy_r18 = 0x0,*_cy_i18 = 0x0,*_cy_r19 = 0x0,*_cy_i19 = 0x0,
		*_cy_r20 = 0x0,*_cy_i20 = 0x0,*_cy_r21 = 0x0,*_cy_i21 = 0x0,*_cy_r22 = 0x0,*_cy_i22 = 0x0,*_cy_r23 = 0x0,*_cy_i23 = 0x0,*_cy_r24 = 0x0,*_cy_i24 = 0x0,*_cy_r25 = 0x0,*_cy_i25 = 0x0,*_cy_r26 = 0x0,*_cy_i26 = 0x0,*_cy_r27 = 0x0,*_cy_i27 = 0x0,*_cy_r28 = 0x0,*_cy_i28 = 0x0,*_cy_r29 = 0x0,*_cy_i29 = 0x0,
		*_cy_r30 = 0x0,*_cy_i30 = 0x0,*_cy_r31 = 0x0,*_cy_i31 = 0x0,*_cy_r32 = 0x0,*_cy_i32 = 0x0,*_cy_r33 = 0x0,*_cy_i33 = 0x0,*_cy_r34 = 0x0,*_cy_i34 = 0x0,*_cy_r35 = 0x0,*_cy_i35 = 0x0,*_cy_r36 = 0x0,*_cy_i36 = 0x0,*_cy_r37 = 0x0,*_cy_i37 = 0x0,*_cy_r38 = 0x0,*_cy_i38 = 0x0,*_cy_r39 = 0x0,*_cy_i39 = 0x0,
		*_cy_r40 = 0x0,*_cy_i40 = 0x0,*_cy_r41 = 0x0,*_cy_i41 = 0x0,*_cy_r42 = 0x0,*_cy_i42 = 0x0,*_cy_r43 = 0x0,*_cy_i43 = 0x0,*_cy_r44 = 0x0,*_cy_i44 = 0x0,*_cy_r45 = 0x0,*_cy_i45 = 0x0,*_cy_r46 = 0x0,*_cy_i46 = 0x0,*_cy_r47 = 0x0,*_cy_i47 = 0x0,*_cy_r48 = 0x0,*_cy_i48 = 0x0,*_cy_r49 = 0x0,*_cy_i49 = 0x0,
		*_cy_r50 = 0x0,*_cy_i50 = 0x0,*_cy_r51 = 0x0,*_cy_i51 = 0x0,*_cy_r52 = 0x0,*_cy_i52 = 0x0,*_cy_r53 = 0x0,*_cy_i53 = 0x0,*_cy_r54 = 0x0,*_cy_i54 = 0x0,*_cy_r55 = 0x0,*_cy_i55 = 0x0,*_cy_r56 = 0x0,*_cy_i56 = 0x0,*_cy_r57 = 0x0,*_cy_i57 = 0x0,*_cy_r58 = 0x0,*_cy_i58 = 0x0,*_cy_r59 = 0x0,*_cy_i59 = 0x0;

	foo_array[0] = base[0];
	foo_array[1] = base[1];
	foo_array[2] = baseinv[0];
	foo_array[3] = baseinv[1];
	wt_arr    = foo_array + 4;
	wtinv_arr = wt_arr    + odd_radix;
	bs_arr    = wtinv_arr + odd_radix;
	bsinv_arr = bs_arr    + odd_radix;

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=ii00=ii01=ii02=ii03=ii04=ii05=ii06=ii07=ii08=ii09=ii10=ii11=ii12=ii13=ii14=-1;

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/RADIX;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/%d in %s.\n", iter,RADIX,func);
		if(INTERACT)fprintf(stderr,"%s",cbuf);
		fp = fopen(   OFILE,"a");
		fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;
		err=ERR_CARRY;
		return(err);
	}

	if(p != psave)
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)RADIX));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;	/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

	#ifdef MULTITHREAD

		/* #Chunks ||ized in carry step is ideally a power of 2, so use the smallest
		power of 2 that is >= the value of the global NTHREADS (but still <= MAX_THREADS):
		*/
		if(isPow2(NTHREADS))
			CY_THREADS = NTHREADS;
		else
		{
			i = leadz32(NTHREADS);
			CY_THREADS = (((uint32)NTHREADS << i) & 0x80000000) >> (i-1);
		}

		if(CY_THREADS > MAX_THREADS)
		{
		//	CY_THREADS = MAX_THREADS;
			fprintf(stderr,"WARN: CY_THREADS = %d exceeds number of cores = %d\n", CY_THREADS, MAX_THREADS);
		}
		ASSERT(HERE, CY_THREADS >= NTHREADS,"CY_THREADS < NTHREADS");
		ASSERT(HERE, isPow2(CY_THREADS)    ,"CY_THREADS not a power of 2!");
		if(CY_THREADS > 1)
		{
			ASSERT(HERE, NDIVR    %CY_THREADS == 0,"NDIVR    %CY_THREADS != 0");
			ASSERT(HERE, n_div_nwt%CY_THREADS == 0,"n_div_nwt%CY_THREADS != 0");
		}

	  #ifdef USE_PTHREAD

		j = (uint32)sizeof(struct cy_thread_data_t);
		tdat = (struct cy_thread_data_t *)calloc(CY_THREADS, j);

		// MacOS does weird things with threading (e.g. Idle" main thread burning 100% of 1 CPU)
		// so on that platform try to be clever and interleave main-thread and threadpool-work processing
		#ifdef OS_TYPE_MACOSX

			if(CY_THREADS > 1) {
				main_work_units = CY_THREADS/2;
				pool_work_units = CY_THREADS - main_work_units;
				ASSERT(HERE, 0x0 != (tpool = threadpool_init(pool_work_units, MAX_THREADS, pool_work_units, &thread_control)), "threadpool_init failed!");
				printf("radix%d_ditN_cy_dif1: Init threadpool of %d threads\n", RADIX, pool_work_units);
			} else {
				main_work_units = 1;
				printf("radix%d_ditN_cy_dif1: CY_THREADS = 1: Using main execution thread, no threadpool needed.\n", RADIX);
			}

		#else

			pool_work_units = CY_THREADS;
			ASSERT(HERE, 0x0 != (tpool = threadpool_init(CY_THREADS, MAX_THREADS, CY_THREADS, &thread_control)), "threadpool_init failed!");

		#endif

		fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);

	  #endif

	#else
		CY_THREADS = 1;
	#endif

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
		// int data:
			tdat[ithread].tid = ithread;
			tdat[ithread].ndivr = NDIVR;

			tdat[ithread].sw  = sw;
			tdat[ithread].nwt = nwt;

		// pointer data:
			tdat[ithread].arrdat = a;			/* Main data array */
			tdat[ithread].wt0 = wt0;
			tdat[ithread].wt1 = wt1;
			tdat[ithread].si  = si;
			tdat[ithread].rn0 = rn0;
			tdat[ithread].rn1 = rn1;
		}
	#endif

	#ifdef USE_SSE2

		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use vector-double type size (16 bytes for SSE2, 32 for AVX) to alloc a block of local storage
		// consisting of 128*2 vec_dbl and (8+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
	#ifdef USE_PTHREAD
		cslots_in_local_store = radix60_creals_in_local_store + (12+RADIX);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

		sm_ptr = (uint64*)(sc_ptr + radix60_creals_in_local_store);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	#elif 0	/*** for some reason the same kind of alloc code as used in radix28 fails here - get memory corruption when init bjmod ptrs ***/

		cslots_in_local_store = radix60_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix60_creals_in_local_store);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	#else

		sc_arr = ALLOC_VEC_DBL(sc_arr, radix60_creals_in_local_store);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

		/* Size here is (8 + radix/2 + 4) 8-byte elements */
		sm_arr = ALLOC_UINT64(sm_arr, (8 + RADIX/2 + 4));	if(!sm_arr){ sprintf(cbuf, "FATAL: unable to allocate sm_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sm_ptr = ALIGN_UINT64(sm_arr);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");
	#endif

	/* Use low 2*RADIX vector-double-sized slots of sc_arr for s1p* temporaries, next 2*RADIX slots for r* temps,
	next RADIX slots for x and y in-place DFT temps, next 7 for the complex root combs needed for the radix-3 and -5 sub-DFTs,
	next 30[avx] | 60[sse2] for the doubled carry pairs, next 2 for ROE and RND_CONST, next 5*RADIX[avx] | 20[sse2] for the
	half_arr table lookup stuff, plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	  #ifdef USE_PTHREAD
		__r0 = sc_ptr;
	  #endif						// Use these as handy address-offset temps:
									tm2 = sc_ptr + 0x1e;    max_err  = tm2 + 0x1e;    sse2_rnd = max_err + 0x1e;
		s1p00r = sc_ptr + 0x00;    s1p10r = tm2 + 0x00;    s1p20r = max_err + 0x00;    s1p30r = sse2_rnd + 0x00;
		s1p01r = sc_ptr + 0x02;    s1p11r = tm2 + 0x02;    s1p21r = max_err + 0x02;    s1p31r = sse2_rnd + 0x02;
		s1p02r = sc_ptr + 0x04;    s1p12r = tm2 + 0x04;    s1p22r = max_err + 0x04;    s1p32r = sse2_rnd + 0x04;
		s1p03r = sc_ptr + 0x06;    s1p13r = tm2 + 0x06;    s1p23r = max_err + 0x06;    s1p33r = sse2_rnd + 0x06;
		s1p04r = sc_ptr + 0x08;    s1p14r = tm2 + 0x08;    s1p24r = max_err + 0x08;    s1p34r = sse2_rnd + 0x08;
		s1p05r = sc_ptr + 0x0a;    s1p15r = tm2 + 0x0a;    s1p25r = max_err + 0x0a;    s1p35r = sse2_rnd + 0x0a;
		s1p06r = sc_ptr + 0x0c;    s1p16r = tm2 + 0x0c;    s1p26r = max_err + 0x0c;    s1p36r = sse2_rnd + 0x0c;
		s1p07r = sc_ptr + 0x0e;    s1p17r = tm2 + 0x0e;    s1p27r = max_err + 0x0e;    s1p37r = sse2_rnd + 0x0e;
		s1p08r = sc_ptr + 0x10;    s1p18r = tm2 + 0x10;    s1p28r = max_err + 0x10;    s1p38r = sse2_rnd + 0x10;
		s1p09r = sc_ptr + 0x12;    s1p19r = tm2 + 0x12;    s1p29r = max_err + 0x12;    s1p39r = sse2_rnd + 0x12;
		s1p0ar = sc_ptr + 0x14;    s1p1ar = tm2 + 0x14;    s1p2ar = max_err + 0x14;    s1p3ar = sse2_rnd + 0x14;
		s1p0br = sc_ptr + 0x16;    s1p1br = tm2 + 0x16;    s1p2br = max_err + 0x16;    s1p3br = sse2_rnd + 0x16;
		s1p0cr = sc_ptr + 0x18;    s1p1cr = tm2 + 0x18;    s1p2cr = max_err + 0x18;    s1p3cr = sse2_rnd + 0x18;
		s1p0dr = sc_ptr + 0x1a;    s1p1dr = tm2 + 0x1a;    s1p2dr = max_err + 0x1a;    s1p3dr = sse2_rnd + 0x1a;
		s1p0er = sc_ptr + 0x1c;    s1p1er = tm2 + 0x1c;    s1p2er = max_err + 0x1c;    s1p3er = sse2_rnd + 0x1c;	// 120 complex
		tmp = sse2_rnd + 0x1e;  tm2 =    tmp + 0x1e;    max_err  = tm2 + 0x1e;    sse2_rnd = max_err + 0x1e;
		r00    = tmp + 0x00;    r10    = tm2 + 0x00;    r20    = max_err + 0x00;    r30    = sse2_rnd + 0x00;
		r01    = tmp + 0x02;    r11    = tm2 + 0x02;    r21    = max_err + 0x02;    r31    = sse2_rnd + 0x02;
		r02    = tmp + 0x04;    r12    = tm2 + 0x04;    r22    = max_err + 0x04;    r32    = sse2_rnd + 0x04;
		r03    = tmp + 0x06;    r13    = tm2 + 0x06;    r23    = max_err + 0x06;    r33    = sse2_rnd + 0x06;
		r04    = tmp + 0x08;    r14    = tm2 + 0x08;    r24    = max_err + 0x08;    r34    = sse2_rnd + 0x08;
		r05    = tmp + 0x0a;    r15    = tm2 + 0x0a;    r25    = max_err + 0x0a;    r35    = sse2_rnd + 0x0a;
		r06    = tmp + 0x0c;    r16    = tm2 + 0x0c;    r26    = max_err + 0x0c;    r36    = sse2_rnd + 0x0c;
		r07    = tmp + 0x0e;    r17    = tm2 + 0x0e;    r27    = max_err + 0x0e;    r37    = sse2_rnd + 0x0e;
		r08    = tmp + 0x10;    r18    = tm2 + 0x10;    r28    = max_err + 0x10;    r38    = sse2_rnd + 0x10;
		r09    = tmp + 0x12;    r19    = tm2 + 0x12;    r29    = max_err + 0x12;    r39    = sse2_rnd + 0x12;
		r0a    = tmp + 0x14;    r1a    = tm2 + 0x14;    r2a    = max_err + 0x14;    r3a    = sse2_rnd + 0x14;
		r0b    = tmp + 0x16;    r1b    = tm2 + 0x16;    r2b    = max_err + 0x16;    r3b    = sse2_rnd + 0x16;
		r0c    = tmp + 0x18;    r1c    = tm2 + 0x18;    r2c    = max_err + 0x18;    r3c    = sse2_rnd + 0x18;
		r0d    = tmp + 0x1a;    r1d    = tm2 + 0x1a;    r2d    = max_err + 0x1a;    r3d    = sse2_rnd + 0x1a;
		r0e    = tmp + 0x1c;    r1e    = tm2 + 0x1c;    r2e    = max_err + 0x1c;    r3e    = sse2_rnd + 0x1c;	// +120 = 240 complex
		tmp = sse2_rnd + 0x1e;
		x00    = tmp + 0x00;
		x01    = tmp + 0x02;
		x02    = tmp + 0x04;
		x03    = tmp + 0x06;
		x04    = tmp + 0x08;
		x05    = tmp + 0x0a;
		x06    = tmp + 0x0c;
		x07    = tmp + 0x0e;
		x08    = tmp + 0x10;
		x09    = tmp + 0x12;
		x0a    = tmp + 0x14;
		x0b    = tmp + 0x16;
		x0c    = tmp + 0x18;
		x0d    = tmp + 0x1a;
		x0e    = tmp + 0x1c;	// +30 = 270 complex
		tmp += 0x1e;
		y00    = tmp + 0x00;
		y01    = tmp + 0x02;
		y02    = tmp + 0x04;
		y03    = tmp + 0x06;
		y04    = tmp + 0x08;
		y05    = tmp + 0x0a;
		y06    = tmp + 0x0c;
		y07    = tmp + 0x0e;
		y08    = tmp + 0x10;
		y09    = tmp + 0x12;
		y0a    = tmp + 0x14;
		y0b    = tmp + 0x16;
		y0c    = tmp + 0x18;
		y0d    = tmp + 0x1a;
		y0e    = tmp + 0x1c;	// +30 = 300 complex
		tmp += 0x1e;
		sse2_c3m1 = tmp + 0x00;
		sse2_s    = tmp + 0x01;
		sse2_cn1  = tmp + 0x02;
		sse2_cn2  = tmp + 0x03;
		sse2_ss3  = tmp + 0x04;
		sse2_sn1  = tmp + 0x05;
		sse2_sn2  = tmp + 0x06;	// +8 = 308 complex
		tmp += 0x08;
	#ifdef USE_AVX
		cy_r00 = tmp + 0x00;
		cy_r04 = tmp + 0x01;
		cy_r08 = tmp + 0x02;
		cy_r12 = tmp + 0x03;
		cy_r16 = tmp + 0x04;
		cy_r20 = tmp + 0x05;
		cy_r24 = tmp + 0x06;
		cy_r28 = tmp + 0x07;
		cy_r32 = tmp + 0x08;
		cy_r36 = tmp + 0x09;
		cy_r40 = tmp + 0x0a;
		cy_r44 = tmp + 0x0b;
		cy_r48 = tmp + 0x0c;
		cy_r52 = tmp + 0x0d;
		cy_r56 = tmp + 0x0e;
		cy_i00 = tmp + 0x0f;
		cy_i04 = tmp + 0x10;
		cy_i08 = tmp + 0x11;
		cy_i12 = tmp + 0x12;
		cy_i16 = tmp + 0x13;
		cy_i20 = tmp + 0x14;
		cy_i24 = tmp + 0x15;
		cy_i28 = tmp + 0x16;
		cy_i32 = tmp + 0x17;
		cy_i36 = tmp + 0x18;
		cy_i40 = tmp + 0x19;
		cy_i44 = tmp + 0x1a;
		cy_i48 = tmp + 0x1b;
		cy_i52 = tmp + 0x1c;
		cy_i56 = tmp + 0x1d;	// +30 = 338 complex
		tmp += 0x1e;
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 340 complex
		// This is where the value of half_arr_offset60 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
								// +20 = 390 complex, round up to nearby multiple of 4
	#else
		cy_r00 = tmp + 0x00;	cy_r02 = tmp + 0x01;
		cy_r04 = tmp + 0x02;	cy_r06 = tmp + 0x03;
		cy_r08 = tmp + 0x04;	cy_r10 = tmp + 0x05;
		cy_r12 = tmp + 0x06;	cy_r14 = tmp + 0x07;
		cy_r16 = tmp + 0x08;	cy_r18 = tmp + 0x09;
		cy_r20 = tmp + 0x0a;	cy_r22 = tmp + 0x0b;
		cy_r24 = tmp + 0x0c;	cy_r26 = tmp + 0x0d;
		cy_r28 = tmp + 0x0e;	cy_r30 = tmp + 0x0f;
		cy_r32 = tmp + 0x10;	cy_r34 = tmp + 0x11;
		cy_r36 = tmp + 0x12;	cy_r38 = tmp + 0x13;
		cy_r40 = tmp + 0x14;	cy_r42 = tmp + 0x15;
		cy_r44 = tmp + 0x16;	cy_r46 = tmp + 0x17;
		cy_r48 = tmp + 0x18;	cy_r50 = tmp + 0x19;
		cy_r52 = tmp + 0x1a;	cy_r54 = tmp + 0x1b;
		cy_r56 = tmp + 0x1c;	cy_r58 = tmp + 0x1d;
		tmp += 0x1e;
		cy_i00 = tmp + 0x00;	cy_i02 = tmp + 0x01;
		cy_i04 = tmp + 0x02;	cy_i06 = tmp + 0x03;
		cy_i08 = tmp + 0x04;	cy_i10 = tmp + 0x05;
		cy_i12 = tmp + 0x06;	cy_i14 = tmp + 0x07;
		cy_i16 = tmp + 0x08;	cy_i18 = tmp + 0x09;
		cy_i20 = tmp + 0x0a;	cy_i22 = tmp + 0x0b;
		cy_i24 = tmp + 0x0c;	cy_i26 = tmp + 0x0d;
		cy_i28 = tmp + 0x0e;	cy_i30 = tmp + 0x0f;
		cy_i32 = tmp + 0x10;	cy_i34 = tmp + 0x11;
		cy_i36 = tmp + 0x12;	cy_i38 = tmp + 0x13;
		cy_i40 = tmp + 0x14;	cy_i42 = tmp + 0x15;
		cy_i44 = tmp + 0x16;	cy_i46 = tmp + 0x17;
		cy_i48 = tmp + 0x18;	cy_i50 = tmp + 0x19;
		cy_i52 = tmp + 0x1a;	cy_i54 = tmp + 0x1b;
		cy_i56 = tmp + 0x1c;	cy_i58 = tmp + 0x1d;	// +60 = 368 complex
		tmp += 0x1e;
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 370 complex
		// This is where the value of half_arr_offset60 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*odd_radix] x 16 for Fermat-mod */
								// +20 = 390 complex, round up to nearby multiple of 4
	#endif
		ASSERT(HERE, half_arr_offset60 == (uint32)(half_arr-sc_ptr), "half_arr_offset60 mismatches actual!");

		ASSERT(HERE, (radix60_creals_in_local_store << l2_sz_vd) >= ((long)half_arr - (long)s1p00r) + (20 << l2_sz_vd), "radix60_creals_in_local_store checksum failed!");
		/* These remain fixed: */
		VEC_DBL_INIT(sse2_c3m1, c3m1);
		VEC_DBL_INIT(sse2_s   , s   );
		VEC_DBL_INIT(sse2_cn1 , cn1 );
		VEC_DBL_INIT(sse2_cn2 , cn2 );
		VEC_DBL_INIT(sse2_ss3 , ss3 );
		VEC_DBL_INIT(sse2_sn1 , sn1 );
		VEC_DBL_INIT(sse2_sn2 , sn2 );

		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
		VEC_DBL_INIT(sse2_rnd, crnd);

		// Propagate the above consts to the remaining threads:
		nbytes = (int)sse2_rnd - (int)sse2_c3m1 + sz_vd;	// #bytes in above block of data, allowing for 'holes' and assuming only that cc0 and sse2_rnd bookend the block
		tmp = sse2_c3m1;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

		/* SSE2 version of the one_half array - we have a 2-bit lookup, low bit is from the low word of the carry pair,
		high bit from the high, i.e. based on this lookup index [listed with LSB at right], we have:

			index	half_lo	half_hi
			00		1.0		1.0
			01		.50		1.0
			10		1.0		.50
			11		.50		.50

		The inverse-weights computation uses a similar table, but with all entries multiplied by .50:

			index2	half_lo	half_hi
			00		.50		.50
			01		.25		.50
			10		.50		.25
			11		.25		.25

		We do similarly for the base[] and baseinv[] table lookups - each of these get 4 further slots in half_arr.
		We also allocate a further 4 16-byte slots [uninitialized] for storage of the wtl,wtn,wtlp1,wtnm1 locals.

		In 4-way SIMD (AVX) mode, we expand this from 2^2 2-vector table entries to 2^4 4-vector entries.
		*/
		tmp = half_arr;

	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
	#ifdef USE_AVX

		base_negacyclic_root = half_arr + RADIX;

	  #if HIACC
		/*
		The pattern of the negacyclic-DWT-weights ("nDWTs") applied to the RADIX complex outputs of the final-DIT-pass is like so:
		The nDWTs multiplying each set of RADIX DIT DFT outputs are simply the product of a single complex-root "base multiplier" rbase
		(separately computed for each batch of DFT outputs), which "base root" multiplies the (0 - RADIX-1)st [4*RADIX]th roots of unity,
		i.e.
			 rbase * (j*I*Pi/2)/RADIX, for j = 0, ..., RADIX-1 .

		See the radix28 version of this routine for additional details.
		*/
		#if 0
		/* Simple qfloat-based loop to crank out the roots: */
			struct qfloat qc,qs,qx,qy,qt,qn,qmul;
			qt = qfdiv(QPIHALF, i64_to_q((int64)RADIX));	// Pi/2/RADIX
			qc = qfcos(qt);	qs = qfsin(qt);
			qx = QONE;		qy = QZRO;
			for(j = 1; j <= RADIX; j++) {
				// Up-multiply the complex exponential:
				qn = qfmul(qx, qc); qt = qfmul(qy, qs); qmul = qfsub(qn, qt);	// Store qxnew in qmul for now.
				qn = qfmul(qx, qs); qt = qfmul(qy, qc); qy   = qfadd(qn, qt); qx = qmul;
				printf("j = %3u: cos = 0x%16llX\n",j,qfdbl_as_uint64(qx));
			}
			exit(0);
		#endif

		tmp = base_negacyclic_root + RADIX*2;	// First 120 = 2*RADIX slots reserved for RADIX/4 copies of the Re/Im parts of the 4 base multipliers
		tm2 = tmp + RADIX/2 - 1;
										tmp->d0 = 1.0;	(tmp+1)->d0 = 0.0;
		tmp64 = 0x3FEFFD315BBF4275ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(01*I*Pi/120) = sin(59*I*Pi/120) */
		tmp64 = 0x3FEFF4C5ED12E61Dull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(02*I*Pi/120) = sin(58*I*Pi/120) */
		tmp64 = 0x3FEFE6BF2E2660AFull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(03*I*Pi/120) = sin(57*I*Pi/120) */	tmp += 2;
		tmp64 = 0x3FEFD31F94F867C6ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(04*I*Pi/120) = sin(56*I*Pi/120) */	tm2 -= 2;
		tmp64 = 0x3FEFB9EA92EC689Bull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(05*I*Pi/120) = sin(55*I*Pi/120) */
		tmp64 = 0x3FEF9B24942FE45Cull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(06*I*Pi/120) = sin(54*I*Pi/120) */
		tmp64 = 0x3FEF76D2FEF3CC4Bull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(07*I*Pi/120) = sin(53*I*Pi/120) */	tmp += 2;
		tmp64 = 0x3FEF4CFC327A0080ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(08*I*Pi/120) = sin(52*I*Pi/120) */	tm2 -= 2;
		tmp64 = 0x3FEF1DA785F71BCEull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(09*I*Pi/120) = sin(51*I*Pi/120) */
		tmp64 = 0x3FEEE8DD4748BF15ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(10*I*Pi/120) = sin(50*I*Pi/120) */
		tmp64 = 0x3FEEAEA6B98095C0ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(11*I*Pi/120) = sin(49*I*Pi/120) */	tmp += 2;
		tmp64 = 0x3FEE6F0E134454FFull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(12*I*Pi/120) = sin(48*I*Pi/120) */	tm2 -= 2;
		tmp64 = 0x3FEE2A1E7D02FE9Full;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(13*I*Pi/120) = sin(47*I*Pi/120) */
		tmp64 = 0x3FEDDFE40EFFB805ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(14*I*Pi/120) = sin(46*I*Pi/120) */
		tmp64 = 0x3FED906BCF328D46ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(15*I*Pi/120) = sin(45*I*Pi/120) */	tmp += 2;
		tmp64 = 0x3FED3BC3AEFF7F95ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(16*I*Pi/120) = sin(44*I*Pi/120) */	tm2 -= 2;
		tmp64 = 0x3FECE1FA88C445BBull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(17*I*Pi/120) = sin(43*I*Pi/120) */
		tmp64 = 0x3FEC83201D3D2C6Dull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(18*I*Pi/120) = sin(42*I*Pi/120) */
		tmp64 = 0x3FEC1F4510C18B95ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(19*I*Pi/120) = sin(41*I*Pi/120) */	tmp += 2;
		tmp64 = 0x3FEBB67AE8584CAAull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(20*I*Pi/120) = sin(40*I*Pi/120) */	tm2 -= 2;
		tmp64 = 0x3FEB48D406A50540ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(21*I*Pi/120) = sin(39*I*Pi/120) */
		tmp64 = 0x3FEAD663A8AE2FDCull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(22*I*Pi/120) = sin(38*I*Pi/120) */
		tmp64 = 0x3FEA5F3DE27D13F2ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(23*I*Pi/120) = sin(37*I*Pi/120) */	tmp += 2;
		tmp64 = 0x3FE9E3779B97F4A8ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(24*I*Pi/120) = sin(36*I*Pi/120) */	tm2 -= 2;
		tmp64 = 0x3FE963268B572492ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(25*I*Pi/120) = sin(35*I*Pi/120) */
		tmp64 = 0x3FE8DE613515A328ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(26*I*Pi/120) = sin(34*I*Pi/120) */
		tmp64 = 0x3FE8553EE43DEF13ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(27*I*Pi/120) = sin(33*I*Pi/120) */	tmp += 2;
		tmp64 = 0x3FE7C7D7A833BEC2ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(28*I*Pi/120) = sin(32*I*Pi/120) */	tm2 -= 2;
		tmp64 = 0x3FE73644501B56CDull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(29*I*Pi/120) = sin(31*I*Pi/120) */
		tmp64 = 0x3FE6A09E667F3BCDull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(30*I*Pi/120) = sin(30*I*Pi/120) */
		tmp64 = 0x3FE607002CD5031Dull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(31*I*Pi/120) = sin(29*I*Pi/120) */	tmp += 2;
		tmp64 = 0x3FE5698496E20BD8ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(32*I*Pi/120) = sin(28*I*Pi/120) */	tm2 -= 2;
		tmp64 = 0x3FE4C8474600EEEEull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(33*I*Pi/120) = sin(27*I*Pi/120) */
		tmp64 = 0x3FE4236484487ABEull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(34*I*Pi/120) = sin(26*I*Pi/120) */
		tmp64 = 0x3FE37AF93F9513EAull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(35*I*Pi/120) = sin(25*I*Pi/120) */	tmp += 2;
		tmp64 = 0x3FE2CF2304755A5Eull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(36*I*Pi/120) = sin(24*I*Pi/120) */	tm2 -= 2;
		tmp64 = 0x3FE21FFFF8FAF674ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(37*I*Pi/120) = sin(23*I*Pi/120) */
		tmp64 = 0x3FE16DAED770771Dull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(38*I*Pi/120) = sin(22*I*Pi/120) */
		tmp64 = 0x3FE0B84EE8F52E9Dull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(39*I*Pi/120) = sin(21*I*Pi/120) */	tmp += 2;
		tmp64 = 0x3FE0000000000000ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(40*I*Pi/120) = sin(20*I*Pi/120) */	tm2 -= 2;
		tmp64 = 0x3FDE89C4E59427B1ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(41*I*Pi/120) = sin(19*I*Pi/120) */
		tmp64 = 0x3FDD0E2E2B44DE01ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(42*I*Pi/120) = sin(18*I*Pi/120) */
		tmp64 = 0x3FDB8D7E6A56D476ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(43*I*Pi/120) = sin(17*I*Pi/120) */	tmp += 2;
		tmp64 = 0x3FDA07F921061AD1ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(44*I*Pi/120) = sin(16*I*Pi/120) */	tm2 -= 2;
		tmp64 = 0x3FD87DE2A6AEA963ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(45*I*Pi/120) = sin(15*I*Pi/120) */
		tmp64 = 0x3FD6EF801FCED33Cull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(46*I*Pi/120) = sin(14*I*Pi/120) */
		tmp64 = 0x3FD55D1771E5BAB9ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(47*I*Pi/120) = sin(13*I*Pi/120) */	tmp += 2;
		tmp64 = 0x3FD3C6EF372FE950ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(48*I*Pi/120) = sin(12*I*Pi/120) */	tm2 -= 2;
		tmp64 = 0x3FD22D4EB2443163ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(49*I*Pi/120) = sin(11*I*Pi/120) */
		tmp64 = 0x3FD0907DC1930690ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(50*I*Pi/120) = sin(10*I*Pi/120) */
		tmp64 = 0x3FCDE189A594FBCCull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(51*I*Pi/120) = sin(09*I*Pi/120) */	tmp += 2;
		tmp64 = 0x3FCA9CD9AC4258F6ull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(52*I*Pi/120) = sin(08*I*Pi/120) */	tm2 -= 2;
		tmp64 = 0x3FC7537E63143E2Eull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(53*I*Pi/120) = sin(07*I*Pi/120) */
		tmp64 = 0x3FC4060B67A85375ull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(54*I*Pi/120) = sin(06*I*Pi/120) */
		tmp64 = 0x3FC0B5150F6DA2D1ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(55*I*Pi/120) = sin(05*I*Pi/120) */	tmp += 2;
		tmp64 = 0x3FBAC2609B3C576Cull;	tmp->d0 = tm2->d0 = *(double *)&tmp64;	/* cos(56*I*Pi/120) = sin(04*I*Pi/120) */	tm2 -= 2;
		tmp64 = 0x3FB415E532398E49ull;	tmp->d1 = tm2->d3 = *(double *)&tmp64;	/* cos(57*I*Pi/120) = sin(03*I*Pi/120) */
		tmp64 = 0x3FAACBC748EFC90Eull;	tmp->d2 = tm2->d2 = *(double *)&tmp64;	/* cos(58*I*Pi/120) = sin(02*I*Pi/120) */
		tmp64 = 0x3F9ACE214390CA91ull;	tmp->d3 = tm2->d1 = *(double *)&tmp64;	/* cos(59*I*Pi/120) = sin(01*I*Pi/120) */	tmp += 2;

		tmp = base_negacyclic_root + RADIX*2;	// reset to point to start of above block
		nbytes = RADIX*sz_vd/2;	// 7 AVX-register-sized complex data

	  #else	// HIACC = false:

		// lower-precision version, which yields slightly more roundoff error, but is simpler and more storage-compact.
		// Init exp(j*I*Pi/2/RADIX), for j = 0-3:
		tmp = base_negacyclic_root + 8;	// First 8 slots reserved for Re/Im parts of the 4 base multipliers
		tmp->d0 = 1.0;
		tmp64 = 0x3FEFFD315BBF4275ull;	tmp->d1 = *(double *)&tmp64;	// cos(01*I*Pi/120)
		tmp64 = 0x3FEFF4C5ED12E61Dull;	tmp->d2 = *(double *)&tmp64;	// cos(02*I*Pi/120)
		tmp64 = 0x3FEFE6BF2E2660AFull;	tmp->d3 = *(double *)&tmp64;	// cos(03*I*Pi/120)

		(++tmp)->d0 = 0.0;
		tmp64 = 0x3FB415E532398E49ull;	tmp->d3 = *(double *)&tmp64;	// sin(03*I*Pi/120)
		tmp64 = 0x3FAACBC748EFC90Eull;	tmp->d2 = *(double *)&tmp64;	// sin(02*I*Pi/120)
		tmp64 = 0x3F9ACE214390CA91ull;	tmp->d1 = *(double *)&tmp64;	// sin(01*I*Pi/120)
		++tmp;
		tmp64 = 0x3FEFD31F94F867C6ull;	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// cos(04*I*Pi/120)
		++tmp;
		tmp64 = 0x3FBAC2609B3C576Cull;	VEC_DBL_INIT(tmp, *(double *)&tmp64);	// sin(04*I*Pi/120)
		tmp = base_negacyclic_root + 8;	// reset to point to start of above block
		nbytes = 4*sz_vd;	// 2 AVX-register-sized complex data

	  #endif	// HIACC toggle

		// Propagate the above consts to the remaining threads:
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	#endif	// AVX?
	}
	else
	{
	#ifdef USE_AVX
		/* Forward-weight multipliers: */
		tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = 1.0;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = 1.0;	tmp->d3 = .50;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = 1.0;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = 1.0;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		/* Inverse-weight multipliers (only needed for mersenne-mod): */
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .25;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .25;	tmp->d2 = .50;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .25;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .50;	tmp->d2 = .25;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .25;	tmp->d2 = .25;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .25;	tmp->d2 = .25;	tmp->d3 = .50;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .50;	tmp->d2 = .50;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .25;	tmp->d2 = .50;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .25;	tmp->d2 = .50;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .50;	tmp->d2 = .25;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .50;	tmp->d2 = .25;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .50;	tmp->d1 = .25;	tmp->d2 = .25;	tmp->d3 = .25;	++tmp;
		tmp->d0 = .25;	tmp->d1 = .25;	tmp->d2 = .25;	tmp->d3 = .25;	++tmp;
		/* Forward-base[] multipliers: */
		tmp->d0 = base   [0];	tmp->d1 = base   [0];	tmp->d2 = base   [0];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [0];	tmp->d2 = base   [0];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [1];	tmp->d2 = base   [0];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [1];	tmp->d2 = base   [0];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [0];	tmp->d2 = base   [1];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [0];	tmp->d2 = base   [1];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [1];	tmp->d2 = base   [1];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [1];	tmp->d2 = base   [1];	tmp->d3 = base   [0];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [0];	tmp->d2 = base   [0];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [0];	tmp->d2 = base   [0];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [1];	tmp->d2 = base   [0];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [1];	tmp->d2 = base   [0];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [0];	tmp->d2 = base   [1];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [0];	tmp->d2 = base   [1];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [0];	tmp->d1 = base   [1];	tmp->d2 = base   [1];	tmp->d3 = base   [1];	++tmp;
		tmp->d0 = base   [1];	tmp->d1 = base   [1];	tmp->d2 = base   [1];	tmp->d3 = base   [1];	++tmp;
		/* Inverse-base[] multipliers: */
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[0];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[0];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[0];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[0];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[1];	++tmp;
		tmp->d0 = baseinv[1];	tmp->d1 = baseinv[1];	tmp->d2 = baseinv[1];	tmp->d3 = baseinv[1];	++tmp;

		nbytes = 64 << l2_sz_vd;

	#elif defined(USE_SSE2)

		ctmp = (struct complex *)tmp;
		/* Forward-weight multipliers: */
		ctmp->re = 1.0;	ctmp->im = 1.0;	++ctmp;
		ctmp->re = .50;	ctmp->im = 1.0;	++ctmp;
		ctmp->re = 1.0;	ctmp->im = .50;	++ctmp;
		ctmp->re = .50;	ctmp->im = .50;	++ctmp;
		/* Inverse-weight multipliers (only needed for mersenne-mod): */
		ctmp->re = .50;	ctmp->im = .50;	++ctmp;
		ctmp->re = .25;	ctmp->im = .50;	++ctmp;
		ctmp->re = .50;	ctmp->im = .25;	++ctmp;
		ctmp->re = .25;	ctmp->im = .25;	++ctmp;
		/* Forward-base[] multipliers: */
		ctmp->re = base   [0];	ctmp->im = base   [0];	++ctmp;
		ctmp->re = base   [1];	ctmp->im = base   [0];	++ctmp;
		ctmp->re = base   [0];	ctmp->im = base   [1];	++ctmp;
		ctmp->re = base   [1];	ctmp->im = base   [1];	++ctmp;
		/* Inverse-base[] multipliers: */
		ctmp->re = baseinv[0];	ctmp->im = baseinv[0];	++ctmp;
		ctmp->re = baseinv[1];	ctmp->im = baseinv[0];	++ctmp;
		ctmp->re = baseinv[0];	ctmp->im = baseinv[1];	++ctmp;
		ctmp->re = baseinv[1];	ctmp->im = baseinv[1];	++ctmp;

		nbytes = 16 << l2_sz_vd;
	#endif

		// Propagate the above consts to the remaining threads:
		tmp = half_arr;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

	}

	/* Floating-point sign mask used for FABS on packed doubles: */
	sign_mask = sm_ptr;
	for(i = 0; i < RE_IM_STRIDE; ++i) {
		*(sign_mask+i) = (uint64)0x7FFFFFFFFFFFFFFFull;
	}
	// Set up the SIMD-tupled-32-bit-int SSE constants used by the carry macros:
	sse_bw  = sm_ptr + RE_IM_STRIDE;	// (#doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
	tmp64 = (uint64)bw;
	tmp64 = tmp64 + (tmp64 << 32);
	for(i = 0; i < RE_IM_STRIDE; ++i) {
		*(sse_bw+i) = tmp64;
	}

	sse_sw  = sse_bw + RE_IM_STRIDE;
	tmp64 = (uint64)sw;
	tmp64 = tmp64 + (tmp64 << 32);
	for(i = 0; i < RE_IM_STRIDE; ++i) {
		*(sse_sw+i) = tmp64;
	}

	sse_n   = sse_sw + RE_IM_STRIDE;
	tmp64 = (uint64)n;
	tmp64 = tmp64 + (tmp64 << 32);
	for(i = 0; i < RE_IM_STRIDE; ++i) {
		*(sse_n +i) = tmp64;
	}

	nbytes = 4 << l2_sz_vd;

#ifdef USE_AVX
	n_minus_sil   = (struct uint32x4 *)sse_n + 1;
	n_minus_silp1 = (struct uint32x4 *)sse_n + 2;
	sinwt         = (struct uint32x4 *)sse_n + 3;
	sinwtm1       = (struct uint32x4 *)sse_n + 4;
	nbytes += 64;;
#endif

	// Propagate the above consts to the remaining threads:
	tmp = (vec_dbl *)sm_ptr;
	tm2 = tmp + cslots_in_local_store;
	for(ithread = 1; ithread < CY_THREADS; ++ithread) {
		memcpy(tm2, tmp, nbytes);
		tmp = tm2;		tm2 += cslots_in_local_store;
	}

	#ifdef USE_AVX
		bjmodn00 = (int*)(sinwtm1 + RE_IM_STRIDE);
	#else
		bjmodn00 = (int*)(sse_n   + RE_IM_STRIDE);
	#endif
		bjmodn01 = bjmodn00 + 1;
		bjmodn02 = bjmodn01 + 1;
		bjmodn03 = bjmodn02 + 1;
		bjmodn04 = bjmodn03 + 1;
		bjmodn05 = bjmodn04 + 1;
		bjmodn06 = bjmodn05 + 1;
		bjmodn07 = bjmodn06 + 1;
		bjmodn08 = bjmodn07 + 1;
		bjmodn09 = bjmodn08 + 1;
		bjmodn10 = bjmodn09 + 1;
		bjmodn11 = bjmodn10 + 1;
		bjmodn12 = bjmodn11 + 1;
		bjmodn13 = bjmodn12 + 1;
		bjmodn14 = bjmodn13 + 1;
		bjmodn15 = bjmodn14 + 1;
		bjmodn16 = bjmodn15 + 1;
		bjmodn17 = bjmodn16 + 1;
		bjmodn18 = bjmodn17 + 1;
		bjmodn19 = bjmodn18 + 1;
		bjmodn20 = bjmodn19 + 1;
		bjmodn21 = bjmodn20 + 1;
		bjmodn22 = bjmodn21 + 1;
		bjmodn23 = bjmodn22 + 1;
		bjmodn24 = bjmodn23 + 1;
		bjmodn25 = bjmodn24 + 1;
		bjmodn26 = bjmodn25 + 1;
		bjmodn27 = bjmodn26 + 1;
		bjmodn28 = bjmodn27 + 1;
		bjmodn29 = bjmodn28 + 1;
		bjmodn30 = bjmodn29 + 1;
		bjmodn31 = bjmodn30 + 1;
		bjmodn32 = bjmodn31 + 1;
		bjmodn33 = bjmodn32 + 1;
		bjmodn34 = bjmodn33 + 1;
		bjmodn35 = bjmodn34 + 1;
		bjmodn36 = bjmodn35 + 1;
		bjmodn37 = bjmodn36 + 1;
		bjmodn38 = bjmodn37 + 1;
		bjmodn39 = bjmodn38 + 1;
		bjmodn40 = bjmodn39 + 1;
		bjmodn41 = bjmodn40 + 1;
		bjmodn42 = bjmodn41 + 1;
		bjmodn43 = bjmodn42 + 1;
		bjmodn44 = bjmodn43 + 1;
		bjmodn45 = bjmodn44 + 1;
		bjmodn46 = bjmodn45 + 1;
		bjmodn47 = bjmodn46 + 1;
		bjmodn48 = bjmodn47 + 1;
		bjmodn49 = bjmodn48 + 1;
		bjmodn50 = bjmodn49 + 1;
		bjmodn51 = bjmodn50 + 1;
		bjmodn52 = bjmodn51 + 1;
		bjmodn53 = bjmodn52 + 1;
		bjmodn54 = bjmodn53 + 1;
		bjmodn55 = bjmodn54 + 1;
		bjmodn56 = bjmodn55 + 1;
		bjmodn57 = bjmodn56 + 1;
		bjmodn58 = bjmodn57 + 1;
		bjmodn59 = bjmodn58 + 1;

	/*********** Defer the per-thread local-mem-block copy until after added wts-index precomputation below ************/
	#endif	/* USE_SSE2 */

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
		#ifdef USE_SSE2
			tdat[ithread].s1p00r   = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].r00      = (long)tdat[ithread].s1p00r + ((long)r00 - (long)s1p00r);
			tdat[ithread].half_arr = (long)tdat[ithread].s1p00r + ((long)half_arr - (long)s1p00r);
		#else
			// In scalar mode use these 2 ptrs to pass wts_idx_incr and the base/baseinv/etc array-ptrs:
			tdat[ithread].s1p00r   = (vec_dbl *)foo_array;
			tdat[ithread].half_arr = (vec_dbl *)&wts_idx_incr;
		#endif	// USE_SSE2
		}
	#endif

		/*   constant index offsets for array load/stores are here.	*/
		pini = NDIVR/CY_THREADS;
		pini += ( (pini >> DAT_BITS) << PAD_BITS );
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;
		p28 = p24 + p04;
		p32 = p28 + p04;
		p36 = p32 + p04;
		p40 = p36 + p04;
		p44 = p40 + p04;
		p48 = p44 + p04;
		p52 = p48 + p04;
		p56 = p52 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p36 = p36 + ( (p36 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		p44 = p44 + ( (p44 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p52 = p52 + ( (p52 >> DAT_BITS) << PAD_BITS );
		p56 = p56 + ( (p56 >> DAT_BITS) << PAD_BITS );

		if(_cy_r00)	/* If it's a new exponent of a range test, need to deallocate these. */
		{
			free((void *)_i     ); _i      = 0x0;

			free((void *)_bjmodn00); _bjmodn00 = 0x0;
			free((void *)_bjmodn01); _bjmodn01 = 0x0;
			free((void *)_bjmodn02); _bjmodn02 = 0x0;
			free((void *)_bjmodn03); _bjmodn03 = 0x0;
			free((void *)_bjmodn04); _bjmodn04 = 0x0;
			free((void *)_bjmodn05); _bjmodn05 = 0x0;
			free((void *)_bjmodn06); _bjmodn06 = 0x0;
			free((void *)_bjmodn07); _bjmodn07 = 0x0;
			free((void *)_bjmodn08); _bjmodn08 = 0x0;
			free((void *)_bjmodn09); _bjmodn09 = 0x0;
			free((void *)_bjmodn10); _bjmodn10 = 0x0;
			free((void *)_bjmodn11); _bjmodn11 = 0x0;
			free((void *)_bjmodn12); _bjmodn12 = 0x0;
			free((void *)_bjmodn13); _bjmodn13 = 0x0;
			free((void *)_bjmodn14); _bjmodn14 = 0x0;
			free((void *)_bjmodn15); _bjmodn15 = 0x0;
			free((void *)_bjmodn16); _bjmodn16 = 0x0;
			free((void *)_bjmodn17); _bjmodn17 = 0x0;
			free((void *)_bjmodn18); _bjmodn18 = 0x0;
			free((void *)_bjmodn19); _bjmodn19 = 0x0;
			free((void *)_bjmodn20); _bjmodn20 = 0x0;
			free((void *)_bjmodn21); _bjmodn21 = 0x0;
			free((void *)_bjmodn22); _bjmodn22 = 0x0;
			free((void *)_bjmodn23); _bjmodn23 = 0x0;
			free((void *)_bjmodn24); _bjmodn24 = 0x0;
			free((void *)_bjmodn25); _bjmodn25 = 0x0;
			free((void *)_bjmodn26); _bjmodn26 = 0x0;
			free((void *)_bjmodn27); _bjmodn27 = 0x0;
			free((void *)_bjmodn28); _bjmodn28 = 0x0;
			free((void *)_bjmodn29); _bjmodn29 = 0x0;
			free((void *)_bjmodn30); _bjmodn30 = 0x0;
			free((void *)_bjmodn31); _bjmodn31 = 0x0;
			free((void *)_bjmodn32); _bjmodn32 = 0x0;
			free((void *)_bjmodn33); _bjmodn33 = 0x0;
			free((void *)_bjmodn34); _bjmodn34 = 0x0;
			free((void *)_bjmodn35); _bjmodn35 = 0x0;
			free((void *)_bjmodn36); _bjmodn36 = 0x0;
			free((void *)_bjmodn37); _bjmodn37 = 0x0;
			free((void *)_bjmodn38); _bjmodn38 = 0x0;
			free((void *)_bjmodn39); _bjmodn39 = 0x0;
			free((void *)_bjmodn40); _bjmodn40 = 0x0;
			free((void *)_bjmodn41); _bjmodn41 = 0x0;
			free((void *)_bjmodn42); _bjmodn42 = 0x0;
			free((void *)_bjmodn43); _bjmodn43 = 0x0;
			free((void *)_bjmodn44); _bjmodn44 = 0x0;
			free((void *)_bjmodn45); _bjmodn45 = 0x0;
			free((void *)_bjmodn46); _bjmodn46 = 0x0;
			free((void *)_bjmodn47); _bjmodn47 = 0x0;
			free((void *)_bjmodn48); _bjmodn48 = 0x0;
			free((void *)_bjmodn49); _bjmodn49 = 0x0;
			free((void *)_bjmodn50); _bjmodn50 = 0x0;
			free((void *)_bjmodn51); _bjmodn51 = 0x0;
			free((void *)_bjmodn52); _bjmodn52 = 0x0;
			free((void *)_bjmodn53); _bjmodn53 = 0x0;
			free((void *)_bjmodn54); _bjmodn54 = 0x0;
			free((void *)_bjmodn55); _bjmodn55 = 0x0;
			free((void *)_bjmodn56); _bjmodn56 = 0x0;
			free((void *)_bjmodn57); _bjmodn57 = 0x0;
			free((void *)_bjmodn58); _bjmodn58 = 0x0;
			free((void *)_bjmodn59); _bjmodn59 = 0x0;

			free((void *)_cy_r00); _cy_r00 = 0x0;		free((void *)_cy_i00); _cy_i00 = 0x0;
			free((void *)_cy_r01); _cy_r01 = 0x0;		free((void *)_cy_i01); _cy_i01 = 0x0;
			free((void *)_cy_r02); _cy_r02 = 0x0;		free((void *)_cy_i02); _cy_i02 = 0x0;
			free((void *)_cy_r03); _cy_r03 = 0x0;		free((void *)_cy_i03); _cy_i03 = 0x0;
			free((void *)_cy_r04); _cy_r04 = 0x0;		free((void *)_cy_i04); _cy_i04 = 0x0;
			free((void *)_cy_r05); _cy_r05 = 0x0;		free((void *)_cy_i05); _cy_i05 = 0x0;
			free((void *)_cy_r06); _cy_r06 = 0x0;		free((void *)_cy_i06); _cy_i06 = 0x0;
			free((void *)_cy_r07); _cy_r07 = 0x0;		free((void *)_cy_i07); _cy_i07 = 0x0;
			free((void *)_cy_r08); _cy_r08 = 0x0;		free((void *)_cy_i08); _cy_i08 = 0x0;
			free((void *)_cy_r09); _cy_r09 = 0x0;		free((void *)_cy_i09); _cy_i09 = 0x0;
			free((void *)_cy_r10); _cy_r10 = 0x0;		free((void *)_cy_i10); _cy_i10 = 0x0;
			free((void *)_cy_r11); _cy_r11 = 0x0;		free((void *)_cy_i11); _cy_i11 = 0x0;
			free((void *)_cy_r12); _cy_r12 = 0x0;		free((void *)_cy_i12); _cy_i12 = 0x0;
			free((void *)_cy_r13); _cy_r13 = 0x0;		free((void *)_cy_i13); _cy_i13 = 0x0;
			free((void *)_cy_r14); _cy_r14 = 0x0;		free((void *)_cy_i14); _cy_i14 = 0x0;
			free((void *)_cy_r15); _cy_r15 = 0x0;		free((void *)_cy_i15); _cy_i15 = 0x0;
			free((void *)_cy_r16); _cy_r16 = 0x0;		free((void *)_cy_i16); _cy_i16 = 0x0;
			free((void *)_cy_r17); _cy_r17 = 0x0;		free((void *)_cy_i17); _cy_i17 = 0x0;
			free((void *)_cy_r18); _cy_r18 = 0x0;		free((void *)_cy_i18); _cy_i18 = 0x0;
			free((void *)_cy_r19); _cy_r19 = 0x0;		free((void *)_cy_i19); _cy_i19 = 0x0;
			free((void *)_cy_r20); _cy_r20 = 0x0;		free((void *)_cy_i20); _cy_i20 = 0x0;
			free((void *)_cy_r21); _cy_r21 = 0x0;		free((void *)_cy_i21); _cy_i21 = 0x0;
			free((void *)_cy_r22); _cy_r22 = 0x0;		free((void *)_cy_i22); _cy_i22 = 0x0;
			free((void *)_cy_r23); _cy_r23 = 0x0;		free((void *)_cy_i23); _cy_i23 = 0x0;
			free((void *)_cy_r24); _cy_r24 = 0x0;		free((void *)_cy_i24); _cy_i24 = 0x0;
			free((void *)_cy_r25); _cy_r25 = 0x0;		free((void *)_cy_i25); _cy_i25 = 0x0;
			free((void *)_cy_r26); _cy_r26 = 0x0;		free((void *)_cy_i26); _cy_i26 = 0x0;
			free((void *)_cy_r27); _cy_r27 = 0x0;		free((void *)_cy_i27); _cy_i27 = 0x0;
			free((void *)_cy_r28); _cy_r28 = 0x0;		free((void *)_cy_i28); _cy_i28 = 0x0;
			free((void *)_cy_r29); _cy_r29 = 0x0;		free((void *)_cy_i29); _cy_i29 = 0x0;
			free((void *)_cy_r30); _cy_r30 = 0x0;		free((void *)_cy_i30); _cy_i30 = 0x0;
			free((void *)_cy_r31); _cy_r31 = 0x0;		free((void *)_cy_i31); _cy_i31 = 0x0;
			free((void *)_cy_r32); _cy_r32 = 0x0;		free((void *)_cy_i32); _cy_i32 = 0x0;
			free((void *)_cy_r33); _cy_r33 = 0x0;		free((void *)_cy_i33); _cy_i33 = 0x0;
			free((void *)_cy_r34); _cy_r34 = 0x0;		free((void *)_cy_i34); _cy_i34 = 0x0;
			free((void *)_cy_r35); _cy_r35 = 0x0;		free((void *)_cy_i35); _cy_i35 = 0x0;
			free((void *)_cy_r36); _cy_r36 = 0x0;		free((void *)_cy_i36); _cy_i36 = 0x0;
			free((void *)_cy_r37); _cy_r37 = 0x0;		free((void *)_cy_i37); _cy_i37 = 0x0;
			free((void *)_cy_r38); _cy_r38 = 0x0;		free((void *)_cy_i38); _cy_i38 = 0x0;
			free((void *)_cy_r39); _cy_r39 = 0x0;		free((void *)_cy_i39); _cy_i39 = 0x0;
			free((void *)_cy_r40); _cy_r40 = 0x0;		free((void *)_cy_i40); _cy_i40 = 0x0;
			free((void *)_cy_r41); _cy_r41 = 0x0;		free((void *)_cy_i41); _cy_i41 = 0x0;
			free((void *)_cy_r42); _cy_r42 = 0x0;		free((void *)_cy_i42); _cy_i42 = 0x0;
			free((void *)_cy_r43); _cy_r43 = 0x0;		free((void *)_cy_i43); _cy_i43 = 0x0;
			free((void *)_cy_r44); _cy_r44 = 0x0;		free((void *)_cy_i44); _cy_i44 = 0x0;
			free((void *)_cy_r45); _cy_r45 = 0x0;		free((void *)_cy_i45); _cy_i45 = 0x0;
			free((void *)_cy_r46); _cy_r46 = 0x0;		free((void *)_cy_i46); _cy_i46 = 0x0;
			free((void *)_cy_r47); _cy_r47 = 0x0;		free((void *)_cy_i47); _cy_i47 = 0x0;
			free((void *)_cy_r48); _cy_r48 = 0x0;		free((void *)_cy_i48); _cy_i48 = 0x0;
			free((void *)_cy_r49); _cy_r49 = 0x0;		free((void *)_cy_i49); _cy_i49 = 0x0;
			free((void *)_cy_r50); _cy_r50 = 0x0;		free((void *)_cy_i50); _cy_i50 = 0x0;
			free((void *)_cy_r51); _cy_r51 = 0x0;		free((void *)_cy_i51); _cy_i51 = 0x0;
			free((void *)_cy_r52); _cy_r52 = 0x0;		free((void *)_cy_i52); _cy_i52 = 0x0;
			free((void *)_cy_r53); _cy_r53 = 0x0;		free((void *)_cy_i53); _cy_i53 = 0x0;
			free((void *)_cy_r54); _cy_r54 = 0x0;		free((void *)_cy_i54); _cy_i54 = 0x0;
			free((void *)_cy_r55); _cy_r55 = 0x0;		free((void *)_cy_i55); _cy_i55 = 0x0;
			free((void *)_cy_r56); _cy_r56 = 0x0;		free((void *)_cy_i56); _cy_i56 = 0x0;
			free((void *)_cy_r57); _cy_r57 = 0x0;		free((void *)_cy_i57); _cy_i57 = 0x0;
			free((void *)_cy_r58); _cy_r58 = 0x0;		free((void *)_cy_i58); _cy_i58 = 0x0;
			free((void *)_cy_r59); _cy_r59 = 0x0;		free((void *)_cy_i59); _cy_i59 = 0x0;

			free((void *)_jstart ); _jstart  = 0x0;
			free((void *)_jhi    ); _jhi     = 0x0;
			free((void *)_maxerr); _maxerr = 0x0;
			free((void *)_col   ); _col    = 0x0;
			free((void *)_co2   ); _co2    = 0x0;
			free((void *)_co3   ); _co3    = 0x0;

			free((void *)_bjmodnini); _bjmodnini = 0x0;
		}

		ptr_prod = (uint32)0;	/* Store bitmask for allocatable-array ptrs here, check vs 0 after all alloc calls finish */
		j = CY_THREADS*sizeof(int);
		_i       	= (int *)malloc(j);	ptr_prod += (uint32)(_i== 0x0);
		_bjmodn00	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn00== 0x0);
		_bjmodn01	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn01== 0x0);
		_bjmodn02	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn02== 0x0);
		_bjmodn03	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn03== 0x0);
		_bjmodn04	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn04== 0x0);
		_bjmodn05	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn05== 0x0);
		_bjmodn06	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn06== 0x0);
		_bjmodn07	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn07== 0x0);
		_bjmodn08	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn08== 0x0);
		_bjmodn09	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn09== 0x0);
		_bjmodn10	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn10== 0x0);
		_bjmodn11	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn11== 0x0);
		_bjmodn12	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn12== 0x0);
		_bjmodn13	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn13== 0x0);
		_bjmodn14	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn14== 0x0);
		_bjmodn15	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn15== 0x0);
		_bjmodn16	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn16== 0x0);
		_bjmodn17	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn17== 0x0);
		_bjmodn18	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn18== 0x0);
		_bjmodn19	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn19== 0x0);
		_bjmodn20	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn20== 0x0);
		_bjmodn21	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn21== 0x0);
		_bjmodn22	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn22== 0x0);
		_bjmodn23	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn23== 0x0);
		_bjmodn24	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn24== 0x0);
		_bjmodn25	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn25== 0x0);
		_bjmodn26	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn26== 0x0);
		_bjmodn27	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn27== 0x0);
		_bjmodn28	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn28== 0x0);
		_bjmodn29	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn29== 0x0);
		_bjmodn30	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn30== 0x0);
		_bjmodn31	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn31== 0x0);
		_bjmodn32	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn32== 0x0);
		_bjmodn33	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn33== 0x0);
		_bjmodn34	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn34== 0x0);
		_bjmodn35	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn35== 0x0);
		_bjmodn36	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn36== 0x0);
		_bjmodn37	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn37== 0x0);
		_bjmodn38	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn38== 0x0);
		_bjmodn39	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn39== 0x0);
		_bjmodn40	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn40== 0x0);
		_bjmodn41	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn41== 0x0);
		_bjmodn42	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn42== 0x0);
		_bjmodn43	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn43== 0x0);
		_bjmodn44	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn44== 0x0);
		_bjmodn45	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn45== 0x0);
		_bjmodn46	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn46== 0x0);
		_bjmodn47	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn47== 0x0);
		_bjmodn48	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn48== 0x0);
		_bjmodn49	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn49== 0x0);
		_bjmodn50	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn50== 0x0);
		_bjmodn51	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn51== 0x0);
		_bjmodn52	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn52== 0x0);
		_bjmodn53	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn53== 0x0);
		_bjmodn54	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn54== 0x0);
		_bjmodn55	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn55== 0x0);
		_bjmodn56	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn56== 0x0);
		_bjmodn57	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn57== 0x0);
		_bjmodn58	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn58== 0x0);
		_bjmodn59	= (int *)malloc(j);	ptr_prod += (uint32)(_bjmodn59== 0x0);
		_jstart  	= (int *)malloc(j);	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(j);	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(j);	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(j);	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(j);	ptr_prod += (uint32)(_co3     == 0x0);

		j = CY_THREADS*sizeof(double);
		_cy_r00	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r00== 0x0);
		_cy_r01	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r01== 0x0);
		_cy_r02	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r02== 0x0);
		_cy_r03	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r03== 0x0);
		_cy_r04	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r04== 0x0);
		_cy_r05	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r05== 0x0);
		_cy_r06	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r06== 0x0);
		_cy_r07	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r07== 0x0);
		_cy_r08	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r08== 0x0);
		_cy_r09	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r09== 0x0);
		_cy_r10	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r10== 0x0);
		_cy_r11	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r11== 0x0);
		_cy_r12	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r12== 0x0);
		_cy_r13	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r13== 0x0);
		_cy_r14	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r14== 0x0);
		_cy_r15	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r15== 0x0);
		_cy_r16	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r16== 0x0);
		_cy_r17	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r17== 0x0);
		_cy_r18	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r18== 0x0);
		_cy_r19	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r19== 0x0);
		_cy_r20	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r20== 0x0);
		_cy_r21	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r21== 0x0);
		_cy_r22	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r22== 0x0);
		_cy_r23	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r23== 0x0);
		_cy_r24	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r24== 0x0);
		_cy_r25	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r25== 0x0);
		_cy_r26	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r26== 0x0);
		_cy_r27	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r27== 0x0);
		_cy_r28	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r28== 0x0);
		_cy_r29	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r29== 0x0);
		_cy_r30	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r30== 0x0);
		_cy_r31	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r31== 0x0);
		_cy_r32	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r32== 0x0);
		_cy_r33	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r33== 0x0);
		_cy_r34	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r34== 0x0);
		_cy_r35	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r35== 0x0);
		_cy_r36	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r36== 0x0);
		_cy_r37	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r37== 0x0);
		_cy_r38	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r38== 0x0);
		_cy_r39	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r39== 0x0);
		_cy_r40	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r40== 0x0);
		_cy_r41	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r41== 0x0);
		_cy_r42	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r42== 0x0);
		_cy_r43	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r43== 0x0);
		_cy_r44	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r44== 0x0);
		_cy_r45	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r45== 0x0);
		_cy_r46	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r46== 0x0);
		_cy_r47	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r47== 0x0);
		_cy_r48	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r48== 0x0);
		_cy_r49	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r49== 0x0);
		_cy_r50	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r50== 0x0);
		_cy_r51	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r51== 0x0);
		_cy_r52	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r52== 0x0);
		_cy_r53	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r53== 0x0);
		_cy_r54	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r54== 0x0);
		_cy_r55	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r55== 0x0);
		_cy_r56	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r56== 0x0);
		_cy_r57	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r57== 0x0);
		_cy_r58	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r58== 0x0);
		_cy_r59	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_r59== 0x0);

		_cy_i00	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i00== 0x0);
		_cy_i01	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i01== 0x0);
		_cy_i02	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i02== 0x0);
		_cy_i03	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i03== 0x0);
		_cy_i04	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i04== 0x0);
		_cy_i05	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i05== 0x0);
		_cy_i06	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i06== 0x0);
		_cy_i07	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i07== 0x0);
		_cy_i08	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i08== 0x0);
		_cy_i09	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i09== 0x0);
		_cy_i10	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i10== 0x0);
		_cy_i11	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i11== 0x0);
		_cy_i12	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i12== 0x0);
		_cy_i13	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i13== 0x0);
		_cy_i14	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i14== 0x0);
		_cy_i15	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i15== 0x0);
		_cy_i16	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i16== 0x0);
		_cy_i17	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i17== 0x0);
		_cy_i18	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i18== 0x0);
		_cy_i19	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i19== 0x0);
		_cy_i20	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i20== 0x0);
		_cy_i21	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i21== 0x0);
		_cy_i22	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i22== 0x0);
		_cy_i23	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i23== 0x0);
		_cy_i24	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i24== 0x0);
		_cy_i25	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i25== 0x0);
		_cy_i26	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i26== 0x0);
		_cy_i27	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i27== 0x0);
		_cy_i28	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i28== 0x0);
		_cy_i29	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i29== 0x0);
		_cy_i30	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i30== 0x0);
		_cy_i31	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i31== 0x0);
		_cy_i32	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i32== 0x0);
		_cy_i33	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i33== 0x0);
		_cy_i34	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i34== 0x0);
		_cy_i35	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i35== 0x0);
		_cy_i36	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i36== 0x0);
		_cy_i37	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i37== 0x0);
		_cy_i38	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i38== 0x0);
		_cy_i39	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i39== 0x0);
		_cy_i40	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i40== 0x0);
		_cy_i41	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i41== 0x0);
		_cy_i42	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i42== 0x0);
		_cy_i43	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i43== 0x0);
		_cy_i44	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i44== 0x0);
		_cy_i45	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i45== 0x0);
		_cy_i46	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i46== 0x0);
		_cy_i47	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i47== 0x0);
		_cy_i48	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i48== 0x0);
		_cy_i49	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i49== 0x0);
		_cy_i50	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i50== 0x0);
		_cy_i51	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i51== 0x0);
		_cy_i52	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i52== 0x0);
		_cy_i53	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i53== 0x0);
		_cy_i54	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i54== 0x0);
		_cy_i55	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i55== 0x0);
		_cy_i56	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i56== 0x0);
		_cy_i57	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i57== 0x0);
		_cy_i58	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i58== 0x0);
		_cy_i59	= (double *)malloc(j);	ptr_prod += (uint32)(_cy_i59== 0x0);

		_maxerr	= (double *)malloc(j);	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/radix-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini in %s.\n", func); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		_bjmodnini[0] = 0;
		_bjmodnini[1] = 0;

		/* For Fermat-mod, since 'adjacent' words are actually stride-2 separated
		in terms of the floating residue array, block boundaries have half the i-index
		(e.g. as in sw*i and bw*i) value they do in the Mersenne-mod case:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			jhi = NDIVR/CY_THREADS;
		}
		else
		{
			jhi = NDIVR/CY_THREADS/2;
		}

		for(j=0; j < jhi; j++)
		{
			_bjmodnini[1] -= sw; _bjmodnini[1] = _bjmodnini[1] + ( (-(int)((uint32)_bjmodnini[1] >> 31)) & n);
		}

		if(CY_THREADS > 1)
		{
			for(ithread = 2; ithread <= CY_THREADS; ithread++)
			{
				_bjmodnini[ithread] = _bjmodnini[ithread-1] + _bjmodnini[1] - n; _bjmodnini[ithread] = _bjmodnini[ithread] + ( (-(int)((uint32)_bjmodnini[ithread] >> 31)) & n);
			}
		}
		/* Check upper element against scalar value, as precomputed in single-thread mode: */
		bjmodnini=0;
		for(j=0; j < jhi*CY_THREADS; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
		ASSERT(HERE, _bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");

	}	/* endif(first_entry) */

/*...The radix-60 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_cy_r00[ithread] = 0;	_cy_i00[ithread] = 0;
		_cy_r01[ithread] = 0;	_cy_i01[ithread] = 0;
		_cy_r02[ithread] = 0;	_cy_i02[ithread] = 0;
		_cy_r03[ithread] = 0;	_cy_i03[ithread] = 0;
		_cy_r04[ithread] = 0;	_cy_i04[ithread] = 0;
		_cy_r05[ithread] = 0;	_cy_i05[ithread] = 0;
		_cy_r06[ithread] = 0;	_cy_i06[ithread] = 0;
		_cy_r07[ithread] = 0;	_cy_i07[ithread] = 0;
		_cy_r08[ithread] = 0;	_cy_i08[ithread] = 0;
		_cy_r09[ithread] = 0;	_cy_i09[ithread] = 0;
		_cy_r10[ithread] = 0;	_cy_i10[ithread] = 0;
		_cy_r11[ithread] = 0;	_cy_i11[ithread] = 0;
		_cy_r12[ithread] = 0;	_cy_i12[ithread] = 0;
		_cy_r13[ithread] = 0;	_cy_i13[ithread] = 0;
		_cy_r14[ithread] = 0;	_cy_i14[ithread] = 0;
		_cy_r15[ithread] = 0;	_cy_i15[ithread] = 0;
		_cy_r16[ithread] = 0;	_cy_i16[ithread] = 0;
		_cy_r17[ithread] = 0;	_cy_i17[ithread] = 0;
		_cy_r18[ithread] = 0;	_cy_i18[ithread] = 0;
		_cy_r19[ithread] = 0;	_cy_i19[ithread] = 0;
		_cy_r20[ithread] = 0;	_cy_i20[ithread] = 0;
		_cy_r21[ithread] = 0;	_cy_i21[ithread] = 0;
		_cy_r22[ithread] = 0;	_cy_i22[ithread] = 0;
		_cy_r23[ithread] = 0;	_cy_i23[ithread] = 0;
		_cy_r24[ithread] = 0;	_cy_i24[ithread] = 0;
		_cy_r25[ithread] = 0;	_cy_i25[ithread] = 0;
		_cy_r26[ithread] = 0;	_cy_i26[ithread] = 0;
		_cy_r27[ithread] = 0;	_cy_i27[ithread] = 0;
		_cy_r28[ithread] = 0;	_cy_i28[ithread] = 0;
		_cy_r29[ithread] = 0;	_cy_i29[ithread] = 0;
		_cy_r30[ithread] = 0;	_cy_i30[ithread] = 0;
		_cy_r31[ithread] = 0;	_cy_i31[ithread] = 0;
		_cy_r32[ithread] = 0;	_cy_i32[ithread] = 0;
		_cy_r33[ithread] = 0;	_cy_i33[ithread] = 0;
		_cy_r34[ithread] = 0;	_cy_i34[ithread] = 0;
		_cy_r35[ithread] = 0;	_cy_i35[ithread] = 0;
		_cy_r36[ithread] = 0;	_cy_i36[ithread] = 0;
		_cy_r37[ithread] = 0;	_cy_i37[ithread] = 0;
		_cy_r38[ithread] = 0;	_cy_i38[ithread] = 0;
		_cy_r39[ithread] = 0;	_cy_i39[ithread] = 0;
		_cy_r40[ithread] = 0;	_cy_i40[ithread] = 0;
		_cy_r41[ithread] = 0;	_cy_i41[ithread] = 0;
		_cy_r42[ithread] = 0;	_cy_i42[ithread] = 0;
		_cy_r43[ithread] = 0;	_cy_i43[ithread] = 0;
		_cy_r44[ithread] = 0;	_cy_i44[ithread] = 0;
		_cy_r45[ithread] = 0;	_cy_i45[ithread] = 0;
		_cy_r46[ithread] = 0;	_cy_i46[ithread] = 0;
		_cy_r47[ithread] = 0;	_cy_i47[ithread] = 0;
		_cy_r48[ithread] = 0;	_cy_i48[ithread] = 0;
		_cy_r49[ithread] = 0;	_cy_i49[ithread] = 0;
		_cy_r50[ithread] = 0;	_cy_i50[ithread] = 0;
		_cy_r51[ithread] = 0;	_cy_i51[ithread] = 0;
		_cy_r52[ithread] = 0;	_cy_i52[ithread] = 0;
		_cy_r53[ithread] = 0;	_cy_i53[ithread] = 0;
		_cy_r54[ithread] = 0;	_cy_i54[ithread] = 0;
		_cy_r55[ithread] = 0;	_cy_i55[ithread] = 0;
		_cy_r56[ithread] = 0;	_cy_i56[ithread] = 0;
		_cy_r57[ithread] = 0;	_cy_i57[ithread] = 0;
		_cy_r58[ithread] = 0;	_cy_i58[ithread] = 0;
		_cy_r59[ithread] = 0;	_cy_i59[ithread] = 0;
	}
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy_r00[      0] = -2;
	}

	*fracmax=0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = 0.0;
	}

for(outer=0; outer <= 1; outer++)
{
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (_i[0] = 1).	*/

		khi = n_div_nwt/CY_THREADS;
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_jstart[ithread] = ithread*NDIVR/CY_THREADS;
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + nwt-1;

			_col[ithread] = ithread*(khi*RADIX);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
			_co2[ithread] = (n>>nwt_bits)-1+RADIX - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
			_co3[ithread] = _co2[ithread]-RADIX;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
		}
	}
	else
	{
		_i[0] = 0;		/* Pointer to the BASE and BASEINV arrays. If n divides p, lowest-order digit is always a smallword (_i[0] = 0).	*/

		khi = 1;
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_jstart[ithread] = ithread*NDIVR/CY_THREADS;
			/*
			For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
			*/
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + 15;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + n_div_nwt/CY_THREADS;
		}

		/* For Fermat-mod, IBDWT access patterns repeat with period NWT = {odd part of radix0},
		so for even radix0 values still only need [radix0 >> trailz(radix0)] bjmodn and ii's:
		*/
		/* indices into IBDWT weights arrays (mod NWT) is here: */
		ii00= 0;
		ii01= (SW_DIV_N*(NDIVR >> 1)) % nwt;	// nwt *not* a power of 2, must use library-mod!
		MOD_ADD32(ii01,ii01,nwt,ii02);
		MOD_ADD32(ii02,ii01,nwt,ii03);
		MOD_ADD32(ii03,ii01,nwt,ii04);
		MOD_ADD32(ii04,ii01,nwt,ii05);
		MOD_ADD32(ii05,ii01,nwt,ii06);
		MOD_ADD32(ii06,ii01,nwt,ii07);
		MOD_ADD32(ii07,ii01,nwt,ii08);
		MOD_ADD32(ii08,ii01,nwt,ii09);
		MOD_ADD32(ii09,ii01,nwt,ii10);
		MOD_ADD32(ii10,ii01,nwt,ii11);
		MOD_ADD32(ii11,ii01,nwt,ii12);
		MOD_ADD32(ii12,ii01,nwt,ii13);
		MOD_ADD32(ii13,ii01,nwt,ii14);
	}

	// In non-power-of-2-runlength case, both Mersenne and Fermat-mod share these next 2 loops:
	if(CY_THREADS > 1)
	{
		for(ithread = 1; ithread < CY_THREADS; ithread++)
		{
			_i[ithread] = ((uint32)(sw - _bjmodnini[ithread]) >> 31);
		}
	}

	// Include 0-thread here ... bjmodn terms all 0 for that, but need jhi computed for all threads:
	j = _bjmodnini[CY_THREADS];
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_bjmodn00[ithread] = _bjmodnini[ithread];
		MOD_ADD32(_bjmodn00[ithread], j, n, _bjmodn01[ithread]);
		MOD_ADD32(_bjmodn01[ithread], j, n, _bjmodn02[ithread]);
		MOD_ADD32(_bjmodn02[ithread], j, n, _bjmodn03[ithread]);
		MOD_ADD32(_bjmodn03[ithread], j, n, _bjmodn04[ithread]);
		MOD_ADD32(_bjmodn04[ithread], j, n, _bjmodn05[ithread]);
		MOD_ADD32(_bjmodn05[ithread], j, n, _bjmodn06[ithread]);
		MOD_ADD32(_bjmodn06[ithread], j, n, _bjmodn07[ithread]);
		MOD_ADD32(_bjmodn07[ithread], j, n, _bjmodn08[ithread]);
		MOD_ADD32(_bjmodn08[ithread], j, n, _bjmodn09[ithread]);
		MOD_ADD32(_bjmodn09[ithread], j, n, _bjmodn10[ithread]);
		MOD_ADD32(_bjmodn10[ithread], j, n, _bjmodn11[ithread]);
		MOD_ADD32(_bjmodn11[ithread], j, n, _bjmodn12[ithread]);
		MOD_ADD32(_bjmodn12[ithread], j, n, _bjmodn13[ithread]);
		MOD_ADD32(_bjmodn13[ithread], j, n, _bjmodn14[ithread]);
		MOD_ADD32(_bjmodn14[ithread], j, n, _bjmodn15[ithread]);
		MOD_ADD32(_bjmodn15[ithread], j, n, _bjmodn16[ithread]);
		MOD_ADD32(_bjmodn16[ithread], j, n, _bjmodn17[ithread]);
		MOD_ADD32(_bjmodn17[ithread], j, n, _bjmodn18[ithread]);
		MOD_ADD32(_bjmodn18[ithread], j, n, _bjmodn19[ithread]);
		MOD_ADD32(_bjmodn19[ithread], j, n, _bjmodn20[ithread]);
		MOD_ADD32(_bjmodn20[ithread], j, n, _bjmodn21[ithread]);
		MOD_ADD32(_bjmodn21[ithread], j, n, _bjmodn22[ithread]);
		MOD_ADD32(_bjmodn22[ithread], j, n, _bjmodn23[ithread]);
		MOD_ADD32(_bjmodn23[ithread], j, n, _bjmodn24[ithread]);
		MOD_ADD32(_bjmodn24[ithread], j, n, _bjmodn25[ithread]);
		MOD_ADD32(_bjmodn25[ithread], j, n, _bjmodn26[ithread]);
		MOD_ADD32(_bjmodn26[ithread], j, n, _bjmodn27[ithread]);
		MOD_ADD32(_bjmodn27[ithread], j, n, _bjmodn28[ithread]);
		MOD_ADD32(_bjmodn28[ithread], j, n, _bjmodn29[ithread]);
		MOD_ADD32(_bjmodn29[ithread], j, n, _bjmodn30[ithread]);
		MOD_ADD32(_bjmodn30[ithread], j, n, _bjmodn31[ithread]);
		MOD_ADD32(_bjmodn31[ithread], j, n, _bjmodn32[ithread]);
		MOD_ADD32(_bjmodn32[ithread], j, n, _bjmodn33[ithread]);
		MOD_ADD32(_bjmodn33[ithread], j, n, _bjmodn34[ithread]);
		MOD_ADD32(_bjmodn34[ithread], j, n, _bjmodn35[ithread]);
		MOD_ADD32(_bjmodn35[ithread], j, n, _bjmodn36[ithread]);
		MOD_ADD32(_bjmodn36[ithread], j, n, _bjmodn37[ithread]);
		MOD_ADD32(_bjmodn37[ithread], j, n, _bjmodn38[ithread]);
		MOD_ADD32(_bjmodn38[ithread], j, n, _bjmodn39[ithread]);
		MOD_ADD32(_bjmodn39[ithread], j, n, _bjmodn40[ithread]);
		MOD_ADD32(_bjmodn40[ithread], j, n, _bjmodn41[ithread]);
		MOD_ADD32(_bjmodn41[ithread], j, n, _bjmodn42[ithread]);
		MOD_ADD32(_bjmodn42[ithread], j, n, _bjmodn43[ithread]);
		MOD_ADD32(_bjmodn43[ithread], j, n, _bjmodn44[ithread]);
		MOD_ADD32(_bjmodn44[ithread], j, n, _bjmodn45[ithread]);
		MOD_ADD32(_bjmodn45[ithread], j, n, _bjmodn46[ithread]);
		MOD_ADD32(_bjmodn46[ithread], j, n, _bjmodn47[ithread]);
		MOD_ADD32(_bjmodn47[ithread], j, n, _bjmodn48[ithread]);
		MOD_ADD32(_bjmodn48[ithread], j, n, _bjmodn49[ithread]);
		MOD_ADD32(_bjmodn49[ithread], j, n, _bjmodn50[ithread]);
		MOD_ADD32(_bjmodn50[ithread], j, n, _bjmodn51[ithread]);
		MOD_ADD32(_bjmodn51[ithread], j, n, _bjmodn52[ithread]);
		MOD_ADD32(_bjmodn52[ithread], j, n, _bjmodn53[ithread]);
		MOD_ADD32(_bjmodn53[ithread], j, n, _bjmodn54[ithread]);
		MOD_ADD32(_bjmodn54[ithread], j, n, _bjmodn55[ithread]);
		MOD_ADD32(_bjmodn55[ithread], j, n, _bjmodn56[ithread]);
		MOD_ADD32(_bjmodn56[ithread], j, n, _bjmodn57[ithread]);
		MOD_ADD32(_bjmodn57[ithread], j, n, _bjmodn58[ithread]);
		MOD_ADD32(_bjmodn58[ithread], j, n, _bjmodn59[ithread]);

		// Every (odd_radix)th bjmodn initializer needs to be forced-to-bigword in fermat-mod DWT case:
		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		{
			/* Start this value off at N in Fermat-mod case, so (bjmodn >= sw) check in
			fermat_carry_norm_errcheck (cf. carry.h) yields a bigword (i == 1) for j= 0:
			*/
			_bjmodn00[ithread] = n;
			_bjmodn15[ithread] = n;
			_bjmodn30[ithread] = n;
			_bjmodn45[ithread] = n;
		}
	}

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		/* Find the circular-index-shift (cf. the head-of-file comments of radix28_ditN_cy_dif1.c) by searching bjmodn01 ... bjmodn[nwt] for the one == bw: */
		if( _bjmodn01[0] == bw ) { wts_idx_incr = 0x1; };
		if( _bjmodn02[0] == bw ) { wts_idx_incr = 0x2; };
		if( _bjmodn03[0] == bw ) { wts_idx_incr = 0x3; };
		if( _bjmodn04[0] == bw ) { wts_idx_incr = 0x4; };
		if( _bjmodn05[0] == bw ) { wts_idx_incr = 0x5; };
		if( _bjmodn06[0] == bw ) { wts_idx_incr = 0x6; };
		if( _bjmodn07[0] == bw ) { wts_idx_incr = 0x7; };
		if( _bjmodn08[0] == bw ) { wts_idx_incr = 0x8; };
		if( _bjmodn09[0] == bw ) { wts_idx_incr = 0x9; };
		if( _bjmodn10[0] == bw ) { wts_idx_incr = 0xa; };
		if( _bjmodn11[0] == bw ) { wts_idx_incr = 0xb; };
		if( _bjmodn12[0] == bw ) { wts_idx_incr = 0xc; };
		if( _bjmodn13[0] == bw ) { wts_idx_incr = 0xd; };
		if( _bjmodn14[0] == bw ) { wts_idx_incr = 0xe; };

		ASSERT(HERE, wts_idx_incr != 0, "wts_idx_incr init failed!");

	#ifdef USE_SSE2
		wts_idx_inc2 = wts_idx_incr << (2*l2_sz_vd - 3);	/* In the SIMD version, use icycle0-6 as actual address
						offsets, so wts_idx_incr includes a *sizeof(vec_dbl) for the array-of-vector-doubles indexing, and another
						doubling|quadrupling|... to reflect the fact that the SIMD version of the loop is equivalent to 2|4|... scalar
						loop executions, i.e. corresponds to [#doubles in each vec_dbl] scalar-code increments of the icycle indices. */
		wts_idx_inc2 %= nwt16;	/* Need an extra mod since [2|4|...]*wts_idx_incr may be >= nwt */
	#endif
		/* Subtract nwt from the increments to ease fast-mod */
		wts_idx_incr -= nwt;
	#ifdef USE_SSE2
		wts_idx_inc2 -= nwt16;
	#endif

		/* Need this both in scalar mode and to ease the SSE2-array init */
		j = _bjmodn00[0] > sw;		bs_arr[0x0] = base[j];	bsinv_arr[0x0] = baseinv[j];
		j = _bjmodn01[0] > sw;		bs_arr[0x1] = base[j];	bsinv_arr[0x1] = baseinv[j];
		j = _bjmodn02[0] > sw;		bs_arr[0x2] = base[j];	bsinv_arr[0x2] = baseinv[j];
		j = _bjmodn03[0] > sw;		bs_arr[0x3] = base[j];	bsinv_arr[0x3] = baseinv[j];
		j = _bjmodn04[0] > sw;		bs_arr[0x4] = base[j];	bsinv_arr[0x4] = baseinv[j];
		j = _bjmodn05[0] > sw;		bs_arr[0x5] = base[j];	bsinv_arr[0x5] = baseinv[j];
		j = _bjmodn06[0] > sw;		bs_arr[0x6] = base[j];	bsinv_arr[0x6] = baseinv[j];
		j = _bjmodn07[0] > sw;		bs_arr[0x7] = base[j];	bsinv_arr[0x7] = baseinv[j];
		j = _bjmodn08[0] > sw;		bs_arr[0x8] = base[j];	bsinv_arr[0x8] = baseinv[j];
		j = _bjmodn09[0] > sw;		bs_arr[0x9] = base[j];	bsinv_arr[0x9] = baseinv[j];
		j = _bjmodn10[0] > sw;		bs_arr[0xa] = base[j];	bsinv_arr[0xa] = baseinv[j];
		j = _bjmodn11[0] > sw;		bs_arr[0xb] = base[j];	bsinv_arr[0xb] = baseinv[j];
		j = _bjmodn12[0] > sw;		bs_arr[0xc] = base[j];	bsinv_arr[0xc] = baseinv[j];
		j = _bjmodn13[0] > sw;		bs_arr[0xd] = base[j];	bsinv_arr[0xd] = baseinv[j];
		j = _bjmodn14[0] > sw;		bs_arr[0xe] = base[j];	bsinv_arr[0xe] = baseinv[j];

		/* Give icycle indices their proper starting values: */
		icycle00 = 0x0;
		icycle01 = 0x1;
		icycle02 = 0x2;
		icycle03 = 0x3;
		icycle04 = 0x4;
		icycle05 = 0x5;
		icycle06 = 0x6;
		icycle07 = 0x7;
		icycle08 = 0x8;
		icycle09 = 0x9;
		icycle10 = 0xa;
		icycle11 = 0xb;
		icycle12 = 0xc;
		icycle13 = 0xd;
		icycle14 = 0xe;

		/* Need this both in scalar mode and to ease the SSE2-array init */
		wt_arr[0x0] = wt0[ii00];	wtinv_arr[0x0] = scale*wt1[ii00];
		wt_arr[0x1] = wt0[ii01];	wtinv_arr[0x1] = scale*wt1[ii01];
		wt_arr[0x2] = wt0[ii02];	wtinv_arr[0x2] = scale*wt1[ii02];
		wt_arr[0x3] = wt0[ii03];	wtinv_arr[0x3] = scale*wt1[ii03];
		wt_arr[0x4] = wt0[ii04];	wtinv_arr[0x4] = scale*wt1[ii04];
		wt_arr[0x5] = wt0[ii05];	wtinv_arr[0x5] = scale*wt1[ii05];
		wt_arr[0x6] = wt0[ii06];	wtinv_arr[0x6] = scale*wt1[ii06];
		wt_arr[0x7] = wt0[ii07];	wtinv_arr[0x7] = scale*wt1[ii07];
		wt_arr[0x8] = wt0[ii08];	wtinv_arr[0x8] = scale*wt1[ii08];
		wt_arr[0x9] = wt0[ii09];	wtinv_arr[0x9] = scale*wt1[ii09];
		wt_arr[0xa] = wt0[ii10];	wtinv_arr[0xa] = scale*wt1[ii10];
		wt_arr[0xb] = wt0[ii11];	wtinv_arr[0xb] = scale*wt1[ii11];
		wt_arr[0xc] = wt0[ii12];	wtinv_arr[0xc] = scale*wt1[ii12];
		wt_arr[0xd] = wt0[ii13];	wtinv_arr[0xd] = scale*wt1[ii13];
		wt_arr[0xe] = wt0[ii14];	wtinv_arr[0xe] = scale*wt1[ii14];

	#ifdef USE_SSE2

		tmp = half_arr;
		tmp->d0 = wt_arr[icycle00];	++tmp;
		tmp->d0 = wt_arr[icycle01];	++tmp;
		tmp->d0 = wt_arr[icycle02];	++tmp;
		tmp->d0 = wt_arr[icycle03];	++tmp;
		tmp->d0 = wt_arr[icycle04];	++tmp;
		tmp->d0 = wt_arr[icycle05];	++tmp;
		tmp->d0 = wt_arr[icycle06];	++tmp;
		tmp->d0 = wt_arr[icycle07];	++tmp;
		tmp->d0 = wt_arr[icycle08];	++tmp;
		tmp->d0 = wt_arr[icycle09];	++tmp;
		tmp->d0 = wt_arr[icycle10];	++tmp;
		tmp->d0 = wt_arr[icycle11];	++tmp;
		tmp->d0 = wt_arr[icycle12];	++tmp;
		tmp->d0 = wt_arr[icycle13];	++tmp;
		tmp->d0 = wt_arr[icycle14];	++tmp;

		tmp->d0 = wtinv_arr[icycle00];	++tmp;
		tmp->d0 = wtinv_arr[icycle01];	++tmp;
		tmp->d0 = wtinv_arr[icycle02];	++tmp;
		tmp->d0 = wtinv_arr[icycle03];	++tmp;
		tmp->d0 = wtinv_arr[icycle04];	++tmp;
		tmp->d0 = wtinv_arr[icycle05];	++tmp;
		tmp->d0 = wtinv_arr[icycle06];	++tmp;
		tmp->d0 = wtinv_arr[icycle07];	++tmp;
		tmp->d0 = wtinv_arr[icycle08];	++tmp;
		tmp->d0 = wtinv_arr[icycle09];	++tmp;
		tmp->d0 = wtinv_arr[icycle10];	++tmp;
		tmp->d0 = wtinv_arr[icycle11];	++tmp;
		tmp->d0 = wtinv_arr[icycle12];	++tmp;
		tmp->d0 = wtinv_arr[icycle13];	++tmp;
		tmp->d0 = wtinv_arr[icycle14];	++tmp;

		/* Now set the imaginary parts to the values corresponding to the 2nd of each pair of scalar-mode loop passes.
		Use this sequence for mod-add, as it is faster than general-mod '% nwt' */
		jcycle00 = icycle00 + wts_idx_incr;		jcycle00 += ( (-(jcycle00 < 0)) & nwt);
		jcycle01 = icycle01 + wts_idx_incr;		jcycle01 += ( (-(jcycle01 < 0)) & nwt);
		jcycle02 = icycle02 + wts_idx_incr;		jcycle02 += ( (-(jcycle02 < 0)) & nwt);
		jcycle03 = icycle03 + wts_idx_incr;		jcycle03 += ( (-(jcycle03 < 0)) & nwt);
		jcycle04 = icycle04 + wts_idx_incr;		jcycle04 += ( (-(jcycle04 < 0)) & nwt);
		jcycle05 = icycle05 + wts_idx_incr;		jcycle05 += ( (-(jcycle05 < 0)) & nwt);
		jcycle06 = icycle06 + wts_idx_incr;		jcycle06 += ( (-(jcycle06 < 0)) & nwt);
		jcycle07 = icycle07 + wts_idx_incr;		jcycle07 += ( (-(jcycle07 < 0)) & nwt);
		jcycle08 = icycle08 + wts_idx_incr;		jcycle08 += ( (-(jcycle08 < 0)) & nwt);
		jcycle09 = icycle09 + wts_idx_incr;		jcycle09 += ( (-(jcycle09 < 0)) & nwt);
		jcycle10 = icycle10 + wts_idx_incr;		jcycle10 += ( (-(jcycle10 < 0)) & nwt);
		jcycle11 = icycle11 + wts_idx_incr;		jcycle11 += ( (-(jcycle11 < 0)) & nwt);
		jcycle12 = icycle12 + wts_idx_incr;		jcycle12 += ( (-(jcycle12 < 0)) & nwt);
		jcycle13 = icycle13 + wts_idx_incr;		jcycle13 += ( (-(jcycle13 < 0)) & nwt);
		jcycle14 = icycle14 + wts_idx_incr;		jcycle14 += ( (-(jcycle14 < 0)) & nwt);

		tmp = half_arr;
		tmp->d1 = wt_arr[jcycle00];	++tmp;
		tmp->d1 = wt_arr[jcycle01];	++tmp;
		tmp->d1 = wt_arr[jcycle02];	++tmp;
		tmp->d1 = wt_arr[jcycle03];	++tmp;
		tmp->d1 = wt_arr[jcycle04];	++tmp;
		tmp->d1 = wt_arr[jcycle05];	++tmp;
		tmp->d1 = wt_arr[jcycle06];	++tmp;
		tmp->d1 = wt_arr[jcycle07];	++tmp;
		tmp->d1 = wt_arr[jcycle08];	++tmp;
		tmp->d1 = wt_arr[jcycle09];	++tmp;
		tmp->d1 = wt_arr[jcycle10];	++tmp;
		tmp->d1 = wt_arr[jcycle11];	++tmp;
		tmp->d1 = wt_arr[jcycle12];	++tmp;
		tmp->d1 = wt_arr[jcycle13];	++tmp;
		tmp->d1 = wt_arr[jcycle14];	++tmp;

		tmp->d1 = wtinv_arr[jcycle00];	++tmp;
		tmp->d1 = wtinv_arr[jcycle01];	++tmp;
		tmp->d1 = wtinv_arr[jcycle02];	++tmp;
		tmp->d1 = wtinv_arr[jcycle03];	++tmp;
		tmp->d1 = wtinv_arr[jcycle04];	++tmp;
		tmp->d1 = wtinv_arr[jcycle05];	++tmp;
		tmp->d1 = wtinv_arr[jcycle06];	++tmp;
		tmp->d1 = wtinv_arr[jcycle07];	++tmp;
		tmp->d1 = wtinv_arr[jcycle08];	++tmp;
		tmp->d1 = wtinv_arr[jcycle09];	++tmp;
		tmp->d1 = wtinv_arr[jcycle10];	++tmp;
		tmp->d1 = wtinv_arr[jcycle11];	++tmp;
		tmp->d1 = wtinv_arr[jcycle12];	++tmp;
		tmp->d1 = wtinv_arr[jcycle13];	++tmp;
		tmp->d1 = wtinv_arr[jcycle14];	++tmp;

	#ifdef USE_AVX
		kcycle00 = jcycle00 + wts_idx_incr;		kcycle00 += ( (-(kcycle00 < 0)) & nwt);
		kcycle01 = jcycle01 + wts_idx_incr;		kcycle01 += ( (-(kcycle01 < 0)) & nwt);
		kcycle02 = jcycle02 + wts_idx_incr;		kcycle02 += ( (-(kcycle02 < 0)) & nwt);
		kcycle03 = jcycle03 + wts_idx_incr;		kcycle03 += ( (-(kcycle03 < 0)) & nwt);
		kcycle04 = jcycle04 + wts_idx_incr;		kcycle04 += ( (-(kcycle04 < 0)) & nwt);
		kcycle05 = jcycle05 + wts_idx_incr;		kcycle05 += ( (-(kcycle05 < 0)) & nwt);
		kcycle06 = jcycle06 + wts_idx_incr;		kcycle06 += ( (-(kcycle06 < 0)) & nwt);
		kcycle07 = jcycle07 + wts_idx_incr;		kcycle07 += ( (-(kcycle07 < 0)) & nwt);
		kcycle08 = jcycle08 + wts_idx_incr;		kcycle08 += ( (-(kcycle08 < 0)) & nwt);
		kcycle09 = jcycle09 + wts_idx_incr;		kcycle09 += ( (-(kcycle09 < 0)) & nwt);
		kcycle10 = jcycle10 + wts_idx_incr;		kcycle10 += ( (-(kcycle10 < 0)) & nwt);
		kcycle11 = jcycle11 + wts_idx_incr;		kcycle11 += ( (-(kcycle11 < 0)) & nwt);
		kcycle12 = jcycle12 + wts_idx_incr;		kcycle12 += ( (-(kcycle12 < 0)) & nwt);
		kcycle13 = jcycle13 + wts_idx_incr;		kcycle13 += ( (-(kcycle13 < 0)) & nwt);
		kcycle14 = jcycle14 + wts_idx_incr;		kcycle14 += ( (-(kcycle14 < 0)) & nwt);

		tmp = half_arr;
		tmp->d2 = wt_arr[kcycle00];	++tmp;
		tmp->d2 = wt_arr[kcycle01];	++tmp;
		tmp->d2 = wt_arr[kcycle02];	++tmp;
		tmp->d2 = wt_arr[kcycle03];	++tmp;
		tmp->d2 = wt_arr[kcycle04];	++tmp;
		tmp->d2 = wt_arr[kcycle05];	++tmp;
		tmp->d2 = wt_arr[kcycle06];	++tmp;
		tmp->d2 = wt_arr[kcycle07];	++tmp;
		tmp->d2 = wt_arr[kcycle08];	++tmp;
		tmp->d2 = wt_arr[kcycle09];	++tmp;
		tmp->d2 = wt_arr[kcycle10];	++tmp;
		tmp->d2 = wt_arr[kcycle11];	++tmp;
		tmp->d2 = wt_arr[kcycle12];	++tmp;
		tmp->d2 = wt_arr[kcycle13];	++tmp;
		tmp->d2 = wt_arr[kcycle14];	++tmp;

		tmp->d2 = wtinv_arr[kcycle00];	++tmp;
		tmp->d2 = wtinv_arr[kcycle01];	++tmp;
		tmp->d2 = wtinv_arr[kcycle02];	++tmp;
		tmp->d2 = wtinv_arr[kcycle03];	++tmp;
		tmp->d2 = wtinv_arr[kcycle04];	++tmp;
		tmp->d2 = wtinv_arr[kcycle05];	++tmp;
		tmp->d2 = wtinv_arr[kcycle06];	++tmp;
		tmp->d2 = wtinv_arr[kcycle07];	++tmp;
		tmp->d2 = wtinv_arr[kcycle08];	++tmp;
		tmp->d2 = wtinv_arr[kcycle09];	++tmp;
		tmp->d2 = wtinv_arr[kcycle10];	++tmp;
		tmp->d2 = wtinv_arr[kcycle11];	++tmp;
		tmp->d2 = wtinv_arr[kcycle12];	++tmp;
		tmp->d2 = wtinv_arr[kcycle13];	++tmp;
		tmp->d2 = wtinv_arr[kcycle14];	++tmp;

		lcycle00 = kcycle00 + wts_idx_incr;		lcycle00 += ( (-(lcycle00 < 0)) & nwt);
		lcycle01 = kcycle01 + wts_idx_incr;		lcycle01 += ( (-(lcycle01 < 0)) & nwt);
		lcycle02 = kcycle02 + wts_idx_incr;		lcycle02 += ( (-(lcycle02 < 0)) & nwt);
		lcycle03 = kcycle03 + wts_idx_incr;		lcycle03 += ( (-(lcycle03 < 0)) & nwt);
		lcycle04 = kcycle04 + wts_idx_incr;		lcycle04 += ( (-(lcycle04 < 0)) & nwt);
		lcycle05 = kcycle05 + wts_idx_incr;		lcycle05 += ( (-(lcycle05 < 0)) & nwt);
		lcycle06 = kcycle06 + wts_idx_incr;		lcycle06 += ( (-(lcycle06 < 0)) & nwt);
		lcycle07 = kcycle07 + wts_idx_incr;		lcycle07 += ( (-(lcycle07 < 0)) & nwt);
		lcycle08 = kcycle08 + wts_idx_incr;		lcycle08 += ( (-(lcycle08 < 0)) & nwt);
		lcycle09 = kcycle09 + wts_idx_incr;		lcycle09 += ( (-(lcycle09 < 0)) & nwt);
		lcycle10 = kcycle10 + wts_idx_incr;		lcycle10 += ( (-(lcycle10 < 0)) & nwt);
		lcycle11 = kcycle11 + wts_idx_incr;		lcycle11 += ( (-(lcycle11 < 0)) & nwt);
		lcycle12 = kcycle12 + wts_idx_incr;		lcycle12 += ( (-(lcycle12 < 0)) & nwt);
		lcycle13 = kcycle13 + wts_idx_incr;		lcycle13 += ( (-(lcycle13 < 0)) & nwt);
		lcycle14 = kcycle14 + wts_idx_incr;		lcycle14 += ( (-(lcycle14 < 0)) & nwt);

		tmp = half_arr;
		tmp->d3 = wt_arr[lcycle00];	++tmp;
		tmp->d3 = wt_arr[lcycle01];	++tmp;
		tmp->d3 = wt_arr[lcycle02];	++tmp;
		tmp->d3 = wt_arr[lcycle03];	++tmp;
		tmp->d3 = wt_arr[lcycle04];	++tmp;
		tmp->d3 = wt_arr[lcycle05];	++tmp;
		tmp->d3 = wt_arr[lcycle06];	++tmp;
		tmp->d3 = wt_arr[lcycle07];	++tmp;
		tmp->d3 = wt_arr[lcycle08];	++tmp;
		tmp->d3 = wt_arr[lcycle09];	++tmp;
		tmp->d3 = wt_arr[lcycle10];	++tmp;
		tmp->d3 = wt_arr[lcycle11];	++tmp;
		tmp->d3 = wt_arr[lcycle12];	++tmp;
		tmp->d3 = wt_arr[lcycle13];	++tmp;
		tmp->d3 = wt_arr[lcycle14];	++tmp;

		tmp->d3 = wtinv_arr[lcycle00];	++tmp;
		tmp->d3 = wtinv_arr[lcycle01];	++tmp;
		tmp->d3 = wtinv_arr[lcycle02];	++tmp;
		tmp->d3 = wtinv_arr[lcycle03];	++tmp;
		tmp->d3 = wtinv_arr[lcycle04];	++tmp;
		tmp->d3 = wtinv_arr[lcycle05];	++tmp;
		tmp->d3 = wtinv_arr[lcycle06];	++tmp;
		tmp->d3 = wtinv_arr[lcycle07];	++tmp;
		tmp->d3 = wtinv_arr[lcycle08];	++tmp;
		tmp->d3 = wtinv_arr[lcycle09];	++tmp;
		tmp->d3 = wtinv_arr[lcycle10];	++tmp;
		tmp->d3 = wtinv_arr[lcycle11];	++tmp;
		tmp->d3 = wtinv_arr[lcycle12];	++tmp;
		tmp->d3 = wtinv_arr[lcycle13];	++tmp;
		tmp->d3 = wtinv_arr[lcycle14];	++tmp;
	#endif

		tmp = half_arr + odd_radix*2;	/* Put the base-mini-arrays right after the weights */

	  #ifdef USE_AVX

		// Each transposed-data quartet in the AVX carry macro needs linearly incrementing bs_arr data (mod odd_radix);
		// Need all [odd_radix] possible such length-4 index subsequences, which will be accessed via their head element
		// by the [ijkl]cycle* index quartets in the respective carry-macro call:
		tmp->d0 = bs_arr[0x0];	tmp->d1 = bs_arr[0x1];	tmp->d2 = bs_arr[0x2];	tmp->d3 = bs_arr[0x3];	++tmp;
		tmp->d0 = bs_arr[0x1];	tmp->d1 = bs_arr[0x2];	tmp->d2 = bs_arr[0x3];	tmp->d3 = bs_arr[0x4];	++tmp;
		tmp->d0 = bs_arr[0x2];	tmp->d1 = bs_arr[0x3];	tmp->d2 = bs_arr[0x4];	tmp->d3 = bs_arr[0x5];	++tmp;
		tmp->d0 = bs_arr[0x3];	tmp->d1 = bs_arr[0x4];	tmp->d2 = bs_arr[0x5];	tmp->d3 = bs_arr[0x6];	++tmp;
		tmp->d0 = bs_arr[0x4];	tmp->d1 = bs_arr[0x5];	tmp->d2 = bs_arr[0x6];	tmp->d3 = bs_arr[0x7];	++tmp;
		tmp->d0 = bs_arr[0x5];	tmp->d1 = bs_arr[0x6];	tmp->d2 = bs_arr[0x7];	tmp->d3 = bs_arr[0x8];	++tmp;
		tmp->d0 = bs_arr[0x6];	tmp->d1 = bs_arr[0x7];	tmp->d2 = bs_arr[0x8];	tmp->d3 = bs_arr[0x9];	++tmp;
		tmp->d0 = bs_arr[0x7];	tmp->d1 = bs_arr[0x8];	tmp->d2 = bs_arr[0x9];	tmp->d3 = bs_arr[0xa];	++tmp;
		tmp->d0 = bs_arr[0x8];	tmp->d1 = bs_arr[0x9];	tmp->d2 = bs_arr[0xa];	tmp->d3 = bs_arr[0xb];	++tmp;
		tmp->d0 = bs_arr[0x9];	tmp->d1 = bs_arr[0xa];	tmp->d2 = bs_arr[0xb];	tmp->d3 = bs_arr[0xc];	++tmp;
		tmp->d0 = bs_arr[0xa];	tmp->d1 = bs_arr[0xb];	tmp->d2 = bs_arr[0xc];	tmp->d3 = bs_arr[0xd];	++tmp;
		tmp->d0 = bs_arr[0xb];	tmp->d1 = bs_arr[0xc];	tmp->d2 = bs_arr[0xd];	tmp->d3 = bs_arr[0xe];	++tmp;
		tmp->d0 = bs_arr[0xc];	tmp->d1 = bs_arr[0xd];	tmp->d2 = bs_arr[0xe];	tmp->d3 = bs_arr[0x0];	++tmp;
		tmp->d0 = bs_arr[0xd];	tmp->d1 = bs_arr[0xe];	tmp->d2 = bs_arr[0x0];	tmp->d3 = bs_arr[0x1];	++tmp;
		tmp->d0 = bs_arr[0xe];	tmp->d1 = bs_arr[0x0];	tmp->d2 = bs_arr[0x1];	tmp->d3 = bs_arr[0x2];	++tmp;

		tmp->d0 = bsinv_arr[0x0];	tmp->d1 = bsinv_arr[0x1];	tmp->d2 = bsinv_arr[0x2];	tmp->d3 = bsinv_arr[0x3];	++tmp;
		tmp->d0 = bsinv_arr[0x1];	tmp->d1 = bsinv_arr[0x2];	tmp->d2 = bsinv_arr[0x3];	tmp->d3 = bsinv_arr[0x4];	++tmp;
		tmp->d0 = bsinv_arr[0x2];	tmp->d1 = bsinv_arr[0x3];	tmp->d2 = bsinv_arr[0x4];	tmp->d3 = bsinv_arr[0x5];	++tmp;
		tmp->d0 = bsinv_arr[0x3];	tmp->d1 = bsinv_arr[0x4];	tmp->d2 = bsinv_arr[0x5];	tmp->d3 = bsinv_arr[0x6];	++tmp;
		tmp->d0 = bsinv_arr[0x4];	tmp->d1 = bsinv_arr[0x5];	tmp->d2 = bsinv_arr[0x6];	tmp->d3 = bsinv_arr[0x7];	++tmp;
		tmp->d0 = bsinv_arr[0x5];	tmp->d1 = bsinv_arr[0x6];	tmp->d2 = bsinv_arr[0x7];	tmp->d3 = bsinv_arr[0x8];	++tmp;
		tmp->d0 = bsinv_arr[0x6];	tmp->d1 = bsinv_arr[0x7];	tmp->d2 = bsinv_arr[0x8];	tmp->d3 = bsinv_arr[0x9];	++tmp;
		tmp->d0 = bsinv_arr[0x7];	tmp->d1 = bsinv_arr[0x8];	tmp->d2 = bsinv_arr[0x9];	tmp->d3 = bsinv_arr[0xa];	++tmp;
		tmp->d0 = bsinv_arr[0x8];	tmp->d1 = bsinv_arr[0x9];	tmp->d2 = bsinv_arr[0xa];	tmp->d3 = bsinv_arr[0xb];	++tmp;
		tmp->d0 = bsinv_arr[0x9];	tmp->d1 = bsinv_arr[0xa];	tmp->d2 = bsinv_arr[0xb];	tmp->d3 = bsinv_arr[0xc];	++tmp;
		tmp->d0 = bsinv_arr[0xa];	tmp->d1 = bsinv_arr[0xb];	tmp->d2 = bsinv_arr[0xc];	tmp->d3 = bsinv_arr[0xd];	++tmp;
		tmp->d0 = bsinv_arr[0xb];	tmp->d1 = bsinv_arr[0xc];	tmp->d2 = bsinv_arr[0xd];	tmp->d3 = bsinv_arr[0xe];	++tmp;
		tmp->d0 = bsinv_arr[0xc];	tmp->d1 = bsinv_arr[0xd];	tmp->d2 = bsinv_arr[0xe];	tmp->d3 = bsinv_arr[0x0];	++tmp;
		tmp->d0 = bsinv_arr[0xd];	tmp->d1 = bsinv_arr[0xe];	tmp->d2 = bsinv_arr[0x0];	tmp->d3 = bsinv_arr[0x1];	++tmp;
		tmp->d0 = bsinv_arr[0xe];	tmp->d1 = bsinv_arr[0x0];	tmp->d2 = bsinv_arr[0x1];	tmp->d3 = bsinv_arr[0x2];	++tmp;

	  #else

		/* In SSE2 mode, because we apply doubled weights to data arranged as [a.re,b.re,...],[a.im,b.im,...] but apply
		doubled base multipliers to shuffled data [a.re,a.im],[b.re,b.im],... (i.e. shuffled to yield same data layout as
		in the scalar case), the weights need to have disparate real and imag parts, but the base/baseinv terms do not: */
		VEC_DBL_INIT(tmp, bs_arr[0x0]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[0x1]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[0x2]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[0x3]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[0x4]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[0x5]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[0x6]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[0x7]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[0x8]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[0x9]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[0xa]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[0xb]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[0xc]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[0xd]);	++tmp;
		VEC_DBL_INIT(tmp, bs_arr[0xe]);	++tmp;

		VEC_DBL_INIT(tmp, bsinv_arr[0x0]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[0x1]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[0x2]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[0x3]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[0x4]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[0x5]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[0x6]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[0x7]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[0x8]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[0x9]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[0xa]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[0xb]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[0xc]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[0xd]);	++tmp;
		VEC_DBL_INIT(tmp, bsinv_arr[0xe]);	++tmp;

	  #endif

		// Propagate the above consts to the remaining threads:
		nbytes = RADIX*sz_vd;
		tmp = half_arr;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}

		icycle00 <<= l2_sz_vd;		jcycle00 <<= l2_sz_vd;
		icycle01 <<= l2_sz_vd;		jcycle01 <<= l2_sz_vd;
		icycle02 <<= l2_sz_vd;		jcycle02 <<= l2_sz_vd;
		icycle03 <<= l2_sz_vd;		jcycle03 <<= l2_sz_vd;
		icycle04 <<= l2_sz_vd;		jcycle04 <<= l2_sz_vd;
		icycle05 <<= l2_sz_vd;		jcycle05 <<= l2_sz_vd;
		icycle06 <<= l2_sz_vd;		jcycle06 <<= l2_sz_vd;
		icycle07 <<= l2_sz_vd;		jcycle07 <<= l2_sz_vd;
		icycle08 <<= l2_sz_vd;		jcycle08 <<= l2_sz_vd;
		icycle09 <<= l2_sz_vd;		jcycle09 <<= l2_sz_vd;
		icycle10 <<= l2_sz_vd;		jcycle10 <<= l2_sz_vd;
		icycle11 <<= l2_sz_vd;		jcycle11 <<= l2_sz_vd;
		icycle12 <<= l2_sz_vd;		jcycle12 <<= l2_sz_vd;
		icycle13 <<= l2_sz_vd;		jcycle13 <<= l2_sz_vd;
		icycle14 <<= l2_sz_vd;		jcycle14 <<= l2_sz_vd;

	  #ifdef USE_AVX
		kcycle00 <<= l2_sz_vd;		lcycle00 <<= l2_sz_vd;
		kcycle01 <<= l2_sz_vd;		lcycle01 <<= l2_sz_vd;
		kcycle02 <<= l2_sz_vd;		lcycle02 <<= l2_sz_vd;
		kcycle03 <<= l2_sz_vd;		lcycle03 <<= l2_sz_vd;
		kcycle04 <<= l2_sz_vd;		lcycle04 <<= l2_sz_vd;
		kcycle05 <<= l2_sz_vd;		lcycle05 <<= l2_sz_vd;
		kcycle06 <<= l2_sz_vd;		lcycle06 <<= l2_sz_vd;
		kcycle07 <<= l2_sz_vd;		lcycle07 <<= l2_sz_vd;
		kcycle08 <<= l2_sz_vd;		lcycle08 <<= l2_sz_vd;
		kcycle09 <<= l2_sz_vd;		lcycle09 <<= l2_sz_vd;
		kcycle10 <<= l2_sz_vd;		lcycle10 <<= l2_sz_vd;
		kcycle11 <<= l2_sz_vd;		lcycle11 <<= l2_sz_vd;
		kcycle12 <<= l2_sz_vd;		lcycle12 <<= l2_sz_vd;
		kcycle13 <<= l2_sz_vd;		lcycle13 <<= l2_sz_vd;
		kcycle14 <<= l2_sz_vd;		lcycle14 <<= l2_sz_vd;
	  #endif

	#endif	/* USE_SSE2 */
	}
/*	fprintf(stderr, "%s: wts_idx_incr = %d\n", func,wts_idx_incr);*/

#ifdef USE_PTHREAD

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart = _jstart[ithread];
		jhi    = _jhi[ithread];
	}
	// These inits must occur just once, in full-pass mode,
	// in order to get the array-index-offset values of the icycle/jcycle indices right:
	else if(full_pass)	// Fermat-mod
	{
		tdat[0].icycle00 = icycle00;
		tdat[0].icycle01 = icycle01;
		tdat[0].icycle02 = icycle02;
		tdat[0].icycle03 = icycle03;
		tdat[0].icycle04 = icycle04;
		tdat[0].icycle05 = icycle05;
		tdat[0].icycle06 = icycle06;
		tdat[0].icycle07 = icycle07;
		tdat[0].icycle08 = icycle08;
		tdat[0].icycle09 = icycle09;
		tdat[0].icycle10 = icycle10;
		tdat[0].icycle11 = icycle11;
		tdat[0].icycle12 = icycle12;
		tdat[0].icycle13 = icycle13;
		tdat[0].icycle14 = icycle14;

	#ifdef USE_SSE2
		tdat[0].wts_idx_inc2 = wts_idx_inc2;
		tdat[0].jcycle00 = jcycle00;
		tdat[0].jcycle01 = jcycle01;
		tdat[0].jcycle02 = jcycle02;
		tdat[0].jcycle03 = jcycle03;
		tdat[0].jcycle04 = jcycle04;
		tdat[0].jcycle05 = jcycle05;
		tdat[0].jcycle06 = jcycle06;
		tdat[0].jcycle07 = jcycle07;
		tdat[0].jcycle08 = jcycle08;
		tdat[0].jcycle09 = jcycle09;
		tdat[0].jcycle10 = jcycle10;
		tdat[0].jcycle11 = jcycle11;
		tdat[0].jcycle12 = jcycle12;
		tdat[0].jcycle13 = jcycle13;
		tdat[0].jcycle14 = jcycle14;

	  #ifdef USE_AVX
		tdat[0].kcycle00 = kcycle00;
		tdat[0].kcycle01 = kcycle01;
		tdat[0].kcycle02 = kcycle02;
		tdat[0].kcycle03 = kcycle03;
		tdat[0].kcycle04 = kcycle04;
		tdat[0].kcycle05 = kcycle05;
		tdat[0].kcycle06 = kcycle06;
		tdat[0].kcycle07 = kcycle07;
		tdat[0].kcycle08 = kcycle08;
		tdat[0].kcycle09 = kcycle09;
		tdat[0].kcycle10 = kcycle10;
		tdat[0].kcycle11 = kcycle11;
		tdat[0].kcycle12 = kcycle12;
		tdat[0].kcycle13 = kcycle13;
		tdat[0].kcycle14 = kcycle14;

		tdat[0].lcycle00 = lcycle00;
		tdat[0].lcycle01 = lcycle01;
		tdat[0].lcycle02 = lcycle02;
		tdat[0].lcycle03 = lcycle03;
		tdat[0].lcycle04 = lcycle04;
		tdat[0].lcycle05 = lcycle05;
		tdat[0].lcycle06 = lcycle06;
		tdat[0].lcycle07 = lcycle07;
		tdat[0].lcycle08 = lcycle08;
		tdat[0].lcycle09 = lcycle09;
		tdat[0].lcycle10 = lcycle10;
		tdat[0].lcycle11 = lcycle11;
		tdat[0].lcycle12 = lcycle12;
		tdat[0].lcycle13 = lcycle13;
		tdat[0].lcycle14 = lcycle14;
	  #endif
	#endif

		// For remaining threads, simulate the loop-evolution of the above indices:
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			jstart = _jstart[ithread];
			jhi    = _jhi[ithread];
			for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
			{
				for(j = jstart; j < jhi; j += stride)
				{
				#ifndef USE_SSE2	// Scalar-double mode uses non-pointerized icycle values:

					icycle00 += wts_idx_incr;		icycle00 += ( (-(int)((uint32)icycle00 >> 31)) & nwt);
					icycle01 += wts_idx_incr;		icycle01 += ( (-(int)((uint32)icycle01 >> 31)) & nwt);
					icycle02 += wts_idx_incr;		icycle02 += ( (-(int)((uint32)icycle02 >> 31)) & nwt);
					icycle03 += wts_idx_incr;		icycle03 += ( (-(int)((uint32)icycle03 >> 31)) & nwt);
					icycle04 += wts_idx_incr;		icycle04 += ( (-(int)((uint32)icycle04 >> 31)) & nwt);
					icycle05 += wts_idx_incr;		icycle05 += ( (-(int)((uint32)icycle05 >> 31)) & nwt);
					icycle06 += wts_idx_incr;		icycle06 += ( (-(int)((uint32)icycle06 >> 31)) & nwt);
					icycle07 += wts_idx_incr;		icycle07 += ( (-(int)((uint32)icycle07 >> 31)) & nwt);
					icycle08 += wts_idx_incr;		icycle08 += ( (-(int)((uint32)icycle08 >> 31)) & nwt);
					icycle09 += wts_idx_incr;		icycle09 += ( (-(int)((uint32)icycle09 >> 31)) & nwt);
					icycle10 += wts_idx_incr;		icycle10 += ( (-(int)((uint32)icycle10 >> 31)) & nwt);
					icycle11 += wts_idx_incr;		icycle11 += ( (-(int)((uint32)icycle11 >> 31)) & nwt);
					icycle12 += wts_idx_incr;		icycle12 += ( (-(int)((uint32)icycle12 >> 31)) & nwt);
					icycle13 += wts_idx_incr;		icycle13 += ( (-(int)((uint32)icycle13 >> 31)) & nwt);
					icycle14 += wts_idx_incr;		icycle14 += ( (-(int)((uint32)icycle14 >> 31)) & nwt);

				#else

					icycle00 += wts_idx_inc2;		icycle00 += ( (-(int)((uint32)icycle00 >> 31)) & nwt16);
					icycle01 += wts_idx_inc2;		icycle01 += ( (-(int)((uint32)icycle01 >> 31)) & nwt16);
					icycle02 += wts_idx_inc2;		icycle02 += ( (-(int)((uint32)icycle02 >> 31)) & nwt16);
					icycle03 += wts_idx_inc2;		icycle03 += ( (-(int)((uint32)icycle03 >> 31)) & nwt16);
					icycle04 += wts_idx_inc2;		icycle04 += ( (-(int)((uint32)icycle04 >> 31)) & nwt16);
					icycle05 += wts_idx_inc2;		icycle05 += ( (-(int)((uint32)icycle05 >> 31)) & nwt16);
					icycle06 += wts_idx_inc2;		icycle06 += ( (-(int)((uint32)icycle06 >> 31)) & nwt16);
					icycle07 += wts_idx_inc2;		icycle07 += ( (-(int)((uint32)icycle07 >> 31)) & nwt16);
					icycle08 += wts_idx_inc2;		icycle08 += ( (-(int)((uint32)icycle08 >> 31)) & nwt16);
					icycle09 += wts_idx_inc2;		icycle09 += ( (-(int)((uint32)icycle09 >> 31)) & nwt16);
					icycle10 += wts_idx_inc2;		icycle10 += ( (-(int)((uint32)icycle10 >> 31)) & nwt16);
					icycle11 += wts_idx_inc2;		icycle11 += ( (-(int)((uint32)icycle11 >> 31)) & nwt16);
					icycle12 += wts_idx_inc2;		icycle12 += ( (-(int)((uint32)icycle12 >> 31)) & nwt16);
					icycle13 += wts_idx_inc2;		icycle13 += ( (-(int)((uint32)icycle13 >> 31)) & nwt16);
					icycle14 += wts_idx_inc2;		icycle14 += ( (-(int)((uint32)icycle14 >> 31)) & nwt16);

					jcycle00 += wts_idx_inc2;		jcycle00 += ( (-(int)((uint32)jcycle00 >> 31)) & nwt16);
					jcycle01 += wts_idx_inc2;		jcycle01 += ( (-(int)((uint32)jcycle01 >> 31)) & nwt16);
					jcycle02 += wts_idx_inc2;		jcycle02 += ( (-(int)((uint32)jcycle02 >> 31)) & nwt16);
					jcycle03 += wts_idx_inc2;		jcycle03 += ( (-(int)((uint32)jcycle03 >> 31)) & nwt16);
					jcycle04 += wts_idx_inc2;		jcycle04 += ( (-(int)((uint32)jcycle04 >> 31)) & nwt16);
					jcycle05 += wts_idx_inc2;		jcycle05 += ( (-(int)((uint32)jcycle05 >> 31)) & nwt16);
					jcycle06 += wts_idx_inc2;		jcycle06 += ( (-(int)((uint32)jcycle06 >> 31)) & nwt16);
					jcycle07 += wts_idx_inc2;		jcycle07 += ( (-(int)((uint32)jcycle07 >> 31)) & nwt16);
					jcycle08 += wts_idx_inc2;		jcycle08 += ( (-(int)((uint32)jcycle08 >> 31)) & nwt16);
					jcycle09 += wts_idx_inc2;		jcycle09 += ( (-(int)((uint32)jcycle09 >> 31)) & nwt16);
					jcycle10 += wts_idx_inc2;		jcycle10 += ( (-(int)((uint32)jcycle10 >> 31)) & nwt16);
					jcycle11 += wts_idx_inc2;		jcycle11 += ( (-(int)((uint32)jcycle11 >> 31)) & nwt16);
					jcycle12 += wts_idx_inc2;		jcycle12 += ( (-(int)((uint32)jcycle12 >> 31)) & nwt16);
					jcycle13 += wts_idx_inc2;		jcycle13 += ( (-(int)((uint32)jcycle13 >> 31)) & nwt16);
					jcycle14 += wts_idx_inc2;		jcycle14 += ( (-(int)((uint32)jcycle14 >> 31)) & nwt16);
				  #ifdef USE_AVX
					kcycle00 += wts_idx_inc2;		kcycle00 += ( (-(int)((uint32)kcycle00 >> 31)) & nwt16);
					kcycle01 += wts_idx_inc2;		kcycle01 += ( (-(int)((uint32)kcycle01 >> 31)) & nwt16);
					kcycle02 += wts_idx_inc2;		kcycle02 += ( (-(int)((uint32)kcycle02 >> 31)) & nwt16);
					kcycle03 += wts_idx_inc2;		kcycle03 += ( (-(int)((uint32)kcycle03 >> 31)) & nwt16);
					kcycle04 += wts_idx_inc2;		kcycle04 += ( (-(int)((uint32)kcycle04 >> 31)) & nwt16);
					kcycle05 += wts_idx_inc2;		kcycle05 += ( (-(int)((uint32)kcycle05 >> 31)) & nwt16);
					kcycle06 += wts_idx_inc2;		kcycle06 += ( (-(int)((uint32)kcycle06 >> 31)) & nwt16);
					kcycle07 += wts_idx_inc2;		kcycle07 += ( (-(int)((uint32)kcycle07 >> 31)) & nwt16);
					kcycle08 += wts_idx_inc2;		kcycle08 += ( (-(int)((uint32)kcycle08 >> 31)) & nwt16);
					kcycle09 += wts_idx_inc2;		kcycle09 += ( (-(int)((uint32)kcycle09 >> 31)) & nwt16);
					kcycle10 += wts_idx_inc2;		kcycle10 += ( (-(int)((uint32)kcycle10 >> 31)) & nwt16);
					kcycle11 += wts_idx_inc2;		kcycle11 += ( (-(int)((uint32)kcycle11 >> 31)) & nwt16);
					kcycle12 += wts_idx_inc2;		kcycle12 += ( (-(int)((uint32)kcycle12 >> 31)) & nwt16);
					kcycle13 += wts_idx_inc2;		kcycle13 += ( (-(int)((uint32)kcycle13 >> 31)) & nwt16);
					kcycle14 += wts_idx_inc2;		kcycle14 += ( (-(int)((uint32)kcycle14 >> 31)) & nwt16);

					lcycle00 += wts_idx_inc2;		lcycle00 += ( (-(int)((uint32)lcycle00 >> 31)) & nwt16);
					lcycle01 += wts_idx_inc2;		lcycle01 += ( (-(int)((uint32)lcycle01 >> 31)) & nwt16);
					lcycle02 += wts_idx_inc2;		lcycle02 += ( (-(int)((uint32)lcycle02 >> 31)) & nwt16);
					lcycle03 += wts_idx_inc2;		lcycle03 += ( (-(int)((uint32)lcycle03 >> 31)) & nwt16);
					lcycle04 += wts_idx_inc2;		lcycle04 += ( (-(int)((uint32)lcycle04 >> 31)) & nwt16);
					lcycle05 += wts_idx_inc2;		lcycle05 += ( (-(int)((uint32)lcycle05 >> 31)) & nwt16);
					lcycle06 += wts_idx_inc2;		lcycle06 += ( (-(int)((uint32)lcycle06 >> 31)) & nwt16);
					lcycle07 += wts_idx_inc2;		lcycle07 += ( (-(int)((uint32)lcycle07 >> 31)) & nwt16);
					lcycle08 += wts_idx_inc2;		lcycle08 += ( (-(int)((uint32)lcycle08 >> 31)) & nwt16);
					lcycle09 += wts_idx_inc2;		lcycle09 += ( (-(int)((uint32)lcycle09 >> 31)) & nwt16);
					lcycle10 += wts_idx_inc2;		lcycle10 += ( (-(int)((uint32)lcycle10 >> 31)) & nwt16);
					lcycle11 += wts_idx_inc2;		lcycle11 += ( (-(int)((uint32)lcycle11 >> 31)) & nwt16);
					lcycle12 += wts_idx_inc2;		lcycle12 += ( (-(int)((uint32)lcycle12 >> 31)) & nwt16);
					lcycle13 += wts_idx_inc2;		lcycle13 += ( (-(int)((uint32)lcycle13 >> 31)) & nwt16);
					lcycle14 += wts_idx_inc2;		lcycle14 += ( (-(int)((uint32)lcycle14 >> 31)) & nwt16);
				  #endif
				#endif
				}
			}
			tdat[ithread].icycle00 = icycle00;
			tdat[ithread].icycle01 = icycle01;
			tdat[ithread].icycle02 = icycle02;
			tdat[ithread].icycle03 = icycle03;
			tdat[ithread].icycle04 = icycle04;
			tdat[ithread].icycle05 = icycle05;
			tdat[ithread].icycle06 = icycle06;
			tdat[ithread].icycle07 = icycle07;
			tdat[ithread].icycle08 = icycle08;
			tdat[ithread].icycle09 = icycle09;
			tdat[ithread].icycle10 = icycle10;
			tdat[ithread].icycle11 = icycle11;
			tdat[ithread].icycle12 = icycle12;
			tdat[ithread].icycle13 = icycle13;
			tdat[ithread].icycle14 = icycle14;

		#ifdef USE_SSE2
			tdat[ithread].wts_idx_inc2 = wts_idx_inc2;
			tdat[ithread].jcycle00 = jcycle00;
			tdat[ithread].jcycle01 = jcycle01;
			tdat[ithread].jcycle02 = jcycle02;
			tdat[ithread].jcycle03 = jcycle03;
			tdat[ithread].jcycle04 = jcycle04;
			tdat[ithread].jcycle05 = jcycle05;
			tdat[ithread].jcycle06 = jcycle06;
			tdat[ithread].jcycle07 = jcycle07;
			tdat[ithread].jcycle08 = jcycle08;
			tdat[ithread].jcycle09 = jcycle09;
			tdat[ithread].jcycle10 = jcycle10;
			tdat[ithread].jcycle11 = jcycle11;
			tdat[ithread].jcycle12 = jcycle12;
			tdat[ithread].jcycle13 = jcycle13;
			tdat[ithread].jcycle14 = jcycle14;

		  #ifdef USE_AVX
			tdat[ithread].kcycle00 = kcycle00;
			tdat[ithread].kcycle01 = kcycle01;
			tdat[ithread].kcycle02 = kcycle02;
			tdat[ithread].kcycle03 = kcycle03;
			tdat[ithread].kcycle04 = kcycle04;
			tdat[ithread].kcycle05 = kcycle05;
			tdat[ithread].kcycle06 = kcycle06;
			tdat[ithread].kcycle07 = kcycle07;
			tdat[ithread].kcycle08 = kcycle08;
			tdat[ithread].kcycle09 = kcycle09;
			tdat[ithread].kcycle10 = kcycle10;
			tdat[ithread].kcycle11 = kcycle11;
			tdat[ithread].kcycle12 = kcycle12;
			tdat[ithread].kcycle13 = kcycle13;
			tdat[ithread].kcycle14 = kcycle14;

			tdat[ithread].lcycle00 = lcycle00;
			tdat[ithread].lcycle01 = lcycle01;
			tdat[ithread].lcycle02 = lcycle02;
			tdat[ithread].lcycle03 = lcycle03;
			tdat[ithread].lcycle04 = lcycle04;
			tdat[ithread].lcycle05 = lcycle05;
			tdat[ithread].lcycle06 = lcycle06;
			tdat[ithread].lcycle07 = lcycle07;
			tdat[ithread].lcycle08 = lcycle08;
			tdat[ithread].lcycle09 = lcycle09;
			tdat[ithread].lcycle10 = lcycle10;
			tdat[ithread].lcycle11 = lcycle11;
			tdat[ithread].lcycle12 = lcycle12;
			tdat[ithread].lcycle13 = lcycle13;
			tdat[ithread].lcycle14 = lcycle14;
		  #endif
		#endif
		}
	}

  #endif

#if defined(USE_SSE2) && defined(USE_PTHREAD)

	tmp = max_err;	VEC_DBL_INIT(tmp, 0.0);
	tm2 = tmp + cslots_in_local_store;
	for(ithread = 1; ithread < CY_THREADS; ++ithread) {
		memcpy(tm2, tmp, sz_vd);
		tmp = tm2;		tm2 += cslots_in_local_store;
	}

#endif	// USE_PTHREAD

	/* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
	if(!full_pass)
	{
		khi = 1;
	}

#ifdef USE_PTHREAD
	/* Populate the thread-specific data structs - use the invariant terms as memchecks: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
	// int data:
		ASSERT(HERE, tdat[ithread].tid == ithread, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].ndivr == NDIVR, "thread-local memcheck fail!");

		tdat[ithread].khi    = khi;
		tdat[ithread].i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		tdat[ithread].jstart = _jstart[ithread];
		tdat[ithread].jhi    = _jhi[ithread];

//	printf("cy60_main: thread %d tdat-init, khi = %d, jlo = %d, jhi = %d\n", ithread,tdat[ithread].khi,tdat[ithread].jstart,tdat[ithread].jhi);

		tdat[ithread].col = _col[ithread];
		tdat[ithread].co2 = _co2[ithread];
		tdat[ithread].co3 = _co3[ithread];
		ASSERT(HERE, tdat[ithread].sw  == sw, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].nwt == nwt, "thread-local memcheck fail!");

	// double data:
		tdat[ithread].maxerr = _maxerr[ithread];
		tdat[ithread].scale = scale;

	// pointer data:
		ASSERT(HERE, tdat[ithread].arrdat == a, "thread-local memcheck fail!");			/* Main data array */
		ASSERT(HERE, tdat[ithread].wt0 == wt0, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].wt1 == wt1, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].si  == si, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].rn0 == rn0, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].rn1 == rn1, "thread-local memcheck fail!");
	#ifdef USE_SSE2
		ASSERT(HERE, tdat[ithread].wts_idx_inc2 == wts_idx_inc2, "thread-local memcheck fail!");
		ASSERT(HERE, tdat[ithread].s1p00r == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].s1p00r;
		ASSERT(HERE, ((tmp + 300)->d0 == c3m1 && (tmp + 300)->d1 == c3m1), "thread-local memcheck failed!");
		tmp = tdat[ithread].half_arr;
		ASSERT(HERE, ((tmp-1)->d0 == crnd && (tmp-1)->d1 == crnd), "thread-local memcheck failed!");
	#endif
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX
			// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
			dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#elif defined(USE_SSE2)
			dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		#endif
			tdat[ithread].bjmodn00 = _bjmodn00[ithread];
			tdat[ithread].bjmodn01 = _bjmodn01[ithread];
			tdat[ithread].bjmodn02 = _bjmodn02[ithread];
			tdat[ithread].bjmodn03 = _bjmodn03[ithread];
			tdat[ithread].bjmodn04 = _bjmodn04[ithread];
			tdat[ithread].bjmodn05 = _bjmodn05[ithread];
			tdat[ithread].bjmodn06 = _bjmodn06[ithread];
			tdat[ithread].bjmodn07 = _bjmodn07[ithread];
			tdat[ithread].bjmodn08 = _bjmodn08[ithread];
			tdat[ithread].bjmodn09 = _bjmodn09[ithread];
			tdat[ithread].bjmodn10 = _bjmodn10[ithread];
			tdat[ithread].bjmodn11 = _bjmodn11[ithread];
			tdat[ithread].bjmodn12 = _bjmodn12[ithread];
			tdat[ithread].bjmodn13 = _bjmodn13[ithread];
			tdat[ithread].bjmodn14 = _bjmodn14[ithread];
			tdat[ithread].bjmodn15 = _bjmodn15[ithread];
			tdat[ithread].bjmodn16 = _bjmodn16[ithread];
			tdat[ithread].bjmodn17 = _bjmodn17[ithread];
			tdat[ithread].bjmodn18 = _bjmodn18[ithread];
			tdat[ithread].bjmodn19 = _bjmodn19[ithread];
			tdat[ithread].bjmodn20 = _bjmodn20[ithread];
			tdat[ithread].bjmodn21 = _bjmodn21[ithread];
			tdat[ithread].bjmodn22 = _bjmodn22[ithread];
			tdat[ithread].bjmodn23 = _bjmodn23[ithread];
			tdat[ithread].bjmodn24 = _bjmodn24[ithread];
			tdat[ithread].bjmodn25 = _bjmodn25[ithread];
			tdat[ithread].bjmodn26 = _bjmodn26[ithread];
			tdat[ithread].bjmodn27 = _bjmodn27[ithread];
			tdat[ithread].bjmodn28 = _bjmodn28[ithread];
			tdat[ithread].bjmodn29 = _bjmodn29[ithread];
			tdat[ithread].bjmodn30 = _bjmodn30[ithread];
			tdat[ithread].bjmodn31 = _bjmodn31[ithread];
			tdat[ithread].bjmodn32 = _bjmodn32[ithread];
			tdat[ithread].bjmodn33 = _bjmodn33[ithread];
			tdat[ithread].bjmodn34 = _bjmodn34[ithread];
			tdat[ithread].bjmodn35 = _bjmodn35[ithread];
			tdat[ithread].bjmodn36 = _bjmodn36[ithread];
			tdat[ithread].bjmodn37 = _bjmodn37[ithread];
			tdat[ithread].bjmodn38 = _bjmodn38[ithread];
			tdat[ithread].bjmodn39 = _bjmodn39[ithread];
			tdat[ithread].bjmodn40 = _bjmodn40[ithread];
			tdat[ithread].bjmodn41 = _bjmodn41[ithread];
			tdat[ithread].bjmodn42 = _bjmodn42[ithread];
			tdat[ithread].bjmodn43 = _bjmodn43[ithread];
			tdat[ithread].bjmodn44 = _bjmodn44[ithread];
			tdat[ithread].bjmodn45 = _bjmodn45[ithread];
			tdat[ithread].bjmodn46 = _bjmodn46[ithread];
			tdat[ithread].bjmodn47 = _bjmodn47[ithread];
			tdat[ithread].bjmodn48 = _bjmodn48[ithread];
			tdat[ithread].bjmodn49 = _bjmodn49[ithread];
			tdat[ithread].bjmodn50 = _bjmodn50[ithread];
			tdat[ithread].bjmodn51 = _bjmodn51[ithread];
			tdat[ithread].bjmodn52 = _bjmodn52[ithread];
			tdat[ithread].bjmodn53 = _bjmodn53[ithread];
			tdat[ithread].bjmodn54 = _bjmodn54[ithread];
			tdat[ithread].bjmodn55 = _bjmodn55[ithread];
			tdat[ithread].bjmodn56 = _bjmodn56[ithread];
			tdat[ithread].bjmodn57 = _bjmodn57[ithread];
			tdat[ithread].bjmodn58 = _bjmodn58[ithread];
			tdat[ithread].bjmodn59 = _bjmodn59[ithread];
			/* init carries	*/
			tdat[ithread].cy_r00 = _cy_r00[ithread];
			tdat[ithread].cy_r01 = _cy_r01[ithread];
			tdat[ithread].cy_r02 = _cy_r02[ithread];
			tdat[ithread].cy_r03 = _cy_r03[ithread];
			tdat[ithread].cy_r04 = _cy_r04[ithread];
			tdat[ithread].cy_r05 = _cy_r05[ithread];
			tdat[ithread].cy_r06 = _cy_r06[ithread];
			tdat[ithread].cy_r07 = _cy_r07[ithread];
			tdat[ithread].cy_r08 = _cy_r08[ithread];
			tdat[ithread].cy_r09 = _cy_r09[ithread];
			tdat[ithread].cy_r10 = _cy_r10[ithread];
			tdat[ithread].cy_r11 = _cy_r11[ithread];
			tdat[ithread].cy_r12 = _cy_r12[ithread];
			tdat[ithread].cy_r13 = _cy_r13[ithread];
			tdat[ithread].cy_r14 = _cy_r14[ithread];
			tdat[ithread].cy_r15 = _cy_r15[ithread];
			tdat[ithread].cy_r16 = _cy_r16[ithread];
			tdat[ithread].cy_r17 = _cy_r17[ithread];
			tdat[ithread].cy_r18 = _cy_r18[ithread];
			tdat[ithread].cy_r19 = _cy_r19[ithread];
			tdat[ithread].cy_r20 = _cy_r20[ithread];
			tdat[ithread].cy_r21 = _cy_r21[ithread];
			tdat[ithread].cy_r22 = _cy_r22[ithread];
			tdat[ithread].cy_r23 = _cy_r23[ithread];
			tdat[ithread].cy_r24 = _cy_r24[ithread];
			tdat[ithread].cy_r25 = _cy_r25[ithread];
			tdat[ithread].cy_r26 = _cy_r26[ithread];
			tdat[ithread].cy_r27 = _cy_r27[ithread];
			tdat[ithread].cy_r28 = _cy_r28[ithread];
			tdat[ithread].cy_r29 = _cy_r29[ithread];
			tdat[ithread].cy_r30 = _cy_r30[ithread];
			tdat[ithread].cy_r31 = _cy_r31[ithread];
			tdat[ithread].cy_r32 = _cy_r32[ithread];
			tdat[ithread].cy_r33 = _cy_r33[ithread];
			tdat[ithread].cy_r34 = _cy_r34[ithread];
			tdat[ithread].cy_r35 = _cy_r35[ithread];
			tdat[ithread].cy_r36 = _cy_r36[ithread];
			tdat[ithread].cy_r37 = _cy_r37[ithread];
			tdat[ithread].cy_r38 = _cy_r38[ithread];
			tdat[ithread].cy_r39 = _cy_r39[ithread];
			tdat[ithread].cy_r40 = _cy_r40[ithread];
			tdat[ithread].cy_r41 = _cy_r41[ithread];
			tdat[ithread].cy_r42 = _cy_r42[ithread];
			tdat[ithread].cy_r43 = _cy_r43[ithread];
			tdat[ithread].cy_r44 = _cy_r44[ithread];
			tdat[ithread].cy_r45 = _cy_r45[ithread];
			tdat[ithread].cy_r46 = _cy_r46[ithread];
			tdat[ithread].cy_r47 = _cy_r47[ithread];
			tdat[ithread].cy_r48 = _cy_r48[ithread];
			tdat[ithread].cy_r49 = _cy_r49[ithread];
			tdat[ithread].cy_r50 = _cy_r50[ithread];
			tdat[ithread].cy_r51 = _cy_r51[ithread];
			tdat[ithread].cy_r52 = _cy_r52[ithread];
			tdat[ithread].cy_r53 = _cy_r53[ithread];
			tdat[ithread].cy_r54 = _cy_r54[ithread];
			tdat[ithread].cy_r55 = _cy_r55[ithread];
			tdat[ithread].cy_r56 = _cy_r56[ithread];
			tdat[ithread].cy_r57 = _cy_r57[ithread];
			tdat[ithread].cy_r58 = _cy_r58[ithread];
			tdat[ithread].cy_r59 = _cy_r59[ithread];
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
		#ifdef USE_SSE2
			dtmp = (tmp)->d0 * (tmp+odd_radix)->d0;	ASSERT(HERE, fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
			dtmp = (tmp)->d1 * (tmp+odd_radix)->d1;	ASSERT(HERE, fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
		#endif
			/* init carries	*/
			tdat[ithread].cy_r00 = _cy_r00[ithread];	tdat[ithread].cy_i00 = _cy_i00[ithread];
			tdat[ithread].cy_r01 = _cy_r01[ithread];	tdat[ithread].cy_i01 = _cy_i01[ithread];
			tdat[ithread].cy_r02 = _cy_r02[ithread];	tdat[ithread].cy_i02 = _cy_i02[ithread];
			tdat[ithread].cy_r03 = _cy_r03[ithread];	tdat[ithread].cy_i03 = _cy_i03[ithread];
			tdat[ithread].cy_r04 = _cy_r04[ithread];	tdat[ithread].cy_i04 = _cy_i04[ithread];
			tdat[ithread].cy_r05 = _cy_r05[ithread];	tdat[ithread].cy_i05 = _cy_i05[ithread];
			tdat[ithread].cy_r06 = _cy_r06[ithread];	tdat[ithread].cy_i06 = _cy_i06[ithread];
			tdat[ithread].cy_r07 = _cy_r07[ithread];	tdat[ithread].cy_i07 = _cy_i07[ithread];
			tdat[ithread].cy_r08 = _cy_r08[ithread];	tdat[ithread].cy_i08 = _cy_i08[ithread];
			tdat[ithread].cy_r09 = _cy_r09[ithread];	tdat[ithread].cy_i09 = _cy_i09[ithread];
			tdat[ithread].cy_r10 = _cy_r10[ithread];	tdat[ithread].cy_i10 = _cy_i10[ithread];
			tdat[ithread].cy_r11 = _cy_r11[ithread];	tdat[ithread].cy_i11 = _cy_i11[ithread];
			tdat[ithread].cy_r12 = _cy_r12[ithread];	tdat[ithread].cy_i12 = _cy_i12[ithread];
			tdat[ithread].cy_r13 = _cy_r13[ithread];	tdat[ithread].cy_i13 = _cy_i13[ithread];
			tdat[ithread].cy_r14 = _cy_r14[ithread];	tdat[ithread].cy_i14 = _cy_i14[ithread];
			tdat[ithread].cy_r15 = _cy_r15[ithread];	tdat[ithread].cy_i15 = _cy_i15[ithread];
			tdat[ithread].cy_r16 = _cy_r16[ithread];	tdat[ithread].cy_i16 = _cy_i16[ithread];
			tdat[ithread].cy_r17 = _cy_r17[ithread];	tdat[ithread].cy_i17 = _cy_i17[ithread];
			tdat[ithread].cy_r18 = _cy_r18[ithread];	tdat[ithread].cy_i18 = _cy_i18[ithread];
			tdat[ithread].cy_r19 = _cy_r19[ithread];	tdat[ithread].cy_i19 = _cy_i19[ithread];
			tdat[ithread].cy_r20 = _cy_r20[ithread];	tdat[ithread].cy_i20 = _cy_i20[ithread];
			tdat[ithread].cy_r21 = _cy_r21[ithread];	tdat[ithread].cy_i21 = _cy_i21[ithread];
			tdat[ithread].cy_r22 = _cy_r22[ithread];	tdat[ithread].cy_i22 = _cy_i22[ithread];
			tdat[ithread].cy_r23 = _cy_r23[ithread];	tdat[ithread].cy_i23 = _cy_i23[ithread];
			tdat[ithread].cy_r24 = _cy_r24[ithread];	tdat[ithread].cy_i24 = _cy_i24[ithread];
			tdat[ithread].cy_r25 = _cy_r25[ithread];	tdat[ithread].cy_i25 = _cy_i25[ithread];
			tdat[ithread].cy_r26 = _cy_r26[ithread];	tdat[ithread].cy_i26 = _cy_i26[ithread];
			tdat[ithread].cy_r27 = _cy_r27[ithread];	tdat[ithread].cy_i27 = _cy_i27[ithread];
			tdat[ithread].cy_r28 = _cy_r28[ithread];	tdat[ithread].cy_i28 = _cy_i28[ithread];
			tdat[ithread].cy_r29 = _cy_r29[ithread];	tdat[ithread].cy_i29 = _cy_i29[ithread];
			tdat[ithread].cy_r30 = _cy_r30[ithread];	tdat[ithread].cy_i30 = _cy_i30[ithread];
			tdat[ithread].cy_r31 = _cy_r31[ithread];	tdat[ithread].cy_i31 = _cy_i31[ithread];
			tdat[ithread].cy_r32 = _cy_r32[ithread];	tdat[ithread].cy_i32 = _cy_i32[ithread];
			tdat[ithread].cy_r33 = _cy_r33[ithread];	tdat[ithread].cy_i33 = _cy_i33[ithread];
			tdat[ithread].cy_r34 = _cy_r34[ithread];	tdat[ithread].cy_i34 = _cy_i34[ithread];
			tdat[ithread].cy_r35 = _cy_r35[ithread];	tdat[ithread].cy_i35 = _cy_i35[ithread];
			tdat[ithread].cy_r36 = _cy_r36[ithread];	tdat[ithread].cy_i36 = _cy_i36[ithread];
			tdat[ithread].cy_r37 = _cy_r37[ithread];	tdat[ithread].cy_i37 = _cy_i37[ithread];
			tdat[ithread].cy_r38 = _cy_r38[ithread];	tdat[ithread].cy_i38 = _cy_i38[ithread];
			tdat[ithread].cy_r39 = _cy_r39[ithread];	tdat[ithread].cy_i39 = _cy_i39[ithread];
			tdat[ithread].cy_r40 = _cy_r40[ithread];	tdat[ithread].cy_i40 = _cy_i40[ithread];
			tdat[ithread].cy_r41 = _cy_r41[ithread];	tdat[ithread].cy_i41 = _cy_i41[ithread];
			tdat[ithread].cy_r42 = _cy_r42[ithread];	tdat[ithread].cy_i42 = _cy_i42[ithread];
			tdat[ithread].cy_r43 = _cy_r43[ithread];	tdat[ithread].cy_i43 = _cy_i43[ithread];
			tdat[ithread].cy_r44 = _cy_r44[ithread];	tdat[ithread].cy_i44 = _cy_i44[ithread];
			tdat[ithread].cy_r45 = _cy_r45[ithread];	tdat[ithread].cy_i45 = _cy_i45[ithread];
			tdat[ithread].cy_r46 = _cy_r46[ithread];	tdat[ithread].cy_i46 = _cy_i46[ithread];
			tdat[ithread].cy_r47 = _cy_r47[ithread];	tdat[ithread].cy_i47 = _cy_i47[ithread];
			tdat[ithread].cy_r48 = _cy_r48[ithread];	tdat[ithread].cy_i48 = _cy_i48[ithread];
			tdat[ithread].cy_r49 = _cy_r49[ithread];	tdat[ithread].cy_i49 = _cy_i49[ithread];
			tdat[ithread].cy_r50 = _cy_r50[ithread];	tdat[ithread].cy_i50 = _cy_i50[ithread];
			tdat[ithread].cy_r51 = _cy_r51[ithread];	tdat[ithread].cy_i51 = _cy_i51[ithread];
			tdat[ithread].cy_r52 = _cy_r52[ithread];	tdat[ithread].cy_i52 = _cy_i52[ithread];
			tdat[ithread].cy_r53 = _cy_r53[ithread];	tdat[ithread].cy_i53 = _cy_i53[ithread];
			tdat[ithread].cy_r54 = _cy_r54[ithread];	tdat[ithread].cy_i54 = _cy_i54[ithread];
			tdat[ithread].cy_r55 = _cy_r55[ithread];	tdat[ithread].cy_i55 = _cy_i55[ithread];
			tdat[ithread].cy_r56 = _cy_r56[ithread];	tdat[ithread].cy_i56 = _cy_i56[ithread];
			tdat[ithread].cy_r57 = _cy_r57[ithread];	tdat[ithread].cy_i57 = _cy_i57[ithread];
			tdat[ithread].cy_r58 = _cy_r58[ithread];	tdat[ithread].cy_i58 = _cy_i58[ithread];
			tdat[ithread].cy_r59 = _cy_r59[ithread];	tdat[ithread].cy_i59 = _cy_i59[ithread];
		}
	}
#endif

#ifdef USE_PTHREAD

	// If also using main thread to do work units, that task-dispatch occurs after all the threadpool-task launches:
	for(ithread = 0; ithread < pool_work_units; ithread++)
	{
		task_control.data = (void*)(&tdat[ithread]);
		threadpool_add_task(tpool, &task_control, task_is_blocking);

#else

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		/***** DEC/HP CC doesn't properly copy init value of maxerr = 0 into threads,
		so need to set once again explicitly for each: *****/
		maxerr = 0.0;
	#ifdef USE_SSE2
	//	VEC_DBL_INIT(max_err, 0.0);	*** must do this in conjunction with thread-local-data-copy
	#endif

	#ifdef USE_SSE2
		*bjmodn00 = _bjmodn00[ithread];
		*bjmodn01 = _bjmodn01[ithread];
		*bjmodn02 = _bjmodn02[ithread];
		*bjmodn03 = _bjmodn03[ithread];
		*bjmodn04 = _bjmodn04[ithread];
		*bjmodn05 = _bjmodn05[ithread];
		*bjmodn06 = _bjmodn06[ithread];
		*bjmodn07 = _bjmodn07[ithread];
		*bjmodn08 = _bjmodn08[ithread];
		*bjmodn09 = _bjmodn09[ithread];
		*bjmodn10 = _bjmodn10[ithread];
		*bjmodn11 = _bjmodn11[ithread];
		*bjmodn12 = _bjmodn12[ithread];
		*bjmodn13 = _bjmodn13[ithread];
		*bjmodn14 = _bjmodn14[ithread];
		*bjmodn15 = _bjmodn15[ithread];
		*bjmodn16 = _bjmodn16[ithread];
		*bjmodn17 = _bjmodn17[ithread];
		*bjmodn18 = _bjmodn18[ithread];
		*bjmodn19 = _bjmodn19[ithread];
		*bjmodn20 = _bjmodn20[ithread];
		*bjmodn21 = _bjmodn21[ithread];
		*bjmodn22 = _bjmodn22[ithread];
		*bjmodn23 = _bjmodn23[ithread];
		*bjmodn24 = _bjmodn24[ithread];
		*bjmodn25 = _bjmodn25[ithread];
		*bjmodn26 = _bjmodn26[ithread];
		*bjmodn27 = _bjmodn27[ithread];
		*bjmodn28 = _bjmodn28[ithread];
		*bjmodn29 = _bjmodn29[ithread];
		*bjmodn30 = _bjmodn30[ithread];
		*bjmodn31 = _bjmodn31[ithread];
		*bjmodn32 = _bjmodn32[ithread];
		*bjmodn33 = _bjmodn33[ithread];
		*bjmodn34 = _bjmodn34[ithread];
		*bjmodn35 = _bjmodn35[ithread];
		*bjmodn36 = _bjmodn36[ithread];
		*bjmodn37 = _bjmodn37[ithread];
		*bjmodn38 = _bjmodn38[ithread];
		*bjmodn39 = _bjmodn39[ithread];
		*bjmodn40 = _bjmodn40[ithread];
		*bjmodn41 = _bjmodn41[ithread];
		*bjmodn42 = _bjmodn42[ithread];
		*bjmodn43 = _bjmodn43[ithread];
		*bjmodn44 = _bjmodn44[ithread];
		*bjmodn45 = _bjmodn45[ithread];
		*bjmodn46 = _bjmodn46[ithread];
		*bjmodn47 = _bjmodn47[ithread];
		*bjmodn48 = _bjmodn48[ithread];
		*bjmodn49 = _bjmodn49[ithread];
		*bjmodn50 = _bjmodn50[ithread];
		*bjmodn51 = _bjmodn51[ithread];
		*bjmodn52 = _bjmodn52[ithread];
		*bjmodn53 = _bjmodn53[ithread];
		*bjmodn54 = _bjmodn54[ithread];
		*bjmodn55 = _bjmodn55[ithread];
		*bjmodn56 = _bjmodn56[ithread];
		*bjmodn57 = _bjmodn57[ithread];
		*bjmodn58 = _bjmodn58[ithread];
		*bjmodn59 = _bjmodn59[ithread];
	#else
		bjmodn00 = _bjmodn00[ithread];
		bjmodn01 = _bjmodn01[ithread];
		bjmodn02 = _bjmodn02[ithread];
		bjmodn03 = _bjmodn03[ithread];
		bjmodn04 = _bjmodn04[ithread];
		bjmodn05 = _bjmodn05[ithread];
		bjmodn06 = _bjmodn06[ithread];
		bjmodn07 = _bjmodn07[ithread];
		bjmodn08 = _bjmodn08[ithread];
		bjmodn09 = _bjmodn09[ithread];
		bjmodn10 = _bjmodn10[ithread];
		bjmodn11 = _bjmodn11[ithread];
		bjmodn12 = _bjmodn12[ithread];
		bjmodn13 = _bjmodn13[ithread];
		bjmodn14 = _bjmodn14[ithread];
		bjmodn15 = _bjmodn15[ithread];
		bjmodn16 = _bjmodn16[ithread];
		bjmodn17 = _bjmodn17[ithread];
		bjmodn18 = _bjmodn18[ithread];
		bjmodn19 = _bjmodn19[ithread];
		bjmodn20 = _bjmodn20[ithread];
		bjmodn21 = _bjmodn21[ithread];
		bjmodn22 = _bjmodn22[ithread];
		bjmodn23 = _bjmodn23[ithread];
		bjmodn24 = _bjmodn24[ithread];
		bjmodn25 = _bjmodn25[ithread];
		bjmodn26 = _bjmodn26[ithread];
		bjmodn27 = _bjmodn27[ithread];
		bjmodn28 = _bjmodn28[ithread];
		bjmodn29 = _bjmodn29[ithread];
		bjmodn30 = _bjmodn30[ithread];
		bjmodn31 = _bjmodn31[ithread];
		bjmodn32 = _bjmodn32[ithread];
		bjmodn33 = _bjmodn33[ithread];
		bjmodn34 = _bjmodn34[ithread];
		bjmodn35 = _bjmodn35[ithread];
		bjmodn36 = _bjmodn36[ithread];
		bjmodn37 = _bjmodn37[ithread];
		bjmodn38 = _bjmodn38[ithread];
		bjmodn39 = _bjmodn39[ithread];
		bjmodn40 = _bjmodn40[ithread];
		bjmodn41 = _bjmodn41[ithread];
		bjmodn42 = _bjmodn42[ithread];
		bjmodn43 = _bjmodn43[ithread];
		bjmodn44 = _bjmodn44[ithread];
		bjmodn45 = _bjmodn45[ithread];
		bjmodn46 = _bjmodn46[ithread];
		bjmodn47 = _bjmodn47[ithread];
		bjmodn48 = _bjmodn48[ithread];
		bjmodn49 = _bjmodn49[ithread];
		bjmodn50 = _bjmodn50[ithread];
		bjmodn51 = _bjmodn51[ithread];
		bjmodn52 = _bjmodn52[ithread];
		bjmodn53 = _bjmodn53[ithread];
		bjmodn54 = _bjmodn54[ithread];
		bjmodn55 = _bjmodn55[ithread];
		bjmodn56 = _bjmodn56[ithread];
		bjmodn57 = _bjmodn57[ithread];
		bjmodn58 = _bjmodn58[ithread];
		bjmodn59 = _bjmodn59[ithread];
	#endif

		i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		jstart = _jstart[ithread];
		jhi    = _jhi[ithread];

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			col = _col[ithread];
			co2 = _co2[ithread];
			co3 = _co3[ithread];

			/* init carries	*/
		#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
			cy_r00->d0 = _cy_r00[ithread];	cy_r00->d1 = _cy_r01[ithread];	cy_r00->d2 = _cy_r02[ithread];	cy_r00->d3 = _cy_r03[ithread];
			cy_r04->d0 = _cy_r04[ithread];	cy_r04->d1 = _cy_r05[ithread];	cy_r04->d2 = _cy_r06[ithread];	cy_r04->d3 = _cy_r07[ithread];
			cy_r08->d0 = _cy_r08[ithread];	cy_r08->d1 = _cy_r09[ithread];	cy_r08->d2 = _cy_r10[ithread];	cy_r08->d3 = _cy_r11[ithread];
			cy_r12->d0 = _cy_r12[ithread];	cy_r12->d1 = _cy_r13[ithread];	cy_r12->d2 = _cy_r14[ithread];	cy_r12->d3 = _cy_r15[ithread];
			cy_r16->d0 = _cy_r16[ithread];	cy_r16->d1 = _cy_r17[ithread];	cy_r16->d2 = _cy_r18[ithread];	cy_r16->d3 = _cy_r19[ithread];
			cy_r20->d0 = _cy_r20[ithread];	cy_r20->d1 = _cy_r21[ithread];	cy_r20->d2 = _cy_r22[ithread];	cy_r20->d3 = _cy_r23[ithread];
			cy_r24->d0 = _cy_r24[ithread];	cy_r24->d1 = _cy_r25[ithread];	cy_r24->d2 = _cy_r26[ithread];	cy_r24->d3 = _cy_r27[ithread];
			cy_r28->d0 = _cy_r28[ithread];	cy_r28->d1 = _cy_r29[ithread];	cy_r28->d2 = _cy_r30[ithread];	cy_r28->d3 = _cy_r31[ithread];
			cy_r32->d0 = _cy_r32[ithread];	cy_r32->d1 = _cy_r33[ithread];	cy_r32->d2 = _cy_r34[ithread];	cy_r32->d3 = _cy_r35[ithread];
			cy_r36->d0 = _cy_r36[ithread];	cy_r36->d1 = _cy_r37[ithread];	cy_r36->d2 = _cy_r38[ithread];	cy_r36->d3 = _cy_r39[ithread];
			cy_r40->d0 = _cy_r40[ithread];	cy_r40->d1 = _cy_r41[ithread];	cy_r40->d2 = _cy_r42[ithread];	cy_r40->d3 = _cy_r43[ithread];
			cy_r44->d0 = _cy_r44[ithread];	cy_r44->d1 = _cy_r45[ithread];	cy_r44->d2 = _cy_r46[ithread];	cy_r44->d3 = _cy_r47[ithread];
			cy_r48->d0 = _cy_r48[ithread];	cy_r48->d1 = _cy_r49[ithread];	cy_r48->d2 = _cy_r50[ithread];	cy_r48->d3 = _cy_r51[ithread];
			cy_r52->d0 = _cy_r52[ithread];	cy_r52->d1 = _cy_r53[ithread];	cy_r52->d2 = _cy_r54[ithread];	cy_r52->d3 = _cy_r55[ithread];
			cy_r56->d0 = _cy_r56[ithread];	cy_r56->d1 = _cy_r57[ithread];	cy_r56->d2 = _cy_r58[ithread];	cy_r56->d3 = _cy_r59[ithread];
		#elif defined(USE_SSE2)
			cy_r00->d0 = _cy_r00[ithread];	cy_r00->d1 = _cy_r01[ithread];
			cy_r02->d0 = _cy_r02[ithread];	cy_r02->d1 = _cy_r03[ithread];
			cy_r04->d0 = _cy_r04[ithread];	cy_r04->d1 = _cy_r05[ithread];
			cy_r06->d0 = _cy_r06[ithread];	cy_r06->d1 = _cy_r07[ithread];
			cy_r08->d0 = _cy_r08[ithread];	cy_r08->d1 = _cy_r09[ithread];
			cy_r10->d0 = _cy_r10[ithread];	cy_r10->d1 = _cy_r11[ithread];
			cy_r12->d0 = _cy_r12[ithread];	cy_r12->d1 = _cy_r13[ithread];
			cy_r14->d0 = _cy_r14[ithread];	cy_r14->d1 = _cy_r15[ithread];
			cy_r16->d0 = _cy_r16[ithread];	cy_r16->d1 = _cy_r17[ithread];
			cy_r18->d0 = _cy_r18[ithread];	cy_r18->d1 = _cy_r19[ithread];
			cy_r20->d0 = _cy_r20[ithread];	cy_r20->d1 = _cy_r21[ithread];
			cy_r22->d0 = _cy_r22[ithread];	cy_r22->d1 = _cy_r23[ithread];
			cy_r24->d0 = _cy_r24[ithread];	cy_r24->d1 = _cy_r25[ithread];
			cy_r26->d0 = _cy_r26[ithread];	cy_r26->d1 = _cy_r27[ithread];
			cy_r28->d0 = _cy_r28[ithread];	cy_r28->d1 = _cy_r29[ithread];
			cy_r30->d0 = _cy_r30[ithread];	cy_r30->d1 = _cy_r31[ithread];
			cy_r32->d0 = _cy_r32[ithread];	cy_r32->d1 = _cy_r33[ithread];
			cy_r34->d0 = _cy_r34[ithread];	cy_r34->d1 = _cy_r35[ithread];
			cy_r36->d0 = _cy_r36[ithread];	cy_r36->d1 = _cy_r37[ithread];
			cy_r38->d0 = _cy_r38[ithread];	cy_r38->d1 = _cy_r39[ithread];
			cy_r40->d0 = _cy_r40[ithread];	cy_r40->d1 = _cy_r41[ithread];
			cy_r42->d0 = _cy_r42[ithread];	cy_r42->d1 = _cy_r43[ithread];
			cy_r44->d0 = _cy_r44[ithread];	cy_r44->d1 = _cy_r45[ithread];
			cy_r46->d0 = _cy_r46[ithread];	cy_r46->d1 = _cy_r47[ithread];
			cy_r48->d0 = _cy_r48[ithread];	cy_r48->d1 = _cy_r49[ithread];
			cy_r50->d0 = _cy_r50[ithread];	cy_r50->d1 = _cy_r51[ithread];
			cy_r52->d0 = _cy_r52[ithread];	cy_r52->d1 = _cy_r53[ithread];
			cy_r54->d0 = _cy_r54[ithread];	cy_r54->d1 = _cy_r55[ithread];
			cy_r56->d0 = _cy_r56[ithread];	cy_r56->d1 = _cy_r57[ithread];
			cy_r58->d0 = _cy_r58[ithread];	cy_r58->d1 = _cy_r59[ithread];
		#else
			/* init carries	*/
			cy_r00 = _cy_r00[ithread];
			cy_r01 = _cy_r01[ithread];
			cy_r02 = _cy_r02[ithread];
			cy_r03 = _cy_r03[ithread];
			cy_r04 = _cy_r04[ithread];
			cy_r05 = _cy_r05[ithread];
			cy_r06 = _cy_r06[ithread];
			cy_r07 = _cy_r07[ithread];
			cy_r08 = _cy_r08[ithread];
			cy_r09 = _cy_r09[ithread];
			cy_r10 = _cy_r10[ithread];
			cy_r11 = _cy_r11[ithread];
			cy_r12 = _cy_r12[ithread];
			cy_r13 = _cy_r13[ithread];
			cy_r14 = _cy_r14[ithread];
			cy_r15 = _cy_r15[ithread];
			cy_r16 = _cy_r16[ithread];
			cy_r17 = _cy_r17[ithread];
			cy_r18 = _cy_r18[ithread];
			cy_r19 = _cy_r19[ithread];
			cy_r20 = _cy_r20[ithread];
			cy_r21 = _cy_r21[ithread];
			cy_r22 = _cy_r22[ithread];
			cy_r23 = _cy_r23[ithread];
			cy_r24 = _cy_r24[ithread];
			cy_r25 = _cy_r25[ithread];
			cy_r26 = _cy_r26[ithread];
			cy_r27 = _cy_r27[ithread];
			cy_r28 = _cy_r28[ithread];
			cy_r29 = _cy_r29[ithread];
			cy_r30 = _cy_r30[ithread];
			cy_r31 = _cy_r31[ithread];
			cy_r32 = _cy_r32[ithread];
			cy_r33 = _cy_r33[ithread];
			cy_r34 = _cy_r34[ithread];
			cy_r35 = _cy_r35[ithread];
			cy_r36 = _cy_r36[ithread];
			cy_r37 = _cy_r37[ithread];
			cy_r38 = _cy_r38[ithread];
			cy_r39 = _cy_r39[ithread];
			cy_r40 = _cy_r40[ithread];
			cy_r41 = _cy_r41[ithread];
			cy_r42 = _cy_r42[ithread];
			cy_r43 = _cy_r43[ithread];
			cy_r44 = _cy_r44[ithread];
			cy_r45 = _cy_r45[ithread];
			cy_r46 = _cy_r46[ithread];
			cy_r47 = _cy_r47[ithread];
			cy_r48 = _cy_r48[ithread];
			cy_r49 = _cy_r49[ithread];
			cy_r50 = _cy_r50[ithread];
			cy_r51 = _cy_r51[ithread];
			cy_r52 = _cy_r52[ithread];
			cy_r53 = _cy_r53[ithread];
			cy_r54 = _cy_r54[ithread];
			cy_r55 = _cy_r55[ithread];
			cy_r56 = _cy_r56[ithread];
			cy_r57 = _cy_r57[ithread];
			cy_r58 = _cy_r58[ithread];
			cy_r59 = _cy_r59[ithread];
		#endif
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
		#ifdef USE_AVX
			cy_r00->d0 = _cy_r00[ithread];	cy_r00->d1 = _cy_r01[ithread];	cy_r00->d2 = _cy_r02[ithread];	cy_r00->d3 = _cy_r03[ithread];
			cy_r04->d0 = _cy_r04[ithread];	cy_r04->d1 = _cy_r05[ithread];	cy_r04->d2 = _cy_r06[ithread];	cy_r04->d3 = _cy_r07[ithread];
			cy_r08->d0 = _cy_r08[ithread];	cy_r08->d1 = _cy_r09[ithread];	cy_r08->d2 = _cy_r10[ithread];	cy_r08->d3 = _cy_r11[ithread];
			cy_r12->d0 = _cy_r12[ithread];	cy_r12->d1 = _cy_r13[ithread];	cy_r12->d2 = _cy_r14[ithread];	cy_r12->d3 = _cy_r15[ithread];
			cy_r16->d0 = _cy_r16[ithread];	cy_r16->d1 = _cy_r17[ithread];	cy_r16->d2 = _cy_r18[ithread];	cy_r16->d3 = _cy_r19[ithread];
			cy_r20->d0 = _cy_r20[ithread];	cy_r20->d1 = _cy_r21[ithread];	cy_r20->d2 = _cy_r22[ithread];	cy_r20->d3 = _cy_r23[ithread];
			cy_r24->d0 = _cy_r24[ithread];	cy_r24->d1 = _cy_r25[ithread];	cy_r24->d2 = _cy_r26[ithread];	cy_r24->d3 = _cy_r27[ithread];
			cy_r28->d0 = _cy_r28[ithread];	cy_r28->d1 = _cy_r29[ithread];	cy_r28->d2 = _cy_r30[ithread];	cy_r28->d3 = _cy_r31[ithread];
			cy_r32->d0 = _cy_r32[ithread];	cy_r32->d1 = _cy_r33[ithread];	cy_r32->d2 = _cy_r34[ithread];	cy_r32->d3 = _cy_r35[ithread];
			cy_r36->d0 = _cy_r36[ithread];	cy_r36->d1 = _cy_r37[ithread];	cy_r36->d2 = _cy_r38[ithread];	cy_r36->d3 = _cy_r39[ithread];
			cy_r40->d0 = _cy_r40[ithread];	cy_r40->d1 = _cy_r41[ithread];	cy_r40->d2 = _cy_r42[ithread];	cy_r40->d3 = _cy_r43[ithread];
			cy_r44->d0 = _cy_r44[ithread];	cy_r44->d1 = _cy_r45[ithread];	cy_r44->d2 = _cy_r46[ithread];	cy_r44->d3 = _cy_r47[ithread];
			cy_r48->d0 = _cy_r48[ithread];	cy_r48->d1 = _cy_r49[ithread];	cy_r48->d2 = _cy_r50[ithread];	cy_r48->d3 = _cy_r51[ithread];
			cy_r52->d0 = _cy_r52[ithread];	cy_r52->d1 = _cy_r53[ithread];	cy_r52->d2 = _cy_r54[ithread];	cy_r52->d3 = _cy_r55[ithread];
			cy_r56->d0 = _cy_r56[ithread];	cy_r56->d1 = _cy_r57[ithread];	cy_r56->d2 = _cy_r58[ithread];	cy_r56->d3 = _cy_r59[ithread];

			cy_i00->d0 = _cy_i00[ithread];	cy_i00->d1 = _cy_i01[ithread];	cy_i00->d2 = _cy_i02[ithread];	cy_i00->d3 = _cy_i03[ithread];
			cy_i04->d0 = _cy_i04[ithread];	cy_i04->d1 = _cy_i05[ithread];	cy_i04->d2 = _cy_i06[ithread];	cy_i04->d3 = _cy_i07[ithread];
			cy_i08->d0 = _cy_i08[ithread];	cy_i08->d1 = _cy_i09[ithread];	cy_i08->d2 = _cy_i10[ithread];	cy_i08->d3 = _cy_i11[ithread];
			cy_i12->d0 = _cy_i12[ithread];	cy_i12->d1 = _cy_i13[ithread];	cy_i12->d2 = _cy_i14[ithread];	cy_i12->d3 = _cy_i15[ithread];
			cy_i16->d0 = _cy_i16[ithread];	cy_i16->d1 = _cy_i17[ithread];	cy_i16->d2 = _cy_i18[ithread];	cy_i16->d3 = _cy_i19[ithread];
			cy_i20->d0 = _cy_i20[ithread];	cy_i20->d1 = _cy_i21[ithread];	cy_i20->d2 = _cy_i22[ithread];	cy_i20->d3 = _cy_i23[ithread];
			cy_i24->d0 = _cy_i24[ithread];	cy_i24->d1 = _cy_i25[ithread];	cy_i24->d2 = _cy_i26[ithread];	cy_i24->d3 = _cy_i27[ithread];
			cy_i28->d0 = _cy_i28[ithread];	cy_i28->d1 = _cy_i29[ithread];	cy_i28->d2 = _cy_i30[ithread];	cy_i28->d3 = _cy_i31[ithread];
			cy_i32->d0 = _cy_i32[ithread];	cy_i32->d1 = _cy_i33[ithread];	cy_i32->d2 = _cy_i34[ithread];	cy_i32->d3 = _cy_i35[ithread];
			cy_i36->d0 = _cy_i36[ithread];	cy_i36->d1 = _cy_i37[ithread];	cy_i36->d2 = _cy_i38[ithread];	cy_i36->d3 = _cy_i39[ithread];
			cy_i40->d0 = _cy_i40[ithread];	cy_i40->d1 = _cy_i41[ithread];	cy_i40->d2 = _cy_i42[ithread];	cy_i40->d3 = _cy_i43[ithread];
			cy_i44->d0 = _cy_i44[ithread];	cy_i44->d1 = _cy_i45[ithread];	cy_i44->d2 = _cy_i46[ithread];	cy_i44->d3 = _cy_i47[ithread];
			cy_i48->d0 = _cy_i48[ithread];	cy_i48->d1 = _cy_i49[ithread];	cy_i48->d2 = _cy_i50[ithread];	cy_i48->d3 = _cy_i51[ithread];
			cy_i52->d0 = _cy_i52[ithread];	cy_i52->d1 = _cy_i53[ithread];	cy_i52->d2 = _cy_i54[ithread];	cy_i52->d3 = _cy_i55[ithread];
			cy_i56->d0 = _cy_i56[ithread];	cy_i56->d1 = _cy_i57[ithread];	cy_i56->d2 = _cy_i58[ithread];	cy_i56->d3 = _cy_i59[ithread];
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			cy_r00->d0 = _cy_r00[ithread];	cy_r00->d1 = _cy_i00[ithread];
			cy_r02->d0 = _cy_r01[ithread];	cy_r02->d1 = _cy_i01[ithread];
			cy_r04->d0 = _cy_r02[ithread];	cy_r04->d1 = _cy_i02[ithread];
			cy_r06->d0 = _cy_r03[ithread];	cy_r06->d1 = _cy_i03[ithread];
			cy_r08->d0 = _cy_r04[ithread];	cy_r08->d1 = _cy_i04[ithread];
			cy_r10->d0 = _cy_r05[ithread];	cy_r10->d1 = _cy_i05[ithread];
			cy_r12->d0 = _cy_r06[ithread];	cy_r12->d1 = _cy_i06[ithread];
			cy_r14->d0 = _cy_r07[ithread];	cy_r14->d1 = _cy_i07[ithread];
			cy_r16->d0 = _cy_r08[ithread];	cy_r16->d1 = _cy_i08[ithread];
			cy_r18->d0 = _cy_r09[ithread];	cy_r18->d1 = _cy_i09[ithread];
			cy_r20->d0 = _cy_r10[ithread];	cy_r20->d1 = _cy_i10[ithread];
			cy_r22->d0 = _cy_r11[ithread];	cy_r22->d1 = _cy_i11[ithread];
			cy_r24->d0 = _cy_r12[ithread];	cy_r24->d1 = _cy_i12[ithread];
			cy_r26->d0 = _cy_r13[ithread];	cy_r26->d1 = _cy_i13[ithread];
			cy_r28->d0 = _cy_r14[ithread];	cy_r28->d1 = _cy_i14[ithread];
			cy_r30->d0 = _cy_r15[ithread];	cy_r30->d1 = _cy_i15[ithread];
			cy_r32->d0 = _cy_r16[ithread];	cy_r32->d1 = _cy_i16[ithread];
			cy_r34->d0 = _cy_r17[ithread];	cy_r34->d1 = _cy_i17[ithread];
			cy_r36->d0 = _cy_r18[ithread];	cy_r36->d1 = _cy_i18[ithread];
			cy_r38->d0 = _cy_r19[ithread];	cy_r38->d1 = _cy_i19[ithread];
			cy_r40->d0 = _cy_r20[ithread];	cy_r40->d1 = _cy_i20[ithread];
			cy_r42->d0 = _cy_r21[ithread];	cy_r42->d1 = _cy_i21[ithread];
			cy_r44->d0 = _cy_r22[ithread];	cy_r44->d1 = _cy_i22[ithread];
			cy_r46->d0 = _cy_r23[ithread];	cy_r46->d1 = _cy_i23[ithread];
			cy_r48->d0 = _cy_r24[ithread];	cy_r48->d1 = _cy_i24[ithread];
			cy_r50->d0 = _cy_r25[ithread];	cy_r50->d1 = _cy_i25[ithread];
			cy_r52->d0 = _cy_r26[ithread];	cy_r52->d1 = _cy_i26[ithread];
			cy_r54->d0 = _cy_r27[ithread];	cy_r54->d1 = _cy_i27[ithread];
			cy_r56->d0 = _cy_r28[ithread];	cy_r56->d1 = _cy_i28[ithread];
			cy_r58->d0 = _cy_r29[ithread];	cy_r58->d1 = _cy_i29[ithread];
			cy_i00->d0 = _cy_r30[ithread];	cy_i00->d1 = _cy_i30[ithread];
			cy_i02->d0 = _cy_r31[ithread];	cy_i02->d1 = _cy_i31[ithread];
			cy_i04->d0 = _cy_r32[ithread];	cy_i04->d1 = _cy_i32[ithread];
			cy_i06->d0 = _cy_r33[ithread];	cy_i06->d1 = _cy_i33[ithread];
			cy_i08->d0 = _cy_r34[ithread];	cy_i08->d1 = _cy_i34[ithread];
			cy_i10->d0 = _cy_r35[ithread];	cy_i10->d1 = _cy_i35[ithread];
			cy_i12->d0 = _cy_r36[ithread];	cy_i12->d1 = _cy_i36[ithread];
			cy_i14->d0 = _cy_r37[ithread];	cy_i14->d1 = _cy_i37[ithread];
			cy_i16->d0 = _cy_r38[ithread];	cy_i16->d1 = _cy_i38[ithread];
			cy_i18->d0 = _cy_r39[ithread];	cy_i18->d1 = _cy_i39[ithread];
			cy_i20->d0 = _cy_r40[ithread];	cy_i20->d1 = _cy_i40[ithread];
			cy_i22->d0 = _cy_r41[ithread];	cy_i22->d1 = _cy_i41[ithread];
			cy_i24->d0 = _cy_r42[ithread];	cy_i24->d1 = _cy_i42[ithread];
			cy_i26->d0 = _cy_r43[ithread];	cy_i26->d1 = _cy_i43[ithread];
			cy_i28->d0 = _cy_r44[ithread];	cy_i28->d1 = _cy_i44[ithread];
			cy_i30->d0 = _cy_r45[ithread];	cy_i30->d1 = _cy_i45[ithread];
			cy_i32->d0 = _cy_r46[ithread];	cy_i32->d1 = _cy_i46[ithread];
			cy_i34->d0 = _cy_r47[ithread];	cy_i34->d1 = _cy_i47[ithread];
			cy_i36->d0 = _cy_r48[ithread];	cy_i36->d1 = _cy_i48[ithread];
			cy_i38->d0 = _cy_r49[ithread];	cy_i38->d1 = _cy_i49[ithread];
			cy_i40->d0 = _cy_r50[ithread];	cy_i40->d1 = _cy_i50[ithread];
			cy_i42->d0 = _cy_r51[ithread];	cy_i42->d1 = _cy_i51[ithread];
			cy_i44->d0 = _cy_r52[ithread];	cy_i44->d1 = _cy_i52[ithread];
			cy_i46->d0 = _cy_r53[ithread];	cy_i46->d1 = _cy_i53[ithread];
			cy_i48->d0 = _cy_r54[ithread];	cy_i48->d1 = _cy_i54[ithread];
			cy_i50->d0 = _cy_r55[ithread];	cy_i50->d1 = _cy_i55[ithread];
			cy_i52->d0 = _cy_r56[ithread];	cy_i52->d1 = _cy_i56[ithread];
			cy_i54->d0 = _cy_r57[ithread];	cy_i54->d1 = _cy_i57[ithread];
			cy_i56->d0 = _cy_r58[ithread];	cy_i56->d1 = _cy_i58[ithread];
			cy_i58->d0 = _cy_r59[ithread];	cy_i58->d1 = _cy_i59[ithread];
		#else
			/* init carries	*/
			cy_r00 = _cy_r00[ithread];	cy_i00 = _cy_i00[ithread];
			cy_r01 = _cy_r01[ithread];	cy_i01 = _cy_i01[ithread];
			cy_r02 = _cy_r02[ithread];	cy_i02 = _cy_i02[ithread];
			cy_r03 = _cy_r03[ithread];	cy_i03 = _cy_i03[ithread];
			cy_r04 = _cy_r04[ithread];	cy_i04 = _cy_i04[ithread];
			cy_r05 = _cy_r05[ithread];	cy_i05 = _cy_i05[ithread];
			cy_r06 = _cy_r06[ithread];	cy_i06 = _cy_i06[ithread];
			cy_r07 = _cy_r07[ithread];	cy_i07 = _cy_i07[ithread];
			cy_r08 = _cy_r08[ithread];	cy_i08 = _cy_i08[ithread];
			cy_r09 = _cy_r09[ithread];	cy_i09 = _cy_i09[ithread];
			cy_r10 = _cy_r10[ithread];	cy_i10 = _cy_i10[ithread];
			cy_r11 = _cy_r11[ithread];	cy_i11 = _cy_i11[ithread];
			cy_r12 = _cy_r12[ithread];	cy_i12 = _cy_i12[ithread];
			cy_r13 = _cy_r13[ithread];	cy_i13 = _cy_i13[ithread];
			cy_r14 = _cy_r14[ithread];	cy_i14 = _cy_i14[ithread];
			cy_r15 = _cy_r15[ithread];	cy_i15 = _cy_i15[ithread];
			cy_r16 = _cy_r16[ithread];	cy_i16 = _cy_i16[ithread];
			cy_r17 = _cy_r17[ithread];	cy_i17 = _cy_i17[ithread];
			cy_r18 = _cy_r18[ithread];	cy_i18 = _cy_i18[ithread];
			cy_r19 = _cy_r19[ithread];	cy_i19 = _cy_i19[ithread];
			cy_r20 = _cy_r20[ithread];	cy_i20 = _cy_i20[ithread];
			cy_r21 = _cy_r21[ithread];	cy_i21 = _cy_i21[ithread];
			cy_r22 = _cy_r22[ithread];	cy_i22 = _cy_i22[ithread];
			cy_r23 = _cy_r23[ithread];	cy_i23 = _cy_i23[ithread];
			cy_r24 = _cy_r24[ithread];	cy_i24 = _cy_i24[ithread];
			cy_r25 = _cy_r25[ithread];	cy_i25 = _cy_i25[ithread];
			cy_r26 = _cy_r26[ithread];	cy_i26 = _cy_i26[ithread];
			cy_r27 = _cy_r27[ithread];	cy_i27 = _cy_i27[ithread];
			cy_r28 = _cy_r28[ithread];	cy_i28 = _cy_i28[ithread];
			cy_r29 = _cy_r29[ithread];	cy_i29 = _cy_i29[ithread];
			cy_r30 = _cy_r30[ithread];	cy_i30 = _cy_i30[ithread];
			cy_r31 = _cy_r31[ithread];	cy_i31 = _cy_i31[ithread];
			cy_r32 = _cy_r32[ithread];	cy_i32 = _cy_i32[ithread];
			cy_r33 = _cy_r33[ithread];	cy_i33 = _cy_i33[ithread];
			cy_r34 = _cy_r34[ithread];	cy_i34 = _cy_i34[ithread];
			cy_r35 = _cy_r35[ithread];	cy_i35 = _cy_i35[ithread];
			cy_r36 = _cy_r36[ithread];	cy_i36 = _cy_i36[ithread];
			cy_r37 = _cy_r37[ithread];	cy_i37 = _cy_i37[ithread];
			cy_r38 = _cy_r38[ithread];	cy_i38 = _cy_i38[ithread];
			cy_r39 = _cy_r39[ithread];	cy_i39 = _cy_i39[ithread];
			cy_r40 = _cy_r40[ithread];	cy_i40 = _cy_i40[ithread];
			cy_r41 = _cy_r41[ithread];	cy_i41 = _cy_i41[ithread];
			cy_r42 = _cy_r42[ithread];	cy_i42 = _cy_i42[ithread];
			cy_r43 = _cy_r43[ithread];	cy_i43 = _cy_i43[ithread];
			cy_r44 = _cy_r44[ithread];	cy_i44 = _cy_i44[ithread];
			cy_r45 = _cy_r45[ithread];	cy_i45 = _cy_i45[ithread];
			cy_r46 = _cy_r46[ithread];	cy_i46 = _cy_i46[ithread];
			cy_r47 = _cy_r47[ithread];	cy_i47 = _cy_i47[ithread];
			cy_r48 = _cy_r48[ithread];	cy_i48 = _cy_i48[ithread];
			cy_r49 = _cy_r49[ithread];	cy_i49 = _cy_i49[ithread];
			cy_r50 = _cy_r50[ithread];	cy_i50 = _cy_i50[ithread];
			cy_r51 = _cy_r51[ithread];	cy_i51 = _cy_i51[ithread];
			cy_r52 = _cy_r52[ithread];	cy_i52 = _cy_i52[ithread];
			cy_r53 = _cy_r53[ithread];	cy_i53 = _cy_i53[ithread];
			cy_r54 = _cy_r54[ithread];	cy_i54 = _cy_i54[ithread];
			cy_r55 = _cy_r55[ithread];	cy_i55 = _cy_i55[ithread];
			cy_r56 = _cy_r56[ithread];	cy_i56 = _cy_i56[ithread];
			cy_r57 = _cy_r57[ithread];	cy_i57 = _cy_i57[ithread];
			cy_r58 = _cy_r58[ithread];	cy_i58 = _cy_i58[ithread];
			cy_r59 = _cy_r59[ithread];	cy_i59 = _cy_i59[ithread];
		#endif
		}

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
			/* In SIMD mode, data are arranged in [re_0,...,re_n-1,im_0,...,im_n-1] groups, not the usual [re_0,im_0],...,[re_n-1,im_n-1] pairs.
			Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
			but to point to the proper real datum, we need to index-map e.g. [0,1,2,3] ==> [0,2,1,3] in 2-way SIMD mode.
			(But only ever need to explicitly do this in debug mode).
			*/
			for(j = jstart; j < jhi; j += stride)
			{
				j1 =  j;
				j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1 + RE_IM_STRIDE;

		/*...The radix-60 DIT pass is here:	*/

			/* The 15 radix-4 subblocks have several different SIMD flow modes: */
			#ifdef USE_AVX
			/* Outputs in AVX mode are temps 2*15*32 = 30*32 = 0x3c0 bytes apart: */
			// Reorder blocks to yield sequentially increasing a-array offsets:
				add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, 0x3c0)
				add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, 0x3c0)
				add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01, 0x3c0)
				add2 = &a[j1+p12];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05, 0x3c0)
				add1 = &a[j1+p16];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, 0x3c0)
				add3 = &a[j1+p20];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03, 0x3c0)
				add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07, 0x3c0)
				add3 = &a[j1+p28];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, 0x3c0)
				add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, 0x3c0)
				add3 = &a[j1+p36];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09, 0x3c0)
				add0 = &a[j1+p40];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0b, 0x3c0)
				add2 = &a[j1+p44];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, 0x3c0)
				add0 = &a[j1+p48];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, 0x3c0)
				add2 = &a[j1+p52];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0d, 0x3c0)
				add1 = &a[j1+p56];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, 0x3c0)

			#elif defined(USE_SSE2)

			/* Outputs in SSE2 mode are temps 2*15*16 = 30*16 = 0x1e0 bytes apart: */
			  #if !GCC_ASM_FULL_INLINE
			// Reorder blocks to yield sequentially increasing a-array offsets:
/* Block 01: */	add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, 0x1e0)
/* Block 03: */	add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, 0x1e0)
/* Block 02: */	add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01, 0x1e0)
/* Block 06: */	add2 = &a[j1+p12];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05, 0x1e0)
/* Block 05: */	add1 = &a[j1+p16];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, 0x1e0)
/* Block 04: */	add3 = &a[j1+p20];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03, 0x1e0)
/* Block 08: */	add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07, 0x1e0)
/* Block 07: */	add3 = &a[j1+p28];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, 0x1e0)
/* Block 09: */	add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, 0x1e0)
/* Block 10: */	add3 = &a[j1+p36];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09, 0x1e0)
/* Block 12: */	add0 = &a[j1+p40];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0b, 0x1e0)
/* Block 11: */	add2 = &a[j1+p44];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, 0x1e0)
/* Block 15: */	add0 = &a[j1+p48];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, 0x1e0)
/* Block 14: */	add2 = &a[j1+p52];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0d, 0x1e0)
/* Block 13: */	add1 = &a[j1+p56];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, 0x1e0)

			  #else

				// The macro name here is misleading (at present): Macro just does 15 x radix-4 DFT:
				add0 = &a[j1    ];
				SSE2_RADIX60_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,r00);

			  #endif

			#endif	// AVX or SSE2?

			#ifdef USE_SSE2	// USE_AVX under same umbrella here

			  #if (OS_BITS == 32) || !defined(USE_64BIT_ASM_STYLE)	// In 64-bit mode, default is to use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just xmm0-7.
				/* Radix-15 DFT uses adjacent temps for inputs, outputs are (cyclic) with pXXr having XX += 12 (24*16 bytes = 0x180) between successive outputs: */
			   #ifndef USE_LITERAL_BYTE_OFFSETS
				SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1     ,  r00,  r01,  r02,  r03,  r04,  r05,  r06,  r07,  r08,  r09,  r0a,  r0b,  r0c,  r0d,  r0e     ,  x00,  x01,  x02,  x03,  x04,  x05,  x06,  x07,  x08,  x09,  x0a,  x0b,  x0c,  x0d,  x0e,        s1p00r,s1p2er,s1p1dr,s1p0cr,s1p3br,s1p2ar,s1p19r,s1p08r,s1p37r,s1p26r,s1p15r,s1p04r,s1p33r,s1p22r,s1p11r);
				SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1     ,  r10,  r11,  r12,  r13,  r14,  r15,  r16,  r17,  r18,  r19,  r1a,  r1b,  r1c,  r1d,  r1e     ,  x00,  x01,  x02,  x03,  x04,  x05,  x06,  x07,  x08,  x09,  x0a,  x0b,  x0c,  x0d,  x0e,        s1p10r,s1p3er,s1p2dr,s1p1cr,s1p0br,s1p3ar,s1p29r,s1p18r,s1p07r,s1p36r,s1p25r,s1p14r,s1p03r,s1p32r,s1p21r);
				SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1     ,  r20,  r21,  r22,  r23,  r24,  r25,  r26,  r27,  r28,  r29,  r2a,  r2b,  r2c,  r2d,  r2e     ,  x00,  x01,  x02,  x03,  x04,  x05,  x06,  x07,  x08,  x09,  x0a,  x0b,  x0c,  x0d,  x0e,        s1p20r,s1p0er,s1p3dr,s1p2cr,s1p1br,s1p0ar,s1p39r,s1p28r,s1p17r,s1p06r,s1p35r,s1p24r,s1p13r,s1p02r,s1p31r);
				SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1     ,  r30,  r31,  r32,  r33,  r34,  r35,  r36,  r37,  r38,  r39,  r3a,  r3b,  r3c,  r3d,  r3e     ,  x00,  x01,  x02,  x03,  x04,  x05,  x06,  x07,  x08,  x09,  x0a,  x0b,  x0c,  x0d,  x0e,        s1p30r,s1p1er,s1p0dr,s1p3cr,s1p2br,s1p1ar,s1p09r,s1p38r,s1p27r,s1p16r,s1p05r,s1p34r,s1p23r,s1p12r,s1p01r);
			   #else	// Versions using base-address-plus-literal-byte-offsets:                                                                                       s1p**r complex-element offsets w.r.to s1p00r (all negative offsets get added to +120):     0    +88    +56    +24    -8     -40    -72    -104   -16    -48    -80    -112   -24    -56    -88 ... when offset < -(row-start offset below), reset by adding 120*16 = 0x780 bytes
				SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1, r00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, s1p00r, 0x000, 0x580, 0x380, 0x180, 0x700, 0x500, 0x300, 0x100, 0x680, 0x480, 0x280, 0x080, 0x600, 0x400, 0x200);	// end offset of each s1p-sequence = start-offset + 0x200
				SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1, r10,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, s1p00r, 0x1e0, 0x760, 0x560, 0x360, 0x160, 0x6e0, 0x4e0, 0x2e0, 0x0e0, 0x660, 0x460, 0x260, 0x060, 0x5e0, 0x3e0);	// s1p10r = +30 complex = +0x1e0 bytes
				SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1, r20,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, s1p00r, 0x3c0, 0x1c0, 0x740, 0x540, 0x340, 0x140, 0x6c0, 0x4c0, 0x2c0, 0x0c0, 0x640, 0x440, 0x240, 0x040, 0x5c0);	// s1p20r = +60 complex = +0x3c0 bytes
				SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1, r30,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, s1p00r, 0x5a0, 0x3a0, 0x1a0, 0x720, 0x520, 0x320, 0x120, 0x6a0, 0x4a0, 0x2a0, 0x0a0, 0x620, 0x420, 0x220, 0x020);	// s1p30r = +90 complex = +0x5a0 bytes
			   #endif
			  #else
				SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,\
					r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r0a,r0b,r0c,r0d,r0e,\
					x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,\
					s1p00r,s1p2er,s1p1dr,s1p0cr,s1p3br,s1p2ar,s1p19r,s1p08r,s1p37r,s1p26r,s1p15r,s1p04r,s1p33r,s1p22r,s1p11r,\
					r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r1a,r1b,r1c,r1d,r1e,\
					y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,\
					s1p10r,s1p3er,s1p2dr,s1p1cr,s1p0br,s1p3ar,s1p29r,s1p18r,s1p07r,s1p36r,s1p25r,s1p14r,s1p03r,s1p32r,s1p21r);
				SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,\
					r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r2a,r2b,r2c,r2d,r2e,\
					x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,\
					s1p20r,s1p0er,s1p3dr,s1p2cr,s1p1br,s1p0ar,s1p39r,s1p28r,s1p17r,s1p06r,s1p35r,s1p24r,s1p13r,s1p02r,s1p31r,\
					r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r3a,r3b,r3c,r3d,r3e,\
					y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,\
					s1p30r,s1p1er,s1p0dr,s1p3cr,s1p2br,s1p1ar,s1p09r,s1p38r,s1p27r,s1p16r,s1p05r,s1p34r,s1p23r,s1p12r,s1p01r);
			  #endif

			#else	// Scalar-double mode:

				jt = j1    ; jp = j2    ;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,rt,it);
				jt = j1+p08; jp = j2+p08;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,rt,it);
				jt = j1+p04; jp = j2+p04;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,rt,it);
				jt = j1+p20; jp = j2+p20;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,rt,it);
				jt = j1+p16; jp = j2+p16;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,rt,it);
				jt = j1+p12; jp = j2+p12;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,rt,it);
				jt = j1+p28; jp = j2+p28;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,rt,it);
				jt = j1+p24; jp = j2+p24;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,rt,it);
				jt = j1+p32; jp = j2+p32;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,rt,it);
				jt = j1+p36; jp = j2+p36;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,rt,it);
				jt = j1+p44; jp = j2+p44;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,rt,it);
				jt = j1+p40; jp = j2+p40;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0br,t0bi,t1br,t1bi,t2br,t2bi,t3br,t3bi,rt,it);
				jt = j1+p56; jp = j2+p56;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t0cr,t0ci,t1cr,t1ci,t2cr,t2ci,t3cr,t3ci,rt,it);
				jt = j1+p52; jp = j2+p52;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t0dr,t0di,t1dr,t1di,t2dr,t2di,t3dr,t3di,rt,it);
				jt = j1+p48; jp = j2+p48;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0er,t0ei,t1er,t1ei,t2er,t2ei,t3er,t3ei,rt,it);

				/*...and now do 4 radix-15 transforms.
				The required output permutation is [in terms of decimal indices...we use base-15 indexing in the actual code below]

					[00,44,28,12,56,40,24,08,52,36,20,04,48,32,16
					 15,59,43,27,11,55,39,23,07,51,35,19,03,47,31
					 30,14,58,42,26,10,54,38,22,06,50,34,18,02,46
					 45,29,13,57,41,25,09,53,37,21,05,49,33,17,01]
				*/
				RADIX_15_DIT(t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a1p00r,a1p00i,a1p2er,a1p2ei,a1p1dr,a1p1di,a1p0cr,a1p0ci,a1p3br,a1p3bi,a1p2ar,a1p2ai,a1p19r,a1p19i,a1p08r,a1p08i,a1p37r,a1p37i,a1p26r,a1p26i,a1p15r,a1p15i,a1p04r,a1p04i,a1p33r,a1p33i,a1p22r,a1p22i,a1p11r,a1p11i,rt,it);
				RADIX_15_DIT(t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a1p10r,a1p10i,a1p3er,a1p3ei,a1p2dr,a1p2di,a1p1cr,a1p1ci,a1p0br,a1p0bi,a1p3ar,a1p3ai,a1p29r,a1p29i,a1p18r,a1p18i,a1p07r,a1p07i,a1p36r,a1p36i,a1p25r,a1p25i,a1p14r,a1p14i,a1p03r,a1p03i,a1p32r,a1p32i,a1p21r,a1p21i,rt,it);
				RADIX_15_DIT(t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a1p20r,a1p20i,a1p0er,a1p0ei,a1p3dr,a1p3di,a1p2cr,a1p2ci,a1p1br,a1p1bi,a1p0ar,a1p0ai,a1p39r,a1p39i,a1p28r,a1p28i,a1p17r,a1p17i,a1p06r,a1p06i,a1p35r,a1p35i,a1p24r,a1p24i,a1p13r,a1p13i,a1p02r,a1p02i,a1p31r,a1p31i,rt,it);
				RADIX_15_DIT(t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a1p30r,a1p30i,a1p1er,a1p1ei,a1p0dr,a1p0di,a1p3cr,a1p3ci,a1p2br,a1p2bi,a1p1ar,a1p1ai,a1p09r,a1p09i,a1p38r,a1p38i,a1p27r,a1p27i,a1p16r,a1p16i,a1p05r,a1p05i,a1p34r,a1p34i,a1p23r,a1p23i,a1p12r,a1p12i,a1p01r,a1p01i,rt,it);

			#endif	// SIMD or not?

		/*...Now do the carries. Since the outputs would
			normally be getting dispatched to 60 separate blocks of the A-array, we need 60 separate carries.	*/

		/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/
			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
			#ifdef USE_AVX

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

				l= j & (nwt-1);						tmp = half_arr + 64;	/* ptr to local storage for the doubled wtl,wtn terms: */
				n_minus_sil  ->d0 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d0 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d0 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d0 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+2) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d1 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d1 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d1 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d1 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+4) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d2 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d2 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d2 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d2 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+6) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d3 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d3 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d3 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d3 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				AVX_cmplx_carry_norm_errcheck0_X4(s1p00r,add1,add2,add3,cy_r00,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p04r,add1,add2,add3,cy_r04,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p08r,add1,add2,add3,cy_r08,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p0cr,add1,add2,add3,cy_r12,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p11r,add1,add2,add3,cy_r16,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p15r,add1,add2,add3,cy_r20,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p19r,add1,add2,add3,cy_r24,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p1dr,add1,add2,add3,cy_r28,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p22r,add1,add2,add3,cy_r32,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p26r,add1,add2,add3,cy_r36,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p2ar,add1,add2,add3,cy_r40,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p2er,add1,add2,add3,cy_r44,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p33r,add1,add2,add3,cy_r48,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p37r,add1,add2,add3,cy_r52,bjmodn52,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p3br,add1,add2,add3,cy_r56,bjmodn56,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);

				co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

			#elif defined(USE_SSE2)

				l= j & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
				ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
				ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
				ctmp->re = wtnm1;	ctmp->im = wtnm1;

				add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

			  #if defined(COMPILER_TYPE_MSVC)

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p0cr,add1,add2,add3,cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p11r,add1,add2,add3,cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p15r,add1,add2,add3,cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p19r,add1,add2,add3,cy_r24,cy_r26,bjmodn24);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p1dr,add1,add2,add3,cy_r28,cy_r30,bjmodn28);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p22r,add1,add2,add3,cy_r32,cy_r34,bjmodn32);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p26r,add1,add2,add3,cy_r36,cy_r38,bjmodn36);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p2ar,add1,add2,add3,cy_r40,cy_r42,bjmodn40);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p2er,add1,add2,add3,cy_r44,cy_r46,bjmodn44);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p33r,add1,add2,add3,cy_r48,cy_r50,bjmodn48);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p37r,add1,add2,add3,cy_r52,cy_r54,bjmodn52);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p3br,add1,add2,add3,cy_r56,cy_r58,bjmodn56);
			   #else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p0cr,add1,add2,add3,cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p11r,add1,add2,add3,cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p15r,add1,add2,add3,cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p19r,add1,add2,add3,cy_r24,cy_r26,bjmodn24);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p1dr,add1,add2,add3,cy_r28,cy_r30,bjmodn28);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p22r,add1,add2,add3,cy_r32,cy_r34,bjmodn32);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p26r,add1,add2,add3,cy_r36,cy_r38,bjmodn36);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p2ar,add1,add2,add3,cy_r40,cy_r42,bjmodn40);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p2er,add1,add2,add3,cy_r44,cy_r46,bjmodn44);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p33r,add1,add2,add3,cy_r48,cy_r50,bjmodn48);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p37r,add1,add2,add3,cy_r52,cy_r54,bjmodn52);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p3br,add1,add2,add3,cy_r56,cy_r58,bjmodn56);
			   #endif

			  #else	/* GCC-style inline ASM: */

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p0cr,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p11r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p15r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p19r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p1dr,add1,add2,add3,cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p22r,add1,add2,add3,cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p26r,add1,add2,add3,cy_r36,cy_r38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p2ar,add1,add2,add3,cy_r40,cy_r42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p2er,add1,add2,add3,cy_r44,cy_r46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p33r,add1,add2,add3,cy_r48,cy_r50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p37r,add1,add2,add3,cy_r52,cy_r54,bjmodn52,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p3br,add1,add2,add3,cy_r56,cy_r58,bjmodn56,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			   #else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p0cr,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p11r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p15r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p19r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p1dr,add1,add2,add3,cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p22r,add1,add2,add3,cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p26r,add1,add2,add3,cy_r36,cy_r38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p2ar,add1,add2,add3,cy_r40,cy_r42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p2er,add1,add2,add3,cy_r44,cy_r46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p33r,add1,add2,add3,cy_r48,cy_r50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p37r,add1,add2,add3,cy_r52,cy_r54,bjmodn52,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p3br,add1,add2,add3,cy_r56,cy_r58,bjmodn56,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			   #endif

			  #endif

				l= (j+2) & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
				ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
				ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
				ctmp->re = wtnm1;	ctmp->im = wtnm1;

			/*	i =((uint32)(sw - *bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];

			  #if defined(COMPILER_TYPE_MSVC)

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p0cr,add1,add2,     cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p11r,add1,add2,     cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p15r,add1,add2,     cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p19r,add1,add2,     cy_r24,cy_r26,bjmodn24);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p1dr,add1,add2,     cy_r28,cy_r30,bjmodn28);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p22r,add1,add2,     cy_r32,cy_r34,bjmodn32);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p26r,add1,add2,     cy_r36,cy_r38,bjmodn36);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p2ar,add1,add2,     cy_r40,cy_r42,bjmodn40);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p2er,add1,add2,     cy_r44,cy_r46,bjmodn44);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p33r,add1,add2,     cy_r48,cy_r50,bjmodn48);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p37r,add1,add2,     cy_r52,cy_r54,bjmodn52);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p3br,add1,add2,     cy_r56,cy_r58,bjmodn56);
			   #else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p0cr,add1,add2,     cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p11r,add1,add2,     cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p15r,add1,add2,     cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p19r,add1,add2,     cy_r24,cy_r26,bjmodn24);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p1dr,add1,add2,     cy_r28,cy_r30,bjmodn28);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p22r,add1,add2,     cy_r32,cy_r34,bjmodn32);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p26r,add1,add2,     cy_r36,cy_r38,bjmodn36);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p2ar,add1,add2,     cy_r40,cy_r42,bjmodn40);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p2er,add1,add2,     cy_r44,cy_r46,bjmodn44);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p33r,add1,add2,     cy_r48,cy_r50,bjmodn48);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p37r,add1,add2,     cy_r52,cy_r54,bjmodn52);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p3br,add1,add2,     cy_r56,cy_r58,bjmodn56);
			   #endif

			  #else	/* GCC-style inline ASM: */

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p0cr,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p11r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p15r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p19r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p1dr,add1,add2,     cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p22r,add1,add2,     cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p26r,add1,add2,     cy_r36,cy_r38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p2ar,add1,add2,     cy_r40,cy_r42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p2er,add1,add2,     cy_r44,cy_r46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p33r,add1,add2,     cy_r48,cy_r50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p37r,add1,add2,     cy_r52,cy_r54,bjmodn52,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p3br,add1,add2,     cy_r56,cy_r58,bjmodn56,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			   #else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p0cr,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p11r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p15r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p19r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p1dr,add1,add2,     cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p22r,add1,add2,     cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p26r,add1,add2,     cy_r36,cy_r38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p2ar,add1,add2,     cy_r40,cy_r42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p2er,add1,add2,     cy_r44,cy_r46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p33r,add1,add2,     cy_r48,cy_r50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p37r,add1,add2,     cy_r52,cy_r54,bjmodn52,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p3br,add1,add2,     cy_r56,cy_r58,bjmodn56,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			   #endif

			  #endif

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

			#else	// Scalar-double mode:

				l= j & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

			/*...set0 is slightly different from others:	*/
			   cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy_r00,bjmodn00   );
				cmplx_carry_norm_errcheck(a1p01r,a1p01i,cy_r01,bjmodn01, 1);
				cmplx_carry_norm_errcheck(a1p02r,a1p02i,cy_r02,bjmodn02, 2);
				cmplx_carry_norm_errcheck(a1p03r,a1p03i,cy_r03,bjmodn03, 3);
				cmplx_carry_norm_errcheck(a1p04r,a1p04i,cy_r04,bjmodn04, 4);
				cmplx_carry_norm_errcheck(a1p05r,a1p05i,cy_r05,bjmodn05, 5);
				cmplx_carry_norm_errcheck(a1p06r,a1p06i,cy_r06,bjmodn06, 6);
				cmplx_carry_norm_errcheck(a1p07r,a1p07i,cy_r07,bjmodn07, 7);
				cmplx_carry_norm_errcheck(a1p08r,a1p08i,cy_r08,bjmodn08, 8);
				cmplx_carry_norm_errcheck(a1p09r,a1p09i,cy_r09,bjmodn09, 9);
				cmplx_carry_norm_errcheck(a1p0ar,a1p0ai,cy_r10,bjmodn10,10);
				cmplx_carry_norm_errcheck(a1p0br,a1p0bi,cy_r11,bjmodn11,11);
				cmplx_carry_norm_errcheck(a1p0cr,a1p0ci,cy_r12,bjmodn12,12);
				cmplx_carry_norm_errcheck(a1p0dr,a1p0di,cy_r13,bjmodn13,13);
				cmplx_carry_norm_errcheck(a1p0er,a1p0ei,cy_r14,bjmodn14,14);
				cmplx_carry_norm_errcheck(a1p10r,a1p10i,cy_r15,bjmodn15,15);
				cmplx_carry_norm_errcheck(a1p11r,a1p11i,cy_r16,bjmodn16,16);
				cmplx_carry_norm_errcheck(a1p12r,a1p12i,cy_r17,bjmodn17,17);
				cmplx_carry_norm_errcheck(a1p13r,a1p13i,cy_r18,bjmodn18,18);
				cmplx_carry_norm_errcheck(a1p14r,a1p14i,cy_r19,bjmodn19,19);
				cmplx_carry_norm_errcheck(a1p15r,a1p15i,cy_r20,bjmodn20,20);
				cmplx_carry_norm_errcheck(a1p16r,a1p16i,cy_r21,bjmodn21,21);
				cmplx_carry_norm_errcheck(a1p17r,a1p17i,cy_r22,bjmodn22,22);
				cmplx_carry_norm_errcheck(a1p18r,a1p18i,cy_r23,bjmodn23,23);
				cmplx_carry_norm_errcheck(a1p19r,a1p19i,cy_r24,bjmodn24,24);
				cmplx_carry_norm_errcheck(a1p1ar,a1p1ai,cy_r25,bjmodn25,25);
				cmplx_carry_norm_errcheck(a1p1br,a1p1bi,cy_r26,bjmodn26,26);
				cmplx_carry_norm_errcheck(a1p1cr,a1p1ci,cy_r27,bjmodn27,27);
				cmplx_carry_norm_errcheck(a1p1dr,a1p1di,cy_r28,bjmodn28,28);
				cmplx_carry_norm_errcheck(a1p1er,a1p1ei,cy_r29,bjmodn29,29);
				cmplx_carry_norm_errcheck(a1p20r,a1p20i,cy_r30,bjmodn30,30);
				cmplx_carry_norm_errcheck(a1p21r,a1p21i,cy_r31,bjmodn31,31);
				cmplx_carry_norm_errcheck(a1p22r,a1p22i,cy_r32,bjmodn32,32);
				cmplx_carry_norm_errcheck(a1p23r,a1p23i,cy_r33,bjmodn33,33);
				cmplx_carry_norm_errcheck(a1p24r,a1p24i,cy_r34,bjmodn34,34);
				cmplx_carry_norm_errcheck(a1p25r,a1p25i,cy_r35,bjmodn35,35);
				cmplx_carry_norm_errcheck(a1p26r,a1p26i,cy_r36,bjmodn36,36);
				cmplx_carry_norm_errcheck(a1p27r,a1p27i,cy_r37,bjmodn37,37);
				cmplx_carry_norm_errcheck(a1p28r,a1p28i,cy_r38,bjmodn38,38);
				cmplx_carry_norm_errcheck(a1p29r,a1p29i,cy_r39,bjmodn39,39);
				cmplx_carry_norm_errcheck(a1p2ar,a1p2ai,cy_r40,bjmodn40,40);
				cmplx_carry_norm_errcheck(a1p2br,a1p2bi,cy_r41,bjmodn41,41);
				cmplx_carry_norm_errcheck(a1p2cr,a1p2ci,cy_r42,bjmodn42,42);
				cmplx_carry_norm_errcheck(a1p2dr,a1p2di,cy_r43,bjmodn43,43);
				cmplx_carry_norm_errcheck(a1p2er,a1p2ei,cy_r44,bjmodn44,44);
				cmplx_carry_norm_errcheck(a1p30r,a1p30i,cy_r45,bjmodn45,45);
				cmplx_carry_norm_errcheck(a1p31r,a1p31i,cy_r46,bjmodn46,46);
				cmplx_carry_norm_errcheck(a1p32r,a1p32i,cy_r47,bjmodn47,47);
				cmplx_carry_norm_errcheck(a1p33r,a1p33i,cy_r48,bjmodn48,48);
				cmplx_carry_norm_errcheck(a1p34r,a1p34i,cy_r49,bjmodn49,49);
				cmplx_carry_norm_errcheck(a1p35r,a1p35i,cy_r50,bjmodn50,50);
				cmplx_carry_norm_errcheck(a1p36r,a1p36i,cy_r51,bjmodn51,51);
				cmplx_carry_norm_errcheck(a1p37r,a1p37i,cy_r52,bjmodn52,52);
				cmplx_carry_norm_errcheck(a1p38r,a1p38i,cy_r53,bjmodn53,53);
				cmplx_carry_norm_errcheck(a1p39r,a1p39i,cy_r54,bjmodn54,54);
				cmplx_carry_norm_errcheck(a1p3ar,a1p3ai,cy_r55,bjmodn55,55);
				cmplx_carry_norm_errcheck(a1p3br,a1p3bi,cy_r56,bjmodn56,56);
				cmplx_carry_norm_errcheck(a1p3cr,a1p3ci,cy_r57,bjmodn57,57);
				cmplx_carry_norm_errcheck(a1p3dr,a1p3di,cy_r58,bjmodn58,58);
				cmplx_carry_norm_errcheck(a1p3er,a1p3ei,cy_r59,bjmodn59,59);

				i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

			#endif	/* #ifdef USE_SSE2 */

			}
			else	/* MODULUS_TYPE_FERMAT */
			{

			#ifdef USE_AVX

				// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.

				/* Get the needed Nth root of -1: */
				add1 = (double *)&rn0[0];
				add2 = (double *)&rn1[0];

				idx_offset = j;
				idx_incr = NDIVR;

				tmp = base_negacyclic_root;	tm2 = tmp+1;

			  #if HIACC
				// Hi-accuracy version needs 7 copies of each base root:
				l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				temp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				VEC_DBL_INIT(tmp+ 32,wt_re);	VEC_DBL_INIT(tm2+ 32,wt_im);
				VEC_DBL_INIT(tmp+ 40,wt_re);	VEC_DBL_INIT(tm2+ 40,wt_im);
				VEC_DBL_INIT(tmp+ 48,wt_re);	VEC_DBL_INIT(tm2+ 48,wt_im);
				VEC_DBL_INIT(tmp+ 56,wt_re);	VEC_DBL_INIT(tm2+ 56,wt_im);
				VEC_DBL_INIT(tmp+ 64,wt_re);	VEC_DBL_INIT(tm2+ 64,wt_im);
				VEC_DBL_INIT(tmp+ 72,wt_re);	VEC_DBL_INIT(tm2+ 72,wt_im);
				VEC_DBL_INIT(tmp+ 80,wt_re);	VEC_DBL_INIT(tm2+ 80,wt_im);
				VEC_DBL_INIT(tmp+ 88,wt_re);	VEC_DBL_INIT(tm2+ 88,wt_im);
				VEC_DBL_INIT(tmp+ 96,wt_re);	VEC_DBL_INIT(tm2+ 96,wt_im);
				VEC_DBL_INIT(tmp+104,wt_re);	VEC_DBL_INIT(tm2+104,wt_im);
				VEC_DBL_INIT(tmp+112,wt_re);	VEC_DBL_INIT(tm2+112,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				VEC_DBL_INIT(tmp+ 32,wt_re);	VEC_DBL_INIT(tm2+ 32,wt_im);
				VEC_DBL_INIT(tmp+ 40,wt_re);	VEC_DBL_INIT(tm2+ 40,wt_im);
				VEC_DBL_INIT(tmp+ 48,wt_re);	VEC_DBL_INIT(tm2+ 48,wt_im);
				VEC_DBL_INIT(tmp+ 56,wt_re);	VEC_DBL_INIT(tm2+ 56,wt_im);
				VEC_DBL_INIT(tmp+ 64,wt_re);	VEC_DBL_INIT(tm2+ 64,wt_im);
				VEC_DBL_INIT(tmp+ 72,wt_re);	VEC_DBL_INIT(tm2+ 72,wt_im);
				VEC_DBL_INIT(tmp+ 80,wt_re);	VEC_DBL_INIT(tm2+ 80,wt_im);
				VEC_DBL_INIT(tmp+ 88,wt_re);	VEC_DBL_INIT(tm2+ 88,wt_im);
				VEC_DBL_INIT(tmp+ 96,wt_re);	VEC_DBL_INIT(tm2+ 96,wt_im);
				VEC_DBL_INIT(tmp+104,wt_re);	VEC_DBL_INIT(tm2+104,wt_im);
				VEC_DBL_INIT(tmp+112,wt_re);	VEC_DBL_INIT(tm2+112,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				VEC_DBL_INIT(tmp+ 32,wt_re);	VEC_DBL_INIT(tm2+ 32,wt_im);
				VEC_DBL_INIT(tmp+ 40,wt_re);	VEC_DBL_INIT(tm2+ 40,wt_im);
				VEC_DBL_INIT(tmp+ 48,wt_re);	VEC_DBL_INIT(tm2+ 48,wt_im);
				VEC_DBL_INIT(tmp+ 56,wt_re);	VEC_DBL_INIT(tm2+ 56,wt_im);
				VEC_DBL_INIT(tmp+ 64,wt_re);	VEC_DBL_INIT(tm2+ 64,wt_im);
				VEC_DBL_INIT(tmp+ 72,wt_re);	VEC_DBL_INIT(tm2+ 72,wt_im);
				VEC_DBL_INIT(tmp+ 80,wt_re);	VEC_DBL_INIT(tm2+ 80,wt_im);
				VEC_DBL_INIT(tmp+ 88,wt_re);	VEC_DBL_INIT(tm2+ 88,wt_im);
				VEC_DBL_INIT(tmp+ 96,wt_re);	VEC_DBL_INIT(tm2+ 96,wt_im);
				VEC_DBL_INIT(tmp+104,wt_re);	VEC_DBL_INIT(tm2+104,wt_im);
				VEC_DBL_INIT(tmp+112,wt_re);	VEC_DBL_INIT(tm2+112,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				VEC_DBL_INIT(tmp+ 32,wt_re);	VEC_DBL_INIT(tm2+ 32,wt_im);
				VEC_DBL_INIT(tmp+ 40,wt_re);	VEC_DBL_INIT(tm2+ 40,wt_im);
				VEC_DBL_INIT(tmp+ 48,wt_re);	VEC_DBL_INIT(tm2+ 48,wt_im);
				VEC_DBL_INIT(tmp+ 56,wt_re);	VEC_DBL_INIT(tm2+ 56,wt_im);
				VEC_DBL_INIT(tmp+ 64,wt_re);	VEC_DBL_INIT(tm2+ 64,wt_im);
				VEC_DBL_INIT(tmp+ 72,wt_re);	VEC_DBL_INIT(tm2+ 72,wt_im);
				VEC_DBL_INIT(tmp+ 80,wt_re);	VEC_DBL_INIT(tm2+ 80,wt_im);
				VEC_DBL_INIT(tmp+ 88,wt_re);	VEC_DBL_INIT(tm2+ 88,wt_im);
				VEC_DBL_INIT(tmp+ 96,wt_re);	VEC_DBL_INIT(tm2+ 96,wt_im);
				VEC_DBL_INIT(tmp+104,wt_re);	VEC_DBL_INIT(tm2+104,wt_im);
				VEC_DBL_INIT(tmp+112,wt_re);	VEC_DBL_INIT(tm2+112,wt_im);

			  #else	// HIACC = false:

				// Get the needed quartet of Nth roots of -1: This is the same code as in the scalar
				// fermat_carry_norm_errcheck() macro, with the single index j replaced by the quartet j,j+2,j+4,j+6:
				l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

				// The above need some inits to prepare for the AVX version of the Fermat-mod carry macro:
				SSE2_fermat_carry_init_loacc(base_negacyclic_root);

			  #endif

				// AVX-custom 4-way carry macro - each contains 4 of the RADIX stride-n/RADIX-separated carries
				// (processed independently in parallel), and steps through sequential-data indicies j,j+2,j+4,j+6:
			  #if HIACC
				// The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl)
				// which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
				// precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
				// uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls:

				tmp = base_negacyclic_root+  0;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p00r,tmp,0xf00,cy_r00,cy_i00,odd_radix,half_arr,sign_mask,icycle00,icycle01,icycle02,icycle03, jcycle00,kcycle00,lcycle00);	// *cycle0 index increments by +4 (mod odd_radix) between macro calls
				tmp = base_negacyclic_root+  8;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p04r,tmp,0xe40,cy_r04,cy_i04,odd_radix,half_arr,sign_mask,icycle04,icycle05,icycle06,icycle07, jcycle04,kcycle04,lcycle04);
				tmp = base_negacyclic_root+ 16;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p08r,tmp,0xd80,cy_r08,cy_i08,odd_radix,half_arr,sign_mask,icycle08,icycle09,icycle10,icycle11, jcycle08,kcycle08,lcycle08);
				tmp = base_negacyclic_root+ 24;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p0cr,tmp,0xcc0,cy_r12,cy_i12,odd_radix,half_arr,sign_mask,icycle12,icycle13,icycle14,icycle00, jcycle12,kcycle12,lcycle12);
				tmp = base_negacyclic_root+ 32;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p11r,tmp,0xc00,cy_r16,cy_i16,odd_radix,half_arr,sign_mask,icycle01,icycle02,icycle03,icycle04, jcycle01,kcycle01,lcycle01);
				tmp = base_negacyclic_root+ 40;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p15r,tmp,0xb40,cy_r20,cy_i20,odd_radix,half_arr,sign_mask,icycle05,icycle06,icycle07,icycle08, jcycle05,kcycle05,lcycle05);
				tmp = base_negacyclic_root+ 48;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p19r,tmp,0xa80,cy_r24,cy_i24,odd_radix,half_arr,sign_mask,icycle09,icycle10,icycle11,icycle12, jcycle09,kcycle09,lcycle09);
				tmp = base_negacyclic_root+ 56;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p1dr,tmp,0x9c0,cy_r28,cy_i28,odd_radix,half_arr,sign_mask,icycle13,icycle14,icycle00,icycle01, jcycle13,kcycle13,lcycle13);
				tmp = base_negacyclic_root+ 64;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p22r,tmp,0x900,cy_r32,cy_i32,odd_radix,half_arr,sign_mask,icycle02,icycle03,icycle04,icycle05, jcycle02,kcycle02,lcycle02);
				tmp = base_negacyclic_root+ 72;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p26r,tmp,0x840,cy_r36,cy_i36,odd_radix,half_arr,sign_mask,icycle06,icycle07,icycle08,icycle09, jcycle06,kcycle06,lcycle06);
				tmp = base_negacyclic_root+ 80;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p2ar,tmp,0x780,cy_r40,cy_i40,odd_radix,half_arr,sign_mask,icycle10,icycle11,icycle12,icycle13, jcycle10,kcycle10,lcycle10);
				tmp = base_negacyclic_root+ 88;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p2er,tmp,0x6c0,cy_r44,cy_i44,odd_radix,half_arr,sign_mask,icycle14,icycle00,icycle01,icycle02, jcycle14,kcycle14,lcycle14);
				tmp = base_negacyclic_root+ 96;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p33r,tmp,0x600,cy_r48,cy_i48,odd_radix,half_arr,sign_mask,icycle03,icycle04,icycle05,icycle06, jcycle03,kcycle03,lcycle03);
				tmp = base_negacyclic_root+104;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p37r,tmp,0x540,cy_r52,cy_i52,odd_radix,half_arr,sign_mask,icycle07,icycle08,icycle09,icycle10, jcycle07,kcycle07,lcycle07);
				tmp = base_negacyclic_root+112;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p3br,tmp,0x480,cy_r56,cy_i56,odd_radix,half_arr,sign_mask,icycle11,icycle12,icycle13,icycle14, jcycle11,kcycle11,lcycle11);

			  #else	// HIACC = false:

				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p00r,base_negacyclic_root,cy_r00,cy_i00,odd_radix,half_arr,sign_mask,                     icycle00,icycle01,icycle02,icycle03, jcycle00,kcycle00,lcycle00);	// *cycle index increments by +4 (mod odd_radix) between macro calls
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p04r,base_negacyclic_root,cy_r04,cy_i04,odd_radix,half_arr,sign_mask,                     icycle04,icycle05,icycle06,icycle07, jcycle04,kcycle04,lcycle04);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p08r,base_negacyclic_root,cy_r08,cy_i08,odd_radix,half_arr,sign_mask,                     icycle08,icycle09,icycle10,icycle11, jcycle08,kcycle08,lcycle08);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p0cr,base_negacyclic_root,cy_r12,cy_i12,odd_radix,half_arr,sign_mask,                     icycle12,icycle13,icycle14,icycle00, jcycle12,kcycle12,lcycle12);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p11r,base_negacyclic_root,cy_r16,cy_i16,odd_radix,half_arr,sign_mask,                     icycle01,icycle02,icycle03,icycle04, jcycle01,kcycle01,lcycle01);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p15r,base_negacyclic_root,cy_r20,cy_i20,odd_radix,half_arr,sign_mask,                     icycle05,icycle06,icycle07,icycle08, jcycle05,kcycle05,lcycle05);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p19r,base_negacyclic_root,cy_r24,cy_i24,odd_radix,half_arr,sign_mask,                     icycle09,icycle10,icycle11,icycle12, jcycle09,kcycle09,lcycle09);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p1dr,base_negacyclic_root,cy_r28,cy_i28,odd_radix,half_arr,sign_mask,                     icycle13,icycle14,icycle00,icycle01, jcycle13,kcycle13,lcycle13);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p22r,base_negacyclic_root,cy_r32,cy_i32,odd_radix,half_arr,sign_mask,                     icycle02,icycle03,icycle04,icycle05, jcycle02,kcycle02,lcycle02);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p26r,base_negacyclic_root,cy_r36,cy_i36,odd_radix,half_arr,sign_mask,                     icycle06,icycle07,icycle08,icycle09, jcycle06,kcycle06,lcycle06);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p2ar,base_negacyclic_root,cy_r40,cy_i40,odd_radix,half_arr,sign_mask,                     icycle10,icycle11,icycle12,icycle13, jcycle10,kcycle10,lcycle10);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p2er,base_negacyclic_root,cy_r44,cy_i44,odd_radix,half_arr,sign_mask,                     icycle14,icycle00,icycle01,icycle02, jcycle14,kcycle14,lcycle14);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p33r,base_negacyclic_root,cy_r48,cy_i48,odd_radix,half_arr,sign_mask,                     icycle03,icycle04,icycle05,icycle06, jcycle03,kcycle03,lcycle03);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p37r,base_negacyclic_root,cy_r52,cy_i52,odd_radix,half_arr,sign_mask,                     icycle07,icycle08,icycle09,icycle10, jcycle07,kcycle07,lcycle07);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p3br,base_negacyclic_root,cy_r56,cy_i56,odd_radix,half_arr,sign_mask,                     icycle11,icycle12,icycle13,icycle14, jcycle11,kcycle11,lcycle11);

			  #endif	// HIACC?

			#elif defined(USE_SSE2)

				// For a description of the data movement for Fermat-mod carries in SSE2 mode, see radix16_ditN_cy_dif1.c.

				/* Get the needed Nth root of -1: */
				add1 = (double *)&rn0[0];
				add2 = (double *)&rn1[0];

				idx_offset = j;
				idx_incr = NDIVR;

			  #if defined(COMPILER_TYPE_MSVC)
			  /* The cy_[r|i]_idx[A|B] names here are not meaningful, each simply stores one [re,im] carry pair,
			  e.g. cy_r01 stores the carries our of [a0.re,a0.im], cy_r23 stores the carries our of [a1.re,a1.im], etc.
			  Here is the actual mapping between these SSE2-mode 2-vector carry pairs and the scalar carries:
														2-vector						 Scalar
														--------	 				   ------------- */
				SSE2_fermat_carry_norm_errcheck(s1p00r,cy_r00,idx_offset,idx_incr,odd_radix,icycle00,jcycle00);	/* cy_r00,cy_i00 */
				SSE2_fermat_carry_norm_errcheck(s1p01r,cy_r02,idx_offset,idx_incr,odd_radix,icycle01,jcycle01);	/* cy_r01,cy_i01 */
				SSE2_fermat_carry_norm_errcheck(s1p02r,cy_r04,idx_offset,idx_incr,odd_radix,icycle02,jcycle02);	/* cy_r02,cy_i02 */
				SSE2_fermat_carry_norm_errcheck(s1p03r,cy_r06,idx_offset,idx_incr,odd_radix,icycle03,jcycle03);	/* cy_r03,cy_i03 */
				SSE2_fermat_carry_norm_errcheck(s1p04r,cy_r08,idx_offset,idx_incr,odd_radix,icycle04,jcycle04);	/* cy_r04,cy_i04 */
				SSE2_fermat_carry_norm_errcheck(s1p05r,cy_r10,idx_offset,idx_incr,odd_radix,icycle05,jcycle05);	/* cy_r05,cy_i05 */
				SSE2_fermat_carry_norm_errcheck(s1p06r,cy_r12,idx_offset,idx_incr,odd_radix,icycle06,jcycle06);	/* cy_r06,cy_i06 */
				SSE2_fermat_carry_norm_errcheck(s1p07r,cy_r14,idx_offset,idx_incr,odd_radix,icycle07,jcycle07);	/* cy_r07,cy_i07 */
				SSE2_fermat_carry_norm_errcheck(s1p08r,cy_r16,idx_offset,idx_incr,odd_radix,icycle08,jcycle08);	/* cy_r08,cy_i08 */
				SSE2_fermat_carry_norm_errcheck(s1p09r,cy_r18,idx_offset,idx_incr,odd_radix,icycle09,jcycle09);	/* cy_r09,cy_i09 */
				SSE2_fermat_carry_norm_errcheck(s1p0ar,cy_r20,idx_offset,idx_incr,odd_radix,icycle10,jcycle10);	/* cy_r10,cy_i10 */
				SSE2_fermat_carry_norm_errcheck(s1p0br,cy_r22,idx_offset,idx_incr,odd_radix,icycle11,jcycle11);	/* cy_r11,cy_i11 */
				SSE2_fermat_carry_norm_errcheck(s1p0cr,cy_r24,idx_offset,idx_incr,odd_radix,icycle12,jcycle12);	/* cy_r12,cy_i12 */
				SSE2_fermat_carry_norm_errcheck(s1p0dr,cy_r26,idx_offset,idx_incr,odd_radix,icycle13,jcycle13);	/* cy_r13,cy_i13 */
				SSE2_fermat_carry_norm_errcheck(s1p0er,cy_r28,idx_offset,idx_incr,odd_radix,icycle14,jcycle14);	/* cy_r14,cy_i14 */
				SSE2_fermat_carry_norm_errcheck(s1p10r,cy_r30,idx_offset,idx_incr,odd_radix,icycle00,jcycle00);	/* cy_r15,cy_i15 */
				SSE2_fermat_carry_norm_errcheck(s1p11r,cy_r32,idx_offset,idx_incr,odd_radix,icycle01,jcycle01);	/* cy_r16,cy_i16 */
				SSE2_fermat_carry_norm_errcheck(s1p12r,cy_r34,idx_offset,idx_incr,odd_radix,icycle02,jcycle02);	/* cy_r17,cy_i17 */
				SSE2_fermat_carry_norm_errcheck(s1p13r,cy_r36,idx_offset,idx_incr,odd_radix,icycle03,jcycle03);	/* cy_r18,cy_i18 */
				SSE2_fermat_carry_norm_errcheck(s1p14r,cy_r38,idx_offset,idx_incr,odd_radix,icycle04,jcycle04);	/* cy_r19,cy_i19 */
				SSE2_fermat_carry_norm_errcheck(s1p15r,cy_r40,idx_offset,idx_incr,odd_radix,icycle05,jcycle05);	/* cy_r20,cy_i20 */
				SSE2_fermat_carry_norm_errcheck(s1p16r,cy_r42,idx_offset,idx_incr,odd_radix,icycle06,jcycle06);	/* cy_r21,cy_i21 */
				SSE2_fermat_carry_norm_errcheck(s1p17r,cy_r44,idx_offset,idx_incr,odd_radix,icycle07,jcycle07);	/* cy_r22,cy_i22 */
				SSE2_fermat_carry_norm_errcheck(s1p18r,cy_r46,idx_offset,idx_incr,odd_radix,icycle08,jcycle08);	/* cy_r23,cy_i23 */
				SSE2_fermat_carry_norm_errcheck(s1p19r,cy_r48,idx_offset,idx_incr,odd_radix,icycle09,jcycle09);	/* cy_r24,cy_i24 */
				SSE2_fermat_carry_norm_errcheck(s1p1ar,cy_r50,idx_offset,idx_incr,odd_radix,icycle10,jcycle10);	/* cy_r25,cy_i25 */
				SSE2_fermat_carry_norm_errcheck(s1p1br,cy_r52,idx_offset,idx_incr,odd_radix,icycle11,jcycle11);	/* cy_r26,cy_i26 */
				SSE2_fermat_carry_norm_errcheck(s1p1cr,cy_r54,idx_offset,idx_incr,odd_radix,icycle12,jcycle12);	/* cy_r27,cy_i27 */
				SSE2_fermat_carry_norm_errcheck(s1p1dr,cy_r56,idx_offset,idx_incr,odd_radix,icycle13,jcycle13);	/* cy_r28,cy_i28 */
				SSE2_fermat_carry_norm_errcheck(s1p1er,cy_r58,idx_offset,idx_incr,odd_radix,icycle14,jcycle14);	/* cy_r29,cy_i29 */
				SSE2_fermat_carry_norm_errcheck(s1p20r,cy_i00,idx_offset,idx_incr,odd_radix,icycle00,jcycle00);	/* cy_r30,cy_i30 */
				SSE2_fermat_carry_norm_errcheck(s1p21r,cy_i02,idx_offset,idx_incr,odd_radix,icycle01,jcycle01);	/* cy_r31,cy_i31 */
				SSE2_fermat_carry_norm_errcheck(s1p22r,cy_i04,idx_offset,idx_incr,odd_radix,icycle02,jcycle02);	/* cy_r32,cy_i32 */
				SSE2_fermat_carry_norm_errcheck(s1p23r,cy_i06,idx_offset,idx_incr,odd_radix,icycle03,jcycle03);	/* cy_r33,cy_i33 */
				SSE2_fermat_carry_norm_errcheck(s1p24r,cy_i08,idx_offset,idx_incr,odd_radix,icycle04,jcycle04);	/* cy_r34,cy_i34 */
				SSE2_fermat_carry_norm_errcheck(s1p25r,cy_i10,idx_offset,idx_incr,odd_radix,icycle05,jcycle05);	/* cy_r35,cy_i35 */
				SSE2_fermat_carry_norm_errcheck(s1p26r,cy_i12,idx_offset,idx_incr,odd_radix,icycle06,jcycle06);	/* cy_r36,cy_i36 */
				SSE2_fermat_carry_norm_errcheck(s1p27r,cy_i14,idx_offset,idx_incr,odd_radix,icycle07,jcycle07);	/* cy_r37,cy_i37 */
				SSE2_fermat_carry_norm_errcheck(s1p28r,cy_i16,idx_offset,idx_incr,odd_radix,icycle08,jcycle08);	/* cy_r38,cy_i38 */
				SSE2_fermat_carry_norm_errcheck(s1p29r,cy_i18,idx_offset,idx_incr,odd_radix,icycle09,jcycle09);	/* cy_r39,cy_i39 */
				SSE2_fermat_carry_norm_errcheck(s1p2ar,cy_i20,idx_offset,idx_incr,odd_radix,icycle10,jcycle10);	/* cy_r40,cy_i40 */
				SSE2_fermat_carry_norm_errcheck(s1p2br,cy_i22,idx_offset,idx_incr,odd_radix,icycle11,jcycle11);	/* cy_r41,cy_i41 */
				SSE2_fermat_carry_norm_errcheck(s1p2cr,cy_i24,idx_offset,idx_incr,odd_radix,icycle12,jcycle12);	/* cy_r42,cy_i42 */
				SSE2_fermat_carry_norm_errcheck(s1p2dr,cy_i26,idx_offset,idx_incr,odd_radix,icycle13,jcycle13);	/* cy_r43,cy_i43 */
				SSE2_fermat_carry_norm_errcheck(s1p2er,cy_i28,idx_offset,idx_incr,odd_radix,icycle14,jcycle14);	/* cy_r44,cy_i44 */
				SSE2_fermat_carry_norm_errcheck(s1p30r,cy_i30,idx_offset,idx_incr,odd_radix,icycle00,jcycle00);	/* cy_r45,cy_i45 */
				SSE2_fermat_carry_norm_errcheck(s1p31r,cy_i32,idx_offset,idx_incr,odd_radix,icycle01,jcycle01);	/* cy_r46,cy_i46 */
				SSE2_fermat_carry_norm_errcheck(s1p32r,cy_i34,idx_offset,idx_incr,odd_radix,icycle02,jcycle02);	/* cy_r47,cy_i47 */
				SSE2_fermat_carry_norm_errcheck(s1p33r,cy_i36,idx_offset,idx_incr,odd_radix,icycle03,jcycle03);	/* cy_r48,cy_i48 */
				SSE2_fermat_carry_norm_errcheck(s1p34r,cy_i38,idx_offset,idx_incr,odd_radix,icycle04,jcycle04);	/* cy_r49,cy_i49 */
				SSE2_fermat_carry_norm_errcheck(s1p35r,cy_i40,idx_offset,idx_incr,odd_radix,icycle05,jcycle05);	/* cy_r50,cy_i50 */
				SSE2_fermat_carry_norm_errcheck(s1p36r,cy_i42,idx_offset,idx_incr,odd_radix,icycle06,jcycle06);	/* cy_r51,cy_i51 */
				SSE2_fermat_carry_norm_errcheck(s1p37r,cy_i44,idx_offset,idx_incr,odd_radix,icycle07,jcycle07);	/* cy_r52,cy_i52 */
				SSE2_fermat_carry_norm_errcheck(s1p38r,cy_i46,idx_offset,idx_incr,odd_radix,icycle08,jcycle08);	/* cy_r53,cy_i53 */
				SSE2_fermat_carry_norm_errcheck(s1p39r,cy_i48,idx_offset,idx_incr,odd_radix,icycle09,jcycle09);	/* cy_r54,cy_i54 */
				SSE2_fermat_carry_norm_errcheck(s1p3ar,cy_i50,idx_offset,idx_incr,odd_radix,icycle10,jcycle10);	/* cy_r55,cy_i55 */
				SSE2_fermat_carry_norm_errcheck(s1p3br,cy_i52,idx_offset,idx_incr,odd_radix,icycle11,jcycle11);	/* cy_r56,cy_i56 */
				SSE2_fermat_carry_norm_errcheck(s1p3cr,cy_i54,idx_offset,idx_incr,odd_radix,icycle12,jcycle12);	/* cy_r57,cy_i57 */
				SSE2_fermat_carry_norm_errcheck(s1p3dr,cy_i56,idx_offset,idx_incr,odd_radix,icycle13,jcycle13);	/* cy_r58,cy_i58 */
				SSE2_fermat_carry_norm_errcheck(s1p3er,cy_i58,idx_offset,idx_incr,odd_radix,icycle14,jcycle14);	/* cy_r59,cy_i59 */

			  #elif (OS_BITS == 32)

				SSE2_fermat_carry_norm_errcheck(s1p00r,cy_r00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle00,jcycle00);
				SSE2_fermat_carry_norm_errcheck(s1p01r,cy_r02,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle01,jcycle01);
				SSE2_fermat_carry_norm_errcheck(s1p02r,cy_r04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle02,jcycle02);
				SSE2_fermat_carry_norm_errcheck(s1p03r,cy_r06,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle03,jcycle03);
				SSE2_fermat_carry_norm_errcheck(s1p04r,cy_r08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle04,jcycle04);
				SSE2_fermat_carry_norm_errcheck(s1p05r,cy_r10,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle05,jcycle05);
				SSE2_fermat_carry_norm_errcheck(s1p06r,cy_r12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle06,jcycle06);
				SSE2_fermat_carry_norm_errcheck(s1p07r,cy_r14,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle07,jcycle07);
				SSE2_fermat_carry_norm_errcheck(s1p08r,cy_r16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle08,jcycle08);
				SSE2_fermat_carry_norm_errcheck(s1p09r,cy_r18,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle09,jcycle09);
				SSE2_fermat_carry_norm_errcheck(s1p0ar,cy_r20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle10,jcycle10);
				SSE2_fermat_carry_norm_errcheck(s1p0br,cy_r22,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle11,jcycle11);
				SSE2_fermat_carry_norm_errcheck(s1p0cr,cy_r24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle12,jcycle12);
				SSE2_fermat_carry_norm_errcheck(s1p0dr,cy_r26,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle13,jcycle13);
				SSE2_fermat_carry_norm_errcheck(s1p0er,cy_r28,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle14,jcycle14);
				SSE2_fermat_carry_norm_errcheck(s1p10r,cy_r30,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle00,jcycle00);
				SSE2_fermat_carry_norm_errcheck(s1p11r,cy_r32,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle01,jcycle01);
				SSE2_fermat_carry_norm_errcheck(s1p12r,cy_r34,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle02,jcycle02);
				SSE2_fermat_carry_norm_errcheck(s1p13r,cy_r36,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle03,jcycle03);
				SSE2_fermat_carry_norm_errcheck(s1p14r,cy_r38,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle04,jcycle04);
				SSE2_fermat_carry_norm_errcheck(s1p15r,cy_r40,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle05,jcycle05);
				SSE2_fermat_carry_norm_errcheck(s1p16r,cy_r42,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle06,jcycle06);
				SSE2_fermat_carry_norm_errcheck(s1p17r,cy_r44,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle07,jcycle07);
				SSE2_fermat_carry_norm_errcheck(s1p18r,cy_r46,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle08,jcycle08);
				SSE2_fermat_carry_norm_errcheck(s1p19r,cy_r48,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle09,jcycle09);
				SSE2_fermat_carry_norm_errcheck(s1p1ar,cy_r50,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle10,jcycle10);
				SSE2_fermat_carry_norm_errcheck(s1p1br,cy_r52,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle11,jcycle11);
				SSE2_fermat_carry_norm_errcheck(s1p1cr,cy_r54,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle12,jcycle12);
				SSE2_fermat_carry_norm_errcheck(s1p1dr,cy_r56,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle13,jcycle13);
				SSE2_fermat_carry_norm_errcheck(s1p1er,cy_r58,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle14,jcycle14);
				SSE2_fermat_carry_norm_errcheck(s1p20r,cy_i00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle00,jcycle00);
				SSE2_fermat_carry_norm_errcheck(s1p21r,cy_i02,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle01,jcycle01);
				SSE2_fermat_carry_norm_errcheck(s1p22r,cy_i04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle02,jcycle02);
				SSE2_fermat_carry_norm_errcheck(s1p23r,cy_i06,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle03,jcycle03);
				SSE2_fermat_carry_norm_errcheck(s1p24r,cy_i08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle04,jcycle04);
				SSE2_fermat_carry_norm_errcheck(s1p25r,cy_i10,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle05,jcycle05);
				SSE2_fermat_carry_norm_errcheck(s1p26r,cy_i12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle06,jcycle06);
				SSE2_fermat_carry_norm_errcheck(s1p27r,cy_i14,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle07,jcycle07);
				SSE2_fermat_carry_norm_errcheck(s1p28r,cy_i16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle08,jcycle08);
				SSE2_fermat_carry_norm_errcheck(s1p29r,cy_i18,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle09,jcycle09);
				SSE2_fermat_carry_norm_errcheck(s1p2ar,cy_i20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle10,jcycle10);
				SSE2_fermat_carry_norm_errcheck(s1p2br,cy_i22,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle11,jcycle11);
				SSE2_fermat_carry_norm_errcheck(s1p2cr,cy_i24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle12,jcycle12);
				SSE2_fermat_carry_norm_errcheck(s1p2dr,cy_i26,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle13,jcycle13);
				SSE2_fermat_carry_norm_errcheck(s1p2er,cy_i28,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle14,jcycle14);
				SSE2_fermat_carry_norm_errcheck(s1p30r,cy_i30,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle00,jcycle00);
				SSE2_fermat_carry_norm_errcheck(s1p31r,cy_i32,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle01,jcycle01);
				SSE2_fermat_carry_norm_errcheck(s1p32r,cy_i34,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle02,jcycle02);
				SSE2_fermat_carry_norm_errcheck(s1p33r,cy_i36,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle03,jcycle03);
				SSE2_fermat_carry_norm_errcheck(s1p34r,cy_i38,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle04,jcycle04);
				SSE2_fermat_carry_norm_errcheck(s1p35r,cy_i40,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle05,jcycle05);
				SSE2_fermat_carry_norm_errcheck(s1p36r,cy_i42,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle06,jcycle06);
				SSE2_fermat_carry_norm_errcheck(s1p37r,cy_i44,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle07,jcycle07);
				SSE2_fermat_carry_norm_errcheck(s1p38r,cy_i46,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle08,jcycle08);
				SSE2_fermat_carry_norm_errcheck(s1p39r,cy_i48,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle09,jcycle09);
				SSE2_fermat_carry_norm_errcheck(s1p3ar,cy_i50,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle10,jcycle10);
				SSE2_fermat_carry_norm_errcheck(s1p3br,cy_i52,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle11,jcycle11);
				SSE2_fermat_carry_norm_errcheck(s1p3cr,cy_i54,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle12,jcycle12);
				SSE2_fermat_carry_norm_errcheck(s1p3dr,cy_i56,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle13,jcycle13);
				SSE2_fermat_carry_norm_errcheck(s1p3er,cy_i58,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle14,jcycle14);

			  #else	// 64-bit SSE2

				SSE2_fermat_carry_norm_errcheck_X2(s1p00r,cy_r00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle00,jcycle00,icycle01,jcycle01);
				SSE2_fermat_carry_norm_errcheck_X2(s1p02r,cy_r04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle02,jcycle02,icycle03,jcycle03);
				SSE2_fermat_carry_norm_errcheck_X2(s1p04r,cy_r08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle04,jcycle04,icycle05,jcycle05);
				SSE2_fermat_carry_norm_errcheck_X2(s1p06r,cy_r12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle06,jcycle06,icycle07,jcycle07);
				SSE2_fermat_carry_norm_errcheck_X2(s1p08r,cy_r16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle08,jcycle08,icycle09,jcycle09);
				SSE2_fermat_carry_norm_errcheck_X2(s1p0ar,cy_r20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle10,jcycle10,icycle11,jcycle11);
				SSE2_fermat_carry_norm_errcheck_X2(s1p0cr,cy_r24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle12,jcycle12,icycle13,jcycle13);
				SSE2_fermat_carry_norm_errcheck_X2(s1p0er,cy_r28,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle14,jcycle14,icycle00,jcycle00);
				SSE2_fermat_carry_norm_errcheck_X2(s1p11r,cy_r32,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle01,jcycle01,icycle02,jcycle02);
				SSE2_fermat_carry_norm_errcheck_X2(s1p13r,cy_r36,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle03,jcycle03,icycle04,jcycle04);
				SSE2_fermat_carry_norm_errcheck_X2(s1p15r,cy_r40,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle05,jcycle05,icycle06,jcycle06);
				SSE2_fermat_carry_norm_errcheck_X2(s1p17r,cy_r44,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle07,jcycle07,icycle08,jcycle08);
				SSE2_fermat_carry_norm_errcheck_X2(s1p19r,cy_r48,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle09,jcycle09,icycle10,jcycle10);
				SSE2_fermat_carry_norm_errcheck_X2(s1p1br,cy_r52,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle11,jcycle11,icycle12,jcycle12);
				SSE2_fermat_carry_norm_errcheck_X2(s1p1dr,cy_r56,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle13,jcycle13,icycle14,jcycle14);
				SSE2_fermat_carry_norm_errcheck_X2(s1p20r,cy_i00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle00,jcycle00,icycle01,jcycle01);
				SSE2_fermat_carry_norm_errcheck_X2(s1p22r,cy_i04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle02,jcycle02,icycle03,jcycle03);
				SSE2_fermat_carry_norm_errcheck_X2(s1p24r,cy_i08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle04,jcycle04,icycle05,jcycle05);
				SSE2_fermat_carry_norm_errcheck_X2(s1p26r,cy_i12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle06,jcycle06,icycle07,jcycle07);
				SSE2_fermat_carry_norm_errcheck_X2(s1p28r,cy_i16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle08,jcycle08,icycle09,jcycle09);
				SSE2_fermat_carry_norm_errcheck_X2(s1p2ar,cy_i20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle10,jcycle10,icycle11,jcycle11);
				SSE2_fermat_carry_norm_errcheck_X2(s1p2cr,cy_i24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle12,jcycle12,icycle13,jcycle13);
				SSE2_fermat_carry_norm_errcheck_X2(s1p2er,cy_i28,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle14,jcycle14,icycle00,jcycle00);
				SSE2_fermat_carry_norm_errcheck_X2(s1p31r,cy_i32,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle01,jcycle01,icycle02,jcycle02);
				SSE2_fermat_carry_norm_errcheck_X2(s1p33r,cy_i36,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle03,jcycle03,icycle04,jcycle04);
				SSE2_fermat_carry_norm_errcheck_X2(s1p35r,cy_i40,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle05,jcycle05,icycle06,jcycle06);
				SSE2_fermat_carry_norm_errcheck_X2(s1p37r,cy_i44,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle07,jcycle07,icycle08,jcycle08);
				SSE2_fermat_carry_norm_errcheck_X2(s1p39r,cy_i48,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle09,jcycle09,icycle10,jcycle10);
				SSE2_fermat_carry_norm_errcheck_X2(s1p3br,cy_i52,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle11,jcycle11,icycle12,jcycle12);
				SSE2_fermat_carry_norm_errcheck_X2(s1p3dr,cy_i56,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle13,jcycle13,icycle14,jcycle14);
			  #endif

			#else	// Scalar-double mode:

				ntmp = 0;
				fermat_carry_norm_errcheckB(a1p00r,a1p00i,cy_r00,cy_i00,icycle00,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p01r,a1p01i,cy_r01,cy_i01,icycle01,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p02r,a1p02i,cy_r02,cy_i02,icycle02,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p03r,a1p03i,cy_r03,cy_i03,icycle03,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p04r,a1p04i,cy_r04,cy_i04,icycle04,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p05r,a1p05i,cy_r05,cy_i05,icycle05,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p06r,a1p06i,cy_r06,cy_i06,icycle06,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p07r,a1p07i,cy_r07,cy_i07,icycle07,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p08r,a1p08i,cy_r08,cy_i08,icycle08,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p09r,a1p09i,cy_r09,cy_i09,icycle09,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p0ar,a1p0ai,cy_r10,cy_i10,icycle10,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p0br,a1p0bi,cy_r11,cy_i11,icycle11,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p0cr,a1p0ci,cy_r12,cy_i12,icycle12,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p0dr,a1p0di,cy_r13,cy_i13,icycle13,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p0er,a1p0ei,cy_r14,cy_i14,icycle14,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p10r,a1p10i,cy_r15,cy_i15,icycle00,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p11r,a1p11i,cy_r16,cy_i16,icycle01,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p12r,a1p12i,cy_r17,cy_i17,icycle02,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p13r,a1p13i,cy_r18,cy_i18,icycle03,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p14r,a1p14i,cy_r19,cy_i19,icycle04,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p15r,a1p15i,cy_r20,cy_i20,icycle05,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p16r,a1p16i,cy_r21,cy_i21,icycle06,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p17r,a1p17i,cy_r22,cy_i22,icycle07,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p18r,a1p18i,cy_r23,cy_i23,icycle08,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p19r,a1p19i,cy_r24,cy_i24,icycle09,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p1ar,a1p1ai,cy_r25,cy_i25,icycle10,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p1br,a1p1bi,cy_r26,cy_i26,icycle11,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p1cr,a1p1ci,cy_r27,cy_i27,icycle12,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p1dr,a1p1di,cy_r28,cy_i28,icycle13,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p1er,a1p1ei,cy_r29,cy_i29,icycle14,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p20r,a1p20i,cy_r30,cy_i30,icycle00,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p21r,a1p21i,cy_r31,cy_i31,icycle01,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p22r,a1p22i,cy_r32,cy_i32,icycle02,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p23r,a1p23i,cy_r33,cy_i33,icycle03,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p24r,a1p24i,cy_r34,cy_i34,icycle04,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p25r,a1p25i,cy_r35,cy_i35,icycle05,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p26r,a1p26i,cy_r36,cy_i36,icycle06,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p27r,a1p27i,cy_r37,cy_i37,icycle07,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p28r,a1p28i,cy_r38,cy_i38,icycle08,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p29r,a1p29i,cy_r39,cy_i39,icycle09,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p2ar,a1p2ai,cy_r40,cy_i40,icycle10,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p2br,a1p2bi,cy_r41,cy_i41,icycle11,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p2cr,a1p2ci,cy_r42,cy_i42,icycle12,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p2dr,a1p2di,cy_r43,cy_i43,icycle13,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p2er,a1p2ei,cy_r44,cy_i44,icycle14,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p30r,a1p30i,cy_r45,cy_i45,icycle00,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p31r,a1p31i,cy_r46,cy_i46,icycle01,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p32r,a1p32i,cy_r47,cy_i47,icycle02,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p33r,a1p33i,cy_r48,cy_i48,icycle03,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p34r,a1p34i,cy_r49,cy_i49,icycle04,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p35r,a1p35i,cy_r50,cy_i50,icycle05,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p36r,a1p36i,cy_r51,cy_i51,icycle06,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p37r,a1p37i,cy_r52,cy_i52,icycle07,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p38r,a1p38i,cy_r53,cy_i53,icycle08,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p39r,a1p39i,cy_r54,cy_i54,icycle09,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p3ar,a1p3ai,cy_r55,cy_i55,icycle10,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p3br,a1p3bi,cy_r56,cy_i56,icycle11,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p3cr,a1p3ci,cy_r57,cy_i57,icycle12,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p3dr,a1p3di,cy_r58,cy_i58,icycle13,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p3er,a1p3ei,cy_r59,cy_i59,icycle14,ntmp,NRTM1,NRT_BITS);

				icycle00 += wts_idx_incr;	/* Inside the loop use this, as it is faster than general-mod '% nwt' */
				icycle01 += wts_idx_incr;
				icycle02 += wts_idx_incr;
				icycle03 += wts_idx_incr;
				icycle04 += wts_idx_incr;
				icycle05 += wts_idx_incr;
				icycle06 += wts_idx_incr;
				icycle07 += wts_idx_incr;
				icycle08 += wts_idx_incr;
				icycle09 += wts_idx_incr;
				icycle10 += wts_idx_incr;
				icycle11 += wts_idx_incr;
				icycle12 += wts_idx_incr;
				icycle13 += wts_idx_incr;
				icycle14 += wts_idx_incr;
				icycle00 += ( (-(int)((uint32)icycle00 >> 31)) & nwt);
				icycle01 += ( (-(int)((uint32)icycle01 >> 31)) & nwt);
				icycle02 += ( (-(int)((uint32)icycle02 >> 31)) & nwt);
				icycle03 += ( (-(int)((uint32)icycle03 >> 31)) & nwt);
				icycle04 += ( (-(int)((uint32)icycle04 >> 31)) & nwt);
				icycle05 += ( (-(int)((uint32)icycle05 >> 31)) & nwt);
				icycle06 += ( (-(int)((uint32)icycle06 >> 31)) & nwt);
				icycle07 += ( (-(int)((uint32)icycle07 >> 31)) & nwt);
				icycle08 += ( (-(int)((uint32)icycle08 >> 31)) & nwt);
				icycle09 += ( (-(int)((uint32)icycle09 >> 31)) & nwt);
				icycle10 += ( (-(int)((uint32)icycle10 >> 31)) & nwt);
				icycle11 += ( (-(int)((uint32)icycle11 >> 31)) & nwt);
				icycle12 += ( (-(int)((uint32)icycle12 >> 31)) & nwt);
				icycle13 += ( (-(int)((uint32)icycle13 >> 31)) & nwt);
				icycle14 += ( (-(int)((uint32)icycle14 >> 31)) & nwt);

			#endif	/* #ifdef USE_SSE2 */

			// Here we nest AVX inside SSE2 since i/jcycle updates are for both, k/l for AVX-only:
			#ifdef USE_SSE2

				icycle00 += wts_idx_inc2;		icycle00 += ( (-(icycle00 < 0)) & nwt16);
				icycle01 += wts_idx_inc2;		icycle01 += ( (-(icycle01 < 0)) & nwt16);
				icycle02 += wts_idx_inc2;		icycle02 += ( (-(icycle02 < 0)) & nwt16);
				icycle03 += wts_idx_inc2;		icycle03 += ( (-(icycle03 < 0)) & nwt16);
				icycle04 += wts_idx_inc2;		icycle04 += ( (-(icycle04 < 0)) & nwt16);
				icycle05 += wts_idx_inc2;		icycle05 += ( (-(icycle05 < 0)) & nwt16);
				icycle06 += wts_idx_inc2;		icycle06 += ( (-(icycle06 < 0)) & nwt16);
				icycle07 += wts_idx_inc2;		icycle07 += ( (-(icycle07 < 0)) & nwt16);
				icycle08 += wts_idx_inc2;		icycle08 += ( (-(icycle08 < 0)) & nwt16);
				icycle09 += wts_idx_inc2;		icycle09 += ( (-(icycle09 < 0)) & nwt16);
				icycle10 += wts_idx_inc2;		icycle10 += ( (-(icycle10 < 0)) & nwt16);
				icycle11 += wts_idx_inc2;		icycle11 += ( (-(icycle11 < 0)) & nwt16);
				icycle12 += wts_idx_inc2;		icycle12 += ( (-(icycle12 < 0)) & nwt16);
				icycle13 += wts_idx_inc2;		icycle13 += ( (-(icycle13 < 0)) & nwt16);
				icycle14 += wts_idx_inc2;		icycle14 += ( (-(icycle14 < 0)) & nwt16);

				jcycle00 += wts_idx_inc2;		jcycle00 += ( (-(jcycle00 < 0)) & nwt16);
				jcycle01 += wts_idx_inc2;		jcycle01 += ( (-(jcycle01 < 0)) & nwt16);
				jcycle02 += wts_idx_inc2;		jcycle02 += ( (-(jcycle02 < 0)) & nwt16);
				jcycle03 += wts_idx_inc2;		jcycle03 += ( (-(jcycle03 < 0)) & nwt16);
				jcycle04 += wts_idx_inc2;		jcycle04 += ( (-(jcycle04 < 0)) & nwt16);
				jcycle05 += wts_idx_inc2;		jcycle05 += ( (-(jcycle05 < 0)) & nwt16);
				jcycle06 += wts_idx_inc2;		jcycle06 += ( (-(jcycle06 < 0)) & nwt16);
				jcycle07 += wts_idx_inc2;		jcycle07 += ( (-(jcycle07 < 0)) & nwt16);
				jcycle08 += wts_idx_inc2;		jcycle08 += ( (-(jcycle08 < 0)) & nwt16);
				jcycle09 += wts_idx_inc2;		jcycle09 += ( (-(jcycle09 < 0)) & nwt16);
				jcycle10 += wts_idx_inc2;		jcycle10 += ( (-(jcycle10 < 0)) & nwt16);
				jcycle11 += wts_idx_inc2;		jcycle11 += ( (-(jcycle11 < 0)) & nwt16);
				jcycle12 += wts_idx_inc2;		jcycle12 += ( (-(jcycle12 < 0)) & nwt16);
				jcycle13 += wts_idx_inc2;		jcycle13 += ( (-(jcycle13 < 0)) & nwt16);
				jcycle14 += wts_idx_inc2;		jcycle14 += ( (-(jcycle14 < 0)) & nwt16);

			  #ifdef USE_AVX
				kcycle00 += wts_idx_inc2;		kcycle00 += ( (-(kcycle00 < 0)) & nwt16);
				kcycle01 += wts_idx_inc2;		kcycle01 += ( (-(kcycle01 < 0)) & nwt16);
				kcycle02 += wts_idx_inc2;		kcycle02 += ( (-(kcycle02 < 0)) & nwt16);
				kcycle03 += wts_idx_inc2;		kcycle03 += ( (-(kcycle03 < 0)) & nwt16);
				kcycle04 += wts_idx_inc2;		kcycle04 += ( (-(kcycle04 < 0)) & nwt16);
				kcycle05 += wts_idx_inc2;		kcycle05 += ( (-(kcycle05 < 0)) & nwt16);
				kcycle06 += wts_idx_inc2;		kcycle06 += ( (-(kcycle06 < 0)) & nwt16);
				kcycle07 += wts_idx_inc2;		kcycle07 += ( (-(kcycle07 < 0)) & nwt16);
				kcycle08 += wts_idx_inc2;		kcycle08 += ( (-(kcycle08 < 0)) & nwt16);
				kcycle09 += wts_idx_inc2;		kcycle09 += ( (-(kcycle09 < 0)) & nwt16);
				kcycle10 += wts_idx_inc2;		kcycle10 += ( (-(kcycle10 < 0)) & nwt16);
				kcycle11 += wts_idx_inc2;		kcycle11 += ( (-(kcycle11 < 0)) & nwt16);
				kcycle12 += wts_idx_inc2;		kcycle12 += ( (-(kcycle12 < 0)) & nwt16);
				kcycle13 += wts_idx_inc2;		kcycle13 += ( (-(kcycle13 < 0)) & nwt16);
				kcycle14 += wts_idx_inc2;		kcycle14 += ( (-(kcycle14 < 0)) & nwt16);

				lcycle00 += wts_idx_inc2;		lcycle00 += ( (-(lcycle00 < 0)) & nwt16);
				lcycle01 += wts_idx_inc2;		lcycle01 += ( (-(lcycle01 < 0)) & nwt16);
				lcycle02 += wts_idx_inc2;		lcycle02 += ( (-(lcycle02 < 0)) & nwt16);
				lcycle03 += wts_idx_inc2;		lcycle03 += ( (-(lcycle03 < 0)) & nwt16);
				lcycle04 += wts_idx_inc2;		lcycle04 += ( (-(lcycle04 < 0)) & nwt16);
				lcycle05 += wts_idx_inc2;		lcycle05 += ( (-(lcycle05 < 0)) & nwt16);
				lcycle06 += wts_idx_inc2;		lcycle06 += ( (-(lcycle06 < 0)) & nwt16);
				lcycle07 += wts_idx_inc2;		lcycle07 += ( (-(lcycle07 < 0)) & nwt16);
				lcycle08 += wts_idx_inc2;		lcycle08 += ( (-(lcycle08 < 0)) & nwt16);
				lcycle09 += wts_idx_inc2;		lcycle09 += ( (-(lcycle09 < 0)) & nwt16);
				lcycle10 += wts_idx_inc2;		lcycle10 += ( (-(lcycle10 < 0)) & nwt16);
				lcycle11 += wts_idx_inc2;		lcycle11 += ( (-(lcycle11 < 0)) & nwt16);
				lcycle12 += wts_idx_inc2;		lcycle12 += ( (-(lcycle12 < 0)) & nwt16);
				lcycle13 += wts_idx_inc2;		lcycle13 += ( (-(lcycle13 < 0)) & nwt16);
				lcycle14 += wts_idx_inc2;		lcycle14 += ( (-(lcycle14 < 0)) & nwt16);
			  #endif
			#endif

			}	/* if(MODULUS_TYPE == ...) */

		/*...The radix-60 DIF pass is here:	*/

			#ifdef USE_SSE2

			/* General indexing for radix-15 done as 3 radix-5 followed by 5 radix-3 is
			RADIX_15_DIF(00,01,02,03,04,05,06,07,08,09,0A,0B,0C,0D,0E) ==>

			RADIX_05_DFT(i0,iC,i9,i6,i3, t0,t1,t2,t3,t4)
			RADIX_05_DFT(iA,i7,i4,i1,iD, t5,t6,t7,t8,t9)
			RADIX_05_DFT(i5,i2,iE,iB,i8, tA,tB,tC,tD,tE)
				RADIX_03_DFT(t0,t5,tA, o0,o1,o2,)
				RADIX_03_DFT(t1,t6,tB, oD,oE,oB,)
				RADIX_03_DFT(t2,t7,tC, o9,oA,oB,)
				RADIX_03_DFT(t3,t8,tD, o8,o6,o7,)
				RADIX_03_DFT(t4,t9,tE, o4,o5,o3,)
			*/
			/*
			RADIX_05_DFT(00,0c,19,26,33, r00,r01,r02,r03,r04);
			RADIX_05_DFT(15,22,2e,3b,08, r05,r06,r07,r08,r09);
			RADIX_05_DFT(2a,37,04,11,1d, r0a,r0b,r0c,r0d,r0e);

			RADIX_03_DFT(r00,r05,r0a,,00,01,02,);
			RADIX_03_DFT(r01,r06,r0b,,13,14,12,);
			RADIX_03_DFT(r02,r07,r0c,,09,10,11,);
			RADIX_03_DFT(r03,r08,r0d,,08,06,07,);
			RADIX_03_DFT(r04,r09,r0e,,04,05,03,);
			*/
			/*
			RADIX_05_DFT(30,3c,09,16,23, r10,r11,r12,r13,r14);
			RADIX_05_DFT(05,12,1e,2b,38, r15,r16,r17,r18,r19);
			RADIX_05_DFT(1a,27,34,01,0d, r1a,r1b,r1c,r1d,r1e);
			*/
			/*
			RADIX_05_DFT(20,2c,39,06,13, r20,r21,r22,r23,r24);
			RADIX_05_DFT(35,02,0e,1b,28, r25,r26,r27,r28,r29);
			RADIX_05_DFT(0a,17,24,31,3d, r2a,r2b,r2c,r2d,r2e);
			*/
			/*
			RADIX_05_DFT(10,1c,29,36,03, r30,r31,r32,r33,r34);
			RADIX_05_DFT(25,32,3e,0b,18, r35,r36,r37,r38,r39);
			RADIX_05_DFT(3a,07,14,21,2d, r3a,r3b,r3c,r3d,r3e);
			*/
			  #if (OS_BITS == 32) || !defined(USE_64BIT_ASM_STYLE)	// In 64-bit mode, default is to use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just xmm0-7.
				/* NOTE that to permit us to re-use the s1p-data as both inputs and outputs of the radix-15 DIF DFT, must swap s1p[0123] -> s1p[0321] in output rows: */
			   #ifndef USE_LITERAL_BYTE_OFFSETS
				SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1,        s1p00r,s1p3br,s1p37r,s1p33r,s1p2er,s1p2ar,s1p26r,s1p22r,s1p1dr,s1p19r,s1p15r,s1p11r,s1p0cr,s1p08r,s1p04r, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r0a,r0b,r0c,r0d,r0e);
				SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1,        s1p30r,s1p2br,s1p27r,s1p23r,s1p1er,s1p1ar,s1p16r,s1p12r,s1p0dr,s1p09r,s1p05r,s1p01r,s1p3cr,s1p38r,s1p34r, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r1a,r1b,r1c,r1d,r1e);
				SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1,        s1p20r,s1p1br,s1p17r,s1p13r,s1p0er,s1p0ar,s1p06r,s1p02r,s1p3dr,s1p39r,s1p35r,s1p31r,s1p2cr,s1p28r,s1p24r, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r2a,r2b,r2c,r2d,r2e);
				SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1,        s1p10r,s1p0br,s1p07r,s1p03r,s1p3er,s1p3ar,s1p36r,s1p32r,s1p2dr,s1p29r,s1p25r,s1p21r,s1p1cr,s1p18r,s1p14r, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r3a,r3b,r3c,r3d,r3e);
			   #else	// Versions using base-address-plus-literal-byte-offsets:
				SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1, s1p00r, 0x000, 0x700, 0x680, 0x600, 0x580, 0x500, 0x480, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, r00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0);
				SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1, s1p00r, 0x5a0, 0x520, 0x4a0, 0x420, 0x3a0, 0x320, 0x2a0, 0x220, 0x1a0, 0x120, 0x0a0, 0x020, 0x720, 0x6a0, 0x620, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, r10,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0);
				SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1, s1p00r, 0x3c0, 0x340, 0x2c0, 0x240, 0x1c0, 0x140, 0x0c0, 0x040, 0x740, 0x6c0, 0x640, 0x5c0, 0x540, 0x4c0, 0x440, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, r20,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0);
				SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1, s1p00r, 0x1e0, 0x160, 0x0e0, 0x060, 0x760, 0x6e0, 0x660, 0x5e0, 0x560, 0x4e0, 0x460, 0x3e0, 0x360, 0x2e0, 0x260, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, r30,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0);
			   #endif
			  #else
				SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1, s1p00r,s1p3br,s1p37r,s1p33r,s1p2er,s1p2ar,s1p26r,s1p22r,s1p1dr,s1p19r,s1p15r,s1p11r,s1p0cr,s1p08r,s1p04r, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r0a,r0b,r0c,r0d,r0e, s1p30r,s1p2br,s1p27r,s1p23r,s1p1er,s1p1ar,s1p16r,s1p12r,s1p0dr,s1p09r,s1p05r,s1p01r,s1p3cr,s1p38r,s1p34r, y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e, r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r1a,r1b,r1c,r1d,r1e);
				SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1, s1p20r,s1p1br,s1p17r,s1p13r,s1p0er,s1p0ar,s1p06r,s1p02r,s1p3dr,s1p39r,s1p35r,s1p31r,s1p2cr,s1p28r,s1p24r, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r2a,r2b,r2c,r2d,r2e, s1p10r,s1p0br,s1p07r,s1p03r,s1p3er,s1p3ar,s1p36r,s1p32r,s1p2dr,s1p29r,s1p25r,s1p21r,s1p1cr,s1p18r,s1p14r, y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e, r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r3a,r3b,r3c,r3d,r3e);
			  #endif

			  #ifdef USE_AVX
			// Reorder blocks to yield sequentially increasing a-array offsets:
				add0 = &a[j1    ];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, 0x3c0)
				add3 = &a[j1+p04];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, 0x3c0)
				add1 = &a[j1+p08];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01, 0x3c0)
				add3 = &a[j1+p12];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, 0x3c0)
				add1 = &a[j1+p16];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0d, 0x3c0)
				add2 = &a[j1+p20];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, 0x3c0)
				add1 = &a[j1+p24];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0b, 0x3c0)
				add2 = &a[j1+p28];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, 0x3c0)
				add0 = &a[j1+p32];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09, 0x3c0)
				add2 = &a[j1+p36];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, 0x3c0)
				add0 = &a[j1+p40];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07, 0x3c0)
				add3 = &a[j1+p44];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, 0x3c0)
				add0 = &a[j1+p48];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05, 0x3c0)
				add3 = &a[j1+p52];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, 0x3c0)
				add1 = &a[j1+p56];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03, 0x3c0)

			  #else

			   #if !GCC_ASM_FULL_INLINE
			// Reorder blocks to yield sequentially increasing a-array offsets:
/* Block 01: */	add0 = &a[j1    ];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, 0x1e0)
/* Block 03: */	add3 = &a[j1+p04];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, 0x1e0)
/* Block 02: */	add1 = &a[j1+p08];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01, 0x1e0)
/* Block 15: */	add3 = &a[j1+p12];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, 0x1e0)
/* Block 14: */	add1 = &a[j1+p16];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0d, 0x1e0)
/* Block 13: */	add2 = &a[j1+p20];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, 0x1e0)
/* Block 12: */	add1 = &a[j1+p24];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0b, 0x1e0)
/* Block 11: */	add2 = &a[j1+p28];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, 0x1e0)
/* Block 10: */	add0 = &a[j1+p32];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09, 0x1e0)
/* Block 09: */	add2 = &a[j1+p36];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, 0x1e0)
/* Block 08: */	add0 = &a[j1+p40];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07, 0x1e0)
/* Block 07: */	add3 = &a[j1+p44];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, 0x1e0)
/* Block 06: */	add0 = &a[j1+p48];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05, 0x1e0)
/* Block 05: */	add3 = &a[j1+p52];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, 0x1e0)
/* Block 04: */	add1 = &a[j1+p56];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03, 0x1e0)

			   #else

				add0 = &a[j1    ];
				SSE2_RADIX60_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,r00);

			   #endif

			  #endif	// AVX or SSE2?

			#else	/* !USE_SSE2 */

				RADIX_15_DIF(a1p00r,a1p00i,a1p3br,a1p3bi,a1p37r,a1p37i,a1p33r,a1p33i,a1p2er,a1p2ei,a1p2ar,a1p2ai,a1p26r,a1p26i,a1p22r,a1p22i,a1p1dr,a1p1di,a1p19r,a1p19i,a1p15r,a1p15i,a1p11r,a1p11i,a1p0cr,a1p0ci,a1p08r,a1p08i,a1p04r,a1p04i,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,rt,it);
				RADIX_15_DIF(a1p30r,a1p30i,a1p2br,a1p2bi,a1p27r,a1p27i,a1p23r,a1p23i,a1p1er,a1p1ei,a1p1ar,a1p1ai,a1p16r,a1p16i,a1p12r,a1p12i,a1p0dr,a1p0di,a1p09r,a1p09i,a1p05r,a1p05i,a1p01r,a1p01i,a1p3cr,a1p3ci,a1p38r,a1p38i,a1p34r,a1p34i,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,rt,it);
				RADIX_15_DIF(a1p20r,a1p20i,a1p1br,a1p1bi,a1p17r,a1p17i,a1p13r,a1p13i,a1p0er,a1p0ei,a1p0ar,a1p0ai,a1p06r,a1p06i,a1p02r,a1p02i,a1p3dr,a1p3di,a1p39r,a1p39i,a1p35r,a1p35i,a1p31r,a1p31i,a1p2cr,a1p2ci,a1p28r,a1p28i,a1p24r,a1p24i,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,rt,it);
				RADIX_15_DIF(a1p10r,a1p10i,a1p0br,a1p0bi,a1p07r,a1p07i,a1p03r,a1p03i,a1p3er,a1p3ei,a1p3ar,a1p3ai,a1p36r,a1p36i,a1p32r,a1p32i,a1p2dr,a1p2di,a1p29r,a1p29i,a1p25r,a1p25i,a1p21r,a1p21i,a1p1cr,a1p1ci,a1p18r,a1p18i,a1p14r,a1p14i,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei,rt,it);

				jt = j1    ; jp = j2    ;	RADIX_04_DIF(t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
				jt = j1+p08; jp = j2+p08;	RADIX_04_DIF(t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
				jt = j1+p04; jp = j2+p04;	RADIX_04_DIF(t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
				jt = j1+p56; jp = j2+p56;	RADIX_04_DIF(t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
				jt = j1+p52; jp = j2+p52;	RADIX_04_DIF(t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
				jt = j1+p48; jp = j2+p48;	RADIX_04_DIF(t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
				jt = j1+p44; jp = j2+p44;	RADIX_04_DIF(t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
				jt = j1+p40; jp = j2+p40;	RADIX_04_DIF(t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
				jt = j1+p36; jp = j2+p36;	RADIX_04_DIF(t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
				jt = j1+p32; jp = j2+p32;	RADIX_04_DIF(t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
				jt = j1+p28; jp = j2+p28;	RADIX_04_DIF(t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
				jt = j1+p24; jp = j2+p24;	RADIX_04_DIF(t0br,t0bi,t1br,t1bi,t2br,t2bi,t3br,t3bi,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
				jt = j1+p20; jp = j2+p20;	RADIX_04_DIF(t0cr,t0ci,t1cr,t1ci,t2cr,t2ci,t3cr,t3ci,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
				jt = j1+p16; jp = j2+p16;	RADIX_04_DIF(t0dr,t0di,t1dr,t1di,t2dr,t2di,t3dr,t3di,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
				jt = j1+p12; jp = j2+p12;	RADIX_04_DIF(t0er,t0ei,t1er,t1ei,t2er,t2ei,t3er,t3ei,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);

			#endif	/* if(USE_SSE2) */

			}

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
				jstart += nwt;
				jhi    += nwt;

				col += RADIX;
				co3 -= RADIX;
			}
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX
			_cy_r00[ithread] = cy_r00->d0;	_cy_r01[ithread] = cy_r00->d1;	_cy_r02[ithread] = cy_r00->d2;	_cy_r03[ithread] = cy_r00->d3;
			_cy_r04[ithread] = cy_r04->d0;	_cy_r05[ithread] = cy_r04->d1;	_cy_r06[ithread] = cy_r04->d2;	_cy_r07[ithread] = cy_r04->d3;
			_cy_r08[ithread] = cy_r08->d0;	_cy_r09[ithread] = cy_r08->d1;	_cy_r10[ithread] = cy_r08->d2;	_cy_r11[ithread] = cy_r08->d3;
			_cy_r12[ithread] = cy_r12->d0;	_cy_r13[ithread] = cy_r12->d1;	_cy_r14[ithread] = cy_r12->d2;	_cy_r15[ithread] = cy_r12->d3;
			_cy_r16[ithread] = cy_r16->d0;	_cy_r17[ithread] = cy_r16->d1;	_cy_r18[ithread] = cy_r16->d2;	_cy_r19[ithread] = cy_r16->d3;
			_cy_r20[ithread] = cy_r20->d0;	_cy_r21[ithread] = cy_r20->d1;	_cy_r22[ithread] = cy_r20->d2;	_cy_r23[ithread] = cy_r20->d3;
			_cy_r24[ithread] = cy_r24->d0;	_cy_r25[ithread] = cy_r24->d1;	_cy_r26[ithread] = cy_r24->d2;	_cy_r27[ithread] = cy_r24->d3;
			_cy_r28[ithread] = cy_r28->d0;	_cy_r29[ithread] = cy_r28->d1;	_cy_r30[ithread] = cy_r28->d2;	_cy_r31[ithread] = cy_r28->d3;
			_cy_r32[ithread] = cy_r32->d0;	_cy_r33[ithread] = cy_r32->d1;	_cy_r34[ithread] = cy_r32->d2;	_cy_r35[ithread] = cy_r32->d3;
			_cy_r36[ithread] = cy_r36->d0;	_cy_r37[ithread] = cy_r36->d1;	_cy_r38[ithread] = cy_r36->d2;	_cy_r39[ithread] = cy_r36->d3;
			_cy_r40[ithread] = cy_r40->d0;	_cy_r41[ithread] = cy_r40->d1;	_cy_r42[ithread] = cy_r40->d2;	_cy_r43[ithread] = cy_r40->d3;
			_cy_r44[ithread] = cy_r44->d0;	_cy_r45[ithread] = cy_r44->d1;	_cy_r46[ithread] = cy_r44->d2;	_cy_r47[ithread] = cy_r44->d3;
			_cy_r48[ithread] = cy_r48->d0;	_cy_r49[ithread] = cy_r48->d1;	_cy_r50[ithread] = cy_r48->d2;	_cy_r51[ithread] = cy_r48->d3;
			_cy_r52[ithread] = cy_r52->d0;	_cy_r53[ithread] = cy_r52->d1;	_cy_r54[ithread] = cy_r52->d2;	_cy_r55[ithread] = cy_r52->d3;
			_cy_r56[ithread] = cy_r56->d0;	_cy_r57[ithread] = cy_r56->d1;	_cy_r58[ithread] = cy_r56->d2;	_cy_r59[ithread] = cy_r56->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			_cy_r00[ithread] = cy_r00->d0;	_cy_r01[ithread] = cy_r00->d1;
			_cy_r02[ithread] = cy_r02->d0;	_cy_r03[ithread] = cy_r02->d1;
			_cy_r04[ithread] = cy_r04->d0;	_cy_r05[ithread] = cy_r04->d1;
			_cy_r06[ithread] = cy_r06->d0;	_cy_r07[ithread] = cy_r06->d1;
			_cy_r08[ithread] = cy_r08->d0;	_cy_r09[ithread] = cy_r08->d1;
			_cy_r10[ithread] = cy_r10->d0;	_cy_r11[ithread] = cy_r10->d1;
			_cy_r12[ithread] = cy_r12->d0;	_cy_r13[ithread] = cy_r12->d1;
			_cy_r14[ithread] = cy_r14->d0;	_cy_r15[ithread] = cy_r14->d1;
			_cy_r16[ithread] = cy_r16->d0;	_cy_r17[ithread] = cy_r16->d1;
			_cy_r18[ithread] = cy_r18->d0;	_cy_r19[ithread] = cy_r18->d1;
			_cy_r20[ithread] = cy_r20->d0;	_cy_r21[ithread] = cy_r20->d1;
			_cy_r22[ithread] = cy_r22->d0;	_cy_r23[ithread] = cy_r22->d1;
			_cy_r24[ithread] = cy_r24->d0;	_cy_r25[ithread] = cy_r24->d1;
			_cy_r26[ithread] = cy_r26->d0;	_cy_r27[ithread] = cy_r26->d1;
			_cy_r28[ithread] = cy_r28->d0;	_cy_r29[ithread] = cy_r28->d1;
			_cy_r30[ithread] = cy_r30->d0;	_cy_r31[ithread] = cy_r30->d1;
			_cy_r32[ithread] = cy_r32->d0;	_cy_r33[ithread] = cy_r32->d1;
			_cy_r34[ithread] = cy_r34->d0;	_cy_r35[ithread] = cy_r34->d1;
			_cy_r36[ithread] = cy_r36->d0;	_cy_r37[ithread] = cy_r36->d1;
			_cy_r38[ithread] = cy_r38->d0;	_cy_r39[ithread] = cy_r38->d1;
			_cy_r40[ithread] = cy_r40->d0;	_cy_r41[ithread] = cy_r40->d1;
			_cy_r42[ithread] = cy_r42->d0;	_cy_r43[ithread] = cy_r42->d1;
			_cy_r44[ithread] = cy_r44->d0;	_cy_r45[ithread] = cy_r44->d1;
			_cy_r46[ithread] = cy_r46->d0;	_cy_r47[ithread] = cy_r46->d1;
			_cy_r48[ithread] = cy_r48->d0;	_cy_r49[ithread] = cy_r48->d1;
			_cy_r50[ithread] = cy_r50->d0;	_cy_r51[ithread] = cy_r50->d1;
			_cy_r52[ithread] = cy_r52->d0;	_cy_r53[ithread] = cy_r52->d1;
			_cy_r54[ithread] = cy_r54->d0;	_cy_r55[ithread] = cy_r54->d1;
			_cy_r56[ithread] = cy_r56->d0;	_cy_r57[ithread] = cy_r56->d1;
			_cy_r58[ithread] = cy_r58->d0;	_cy_r59[ithread] = cy_r58->d1;
			maxerr = MAX(max_err->d0,max_err->d1);
		#else
			_cy_r00[ithread] = cy_r00;
			_cy_r01[ithread] = cy_r01;
			_cy_r02[ithread] = cy_r02;
			_cy_r03[ithread] = cy_r03;
			_cy_r04[ithread] = cy_r04;
			_cy_r05[ithread] = cy_r05;
			_cy_r06[ithread] = cy_r06;
			_cy_r07[ithread] = cy_r07;
			_cy_r08[ithread] = cy_r08;
			_cy_r09[ithread] = cy_r09;
			_cy_r10[ithread] = cy_r10;
			_cy_r11[ithread] = cy_r11;
			_cy_r12[ithread] = cy_r12;
			_cy_r13[ithread] = cy_r13;
			_cy_r14[ithread] = cy_r14;
			_cy_r15[ithread] = cy_r15;
			_cy_r16[ithread] = cy_r16;
			_cy_r17[ithread] = cy_r17;
			_cy_r18[ithread] = cy_r18;
			_cy_r19[ithread] = cy_r19;
			_cy_r20[ithread] = cy_r20;
			_cy_r21[ithread] = cy_r21;
			_cy_r22[ithread] = cy_r22;
			_cy_r23[ithread] = cy_r23;
			_cy_r24[ithread] = cy_r24;
			_cy_r25[ithread] = cy_r25;
			_cy_r26[ithread] = cy_r26;
			_cy_r27[ithread] = cy_r27;
			_cy_r28[ithread] = cy_r28;
			_cy_r29[ithread] = cy_r29;
			_cy_r30[ithread] = cy_r30;
			_cy_r31[ithread] = cy_r31;
			_cy_r32[ithread] = cy_r32;
			_cy_r33[ithread] = cy_r33;
			_cy_r34[ithread] = cy_r34;
			_cy_r35[ithread] = cy_r35;
			_cy_r36[ithread] = cy_r36;
			_cy_r37[ithread] = cy_r37;
			_cy_r38[ithread] = cy_r38;
			_cy_r39[ithread] = cy_r39;
			_cy_r40[ithread] = cy_r40;
			_cy_r41[ithread] = cy_r41;
			_cy_r42[ithread] = cy_r42;
			_cy_r43[ithread] = cy_r43;
			_cy_r44[ithread] = cy_r44;
			_cy_r45[ithread] = cy_r45;
			_cy_r46[ithread] = cy_r46;
			_cy_r47[ithread] = cy_r47;
			_cy_r48[ithread] = cy_r48;
			_cy_r49[ithread] = cy_r49;
			_cy_r50[ithread] = cy_r50;
			_cy_r51[ithread] = cy_r51;
			_cy_r52[ithread] = cy_r52;
			_cy_r53[ithread] = cy_r53;
			_cy_r54[ithread] = cy_r54;
			_cy_r55[ithread] = cy_r55;
			_cy_r56[ithread] = cy_r56;
			_cy_r57[ithread] = cy_r57;
			_cy_r58[ithread] = cy_r58;
			_cy_r59[ithread] = cy_r59;
		#endif
		}
		else
		{
		#ifdef USE_AVX
			_cy_r00[ithread] = cy_r00->d0;	_cy_r01[ithread] = cy_r00->d1;	_cy_r02[ithread] = cy_r00->d2;	_cy_r03[ithread] = cy_r00->d3;
			_cy_r04[ithread] = cy_r04->d0;	_cy_r05[ithread] = cy_r04->d1;	_cy_r06[ithread] = cy_r04->d2;	_cy_r07[ithread] = cy_r04->d3;
			_cy_r08[ithread] = cy_r08->d0;	_cy_r09[ithread] = cy_r08->d1;	_cy_r10[ithread] = cy_r08->d2;	_cy_r11[ithread] = cy_r08->d3;
			_cy_r12[ithread] = cy_r12->d0;	_cy_r13[ithread] = cy_r12->d1;	_cy_r14[ithread] = cy_r12->d2;	_cy_r15[ithread] = cy_r12->d3;
			_cy_r16[ithread] = cy_r16->d0;	_cy_r17[ithread] = cy_r16->d1;	_cy_r18[ithread] = cy_r16->d2;	_cy_r19[ithread] = cy_r16->d3;
			_cy_r20[ithread] = cy_r20->d0;	_cy_r21[ithread] = cy_r20->d1;	_cy_r22[ithread] = cy_r20->d2;	_cy_r23[ithread] = cy_r20->d3;
			_cy_r24[ithread] = cy_r24->d0;	_cy_r25[ithread] = cy_r24->d1;	_cy_r26[ithread] = cy_r24->d2;	_cy_r27[ithread] = cy_r24->d3;
			_cy_r28[ithread] = cy_r28->d0;	_cy_r29[ithread] = cy_r28->d1;	_cy_r30[ithread] = cy_r28->d2;	_cy_r31[ithread] = cy_r28->d3;
			_cy_r32[ithread] = cy_r32->d0;	_cy_r33[ithread] = cy_r32->d1;	_cy_r34[ithread] = cy_r32->d2;	_cy_r35[ithread] = cy_r32->d3;
			_cy_r36[ithread] = cy_r36->d0;	_cy_r37[ithread] = cy_r36->d1;	_cy_r38[ithread] = cy_r36->d2;	_cy_r39[ithread] = cy_r36->d3;
			_cy_r40[ithread] = cy_r40->d0;	_cy_r41[ithread] = cy_r40->d1;	_cy_r42[ithread] = cy_r40->d2;	_cy_r43[ithread] = cy_r40->d3;
			_cy_r44[ithread] = cy_r44->d0;	_cy_r45[ithread] = cy_r44->d1;	_cy_r46[ithread] = cy_r44->d2;	_cy_r47[ithread] = cy_r44->d3;
			_cy_r48[ithread] = cy_r48->d0;	_cy_r49[ithread] = cy_r48->d1;	_cy_r50[ithread] = cy_r48->d2;	_cy_r51[ithread] = cy_r48->d3;
			_cy_r52[ithread] = cy_r52->d0;	_cy_r53[ithread] = cy_r52->d1;	_cy_r54[ithread] = cy_r52->d2;	_cy_r55[ithread] = cy_r52->d3;
			_cy_r56[ithread] = cy_r56->d0;	_cy_r57[ithread] = cy_r56->d1;	_cy_r58[ithread] = cy_r56->d2;	_cy_r59[ithread] = cy_r56->d3;

			_cy_i00[ithread] = cy_i00->d0;	_cy_i01[ithread] = cy_i00->d1;	_cy_i02[ithread] = cy_i00->d2;	_cy_i03[ithread] = cy_i00->d3;
			_cy_i04[ithread] = cy_i04->d0;	_cy_i05[ithread] = cy_i04->d1;	_cy_i06[ithread] = cy_i04->d2;	_cy_i07[ithread] = cy_i04->d3;
			_cy_i08[ithread] = cy_i08->d0;	_cy_i09[ithread] = cy_i08->d1;	_cy_i10[ithread] = cy_i08->d2;	_cy_i11[ithread] = cy_i08->d3;
			_cy_i12[ithread] = cy_i12->d0;	_cy_i13[ithread] = cy_i12->d1;	_cy_i14[ithread] = cy_i12->d2;	_cy_i15[ithread] = cy_i12->d3;
			_cy_i16[ithread] = cy_i16->d0;	_cy_i17[ithread] = cy_i16->d1;	_cy_i18[ithread] = cy_i16->d2;	_cy_i19[ithread] = cy_i16->d3;
			_cy_i20[ithread] = cy_i20->d0;	_cy_i21[ithread] = cy_i20->d1;	_cy_i22[ithread] = cy_i20->d2;	_cy_i23[ithread] = cy_i20->d3;
			_cy_i24[ithread] = cy_i24->d0;	_cy_i25[ithread] = cy_i24->d1;	_cy_i26[ithread] = cy_i24->d2;	_cy_i27[ithread] = cy_i24->d3;
			_cy_i28[ithread] = cy_i28->d0;	_cy_i29[ithread] = cy_i28->d1;	_cy_i30[ithread] = cy_i28->d2;	_cy_i31[ithread] = cy_i28->d3;
			_cy_i32[ithread] = cy_i32->d0;	_cy_i33[ithread] = cy_i32->d1;	_cy_i34[ithread] = cy_i32->d2;	_cy_i35[ithread] = cy_i32->d3;
			_cy_i36[ithread] = cy_i36->d0;	_cy_i37[ithread] = cy_i36->d1;	_cy_i38[ithread] = cy_i36->d2;	_cy_i39[ithread] = cy_i36->d3;
			_cy_i40[ithread] = cy_i40->d0;	_cy_i41[ithread] = cy_i40->d1;	_cy_i42[ithread] = cy_i40->d2;	_cy_i43[ithread] = cy_i40->d3;
			_cy_i44[ithread] = cy_i44->d0;	_cy_i45[ithread] = cy_i44->d1;	_cy_i46[ithread] = cy_i44->d2;	_cy_i47[ithread] = cy_i44->d3;
			_cy_i48[ithread] = cy_i48->d0;	_cy_i49[ithread] = cy_i48->d1;	_cy_i50[ithread] = cy_i48->d2;	_cy_i51[ithread] = cy_i48->d3;
			_cy_i52[ithread] = cy_i52->d0;	_cy_i53[ithread] = cy_i52->d1;	_cy_i54[ithread] = cy_i52->d2;	_cy_i55[ithread] = cy_i52->d3;
			_cy_i56[ithread] = cy_i56->d0;	_cy_i57[ithread] = cy_i56->d1;	_cy_i58[ithread] = cy_i56->d2;	_cy_i59[ithread] = cy_i56->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
		#elif defined(USE_SSE2)
			// Carry pattern for Fermat-mod in SSE2 mode is kinda funky:
			_cy_r00[ithread] = cy_r00->d0;	_cy_i00[ithread] = cy_r00->d1;
			_cy_r01[ithread] = cy_r02->d0;	_cy_i01[ithread] = cy_r02->d1;
			_cy_r02[ithread] = cy_r04->d0;	_cy_i02[ithread] = cy_r04->d1;
			_cy_r03[ithread] = cy_r06->d0;	_cy_i03[ithread] = cy_r06->d1;
			_cy_r04[ithread] = cy_r08->d0;	_cy_i04[ithread] = cy_r08->d1;
			_cy_r05[ithread] = cy_r10->d0;	_cy_i05[ithread] = cy_r10->d1;
			_cy_r06[ithread] = cy_r12->d0;	_cy_i06[ithread] = cy_r12->d1;
			_cy_r07[ithread] = cy_r14->d0;	_cy_i07[ithread] = cy_r14->d1;
			_cy_r08[ithread] = cy_r16->d0;	_cy_i08[ithread] = cy_r16->d1;
			_cy_r09[ithread] = cy_r18->d0;	_cy_i09[ithread] = cy_r18->d1;
			_cy_r10[ithread] = cy_r20->d0;	_cy_i10[ithread] = cy_r20->d1;
			_cy_r11[ithread] = cy_r22->d0;	_cy_i11[ithread] = cy_r22->d1;
			_cy_r12[ithread] = cy_r24->d0;	_cy_i12[ithread] = cy_r24->d1;
			_cy_r13[ithread] = cy_r26->d0;	_cy_i13[ithread] = cy_r26->d1;
			_cy_r14[ithread] = cy_r28->d0;	_cy_i14[ithread] = cy_r28->d1;
			_cy_r15[ithread] = cy_r30->d0;	_cy_i15[ithread] = cy_r30->d1;
			_cy_r16[ithread] = cy_r32->d0;	_cy_i16[ithread] = cy_r32->d1;
			_cy_r17[ithread] = cy_r34->d0;	_cy_i17[ithread] = cy_r34->d1;
			_cy_r18[ithread] = cy_r36->d0;	_cy_i18[ithread] = cy_r36->d1;
			_cy_r19[ithread] = cy_r38->d0;	_cy_i19[ithread] = cy_r38->d1;
			_cy_r20[ithread] = cy_r40->d0;	_cy_i20[ithread] = cy_r40->d1;
			_cy_r21[ithread] = cy_r42->d0;	_cy_i21[ithread] = cy_r42->d1;
			_cy_r22[ithread] = cy_r44->d0;	_cy_i22[ithread] = cy_r44->d1;
			_cy_r23[ithread] = cy_r46->d0;	_cy_i23[ithread] = cy_r46->d1;
			_cy_r24[ithread] = cy_r48->d0;	_cy_i24[ithread] = cy_r48->d1;
			_cy_r25[ithread] = cy_r50->d0;	_cy_i25[ithread] = cy_r50->d1;
			_cy_r26[ithread] = cy_r52->d0;	_cy_i26[ithread] = cy_r52->d1;
			_cy_r27[ithread] = cy_r54->d0;	_cy_i27[ithread] = cy_r54->d1;
			_cy_r28[ithread] = cy_r56->d0;	_cy_i28[ithread] = cy_r56->d1;
			_cy_r29[ithread] = cy_r58->d0;	_cy_i29[ithread] = cy_r58->d1;
			_cy_r30[ithread] = cy_i00->d0;	_cy_i30[ithread] = cy_i00->d1;
			_cy_r31[ithread] = cy_i02->d0;	_cy_i31[ithread] = cy_i02->d1;
			_cy_r32[ithread] = cy_i04->d0;	_cy_i32[ithread] = cy_i04->d1;
			_cy_r33[ithread] = cy_i06->d0;	_cy_i33[ithread] = cy_i06->d1;
			_cy_r34[ithread] = cy_i08->d0;	_cy_i34[ithread] = cy_i08->d1;
			_cy_r35[ithread] = cy_i10->d0;	_cy_i35[ithread] = cy_i10->d1;
			_cy_r36[ithread] = cy_i12->d0;	_cy_i36[ithread] = cy_i12->d1;
			_cy_r37[ithread] = cy_i14->d0;	_cy_i37[ithread] = cy_i14->d1;
			_cy_r38[ithread] = cy_i16->d0;	_cy_i38[ithread] = cy_i16->d1;
			_cy_r39[ithread] = cy_i18->d0;	_cy_i39[ithread] = cy_i18->d1;
			_cy_r40[ithread] = cy_i20->d0;	_cy_i40[ithread] = cy_i20->d1;
			_cy_r41[ithread] = cy_i22->d0;	_cy_i41[ithread] = cy_i22->d1;
			_cy_r42[ithread] = cy_i24->d0;	_cy_i42[ithread] = cy_i24->d1;
			_cy_r43[ithread] = cy_i26->d0;	_cy_i43[ithread] = cy_i26->d1;
			_cy_r44[ithread] = cy_i28->d0;	_cy_i44[ithread] = cy_i28->d1;
			_cy_r45[ithread] = cy_i30->d0;	_cy_i45[ithread] = cy_i30->d1;
			_cy_r46[ithread] = cy_i32->d0;	_cy_i46[ithread] = cy_i32->d1;
			_cy_r47[ithread] = cy_i34->d0;	_cy_i47[ithread] = cy_i34->d1;
			_cy_r48[ithread] = cy_i36->d0;	_cy_i48[ithread] = cy_i36->d1;
			_cy_r49[ithread] = cy_i38->d0;	_cy_i49[ithread] = cy_i38->d1;
			_cy_r50[ithread] = cy_i40->d0;	_cy_i50[ithread] = cy_i40->d1;
			_cy_r51[ithread] = cy_i42->d0;	_cy_i51[ithread] = cy_i42->d1;
			_cy_r52[ithread] = cy_i44->d0;	_cy_i52[ithread] = cy_i44->d1;
			_cy_r53[ithread] = cy_i46->d0;	_cy_i53[ithread] = cy_i46->d1;
			_cy_r54[ithread] = cy_i48->d0;	_cy_i54[ithread] = cy_i48->d1;
			_cy_r55[ithread] = cy_i50->d0;	_cy_i55[ithread] = cy_i50->d1;
			_cy_r56[ithread] = cy_i52->d0;	_cy_i56[ithread] = cy_i52->d1;
			_cy_r57[ithread] = cy_i54->d0;	_cy_i57[ithread] = cy_i54->d1;
			_cy_r58[ithread] = cy_i56->d0;	_cy_i58[ithread] = cy_i56->d1;
			_cy_r59[ithread] = cy_i58->d0;	_cy_i59[ithread] = cy_i58->d1;
			maxerr = MAX(max_err->d0,max_err->d1);
		#else
			_cy_r00[ithread] = cy_r00;	_cy_i00[ithread] = cy_i00;
			_cy_r01[ithread] = cy_r01;	_cy_i01[ithread] = cy_i01;
			_cy_r02[ithread] = cy_r02;	_cy_i02[ithread] = cy_i02;
			_cy_r03[ithread] = cy_r03;	_cy_i03[ithread] = cy_i03;
			_cy_r04[ithread] = cy_r04;	_cy_i04[ithread] = cy_i04;
			_cy_r05[ithread] = cy_r05;	_cy_i05[ithread] = cy_i05;
			_cy_r06[ithread] = cy_r06;	_cy_i06[ithread] = cy_i06;
			_cy_r07[ithread] = cy_r07;	_cy_i07[ithread] = cy_i07;
			_cy_r08[ithread] = cy_r08;	_cy_i08[ithread] = cy_i08;
			_cy_r09[ithread] = cy_r09;	_cy_i09[ithread] = cy_i09;
			_cy_r10[ithread] = cy_r10;	_cy_i10[ithread] = cy_i10;
			_cy_r11[ithread] = cy_r11;	_cy_i11[ithread] = cy_i11;
			_cy_r12[ithread] = cy_r12;	_cy_i12[ithread] = cy_i12;
			_cy_r13[ithread] = cy_r13;	_cy_i13[ithread] = cy_i13;
			_cy_r14[ithread] = cy_r14;	_cy_i14[ithread] = cy_i14;
			_cy_r15[ithread] = cy_r15;	_cy_i15[ithread] = cy_i15;
			_cy_r16[ithread] = cy_r16;	_cy_i16[ithread] = cy_i16;
			_cy_r17[ithread] = cy_r17;	_cy_i17[ithread] = cy_i17;
			_cy_r18[ithread] = cy_r18;	_cy_i18[ithread] = cy_i18;
			_cy_r19[ithread] = cy_r19;	_cy_i19[ithread] = cy_i19;
			_cy_r20[ithread] = cy_r20;	_cy_i20[ithread] = cy_i20;
			_cy_r21[ithread] = cy_r21;	_cy_i21[ithread] = cy_i21;
			_cy_r22[ithread] = cy_r22;	_cy_i22[ithread] = cy_i22;
			_cy_r23[ithread] = cy_r23;	_cy_i23[ithread] = cy_i23;
			_cy_r24[ithread] = cy_r24;	_cy_i24[ithread] = cy_i24;
			_cy_r25[ithread] = cy_r25;	_cy_i25[ithread] = cy_i25;
			_cy_r26[ithread] = cy_r26;	_cy_i26[ithread] = cy_i26;
			_cy_r27[ithread] = cy_r27;	_cy_i27[ithread] = cy_i27;
			_cy_r28[ithread] = cy_r28;	_cy_i28[ithread] = cy_i28;
			_cy_r29[ithread] = cy_r29;	_cy_i29[ithread] = cy_i29;
			_cy_r30[ithread] = cy_r30;	_cy_i30[ithread] = cy_i30;
			_cy_r31[ithread] = cy_r31;	_cy_i31[ithread] = cy_i31;
			_cy_r32[ithread] = cy_r32;	_cy_i32[ithread] = cy_i32;
			_cy_r33[ithread] = cy_r33;	_cy_i33[ithread] = cy_i33;
			_cy_r34[ithread] = cy_r34;	_cy_i34[ithread] = cy_i34;
			_cy_r35[ithread] = cy_r35;	_cy_i35[ithread] = cy_i35;
			_cy_r36[ithread] = cy_r36;	_cy_i36[ithread] = cy_i36;
			_cy_r37[ithread] = cy_r37;	_cy_i37[ithread] = cy_i37;
			_cy_r38[ithread] = cy_r38;	_cy_i38[ithread] = cy_i38;
			_cy_r39[ithread] = cy_r39;	_cy_i39[ithread] = cy_i39;
			_cy_r40[ithread] = cy_r40;	_cy_i40[ithread] = cy_i40;
			_cy_r41[ithread] = cy_r41;	_cy_i41[ithread] = cy_i41;
			_cy_r42[ithread] = cy_r42;	_cy_i42[ithread] = cy_i42;
			_cy_r43[ithread] = cy_r43;	_cy_i43[ithread] = cy_i43;
			_cy_r44[ithread] = cy_r44;	_cy_i44[ithread] = cy_i44;
			_cy_r45[ithread] = cy_r45;	_cy_i45[ithread] = cy_i45;
			_cy_r46[ithread] = cy_r46;	_cy_i46[ithread] = cy_i46;
			_cy_r47[ithread] = cy_r47;	_cy_i47[ithread] = cy_i47;
			_cy_r48[ithread] = cy_r48;	_cy_i48[ithread] = cy_i48;
			_cy_r49[ithread] = cy_r49;	_cy_i49[ithread] = cy_i49;
			_cy_r50[ithread] = cy_r50;	_cy_i50[ithread] = cy_i50;
			_cy_r51[ithread] = cy_r51;	_cy_i51[ithread] = cy_i51;
			_cy_r52[ithread] = cy_r52;	_cy_i52[ithread] = cy_i52;
			_cy_r53[ithread] = cy_r53;	_cy_i53[ithread] = cy_i53;
			_cy_r54[ithread] = cy_r54;	_cy_i54[ithread] = cy_i54;
			_cy_r55[ithread] = cy_r55;	_cy_i55[ithread] = cy_i55;
			_cy_r56[ithread] = cy_r56;	_cy_i56[ithread] = cy_i56;
			_cy_r57[ithread] = cy_r57;	_cy_i57[ithread] = cy_i57;
			_cy_r58[ithread] = cy_r58;	_cy_i58[ithread] = cy_i58;
			_cy_r59[ithread] = cy_r59;	_cy_i59[ithread] = cy_i59;
		#endif
		}

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(_maxerr[ithread] < maxerr)
		{
			_maxerr[ithread] = maxerr;
		}

  #endif	// #ifdef USE_PTHREAD

	}	/******* END OF PARALLEL FOR-LOOP ********/

#ifdef USE_PTHREAD	// End of threadpool-based dispatch: Add a small wait-loop to ensure all threads complete

  #ifdef OS_TYPE_MACOSX

	/*** Main execution thread executes remaining chunks in serial fashion (but in || with the pool threads): ***/
	for(j = 0; j < main_work_units; ++j)
	{
	//	printf("adding main task %d\n",j + pool_work_units);
		ASSERT(HERE, 0x0 == cy60_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;
	ns_time.tv_sec  = 0.0001;// (time_t)seconds
	ns_time.tv_nsec = 0;	// (long)nanoseconds - At least allegedly, but under OS X it seems to be finer-grained than that

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(HERE, 0 == nanosleep(&ns_time, 0x0), "nanosleep fail!");
	}
//	printf("radix32_ditN_cy_dif1 end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = tdat[ithread].maxerr;
		if(maxerr < _maxerr[ithread]) {
			maxerr = _maxerr[ithread];
		}

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			_cy_r00[ithread] = tdat[ithread].cy_r00;
			_cy_r01[ithread] = tdat[ithread].cy_r01;
			_cy_r02[ithread] = tdat[ithread].cy_r02;
			_cy_r03[ithread] = tdat[ithread].cy_r03;
			_cy_r04[ithread] = tdat[ithread].cy_r04;
			_cy_r05[ithread] = tdat[ithread].cy_r05;
			_cy_r06[ithread] = tdat[ithread].cy_r06;
			_cy_r07[ithread] = tdat[ithread].cy_r07;
			_cy_r08[ithread] = tdat[ithread].cy_r08;
			_cy_r09[ithread] = tdat[ithread].cy_r09;
			_cy_r10[ithread] = tdat[ithread].cy_r10;
			_cy_r11[ithread] = tdat[ithread].cy_r11;
			_cy_r12[ithread] = tdat[ithread].cy_r12;
			_cy_r13[ithread] = tdat[ithread].cy_r13;
			_cy_r14[ithread] = tdat[ithread].cy_r14;
			_cy_r15[ithread] = tdat[ithread].cy_r15;
			_cy_r16[ithread] = tdat[ithread].cy_r16;
			_cy_r17[ithread] = tdat[ithread].cy_r17;
			_cy_r18[ithread] = tdat[ithread].cy_r18;
			_cy_r19[ithread] = tdat[ithread].cy_r19;
			_cy_r20[ithread] = tdat[ithread].cy_r20;
			_cy_r21[ithread] = tdat[ithread].cy_r21;
			_cy_r22[ithread] = tdat[ithread].cy_r22;
			_cy_r23[ithread] = tdat[ithread].cy_r23;
			_cy_r24[ithread] = tdat[ithread].cy_r24;
			_cy_r25[ithread] = tdat[ithread].cy_r25;
			_cy_r26[ithread] = tdat[ithread].cy_r26;
			_cy_r27[ithread] = tdat[ithread].cy_r27;
			_cy_r28[ithread] = tdat[ithread].cy_r28;
			_cy_r29[ithread] = tdat[ithread].cy_r29;
			_cy_r30[ithread] = tdat[ithread].cy_r30;
			_cy_r31[ithread] = tdat[ithread].cy_r31;
			_cy_r32[ithread] = tdat[ithread].cy_r32;
			_cy_r33[ithread] = tdat[ithread].cy_r33;
			_cy_r34[ithread] = tdat[ithread].cy_r34;
			_cy_r35[ithread] = tdat[ithread].cy_r35;
			_cy_r36[ithread] = tdat[ithread].cy_r36;
			_cy_r37[ithread] = tdat[ithread].cy_r37;
			_cy_r38[ithread] = tdat[ithread].cy_r38;
			_cy_r39[ithread] = tdat[ithread].cy_r39;
			_cy_r40[ithread] = tdat[ithread].cy_r40;
			_cy_r41[ithread] = tdat[ithread].cy_r41;
			_cy_r42[ithread] = tdat[ithread].cy_r42;
			_cy_r43[ithread] = tdat[ithread].cy_r43;
			_cy_r44[ithread] = tdat[ithread].cy_r44;
			_cy_r45[ithread] = tdat[ithread].cy_r45;
			_cy_r46[ithread] = tdat[ithread].cy_r46;
			_cy_r47[ithread] = tdat[ithread].cy_r47;
			_cy_r48[ithread] = tdat[ithread].cy_r48;
			_cy_r49[ithread] = tdat[ithread].cy_r49;
			_cy_r50[ithread] = tdat[ithread].cy_r50;
			_cy_r51[ithread] = tdat[ithread].cy_r51;
			_cy_r52[ithread] = tdat[ithread].cy_r52;
			_cy_r53[ithread] = tdat[ithread].cy_r53;
			_cy_r54[ithread] = tdat[ithread].cy_r54;
			_cy_r55[ithread] = tdat[ithread].cy_r55;
			_cy_r56[ithread] = tdat[ithread].cy_r56;
			_cy_r57[ithread] = tdat[ithread].cy_r57;
			_cy_r58[ithread] = tdat[ithread].cy_r58;
			_cy_r59[ithread] = tdat[ithread].cy_r59;
		}
		else
		{
			_cy_r00[ithread] = tdat[ithread].cy_r00;	_cy_i00[ithread] = tdat[ithread].cy_i00;
			_cy_r01[ithread] = tdat[ithread].cy_r01;	_cy_i01[ithread] = tdat[ithread].cy_i01;
			_cy_r02[ithread] = tdat[ithread].cy_r02;	_cy_i02[ithread] = tdat[ithread].cy_i02;
			_cy_r03[ithread] = tdat[ithread].cy_r03;	_cy_i03[ithread] = tdat[ithread].cy_i03;
			_cy_r04[ithread] = tdat[ithread].cy_r04;	_cy_i04[ithread] = tdat[ithread].cy_i04;
			_cy_r05[ithread] = tdat[ithread].cy_r05;	_cy_i05[ithread] = tdat[ithread].cy_i05;
			_cy_r06[ithread] = tdat[ithread].cy_r06;	_cy_i06[ithread] = tdat[ithread].cy_i06;
			_cy_r07[ithread] = tdat[ithread].cy_r07;	_cy_i07[ithread] = tdat[ithread].cy_i07;
			_cy_r08[ithread] = tdat[ithread].cy_r08;	_cy_i08[ithread] = tdat[ithread].cy_i08;
			_cy_r09[ithread] = tdat[ithread].cy_r09;	_cy_i09[ithread] = tdat[ithread].cy_i09;
			_cy_r10[ithread] = tdat[ithread].cy_r10;	_cy_i10[ithread] = tdat[ithread].cy_i10;
			_cy_r11[ithread] = tdat[ithread].cy_r11;	_cy_i11[ithread] = tdat[ithread].cy_i11;
			_cy_r12[ithread] = tdat[ithread].cy_r12;	_cy_i12[ithread] = tdat[ithread].cy_i12;
			_cy_r13[ithread] = tdat[ithread].cy_r13;	_cy_i13[ithread] = tdat[ithread].cy_i13;
			_cy_r14[ithread] = tdat[ithread].cy_r14;	_cy_i14[ithread] = tdat[ithread].cy_i14;
			_cy_r15[ithread] = tdat[ithread].cy_r15;	_cy_i15[ithread] = tdat[ithread].cy_i15;
			_cy_r16[ithread] = tdat[ithread].cy_r16;	_cy_i16[ithread] = tdat[ithread].cy_i16;
			_cy_r17[ithread] = tdat[ithread].cy_r17;	_cy_i17[ithread] = tdat[ithread].cy_i17;
			_cy_r18[ithread] = tdat[ithread].cy_r18;	_cy_i18[ithread] = tdat[ithread].cy_i18;
			_cy_r19[ithread] = tdat[ithread].cy_r19;	_cy_i19[ithread] = tdat[ithread].cy_i19;
			_cy_r20[ithread] = tdat[ithread].cy_r20;	_cy_i20[ithread] = tdat[ithread].cy_i20;
			_cy_r21[ithread] = tdat[ithread].cy_r21;	_cy_i21[ithread] = tdat[ithread].cy_i21;
			_cy_r22[ithread] = tdat[ithread].cy_r22;	_cy_i22[ithread] = tdat[ithread].cy_i22;
			_cy_r23[ithread] = tdat[ithread].cy_r23;	_cy_i23[ithread] = tdat[ithread].cy_i23;
			_cy_r24[ithread] = tdat[ithread].cy_r24;	_cy_i24[ithread] = tdat[ithread].cy_i24;
			_cy_r25[ithread] = tdat[ithread].cy_r25;	_cy_i25[ithread] = tdat[ithread].cy_i25;
			_cy_r26[ithread] = tdat[ithread].cy_r26;	_cy_i26[ithread] = tdat[ithread].cy_i26;
			_cy_r27[ithread] = tdat[ithread].cy_r27;	_cy_i27[ithread] = tdat[ithread].cy_i27;
			_cy_r28[ithread] = tdat[ithread].cy_r28;	_cy_i28[ithread] = tdat[ithread].cy_i28;
			_cy_r29[ithread] = tdat[ithread].cy_r29;	_cy_i29[ithread] = tdat[ithread].cy_i29;
			_cy_r30[ithread] = tdat[ithread].cy_r30;	_cy_i30[ithread] = tdat[ithread].cy_i30;
			_cy_r31[ithread] = tdat[ithread].cy_r31;	_cy_i31[ithread] = tdat[ithread].cy_i31;
			_cy_r32[ithread] = tdat[ithread].cy_r32;	_cy_i32[ithread] = tdat[ithread].cy_i32;
			_cy_r33[ithread] = tdat[ithread].cy_r33;	_cy_i33[ithread] = tdat[ithread].cy_i33;
			_cy_r34[ithread] = tdat[ithread].cy_r34;	_cy_i34[ithread] = tdat[ithread].cy_i34;
			_cy_r35[ithread] = tdat[ithread].cy_r35;	_cy_i35[ithread] = tdat[ithread].cy_i35;
			_cy_r36[ithread] = tdat[ithread].cy_r36;	_cy_i36[ithread] = tdat[ithread].cy_i36;
			_cy_r37[ithread] = tdat[ithread].cy_r37;	_cy_i37[ithread] = tdat[ithread].cy_i37;
			_cy_r38[ithread] = tdat[ithread].cy_r38;	_cy_i38[ithread] = tdat[ithread].cy_i38;
			_cy_r39[ithread] = tdat[ithread].cy_r39;	_cy_i39[ithread] = tdat[ithread].cy_i39;
			_cy_r40[ithread] = tdat[ithread].cy_r40;	_cy_i40[ithread] = tdat[ithread].cy_i40;
			_cy_r41[ithread] = tdat[ithread].cy_r41;	_cy_i41[ithread] = tdat[ithread].cy_i41;
			_cy_r42[ithread] = tdat[ithread].cy_r42;	_cy_i42[ithread] = tdat[ithread].cy_i42;
			_cy_r43[ithread] = tdat[ithread].cy_r43;	_cy_i43[ithread] = tdat[ithread].cy_i43;
			_cy_r44[ithread] = tdat[ithread].cy_r44;	_cy_i44[ithread] = tdat[ithread].cy_i44;
			_cy_r45[ithread] = tdat[ithread].cy_r45;	_cy_i45[ithread] = tdat[ithread].cy_i45;
			_cy_r46[ithread] = tdat[ithread].cy_r46;	_cy_i46[ithread] = tdat[ithread].cy_i46;
			_cy_r47[ithread] = tdat[ithread].cy_r47;	_cy_i47[ithread] = tdat[ithread].cy_i47;
			_cy_r48[ithread] = tdat[ithread].cy_r48;	_cy_i48[ithread] = tdat[ithread].cy_i48;
			_cy_r49[ithread] = tdat[ithread].cy_r49;	_cy_i49[ithread] = tdat[ithread].cy_i49;
			_cy_r50[ithread] = tdat[ithread].cy_r50;	_cy_i50[ithread] = tdat[ithread].cy_i50;
			_cy_r51[ithread] = tdat[ithread].cy_r51;	_cy_i51[ithread] = tdat[ithread].cy_i51;
			_cy_r52[ithread] = tdat[ithread].cy_r52;	_cy_i52[ithread] = tdat[ithread].cy_i52;
			_cy_r53[ithread] = tdat[ithread].cy_r53;	_cy_i53[ithread] = tdat[ithread].cy_i53;
			_cy_r54[ithread] = tdat[ithread].cy_r54;	_cy_i54[ithread] = tdat[ithread].cy_i54;
			_cy_r55[ithread] = tdat[ithread].cy_r55;	_cy_i55[ithread] = tdat[ithread].cy_i55;
			_cy_r56[ithread] = tdat[ithread].cy_r56;	_cy_i56[ithread] = tdat[ithread].cy_i56;
			_cy_r57[ithread] = tdat[ithread].cy_r57;	_cy_i57[ithread] = tdat[ithread].cy_i57;
			_cy_r58[ithread] = tdat[ithread].cy_r58;	_cy_i58[ithread] = tdat[ithread].cy_i58;
			_cy_r59[ithread] = tdat[ithread].cy_r59;	_cy_i59[ithread] = tdat[ithread].cy_i59;
		}

	}
#endif

	if(full_pass) {
	//	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

	/*   Wraparound carry cleanup loop is here: ***
	!
	!   (1) Invert the radix-28 forward DIF FFT of the first block of 28 complex elements in A and unweight;
	!   (2) Propagate cleanup carries among the real and imaginary parts of the 28 outputs of (1);
	!   (3) Reweight and perform a radix-28 forward DIF FFT on the result of (2);
	!   (4) If any of the exit carries from (2) are nonzero, advance to the next 28 elements and repeat (1-4).
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t00r = _cy_r00[CY_THREADS - 1];
		t01r = _cy_r01[CY_THREADS - 1];
		t02r = _cy_r02[CY_THREADS - 1];
		t03r = _cy_r03[CY_THREADS - 1];
		t04r = _cy_r04[CY_THREADS - 1];
		t05r = _cy_r05[CY_THREADS - 1];
		t06r = _cy_r06[CY_THREADS - 1];
		t07r = _cy_r07[CY_THREADS - 1];
		t08r = _cy_r08[CY_THREADS - 1];
		t09r = _cy_r09[CY_THREADS - 1];
		t0ar = _cy_r10[CY_THREADS - 1];
		t0br = _cy_r11[CY_THREADS - 1];
		t0cr = _cy_r12[CY_THREADS - 1];
		t0dr = _cy_r13[CY_THREADS - 1];
		t0er = _cy_r14[CY_THREADS - 1];
		t10r = _cy_r15[CY_THREADS - 1];
		t11r = _cy_r16[CY_THREADS - 1];
		t12r = _cy_r17[CY_THREADS - 1];
		t13r = _cy_r18[CY_THREADS - 1];
		t14r = _cy_r19[CY_THREADS - 1];
		t15r = _cy_r20[CY_THREADS - 1];
		t16r = _cy_r21[CY_THREADS - 1];
		t17r = _cy_r22[CY_THREADS - 1];
		t18r = _cy_r23[CY_THREADS - 1];
		t19r = _cy_r24[CY_THREADS - 1];
		t1ar = _cy_r25[CY_THREADS - 1];
		t1br = _cy_r26[CY_THREADS - 1];
		t1cr = _cy_r27[CY_THREADS - 1];
		t1dr = _cy_r28[CY_THREADS - 1];
		t1er = _cy_r29[CY_THREADS - 1];
		t20r = _cy_r30[CY_THREADS - 1];
		t21r = _cy_r31[CY_THREADS - 1];
		t22r = _cy_r32[CY_THREADS - 1];
		t23r = _cy_r33[CY_THREADS - 1];
		t24r = _cy_r34[CY_THREADS - 1];
		t25r = _cy_r35[CY_THREADS - 1];
		t26r = _cy_r36[CY_THREADS - 1];
		t27r = _cy_r37[CY_THREADS - 1];
		t28r = _cy_r38[CY_THREADS - 1];
		t29r = _cy_r39[CY_THREADS - 1];
		t2ar = _cy_r40[CY_THREADS - 1];
		t2br = _cy_r41[CY_THREADS - 1];
		t2cr = _cy_r42[CY_THREADS - 1];
		t2dr = _cy_r43[CY_THREADS - 1];
		t2er = _cy_r44[CY_THREADS - 1];
		t30r = _cy_r45[CY_THREADS - 1];
		t31r = _cy_r46[CY_THREADS - 1];
		t32r = _cy_r47[CY_THREADS - 1];
		t33r = _cy_r48[CY_THREADS - 1];
		t34r = _cy_r49[CY_THREADS - 1];
		t35r = _cy_r50[CY_THREADS - 1];
		t36r = _cy_r51[CY_THREADS - 1];
		t37r = _cy_r52[CY_THREADS - 1];
		t38r = _cy_r53[CY_THREADS - 1];
		t39r = _cy_r54[CY_THREADS - 1];
		t3ar = _cy_r55[CY_THREADS - 1];
		t3br = _cy_r56[CY_THREADS - 1];
		t3cr = _cy_r57[CY_THREADS - 1];
		t3dr = _cy_r58[CY_THREADS - 1];
		t3er = _cy_r59[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1, "CY_THREADS must be > 1!");	/* Make sure loop only gets executed if multiple threads */
			_cy_r00[ithread] = _cy_r00[ithread-1];
			_cy_r01[ithread] = _cy_r01[ithread-1];
			_cy_r02[ithread] = _cy_r02[ithread-1];
			_cy_r03[ithread] = _cy_r03[ithread-1];
			_cy_r04[ithread] = _cy_r04[ithread-1];
			_cy_r05[ithread] = _cy_r05[ithread-1];
			_cy_r06[ithread] = _cy_r06[ithread-1];
			_cy_r07[ithread] = _cy_r07[ithread-1];
			_cy_r08[ithread] = _cy_r08[ithread-1];
			_cy_r09[ithread] = _cy_r09[ithread-1];
			_cy_r10[ithread] = _cy_r10[ithread-1];
			_cy_r11[ithread] = _cy_r11[ithread-1];
			_cy_r12[ithread] = _cy_r12[ithread-1];
			_cy_r13[ithread] = _cy_r13[ithread-1];
			_cy_r14[ithread] = _cy_r14[ithread-1];
			_cy_r15[ithread] = _cy_r15[ithread-1];
			_cy_r16[ithread] = _cy_r16[ithread-1];
			_cy_r17[ithread] = _cy_r17[ithread-1];
			_cy_r18[ithread] = _cy_r18[ithread-1];
			_cy_r19[ithread] = _cy_r19[ithread-1];
			_cy_r20[ithread] = _cy_r20[ithread-1];
			_cy_r21[ithread] = _cy_r21[ithread-1];
			_cy_r22[ithread] = _cy_r22[ithread-1];
			_cy_r23[ithread] = _cy_r23[ithread-1];
			_cy_r24[ithread] = _cy_r24[ithread-1];
			_cy_r25[ithread] = _cy_r25[ithread-1];
			_cy_r26[ithread] = _cy_r26[ithread-1];
			_cy_r27[ithread] = _cy_r27[ithread-1];
			_cy_r28[ithread] = _cy_r28[ithread-1];
			_cy_r29[ithread] = _cy_r29[ithread-1];
			_cy_r30[ithread] = _cy_r30[ithread-1];
			_cy_r31[ithread] = _cy_r31[ithread-1];
			_cy_r32[ithread] = _cy_r32[ithread-1];
			_cy_r33[ithread] = _cy_r33[ithread-1];
			_cy_r34[ithread] = _cy_r34[ithread-1];
			_cy_r35[ithread] = _cy_r35[ithread-1];
			_cy_r36[ithread] = _cy_r36[ithread-1];
			_cy_r37[ithread] = _cy_r37[ithread-1];
			_cy_r38[ithread] = _cy_r38[ithread-1];
			_cy_r39[ithread] = _cy_r39[ithread-1];
			_cy_r40[ithread] = _cy_r40[ithread-1];
			_cy_r41[ithread] = _cy_r41[ithread-1];
			_cy_r42[ithread] = _cy_r42[ithread-1];
			_cy_r43[ithread] = _cy_r43[ithread-1];
			_cy_r44[ithread] = _cy_r44[ithread-1];
			_cy_r45[ithread] = _cy_r45[ithread-1];
			_cy_r46[ithread] = _cy_r46[ithread-1];
			_cy_r47[ithread] = _cy_r47[ithread-1];
			_cy_r48[ithread] = _cy_r48[ithread-1];
			_cy_r49[ithread] = _cy_r49[ithread-1];
			_cy_r50[ithread] = _cy_r50[ithread-1];
			_cy_r51[ithread] = _cy_r51[ithread-1];
			_cy_r52[ithread] = _cy_r52[ithread-1];
			_cy_r53[ithread] = _cy_r53[ithread-1];
			_cy_r54[ithread] = _cy_r54[ithread-1];
			_cy_r55[ithread] = _cy_r55[ithread-1];
			_cy_r56[ithread] = _cy_r56[ithread-1];
			_cy_r57[ithread] = _cy_r57[ithread-1];
			_cy_r58[ithread] = _cy_r58[ithread-1];
			_cy_r59[ithread] = _cy_r59[ithread-1];
		}

		_cy_r00[0] =+t3er;	/* ...The wraparound carry is here: */
		_cy_r01[0] = t00r;
		_cy_r02[0] = t01r;
		_cy_r03[0] = t02r;
		_cy_r04[0] = t03r;
		_cy_r05[0] = t04r;
		_cy_r06[0] = t05r;
		_cy_r07[0] = t06r;
		_cy_r08[0] = t07r;
		_cy_r09[0] = t08r;
		_cy_r10[0] = t09r;
		_cy_r11[0] = t0ar;
		_cy_r12[0] = t0br;
		_cy_r13[0] = t0cr;
		_cy_r14[0] = t0dr;
		_cy_r15[0] = t0er;
		_cy_r16[0] = t10r;
		_cy_r17[0] = t11r;
		_cy_r18[0] = t12r;
		_cy_r19[0] = t13r;
		_cy_r20[0] = t14r;
		_cy_r21[0] = t15r;
		_cy_r22[0] = t16r;
		_cy_r23[0] = t17r;
		_cy_r24[0] = t18r;
		_cy_r25[0] = t19r;
		_cy_r26[0] = t1ar;
		_cy_r27[0] = t1br;
		_cy_r28[0] = t1cr;
		_cy_r29[0] = t1dr;
		_cy_r30[0] = t1er;
		_cy_r31[0] = t20r;
		_cy_r32[0] = t21r;
		_cy_r33[0] = t22r;
		_cy_r34[0] = t23r;
		_cy_r35[0] = t24r;
		_cy_r36[0] = t25r;
		_cy_r37[0] = t26r;
		_cy_r38[0] = t27r;
		_cy_r39[0] = t28r;
		_cy_r40[0] = t29r;
		_cy_r41[0] = t2ar;
		_cy_r42[0] = t2br;
		_cy_r43[0] = t2cr;
		_cy_r44[0] = t2dr;
		_cy_r45[0] = t2er;
		_cy_r46[0] = t30r;
		_cy_r47[0] = t31r;
		_cy_r48[0] = t32r;
		_cy_r49[0] = t33r;
		_cy_r50[0] = t34r;
		_cy_r51[0] = t35r;
		_cy_r52[0] = t36r;
		_cy_r53[0] = t37r;
		_cy_r54[0] = t38r;
		_cy_r55[0] = t39r;
		_cy_r56[0] = t3ar;
		_cy_r57[0] = t3br;
		_cy_r58[0] = t3cr;
		_cy_r59[0] = t3dr;
	}
	else
	{
		t00r = _cy_r00[CY_THREADS - 1];	t00i = _cy_i00[CY_THREADS - 1];
		t01r = _cy_r01[CY_THREADS - 1];	t01i = _cy_i01[CY_THREADS - 1];
		t02r = _cy_r02[CY_THREADS - 1];	t02i = _cy_i02[CY_THREADS - 1];
		t03r = _cy_r03[CY_THREADS - 1];	t03i = _cy_i03[CY_THREADS - 1];
		t04r = _cy_r04[CY_THREADS - 1];	t04i = _cy_i04[CY_THREADS - 1];
		t05r = _cy_r05[CY_THREADS - 1];	t05i = _cy_i05[CY_THREADS - 1];
		t06r = _cy_r06[CY_THREADS - 1];	t06i = _cy_i06[CY_THREADS - 1];
		t07r = _cy_r07[CY_THREADS - 1];	t07i = _cy_i07[CY_THREADS - 1];
		t08r = _cy_r08[CY_THREADS - 1];	t08i = _cy_i08[CY_THREADS - 1];
		t09r = _cy_r09[CY_THREADS - 1];	t09i = _cy_i09[CY_THREADS - 1];
		t0ar = _cy_r10[CY_THREADS - 1];	t0ai = _cy_i10[CY_THREADS - 1];
		t0br = _cy_r11[CY_THREADS - 1];	t0bi = _cy_i11[CY_THREADS - 1];
		t0cr = _cy_r12[CY_THREADS - 1];	t0ci = _cy_i12[CY_THREADS - 1];
		t0dr = _cy_r13[CY_THREADS - 1];	t0di = _cy_i13[CY_THREADS - 1];
		t0er = _cy_r14[CY_THREADS - 1];	t0ei = _cy_i14[CY_THREADS - 1];
		t10r = _cy_r15[CY_THREADS - 1];	t10i = _cy_i15[CY_THREADS - 1];
		t11r = _cy_r16[CY_THREADS - 1];	t11i = _cy_i16[CY_THREADS - 1];
		t12r = _cy_r17[CY_THREADS - 1];	t12i = _cy_i17[CY_THREADS - 1];
		t13r = _cy_r18[CY_THREADS - 1];	t13i = _cy_i18[CY_THREADS - 1];
		t14r = _cy_r19[CY_THREADS - 1];	t14i = _cy_i19[CY_THREADS - 1];
		t15r = _cy_r20[CY_THREADS - 1];	t15i = _cy_i20[CY_THREADS - 1];
		t16r = _cy_r21[CY_THREADS - 1];	t16i = _cy_i21[CY_THREADS - 1];
		t17r = _cy_r22[CY_THREADS - 1];	t17i = _cy_i22[CY_THREADS - 1];
		t18r = _cy_r23[CY_THREADS - 1];	t18i = _cy_i23[CY_THREADS - 1];
		t19r = _cy_r24[CY_THREADS - 1];	t19i = _cy_i24[CY_THREADS - 1];
		t1ar = _cy_r25[CY_THREADS - 1];	t1ai = _cy_i25[CY_THREADS - 1];
		t1br = _cy_r26[CY_THREADS - 1];	t1bi = _cy_i26[CY_THREADS - 1];
		t1cr = _cy_r27[CY_THREADS - 1];	t1ci = _cy_i27[CY_THREADS - 1];
		t1dr = _cy_r28[CY_THREADS - 1];	t1di = _cy_i28[CY_THREADS - 1];
		t1er = _cy_r29[CY_THREADS - 1];	t1ei = _cy_i29[CY_THREADS - 1];
		t20r = _cy_r30[CY_THREADS - 1];	t20i = _cy_i30[CY_THREADS - 1];
		t21r = _cy_r31[CY_THREADS - 1];	t21i = _cy_i31[CY_THREADS - 1];
		t22r = _cy_r32[CY_THREADS - 1];	t22i = _cy_i32[CY_THREADS - 1];
		t23r = _cy_r33[CY_THREADS - 1];	t23i = _cy_i33[CY_THREADS - 1];
		t24r = _cy_r34[CY_THREADS - 1];	t24i = _cy_i34[CY_THREADS - 1];
		t25r = _cy_r35[CY_THREADS - 1];	t25i = _cy_i35[CY_THREADS - 1];
		t26r = _cy_r36[CY_THREADS - 1];	t26i = _cy_i36[CY_THREADS - 1];
		t27r = _cy_r37[CY_THREADS - 1];	t27i = _cy_i37[CY_THREADS - 1];
		t28r = _cy_r38[CY_THREADS - 1];	t28i = _cy_i38[CY_THREADS - 1];
		t29r = _cy_r39[CY_THREADS - 1];	t29i = _cy_i39[CY_THREADS - 1];
		t2ar = _cy_r40[CY_THREADS - 1];	t2ai = _cy_i40[CY_THREADS - 1];
		t2br = _cy_r41[CY_THREADS - 1];	t2bi = _cy_i41[CY_THREADS - 1];
		t2cr = _cy_r42[CY_THREADS - 1];	t2ci = _cy_i42[CY_THREADS - 1];
		t2dr = _cy_r43[CY_THREADS - 1];	t2di = _cy_i43[CY_THREADS - 1];
		t2er = _cy_r44[CY_THREADS - 1];	t2ei = _cy_i44[CY_THREADS - 1];
		t30r = _cy_r45[CY_THREADS - 1];	t30i = _cy_i45[CY_THREADS - 1];
		t31r = _cy_r46[CY_THREADS - 1];	t31i = _cy_i46[CY_THREADS - 1];
		t32r = _cy_r47[CY_THREADS - 1];	t32i = _cy_i47[CY_THREADS - 1];
		t33r = _cy_r48[CY_THREADS - 1];	t33i = _cy_i48[CY_THREADS - 1];
		t34r = _cy_r49[CY_THREADS - 1];	t34i = _cy_i49[CY_THREADS - 1];
		t35r = _cy_r50[CY_THREADS - 1];	t35i = _cy_i50[CY_THREADS - 1];
		t36r = _cy_r51[CY_THREADS - 1];	t36i = _cy_i51[CY_THREADS - 1];
		t37r = _cy_r52[CY_THREADS - 1];	t37i = _cy_i52[CY_THREADS - 1];
		t38r = _cy_r53[CY_THREADS - 1];	t38i = _cy_i53[CY_THREADS - 1];
		t39r = _cy_r54[CY_THREADS - 1];	t39i = _cy_i54[CY_THREADS - 1];
		t3ar = _cy_r55[CY_THREADS - 1];	t3ai = _cy_i55[CY_THREADS - 1];
		t3br = _cy_r56[CY_THREADS - 1];	t3bi = _cy_i56[CY_THREADS - 1];
		t3cr = _cy_r57[CY_THREADS - 1];	t3ci = _cy_i57[CY_THREADS - 1];
		t3dr = _cy_r58[CY_THREADS - 1];	t3di = _cy_i58[CY_THREADS - 1];
		t3er = _cy_r59[CY_THREADS - 1];	t3ei = _cy_i59[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1, "CY_THREADS must be > 1!");	/* Make sure loop only gets executed if multiple threads */
			_cy_r00[ithread] = _cy_r00[ithread-1];		_cy_i00[ithread] = _cy_i00[ithread-1];
			_cy_r01[ithread] = _cy_r01[ithread-1];		_cy_i01[ithread] = _cy_i01[ithread-1];
			_cy_r02[ithread] = _cy_r02[ithread-1];		_cy_i02[ithread] = _cy_i02[ithread-1];
			_cy_r03[ithread] = _cy_r03[ithread-1];		_cy_i03[ithread] = _cy_i03[ithread-1];
			_cy_r04[ithread] = _cy_r04[ithread-1];		_cy_i04[ithread] = _cy_i04[ithread-1];
			_cy_r05[ithread] = _cy_r05[ithread-1];		_cy_i05[ithread] = _cy_i05[ithread-1];
			_cy_r06[ithread] = _cy_r06[ithread-1];		_cy_i06[ithread] = _cy_i06[ithread-1];
			_cy_r07[ithread] = _cy_r07[ithread-1];		_cy_i07[ithread] = _cy_i07[ithread-1];
			_cy_r08[ithread] = _cy_r08[ithread-1];		_cy_i08[ithread] = _cy_i08[ithread-1];
			_cy_r09[ithread] = _cy_r09[ithread-1];		_cy_i09[ithread] = _cy_i09[ithread-1];
			_cy_r10[ithread] = _cy_r10[ithread-1];		_cy_i10[ithread] = _cy_i10[ithread-1];
			_cy_r11[ithread] = _cy_r11[ithread-1];		_cy_i11[ithread] = _cy_i11[ithread-1];
			_cy_r12[ithread] = _cy_r12[ithread-1];		_cy_i12[ithread] = _cy_i12[ithread-1];
			_cy_r13[ithread] = _cy_r13[ithread-1];		_cy_i13[ithread] = _cy_i13[ithread-1];
			_cy_r14[ithread] = _cy_r14[ithread-1];		_cy_i14[ithread] = _cy_i14[ithread-1];
			_cy_r15[ithread] = _cy_r15[ithread-1];		_cy_i15[ithread] = _cy_i15[ithread-1];
			_cy_r16[ithread] = _cy_r16[ithread-1];		_cy_i16[ithread] = _cy_i16[ithread-1];
			_cy_r17[ithread] = _cy_r17[ithread-1];		_cy_i17[ithread] = _cy_i17[ithread-1];
			_cy_r18[ithread] = _cy_r18[ithread-1];		_cy_i18[ithread] = _cy_i18[ithread-1];
			_cy_r19[ithread] = _cy_r19[ithread-1];		_cy_i19[ithread] = _cy_i19[ithread-1];
			_cy_r20[ithread] = _cy_r20[ithread-1];		_cy_i20[ithread] = _cy_i20[ithread-1];
			_cy_r21[ithread] = _cy_r21[ithread-1];		_cy_i21[ithread] = _cy_i21[ithread-1];
			_cy_r22[ithread] = _cy_r22[ithread-1];		_cy_i22[ithread] = _cy_i22[ithread-1];
			_cy_r23[ithread] = _cy_r23[ithread-1];		_cy_i23[ithread] = _cy_i23[ithread-1];
			_cy_r24[ithread] = _cy_r24[ithread-1];		_cy_i24[ithread] = _cy_i24[ithread-1];
			_cy_r25[ithread] = _cy_r25[ithread-1];		_cy_i25[ithread] = _cy_i25[ithread-1];
			_cy_r26[ithread] = _cy_r26[ithread-1];		_cy_i26[ithread] = _cy_i26[ithread-1];
			_cy_r27[ithread] = _cy_r27[ithread-1];		_cy_i27[ithread] = _cy_i27[ithread-1];
			_cy_r28[ithread] = _cy_r28[ithread-1];		_cy_i28[ithread] = _cy_i28[ithread-1];
			_cy_r29[ithread] = _cy_r29[ithread-1];		_cy_i29[ithread] = _cy_i29[ithread-1];
			_cy_r30[ithread] = _cy_r30[ithread-1];		_cy_i30[ithread] = _cy_i30[ithread-1];
			_cy_r31[ithread] = _cy_r31[ithread-1];		_cy_i31[ithread] = _cy_i31[ithread-1];
			_cy_r32[ithread] = _cy_r32[ithread-1];		_cy_i32[ithread] = _cy_i32[ithread-1];
			_cy_r33[ithread] = _cy_r33[ithread-1];		_cy_i33[ithread] = _cy_i33[ithread-1];
			_cy_r34[ithread] = _cy_r34[ithread-1];		_cy_i34[ithread] = _cy_i34[ithread-1];
			_cy_r35[ithread] = _cy_r35[ithread-1];		_cy_i35[ithread] = _cy_i35[ithread-1];
			_cy_r36[ithread] = _cy_r36[ithread-1];		_cy_i36[ithread] = _cy_i36[ithread-1];
			_cy_r37[ithread] = _cy_r37[ithread-1];		_cy_i37[ithread] = _cy_i37[ithread-1];
			_cy_r38[ithread] = _cy_r38[ithread-1];		_cy_i38[ithread] = _cy_i38[ithread-1];
			_cy_r39[ithread] = _cy_r39[ithread-1];		_cy_i39[ithread] = _cy_i39[ithread-1];
			_cy_r40[ithread] = _cy_r40[ithread-1];		_cy_i40[ithread] = _cy_i40[ithread-1];
			_cy_r41[ithread] = _cy_r41[ithread-1];		_cy_i41[ithread] = _cy_i41[ithread-1];
			_cy_r42[ithread] = _cy_r42[ithread-1];		_cy_i42[ithread] = _cy_i42[ithread-1];
			_cy_r43[ithread] = _cy_r43[ithread-1];		_cy_i43[ithread] = _cy_i43[ithread-1];
			_cy_r44[ithread] = _cy_r44[ithread-1];		_cy_i44[ithread] = _cy_i44[ithread-1];
			_cy_r45[ithread] = _cy_r45[ithread-1];		_cy_i45[ithread] = _cy_i45[ithread-1];
			_cy_r46[ithread] = _cy_r46[ithread-1];		_cy_i46[ithread] = _cy_i46[ithread-1];
			_cy_r47[ithread] = _cy_r47[ithread-1];		_cy_i47[ithread] = _cy_i47[ithread-1];
			_cy_r48[ithread] = _cy_r48[ithread-1];		_cy_i48[ithread] = _cy_i48[ithread-1];
			_cy_r49[ithread] = _cy_r49[ithread-1];		_cy_i49[ithread] = _cy_i49[ithread-1];
			_cy_r50[ithread] = _cy_r50[ithread-1];		_cy_i50[ithread] = _cy_i50[ithread-1];
			_cy_r51[ithread] = _cy_r51[ithread-1];		_cy_i51[ithread] = _cy_i51[ithread-1];
			_cy_r52[ithread] = _cy_r52[ithread-1];		_cy_i52[ithread] = _cy_i52[ithread-1];
			_cy_r53[ithread] = _cy_r53[ithread-1];		_cy_i53[ithread] = _cy_i53[ithread-1];
			_cy_r54[ithread] = _cy_r54[ithread-1];		_cy_i54[ithread] = _cy_i54[ithread-1];
			_cy_r55[ithread] = _cy_r55[ithread-1];		_cy_i55[ithread] = _cy_i55[ithread-1];
			_cy_r56[ithread] = _cy_r56[ithread-1];		_cy_i56[ithread] = _cy_i56[ithread-1];
			_cy_r57[ithread] = _cy_r57[ithread-1];		_cy_i57[ithread] = _cy_i57[ithread-1];
			_cy_r58[ithread] = _cy_r58[ithread-1];		_cy_i58[ithread] = _cy_i58[ithread-1];
			_cy_r59[ithread] = _cy_r59[ithread-1];		_cy_i59[ithread] = _cy_i59[ithread-1];
		}

		_cy_r00[0] =-t3ei;	_cy_i00[0] =+t3er;	/* ...The 2 Mo"bius carries are here: */
		_cy_r01[0] = t00r;	_cy_i01[0] = t00i;
		_cy_r02[0] = t01r;	_cy_i02[0] = t01i;
		_cy_r03[0] = t02r;	_cy_i03[0] = t02i;
		_cy_r04[0] = t03r;	_cy_i04[0] = t03i;
		_cy_r05[0] = t04r;	_cy_i05[0] = t04i;
		_cy_r06[0] = t05r;	_cy_i06[0] = t05i;
		_cy_r07[0] = t06r;	_cy_i07[0] = t06i;
		_cy_r08[0] = t07r;	_cy_i08[0] = t07i;
		_cy_r09[0] = t08r;	_cy_i09[0] = t08i;
		_cy_r10[0] = t09r;	_cy_i10[0] = t09i;
		_cy_r11[0] = t0ar;	_cy_i11[0] = t0ai;
		_cy_r12[0] = t0br;	_cy_i12[0] = t0bi;
		_cy_r13[0] = t0cr;	_cy_i13[0] = t0ci;
		_cy_r14[0] = t0dr;	_cy_i14[0] = t0di;
		_cy_r15[0] = t0er;	_cy_i15[0] = t0ei;
		_cy_r16[0] = t10r;	_cy_i16[0] = t10i;
		_cy_r17[0] = t11r;	_cy_i17[0] = t11i;
		_cy_r18[0] = t12r;	_cy_i18[0] = t12i;
		_cy_r19[0] = t13r;	_cy_i19[0] = t13i;
		_cy_r20[0] = t14r;	_cy_i20[0] = t14i;
		_cy_r21[0] = t15r;	_cy_i21[0] = t15i;
		_cy_r22[0] = t16r;	_cy_i22[0] = t16i;
		_cy_r23[0] = t17r;	_cy_i23[0] = t17i;
		_cy_r24[0] = t18r;	_cy_i24[0] = t18i;
		_cy_r25[0] = t19r;	_cy_i25[0] = t19i;
		_cy_r26[0] = t1ar;	_cy_i26[0] = t1ai;
		_cy_r27[0] = t1br;	_cy_i27[0] = t1bi;
		_cy_r28[0] = t1cr;	_cy_i28[0] = t1ci;
		_cy_r29[0] = t1dr;	_cy_i29[0] = t1di;
		_cy_r30[0] = t1er;	_cy_i30[0] = t1ei;
		_cy_r31[0] = t20r;	_cy_i31[0] = t20i;
		_cy_r32[0] = t21r;	_cy_i32[0] = t21i;
		_cy_r33[0] = t22r;	_cy_i33[0] = t22i;
		_cy_r34[0] = t23r;	_cy_i34[0] = t23i;
		_cy_r35[0] = t24r;	_cy_i35[0] = t24i;
		_cy_r36[0] = t25r;	_cy_i36[0] = t25i;
		_cy_r37[0] = t26r;	_cy_i37[0] = t26i;
		_cy_r38[0] = t27r;	_cy_i38[0] = t27i;
		_cy_r39[0] = t28r;	_cy_i39[0] = t28i;
		_cy_r40[0] = t29r;	_cy_i40[0] = t29i;
		_cy_r41[0] = t2ar;	_cy_i41[0] = t2ai;
		_cy_r42[0] = t2br;	_cy_i42[0] = t2bi;
		_cy_r43[0] = t2cr;	_cy_i43[0] = t2ci;
		_cy_r44[0] = t2dr;	_cy_i44[0] = t2di;
		_cy_r45[0] = t2er;	_cy_i45[0] = t2ei;
		_cy_r46[0] = t30r;	_cy_i46[0] = t30i;
		_cy_r47[0] = t31r;	_cy_i47[0] = t31i;
		_cy_r48[0] = t32r;	_cy_i48[0] = t32i;
		_cy_r49[0] = t33r;	_cy_i49[0] = t33i;
		_cy_r50[0] = t34r;	_cy_i50[0] = t34i;
		_cy_r51[0] = t35r;	_cy_i51[0] = t35i;
		_cy_r52[0] = t36r;	_cy_i52[0] = t36i;
		_cy_r53[0] = t37r;	_cy_i53[0] = t37i;
		_cy_r54[0] = t38r;	_cy_i54[0] = t38i;
		_cy_r55[0] = t39r;	_cy_i55[0] = t39i;
		_cy_r56[0] = t3ar;	_cy_i56[0] = t3ai;
		_cy_r57[0] = t3br;	_cy_i57[0] = t3bi;
		_cy_r58[0] = t3cr;	_cy_i58[0] = t3ci;
		_cy_r59[0] = t3dr;	_cy_i59[0] = t3di;
	}

	full_pass = 0;
	scale = 1;

	/*
	For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
	*/
	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		j_jhi =15;
	}
	else
	{
		j_jhi = 7;
	}

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		for(j = ithread*pini; j <= ithread*pini + j_jhi; j++)
		{
			jt = j;	// NB: pini *already* padded, so do not additionally pad index here!
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p04;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p08;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p12;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p16;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p20;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p24;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p28;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p32;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p36;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p40;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p44;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p48;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p52;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
			jt = j + p56;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
		}
	}
}	/* endfor(outer) */

	t00r = 0;
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		t00r += fabs(_cy_r00[0])+fabs(_cy_i00[0])+fabs(_cy_r01[0])+fabs(_cy_i01[0])+fabs(_cy_r02[0])+fabs(_cy_i02[0])+fabs(_cy_r03[0])+fabs(_cy_i03[0])+fabs(_cy_r04[0])+fabs(_cy_i04[0])+fabs(_cy_r05[0])+fabs(_cy_i05[0])+fabs(_cy_r06[0])+fabs(_cy_i06[0])+fabs(_cy_r07[0])+fabs(_cy_i07[0])+fabs(_cy_r08[0])+fabs(_cy_i08[0])+fabs(_cy_r09[0])+fabs(_cy_i09[0]);
		t00r += fabs(_cy_r10[0])+fabs(_cy_i10[0])+fabs(_cy_r11[0])+fabs(_cy_i11[0])+fabs(_cy_r12[0])+fabs(_cy_i12[0])+fabs(_cy_r13[0])+fabs(_cy_i13[0])+fabs(_cy_r14[0])+fabs(_cy_i14[0])+fabs(_cy_r15[0])+fabs(_cy_i15[0])+fabs(_cy_r16[0])+fabs(_cy_i16[0])+fabs(_cy_r17[0])+fabs(_cy_i17[0])+fabs(_cy_r18[0])+fabs(_cy_i18[0])+fabs(_cy_r19[0])+fabs(_cy_i19[0]);
		t00r += fabs(_cy_r20[0])+fabs(_cy_i20[0])+fabs(_cy_r21[0])+fabs(_cy_i21[0])+fabs(_cy_r22[0])+fabs(_cy_i22[0])+fabs(_cy_r23[0])+fabs(_cy_i23[0])+fabs(_cy_r24[0])+fabs(_cy_i24[0])+fabs(_cy_r25[0])+fabs(_cy_i25[0])+fabs(_cy_r26[0])+fabs(_cy_i26[0])+fabs(_cy_r27[0])+fabs(_cy_i27[0])+fabs(_cy_r28[0])+fabs(_cy_i28[0])+fabs(_cy_r29[0])+fabs(_cy_i29[0]);
		t00r += fabs(_cy_r30[0])+fabs(_cy_i30[0])+fabs(_cy_r31[0])+fabs(_cy_i31[0])+fabs(_cy_r32[0])+fabs(_cy_i32[0])+fabs(_cy_r33[0])+fabs(_cy_i33[0])+fabs(_cy_r34[0])+fabs(_cy_i34[0])+fabs(_cy_r35[0])+fabs(_cy_i35[0])+fabs(_cy_r36[0])+fabs(_cy_i36[0])+fabs(_cy_r37[0])+fabs(_cy_i37[0])+fabs(_cy_r38[0])+fabs(_cy_i38[0])+fabs(_cy_r39[0])+fabs(_cy_i39[0]);
		t00r += fabs(_cy_r40[0])+fabs(_cy_i40[0])+fabs(_cy_r41[0])+fabs(_cy_i41[0])+fabs(_cy_r42[0])+fabs(_cy_i42[0])+fabs(_cy_r43[0])+fabs(_cy_i43[0])+fabs(_cy_r44[0])+fabs(_cy_i44[0])+fabs(_cy_r45[0])+fabs(_cy_i45[0])+fabs(_cy_r46[0])+fabs(_cy_i46[0])+fabs(_cy_r47[0])+fabs(_cy_i47[0])+fabs(_cy_r48[0])+fabs(_cy_i48[0])+fabs(_cy_r49[0])+fabs(_cy_i49[0]);
		t00r += fabs(_cy_r50[0])+fabs(_cy_i50[0])+fabs(_cy_r51[0])+fabs(_cy_i51[0])+fabs(_cy_r52[0])+fabs(_cy_i52[0])+fabs(_cy_r53[0])+fabs(_cy_i53[0])+fabs(_cy_r54[0])+fabs(_cy_i54[0])+fabs(_cy_r55[0])+fabs(_cy_i55[0])+fabs(_cy_r56[0])+fabs(_cy_i56[0])+fabs(_cy_r57[0])+fabs(_cy_i57[0])+fabs(_cy_r58[0])+fabs(_cy_i58[0])+fabs(_cy_r59[0])+fabs(_cy_i59[0]);

		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
	}

	if(t00r != 0.0)
	{
		sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in %s - input wordsize may be too small.\n",iter,func);
		if(INTERACT)fprintf(stderr,"%s",cbuf);
		fp = fopen(   OFILE,"a");
		fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;
		err=ERR_CARRY;
		return(err);
	}
	return(0);
}


/***************/

void radix60_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-60 complex DIF FFT pass on the data in the length-N real vector A.
*/
	int j,jt,jp,j1,j2;
	static int n60,first_entry=TRUE,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40,p44,p48,p52,p56;
	static double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/

					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it,
		x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,
		t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,
		t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,
		t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,
		t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei;

	if(!first_entry && (n/60) != n60)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n60=n/60;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n60;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;
		p28 = p24 + p04;
		p32 = p28 + p04;
		p36 = p32 + p04;
		p40 = p36 + p04;
		p44 = p40 + p04;
		p48 = p44 + p04;
		p52 = p48 + p04;
		p56 = p52 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p36 = p36 + ( (p36 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		p44 = p44 + ( (p44 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p52 = p52 + ( (p52 >> DAT_BITS) << PAD_BITS );
		p56 = p56 + ( (p56 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-60 pass is here.	*/

	for(j = 0; j < n60; j += 2)
	{
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;

		/*...gather the needed data (60 64-bit complex, i.e.120 64-bit reals) and do 4 radix-15 transforms...*/
	/*
	Twiddleless version arranges 4 sets of radix-15 DFT inputs as follows:
	0 in upper left corner, decrement 4 horizontally and 15 vertically, indexing modulo 60:

		RADIX_15_DFT(00,56,52,48,44,40,36,32,28,24,20,16,12,08,04)
		RADIX_15_DFT(45,41,37,33,29,25,21,17,13,09,05,01,57,53,49)
		RADIX_15_DFT(30,26,22,18,14,10,06,02,58,54,50,46,42,38,34)
		RADIX_15_DFT(15,11,07,03,59,55,51,47,43,39,35,31,27,23,19)

	If we subtract costant offset 1,2,3 from the last 3 rows as we do in the implementation to reduce the number
	of index offsets needing to be stored, we decrement 4 horizontally and 16 vertically.
	*/
	#if 0
		jt = j1    ; jp = j2    ;	RADIX_15_DIF(a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,rt,it);
		jt = j1+p01; jp = j2+p01;	RADIX_15_DIF(a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,rt,it);
		jt = j1+p02; jp = j2+p02;	RADIX_15_DIF(a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,rt,it);
		jt = j1+p03; jp = j2+p03;	RADIX_15_DIF(a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei,rt,it);
	#else
		jt = j1    ; jp = j2    ;	RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,rt,it);
		jt = j1+p01; jp = j2+p01;	RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,rt,it);
		jt = j1+p02; jp = j2+p02;	RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,rt,it);
		jt = j1+p03; jp = j2+p03;	RADIX_15_DIF_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p56],a[jp+p56],a[jt+p52],a[jp+p52],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei,rt,it);
	#endif
		/*...and now do 15 radix-4 transforms.
		The required output permutation is:

			[ 0, 1, 2, 3
			, 9, 8,11,10
			, 6, 7, 5, 4
			,57,56,59,58
			,54,55,53,52
			,48,49,50,51
			,46,47,45,44
			,40,41,42,43
			,39,38,36,37
			,32,33,34,35
			,31,30,28,29
			,25,24,27,26
			,23,22,20,21
			,17,16,19,18
			,14,15,13,12].
		*/
												/*          inputs                    */ /*                                      outputs                              */
		jt = j1    ; jp = j2    ;	RADIX_04_DIF(t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIF(t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIF(t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
		jt = j1+p56; jp = j2+p56;	RADIX_04_DIF(t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
		jt = j1+p52; jp = j2+p52;	RADIX_04_DIF(t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
		jt = j1+p48; jp = j2+p48;	RADIX_04_DIF(t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
		jt = j1+p44; jp = j2+p44;	RADIX_04_DIF(t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
		jt = j1+p40; jp = j2+p40;	RADIX_04_DIF(t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
		jt = j1+p36; jp = j2+p36;	RADIX_04_DIF(t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIF(t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIF(t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIF(t0br,t0bi,t1br,t1bi,t2br,t2bi,t3br,t3bi,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIF(t0cr,t0ci,t1cr,t1ci,t2cr,t2ci,t3cr,t3ci,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIF(t0dr,t0di,t1dr,t1di,t2dr,t2di,t3dr,t3di,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIF(t0er,t0ei,t1er,t1ei,t2er,t2ei,t3er,t3ei,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
		/* Totals: 4*radix15 + 15*radix04 = 4*(162 FADD, 50 FMUL) + 15*(16 FADD, 0 FMUL) = 888 FADD, 200 FMUL	*/
	}
}

/***************/

void radix60_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Time
!
!...Subroutine to perform an initial radix-60 complex DIT FFT pass on the data in the length-N real vector A.
*/
	int j,jt,jp,j1,j2;
	static int n60,first_entry=TRUE,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40,p44,p48,p52,p56;
	static double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/

					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it,
		x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,
		t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,
		t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,
		t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,
		t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei;

	if(!first_entry && (n/60) != n60)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n60=n/60;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n60;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;
		p28 = p24 + p04;
		p32 = p28 + p04;
		p36 = p32 + p04;
		p40 = p36 + p04;
		p44 = p40 + p04;
		p48 = p44 + p04;
		p52 = p48 + p04;
		p56 = p52 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p36 = p36 + ( (p36 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		p44 = p44 + ( (p44 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p52 = p52 + ( (p52 >> DAT_BITS) << PAD_BITS );
		p56 = p56 + ( (p56 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-60 pass is here.	*/

	for(j = 0; j < n60; j += 2)
	{
	#ifdef USE_AVX
		j1 = (j & mask02) + br8[j&7];
	#elif defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
	#else
		j1 = j;
	#endif
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;

	/*
		Twiddleless version requires us to swap inputs as follows:
		indices x[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59
			  ->x[0 56,52,48,44,40,36,32,28,24,20,16,12,08,04|45,41,37,33,29,25,21,17,13,09,05,01,57,53,49|30,26,22,18,14,10,06,02,58,54,50,46,42,38,34|15,11,07,03,59,55,51,47,43,39,35,31,27,23,19

		I.e. start out with first quartet of indices {0,15,30,45}, permute those according to
		(0,15,30,45}*59%60 = {0,45,30,15), then each is head of a length-15 list of indices with decrement 4.

		Remember, inputs to DIT are bit-reversed, so
		a[ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59] contain
		x[ 0 30 15 45 05 35 20 50 10 40 25 55 01 31 16 46 06 36 21 51 11 41 26 56 02 32 17 47 07 37 22 52 12 42 27 57 03 33 18 48 08 38 23 53 13 43 28 58 04 34 19 49 09 39 24 54 14 44 29 59], which get swapped to
	(Look at first nontrivial one, i.e. x[1 -> 56] ... in terms of a[] this translates to a[12 -> 23])
	-->	x[ 0 30 45 15 40 10 25 55 20 50 05 35 56 26 41 11 36 06 21 51 16 46 01 31 52 22 37 07 32 02 17 47 12 42 57 27 48 18 33 03 28 58 13 43 08 38 53 23 44 14 29 59 24 54 09 39 04 34 49 19]
		a[ 0  1  3  2| 9  8 10 11| 6  7  4  5|23 22 21 20|17 16 18 19|14 15 12 13|31 30 29 28|25 24 26 27|32 33 35 34|39 38 37 36|46 47 44 45|40 41 43 42|57 56 58 59|54 55 52 53|48 49 51 50].
	*/
		/*...gather the needed data (60 64-bit complex, i.e.120 64-bit reals) and do 15 radix-4 transforms...*/
												/*                                   inputs                                   */ /*              outputs              */
		jt = j1    ; jp = j2    ;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,rt,it);
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,rt,it);
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,rt,it);
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,rt,it);
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,rt,it);
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,rt,it);
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,rt,it);
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,rt,it);
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,rt,it);
		jt = j1+p36; jp = j2+p36;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,rt,it);
		jt = j1+p44; jp = j2+p44;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,rt,it);
		jt = j1+p40; jp = j2+p40;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0br,t0bi,t1br,t1bi,t2br,t2bi,t3br,t3bi,rt,it);
		jt = j1+p56; jp = j2+p56;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t0cr,t0ci,t1cr,t1ci,t2cr,t2ci,t3cr,t3ci,rt,it);
		jt = j1+p52; jp = j2+p52;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t0dr,t0di,t1dr,t1di,t2dr,t2di,t3dr,t3di,rt,it);
		jt = j1+p48; jp = j2+p48;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0er,t0ei,t1er,t1ei,t2er,t2ei,t3er,t3ei,rt,it);

		/*...and now do 4 radix-15 transforms.
		The required output permutation is:

			[ 0,44,28,12,56,40,24, 8,52,36,20, 4,48,32,16
			 15,59,43,27,11,55,39,23, 7,51,35,19, 3,47,31
			 30,14,58,42,26,10,54,38,22, 6,50,34,18, 2,46
			 45,29,13,57,41,25, 9,53,37,21, 5,49,33,17, 1]
		*/
												/*                                                  inputs                                                                                          */ /*                                                                    temps                                                                        */ /* outputs: --> */
	#if 0
		jt = j1    ; jp = j2    ;	RADIX_15_DIT_C(t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a[jt    ],a[jp    ],a[jt+p44],a[jp+p44],a[jt+p28],a[jp+p28],a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],rt,it);
		jt = j1+p03; jp = j2+p03;	RADIX_15_DIT_C(t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p44],a[jp+p44],a[jt+p28],a[jp+p28],rt,it);
		jt = j1+p02; jp = j2+p02;	RADIX_15_DIT_C(t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a[jt+p28],a[jp+p28],a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p44],a[jp+p44],rt,it);
		jt = j1+p01; jp = j2+p01;	RADIX_15_DIT_C(t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a[jt+p44],a[jp+p44],a[jt+p28],a[jp+p28],a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],rt,it);
	#else
		jt = j1    ; jp = j2    ;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a[jt    ],a[jp    ],a[jt+p44],a[jp+p44],a[jt+p28],a[jp+p28],a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],rt,it);
		jt = j1+p03; jp = j2+p03;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p44],a[jp+p44],a[jt+p28],a[jp+p28],rt,it);
		jt = j1+p02; jp = j2+p02;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a[jt+p28],a[jp+p28],a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],a[jt+p44],a[jp+p44],rt,it);
		jt = j1+p01; jp = j2+p01;	RADIX_15_DIT_B(s,c3m1,cn1,cn2,ss3,sn1,sn2,t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a[jt+p44],a[jp+p44],a[jt+p28],a[jp+p28],a[jt+p12],a[jp+p12],a[jt+p56],a[jp+p56],a[jt+p40],a[jp+p40],a[jt+p24],a[jp+p24],a[jt+p08],a[jp+p08],a[jt+p52],a[jp+p52],a[jt+p36],a[jp+p36],a[jt+p20],a[jp+p20],a[jt+p04],a[jp+p04],a[jt+p48],a[jp+p48],a[jt+p32],a[jp+p32],a[jt+p16],a[jp+p16],a[jt    ],a[jp    ],rt,it);
	#endif
	}
}

/******************** Multithreaded function body: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void* 
	cy60_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		const uint32 RADIX = 60, odd_radix = 15;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40,p44,p48,p52,p56;
		int j,j1,j2,k,l;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
	#ifdef USE_AVX
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#else
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	#endif
		double rt,it,wt_re,wt_im;	/* Fermat-mod weights stuff, used in both scalar and AVX mode */
		int k1,k2;

	#ifdef USE_SSE2

		const double c3m1= -1.50000000000000000000;	/* cos(twopi/3)-1	*/
	  #ifdef USE_AVX
		const int l2_sz_vd = 5;
	  #else
		const int l2_sz_vd = 4;
	  #endif
		double *add0, *add1, *add2, *add3;
		const double crnd = 3.0*0x4000000*0x2000000;
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,
			*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,
			*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27,*bjmodn28,*bjmodn29,
			*bjmodn30,*bjmodn31,*bjmodn32,*bjmodn33,*bjmodn34,*bjmodn35,*bjmodn36,*bjmodn37,*bjmodn38,*bjmodn39,
			*bjmodn40,*bjmodn41,*bjmodn42,*bjmodn43,*bjmodn44,*bjmodn45,*bjmodn46,*bjmodn47,*bjmodn48,*bjmodn49,
			*bjmodn50,*bjmodn51,*bjmodn52,*bjmodn53,*bjmodn54,*bjmodn55,*bjmodn56,*bjmodn57,*bjmodn58,*bjmodn59;
		vec_dbl *tmp, *tm2, *max_err, *sse2_rnd, *half_arr, *sse2_c3m1, *sse2_s, *sse2_cn1, *sse2_cn2, *sse2_ss3, *sse2_sn1, *sse2_sn2,
			*s1p00r,*s1p00i,*s1p01r,*s1p01i,*s1p02r,*s1p02i,*s1p03r,*s1p03i,*s1p04r,*s1p04i,*s1p05r,*s1p05i,*s1p06r,*s1p06i,*s1p07r,*s1p07i,*s1p08r,*s1p08i,*s1p09r,*s1p09i,*s1p0ar,*s1p0ai,*s1p0br,*s1p0bi,*s1p0cr,*s1p0ci,*s1p0dr,*s1p0di,*s1p0er,*s1p0ei,
			*s1p10r,*s1p10i,*s1p11r,*s1p11i,*s1p12r,*s1p12i,*s1p13r,*s1p13i,*s1p14r,*s1p14i,*s1p15r,*s1p15i,*s1p16r,*s1p16i,*s1p17r,*s1p17i,*s1p18r,*s1p18i,*s1p19r,*s1p19i,*s1p1ar,*s1p1ai,*s1p1br,*s1p1bi,*s1p1cr,*s1p1ci,*s1p1dr,*s1p1di,*s1p1er,*s1p1ei,
			*s1p20r,*s1p20i,*s1p21r,*s1p21i,*s1p22r,*s1p22i,*s1p23r,*s1p23i,*s1p24r,*s1p24i,*s1p25r,*s1p25i,*s1p26r,*s1p26i,*s1p27r,*s1p27i,*s1p28r,*s1p28i,*s1p29r,*s1p29i,*s1p2ar,*s1p2ai,*s1p2br,*s1p2bi,*s1p2cr,*s1p2ci,*s1p2dr,*s1p2di,*s1p2er,*s1p2ei,
			*s1p30r,*s1p30i,*s1p31r,*s1p31i,*s1p32r,*s1p32i,*s1p33r,*s1p33i,*s1p34r,*s1p34i,*s1p35r,*s1p35i,*s1p36r,*s1p36i,*s1p37r,*s1p37i,*s1p38r,*s1p38i,*s1p39r,*s1p39i,*s1p3ar,*s1p3ai,*s1p3br,*s1p3bi,*s1p3cr,*s1p3ci,*s1p3dr,*s1p3di,*s1p3er,*s1p3ei,
			*r00,*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09,*r0a,*r0b,*r0c,*r0d,*r0e,
			*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r1a,*r1b,*r1c,*r1d,*r1e,
			*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r2a,*r2b,*r2c,*r2d,*r2e,
			*r30,*r31,*r32,*r33,*r34,*r35,*r36,*r37,*r38,*r39,*r3a,*r3b,*r3c,*r3d,*r3e,
			*x00,*x01,*x02,*x03,*x04,*x05,*x06,*x07,*x08,*x09,*x0a,*x0b,*x0c,*x0d,*x0e,
			*y00,*y01,*y02,*y03,*y04,*y05,*y06,*y07,*y08,*y09,*y0a,*y0b,*y0c,*y0d,*y0e,
			*cy_r00,*cy_r04,*cy_r08,*cy_r12,*cy_r16,*cy_r20,*cy_r24,*cy_r28,*cy_r32,*cy_r36,*cy_r40,*cy_r44,*cy_r48,*cy_r52,*cy_r56,
			*cy_i00,*cy_i04,*cy_i08,*cy_i12,*cy_i16,*cy_i20,*cy_i24,*cy_i28,*cy_i32,*cy_i36,*cy_i40,*cy_i44,*cy_i48,*cy_i52,*cy_i56;
	  #ifndef USE_AVX
		vec_dbl
			*cy_r02,*cy_r06,*cy_r10,*cy_r14,*cy_r18,*cy_r22,*cy_r26,*cy_r30,*cy_r34,*cy_r38,*cy_r42,*cy_r46,*cy_r50,*cy_r54,*cy_r58,
			*cy_i02,*cy_i06,*cy_i10,*cy_i14,*cy_i18,*cy_i22,*cy_i26,*cy_i30,*cy_i34,*cy_i38,*cy_i42,*cy_i46,*cy_i50,*cy_i54,*cy_i58;
	  #else
			vec_dbl *base_negacyclic_root;
	  #endif

		/* These are used in conjunction with the langth-7 arrays in the USE_SCALAR_CARRY #define below;
		In SSE2 mode store doubled versions of these data in the scratch storage accessed via the half_arr pointer: */
		int idx_offset, idx_incr;
		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else

		const double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
						s   =  0.86602540378443864675,	/* sin(twopi/3)		*/
						cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
						cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
						ss3 =  0.95105651629515357210,	/*  sin(u) */
						sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
						sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
		double *base, *baseinv, *wt_arr, *wtinv_arr, *bs_arr, *bsinv_arr, bs,bsinv, wts_idx_incr;
		const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int jt,jp,m,m2,ntmp;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		double temp,frac,
			a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,a1p08r,a1p08i,a1p09r,a1p09i,a1p0ar,a1p0ai,a1p0br,a1p0bi,a1p0cr,a1p0ci,a1p0dr,a1p0di,a1p0er,a1p0ei,
			a1p10r,a1p10i,a1p11r,a1p11i,a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a1p1ar,a1p1ai,a1p1br,a1p1bi,a1p1cr,a1p1ci,a1p1dr,a1p1di,a1p1er,a1p1ei,
			a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,a1p24r,a1p24i,a1p25r,a1p25i,a1p26r,a1p26i,a1p27r,a1p27i,a1p28r,a1p28i,a1p29r,a1p29i,a1p2ar,a1p2ai,a1p2br,a1p2bi,a1p2cr,a1p2ci,a1p2dr,a1p2di,a1p2er,a1p2ei,
			a1p30r,a1p30i,a1p31r,a1p31i,a1p32r,a1p32i,a1p33r,a1p33i,a1p34r,a1p34i,a1p35r,a1p35i,a1p36r,a1p36i,a1p37r,a1p37i,a1p38r,a1p38i,a1p39r,a1p39i,a1p3ar,a1p3ai,a1p3br,a1p3bi,a1p3cr,a1p3ci,a1p3dr,a1p3di,a1p3er,a1p3ei,
			cy_r00,cy_i00,cy_r01,cy_i01,cy_r02,cy_i02,cy_r03,cy_i03,cy_r04,cy_i04,cy_r05,cy_i05,cy_r06,cy_i06,cy_r07,cy_i07,cy_r08,cy_i08,cy_r09,cy_i09,
			cy_r10,cy_i10,cy_r11,cy_i11,cy_r12,cy_i12,cy_r13,cy_i13,cy_r14,cy_i14,cy_r15,cy_i15,cy_r16,cy_i16,cy_r17,cy_i17,cy_r18,cy_i18,cy_r19,cy_i19,
			cy_r20,cy_i20,cy_r21,cy_i21,cy_r22,cy_i22,cy_r23,cy_i23,cy_r24,cy_i24,cy_r25,cy_i25,cy_r26,cy_i26,cy_r27,cy_i27,cy_r28,cy_i28,cy_r29,cy_i29,
			cy_r30,cy_i30,cy_r31,cy_i31,cy_r32,cy_i32,cy_r33,cy_i33,cy_r34,cy_i34,cy_r35,cy_i35,cy_r36,cy_i36,cy_r37,cy_i37,cy_r38,cy_i38,cy_r39,cy_i39,
			cy_r40,cy_i40,cy_r41,cy_i41,cy_r42,cy_i42,cy_r43,cy_i43,cy_r44,cy_i44,cy_r45,cy_i45,cy_r46,cy_i46,cy_r47,cy_i47,cy_r48,cy_i48,cy_r49,cy_i49,
			cy_r50,cy_i50,cy_r51,cy_i51,cy_r52,cy_i52,cy_r53,cy_i53,cy_r54,cy_i54,cy_r55,cy_i55,cy_r56,cy_i56,cy_r57,cy_i57,cy_r58,cy_i58,cy_r59,cy_i59,
			x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,
			t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,
			t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,
			t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,
			t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei;
		int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,
			bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,
			bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,bjmodn28,bjmodn29,
			bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,bjmodn36,bjmodn37,bjmodn38,bjmodn39,
			bjmodn40,bjmodn41,bjmodn42,bjmodn43,bjmodn44,bjmodn45,bjmodn46,bjmodn47,bjmodn48,bjmodn49,
			bjmodn50,bjmodn51,bjmodn52,bjmodn53,bjmodn54,bjmodn55,bjmodn56,bjmodn57,bjmodn58,bjmodn59;

	#endif

		struct cy_thread_data_t* thread_arg = targ;
	// int data:
		int NDIVR = thread_arg->ndivr;
		int n = NDIVR*RADIX;
		int khi    = thread_arg->khi;
		int i      = thread_arg->i;	/* Pointer to the BASE and BASEINV arrays.	*/
		int jstart = thread_arg->jstart;
		int jhi    = thread_arg->jhi;
		int col = thread_arg->col;
		int co2 = thread_arg->co2;
		int co3 = thread_arg->co3;
		int sw  = thread_arg->sw;
		int nwt = thread_arg->nwt;
		int icycle00 = thread_arg->icycle00;
		int icycle01 = thread_arg->icycle01;
		int icycle02 = thread_arg->icycle02;
		int icycle03 = thread_arg->icycle03;
		int icycle04 = thread_arg->icycle04;
		int icycle05 = thread_arg->icycle05;
		int icycle06 = thread_arg->icycle06;
		int icycle07 = thread_arg->icycle07;
		int icycle08 = thread_arg->icycle08;
		int icycle09 = thread_arg->icycle09;
		int icycle10 = thread_arg->icycle10;
		int icycle11 = thread_arg->icycle11;
		int icycle12 = thread_arg->icycle12;
		int icycle13 = thread_arg->icycle13;
		int icycle14 = thread_arg->icycle14;
	#ifdef USE_SSE2
		int wts_idx_inc2 = thread_arg->wts_idx_inc2;
		int jcycle00 = thread_arg->jcycle00;
		int jcycle01 = thread_arg->jcycle01;
		int jcycle02 = thread_arg->jcycle02;
		int jcycle03 = thread_arg->jcycle03;
		int jcycle04 = thread_arg->jcycle04;
		int jcycle05 = thread_arg->jcycle05;
		int jcycle06 = thread_arg->jcycle06;
		int jcycle07 = thread_arg->jcycle07;
		int jcycle08 = thread_arg->jcycle08;
		int jcycle09 = thread_arg->jcycle09;
		int jcycle10 = thread_arg->jcycle10;
		int jcycle11 = thread_arg->jcycle11;
		int jcycle12 = thread_arg->jcycle12;
		int jcycle13 = thread_arg->jcycle13;
		int jcycle14 = thread_arg->jcycle14;
	  #ifdef USE_AVX
		int kcycle00 = thread_arg->kcycle00;
		int kcycle01 = thread_arg->kcycle01;
		int kcycle02 = thread_arg->kcycle02;
		int kcycle03 = thread_arg->kcycle03;
		int kcycle04 = thread_arg->kcycle04;
		int kcycle05 = thread_arg->kcycle05;
		int kcycle06 = thread_arg->kcycle06;
		int kcycle07 = thread_arg->kcycle07;
		int kcycle08 = thread_arg->kcycle08;
		int kcycle09 = thread_arg->kcycle09;
		int kcycle10 = thread_arg->kcycle10;
		int kcycle11 = thread_arg->kcycle11;
		int kcycle12 = thread_arg->kcycle12;
		int kcycle13 = thread_arg->kcycle13;
		int kcycle14 = thread_arg->kcycle14;

		int lcycle00 = thread_arg->lcycle00;
		int lcycle01 = thread_arg->lcycle01;
		int lcycle02 = thread_arg->lcycle02;
		int lcycle03 = thread_arg->lcycle03;
		int lcycle04 = thread_arg->lcycle04;
		int lcycle05 = thread_arg->lcycle05;
		int lcycle06 = thread_arg->lcycle06;
		int lcycle07 = thread_arg->lcycle07;
		int lcycle08 = thread_arg->lcycle08;
		int lcycle09 = thread_arg->lcycle09;
		int lcycle10 = thread_arg->lcycle10;
		int lcycle11 = thread_arg->lcycle11;
		int lcycle12 = thread_arg->lcycle12;
		int lcycle13 = thread_arg->lcycle13;
		int lcycle14 = thread_arg->lcycle14;
	  #endif
	#endif

	// double data:
		double maxerr = thread_arg->maxerr;
		double scale = thread_arg->scale;

	// pointer data:
		double *a = thread_arg->arrdat;
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
		int *si = thread_arg->si;
		struct complex *rn0 = thread_arg->rn0;
		struct complex *rn1 = thread_arg->rn1;

		/*   constant index offsets for array load/stores are here.	*/
		p01 = NDIVR;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;
		p28 = p24 + p04;
		p32 = p28 + p04;
		p36 = p32 + p04;
		p40 = p36 + p04;
		p44 = p40 + p04;
		p48 = p44 + p04;
		p52 = p48 + p04;
		p56 = p52 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p32 = p32 + ( (p32 >> DAT_BITS) << PAD_BITS );
		p36 = p36 + ( (p36 >> DAT_BITS) << PAD_BITS );
		p40 = p40 + ( (p40 >> DAT_BITS) << PAD_BITS );
		p44 = p44 + ( (p44 >> DAT_BITS) << PAD_BITS );
		p48 = p48 + ( (p48 >> DAT_BITS) << PAD_BITS );
		p52 = p52 + ( (p52 >> DAT_BITS) << PAD_BITS );
		p56 = p56 + ( (p56 >> DAT_BITS) << PAD_BITS );

	#ifdef USE_SSE2
		uint32 nwt16 = nwt << l2_sz_vd;	// nwt*sizeof(vec_dbl); the '16' is a historical naming artifact dating to first SSE2 code

		tmp = thread_arg->s1p00r;
									// Use these as handy address-offset temps:
									tm2 = tmp + 0x1e;    max_err  = tm2 + 0x1e;    sse2_rnd = max_err + 0x1e;
		s1p00r = tmp + 0x00;    s1p10r = tm2 + 0x00;    s1p20r = max_err + 0x00;    s1p30r = sse2_rnd + 0x00;
		s1p00i = tmp + 0x01;    s1p10i = tm2 + 0x01;    s1p20i = max_err + 0x01;    s1p30i = sse2_rnd + 0x01;
		s1p01r = tmp + 0x02;    s1p11r = tm2 + 0x02;    s1p21r = max_err + 0x02;    s1p31r = sse2_rnd + 0x02;
		s1p01i = tmp + 0x03;    s1p11i = tm2 + 0x03;    s1p21i = max_err + 0x03;    s1p31i = sse2_rnd + 0x03;
		s1p02r = tmp + 0x04;    s1p12r = tm2 + 0x04;    s1p22r = max_err + 0x04;    s1p32r = sse2_rnd + 0x04;
		s1p02i = tmp + 0x05;    s1p12i = tm2 + 0x05;    s1p22i = max_err + 0x05;    s1p32i = sse2_rnd + 0x05;
		s1p03r = tmp + 0x06;    s1p13r = tm2 + 0x06;    s1p23r = max_err + 0x06;    s1p33r = sse2_rnd + 0x06;
		s1p03i = tmp + 0x07;    s1p13i = tm2 + 0x07;    s1p23i = max_err + 0x07;    s1p33i = sse2_rnd + 0x07;
		s1p04r = tmp + 0x08;    s1p14r = tm2 + 0x08;    s1p24r = max_err + 0x08;    s1p34r = sse2_rnd + 0x08;
		s1p04i = tmp + 0x09;    s1p14i = tm2 + 0x09;    s1p24i = max_err + 0x09;    s1p34i = sse2_rnd + 0x09;
		s1p05r = tmp + 0x0a;    s1p15r = tm2 + 0x0a;    s1p25r = max_err + 0x0a;    s1p35r = sse2_rnd + 0x0a;
		s1p05i = tmp + 0x0b;    s1p15i = tm2 + 0x0b;    s1p25i = max_err + 0x0b;    s1p35i = sse2_rnd + 0x0b;
		s1p06r = tmp + 0x0c;    s1p16r = tm2 + 0x0c;    s1p26r = max_err + 0x0c;    s1p36r = sse2_rnd + 0x0c;
		s1p06i = tmp + 0x0d;    s1p16i = tm2 + 0x0d;    s1p26i = max_err + 0x0d;    s1p36i = sse2_rnd + 0x0d;
		s1p07r = tmp + 0x0e;    s1p17r = tm2 + 0x0e;    s1p27r = max_err + 0x0e;    s1p37r = sse2_rnd + 0x0e;
		s1p07i = tmp + 0x0f;    s1p17i = tm2 + 0x0f;    s1p27i = max_err + 0x0f;    s1p37i = sse2_rnd + 0x0f;
		s1p08r = tmp + 0x10;    s1p18r = tm2 + 0x10;    s1p28r = max_err + 0x10;    s1p38r = sse2_rnd + 0x10;
		s1p08i = tmp + 0x11;    s1p18i = tm2 + 0x11;    s1p28i = max_err + 0x11;    s1p38i = sse2_rnd + 0x11;
		s1p09r = tmp + 0x12;    s1p19r = tm2 + 0x12;    s1p29r = max_err + 0x12;    s1p39r = sse2_rnd + 0x12;
		s1p09i = tmp + 0x13;    s1p19i = tm2 + 0x13;    s1p29i = max_err + 0x13;    s1p39i = sse2_rnd + 0x13;
		s1p0ar = tmp + 0x14;    s1p1ar = tm2 + 0x14;    s1p2ar = max_err + 0x14;    s1p3ar = sse2_rnd + 0x14;
		s1p0ai = tmp + 0x15;    s1p1ai = tm2 + 0x15;    s1p2ai = max_err + 0x15;    s1p3ai = sse2_rnd + 0x15;
		s1p0br = tmp + 0x16;    s1p1br = tm2 + 0x16;    s1p2br = max_err + 0x16;    s1p3br = sse2_rnd + 0x16;
		s1p0bi = tmp + 0x17;    s1p1bi = tm2 + 0x17;    s1p2bi = max_err + 0x17;    s1p3bi = sse2_rnd + 0x17;
		s1p0cr = tmp + 0x18;    s1p1cr = tm2 + 0x18;    s1p2cr = max_err + 0x18;    s1p3cr = sse2_rnd + 0x18;
		s1p0ci = tmp + 0x19;    s1p1ci = tm2 + 0x19;    s1p2ci = max_err + 0x19;    s1p3ci = sse2_rnd + 0x19;
		s1p0dr = tmp + 0x1a;    s1p1dr = tm2 + 0x1a;    s1p2dr = max_err + 0x1a;    s1p3dr = sse2_rnd + 0x1a;
		s1p0di = tmp + 0x1b;    s1p1di = tm2 + 0x1b;    s1p2di = max_err + 0x1b;    s1p3di = sse2_rnd + 0x1b;
		s1p0er = tmp + 0x1c;    s1p1er = tm2 + 0x1c;    s1p2er = max_err + 0x1c;    s1p3er = sse2_rnd + 0x1c;
		s1p0ei = tmp + 0x1d;    s1p1ei = tm2 + 0x1d;    s1p2ei = max_err + 0x1d;    s1p3ei = sse2_rnd + 0x1d;
		tmp = sse2_rnd + 0x1e;  tm2 =    tmp + 0x1e;    max_err  = tm2 + 0x1e;    sse2_rnd = max_err + 0x1e;
		r00    = tmp + 0x00;    r10    = tm2 + 0x00;    r20    = max_err + 0x00;    r30    = sse2_rnd + 0x00;
		r01    = tmp + 0x02;    r11    = tm2 + 0x02;    r21    = max_err + 0x02;    r31    = sse2_rnd + 0x02;
		r02    = tmp + 0x04;    r12    = tm2 + 0x04;    r22    = max_err + 0x04;    r32    = sse2_rnd + 0x04;
		r03    = tmp + 0x06;    r13    = tm2 + 0x06;    r23    = max_err + 0x06;    r33    = sse2_rnd + 0x06;
		r04    = tmp + 0x08;    r14    = tm2 + 0x08;    r24    = max_err + 0x08;    r34    = sse2_rnd + 0x08;
		r05    = tmp + 0x0a;    r15    = tm2 + 0x0a;    r25    = max_err + 0x0a;    r35    = sse2_rnd + 0x0a;
		r06    = tmp + 0x0c;    r16    = tm2 + 0x0c;    r26    = max_err + 0x0c;    r36    = sse2_rnd + 0x0c;
		r07    = tmp + 0x0e;    r17    = tm2 + 0x0e;    r27    = max_err + 0x0e;    r37    = sse2_rnd + 0x0e;
		r08    = tmp + 0x10;    r18    = tm2 + 0x10;    r28    = max_err + 0x10;    r38    = sse2_rnd + 0x10;
		r09    = tmp + 0x12;    r19    = tm2 + 0x12;    r29    = max_err + 0x12;    r39    = sse2_rnd + 0x12;
		r0a    = tmp + 0x14;    r1a    = tm2 + 0x14;    r2a    = max_err + 0x14;    r3a    = sse2_rnd + 0x14;
		r0b    = tmp + 0x16;    r1b    = tm2 + 0x16;    r2b    = max_err + 0x16;    r3b    = sse2_rnd + 0x16;
		r0c    = tmp + 0x18;    r1c    = tm2 + 0x18;    r2c    = max_err + 0x18;    r3c    = sse2_rnd + 0x18;
		r0d    = tmp + 0x1a;    r1d    = tm2 + 0x1a;    r2d    = max_err + 0x1a;    r3d    = sse2_rnd + 0x1a;
		r0e    = tmp + 0x1c;    r1e    = tm2 + 0x1c;    r2e    = max_err + 0x1c;    r3e    = sse2_rnd + 0x1c;
		tmp = sse2_rnd + 0x1e;
		x00    = tmp + 0x00;
		x01    = tmp + 0x02;
		x02    = tmp + 0x04;
		x03    = tmp + 0x06;
		x04    = tmp + 0x08;
		x05    = tmp + 0x0a;
		x06    = tmp + 0x0c;
		x07    = tmp + 0x0e;
		x08    = tmp + 0x10;
		x09    = tmp + 0x12;
		x0a    = tmp + 0x14;
		x0b    = tmp + 0x16;
		x0c    = tmp + 0x18;
		x0d    = tmp + 0x1a;
		x0e    = tmp + 0x1c;
		tmp += 0x1e;
		y00    = tmp + 0x00;
		y01    = tmp + 0x02;
		y02    = tmp + 0x04;
		y03    = tmp + 0x06;
		y04    = tmp + 0x08;
		y05    = tmp + 0x0a;
		y06    = tmp + 0x0c;
		y07    = tmp + 0x0e;
		y08    = tmp + 0x10;
		y09    = tmp + 0x12;
		y0a    = tmp + 0x14;
		y0b    = tmp + 0x16;
		y0c    = tmp + 0x18;
		y0d    = tmp + 0x1a;
		y0e    = tmp + 0x1c;
		tmp += 0x1e;

		ASSERT(HERE, (tmp->d0 == tmp->d1) && (tmp->d0 == c3m1), "thread-local memcheck failed!");
		sse2_c3m1 = tmp + 0x00;
		sse2_s    = tmp + 0x01;
		sse2_cn1  = tmp + 0x02;
		sse2_cn2  = tmp + 0x03;
		sse2_ss3  = tmp + 0x04;
		sse2_sn1  = tmp + 0x05;
		sse2_sn2  = tmp + 0x06;
		tmp += 0x08;
	   	tm2 = tmp + 0x1e;	/* Use these as handy temps */

	  #ifdef USE_AVX
		cy_r00 = tmp + 0x00;
		cy_r04 = tmp + 0x01;
		cy_r08 = tmp + 0x02;
		cy_r12 = tmp + 0x03;
		cy_r16 = tmp + 0x04;
		cy_r20 = tmp + 0x05;
		cy_r24 = tmp + 0x06;
		cy_r28 = tmp + 0x07;
		cy_r32 = tmp + 0x08;
		cy_r36 = tmp + 0x09;
		cy_r40 = tmp + 0x0a;
		cy_r44 = tmp + 0x0b;
		cy_r48 = tmp + 0x0c;
		cy_r52 = tmp + 0x0d;
		cy_r56 = tmp + 0x0e;
		cy_i00 = tmp + 0x0f;
		cy_i04 = tmp + 0x10;
		cy_i08 = tmp + 0x11;
		cy_i12 = tmp + 0x12;
		cy_i16 = tmp + 0x13;
		cy_i20 = tmp + 0x14;
		cy_i24 = tmp + 0x15;
		cy_i28 = tmp + 0x16;
		cy_i32 = tmp + 0x17;
		cy_i36 = tmp + 0x18;
		cy_i40 = tmp + 0x19;
		cy_i44 = tmp + 0x1a;
		cy_i48 = tmp + 0x1b;
		cy_i52 = tmp + 0x1c;
		cy_i56 = tmp + 0x1d;	// +30 = 338 complex
		tmp += 0x1e;
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 340 complex
		// This is where the value of half_arr_offset60 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 vec_dbl for Mersenne-mod, and 3.5*RADIX[avx] | RADIX[sse2] for Fermat-mod */
								// +20 = 390 complex, round up to nearby multiple of 4
		base_negacyclic_root = half_arr + RADIX;	// Only used for Fermat-mod
	  #else
		cy_r00 = tmp + 0x00;	cy_r02 = tmp + 0x01;
		cy_r04 = tmp + 0x02;	cy_r06 = tmp + 0x03;
		cy_r08 = tmp + 0x04;	cy_r10 = tmp + 0x05;
		cy_r12 = tmp + 0x06;	cy_r14 = tmp + 0x07;
		cy_r16 = tmp + 0x08;	cy_r18 = tmp + 0x09;
		cy_r20 = tmp + 0x0a;	cy_r22 = tmp + 0x0b;
		cy_r24 = tmp + 0x0c;	cy_r26 = tmp + 0x0d;
		cy_r28 = tmp + 0x0e;	cy_r30 = tmp + 0x0f;
		cy_r32 = tmp + 0x10;	cy_r34 = tmp + 0x11;
		cy_r36 = tmp + 0x12;	cy_r38 = tmp + 0x13;
		cy_r40 = tmp + 0x14;	cy_r42 = tmp + 0x15;
		cy_r44 = tmp + 0x16;	cy_r46 = tmp + 0x17;
		cy_r48 = tmp + 0x18;	cy_r50 = tmp + 0x19;
		cy_r52 = tmp + 0x1a;	cy_r54 = tmp + 0x1b;
		cy_r56 = tmp + 0x1c;	cy_r58 = tmp + 0x1d;
		tmp += 0x1e;
		cy_i00 = tmp + 0x00;	cy_i02 = tmp + 0x01;
		cy_i04 = tmp + 0x02;	cy_i06 = tmp + 0x03;
		cy_i08 = tmp + 0x04;	cy_i10 = tmp + 0x05;
		cy_i12 = tmp + 0x06;	cy_i14 = tmp + 0x07;
		cy_i16 = tmp + 0x08;	cy_i18 = tmp + 0x09;
		cy_i20 = tmp + 0x0a;	cy_i22 = tmp + 0x0b;
		cy_i24 = tmp + 0x0c;	cy_i26 = tmp + 0x0d;
		cy_i28 = tmp + 0x0e;	cy_i30 = tmp + 0x0f;
		cy_i32 = tmp + 0x10;	cy_i34 = tmp + 0x11;
		cy_i36 = tmp + 0x12;	cy_i38 = tmp + 0x13;
		cy_i40 = tmp + 0x14;	cy_i42 = tmp + 0x15;
		cy_i44 = tmp + 0x16;	cy_i46 = tmp + 0x17;
		cy_i48 = tmp + 0x18;	cy_i50 = tmp + 0x19;
		cy_i52 = tmp + 0x1a;	cy_i54 = tmp + 0x1b;
		cy_i56 = tmp + 0x1c;	cy_i58 = tmp + 0x1d;	// +60 = 368 complex
		tmp += 0x1e;
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;	// +2 = 370 complex
		// This is where the value of half_arr_offset60 comes from
		half_arr= tmp + 0x02;	/* This table needs 20 x 16 bytes for Mersenne-mod, and [4*odd_radix] x 16 for Fermat-mod */
								// +20 = 390 complex, round up to nearby multiple of 4
	  #endif

		ASSERT(HERE, (s1p00r == thread_arg->s1p00r), "thread-local memcheck failed!");
		ASSERT(HERE, (r00 == thread_arg->r00), "thread-local memcheck failed!");
		ASSERT(HERE, (half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
		ASSERT(HERE, (sse2_rnd->d0 == crnd && sse2_rnd->d1 == crnd), "thread-local memcheck failed!");
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		tmp = half_arr;
	  #ifdef USE_AVX
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif
	} else {
		dtmp = (half_arr)->d0 * (half_arr+odd_radix)->d0;	ASSERT(HERE, fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
		dtmp = (half_arr)->d1 * (half_arr+odd_radix)->d1;	ASSERT(HERE, fabs(dtmp - scale) < EPS, "thread-local memcheck failed!");
	}

		VEC_DBL_INIT(max_err, 0.0);

		sign_mask = (uint64*)(s1p00r + radix60_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (#doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
		sse_sw  = sse_bw    + RE_IM_STRIDE;
		sse_n   = sse_sw    + RE_IM_STRIDE;
	  #ifdef USE_AVX
		n_minus_sil   = (struct uint32x4 *)sse_n + 1;
		n_minus_silp1 = (struct uint32x4 *)sse_n + 2;
		sinwt         = (struct uint32x4 *)sse_n + 3;
		sinwtm1       = (struct uint32x4 *)sse_n + 4;

		bjmodn00 = (int*)(sinwtm1 + RE_IM_STRIDE);
	  #else
		bjmodn00 = (int*)(sse_n + RE_IM_STRIDE);
	  #endif
		bjmodn01 = bjmodn00 +  1;
		bjmodn02 = bjmodn00 +  2;
		bjmodn03 = bjmodn00 +  3;
		bjmodn04 = bjmodn00 +  4;
		bjmodn05 = bjmodn00 +  5;
		bjmodn06 = bjmodn00 +  6;
		bjmodn07 = bjmodn00 +  7;
		bjmodn08 = bjmodn00 +  8;
		bjmodn09 = bjmodn00 +  9;
		bjmodn10 = bjmodn00 + 10;
		bjmodn11 = bjmodn00 + 11;
		bjmodn12 = bjmodn00 + 12;
		bjmodn13 = bjmodn00 + 13;
		bjmodn14 = bjmodn00 + 14;
		bjmodn15 = bjmodn00 + 15;
		bjmodn16 = bjmodn00 + 16;
		bjmodn17 = bjmodn00 + 17;
		bjmodn18 = bjmodn00 + 18;
		bjmodn19 = bjmodn00 + 19;
		bjmodn20 = bjmodn00 + 20;
		bjmodn21 = bjmodn00 + 21;
		bjmodn22 = bjmodn00 + 22;
		bjmodn23 = bjmodn00 + 23;
		bjmodn24 = bjmodn00 + 24;
		bjmodn25 = bjmodn00 + 25;
		bjmodn26 = bjmodn00 + 26;
		bjmodn27 = bjmodn00 + 27;
		bjmodn28 = bjmodn00 + 28;
		bjmodn29 = bjmodn00 + 29;
		bjmodn30 = bjmodn00 + 30;
		bjmodn31 = bjmodn00 + 31;
		bjmodn32 = bjmodn00 + 32;
		bjmodn33 = bjmodn00 + 33;
		bjmodn34 = bjmodn00 + 34;
		bjmodn35 = bjmodn00 + 35;
		bjmodn36 = bjmodn00 + 36;
		bjmodn37 = bjmodn00 + 37;
		bjmodn38 = bjmodn00 + 38;
		bjmodn39 = bjmodn00 + 39;
		bjmodn40 = bjmodn00 + 40;
		bjmodn41 = bjmodn00 + 41;
		bjmodn42 = bjmodn00 + 42;
		bjmodn43 = bjmodn00 + 43;
		bjmodn44 = bjmodn00 + 44;
		bjmodn45 = bjmodn00 + 45;
		bjmodn46 = bjmodn00 + 46;
		bjmodn47 = bjmodn00 + 47;
		bjmodn48 = bjmodn00 + 48;
		bjmodn49 = bjmodn00 + 49;
		bjmodn50 = bjmodn00 + 50;
		bjmodn51 = bjmodn00 + 51;
		bjmodn52 = bjmodn00 + 52;
		bjmodn53 = bjmodn00 + 53;
		bjmodn54 = bjmodn00 + 54;
		bjmodn55 = bjmodn00 + 55;
		bjmodn56 = bjmodn00 + 56;
		bjmodn57 = bjmodn00 + 57;
		bjmodn58 = bjmodn00 + 58;
		bjmodn59 = bjmodn00 + 59;
	#else

		// In scalar mode use these 2 ptrs to pass wts_idx_incr and the base/baseinv/etc array-ptrs:
		wts_idx_incr = *(int *)thread_arg->half_arr;
		base      = (double *)thread_arg->s1p00r;
		baseinv   = base + 2;
		wt_arr    = base + 4;
		wtinv_arr = wt_arr    + odd_radix;
		bs_arr    = wtinv_arr + odd_radix;
		bsinv_arr = bs_arr    + odd_radix;

	#endif	// USE_SSE2 ?

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			/* Init DWT-indices: */						/* init carries	*/
		#ifdef USE_AVX
			*bjmodn00 = thread_arg->bjmodn00;			cy_r00->d0 = thread_arg->cy_r00;
			*bjmodn01 = thread_arg->bjmodn01;			cy_r00->d1 = thread_arg->cy_r01;
			*bjmodn02 = thread_arg->bjmodn02;			cy_r00->d2 = thread_arg->cy_r02;
			*bjmodn03 = thread_arg->bjmodn03;			cy_r00->d3 = thread_arg->cy_r03;
			*bjmodn04 = thread_arg->bjmodn04;			cy_r04->d0 = thread_arg->cy_r04;
			*bjmodn05 = thread_arg->bjmodn05;			cy_r04->d1 = thread_arg->cy_r05;
			*bjmodn06 = thread_arg->bjmodn06;			cy_r04->d2 = thread_arg->cy_r06;
			*bjmodn07 = thread_arg->bjmodn07;			cy_r04->d3 = thread_arg->cy_r07;
			*bjmodn08 = thread_arg->bjmodn08;			cy_r08->d0 = thread_arg->cy_r08;
			*bjmodn09 = thread_arg->bjmodn09;			cy_r08->d1 = thread_arg->cy_r09;
			*bjmodn10 = thread_arg->bjmodn10;			cy_r08->d2 = thread_arg->cy_r10;
			*bjmodn11 = thread_arg->bjmodn11;			cy_r08->d3 = thread_arg->cy_r11;
			*bjmodn12 = thread_arg->bjmodn12;			cy_r12->d0 = thread_arg->cy_r12;
			*bjmodn13 = thread_arg->bjmodn13;			cy_r12->d1 = thread_arg->cy_r13;
			*bjmodn14 = thread_arg->bjmodn14;			cy_r12->d2 = thread_arg->cy_r14;
			*bjmodn15 = thread_arg->bjmodn15;			cy_r12->d3 = thread_arg->cy_r15;
			*bjmodn16 = thread_arg->bjmodn16;			cy_r16->d0 = thread_arg->cy_r16;
			*bjmodn17 = thread_arg->bjmodn17;			cy_r16->d1 = thread_arg->cy_r17;
			*bjmodn18 = thread_arg->bjmodn18;			cy_r16->d2 = thread_arg->cy_r18;
			*bjmodn19 = thread_arg->bjmodn19;			cy_r16->d3 = thread_arg->cy_r19;
			*bjmodn20 = thread_arg->bjmodn20;			cy_r20->d0 = thread_arg->cy_r20;
			*bjmodn21 = thread_arg->bjmodn21;			cy_r20->d1 = thread_arg->cy_r21;
			*bjmodn22 = thread_arg->bjmodn22;			cy_r20->d2 = thread_arg->cy_r22;
			*bjmodn23 = thread_arg->bjmodn23;			cy_r20->d3 = thread_arg->cy_r23;
			*bjmodn24 = thread_arg->bjmodn24;			cy_r24->d0 = thread_arg->cy_r24;
			*bjmodn25 = thread_arg->bjmodn25;			cy_r24->d1 = thread_arg->cy_r25;
			*bjmodn26 = thread_arg->bjmodn26;			cy_r24->d2 = thread_arg->cy_r26;
			*bjmodn27 = thread_arg->bjmodn27;			cy_r24->d3 = thread_arg->cy_r27;
			*bjmodn28 = thread_arg->bjmodn28;			cy_r28->d0 = thread_arg->cy_r28;
			*bjmodn29 = thread_arg->bjmodn29;			cy_r28->d1 = thread_arg->cy_r29;
			*bjmodn30 = thread_arg->bjmodn30;			cy_r28->d2 = thread_arg->cy_r30;
			*bjmodn31 = thread_arg->bjmodn31;			cy_r28->d3 = thread_arg->cy_r31;
			*bjmodn32 = thread_arg->bjmodn32;			cy_r32->d0 = thread_arg->cy_r32;
			*bjmodn33 = thread_arg->bjmodn33;			cy_r32->d1 = thread_arg->cy_r33;
			*bjmodn34 = thread_arg->bjmodn34;			cy_r32->d2 = thread_arg->cy_r34;
			*bjmodn35 = thread_arg->bjmodn35;			cy_r32->d3 = thread_arg->cy_r35;
			*bjmodn36 = thread_arg->bjmodn36;			cy_r36->d0 = thread_arg->cy_r36;
			*bjmodn37 = thread_arg->bjmodn37;			cy_r36->d1 = thread_arg->cy_r37;
			*bjmodn38 = thread_arg->bjmodn38;			cy_r36->d2 = thread_arg->cy_r38;
			*bjmodn39 = thread_arg->bjmodn39;			cy_r36->d3 = thread_arg->cy_r39;
			*bjmodn40 = thread_arg->bjmodn40;			cy_r40->d0 = thread_arg->cy_r40;
			*bjmodn41 = thread_arg->bjmodn41;			cy_r40->d1 = thread_arg->cy_r41;
			*bjmodn42 = thread_arg->bjmodn42;			cy_r40->d2 = thread_arg->cy_r42;
			*bjmodn43 = thread_arg->bjmodn43;			cy_r40->d3 = thread_arg->cy_r43;
			*bjmodn44 = thread_arg->bjmodn44;			cy_r44->d0 = thread_arg->cy_r44;
			*bjmodn45 = thread_arg->bjmodn45;			cy_r44->d1 = thread_arg->cy_r45;
			*bjmodn46 = thread_arg->bjmodn46;			cy_r44->d2 = thread_arg->cy_r46;
			*bjmodn47 = thread_arg->bjmodn47;			cy_r44->d3 = thread_arg->cy_r47;
			*bjmodn48 = thread_arg->bjmodn48;			cy_r48->d0 = thread_arg->cy_r48;
			*bjmodn49 = thread_arg->bjmodn49;			cy_r48->d1 = thread_arg->cy_r49;
			*bjmodn50 = thread_arg->bjmodn50;			cy_r48->d2 = thread_arg->cy_r50;
			*bjmodn51 = thread_arg->bjmodn51;			cy_r48->d3 = thread_arg->cy_r51;
			*bjmodn52 = thread_arg->bjmodn52;			cy_r52->d0 = thread_arg->cy_r52;
			*bjmodn53 = thread_arg->bjmodn53;			cy_r52->d1 = thread_arg->cy_r53;
			*bjmodn54 = thread_arg->bjmodn54;			cy_r52->d2 = thread_arg->cy_r54;
			*bjmodn55 = thread_arg->bjmodn55;			cy_r52->d3 = thread_arg->cy_r55;
			*bjmodn56 = thread_arg->bjmodn56;			cy_r56->d0 = thread_arg->cy_r56;
			*bjmodn57 = thread_arg->bjmodn57;			cy_r56->d1 = thread_arg->cy_r57;
			*bjmodn58 = thread_arg->bjmodn58;			cy_r56->d2 = thread_arg->cy_r58;
			*bjmodn59 = thread_arg->bjmodn59;			cy_r56->d3 = thread_arg->cy_r59;

		#elif defined(USE_SSE2)

			*bjmodn00 = thread_arg->bjmodn00;			cy_r00->d0 = thread_arg->cy_r00;
			*bjmodn01 = thread_arg->bjmodn01;			cy_r00->d1 = thread_arg->cy_r01;
			*bjmodn02 = thread_arg->bjmodn02;			cy_r02->d0 = thread_arg->cy_r02;
			*bjmodn03 = thread_arg->bjmodn03;			cy_r02->d1 = thread_arg->cy_r03;
			*bjmodn04 = thread_arg->bjmodn04;			cy_r04->d0 = thread_arg->cy_r04;
			*bjmodn05 = thread_arg->bjmodn05;			cy_r04->d1 = thread_arg->cy_r05;
			*bjmodn06 = thread_arg->bjmodn06;			cy_r06->d0 = thread_arg->cy_r06;
			*bjmodn07 = thread_arg->bjmodn07;			cy_r06->d1 = thread_arg->cy_r07;
			*bjmodn08 = thread_arg->bjmodn08;			cy_r08->d0 = thread_arg->cy_r08;
			*bjmodn09 = thread_arg->bjmodn09;			cy_r08->d1 = thread_arg->cy_r09;
			*bjmodn10 = thread_arg->bjmodn10;			cy_r10->d0 = thread_arg->cy_r10;
			*bjmodn11 = thread_arg->bjmodn11;			cy_r10->d1 = thread_arg->cy_r11;
			*bjmodn12 = thread_arg->bjmodn12;			cy_r12->d0 = thread_arg->cy_r12;
			*bjmodn13 = thread_arg->bjmodn13;			cy_r12->d1 = thread_arg->cy_r13;
			*bjmodn14 = thread_arg->bjmodn14;			cy_r14->d0 = thread_arg->cy_r14;
			*bjmodn15 = thread_arg->bjmodn15;			cy_r14->d1 = thread_arg->cy_r15;
			*bjmodn16 = thread_arg->bjmodn16;			cy_r16->d0 = thread_arg->cy_r16;
			*bjmodn17 = thread_arg->bjmodn17;			cy_r16->d1 = thread_arg->cy_r17;
			*bjmodn18 = thread_arg->bjmodn18;			cy_r18->d0 = thread_arg->cy_r18;
			*bjmodn19 = thread_arg->bjmodn19;			cy_r18->d1 = thread_arg->cy_r19;
			*bjmodn20 = thread_arg->bjmodn20;			cy_r20->d0 = thread_arg->cy_r20;
			*bjmodn21 = thread_arg->bjmodn21;			cy_r20->d1 = thread_arg->cy_r21;
			*bjmodn22 = thread_arg->bjmodn22;			cy_r22->d0 = thread_arg->cy_r22;
			*bjmodn23 = thread_arg->bjmodn23;			cy_r22->d1 = thread_arg->cy_r23;
			*bjmodn24 = thread_arg->bjmodn24;			cy_r24->d0 = thread_arg->cy_r24;
			*bjmodn25 = thread_arg->bjmodn25;			cy_r24->d1 = thread_arg->cy_r25;
			*bjmodn26 = thread_arg->bjmodn26;			cy_r26->d0 = thread_arg->cy_r26;
			*bjmodn27 = thread_arg->bjmodn27;			cy_r26->d1 = thread_arg->cy_r27;
			*bjmodn28 = thread_arg->bjmodn28;			cy_r28->d0 = thread_arg->cy_r28;
			*bjmodn29 = thread_arg->bjmodn29;			cy_r28->d1 = thread_arg->cy_r29;
			*bjmodn30 = thread_arg->bjmodn30;			cy_r30->d0 = thread_arg->cy_r30;
			*bjmodn31 = thread_arg->bjmodn31;			cy_r30->d1 = thread_arg->cy_r31;
			*bjmodn32 = thread_arg->bjmodn32;			cy_r32->d0 = thread_arg->cy_r32;
			*bjmodn33 = thread_arg->bjmodn33;			cy_r32->d1 = thread_arg->cy_r33;
			*bjmodn34 = thread_arg->bjmodn34;			cy_r34->d0 = thread_arg->cy_r34;
			*bjmodn35 = thread_arg->bjmodn35;			cy_r34->d1 = thread_arg->cy_r35;
			*bjmodn36 = thread_arg->bjmodn36;			cy_r36->d0 = thread_arg->cy_r36;
			*bjmodn37 = thread_arg->bjmodn37;			cy_r36->d1 = thread_arg->cy_r37;
			*bjmodn38 = thread_arg->bjmodn38;			cy_r38->d0 = thread_arg->cy_r38;
			*bjmodn39 = thread_arg->bjmodn39;			cy_r38->d1 = thread_arg->cy_r39;
			*bjmodn40 = thread_arg->bjmodn40;			cy_r40->d0 = thread_arg->cy_r40;
			*bjmodn41 = thread_arg->bjmodn41;			cy_r40->d1 = thread_arg->cy_r41;
			*bjmodn42 = thread_arg->bjmodn42;			cy_r42->d0 = thread_arg->cy_r42;
			*bjmodn43 = thread_arg->bjmodn43;			cy_r42->d1 = thread_arg->cy_r43;
			*bjmodn44 = thread_arg->bjmodn44;			cy_r44->d0 = thread_arg->cy_r44;
			*bjmodn45 = thread_arg->bjmodn45;			cy_r44->d1 = thread_arg->cy_r45;
			*bjmodn46 = thread_arg->bjmodn46;			cy_r46->d0 = thread_arg->cy_r46;
			*bjmodn47 = thread_arg->bjmodn47;			cy_r46->d1 = thread_arg->cy_r47;
			*bjmodn48 = thread_arg->bjmodn48;			cy_r48->d0 = thread_arg->cy_r48;
			*bjmodn49 = thread_arg->bjmodn49;			cy_r48->d1 = thread_arg->cy_r49;
			*bjmodn50 = thread_arg->bjmodn50;			cy_r50->d0 = thread_arg->cy_r50;
			*bjmodn51 = thread_arg->bjmodn51;			cy_r50->d1 = thread_arg->cy_r51;
			*bjmodn52 = thread_arg->bjmodn52;			cy_r52->d0 = thread_arg->cy_r52;
			*bjmodn53 = thread_arg->bjmodn53;			cy_r52->d1 = thread_arg->cy_r53;
			*bjmodn54 = thread_arg->bjmodn54;			cy_r54->d0 = thread_arg->cy_r54;
			*bjmodn55 = thread_arg->bjmodn55;			cy_r54->d1 = thread_arg->cy_r55;
			*bjmodn56 = thread_arg->bjmodn56;			cy_r56->d0 = thread_arg->cy_r56;
			*bjmodn57 = thread_arg->bjmodn57;			cy_r56->d1 = thread_arg->cy_r57;
			*bjmodn58 = thread_arg->bjmodn58;			cy_r58->d0 = thread_arg->cy_r58;
			*bjmodn59 = thread_arg->bjmodn59;			cy_r58->d1 = thread_arg->cy_r59;

		#else

			bjmodn00 = thread_arg->bjmodn00;			cy_r00 = thread_arg->cy_r00;
			bjmodn01 = thread_arg->bjmodn01;			cy_r01 = thread_arg->cy_r01;
			bjmodn02 = thread_arg->bjmodn02;			cy_r02 = thread_arg->cy_r02;
			bjmodn03 = thread_arg->bjmodn03;			cy_r03 = thread_arg->cy_r03;
			bjmodn04 = thread_arg->bjmodn04;			cy_r04 = thread_arg->cy_r04;
			bjmodn05 = thread_arg->bjmodn05;			cy_r05 = thread_arg->cy_r05;
			bjmodn06 = thread_arg->bjmodn06;			cy_r06 = thread_arg->cy_r06;
			bjmodn07 = thread_arg->bjmodn07;			cy_r07 = thread_arg->cy_r07;
			bjmodn08 = thread_arg->bjmodn08;			cy_r08 = thread_arg->cy_r08;
			bjmodn09 = thread_arg->bjmodn09;			cy_r09 = thread_arg->cy_r09;
			bjmodn10 = thread_arg->bjmodn10;			cy_r10 = thread_arg->cy_r10;
			bjmodn11 = thread_arg->bjmodn11;			cy_r11 = thread_arg->cy_r11;
			bjmodn12 = thread_arg->bjmodn12;			cy_r12 = thread_arg->cy_r12;
			bjmodn13 = thread_arg->bjmodn13;			cy_r13 = thread_arg->cy_r13;
			bjmodn14 = thread_arg->bjmodn14;			cy_r14 = thread_arg->cy_r14;
			bjmodn15 = thread_arg->bjmodn15;			cy_r15 = thread_arg->cy_r15;
			bjmodn16 = thread_arg->bjmodn16;			cy_r16 = thread_arg->cy_r16;
			bjmodn17 = thread_arg->bjmodn17;			cy_r17 = thread_arg->cy_r17;
			bjmodn18 = thread_arg->bjmodn18;			cy_r18 = thread_arg->cy_r18;
			bjmodn19 = thread_arg->bjmodn19;			cy_r19 = thread_arg->cy_r19;
			bjmodn20 = thread_arg->bjmodn20;			cy_r20 = thread_arg->cy_r20;
			bjmodn21 = thread_arg->bjmodn21;			cy_r21 = thread_arg->cy_r21;
			bjmodn22 = thread_arg->bjmodn22;			cy_r22 = thread_arg->cy_r22;
			bjmodn23 = thread_arg->bjmodn23;			cy_r23 = thread_arg->cy_r23;
			bjmodn24 = thread_arg->bjmodn24;			cy_r24 = thread_arg->cy_r24;
			bjmodn25 = thread_arg->bjmodn25;			cy_r25 = thread_arg->cy_r25;
			bjmodn26 = thread_arg->bjmodn26;			cy_r26 = thread_arg->cy_r26;
			bjmodn27 = thread_arg->bjmodn27;			cy_r27 = thread_arg->cy_r27;
			bjmodn28 = thread_arg->bjmodn28;			cy_r28 = thread_arg->cy_r28;
			bjmodn29 = thread_arg->bjmodn29;			cy_r29 = thread_arg->cy_r29;
			bjmodn30 = thread_arg->bjmodn30;			cy_r30 = thread_arg->cy_r30;
			bjmodn31 = thread_arg->bjmodn31;			cy_r31 = thread_arg->cy_r31;
			bjmodn32 = thread_arg->bjmodn32;			cy_r32 = thread_arg->cy_r32;
			bjmodn33 = thread_arg->bjmodn33;			cy_r33 = thread_arg->cy_r33;
			bjmodn34 = thread_arg->bjmodn34;			cy_r34 = thread_arg->cy_r34;
			bjmodn35 = thread_arg->bjmodn35;			cy_r35 = thread_arg->cy_r35;
			bjmodn36 = thread_arg->bjmodn36;			cy_r36 = thread_arg->cy_r36;
			bjmodn37 = thread_arg->bjmodn37;			cy_r37 = thread_arg->cy_r37;
			bjmodn38 = thread_arg->bjmodn38;			cy_r38 = thread_arg->cy_r38;
			bjmodn39 = thread_arg->bjmodn39;			cy_r39 = thread_arg->cy_r39;
			bjmodn40 = thread_arg->bjmodn40;			cy_r40 = thread_arg->cy_r40;
			bjmodn41 = thread_arg->bjmodn41;			cy_r41 = thread_arg->cy_r41;
			bjmodn42 = thread_arg->bjmodn42;			cy_r42 = thread_arg->cy_r42;
			bjmodn43 = thread_arg->bjmodn43;			cy_r43 = thread_arg->cy_r43;
			bjmodn44 = thread_arg->bjmodn44;			cy_r44 = thread_arg->cy_r44;
			bjmodn45 = thread_arg->bjmodn45;			cy_r45 = thread_arg->cy_r45;
			bjmodn46 = thread_arg->bjmodn46;			cy_r46 = thread_arg->cy_r46;
			bjmodn47 = thread_arg->bjmodn47;			cy_r47 = thread_arg->cy_r47;
			bjmodn48 = thread_arg->bjmodn48;			cy_r48 = thread_arg->cy_r48;
			bjmodn49 = thread_arg->bjmodn49;			cy_r49 = thread_arg->cy_r49;
			bjmodn50 = thread_arg->bjmodn50;			cy_r50 = thread_arg->cy_r50;
			bjmodn51 = thread_arg->bjmodn51;			cy_r51 = thread_arg->cy_r51;
			bjmodn52 = thread_arg->bjmodn52;			cy_r52 = thread_arg->cy_r52;
			bjmodn53 = thread_arg->bjmodn53;			cy_r53 = thread_arg->cy_r53;
			bjmodn54 = thread_arg->bjmodn54;			cy_r54 = thread_arg->cy_r54;
			bjmodn55 = thread_arg->bjmodn55;			cy_r55 = thread_arg->cy_r55;
			bjmodn56 = thread_arg->bjmodn56;			cy_r56 = thread_arg->cy_r56;
			bjmodn57 = thread_arg->bjmodn57;			cy_r57 = thread_arg->cy_r57;
			bjmodn58 = thread_arg->bjmodn58;			cy_r58 = thread_arg->cy_r58;
			bjmodn59 = thread_arg->bjmodn59;			cy_r59 = thread_arg->cy_r59;

		#endif
		}
		else	/* Fermat-mod uses "double helix" carry scheme - 2 separate sets of real/imaginary carries for right-angle transform, plus "twisted" wraparound step. */
		{
		#ifdef USE_AVX

			cy_r00->d0 = thread_arg->cy_r00;	cy_i00->d0 = thread_arg->cy_i00;
			cy_r00->d1 = thread_arg->cy_r01;	cy_i00->d1 = thread_arg->cy_i01;
			cy_r00->d2 = thread_arg->cy_r02;	cy_i00->d2 = thread_arg->cy_i02;
			cy_r00->d3 = thread_arg->cy_r03;	cy_i00->d3 = thread_arg->cy_i03;
			cy_r04->d0 = thread_arg->cy_r04;	cy_i04->d0 = thread_arg->cy_i04;
			cy_r04->d1 = thread_arg->cy_r05;	cy_i04->d1 = thread_arg->cy_i05;
			cy_r04->d2 = thread_arg->cy_r06;	cy_i04->d2 = thread_arg->cy_i06;
			cy_r04->d3 = thread_arg->cy_r07;	cy_i04->d3 = thread_arg->cy_i07;
			cy_r08->d0 = thread_arg->cy_r08;	cy_i08->d0 = thread_arg->cy_i08;
			cy_r08->d1 = thread_arg->cy_r09;	cy_i08->d1 = thread_arg->cy_i09;
			cy_r08->d2 = thread_arg->cy_r10;	cy_i08->d2 = thread_arg->cy_i10;
			cy_r08->d3 = thread_arg->cy_r11;	cy_i08->d3 = thread_arg->cy_i11;
			cy_r12->d0 = thread_arg->cy_r12;	cy_i12->d0 = thread_arg->cy_i12;
			cy_r12->d1 = thread_arg->cy_r13;	cy_i12->d1 = thread_arg->cy_i13;
			cy_r12->d2 = thread_arg->cy_r14;	cy_i12->d2 = thread_arg->cy_i14;
			cy_r12->d3 = thread_arg->cy_r15;	cy_i12->d3 = thread_arg->cy_i15;
			cy_r16->d0 = thread_arg->cy_r16;	cy_i16->d0 = thread_arg->cy_i16;
			cy_r16->d1 = thread_arg->cy_r17;	cy_i16->d1 = thread_arg->cy_i17;
			cy_r16->d2 = thread_arg->cy_r18;	cy_i16->d2 = thread_arg->cy_i18;
			cy_r16->d3 = thread_arg->cy_r19;	cy_i16->d3 = thread_arg->cy_i19;
			cy_r20->d0 = thread_arg->cy_r20;	cy_i20->d0 = thread_arg->cy_i20;
			cy_r20->d1 = thread_arg->cy_r21;	cy_i20->d1 = thread_arg->cy_i21;
			cy_r20->d2 = thread_arg->cy_r22;	cy_i20->d2 = thread_arg->cy_i22;
			cy_r20->d3 = thread_arg->cy_r23;	cy_i20->d3 = thread_arg->cy_i23;
			cy_r24->d0 = thread_arg->cy_r24;	cy_i24->d0 = thread_arg->cy_i24;
			cy_r24->d1 = thread_arg->cy_r25;	cy_i24->d1 = thread_arg->cy_i25;
			cy_r24->d2 = thread_arg->cy_r26;	cy_i24->d2 = thread_arg->cy_i26;
			cy_r24->d3 = thread_arg->cy_r27;	cy_i24->d3 = thread_arg->cy_i27;
			cy_r28->d0 = thread_arg->cy_r28;	cy_i28->d0 = thread_arg->cy_i28;
			cy_r28->d1 = thread_arg->cy_r29;	cy_i28->d1 = thread_arg->cy_i29;
			cy_r28->d2 = thread_arg->cy_r30;	cy_i28->d2 = thread_arg->cy_i30;
			cy_r28->d3 = thread_arg->cy_r31;	cy_i28->d3 = thread_arg->cy_i31;
			cy_r32->d0 = thread_arg->cy_r32;	cy_i32->d0 = thread_arg->cy_i32;
			cy_r32->d1 = thread_arg->cy_r33;	cy_i32->d1 = thread_arg->cy_i33;
			cy_r32->d2 = thread_arg->cy_r34;	cy_i32->d2 = thread_arg->cy_i34;
			cy_r32->d3 = thread_arg->cy_r35;	cy_i32->d3 = thread_arg->cy_i35;
			cy_r36->d0 = thread_arg->cy_r36;	cy_i36->d0 = thread_arg->cy_i36;
			cy_r36->d1 = thread_arg->cy_r37;	cy_i36->d1 = thread_arg->cy_i37;
			cy_r36->d2 = thread_arg->cy_r38;	cy_i36->d2 = thread_arg->cy_i38;
			cy_r36->d3 = thread_arg->cy_r39;	cy_i36->d3 = thread_arg->cy_i39;
			cy_r40->d0 = thread_arg->cy_r40;	cy_i40->d0 = thread_arg->cy_i40;
			cy_r40->d1 = thread_arg->cy_r41;	cy_i40->d1 = thread_arg->cy_i41;
			cy_r40->d2 = thread_arg->cy_r42;	cy_i40->d2 = thread_arg->cy_i42;
			cy_r40->d3 = thread_arg->cy_r43;	cy_i40->d3 = thread_arg->cy_i43;
			cy_r44->d0 = thread_arg->cy_r44;	cy_i44->d0 = thread_arg->cy_i44;
			cy_r44->d1 = thread_arg->cy_r45;	cy_i44->d1 = thread_arg->cy_i45;
			cy_r44->d2 = thread_arg->cy_r46;	cy_i44->d2 = thread_arg->cy_i46;
			cy_r44->d3 = thread_arg->cy_r47;	cy_i44->d3 = thread_arg->cy_i47;
			cy_r48->d0 = thread_arg->cy_r48;	cy_i48->d0 = thread_arg->cy_i48;
			cy_r48->d1 = thread_arg->cy_r49;	cy_i48->d1 = thread_arg->cy_i49;
			cy_r48->d2 = thread_arg->cy_r50;	cy_i48->d2 = thread_arg->cy_i50;
			cy_r48->d3 = thread_arg->cy_r51;	cy_i48->d3 = thread_arg->cy_i51;
			cy_r52->d0 = thread_arg->cy_r52;	cy_i52->d0 = thread_arg->cy_i52;
			cy_r52->d1 = thread_arg->cy_r53;	cy_i52->d1 = thread_arg->cy_i53;
			cy_r52->d2 = thread_arg->cy_r54;	cy_i52->d2 = thread_arg->cy_i54;
			cy_r52->d3 = thread_arg->cy_r55;	cy_i52->d3 = thread_arg->cy_i55;
			cy_r56->d0 = thread_arg->cy_r56;	cy_i56->d0 = thread_arg->cy_i56;
			cy_r56->d1 = thread_arg->cy_r57;	cy_i56->d1 = thread_arg->cy_i57;
			cy_r56->d2 = thread_arg->cy_r58;	cy_i56->d2 = thread_arg->cy_i58;
			cy_r56->d3 = thread_arg->cy_r59;	cy_i56->d3 = thread_arg->cy_i59;

		#elif defined(USE_SSE2)

			cy_r00->d0 = thread_arg->cy_r00;	cy_r00->d1 = thread_arg->cy_i00;
			cy_r02->d0 = thread_arg->cy_r01;	cy_r02->d1 = thread_arg->cy_i01;
			cy_r04->d0 = thread_arg->cy_r02;	cy_r04->d1 = thread_arg->cy_i02;
			cy_r06->d0 = thread_arg->cy_r03;	cy_r06->d1 = thread_arg->cy_i03;
			cy_r08->d0 = thread_arg->cy_r04;	cy_r08->d1 = thread_arg->cy_i04;
			cy_r10->d0 = thread_arg->cy_r05;	cy_r10->d1 = thread_arg->cy_i05;
			cy_r12->d0 = thread_arg->cy_r06;	cy_r12->d1 = thread_arg->cy_i06;
			cy_r14->d0 = thread_arg->cy_r07;	cy_r14->d1 = thread_arg->cy_i07;
			cy_r16->d0 = thread_arg->cy_r08;	cy_r16->d1 = thread_arg->cy_i08;
			cy_r18->d0 = thread_arg->cy_r09;	cy_r18->d1 = thread_arg->cy_i09;
			cy_r20->d0 = thread_arg->cy_r10;	cy_r20->d1 = thread_arg->cy_i10;
			cy_r22->d0 = thread_arg->cy_r11;	cy_r22->d1 = thread_arg->cy_i11;
			cy_r24->d0 = thread_arg->cy_r12;	cy_r24->d1 = thread_arg->cy_i12;
			cy_r26->d0 = thread_arg->cy_r13;	cy_r26->d1 = thread_arg->cy_i13;
			cy_r28->d0 = thread_arg->cy_r14;	cy_r28->d1 = thread_arg->cy_i14;
			cy_r30->d0 = thread_arg->cy_r15;	cy_r30->d1 = thread_arg->cy_i15;
			cy_r32->d0 = thread_arg->cy_r16;	cy_r32->d1 = thread_arg->cy_i16;
			cy_r34->d0 = thread_arg->cy_r17;	cy_r34->d1 = thread_arg->cy_i17;
			cy_r36->d0 = thread_arg->cy_r18;	cy_r36->d1 = thread_arg->cy_i18;
			cy_r38->d0 = thread_arg->cy_r19;	cy_r38->d1 = thread_arg->cy_i19;
			cy_r40->d0 = thread_arg->cy_r20;	cy_r40->d1 = thread_arg->cy_i20;
			cy_r42->d0 = thread_arg->cy_r21;	cy_r42->d1 = thread_arg->cy_i21;
			cy_r44->d0 = thread_arg->cy_r22;	cy_r44->d1 = thread_arg->cy_i22;
			cy_r46->d0 = thread_arg->cy_r23;	cy_r46->d1 = thread_arg->cy_i23;
			cy_r48->d0 = thread_arg->cy_r24;	cy_r48->d1 = thread_arg->cy_i24;
			cy_r50->d0 = thread_arg->cy_r25;	cy_r50->d1 = thread_arg->cy_i25;
			cy_r52->d0 = thread_arg->cy_r26;	cy_r52->d1 = thread_arg->cy_i26;
			cy_r54->d0 = thread_arg->cy_r27;	cy_r54->d1 = thread_arg->cy_i27;
			cy_r56->d0 = thread_arg->cy_r28;	cy_r56->d1 = thread_arg->cy_i28;
			cy_r58->d0 = thread_arg->cy_r29;	cy_r58->d1 = thread_arg->cy_i29;
			cy_i00->d0 = thread_arg->cy_r30;	cy_i00->d1 = thread_arg->cy_i30;
			cy_i02->d0 = thread_arg->cy_r31;	cy_i02->d1 = thread_arg->cy_i31;
			cy_i04->d0 = thread_arg->cy_r32;	cy_i04->d1 = thread_arg->cy_i32;
			cy_i06->d0 = thread_arg->cy_r33;	cy_i06->d1 = thread_arg->cy_i33;
			cy_i08->d0 = thread_arg->cy_r34;	cy_i08->d1 = thread_arg->cy_i34;
			cy_i10->d0 = thread_arg->cy_r35;	cy_i10->d1 = thread_arg->cy_i35;
			cy_i12->d0 = thread_arg->cy_r36;	cy_i12->d1 = thread_arg->cy_i36;
			cy_i14->d0 = thread_arg->cy_r37;	cy_i14->d1 = thread_arg->cy_i37;
			cy_i16->d0 = thread_arg->cy_r38;	cy_i16->d1 = thread_arg->cy_i38;
			cy_i18->d0 = thread_arg->cy_r39;	cy_i18->d1 = thread_arg->cy_i39;
			cy_i20->d0 = thread_arg->cy_r40;	cy_i20->d1 = thread_arg->cy_i40;
			cy_i22->d0 = thread_arg->cy_r41;	cy_i22->d1 = thread_arg->cy_i41;
			cy_i24->d0 = thread_arg->cy_r42;	cy_i24->d1 = thread_arg->cy_i42;
			cy_i26->d0 = thread_arg->cy_r43;	cy_i26->d1 = thread_arg->cy_i43;
			cy_i28->d0 = thread_arg->cy_r44;	cy_i28->d1 = thread_arg->cy_i44;
			cy_i30->d0 = thread_arg->cy_r45;	cy_i30->d1 = thread_arg->cy_i45;
			cy_i32->d0 = thread_arg->cy_r46;	cy_i32->d1 = thread_arg->cy_i46;
			cy_i34->d0 = thread_arg->cy_r47;	cy_i34->d1 = thread_arg->cy_i47;
			cy_i36->d0 = thread_arg->cy_r48;	cy_i36->d1 = thread_arg->cy_i48;
			cy_i38->d0 = thread_arg->cy_r49;	cy_i38->d1 = thread_arg->cy_i49;
			cy_i40->d0 = thread_arg->cy_r50;	cy_i40->d1 = thread_arg->cy_i50;
			cy_i42->d0 = thread_arg->cy_r51;	cy_i42->d1 = thread_arg->cy_i51;
			cy_i44->d0 = thread_arg->cy_r52;	cy_i44->d1 = thread_arg->cy_i52;
			cy_i46->d0 = thread_arg->cy_r53;	cy_i46->d1 = thread_arg->cy_i53;
			cy_i48->d0 = thread_arg->cy_r54;	cy_i48->d1 = thread_arg->cy_i54;
			cy_i50->d0 = thread_arg->cy_r55;	cy_i50->d1 = thread_arg->cy_i55;
			cy_i52->d0 = thread_arg->cy_r56;	cy_i52->d1 = thread_arg->cy_i56;
			cy_i54->d0 = thread_arg->cy_r57;	cy_i54->d1 = thread_arg->cy_i57;
			cy_i56->d0 = thread_arg->cy_r58;	cy_i56->d1 = thread_arg->cy_i58;
			cy_i58->d0 = thread_arg->cy_r59;	cy_i58->d1 = thread_arg->cy_i59;

		#else

			cy_r00 = thread_arg->cy_r00;		cy_i00 = thread_arg->cy_i00;
			cy_r01 = thread_arg->cy_r01;		cy_i01 = thread_arg->cy_i01;
			cy_r02 = thread_arg->cy_r02;		cy_i02 = thread_arg->cy_i02;
			cy_r03 = thread_arg->cy_r03;		cy_i03 = thread_arg->cy_i03;
			cy_r04 = thread_arg->cy_r04;		cy_i04 = thread_arg->cy_i04;
			cy_r05 = thread_arg->cy_r05;		cy_i05 = thread_arg->cy_i05;
			cy_r06 = thread_arg->cy_r06;		cy_i06 = thread_arg->cy_i06;
			cy_r07 = thread_arg->cy_r07;		cy_i07 = thread_arg->cy_i07;
			cy_r08 = thread_arg->cy_r08;		cy_i08 = thread_arg->cy_i08;
			cy_r09 = thread_arg->cy_r09;		cy_i09 = thread_arg->cy_i09;
			cy_r10 = thread_arg->cy_r10;		cy_i10 = thread_arg->cy_i10;
			cy_r11 = thread_arg->cy_r11;		cy_i11 = thread_arg->cy_i11;
			cy_r12 = thread_arg->cy_r12;		cy_i12 = thread_arg->cy_i12;
			cy_r13 = thread_arg->cy_r13;		cy_i13 = thread_arg->cy_i13;
			cy_r14 = thread_arg->cy_r14;		cy_i14 = thread_arg->cy_i14;
			cy_r15 = thread_arg->cy_r15;		cy_i15 = thread_arg->cy_i15;
			cy_r16 = thread_arg->cy_r16;		cy_i16 = thread_arg->cy_i16;
			cy_r17 = thread_arg->cy_r17;		cy_i17 = thread_arg->cy_i17;
			cy_r18 = thread_arg->cy_r18;		cy_i18 = thread_arg->cy_i18;
			cy_r19 = thread_arg->cy_r19;		cy_i19 = thread_arg->cy_i19;
			cy_r20 = thread_arg->cy_r20;		cy_i20 = thread_arg->cy_i20;
			cy_r21 = thread_arg->cy_r21;		cy_i21 = thread_arg->cy_i21;
			cy_r22 = thread_arg->cy_r22;		cy_i22 = thread_arg->cy_i22;
			cy_r23 = thread_arg->cy_r23;		cy_i23 = thread_arg->cy_i23;
			cy_r24 = thread_arg->cy_r24;		cy_i24 = thread_arg->cy_i24;
			cy_r25 = thread_arg->cy_r25;		cy_i25 = thread_arg->cy_i25;
			cy_r26 = thread_arg->cy_r26;		cy_i26 = thread_arg->cy_i26;
			cy_r27 = thread_arg->cy_r27;		cy_i27 = thread_arg->cy_i27;
			cy_r28 = thread_arg->cy_r28;		cy_i28 = thread_arg->cy_i28;
			cy_r29 = thread_arg->cy_r29;		cy_i29 = thread_arg->cy_i29;
			cy_r30 = thread_arg->cy_r30;		cy_i30 = thread_arg->cy_i30;
			cy_r31 = thread_arg->cy_r31;		cy_i31 = thread_arg->cy_i31;
			cy_r32 = thread_arg->cy_r32;		cy_i32 = thread_arg->cy_i32;
			cy_r33 = thread_arg->cy_r33;		cy_i33 = thread_arg->cy_i33;
			cy_r34 = thread_arg->cy_r34;		cy_i34 = thread_arg->cy_i34;
			cy_r35 = thread_arg->cy_r35;		cy_i35 = thread_arg->cy_i35;
			cy_r36 = thread_arg->cy_r36;		cy_i36 = thread_arg->cy_i36;
			cy_r37 = thread_arg->cy_r37;		cy_i37 = thread_arg->cy_i37;
			cy_r38 = thread_arg->cy_r38;		cy_i38 = thread_arg->cy_i38;
			cy_r39 = thread_arg->cy_r39;		cy_i39 = thread_arg->cy_i39;
			cy_r40 = thread_arg->cy_r40;		cy_i40 = thread_arg->cy_i40;
			cy_r41 = thread_arg->cy_r41;		cy_i41 = thread_arg->cy_i41;
			cy_r42 = thread_arg->cy_r42;		cy_i42 = thread_arg->cy_i42;
			cy_r43 = thread_arg->cy_r43;		cy_i43 = thread_arg->cy_i43;
			cy_r44 = thread_arg->cy_r44;		cy_i44 = thread_arg->cy_i44;
			cy_r45 = thread_arg->cy_r45;		cy_i45 = thread_arg->cy_i45;
			cy_r46 = thread_arg->cy_r46;		cy_i46 = thread_arg->cy_i46;
			cy_r47 = thread_arg->cy_r47;		cy_i47 = thread_arg->cy_i47;
			cy_r48 = thread_arg->cy_r48;		cy_i48 = thread_arg->cy_i48;
			cy_r49 = thread_arg->cy_r49;		cy_i49 = thread_arg->cy_i49;
			cy_r50 = thread_arg->cy_r50;		cy_i50 = thread_arg->cy_i50;
			cy_r51 = thread_arg->cy_r51;		cy_i51 = thread_arg->cy_i51;
			cy_r52 = thread_arg->cy_r52;		cy_i52 = thread_arg->cy_i52;
			cy_r53 = thread_arg->cy_r53;		cy_i53 = thread_arg->cy_i53;
			cy_r54 = thread_arg->cy_r54;		cy_i54 = thread_arg->cy_i54;
			cy_r55 = thread_arg->cy_r55;		cy_i55 = thread_arg->cy_i55;
			cy_r56 = thread_arg->cy_r56;		cy_i56 = thread_arg->cy_i56;
			cy_r57 = thread_arg->cy_r57;		cy_i57 = thread_arg->cy_i57;
			cy_r58 = thread_arg->cy_r58;		cy_i58 = thread_arg->cy_i58;
			cy_r59 = thread_arg->cy_r59;		cy_i59 = thread_arg->cy_i59;

		#endif
		}

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
			for(j = jstart; j < jhi; j += stride)	// Stride = 4 reals for SSE2, 8 for AVX
			{
				j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1+RE_IM_STRIDE;

			#ifdef USE_AVX
			/* Outputs in SSE2 modes are temps 2*15*32 = 30*32 = 0x3c0 bytes apart: */
			// Reorder blocks to yield sequentially increasing a-array offsets:
				add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, 0x3c0)
				add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, 0x3c0)
				add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01, 0x3c0)
				add2 = &a[j1+p12];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05, 0x3c0)
				add1 = &a[j1+p16];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, 0x3c0)
				add3 = &a[j1+p20];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03, 0x3c0)
				add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07, 0x3c0)
				add3 = &a[j1+p28];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, 0x3c0)
				add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, 0x3c0)
				add3 = &a[j1+p36];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09, 0x3c0)
				add0 = &a[j1+p40];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0b, 0x3c0)
				add2 = &a[j1+p44];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, 0x3c0)
				add0 = &a[j1+p48];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, 0x3c0)
				add2 = &a[j1+p52];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0d, 0x3c0)
				add1 = &a[j1+p56];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, 0x3c0)

			#elif defined(USE_SSE2)

			/* Outputs in SSE2 modes are temps 2*15*16 = 30*16 = 0x1e0 bytes apart: */
			  #if !GCC_ASM_FULL_INLINE
			// Reorder blocks to yield sequentially increasing a-array offsets:
/* Block 01: */	add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, 0x1e0)
/* Block 03: */	add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, 0x1e0)
/* Block 02: */	add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01, 0x1e0)
/* Block 06: */	add2 = &a[j1+p12];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05, 0x1e0)
/* Block 05: */	add1 = &a[j1+p16];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, 0x1e0)
/* Block 04: */	add3 = &a[j1+p20];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03, 0x1e0)
/* Block 08: */	add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07, 0x1e0)
/* Block 07: */	add3 = &a[j1+p28];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, 0x1e0)
/* Block 09: */	add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, 0x1e0)
/* Block 10: */	add3 = &a[j1+p36];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09, 0x1e0)
/* Block 12: */	add0 = &a[j1+p40];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0b, 0x1e0)
/* Block 11: */	add2 = &a[j1+p44];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, 0x1e0)
/* Block 15: */	add0 = &a[j1+p48];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, 0x1e0)
/* Block 14: */	add2 = &a[j1+p52];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0d, 0x1e0)
/* Block 13: */	add1 = &a[j1+p56];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, 0x1e0)

			  #else

				add0 = &a[j1    ];
				SSE2_RADIX60_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,r00);

			  #endif

			  #if (OS_BITS == 32) || !defined(USE_64BIT_ASM_STYLE)	// In 64-bit mode, default is to use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just xmm0-7.
				/* Radix-15 DFT uses adjacent temps for inputs, outputs are (cyclic) with pXXr having XX += 12 (24*16 bytes = 0x180) between successive outputs: */
			   #ifndef USE_LITERAL_BYTE_OFFSETS
				SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1     ,  r00,  r01,  r02,  r03,  r04,  r05,  r06,  r07,  r08,  r09,  r0a,  r0b,  r0c,  r0d,  r0e     ,  x00,  x01,  x02,  x03,  x04,  x05,  x06,  x07,  x08,  x09,  x0a,  x0b,  x0c,  x0d,  x0e,        s1p00r,s1p2er,s1p1dr,s1p0cr,s1p3br,s1p2ar,s1p19r,s1p08r,s1p37r,s1p26r,s1p15r,s1p04r,s1p33r,s1p22r,s1p11r);
				SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1     ,  r10,  r11,  r12,  r13,  r14,  r15,  r16,  r17,  r18,  r19,  r1a,  r1b,  r1c,  r1d,  r1e     ,  x00,  x01,  x02,  x03,  x04,  x05,  x06,  x07,  x08,  x09,  x0a,  x0b,  x0c,  x0d,  x0e,        s1p10r,s1p3er,s1p2dr,s1p1cr,s1p0br,s1p3ar,s1p29r,s1p18r,s1p07r,s1p36r,s1p25r,s1p14r,s1p03r,s1p32r,s1p21r);
				SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1     ,  r20,  r21,  r22,  r23,  r24,  r25,  r26,  r27,  r28,  r29,  r2a,  r2b,  r2c,  r2d,  r2e     ,  x00,  x01,  x02,  x03,  x04,  x05,  x06,  x07,  x08,  x09,  x0a,  x0b,  x0c,  x0d,  x0e,        s1p20r,s1p0er,s1p3dr,s1p2cr,s1p1br,s1p0ar,s1p39r,s1p28r,s1p17r,s1p06r,s1p35r,s1p24r,s1p13r,s1p02r,s1p31r);
				SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1     ,  r30,  r31,  r32,  r33,  r34,  r35,  r36,  r37,  r38,  r39,  r3a,  r3b,  r3c,  r3d,  r3e     ,  x00,  x01,  x02,  x03,  x04,  x05,  x06,  x07,  x08,  x09,  x0a,  x0b,  x0c,  x0d,  x0e,        s1p30r,s1p1er,s1p0dr,s1p3cr,s1p2br,s1p1ar,s1p09r,s1p38r,s1p27r,s1p16r,s1p05r,s1p34r,s1p23r,s1p12r,s1p01r);
			   #else	// Versions using base-address-plus-literal-byte-offsets:                                                                                       s1p**r complex-element offsets w.r.to s1p00r (all negative offsets get added to +120):     0    +88    +56    +24    -8     -40    -72    -104   -16    -48    -80    -112   -24    -56    -88 ... when offset < -(row-start offset below), reset by adding 120*16 = 0x780 bytes
				SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1, r00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, s1p00r, 0x000, 0x580, 0x380, 0x180, 0x700, 0x500, 0x300, 0x100, 0x680, 0x480, 0x280, 0x080, 0x600, 0x400, 0x200);	// end offset of each s1p-sequence = start-offset + 0x200
				SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1, r10,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, s1p00r, 0x1e0, 0x760, 0x560, 0x360, 0x160, 0x6e0, 0x4e0, 0x2e0, 0x0e0, 0x660, 0x460, 0x260, 0x060, 0x5e0, 0x3e0);	// s1p10r = +30 complex = +0x1e0 bytes
				SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1, r20,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, s1p00r, 0x3c0, 0x1c0, 0x740, 0x540, 0x340, 0x140, 0x6c0, 0x4c0, 0x2c0, 0x0c0, 0x640, 0x440, 0x240, 0x040, 0x5c0);	// s1p20r = +60 complex = +0x3c0 bytes
				SSE2_RADIX_15_DIT(sse2_c3m1,sse2_cn1, r30,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, s1p00r, 0x5a0, 0x3a0, 0x1a0, 0x720, 0x520, 0x320, 0x120, 0x6a0, 0x4a0, 0x2a0, 0x0a0, 0x620, 0x420, 0x220, 0x020);	// s1p30r = +90 complex = +0x5a0 bytes
			   #endif
			  #else
				SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,\
					r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r0a,r0b,r0c,r0d,r0e,\
					x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,\
					s1p00r,s1p2er,s1p1dr,s1p0cr,s1p3br,s1p2ar,s1p19r,s1p08r,s1p37r,s1p26r,s1p15r,s1p04r,s1p33r,s1p22r,s1p11r,\
					r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r1a,r1b,r1c,r1d,r1e,\
					y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,\
					s1p10r,s1p3er,s1p2dr,s1p1cr,s1p0br,s1p3ar,s1p29r,s1p18r,s1p07r,s1p36r,s1p25r,s1p14r,s1p03r,s1p32r,s1p21r);
				SSE2_RADIX_15_DIT_X2(sse2_c3m1,sse2_cn1,\
					r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r2a,r2b,r2c,r2d,r2e,\
					x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e,\
					s1p20r,s1p0er,s1p3dr,s1p2cr,s1p1br,s1p0ar,s1p39r,s1p28r,s1p17r,s1p06r,s1p35r,s1p24r,s1p13r,s1p02r,s1p31r,\
					r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r3a,r3b,r3c,r3d,r3e,\
					y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e,\
					s1p30r,s1p1er,s1p0dr,s1p3cr,s1p2br,s1p1ar,s1p09r,s1p38r,s1p27r,s1p16r,s1p05r,s1p34r,s1p23r,s1p12r,s1p01r);
			  #endif

			#else	// Scalar-double mode:

				jt = j1    ; jp = j2    ;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,rt,it);
				jt = j1+p08; jp = j2+p08;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,rt,it);
				jt = j1+p04; jp = j2+p04;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,rt,it);
				jt = j1+p20; jp = j2+p20;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,rt,it);
				jt = j1+p16; jp = j2+p16;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,rt,it);
				jt = j1+p12; jp = j2+p12;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,rt,it);
				jt = j1+p28; jp = j2+p28;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,rt,it);
				jt = j1+p24; jp = j2+p24;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,rt,it);
				jt = j1+p32; jp = j2+p32;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,rt,it);
				jt = j1+p36; jp = j2+p36;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,rt,it);
				jt = j1+p44; jp = j2+p44;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,rt,it);
				jt = j1+p40; jp = j2+p40;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0br,t0bi,t1br,t1bi,t2br,t2bi,t3br,t3bi,rt,it);
				jt = j1+p56; jp = j2+p56;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t0cr,t0ci,t1cr,t1ci,t2cr,t2ci,t3cr,t3ci,rt,it);
				jt = j1+p52; jp = j2+p52;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t0dr,t0di,t1dr,t1di,t2dr,t2di,t3dr,t3di,rt,it);
				jt = j1+p48; jp = j2+p48;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0er,t0ei,t1er,t1ei,t2er,t2ei,t3er,t3ei,rt,it);

				RADIX_15_DIT(t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a1p00r,a1p00i,a1p2er,a1p2ei,a1p1dr,a1p1di,a1p0cr,a1p0ci,a1p3br,a1p3bi,a1p2ar,a1p2ai,a1p19r,a1p19i,a1p08r,a1p08i,a1p37r,a1p37i,a1p26r,a1p26i,a1p15r,a1p15i,a1p04r,a1p04i,a1p33r,a1p33i,a1p22r,a1p22i,a1p11r,a1p11i,rt,it);
				RADIX_15_DIT(t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a1p10r,a1p10i,a1p3er,a1p3ei,a1p2dr,a1p2di,a1p1cr,a1p1ci,a1p0br,a1p0bi,a1p3ar,a1p3ai,a1p29r,a1p29i,a1p18r,a1p18i,a1p07r,a1p07i,a1p36r,a1p36i,a1p25r,a1p25i,a1p14r,a1p14i,a1p03r,a1p03i,a1p32r,a1p32i,a1p21r,a1p21i,rt,it);
				RADIX_15_DIT(t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a1p20r,a1p20i,a1p0er,a1p0ei,a1p3dr,a1p3di,a1p2cr,a1p2ci,a1p1br,a1p1bi,a1p0ar,a1p0ai,a1p39r,a1p39i,a1p28r,a1p28i,a1p17r,a1p17i,a1p06r,a1p06i,a1p35r,a1p35i,a1p24r,a1p24i,a1p13r,a1p13i,a1p02r,a1p02i,a1p31r,a1p31i,rt,it);
				RADIX_15_DIT(t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,a1p30r,a1p30i,a1p1er,a1p1ei,a1p0dr,a1p0di,a1p3cr,a1p3ci,a1p2br,a1p2bi,a1p1ar,a1p1ai,a1p09r,a1p09i,a1p38r,a1p38i,a1p27r,a1p27i,a1p16r,a1p16i,a1p05r,a1p05i,a1p34r,a1p34i,a1p23r,a1p23i,a1p12r,a1p12i,a1p01r,a1p01i,rt,it);

			#endif	// SIMD or not?

			/*...Now do the carries. Since the outputs would
				normally be getting dispatched to 60 separate blocks of the A-array, we need 60 separate carries.	*/

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
			#ifdef USE_AVX

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

				l= j & (nwt-1);						tmp = half_arr + 64;	/* ptr to local storage for the doubled wtl,wtn terms: */
				n_minus_sil  ->d0 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d0 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d0 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d0 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+2) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d1 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d1 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d1 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d1 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+4) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d2 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d2 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d2 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d2 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				l= (j+6) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
				n_minus_sil  ->d3 = n-si[l  ];		tmp->d0 = wt0[    l  ];
				n_minus_silp1->d3 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
				sinwt        ->d3 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
				sinwtm1      ->d3 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

				AVX_cmplx_carry_norm_errcheck0_X4(s1p00r,add1,add2,add3,cy_r00,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p04r,add1,add2,add3,cy_r04,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p08r,add1,add2,add3,cy_r08,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p0cr,add1,add2,add3,cy_r12,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p11r,add1,add2,add3,cy_r16,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p15r,add1,add2,add3,cy_r20,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p19r,add1,add2,add3,cy_r24,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p1dr,add1,add2,add3,cy_r28,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p22r,add1,add2,add3,cy_r32,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p26r,add1,add2,add3,cy_r36,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p2ar,add1,add2,add3,cy_r40,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p2er,add1,add2,add3,cy_r44,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p33r,add1,add2,add3,cy_r48,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p37r,add1,add2,add3,cy_r52,bjmodn52,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p3br,add1,add2,add3,cy_r56,bjmodn56,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);

				co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

			#elif defined(USE_SSE2)

				l= j & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
				ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
				ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
				ctmp->re = wtnm1;	ctmp->im = wtnm1;

				add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p0cr,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p11r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p15r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p19r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p1dr,add1,add2,add3,cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p22r,add1,add2,add3,cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p26r,add1,add2,add3,cy_r36,cy_r38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p2ar,add1,add2,add3,cy_r40,cy_r42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p2er,add1,add2,add3,cy_r44,cy_r46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p33r,add1,add2,add3,cy_r48,cy_r50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p37r,add1,add2,add3,cy_r52,cy_r54,bjmodn52,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p3br,add1,add2,add3,cy_r56,cy_r58,bjmodn56,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			   #else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p0cr,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p11r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p15r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p19r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p1dr,add1,add2,add3,cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p22r,add1,add2,add3,cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p26r,add1,add2,add3,cy_r36,cy_r38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p2ar,add1,add2,add3,cy_r40,cy_r42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p2er,add1,add2,add3,cy_r44,cy_r46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p33r,add1,add2,add3,cy_r48,cy_r50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p37r,add1,add2,add3,cy_r52,cy_r54,bjmodn52,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p3br,add1,add2,add3,cy_r56,cy_r58,bjmodn56,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			   #endif

				l= (j+2) & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
				n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
				n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
				sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
				ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
				ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
				ctmp->re = wtnm1;	ctmp->im = wtnm1;

			/*	i =((uint32)(sw - *bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];

			   #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p0cr,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p11r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p15r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p19r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p1dr,add1,add2,     cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p22r,add1,add2,     cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p26r,add1,add2,     cy_r36,cy_r38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p2ar,add1,add2,     cy_r40,cy_r42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p2er,add1,add2,     cy_r44,cy_r46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p33r,add1,add2,     cy_r48,cy_r50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p37r,add1,add2,     cy_r52,cy_r54,bjmodn52,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p3br,add1,add2,     cy_r56,cy_r58,bjmodn56,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			   #else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p0cr,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p11r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p15r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p19r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p1dr,add1,add2,     cy_r28,cy_r30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p22r,add1,add2,     cy_r32,cy_r34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p26r,add1,add2,     cy_r36,cy_r38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p2ar,add1,add2,     cy_r40,cy_r42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p2er,add1,add2,     cy_r44,cy_r46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p33r,add1,add2,     cy_r48,cy_r50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p37r,add1,add2,     cy_r52,cy_r54,bjmodn52,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p3br,add1,add2,     cy_r56,cy_r58,bjmodn56,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			   #endif

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

			#else	// Scalar-double mode:

				l= j & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

			/*...set0 is slightly different from others:	*/
			   cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy_r00,bjmodn00   );
				cmplx_carry_norm_errcheck(a1p01r,a1p01i,cy_r01,bjmodn01, 1);
				cmplx_carry_norm_errcheck(a1p02r,a1p02i,cy_r02,bjmodn02, 2);
				cmplx_carry_norm_errcheck(a1p03r,a1p03i,cy_r03,bjmodn03, 3);
				cmplx_carry_norm_errcheck(a1p04r,a1p04i,cy_r04,bjmodn04, 4);
				cmplx_carry_norm_errcheck(a1p05r,a1p05i,cy_r05,bjmodn05, 5);
				cmplx_carry_norm_errcheck(a1p06r,a1p06i,cy_r06,bjmodn06, 6);
				cmplx_carry_norm_errcheck(a1p07r,a1p07i,cy_r07,bjmodn07, 7);
				cmplx_carry_norm_errcheck(a1p08r,a1p08i,cy_r08,bjmodn08, 8);
				cmplx_carry_norm_errcheck(a1p09r,a1p09i,cy_r09,bjmodn09, 9);
				cmplx_carry_norm_errcheck(a1p0ar,a1p0ai,cy_r10,bjmodn10,10);
				cmplx_carry_norm_errcheck(a1p0br,a1p0bi,cy_r11,bjmodn11,11);
				cmplx_carry_norm_errcheck(a1p0cr,a1p0ci,cy_r12,bjmodn12,12);
				cmplx_carry_norm_errcheck(a1p0dr,a1p0di,cy_r13,bjmodn13,13);
				cmplx_carry_norm_errcheck(a1p0er,a1p0ei,cy_r14,bjmodn14,14);
				cmplx_carry_norm_errcheck(a1p10r,a1p10i,cy_r15,bjmodn15,15);
				cmplx_carry_norm_errcheck(a1p11r,a1p11i,cy_r16,bjmodn16,16);
				cmplx_carry_norm_errcheck(a1p12r,a1p12i,cy_r17,bjmodn17,17);
				cmplx_carry_norm_errcheck(a1p13r,a1p13i,cy_r18,bjmodn18,18);
				cmplx_carry_norm_errcheck(a1p14r,a1p14i,cy_r19,bjmodn19,19);
				cmplx_carry_norm_errcheck(a1p15r,a1p15i,cy_r20,bjmodn20,20);
				cmplx_carry_norm_errcheck(a1p16r,a1p16i,cy_r21,bjmodn21,21);
				cmplx_carry_norm_errcheck(a1p17r,a1p17i,cy_r22,bjmodn22,22);
				cmplx_carry_norm_errcheck(a1p18r,a1p18i,cy_r23,bjmodn23,23);
				cmplx_carry_norm_errcheck(a1p19r,a1p19i,cy_r24,bjmodn24,24);
				cmplx_carry_norm_errcheck(a1p1ar,a1p1ai,cy_r25,bjmodn25,25);
				cmplx_carry_norm_errcheck(a1p1br,a1p1bi,cy_r26,bjmodn26,26);
				cmplx_carry_norm_errcheck(a1p1cr,a1p1ci,cy_r27,bjmodn27,27);
				cmplx_carry_norm_errcheck(a1p1dr,a1p1di,cy_r28,bjmodn28,28);
				cmplx_carry_norm_errcheck(a1p1er,a1p1ei,cy_r29,bjmodn29,29);
				cmplx_carry_norm_errcheck(a1p20r,a1p20i,cy_r30,bjmodn30,30);
				cmplx_carry_norm_errcheck(a1p21r,a1p21i,cy_r31,bjmodn31,31);
				cmplx_carry_norm_errcheck(a1p22r,a1p22i,cy_r32,bjmodn32,32);
				cmplx_carry_norm_errcheck(a1p23r,a1p23i,cy_r33,bjmodn33,33);
				cmplx_carry_norm_errcheck(a1p24r,a1p24i,cy_r34,bjmodn34,34);
				cmplx_carry_norm_errcheck(a1p25r,a1p25i,cy_r35,bjmodn35,35);
				cmplx_carry_norm_errcheck(a1p26r,a1p26i,cy_r36,bjmodn36,36);
				cmplx_carry_norm_errcheck(a1p27r,a1p27i,cy_r37,bjmodn37,37);
				cmplx_carry_norm_errcheck(a1p28r,a1p28i,cy_r38,bjmodn38,38);
				cmplx_carry_norm_errcheck(a1p29r,a1p29i,cy_r39,bjmodn39,39);
				cmplx_carry_norm_errcheck(a1p2ar,a1p2ai,cy_r40,bjmodn40,40);
				cmplx_carry_norm_errcheck(a1p2br,a1p2bi,cy_r41,bjmodn41,41);
				cmplx_carry_norm_errcheck(a1p2cr,a1p2ci,cy_r42,bjmodn42,42);
				cmplx_carry_norm_errcheck(a1p2dr,a1p2di,cy_r43,bjmodn43,43);
				cmplx_carry_norm_errcheck(a1p2er,a1p2ei,cy_r44,bjmodn44,44);
				cmplx_carry_norm_errcheck(a1p30r,a1p30i,cy_r45,bjmodn45,45);
				cmplx_carry_norm_errcheck(a1p31r,a1p31i,cy_r46,bjmodn46,46);
				cmplx_carry_norm_errcheck(a1p32r,a1p32i,cy_r47,bjmodn47,47);
				cmplx_carry_norm_errcheck(a1p33r,a1p33i,cy_r48,bjmodn48,48);
				cmplx_carry_norm_errcheck(a1p34r,a1p34i,cy_r49,bjmodn49,49);
				cmplx_carry_norm_errcheck(a1p35r,a1p35i,cy_r50,bjmodn50,50);
				cmplx_carry_norm_errcheck(a1p36r,a1p36i,cy_r51,bjmodn51,51);
				cmplx_carry_norm_errcheck(a1p37r,a1p37i,cy_r52,bjmodn52,52);
				cmplx_carry_norm_errcheck(a1p38r,a1p38i,cy_r53,bjmodn53,53);
				cmplx_carry_norm_errcheck(a1p39r,a1p39i,cy_r54,bjmodn54,54);
				cmplx_carry_norm_errcheck(a1p3ar,a1p3ai,cy_r55,bjmodn55,55);
				cmplx_carry_norm_errcheck(a1p3br,a1p3bi,cy_r56,bjmodn56,56);
				cmplx_carry_norm_errcheck(a1p3cr,a1p3ci,cy_r57,bjmodn57,57);
				cmplx_carry_norm_errcheck(a1p3dr,a1p3di,cy_r58,bjmodn58,58);
				cmplx_carry_norm_errcheck(a1p3er,a1p3ei,cy_r59,bjmodn59,59);

				i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

			#endif	// SIMD or not?

			}
			else	/* MODULUS_TYPE_FERMAT */
			{

			#ifdef USE_AVX

				// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.

				/* Get the needed Nth root of -1: */
				add1 = (double *)&rn0[0];
				add2 = (double *)&rn1[0];

				idx_offset = j;
				idx_incr = NDIVR;

				tmp = base_negacyclic_root;	tm2 = tmp+1;

			  #if HIACC
				// Hi-accuracy version needs 7 copies of each base root:
				l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				VEC_DBL_INIT(tmp+ 32,wt_re);	VEC_DBL_INIT(tm2+ 32,wt_im);
				VEC_DBL_INIT(tmp+ 40,wt_re);	VEC_DBL_INIT(tm2+ 40,wt_im);
				VEC_DBL_INIT(tmp+ 48,wt_re);	VEC_DBL_INIT(tm2+ 48,wt_im);
				VEC_DBL_INIT(tmp+ 56,wt_re);	VEC_DBL_INIT(tm2+ 56,wt_im);
				VEC_DBL_INIT(tmp+ 64,wt_re);	VEC_DBL_INIT(tm2+ 64,wt_im);
				VEC_DBL_INIT(tmp+ 72,wt_re);	VEC_DBL_INIT(tm2+ 72,wt_im);
				VEC_DBL_INIT(tmp+ 80,wt_re);	VEC_DBL_INIT(tm2+ 80,wt_im);
				VEC_DBL_INIT(tmp+ 88,wt_re);	VEC_DBL_INIT(tm2+ 88,wt_im);
				VEC_DBL_INIT(tmp+ 96,wt_re);	VEC_DBL_INIT(tm2+ 96,wt_im);
				VEC_DBL_INIT(tmp+104,wt_re);	VEC_DBL_INIT(tm2+104,wt_im);
				VEC_DBL_INIT(tmp+112,wt_re);	VEC_DBL_INIT(tm2+112,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				VEC_DBL_INIT(tmp+ 32,wt_re);	VEC_DBL_INIT(tm2+ 32,wt_im);
				VEC_DBL_INIT(tmp+ 40,wt_re);	VEC_DBL_INIT(tm2+ 40,wt_im);
				VEC_DBL_INIT(tmp+ 48,wt_re);	VEC_DBL_INIT(tm2+ 48,wt_im);
				VEC_DBL_INIT(tmp+ 56,wt_re);	VEC_DBL_INIT(tm2+ 56,wt_im);
				VEC_DBL_INIT(tmp+ 64,wt_re);	VEC_DBL_INIT(tm2+ 64,wt_im);
				VEC_DBL_INIT(tmp+ 72,wt_re);	VEC_DBL_INIT(tm2+ 72,wt_im);
				VEC_DBL_INIT(tmp+ 80,wt_re);	VEC_DBL_INIT(tm2+ 80,wt_im);
				VEC_DBL_INIT(tmp+ 88,wt_re);	VEC_DBL_INIT(tm2+ 88,wt_im);
				VEC_DBL_INIT(tmp+ 96,wt_re);	VEC_DBL_INIT(tm2+ 96,wt_im);
				VEC_DBL_INIT(tmp+104,wt_re);	VEC_DBL_INIT(tm2+104,wt_im);
				VEC_DBL_INIT(tmp+112,wt_re);	VEC_DBL_INIT(tm2+112,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				VEC_DBL_INIT(tmp+ 32,wt_re);	VEC_DBL_INIT(tm2+ 32,wt_im);
				VEC_DBL_INIT(tmp+ 40,wt_re);	VEC_DBL_INIT(tm2+ 40,wt_im);
				VEC_DBL_INIT(tmp+ 48,wt_re);	VEC_DBL_INIT(tm2+ 48,wt_im);
				VEC_DBL_INIT(tmp+ 56,wt_re);	VEC_DBL_INIT(tm2+ 56,wt_im);
				VEC_DBL_INIT(tmp+ 64,wt_re);	VEC_DBL_INIT(tm2+ 64,wt_im);
				VEC_DBL_INIT(tmp+ 72,wt_re);	VEC_DBL_INIT(tm2+ 72,wt_im);
				VEC_DBL_INIT(tmp+ 80,wt_re);	VEC_DBL_INIT(tm2+ 80,wt_im);
				VEC_DBL_INIT(tmp+ 88,wt_re);	VEC_DBL_INIT(tm2+ 88,wt_im);
				VEC_DBL_INIT(tmp+ 96,wt_re);	VEC_DBL_INIT(tm2+ 96,wt_im);
				VEC_DBL_INIT(tmp+104,wt_re);	VEC_DBL_INIT(tm2+104,wt_im);
				VEC_DBL_INIT(tmp+112,wt_re);	VEC_DBL_INIT(tm2+112,wt_im);
				tmp += 2;	tm2 += 2;
				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
				VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
				VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
				VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
				VEC_DBL_INIT(tmp+ 32,wt_re);	VEC_DBL_INIT(tm2+ 32,wt_im);
				VEC_DBL_INIT(tmp+ 40,wt_re);	VEC_DBL_INIT(tm2+ 40,wt_im);
				VEC_DBL_INIT(tmp+ 48,wt_re);	VEC_DBL_INIT(tm2+ 48,wt_im);
				VEC_DBL_INIT(tmp+ 56,wt_re);	VEC_DBL_INIT(tm2+ 56,wt_im);
				VEC_DBL_INIT(tmp+ 64,wt_re);	VEC_DBL_INIT(tm2+ 64,wt_im);
				VEC_DBL_INIT(tmp+ 72,wt_re);	VEC_DBL_INIT(tm2+ 72,wt_im);
				VEC_DBL_INIT(tmp+ 80,wt_re);	VEC_DBL_INIT(tm2+ 80,wt_im);
				VEC_DBL_INIT(tmp+ 88,wt_re);	VEC_DBL_INIT(tm2+ 88,wt_im);
				VEC_DBL_INIT(tmp+ 96,wt_re);	VEC_DBL_INIT(tm2+ 96,wt_im);
				VEC_DBL_INIT(tmp+104,wt_re);	VEC_DBL_INIT(tm2+104,wt_im);
				VEC_DBL_INIT(tmp+112,wt_re);	VEC_DBL_INIT(tm2+112,wt_im);

			  #else	// HIACC = false:

				// Get the needed quartet of Nth roots of -1: This is the same code as in the scalar
				// fermat_carry_norm_errcheck() macro, with the single index j replaced by the quartet j,j+2,j+4,j+6:
				l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

				l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
				dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
				rt  =rn1[k2].re;			it   =rn1[k2].im;
				wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
				VEC_DBL_INIT(tmp,wt_re);	++tmp;	VEC_DBL_INIT(tmp,wt_im);	++tmp;

				// The above need some inits to prepare for the AVX version of the Fermat-mod carry macro:
				SSE2_fermat_carry_init_loacc(base_negacyclic_root);

			  #endif

				// AVX-custom 4-way carry macro - each contains 4 of the RADIX stride-n/RADIX-separated carries
				// (processed independently in parallel), and steps through sequential-data indicies j,j+2,j+4,j+6:
			  #if HIACC
				// The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl)
				// which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
				// precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
				// uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls:

				tmp = base_negacyclic_root+  0;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p00r,tmp,0xf00,cy_r00,cy_i00,odd_radix,half_arr,sign_mask,icycle00,icycle01,icycle02,icycle03, jcycle00,kcycle00,lcycle00);	// *cycle0 index increments by +4 (mod odd_radix) between macro calls
				tmp = base_negacyclic_root+  8;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p04r,tmp,0xe40,cy_r04,cy_i04,odd_radix,half_arr,sign_mask,icycle04,icycle05,icycle06,icycle07, jcycle04,kcycle04,lcycle04);
				tmp = base_negacyclic_root+ 16;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p08r,tmp,0xd80,cy_r08,cy_i08,odd_radix,half_arr,sign_mask,icycle08,icycle09,icycle10,icycle11, jcycle08,kcycle08,lcycle08);
				tmp = base_negacyclic_root+ 24;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p0cr,tmp,0xcc0,cy_r12,cy_i12,odd_radix,half_arr,sign_mask,icycle12,icycle13,icycle14,icycle00, jcycle12,kcycle12,lcycle12);
				tmp = base_negacyclic_root+ 32;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p11r,tmp,0xc00,cy_r16,cy_i16,odd_radix,half_arr,sign_mask,icycle01,icycle02,icycle03,icycle04, jcycle01,kcycle01,lcycle01);
				tmp = base_negacyclic_root+ 40;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p15r,tmp,0xb40,cy_r20,cy_i20,odd_radix,half_arr,sign_mask,icycle05,icycle06,icycle07,icycle08, jcycle05,kcycle05,lcycle05);
				tmp = base_negacyclic_root+ 48;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p19r,tmp,0xa80,cy_r24,cy_i24,odd_radix,half_arr,sign_mask,icycle09,icycle10,icycle11,icycle12, jcycle09,kcycle09,lcycle09);
				tmp = base_negacyclic_root+ 56;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p1dr,tmp,0x9c0,cy_r28,cy_i28,odd_radix,half_arr,sign_mask,icycle13,icycle14,icycle00,icycle01, jcycle13,kcycle13,lcycle13);
				tmp = base_negacyclic_root+ 64;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p22r,tmp,0x900,cy_r32,cy_i32,odd_radix,half_arr,sign_mask,icycle02,icycle03,icycle04,icycle05, jcycle02,kcycle02,lcycle02);
				tmp = base_negacyclic_root+ 72;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p26r,tmp,0x840,cy_r36,cy_i36,odd_radix,half_arr,sign_mask,icycle06,icycle07,icycle08,icycle09, jcycle06,kcycle06,lcycle06);
				tmp = base_negacyclic_root+ 80;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p2ar,tmp,0x780,cy_r40,cy_i40,odd_radix,half_arr,sign_mask,icycle10,icycle11,icycle12,icycle13, jcycle10,kcycle10,lcycle10);
				tmp = base_negacyclic_root+ 88;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p2er,tmp,0x6c0,cy_r44,cy_i44,odd_radix,half_arr,sign_mask,icycle14,icycle00,icycle01,icycle02, jcycle14,kcycle14,lcycle14);
				tmp = base_negacyclic_root+ 96;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p33r,tmp,0x600,cy_r48,cy_i48,odd_radix,half_arr,sign_mask,icycle03,icycle04,icycle05,icycle06, jcycle03,kcycle03,lcycle03);
				tmp = base_negacyclic_root+104;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p37r,tmp,0x540,cy_r52,cy_i52,odd_radix,half_arr,sign_mask,icycle07,icycle08,icycle09,icycle10, jcycle07,kcycle07,lcycle07);
				tmp = base_negacyclic_root+112;	SSE2_fermat_carry_norm_errcheck_X4_hiacc(s1p3br,tmp,0x480,cy_r56,cy_i56,odd_radix,half_arr,sign_mask,icycle11,icycle12,icycle13,icycle14, jcycle11,kcycle11,lcycle11);

			  #else	// HIACC = false:

				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p00r,base_negacyclic_root,cy_r00,cy_i00,odd_radix,half_arr,sign_mask,                     icycle00,icycle01,icycle02,icycle03, jcycle00,kcycle00,lcycle00);	// *cycle index increments by +4 (mod odd_radix) between macro calls
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p04r,base_negacyclic_root,cy_r04,cy_i04,odd_radix,half_arr,sign_mask,                     icycle04,icycle05,icycle06,icycle07, jcycle04,kcycle04,lcycle04);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p08r,base_negacyclic_root,cy_r08,cy_i08,odd_radix,half_arr,sign_mask,                     icycle08,icycle09,icycle10,icycle11, jcycle08,kcycle08,lcycle08);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p0cr,base_negacyclic_root,cy_r12,cy_i12,odd_radix,half_arr,sign_mask,                     icycle12,icycle13,icycle14,icycle00, jcycle12,kcycle12,lcycle12);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p11r,base_negacyclic_root,cy_r16,cy_i16,odd_radix,half_arr,sign_mask,                     icycle01,icycle02,icycle03,icycle04, jcycle01,kcycle01,lcycle01);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p15r,base_negacyclic_root,cy_r20,cy_i20,odd_radix,half_arr,sign_mask,                     icycle05,icycle06,icycle07,icycle08, jcycle05,kcycle05,lcycle05);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p19r,base_negacyclic_root,cy_r24,cy_i24,odd_radix,half_arr,sign_mask,                     icycle09,icycle10,icycle11,icycle12, jcycle09,kcycle09,lcycle09);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p1dr,base_negacyclic_root,cy_r28,cy_i28,odd_radix,half_arr,sign_mask,                     icycle13,icycle14,icycle00,icycle01, jcycle13,kcycle13,lcycle13);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p22r,base_negacyclic_root,cy_r32,cy_i32,odd_radix,half_arr,sign_mask,                     icycle02,icycle03,icycle04,icycle05, jcycle02,kcycle02,lcycle02);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p26r,base_negacyclic_root,cy_r36,cy_i36,odd_radix,half_arr,sign_mask,                     icycle06,icycle07,icycle08,icycle09, jcycle06,kcycle06,lcycle06);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p2ar,base_negacyclic_root,cy_r40,cy_i40,odd_radix,half_arr,sign_mask,                     icycle10,icycle11,icycle12,icycle13, jcycle10,kcycle10,lcycle10);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p2er,base_negacyclic_root,cy_r44,cy_i44,odd_radix,half_arr,sign_mask,                     icycle14,icycle00,icycle01,icycle02, jcycle14,kcycle14,lcycle14);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p33r,base_negacyclic_root,cy_r48,cy_i48,odd_radix,half_arr,sign_mask,                     icycle03,icycle04,icycle05,icycle06, jcycle03,kcycle03,lcycle03);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p37r,base_negacyclic_root,cy_r52,cy_i52,odd_radix,half_arr,sign_mask,                     icycle07,icycle08,icycle09,icycle10, jcycle07,kcycle07,lcycle07);
				SSE2_fermat_carry_norm_errcheck_X4_loacc(s1p3br,base_negacyclic_root,cy_r56,cy_i56,odd_radix,half_arr,sign_mask,                     icycle11,icycle12,icycle13,icycle14, jcycle11,kcycle11,lcycle11);

			  #endif	// HIACC?

			#elif defined(USE_SSE2)

				// For a description of the data movement for Fermat-mod carries in SSE2 mode, see radix16_ditN_cy_dif1.c.

				/* Get the needed Nth root of -1: */
				add1 = (double *)&rn0[0];
				add2 = (double *)&rn1[0];

				idx_offset = j;
				idx_incr = NDIVR;

				SSE2_fermat_carry_norm_errcheck_X2(s1p00r,cy_r00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle00,jcycle00,icycle01,jcycle01);
				SSE2_fermat_carry_norm_errcheck_X2(s1p02r,cy_r04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle02,jcycle02,icycle03,jcycle03);
				SSE2_fermat_carry_norm_errcheck_X2(s1p04r,cy_r08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle04,jcycle04,icycle05,jcycle05);
				SSE2_fermat_carry_norm_errcheck_X2(s1p06r,cy_r12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle06,jcycle06,icycle07,jcycle07);
				SSE2_fermat_carry_norm_errcheck_X2(s1p08r,cy_r16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle08,jcycle08,icycle09,jcycle09);
				SSE2_fermat_carry_norm_errcheck_X2(s1p0ar,cy_r20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle10,jcycle10,icycle11,jcycle11);
				SSE2_fermat_carry_norm_errcheck_X2(s1p0cr,cy_r24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle12,jcycle12,icycle13,jcycle13);
				SSE2_fermat_carry_norm_errcheck_X2(s1p0er,cy_r28,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle14,jcycle14,icycle00,jcycle00);
				SSE2_fermat_carry_norm_errcheck_X2(s1p11r,cy_r32,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle01,jcycle01,icycle02,jcycle02);
				SSE2_fermat_carry_norm_errcheck_X2(s1p13r,cy_r36,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle03,jcycle03,icycle04,jcycle04);
				SSE2_fermat_carry_norm_errcheck_X2(s1p15r,cy_r40,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle05,jcycle05,icycle06,jcycle06);
				SSE2_fermat_carry_norm_errcheck_X2(s1p17r,cy_r44,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle07,jcycle07,icycle08,jcycle08);
				SSE2_fermat_carry_norm_errcheck_X2(s1p19r,cy_r48,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle09,jcycle09,icycle10,jcycle10);
				SSE2_fermat_carry_norm_errcheck_X2(s1p1br,cy_r52,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle11,jcycle11,icycle12,jcycle12);
				SSE2_fermat_carry_norm_errcheck_X2(s1p1dr,cy_r56,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle13,jcycle13,icycle14,jcycle14);
				SSE2_fermat_carry_norm_errcheck_X2(s1p20r,cy_i00,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle00,jcycle00,icycle01,jcycle01);
				SSE2_fermat_carry_norm_errcheck_X2(s1p22r,cy_i04,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle02,jcycle02,icycle03,jcycle03);
				SSE2_fermat_carry_norm_errcheck_X2(s1p24r,cy_i08,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle04,jcycle04,icycle05,jcycle05);
				SSE2_fermat_carry_norm_errcheck_X2(s1p26r,cy_i12,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle06,jcycle06,icycle07,jcycle07);
				SSE2_fermat_carry_norm_errcheck_X2(s1p28r,cy_i16,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle08,jcycle08,icycle09,jcycle09);
				SSE2_fermat_carry_norm_errcheck_X2(s1p2ar,cy_i20,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle10,jcycle10,icycle11,jcycle11);
				SSE2_fermat_carry_norm_errcheck_X2(s1p2cr,cy_i24,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle12,jcycle12,icycle13,jcycle13);
				SSE2_fermat_carry_norm_errcheck_X2(s1p2er,cy_i28,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle14,jcycle14,icycle00,jcycle00);
				SSE2_fermat_carry_norm_errcheck_X2(s1p31r,cy_i32,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle01,jcycle01,icycle02,jcycle02);
				SSE2_fermat_carry_norm_errcheck_X2(s1p33r,cy_i36,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle03,jcycle03,icycle04,jcycle04);
				SSE2_fermat_carry_norm_errcheck_X2(s1p35r,cy_i40,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle05,jcycle05,icycle06,jcycle06);
				SSE2_fermat_carry_norm_errcheck_X2(s1p37r,cy_i44,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle07,jcycle07,icycle08,jcycle08);
				SSE2_fermat_carry_norm_errcheck_X2(s1p39r,cy_i48,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle09,jcycle09,icycle10,jcycle10);
				SSE2_fermat_carry_norm_errcheck_X2(s1p3br,cy_i52,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle11,jcycle11,icycle12,jcycle12);
				SSE2_fermat_carry_norm_errcheck_X2(s1p3dr,cy_i56,NRT_BITS,NRTM1,idx_offset,idx_incr,odd_radix,half_arr,sign_mask,add1,add2,icycle13,jcycle13,icycle14,jcycle14);

			#else	// Scalar-double mode:

				ntmp = 0;
				fermat_carry_norm_errcheckB(a1p00r,a1p00i,cy_r00,cy_i00,icycle00,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p01r,a1p01i,cy_r01,cy_i01,icycle01,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p02r,a1p02i,cy_r02,cy_i02,icycle02,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p03r,a1p03i,cy_r03,cy_i03,icycle03,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p04r,a1p04i,cy_r04,cy_i04,icycle04,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p05r,a1p05i,cy_r05,cy_i05,icycle05,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p06r,a1p06i,cy_r06,cy_i06,icycle06,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p07r,a1p07i,cy_r07,cy_i07,icycle07,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p08r,a1p08i,cy_r08,cy_i08,icycle08,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p09r,a1p09i,cy_r09,cy_i09,icycle09,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p0ar,a1p0ai,cy_r10,cy_i10,icycle10,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p0br,a1p0bi,cy_r11,cy_i11,icycle11,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p0cr,a1p0ci,cy_r12,cy_i12,icycle12,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p0dr,a1p0di,cy_r13,cy_i13,icycle13,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p0er,a1p0ei,cy_r14,cy_i14,icycle14,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p10r,a1p10i,cy_r15,cy_i15,icycle00,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p11r,a1p11i,cy_r16,cy_i16,icycle01,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p12r,a1p12i,cy_r17,cy_i17,icycle02,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p13r,a1p13i,cy_r18,cy_i18,icycle03,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p14r,a1p14i,cy_r19,cy_i19,icycle04,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p15r,a1p15i,cy_r20,cy_i20,icycle05,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p16r,a1p16i,cy_r21,cy_i21,icycle06,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p17r,a1p17i,cy_r22,cy_i22,icycle07,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p18r,a1p18i,cy_r23,cy_i23,icycle08,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p19r,a1p19i,cy_r24,cy_i24,icycle09,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p1ar,a1p1ai,cy_r25,cy_i25,icycle10,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p1br,a1p1bi,cy_r26,cy_i26,icycle11,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p1cr,a1p1ci,cy_r27,cy_i27,icycle12,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p1dr,a1p1di,cy_r28,cy_i28,icycle13,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p1er,a1p1ei,cy_r29,cy_i29,icycle14,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p20r,a1p20i,cy_r30,cy_i30,icycle00,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p21r,a1p21i,cy_r31,cy_i31,icycle01,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p22r,a1p22i,cy_r32,cy_i32,icycle02,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p23r,a1p23i,cy_r33,cy_i33,icycle03,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p24r,a1p24i,cy_r34,cy_i34,icycle04,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p25r,a1p25i,cy_r35,cy_i35,icycle05,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p26r,a1p26i,cy_r36,cy_i36,icycle06,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p27r,a1p27i,cy_r37,cy_i37,icycle07,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p28r,a1p28i,cy_r38,cy_i38,icycle08,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p29r,a1p29i,cy_r39,cy_i39,icycle09,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p2ar,a1p2ai,cy_r40,cy_i40,icycle10,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p2br,a1p2bi,cy_r41,cy_i41,icycle11,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p2cr,a1p2ci,cy_r42,cy_i42,icycle12,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p2dr,a1p2di,cy_r43,cy_i43,icycle13,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p2er,a1p2ei,cy_r44,cy_i44,icycle14,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p30r,a1p30i,cy_r45,cy_i45,icycle00,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p31r,a1p31i,cy_r46,cy_i46,icycle01,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p32r,a1p32i,cy_r47,cy_i47,icycle02,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p33r,a1p33i,cy_r48,cy_i48,icycle03,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p34r,a1p34i,cy_r49,cy_i49,icycle04,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p35r,a1p35i,cy_r50,cy_i50,icycle05,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p36r,a1p36i,cy_r51,cy_i51,icycle06,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p37r,a1p37i,cy_r52,cy_i52,icycle07,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p38r,a1p38i,cy_r53,cy_i53,icycle08,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p39r,a1p39i,cy_r54,cy_i54,icycle09,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p3ar,a1p3ai,cy_r55,cy_i55,icycle10,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p3br,a1p3bi,cy_r56,cy_i56,icycle11,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p3cr,a1p3ci,cy_r57,cy_i57,icycle12,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p3dr,a1p3di,cy_r58,cy_i58,icycle13,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
				fermat_carry_norm_errcheckB(a1p3er,a1p3ei,cy_r59,cy_i59,icycle14,ntmp,NRTM1,NRT_BITS);

				icycle00 += wts_idx_incr;	/* Inside the loop use this, as it is faster than general-mod '% nwt' */
				icycle01 += wts_idx_incr;
				icycle02 += wts_idx_incr;
				icycle03 += wts_idx_incr;
				icycle04 += wts_idx_incr;
				icycle05 += wts_idx_incr;
				icycle06 += wts_idx_incr;
				icycle07 += wts_idx_incr;
				icycle08 += wts_idx_incr;
				icycle09 += wts_idx_incr;
				icycle10 += wts_idx_incr;
				icycle11 += wts_idx_incr;
				icycle12 += wts_idx_incr;
				icycle13 += wts_idx_incr;
				icycle14 += wts_idx_incr;
				icycle00 += ( (-(int)((uint32)icycle00 >> 31)) & nwt);
				icycle01 += ( (-(int)((uint32)icycle01 >> 31)) & nwt);
				icycle02 += ( (-(int)((uint32)icycle02 >> 31)) & nwt);
				icycle03 += ( (-(int)((uint32)icycle03 >> 31)) & nwt);
				icycle04 += ( (-(int)((uint32)icycle04 >> 31)) & nwt);
				icycle05 += ( (-(int)((uint32)icycle05 >> 31)) & nwt);
				icycle06 += ( (-(int)((uint32)icycle06 >> 31)) & nwt);
				icycle07 += ( (-(int)((uint32)icycle07 >> 31)) & nwt);
				icycle08 += ( (-(int)((uint32)icycle08 >> 31)) & nwt);
				icycle09 += ( (-(int)((uint32)icycle09 >> 31)) & nwt);
				icycle10 += ( (-(int)((uint32)icycle10 >> 31)) & nwt);
				icycle11 += ( (-(int)((uint32)icycle11 >> 31)) & nwt);
				icycle12 += ( (-(int)((uint32)icycle12 >> 31)) & nwt);
				icycle13 += ( (-(int)((uint32)icycle13 >> 31)) & nwt);
				icycle14 += ( (-(int)((uint32)icycle14 >> 31)) & nwt);

			#endif	/* #ifdef USE_SSE2 */

			// Here we nest AVX inside SSE2 since i/jcycle updates are for both, k/l for AVX-only:
			#ifdef USE_SSE2

				icycle00 += wts_idx_inc2;		icycle00 += ( (-(icycle00 < 0)) & nwt16);
				icycle01 += wts_idx_inc2;		icycle01 += ( (-(icycle01 < 0)) & nwt16);
				icycle02 += wts_idx_inc2;		icycle02 += ( (-(icycle02 < 0)) & nwt16);
				icycle03 += wts_idx_inc2;		icycle03 += ( (-(icycle03 < 0)) & nwt16);
				icycle04 += wts_idx_inc2;		icycle04 += ( (-(icycle04 < 0)) & nwt16);
				icycle05 += wts_idx_inc2;		icycle05 += ( (-(icycle05 < 0)) & nwt16);
				icycle06 += wts_idx_inc2;		icycle06 += ( (-(icycle06 < 0)) & nwt16);
				icycle07 += wts_idx_inc2;		icycle07 += ( (-(icycle07 < 0)) & nwt16);
				icycle08 += wts_idx_inc2;		icycle08 += ( (-(icycle08 < 0)) & nwt16);
				icycle09 += wts_idx_inc2;		icycle09 += ( (-(icycle09 < 0)) & nwt16);
				icycle10 += wts_idx_inc2;		icycle10 += ( (-(icycle10 < 0)) & nwt16);
				icycle11 += wts_idx_inc2;		icycle11 += ( (-(icycle11 < 0)) & nwt16);
				icycle12 += wts_idx_inc2;		icycle12 += ( (-(icycle12 < 0)) & nwt16);
				icycle13 += wts_idx_inc2;		icycle13 += ( (-(icycle13 < 0)) & nwt16);
				icycle14 += wts_idx_inc2;		icycle14 += ( (-(icycle14 < 0)) & nwt16);

				jcycle00 += wts_idx_inc2;		jcycle00 += ( (-(jcycle00 < 0)) & nwt16);
				jcycle01 += wts_idx_inc2;		jcycle01 += ( (-(jcycle01 < 0)) & nwt16);
				jcycle02 += wts_idx_inc2;		jcycle02 += ( (-(jcycle02 < 0)) & nwt16);
				jcycle03 += wts_idx_inc2;		jcycle03 += ( (-(jcycle03 < 0)) & nwt16);
				jcycle04 += wts_idx_inc2;		jcycle04 += ( (-(jcycle04 < 0)) & nwt16);
				jcycle05 += wts_idx_inc2;		jcycle05 += ( (-(jcycle05 < 0)) & nwt16);
				jcycle06 += wts_idx_inc2;		jcycle06 += ( (-(jcycle06 < 0)) & nwt16);
				jcycle07 += wts_idx_inc2;		jcycle07 += ( (-(jcycle07 < 0)) & nwt16);
				jcycle08 += wts_idx_inc2;		jcycle08 += ( (-(jcycle08 < 0)) & nwt16);
				jcycle09 += wts_idx_inc2;		jcycle09 += ( (-(jcycle09 < 0)) & nwt16);
				jcycle10 += wts_idx_inc2;		jcycle10 += ( (-(jcycle10 < 0)) & nwt16);
				jcycle11 += wts_idx_inc2;		jcycle11 += ( (-(jcycle11 < 0)) & nwt16);
				jcycle12 += wts_idx_inc2;		jcycle12 += ( (-(jcycle12 < 0)) & nwt16);
				jcycle13 += wts_idx_inc2;		jcycle13 += ( (-(jcycle13 < 0)) & nwt16);
				jcycle14 += wts_idx_inc2;		jcycle14 += ( (-(jcycle14 < 0)) & nwt16);

			  #ifdef USE_AVX
				kcycle00 += wts_idx_inc2;		kcycle00 += ( (-(kcycle00 < 0)) & nwt16);
				kcycle01 += wts_idx_inc2;		kcycle01 += ( (-(kcycle01 < 0)) & nwt16);
				kcycle02 += wts_idx_inc2;		kcycle02 += ( (-(kcycle02 < 0)) & nwt16);
				kcycle03 += wts_idx_inc2;		kcycle03 += ( (-(kcycle03 < 0)) & nwt16);
				kcycle04 += wts_idx_inc2;		kcycle04 += ( (-(kcycle04 < 0)) & nwt16);
				kcycle05 += wts_idx_inc2;		kcycle05 += ( (-(kcycle05 < 0)) & nwt16);
				kcycle06 += wts_idx_inc2;		kcycle06 += ( (-(kcycle06 < 0)) & nwt16);
				kcycle07 += wts_idx_inc2;		kcycle07 += ( (-(kcycle07 < 0)) & nwt16);
				kcycle08 += wts_idx_inc2;		kcycle08 += ( (-(kcycle08 < 0)) & nwt16);
				kcycle09 += wts_idx_inc2;		kcycle09 += ( (-(kcycle09 < 0)) & nwt16);
				kcycle10 += wts_idx_inc2;		kcycle10 += ( (-(kcycle10 < 0)) & nwt16);
				kcycle11 += wts_idx_inc2;		kcycle11 += ( (-(kcycle11 < 0)) & nwt16);
				kcycle12 += wts_idx_inc2;		kcycle12 += ( (-(kcycle12 < 0)) & nwt16);
				kcycle13 += wts_idx_inc2;		kcycle13 += ( (-(kcycle13 < 0)) & nwt16);
				kcycle14 += wts_idx_inc2;		kcycle14 += ( (-(kcycle14 < 0)) & nwt16);

				lcycle00 += wts_idx_inc2;		lcycle00 += ( (-(lcycle00 < 0)) & nwt16);
				lcycle01 += wts_idx_inc2;		lcycle01 += ( (-(lcycle01 < 0)) & nwt16);
				lcycle02 += wts_idx_inc2;		lcycle02 += ( (-(lcycle02 < 0)) & nwt16);
				lcycle03 += wts_idx_inc2;		lcycle03 += ( (-(lcycle03 < 0)) & nwt16);
				lcycle04 += wts_idx_inc2;		lcycle04 += ( (-(lcycle04 < 0)) & nwt16);
				lcycle05 += wts_idx_inc2;		lcycle05 += ( (-(lcycle05 < 0)) & nwt16);
				lcycle06 += wts_idx_inc2;		lcycle06 += ( (-(lcycle06 < 0)) & nwt16);
				lcycle07 += wts_idx_inc2;		lcycle07 += ( (-(lcycle07 < 0)) & nwt16);
				lcycle08 += wts_idx_inc2;		lcycle08 += ( (-(lcycle08 < 0)) & nwt16);
				lcycle09 += wts_idx_inc2;		lcycle09 += ( (-(lcycle09 < 0)) & nwt16);
				lcycle10 += wts_idx_inc2;		lcycle10 += ( (-(lcycle10 < 0)) & nwt16);
				lcycle11 += wts_idx_inc2;		lcycle11 += ( (-(lcycle11 < 0)) & nwt16);
				lcycle12 += wts_idx_inc2;		lcycle12 += ( (-(lcycle12 < 0)) & nwt16);
				lcycle13 += wts_idx_inc2;		lcycle13 += ( (-(lcycle13 < 0)) & nwt16);
				lcycle14 += wts_idx_inc2;		lcycle14 += ( (-(lcycle14 < 0)) & nwt16);
			  #endif
			#endif

			}	/* if(MODULUS_TYPE == ...) */

			/*...The radix-60 DIF pass is here:	*/
			#ifdef USE_SSE2
			  #if (OS_BITS == 32) || !defined(USE_64BIT_ASM_STYLE)	// In 64-bit mode, default is to use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just xmm0-7.
				/* NOTE that to permit us to re-use the s1p-data as both inputs and outputs of the radix-15 DIF DFT, must swap s1p[0123] -> s1p[0321] in output rows: */
			   #ifndef USE_LITERAL_BYTE_OFFSETS
				SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1,        s1p00r,s1p3br,s1p37r,s1p33r,s1p2er,s1p2ar,s1p26r,s1p22r,s1p1dr,s1p19r,s1p15r,s1p11r,s1p0cr,s1p08r,s1p04r, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r0a,r0b,r0c,r0d,r0e);
				SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1,        s1p30r,s1p2br,s1p27r,s1p23r,s1p1er,s1p1ar,s1p16r,s1p12r,s1p0dr,s1p09r,s1p05r,s1p01r,s1p3cr,s1p38r,s1p34r, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r1a,r1b,r1c,r1d,r1e);
				SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1,        s1p20r,s1p1br,s1p17r,s1p13r,s1p0er,s1p0ar,s1p06r,s1p02r,s1p3dr,s1p39r,s1p35r,s1p31r,s1p2cr,s1p28r,s1p24r, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r2a,r2b,r2c,r2d,r2e);
				SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1,        s1p10r,s1p0br,s1p07r,s1p03r,s1p3er,s1p3ar,s1p36r,s1p32r,s1p2dr,s1p29r,s1p25r,s1p21r,s1p1cr,s1p18r,s1p14r, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r3a,r3b,r3c,r3d,r3e);
			   #else	// Versions using base-address-plus-literal-byte-offsets:
				SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1, s1p00r, 0x000, 0x700, 0x680, 0x600, 0x580, 0x500, 0x480, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, r00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0);
				SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1, s1p00r, 0x5a0, 0x520, 0x4a0, 0x420, 0x3a0, 0x320, 0x2a0, 0x220, 0x1a0, 0x120, 0x0a0, 0x020, 0x720, 0x6a0, 0x620, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, r10,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0);
				SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1, s1p00r, 0x3c0, 0x340, 0x2c0, 0x240, 0x1c0, 0x140, 0x0c0, 0x040, 0x740, 0x6c0, 0x640, 0x5c0, 0x540, 0x4c0, 0x440, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, r20,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0);
				SSE2_RADIX_15_DIF(sse2_c3m1,sse2_cn1, s1p00r, 0x1e0, 0x160, 0x0e0, 0x060, 0x760, 0x6e0, 0x660, 0x5e0, 0x560, 0x4e0, 0x460, 0x3e0, 0x360, 0x2e0, 0x260, x00,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0, r30,0x000,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180,0x1a0,0x1c0);
			   #endif
			  #else
				SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1, s1p00r,s1p3br,s1p37r,s1p33r,s1p2er,s1p2ar,s1p26r,s1p22r,s1p1dr,s1p19r,s1p15r,s1p11r,s1p0cr,s1p08r,s1p04r, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r00,r01,r02,r03,r04,r05,r06,r07,r08,r09,r0a,r0b,r0c,r0d,r0e, s1p30r,s1p2br,s1p27r,s1p23r,s1p1er,s1p1ar,s1p16r,s1p12r,s1p0dr,s1p09r,s1p05r,s1p01r,s1p3cr,s1p38r,s1p34r, y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e, r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r1a,r1b,r1c,r1d,r1e);
				SSE2_RADIX_15_DIF_X2(sse2_c3m1,sse2_cn1, s1p20r,s1p1br,s1p17r,s1p13r,s1p0er,s1p0ar,s1p06r,s1p02r,s1p3dr,s1p39r,s1p35r,s1p31r,s1p2cr,s1p28r,s1p24r, x00,x01,x02,x03,x04,x05,x06,x07,x08,x09,x0a,x0b,x0c,x0d,x0e, r20,r21,r22,r23,r24,r25,r26,r27,r28,r29,r2a,r2b,r2c,r2d,r2e, s1p10r,s1p0br,s1p07r,s1p03r,s1p3er,s1p3ar,s1p36r,s1p32r,s1p2dr,s1p29r,s1p25r,s1p21r,s1p1cr,s1p18r,s1p14r, y00,y01,y02,y03,y04,y05,y06,y07,y08,y09,y0a,y0b,y0c,y0d,y0e, r30,r31,r32,r33,r34,r35,r36,r37,r38,r39,r3a,r3b,r3c,r3d,r3e);
			  #endif

			  #ifdef USE_AVX
			// Reorder blocks to yield sequentially increasing a-array offsets:
				add0 = &a[j1    ];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, 0x3c0)
				add3 = &a[j1+p04];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, 0x3c0)
				add1 = &a[j1+p08];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01, 0x3c0)
				add3 = &a[j1+p12];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, 0x3c0)
				add1 = &a[j1+p16];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0d, 0x3c0)
				add2 = &a[j1+p20];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, 0x3c0)
				add1 = &a[j1+p24];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0b, 0x3c0)
				add2 = &a[j1+p28];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, 0x3c0)
				add0 = &a[j1+p32];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09, 0x3c0)
				add2 = &a[j1+p36];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, 0x3c0)
				add0 = &a[j1+p40];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07, 0x3c0)
				add3 = &a[j1+p44];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, 0x3c0)
				add0 = &a[j1+p48];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05, 0x3c0)
				add3 = &a[j1+p52];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, 0x3c0)
				add1 = &a[j1+p56];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03, 0x3c0)

			  #else

			   #if !GCC_ASM_FULL_INLINE
				add0 = &a[j1    ];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00, 0x1e0)
				add1 = &a[j1+p08];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01, 0x1e0)
				add3 = &a[j1+p04];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02, 0x1e0)
				add1 = &a[j1+p56];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03, 0x1e0)
				add3 = &a[j1+p52];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04, 0x1e0)
				add0 = &a[j1+p48];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05, 0x1e0)
				add3 = &a[j1+p44];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06, 0x1e0)
				add0 = &a[j1+p40];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07, 0x1e0)
				add2 = &a[j1+p36];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08, 0x1e0)
				add0 = &a[j1+p32];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09, 0x1e0)
				add2 = &a[j1+p28];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0a, 0x1e0)
				add1 = &a[j1+p24];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0b, 0x1e0)
				add2 = &a[j1+p20];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0c, 0x1e0)
				add1 = &a[j1+p16];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0d, 0x1e0)
				add3 = &a[j1+p12];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0e, 0x1e0)
			   #else
				add0 = &a[j1    ];
				SSE2_RADIX60_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,r00);
			   #endif

			  #endif	// AVX or SSE2?

			#else	/* !USE_SSE2 */

				RADIX_15_DIF(a1p00r,a1p00i,a1p3br,a1p3bi,a1p37r,a1p37i,a1p33r,a1p33i,a1p2er,a1p2ei,a1p2ar,a1p2ai,a1p26r,a1p26i,a1p22r,a1p22i,a1p1dr,a1p1di,a1p19r,a1p19i,a1p15r,a1p15i,a1p11r,a1p11i,a1p0cr,a1p0ci,a1p08r,a1p08i,a1p04r,a1p04i,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,t0dr,t0di,t0er,t0ei,rt,it);
				RADIX_15_DIF(a1p30r,a1p30i,a1p2br,a1p2bi,a1p27r,a1p27i,a1p23r,a1p23i,a1p1er,a1p1ei,a1p1ar,a1p1ai,a1p16r,a1p16i,a1p12r,a1p12i,a1p0dr,a1p0di,a1p09r,a1p09i,a1p05r,a1p05i,a1p01r,a1p01i,a1p3cr,a1p3ci,a1p38r,a1p38i,a1p34r,a1p34i,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,t1dr,t1di,t1er,t1ei,rt,it);
				RADIX_15_DIF(a1p20r,a1p20i,a1p1br,a1p1bi,a1p17r,a1p17i,a1p13r,a1p13i,a1p0er,a1p0ei,a1p0ar,a1p0ai,a1p06r,a1p06i,a1p02r,a1p02i,a1p3dr,a1p3di,a1p39r,a1p39i,a1p35r,a1p35i,a1p31r,a1p31i,a1p2cr,a1p2ci,a1p28r,a1p28i,a1p24r,a1p24i,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,t2dr,t2di,t2er,t2ei,rt,it);
				RADIX_15_DIF(a1p10r,a1p10i,a1p0br,a1p0bi,a1p07r,a1p07i,a1p03r,a1p03i,a1p3er,a1p3ei,a1p3ar,a1p3ai,a1p36r,a1p36i,a1p32r,a1p32i,a1p2dr,a1p2di,a1p29r,a1p29i,a1p25r,a1p25i,a1p21r,a1p21i,a1p1cr,a1p1ci,a1p18r,a1p18i,a1p14r,a1p14i,x0r,x0i,x1r,x1i,x2r,x2i,x3r,x3i,x4r,x4i,x5r,x5i,x6r,x6i,x7r,x7i,x8r,x8i,x9r,x9i,xar,xai,xbr,xbi,xcr,xci,xdr,xdi,xer,xei,t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,t3dr,t3di,t3er,t3ei,rt,it);

				jt = j1    ; jp = j2    ;	RADIX_04_DIF(t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
				jt = j1+p08; jp = j2+p08;	RADIX_04_DIF(t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
				jt = j1+p04; jp = j2+p04;	RADIX_04_DIF(t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
				jt = j1+p56; jp = j2+p56;	RADIX_04_DIF(t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
				jt = j1+p52; jp = j2+p52;	RADIX_04_DIF(t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
				jt = j1+p48; jp = j2+p48;	RADIX_04_DIF(t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
				jt = j1+p44; jp = j2+p44;	RADIX_04_DIF(t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
				jt = j1+p40; jp = j2+p40;	RADIX_04_DIF(t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
				jt = j1+p36; jp = j2+p36;	RADIX_04_DIF(t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
				jt = j1+p32; jp = j2+p32;	RADIX_04_DIF(t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
				jt = j1+p28; jp = j2+p28;	RADIX_04_DIF(t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
				jt = j1+p24; jp = j2+p24;	RADIX_04_DIF(t0br,t0bi,t1br,t1bi,t2br,t2bi,t3br,t3bi,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
				jt = j1+p20; jp = j2+p20;	RADIX_04_DIF(t0cr,t0ci,t1cr,t1ci,t2cr,t2ci,t3cr,t3ci,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
				jt = j1+p16; jp = j2+p16;	RADIX_04_DIF(t0dr,t0di,t1dr,t1di,t2dr,t2di,t3dr,t3di,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
				jt = j1+p12; jp = j2+p12;	RADIX_04_DIF(t0er,t0ei,t1er,t1ei,t2er,t2ei,t3er,t3ei,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);

			#endif	// SSE2 or AVX?

			}	/* end for(j=_jstart; j < _jhi; j += 2) */

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
				jstart += nwt;
				jhi    += nwt;

				col += RADIX;
				co3 -= RADIX;
			}
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
		#ifdef USE_AVX

			thread_arg->cy_r00 = cy_r00->d0;
			thread_arg->cy_r01 = cy_r00->d1;
			thread_arg->cy_r02 = cy_r00->d2;
			thread_arg->cy_r03 = cy_r00->d3;
			thread_arg->cy_r04 = cy_r04->d0;
			thread_arg->cy_r05 = cy_r04->d1;
			thread_arg->cy_r06 = cy_r04->d2;
			thread_arg->cy_r07 = cy_r04->d3;
			thread_arg->cy_r08 = cy_r08->d0;
			thread_arg->cy_r09 = cy_r08->d1;
			thread_arg->cy_r10 = cy_r08->d2;
			thread_arg->cy_r11 = cy_r08->d3;
			thread_arg->cy_r12 = cy_r12->d0;
			thread_arg->cy_r13 = cy_r12->d1;
			thread_arg->cy_r14 = cy_r12->d2;
			thread_arg->cy_r15 = cy_r12->d3;
			thread_arg->cy_r16 = cy_r16->d0;
			thread_arg->cy_r17 = cy_r16->d1;
			thread_arg->cy_r18 = cy_r16->d2;
			thread_arg->cy_r19 = cy_r16->d3;
			thread_arg->cy_r20 = cy_r20->d0;
			thread_arg->cy_r21 = cy_r20->d1;
			thread_arg->cy_r22 = cy_r20->d2;
			thread_arg->cy_r23 = cy_r20->d3;
			thread_arg->cy_r24 = cy_r24->d0;
			thread_arg->cy_r25 = cy_r24->d1;
			thread_arg->cy_r26 = cy_r24->d2;
			thread_arg->cy_r27 = cy_r24->d3;
			thread_arg->cy_r28 = cy_r28->d0;
			thread_arg->cy_r29 = cy_r28->d1;
			thread_arg->cy_r30 = cy_r28->d2;
			thread_arg->cy_r31 = cy_r28->d3;
			thread_arg->cy_r32 = cy_r32->d0;
			thread_arg->cy_r33 = cy_r32->d1;
			thread_arg->cy_r34 = cy_r32->d2;
			thread_arg->cy_r35 = cy_r32->d3;
			thread_arg->cy_r36 = cy_r36->d0;
			thread_arg->cy_r37 = cy_r36->d1;
			thread_arg->cy_r38 = cy_r36->d2;
			thread_arg->cy_r39 = cy_r36->d3;
			thread_arg->cy_r40 = cy_r40->d0;
			thread_arg->cy_r41 = cy_r40->d1;
			thread_arg->cy_r42 = cy_r40->d2;
			thread_arg->cy_r43 = cy_r40->d3;
			thread_arg->cy_r44 = cy_r44->d0;
			thread_arg->cy_r45 = cy_r44->d1;
			thread_arg->cy_r46 = cy_r44->d2;
			thread_arg->cy_r47 = cy_r44->d3;
			thread_arg->cy_r48 = cy_r48->d0;
			thread_arg->cy_r49 = cy_r48->d1;
			thread_arg->cy_r50 = cy_r48->d2;
			thread_arg->cy_r51 = cy_r48->d3;
			thread_arg->cy_r52 = cy_r52->d0;
			thread_arg->cy_r53 = cy_r52->d1;
			thread_arg->cy_r54 = cy_r52->d2;
			thread_arg->cy_r55 = cy_r52->d3;
			thread_arg->cy_r56 = cy_r56->d0;
			thread_arg->cy_r57 = cy_r56->d1;
			thread_arg->cy_r58 = cy_r56->d2;
			thread_arg->cy_r59 = cy_r56->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );

		#elif defined(USE_SSE2)

			thread_arg->cy_r00 = cy_r00->d0;
			thread_arg->cy_r01 = cy_r00->d1;
			thread_arg->cy_r02 = cy_r02->d0;
			thread_arg->cy_r03 = cy_r02->d1;
			thread_arg->cy_r04 = cy_r04->d0;
			thread_arg->cy_r05 = cy_r04->d1;
			thread_arg->cy_r06 = cy_r06->d0;
			thread_arg->cy_r07 = cy_r06->d1;
			thread_arg->cy_r08 = cy_r08->d0;
			thread_arg->cy_r09 = cy_r08->d1;
			thread_arg->cy_r10 = cy_r10->d0;
			thread_arg->cy_r11 = cy_r10->d1;
			thread_arg->cy_r12 = cy_r12->d0;
			thread_arg->cy_r13 = cy_r12->d1;
			thread_arg->cy_r14 = cy_r14->d0;
			thread_arg->cy_r15 = cy_r14->d1;
			thread_arg->cy_r16 = cy_r16->d0;
			thread_arg->cy_r17 = cy_r16->d1;
			thread_arg->cy_r18 = cy_r18->d0;
			thread_arg->cy_r19 = cy_r18->d1;
			thread_arg->cy_r20 = cy_r20->d0;
			thread_arg->cy_r21 = cy_r20->d1;
			thread_arg->cy_r22 = cy_r22->d0;
			thread_arg->cy_r23 = cy_r22->d1;
			thread_arg->cy_r24 = cy_r24->d0;
			thread_arg->cy_r25 = cy_r24->d1;
			thread_arg->cy_r26 = cy_r26->d0;
			thread_arg->cy_r27 = cy_r26->d1;
			thread_arg->cy_r28 = cy_r28->d0;
			thread_arg->cy_r29 = cy_r28->d1;
			thread_arg->cy_r30 = cy_r30->d0;
			thread_arg->cy_r31 = cy_r30->d1;
			thread_arg->cy_r32 = cy_r32->d0;
			thread_arg->cy_r33 = cy_r32->d1;
			thread_arg->cy_r34 = cy_r34->d0;
			thread_arg->cy_r35 = cy_r34->d1;
			thread_arg->cy_r36 = cy_r36->d0;
			thread_arg->cy_r37 = cy_r36->d1;
			thread_arg->cy_r38 = cy_r38->d0;
			thread_arg->cy_r39 = cy_r38->d1;
			thread_arg->cy_r40 = cy_r40->d0;
			thread_arg->cy_r41 = cy_r40->d1;
			thread_arg->cy_r42 = cy_r42->d0;
			thread_arg->cy_r43 = cy_r42->d1;
			thread_arg->cy_r44 = cy_r44->d0;
			thread_arg->cy_r45 = cy_r44->d1;
			thread_arg->cy_r46 = cy_r46->d0;
			thread_arg->cy_r47 = cy_r46->d1;
			thread_arg->cy_r48 = cy_r48->d0;
			thread_arg->cy_r49 = cy_r48->d1;
			thread_arg->cy_r50 = cy_r50->d0;
			thread_arg->cy_r51 = cy_r50->d1;
			thread_arg->cy_r52 = cy_r52->d0;
			thread_arg->cy_r53 = cy_r52->d1;
			thread_arg->cy_r54 = cy_r54->d0;
			thread_arg->cy_r55 = cy_r54->d1;
			thread_arg->cy_r56 = cy_r56->d0;
			thread_arg->cy_r57 = cy_r56->d1;
			thread_arg->cy_r58 = cy_r58->d0;
			thread_arg->cy_r59 = cy_r58->d1;
			maxerr = MAX(max_err->d0,max_err->d1);

		#else

			thread_arg->cy_r00 = cy_r00;
			thread_arg->cy_r01 = cy_r01;
			thread_arg->cy_r02 = cy_r02;
			thread_arg->cy_r03 = cy_r03;
			thread_arg->cy_r04 = cy_r04;
			thread_arg->cy_r05 = cy_r05;
			thread_arg->cy_r06 = cy_r06;
			thread_arg->cy_r07 = cy_r07;
			thread_arg->cy_r08 = cy_r08;
			thread_arg->cy_r09 = cy_r09;
			thread_arg->cy_r10 = cy_r10;
			thread_arg->cy_r11 = cy_r11;
			thread_arg->cy_r12 = cy_r12;
			thread_arg->cy_r13 = cy_r13;
			thread_arg->cy_r14 = cy_r14;
			thread_arg->cy_r15 = cy_r15;
			thread_arg->cy_r16 = cy_r16;
			thread_arg->cy_r17 = cy_r17;
			thread_arg->cy_r18 = cy_r18;
			thread_arg->cy_r19 = cy_r19;
			thread_arg->cy_r20 = cy_r20;
			thread_arg->cy_r21 = cy_r21;
			thread_arg->cy_r22 = cy_r22;
			thread_arg->cy_r23 = cy_r23;
			thread_arg->cy_r24 = cy_r24;
			thread_arg->cy_r25 = cy_r25;
			thread_arg->cy_r26 = cy_r26;
			thread_arg->cy_r27 = cy_r27;
			thread_arg->cy_r28 = cy_r28;
			thread_arg->cy_r29 = cy_r29;
			thread_arg->cy_r30 = cy_r30;
			thread_arg->cy_r31 = cy_r31;
			thread_arg->cy_r32 = cy_r32;
			thread_arg->cy_r33 = cy_r33;
			thread_arg->cy_r34 = cy_r34;
			thread_arg->cy_r35 = cy_r35;
			thread_arg->cy_r36 = cy_r36;
			thread_arg->cy_r37 = cy_r37;
			thread_arg->cy_r38 = cy_r38;
			thread_arg->cy_r39 = cy_r39;
			thread_arg->cy_r40 = cy_r40;
			thread_arg->cy_r41 = cy_r41;
			thread_arg->cy_r42 = cy_r42;
			thread_arg->cy_r43 = cy_r43;
			thread_arg->cy_r44 = cy_r44;
			thread_arg->cy_r45 = cy_r45;
			thread_arg->cy_r46 = cy_r46;
			thread_arg->cy_r47 = cy_r47;
			thread_arg->cy_r48 = cy_r48;
			thread_arg->cy_r49 = cy_r49;
			thread_arg->cy_r50 = cy_r50;
			thread_arg->cy_r51 = cy_r51;
			thread_arg->cy_r52 = cy_r52;
			thread_arg->cy_r53 = cy_r53;
			thread_arg->cy_r54 = cy_r54;
			thread_arg->cy_r55 = cy_r55;
			thread_arg->cy_r56 = cy_r56;
			thread_arg->cy_r57 = cy_r57;
			thread_arg->cy_r58 = cy_r58;
			thread_arg->cy_r59 = cy_r59;

		#endif	// SSE2 or AVX?
		}
		else
		{
		#ifdef USE_AVX

			thread_arg->cy_r00 = cy_r00->d0;	thread_arg->cy_i00 = cy_i00->d0;
			thread_arg->cy_r01 = cy_r00->d1;	thread_arg->cy_i01 = cy_i00->d1;
			thread_arg->cy_r02 = cy_r00->d2;	thread_arg->cy_i02 = cy_i00->d2;
			thread_arg->cy_r03 = cy_r00->d3;	thread_arg->cy_i03 = cy_i00->d3;
			thread_arg->cy_r04 = cy_r04->d0;	thread_arg->cy_i04 = cy_i04->d0;
			thread_arg->cy_r05 = cy_r04->d1;	thread_arg->cy_i05 = cy_i04->d1;
			thread_arg->cy_r06 = cy_r04->d2;	thread_arg->cy_i06 = cy_i04->d2;
			thread_arg->cy_r07 = cy_r04->d3;	thread_arg->cy_i07 = cy_i04->d3;
			thread_arg->cy_r08 = cy_r08->d0;	thread_arg->cy_i08 = cy_i08->d0;
			thread_arg->cy_r09 = cy_r08->d1;	thread_arg->cy_i09 = cy_i08->d1;
			thread_arg->cy_r10 = cy_r08->d2;	thread_arg->cy_i10 = cy_i08->d2;
			thread_arg->cy_r11 = cy_r08->d3;	thread_arg->cy_i11 = cy_i08->d3;
			thread_arg->cy_r12 = cy_r12->d0;	thread_arg->cy_i12 = cy_i12->d0;
			thread_arg->cy_r13 = cy_r12->d1;	thread_arg->cy_i13 = cy_i12->d1;
			thread_arg->cy_r14 = cy_r12->d2;	thread_arg->cy_i14 = cy_i12->d2;
			thread_arg->cy_r15 = cy_r12->d3;	thread_arg->cy_i15 = cy_i12->d3;
			thread_arg->cy_r16 = cy_r16->d0;	thread_arg->cy_i16 = cy_i16->d0;
			thread_arg->cy_r17 = cy_r16->d1;	thread_arg->cy_i17 = cy_i16->d1;
			thread_arg->cy_r18 = cy_r16->d2;	thread_arg->cy_i18 = cy_i16->d2;
			thread_arg->cy_r19 = cy_r16->d3;	thread_arg->cy_i19 = cy_i16->d3;
			thread_arg->cy_r20 = cy_r20->d0;	thread_arg->cy_i20 = cy_i20->d0;
			thread_arg->cy_r21 = cy_r20->d1;	thread_arg->cy_i21 = cy_i20->d1;
			thread_arg->cy_r22 = cy_r20->d2;	thread_arg->cy_i22 = cy_i20->d2;
			thread_arg->cy_r23 = cy_r20->d3;	thread_arg->cy_i23 = cy_i20->d3;
			thread_arg->cy_r24 = cy_r24->d0;	thread_arg->cy_i24 = cy_i24->d0;
			thread_arg->cy_r25 = cy_r24->d1;	thread_arg->cy_i25 = cy_i24->d1;
			thread_arg->cy_r26 = cy_r24->d2;	thread_arg->cy_i26 = cy_i24->d2;
			thread_arg->cy_r27 = cy_r24->d3;	thread_arg->cy_i27 = cy_i24->d3;
			thread_arg->cy_r28 = cy_r28->d0;	thread_arg->cy_i28 = cy_i28->d0;
			thread_arg->cy_r29 = cy_r28->d1;	thread_arg->cy_i29 = cy_i28->d1;
			thread_arg->cy_r30 = cy_r28->d2;	thread_arg->cy_i30 = cy_i28->d2;
			thread_arg->cy_r31 = cy_r28->d3;	thread_arg->cy_i31 = cy_i28->d3;
			thread_arg->cy_r32 = cy_r32->d0;	thread_arg->cy_i32 = cy_i32->d0;
			thread_arg->cy_r33 = cy_r32->d1;	thread_arg->cy_i33 = cy_i32->d1;
			thread_arg->cy_r34 = cy_r32->d2;	thread_arg->cy_i34 = cy_i32->d2;
			thread_arg->cy_r35 = cy_r32->d3;	thread_arg->cy_i35 = cy_i32->d3;
			thread_arg->cy_r36 = cy_r36->d0;	thread_arg->cy_i36 = cy_i36->d0;
			thread_arg->cy_r37 = cy_r36->d1;	thread_arg->cy_i37 = cy_i36->d1;
			thread_arg->cy_r38 = cy_r36->d2;	thread_arg->cy_i38 = cy_i36->d2;
			thread_arg->cy_r39 = cy_r36->d3;	thread_arg->cy_i39 = cy_i36->d3;
			thread_arg->cy_r40 = cy_r40->d0;	thread_arg->cy_i40 = cy_i40->d0;
			thread_arg->cy_r41 = cy_r40->d1;	thread_arg->cy_i41 = cy_i40->d1;
			thread_arg->cy_r42 = cy_r40->d2;	thread_arg->cy_i42 = cy_i40->d2;
			thread_arg->cy_r43 = cy_r40->d3;	thread_arg->cy_i43 = cy_i40->d3;
			thread_arg->cy_r44 = cy_r44->d0;	thread_arg->cy_i44 = cy_i44->d0;
			thread_arg->cy_r45 = cy_r44->d1;	thread_arg->cy_i45 = cy_i44->d1;
			thread_arg->cy_r46 = cy_r44->d2;	thread_arg->cy_i46 = cy_i44->d2;
			thread_arg->cy_r47 = cy_r44->d3;	thread_arg->cy_i47 = cy_i44->d3;
			thread_arg->cy_r48 = cy_r48->d0;	thread_arg->cy_i48 = cy_i48->d0;
			thread_arg->cy_r49 = cy_r48->d1;	thread_arg->cy_i49 = cy_i48->d1;
			thread_arg->cy_r50 = cy_r48->d2;	thread_arg->cy_i50 = cy_i48->d2;
			thread_arg->cy_r51 = cy_r48->d3;	thread_arg->cy_i51 = cy_i48->d3;
			thread_arg->cy_r52 = cy_r52->d0;	thread_arg->cy_i52 = cy_i52->d0;
			thread_arg->cy_r53 = cy_r52->d1;	thread_arg->cy_i53 = cy_i52->d1;
			thread_arg->cy_r54 = cy_r52->d2;	thread_arg->cy_i54 = cy_i52->d2;
			thread_arg->cy_r55 = cy_r52->d3;	thread_arg->cy_i55 = cy_i52->d3;
			thread_arg->cy_r56 = cy_r56->d0;	thread_arg->cy_i56 = cy_i56->d0;
			thread_arg->cy_r57 = cy_r56->d1;	thread_arg->cy_i57 = cy_i56->d1;
			thread_arg->cy_r58 = cy_r56->d2;	thread_arg->cy_i58 = cy_i56->d2;
			thread_arg->cy_r59 = cy_r56->d3;	thread_arg->cy_i59 = cy_i56->d3;
			maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );

		#elif defined(USE_SSE2)

			thread_arg->cy_r00 = cy_r00->d0;	thread_arg->cy_i00 = cy_r00->d1;
			thread_arg->cy_r01 = cy_r02->d0;	thread_arg->cy_i01 = cy_r02->d1;
			thread_arg->cy_r02 = cy_r04->d0;	thread_arg->cy_i02 = cy_r04->d1;
			thread_arg->cy_r03 = cy_r06->d0;	thread_arg->cy_i03 = cy_r06->d1;
			thread_arg->cy_r04 = cy_r08->d0;	thread_arg->cy_i04 = cy_r08->d1;
			thread_arg->cy_r05 = cy_r10->d0;	thread_arg->cy_i05 = cy_r10->d1;
			thread_arg->cy_r06 = cy_r12->d0;	thread_arg->cy_i06 = cy_r12->d1;
			thread_arg->cy_r07 = cy_r14->d0;	thread_arg->cy_i07 = cy_r14->d1;
			thread_arg->cy_r08 = cy_r16->d0;	thread_arg->cy_i08 = cy_r16->d1;
			thread_arg->cy_r09 = cy_r18->d0;	thread_arg->cy_i09 = cy_r18->d1;
			thread_arg->cy_r10 = cy_r20->d0;	thread_arg->cy_i10 = cy_r20->d1;
			thread_arg->cy_r11 = cy_r22->d0;	thread_arg->cy_i11 = cy_r22->d1;
			thread_arg->cy_r12 = cy_r24->d0;	thread_arg->cy_i12 = cy_r24->d1;
			thread_arg->cy_r13 = cy_r26->d0;	thread_arg->cy_i13 = cy_r26->d1;
			thread_arg->cy_r14 = cy_r28->d0;	thread_arg->cy_i14 = cy_r28->d1;
			thread_arg->cy_r15 = cy_r30->d0;	thread_arg->cy_i15 = cy_r30->d1;
			thread_arg->cy_r16 = cy_r32->d0;	thread_arg->cy_i16 = cy_r32->d1;
			thread_arg->cy_r17 = cy_r34->d0;	thread_arg->cy_i17 = cy_r34->d1;
			thread_arg->cy_r18 = cy_r36->d0;	thread_arg->cy_i18 = cy_r36->d1;
			thread_arg->cy_r19 = cy_r38->d0;	thread_arg->cy_i19 = cy_r38->d1;
			thread_arg->cy_r20 = cy_r40->d0;	thread_arg->cy_i20 = cy_r40->d1;
			thread_arg->cy_r21 = cy_r42->d0;	thread_arg->cy_i21 = cy_r42->d1;
			thread_arg->cy_r22 = cy_r44->d0;	thread_arg->cy_i22 = cy_r44->d1;
			thread_arg->cy_r23 = cy_r46->d0;	thread_arg->cy_i23 = cy_r46->d1;
			thread_arg->cy_r24 = cy_r48->d0;	thread_arg->cy_i24 = cy_r48->d1;
			thread_arg->cy_r25 = cy_r50->d0;	thread_arg->cy_i25 = cy_r50->d1;
			thread_arg->cy_r26 = cy_r52->d0;	thread_arg->cy_i26 = cy_r52->d1;
			thread_arg->cy_r27 = cy_r54->d0;	thread_arg->cy_i27 = cy_r54->d1;
			thread_arg->cy_r28 = cy_r56->d0;	thread_arg->cy_i28 = cy_r56->d1;
			thread_arg->cy_r29 = cy_r58->d0;	thread_arg->cy_i29 = cy_r58->d1;
			thread_arg->cy_r30 = cy_i00->d0;	thread_arg->cy_i30 = cy_i00->d1;
			thread_arg->cy_r31 = cy_i02->d0;	thread_arg->cy_i31 = cy_i02->d1;
			thread_arg->cy_r32 = cy_i04->d0;	thread_arg->cy_i32 = cy_i04->d1;
			thread_arg->cy_r33 = cy_i06->d0;	thread_arg->cy_i33 = cy_i06->d1;
			thread_arg->cy_r34 = cy_i08->d0;	thread_arg->cy_i34 = cy_i08->d1;
			thread_arg->cy_r35 = cy_i10->d0;	thread_arg->cy_i35 = cy_i10->d1;
			thread_arg->cy_r36 = cy_i12->d0;	thread_arg->cy_i36 = cy_i12->d1;
			thread_arg->cy_r37 = cy_i14->d0;	thread_arg->cy_i37 = cy_i14->d1;
			thread_arg->cy_r38 = cy_i16->d0;	thread_arg->cy_i38 = cy_i16->d1;
			thread_arg->cy_r39 = cy_i18->d0;	thread_arg->cy_i39 = cy_i18->d1;
			thread_arg->cy_r40 = cy_i20->d0;	thread_arg->cy_i40 = cy_i20->d1;
			thread_arg->cy_r41 = cy_i22->d0;	thread_arg->cy_i41 = cy_i22->d1;
			thread_arg->cy_r42 = cy_i24->d0;	thread_arg->cy_i42 = cy_i24->d1;
			thread_arg->cy_r43 = cy_i26->d0;	thread_arg->cy_i43 = cy_i26->d1;
			thread_arg->cy_r44 = cy_i28->d0;	thread_arg->cy_i44 = cy_i28->d1;
			thread_arg->cy_r45 = cy_i30->d0;	thread_arg->cy_i45 = cy_i30->d1;
			thread_arg->cy_r46 = cy_i32->d0;	thread_arg->cy_i46 = cy_i32->d1;
			thread_arg->cy_r47 = cy_i34->d0;	thread_arg->cy_i47 = cy_i34->d1;
			thread_arg->cy_r48 = cy_i36->d0;	thread_arg->cy_i48 = cy_i36->d1;
			thread_arg->cy_r49 = cy_i38->d0;	thread_arg->cy_i49 = cy_i38->d1;
			thread_arg->cy_r50 = cy_i40->d0;	thread_arg->cy_i50 = cy_i40->d1;
			thread_arg->cy_r51 = cy_i42->d0;	thread_arg->cy_i51 = cy_i42->d1;
			thread_arg->cy_r52 = cy_i44->d0;	thread_arg->cy_i52 = cy_i44->d1;
			thread_arg->cy_r53 = cy_i46->d0;	thread_arg->cy_i53 = cy_i46->d1;
			thread_arg->cy_r54 = cy_i48->d0;	thread_arg->cy_i54 = cy_i48->d1;
			thread_arg->cy_r55 = cy_i50->d0;	thread_arg->cy_i55 = cy_i50->d1;
			thread_arg->cy_r56 = cy_i52->d0;	thread_arg->cy_i56 = cy_i52->d1;
			thread_arg->cy_r57 = cy_i54->d0;	thread_arg->cy_i57 = cy_i54->d1;
			thread_arg->cy_r58 = cy_i56->d0;	thread_arg->cy_i58 = cy_i56->d1;
			thread_arg->cy_r59 = cy_i58->d0;	thread_arg->cy_i59 = cy_i58->d1;
			maxerr = MAX(max_err->d0,max_err->d1);

		#else

			thread_arg->cy_r00 = cy_r00;		thread_arg->cy_i00 = cy_i00;
			thread_arg->cy_r01 = cy_r01;		thread_arg->cy_i01 = cy_i01;
			thread_arg->cy_r02 = cy_r02;		thread_arg->cy_i02 = cy_i02;
			thread_arg->cy_r03 = cy_r03;		thread_arg->cy_i03 = cy_i03;
			thread_arg->cy_r04 = cy_r04;		thread_arg->cy_i04 = cy_i04;
			thread_arg->cy_r05 = cy_r05;		thread_arg->cy_i05 = cy_i05;
			thread_arg->cy_r06 = cy_r06;		thread_arg->cy_i06 = cy_i06;
			thread_arg->cy_r07 = cy_r07;		thread_arg->cy_i07 = cy_i07;
			thread_arg->cy_r08 = cy_r08;		thread_arg->cy_i08 = cy_i08;
			thread_arg->cy_r09 = cy_r09;		thread_arg->cy_i09 = cy_i09;
			thread_arg->cy_r10 = cy_r10;		thread_arg->cy_i10 = cy_i10;
			thread_arg->cy_r11 = cy_r11;		thread_arg->cy_i11 = cy_i11;
			thread_arg->cy_r12 = cy_r12;		thread_arg->cy_i12 = cy_i12;
			thread_arg->cy_r13 = cy_r13;		thread_arg->cy_i13 = cy_i13;
			thread_arg->cy_r14 = cy_r14;		thread_arg->cy_i14 = cy_i14;
			thread_arg->cy_r15 = cy_r15;		thread_arg->cy_i15 = cy_i15;
			thread_arg->cy_r16 = cy_r16;		thread_arg->cy_i16 = cy_i16;
			thread_arg->cy_r17 = cy_r17;		thread_arg->cy_i17 = cy_i17;
			thread_arg->cy_r18 = cy_r18;		thread_arg->cy_i18 = cy_i18;
			thread_arg->cy_r19 = cy_r19;		thread_arg->cy_i19 = cy_i19;
			thread_arg->cy_r20 = cy_r20;		thread_arg->cy_i20 = cy_i20;
			thread_arg->cy_r21 = cy_r21;		thread_arg->cy_i21 = cy_i21;
			thread_arg->cy_r22 = cy_r22;		thread_arg->cy_i22 = cy_i22;
			thread_arg->cy_r23 = cy_r23;		thread_arg->cy_i23 = cy_i23;
			thread_arg->cy_r24 = cy_r24;		thread_arg->cy_i24 = cy_i24;
			thread_arg->cy_r25 = cy_r25;		thread_arg->cy_i25 = cy_i25;
			thread_arg->cy_r26 = cy_r26;		thread_arg->cy_i26 = cy_i26;
			thread_arg->cy_r27 = cy_r27;		thread_arg->cy_i27 = cy_i27;
			thread_arg->cy_r28 = cy_r28;		thread_arg->cy_i28 = cy_i28;
			thread_arg->cy_r29 = cy_r29;		thread_arg->cy_i29 = cy_i29;
			thread_arg->cy_r30 = cy_r30;		thread_arg->cy_i30 = cy_i30;
			thread_arg->cy_r31 = cy_r31;		thread_arg->cy_i31 = cy_i31;
			thread_arg->cy_r32 = cy_r32;		thread_arg->cy_i32 = cy_i32;
			thread_arg->cy_r33 = cy_r33;		thread_arg->cy_i33 = cy_i33;
			thread_arg->cy_r34 = cy_r34;		thread_arg->cy_i34 = cy_i34;
			thread_arg->cy_r35 = cy_r35;		thread_arg->cy_i35 = cy_i35;
			thread_arg->cy_r36 = cy_r36;		thread_arg->cy_i36 = cy_i36;
			thread_arg->cy_r37 = cy_r37;		thread_arg->cy_i37 = cy_i37;
			thread_arg->cy_r38 = cy_r38;		thread_arg->cy_i38 = cy_i38;
			thread_arg->cy_r39 = cy_r39;		thread_arg->cy_i39 = cy_i39;
			thread_arg->cy_r40 = cy_r40;		thread_arg->cy_i40 = cy_i40;
			thread_arg->cy_r41 = cy_r41;		thread_arg->cy_i41 = cy_i41;
			thread_arg->cy_r42 = cy_r42;		thread_arg->cy_i42 = cy_i42;
			thread_arg->cy_r43 = cy_r43;		thread_arg->cy_i43 = cy_i43;
			thread_arg->cy_r44 = cy_r44;		thread_arg->cy_i44 = cy_i44;
			thread_arg->cy_r45 = cy_r45;		thread_arg->cy_i45 = cy_i45;
			thread_arg->cy_r46 = cy_r46;		thread_arg->cy_i46 = cy_i46;
			thread_arg->cy_r47 = cy_r47;		thread_arg->cy_i47 = cy_i47;
			thread_arg->cy_r48 = cy_r48;		thread_arg->cy_i48 = cy_i48;
			thread_arg->cy_r49 = cy_r49;		thread_arg->cy_i49 = cy_i49;
			thread_arg->cy_r50 = cy_r50;		thread_arg->cy_i50 = cy_i50;
			thread_arg->cy_r51 = cy_r51;		thread_arg->cy_i51 = cy_i51;
			thread_arg->cy_r52 = cy_r52;		thread_arg->cy_i52 = cy_i52;
			thread_arg->cy_r53 = cy_r53;		thread_arg->cy_i53 = cy_i53;
			thread_arg->cy_r54 = cy_r54;		thread_arg->cy_i54 = cy_i54;
			thread_arg->cy_r55 = cy_r55;		thread_arg->cy_i55 = cy_i55;
			thread_arg->cy_r56 = cy_r56;		thread_arg->cy_i56 = cy_i56;
			thread_arg->cy_r57 = cy_r57;		thread_arg->cy_i57 = cy_i57;
			thread_arg->cy_r58 = cy_r58;		thread_arg->cy_i58 = cy_i58;
			thread_arg->cy_r59 = cy_r59;		thread_arg->cy_i59 = cy_i59;

		#endif

		}

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(thread_arg->maxerr < maxerr)
		{
			thread_arg->maxerr = maxerr;
		}
		return 0x0;
	}
#endif

