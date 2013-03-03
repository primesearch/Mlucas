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
#include "radix13.h"

	#undef DEBUG_SSE2
//	#define DEBUG_SSE2

#ifdef USE_SSE2

	#if defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
		#define GCC_ASM_FULL_INLINE  1	// 0 to use small-macros to assemble radix-52 DFTs, 1 to inline fuse macros as a few big blobs of asm (64-bit only)
	#else
		#undef GCC_ASM_FULL_INLINE
	#endif

	const int radix52_creals_in_local_store = 276;

  #ifdef USE_PTHREAD

	#ifndef USE_SSE2
		#error Pthreading only available in SSE2 mode!
	#endif

	// Use non-pooled simple spawn/rejoin thread-team model
	#include "threadpool.h"

	struct cy_thread_data_t{
	// int data - if needed, pad to yield an even number of these:
		int tid;
		int ndivr;
	
		int khi;
		int i;
		int jstart;
		int jhi;
		int col;
		int co2;
		int co3;
		int sw;
		int nwt;

	// double data:
		double maxerr;
		double scale;

	// pointer data:
		double *arrdat;			/* Main data array */
		double *wt0;
		double *wt1;
		int *si;
		struct complex *s1p00r;
		struct complex *half_arr;

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
		/* carries: */
		double cy00;
		double cy01;
		double cy02;
		double cy03;
		double cy04;
		double cy05;
		double cy06;
		double cy07;
		double cy08;
		double cy09;
		double cy10;
		double cy11;
		double cy12;
		double cy13;
		double cy14;
		double cy15;
		double cy16;
		double cy17;
		double cy18;
		double cy19;
		double cy20;
		double cy21;
		double cy22;
		double cy23;
		double cy24;
		double cy25;
		double cy26;
		double cy27;
		double cy28;
		double cy29;
		double cy30;
		double cy31;
		double cy32;
		double cy33;
		double cy34;
		double cy35;
		double cy36;
		double cy37;
		double cy38;
		double cy39;
		double cy40;
		double cy41;
		double cy42;
		double cy43;
		double cy44;
		double cy45;
		double cy46;
		double cy47;
		double cy48;
		double cy49;
		double cy50;
		double cy51;
	};

  #endif

//	#define ERR_CHECK_ALL	/* #define this to do ROE checking of all convolution outputs, rather than just every 36th one */

	#include "sse2_macro.h"

	#if defined(COMPILER_TYPE_MSVC)
		/*
		__cc pointer offsets:
		-0x010 = 2.0
		0x000 =  DC1
		0x010 =  DC3
		0x020 =  DC4
		0x030 =  DS1
		0x040 =  DS2
		0x050 =  DS3
		0x060 =  DS4
		0x070 =  DS5
		0x080 = DC23
		0x090 = DC54
		0x0a0 = DC65
		0x0b0 = DS63
		0x0c0 = DS74
		0x0d0 = DS85
		0x0e0 = DS93
		0x0f0 = DSa4
		0x100 = DSb5
		*/
		/*...Radix-13 DFT: Inputs in memory locations __I0 + [__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA,__iB,__iC],\
		where r0 is a memory address and the i's are literal [byte] offsets. Outputs similarly go into memory locations\
		__O0 + [__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA,__oB,__oC], assumed disjoint with inputs:\

		Total SSE opcounts:
				32-bit:		148 MOVAPS, 164 ADD/SUBPD, 76 MULPD
				64-bit:		124 MOVAPS, 164 ADD/SUBPD, 76 MULPD
		*/\
		#define SSE2_RADIX_13_DFT(__cc,\
			__I0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8,__i9,__iA,__iB,__iC,\
			__O0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8,__o9,__oA,__oB,__oC)\
		{\
		/***************/\
		/* REAL PARTS: */\
		/***************/\
			__asm	mov	eax, __I0	\
			__asm	mov	ebx, __cc	\
			__asm	mov	ecx, __O0	\
			__asm	movaps	xmm7,[eax+__i6]		/* xr7 = __A6r */\
			__asm	movaps	xmm5,[eax+__i5]		/* xr5 = __A5r */\
			__asm	movaps	xmm4,[eax+__i4]		/* xr4 = __A4r */\
			__asm	movaps	xmm6,[eax+__i3]		/* xr6 = __A3r */\
			__asm	movaps	xmm3,[eax+__i2]		/* xr3 = __A2r */\
			__asm	movaps	xmm1,[eax+__i1]		/* xr1 = __A1r */\
		/* xr-terms need 8 registers for each side: */\
			__asm	movaps	xmm0,[ebx-0x10]		/* 2.0 */\
			__asm	addpd	xmm7,[eax+__i7]		/* xr7 += __A7r */\
			__asm	addpd	xmm5,[eax+__i8]		/* xr5 += __A8r */\
			__asm	addpd	xmm4,[eax+__i9]		/* xr4 += __A9r */\
			__asm	addpd	xmm6,[eax+__iA]		/* xr6 += __Aar */\
			__asm	addpd	xmm3,[eax+__iB]		/* xr3 += __Abr */\
			__asm	addpd	xmm1,[eax+__iC]		/* xr1 += __Acr */\
			__asm	subpd	xmm1,xmm5			__asm	mulpd	xmm5,xmm0	/* xr1 -= xr5;	xr5 *= 2.0 */\
			__asm	subpd	xmm3,xmm6			__asm	mulpd	xmm6,xmm0	/* xr3 -= xr6;	xr6 *= 2.0 */\
			__asm	subpd	xmm4,xmm7			__asm	mulpd	xmm7,xmm0	/* xr4 -= xr7;	xr7 *= 2.0 */\
			__asm	addpd	xmm5,xmm1			/* xr5 += xr1 */\
			__asm	addpd	xmm6,xmm3			/* xr6 += xr3 */\
			__asm	addpd	xmm7,xmm4			/* xr7 += xr4 */\
			__asm	movaps	xmm2,xmm5			/* xr2  = xr5 */\
			__asm	addpd	xmm2,xmm6			/* xr2 += xr6 */\
			__asm	addpd	xmm2,xmm7			/* xr2 += xr7 */\
			__asm	movaps	xmm0,[eax]			/* xr0 = __A0r */\
			__asm	addpd	xmm0,xmm2			/* xr0 += xr2 */\
			__asm	mulpd	xmm2,[ebx]			/* xr2 *= DC1 */\
			__asm	movaps	[ecx],xmm0			/* __B0r = xr0 */\
			__asm	addpd	xmm0,xmm2			/* xr0 += xr2 */\
			__asm	movaps	xmm2,[ebx+0x010]	/* xr2 = DC3 */\
			__asm	mulpd	xmm6,xmm2			/* xr6 *= xr2 */\
			__asm	mulpd	xmm5,xmm2			/* xr5 *= xr2 */\
			__asm	mulpd	xmm7,xmm2			/* xr7 *= xr2 */\
			__asm	movaps	[ecx+__o1],xmm6		/* __B1r = xr6 */\
			__asm	movaps	[ecx+__o2],xmm5		/* __B2r = xr5 */\
			__asm	movaps	[ecx+__o3],xmm7		/* __B3r = xr7 */\
			__asm	movaps	xmm2,[ebx+0x080]	/* xr2 = DC23 */\
			__asm	mulpd	xmm5,xmm2			/* xr5 *= xr2 */\
			__asm	mulpd	xmm7,xmm2			/* xr7 *= xr2 */\
			__asm	mulpd	xmm6,xmm2			/* xr6 *= xr2 */\
			__asm	addpd	xmm5,xmm0			/* xr5 += xr0 */\
			__asm	addpd	xmm7,xmm0			/* xr7 += xr0 */\
			__asm	addpd	xmm6,xmm0			/* xr6 += xr0 */\
			__asm	addpd	xmm5,[ecx+__o1]		/* xr5 += __B1r */\
			__asm	addpd	xmm7,[ecx+__o2]		/* xr7 += __B2r */\
			__asm	addpd	xmm6,[ecx+__o3]		/* xr6 += __B3r */\
		/* Rearrange the 3x3 multiply block to ease in-place computation: */\
			__asm	movaps	xmm2,[ebx+0x0a0]	/* xr2 = DC65 */\
			__asm	movaps	xmm0,[ebx+0x090]	/* xr0 = DC54 */\
			__asm	movaps	[ecx+__o1],xmm4		/* __B1r = xr4 */\
			__asm	movaps	[ecx+__o2],xmm1		/* __B2r = xr1 */\
			__asm	movaps	[ecx+__o3],xmm3		/* __B3r = xr3 */\
			__asm	mulpd	xmm4,xmm2			/* xr4 *= DC65 */\
			__asm	mulpd	xmm1,xmm2			/* xr1 *= DC65 */\
			__asm	mulpd	xmm3,xmm2			/* xr3 *= DC65 */\
			__asm	addpd	xmm4,[ecx+__o3]		/* xr4 += __B3r */\
			__asm	subpd	xmm1,[ecx+__o1]		/* xr1 -= __B1r */\
			__asm	addpd	xmm3,[ecx+__o2]		/* xr3 += __B2r */\
			__asm	mulpd	xmm4,xmm0			/* xr4 *= DC54 */\
			__asm	mulpd	xmm1,xmm0			/* xr1 *= DC54 */\
			__asm	mulpd	xmm3,xmm0			/* xr3 *= DC54 */\
			__asm	movaps	xmm2,[ebx+0x020]	/* xr2 = DC4 */\
			__asm	addpd	xmm4,[ecx+__o2]		/* xr4 += __B2r */\
			__asm	subpd	xmm1,[ecx+__o3]		/* xr1 -= __B3r */\
			__asm	subpd	xmm3,[ecx+__o1]		/* xr3 -= __B1r */\
			__asm	mulpd	xmm4,xmm2			/* xr4 *= DC4 */\
			__asm	mulpd	xmm1,xmm2			/* xr1 *= DC4 */\
			__asm	mulpd	xmm3,xmm2			/* xr3 *= DC4 */\
			__asm	movaps	xmm0,[ebx-0x010]	/* 2.0 */\
		/* Spill into destination outputs: */\
			__asm	subpd	xmm6,xmm4			__asm	mulpd	xmm4,xmm0	/* xr6 -= xr1;	xr5 *= 2.0 */\
			__asm	subpd	xmm7,xmm1			__asm	mulpd	xmm1,xmm0	/* xr7 -= xr3;	xr6 *= 2.0 */\
			__asm	subpd	xmm5,xmm3			__asm	mulpd	xmm3,xmm0	/* xr5 -= xr4;	xr7 *= 2.0 */\
			__asm	movaps	[ecx+__o3],xmm7		/* __B3r = xr7 */\
			__asm	movaps	[ecx+__o8],xmm5		/* __B8r = xr5 */\
			__asm	addpd	xmm4,xmm6			/* xr4 += xr6 */\
			__asm	addpd	xmm1,xmm7			/* xr7 += xr1 */\
			__asm	addpd	xmm3,xmm5			/* xr5 += xr3 */\
			__asm	movaps	[ecx+__o4],xmm6		/* __B4r = xr6 */\
			__asm	movaps	[ecx+__o6],xmm4		/* __B6r = xr4 */\
			__asm	movaps	[ecx+__o2],xmm1		/* __B2r = xr7 */\
			__asm	movaps	[ecx+__o1],xmm3		/* __B1r = xr5 */\
		/* yi-terms: */\
			__asm	add	eax, 0x10	\
			__asm	movaps	xmm1,[eax+__i1]		/* xr1 = __A1i */\
			__asm	movaps	xmm2,[eax+__i2]		/* xr2 = __A2i */\
			__asm	movaps	xmm5,[eax+__i3]		/* xr5 = __A3i */\
			__asm	movaps	xmm3,[eax+__i4]		/* xr3 = __A4i */\
			__asm	movaps	xmm4,[eax+__i5]		/* xr4 = __A5i */\
			__asm	movaps	xmm6,[eax+__i6]		/* xr6 = __A6i */\
			__asm	subpd	xmm1,[eax+__iC]		/* xr1 -= __Aci */\
			__asm	subpd	xmm2,[eax+__iB]		/* xr2 -= __Abi */\
			__asm	subpd	xmm5,[eax+__iA]		/* xr5 -= __Aai */\
			__asm	subpd	xmm3,[eax+__i9]		/* xr3 -= __A9i */\
			__asm	subpd	xmm4,[eax+__i8]		/* xr4 -= __A8i */\
			__asm	subpd	xmm6,[eax+__i7]		/* xr6 -= __A7i */\
			__asm	movaps	xmm7,xmm1			/* xr7 = xr1 */\
			__asm	movaps	xmm0,xmm2			/* xr0 = xr2 */\
			__asm	subpd	xmm7,xmm3			/* xr7 -= xr3 */\
			__asm	addpd	xmm0,xmm4			/* xr0 += xr4 */\
			__asm	addpd	xmm7,xmm5			/* xr7 += xr5 */\
			__asm	addpd	xmm0,xmm6			/* xr0 += xr6 */\
			__asm	addpd	xmm1,xmm3			/* xr1 += xr3 */\
			__asm	addpd	xmm5,xmm3			/* xr5 += xr3 */\
			__asm	subpd	xmm2,xmm4			/* xr2 -= xr4 */\
			__asm	subpd	xmm6,xmm4			/* xr6 -= xr4 */\
			__asm	movaps	xmm4,xmm0			/* xr4 = xr0 */\
			__asm	movaps	xmm3,xmm7			/* xr3 = xr7 */\
			__asm	mulpd	xmm0,[ebx+0x040]	/* xr0 *= DS2 */\
			__asm	mulpd	xmm7,[ebx+0x040]	/* xr7 *= DS2 */\
			__asm	mulpd	xmm4,[ebx+0x030]	/* xr4 *= DS1 */\
			__asm	mulpd	xmm3,[ebx+0x030]	/* xr3 *= DS1 */\
			__asm	subpd	xmm4,xmm7			/* xr4 -= xr7 */\
			__asm	addpd	xmm3,xmm0			/* xr3 += xr0 */\
			__asm	movaps  [ecx+__oC],xmm4		/* __Bcr = xr4; tmp-store in Bcr */\
			__asm	movaps	xmm0,xmm1			/* xr0 = xr1 */\
			__asm	movaps	xmm4,xmm5			/* xr4 = xr5 */\
			__asm	addpd	xmm0,xmm2			/* xr0 += xr2 */\
			__asm	addpd	xmm4,xmm6			/* xr4 += xr6 */\
		/*\
		xr7 = DS3*xr0-DS6*xr4;\
		xr0 = DS3*xr4-DS9*xr0;\
		*/\
			__asm	movaps  [ecx+__oB],xmm0		/* __Bbr = xr0; tmp-store in Bbr */\
			__asm	movaps	xmm7,xmm4			/* xr7 = xr4 */\
			__asm	mulpd	xmm4,[ebx+0x0b0]	/* xr4 *= DS63 */\
			__asm	mulpd	xmm0,[ebx+0x0e0]	/* xr0 *= DS93 */\
			__asm	addpd	xmm4,[ecx+__oB]		/* xr4 += __Bbr */\
			__asm	addpd	xmm0,xmm7			/* xr0 += xr7 */\
			__asm	mulpd	xmm4,[ebx+0x050]	/* xr4 *= DS3 */\
			__asm	mulpd	xmm0,[ebx+0x050]	/* xr0 *= DS3 */\
		/*\
		xmm7,xmm4+DS4*xmm2-DS7*xr6;\
		xmm2,xmm0+DS4*xmm6-DSa*xr2;\
		*/\
			__asm	movaps  [ecx+__oB],xmm2		/* __Bbr = xr2 */\
			__asm	movaps	xmm7,xmm6			/* xr7 = xr6 */\
			__asm	mulpd	xmm6,[ebx+0x0c0]	/* xr6 *= DS74 */\
			__asm	mulpd	xmm2,[ebx+0x0f0]	/* xr2 *= DSa4 */\
			__asm	addpd	xmm6,[ecx+__oB]		/* xr6 += __Bbr */\
			__asm	addpd	xmm2,xmm7			/* xr2 += xr7 */\
			__asm	mulpd	xmm6,[ebx+0x060]	/* xr6 *= DS4 */\
			__asm	mulpd	xmm2,[ebx+0x060]	/* xr2 *= DS4 */\
			__asm	addpd	xmm6,xmm4			/* xr6 += xr4 */\
			__asm	addpd	xmm2,xmm0			/* xr2 += xr0 */\
		/*\
		xmm7,xmm4+DS5*xmm1-DS8*xr5;\
		xmm0,xmm0+DS5*xmm5-DSb*xr1;\
		*/\
			__asm	movaps  [ecx+__oB],xmm1		/* __Bbr = xr1 */\
			__asm	movaps	xmm7,xmm5			/* xr7 = xr5 */\
			__asm	mulpd	xmm5,[ebx+0x0d0]	/* xr5 *= DS85 */\
			__asm	mulpd	xmm1,[ebx+0x100]	/* xr1 *= DSb5 */\
			__asm	addpd	xmm5,[ecx+__oB]		/* xr5 += __Bbr */\
			__asm	addpd	xmm1,xmm7			/* xr1 += xr7 */\
			__asm	mulpd	xmm5,[ebx+0x070]	/* xr5 *= DS5 */\
			__asm	mulpd	xmm1,[ebx+0x070]	/* xr1 *= DS5 */\
			__asm	addpd	xmm5,xmm4			/* xr5 += xr4 */\
			__asm	addpd	xmm1,xmm0			/* xr1 += xr0 */\
			__asm	movaps  xmm4,[ecx+__oC]		/* xr4 = __Bcr */\
			__asm	movaps	xmm7,xmm3			/* xr7  = xr3 */\
			__asm	movaps	xmm0,xmm4			/* xr0  = xr4 */\
			__asm	addpd	xmm7,xmm6			/* xr7 += xr6 */\
			__asm	addpd	xmm0,xmm5			/* xr0 += xr5 */\
			__asm	subpd	xmm6,xmm3			/* xr6 -= xr3 */\
			__asm	subpd	xmm5,xmm4			/* xr5 -= xr4 */\
			__asm	addpd	xmm6,xmm2			/* xr6 += xr2 */\
			__asm	addpd	xmm5,xmm1			/* xr5 += xr1 */\
			__asm	addpd	xmm2,xmm3			/* xr2 += xr3 */\
			__asm	addpd	xmm4,xmm1			/* xr4 += xr1 */\
		/* Combine xmm and yi-terms to get real parts of outputs: */\
			/* __B6r -= xr7; xr7 *= xr1; xr7 += __B6r; __B7r = xr7 */\
			__asm movaps xmm3,[ecx+__o6]	\
			__asm movaps xmm1,xmm3			\
			__asm subpd xmm3,xmm7			\
			__asm addpd xmm1,xmm7			\
			__asm movaps [ecx+__o6],xmm3	\
			__asm movaps [ecx+__o7],xmm1	\
			/* __B8r -= xr6;	xr6 *= xr1; xr6 += __B8r;   __B5r = xr6 */\
			__asm movaps xmm3,[ecx+__o8]	\
			__asm movaps xmm1,xmm3			\
			__asm subpd xmm3,xmm6			\
			__asm addpd xmm1,xmm6			\
			__asm movaps [ecx+__o8],xmm3	\
			__asm movaps [ecx+__o5],xmm1	\
			/* __B2r -= xr2;	xr2 *= xr1; xr2 += __B2r;   __Bbr = xr2 */\
			__asm movaps xmm3,[ecx+__o2]	\
			__asm movaps xmm1,xmm3			\
			__asm subpd xmm3,xmm2			\
			__asm addpd xmm1,xmm2			\
			__asm movaps [ecx+__o2],xmm3	\
			__asm movaps [ecx+__oB],xmm1	\
			/* __B3r -= xr0;	xr0 *= xr1; xr0 += __B3r;   __Bar = xr0 */\
			__asm movaps xmm3,[ecx+__o3]	\
			__asm movaps xmm1,xmm3			\
			__asm subpd xmm3,xmm0			\
			__asm addpd xmm1,xmm0			\
			__asm movaps [ecx+__o3],xmm3	\
			__asm movaps [ecx+__oA],xmm1	\
			/* __B4r -= xr5;	xr5 *= xr1; xr5 += __B4r;   __B9r = xr5 */\
			__asm movaps xmm3,[ecx+__o4]	\
			__asm movaps xmm1,xmm3			\
			__asm subpd xmm3,xmm5			\
			__asm addpd xmm1,xmm5			\
			__asm movaps [ecx+__o4],xmm3	\
			__asm movaps [ecx+__o9],xmm1	\
			/* __B1r -= xr4;	xr4 *= xr1; xr4 += __B1r;   __Bcr = xr4 */\
			__asm movaps xmm3,[ecx+__o1]	\
			__asm movaps xmm1,xmm3			\
			__asm subpd xmm3,xmm4			\
			__asm addpd xmm1,xmm4			\
			__asm movaps [ecx+__o1],xmm3	\
			__asm movaps [ecx+__oC],xmm1	\
		/***************/\
		/* IMAG PARTS: Swap __**r <--> __**i, Replace __B[j] with __B[13-j] for j > 0: */\
		/***************/\
			__asm	add	ecx, 0x10	\
			__asm	movaps	xmm7,[eax+__i6]		/* xr7 = __A6i */\
			__asm	movaps	xmm5,[eax+__i5]		/* xr5 = __A5i */\
			__asm	movaps	xmm4,[eax+__i4]		/* xr4 = __A4i */\
			__asm	movaps	xmm6,[eax+__i3]		/* xr6 = __A3i */\
			__asm	movaps	xmm3,[eax+__i2]		/* xr3 = __A2i */\
			__asm	movaps	xmm1,[eax+__i1]		/* xr1 = __A1i */\
		/* xi-terms need 8 registers for each side: */\
			__asm	movaps	xmm0,[ebx-0x10]		/* 2.0 */\
			__asm	addpd	xmm7,[eax+__i7]		/* xr7 += __A7i */\
			__asm	addpd	xmm5,[eax+__i8]		/* xr5 += __A8i */\
			__asm	addpd	xmm4,[eax+__i9]		/* xr4 += __A9i */\
			__asm	addpd	xmm6,[eax+__iA]		/* xr6 += __Aai */\
			__asm	addpd	xmm3,[eax+__iB]		/* xr3 += __Abi */\
			__asm	addpd	xmm1,[eax+__iC]		/* xr1 += __Aci */\
			__asm	subpd	xmm1,xmm5			__asm	mulpd	xmm5,xmm0	/* xr1 -= xr5;	xr5 *= 2.0 */\
			__asm	subpd	xmm3,xmm6			__asm	mulpd	xmm6,xmm0	/* xr3 -= xr6;	xr6 *= 2.0 */\
			__asm	subpd	xmm4,xmm7			__asm	mulpd	xmm7,xmm0	/* xr4 -= xr7;	xr7 *= 2.0 */\
			__asm	addpd	xmm5,xmm1			/* xr5 += xr1 */\
			__asm	addpd	xmm6,xmm3			/* xr6 += xr3 */\
			__asm	addpd	xmm7,xmm4			/* xr7 += xr4 */\
			__asm	movaps	xmm2,xmm5			/* xr2  = xr5 */\
			__asm	addpd	xmm2,xmm6			/* xr2 += xr6 */\
			__asm	addpd	xmm2,xmm7			/* xr2 += xr7 */\
			__asm	movaps	xmm0,[eax]			/* xr0 = __A0i */\
			__asm	addpd	xmm0,xmm2			/* xr0 += xr2 */\
			__asm	mulpd	xmm2,[ebx]			/* xr2 *= DC1 */\
			__asm	movaps	[ecx],xmm0			/* __B0i = xr0 */\
			__asm	addpd	xmm0,xmm2			/* xr0 += xr2 */\
			__asm	movaps	xmm2,[ebx+0x010]	/* xr2 = DC3 */\
			__asm	mulpd	xmm6,xmm2			/* xr6 *= xr2 */\
			__asm	mulpd	xmm5,xmm2			/* xr5 *= xr2 */\
			__asm	mulpd	xmm7,xmm2			/* xr7 *= xr2 */\
			__asm	movaps  [ecx+__oC],xmm6		/* __Bci = xr6 */\
			__asm	movaps  [ecx+__oB],xmm5		/* __Bbi = xr5 */\
			__asm	movaps  [ecx+__oA],xmm7		/* __Bai = xr7 */\
			__asm	movaps	xmm2,[ebx+0x080]	/* xr2 = DC23 */\
			__asm	mulpd	xmm5,xmm2			/* xr5 *= xr2 */\
			__asm	mulpd	xmm7,xmm2			/* xr7 *= xr2 */\
			__asm	mulpd	xmm6,xmm2			/* xr6 *= xr2 */\
			__asm	addpd	xmm5,xmm0			/* xr5 += xr0 */\
			__asm	addpd	xmm7,xmm0			/* xr7 += xr0 */\
			__asm	addpd	xmm6,xmm0			/* xr6 += xr0 */\
			__asm	addpd	xmm5,[ecx+__oC]		/* xr5 += __Bci */\
			__asm	addpd	xmm7,[ecx+__oB]		/* xr7 += __Bbi */\
			__asm	addpd	xmm6,[ecx+__oA]		/* xr6 += __Bai */\
		/* Rearrange the 3x3 multiply block to ease in-place computation: */\
			__asm	movaps	xmm2,[ebx+0x0a0]	/* xr2 = DC65 */\
			__asm	movaps	xmm0,[ebx+0x090]	/* xr0 = DC54 */\
			__asm	movaps  [ecx+__oC],xmm4		/* __Bci = xr4 */\
			__asm	movaps  [ecx+__oB],xmm1		/* __Bbi = xr1 */\
			__asm	movaps  [ecx+__oA],xmm3		/* __Bai = xr3 */\
			__asm	mulpd	xmm4,xmm2			/* xr4 *= DC65 */\
			__asm	mulpd	xmm1,xmm2			/* xr1 *= DC65 */\
			__asm	mulpd	xmm3,xmm2			/* xr3 *= DC65 */\
			__asm	addpd	xmm4,[ecx+__oA]		/* xr4 += __Bai */\
			__asm	subpd	xmm1,[ecx+__oC]		/* xr1 -= __Bci */\
			__asm	addpd	xmm3,[ecx+__oB]		/* xr3 += __Bbi */\
			__asm	mulpd	xmm4,xmm0			/* xr4 *= DC54 */\
			__asm	mulpd	xmm1,xmm0			/* xr1 *= DC54 */\
			__asm	mulpd	xmm3,xmm0			/* xr3 *= DC54 */\
			__asm	movaps	xmm2,[ebx+0x020]	/* xr2 = DC4 */\
			__asm	addpd	xmm4,[ecx+__oB]		/* xr4 += __Bbi */\
			__asm	subpd	xmm1,[ecx+__oA]		/* xr1 -= __Bai */\
			__asm	subpd	xmm3,[ecx+__oC]		/* xr3 -= __Bci */\
			__asm	mulpd	xmm4,xmm2			/* xr4 *= DC4 */\
			__asm	mulpd	xmm1,xmm2			/* xr1 *= DC4 */\
			__asm	mulpd	xmm3,xmm2			/* xr3 *= DC4 */\
			__asm	movaps	xmm0,[ebx-0x010]	/* 2.0 */\
		/* Spill into destination outputs: */\
			__asm	subpd	xmm6,xmm4			__asm	mulpd	xmm4,xmm0	/* xr6 -= xr1;	xr5 *= 2.0 */\
			__asm	subpd	xmm7,xmm1			__asm	mulpd	xmm1,xmm0	/* xr7 -= xr3;	xr6 *= 2.0 */\
			__asm	subpd	xmm5,xmm3			__asm	mulpd	xmm3,xmm0	/* xr5 -= xr4;	xr7 *= 2.0 */\
			__asm	movaps  [ecx+__oA],xmm7		/* __Bai = xr7 */\
			__asm	movaps	[ecx+__o5],xmm5		/* __B5i = xr5 */\
			__asm	addpd	xmm4,xmm6			/* xr4 += xr6 */\
			__asm	addpd	xmm1,xmm7			/* xr7 += xr1 */\
			__asm	addpd	xmm3,xmm5			/* xr5 += xr3 */\
			__asm	movaps	[ecx+__o9],xmm6		/* __B9i = xr6 */\
			__asm	movaps	[ecx+__o7],xmm4		/* __B7i = xr4 */\
			__asm	movaps  [ecx+__oB],xmm1		/* __Bbi = xr7 */\
			__asm	movaps  [ecx+__oC],xmm3		/* __Bci = xr5 */\
		/* yr-terms: */\
			__asm	sub	eax, 0x10	\
			__asm	movaps	xmm1,[eax+__i1]		/* xr1 = __A1r */\
			__asm	movaps	xmm2,[eax+__i2]		/* xr2 = __A2r */\
			__asm	movaps	xmm5,[eax+__i3]		/* xr5 = __A3r */\
			__asm	movaps	xmm3,[eax+__i4]		/* xr3 = __A4r */\
			__asm	movaps	xmm4,[eax+__i5]		/* xr4 = __A5r */\
			__asm	movaps	xmm6,[eax+__i6]		/* xr6 = __A6r */\
			__asm	subpd	xmm1,[eax+__iC]		/* xr1 -= __Acr */\
			__asm	subpd	xmm2,[eax+__iB]		/* xr2 -= __Abr */\
			__asm	subpd	xmm5,[eax+__iA]		/* xr5 -= __Aar */\
			__asm	subpd	xmm3,[eax+__i9]		/* xr3 -= __A9r */\
			__asm	subpd	xmm4,[eax+__i8]		/* xr4 -= __A8r */\
			__asm	subpd	xmm6,[eax+__i7]		/* xr6 -= __A7r */\
			__asm	movaps	xmm7,xmm1			/* xr7 = xr1 */\
			__asm	movaps	xmm0,xmm2			/* xr0 = xr2 */\
			__asm	subpd	xmm7,xmm3			/* xr7 -= xr3 */\
			__asm	addpd	xmm0,xmm4			/* xr0 += xr4 */\
			__asm	addpd	xmm7,xmm5			/* xr7 += xr5 */\
			__asm	addpd	xmm0,xmm6			/* xr0 += xr6 */\
			__asm	addpd	xmm1,xmm3			/* xr1 += xr3 */\
			__asm	addpd	xmm5,xmm3			/* xr5 += xr3 */\
			__asm	subpd	xmm2,xmm4			/* xr2 -= xr4 */\
			__asm	subpd	xmm6,xmm4			/* xr6 -= xr4 */\
			__asm	movaps	xmm4,xmm0			/* xr4 = xr0 */\
			__asm	movaps	xmm3,xmm7			/* xr3	= xr7 */\
			__asm	mulpd	xmm0,[ebx+0x040]	/* xr0 *= DS2 */\
			__asm	mulpd	xmm7,[ebx+0x040]	/* xr7 *= DS2 */\
			__asm	mulpd	xmm4,[ebx+0x030]	/* xr4 *= DS1 */\
			__asm	mulpd	xmm3,[ebx+0x030]	/* xr3	*= DS1 */\
			__asm	subpd	xmm4,xmm7			/* xr4 -= xr7 */\
			__asm	addpd	xmm3,xmm0			/* xr3	+= xr0 */\
			__asm	movaps	[ecx+__o1],xmm4		/* __B1i = xr4; tmp-store in B1i */\
			__asm	movaps	xmm0,xmm1			/* xr0 = xr1 */\
			__asm	movaps	xmm4,xmm5			/* xr4 = xr5 */\
			__asm	addpd	xmm0,xmm2			/* xr0 += xr2 */\
			__asm	addpd	xmm4,xmm6			/* xr4 += xr6 */\
		/*\
		xr7 = DS3*xr0-DS6*xr4;\
		xr0 = DS3*xr4-DS9*xr0;\
		*/\
			__asm	movaps	[ecx+__o2],xmm0		/*	__B2i = xr0; tmp-store in B2i */\
			__asm	movaps	xmm7,xmm4			/*	  xr7 = xr4 */\
			__asm	mulpd	xmm4,[ebx+0x0b0]	/*	xr4 *= DS63 */\
			__asm	mulpd	xmm0,[ebx+0x0e0]	/*	xr0 *= DS93 */\
			__asm	addpd	xmm4,[ecx+__o2]		/*	xr4 += __B2i */\
			__asm	addpd	xmm0,xmm7			/*	xr0 +=	xr7 */\
			__asm	mulpd	xmm4,[ebx+0x050]	/*	xr4 *= DS3 */\
			__asm	mulpd	xmm0,[ebx+0x050]	/*	xr0 *= DS3 */\
		/*\
		xmm7,xmm4+DS4*xmm2-DS7*xr6;\
		xmm2,xmm0+DS4*xmm6-DSa*xr2;\
		*/\
			__asm	movaps	[ecx+__o2],xmm2		/*	__B2i = xr2 */\
			__asm	movaps	xmm7,xmm6			/*	  xr7 = xr6 */\
			__asm	mulpd	xmm6,[ebx+0x0c0]	/*	xr6 *= DS74 */\
			__asm	mulpd	xmm2,[ebx+0x0f0]	/*	xr2 *= DSa4 */\
			__asm	addpd	xmm6,[ecx+__o2]		/*	xr6 += __B2i */\
			__asm	addpd	xmm2,xmm7			/*	xr2 +=	xr7 */\
			__asm	mulpd	xmm6,[ebx+0x060]	/*	xr6 *= DS4 */\
			__asm	mulpd	xmm2,[ebx+0x060]	/*	xr2 *= DS4 */\
			__asm	addpd	xmm6,xmm4			/*	xr6 += xr4 */\
			__asm	addpd	xmm2,xmm0			/*	xr2 += xr0 */\
		/*\
		xmm7,xmm4+DS5*xmm1-DS8*xr5;\
		xmm0,xmm0+DS5*xmm5-DSb*xr1;\
		*/\
			__asm	movaps	[ecx+__o2],xmm1		/*	__B2i = xr1 */\
			__asm	movaps	xmm7,xmm5			/*	  xr7 = xr5 */\
			__asm	mulpd	xmm5,[ebx+0x0d0]	/*	xr5 *= DS85 */\
			__asm	mulpd	xmm1,[ebx+0x100]	/*	xr1 *= DSb5 */\
			__asm	addpd	xmm5,[ecx+__o2]		/*	xr5 += __B2i */\
			__asm	addpd	xmm1,xmm7			/*	xr1 +=	xr7 */\
			__asm	mulpd	xmm5,[ebx+0x070]	/*	xr5 *= DS5 */\
			__asm	mulpd	xmm1,[ebx+0x070]	/*	xr1 *= DS5 */\
			__asm	addpd	xmm5,xmm4			/*	xr5 += xr4 */\
			__asm	addpd	xmm1,xmm0			/*	xr1 += xr0 */\
			__asm	movaps	xmm4,[ecx+__o1]		/*	xr4 = __B1i */\
			__asm	movaps	xmm7,xmm3			/*	xr7  = xr3 */\
			__asm	movaps	xmm0,xmm4			/*	xr0  = xr4 */\
			__asm	addpd	xmm7,xmm6			/*	xr7 += xr6 */\
			__asm	addpd	xmm0,xmm5			/*	xr0 += xr5 */\
			__asm	subpd	xmm6,xmm3			/*	xr6 -= xr3 */\
			__asm	subpd	xmm5,xmm4			/*	xr5 -= xr4 */\
			__asm	addpd	xmm6,xmm2			/*	xr6 += xr2 */\
			__asm	addpd	xmm5,xmm1			/*	xr5 += xr1 */\
			__asm	addpd	xmm2,xmm3			/*	xr2 += xr3 */\
			__asm	addpd	xmm4,xmm1			/*	xr4 += xr1 */\
		/* Combine xmm and yi-terms to get real parts of outputs: */\
			/* __B7i -= xr7; xr7 *= xr1; xr7 += __B7i; __B6i = xr7 */\
			__asm movaps xmm3,[ecx+__o7]	\
			__asm movaps xmm1,xmm3			\
			__asm subpd xmm3,xmm7			\
			__asm addpd xmm1,xmm7			\
			__asm movaps [ecx+__o7],xmm3	\
			__asm movaps [ecx+__o6],xmm1	\
			/* __B5i -= xr6; xr6 *= xr1; xr6+= __B5i; __B8i = xr6 */\
			__asm movaps xmm3,[ecx+__o5]	\
			__asm movaps xmm1,xmm3			\
			__asm subpd xmm3,xmm6			\
			__asm addpd xmm1,xmm6			\
			__asm movaps [ecx+__o5],xmm3	\
			__asm movaps [ecx+__o8],xmm1	\
			/* __Bbi -= xr2; xr2 *= xr1; xr2+= __Bbi; __B2i = xr2 */\
			__asm movaps xmm3,[ecx+__oB]	\
			__asm movaps xmm1,xmm3			\
			__asm subpd xmm3,xmm2			\
			__asm addpd xmm1,xmm2			\
			__asm movaps [ecx+__oB],xmm3	\
			__asm movaps [ecx+__o2],xmm1	\
			/* __Bai -= xr0; xr0 *= xr1; xr0+= __Bai; __B3i = xr0 */\
			__asm movaps xmm3,[ecx+__oA]	\
			__asm movaps xmm1,xmm3			\
			__asm subpd xmm3,xmm0			\
			__asm addpd xmm1,xmm0			\
			__asm movaps [ecx+__oA],xmm3	\
			__asm movaps [ecx+__o3],xmm1	\
			/* __B9i -= xr5; xr5 *= xr1; xr5+= __B9i; __B4i = xr5 */\
			__asm movaps xmm3,[ecx+__o9]	\
			__asm movaps xmm1,xmm3			\
			__asm subpd xmm3,xmm5			\
			__asm addpd xmm1,xmm5			\
			__asm movaps [ecx+__o9],xmm3	\
			__asm movaps [ecx+__o4],xmm1	\
			/* __Bci -= xr4; xr4 *= xr1; xr4+= __Bci; __B1i = xr4 */\
			__asm movaps xmm3,[ecx+__oC]	\
			__asm movaps xmm1,xmm3			\
			__asm subpd xmm3,xmm4			\
			__asm addpd xmm1,xmm4			\
			__asm movaps [ecx+__oC],xmm3	\
			__asm movaps [ecx+__o1],xmm1	\
		}

	#else	/* GCC-style inline ASM: */

	  #if OS_BITS == 32

		#define SSE2_RADIX_13_DFT(Xcc, XI0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8,Xi9,XiA,XiB,XiC, XO0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8,Xo9,XoA,XoB,XoC)\
		{\
		__asm__ volatile (\
			"movl	%[__I0],%%eax			\n\t"\
			"movl	%[__cc],%%ebx			\n\t"\
			"movl	%[__O0],%%ecx			\n\t"\
		"/* xr-terms need 8 registers for each side: */\n\t"\
			"movaps	%c[__i6](%%eax),%%xmm7	\n\t"\
			"movaps	%c[__i5](%%eax),%%xmm5	\n\t"\
			"movaps	%c[__i4](%%eax),%%xmm4	\n\t"\
			"movaps	%c[__i3](%%eax),%%xmm6	\n\t"\
			"movaps	%c[__i2](%%eax),%%xmm3	\n\t"\
			"movaps	%c[__i1](%%eax),%%xmm1	\n\t"\
			"movaps	-0x010(%%ebx),%%xmm0	\n\t"\
			"addpd	%c[__i7](%%eax),%%xmm7	\n\t"\
			"addpd	%c[__i8](%%eax),%%xmm5	\n\t"\
			"addpd	%c[__i9](%%eax),%%xmm4	\n\t"\
			"addpd	%c[__iA](%%eax),%%xmm6	\n\t"\
			"addpd	%c[__iB](%%eax),%%xmm3	\n\t"\
			"addpd	%c[__iC](%%eax),%%xmm1	\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t"\
			"mulpd	%%xmm0,%%xmm5			\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t"\
			"mulpd	%%xmm0,%%xmm6			\n\t"\
			"subpd	%%xmm7,%%xmm4			\n\t"\
			"mulpd	%%xmm0,%%xmm7			\n\t"\
			"addpd	%%xmm1,%%xmm5			\n\t"\
			"addpd	%%xmm3,%%xmm6			\n\t"\
			"addpd	%%xmm4,%%xmm7			\n\t"\
			"movaps	%%xmm5,%%xmm2			\n\t"\
			"addpd	%%xmm6,%%xmm2			\n\t"\
			"addpd	%%xmm7,%%xmm2			\n\t"\
			"movaps	(%%eax),%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm0			\n\t"\
			"mulpd	(%%ebx),%%xmm2			\n\t"\
			"movaps	%%xmm0,(%%ecx)			\n\t"\
			"addpd	%%xmm2,%%xmm0			\n\t"\
			"movaps	0x010(%%ebx),%%xmm2		\n\t"\
			"mulpd	%%xmm2,%%xmm6			\n\t"\
			"mulpd	%%xmm2,%%xmm5			\n\t"\
			"mulpd	%%xmm2,%%xmm7			\n\t"\
			"movaps	%%xmm6,%c[__o1](%%ecx)	\n\t"\
			"movaps	%%xmm5,%c[__o2](%%ecx)	\n\t"\
			"movaps	%%xmm7,%c[__o3](%%ecx)	\n\t"\
			"movaps	0x080(%%ebx),%%xmm2		\n\t"\
			"mulpd	%%xmm2,%%xmm5			\n\t"\
			"mulpd	%%xmm2,%%xmm7			\n\t"\
			"mulpd	%%xmm2,%%xmm6			\n\t"\
			"addpd	%%xmm0,%%xmm5			\n\t"\
			"addpd	%%xmm0,%%xmm7			\n\t"\
			"addpd	%%xmm0,%%xmm6			\n\t"\
			"addpd	%c[__o1](%%ecx),%%xmm5	\n\t"\
			"addpd	%c[__o2](%%ecx),%%xmm7	\n\t"\
			"addpd	%c[__o3](%%ecx),%%xmm6	\n\t"\
			"movaps	0x0a0(%%ebx),%%xmm2		\n\t"\
			"movaps	0x090(%%ebx),%%xmm0		\n\t"\
			"movaps	%%xmm4,%c[__o1](%%ecx)	\n\t"\
			"movaps	%%xmm1,%c[__o2](%%ecx)	\n\t"\
			"movaps	%%xmm3,%c[__o3](%%ecx)	\n\t"\
			"mulpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	%%xmm2,%%xmm1			\n\t"\
			"mulpd	%%xmm2,%%xmm3			\n\t"\
			"addpd	%c[__o3](%%ecx),%%xmm4	\n\t"\
			"subpd	%c[__o1](%%ecx),%%xmm1	\n\t"\
			"addpd	%c[__o2](%%ecx),%%xmm3	\n\t"\
			"mulpd	%%xmm0,%%xmm4			\n\t"\
			"mulpd	%%xmm0,%%xmm1			\n\t"\
			"mulpd	%%xmm0,%%xmm3			\n\t"\
			"movaps	0x020(%%ebx),%%xmm2		\n\t"\
			"addpd	%c[__o2](%%ecx),%%xmm4	\n\t"\
			"subpd	%c[__o3](%%ecx),%%xmm1	\n\t"\
			"subpd	%c[__o1](%%ecx),%%xmm3	\n\t"\
			"mulpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	%%xmm2,%%xmm1			\n\t"\
			"mulpd	%%xmm2,%%xmm3			\n\t"\
			"movaps	-0x010(%%ebx),%%xmm0	\n\t"\
			"subpd	%%xmm4,%%xmm6			\n\t"\
			"mulpd	%%xmm0,%%xmm4			\n\t"\
			"subpd	%%xmm1,%%xmm7			\n\t"\
			"mulpd	%%xmm0,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm5			\n\t"\
			"mulpd	%%xmm0,%%xmm3			\n\t"\
			"movaps	%%xmm7,%c[__o3](%%ecx)	\n\t"\
			"movaps	%%xmm5,%c[__o8](%%ecx)	\n\t"\
			"addpd	%%xmm6,%%xmm4			\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"addpd	%%xmm3,%%xmm5			\n\t"\
			"movaps	%%xmm5,%c[__o1](%%ecx)	\n\t"\
			"movaps	%%xmm7,%c[__o2](%%ecx)	\n\t"\
			"movaps	%%xmm6,%c[__o4](%%ecx)	\n\t"\
			"movaps	%%xmm4,%c[__o6](%%ecx)	\n\t"\
		"/* yi-terms: */				\n\t"\
			"addl	$0x10,%%eax				\n\t"\
			"movaps	%c[__i1](%%eax),%%xmm1	\n\t"\
			"movaps	%c[__i2](%%eax),%%xmm2	\n\t"\
			"movaps	%c[__i3](%%eax),%%xmm5	\n\t"\
			"movaps	%c[__i4](%%eax),%%xmm3	\n\t"\
			"movaps	%c[__i5](%%eax),%%xmm4	\n\t"\
			"movaps	%c[__i6](%%eax),%%xmm6	\n\t"\
			"subpd	%c[__iC](%%eax),%%xmm1	\n\t"\
			"subpd	%c[__iB](%%eax),%%xmm2	\n\t"\
			"subpd	%c[__iA](%%eax),%%xmm5	\n\t"\
			"subpd	%c[__i9](%%eax),%%xmm3	\n\t"\
			"subpd	%c[__i8](%%eax),%%xmm4	\n\t"\
			"subpd	%c[__i7](%%eax),%%xmm6	\n\t"\
			"movaps	%%xmm1,%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm0			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"addpd	%%xmm4,%%xmm0			\n\t"\
			"addpd	%%xmm5,%%xmm7			\n\t"\
			"addpd	%%xmm6,%%xmm0			\n\t"\
			"addpd	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm3,%%xmm5			\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"subpd	%%xmm4,%%xmm6			\n\t"\
			"movaps	%%xmm0,%%xmm4			\n\t"\
			"movaps	%%xmm7,%%xmm3			\n\t"\
			"mulpd	0x040(%%ebx),%%xmm0		\n\t"\
			"mulpd	0x040(%%ebx),%%xmm7		\n\t"\
			"mulpd	0x030(%%ebx),%%xmm4		\n\t"\
			"mulpd	0x030(%%ebx),%%xmm3		\n\t"\
			"subpd	%%xmm7,%%xmm4			\n\t"\
			"addpd	%%xmm0,%%xmm3			\n\t"\
			"movaps	%%xmm4,%c[__oC](%%ecx)	\n\t"\
			"movaps	%%xmm1,%%xmm0			\n\t"\
			"movaps	%%xmm5,%%xmm4			\n\t"\
			"addpd	%%xmm2,%%xmm0			\n\t"\
			"addpd	%%xmm6,%%xmm4			\n\t"\
			"movaps	%%xmm0,%c[__oB](%%ecx)	\n\t"\
			"movaps	%%xmm4,%%xmm7			\n\t"\
			"mulpd	0x0b0(%%ebx),%%xmm4		\n\t"\
			"mulpd	0x0e0(%%ebx),%%xmm0		\n\t"\
			"addpd	%c[__oB](%%ecx),%%xmm4	\n\t"\
			"addpd	%%xmm7,%%xmm0			\n\t"\
			"mulpd	0x050(%%ebx),%%xmm4		\n\t"\
			"mulpd	0x050(%%ebx),%%xmm0		\n\t"\
			"movaps	%%xmm2,%c[__oB](%%ecx)	\n\t"\
			"movaps	%%xmm6,%%xmm7			\n\t"\
			"mulpd	0x0c0(%%ebx),%%xmm6		\n\t"\
			"mulpd	0x0f0(%%ebx),%%xmm2		\n\t"\
			"addpd	%c[__oB](%%ecx),%%xmm6	\n\t"\
			"addpd	%%xmm7,%%xmm2			\n\t"\
			"mulpd	0x060(%%ebx),%%xmm6		\n\t"\
			"mulpd	0x060(%%ebx),%%xmm2		\n\t"\
			"addpd	%%xmm4,%%xmm6			\n\t"\
			"addpd	%%xmm0,%%xmm2			\n\t"\
			"movaps	%%xmm1,%c[__oB](%%ecx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"mulpd	0x0d0(%%ebx),%%xmm5		\n\t"\
			"mulpd	0x100(%%ebx),%%xmm1		\n\t"\
			"addpd	%c[__oB](%%ecx),%%xmm5	\n\t"\
			"addpd	%%xmm7,%%xmm1			\n\t"\
			"mulpd	0x070(%%ebx),%%xmm5		\n\t"\
			"mulpd	0x070(%%ebx),%%xmm1		\n\t"\
			"addpd	%%xmm4,%%xmm5			\n\t"\
			"addpd	%%xmm0,%%xmm1			\n\t"\
			"movaps	%c[__oC](%%ecx),%%xmm4	\n\t"\
			"movaps	%%xmm3,%%xmm7			\n\t"\
			"movaps	%%xmm4,%%xmm0			\n\t"\
			"addpd	%%xmm6,%%xmm7			\n\t"\
			"addpd	%%xmm5,%%xmm0			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"addpd	%%xmm2,%%xmm6			\n\t"\
			"addpd	%%xmm1,%%xmm5			\n\t"\
			"addpd	%%xmm3,%%xmm2			\n\t"\
			"addpd	%%xmm1,%%xmm4			\n\t"\
			"movaps	%c[__o6](%%ecx),%%xmm3	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"subpd	%%xmm7,%%xmm3			\n\t"\
			"addpd	%%xmm7,%%xmm1			\n\t"\
			"movaps	%%xmm3,%c[__o6](%%ecx)	\n\t"\
			"movaps	%%xmm1,%c[__o7](%%ecx)	\n\t"\
			"movaps	%c[__o8](%%ecx),%%xmm3	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t"\
			"addpd	%%xmm6,%%xmm1			\n\t"\
			"movaps	%%xmm3,%c[__o8](%%ecx)	\n\t"\
			"movaps	%%xmm1,%c[__o5](%%ecx)	\n\t"\
			"movaps	%c[__o2](%%ecx),%%xmm3	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"subpd	%%xmm2,%%xmm3			\n\t"\
			"addpd	%%xmm2,%%xmm1			\n\t"\
			"movaps	%%xmm3,%c[__o2](%%ecx)	\n\t"\
			"movaps	%%xmm1,%c[__oB](%%ecx)	\n\t"\
			"movaps	%c[__o3](%%ecx),%%xmm3	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"subpd	%%xmm0,%%xmm3			\n\t"\
			"addpd	%%xmm0,%%xmm1			\n\t"\
			"movaps	%%xmm3,%c[__o3](%%ecx)	\n\t"\
			"movaps	%%xmm1,%c[__oA](%%ecx)	\n\t"\
			"movaps	%c[__o4](%%ecx),%%xmm3	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"subpd	%%xmm5,%%xmm3			\n\t"\
			"addpd	%%xmm5,%%xmm1			\n\t"\
			"movaps	%%xmm3,%c[__o4](%%ecx)	\n\t"\
			"movaps	%%xmm1,%c[__o9](%%ecx)	\n\t"\
			"movaps	%c[__o1](%%ecx),%%xmm3	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"subpd	%%xmm4,%%xmm3			\n\t"\
			"addpd	%%xmm4,%%xmm1			\n\t"\
			"movaps	%%xmm3,%c[__o1](%%ecx)	\n\t"\
			"movaps	%%xmm1,%c[__oC](%%ecx)	\n\t"\
			"/***************/			\n\t"\
			"/* IMAG PARTS: */			\n\t"\
			"/***************/			\n\t"\
		"/* xi-terms need 8 registers for each	side: */\n\t"\
			"addl	$0x10,%%ecx				\n\t"\
			"movaps	%c[__i6](%%eax),%%xmm7	\n\t"\
			"movaps	%c[__i5](%%eax),%%xmm5	\n\t"\
			"movaps	%c[__i4](%%eax),%%xmm4	\n\t"\
			"movaps	%c[__i3](%%eax),%%xmm6	\n\t"\
			"movaps	%c[__i2](%%eax),%%xmm3	\n\t"\
			"movaps	%c[__i1](%%eax),%%xmm1	\n\t"\
			"movaps	-0x10(%%ebx),%%xmm0		\n\t"\
			"addpd	%c[__i7](%%eax),%%xmm7	\n\t"\
			"addpd	%c[__i8](%%eax),%%xmm5	\n\t"\
			"addpd	%c[__i9](%%eax),%%xmm4	\n\t"\
			"addpd	%c[__iA](%%eax),%%xmm6	\n\t"\
			"addpd	%c[__iB](%%eax),%%xmm3	\n\t"\
			"addpd	%c[__iC](%%eax),%%xmm1	\n\t"\
			"subpd	%%xmm5,%%xmm1			\n\t"\
			"mulpd	%%xmm0,%%xmm5			\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t"\
			"mulpd	%%xmm0,%%xmm6			\n\t"\
			"subpd	%%xmm7,%%xmm4			\n\t"\
			"mulpd	%%xmm0,%%xmm7			\n\t"\
			"addpd	%%xmm1,%%xmm5			\n\t"\
			"addpd	%%xmm3,%%xmm6			\n\t"\
			"addpd	%%xmm4,%%xmm7			\n\t"\
			"movaps	%%xmm5,%%xmm2			\n\t"\
			"addpd	%%xmm6,%%xmm2			\n\t"\
			"addpd	%%xmm7,%%xmm2			\n\t"\
			"movaps	(%%eax),%%xmm0			\n\t"\
			"addpd	%%xmm2,%%xmm0			\n\t"\
			"mulpd	(%%ebx),%%xmm2			\n\t"\
			"movaps	%%xmm0,(%%ecx)			\n\t"\
			"addpd	%%xmm2,%%xmm0			\n\t"\
			"movaps	0x010(%%ebx),%%xmm2		\n\t"\
			"mulpd	%%xmm2,%%xmm6			\n\t"\
			"mulpd	%%xmm2,%%xmm5			\n\t"\
			"mulpd	%%xmm2,%%xmm7			\n\t"\
			"movaps	%%xmm6,%c[__oC](%%ecx)	\n\t"\
			"movaps	%%xmm5,%c[__oB](%%ecx)	\n\t"\
			"movaps	%%xmm7,%c[__oA](%%ecx)	\n\t"\
			"movaps	0x080(%%ebx),%%xmm2		\n\t"\
			"mulpd	%%xmm2,%%xmm5			\n\t"\
			"mulpd	%%xmm2,%%xmm7			\n\t"\
			"mulpd	%%xmm2,%%xmm6			\n\t"\
			"addpd	%%xmm0,%%xmm5			\n\t"\
			"addpd	%%xmm0,%%xmm7			\n\t"\
			"addpd	%%xmm0,%%xmm6			\n\t"\
			"addpd	%c[__oC](%%ecx),%%xmm5	\n\t"\
			"addpd	%c[__oB](%%ecx),%%xmm7	\n\t"\
			"addpd	%c[__oA](%%ecx),%%xmm6	\n\t"\
			"movaps	0x0a0(%%ebx),%%xmm2		\n\t"\
			"movaps	0x090(%%ebx),%%xmm0		\n\t"\
			"movaps	%%xmm4,%c[__oC](%%ecx)	\n\t"\
			"movaps	%%xmm1,%c[__oB](%%ecx)	\n\t"\
			"movaps	%%xmm3,%c[__oA](%%ecx)	\n\t"\
			"mulpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	%%xmm2,%%xmm1			\n\t"\
			"mulpd	%%xmm2,%%xmm3			\n\t"\
			"addpd	%c[__oA](%%ecx),%%xmm4	\n\t"\
			"subpd	%c[__oC](%%ecx),%%xmm1	\n\t"\
			"addpd	%c[__oB](%%ecx),%%xmm3	\n\t"\
			"mulpd	%%xmm0,%%xmm4			\n\t"\
			"mulpd	%%xmm0,%%xmm1			\n\t"\
			"mulpd	%%xmm0,%%xmm3			\n\t"\
			"movaps	0x020(%%ebx),%%xmm2		\n\t"\
			"addpd	%c[__oB](%%ecx),%%xmm4	\n\t"\
			"subpd	%c[__oA](%%ecx),%%xmm1	\n\t"\
			"subpd	%c[__oC](%%ecx),%%xmm3	\n\t"\
			"mulpd	%%xmm2,%%xmm4			\n\t"\
			"mulpd	%%xmm2,%%xmm1			\n\t"\
			"mulpd	%%xmm2,%%xmm3			\n\t"\
			"movaps	-0x010(%%ebx),%%xmm0	\n\t"\
			"subpd	%%xmm4,%%xmm6			\n\t"\
			"mulpd	%%xmm0,%%xmm4			\n\t"\
			"subpd	%%xmm1,%%xmm7			\n\t"\
			"mulpd	%%xmm0,%%xmm1			\n\t"\
			"subpd	%%xmm3,%%xmm5			\n\t"\
			"mulpd	%%xmm0,%%xmm3			\n\t"\
			"movaps	%%xmm7,%c[__oA](%%ecx)	\n\t"\
			"movaps	%%xmm5,%c[__o5](%%ecx)	\n\t"\
			"addpd	%%xmm6,%%xmm4			\n\t"\
			"addpd	%%xmm1,%%xmm7			\n\t"\
			"addpd	%%xmm3,%%xmm5			\n\t"\
			"movaps	%%xmm5,%c[__oC](%%ecx)	\n\t"\
			"movaps	%%xmm7,%c[__oB](%%ecx)	\n\t"\
			"movaps	%%xmm6,%c[__o9](%%ecx)	\n\t"\
			"movaps	%%xmm4,%c[__o7](%%ecx)	\n\t"\
		"/* yr-terms: */				\n\t"\
			"subl	$0x10,%%eax				\n\t"\
			"movaps	%c[__i1](%%eax),%%xmm1	\n\t"\
			"movaps	%c[__i2](%%eax),%%xmm2	\n\t"\
			"movaps	%c[__i3](%%eax),%%xmm5	\n\t"\
			"movaps	%c[__i4](%%eax),%%xmm3	\n\t"\
			"movaps	%c[__i5](%%eax),%%xmm4	\n\t"\
			"movaps	%c[__i6](%%eax),%%xmm6	\n\t"\
			"subpd	%c[__iC](%%eax),%%xmm1	\n\t"\
			"subpd	%c[__iB](%%eax),%%xmm2	\n\t"\
			"subpd	%c[__iA](%%eax),%%xmm5	\n\t"\
			"subpd	%c[__i9](%%eax),%%xmm3	\n\t"\
			"subpd	%c[__i8](%%eax),%%xmm4	\n\t"\
			"subpd	%c[__i7](%%eax),%%xmm6	\n\t"\
			"movaps	%%xmm1,%%xmm7			\n\t"\
			"movaps	%%xmm2,%%xmm0			\n\t"\
			"subpd	%%xmm3,%%xmm7			\n\t"\
			"addpd	%%xmm4,%%xmm0			\n\t"\
			"addpd	%%xmm5,%%xmm7			\n\t"\
			"addpd	%%xmm6,%%xmm0			\n\t"\
			"addpd	%%xmm3,%%xmm1			\n\t"\
			"addpd	%%xmm3,%%xmm5			\n\t"\
			"subpd	%%xmm4,%%xmm2			\n\t"\
			"subpd	%%xmm4,%%xmm6			\n\t"\
			"movaps	%%xmm0,%%xmm4			\n\t"\
			"movaps	%%xmm7,%%xmm3			\n\t"\
			"mulpd	0x040(%%ebx),%%xmm0		\n\t"\
			"mulpd	0x040(%%ebx),%%xmm7		\n\t"\
			"mulpd	0x030(%%ebx),%%xmm4		\n\t"\
			"mulpd	0x030(%%ebx),%%xmm3		\n\t"\
			"subpd	%%xmm7,%%xmm4			\n\t"\
			"addpd	%%xmm0,%%xmm3			\n\t"\
			"movaps	%%xmm4,%c[__o1](%%ecx)	\n\t"\
			"movaps	%%xmm1,%%xmm0			\n\t"\
			"movaps	%%xmm5,%%xmm4			\n\t"\
			"addpd	%%xmm2,%%xmm0			\n\t"\
			"addpd	%%xmm6,%%xmm4			\n\t"\
			"movaps	%%xmm0,%c[__o2](%%ecx)	\n\t"\
			"movaps	%%xmm4,%%xmm7			\n\t"\
			"mulpd	0x0b0(%%ebx),%%xmm4		\n\t"\
			"mulpd	0x0e0(%%ebx),%%xmm0		\n\t"\
			"addpd	%c[__o2](%%ecx),%%xmm4	\n\t"\
			"addpd	%%xmm7,%%xmm0			\n\t"\
			"mulpd	0x050(%%ebx),%%xmm4		\n\t"\
			"mulpd	0x050(%%ebx),%%xmm0		\n\t"\
			"movaps	%%xmm2,%c[__o2](%%ecx)	\n\t"\
			"movaps	%%xmm6,%%xmm7			\n\t"\
			"mulpd	0x0c0(%%ebx),%%xmm6		\n\t"\
			"mulpd	0x0f0(%%ebx),%%xmm2		\n\t"\
			"addpd	%c[__o2](%%ecx),%%xmm6	\n\t"\
			"addpd	%%xmm7,%%xmm2			\n\t"\
			"mulpd	0x060(%%ebx),%%xmm6		\n\t"\
			"mulpd	0x060(%%ebx),%%xmm2		\n\t"\
			"addpd	%%xmm4,%%xmm6			\n\t"\
			"addpd	%%xmm0,%%xmm2			\n\t"\
			"movaps	%%xmm1,%c[__o2](%%ecx)	\n\t"\
			"movaps	%%xmm5,%%xmm7			\n\t"\
			"mulpd	0x0d0(%%ebx),%%xmm5		\n\t"\
			"mulpd	0x100(%%ebx),%%xmm1		\n\t"\
			"addpd	%c[__o2](%%ecx),%%xmm5	\n\t"\
			"addpd	%%xmm7,%%xmm1			\n\t"\
			"mulpd	0x070(%%ebx),%%xmm5		\n\t"\
			"mulpd	0x070(%%ebx),%%xmm1		\n\t"\
			"addpd	%%xmm4,%%xmm5			\n\t"\
			"addpd	%%xmm0,%%xmm1			\n\t"\
			"movaps	%c[__o1](%%ecx),%%xmm4	\n\t"\
			"movaps	%%xmm3,%%xmm7			\n\t"\
			"movaps	%%xmm4,%%xmm0			\n\t"\
			"addpd	%%xmm6,%%xmm7			\n\t"\
			"addpd	%%xmm5,%%xmm0			\n\t"\
			"subpd	%%xmm3,%%xmm6			\n\t"\
			"subpd	%%xmm4,%%xmm5			\n\t"\
			"addpd	%%xmm2,%%xmm6			\n\t"\
			"addpd	%%xmm1,%%xmm5			\n\t"\
			"addpd	%%xmm3,%%xmm2			\n\t"\
			"addpd	%%xmm1,%%xmm4			\n\t"\
			"movaps	%c[__o7](%%ecx),%%xmm3	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"subpd	%%xmm7,%%xmm3			\n\t"\
			"addpd	%%xmm7,%%xmm1			\n\t"\
			"movaps	%%xmm3,%c[__o7](%%ecx)	\n\t"\
			"movaps	%%xmm1,%c[__o6](%%ecx)	\n\t"\
			"movaps	%c[__o5](%%ecx),%%xmm3	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"subpd	%%xmm6,%%xmm3			\n\t"\
			"addpd	%%xmm6,%%xmm1			\n\t"\
			"movaps	%%xmm3,%c[__o5](%%ecx)	\n\t"\
			"movaps	%%xmm1,%c[__o8](%%ecx)	\n\t"\
			"movaps	%c[__oB](%%ecx),%%xmm3	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"subpd	%%xmm2,%%xmm3			\n\t"\
			"addpd	%%xmm2,%%xmm1			\n\t"\
			"movaps	%%xmm3,%c[__oB](%%ecx)	\n\t"\
			"movaps	%%xmm1,%c[__o2](%%ecx)	\n\t"\
			"movaps	%c[__oA](%%ecx),%%xmm3	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"subpd	%%xmm0,%%xmm3			\n\t"\
			"addpd	%%xmm0,%%xmm1			\n\t"\
			"movaps	%%xmm3,%c[__oA](%%ecx)	\n\t"\
			"movaps	%%xmm1,%c[__o3](%%ecx)	\n\t"\
			"movaps	%c[__o9](%%ecx),%%xmm3	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"subpd	%%xmm5,%%xmm3			\n\t"\
			"addpd	%%xmm5,%%xmm1			\n\t"\
			"movaps	%%xmm3,%c[__o9](%%ecx)	\n\t"\
			"movaps	%%xmm1,%c[__o4](%%ecx)	\n\t"\
			"movaps	%c[__oC](%%ecx),%%xmm3	\n\t"\
			"movaps	%%xmm3,%%xmm1			\n\t"\
			"subpd	%%xmm4,%%xmm3			\n\t"\
			"addpd	%%xmm4,%%xmm1			\n\t"\
			"movaps	%%xmm3,%c[__oC](%%ecx)	\n\t"\
			"movaps	%%xmm1,%c[__o1](%%ecx)	\n\t"\
			:					/* outputs: none */\
			: [__cc] "m" (Xcc)	/* All inputs from memory addresses here */\
			 ,[__I0] "m" (XI0)\
			 ,[__i1] "e" (Xi1)\
			 ,[__i2] "e" (Xi2)\
			 ,[__i3] "e" (Xi3)\
			 ,[__i4] "e" (Xi4)\
			 ,[__i5] "e" (Xi5)\
			 ,[__i6] "e" (Xi6)\
			 ,[__i7] "e" (Xi7)\
			 ,[__i8] "e" (Xi8)\
			 ,[__i9] "e" (Xi9)\
			 ,[__iA] "e" (XiA)\
			 ,[__iB] "e" (XiB)\
			 ,[__iC] "e" (XiC)\
			 ,[__O0] "m" (XO0)\
			 ,[__o1] "e" (Xo1)\
			 ,[__o2] "e" (Xo2)\
			 ,[__o3] "e" (Xo3)\
			 ,[__o4] "e" (Xo4)\
			 ,[__o5] "e" (Xo5)\
			 ,[__o6] "e" (Xo6)\
			 ,[__o7] "e" (Xo7)\
			 ,[__o8] "e" (Xo8)\
			 ,[__o9] "e" (Xo9)\
			 ,[__oA] "e" (XoA)\
			 ,[__oB] "e" (XoB)\
			 ,[__oC] "e" (XoC)\
			: "eax","ebx","ecx"		 /* Clobbered registers */\
			);\
		}

	  #else

		#define USE_64BIT_ASM_STYLE	1	// My x86 timings indicate xmm0-15-using version of radix-13 DFT is about the same speed as 8-register, but keep both along with simple toggle

		#if !USE_64BIT_ASM_STYLE

		// Simple 64-bit-ified version of the above 32-bit ASM macro, using just xmm0-7
		#define SSE2_RADIX_13_DFT(Xcc, XI0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8,Xi9,XiA,XiB,XiC, XO0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8,Xo9,XoA,XoB,XoC)\
		{\
		__asm__ volatile (\
			"movq	%[__I0],%%rax		\n\t"\
			"movq	%[__cc],%%rbx		\n\t"\
			"movq	%[__O0],%%rcx		\n\t"\
		"/* xr-terms need 8 registers for each side: */\n\t"\
			"movaps %c[__i6](%%rax),%%xmm7	\n\t"\
			"movaps %c[__i5](%%rax),%%xmm5	\n\t"\
			"movaps %c[__i4](%%rax),%%xmm4	\n\t"\
			"movaps %c[__i3](%%rax),%%xmm6	\n\t"\
			"movaps %c[__i2](%%rax),%%xmm3	\n\t"\
			"movaps %c[__i1](%%rax),%%xmm1	\n\t"\
			"movaps -0x010(%%rbx),%%xmm0	\n\t"\
			"addpd	%c[__i7](%%rax),%%xmm7	\n\t"\
			"addpd	%c[__i8](%%rax),%%xmm5	\n\t"\
			"addpd	%c[__i9](%%rax),%%xmm4	\n\t"\
			"addpd	%c[__iA](%%rax),%%xmm6	\n\t"\
			"addpd	%c[__iB](%%rax),%%xmm3	\n\t"\
			"addpd	%c[__iC](%%rax),%%xmm1	\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"mulpd	%%xmm0,%%xmm5		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"mulpd	%%xmm0,%%xmm6		\n\t"\
			"subpd	%%xmm7,%%xmm4		\n\t"\
			"mulpd	%%xmm0,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"addpd	%%xmm4,%%xmm7		\n\t"\
			"movaps %%xmm5,%%xmm2		\n\t"\
			"addpd	%%xmm6,%%xmm2		\n\t"\
			"addpd	%%xmm7,%%xmm2		\n\t"\
			"movaps (%%rax),%%xmm0		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"mulpd	(%%rbx),%%xmm2		\n\t"\
			"movaps %%xmm0,(%%rcx)	\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"movaps 0x010(%%rbx),%%xmm2	 \n\t"\
			"mulpd	%%xmm2,%%xmm6		\n\t"\
			"mulpd	%%xmm2,%%xmm5		\n\t"\
			"mulpd	%%xmm2,%%xmm7		\n\t"\
			"movaps %%xmm6,%c[__o1](%%rcx)	\n\t"\
			"movaps %%xmm5,%c[__o2](%%rcx)	\n\t"\
			"movaps %%xmm7,%c[__o3](%%rcx)	\n\t"\
			"movaps 0x080(%%rbx),%%xmm2	 \n\t"\
			"mulpd	%%xmm2,%%xmm5		\n\t"\
			"mulpd	%%xmm2,%%xmm7		\n\t"\
			"mulpd	%%xmm2,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm5		\n\t"\
			"addpd	%%xmm0,%%xmm7		\n\t"\
			"addpd	%%xmm0,%%xmm6		\n\t"\
			"addpd	%c[__o1](%%rcx),%%xmm5	\n\t"\
			"addpd	%c[__o2](%%rcx),%%xmm7	\n\t"\
			"addpd	%c[__o3](%%rcx),%%xmm6	\n\t"\
			"movaps 0x0a0(%%rbx),%%xmm2	 \n\t"\
			"movaps 0x090(%%rbx),%%xmm0	 \n\t"\
			"movaps %%xmm4,%c[__o1](%%rcx)	\n\t"\
			"movaps %%xmm1,%c[__o2](%%rcx)	\n\t"\
			"movaps %%xmm3,%c[__o3](%%rcx)	\n\t"\
			"mulpd	%%xmm2,%%xmm4		\n\t"\
			"mulpd	%%xmm2,%%xmm1		\n\t"\
			"mulpd	%%xmm2,%%xmm3		\n\t"\
			"addpd	%c[__o3](%%rcx),%%xmm4	\n\t"\
			"subpd	%c[__o1](%%rcx),%%xmm1	\n\t"\
			"addpd	%c[__o2](%%rcx),%%xmm3	\n\t"\
			"mulpd	%%xmm0,%%xmm4		\n\t"\
			"mulpd	%%xmm0,%%xmm1		\n\t"\
			"mulpd	%%xmm0,%%xmm3		\n\t"\
			"movaps 0x020(%%rbx),%%xmm2	 \n\t"\
			"addpd	%c[__o2](%%rcx),%%xmm4	\n\t"\
			"subpd	%c[__o3](%%rcx),%%xmm1	\n\t"\
			"subpd	%c[__o1](%%rcx),%%xmm3	\n\t"\
			"mulpd	%%xmm2,%%xmm4		\n\t"\
			"mulpd	%%xmm2,%%xmm1		\n\t"\
			"mulpd	%%xmm2,%%xmm3		\n\t"\
			"movaps -0x010(%%rbx),%%xmm0	\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t"\
			"mulpd	%%xmm0,%%xmm4		\n\t"\
			"subpd	%%xmm1,%%xmm7		\n\t"\
			"mulpd	%%xmm0,%%xmm1		\n\t"\
			"subpd	%%xmm3,%%xmm5		\n\t"\
			"mulpd	%%xmm0,%%xmm3		\n\t"\
			"movaps %%xmm7,%c[__o3](%%rcx)	\n\t"\
			"movaps %%xmm5,%c[__o8](%%rcx)	\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t"\
			"addpd	%%xmm1,%%xmm7		\n\t"\
			"addpd	%%xmm3,%%xmm5		\n\t"\
			"movaps %%xmm5,%c[__o1](%%rcx)	\n\t"\
			"movaps %%xmm7,%c[__o2](%%rcx)	\n\t"\
			"movaps %%xmm6,%c[__o4](%%rcx)	\n\t"\
			"movaps %%xmm4,%c[__o6](%%rcx)	\n\t"\
		"/* yi-terms: */			\n\t"\
			"addq	$0x10,%%rax		 \n\t"\
			"movaps %c[__i1](%%rax),%%xmm1	\n\t"\
			"movaps %c[__i2](%%rax),%%xmm2	\n\t"\
			"movaps %c[__i3](%%rax),%%xmm5	\n\t"\
			"movaps %c[__i4](%%rax),%%xmm3	\n\t"\
			"movaps %c[__i5](%%rax),%%xmm4	\n\t"\
			"movaps %c[__i6](%%rax),%%xmm6	\n\t"\
			"subpd	%c[__iC](%%rax),%%xmm1	\n\t"\
			"subpd	%c[__iB](%%rax),%%xmm2	\n\t"\
			"subpd	%c[__iA](%%rax),%%xmm5	\n\t"\
			"subpd	%c[__i9](%%rax),%%xmm3	\n\t"\
			"subpd	%c[__i8](%%rax),%%xmm4	\n\t"\
			"subpd	%c[__i7](%%rax),%%xmm6	\n\t"\
			"movaps %%xmm1,%%xmm7		\n\t"\
			"movaps %%xmm2,%%xmm0		\n\t"\
			"subpd	%%xmm3,%%xmm7		\n\t"\
			"addpd	%%xmm4,%%xmm0		\n\t"\
			"addpd	%%xmm5,%%xmm7		\n\t"\
			"addpd	%%xmm6,%%xmm0		\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t"\
			"addpd	%%xmm3,%%xmm5		\n\t"\
			"subpd	%%xmm4,%%xmm2		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t"\
			"movaps %%xmm0,%%xmm4		\n\t"\
			"movaps %%xmm7,%%xmm3		\n\t"\
			"mulpd	0x040(%%rbx),%%xmm0	 \n\t"\
			"mulpd	0x040(%%rbx),%%xmm7	 \n\t"\
			"mulpd	0x030(%%rbx),%%xmm4	 \n\t"\
			"mulpd	0x030(%%rbx),%%xmm3	 \n\t"\
			"subpd	%%xmm7,%%xmm4		\n\t"\
			"addpd	%%xmm0,%%xmm3		\n\t"\
			"movaps %%xmm4,%c[__oC](%%rcx)	\n\t"\
			"movaps %%xmm1,%%xmm0		\n\t"\
			"movaps %%xmm5,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t"\
			"movaps %%xmm0,%c[__oB](%%rcx)	\n\t"\
			"movaps %%xmm4,%%xmm7		\n\t"\
			"mulpd	0x0b0(%%rbx),%%xmm4	 \n\t"\
			"mulpd	0x0e0(%%rbx),%%xmm0	 \n\t"\
			"addpd	%c[__oB](%%rcx),%%xmm4	\n\t"\
			"addpd	%%xmm7,%%xmm0		\n\t"\
			"mulpd	0x050(%%rbx),%%xmm4	 \n\t"\
			"mulpd	0x050(%%rbx),%%xmm0	 \n\t"\
			"movaps %%xmm2,%c[__oB](%%rcx)	\n\t"\
			"movaps %%xmm6,%%xmm7		\n\t"\
			"mulpd	0x0c0(%%rbx),%%xmm6	 \n\t"\
			"mulpd	0x0f0(%%rbx),%%xmm2	 \n\t"\
			"addpd	%c[__oB](%%rcx),%%xmm6	\n\t"\
			"addpd	%%xmm7,%%xmm2		\n\t"\
			"mulpd	0x060(%%rbx),%%xmm6	 \n\t"\
			"mulpd	0x060(%%rbx),%%xmm2	 \n\t"\
			"addpd	%%xmm4,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm2		\n\t"\
			"movaps %%xmm1,%c[__oB](%%rcx)	\n\t"\
			"movaps %%xmm5,%%xmm7		\n\t"\
			"mulpd	0x0d0(%%rbx),%%xmm5	 \n\t"\
			"mulpd	0x100(%%rbx),%%xmm1	 \n\t"\
			"addpd	%c[__oB](%%rcx),%%xmm5	\n\t"\
			"addpd	%%xmm7,%%xmm1		\n\t"\
			"mulpd	0x070(%%rbx),%%xmm5	 \n\t"\
			"mulpd	0x070(%%rbx),%%xmm1	 \n\t"\
			"addpd	%%xmm4,%%xmm5		\n\t"\
			"addpd	%%xmm0,%%xmm1		\n\t"\
			"movaps %c[__oC](%%rcx),%%xmm4	\n\t"\
			"movaps %%xmm3,%%xmm7		\n\t"\
			"movaps %%xmm4,%%xmm0		\n\t"\
			"addpd	%%xmm6,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm0		\n\t"\
			"subpd	%%xmm3,%%xmm6		\n\t"\
			"subpd	%%xmm4,%%xmm5		\n\t"\
			"addpd	%%xmm2,%%xmm6		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm2		\n\t"\
			"addpd	%%xmm1,%%xmm4		\n\t"\
			"movaps %c[__o6](%%rcx),%%xmm3	\n\t"\
			"movaps %%xmm3,%%xmm1		\n\t"\
			"subpd	%%xmm7,%%xmm3		\n\t"\
			"addpd	%%xmm7,%%xmm1		\n\t"\
			"movaps %%xmm3,%c[__o6](%%rcx)	\n\t"\
			"movaps %%xmm1,%c[__o7](%%rcx)	\n\t"\
			"movaps %c[__o8](%%rcx),%%xmm3	\n\t"\
			"movaps %%xmm3,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"addpd	%%xmm6,%%xmm1		\n\t"\
			"movaps %%xmm3,%c[__o8](%%rcx)	\n\t"\
			"movaps %%xmm1,%c[__o5](%%rcx)	\n\t"\
			"movaps %c[__o2](%%rcx),%%xmm3	\n\t"\
			"movaps %%xmm3,%%xmm1		\n\t"\
			"subpd	%%xmm2,%%xmm3		\n\t"\
			"addpd	%%xmm2,%%xmm1		\n\t"\
			"movaps %%xmm3,%c[__o2](%%rcx)	\n\t"\
			"movaps %%xmm1,%c[__oB](%%rcx)	\n\t"\
			"movaps %c[__o3](%%rcx),%%xmm3	\n\t"\
			"movaps %%xmm3,%%xmm1		\n\t"\
			"subpd	%%xmm0,%%xmm3		\n\t"\
			"addpd	%%xmm0,%%xmm1		\n\t"\
			"movaps %%xmm3,%c[__o3](%%rcx)	\n\t"\
			"movaps %%xmm1,%c[__oA](%%rcx)	\n\t"\
			"movaps %c[__o4](%%rcx),%%xmm3	\n\t"\
			"movaps %%xmm3,%%xmm1		\n\t"\
			"subpd	%%xmm5,%%xmm3		\n\t"\
			"addpd	%%xmm5,%%xmm1		\n\t"\
			"movaps %%xmm3,%c[__o4](%%rcx)	\n\t"\
			"movaps %%xmm1,%c[__o9](%%rcx)	\n\t"\
			"movaps %c[__o1](%%rcx),%%xmm3	\n\t"\
			"movaps %%xmm3,%%xmm1		\n\t"\
			"subpd	%%xmm4,%%xmm3		\n\t"\
			"addpd	%%xmm4,%%xmm1		\n\t"\
			"movaps %%xmm3,%c[__o1](%%rcx)	\n\t"\
			"movaps %%xmm1,%c[__oC](%%rcx)	\n\t"\
			"/***************/			\n\t"\
			"/* IMAG PARTS: */			\n\t"\
			"/***************/			\n\t"\
		"/* xi-terms need 8 registers for each	side: */\n\t"\
			"addq	$0x10,%%rcx		 \n\t"\
			"movaps %c[__i6](%%rax),%%xmm7	\n\t"\
			"movaps %c[__i5](%%rax),%%xmm5	\n\t"\
			"movaps %c[__i4](%%rax),%%xmm4	\n\t"\
			"movaps %c[__i3](%%rax),%%xmm6	\n\t"\
			"movaps %c[__i2](%%rax),%%xmm3	\n\t"\
			"movaps %c[__i1](%%rax),%%xmm1	\n\t"\
			"movaps -0x10(%%rbx),%%xmm0	 \n\t"\
			"addpd	%c[__i7](%%rax),%%xmm7	\n\t"\
			"addpd	%c[__i8](%%rax),%%xmm5	\n\t"\
			"addpd	%c[__i9](%%rax),%%xmm4	\n\t"\
			"addpd	%c[__iA](%%rax),%%xmm6	\n\t"\
			"addpd	%c[__iB](%%rax),%%xmm3	\n\t"\
			"addpd	%c[__iC](%%rax),%%xmm1	\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"mulpd	%%xmm0,%%xmm5		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"mulpd	%%xmm0,%%xmm6		\n\t"\
			"subpd	%%xmm7,%%xmm4		\n\t"\
			"mulpd	%%xmm0,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"addpd	%%xmm4,%%xmm7		\n\t"\
			"movaps %%xmm5,%%xmm2		\n\t"\
			"addpd	%%xmm6,%%xmm2		\n\t"\
			"addpd	%%xmm7,%%xmm2		\n\t"\
			"movaps (%%rax),%%xmm0		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"mulpd	(%%rbx),%%xmm2		\n\t"\
			"movaps %%xmm0,(%%rcx)	\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"movaps 0x010(%%rbx),%%xmm2	 \n\t"\
			"mulpd	%%xmm2,%%xmm6		\n\t"\
			"mulpd	%%xmm2,%%xmm5		\n\t"\
			"mulpd	%%xmm2,%%xmm7		\n\t"\
			"movaps %%xmm6,%c[__oC](%%rcx)	\n\t"\
			"movaps %%xmm5,%c[__oB](%%rcx)	\n\t"\
			"movaps %%xmm7,%c[__oA](%%rcx)	\n\t"\
			"movaps 0x080(%%rbx),%%xmm2	 \n\t"\
			"mulpd	%%xmm2,%%xmm5		\n\t"\
			"mulpd	%%xmm2,%%xmm7		\n\t"\
			"mulpd	%%xmm2,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm5		\n\t"\
			"addpd	%%xmm0,%%xmm7		\n\t"\
			"addpd	%%xmm0,%%xmm6		\n\t"\
			"addpd	%c[__oC](%%rcx),%%xmm5	\n\t"\
			"addpd	%c[__oB](%%rcx),%%xmm7	\n\t"\
			"addpd	%c[__oA](%%rcx),%%xmm6	\n\t"\
			"movaps 0x0a0(%%rbx),%%xmm2	 \n\t"\
			"movaps 0x090(%%rbx),%%xmm0	 \n\t"\
			"movaps %%xmm4,%c[__oC](%%rcx)	\n\t"\
			"movaps %%xmm1,%c[__oB](%%rcx)	\n\t"\
			"movaps %%xmm3,%c[__oA](%%rcx)	\n\t"\
			"mulpd	%%xmm2,%%xmm4		\n\t"\
			"mulpd	%%xmm2,%%xmm1		\n\t"\
			"mulpd	%%xmm2,%%xmm3		\n\t"\
			"addpd	%c[__oA](%%rcx),%%xmm4	\n\t"\
			"subpd	%c[__oC](%%rcx),%%xmm1	\n\t"\
			"addpd	%c[__oB](%%rcx),%%xmm3	\n\t"\
			"mulpd	%%xmm0,%%xmm4		\n\t"\
			"mulpd	%%xmm0,%%xmm1		\n\t"\
			"mulpd	%%xmm0,%%xmm3		\n\t"\
			"movaps 0x020(%%rbx),%%xmm2	 \n\t"\
			"addpd	%c[__oB](%%rcx),%%xmm4	\n\t"\
			"subpd	%c[__oA](%%rcx),%%xmm1	\n\t"\
			"subpd	%c[__oC](%%rcx),%%xmm3	\n\t"\
			"mulpd	%%xmm2,%%xmm4		\n\t"\
			"mulpd	%%xmm2,%%xmm1		\n\t"\
			"mulpd	%%xmm2,%%xmm3		\n\t"\
			"movaps -0x010(%%rbx),%%xmm0	\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t"\
			"mulpd	%%xmm0,%%xmm4		\n\t"\
			"subpd	%%xmm1,%%xmm7		\n\t"\
			"mulpd	%%xmm0,%%xmm1		\n\t"\
			"subpd	%%xmm3,%%xmm5		\n\t"\
			"mulpd	%%xmm0,%%xmm3		\n\t"\
			"movaps %%xmm7,%c[__oA](%%rcx)	\n\t"\
			"movaps %%xmm5,%c[__o5](%%rcx)	\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t"\
			"addpd	%%xmm1,%%xmm7		\n\t"\
			"addpd	%%xmm3,%%xmm5		\n\t"\
			"movaps %%xmm5,%c[__oC](%%rcx)	\n\t"\
			"movaps %%xmm7,%c[__oB](%%rcx)	\n\t"\
			"movaps %%xmm6,%c[__o9](%%rcx)	\n\t"\
			"movaps %%xmm4,%c[__o7](%%rcx)	\n\t"\
		"/* yr-terms: */			\n\t"\
			"subq	$0x10,%%rax		 \n\t"\
			"movaps %c[__i1](%%rax),%%xmm1	\n\t"\
			"movaps %c[__i2](%%rax),%%xmm2	\n\t"\
			"movaps %c[__i3](%%rax),%%xmm5	\n\t"\
			"movaps %c[__i4](%%rax),%%xmm3	\n\t"\
			"movaps %c[__i5](%%rax),%%xmm4	\n\t"\
			"movaps %c[__i6](%%rax),%%xmm6	\n\t"\
			"subpd	%c[__iC](%%rax),%%xmm1	\n\t"\
			"subpd	%c[__iB](%%rax),%%xmm2	\n\t"\
			"subpd	%c[__iA](%%rax),%%xmm5	\n\t"\
			"subpd	%c[__i9](%%rax),%%xmm3	\n\t"\
			"subpd	%c[__i8](%%rax),%%xmm4	\n\t"\
			"subpd	%c[__i7](%%rax),%%xmm6	\n\t"\
			"movaps %%xmm1,%%xmm7		\n\t"\
			"movaps %%xmm2,%%xmm0		\n\t"\
			"subpd	%%xmm3,%%xmm7		\n\t"\
			"addpd	%%xmm4,%%xmm0		\n\t"\
			"addpd	%%xmm5,%%xmm7		\n\t"\
			"addpd	%%xmm6,%%xmm0		\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t"\
			"addpd	%%xmm3,%%xmm5		\n\t"\
			"subpd	%%xmm4,%%xmm2		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t"\
			"movaps %%xmm0,%%xmm4		\n\t"\
			"movaps %%xmm7,%%xmm3		\n\t"\
			"mulpd	0x040(%%rbx),%%xmm0	 \n\t"\
			"mulpd	0x040(%%rbx),%%xmm7	 \n\t"\
			"mulpd	0x030(%%rbx),%%xmm4	 \n\t"\
			"mulpd	0x030(%%rbx),%%xmm3	 \n\t"\
			"subpd	%%xmm7,%%xmm4		\n\t"\
			"addpd	%%xmm0,%%xmm3		\n\t"\
			"movaps %%xmm4,%c[__o1](%%rcx)	\n\t"\
			"movaps %%xmm1,%%xmm0		\n\t"\
			"movaps %%xmm5,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t"\
			"movaps %%xmm0,%c[__o2](%%rcx)	\n\t"\
			"movaps %%xmm4,%%xmm7		\n\t"\
			"mulpd	0x0b0(%%rbx),%%xmm4	 \n\t"\
			"mulpd	0x0e0(%%rbx),%%xmm0	 \n\t"\
			"addpd	%c[__o2](%%rcx),%%xmm4	\n\t"\
			"addpd	%%xmm7,%%xmm0		\n\t"\
			"mulpd	0x050(%%rbx),%%xmm4	 \n\t"\
			"mulpd	0x050(%%rbx),%%xmm0	 \n\t"\
			"movaps %%xmm2,%c[__o2](%%rcx)	\n\t"\
			"movaps %%xmm6,%%xmm7		\n\t"\
			"mulpd	0x0c0(%%rbx),%%xmm6	 \n\t"\
			"mulpd	0x0f0(%%rbx),%%xmm2	 \n\t"\
			"addpd	%c[__o2](%%rcx),%%xmm6	\n\t"\
			"addpd	%%xmm7,%%xmm2		\n\t"\
			"mulpd	0x060(%%rbx),%%xmm6	 \n\t"\
			"mulpd	0x060(%%rbx),%%xmm2	 \n\t"\
			"addpd	%%xmm4,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm2		\n\t"\
			"movaps %%xmm1,%c[__o2](%%rcx)	\n\t"\
			"movaps %%xmm5,%%xmm7		\n\t"\
			"mulpd	0x0d0(%%rbx),%%xmm5	 \n\t"\
			"mulpd	0x100(%%rbx),%%xmm1	 \n\t"\
			"addpd	%c[__o2](%%rcx),%%xmm5	\n\t"\
			"addpd	%%xmm7,%%xmm1		\n\t"\
			"mulpd	0x070(%%rbx),%%xmm5	 \n\t"\
			"mulpd	0x070(%%rbx),%%xmm1	 \n\t"\
			"addpd	%%xmm4,%%xmm5		\n\t"\
			"addpd	%%xmm0,%%xmm1		\n\t"\
			"movaps %c[__o1](%%rcx),%%xmm4	\n\t"\
			"movaps %%xmm3,%%xmm7		\n\t"\
			"movaps %%xmm4,%%xmm0		\n\t"\
			"addpd	%%xmm6,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm0		\n\t"\
			"subpd	%%xmm3,%%xmm6		\n\t"\
			"subpd	%%xmm4,%%xmm5		\n\t"\
			"addpd	%%xmm2,%%xmm6		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm2		\n\t"\
			"addpd	%%xmm1,%%xmm4		\n\t"\
			"movaps %c[__o7](%%rcx),%%xmm3	\n\t"\
			"movaps %%xmm3,%%xmm1		\n\t"\
			"subpd	%%xmm7,%%xmm3		\n\t"\
			"addpd	%%xmm7,%%xmm1		\n\t"\
			"movaps %%xmm3,%c[__o7](%%rcx)	\n\t"\
			"movaps %%xmm1,%c[__o6](%%rcx)	\n\t"\
			"movaps %c[__o5](%%rcx),%%xmm3	\n\t"\
			"movaps %%xmm3,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"addpd	%%xmm6,%%xmm1		\n\t"\
			"movaps %%xmm3,%c[__o5](%%rcx)	\n\t"\
			"movaps %%xmm1,%c[__o8](%%rcx)	\n\t"\
			"movaps %c[__oB](%%rcx),%%xmm3	\n\t"\
			"movaps %%xmm3,%%xmm1		\n\t"\
			"subpd	%%xmm2,%%xmm3		\n\t"\
			"addpd	%%xmm2,%%xmm1		\n\t"\
			"movaps %%xmm3,%c[__oB](%%rcx)	\n\t"\
			"movaps %%xmm1,%c[__o2](%%rcx)	\n\t"\
			"movaps %c[__oA](%%rcx),%%xmm3	\n\t"\
			"movaps %%xmm3,%%xmm1		\n\t"\
			"subpd	%%xmm0,%%xmm3		\n\t"\
			"addpd	%%xmm0,%%xmm1		\n\t"\
			"movaps %%xmm3,%c[__oA](%%rcx)	\n\t"\
			"movaps %%xmm1,%c[__o3](%%rcx)	\n\t"\
			"movaps %c[__o9](%%rcx),%%xmm3	\n\t"\
			"movaps %%xmm3,%%xmm1		\n\t"\
			"subpd	%%xmm5,%%xmm3		\n\t"\
			"addpd	%%xmm5,%%xmm1		\n\t"\
			"movaps %%xmm3,%c[__o9](%%rcx)	\n\t"\
			"movaps %%xmm1,%c[__o4](%%rcx)	\n\t"\
			"movaps %c[__oC](%%rcx),%%xmm3	\n\t"\
			"movaps %%xmm3,%%xmm1		\n\t"\
			"subpd	%%xmm4,%%xmm3		\n\t"\
			"addpd	%%xmm4,%%xmm1		\n\t"\
			"movaps %%xmm3,%c[__oC](%%rcx)	\n\t"\
			"movaps %%xmm1,%c[__o1](%%rcx)	\n\t"\
			:					/* outputs: none */\
			: [__cc] "m" (Xcc)	/* All inputs from memory addresses here */\
			 ,[__I0] "m" (XI0)\
			 ,[__i1] "e" (Xi1)\
			 ,[__i2] "e" (Xi2)\
			 ,[__i3] "e" (Xi3)\
			 ,[__i4] "e" (Xi4)\
			 ,[__i5] "e" (Xi5)\
			 ,[__i6] "e" (Xi6)\
			 ,[__i7] "e" (Xi7)\
			 ,[__i8] "e" (Xi8)\
			 ,[__i9] "e" (Xi9)\
			 ,[__iA] "e" (XiA)\
			 ,[__iB] "e" (XiB)\
			 ,[__iC] "e" (XiC)\
			 ,[__O0] "m" (XO0)\
			 ,[__o1] "e" (Xo1)\
			 ,[__o2] "e" (Xo2)\
			 ,[__o3] "e" (Xo3)\
			 ,[__o4] "e" (Xo4)\
			 ,[__o5] "e" (Xo5)\
			 ,[__o6] "e" (Xo6)\
			 ,[__o7] "e" (Xo7)\
			 ,[__o8] "e" (Xo8)\
			 ,[__o9] "e" (Xo9)\
			 ,[__oA] "e" (XoA)\
			 ,[__oB] "e" (XoB)\
			 ,[__oC] "e" (XoC)\
			: "rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	 /* Clobbered registers */\
			);\
		}

		#else // USE_64BIT_ASM_STYLE - 2-Instruction-stream-overlapped 64-bit-ified version of the above 32-bit ASM macro, using all of xmm0-15

		#define SSE2_RADIX_13_DFT(Xcc, XI0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8,Xi9,XiA,XiB,XiC, XO0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8,Xo9,XoA,XoB,XoC)\
		{\
		__asm__ volatile (\
			"movq	%[__I0],%%rax							\n\t"\
			"movq	%[__cc],%%rbx							\n\t"\
			"movq	%[__O0],%%rcx							\n\t"\
		"/* xr-terms need 8 registers for each side: */		\n\t"\
			"movaps	%c[__i6](%%rax),%%xmm7					\n\t"\
			"movaps	%c[__i5](%%rax),%%xmm5					\n\t"\
			"movaps	%c[__i4](%%rax),%%xmm4					\n\t"\
			"movaps	%c[__i3](%%rax),%%xmm6					\n\t		/* yi-terms: */									\n\t"\
			"movaps	%c[__i2](%%rax),%%xmm3					\n\t			movq	%%rax,%%rdx							\n\t"\
			"movaps	%c[__i1](%%rax),%%xmm1					\n\t			addq	$0x10,%%rdx							\n\t"\
			"movaps	-0x010(%%rbx),%%xmm0					\n\t"\
			"addpd	%c[__i7](%%rax),%%xmm7					\n\t			movaps	%c[__i1](%%rdx),%%xmm9				\n\t"\
			"addpd	%c[__i8](%%rax),%%xmm5					\n\t			movaps	%c[__i2](%%rdx),%%xmm10				\n\t"\
			"addpd	%c[__i9](%%rax),%%xmm4					\n\t			movaps	%c[__i3](%%rdx),%%xmm13				\n\t"\
			"addpd	%c[__iA](%%rax),%%xmm6					\n\t			movaps	%c[__i4](%%rdx),%%xmm11				\n\t"\
			"addpd	%c[__iB](%%rax),%%xmm3					\n\t			movaps	%c[__i5](%%rdx),%%xmm12				\n\t"\
			"addpd	%c[__iC](%%rax),%%xmm1					\n\t			movaps	%c[__i6](%%rdx),%%xmm14				\n\t"\
			"subpd	%%xmm5,%%xmm1							\n\t			subpd	%c[__iC](%%rdx),%%xmm9				\n\t"\
			"mulpd	%%xmm0,%%xmm5							\n\t			subpd	%c[__iB](%%rdx),%%xmm10				\n\t"\
			"subpd	%%xmm6,%%xmm3							\n\t			subpd	%c[__iA](%%rdx),%%xmm13				\n\t"\
			"mulpd	%%xmm0,%%xmm6							\n\t			subpd	%c[__i9](%%rdx),%%xmm11				\n\t"\
			"subpd	%%xmm7,%%xmm4							\n\t			subpd	%c[__i8](%%rdx),%%xmm12				\n\t"\
			"mulpd	%%xmm0,%%xmm7							\n\t			subpd	%c[__i7](%%rdx),%%xmm14				\n\t"\
			"addpd	%%xmm1,%%xmm5							\n\t			movaps	%%xmm9 ,%%xmm15						\n\t"\
			"addpd	%%xmm3,%%xmm6							\n\t			movaps	%%xmm10,%%xmm8						\n\t"\
			"addpd	%%xmm4,%%xmm7							\n\t			subpd	%%xmm11,%%xmm15						\n\t"\
			"movaps	%%xmm5,%%xmm2							\n\t			addpd	%%xmm12,%%xmm8						\n\t"\
			"addpd	%%xmm6,%%xmm2							\n\t			addpd	%%xmm13,%%xmm15						\n\t"\
			"addpd	%%xmm7,%%xmm2							\n\t			addpd	%%xmm14,%%xmm8						\n\t"\
			"movaps	(%%rax),%%xmm0							\n\t			addpd	%%xmm11,%%xmm9						\n\t"\
			"addpd	%%xmm2,%%xmm0							\n\t			addpd	%%xmm11,%%xmm13						\n\t"\
			"mulpd	(%%rbx),%%xmm2							\n\t			subpd	%%xmm12,%%xmm10						\n\t"\
			"movaps	%%xmm0,(%%rcx)							\n\t			subpd	%%xmm12,%%xmm14						\n\t"\
			"addpd	%%xmm2,%%xmm0							\n\t			movaps	%%xmm8 ,%%xmm12						\n\t"\
			"movaps	0x010(%%rbx),%%xmm2						\n\t			movaps	%%xmm15,%%xmm11						\n\t"\
			"mulpd	%%xmm2,%%xmm6							\n\t			mulpd	0x040(%%rbx),%%xmm8					\n\t"\
			"mulpd	%%xmm2,%%xmm5							\n\t			mulpd	0x040(%%rbx),%%xmm15				\n\t"\
			"mulpd	%%xmm2,%%xmm7							\n\t			mulpd	0x030(%%rbx),%%xmm12				\n\t"\
			"movaps	%%xmm6,%c[__o1](%%rcx)					\n\t			mulpd	0x030(%%rbx),%%xmm11				\n\t"\
			"movaps	%%xmm5,%c[__o2](%%rcx)					\n\t			subpd	%%xmm15,%%xmm12						\n\t"\
			"movaps	%%xmm7,%c[__o3](%%rcx)					\n\t			addpd	%%xmm8 ,%%xmm11						\n\t"\
			"movaps	0x080(%%rbx),%%xmm2						\n\t			movaps	%%xmm12,%c[__oC](%%rcx)				\n\t"\
			"mulpd	%%xmm2,%%xmm5							\n\t			movaps	%%xmm9 ,%%xmm8						\n\t"\
			"mulpd	%%xmm2,%%xmm7							\n\t			movaps	%%xmm13,%%xmm12						\n\t"\
			"mulpd	%%xmm2,%%xmm6							\n\t			addpd	%%xmm10,%%xmm8						\n\t"\
			"addpd	%%xmm0,%%xmm5							\n\t			addpd	%%xmm14,%%xmm12						\n\t"\
			"addpd	%%xmm0,%%xmm7							\n\t			movaps	%%xmm8 ,%c[__oB](%%rcx)				\n\t"\
			"addpd	%%xmm0,%%xmm6							\n\t			movaps	%%xmm12,%%xmm15						\n\t"\
			"addpd	%c[__o1](%%rcx),%%xmm5					\n\t			mulpd	0x0b0(%%rbx),%%xmm12				\n\t"\
			"addpd	%c[__o2](%%rcx),%%xmm7					\n\t			mulpd	0x0e0(%%rbx),%%xmm8					\n\t"\
			"addpd	%c[__o3](%%rcx),%%xmm6					\n\t			addpd	%c[__oB](%%rcx),%%xmm12				\n\t"\
			"movaps	0x0a0(%%rbx),%%xmm2						\n\t			addpd	%%xmm15,%%xmm8						\n\t"\
			"movaps	0x090(%%rbx),%%xmm0						\n\t			mulpd	0x050(%%rbx),%%xmm12				\n\t"\
			"movaps	%%xmm4,%c[__o1](%%rcx)					\n\t			mulpd	0x050(%%rbx),%%xmm8					\n\t"\
			"movaps	%%xmm1,%c[__o2](%%rcx)					\n\t			movaps	%%xmm10,%c[__oB](%%rcx)				\n\t"\
			"movaps	%%xmm3,%c[__o3](%%rcx)					\n\t			movaps	%%xmm14,%%xmm15						\n\t"\
			"mulpd	%%xmm2,%%xmm4							\n\t			mulpd	0x0c0(%%rbx),%%xmm14				\n\t"\
			"mulpd	%%xmm2,%%xmm1							\n\t			mulpd	0x0f0(%%rbx),%%xmm10				\n\t"\
			"mulpd	%%xmm2,%%xmm3							\n\t			addpd	%c[__oB](%%rcx),%%xmm14				\n\t"\
			"addpd	%c[__o3](%%rcx),%%xmm4					\n\t			addpd	%%xmm15,%%xmm10						\n\t"\
			"subpd	%c[__o1](%%rcx),%%xmm1					\n\t			mulpd	0x060(%%rbx),%%xmm14				\n\t"\
			"addpd	%c[__o2](%%rcx),%%xmm3					\n\t			mulpd	0x060(%%rbx),%%xmm10				\n\t"\
			"mulpd	%%xmm0,%%xmm4							\n\t			addpd	%%xmm12,%%xmm14						\n\t"\
			"mulpd	%%xmm0,%%xmm1							\n\t			addpd	%%xmm8 ,%%xmm10						\n\t"\
			"mulpd	%%xmm0,%%xmm3							\n\t			movaps	%%xmm9 ,%c[__oB](%%rcx)				\n\t"\
			"movaps	0x020(%%rbx),%%xmm2						\n\t			movaps	%%xmm13,%%xmm15						\n\t"\
			"addpd	%c[__o2](%%rcx),%%xmm4					\n\t			mulpd	0x0d0(%%rbx),%%xmm13				\n\t"\
			"subpd	%c[__o3](%%rcx),%%xmm1					\n\t			mulpd	0x100(%%rbx),%%xmm9					\n\t"\
			"subpd	%c[__o1](%%rcx),%%xmm3					\n\t			addpd	%c[__oB](%%rcx),%%xmm13				\n\t"\
			"mulpd	%%xmm2,%%xmm4							\n\t			addpd	%%xmm15,%%xmm9						\n\t"\
			"mulpd	%%xmm2,%%xmm1							\n\t			mulpd	0x070(%%rbx),%%xmm13				\n\t"\
			"mulpd	%%xmm2,%%xmm3							\n\t			mulpd	0x070(%%rbx),%%xmm9					\n\t"\
			"movaps	-0x010(%%rbx),%%xmm0					\n\t			addpd	%%xmm12,%%xmm13						\n\t"\
		"/* Spill into destination outputs: */				\n\t			addpd	%%xmm8 ,%%xmm9						\n\t"\
			"subpd	%%xmm4,%%xmm6							\n\t			movaps	%c[__oC](%%rcx),%%xmm12				\n\t"\
			"mulpd	%%xmm0,%%xmm4							\n\t			movaps	%%xmm11,%%xmm15						\n\t"\
			"subpd	%%xmm1,%%xmm7							\n\t			movaps	%%xmm12,%%xmm8						\n\t"\
			"mulpd	%%xmm0,%%xmm1							\n\t			addpd	%%xmm14,%%xmm15						\n\t"\
			"subpd	%%xmm3,%%xmm5							\n\t			addpd	%%xmm13,%%xmm8						\n\t"\
			"mulpd	%%xmm0,%%xmm3							\n\t			subpd	%%xmm11,%%xmm14						\n\t"\
		"/*	movaps	%%xmm7,%c[__o3](%%rcx)	*/				\n\t			subpd	%%xmm12,%%xmm13						\n\t"\
		"/*	movaps	%%xmm5,%c[__o8](%%rcx)	*/				\n\t			addpd	%%xmm10,%%xmm14						\n\t"\
			"addpd	%%xmm6,%%xmm4							\n\t			addpd	%%xmm9 ,%%xmm13						\n\t"\
			"addpd	%%xmm7,%%xmm1							\n\t			addpd	%%xmm11,%%xmm10						\n\t"\
			"addpd	%%xmm5,%%xmm3							\n\t			addpd	%%xmm9 ,%%xmm12						\n\t"\
		"/*	movaps	%%xmm6,%c[__o4](%%rcx)	*/				\n\t"\
		"/*	movaps	%%xmm4,%c[__o6](%%rcx)	*/				\n\t"\
		"/*	movaps	%%xmm1,%c[__o2](%%rcx)	*/				\n\t"\
		"/*	movaps	%%xmm3,%c[__o1](%%rcx)	*/				\n\t"\
	  "/* yi-data in xmm8,10,12,13,14,15; xmm0,2,9,11 free: */	\n\t"\
		"/*	movaps	%c[__o6](%%rcx),%%xmm4	*/				\n\t"\
			"movaps	%%xmm4 ,%%xmm0							\n\t"\
			"subpd	%%xmm15,%%xmm4							\n\t		/*	movaps	%c[__o3](%%rcx),%%xmm7	*/			\n\t"\
			"addpd	%%xmm15,%%xmm0							\n\t			movaps	%%xmm7 ,%%xmm9						\n\t"\
			"movaps	%%xmm4 ,%c[__o6](%%rcx)					\n\t			subpd	%%xmm8 ,%%xmm7						\n\t"\
			"movaps	%%xmm0 ,%c[__o7](%%rcx)					\n\t			addpd	%%xmm8 ,%%xmm9						\n\t"\
		"/*	movaps	%c[__o8](%%rcx),%%xmm5	*/				\n\t			movaps	%%xmm7 ,%c[__o3](%%rcx)				\n\t"\
			"movaps	%%xmm5 ,%%xmm0							\n\t			movaps	%%xmm9 ,%c[__oA](%%rcx)				\n\t"\
			"subpd	%%xmm14,%%xmm5							\n\t		/*	movaps	%c[__o4](%%rcx),%%xmm6	*/			\n\t"\
			"addpd	%%xmm14,%%xmm0							\n\t			movaps	%%xmm6 ,%%xmm9						\n\t"\
			"movaps	%%xmm5 ,%c[__o8](%%rcx)					\n\t			subpd	%%xmm13,%%xmm6						\n\t"\
			"movaps	%%xmm0 ,%c[__o5](%%rcx)					\n\t			addpd	%%xmm13,%%xmm9						\n\t"\
		"/*	movaps	%c[__o2](%%rcx),%%xmm1	*/				\n\t			movaps	%%xmm6 ,%c[__o4](%%rcx)				\n\t"\
			"movaps	%%xmm1 ,%%xmm0							\n\t			movaps	%%xmm9 ,%c[__o9](%%rcx)				\n\t"\
			"subpd	%%xmm10,%%xmm1							\n\t		/*	movaps	%c[__o1](%%rcx),%%xmm3	*/			\n\t"\
			"addpd	%%xmm10,%%xmm0							\n\t			movaps	%%xmm3 ,%%xmm9						\n\t"\
			"movaps	%%xmm1 ,%c[__o2](%%rcx)					\n\t			subpd	%%xmm12,%%xmm3						\n\t"\
			"movaps	%%xmm0 ,%c[__oB](%%rcx)					\n\t			addpd	%%xmm12,%%xmm9						\n\t"\
			"/***************/								\n\t			movaps	%%xmm3 ,%c[__o1](%%rcx)				\n\t"\
			"/* IMAG PARTS: */								\n\t			movaps	%%xmm9 ,%c[__oC](%%rcx)				\n\t"\
			"/***************/								\n\t"\
		"/* xi-terms need 8 registers for each	side: */	\n\t"\
			"addq	$0x10,%%rax		 						\n\t"\
			"addq	$0x10,%%rcx								\n\t"\
			"movaps	%c[__i6](%%rax),%%xmm7					\n\t"\
			"movaps	%c[__i5](%%rax),%%xmm5					\n\t"\
			"movaps	%c[__i4](%%rax),%%xmm4					\n\t"\
			"movaps	%c[__i3](%%rax),%%xmm6					\n\t		/* yi-terms: */									\n\t"\
			"movaps	%c[__i2](%%rax),%%xmm3					\n\t"\
			"movaps	%c[__i1](%%rax),%%xmm1					\n\t			subq	$0x10,%%rdx							\n\t"\
			"movaps	-0x010(%%rbx),%%xmm0					\n\t"\
			"addpd	%c[__i7](%%rax),%%xmm7					\n\t			movaps	%c[__i1](%%rdx),%%xmm9				\n\t"\
			"addpd	%c[__i8](%%rax),%%xmm5					\n\t			movaps	%c[__i2](%%rdx),%%xmm10				\n\t"\
			"addpd	%c[__i9](%%rax),%%xmm4					\n\t			movaps	%c[__i3](%%rdx),%%xmm13				\n\t"\
			"addpd	%c[__iA](%%rax),%%xmm6					\n\t			movaps	%c[__i4](%%rdx),%%xmm11				\n\t"\
			"addpd	%c[__iB](%%rax),%%xmm3					\n\t			movaps	%c[__i5](%%rdx),%%xmm12				\n\t"\
			"addpd	%c[__iC](%%rax),%%xmm1					\n\t			movaps	%c[__i6](%%rdx),%%xmm14				\n\t"\
			"subpd	%%xmm5,%%xmm1							\n\t			subpd	%c[__iC](%%rdx),%%xmm9				\n\t"\
			"mulpd	%%xmm0,%%xmm5							\n\t			subpd	%c[__iB](%%rdx),%%xmm10				\n\t"\
			"subpd	%%xmm6,%%xmm3							\n\t			subpd	%c[__iA](%%rdx),%%xmm13				\n\t"\
			"mulpd	%%xmm0,%%xmm6							\n\t			subpd	%c[__i9](%%rdx),%%xmm11				\n\t"\
			"subpd	%%xmm7,%%xmm4							\n\t			subpd	%c[__i8](%%rdx),%%xmm12				\n\t"\
			"mulpd	%%xmm0,%%xmm7							\n\t			subpd	%c[__i7](%%rdx),%%xmm14				\n\t"\
			"addpd	%%xmm1,%%xmm5							\n\t			movaps	%%xmm9 ,%%xmm15						\n\t"\
			"addpd	%%xmm3,%%xmm6							\n\t			movaps	%%xmm10,%%xmm8						\n\t"\
			"addpd	%%xmm4,%%xmm7							\n\t			subpd	%%xmm11,%%xmm15						\n\t"\
			"movaps	%%xmm5,%%xmm2							\n\t			addpd	%%xmm12,%%xmm8						\n\t"\
			"addpd	%%xmm6,%%xmm2							\n\t			addpd	%%xmm13,%%xmm15						\n\t"\
			"addpd	%%xmm7,%%xmm2							\n\t			addpd	%%xmm14,%%xmm8						\n\t"\
			"movaps	(%%rax),%%xmm0							\n\t			addpd	%%xmm11,%%xmm9						\n\t"\
			"addpd	%%xmm2,%%xmm0							\n\t			addpd	%%xmm11,%%xmm13						\n\t"\
			"mulpd	(%%rbx),%%xmm2							\n\t			subpd	%%xmm12,%%xmm10						\n\t"\
			"movaps	%%xmm0,(%%rcx)							\n\t			subpd	%%xmm12,%%xmm14						\n\t"\
			"addpd	%%xmm2,%%xmm0							\n\t			movaps	%%xmm8 ,%%xmm12						\n\t"\
			"movaps	0x010(%%rbx),%%xmm2						\n\t			movaps	%%xmm15,%%xmm11						\n\t"\
			"mulpd	%%xmm2,%%xmm6							\n\t			mulpd	0x040(%%rbx),%%xmm8					\n\t"\
			"mulpd	%%xmm2,%%xmm5							\n\t			mulpd	0x040(%%rbx),%%xmm15				\n\t"\
			"mulpd	%%xmm2,%%xmm7							\n\t			mulpd	0x030(%%rbx),%%xmm12				\n\t"\
			"movaps	%%xmm6,%c[__oC](%%rcx)					\n\t			mulpd	0x030(%%rbx),%%xmm11				\n\t"\
			"movaps	%%xmm5,%c[__oB](%%rcx)					\n\t			subpd	%%xmm15,%%xmm12						\n\t"\
			"movaps	%%xmm7,%c[__oA](%%rcx)					\n\t			addpd	%%xmm8 ,%%xmm11						\n\t"\
			"movaps	0x080(%%rbx),%%xmm2						\n\t			movaps	%%xmm12,%c[__o1](%%rcx)				\n\t"\
			"mulpd	%%xmm2,%%xmm5							\n\t			movaps	%%xmm9 ,%%xmm8						\n\t"\
			"mulpd	%%xmm2,%%xmm7							\n\t			movaps	%%xmm13,%%xmm12						\n\t"\
			"mulpd	%%xmm2,%%xmm6							\n\t			addpd	%%xmm10,%%xmm8						\n\t"\
			"addpd	%%xmm0,%%xmm5							\n\t			addpd	%%xmm14,%%xmm12						\n\t"\
			"addpd	%%xmm0,%%xmm7							\n\t			movaps	%%xmm8 ,%c[__o2](%%rcx)				\n\t"\
			"addpd	%%xmm0,%%xmm6							\n\t			movaps	%%xmm12,%%xmm15						\n\t"\
			"addpd	%c[__oC](%%rcx),%%xmm5					\n\t			mulpd	0x0b0(%%rbx),%%xmm12				\n\t"\
			"addpd	%c[__oB](%%rcx),%%xmm7					\n\t			mulpd	0x0e0(%%rbx),%%xmm8					\n\t"\
			"addpd	%c[__oA](%%rcx),%%xmm6					\n\t			addpd	%c[__o2](%%rcx),%%xmm12				\n\t"\
			"movaps	0x0a0(%%rbx),%%xmm2						\n\t			addpd	%%xmm15,%%xmm8						\n\t"\
			"movaps	0x090(%%rbx),%%xmm0						\n\t			mulpd	0x050(%%rbx),%%xmm12				\n\t"\
			"movaps	%%xmm4,%c[__oC](%%rcx)					\n\t			mulpd	0x050(%%rbx),%%xmm8					\n\t"\
			"movaps	%%xmm1,%c[__oB](%%rcx)					\n\t			movaps	%%xmm10,%c[__o2](%%rcx)				\n\t"\
			"movaps	%%xmm3,%c[__oA](%%rcx)					\n\t			movaps	%%xmm14,%%xmm15						\n\t"\
			"mulpd	%%xmm2,%%xmm4							\n\t			mulpd	0x0c0(%%rbx),%%xmm14				\n\t"\
			"mulpd	%%xmm2,%%xmm1							\n\t			mulpd	0x0f0(%%rbx),%%xmm10				\n\t"\
			"mulpd	%%xmm2,%%xmm3							\n\t			addpd	%c[__o2](%%rcx),%%xmm14				\n\t"\
			"addpd	%c[__oA](%%rcx),%%xmm4					\n\t			addpd	%%xmm15,%%xmm10						\n\t"\
			"subpd	%c[__oC](%%rcx),%%xmm1					\n\t			mulpd	0x060(%%rbx),%%xmm14				\n\t"\
			"addpd	%c[__oB](%%rcx),%%xmm3					\n\t			mulpd	0x060(%%rbx),%%xmm10				\n\t"\
			"mulpd	%%xmm0,%%xmm4							\n\t			addpd	%%xmm12,%%xmm14						\n\t"\
			"mulpd	%%xmm0,%%xmm1							\n\t			addpd	%%xmm8 ,%%xmm10						\n\t"\
			"mulpd	%%xmm0,%%xmm3							\n\t			movaps	%%xmm9 ,%c[__o2](%%rcx)				\n\t"\
			"movaps	0x020(%%rbx),%%xmm2						\n\t			movaps	%%xmm13,%%xmm15						\n\t"\
			"addpd	%c[__oB](%%rcx),%%xmm4					\n\t			mulpd	0x0d0(%%rbx),%%xmm13				\n\t"\
			"subpd	%c[__oA](%%rcx),%%xmm1					\n\t			mulpd	0x100(%%rbx),%%xmm9					\n\t"\
			"subpd	%c[__oC](%%rcx),%%xmm3					\n\t			addpd	%c[__o2](%%rcx),%%xmm13				\n\t"\
			"mulpd	%%xmm2,%%xmm4							\n\t			addpd	%%xmm15,%%xmm9						\n\t"\
			"mulpd	%%xmm2,%%xmm1							\n\t			mulpd	0x070(%%rbx),%%xmm13				\n\t"\
			"mulpd	%%xmm2,%%xmm3							\n\t			mulpd	0x070(%%rbx),%%xmm9					\n\t"\
			"movaps	-0x010(%%rbx),%%xmm0					\n\t			addpd	%%xmm12,%%xmm13						\n\t"\
		"/* Spill into destination outputs: */				\n\t			addpd	%%xmm8 ,%%xmm9						\n\t"\
			"subpd	%%xmm4,%%xmm6							\n\t			movaps	%c[__o1](%%rcx),%%xmm12				\n\t"\
			"mulpd	%%xmm0,%%xmm4							\n\t			movaps	%%xmm11,%%xmm15						\n\t"\
			"subpd	%%xmm1,%%xmm7							\n\t			movaps	%%xmm12,%%xmm8						\n\t"\
			"mulpd	%%xmm0,%%xmm1							\n\t			addpd	%%xmm14,%%xmm15						\n\t"\
			"subpd	%%xmm3,%%xmm5							\n\t			addpd	%%xmm13,%%xmm8						\n\t"\
			"mulpd	%%xmm0,%%xmm3							\n\t			subpd	%%xmm11,%%xmm14						\n\t"\
		"/*	movaps	%%xmm7,%c[__oA](%%rcx)	*/				\n\t			subpd	%%xmm12,%%xmm13						\n\t"\
		"/*	movaps	%%xmm5,%c[__o5](%%rcx)	*/				\n\t			addpd	%%xmm10,%%xmm14						\n\t"\
			"addpd	%%xmm6,%%xmm4							\n\t			addpd	%%xmm9 ,%%xmm13						\n\t"\
			"addpd	%%xmm7,%%xmm1							\n\t			addpd	%%xmm11,%%xmm10						\n\t"\
			"addpd	%%xmm5,%%xmm3							\n\t			addpd	%%xmm9 ,%%xmm12						\n\t"\
		"/*	movaps	%%xmm6,%c[__o9](%%rcx)	*/				\n\t"\
		"/*	movaps	%%xmm4,%c[__o7](%%rcx)	*/				\n\t"\
		"/*	movaps	%%xmm1,%c[__oB](%%rcx)	*/				\n\t"\
		"/*	movaps	%%xmm3,%c[__oC](%%rcx)	*/				\n\t"\
	  "/* yi-data in xmm8,10,12,13,14,15; xmm0,2,9,11 free: */	\n\t"\
		"/*	movaps	%c[__o7](%%rcx),%%xmm4	*/				\n\t"\
			"movaps	%%xmm4 ,%%xmm0							\n\t"\
			"subpd	%%xmm15,%%xmm4							\n\t		/*	movaps	%c[__oA](%%rcx),%%xmm7	*/			\n\t"\
			"addpd	%%xmm15,%%xmm0							\n\t			movaps	%%xmm7 ,%%xmm9						\n\t"\
			"movaps	%%xmm4 ,%c[__o7](%%rcx)					\n\t			subpd	%%xmm8 ,%%xmm7						\n\t"\
			"movaps	%%xmm0 ,%c[__o6](%%rcx)					\n\t			addpd	%%xmm8 ,%%xmm9						\n\t"\
		"/*	movaps	%c[__o5](%%rcx),%%xmm5	*/				\n\t			movaps	%%xmm7 ,%c[__oA](%%rcx)				\n\t"\
			"movaps	%%xmm5 ,%%xmm0							\n\t			movaps	%%xmm9 ,%c[__o3](%%rcx)				\n\t"\
			"subpd	%%xmm14,%%xmm5							\n\t		/*	movaps	%c[__o9](%%rcx),%%xmm6	*/			\n\t"\
			"addpd	%%xmm14,%%xmm0							\n\t			movaps	%%xmm6 ,%%xmm9						\n\t"\
			"movaps	%%xmm5 ,%c[__o5](%%rcx)					\n\t			subpd	%%xmm13,%%xmm6						\n\t"\
			"movaps	%%xmm0 ,%c[__o8](%%rcx)					\n\t			addpd	%%xmm13,%%xmm9						\n\t"\
		"/*	movaps	%c[__oB](%%rcx),%%xmm1	*/				\n\t			movaps	%%xmm6 ,%c[__o9](%%rcx)				\n\t"\
			"movaps	%%xmm1 ,%%xmm0							\n\t			movaps	%%xmm9 ,%c[__o4](%%rcx)				\n\t"\
			"subpd	%%xmm10,%%xmm1							\n\t		/*	movaps	%c[__oC](%%rcx),%%xmm3	*/			\n\t"\
			"addpd	%%xmm10,%%xmm0							\n\t			movaps	%%xmm3 ,%%xmm9						\n\t"\
			"movaps	%%xmm1 ,%c[__oB](%%rcx)					\n\t			subpd	%%xmm12,%%xmm3						\n\t"\
			"movaps	%%xmm0 ,%c[__o2](%%rcx)					\n\t			addpd	%%xmm12,%%xmm9						\n\t"\
			"																movaps	%%xmm3 ,%c[__oC](%%rcx)				\n\t"\
			"																movaps	%%xmm9 ,%c[__o1](%%rcx)				\n\t"\
			:					/* outputs: none */\
			: [__cc] "m" (Xcc)	/* All inputs from memory addresses here */\
			 ,[__I0] "m" (XI0)\
			 ,[__i1] "e" (Xi1)\
			 ,[__i2] "e" (Xi2)\
			 ,[__i3] "e" (Xi3)\
			 ,[__i4] "e" (Xi4)\
			 ,[__i5] "e" (Xi5)\
			 ,[__i6] "e" (Xi6)\
			 ,[__i7] "e" (Xi7)\
			 ,[__i8] "e" (Xi8)\
			 ,[__i9] "e" (Xi9)\
			 ,[__iA] "e" (XiA)\
			 ,[__iB] "e" (XiB)\
			 ,[__iC] "e" (XiC)\
			 ,[__O0] "m" (XO0)\
			 ,[__o1] "e" (Xo1)\
			 ,[__o2] "e" (Xo2)\
			 ,[__o3] "e" (Xo3)\
			 ,[__o4] "e" (Xo4)\
			 ,[__o5] "e" (Xo5)\
			 ,[__o6] "e" (Xo6)\
			 ,[__o7] "e" (Xo7)\
			 ,[__o8] "e" (Xo8)\
			 ,[__o9] "e" (Xo9)\
			 ,[__oA] "e" (XoA)\
			 ,[__oB] "e" (XoB)\
			 ,[__oC] "e" (XoC)\
			: "rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
			);\
		}

		#endif // USE_32BIT_ASM_STYLE

	  #endif	// IF(GCC), USE 32/64-BIT ASM STYLE

	#endif	// MSVC / GCC

	#ifdef COMPILER_TYPE_GCC

		#if OS_BITS == 32

		//	#include "radix52_ditN_cy_dif1_gcc32.h"

		#else

			#include "radix52_ditN_cy_dif1_gcc64.h"

		#endif

	#endif

#endif	// SSE2

/**************/

int radix52_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-52 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-52 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	const uint32 RADIX = 52;
	const double crnd = 3.0*0x4000000*0x2000000;
	int NDIVR,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40,p44,p48;
	static double radix_inv, n2inv;
	double scale, maxerr = 0.0;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1,
	t00r,t01r,t02r,t03r,t04r,t05r,t06r,t07r,t08r,t09r,t0ar,t0br,t0cr,
	t10r,t11r,t12r,t13r,t14r,t15r,t16r,t17r,t18r,t19r,t1ar,t1br,t1cr,
	t20r,t21r,t22r,t23r,t24r,t25r,t26r,t27r,t28r,t29r,t2ar,t2br,t2cr,
	t30r,t31r,t32r,t33r,t34r,t35r,t36r,t37r,t38r,t39r,t3ar,t3br,t3cr;

#ifdef USE_SSE2

  #if !(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))
	#error SSE2 code not supported for this compiler!
  #endif

	static int cslots_in_local_store;
	static struct complex *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD

	#ifdef USE_PTHREAD
		static struct complex *__r0;	// Base address for discrete per-thread local stores
		static struct cy_thread_data_t *tdat = 0x0;
		// Threadpool-based dispatch stuff:
		static int main_work_units = 0, pool_work_units = 0;
		static struct threadpool *tpool = 0x0;
		static int task_is_blocking = TRUE;
		static thread_control_t thread_control = {0,0,0};
		// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
		static task_control_t   task_control = {NULL, (void*)cy52_process_chunk, NULL, 0x0};
	#endif

  #else
	double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
  #endif

	static int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09
		,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19
		,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27,*bjmodn28,*bjmodn29
		,*bjmodn30,*bjmodn31,*bjmodn32,*bjmodn33,*bjmodn34,*bjmodn35,*bjmodn36,*bjmodn37,*bjmodn38,*bjmodn39
		,*bjmodn40,*bjmodn41,*bjmodn42,*bjmodn43,*bjmodn44,*bjmodn45,*bjmodn46,*bjmodn47,*bjmodn48,*bjmodn49
		,*bjmodn50,*bjmodn51;
	static struct complex *two,*rad13_const, *max_err, *sse2_rnd, *half_arr, *tmp     /* rad13_const needs 18*16 bytes allocated */
		,*r00r,*r00i,*r01r,*r01i,*r02r,*r02i,*r03r,*r03i,*r04r,*r04i,*r05r,*r05i,*r06r,*r06i,*r07r,*r07i,*r08r,*r08i,*r09r,*r09i,*r0ar,*r0ai,*r0br,*r0bi,*r0cr,*r0ci
		,*r10r,*r10i,*r11r,*r11i,*r12r,*r12i,*r13r,*r13i,*r14r,*r14i,*r15r,*r15i,*r16r,*r16i,*r17r,*r17i,*r18r,*r18i,*r19r,*r19i,*r1ar,*r1ai,*r1br,*r1bi,*r1cr,*r1ci
		,*r20r,*r20i,*r21r,*r21i,*r22r,*r22i,*r23r,*r23i,*r24r,*r24i,*r25r,*r25i,*r26r,*r26i,*r27r,*r27i,*r28r,*r28i,*r29r,*r29i,*r2ar,*r2ai,*r2br,*r2bi,*r2cr,*r2ci
		,*r30r,*r30i,*r31r,*r31i,*r32r,*r32i,*r33r,*r33i,*r34r,*r34i,*r35r,*r35i,*r36r,*r36i,*r37r,*r37i,*r38r,*r38i,*r39r,*r39i,*r3ar,*r3ai,*r3br,*r3bi,*r3cr,*r3ci
		,*s1p00r,*s1p00i,*s1p01r,*s1p01i,*s1p02r,*s1p02i,*s1p03r,*s1p03i,*s1p04r,*s1p04i,*s1p05r,*s1p05i,*s1p06r,*s1p06i,*s1p07r,*s1p07i,*s1p08r,*s1p08i,*s1p09r,*s1p09i,*s1p0ar,*s1p0ai,*s1p0br,*s1p0bi,*s1p0cr,*s1p0ci
		,*s1p10r,*s1p10i,*s1p11r,*s1p11i,*s1p12r,*s1p12i,*s1p13r,*s1p13i,*s1p14r,*s1p14i,*s1p15r,*s1p15i,*s1p16r,*s1p16i,*s1p17r,*s1p17i,*s1p18r,*s1p18i,*s1p19r,*s1p19i,*s1p1ar,*s1p1ai,*s1p1br,*s1p1bi,*s1p1cr,*s1p1ci
		,*s1p20r,*s1p20i,*s1p21r,*s1p21i,*s1p22r,*s1p22i,*s1p23r,*s1p23i,*s1p24r,*s1p24i,*s1p25r,*s1p25i,*s1p26r,*s1p26i,*s1p27r,*s1p27i,*s1p28r,*s1p28i,*s1p29r,*s1p29i,*s1p2ar,*s1p2ai,*s1p2br,*s1p2bi,*s1p2cr,*s1p2ci
		,*s1p30r,*s1p30i,*s1p31r,*s1p31i,*s1p32r,*s1p32i,*s1p33r,*s1p33i,*s1p34r,*s1p34i,*s1p35r,*s1p35i,*s1p36r,*s1p36i,*s1p37r,*s1p37i,*s1p38r,*s1p38i,*s1p39r,*s1p39i,*s1p3ar,*s1p3ai,*s1p3br,*s1p3bi,*s1p3cr,*s1p3ci
		,*cy00,*cy02,*cy04,*cy06,*cy08
		,*cy10,*cy12,*cy14,*cy16,*cy18
		,*cy20,*cy22,*cy24,*cy26,*cy28
		,*cy30,*cy32,*cy34,*cy36,*cy38
		,*cy40,*cy42,*cy44,*cy46,*cy48,*cy50;

#else

	int m,m2,
	bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,
	bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,
	bjmodn26,bjmodn27,bjmodn28,bjmodn29,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,bjmodn36,bjmodn37,bjmodn38,
	bjmodn39,bjmodn40,bjmodn41,bjmodn42,bjmodn43,bjmodn44,bjmodn45,bjmodn46,bjmodn47,bjmodn48,bjmodn49,bjmodn50,bjmodn51;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	double rt,it,temp,frac;
	double wt,wtinv,wtA,wtB,wtC,
	t00i,t01i,t02i,t03i,t04i,t05i,t06i,t07i,t08i,t09i,t0ai,t0bi,t0ci,
	t10i,t11i,t12i,t13i,t14i,t15i,t16i,t17i,t18i,t19i,t1ai,t1bi,t1ci,
	t20i,t21i,t22i,t23i,t24i,t25i,t26i,t27i,t28i,t29i,t2ai,t2bi,t2ci,
	t30i,t31i,t32i,t33i,t34i,t35i,t36i,t37i,t38i,t39i,t3ai,t3bi,t3ci,
	a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,a1p08r,a1p08i,a1p09r,a1p09i,a1p0ar,a1p0ai,a1p0br,a1p0bi,a1p0cr,a1p0ci,
	a1p10r,a1p10i,a1p11r,a1p11i,a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a1p1ar,a1p1ai,a1p1br,a1p1bi,a1p1cr,a1p1ci,
	a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,a1p24r,a1p24i,a1p25r,a1p25i,a1p26r,a1p26i,a1p27r,a1p27i,a1p28r,a1p28i,a1p29r,a1p29i,a1p2ar,a1p2ai,a1p2br,a1p2bi,a1p2cr,a1p2ci,
	a1p30r,a1p30i,a1p31r,a1p31i,a1p32r,a1p32i,a1p33r,a1p33i,a1p34r,a1p34i,a1p35r,a1p35i,a1p36r,a1p36i,a1p37r,a1p37i,a1p38r,a1p38i,a1p39r,a1p39i,a1p3ar,a1p3ai,a1p3br,a1p3bi,a1p3cr,a1p3ci,
	cy00,cy01,cy02,cy03,cy04,cy05,cy06,cy07,cy08,cy09,cy10,cy11,cy12,
	cy13,cy14,cy15,cy16,cy17,cy18,cy19,cy20,cy21,cy22,cy23,cy24,cy25,
	cy26,cy27,cy28,cy29,cy30,cy31,cy32,cy33,cy34,cy35,cy36,cy37,cy38,
	cy39,cy40,cy41,cy42,cy43,cy44,cy45,cy46,cy47,cy48,cy49,cy50,cy51;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_bjmodnini = 0x0,
	*_bjmodn00 = 0x0,*_bjmodn01 = 0x0,*_bjmodn02 = 0x0,*_bjmodn03 = 0x0,*_bjmodn04 = 0x0,*_bjmodn05 = 0x0,*_bjmodn06 = 0x0,*_bjmodn07 = 0x0,*_bjmodn08 = 0x0,*_bjmodn09 = 0x0,
	*_bjmodn10 = 0x0,*_bjmodn11 = 0x0,*_bjmodn12 = 0x0,*_bjmodn13 = 0x0,*_bjmodn14 = 0x0,*_bjmodn15 = 0x0,*_bjmodn16 = 0x0,*_bjmodn17 = 0x0,*_bjmodn18 = 0x0,*_bjmodn19 = 0x0,
	*_bjmodn20 = 0x0,*_bjmodn21 = 0x0,*_bjmodn22 = 0x0,*_bjmodn23 = 0x0,*_bjmodn24 = 0x0,*_bjmodn25 = 0x0,*_bjmodn26 = 0x0,*_bjmodn27 = 0x0,*_bjmodn28 = 0x0,*_bjmodn29 = 0x0,
	*_bjmodn30 = 0x0,*_bjmodn31 = 0x0,*_bjmodn32 = 0x0,*_bjmodn33 = 0x0,*_bjmodn34 = 0x0,*_bjmodn35 = 0x0,*_bjmodn36 = 0x0,*_bjmodn37 = 0x0,*_bjmodn38 = 0x0,*_bjmodn39 = 0x0,
	*_bjmodn40 = 0x0,*_bjmodn41 = 0x0,*_bjmodn42 = 0x0,*_bjmodn43 = 0x0,*_bjmodn44 = 0x0,*_bjmodn45 = 0x0,*_bjmodn46 = 0x0,*_bjmodn47 = 0x0,*_bjmodn48 = 0x0,*_bjmodn49 = 0x0,
	*_bjmodn50 = 0x0,*_bjmodn51 = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_maxerr = 0x0,
	*_cy00 = 0x0,*_cy01 = 0x0,*_cy02 = 0x0,*_cy03 = 0x0,*_cy04 = 0x0,*_cy05 = 0x0,*_cy06 = 0x0,*_cy07 = 0x0,*_cy08 = 0x0,*_cy09 = 0x0,
	*_cy10 = 0x0,*_cy11 = 0x0,*_cy12 = 0x0,*_cy13 = 0x0,*_cy14 = 0x0,*_cy15 = 0x0,*_cy16 = 0x0,*_cy17 = 0x0,*_cy18 = 0x0,*_cy19 = 0x0,
	*_cy20 = 0x0,*_cy21 = 0x0,*_cy22 = 0x0,*_cy23 = 0x0,*_cy24 = 0x0,*_cy25 = 0x0,*_cy26 = 0x0,*_cy27 = 0x0,*_cy28 = 0x0,*_cy29 = 0x0,
	*_cy30 = 0x0,*_cy31 = 0x0,*_cy32 = 0x0,*_cy33 = 0x0,*_cy34 = 0x0,*_cy35 = 0x0,*_cy36 = 0x0,*_cy37 = 0x0,*_cy38 = 0x0,*_cy39 = 0x0,
	*_cy40 = 0x0,*_cy41 = 0x0,*_cy42 = 0x0,*_cy43 = 0x0,*_cy44 = 0x0,*_cy45 = 0x0,*_cy46 = 0x0,*_cy47 = 0x0,*_cy48 = 0x0,*_cy49 = 0x0,
	*_cy50 = 0x0,*_cy51 = 0x0;

	if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
	{
		ASSERT(HERE, 0, "radix24_ditN_cy_dif1: Fermat-mod only available for radices 7,8,9,15 and their multiples!");
	}

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/RADIX;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/RADIX in radix24_ditN_cy_dif1.\n",iter);
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

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
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
		if(0 != (j & 0xf)) {
			printf("sizeof(cy_thread_data_t) = %x\n",j);
			ASSERT(HERE, 0, "struct cy_thread_data_t not 16-byte size multiple!");
		}
		tdat = (struct cy_thread_data_t *)calloc(CY_THREADS, sizeof(struct cy_thread_data_t));

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

	#ifdef USE_SSE2

		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use double-complex type size (16 bytes) to alloc a block of local storage
		// consisting of 88 dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix52_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_COMPLEX(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		sm_ptr = (uint64*)(sc_ptr + radix52_creals_in_local_store);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 48 16-byte slots of sc_arr for temporaries, next 2 for the doubled cos and c3m1 terms,
	next 52/2 = 26 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
	#ifdef USE_PTHREAD
		__r0 = sc_ptr;
	#endif
									two = sc_ptr + 0x1a;	max_err  = two + 0x1a;		sse2_rnd = max_err + 0x1a;	/* Use these as handy temps */
		s1p00r = sc_ptr + 0x00;    s1p10r = two + 0x00;    s1p20r = max_err + 0x00;    s1p30r = sse2_rnd + 0x00;
		s1p00i = sc_ptr + 0x01;    s1p10i = two + 0x01;    s1p20i = max_err + 0x01;    s1p30i = sse2_rnd + 0x01;
		s1p01r = sc_ptr + 0x02;    s1p11r = two + 0x02;    s1p21r = max_err + 0x02;    s1p31r = sse2_rnd + 0x02;
		s1p01i = sc_ptr + 0x03;    s1p11i = two + 0x03;    s1p21i = max_err + 0x03;    s1p31i = sse2_rnd + 0x03;
		s1p02r = sc_ptr + 0x04;    s1p12r = two + 0x04;    s1p22r = max_err + 0x04;    s1p32r = sse2_rnd + 0x04;
		s1p02i = sc_ptr + 0x05;    s1p12i = two + 0x05;    s1p22i = max_err + 0x05;    s1p32i = sse2_rnd + 0x05;
		s1p03r = sc_ptr + 0x06;    s1p13r = two + 0x06;    s1p23r = max_err + 0x06;    s1p33r = sse2_rnd + 0x06;
		s1p03i = sc_ptr + 0x07;    s1p13i = two + 0x07;    s1p23i = max_err + 0x07;    s1p33i = sse2_rnd + 0x07;
		s1p04r = sc_ptr + 0x08;    s1p14r = two + 0x08;    s1p24r = max_err + 0x08;    s1p34r = sse2_rnd + 0x08;
		s1p04i = sc_ptr + 0x09;    s1p14i = two + 0x09;    s1p24i = max_err + 0x09;    s1p34i = sse2_rnd + 0x09;
		s1p05r = sc_ptr + 0x0a;    s1p15r = two + 0x0a;    s1p25r = max_err + 0x0a;    s1p35r = sse2_rnd + 0x0a;
		s1p05i = sc_ptr + 0x0b;    s1p15i = two + 0x0b;    s1p25i = max_err + 0x0b;    s1p35i = sse2_rnd + 0x0b;
		s1p06r = sc_ptr + 0x0c;    s1p16r = two + 0x0c;    s1p26r = max_err + 0x0c;    s1p36r = sse2_rnd + 0x0c;
		s1p06i = sc_ptr + 0x0d;    s1p16i = two + 0x0d;    s1p26i = max_err + 0x0d;    s1p36i = sse2_rnd + 0x0d;
		s1p07r = sc_ptr + 0x0e;    s1p17r = two + 0x0e;    s1p27r = max_err + 0x0e;    s1p37r = sse2_rnd + 0x0e;
		s1p07i = sc_ptr + 0x0f;    s1p17i = two + 0x0f;    s1p27i = max_err + 0x0f;    s1p37i = sse2_rnd + 0x0f;
		s1p08r = sc_ptr + 0x10;    s1p18r = two + 0x10;    s1p28r = max_err + 0x10;    s1p38r = sse2_rnd + 0x10;
		s1p08i = sc_ptr + 0x11;    s1p18i = two + 0x11;    s1p28i = max_err + 0x11;    s1p38i = sse2_rnd + 0x11;
		s1p09r = sc_ptr + 0x12;    s1p19r = two + 0x12;    s1p29r = max_err + 0x12;    s1p39r = sse2_rnd + 0x12;
		s1p09i = sc_ptr + 0x13;    s1p19i = two + 0x13;    s1p29i = max_err + 0x13;    s1p39i = sse2_rnd + 0x13;
		s1p0ar = sc_ptr + 0x14;    s1p1ar = two + 0x14;    s1p2ar = max_err + 0x14;    s1p3ar = sse2_rnd + 0x14;
		s1p0ai = sc_ptr + 0x15;    s1p1ai = two + 0x15;    s1p2ai = max_err + 0x15;    s1p3ai = sse2_rnd + 0x15;
		s1p0br = sc_ptr + 0x16;    s1p1br = two + 0x16;    s1p2br = max_err + 0x16;    s1p3br = sse2_rnd + 0x16;
		s1p0bi = sc_ptr + 0x17;    s1p1bi = two + 0x17;    s1p2bi = max_err + 0x17;    s1p3bi = sse2_rnd + 0x17;
		s1p0cr = sc_ptr + 0x18;    s1p1cr = two + 0x18;    s1p2cr = max_err + 0x18;    s1p3cr = sse2_rnd + 0x18;
		s1p0ci = sc_ptr + 0x19;    s1p1ci = two + 0x19;    s1p2ci = max_err + 0x19;    s1p3ci = sse2_rnd + 0x19;

		tmp = sse2_rnd + 0x1a; two = tmp + 0x1a;    max_err = two + 0x1a;   sse2_rnd = max_err + 0x1a;	/* Use these as handy temps */
		r00r = tmp + 0x00;    r10r = two + 0x00;    r20r = max_err + 0x00;    r30r = sse2_rnd + 0x00;
		r00i = tmp + 0x01;    r10i = two + 0x01;    r20i = max_err + 0x01;    r30i = sse2_rnd + 0x01;
		r01r = tmp + 0x02;    r11r = two + 0x02;    r21r = max_err + 0x02;    r31r = sse2_rnd + 0x02;
		r01i = tmp + 0x03;    r11i = two + 0x03;    r21i = max_err + 0x03;    r31i = sse2_rnd + 0x03;
		r02r = tmp + 0x04;    r12r = two + 0x04;    r22r = max_err + 0x04;    r32r = sse2_rnd + 0x04;
		r02i = tmp + 0x05;    r12i = two + 0x05;    r22i = max_err + 0x05;    r32i = sse2_rnd + 0x05;
		r03r = tmp + 0x06;    r13r = two + 0x06;    r23r = max_err + 0x06;    r33r = sse2_rnd + 0x06;
		r03i = tmp + 0x07;    r13i = two + 0x07;    r23i = max_err + 0x07;    r33i = sse2_rnd + 0x07;
		r04r = tmp + 0x08;    r14r = two + 0x08;    r24r = max_err + 0x08;    r34r = sse2_rnd + 0x08;
		r04i = tmp + 0x09;    r14i = two + 0x09;    r24i = max_err + 0x09;    r34i = sse2_rnd + 0x09;
		r05r = tmp + 0x0a;    r15r = two + 0x0a;    r25r = max_err + 0x0a;    r35r = sse2_rnd + 0x0a;
		r05i = tmp + 0x0b;    r15i = two + 0x0b;    r25i = max_err + 0x0b;    r35i = sse2_rnd + 0x0b;
		r06r = tmp + 0x0c;    r16r = two + 0x0c;    r26r = max_err + 0x0c;    r36r = sse2_rnd + 0x0c;
		r06i = tmp + 0x0d;    r16i = two + 0x0d;    r26i = max_err + 0x0d;    r36i = sse2_rnd + 0x0d;
		r07r = tmp + 0x0e;    r17r = two + 0x0e;    r27r = max_err + 0x0e;    r37r = sse2_rnd + 0x0e;
		r07i = tmp + 0x0f;    r17i = two + 0x0f;    r27i = max_err + 0x0f;    r37i = sse2_rnd + 0x0f;
		r08r = tmp + 0x10;    r18r = two + 0x10;    r28r = max_err + 0x10;    r38r = sse2_rnd + 0x10;
		r08i = tmp + 0x11;    r18i = two + 0x11;    r28i = max_err + 0x11;    r38i = sse2_rnd + 0x11;
		r09r = tmp + 0x12;    r19r = two + 0x12;    r29r = max_err + 0x12;    r39r = sse2_rnd + 0x12;
		r09i = tmp + 0x13;    r19i = two + 0x13;    r29i = max_err + 0x13;    r39i = sse2_rnd + 0x13;
		r0ar = tmp + 0x14;    r1ar = two + 0x14;    r2ar = max_err + 0x14;    r3ar = sse2_rnd + 0x14;
		r0ai = tmp + 0x15;    r1ai = two + 0x15;    r2ai = max_err + 0x15;    r3ai = sse2_rnd + 0x15;
		r0br = tmp + 0x16;    r1br = two + 0x16;    r2br = max_err + 0x16;    r3br = sse2_rnd + 0x16;
		r0bi = tmp + 0x17;    r1bi = two + 0x17;    r2bi = max_err + 0x17;    r3bi = sse2_rnd + 0x17;
		r0cr = tmp + 0x18;    r1cr = two + 0x18;    r2cr = max_err + 0x18;    r3cr = sse2_rnd + 0x18;
		r0ci = tmp + 0x19;    r1ci = two + 0x19;    r2ci = max_err + 0x19;    r3ci = sse2_rnd + 0x19;
		tmp = sse2_rnd + 0x1a;

		rad13_const = tmp + 0x01;	/* Leave an extra slot at radix13_const-1 for the constant two = 2.0: */
		tmp += 0x14;	/* Need 20 16-byte slots for two+sincos, but offset the carry slots by the next-larger multiple of 4 */

		cy00 = tmp + 0x00;
		cy02 = tmp + 0x01;
		cy04 = tmp + 0x02;
		cy06 = tmp + 0x03;
		cy08 = tmp + 0x04;
		cy10 = tmp + 0x05;
		cy12 = tmp + 0x06;
		cy14 = tmp + 0x07;
		cy16 = tmp + 0x08;
		cy18 = tmp + 0x09;
		cy20 = tmp + 0x0a;
		cy22 = tmp + 0x0b;
		cy24 = tmp + 0x0c;
		cy26 = tmp + 0x0d;
		cy28 = tmp + 0x0e;
		cy30 = tmp + 0x0f;
		cy32 = tmp + 0x10;
		cy34 = tmp + 0x11;
		cy36 = tmp + 0x12;
		cy38 = tmp + 0x13;
		cy40 = tmp + 0x14;
		cy42 = tmp + 0x15;
		cy44 = tmp + 0x16;
		cy46 = tmp + 0x17;
		cy48 = tmp + 0x18;
		cy50 = tmp + 0x19;
		tmp += 0x1a;

		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;	/* This table needs 20x16 bytes for Mersenne-mod, and radixx16 for Fermat-mod */
		ASSERT(HERE, (radix52_creals_in_local_store << 4) >= ((long)half_arr - (long)s1p00r) + (20 << 4), "radix52_creals_in_local_store checksum failed!");
		/* These remain fixed: */
		tmp = rad13_const-1;							/* __cc pointer offsets: */
		tmp->re = tmp->im =  2.0;		++tmp;	/*	-0x010 = 2.0 */
		tmp->re = tmp->im =  DC1;		++tmp;	/*	0x000 =  DC1 */
		tmp->re = tmp->im =  DC3;		++tmp;	/*	0x010 =  DC3 */
		tmp->re = tmp->im =  DC4;		++tmp;	/*	0x020 =  DC4 */
		tmp->re = tmp->im =  DS1;		++tmp;	/*	0x030 =  DS1 */
		tmp->re = tmp->im =  DS2;		++tmp;	/*	0x040 =  DS2 */
		tmp->re = tmp->im =  DS3;		++tmp;	/*	0x050 =  DS3 */
		tmp->re = tmp->im =  DS4;		++tmp;	/*	0x060 =  DS4 */
		tmp->re = tmp->im =  DS5;		++tmp;	/*	0x070 =  DS5 */
		tmp->re = tmp->im = DC23;		++tmp;	/*	0x080 = DC23 */
		tmp->re = tmp->im = DC54;		++tmp;	/*	0x090 = DC54 */
		tmp->re = tmp->im = DC65;		++tmp;	/*	0x0a0 = DC65 */
		tmp->re = tmp->im = DS63;		++tmp;	/*	0x0b0 = DS63 */
		tmp->re = tmp->im = DS74;		++tmp;	/*	0x0c0 = DS74 */
		tmp->re = tmp->im = DS85;		++tmp;	/*	0x0d0 = DS85 */
		tmp->re = tmp->im = DS93;		++tmp;	/*	0x0e0 = DS93 */
		tmp->re = tmp->im = DSa4;		++tmp;	/*	0x0f0 = DSa4 */
		tmp->re = tmp->im = DSb5;		++tmp;	/*	0x100 = DSb5 */

		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
		sse2_rnd->re = sse2_rnd->im = crnd;

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
		*/
		tmp = half_arr;
		/* Forward-weight multipliers: */
		tmp->re = 1.0;	tmp->im = 1.0;	++tmp;
		tmp->re = .50;	tmp->im = 1.0;	++tmp;
		tmp->re = 1.0;	tmp->im = .50;	++tmp;
		tmp->re = .50;	tmp->im = .50;	++tmp;
		/* Inverse-weight multipliers: */
		tmp->re = .50;	tmp->im = .50;	++tmp;
		tmp->re = .25;	tmp->im = .50;	++tmp;
		tmp->re = .50;	tmp->im = .25;	++tmp;
		tmp->re = .25;	tmp->im = .25;	++tmp;
		/* Forward-base[] multipliers: */
		tmp->re = base   [0];	tmp->im = base   [0];	++tmp;
		tmp->re = base   [1];	tmp->im = base   [0];	++tmp;
		tmp->re = base   [0];	tmp->im = base   [1];	++tmp;
		tmp->re = base   [1];	tmp->im = base   [1];	++tmp;
		/* Inverse-base[] multipliers: */
		tmp->re = baseinv[0];	tmp->im = baseinv[0];	++tmp;
		tmp->re = baseinv[1];	tmp->im = baseinv[0];	++tmp;
		tmp->re = baseinv[0];	tmp->im = baseinv[1];	++tmp;
		tmp->re = baseinv[1];	tmp->im = baseinv[1];	++tmp;

		/* Floating-point sign mask used for FABS on packed doubles: */
		sign_mask = sm_ptr;
		*sign_mask++ = (uint64)0x7FFFFFFFFFFFFFFFull;
		*sign_mask-- = (uint64)0x7FFFFFFFFFFFFFFFull;

		sse_bw  = sm_ptr + 2;
		tmp64 = (uint64)bw;
		tmp64 = tmp64 + (tmp64 << 32);
		*sse_bw++ = tmp64;
		*sse_bw-- = tmp64;

		sse_sw  = sm_ptr + 4;
		tmp64 = (uint64)sw;
		tmp64 = tmp64 + (tmp64 << 32);
		*sse_sw++ = tmp64;
		*sse_sw-- = tmp64;

		sse_n   = sm_ptr + 6;
		tmp64 = (uint64)n;
		tmp64 = tmp64 + (tmp64 << 32);
		*sse_n++ = tmp64;
		*sse_n-- = tmp64;

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
		tdat[ithread].s1p00r = __r0 + ithread*cslots_in_local_store;
		tdat[ithread].half_arr = (long)tdat[ithread].s1p00r + ((long)half_arr - (long)s1p00r);
	}
#endif

		bjmodn00 = (uint32*)(sm_ptr + 8);
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

	  #ifdef USE_PTHREAD
		tmp = __r0 + cslots_in_local_store;
		/* Init thread 1-CY_THREADS's local stores and pointers: */
		for(i = 1; i < CY_THREADS; ++i) {
			/* Only care about the constants for each thread here, but easier to just copy the entire thread0 local store: */
			memcpy(tmp, __r0, cslots_in_local_store<<4);	// bytewise copy treats complex and uint64 subdata the same
			tmp += cslots_in_local_store;
		}
	  #endif

	#endif	// USE_SSE2

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

		if(_cy00)	/* If it's a new exponent of a range test, need to deallocate these. */
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

			free((void *)_cy00); _cy00 = 0x0;
			free((void *)_cy01); _cy01 = 0x0;
			free((void *)_cy02); _cy02 = 0x0;
			free((void *)_cy03); _cy03 = 0x0;
			free((void *)_cy04); _cy04 = 0x0;
			free((void *)_cy05); _cy05 = 0x0;
			free((void *)_cy06); _cy06 = 0x0;
			free((void *)_cy07); _cy07 = 0x0;
			free((void *)_cy08); _cy08 = 0x0;
			free((void *)_cy09); _cy09 = 0x0;
			free((void *)_cy10); _cy10 = 0x0;
			free((void *)_cy11); _cy11 = 0x0;
			free((void *)_cy12); _cy12 = 0x0;
			free((void *)_cy13); _cy13 = 0x0;
			free((void *)_cy14); _cy14 = 0x0;
			free((void *)_cy15); _cy15 = 0x0;
			free((void *)_cy16); _cy16 = 0x0;
			free((void *)_cy17); _cy17 = 0x0;
			free((void *)_cy18); _cy18 = 0x0;
			free((void *)_cy19); _cy19 = 0x0;
			free((void *)_cy20); _cy20 = 0x0;
			free((void *)_cy21); _cy21 = 0x0;
			free((void *)_cy22); _cy22 = 0x0;
			free((void *)_cy23); _cy23 = 0x0;
			free((void *)_cy24); _cy24 = 0x0;
			free((void *)_cy25); _cy25 = 0x0;
			free((void *)_cy26); _cy26 = 0x0;
			free((void *)_cy27); _cy27 = 0x0;
			free((void *)_cy28); _cy28 = 0x0;
			free((void *)_cy29); _cy29 = 0x0;
			free((void *)_cy30); _cy30 = 0x0;
			free((void *)_cy31); _cy31 = 0x0;
			free((void *)_cy32); _cy32 = 0x0;
			free((void *)_cy33); _cy33 = 0x0;
			free((void *)_cy34); _cy34 = 0x0;
			free((void *)_cy35); _cy35 = 0x0;
			free((void *)_cy36); _cy36 = 0x0;
			free((void *)_cy37); _cy37 = 0x0;
			free((void *)_cy38); _cy38 = 0x0;
			free((void *)_cy39); _cy39 = 0x0;
			free((void *)_cy40); _cy40 = 0x0;
			free((void *)_cy41); _cy41 = 0x0;
			free((void *)_cy42); _cy42 = 0x0;
			free((void *)_cy43); _cy43 = 0x0;
			free((void *)_cy44); _cy44 = 0x0;
			free((void *)_cy45); _cy45 = 0x0;
			free((void *)_cy46); _cy46 = 0x0;
			free((void *)_cy47); _cy47 = 0x0;
			free((void *)_cy48); _cy48 = 0x0;
			free((void *)_cy49); _cy49 = 0x0;
			free((void *)_cy50); _cy50 = 0x0;
			free((void *)_cy51); _cy51 = 0x0;

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
		_jstart  	= (int *)malloc(j);	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(j);	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(j);	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(j);	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(j);	ptr_prod += (uint32)(_co3     == 0x0);

		j = CY_THREADS*sizeof(double);
		_cy00	= (double *)malloc(j);	ptr_prod += (uint32)(_cy00== 0x0);
		_cy01	= (double *)malloc(j);	ptr_prod += (uint32)(_cy01== 0x0);
		_cy02	= (double *)malloc(j);	ptr_prod += (uint32)(_cy02== 0x0);
		_cy03	= (double *)malloc(j);	ptr_prod += (uint32)(_cy03== 0x0);
		_cy04	= (double *)malloc(j);	ptr_prod += (uint32)(_cy04== 0x0);
		_cy05	= (double *)malloc(j);	ptr_prod += (uint32)(_cy05== 0x0);
		_cy06	= (double *)malloc(j);	ptr_prod += (uint32)(_cy06== 0x0);
		_cy07	= (double *)malloc(j);	ptr_prod += (uint32)(_cy07== 0x0);
		_cy08	= (double *)malloc(j);	ptr_prod += (uint32)(_cy08== 0x0);
		_cy09	= (double *)malloc(j);	ptr_prod += (uint32)(_cy09== 0x0);
		_cy10	= (double *)malloc(j);	ptr_prod += (uint32)(_cy10== 0x0);
		_cy11	= (double *)malloc(j);	ptr_prod += (uint32)(_cy11== 0x0);
		_cy12	= (double *)malloc(j);	ptr_prod += (uint32)(_cy12== 0x0);
		_cy13	= (double *)malloc(j);	ptr_prod += (uint32)(_cy13== 0x0);
		_cy14	= (double *)malloc(j);	ptr_prod += (uint32)(_cy14== 0x0);
		_cy15	= (double *)malloc(j);	ptr_prod += (uint32)(_cy15== 0x0);
		_cy16	= (double *)malloc(j);	ptr_prod += (uint32)(_cy16== 0x0);
		_cy17	= (double *)malloc(j);	ptr_prod += (uint32)(_cy17== 0x0);
		_cy18	= (double *)malloc(j);	ptr_prod += (uint32)(_cy18== 0x0);
		_cy19	= (double *)malloc(j);	ptr_prod += (uint32)(_cy19== 0x0);
		_cy20	= (double *)malloc(j);	ptr_prod += (uint32)(_cy20== 0x0);
		_cy21	= (double *)malloc(j);	ptr_prod += (uint32)(_cy21== 0x0);
		_cy22	= (double *)malloc(j);	ptr_prod += (uint32)(_cy22== 0x0);
		_cy23	= (double *)malloc(j);	ptr_prod += (uint32)(_cy23== 0x0);
		_cy24	= (double *)malloc(j);	ptr_prod += (uint32)(_cy24== 0x0);
		_cy25	= (double *)malloc(j);	ptr_prod += (uint32)(_cy25== 0x0);
		_cy26	= (double *)malloc(j);	ptr_prod += (uint32)(_cy26== 0x0);
		_cy27	= (double *)malloc(j);	ptr_prod += (uint32)(_cy27== 0x0);
		_cy28	= (double *)malloc(j);	ptr_prod += (uint32)(_cy28== 0x0);
		_cy29	= (double *)malloc(j);	ptr_prod += (uint32)(_cy29== 0x0);
		_cy30	= (double *)malloc(j);	ptr_prod += (uint32)(_cy30== 0x0);
		_cy31	= (double *)malloc(j);	ptr_prod += (uint32)(_cy31== 0x0);
		_cy32	= (double *)malloc(j);	ptr_prod += (uint32)(_cy32== 0x0);
		_cy33	= (double *)malloc(j);	ptr_prod += (uint32)(_cy33== 0x0);
		_cy34	= (double *)malloc(j);	ptr_prod += (uint32)(_cy34== 0x0);
		_cy35	= (double *)malloc(j);	ptr_prod += (uint32)(_cy35== 0x0);
		_cy36	= (double *)malloc(j);	ptr_prod += (uint32)(_cy36== 0x0);
		_cy37	= (double *)malloc(j);	ptr_prod += (uint32)(_cy37== 0x0);
		_cy38	= (double *)malloc(j);	ptr_prod += (uint32)(_cy38== 0x0);
		_cy39	= (double *)malloc(j);	ptr_prod += (uint32)(_cy39== 0x0);
		_cy40	= (double *)malloc(j);	ptr_prod += (uint32)(_cy40== 0x0);
		_cy41	= (double *)malloc(j);	ptr_prod += (uint32)(_cy41== 0x0);
		_cy42	= (double *)malloc(j);	ptr_prod += (uint32)(_cy42== 0x0);
		_cy43	= (double *)malloc(j);	ptr_prod += (uint32)(_cy43== 0x0);
		_cy44	= (double *)malloc(j);	ptr_prod += (uint32)(_cy44== 0x0);
		_cy45	= (double *)malloc(j);	ptr_prod += (uint32)(_cy45== 0x0);
		_cy46	= (double *)malloc(j);	ptr_prod += (uint32)(_cy46== 0x0);
		_cy47	= (double *)malloc(j);	ptr_prod += (uint32)(_cy47== 0x0);
		_cy48	= (double *)malloc(j);	ptr_prod += (uint32)(_cy48== 0x0);
		_cy49	= (double *)malloc(j);	ptr_prod += (uint32)(_cy49== 0x0);
		_cy50	= (double *)malloc(j);	ptr_prod += (uint32)(_cy50== 0x0);
		_cy51	= (double *)malloc(j);	ptr_prod += (uint32)(_cy51== 0x0);

		_maxerr	= (double *)malloc(j);	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays in radix52_ditN_cy_dif1.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/52-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini in radix52_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		_bjmodnini[0] = 0;
		_bjmodnini[1] = 0;

		/* For Fermat-mod, since 'adjacent' words are actually stride-2 separated
		in terms of the floating residue array, block boundaries have half the i-index
		(e.g. as in sw*i and bw*i) value they do in the Mersenne-mod case:
		*/
		jhi = NDIVR/CY_THREADS;

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

/*...The radix-52 final DIT pass is here.	*/

	/* init carries	*/
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_cy00[ithread] = 0;
		_cy01[ithread] = 0;
		_cy02[ithread] = 0;
		_cy03[ithread] = 0;
		_cy04[ithread] = 0;
		_cy05[ithread] = 0;
		_cy06[ithread] = 0;
		_cy07[ithread] = 0;
		_cy08[ithread] = 0;
		_cy09[ithread] = 0;
		_cy10[ithread] = 0;
		_cy11[ithread] = 0;
		_cy12[ithread] = 0;
		_cy13[ithread] = 0;
		_cy14[ithread] = 0;
		_cy15[ithread] = 0;
		_cy16[ithread] = 0;
		_cy17[ithread] = 0;
		_cy18[ithread] = 0;
		_cy19[ithread] = 0;
		_cy20[ithread] = 0;
		_cy21[ithread] = 0;
		_cy22[ithread] = 0;
		_cy23[ithread] = 0;
		_cy24[ithread] = 0;
		_cy25[ithread] = 0;
		_cy26[ithread] = 0;
		_cy27[ithread] = 0;
		_cy28[ithread] = 0;
		_cy29[ithread] = 0;
		_cy30[ithread] = 0;
		_cy31[ithread] = 0;
		_cy32[ithread] = 0;
		_cy33[ithread] = 0;
		_cy34[ithread] = 0;
		_cy35[ithread] = 0;
		_cy36[ithread] = 0;
		_cy37[ithread] = 0;
		_cy38[ithread] = 0;
		_cy39[ithread] = 0;
		_cy40[ithread] = 0;
		_cy41[ithread] = 0;
		_cy42[ithread] = 0;
		_cy43[ithread] = 0;
		_cy44[ithread] = 0;
		_cy45[ithread] = 0;
		_cy46[ithread] = 0;
		_cy47[ithread] = 0;
		_cy48[ithread] = 0;
		_cy49[ithread] = 0;
		_cy50[ithread] = 0;
		_cy51[ithread] = 0;
	}
	/* If an LL test, init the subtract-2: */
	if(TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		_cy00[      0] = -2;
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
	_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. If n does not divide p, lowest-order digit is always a bigword (_i[0] = 1).	*/

	if(CY_THREADS > 1)
	{
		for(ithread = 1; ithread < CY_THREADS; ithread++)
		{
			_i[ithread] = ((uint32)(sw - _bjmodnini[ithread]) >> 31);
		}
	}

	/*
	Moved this inside the outer-loop, so on cleanup pass can use it to reset _col,_co2,_co3 starting values,
	then simply overwrite it with 1 prior to starting the k-loop.
	*/
	khi = n_div_nwt/CY_THREADS;
	j = _bjmodnini[CY_THREADS];
	// Include 0-thread here ... bjmodn terms all 0 for that, but need jhi computed for all threads:
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

		_jstart[ithread] = ithread*NDIVR/CY_THREADS;
		if(!full_pass)
			_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
		else
			_jhi[ithread] = _jstart[ithread] + nwt-1;

		_col[ithread] = ithread*(khi*RADIX);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
		_co2[ithread] = (n>>nwt_bits)-1+RADIX - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
		_co3[ithread] = _co2[ithread]-RADIX;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
	}

    /* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
    if(!full_pass)
    {
		khi = 1;
	#if FFT_DEBUG
		fprintf(dbg_file, "radix52_ditN_cy_dif1: Cleanup Pass:\n", n);
	#endif
	}

/* Needed to remove the prefetch-address vars add0 & add for this to compile properly: */
#ifdef USE_OMP
	#error radix52_ditN_cy_dif1: OpenMP private list needs updating!
	omp_set_num_threads(CY_THREADS);
//#undef PFETCH
	#pragma omp parallel for private(temp,frac,maxerr,i,j,j1,j2,jt,jp,jstart,jhi,k,khi,l,col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,wtl,wtlp1,wtn,wtnm1,wt,wtinv,wtA,wtB,wtC,wt_re,wt_im,rt,it,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0a,t0b,t0c,t0d,t0e,t0f,t0g,t0h,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1a,t1b,t1c,t1d,t1e,t1f,t1g,t1h,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2a,t2b,t2c,t2d,t2e,t2f,t2g,t2h,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3a,t3b,t3c,t3d,t3e,t3f,t3g,t3h,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r,a1p28r,a1p29r,a1p30r,a1p31r,a1p32r,a1p33r,a1p34r,a1p35r,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i,a1p28i,a1p29i,a1p30i,a1p31i,a1p32i,a1p33i,a1p34i,a1p35i,bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,bjmodn28,bjmodn29,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,cy00,cy01,cy02,cy03,cy04,cy05,cy06,cy07,cy08,cy09,cy10,cy11,cy12,cy13,cy14,cy15,cy16,cy17,cy18,cy19,cy20,cy21,cy22,cy23,cy24,cy25,cy26,cy27,cy28,cy29,cy30,cy31,cy32,cy33,cy34,cy35,cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27,cy_i28,cy_i29,cy_i30,cy_i31,cy_i32,cy_i33,cy_i34,cy_i35) default(shared) schedule(static)
*******************to-do update private-data list!*******************
#endif

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
		ASSERT(HERE, tdat[ithread].s1p00r == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].half_arr;
		ASSERT(HERE, ((tmp-1)->re == crnd && (tmp-1)->im == crnd), "thread-local memcheck failed!");
		ASSERT(HERE, (tmp+10)->re * (tmp+14)->re == 1.0 && (tmp+10)->im * (tmp+14)->im == 1.0, "thread-local memcheck failed!");

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
		/* init carries	*/
		tdat[ithread].cy00 = _cy00[ithread];
		tdat[ithread].cy01 = _cy01[ithread];
		tdat[ithread].cy02 = _cy02[ithread];
		tdat[ithread].cy03 = _cy03[ithread];
		tdat[ithread].cy04 = _cy04[ithread];
		tdat[ithread].cy05 = _cy05[ithread];
		tdat[ithread].cy06 = _cy06[ithread];
		tdat[ithread].cy07 = _cy07[ithread];
		tdat[ithread].cy08 = _cy08[ithread];
		tdat[ithread].cy09 = _cy09[ithread];
		tdat[ithread].cy10 = _cy10[ithread];
		tdat[ithread].cy11 = _cy11[ithread];
		tdat[ithread].cy12 = _cy12[ithread];
		tdat[ithread].cy13 = _cy13[ithread];
		tdat[ithread].cy14 = _cy14[ithread];
		tdat[ithread].cy15 = _cy15[ithread];
		tdat[ithread].cy16 = _cy16[ithread];
		tdat[ithread].cy17 = _cy17[ithread];
		tdat[ithread].cy18 = _cy18[ithread];
		tdat[ithread].cy19 = _cy19[ithread];
		tdat[ithread].cy20 = _cy20[ithread];
		tdat[ithread].cy21 = _cy21[ithread];
		tdat[ithread].cy22 = _cy22[ithread];
		tdat[ithread].cy23 = _cy23[ithread];
		tdat[ithread].cy24 = _cy24[ithread];
		tdat[ithread].cy25 = _cy25[ithread];
		tdat[ithread].cy26 = _cy26[ithread];
		tdat[ithread].cy27 = _cy27[ithread];
		tdat[ithread].cy28 = _cy28[ithread];
		tdat[ithread].cy29 = _cy29[ithread];
		tdat[ithread].cy30 = _cy30[ithread];
		tdat[ithread].cy31 = _cy31[ithread];
		tdat[ithread].cy32 = _cy32[ithread];
		tdat[ithread].cy33 = _cy33[ithread];
		tdat[ithread].cy34 = _cy34[ithread];
		tdat[ithread].cy35 = _cy35[ithread];
		tdat[ithread].cy36 = _cy36[ithread];
		tdat[ithread].cy37 = _cy37[ithread];
		tdat[ithread].cy38 = _cy38[ithread];
		tdat[ithread].cy39 = _cy39[ithread];
		tdat[ithread].cy40 = _cy40[ithread];
		tdat[ithread].cy41 = _cy41[ithread];
		tdat[ithread].cy42 = _cy42[ithread];
		tdat[ithread].cy43 = _cy43[ithread];
		tdat[ithread].cy44 = _cy44[ithread];
		tdat[ithread].cy45 = _cy45[ithread];
		tdat[ithread].cy46 = _cy46[ithread];
		tdat[ithread].cy47 = _cy47[ithread];
		tdat[ithread].cy48 = _cy48[ithread];
		tdat[ithread].cy49 = _cy49[ithread];
		tdat[ithread].cy50 = _cy50[ithread];
		tdat[ithread].cy51 = _cy51[ithread];
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
		max_err->re = 0.0;	max_err->im = 0.0;
	#endif

		i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		jstart = _jstart[ithread];
		jhi    = _jhi[ithread];

		col = _col[ithread];
		co2 = _co2[ithread];
		co3 = _co3[ithread];

	#ifdef USE_SSE2								/* init carries	*/
		*bjmodn00 = _bjmodn00[ithread];			cy00->re = _cy00[ithread];
		*bjmodn01 = _bjmodn01[ithread];			cy00->im = _cy01[ithread];
		*bjmodn02 = _bjmodn02[ithread];			cy02->re = _cy02[ithread];
		*bjmodn03 = _bjmodn03[ithread];			cy02->im = _cy03[ithread];
		*bjmodn04 = _bjmodn04[ithread];			cy04->re = _cy04[ithread];
		*bjmodn05 = _bjmodn05[ithread];			cy04->im = _cy05[ithread];
		*bjmodn06 = _bjmodn06[ithread];			cy06->re = _cy06[ithread];
		*bjmodn07 = _bjmodn07[ithread];			cy06->im = _cy07[ithread];
		*bjmodn08 = _bjmodn08[ithread];			cy08->re = _cy08[ithread];
		*bjmodn09 = _bjmodn09[ithread];			cy08->im = _cy09[ithread];
		*bjmodn10 = _bjmodn10[ithread];			cy10->re = _cy10[ithread];
		*bjmodn11 = _bjmodn11[ithread];			cy10->im = _cy11[ithread];
		*bjmodn12 = _bjmodn12[ithread];			cy12->re = _cy12[ithread];
		*bjmodn13 = _bjmodn13[ithread];			cy12->im = _cy13[ithread];
		*bjmodn14 = _bjmodn14[ithread];			cy14->re = _cy14[ithread];
		*bjmodn15 = _bjmodn15[ithread];			cy14->im = _cy15[ithread];
		*bjmodn16 = _bjmodn16[ithread];			cy16->re = _cy16[ithread];
		*bjmodn17 = _bjmodn17[ithread];			cy16->im = _cy17[ithread];
		*bjmodn18 = _bjmodn18[ithread];			cy18->re = _cy18[ithread];
		*bjmodn19 = _bjmodn19[ithread];			cy18->im = _cy19[ithread];
		*bjmodn20 = _bjmodn20[ithread];			cy20->re = _cy20[ithread];
		*bjmodn21 = _bjmodn21[ithread];			cy20->im = _cy21[ithread];
		*bjmodn22 = _bjmodn22[ithread];			cy22->re = _cy22[ithread];
		*bjmodn23 = _bjmodn23[ithread];			cy22->im = _cy23[ithread];
		*bjmodn24 = _bjmodn24[ithread];			cy24->re = _cy24[ithread];
		*bjmodn25 = _bjmodn25[ithread];			cy24->im = _cy25[ithread];
		*bjmodn26 = _bjmodn26[ithread];			cy26->re = _cy26[ithread];
		*bjmodn27 = _bjmodn27[ithread];			cy26->im = _cy27[ithread];
		*bjmodn28 = _bjmodn28[ithread];			cy28->re = _cy28[ithread];
		*bjmodn29 = _bjmodn29[ithread];			cy28->im = _cy29[ithread];
		*bjmodn30 = _bjmodn30[ithread];			cy30->re = _cy30[ithread];
		*bjmodn31 = _bjmodn31[ithread];			cy30->im = _cy31[ithread];
		*bjmodn32 = _bjmodn32[ithread];			cy32->re = _cy32[ithread];
		*bjmodn33 = _bjmodn33[ithread];			cy32->im = _cy33[ithread];
		*bjmodn34 = _bjmodn34[ithread];			cy34->re = _cy34[ithread];
		*bjmodn35 = _bjmodn35[ithread];			cy34->im = _cy35[ithread];
		*bjmodn36 = _bjmodn36[ithread];			cy36->re = _cy36[ithread];
		*bjmodn37 = _bjmodn37[ithread];			cy36->im = _cy37[ithread];
		*bjmodn38 = _bjmodn38[ithread];			cy38->re = _cy38[ithread];
		*bjmodn39 = _bjmodn39[ithread];			cy38->im = _cy39[ithread];
		*bjmodn40 = _bjmodn40[ithread];			cy40->re = _cy40[ithread];
		*bjmodn41 = _bjmodn41[ithread];			cy40->im = _cy41[ithread];
		*bjmodn42 = _bjmodn42[ithread];			cy42->re = _cy42[ithread];
		*bjmodn43 = _bjmodn43[ithread];			cy42->im = _cy43[ithread];
		*bjmodn44 = _bjmodn44[ithread];			cy44->re = _cy44[ithread];
		*bjmodn45 = _bjmodn45[ithread];			cy44->im = _cy45[ithread];
		*bjmodn46 = _bjmodn46[ithread];			cy46->re = _cy46[ithread];
		*bjmodn47 = _bjmodn47[ithread];			cy46->im = _cy47[ithread];
		*bjmodn48 = _bjmodn48[ithread];			cy48->re = _cy48[ithread];
		*bjmodn49 = _bjmodn49[ithread];			cy48->im = _cy49[ithread];
		*bjmodn50 = _bjmodn50[ithread];			cy50->re = _cy50[ithread];
		*bjmodn51 = _bjmodn51[ithread];			cy50->im = _cy51[ithread];
	#else
		bjmodn00 = _bjmodn00[ithread];			cy00 = _cy00[ithread];
		bjmodn01 = _bjmodn01[ithread];			cy01 = _cy01[ithread];
		bjmodn02 = _bjmodn02[ithread];			cy02 = _cy02[ithread];
		bjmodn03 = _bjmodn03[ithread];			cy03 = _cy03[ithread];
		bjmodn04 = _bjmodn04[ithread];			cy04 = _cy04[ithread];
		bjmodn05 = _bjmodn05[ithread];			cy05 = _cy05[ithread];
		bjmodn06 = _bjmodn06[ithread];			cy06 = _cy06[ithread];
		bjmodn07 = _bjmodn07[ithread];			cy07 = _cy07[ithread];
		bjmodn08 = _bjmodn08[ithread];			cy08 = _cy08[ithread];
		bjmodn09 = _bjmodn09[ithread];			cy09 = _cy09[ithread];
		bjmodn10 = _bjmodn10[ithread];			cy10 = _cy10[ithread];
		bjmodn11 = _bjmodn11[ithread];			cy11 = _cy11[ithread];
		bjmodn12 = _bjmodn12[ithread];			cy12 = _cy12[ithread];
		bjmodn13 = _bjmodn13[ithread];			cy13 = _cy13[ithread];
		bjmodn14 = _bjmodn14[ithread];			cy14 = _cy14[ithread];
		bjmodn15 = _bjmodn15[ithread];			cy15 = _cy15[ithread];
		bjmodn16 = _bjmodn16[ithread];			cy16 = _cy16[ithread];
		bjmodn17 = _bjmodn17[ithread];			cy17 = _cy17[ithread];
		bjmodn18 = _bjmodn18[ithread];			cy18 = _cy18[ithread];
		bjmodn19 = _bjmodn19[ithread];			cy19 = _cy19[ithread];
		bjmodn20 = _bjmodn20[ithread];			cy20 = _cy20[ithread];
		bjmodn21 = _bjmodn21[ithread];			cy21 = _cy21[ithread];
		bjmodn22 = _bjmodn22[ithread];			cy22 = _cy22[ithread];
		bjmodn23 = _bjmodn23[ithread];			cy23 = _cy23[ithread];
		bjmodn24 = _bjmodn24[ithread];			cy24 = _cy24[ithread];
		bjmodn25 = _bjmodn25[ithread];			cy25 = _cy25[ithread];
		bjmodn26 = _bjmodn26[ithread];			cy26 = _cy26[ithread];
		bjmodn27 = _bjmodn27[ithread];			cy27 = _cy27[ithread];
		bjmodn28 = _bjmodn28[ithread];			cy28 = _cy28[ithread];
		bjmodn29 = _bjmodn29[ithread];			cy29 = _cy29[ithread];
		bjmodn30 = _bjmodn30[ithread];			cy30 = _cy30[ithread];
		bjmodn31 = _bjmodn31[ithread];			cy31 = _cy31[ithread];
		bjmodn32 = _bjmodn32[ithread];			cy32 = _cy32[ithread];
		bjmodn33 = _bjmodn33[ithread];			cy33 = _cy33[ithread];
		bjmodn34 = _bjmodn34[ithread];			cy34 = _cy34[ithread];
		bjmodn35 = _bjmodn35[ithread];			cy35 = _cy35[ithread];
		bjmodn36 = _bjmodn36[ithread];			cy36 = _cy36[ithread];
		bjmodn37 = _bjmodn37[ithread];			cy37 = _cy37[ithread];
		bjmodn38 = _bjmodn38[ithread];			cy38 = _cy38[ithread];
		bjmodn39 = _bjmodn39[ithread];			cy39 = _cy39[ithread];
		bjmodn40 = _bjmodn40[ithread];			cy40 = _cy40[ithread];
		bjmodn41 = _bjmodn41[ithread];			cy41 = _cy41[ithread];
		bjmodn42 = _bjmodn42[ithread];			cy42 = _cy42[ithread];
		bjmodn43 = _bjmodn43[ithread];			cy43 = _cy43[ithread];
		bjmodn44 = _bjmodn44[ithread];			cy44 = _cy44[ithread];
		bjmodn45 = _bjmodn45[ithread];			cy45 = _cy45[ithread];
		bjmodn46 = _bjmodn46[ithread];			cy46 = _cy46[ithread];
		bjmodn47 = _bjmodn47[ithread];			cy47 = _cy47[ithread];
		bjmodn48 = _bjmodn48[ithread];			cy48 = _cy48[ithread];
		bjmodn49 = _bjmodn49[ithread];			cy49 = _cy49[ithread];
		bjmodn50 = _bjmodn50[ithread];			cy50 = _cy50[ithread];
		bjmodn51 = _bjmodn51[ithread];			cy51 = _cy51[ithread];
	#endif

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
		#ifdef USE_SSE2
			for(j = jstart; j < jhi; j += 4)
			{
			/* In SSE2 mode, data are arranged in [re0,re1,im0,im1] quartets, not the usual [re0,im0],[re1,im1] pairs.
			Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
			but to point to the proper real datum, we need to bit-reverse bits <0:1> of j, i.e. [0,1,2,3] ==> [0,2,1,3].
			*/
				j1 = (j & mask01) + br4[j&3];
		#elif defined(USE_SSE2)	/* This allows us to use #if 0 above and disable sse2-based *computation*, while still using sse2-style data layout */
			for(j = jstart; j < jhi; j += 2)	/* Each inner loop execution processes (radix(1)*nwt) array data.	*/
			{
				j1 = (j & mask01) + br4[j&3];
		#else
			for(j = jstart; j < jhi; j += 2)	/* Each inner loop execution processes (radix(1)*nwt) array data.	*/
			{
				j1 =  j;
		#endif
				j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1+RE_IM_STRIDE;

		#ifdef USE_SSE2

			/* Outputs in SSE2 modes are temps 2*13*16 = 26*16 = 0x1a0 bytes apart: */
		  #if !GCC_ASM_FULL_INLINE
			// Reorder blocks to yield sequentially increasing a-array offsets:
/* Block 1 : */	add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00r, 0x1a0)
/* Block 4 : */	add3 = &a[j1+p04];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03r, 0x1a0)
/* Block 7 : */	add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06r, 0x1a0)
/* Block 10: */	add2 = &a[j1+p12];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09r, 0x1a0)
/* Block 13: */	add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0cr, 0x1a0)
/* Block 3 : */	add3 = &a[j1+p20];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02r, 0x1a0)
/* Block 6 : */	add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05r, 0x1a0)
/* Block 9 : */	add2 = &a[j1+p28];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08r, 0x1a0)
/* Block 12: */	add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0br, 0x1a0)
/* Block 2 : */	add3 = &a[j1+p36];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01r, 0x1a0)
/* Block 5 : */	add1 = &a[j1+p40];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04r, 0x1a0)
/* Block 8 : */	add2 = &a[j1+p44];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07r, 0x1a0)
/* Block 11: */	add0 = &a[j1+p48];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0ar, 0x1a0)
		  #else
			add0 = &a[j1    ];
			SSE2_RADIX52_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,r00r);
		  #endif

			/* Radix-13 DFT uses adjacent temps for inputs, outputs are (cyclic) with pXXr having XX += 40 (40*32 bytes = 0x500) or XX -= 12 (-12*32 bytes = -0x180) between successive outputs: */
																													/*  a1p00r,a1p31r,a1p22r,a1p13r,a1p04r,a1p35r,a1p26r,a1p17r,a1p08r,a1p39r,a1p2ar,a1p1br,a1p0cr */
			SSE2_RADIX_13_DFT(rad13_const, r00r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180, s1p00r, 0x500, 0x380, 0x200, 0x080, 0x580, 0x400, 0x280, 0x100, 0x600, 0x480, 0x300, 0x180)
																													/*  a1p30r,a1p21r,a1p12r,a1p03r,a1p34r,a1p25r,a1p16r,a1p07r,a1p38r,a1p29r,a1p1ar,a1p0br,a1p3cr */
			SSE2_RADIX_13_DFT(rad13_const, r10r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180, s1p30r,-0x180,-0x300,-0x480, 0x080,-0x100,-0x280,-0x400, 0x100,-0x080,-0x200,-0x380, 0x180)
																													/*  a1p20r,a1p11r,a1p02r,a1p33r,a1p24r,a1p15r,a1p06r,a1p37r,a1p28r,a1p19r,a1p0ar,a1p3br,a1p2cr */
			SSE2_RADIX_13_DFT(rad13_const, r20r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180, s1p20r,-0x180,-0x300, 0x200, 0x080,-0x100,-0x280, 0x280, 0x100,-0x080,-0x200, 0x300, 0x180)
																													/*  a1p10r,a1p01r,a1p32r,a1p23r,a1p14r,a1p05r,a1p36r,a1p27r,a1p18r,a1p09r,a1p3ar,a1p2br,a1p1cr */
			SSE2_RADIX_13_DFT(rad13_const, r30r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180, s1p10r,-0x180, 0x380, 0x200, 0x080,-0x100, 0x400, 0x280, 0x100,-0x080, 0x480, 0x300, 0x180)

		#else	/* !USE_SSE2 */

		/*...gather the needed data (52 64-bit complex, i.e. 104 64-bit reals) and do 13 radix-4 transforms...*/
						/*                                   inputs                                   */ /*              outputs              */
			jt = j1    ; jp = j2    ;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,rt,it);
			jt = j1+p36; jp = j2+p36;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,rt,it);
			jt = j1+p20; jp = j2+p20;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,rt,it);
			jt = j1+p04; jp = j2+p04;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,rt,it);
			jt = j1+p40; jp = j2+p40;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,rt,it);
			jt = j1+p24; jp = j2+p24;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,rt,it);
			jt = j1+p08; jp = j2+p08;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,rt,it);
			jt = j1+p44; jp = j2+p44;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,rt,it);
			jt = j1+p28; jp = j2+p28;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,rt,it);
			jt = j1+p12; jp = j2+p12;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,rt,it);
			jt = j1+p48; jp = j2+p48;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,rt,it);
			jt = j1+p32; jp = j2+p32;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0br,t0bi,t1br,t1bi,t2br,t2bi,t3br,t3bi,rt,it);
			jt = j1+p16; jp = j2+p16;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0cr,t0ci,t1cr,t1ci,t2cr,t2ci,t3cr,t3ci,rt,it);
			/*...and now do 4 radix-13 transforms. Convert output p-indices from decimal to base-13: */
													/*                                                  inputs                                                                      */ /* outputs: --> */
			RADIX_13_DFT(t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,a1p00r,a1p00i,a1p31r,a1p31i,a1p22r,a1p22i,a1p13r,a1p13i,a1p04r,a1p04i,a1p35r,a1p35i,a1p26r,a1p26i,a1p17r,a1p17i,a1p08r,a1p08i,a1p39r,a1p39i,a1p2ar,a1p2ai,a1p1br,a1p1bi,a1p0cr,a1p0ci);
			RADIX_13_DFT(t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,a1p30r,a1p30i,a1p21r,a1p21i,a1p12r,a1p12i,a1p03r,a1p03i,a1p34r,a1p34i,a1p25r,a1p25i,a1p16r,a1p16i,a1p07r,a1p07i,a1p38r,a1p38i,a1p29r,a1p29i,a1p1ar,a1p1ai,a1p0br,a1p0bi,a1p3cr,a1p3ci);
			RADIX_13_DFT(t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,a1p20r,a1p20i,a1p11r,a1p11i,a1p02r,a1p02i,a1p33r,a1p33i,a1p24r,a1p24i,a1p15r,a1p15i,a1p06r,a1p06i,a1p37r,a1p37i,a1p28r,a1p28i,a1p19r,a1p19i,a1p0ar,a1p0ai,a1p3br,a1p3bi,a1p2cr,a1p2ci);
			RADIX_13_DFT(t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,a1p10r,a1p10i,a1p01r,a1p01i,a1p32r,a1p32i,a1p23r,a1p23i,a1p14r,a1p14i,a1p05r,a1p05i,a1p36r,a1p36i,a1p27r,a1p27i,a1p18r,a1p18i,a1p09r,a1p09i,a1p3ar,a1p3ai,a1p2br,a1p2bi,a1p1cr,a1p1ci);

		#endif

	/*...Now do the carries. Since the outputs would
		normally be getting dispatched to 52 separate blocks of the A-array, we need 52 separate carries.	*/

				l= j & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/
		#ifdef USE_SSE2

				tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				tmp->re = wtl;		tmp->im = wtl;	++tmp;
				tmp->re = wtn;		tmp->im = wtn;	++tmp;
				tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
				tmp->re = wtnm1;	tmp->im = wtnm1;

				add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

		  #if defined(COMPILER_TYPE_MSVC)

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy04,cy06,bjmodn04);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy08,cy10,bjmodn08);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p0cr,add1,add2,add3,cy12,cy14,bjmodn12);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p13r,add1,add2,add3,cy16,cy18,bjmodn16);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p17r,add1,add2,add3,cy20,cy22,bjmodn20);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p1br,add1,add2,add3,cy24,cy26,bjmodn24);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p22r,add1,add2,add3,cy28,cy30,bjmodn28);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p26r,add1,add2,add3,cy32,cy34,bjmodn32);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p2ar,add1,add2,add3,cy36,cy38,bjmodn36);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p31r,add1,add2,add3,cy40,cy42,bjmodn40);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p35r,add1,add2,add3,cy44,cy46,bjmodn44);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p39r,add1,add2,add3,cy48,cy50,bjmodn48);
			  #else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy04,cy06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy08,cy10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p0cr,add1,add2,add3,cy12,cy14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p13r,add1,add2,add3,cy16,cy18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p17r,add1,add2,add3,cy20,cy22,bjmodn20);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p1br,add1,add2,add3,cy24,cy26,bjmodn24);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p22r,add1,add2,add3,cy28,cy30,bjmodn28);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p26r,add1,add2,add3,cy32,cy34,bjmodn32);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p2ar,add1,add2,add3,cy36,cy38,bjmodn36);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p31r,add1,add2,add3,cy40,cy42,bjmodn40);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p35r,add1,add2,add3,cy44,cy46,bjmodn44);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p39r,add1,add2,add3,cy48,cy50,bjmodn48);
			  #endif

		  #else	/* GCC-style inline ASM: */

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p0cr,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p13r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p17r,add1,add2,add3,cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p1br,add1,add2,add3,cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p22r,add1,add2,add3,cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p26r,add1,add2,add3,cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p2ar,add1,add2,add3,cy36,cy38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p31r,add1,add2,add3,cy40,cy42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p35r,add1,add2,add3,cy44,cy46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p39r,add1,add2,add3,cy48,cy50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p0cr,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p13r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p17r,add1,add2,add3,cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p1br,add1,add2,add3,cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p22r,add1,add2,add3,cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p26r,add1,add2,add3,cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p2ar,add1,add2,add3,cy36,cy38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p31r,add1,add2,add3,cy40,cy42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p35r,add1,add2,add3,cy44,cy46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p39r,add1,add2,add3,cy48,cy50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #endif

				/* Bizarre - when I disabled the diagnostic prints above and below, the resulting GCC build immediately gave
					fatal roundoff errors starting on iteration #5 - so insert the bogus [never taken] if() here as a workaround.
					Equally bizarre, inserting the bogus if() *before* the 4 carry-macro calls above gave the correct result as well,
					but ran fully 10% slower. Good old GCC...
				*/
				if(j < 0)
				{
					fprintf(stderr, "Iter %3d\n",iter);
				}

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

				tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				tmp->re = wtl;		tmp->im = wtl;	++tmp;
				tmp->re = wtn;		tmp->im = wtn;	++tmp;
				tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
				tmp->re = wtnm1;	tmp->im = wtnm1;

			/*	i =((uint32)(sw - *bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];

		  #if defined(COMPILER_TYPE_MSVC)

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy00,cy02,bjmodn00);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy04,cy06,bjmodn04);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy08,cy10,bjmodn08);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p0cr,add1,add2,     cy12,cy14,bjmodn12);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p13r,add1,add2,     cy16,cy18,bjmodn16);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p17r,add1,add2,     cy20,cy22,bjmodn20);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p1br,add1,add2,     cy24,cy26,bjmodn24);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p22r,add1,add2,     cy28,cy30,bjmodn28);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p26r,add1,add2,     cy32,cy34,bjmodn32);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p2ar,add1,add2,     cy36,cy38,bjmodn36);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p31r,add1,add2,     cy40,cy42,bjmodn40);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p35r,add1,add2,     cy44,cy46,bjmodn44);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p39r,add1,add2,     cy48,cy50,bjmodn48);
			  #else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy00,cy02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy04,cy06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy08,cy10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p0cr,add1,add2,     cy12,cy14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p13r,add1,add2,     cy16,cy18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p17r,add1,add2,     cy20,cy22,bjmodn20);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p1br,add1,add2,     cy24,cy26,bjmodn24);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p22r,add1,add2,     cy28,cy30,bjmodn28);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p26r,add1,add2,     cy32,cy34,bjmodn32);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p2ar,add1,add2,     cy36,cy38,bjmodn36);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p31r,add1,add2,     cy40,cy42,bjmodn40);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p35r,add1,add2,     cy44,cy46,bjmodn44);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p39r,add1,add2,     cy48,cy50,bjmodn48);
			  #endif

		  #else	/* GCC-style inline ASM: */

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p0cr,add1,add2,     cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p13r,add1,add2,     cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p17r,add1,add2,     cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p1br,add1,add2,     cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p22r,add1,add2,     cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p26r,add1,add2,     cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p2ar,add1,add2,     cy36,cy38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p31r,add1,add2,     cy40,cy42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p35r,add1,add2,     cy44,cy46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p39r,add1,add2,     cy48,cy50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p0cr,add1,add2,     cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p13r,add1,add2,     cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p17r,add1,add2,     cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p1br,add1,add2,     cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p22r,add1,add2,     cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p26r,add1,add2,     cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p2ar,add1,add2,     cy36,cy38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p31r,add1,add2,     cy40,cy42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p35r,add1,add2,     cy44,cy46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p39r,add1,add2,     cy48,cy50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #endif

		  #endif

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

		#else	/* #ifdef USE_SSE2 */

				/*...set0 is slightly different from others:	*/
			   cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy00,bjmodn00   );
				cmplx_carry_norm_errcheck(a1p01r,a1p01i,cy01,bjmodn01,1 );
				cmplx_carry_norm_errcheck(a1p02r,a1p02i,cy02,bjmodn02,2 );
				cmplx_carry_norm_errcheck(a1p03r,a1p03i,cy03,bjmodn03,3 );
				cmplx_carry_norm_errcheck(a1p04r,a1p04i,cy04,bjmodn04,4 );
				cmplx_carry_norm_errcheck(a1p05r,a1p05i,cy05,bjmodn05,5 );
				cmplx_carry_norm_errcheck(a1p06r,a1p06i,cy06,bjmodn06,6 );
				cmplx_carry_norm_errcheck(a1p07r,a1p07i,cy07,bjmodn07,7 );
				cmplx_carry_norm_errcheck(a1p08r,a1p08i,cy08,bjmodn08,8 );
				cmplx_carry_norm_errcheck(a1p09r,a1p09i,cy09,bjmodn09,9 );
				cmplx_carry_norm_errcheck(a1p0ar,a1p0ai,cy10,bjmodn10,10);
				cmplx_carry_norm_errcheck(a1p0br,a1p0bi,cy11,bjmodn11,11);
				cmplx_carry_norm_errcheck(a1p0cr,a1p0ci,cy12,bjmodn12,12);
				cmplx_carry_norm_errcheck(a1p10r,a1p10i,cy13,bjmodn13,13);
				cmplx_carry_norm_errcheck(a1p11r,a1p11i,cy14,bjmodn14,14);
				cmplx_carry_norm_errcheck(a1p12r,a1p12i,cy15,bjmodn15,15);
				cmplx_carry_norm_errcheck(a1p13r,a1p13i,cy16,bjmodn16,16);
				cmplx_carry_norm_errcheck(a1p14r,a1p14i,cy17,bjmodn17,17);
				cmplx_carry_norm_errcheck(a1p15r,a1p15i,cy18,bjmodn18,18);
				cmplx_carry_norm_errcheck(a1p16r,a1p16i,cy19,bjmodn19,19);
				cmplx_carry_norm_errcheck(a1p17r,a1p17i,cy20,bjmodn20,20);
				cmplx_carry_norm_errcheck(a1p18r,a1p18i,cy21,bjmodn21,21);
				cmplx_carry_norm_errcheck(a1p19r,a1p19i,cy22,bjmodn22,22);
				cmplx_carry_norm_errcheck(a1p1ar,a1p1ai,cy23,bjmodn23,23);
				cmplx_carry_norm_errcheck(a1p1br,a1p1bi,cy24,bjmodn24,24);
				cmplx_carry_norm_errcheck(a1p1cr,a1p1ci,cy25,bjmodn25,25);
				cmplx_carry_norm_errcheck(a1p20r,a1p20i,cy26,bjmodn26,26);
				cmplx_carry_norm_errcheck(a1p21r,a1p21i,cy27,bjmodn27,27);
				cmplx_carry_norm_errcheck(a1p22r,a1p22i,cy28,bjmodn28,28);
				cmplx_carry_norm_errcheck(a1p23r,a1p23i,cy29,bjmodn29,29);
				cmplx_carry_norm_errcheck(a1p24r,a1p24i,cy30,bjmodn30,30);
				cmplx_carry_norm_errcheck(a1p25r,a1p25i,cy31,bjmodn31,31);
				cmplx_carry_norm_errcheck(a1p26r,a1p26i,cy32,bjmodn32,32);
				cmplx_carry_norm_errcheck(a1p27r,a1p27i,cy33,bjmodn33,33);
				cmplx_carry_norm_errcheck(a1p28r,a1p28i,cy34,bjmodn34,34);
				cmplx_carry_norm_errcheck(a1p29r,a1p29i,cy35,bjmodn35,35);
				cmplx_carry_norm_errcheck(a1p2ar,a1p2ai,cy36,bjmodn36,36);
				cmplx_carry_norm_errcheck(a1p2br,a1p2bi,cy37,bjmodn37,37);
				cmplx_carry_norm_errcheck(a1p2cr,a1p2ci,cy38,bjmodn38,38);
				cmplx_carry_norm_errcheck(a1p30r,a1p30i,cy39,bjmodn39,39);
				cmplx_carry_norm_errcheck(a1p31r,a1p31i,cy40,bjmodn40,40);
				cmplx_carry_norm_errcheck(a1p32r,a1p32i,cy41,bjmodn41,41);
				cmplx_carry_norm_errcheck(a1p33r,a1p33i,cy42,bjmodn42,42);
				cmplx_carry_norm_errcheck(a1p34r,a1p34i,cy43,bjmodn43,43);
				cmplx_carry_norm_errcheck(a1p35r,a1p35i,cy44,bjmodn44,44);
				cmplx_carry_norm_errcheck(a1p36r,a1p36i,cy45,bjmodn45,45);
				cmplx_carry_norm_errcheck(a1p37r,a1p37i,cy46,bjmodn46,46);
				cmplx_carry_norm_errcheck(a1p38r,a1p38i,cy47,bjmodn47,47);
				cmplx_carry_norm_errcheck(a1p39r,a1p39i,cy48,bjmodn48,48);
				cmplx_carry_norm_errcheck(a1p3ar,a1p3ai,cy49,bjmodn49,49);
				cmplx_carry_norm_errcheck(a1p3br,a1p3bi,cy50,bjmodn50,50);
				cmplx_carry_norm_errcheck(a1p3cr,a1p3ci,cy51,bjmodn51,51);
/*
if(j==jstart && !j && full_pass && (a1p01r > 0.0)) {
	fprintf(stderr, "maxerr = %20.15f!\n",maxerr);
}
*/
				i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
				co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		#endif	/* #ifdef USE_SSE2 */

			/*...The radix-52 DIF pass is here:	*/

		#ifdef USE_SSE2

			/*...gather the needed data (52 64-bit complex, i.e. 104 64-bit reals) and do 4 radix-13 transforms...*/
			/* Radix-13 DFT inputs are (cyclic) with pXXr having XX += 48 (48*32 bytes = +0x600) or XX -= 4 (4*32 bytes = -0x080), outputs are adjacent 32-byte-separated temps: */
			SSE2_RADIX_13_DFT(rad13_const, s1p00r, 0x600, 0x580, 0x500, 0x480, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, r00r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180)
			SSE2_RADIX_13_DFT(rad13_const, s1p30r,-0x080,-0x100,-0x180,-0x200,-0x280,-0x300,-0x380,-0x400,-0x480, 0x180, 0x100, 0x080, r10r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180)
			SSE2_RADIX_13_DFT(rad13_const, s1p20r,-0x080,-0x100,-0x180,-0x200,-0x280,-0x300, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, r20r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180)
			SSE2_RADIX_13_DFT(rad13_const, s1p10r,-0x080,-0x100,-0x180, 0x480, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, r30r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180)

			/*...and now do 13 radix-4 transforms...*/
			/* Inputs in SSE2 modes are temps 2*13*16 = 26*16 = 0x1a0 bytes apart. Notice how the add* indices after the first row repeat with period 4: */
		  #if !GCC_ASM_FULL_INLINE
			// Reorder blocks to yield sequentially increasing a-array offsets:
/* Block 01 : */	add0 = &a[j1    ];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00r, 0x1a0)
/* Block 13 : */	add3 = &a[j1+p04];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0cr, 0x1a0)
/* Block 12 : */	add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0br, 0x1a0)
/* Block 11 : */	add2 = &a[j1+p12];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0ar, 0x1a0)
/* Block 10 : */	add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09r, 0x1a0)
/* Block 09 : */	add3 = &a[j1+p20];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08r, 0x1a0)
/* Block 08 : */	add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07r, 0x1a0)
/* Block 07 : */	add2 = &a[j1+p28];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06r, 0x1a0)
/* Block 06 : */	add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05r, 0x1a0)
/* Block 05 : */	add3 = &a[j1+p36];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04r, 0x1a0)
/* Block 04 : */	add1 = &a[j1+p40];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03r, 0x1a0)
/* Block 03 : */	add2 = &a[j1+p44];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02r, 0x1a0)
/* Block 02 : */	add0 = &a[j1+p48];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01r, 0x1a0)
		  #else
			add0 = &a[j1    ];
			SSE2_RADIX52_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,r00r);
		  #endif

		#else	/* !USE_SSE2 */

		/*...gather the needed data (52 64-bit complex, i.e 104 64-bit reals) and do 4 radix-13 transforms...
		Twiddleless version arranges 4 sets of radix-13 DFT inputs as follows: 0 in upper left corner, decrement 4 horizontally and 13 vertically:
									decimal:											base-13:
			RADIX_13_DFT(00,48,44,40,36,32,28,24,20,16,12,08,04)		00,39,35,31,2a,26,22,1b,17,13,0c,08,04
			RADIX_13_DFT(39,35,31,27,23,19,15,11,07,03,51,47,43)		30,29,25,21,1a,16,12,0b,07,03,3c,38,34
			RADIX_13_DFT(26,22,18,14,10,06,02,50,46,42,38,34,30)		20,19,15,11,0a,06,02,3b,37,33,2c,28,24
			RADIX_13_DFT(13,09,05,01,49,45,41,37,33,29,25,21,17)		10,09,05,01,3a,36,32,2b,27,23,1c,18,14

		Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF] to properly permute the radix-4 DFT outputs.
		*/
			RADIX_13_DFT(a1p00r,a1p00i,a1p39r,a1p39i,a1p35r,a1p35i,a1p31r,a1p31i,a1p2ar,a1p2ai,a1p26r,a1p26i,a1p22r,a1p22i,a1p1br,a1p1bi,a1p17r,a1p17i,a1p13r,a1p13i,a1p0cr,a1p0ci,a1p08r,a1p08i,a1p04r,a1p04i,t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci);
			RADIX_13_DFT(a1p30r,a1p30i,a1p29r,a1p29i,a1p25r,a1p25i,a1p21r,a1p21i,a1p1ar,a1p1ai,a1p16r,a1p16i,a1p12r,a1p12i,a1p0br,a1p0bi,a1p07r,a1p07i,a1p03r,a1p03i,a1p3cr,a1p3ci,a1p38r,a1p38i,a1p34r,a1p34i,t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci);
			RADIX_13_DFT(a1p20r,a1p20i,a1p19r,a1p19i,a1p15r,a1p15i,a1p11r,a1p11i,a1p0ar,a1p0ai,a1p06r,a1p06i,a1p02r,a1p02i,a1p3br,a1p3bi,a1p37r,a1p37i,a1p33r,a1p33i,a1p2cr,a1p2ci,a1p28r,a1p28i,a1p24r,a1p24i,t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci);
			RADIX_13_DFT(a1p10r,a1p10i,a1p09r,a1p09i,a1p05r,a1p05i,a1p01r,a1p01i,a1p3ar,a1p3ai,a1p36r,a1p36i,a1p32r,a1p32i,a1p2br,a1p2bi,a1p27r,a1p27i,a1p23r,a1p23i,a1p1cr,a1p1ci,a1p18r,a1p18i,a1p14r,a1p14i,t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci);
			/*...and now do 13 radix-4 transforms.*/
			jt = j1    ; jp = j2    ;	RADIX_04_DIF(t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
			jt = j1+p48; jp = j2+p48;	RADIX_04_DIF(t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
			jt = j1+p44; jp = j2+p44;	RADIX_04_DIF(t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
			jt = j1+p40; jp = j2+p40;	RADIX_04_DIF(t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
			jt = j1+p36; jp = j2+p36;	RADIX_04_DIF(t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
			jt = j1+p32; jp = j2+p32;	RADIX_04_DIF(t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
			jt = j1+p28; jp = j2+p28;	RADIX_04_DIF(t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
			jt = j1+p24; jp = j2+p24;	RADIX_04_DIF(t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
			jt = j1+p20; jp = j2+p20;	RADIX_04_DIF(t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
			jt = j1+p16; jp = j2+p16;	RADIX_04_DIF(t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
			jt = j1+p12; jp = j2+p12;	RADIX_04_DIF(t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
			jt = j1+p08; jp = j2+p08;	RADIX_04_DIF(t0br,t0bi,t1br,t1bi,t2br,t2bi,t3br,t3bi,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
			jt = j1+p04; jp = j2+p04;	RADIX_04_DIF(t0cr,t0ci,t1cr,t1ci,t2cr,t2ci,t3cr,t3ci,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
		#endif

			}

			jstart += nwt;
			jhi    += nwt;

			col += RADIX;
			co3 -= RADIX;
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
	#ifdef USE_SSE2
		_cy00[ithread] = cy00->re;
		_cy01[ithread] = cy00->im;
		_cy02[ithread] = cy02->re;
		_cy03[ithread] = cy02->im;
		_cy04[ithread] = cy04->re;
		_cy05[ithread] = cy04->im;
		_cy06[ithread] = cy06->re;
		_cy07[ithread] = cy06->im;
		_cy08[ithread] = cy08->re;
		_cy09[ithread] = cy08->im;
		_cy10[ithread] = cy10->re;
		_cy11[ithread] = cy10->im;
		_cy12[ithread] = cy12->re;
		_cy13[ithread] = cy12->im;
		_cy14[ithread] = cy14->re;
		_cy15[ithread] = cy14->im;
		_cy16[ithread] = cy16->re;
		_cy17[ithread] = cy16->im;
		_cy18[ithread] = cy18->re;
		_cy19[ithread] = cy18->im;
		_cy20[ithread] = cy20->re;
		_cy21[ithread] = cy20->im;
		_cy22[ithread] = cy22->re;
		_cy23[ithread] = cy22->im;
		_cy24[ithread] = cy24->re;
		_cy25[ithread] = cy24->im;
		_cy26[ithread] = cy26->re;
		_cy27[ithread] = cy26->im;
		_cy28[ithread] = cy28->re;
		_cy29[ithread] = cy28->im;
		_cy30[ithread] = cy30->re;
		_cy31[ithread] = cy30->im;
		_cy32[ithread] = cy32->re;
		_cy33[ithread] = cy32->im;
		_cy34[ithread] = cy34->re;
		_cy35[ithread] = cy34->im;
		_cy36[ithread] = cy36->re;
		_cy37[ithread] = cy36->im;
		_cy38[ithread] = cy38->re;
		_cy39[ithread] = cy38->im;
		_cy40[ithread] = cy40->re;
		_cy41[ithread] = cy40->im;
		_cy42[ithread] = cy42->re;
		_cy43[ithread] = cy42->im;
		_cy44[ithread] = cy44->re;
		_cy45[ithread] = cy44->im;
		_cy46[ithread] = cy46->re;
		_cy47[ithread] = cy46->im;
		_cy48[ithread] = cy48->re;
		_cy49[ithread] = cy48->im;
		_cy50[ithread] = cy50->re;
		_cy51[ithread] = cy50->im;

		if(max_err->re > max_err->im)
			maxerr = max_err->re;
		else
			maxerr = max_err->im;
	#else
		_cy00[ithread] = cy00;
		_cy01[ithread] = cy01;
		_cy02[ithread] = cy02;
		_cy03[ithread] = cy03;
		_cy04[ithread] = cy04;
		_cy05[ithread] = cy05;
		_cy06[ithread] = cy06;
		_cy07[ithread] = cy07;
		_cy08[ithread] = cy08;
		_cy09[ithread] = cy09;
		_cy10[ithread] = cy10;
		_cy11[ithread] = cy11;
		_cy12[ithread] = cy12;
		_cy13[ithread] = cy13;
		_cy14[ithread] = cy14;
		_cy15[ithread] = cy15;
		_cy16[ithread] = cy16;
		_cy17[ithread] = cy17;
		_cy18[ithread] = cy18;
		_cy19[ithread] = cy19;
		_cy20[ithread] = cy20;
		_cy21[ithread] = cy21;
		_cy22[ithread] = cy22;
		_cy23[ithread] = cy23;
		_cy24[ithread] = cy24;
		_cy25[ithread] = cy25;
		_cy26[ithread] = cy26;
		_cy27[ithread] = cy27;
		_cy28[ithread] = cy28;
		_cy29[ithread] = cy29;
		_cy30[ithread] = cy30;
		_cy31[ithread] = cy31;
		_cy32[ithread] = cy32;
		_cy33[ithread] = cy33;
		_cy34[ithread] = cy34;
		_cy35[ithread] = cy35;
		_cy36[ithread] = cy36;
		_cy37[ithread] = cy37;
		_cy38[ithread] = cy38;
		_cy39[ithread] = cy39;
		_cy40[ithread] = cy40;
		_cy41[ithread] = cy41;
		_cy42[ithread] = cy42;
		_cy43[ithread] = cy43;
		_cy44[ithread] = cy44;
		_cy45[ithread] = cy45;
		_cy46[ithread] = cy46;
		_cy47[ithread] = cy47;
		_cy48[ithread] = cy48;
		_cy49[ithread] = cy49;
		_cy50[ithread] = cy50;
		_cy51[ithread] = cy51;
	#endif

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
		ASSERT(HERE, 0x0 == cy52_process_chunk( (void*)(&tdat[j + pool_work_units]) ), "Main-thread task failure!");
	}

  #endif

	struct timespec ns_time;
	ns_time.tv_sec  = 0.0001;// (time_t)seconds
	ns_time.tv_nsec = 0;	// (long)nanoseconds - At least allegedly, but under OS X it seems to be finer-grained than that

	while(tpool && tpool->free_tasks_queue.num_tasks != pool_work_units) {
		ASSERT(HERE, 0 == nanosleep(&ns_time, 0x0), "nanosleep fail!");
	}
//	printf("radix52_ditN_cy_dif1 end  ; #tasks = %d, #free_tasks = %d\n", tpool->tasks_queue.num_tasks, tpool->free_tasks_queue.num_tasks);

	/* Copy the thread-specific output carry data back to shared memory: */
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_maxerr[ithread] = tdat[ithread].maxerr;
		if(maxerr < _maxerr[ithread]) {
			maxerr = _maxerr[ithread];
		}

		_cy00[ithread] = tdat[ithread].cy00;
		_cy01[ithread] = tdat[ithread].cy01;
		_cy02[ithread] = tdat[ithread].cy02;
		_cy03[ithread] = tdat[ithread].cy03;
		_cy04[ithread] = tdat[ithread].cy04;
		_cy05[ithread] = tdat[ithread].cy05;
		_cy06[ithread] = tdat[ithread].cy06;
		_cy07[ithread] = tdat[ithread].cy07;
		_cy08[ithread] = tdat[ithread].cy08;
		_cy09[ithread] = tdat[ithread].cy09;
		_cy10[ithread] = tdat[ithread].cy10;
		_cy11[ithread] = tdat[ithread].cy11;
		_cy12[ithread] = tdat[ithread].cy12;
		_cy13[ithread] = tdat[ithread].cy13;
		_cy14[ithread] = tdat[ithread].cy14;
		_cy15[ithread] = tdat[ithread].cy15;
		_cy16[ithread] = tdat[ithread].cy16;
		_cy17[ithread] = tdat[ithread].cy17;
		_cy18[ithread] = tdat[ithread].cy18;
		_cy19[ithread] = tdat[ithread].cy19;
		_cy20[ithread] = tdat[ithread].cy20;
		_cy21[ithread] = tdat[ithread].cy21;
		_cy22[ithread] = tdat[ithread].cy22;
		_cy23[ithread] = tdat[ithread].cy23;
		_cy24[ithread] = tdat[ithread].cy24;
		_cy25[ithread] = tdat[ithread].cy25;
		_cy26[ithread] = tdat[ithread].cy26;
		_cy27[ithread] = tdat[ithread].cy27;
		_cy28[ithread] = tdat[ithread].cy28;
		_cy29[ithread] = tdat[ithread].cy29;
		_cy30[ithread] = tdat[ithread].cy30;
		_cy31[ithread] = tdat[ithread].cy31;
		_cy32[ithread] = tdat[ithread].cy32;
		_cy33[ithread] = tdat[ithread].cy33;
		_cy34[ithread] = tdat[ithread].cy34;
		_cy35[ithread] = tdat[ithread].cy35;
		_cy36[ithread] = tdat[ithread].cy36;
		_cy37[ithread] = tdat[ithread].cy37;
		_cy38[ithread] = tdat[ithread].cy38;
		_cy39[ithread] = tdat[ithread].cy39;
		_cy40[ithread] = tdat[ithread].cy40;
		_cy41[ithread] = tdat[ithread].cy41;
		_cy42[ithread] = tdat[ithread].cy42;
		_cy43[ithread] = tdat[ithread].cy43;
		_cy44[ithread] = tdat[ithread].cy44;
		_cy45[ithread] = tdat[ithread].cy45;
		_cy46[ithread] = tdat[ithread].cy46;
		_cy47[ithread] = tdat[ithread].cy47;
		_cy48[ithread] = tdat[ithread].cy48;
		_cy49[ithread] = tdat[ithread].cy49;
		_cy50[ithread] = tdat[ithread].cy50;
		_cy51[ithread] = tdat[ithread].cy51;
	}
#endif

	if(full_pass) {
	//≈	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

	/*   Wraparound carry cleanup loop is here: ***
	!
	!   (1) Invert the radix-52 forward DIF FFT of the first block of 52 complex elements in A and unweight;
	!   (2) Propagate cleanup carries among the real and imaginary parts of the 52 outputs of (1);
	!   (3) Reweight and perform a radix-52 forward DIF FFT on the result of (2);
	!   (4) If any of the exit carries from (2) are nonzero, advance to the next 52 elements and repeat (1-4).
	*/
	t00r = _cy00[CY_THREADS - 1];
	t01r = _cy01[CY_THREADS - 1];
	t02r = _cy02[CY_THREADS - 1];
	t03r = _cy03[CY_THREADS - 1];
	t04r = _cy04[CY_THREADS - 1];
	t05r = _cy05[CY_THREADS - 1];
	t06r = _cy06[CY_THREADS - 1];
	t07r = _cy07[CY_THREADS - 1];
	t08r = _cy08[CY_THREADS - 1];
	t09r = _cy09[CY_THREADS - 1];
	t0ar = _cy10[CY_THREADS - 1];
	t0br = _cy11[CY_THREADS - 1];
	t0cr = _cy12[CY_THREADS - 1];
	t10r = _cy13[CY_THREADS - 1];
	t11r = _cy14[CY_THREADS - 1];
	t12r = _cy15[CY_THREADS - 1];
	t13r = _cy16[CY_THREADS - 1];
	t14r = _cy17[CY_THREADS - 1];
	t15r = _cy18[CY_THREADS - 1];
	t16r = _cy19[CY_THREADS - 1];
	t17r = _cy20[CY_THREADS - 1];
	t18r = _cy21[CY_THREADS - 1];
	t19r = _cy22[CY_THREADS - 1];
	t1ar = _cy23[CY_THREADS - 1];
	t1br = _cy24[CY_THREADS - 1];
	t1cr = _cy25[CY_THREADS - 1];
	t20r = _cy26[CY_THREADS - 1];
	t21r = _cy27[CY_THREADS - 1];
	t22r = _cy28[CY_THREADS - 1];
	t23r = _cy29[CY_THREADS - 1];
	t24r = _cy30[CY_THREADS - 1];
	t25r = _cy31[CY_THREADS - 1];
	t26r = _cy32[CY_THREADS - 1];
	t27r = _cy33[CY_THREADS - 1];
	t28r = _cy34[CY_THREADS - 1];
	t29r = _cy35[CY_THREADS - 1];
	t2ar = _cy36[CY_THREADS - 1];
	t2br = _cy37[CY_THREADS - 1];
	t2cr = _cy38[CY_THREADS - 1];
	t30r = _cy39[CY_THREADS - 1];
	t31r = _cy40[CY_THREADS - 1];
	t32r = _cy41[CY_THREADS - 1];
	t33r = _cy42[CY_THREADS - 1];
	t34r = _cy43[CY_THREADS - 1];
	t35r = _cy44[CY_THREADS - 1];
	t36r = _cy45[CY_THREADS - 1];
	t37r = _cy46[CY_THREADS - 1];
	t38r = _cy47[CY_THREADS - 1];
	t39r = _cy48[CY_THREADS - 1];
	t3ar = _cy49[CY_THREADS - 1];
	t3br = _cy50[CY_THREADS - 1];
	t3cr = _cy51[CY_THREADS - 1];

	for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
	{
		ASSERT(HERE, CY_THREADS > 1,"radix52_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
		_cy00[ithread] = _cy00[ithread-1];
		_cy01[ithread] = _cy01[ithread-1];
		_cy02[ithread] = _cy02[ithread-1];
		_cy03[ithread] = _cy03[ithread-1];
		_cy04[ithread] = _cy04[ithread-1];
		_cy05[ithread] = _cy05[ithread-1];
		_cy06[ithread] = _cy06[ithread-1];
		_cy07[ithread] = _cy07[ithread-1];
		_cy08[ithread] = _cy08[ithread-1];
		_cy09[ithread] = _cy09[ithread-1];
		_cy10[ithread] = _cy10[ithread-1];
		_cy11[ithread] = _cy11[ithread-1];
		_cy12[ithread] = _cy12[ithread-1];
		_cy13[ithread] = _cy13[ithread-1];
		_cy14[ithread] = _cy14[ithread-1];
		_cy15[ithread] = _cy15[ithread-1];
		_cy16[ithread] = _cy16[ithread-1];
		_cy17[ithread] = _cy17[ithread-1];
		_cy18[ithread] = _cy18[ithread-1];
		_cy19[ithread] = _cy19[ithread-1];
		_cy20[ithread] = _cy20[ithread-1];
		_cy21[ithread] = _cy21[ithread-1];
		_cy22[ithread] = _cy22[ithread-1];
		_cy23[ithread] = _cy23[ithread-1];
		_cy24[ithread] = _cy24[ithread-1];
		_cy25[ithread] = _cy25[ithread-1];
		_cy26[ithread] = _cy26[ithread-1];
		_cy27[ithread] = _cy27[ithread-1];
		_cy28[ithread] = _cy28[ithread-1];
		_cy29[ithread] = _cy29[ithread-1];
		_cy30[ithread] = _cy30[ithread-1];
		_cy31[ithread] = _cy31[ithread-1];
		_cy32[ithread] = _cy32[ithread-1];
		_cy33[ithread] = _cy33[ithread-1];
		_cy34[ithread] = _cy34[ithread-1];
		_cy35[ithread] = _cy35[ithread-1];
		_cy36[ithread] = _cy36[ithread-1];
		_cy37[ithread] = _cy37[ithread-1];
		_cy38[ithread] = _cy38[ithread-1];
		_cy39[ithread] = _cy39[ithread-1];
		_cy40[ithread] = _cy40[ithread-1];
		_cy41[ithread] = _cy41[ithread-1];
		_cy42[ithread] = _cy42[ithread-1];
		_cy43[ithread] = _cy43[ithread-1];
		_cy44[ithread] = _cy44[ithread-1];
		_cy45[ithread] = _cy45[ithread-1];
		_cy46[ithread] = _cy46[ithread-1];
		_cy47[ithread] = _cy47[ithread-1];
		_cy48[ithread] = _cy48[ithread-1];
		_cy49[ithread] = _cy49[ithread-1];
		_cy50[ithread] = _cy50[ithread-1];
		_cy51[ithread] = _cy51[ithread-1];
	}

	_cy00[0] =+t3cr;	/* ...The wraparound carry is here: */
	_cy01[0] = t00r;
	_cy02[0] = t01r;
	_cy03[0] = t02r;
	_cy04[0] = t03r;
	_cy05[0] = t04r;
	_cy06[0] = t05r;
	_cy07[0] = t06r;
	_cy08[0] = t07r;
	_cy09[0] = t08r;
	_cy10[0] = t09r;
	_cy11[0] = t0ar;
	_cy12[0] = t0br;
	_cy13[0] = t0cr;
	_cy14[0] = t10r;
	_cy15[0] = t11r;
	_cy16[0] = t12r;
	_cy17[0] = t13r;
	_cy18[0] = t14r;
	_cy19[0] = t15r;
	_cy20[0] = t16r;
	_cy21[0] = t17r;
	_cy22[0] = t18r;
	_cy23[0] = t19r;
	_cy24[0] = t1ar;
	_cy25[0] = t1br;
	_cy26[0] = t1cr;
	_cy27[0] = t20r;
	_cy28[0] = t21r;
	_cy29[0] = t22r;
	_cy30[0] = t23r;
	_cy31[0] = t24r;
	_cy32[0] = t25r;
	_cy33[0] = t26r;
	_cy34[0] = t27r;
	_cy35[0] = t28r;
	_cy36[0] = t29r;
	_cy37[0] = t2ar;
	_cy38[0] = t2br;
	_cy39[0] = t2cr;
	_cy40[0] = t30r;
	_cy41[0] = t31r;
	_cy42[0] = t32r;
	_cy43[0] = t33r;
	_cy44[0] = t34r;
	_cy45[0] = t35r;
	_cy46[0] = t36r;
	_cy47[0] = t37r;
	_cy48[0] = t38r;
	_cy49[0] = t39r;
	_cy50[0] = t3ar;
	_cy51[0] = t3br;

	full_pass = 0;
	scale = 1;
	j_jhi = 7;

    for(ithread = 0; ithread < CY_THREADS; ithread++)
    {
		for(j = ithread*pini; j <= ithread*pini + j_jhi; j++)
		{
			jt = j;
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
		}
	}
}	/* endfor(outer) */

	t00r = 0;
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		t00r += fabs(_cy00[0])+fabs(_cy01[0])+fabs(_cy02[0])+fabs(_cy03[0])+fabs(_cy04[0])+fabs(_cy05[0])+fabs(_cy06[0])+fabs(_cy07[0])+fabs(_cy08[0])+fabs(_cy09[0]);
		t00r += fabs(_cy10[0])+fabs(_cy11[0])+fabs(_cy12[0])+fabs(_cy13[0])+fabs(_cy14[0])+fabs(_cy15[0])+fabs(_cy16[0])+fabs(_cy17[0])+fabs(_cy18[0])+fabs(_cy19[0]);
		t00r += fabs(_cy20[0])+fabs(_cy21[0])+fabs(_cy22[0])+fabs(_cy23[0])+fabs(_cy24[0])+fabs(_cy25[0])+fabs(_cy26[0])+fabs(_cy27[0])+fabs(_cy28[0])+fabs(_cy29[0]);
		t00r += fabs(_cy30[0])+fabs(_cy31[0])+fabs(_cy32[0])+fabs(_cy33[0])+fabs(_cy34[0])+fabs(_cy35[0])+fabs(_cy36[0])+fabs(_cy37[0])+fabs(_cy38[0])+fabs(_cy39[0]);
		t00r += fabs(_cy40[0])+fabs(_cy41[0])+fabs(_cy42[0])+fabs(_cy43[0])+fabs(_cy44[0])+fabs(_cy45[0])+fabs(_cy46[0])+fabs(_cy47[0])+fabs(_cy48[0])+fabs(_cy49[0]);
		t00r += fabs(_cy50[0])+fabs(_cy51[0]);

		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
    }

	if(t00r != 0.0)
	{
		sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix52_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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

void radix52_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-52 complex DIF FFT pass on the data in the length-N real vector A.
*/
	int j,jt,jp,j1,j2;
	static int n52,first_entry=TRUE,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40,p44,p48;
	double rt,it,
		t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,
		t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,
		t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,
		t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci;

	if(!first_entry && (n/52) != n52)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n52=n/52;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n52;
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
	}

/*...The radix-52 pass is here.	*/

    for(j=0; j < n52; j += 2)
    {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

		/*...gather the needed data (52 64-bit complex, i.e 104 64-bit reals) and do 4 radix-13 transforms...*/
	/*
	Twiddleless version arranges 4 sets of radix-13 DFT inputs as follows: 0 in upper left corner, decrement 4 horizontally and 13 vertically:

		RADIX_13_DFT(00,48,44,40,36,32,28,24,20,16,12,08,04)
		RADIX_13_DFT(39,35,31,27,23,19,15,11,07,03,51,47,43)
		RADIX_13_DFT(26,22,18,14,10,06,02,50,46,42,38,34,30)
		RADIX_13_DFT(13,09,05,01,49,45,41,37,33,29,25,21,17)

	Use the supercalafragalistic Ancient Chinese Secret index-munging formula [SACSIMPF] to properly permute the radix-4 DFT outputs.
	*/
		jt = j1    ; jp = j2    ;	RADIX_13_DFT(a[jt    ],a[jp    ],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci);
		jt = j1+p03; jp = j2+p03;	RADIX_13_DFT(a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci);
		jt = j1+p02; jp = j2+p02;	RADIX_13_DFT(a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci);
		jt = j1+p01; jp = j2+p01;	RADIX_13_DFT(a[jt+p12],a[jp+p12],a[jt+p08],a[jp+p08],a[jt+p04],a[jp+p04],a[jt    ],a[jp    ],a[jt+p48],a[jp+p48],a[jt+p44],a[jp+p44],a[jt+p40],a[jp+p40],a[jt+p36],a[jp+p36],a[jt+p32],a[jp+p32],a[jt+p28],a[jp+p28],a[jt+p24],a[jp+p24],a[jt+p20],a[jp+p20],a[jt+p16],a[jp+p16],t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci);

		/*...and now do 13 radix-4 transforms.
		The required output permutation is:

			[0, 1, 3, 2,
			48,49,51,50,
			46,47,44,45,
			41,40,42,43,
			39,38,37,36,
			32,33,35,34,
			30,31,28,29,
			25,24,26,27,
			23,22,21,20,
			16,17,19,18,
			14,15,12,13,
			 9, 8,10,11,
			 7, 6, 5, 4].
		*/
												/*          inputs                    */ /*                                      outputs                              */
		jt = j1    ; jp = j2    ;	RADIX_04_DIF(t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
		jt = j1+p48; jp = j2+p48;	RADIX_04_DIF(t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
		jt = j1+p44; jp = j2+p44;	RADIX_04_DIF(t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
		jt = j1+p40; jp = j2+p40;	RADIX_04_DIF(t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
		jt = j1+p36; jp = j2+p36;	RADIX_04_DIF(t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIF(t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIF(t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIF(t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIF(t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIF(t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIF(t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIF(t0br,t0bi,t1br,t1bi,t2br,t2bi,t3br,t3bi,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIF(t0cr,t0ci,t1cr,t1ci,t2cr,t2ci,t3cr,t3ci,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
		/* Totals: 4*radix13 + 13*radix04 = 4*(164 FADD, 64 FMUL) + 13*(16 FADD, 0 FMUL) = 864 FADD, 256 FMUL	*/
	}
}

/***************/

/* Aug 2011: Modify just the pass-1 DIT routine to actually use the radix-52 SSE2 code, in order to test SSE2 support in test_fft_radix: */
void radix52_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-52 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix52_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
*/
	int i,j,j1,j2;
	static int n52,first_entry=TRUE,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40,p44,p48;
#ifdef USE_SSE2
	double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
	static struct complex *sc_arr = 0x0, *sc_ptr;
	static struct complex *two,*rad13_const,*tmp0,*tmp1,*tmp2,*tmp3     /* rad13_const needs 18*16 bytes allocated */
		,*r00r,*r00i,*r01r,*r01i,*r02r,*r02i,*r03r,*r03i,*r04r,*r04i,*r05r,*r05i,*r06r,*r06i,*r07r,*r07i,*r08r,*r08i,*r09r,*r09i,*r0ar,*r0ai,*r0br,*r0bi,*r0cr,*r0ci
		,*r10r,*r10i,*r11r,*r11i,*r12r,*r12i,*r13r,*r13i,*r14r,*r14i,*r15r,*r15i,*r16r,*r16i,*r17r,*r17i,*r18r,*r18i,*r19r,*r19i,*r1ar,*r1ai,*r1br,*r1bi,*r1cr,*r1ci
		,*r20r,*r20i,*r21r,*r21i,*r22r,*r22i,*r23r,*r23i,*r24r,*r24i,*r25r,*r25i,*r26r,*r26i,*r27r,*r27i,*r28r,*r28i,*r29r,*r29i,*r2ar,*r2ai,*r2br,*r2bi,*r2cr,*r2ci
		,*r30r,*r30i,*r31r,*r31i,*r32r,*r32i,*r33r,*r33i,*r34r,*r34i,*r35r,*r35i,*r36r,*r36i,*r37r,*r37i,*r38r,*r38i,*r39r,*r39i,*r3ar,*r3ai,*r3br,*r3bi,*r3cr,*r3ci
		,*s1p00r,*s1p00i,*s1p01r,*s1p01i,*s1p02r,*s1p02i,*s1p03r,*s1p03i,*s1p04r,*s1p04i,*s1p05r,*s1p05i,*s1p06r,*s1p06i,*s1p07r,*s1p07i,*s1p08r,*s1p08i,*s1p09r,*s1p09i,*s1p0ar,*s1p0ai,*s1p0br,*s1p0bi,*s1p0cr,*s1p0ci
		,*s1p10r,*s1p10i,*s1p11r,*s1p11i,*s1p12r,*s1p12i,*s1p13r,*s1p13i,*s1p14r,*s1p14i,*s1p15r,*s1p15i,*s1p16r,*s1p16i,*s1p17r,*s1p17i,*s1p18r,*s1p18i,*s1p19r,*s1p19i,*s1p1ar,*s1p1ai,*s1p1br,*s1p1bi,*s1p1cr,*s1p1ci
		,*s1p20r,*s1p20i,*s1p21r,*s1p21i,*s1p22r,*s1p22i,*s1p23r,*s1p23i,*s1p24r,*s1p24i,*s1p25r,*s1p25i,*s1p26r,*s1p26i,*s1p27r,*s1p27i,*s1p28r,*s1p28i,*s1p29r,*s1p29i,*s1p2ar,*s1p2ai,*s1p2br,*s1p2bi,*s1p2cr,*s1p2ci
		,*s1p30r,*s1p30i,*s1p31r,*s1p31i,*s1p32r,*s1p32i,*s1p33r,*s1p33i,*s1p34r,*s1p34i,*s1p35r,*s1p35i,*s1p36r,*s1p36i,*s1p37r,*s1p37i,*s1p38r,*s1p38i,*s1p39r,*s1p39i,*s1p3ar,*s1p3ai,*s1p3br,*s1p3bi,*s1p3cr,*s1p3ci;

		sc_arr = ALLOC_COMPLEX(sc_arr, 248);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 104x2 16-byte slots of sc_arr for temporaries, next 21 for the constants needed by the radix-13 DFT,
	next 22 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
		tmp0 = sc_ptr;           tmp1 = tmp0 + 0x1a;      tmp2 = tmp1 + 0x1a;      tmp3 = tmp2 + 0x1a;	/* Temps needed so we can do these inits in || */
		s1p00r = tmp0 + 0x00;    s1p10r = tmp1 + 0x00;    s1p20r = tmp2 + 0x00;    s1p30r = tmp3 + 0x00;
		s1p00i = tmp0 + 0x01;    s1p10i = tmp1 + 0x01;    s1p20i = tmp2 + 0x01;    s1p30i = tmp3 + 0x01;
		s1p01r = tmp0 + 0x02;    s1p11r = tmp1 + 0x02;    s1p21r = tmp2 + 0x02;    s1p31r = tmp3 + 0x02;
		s1p01i = tmp0 + 0x03;    s1p11i = tmp1 + 0x03;    s1p21i = tmp2 + 0x03;    s1p31i = tmp3 + 0x03;
		s1p02r = tmp0 + 0x04;    s1p12r = tmp1 + 0x04;    s1p22r = tmp2 + 0x04;    s1p32r = tmp3 + 0x04;
		s1p02i = tmp0 + 0x05;    s1p12i = tmp1 + 0x05;    s1p22i = tmp2 + 0x05;    s1p32i = tmp3 + 0x05;
		s1p03r = tmp0 + 0x06;    s1p13r = tmp1 + 0x06;    s1p23r = tmp2 + 0x06;    s1p33r = tmp3 + 0x06;
		s1p03i = tmp0 + 0x07;    s1p13i = tmp1 + 0x07;    s1p23i = tmp2 + 0x07;    s1p33i = tmp3 + 0x07;
		s1p04r = tmp0 + 0x08;    s1p14r = tmp1 + 0x08;    s1p24r = tmp2 + 0x08;    s1p34r = tmp3 + 0x08;
		s1p04i = tmp0 + 0x09;    s1p14i = tmp1 + 0x09;    s1p24i = tmp2 + 0x09;    s1p34i = tmp3 + 0x09;
		s1p05r = tmp0 + 0x0a;    s1p15r = tmp1 + 0x0a;    s1p25r = tmp2 + 0x0a;    s1p35r = tmp3 + 0x0a;
		s1p05i = tmp0 + 0x0b;    s1p15i = tmp1 + 0x0b;    s1p25i = tmp2 + 0x0b;    s1p35i = tmp3 + 0x0b;
		s1p06r = tmp0 + 0x0c;    s1p16r = tmp1 + 0x0c;    s1p26r = tmp2 + 0x0c;    s1p36r = tmp3 + 0x0c;
		s1p06i = tmp0 + 0x0d;    s1p16i = tmp1 + 0x0d;    s1p26i = tmp2 + 0x0d;    s1p36i = tmp3 + 0x0d;
		s1p07r = tmp0 + 0x0e;    s1p17r = tmp1 + 0x0e;    s1p27r = tmp2 + 0x0e;    s1p37r = tmp3 + 0x0e;
		s1p07i = tmp0 + 0x0f;    s1p17i = tmp1 + 0x0f;    s1p27i = tmp2 + 0x0f;    s1p37i = tmp3 + 0x0f;
		s1p08r = tmp0 + 0x10;    s1p18r = tmp1 + 0x10;    s1p28r = tmp2 + 0x10;    s1p38r = tmp3 + 0x10;
		s1p08i = tmp0 + 0x11;    s1p18i = tmp1 + 0x11;    s1p28i = tmp2 + 0x11;    s1p38i = tmp3 + 0x11;
		s1p09r = tmp0 + 0x12;    s1p19r = tmp1 + 0x12;    s1p29r = tmp2 + 0x12;    s1p39r = tmp3 + 0x12;
		s1p09i = tmp0 + 0x13;    s1p19i = tmp1 + 0x13;    s1p29i = tmp2 + 0x13;    s1p39i = tmp3 + 0x13;
		s1p0ar = tmp0 + 0x14;    s1p1ar = tmp1 + 0x14;    s1p2ar = tmp2 + 0x14;    s1p3ar = tmp3 + 0x14;
		s1p0ai = tmp0 + 0x15;    s1p1ai = tmp1 + 0x15;    s1p2ai = tmp2 + 0x15;    s1p3ai = tmp3 + 0x15;
		s1p0br = tmp0 + 0x16;    s1p1br = tmp1 + 0x16;    s1p2br = tmp2 + 0x16;    s1p3br = tmp3 + 0x16;
		s1p0bi = tmp0 + 0x17;    s1p1bi = tmp1 + 0x17;    s1p2bi = tmp2 + 0x17;    s1p3bi = tmp3 + 0x17;
		s1p0cr = tmp0 + 0x18;    s1p1cr = tmp1 + 0x18;    s1p2cr = tmp2 + 0x18;    s1p3cr = tmp3 + 0x18;
		s1p0ci = tmp0 + 0x19;    s1p1ci = tmp1 + 0x19;    s1p2ci = tmp2 + 0x19;    s1p3ci = tmp3 + 0x19;

		tmp0 = tmp3 + 0x1a;    tmp1 = tmp0 + 0x1a;    tmp2 = tmp1 + 0x1a;    tmp3 = tmp2 + 0x1a;	/* Use these as handy temps */
		r00r = tmp0 + 0x00;    r10r = tmp1 + 0x00;    r20r = tmp2 + 0x00;    r30r = tmp3 + 0x00;
		r00i = tmp0 + 0x01;    r10i = tmp1 + 0x01;    r20i = tmp2 + 0x01;    r30i = tmp3 + 0x01;
		r01r = tmp0 + 0x02;    r11r = tmp1 + 0x02;    r21r = tmp2 + 0x02;    r31r = tmp3 + 0x02;
		r01i = tmp0 + 0x03;    r11i = tmp1 + 0x03;    r21i = tmp2 + 0x03;    r31i = tmp3 + 0x03;
		r02r = tmp0 + 0x04;    r12r = tmp1 + 0x04;    r22r = tmp2 + 0x04;    r32r = tmp3 + 0x04;
		r02i = tmp0 + 0x05;    r12i = tmp1 + 0x05;    r22i = tmp2 + 0x05;    r32i = tmp3 + 0x05;
		r03r = tmp0 + 0x06;    r13r = tmp1 + 0x06;    r23r = tmp2 + 0x06;    r33r = tmp3 + 0x06;
		r03i = tmp0 + 0x07;    r13i = tmp1 + 0x07;    r23i = tmp2 + 0x07;    r33i = tmp3 + 0x07;
		r04r = tmp0 + 0x08;    r14r = tmp1 + 0x08;    r24r = tmp2 + 0x08;    r34r = tmp3 + 0x08;
		r04i = tmp0 + 0x09;    r14i = tmp1 + 0x09;    r24i = tmp2 + 0x09;    r34i = tmp3 + 0x09;
		r05r = tmp0 + 0x0a;    r15r = tmp1 + 0x0a;    r25r = tmp2 + 0x0a;    r35r = tmp3 + 0x0a;
		r05i = tmp0 + 0x0b;    r15i = tmp1 + 0x0b;    r25i = tmp2 + 0x0b;    r35i = tmp3 + 0x0b;
		r06r = tmp0 + 0x0c;    r16r = tmp1 + 0x0c;    r26r = tmp2 + 0x0c;    r36r = tmp3 + 0x0c;
		r06i = tmp0 + 0x0d;    r16i = tmp1 + 0x0d;    r26i = tmp2 + 0x0d;    r36i = tmp3 + 0x0d;
		r07r = tmp0 + 0x0e;    r17r = tmp1 + 0x0e;    r27r = tmp2 + 0x0e;    r37r = tmp3 + 0x0e;
		r07i = tmp0 + 0x0f;    r17i = tmp1 + 0x0f;    r27i = tmp2 + 0x0f;    r37i = tmp3 + 0x0f;
		r08r = tmp0 + 0x10;    r18r = tmp1 + 0x10;    r28r = tmp2 + 0x10;    r38r = tmp3 + 0x10;
		r08i = tmp0 + 0x11;    r18i = tmp1 + 0x11;    r28i = tmp2 + 0x11;    r38i = tmp3 + 0x11;
		r09r = tmp0 + 0x12;    r19r = tmp1 + 0x12;    r29r = tmp2 + 0x12;    r39r = tmp3 + 0x12;
		r09i = tmp0 + 0x13;    r19i = tmp1 + 0x13;    r29i = tmp2 + 0x13;    r39i = tmp3 + 0x13;
		r0ar = tmp0 + 0x14;    r1ar = tmp1 + 0x14;    r2ar = tmp2 + 0x14;    r3ar = tmp3 + 0x14;
		r0ai = tmp0 + 0x15;    r1ai = tmp1 + 0x15;    r2ai = tmp2 + 0x15;    r3ai = tmp3 + 0x15;
		r0br = tmp0 + 0x16;    r1br = tmp1 + 0x16;    r2br = tmp2 + 0x16;    r3br = tmp3 + 0x16;
		r0bi = tmp0 + 0x17;    r1bi = tmp1 + 0x17;    r2bi = tmp2 + 0x17;    r3bi = tmp3 + 0x17;
		r0cr = tmp0 + 0x18;    r1cr = tmp1 + 0x18;    r2cr = tmp2 + 0x18;    r3cr = tmp3 + 0x18;
		r0ci = tmp0 + 0x19;    r1ci = tmp1 + 0x19;    r2ci = tmp2 + 0x19;    r3ci = tmp3 + 0x19;
		/* __cc pointer offsets: */
		tmp0 = tmp3 + 0x1a;
		rad13_const = tmp0 + 0x01;	/* Leave an extra slot at radix13_const-1 for the constant two = 2.0: */
					tmp0->re =  2.0;		tmp0->im =  2.0;	/*	-0x010 = 2.0 */
		++tmp0;		tmp0->re =  DC1;		tmp0->im =  DC1;	/*	0x000 =  DC1 */
		++tmp0;		tmp0->re =  DC3;		tmp0->im =  DC3;	/*	0x010 =  DC3 */
		++tmp0;		tmp0->re =  DC4;		tmp0->im =  DC4;	/*	0x020 =  DC4 */
		++tmp0;		tmp0->re =  DS1;		tmp0->im =  DS1;	/*	0x030 =  DS1 */
		++tmp0;		tmp0->re =  DS2;		tmp0->im =  DS2;	/*	0x040 =  DS2 */
		++tmp0;		tmp0->re =  DS3;		tmp0->im =  DS3;	/*	0x050 =  DS3 */
		++tmp0;		tmp0->re =  DS4;		tmp0->im =  DS4;	/*	0x060 =  DS4 */
		++tmp0;		tmp0->re =  DS5;		tmp0->im =  DS5;	/*	0x070 =  DS5 */
		++tmp0;		tmp0->re = DC23;		tmp0->im = DC23;	/*	0x080 = DC23 */
		++tmp0;		tmp0->re = DC54;		tmp0->im = DC54;	/*	0x090 = DC54 */
		++tmp0;		tmp0->re = DC65;		tmp0->im = DC65;	/*	0x0a0 = DC65 */
		++tmp0;		tmp0->re = DS63;		tmp0->im = DS63;	/*	0x0b0 = DS63 */
		++tmp0;		tmp0->re = DS74;		tmp0->im = DS74;	/*	0x0c0 = DS74 */
		++tmp0;		tmp0->re = DS85;		tmp0->im = DS85;	/*	0x0d0 = DS85 */
		++tmp0;		tmp0->re = DS93;		tmp0->im = DS93;	/*	0x0e0 = DS93 */
		++tmp0;		tmp0->re = DSa4;		tmp0->im = DSa4;	/*	0x0f0 = DSa4 */
		++tmp0;		tmp0->re = DSb5;		tmp0->im = DSb5;	/*	0x100 = DSb5 */

#else
	int jt,jp;
	double rt,it,
		t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,
		t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,
		t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,
		t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci;
#endif

	if(!first_entry && (n/52) != n52)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n52=n/52;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n52;
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
	}

/*...The radix-52 pass is here.	*/

#ifdef USE_SSE2
	for(j = 0; j < n52; j += 4)
	{
	/* In SSE2 mode, data are arranged in [re0,re1,im0,im1] quartets, not the usual [re0,im0],[re1,im1] pairs.
	Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
	but to point to the proper real datum, we need to bit-reverse bits <0:1> of j, i.e. [0,1,2,3] ==> [0,2,1,3].
	*/
		j1 = (j & mask01) + br4[j&3];
#else
	for(j = 0; j < n52; j += 2)	/* Each inner loop execution processes (radix(1)*nwt) array data.	*/
	{
		j1 =  j;
#endif
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;

	/*
		Twiddleless version requires us to swap inputs as follows:

		bit-reversal index array =
		[ 0 26 13 39 1 27 14 40 2 28 15 41 3 29 16 42 4 30 17 43 5 31 18 44 6 32 19 45 7 33 20 46 8 34 21 47 9 35 22 48 10 36 23 49 11 37 24 50 12 38 25 51]

		DIT input-scramble array =
		[ 0 48 44 40 36 32 28 24 20 16 12 8 4 39 35 31 27 23 19 15 11 7 3 51 47 43 26 22 18 14 10 6 2 50 46 42 38 34 30 13 9 5 1 49 45 41 37 33 29 25 21 17]

		DIT input-scramble + bit-reversal array =
		[ 0 26 39 13 48 22 35 9 44 18 31 5 40 14 27 1 36 10 23 49 32 6 19 45 28 2 15 41 24 50 11 37 20 46 7 33 16 42 3 29 12 3 8 51 25 8 34 47 21 4 30 43 17]

		Combined DIT input-scramble array =
		[0  1  3  2
		39 38 37 36
		23 22 21 20
		 7  6  5  4
		41 40 42 43
		25 24 26 27
		 9  8 10 11
		46 47 44 45
		30 31 28 29
		14 15 12 13
		48 49 51 50
		32 33 35 34
		16 17 19 18]
	*/
	#ifdef USE_SSE2

		/* Outputs in SSE2 modes are temps 2*13*16 = 26*16 = 0x1a0 bytes apart: */
		add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00r, 0x1a0)
		add3 = &a[j1+p04];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03r, 0x1a0)
		add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06r, 0x1a0)
		add2 = &a[j1+p12];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09r, 0x1a0)
		add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0cr, 0x1a0)
		add3 = &a[j1+p20];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02r, 0x1a0)
		add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05r, 0x1a0)
		add2 = &a[j1+p28];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08r, 0x1a0)
		add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0br, 0x1a0)
		add3 = &a[j1+p36];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01r, 0x1a0)
		add1 = &a[j1+p40];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04r, 0x1a0)
		add2 = &a[j1+p44];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07r, 0x1a0)
		add0 = &a[j1+p48];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0ar, 0x1a0)

		/* Radix-13 DFT uses adjacent temps for inputs, outputs are (cyclic) with pXXr having XX += 40 (40*32 bytes = 0x500) or XX -= 12 (-12*32 bytes = -0x180) between successive outputs: */
																												/*  a1p00r,a1p31r,a1p22r,a1p13r,a1p04r,a1p35r,a1p26r,a1p17r,a1p08r,a1p39r,a1p2ar,a1p1br,a1p0cr */
		SSE2_RADIX_13_DFT(rad13_const, r00r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180, s1p00r, 0x500, 0x380, 0x200, 0x080, 0x580, 0x400, 0x280, 0x100, 0x600, 0x480, 0x300, 0x180)
																												/*  a1p30r,a1p21r,a1p12r,a1p03r,a1p34r,a1p25r,a1p16r,a1p07r,a1p38r,a1p29r,a1p1ar,a1p0br,a1p3cr */
		SSE2_RADIX_13_DFT(rad13_const, r10r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180, s1p30r,-0x180,-0x300,-0x480, 0x080,-0x100,-0x280,-0x400, 0x100,-0x080,-0x200,-0x380, 0x180)
																												/*  a1p20r,a1p11r,a1p02r,a1p33r,a1p24r,a1p15r,a1p06r,a1p37r,a1p28r,a1p19r,a1p0ar,a1p3br,a1p2cr */
		SSE2_RADIX_13_DFT(rad13_const, r20r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180, s1p20r,-0x180,-0x300, 0x200, 0x080,-0x100,-0x280, 0x280, 0x100,-0x080,-0x200, 0x300, 0x180)
																												/*  a1p10r,a1p01r,a1p32r,a1p23r,a1p14r,a1p05r,a1p36r,a1p27r,a1p18r,a1p09r,a1p3ar,a1p2br,a1p1cr */
		SSE2_RADIX_13_DFT(rad13_const, r30r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180, s1p10r,-0x180, 0x380, 0x200, 0x080,-0x100, 0x400, 0x280, 0x100,-0x080, 0x480, 0x300, 0x180)

		/* Now copy outputs back into main array: */
		tmp0 = s1p00r;	a[j1        ] = tmp0->re;	a[j1        +1] = tmp0->im;	++tmp0;	a[j1        +RE_IM_STRIDE] = tmp0->re;	a[j1        +RE_IM_STRIDE+1] = tmp0->im;
		tmp1 = tmp0+1;	a[j1    +p01] = tmp1->re;	a[j1    +p01+1] = tmp1->im;	++tmp1;	a[j1    +p01+RE_IM_STRIDE] = tmp1->re;	a[j1    +p01+RE_IM_STRIDE+1] = tmp1->im;
		tmp2 = tmp1+1;	a[j1    +p02] = tmp2->re;	a[j1    +p02+1] = tmp2->im;	++tmp2;	a[j1    +p02+RE_IM_STRIDE] = tmp2->re;	a[j1    +p02+RE_IM_STRIDE+1] = tmp2->im;
		tmp3 = tmp2+1;	a[j1    +p03] = tmp3->re;	a[j1    +p03+1] = tmp3->im;	++tmp3;	a[j1    +p03+RE_IM_STRIDE] = tmp3->re;	a[j1    +p03+RE_IM_STRIDE+1] = tmp3->im;

		tmp0 = tmp3+1;	a[j1+p04    ] = tmp0->re;	a[j1+p04    +1] = tmp0->im;	++tmp0;	a[j1+p04    +RE_IM_STRIDE] = tmp0->re;	a[j1+p04    +RE_IM_STRIDE+1] = tmp0->im;
		tmp1 = tmp0+1;	a[j1+p04+p01] = tmp1->re;	a[j1+p04+p01+1] = tmp1->im;	++tmp1;	a[j1+p04+p01+RE_IM_STRIDE] = tmp1->re;	a[j1+p04+p01+RE_IM_STRIDE+1] = tmp1->im;
		tmp2 = tmp1+1;	a[j1+p04+p02] = tmp2->re;	a[j1+p04+p02+1] = tmp2->im;	++tmp2;	a[j1+p04+p02+RE_IM_STRIDE] = tmp2->re;	a[j1+p04+p02+RE_IM_STRIDE+1] = tmp2->im;
		tmp3 = tmp2+1;	a[j1+p04+p03] = tmp3->re;	a[j1+p04+p03+1] = tmp3->im;	++tmp3;	a[j1+p04+p03+RE_IM_STRIDE] = tmp3->re;	a[j1+p04+p03+RE_IM_STRIDE+1] = tmp3->im;

		tmp0 = tmp3+1;	a[j1+p08    ] = tmp0->re;	a[j1+p08    +1] = tmp0->im;	++tmp0;	a[j1+p08    +RE_IM_STRIDE] = tmp0->re;	a[j1+p08    +RE_IM_STRIDE+1] = tmp0->im;
		tmp1 = tmp0+1;	a[j1+p08+p01] = tmp1->re;	a[j1+p08+p01+1] = tmp1->im;	++tmp1;	a[j1+p08+p01+RE_IM_STRIDE] = tmp1->re;	a[j1+p08+p01+RE_IM_STRIDE+1] = tmp1->im;
		tmp2 = tmp1+1;	a[j1+p08+p02] = tmp2->re;	a[j1+p08+p02+1] = tmp2->im;	++tmp2;	a[j1+p08+p02+RE_IM_STRIDE] = tmp2->re;	a[j1+p08+p02+RE_IM_STRIDE+1] = tmp2->im;
		tmp3 = tmp2+1;	a[j1+p08+p03] = tmp3->re;	a[j1+p08+p03+1] = tmp3->im;	++tmp3;	a[j1+p08+p03+RE_IM_STRIDE] = tmp3->re;	a[j1+p08+p03+RE_IM_STRIDE+1] = tmp3->im;

		tmp0 = tmp3+1;	a[j1+p12    ] = tmp0->re;	a[j1+p12    +1] = tmp0->im;	++tmp0;	a[j1+p12    +RE_IM_STRIDE] = tmp0->re;	a[j1+p12    +RE_IM_STRIDE+1] = tmp0->im;
		tmp1 = tmp0+1;	a[j1+p12+p01] = tmp1->re;	a[j1+p12+p01+1] = tmp1->im;	++tmp1;	a[j1+p12+p01+RE_IM_STRIDE] = tmp1->re;	a[j1+p12+p01+RE_IM_STRIDE+1] = tmp1->im;
		tmp2 = tmp1+1;	a[j1+p12+p02] = tmp2->re;	a[j1+p12+p02+1] = tmp2->im;	++tmp2;	a[j1+p12+p02+RE_IM_STRIDE] = tmp2->re;	a[j1+p12+p02+RE_IM_STRIDE+1] = tmp2->im;
		tmp3 = tmp2+1;	a[j1+p12+p03] = tmp3->re;	a[j1+p12+p03+1] = tmp3->im;	++tmp3;	a[j1+p12+p03+RE_IM_STRIDE] = tmp3->re;	a[j1+p12+p03+RE_IM_STRIDE+1] = tmp3->im;

		tmp0 = tmp3+1;	a[j1+p16    ] = tmp0->re;	a[j1+p16    +1] = tmp0->im;	++tmp0;	a[j1+p16    +RE_IM_STRIDE] = tmp0->re;	a[j1+p16    +RE_IM_STRIDE+1] = tmp0->im;
		tmp1 = tmp0+1;	a[j1+p16+p01] = tmp1->re;	a[j1+p16+p01+1] = tmp1->im;	++tmp1;	a[j1+p16+p01+RE_IM_STRIDE] = tmp1->re;	a[j1+p16+p01+RE_IM_STRIDE+1] = tmp1->im;
		tmp2 = tmp1+1;	a[j1+p16+p02] = tmp2->re;	a[j1+p16+p02+1] = tmp2->im;	++tmp2;	a[j1+p16+p02+RE_IM_STRIDE] = tmp2->re;	a[j1+p16+p02+RE_IM_STRIDE+1] = tmp2->im;
		tmp3 = tmp2+1;	a[j1+p16+p03] = tmp3->re;	a[j1+p16+p03+1] = tmp3->im;	++tmp3;	a[j1+p16+p03+RE_IM_STRIDE] = tmp3->re;	a[j1+p16+p03+RE_IM_STRIDE+1] = tmp3->im;

		tmp0 = tmp3+1;	a[j1+p20    ] = tmp0->re;	a[j1+p20    +1] = tmp0->im;	++tmp0;	a[j1+p20    +RE_IM_STRIDE] = tmp0->re;	a[j1+p20    +RE_IM_STRIDE+1] = tmp0->im;
		tmp1 = tmp0+1;	a[j1+p20+p01] = tmp1->re;	a[j1+p20+p01+1] = tmp1->im;	++tmp1;	a[j1+p20+p01+RE_IM_STRIDE] = tmp1->re;	a[j1+p20+p01+RE_IM_STRIDE+1] = tmp1->im;
		tmp2 = tmp1+1;	a[j1+p20+p02] = tmp2->re;	a[j1+p20+p02+1] = tmp2->im;	++tmp2;	a[j1+p20+p02+RE_IM_STRIDE] = tmp2->re;	a[j1+p20+p02+RE_IM_STRIDE+1] = tmp2->im;
		tmp3 = tmp2+1;	a[j1+p20+p03] = tmp3->re;	a[j1+p20+p03+1] = tmp3->im;	++tmp3;	a[j1+p20+p03+RE_IM_STRIDE] = tmp3->re;	a[j1+p20+p03+RE_IM_STRIDE+1] = tmp3->im;

		tmp0 = tmp3+1;	a[j1+p24    ] = tmp0->re;	a[j1+p24    +1] = tmp0->im;	++tmp0;	a[j1+p24    +RE_IM_STRIDE] = tmp0->re;	a[j1+p24    +RE_IM_STRIDE+1] = tmp0->im;
		tmp1 = tmp0+1;	a[j1+p24+p01] = tmp1->re;	a[j1+p24+p01+1] = tmp1->im;	++tmp1;	a[j1+p24+p01+RE_IM_STRIDE] = tmp1->re;	a[j1+p24+p01+RE_IM_STRIDE+1] = tmp1->im;
		tmp2 = tmp1+1;	a[j1+p24+p02] = tmp2->re;	a[j1+p24+p02+1] = tmp2->im;	++tmp2;	a[j1+p24+p02+RE_IM_STRIDE] = tmp2->re;	a[j1+p24+p02+RE_IM_STRIDE+1] = tmp2->im;
		tmp3 = tmp2+1;	a[j1+p24+p03] = tmp3->re;	a[j1+p24+p03+1] = tmp3->im;	++tmp3;	a[j1+p24+p03+RE_IM_STRIDE] = tmp3->re;	a[j1+p24+p03+RE_IM_STRIDE+1] = tmp3->im;

		tmp0 = tmp3+1;	a[j1+p28    ] = tmp0->re;	a[j1+p28    +1] = tmp0->im;	++tmp0;	a[j1+p28    +RE_IM_STRIDE] = tmp0->re;	a[j1+p28    +RE_IM_STRIDE+1] = tmp0->im;
		tmp1 = tmp0+1;	a[j1+p28+p01] = tmp1->re;	a[j1+p28+p01+1] = tmp1->im;	++tmp1;	a[j1+p28+p01+RE_IM_STRIDE] = tmp1->re;	a[j1+p28+p01+RE_IM_STRIDE+1] = tmp1->im;
		tmp2 = tmp1+1;	a[j1+p28+p02] = tmp2->re;	a[j1+p28+p02+1] = tmp2->im;	++tmp2;	a[j1+p28+p02+RE_IM_STRIDE] = tmp2->re;	a[j1+p28+p02+RE_IM_STRIDE+1] = tmp2->im;
		tmp3 = tmp2+1;	a[j1+p28+p03] = tmp3->re;	a[j1+p28+p03+1] = tmp3->im;	++tmp3;	a[j1+p28+p03+RE_IM_STRIDE] = tmp3->re;	a[j1+p28+p03+RE_IM_STRIDE+1] = tmp3->im;

		tmp0 = tmp3+1;	a[j1+p32    ] = tmp0->re;	a[j1+p32    +1] = tmp0->im;	++tmp0;	a[j1+p32    +RE_IM_STRIDE] = tmp0->re;	a[j1+p32    +RE_IM_STRIDE+1] = tmp0->im;
		tmp1 = tmp0+1;	a[j1+p32+p01] = tmp1->re;	a[j1+p32+p01+1] = tmp1->im;	++tmp1;	a[j1+p32+p01+RE_IM_STRIDE] = tmp1->re;	a[j1+p32+p01+RE_IM_STRIDE+1] = tmp1->im;
		tmp2 = tmp1+1;	a[j1+p32+p02] = tmp2->re;	a[j1+p32+p02+1] = tmp2->im;	++tmp2;	a[j1+p32+p02+RE_IM_STRIDE] = tmp2->re;	a[j1+p32+p02+RE_IM_STRIDE+1] = tmp2->im;
		tmp3 = tmp2+1;	a[j1+p32+p03] = tmp3->re;	a[j1+p32+p03+1] = tmp3->im;	++tmp3;	a[j1+p32+p03+RE_IM_STRIDE] = tmp3->re;	a[j1+p32+p03+RE_IM_STRIDE+1] = tmp3->im;

		tmp0 = tmp3+1;	a[j1+p36    ] = tmp0->re;	a[j1+p36    +1] = tmp0->im;	++tmp0;	a[j1+p36    +RE_IM_STRIDE] = tmp0->re;	a[j1+p36    +RE_IM_STRIDE+1] = tmp0->im;
		tmp1 = tmp0+1;	a[j1+p36+p01] = tmp1->re;	a[j1+p36+p01+1] = tmp1->im;	++tmp1;	a[j1+p36+p01+RE_IM_STRIDE] = tmp1->re;	a[j1+p36+p01+RE_IM_STRIDE+1] = tmp1->im;
		tmp2 = tmp1+1;	a[j1+p36+p02] = tmp2->re;	a[j1+p36+p02+1] = tmp2->im;	++tmp2;	a[j1+p36+p02+RE_IM_STRIDE] = tmp2->re;	a[j1+p36+p02+RE_IM_STRIDE+1] = tmp2->im;
		tmp3 = tmp2+1;	a[j1+p36+p03] = tmp3->re;	a[j1+p36+p03+1] = tmp3->im;	++tmp3;	a[j1+p36+p03+RE_IM_STRIDE] = tmp3->re;	a[j1+p36+p03+RE_IM_STRIDE+1] = tmp3->im;

		tmp0 = tmp3+1;	a[j1+p40    ] = tmp0->re;	a[j1+p40    +1] = tmp0->im;	++tmp0;	a[j1+p40    +RE_IM_STRIDE] = tmp0->re;	a[j1+p40    +RE_IM_STRIDE+1] = tmp0->im;
		tmp1 = tmp0+1;	a[j1+p40+p01] = tmp1->re;	a[j1+p40+p01+1] = tmp1->im;	++tmp1;	a[j1+p40+p01+RE_IM_STRIDE] = tmp1->re;	a[j1+p40+p01+RE_IM_STRIDE+1] = tmp1->im;
		tmp2 = tmp1+1;	a[j1+p40+p02] = tmp2->re;	a[j1+p40+p02+1] = tmp2->im;	++tmp2;	a[j1+p40+p02+RE_IM_STRIDE] = tmp2->re;	a[j1+p40+p02+RE_IM_STRIDE+1] = tmp2->im;
		tmp3 = tmp2+1;	a[j1+p40+p03] = tmp3->re;	a[j1+p40+p03+1] = tmp3->im;	++tmp3;	a[j1+p40+p03+RE_IM_STRIDE] = tmp3->re;	a[j1+p40+p03+RE_IM_STRIDE+1] = tmp3->im;

		tmp0 = tmp3+1;	a[j1+p44    ] = tmp0->re;	a[j1+p44    +1] = tmp0->im;	++tmp0;	a[j1+p44    +RE_IM_STRIDE] = tmp0->re;	a[j1+p44    +RE_IM_STRIDE+1] = tmp0->im;
		tmp1 = tmp0+1;	a[j1+p44+p01] = tmp1->re;	a[j1+p44+p01+1] = tmp1->im;	++tmp1;	a[j1+p44+p01+RE_IM_STRIDE] = tmp1->re;	a[j1+p44+p01+RE_IM_STRIDE+1] = tmp1->im;
		tmp2 = tmp1+1;	a[j1+p44+p02] = tmp2->re;	a[j1+p44+p02+1] = tmp2->im;	++tmp2;	a[j1+p44+p02+RE_IM_STRIDE] = tmp2->re;	a[j1+p44+p02+RE_IM_STRIDE+1] = tmp2->im;
		tmp3 = tmp2+1;	a[j1+p44+p03] = tmp3->re;	a[j1+p44+p03+1] = tmp3->im;	++tmp3;	a[j1+p44+p03+RE_IM_STRIDE] = tmp3->re;	a[j1+p44+p03+RE_IM_STRIDE+1] = tmp3->im;

		tmp0 = tmp3+1;	a[j1+p48    ] = tmp0->re;	a[j1+p48    +1] = tmp0->im;	++tmp0;	a[j1+p48    +RE_IM_STRIDE] = tmp0->re;	a[j1+p48    +RE_IM_STRIDE+1] = tmp0->im;
		tmp1 = tmp0+1;	a[j1+p48+p01] = tmp1->re;	a[j1+p48+p01+1] = tmp1->im;	++tmp1;	a[j1+p48+p01+RE_IM_STRIDE] = tmp1->re;	a[j1+p48+p01+RE_IM_STRIDE+1] = tmp1->im;
		tmp2 = tmp1+1;	a[j1+p48+p02] = tmp2->re;	a[j1+p48+p02+1] = tmp2->im;	++tmp2;	a[j1+p48+p02+RE_IM_STRIDE] = tmp2->re;	a[j1+p48+p02+RE_IM_STRIDE+1] = tmp2->im;
		tmp3 = tmp2+1;	a[j1+p48+p03] = tmp3->re;	a[j1+p48+p03+1] = tmp3->im;	++tmp3;	a[j1+p48+p03+RE_IM_STRIDE] = tmp3->re;	a[j1+p48+p03+RE_IM_STRIDE+1] = tmp3->im;

	#else	/* !USE_SSE2 */

		/*...gather the needed data (52 64-bit complex, i.e. 104 64-bit reals) and do 13 radix-4 transforms...*/
					/*                                   inputs                                   */ /*              outputs              */
		jt = j1    ; jp = j2    ;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t00r,t00i,t10r,t10i,t20r,t20i,t30r,t30i,rt,it);
		jt = j1+p36; jp = j2+p36;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t01r,t01i,t11r,t11i,t21r,t21i,t31r,t31i,rt,it);
		jt = j1+p20; jp = j2+p20;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t02r,t02i,t12r,t12i,t22r,t22i,t32r,t32i,rt,it);
		jt = j1+p04; jp = j2+p04;	RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],t03r,t03i,t13r,t13i,t23r,t23i,t33r,t33i,rt,it);
		jt = j1+p40; jp = j2+p40;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t04r,t04i,t14r,t14i,t24r,t24i,t34r,t34i,rt,it);
		jt = j1+p24; jp = j2+p24;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t05r,t05i,t15r,t15i,t25r,t25i,t35r,t35i,rt,it);
		jt = j1+p08; jp = j2+p08;	RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],t06r,t06i,t16r,t16i,t26r,t26i,t36r,t36i,rt,it);
		jt = j1+p44; jp = j2+p44;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t07r,t07i,t17r,t17i,t27r,t27i,t37r,t37i,rt,it);
		jt = j1+p28; jp = j2+p28;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t08r,t08i,t18r,t18i,t28r,t28i,t38r,t38i,rt,it);
		jt = j1+p12; jp = j2+p12;	RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],t09r,t09i,t19r,t19i,t29r,t29i,t39r,t39i,rt,it);
		jt = j1+p48; jp = j2+p48;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0ar,t0ai,t1ar,t1ai,t2ar,t2ai,t3ar,t3ai,rt,it);
		jt = j1+p32; jp = j2+p32;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0br,t0bi,t1br,t1bi,t2br,t2bi,t3br,t3bi,rt,it);
		jt = j1+p16; jp = j2+p16;	RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],t0cr,t0ci,t1cr,t1ci,t2cr,t2ci,t3cr,t3ci,rt,it);

		/*...and now do 4 radix-13 transforms.
		The required output permutation is:

			[ 0,40,28,16, 4,44,32,20, 8,48,36,24,12,
			 39,27,15, 3,43,31,19, 7,47,35,23,11,51,
			 26,14, 2,42,30,18, 6,46,34,22,10,50,38,
			 13, 1,41,29,17, 5,45,33,21, 9,49,37,25]
		*/
												/*                                                  inputs                                                                      */ /* outputs: --> */
		jt = j1    ; jp = j2    ;	RADIX_13_DFT(t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p28],a[jp+p28],a[jt+p16],a[jp+p16],a[jt+p04],a[jp+p04],a[jt+p44],a[jp+p44],a[jt+p32],a[jp+p32],a[jt+p20],a[jp+p20],a[jt+p08],a[jp+p08],a[jt+p48],a[jp+p48],a[jt+p36],a[jp+p36],a[jt+p24],a[jp+p24],a[jt+p12],a[jp+p12]);
		jt = j1+p03; jp = j2+p03;	RADIX_13_DFT(t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,a[jt+p36],a[jp+p36],a[jt+p24],a[jp+p24],a[jt+p12],a[jp+p12],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p28],a[jp+p28],a[jt+p16],a[jp+p16],a[jt+p04],a[jp+p04],a[jt+p44],a[jp+p44],a[jt+p32],a[jp+p32],a[jt+p20],a[jp+p20],a[jt+p08],a[jp+p08],a[jt+p48],a[jp+p48]);
		jt = j1+p02; jp = j2+p02;	RADIX_13_DFT(t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,a[jt+p24],a[jp+p24],a[jt+p12],a[jp+p12],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p28],a[jp+p28],a[jt+p16],a[jp+p16],a[jt+p04],a[jp+p04],a[jt+p44],a[jp+p44],a[jt+p32],a[jp+p32],a[jt+p20],a[jp+p20],a[jt+p08],a[jp+p08],a[jt+p48],a[jp+p48],a[jt+p36],a[jp+p36]);
		jt = j1+p01; jp = j2+p01;	RADIX_13_DFT(t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci,a[jt+p12],a[jp+p12],a[jt    ],a[jp    ],a[jt+p40],a[jp+p40],a[jt+p28],a[jp+p28],a[jt+p16],a[jp+p16],a[jt+p04],a[jp+p04],a[jt+p44],a[jp+p44],a[jt+p32],a[jp+p32],a[jt+p20],a[jp+p20],a[jt+p08],a[jp+p08],a[jt+p48],a[jp+p48],a[jt+p36],a[jp+p36],a[jt+p24],a[jp+p24]);

	#endif	/* USE_SSE2 */
	}
}

/******************** Multithreaded function body: ***************************/

#ifdef USE_PTHREAD

	#ifndef USE_SSE2
		#error pthreaded carry code requires SSE2-enabled build!
	#endif
	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void* 
	cy52_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		const uint32 RADIX = 52;
		const double crnd = 3.0*0x4000000*0x2000000;
		int j,j1,j2,k;
		int l,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
		uint32 p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40,p44,p48;
		double *add0, *add1, *add2, *add3;
		int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09
			,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19
			,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27,*bjmodn28,*bjmodn29
			,*bjmodn30,*bjmodn31,*bjmodn32,*bjmodn33,*bjmodn34,*bjmodn35,*bjmodn36,*bjmodn37,*bjmodn38,*bjmodn39
			,*bjmodn40,*bjmodn41,*bjmodn42,*bjmodn43,*bjmodn44,*bjmodn45,*bjmodn46,*bjmodn47,*bjmodn48,*bjmodn49
			,*bjmodn50,*bjmodn51;
		struct complex *two,*rad13_const, *max_err, *sse2_rnd, *half_arr, *tmp     /* rad13_const needs 18*16 bytes allocated */
			,*r00r,*r00i,*r01r,*r01i,*r02r,*r02i,*r03r,*r03i,*r04r,*r04i,*r05r,*r05i,*r06r,*r06i,*r07r,*r07i,*r08r,*r08i,*r09r,*r09i,*r0ar,*r0ai,*r0br,*r0bi,*r0cr,*r0ci
			,*r10r,*r10i,*r11r,*r11i,*r12r,*r12i,*r13r,*r13i,*r14r,*r14i,*r15r,*r15i,*r16r,*r16i,*r17r,*r17i,*r18r,*r18i,*r19r,*r19i,*r1ar,*r1ai,*r1br,*r1bi,*r1cr,*r1ci
			,*r20r,*r20i,*r21r,*r21i,*r22r,*r22i,*r23r,*r23i,*r24r,*r24i,*r25r,*r25i,*r26r,*r26i,*r27r,*r27i,*r28r,*r28i,*r29r,*r29i,*r2ar,*r2ai,*r2br,*r2bi,*r2cr,*r2ci
			,*r30r,*r30i,*r31r,*r31i,*r32r,*r32i,*r33r,*r33i,*r34r,*r34i,*r35r,*r35i,*r36r,*r36i,*r37r,*r37i,*r38r,*r38i,*r39r,*r39i,*r3ar,*r3ai,*r3br,*r3bi,*r3cr,*r3ci
			,*s1p00r,*s1p00i,*s1p01r,*s1p01i,*s1p02r,*s1p02i,*s1p03r,*s1p03i,*s1p04r,*s1p04i,*s1p05r,*s1p05i,*s1p06r,*s1p06i,*s1p07r,*s1p07i,*s1p08r,*s1p08i,*s1p09r,*s1p09i,*s1p0ar,*s1p0ai,*s1p0br,*s1p0bi,*s1p0cr,*s1p0ci
			,*s1p10r,*s1p10i,*s1p11r,*s1p11i,*s1p12r,*s1p12i,*s1p13r,*s1p13i,*s1p14r,*s1p14i,*s1p15r,*s1p15i,*s1p16r,*s1p16i,*s1p17r,*s1p17i,*s1p18r,*s1p18i,*s1p19r,*s1p19i,*s1p1ar,*s1p1ai,*s1p1br,*s1p1bi,*s1p1cr,*s1p1ci
			,*s1p20r,*s1p20i,*s1p21r,*s1p21i,*s1p22r,*s1p22i,*s1p23r,*s1p23i,*s1p24r,*s1p24i,*s1p25r,*s1p25i,*s1p26r,*s1p26i,*s1p27r,*s1p27i,*s1p28r,*s1p28i,*s1p29r,*s1p29i,*s1p2ar,*s1p2ai,*s1p2br,*s1p2bi,*s1p2cr,*s1p2ci
			,*s1p30r,*s1p30i,*s1p31r,*s1p31i,*s1p32r,*s1p32i,*s1p33r,*s1p33i,*s1p34r,*s1p34i,*s1p35r,*s1p35i,*s1p36r,*s1p36i,*s1p37r,*s1p37i,*s1p38r,*s1p38i,*s1p39r,*s1p39i,*s1p3ar,*s1p3ai,*s1p3br,*s1p3bi,*s1p3cr,*s1p3ci
			,*cy00,*cy02,*cy04,*cy06,*cy08
			,*cy10,*cy12,*cy14,*cy16,*cy18
			,*cy20,*cy22,*cy24,*cy26,*cy28
			,*cy30,*cy32,*cy34,*cy36,*cy38
			,*cy40,*cy42,*cy44,*cy46,*cy48,*cy50;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

		struct cy_thread_data_t* thread_arg = targ;
	// int data:
#if FFT_DEBUG
	int ithread = thread_arg->tid;	/* unique thread index (use for debug) */
	fprintf(dbg_file,"cy52_process_chunk: thread %d, NDIVR = %d, NWT = %d, &wt0,1 = %llx %llx\n"\
		, ithread, thread_arg->ndivr, thread_arg->nwt, (uint64)thread_arg->wt0, (uint64)thread_arg->wt1);
#endif
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

	// double data:
		double maxerr = thread_arg->maxerr;
		double scale = thread_arg->scale;

	// pointer data:
		double *a = thread_arg->arrdat;
		double *wt0 = thread_arg->wt0;
		double *wt1 = thread_arg->wt1;
		int *si = thread_arg->si;

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

		tmp	= thread_arg->s1p00r;
								    two = tmp + 0x1a;	max_err  = two + 0x1a;		sse2_rnd = max_err + 0x1a;	/* Use these as handy temps */
		s1p00r = tmp + 0x00;    s1p10r = two + 0x00;    s1p20r = max_err + 0x00;    s1p30r = sse2_rnd + 0x00;
		s1p00i = tmp + 0x01;    s1p10i = two + 0x01;    s1p20i = max_err + 0x01;    s1p30i = sse2_rnd + 0x01;
		s1p01r = tmp + 0x02;    s1p11r = two + 0x02;    s1p21r = max_err + 0x02;    s1p31r = sse2_rnd + 0x02;
		s1p01i = tmp + 0x03;    s1p11i = two + 0x03;    s1p21i = max_err + 0x03;    s1p31i = sse2_rnd + 0x03;
		s1p02r = tmp + 0x04;    s1p12r = two + 0x04;    s1p22r = max_err + 0x04;    s1p32r = sse2_rnd + 0x04;
		s1p02i = tmp + 0x05;    s1p12i = two + 0x05;    s1p22i = max_err + 0x05;    s1p32i = sse2_rnd + 0x05;
		s1p03r = tmp + 0x06;    s1p13r = two + 0x06;    s1p23r = max_err + 0x06;    s1p33r = sse2_rnd + 0x06;
		s1p03i = tmp + 0x07;    s1p13i = two + 0x07;    s1p23i = max_err + 0x07;    s1p33i = sse2_rnd + 0x07;
		s1p04r = tmp + 0x08;    s1p14r = two + 0x08;    s1p24r = max_err + 0x08;    s1p34r = sse2_rnd + 0x08;
		s1p04i = tmp + 0x09;    s1p14i = two + 0x09;    s1p24i = max_err + 0x09;    s1p34i = sse2_rnd + 0x09;
		s1p05r = tmp + 0x0a;    s1p15r = two + 0x0a;    s1p25r = max_err + 0x0a;    s1p35r = sse2_rnd + 0x0a;
		s1p05i = tmp + 0x0b;    s1p15i = two + 0x0b;    s1p25i = max_err + 0x0b;    s1p35i = sse2_rnd + 0x0b;
		s1p06r = tmp + 0x0c;    s1p16r = two + 0x0c;    s1p26r = max_err + 0x0c;    s1p36r = sse2_rnd + 0x0c;
		s1p06i = tmp + 0x0d;    s1p16i = two + 0x0d;    s1p26i = max_err + 0x0d;    s1p36i = sse2_rnd + 0x0d;
		s1p07r = tmp + 0x0e;    s1p17r = two + 0x0e;    s1p27r = max_err + 0x0e;    s1p37r = sse2_rnd + 0x0e;
		s1p07i = tmp + 0x0f;    s1p17i = two + 0x0f;    s1p27i = max_err + 0x0f;    s1p37i = sse2_rnd + 0x0f;
		s1p08r = tmp + 0x10;    s1p18r = two + 0x10;    s1p28r = max_err + 0x10;    s1p38r = sse2_rnd + 0x10;
		s1p08i = tmp + 0x11;    s1p18i = two + 0x11;    s1p28i = max_err + 0x11;    s1p38i = sse2_rnd + 0x11;
		s1p09r = tmp + 0x12;    s1p19r = two + 0x12;    s1p29r = max_err + 0x12;    s1p39r = sse2_rnd + 0x12;
		s1p09i = tmp + 0x13;    s1p19i = two + 0x13;    s1p29i = max_err + 0x13;    s1p39i = sse2_rnd + 0x13;
		s1p0ar = tmp + 0x14;    s1p1ar = two + 0x14;    s1p2ar = max_err + 0x14;    s1p3ar = sse2_rnd + 0x14;
		s1p0ai = tmp + 0x15;    s1p1ai = two + 0x15;    s1p2ai = max_err + 0x15;    s1p3ai = sse2_rnd + 0x15;
		s1p0br = tmp + 0x16;    s1p1br = two + 0x16;    s1p2br = max_err + 0x16;    s1p3br = sse2_rnd + 0x16;
		s1p0bi = tmp + 0x17;    s1p1bi = two + 0x17;    s1p2bi = max_err + 0x17;    s1p3bi = sse2_rnd + 0x17;
		s1p0cr = tmp + 0x18;    s1p1cr = two + 0x18;    s1p2cr = max_err + 0x18;    s1p3cr = sse2_rnd + 0x18;
		s1p0ci = tmp + 0x19;    s1p1ci = two + 0x19;    s1p2ci = max_err + 0x19;    s1p3ci = sse2_rnd + 0x19;

		tmp = sse2_rnd + 0x1a; two = tmp + 0x1a;    max_err = two + 0x1a;   sse2_rnd = max_err + 0x1a;	/* Use these as handy temps */
		r00r = tmp + 0x00;    r10r = two + 0x00;    r20r = max_err + 0x00;    r30r = sse2_rnd + 0x00;
		r00i = tmp + 0x01;    r10i = two + 0x01;    r20i = max_err + 0x01;    r30i = sse2_rnd + 0x01;
		r01r = tmp + 0x02;    r11r = two + 0x02;    r21r = max_err + 0x02;    r31r = sse2_rnd + 0x02;
		r01i = tmp + 0x03;    r11i = two + 0x03;    r21i = max_err + 0x03;    r31i = sse2_rnd + 0x03;
		r02r = tmp + 0x04;    r12r = two + 0x04;    r22r = max_err + 0x04;    r32r = sse2_rnd + 0x04;
		r02i = tmp + 0x05;    r12i = two + 0x05;    r22i = max_err + 0x05;    r32i = sse2_rnd + 0x05;
		r03r = tmp + 0x06;    r13r = two + 0x06;    r23r = max_err + 0x06;    r33r = sse2_rnd + 0x06;
		r03i = tmp + 0x07;    r13i = two + 0x07;    r23i = max_err + 0x07;    r33i = sse2_rnd + 0x07;
		r04r = tmp + 0x08;    r14r = two + 0x08;    r24r = max_err + 0x08;    r34r = sse2_rnd + 0x08;
		r04i = tmp + 0x09;    r14i = two + 0x09;    r24i = max_err + 0x09;    r34i = sse2_rnd + 0x09;
		r05r = tmp + 0x0a;    r15r = two + 0x0a;    r25r = max_err + 0x0a;    r35r = sse2_rnd + 0x0a;
		r05i = tmp + 0x0b;    r15i = two + 0x0b;    r25i = max_err + 0x0b;    r35i = sse2_rnd + 0x0b;
		r06r = tmp + 0x0c;    r16r = two + 0x0c;    r26r = max_err + 0x0c;    r36r = sse2_rnd + 0x0c;
		r06i = tmp + 0x0d;    r16i = two + 0x0d;    r26i = max_err + 0x0d;    r36i = sse2_rnd + 0x0d;
		r07r = tmp + 0x0e;    r17r = two + 0x0e;    r27r = max_err + 0x0e;    r37r = sse2_rnd + 0x0e;
		r07i = tmp + 0x0f;    r17i = two + 0x0f;    r27i = max_err + 0x0f;    r37i = sse2_rnd + 0x0f;
		r08r = tmp + 0x10;    r18r = two + 0x10;    r28r = max_err + 0x10;    r38r = sse2_rnd + 0x10;
		r08i = tmp + 0x11;    r18i = two + 0x11;    r28i = max_err + 0x11;    r38i = sse2_rnd + 0x11;
		r09r = tmp + 0x12;    r19r = two + 0x12;    r29r = max_err + 0x12;    r39r = sse2_rnd + 0x12;
		r09i = tmp + 0x13;    r19i = two + 0x13;    r29i = max_err + 0x13;    r39i = sse2_rnd + 0x13;
		r0ar = tmp + 0x14;    r1ar = two + 0x14;    r2ar = max_err + 0x14;    r3ar = sse2_rnd + 0x14;
		r0ai = tmp + 0x15;    r1ai = two + 0x15;    r2ai = max_err + 0x15;    r3ai = sse2_rnd + 0x15;
		r0br = tmp + 0x16;    r1br = two + 0x16;    r2br = max_err + 0x16;    r3br = sse2_rnd + 0x16;
		r0bi = tmp + 0x17;    r1bi = two + 0x17;    r2bi = max_err + 0x17;    r3bi = sse2_rnd + 0x17;
		r0cr = tmp + 0x18;    r1cr = two + 0x18;    r2cr = max_err + 0x18;    r3cr = sse2_rnd + 0x18;
		r0ci = tmp + 0x19;    r1ci = two + 0x19;    r2ci = max_err + 0x19;    r3ci = sse2_rnd + 0x19;
		tmp = sse2_rnd + 0x1a;

		rad13_const = tmp + 0x01;	/* Leave an extra slot at radix13_const-1 for the constant two = 2.0: */
		tmp += 0x14;	/* Need 20 16-byte slots for two+sincos, but offset the carry slots by the next-larger multiple of 4 */

		cy00 = tmp + 0x00;
		cy02 = tmp + 0x01;
		cy04 = tmp + 0x02;
		cy06 = tmp + 0x03;
		cy08 = tmp + 0x04;
		cy10 = tmp + 0x05;
		cy12 = tmp + 0x06;
		cy14 = tmp + 0x07;
		cy16 = tmp + 0x08;
		cy18 = tmp + 0x09;
		cy20 = tmp + 0x0a;
		cy22 = tmp + 0x0b;
		cy24 = tmp + 0x0c;
		cy26 = tmp + 0x0d;
		cy28 = tmp + 0x0e;
		cy30 = tmp + 0x0f;
		cy32 = tmp + 0x10;
		cy34 = tmp + 0x11;
		cy36 = tmp + 0x12;
		cy38 = tmp + 0x13;
		cy40 = tmp + 0x14;
		cy42 = tmp + 0x15;
		cy44 = tmp + 0x16;
		cy46 = tmp + 0x17;
		cy48 = tmp + 0x18;
		cy50 = tmp + 0x19;
		tmp += 0x1a;

		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;	/* This table needs 20x16 bytes for Mersenne-mod, and radixx16 for Fermat-mod */

		ASSERT(HERE, (s1p00r == thread_arg->s1p00r), "thread-local memcheck failed!");
		ASSERT(HERE, (half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
		ASSERT(HERE, (sse2_rnd->re == crnd && sse2_rnd->im == crnd), "thread-local memcheck failed!");
		ASSERT(HERE, (half_arr+10)->re * (half_arr+14)->re == 1.0 && (half_arr+10)->im * (half_arr+14)->im == 1.0, "thread-local memcheck failed!");

		max_err->re = 0.0;	max_err->im = 0.0;

		sign_mask = (uint64*)(s1p00r + radix52_creals_in_local_store);
		sse_bw  = sign_mask + 2;
		sse_sw  = sign_mask + 4;
		sse_n   = sign_mask + 6;
		bjmodn00 = (int*)(sign_mask + 8);
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

		/* Init DWT-indices: */	/* init carries	*/
		*bjmodn00 = thread_arg->bjmodn00;	cy00->re = thread_arg->cy00;
		*bjmodn01 = thread_arg->bjmodn01;	cy00->im = thread_arg->cy01;
		*bjmodn02 = thread_arg->bjmodn02;	cy02->re = thread_arg->cy02;
		*bjmodn03 = thread_arg->bjmodn03;	cy02->im = thread_arg->cy03;
		*bjmodn04 = thread_arg->bjmodn04;	cy04->re = thread_arg->cy04;
		*bjmodn05 = thread_arg->bjmodn05;	cy04->im = thread_arg->cy05;
		*bjmodn06 = thread_arg->bjmodn06;	cy06->re = thread_arg->cy06;
		*bjmodn07 = thread_arg->bjmodn07;	cy06->im = thread_arg->cy07;
		*bjmodn08 = thread_arg->bjmodn08;	cy08->re = thread_arg->cy08;
		*bjmodn09 = thread_arg->bjmodn09;	cy08->im = thread_arg->cy09;
		*bjmodn10 = thread_arg->bjmodn10;	cy10->re = thread_arg->cy10;
		*bjmodn11 = thread_arg->bjmodn11;	cy10->im = thread_arg->cy11;
		*bjmodn12 = thread_arg->bjmodn12;	cy12->re = thread_arg->cy12;
		*bjmodn13 = thread_arg->bjmodn13;	cy12->im = thread_arg->cy13;
		*bjmodn14 = thread_arg->bjmodn14;	cy14->re = thread_arg->cy14;
		*bjmodn15 = thread_arg->bjmodn15;	cy14->im = thread_arg->cy15;
		*bjmodn16 = thread_arg->bjmodn16;	cy16->re = thread_arg->cy16;
		*bjmodn17 = thread_arg->bjmodn17;	cy16->im = thread_arg->cy17;
		*bjmodn18 = thread_arg->bjmodn18;	cy18->re = thread_arg->cy18;
		*bjmodn19 = thread_arg->bjmodn19;	cy18->im = thread_arg->cy19;
		*bjmodn20 = thread_arg->bjmodn20;	cy20->re = thread_arg->cy20;
		*bjmodn21 = thread_arg->bjmodn21;	cy20->im = thread_arg->cy21;
		*bjmodn22 = thread_arg->bjmodn22;	cy22->re = thread_arg->cy22;
		*bjmodn23 = thread_arg->bjmodn23;	cy22->im = thread_arg->cy23;
		*bjmodn24 = thread_arg->bjmodn24;	cy24->re = thread_arg->cy24;
		*bjmodn25 = thread_arg->bjmodn25;	cy24->im = thread_arg->cy25;
		*bjmodn26 = thread_arg->bjmodn26;	cy26->re = thread_arg->cy26;
		*bjmodn27 = thread_arg->bjmodn27;	cy26->im = thread_arg->cy27;
		*bjmodn28 = thread_arg->bjmodn28;	cy28->re = thread_arg->cy28;
		*bjmodn29 = thread_arg->bjmodn29;	cy28->im = thread_arg->cy29;
		*bjmodn30 = thread_arg->bjmodn30;	cy30->re = thread_arg->cy30;
		*bjmodn31 = thread_arg->bjmodn31;	cy30->im = thread_arg->cy31;
		*bjmodn32 = thread_arg->bjmodn32;	cy32->re = thread_arg->cy32;
		*bjmodn33 = thread_arg->bjmodn33;	cy32->im = thread_arg->cy33;
		*bjmodn34 = thread_arg->bjmodn34;	cy34->re = thread_arg->cy34;
		*bjmodn35 = thread_arg->bjmodn35;	cy34->im = thread_arg->cy35;
		*bjmodn36 = thread_arg->bjmodn36;	cy36->re = thread_arg->cy36;
		*bjmodn37 = thread_arg->bjmodn37;	cy36->im = thread_arg->cy37;
		*bjmodn38 = thread_arg->bjmodn38;	cy38->re = thread_arg->cy38;
		*bjmodn39 = thread_arg->bjmodn39;	cy38->im = thread_arg->cy39;
		*bjmodn40 = thread_arg->bjmodn40;	cy40->re = thread_arg->cy40;
		*bjmodn41 = thread_arg->bjmodn41;	cy40->im = thread_arg->cy41;
		*bjmodn42 = thread_arg->bjmodn42;	cy42->re = thread_arg->cy42;
		*bjmodn43 = thread_arg->bjmodn43;	cy42->im = thread_arg->cy43;
		*bjmodn44 = thread_arg->bjmodn44;	cy44->re = thread_arg->cy44;
		*bjmodn45 = thread_arg->bjmodn45;	cy44->im = thread_arg->cy45;
		*bjmodn46 = thread_arg->bjmodn46;	cy46->re = thread_arg->cy46;
		*bjmodn47 = thread_arg->bjmodn47;	cy46->im = thread_arg->cy47;
		*bjmodn48 = thread_arg->bjmodn48;	cy48->re = thread_arg->cy48;
		*bjmodn49 = thread_arg->bjmodn49;	cy48->im = thread_arg->cy49;
		*bjmodn50 = thread_arg->bjmodn50;	cy50->re = thread_arg->cy50;
		*bjmodn51 = thread_arg->bjmodn51;	cy50->im = thread_arg->cy51;

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
			for(j = jstart; j < jhi; j += 4)
			{
			/* In SSE2 mode, data are arranged in [re0,re1,im0,im1] quartets, not the usual [re0,im0],[re1,im1] pairs.
			Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
			but to point to the proper real datum, we need to bit-reverse bits <0:1> of j, i.e. [0,1,2,3] ==> [0,2,1,3].
			*/
				j1 = (j & mask01) + br4[j&3];
				j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1+RE_IM_STRIDE;

				/* Outputs in SSE2 modes are temps 2*13*16 = 26*16 = 0x1a0 bytes apart: */
			  #if !GCC_ASM_FULL_INLINE
				// Reorder blocks to yield sequentially increasing a-array offsets:
	/* Block 1 : */	add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00r, 0x1a0)
	/* Block 4 : */	add3 = &a[j1+p04];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03r, 0x1a0)
	/* Block 7 : */	add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06r, 0x1a0)
	/* Block 10: */	add2 = &a[j1+p12];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09r, 0x1a0)
	/* Block 13: */	add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0cr, 0x1a0)
	/* Block 3 : */	add3 = &a[j1+p20];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02r, 0x1a0)
	/* Block 6 : */	add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05r, 0x1a0)
	/* Block 9 : */	add2 = &a[j1+p28];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08r, 0x1a0)
	/* Block 12: */	add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0br, 0x1a0)
	/* Block 2 : */	add3 = &a[j1+p36];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01r, 0x1a0)
	/* Block 5 : */	add1 = &a[j1+p40];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04r, 0x1a0)
	/* Block 8 : */	add2 = &a[j1+p44];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07r, 0x1a0)
	/* Block 11: */	add0 = &a[j1+p48];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0ar, 0x1a0)
			  #else
				add0 = &a[j1    ];
				SSE2_RADIX52_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,r00r);
			  #endif

				/* Radix-13 DFT uses adjacent temps for inputs, outputs are (cyclic) with pXXr having XX += 40 (40*32 bytes = 0x500) or XX -= 12 (-12*32 bytes = -0x180) between successive outputs: */
																														/*  a1p00r,a1p31r,a1p22r,a1p13r,a1p04r,a1p35r,a1p26r,a1p17r,a1p08r,a1p39r,a1p2ar,a1p1br,a1p0cr */
				SSE2_RADIX_13_DFT(rad13_const, r00r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180, s1p00r, 0x500, 0x380, 0x200, 0x080, 0x580, 0x400, 0x280, 0x100, 0x600, 0x480, 0x300, 0x180)
																														/*  a1p30r,a1p21r,a1p12r,a1p03r,a1p34r,a1p25r,a1p16r,a1p07r,a1p38r,a1p29r,a1p1ar,a1p0br,a1p3cr */
				SSE2_RADIX_13_DFT(rad13_const, r10r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180, s1p30r,-0x180,-0x300,-0x480, 0x080,-0x100,-0x280,-0x400, 0x100,-0x080,-0x200,-0x380, 0x180)
																														/*  a1p20r,a1p11r,a1p02r,a1p33r,a1p24r,a1p15r,a1p06r,a1p37r,a1p28r,a1p19r,a1p0ar,a1p3br,a1p2cr */
				SSE2_RADIX_13_DFT(rad13_const, r20r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180, s1p20r,-0x180,-0x300, 0x200, 0x080,-0x100,-0x280, 0x280, 0x100,-0x080,-0x200, 0x300, 0x180)
																														/*  a1p10r,a1p01r,a1p32r,a1p23r,a1p14r,a1p05r,a1p36r,a1p27r,a1p18r,a1p09r,a1p3ar,a1p2br,a1p1cr */
				SSE2_RADIX_13_DFT(rad13_const, r30r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180, s1p10r,-0x180, 0x380, 0x200, 0x080,-0x100, 0x400, 0x280, 0x100,-0x080, 0x480, 0x300, 0x180)

				l= j & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				tmp->re = wtl;		tmp->im = wtl;	++tmp;
				tmp->re = wtn;		tmp->im = wtn;	++tmp;
				tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
				tmp->re = wtnm1;	tmp->im = wtnm1;

				add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
				add2 = &wt1[co2-1];
				add3 = &wt1[co3-1];

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p0cr,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p13r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p17r,add1,add2,add3,cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p1br,add1,add2,add3,cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p22r,add1,add2,add3,cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p26r,add1,add2,add3,cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p2ar,add1,add2,add3,cy36,cy38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p31r,add1,add2,add3,cy40,cy42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p35r,add1,add2,add3,cy44,cy46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck1_2B(s1p39r,add1,add2,add3,cy48,cy50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy00,cy02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p0cr,add1,add2,add3,cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p13r,add1,add2,add3,cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p17r,add1,add2,add3,cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p1br,add1,add2,add3,cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p22r,add1,add2,add3,cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p26r,add1,add2,add3,cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p2ar,add1,add2,add3,cy36,cy38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p31r,add1,add2,add3,cy40,cy42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p35r,add1,add2,add3,cy44,cy46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p39r,add1,add2,add3,cy48,cy50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #endif

				l= (j+2) & (nwt-1);
				n_minus_sil   = n-si[l  ];
				n_minus_silp1 = n-si[l+1];
				sinwt   = si[nwt-l  ];
				sinwtm1 = si[nwt-l-1];

				wtl     =wt0[    l  ];
				wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
				wtlp1   =wt0[    l+1];
				wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

				tmp = half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
				tmp->re = wtl;		tmp->im = wtl;	++tmp;
				tmp->re = wtn;		tmp->im = wtn;	++tmp;
				tmp->re = wtlp1;	tmp->im = wtlp1;++tmp;
				tmp->re = wtnm1;	tmp->im = wtnm1;

			/*	i =((uint32)(sw - *bjmodn0) >> 31);	Don't need this here, since no special index-0 macro in the set below */

				co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
							(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

				add1 = &wt1[col  ];
				add2 = &wt1[co2-1];

			  #ifdef ERR_CHECK_ALL
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p0cr,add1,add2,     cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p13r,add1,add2,     cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p17r,add1,add2,     cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p1br,add1,add2,     cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p22r,add1,add2,     cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p26r,add1,add2,     cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p2ar,add1,add2,     cy36,cy38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p31r,add1,add2,     cy40,cy42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p35r,add1,add2,     cy44,cy46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_errcheck2_2B(s1p39r,add1,add2,     cy48,cy50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #else
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p00r,add1,add2,     cy00,cy02,bjmodn00,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p04r,add1,add2,     cy04,cy06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p08r,add1,add2,     cy08,cy10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p0cr,add1,add2,     cy12,cy14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p13r,add1,add2,     cy16,cy18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p17r,add1,add2,     cy20,cy22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p1br,add1,add2,     cy24,cy26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p22r,add1,add2,     cy28,cy30,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p26r,add1,add2,     cy32,cy34,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p2ar,add1,add2,     cy36,cy38,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p31r,add1,add2,     cy40,cy42,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p35r,add1,add2,     cy44,cy46,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B (s1p39r,add1,add2,     cy48,cy50,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			  #endif

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

				/*...gather the needed data (52 64-bit complex, i.e. 104 64-bit reals) and do 4 radix-13 transforms...*/
				/* Radix-13 DFT inputs are (cyclic) with pXXr having XX += 48 (48*32 bytes = +0x600) or XX -= 4 (4*32 bytes = -0x080), outputs are adjacent 32-byte-separated temps: */
				SSE2_RADIX_13_DFT(rad13_const, s1p00r, 0x600, 0x580, 0x500, 0x480, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, r00r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180)
				SSE2_RADIX_13_DFT(rad13_const, s1p30r,-0x080,-0x100,-0x180,-0x200,-0x280,-0x300,-0x380,-0x400,-0x480, 0x180, 0x100, 0x080, r10r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180)
				SSE2_RADIX_13_DFT(rad13_const, s1p20r,-0x080,-0x100,-0x180,-0x200,-0x280,-0x300, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, r20r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180)
				SSE2_RADIX_13_DFT(rad13_const, s1p10r,-0x080,-0x100,-0x180, 0x480, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, r30r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180)

				/*...and now do 13 radix-4 transforms...*/
				/* Inputs in SSE2 modes are temps 2*13*16 = 26*16 = 0x1a0 bytes apart. Notice how the add* indices after the first row repeat with period 4: */
		  #if !GCC_ASM_FULL_INLINE
			// Reorder blocks to yield sequentially increasing a-array offsets:
/* Block 01: */	add0 = &a[j1    ];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00r, 0x1a0)
/* Block 13: */	add3 = &a[j1+p04];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0cr, 0x1a0)
/* Block 12: */	add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0br, 0x1a0)
/* Block 11: */	add2 = &a[j1+p12];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0ar, 0x1a0)
/* Block 10: */	add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09r, 0x1a0)
/* Block 09: */	add3 = &a[j1+p20];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08r, 0x1a0)
/* Block 08: */	add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07r, 0x1a0)
/* Block 07: */	add2 = &a[j1+p28];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06r, 0x1a0)
/* Block 06: */	add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05r, 0x1a0)
/* Block 05: */	add3 = &a[j1+p36];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04r, 0x1a0)
/* Block 04: */	add1 = &a[j1+p40];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03r, 0x1a0)
/* Block 03: */	add2 = &a[j1+p44];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02r, 0x1a0)
/* Block 02: */	add0 = &a[j1+p48];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01r, 0x1a0)
		  #else
				add0 = &a[j1    ];
				SSE2_RADIX52_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,r00r);
		  #endif

	
			}	/* end for(j=_jstart; j < _jhi; j += 2) */

			jstart += nwt;
			jhi    += nwt;

			col += RADIX;
			co3 -= RADIX;
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
		thread_arg->cy00 = cy00->re;
		thread_arg->cy01 = cy00->im;
		thread_arg->cy02 = cy02->re;
		thread_arg->cy03 = cy02->im;
		thread_arg->cy04 = cy04->re;
		thread_arg->cy05 = cy04->im;
		thread_arg->cy06 = cy06->re;
		thread_arg->cy07 = cy06->im;
		thread_arg->cy08 = cy08->re;
		thread_arg->cy09 = cy08->im;
		thread_arg->cy10 = cy10->re;
		thread_arg->cy11 = cy10->im;
		thread_arg->cy12 = cy12->re;
		thread_arg->cy13 = cy12->im;
		thread_arg->cy14 = cy14->re;
		thread_arg->cy15 = cy14->im;
		thread_arg->cy16 = cy16->re;
		thread_arg->cy17 = cy16->im;
		thread_arg->cy18 = cy18->re;
		thread_arg->cy19 = cy18->im;
		thread_arg->cy20 = cy20->re;
		thread_arg->cy21 = cy20->im;
		thread_arg->cy22 = cy22->re;
		thread_arg->cy23 = cy22->im;
		thread_arg->cy24 = cy24->re;
		thread_arg->cy25 = cy24->im;
		thread_arg->cy26 = cy26->re;
		thread_arg->cy27 = cy26->im;
		thread_arg->cy28 = cy28->re;
		thread_arg->cy29 = cy28->im;
		thread_arg->cy30 = cy30->re;
		thread_arg->cy31 = cy30->im;
		thread_arg->cy32 = cy32->re;
		thread_arg->cy33 = cy32->im;
		thread_arg->cy34 = cy34->re;
		thread_arg->cy35 = cy34->im;
		thread_arg->cy36 = cy36->re;
		thread_arg->cy37 = cy36->im;
		thread_arg->cy38 = cy38->re;
		thread_arg->cy39 = cy38->im;
		thread_arg->cy40 = cy40->re;
		thread_arg->cy41 = cy40->im;
		thread_arg->cy42 = cy42->re;
		thread_arg->cy43 = cy42->im;
		thread_arg->cy44 = cy44->re;
		thread_arg->cy45 = cy44->im;
		thread_arg->cy46 = cy46->re;
		thread_arg->cy47 = cy46->im;
		thread_arg->cy48 = cy48->re;
		thread_arg->cy49 = cy48->im;
		thread_arg->cy50 = cy50->re;
		thread_arg->cy51 = cy50->im;

		if(max_err->re > max_err->im)
			maxerr = max_err->re;
		else
			maxerr = max_err->im;

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(thread_arg->maxerr < maxerr)
		{
			thread_arg->maxerr = maxerr;
		}

		return 0x0;
	}
#endif

