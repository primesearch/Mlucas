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

#ifdef MULTITHREAD
	#ifndef USE_PTHREAD
		#error Pthreads is only thread model currently supported!
	#endif
#endif

#ifdef USE_SSE2

	#define ERR_CHECK_ALL	/* #define this to do ROE checking of all convolution outputs, rather than just every 36th one */
	#if defined(USE_AVX) && !defined(ERR_CHECK_ALL)
		#error ERR_CHECK_ALL *required* for AVX-mode builds!
	#endif

	#define EPS 1e-10

	#if defined(COMPILER_TYPE_GCC) && (OS_BITS == 64)
		#define GCC_ASM_FULL_INLINE  1	// 0 to use small-macros to assemble radix-52 DFTs, 1 to inline fuse macros as a few big blobs of asm (64-bit only)
	#else
		#undef GCC_ASM_FULL_INLINE
	#endif

  // For Mersenne-mod we need (16 [SSE2] or 64 [AVX]) + 4 added slots for the half_arr lookup tables.
  // Add relevant number (half_arr_offset52 + RADIX) to get required value of radix52_creals_in_local_store:
  #ifdef USE_AVX
	const int half_arr_offset52 = 243;	// + RADIX = 295; Used for thread local-storage-integrity checking
	const int radix52_creals_in_local_store = 364;	// (half_arr_offset52 + RADIX) + 68 and round up to nearest multiple of 4
  #else
	const int half_arr_offset52 = 256;	// + RADIX = 308; Used for thread local-storage-integrity checking
	const int radix52_creals_in_local_store = 328;	// (half_arr_offset52 + RADIX) = 20 and round up to nearest multiple of 4
  #endif

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
			: "cc","memory","eax","ebx","ecx"		 /* Clobbered registers */\
			);\
		}

	  #else

		#ifdef USE_AVX

			#define SSE2_RADIX_13_DFT(Xcc, XI0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8,Xi9,XiA,XiB,XiC, XO0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8,Xo9,XoA,XoB,XoC)\
			{\
			__asm__ volatile (\
				"movq	%[__I0],%%rax							\n\t"\
				"movq	%[__cc],%%rbx							\n\t"\
				"movq	%[__O0],%%rcx							\n\t"\
			"/* xr-terms need 8 registers for each side: */		\n\t"\
				"vmovaps	%c[__i6](%%rax),%%ymm7				\n\t"\
				"vmovaps	%c[__i5](%%rax),%%ymm5				\n\t"\
				"vmovaps	%c[__i4](%%rax),%%ymm4				\n\t"\
				"vmovaps	%c[__i3](%%rax),%%ymm6				\n\t		/* yi-terms: */									\n\t"\
				"vmovaps	%c[__i2](%%rax),%%ymm3				\n\t			movq	%%rax,%%rdx							\n\t"\
				"vmovaps	%c[__i1](%%rax),%%ymm1				\n\t			addq	$0x020,%%rdx						\n\t"\
				"vmovaps	-0x020(%%rbx),%%ymm0				\n\t"\
				"vaddpd	%c[__i7](%%rax),%%ymm7,%%ymm7			\n\t			vmovaps	%c[__i1](%%rdx),%%ymm9 				\n\t"\
				"vaddpd	%c[__i8](%%rax),%%ymm5,%%ymm5			\n\t			vmovaps	%c[__i2](%%rdx),%%ymm10				\n\t"\
				"vaddpd	%c[__i9](%%rax),%%ymm4,%%ymm4			\n\t			vmovaps	%c[__i3](%%rdx),%%ymm13				\n\t"\
				"vaddpd	%c[__iA](%%rax),%%ymm6,%%ymm6			\n\t			vmovaps	%c[__i4](%%rdx),%%ymm11				\n\t"\
				"vaddpd	%c[__iB](%%rax),%%ymm3,%%ymm3			\n\t			vmovaps	%c[__i5](%%rdx),%%ymm12				\n\t"\
				"vaddpd	%c[__iC](%%rax),%%ymm1,%%ymm1			\n\t			vmovaps	%c[__i6](%%rdx),%%ymm14				\n\t"\
				"vsubpd	%%ymm5,%%ymm1,%%ymm1					\n\t			vsubpd	%c[__iC](%%rdx),%%ymm9 ,%%ymm9 		\n\t"\
				"vmulpd	%%ymm0,%%ymm5,%%ymm5					\n\t			vsubpd	%c[__iB](%%rdx),%%ymm10,%%ymm10		\n\t"\
				"vsubpd	%%ymm6,%%ymm3,%%ymm3					\n\t			vsubpd	%c[__iA](%%rdx),%%ymm13,%%ymm13		\n\t"\
				"vmulpd	%%ymm0,%%ymm6,%%ymm6					\n\t			vsubpd	%c[__i9](%%rdx),%%ymm11,%%ymm11		\n\t"\
				"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t			vsubpd	%c[__i8](%%rdx),%%ymm12,%%ymm12		\n\t"\
				"vmulpd	%%ymm0,%%ymm7,%%ymm7					\n\t			vsubpd	%c[__i7](%%rdx),%%ymm14,%%ymm14		\n\t"\
				"vaddpd	%%ymm1,%%ymm5,%%ymm5					\n\t			vmovaps	%%ymm9 ,%%ymm15						\n\t"\
				"vaddpd	%%ymm3,%%ymm6,%%ymm6					\n\t			vmovaps	%%ymm10,%%ymm8 						\n\t"\
				"vaddpd	%%ymm4,%%ymm7,%%ymm7					\n\t			vsubpd	%%ymm11,%%ymm15,%%ymm15				\n\t"\
				"vmovaps	%%ymm5,%%ymm2						\n\t			vaddpd	%%ymm12,%%ymm8 ,%%ymm8 				\n\t"\
				"vaddpd	%%ymm6,%%ymm2,%%ymm2					\n\t			vaddpd	%%ymm13,%%ymm15,%%ymm15				\n\t"\
				"vaddpd	%%ymm7,%%ymm2,%%ymm2					\n\t			vaddpd	%%ymm14,%%ymm8 ,%%ymm8 				\n\t"\
				"vmovaps	(%%rax),%%ymm0						\n\t			vaddpd	%%ymm11,%%ymm9 ,%%ymm9 				\n\t"\
				"vaddpd	%%ymm2,%%ymm0,%%ymm0					\n\t			vaddpd	%%ymm11,%%ymm13,%%ymm13				\n\t"\
				"vmulpd	(%%rbx),%%ymm2,%%ymm2					\n\t			vsubpd	%%ymm12,%%ymm10,%%ymm10				\n\t"\
				"vmovaps	%%ymm0,(%%rcx)						\n\t			vsubpd	%%ymm12,%%ymm14,%%ymm14				\n\t"\
				"vaddpd	%%ymm2,%%ymm0,%%ymm0					\n\t			vmovaps	%%ymm8 ,%%ymm12						\n\t"\
				"vmovaps	0x020(%%rbx),%%ymm2					\n\t			vmovaps	%%ymm15,%%ymm11						\n\t"\
				"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t			vmulpd	0x080(%%rbx),%%ymm8 ,%%ymm8 		\n\t"\
				"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t			vmulpd	0x080(%%rbx),%%ymm15,%%ymm15		\n\t"\
				"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t			vmulpd	0x060(%%rbx),%%ymm12,%%ymm12		\n\t"\
				"vmovaps	%%ymm6,%c[__o1](%%rcx)				\n\t			vmulpd	0x060(%%rbx),%%ymm11,%%ymm11		\n\t"\
				"vmovaps	%%ymm5,%c[__o2](%%rcx)				\n\t			vsubpd	%%ymm15,%%ymm12,%%ymm12				\n\t"\
				"vmovaps	%%ymm7,%c[__o3](%%rcx)				\n\t			vaddpd	%%ymm8 ,%%ymm11,%%ymm11				\n\t"\
				"vmovaps	0x100(%%rbx),%%ymm2					\n\t			vmovaps	%%ymm12,%c[__oC](%%rcx)				\n\t"\
				"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t			vmovaps	%%ymm9 ,%%ymm8 						\n\t"\
				"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t			vmovaps	%%ymm13,%%ymm12						\n\t"\
				"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t			vaddpd	%%ymm10,%%ymm8 ,%%ymm8 				\n\t"\
				"vaddpd	%%ymm0,%%ymm5,%%ymm5					\n\t			vaddpd	%%ymm14,%%ymm12,%%ymm12				\n\t"\
				"vaddpd	%%ymm0,%%ymm7,%%ymm7					\n\t			vmovaps	%%ymm8 ,%c[__oB](%%rcx)				\n\t"\
				"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t			vmovaps	%%ymm12,%%ymm15						\n\t"\
				"vaddpd	%c[__o1](%%rcx),%%ymm5,%%ymm5			\n\t			vmulpd	0x160(%%rbx),%%ymm12,%%ymm12		\n\t"\
				"vaddpd	%c[__o2](%%rcx),%%ymm7,%%ymm7			\n\t			vmulpd	0x1c0(%%rbx),%%ymm8 ,%%ymm8 		\n\t"\
				"vaddpd	%c[__o3](%%rcx),%%ymm6,%%ymm6			\n\t			vaddpd	%c[__oB](%%rcx),%%ymm12,%%ymm12		\n\t"\
				"vmovaps	0x140(%%rbx),%%ymm2					\n\t			vaddpd	%%ymm15,%%ymm8 ,%%ymm8 				\n\t"\
				"vmovaps	0x120(%%rbx),%%ymm0					\n\t			vmulpd	0x0a0(%%rbx),%%ymm12,%%ymm12		\n\t"\
				"vmovaps	%%ymm4,%c[__o1](%%rcx)				\n\t			vmulpd	0x0a0(%%rbx),%%ymm8 ,%%ymm8 		\n\t"\
				"vmovaps	%%ymm1,%c[__o2](%%rcx)				\n\t			vmovaps	%%ymm10,%c[__oB](%%rcx)				\n\t"\
				"vmovaps	%%ymm3,%c[__o3](%%rcx)				\n\t			vmovaps	%%ymm14,%%ymm15						\n\t"\
				"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t			vmulpd	0x180(%%rbx),%%ymm14,%%ymm14		\n\t"\
				"vmulpd	%%ymm2,%%ymm1,%%ymm1					\n\t			vmulpd	0x1e0(%%rbx),%%ymm10,%%ymm10		\n\t"\
				"vmulpd	%%ymm2,%%ymm3,%%ymm3					\n\t			vaddpd	%c[__oB](%%rcx),%%ymm14,%%ymm14		\n\t"\
				"vaddpd	%c[__o3](%%rcx),%%ymm4,%%ymm4			\n\t			vaddpd	%%ymm15,%%ymm10,%%ymm10				\n\t"\
				"vsubpd	%c[__o1](%%rcx),%%ymm1,%%ymm1			\n\t			vmulpd	0x0c0(%%rbx),%%ymm14,%%ymm14		\n\t"\
				"vaddpd	%c[__o2](%%rcx),%%ymm3,%%ymm3			\n\t			vmulpd	0x0c0(%%rbx),%%ymm10,%%ymm10		\n\t"\
				"vmulpd	%%ymm0,%%ymm4,%%ymm4					\n\t			vaddpd	%%ymm12,%%ymm14,%%ymm14				\n\t"\
				"vmulpd	%%ymm0,%%ymm1,%%ymm1					\n\t			vaddpd	%%ymm8 ,%%ymm10,%%ymm10				\n\t"\
				"vmulpd	%%ymm0,%%ymm3,%%ymm3					\n\t			vmovaps	%%ymm9 ,%c[__oB](%%rcx)				\n\t"\
				"vmovaps	0x040(%%rbx),%%ymm2					\n\t			vmovaps	%%ymm13,%%ymm15						\n\t"\
				"vaddpd	%c[__o2](%%rcx),%%ymm4,%%ymm4			\n\t			vmulpd	0x1a0(%%rbx),%%ymm13,%%ymm13		\n\t"\
				"vsubpd	%c[__o3](%%rcx),%%ymm1,%%ymm1			\n\t			vmulpd	0x200(%%rbx),%%ymm9 ,%%ymm9 		\n\t"\
				"vsubpd	%c[__o1](%%rcx),%%ymm3,%%ymm3			\n\t			vaddpd	%c[__oB](%%rcx),%%ymm13,%%ymm13		\n\t"\
				"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t			vaddpd	%%ymm15,%%ymm9 ,%%ymm9 				\n\t"\
				"vmulpd	%%ymm2,%%ymm1,%%ymm1					\n\t			vmulpd	0x0e0(%%rbx),%%ymm13,%%ymm13		\n\t"\
				"vmulpd	%%ymm2,%%ymm3,%%ymm3					\n\t			vmulpd	0x0e0(%%rbx),%%ymm9 ,%%ymm9 		\n\t"\
				"vmovaps	-0x020(%%rbx),%%ymm0				\n\t			vaddpd	%%ymm12,%%ymm13,%%ymm13				\n\t"\
			"/* Spill into destination outputs: */				\n\t			vaddpd	%%ymm8 ,%%ymm9 ,%%ymm9 				\n\t"\
				"vsubpd	%%ymm4,%%ymm6,%%ymm6					\n\t			vmovaps	%c[__oC](%%rcx),%%ymm12				\n\t"\
				"vmulpd	%%ymm0,%%ymm4,%%ymm4					\n\t			vmovaps	%%ymm11,%%ymm15						\n\t"\
				"vsubpd	%%ymm1,%%ymm7,%%ymm7					\n\t			vmovaps	%%ymm12,%%ymm8 						\n\t"\
				"vmulpd	%%ymm0,%%ymm1,%%ymm1					\n\t			vaddpd	%%ymm14,%%ymm15,%%ymm15				\n\t"\
				"vsubpd	%%ymm3,%%ymm5,%%ymm5					\n\t			vaddpd	%%ymm13,%%ymm8 ,%%ymm8 				\n\t"\
				"vmulpd	%%ymm0,%%ymm3,%%ymm3					\n\t			vsubpd	%%ymm11,%%ymm14,%%ymm14				\n\t"\
			"/*	vmovaps	%%ymm7,%c[__o3](%%rcx)	*/				\n\t			vsubpd	%%ymm12,%%ymm13,%%ymm13				\n\t"\
			"/*	vmovaps	%%ymm5,%c[__o8](%%rcx)	*/				\n\t			vaddpd	%%ymm10,%%ymm14,%%ymm14				\n\t"\
				"vaddpd	%%ymm6,%%ymm4,%%ymm4					\n\t			vaddpd	%%ymm9 ,%%ymm13,%%ymm13				\n\t"\
				"vaddpd	%%ymm7,%%ymm1,%%ymm1					\n\t			vaddpd	%%ymm11,%%ymm10,%%ymm10				\n\t"\
				"vaddpd	%%ymm5,%%ymm3,%%ymm3					\n\t			vaddpd	%%ymm9 ,%%ymm12,%%ymm12				\n\t"\
			"/*	vmovaps	%%ymm6,%c[__o4](%%rcx)	*/				\n\t"\
			"/*	vmovaps	%%ymm4,%c[__o6](%%rcx)	*/				\n\t"\
			"/*	vmovaps	%%ymm1,%c[__o2](%%rcx)	*/				\n\t"\
			"/*	vmovaps	%%ymm3,%c[__o1](%%rcx)	*/				\n\t"\
		  "/* yi-data in ymm8,10,12,13,14,15; ymm0,2,9,11 free: */	\n\t"\
			"/*	vmovaps	%c[__o6](%%rcx),%%ymm4	*/				\n\t"\
				"vmovaps	%%ymm4 ,%%ymm0						\n\t"\
				"vsubpd	%%ymm15,%%ymm4,%%ymm4					\n\t		/*	vmovaps	%c[__o3](%%rcx),%%ymm7	*/			\n\t"\
				"vaddpd	%%ymm15,%%ymm0,%%ymm0					\n\t			vmovaps	%%ymm7 ,%%ymm9 						\n\t"\
				"vmovaps	%%ymm4 ,%c[__o6](%%rcx)				\n\t			vsubpd	%%ymm8 ,%%ymm7 ,%%ymm7 				\n\t"\
				"vmovaps	%%ymm0 ,%c[__o7](%%rcx)				\n\t			vaddpd	%%ymm8 ,%%ymm9 ,%%ymm9 				\n\t"\
			"/*	vmovaps	%c[__o8](%%rcx),%%ymm5	*/				\n\t			vmovaps	%%ymm7 ,%c[__o3](%%rcx)				\n\t"\
				"vmovaps	%%ymm5 ,%%ymm0						\n\t			vmovaps	%%ymm9 ,%c[__oA](%%rcx)				\n\t"\
				"vsubpd	%%ymm14,%%ymm5,%%ymm5					\n\t		/*	vmovaps	%c[__o4](%%rcx),%%ymm6	*/			\n\t"\
				"vaddpd	%%ymm14,%%ymm0,%%ymm0					\n\t			vmovaps	%%ymm6 ,%%ymm9 						\n\t"\
				"vmovaps	%%ymm5 ,%c[__o8](%%rcx)				\n\t			vsubpd	%%ymm13,%%ymm6 ,%%ymm6 				\n\t"\
				"vmovaps	%%ymm0 ,%c[__o5](%%rcx)				\n\t			vaddpd	%%ymm13,%%ymm9 ,%%ymm9 				\n\t"\
			"/*	vmovaps	%c[__o2](%%rcx),%%ymm1	*/				\n\t			vmovaps	%%ymm6 ,%c[__o4](%%rcx)				\n\t"\
				"vmovaps	%%ymm1 ,%%ymm0						\n\t			vmovaps	%%ymm9 ,%c[__o9](%%rcx)				\n\t"\
				"vsubpd	%%ymm10,%%ymm1,%%ymm1					\n\t		/*	vmovaps	%c[__o1](%%rcx),%%ymm3	*/			\n\t"\
				"vaddpd	%%ymm10,%%ymm0,%%ymm0					\n\t			vmovaps	%%ymm3 ,%%ymm9 						\n\t"\
				"vmovaps	%%ymm1 ,%c[__o2](%%rcx)				\n\t			vsubpd	%%ymm12,%%ymm3 ,%%ymm3 				\n\t"\
				"vmovaps	%%ymm0 ,%c[__oB](%%rcx)				\n\t			vaddpd	%%ymm12,%%ymm9 ,%%ymm9 				\n\t"\
				"/***************/								\n\t			vmovaps	%%ymm3 ,%c[__o1](%%rcx)				\n\t"\
				"/* IMAG PARTS: */								\n\t			vmovaps	%%ymm9 ,%c[__oC](%%rcx)				\n\t"\
				"/***************/								\n\t"\
			"/* xi-terms need 8 registers for each	side: */	\n\t"\
				"addq	$0x020,%%rax		 					\n\t"\
				"addq	$0x020,%%rcx							\n\t"\
				"vmovaps	%c[__i6](%%rax),%%ymm7				\n\t"\
				"vmovaps	%c[__i5](%%rax),%%ymm5				\n\t"\
				"vmovaps	%c[__i4](%%rax),%%ymm4				\n\t"\
				"vmovaps	%c[__i3](%%rax),%%ymm6				\n\t		/* yi-terms: */									\n\t"\
				"vmovaps	%c[__i2](%%rax),%%ymm3				\n\t"\
				"vmovaps	%c[__i1](%%rax),%%ymm1				\n\t			subq	$0x020,%%rdx						\n\t"\
				"vmovaps	-0x020(%%rbx),%%ymm0				\n\t"\
				"vaddpd	%c[__i7](%%rax),%%ymm7,%%ymm7			\n\t			vmovaps	%c[__i1](%%rdx),%%ymm9 				\n\t"\
				"vaddpd	%c[__i8](%%rax),%%ymm5,%%ymm5			\n\t			vmovaps	%c[__i2](%%rdx),%%ymm10				\n\t"\
				"vaddpd	%c[__i9](%%rax),%%ymm4,%%ymm4			\n\t			vmovaps	%c[__i3](%%rdx),%%ymm13				\n\t"\
				"vaddpd	%c[__iA](%%rax),%%ymm6,%%ymm6			\n\t			vmovaps	%c[__i4](%%rdx),%%ymm11				\n\t"\
				"vaddpd	%c[__iB](%%rax),%%ymm3,%%ymm3			\n\t			vmovaps	%c[__i5](%%rdx),%%ymm12				\n\t"\
				"vaddpd	%c[__iC](%%rax),%%ymm1,%%ymm1			\n\t			vmovaps	%c[__i6](%%rdx),%%ymm14				\n\t"\
				"vsubpd	%%ymm5,%%ymm1,%%ymm1					\n\t			vsubpd	%c[__iC](%%rdx),%%ymm9 ,%%ymm9 		\n\t"\
				"vmulpd	%%ymm0,%%ymm5,%%ymm5					\n\t			vsubpd	%c[__iB](%%rdx),%%ymm10,%%ymm10		\n\t"\
				"vsubpd	%%ymm6,%%ymm3,%%ymm3					\n\t			vsubpd	%c[__iA](%%rdx),%%ymm13,%%ymm13		\n\t"\
				"vmulpd	%%ymm0,%%ymm6,%%ymm6					\n\t			vsubpd	%c[__i9](%%rdx),%%ymm11,%%ymm11		\n\t"\
				"vsubpd	%%ymm7,%%ymm4,%%ymm4					\n\t			vsubpd	%c[__i8](%%rdx),%%ymm12,%%ymm12		\n\t"\
				"vmulpd	%%ymm0,%%ymm7,%%ymm7					\n\t			vsubpd	%c[__i7](%%rdx),%%ymm14,%%ymm14		\n\t"\
				"vaddpd	%%ymm1,%%ymm5,%%ymm5					\n\t			vmovaps	%%ymm9 ,%%ymm15						\n\t"\
				"vaddpd	%%ymm3,%%ymm6,%%ymm6					\n\t			vmovaps	%%ymm10,%%ymm8 						\n\t"\
				"vaddpd	%%ymm4,%%ymm7,%%ymm7					\n\t			vsubpd	%%ymm11,%%ymm15,%%ymm15				\n\t"\
				"vmovaps	%%ymm5,%%ymm2						\n\t			vaddpd	%%ymm12,%%ymm8 ,%%ymm8 				\n\t"\
				"vaddpd	%%ymm6,%%ymm2,%%ymm2					\n\t			vaddpd	%%ymm13,%%ymm15,%%ymm15				\n\t"\
				"vaddpd	%%ymm7,%%ymm2,%%ymm2					\n\t			vaddpd	%%ymm14,%%ymm8 ,%%ymm8 				\n\t"\
				"vmovaps	(%%rax),%%ymm0						\n\t			vaddpd	%%ymm11,%%ymm9 ,%%ymm9 				\n\t"\
				"vaddpd	%%ymm2,%%ymm0,%%ymm0					\n\t			vaddpd	%%ymm11,%%ymm13,%%ymm13				\n\t"\
				"vmulpd	(%%rbx),%%ymm2,%%ymm2					\n\t			vsubpd	%%ymm12,%%ymm10,%%ymm10				\n\t"\
				"vmovaps	%%ymm0,(%%rcx)						\n\t			vsubpd	%%ymm12,%%ymm14,%%ymm14				\n\t"\
				"vaddpd	%%ymm2,%%ymm0,%%ymm0					\n\t			vmovaps	%%ymm8 ,%%ymm12						\n\t"\
				"vmovaps	0x020(%%rbx),%%ymm2					\n\t			vmovaps	%%ymm15,%%ymm11						\n\t"\
				"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t			vmulpd	0x080(%%rbx),%%ymm8 ,%%ymm8 		\n\t"\
				"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t			vmulpd	0x080(%%rbx),%%ymm15,%%ymm15		\n\t"\
				"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t			vmulpd	0x060(%%rbx),%%ymm12,%%ymm12		\n\t"\
				"vmovaps	%%ymm6,%c[__oC](%%rcx)				\n\t			vmulpd	0x060(%%rbx),%%ymm11,%%ymm11		\n\t"\
				"vmovaps	%%ymm5,%c[__oB](%%rcx)				\n\t			vsubpd	%%ymm15,%%ymm12,%%ymm12				\n\t"\
				"vmovaps	%%ymm7,%c[__oA](%%rcx)				\n\t			vaddpd	%%ymm8 ,%%ymm11,%%ymm11				\n\t"\
				"vmovaps	0x100(%%rbx),%%ymm2					\n\t			vmovaps	%%ymm12,%c[__o1](%%rcx)				\n\t"\
				"vmulpd	%%ymm2,%%ymm5,%%ymm5					\n\t			vmovaps	%%ymm9 ,%%ymm8 						\n\t"\
				"vmulpd	%%ymm2,%%ymm7,%%ymm7					\n\t			vmovaps	%%ymm13,%%ymm12						\n\t"\
				"vmulpd	%%ymm2,%%ymm6,%%ymm6					\n\t			vaddpd	%%ymm10,%%ymm8 ,%%ymm8 				\n\t"\
				"vaddpd	%%ymm0,%%ymm5,%%ymm5					\n\t			vaddpd	%%ymm14,%%ymm12,%%ymm12				\n\t"\
				"vaddpd	%%ymm0,%%ymm7,%%ymm7					\n\t			vmovaps	%%ymm8 ,%c[__o2](%%rcx)				\n\t"\
				"vaddpd	%%ymm0,%%ymm6,%%ymm6					\n\t			vmovaps	%%ymm12,%%ymm15						\n\t"\
				"vaddpd	%c[__oC](%%rcx),%%ymm5,%%ymm5			\n\t			vmulpd	0x160(%%rbx),%%ymm12,%%ymm12		\n\t"\
				"vaddpd	%c[__oB](%%rcx),%%ymm7,%%ymm7			\n\t			vmulpd	0x1c0(%%rbx),%%ymm8 ,%%ymm8 		\n\t"\
				"vaddpd	%c[__oA](%%rcx),%%ymm6,%%ymm6			\n\t			vaddpd	%c[__o2](%%rcx),%%ymm12,%%ymm12		\n\t"\
				"vmovaps	0x140(%%rbx),%%ymm2					\n\t			vaddpd	%%ymm15,%%ymm8 ,%%ymm8 				\n\t"\
				"vmovaps	0x120(%%rbx),%%ymm0					\n\t			vmulpd	0x0a0(%%rbx),%%ymm12,%%ymm12		\n\t"\
				"vmovaps	%%ymm4,%c[__oC](%%rcx)				\n\t			vmulpd	0x0a0(%%rbx),%%ymm8 ,%%ymm8 		\n\t"\
				"vmovaps	%%ymm1,%c[__oB](%%rcx)				\n\t			vmovaps	%%ymm10,%c[__o2](%%rcx)				\n\t"\
				"vmovaps	%%ymm3,%c[__oA](%%rcx)				\n\t			vmovaps	%%ymm14,%%ymm15						\n\t"\
				"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t			vmulpd	0x180(%%rbx),%%ymm14,%%ymm14		\n\t"\
				"vmulpd	%%ymm2,%%ymm1,%%ymm1					\n\t			vmulpd	0x1e0(%%rbx),%%ymm10,%%ymm10		\n\t"\
				"vmulpd	%%ymm2,%%ymm3,%%ymm3					\n\t			vaddpd	%c[__o2](%%rcx),%%ymm14,%%ymm14		\n\t"\
				"vaddpd	%c[__oA](%%rcx),%%ymm4,%%ymm4			\n\t			vaddpd	%%ymm15,%%ymm10,%%ymm10				\n\t"\
				"vsubpd	%c[__oC](%%rcx),%%ymm1,%%ymm1			\n\t			vmulpd	0x0c0(%%rbx),%%ymm14,%%ymm14		\n\t"\
				"vaddpd	%c[__oB](%%rcx),%%ymm3,%%ymm3			\n\t			vmulpd	0x0c0(%%rbx),%%ymm10,%%ymm10		\n\t"\
				"vmulpd	%%ymm0,%%ymm4,%%ymm4					\n\t			vaddpd	%%ymm12,%%ymm14,%%ymm14				\n\t"\
				"vmulpd	%%ymm0,%%ymm1,%%ymm1					\n\t			vaddpd	%%ymm8 ,%%ymm10,%%ymm10				\n\t"\
				"vmulpd	%%ymm0,%%ymm3,%%ymm3					\n\t			vmovaps	%%ymm9 ,%c[__o2](%%rcx)				\n\t"\
				"vmovaps	0x040(%%rbx),%%ymm2					\n\t			vmovaps	%%ymm13,%%ymm15						\n\t"\
				"vaddpd	%c[__oB](%%rcx),%%ymm4,%%ymm4			\n\t			vmulpd	0x1a0(%%rbx),%%ymm13,%%ymm13		\n\t"\
				"vsubpd	%c[__oA](%%rcx),%%ymm1,%%ymm1			\n\t			vmulpd	0x200(%%rbx),%%ymm9 ,%%ymm9 		\n\t"\
				"vsubpd	%c[__oC](%%rcx),%%ymm3,%%ymm3			\n\t			vaddpd	%c[__o2](%%rcx),%%ymm13,%%ymm13		\n\t"\
				"vmulpd	%%ymm2,%%ymm4,%%ymm4					\n\t			vaddpd	%%ymm15,%%ymm9 ,%%ymm9 				\n\t"\
				"vmulpd	%%ymm2,%%ymm1,%%ymm1					\n\t			vmulpd	0x0e0(%%rbx),%%ymm13,%%ymm13		\n\t"\
				"vmulpd	%%ymm2,%%ymm3,%%ymm3					\n\t			vmulpd	0x0e0(%%rbx),%%ymm9 ,%%ymm9 		\n\t"\
				"vmovaps	-0x020(%%rbx),%%ymm0				\n\t			vaddpd	%%ymm12,%%ymm13,%%ymm13				\n\t"\
			"/* Spill into destination outputs: */				\n\t			vaddpd	%%ymm8 ,%%ymm9 ,%%ymm9 				\n\t"\
				"vsubpd	%%ymm4,%%ymm6,%%ymm6					\n\t			vmovaps	%c[__o1](%%rcx),%%ymm12				\n\t"\
				"vmulpd	%%ymm0,%%ymm4,%%ymm4					\n\t			vmovaps	%%ymm11,%%ymm15						\n\t"\
				"vsubpd	%%ymm1,%%ymm7,%%ymm7					\n\t			vmovaps	%%ymm12,%%ymm8 						\n\t"\
				"vmulpd	%%ymm0,%%ymm1,%%ymm1					\n\t			vaddpd	%%ymm14,%%ymm15,%%ymm15				\n\t"\
				"vsubpd	%%ymm3,%%ymm5,%%ymm5					\n\t			vaddpd	%%ymm13,%%ymm8 ,%%ymm8 				\n\t"\
				"vmulpd	%%ymm0,%%ymm3,%%ymm3					\n\t			vsubpd	%%ymm11,%%ymm14,%%ymm14				\n\t"\
			"/*	vmovaps	%%ymm7,%c[__oA](%%rcx)	*/				\n\t			vsubpd	%%ymm12,%%ymm13,%%ymm13				\n\t"\
			"/*	vmovaps	%%ymm5,%c[__o5](%%rcx)	*/				\n\t			vaddpd	%%ymm10,%%ymm14,%%ymm14				\n\t"\
				"vaddpd	%%ymm6,%%ymm4,%%ymm4					\n\t			vaddpd	%%ymm9 ,%%ymm13,%%ymm13				\n\t"\
				"vaddpd	%%ymm7,%%ymm1,%%ymm1					\n\t			vaddpd	%%ymm11,%%ymm10,%%ymm10				\n\t"\
				"vaddpd	%%ymm5,%%ymm3,%%ymm3					\n\t			vaddpd	%%ymm9 ,%%ymm12,%%ymm12				\n\t"\
			"/*	vmovaps	%%ymm6,%c[__o9](%%rcx)	*/				\n\t"\
			"/*	vmovaps	%%ymm4,%c[__o7](%%rcx)	*/				\n\t"\
			"/*	vmovaps	%%ymm1,%c[__oB](%%rcx)	*/				\n\t"\
			"/*	vmovaps	%%ymm3,%c[__oC](%%rcx)	*/				\n\t"\
		  "/* yi-data in ymm8,10,12,13,14,15; ymm0,2,9,11 free: */	\n\t"\
			"/*	vmovaps	%c[__o7](%%rcx),%%ymm4	*/				\n\t"\
				"vmovaps	%%ymm4 ,%%ymm0						\n\t"\
				"vsubpd	%%ymm15,%%ymm4,%%ymm4					\n\t		/*	vmovaps	%c[__oA](%%rcx),%%ymm7	*/			\n\t"\
				"vaddpd	%%ymm15,%%ymm0,%%ymm0					\n\t			vmovaps	%%ymm7 ,%%ymm9 						\n\t"\
				"vmovaps	%%ymm4 ,%c[__o7](%%rcx)				\n\t			vsubpd	%%ymm8 ,%%ymm7 ,%%ymm7 				\n\t"\
				"vmovaps	%%ymm0 ,%c[__o6](%%rcx)				\n\t			vaddpd	%%ymm8 ,%%ymm9 ,%%ymm9 				\n\t"\
			"/*	vmovaps	%c[__o5](%%rcx),%%ymm5	*/				\n\t			vmovaps	%%ymm7 ,%c[__oA](%%rcx)				\n\t"\
				"vmovaps	%%ymm5 ,%%ymm0						\n\t			vmovaps	%%ymm9 ,%c[__o3](%%rcx)				\n\t"\
				"vsubpd	%%ymm14,%%ymm5,%%ymm5					\n\t		/*	vmovaps	%c[__o9](%%rcx),%%ymm6	*/			\n\t"\
				"vaddpd	%%ymm14,%%ymm0,%%ymm0					\n\t			vmovaps	%%ymm6 ,%%ymm9 						\n\t"\
				"vmovaps	%%ymm5 ,%c[__o5](%%rcx)				\n\t			vsubpd	%%ymm13,%%ymm6 ,%%ymm6 				\n\t"\
				"vmovaps	%%ymm0 ,%c[__o8](%%rcx)				\n\t			vaddpd	%%ymm13,%%ymm9 ,%%ymm9 				\n\t"\
			"/*	vmovaps	%c[__oB](%%rcx),%%ymm1	*/				\n\t			vmovaps	%%ymm6 ,%c[__o9](%%rcx)				\n\t"\
				"vmovaps	%%ymm1 ,%%ymm0						\n\t			vmovaps	%%ymm9 ,%c[__o4](%%rcx)				\n\t"\
				"vsubpd	%%ymm10,%%ymm1,%%ymm1					\n\t		/*	vmovaps	%c[__oC](%%rcx),%%ymm3	*/			\n\t"\
				"vaddpd	%%ymm10,%%ymm0,%%ymm0					\n\t			vmovaps	%%ymm3 ,%%ymm9 						\n\t"\
				"vmovaps	%%ymm1 ,%c[__oB](%%rcx)				\n\t			vsubpd	%%ymm12,%%ymm3 ,%%ymm3 				\n\t"\
				"vmovaps	%%ymm0 ,%c[__o2](%%rcx)				\n\t			vaddpd	%%ymm12,%%ymm9 ,%%ymm9 				\n\t"\
				"																vmovaps	%%ymm3 ,%c[__oC](%%rcx)				\n\t"\
				"																vmovaps	%%ymm9 ,%c[__o1](%%rcx)				\n\t"\
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
				: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
				);\
			}

		#else	// 64-bit SSE2:

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
				: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	 /* Clobbered registers */\
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
				: "cc","memory","rax","rbx","rcx","rdx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
				);\
			}

		  #endif // USE_64BIT_ASM_STYLE

		#endif	// AVX and 64-bit SSE2

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

#ifdef USE_PTHREAD

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
		vec_dbl *s1p00r;
		vec_dbl *half_arr;

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
	int NDIVR,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k,khi,l,outer,nbytes;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40,p44,p48;
	static double radix_inv, n2inv;
	double scale, dtmp, maxerr = 0.0;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3;
  #ifdef USE_AVX
	static struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
  #else
	int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
  #endif
	double
	t00r,t01r,t02r,t03r,t04r,t05r,t06r,t07r,t08r,t09r,t0ar,t0br,t0cr,
	t10r,t11r,t12r,t13r,t14r,t15r,t16r,t17r,t18r,t19r,t1ar,t1br,t1cr,
	t20r,t21r,t22r,t23r,t24r,t25r,t26r,t27r,t28r,t29r,t2ar,t2br,t2cr,
	t30r,t31r,t32r,t33r,t34r,t35r,t36r,t37r,t38r,t39r,t3ar,t3br,t3cr;

#ifdef USE_SSE2

  #if !(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))
	#error SSE2 code not supported for this compiler!
  #endif

	static int cslots_in_local_store;
	static vec_dbl *sc_arr = 0x0, *sc_ptr;
	static uint64 *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;

  #ifdef MULTITHREAD
	static vec_dbl *__r0;	// Base address for discrete per-thread local stores
  #else
	double *add0, *add1, *add2, *add3;	/* Addresses into array sections */
  #endif

	static int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09
		,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19
		,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27,*bjmodn28,*bjmodn29
		,*bjmodn30,*bjmodn31,*bjmodn32,*bjmodn33,*bjmodn34,*bjmodn35,*bjmodn36,*bjmodn37,*bjmodn38,*bjmodn39
		,*bjmodn40,*bjmodn41,*bjmodn42,*bjmodn43,*bjmodn44,*bjmodn45,*bjmodn46,*bjmodn47,*bjmodn48,*bjmodn49
		,*bjmodn50,*bjmodn51;
	const double crnd = 3.0*0x4000000*0x2000000;
	struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
	static vec_dbl *two,*rad13_const, *max_err, *sse2_rnd, *half_arr, *tmp,*tm2     /* rad13_const needs 18*16 bytes allocated */
		,*r00r,*r01r,*r02r,*r03r,*r04r,*r05r,*r06r,*r07r,*r08r,*r09r,*r0ar,*r0br,*r0cr
		,*r10r,*r11r,*r12r,*r13r,*r14r,*r15r,*r16r,*r17r,*r18r,*r19r,*r1ar,*r1br,*r1cr
		,*r20r,*r21r,*r22r,*r23r,*r24r,*r25r,*r26r,*r27r,*r28r,*r29r,*r2ar,*r2br,*r2cr
		,*r30r,*r31r,*r32r,*r33r,*r34r,*r35r,*r36r,*r37r,*r38r,*r39r,*r3ar,*r3br,*r3cr
		,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p0ar,*s1p0br,*s1p0cr
		,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p1ar,*s1p1br,*s1p1cr
		,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r,*s1p28r,*s1p29r,*s1p2ar,*s1p2br,*s1p2cr
		,*s1p30r,*s1p31r,*s1p32r,*s1p33r,*s1p34r,*s1p35r,*s1p36r,*s1p37r,*s1p38r,*s1p39r,*s1p3ar,*s1p3br,*s1p3cr;
	static vec_dbl
		*cy00,*cy04,*cy08,*cy12,*cy16,*cy20,*cy24,*cy28,*cy32,*cy36,*cy40,*cy44,*cy48;
  #ifndef USE_AVX
	static vec_dbl
		*cy02,*cy06,*cy10,*cy14,*cy18,*cy22,*cy26,*cy30,*cy34,*cy38,*cy42,*cy46,*cy50;
  #endif

#endif

#ifdef MULTITHREAD

	static struct cy_thread_data_t *tdat = 0x0;
	// Threadpool-based dispatch stuff:
	static int main_work_units = 0, pool_work_units = 0;
	static struct threadpool *tpool = 0x0;
	static int task_is_blocking = TRUE;
	static thread_control_t thread_control = {0,0,0};
	// First 3 subfields same for all threads, 4th provides thread-specifc data, will be inited at thread dispatch:
	static task_control_t   task_control = {NULL, (void*)cy52_process_chunk, NULL, 0x0};

#elif !defined(USE_SSE2)

	// Vars needed in scalar mode only:
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	int m,m2,
		bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,
		bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,
		bjmodn26,bjmodn27,bjmodn28,bjmodn29,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,bjmodn36,bjmodn37,bjmodn38,
		bjmodn39,bjmodn40,bjmodn41,bjmodn42,bjmodn43,bjmodn44,bjmodn45,bjmodn46,bjmodn47,bjmodn48,bjmodn49,bjmodn50,bjmodn51;
	double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	double rt,it,temp,frac,
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
		}
	#endif

	#ifdef USE_SSE2

		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		// Use double-complex type size (16 bytes) to alloc a block of local storage
		// consisting of 88 dcomplex and (12+RADIX/2) uint64 element slots per thread
		// (Add as many padding elts to the latter as needed to make it a multiple of 4):
		cslots_in_local_store = radix52_creals_in_local_store + (((12+RADIX/2)/2 + 3) & ~0x3);
		sc_arr = ALLOC_VEC_DBL(sc_arr, cslots_in_local_store*CY_THREADS);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_VEC_DBL(sc_arr);
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
		s1p01r = sc_ptr + 0x02;    s1p11r = two + 0x02;    s1p21r = max_err + 0x02;    s1p31r = sse2_rnd + 0x02;
		s1p02r = sc_ptr + 0x04;    s1p12r = two + 0x04;    s1p22r = max_err + 0x04;    s1p32r = sse2_rnd + 0x04;
		s1p03r = sc_ptr + 0x06;    s1p13r = two + 0x06;    s1p23r = max_err + 0x06;    s1p33r = sse2_rnd + 0x06;
		s1p04r = sc_ptr + 0x08;    s1p14r = two + 0x08;    s1p24r = max_err + 0x08;    s1p34r = sse2_rnd + 0x08;
		s1p05r = sc_ptr + 0x0a;    s1p15r = two + 0x0a;    s1p25r = max_err + 0x0a;    s1p35r = sse2_rnd + 0x0a;
		s1p06r = sc_ptr + 0x0c;    s1p16r = two + 0x0c;    s1p26r = max_err + 0x0c;    s1p36r = sse2_rnd + 0x0c;
		s1p07r = sc_ptr + 0x0e;    s1p17r = two + 0x0e;    s1p27r = max_err + 0x0e;    s1p37r = sse2_rnd + 0x0e;
		s1p08r = sc_ptr + 0x10;    s1p18r = two + 0x10;    s1p28r = max_err + 0x10;    s1p38r = sse2_rnd + 0x10;
		s1p09r = sc_ptr + 0x12;    s1p19r = two + 0x12;    s1p29r = max_err + 0x12;    s1p39r = sse2_rnd + 0x12;
		s1p0ar = sc_ptr + 0x14;    s1p1ar = two + 0x14;    s1p2ar = max_err + 0x14;    s1p3ar = sse2_rnd + 0x14;
		s1p0br = sc_ptr + 0x16;    s1p1br = two + 0x16;    s1p2br = max_err + 0x16;    s1p3br = sse2_rnd + 0x16;
		s1p0cr = sc_ptr + 0x18;    s1p1cr = two + 0x18;    s1p2cr = max_err + 0x18;    s1p3cr = sse2_rnd + 0x18;
	// sc_ptr += 104
		tmp = sse2_rnd + 0x1a; two = tmp + 0x1a;    max_err = two + 0x1a;   sse2_rnd = max_err + 0x1a;	/* Use these as handy temps */
		r00r = tmp + 0x00;    r10r = two + 0x00;    r20r = max_err + 0x00;    r30r = sse2_rnd + 0x00;
		r01r = tmp + 0x02;    r11r = two + 0x02;    r21r = max_err + 0x02;    r31r = sse2_rnd + 0x02;
		r02r = tmp + 0x04;    r12r = two + 0x04;    r22r = max_err + 0x04;    r32r = sse2_rnd + 0x04;
		r03r = tmp + 0x06;    r13r = two + 0x06;    r23r = max_err + 0x06;    r33r = sse2_rnd + 0x06;
		r04r = tmp + 0x08;    r14r = two + 0x08;    r24r = max_err + 0x08;    r34r = sse2_rnd + 0x08;
		r05r = tmp + 0x0a;    r15r = two + 0x0a;    r25r = max_err + 0x0a;    r35r = sse2_rnd + 0x0a;
		r06r = tmp + 0x0c;    r16r = two + 0x0c;    r26r = max_err + 0x0c;    r36r = sse2_rnd + 0x0c;
		r07r = tmp + 0x0e;    r17r = two + 0x0e;    r27r = max_err + 0x0e;    r37r = sse2_rnd + 0x0e;
		r08r = tmp + 0x10;    r18r = two + 0x10;    r28r = max_err + 0x10;    r38r = sse2_rnd + 0x10;
		r09r = tmp + 0x12;    r19r = two + 0x12;    r29r = max_err + 0x12;    r39r = sse2_rnd + 0x12;
		r0ar = tmp + 0x14;    r1ar = two + 0x14;    r2ar = max_err + 0x14;    r3ar = sse2_rnd + 0x14;
		r0br = tmp + 0x16;    r1br = two + 0x16;    r2br = max_err + 0x16;    r3br = sse2_rnd + 0x16;
		r0cr = tmp + 0x18;    r1cr = two + 0x18;    r2cr = max_err + 0x18;    r3cr = sse2_rnd + 0x18;
	// sc_ptr += 208
		tmp = sse2_rnd + 0x1a;
		two = tmp;
		rad13_const = tmp + 0x01;	/* Leave an extra slot at radix13_const-1 for the constant two = 2.0: */
		tmp += 0x14;	/* Need 20 16-byte slots for two+sincos, but offset the carry slots by the next-larger multiple of 4 */
	// sc_ptr += 228
	#ifdef USE_AVX
		cy00 = tmp + 0x00;
		cy04 = tmp + 0x01;
		cy08 = tmp + 0x02;
		cy12 = tmp + 0x03;
		cy16 = tmp + 0x04;
		cy20 = tmp + 0x05;
		cy24 = tmp + 0x06;
		cy28 = tmp + 0x07;
		cy32 = tmp + 0x08;
		cy36 = tmp + 0x09;
		cy40 = tmp + 0x0a;
		cy44 = tmp + 0x0b;
		cy48 = tmp + 0x0c;
		tmp += 0x0d;
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
	// sc_ptr += 243; This is where the value of half_arr_offset52 comes from
		half_arr= tmp + 0x02;	/* This table needs 20x16 bytes for Mersenne-mod, and radixx16 for Fermat-mod */
	#else
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
	// sc_ptr += 256; This is where the value of half_arr_offset52 comes from
		half_arr= tmp + 0x02;	/* This table needs 20x16 bytes for Mersenne-mod, and radixx16 for Fermat-mod */
	#endif
		ASSERT(HERE, (radix52_creals_in_local_store << l2_sz_vd) >= ((long)half_arr - (long)s1p00r) + (20 << l2_sz_vd), "radix52_creals_in_local_store checksum failed!");
		/* These remain fixed: */
		tmp = two;		/* __cc pointer offsets: */
		VEC_DBL_INIT(tmp,  2.0);	++tmp;	/*	-0x010 = 2.0 */
		VEC_DBL_INIT(tmp,  DC1);	++tmp;	/*	0x000 =  DC1 */
		VEC_DBL_INIT(tmp,  DC3);	++tmp;	/*	0x010 =  DC3 */
		VEC_DBL_INIT(tmp,  DC4);	++tmp;	/*	0x020 =  DC4 */
		VEC_DBL_INIT(tmp,  DS1);	++tmp;	/*	0x030 =  DS1 */
		VEC_DBL_INIT(tmp,  DS2);	++tmp;	/*	0x040 =  DS2 */
		VEC_DBL_INIT(tmp,  DS3);	++tmp;	/*	0x050 =  DS3 */
		VEC_DBL_INIT(tmp,  DS4);	++tmp;	/*	0x060 =  DS4 */
		VEC_DBL_INIT(tmp,  DS5);	++tmp;	/*	0x070 =  DS5 */
		VEC_DBL_INIT(tmp, DC23);	++tmp;	/*	0x080 = DC23 */
		VEC_DBL_INIT(tmp, DC54);	++tmp;	/*	0x090 = DC54 */
		VEC_DBL_INIT(tmp, DC65);	++tmp;	/*	0x0a0 = DC65 */
		VEC_DBL_INIT(tmp, DS63);	++tmp;	/*	0x0b0 = DS63 */
		VEC_DBL_INIT(tmp, DS74);	++tmp;	/*	0x0c0 = DS74 */
		VEC_DBL_INIT(tmp, DS85);	++tmp;	/*	0x0d0 = DS85 */
		VEC_DBL_INIT(tmp, DS93);	++tmp;	/*	0x0e0 = DS93 */
		VEC_DBL_INIT(tmp, DSa4);	++tmp;	/*	0x0f0 = DSa4 */
		VEC_DBL_INIT(tmp, DSb5);	++tmp;	/*	0x100 = DSb5 */
		VEC_DBL_INIT(sse2_rnd,crnd);

		// Propagate the above consts to the remaining threads:
		nbytes = (int)tmp - (int)two;	// #bytes in above sincos block of data
		tmp = two;
		tm2 = tmp + cslots_in_local_store;
		for(ithread = 1; ithread < CY_THREADS; ++ithread) {
			memcpy(tm2, tmp, nbytes);
			tmp = tm2;		tm2 += cslots_in_local_store;
		}
		nbytes = sz_vd;	// sse2_rnd is a solo (in the SIMD-vector) datum
		tmp = sse2_rnd;
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
		*/
		tmp = half_arr;

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

	#endif	// USE_SSE2

	#ifdef USE_PTHREAD
		/* Populate the elements of the thread-specific data structs which don't change after init: */
		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
		#ifdef USE_SSE2
			tdat[ithread].s1p00r = __r0 + ithread*cslots_in_local_store;
			tdat[ithread].half_arr = (long)tdat[ithread].s1p00r + ((long)half_arr - (long)s1p00r);
		#else	// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
			tdat[ithread].s1p00r      = (vec_dbl *)base;
			tdat[ithread].half_arr = (vec_dbl *)baseinv;
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
	#ifdef USE_SSE2
		ASSERT(HERE, tdat[ithread].s1p00r == __r0 + ithread*cslots_in_local_store, "thread-local memcheck fail!");
		tmp = tdat[ithread].half_arr;
		ASSERT(HERE, ((tmp-1)->d0 == crnd && (tmp-1)->d1 == crnd), "thread-local memcheck failed!");
	  #ifdef USE_AVX
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif
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
	//	VEC_DBL_INIT(max_err, 0.0);	*** must do this in conjunction with thread-local-data-copy
	#endif

		i      = _i[ithread];	/* Pointer to the BASE and BASEINV arrays.	*/
		jstart = _jstart[ithread];
		jhi    = _jhi[ithread];

		col = _col[ithread];
		co2 = _co2[ithread];
		co3 = _co3[ithread];

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
	#endif
		/* init carries	*/
	#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
		cy00->d0 = _cy00[ithread];	cy00->d1 = _cy01[ithread];	cy00->d2 = _cy02[ithread];	cy00->d3 = _cy03[ithread];
		cy04->d0 = _cy04[ithread];	cy04->d1 = _cy05[ithread];	cy04->d2 = _cy06[ithread];	cy04->d3 = _cy07[ithread];
		cy08->d0 = _cy08[ithread];	cy08->d1 = _cy09[ithread];	cy08->d2 = _cy10[ithread];	cy08->d3 = _cy11[ithread];
		cy12->d0 = _cy12[ithread];	cy12->d1 = _cy13[ithread];	cy12->d2 = _cy14[ithread];	cy12->d3 = _cy15[ithread];
		cy16->d0 = _cy16[ithread];	cy16->d1 = _cy17[ithread];	cy16->d2 = _cy18[ithread];	cy16->d3 = _cy19[ithread];
		cy20->d0 = _cy20[ithread];	cy20->d1 = _cy21[ithread];	cy20->d2 = _cy22[ithread];	cy20->d3 = _cy23[ithread];
		cy24->d0 = _cy24[ithread];	cy24->d1 = _cy25[ithread];	cy24->d2 = _cy26[ithread];	cy24->d3 = _cy27[ithread];
		cy28->d0 = _cy28[ithread];	cy28->d1 = _cy29[ithread];	cy28->d2 = _cy30[ithread];	cy28->d3 = _cy31[ithread];
		cy32->d0 = _cy32[ithread];	cy32->d1 = _cy33[ithread];	cy32->d2 = _cy34[ithread];	cy32->d3 = _cy35[ithread];
		cy36->d0 = _cy36[ithread];	cy36->d1 = _cy37[ithread];	cy36->d2 = _cy38[ithread];	cy36->d3 = _cy39[ithread];
		cy40->d0 = _cy40[ithread];	cy40->d1 = _cy41[ithread];	cy40->d2 = _cy42[ithread];	cy40->d3 = _cy43[ithread];
		cy44->d0 = _cy44[ithread];	cy44->d1 = _cy45[ithread];	cy44->d2 = _cy46[ithread];	cy44->d3 = _cy47[ithread];
		cy48->d0 = _cy48[ithread];	cy48->d1 = _cy49[ithread];	cy48->d2 = _cy50[ithread];	cy48->d3 = _cy51[ithread];
	#elif defined(USE_SSE2)
		cy00->d0 = _cy00[ithread];	cy00->d1 = _cy01[ithread];
		cy02->d0 = _cy02[ithread];	cy02->d1 = _cy03[ithread];
		cy04->d0 = _cy04[ithread];	cy04->d1 = _cy05[ithread];
		cy06->d0 = _cy06[ithread];	cy06->d1 = _cy07[ithread];
		cy08->d0 = _cy08[ithread];	cy08->d1 = _cy09[ithread];
		cy10->d0 = _cy10[ithread];	cy10->d1 = _cy11[ithread];
		cy12->d0 = _cy12[ithread];	cy12->d1 = _cy13[ithread];
		cy14->d0 = _cy14[ithread];	cy14->d1 = _cy15[ithread];
		cy16->d0 = _cy16[ithread];	cy16->d1 = _cy17[ithread];
		cy18->d0 = _cy18[ithread];	cy18->d1 = _cy19[ithread];
		cy20->d0 = _cy20[ithread];	cy20->d1 = _cy21[ithread];
		cy22->d0 = _cy22[ithread];	cy22->d1 = _cy23[ithread];
		cy24->d0 = _cy24[ithread];	cy24->d1 = _cy25[ithread];
		cy26->d0 = _cy26[ithread];	cy26->d1 = _cy27[ithread];
		cy28->d0 = _cy28[ithread];	cy28->d1 = _cy29[ithread];
		cy30->d0 = _cy30[ithread];	cy30->d1 = _cy31[ithread];
		cy32->d0 = _cy32[ithread];	cy32->d1 = _cy33[ithread];
		cy34->d0 = _cy34[ithread];	cy34->d1 = _cy35[ithread];
		cy36->d0 = _cy36[ithread];	cy36->d1 = _cy37[ithread];
		cy38->d0 = _cy38[ithread];	cy38->d1 = _cy39[ithread];
		cy40->d0 = _cy40[ithread];	cy40->d1 = _cy41[ithread];
		cy42->d0 = _cy42[ithread];	cy42->d1 = _cy43[ithread];
		cy44->d0 = _cy44[ithread];	cy44->d1 = _cy45[ithread];
		cy46->d0 = _cy46[ithread];	cy46->d1 = _cy47[ithread];
		cy48->d0 = _cy48[ithread];	cy48->d1 = _cy49[ithread];
		cy50->d0 = _cy50[ithread];	cy50->d1 = _cy51[ithread];
	#else
		cy00 = _cy00[ithread];
		cy01 = _cy01[ithread];
		cy02 = _cy02[ithread];
		cy03 = _cy03[ithread];
		cy04 = _cy04[ithread];
		cy05 = _cy05[ithread];
		cy06 = _cy06[ithread];
		cy07 = _cy07[ithread];
		cy08 = _cy08[ithread];
		cy09 = _cy09[ithread];
		cy10 = _cy10[ithread];
		cy11 = _cy11[ithread];
		cy12 = _cy12[ithread];
		cy13 = _cy13[ithread];
		cy14 = _cy14[ithread];
		cy15 = _cy15[ithread];
		cy16 = _cy16[ithread];
		cy17 = _cy17[ithread];
		cy18 = _cy18[ithread];
		cy19 = _cy19[ithread];
		cy20 = _cy20[ithread];
		cy21 = _cy21[ithread];
		cy22 = _cy22[ithread];
		cy23 = _cy23[ithread];
		cy24 = _cy24[ithread];
		cy25 = _cy25[ithread];
		cy26 = _cy26[ithread];
		cy27 = _cy27[ithread];
		cy28 = _cy28[ithread];
		cy29 = _cy29[ithread];
		cy30 = _cy30[ithread];
		cy31 = _cy31[ithread];
		cy32 = _cy32[ithread];
		cy33 = _cy33[ithread];
		cy34 = _cy34[ithread];
		cy35 = _cy35[ithread];
		cy36 = _cy36[ithread];
		cy37 = _cy37[ithread];
		cy38 = _cy38[ithread];
		cy39 = _cy39[ithread];
		cy40 = _cy40[ithread];
		cy41 = _cy41[ithread];
		cy42 = _cy42[ithread];
		cy43 = _cy43[ithread];
		cy44 = _cy44[ithread];
		cy45 = _cy45[ithread];
		cy46 = _cy46[ithread];
		cy47 = _cy47[ithread];
		cy48 = _cy48[ithread];
		cy49 = _cy49[ithread];
		cy50 = _cy50[ithread];
		cy51 = _cy51[ithread];
	#endif

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

		#ifdef USE_SSE2

		  #ifdef USE_AVX
			/* Outputs in AVX mode are temps 2*13*32 = 26*32 = 0x340 bytes apart: */
			// Reorder blocks to yield sequentially increasing a-array offsets:
/* Block 1 : */	add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00r, 0x340)
/* Block 4 : */	add3 = &a[j1+p04];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03r, 0x340)
/* Block 7 : */	add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06r, 0x340)
/* Block 10: */	add2 = &a[j1+p12];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09r, 0x340)
/* Block 13: */	add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0cr, 0x340)
/* Block 3 : */	add3 = &a[j1+p20];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02r, 0x340)
/* Block 6 : */	add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05r, 0x340)
/* Block 9 : */	add2 = &a[j1+p28];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08r, 0x340)
/* Block 12: */	add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0br, 0x340)
/* Block 2 : */	add3 = &a[j1+p36];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01r, 0x340)
/* Block 5 : */	add1 = &a[j1+p40];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04r, 0x340)
/* Block 8 : */	add2 = &a[j1+p44];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07r, 0x340)
/* Block 11: */	add0 = &a[j1+p48];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0ar, 0x340)

			/* Radix-13 DFT uses adjacent temps for inputs, outputs are (cyclic) with pXXr having XX += 40 (40*32 bytes = 0x500) or XX -= 12 (-12*32 bytes = -0x180) between successive outputs: */
																													/*  a1p00r,a1p31r,a1p22r,a1p13r,a1p04r,a1p35r,a1p26r,a1p17r,a1p08r,a1p39r,a1p2ar,a1p1br,a1p0cr */
			SSE2_RADIX_13_DFT(rad13_const, r00r,0x040,0x080,0x0c0,0x100,0x140,0x180,0x1c0,0x200,0x240,0x280,0x2c0,0x300, s1p00r, 0xa00, 0x700, 0x400, 0x100, 0xb00, 0x800, 0x500, 0x200, 0xc00, 0x900, 0x600, 0x300)
																													/*  a1p30r,a1p21r,a1p12r,a1p03r,a1p34r,a1p25r,a1p16r,a1p07r,a1p38r,a1p29r,a1p1ar,a1p0br,a1p3cr */
			SSE2_RADIX_13_DFT(rad13_const, r10r,0x040,0x080,0x0c0,0x100,0x140,0x180,0x1c0,0x200,0x240,0x280,0x2c0,0x300, s1p30r,-0x300,-0x600,-0x900, 0x100,-0x200,-0x500,-0x800, 0x200,-0x100,-0x400,-0x700, 0x300)
																													/*  a1p20r,a1p11r,a1p02r,a1p33r,a1p24r,a1p15r,a1p06r,a1p37r,a1p28r,a1p19r,a1p0ar,a1p3br,a1p2cr */
			SSE2_RADIX_13_DFT(rad13_const, r20r,0x040,0x080,0x0c0,0x100,0x140,0x180,0x1c0,0x200,0x240,0x280,0x2c0,0x300, s1p20r,-0x300,-0x600, 0x400, 0x100,-0x200,-0x500, 0x500, 0x200,-0x100,-0x400, 0x600, 0x300)
																													/*  a1p10r,a1p01r,a1p32r,a1p23r,a1p14r,a1p05r,a1p36r,a1p27r,a1p18r,a1p09r,a1p3ar,a1p2br,a1p1cr */
			SSE2_RADIX_13_DFT(rad13_const, r30r,0x040,0x080,0x0c0,0x100,0x140,0x180,0x1c0,0x200,0x240,0x280,0x2c0,0x300, s1p10r,-0x300, 0x700, 0x400, 0x100,-0x200, 0x800, 0x500, 0x200,-0x100, 0x900, 0x600, 0x300)

		  #else

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

		  #endif	// AVX or SSE2?

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

			AVX_cmplx_carry_norm_errcheck0_X4(s1p00r,add1,add2,add3,cy00,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			AVX_cmplx_carry_norm_errcheck1_X4(s1p04r,add1,add2,add3,cy04,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			AVX_cmplx_carry_norm_errcheck1_X4(s1p08r,add1,add2,add3,cy08,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			AVX_cmplx_carry_norm_errcheck1_X4(s1p0cr,add1,add2,add3,cy12,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			AVX_cmplx_carry_norm_errcheck1_X4(s1p13r,add1,add2,add3,cy16,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			AVX_cmplx_carry_norm_errcheck1_X4(s1p17r,add1,add2,add3,cy20,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			AVX_cmplx_carry_norm_errcheck1_X4(s1p1br,add1,add2,add3,cy24,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			AVX_cmplx_carry_norm_errcheck1_X4(s1p22r,add1,add2,add3,cy28,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			AVX_cmplx_carry_norm_errcheck1_X4(s1p26r,add1,add2,add3,cy32,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			AVX_cmplx_carry_norm_errcheck1_X4(s1p2ar,add1,add2,add3,cy36,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			AVX_cmplx_carry_norm_errcheck1_X4(s1p31r,add1,add2,add3,cy40,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			AVX_cmplx_carry_norm_errcheck1_X4(s1p35r,add1,add2,add3,cy44,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
			AVX_cmplx_carry_norm_errcheck1_X4(s1p39r,add1,add2,add3,cy48,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);

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

		/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/

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

			i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
				 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		#endif	// USE_AVX?

		/*...The radix-52 DIF pass is here:	*/

		#ifdef USE_SSE2

		  #ifdef USE_AVX
			/*...gather the needed data (52 64-bit complex, i.e. 104 64-bit reals) and do 4 radix-13 transforms...*/
			/* Radix-13 DFT inputs are (cyclic) with pXXr having XX += 48 (48*32 bytes = +0x600) or XX -= 4 (4*32 bytes = -0x080), outputs are adjacent 32-byte-separated temps: */
			SSE2_RADIX_13_DFT(rad13_const, s1p00r, 0xc00, 0xb00, 0xa00, 0x900, 0x800, 0x700, 0x600, 0x500, 0x400, 0x300, 0x200, 0x100, r00r,0x040,0x080,0x0c0,0x100,0x140,0x180,0x1c0,0x200,0x240,0x280,0x2c0,0x300)
			SSE2_RADIX_13_DFT(rad13_const, s1p30r,-0x100,-0x200,-0x300,-0x400,-0x500,-0x600,-0x700,-0x800,-0x900, 0x300, 0x200, 0x100, r10r,0x040,0x080,0x0c0,0x100,0x140,0x180,0x1c0,0x200,0x240,0x280,0x2c0,0x300)
			SSE2_RADIX_13_DFT(rad13_const, s1p20r,-0x100,-0x200,-0x300,-0x400,-0x500,-0x600, 0x600, 0x500, 0x400, 0x300, 0x200, 0x100, r20r,0x040,0x080,0x0c0,0x100,0x140,0x180,0x1c0,0x200,0x240,0x280,0x2c0,0x300)
			SSE2_RADIX_13_DFT(rad13_const, s1p10r,-0x100,-0x200,-0x300, 0x900, 0x800, 0x700, 0x600, 0x500, 0x400, 0x300, 0x200, 0x100, r30r,0x040,0x080,0x0c0,0x100,0x140,0x180,0x1c0,0x200,0x240,0x280,0x2c0,0x300)

			/*...and now do 13 radix-4 transforms...*/
			/* Inputs in AVX mode are temps 2*13*32 = 26*32 = 0x340 bytes apart. Notice how the add* indices after the first row repeat with period 4: */
			// Reorder blocks to yield sequentially increasing a-array offsets:
/* Block 01 : */	add0 = &a[j1    ];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00r, 0x340)
/* Block 13 : */	add3 = &a[j1+p04];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0cr, 0x340)
/* Block 12 : */	add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0br, 0x340)
/* Block 11 : */	add2 = &a[j1+p12];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0ar, 0x340)
/* Block 10 : */	add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09r, 0x340)
/* Block 09 : */	add3 = &a[j1+p20];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08r, 0x340)
/* Block 08 : */	add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07r, 0x340)
/* Block 07 : */	add2 = &a[j1+p28];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06r, 0x340)
/* Block 06 : */	add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05r, 0x340)
/* Block 05 : */	add3 = &a[j1+p36];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04r, 0x340)
/* Block 04 : */	add1 = &a[j1+p40];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03r, 0x340)
/* Block 03 : */	add2 = &a[j1+p44];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02r, 0x340)
/* Block 02 : */	add0 = &a[j1+p48];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01r, 0x340)

		  #else

			/*...gather the needed data (52 64-bit complex, i.e. 104 64-bit reals) and do 4 radix-13 transforms...*/
			/* Radix-13 DFT inputs are (cyclic) with pXXr having XX += 48 (48*32 bytes = +0x600) or XX -= 4 (4*32 bytes = -0x080), outputs are adjacent 32-byte-separated temps: */
			SSE2_RADIX_13_DFT(rad13_const, s1p00r, 0x600, 0x580, 0x500, 0x480, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, r00r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180)
			SSE2_RADIX_13_DFT(rad13_const, s1p30r,-0x080,-0x100,-0x180,-0x200,-0x280,-0x300,-0x380,-0x400,-0x480, 0x180, 0x100, 0x080, r10r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180)
			SSE2_RADIX_13_DFT(rad13_const, s1p20r,-0x080,-0x100,-0x180,-0x200,-0x280,-0x300, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, r20r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180)
			SSE2_RADIX_13_DFT(rad13_const, s1p10r,-0x080,-0x100,-0x180, 0x480, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, r30r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180)

			/*...and now do 13 radix-4 transforms...*/
			/* Inputs in SSE2 mode are temps 2*13*16 = 26*16 = 0x1a0 bytes apart. Notice how the add* indices after the first row repeat with period 4: */
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

		  #endif	// AVX or SSE2?

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
	#ifdef USE_AVX	// AVX and AVX2 both use 256-bit registers
		_cy00[ithread] = cy00->d0;	_cy01[ithread] = cy00->d1;	_cy02[ithread] = cy00->d2;	_cy03[ithread] = cy00->d3;
		_cy04[ithread] = cy04->d0;	_cy05[ithread] = cy04->d1;	_cy06[ithread] = cy04->d2;	_cy07[ithread] = cy04->d3;
		_cy08[ithread] = cy08->d0;	_cy09[ithread] = cy08->d1;	_cy10[ithread] = cy08->d2;	_cy11[ithread] = cy08->d3;
		_cy12[ithread] = cy12->d0;	_cy13[ithread] = cy12->d1;	_cy14[ithread] = cy12->d2;	_cy15[ithread] = cy12->d3;
		_cy16[ithread] = cy16->d0;	_cy17[ithread] = cy16->d1;	_cy18[ithread] = cy16->d2;	_cy19[ithread] = cy16->d3;
		_cy20[ithread] = cy20->d0;	_cy21[ithread] = cy20->d1;	_cy22[ithread] = cy20->d2;	_cy23[ithread] = cy20->d3;
		_cy24[ithread] = cy24->d0;	_cy25[ithread] = cy24->d1;	_cy26[ithread] = cy24->d2;	_cy27[ithread] = cy24->d3;
		_cy28[ithread] = cy28->d0;	_cy29[ithread] = cy28->d1;	_cy30[ithread] = cy28->d2;	_cy31[ithread] = cy28->d3;
		_cy32[ithread] = cy32->d0;	_cy33[ithread] = cy32->d1;	_cy34[ithread] = cy32->d2;	_cy35[ithread] = cy32->d3;
		_cy36[ithread] = cy36->d0;	_cy37[ithread] = cy36->d1;	_cy38[ithread] = cy36->d2;	_cy39[ithread] = cy36->d3;
		_cy40[ithread] = cy40->d0;	_cy41[ithread] = cy40->d1;	_cy42[ithread] = cy40->d2;	_cy43[ithread] = cy40->d3;
		_cy44[ithread] = cy44->d0;	_cy45[ithread] = cy44->d1;	_cy46[ithread] = cy44->d2;	_cy47[ithread] = cy44->d3;
		_cy48[ithread] = cy48->d0;	_cy49[ithread] = cy48->d1;	_cy50[ithread] = cy48->d2;	_cy51[ithread] = cy48->d3;
		maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );
	#elif defined(USE_SSE2)
		_cy00[ithread] = cy00->d0;	_cy01[ithread] = cy00->d1;
		_cy02[ithread] = cy02->d0;	_cy03[ithread] = cy02->d1;
		_cy04[ithread] = cy04->d0;	_cy05[ithread] = cy04->d1;
		_cy06[ithread] = cy06->d0;	_cy07[ithread] = cy06->d1;
		_cy08[ithread] = cy08->d0;	_cy09[ithread] = cy08->d1;
		_cy10[ithread] = cy10->d0;	_cy11[ithread] = cy10->d1;
		_cy12[ithread] = cy12->d0;	_cy13[ithread] = cy12->d1;
		_cy14[ithread] = cy14->d0;	_cy15[ithread] = cy14->d1;
		_cy16[ithread] = cy16->d0;	_cy17[ithread] = cy16->d1;
		_cy18[ithread] = cy18->d0;	_cy19[ithread] = cy18->d1;
		_cy20[ithread] = cy20->d0;	_cy21[ithread] = cy20->d1;
		_cy22[ithread] = cy22->d0;	_cy23[ithread] = cy22->d1;
		_cy24[ithread] = cy24->d0;	_cy25[ithread] = cy24->d1;
		_cy26[ithread] = cy26->d0;	_cy27[ithread] = cy26->d1;
		_cy28[ithread] = cy28->d0;	_cy29[ithread] = cy28->d1;
		_cy30[ithread] = cy30->d0;	_cy31[ithread] = cy30->d1;
		_cy32[ithread] = cy32->d0;	_cy33[ithread] = cy32->d1;
		_cy34[ithread] = cy34->d0;	_cy35[ithread] = cy34->d1;
		_cy36[ithread] = cy36->d0;	_cy37[ithread] = cy36->d1;
		_cy38[ithread] = cy38->d0;	_cy39[ithread] = cy38->d1;
		_cy40[ithread] = cy40->d0;	_cy41[ithread] = cy40->d1;
		_cy42[ithread] = cy42->d0;	_cy43[ithread] = cy42->d1;
		_cy44[ithread] = cy44->d0;	_cy45[ithread] = cy44->d1;
		_cy46[ithread] = cy46->d0;	_cy47[ithread] = cy46->d1;
		_cy48[ithread] = cy48->d0;	_cy49[ithread] = cy48->d1;
		_cy50[ithread] = cy50->d0;	_cy51[ithread] = cy50->d1;
		maxerr = MAX(max_err->d0,max_err->d1);
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
	static int NDIVR,first_entry=TRUE,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40,p44,p48;
	double rt,it,
		t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,
		t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,
		t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,
		t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci;

	if(!first_entry && (n/52) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR=n/52;

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
	}

/*...The radix-52 pass is here.	*/

	for(j=0; j < NDIVR; j += 2)
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
	int j,j1,j2;
	static int NDIVR,first_entry=TRUE,p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40,p44,p48;
	int jt,jp;
	double rt,it,
		t00r,t00i,t01r,t01i,t02r,t02i,t03r,t03i,t04r,t04i,t05r,t05i,t06r,t06i,t07r,t07i,t08r,t08i,t09r,t09i,t0ar,t0ai,t0br,t0bi,t0cr,t0ci,
		t10r,t10i,t11r,t11i,t12r,t12i,t13r,t13i,t14r,t14i,t15r,t15i,t16r,t16i,t17r,t17i,t18r,t18i,t19r,t19i,t1ar,t1ai,t1br,t1bi,t1cr,t1ci,
		t20r,t20i,t21r,t21i,t22r,t22i,t23r,t23i,t24r,t24i,t25r,t25i,t26r,t26i,t27r,t27i,t28r,t28i,t29r,t29i,t2ar,t2ai,t2br,t2bi,t2cr,t2ci,
		t30r,t30i,t31r,t31i,t32r,t32i,t33r,t33i,t34r,t34i,t35r,t35i,t36r,t36i,t37r,t37i,t38r,t38i,t39r,t39i,t3ar,t3ai,t3br,t3bi,t3cr,t3ci;

	if(!first_entry && (n/52) != NDIVR)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		NDIVR=n/52;

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
	}

/*...The radix-52 pass is here.	*/

	for(j=0; j < NDIVR; j += 2)
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
	}
}

/******************** Multithreaded function body: ***************************/

#ifdef USE_PTHREAD

	#ifndef COMPILER_TYPE_GCC
		#error pthreaded carry code requires GCC build!
	#endif

	void* 
	cy52_process_chunk(void*targ)	// Thread-arg pointer *must* be cast to void and specialized inside the function
	{
		const uint32 RADIX = 52;
		const int stride = (int)RE_IM_STRIDE << 1;	// main-array loop stride = 2*RE_IM_STRIDE
		uint32 p01,p02,p03,p04,p08,p12,p16,p20,p24,p28,p32,p36,p40,p44,p48;
		int j,j1,j2,k,l;
		double wtl,wtlp1,wtn,wtnm1;	/* Mersenne-mod weights stuff */
	#ifdef USE_AVX
		struct uint32x4 *n_minus_sil,*n_minus_silp1,*sinwt,*sinwtm1;
	#else
		int n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	#endif

	#ifdef USE_SSE2

		double *add0, *add1, *add2, *add3;
		const double crnd = 3.0*0x4000000*0x2000000;
		struct complex *ctmp;	// Hybrid AVX-DFT/SSE2-carry scheme used for Mersenne-mod needs a 2-word-double pointer
		int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09
			,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19
			,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27,*bjmodn28,*bjmodn29
			,*bjmodn30,*bjmodn31,*bjmodn32,*bjmodn33,*bjmodn34,*bjmodn35,*bjmodn36,*bjmodn37,*bjmodn38,*bjmodn39
			,*bjmodn40,*bjmodn41,*bjmodn42,*bjmodn43,*bjmodn44,*bjmodn45,*bjmodn46,*bjmodn47,*bjmodn48,*bjmodn49
			,*bjmodn50,*bjmodn51;
		vec_dbl *two,*rad13_const, *max_err, *sse2_rnd, *half_arr, *tmp     /* rad13_const needs 18*16 bytes allocated */
			,*r00r,*r01r,*r02r,*r03r,*r04r,*r05r,*r06r,*r07r,*r08r,*r09r,*r0ar,*r0br,*r0cr
			,*r10r,*r11r,*r12r,*r13r,*r14r,*r15r,*r16r,*r17r,*r18r,*r19r,*r1ar,*r1br,*r1cr
			,*r20r,*r21r,*r22r,*r23r,*r24r,*r25r,*r26r,*r27r,*r28r,*r29r,*r2ar,*r2br,*r2cr
			,*r30r,*r31r,*r32r,*r33r,*r34r,*r35r,*r36r,*r37r,*r38r,*r39r,*r3ar,*r3br,*r3cr
			,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p0ar,*s1p0br,*s1p0cr
			,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p1ar,*s1p1br,*s1p1cr
			,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r,*s1p28r,*s1p29r,*s1p2ar,*s1p2br,*s1p2cr
			,*s1p30r,*s1p31r,*s1p32r,*s1p33r,*s1p34r,*s1p35r,*s1p36r,*s1p37r,*s1p38r,*s1p39r,*s1p3ar,*s1p3br,*s1p3cr;
		vec_dbl
			*cy00,*cy04,*cy08,*cy12,*cy16,*cy20,*cy24,*cy28,*cy32,*cy36,*cy40,*cy44,*cy48;
	  #ifndef USE_AVX
		vec_dbl
			*cy02,*cy06,*cy10,*cy14,*cy18,*cy22,*cy26,*cy30,*cy34,*cy38,*cy42,*cy46,*cy50;
	  #endif
		double dtmp;
		uint64 *sign_mask, *sse_bw, *sse_sw, *sse_n;

	#else

		double *base, *baseinv;
		const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
		int jt,jp,m,m2;
		double wt,wtinv,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
		double rt,it,temp,frac,
			t00r,t01r,t02r,t03r,t04r,t05r,t06r,t07r,t08r,t09r,t0ar,t0br,t0cr,
			t10r,t11r,t12r,t13r,t14r,t15r,t16r,t17r,t18r,t19r,t1ar,t1br,t1cr,
			t20r,t21r,t22r,t23r,t24r,t25r,t26r,t27r,t28r,t29r,t2ar,t2br,t2cr,
			t30r,t31r,t32r,t33r,t34r,t35r,t36r,t37r,t38r,t39r,t3ar,t3br,t3cr,
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
		int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,
			bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,
			bjmodn26,bjmodn27,bjmodn28,bjmodn29,bjmodn30,bjmodn31,bjmodn32,bjmodn33,bjmodn34,bjmodn35,bjmodn36,bjmodn37,bjmodn38,
			bjmodn39,bjmodn40,bjmodn41,bjmodn42,bjmodn43,bjmodn44,bjmodn45,bjmodn46,bjmodn47,bjmodn48,bjmodn49,bjmodn50,bjmodn51;

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

	#ifdef USE_SSE2
		tmp	= thread_arg->s1p00r;
								    two = tmp + 0x1a;	max_err  = two + 0x1a;		sse2_rnd = max_err + 0x1a;	/* Use these as handy temps */
		s1p00r = tmp + 0x00;    s1p10r = two + 0x00;    s1p20r = max_err + 0x00;    s1p30r = sse2_rnd + 0x00;
		s1p01r = tmp + 0x02;    s1p11r = two + 0x02;    s1p21r = max_err + 0x02;    s1p31r = sse2_rnd + 0x02;
		s1p02r = tmp + 0x04;    s1p12r = two + 0x04;    s1p22r = max_err + 0x04;    s1p32r = sse2_rnd + 0x04;
		s1p03r = tmp + 0x06;    s1p13r = two + 0x06;    s1p23r = max_err + 0x06;    s1p33r = sse2_rnd + 0x06;
		s1p04r = tmp + 0x08;    s1p14r = two + 0x08;    s1p24r = max_err + 0x08;    s1p34r = sse2_rnd + 0x08;
		s1p05r = tmp + 0x0a;    s1p15r = two + 0x0a;    s1p25r = max_err + 0x0a;    s1p35r = sse2_rnd + 0x0a;
		s1p06r = tmp + 0x0c;    s1p16r = two + 0x0c;    s1p26r = max_err + 0x0c;    s1p36r = sse2_rnd + 0x0c;
		s1p07r = tmp + 0x0e;    s1p17r = two + 0x0e;    s1p27r = max_err + 0x0e;    s1p37r = sse2_rnd + 0x0e;
		s1p08r = tmp + 0x10;    s1p18r = two + 0x10;    s1p28r = max_err + 0x10;    s1p38r = sse2_rnd + 0x10;
		s1p09r = tmp + 0x12;    s1p19r = two + 0x12;    s1p29r = max_err + 0x12;    s1p39r = sse2_rnd + 0x12;
		s1p0ar = tmp + 0x14;    s1p1ar = two + 0x14;    s1p2ar = max_err + 0x14;    s1p3ar = sse2_rnd + 0x14;
		s1p0br = tmp + 0x16;    s1p1br = two + 0x16;    s1p2br = max_err + 0x16;    s1p3br = sse2_rnd + 0x16;
		s1p0cr = tmp + 0x18;    s1p1cr = two + 0x18;    s1p2cr = max_err + 0x18;    s1p3cr = sse2_rnd + 0x18;

		tmp = sse2_rnd + 0x1a; two = tmp + 0x1a;    max_err = two + 0x1a;   sse2_rnd = max_err + 0x1a;	/* Use these as handy temps */
		r00r = tmp + 0x00;    r10r = two + 0x00;    r20r = max_err + 0x00;    r30r = sse2_rnd + 0x00;
		r01r = tmp + 0x02;    r11r = two + 0x02;    r21r = max_err + 0x02;    r31r = sse2_rnd + 0x02;
		r02r = tmp + 0x04;    r12r = two + 0x04;    r22r = max_err + 0x04;    r32r = sse2_rnd + 0x04;
		r03r = tmp + 0x06;    r13r = two + 0x06;    r23r = max_err + 0x06;    r33r = sse2_rnd + 0x06;
		r04r = tmp + 0x08;    r14r = two + 0x08;    r24r = max_err + 0x08;    r34r = sse2_rnd + 0x08;
		r05r = tmp + 0x0a;    r15r = two + 0x0a;    r25r = max_err + 0x0a;    r35r = sse2_rnd + 0x0a;
		r06r = tmp + 0x0c;    r16r = two + 0x0c;    r26r = max_err + 0x0c;    r36r = sse2_rnd + 0x0c;
		r07r = tmp + 0x0e;    r17r = two + 0x0e;    r27r = max_err + 0x0e;    r37r = sse2_rnd + 0x0e;
		r08r = tmp + 0x10;    r18r = two + 0x10;    r28r = max_err + 0x10;    r38r = sse2_rnd + 0x10;
		r09r = tmp + 0x12;    r19r = two + 0x12;    r29r = max_err + 0x12;    r39r = sse2_rnd + 0x12;
		r0ar = tmp + 0x14;    r1ar = two + 0x14;    r2ar = max_err + 0x14;    r3ar = sse2_rnd + 0x14;
		r0br = tmp + 0x16;    r1br = two + 0x16;    r2br = max_err + 0x16;    r3br = sse2_rnd + 0x16;
		r0cr = tmp + 0x18;    r1cr = two + 0x18;    r2cr = max_err + 0x18;    r3cr = sse2_rnd + 0x18;
		tmp = sse2_rnd + 0x1a;
		two = tmp;
		rad13_const = tmp + 0x01;	/* Leave an extra slot at radix13_const-1 for the constant two = 2.0: */
		tmp += 0x14;	/* Need 20 16-byte slots for two+sincos, but offset the carry slots by the next-larger multiple of 4 */

	  #ifdef USE_AVX
		cy00 = tmp + 0x00;
		cy04 = tmp + 0x01;
		cy08 = tmp + 0x02;
		cy12 = tmp + 0x03;
		cy16 = tmp + 0x04;
		cy20 = tmp + 0x05;
		cy24 = tmp + 0x06;
		cy28 = tmp + 0x07;
		cy32 = tmp + 0x08;
		cy36 = tmp + 0x09;
		cy40 = tmp + 0x0a;
		cy44 = tmp + 0x0b;
		cy48 = tmp + 0x0c;
		tmp += 0x0d;
		max_err = tmp + 0x00;
		sse2_rnd= tmp + 0x01;
		half_arr= tmp + 0x02;	/* This table needs 20x16 bytes for Mersenne-mod, and radixx16 for Fermat-mod */
	  #else
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
	  #endif
		ASSERT(HERE, (s1p00r == thread_arg->s1p00r), "thread-local memcheck failed!");
		ASSERT(HERE, (half_arr == thread_arg->half_arr), "thread-local memcheck failed!");
		ASSERT(HERE, (sse2_rnd->d0 == crnd && sse2_rnd->d1 == crnd), "thread-local memcheck failed!");
		tmp = half_arr;
	  #ifdef USE_AVX
		// Grab some elt of base-data [offset by, say, +32] and mpy by its inverse [+16 further]
		dtmp = (tmp+40)->d0 * (tmp+56)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+40)->d1 * (tmp+56)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #else	// SSE2:
		dtmp = (tmp+10)->d0 * (tmp+14)->d0;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
		dtmp = (tmp+10)->d1 * (tmp+14)->d1;	ASSERT(HERE, fabs(dtmp - 1.0) < EPS, "thread-local memcheck failed!");
	  #endif

		VEC_DBL_INIT(max_err, 0.0);

		sign_mask = (uint64*)(s1p00r + radix52_creals_in_local_store);
		sse_bw  = sign_mask + RE_IM_STRIDE;	// (  #doubles in a SIMD complex) x 32-bits = RE_IM_STRIDE x 64-bits
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
	#else

		// In scalar mode use these 2 ptrs to pass the base & baseinv arrays:
		base    = (double *)thread_arg->s1p00r  ;
		baseinv = (double *)thread_arg->half_arr;

	#endif	// USE_SSE2 ?

		/* Init DWT-indices: */						/* init carries	*/
	#ifdef USE_AVX
		*bjmodn00 = thread_arg->bjmodn00;		cy00->d0 = thread_arg->cy00;
		*bjmodn01 = thread_arg->bjmodn01;		cy00->d1 = thread_arg->cy01;
		*bjmodn02 = thread_arg->bjmodn02;		cy00->d2 = thread_arg->cy02;
		*bjmodn03 = thread_arg->bjmodn03;		cy00->d3 = thread_arg->cy03;
		*bjmodn04 = thread_arg->bjmodn04;		cy04->d0 = thread_arg->cy04;
		*bjmodn05 = thread_arg->bjmodn05;		cy04->d1 = thread_arg->cy05;
		*bjmodn06 = thread_arg->bjmodn06;		cy04->d2 = thread_arg->cy06;
		*bjmodn07 = thread_arg->bjmodn07;		cy04->d3 = thread_arg->cy07;
		*bjmodn08 = thread_arg->bjmodn08;		cy08->d0 = thread_arg->cy08;
		*bjmodn09 = thread_arg->bjmodn09;		cy08->d1 = thread_arg->cy09;
		*bjmodn10 = thread_arg->bjmodn10;		cy08->d2 = thread_arg->cy10;
		*bjmodn11 = thread_arg->bjmodn11;		cy08->d3 = thread_arg->cy11;
		*bjmodn12 = thread_arg->bjmodn12;		cy12->d0 = thread_arg->cy12;
		*bjmodn13 = thread_arg->bjmodn13;		cy12->d1 = thread_arg->cy13;
		*bjmodn14 = thread_arg->bjmodn14;		cy12->d2 = thread_arg->cy14;
		*bjmodn15 = thread_arg->bjmodn15;		cy12->d3 = thread_arg->cy15;
		*bjmodn16 = thread_arg->bjmodn16;		cy16->d0 = thread_arg->cy16;
		*bjmodn17 = thread_arg->bjmodn17;		cy16->d1 = thread_arg->cy17;
		*bjmodn18 = thread_arg->bjmodn18;		cy16->d2 = thread_arg->cy18;
		*bjmodn19 = thread_arg->bjmodn19;		cy16->d3 = thread_arg->cy19;
		*bjmodn20 = thread_arg->bjmodn20;		cy20->d0 = thread_arg->cy20;
		*bjmodn21 = thread_arg->bjmodn21;		cy20->d1 = thread_arg->cy21;
		*bjmodn22 = thread_arg->bjmodn22;		cy20->d2 = thread_arg->cy22;
		*bjmodn23 = thread_arg->bjmodn23;		cy20->d3 = thread_arg->cy23;
		*bjmodn24 = thread_arg->bjmodn24;		cy24->d0 = thread_arg->cy24;
		*bjmodn25 = thread_arg->bjmodn25;		cy24->d1 = thread_arg->cy25;
		*bjmodn26 = thread_arg->bjmodn26;		cy24->d2 = thread_arg->cy26;
		*bjmodn27 = thread_arg->bjmodn27;		cy24->d3 = thread_arg->cy27;
		*bjmodn28 = thread_arg->bjmodn28;		cy28->d0 = thread_arg->cy28;
		*bjmodn29 = thread_arg->bjmodn29;		cy28->d1 = thread_arg->cy29;
		*bjmodn30 = thread_arg->bjmodn30;		cy28->d2 = thread_arg->cy30;
		*bjmodn31 = thread_arg->bjmodn31;		cy28->d3 = thread_arg->cy31;
		*bjmodn32 = thread_arg->bjmodn32;		cy32->d0 = thread_arg->cy32;
		*bjmodn33 = thread_arg->bjmodn33;		cy32->d1 = thread_arg->cy33;
		*bjmodn34 = thread_arg->bjmodn34;		cy32->d2 = thread_arg->cy34;
		*bjmodn35 = thread_arg->bjmodn35;		cy32->d3 = thread_arg->cy35;
		*bjmodn36 = thread_arg->bjmodn36;		cy36->d0 = thread_arg->cy36;
		*bjmodn37 = thread_arg->bjmodn37;		cy36->d1 = thread_arg->cy37;
		*bjmodn38 = thread_arg->bjmodn38;		cy36->d2 = thread_arg->cy38;
		*bjmodn39 = thread_arg->bjmodn39;		cy36->d3 = thread_arg->cy39;
		*bjmodn40 = thread_arg->bjmodn40;		cy40->d0 = thread_arg->cy40;
		*bjmodn41 = thread_arg->bjmodn41;		cy40->d1 = thread_arg->cy41;
		*bjmodn42 = thread_arg->bjmodn42;		cy40->d2 = thread_arg->cy42;
		*bjmodn43 = thread_arg->bjmodn43;		cy40->d3 = thread_arg->cy43;
		*bjmodn44 = thread_arg->bjmodn44;		cy44->d0 = thread_arg->cy44;
		*bjmodn45 = thread_arg->bjmodn45;		cy44->d1 = thread_arg->cy45;
		*bjmodn46 = thread_arg->bjmodn46;		cy44->d2 = thread_arg->cy46;
		*bjmodn47 = thread_arg->bjmodn47;		cy44->d3 = thread_arg->cy47;
		*bjmodn48 = thread_arg->bjmodn48;		cy48->d0 = thread_arg->cy48;
		*bjmodn49 = thread_arg->bjmodn49;		cy48->d1 = thread_arg->cy49;
		*bjmodn50 = thread_arg->bjmodn50;		cy48->d2 = thread_arg->cy50;
		*bjmodn51 = thread_arg->bjmodn51;		cy48->d3 = thread_arg->cy51;

	#elif defined(USE_SSE2)

		*bjmodn00 = thread_arg->bjmodn00;		cy00->d0 = thread_arg->cy00;
		*bjmodn01 = thread_arg->bjmodn01;		cy00->d1 = thread_arg->cy01;
		*bjmodn02 = thread_arg->bjmodn02;		cy02->d0 = thread_arg->cy02;
		*bjmodn03 = thread_arg->bjmodn03;		cy02->d1 = thread_arg->cy03;
		*bjmodn04 = thread_arg->bjmodn04;		cy04->d0 = thread_arg->cy04;
		*bjmodn05 = thread_arg->bjmodn05;		cy04->d1 = thread_arg->cy05;
		*bjmodn06 = thread_arg->bjmodn06;		cy06->d0 = thread_arg->cy06;
		*bjmodn07 = thread_arg->bjmodn07;		cy06->d1 = thread_arg->cy07;
		*bjmodn08 = thread_arg->bjmodn08;		cy08->d0 = thread_arg->cy08;
		*bjmodn09 = thread_arg->bjmodn09;		cy08->d1 = thread_arg->cy09;
		*bjmodn10 = thread_arg->bjmodn10;		cy10->d0 = thread_arg->cy10;
		*bjmodn11 = thread_arg->bjmodn11;		cy10->d1 = thread_arg->cy11;
		*bjmodn12 = thread_arg->bjmodn12;		cy12->d0 = thread_arg->cy12;
		*bjmodn13 = thread_arg->bjmodn13;		cy12->d1 = thread_arg->cy13;
		*bjmodn14 = thread_arg->bjmodn14;		cy14->d0 = thread_arg->cy14;
		*bjmodn15 = thread_arg->bjmodn15;		cy14->d1 = thread_arg->cy15;
		*bjmodn16 = thread_arg->bjmodn16;		cy16->d0 = thread_arg->cy16;
		*bjmodn17 = thread_arg->bjmodn17;		cy16->d1 = thread_arg->cy17;
		*bjmodn18 = thread_arg->bjmodn18;		cy18->d0 = thread_arg->cy18;
		*bjmodn19 = thread_arg->bjmodn19;		cy18->d1 = thread_arg->cy19;
		*bjmodn20 = thread_arg->bjmodn20;		cy20->d0 = thread_arg->cy20;
		*bjmodn21 = thread_arg->bjmodn21;		cy20->d1 = thread_arg->cy21;
		*bjmodn22 = thread_arg->bjmodn22;		cy22->d0 = thread_arg->cy22;
		*bjmodn23 = thread_arg->bjmodn23;		cy22->d1 = thread_arg->cy23;
		*bjmodn24 = thread_arg->bjmodn24;		cy24->d0 = thread_arg->cy24;
		*bjmodn25 = thread_arg->bjmodn25;		cy24->d1 = thread_arg->cy25;
		*bjmodn26 = thread_arg->bjmodn26;		cy26->d0 = thread_arg->cy26;
		*bjmodn27 = thread_arg->bjmodn27;		cy26->d1 = thread_arg->cy27;
		*bjmodn28 = thread_arg->bjmodn28;		cy28->d0 = thread_arg->cy28;
		*bjmodn29 = thread_arg->bjmodn29;		cy28->d1 = thread_arg->cy29;
		*bjmodn30 = thread_arg->bjmodn30;		cy30->d0 = thread_arg->cy30;
		*bjmodn31 = thread_arg->bjmodn31;		cy30->d1 = thread_arg->cy31;
		*bjmodn32 = thread_arg->bjmodn32;		cy32->d0 = thread_arg->cy32;
		*bjmodn33 = thread_arg->bjmodn33;		cy32->d1 = thread_arg->cy33;
		*bjmodn34 = thread_arg->bjmodn34;		cy34->d0 = thread_arg->cy34;
		*bjmodn35 = thread_arg->bjmodn35;		cy34->d1 = thread_arg->cy35;
		*bjmodn36 = thread_arg->bjmodn36;		cy36->d0 = thread_arg->cy36;
		*bjmodn37 = thread_arg->bjmodn37;		cy36->d1 = thread_arg->cy37;
		*bjmodn38 = thread_arg->bjmodn38;		cy38->d0 = thread_arg->cy38;
		*bjmodn39 = thread_arg->bjmodn39;		cy38->d1 = thread_arg->cy39;
		*bjmodn40 = thread_arg->bjmodn40;		cy40->d0 = thread_arg->cy40;
		*bjmodn41 = thread_arg->bjmodn41;		cy40->d1 = thread_arg->cy41;
		*bjmodn42 = thread_arg->bjmodn42;		cy42->d0 = thread_arg->cy42;
		*bjmodn43 = thread_arg->bjmodn43;		cy42->d1 = thread_arg->cy43;
		*bjmodn44 = thread_arg->bjmodn44;		cy44->d0 = thread_arg->cy44;
		*bjmodn45 = thread_arg->bjmodn45;		cy44->d1 = thread_arg->cy45;
		*bjmodn46 = thread_arg->bjmodn46;		cy46->d0 = thread_arg->cy46;
		*bjmodn47 = thread_arg->bjmodn47;		cy46->d1 = thread_arg->cy47;
		*bjmodn48 = thread_arg->bjmodn48;		cy48->d0 = thread_arg->cy48;
		*bjmodn49 = thread_arg->bjmodn49;		cy48->d1 = thread_arg->cy49;
		*bjmodn50 = thread_arg->bjmodn50;		cy50->d0 = thread_arg->cy50;
		*bjmodn51 = thread_arg->bjmodn51;		cy50->d1 = thread_arg->cy51;

	#else

		bjmodn00 = thread_arg->bjmodn00;		cy00 = thread_arg->cy00;
		bjmodn01 = thread_arg->bjmodn01;		cy01 = thread_arg->cy01;
		bjmodn02 = thread_arg->bjmodn02;		cy02 = thread_arg->cy02;
		bjmodn03 = thread_arg->bjmodn03;		cy03 = thread_arg->cy03;
		bjmodn04 = thread_arg->bjmodn04;		cy04 = thread_arg->cy04;
		bjmodn05 = thread_arg->bjmodn05;		cy05 = thread_arg->cy05;
		bjmodn06 = thread_arg->bjmodn06;		cy06 = thread_arg->cy06;
		bjmodn07 = thread_arg->bjmodn07;		cy07 = thread_arg->cy07;
		bjmodn08 = thread_arg->bjmodn08;		cy08 = thread_arg->cy08;
		bjmodn09 = thread_arg->bjmodn09;		cy09 = thread_arg->cy09;
		bjmodn10 = thread_arg->bjmodn10;		cy10 = thread_arg->cy10;
		bjmodn11 = thread_arg->bjmodn11;		cy11 = thread_arg->cy11;
		bjmodn12 = thread_arg->bjmodn12;		cy12 = thread_arg->cy12;
		bjmodn13 = thread_arg->bjmodn13;		cy13 = thread_arg->cy13;
		bjmodn14 = thread_arg->bjmodn14;		cy14 = thread_arg->cy14;
		bjmodn15 = thread_arg->bjmodn15;		cy15 = thread_arg->cy15;
		bjmodn16 = thread_arg->bjmodn16;		cy16 = thread_arg->cy16;
		bjmodn17 = thread_arg->bjmodn17;		cy17 = thread_arg->cy17;
		bjmodn18 = thread_arg->bjmodn18;		cy18 = thread_arg->cy18;
		bjmodn19 = thread_arg->bjmodn19;		cy19 = thread_arg->cy19;
		bjmodn20 = thread_arg->bjmodn20;		cy20 = thread_arg->cy20;
		bjmodn21 = thread_arg->bjmodn21;		cy21 = thread_arg->cy21;
		bjmodn22 = thread_arg->bjmodn22;		cy22 = thread_arg->cy22;
		bjmodn23 = thread_arg->bjmodn23;		cy23 = thread_arg->cy23;
		bjmodn24 = thread_arg->bjmodn24;		cy24 = thread_arg->cy24;
		bjmodn25 = thread_arg->bjmodn25;		cy25 = thread_arg->cy25;
		bjmodn26 = thread_arg->bjmodn26;		cy26 = thread_arg->cy26;
		bjmodn27 = thread_arg->bjmodn27;		cy27 = thread_arg->cy27;
		bjmodn28 = thread_arg->bjmodn28;		cy28 = thread_arg->cy28;
		bjmodn29 = thread_arg->bjmodn29;		cy29 = thread_arg->cy29;
		bjmodn30 = thread_arg->bjmodn30;		cy30 = thread_arg->cy30;
		bjmodn31 = thread_arg->bjmodn31;		cy31 = thread_arg->cy31;
		bjmodn32 = thread_arg->bjmodn32;		cy32 = thread_arg->cy32;
		bjmodn33 = thread_arg->bjmodn33;		cy33 = thread_arg->cy33;
		bjmodn34 = thread_arg->bjmodn34;		cy34 = thread_arg->cy34;
		bjmodn35 = thread_arg->bjmodn35;		cy35 = thread_arg->cy35;
		bjmodn36 = thread_arg->bjmodn36;		cy36 = thread_arg->cy36;
		bjmodn37 = thread_arg->bjmodn37;		cy37 = thread_arg->cy37;
		bjmodn38 = thread_arg->bjmodn38;		cy38 = thread_arg->cy38;
		bjmodn39 = thread_arg->bjmodn39;		cy39 = thread_arg->cy39;
		bjmodn40 = thread_arg->bjmodn40;		cy40 = thread_arg->cy40;
		bjmodn41 = thread_arg->bjmodn41;		cy41 = thread_arg->cy41;
		bjmodn42 = thread_arg->bjmodn42;		cy42 = thread_arg->cy42;
		bjmodn43 = thread_arg->bjmodn43;		cy43 = thread_arg->cy43;
		bjmodn44 = thread_arg->bjmodn44;		cy44 = thread_arg->cy44;
		bjmodn45 = thread_arg->bjmodn45;		cy45 = thread_arg->cy45;
		bjmodn46 = thread_arg->bjmodn46;		cy46 = thread_arg->cy46;
		bjmodn47 = thread_arg->bjmodn47;		cy47 = thread_arg->cy47;
		bjmodn48 = thread_arg->bjmodn48;		cy48 = thread_arg->cy48;
		bjmodn49 = thread_arg->bjmodn49;		cy49 = thread_arg->cy49;
		bjmodn50 = thread_arg->bjmodn50;		cy50 = thread_arg->cy50;
		bjmodn51 = thread_arg->bjmodn51;		cy51 = thread_arg->cy51;

	#endif

		for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
		{
			for(j = jstart; j < jhi; j += stride)	// Stride = 4 reals for SSE2, 8 for AVX
			{
				j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
				j2 = j1+RE_IM_STRIDE;

			#ifdef USE_AVX
				/* Outputs in AVX mode are temps 2*13*32 = 26*32 = 0x340 bytes apart: */
				// Reorder blocks to yield sequentially increasing a-array offsets:
	/* Block 1 : */	add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00r, 0x340)
	/* Block 4 : */	add3 = &a[j1+p04];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03r, 0x340)
	/* Block 7 : */	add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06r, 0x340)
	/* Block 10: */	add2 = &a[j1+p12];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09r, 0x340)
	/* Block 13: */	add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0cr, 0x340)
	/* Block 3 : */	add3 = &a[j1+p20];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02r, 0x340)
	/* Block 6 : */	add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05r, 0x340)
	/* Block 9 : */	add2 = &a[j1+p28];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08r, 0x340)
	/* Block 12: */	add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0br, 0x340)
	/* Block 2 : */	add3 = &a[j1+p36];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01r, 0x340)
	/* Block 5 : */	add1 = &a[j1+p40];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04r, 0x340)
	/* Block 8 : */	add2 = &a[j1+p44];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07r, 0x340)
	/* Block 11: */	add0 = &a[j1+p48];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0ar, 0x340)

				/* Radix-13 DFT uses adjacent temps for inputs, outputs are (cyclic) with pXXr having XX += 40 (40*32 bytes = 0x500) or XX -= 12 (-12*32 bytes = -0x180) between successive outputs: */
																														/*  a1p00r,a1p31r,a1p22r,a1p13r,a1p04r,a1p35r,a1p26r,a1p17r,a1p08r,a1p39r,a1p2ar,a1p1br,a1p0cr */
				SSE2_RADIX_13_DFT(rad13_const, r00r,0x040,0x080,0x0c0,0x100,0x140,0x180,0x1c0,0x200,0x240,0x280,0x2c0,0x300, s1p00r, 0xa00, 0x700, 0x400, 0x100, 0xb00, 0x800, 0x500, 0x200, 0xc00, 0x900, 0x600, 0x300)
																														/*  a1p30r,a1p21r,a1p12r,a1p03r,a1p34r,a1p25r,a1p16r,a1p07r,a1p38r,a1p29r,a1p1ar,a1p0br,a1p3cr */
				SSE2_RADIX_13_DFT(rad13_const, r10r,0x040,0x080,0x0c0,0x100,0x140,0x180,0x1c0,0x200,0x240,0x280,0x2c0,0x300, s1p30r,-0x300,-0x600,-0x900, 0x100,-0x200,-0x500,-0x800, 0x200,-0x100,-0x400,-0x700, 0x300)
																														/*  a1p20r,a1p11r,a1p02r,a1p33r,a1p24r,a1p15r,a1p06r,a1p37r,a1p28r,a1p19r,a1p0ar,a1p3br,a1p2cr */
				SSE2_RADIX_13_DFT(rad13_const, r20r,0x040,0x080,0x0c0,0x100,0x140,0x180,0x1c0,0x200,0x240,0x280,0x2c0,0x300, s1p20r,-0x300,-0x600, 0x400, 0x100,-0x200,-0x500, 0x500, 0x200,-0x100,-0x400, 0x600, 0x300)
																														/*  a1p10r,a1p01r,a1p32r,a1p23r,a1p14r,a1p05r,a1p36r,a1p27r,a1p18r,a1p09r,a1p3ar,a1p2br,a1p1cr */
				SSE2_RADIX_13_DFT(rad13_const, r30r,0x040,0x080,0x0c0,0x100,0x140,0x180,0x1c0,0x200,0x240,0x280,0x2c0,0x300, s1p10r,-0x300, 0x700, 0x400, 0x100,-0x200, 0x800, 0x500, 0x200,-0x100, 0x900, 0x600, 0x300)

			#elif defined(USE_SSE2)

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

				AVX_cmplx_carry_norm_errcheck0_X4(s1p00r,add1,add2,add3,cy00,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p04r,add1,add2,add3,cy04,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p08r,add1,add2,add3,cy08,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p0cr,add1,add2,add3,cy12,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p13r,add1,add2,add3,cy16,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p17r,add1,add2,add3,cy20,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p1br,add1,add2,add3,cy24,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p22r,add1,add2,add3,cy28,bjmodn28,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p26r,add1,add2,add3,cy32,bjmodn32,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p2ar,add1,add2,add3,cy36,bjmodn36,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p31r,add1,add2,add3,cy40,bjmodn40,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p35r,add1,add2,add3,cy44,bjmodn44,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				AVX_cmplx_carry_norm_errcheck1_X4(s1p39r,add1,add2,add3,cy48,bjmodn48,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);

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

			/************ See the radix16_ditN_cy_dif1 routine for details on how the SSE2 carry stuff works **********/

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

				i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
				co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

			#endif	// USE_AVX?

			/*...The radix-52 DIF pass is here:	*/
			#ifdef USE_AVX
				/*...gather the needed data (52 64-bit complex, i.e. 104 64-bit reals) and do 4 radix-13 transforms...*/
				/* Radix-13 DFT inputs are (cyclic) with pXXr having XX += 48 (48*32 bytes = +0x600) or XX -= 4 (4*32 bytes = -0x080), outputs are adjacent 32-byte-separated temps: */
				SSE2_RADIX_13_DFT(rad13_const, s1p00r, 0xc00, 0xb00, 0xa00, 0x900, 0x800, 0x700, 0x600, 0x500, 0x400, 0x300, 0x200, 0x100, r00r,0x040,0x080,0x0c0,0x100,0x140,0x180,0x1c0,0x200,0x240,0x280,0x2c0,0x300)
				SSE2_RADIX_13_DFT(rad13_const, s1p30r,-0x100,-0x200,-0x300,-0x400,-0x500,-0x600,-0x700,-0x800,-0x900, 0x300, 0x200, 0x100, r10r,0x040,0x080,0x0c0,0x100,0x140,0x180,0x1c0,0x200,0x240,0x280,0x2c0,0x300)
				SSE2_RADIX_13_DFT(rad13_const, s1p20r,-0x100,-0x200,-0x300,-0x400,-0x500,-0x600, 0x600, 0x500, 0x400, 0x300, 0x200, 0x100, r20r,0x040,0x080,0x0c0,0x100,0x140,0x180,0x1c0,0x200,0x240,0x280,0x2c0,0x300)
				SSE2_RADIX_13_DFT(rad13_const, s1p10r,-0x100,-0x200,-0x300, 0x900, 0x800, 0x700, 0x600, 0x500, 0x400, 0x300, 0x200, 0x100, r30r,0x040,0x080,0x0c0,0x100,0x140,0x180,0x1c0,0x200,0x240,0x280,0x2c0,0x300)

				/*...and now do 13 radix-4 transforms...*/
				/* Inputs in AVX mode are temps 2*13*32 = 26*32 = 0x340 bytes apart. Notice how the add* indices after the first row repeat with period 4: */
				// Reorder blocks to yield sequentially increasing a-array offsets:
	/* Block 01 : */	add0 = &a[j1    ];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r00r, 0x340)
	/* Block 13 : */	add3 = &a[j1+p04];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0cr, 0x340)
	/* Block 12 : */	add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0br, 0x340)
	/* Block 11 : */	add2 = &a[j1+p12];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r0ar, 0x340)
	/* Block 10 : */	add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r09r, 0x340)
	/* Block 09 : */	add3 = &a[j1+p20];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r08r, 0x340)
	/* Block 08 : */	add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r07r, 0x340)
	/* Block 07 : */	add2 = &a[j1+p28];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r06r, 0x340)
	/* Block 06 : */	add0 = &a[j1+p32];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r05r, 0x340)
	/* Block 05 : */	add3 = &a[j1+p36];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r04r, 0x340)
	/* Block 04 : */	add1 = &a[j1+p40];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r03r, 0x340)
	/* Block 03 : */	add2 = &a[j1+p44];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r02r, 0x340)
	/* Block 02 : */	add0 = &a[j1+p48];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(add0, add1, add2, add3, r01r, 0x340)

			#elif defined(USE_SSE2)

				/*...gather the needed data (52 64-bit complex, i.e. 104 64-bit reals) and do 4 radix-13 transforms...*/
				/* Radix-13 DFT inputs are (cyclic) with pXXr having XX += 48 (48*32 bytes = +0x600) or XX -= 4 (4*32 bytes = -0x080), outputs are adjacent 32-byte-separated temps: */
				SSE2_RADIX_13_DFT(rad13_const, s1p00r, 0x600, 0x580, 0x500, 0x480, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, r00r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180)
				SSE2_RADIX_13_DFT(rad13_const, s1p30r,-0x080,-0x100,-0x180,-0x200,-0x280,-0x300,-0x380,-0x400,-0x480, 0x180, 0x100, 0x080, r10r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180)
				SSE2_RADIX_13_DFT(rad13_const, s1p20r,-0x080,-0x100,-0x180,-0x200,-0x280,-0x300, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, r20r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180)
				SSE2_RADIX_13_DFT(rad13_const, s1p10r,-0x080,-0x100,-0x180, 0x480, 0x400, 0x380, 0x300, 0x280, 0x200, 0x180, 0x100, 0x080, r30r,0x020,0x040,0x060,0x080,0x0a0,0x0c0,0x0e0,0x100,0x120,0x140,0x160,0x180)

				/*...and now do 13 radix-4 transforms...*/
				/* Inputs in SSE2 mode are temps 2*13*16 = 26*16 = 0x1a0 bytes apart. Notice how the add* indices after the first row repeat with period 4: */
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

			#else

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

			#endif	// AVX or SSE2?

			}	/* end for(j=_jstart; j < _jhi; j += 2) */

			jstart += nwt;
			jhi    += nwt;

			col += RADIX;
			co3 -= RADIX;
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
	#ifdef USE_AVX

		thread_arg->cy00 = cy00->d0;
		thread_arg->cy01 = cy00->d1;
		thread_arg->cy02 = cy00->d2;
		thread_arg->cy03 = cy00->d3;
		thread_arg->cy04 = cy04->d0;
		thread_arg->cy05 = cy04->d1;
		thread_arg->cy06 = cy04->d2;
		thread_arg->cy07 = cy04->d3;
		thread_arg->cy08 = cy08->d0;
		thread_arg->cy09 = cy08->d1;
		thread_arg->cy10 = cy08->d2;
		thread_arg->cy11 = cy08->d3;
		thread_arg->cy12 = cy12->d0;
		thread_arg->cy13 = cy12->d1;
		thread_arg->cy14 = cy12->d2;
		thread_arg->cy15 = cy12->d3;
		thread_arg->cy16 = cy16->d0;
		thread_arg->cy17 = cy16->d1;
		thread_arg->cy18 = cy16->d2;
		thread_arg->cy19 = cy16->d3;
		thread_arg->cy20 = cy20->d0;
		thread_arg->cy21 = cy20->d1;
		thread_arg->cy22 = cy20->d2;
		thread_arg->cy23 = cy20->d3;
		thread_arg->cy24 = cy24->d0;
		thread_arg->cy25 = cy24->d1;
		thread_arg->cy26 = cy24->d2;
		thread_arg->cy27 = cy24->d3;
		thread_arg->cy28 = cy28->d0;
		thread_arg->cy29 = cy28->d1;
		thread_arg->cy30 = cy28->d2;
		thread_arg->cy31 = cy28->d3;
		thread_arg->cy32 = cy32->d0;
		thread_arg->cy33 = cy32->d1;
		thread_arg->cy34 = cy32->d2;
		thread_arg->cy35 = cy32->d3;
		thread_arg->cy36 = cy36->d0;
		thread_arg->cy37 = cy36->d1;
		thread_arg->cy38 = cy36->d2;
		thread_arg->cy39 = cy36->d3;
		thread_arg->cy40 = cy40->d0;
		thread_arg->cy41 = cy40->d1;
		thread_arg->cy42 = cy40->d2;
		thread_arg->cy43 = cy40->d3;
		thread_arg->cy44 = cy44->d0;
		thread_arg->cy45 = cy44->d1;
		thread_arg->cy46 = cy44->d2;
		thread_arg->cy47 = cy44->d3;
		thread_arg->cy48 = cy48->d0;
		thread_arg->cy49 = cy48->d1;
		thread_arg->cy50 = cy48->d2;
		thread_arg->cy51 = cy48->d3;
		maxerr = MAX( MAX(max_err->d0,max_err->d1) , MAX(max_err->d2,max_err->d3) );

	#elif defined(USE_SSE2)

		thread_arg->cy00 = cy00->d0;
		thread_arg->cy01 = cy00->d1;
		thread_arg->cy02 = cy02->d0;
		thread_arg->cy03 = cy02->d1;
		thread_arg->cy04 = cy04->d0;
		thread_arg->cy05 = cy04->d1;
		thread_arg->cy06 = cy06->d0;
		thread_arg->cy07 = cy06->d1;
		thread_arg->cy08 = cy08->d0;
		thread_arg->cy09 = cy08->d1;
		thread_arg->cy10 = cy10->d0;
		thread_arg->cy11 = cy10->d1;
		thread_arg->cy12 = cy12->d0;
		thread_arg->cy13 = cy12->d1;
		thread_arg->cy14 = cy14->d0;
		thread_arg->cy15 = cy14->d1;
		thread_arg->cy16 = cy16->d0;
		thread_arg->cy17 = cy16->d1;
		thread_arg->cy18 = cy18->d0;
		thread_arg->cy19 = cy18->d1;
		thread_arg->cy20 = cy20->d0;
		thread_arg->cy21 = cy20->d1;
		thread_arg->cy22 = cy22->d0;
		thread_arg->cy23 = cy22->d1;
		thread_arg->cy24 = cy24->d0;
		thread_arg->cy25 = cy24->d1;
		thread_arg->cy26 = cy26->d0;
		thread_arg->cy27 = cy26->d1;
		thread_arg->cy28 = cy28->d0;
		thread_arg->cy29 = cy28->d1;
		thread_arg->cy30 = cy30->d0;
		thread_arg->cy31 = cy30->d1;
		thread_arg->cy32 = cy32->d0;
		thread_arg->cy33 = cy32->d1;
		thread_arg->cy34 = cy34->d0;
		thread_arg->cy35 = cy34->d1;
		thread_arg->cy36 = cy36->d0;
		thread_arg->cy37 = cy36->d1;
		thread_arg->cy38 = cy38->d0;
		thread_arg->cy39 = cy38->d1;
		thread_arg->cy40 = cy40->d0;
		thread_arg->cy41 = cy40->d1;
		thread_arg->cy42 = cy42->d0;
		thread_arg->cy43 = cy42->d1;
		thread_arg->cy44 = cy44->d0;
		thread_arg->cy45 = cy44->d1;
		thread_arg->cy46 = cy46->d0;
		thread_arg->cy47 = cy46->d1;
		thread_arg->cy48 = cy48->d0;
		thread_arg->cy49 = cy48->d1;
		thread_arg->cy50 = cy50->d0;
		thread_arg->cy51 = cy50->d1;
		maxerr = MAX(max_err->d0,max_err->d1);

	#else

		thread_arg->cy00 = cy00;
		thread_arg->cy01 = cy01;
		thread_arg->cy02 = cy02;
		thread_arg->cy03 = cy03;
		thread_arg->cy04 = cy04;
		thread_arg->cy05 = cy05;
		thread_arg->cy06 = cy06;
		thread_arg->cy07 = cy07;
		thread_arg->cy08 = cy08;
		thread_arg->cy09 = cy09;
		thread_arg->cy10 = cy10;
		thread_arg->cy11 = cy11;
		thread_arg->cy12 = cy12;
		thread_arg->cy13 = cy13;
		thread_arg->cy14 = cy14;
		thread_arg->cy15 = cy15;
		thread_arg->cy16 = cy16;
		thread_arg->cy17 = cy17;
		thread_arg->cy18 = cy18;
		thread_arg->cy19 = cy19;
		thread_arg->cy20 = cy20;
		thread_arg->cy21 = cy21;
		thread_arg->cy22 = cy22;
		thread_arg->cy23 = cy23;
		thread_arg->cy24 = cy24;
		thread_arg->cy25 = cy25;
		thread_arg->cy26 = cy26;
		thread_arg->cy27 = cy27;
		thread_arg->cy28 = cy28;
		thread_arg->cy29 = cy29;
		thread_arg->cy30 = cy30;
		thread_arg->cy31 = cy31;
		thread_arg->cy32 = cy32;
		thread_arg->cy33 = cy33;
		thread_arg->cy34 = cy34;
		thread_arg->cy35 = cy35;
		thread_arg->cy36 = cy36;
		thread_arg->cy37 = cy37;
		thread_arg->cy38 = cy38;
		thread_arg->cy39 = cy39;
		thread_arg->cy40 = cy40;
		thread_arg->cy41 = cy41;
		thread_arg->cy42 = cy42;
		thread_arg->cy43 = cy43;
		thread_arg->cy44 = cy44;
		thread_arg->cy45 = cy45;
		thread_arg->cy46 = cy46;
		thread_arg->cy47 = cy47;
		thread_arg->cy48 = cy48;
		thread_arg->cy49 = cy49;
		thread_arg->cy50 = cy50;
		thread_arg->cy51 = cy51;

	#endif	// SSE2 or AVX?

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(thread_arg->maxerr < maxerr)
		{
			thread_arg->maxerr = maxerr;
		}

		return 0x0;
	}
#endif

