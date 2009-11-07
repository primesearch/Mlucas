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

#ifdef USE_SSE2

	#undef DEBUG_SSE2
//	#define DEBUG_SSE2

	#if defined(COMPILER_TYPE_MSVC)

		#include "sse2_macro.h"

		/* DIT radix-4 subconvolution, sans twiddles, inputs in __i0-3, outputs in __o0-3, possibly coincident with inputs: */
		#define SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(__i0,__i1,__i2,__i3, __o0,__o1,__o2,__o3)\
		{\
			__asm	mov	eax, __i0\
			__asm	mov	ebx, __i1\
			__asm	mov	ecx, __i2\
			__asm	mov	edx, __i3\
			\
			__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */				__asm	movaps	xmm4,[ecx     ]	/* a[jt+p2] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p2] */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[jt   ] */	__asm	movaps	xmm6,xmm4		/* xmm4 <- cpy a[jt+p2] */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[jp   ] */	__asm	movaps	xmm7,xmm5		/* xmm5 <- cpy a[jp+p2] */\
			\
			__asm	addpd	xmm0,[ebx     ]	/* t1 */					__asm	addpd	xmm4,[edx     ]	/* t5 */\
			__asm	addpd	xmm1,[ebx+0x10]	/* t2 */					__asm	addpd	xmm5,[edx+0x10]	/* t6 */\
			__asm	subpd	xmm2,[ebx     ]	/* t3 */					__asm	subpd	xmm6,[edx     ]	/* t7 */\
			__asm	subpd	xmm3,[ebx+0x10]	/* t4 */					__asm	subpd	xmm7,[edx+0x10]	/* t8 */\
			\
			/* Finish radix-4 butterfly and store results into output-array slots: */\
			\
			__asm	mov	eax, __o0\
			__asm	mov	ebx, __o1\
			__asm	mov	ecx, __o2\
			__asm	mov	edx, __o3\
			\
			__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */				__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */\
			__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */				__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */\
			__asm	movaps	[ecx      ],xmm0	/* <- ~t5 */			__asm	movaps	[edx      ],xmm2	/* <- ~t7 */\
			__asm	movaps	[ecx+0x010],xmm1	/* <- ~t6 */			__asm	movaps	[ebx+0x010],xmm3	/* <- ~t4 */\
			__asm	addpd	xmm4,xmm4	/*          2*t5 */				__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm5,xmm5	/*          2*t6 */				__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */				__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */\
			__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */				__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */\
			__asm	movaps	[eax      ],xmm4	/* <- ~t1 */			__asm	movaps	[ebx      ],xmm7	/* <- ~t3 */\
			__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */			__asm	movaps	[edx+0x010],xmm6	/* <- ~t8 */\
		}

		/* DIT radix-4 subconvolution, sans twiddles, inputs in __i0-3, outputs in __o0-3, possibly coincident with inputs: */
		#define SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(__i0,__i1,__i2,__i3, __o0,__o1,__o2,__o3)\
		{\
			__asm	mov	eax, __i0\
			__asm	mov	ebx, __i1\
			__asm	mov	ecx, __i2\
			__asm	mov	edx, __i3\
			\
			__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */				__asm	movaps	xmm4,[ebx     ]	/* a[jt+p1] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */				__asm	movaps	xmm5,[ebx+0x10]	/* a[jp+p1] */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[jt   ] */	__asm	movaps	xmm6,xmm4		/* xmm4 <- cpy a[jt+p1] */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[jp   ] */	__asm	movaps	xmm7,xmm5		/* xmm5 <- cpy a[jp+p1] */\
			\
			__asm	addpd	xmm0,[ecx     ]	/* t1 */					__asm	addpd	xmm4,[edx     ]	/* t5 */\
			__asm	addpd	xmm1,[ecx+0x10]	/* t2 */					__asm	addpd	xmm5,[edx+0x10]	/* t6 */\
			__asm	subpd	xmm2,[ecx     ]	/* t3 */					__asm	subpd	xmm6,[edx     ]	/* t7 */\
			__asm	subpd	xmm3,[ecx+0x10]	/* t4 */					__asm	subpd	xmm7,[edx+0x10]	/* t8 */\
			\
			/* Finish radix-4 butterfly and store results into output-array slots: */\
			\
			__asm	mov	eax, __o0\
			__asm	mov	ebx, __o1\
			__asm	mov	ecx, __o2\
			__asm	mov	edx, __o3\
			\
			__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */				__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */\
			__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */				__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */\
			__asm	movaps	[ebx      ],xmm0	/* <- ~t5 */			__asm	movaps	[ecx      ],xmm2	/* <- ~t7 */\
			__asm	movaps	[ebx+0x010],xmm1	/* <- ~t6 */			__asm	movaps	[edx+0x010],xmm3	/* <- ~t4 */\
			__asm	addpd	xmm4,xmm4	/*          2*t5 */				__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm5,xmm5	/*          2*t6 */				__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */				__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */\
			__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */				__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */\
			__asm	movaps	[eax      ],xmm4	/* <- ~t1 */			__asm	movaps	[edx      ],xmm7	/* <- ~t3 */\
			__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */			__asm	movaps	[ecx+0x010],xmm6	/* <- ~t8 */\
		}

		/*...Radix-7 DFT: Inputs in memory locations __i0-6, outputs go into memory locations __o0-6, possibly coincident with inputs:\
		*/\
		#define SSE2_RADIX_07_DFT(__i0,__i1,__i2,__i3,__i4,__i5,__i6, __cc, __o0,__o1,__o2,__o3,__o4,__o5,__o6)\
		{\
		/*\
			t1r=A1r+A6r;	\
			t6r=A1r-A6r;	\
							\
			t2r=A2r+A5r;	\
			t5r=A2r-A5r;	\
							\
			t3r=A3r+A4r;	\
			t4r=A3r-A4r;	\
		*/\
			__asm	mov	eax, __i1	\
			__asm	mov	ebx, __i2	\
			__asm	mov	ecx, __i3	\
			__asm	mov	edx, __i4	\
			__asm	mov	esi, __i5	\
			__asm	mov	edi, __i6	\
			__asm	movaps	xmm6,[eax     ]	/* A1r */\
			__asm	movaps	xmm1,[edi     ]	/* A6r */\
			__asm	movaps	xmm5,[ebx     ]	/* A2r */\
			__asm	movaps	xmm2,[esi     ]	/* A5r */\
			__asm	movaps	xmm4,[ecx     ]	/* A3r */\
			__asm	movaps	xmm3,[edx     ]	/* A4r */\
			\
			__asm	mov	ebx, __i0	\
			__asm	subpd	xmm6,xmm1	/* t6r = A1r-A6r */\
			__asm	addpd	xmm1,xmm1	/*         2*A6r */\
			__asm	addpd	xmm1,xmm6	/* t1r = A1r+A6r */\
			\
			__asm	subpd	xmm5,xmm2	/* t5r = A2r-A5r */\
			__asm	addpd	xmm2,xmm2	/*         2*A5r */\
			__asm	addpd	xmm2,xmm5	/* t2r = A2r+A5r */\
			\
			__asm	movaps	xmm0,[ebx     ]	/* Ar0 */\
			__asm	subpd	xmm4,xmm3	/* t4r = A3r-A4r */\
			__asm	addpd	xmm3,xmm3	/*         2*A4r */\
			__asm	addpd	xmm3,xmm4	/* t3r = A3r+A4r */\
		/*\
			rt  = t1r+t2r+t3r;	\
			B0r = rt + A0r;		\
			t0r = rt*cx0 + A0r;			t3r=(t6r-t4r+t5r)*sx0;	\
			t1r = t1r-t2r;				t6r= t6r-t5r;			\
			t2r = t3r-t2r;				t5r= t4r+t5r;			\
			t3r =(t1r+t2r)*cx3;			t4r=(t5r-t6r)*sx3;		\
			t1r = t1r*cx1;				t6r= t6r*sx1;			\
			t2r = t2r*cx2;				t5r= t5r*sx2;			\
			tt  = t1r-t3r;				t6r= t4r+t6r;			\
			t2r = t2r-t3r;				t5r= t4r-t5r;			\
																\
			t1r= t0r- tt-t2r;			t4r= t3r-t6r-t5r;		\
			t2r= t0r+t2r;				t5r= t3r+t5r;			\
			t0r= t0r+ tt;				t3r= t3r+t6r;			\
		*/\
			__asm	mov	ecx, __o0	/* Assume that this might be the same address as any of i0-i6 */\
			__asm	mov	esi, __cc	\
			__asm	movaps	[esi+0x80],xmm0	/* cpy t0 into scratch sincos slot */	__asm	movaps	[esi+0x90],xmm6	/* cpy t6 into scratch sincos slot */	\
			__asm	addpd	xmm0,xmm1	/*~A0 = A0+t1 */							__asm	movaps	xmm7,xmm5	/* cpy t5 */			\
			__asm	addpd	xmm3,xmm2	/*~t3 = t3+t2 */							__asm	subpd	xmm5,xmm4	/*~t5 = t5-t4 */		\
			__asm	subpd	xmm1,xmm2	/*~t1 = t1-t2 */							__asm	subpd	xmm6,xmm7	/*~t6 = t6-t5 */		\
			__asm	addpd	xmm2,xmm2	/* 2*t2 */									__asm	addpd	xmm4,xmm7	/*~t5 = t4+t5 */		\
			__asm	addpd	xmm0,xmm3	/* B0 */									__asm	addpd	xmm5,[esi+0x90]	/* t3 = [t5-t4]+t6 */	\
			__asm	subpd	xmm3,xmm2	/*~t2 =  [t2+t3] - 2*t2 = t3-t2 */			__asm	movaps	xmm7,xmm4	/* cpy t5 */			\
			__asm	movaps	[ecx     ],xmm0	/* <-B0, xmm0 FREE */					__asm	subpd	xmm4,xmm6	/* t4 = ~t5-~t6 */		\
			__asm	movaps	xmm2,xmm1	/* cpy ~t1 */																	\
			__asm	subpd	xmm0,[esi+0x80]	/* r = B0 - t0 */						__asm	mulpd	xmm5,[esi+0x10]	/*~t3 = t3*sx0 */	\
			__asm	addpd	xmm2,xmm3	/* ~t1+~t2 */																	\
			__asm	mulpd	xmm3,[esi+0x40]	/* t2 = t2*cx2 */						__asm	mulpd	xmm4,[esi+0x70]	/*~t4 = t4*sx3 */	\
			__asm	mulpd	xmm1,[esi+0x20]	/* t1 = t1*cx1 */						__asm	mulpd	xmm6,[esi+0x30]	/*~t6 = t6*sx1 */	\
			__asm	mulpd	xmm0,[esi]     	/* ~r = r*(cx0-1) */					__asm	mulpd	xmm7,[esi+0x50]	/*~t5 = t5*sx2 */	\
			__asm	mulpd	xmm2,[esi+0x60]	/* t3 =(t1+t2)*cx3 */																		\
			__asm	addpd	xmm0,[ecx     ]	/* t0 =~r + B0 */						__asm	addpd	xmm6,xmm4	/*~t6 = t4+t6 */		\
			__asm	subpd	xmm1,xmm2	/* tt = t1-t3 */							__asm	subpd	xmm4,xmm7	/*~t5 = t4-t5, xmm7 FREE */\
			__asm	subpd	xmm3,xmm2	/* t2 = t2-t3, xmm2 FREE */					\
			__asm	mov	eax, __o1													\
			__asm	mov	ebx, __o2													\
			__asm	mov	ecx, __o3													\
			__asm	mov	edx, __o4													\
			__asm	mov	esi, __o5													\
			__asm	mov	edi, __o6													\
			__asm	movaps	xmm2,xmm0	/* cpy t0 */								__asm	movaps	xmm7,xmm5	/* cpy t3 */		\
			__asm	addpd	xmm0,xmm1	/*~t0 = t0+tt */							__asm	addpd	xmm5,xmm6	/*~t3 = t3+t6 */	\
			__asm	addpd	xmm1,xmm3	/*~tt = tt+t2 */							__asm	addpd	xmm6,xmm4	/*      t6+t5 */	\
			__asm	addpd	xmm3,xmm2	/*~t2 = t2+t0 */							__asm	addpd	xmm4,xmm7	/*~t5 = t5+t3 */	\
			__asm	subpd	xmm2,xmm1	/*~t1 = t0-tt-t2 */							__asm	subpd	xmm7,xmm6	/*~t4 = t3-t6-t5 */	\
			__asm	movaps	[eax     ],xmm0	/* B1 <- t0 */							__asm	movaps	[edi     ],xmm5	/* B6 <- t3 */	\
			__asm	movaps	[ebx     ],xmm2	/* B2 <- t1 */							__asm	movaps	[esi     ],xmm7	/* B5 <- t4 */	\
			__asm	movaps	[ecx     ],xmm3	/* B3 <- t2 */							__asm	movaps	[edx     ],xmm4	/* B4 <- t5 */	\
			\
		/************************** Imaginary Parts: ******************************************/\
			\
			__asm	mov	eax, __i1	\
			__asm	mov	ebx, __i2	\
			__asm	mov	ecx, __i3	\
			__asm	mov	edx, __i4	\
			__asm	mov	esi, __i5	\
			__asm	mov	edi, __i6	\
			__asm	movaps	xmm6,[eax+0x10]	/* A1i */\
			__asm	movaps	xmm1,[edi+0x10]	/* A6i */\
			__asm	movaps	xmm5,[ebx+0x10]	/* A2i */\
			__asm	movaps	xmm2,[esi+0x10]	/* A5i */\
			__asm	movaps	xmm4,[ecx+0x10]	/* A3i */\
			__asm	movaps	xmm3,[edx+0x10]	/* A4i */\
			\
			__asm	mov	ebx, __i0	\
			__asm	subpd	xmm6,xmm1	/* t6i = A1i-A6i */\
			__asm	addpd	xmm1,xmm1	/*         2*A6i */\
			__asm	addpd	xmm1,xmm6	/* t1i = A1i+A6i */\
			\
			__asm	subpd	xmm5,xmm2	/* t5i = A2i-A5i */\
			__asm	addpd	xmm2,xmm2	/*         2*A5i */\
			__asm	addpd	xmm2,xmm5	/* t2i = A2i+A5i */\
			\
			__asm	movaps	xmm0,[ebx+0x10]	/* Ai0 */\
			__asm	subpd	xmm4,xmm3	/* t4i = A3i-A4i */\
			__asm	addpd	xmm3,xmm3	/*         2*A4i */\
			__asm	addpd	xmm3,xmm4	/* t3i = A3i+A4i */\
		/*\
			it  = t1i+t2i+t3i;	\
			B0i = it + A0i;		\
			t0i = it*cx0 + A0i;			t3i=(t6i-t4i+t5i)*sx0;	\
			t1i = t1i-t2i;				t6i= t6i-t5i;			\
			t2i = t2i-t3i;				t5i= t4i+t5i;			\
			t3i =(t1i-t2i)*cx3;			t4i=(t5i-t6i)*sx3;		\
			t1i = t1i*cx1;				t6i= t6i*sx1;			\
			t2i = t2i*cx2;				t5i= t5i*sx2;			\
			it  = t1i-t3i;				t6i= t4i+t6i;			\
			t2i = t2i-t3i;				t5i= t4i-t5i;			\
																\
			t1i= t0i- it-t2i;			t4i= t3i-t6i-t5i;		\
			t2i= t0i+t2i;				t5i= t3i+t5i;			\
			t0i= t0i+ it;				t3i= t3i+t6i;			\
		*/\
			__asm	mov	ecx, __o0	\
			__asm	mov	esi, __cc	\
			__asm	movaps	[esi+0x80],xmm0	/* cpy t0 into scratch sincos slot */	__asm	movaps	[esi+0x90],xmm6	/* cpy t6 into scratch sincos slot */	\
			__asm	addpd	xmm0,xmm1	/*~A0 = A0+t1 */							__asm	movaps	xmm7,xmm5	/* cpy t5 */			\
			__asm	addpd	xmm3,xmm2	/*~t3 = t3+t2 */							__asm	subpd	xmm5,xmm4	/*~t5 = t5-t4 */		\
			__asm	subpd	xmm1,xmm2	/*~t1 = t1-t2 */							__asm	subpd	xmm6,xmm7	/*~t6 = t6-t5 */		\
			__asm	addpd	xmm2,xmm2	/* 2*t2 */									__asm	addpd	xmm4,xmm7	/*~t5 = t4+t5 */		\
			__asm	addpd	xmm0,xmm3	/* B0 */									__asm	addpd	xmm5,[esi+0x90]	/* t3 = [t5-t4]+t6 */	\
			__asm	subpd	xmm3,xmm2	/*~t2 =  [t2+t3] - 2*t2 = t3-t2 */			__asm	movaps	xmm7,xmm4	/* cpy t5 */			\
			__asm	movaps	[ecx+0x10],xmm0	/* <-B0, xmm0 FREE */					__asm	subpd	xmm4,xmm6	/* t4 = ~t5-~t6 */		\
			__asm	movaps	xmm2,xmm1	/* cpy ~t1 */																	\
			__asm	subpd	xmm0,[esi+0x80]	/* r = B0 - t0 */						__asm	mulpd	xmm5,[esi+0x10]	/*~t3 = t3*sx0 */	\
			__asm	addpd	xmm2,xmm3	/* ~t1+~t2 */																	\
			__asm	mulpd	xmm3,[esi+0x40]	/* t2 = t2*cx2 */						__asm	mulpd	xmm4,[esi+0x70]	/*~t4 = t4*sx3 */	\
			__asm	mulpd	xmm1,[esi+0x20]	/* t1 = t1*cx1 */						__asm	mulpd	xmm6,[esi+0x30]	/*~t6 = t6*sx1 */	\
			__asm	mulpd	xmm0,[esi]     	/* ~r = r*(cx0-1) */					__asm	mulpd	xmm7,[esi+0x50]	/*~t5 = t5*sx2 */	\
			__asm	mulpd	xmm2,[esi+0x60]	/* t3 =(t1+t2)*cx3 */																		\
			__asm	addpd	xmm0,[ecx+0x10]	/* t0 =~r + B0 */						__asm	addpd	xmm6,xmm4	/*~t6 = t4+t6 */		\
			__asm	subpd	xmm1,xmm2	/* tt = t1-t3 */							__asm	subpd	xmm4,xmm7	/*~t5 = t4-t5, xmm7 FREE */\
			__asm	subpd	xmm3,xmm2	/* t2 = t2-t3, xmm2 FREE */					\
			__asm	mov	eax, __o1													\
			__asm	mov	ebx, __o2													\
			__asm	mov	ecx, __o3													\
			__asm	movaps	xmm2,xmm0	/* cpy t0 */								__asm	movaps	xmm7,xmm5	/* cpy t3 */		\
			__asm	addpd	xmm0,xmm1	/*~t0 = t0+tt */							__asm	addpd	xmm5,xmm6	/*~t3 = t3+t6 */	\
			__asm	addpd	xmm1,xmm3	/*~tt = tt+t2 */							__asm	addpd	xmm6,xmm4	/*      t6+t5 */	\
			__asm	addpd	xmm3,xmm2	/*~t2 = t2+t0 */							__asm	addpd	xmm4,xmm7	/*~t5 = t5+t3 */	\
			__asm	subpd	xmm2,xmm1	/*~t1 = t0-tt-t2, xmm1 FREE */				__asm	subpd	xmm7,xmm6	/*~t4 = t3-t6-t5, xmm6 FREE */	\
		/*\
			B1r =t0r-t3i;					B1i*=t0i+t3r;\
			B2r =t1r-t4i;					B2i*=t1i+t4r;\
			B3r*=t2r+t5i;					B3i =t2i-t5r;\
			B4r*=t2r-t5i;					B4i =t2i+t5r;\
			B5r =t1r+t4i;					B5i*=t1i-t4r;\
			B6r =t0r+t3i;					B6i*=t0i-t3r;\
		*/\
			__asm	mov	edx, __o4\
			__asm	mov	esi, __o5\
			__asm	mov	edi, __o6\
			/* xmm1,6 FREE */\
			__asm	movaps	xmm1,[eax     ]	/* t0r */					__asm	movaps	xmm6,[edi     ]	/* t3r */					\
			__asm	subpd	xmm1,xmm5	/* B1r =t0r-t3i */				__asm	subpd	xmm0,xmm6	/* B6i =t0i-t3r */				\
			__asm	addpd	xmm5,xmm5	/*        2*t3i */				__asm	addpd	xmm6,xmm6	/*        2*t3r */				\
			__asm	addpd	xmm5,xmm1	/* B6r =t0r+t3i */				__asm	addpd	xmm6,xmm0	/* B1i =t0i+t3r */				\
			__asm	movaps	[eax     ],xmm1	/* <-B1r */					__asm	movaps	[edi+0x10],xmm0	/* <-B6i */		\
			__asm	movaps	[edi     ],xmm5	/* <-B6r */					__asm	movaps	[eax+0x10],xmm6	/* <-B1i */		\
			\
			__asm	movaps	xmm1,[ebx     ]	/* t1r */					__asm	movaps	xmm6,[esi     ]	/* t4r */					\
			__asm	subpd	xmm1,xmm7	/* B2r =t1r-t4i */				__asm	subpd	xmm2,xmm6	/* B5i =t1i-t4r */				\
			__asm	addpd	xmm7,xmm7	/*        2*t4i */				__asm	addpd	xmm6,xmm6	/*        2*t4r */				\
			__asm	addpd	xmm7,xmm1	/* B5r =t1r+t4i */				__asm	addpd	xmm6,xmm2	/* B2i =t1i+t4r */				\
			__asm	movaps	[ebx     ],xmm1	/* <-B2r */					__asm	movaps	[esi+0x10],xmm2	/* <-B5i */		\
			__asm	movaps	[esi     ],xmm7	/* <-B5r*/					__asm	movaps	[ebx+0x10],xmm6	/* <-B2i */		\
			\
			/* Note the order reversal on this pair of outputs: */\
			__asm	movaps	xmm0,[ecx     ]	/* t2r */					__asm	movaps	xmm5,[edx     ]	/* t5r */					\
			__asm	subpd	xmm0,xmm4	/* B4r =t2r-t5i */				__asm	subpd	xmm3,xmm5	/* B3i =t2i-t5r */				\
			__asm	addpd	xmm4,xmm4	/*        2*t5i */				__asm	addpd	xmm5,xmm5	/*        2*t5r */				\
			__asm	addpd	xmm4,xmm0	/* B3r =t2r+t5i */				__asm	addpd	xmm5,xmm3	/* B4i =t2i+t5r */				\
			__asm	movaps	[edx     ],xmm0	/* <-B4r */					__asm	movaps	[ecx+0x10],xmm3	/* <-B3i */		\
			__asm	movaps	[ecx     ],xmm4	/* <-B3r*/					__asm	movaps	[edx+0x10],xmm5	/* <-B4i */		\
		}

	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "radix28_ditN_cy_dif1_gcc32.h"

		#else

			#include "radix28_ditN_cy_dif1_gcc64.h"

		#endif

	#endif

#endif

#define RADIX_07_DFT_NUSS(\
	__A0r,__A0i,\
	__A1r,__A1i,\
	__A2r,__A2i,\
	__A3r,__A3i,\
	__A4r,__A4i,\
	__A5r,__A5i,\
	__A6r,__A6i,\
	__t0r,__t0i,\
	__t1r,__t1i,\
	__t2r,__t2i,\
	__t3r,__t3i,\
	__t4r,__t4i,\
	__t5r,__t5i,\
	__t6r,__t6i,\
	__B0r,__B0i,\
	__B1r,__B1i,\
	__B2r,__B2i,\
	__B3r,__B3i,\
	__B4r,__B4i,\
	__B5r,__B5i,\
	__B6r,__B6i,\
	__cx0,__sx0,\
	__cx1,__sx1,\
	__cx2,__sx2,\
	__cx3,__sx3,\
	__rt,__it)\
{\
	__t0r = __A0r;						__t0i = __A0i;				\
	__t6r = __A1r - __A6r;				__t6i = __A1i - __A6i;	/* x1 - x6 */	\
	__t1r = __A1r + __A6r;				__t1i = __A1i + __A6i;	/* x1 + x6 */	\
	\
	__t5r = __A2r - __A5r;				__t5i = __A2i - __A5i;	/* x2 - x5 */	\
	__t2r = __A2r + __A5r;				__t2i = __A2i + __A5i;	/* x2 + x5 */	\
	\
	__t4r = __A3r - __A4r;				__t4i = __A3i - __A4i;	/* x3 - x4 */	\
	__t3r = __A3r + __A4r;				__t3i = __A3i + __A4i;	/* x3 + x4 */	\
	\
	__rt = __t1r+__t2r+__t3r;			__it = __t1i+__t2i+__t3i;		\
	__B0r= __rt+__t0r;					__B0i= __it+__t0i;				\
	__t0r= __rt*__cx0+__t0r;			__t0i= __it*__cx0+__t0i;		\
	__t1r= __t1r-__t2r;					__t1i= __t1i-__t2i;				\
	__t2r= __t3r-__t2r;					__t2i= __t3i-__t2i;				\
	__t3r=(__t1r+__t2r)*__cx3;			__t3i=(__t1i+__t2i)*__cx3;		\
	__t1r= __t1r*__cx1;					__t1i= __t1i*__cx1;				\
	__t2r= __t2r*__cx2;					__t2i= __t2i*__cx2;				\
	__rt = __t1r-__t3r;					__it = __t1i-__t3i;				\
	__t2r= __t2r-__t3r;					__t2i= __t2i-__t3i;				\
																		\
	__t1r= __t0r-__rt-__t2r;			__t1i= __t0i-__it-__t2i;		\
	__t2r= __t0r+__t2r;					__t2i= __t0i+__t2i;				\
	__t0r= __t0r+__rt;					__t0i= __t0i+__it;				\
																		\
	__t3r=(__t6r-__t4r+__t5r)*__sx0;	__t3i=(__t6i-__t4i+__t5i)*__sx0;\
	__t6r= __t6r-__t5r;					__t6i= __t6i-__t5i;				\
	__t5r= __t4r+__t5r;					__t5i= __t4i+__t5i;				\
	__t4r=(__t5r-__t6r)*__sx3;			__t4i=(__t5i-__t6i)*__sx3;		\
	__t6r= __t6r*__sx1;					__t6i= __t6i*__sx1;				\
	__t5r= __t5r*__sx2;					__t5i= __t5i*__sx2;				\
	__t6r= __t4r+__t6r;					__t6i= __t4i+__t6i;				\
	__t5r= __t4r-__t5r;					__t5i= __t4i-__t5i;				\
																		\
	__t4r= __t3r-__t6r-__t5r;			__t4i= __t3i-__t6i-__t5i;		\
	__t5r= __t3r+__t5r;					__t5i= __t3i+__t5i;				\
	__t3r= __t3r+__t6r;					__t3i= __t3i+__t6i;				\
																		\
	__B1r =__t0r-__t3i;					__B1i =__t0i+__t3r;				\
	__B2r =__t1r-__t4i;					__B2i =__t1i+__t4r;				\
	__B3r =__t2r+__t5i;					__B3i =__t2i-__t5r;				\
	__B4r =__t2r-__t5i;					__B4i =__t2i+__t5r;				\
	__B5r =__t1r+__t4i;					__B5i =__t1i-__t4r;				\
	__B6r =__t0r+__t3i;					__B6i =__t0i-__t3r;				\
}

/**************/

int radix28_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-28 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-28 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix7/8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n28,i,j,j1,j2,jt,jp,jstart,jhi,full_pass,k1,k2,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p08,p12,p16,p20,p24;
	static double radix_inv, n2inv;

#if defined(USE_SSE2) || !defined (LO_ADD)
	/* SSE2 version assumes LO_ADD = 0, i.e. the low-mul Nussbaumer-style DFT implementation: */
	static double	cx0 =-0.16666666666666666667,	/* (cc1+cc2+cc3)/3 */
				 	cx1 = 1.52445866976115265675, 	/*  cc1-cc3		*/
				 	cx2 = 0.67844793394610472196, 	/*  cc2-cc3		*/
				 	cx3 = 0.73430220123575245957,	/* (cc1+cc2-2*cc3)/3	*/
				/* Switch the sign of ss3 in these: */
				 	sx0 = 0.44095855184409843174,	/* (ss1+ss2-ss3)/3	*/
				 	sx1 = 1.21571522158558792920, 	/*  ss1+ss3		*/
				 	sx2 = 1.40881165129938172752, 	/*  ss2+ss3		*/
				 	sx3 = 0.87484229096165655224;	/* (ss1+ss2+2*ss3)/3	*/
#else
	/* Non-SSE2 version assumes LO_ADD = 1 and uses the corresponding versions of the sincos constants: */
	static double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
#endif

	double scale
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t54,t55;
	double maxerr = 0.0;
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wtl,wtlp1,wtn,wtnm1;
	double wt_re,wt_im;									/* Fermat-mod weights stuff */
//	int ii00,ii01,ii02,ii03,ii04,ii05,ii06,ii07,ii08,ii09,ii10,ii11,ii12,ii13,ii14,ii15,ii16,ii17,ii18,ii19,ii20,ii21,ii22,ii23,ii24,ii25,ii26,ii27;	/* indices into weights arrays (mod NWT) */

#ifdef USE_SSE2

	int i0,i1,m0,m1,m3;	/* m2 already def'd for regular carry sequence */
	double *add0, *add1, *add2, *add3, *add4, *add5, *add6;	/* Addresses into array sections */
	static struct complex *sc_arr = 0x0, *sc_ptr;
	static struct complex *isrt2, *cc0, *ss0, *cc1, *ss1, *cc2, *ss2, *cc3, *ss3, *max_err, *sse2_rnd, *half_arr, *tmp
	,*s1p00r,*s1p01r,*s1p02r,*s1p03r,*s1p04r,*s1p05r,*s1p06r,*s1p07r,*s1p08r,*s1p09r,*s1p10r,*s1p11r,*s1p12r,*s1p13r,*s1p14r,*s1p15r,*s1p16r,*s1p17r,*s1p18r,*s1p19r,*s1p20r,*s1p21r,*s1p22r,*s1p23r,*s1p24r,*s1p25r,*s1p26r,*s1p27r
	,*s1p00i,*s1p01i,*s1p02i,*s1p03i,*s1p04i,*s1p05i,*s1p06i,*s1p07i,*s1p08i,*s1p09i,*s1p10i,*s1p11i,*s1p12i,*s1p13i,*s1p14i,*s1p15i,*s1p16i,*s1p17i,*s1p18i,*s1p19i,*s1p20i,*s1p21i,*s1p22i,*s1p23i,*s1p24i,*s1p25i,*s1p26i,*s1p27i;

	static struct complex *cy_r00,*cy_r02,*cy_r04,*cy_r06,*cy_r08,*cy_r10,*cy_r12,*cy_r14,*cy_r16,*cy_r18,*cy_r20,*cy_r22,*cy_r24,*cy_r26;
	static struct complex *cy_i00,*cy_i02,*cy_i04,*cy_i06,*cy_i08,*cy_i10,*cy_i12,*cy_i14,*cy_i16,*cy_i18,*cy_i20,*cy_i22,*cy_i24,*cy_i26;
	static uint64 *sm_arr = 0x0, *sm_ptr, *sign_mask, *sse_bw, *sse_sw, *sse_n;
	uint64 tmp64;
	static int *bjmodn00,*bjmodn01,*bjmodn02,*bjmodn03,*bjmodn04,*bjmodn05,*bjmodn06,*bjmodn07,*bjmodn08,*bjmodn09,*bjmodn10,*bjmodn11,*bjmodn12,*bjmodn13,*bjmodn14,*bjmodn15,*bjmodn16,*bjmodn17,*bjmodn18,*bjmodn19,*bjmodn20,*bjmodn21,*bjmodn22,*bjmodn23,*bjmodn24,*bjmodn25,*bjmodn26,*bjmodn27;

#else
  #if PFETCH
	double *addr, *addp;
  #endif
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	double rt,it,temp,frac;
	double wt,wtinv,wtA,wtB,wtC;
	int m,m2;
	int bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27;
	double re,im
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i
	,cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27
	,cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27;

#endif

/*...stuff for the multithreaded implementation is here:	*/
	static uint32 CY_THREADS,pini;
	int ithread,j_jhi;
	uint32 ptr_prod;
	static int *_bjmodn00 = 0x0,*_bjmodn01 = 0x0,*_bjmodn02 = 0x0,*_bjmodn03 = 0x0,*_bjmodn04 = 0x0,*_bjmodn05 = 0x0,*_bjmodn06 = 0x0,*_bjmodn07 = 0x0,*_bjmodn08 = 0x0,*_bjmodn09 = 0x0,*_bjmodn10 = 0x0,*_bjmodn11 = 0x0,*_bjmodn12 = 0x0,*_bjmodn13 = 0x0,*_bjmodn14 = 0x0,*_bjmodn15 = 0x0,*_bjmodn16 = 0x0,*_bjmodn17 = 0x0,*_bjmodn18 = 0x0,*_bjmodn19 = 0x0,*_bjmodn20 = 0x0,*_bjmodn21 = 0x0,*_bjmodn22 = 0x0,*_bjmodn23 = 0x0,*_bjmodn24 = 0x0,*_bjmodn25 = 0x0,*_bjmodn26 = 0x0,*_bjmodn27 = 0x0;
	static int *_bjmodnini = 0x0;
	static int *_i, *_jstart = 0x0, *_jhi = 0x0, *_col = 0x0, *_co2 = 0x0, *_co3 = 0x0;
	static double *_maxerr = 0x0,
	*_cy_r00 = 0x0,*_cy_r01 = 0x0,*_cy_r02 = 0x0,*_cy_r03 = 0x0,*_cy_r04 = 0x0,*_cy_r05 = 0x0,*_cy_r06 = 0x0,*_cy_r07 = 0x0,*_cy_r08 = 0x0,*_cy_r09 = 0x0,*_cy_r10 = 0x0,*_cy_r11 = 0x0,*_cy_r12 = 0x0,*_cy_r13 = 0x0,*_cy_r14 = 0x0,*_cy_r15 = 0x0,*_cy_r16 = 0x0,*_cy_r17 = 0x0,*_cy_r18 = 0x0,*_cy_r19 = 0x0,*_cy_r20 = 0x0,*_cy_r21 = 0x0,*_cy_r22 = 0x0,*_cy_r23 = 0x0,*_cy_r24 = 0x0,*_cy_r25 = 0x0,*_cy_r26 = 0x0,*_cy_r27 = 0x0,
	*_cy_i00 = 0x0,*_cy_i01 = 0x0,*_cy_i02 = 0x0,*_cy_i03 = 0x0,*_cy_i04 = 0x0,*_cy_i05 = 0x0,*_cy_i06 = 0x0,*_cy_i07 = 0x0,*_cy_i08 = 0x0,*_cy_i09 = 0x0,*_cy_i10 = 0x0,*_cy_i11 = 0x0,*_cy_i12 = 0x0,*_cy_i13 = 0x0,*_cy_i14 = 0x0,*_cy_i15 = 0x0,*_cy_i16 = 0x0,*_cy_i17 = 0x0,*_cy_i18 = 0x0,*_cy_i19 = 0x0,*_cy_i20 = 0x0,*_cy_i21 = 0x0,*_cy_i22 = 0x0,*_cy_i23 = 0x0,*_cy_i24 = 0x0,*_cy_i25 = 0x0,*_cy_i26 = 0x0,*_cy_i27 = 0x0;

/*...change n28 and n_div_wt to non-static to work around a gcc compiler bug. */
	n28   = n/28;
	n_div_nwt = n28 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n28)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/28 in radix28_ditN_cy_dif1.\n",iter);
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
			CY_THREADS = MAX_THREADS;

		ASSERT(HERE, CY_THREADS >= NTHREADS,"radix28_ditN_cy_dif1.c: CY_THREADS < NTHREADS");
		ASSERT(HERE, isPow2(CY_THREADS)    ,"radix28_ditN_cy_dif1.c: CY_THREADS not a power of 2!");
		if(CY_THREADS > 1)
		{
			ASSERT(HERE, n28      %CY_THREADS == 0,"radix28_ditN_cy_dif1.c: n28      %CY_THREADS != 0");
			ASSERT(HERE, n_div_nwt%CY_THREADS == 0,"radix28_ditN_cy_dif1.c: n_div_nwt%CY_THREADS != 0");
		}

	#ifdef MULTITHREAD
		fprintf(stderr,"Using %d threads in carry step\n", CY_THREADS);
	#endif

		ASSERT(HERE, LO_ADD,"radix28_ditN_cy_dif1.c: LO_ADD");
		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)28));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;	/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

	#ifdef USE_SSE2

		ASSERT(HERE, (MODULUS_TYPE == MODULUS_TYPE_MERSENNE), "SSE2 currently only supports Mersenne-mod!");

		sc_arr = ALLOC_COMPLEX(sc_arr, 128);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt0    & 0x3f) == 0, "wt0[]  not 64-byte aligned!");
		ASSERT(HERE, ((uint32)wt1    & 0x3f) == 0, "wt1[]  not 64-byte aligned!");

		/* Size here is [8 + radix/2 + 4] 8-byte elements */
		sm_arr = ALLOC_UINT64(sm_arr, 26);	if(!sm_arr){ sprintf(cbuf, "FATAL: unable to allocate sm_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sm_ptr = ALIGN_UINT64(sm_arr);
		ASSERT(HERE, ((uint32)sm_ptr & 0x3f) == 0, "sm_ptr not 64-byte aligned!");

	/* Use low 64 16-byte slots of sc_arr for temporaries, next 7 for the nontrivial complex 16th roots,
	next 32 for the doubled carry pairs, next 2 for ROE and RND_CONST, next 20 for the half_arr table lookup stuff,
	plus at least 3 more slots to allow for 64-byte alignment of the array:
	*/
		s1p00r = sc_ptr + 0x00;		cc0		= sc_ptr + 0x38;
		s1p00i = sc_ptr + 0x01;		ss0		= sc_ptr + 0x39;
		s1p01r = sc_ptr + 0x02;		cc1		= sc_ptr + 0x3a;
		s1p01i = sc_ptr + 0x03;		ss1		= sc_ptr + 0x3b;
		s1p02r = sc_ptr + 0x04;		cc2		= sc_ptr + 0x3c;
		s1p02i = sc_ptr + 0x05;		ss2		= sc_ptr + 0x3d;
		s1p03r = sc_ptr + 0x06;		cc3  	= sc_ptr + 0x3e;
		s1p03i = sc_ptr + 0x07;		ss3		= sc_ptr + 0x3f;	/* Extra 2 slots for scratch storage here */
		s1p04r = sc_ptr + 0x08;		cy_r00	= sc_ptr + 0x42;
		s1p04i = sc_ptr + 0x09;		cy_r02	= sc_ptr + 0x43;
		s1p05r = sc_ptr + 0x0a;		cy_r04	= sc_ptr + 0x44;
		s1p05i = sc_ptr + 0x0b;		cy_r06	= sc_ptr + 0x45;
		s1p06r = sc_ptr + 0x0c;		cy_r08	= sc_ptr + 0x46;
		s1p06i = sc_ptr + 0x0d;		cy_r10	= sc_ptr + 0x47;
		s1p07r = sc_ptr + 0x0e;		cy_r12	= sc_ptr + 0x48;
		s1p07i = sc_ptr + 0x0f;		cy_r14	= sc_ptr + 0x49;
		s1p08r = sc_ptr + 0x10;		cy_r16	= sc_ptr + 0x4a;
		s1p08i = sc_ptr + 0x11;		cy_r18	= sc_ptr + 0x4b;
		s1p09r = sc_ptr + 0x12;		cy_r20	= sc_ptr + 0x4c;
		s1p09i = sc_ptr + 0x13;		cy_r22	= sc_ptr + 0x4d;
		s1p10r = sc_ptr + 0x14;		cy_r24	= sc_ptr + 0x4e;
		s1p10i = sc_ptr + 0x15;		cy_r26	= sc_ptr + 0x4f;
		s1p11r = sc_ptr + 0x16;		cy_i00	= sc_ptr + 0x50;
		s1p11i = sc_ptr + 0x17;		cy_i02	= sc_ptr + 0x51;
		s1p12r = sc_ptr + 0x18;		cy_i04	= sc_ptr + 0x52;
		s1p12i = sc_ptr + 0x19;		cy_i06	= sc_ptr + 0x53;
		s1p13r = sc_ptr + 0x1a;		cy_i08	= sc_ptr + 0x54;
		s1p13i = sc_ptr + 0x1b;		cy_i10	= sc_ptr + 0x55;
		s1p14r = sc_ptr + 0x1c;		cy_i12	= sc_ptr + 0x56;
		s1p14i = sc_ptr + 0x1d;		cy_i14	= sc_ptr + 0x57;
		s1p15r = sc_ptr + 0x1e;		cy_i16	= sc_ptr + 0x58;
		s1p15i = sc_ptr + 0x1f;		cy_i18	= sc_ptr + 0x59;
		s1p16r = sc_ptr + 0x20;		cy_i20	= sc_ptr + 0x5a;
		s1p16i = sc_ptr + 0x21;		cy_i22	= sc_ptr + 0x5b;
		s1p17r = sc_ptr + 0x22;		cy_i24	= sc_ptr + 0x5c;
		s1p17i = sc_ptr + 0x23;		cy_i26	= sc_ptr + 0x5d;
		s1p18r = sc_ptr + 0x24;		max_err = sc_ptr + 0x5e;
		s1p18i = sc_ptr + 0x25;		sse2_rnd= sc_ptr + 0x5f;
		s1p19r = sc_ptr + 0x26;		half_arr= sc_ptr + 0x60;	/* This table needs 20x16 bytes */
		s1p19i = sc_ptr + 0x27;
		s1p20r = sc_ptr + 0x28;
		s1p20i = sc_ptr + 0x29;
		s1p21r = sc_ptr + 0x2a;
		s1p21i = sc_ptr + 0x2b;
		s1p22r = sc_ptr + 0x2c;
		s1p22i = sc_ptr + 0x2d;
		s1p23r = sc_ptr + 0x2e;
		s1p23i = sc_ptr + 0x2f;
		s1p24r = sc_ptr + 0x30;
		s1p24i = sc_ptr + 0x31;
		s1p25r = sc_ptr + 0x32;
		s1p25i = sc_ptr + 0x33;
		s1p26r = sc_ptr + 0x34;
		s1p26i = sc_ptr + 0x35;
		s1p27r = sc_ptr + 0x36;
		s1p27i = sc_ptr + 0x37;

		/* These remain fixed: */
		/* cc0 = (cc1+cc2+cc3)/3 - 1; subtract 1 from Nussbaumer's definition in order to ease in-place computation */
		cc0->re = cx0-1;cc0->im = cx0-1;	ss0->re = sx0;	ss0->im = sx0;
		cc1->re = cx1;	cc1->im = cx1;		ss1->re = sx1;	ss1->im = sx1;
		cc2->re = cx2;	cc2->im = cx2;		ss2->re = sx2;	ss2->im = sx2;
		cc3->re = cx3;	cc3->im = cx3;		ss3->re = sx3;	ss3->im = sx3;

		/* SSE2 math = 53-mantissa-bit IEEE double-float: */
		sse2_rnd->re = 3.0*0x4000000*0x2000000;
		sse2_rnd->im = 3.0*0x4000000*0x2000000;

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

	#if 0
		// Set up the quadrupled-32-bit-int SSE constants used by the carry macros:
		sse_bw  = sm_ptr + 2;
		__asm	mov	eax, bw
		__asm	mov	ebx, sse_bw
		__asm	movd	xmm0,eax	/* Move actual *value* of reg eax into low 32 bits of xmm0 */
		__asm	pshufd	xmm0,xmm0,0	/* Broadcast low 32 bits of xmm0 to all 4 slots of xmm0 */
		__asm	movaps	[ebx],xmm0

		sse_sw  = sm_ptr + 4;
		__asm	lea	eax, sw
		__asm	mov	ebx, sse_sw
		__asm	movd	xmm0,[eax]	/* Variant 2: Move contents of address pointed to by reg eax into low 32 bits of xmm0 */
		__asm	pshufd	xmm0,xmm0,0	/* Broadcast low 32 bits of xmm0 to all 4 slots of xmm0 */
		__asm	movaps	[ebx],xmm0

		sse_n   = sm_ptr + 6;
		__asm	lea	eax, n
		__asm	mov	ebx, sse_n
		__asm	movd	xmm0,[eax]
		__asm	pshufd	xmm0,xmm0,0	// Broadcast low 32 bits of xmm0 to all 4 slots of xmm0
		__asm	movaps	[ebx],xmm0
	#else
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

	#endif

		/*   constant index offsets for array load/stores are here.	*/

		p01 = n28;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p08 = p04 + p04;
		p12 = p08 + p04;
		p16 = p12 + p04;
		p20 = p16 + p04;
		p24 = p20 + p04;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );

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

			free((void *)_jstart ); _jstart  = 0x0;
			free((void *)_jhi    ); _jhi     = 0x0;
			free((void *)_maxerr); _maxerr = 0x0;
			free((void *)_col   ); _col    = 0x0;
			free((void *)_co2   ); _co2    = 0x0;
			free((void *)_co3   ); _co3    = 0x0;

			free((void *)_bjmodnini); _bjmodnini = 0x0;
		}

		ptr_prod = (uint32)0;	/* Store bitmask for allocatable-array ptrs here, check vs 0 after all alloc calls finish */
		_i       	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_i== 0x0);
		_bjmodn00	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn00== 0x0);
		_bjmodn01	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn01== 0x0);
		_bjmodn02	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn02== 0x0);
		_bjmodn03	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn03== 0x0);
		_bjmodn04	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn04== 0x0);
		_bjmodn05	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn05== 0x0);
		_bjmodn06	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn06== 0x0);
		_bjmodn07	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn07== 0x0);
		_bjmodn08	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn08== 0x0);
		_bjmodn09	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn09== 0x0);
		_bjmodn10	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn10== 0x0);
		_bjmodn11	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn11== 0x0);
		_bjmodn12	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn12== 0x0);
		_bjmodn13	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn13== 0x0);
		_bjmodn14	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn14== 0x0);
		_bjmodn15	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn15== 0x0);
		_bjmodn16	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn16== 0x0);
		_bjmodn17	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn17== 0x0);
		_bjmodn18	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn18== 0x0);
		_bjmodn19	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn19== 0x0);
		_bjmodn20	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn20== 0x0);
		_bjmodn21	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn21== 0x0);
		_bjmodn22	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn22== 0x0);
		_bjmodn23	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn23== 0x0);
		_bjmodn24	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn24== 0x0);
		_bjmodn25	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn25== 0x0);
		_bjmodn26	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn26== 0x0);
		_bjmodn27	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_bjmodn27== 0x0);
		_jstart  	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_jstart  == 0x0);
		_jhi     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_jhi     == 0x0);
		_col     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_col     == 0x0);
		_co2     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_co2     == 0x0);
		_co3     	= (int *)malloc(CY_THREADS*sizeof(int));	ptr_prod += (uint32)(_co3     == 0x0);

		_cy_r00	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r00== 0x0);
		_cy_r01	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r01== 0x0);
		_cy_r02	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r02== 0x0);
		_cy_r03	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r03== 0x0);
		_cy_r04	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r04== 0x0);
		_cy_r05	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r05== 0x0);
		_cy_r06	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r06== 0x0);
		_cy_r07	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r07== 0x0);
		_cy_r08	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r08== 0x0);
		_cy_r09	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r09== 0x0);
		_cy_r10	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r10== 0x0);
		_cy_r11	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r11== 0x0);
		_cy_r12	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r12== 0x0);
		_cy_r13	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r13== 0x0);
		_cy_r14	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r14== 0x0);
		_cy_r15	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r15== 0x0);
		_cy_r16	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r16== 0x0);
		_cy_r17	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r17== 0x0);
		_cy_r18	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r18== 0x0);
		_cy_r19	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r19== 0x0);
		_cy_r20	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r20== 0x0);
		_cy_r21	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r21== 0x0);
		_cy_r22	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r22== 0x0);
		_cy_r23	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r23== 0x0);
		_cy_r24	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r24== 0x0);
		_cy_r25	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r25== 0x0);
		_cy_r26	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r26== 0x0);
		_cy_r27	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_r27== 0x0);

		_cy_i00	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i00== 0x0);
		_cy_i01	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i01== 0x0);
		_cy_i02	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i02== 0x0);
		_cy_i03	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i03== 0x0);
		_cy_i04	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i04== 0x0);
		_cy_i05	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i05== 0x0);
		_cy_i06	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i06== 0x0);
		_cy_i07	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i07== 0x0);
		_cy_i08	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i08== 0x0);
		_cy_i09	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i09== 0x0);
		_cy_i10	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i10== 0x0);
		_cy_i11	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i11== 0x0);
		_cy_i12	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i12== 0x0);
		_cy_i13	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i13== 0x0);
		_cy_i14	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i14== 0x0);
		_cy_i15	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i15== 0x0);
		_cy_i16	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i16== 0x0);
		_cy_i17	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i17== 0x0);
		_cy_i18	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i18== 0x0);
		_cy_i19	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i19== 0x0);
		_cy_i20	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i20== 0x0);
		_cy_i21	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i21== 0x0);
		_cy_i22	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i22== 0x0);
		_cy_i23	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i23== 0x0);
		_cy_i24	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i24== 0x0);
		_cy_i25	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i25== 0x0);
		_cy_i26	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i26== 0x0);
		_cy_i27	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_cy_i27== 0x0);

		_maxerr	= (double *)malloc(CY_THREADS*sizeof(double));	ptr_prod += (uint32)(_maxerr== 0x0);

		ASSERT(HERE, ptr_prod == 0, "FATAL: unable to allocate one or more auxiliary arrays in radix28_ditN_cy_dif1.");

		/* Create (THREADS + 1) copies of _bjmodnini and use the extra (uppermost) one to store the "master" increment,
		i.e. the one that n2/28-separated FFT outputs need:
		*/
		_bjmodnini = (int *)malloc((CY_THREADS + 1)*sizeof(int));	if(!_bjmodnini){ sprintf(cbuf,"FATAL: unable to allocate array _bjmodnini in radix28_ditN_cy_dif1.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		_bjmodnini[0] = 0;
		_bjmodnini[1] = 0;

		/* For Fermat-mod, since 'adjacent' words are actually stride-2 separated
		in terms of the floating residue array, block boundaries have half the i-index
		(e.g. as in sw*i and bw*i) value they do in the Mersenne-mod case:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			jhi = n28/CY_THREADS;
		}
		else
		{
			jhi = n28/CY_THREADS/2;
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
		for(j=0; j < n28; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
		ASSERT(HERE, _bjmodnini[CY_THREADS] == bjmodnini,"_bjmodnini[CY_THREADS] != bjmodnini");

	}	/* endif(first_entry) */

/*...The radix-28 final DIT pass is here.	*/

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
		_i[0] = 1;		/* Pointer to the BASE and BASEINV arrays. If n does not divide p, lowest-order digit is always a bigword (_i[0] = 1).	*/

		/*
		Moved this inside the outer-loop, so on cleanup pass can use it to reset _col,_co2,_co3 starting values,
		then simply overwrite it with 1 prior to starting the k-loop.
		*/
		khi = n_div_nwt/CY_THREADS;

		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_jstart[ithread] = ithread*n28/CY_THREADS;
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + 7;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + nwt-1;

			_col[ithread] = ithread*(khi*28);			/* col gets incremented by RADIX_VEC[0] on every pass through the k-loop */
			_co2[ithread] = (n>>nwt_bits)-1+28 - _col[ithread];	/* co2 gets decremented by RADIX_VEC[0] on every pass through the k-loop */
			_co3[ithread] = _co2[ithread]-28;			/* At the start of each new j-loop, co3=co2-RADIX_VEC[0]	*/
		}
	}
	else
	{
		_i[0] = 0;		/* Pointer to the BASE and BASEINV arrays. If n divides p, lowest-order digit is always a smallword (_i[0] = 0).	*/

		khi = 1;

		for(ithread = 0; ithread < CY_THREADS; ithread++)
		{
			_jstart[ithread] = ithread*n28/CY_THREADS;
			/*
			For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
			*/
			if(!full_pass)
				_jhi[ithread] = _jstart[ithread] + 15;		/* Cleanup loop assumes carryins propagate at most 4 words up. */
			else
				_jhi[ithread] = _jstart[ithread] + n_div_nwt/CY_THREADS;
		}
#if 0
		/* For Fermat-mod, IBDWT access patterns repeat with period NWT = {odd part of radix0},
		so for even radix0's only really need that many bjmodn and ii's, but that would require
		specialized carry macros that don't update ii and bjmodn - not worth the trouble.
		*/
		/* indices into IBDWT weights arrays (mod NWT) is here: */
		ii00= 0;
		ii01= (SW_DIV_N*n28/2) % nwt;
		ii02= (ii01+ ii01) % nwt;
		ii03= (ii02+ ii01) % nwt;
		ii04= (ii03+ ii01) % nwt;
		ii05= (ii04+ ii01) % nwt;
		ii06= (ii05+ ii01) % nwt;
		ii07= (ii06+ ii01) % nwt;
		ii08= (ii07+ ii01) % nwt;
		ii09= (ii08+ ii01) % nwt;
		ii10= (ii09+ ii01) % nwt;
		ii11= (ii10+ ii01) % nwt;
		ii12= (ii11+ ii01) % nwt;
		ii13= (ii12+ ii01) % nwt;
		ii14= (ii13+ ii01) % nwt;
		ii15= (ii14+ ii01) % nwt;
		ii16= (ii15+ ii01) % nwt;
		ii17= (ii16+ ii01) % nwt;
		ii18= (ii17+ ii01) % nwt;
		ii19= (ii18+ ii01) % nwt;
		ii20= (ii19+ ii01) % nwt;
		ii21= (ii20+ ii01) % nwt;
		ii22= (ii21+ ii01) % nwt;
		ii23= (ii22+ ii01) % nwt;
		ii24= (ii23+ ii01) % nwt;
		ii25= (ii24+ ii01) % nwt;
		ii26= (ii25+ ii01) % nwt;
		ii27= (ii26+ ii01) % nwt;
#endif
	}

	if(CY_THREADS > 1)
	{
		for(ithread = 1; ithread < CY_THREADS; ithread++)
		{
			_i[ithread] = ((uint32)(sw - _bjmodnini[ithread]) >> 31);
		}
	}

	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		_bjmodn00[ithread] = _bjmodnini[ithread];
		_bjmodn01[ithread] = _bjmodn00[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn01[ithread] = _bjmodn01[ithread] + ( (-(int)((uint32)_bjmodn01[ithread] >> 31)) & n);
		_bjmodn02[ithread] = _bjmodn01[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn02[ithread] = _bjmodn02[ithread] + ( (-(int)((uint32)_bjmodn02[ithread] >> 31)) & n);
		_bjmodn03[ithread] = _bjmodn02[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn03[ithread] = _bjmodn03[ithread] + ( (-(int)((uint32)_bjmodn03[ithread] >> 31)) & n);
		_bjmodn04[ithread] = _bjmodn03[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn04[ithread] = _bjmodn04[ithread] + ( (-(int)((uint32)_bjmodn04[ithread] >> 31)) & n);
		_bjmodn05[ithread] = _bjmodn04[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn05[ithread] = _bjmodn05[ithread] + ( (-(int)((uint32)_bjmodn05[ithread] >> 31)) & n);
		_bjmodn06[ithread] = _bjmodn05[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn06[ithread] = _bjmodn06[ithread] + ( (-(int)((uint32)_bjmodn06[ithread] >> 31)) & n);
		_bjmodn07[ithread] = _bjmodn06[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn07[ithread] = _bjmodn07[ithread] + ( (-(int)((uint32)_bjmodn07[ithread] >> 31)) & n);
		_bjmodn08[ithread] = _bjmodn07[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn08[ithread] = _bjmodn08[ithread] + ( (-(int)((uint32)_bjmodn08[ithread] >> 31)) & n);
		_bjmodn09[ithread] = _bjmodn08[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn09[ithread] = _bjmodn09[ithread] + ( (-(int)((uint32)_bjmodn09[ithread] >> 31)) & n);
		_bjmodn10[ithread] = _bjmodn09[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn10[ithread] = _bjmodn10[ithread] + ( (-(int)((uint32)_bjmodn10[ithread] >> 31)) & n);
		_bjmodn11[ithread] = _bjmodn10[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn11[ithread] = _bjmodn11[ithread] + ( (-(int)((uint32)_bjmodn11[ithread] >> 31)) & n);
		_bjmodn12[ithread] = _bjmodn11[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn12[ithread] = _bjmodn12[ithread] + ( (-(int)((uint32)_bjmodn12[ithread] >> 31)) & n);
		_bjmodn13[ithread] = _bjmodn12[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn13[ithread] = _bjmodn13[ithread] + ( (-(int)((uint32)_bjmodn13[ithread] >> 31)) & n);
		_bjmodn14[ithread] = _bjmodn13[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn14[ithread] = _bjmodn14[ithread] + ( (-(int)((uint32)_bjmodn14[ithread] >> 31)) & n);
		_bjmodn15[ithread] = _bjmodn14[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn15[ithread] = _bjmodn15[ithread] + ( (-(int)((uint32)_bjmodn15[ithread] >> 31)) & n);
		_bjmodn16[ithread] = _bjmodn15[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn16[ithread] = _bjmodn16[ithread] + ( (-(int)((uint32)_bjmodn16[ithread] >> 31)) & n);
		_bjmodn17[ithread] = _bjmodn16[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn17[ithread] = _bjmodn17[ithread] + ( (-(int)((uint32)_bjmodn17[ithread] >> 31)) & n);
		_bjmodn18[ithread] = _bjmodn17[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn18[ithread] = _bjmodn18[ithread] + ( (-(int)((uint32)_bjmodn18[ithread] >> 31)) & n);
		_bjmodn19[ithread] = _bjmodn18[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn19[ithread] = _bjmodn19[ithread] + ( (-(int)((uint32)_bjmodn19[ithread] >> 31)) & n);
		_bjmodn20[ithread] = _bjmodn19[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn20[ithread] = _bjmodn20[ithread] + ( (-(int)((uint32)_bjmodn20[ithread] >> 31)) & n);
		_bjmodn21[ithread] = _bjmodn20[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn21[ithread] = _bjmodn21[ithread] + ( (-(int)((uint32)_bjmodn21[ithread] >> 31)) & n);
		_bjmodn22[ithread] = _bjmodn21[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn22[ithread] = _bjmodn22[ithread] + ( (-(int)((uint32)_bjmodn22[ithread] >> 31)) & n);
		_bjmodn23[ithread] = _bjmodn22[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn23[ithread] = _bjmodn23[ithread] + ( (-(int)((uint32)_bjmodn23[ithread] >> 31)) & n);
		_bjmodn24[ithread] = _bjmodn23[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn24[ithread] = _bjmodn24[ithread] + ( (-(int)((uint32)_bjmodn24[ithread] >> 31)) & n);
		_bjmodn25[ithread] = _bjmodn24[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn25[ithread] = _bjmodn25[ithread] + ( (-(int)((uint32)_bjmodn25[ithread] >> 31)) & n);
		_bjmodn26[ithread] = _bjmodn25[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn26[ithread] = _bjmodn26[ithread] + ( (-(int)((uint32)_bjmodn26[ithread] >> 31)) & n);
		_bjmodn27[ithread] = _bjmodn26[ithread] + _bjmodnini[CY_THREADS] - n; _bjmodn27[ithread] = _bjmodn27[ithread] + ( (-(int)((uint32)_bjmodn27[ithread] >> 31)) & n);

		if(MODULUS_TYPE == MODULUS_TYPE_FERMAT)
		{
			/* Start this value off at N in Fermat-mod case, so (bjmodn >= sw) check in
			fermat_carry_norm_errcheck (cf. carry.h) yields a bigword (i == 1) for j= 0:
			*/
			_bjmodn00[ithread] = n;
			_bjmodn07[ithread] = n;
			_bjmodn14[ithread] = n;
			_bjmodn21[ithread] = n;

		}
	}

	/* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
	if(!full_pass)
	{
		khi = 1;
	#if FFT_DEBUG
		fprintf(dbg_file, "radix28_ditN_cy_dif1: Cleanup Pass:\n", n);
	#endif
	}

/* Needed to remove the prefetch-address vars add0 & add for this to compile properly: */
#ifdef MULTITHREAD
	omp_set_num_threads(CY_THREADS);
//#undef PFETCH
	#pragma omp parallel for private(temp,frac,maxerr,i,j,j1,jstart,jhi,k,k1,k2,l,col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,wtl,wtlp1,wtn,wtnm1,wt,wtinv,wtA,wtB,wtC,wt_re,wt_im,rt,it,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i,bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27,cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27) default(shared) schedule(static)
#endif

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

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
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
		#endif
		}

	#ifdef USE_SSE2
		/* init carries	*/
		cy_r00->re = _cy_r00[ithread];	cy_i00->re = _cy_i00[ithread];
		cy_r00->im = _cy_r01[ithread];	cy_i00->im = _cy_i01[ithread];
		cy_r02->re = _cy_r02[ithread];	cy_i02->re = _cy_i02[ithread];
		cy_r02->im = _cy_r03[ithread];	cy_i02->im = _cy_i03[ithread];
		cy_r04->re = _cy_r04[ithread];	cy_i04->re = _cy_i04[ithread];
		cy_r04->im = _cy_r05[ithread];	cy_i04->im = _cy_i05[ithread];
		cy_r06->re = _cy_r06[ithread];	cy_i06->re = _cy_i06[ithread];
		cy_r06->im = _cy_r07[ithread];	cy_i06->im = _cy_i07[ithread];
		cy_r08->re = _cy_r08[ithread];	cy_i08->re = _cy_i08[ithread];
		cy_r08->im = _cy_r09[ithread];	cy_i08->im = _cy_i09[ithread];
		cy_r10->re = _cy_r10[ithread];	cy_i10->re = _cy_i10[ithread];
		cy_r10->im = _cy_r11[ithread];	cy_i10->im = _cy_i11[ithread];
		cy_r12->re = _cy_r12[ithread];	cy_i12->re = _cy_i12[ithread];
		cy_r12->im = _cy_r13[ithread];	cy_i12->im = _cy_i13[ithread];
		cy_r14->re = _cy_r14[ithread];	cy_i14->re = _cy_i14[ithread];
		cy_r14->im = _cy_r15[ithread];	cy_i14->im = _cy_i15[ithread];
		cy_r16->re = _cy_r16[ithread];	cy_i16->re = _cy_i16[ithread];
		cy_r16->im = _cy_r17[ithread];	cy_i16->im = _cy_i17[ithread];
		cy_r18->re = _cy_r18[ithread];	cy_i18->re = _cy_i18[ithread];
		cy_r18->im = _cy_r19[ithread];	cy_i18->im = _cy_i19[ithread];
		cy_r20->re = _cy_r20[ithread];	cy_i20->re = _cy_i20[ithread];
		cy_r20->im = _cy_r21[ithread];	cy_i20->im = _cy_i21[ithread];
		cy_r22->re = _cy_r22[ithread];	cy_i22->re = _cy_i22[ithread];
		cy_r22->im = _cy_r23[ithread];	cy_i22->im = _cy_i23[ithread];
		cy_r24->re = _cy_r24[ithread];	cy_i24->re = _cy_i24[ithread];
		cy_r24->im = _cy_r25[ithread];	cy_i24->im = _cy_i25[ithread];
		cy_r26->re = _cy_r26[ithread];	cy_i26->re = _cy_i26[ithread];
		cy_r26->im = _cy_r27[ithread];	cy_i26->im = _cy_i27[ithread];
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

		#ifdef DEBUG_SSE2
			rng_isaac_init(TRUE);
			jt = j1;		jp = j2;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			fprintf(stderr, "radix28_wrapper: A_in[00] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_wrapper: A_in[01] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_wrapper: A_in[02] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_wrapper: A_in[03] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);
			jt = j1+p04;	jp = j2+p04;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			fprintf(stderr, "radix28_wrapper: A_in[04] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_wrapper: A_in[05] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_wrapper: A_in[06] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_wrapper: A_in[07] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);
			jt = j1+p08;	jp = j2+p08;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			fprintf(stderr, "radix28_wrapper: A_in[08] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_wrapper: A_in[09] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_wrapper: A_in[10] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_wrapper: A_in[11] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);
			jt = j1+p12;	jp = j2+p12;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			fprintf(stderr, "radix28_wrapper: A_in[12] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_wrapper: A_in[13] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_wrapper: A_in[14] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_wrapper: A_in[15] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);
			jt = j1+p16;	jp = j2+p16;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			fprintf(stderr, "radix28_wrapper: A_in[16] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_wrapper: A_in[17] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_wrapper: A_in[18] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_wrapper: A_in[29] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);
			jt = j1+p20;	jp = j2+p20;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			fprintf(stderr, "radix28_wrapper: A_in[20] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_wrapper: A_in[21] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_wrapper: A_in[22] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_wrapper: A_in[23] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);
			jt = j1+p24;	jp = j2+p24;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			fprintf(stderr, "radix28_wrapper: A_in[24] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_wrapper: A_in[25] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_wrapper: A_in[26] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_wrapper: A_in[27] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);
			fprintf(stderr, "\n");
		#endif

		/*...The radix-28 DIT pass is here:	*/

		/* EWM: 10/18/04: We swap the odd-index outputs of each of the radix-4 DIT transforms (1<=>3, 5<=>7, etc.) so that the indexing
						  of all the radix-7 transforms (really just the 2nd and 4th of these) winds up being in-place. This allows us
						  to properly re-use the ajp1 variables in the carry-pass version of this routine.
		*/
		#ifdef USE_SSE2

		  #if defined(COMPILER_TYPE_MSVC)

			/* Since doing radix-7 in-place here, outputs of radix-4 are in consecutive memory locs, i.e. 0x20 bytes apart, e.g. the distance between s1p00r and s1p01r: */

			add0 = &a[j1    ]; 	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p00r,s1p03r,s1p02r,s1p01r)
			add3 = &a[j1+p12];	add2 = add3+p01;	add1 = add3+p02;	add0 = add3+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p04r,s1p07r,s1p06r,s1p05r)
			add1 = &a[j1+p24];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p08r,s1p11r,s1p10r,s1p09r)
			add1 = &a[j1+p08];	add0 = add1+p01;	add2 = add1+p02;	add3 = add1+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p12r,s1p15r,s1p14r,s1p13r)
			add2 = &a[j1+p20];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p16r,s1p19r,s1p18r,s1p17r)
			add2 = &a[j1+p04];	add3 = add2+p01;	add0 = add2+p02;	add1 = add2+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p20r,s1p23r,s1p22r,s1p21r)
			add0 = &a[j1+p16];	add1 = add0+p01;	add3 = add0+p02;	add2 = add0+p03;	SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_E(add0,add1,add2,add3, s1p24r,s1p27r,s1p26r,s1p25r)

			/*...and now do 4 radix-7 transforms...*/

			SSE2_RADIX_07_DFT(s1p00r,s1p04r,s1p08r,s1p12r,s1p16r,s1p20r,s1p24r,cc0,s1p00r,s1p08r,s1p16r,s1p24r,s1p04r,s1p12r,s1p20r);
			SSE2_RADIX_07_DFT(s1p03r,s1p07r,s1p11r,s1p15r,s1p19r,s1p23r,s1p27r,cc0,s1p07r,s1p15r,s1p23r,s1p03r,s1p11r,s1p19r,s1p27r);
			SSE2_RADIX_07_DFT(s1p02r,s1p06r,s1p10r,s1p14r,s1p18r,s1p22r,s1p26r,cc0,s1p14r,s1p22r,s1p02r,s1p10r,s1p18r,s1p26r,s1p06r);
			SSE2_RADIX_07_DFT(s1p01r,s1p05r,s1p09r,s1p13r,s1p17r,s1p21r,s1p25r,cc0,s1p21r,s1p01r,s1p09r,s1p17r,s1p25r,s1p05r,s1p13r);

		  #else	/* GCC-style inline ASM: */

			add0 = &a[j1    ];
			SSE2_RADIX28_DIT_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p12,p16,p20,p24,s1p00r,cc0);

		  #endif

		#else	/* USE_SSE2 */

		  /*...gather the needed data (28 64-bit complex, i.e. 56 64-bit reals) and do 7 radix-4 transforms...*/
							 /*                                      outputs                                      */ /*                          inputs                           */
			RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a1p00r,a1p00i,a1p03r,a1p03i,a1p02r,a1p02i,a1p01r,a1p01i,rt,it);	jt = j1+p12; jp = j2+p12;
			RADIX_04_DIT(a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a1p04r,a1p04i,a1p07r,a1p07i,a1p06r,a1p06i,a1p05r,a1p05i,rt,it);	jt = j1+p24; jp = j2+p24;
			RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a1p08r,a1p08i,a1p11r,a1p11i,a1p10r,a1p10i,a1p09r,a1p09i,rt,it);	jt = j1+p08; jp = j2+p08;
			RADIX_04_DIT(a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a1p12r,a1p12i,a1p15r,a1p15i,a1p14r,a1p14i,a1p13r,a1p13i,rt,it);	jt = j1+p20; jp = j2+p20;
			RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a1p16r,a1p16i,a1p19r,a1p19i,a1p18r,a1p18i,a1p17r,a1p17i,rt,it);	jt = j1+p04; jp = j2+p04;
			RADIX_04_DIT(a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a1p20r,a1p20i,a1p23r,a1p23i,a1p22r,a1p22i,a1p21r,a1p21i,rt,it);	jt = j1+p16; jp = j2+p16;
			RADIX_04_DIT(a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a1p24r,a1p24i,a1p27r,a1p27i,a1p26r,a1p26i,a1p25r,a1p25i,rt,it);

		  /*...and now do 4 radix-7 transforms...*/
		  #if LO_ADD
							 /*                                                   inputs                                                  */ /*               intermediates              */ /*                                                                     outputs                                                                         */
			RADIX_07_DFT(a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,a1p24r,a1p24i,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
			RADIX_07_DFT(a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,a1p27r,a1p27i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
			RADIX_07_DFT(a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p14r,a1p14i,a1p22r,a1p22i,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,a1p26r,a1p26i,a1p06r,a1p06i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
			RADIX_07_DFT(a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p21r,a1p21i,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,a1p25r,a1p25i,a1p05r,a1p05i,a1p13r,a1p13i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		  #else
			RADIX_07_DFT_NUSS(a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,a1p24r,a1p24i,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3,rt,it);
			RADIX_07_DFT_NUSS(a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,a1p27r,a1p27i,cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3,rt,it);
			RADIX_07_DFT_NUSS(a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p14r,a1p14i,a1p22r,a1p22i,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,a1p26r,a1p26i,a1p06r,a1p06i,cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3,rt,it);
			RADIX_07_DFT_NUSS(a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p21r,a1p21i,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,a1p25r,a1p25i,a1p05r,a1p05i,a1p13r,a1p13i,cx0,sx0,cx1,sx1,cx2,sx2,cx3,sx3,rt,it);
		  #endif

		#endif	/* USE_SSE2 */

		/*...Now do the carries. Since the outputs would
			normally be getting dispatched to 28 separate blocks of the A-array, we need 28 separate carries.	*/

//		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
//		{
#if 1
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

				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24);

		  #else	/* GCC-style inline ASM: */

				SSE2_cmplx_carry_norm_errcheck0_2B(s1p00r,add1,add2,add3,cy_r00,cy_r02,bjmodn00,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p04r,add1,add2,add3,cy_r04,cy_r06,bjmodn04,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p08r,add1,add2,add3,cy_r08,cy_r10,bjmodn08,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p12r,add1,add2,add3,cy_r12,cy_r14,bjmodn12,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p16r,add1,add2,add3,cy_r16,cy_r18,bjmodn16,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p20r,add1,add2,add3,cy_r20,cy_r22,bjmodn20,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck1_2B (s1p24r,add1,add2,add3,cy_r24,cy_r26,bjmodn24,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);

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

				SSE2_cmplx_carry_norm_nocheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24);

		  #else	/* GCC-style inline ASM: */

				SSE2_cmplx_carry_norm_nocheck2_2B(s1p00r,add1,add2,     cy_r00,cy_r02,bjmodn00,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p04r,add1,add2,     cy_r04,cy_r06,bjmodn04,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p08r,add1,add2,     cy_r08,cy_r10,bjmodn08,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p12r,add1,add2,     cy_r12,cy_r14,bjmodn12,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p16r,add1,add2,     cy_r16,cy_r18,bjmodn16,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p20r,add1,add2,     cy_r20,cy_r22,bjmodn20,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);
				SSE2_cmplx_carry_norm_nocheck2_2B(s1p24r,add1,add2,     cy_r24,cy_r26,bjmodn24,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_n,sse_sw);

		  #endif

				i =((uint32)(sw - *bjmodn00) >> 31);	/* get ready for the next set...	*/

		#else	/* USE_SSE2 */

			/*...set0 is slightly different from others:	*/
		   cmplx_carry_norm_errcheck0(a1p00r,a1p00i,cy_r00,bjmodn00   );
			cmplx_carry_norm_errcheck(a1p01r,a1p01i,cy_r01,bjmodn01,1 );
			cmplx_carry_norm_errcheck(a1p02r,a1p02i,cy_r02,bjmodn02,2 );
			cmplx_carry_norm_errcheck(a1p03r,a1p03i,cy_r03,bjmodn03,3 );
			cmplx_carry_norm_errcheck(a1p04r,a1p04i,cy_r04,bjmodn04,4 );
			cmplx_carry_norm_errcheck(a1p05r,a1p05i,cy_r05,bjmodn05,5 );
			cmplx_carry_norm_errcheck(a1p06r,a1p06i,cy_r06,bjmodn06,6 );
			cmplx_carry_norm_errcheck(a1p07r,a1p07i,cy_r07,bjmodn07,7 );
			cmplx_carry_norm_errcheck(a1p08r,a1p08i,cy_r08,bjmodn08,8 );
			cmplx_carry_norm_errcheck(a1p09r,a1p09i,cy_r09,bjmodn09,9 );
			cmplx_carry_norm_errcheck(a1p10r,a1p10i,cy_r10,bjmodn10,10);
			cmplx_carry_norm_errcheck(a1p11r,a1p11i,cy_r11,bjmodn11,11);
			cmplx_carry_norm_errcheck(a1p12r,a1p12i,cy_r12,bjmodn12,12);
			cmplx_carry_norm_errcheck(a1p13r,a1p13i,cy_r13,bjmodn13,13);
			cmplx_carry_norm_errcheck(a1p14r,a1p14i,cy_r14,bjmodn14,14);
			cmplx_carry_norm_errcheck(a1p15r,a1p15i,cy_r15,bjmodn15,15);
			cmplx_carry_norm_errcheck(a1p16r,a1p16i,cy_r16,bjmodn16,16);
			cmplx_carry_norm_errcheck(a1p17r,a1p17i,cy_r17,bjmodn17,17);
			cmplx_carry_norm_errcheck(a1p18r,a1p18i,cy_r18,bjmodn18,18);
			cmplx_carry_norm_errcheck(a1p19r,a1p19i,cy_r19,bjmodn19,19);
			cmplx_carry_norm_errcheck(a1p20r,a1p20i,cy_r20,bjmodn20,20);
			cmplx_carry_norm_errcheck(a1p21r,a1p21i,cy_r21,bjmodn21,21);
			cmplx_carry_norm_errcheck(a1p22r,a1p22i,cy_r22,bjmodn22,22);
			cmplx_carry_norm_errcheck(a1p23r,a1p23i,cy_r23,bjmodn23,23);
			cmplx_carry_norm_errcheck(a1p24r,a1p24i,cy_r24,bjmodn24,24);
			cmplx_carry_norm_errcheck(a1p25r,a1p25i,cy_r25,bjmodn25,25);
			cmplx_carry_norm_errcheck(a1p26r,a1p26i,cy_r26,bjmodn26,26);
			cmplx_carry_norm_errcheck(a1p27r,a1p27i,cy_r27,bjmodn27,27);

			i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
						and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

		#endif	/* #ifdef USE_SSE2 */
#elif 0
			}
			else
			{
				ASSERT(HERE, 0, "Fermat-mod carries not yet supported for SSE2!");
		}
		else
		{
			fermat_carry_norm_errcheck(a1p00r,a1p00i,cy_r00,cy_i00,ii00,bjmodn00,0 *n28);
			fermat_carry_norm_errcheck(a1p01r,a1p01i,cy_r01,cy_i01,ii01,bjmodn01,1 *n28);
			fermat_carry_norm_errcheck(a1p02r,a1p02i,cy_r02,cy_i02,ii02,bjmodn02,2 *n28);
			fermat_carry_norm_errcheck(a1p03r,a1p03i,cy_r03,cy_i03,ii03,bjmodn03,3 *n28);
			fermat_carry_norm_errcheck(a1p04r,a1p04i,cy_r04,cy_i04,ii04,bjmodn04,4 *n28);
			fermat_carry_norm_errcheck(a1p05r,a1p05i,cy_r05,cy_i05,ii05,bjmodn05,5 *n28);
			fermat_carry_norm_errcheck(a1p06r,a1p06i,cy_r06,cy_i06,ii06,bjmodn06,6 *n28);
			fermat_carry_norm_errcheck(a1p07r,a1p07i,cy_r07,cy_i07,ii07,bjmodn07,7 *n28);
			fermat_carry_norm_errcheck(a1p08r,a1p08i,cy_r08,cy_i08,ii08,bjmodn08,8 *n28);
			fermat_carry_norm_errcheck(a1p09r,a1p09i,cy_r09,cy_i09,ii09,bjmodn09,9 *n28);
			fermat_carry_norm_errcheck(a1p10r,a1p10i,cy_r10,cy_i10,ii10,bjmodn10,10*n28);
			fermat_carry_norm_errcheck(a1p11r,a1p11i,cy_r11,cy_i11,ii11,bjmodn11,11*n28);
			fermat_carry_norm_errcheck(a1p12r,a1p12i,cy_r12,cy_i12,ii12,bjmodn12,12*n28);
			fermat_carry_norm_errcheck(a1p13r,a1p13i,cy_r13,cy_i13,ii13,bjmodn13,13*n28);
			fermat_carry_norm_errcheck(a1p14r,a1p14i,cy_r14,cy_i14,ii14,bjmodn14,14*n28);
			fermat_carry_norm_errcheck(a1p15r,a1p15i,cy_r15,cy_i15,ii15,bjmodn15,15*n28);
			fermat_carry_norm_errcheck(a1p16r,a1p16i,cy_r16,cy_i16,ii16,bjmodn16,16*n28);
			fermat_carry_norm_errcheck(a1p17r,a1p17i,cy_r17,cy_i17,ii17,bjmodn17,17*n28);
			fermat_carry_norm_errcheck(a1p18r,a1p18i,cy_r18,cy_i18,ii18,bjmodn18,18*n28);
			fermat_carry_norm_errcheck(a1p19r,a1p19i,cy_r19,cy_i19,ii19,bjmodn19,19*n28);
			fermat_carry_norm_errcheck(a1p20r,a1p20i,cy_r20,cy_i20,ii20,bjmodn20,20*n28);
			fermat_carry_norm_errcheck(a1p21r,a1p21i,cy_r21,cy_i21,ii21,bjmodn21,21*n28);
			fermat_carry_norm_errcheck(a1p22r,a1p22i,cy_r22,cy_i22,ii22,bjmodn22,22*n28);
			fermat_carry_norm_errcheck(a1p23r,a1p23i,cy_r23,cy_i23,ii23,bjmodn23,23*n28);
			fermat_carry_norm_errcheck(a1p24r,a1p24i,cy_r24,cy_i24,ii24,bjmodn24,24*n28);
			fermat_carry_norm_errcheck(a1p25r,a1p25i,cy_r25,cy_i25,ii25,bjmodn25,25*n28);
			fermat_carry_norm_errcheck(a1p26r,a1p26i,cy_r26,cy_i26,ii26,bjmodn26,26*n28);
			fermat_carry_norm_errcheck(a1p27r,a1p27i,cy_r27,cy_i27,ii27,bjmodn27,27*n28);
		}
#endif

/*...The radix-28 DIF pass is here:	*/

		#ifdef USE_SSE2

		  #if defined(COMPILER_TYPE_MSVC)

			SSE2_RADIX_07_DFT(s1p00r,s1p24r,s1p20r,s1p16r,s1p12r,s1p08r,s1p04r,cc0,s1p00r,s1p04r,s1p08r,s1p12r,s1p16r,s1p20r,s1p24r);
			SSE2_RADIX_07_DFT(s1p21r,s1p17r,s1p13r,s1p09r,s1p05r,s1p01r,s1p25r,cc0,s1p01r,s1p05r,s1p09r,s1p13r,s1p17r,s1p21r,s1p25r);
			SSE2_RADIX_07_DFT(s1p14r,s1p10r,s1p06r,s1p02r,s1p26r,s1p22r,s1p18r,cc0,s1p02r,s1p06r,s1p10r,s1p14r,s1p18r,s1p22r,s1p26r);
			SSE2_RADIX_07_DFT(s1p07r,s1p03r,s1p27r,s1p23r,s1p19r,s1p15r,s1p11r,cc0,s1p03r,s1p07r,s1p11r,s1p15r,s1p19r,s1p23r,s1p27r);

			/* Since doing radix-7 in-place here, inputs of radix-4 are in consecutive memory locs, i.e. 0x20 bytes apart, e.g. the distance between s1p00r and s1p01r: */

			add0 = &a[j1    ]; 	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p00r,s1p01r,s1p02r,s1p03r, add0,add1,add2,add3)
			add1 = &a[j1+p24];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p04r,s1p05r,s1p06r,s1p07r, add0,add1,add2,add3)
			add3 = &a[j1+p20];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p08r,s1p09r,s1p10r,s1p11r, add0,add1,add2,add3)
			add0 = &a[j1+p16];	add1 = add0+p01;	add2 = add0+p02;	add3 = add0+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p12r,s1p13r,s1p14r,s1p15r, add0,add1,add2,add3)
			add2 = &a[j1+p12];	add3 = add2+p01;	add1 = add2+p02;	add0 = add2+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p16r,s1p17r,s1p18r,s1p19r, add0,add1,add2,add3)
			add1 = &a[j1+p08];	add0 = add1+p01;	add3 = add1+p02;	add2 = add1+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p20r,s1p21r,s1p22r,s1p23r, add0,add1,add2,add3)
			add3 = &a[j1+p04];	add2 = add3+p01;	add0 = add3+p02;	add1 = add3+p03;	SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_E(s1p24r,s1p25r,s1p26r,s1p27r, add0,add1,add2,add3)

		  #else	/* GCC-style inline ASM: */

			add0 = &a[j1    ];
			SSE2_RADIX28_DIF_NOTWIDDLE(add0,p01,p02,p03,p04,p08,p12,p16,p20,p24,s1p00r,cc0);

		  #endif

		#if 0//DEBUG_SSE2
			fprintf(stderr, "radix28_carry_out: s00= %20.5f, %20.5f\n",s1p00r->re,s1p00i->re);
			fprintf(stderr, "radix28_carry_out: s01= %20.5f, %20.5f\n",s1p01r->re,s1p01i->re);
			fprintf(stderr, "radix28_carry_out: s02= %20.5f, %20.5f\n",s1p02r->re,s1p02i->re);
			fprintf(stderr, "radix28_carry_out: s03= %20.5f, %20.5f\n",s1p03r->re,s1p03i->re);
			fprintf(stderr, "radix28_carry_out: s04= %20.5f, %20.5f\n",s1p04r->re,s1p04i->re);
			fprintf(stderr, "radix28_carry_out: s05= %20.5f, %20.5f\n",s1p05r->re,s1p05i->re);
			fprintf(stderr, "radix28_carry_out: s06= %20.5f, %20.5f\n",s1p06r->re,s1p06i->re);
			fprintf(stderr, "radix28_carry_out: s07= %20.5f, %20.5f\n",s1p07r->re,s1p07i->re);
			fprintf(stderr, "radix28_carry_out: s08= %20.5f, %20.5f\n",s1p08r->re,s1p08i->re);
			fprintf(stderr, "radix28_carry_out: s09= %20.5f, %20.5f\n",s1p09r->re,s1p09i->re);
			fprintf(stderr, "radix28_carry_out: s10= %20.5f, %20.5f\n",s1p10r->re,s1p10i->re);
			fprintf(stderr, "radix28_carry_out: s11= %20.5f, %20.5f\n",s1p11r->re,s1p11i->re);
			fprintf(stderr, "radix28_carry_out: s12= %20.5f, %20.5f\n",s1p12r->re,s1p12i->re);
			fprintf(stderr, "radix28_carry_out: s13= %20.5f, %20.5f\n",s1p13r->re,s1p13i->re);
			fprintf(stderr, "radix28_carry_out: s14= %20.5f, %20.5f\n",s1p14r->re,s1p14i->re);
			fprintf(stderr, "radix28_carry_out: s15= %20.5f, %20.5f\n",s1p15r->re,s1p15i->re);
			fprintf(stderr, "radix28_carry_out: s16= %20.5f, %20.5f\n",s1p16r->re,s1p16i->re);
			fprintf(stderr, "radix28_carry_out: s17= %20.5f, %20.5f\n",s1p17r->re,s1p17i->re);
			fprintf(stderr, "radix28_carry_out: s18= %20.5f, %20.5f\n",s1p18r->re,s1p18i->re);
			fprintf(stderr, "radix28_carry_out: s19= %20.5f, %20.5f\n",s1p19r->re,s1p19i->re);
			fprintf(stderr, "radix28_carry_out: s20= %20.5f, %20.5f\n",s1p20r->re,s1p20i->re);
			fprintf(stderr, "radix28_carry_out: s21= %20.5f, %20.5f\n",s1p21r->re,s1p21i->re);
			fprintf(stderr, "radix28_carry_out: s22= %20.5f, %20.5f\n",s1p22r->re,s1p22i->re);
			fprintf(stderr, "radix28_carry_out: s23= %20.5f, %20.5f\n",s1p23r->re,s1p23i->re);
			fprintf(stderr, "radix28_carry_out: s24= %20.5f, %20.5f\n",s1p24r->re,s1p24i->re);
			fprintf(stderr, "radix28_carry_out: s25= %20.5f, %20.5f\n",s1p25r->re,s1p25i->re);
			fprintf(stderr, "radix28_carry_out: s26= %20.5f, %20.5f\n",s1p26r->re,s1p26i->re);
			fprintf(stderr, "radix28_carry_out: s27= %20.5f, %20.5f\n",s1p27r->re,s1p27i->re);
			exit(0);
		#endif

		#ifdef DEBUG_SSE2
			jt = j1;		jp = j2;
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[00] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[01] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[02] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[03] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);			jt = j1+p04;	jp = j2+p04;
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[04] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[05] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[06] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[07] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);			jt = j1+p08;	jp = j2+p08;
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[08] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[09] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[10] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[11] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);			jt = j1+p12;	jp = j2+p12;
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[12] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[13] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[14] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[15] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);			jt = j1+p16;	jp = j2+p16;
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[16] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[17] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[18] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[29] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);			jt = j1+p20;	jp = j2+p20;
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[20] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[21] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[22] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[23] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);			jt = j1+p24;	jp = j2+p24;
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[24] = %24.5f, %24.5f\n",a[jt    ],a[jp    ]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[25] = %24.5f, %24.5f\n",a[jt+p01],a[jp+p01]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[26] = %24.5f, %24.5f\n",a[jt+p02],a[jp+p02]);
			fprintf(stderr, "radix28_ditN_cy_dif1: A_out[27] = %24.5f, %24.5f\n",a[jt+p03],a[jp+p03]);
			exit(0);
		#endif

		#else	/* !USE_SSE2 */

		  #if PFETCH
			addr = &a[j1];
			prefetch_p_doubles(addr);
		  #endif

		/*...gather the needed data (28 64-bit complex, i.e. 56 64-bit reals) and do 4 radix-7 transforms...*/
							 /*                                                                      inputs                                                                         */ /*               intermediates              */ /*                                                  outputs                                                  */
		  #if PFETCH
			RADIX_07_DFT_PFETCH(a1p00r,a1p00i,a1p24r,a1p24i,a1p20r,a1p20i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im,addr,addp,p01,p02,p03);	jt=p04+p01;	jp=p04+p02;
			RADIX_07_DFT_PFETCH(a1p21r,a1p21i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,a1p05r,a1p05i,a1p01r,a1p01i,a1p25r,a1p25i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im,addr,addp,p04, jt, jp);	jt=p08-p01;	jp=p08+p01;
			RADIX_07_DFT_PFETCH(a1p14r,a1p14i,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p26r,a1p26i,a1p22r,a1p22i,a1p18r,a1p18i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im,addr,addp, jt,p08, jp);	jt=p08+p02;	jp=p08+p03;
			RADIX_07_DFT_PFETCH(a1p07r,a1p07i,a1p03r,a1p03i,a1p27r,a1p27i,a1p23r,a1p23i,a1p19r,a1p19i,a1p15r,a1p15i,a1p11r,a1p11i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im,addr,addp, jt, jp,p12);
		  #else
			RADIX_07_DFT       (a1p00r,a1p00i,a1p24r,a1p24i,a1p20r,a1p20i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
			RADIX_07_DFT       (a1p21r,a1p21i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,a1p05r,a1p05i,a1p01r,a1p01i,a1p25r,a1p25i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
			RADIX_07_DFT       (a1p14r,a1p14i,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p26r,a1p26i,a1p22r,a1p22i,a1p18r,a1p18i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
			RADIX_07_DFT       (a1p07r,a1p07i,a1p03r,a1p03i,a1p27r,a1p27i,a1p23r,a1p23i,a1p19r,a1p19i,a1p15r,a1p15i,a1p11r,a1p11i,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		  #endif

		/*...and now do 7 radix-4 transforms...*/
							 /*                          inputs                           */ /*                                      outputs                                      */
		  #if PFETCH
			addp = addr+p12+p01;
			prefetch_p_doubles(addp);

			RADIX_04_DIF_PFETCH(a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],rt,it,addr,addp,p12+p02,p12+p03);	jt = j1+p24; jp = j2+p24;
			RADIX_04_DIF_PFETCH(a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it,addr,addp,p16    ,p16+p01);	jt = j1+p20; jp = j2+p20;
			RADIX_04_DIF_PFETCH(a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it,addr,addp,p16+p02,p16+p03);	jt = j1+p16; jp = j2+p16;
			RADIX_04_DIF_PFETCH(a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it,addr,addp,p20    ,p20+p01);	jt = j1+p12; jp = j2+p12;
			RADIX_04_DIF_PFETCH(a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it,addr,addp,p20+p02,p20+p03);	jt = j1+p08; jp = j2+p08;
			RADIX_04_DIF_PFETCH(a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it,addr,addp,p24    ,p20+p01);	jt = j1+p04; jp = j2+p04;
			RADIX_04_DIF_PFETCH(a1p24r,a1p24i,a1p25r,a1p25i,a1p26r,a1p26i,a1p27r,a1p27i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it,addr,addp,p20+p02,p20+p03);
		  #else
			RADIX_04_DIF       (a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],rt,it);	jt = j1+p24; jp = j2+p24;
			RADIX_04_DIF       (a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p20; jp = j2+p20;
			RADIX_04_DIF       (a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);	jt = j1+p16; jp = j2+p16;
			RADIX_04_DIF       (a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],rt,it);	jt = j1+p12; jp = j2+p12;
			RADIX_04_DIF       (a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],a[jt    ],a[jp    ],a[jt+p01],a[jp+p01],rt,it);	jt = j1+p08; jp = j2+p08;
			RADIX_04_DIF       (a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],a[jt+p03],a[jp+p03],a[jt+p02],a[jp+p02],rt,it);	jt = j1+p04; jp = j2+p04;
			RADIX_04_DIF       (a1p24r,a1p24i,a1p25r,a1p25i,a1p26r,a1p26i,a1p27r,a1p27i,a[jt+p02],a[jp+p02],a[jt+p03],a[jp+p03],a[jt+p01],a[jp+p01],a[jt    ],a[jp    ],rt,it);
		  #endif

		#endif	/* !USE_SSE2 */

			}

			if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
			{
				jstart += nwt;
				jhi    += nwt;

				col += 28;
				co3 -= 28;
			}
		}	/* end for(k=1; k <= khi; k++) */

		/* At end of each thread-processed work chunk, dump the
		carryouts into their non-thread-private array slots:
		*/
	#ifdef USE_SSE2
		_cy_r00[ithread] = cy_r00->re;	_cy_i00[ithread] = cy_i00->re;
		_cy_r01[ithread] = cy_r00->im;	_cy_i01[ithread] = cy_i00->im;
		_cy_r02[ithread] = cy_r02->re;	_cy_i02[ithread] = cy_i02->re;
		_cy_r03[ithread] = cy_r02->im;	_cy_i03[ithread] = cy_i02->im;
		_cy_r04[ithread] = cy_r04->re;	_cy_i04[ithread] = cy_i04->re;
		_cy_r05[ithread] = cy_r04->im;	_cy_i05[ithread] = cy_i04->im;
		_cy_r06[ithread] = cy_r06->re;	_cy_i06[ithread] = cy_i06->re;
		_cy_r07[ithread] = cy_r06->im;	_cy_i07[ithread] = cy_i06->im;
		_cy_r08[ithread] = cy_r08->re;	_cy_i08[ithread] = cy_i08->re;
		_cy_r09[ithread] = cy_r08->im;	_cy_i09[ithread] = cy_i08->im;
		_cy_r10[ithread] = cy_r10->re;	_cy_i10[ithread] = cy_i10->re;
		_cy_r11[ithread] = cy_r10->im;	_cy_i11[ithread] = cy_i10->im;
		_cy_r12[ithread] = cy_r12->re;	_cy_i12[ithread] = cy_i12->re;
		_cy_r13[ithread] = cy_r12->im;	_cy_i13[ithread] = cy_i12->im;
		_cy_r14[ithread] = cy_r14->re;	_cy_i14[ithread] = cy_i14->re;
		_cy_r15[ithread] = cy_r14->im;	_cy_i15[ithread] = cy_i14->im;
		_cy_r16[ithread] = cy_r16->re;	_cy_i16[ithread] = cy_i16->re;
		_cy_r17[ithread] = cy_r16->im;	_cy_i17[ithread] = cy_i16->im;
		_cy_r18[ithread] = cy_r18->re;	_cy_i18[ithread] = cy_i18->re;
		_cy_r19[ithread] = cy_r18->im;	_cy_i19[ithread] = cy_i18->im;
		_cy_r20[ithread] = cy_r20->re;	_cy_i20[ithread] = cy_i20->re;
		_cy_r21[ithread] = cy_r20->im;	_cy_i21[ithread] = cy_i20->im;
		_cy_r22[ithread] = cy_r22->re;	_cy_i22[ithread] = cy_i22->re;
		_cy_r23[ithread] = cy_r22->im;	_cy_i23[ithread] = cy_i22->im;
		_cy_r24[ithread] = cy_r24->re;	_cy_i24[ithread] = cy_i24->re;
		_cy_r25[ithread] = cy_r24->im;	_cy_i25[ithread] = cy_i24->im;
		_cy_r26[ithread] = cy_r26->re;	_cy_i26[ithread] = cy_i26->re;
		_cy_r27[ithread] = cy_r26->im;	_cy_i27[ithread] = cy_i26->im;

		if(max_err->re > max_err->im)
			maxerr = max_err->re;
		else
			maxerr = max_err->im;
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
	#endif

		/* Since will lose separate maxerr values when threads are merged, save them after each pass. */
		if(_maxerr[ithread] < maxerr)
		{
			_maxerr[ithread] = maxerr;
		}

	}	/******* END OF PARALLEL FOR-LOOP ********/

	if(!full_pass)break;

	/*   Wraparound carry cleanup loop is here: ***
	!
	!   (1) Invert the radix-28 forward DIF FFT of the first block of 28 complex elements in A and unweight;
	!   (2) Propagate cleanup carries among the real and imaginary parts of the 28 outputs of (1);
	!   (3) Reweight and perform a radix-28 forward DIF FFT on the result of (2);
	!   (4) If any of the exit carries from (2) are nonzero, advance to the next 28 elements and repeat (1-4).
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t00= _cy_r00[CY_THREADS - 1];
		t02= _cy_r01[CY_THREADS - 1];
		t04= _cy_r02[CY_THREADS - 1];
		t06= _cy_r03[CY_THREADS - 1];
		t08= _cy_r04[CY_THREADS - 1];
		t10= _cy_r05[CY_THREADS - 1];
		t12= _cy_r06[CY_THREADS - 1];
		t14= _cy_r07[CY_THREADS - 1];
		t16= _cy_r08[CY_THREADS - 1];
		t18= _cy_r09[CY_THREADS - 1];
		t20= _cy_r10[CY_THREADS - 1];
		t22= _cy_r11[CY_THREADS - 1];
		t24= _cy_r12[CY_THREADS - 1];
		t26= _cy_r13[CY_THREADS - 1];
		t28= _cy_r14[CY_THREADS - 1];
		t30= _cy_r15[CY_THREADS - 1];
		t32= _cy_r16[CY_THREADS - 1];
		t34= _cy_r17[CY_THREADS - 1];
		t36= _cy_r18[CY_THREADS - 1];
		t38= _cy_r19[CY_THREADS - 1];
		t40= _cy_r20[CY_THREADS - 1];
		t42= _cy_r21[CY_THREADS - 1];
		t44= _cy_r22[CY_THREADS - 1];
		t46= _cy_r23[CY_THREADS - 1];
		t48= _cy_r24[CY_THREADS - 1];
		t50= _cy_r25[CY_THREADS - 1];
		t52= _cy_r26[CY_THREADS - 1];
		t54= _cy_r27[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1,"radix28_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
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
		}

		_cy_r00[0] =+t54;	/* ...The wraparound carry is here: */
		_cy_r01[0] = t00;
		_cy_r02[0] = t02;
		_cy_r03[0] = t04;
		_cy_r04[0] = t06;
		_cy_r05[0] = t08;
		_cy_r06[0] = t10;
		_cy_r07[0] = t12;
		_cy_r08[0] = t14;
		_cy_r09[0] = t16;
		_cy_r10[0] = t18;
		_cy_r11[0] = t20;
		_cy_r12[0] = t22;
		_cy_r13[0] = t24;
		_cy_r14[0] = t26;
		_cy_r15[0] = t28;
		_cy_r16[0] = t30;
		_cy_r17[0] = t32;
		_cy_r18[0] = t34;
		_cy_r19[0] = t36;
		_cy_r20[0] = t38;
		_cy_r21[0] = t40;
		_cy_r22[0] = t42;
		_cy_r23[0] = t44;
		_cy_r24[0] = t46;
		_cy_r25[0] = t48;
		_cy_r26[0] = t50;
		_cy_r27[0] = t52;
	}
	else
	{
		t00= _cy_r00[CY_THREADS - 1];	t01= _cy_i00[CY_THREADS - 1];
		t02= _cy_r01[CY_THREADS - 1];	t03= _cy_i01[CY_THREADS - 1];
		t04= _cy_r02[CY_THREADS - 1];	t05= _cy_i02[CY_THREADS - 1];
		t06= _cy_r03[CY_THREADS - 1];	t07= _cy_i03[CY_THREADS - 1];
		t08= _cy_r04[CY_THREADS - 1];	t09= _cy_i04[CY_THREADS - 1];
		t10= _cy_r05[CY_THREADS - 1];	t11= _cy_i05[CY_THREADS - 1];
		t12= _cy_r06[CY_THREADS - 1];	t13= _cy_i06[CY_THREADS - 1];
		t14= _cy_r07[CY_THREADS - 1];	t15= _cy_i07[CY_THREADS - 1];
		t16= _cy_r08[CY_THREADS - 1];	t17= _cy_i08[CY_THREADS - 1];
		t18= _cy_r09[CY_THREADS - 1];	t19= _cy_i09[CY_THREADS - 1];
		t20= _cy_r10[CY_THREADS - 1];	t21= _cy_i10[CY_THREADS - 1];
		t22= _cy_r11[CY_THREADS - 1];	t23= _cy_i11[CY_THREADS - 1];
		t24= _cy_r12[CY_THREADS - 1];	t25= _cy_i12[CY_THREADS - 1];
		t26= _cy_r13[CY_THREADS - 1];	t27= _cy_i13[CY_THREADS - 1];
		t28= _cy_r14[CY_THREADS - 1];	t29= _cy_i14[CY_THREADS - 1];
		t30= _cy_r15[CY_THREADS - 1];	t31= _cy_i15[CY_THREADS - 1];
		t32= _cy_r16[CY_THREADS - 1];	t33= _cy_i16[CY_THREADS - 1];
		t34= _cy_r17[CY_THREADS - 1];	t35= _cy_i17[CY_THREADS - 1];
		t36= _cy_r18[CY_THREADS - 1];	t37= _cy_i18[CY_THREADS - 1];
		t38= _cy_r19[CY_THREADS - 1];	t39= _cy_i19[CY_THREADS - 1];
		t40= _cy_r20[CY_THREADS - 1];	t41= _cy_i20[CY_THREADS - 1];
		t42= _cy_r21[CY_THREADS - 1];	t43= _cy_i21[CY_THREADS - 1];
		t44= _cy_r22[CY_THREADS - 1];	t45= _cy_i22[CY_THREADS - 1];
		t46= _cy_r23[CY_THREADS - 1];	t47= _cy_i23[CY_THREADS - 1];
		t48= _cy_r24[CY_THREADS - 1];	t49= _cy_i24[CY_THREADS - 1];
		t50= _cy_r25[CY_THREADS - 1];	t51= _cy_i25[CY_THREADS - 1];
		t52= _cy_r26[CY_THREADS - 1];	t53= _cy_i26[CY_THREADS - 1];
		t54= _cy_r27[CY_THREADS - 1];	t55= _cy_i27[CY_THREADS - 1];

		for(ithread = CY_THREADS - 1; ithread > 0; ithread--)
		{
			ASSERT(HERE, CY_THREADS > 1,"radix28_ditN_cy_dif1.c: ");	/* Make sure loop only gets executed if multiple threads */
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
		}

		_cy_r00[0] =-t54;	_cy_i00[0] =+t55;	/* ...The 2 Mo"bius carries are here: */
		_cy_r01[0] = t00;	_cy_i01[0] = t01;
		_cy_r02[0] = t02;	_cy_i02[0] = t03;
		_cy_r03[0] = t04;	_cy_i03[0] = t05;
		_cy_r04[0] = t06;	_cy_i04[0] = t07;
		_cy_r05[0] = t08;	_cy_i05[0] = t09;
		_cy_r06[0] = t10;	_cy_i06[0] = t11;
		_cy_r07[0] = t12;	_cy_i07[0] = t13;
		_cy_r08[0] = t14;	_cy_i08[0] = t15;
		_cy_r09[0] = t16;	_cy_i09[0] = t17;
		_cy_r10[0] = t18;	_cy_i10[0] = t19;
		_cy_r11[0] = t20;	_cy_i11[0] = t21;
		_cy_r12[0] = t22;	_cy_i12[0] = t23;
		_cy_r13[0] = t24;	_cy_i13[0] = t25;
		_cy_r14[0] = t26;	_cy_i14[0] = t27;
		_cy_r15[0] = t28;	_cy_i15[0] = t29;
		_cy_r16[0] = t30;	_cy_i16[0] = t31;
		_cy_r17[0] = t32;	_cy_i17[0] = t33;
		_cy_r18[0] = t34;	_cy_i18[0] = t35;
		_cy_r19[0] = t36;	_cy_i19[0] = t37;
		_cy_r20[0] = t38;	_cy_i20[0] = t39;
		_cy_r21[0] = t40;	_cy_i21[0] = t41;
		_cy_r22[0] = t42;	_cy_i22[0] = t43;
		_cy_r23[0] = t44;	_cy_i23[0] = t45;
		_cy_r24[0] = t46;	_cy_i24[0] = t47;
		_cy_r25[0] = t48;	_cy_i25[0] = t49;
		_cy_r26[0] = t50;	_cy_i26[0] = t51;
		_cy_r27[0] = t52;	_cy_i27[0] = t53;
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
			a[j     ] *= radix_inv;
			a[j +p01] *= radix_inv;
			a[j +p02] *= radix_inv;
			a[j +p03] *= radix_inv;	jt = j + p04;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;	jt = j + p08;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;	jt = j + p12;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;	jt = j + p16;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;	jt = j + p20;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;	jt = j + p24;
			a[jt    ] *= radix_inv;
			a[jt+p01] *= radix_inv;
			a[jt+p02] *= radix_inv;
			a[jt+p03] *= radix_inv;
		}
	}
}	/* endfor(outer) */

	t00 = 0;
	for(ithread = 0; ithread < CY_THREADS; ithread++)
	{
		t00 += fabs(_cy_r00[0])+fabs(_cy_r01[0])+fabs(_cy_r02[0])+fabs(_cy_r03[0])+fabs(_cy_r04[0])+fabs(_cy_r05[0])+fabs(_cy_r06[0])+fabs(_cy_r07[0])+fabs(_cy_r08[0])+fabs(_cy_r09[0])+fabs(_cy_r10[0])+fabs(_cy_r11[0])+fabs(_cy_r12[0])+fabs(_cy_r13[0])+fabs(_cy_r14[0])+fabs(_cy_r15[0])+fabs(_cy_r16[0])+fabs(_cy_r17[0])+fabs(_cy_r18[0])+fabs(_cy_r19[0])+fabs(_cy_r20[0])+fabs(_cy_r21[0])+fabs(_cy_r22[0])+fabs(_cy_r23[0])+fabs(_cy_r24[0])+fabs(_cy_r25[0])+fabs(_cy_r26[0])+fabs(_cy_r27[0]);
		t00 += fabs(_cy_i00[0])+fabs(_cy_i01[0])+fabs(_cy_i02[0])+fabs(_cy_i03[0])+fabs(_cy_i04[0])+fabs(_cy_i05[0])+fabs(_cy_i06[0])+fabs(_cy_i07[0])+fabs(_cy_i08[0])+fabs(_cy_i09[0])+fabs(_cy_i10[0])+fabs(_cy_i11[0])+fabs(_cy_i12[0])+fabs(_cy_i13[0])+fabs(_cy_i14[0])+fabs(_cy_i15[0])+fabs(_cy_i16[0])+fabs(_cy_i17[0])+fabs(_cy_i18[0])+fabs(_cy_i19[0])+fabs(_cy_i20[0])+fabs(_cy_i21[0])+fabs(_cy_i22[0])+fabs(_cy_i23[0])+fabs(_cy_i24[0])+fabs(_cy_i25[0])+fabs(_cy_i26[0])+fabs(_cy_i27[0]);

		if(*fracmax < _maxerr[ithread])
			*fracmax = _maxerr[ithread];
	}

	if(t00 != 0.0)
	{
			sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix28_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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

int radix28_ditN_cy_dif1_nochk(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter                 , uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-28 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-28 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix7/8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n28, bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27
		,i,j,j1,j2,jstart,jhi,iroot,root_incr,k1,k2,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	/* The current versions of the macros in dft_macro.h doesn't allow a #define LO_ADD;
	it simply assumes LO_ADD = 1 (and this is explicitly checked at runtime),
	so must use the corresponding versions of the sincos constants :
	*/
	static double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	static double radix_inv, n2inv;
	double re,im,rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i
	,cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27
	,cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27
	,temp,scale;
#if PFETCH
	double *addr, *addp;
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	int ii00,ii01,ii02,ii03,ii04,ii05,ii06,ii07,ii08,ii09,ii10,ii11,ii12,ii13,ii14,ii15,ii16,ii17,ii18,ii19,ii20,ii21,ii22,ii23,ii24,ii25,ii26,ii27;	/* indices into weights arrays (mod NWT) */
	double wt_re,wt_im;									/* Fermat-mod weights stuff */

/*...change n28 and n_div_wt to non-static to work around a gcc compiler bug. */
	n28   = n/28;
	n_div_nwt = n28 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n28)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/28 in radix28_ditN_cy_dif1.\n",iter);
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
		ASSERT(HERE, LO_ADD,"radix28_ditN_cy_dif1.c: LO_ADD");
		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)28));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

		/*   constant index offsets for array load/stores are here.	*/

		p01 = n28;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p10 = p09 + p01;
		p11 = p10 + p01;
		p12 = p11 + p01;
		p13 = p12 + p01;
		p14 = p13 + p01;
		p15 = p14 + p01;
		p16 = p15 + p01;
		p17 = p16 + p01;
		p18 = p17 + p01;
		p19 = p18 + p01;
		p20 = p19 + p01;
		p21 = p20 + p01;
		p22 = p21 + p01;
		p23 = p22 + p01;
		p24 = p23 + p01;
		p25 = p24 + p01;
		p26 = p25 + p01;
		p27 = p26 + p01;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p17 = p17 + ( (p17 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p19 = p19 + ( (p19 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p21 = p21 + ( (p21 >> DAT_BITS) << PAD_BITS );
		p22 = p22 + ( (p22 >> DAT_BITS) << PAD_BITS );
		p23 = p23 + ( (p23 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p25 = p25 + ( (p25 >> DAT_BITS) << PAD_BITS );
		p26 = p26 + ( (p26 >> DAT_BITS) << PAD_BITS );
		p27 = p27 + ( (p27 >> DAT_BITS) << PAD_BITS );

		/* For Fermat-mod, since 'adjacent' words are actually stride-2 separated
		in terms of the floating residue array, block boundaries have half the i-index
		(e.g. as in sw*i and bw*i) value they do in the Mersenne-mod case:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			bjmodnini=0;
			for(j=0; j < n28; j++)
			{
				bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
			}
		}
		else
		{
			bjmodnini=0;
			for(j=0; j < n28/2; j++)
			{
				bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
			}
		}
	}

/*...The radix-28 final DIT pass is here.	*/

	/* init carries	*/
	cy_r00= 0;	cy_i00= 0;
	cy_r01= 0;	cy_i01= 0;
	cy_r02= 0;	cy_i02= 0;
	cy_r03= 0;	cy_i03= 0;
	cy_r04= 0;	cy_i04= 0;
	cy_r05= 0;	cy_i05= 0;
	cy_r06= 0;	cy_i06= 0;
	cy_r07= 0;	cy_i07= 0;
	cy_r08= 0;	cy_i08= 0;
	cy_r09= 0;	cy_i09= 0;
	cy_r10= 0;	cy_i10= 0;
	cy_r11= 0;	cy_i11= 0;
	cy_r12= 0;	cy_i12= 0;
	cy_r13= 0;	cy_i13= 0;
	cy_r14= 0;	cy_i14= 0;
	cy_r15= 0;	cy_i15= 0;
	cy_r16= 0;	cy_i16= 0;
	cy_r17= 0;	cy_i17= 0;
	cy_r18= 0;	cy_i18= 0;
	cy_r19= 0;	cy_i19= 0;
	cy_r20= 0;	cy_i20= 0;
	cy_r21= 0;	cy_i21= 0;
	cy_r22= 0;	cy_i22= 0;
	cy_r23= 0;	cy_i23= 0;
	cy_r24= 0;	cy_i24= 0;
	cy_r25= 0;	cy_i25= 0;
	cy_r26= 0;	cy_i26= 0;
	cy_r27= 0;	cy_i27= 0;

	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		cy_r00 = -2;
	}

	iroot = 0;	/* init sincos array index	*/
	root_incr = 1;	/* init sincos array index increment (set = 1 for normal carry pass, = 0 for wrapper pass)	*/

	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart = 0;
		jhi = jstart+nwt-1;
		khi = n_div_nwt;
	}
	else
	{
		jstart = 0;
		jhi = n_div_nwt;
		khi = 1;
	}

for(outer=0; outer <= 1; outer++)
{
	i = 0;		/* Index into the BASE and BASEINV arrays. */
	/* If bw > 0 (i.e. n does not divide p), lowest-order digit always a bigword:	*/
	if(bw > 0)
		i = 1;

	bjmodn00= 0;
	bjmodn01= bjmodnini;
	bjmodn02= bjmodn01+bjmodnini-n; bjmodn02= bjmodn02+ ( (-(int)((uint32)bjmodn02>> 31)) & n);
	bjmodn03= bjmodn02+bjmodnini-n; bjmodn03= bjmodn03+ ( (-(int)((uint32)bjmodn03>> 31)) & n);
	bjmodn04= bjmodn03+bjmodnini-n; bjmodn04= bjmodn04+ ( (-(int)((uint32)bjmodn04>> 31)) & n);
	bjmodn05= bjmodn04+bjmodnini-n; bjmodn05= bjmodn05+ ( (-(int)((uint32)bjmodn05>> 31)) & n);
	bjmodn06= bjmodn05+bjmodnini-n; bjmodn06= bjmodn06+ ( (-(int)((uint32)bjmodn06>> 31)) & n);
	bjmodn07= bjmodn06+bjmodnini-n; bjmodn07= bjmodn07+ ( (-(int)((uint32)bjmodn07>> 31)) & n);
	bjmodn08= bjmodn07+bjmodnini-n; bjmodn08= bjmodn08+ ( (-(int)((uint32)bjmodn08>> 31)) & n);
	bjmodn09= bjmodn08+bjmodnini-n; bjmodn09= bjmodn09+ ( (-(int)((uint32)bjmodn09>> 31)) & n);
	bjmodn10= bjmodn09+bjmodnini-n; bjmodn10= bjmodn10+ ( (-(int)((uint32)bjmodn10>> 31)) & n);
	bjmodn11= bjmodn10+bjmodnini-n; bjmodn11= bjmodn11+ ( (-(int)((uint32)bjmodn11>> 31)) & n);
	bjmodn12= bjmodn11+bjmodnini-n; bjmodn12= bjmodn12+ ( (-(int)((uint32)bjmodn12>> 31)) & n);
	bjmodn13= bjmodn12+bjmodnini-n; bjmodn13= bjmodn13+ ( (-(int)((uint32)bjmodn13>> 31)) & n);
	bjmodn14= bjmodn13+bjmodnini-n; bjmodn14= bjmodn14+ ( (-(int)((uint32)bjmodn14>> 31)) & n);
	bjmodn15= bjmodn14+bjmodnini-n; bjmodn15= bjmodn15+ ( (-(int)((uint32)bjmodn15>> 31)) & n);
	bjmodn16= bjmodn15+bjmodnini-n; bjmodn16= bjmodn16+ ( (-(int)((uint32)bjmodn16>> 31)) & n);
	bjmodn17= bjmodn16+bjmodnini-n; bjmodn17= bjmodn17+ ( (-(int)((uint32)bjmodn17>> 31)) & n);
	bjmodn18= bjmodn17+bjmodnini-n; bjmodn18= bjmodn18+ ( (-(int)((uint32)bjmodn18>> 31)) & n);
	bjmodn19= bjmodn18+bjmodnini-n; bjmodn19= bjmodn19+ ( (-(int)((uint32)bjmodn19>> 31)) & n);
	bjmodn20= bjmodn19+bjmodnini-n; bjmodn20= bjmodn20+ ( (-(int)((uint32)bjmodn20>> 31)) & n);
	bjmodn21= bjmodn20+bjmodnini-n; bjmodn21= bjmodn21+ ( (-(int)((uint32)bjmodn21>> 31)) & n);
	bjmodn22= bjmodn21+bjmodnini-n; bjmodn22= bjmodn22+ ( (-(int)((uint32)bjmodn22>> 31)) & n);
	bjmodn23= bjmodn22+bjmodnini-n; bjmodn23= bjmodn23+ ( (-(int)((uint32)bjmodn23>> 31)) & n);
	bjmodn24= bjmodn23+bjmodnini-n; bjmodn24= bjmodn24+ ( (-(int)((uint32)bjmodn24>> 31)) & n);
	bjmodn25= bjmodn24+bjmodnini-n; bjmodn25= bjmodn25+ ( (-(int)((uint32)bjmodn25>> 31)) & n);
	bjmodn26= bjmodn25+bjmodnini-n; bjmodn26= bjmodn26+ ( (-(int)((uint32)bjmodn26>> 31)) & n);
	bjmodn27= bjmodn26+bjmodnini-n; bjmodn27= bjmodn27+ ( (-(int)((uint32)bjmodn27>> 31)) & n);

	/* For Fermat-mod, IBDWT access patterns repeat with period NWT = {odd part of radix0},
	so for even radix0's only really need that many bjmodn and ii's, but that would require
	specialized carry macros that don't update ii and bjmodn - not worth the trouble.
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		col=0;
		co2=(n >> nwt_bits)-1+28;
		co3=co2-28;		/* At the start of each new j-loop, co3=co2-radix(1)	*/
	}
	else
	{
		/* indices into IBDWT weights arrays (mod NWT) is here: */
		ii00= 0;
		ii01= (SW_DIV_N*n28/2) % nwt;
		ii02= (ii01+ ii01) % nwt;
		ii03= (ii02+ ii01) % nwt;
		ii04= (ii03+ ii01) % nwt;
		ii05= (ii04+ ii01) % nwt;
		ii06= (ii05+ ii01) % nwt;
		ii07= (ii06+ ii01) % nwt;
		ii08= (ii07+ ii01) % nwt;
		ii09= (ii08+ ii01) % nwt;
		ii10= (ii09+ ii01) % nwt;
		ii11= (ii10+ ii01) % nwt;
		ii12= (ii11+ ii01) % nwt;
		ii13= (ii12+ ii01) % nwt;
		ii14= (ii13+ ii01) % nwt;
		ii15= (ii14+ ii01) % nwt;
		ii16= (ii15+ ii01) % nwt;
		ii17= (ii16+ ii01) % nwt;
		ii18= (ii17+ ii01) % nwt;
		ii19= (ii18+ ii01) % nwt;
		ii20= (ii19+ ii01) % nwt;
		ii21= (ii20+ ii01) % nwt;
		ii22= (ii21+ ii01) % nwt;
		ii23= (ii22+ ii01) % nwt;
		ii24= (ii23+ ii01) % nwt;
		ii25= (ii24+ ii01) % nwt;
		ii26= (ii25+ ii01) % nwt;
		ii27= (ii26+ ii01) % nwt;

		/* Start this value off at N in Fermat-mod case, so (bjmodn >= sw) check in
		fermat_carry_norm_errcheck (cf. carry.h) yields a bigword (i == 1) for j= 0:
		*/
		bjmodn00= n;
		bjmodn07= n;
		bjmodn14= n;
		bjmodn21= n;
	}

	for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
	{
		for(j=jstart; j<jhi; j += 2)	/* Each inner loop execution processes (radix(1)*nwt) array data.	*/
		{
		#ifdef USE_SSE2
			j1 = (j & mask01) + br4[j&3];
			j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
		#else
			j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		#endif
			j2 = j1+RE_IM_STRIDE;

/*...The radix-28 DIT pass is here:	*/

/* EWM: 10/18/04: We swap the odd-index outputs of each of the radix-4 DIT transforms (1<=>3, 5<=>7, etc.) so that the indexing
                  of all the radix-7 transforms (really just the 2nd and 4th of these) winds up being in-place. This allows us
                  to properly re-use the ajp1 variables in the carry-pass version of this routine.
*/

/*...gather the needed data (28 64-bit complex, i.e. 56 64-bit reals) and do 7 radix-4 transforms...*/
                     /*                                      outputs                                      */ /*                          inputs                           */
	RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a1p00r,a1p00i,a1p03r,a1p03i,a1p02r,a1p02i,a1p01r,a1p01i,rt,it);
	RADIX_04_DIT(a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p13],a[j2+p13],a[j1+p12],a[j2+p12],a1p04r,a1p04i,a1p07r,a1p07i,a1p06r,a1p06i,a1p05r,a1p05i,rt,it);
	RADIX_04_DIT(a[j1+p25],a[j2+p25],a[j1+p24],a[j2+p24],a[j1+p26],a[j2+p26],a[j1+p27],a[j2+p27],a1p08r,a1p08i,a1p11r,a1p11i,a1p10r,a1p10i,a1p09r,a1p09i,rt,it);
	RADIX_04_DIT(a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a1p12r,a1p12i,a1p15r,a1p15i,a1p14r,a1p14i,a1p13r,a1p13i,rt,it);
	RADIX_04_DIT(a[j1+p22],a[j2+p22],a[j1+p23],a[j2+p23],a[j1+p20],a[j2+p20],a[j1+p21],a[j2+p21],a1p16r,a1p16i,a1p19r,a1p19i,a1p18r,a1p18i,a1p17r,a1p17i,rt,it);
	RADIX_04_DIT(a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05],a1p20r,a1p20i,a1p23r,a1p23i,a1p22r,a1p22i,a1p21r,a1p21i,rt,it);
	RADIX_04_DIT(a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],a1p24r,a1p24i,a1p27r,a1p27i,a1p26r,a1p26i,a1p25r,a1p25i,rt,it);

/*...and now do 4 radix-4 transforms...*/
                     /*                                                   inputs                                                  */ /*               intermediates              */ /*                                                                     outputs                                                                         */
	RADIX_07_DFT(a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p00r,a1p00i,a1p08r,a1p08i,a1p16r,a1p16i,a1p24r,a1p24i,a1p04r,a1p04i,a1p12r,a1p12i,a1p20r,a1p20i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	RADIX_07_DFT(a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p07r,a1p07i,a1p15r,a1p15i,a1p23r,a1p23i,a1p03r,a1p03i,a1p11r,a1p11i,a1p19r,a1p19i,a1p27r,a1p27i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	RADIX_07_DFT(a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p14r,a1p14i,a1p22r,a1p22i,a1p02r,a1p02i,a1p10r,a1p10i,a1p18r,a1p18i,a1p26r,a1p26i,a1p06r,a1p06i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	RADIX_07_DFT(a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p21r,a1p21i,a1p01r,a1p01i,a1p09r,a1p09i,a1p17r,a1p17i,a1p25r,a1p25i,a1p05r,a1p05i,a1p13r,a1p13i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);

/*...Now do the carries. Since the outputs would
	normally be getting dispatched to 28 separate blocks of the A-array, we need 28 separate carries.	*/

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
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
			 cmplx_carry_norm_nocheck0(a1p00r,a1p00i,cy_r00,bjmodn00   );
			cmplx_carry_norm_nocheck(a1p01r,a1p01i,cy_r01,bjmodn01,1 );
			cmplx_carry_norm_nocheck(a1p02r,a1p02i,cy_r02,bjmodn02,2 );
			cmplx_carry_norm_nocheck(a1p03r,a1p03i,cy_r03,bjmodn03,3 );
			cmplx_carry_norm_nocheck(a1p04r,a1p04i,cy_r04,bjmodn04,4 );
			cmplx_carry_norm_nocheck(a1p05r,a1p05i,cy_r05,bjmodn05,5 );
			cmplx_carry_norm_nocheck(a1p06r,a1p06i,cy_r06,bjmodn06,6 );
			cmplx_carry_norm_nocheck(a1p07r,a1p07i,cy_r07,bjmodn07,7 );
			cmplx_carry_norm_nocheck(a1p08r,a1p08i,cy_r08,bjmodn08,8 );
			cmplx_carry_norm_nocheck(a1p09r,a1p09i,cy_r09,bjmodn09,9 );
			cmplx_carry_norm_nocheck(a1p10r,a1p10i,cy_r10,bjmodn10,10);
			cmplx_carry_norm_nocheck(a1p11r,a1p11i,cy_r11,bjmodn11,11);
			cmplx_carry_norm_nocheck(a1p12r,a1p12i,cy_r12,bjmodn12,12);
			cmplx_carry_norm_nocheck(a1p13r,a1p13i,cy_r13,bjmodn13,13);
			cmplx_carry_norm_nocheck(a1p14r,a1p14i,cy_r14,bjmodn14,14);
			cmplx_carry_norm_nocheck(a1p15r,a1p15i,cy_r15,bjmodn15,15);
			cmplx_carry_norm_nocheck(a1p16r,a1p16i,cy_r16,bjmodn16,16);
			cmplx_carry_norm_nocheck(a1p17r,a1p17i,cy_r17,bjmodn17,17);
			cmplx_carry_norm_nocheck(a1p18r,a1p18i,cy_r18,bjmodn18,18);
			cmplx_carry_norm_nocheck(a1p19r,a1p19i,cy_r19,bjmodn19,19);
			cmplx_carry_norm_nocheck(a1p20r,a1p20i,cy_r20,bjmodn20,20);
			cmplx_carry_norm_nocheck(a1p21r,a1p21i,cy_r21,bjmodn21,21);
			cmplx_carry_norm_nocheck(a1p22r,a1p22i,cy_r22,bjmodn22,22);
			cmplx_carry_norm_nocheck(a1p23r,a1p23i,cy_r23,bjmodn23,23);
			cmplx_carry_norm_nocheck(a1p24r,a1p24i,cy_r24,bjmodn24,24);
			cmplx_carry_norm_nocheck(a1p25r,a1p25i,cy_r25,bjmodn25,25);
			cmplx_carry_norm_nocheck(a1p26r,a1p26i,cy_r26,bjmodn26,26);
			cmplx_carry_norm_nocheck(a1p27r,a1p27i,cy_r27,bjmodn27,27);

			i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
						and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/
		}
		else
		{
			fermat_carry_norm_nocheck(a1p00r,a1p00i,cy_r00,cy_i00,ii00,bjmodn00,0 *n28);
			fermat_carry_norm_nocheck(a1p01r,a1p01i,cy_r01,cy_i01,ii01,bjmodn01,1 *n28);
			fermat_carry_norm_nocheck(a1p02r,a1p02i,cy_r02,cy_i02,ii02,bjmodn02,2 *n28);
			fermat_carry_norm_nocheck(a1p03r,a1p03i,cy_r03,cy_i03,ii03,bjmodn03,3 *n28);
			fermat_carry_norm_nocheck(a1p04r,a1p04i,cy_r04,cy_i04,ii04,bjmodn04,4 *n28);
			fermat_carry_norm_nocheck(a1p05r,a1p05i,cy_r05,cy_i05,ii05,bjmodn05,5 *n28);
			fermat_carry_norm_nocheck(a1p06r,a1p06i,cy_r06,cy_i06,ii06,bjmodn06,6 *n28);
			fermat_carry_norm_nocheck(a1p07r,a1p07i,cy_r07,cy_i07,ii07,bjmodn07,7 *n28);
			fermat_carry_norm_nocheck(a1p08r,a1p08i,cy_r08,cy_i08,ii08,bjmodn08,8 *n28);
			fermat_carry_norm_nocheck(a1p09r,a1p09i,cy_r09,cy_i09,ii09,bjmodn09,9 *n28);
			fermat_carry_norm_nocheck(a1p10r,a1p10i,cy_r10,cy_i10,ii10,bjmodn10,10*n28);
			fermat_carry_norm_nocheck(a1p11r,a1p11i,cy_r11,cy_i11,ii11,bjmodn11,11*n28);
			fermat_carry_norm_nocheck(a1p12r,a1p12i,cy_r12,cy_i12,ii12,bjmodn12,12*n28);
			fermat_carry_norm_nocheck(a1p13r,a1p13i,cy_r13,cy_i13,ii13,bjmodn13,13*n28);
			fermat_carry_norm_nocheck(a1p14r,a1p14i,cy_r14,cy_i14,ii14,bjmodn14,14*n28);
			fermat_carry_norm_nocheck(a1p15r,a1p15i,cy_r15,cy_i15,ii15,bjmodn15,15*n28);
			fermat_carry_norm_nocheck(a1p16r,a1p16i,cy_r16,cy_i16,ii16,bjmodn16,16*n28);
			fermat_carry_norm_nocheck(a1p17r,a1p17i,cy_r17,cy_i17,ii17,bjmodn17,17*n28);
			fermat_carry_norm_nocheck(a1p18r,a1p18i,cy_r18,cy_i18,ii18,bjmodn18,18*n28);
			fermat_carry_norm_nocheck(a1p19r,a1p19i,cy_r19,cy_i19,ii19,bjmodn19,19*n28);
			fermat_carry_norm_nocheck(a1p20r,a1p20i,cy_r20,cy_i20,ii20,bjmodn20,20*n28);
			fermat_carry_norm_nocheck(a1p21r,a1p21i,cy_r21,cy_i21,ii21,bjmodn21,21*n28);
			fermat_carry_norm_nocheck(a1p22r,a1p22i,cy_r22,cy_i22,ii22,bjmodn22,22*n28);
			fermat_carry_norm_nocheck(a1p23r,a1p23i,cy_r23,cy_i23,ii23,bjmodn23,23*n28);
			fermat_carry_norm_nocheck(a1p24r,a1p24i,cy_r24,cy_i24,ii24,bjmodn24,24*n28);
			fermat_carry_norm_nocheck(a1p25r,a1p25i,cy_r25,cy_i25,ii25,bjmodn25,25*n28);
			fermat_carry_norm_nocheck(a1p26r,a1p26i,cy_r26,cy_i26,ii26,bjmodn26,26*n28);
			fermat_carry_norm_nocheck(a1p27r,a1p27i,cy_r27,cy_i27,ii27,bjmodn27,27*n28);
		}

/*...The radix-28 DIF pass is here:	*/
#if PFETCH
addr = &a[j1];
prefetch_p_doubles(addr);
#endif

/*...gather the needed data (28 64-bit complex, i.e. 56 64-bit reals) and do 4 radix-7 transforms...*/
                     /*                                                                      inputs                                                                         */ /*               intermediates              */ /*                                                  outputs                                                  */
#if PFETCH
	RADIX_07_DFT_PFETCH(a1p00r,a1p00i,a1p24r,a1p24i,a1p20r,a1p20i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im,addr,addp,p01,p02,p03);
	RADIX_07_DFT_PFETCH(a1p21r,a1p21i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,a1p05r,a1p05i,a1p01r,a1p01i,a1p25r,a1p25i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im,addr,addp,p04,p05,p06);
	RADIX_07_DFT_PFETCH(a1p14r,a1p14i,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p26r,a1p26i,a1p22r,a1p22i,a1p18r,a1p18i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im,addr,addp,p07,p08,p09);
	RADIX_07_DFT_PFETCH(a1p07r,a1p07i,a1p03r,a1p03i,a1p27r,a1p27i,a1p23r,a1p23i,a1p19r,a1p19i,a1p15r,a1p15i,a1p11r,a1p11i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im,addr,addp,p10,p11,p12);
#else
	RADIX_07_DFT       (a1p00r,a1p00i,a1p24r,a1p24i,a1p20r,a1p20i,a1p16r,a1p16i,a1p12r,a1p12i,a1p08r,a1p08i,a1p04r,a1p04i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	RADIX_07_DFT       (a1p21r,a1p21i,a1p17r,a1p17i,a1p13r,a1p13i,a1p09r,a1p09i,a1p05r,a1p05i,a1p01r,a1p01i,a1p25r,a1p25i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	RADIX_07_DFT       (a1p14r,a1p14i,a1p10r,a1p10i,a1p06r,a1p06i,a1p02r,a1p02i,a1p26r,a1p26i,a1p22r,a1p22i,a1p18r,a1p18i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	RADIX_07_DFT       (a1p07r,a1p07i,a1p03r,a1p03i,a1p27r,a1p27i,a1p23r,a1p23i,a1p19r,a1p19i,a1p15r,a1p15i,a1p11r,a1p11i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
#endif

/*...and now do 7 radix-4 transforms...*/
                     /*                          inputs                           */ /*                                      outputs                                      */
#if PFETCH
	addp = addr+p13;
	prefetch_p_doubles(addp);

	RADIX_04_DIF_PFETCH(a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],rt,it,addr,addp,p14,p15);
	RADIX_04_DIF_PFETCH(a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,a[j1+p25],a[j2+p25],a[j1+p24],a[j2+p24],a[j1+p27],a[j2+p27],a[j1+p26],a[j2+p26],rt,it,addr,addp,p16,p17);
	RADIX_04_DIF_PFETCH(a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a[j1+p22],a[j2+p22],a[j1+p23],a[j2+p23],a[j1+p21],a[j2+p21],a[j1+p20],a[j2+p20],rt,it,addr,addp,p18,p19);
	RADIX_04_DIF_PFETCH(a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p18],a[j2+p18],a[j1+p19],a[j2+p19],rt,it,addr,addp,p20,p21);
	RADIX_04_DIF_PFETCH(a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],rt,it,addr,addp,p22,p23);
	RADIX_04_DIF_PFETCH(a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p11],a[j2+p11],a[j1+p10],a[j2+p10],rt,it,addr,addp,p24,p25);
	RADIX_04_DIF_PFETCH(a1p24r,a1p24i,a1p25r,a1p25i,a1p26r,a1p26i,a1p27r,a1p27i,a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],rt,it,addr,addp,p26,p27);
#else
	RADIX_04_DIF       (a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],rt,it);
	RADIX_04_DIF       (a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,a[j1+p25],a[j2+p25],a[j1+p24],a[j2+p24],a[j1+p27],a[j2+p27],a[j1+p26],a[j2+p26],rt,it);
	RADIX_04_DIF       (a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a[j1+p22],a[j2+p22],a[j1+p23],a[j2+p23],a[j1+p21],a[j2+p21],a[j1+p20],a[j2+p20],rt,it);
	RADIX_04_DIF       (a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p18],a[j2+p18],a[j1+p19],a[j2+p19],rt,it);
	RADIX_04_DIF       (a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],rt,it);
	RADIX_04_DIF       (a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p11],a[j2+p11],a[j1+p10],a[j2+p10],rt,it);
	RADIX_04_DIF       (a1p24r,a1p24i,a1p25r,a1p25i,a1p26r,a1p26i,a1p27r,a1p27i,a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],rt,it);
#endif

		iroot += root_incr;		/* increment sincos index.	*/

		}

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			jstart += nwt;
			jhi    += nwt;

			col += 28;
			co3 -= 28;
		}
	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-28 forward DIF FFT of the first block of 28 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 28 outputs of (1);
!   (3) Reweight and perform a radix-28 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 28 elements and repeat (1-4).
*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t1    = cy_r27;
		cy_r27= cy_r26;
		cy_r26= cy_r25;
		cy_r25= cy_r24;
		cy_r24= cy_r23;
		cy_r23= cy_r22;
		cy_r22= cy_r21;
		cy_r21= cy_r20;
		cy_r20= cy_r19;
		cy_r19= cy_r18;
		cy_r18= cy_r17;
		cy_r17= cy_r16;
		cy_r16= cy_r15;
		cy_r15= cy_r14;
		cy_r14= cy_r13;
		cy_r13= cy_r12;
		cy_r12= cy_r11;
		cy_r11= cy_r10;
		cy_r10= cy_r09;
		cy_r09= cy_r08;
		cy_r08= cy_r07;
		cy_r07= cy_r06;
		cy_r06= cy_r05;
		cy_r05= cy_r04;
		cy_r04= cy_r03;
		cy_r03= cy_r02;
		cy_r02= cy_r01;
		cy_r01= cy_r00;
		cy_r00=    t1 ;
	}
	else
	{
		/* ...The 2 Mo"bius carries are here: */
		t1    = cy_r27;	t2    = cy_i27;
		cy_r27= cy_r26;	cy_i27= cy_i26;
		cy_r26= cy_r25;	cy_i26= cy_i25;
		cy_r25= cy_r24;	cy_i25= cy_i24;
		cy_r24= cy_r23;	cy_i24= cy_i23;
		cy_r23= cy_r22;	cy_i23= cy_i22;
		cy_r22= cy_r21;	cy_i22= cy_i21;
		cy_r21= cy_r20;	cy_i21= cy_i20;
		cy_r20= cy_r19;	cy_i20= cy_i19;
		cy_r19= cy_r18;	cy_i19= cy_i18;
		cy_r18= cy_r17;	cy_i18= cy_i17;
		cy_r17= cy_r16;	cy_i17= cy_i16;
		cy_r16= cy_r15;	cy_i16= cy_i15;
		cy_r15= cy_r14;	cy_i15= cy_i14;
		cy_r14= cy_r13;	cy_i14= cy_i13;
		cy_r13= cy_r12;	cy_i13= cy_i12;
		cy_r12= cy_r11;	cy_i12= cy_i11;
		cy_r11= cy_r10;	cy_i11= cy_i10;
		cy_r10= cy_r09;	cy_i10= cy_i09;
		cy_r09= cy_r08;	cy_i09= cy_i08;
		cy_r08= cy_r07;	cy_i08= cy_i07;
		cy_r07= cy_r06;	cy_i07= cy_i06;
		cy_r06= cy_r05;	cy_i06= cy_i05;
		cy_r05= cy_r04;	cy_i05= cy_i04;
		cy_r04= cy_r03;	cy_i04= cy_i03;
		cy_r03= cy_r02;	cy_i03= cy_i02;
		cy_r02= cy_r01;	cy_i02= cy_i01;
		cy_r01= cy_r00;	cy_i01= cy_i00;
		cy_r00=   -t2 ;	cy_i00=   +t1 ;
	}

	iroot = 0;
	root_incr = 0;
	scale = 1;

	jstart = 0;
	/*
	For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
	*/
	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		jhi =15;
	}
	else
	{
		jhi = 7;
	}
	khi = 1;

	for(j=0; j<=jhi; j++)
	{
		a[j    ] *= radix_inv;
		a[j+p01] *= radix_inv;
		a[j+p02] *= radix_inv;
		a[j+p03] *= radix_inv;
		a[j+p04] *= radix_inv;
		a[j+p05] *= radix_inv;
		a[j+p06] *= radix_inv;
		a[j+p07] *= radix_inv;
		a[j+p08] *= radix_inv;
		a[j+p09] *= radix_inv;
		a[j+p10] *= radix_inv;
		a[j+p11] *= radix_inv;
		a[j+p12] *= radix_inv;
		a[j+p13] *= radix_inv;
		a[j+p14] *= radix_inv;
		a[j+p15] *= radix_inv;
		a[j+p16] *= radix_inv;
		a[j+p17] *= radix_inv;
		a[j+p18] *= radix_inv;
		a[j+p19] *= radix_inv;
		a[j+p20] *= radix_inv;
		a[j+p21] *= radix_inv;
		a[j+p22] *= radix_inv;
		a[j+p23] *= radix_inv;
		a[j+p24] *= radix_inv;
		a[j+p25] *= radix_inv;
		a[j+p26] *= radix_inv;
		a[j+p27] *= radix_inv;
	}
}

	if(fabs(cy_r00)+fabs(cy_r01)+fabs(cy_r02)+fabs(cy_r03)+fabs(cy_r04)+fabs(cy_r05)+fabs(cy_r06)+fabs(cy_r07)+fabs(cy_r08)+fabs(cy_r09)+fabs(cy_r10)+fabs(cy_r11)+fabs(cy_r12)+fabs(cy_r13)+fabs(cy_r14)+fabs(cy_r15)+fabs(cy_r16)+fabs(cy_r17)+fabs(cy_r18)+fabs(cy_r19)+fabs(cy_r20)+fabs(cy_r21)+fabs(cy_r22)+fabs(cy_r23)+fabs(cy_r24)+fabs(cy_r25)+fabs(cy_r26)+fabs(cy_r27)
		+fabs(cy_i00)+fabs(cy_i01)+fabs(cy_i02)+fabs(cy_i03)+fabs(cy_i04)+fabs(cy_i05)+fabs(cy_i06)+fabs(cy_i07)+fabs(cy_i08)+fabs(cy_i09)+fabs(cy_i10)+fabs(cy_i11)+fabs(cy_i12)+fabs(cy_i13)+fabs(cy_i14)+fabs(cy_i15)+fabs(cy_i16)+fabs(cy_i17)+fabs(cy_i18)+fabs(cy_i19)+fabs(cy_i20)+fabs(cy_i21)+fabs(cy_i22)+fabs(cy_i23)+fabs(cy_i24)+fabs(cy_i25)+fabs(cy_i26)+fabs(cy_i27) != 0.0)
	{
			sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix28_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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

void radix28_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-28 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix7_dif_pass for details on the radix-15 subtransforms.
*/
	int j,j1,j2;
	static int n28,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27, first_entry=TRUE;
	static double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	double re,im,rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i;

	if(!first_entry && (n/28) != n28)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n28=n/28;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n28;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p10 = p09 + p01;
		p11 = p10 + p01;
		p12 = p11 + p01;
		p13 = p12 + p01;
		p14 = p13 + p01;
		p15 = p14 + p01;
		p16 = p15 + p01;
		p17 = p16 + p01;
		p18 = p17 + p01;
		p19 = p18 + p01;
		p20 = p19 + p01;
		p21 = p20 + p01;
		p22 = p21 + p01;
		p23 = p22 + p01;
		p24 = p23 + p01;
		p25 = p24 + p01;
		p26 = p25 + p01;
		p27 = p26 + p01;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p17 = p17 + ( (p17 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p19 = p19 + ( (p19 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p21 = p21 + ( (p21 >> DAT_BITS) << PAD_BITS );
		p22 = p22 + ( (p22 >> DAT_BITS) << PAD_BITS );
		p23 = p23 + ( (p23 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p25 = p25 + ( (p25 >> DAT_BITS) << PAD_BITS );
		p26 = p26 + ( (p26 >> DAT_BITS) << PAD_BITS );
		p27 = p27 + ( (p27 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-28 pass is here.	*/

	for(j=0; j < n28; j += 2)
	{
	#if defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*...gather the needed data (28 64-bit complex, i.e. 56 64-bit reals) and do 4 radix-7 transforms...*/

	/*
		Twiddleless version requires us to swap inputs as follows:
		indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27
			  -> 0,21,14, 7,24,17,10, 3,20,13, 6,27,16, 9, 2,23,12, 5,26,19, 8, 1,22,15, 4,25,18,11
		I.e. start out with first septet of indices {0,4,8,12,16,20,24}, permute those according to
		{0,4,8,12,16,20,24}*27%28 = {0,24,20,16,12,8,4}, then each is head of a length-4 list of indices with decrement 7.
	*/
						 /*                                                                      inputs                                                                         */ /*               intermediates              */ /*                                                  outputs                                                  */
		RADIX_07_DFT(a[j1    ],a[j2    ],a[j1+p24],a[j2+p24],a[j1+p20],a[j2+p20],a[j1+p16],a[j2+p16],a[j1+p12],a[j2+p12],a[j1+p08],a[j2+p08],a[j1+p04],a[j2+p04],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a[j1+p21],a[j2+p21],a[j1+p17],a[j2+p17],a[j1+p13],a[j2+p13],a[j1+p09],a[j2+p09],a[j1+p05],a[j2+p05],a[j1+p01],a[j2+p01],a[j1+p25],a[j2+p25],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a[j1+p14],a[j2+p14],a[j1+p10],a[j2+p10],a[j1+p06],a[j2+p06],a[j1+p02],a[j2+p02],a[j1+p26],a[j2+p26],a[j1+p22],a[j2+p22],a[j1+p18],a[j2+p18],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a[j1+p07],a[j2+p07],a[j1+p03],a[j2+p03],a[j1+p27],a[j2+p27],a[j1+p23],a[j2+p23],a[j1+p19],a[j2+p19],a[j1+p15],a[j2+p15],a[j1+p11],a[j2+p11],t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);

	/*...and now do 7 radix-4 transforms...*/
						 /*                          inputs                           */ /*                                      outputs                                      */
		RADIX_04_DIF(a1p00r,a1p00i,a1p01r,a1p01i,a1p02r,a1p02i,a1p03r,a1p03i,a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p02],a[j2+p02],a[j1+p03],a[j2+p03],rt,it);
		RADIX_04_DIF(a1p04r,a1p04i,a1p05r,a1p05i,a1p06r,a1p06i,a1p07r,a1p07i,a[j1+p25],a[j2+p25],a[j1+p24],a[j2+p24],a[j1+p27],a[j2+p27],a[j1+p26],a[j2+p26],rt,it);
		RADIX_04_DIF(a1p08r,a1p08i,a1p09r,a1p09i,a1p10r,a1p10i,a1p11r,a1p11i,a[j1+p22],a[j2+p22],a[j1+p23],a[j2+p23],a[j1+p21],a[j2+p21],a[j1+p20],a[j2+p20],rt,it);
		RADIX_04_DIF(a1p12r,a1p12i,a1p13r,a1p13i,a1p14r,a1p14i,a1p15r,a1p15i,a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p18],a[j2+p18],a[j1+p19],a[j2+p19],rt,it);
		RADIX_04_DIF(a1p16r,a1p16i,a1p17r,a1p17i,a1p18r,a1p18i,a1p19r,a1p19i,a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],rt,it);
		RADIX_04_DIF(a1p20r,a1p20i,a1p21r,a1p21i,a1p22r,a1p22i,a1p23r,a1p23i,a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p11],a[j2+p11],a[j1+p10],a[j2+p10],rt,it);
		RADIX_04_DIF(a1p24r,a1p24i,a1p25r,a1p25i,a1p26r,a1p26i,a1p27r,a1p27i,a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p05],a[j2+p05],a[j1+p04],a[j2+p04],rt,it);
	}
}

/*
   0 (   0)  142.0000000000   129.0000000000   142.0000000000   129.0000000000
   1 (  14)  -18.0000000000   -17.0000000000   -18.0000000000   -17.0000000000
   2 (   7)  -33.0000000000   -22.0000000000   -15.0000000000   -14.0000000000
   3 (  21)  -15.0000000000   -14.0000000000   -33.0000000000   -22.0000000000
   4 (   1)  -11.6595587366   -19.7679307866    -7.6577999302     8.5976313922
   5 (  15)    2.8294700831     0.5106027780    32.4689557080    -0.3465505441
   6 (   8)   46.1480382688    -2.2866689587     9.9796188501    -9.6408663259
   7 (  22)  -17.5691832932     0.3833862984    -7.6241120460    -1.5581921548
   8 (   2)    3.4127732919    10.2586355896   -30.9999555547    -5.6414488172
   9 (  16)   -8.0030283504    -8.3212414840     2.0104983118   -20.7200709104
  10 (   9)   13.2364948412     1.8767883671    17.5130677056    -5.0313101030
  11 (  23)    1.0018763035   -12.2220878772     1.3671011496     0.7498989339
  12 (   3)    9.4157637378     2.9550563366    -6.8274600541   -10.0985523545
  13 (  17)  -25.6563211655    16.0852626944   -31.4472791768    -9.5901194084
  14 (  10)    2.2074782914   -13.5663212892    -3.5460529665    29.9250169837
  15 (  24)   11.5587101862    -6.4660944566     5.8409045449     5.9151760968
  16 (   4)    5.8409045449     5.9151760968     9.4157637378     2.9550563366
  17 (  18)   -3.5460529665    29.9250169837   -25.6563211655    16.0852626944
  18 (  11)  -31.4472791768    -9.5901194084    11.5587101862    -6.4660944566
  19 (  25)   -6.8274600541   -10.0985523545     2.2074782914   -13.5663212892
  20 (   5)    2.0104983118   -20.7200709104    13.2364948412     1.8767883671
  21 (  19)  -30.9999555547    -5.6414488172     1.0018763035   -12.2220878772
  22 (  12)   17.5130677056    -5.0313101030     3.4127732919    10.2586355896
  23 (  26)    1.3671011496     0.7498989339    -8.0030283504    -8.3212414840
  24 (   6)    9.9796188501    -9.6408663259     2.8294700831     0.5106027780
  25 (  20)   -7.6241120460    -1.5581921548   -11.6595587366   -19.7679307866
  26 (  13)   -7.6577999302     8.5976313922    46.1480382688    -2.2866689587
  27 (  27)   32.4689557080    -0.3465505441   -17.5691832932     0.3833862984
*/

/***************/

void radix28_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-28 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix7_dif_pass for details on the radix-15 subtransforms.
*/
	int j,j1,j2;
	static int n28,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27, first_entry=TRUE;
	static double	uc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					us1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					uc2 =-.22252093395631440426,	 /* cos(2u)	*/
					us2 = .97492791218182360702,	 /* sin(2u)	*/
					uc3 =-.90096886790241912622,	 /* cos(3u)	*/
					us3 = .43388373911755812050;	 /* sin(3u)	*/
	double re,im,rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i;

	if(!first_entry && (n/28) != n28)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n28=n/28;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n28;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p10 = p09 + p01;
		p11 = p10 + p01;
		p12 = p11 + p01;
		p13 = p12 + p01;
		p14 = p13 + p01;
		p15 = p14 + p01;
		p16 = p15 + p01;
		p17 = p16 + p01;
		p18 = p17 + p01;
		p19 = p18 + p01;
		p20 = p19 + p01;
		p21 = p20 + p01;
		p22 = p21 + p01;
		p23 = p22 + p01;
		p24 = p23 + p01;
		p25 = p24 + p01;
		p26 = p25 + p01;
		p27 = p26 + p01;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p17 = p17 + ( (p17 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p19 = p19 + ( (p19 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p21 = p21 + ( (p21 >> DAT_BITS) << PAD_BITS );
		p22 = p22 + ( (p22 >> DAT_BITS) << PAD_BITS );
		p23 = p23 + ( (p23 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p25 = p25 + ( (p25 >> DAT_BITS) << PAD_BITS );
		p26 = p26 + ( (p26 >> DAT_BITS) << PAD_BITS );
		p27 = p27 + ( (p27 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-28 pass is here.	*/

	for(j=0; j < n28; j += 2)
	{
	#if defined(USE_SSE2)
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*
		Twiddleless version requires us to swap inputs as follows:
		indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27
			  -> 0,24,20,16,12, 8, 4,21,17,13, 9, 5, 1,25,14,10, 6, 2,26,22,18, 7, 3,27,23,19,15,11

		I.e. start out with first quartet of indices {0,7,14,21}, permute those according to
		  {0,7,14,21}*27%28 = {0,21,14,7}, then each is head of a length-7 list of indices with decrement 4

		Remember, inputs to DIT are bit-reversed, so
		a[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27] contain
		x[0,14, 7,21, 1,15, 8,22, 2,16, 9,23, 3,17,10,24, 4,18,11,25, 5,19,12,26, 6,20,13,27], which get swapped to
		x[0,14,21, 7,24,10,17, 3,20, 6,13,27,16, 2, 9,23,12,26, 5,19, 8,22, 1,15, 4,18,25,11], which means the a-indices get swapped as
		a[0, 1, 3, 2,15,14,13,12,25,24,26,27, 9, 8,10,11,22,23,20,21, 6, 7, 4, 5,16,17,19,18].
	*/
	/* EWM: 10/18/04: We swap the odd-index outputs of each of the radix-4 DIT transforms (1<=>3, 5<=>7, etc.) so that the indexing
					  of all the radix-7 transforms (really just the 2nd and 4th of these) winds up being in-place. This allows us
					  to properly re-use the ajp1 variables in the carry-pass version of this routine.
					  (NOTE: Didn't bother to do the swapping in the macro-less version of this routine - would only need to
					   do the swapping there if for some reason we wanted the macro-less DIT code in the carry routines.)
	*/
	/*...gather the needed data (28 64-bit complex, i.e. 56 64-bit reals) and do 7 radix-4 transforms...*/
					/*                                    inputs                                  */ /*                         outputs                   */
		RADIX_04_DIT(a[j1    ],a[j2    ],a[j1+p01],a[j2+p01],a[j1+p03],a[j2+p03],a[j1+p02],a[j2+p02],a1p00r,a1p00i,a1p03r,a1p03i,a1p02r,a1p02i,a1p01r,a1p01i,rt,it);
		RADIX_04_DIT(a[j1+p15],a[j2+p15],a[j1+p14],a[j2+p14],a[j1+p13],a[j2+p13],a[j1+p12],a[j2+p12],a1p04r,a1p04i,a1p07r,a1p07i,a1p06r,a1p06i,a1p05r,a1p05i,rt,it);
		RADIX_04_DIT(a[j1+p25],a[j2+p25],a[j1+p24],a[j2+p24],a[j1+p26],a[j2+p26],a[j1+p27],a[j2+p27],a1p08r,a1p08i,a1p11r,a1p11i,a1p10r,a1p10i,a1p09r,a1p09i,rt,it);
		RADIX_04_DIT(a[j1+p09],a[j2+p09],a[j1+p08],a[j2+p08],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a1p12r,a1p12i,a1p15r,a1p15i,a1p14r,a1p14i,a1p13r,a1p13i,rt,it);
		RADIX_04_DIT(a[j1+p22],a[j2+p22],a[j1+p23],a[j2+p23],a[j1+p20],a[j2+p20],a[j1+p21],a[j2+p21],a1p16r,a1p16i,a1p19r,a1p19i,a1p18r,a1p18i,a1p17r,a1p17i,rt,it);
		RADIX_04_DIT(a[j1+p06],a[j2+p06],a[j1+p07],a[j2+p07],a[j1+p04],a[j2+p04],a[j1+p05],a[j2+p05],a1p20r,a1p20i,a1p23r,a1p23i,a1p22r,a1p22i,a1p21r,a1p21i,rt,it);
		RADIX_04_DIT(a[j1+p16],a[j2+p16],a[j1+p17],a[j2+p17],a[j1+p19],a[j2+p19],a[j1+p18],a[j2+p18],a1p24r,a1p24i,a1p27r,a1p27i,a1p26r,a1p26i,a1p25r,a1p25i,rt,it);

	/*...and now do 4 radix-4 transforms...*/
					/*                                              inputs                                          */ /*               intermediates              */ /*                                                                     outputs                                                           */
		RADIX_07_DFT(a1p00r,a1p00i,a1p04r,a1p04i,a1p08r,a1p08i,a1p12r,a1p12i,a1p16r,a1p16i,a1p20r,a1p20i,a1p24r,a1p24i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a[j1    ],a[j2    ],a[j1+p08],a[j2+p08],a[j1+p16],a[j2+p16],a[j1+p24],a[j2+p24],a[j1+p04],a[j2+p04],a[j1+p12],a[j2+p12],a[j1+p20],a[j2+p20],uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a1p03r,a1p03i,a1p07r,a1p07i,a1p11r,a1p11i,a1p15r,a1p15i,a1p19r,a1p19i,a1p23r,a1p23i,a1p27r,a1p27i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a[j1+p07],a[j2+p07],a[j1+p15],a[j2+p15],a[j1+p23],a[j2+p23],a[j1+p03],a[j2+p03],a[j1+p11],a[j2+p11],a[j1+p19],a[j2+p19],a[j1+p27],a[j2+p27],uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a1p02r,a1p02i,a1p06r,a1p06i,a1p10r,a1p10i,a1p14r,a1p14i,a1p18r,a1p18i,a1p22r,a1p22i,a1p26r,a1p26i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a[j1+p14],a[j2+p14],a[j1+p22],a[j2+p22],a[j1+p02],a[j2+p02],a[j1+p10],a[j2+p10],a[j1+p18],a[j2+p18],a[j1+p26],a[j2+p26],a[j1+p06],a[j2+p06],uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
		RADIX_07_DFT(a1p01r,a1p01i,a1p05r,a1p05i,a1p09r,a1p09i,a1p13r,a1p13i,a1p17r,a1p17i,a1p21r,a1p21i,a1p25r,a1p25i,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,a[j1+p21],a[j2+p21],a[j1+p01],a[j2+p01],a[j1+p09],a[j2+p09],a[j1+p17],a[j2+p17],a[j1+p25],a[j2+p25],a[j1+p05],a[j2+p05],a[j1+p13],a[j2+p13],uc1,us1,uc2,us2,uc3,us3,rt,it,re,im);
	}
}

