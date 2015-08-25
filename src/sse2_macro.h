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

/****************************************************************************
 * We now include this header file if it was not included before.
 ****************************************************************************/
#ifndef sse2_macro_included
#define sse2_macro_included

#ifdef USE_SSE2

	#ifdef COMPILER_TYPE_MSVC

		/* ...Special 2-register, Karatsuba-style CMUL for the last 4 of each radix-4-with-twiddles subtransform.
		   Address of head element of the 4 neighboring main-array elts in addr [assumed already loaded into one of eax,abx,ecx,edx],
		   sincos data triplet starts at cddr [similarly loaded into GPR], result returned as [re,im] in registers [r0,r1], respectively.

		"Normal" complex MUL:

			(a+bI)*(c+dI) = x+yI,

			x = ac-bd
			y = ad+bc, 2 ADD, 4 MUL.

		Define w, z such that

			x = v+w
			y = v-w

		so

			v = (x+y)/2 = [ac+ad+bc-bd]/2 = [1] - bd,
			w = (x-y)/2 = [ac-ad-bc-bd]/2 = [4] - bd.

		So far, so good. The problem is that because of the sign-flip on the bd term in the x-output of CMUL,
		no simple Karatsuba-style combo of form (a+-b)*(c+-d) gives the right sign combination we need:

		1.	(a+b)*(c+d)/2 = [ac+ad+bc+bd]/2
		2.	(a+b)*(c-d)/2 = [ac-ad+bc-bd]/2
		3.	(a-b)*(c+d)/2 = [ac+ad-bc-bd]/2
		4.	(a-b)*(c-d)/2 = [ac-ad-bc+bd]/2

		E.g. 1+-4 gives the proper ac,ad,bc terms, but each needs the bd term flipped. If we e.g. change the sign on the "d" we fix the bd terms but the ad term is now flipped.

		In terms of [1-4], we have

			v = [ac+ad+bc-bd]/2 = [1] - bd,
			w = [ac-ad-bc-bd]/2 = [4] - bd.

		and

			x = v+w = [1]+[4]-2*bd,	 [Can also do as x = y + 2*([4]-bd), if that proves favorable for some reason.]
			y = v-w = [1]-[4] .

		Let's say (c+d)/2,(c-d)/2,d are in memlocs mem0,1,2, resp, and a,b are in memlocs [a],[b].

		Then [1]:could be done in SSE2 as shown in the macro below.
		Needs just 2 registers, plus one intermediate store. Total cost:

			4 MEM [+9 implicit], 6 ADD, 3 MUL        [variant A]
			4 MEM [+7 implicit], 7 ADD, 3 MUL        [variant B]
					 [or 6 ADD, 4 MUL].
		*/
		#define SSE2_CMUL_2REG_A(__addr,__cddr,__r0,__r1)\
		{\
			/* [variant A]:														*/	\
			__asm	movaps	__r1,[__addr+0x10]	/* b							*/	\
			__asm	movaps	__r0,[__addr     ]	/* a							*/	\
			__asm	mulpd	__r1,[__cddr+0x20]	/* b*d							*/	\
			__asm	movaps	[__cddr+0x30],__r1	/* store b*d to tmp	[variant B]:*/	\
			__asm	movaps	__r1,[__addr     ]	/* a							*/	\
			__asm	subpd	__r0,[__addr+0x10]	/* r0: a-b						*/	\
			__asm	addpd	__r1,[__addr+0x10]	/* r1: a+b						*/	\
			__asm	mulpd	__r1,[__cddr     ]	/* [1]							*/	\
			__asm	mulpd	__r0,[__cddr+0x10]	/* [4]							*/	\
			__asm	subpd	__r1,__r0			/* y = [1]-[4]					*/	\
			__asm	subpd	__r0,[__cddr+0x30]	/* [4]-b*d						*/	\
			__asm	addpd	__r0,__r0			/* 2*([4]-b*d)					*/	\
			__asm	addpd	__r0,__r1			/* x = y + 2*([4]-b*d).			*/	\
		}

		/* Complex multiply of 2 roots of unity - use e.g. for "multiply up" of sincos twiddles. */\
		#define SSE2_CMUL_EXPO(__cA,__cB,__cAmB,__cApB)\
		{\
			__asm	mov	eax, __cA	\
			__asm	mov	ebx, __cB	\
			__asm	mov	ecx, __cAmB	\
			__asm	mov	edx, __cApB	\
			\
			__asm	movaps	xmm0,[eax     ]	/* __cA */		__asm	movaps	xmm4,[ebx     ]	/* __cB */	\
			__asm	movaps	xmm2,[eax+0x10]	/* __sA */		__asm	movaps	xmm5,[ebx+0x10]	/* __sB */	\
			__asm	movaps	xmm1,xmm0	/* cpy __cA */		\
			__asm	movaps	xmm3,xmm2	/* cpy __sA */		\
			\
			__asm	mulpd	xmm0,xmm4	/* t1 = __cA*__cB */\
			__asm	mulpd	xmm1,xmm5	/* t2 = __cA*__sB */\
			__asm	mulpd	xmm2,xmm4	/* rt = __sA*__cB */\
			__asm	mulpd	xmm3,xmm5	/* it = __sA*__sB */\
			__asm	movaps	xmm4,xmm0	/* t1 copy */		\
			__asm	movaps	xmm5,xmm1	/* t2 copy */		\
			__asm	addpd	xmm0,xmm3	/* __cA-B */		\
			__asm	subpd	xmm1,xmm2	/* __sA-B */		\
			__asm	subpd	xmm4,xmm3	/* __cA+B */		\
			__asm	addpd	xmm5,xmm2	/* __sA+B */		\
			__asm	movaps	[ecx     ],xmm0	/* __cA-B */	\
			__asm	movaps	[ecx+0x10],xmm1	/* __sA-B */	\
			__asm	movaps	[edx     ],xmm4	/* __cA+B */	\
			__asm	movaps	[edx+0x10],xmm5	/* __sA+B */	\
		}

		/*...Radix-3 DFT: Inputs in memory locations __i0,__i1,__i2,\
				outputs go into memory locations __o0,__o1,__o2, possibly coincident with inputs:\
		*/\
		#define SSE2_RADIX_03_DFT(__i0,__i1,__i2, __c3m1, __o0,__o1,__o2)\
		{\
		/*\
			t2r =A1r - A2r;			t2i =A1i - A2i;	\
			t1r =A1r + A2r;			t1i =A1i + A2i;	\
			B0r =A0r + t1r;			B0i =A0i + t1i;	\
			B1r =B0r + t1r*c3m1;	B1i =B0i + t1i*c3m1;\
			rt  =c*t2r;				it  =c*t2i;		\
			B2r =B1r + it;			B2i =B1i - rt;	\
			B1r =B1r - it;			B1i =B1i + rt;	\
		*/\
			__asm	mov	eax, __i0\
			__asm	mov	ebx, __i1\
			__asm	mov	ecx, __i2\
			__asm	mov	edx, __c3m1\
			\
			__asm	movaps	xmm2,[ebx     ]	/* A1r */\
			__asm	movaps	xmm3,[ebx+0x10]	/* A1i */\
			__asm	movaps	xmm0,[eax     ]	/* A0r */\
			__asm	movaps	xmm1,[eax+0x10]	/* A0i */\
			__asm	movaps	xmm6,[ecx     ]	/* A2r */\
			__asm	movaps	xmm7,[ecx+0x10]	/* A2i */\
			__asm	movaps	xmm4,xmm2		/* t2r: init to copy of A1r */\
			__asm	movaps	xmm5,xmm3		/* t2i: init to copy of A1i */\
			\
			__asm	mov	eax, __o0\
			__asm	mov	ebx, __o1\
			__asm	mov	ecx, __o2\
			__asm	addpd	xmm2,xmm6		/* t1r=A1r+A2r */\
			__asm	addpd	xmm3,xmm7		/* t1i=A1i+A2i */\
			__asm	subpd	xmm4,xmm6		/* t2r=A1r-A2r */\
			__asm	subpd	xmm5,xmm7		/* t2i=A1i-A2i */\
			__asm	addpd	xmm0,xmm2	/* B0r = A0r+t1r */\
			__asm	addpd	xmm1,xmm3	/* B0i = A0i+t1i */\
			__asm	movaps	xmm6,[edx     ]	/* c3m1*/\
			__asm	movaps	xmm7,[edx+0x10]	/* c   */\
			__asm	movaps	[eax      ],xmm0	/* <- B0r */\
			__asm	movaps	[eax+0x010],xmm1	/* <- B0i */\
			\
			__asm	mulpd	xmm2,xmm6		/* t1r *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t1i *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = c*t2r   */\
			__asm	mulpd	xmm5,xmm7		/* it = c*t2i   */\
			__asm	addpd	xmm2,xmm0	/* B1r = B0r+c3m1*t1r */\
			__asm	addpd	xmm3,xmm1	/* B1i = B0i+c3m1*t1i */\
			\
			__asm	movaps	xmm0,xmm2	/* copy B1r */\
			__asm	movaps	xmm1,xmm3	/* copy B1i */\
			\
			__asm	subpd	xmm2,xmm5	/* B1r = B1r-it */\
			__asm	addpd	xmm3,xmm4	/* B1i = B1i+rt */\
			__asm	addpd	xmm0,xmm5	/* B2r = B1r+it */\
			__asm	subpd	xmm1,xmm4	/* B2i = B1i-rt */\
			\
			__asm	movaps	[ebx      ],xmm2	/* <- B1r */\
			__asm	movaps	[ebx+0x010],xmm3	/* <- B1i */\
			__asm	movaps	[ecx      ],xmm0	/* <- B2r */\
			__asm	movaps	[ecx+0x010],xmm1	/* <- B2i */\
		}

		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		#define SSE2_RADIX4_DIT_0TWIDDLE(__add0, __add1, __add2, __add3, __tmp)\
		{\
			__asm	mov	eax, __add0\
			__asm	mov	ebx, __add1\
			__asm	mov	ecx, __add2\
			__asm	mov	edx, __add3\
			\
			__asm	movaps	xmm2,[eax     ]	/* t1 */				__asm	movaps	xmm6,[ecx     ]	/* t5 */\
			__asm	movaps	xmm3,[eax+0x10]	/* t2 */				__asm	movaps	xmm7,[ecx+0x10]	/* t6 */\
			__asm	movaps	xmm0,[ebx     ]	/* t3 */				__asm	movaps	xmm4,[edx     ]	/* t7 */\
			__asm	movaps	xmm1,[ebx+0x10]	/* t4 */				__asm	movaps	xmm5,[edx+0x10]	/* t8 */\
			\
			__asm	subpd	xmm2,xmm0		/* ~t3 = t1-t3 */		__asm	subpd	xmm6,xmm4		/* ~t7 = t5-t7 */\
			__asm	subpd	xmm3,xmm1		/* ~t4 = t2-t4 */		__asm	subpd	xmm7,xmm5		/* ~t8 = t6-t8 */\
			__asm	addpd	xmm0,xmm0		/*        2*t3 */		__asm	addpd	xmm4,xmm4		/*        2*t7 */\
			__asm	addpd	xmm1,xmm1		/*        2*t4 */		__asm	addpd	xmm5,xmm5		/*        2*t8 */\
			__asm	addpd	xmm0,xmm2		/* ~t1 = t1+t3 */		__asm	addpd	xmm4,xmm6		/* ~t5 = t5+t7 */\
			__asm	addpd	xmm1,xmm3		/* ~t2 = t2+t4 */		__asm	addpd	xmm5,xmm7		/* ~t6 = t6+t8 */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			__asm	mov	eax, __tmp	\
			__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */			__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */\
			__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */			__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */\
			__asm	movaps	[eax+0x040],xmm0	/* <- ~t5 */		__asm	movaps	[eax+0x060],xmm2	/* <- ~t7 */\
			__asm	movaps	[eax+0x050],xmm1	/* <- ~t6 */		__asm	movaps	[eax+0x030],xmm3	/* <- ~t4 */\
			__asm	addpd	xmm4,xmm4	/*          2*t5 */			__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm5,xmm5	/*          2*t6 */			__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */			__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */\
			__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */			__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */\
			__asm	movaps	[eax      ],xmm4	/* <- ~t1 */		__asm	movaps	[eax+0x020],xmm7	/* <- ~t3 */\
			__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */		__asm	movaps	[eax+0x070],xmm6	/* <- ~t8 */\
		}

		/* Full-inline-asm version of SSE2_RADIX4_DIT_0TWIDDLE. Assumes base addresses add0-3 in e[a-d]x,
		base offset in edi, these must remain unchanged, so esi is available for temporary storage. */
		#define SSE2_RADIX4_DIT_0TWIDDLE_B(__tmp)\
		{\
			__asm	movaps	xmm2,[eax     ]	/* t1 */				__asm	movaps	xmm6,[ecx     ]	/* t5 */\
			__asm	movaps	xmm3,[eax+0x10]	/* t2 */				__asm	movaps	xmm7,[ecx+0x10]	/* t6 */\
			__asm	movaps	xmm0,[ebx     ]	/* t3 */				__asm	movaps	xmm4,[edx     ]	/* t7 */\
			__asm	movaps	xmm1,[ebx+0x10]	/* t4 */				__asm	movaps	xmm5,[edx+0x10]	/* t8 */\
			\
			__asm	subpd	xmm2,xmm0		/* ~t3 = t1-t3 */		__asm	subpd	xmm6,xmm4		/* ~t7 = t5-t7 */\
			__asm	subpd	xmm3,xmm1		/* ~t4 = t2-t4 */		__asm	subpd	xmm7,xmm5		/* ~t8 = t6-t8 */\
			__asm	addpd	xmm0,xmm0		/*        2*t3 */		__asm	addpd	xmm4,xmm4		/*        2*t7 */\
			__asm	addpd	xmm1,xmm1		/*        2*t4 */		__asm	addpd	xmm5,xmm5		/*        2*t8 */\
			__asm	addpd	xmm0,xmm2		/* ~t1 = t1+t3 */		__asm	addpd	xmm4,xmm6		/* ~t5 = t5+t7 */\
			__asm	addpd	xmm1,xmm3		/* ~t2 = t2+t4 */		__asm	addpd	xmm5,xmm7		/* ~t6 = t6+t8 */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			__asm	mov	esi, __tmp	\
			__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */			__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */\
			__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */			__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */\
			__asm	movaps	[esi+0x040],xmm0	/* <- ~t5 */		__asm	movaps	[esi+0x060],xmm2	/* <- ~t7 */\
			__asm	movaps	[esi+0x050],xmm1	/* <- ~t6 */		__asm	movaps	[esi+0x030],xmm3	/* <- ~t4 */\
			__asm	addpd	xmm4,xmm4	/*          2*t5 */			__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm5,xmm5	/*          2*t6 */			__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */			__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */\
			__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */			__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */\
			__asm	movaps	[esi      ],xmm4	/* <- ~t1 */		__asm	movaps	[esi+0x020],xmm7	/* <- ~t3 */\
			__asm	movaps	[esi+0x010],xmm5	/* <- ~t2 */		__asm	movaps	[esi+0x070],xmm6	/* <- ~t8 */\
		}

		/* DIT radix-4 subconvolution, sans twiddles, inputs in add0,1,2,3, outputs written 16*__stride bytes apart: */
		#define SSE2_RADIX4_DIT_0TWIDDLE_STRIDE(__add0, __add1, __add2, __add3, __tmp, __stride)\
		{\
			__asm	mov	eax, __add0\
			__asm	mov	ebx, __add1\
			__asm	mov	ecx, __add2\
			__asm	mov	edx, __add3\
			\
			__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */				__asm	movaps	xmm4,[ecx     ]	/* a[jt+p2] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p2] */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[jt   ] */	__asm	movaps	xmm6,xmm4		/* xmm4 <- cpy a[jt+p2] */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[jp   ] */	__asm	movaps	xmm7,xmm5		/* xmm5 <- cpy a[jp+p2] */\
			\
			__asm	mov	eax, __tmp		/* addout_0 */\
			__asm	mov	ecx, __stride\
			\
			__asm	addpd	xmm0,[ebx     ]	/* t1 */					__asm	addpd	xmm4,[edx     ]	/* t5 */\
			__asm	addpd	xmm1,[ebx+0x10]	/* t2 */					__asm	addpd	xmm5,[edx+0x10]	/* t6 */\
			__asm	subpd	xmm2,[ebx     ]	/* t3 */					__asm	subpd	xmm6,[edx     ]	/* t7 */\
			__asm	subpd	xmm3,[ebx+0x10]	/* t4 */					__asm	subpd	xmm7,[edx+0x10]	/* t8 */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			__asm	mov	ebx, eax	/* addout_0; need to add 1*stride */\
			__asm	add	ebx, ecx	/* addout_1 */\
			__asm	mov	edx, ebx	/* addout_1; need to add 2*stride */\
			__asm	add	ecx, ecx	/* 2*stride */\
			__asm	add	edx, ecx	/* addout_3 */\
			__asm	add	ecx, eax	/* addout_2 */\
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

		/* Variant of SSE2_RADIX4_DIT_0TWIDDLE_STRIDE, which assumes __add0,1,2,3 enter in eax,ebx,ecx,edx,
		and that the contents of eax and esi must not be changed on return.
		*/
		#define SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_B(__tmp, __stride)\
		{\
			__asm	movaps	xmm2,[eax     ]	/* a[jt   ] */				__asm	movaps	xmm6,[ecx     ]	/* a[jt+p2] */\
			__asm	movaps	xmm3,[eax+0x10]	/* a[jp   ] */				__asm	movaps	xmm7,[ecx+0x10]	/* a[jp+p2] */\
			__asm	movaps	xmm0,[ebx     ]	/* a[jt+p1] */				__asm	movaps	xmm4,[edx     ]	/* a[jt+p3] */\
			__asm	movaps	xmm1,[ebx+0x10]	/* a[jp+p1] */				__asm	movaps	xmm5,[edx+0x10]	/* a[jp+p3] */\
			\
			__asm	mov	edi, __tmp		/* addout_0 */\
			__asm	mov	ecx, __stride\
			\
			__asm	subpd	xmm2,xmm0		/*~t3 = t1-t3 */			__asm	subpd	xmm6,xmm4		/* ~t7 = t5-t7 */\
			__asm	subpd	xmm3,xmm1		/*~t4 = t2-t4 */			__asm	subpd	xmm7,xmm5		/* ~t8 = t6-t8 */\
			__asm	addpd	xmm0,xmm0		/*       2*t3 */			__asm	addpd	xmm4,xmm4		/*        2*t7 */\
			__asm	addpd	xmm1,xmm1		/*       2*t4 */			__asm	addpd	xmm5,xmm5		/*        2*t8 */\
			__asm	addpd	xmm0,xmm2		/*~t1 = t1+t3 */			__asm	addpd	xmm4,xmm6		/* ~t5 = t5+t7 */\
			__asm	addpd	xmm1,xmm3		/*~t2 = t2+t4 */			__asm	addpd	xmm5,xmm7		/* ~t6 = t6+t8 */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			__asm	mov	ebx, edi	/* addout_0; need to add 1*stride */\
			__asm	add	ebx, ecx	/* addout_1 */\
			__asm	mov	edx, ebx	/* addout_1; need to add 2*stride */\
			__asm	add	ecx, ecx	/* 2*stride */\
			__asm	add	edx, ecx	/* addout_3 */\
			__asm	add	ecx, edi	/* addout_2 */\
			\
			__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */				__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */\
			__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */				__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */\
			__asm	movaps	[ecx      ],xmm0	/* <- ~t5 */			__asm	movaps	[edx      ],xmm2	/* <- ~t7 */\
			__asm	movaps	[ecx+0x010],xmm1	/* <- ~t6 */			__asm	movaps	[ebx+0x010],xmm3	/* <- ~t4 */\
			__asm	addpd	xmm4,xmm4	/*          2*t5 */				__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm5,xmm5	/*          2*t6 */				__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */				__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */\
			__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */				__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */\
			__asm	movaps	[edi      ],xmm4	/* <- ~t1 */			__asm	movaps	[ebx      ],xmm7	/* <- ~t3 */\
			__asm	movaps	[edi+0x010],xmm5	/* <- ~t2 */			__asm	movaps	[edx+0x010],xmm6	/* <- ~t8 */\
		}

		/* Variant of SSE2_RADIX4_DIT_0TWIDDLE_STRIDE, which assumes __add0,1,2,3 enter in in0,in1,in2,in3,
		and that the contents of these 4 register aliases must not be changed on return.
		*/
		#define SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(in0,in1,in2,in3, __stride1,__stride2, __out)\
		{\
			__asm	movaps	xmm2,[in0     ]	/* a[jt   ] */				__asm	movaps	xmm6,[in2     ]	/* a[jt+p2] */\
			__asm	movaps	xmm3,[in0+0x10]	/* a[jp   ] */				__asm	movaps	xmm7,[in2+0x10]	/* a[jp+p2] */\
			__asm	movaps	xmm0,[in1     ]	/* a[jt+p1] */				__asm	movaps	xmm4,[in3     ]	/* a[jt+p3] */\
			__asm	movaps	xmm1,[in1+0x10]	/* a[jp+p1] */				__asm	movaps	xmm5,[in3+0x10]	/* a[jp+p3] */\
			\
			__asm	mov	esi, __out		/* addout_0 */\
			__asm	mov	edi, __stride1\
			__asm	add	edi, esi		/* addout_1 */\
			\
			__asm	subpd	xmm2,xmm0		/*~t3 = t1-t3 */			__asm	subpd	xmm6,xmm4		/* ~t7 = t5-t7 */\
			__asm	subpd	xmm3,xmm1		/*~t4 = t2-t4 */			__asm	subpd	xmm7,xmm5		/* ~t8 = t6-t8 */\
			__asm	addpd	xmm0,xmm0		/*       2*t3 */			__asm	addpd	xmm4,xmm4		/*        2*t7 */\
			__asm	addpd	xmm1,xmm1		/*       2*t4 */			__asm	addpd	xmm5,xmm5		/*        2*t8 */\
			__asm	addpd	xmm0,xmm2		/*~t1 = t1+t3 */			__asm	addpd	xmm4,xmm6		/* ~t5 = t5+t7 */\
			__asm	addpd	xmm1,xmm3		/*~t2 = t2+t4 */			__asm	addpd	xmm5,xmm7		/* ~t6 = t6+t8 */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			\
			__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */				__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */\
			__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */				__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */\
			__asm	movaps	[esi+__stride2      ],xmm0	/* <- ~t5 */	__asm	movaps	[edi+__stride2      ],xmm2	/* <- ~t7 */\
			__asm	movaps	[esi+__stride2+0x010],xmm1	/* <- ~t6 */	__asm	movaps	[edi          +0x010],xmm3	/* <- ~t4 */\
			__asm	addpd	xmm4,xmm4	/*          2*t5 */				__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm5,xmm5	/*          2*t6 */				__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */				__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */\
			__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */				__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */\
			__asm	movaps	[esi      ],xmm4	/* <- ~t1 */			__asm	movaps	[edi                ],xmm7	/* <- ~t3 */\
			__asm	movaps	[esi+0x010],xmm5	/* <- ~t2 */			__asm	movaps	[edi+__stride2+0x010],xmm6	/* <- ~t8 */\
		}

		/* In-place Variant of SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C. We use separate pseudonyms for the in-and-out-address
		registers to allow for any desired index permutation here: */
		#define SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_D(in0,in1,in2,in3, out0,out1,out2,out3)\
		{\
			__asm	movaps	xmm2,[in0     ]	/* a[jt   ] */				__asm	movaps	xmm6,[in2     ]	/* a[jt+p2] */\
			__asm	movaps	xmm3,[in0+0x10]	/* a[jp   ] */				__asm	movaps	xmm7,[in2+0x10]	/* a[jp+p2] */\
			__asm	movaps	xmm0,[in1     ]	/* a[jt+p1] */				__asm	movaps	xmm4,[in3     ]	/* a[jt+p3] */\
			__asm	movaps	xmm1,[in1+0x10]	/* a[jp+p1] */				__asm	movaps	xmm5,[in3+0x10]	/* a[jp+p3] */\
			\
			__asm	subpd	xmm2,xmm0		/*~t3 = t1-t3 */			__asm	subpd	xmm6,xmm4		/* ~t7 = t5-t7 */\
			__asm	subpd	xmm3,xmm1		/*~t4 = t2-t4 */			__asm	subpd	xmm7,xmm5		/* ~t8 = t6-t8 */\
			__asm	addpd	xmm0,xmm0		/*       2*t3 */			__asm	addpd	xmm4,xmm4		/*        2*t7 */\
			__asm	addpd	xmm1,xmm1		/*       2*t4 */			__asm	addpd	xmm5,xmm5		/*        2*t8 */\
			__asm	addpd	xmm0,xmm2		/*~t1 = t1+t3 */			__asm	addpd	xmm4,xmm6		/* ~t5 = t5+t7 */\
			__asm	addpd	xmm1,xmm3		/*~t2 = t2+t4 */			__asm	addpd	xmm5,xmm7		/* ~t6 = t6+t8 */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			\
			__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */				__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */\
			__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */				__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */\
			__asm	movaps	[out2      ],xmm0	/* <- ~t5 */			__asm	movaps	[out3      ],xmm2	/* <- ~t7 */\
			__asm	movaps	[out2+0x010],xmm1	/* <- ~t6 */			__asm	movaps	[out1+0x010],xmm3	/* <- ~t4 */\
			__asm	addpd	xmm4,xmm4	/*          2*t5 */				__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm5,xmm5	/*          2*t6 */				__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */				__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */\
			__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */				__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */\
			__asm	movaps	[out0      ],xmm4	/* <- ~t1 */			__asm	movaps	[out1      ],xmm7	/* <- ~t3 */\
			__asm	movaps	[out0+0x010],xmm5	/* <- ~t2 */			__asm	movaps	[out3+0x010],xmm6	/* <- ~t8 */\
		}

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

		/*******************************************************************************************************/
		/***************************** DIF versions of the radix-4 DFT macros: *********************************/
		/*******************************************************************************************************/

		/* DIF radix-4 subconvolution, sans twiddles, inputs 16*__stride bytes apart, outputs written to add0,1,2,3: */
		#define SSE2_RADIX4_DIF_0TWIDDLE_STRIDE(__add0, __add1, __add2, __add3, __tmp, __stride)\
		{\
			__asm	mov	eax, __tmp	\
			__asm	mov	esi, __stride\
			__asm	mov	ebx, eax\
			__asm	add	ebx, esi	/* add_in1 */\
			__asm	shl	esi, 1		/* stride*2 */\
			\
			__asm	movaps	xmm0,[eax     ]	/* t1 */					__asm	movaps	xmm2,[ebx     ]	/* t5 */\
			__asm	movaps	xmm1,[eax+0x10]	/* t2 */					__asm	movaps	xmm3,[ebx+0x10]	/* t6 */\
			__asm	movaps	xmm4,[eax     ]	/* xmm2 <- cpy t1 */		__asm	movaps	xmm6,[ebx     ]	/* xmm2 <- cpy t5 */\
			__asm	movaps	xmm5,[eax+0x10]	/* xmm3 <- cpy t2 */		__asm	movaps	xmm7,[ebx+0x10]	/* xmm3 <- cpy t6 */\
			\
			__asm	add	eax, esi	/* add_in2 */\
			__asm	add	ebx, esi	/* add_in3 */\
			\
			__asm	addpd	xmm0,[eax     ]	/* ~t1 = t1 + t5 */			__asm	addpd	xmm2,[ebx     ]	/* ~t3 = t3 + t7 */\
			__asm	addpd	xmm1,[eax+0x10]	/* ~t2 = t2 + t6 */			__asm	addpd	xmm3,[ebx+0x10]	/* ~t4 = t4 + t8 */\
			__asm	subpd	xmm4,[eax     ]	/* ~t5 = t1 - t5 */			__asm	subpd	xmm6,[ebx     ]	/* ~t7 = t3 - t7 */\
			__asm	subpd	xmm5,[eax+0x10]	/* ~t6 = t2 - t6 */			__asm	subpd	xmm7,[ebx+0x10]	/* ~t8 = t4 - t8 */\
			\
			/* Finish radix-4 butterfly and store results into main-array slots: */\
			__asm	mov	eax, __add0\
			__asm	mov	ebx, __add1\
			__asm	mov	ecx, __add2\
			__asm	mov	edx, __add3\
			\
			__asm	subpd	xmm0,xmm2	/* ~t3 <- t1 -t3 */				__asm	subpd	xmm4,xmm7	/* ~t5 <- t5 -t8 */\
			__asm	subpd	xmm1,xmm3	/* ~t4 <- t2 -t4 */				__asm	subpd	xmm5,xmm6	/* ~t8 <- t6 -t7 */\
			__asm	movaps	[ebx      ],xmm0	/* <- ~t3 */			__asm	movaps	[ecx      ],xmm4	/* <- ~t5 */\
			__asm	movaps	[ebx+0x010],xmm1	/* <- ~t4 */			__asm	movaps	[edx+0x010],xmm5	/* <- ~t8 */\
			__asm	addpd	xmm2,xmm2	/*          2*t3 */				__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm3,xmm3	/*          2*t4 */				__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm2,xmm0	/* ~t1 <- t1 +t3 */				__asm	addpd	xmm7,xmm4	/* ~t7 <- t5 +t8 */\
			__asm	addpd	xmm3,xmm1	/* ~t2 <- t2 +t4 */				__asm	addpd	xmm6,xmm5	/* ~t6 <- t6 +t7 */\
			__asm	movaps	[eax      ],xmm2	/* <- ~t1 */			__asm	movaps	[edx      ],xmm7	/* <- ~t7 */\
			__asm	movaps	[eax+0x010],xmm3	/* <- ~t2 */			__asm	movaps	[ecx+0x010],xmm6	/* <- ~t6 */\
		}

		/* Variant of SSE2_RADIX4_DIF_0TWIDDLE_STRIDE, which
		Assumes __add0,1,2,3 enter in eax,ebx,ecx,edx, and that the contents of eax must not be changed on return.
		*/
		#define SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_B(__tmp, __stride)\
		{\
			__asm	mov	edi, __tmp	\
			__asm	mov	esi, __stride\
			__asm	add	esi, edi	/* add_in1 */\
			\
			__asm	movaps	xmm0,[edi     ]	/* t1 */					__asm	movaps	xmm2,[esi     ]	/* t5 */\
			__asm	movaps	xmm1,[edi+0x10]	/* t2 */					__asm	movaps	xmm3,[esi+0x10]	/* t6 */\
			__asm	movaps	xmm4,[edi     ]	/* xmm2 <- cpy t1 */		__asm	movaps	xmm6,[esi     ]	/* xmm2 <- cpy t5 */\
			__asm	movaps	xmm5,[edi+0x10]	/* xmm3 <- cpy t2 */		__asm	movaps	xmm7,[esi+0x10]	/* xmm3 <- cpy t6 */\
			\
			__asm	sub	esi, edi	/* stride */\
			__asm	shl	esi, 1		/* stride*2 */\
			__asm	add	edi, esi	/* add_in2 */\
			__asm	shr	esi, 1		/* stride */\
			__asm	add	esi, edi	/* add_in3 */\
			\
			__asm	addpd	xmm0,[edi     ]	/* ~t1 = t1 + t5 */			__asm	addpd	xmm2,[esi     ]	/* ~t3 = t3 + t7 */\
			__asm	addpd	xmm1,[edi+0x10]	/* ~t2 = t2 + t6 */			__asm	addpd	xmm3,[esi+0x10]	/* ~t4 = t4 + t8 */\
			__asm	subpd	xmm4,[edi     ]	/* ~t5 = t1 - t5 */			__asm	subpd	xmm6,[esi     ]	/* ~t7 = t3 - t7 */\
			__asm	subpd	xmm5,[edi+0x10]	/* ~t6 = t2 - t6 */			__asm	subpd	xmm7,[esi+0x10]	/* ~t8 = t4 - t8 */\
			\
			/* Finish radix-4 butterfly and store results into main-array slots: */\
			\
			__asm	subpd	xmm0,xmm2	/* ~t3 <- t1 -t3 */				__asm	subpd	xmm4,xmm7	/* ~t5 <- t5 -t8 */\
			__asm	subpd	xmm1,xmm3	/* ~t4 <- t2 -t4 */				__asm	subpd	xmm5,xmm6	/* ~t8 <- t6 -t7 */\
			__asm	movaps	[ebx      ],xmm0	/* <- ~t3 */			__asm	movaps	[ecx      ],xmm4	/* <- ~t5 */\
			__asm	movaps	[ebx+0x010],xmm1	/* <- ~t4 */			__asm	movaps	[edx+0x010],xmm5	/* <- ~t8 */\
			__asm	addpd	xmm2,xmm2	/*          2*t3 */				__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm3,xmm3	/*          2*t4 */				__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm2,xmm0	/* ~t1 <- t1 +t3 */				__asm	addpd	xmm7,xmm4	/* ~t7 <- t5 +t8 */\
			__asm	addpd	xmm3,xmm1	/* ~t2 <- t2 +t4 */				__asm	addpd	xmm6,xmm5	/* ~t6 <- t6 +t7 */\
			__asm	movaps	[eax      ],xmm2	/* <- ~t1 */			__asm	movaps	[edx      ],xmm7	/* <- ~t7 */\
			__asm	movaps	[eax+0x010],xmm3	/* <- ~t2 */			__asm	movaps	[ecx+0x010],xmm6	/* <- ~t6 */\
		}

		/* Variant of SSE2_RADIX4_DIF_0TWIDDLE_STRIDE, which assumes __add0,1,2,3 enter in out0,out1,out2,out3,
		and that the contents of these 4 register aliases must not be changed on return.
		*/
		#define SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(__in, __stride1,__stride2, out0,out1,out2,out3)\
		{\
			__asm	mov	edi, __in	\
			__asm	mov	esi, __stride1\
			__asm	add	esi, edi	/* add_in1 */\
			\
			__asm	movaps	xmm4,[edi     ]	/* t1 */					__asm	movaps	xmm6,[esi     ]	/* t3 */\
			__asm	movaps	xmm5,[edi+0x10]	/* t2 */					__asm	movaps	xmm7,[esi+0x10]	/* t4 */\
			__asm	add	edi, __stride2\
			__asm	add	esi, __stride2\
			__asm	movaps	xmm0,[edi     ]	/* t5 */					__asm	movaps	xmm2,[esi     ]	/* t7 */\
			__asm	movaps	xmm1,[edi+0x10]	/* t6 */					__asm	movaps	xmm3,[esi+0x10]	/* t8 */\
			\
			__asm	subpd	xmm4,xmm0		/*~t5 = t1-t5 */			__asm	subpd	xmm6,xmm2		/* ~t7 = t3-t7 */\
			__asm	subpd	xmm5,xmm1		/*~t6 = t2-t6 */			__asm	subpd	xmm7,xmm3		/* ~t8 = t4-t8 */\
			__asm	addpd	xmm0,xmm0		/*       2*t5 */			__asm	addpd	xmm2,xmm2		/*        2*t7 */\
			__asm	addpd	xmm1,xmm1		/*       2*t6 */			__asm	addpd	xmm3,xmm3		/*        2*t8 */\
			__asm	addpd	xmm0,xmm4		/*~t1 = t1+t5 */			__asm	addpd	xmm2,xmm6		/* ~t3 = t3+t7 */\
			__asm	addpd	xmm1,xmm5		/*~t2 = t2+t6 */			__asm	addpd	xmm3,xmm7		/* ~t4 = t4+t8 */\
			/* Finish radix-4 butterfly and store results into main-array slots: */\
			\
			__asm	subpd	xmm0,xmm2	/* ~t3 <- t1 -t3 */				__asm	subpd	xmm4,xmm7	/* ~t5 <- t5 -t8 */\
			__asm	subpd	xmm1,xmm3	/* ~t4 <- t2 -t4 */				__asm	subpd	xmm5,xmm6	/* ~t8 <- t6 -t7 */\
			__asm	movaps	[out1      ],xmm0	/* <- ~t3 */			__asm	movaps	[out2      ],xmm4	/* <- ~t5 */\
			__asm	movaps	[out1+0x010],xmm1	/* <- ~t4 */			__asm	movaps	[out3+0x010],xmm5	/* <- ~t8 */\
			__asm	addpd	xmm2,xmm2	/*          2*t3 */				__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm3,xmm3	/*          2*t4 */				__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm2,xmm0	/* ~t1 <- t1 +t3 */				__asm	addpd	xmm7,xmm4	/* ~t7 <- t5 +t8 */\
			__asm	addpd	xmm3,xmm1	/* ~t2 <- t2 +t4 */				__asm	addpd	xmm6,xmm5	/* ~t6 <- t6 +t7 */\
			__asm	movaps	[out0      ],xmm2	/* <- ~t1 */			__asm	movaps	[out3      ],xmm7	/* <- ~t7 */\
			__asm	movaps	[out0+0x010],xmm3	/* <- ~t2 */			__asm	movaps	[out2+0x010],xmm6	/* <- ~t6 */\
		}

		/* In-place Variant of SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C. We use separate pseudonyms for the in-and-out-address
		registers to allow for any desired index permutation here: */
		#define SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_D(in0,in1,in2,in3, out0,out1,out2,out3)\
		{\
			__asm	movaps	xmm4,[in0     ]	/* t1 */					__asm	movaps	xmm6,[in1     ]	/* t3 */\
			__asm	movaps	xmm5,[in0+0x10]	/* t2 */					__asm	movaps	xmm7,[in1+0x10]	/* t4 */\
			__asm	movaps	xmm0,[in2     ]	/* t5 */					__asm	movaps	xmm2,[in3     ]	/* t7 */\
			__asm	movaps	xmm1,[in2+0x10]	/* t6 */					__asm	movaps	xmm3,[in3+0x10]	/* t8 */\
			\
			__asm	subpd	xmm4,xmm0		/*~t5 = t1-t5 */			__asm	subpd	xmm6,xmm2		/* ~t7 = t3-t7 */\
			__asm	subpd	xmm5,xmm1		/*~t6 = t2-t6 */			__asm	subpd	xmm7,xmm3		/* ~t8 = t4-t8 */\
			__asm	addpd	xmm0,xmm0		/*       2*t5 */			__asm	addpd	xmm2,xmm2		/*        2*t7 */\
			__asm	addpd	xmm1,xmm1		/*       2*t6 */			__asm	addpd	xmm3,xmm3		/*        2*t8 */\
			__asm	addpd	xmm0,xmm4		/*~t1 = t1+t5 */			__asm	addpd	xmm2,xmm6		/* ~t3 = t3+t7 */\
			__asm	addpd	xmm1,xmm5		/*~t2 = t2+t6 */			__asm	addpd	xmm3,xmm7		/* ~t4 = t4+t8 */\
			/* Finish radix-4 butterfly and store results into main-array slots: */\
			\
			__asm	subpd	xmm0,xmm2	/* ~t3 <- t1 -t3 */				__asm	subpd	xmm4,xmm7	/* ~t5 <- t5 -t8 */\
			__asm	subpd	xmm1,xmm3	/* ~t4 <- t2 -t4 */				__asm	subpd	xmm5,xmm6	/* ~t8 <- t6 -t7 */\
			__asm	movaps	[out1      ],xmm0	/* <- ~t3 */			__asm	movaps	[out2      ],xmm4	/* <- ~t5 */\
			__asm	movaps	[out1+0x010],xmm1	/* <- ~t4 */			__asm	movaps	[out3+0x010],xmm5	/* <- ~t8 */\
			__asm	addpd	xmm2,xmm2	/*          2*t3 */				__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm3,xmm3	/*          2*t4 */				__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm2,xmm0	/* ~t1 <- t1 +t3 */				__asm	addpd	xmm7,xmm4	/* ~t7 <- t5 +t8 */\
			__asm	addpd	xmm3,xmm1	/* ~t2 <- t2 +t4 */				__asm	addpd	xmm6,xmm5	/* ~t6 <- t6 +t7 */\
			__asm	movaps	[out0      ],xmm2	/* <- ~t1 */			__asm	movaps	[out3      ],xmm7	/* <- ~t7 */\
			__asm	movaps	[out0+0x010],xmm3	/* <- ~t2 */			__asm	movaps	[out2+0x010],xmm6	/* <- ~t6 */\
		}

		/* DIF radix-4 subconvolution, sans twiddles, inputs in __i0-3, outputs in __o0-3, possibly coincident with inputs: */
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

		/* The SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO variant means that we munge the second set of 4 output as follows:

			t3 =(t3+t4)*ISRT2;t4 =(t4-t3)*ISRT2;
			t7 =(t7+t8)*ISRT2;t8 =(t8-t7)*ISRT2;
		*/
		/* Full-inline-asm version of SSE2_RADIX4_DIT_0TWIDDLE. Assumes base addresses add0-3 in e[a-d]x,
		base offset in edi, these must remain unchanged, so esi is available for temporary storage. */
		#define SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(__add0, __add1, __add2, __add3, __tmp)\
		{\
			__asm	mov	eax, __add0\
			__asm	mov	ebx, __add1\
			__asm	mov	ecx, __add2\
			__asm	mov	edx, __add3\
			\
			__asm	movaps	xmm2,[eax     ]	/* t1 */				__asm	movaps	xmm6,[ecx     ]	/* t5 */\
			__asm	movaps	xmm3,[eax+0x10]	/* t2 */				__asm	movaps	xmm7,[ecx+0x10]	/* t6 */\
			__asm	movaps	xmm0,[ebx     ]	/* t3 */				__asm	movaps	xmm4,[edx     ]	/* t7 */\
			__asm	movaps	xmm1,[ebx+0x10]	/* t4 */				__asm	movaps	xmm5,[edx+0x10]	/* t8 */\
			\
			__asm	subpd	xmm2,xmm0		/* ~t3 = t1-t3 */		__asm	subpd	xmm6,xmm4		/* ~t7 = t5-t7 */\
			__asm	subpd	xmm3,xmm1		/* ~t4 = t2-t4 */		__asm	subpd	xmm7,xmm5		/* ~t8 = t6-t8 */\
			__asm	addpd	xmm0,xmm0		/*        2*t3 */		__asm	addpd	xmm4,xmm4		/*        2*t7 */\
			__asm	addpd	xmm1,xmm1		/*        2*t4 */		__asm	addpd	xmm5,xmm5		/*        2*t8 */\
			__asm	addpd	xmm0,xmm2		/* ~t1 = t1+t3 */		__asm	addpd	xmm4,xmm6		/* ~t5 = t5+t7 */\
			__asm	addpd	xmm1,xmm3		/* ~t2 = t2+t4 */		__asm	addpd	xmm5,xmm7		/* ~t6 = t6+t8 */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			__asm	mov	eax, __tmp										__asm	mov	edx, isrt2\
			\
			__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */		\
			__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */		\
			__asm	movaps	[eax+0x040],xmm0	/* <- ~t5 */	\
			__asm	movaps	[eax+0x050],xmm1	/* <- ~t6 */			__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */\
			__asm	addpd	xmm4,xmm4	/*          2*t5 */				__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */\
			__asm	addpd	xmm5,xmm5	/*          2*t6 */				__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */				__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */				__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */\
			__asm	movaps	[eax      ],xmm4	/* <- ~t1 */			__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */\
			__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */	\
			\
			/*\
			t3 =(t3+t4)*ISRT2;t4 =(t4-t3)*ISRT2;\
			t7 =(t7+t8)*ISRT2;t8 =(t8-t7)*ISRT2;\
			*/\
			__asm	movaps	xmm5,[edx]	/* ISRT2 */\
			__asm	movaps	xmm0,xmm3	/* cpy t4 */\
			__asm	movaps	xmm1,xmm6	/* cpy t8 */\
			__asm	subpd	xmm3,xmm7	/* 4-3*/\
			__asm	subpd	xmm6,xmm2	/* 8-7*/\
			__asm	addpd	xmm0,xmm7	/* 4+3*/\
			__asm	addpd	xmm1,xmm2	/* 8+7*/\
			__asm	mulpd	xmm3,xmm5	/* (4-3)*ISRT2 */\
			__asm	mulpd	xmm6,xmm5	/* (8-7)*ISRT2 */\
			__asm	mulpd	xmm0,xmm5	/* (4+3)*ISRT2 */\
			__asm	mulpd	xmm1,xmm5	/* (8+7)*ISRT2 */\
			__asm	movaps	[eax+0x030],xmm3	/* a[jp+p4 ] <- ~t4 */\
			__asm	movaps	[eax+0x070],xmm6	/* a[jp+p12] <- ~t8 */\
			__asm	movaps	[eax+0x020],xmm0	/* a[jt+p4 ] <- ~t3 */\
			__asm	movaps	[eax+0x060],xmm1	/* a[jt+p12] <- ~t7 */\
		}

		/* Full-inline-asm version of SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO. Assumes base addresses add0-3 in e[a-d]x,
		base offset in edi, these must remain unchanged, so esi is available for temporary storage. */
		#define SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(__tmp)\
		{\
			__asm	movaps	xmm2,[eax     ]	/* t1 */				__asm	movaps	xmm6,[ecx     ]	/* t5 */\
			__asm	movaps	xmm3,[eax+0x10]	/* t2 */				__asm	movaps	xmm7,[ecx+0x10]	/* t6 */\
			__asm	movaps	xmm0,[ebx     ]	/* t3 */				__asm	movaps	xmm4,[edx     ]	/* t7 */\
			__asm	movaps	xmm1,[ebx+0x10]	/* t4 */				__asm	movaps	xmm5,[edx+0x10]	/* t8 */\
			\
			__asm	subpd	xmm2,xmm0		/* ~t3 = t1-t3 */		__asm	subpd	xmm6,xmm4		/* ~t7 = t5-t7 */\
			__asm	subpd	xmm3,xmm1		/* ~t4 = t2-t4 */		__asm	subpd	xmm7,xmm5		/* ~t8 = t6-t8 */\
			__asm	addpd	xmm0,xmm0		/*        2*t3 */		__asm	addpd	xmm4,xmm4		/*        2*t7 */\
			__asm	addpd	xmm1,xmm1		/*        2*t4 */		__asm	addpd	xmm5,xmm5		/*        2*t8 */\
			__asm	addpd	xmm0,xmm2		/* ~t1 = t1+t3 */		__asm	addpd	xmm4,xmm6		/* ~t5 = t5+t7 */\
			__asm	addpd	xmm1,xmm3		/* ~t2 = t2+t4 */		__asm	addpd	xmm5,xmm7		/* ~t6 = t6+t8 */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			__asm	mov	esi, __tmp\
			\
			__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */		\
			__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */		\
			__asm	movaps	[esi+0x040],xmm0	/* <- ~t5 */	\
			__asm	movaps	[esi+0x050],xmm1	/* <- ~t6 */			__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */\
			__asm	addpd	xmm4,xmm4	/*          2*t5 */				__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */\
			__asm	addpd	xmm5,xmm5	/*          2*t6 */				__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */				__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */				__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */\
			__asm	movaps	[esi      ],xmm4	/* <- ~t1 */			__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */\
			__asm	movaps	[esi+0x010],xmm5	/* <- ~t2 */	\
			\
			/*\
			t3 =(t3+t4)*ISRT2;t4 =(t4-t3)*ISRT2;\
			t7 =(t7+t8)*ISRT2;t8 =(t8-t7)*ISRT2;\
			*/\
			__asm	mov	esi, isrt2\
			__asm	movaps	xmm5,[esi]	/* ISRT2 */\
			__asm	movaps	xmm0,xmm3	/* cpy t4 */\
			__asm	movaps	xmm1,xmm6	/* cpy t8 */\
			__asm	subpd	xmm3,xmm7	/* 4-3*/\
			__asm	subpd	xmm6,xmm2	/* 8-7*/\
			__asm	addpd	xmm0,xmm7	/* 4+3*/\
			__asm	addpd	xmm1,xmm2	/* 8+7*/\
			__asm	mov	esi, __tmp\
			__asm	mulpd	xmm3,xmm5	/* (4-3)*ISRT2 */\
			__asm	mulpd	xmm6,xmm5	/* (8-7)*ISRT2 */\
			__asm	mulpd	xmm0,xmm5	/* (4+3)*ISRT2 */\
			__asm	mulpd	xmm1,xmm5	/* (8+7)*ISRT2 */\
			__asm	movaps	[esi+0x030],xmm3	/* a[jp+p4 ] <- ~t4 */\
			__asm	movaps	[esi+0x070],xmm6	/* a[jp+p12] <- ~t8 */\
			__asm	movaps	[esi+0x020],xmm0	/* a[jt+p4 ] <- ~t3 */\
			__asm	movaps	[esi+0x060],xmm1	/* a[jt+p12] <- ~t7 */\
		}

		/* DIT radix-4 subconvolution, sans twiddles - this is the in-place version needed by the wrapper_square routines.
		Assumes the 4 addresses in question, __add0, __add1, __add2, __add3, enter in eax,ebx,ecx,edx, respectively:
		*/
		#define SSE2_RADIX4_DIT_IN_PLACE()\
		{\
			__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */\
			__asm	movaps	xmm2,[eax     ]	/* xmm2 <- cpy a[jt   ] */\
			__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[jp   ] */\
			\
			__asm	addpd	xmm0,[ebx     ]	/* t1 */\
			__asm	addpd	xmm1,[ebx+0x10]	/* t2 */\
			__asm	subpd	xmm2,[ebx     ]	/* t3 */\
			__asm	subpd	xmm3,[ebx+0x10]	/* t4 */\
			\
			__asm	movaps	xmm4,[ecx     ]	/* a[jt+p2] */\
			__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p2] */\
			__asm	movaps	xmm6,[ecx     ]	/* xmm2 <- cpy a[jt+p2] */\
			__asm	movaps	xmm7,[ecx+0x10]	/* xmm3 <- cpy a[jp+p2] */\
			\
			__asm	addpd	xmm4,[edx     ]	/* t5 */\
			__asm	addpd	xmm5,[edx+0x10]	/* t6 */\
			__asm	subpd	xmm6,[edx     ]	/* t7 */\
			__asm	subpd	xmm7,[edx+0x10]	/* t8 */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */\
			__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */\
			__asm	movaps	[ebx      ],xmm0	/* <- ~t5 */\
			__asm	movaps	[ebx+0x010],xmm1	/* <- ~t6 */\
			__asm	addpd	xmm4,xmm4	/*          2*t5 */\
			__asm	addpd	xmm5,xmm5	/*          2*t6 */\
			__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */\
			__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */\
			__asm	movaps	[eax      ],xmm4	/* <- ~t1 */\
			__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */\
			\
			__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */\
			__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */\
			__asm	movaps	[edx      ],xmm2	/* <- ~t7 */\
			__asm	movaps	[ecx+0x010],xmm3	/* <- ~t4 */\
			__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */\
			__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */\
			__asm	movaps	[ecx      ],xmm7	/* <- ~t3 */\
			__asm	movaps	[edx+0x010],xmm6	/* <- ~t8 */\
		}

		/* DIT radix-4 subconvolution, sans twiddles, with low address and stride as arguments. Stride must be a numeric [decimal or hex] literal, e.g. 16, 0x100 or 100H: */
		#define SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(__reg0, __reg1, __reg2, __reg3, __reg4, __reg5, __reg6, __reg7)\
		{\
		/*\
			a[jt    ]=t00+t20;			a[jp    ]=t01+t21;\
			a[jt+p10]=t00-t20;			a[jp+p10]=t01-t21;\
		*/\
			__asm	addpd	__reg0,__reg4		/* t1 +t17 */\
			__asm	addpd	__reg1,__reg5		/* t2 +t18 */\
			__asm	movaps	[eax     ],__reg0	/* a[jt+p0 ] */\
			__asm	movaps	[eax+0x10],__reg1	/* a[jp+p0 ] */\
			__asm	addpd	__reg4,__reg4		/*   2*t17 */\
			__asm	addpd	__reg5,__reg5		/*   2*t18 */\
			__asm	subpd	__reg0,__reg4		/*~t17 <- t1 -t17 */\
			__asm	subpd	__reg1,__reg5		/*~t18 <- t2 -t18 */\
			__asm	movaps	[ecx     ],__reg0	/* a[jt+p8 ] */\
			__asm	movaps	[ecx+0x10],__reg1	/* a[jp+p8 ] */\
		/*\
			mpy by E^-4 = -I is inlined here...\
			a[jt+p08]=t10+t31;			a[jp+p08]=t11-t30;\
			a[jt+p18]=t10-t31;			a[jp+p18]=t11+t30;\
		*/\
			__asm	addpd	__reg2,__reg7		/* rt  <- t9 +t26 */\
			__asm	subpd	__reg3,__reg6		/* it  <- t10-t25 */\
			__asm	movaps	[ebx     ],__reg2	/* a[jt+p4 ] */\
			__asm	movaps	[ebx+0x10],__reg3	/* a[jp+p4 ] */\
			__asm	addpd	__reg7,__reg7		/*          2*t26 */\
			__asm	addpd	__reg6,__reg6		/*          2*t25 */\
			__asm	subpd	__reg2,__reg7		/*~t26 <- t9 -t26 */\
			__asm	addpd	__reg3,__reg6		/*~t25 <- t10+t25 */\
			__asm	movaps	[edx     ],__reg2	/* a[jt+p12] */\
			__asm	movaps	[edx+0x10],__reg3	/* a[jp+p12] */\
		}

		/* DIF radix-4 subconvolution, sans twiddles - this is the in-place version needed by the carry routines */
		#define SSE2_RADIX4_DIF_IN_PLACE(__add0, __add1, __add2, __add3)\
		{\
			__asm	mov	eax, __add0\
			__asm	mov	ebx, __add1\
			__asm	mov	ecx, __add2\
			__asm	mov	edx, __add3\
			\
			__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */\
			__asm	movaps	xmm2,[eax     ]	/* xmm2 <- cpy a[jt   ] */\
			__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[jp   ] */\
			\
			__asm	addpd	xmm0,[ebx     ]	/* t1 */\
			__asm	addpd	xmm1,[ebx+0x10]	/* t2 */\
			__asm	subpd	xmm2,[ebx     ]	/* t3 */\
			__asm	subpd	xmm3,[ebx+0x10]	/* t4 */\
			\
			__asm	movaps	xmm4,[ecx     ]	/* a[jt+p2] */\
			__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p2] */\
			__asm	movaps	xmm6,[ecx     ]	/* xmm2 <- cpy a[jt+p2] */\
			__asm	movaps	xmm7,[ecx+0x10]	/* xmm3 <- cpy a[jp+p2] */\
			\
			__asm	addpd	xmm4,[edx     ]	/* t5 */\
			__asm	addpd	xmm5,[edx+0x10]	/* t6 */\
			__asm	subpd	xmm6,[edx     ]	/* t7 */\
			__asm	subpd	xmm7,[edx+0x10]	/* t8 */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */\
			__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */\
			__asm	movaps	[ebx      ],xmm0	/* <- ~t5 */\
			__asm	movaps	[ebx+0x010],xmm1	/* <- ~t6 */\
			__asm	addpd	xmm4,xmm4	/*          2*t5 */\
			__asm	addpd	xmm5,xmm5	/*          2*t6 */\
			__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */\
			__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */\
			__asm	movaps	[eax      ],xmm4	/* <- ~t1 */\
			__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */\
			\
			__asm	subpd	xmm2,xmm7	/* ~t3 <- t3 -t8 */\
			__asm	subpd	xmm3,xmm6	/* ~t8 <- t4 -t7 */\
			__asm	movaps	[ecx      ],xmm2	/* <- ~t3 */\
			__asm	movaps	[edx+0x010],xmm3	/* <- ~t8 */\
			__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm7,xmm2	/* ~t7 <- t3 +t8 */\
			__asm	addpd	xmm6,xmm3	/* ~t4 <- t4 +t7 */\
			__asm	movaps	[edx      ],xmm7	/* <- ~t7 */\
			__asm	movaps	[ecx+0x010],xmm6	/* <- ~t4 */\
		}

		/* The SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO variant means that we munge the second set of 4 output as follows:

			t3 =(t3-t4)*ISRT2;	t4 =(t3+t4)*ISRT2;
			t7 =(t7+t8)*ISRT2;	t8 =(t7-t8)*ISRT2;
		*/
		#define SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(__add0, __add1, __add2, __add3)\
		{\
			__asm	mov	eax, __add0\
			__asm	mov	ebx, __add1\
			__asm	mov	ecx, __add2\
			__asm	mov	edx, __add3\
			\
			__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */\
			__asm	movaps	xmm2,[eax     ]	/* xmm2 <- cpy a[jt   ] */\
			__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[jp   ] */\
			\
			__asm	addpd	xmm0,[ebx     ]	/* t1 */\
			__asm	addpd	xmm1,[ebx+0x10]	/* t2 */\
			__asm	subpd	xmm2,[ebx     ]	/* t3 */\
			__asm	subpd	xmm3,[ebx+0x10]	/* t4 */\
			\
			__asm	movaps	xmm4,[ecx     ]	/* a[jt+p2] */\
			__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p2] */\
			__asm	movaps	xmm6,[ecx     ]	/* xmm2 <- cpy a[jt+p2] */\
			__asm	movaps	xmm7,[ecx+0x10]	/* xmm3 <- cpy a[jp+p2] */\
			\
			__asm	addpd	xmm6,[edx     ]	/* t5 */\
			__asm	addpd	xmm7,[edx+0x10]	/* t6 */\
			__asm	subpd	xmm4,[edx     ]	/* t7 */\
			__asm	subpd	xmm5,[edx+0x10]	/* t8 */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			/*\
			~t5 =t1 -t5;		~t1 =t1 +t5;										~t7 =t3 +t8;		~t3 =t3 -t8;\
			~t6 =t2 -t6;		~t2 =t2 +t6;										~t8 =t4 -t7;		~t4 =t4 +t7;\
			*/\
			__asm	subpd	xmm0,xmm6	/*~t5 */									__asm	mov	edi, isrt2\
			__asm	subpd	xmm1,xmm7	/*~t6 */									\
			__asm	movaps	[ebx      ],xmm0	/* a[jt+p8 ] <- ~t5 */				__asm	subpd	xmm2,xmm5	/*~t3 */\
			__asm	movaps	[ebx+0x010],xmm1	/* a[jp+p8 ] <- ~t6 */				__asm	subpd	xmm3,xmm4	/*~t8 */\
			__asm	addpd	xmm6,xmm6	/* 2*t5 */									__asm	addpd	xmm5,xmm5	/* 2*t8 */\
			__asm	addpd	xmm7,xmm7	/* 2*t6 */									__asm	addpd	xmm4,xmm4	/* 2*t7 */\
			__asm	addpd	xmm6,xmm0	/*~t1 */									__asm	addpd	xmm5,xmm2	/*~t7 */\
			__asm	addpd	xmm7,xmm1	/*~t2 */									__asm	addpd	xmm4,xmm3	/*~t4 */\
			__asm	movaps	[eax      ],xmm6	/* a[jt    ] <- ~t1 */\
			__asm	movaps	[eax+0x010],xmm7	/* a[jp    ] <- ~t2 */\
			/*\
			t3 =(t3-t4)*ISRT2;	t4 =(t3+t4)*ISRT2;\
			t7 =(t7-t8)*ISRT2;	t8 =(t7+t8)*ISRT2;\
			*/\
			__asm	movaps	xmm6,xmm2	/* cpy t3 */\
			__asm	movaps	xmm7,xmm5	/* cpy t7 */\
			__asm	subpd	xmm2,xmm4	/* 3-4*/\
			__asm	subpd	xmm5,xmm3	/* 7-8*/\
			__asm	addpd	xmm6,xmm4	/* 3+4*/\
			__asm	addpd	xmm7,xmm3	/* 7+8*/\
			__asm	mulpd	xmm2,[edi]	/* (3-4)*ISRT2 */\
			__asm	mulpd	xmm5,[edi]	/* (7-8)*ISRT2 */\
			__asm	mulpd	xmm6,[edi]	/* (3+4)*ISRT2 */\
			__asm	mulpd	xmm7,[edi]	/* (7+8)*ISRT2 */\
			__asm	movaps	[ecx      ],xmm2	/* a[jt+p4 ] <- ~t3 */\
			__asm	movaps	[edx      ],xmm5	/* a[jp+p12] <- ~t7 */\
			__asm	movaps	[ecx+0x010],xmm6	/* a[jp+p4 ] <- ~t4 */\
			__asm	movaps	[edx+0x010],xmm7	/* a[jt+p12] <- ~t8 */\
		}

		/* In-place version of SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO.
		Assumes the 4 addresses in question, __add0, __add1, __add2, __add3, enter in eax,ebx,ecx,edx, respectively:
		*/
		#define SSE2_RADIX4_DIT_IN_PLACE_2NDOFTWO()\
		{\
			__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */\
			__asm	movaps	xmm2,[eax     ]	/* xmm2 <- cpy a[jt   ] */\
			__asm	movaps	xmm3,[eax+0x10]	/* xmm3 <- cpy a[jp   ] */\
			\
			__asm	addpd	xmm0,[ebx     ]	/* t1 */\
			__asm	addpd	xmm1,[ebx+0x10]	/* t2 */\
			__asm	subpd	xmm2,[ebx     ]	/* t3 */\
			__asm	subpd	xmm3,[ebx+0x10]	/* t4 */\
			\
			__asm	movaps	xmm4,[ecx     ]	/* a[jt+p2] */\
			__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p2] */\
			__asm	movaps	xmm6,[ecx     ]	/* xmm2 <- cpy a[jt+p2] */\
			__asm	movaps	xmm7,[ecx+0x10]	/* xmm3 <- cpy a[jp+p2] */\
			\
			__asm	addpd	xmm4,[edx     ]	/* t5 */\
			__asm	addpd	xmm5,[edx+0x10]	/* t6 */\
			__asm	subpd	xmm6,[edx     ]	/* t7 */\
			__asm	subpd	xmm7,[edx+0x10]	/* t8 */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			\
			__asm	subpd	xmm0,xmm4	/* ~t5 <- t1 -t5 */\
			__asm	subpd	xmm1,xmm5	/* ~t6 <- t2 -t6 */\
			__asm	movaps	[ebx      ],xmm0	/* <- ~t5 */\
			__asm	movaps	[ebx+0x010],xmm1	/* <- ~t6 */\
			__asm	addpd	xmm4,xmm4	/*          2*t5 */\
			__asm	addpd	xmm5,xmm5	/*          2*t6 */\
			__asm	addpd	xmm4,xmm0	/* ~t1 <- t1 +t5 */\
			__asm	addpd	xmm5,xmm1	/* ~t2 <- t2 +t6 */\
			__asm	movaps	[eax      ],xmm4	/* <- ~t1 */\
			__asm	movaps	[eax+0x010],xmm5	/* <- ~t2 */\
			\
			__asm	mov	ebx, isrt2\
			\
			__asm	subpd	xmm2,xmm7	/* ~t7 <- t3 -t8 */\
			__asm	subpd	xmm3,xmm6	/* ~t4 <- t4 -t7 */\
			__asm	addpd	xmm7,xmm7	/*          2*t8 */\
			__asm	addpd	xmm6,xmm6	/*          2*t7 */\
			__asm	addpd	xmm7,xmm2	/* ~t3 <- t3 +t8 */\
			__asm	addpd	xmm6,xmm3	/* ~t8 <- t4 +t7 */\
			\
			/*\
			t3 =(t3+t4)*ISRT2;t4 =(t4-t3)*ISRT2;\
			t7 =(t7+t8)*ISRT2;t8 =(t8-t7)*ISRT2;\
			*/\
			__asm	movaps	xmm0,xmm3	/* cpy t4 */\
			__asm	movaps	xmm1,xmm6	/* cpy t8 */\
			__asm	subpd	xmm3,xmm7	/* 4-3*/\
			__asm	subpd	xmm6,xmm2	/* 8-7*/\
			__asm	addpd	xmm0,xmm7	/* 4+3*/\
			__asm	addpd	xmm1,xmm2	/* 8+7*/\
			__asm	mulpd	xmm3,[ebx]	/* (4-3)*ISRT2 */\
			__asm	mulpd	xmm6,[ebx]	/* (8-7)*ISRT2 */\
			__asm	mulpd	xmm0,[ebx]	/* (4+3)*ISRT2 */\
			__asm	mulpd	xmm1,[ebx]	/* (8+7)*ISRT2 */\
			__asm	movaps	[ecx+0x010],xmm3	/* a[jp+p4 ] <- ~t4 */\
			__asm	movaps	[edx+0x010],xmm6	/* a[jp+p12] <- ~t8 */\
			__asm	movaps	[ecx      ],xmm0	/* a[jt+p4 ] <- ~t3 */\
			__asm	movaps	[edx      ],xmm1	/* a[jt+p12] <- ~t7 */\
		}

		/* Combine the 2 radix-4 subtransforms: r00,01,02,03,04,05,06,07 in registers xmm6,7,2,4,0,1,5,3
		in addition to the r-temps, if needed.

		The SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO variant means the ISRT2 stuff is already done
		on the t07,07,0e,0f outputs, i.e. both radix-2 combos here look identical.
		*/
		/* DIF radix-4 subconvolution, with twiddles.	Cost: 30 MOVapd, 28 ADD/SUBpd, 16 MULpd.
		Assumes the 4 twiddles c0,4,8,12 are bit-reversed in memory, i.e. c0,8,4,12 are adjacent in memory starting at c0.
		*/
		#define SSE2_RADIX4_DIF_4TWIDDLE(__add0, __add1, __add2, __add3, __tmp, __c0)\
		{\
			/* Do the p0,p8 combo: */\
			__asm	mov	eax, __add0\
			__asm	mov	ebx, __c0\
			__asm	mov	ecx, __add2\
			\
			__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */				__asm	movaps	xmm4,[ecx     ]	/* a[jt+p8 ] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p8 ] */\
			__asm	movaps	xmm6,[ebx     ]	/* c0 */\
			__asm	movaps	xmm7,[ebx+0x10]	/* s0 */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[jt   ] */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[jp   ] */\
			\
			__asm	mulpd	xmm0,xmm6		/* a[jt   ]*c0 */\
			__asm	mulpd	xmm1,xmm6		/* a[jp   ]*c0 */\
			__asm	mulpd	xmm2,xmm7		/* a[jt   ]*s0 */\
			__asm	mulpd	xmm3,xmm7		/* a[jp   ]*s0	xmm6,7 free */\
																		__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[jt+p8 ] */\
			__asm	addpd	xmm1,xmm2		/* xmm1 <- t2 */			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[jp+p8 ] */\
																		__asm	mulpd	xmm4,[ebx+0x20] /* a[jt+p8 ]*c8 */\
			__asm	subpd	xmm0,xmm3		/* xmm0 <- t1 */			__asm	mulpd	xmm5,[ebx+0x20]	/* a[jp+p8 ]*c8 */\
																		__asm	mulpd	xmm6,[ebx+0x30]	/* a[jt+p8 ]*s8 */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy t1 */		__asm	mulpd	xmm7,[ebx+0x30]	/* a[jp+p8 ]*s8 */\
																		__asm	addpd	xmm5,xmm6	    /* xmm5 <- it */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy t2 */		__asm	subpd	xmm4,xmm7		/* xmm4 <- rt    xmm6,7 free */\
			\
																		__asm	mov	ecx, __add3\
																		__asm	add	ebx, 0x60\
																		__asm	movaps	xmm6,[ecx     ]	/* a[jt+p12] */\
																		__asm	movaps	xmm7,[ecx+0x10]	/* a[jp+p12] */\
			\
			__asm	addpd	xmm0,xmm4		/* ~t1 <- t1 +rt */\
			__asm	addpd	xmm1,xmm5		/* ~t2 <- t2 +it */\
			__asm	subpd	xmm2,xmm4		/* ~t3 <- t1 -rt */\
			__asm	subpd	xmm3,xmm5		/* ~t4 <- t2 -it	xmm4,5 free */\
			\
			/* Do the p4,12 combo: */\
			__asm	movaps	xmm4,xmm6		/* xmm4 <- cpy a[jt+p12] */\
			__asm	movaps	xmm5,xmm7		/* xmm5 <- cpy a[jp+p12] */\
			\
			__asm	mulpd	xmm4,[ebx     ]	/* a[jt+p12]*c12 */\
			__asm	mulpd	xmm5,[ebx     ]	/* a[jp+p12]*c12 */\
			__asm	mulpd	xmm6,[ebx+0x10]	/* a[jt+p12]*s12 */\
			__asm	mulpd	xmm7,[ebx+0x10]	/* a[jp+p12]*s12 */\
			__asm	mov	edx, __tmp	\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- it */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- rt */\
			__asm	movaps	[edx+0x010],xmm5	/* store it */\
			__asm	movaps	[edx      ],xmm4	/* store rt */\
			\
			__asm	mov	eax, __add1\
			__asm	sub	ebx, 0x20\
			__asm	movaps	xmm4,[eax     ]	/* a[jt+p4] */\
			__asm	movaps	xmm5,[eax+0x10]	/* a[jp+p4] */\
			__asm	movaps	xmm6,xmm4		/* xmm4 <- cpy a[jt+p4] */\
			__asm	movaps	xmm7,xmm5		/* xmm5 <- cpy a[jp+p4] */\
			\
			__asm	mulpd	xmm4,[ebx     ]	/* a[jt+p4]*c4 */\
			__asm	mulpd	xmm5,[ebx     ]	/* a[jp+p4]*c4 */\
			__asm	mulpd	xmm6,[ebx+0x10]	/* a[jt+p4]*s4 */\
			__asm	mulpd	xmm7,[ebx+0x10]	/* a[jp+p4]*s4 */\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- t6 */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- t5 	xmm6,7 free */\
			__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy t6 */\
			__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy t5 */\
			\
			__asm	subpd	xmm4,[edx      ]	/* ~t7 <- t5 -rt */\
			__asm	subpd	xmm5,[edx+0x010]	/* ~t8 <- t6 -it */\
			__asm	addpd	xmm6,[edx      ]	/* ~t5 <- t5 +rt */\
			__asm	addpd	xmm7,[edx+0x010]	/* ~t6 <- t6 +it */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			/*\
			~t5 =t1 -t5;		~t1 =t1 +t5;										~t7 =t3 +t8;		~t3 =t3 -t8;\
			~t6 =t2 -t6;		~t2 =t2 +t6;										~t8 =t4 -t7;		~t4 =t4 +t7;\
			*/\
			__asm	subpd	xmm0,xmm6	/*~t5 */									__asm	subpd	xmm2,xmm5	/*~t3 */\
			__asm	subpd	xmm1,xmm7	/*~t6 */									__asm	subpd	xmm3,xmm4	/*~t8 */\
			__asm	movaps	[edx+0x040],xmm0	/* a[jt+p8 ] <- ~t5 */				__asm	movaps	[edx+0x020],xmm2	/* a[jt+p4 ] <- ~t3 */\
			__asm	movaps	[edx+0x050],xmm1	/* a[jp+p8 ] <- ~t6 */				__asm	movaps	[edx+0x070],xmm3	/* a[jp+p12] <- ~t8 */\
			__asm	addpd	xmm6,xmm6	/* 2*t5 */									__asm	addpd	xmm5,xmm5	/* 2*t8 */\
			__asm	addpd	xmm7,xmm7	/* 2*t6 */									__asm	addpd	xmm4,xmm4	/* 2*t7 */\
			__asm	addpd	xmm6,xmm0	/*~t1 */									__asm	addpd	xmm5,xmm2	/*~t7 */\
			__asm	addpd	xmm7,xmm1	/*~t2 */									__asm	addpd	xmm4,xmm3	/*~t4 */\
			__asm	movaps	[edx      ],xmm6	/* a[jt    ] <- ~t1 */				__asm	movaps	[edx+0x060],xmm5	/* a[jt+p12] <- ~t7 */\
			__asm	movaps	[edx+0x010],xmm7	/* a[jp    ] <- ~t2 */				__asm	movaps	[edx+0x030],xmm4	/* a[jp+p4 ] <- ~t4 */\
			\
		}

		/* 2 Variants of SSE2_RADIX4_DIT_0TWIDDLE, which assume that the base offset [not used inside the macro] p4 enters in edi,
		addresses __add0,1,2,3 enter in eax,ebx,ecx,edx, and that the contents of these 5 regs must not be changed on return.
		Also assumes the 4 twiddles c0,4,8,12 are bit-reversed in memory, i.e. c0,8,4,12 are adjacent in memory starting at c0.
		[The A-version assumes c0 is the trivial 0th root of unity, so instead passes sincos data starting with c8].
		*/
		#define SSE2_RADIX4_DIF_4TWIDDLE_A(__tmp,__c8)\
		{\
			/* Do the p0,p8 combo: */\
			__asm	mov	esi, __c8\
			\
			__asm	movaps	xmm0,[eax     ]	/* t1 = a[jt   ] */			__asm	movaps	xmm4,[ecx     ]	/* a[jt+p8 ] */\
			__asm	movaps	xmm1,[eax+0x10]	/* t2 = a[jp   ] */			__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p8 ] */\
																		__asm	movaps	xmm2,[esi     ]	/* c8 */\
																		__asm	movaps	xmm3,[esi+0x10]	/* s8 */\
																		__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[jt+p8 ] */\
																		__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[jp+p8 ] */\
																		__asm	mulpd	xmm4,xmm2		/* a[jt+p8 ]*c8 */\
																		__asm	mulpd	xmm5,xmm2		/* a[jp+p8 ]*c8 */\
																		__asm	mulpd	xmm6,xmm3		/* a[jt+p8 ]*s8 */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy t1 */		__asm	mulpd	xmm7,xmm3		/* a[jp+p8 ]*s8 */\
																		__asm	addpd	xmm5,xmm6	    /* xmm5 <- it */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy t2 */		__asm	subpd	xmm4,xmm7		/* xmm4 <- rt    xmm6,7 free */\
			\
			__asm	addpd	xmm0,xmm4		/* ~t1 <- t1 +rt */\
			__asm	addpd	xmm1,xmm5		/* ~t2 <- t2 +it */\
			__asm	subpd	xmm2,xmm4		/* ~t3 <- t1 -rt */\
			__asm	subpd	xmm3,xmm5		/* ~t4 <- t2 -it	xmm4,5 free */\
			\
			/* Do the p4,12 combo: */\
			__asm	add	esi, 0x40\
			\
			__asm	movaps	xmm4,[edx     ]	/* a[jt+p12] */\
			__asm	movaps	xmm5,[edx+0x10]	/* a[jp+p12] */\
			__asm	movaps	xmm6,[edx     ]	/* xmm2 <- cpy a[jt+p12] */\
			__asm	movaps	xmm7,[edx+0x10]	/* xmm3 <- cpy a[jp+p12] */\
			\
			__asm	mulpd	xmm4,[esi     ]	/* a[jt+p12]*c12 */\
			__asm	mulpd	xmm5,[esi     ]	/* a[jp+p12]*c12 */\
			__asm	mulpd	xmm6,[esi+0x10]	/* a[jt+p12]*s12 */\
			__asm	mulpd	xmm7,[esi+0x10]	/* a[jp+p12]*s12 */\
			\
			__asm	mov	esi, __tmp	\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- it */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- rt */\
			__asm	movaps	[esi+0x010],xmm5	/* store it */\
			__asm	movaps	[esi      ],xmm4	/* store rt */\
			\
			__asm	mov	esi, __c8\
			__asm	add	esi, 0x20\
			__asm	movaps	xmm4,[ebx     ]	/* a[jt+p4] */\
			__asm	movaps	xmm5,[ebx+0x10]	/* a[jp+p4] */\
			__asm	movaps	xmm6,[ebx     ]	/* xmm2 <- cpy a[jt+p4] */\
			__asm	movaps	xmm7,[ebx+0x10]	/* xmm3 <- cpy a[jp+p4] */\
			\
			__asm	mulpd	xmm4,[esi     ]	/* a[jt+p4]*c4 */\
			__asm	mulpd	xmm5,[esi     ]	/* a[jp+p4]*c4 */\
			__asm	mulpd	xmm6,[esi+0x10]	/* a[jt+p4]*s4 */\
			__asm	mulpd	xmm7,[esi+0x10]	/* a[jp+p4]*s4 */\
			__asm	mov	esi, __tmp	\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- t6 */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- t5 	xmm6,7 free */\
			__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy t6 */\
			__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy t5 */\
			\
			__asm	subpd	xmm4,[esi      ]	/* ~t7 <- t5 -rt */\
			__asm	subpd	xmm5,[esi+0x010]	/* ~t8 <- t6 -it */\
			__asm	addpd	xmm6,[esi      ]	/* ~t5 <- t5 +rt */\
			__asm	addpd	xmm7,[esi+0x010]	/* ~t6 <- t6 +it */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			/*\
			~t5 =t1 -t5;		~t1 =t1 +t5;										~t7 =t3 +t8;		~t3 =t3 -t8;\
			~t6 =t2 -t6;		~t2 =t2 +t6;										~t8 =t4 -t7;		~t4 =t4 +t7;\
			*/\
			__asm	subpd	xmm0,xmm6	/*~t5 */									__asm	subpd	xmm2,xmm5	/*~t3 */\
			__asm	subpd	xmm1,xmm7	/*~t6 */									__asm	subpd	xmm3,xmm4	/*~t8 */\
			__asm	movaps	[esi+0x040],xmm0	/* a[jt+p8 ] <- ~t5 */				__asm	movaps	[esi+0x020],xmm2	/* a[jt+p4 ] <- ~t3 */\
			__asm	movaps	[esi+0x050],xmm1	/* a[jp+p8 ] <- ~t6 */				__asm	movaps	[esi+0x070],xmm3	/* a[jp+p12] <- ~t8 */\
			__asm	addpd	xmm6,xmm6	/* 2*t5 */									__asm	addpd	xmm5,xmm5	/* 2*t8 */\
			__asm	addpd	xmm7,xmm7	/* 2*t6 */									__asm	addpd	xmm4,xmm4	/* 2*t7 */\
			__asm	addpd	xmm6,xmm0	/*~t1 */									__asm	addpd	xmm5,xmm2	/*~t7 */\
			__asm	addpd	xmm7,xmm1	/*~t2 */									__asm	addpd	xmm4,xmm3	/*~t4 */\
			__asm	movaps	[esi      ],xmm6	/* a[jt    ] <- ~t1 */				__asm	movaps	[esi+0x060],xmm5	/* a[jt+p12] <- ~t7 */\
			__asm	movaps	[esi+0x010],xmm7	/* a[jp    ] <- ~t2 */				__asm	movaps	[esi+0x030],xmm4	/* a[jp+p4 ] <- ~t4 */\
			\
		}

		#define SSE2_RADIX4_DIF_4TWIDDLE_B(__tmp,__c0)\
		{\
			/* Do the p0,p8 combo: */\
			__asm	mov	esi, __c0\
			\
			__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */				__asm	movaps	xmm4,[ecx     ]	/* a[jt+p8 ] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p8 ] */\
			__asm	movaps	xmm6,[esi     ]	/* c0 */\
			__asm	movaps	xmm7,[esi+0x10]	/* s0 */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[jt   ] */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[jp   ] */\
			\
			__asm	mulpd	xmm0,xmm6		/* a[jt   ]*c0 */\
			__asm	mulpd	xmm1,xmm6		/* a[jp   ]*c0 */\
			__asm	mulpd	xmm2,xmm7		/* a[jt   ]*s0 */\
			__asm	mulpd	xmm3,xmm7		/* a[jp   ]*s0	xmm6,7 free */\
																		__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[jt+p8 ] */\
			__asm	addpd	xmm1,xmm2		/* xmm1 <- t2 */			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[jp+p8 ] */\
																		__asm	mulpd	xmm4,[esi+0x20] /* a[jt+p8 ]*c8 */\
			__asm	subpd	xmm0,xmm3		/* xmm0 <- t1 */			__asm	mulpd	xmm5,[esi+0x20]	/* a[jp+p8 ]*c8 */\
																		__asm	mulpd	xmm6,[esi+0x30]	/* a[jt+p8 ]*s8 */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy t1 */		__asm	mulpd	xmm7,[esi+0x30]	/* a[jp+p8 ]*s8 */\
																		__asm	addpd	xmm5,xmm6	    /* xmm5 <- it */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy t2 */		__asm	subpd	xmm4,xmm7		/* xmm4 <- rt    xmm6,7 free */\
			\
			__asm	addpd	xmm0,xmm4		/* ~t1 <- t1 +rt */\
			__asm	addpd	xmm1,xmm5		/* ~t2 <- t2 +it */\
			__asm	subpd	xmm2,xmm4		/* ~t3 <- t1 -rt */\
			__asm	subpd	xmm3,xmm5		/* ~t4 <- t2 -it	xmm4,5 free */\
			\
			/* Do the p4,12 combo: */\
			__asm	add	esi, 0x60\
			\
			__asm	movaps	xmm4,[edx     ]	/* a[jt+p12] */\
			__asm	movaps	xmm5,[edx+0x10]	/* a[jp+p12] */\
			__asm	movaps	xmm6,[edx     ]	/* xmm2 <- cpy a[jt+p12] */\
			__asm	movaps	xmm7,[edx+0x10]	/* xmm3 <- cpy a[jp+p12] */\
			\
			__asm	mulpd	xmm4,[esi     ]	/* a[jt+p12]*c12 */\
			__asm	mulpd	xmm5,[esi     ]	/* a[jp+p12]*c12 */\
			__asm	mulpd	xmm6,[esi+0x10]	/* a[jt+p12]*s12 */\
			__asm	mulpd	xmm7,[esi+0x10]	/* a[jp+p12]*s12 */\
			\
			__asm	mov	esi, __tmp	\
			__asm	addpd	xmm5,xmm6		/* xmm5 <- it */\
			__asm	subpd	xmm4,xmm7		/* xmm4 <- rt */\
			__asm	movaps	[esi+0x010],xmm5	/* store it */\
			__asm	movaps	[esi      ],xmm4	/* store rt */\
			\
			__asm	mov	esi, __c0\
			__asm	add	esi, 0x40\
			__asm	movaps	xmm4,[ebx     ]	/* a[jt+p4] */\
			__asm	movaps	xmm5,[ebx+0x10]	/* a[jp+p4] */\
			__asm	movaps	xmm6,[ebx     ]	/* xmm2 <- cpy a[jt+p4] */\
			__asm	movaps	xmm7,[ebx+0x10]	/* xmm3 <- cpy a[jp+p4] */\
			\
			__asm	mulpd	xmm4,[esi     ]	/* a[jt+p4]*c4 */\
			__asm	mulpd	xmm5,[esi     ]	/* a[jp+p4]*c4 */\
			__asm	mulpd	xmm6,[esi+0x10]	/* a[jt+p4]*s4 */\
			__asm	mulpd	xmm7,[esi+0x10]	/* a[jp+p4]*s4 */\
			__asm	mov	esi, __tmp	\
			__asm	addpd	xmm5,xmm6		/* xmm5 <- t6 */\
			__asm	subpd	xmm4,xmm7		/* xmm4 <- t5 	xmm6,7 free */\
			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy t6 */\
			__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy t5 */\
			\
			__asm	subpd	xmm4,[esi      ]	/* ~t7 <- t5 -rt */\
			__asm	subpd	xmm5,[esi+0x010]	/* ~t8 <- t6 -it */\
			__asm	addpd	xmm6,[esi      ]	/* ~t5 <- t5 +rt */\
			__asm	addpd	xmm7,[esi+0x010]	/* ~t6 <- t6 +it */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			/*\
			~t5 =t1 -t5;		~t1 =t1 +t5;										~t7 =t3 +t8;		~t3 =t3 -t8;\
			~t6 =t2 -t6;		~t2 =t2 +t6;										~t8 =t4 -t7;		~t4 =t4 +t7;\
			*/\
			__asm	subpd	xmm0,xmm6			/*~t5 */							__asm	subpd	xmm2,xmm5			/*~t3 */\
			__asm	subpd	xmm1,xmm7			/*~t6 */							__asm	subpd	xmm3,xmm4			/*~t8 */\
			__asm	movaps	[esi+0x040],xmm0	/* a[jt+p8 ] <- ~t5 */				__asm	movaps	[esi+0x020],xmm2	/* a[jt+p4 ] <- ~t3 */\
			__asm	movaps	[esi+0x050],xmm1	/* a[jp+p8 ] <- ~t6 */				__asm	movaps	[esi+0x070],xmm3	/* a[jp+p12] <- ~t8 */\
			__asm	addpd	xmm6,xmm6			/* 2*t5 */							__asm	addpd	xmm5,xmm5			/* 2*t8 */\
			__asm	addpd	xmm7,xmm7			/* 2*t6 */							__asm	addpd	xmm4,xmm4			/* 2*t7 */\
			__asm	addpd	xmm6,xmm0			/*~t1 */							__asm	addpd	xmm5,xmm2			/*~t7 */\
			__asm	addpd	xmm7,xmm1			/*~t2 */							__asm	addpd	xmm4,xmm3			/*~t4 */\
			__asm	movaps	[esi      ],xmm6	/* a[jt    ] <- ~t1 */				__asm	movaps	[esi+0x060],xmm5	/* a[jt+p12] <- ~t7 */\
			__asm	movaps	[esi+0x010],xmm7	/* a[jp    ] <- ~t2 */				__asm	movaps	[esi+0x030],xmm4	/* a[jp+p4 ] <- ~t4 */\
			\
		}

		/* DIF radix-4 subconvolution, with twiddles.
		The SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO variant means that we munge the second set of 4 output as follows:

			t3 =(t3-t4)*ISRT2;	t4 =(t3+t4)*ISRT2;
			t7 =(t7+t8)*ISRT2;	t8 =(t7-t8)*ISRT2;

		Assumes the 4 twiddles c0,4,8,12 are bit-reversed in memory, i.e. c0,8,4,12 are adjacent in memory starting at c0.
		*/
		/* Cost: 30 MOVapd, 28 ADD/SUBpd, 16 MULpd */
		#define SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(__add0, __add1, __add2, __add3, __tmp, __c0)\
		{\
			/* Do the p0,p8 combo: */\
			__asm	mov	eax, __add0\
			__asm	mov	ebx, __c0\
			__asm	mov	ecx, __add2\
			\
			__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */				__asm	movaps	xmm4,[ecx     ]	/* a[jt+p8 ] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p8 ] */\
			__asm	movaps	xmm6,[ebx     ]	/* c0 */\
			__asm	movaps	xmm7,[ebx+0x10]	/* s0 */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[jt   ] */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[jp   ] */\
			\
			__asm	mulpd	xmm0,xmm6		/* a[jt   ]*c0 */\
			__asm	mulpd	xmm1,xmm6		/* a[jp   ]*c0 */\
			__asm	mulpd	xmm2,xmm7		/* a[jt   ]*s0 */\
			__asm	mulpd	xmm3,xmm7		/* a[jp   ]*s0	xmm6,7 free */\
																		__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[jt+p8 ] */\
			__asm	addpd	xmm1,xmm2		/* xmm1 <- t2 */			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[jp+p8 ] */\
																		__asm	mulpd	xmm4,[ebx+0x20] /* a[jt+p8 ]*c8 */\
			__asm	subpd	xmm0,xmm3		/* xmm0 <- t1 */			__asm	mulpd	xmm5,[ebx+0x20]	/* a[jp+p8 ]*c8 */\
																		__asm	mulpd	xmm6,[ebx+0x30]	/* a[jt+p8 ]*s8 */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy t1 */		__asm	mulpd	xmm7,[ebx+0x30]	/* a[jp+p8 ]*s8 */\
																		__asm	addpd	xmm5,xmm6	    /* xmm5 <- it */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy t2 */		__asm	subpd	xmm4,xmm7		/* xmm4 <- rt    xmm6,7 free */\
			\
			__asm	addpd	xmm0,xmm4	/* ~t1 <- t1 +rt */\
			__asm	addpd	xmm1,xmm5	/* ~t2 <- t2 +it */\
			__asm	subpd	xmm2,xmm4	/* ~t3 <- t1 -rt */\
			__asm	subpd	xmm3,xmm5	/* ~t4 <- t2 -it	xmm4,5 free */\
			\
			/* Do the p4,12 combo: */\
			__asm	mov	eax, __add1\
			__asm	add	ebx, 0x60\
			__asm	mov	ecx, __add3\
			\
			__asm	movaps	xmm4,[ecx     ]	/* a[jt+p12] */\
			__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p12] */\
			__asm	movaps	xmm6,[ecx     ]	/* xmm2 <- cpy a[jt+p12] */\
			__asm	movaps	xmm7,[ecx+0x10]	/* xmm3 <- cpy a[jp+p12] */\
			\
			__asm	mulpd	xmm4,[ebx     ]	/* a[jt+p12]*c12 */\
			__asm	mulpd	xmm5,[ebx     ]	/* a[jp+p12]*c12 */\
			__asm	mulpd	xmm6,[ebx+0x10]	/* a[jt+p12]*s12 */\
			__asm	mulpd	xmm7,[ebx+0x10]	/* a[jp+p12]*s12 */\
			\
			__asm	mov	edx, __tmp	\
			\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- it */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- rt */\
			__asm	movaps	[edx+0x010],xmm5	/* tmp store it */\
			__asm	movaps	[edx      ],xmm4	/* tmp store rt */\
			\
			__asm	sub	ebx, 0x20\
			__asm	movaps	xmm4,[eax     ]	/* a[jt+p4] */\
			__asm	movaps	xmm5,[eax+0x10]	/* a[jp+p4] */\
			__asm	movaps	xmm6,[eax     ]	/* xmm2 <- cpy a[jt+p4] */\
			__asm	movaps	xmm7,[eax+0x10]	/* xmm3 <- cpy a[jp+p4] */\
			\
			__asm	mulpd	xmm4,[ebx     ]	/* a[jt+p4]*c4 */\
			__asm	mulpd	xmm5,[ebx     ]	/* a[jp+p4]*c4 */\
			__asm	mulpd	xmm6,[ebx+0x10]	/* a[jt+p4]*s4 */\
			__asm	mulpd	xmm7,[ebx+0x10]	/* a[jp+p4]*s4 */\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- t6 */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- t5 	xmm6,7 free */\
			__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy t6 */\
			__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy t5 */\
			\
			__asm	subpd	xmm4,[edx      ]	/* ~t7 <- t5 -rt */\
			__asm	subpd	xmm5,[edx+0x010]	/* ~t8 <- t6 -it */\
			__asm	addpd	xmm6,[edx      ]	/* ~t5 <- t5 +rt */\
			__asm	addpd	xmm7,[edx+0x010]	/* ~t6 <- t6 +it */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			/*\
			~t5 =t1 -t5;		~t1 =t1 +t5;										~t7 =t3 +t8;		~t3 =t3 -t8;\
			~t6 =t2 -t6;		~t2 =t2 +t6;										~t8 =t4 -t7;		~t4 =t4 +t7;\
			*/\
			__asm	subpd	xmm0,xmm6	/*~t5 */									__asm	mov	ecx, isrt2\
			__asm	subpd	xmm1,xmm7	/*~t6 */									\
			__asm	movaps	[edx+0x040],xmm0	/* a[jt+p8 ] <- ~t5 */				__asm	subpd	xmm2,xmm5	/*~t3 */\
			__asm	movaps	[edx+0x050],xmm1	/* a[jp+p8 ] <- ~t6 */				__asm	subpd	xmm3,xmm4	/*~t8 */\
			__asm	addpd	xmm6,xmm6	/* 2*t5 */									__asm	addpd	xmm5,xmm5	/* 2*t8 */\
			__asm	addpd	xmm7,xmm7	/* 2*t6 */									__asm	addpd	xmm4,xmm4	/* 2*t7 */\
			__asm	addpd	xmm6,xmm0	/*~t1 */									__asm	addpd	xmm5,xmm2	/*~t7 */\
			__asm	addpd	xmm7,xmm1	/*~t2 */									__asm	addpd	xmm4,xmm3	/*~t4 */\
			__asm	movaps	[edx      ],xmm6	/* a[jt    ] <- ~t1 */\
			__asm	movaps	[edx+0x010],xmm7	/* a[jp    ] <- ~t2 */\
			/*\
			t3 =(t3-t4)*ISRT2;	t4 =(t3+t4)*ISRT2;\
			t7 =(t7-t8)*ISRT2;	t8 =(t7+t8)*ISRT2;\
			*/\
			__asm	movaps	xmm0,[ecx]	/* ISRT2 */\
			__asm	movaps	xmm6,xmm2	/* cpy t3 */\
			__asm	movaps	xmm7,xmm5	/* cpy t7 */\
			__asm	subpd	xmm2,xmm4	/* 3-4*/\
			__asm	subpd	xmm5,xmm3	/* 7-8*/\
			__asm	addpd	xmm6,xmm4	/* 3+4*/\
			__asm	addpd	xmm7,xmm3	/* 7+8*/\
			__asm	mulpd	xmm2,xmm0	/* (3-4)*ISRT2 */\
			__asm	mulpd	xmm5,xmm0	/* (7-8)*ISRT2 */\
			__asm	mulpd	xmm6,xmm0	/* (3+4)*ISRT2 */\
			__asm	mulpd	xmm7,xmm0	/* (7+8)*ISRT2 */\
			__asm	movaps	[edx+0x020],xmm2	/* a[jt+p4 ] <- ~t3 */\
			__asm	movaps	[edx+0x060],xmm5	/* a[jp+p12] <- ~t7 */\
			__asm	movaps	[edx+0x030],xmm6	/* a[jp+p4 ] <- ~t4 */\
			__asm	movaps	[edx+0x070],xmm7	/* a[jt+p12] <- ~t8 */\
		}

		/* Variant of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO, which assumes that the base offset [not used inside the macro] is in edi,
		addresses __add0,1,2,3 enter in eax,ebx,ecx,edx, and that the contents of these 5 regs must not be changed on return.
		Also assumes the 4 twiddles c0,4,8,12 are bit-reversed in memory, i.e. c0,8,4,12 are adjacent in memory starting at c0.
		*/
		#define SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(__tmp, __c0)\
		{\
			/* Do the p0,p8 combo: */\
			__asm	mov	esi, __c0\
			\
			__asm	movaps	xmm0,[eax     ]	/* a[jt   ] */				__asm	movaps	xmm4,[ecx     ]	/* a[jt+p8 ] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[jp   ] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p8 ] */\
			__asm	movaps	xmm6,[esi     ]	/* c0 */\
			__asm	movaps	xmm7,[esi+0x10]	/* s0 */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[jt   ] */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[jp   ] */\
			\
			__asm	mulpd	xmm0,xmm6		/* a[jt   ]*c0 */\
			__asm	mulpd	xmm1,xmm6		/* a[jp   ]*c0 */\
			__asm	mulpd	xmm2,xmm7		/* a[jt   ]*s0 */\
			__asm	mulpd	xmm3,xmm7		/* a[jp   ]*s0	xmm6,7 free */\
																		__asm	add	esi, 0x20\
																		__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[jt+p8 ] */\
			__asm	addpd	xmm1,xmm2		/* xmm1 <- t2 */			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[jp+p8 ] */\
																		__asm	mulpd	xmm4,[esi     ] /* a[jt+p8 ]*c8 */\
			__asm	subpd	xmm0,xmm3		/* xmm0 <- t1 */			__asm	mulpd	xmm5,[esi     ]	/* a[jp+p8 ]*c8 */\
																		__asm	mulpd	xmm6,[esi+0x10]	/* a[jt+p8 ]*s8 */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy t1 */		__asm	mulpd	xmm7,[esi+0x10]	/* a[jp+p8 ]*s8 */\
																		__asm	addpd	xmm5,xmm6	    /* xmm5 <- it */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy t2 */		__asm	subpd	xmm4,xmm7		/* xmm4 <- rt    xmm6,7 free */\
			\
			__asm	addpd	xmm0,xmm4	/* ~t1 <- t1 +rt */\
			__asm	addpd	xmm1,xmm5	/* ~t2 <- t2 +it */\
			__asm	subpd	xmm2,xmm4	/* ~t3 <- t1 -rt */\
			__asm	subpd	xmm3,xmm5	/* ~t4 <- t2 -it	xmm4,5 free */\
			\
			/* Do the p4,12 combo: */\
			__asm	add	esi, 0x40\
			\
			__asm	movaps	xmm4,[edx     ]	/* a[jt+p12] */\
			__asm	movaps	xmm5,[edx+0x10]	/* a[jp+p12] */\
			__asm	movaps	xmm6,[edx     ]	/* xmm2 <- cpy a[jt+p12] */\
			__asm	movaps	xmm7,[edx+0x10]	/* xmm3 <- cpy a[jp+p12] */\
			\
			__asm	mulpd	xmm4,[esi     ]	/* a[jt+p12]*c12 */\
			__asm	mulpd	xmm5,[esi     ]	/* a[jp+p12]*c12 */\
			__asm	mulpd	xmm6,[esi+0x10]	/* a[jt+p12]*s12 */\
			__asm	mulpd	xmm7,[esi+0x10]	/* a[jp+p12]*s12 */\
			\
			__asm	mov	esi, __tmp	\
			\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- it */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- rt */\
			__asm	movaps	[esi+0x010],xmm5	/* tmp store it */\
			__asm	movaps	[esi      ],xmm4	/* tmp store rt */\
			\
			__asm	mov	esi, __c0\
			__asm	add	esi, 0x40\
			__asm	movaps	xmm4,[ebx     ]	/* a[jt+p4] */\
			__asm	movaps	xmm5,[ebx+0x10]	/* a[jp+p4] */\
			__asm	movaps	xmm6,[ebx     ]	/* xmm2 <- cpy a[jt+p4] */\
			__asm	movaps	xmm7,[ebx+0x10]	/* xmm3 <- cpy a[jp+p4] */\
			\
			__asm	mulpd	xmm4,[esi     ]	/* a[jt+p4]*c4 */\
			__asm	mulpd	xmm5,[esi     ]	/* a[jp+p4]*c4 */\
			__asm	mulpd	xmm6,[esi+0x10]	/* a[jt+p4]*s4 */\
			__asm	mulpd	xmm7,[esi+0x10]	/* a[jp+p4]*s4 */\
			__asm	mov	esi, __tmp	\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- t6 */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- t5 	xmm6,7 free */\
			__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy t6 */\
			__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy t5 */\
			\
			__asm	subpd	xmm4,[esi      ]	/* ~t7 <- t5 -rt */\
			__asm	subpd	xmm5,[esi+0x010]	/* ~t8 <- t6 -it */\
			__asm	addpd	xmm6,[esi      ]	/* ~t5 <- t5 +rt */\
			__asm	addpd	xmm7,[esi+0x010]	/* ~t6 <- t6 +it */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			/*\
			~t5 =t1 -t5;		~t1 =t1 +t5;										~t7 =t3 +t8;		~t3 =t3 -t8;\
			~t6 =t2 -t6;		~t2 =t2 +t6;										~t8 =t4 -t7;		~t4 =t4 +t7;\
			*/\
			__asm	subpd	xmm0,xmm6	/*~t5 */\
			__asm	subpd	xmm1,xmm7	/*~t6 */									\
			__asm	movaps	[esi+0x040],xmm0	/* a[jt+p8 ] <- ~t5 */				__asm	subpd	xmm2,xmm5	/*~t3 */\
			__asm	movaps	[esi+0x050],xmm1	/* a[jp+p8 ] <- ~t6 */				__asm	subpd	xmm3,xmm4	/*~t8 */\
			__asm	addpd	xmm6,xmm6	/* 2*t5 */									__asm	addpd	xmm5,xmm5	/* 2*t8 */\
			__asm	addpd	xmm7,xmm7	/* 2*t6 */									__asm	addpd	xmm4,xmm4	/* 2*t7 */\
			__asm	addpd	xmm6,xmm0	/*~t1 */									__asm	addpd	xmm5,xmm2	/*~t7 */\
			__asm	addpd	xmm7,xmm1	/*~t2 */									__asm	addpd	xmm4,xmm3	/*~t4 */\
			__asm	movaps	[esi      ],xmm6	/* a[jt    ] <- ~t1 */\
			__asm	movaps	[esi+0x010],xmm7	/* a[jp    ] <- ~t2 */\
			/*\
			t3 =(t3-t4)*ISRT2;	t4 =(t3+t4)*ISRT2;\
			t7 =(t7-t8)*ISRT2;	t8 =(t7+t8)*ISRT2;\
			*/\
			__asm	mov	esi, isrt2\
			__asm	movaps	xmm0,[esi]	/* ISRT2 */\
			__asm	movaps	xmm6,xmm2	/* cpy t3 */\
			__asm	movaps	xmm7,xmm5	/* cpy t7 */\
			__asm	subpd	xmm2,xmm4	/* 3-4*/\
			__asm	subpd	xmm5,xmm3	/* 7-8*/\
			__asm	addpd	xmm6,xmm4	/* 3+4*/\
			__asm	addpd	xmm7,xmm3	/* 7+8*/\
			__asm	mov	esi, __tmp	\
			__asm	mulpd	xmm2,xmm0	/* (3-4)*ISRT2 */\
			__asm	mulpd	xmm5,xmm0	/* (7-8)*ISRT2 */\
			__asm	mulpd	xmm6,xmm0	/* (3+4)*ISRT2 */\
			__asm	mulpd	xmm7,xmm0	/* (7+8)*ISRT2 */\
			__asm	movaps	[esi+0x020],xmm2	/* a[jt+p4 ] <- ~t3 */\
			__asm	movaps	[esi+0x060],xmm5	/* a[jp+p12] <- ~t7 */\
			__asm	movaps	[esi+0x030],xmm6	/* a[jp+p4 ] <- ~t4 */\
			__asm	movaps	[esi+0x070],xmm7	/* a[jt+p12] <- ~t8 */\
		}

		/* Combine two radix-4 DIF subtransforms. The SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO variant means the
		ISRT2 stuff is already done on the second input quartet, i.e. both radix-2 combos here look identical.
		*/
		#define SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(__r0,__r1,__r2,__r3,__r4,__r5,__r6,__r7)\
		{\
			__asm	mov	eax, __r0\
			__asm	mov	ebx, __r4\
			__asm	mov	ecx, __r2\
			__asm	mov	edx, __r6\
			\
			__asm	movaps	xmm0,[eax     ]	/* t0 */			__asm	movaps	xmm4,[ecx     ]	/* t4 */\
			__asm	movaps	xmm1,[eax+0x10]	/* t1 */			__asm	movaps	xmm5,[ecx+0x10]	/* t5 */\
			__asm	movaps	xmm2,[ebx     ]	/* cpy t8 */		__asm	movaps	xmm7,[edx+0x10]	/* td */\
			__asm	movaps	xmm3,[ebx+0x10]	/* cpy t9 */		__asm	movaps	xmm6,[edx     ]	/* tc */\
			__asm	subpd	xmm0,xmm2		/* t8 = 0-8 */		__asm	subpd	xmm4,xmm7		/* t4 = 4-d */\
			__asm	subpd	xmm1,xmm3		/* t9 = 1-9 */		__asm	subpd	xmm5,xmm6		/* td = 5-c */\
			__asm	addpd	xmm2,xmm2		/*      2*8 */		__asm	addpd	xmm7,xmm7		/*      2*d */\
			__asm	addpd	xmm3,xmm3		/*      2*9 */		__asm	addpd	xmm6,xmm6		/*      2*c */\
			__asm	addpd	xmm2,xmm0		/* t0 = 0+8 */		__asm	addpd	xmm7,xmm4		/* tc = 4+d */\
			__asm	addpd	xmm3,xmm1		/* t1 = 1+9 */		__asm	addpd	xmm6,xmm5		/* t5 = 5+c */\
			\
			__asm	movaps	[ebx     ],xmm0	/* t8 */			__asm	movaps	[ecx     ],xmm4	/* t4 */\
			__asm	movaps	[ebx+0x10],xmm1	/* t9 */			__asm	movaps	[edx+0x10],xmm5	/* td */\
			__asm	movaps	[eax     ],xmm2	/* t0 */			__asm	movaps	[edx     ],xmm7	/* tc */\
			__asm	movaps	[eax+0x10],xmm3	/* t1 */			__asm	movaps	[ecx+0x10],xmm6	/* t5 */\
			\
			__asm	mov	eax, __r1\
			__asm	mov	ebx, __r5\
			__asm	mov	ecx, __r3\
			__asm	mov	edx, __r7\
			\
			__asm	movaps	xmm0,[eax     ]	/* t2 */			__asm	movaps	xmm4,[ecx     ]	/* t6 */\
			__asm	movaps	xmm1,[eax+0x10]	/* t3 */			__asm	movaps	xmm5,[ecx+0x10]	/* t7 */\
			__asm	movaps	xmm2,[ebx     ]	/* cpy ta */		__asm	movaps	xmm7,[edx+0x10]	/* tf */\
			__asm	movaps	xmm3,[ebx+0x10]	/* cpy tb */		__asm	movaps	xmm6,[edx     ]	/* te */\
			__asm	subpd	xmm0,xmm2		/* ta = 2-a */		__asm	subpd	xmm4,xmm7		/* t6 = 6-f */\
			__asm	subpd	xmm1,xmm3		/* tb = 3-b */		__asm	subpd	xmm5,xmm6		/* tf = 7-e */\
			__asm	addpd	xmm2,xmm2		/*      2*a */		__asm	addpd	xmm7,xmm7		/*      2*f */\
			__asm	addpd	xmm3,xmm3		/*      2*b */		__asm	addpd	xmm6,xmm6		/*      2*e */\
			__asm	addpd	xmm2,xmm0		/* t2 = 2+a */		__asm	addpd	xmm7,xmm4		/* te = 6+f */\
			__asm	addpd	xmm3,xmm1		/* t3 = 3+b */		__asm	addpd	xmm6,xmm5		/* t7 = 7+e */\
			\
			__asm	movaps	[ebx     ],xmm0	/* ta */			__asm	movaps	[ecx     ],xmm4	/* t6 */\
			__asm	movaps	[ebx+0x10],xmm1	/* tb */			__asm	movaps	[edx+0x10],xmm5	/* tf */\
			__asm	movaps	[eax     ],xmm2	/* t2 */			__asm	movaps	[edx     ],xmm7	/* te */\
			__asm	movaps	[eax+0x10],xmm3	/* t3 */			__asm	movaps	[ecx+0x10],xmm6	/* t7 */\
		}

		/* 2 Variants of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO, which assume that the base offset [not used inside the macro] is in edi,
		main-array addresses __add1,3,5,7 are in eax,ebx,ecx,edx, and that the contents of these 5 regs must not be changed on return.
		This leaves esi for local address computation, requiring us to assume that the local addresses r0,1,2,3,4,5,6,7
		occupy consecutive memory locations.
		A-variant writes outputs to main-array addresses __add0-7; B-variant writes back to local addresses r0,1,2,3,4,5,6,7.
		*/
		#define SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(__r0)\
		{\
			__asm	mov	esi, __r0\
			\
			__asm	movaps	xmm0,[esi+0x20]	/* t2 */			__asm	movaps	xmm4,[esi+0x60]	/* t6 */\
			__asm	movaps	xmm1,[esi+0x30]	/* t3 */			__asm	movaps	xmm5,[esi+0x70]	/* t7 */\
			__asm	movaps	xmm2,[esi+0xa0]	/* cpy ta */		__asm	movaps	xmm7,[esi+0xf0]	/* tf */\
			__asm	movaps	xmm3,[esi+0xb0]	/* cpy tb */		__asm	movaps	xmm6,[esi+0xe0]	/* te */\
			__asm	subpd	xmm0,xmm2		/* ta = 2-a */		__asm	subpd	xmm4,xmm7		/* t6 = 6-f */\
			__asm	subpd	xmm1,xmm3		/* tb = 3-b */		__asm	subpd	xmm5,xmm6		/* tf = 7-e */\
			__asm	addpd	xmm2,xmm2		/*      2*a */		__asm	addpd	xmm7,xmm7		/*      2*f */\
			__asm	addpd	xmm3,xmm3		/*      2*b */		__asm	addpd	xmm6,xmm6		/*      2*e */\
			__asm	addpd	xmm2,xmm0		/* t2 = 2+a */		__asm	addpd	xmm7,xmm4		/* te = 6+f */\
			__asm	addpd	xmm3,xmm1		/* t3 = 3+b */		__asm	addpd	xmm6,xmm5		/* t7 = 7+e */\
			\
			__asm	mov	eax, ecx	/* &a[j1+p5] */\
			__asm	sub	eax, edi	/* &a[j1+p4] */\
			__asm	mov	ebx, edx	/* &a[j1+p7] */\
			__asm	sub	ebx, edi	/* &a[j1+p6] */\
			\
			__asm	movaps	[ecx     ],xmm0	/* ta */			__asm	movaps	[ebx     ],xmm4	/* t6 */\
			__asm	movaps	[ecx+0x10],xmm1	/* tb */			__asm	movaps	[edx+0x10],xmm5	/* tf */\
			__asm	movaps	[eax     ],xmm2	/* t2 */			__asm	movaps	[edx     ],xmm7	/* te */\
			__asm	movaps	[eax+0x10],xmm3	/* t3 */			__asm	movaps	[ebx+0x10],xmm6	/* t7 */\
			\
			__asm	mov	edi, p4	/* Can't simply do p4 = (p1<<2) in-place, due to array-padding scheme */\
			__asm	shl	edi,  3		/* 8-bytes for array-of-doubles */\
			__asm	sub	eax, edi	/* &a[j1+p4] -> &a[j1+p0] */\
			__asm	sub	ebx, edi	/* &a[j1+p6] -> &a[j1+p2] */\
			__asm	sub	ecx, edi	/* &a[j1+p5] -> &a[j1+p1] */\
			__asm	sub	edx, edi	/* &a[j1+p7] -> &a[j1+p3] */\
			__asm	shr	edi,  2	/* p1 = (p4>>2) works, because any excess padding [i.e. present in p4, but not p1] gets shifted off */\
			\
			__asm	movaps	xmm0,[esi     ]	/* t0 */			__asm	movaps	xmm4,[esi+0x40]	/* t4 */\
			__asm	movaps	xmm1,[esi+0x10]	/* t1 */			__asm	movaps	xmm5,[esi+0x50]	/* t5 */\
			__asm	movaps	xmm2,[esi+0x80]	/* cpy t8 */		__asm	movaps	xmm7,[esi+0xd0]	/* td */\
			__asm	movaps	xmm3,[esi+0x90]	/* cpy t9 */		__asm	movaps	xmm6,[esi+0xc0]	/* tc */\
			__asm	subpd	xmm0,xmm2		/* t8 = 0-8 */		__asm	subpd	xmm4,xmm7		/* t4 = 4-d */\
			__asm	subpd	xmm1,xmm3		/* t9 = 1-9 */		__asm	subpd	xmm5,xmm6		/* td = 5-c */\
			__asm	addpd	xmm2,xmm2		/*      2*8 */		__asm	addpd	xmm7,xmm7		/*      2*d */\
			__asm	addpd	xmm3,xmm3		/*      2*9 */		__asm	addpd	xmm6,xmm6		/*      2*c */\
			__asm	addpd	xmm2,xmm0		/* t0 = 0+8 */		__asm	addpd	xmm7,xmm4		/* tc = 4+d */\
			__asm	addpd	xmm3,xmm1		/* t1 = 1+9 */		__asm	addpd	xmm6,xmm5		/* t5 = 5+c */\
			\
			__asm	movaps	[ecx     ],xmm0	/* t8 */			__asm	movaps	[ebx     ],xmm4	/* t4 */\
			__asm	movaps	[ecx+0x10],xmm1	/* t9 */			__asm	movaps	[edx+0x10],xmm5	/* td */\
			__asm	movaps	[eax     ],xmm2	/* t0 */			__asm	movaps	[edx     ],xmm7	/* tc */\
			__asm	movaps	[eax+0x10],xmm3	/* t1 */			__asm	movaps	[ebx+0x10],xmm6	/* t5 */\
		}

		#define SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(__r0)\
		{\
			__asm	mov	esi, __r0\
			\
			__asm	movaps	xmm0,[esi     ]	/* t0 */			__asm	movaps	xmm4,[esi+0x40]	/* t4 */\
			__asm	movaps	xmm1,[esi+0x10]	/* t1 */			__asm	movaps	xmm5,[esi+0x50]	/* t5 */\
			__asm	movaps	xmm2,[esi+0x80]	/* cpy t8 */		__asm	movaps	xmm7,[esi+0xd0]	/* td */\
			__asm	movaps	xmm3,[esi+0x90]	/* cpy t9 */		__asm	movaps	xmm6,[esi+0xc0]	/* tc */\
			__asm	subpd	xmm0,xmm2		/* t8 = 0-8 */		__asm	subpd	xmm4,xmm7		/* t4 = 4-d */\
			__asm	subpd	xmm1,xmm3		/* t9 = 1-9 */		__asm	subpd	xmm5,xmm6		/* td = 5-c */\
			__asm	addpd	xmm2,xmm2		/*      2*8 */		__asm	addpd	xmm7,xmm7		/*      2*d */\
			__asm	addpd	xmm3,xmm3		/*      2*9 */		__asm	addpd	xmm6,xmm6		/*      2*c */\
			__asm	addpd	xmm2,xmm0		/* t0 = 0+8 */		__asm	addpd	xmm7,xmm4		/* tc = 4+d */\
			__asm	addpd	xmm3,xmm1		/* t1 = 1+9 */		__asm	addpd	xmm6,xmm5		/* t5 = 5+c */\
			\
			__asm	movaps	[esi+0x80],xmm0	/* t8 */			__asm	movaps	[esi+0x40],xmm4	/* t4 */\
			__asm	movaps	[esi+0x90],xmm1	/* t9 */			__asm	movaps	[esi+0xd0],xmm5	/* td */\
			__asm	movaps	[esi     ],xmm2	/* t0 */			__asm	movaps	[esi+0xc0],xmm7	/* tc */\
			__asm	movaps	[esi+0x10],xmm3	/* t1 */			__asm	movaps	[esi+0x50],xmm6	/* t5 */\
			\
			__asm	movaps	xmm0,[esi+0x20]	/* t2 */			__asm	movaps	xmm4,[esi+0x60]	/* t6 */\
			__asm	movaps	xmm1,[esi+0x30]	/* t3 */			__asm	movaps	xmm5,[esi+0x70]	/* t7 */\
			__asm	movaps	xmm2,[esi+0xa0]	/* cpy ta */		__asm	movaps	xmm7,[esi+0xf0]	/* tf */\
			__asm	movaps	xmm3,[esi+0xb0]	/* cpy tb */		__asm	movaps	xmm6,[esi+0xe0]	/* te */\
			__asm	subpd	xmm0,xmm2		/* ta = 2-a */		__asm	subpd	xmm4,xmm7		/* t6 = 6-f */\
			__asm	subpd	xmm1,xmm3		/* tb = 3-b */		__asm	subpd	xmm5,xmm6		/* tf = 7-e */\
			__asm	addpd	xmm2,xmm2		/*      2*a */		__asm	addpd	xmm7,xmm7		/*      2*f */\
			__asm	addpd	xmm3,xmm3		/*      2*b */		__asm	addpd	xmm6,xmm6		/*      2*e */\
			__asm	addpd	xmm2,xmm0		/* t2 = 2+a */		__asm	addpd	xmm7,xmm4		/* te = 6+f */\
			__asm	addpd	xmm3,xmm1		/* t3 = 3+b */		__asm	addpd	xmm6,xmm5		/* t7 = 7+e */\
			\
			__asm	movaps	[esi+0xa0],xmm0	/* ta */			__asm	movaps	[esi+0x60],xmm4	/* t6 */\
			__asm	movaps	[esi+0xb0],xmm1	/* tb */			__asm	movaps	[esi+0xf0],xmm5	/* tf */\
			__asm	movaps	[esi+0x20],xmm2	/* t2 */			__asm	movaps	[esi+0xe0],xmm7	/* te */\
			__asm	movaps	[esi+0x30],xmm3	/* t3 */			__asm	movaps	[esi+0x70],xmm6	/* t7 */\
		}

		/* Twiddleless version of SSE2_RADIX8_DIF_TWIDDLE. Inputs enter in memory locations __r0 + [__i1,__i2,__i3,__i4,__i5,__i6,__i7],;
		where r0 is a memory address and the i's are LITERAL [BYTE] OFFSETS. Outputs go into memory locations __o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7, assumed disjoint with inputs:\
		*/
		#define SSE2_RADIX8_DIF_0TWIDDLE(__r0, __i1,__i2,__i3,__i4,__i5,__i6,__i7, __o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7, __isrt2)\
		{\
			/* Can't get these via simple load-one-and-shift-as-needed due to array padding scheme */\
			__asm	mov	eax, __r0\
			__asm	mov	ebx, __i2\
			__asm	mov	ecx, __i4\
			__asm	mov	edx, __i6\
			__asm	add ebx, eax\
			__asm	add	ecx, eax\
			__asm	add	edx, eax\
			\
			/* Do the p0,p4 combo: */\
			\
			__asm	movaps	xmm0,[eax     ]	/* a[j1+p0] */				__asm	movaps	xmm4,[ecx     ]	/* a[j1+p4] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[j2+p0] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[j2+p4] */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[j1   ] */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[j2   ] */\
			\
			__asm	addpd	xmm0,xmm4		/* ~t1 <- t1 +rt */\
			__asm	addpd	xmm1,xmm5		/* ~t2 <- t2 +it */\
			__asm	subpd	xmm2,xmm4		/* ~t3 <- t1 -rt */\
			__asm	subpd	xmm3,xmm5		/* ~t4 <- t2 -it	xmm4,5 free */\
			\
			/* Do the p2,6 combo: */\
			\
			__asm	movaps	xmm4,[ebx     ]	/* a[jt+p2] */\
			__asm	movaps	xmm5,[ebx+0x10]	/* a[jp+p2] */\
			__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[jt+p2] */\
			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[jp+p2] */\
			\
			__asm	addpd	xmm4,[edx      ]	/* ~t4 <- t4 +rt */\
			__asm	addpd	xmm5,[edx+0x010]	/* ~t5 <- t5 +it */\
			__asm	subpd	xmm6,[edx      ]	/* ~t6 <- t4 -rt */\
			__asm	subpd	xmm7,[edx+0x010]	/* ~t7 <- t5 -it */\
			\
			/* Finish radix-4 butterfly and store results back into input-array slots: */\
			/*\
			~t4 =t0 -t4;		~t0 =t0 +t4;							~t6 =t2 +t7;		~t2 =t2 -t7;\
			~t5 =t1 -t5;		~t1 =t1 +t5;							~t7 =t3 -t6;		~t3 =t3 +t6;\
			*/\
			__asm	subpd	xmm0,xmm4	/*~t4 */						__asm	subpd	xmm2,xmm7	/*~t2 */\
			__asm	subpd	xmm1,xmm5	/*~t5 */						__asm	subpd	xmm3,xmm6	/*~t7 */\
			__asm	movaps	[ecx      ],xmm0	/* a[jt+p4] <- ~t4 */	__asm	movaps	[ebx      ],xmm2	/* a[jt+p2] <- ~t2 */\
			__asm	movaps	[ecx+0x010],xmm1	/* a[jp+p4] <- ~t5 */	__asm	movaps	[edx+0x010],xmm3	/* a[jp+p6] <- ~t7 */\
			__asm	addpd	xmm4,xmm4	/* 2*t4 */						__asm	addpd	xmm7,xmm7	/* 2*t7 */\
			__asm	addpd	xmm5,xmm5	/* 2*t5 */						__asm	addpd	xmm6,xmm6	/* 2*t6 */\
			__asm	addpd	xmm4,xmm0	/*~t0 */						__asm	addpd	xmm7,xmm2	/*~t6 */\
			__asm	addpd	xmm5,xmm1	/*~t1 */						__asm	addpd	xmm6,xmm3	/*~t3 */\
			__asm	movaps	[eax      ],xmm4	/* a[jt+p0] <- ~t0 */	__asm	movaps	[edx      ],xmm7	/* a[jt+p6] <- ~t6 */\
			__asm	movaps	[eax+0x010],xmm5	/* a[jp+p0] <- ~t1 */	__asm	movaps	[ebx+0x010],xmm6	/* a[jp+p2] <- ~t3 */\
			\
			/*******************************************************/\
			/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1): Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\
			/*******************************************************/\
			__asm	mov	ebx, __i3\
			__asm	mov	ecx, __i5\
			__asm	mov	edx, __i7\
			__asm	add ebx, eax\
			__asm	add	ecx, eax\
			__asm	add	edx, eax\
			__asm	add	eax, __i1\
			\
			/* Do the p1,p5 combo: */\
			\
			__asm	movaps	xmm0,[eax     ]	/* a[jt+p1] */				__asm	movaps	xmm4,[ecx     ]	/* a[jt+p5] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[jp+p1] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p5] */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[jt+p1] */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[jp+p1] */\
			\
			__asm	addpd	xmm0,xmm4	/* ~t8 <- t8 +rt */\
			__asm	addpd	xmm1,xmm5	/* ~t9 <- t9 +it */\
			__asm	subpd	xmm2,xmm4	/* ~ta <- t8 -rt */\
			__asm	subpd	xmm3,xmm5	/* ~tb <- t9 -it	xmm4,5 free */\
			\
			/* Do the p3,7 combo: */\
			\
			__asm	movaps	xmm4,[ebx     ]	/* a[jt+p3] */\
			__asm	movaps	xmm5,[ebx+0x10]	/* a[jp+p3] */\
			__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[jt+p3] */\
			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[jp+p3] */\
			\
			__asm	subpd	xmm4,[edx      ]	/* ~te <- tc -rt */\
			__asm	subpd	xmm5,[edx+0x010]	/* ~tf <- td -it */\
			__asm	addpd	xmm6,[edx      ]	/* ~tc <- tc +rt */\
			__asm	addpd	xmm7,[edx+0x010]	/* ~td <- td +it */\
			\
			/* Finish radix-4 butterfly and store just the 1st of the 8 outputs into output-array slots: */\
			/*\
			~tc =t8 -tc;		~t8 =t8 +tc;							~te =ta +tf;		~ta =ta -tf;\
			~td =t9 -td;		~t9 =t9 +td;							~tf =tb -te;		~tb =tb +te;\
			*/\
			__asm	subpd	xmm0,xmm6	/*~tc */\
			__asm	subpd	xmm1,xmm7	/*~td */\
		/*	__asm	movaps	[esi+0x040],xmm0	// a[jt+p5] <-~tc */	__asm	subpd	xmm2,xmm5	/*~ta */\
		/*	__asm	movaps	[esi+0x050],xmm1	// a[jp+p5] <-~td */	__asm	subpd	xmm3,xmm4	/*~tf */\
			__asm	addpd	xmm6,xmm6	/* 2*tc */						__asm	addpd	xmm5,xmm5	/* 2*tf */\
			__asm	addpd	xmm7,xmm7	/* 2*td */						__asm	addpd	xmm4,xmm4	/* 2*te */\
			__asm	addpd	xmm6,xmm0	/*~t8 */						__asm	addpd	xmm5,xmm2	/*~te */\
			__asm	addpd	xmm7,xmm1	/*~t9 */						__asm	addpd	xmm4,xmm3	/*~tb */\
			__asm	movaps	[eax     ],xmm6	/* a[jt+p1] <-~t8 */\
			__asm	movaps	[eax+0x10],xmm7	/* a[jp+p1] <-~t9; xmm6,7 FREE */\
			/*\
			ta =(ta-tb)*ISRT2;	tb =(ta+tb)*ISRT2;\
			te =(te-tf)*ISRT2;	tf =(te+tf)*ISRT2;\
			*/\
			__asm	movaps	xmm6,xmm2	/* xmm6 <- cpy ta */\
			__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy te */\
			__asm	subpd	xmm2,xmm4	/* ta-tb*/\
			__asm	subpd	xmm5,xmm3	/* te-tf*/\
			__asm	addpd	xmm4,xmm6	/* ta+tb*/\
			__asm	addpd	xmm3,xmm7	/* te+tf*/\
			\
			__asm	mov	esi, __isrt2\
			__asm	movaps	xmm6,[esi]	/* ISRT2 */\
			__asm	mulpd	xmm2,xmm6	/* ~ta = (ta-tb)*ISRT2 */\
			__asm	mulpd	xmm5,xmm6	/* ~te = (te-tf)*ISRT2 */\
			__asm	mulpd	xmm4,xmm6	/* ~tb = (ta+tb)*ISRT2 */\
			__asm	mulpd	xmm3,xmm6	/* ~tf = (te+tf)*ISRT2 */\
		/*	__asm	movaps	[esi+0x020],xmm2	// a[jt+p3] <-~ta */\
		/*	__asm	movaps	[esi+0x060],xmm5	// a[jp+p7] <-~te */\
		/*	__asm	movaps	[esi+0x030],xmm4	// a[jp+p3] <-~tb */\
		/*	__asm	movaps	[esi+0x070],xmm3	// a[jt+p7] <-~tf */\
			\
			/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs to temp array: */\
		/***** t0,1,2,3,4,5,6,7 in __i0,2,4,6 [eax,ebx,ecx,edx] *****/\
		/***** t8,9,a,b,c,d,e,f in __i1,xmm2,4,0,1,5,3; xmm6,7 FREE */\
			__asm	mov	edi, __r0\
			__asm	mov	esi, __r0\
			__asm	add esi, 0x10\
			\
			__asm	mov	eax, __o4\
			__asm	mov	ebx, __o5\
			__asm	mov	ecx, __o6\
			__asm	mov	edx, __o7\
			\
			__asm	movaps	xmm6,[edi+__i2]	/* t2 */\
			__asm	movaps	xmm7,[esi+__i2]	/* t3 */\
			__asm	subpd   xmm6,xmm2		/* ta = 2-a */\
			__asm	subpd   xmm7,xmm4		/* tb = 3-b */\
			__asm	addpd   xmm2,xmm2		/*      2*a */\
			__asm	addpd   xmm4,xmm4		/*      2*b */\
			__asm	addpd   xmm2,xmm6		/* t2 = 2+a */\
			__asm	addpd   xmm4,xmm7		/* t3 = 3+b */\
			\
			__asm	movaps	[ebx     ],xmm6	/* a[j1+p5] = ta */\
			__asm	movaps	[ebx+0x10],xmm7	/* a[j2+p5] = tb */\
			__asm	movaps	[eax     ],xmm2	/* a[j1+p4] = t2 */\
			__asm	movaps	[eax+0x10],xmm4	/* a[j2+p4] = t3 */\
			\
			/* Can't do 2nd set side-by-side with 1st, since only 2 registers [xmm6,7] available for temp storage: */\
			__asm	movaps	xmm6,[edi+__i6]	/* t6 */\
			__asm	movaps	xmm7,[esi+__i6]	/* t7 */\
			__asm	subpd   xmm6,xmm3		/* t6 = 6-f */\
			__asm	subpd   xmm7,xmm5		/* tf = 7-e */\
			__asm	addpd   xmm3,xmm3		/*      2*f */\
			__asm	addpd   xmm5,xmm5		/*      2*e */\
			__asm	addpd   xmm3,xmm6		/* te = 6+f */\
			__asm	addpd   xmm5,xmm7		/* t7 = 7+e */\
			\
			__asm	movaps	[ecx     ],xmm6	/* a[j1+p6] = t6 */\
			__asm	movaps	[edx+0x10],xmm7	/* a[j2+p7] = tf */\
			__asm	movaps	[edx     ],xmm3	/* a[j1+p7] = te */\
			__asm	movaps	[ecx+0x10],xmm5	/* a[j2+p6] = t7 */\
			\
			__asm	mov	eax, __o0\
			__asm	mov	ebx, __o1\
			__asm	mov	ecx, __o2\
			__asm	mov	edx, __o3\
			\
			__asm	movaps	xmm6,[edi     ]	/* t0 */			__asm	movaps	xmm4,[edi+__i4]	/* t4 */\
			__asm	movaps	xmm7,[esi     ]	/* t1 */			__asm	movaps	xmm5,[esi+__i4]	/* t5 */\
			__asm	movaps	xmm2,[edi+__i1]	/* t8 */\
			__asm	movaps	xmm3,[esi+__i1]	/* t9 */\
			__asm	subpd   xmm6,xmm2		/* t8 = 0-8 */		__asm	subpd   xmm4,xmm1		/* t4 = 4-d */\
			__asm	subpd   xmm7,xmm3		/* t9 = 1-9 */		__asm	subpd   xmm5,xmm0		/* td = 5-c */\
			__asm	addpd   xmm2,xmm2		/*      2*8 */		__asm	addpd   xmm1,xmm1		/*      2*d */\
			__asm	addpd   xmm3,xmm3		/*      2*9 */		__asm	addpd   xmm0,xmm0		/*      2*c */\
			__asm	addpd   xmm2,xmm6		/* t0 = 0+8 */		__asm	addpd   xmm1,xmm4		/* tc = 4+d */\
			__asm	addpd   xmm3,xmm7		/* t1 = 1+9 */		__asm	addpd   xmm0,xmm5		/* t5 = 5+c */\
			\
			__asm	movaps	[ebx     ],xmm6	/* t8 */			__asm	movaps	[ecx     ],xmm4	/* t4 */\
			__asm	movaps	[ebx+0x10],xmm7	/* t9 */			__asm	movaps	[edx+0x10],xmm5	/* td */\
			__asm	movaps	[eax     ],xmm2	/* t0 */			__asm	movaps	[edx     ],xmm1	/* tc */\
			__asm	movaps	[eax+0x10],xmm3	/* t1 */			__asm	movaps	[ecx+0x10],xmm0	/* t5 */\
		}
												/* Totals: 140 load/store [63 movaps], 54 add/subpd, 4 mulpd */

		/* NOTE: Must use out-array for temp storage here, otherwise input-index permutations hose outputs */
		#define SSE2_RADIX8_DIT_0TWIDDLE(__in0,__in1,__in2,__in3,__in4,__in5,__in6,__in7, __out, __isrt2)\
		{\
			__asm	mov	eax, __in4\
			__asm	mov	ebx, __in5\
			__asm	mov	ecx, __in6\
			__asm	mov	edx, __in7\
			__asm	mov	esi, __out\
			\
			/*** 2nd of the 2 length-4 subtransforms gets done first, due to e.g. t1-+t9 combos in final step: ***/\
			/*\
			t4r = a[j1+p4];		t4i = a[j2+p4];\
			t5r = a[j1+p5];		t5i = a[j2+p5];\
			~t5r = t4r-t5r;		~t5i = t4i-t5i;\
			~t4r = t4r+t5r;		~t4i = t4i+t5i;\
			*/\
			__asm	movaps	xmm0,[eax     ]	/* xmm0 <- a[j1+p4] = t4r*/\
			__asm	movaps	xmm1,[eax+0x10]	/* xmm1 <- a[j2+p4] = t4i*/\
			__asm	movaps	xmm2,xmm0		/* xmm4 <- copy of t4r   */\
			__asm	movaps	xmm3,xmm1		/* xmm5 <- copy of t4i   */\
			__asm	addpd	xmm2,[ebx     ]	/* xmm2 <- t4r */\
			__asm	addpd	xmm3,[ebx+0x10]	/* xmm3 <- t4i */\
			__asm	subpd	xmm0,[ebx     ]	/* xmm0 <- t5r */\
			__asm	subpd	xmm1,[ebx+0x10]	/* xmm1 <- t5i */\
			/*\
			t6r = a[j1+p6];		t6i = a[j2+p6];\
			rt  = a[j1+p7];		it  = a[j2+p7];\
			~t7r = t6r-t7r;		~t7i = t6i-t7i\
			~t6r = t6r+t7r;		~t6i = t6i+t7i\
			*/\
			__asm	movaps	xmm4,[ecx     ]	/* xmm4 <- a[j1+p6] = t6r*/\
			__asm	movaps	xmm5,[ecx+0x10]	/* xmm5 <- a[j2+p6] = t6i*/\
			__asm	movaps	xmm6,xmm4		/* xmm6 <- copy of t6r   */\
			__asm	movaps	xmm7,xmm5		/* xmm7 <- copy of t6i   */\
			__asm	addpd	xmm6,[edx     ]	/* xmm6 <- t6r */\
			__asm	addpd	xmm7,[edx+0x10]	/* xmm7 <- t6i */\
			__asm	subpd	xmm4,[edx     ]	/* xmm4 <- t7r */\
			__asm	subpd	xmm5,[edx+0x10]	/* xmm5 <- t7i */\
			/* Copy t6r,i into main-array slot add6 */\
			__asm	movaps	[esi+0xc0],xmm6\
			__asm	movaps	[esi+0xd0],xmm7\
			/* Copy t7r,i into main-array slot add7 */\
			__asm	movaps	[esi+0xe0],xmm4\
			__asm	movaps	[esi+0xf0],xmm5\
			\
			/** GPRs: ***** SSE Regs: ***** Temp array: ***\
			* eax, add4    xmm0 <- t5r    add0 <- unused  *\
			* ebx, add5    xmm1 <- t5i    add1 <- unused  *\
			* ecx, add6    xmm2 <- t4r    add2 <- unused  *\
			* edx, add7    xmm3 <- t4i    add3 <- unused  *\
			* esi,   p4    xmm4 <- t7r    add4 <- unused  *\
			* edi, ----    xmm5 <- t7i    add5 <- unused  *\
			*              xmm6 <- t6r    add6 <- t6r,i   *\
			*              xmm7 <- t6i    add7 <- t7r,i   *\
			**********************************************/\
			\
			/*\
			~t6r = t4r - t6r;	~t4r = t4r + t6r;	copies of t6r in add6     , xmm6\
			~t6i = t4i - t6i;	~t4i = t4i + t6i;	copies of t6i in add6+0x10, xmm7\
			\
			~t7r = t5r - t7i;	~t5r = t5r + t7i;	copies of t7r in add7     , xmm4\
			~t7i = t5i + t7r;	~t5i = t5i - t7r;	copies of t7i in add7+0x10, xmm5\
			*/\
			/* Move outputs t5r,i into a[j1,j2+p5], first doing the addsub and mul by ISRT2: */\
			/* Move outputs t7r,i into a[j1,j2+p7], first doing the addsub and mul by ISRT2: */\
			\
			__asm	addpd	xmm6,xmm2		/* xmm6 <- ~t4r */\
			__asm	addpd	xmm7,xmm3		/* xmm7 <- ~t4i */\
			__asm	subpd	xmm2,[esi+0xc0]	/* xmm2 <- ~t6r */\
			__asm	subpd	xmm3,[esi+0xd0]	/* xmm3 <- ~t6i */\
			/* Move t6r,i into a[j1,j2+p6] */\
			__asm	movaps	[esi+0xc0],xmm2	/* add6r <- ~t6r */\
			__asm	movaps	[esi+0xd0],xmm3	/* add6i <- ~t6i */\
			\
			__asm	movaps	xmm2,xmm4	/* xmm2 <- copy of t7r */\
			__asm	movaps	xmm3,xmm5	/* xmm3 <- copy of t7i */\
			__asm	addpd	xmm5,xmm0	/* xmm5 <-~t5r */\
			__asm	subpd	xmm0,xmm3	/* xmm0 <-~t7r */\
			__asm	addpd	xmm4,xmm1	/* xmm4 <-~t7i */\
			__asm	subpd	xmm1,xmm2	/* xmm1 <-~t5i */\
			\
			__asm	movaps	xmm2,xmm5	/* xmm2 <- copy of~t5r */\
			__asm	movaps	xmm3,xmm1	/* xmm3 <- copy of~t5i */\
			__asm	addpd	xmm5,xmm1	/* xmm5 <-~(t5r+t5i), xmm1 FREE */\
			__asm	mov	edi, __isrt2\
			__asm	movaps	xmm1,[edi]	/* xmm1 <- ISRT2 */\
			__asm	subpd	xmm2,xmm3	/* xmm2 <-~(t5r-t5i), xmm3 FREE */\
			__asm	mulpd	xmm5,xmm1	/* xmm5 <- (t5r+t5i)*ISRT2 */\
			__asm	mulpd	xmm2,xmm1	/* xmm2 <- (t5r-t5i)*ISRT2 */\
			__asm	movaps	xmm3,xmm0	/* xmm3 <- copy of~t7r */\
			\
			__asm	movaps	[esi+0xa0],xmm5	/* add5r<- (t5r+t5i)*ISRT2, xmm5 FREE */\
			__asm	movaps	xmm5,xmm4	/* xmm5 <- copy of~t7i */\
			__asm	addpd	xmm0,xmm4	/* xmm0 <-~(t7r+t7i) */\
			__asm	movaps	[esi+0xb0],xmm2	/* add5i<- (t5r-t5i)*ISRT2 */\
			__asm	subpd	xmm3,xmm5	/* xmm3 <-~(t7r-t7i) */\
			__asm	mulpd	xmm0,xmm1	/* xmm0 <- (t7r+t7i)*ISRT2 */\
			__asm	mulpd	xmm3,xmm1	/* xmm3 <- (t7r-t7i)*ISRT2 */\
			__asm	movaps	[esi+0xe0],xmm0	/* add7r<- (t7r+t7i)*ISRT2 */\
			__asm	movaps	[esi+0xf0],xmm3	/* add7i<- (t7r-t7i)*ISRT2 */\
			\
			/** GPRs: ***** SSE Regs: ******** Temp array: *************\
			* eax, add4    xmm0 <- unused    add0 <- unused            *\
			* ebx, add5    xmm1 <- unused    add1 <- unused            *\
			* ecx, add6    xmm2 <- unused    add2 <- unused            *\
			* edx, add7    xmm3 <- unused    add3 <- unused            *\
			* esi,   p4    xmm4 <- unused    add4 <- unused            *\
			* edi,isrt2    xmm5 <- unused    add5 <- (t5r+-t5i)*ISRT2  *\
			*              xmm6 <- t4r       add6 <-  t6r,i            *\
			*              xmm7 <- t4i       add7 <- (t7r+-t7i)*ISRT2  *\
			***********************************************************/\
			\
			/**************** 1st of the 2 length-4 subtransforms... **************/\
			/*\
			t0r = a[j1   ];		t0i = a[j2   ];\
			rt  = a[j1+p1];		it  = a[j2+p1];\
			t1r = t0r-rt;		t1i = t0i-it;\
			t0r = t0r+rt;		t0i = t0i+it;\
			*/\
			__asm	mov	eax, __in0\
			__asm	mov	ebx, __in1\
			__asm	mov	ecx, __in2\
			__asm	mov	edx, __in3\
			\
			__asm	movaps	xmm0,[eax     ]	/* xmm0 <- a[j1   ] = t0r*/\
			__asm	movaps	xmm1,[eax+0x10]	/* xmm1 <- a[j2   ] = t0i*/\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- copy of t0r   */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- copy of t0i   */\
			__asm	addpd	xmm2,[ebx     ]	/*~xmm2 <- t0r */\
			__asm	addpd	xmm3,[ebx+0x10]	/*~xmm3 <- t0i */\
			__asm	subpd	xmm0,[ebx     ]	/*~xmm0 <- t1r */\
			__asm	subpd	xmm1,[ebx+0x10]	/*~xmm1 <- t1i */\
			\
			/* Move   t4r,i into temp[j1+p0] in anticipation of final outputs [t0+t4]r,i which will go there: */\
			__asm	movaps	[esi     ],xmm6\
			__asm	movaps	[esi+0x10],xmm7	/* add0 <-  t4r,i */\
			__asm	addpd	xmm6,xmm6		/* xmm6 <- 2*t4r */\
			__asm	addpd	xmm7,xmm7		/* xmm7 <- 2*t4i */\
			/* Move 2*t4r,i into temp[j1+p4] in anticipation of final outputs [t0-t4]r,i which will go there: */\
			__asm	movaps	[esi+0x80],xmm6	/* add4 <- 2*t4r,i */\
			__asm	movaps	[esi+0x90],xmm7	/* xmm4-7 FREE */\
			/*\
			t2r = a[j1+p2];		t2i = a[j2+p2];\
			rt  = a[j1+p3];		it  = a[j2+p3];\
			t3r = t2r-rt;		t3i = t2i-it;\
			t2r = t2r+rt;		t2i = t2i+it;\
			*/\
			__asm	movaps	xmm4,[ecx     ]	/* xmm4 <- a[j1+p2] = t2r*/\
			__asm	movaps	xmm5,[ecx+0x10]	/* xmm5 <- a[j2+p2] = t2i*/\
			__asm	movaps	xmm6,xmm4		/* xmm6 <- copy of t2r   */\
			__asm	movaps	xmm7,xmm5		/* xmm7 <- copy of t2i   */\
			__asm	addpd	xmm6,[edx     ]	/*~xmm6 <- t2r */\
			__asm	addpd	xmm7,[edx+0x10]	/*~xmm7 <- t2i */\
			__asm	subpd	xmm4,[edx     ]	/*~xmm4 <- t3r */\
			__asm	subpd	xmm5,[edx+0x10]	/*~xmm5 <- t3i */\
			/* Copy t5,6 into temp[j1+p2] */\
			__asm	movaps	[esi+0x40],xmm6\
			__asm	movaps	[esi+0x50],xmm7	/* add2 <- t2r,i*/\
			\
			/** GPRs: ***** SSE Regs: ******** Main array: *************\
			* eax, add0    xmm0 <- t1r       add0 <-   t4r,i           *\
			* ebx, add1    xmm1 <- t1i       add1 <- unused            *\
			* ecx, add2    xmm2 <- t0r       add2 <-   t2r,i           *\
			* edx, add3    xmm3 <- t0i       add3 <- unused            *\
			* esi,   p4    xmm4 <- t3r       add4 <- 2*t4r,i           *\
			* edi, add4    xmm5 <- t3i       add5 <- (t5r+-t5i)*ISRT2  *\
			*              xmm6 <- t2r       add6 <-   t6r,i           *\
			*              xmm7 <- t2i       add7 <- (t7r+-t7i)*ISRT2  *\
			***********************************************************/\
			\
			/*\
			~t2r = t0r - t2r;	~t0r = t0r + t2r;	copies of t2r in add2     , xmm6\
			~t2i = t0i - t2i;	~t0i = t0i + t2i;	copies of t2i in add2+0x10, xmm7\
			\
			~t3r = t1r - t3i;	~t1r = t1r + t3i;\
			~t3i = t1i + t3r;	~t1i = t1i - t3r;\
			*/\
			__asm	addpd	xmm6,xmm2		/* xmm6 <- ~t0r*/\
			__asm	addpd	xmm7,xmm3		/* xmm7 <- ~t0i*/\
			__asm	subpd	xmm2,[esi+0x40]	/* xmm2 <- ~t2r*/\
			__asm	subpd	xmm3,[esi+0x50]	/* xmm3 <- ~t2i*/\
			\
			/* Compute and dump first 2 outputs now, in order to free up 2 registers: */\
			__asm	addpd	xmm6,[esi     ]	/* t1+t4r*/\
			__asm	addpd	xmm7,[esi+0x10]	/* t2+t4i*/\
			__asm	movaps	[esi     ],xmm6	/* a[j1   ], DONE. */\
			__asm	movaps	[esi+0x10],xmm7	/* a[j2   ], DONE. */\
			\
			__asm	subpd	xmm6,[esi+0x80]	/* t1-t4r = [t1+t4r] - 2*t4r */\
			__asm	subpd	xmm7,[esi+0x90]	/* t2-t4i = [t2+t4i] - 2*t4i */\
			__asm	movaps	[esi+0x80],xmm6	/* a[j1+p4], DONE. */\
			__asm	movaps	[esi+0x90],xmm7	/* a[j2+p4], DONE. */\
			\
			__asm	movaps	xmm6,xmm4		/* xmm6 <- copy of t3r*/\
			__asm	movaps	xmm7,xmm5		/* xmm7 <- copy of t3i*/\
			__asm	addpd	xmm5,xmm0	/* xmm5 <- ~t1r*/\
			__asm	subpd	xmm0,xmm7	/* xmm0 <- ~t3r*/\
			__asm	addpd	xmm4,xmm1	/* xmm4 <- ~t3i*/\
			__asm	subpd	xmm1,xmm6	/* xmm1 <- ~t1i*/\
			\
			/** GPRs: ***** SSE Regs: ******** Main array: *************\
			* eax, ----    xmm0 <- t3r       add0 <- DONE              *\
			* ebx, add1    xmm1 <- t1i       add1 <- unused            *\
			* ecx, add2    xmm2 <- t2r       add2 <- unused            *\
			* edx, add3    xmm3 <- t2i       add3 <- unused            *\
			* esi,   p4    xmm4 <- t3i       add4 <- DONE              *\
			* edi, add4    xmm5 <- t1r       add5 <- (t5r+-t5i)*ISRT2  *\
			*              xmm6 <- unused    add6 <-  t6r,14           *\
			*              xmm7 <- unused    add7 <- (t7r+-t7i)*ISRT2  *\
			***********************************************************/\
			\
			/* Now combine the two half-transforms & store outputs back into original array slots: */\
			/*\
			a[j1   ] = t0r+t4r;			a[j2   ] = t0i+t4i;	already done\
			a[j1+p4] = t0r-t4r;			a[j2+p4] = t0i-t4i;	already done\
			*/\
			\
			/*\
			a[j1+p2] = t2r+t6i;			a[j2+p2] = t2i-t6r;\
			a[j1+p6] = t2r-t6i;			a[j2+p6] = t2i+t6r;\
			*/\
			__asm	movaps	xmm6,xmm2		/* xmm6 <- copy of t2r*/\
			__asm	movaps	xmm7,xmm3		/* xmm7 <- copy of t2i*/\
			__asm	addpd	xmm2,[esi+0xd0]	/*  rt */\
			__asm	subpd	xmm3,[esi+0xc0]	/*  it */\
			__asm	subpd	xmm6,[esi+0xd0]	/* ~t2r */\
			__asm	addpd	xmm7,[esi+0xc0]	/* ~t2i */\
			__asm	movaps	[esi+0xc0],xmm2	/* a[j1+p6] */\
			__asm	movaps	[esi+0xd0],xmm3	/* a[j2+p6] */\
			__asm	movaps	[esi+0x40],xmm6	/* a[j1+p2] */\
			__asm	movaps	[esi+0x50],xmm7	/* a[j2+p2] xmm2,3,6,7 FREE ***/\
			\
			/* Move these 2 blocks together to enable output-index-swapping 1<>7, 3<>5 required to match scalar RADIX_08_DIT output oredring: */\
			/*\
			rt = (t5r+t5i)*ISRT2;	it = (t5r-t5i)*ISRT2	[rt,it] in add5			rt = (t7r-t7i)*ISRT2;		it = (t7r+t7i)*ISRT2;	[it,rt] in add7; NOTE reversed order!\
			a[j1+p1] = t1r+rt;		a[j1+p1] = t1i-it;								a[j1+p3] = t3r-rt;			a[j2+p3] = t3i-it;\
			a[j1+p5] = t1r-rt;		a[j1+p5] = t1i+it;								a[j1+p7] = t3r+rt;			a[j2+p7] = t3i+it;\
			*/\
			__asm	movaps	xmm2,xmm5		/* xmm6 <- copy of t1r*/				__asm	movaps	xmm6,xmm0		/* xmm6 <- copy of t3r*/\
			__asm	movaps	xmm3,xmm1		/* xmm7 <- copy of t1i*/				__asm	movaps	xmm7,xmm4		/* xmm7 <- copy of t3i*/\
			__asm	addpd	xmm5,[esi+0xa0]	/* ~t5r */								__asm	subpd	xmm0,[esi+0xf0]	/* ~t7r */\
			__asm	subpd	xmm1,[esi+0xb0]	/* ~t5i */								__asm	subpd	xmm4,[esi+0xe0]	/* ~t7i */\
			__asm	subpd	xmm2,[esi+0xa0]	/* ~t1r */								__asm	addpd	xmm6,[esi+0xf0]	/* ~t3r */\
			__asm	addpd	xmm3,[esi+0xb0]	/* ~t1i */								__asm	addpd	xmm7,[esi+0xe0]	/* ~t3i */\
			__asm	movaps	[esi+0xe0],xmm5	/* a[j1+p7] */							__asm	movaps	[esi+0xa0],xmm0	/* a[j1+p5] */\
			__asm	movaps	[esi+0xf0],xmm1	/* a[j2+p7] */							__asm	movaps	[esi+0xb0],xmm4	/* a[j2+p5] */\
			__asm	movaps	[esi+0x60],xmm2	/* a[j1+p3] */							__asm	movaps	[esi+0x20],xmm6	/* a[j1+p1] */\
			__asm	movaps	[esi+0x70],xmm3	/* a[j2+p3] */							__asm	movaps	[esi+0x30],xmm7	/* a[j2+p1] */\
		}
												/* Totals: 140 load/store [63 movaps], 54 add/subpd, 4 mulpd */

		#define SSE2_RADIX8_DIF_TWIDDLE(__add0,__p1,__p2,__p4,__p6,__r0,__c0)\
		{\
			/*******************************************************/\
			/* Inline contents of SSE2_RADIX4_DIF_4TWIDDLE_B(r0,c0): Do 1st of 2 radix-4 subtransforms and store outputs to temps starting at r0: */\
			/*******************************************************/\
			/* Can't get these via simple load-one-and-shift-as-needed due to array padding scheme */\
			__asm	mov	eax, __add0	/* &a[j1+p0] */\
			__asm	mov	ebx, __p2	/* &a[j1+p2] */\
			__asm	mov	ecx, __p4	/* &a[j1+p4] */\
			__asm	mov	edx, __p6	/* &a[j1+p6] */\
			__asm	shl	ebx,  3		/* 8-bytes for array-of-doubles */\
			__asm	shl	ecx,  3\
			__asm	shl	edx,  3\
			__asm	add	ebx, eax\
			__asm	add	ecx, eax\
			__asm	add	edx, eax\
			\
			/* Do the p0,p4 combo: */\
			__asm	mov	esi, __c0\
			\
			__asm	movaps	xmm0,[eax     ]	/* a[j1+p0] */				__asm	movaps	xmm4,[ecx     ]	/* a[j1+p4] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[j2+p0] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[j2+p4] */\
			__asm	movaps	xmm6,[esi     ]	/* c0 */\
			__asm	movaps	xmm7,[esi+0x10]	/* s0 */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[j1   ] */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[j2   ] */\
			\
			__asm	mulpd	xmm0,xmm6		/* a[j1   ]*c0 */\
			__asm	mulpd	xmm1,xmm6		/* a[j2   ]*c0 */\
			__asm	mulpd	xmm2,xmm7		/* a[j1   ]*s0 */\
			__asm	mulpd	xmm3,xmm7		/* a[j2   ]*s0	xmm6,7 free */\
																		__asm	add	esi, 0x20	/* c4 */\
																		__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[j1+p4] */\
			__asm	addpd	xmm1,xmm2		/* xmm1 <- t2 */			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[j2+p4] */\
																		__asm	mulpd	xmm4,[esi     ] /* a[j1+p4]*c4 */\
			__asm	subpd	xmm0,xmm3		/* xmm0 <- t1 */			__asm	mulpd	xmm5,[esi     ]	/* a[j2+p4]*c4 */\
																		__asm	mulpd	xmm6,[esi+0x10]	/* a[j1+p4]*s4 */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy t1 */		__asm	mulpd	xmm7,[esi+0x10]	/* a[j2+p4]*s4 */\
																		__asm	addpd	xmm5,xmm6	    /* xmm5 <- it */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy t2 */		__asm	subpd	xmm4,xmm7		/* xmm4 <- rt    xmm6,7 free */\
			\
			__asm	addpd	xmm0,xmm4		/* ~t1 <- t1 +rt */\
			__asm	addpd	xmm1,xmm5		/* ~t2 <- t2 +it */\
			__asm	subpd	xmm2,xmm4		/* ~t3 <- t1 -rt */\
			__asm	subpd	xmm3,xmm5		/* ~t4 <- t2 -it	xmm4,5 free */\
			\
			/* Do the p2,6 combo: */\
			__asm	add	esi, 0x40	/* c6 */\
			\
			__asm	movaps	xmm4,[edx     ]	/* a[jt+p6] */\
			__asm	movaps	xmm5,[edx+0x10]	/* a[jp+p6] */\
			__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[jt+p6] */\
			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[jp+p6] */\
			\
			__asm	mulpd	xmm4,[esi     ]	/* a[jt+p6]*c6 */\
			__asm	mulpd	xmm5,[esi     ]	/* a[jp+p6]*c6 */\
			__asm	mulpd	xmm6,[esi+0x10]	/* a[jt+p6]*s6 */\
			__asm	mulpd	xmm7,[esi+0x10]	/* a[jp+p6]*s6 */\
			\
			__asm	mov	esi, __r0\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- it */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- rt */\
			__asm	movaps	[esi+0x010],xmm5	/* store it in a[j2+p0] */\
			__asm	movaps	[esi      ],xmm4	/* store rt in a[j1+p0] */\
			\
			__asm	mov	esi, __c0\
			__asm	add	esi, 0x40	/* c2 */\
			__asm	movaps	xmm4,[ebx     ]	/* a[jt+p2] */\
			__asm	movaps	xmm5,[ebx+0x10]	/* a[jp+p2] */\
			__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[jt+p2] */\
			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[jp+p2] */\
			\
			__asm	mulpd	xmm4,[esi     ]	/* a[jt+p2]*c2 */\
			__asm	mulpd	xmm5,[esi     ]	/* a[jp+p2]*c2 */\
			__asm	mulpd	xmm6,[esi+0x10]	/* a[jt+p2]*s2 */\
			__asm	mulpd	xmm7,[esi+0x10]	/* a[jp+p2]*s2 */\
			__asm	mov	esi, __r0\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- t5 */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- t4 	xmm6,7 free */\
			__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy t5 */\
			__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy t4 */\
			\
			__asm	subpd	xmm4,[esi      ]	/* ~t6 <- t4 -rt */\
			__asm	subpd	xmm5,[esi+0x010]	/* ~t7 <- t5 -it */\
			__asm	addpd	xmm6,[esi      ]	/* ~t4 <- t4 +rt */\
			__asm	addpd	xmm7,[esi+0x010]	/* ~t5 <- t5 +it */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			/*\
			~t4 =t0 -t4;		~t0 =t0 +t4;							~t6 =t2 +t7;		~t2 =t2 -t7;\
			~t5 =t1 -t5;		~t1 =t1 +t5;							~t7 =t3 -t6;		~t3 =t3 +t6;\
			*/\
			__asm	subpd	xmm0,xmm6	/*~t4 */						__asm	subpd	xmm2,xmm5	/*~t2 */\
			__asm	subpd	xmm1,xmm7	/*~t5 */						__asm	subpd	xmm3,xmm4	/*~t7 */\
			__asm	movaps	[esi+0x040],xmm0	/* a[jt+p4] <- ~t4 */	__asm	movaps	[esi+0x020],xmm2	/* a[jt+p2] <- ~t2 */\
			__asm	movaps	[esi+0x050],xmm1	/* a[jp+p4] <- ~t5 */	__asm	movaps	[esi+0x070],xmm3	/* a[jp+p6] <- ~t7 */\
			__asm	addpd	xmm6,xmm6	/* 2*t4 */						__asm	addpd	xmm5,xmm5	/* 2*t7 */\
			__asm	addpd	xmm7,xmm7	/* 2*t5 */						__asm	addpd	xmm4,xmm4	/* 2*t6 */\
			__asm	addpd	xmm6,xmm0	/*~t0 */						__asm	addpd	xmm5,xmm2	/*~t6 */\
			__asm	addpd	xmm7,xmm1	/*~t1 */						__asm	addpd	xmm4,xmm3	/*~t3 */\
			__asm	movaps	[esi      ],xmm6	/* a[jt+p0] <- ~t0 */	__asm	movaps	[esi+0x060],xmm5	/* a[jt+p6] <- ~t6 */\
			__asm	movaps	[esi+0x010],xmm7	/* a[jp+p0] <- ~t1 */	__asm	movaps	[esi+0x030],xmm4	/* a[jp+p2] <- ~t3 */\
			\
			/*******************************************************/\
			/* Inline of SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r8,c1): Do 2nd of 2 radix-4 subtransforms, but now keep outputs in registers: */\
			/*******************************************************/\
			__asm	mov	edi, __p1\
			__asm	shl	edi,  3		/* 8-bytes for array-of-doubles */\
			__asm	add	eax, edi	/* &a[j1+p1] */\
			__asm	add	ebx, edi	/* &a[j1+p3] */\
			__asm	add	ecx, edi	/* &a[j1+p5] */\
			__asm	add	edx, edi	/* &a[j1+p7] */\
			\
			/* Do the p1,p5 combo: */\
			__asm	mov	esi, __c0\
			__asm	add	esi, 0x80	/* c1 */\
			__asm	movaps	xmm0,[eax     ]	/* a[jt+p1] */				__asm	movaps	xmm4,[ecx     ]	/* a[jt+p5] */\
			__asm	movaps	xmm1,[eax+0x10]	/* a[jp+p1] */				__asm	movaps	xmm5,[ecx+0x10]	/* a[jp+p5] */\
			__asm	movaps	xmm6,[esi     ]	/* c1 */\
			__asm	movaps	xmm7,[esi+0x10]	/* s1 */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy a[jt+p1] */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy a[jp+p1] */\
			\
			__asm	mulpd	xmm0,xmm6		/* a[jt+p1]*c1 */\
			__asm	mulpd	xmm1,xmm6		/* a[jp+p1]*c1 */\
			__asm	mulpd	xmm2,xmm7		/* a[jt+p1]*s1 */\
			__asm	mulpd	xmm3,xmm7		/* a[jp+p1]*s1	xmm6,7 free */\
																		__asm	add	esi, 0x20	/* c5 */\
																		__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[jt+p5] */\
			__asm	addpd	xmm1,xmm2		/* xmm1 <- t9 */			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[jp+p5] */\
																		__asm	mulpd	xmm4,[esi     ] /* a[jt+p5]*c5 */\
			__asm	subpd	xmm0,xmm3		/* xmm0 <- t8 */			__asm	mulpd	xmm5,[esi     ]	/* a[jp+p5]*c5 */\
																		__asm	mulpd	xmm6,[esi+0x10]	/* a[jt+p5]*s5 */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- cpy t8 */		__asm	mulpd	xmm7,[esi+0x10]	/* a[jp+p5]*s5 */\
																		__asm	addpd	xmm5,xmm6	    /* xmm5 <- it */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- cpy t9 */		__asm	subpd	xmm4,xmm7		/* xmm4 <- rt    xmm6,7 free */\
			\
			__asm	addpd	xmm0,xmm4	/* ~t8 <- t8 +rt */\
			__asm	addpd	xmm1,xmm5	/* ~t9 <- t9 +it */\
			__asm	subpd	xmm2,xmm4	/* ~ta <- t8 -rt */\
			__asm	subpd	xmm3,xmm5	/* ~tb <- t9 -it	xmm4,5 free */\
			\
			/* Do the p3,7 combo: */\
			__asm	add	esi, 0x40	/* c7 */\
			\
			__asm	movaps	xmm4,[edx     ]	/* a[jt+p7] */\
			__asm	movaps	xmm5,[edx+0x10]	/* a[jp+p7] */\
			__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[jt+p7] */\
			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[jp+p7] */\
			\
			__asm	mulpd	xmm4,[esi     ]	/* a[jt+p7]*c7 */\
			__asm	mulpd	xmm5,[esi     ]	/* a[jp+p7]*c7 */\
			__asm	mulpd	xmm6,[esi+0x10]	/* a[jt+p7]*s7 */\
			__asm	mulpd	xmm7,[esi+0x10]	/* a[jp+p7]*s7 */\
			\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- it */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- rt */\
			__asm	movaps	[eax+0x010],xmm5	/* store it in a[j2+p1] */\
			__asm	movaps	[eax      ],xmm4	/* store rt in a[j1+p1] */\
			\
			__asm	sub	esi, 0x20	/* c3 */\
			__asm	movaps	xmm4,[ebx     ]	/* a[jt+p3] */\
			__asm	movaps	xmm5,[ebx+0x10]	/* a[jp+p3] */\
			__asm	movaps	xmm6,xmm4		/* xmm6 <- cpy a[jt+p3] */\
			__asm	movaps	xmm7,xmm5		/* xmm7 <- cpy a[jp+p3] */\
			\
			__asm	mulpd	xmm4,[esi     ]	/* a[jt+p3]*c3 */\
			__asm	mulpd	xmm5,[esi     ]	/* a[jp+p3]*c3 */\
			__asm	mulpd	xmm6,[esi+0x10]	/* a[jt+p3]*s3 */\
			__asm	mulpd	xmm7,[esi+0x10]	/* a[jp+p3]*s3 */\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- td */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- tc 	xmm6,7 free */\
			__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy td */\
			__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy tc */\
			\
			__asm	subpd	xmm4,[eax      ]	/* ~te <- tc -rt */\
			__asm	subpd	xmm5,[eax+0x010]	/* ~tf <- td -it */\
			__asm	addpd	xmm6,[eax      ]	/* ~tc <- tc +rt */\
			__asm	addpd	xmm7,[eax+0x010]	/* ~td <- td +it */\
			\
			/* Finish radix-4 butterfly and store just the 1st of the 8 outputs into temporary-array slots: */\
			/*\
			~tc =t8 -tc;		~t8 =t8 +tc;							~te =ta +tf;		~ta =ta -tf;\
			~td =t9 -td;		~t9 =t9 +td;							~tf =tb -te;		~tb =tb +te;\
			*/\
			__asm	subpd	xmm0,xmm6	/*~tc */\
			__asm	subpd	xmm1,xmm7	/*~td */\
		/*	__asm	movaps	[esi+0x040],xmm0	// a[jt+p5] <-~tc */	__asm	subpd	xmm2,xmm5	/*~ta */\
		/*	__asm	movaps	[esi+0x050],xmm1	// a[jp+p5] <-~td */	__asm	subpd	xmm3,xmm4	/*~tf */\
			__asm	addpd	xmm6,xmm6	/* 2*tc */						__asm	addpd	xmm5,xmm5	/* 2*tf */\
			__asm	addpd	xmm7,xmm7	/* 2*td */						__asm	addpd	xmm4,xmm4	/* 2*te */\
			__asm	addpd	xmm6,xmm0	/*~t8 */						__asm	addpd	xmm5,xmm2	/*~te */\
			__asm	addpd	xmm7,xmm1	/*~t9 */						__asm	addpd	xmm4,xmm3	/*~tb */\
			__asm	movaps	[eax      ],xmm6	/* a[jt+p1] <-~t8 */\
			__asm	movaps	[eax+0x010],xmm7	/* a[jp+p1] <-~t9; xmm6,7 FREE */\
			/*\
			ta =(ta-tb)*ISRT2;	tb =(ta+tb)*ISRT2;\
			te =(te-tf)*ISRT2;	tf =(te+tf)*ISRT2;\
			*/\
			__asm	movaps	xmm6,xmm2	/* xmm6 <- cpy ta */\
			__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy te */\
			__asm	subpd	xmm2,xmm4	/* ta-tb*/\
			__asm	subpd	xmm5,xmm3	/* te-tf*/\
			__asm	addpd	xmm4,xmm6	/* ta+tb*/\
			__asm	addpd	xmm3,xmm7	/* te+tf*/\
			\
			__asm	mov	esi, isrt2\
			__asm	movaps	xmm6,[esi]	/* ISRT2 */\
			__asm	mulpd	xmm2,xmm6	/* ~ta = (ta-tb)*ISRT2 */\
			__asm	mulpd	xmm5,xmm6	/* ~te = (te-tf)*ISRT2 */\
			__asm	mulpd	xmm4,xmm6	/* ~tb = (ta+tb)*ISRT2 */\
			__asm	mulpd	xmm3,xmm6	/* ~tf = (te+tf)*ISRT2 */\
		/*	__asm	movaps	[esi+0x020],xmm2	// a[jt+p3] <-~ta */\
		/*	__asm	movaps	[esi+0x060],xmm5	// a[jp+p7] <-~te */\
		/*	__asm	movaps	[esi+0x030],xmm4	// a[jp+p3] <-~tb */\
		/*	__asm	movaps	[esi+0x070],xmm3	// a[jt+p7] <-~tf */\
			\
			/* Inline the contents of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_A(r0): Combine the 2 radix-4 subtransforms and write outputs back to main array: */\
		/***** t0,1,2,3,4,5,6,7 in r0+0x00,10,20,30,40,50,60,70 **************/\
		/***** t8,9,a,b,c,d,e,f in a[j1,j2 + p1],xmm2,4,0,1,5,3; xmm6,7 FREE */\
			__asm	mov	esi, __r0\
			\
			__asm	mov	eax, ecx	/* ecx = a[j1+p5] */\
			__asm	sub	eax, edi	/* eax = a[j1+p4] */\
			__asm	mov	ebx, edx	/* edx = a[j1+p7] */\
			__asm	sub	ebx, edi	/* ebx = a[j1+p6] */\
			\
			__asm	movaps	xmm6,[esi+0x20]	/* t2 */\
			__asm	movaps	xmm7,[esi+0x30]	/* t3 */\
			__asm	subpd	xmm6,xmm2		/* ta = 2-a */\
			__asm	subpd	xmm7,xmm4		/* tb = 3-b */\
			__asm	addpd	xmm2,xmm2		/*      2*a */\
			__asm	addpd	xmm4,xmm4		/*      2*b */\
			__asm	addpd	xmm2,xmm6		/* t2 = 2+a */\
			__asm	addpd	xmm4,xmm7		/* t3 = 3+b */\
			\
			__asm	movaps	[ecx     ],xmm6	/* a[j1+p5] = ta */\
			__asm	movaps	[ecx+0x10],xmm7	/* a[j2+p5] = tb */\
			__asm	movaps	[eax     ],xmm2	/* a[j1+p4] = t2 */\
			__asm	movaps	[eax+0x10],xmm4	/* a[j2+p4] = t3 */\
			\
			/* Can't do 2nd set side-by-side with 1st, since only 2 registers [xmm6,7] available for temp storage: */\
			__asm	movaps	xmm6,[esi+0x60]	/* t6 */\
			__asm	movaps	xmm7,[esi+0x70]	/* t7 */\
			__asm	subpd	xmm6,xmm3		/* t6 = 6-f */\
			__asm	subpd	xmm7,xmm5		/* tf = 7-e */\
			__asm	addpd	xmm3,xmm3		/*      2*f */\
			__asm	addpd	xmm5,xmm5		/*      2*e */\
			__asm	addpd	xmm3,xmm6		/* te = 6+f */\
			__asm	addpd	xmm5,xmm7		/* t7 = 7+e */\
			\
			__asm	movaps	[ebx     ],xmm6	/* a[j1+p6] = t6 */\
			__asm	movaps	[edx+0x10],xmm7	/* a[j2+p7] = tf */\
			__asm	movaps	[edx     ],xmm3	/* a[j1+p7] = te */\
			__asm	movaps	[ebx+0x10],xmm5	/* a[j2+p6] = t7 */\
			\
			__asm	mov	edi, __p4	/* Can't simply do p4 = (p1<<2) in-place, due to array-padding scheme */\
			__asm	shl	edi,  3		/* 8-bytes for array-of-doubles */\
			__asm	sub	eax, edi	/* &a[j1+p0] */\
			__asm	sub	ebx, edi	/* &a[j1+p2] */\
			__asm	sub	ecx, edi	/* &a[j1+p1] */\
			__asm	sub	edx, edi	/* &a[j1+p3] */\
			__asm	shr	edi,  2	/* p1 = (p4>>2) works, because any excess padding [i.e. present in p4, but not p1] gets shifted off */\
			\
			__asm	movaps	xmm6,[esi     ]	/* t0 */			__asm	movaps	xmm4,[esi+0x40]	/* t4 */\
			__asm	movaps	xmm7,[esi+0x10]	/* t1 */			__asm	movaps	xmm5,[esi+0x50]	/* t5 */\
			__asm	movaps	xmm2,[ecx     ]	/* t8 in a[j1+p1] */\
			__asm	movaps	xmm3,[ecx+0x10]	/* t9 in a[j2+p1] */\
			__asm	subpd	xmm6,xmm2		/* t8 = 0-8 */		__asm	subpd	xmm4,xmm1		/* t4 = 4-d */\
			__asm	subpd	xmm7,xmm3		/* t9 = 1-9 */		__asm	subpd	xmm5,xmm0		/* td = 5-c */\
			__asm	addpd	xmm2,xmm2		/*      2*8 */		__asm	addpd	xmm1,xmm1		/*      2*d */\
			__asm	addpd	xmm3,xmm3		/*      2*9 */		__asm	addpd	xmm0,xmm0		/*      2*c */\
			__asm	addpd	xmm2,xmm6		/* t0 = 0+8 */		__asm	addpd	xmm1,xmm4		/* tc = 4+d */\
			__asm	addpd	xmm3,xmm7		/* t1 = 1+9 */		__asm	addpd	xmm0,xmm5		/* t5 = 5+c */\
			\
			__asm	movaps	[ecx     ],xmm6	/* t8 */			__asm	movaps	[ebx     ],xmm4	/* t4 */\
			__asm	movaps	[ecx+0x10],xmm7	/* t9 */			__asm	movaps	[edx+0x10],xmm5	/* td */\
			__asm	movaps	[eax     ],xmm2	/* t0 */			__asm	movaps	[edx     ],xmm1	/* tc */\
			__asm	movaps	[eax+0x10],xmm3	/* t1 */			__asm	movaps	[ebx+0x10],xmm0	/* t5 */\
		}
												/* Totals: 185 load/store [93 movaps], 84 add/subpd, 36 mulpd */

		#define SSE2_RADIX8_DIT_TWIDDLE(__add0,__p4,__p5,__p6,__p7,__r0,__c4)\
		{\
			/* Can't get these via simple load-one-and-shift-as-needed due to array padding scheme */\
			__asm	mov	esi, __p4	/* esi will store ptr offset for p4 throughout */\
			__asm	shl	esi,  3		/* 8-bytes for array-of-doubles */\
			__asm	mov	eax, __add0	/* &a[j1+p0] */\
			__asm	mov	ebx, __p5\
			__asm	mov	ecx, __p6\
			__asm	mov	edx, __p7\
			__asm	shl	ebx,  3\
			__asm	shl	ecx,  3\
			__asm	shl	edx,  3\
			__asm	add	ebx, eax	/* &a[j1+p5] */\
			__asm	add	ecx, eax	/* &a[j1+p6] */\
			__asm	add	edx, eax	/* &a[j1+p7] */\
			__asm	add	eax, esi	/* &a[j1+p4] */\
			\
			/*** 2nd of the 2 length-4 subtransforms gets done first, due to e.g. t1-+t9 combos in final step: ***/\
			/*\
			t9 =a[j1+p4];	t10=a[j2+p4];\
			t11=a[j1+p5];	t12=a[j2+p5];\
			~t11=t9 -t11;	~t12=t10-t12;\
			~t9 =t9 +t11;	~t10=t10+t12;\
			*/\
			__asm	movaps	xmm0,[eax]		/* xmm0 <- a[j1+p4] = t9 */\
			__asm	movaps	xmm1,[eax+0x10]	/* xmm1 <- a[j2+p4] = t10*/\
			__asm	movaps	xmm2,xmm0		/* xmm4 <- copy of t9    */\
			__asm	movaps	xmm3,xmm1		/* xmm5 <- copy of t10   */\
			__asm	addpd	xmm2,[ebx]		/* xmm2 <- t9  */\
			__asm	addpd	xmm3,[ebx+0x10]	/* xmm3 <- t10 */\
			__asm	subpd	xmm0,[ebx]		/* xmm0 <- t11 */\
			__asm	subpd	xmm1,[ebx+0x10]	/* xmm1 <- t12 */\
			/*\
			t13=a[j1+p6];	t14=a[j2+p6];\
			rt =a[j1+p7];	it =a[j2+p7];\
			~t15=t13-t15	~t16=t14-t16\
			~t13=t13+t15	~t14=t14+t16\
			*/\
			__asm	movaps	xmm4,[ecx]		/* xmm4 <- a[j1+p6] = t13*/\
			__asm	movaps	xmm5,[ecx+0x10]	/* xmm5 <- a[j2+p6] = t14*/\
			__asm	movaps	xmm6,xmm4		/* xmm6 <- copy of t13   */\
			__asm	movaps	xmm7,xmm5		/* xmm7 <- copy of t14   */\
			__asm	addpd	xmm6,[edx]		/* xmm6 <- t13 */\
			__asm	addpd	xmm7,[edx+0x10]	/* xmm7 <- t14 */\
			__asm	subpd	xmm4,[edx]		/* xmm4 <- t15 */\
			__asm	subpd	xmm5,[edx+0x10]	/* xmm5 <- t16 */\
			/* Copy t13,14 into main-array slot add6 */\
			__asm	movaps	[ecx     ],xmm6\
			__asm	movaps	[ecx+0x10],xmm7\
			/* Copy t15,16 into main-array slot add7 */\
			__asm	movaps	[edx     ],xmm4\
			__asm	movaps	[edx+0x10],xmm5\
			\
			/** GPRs: ***** SSE Regs: ***** Main array: ***\
			* eax, add4    xmm0 <- t11    add0 <- unused  *\
			* ebx, add5    xmm1 <- t12    add1 <- unused  *\
			* ecx, add6    xmm2 <- t9     add2 <- unused  *\
			* edx, add7    xmm3 <- t10    add3 <- unused  *\
			* esi,   p4    xmm4 <- t15    add4 <- unused  *\
			* edi, ----    xmm5 <- t16    add5 <- unused  *\
			*              xmm6 <- t13    add6 <- t13,14  *\
			*              xmm7 <- t14    add7 <- t15,16  *\
			**********************************************/\
			\
			/*\
			rt =t13;	t13=t9 -rt ;	t9 =t9 +rt ;	copies of t13 in add6     , xmm6\
			it =t14;	t14=t10-it ;	t10=t10+it ;	copies of t14 in add6+0x10, xmm7\
			\
			rt =t15;	t15=t11-t16;	t11=t11+t16;	copies of t15 in add7     , xmm4\
						t16=t12+rt ;	t12=t12-rt ;	copies of t16 in add7+0x10, xmm5\
			*/\
			/* Move outputs t11,12 into a[j1,j2+p5], first doing the addsub and mul by ISRT2: */\
			/* Move outputs t15,16 into a[j1,j2+p7], first doing the addsub and mul by ISRT2: */\
			\
			__asm	addpd	xmm6,xmm2		/* xmm6 <- ~t9  */\
			__asm	addpd	xmm7,xmm3		/* xmm7 <- ~t10 */\
			__asm	subpd	xmm2,[ecx     ]	/* xmm2 <- ~t13 */\
			__asm	subpd	xmm3,[ecx+0x10]	/* xmm3 <- ~t14 */\
			/* Move t13,14 into a[j1,j2+p6] */\
			__asm	movaps	[ecx     ],xmm2	/* add6r <- ~t13 */\
			__asm	movaps	[ecx+0x10],xmm3	/* add6i <- ~t14 */\
			\
			__asm	movaps	xmm2,xmm4	/* xmm2 <- copy of t15 */\
			__asm	movaps	xmm3,xmm5	/* xmm3 <- copy of t16 */\
			__asm	addpd	xmm5,xmm0	/* xmm5 <-~t11 */\
			__asm	subpd	xmm0,xmm3	/* xmm0 <-~t15 */\
			__asm	addpd	xmm4,xmm1	/* xmm4 <-~t16 */\
			__asm	subpd	xmm1,xmm2	/* xmm1 <-~t12 */\
			\
			__asm	movaps	xmm2,xmm5	/* xmm2 <- copy of~t11 */\
			__asm	movaps	xmm3,xmm1	/* xmm3 <- copy of~t12 */\
			__asm	addpd	xmm5,xmm1	/* xmm5 <-~(t11+t12), xmm1 FREE */\
			__asm	mov	edi,isrt2\
			__asm	movaps	xmm1,[edi]	/* xmm1 <- ISRT2 */\
			__asm	subpd	xmm2,xmm3	/* xmm2 <-~(t11-t12), xmm3 FREE */\
			__asm	mulpd	xmm5,xmm1	/* xmm5 <- (t11+t12)*ISRT2 */\
			__asm	mulpd	xmm2,xmm1	/* xmm2 <- (t11-t12)*ISRT2 */\
			__asm	movaps	xmm3,xmm0	/* xmm3 <- copy of~t15 */\
			\
			__asm	movaps	[ebx     ],xmm5	/* add5r<- (t11+t12)*ISRT2, xmm5 FREE */\
			__asm	movaps	xmm5,xmm4	/* xmm5 <- copy of~t16 */\
			__asm	addpd	xmm0,xmm4	/* xmm0 <-~(t15+t16) */\
			__asm	movaps	[ebx+0x10],xmm2	/* add5i<- (t11-t12)*ISRT2 */\
			__asm	subpd	xmm3,xmm5	/* xmm3 <-~(t15-t16) */\
			__asm	mulpd	xmm0,xmm1	/* xmm0 <- (t15+t16)*ISRT2 */\
			__asm	mulpd	xmm3,xmm1	/* xmm3 <- (t15-t16)*ISRT2 */\
			__asm	movaps	[edx     ],xmm0	/* add7r<- (t15+t16)*ISRT2 */\
			__asm	movaps	[edx+0x10],xmm3	/* add7i<- (t15-t16)*ISRT2 */\
			\
			/** GPRs: ***** SSE Regs: ******** Main array: *************\
			* eax, add4    xmm0 <- unused    add0 <- unused            *\
			* ebx, add5    xmm1 <- unused    add1 <- unused            *\
			* ecx, add6    xmm2 <- unused    add2 <- unused            *\
			* edx, add7    xmm3 <- unused    add3 <- unused            *\
			* esi,   p4    xmm4 <- unused    add4 <- unused            *\
			* edi,isrt2    xmm5 <- unused    add5 <- (t11+-t12)*ISRT2  *\
			*              xmm6 <- t9        add6 <-  t13,14           *\
			*              xmm7 <- t10       add7 <- (t15+-t16)*ISRT2  *\
			***********************************************************/\
			\
			/**************** 1st of the 2 length-4 subtransforms... **************/\
			/*\
			t1 =a[j1   ];	t2 =a[j2   ];\
			rt =a[j1+p1];	it =a[j2+p1];\
			t3 =t1 -rt;		t4 =t2 -it;\
			t1 =t1 +rt;		t2 =t2 +it;\
			*/\
			__asm	mov	edi, eax	/* &a[j1+p4] */\
			__asm	sub	eax, esi	/* &a[j1+p0] */\
			__asm	sub	ebx, esi	/* &a[j1+p1] */\
			__asm	sub	ecx, esi	/* &a[j1+p2] */\
			__asm	sub	edx, esi	/* &a[j1+p3] */\
			\
			__asm	movaps	xmm0,[eax]		/* xmm0 <- a[j1   ] = t1 */\
			__asm	movaps	xmm1,[eax+0x10]	/* xmm1 <- a[j2   ] = t2 */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- copy of t1    */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- copy of t2    */\
			__asm	addpd	xmm2,[ebx]		/*~xmm2 <- t1  */\
			__asm	addpd	xmm3,[ebx+0x10]	/*~xmm3 <- t2  */\
			__asm	subpd	xmm0,[ebx]		/*~xmm0 <- t3  */\
			__asm	subpd	xmm1,[ebx+0x10]	/*~xmm1 <- t4  */\
			\
			/* Move   t9,  t10 into a[j1,j2   ] in anticipation of final outputs t1+-t9, t2+-t10 which will go there: */\
			__asm	movaps	[eax     ],xmm6\
			__asm	movaps	[eax+0x10],xmm7	/* add0 <-  t9,t10 */\
			__asm	addpd	xmm6,xmm6		/* xmm6 <- 2*t9  */\
			__asm	addpd	xmm7,xmm7		/* xmm7 <- 2*t10 */\
			/* Move 2*t9,2*t10 into a[j1,j2+p4] in anticipation of final outputs t1+-t9, t2+-t10 which will go there: */\
			__asm	movaps	[edi     ],xmm6	/* add4 <- 2*t9,2*t10  */\
			__asm	movaps	[edi+0x10],xmm7	/* xmm4-7 FREE */\
			/*\
			t5 =a[j1+p2];	t6 =a[j2+p2];\
			rt =a[j1+p3];	it =a[j2+p3];\
			t7 =t5 -rt;		t8 =t6 -it;\
			t5 =t5 +rt;		t6 =t6 +it;\
			*/\
			__asm	movaps	xmm4,[ecx]		/* xmm4 <- a[j1+p2] = t5 */\
			__asm	movaps	xmm5,[ecx+0x10]	/* xmm5 <- a[j2+p2] = t6 */\
			__asm	movaps	xmm6,xmm4		/* xmm6 <- copy of t5    */\
			__asm	movaps	xmm7,xmm5		/* xmm7 <- copy of t6    */\
			__asm	addpd	xmm6,[edx]		/*~xmm6 <- t5  */\
			__asm	addpd	xmm7,[edx+0x10]	/*~xmm7 <- t6  */\
			__asm	subpd	xmm4,[edx]		/*~xmm4 <- t7  */\
			__asm	subpd	xmm5,[edx+0x10]	/*~xmm5 <- t8  */\
			/* Copy t5,6 into main-array slots a[j1,j2+p2] */\
			__asm	movaps	[ecx     ],xmm6\
			__asm	movaps	[ecx+0x10],xmm7	/* add2 <-  t5,t6 */\
			\
			/** GPRs: ***** SSE Regs: ******** Main array: *************\
			* eax, add0    xmm0 <- t3        add0 <-  t9,t10           *\
			* ebx, add1    xmm1 <- t4        add1 <- unused            *\
			* ecx, add2    xmm2 <- t1        add2 <-  t5,t6            *\
			* edx, add3    xmm3 <- t2        add3 <- unused            *\
			* esi,   p4    xmm4 <- t7        add4 <- 2*t9,2*t10        *\
			* edi, add4    xmm5 <- t8        add5 <- (t11+-t12)*ISRT2  *\
			*              xmm6 <- t5        add6 <-  t13,14           *\
			*              xmm7 <- t6        add7 <- (t15+-t16)*ISRT2  *\
			***********************************************************/\
			\
			/*\
			rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;	copies of t5 in add2     , xmm6\
			it =t6;	t6 =t2 -it;	t2 =t2 +it;	copies of t6 in add2+0x10, xmm7\
			\
			rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;\
					t8 =t4 +rt;	t4 =t4 -rt;\
			*/\
			__asm	addpd	xmm6,xmm2		/* xmm6 <- ~t1 */\
			__asm	addpd	xmm7,xmm3		/* xmm7 <- ~t2 */\
			__asm	subpd	xmm2,[ecx     ]	/* xmm2 <- ~t5 */\
			__asm	subpd	xmm3,[ecx+0x10]	/* xmm3 <- ~t6 */\
			\
			/* Compute and dump first 2 outputs now, in order to free up 2 registers: */\
			__asm	addpd	xmm6,[eax     ]	/* t1+t9 */\
			__asm	addpd	xmm7,[eax+0x10]	/* t2+t10*/\
			__asm	movaps	[eax     ],xmm6	/* a[j1   ], DONE. */\
			__asm	movaps	[eax+0x10],xmm7	/* a[j2   ], DONE. */\
			\
			__asm	subpd	xmm6,[edi     ]	/* t1-t9  = [t1+t9 ] - 2*t9  */\
			__asm	subpd	xmm7,[edi+0x10]	/* t2-t10 = [t2+t10] - 2*t10 */\
			__asm	movaps	[edi     ],xmm6	/* Spill t1-t9  to a[j1+p4]. */\
			__asm	movaps	[edi+0x10],xmm7	/* Spill t2-t10 to a[j2+p4]. */\
			__asm	movaps	[edx     ],xmm6	/* Copy: t1-t9  -> a[j1+p3]. */\
			__asm	movaps	[edx+0x10],xmm7	/* Copy: t2-t10 -> a[j2+p3]. */\
			\
			__asm	movaps	xmm6,xmm4		/* xmm6 <- copy of t7 */\
			__asm	movaps	xmm7,xmm5		/* xmm7 <- copy of t8 */\
			__asm	addpd	xmm5,xmm0	/* xmm5 <- ~t3 */\
			__asm	subpd	xmm0,xmm7	/* xmm0 <- ~t7 */\
			__asm	addpd	xmm4,xmm1	/* xmm4 <- ~t8 */\
			__asm	subpd	xmm1,xmm6	/* xmm1 <- ~t4 */\
			\
			/** GPRs: ***** SSE Regs: ******** Main array: *************\
			* eax, ----    xmm0 <- t7        add0 <- DONE              *\
			* ebx, add1    xmm1 <- t4        add1 <- unused            *\
			* ecx, add2    xmm2 <- t5        add2 <- unused            *\
			* edx, add3    xmm3 <- t6        add3 <- t1-t9,t2-t10      *\
			* esi,   p4    xmm4 <- t8        add4 <- t1-t9,t2-t10      *\
			* edi, add4    xmm5 <- t3        add5 <- (t11+-t12)*ISRT2  *\
			*              xmm6 <- unused    add6 <-  t13,14           *\
			*              xmm7 <- unused    add7 <- (t15+-t16)*ISRT2  *\
			***********************************************************/\
			\
			/* Now combine the two half-transforms & store outputs back into original array slots: */\
			/*\
			a[j1   ]=t1+t9;				a[j2   ]=t2+t10;	already done\
			~t1     =t1-t9;				~t2     =t2-t10;	copies in add3[edx],add4[edi]\
			a[j1+p4]=~t1*c4 +~t2*s4;	a[j2+p4]=~t2*c4 -~t1*s4;\
			*/\
			__asm	mov	eax,   c4\
			__asm	movaps	xmm6,[edi     ]	/* t1-t9  <- a[j1+p4]. */\
			__asm	movaps	xmm7,[edi+0x10]	/* t2-t10 <- a[j2+p4]. */\
			__asm	mulpd	xmm6,[eax+0x10]	/*~t1*s4 */\
			__asm	mulpd	xmm7,[eax+0x10]	/*~t2*s4 */\
			__asm	movaps	[edi     ],xmm6	/* a[j1+p4] <- ~t1*s4 */\
			__asm	movaps	[edi+0x10],xmm7	/* a[j2+p4] <- ~t2*s4 */\
			\
			__asm	movaps	xmm6,[edx     ]	/* xmm6 <- cpy ~t1 */\
			__asm	movaps	xmm7,[edx+0x10]	/* xmm7 <- cpy ~t2 */\
			__asm	mulpd	xmm6,[eax     ]	/*~t1*c4 */\
			__asm	mulpd	xmm7,[eax     ]	/*~t2*c4 */\
			__asm	addpd	xmm6,[edi+0x10]	/* ~t1*c4 +~t2*s4 */\
			__asm	subpd	xmm7,[edi     ]	/* ~t2*c4 -~t1*s4 */\
			\
			__asm	movaps	[edi     ],xmm6	/* a[j1+p4] <- ~t1*c4 +~t2*s4 */\
			__asm	movaps	[edi+0x10],xmm7	/* a[j2+p4] <- ~t2*c4 -~t1*s4 */\
			\
			/*\
			rt=(t11+t12)*ISRT2;			it=(t11-t12)*ISRT2;	precomputed, in add5\
			t11     =t3+rt;				t12       =t4-it;\
			t3      =t3-rt;				t4        =t4+it;\
			a[j1+p1]=t11*c1 +t12*s1;	a[j2+p1]=t12*c1 -t11*s1;\
			a[j1+p5]=t3 *c5 +t4 *s5;	a[j2+p5]=t4 *c5 -t3 *s5;\
			*/\
			__asm	mov	edi, ebx	/* add1 cpy */\
			__asm	add	edi, esi	/* add5 */\
			__asm	add	eax,0x60	/* c1 */\
			__asm	movaps	xmm6,xmm5		/* xmm6 <- copy of t3 */\
			__asm	movaps	xmm7,xmm1		/* xmm7 <- copy of t4 */\
			__asm	addpd	xmm5,[edi     ]	/* ~t11 */\
			__asm	subpd	xmm1,[edi+0x10]	/* ~t12 */\
			__asm	subpd	xmm6,[edi     ]	/* ~t3  */\
			__asm	addpd	xmm7,[edi+0x10]	/* ~t4  */\
			__asm	movaps	[edi     ],xmm6	/* spill ~t3,t4 to add5 */\
			__asm	movaps	[edi+0x10],xmm7\
			__asm	movaps	xmm6,xmm5	/* xmm6 <- cpy ~t11*/\
			__asm	movaps	xmm7,xmm1	/* xmm7 <- cpy ~t12*/\
			\
			__asm	mulpd	xmm5,[eax     ]	/* t11*c1 */\
			__asm	mulpd	xmm1,[eax     ]	/* t12*c1 */\
			__asm	mulpd	xmm6,[eax+0x10]	/* t11*s1 */\
			__asm	mulpd	xmm7,[eax+0x10]	/* t12*s1 */\
			__asm	subpd	xmm1,xmm6	/* t12*c1 - t11*s1 */\
			__asm	addpd	xmm5,xmm7	/* t11*c1 + t12*s1 */\
			__asm	movaps	[ebx+0x10],xmm1	/* a[j2+p1] */\
			__asm	movaps	[ebx     ],xmm5	/* a[j1+p1] */\
			\
			__asm	add	eax,0x20	/* c5 */\
			__asm	movaps	xmm5,[edi     ]	/* t3 */\
			__asm	movaps	xmm1,[edi+0x10]	/* t4 */\
			__asm	movaps	xmm6,xmm5	/* xmm6 <- copy t3 */\
			__asm	movaps	xmm7,xmm1	/* xmm7 <- copy t4 */\
			\
			__asm	mulpd	xmm5,[eax     ]	/* t3 *c5 */\
			__asm	mulpd	xmm1,[eax     ]	/* t4 *c5 */\
			__asm	mulpd	xmm6,[eax+0x10]	/* t3 *s5 */\
			__asm	mulpd	xmm7,[eax+0x10]	/* t4 *s5 */\
			__asm	subpd	xmm1,xmm6	/* t4*c5 - t3*s5 */\
			__asm	addpd	xmm5,xmm7	/* t3*c5 + t4*s5 */\
			__asm	movaps	[edi+0x10],xmm1	/* a[j2+p5] */\
			__asm	movaps	[edi     ],xmm5	/* a[j1+p5] */	/*** xmm1,5,6,7 FREE ***/\
			\
			/*\
			rt      =t5+t14;			it        =t6-t13;\
			t5      =t5-t14;			t6        =t6+t13;\
			a[j1+p2]=rt *c2 +it *s2;	a[j2+p2]=it *c2 -rt *s2;\
			a[j1+p6]=t5 *c6 +t6 *s6;	a[j2+p6]=t6 *c6 -t5 *s6;\
			*/\
			__asm	mov	edi, ecx	/* add2 cpy */\
			__asm	add	edi, esi	/* add6 */\
			__asm	sub	eax,0x60	/* c2 */\
			__asm	movaps	xmm6,xmm2		/* xmm6 <- copy of t5 */\
			__asm	movaps	xmm7,xmm3		/* xmm7 <- copy of t6 */\
			__asm	addpd	xmm2,[edi+0x10]	/*  rt */\
			__asm	subpd	xmm3,[edi     ]	/*  it */\
			__asm	subpd	xmm6,[edi+0x10]	/* ~t5  */\
			__asm	addpd	xmm7,[edi     ]	/* ~t6  */\
			__asm	movaps	xmm1,xmm2	/* xmm1 <- cpy rt */\
			__asm	movaps	xmm5,xmm3	/* xmm5 <- cpy it */\
			\
			__asm	mulpd	xmm2,[eax     ]	/* rt*c2 */\
			__asm	mulpd	xmm3,[eax     ]	/* it*c2 */\
			__asm	mulpd	xmm1,[eax+0x10]	/* rt*s2 */\
			__asm	mulpd	xmm5,[eax+0x10]	/* it*s2 */\
			__asm	subpd	xmm3,xmm1	/* it*c2 - rt*s2 */\
			__asm	addpd	xmm2,xmm5	/* rt*c2 + it*s2 */\
			__asm	movaps	[ecx+0x10],xmm3	/* a[j2+p2] */\
			__asm	movaps	[ecx     ],xmm2	/* a[j1+p2] */\
			\
			__asm	add	eax,0x20	/* c6 */\
			__asm	movaps	xmm1,xmm6	/* xmm1 <- cpy t5 */\
			__asm	movaps	xmm5,xmm7	/* xmm5 <- cpy t6 */\
			\
			__asm	mulpd	xmm6,[eax     ]	/* t5*c6 */\
			__asm	mulpd	xmm7,[eax     ]	/* t6*c6 */\
			__asm	mulpd	xmm1,[eax+0x10]	/* t5*s6 */\
			__asm	mulpd	xmm5,[eax+0x10]	/* t6*s6 */\
			__asm	subpd	xmm7,xmm1	/* t6*c6 - t5*s6 */\
			__asm	addpd	xmm6,xmm5	/* t5*c6 + t6*s6 */\
			__asm	movaps	[edi+0x10],xmm7	/* a[j2+p6] */\
			__asm	movaps	[edi     ],xmm6	/* a[j1+p6] */\
			\
			/*\
			rt=(t15-t16)*ISRT2;			it=(t15+t16)*ISRT2;	precomputed, it,rt] in add7; NOTE reversed order!\
			t15     =t7-rt;				t16       =t8-it;\
			t7      =t7+rt;				t8        =t8+it;\
			a[j1+p3]=t15*c3 +t16*s3;	a[j2+p3]=t16*c3 -t15*s3;\
			a[j1+p7]=t7 *c7 +t8 *s7;	a[j2+p7]=t8 *c7 -t7 *s7;\
			*/\
			__asm	mov	edi, edx	/* add3 cpy */\
			__asm	add	edi, esi	/* add7 */\
			__asm	add	eax,0x60	/* c3 */\
			__asm	movaps	xmm6,xmm0		/* xmm6 <- copy of t7 */\
			__asm	movaps	xmm7,xmm4		/* xmm7 <- copy of t8 */\
			__asm	subpd	xmm0,[edi+0x10]	/* ~t15 */\
			__asm	subpd	xmm4,[edi     ]	/* ~t16 */\
			__asm	addpd	xmm6,[edi+0x10]	/* ~t7  */\
			__asm	addpd	xmm7,[edi     ]	/* ~t8  */\
			\
			__asm	movaps	xmm1,xmm0	/* xmm6 <- cpy t15*/\
			__asm	movaps	xmm5,xmm4	/* xmm7 <- cpy t16*/\
			\
			__asm	mulpd	xmm0,[eax     ]	/* t15*c3 */\
			__asm	mulpd	xmm4,[eax     ]	/* t16*c3 */\
			__asm	mulpd	xmm1,[eax+0x10]	/* t15*s3 */\
			__asm	mulpd	xmm5,[eax+0x10]	/* t16*s3 */\
			__asm	subpd	xmm4,xmm1	/* t16*c3 -t15*s3 */\
			__asm	addpd	xmm0,xmm5	/* t15*c3 +t16*s3 */\
			__asm	movaps	[edx+0x10],xmm4	/* a[j2+p3] */\
			__asm	movaps	[edx     ],xmm0	/* a[j1+p3] */\
			\
			__asm	add	eax,0x20	/* c7 */\
			__asm	movaps	xmm1,xmm6	/* xmm1 <- cpy t7 */\
			__asm	movaps	xmm5,xmm7	/* xmm5 <- cpy t8 */\
			\
			__asm	mulpd	xmm6,[eax     ]	/* t7*c7 */\
			__asm	mulpd	xmm7,[eax     ]	/* t8*c7 */\
			__asm	mulpd	xmm1,[eax+0x10]	/* t7*s7 */\
			__asm	mulpd	xmm5,[eax+0x10]	/* t8*s7 */\
			__asm	subpd	xmm7,xmm1	/* t8*c7 - t7*s7 */\
			__asm	addpd	xmm6,xmm5	/* t7*c7 + t8*s7 */\
			__asm	movaps	[edi+0x10],xmm7	/* a[j2+p7] */\
			__asm	movaps	[edi     ],xmm6	/* a[j1+p7] */\
		}
												/* Totals: 219 load/store [89 movaps], 68 add/subpd, 32 mulpd */

/******************************************************************************************************************/

		/* DIT version of SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS: switch r4<->rc, r5<->rd, r6<->re, r7<->rf. */
		#define SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(__r0,__r1,__r2,__r3,__r4,__r5,__r6,__r7)\
		{\
			__asm	mov	eax, __r0\
			__asm	mov	ebx, __r4\
			__asm	mov	ecx, __r2\
			__asm	mov	edx, __r6\
			\
			__asm	movaps	xmm0,[eax     ]	/* t0 */			__asm	movaps	xmm4,[ecx     ]	/* t4 */\
			__asm	movaps	xmm1,[eax+0x10]	/* t1 */			__asm	movaps	xmm5,[ecx+0x10]	/* t5 */\
			__asm	movaps	xmm2,[ebx     ]	/* cpy t8 */		__asm	movaps	xmm7,[edx+0x10]	/* td */\
			__asm	movaps	xmm3,[ebx+0x10]	/* cpy t9 */		__asm	movaps	xmm6,[edx     ]	/* tc */\
			__asm	subpd	xmm0,xmm2		/* t8 = 0-8 */		__asm	subpd	xmm4,xmm7		/* tc = 4-d */\
			__asm	subpd	xmm1,xmm3		/* t9 = 1-9 */		__asm	subpd	xmm5,xmm6		/* t5 = 5-c */\
			__asm	addpd	xmm2,xmm2		/*      2*8 */		__asm	addpd	xmm7,xmm7		/*      2*d */\
			__asm	addpd	xmm3,xmm3		/*      2*9 */		__asm	addpd	xmm6,xmm6		/*      2*c */\
			__asm	addpd	xmm2,xmm0		/* t0 = 0+8 */		__asm	addpd	xmm7,xmm4		/* t4 = 4+d */\
			__asm	addpd	xmm3,xmm1		/* t1 = 1+9 */		__asm	addpd	xmm6,xmm5		/* td = 5+c */\
			\
			__asm	movaps	[ebx     ],xmm0	/* t8 */			__asm	movaps	[edx     ],xmm4	/* t4 */\
			__asm	movaps	[ebx+0x10],xmm1	/* t9 */			__asm	movaps	[ecx+0x10],xmm5	/* td */\
			__asm	movaps	[eax     ],xmm2	/* t0 */			__asm	movaps	[ecx     ],xmm7	/* tc */\
			__asm	movaps	[eax+0x10],xmm3	/* t1 */			__asm	movaps	[edx+0x10],xmm6	/* t5 */\
			\
			__asm	mov	eax, __r1\
			__asm	mov	ebx, __r5\
			__asm	mov	ecx, __r3\
			__asm	mov	edx, __r7\
			\
			__asm	movaps	xmm0,[eax     ]	/* t2 */			__asm	movaps	xmm4,[ecx     ]	/* t6 */\
			__asm	movaps	xmm1,[eax+0x10]	/* t3 */			__asm	movaps	xmm5,[ecx+0x10]	/* t7 */\
			__asm	movaps	xmm2,[ebx     ]	/* cpy ta */		__asm	movaps	xmm7,[edx+0x10]	/* tf */\
			__asm	movaps	xmm3,[ebx+0x10]	/* cpy tb */		__asm	movaps	xmm6,[edx     ]	/* te */\
			__asm	subpd	xmm0,xmm2		/* ta = 2-a */		__asm	subpd	xmm4,xmm7		/* te = 6-f */\
			__asm	subpd	xmm1,xmm3		/* tb = 3-b */		__asm	subpd	xmm5,xmm6		/* t7 = 7-e */\
			__asm	addpd	xmm2,xmm2		/*      2*a */		__asm	addpd	xmm7,xmm7		/*      2*f */\
			__asm	addpd	xmm3,xmm3		/*      2*b */		__asm	addpd	xmm6,xmm6		/*      2*e */\
			__asm	addpd	xmm2,xmm0		/* t2 = 2+a */		__asm	addpd	xmm7,xmm4		/* t6 = 6+f */\
			__asm	addpd	xmm3,xmm1		/* t3 = 3+b */		__asm	addpd	xmm6,xmm5		/* tf = 7+e */\
			\
			__asm	movaps	[ebx     ],xmm0	/* ta */			__asm	movaps	[edx     ],xmm4	/* t6 */\
			__asm	movaps	[ebx+0x10],xmm1	/* tb */			__asm	movaps	[ecx+0x10],xmm5	/* tf */\
			__asm	movaps	[eax     ],xmm2	/* t2 */			__asm	movaps	[ecx     ],xmm7	/* te */\
			__asm	movaps	[eax+0x10],xmm3	/* t3 */			__asm	movaps	[edx+0x10],xmm6	/* t7 */\
		}

		/* Does the 2nd pair of radix-2 butterflies-with-twiddles of a radix-4 DIT, e.g of form like so:

			rt       =t0r+t2r;			it       =t0i+t2i;					rt       =t1r+t3i;		it         =t1i-t3r;
			t2r      =t0r-t2r;			t2i      =t0i-t2i;					t3i      =t1r-t3i;		t3r        =t1i+t3r;
			a[jt+p0]=rt *c0+it *s0;	a[jp+p0]=it *c0-rt *s0;			a[jt+p1]=rt *c1+it *s1;	a[jp+p1]=it *c1-rt *s1;
			a[jt+p2]=t2r*c2+t2i*s2;	a[jp+p8]=t2r*c2-t2i*s2;			a[jt+p3]=t3i*c3+t3r*s3;	a[jp+p3]=t3r*c3-t3i*s3;

		Assumes:

			- Addresses of the 4 real outputs a[jt+p0,1,2,3] enter in add0,1,2,3;
			- t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i enter in SSE registers xmm2,3,0,1,4,5,6,7 [don't ask ;)];
			- Address of t0r is in eax, t1r is in edx.
			- Address of c0,c1,c2,c3 are adjacent in memory.
		*/
		#define SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF(__add0, __add1, __add2, __add3, __c0)\
		{\
			__asm	mov	ebx, add0											__asm	mov	edi, add1\
			__asm	mov	esi, __c0\
			\
			__asm	subpd	xmm2,xmm4		/*~t2r <- t0r-t2r */			__asm	subpd	xmm0,xmm7		/*~t3i <- t1r-t3i */\
			__asm	subpd	xmm3,xmm5		/*~t2i <- t0i-t2i */			__asm	subpd	xmm1,xmm6		/* it  <- t1i-t3r */\
			__asm	addpd	xmm4,xmm4		/*          2*t2r */			__asm	addpd	xmm7,xmm7		/*          2*t3i */\
			__asm	addpd	xmm5,xmm5		/*          2*t2i */			__asm	addpd	xmm6,xmm6		/*          2*t3r */\
			__asm	addpd	xmm4,xmm2		/* rt  <- t0r+t2r */			__asm	addpd	xmm7,xmm0		/* rt  <- t1r+t3i */\
			__asm	addpd	xmm5,xmm3		/* it  <- t0i+t2i */			__asm	addpd	xmm6,xmm1		/*~t3r <- t1i+t3r */\
			__asm	movaps	[eax      ],xmm2	/* store ~t2r */			__asm	movaps	[edx+0x010],xmm0	/* store ~t3i */\
			__asm	movaps	[eax+0x010],xmm3	/* store ~t2i */			__asm	movaps	[edx      ],xmm6	/* store ~t3r */\
			__asm	movaps	xmm2,xmm4		/* rt copy */					__asm	movaps	xmm0,xmm7		/* rt copy */\
			__asm	movaps	xmm3,xmm5		/* it copy */					__asm	movaps	xmm6,xmm1		/* it copy */\
			__asm	mulpd	xmm4,[esi     ]	/* rt * c0 */					__asm	mulpd	xmm7,[esi+0x20]	/* rt* c1 */\
			__asm	mulpd	xmm5,[esi     ]	/* it * c0 */					__asm	mulpd	xmm1,[esi+0x20]	/* it* c1 */\
			__asm	mulpd	xmm2,[esi+0x10]	/* rt * s0 */					__asm	mulpd	xmm0,[esi+0x30]	/* rt* s1 */\
			__asm	mulpd	xmm3,[esi+0x10]	/* it * s0 */					__asm	mulpd	xmm6,[esi+0x30]	/* it* s1 */\
			__asm	subpd	xmm5,xmm2	/* xmm5 <- im */					__asm	subpd	xmm1,xmm0	/* xmm1 <- im */\
			__asm	addpd	xmm4,xmm3	/* xmm4 <- re */					__asm	addpd	xmm7,xmm6	/* xmm7 <- re */\
			\
			__asm	movaps	[ebx+0x10],xmm5	/* a[jp+p0] */					__asm	movaps	[edi+0x10],xmm1	/* a[jp+p1] */\
			__asm	movaps	[ebx     ],xmm4	/* a[jt+p0] */					__asm	movaps	[edi     ],xmm7	/* a[jt+p1] */\
			\
			__asm	mov	ebx, add2											__asm	mov	edi, add3\
			__asm	add	esi, 0x40	/* c2 */\
			\
			__asm	movaps	xmm4,[eax      ]	/* load ~t2r */				__asm	movaps	xmm0,[edx+0x010]	/* load ~t3i */\
			__asm	movaps	xmm5,[eax+0x010]	/* load ~t2i */				__asm	movaps	xmm6,[edx      ]	/* load ~t3r */\
			__asm	movaps	xmm2,xmm4		/* re copy */					__asm	movaps	xmm1,xmm0		/*t3i copy */\
			__asm	movaps	xmm3,xmm5		/* im copy */					__asm	movaps	xmm7,xmm6		/*t3r copy */\
			__asm	mulpd	xmm4,[esi     ]	/* re * c2 */					__asm	mulpd	xmm0,[esi+0x20]	/*t3i* c3 */\
			__asm	mulpd	xmm5,[esi     ]	/* im * c2 */					__asm	mulpd	xmm6,[esi+0x20]	/*t3r* c3 */\
			__asm	mulpd	xmm2,[esi+0x10]	/* re * c2 */					__asm	mulpd	xmm1,[esi+0x30]	/*t3i* s3 */\
			__asm	mulpd	xmm3,[esi+0x10]	/* im * c2 */					__asm	mulpd	xmm7,[esi+0x30]	/*t3r* s3 */\
			__asm	subpd	xmm5,xmm2	/* xmm5 <- im */					__asm	subpd	xmm6,xmm1	/* xmm6 <- im */\
			__asm	addpd	xmm4,xmm3	/* xmm4 <- re */					__asm	addpd	xmm0,xmm7	/* xmm7 <- re */\
			\
			__asm	movaps	[ebx+0x10],xmm5	/* a[jp+p2] */					__asm	movaps	[edi+0x10],xmm6	/* a[jp+p3] */\
			__asm	movaps	[ebx     ],xmm4	/* a[jt+p2] */					__asm	movaps	[edi     ],xmm0	/* a[jt+p3] */\
		}

		/* Full-inline-asm version of SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF.

		Assumes:
			- Addresses of the real outputs a[jt+p0,1] enter in eax,ebx;
			- Address of t0r is in ecx, t1r is in edx.
			- Base offset p2 in edi, these must remain unchanged, so esi is available for temporary storage
			- t0r,t0i,t1r,t1i,t2r,t2i,t3r,t3i enter in SSE registers xmm2,3,0,1,4,5,6,7 [don't ask ;)];
			- Address of c0,c1,c2,c3 are adjacent in memory.
		*/
		#define SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(__c0)\
		{\
			__asm	mov	esi, __c0\
			\
			__asm	subpd	xmm2,xmm4		/*~t2r <- t0r-t2r */			__asm	subpd	xmm0,xmm7		/*~t3i <- t1r-t3i */\
			__asm	subpd	xmm3,xmm5		/*~t2i <- t0i-t2i */			__asm	subpd	xmm1,xmm6		/* it  <- t1i-t3r */\
			__asm	addpd	xmm4,xmm4		/*          2*t2r */			__asm	addpd	xmm7,xmm7		/*          2*t3i */\
			__asm	addpd	xmm5,xmm5		/*          2*t2i */			__asm	addpd	xmm6,xmm6		/*          2*t3r */\
			__asm	addpd	xmm4,xmm2		/* rt  <- t0r+t2r */			__asm	addpd	xmm7,xmm0		/* rt  <- t1r+t3i */\
			__asm	addpd	xmm5,xmm3		/* it  <- t0i+t2i */			__asm	addpd	xmm6,xmm1		/*~t3r <- t1i+t3r */\
			__asm	movaps	[ecx      ],xmm2	/* store ~t2r */			__asm	movaps	[edx+0x010],xmm0	/* store ~t3i */\
			__asm	movaps	[ecx+0x010],xmm3	/* store ~t2i */			__asm	movaps	[edx      ],xmm6	/* store ~t3r */\
			__asm	movaps	xmm2,xmm4		/* rt copy */					__asm	movaps	xmm0,xmm7		/* rt copy */\
			__asm	movaps	xmm3,xmm5		/* it copy */					__asm	movaps	xmm6,xmm1		/* it copy */\
			__asm	mulpd	xmm4,[esi     ]	/* rt * c0 */					__asm	mulpd	xmm7,[esi+0x20]	/* rt* c1 */\
			__asm	mulpd	xmm5,[esi     ]	/* it * c0 */					__asm	mulpd	xmm1,[esi+0x20]	/* it* c1 */\
			__asm	mulpd	xmm2,[esi+0x10]	/* rt * s0 */					__asm	mulpd	xmm0,[esi+0x30]	/* rt* s1 */\
			__asm	mulpd	xmm3,[esi+0x10]	/* it * s0 */					__asm	mulpd	xmm6,[esi+0x30]	/* it* s1 */\
			__asm	subpd	xmm5,xmm2	/* xmm5 <- im */					__asm	subpd	xmm1,xmm0	/* xmm1 <- im */\
			__asm	addpd	xmm4,xmm3	/* xmm4 <- re */					__asm	addpd	xmm7,xmm6	/* xmm7 <- re */\
			\
			__asm	movaps	[eax+0x10],xmm5	/* a[jp+p0] */					__asm	movaps	[ebx+0x10],xmm1	/* a[jp+p1] */\
			__asm	movaps	[eax     ],xmm4	/* a[jt+p0] */					__asm	movaps	[ebx     ],xmm7	/* a[jt+p1] */\
			\
			__asm	add	eax, edi											__asm	add	ebx, edi\
			__asm	add	esi, 0x40	/* c2 */\
			\
			__asm	movaps	xmm4,[ecx      ]	/* load ~t2r */				__asm	movaps	xmm0,[edx+0x010]	/* load ~t3i */\
			__asm	movaps	xmm5,[ecx+0x010]	/* load ~t2i */				__asm	movaps	xmm6,[edx      ]	/* load ~t3r */\
			__asm	movaps	xmm2,xmm4		/* re copy */					__asm	movaps	xmm1,xmm0		/*t3i copy */\
			__asm	movaps	xmm3,xmm5		/* im copy */					__asm	movaps	xmm7,xmm6		/*t3r copy */\
			__asm	mulpd	xmm4,[esi     ]	/* re * c2 */					__asm	mulpd	xmm0,[esi+0x20]	/*t3i* c3 */\
			__asm	mulpd	xmm5,[esi     ]	/* im * c2 */					__asm	mulpd	xmm6,[esi+0x20]	/*t3r* c3 */\
			__asm	mulpd	xmm2,[esi+0x10]	/* re * c2 */					__asm	mulpd	xmm1,[esi+0x30]	/*t3i* s3 */\
			__asm	mulpd	xmm3,[esi+0x10]	/* im * c2 */					__asm	mulpd	xmm7,[esi+0x30]	/*t3r* s3 */\
			__asm	subpd	xmm5,xmm2	/* xmm5 <- im */					__asm	subpd	xmm6,xmm1	/* xmm6 <- im */\
			__asm	addpd	xmm4,xmm3	/* xmm4 <- re */					__asm	addpd	xmm0,xmm7	/* xmm7 <- re */\
			\
			__asm	movaps	[eax+0x10],xmm5	/* a[jp+p2] */					__asm	movaps	[ebx+0x10],xmm6	/* a[jp+p3] */\
			__asm	movaps	[eax     ],xmm4	/* a[jt+p2] */					__asm	movaps	[ebx     ],xmm0	/* a[jt+p3] */\
		}

		/* Variant of SSE2_RADIX4_DIF_4TWIDDLE needed by the radix32_wrapper_square routine */
		/* Assumes addresses __add0, __add1 enter in eax, ebx, respectively: */
		#define SSE2_RADIX4_DIF_4WRAPPER(__c00,__c08,__c10,__c18, __tmp)\
		{\
			__asm	mov	ecx,__tmp\
			\
			/* Do the p00,p10 combo: */\
			__asm	mov	edx,__c00\
			\
			/* For interleaved [j1,j2] version, replace e.g.\
			\
				__asm	movaps	xmm0,[eax+0x40]	// a[jt+p8]\
				__asm	movaps	xmm1,[eax+0x50]	// a[jp+p8]\
				__asm	movaps	xmm2,[eax+0x40]	// xmm2 <- cpy a[jt+p8]\
				__asm	movaps	xmm3,[eax+0x50]	// xmm3 <- cpy a[jp+p8]\
			\
			by the following:\
			*/\
			/* Real parts: */\
			__asm	movaps		xmm6,[eax     ]	/* a[j1+p0 ], this is the scratch xmm register */\
			__asm	movaps		xmm0,[eax     ]	/* a[j1+p0 ], this is the active  xmm register */\
			__asm	unpckhpd	xmm6,[ebx     ]	/* a[j2+p0 ] gets read twice */\
			__asm	unpcklpd	xmm0,[ebx     ]	/* a[jt+p0 ] */\
			__asm	movaps	[ecx+0x200],xmm6	/* Store hi real in t00+32 */\
			/* Imag parts: */\
			__asm	movaps		xmm7,[eax+0x10]\
			__asm	movaps		xmm1,[eax+0x10]\
			__asm	unpckhpd	xmm7,[ebx+0x10]\
			__asm	unpcklpd	xmm1,[ebx+0x10]	/* a[jp+p0 ] */\
			__asm	movaps	[ecx+0x210],xmm7	/* Store hi imag in t01+32 */\
			\
			__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy a[jt+p0] */\
			__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy a[jp+p0] */\
			\
			/* From here on, things are identical to the code in radix32_dif_pass: */\
			__asm	mulpd	xmm0,[edx     ]	/* a[jt+p0]*c0 */\
			__asm	mulpd	xmm1,[edx     ]	/* a[jp+p0]*c0 */\
			__asm	mulpd	xmm2,[edx+0x10]	/* a[jt+p0]*s0 */\
			__asm	mulpd	xmm3,[edx+0x10]	/* a[jp+p0]*s0 */\
			__asm	addpd	xmm1,xmm2	/* xmm1 <- t01*/\
			__asm	subpd	xmm0,xmm3	/* xmm0 <- t00*/\
			__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy t00*/\
			__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy t01*/\
			\
			__asm	mov	edx,__c10\
			/* Real parts: */\
			__asm	movaps		xmm6,[eax+0x100]	/* a[j1+p10], this is the scratch xmm register  */\
			__asm	movaps		xmm4,[eax+0x100]	/* a[j1+p10], this is the active  xmm register */\
			__asm	unpckhpd	xmm6,[ebx+0x100]	/* a[j2+p10] gets read twice */\
			__asm	unpcklpd	xmm4,[ebx+0x100]	/* a[jt+p10] */\
			__asm	movaps	[ecx+0x220],xmm6	/* Store hi real in t11+32 */\
			/* Imag parts: */\
			__asm	movaps		xmm7,[eax+0x110]\
			__asm	movaps		xmm5,[eax+0x110]\
			__asm	unpckhpd	xmm7,[ebx+0x110]\
			__asm	unpcklpd	xmm5,[ebx+0x110]	/* a[jp+p10] */\
			__asm	movaps	[ecx+0x230],xmm7	/* Store hi imag in t12+32 */\
			\
			__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p10] */\
			__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p10] */\
			\
			__asm	mulpd	xmm4,[edx     ]	/* a[jt+p10]*c10 */\
			__asm	mulpd	xmm5,[edx     ]	/* a[jp+p10]*c10 */\
			__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p10]*s10 */\
			__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p10]*s10 */\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- it */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- rt		xmm6,7 free */\
			\
			__asm	addpd	xmm0,xmm4	/* ~t1 <- t1 +rt */\
			__asm	addpd	xmm1,xmm5	/* ~t2 <- t2 +it */\
			__asm	subpd	xmm2,xmm4	/* ~t3 <- t1 -rt */\
			__asm	subpd	xmm3,xmm5	/* ~t4 <- t2 -it	xmm4,5 free */\
			\
			/* Do the p08,p18 [hexadecimal here] combo - do p18 first so register assignments come out in same relative order as for p0,8 */\
			__asm	mov	edx,__c18\
			/* Real parts: */\
			__asm	movaps		xmm6,[eax+0x180]	/* a[j1+p18], this is the scratch xmm register  */\
			__asm	movaps		xmm4,[eax+0x180]	/* a[j1+p18], this is the active  xmm register */\
			__asm	unpckhpd	xmm6,[ebx+0x180]	/* a[j2+p18] gets read twice */\
			__asm	unpcklpd	xmm4,[ebx+0x180]	/* a[jt+p18] */\
			__asm	movaps	[ecx+0x260],xmm6	/* Store hi real in t15+32 */\
			/* Imag parts: */\
			__asm	movaps		xmm7,[eax+0x190]\
			__asm	movaps		xmm5,[eax+0x190]\
			__asm	unpckhpd	xmm7,[ebx+0x190]\
			__asm	unpcklpd	xmm5,[ebx+0x190]	/* a[jp+p18] */\
			__asm	movaps	[ecx+0x270],xmm7	/* Store hi imag in t16+32 */\
			\
			__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p18] */\
			__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p18] */\
			\
			__asm	mulpd	xmm4,[edx     ]	/* a[jt+p18]*c18 */\
			__asm	mulpd	xmm5,[edx     ]	/* a[jp+p18]*c18 */\
			__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p18]*s18 */\
			__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p18]*s18 */\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- it */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- rt		xmm6,7 free */\
			__asm	movaps	[ecx+0x010],xmm5	/* Store it */\
			__asm	movaps	[ecx      ],xmm4	/* Store rt */\
			\
			__asm	mov	edx,__c08\
			/* Real parts: */\
			__asm	movaps		xmm6,[eax+0x080]	/* a[j1+p08], this is the scratch xmm register  */\
			__asm	movaps		xmm4,[eax+0x080]	/* a[j1+p08], this is the active  xmm register */\
			__asm	unpckhpd	xmm6,[ebx+0x080]	/* a[j2+p08] gets read twice */\
			__asm	unpcklpd	xmm4,[ebx+0x080]	/* a[jt+p08] */\
			__asm	movaps	[ecx+0x240],xmm6	/* Store hi real in t13+32 */\
			/* Imag parts: */\
			__asm	movaps		xmm7,[eax+0x090]\
			__asm	movaps		xmm5,[eax+0x090]\
			__asm	unpckhpd	xmm7,[ebx+0x090]\
			__asm	unpcklpd	xmm5,[ebx+0x090]	/* a[jp+p08] */\
			__asm	movaps	[ecx+0x250],xmm7	/* Store hi imag in t14+32 */\
			\
			__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p08] */\
			__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p08] */\
			\
			__asm	mulpd	xmm4,[edx     ]	/* a[jt+p08]*c08 */\
			__asm	mulpd	xmm5,[edx     ]	/* a[jp+p08]*c08 */\
			__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p08]*s08 */\
			__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p08]*s08 */\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- t6 */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- t5 	xmm6,7 free */\
			__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy t6 */\
			__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy t5 */\
			\
			__asm	subpd	xmm4,[ecx      ]	/* ~t7 <- t5 -rt */\
			__asm	subpd	xmm5,[ecx+0x010]	/* ~t8 <- t6 -it */\
			__asm	addpd	xmm6,[ecx      ]	/* ~t5 <- t5 +rt */\
			__asm	addpd	xmm7,[ecx+0x010]	/* ~t6 <- t6 +it */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			/*\
			~t5 =t1 -t5;		~t1 =t1 +t5;\
			~t6 =t2 -t6;		~t2 =t2 +t6;\
			*/\
			__asm	subpd	xmm0,xmm6	/*~t5 */\
			__asm	subpd	xmm1,xmm7	/*~t6 */\
			__asm	movaps	[ecx+0x040],xmm0	/* a[jt+p8 ] <- ~t5 */\
			__asm	movaps	[ecx+0x050],xmm1	/* a[jp+p8 ] <- ~t6 */\
			__asm	addpd	xmm6,xmm6	/* 2*t5 */\
			__asm	addpd	xmm7,xmm7	/* 2*t6 */\
			__asm	addpd	xmm6,xmm0	/*~t1 */\
			__asm	addpd	xmm7,xmm1	/*~t2 */\
			__asm	movaps	[ecx      ],xmm6	/* a[jt    ] <- ~t1 */\
			__asm	movaps	[ecx+0x010],xmm7	/* a[jp    ] <- ~t2 */\
			\
			/*\
			~t7 =t3 +t8;		~t3 =t3 -t8;\
			~t8 =t4 -t7;		~t4 =t4 +t7;\
			*/\
			__asm	subpd	xmm2,xmm5	/*~t3 */\
			__asm	subpd	xmm3,xmm4	/*~t8 */\
			__asm	movaps	[ecx+0x020],xmm2	/* a[jt+p4 ] <- ~t3 */\
			__asm	movaps	[ecx+0x070],xmm3	/* a[jp+p12] <- ~t8 */\
			__asm	addpd	xmm5,xmm5	/* 2*t8 */\
			__asm	addpd	xmm4,xmm4	/* 2*t7 */\
			__asm	addpd	xmm5,xmm2	/*~t7 */\
			__asm	addpd	xmm4,xmm3	/*~t4 */\
			__asm	movaps	[ecx+0x060],xmm5	/* a[jt+p12] <- ~t7 */\
			__asm	movaps	[ecx+0x030],xmm4	/* a[jp+p4 ] <- ~t4 */\
		}

		/* Assumes addresses __add0, __add1 enter in eax, ebx, respectively: */
		#define SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(__c04,__c0C,__c14,__c1C, __tmp)\
		{\
			__asm	mov	ecx,__tmp\
			\
			/* Do the p00,p10 combo: */\
			__asm	mov	edx,__c04\
			\
			/* For interleaved [j1,j2] version, replace e.g.\
			\
				__asm	movaps	xmm0,[eax+0x40]	// a[jt+p8]\
				__asm	movaps	xmm1,[eax+0x50]	// a[jp+p8]\
				__asm	movaps	xmm2,[eax+0x40]	// xmm2 <- cpy a[jt+p8]\
				__asm	movaps	xmm3,[eax+0x50]	// xmm3 <- cpy a[jp+p8]\
			\
			by the following:\
			*/\
			/* Real parts: */\
			__asm	movaps		xmm6,[eax+0x40]	/* a[j1+p4 ], this is the scratch xmm register */\
			__asm	movaps		xmm0,[eax+0x40]	/* a[j1+p4 ], this is the active  xmm register */\
			__asm	unpckhpd	xmm6,[ebx+0x40]	/* a[j2+p4 ] gets read twice */\
			__asm	unpcklpd	xmm0,[ebx+0x40]	/* a[jt+p4 ] */\
			__asm	movaps	[ecx+0x200],xmm6	/* Store hi real in t00+32 */\
			/* Imag parts: */\
			__asm	movaps		xmm7,[eax+0x50]\
			__asm	movaps		xmm1,[eax+0x50]\
			__asm	unpckhpd	xmm7,[ebx+0x50]\
			__asm	unpcklpd	xmm1,[ebx+0x50]	/* a[jp+p4 ] */\
			__asm	movaps	[ecx+0x210],xmm7	/* Store hi imag in t01+32 */\
			\
			__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy a[jt+p4] */\
			__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy a[jp+p4] */\
			\
			/* From here on, things are identical to the code in radix32_dif_pass: */\
			__asm	mulpd	xmm0,[edx     ]	/* a[jt+p4]*c4 */\
			__asm	mulpd	xmm1,[edx     ]	/* a[jp+p4]*c4 */\
			__asm	mulpd	xmm2,[edx+0x10]	/* a[jt+p4]*s4 */\
			__asm	mulpd	xmm3,[edx+0x10]	/* a[jp+p4]*s4 */\
			__asm	addpd	xmm1,xmm2	/* xmm1 <- t01*/\
			__asm	subpd	xmm0,xmm3	/* xmm0 <- t00*/\
			__asm	movaps	xmm2,xmm0	/* xmm2 <- cpy t00*/\
			__asm	movaps	xmm3,xmm1	/* xmm3 <- cpy t01*/\
			\
			__asm	mov	edx,__c14\
			/* Real parts: */\
			__asm	movaps		xmm6,[eax+0x140]	/* a[j1+p14], this is the scratch xmm register  */\
			__asm	movaps		xmm4,[eax+0x140]	/* a[j1+p14], this is the active  xmm register */\
			__asm	unpckhpd	xmm6,[ebx+0x140]	/* a[j2+p14] gets read twice */\
			__asm	unpcklpd	xmm4,[ebx+0x140]	/* a[jt+p14] */\
			__asm	movaps	[ecx+0x220],xmm6	/* Store hi real in t11+32 */\
			/* Imag parts: */\
			__asm	movaps		xmm7,[eax+0x150]\
			__asm	movaps		xmm5,[eax+0x150]\
			__asm	unpckhpd	xmm7,[ebx+0x150]\
			__asm	unpcklpd	xmm5,[ebx+0x150]	/* a[jp+p14] */\
			__asm	movaps	[ecx+0x230],xmm7	/* Store hi imag in t12+32 */\
			\
			__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p14] */\
			__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p14] */\
			\
			__asm	mulpd	xmm4,[edx     ]	/* a[jt+p14]*c14 */\
			__asm	mulpd	xmm5,[edx     ]	/* a[jp+p14]*c14 */\
			__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p14]*s14 */\
			__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p14]*s14 */\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- it */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- rt		xmm6,7 free */\
			\
			__asm	addpd	xmm0,xmm4	/* ~t13<- t13+rt */\
			__asm	addpd	xmm1,xmm5	/* ~t14<- t14+it */\
			__asm	subpd	xmm2,xmm4	/* ~t15<- t13-rt */\
			__asm	subpd	xmm3,xmm5	/* ~t16<- t14-it	xmm4,5 free */\
			\
			/* Do the p08,p18 [hexadecimal here] combo - do p18 first so register assignments come out in same relative order as for p0,8 */\
			__asm	mov	edx,__c1C\
			/* Real parts: */\
			__asm	movaps		xmm6,[eax+0x1c0]	/* a[j1+p1C], this is the scratch xmm register  */\
			__asm	movaps		xmm4,[eax+0x1c0]	/* a[j1+p1C], this is the active  xmm register */\
			__asm	unpckhpd	xmm6,[ebx+0x1c0]	/* a[j2+p1C] gets read twice */\
			__asm	unpcklpd	xmm4,[ebx+0x1c0]	/* a[jt+p1C] */\
			__asm	movaps	[ecx+0x260],xmm6	/* Store hi real in t15+32 */\
			/* Imag parts: */\
			__asm	movaps		xmm7,[eax+0x1d0]\
			__asm	movaps		xmm5,[eax+0x1d0]\
			__asm	unpckhpd	xmm7,[ebx+0x1d0]\
			__asm	unpcklpd	xmm5,[ebx+0x1d0]	/* a[jp+p1C] */\
			__asm	movaps	[ecx+0x270],xmm7	/* Store hi imag in t16+32 */\
			\
			__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p1C] */\
			__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p1C] */\
			\
			__asm	mulpd	xmm4,[edx     ]	/* a[jt+p1C]*c1C */\
			__asm	mulpd	xmm5,[edx     ]	/* a[jp+p1C]*c1C */\
			__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p1C]*s1C */\
			__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p1C]*s1C */\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- it */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- rt		xmm6,7 free */\
			__asm	movaps	[ecx+0x010],xmm5	/* Store it */\
			__asm	movaps	[ecx      ],xmm4	/* Store rt */\
			\
			__asm	mov	edx,__c0C\
			/* Real parts: */\
			__asm	movaps		xmm6,[eax+0x0c0]	/* a[j1+p0C], this is the scratch xmm register  */\
			__asm	movaps		xmm4,[eax+0x0c0]	/* a[j1+p0C], this is the active  xmm register */\
			__asm	unpckhpd	xmm6,[ebx+0x0c0]	/* a[j2+p0C] gets read twice */\
			__asm	unpcklpd	xmm4,[ebx+0x0c0]	/* a[jt+p0C] */\
			__asm	movaps	[ecx+0x240],xmm6	/* Store hi real in t13+32 */\
			/* Imag parts: */\
			__asm	movaps		xmm7,[eax+0x0d0]\
			__asm	movaps		xmm5,[eax+0x0d0]\
			__asm	unpckhpd	xmm7,[ebx+0x0d0]\
			__asm	unpcklpd	xmm5,[ebx+0x0d0]	/* a[jp+p0C] */\
			__asm	movaps	[ecx+0x250],xmm7	/* Store hi imag in t14+32 */\
			\
			__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy a[jt+p0C] */\
			__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy a[jp+p0C] */\
			\
			__asm	mulpd	xmm4,[edx     ]	/* a[jt+p0C]*c0C*/\
			__asm	mulpd	xmm5,[edx     ]	/* a[jp+p0C]*c0C*/\
			__asm	mulpd	xmm6,[edx+0x10]	/* a[jt+p0C]*s0C*/\
			__asm	mulpd	xmm7,[edx+0x10]	/* a[jp+p0C]*s0C*/\
			__asm	addpd	xmm5,xmm6	/* xmm5 <- t14 */\
			__asm	subpd	xmm4,xmm7	/* xmm4 <- t13 */\
			__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy t14 */\
			__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy t13 */\
			\
			__asm	subpd	xmm4,[ecx      ]	/* ~t15<- t13-rt */\
			__asm	subpd	xmm5,[ecx+0x010]	/* ~t16<- t14-it */\
			__asm	addpd	xmm6,[ecx      ]	/* ~t13<- t13+rt */\
			__asm	addpd	xmm7,[ecx+0x010]	/* ~t14<- t14+it */\
			\
			/* Finish radix-4 butterfly and store results into temporary-array slots: */\
			/*\
			~t5 =t1 -t5;		~t1 =t1 +t5;\
			~t6 =t2 -t6;		~t2 =t2 +t6;\
			*/\
			__asm	subpd	xmm0,xmm6	/*~t5 */\
			__asm	subpd	xmm1,xmm7	/*~t6 */\
			__asm	movaps	[ecx+0x040],xmm0	/* a[jt+p8 ] <- ~t5 */\
			__asm	movaps	[ecx+0x050],xmm1	/* a[jp+p8 ] <- ~t6 */\
			__asm	addpd	xmm6,xmm6	/* 2*t5 */\
			__asm	addpd	xmm7,xmm7	/* 2*t6 */\
			__asm	addpd	xmm6,xmm0	/*~t1 */\
			__asm	addpd	xmm7,xmm1	/*~t2 */\
			__asm	movaps	[ecx      ],xmm6	/* a[jt    ] <- ~t1 */\
			__asm	movaps	[ecx+0x010],xmm7	/* a[jp    ] <- ~t2 */\
			\
			/*\
			~t7 =t3 +t8;		~t3 =t3 -t8;\
			~t8 =t4 -t7;		~t4 =t4 +t7;\
			*/\
			__asm	mov	edx, isrt2\
			\
			__asm	subpd	xmm2,xmm5	/*~t3 */\
			__asm	subpd	xmm3,xmm4	/*~t8 */\
			__asm	addpd	xmm5,xmm5	/* 2*t8 */\
			__asm	addpd	xmm4,xmm4	/* 2*t7 */\
			__asm	addpd	xmm5,xmm2	/*~t7 */\
			__asm	addpd	xmm4,xmm3	/*~t4 */\
			/*\
			t3 =(t3-t4)*ISRT2;	t4 =(t3+t4)*ISRT2;\
			t7 =(t7-t8)*ISRT2;	t8 =(t7+t8)*ISRT2;\
			*/\
			__asm	movaps	xmm6,xmm2	/* cpy t3 */\
			__asm	movaps	xmm7,xmm5	/* cpy t7 */\
			__asm	subpd	xmm2,xmm4	/* 3-4*/\
			__asm	subpd	xmm5,xmm3	/* 7-8*/\
			__asm	addpd	xmm6,xmm4	/* 3+4*/\
			__asm	addpd	xmm7,xmm3	/* 7+8*/\
			__asm	mulpd	xmm2,[edx]	/* (3-4)*ISRT2 */\
			__asm	mulpd	xmm5,[edx]	/* (7-8)*ISRT2 */\
			__asm	mulpd	xmm6,[edx]	/* (3+4)*ISRT2 */\
			__asm	mulpd	xmm7,[edx]	/* (7+8)*ISRT2 */\
			__asm	movaps	[ecx+0x020],xmm2	/* a[jt+p4 ] <- ~t3 */\
			__asm	movaps	[ecx+0x060],xmm5	/* a[jp+p12] <- ~t7 */\
			__asm	movaps	[ecx+0x030],xmm6	/* a[jp+p4 ] <- ~t4 */\
			__asm	movaps	[ecx+0x070],xmm7	/* a[jt+p12] <- ~t8 */\
		}

		/*
		SSE2-ified version of PAIR_SQUARE_4. Data enter in [__tAr, ~tAr], [__tAi, ~tAi], pairs, where the imaginary part
		of each input pair is assumed offset +0x10 in memory from the real part, i.e. needs no explicit pointer reference.

		For the sincos twiddles: using the notation of the scalar PAIR_SQUARE_4() macro,"__c" means [c0,s1], "__s" means [s0,c1].
		For these, due to the buterfly indexing pattern, we cannot assume that __s = __c + 0x10, so feed both pointers explicitly.

		We use shufpd xmm, xmm, 1 to swap lo and hi doubles of an xmm register for the various operations with one swapped input.
		*/
		#define PAIR_SQUARE_4_SSE2(__tAr, __tBr, __tCr, __tDr, __c, __s, __forth)\
		{\
		/*   calculate cross-product terms...\
			__rt=__tAr* ~tDr+__tAi* ~tDi; __rt=__rt+__rt;\
			__it=__tAi* ~tDr-__tAr* ~tDi; __it=__it+__it;\
		*/\
			__asm	mov	edx, __tDr\
			__asm	mov	eax, __tAr\
			\
			__asm	movaps	xmm6,[edx     ]	/* tDr */	\
			__asm	movaps	xmm7,[edx+0x10]	/* tDi */	\
			__asm	movaps	xmm0,[eax     ]	/* tAr */	\
			__asm	movaps	xmm3,[eax+0x10]	/* tAi */	\
			__asm	shufpd	xmm6,xmm6,1	/*~tDr */	\
			__asm	shufpd	xmm7,xmm7,1	/*~tDi */	\
			__asm	movaps	xmm2,[eax     ]	/* cpy tAr */	\
			__asm	movaps	xmm1,[eax+0x10]	/* cpy tAi */	\
			\
			__asm	mulpd	xmm0,xmm6	/* tAr*~tDr */	\
			__asm	mulpd	xmm3,xmm7	/* tAi*~tDi */	\
			__asm	mulpd	xmm1,xmm6	/* tAi*~tDr */	\
			__asm	mulpd	xmm2,xmm7	/* tAr*~tDi */	\
			__asm	addpd	xmm0,xmm3	/* rt */	\
			__asm	subpd	xmm1,xmm2	/* it */	\
			__asm	addpd	xmm0,xmm0	/* rt=rt+rt */	\
			__asm	addpd	xmm1,xmm1	/* it=it+it; xmm2-7 free */	\
			\
		/*\
			__st=__tBr* ~tCr+__tBi* ~tCi; __st=__st+__st;\
			__jt=__tBi* ~tCr-__tBr* ~tCi; __jt=__jt+__jt;\
		*/\
			__asm	mov	ecx, __tCr\
			__asm	mov	ebx, __tBr\
			\
			__asm	movaps	xmm6,[ecx     ]	/* tCr */	\
			__asm	movaps	xmm7,[ecx+0x10]	/* tCi */	\
			__asm	movaps	xmm2,[ebx     ]	/* tBr */	\
			__asm	movaps	xmm5,[ebx+0x10]	/* tBi */	\
			__asm	shufpd	xmm6,xmm6,1	/*~tCr */	\
			__asm	shufpd	xmm7,xmm7,1	/*~tCi */	\
			__asm	movaps	xmm4,[ebx     ]	/* cpy tBr */	\
			__asm	movaps	xmm3,[ebx+0x10]	/* cpy tBi */	\
			\
			__asm	mulpd	xmm2,xmm6	/* tBr*~tCr */	\
			__asm	mulpd	xmm5,xmm7	/* tBi*~tCi */	\
			__asm	mulpd	xmm3,xmm6	/* tBi*~tCr */	\
			__asm	mulpd	xmm4,xmm7	/* tBr*~tCi */	\
			__asm	addpd	xmm2,xmm5	/* st */	\
			__asm	subpd	xmm3,xmm4	/* jt */	\
			__asm	addpd	xmm2,xmm2	/* st=st+st */	\
			__asm	addpd	xmm3,xmm3	/* jt=jt+jt; xmm4-7 free */	\
		\
		/*   now calculate square terms and __store back in the same temporaries:	*/					\
		/*	__tmp=(__tAr+__tAi)*(__tAr-__tAi); __tAi=__tAr*__tAi; __tAi=__tAi+__tAi; __tAr=__tmp;	*/	\
		\
			__asm	movaps	xmm4,[eax     ]	/* __tAr */	\
			__asm	movaps	xmm5,[eax+0x10]	/* __tAi */	\
			__asm	subpd	xmm4,xmm5		/* (__tAr-__tAi) */	\
			__asm	addpd	xmm5,xmm5		/*      2*__tAi  */	\
			__asm	addpd	xmm5,xmm4		/* (__tAr+__tAi) */	\
			__asm	mulpd	xmm4,xmm5		/*>__tAr */	\
			\
			__asm	movaps	xmm5,[eax     ]	/* __tAr */	\
			__asm	mulpd	xmm5,[eax+0x10]	/* __tAr*__tAi */	\
			__asm	addpd	xmm5,xmm5		/*>__tAi */	\
			__asm	movaps	[eax     ],xmm4	/* tmp store >__tAr */	\
			__asm	movaps	[eax+0x10],xmm5	/* tmp store >__tAi */	\
			\
			__asm	subpd	xmm0,xmm4	/* rt-__tAr */	\
			__asm	subpd	xmm1,xmm5	/* it-__tAi; xmm4-7 free */	\
		\
		/*	__tmp=(__tBr+__tBi)*(__tBr-__tBi); __tBi=__tBr*__tBi; __tBi=__tBi+__tBi; __tBr=__tmp;	*/	\
		/*** [Can be done in parallel with above segment] ***/											\
		\
			__asm	movaps	xmm6,[ebx     ]	/* __tBr */	\
			__asm	movaps	xmm7,[ebx+0x10]	/* __tBi */	\
			__asm	subpd	xmm6,xmm7		/* (__tBr-__tBi) */	\
			__asm	addpd	xmm7,xmm7		/*      2*__tBi  */	\
			__asm	addpd	xmm7,xmm6		/* (__tBr+__tBi) */	\
			__asm	mulpd	xmm6,xmm7		/*>__tBr */	\
			\
			__asm	movaps	xmm7,[ebx     ]	/* __tBr */	\
			__asm	mulpd	xmm7,[ebx+0x10]	/* __tBr*__tBi */	\
			__asm	addpd	xmm7,xmm7		/*>__tBi */	\
			__asm	movaps	[ebx     ],xmm6	/* tmp store >__tBr */	\
			__asm	movaps	[ebx+0x10],xmm7	/* tmp store >__tBi */	\
			\
			__asm	subpd	xmm2,xmm6	/* st-__tBr */	\
			__asm	subpd	xmm3,xmm7	/* jt-__tBi; xmm4-7 free */	\
		\
		/*	__tmp=(__tDr+__tDi)*(__tDr-__tDi); __tDi=__tDr*__tDi; __tDi=__tDi+__tDi; __tDr=__tmp;	*/	\
		\
			__asm	movaps	xmm4,[edx     ]	/* __tDr */	\
			__asm	movaps	xmm5,[edx+0x10]	/* __tDi */	\
			__asm	subpd	xmm4,xmm5		/* (__tDr-__tDi) */	\
			__asm	addpd	xmm5,xmm5		/*      2*__tDi  */	\
			__asm	addpd	xmm5,xmm4		/* (__tDr+__tDi) */	\
			__asm	mulpd	xmm4,xmm5		/*>__tDr */	\
			\
			__asm	movaps	xmm5,[edx     ]	/* __tDr */	\
			__asm	mulpd	xmm5,[edx+0x10]	/* __tDr*__tDi */	\
			__asm	addpd	xmm5,xmm5		/*>__tDi */	\
			__asm	movaps	[edx     ],xmm4	/* tmp store ~tDr */	\
			__asm	movaps	[edx+0x10],xmm5	/* tmp store ~tDi */	\
			__asm	shufpd	xmm4,xmm4,1	/*~tDr */	\
			__asm	shufpd	xmm5,xmm5,1	/*~tDi */	\
			\
			__asm	subpd	xmm0,xmm4	/* rt-__tAr- ~tDr */	\
			__asm	addpd	xmm1,xmm5	/* it-__tAi+ ~tDi; xmm4-7 free */	\
		\
		/*	__tmp=(__tCr+__tCi)*(__tCr-__tCi); __tCi=__tCr*__tCi; __tCi=__tCi+__tCi; __tCr=__tmp;	*/	\
		/*** [Can be done in parallel with above segment] ***/											\
		\
			__asm	movaps	xmm6,[ecx     ]	/* __tCr */	\
			__asm	movaps	xmm7,[ecx+0x10]	/* __tCi */	\
			__asm	subpd	xmm6,xmm7		/* (__tCr-__tCi) */	\
			__asm	addpd	xmm7,xmm7		/*      2*__tCi  */	\
			__asm	addpd	xmm7,xmm6		/* (__tCr+__tCi) */	\
			__asm	mulpd	xmm6,xmm7		/*>__tCr */	\
			\
			__asm	movaps	xmm7,[ecx     ]	/* __tCr */	\
			__asm	mulpd	xmm7,[ecx+0x10]	/* __tCr*__tCi */	\
			__asm	addpd	xmm7,xmm7		/*>__tCi */	\
			__asm	movaps	[ecx     ],xmm6	/* tmp store ~tCr */	\
			__asm	movaps	[ecx+0x10],xmm7	/* tmp store ~tCi */	\
			__asm	shufpd	xmm6,xmm6,1	/*~tCr */	\
			__asm	shufpd	xmm7,xmm7,1	/*~tCi */	\
			\
			__asm	subpd	xmm2,xmm6	/* st-__tBr- ~tCr */	\
			__asm	addpd	xmm3,xmm7	/* jt-__tBi+ ~tCi; xmm4-7 free */	\
		\
		/*	__tmp=((1.0+__c)*__rt-__s*__it)*0.25;					\
			__it =((1.0+__c)*__it+__s*__rt)*0.25;	__rt=__tmp;	*/	\
		/*** [Can be done in parallel with above segment] ***/		\
		\
			__asm	mov	eax, __c	\
			__asm	mov	ebx, __s	\
			__asm	mov	edx, __forth	\
			__asm	movaps	xmm4,xmm0	/* cpy rt */	\
			__asm	movaps	xmm5,xmm1	/* cpy it */	\
			__asm	mulpd	xmm0,[eax]	/* c*rt */	\
			__asm	mulpd	xmm1,[eax]	/* c*it */	\
			__asm	addpd	xmm0,xmm4	/* (c+1.0)*rt */	\
			__asm	addpd	xmm1,xmm5	/* (c+1.0)*it */	\
			__asm	mulpd	xmm4,[ebx]	/* s*rt */	\
			__asm	mulpd	xmm5,[ebx]	/* s*it */	\
			__asm	subpd	xmm0,xmm5	/* (c+1.0)*rt-s*it */				\
			__asm	addpd	xmm1,xmm4	/* (c+1.0)*it+s*rt; xmm4,5 free */	\
			__asm	mulpd	xmm0,[edx]	/* -rt Both of these inherit the sign flip [w.r.to the non-SSE2 PAIR_SQUARE_4 macro] */	\
			__asm	mulpd	xmm1,[edx]	/* -it that resulted from the in-place-friendlier (rt-__tAr- ~tDr) reordering above. */	\
		\
		/*	__tmp=((1.0-__s)*__st-__c*__jt)*0.25;					\
			__jt =((1.0-__s)*__jt+__c*__st)*0.25	__st=__tmp;	*/	\
		/*** [Can be done in parallel wjth above segment] ***/		\
		\
			__asm	movaps	xmm6,xmm2	/* cpy st */	\
			__asm	movaps	xmm7,xmm3	/* cpy jt */	\
			__asm	mulpd	xmm2,[ebx]	/* s*st */	\
			__asm	mulpd	xmm3,[ebx]	/* s*jt */	\
			__asm	subpd	xmm2,xmm6	/* (s-1.0)*st, note sign flip! */	\
			__asm	subpd	xmm3,xmm7	/* (s-1.0)*jt, note sign flip! */	\
			__asm	mulpd	xmm6,[eax]	/* c*st */	\
			__asm	mulpd	xmm7,[eax]	/* c*jt */	\
			__asm	addpd	xmm2,xmm7	/* -[(1.0-s)*st-c*jt] */				\
			__asm	subpd	xmm3,xmm6	/* -[(1.0-s)*jt+c*st]; xmm6,7 free */	\
			__asm	mulpd	xmm2,[edx]	/* +st Sign flip due to (s-1.0) reordering here */	\
			__asm	mulpd	xmm3,[edx]	/* +jt cancels earlier one due to in-place-friendlier (st-__tBr- ~tCr) reordering above. */	\
		\
		/*...and now complete and store the results. We flip the signs on st and jt here to undo the above -st,-jt negations. */\
		/*	__tAr = (__tAr+__rt);	\
			__tAi = (__tAi+__it);	\
			__tBr = (__tBr-__st);	\
			__tBi = (__tBi-__jt); */\
		\
			__asm	mov	eax, __tAr	\
			__asm	mov	ebx, __tBr	\
			\
			__asm	movaps	xmm4,[eax     ]	/* __tAr */	\
			__asm	movaps	xmm5,[eax+0x10]	/* __tAi */	\
			__asm	movaps	xmm6,[ebx     ]	/* __tBr */	\
			__asm	movaps	xmm7,[ebx+0x10]	/* __tBi */	\
			__asm	addpd	xmm4,xmm0		/* (__tAr+__rt) */	\
			__asm	addpd	xmm5,xmm1		/* (__tAi+__it) */	\
			__asm	subpd	xmm6,xmm2		/* (__tBr-__st) */	\
			__asm	subpd	xmm7,xmm3		/* (__tBi-__jt) */	\
			__asm	movaps	[eax     ],xmm4	/* store >__tAr */	\
			__asm	movaps	[eax+0x10],xmm5	/* store >__tAi */	\
			__asm	movaps	[ebx     ],xmm6	/* store >__tBr */	\
			__asm	movaps	[ebx+0x10],xmm7	/* store >__tBi */	\
			\
		/*...N-j terms are as above, but with the replacements: __tAr<--> ~tDr, __tAi<--> ~tDi, __it|-->-__it. */\
		/*	__tDr = (__tDr+ ~rt);	\
			__tDi = (__tDi- ~it);	\
			__tCr = (__tCr- ~st);	\
			__tCi = (__tCi+ ~jt); */\
			\
			__asm	mov	ecx, __tCr	\
			__asm	mov	edx, __tDr	\
			\
			__asm	shufpd	xmm0,xmm0,1		/* ~rt */	\
			__asm	shufpd	xmm1,xmm1,1		/* ~it */	\
			__asm	shufpd	xmm2,xmm2,1		/* ~st */	\
			__asm	shufpd	xmm3,xmm3,1		/* ~jt */	\
			\
			__asm	movaps	xmm4,[edx     ]	/* __tDr */	\
			__asm	movaps	xmm5,[edx+0x10]	/* __tDi */	\
			__asm	movaps	xmm6,[ecx     ]	/* __tCr */	\
			__asm	movaps	xmm7,[ecx+0x10]	/* __tCi */	\
			__asm	addpd	xmm4,xmm0		/* (__tDr+ ~rt) */	\
			__asm	subpd	xmm5,xmm1		/* (__tDi- ~it) */	\
			__asm	subpd	xmm6,xmm2		/* (__tCr- ~st) */	\
			__asm	addpd	xmm7,xmm3		/* (__tCi+ ~jt) */	\
			__asm	movaps	[edx     ],xmm4	/* store >__tDr */	\
			__asm	movaps	[edx+0x10],xmm5	/* store >__tDi */	\
			__asm	movaps	[ecx     ],xmm6	/* store >__tCr */	\
			__asm	movaps	[ecx+0x10],xmm7	/* store >__tCi */	\
			\
			/**************************************************************/\
			/* Total Cost:  52 MOVapd,  48 ADD/SUBpd, 28 MULpd, 12 SHUFpd */\
			/******** [~1/4 the cost of a radix-16 DIF or DIT pass] *******/\
		}

	/******************************************************************************************************************************************/

		/*...Radix-5 DFT: Inputs enter in memory locations __i0,__i1,__i2,__i3,__i4, assumed disjoint with inputs: */\
		#define SSE2_RADIX_05_DFT_0TWIDDLE(__i0,__i1,__i2,__i3,__i4, __c1, __o0,__o1,__o2,__o3,__o4)\
		{\
			__asm	mov	eax, __i1\
			__asm	mov	ebx, __i2\
			__asm	mov	ecx, __i3\
			__asm	mov	edx, __i4\
			__asm	mov	esi, __i0\
			__asm	mov	edi, __o0\
			\
			__asm	movaps	xmm0,[eax     ]	/* Ar1 */\
			__asm	movaps	xmm1,[eax+0x10]	/* Ai1 */\
			__asm	movaps	xmm2,[ebx     ]	/* Ar2 */\
			__asm	movaps	xmm3,[ebx+0x10]	/* Ai2 */\
			__asm	movaps	xmm4,[ecx     ]	/* Ar3 */\
			__asm	movaps	xmm5,[ecx+0x10]	/* Ai3 */\
			__asm	movaps	xmm6,[edx     ]	/* Ar4 */\
			__asm	movaps	xmm7,[edx+0x10]	/* Ai4 */\
		/*\
			Br4 = Ar1 - Ar4;		Bi4 = Ai1 - Ai4;			\
			Br1 = Ar1 + Ar4;		Bi1 = Ai1 + Ai4;			\
		*/\
			__asm	subpd	xmm0,xmm6		/* Br4 = Ar1 - Ar4 */\
			__asm	subpd	xmm1,xmm7		/* Bi4 = Ai1 - Ai4 */\
			__asm	addpd	xmm6,xmm6		/*           2*Ar4 */\
			__asm	addpd	xmm7,xmm7		/*           2*Ai4 */\
			__asm	addpd	xmm6,xmm0		/* Br1 = Ar1 + Ar4 */\
			__asm	addpd	xmm7,xmm1		/* Bi1 = Ai1 + Ai4 */\
		/*\
			Br3 = Ar2 - Ar3;		Bi3 = Ai2 - Ai3;			\
			Br2 = Ar2 + Ar3;		Bi2 = Ai2 + Ai3;			\
		*/\
			__asm	subpd	xmm2,xmm4		/* Br3 = Ar2 - Ar3 */\
			__asm	subpd	xmm3,xmm5		/* Bi3 = Ai2 - Ai3 */\
			__asm	addpd	xmm4,xmm4		/*           2*Ar3 */\
			__asm	addpd	xmm5,xmm5		/*           2*Ai3 */\
			__asm	addpd	xmm4,xmm2		/* Br2 = Ar2 + Ar3 */\
			__asm	addpd	xmm5,xmm3		/* Bi2 = Ai2 + Ai3 */\
		/*
			rt  = Br1 + Br2;			it  = Bi1 + Bi2;		\
			Br0 = Ar0 + rt;			Bi0 = Ai0 + it;		// y0\
			rt  = Br0 + rt  *cc1;		it  = Bi0 + it  *cc1;	\
			Br2 =(Br1 - Br2)*cc2;		Bi2 =(Bi1 - Bi2)*cc2;	\
		*/\
			__asm	mov	eax, __c1\
			__asm	subpd	xmm6,xmm4		/*~Br2 = Br1 - Br2 */\
			__asm	subpd	xmm7,xmm5		/*~Bi2 = Bi1 - Bi2 */\
			__asm	addpd	xmm4,xmm4		/*           2*Br2 */\
			__asm	addpd	xmm5,xmm5		/*           2*Bi2 */\
			__asm	addpd	xmm4,xmm6		/*  rt = Br1 + Br2 */\
			__asm	addpd	xmm5,xmm7		/*  it = Bi1 + Bi2 */\
			__asm	addpd	xmm4,[esi     ]	/* Br0,   rt lost */\
			__asm	addpd	xmm5,[esi+0x10]	/* Bi0,   it lost */\
			__asm	movaps	[edi     ],xmm4	/* Write Br0 */\
			__asm	movaps	[edi+0x10],xmm5	/* Write Bi0 */\
			__asm	mulpd	xmm6,[eax+0x10]	/* Br2 =(Br1 - Br2)*cc2 */\
			__asm	mulpd	xmm7,[eax+0x10]	/* Bi2 =(Bi1 - Bi2)*cc2 */\
			__asm	subpd	xmm4,[esi     ]	/* Restore rt */\
			__asm	subpd	xmm5,[esi+0x10]	/* Restore it */\
			__asm	mulpd	xmm4,[eax     ]	/* rt *= cc1 */\
			__asm	mulpd	xmm5,[eax     ]	/* it *= cc1 */\
			__asm	addpd	xmm4,[edi     ]	/* rt  = Br0 + rt*cc1*/\
			__asm	addpd	xmm5,[edi+0x10]	/* it  = Bi0 + it*cc1 */\
		/*\
			Br1 = rt  + Br2;		Bi1 = it  + Bi2;			\
			Br2 = rt  - Br2;		Bi2 = it  - Bi2;			\
		*/\
			__asm	subpd	xmm4,xmm6		/*~Br2 =  rt - Br2 */\
			__asm	subpd	xmm5,xmm7		/*~Bi2 =  it - Bi2 */\
			__asm	addpd	xmm6,xmm6		/*           2*Br2 */\
			__asm	addpd	xmm7,xmm7		/*           2*Bi2 */\
			__asm	addpd	xmm6,xmm4		/*~Br1 =  rt + Br2 */\
			__asm	addpd	xmm7,xmm5		/*~Bi1 =  it + Bi2 */\
			__asm	movaps	[esi     ],xmm4	/* tmpstr Br2 */\
			__asm	movaps	[esi+0x10],xmm5	/* tmpstr Bi2, xmm4,5 FREE */\
		/*\
			rt  =(Br4 - Br3)* s2;	it  =(Bi4 - Bi3)* s2;	\
		*/\
			__asm	movaps	xmm4,xmm0	/* copy Br4 */\
			__asm	movaps	xmm5,xmm1	/* copy Bi4 */\
			__asm	subpd	xmm0,xmm2		/*~Br4 = Br4 - Br3 */\
			__asm	subpd	xmm1,xmm3		/*~Bi4 = Bi4 - Bi3 */\
			__asm	mulpd	xmm0,[eax+0x20]	/*  rt =(Br4 - Br3)* s2 */\
			__asm	mulpd	xmm1,[eax+0x20]	/*  it =(Bi4 - Bi3)* s2 */\
			__asm	mulpd	xmm2,[eax+0x30]	/*~Br3 =       Br3 *ss1 */\
			__asm	mulpd	xmm3,[eax+0x30]	/*~Bi3 =       Bi3 *ss1 */\
			__asm	mulpd	xmm4,[eax+0x40]	/*~Br4 =       Br4 *ss2 */\
			__asm	mulpd	xmm5,[eax+0x40]	/*~Bi4 =       Bi4 *ss2 */\
		/*\
			Br3 = rt  + Br3 *ss1;	Bi3 = it  + Bi3 *ss1;	\
			Br4 = rt  - Br4 *ss2;	Bi4 = it  - Bi4 *ss2;	\
		*/\
			__asm	addpd	xmm2,xmm0		/* Br3 = rt  + Br3 *ss1 */\
			__asm	addpd	xmm3,xmm1		/* Bi3 = it  + Bi3 *ss1 */\
			__asm	subpd	xmm0,xmm4		/* Br4 = rt  - Br4 *ss2 */\
			__asm	subpd	xmm1,xmm5		/* Bi4 = it  - Bi4 *ss2, xmm4,5 FREE */\
			__asm	movaps	xmm4,[esi     ]	/* Br2 */\
			__asm	movaps	xmm5,[esi+0x10]	/* Bi2 */\
		/*\
			y4r = Br1 + Bi3;			y4i = Bi1 - Br3;\
			y1r = Br1 - Bi3;			y1i = Bi1 + Br3;\
		*/\
			__asm	mov	eax, __o1\
			__asm	mov	edx, __o4\
			__asm	subpd	xmm6,xmm3		/* y1r = Br1 - Bi3 */\
			__asm	subpd	xmm7,xmm2		/* y4i = Bi1 - Br3 */\
			__asm	addpd	xmm3,xmm3		/*           2*Bi3 */\
			__asm	addpd	xmm2,xmm2		/*           2*Br3 */\
			__asm	movaps	[eax     ],xmm6	/* Write y1r */\
			__asm	movaps	[edx+0x10],xmm7	/* Write y4i */\
			__asm	addpd	xmm3,xmm6		/* y4r = Br1 + Bi3 */\
			__asm	addpd	xmm2,xmm7		/* y1i = Bi1 + Br3 */\
			__asm	movaps	[edx     ],xmm3	/* Write y4r */\
			__asm	movaps	[eax+0x10],xmm2	/* Write y1i */\
		/*\
			y3r = Br2 + Bi4;			y3i = Bi2 - Br4;\
			y2r = Br2 - Bi4;			y2i = Bi2 + Br4;\
		*/\
			__asm	mov	ebx, __o2\
			__asm	mov	ecx, __o3\
			__asm	subpd	xmm4,xmm1		/* y2r = Br2 - Bi4 */\
			__asm	subpd	xmm5,xmm0		/* y3i = Bi2 - Br4 */\
			__asm	addpd	xmm1,xmm1		/*           2*Bi4 */\
			__asm	addpd	xmm0,xmm0		/*           2*Br4 */\
			__asm	movaps	[ebx     ],xmm4	/* Write y2r */\
			__asm	movaps	[ecx+0x10],xmm5	/* Write y3i */\
			__asm	addpd	xmm1,xmm4		/* y3r = Br2 + Bi4 */\
			__asm	addpd	xmm0,xmm5		/* y2i = Bi2 + Br4 */\
			__asm	movaps	[ecx     ],xmm1	/* Write y3r */\
			__asm	movaps	[ebx+0x10],xmm0	/* Write y2i */\
		}

	/******************************************************************************************************************************************/

		/*...Radix-7 DFT: Inputs in memory locations __i0-6, outputs go into memory locations __o0-6, possibly coincident with inputs:\ */\
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

	/******************************************************************************************************************************************/

	#elif defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC)

		#if OS_BITS == 32

			#include "sse2_macro_gcc32.h"

		#else

			#include "sse2_macro_gcc64.h"

		#endif

	#else

		#error sse2_macro.h: No implementation for SSE2 macros on this platform!

	#endif

#endif

#endif	/* #ifndef sse2_macro_included */
