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
#ifndef radix09_sse_macro_h_included
#define radix09_sse_macro_h_included

	#include "sse2_macro.h"

	#if defined(COMPILER_TYPE_MSVC)
	#error Unsupported - Need to modify maccro to base all sincos addresses off of c1 [c2 = c1+0x20, c3m1 = c1+0x40, c4 = c1+0x60]
		/*...Radix-9 DIF: Outputs [t's] are temporaries assumed to be adjacent in memory, i.e. reals separated by 0x20 bytes.\
				Inputs enter in memory locations __i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8, assumed disjoint with outputs:\
		*/\
		#define SSE2_RADIX_09_DIF(__i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8, __c1,__c2,__c3m1,__c4, __t00)\
		{\
		/*...gather the needed data (9 64-b__it complex, i.e. 18 64-b__it reals) and do three radix-3 transforms: */\
		/*\
			t00 =A0r;				t01 =A0i;\
			t02 =A3r+A6r;			t03 =A3i+A6i;\
			t04 =A3r-A6r;			t05 =A3i-A6i;\
			t00 =t00+t02;			t01 =t01+t03;\
			t02 =t00+c3m1*t02;		t03 =t01+c3m1*t03;\
			rt  =s3*t04;			it  =s3*t05;\
			t04 =t02+it;			t05 =t03-rt;\
			t02 =t02-it;			t03 =t03+rt;\
		*/\
			__asm	mov	eax, __i0\
			__asm	mov	ebx, __i3\
			__asm	mov	ecx, __i6\
			__asm	mov	edx, __c3m1\
			\
			__asm	movaps	xmm2,[ebx     ]	/* A3r */\
			__asm	movaps	xmm3,[ebx+0x10]	/* A3i */\
			__asm	movaps	xmm0,[eax     ]	/* A0r */\
			__asm	movaps	xmm1,[eax+0x10]	/* A0i */\
			__asm	movaps	xmm6,[ecx     ]	/* A6r */\
			__asm	movaps	xmm7,[ecx+0x10]	/* A6i */\
			__asm	movaps	xmm4,xmm2		/* t04: init to copy of A3r */\
			__asm	movaps	xmm5,xmm3		/* t05: init to copy of A3i */\
			\
		__asm	mov	esi, __t00\
		__asm	mov	edi, esi\
		__asm	add	edi, 0x20\
			\
			__asm	addpd	xmm2,xmm6		/* t02=A3r+A6r */\
			__asm	addpd	xmm3,xmm7		/* t03=A3i+A6i */\
			__asm	subpd	xmm4,xmm6		/* t04=A3r-A6r */\
			__asm	subpd	xmm5,xmm7		/* t05=A3i-A6i */\
			__asm	addpd	xmm0,xmm2	/* t00 = t00+t02 */\
			__asm	addpd	xmm1,xmm3	/* t01 = t01+t03 */\
			__asm	movaps	xmm6,[edx     ]	/* c3m1*/\
			__asm	movaps	xmm7,[edx+0x10]	/* s3  */\
			__asm	movaps	[esi      ],xmm0	/* <- t00 */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t01 */\
			\
		__asm	mov	edx, esi\
		__asm	add	edx, 0x40\
			\
			__asm	mulpd	xmm2,xmm6		/* t02 *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t03 *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = t04*s3   */\
			__asm	mulpd	xmm5,xmm7		/* it = t05*s3   */\
			__asm	addpd	xmm2,xmm0	/* t02 = t00+c3m1*t02 */\
			__asm	addpd	xmm3,xmm1	/* t03 = t01+c3m1*t03 */\
			\
			__asm	movaps	xmm0,xmm2	/* copy t02 */\
			__asm	movaps	xmm1,xmm3	/* copy t03 */\
			\
			__asm	subpd	xmm2,xmm5	/* t02 = t02-it */\
			__asm	addpd	xmm3,xmm4	/* t03 = t03+rt */\
			__asm	addpd	xmm0,xmm5	/* t04 = t02+it */\
			__asm	subpd	xmm1,xmm4	/* t05 = t03-rt */\
			\
			__asm	movaps	[edi      ],xmm2	/* <- t02 */\
			__asm	movaps	[edi+0x010],xmm3	/* <- t03 */\
			__asm	movaps	[edx      ],xmm0	/* <- t04 */\
			__asm	movaps	[edx+0x010],xmm1	/* <- t05 */\
			\
		/*\
			t06 =A1r;				t07 =A1i;\
			t08 =A4r+A7r;			t09 =A4i+A7i;\
			t0a =A4r-A7r;			t0b =A4i-A7i;\
			t06 =t06+t08;			t07 =t07+t09;\
			t08 =t06+c3m1*t08;		t09 =t07+c3m1*t09;\
			rt  =s3*t0a;			it  =s3*t0b;\
			t0a =t08+it;			t0b =t09-rt;\
			t08 =t08-it;			t09 =t09+rt;\
		*/\
			__asm	mov	eax, __i1\
			__asm	mov	ebx, __i4\
			__asm	mov	ecx, __i7\
			\
			__asm	movaps	xmm4,[ebx     ]	/* A4r */\
			__asm	movaps	xmm5,[ebx+0x10]	/* A4i */\
			__asm	movaps	xmm0,[eax     ]	/* A1r */\
			__asm	movaps	xmm1,[eax+0x10]	/* A1i */\
			__asm	movaps	xmm2,[ecx     ]	/* A7r */\
			__asm	movaps	xmm3,[ecx+0x10]	/* A7i */\
			\
		__asm	add	esi, 0x60\
			\
			__asm	subpd	xmm4,xmm2		/* t0a=A4r-A7r */\
			__asm	subpd	xmm5,xmm3		/* t0b=A4i-A7i */\
			__asm	addpd	xmm2,xmm2		/*       2*A7r */\
			__asm	addpd	xmm3,xmm3		/*       2*A7i */\
			__asm	addpd	xmm2,xmm4		/* t08=A4r+A7r */\
			__asm	addpd	xmm3,xmm5		/* t09=A4i+A7i */\
			__asm	addpd	xmm0,xmm2	/* t06 = t06+t08 */\
			__asm	addpd	xmm1,xmm3	/* t07 = t07+t09 */\
			__asm	movaps	[esi      ],xmm0	/* <- t06 */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t07 */\
			\
		__asm	add	edi, 0x60\
		__asm	add	edx, 0x60\
			\
			__asm	mulpd	xmm2,xmm6		/* t08 *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t09 *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = t0a*s3   */\
			__asm	mulpd	xmm5,xmm7		/* it = t0b*s3   */\
			__asm	addpd	xmm2,xmm0	/* t08 = t06+c3m1*t08 */\
			__asm	addpd	xmm3,xmm1	/* t09 = t07+c3m1*t09 */\
			\
			__asm	movaps	xmm0,xmm2	/* copy t08 */\
			__asm	movaps	xmm1,xmm3	/* copy t09 */\
			\
			__asm	subpd	xmm2,xmm5	/* t08 = t08-it */\
			__asm	addpd	xmm3,xmm4	/* t09 = t09+rt */\
			__asm	addpd	xmm0,xmm5	/* t0a = t08+it */\
			__asm	subpd	xmm1,xmm4	/* t0b = t09-rt */\
			\
			__asm	movaps	[edi      ],xmm2	/* <- t08 */\
			__asm	movaps	[edi+0x010],xmm3	/* <- t09 */\
			__asm	movaps	[edx      ],xmm0	/* <- t0a */\
			__asm	movaps	[edx+0x010],xmm1	/* <- t0b */\
			\
		/*\
			t0c =A2r;				t0d =A2i;\
			t0e =A5r+A8r;			t0f =A5i+A8i;\
			t0g =A5r-A8r;			t0h =A5i-A8i;\
			t0c =t0c+t0e;			t0d =t0d+t0f;\
			t0e =t0c+c3m1*t0e;		t0f =t0d+c3m1*t0f;\
			rt  =s3*t0g;			it  =s3*t0h;\
			t0g =t0e+it;			t0h =t0f-rt;\
			t0e =t0e-it;			t0f =t0f+rt;\
		*/\
			__asm	mov	eax, __i2\
			__asm	mov	ebx, __i5\
			__asm	mov	ecx, __i8\
			\
			__asm	movaps	xmm4,[ebx     ]	/* A5r */\
			__asm	movaps	xmm5,[ebx+0x10]	/* A5i */\
			__asm	movaps	xmm0,[eax     ]	/* A2r */\
			__asm	movaps	xmm1,[eax+0x10]	/* A2i */\
			__asm	movaps	xmm2,[ecx     ]	/* A8r */\
			__asm	movaps	xmm3,[ecx+0x10]	/* A8i */\
			\
		__asm	add	esi, 0x60\
			\
			__asm	subpd	xmm4,xmm2		/* t0g=A5r-A8r */\
			__asm	subpd	xmm5,xmm3		/* t0h=A5i-A8i */\
			__asm	addpd	xmm2,xmm2		/*       2*A8r */\
			__asm	addpd	xmm3,xmm3		/*       2*A8i */\
			__asm	addpd	xmm2,xmm4		/* t0e=A5r+A8r */\
			__asm	addpd	xmm3,xmm5		/* t0f=A5i+A8i */\
			__asm	addpd	xmm0,xmm2	/* t0c = t0c+t0e */\
			__asm	addpd	xmm1,xmm3	/* t0d = t0d+t0f */\
			__asm	movaps	[esi      ],xmm0	/* <- t0c */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t0d */\
			\
		__asm	add	edi, 0x60\
		__asm	add	edx, 0x60\
			\
			__asm	mulpd	xmm2,xmm6		/* t0e *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t0f *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = t0g*s3   */\
			__asm	mulpd	xmm5,xmm7		/* it = t0h*s3   */\
			__asm	addpd	xmm2,xmm0	/* t0e = t0c+c3m1*t0e */\
			__asm	addpd	xmm3,xmm1	/* t0f = t0d+c3m1*t0f */\
			\
			__asm	movaps	xmm0,xmm2	/* copy t0e */\
			__asm	movaps	xmm1,xmm3	/* copy t0f */\
			\
			__asm	subpd	xmm2,xmm5	/* t0e = t0e-it */\
			__asm	addpd	xmm3,xmm4	/* t0f = t0f+rt */\
			__asm	addpd	xmm0,xmm5	/* t0g = t0e+it */\
			__asm	subpd	xmm1,xmm4	/* t0h = t0f-rt */\
			\
			__asm	movaps	[edi      ],xmm2	/* <- t0e */\
			__asm	movaps	[edi+0x010],xmm3	/* <- t0f */\
			__asm	movaps	[edx      ],xmm0	/* <- t0g */\
			__asm	movaps	[edx+0x010],xmm1	/* <- t0h */\
			\
		/******************************************************************************/\
		/*...and now do three more radix-3 transforms, including the twiddle factors: */\
		/******************************************************************************/\
		/*\
			rt  =t06;				it  =t07;\
			t06 =rt+t0c;			t07 =it+t0d;\
			t0c =rt-t0c;			t0d =it-t0d;\
			t00 =t00+t06;			t01 =t01+t07;\
			t06 =t00+c3m1*t06;		t07 =t01+c3m1*t07;\
			rt  =s3*t0c;			it  =s3*t0d;\
			t0c =t06+it;			t0d =t07-rt;\
			t06 =t06-it;			t07 =t07+rt;\
		*/\
			__asm	mov	eax, __t00	/* t00 */\
			__asm	mov	ebx, eax\
			__asm	mov	ecx, eax\
			__asm	add	ebx, 0x60	/* t06 */\
			__asm	add	ecx, 0xc0	/* t0c */\
			\
			__asm	movaps	xmm4,[ebx     ]	/* t06 */\
			__asm	movaps	xmm5,[ebx+0x10]	/* t07 */\
			__asm	movaps	xmm2,[ecx     ]	/* t0c */\
			__asm	movaps	xmm3,[ecx+0x10]	/* t0d */\
			__asm	movaps	xmm0,[eax     ]	/* t00 */\
			__asm	movaps	xmm1,[eax+0x10]	/* t01 */\
			\
			__asm	subpd	xmm4,xmm2		/* t0c=t06-t0c */\
			__asm	subpd	xmm5,xmm3		/* t0d=t07-t0d */\
			__asm	addpd	xmm2,xmm2		/*       2*t0c */\
			__asm	addpd	xmm3,xmm3		/*       2*t0d */\
			__asm	addpd	xmm2,xmm4		/* t06=t06+t0c */\
			__asm	addpd	xmm3,xmm5		/* t07=t07+t0d */\
			__asm	addpd	xmm0,xmm2	/* t00 = t00+t06 */\
			__asm	addpd	xmm1,xmm3	/* t01 = t01+t07 */\
			__asm	movaps	[eax      ],xmm0	/* <- t00 */\
			__asm	movaps	[eax+0x010],xmm1	/* <- t01 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t06 *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t07 *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = t0c*s3   */\
			__asm	mulpd	xmm5,xmm7		/* it = t0d*s3   */\
			__asm	addpd	xmm2,xmm0	/* t06 = t00+c3m1*t06 */\
			__asm	addpd	xmm3,xmm1	/* t07 = t01+c3m1*t07 */\
			\
			__asm	movaps	xmm0,xmm2	/* copy t06 */\
			__asm	movaps	xmm1,xmm3	/* copy t07 */\
			\
			__asm	subpd	xmm2,xmm5	/* t06 = t06-it */\
			__asm	addpd	xmm3,xmm4	/* t07 = t07+rt */\
			__asm	addpd	xmm0,xmm5	/* t0c = t06+it */\
			__asm	subpd	xmm1,xmm4	/* t0d = t07-rt */\
			\
			__asm	movaps	[ebx      ],xmm2	/* <- t06 */\
			__asm	movaps	[ebx+0x010],xmm3	/* <- t07 */\
			__asm	movaps	[ecx      ],xmm0	/* <- t0c */\
			__asm	movaps	[ecx+0x010],xmm1	/* <- t0d */\
			\
		/*\
			rt  =t08*c -t09*s;		it  =t08*s +t09*c;\
			tt  =t0e*c2-t0f*s2;		t0f =t0e*s2+t0f*c2;	t0e=tt;\
			t08 =rt+t0e;			t09 =it+t0f;\
			t0e =rt-t0e;			t0f =it-t0f;\
			t02 =t02+t08;			t03 =t03+t09;\
			t08 =t02+c3m1*t08;		t09 =t03+c3m1*t09;\
			rt  =s3*t0e;			it  =s3*t0f;\
			t0e =t08+it;			t0f =t09-rt;\
			t08 =t08-it;			t09 =t09+rt;\
		*/\
			__asm	add	eax, 0x20	/* t02 */\
			__asm	add	ebx, 0x20	/* t08 */\
			__asm	add	ecx, 0x20	/* t0e */\
			__asm	mov	edi, __c1\
			\
			__asm	movaps	xmm2,[ebx     ]	/* t08 */\
			__asm	movaps	xmm3,[ebx+0x10]	/* t09 */\
			__asm	movaps	xmm6,[edi     ]	/* c */\
			__asm	movaps	xmm7,[edi+0x10]	/* s */\
			__asm	movaps	xmm0,xmm2		/* copy of t08 */\
			__asm	movaps	xmm1,xmm3		/* copy of t09 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t08*c  */\
			__asm	mulpd	xmm3,xmm6		/* t09*c  */\
			__asm	mulpd	xmm0,xmm7		/* t08*s  */\
			__asm	mulpd	xmm1,xmm7		/* t09*s  */\
			__asm	subpd	xmm2,xmm1	/* xmm2 <- rt */\
			__asm	addpd	xmm3,xmm0	/* xmm3 <- it 	xmm6,7 free */\
			\
		__asm	add	edi, 0x20	/* c2 */\
			__asm	movaps	xmm4,[ecx     ]	/* t0e */\
			__asm	movaps	xmm5,[ecx+0x10]	/* t0f */\
			__asm	movaps	xmm6,[edi     ]	/* c2*/\
			__asm	movaps	xmm7,[edi+0x10]	/* s2*/\
			__asm	movaps	xmm0,xmm4		/* copy of t0e */\
			__asm	movaps	xmm1,xmm5		/* copy of t0f */\
			\
			__asm	mov	edx, __c3m1\
			__asm	mulpd	xmm4,xmm6		/* t0e*c2 */\
			__asm	mulpd	xmm5,xmm6		/* t0f*c2 */\
			__asm	mulpd	xmm0,xmm7		/* t0e*s2 */\
			__asm	mulpd	xmm1,xmm7		/* t0f*s2 */\
			__asm	subpd	xmm4,xmm1	/* ~t0e */\
			__asm	addpd	xmm5,xmm0	/* ~t0f 	xmm6,7 free */\
			\
			__asm	movaps	xmm6,[edx     ]	/* c3m1*/\
			__asm	movaps	xmm7,[edx+0x10]	/* s3  */\
			__asm	movaps	xmm0,[eax     ]	/* t02 */\
			__asm	movaps	xmm1,[eax+0x10]	/* t03 */\
			\
			__asm	subpd	xmm2,xmm4		/* t0e =rt-t0e */\
			__asm	subpd	xmm3,xmm5		/* t0f =it-t0f */\
			__asm	addpd	xmm4,xmm4		/*       2*t0e */\
			__asm	addpd	xmm5,xmm5		/*       2*t0f */\
			__asm	addpd	xmm4,xmm2		/* t08 =rt+t0e */\
			__asm	addpd	xmm5,xmm3		/* t09 =it+t0f */\
			__asm	addpd	xmm0,xmm4	/* t02 = t02+t08 */\
			__asm	addpd	xmm1,xmm5	/* t03 = t03+t09 */\
			__asm	movaps	[eax      ],xmm0	/* <- t02 */\
			__asm	movaps	[eax+0x010],xmm1	/* <- t03 */\
			\
			__asm	mulpd	xmm4,xmm6		/* t08 *= c3m1 */\
			__asm	mulpd	xmm5,xmm6		/* t09 *= c3m1 */\
			__asm	mulpd	xmm2,xmm7		/* rt = t0e*s3   */\
			__asm	mulpd	xmm3,xmm7		/* it = t0f*s3   */\
			__asm	addpd	xmm4,xmm0	/* t08 = t02+c3m1*t08 */\
			__asm	addpd	xmm5,xmm1	/* t09 = t03+c3m1*t09 */\
			\
			__asm	movaps	xmm0,xmm4	/* copy t08 */\
			__asm	movaps	xmm1,xmm5	/* copy t09 */\
			\
			__asm	subpd	xmm4,xmm3	/* t08 = t08-it */\
			__asm	addpd	xmm5,xmm2	/* t09 = t09+rt */\
			__asm	addpd	xmm0,xmm3	/* t0e = t08+it */\
			__asm	subpd	xmm1,xmm2	/* t0f = t09-rt */\
			\
			__asm	movaps	[ebx      ],xmm4	/* <- t08 */\
			__asm	movaps	[ebx+0x010],xmm5	/* <- t09 */\
			__asm	movaps	[ecx      ],xmm0	/* <- t0e */\
			__asm	movaps	[ecx+0x010],xmm1	/* <- t0f */\
			\
		/*\
			rt  =t0a*c2-t0b*s2;	it  =t0a*s2+t0b*c2;\
			tt  =t0g*c4-t0h*s4;	t0h =t0g*s4+t0h*c4;	t0g=tt;\
			t0a =rt+t0g;			t0b =it+t0h;\
			t0g =rt-t0g;			t0h =it-t0h;\
			t04 =t04+t0a;			t05 =t05+t0b;\
			t0a =t04+c3m1*t0a;	t0b =t05+c3m1*t0b;\
			rt  =s3*t0g;			it  =s3*t0h;\
			t0g =t0a+it;			t0h =t0b-rt;\
			t0a =t0a-it;			t0b =t0b+rt;\
		*/\
			__asm	add	eax, 0x20	/* t04 */\
			__asm	add	ebx, 0x20	/* t0a */\
			__asm	add	ecx, 0x20	/* t0g */\
			\
			__asm	movaps	xmm2,[ebx     ]	/* t0a */\
			__asm	movaps	xmm3,[ebx+0x10]	/* t0b */\
			__asm	movaps	xmm6,[edi     ]	/* c2*/\
			__asm	movaps	xmm7,[edi+0x10]	/* s2*/\
			__asm	movaps	xmm0,xmm2		/* copy of t0a */\
			__asm	movaps	xmm1,xmm3		/* copy of t0b */\
			\
			__asm	mulpd	xmm2,xmm6		/* t0a*c  */\
			__asm	mulpd	xmm3,xmm6		/* t0b*c  */\
			__asm	mulpd	xmm0,xmm7		/* t0a*s  */\
			__asm	mulpd	xmm1,xmm7		/* t0b*s  */\
			__asm	subpd	xmm2,xmm1	/* xmm2 <- rt */\
			__asm	addpd	xmm3,xmm0	/* xmm3 <- it 	xmm6,7 free */\
			\
		__asm	add	edi, 0x40	/* c4 */\
			__asm	movaps	xmm4,[ecx     ]	/* t0g */\
			__asm	movaps	xmm5,[ecx+0x10]	/* t0h */\
			__asm	movaps	xmm6,[edi     ]	/* c4*/\
			__asm	movaps	xmm7,[edi+0x10]	/* s4*/\
			__asm	movaps	xmm0,xmm4		/* copy of t0g */\
			__asm	movaps	xmm1,xmm5		/* copy of t0h */\
			\
			__asm	mulpd	xmm4,xmm6		/* t0g*c4 */\
			__asm	mulpd	xmm5,xmm6		/* t0h*c4 */\
			__asm	mulpd	xmm0,xmm7		/* t0g*s4 */\
			__asm	mulpd	xmm1,xmm7		/* t0h*s4 */\
			__asm	subpd	xmm4,xmm1	/* ~t0g */\
			__asm	addpd	xmm5,xmm0	/* ~t0h 	xmm6,7 free */\
			\
			__asm	movaps	xmm6,[edx     ]	/* c3m1*/\
			__asm	movaps	xmm7,[edx+0x10]	/* s3  */\
			__asm	movaps	xmm0,[eax     ]	/* t04 */\
			__asm	movaps	xmm1,[eax+0x10]	/* t05 */\
			\
			__asm	subpd	xmm2,xmm4		/* t0g =rt-t0g */\
			__asm	subpd	xmm3,xmm5		/* t0h =it-t0h */\
			__asm	addpd	xmm4,xmm4		/*       2*t0g */\
			__asm	addpd	xmm5,xmm5		/*       2*t0h */\
			__asm	addpd	xmm4,xmm2		/* t0a =rt+t0g */\
			__asm	addpd	xmm5,xmm3		/* t0b =it+t0h */\
			__asm	addpd	xmm0,xmm4	/* t04 = t04+t0a */\
			__asm	addpd	xmm1,xmm5	/* t05 = t05+t0b */\
			__asm	movaps	[eax      ],xmm0	/* <- t04 */\
			__asm	movaps	[eax+0x010],xmm1	/* <- t05 */\
			\
			__asm	mulpd	xmm4,xmm6		/* t0a *= c3m1 */\
			__asm	mulpd	xmm5,xmm6		/* t0b *= c3m1 */\
			__asm	mulpd	xmm2,xmm7		/* rt = t0g*s3   */\
			__asm	mulpd	xmm3,xmm7		/* it = t0h*s3   */\
			__asm	addpd	xmm4,xmm0	/* t0a = t04+c3m1*t0a */\
			__asm	addpd	xmm5,xmm1	/* t0b = t05+c3m1*t0b */\
			\
			__asm	movaps	xmm0,xmm4	/* copy t0a */\
			__asm	movaps	xmm1,xmm5	/* copy t0b */\
			\
			__asm	subpd	xmm4,xmm3	/* t0a = t0a-it */\
			__asm	addpd	xmm5,xmm2	/* t0b = t0b+rt */\
			__asm	addpd	xmm0,xmm3	/* t0g = t0a+it */\
			__asm	subpd	xmm1,xmm2	/* t0h = t0b-rt */\
			\
			__asm	movaps	[ebx      ],xmm4	/* <- t0a */\
			__asm	movaps	[ebx+0x010],xmm5	/* <- t0b */\
			__asm	movaps	[ecx      ],xmm0	/* <- t0g */\
			__asm	movaps	[ecx+0x010],xmm1	/* <- t0h */\
		}

		/*...Radix-9 DIF: Inputs enter in memory locations __i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8.\
				Outputs go into memory locations __o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8, which may coincide with inputs:\
		*/\
		#define SSE2_RADIX_09_DIF_B(__i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8, __c1,__c2,__c3m1,__c4, __o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8)\
		{\
		/*...gather the needed data (9 64-b__it complex, i.e. 18 64-b__it reals) and do three radix-3 transforms: */\
		/*\
			t00 =A0r;				t01 =A0i;\
			t02 =A3r+A6r;			t03 =A3i+A6i;\
			t04 =A3r-A6r;			t05 =A3i-A6i;\
			t00 =t00+t02;			t01 =t01+t03;\
			t02 =t00+c3m1*t02;		t03 =t01+c3m1*t03;\
			rt  =s3*t04;			it  =s3*t05;\
			t04 =t02+it;			t05 =t03-rt;\
			t02 =t02-it;			t03 =t03+rt;\
		*/\
			__asm	mov	eax, __i0\
			__asm	mov	ebx, __i3\
			__asm	mov	ecx, __i6\
			__asm	mov	edx, __c3m1\
			\
			__asm	movaps	xmm2,[ebx     ]	/* A3r */\
			__asm	movaps	xmm3,[ebx+0x10]	/* A3i */\
			__asm	movaps	xmm0,[eax     ]	/* A0r */\
			__asm	movaps	xmm1,[eax+0x10]	/* A0i */\
			__asm	movaps	xmm6,[ecx     ]	/* A6r */\
			__asm	movaps	xmm7,[ecx+0x10]	/* A6i */\
			__asm	movaps	xmm4,xmm2		/* t04: init to copy of A3r */\
			__asm	movaps	xmm5,xmm3		/* t05: init to copy of A3i */\
			\
			__asm	addpd	xmm2,xmm6		/* t02=A3r+A6r */\
			__asm	addpd	xmm3,xmm7		/* t03=A3i+A6i */\
			__asm	subpd	xmm4,xmm6		/* t04=A3r-A6r */\
			__asm	subpd	xmm5,xmm7		/* t05=A3i-A6i */\
			__asm	addpd	xmm0,xmm2	/* t00 = t00+t02 */\
			__asm	addpd	xmm1,xmm3	/* t01 = t01+t03 */\
			__asm	movaps	xmm6,[edx     ]	/* c3m1*/\
			__asm	movaps	xmm7,[edx+0x10]	/* s3  */\
			__asm	movaps	[eax      ],xmm0	/* <- t00 */\
			__asm	movaps	[eax+0x010],xmm1	/* <- t01 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t02 *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t03 *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = t04*s3   */\
			__asm	mulpd	xmm5,xmm7		/* it = t05*s3   */\
			__asm	addpd	xmm2,xmm0	/* t02 = t00+c3m1*t02 */\
			__asm	addpd	xmm3,xmm1	/* t03 = t01+c3m1*t03 */\
			\
			__asm	movaps	xmm0,xmm2	/* copy t02 */\
			__asm	movaps	xmm1,xmm3	/* copy t03 */\
			\
			__asm	subpd	xmm2,xmm5	/* t02 = t02-it */\
			__asm	addpd	xmm3,xmm4	/* t03 = t03+rt */\
			__asm	addpd	xmm0,xmm5	/* t04 = t02+it */\
			__asm	subpd	xmm1,xmm4	/* t05 = t03-rt */\
			\
			__asm	movaps	[ebx      ],xmm2	/* <- t02 */\
			__asm	movaps	[ebx+0x010],xmm3	/* <- t03 */\
			__asm	movaps	[ecx      ],xmm0	/* <- t04 */\
			__asm	movaps	[ecx+0x010],xmm1	/* <- t05 */\
			\
		/*\
			t06 =A1r;				t07 =A1i;\
			t08 =A4r+A7r;			t09 =A4i+A7i;\
			t0a =A4r-A7r;			t0b =A4i-A7i;\
			t06 =t06+t08;			t07 =t07+t09;\
			t08 =t06+c3m1*t08;		t09 =t07+c3m1*t09;\
			rt  =s3*t0a;			it  =s3*t0b;\
			t0a =t08+it;			t0b =t09-rt;\
			t08 =t08-it;			t09 =t09+rt;\
		*/\
			__asm	mov	eax, __i1\
			__asm	mov	ebx, __i4\
			__asm	mov	ecx, __i7\
			\
			__asm	movaps	xmm4,[ebx     ]	/* A4r */\
			__asm	movaps	xmm5,[ebx+0x10]	/* A4i */\
			__asm	movaps	xmm0,[eax     ]	/* A1r */\
			__asm	movaps	xmm1,[eax+0x10]	/* A1i */\
			__asm	movaps	xmm2,[ecx     ]	/* A7r */\
			__asm	movaps	xmm3,[ecx+0x10]	/* A7i */\
			\
			__asm	subpd	xmm4,xmm2		/* t0a=A4r-A7r */\
			__asm	subpd	xmm5,xmm3		/* t0b=A4i-A7i */\
			__asm	addpd	xmm2,xmm2		/*       2*A7r */\
			__asm	addpd	xmm3,xmm3		/*       2*A7i */\
			__asm	addpd	xmm2,xmm4		/* t08=A4r+A7r */\
			__asm	addpd	xmm3,xmm5		/* t09=A4i+A7i */\
			__asm	addpd	xmm0,xmm2	/* t06 = t06+t08 */\
			__asm	addpd	xmm1,xmm3	/* t07 = t07+t09 */\
			__asm	movaps	[eax      ],xmm0	/* <- t06 */\
			__asm	movaps	[eax+0x010],xmm1	/* <- t07 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t08 *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t09 *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = t0a*s3   */\
			__asm	mulpd	xmm5,xmm7		/* it = t0b*s3   */\
			__asm	addpd	xmm2,xmm0	/* t08 = t06+c3m1*t08 */\
			__asm	addpd	xmm3,xmm1	/* t09 = t07+c3m1*t09 */\
			\
			__asm	movaps	xmm0,xmm2	/* copy t08 */\
			__asm	movaps	xmm1,xmm3	/* copy t09 */\
			\
			__asm	subpd	xmm2,xmm5	/* t08 = t08-it */\
			__asm	addpd	xmm3,xmm4	/* t09 = t09+rt */\
			__asm	addpd	xmm0,xmm5	/* t0a = t08+it */\
			__asm	subpd	xmm1,xmm4	/* t0b = t09-rt */\
			\
			__asm	movaps	[ebx      ],xmm2	/* <- t08 */\
			__asm	movaps	[ebx+0x010],xmm3	/* <- t09 */\
			__asm	movaps	[ecx      ],xmm0	/* <- t0a */\
			__asm	movaps	[ecx+0x010],xmm1	/* <- t0b */\
			\
		/*\
			t0c =A2r;				t0d =A2i;\
			t0e =A5r+A8r;			t0f =A5i+A8i;\
			t0g =A5r-A8r;			t0h =A5i-A8i;\
			t0c =t0c+t0e;			t0d =t0d+t0f;\
			t0e =t0c+c3m1*t0e;		t0f =t0d+c3m1*t0f;\
			rt  =s3*t0g;			it  =s3*t0h;\
			t0g =t0e+it;			t0h =t0f-rt;\
			t0e =t0e-it;			t0f =t0f+rt;\
		*/\
			__asm	mov	eax, __i2\
			__asm	mov	ebx, __i5\
			__asm	mov	ecx, __i8\
			\
			__asm	movaps	xmm4,[ebx     ]	/* A5r */\
			__asm	movaps	xmm5,[ebx+0x10]	/* A5i */\
			__asm	movaps	xmm0,[eax     ]	/* A2r */\
			__asm	movaps	xmm1,[eax+0x10]	/* A2i */\
			__asm	movaps	xmm2,[ecx     ]	/* A8r */\
			__asm	movaps	xmm3,[ecx+0x10]	/* A8i */\
			\
			__asm	subpd	xmm4,xmm2		/* t0g=A5r-A8r */\
			__asm	subpd	xmm5,xmm3		/* t0h=A5i-A8i */\
			__asm	addpd	xmm2,xmm2		/*       2*A8r */\
			__asm	addpd	xmm3,xmm3		/*       2*A8i */\
			__asm	addpd	xmm2,xmm4		/* t0e=A5r+A8r */\
			__asm	addpd	xmm3,xmm5		/* t0f=A5i+A8i */\
			__asm	addpd	xmm0,xmm2	/* t0c = t0c+t0e */\
			__asm	addpd	xmm1,xmm3	/* t0d = t0d+t0f */\
			__asm	movaps	[eax      ],xmm0	/* <- t0c */\
			__asm	movaps	[eax+0x010],xmm1	/* <- t0d */\
			\
			__asm	mulpd	xmm2,xmm6		/* t0e *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t0f *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = t0g*s3   */\
			__asm	mulpd	xmm5,xmm7		/* it = t0h*s3   */\
			__asm	addpd	xmm2,xmm0	/* t0e = t0c+c3m1*t0e */\
			__asm	addpd	xmm3,xmm1	/* t0f = t0d+c3m1*t0f */\
			\
			__asm	movaps	xmm0,xmm2	/* copy t0e */\
			__asm	movaps	xmm1,xmm3	/* copy t0f */\
			\
			__asm	subpd	xmm2,xmm5	/* t0e = t0e-it */\
			__asm	addpd	xmm3,xmm4	/* t0f = t0f+rt */\
			__asm	addpd	xmm0,xmm5	/* t0g = t0e+it */\
			__asm	subpd	xmm1,xmm4	/* t0h = t0f-rt */\
			\
			__asm	movaps	[ebx      ],xmm2	/* <- t0e */\
			__asm	movaps	[ebx+0x010],xmm3	/* <- t0f */\
			__asm	movaps	[ecx      ],xmm0	/* <- t0g */\
			__asm	movaps	[ecx+0x010],xmm1	/* <- t0h */\
			\
		/******************************************************************************/\
		/*...and now do three more radix-3 transforms, including the twiddle factors: */\
		/******************************************************************************/\
		/*\
			rt  =t06;				it  =t07;\
			t06 =rt+t0c;			t07 =it+t0d;\
			t0c =rt-t0c;			t0d =it-t0d;\
			t00 =t00+t06;			t01 =t01+t07;\
			t06 =t00+c3m1*t06;		t07 =t01+c3m1*t07;\
			rt  =s3*t0c;			it  =s3*t0d;\
			t0c =t06+it;			t0d =t07-rt;\
			t06 =t06-it;			t07 =t07+rt;\
		*/\
		__asm	mov	eax, __i0\
		__asm	mov	ebx, __i1\
		__asm	mov	ecx, __i2\
			\
			__asm	movaps	xmm4,[ebx     ]	/* t06 */\
			__asm	movaps	xmm5,[ebx+0x10]	/* t07 */\
			__asm	movaps	xmm2,[ecx     ]	/* t0c */\
			__asm	movaps	xmm3,[ecx+0x10]	/* t0d */\
			__asm	movaps	xmm0,[eax     ]	/* t00 */\
			__asm	movaps	xmm1,[eax+0x10]	/* t01 */\
			\
		__asm	mov	esi, __o0\
		__asm	mov	edi, __o1\
		__asm	mov	edx, __o2\
			\
			__asm	subpd	xmm4,xmm2		/* t0c=t06-t0c */\
			__asm	subpd	xmm5,xmm3		/* t0d=t07-t0d */\
			__asm	addpd	xmm2,xmm2		/*       2*t0c */\
			__asm	addpd	xmm3,xmm3		/*       2*t0d */\
			__asm	addpd	xmm2,xmm4		/* t06=t06+t0c */\
			__asm	addpd	xmm3,xmm5		/* t07=t07+t0d */\
			__asm	addpd	xmm0,xmm2	/* t00 = t00+t06 */\
			__asm	addpd	xmm1,xmm3	/* t01 = t01+t07 */\
			__asm	movaps	[esi      ],xmm0	/* <- t00 */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t01 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t06 *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t07 *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = t0c*s3   */\
			__asm	mulpd	xmm5,xmm7		/* it = t0d*s3   */\
			__asm	addpd	xmm2,xmm0	/* t06 = t00+c3m1*t06 */\
			__asm	addpd	xmm3,xmm1	/* t07 = t01+c3m1*t07 */\
			\
			__asm	movaps	xmm0,xmm2	/* copy t06 */\
			__asm	movaps	xmm1,xmm3	/* copy t07 */\
			\
			__asm	subpd	xmm2,xmm5	/* t06 = t06-it */\
			__asm	addpd	xmm3,xmm4	/* t07 = t07+rt */\
			__asm	addpd	xmm0,xmm5	/* t0c = t06+it */\
			__asm	subpd	xmm1,xmm4	/* t0d = t07-rt */\
			\
			__asm	movaps	[edi      ],xmm2	/* <- t06 */\
			__asm	movaps	[edi+0x010],xmm3	/* <- t07 */\
			__asm	movaps	[edx      ],xmm0	/* <- t0c */\
			__asm	movaps	[edx+0x010],xmm1	/* <- t0d */\
			\
		/*\
			rt  =t08*c -t09*s;		it  =t08*s +t09*c;\
			tt  =t0e*c2-t0f*s2;		t0f =t0e*s2+t0f*c2;	t0e=tt;\
			t08 =rt+t0e;			t09 =it+t0f;\
			t0e =rt-t0e;			t0f =it-t0f;\
			t02 =t02+t08;			t03 =t03+t09;\
			t08 =t02+c3m1*t08;		t09 =t03+c3m1*t09;\
			rt  =s3*t0e;			it  =s3*t0f;\
			t0e =t08+it;			t0f =t09-rt;\
			t08 =t08-it;			t09 =t09+rt;\
		*/\
		__asm	mov	eax, __i3\
		__asm	mov	ebx, __i4\
		__asm	mov	ecx, __i5\
			__asm	mov	edi, __c1\
			\
			__asm	movaps	xmm2,[ebx     ]	/* t08 */\
			__asm	movaps	xmm3,[ebx+0x10]	/* t09 */\
			__asm	movaps	xmm6,[edi     ]	/* c */\
			__asm	movaps	xmm7,[edi+0x10]	/* s */\
			__asm	movaps	xmm0,xmm2		/* copy of t08 */\
			__asm	movaps	xmm1,xmm3		/* copy of t09 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t08*c  */\
			__asm	mulpd	xmm3,xmm6		/* t09*c  */\
			__asm	mulpd	xmm0,xmm7		/* t08*s  */\
			__asm	mulpd	xmm1,xmm7		/* t09*s  */\
			__asm	subpd	xmm2,xmm1	/* xmm2 <- rt */\
			__asm	addpd	xmm3,xmm0	/* xmm3 <- it 	xmm6,7 free */\
			\
		__asm	add	edi, 0x20	/* c2 */\
			__asm	movaps	xmm4,[ecx     ]	/* t0e */\
			__asm	movaps	xmm5,[ecx+0x10]	/* t0f */\
			__asm	movaps	xmm6,[edi     ]	/* c2*/\
			__asm	movaps	xmm7,[edi+0x10]	/* s2*/\
			__asm	movaps	xmm0,xmm4		/* copy of t0e */\
			__asm	movaps	xmm1,xmm5		/* copy of t0f */\
			\
			__asm	mov	edx, __c3m1\
			__asm	mulpd	xmm4,xmm6		/* t0e*c2 */\
			__asm	mulpd	xmm5,xmm6		/* t0f*c2 */\
			__asm	mulpd	xmm0,xmm7		/* t0e*s2 */\
			__asm	mulpd	xmm1,xmm7		/* t0f*s2 */\
			__asm	subpd	xmm4,xmm1	/* ~t0e */\
			__asm	addpd	xmm5,xmm0	/* ~t0f 	xmm6,7 free */\
			\
			__asm	movaps	xmm6,[edx     ]	/* c3m1*/\
			__asm	movaps	xmm7,[edx+0x10]	/* s3  */\
			__asm	movaps	xmm0,[eax     ]	/* t02 */\
			__asm	movaps	xmm1,[eax+0x10]	/* t03 */\
			\
		__asm	mov	esi, __o3\
		__asm	mov	edi, __o4\
		__asm	mov	edx, __o5\
			\
			__asm	subpd	xmm2,xmm4		/* t0e =rt-t0e */\
			__asm	subpd	xmm3,xmm5		/* t0f =it-t0f */\
			__asm	addpd	xmm4,xmm4		/*       2*t0e */\
			__asm	addpd	xmm5,xmm5		/*       2*t0f */\
			__asm	addpd	xmm4,xmm2		/* t08 =rt+t0e */\
			__asm	addpd	xmm5,xmm3		/* t09 =it+t0f */\
			__asm	addpd	xmm0,xmm4	/* t02 = t02+t08 */\
			__asm	addpd	xmm1,xmm5	/* t03 = t03+t09 */\
			__asm	movaps	[esi      ],xmm0	/* <- t02 */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t03 */\
			\
			__asm	mulpd	xmm4,xmm6		/* t08 *= c3m1 */\
			__asm	mulpd	xmm5,xmm6		/* t09 *= c3m1 */\
			__asm	mulpd	xmm2,xmm7		/* rt = t0e*s3   */\
			__asm	mulpd	xmm3,xmm7		/* it = t0f*s3   */\
			__asm	addpd	xmm4,xmm0	/* t08 = t02+c3m1*t08 */\
			__asm	addpd	xmm5,xmm1	/* t09 = t03+c3m1*t09 */\
			\
			__asm	movaps	xmm0,xmm4	/* copy t08 */\
			__asm	movaps	xmm1,xmm5	/* copy t09 */\
			\
			__asm	subpd	xmm4,xmm3	/* t08 = t08-it */\
			__asm	addpd	xmm5,xmm2	/* t09 = t09+rt */\
			__asm	addpd	xmm0,xmm3	/* t0e = t08+it */\
			__asm	subpd	xmm1,xmm2	/* t0f = t09-rt */\
			\
			__asm	movaps	[edi      ],xmm4	/* <- t08 */\
			__asm	movaps	[edi+0x010],xmm5	/* <- t09 */\
			__asm	movaps	[edx      ],xmm0	/* <- t0e */\
			__asm	movaps	[edx+0x010],xmm1	/* <- t0f */\
			\
		/*\
			rt  =t0a*c2-t0b*s2;	it  =t0a*s2+t0b*c2;\
			tt  =t0g*c4-t0h*s4;	t0h =t0g*s4+t0h*c4;	t0g=tt;\
			t0a =rt+t0g;			t0b =it+t0h;\
			t0g =rt-t0g;			t0h =it-t0h;\
			t04 =t04+t0a;			t05 =t05+t0b;\
			t0a =t04+c3m1*t0a;	t0b =t05+c3m1*t0b;\
			rt  =s3*t0g;			it  =s3*t0h;\
			t0g =t0a+it;			t0h =t0b-rt;\
			t0a =t0a-it;			t0b =t0b+rt;\
		*/\
		__asm	mov	eax, __i6\
		__asm	mov	ebx, __i7\
		__asm	mov	ecx, __i8\
			__asm	mov	edi, __c2\
			__asm	mov	edx, __c3m1\
			\
			__asm	movaps	xmm2,[ebx     ]	/* t0a */\
			__asm	movaps	xmm3,[ebx+0x10]	/* t0b */\
			__asm	movaps	xmm6,[edi     ]	/* c2*/\
			__asm	movaps	xmm7,[edi+0x10]	/* s2*/\
			__asm	movaps	xmm0,xmm2		/* copy of t0a */\
			__asm	movaps	xmm1,xmm3		/* copy of t0b */\
			\
			__asm	mulpd	xmm2,xmm6		/* t0a*c  */\
			__asm	mulpd	xmm3,xmm6		/* t0b*c  */\
			__asm	mulpd	xmm0,xmm7		/* t0a*s  */\
			__asm	mulpd	xmm1,xmm7		/* t0b*s  */\
			__asm	subpd	xmm2,xmm1	/* xmm2 <- rt */\
			__asm	addpd	xmm3,xmm0	/* xmm3 <- it 	xmm6,7 free */\
			\
		__asm	add	edi, 0x40	/* c4 */\
			__asm	movaps	xmm4,[ecx     ]	/* t0g */\
			__asm	movaps	xmm5,[ecx+0x10]	/* t0h */\
			__asm	movaps	xmm6,[edi     ]	/* c4*/\
			__asm	movaps	xmm7,[edi+0x10]	/* s4*/\
			__asm	movaps	xmm0,xmm4		/* copy of t0g */\
			__asm	movaps	xmm1,xmm5		/* copy of t0h */\
			\
			__asm	mulpd	xmm4,xmm6		/* t0g*c4 */\
			__asm	mulpd	xmm5,xmm6		/* t0h*c4 */\
			__asm	mulpd	xmm0,xmm7		/* t0g*s4 */\
			__asm	mulpd	xmm1,xmm7		/* t0h*s4 */\
			__asm	subpd	xmm4,xmm1	/* ~t0g */\
			__asm	addpd	xmm5,xmm0	/* ~t0h 	xmm6,7 free */\
			\
			__asm	movaps	xmm6,[edx     ]	/* c3m1*/\
			__asm	movaps	xmm7,[edx+0x10]	/* s3  */\
			__asm	movaps	xmm0,[eax     ]	/* t04 */\
			__asm	movaps	xmm1,[eax+0x10]	/* t05 */\
			\
		__asm	mov	esi, __o6\
		__asm	mov	edi, __o7\
		__asm	mov	edx, __o8\
			\
			__asm	subpd	xmm2,xmm4		/* t0g =rt-t0g */\
			__asm	subpd	xmm3,xmm5		/* t0h =it-t0h */\
			__asm	addpd	xmm4,xmm4		/*       2*t0g */\
			__asm	addpd	xmm5,xmm5		/*       2*t0h */\
			__asm	addpd	xmm4,xmm2		/* t0a =rt+t0g */\
			__asm	addpd	xmm5,xmm3		/* t0b =it+t0h */\
			__asm	addpd	xmm0,xmm4	/* t04 = t04+t0a */\
			__asm	addpd	xmm1,xmm5	/* t05 = t05+t0b */\
			__asm	movaps	[esi      ],xmm0	/* <- t04 */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t05 */\
			\
			__asm	mulpd	xmm4,xmm6		/* t0a *= c3m1 */\
			__asm	mulpd	xmm5,xmm6		/* t0b *= c3m1 */\
			__asm	mulpd	xmm2,xmm7		/* rt = t0g*s3   */\
			__asm	mulpd	xmm3,xmm7		/* it = t0h*s3   */\
			__asm	addpd	xmm4,xmm0	/* t0a = t04+c3m1*t0a */\
			__asm	addpd	xmm5,xmm1	/* t0b = t05+c3m1*t0b */\
			\
			__asm	movaps	xmm0,xmm4	/* copy t0a */\
			__asm	movaps	xmm1,xmm5	/* copy t0b */\
			\
			__asm	subpd	xmm4,xmm3	/* t0a = t0a-it */\
			__asm	addpd	xmm5,xmm2	/* t0b = t0b+rt */\
			__asm	addpd	xmm0,xmm3	/* t0g = t0a+it */\
			__asm	subpd	xmm1,xmm2	/* t0h = t0b-rt */\
			\
			__asm	movaps	[edi      ],xmm4	/* <- t0a */\
			__asm	movaps	[edi+0x010],xmm5	/* <- t0b */\
			__asm	movaps	[edx      ],xmm0	/* <- t0g */\
			__asm	movaps	[edx+0x010],xmm1	/* <- t0h */\
		}

		/*...Radix-9 DIT: t's are temporaries assumed to be adjacent in memory, i.e. reals separated by 0x20 bytes.\
				outputs go into memory locations __o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8, assumed disjoint with inputs:\
		*/\
		#define SSE2_RADIX_09_DIT(__t00, __c1,__c2,__c3m1,__c4, __o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8)\
		{\
		/*...gather the needed data (9 64-b__it complex, i.e. 18 64-b__it reals) and do three radix-3 transforms: */\
		/*\
			rt  =t02;				it  =t03;\
			t02 =rt+t04;			t03 =it+t05;\
			t04 =rt-t04;			t05 =it-t05;\
			t00 =t00+t02;			t01 =t01+t03;\
			t02 =t00+c3m1*t02;	t03 =t01+c3m1*t03;\
			rt  =s3*t04;			it  =s3*t05;\
			t04 =t02-it;			t05 =t03+rt;\
			t02 =t02+it;			t03 =t03-rt;\
		*/\
			__asm	mov	eax, __t00\
			__asm	mov	edx, __c3m1\
			__asm	mov	ebx, eax\
			__asm	mov	ecx, eax\
			__asm	add	ebx, 0x20\
			__asm	add	ecx, 0x40\
			\
			__asm	movaps	xmm2,[ebx     ]	/* A3r */\
			__asm	movaps	xmm3,[ebx+0x10]	/* A3i */\
			__asm	movaps	xmm0,[eax     ]	/* A0r */\
			__asm	movaps	xmm1,[eax+0x10]	/* A0i */\
			__asm	movaps	xmm6,[ecx     ]	/* A6r */\
			__asm	movaps	xmm7,[ecx+0x10]	/* A6i */\
			__asm	movaps	xmm4,xmm2		/* t04: init to copy of A3r */\
			__asm	movaps	xmm5,xmm3		/* t05: init to copy of A3i */\
			\
			__asm	addpd	xmm2,xmm6		/* t02=A3r+A6r */\
			__asm	addpd	xmm3,xmm7		/* t03=A3i+A6i */\
			__asm	subpd	xmm4,xmm6		/* t04=A3r-A6r */\
			__asm	subpd	xmm5,xmm7		/* t05=A3i-A6i */\
			__asm	addpd	xmm0,xmm2	/* t00 = t00+t02 */\
			__asm	addpd	xmm1,xmm3	/* t01 = t01+t03 */\
			__asm	movaps	xmm6,[edx     ]	/* c3m1*/\
			__asm	movaps	xmm7,[edx+0x10]	/* s3  */\
			__asm	movaps	[eax      ],xmm0	/* <- t00 */\
			__asm	movaps	[eax+0x010],xmm1	/* <- t01 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t02 *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t03 *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = t04*s3   */\
			__asm	mulpd	xmm5,xmm7		/* it = t05*s3   */\
			__asm	addpd	xmm2,xmm0	/* t02 = t00+c3m1*t02 */\
			__asm	addpd	xmm3,xmm1	/* t03 = t01+c3m1*t03 */\
			\
			__asm	movaps	xmm0,xmm2	/* copy t02 */\
			__asm	movaps	xmm1,xmm3	/* copy t03 */\
			\
			__asm	addpd	xmm2,xmm5	/* t02 = t02-it */\
			__asm	subpd	xmm3,xmm4	/* t03 = t03+rt */\
			__asm	subpd	xmm0,xmm5	/* t04 = t02+it */\
			__asm	addpd	xmm1,xmm4	/* t05 = t03-rt */\
			\
			__asm	movaps	[ebx      ],xmm2	/* <- t02 */\
			__asm	movaps	[ebx+0x010],xmm3	/* <- t03 */\
			__asm	movaps	[ecx      ],xmm0	/* <- t04 */\
			__asm	movaps	[ecx+0x010],xmm1	/* <- t05 */\
			\
		/*\
			rt  =t08;				it  =t09;\
			t08 =rt+t0a;			t09 =it+t0b;\
			t0a =rt-t0a;			t0b =it-t0b;\
			t06 =t06+t08;			t07 =t07+t09;\
			t08 =t06+c3m1*t08;	t09 =t07+c3m1*t09;\
			rt  =s3*t0a;			it  =s3*t0b;\
			t0a =t08-it;			t0b =t09+rt;\
			t08 =t08+it;			t09 =t09-rt;\
		*/\
			__asm	add	eax, 0x60	/* t06 */\
			__asm	add	ebx, 0x60	/* t08 */\
			__asm	add	ecx, 0x60	/* t0a */\
			\
			__asm	movaps	xmm4,[ebx     ]	/* A3r */\
			__asm	movaps	xmm5,[ebx+0x10]	/* A3i */\
			__asm	movaps	xmm2,[ecx     ]	/* A6r */\
			__asm	movaps	xmm3,[ecx+0x10]	/* A6i */\
			__asm	movaps	xmm0,[eax     ]	/* A0r */\
			__asm	movaps	xmm1,[eax+0x10]	/* A0i */\
			\
			__asm	subpd	xmm4,xmm2		/* t04=A3r-A6r */\
			__asm	subpd	xmm5,xmm3		/* t05=A3i-A6i */\
			__asm	addpd	xmm2,xmm2		/*       2*A6r */\
			__asm	addpd	xmm3,xmm3		/*       2*A6i */\
			__asm	addpd	xmm2,xmm4		/* t02=A3r+A6r */\
			__asm	addpd	xmm3,xmm5		/* t03=A3i+A6i */\
			__asm	addpd	xmm0,xmm2	/* t00 = t00+t02 */\
			__asm	addpd	xmm1,xmm3	/* t01 = t01+t03 */\
			__asm	movaps	[eax      ],xmm0	/* <- t00 */\
			__asm	movaps	[eax+0x010],xmm1	/* <- t01 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t02 *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t03 *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = t04*s3   */\
			__asm	mulpd	xmm5,xmm7		/* it = t05*s3   */\
			__asm	addpd	xmm2,xmm0	/* t02 = t00+c3m1*t02 */\
			__asm	addpd	xmm3,xmm1	/* t03 = t01+c3m1*t03 */\
			\
			__asm	movaps	xmm0,xmm2	/* copy t02 */\
			__asm	movaps	xmm1,xmm3	/* copy t03 */\
			\
			__asm	addpd	xmm2,xmm5	/* t02 = t02-it */\
			__asm	subpd	xmm3,xmm4	/* t03 = t03+rt */\
			__asm	subpd	xmm0,xmm5	/* t04 = t02+it */\
			__asm	addpd	xmm1,xmm4	/* t05 = t03-rt */\
			\
			__asm	movaps	[ebx      ],xmm2	/* <- t02 */\
			__asm	movaps	[ebx+0x010],xmm3	/* <- t03 */\
			__asm	movaps	[ecx      ],xmm0	/* <- t04 */\
			__asm	movaps	[ecx+0x010],xmm1	/* <- t05 */\
			\
		/*\
			rt  =t0e;				it  =t0f;\
			t0e =rt+t0g;			t0f =it+t0h;\
			t0g =rt-t0g;			t0h =it-t0h;\
			t0c =t0c+t0e;			t0d =t0d+t0f;\
			t0e =t0c+c3m1*t0e;	t0f =t0d+c3m1*t0f;\
			rt  =s3*t0g;			it  =s3*t0h;\
			t0g =t0e-it;			t0h =t0f+rt;\
			t0e =t0e+it;			t0f =t0f-rt;\
		*/\
			__asm	add	eax, 0x60	/* t06 */\
			__asm	add	ebx, 0x60	/* t08 */\
			__asm	add	ecx, 0x60	/* t0a */\
			\
			__asm	movaps	xmm4,[ebx     ]	/* A3r */\
			__asm	movaps	xmm5,[ebx+0x10]	/* A3i */\
			__asm	movaps	xmm2,[ecx     ]	/* A6r */\
			__asm	movaps	xmm3,[ecx+0x10]	/* A6i */\
			__asm	movaps	xmm0,[eax     ]	/* A0r */\
			__asm	movaps	xmm1,[eax+0x10]	/* A0i */\
			\
			__asm	subpd	xmm4,xmm2		/* t04=A3r-A6r */\
			__asm	subpd	xmm5,xmm3		/* t05=A3i-A6i */\
			__asm	addpd	xmm2,xmm2		/*       2*A6r */\
			__asm	addpd	xmm3,xmm3		/*       2*A6i */\
			__asm	addpd	xmm2,xmm4		/* t02=A3r+A6r */\
			__asm	addpd	xmm3,xmm5		/* t03=A3i+A6i */\
			__asm	addpd	xmm0,xmm2	/* t00 = t00+t02 */\
			__asm	addpd	xmm1,xmm3	/* t01 = t01+t03 */\
			__asm	movaps	[eax      ],xmm0	/* <- t00 */\
			__asm	movaps	[eax+0x010],xmm1	/* <- t01 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t02 *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t03 *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = t04*s3   */\
			__asm	mulpd	xmm5,xmm7		/* it = t05*s3   */\
			__asm	addpd	xmm2,xmm0	/* t02 = t00+c3m1*t02 */\
			__asm	addpd	xmm3,xmm1	/* t03 = t01+c3m1*t03 */\
			\
			__asm	movaps	xmm0,xmm2	/* copy t02 */\
			__asm	movaps	xmm1,xmm3	/* copy t03 */\
			\
			__asm	addpd	xmm2,xmm5	/* t02 = t02-it */\
			__asm	subpd	xmm3,xmm4	/* t03 = t03+rt */\
			__asm	subpd	xmm0,xmm5	/* t04 = t02+it */\
			__asm	addpd	xmm1,xmm4	/* t05 = t03-rt */\
			\
			__asm	movaps	[ebx      ],xmm2	/* <- t02 */\
			__asm	movaps	[ebx+0x010],xmm3	/* <- t03 */\
			__asm	movaps	[ecx      ],xmm0	/* <- t04 */\
			__asm	movaps	[ecx+0x010],xmm1	/* <- t05 */\
			\
		/******************************************************************************/\
		/*...and now do three more radix-3 transforms, including the twiddle factors: */\
		/******************************************************************************/\
		/*\
			rt  =t06;				it  =t07;\
			t06 =rt+t0c;			t07 =it+t0d;\
			t0c =rt-t0c;			t0d =it-t0d;\
			t00 =t00+t06;			t01 =t01+t07;\
			A0r =t00;				A0i =t01;\
			t06 =t00+c3m1*t06;		t07 =t01+c3m1*t07;\
			rt  =s3*t0c;			it  =s3*t0d;\
			A3r =t06+it;			A3i =t07-rt;\
			A6r =t06-it;			A6i =t07+rt;\
		*/\
			__asm	mov	eax, __t00	/* t00 */\
			__asm	mov	ebx, eax\
			__asm	mov	ecx, eax\
			__asm	add	ebx, 0x60	/* t06 */\
			__asm	add	ecx, 0xc0	/* t0c */\
			\
		__asm	mov	esi, __o0\
			\
			__asm	movaps	xmm4,[ebx     ]	/* A3r */\
			__asm	movaps	xmm5,[ebx+0x10]	/* A3i */\
			__asm	movaps	xmm2,[ecx     ]	/* A6r */\
			__asm	movaps	xmm3,[ecx+0x10]	/* A6i */\
			__asm	movaps	xmm0,[eax     ]	/* A0r */\
			__asm	movaps	xmm1,[eax+0x10]	/* A0i */\
			\
			__asm	subpd	xmm4,xmm2		/* t04=A3r-A6r */\
			__asm	subpd	xmm5,xmm3		/* t05=A3i-A6i */\
			__asm	addpd	xmm2,xmm2		/*       2*A6r */\
			__asm	addpd	xmm3,xmm3		/*       2*A6i */\
			__asm	addpd	xmm2,xmm4		/* t02=A3r+A6r */\
			__asm	addpd	xmm3,xmm5		/* t03=A3i+A6i */\
			__asm	addpd	xmm0,xmm2	/* t00 = t00+t02 */\
			__asm	addpd	xmm1,xmm3	/* t01 = t01+t03 */\
			__asm	movaps	[esi      ],xmm0	/* <- t00 */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t01 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t02 *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t03 *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = t04*s3   */\
			__asm	mulpd	xmm5,xmm7		/* it = t05*s3   */\
			__asm	addpd	xmm2,xmm0	/* t02 = t00+c3m1*t02 */\
			__asm	addpd	xmm3,xmm1	/* t03 = t01+c3m1*t03 */\
			\
			__asm	movaps	xmm0,xmm2	/* copy t02 */\
			__asm	movaps	xmm1,xmm3	/* copy t03 */\
			\
		__asm	mov	edi, __o3\
		__asm	mov	esi, __o6\
			\
			__asm	addpd	xmm2,xmm5	/* t06 = t06+it */\
			__asm	subpd	xmm3,xmm4	/* t07 = t07-rt */\
			__asm	subpd	xmm0,xmm5	/* t0c = t06-it */\
			__asm	addpd	xmm1,xmm4	/* t0d = t07+rt */\
			\
			__asm	movaps	[edi      ],xmm2	/* <- t06 */\
			__asm	movaps	[edi+0x010],xmm3	/* <- t07 */\
			__asm	movaps	[esi      ],xmm0	/* <- t0c */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t0d */\
			\
		/*\
			rt  =t08*c +t09*s;		it  =t09*c -t08*s;\
			tt  =t0e*c2+t0f*s2;		t0f =t0f*c2-t0e*s2;	t0e=tt;\
			t08 =rt+t0e;			t09 =it+t0f;\
			t0e =rt-t0e;			t0f =it-t0f;\
			t02 =t02+t08;			t03 =t03+t09;\
			A7r =t02;				A7i =t03;\
			t08 =t02+c3m1*t08;		t09 =t03+c3m1*t09;\
			rt  =s3*t0e;			it  =s3*t0f;\
			A1r =t08+it;			A1i =t09-rt;\
			A4r =t08-it;			A4i =t09+rt;\
		*/\
			__asm	add	eax, 0x20	/* t02 */\
			__asm	add	ebx, 0x20	/* t08 */\
			__asm	add	ecx, 0x20	/* t0e */\
			__asm	mov	edi, __c1\
			\
			__asm	movaps	xmm2,[ebx     ]	/* t08 */\
			__asm	movaps	xmm3,[ebx+0x10]	/* t09 */\
			__asm	movaps	xmm6,[edi     ]	/* c */\
			__asm	movaps	xmm7,[edi+0x10]	/* s */\
			__asm	movaps	xmm0,xmm2		/* copy of t08 */\
			__asm	movaps	xmm1,xmm3		/* copy of t09 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t08*c  */\
			__asm	mulpd	xmm3,xmm6		/* t09*c  */\
			__asm	mulpd	xmm0,xmm7		/* t08*s  */\
			__asm	mulpd	xmm1,xmm7		/* t09*s  */\
			__asm	addpd	xmm2,xmm1	/* xmm2 <- rt */\
			__asm	subpd	xmm3,xmm0	/* xmm3 <- it 	xmm6,7 free */\
			\
		__asm	add	edi, 0x20	/* c2 */\
			__asm	movaps	xmm4,[ecx     ]	/* t0e */\
			__asm	movaps	xmm5,[ecx+0x10]	/* t0f */\
			__asm	movaps	xmm6,[edi     ]	/* c2*/\
			__asm	movaps	xmm7,[edi+0x10]	/* s2*/\
			__asm	movaps	xmm0,xmm4		/* copy of t0e */\
			__asm	movaps	xmm1,xmm5		/* copy of t0f */\
			\
			__asm	mulpd	xmm4,xmm6		/* t0e*c2 */\
			__asm	mulpd	xmm5,xmm6		/* t0f*c2 */\
			__asm	mulpd	xmm0,xmm7		/* t0e*s2 */\
			__asm	mulpd	xmm1,xmm7		/* t0f*s2 */\
			__asm	addpd	xmm4,xmm1	/* ~t0e */\
			__asm	subpd	xmm5,xmm0	/* ~t0f 	xmm6,7 free */\
			\
			__asm	movaps	xmm6,[edx     ]	/* c3m1*/\
			__asm	movaps	xmm7,[edx+0x10]	/* s3  */\
			__asm	movaps	xmm0,[eax     ]	/* t02 */\
			__asm	movaps	xmm1,[eax+0x10]	/* t03 */\
			\
		__asm	mov	esi, __o7\
			\
			__asm	subpd	xmm2,xmm4		/* t0e =rt-t0e */\
			__asm	subpd	xmm3,xmm5		/* t0f =it-t0f */\
			__asm	addpd	xmm4,xmm4		/*       2*t0e */\
			__asm	addpd	xmm5,xmm5		/*       2*t0f */\
			__asm	addpd	xmm4,xmm2		/* t08 =rt+t0e */\
			__asm	addpd	xmm5,xmm3		/* t09 =it+t0f */\
			__asm	addpd	xmm0,xmm4	/* t02 = t02+t08 */\
			__asm	addpd	xmm1,xmm5	/* t03 = t03+t09 */\
			__asm	movaps	[esi      ],xmm0	/* <- t02, put in slot 7 of 1,4,7 */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t03 */\
			\
			__asm	mulpd	xmm4,xmm6		/* t08 *= c3m1 */\
			__asm	mulpd	xmm5,xmm6		/* t09 *= c3m1 */\
			__asm	mulpd	xmm2,xmm7		/* rt = t0e*s3   */\
			__asm	mulpd	xmm3,xmm7		/* it = t0f*s3   */\
			__asm	addpd	xmm4,xmm0	/* t08 = t02+c3m1*t08 */\
			__asm	addpd	xmm5,xmm1	/* t09 = t03+c3m1*t09 */\
			\
		__asm	mov	edi, __o1\
		__asm	mov	esi, __o4\
			\
			__asm	movaps	xmm0,xmm4	/* copy t08 */\
			__asm	movaps	xmm1,xmm5	/* copy t09 */\
			\
			__asm	addpd	xmm4,xmm3	/* t08 = t08+it */\
			__asm	subpd	xmm5,xmm2	/* t09 = t09-rt */\
			__asm	subpd	xmm0,xmm3	/* t0e = t08-it */\
			__asm	addpd	xmm1,xmm2	/* t0f = t09+rt */\
			\
			__asm	movaps	[edi      ],xmm4	/* <- t08, put in slot 1 of 1,4,7 */\
			__asm	movaps	[edi+0x010],xmm5	/* <- t09 */\
			__asm	movaps	[esi      ],xmm0	/* <- t0e, put in slot 4 of 1,4,7 */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t0f */\
		/*\
			rt  =t0a*c2+t0b*s2;		it  =t0b*c2-t0a*s2;\
			tt  =t0g*c4+t0h*s4;		t0h =t0h*c4-t0g*s4;	t0g=tt;\
			t0a =rt+t0g;			t0b =it+t0h;\
			t0g =rt-t0g;			t0h =it-t0h;\
			t04 =t04+t0a;			t05 =t05+t0b;\
			A5r =t04;				A5i =t05;\
			t0a =t04+c3m1*t0a;		t0b =t05+c3m1*t0b;\
			rt  =s3*t0g;			it  =s3*t0h;\
			A8r =t0a+it;			A8i =t0b-rt;\
			A2r =t0a-it;			A2i =t0b+rt;\
		*/\
			__asm	add	eax, 0x20	/* t04 */\
			__asm	add	ebx, 0x20	/* t0a */\
			__asm	add	ecx, 0x20	/* t0g */\
			__asm	mov	edi, __c2\
			\
			__asm	movaps	xmm2,[ebx     ]	/* t08 */\
			__asm	movaps	xmm3,[ebx+0x10]	/* t09 */\
			__asm	movaps	xmm6,[edi     ]	/* c */\
			__asm	movaps	xmm7,[edi+0x10]	/* s */\
			__asm	movaps	xmm0,xmm2		/* copy of t08 */\
			__asm	movaps	xmm1,xmm3		/* copy of t09 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t08*c  */\
			__asm	mulpd	xmm3,xmm6		/* t09*c  */\
			__asm	mulpd	xmm0,xmm7		/* t08*s  */\
			__asm	mulpd	xmm1,xmm7		/* t09*s  */\
			__asm	addpd	xmm2,xmm1	/* xmm2 <- rt */\
			__asm	subpd	xmm3,xmm0	/* xmm3 <- it 	xmm6,7 free */\
			\
		__asm	add	edi, 0x40	/* c4 */\
			__asm	movaps	xmm4,[ecx     ]	/* t0e */\
			__asm	movaps	xmm5,[ecx+0x10]	/* t0f */\
			__asm	movaps	xmm6,[edi     ]	/* c2*/\
			__asm	movaps	xmm7,[edi+0x10]	/* s2*/\
			__asm	movaps	xmm0,xmm4		/* copy of t0e */\
			__asm	movaps	xmm1,xmm5		/* copy of t0f */\
			\
			__asm	mulpd	xmm4,xmm6		/* t0e*c2 */\
			__asm	mulpd	xmm5,xmm6		/* t0f*c2 */\
			__asm	mulpd	xmm0,xmm7		/* t0e*s2 */\
			__asm	mulpd	xmm1,xmm7		/* t0f*s2 */\
			__asm	addpd	xmm4,xmm1	/* ~t0e */\
			__asm	subpd	xmm5,xmm0	/* ~t0f 	xmm6,7 free */\
			\
			__asm	movaps	xmm6,[edx     ]	/* c3m1*/\
			__asm	movaps	xmm7,[edx+0x10]	/* s3  */\
			__asm	movaps	xmm0,[eax     ]	/* t04 */\
			__asm	movaps	xmm1,[eax+0x10]	/* t05 */\
			\
		__asm	mov	esi, __o5\
			\
			__asm	subpd	xmm2,xmm4		/* t0g =rt-t0g */\
			__asm	subpd	xmm3,xmm5		/* t0h =it-t0h */\
			__asm	addpd	xmm4,xmm4		/*       2*t0g */\
			__asm	addpd	xmm5,xmm5		/*       2*t0h */\
			__asm	addpd	xmm4,xmm2		/* t0a =rt+t0g */\
			__asm	addpd	xmm5,xmm3		/* t0b =it+t0h */\
			__asm	addpd	xmm0,xmm4	/* t04 = t04+t0a */\
			__asm	addpd	xmm1,xmm5	/* t05 = t05+t0b */\
			__asm	movaps	[esi      ],xmm0	/* <- t04, put in slot 5 of 2,5,8 */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t05 */\
			\
			__asm	mulpd	xmm4,xmm6		/* t0a *= c3m1 */\
			__asm	mulpd	xmm5,xmm6		/* t0b *= c3m1 */\
			__asm	mulpd	xmm2,xmm7		/* rt = t0g*s3   */\
			__asm	mulpd	xmm3,xmm7		/* it = t0h*s3   */\
			__asm	addpd	xmm4,xmm0	/* t0a = t04+c3m1*t0a */\
			__asm	addpd	xmm5,xmm1	/* t0b = t05+c3m1*t0b */\
			\
		__asm	mov	edi, __o8\
		__asm	mov	esi, __o2\
			\
			__asm	movaps	xmm0,xmm4	/* copy t0a */\
			__asm	movaps	xmm1,xmm5	/* copy t0b */\
			\
			__asm	addpd	xmm4,xmm3	/* t0a = t0a+it */\
			__asm	subpd	xmm5,xmm2	/* t0b = t0b-rt */\
			__asm	subpd	xmm0,xmm3	/* t0g = t0a-it */\
			__asm	addpd	xmm1,xmm2	/* t0h = t0b+rt */\
			\
			__asm	movaps	[edi      ],xmm4	/* <- t0a, put in slot 8 of 2,5,8 */\
			__asm	movaps	[edi+0x010],xmm5	/* <- t0b */\
			__asm	movaps	[esi      ],xmm0	/* <- t0g, put in slot 2 of 2,5,8 */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t0h */\
		}

		/*...Radix-9 DIT: Inputs enter in memory locations __i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8.\
				Outputs go into memory locations __o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8, which may coincide with inputs:\
		*/\
		#define SSE2_RADIX_09_DIT_B(__i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8, __c1,__c2,__c3m1,__c4, __o0,__o1,__o2,__o3,__o4,__o5,__o6,__o7,__o8)\
		{\
		/*...gather the needed data (9 64-b__it complex, i.e. 18 64-b__it reals) and do three radix-3 transforms: */\
		/*\
			rt  =t02;				it  =t03;\
			t02 =rt+t04;			t03 =it+t05;\
			t04 =rt-t04;			t05 =it-t05;\
			t00 =t00+t02;			t01 =t01+t03;\
			t02 =t00+c3m1*t02;	t03 =t01+c3m1*t03;\
			rt  =s3*t04;			it  =s3*t05;\
			t04 =t02-it;			t05 =t03+rt;\
			t02 =t02+it;			t03 =t03-rt;\
		*/\
			__asm	mov	eax, __i0\
			__asm	mov	ebx, __i1\
			__asm	mov	ecx, __i2\
			__asm	mov	edx, __c3m1\
			\
			__asm	movaps	xmm2,[ebx     ]	/* A3r */\
			__asm	movaps	xmm3,[ebx+0x10]	/* A3i */\
			__asm	movaps	xmm0,[eax     ]	/* A0r */\
			__asm	movaps	xmm1,[eax+0x10]	/* A0i */\
			__asm	movaps	xmm6,[ecx     ]	/* A6r */\
			__asm	movaps	xmm7,[ecx+0x10]	/* A6i */\
			__asm	movaps	xmm4,xmm2		/* t04: init to copy of A3r */\
			__asm	movaps	xmm5,xmm3		/* t05: init to copy of A3i */\
			\
			__asm	addpd	xmm2,xmm6		/* t02=A3r+A6r */\
			__asm	addpd	xmm3,xmm7		/* t03=A3i+A6i */\
			__asm	subpd	xmm4,xmm6		/* t04=A3r-A6r */\
			__asm	subpd	xmm5,xmm7		/* t05=A3i-A6i */\
			__asm	addpd	xmm0,xmm2	/* t00 = t00+t02 */\
			__asm	addpd	xmm1,xmm3	/* t01 = t01+t03 */\
			__asm	movaps	xmm6,[edx     ]	/* c3m1*/\
			__asm	movaps	xmm7,[edx+0x10]	/* s3  */\
			__asm	movaps	[eax      ],xmm0	/* <- t00 */\
			__asm	movaps	[eax+0x010],xmm1	/* <- t01 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t02 *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t03 *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = t04*s3   */\
			__asm	mulpd	xmm5,xmm7		/* it = t05*s3   */\
			__asm	addpd	xmm2,xmm0	/* t02 = t00+c3m1*t02 */\
			__asm	addpd	xmm3,xmm1	/* t03 = t01+c3m1*t03 */\
			\
			__asm	movaps	xmm0,xmm2	/* copy t02 */\
			__asm	movaps	xmm1,xmm3	/* copy t03 */\
			\
			__asm	addpd	xmm2,xmm5	/* t02 = t02-it */\
			__asm	subpd	xmm3,xmm4	/* t03 = t03+rt */\
			__asm	subpd	xmm0,xmm5	/* t04 = t02+it */\
			__asm	addpd	xmm1,xmm4	/* t05 = t03-rt */\
			\
			__asm	movaps	[ebx      ],xmm2	/* <- t02 */\
			__asm	movaps	[ebx+0x010],xmm3	/* <- t03 */\
			__asm	movaps	[ecx      ],xmm0	/* <- t04 */\
			__asm	movaps	[ecx+0x010],xmm1	/* <- t05 */\
			\
		/*\
			rt  =t08;				it  =t09;\
			t08 =rt+t0a;			t09 =it+t0b;\
			t0a =rt-t0a;			t0b =it-t0b;\
			t06 =t06+t08;			t07 =t07+t09;\
			t08 =t06+c3m1*t08;	t09 =t07+c3m1*t09;\
			rt  =s3*t0a;			it  =s3*t0b;\
			t0a =t08-it;			t0b =t09+rt;\
			t08 =t08+it;			t09 =t09-rt;\
		*/\
			__asm	mov	eax, __i3\
			__asm	mov	ebx, __i4\
			__asm	mov	ecx, __i5\
			\
			__asm	movaps	xmm4,[ebx     ]	/* A3r */\
			__asm	movaps	xmm5,[ebx+0x10]	/* A3i */\
			__asm	movaps	xmm2,[ecx     ]	/* A6r */\
			__asm	movaps	xmm3,[ecx+0x10]	/* A6i */\
			__asm	movaps	xmm0,[eax     ]	/* A0r */\
			__asm	movaps	xmm1,[eax+0x10]	/* A0i */\
			\
			__asm	subpd	xmm4,xmm2		/* t04=A3r-A6r */\
			__asm	subpd	xmm5,xmm3		/* t05=A3i-A6i */\
			__asm	addpd	xmm2,xmm2		/*       2*A6r */\
			__asm	addpd	xmm3,xmm3		/*       2*A6i */\
			__asm	addpd	xmm2,xmm4		/* t02=A3r+A6r */\
			__asm	addpd	xmm3,xmm5		/* t03=A3i+A6i */\
			__asm	addpd	xmm0,xmm2	/* t00 = t00+t02 */\
			__asm	addpd	xmm1,xmm3	/* t01 = t01+t03 */\
			__asm	movaps	[eax      ],xmm0	/* <- t00 */\
			__asm	movaps	[eax+0x010],xmm1	/* <- t01 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t02 *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t03 *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = t04*s3   */\
			__asm	mulpd	xmm5,xmm7		/* it = t05*s3   */\
			__asm	addpd	xmm2,xmm0	/* t02 = t00+c3m1*t02 */\
			__asm	addpd	xmm3,xmm1	/* t03 = t01+c3m1*t03 */\
			\
			__asm	movaps	xmm0,xmm2	/* copy t02 */\
			__asm	movaps	xmm1,xmm3	/* copy t03 */\
			\
			__asm	addpd	xmm2,xmm5	/* t02 = t02-it */\
			__asm	subpd	xmm3,xmm4	/* t03 = t03+rt */\
			__asm	subpd	xmm0,xmm5	/* t04 = t02+it */\
			__asm	addpd	xmm1,xmm4	/* t05 = t03-rt */\
			\
			__asm	movaps	[ebx      ],xmm2	/* <- t02 */\
			__asm	movaps	[ebx+0x010],xmm3	/* <- t03 */\
			__asm	movaps	[ecx      ],xmm0	/* <- t04 */\
			__asm	movaps	[ecx+0x010],xmm1	/* <- t05 */\
			\
		/*\
			rt  =t0e;				it  =t0f;\
			t0e =rt+t0g;			t0f =it+t0h;\
			t0g =rt-t0g;			t0h =it-t0h;\
			t0c =t0c+t0e;			t0d =t0d+t0f;\
			t0e =t0c+c3m1*t0e;	t0f =t0d+c3m1*t0f;\
			rt  =s3*t0g;			it  =s3*t0h;\
			t0g =t0e-it;			t0h =t0f+rt;\
			t0e =t0e+it;			t0f =t0f-rt;\
		*/\
			__asm	mov	eax, __i6\
			__asm	mov	ebx, __i7\
			__asm	mov	ecx, __i8\
			\
			__asm	movaps	xmm4,[ebx     ]	/* A3r */\
			__asm	movaps	xmm5,[ebx+0x10]	/* A3i */\
			__asm	movaps	xmm2,[ecx     ]	/* A6r */\
			__asm	movaps	xmm3,[ecx+0x10]	/* A6i */\
			__asm	movaps	xmm0,[eax     ]	/* A0r */\
			__asm	movaps	xmm1,[eax+0x10]	/* A0i */\
			\
			__asm	subpd	xmm4,xmm2		/* t04=A3r-A6r */\
			__asm	subpd	xmm5,xmm3		/* t05=A3i-A6i */\
			__asm	addpd	xmm2,xmm2		/*       2*A6r */\
			__asm	addpd	xmm3,xmm3		/*       2*A6i */\
			__asm	addpd	xmm2,xmm4		/* t02=A3r+A6r */\
			__asm	addpd	xmm3,xmm5		/* t03=A3i+A6i */\
			__asm	addpd	xmm0,xmm2	/* t00 = t00+t02 */\
			__asm	addpd	xmm1,xmm3	/* t01 = t01+t03 */\
			__asm	movaps	[eax      ],xmm0	/* <- t00 */\
			__asm	movaps	[eax+0x010],xmm1	/* <- t01 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t02 *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t03 *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = t04*s3   */\
			__asm	mulpd	xmm5,xmm7		/* it = t05*s3   */\
			__asm	addpd	xmm2,xmm0	/* t02 = t00+c3m1*t02 */\
			__asm	addpd	xmm3,xmm1	/* t03 = t01+c3m1*t03 */\
			\
			__asm	movaps	xmm0,xmm2	/* copy t02 */\
			__asm	movaps	xmm1,xmm3	/* copy t03 */\
			\
			__asm	addpd	xmm2,xmm5	/* t02 = t02-it */\
			__asm	subpd	xmm3,xmm4	/* t03 = t03+rt */\
			__asm	subpd	xmm0,xmm5	/* t04 = t02+it */\
			__asm	addpd	xmm1,xmm4	/* t05 = t03-rt */\
			\
			__asm	movaps	[ebx      ],xmm2	/* <- t02 */\
			__asm	movaps	[ebx+0x010],xmm3	/* <- t03 */\
			__asm	movaps	[ecx      ],xmm0	/* <- t04 */\
			__asm	movaps	[ecx+0x010],xmm1	/* <- t05 */\
			\
		/******************************************************************************/\
		/*...and now do three more radix-3 transforms, including the twiddle factors: */\
		/******************************************************************************/\
		/*\
			rt  =t06;				it  =t07;\
			t06 =rt+t0c;			t07 =it+t0d;\
			t0c =rt-t0c;			t0d =it-t0d;\
			t00 =t00+t06;			t01 =t01+t07;\
			A0r =t00;				A0i =t01;\
			t06 =t00+c3m1*t06;		t07 =t01+c3m1*t07;\
			rt  =s3*t0c;			it  =s3*t0d;\
			A3r =t06+it;			A3i =t07-rt;\
			A6r =t06-it;			A6i =t07+rt;\
		*/\
			__asm	mov	eax, __i0\
			__asm	mov	ebx, __i3\
			__asm	mov	ecx, __i6\
			\
		__asm	mov	esi, __o0\
			\
			__asm	movaps	xmm4,[ebx     ]	/* A3r */\
			__asm	movaps	xmm5,[ebx+0x10]	/* A3i */\
			__asm	movaps	xmm2,[ecx     ]	/* A6r */\
			__asm	movaps	xmm3,[ecx+0x10]	/* A6i */\
			__asm	movaps	xmm0,[eax     ]	/* A0r */\
			__asm	movaps	xmm1,[eax+0x10]	/* A0i */\
			\
			__asm	subpd	xmm4,xmm2		/* t04=A3r-A6r */\
			__asm	subpd	xmm5,xmm3		/* t05=A3i-A6i */\
			__asm	addpd	xmm2,xmm2		/*       2*A6r */\
			__asm	addpd	xmm3,xmm3		/*       2*A6i */\
			__asm	addpd	xmm2,xmm4		/* t02=A3r+A6r */\
			__asm	addpd	xmm3,xmm5		/* t03=A3i+A6i */\
			__asm	addpd	xmm0,xmm2	/* t00 = t00+t02 */\
			__asm	addpd	xmm1,xmm3	/* t01 = t01+t03 */\
			__asm	movaps	[esi      ],xmm0	/* <- t00 */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t01 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t02 *= c3m1 */\
			__asm	mulpd	xmm3,xmm6		/* t03 *= c3m1 */\
			__asm	mulpd	xmm4,xmm7		/* rt = t04*s3   */\
			__asm	mulpd	xmm5,xmm7		/* it = t05*s3   */\
			__asm	addpd	xmm2,xmm0	/* t02 = t00+c3m1*t02 */\
			__asm	addpd	xmm3,xmm1	/* t03 = t01+c3m1*t03 */\
			\
			__asm	movaps	xmm0,xmm2	/* copy t02 */\
			__asm	movaps	xmm1,xmm3	/* copy t03 */\
			\
		__asm	mov	edi, __o3\
		__asm	mov	esi, __o6\
			\
			__asm	addpd	xmm2,xmm5	/* t06 = t06+it */\
			__asm	subpd	xmm3,xmm4	/* t07 = t07-rt */\
			__asm	subpd	xmm0,xmm5	/* t0c = t06-it */\
			__asm	addpd	xmm1,xmm4	/* t0d = t07+rt */\
			\
			__asm	movaps	[edi      ],xmm2	/* <- t06 */\
			__asm	movaps	[edi+0x010],xmm3	/* <- t07 */\
			__asm	movaps	[esi      ],xmm0	/* <- t0c */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t0d */\
			\
		/*\
			rt  =t08*c +t09*s;		it  =t09*c -t08*s;\
			tt  =t0e*c2+t0f*s2;		t0f =t0f*c2-t0e*s2;	t0e=tt;\
			t08 =rt+t0e;			t09 =it+t0f;\
			t0e =rt-t0e;			t0f =it-t0f;\
			t02 =t02+t08;			t03 =t03+t09;\
			A7r =t02;				A7i =t03;\
			t08 =t02+c3m1*t08;		t09 =t03+c3m1*t09;\
			rt  =s3*t0e;			it  =s3*t0f;\
			A1r =t08+it;			A1i =t09-rt;\
			A4r =t08-it;			A4i =t09+rt;\
		*/\
			__asm	mov	eax, __i1\
			__asm	mov	ebx, __i4\
			__asm	mov	ecx, __i7\
			__asm	mov	edi, __c1\
			\
			__asm	movaps	xmm2,[ebx     ]	/* t08 */\
			__asm	movaps	xmm3,[ebx+0x10]	/* t09 */\
			__asm	movaps	xmm6,[edi     ]	/* c */\
			__asm	movaps	xmm7,[edi+0x10]	/* s */\
			__asm	movaps	xmm0,xmm2		/* copy of t08 */\
			__asm	movaps	xmm1,xmm3		/* copy of t09 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t08*c  */\
			__asm	mulpd	xmm3,xmm6		/* t09*c  */\
			__asm	mulpd	xmm0,xmm7		/* t08*s  */\
			__asm	mulpd	xmm1,xmm7		/* t09*s  */\
			__asm	addpd	xmm2,xmm1	/* xmm2 <- rt */\
			__asm	subpd	xmm3,xmm0	/* xmm3 <- it 	xmm6,7 free */\
			\
		__asm	add	edi, 0x20	/* c2 */\
			__asm	movaps	xmm4,[ecx     ]	/* t0e */\
			__asm	movaps	xmm5,[ecx+0x10]	/* t0f */\
			__asm	movaps	xmm6,[edi     ]	/* c2*/\
			__asm	movaps	xmm7,[edi+0x10]	/* s2*/\
			__asm	movaps	xmm0,xmm4		/* copy of t0e */\
			__asm	movaps	xmm1,xmm5		/* copy of t0f */\
			\
			__asm	mulpd	xmm4,xmm6		/* t0e*c2 */\
			__asm	mulpd	xmm5,xmm6		/* t0f*c2 */\
			__asm	mulpd	xmm0,xmm7		/* t0e*s2 */\
			__asm	mulpd	xmm1,xmm7		/* t0f*s2 */\
			__asm	addpd	xmm4,xmm1	/* ~t0e */\
			__asm	subpd	xmm5,xmm0	/* ~t0f 	xmm6,7 free */\
			\
			__asm	movaps	xmm6,[edx     ]	/* c3m1*/\
			__asm	movaps	xmm7,[edx+0x10]	/* s3  */\
			__asm	movaps	xmm0,[eax     ]	/* t02 */\
			__asm	movaps	xmm1,[eax+0x10]	/* t03 */\
			\
		__asm	mov	esi, __o7\
			\
			__asm	subpd	xmm2,xmm4		/* t0e =rt-t0e */\
			__asm	subpd	xmm3,xmm5		/* t0f =it-t0f */\
			__asm	addpd	xmm4,xmm4		/*       2*t0e */\
			__asm	addpd	xmm5,xmm5		/*       2*t0f */\
			__asm	addpd	xmm4,xmm2		/* t08 =rt+t0e */\
			__asm	addpd	xmm5,xmm3		/* t09 =it+t0f */\
			__asm	addpd	xmm0,xmm4	/* t02 = t02+t08 */\
			__asm	addpd	xmm1,xmm5	/* t03 = t03+t09 */\
			__asm	movaps	[esi      ],xmm0	/* <- t02, put in slot 7 of 1,4,7 */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t03 */\
			\
			__asm	mulpd	xmm4,xmm6		/* t08 *= c3m1 */\
			__asm	mulpd	xmm5,xmm6		/* t09 *= c3m1 */\
			__asm	mulpd	xmm2,xmm7		/* rt = t0e*s3   */\
			__asm	mulpd	xmm3,xmm7		/* it = t0f*s3   */\
			__asm	addpd	xmm4,xmm0	/* t08 = t02+c3m1*t08 */\
			__asm	addpd	xmm5,xmm1	/* t09 = t03+c3m1*t09 */\
			\
		__asm	mov	edi, __o1\
		__asm	mov	esi, __o4\
			\
			__asm	movaps	xmm0,xmm4	/* copy t08 */\
			__asm	movaps	xmm1,xmm5	/* copy t09 */\
			\
			__asm	addpd	xmm4,xmm3	/* t08 = t08+it */\
			__asm	subpd	xmm5,xmm2	/* t09 = t09-rt */\
			__asm	subpd	xmm0,xmm3	/* t0e = t08-it */\
			__asm	addpd	xmm1,xmm2	/* t0f = t09+rt */\
			\
			__asm	movaps	[edi      ],xmm4	/* <- t08, put in slot 1 of 1,4,7 */\
			__asm	movaps	[edi+0x010],xmm5	/* <- t09 */\
			__asm	movaps	[esi      ],xmm0	/* <- t0e, put in slot 4 of 1,4,7 */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t0f */\
		/*\
			rt  =t0a*c2+t0b*s2;		it  =t0b*c2-t0a*s2;\
			tt  =t0g*c4+t0h*s4;		t0h =t0h*c4-t0g*s4;	t0g=tt;\
			t0a =rt+t0g;			t0b =it+t0h;\
			t0g =rt-t0g;			t0h =it-t0h;\
			t04 =t04+t0a;			t05 =t05+t0b;\
			A5r =t04;				A5i =t05;\
			t0a =t04+c3m1*t0a;		t0b =t05+c3m1*t0b;\
			rt  =s3*t0g;			it  =s3*t0h;\
			A8r =t0a+it;			A8i =t0b-rt;\
			A2r =t0a-it;			A2i =t0b+rt;\
		*/\
			__asm	mov	eax, __i2\
			__asm	mov	ebx, __i5\
			__asm	mov	ecx, __i8\
			__asm	mov	edi, __c2\
			\
			__asm	movaps	xmm2,[ebx     ]	/* t08 */\
			__asm	movaps	xmm3,[ebx+0x10]	/* t09 */\
			__asm	movaps	xmm6,[edi     ]	/* c */\
			__asm	movaps	xmm7,[edi+0x10]	/* s */\
			__asm	movaps	xmm0,xmm2		/* copy of t08 */\
			__asm	movaps	xmm1,xmm3		/* copy of t09 */\
			\
			__asm	mulpd	xmm2,xmm6		/* t08*c  */\
			__asm	mulpd	xmm3,xmm6		/* t09*c  */\
			__asm	mulpd	xmm0,xmm7		/* t08*s  */\
			__asm	mulpd	xmm1,xmm7		/* t09*s  */\
			__asm	addpd	xmm2,xmm1	/* xmm2 <- rt */\
			__asm	subpd	xmm3,xmm0	/* xmm3 <- it 	xmm6,7 free */\
			\
		__asm	add	edi, 0x40	/* c4 */\
			__asm	movaps	xmm4,[ecx     ]	/* t0e */\
			__asm	movaps	xmm5,[ecx+0x10]	/* t0f */\
			__asm	movaps	xmm6,[edi     ]	/* c2*/\
			__asm	movaps	xmm7,[edi+0x10]	/* s2*/\
			__asm	movaps	xmm0,xmm4		/* copy of t0e */\
			__asm	movaps	xmm1,xmm5		/* copy of t0f */\
			\
			__asm	mulpd	xmm4,xmm6		/* t0e*c2 */\
			__asm	mulpd	xmm5,xmm6		/* t0f*c2 */\
			__asm	mulpd	xmm0,xmm7		/* t0e*s2 */\
			__asm	mulpd	xmm1,xmm7		/* t0f*s2 */\
			__asm	addpd	xmm4,xmm1	/* ~t0e */\
			__asm	subpd	xmm5,xmm0	/* ~t0f 	xmm6,7 free */\
			\
			__asm	movaps	xmm6,[edx     ]	/* c3m1*/\
			__asm	movaps	xmm7,[edx+0x10]	/* s3  */\
			__asm	movaps	xmm0,[eax     ]	/* t04 */\
			__asm	movaps	xmm1,[eax+0x10]	/* t05 */\
			\
		__asm	mov	esi, __o5\
			\
			__asm	subpd	xmm2,xmm4		/* t0g =rt-t0g */\
			__asm	subpd	xmm3,xmm5		/* t0h =it-t0h */\
			__asm	addpd	xmm4,xmm4		/*       2*t0g */\
			__asm	addpd	xmm5,xmm5		/*       2*t0h */\
			__asm	addpd	xmm4,xmm2		/* t0a =rt+t0g */\
			__asm	addpd	xmm5,xmm3		/* t0b =it+t0h */\
			__asm	addpd	xmm0,xmm4	/* t04 = t04+t0a */\
			__asm	addpd	xmm1,xmm5	/* t05 = t05+t0b */\
			__asm	movaps	[esi      ],xmm0	/* <- t04, put in slot 5 of 2,5,8 */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t05 */\
			\
			__asm	mulpd	xmm4,xmm6		/* t0a *= c3m1 */\
			__asm	mulpd	xmm5,xmm6		/* t0b *= c3m1 */\
			__asm	mulpd	xmm2,xmm7		/* rt = t0g*s3   */\
			__asm	mulpd	xmm3,xmm7		/* it = t0h*s3   */\
			__asm	addpd	xmm4,xmm0	/* t0a = t04+c3m1*t0a */\
			__asm	addpd	xmm5,xmm1	/* t0b = t05+c3m1*t0b */\
			\
		__asm	mov	edi, __o8\
		__asm	mov	esi, __o2\
			\
			__asm	movaps	xmm0,xmm4	/* copy t0a */\
			__asm	movaps	xmm1,xmm5	/* copy t0b */\
			\
			__asm	addpd	xmm4,xmm3	/* t0a = t0a+it */\
			__asm	subpd	xmm5,xmm2	/* t0b = t0b-rt */\
			__asm	subpd	xmm0,xmm3	/* t0g = t0a-it */\
			__asm	addpd	xmm1,xmm2	/* t0h = t0b+rt */\
			\
			__asm	movaps	[edi      ],xmm4	/* <- t0a, put in slot 8 of 2,5,8 */\
			__asm	movaps	[edi+0x010],xmm5	/* <- t0b */\
			__asm	movaps	[esi      ],xmm0	/* <- t0g, put in slot 2 of 2,5,8 */\
			__asm	movaps	[esi+0x010],xmm1	/* <- t0h */\
		}

	#else	/* GCC-style inline ASM: */

	  #if OS_BITS == 32

			/*...Radix-9 DIF: Outs adjacent in memory starting at address __o0, i.e. reals separated by 0x20 bytes.\
				Inputs enter in memory locations __i0,__i1,__i2,__i3,__i4,__i5,__i6,__i7,__i8, assumed disjoint with outputs:\
			*/\
			#define SSE2_RADIX_09_DIF(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8)\
			{\
			__asm__ volatile (\
			"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
			"movl	%[__i0]	,%%eax 	\n\t"/* __i0-8; e[abc]x store input addresses */\
			"movl	%[__o0]	,%%esi	\n\t"\
			"movl	%[__c1]	,%%edx 	\n\t"/* edx stores trig addresses throughout */\
			/* Block 1: */\
				"movl	%%eax		,%%ebx\n\t"\
				"movl	%%eax		,%%ecx\n\t"\
				"movl	%[__i3]	,%%ebx 	\n\t"\
				"movl	%[__i6]	,%%ecx 	\n\t"\
				"movaps	    (%%ebx)	,%%xmm2\n\t"\
				"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
				"movaps	    (%%eax)	,%%xmm0\n\t"\
				"movaps	0x10(%%eax)	,%%xmm1\n\t"\
				"movaps	    (%%ecx)	,%%xmm6\n\t"\
				"movaps	0x10(%%ecx)	,%%xmm7\n\t"\
				"movaps	%%xmm2		,%%xmm4\n\t"\
				"movaps	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm6		,%%xmm2\n\t"\
				"addpd	%%xmm7		,%%xmm3\n\t"\
				"subpd	%%xmm6		,%%xmm4\n\t"\
				"subpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm0\n\t"\
				"addpd	%%xmm3		,%%xmm1\n\t"\
				"movaps	0x40(%%edx)	,%%xmm6\n\t"/* c3m1 */\
				"movaps	0x50(%%edx)	,%%xmm7\n\t"\
				"movaps	%%xmm0		,    (%%esi)	/* <- t00 */\n\t"\
				"movaps	%%xmm1		,0x10(%%esi)	/* <- t01 */\n\t"\
			"movl	%[__o1]	,%%edi	\n\t"\
			"movl	%[__o2]	,%%esi	\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm4\n\t"\
				"mulpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm0		,%%xmm2\n\t"\
				"addpd	%%xmm1		,%%xmm3\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"subpd	%%xmm5		,%%xmm2\n\t"\
				"addpd	%%xmm4		,%%xmm3\n\t"\
				"addpd	%%xmm5		,%%xmm0\n\t"\
				"subpd	%%xmm4		,%%xmm1\n\t"\
				"movaps	%%xmm2		,    (%%edi)	/* <- t02 */\n\t"\
				"movaps	%%xmm3		,0x10(%%edi)	/* <- t03 */\n\t"\
				"movaps	%%xmm0		,    (%%esi)	/* <- t04 */\n\t"\
				"movaps	%%xmm1		,0x10(%%esi)	/* <- t05 */\n\t"\
			/* Block 2: */\
				"movl	%[__i1]	,%%eax 	\n\t"\
				"movl	%[__i4]	,%%ebx 	\n\t"\
				"movl	%[__i7]	,%%ecx 	\n\t"\
				"movaps	    (%%ebx)	,%%xmm4\n\t"\
				"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
				"movaps	    (%%eax)	,%%xmm0\n\t"\
				"movaps	0x10(%%eax)	,%%xmm1\n\t"\
				"movaps	    (%%ecx)	,%%xmm2\n\t"\
				"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movl	%[__o3]	,%%edi	\n\t"\
				"subpd	%%xmm2		,%%xmm4\n\t"\
				"subpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm2\n\t"\
				"addpd	%%xmm3		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm2\n\t"\
				"addpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm2		,%%xmm0\n\t"\
				"addpd	%%xmm3		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%edi)	/* <- t06 */\n\t"\
				"movaps	%%xmm1		,0x10(%%edi)	/* <- t07 */\n\t"\
			"movl	%[__o4]	,%%esi	\n\t"\
			"movl	%[__o5]	,%%edi	\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm4\n\t"\
				"mulpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm0		,%%xmm2\n\t"\
				"addpd	%%xmm1		,%%xmm3\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"subpd	%%xmm5		,%%xmm2\n\t"\
				"addpd	%%xmm4		,%%xmm3\n\t"\
				"addpd	%%xmm5		,%%xmm0\n\t"\
				"subpd	%%xmm4		,%%xmm1\n\t"\
				"movaps	%%xmm2		,    (%%esi)	/* <- t08 */\n\t"\
				"movaps	%%xmm3		,0x10(%%esi)	/* <- t09 */\n\t"\
				"movaps	%%xmm0		,    (%%edi)	/* <- t0a */\n\t"\
				"movaps	%%xmm1		,0x10(%%edi)	/* <- t0b */\n\t"\
			/* Block 3: */\
				"movl	%[__i2]	,%%eax 	\n\t"\
				"movl	%[__i5]	,%%ebx 	\n\t"\
				"movl	%[__i8]	,%%ecx 	\n\t"\
				"movaps	    (%%ebx)	,%%xmm4\n\t"\
				"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
				"movaps	    (%%eax)	,%%xmm0\n\t"\
				"movaps	0x10(%%eax)	,%%xmm1\n\t"\
				"movaps	    (%%ecx)	,%%xmm2\n\t"\
				"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
			"movl	%[__o6]	,%%esi	\n\t"\
				"subpd	%%xmm2		,%%xmm4\n\t"\
				"subpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm2\n\t"\
				"addpd	%%xmm3		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm2\n\t"\
				"addpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm2		,%%xmm0\n\t"\
				"addpd	%%xmm3		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%esi)	/* <- t0c */\n\t"\
				"movaps	%%xmm1		,0x10(%%esi)	/* <- t0d */\n\t"\
			"movl	%[__o7]	,%%edi	\n\t"\
			"movl	%[__o8]	,%%esi	\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm4\n\t"\
				"mulpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm0		,%%xmm2\n\t"\
				"addpd	%%xmm1		,%%xmm3\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"subpd	%%xmm5		,%%xmm2\n\t"\
				"addpd	%%xmm4		,%%xmm3\n\t"\
				"addpd	%%xmm5		,%%xmm0\n\t"\
				"subpd	%%xmm4		,%%xmm1\n\t"\
				"movaps	%%xmm2		,    (%%edi)	/* <- t0e */\n\t"\
				"movaps	%%xmm3		,0x10(%%edi)	/* <- t0f */\n\t"\
				"movaps	%%xmm0		,    (%%esi)	/* <- t0g */\n\t"\
				"movaps	%%xmm1		,0x10(%%esi)	/* <- t0h */\n\t"\
			/*****************************************/\
			/* now do three more radix-3 transforms: */\
			/*****************************************/\
			/* Block 1: */\
			"movl	%[__o0]	,%%eax	\n\t"\
			"movl	%[__o3]	,%%ebx	\n\t"\
			"movl	%[__o6]	,%%ecx	\n\t"\
				"movaps	    (%%ebx)	,%%xmm4\n\t"\
				"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
				"movaps	    (%%ecx)	,%%xmm2\n\t"\
				"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
				"movaps	    (%%eax)	,%%xmm0\n\t"\
				"movaps	0x10(%%eax)	,%%xmm1\n\t"\
				"subpd	%%xmm2		,%%xmm4\n\t"\
				"subpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm2\n\t"\
				"addpd	%%xmm3		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm2\n\t"\
				"addpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm2		,%%xmm0\n\t"\
				"addpd	%%xmm3		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%eax)\n\t"\
				"movaps	%%xmm1		,0x10(%%eax)\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm4\n\t"\
				"mulpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm0		,%%xmm2\n\t"\
				"addpd	%%xmm1		,%%xmm3\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"subpd	%%xmm5		,%%xmm2\n\t"\
				"addpd	%%xmm4		,%%xmm3\n\t"\
				"addpd	%%xmm5		,%%xmm0\n\t"\
				"subpd	%%xmm4		,%%xmm1\n\t"\
				"movaps	%%xmm2		,    (%%ebx)\n\t"\
				"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
				"movaps	%%xmm0		,    (%%ecx)\n\t"\
				"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			/* Block 2: */\
			"movl	%[__o1]	,%%eax	\n\t"\
			"movl	%[__o4]	,%%ebx	\n\t"\
			"movl	%[__o7]	,%%ecx	\n\t"\
				"movaps	    (%%ebx)	,%%xmm2\n\t"\
				"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
				"movaps	    (%%edx)	,%%xmm6\n\t"/* c1 */\
				"movaps	0x10(%%edx)	,%%xmm7\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm0\n\t"\
				"mulpd	%%xmm7		,%%xmm1\n\t"\
				"subpd	%%xmm1		,%%xmm2\n\t"\
				"addpd	%%xmm0		,%%xmm3\n\t"\
				"movaps	    (%%ecx)	,%%xmm4\n\t"\
				"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
				"movaps	0x20(%%edx)	,%%xmm6\n\t"/* c2 */\
				"movaps	0x30(%%edx)	,%%xmm7\n\t"\
				"movaps	%%xmm4		,%%xmm0\n\t"\
				"movaps	%%xmm5		,%%xmm1\n\t"\
				"mulpd	%%xmm6		,%%xmm4\n\t"\
				"mulpd	%%xmm6		,%%xmm5\n\t"\
				"mulpd	%%xmm7		,%%xmm0\n\t"\
				"mulpd	%%xmm7		,%%xmm1\n\t"\
				"subpd	%%xmm1		,%%xmm4\n\t"\
				"addpd	%%xmm0		,%%xmm5\n\t"\
				"movaps	0x40(%%edx)	,%%xmm6\n\t"/* c3m1 */\
				"movaps	0x50(%%edx)	,%%xmm7\n\t"\
				"movaps	    (%%eax)	,%%xmm0\n\t"\
				"movaps	0x10(%%eax)	,%%xmm1\n\t"\
				"subpd	%%xmm4		,%%xmm2\n\t"\
				"subpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm4\n\t"\
				"addpd	%%xmm5		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm4\n\t"\
				"addpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm4		,%%xmm0\n\t"\
				"addpd	%%xmm5		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%eax)\n\t"\
				"movaps	%%xmm1		,0x10(%%eax)\n\t"\
				"mulpd	%%xmm6		,%%xmm4\n\t"\
				"mulpd	%%xmm6		,%%xmm5\n\t"\
				"mulpd	%%xmm7		,%%xmm2\n\t"\
				"mulpd	%%xmm7		,%%xmm3\n\t"\
				"addpd	%%xmm0		,%%xmm4\n\t"\
				"addpd	%%xmm1		,%%xmm5\n\t"\
				"movaps	%%xmm4		,%%xmm0\n\t"\
				"movaps	%%xmm5		,%%xmm1\n\t"\
				"subpd	%%xmm3		,%%xmm4\n\t"\
				"addpd	%%xmm2		,%%xmm5\n\t"\
				"addpd	%%xmm3		,%%xmm0\n\t"\
				"subpd	%%xmm2		,%%xmm1\n\t"\
				"movaps	%%xmm4		,    (%%ebx)\n\t"\
				"movaps	%%xmm5		,0x10(%%ebx)\n\t"\
				"movaps	%%xmm0		,    (%%ecx)\n\t"\
				"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			/* Block 3: */\
			"movl	%[__o2]	,%%eax	\n\t"\
			"movl	%[__o5]	,%%ebx	\n\t"\
			"movl	%[__o8]	,%%ecx	\n\t"\
				"movaps	    (%%ebx)	,%%xmm2\n\t"\
				"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
				"movaps	0x20(%%edx)	,%%xmm6\n\t"/* c2 */\
				"movaps	0x30(%%edx)	,%%xmm7\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm0\n\t"\
				"mulpd	%%xmm7		,%%xmm1\n\t"\
				"subpd	%%xmm1		,%%xmm2\n\t"\
				"addpd	%%xmm0		,%%xmm3\n\t"\
				"movaps	    (%%ecx)	,%%xmm4\n\t"\
				"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
				"movaps	0x60(%%edx)	,%%xmm6\n\t"/* c4 */\
				"movaps	0x70(%%edx)	,%%xmm7\n\t"\
				"movaps	%%xmm4		,%%xmm0\n\t"\
				"movaps	%%xmm5		,%%xmm1\n\t"\
				"mulpd	%%xmm6		,%%xmm4\n\t"\
				"mulpd	%%xmm6		,%%xmm5\n\t"\
				"mulpd	%%xmm7		,%%xmm0\n\t"\
				"mulpd	%%xmm7		,%%xmm1\n\t"\
				"subpd	%%xmm1		,%%xmm4\n\t"\
				"addpd	%%xmm0		,%%xmm5\n\t"\
				"movaps	0x40(%%edx)	,%%xmm6\n\t"/* c3m1 */\
				"movaps	0x50(%%edx)	,%%xmm7\n\t"\
				"movaps	    (%%eax)	,%%xmm0\n\t"\
				"movaps	0x10(%%eax)	,%%xmm1\n\t"\
				"\n\t"\
				"subpd	%%xmm4		,%%xmm2\n\t"\
				"subpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm4\n\t"\
				"addpd	%%xmm5		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm4\n\t"\
				"addpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm4		,%%xmm0\n\t"\
				"addpd	%%xmm5		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%eax)\n\t"\
				"movaps	%%xmm1		,0x10(%%eax)\n\t"\
				"mulpd	%%xmm6		,%%xmm4\n\t"\
				"mulpd	%%xmm6		,%%xmm5\n\t"\
				"mulpd	%%xmm7		,%%xmm2\n\t"\
				"mulpd	%%xmm7		,%%xmm3\n\t"\
				"addpd	%%xmm0		,%%xmm4\n\t"\
				"addpd	%%xmm1		,%%xmm5\n\t"\
				"\n\t"\
				"movaps	%%xmm4		,%%xmm0\n\t"\
				"movaps	%%xmm5		,%%xmm1\n\t"\
				"subpd	%%xmm3		,%%xmm4\n\t"\
				"addpd	%%xmm2		,%%xmm5\n\t"\
				"addpd	%%xmm3		,%%xmm0\n\t"\
				"subpd	%%xmm2		,%%xmm1\n\t"\
				"movaps	%%xmm4		,    (%%ebx)\n\t"\
				"movaps	%%xmm5		,0x10(%%ebx)\n\t"\
				"movaps	%%xmm0		,    (%%ecx)\n\t"\
				"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			"popl %%ebx	\n\t"\
				:					/* outputs: none */\
				: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
				 ,[__i1] "m" (Xi1)\
				 ,[__i2] "m" (Xi2)\
				 ,[__i3] "m" (Xi3)\
				 ,[__i4] "m" (Xi4)\
				 ,[__i5] "m" (Xi5)\
				 ,[__i6] "m" (Xi6)\
				 ,[__i7] "m" (Xi7)\
				 ,[__i8] "m" (Xi8)\
				 ,[__c1] "m" (Xc1)\
				 ,[__o0] "m" (Xo0)\
				 ,[__o1] "m" (Xo1)\
				 ,[__o2] "m" (Xo2)\
				 ,[__o3] "m" (Xo3)\
				 ,[__o4] "m" (Xo4)\
				 ,[__o5] "m" (Xo5)\
				 ,[__o6] "m" (Xo6)\
				 ,[__o7] "m" (Xo7)\
				 ,[__o8] "m" (Xo8)\
				: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
				);\
			}

			/*...Radix-9 DIT: Ins in memory locations __i0-8.\
				Outs at memory locations __o0-8, assumed disjoint with inputs:\
			*/\
			#define SSE2_RADIX_09_DIT(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8)\
			{\
			__asm__ volatile (\
			"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
			"movl	%[__i0]	,%%eax\n\t"\
			"movl	%[__i1]	,%%ebx\n\t"\
			"movl	%[__i2]	,%%ecx\n\t"\
			"movl	%[__c1]	,%%edx 	\n\t"/* edx stores trig addresses throughout */\
			/* Block 1: */\
				"movaps	    (%%ebx)	,%%xmm2\n\t"\
				"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
				"movaps	    (%%eax)	,%%xmm0\n\t"\
				"movaps	0x10(%%eax)	,%%xmm1\n\t"\
				"movaps	    (%%ecx)	,%%xmm6\n\t"\
				"movaps	0x10(%%ecx)	,%%xmm7\n\t"\
				"movaps	%%xmm2		,%%xmm4\n\t"\
				"movaps	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm6		,%%xmm2\n\t"\
				"addpd	%%xmm7		,%%xmm3\n\t"\
				"subpd	%%xmm6		,%%xmm4\n\t"\
				"subpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm0\n\t"\
				"addpd	%%xmm3		,%%xmm1\n\t"\
				"movaps	0x40(%%edx)	,%%xmm6\n\t"/* c3m1 */\
				"movaps	0x50(%%edx)	,%%xmm7\n\t"\
				"movaps	%%xmm0		,    (%%eax)\n\t"\
				"movaps	%%xmm1		,0x10(%%eax)\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm4\n\t"\
				"mulpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm0		,%%xmm2\n\t"\
				"addpd	%%xmm1		,%%xmm3\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"addpd	%%xmm5		,%%xmm2\n\t"\
				"subpd	%%xmm4		,%%xmm3\n\t"\
				"subpd	%%xmm5		,%%xmm0\n\t"\
				"addpd	%%xmm4		,%%xmm1\n\t"\
				"movaps	%%xmm2		,    (%%ebx)\n\t"\
				"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
				"movaps	%%xmm0		,    (%%ecx)\n\t"\
				"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			/* Block 2: */\
			"movl	%[__i3]	,%%eax\n\t"\
			"movl	%[__i4]	,%%ebx\n\t"\
			"movl	%[__i5]	,%%ecx\n\t"\
				"movaps	    (%%ebx)	,%%xmm4\n\t"\
				"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
				"movaps	    (%%ecx)	,%%xmm2\n\t"\
				"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
				"movaps	    (%%eax)	,%%xmm0\n\t"\
				"movaps	0x10(%%eax)	,%%xmm1\n\t"\
				"subpd	%%xmm2		,%%xmm4\n\t"\
				"subpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm2\n\t"\
				"addpd	%%xmm3		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm2\n\t"\
				"addpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm2		,%%xmm0\n\t"\
				"addpd	%%xmm3		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%eax)\n\t"\
				"movaps	%%xmm1		,0x10(%%eax)\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm4\n\t"\
				"mulpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm0		,%%xmm2\n\t"\
				"addpd	%%xmm1		,%%xmm3\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"addpd	%%xmm5		,%%xmm2\n\t"\
				"subpd	%%xmm4		,%%xmm3\n\t"\
				"subpd	%%xmm5		,%%xmm0\n\t"\
				"addpd	%%xmm4		,%%xmm1\n\t"\
				"movaps	%%xmm2		,    (%%ebx)\n\t"\
				"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
				"movaps	%%xmm0		,    (%%ecx)\n\t"\
				"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			/* Block 3: */\
			"movl	%[__i6]	,%%eax\n\t"\
			"movl	%[__i7]	,%%ebx\n\t"\
			"movl	%[__i8]	,%%ecx\n\t"\
				"movaps	    (%%ebx)	,%%xmm4\n\t"\
				"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
				"movaps	    (%%ecx)	,%%xmm2\n\t"\
				"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
				"movaps	    (%%eax)	,%%xmm0\n\t"\
				"movaps	0x10(%%eax)	,%%xmm1\n\t"\
				"subpd	%%xmm2		,%%xmm4\n\t"\
				"subpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm2\n\t"\
				"addpd	%%xmm3		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm2\n\t"\
				"addpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm2		,%%xmm0\n\t"\
				"addpd	%%xmm3		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%eax)\n\t"\
				"movaps	%%xmm1		,0x10(%%eax)\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm4\n\t"\
				"mulpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm0		,%%xmm2\n\t"\
				"addpd	%%xmm1		,%%xmm3\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"addpd	%%xmm5		,%%xmm2\n\t"\
				"subpd	%%xmm4		,%%xmm3\n\t"\
				"subpd	%%xmm5		,%%xmm0\n\t"\
				"addpd	%%xmm4		,%%xmm1\n\t"\
				"movaps	%%xmm2		,    (%%ebx)\n\t"\
				"movaps	%%xmm3		,0x10(%%ebx)\n\t"\
				"movaps	%%xmm0		,    (%%ecx)\n\t"\
				"movaps	%%xmm1		,0x10(%%ecx)\n\t"\
			/*****************************************/\
			/* now do three more radix-3 transforms: */\
			/*****************************************/\
			/* Block 1: */\
			"movl	%[__o0]	,%%esi 	\n\t"/* __o0-8: esi,edi store output addresses throughout */\
			"movl	%[__i0]	,%%eax\n\t"\
			"movl	%[__i3]	,%%ebx\n\t"\
			"movl	%[__i6]	,%%ecx\n\t"\
				"movaps	    (%%ebx)	,%%xmm4\n\t"\
				"movaps	0x10(%%ebx)	,%%xmm5\n\t"\
				"movaps	    (%%ecx)	,%%xmm2\n\t"\
				"movaps	0x10(%%ecx)	,%%xmm3\n\t"\
				"movaps	    (%%eax)	,%%xmm0\n\t"\
				"movaps	0x10(%%eax)	,%%xmm1\n\t"\
				"subpd	%%xmm2		,%%xmm4\n\t"\
				"subpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm2\n\t"\
				"addpd	%%xmm3		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm2\n\t"\
				"addpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm2		,%%xmm0\n\t"\
				"addpd	%%xmm3		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%esi)\n\t"\
				"movaps	%%xmm1		,0x10(%%esi)\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm4\n\t"\
				"mulpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm0		,%%xmm2\n\t"\
				"addpd	%%xmm1		,%%xmm3\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
			"movl	%[__o3]	,%%edi 	\n\t"\
			"movl	%[__o6]	,%%esi 	\n\t"\
				"addpd	%%xmm5		,%%xmm2\n\t"\
				"subpd	%%xmm4		,%%xmm3\n\t"\
				"subpd	%%xmm5		,%%xmm0\n\t"\
				"addpd	%%xmm4		,%%xmm1\n\t"\
				"movaps	%%xmm2		,    (%%edi)\n\t"\
				"movaps	%%xmm3		,0x10(%%edi)\n\t"\
				"movaps	%%xmm0		,    (%%esi)\n\t"\
				"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			/* Block 2: */\
			"movl	%[__i1]	,%%eax\n\t"\
			"movl	%[__i4]	,%%ebx\n\t"\
			"movl	%[__i7]	,%%ecx\n\t"\
				"movaps	    (%%ebx)	,%%xmm2\n\t"\
				"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
				"movaps	    (%%edx)	,%%xmm6\n\t"/* c1 */\
				"movaps	0x10(%%edx)	,%%xmm7\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm0\n\t"\
				"mulpd	%%xmm7		,%%xmm1\n\t"\
				"addpd	%%xmm1		,%%xmm2\n\t"\
				"subpd	%%xmm0		,%%xmm3\n\t"\
				"movaps	    (%%ecx)	,%%xmm4\n\t"\
				"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
				"movaps	0x20(%%edx)	,%%xmm6\n\t"/* c2 */\
				"movaps	0x30(%%edx)	,%%xmm7\n\t"\
				"movaps	%%xmm4		,%%xmm0\n\t"\
				"movaps	%%xmm5		,%%xmm1\n\t"\
				"mulpd	%%xmm6		,%%xmm4\n\t"\
				"mulpd	%%xmm6		,%%xmm5\n\t"\
				"mulpd	%%xmm7		,%%xmm0\n\t"\
				"mulpd	%%xmm7		,%%xmm1\n\t"\
				"addpd	%%xmm1		,%%xmm4\n\t"\
				"subpd	%%xmm0		,%%xmm5\n\t"\
				"movaps	0x40(%%edx)	,%%xmm6\n\t"/* c3m1 */\
				"movaps	0x50(%%edx)	,%%xmm7\n\t"\
				"movaps	    (%%eax)	,%%xmm0\n\t"\
				"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movl	%[__o7]	,%%esi 	\n\t"\
				"subpd	%%xmm4		,%%xmm2\n\t"\
				"subpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm4\n\t"\
				"addpd	%%xmm5		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm4\n\t"\
				"addpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm4		,%%xmm0\n\t"\
				"addpd	%%xmm5		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%esi)\n\t"\
				"movaps	%%xmm1		,0x10(%%esi)\n\t"\
				"mulpd	%%xmm6		,%%xmm4\n\t"\
				"mulpd	%%xmm6		,%%xmm5\n\t"\
				"mulpd	%%xmm7		,%%xmm2\n\t"\
				"mulpd	%%xmm7		,%%xmm3\n\t"\
				"addpd	%%xmm0		,%%xmm4\n\t"\
				"addpd	%%xmm1		,%%xmm5\n\t"\
			"movl	%[__o1]	,%%edi 	\n\t"\
			"movl	%[__o4]	,%%esi 	\n\t"\
				"movaps	%%xmm4		,%%xmm0\n\t"\
				"movaps	%%xmm5		,%%xmm1\n\t"\
				"addpd	%%xmm3		,%%xmm4\n\t"\
				"subpd	%%xmm2		,%%xmm5\n\t"\
				"subpd	%%xmm3		,%%xmm0\n\t"\
				"addpd	%%xmm2		,%%xmm1\n\t"\
				"movaps	%%xmm4		,    (%%edi)\n\t"\
				"movaps	%%xmm5		,0x10(%%edi)\n\t"\
				"movaps	%%xmm0		,    (%%esi)\n\t"\
				"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			/* Block 3: */\
			"movl	%[__i2]	,%%eax\n\t"\
			"movl	%[__i5]	,%%ebx\n\t"\
			"movl	%[__i8]	,%%ecx\n\t"\
				"movaps	    (%%ebx)	,%%xmm2\n\t"\
				"movaps	0x10(%%ebx)	,%%xmm3\n\t"\
				"movaps	0x20(%%edx)	,%%xmm6\n\t"/* c2 */\
				"movaps	0x30(%%edx)	,%%xmm7\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm0\n\t"\
				"mulpd	%%xmm7		,%%xmm1\n\t"\
				"addpd	%%xmm1		,%%xmm2\n\t"\
				"subpd	%%xmm0		,%%xmm3\n\t"\
				"movaps	    (%%ecx)	,%%xmm4\n\t"\
				"movaps	0x10(%%ecx)	,%%xmm5\n\t"\
				"movaps	0x60(%%edx)	,%%xmm6\n\t"/* c4 */\
				"movaps	0x70(%%edx)	,%%xmm7\n\t"\
				"movaps	%%xmm4		,%%xmm0\n\t"\
				"movaps	%%xmm5		,%%xmm1\n\t"\
				"mulpd	%%xmm6		,%%xmm4\n\t"\
				"mulpd	%%xmm6		,%%xmm5\n\t"\
				"mulpd	%%xmm7		,%%xmm0\n\t"\
				"mulpd	%%xmm7		,%%xmm1\n\t"\
				"addpd	%%xmm1		,%%xmm4\n\t"\
				"subpd	%%xmm0		,%%xmm5\n\t"\
				"movaps	0x40(%%edx)	,%%xmm6\n\t"\
				"movaps	0x50(%%edx)	,%%xmm7\n\t"\
				"movaps	    (%%eax)	,%%xmm0\n\t"\
				"movaps	0x10(%%eax)	,%%xmm1\n\t"\
			"movl	%[__o5]	,%%esi 	\n\t"\
				"subpd	%%xmm4		,%%xmm2\n\t"\
				"subpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm4\n\t"\
				"addpd	%%xmm5		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm4\n\t"\
				"addpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm4		,%%xmm0\n\t"\
				"addpd	%%xmm5		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%esi)\n\t"\
				"movaps	%%xmm1		,0x10(%%esi)\n\t"\
				"mulpd	%%xmm6		,%%xmm4\n\t"\
				"mulpd	%%xmm6		,%%xmm5\n\t"\
				"mulpd	%%xmm7		,%%xmm2\n\t"\
				"mulpd	%%xmm7		,%%xmm3\n\t"\
				"addpd	%%xmm0		,%%xmm4\n\t"\
				"addpd	%%xmm1		,%%xmm5\n\t"\
			"movl	%[__o8]	,%%edi 	\n\t"\
			"movl	%[__o2]	,%%esi 	\n\t"\
				"movaps	%%xmm4		,%%xmm0\n\t"\
				"movaps	%%xmm5		,%%xmm1\n\t"\
				"addpd	%%xmm3		,%%xmm4\n\t"\
				"subpd	%%xmm2		,%%xmm5\n\t"\
				"subpd	%%xmm3		,%%xmm0\n\t"\
				"addpd	%%xmm2		,%%xmm1\n\t"\
				"movaps	%%xmm4		,    (%%edi)\n\t"\
				"movaps	%%xmm5		,0x10(%%edi)\n\t"\
				"movaps	%%xmm0		,    (%%esi)\n\t"\
				"movaps	%%xmm1		,0x10(%%esi)\n\t"\
			"popl %%ebx	\n\t"\
				:					/* outputs: none */\
				: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
				 ,[__i1] "m" (Xi1)\
				 ,[__i2] "m" (Xi2)\
				 ,[__i3] "m" (Xi3)\
				 ,[__i4] "m" (Xi4)\
				 ,[__i5] "m" (Xi5)\
				 ,[__i6] "m" (Xi6)\
				 ,[__i7] "m" (Xi7)\
				 ,[__i8] "m" (Xi8)\
				 ,[__c1] "m" (Xc1)\
				 ,[__o0] "m" (Xo0)\
				 ,[__o1] "m" (Xo1)\
				 ,[__o2] "m" (Xo2)\
				 ,[__o3] "m" (Xo3)\
				 ,[__o4] "m" (Xo4)\
				 ,[__o5] "m" (Xo5)\
				 ,[__o6] "m" (Xo6)\
				 ,[__o7] "m" (Xo7)\
				 ,[__o8] "m" (Xo8)\
				: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
				);\
			}

	  #else

		#ifdef USE_AVX

			/*...Radix-9 DIF: Ins in memory locations __i0-8.\
				Outs at memory locations __o0-8, assumed disjoint with inputs:\
			*/\
			#define SSE2_RADIX_09_DIF(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8)\
			{\
			__asm__ volatile (\
			"movq	%[__i0]	,%%rax 	\n\t"/* __i0-8; e[abc]x store input addresses */\
			"movq	%[__o0]	,%%rsi	\n\t"\
			"movq	%[__c1]	,%%rdx 	\n\t"/* edx stores trig addresses throughout */\
			/* Block 1: */\
				"movq	%%rax		,%%rbx\n\t"\
				"movq	%%rax		,%%rcx\n\t"\
				"movq	%[__i3]	,%%rbx 	\n\t"\
				"movq	%[__i6]	,%%rcx 	\n\t"\
				"vmovaps	    (%%rbx)	,%%ymm2\n\t"\
				"vmovaps	0x20(%%rbx)	,%%ymm3\n\t"\
				"vmovaps	    (%%rax)	,%%ymm0\n\t"\
				"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
				"vmovaps	    (%%rcx)	,%%ymm6\n\t"\
				"vmovaps	0x20(%%rcx)	,%%ymm7\n\t"\
				"vmovaps	%%ymm2		,%%ymm4\n\t"\
				"vmovaps	%%ymm3		,%%ymm5\n\t"\
				"vaddpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm7,%%ymm3,%%ymm3\n\t"\
				"vsubpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
				"vsubpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t"\
				"vmovaps	0x80(%%rdx)	,%%ymm6\n\t"/* c3m1 */\
				"vmovaps	0xa0(%%rdx)	,%%ymm7\n\t"\
				"vmovaps	%%ymm0		,    (%%rsi)	/* <- t00 */\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t01 */\n\t"\
			"movq	%[__o1]	,%%rdi	\n\t"\
			"movq	%[__o2]	,%%rsi	\n\t"\
				"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
				"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
				"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t"\
				"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t"\
				"vmovaps	%%ymm2		,%%ymm0\n\t"\
				"vmovaps	%%ymm3		,%%ymm1\n\t"\
				"vsubpd	%%ymm5,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm4,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm5,%%ymm0,%%ymm0\n\t"\
				"vsubpd	%%ymm4,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm2		,    (%%rdi)	/* <- t02 */\n\t"\
				"vmovaps	%%ymm3		,0x20(%%rdi)	/* <- t03 */\n\t"\
				"vmovaps	%%ymm0		,    (%%rsi)	/* <- t04 */\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t05 */\n\t"\
			/* Block 2: */\
				"movq	%[__i1]	,%%rax 	\n\t"\
				"movq	%[__i4]	,%%rbx 	\n\t"\
				"movq	%[__i7]	,%%rcx 	\n\t"\
				"vmovaps	    (%%rbx)	,%%ymm4\n\t"\
				"vmovaps	0x20(%%rbx)	,%%ymm5\n\t"\
				"vmovaps	    (%%rax)	,%%ymm0\n\t"\
				"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
				"vmovaps	    (%%rcx)	,%%ymm2\n\t"\
				"vmovaps	0x20(%%rcx)	,%%ymm3\n\t"\
			"movq	%[__o3]	,%%rdi	\n\t"\
				"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
				"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm0		,    (%%rdi)	/* <- t06 */\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rdi)	/* <- t07 */\n\t"\
			"movq	%[__o4]	,%%rsi	\n\t"\
			"movq	%[__o5]	,%%rdi	\n\t"\
				"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
				"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
				"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t"\
				"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t"\
				"vmovaps	%%ymm2		,%%ymm0\n\t"\
				"vmovaps	%%ymm3		,%%ymm1\n\t"\
				"vsubpd	%%ymm5,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm4,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm5,%%ymm0,%%ymm0\n\t"\
				"vsubpd	%%ymm4,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm2		,    (%%rsi)	/* <- t08 */\n\t"\
				"vmovaps	%%ymm3		,0x20(%%rsi)	/* <- t09 */\n\t"\
				"vmovaps	%%ymm0		,    (%%rdi)	/* <- t0a */\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rdi)	/* <- t0b */\n\t"\
			/* Block 3: */\
				"movq	%[__i2]	,%%rax 	\n\t"\
				"movq	%[__i5]	,%%rbx 	\n\t"\
				"movq	%[__i8]	,%%rcx 	\n\t"\
				"vmovaps	    (%%rbx)	,%%ymm4\n\t"\
				"vmovaps	0x20(%%rbx)	,%%ymm5\n\t"\
				"vmovaps	    (%%rax)	,%%ymm0\n\t"\
				"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
				"vmovaps	    (%%rcx)	,%%ymm2\n\t"\
				"vmovaps	0x20(%%rcx)	,%%ymm3\n\t"\
			"movq	%[__o6]	,%%rsi	\n\t"\
				"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
				"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm0		,    (%%rsi)	/* <- t0c */\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t0d */\n\t"\
			"movq	%[__o7]	,%%rdi	\n\t"\
			"movq	%[__o8]	,%%rsi	\n\t"\
				"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
				"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
				"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t"\
				"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t"\
				"vmovaps	%%ymm2		,%%ymm0\n\t"\
				"vmovaps	%%ymm3		,%%ymm1\n\t"\
				"vsubpd	%%ymm5,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm4,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm5,%%ymm0,%%ymm0\n\t"\
				"vsubpd	%%ymm4,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm2		,    (%%rdi)	/* <- t0e */\n\t"\
				"vmovaps	%%ymm3		,0x20(%%rdi)	/* <- t0f */\n\t"\
				"vmovaps	%%ymm0		,    (%%rsi)	/* <- t0g */\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t0h */\n\t"\
			/*****************************************/\
			/* now do three more radix-3 transforms: */\
			/*****************************************/\
			/* Block 1: */\
			"movq	%[__o0]	,%%rax	\n\t"\
			"movq	%[__o3]	,%%rbx	\n\t"\
			"movq	%[__o6]	,%%rcx	\n\t"\
				"vmovaps	    (%%rbx)	,%%ymm4\n\t"\
				"vmovaps	0x20(%%rbx)	,%%ymm5\n\t"\
				"vmovaps	    (%%rcx)	,%%ymm2\n\t"\
				"vmovaps	0x20(%%rcx)	,%%ymm3\n\t"\
				"vmovaps	    (%%rax)	,%%ymm0\n\t"\
				"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
				"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
				"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm0		,    (%%rax)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rax)\n\t"\
				"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
				"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
				"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t"\
				"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t"\
				"vmovaps	%%ymm2		,%%ymm0\n\t"\
				"vmovaps	%%ymm3		,%%ymm1\n\t"\
				"vsubpd	%%ymm5,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm4,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm5,%%ymm0,%%ymm0\n\t"\
				"vsubpd	%%ymm4,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm2		,    (%%rbx)\n\t"\
				"vmovaps	%%ymm3		,0x20(%%rbx)\n\t"\
				"vmovaps	%%ymm0		,    (%%rcx)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rcx)\n\t"\
			/* Block 2: */\
			"movq	%[__o1]	,%%rax	\n\t"\
			"movq	%[__o4]	,%%rbx	\n\t"\
			"movq	%[__o7]	,%%rcx	\n\t"\
				"vmovaps	    (%%rbx)	,%%ymm2\n\t"\
				"vmovaps	0x20(%%rbx)	,%%ymm3\n\t"\
				"vmovaps	    (%%rdx)	,%%ymm6\n\t"/* c1 */\
				"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
				"vmovaps	%%ymm2		,%%ymm0\n\t"\
				"vmovaps	%%ymm3		,%%ymm1\n\t"\
				"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
				"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
				"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
				"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t"\
				"vsubpd	%%ymm1,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm0,%%ymm3,%%ymm3\n\t"\
				"vmovaps	    (%%rcx)	,%%ymm4\n\t"\
				"vmovaps	0x20(%%rcx)	,%%ymm5\n\t"\
				"vmovaps	0x40(%%rdx)	,%%ymm6\n\t"/* c2 */\
				"vmovaps	0x60(%%rdx)	,%%ymm7\n\t"\
				"vmovaps	%%ymm4		,%%ymm0\n\t"\
				"vmovaps	%%ymm5		,%%ymm1\n\t"\
				"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
				"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t"\
				"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
				"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t"\
				"vsubpd	%%ymm1,%%ymm4,%%ymm4\n\t"\
				"vaddpd	%%ymm0,%%ymm5,%%ymm5\n\t"\
				"vmovaps	0x80(%%rdx)	,%%ymm6\n\t"/* c3m1 */\
				"vmovaps	0xa0(%%rdx)	,%%ymm7\n\t"\
				"vmovaps	    (%%rax)	,%%ymm0\n\t"\
				"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
				"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
				"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t"\
				"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
				"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm4,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm5,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm0		,    (%%rax)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rax)\n\t"\
				"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
				"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t"\
				"vmulpd	%%ymm7,%%ymm2,%%ymm2\n\t"\
				"vmulpd	%%ymm7,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm0,%%ymm4,%%ymm4\n\t"\
				"vaddpd	%%ymm1,%%ymm5,%%ymm5\n\t"\
				"vmovaps	%%ymm4		,%%ymm0\n\t"\
				"vmovaps	%%ymm5		,%%ymm1\n\t"\
				"vsubpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
				"vaddpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm3,%%ymm0,%%ymm0\n\t"\
				"vsubpd	%%ymm2,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm4		,    (%%rbx)\n\t"\
				"vmovaps	%%ymm5		,0x20(%%rbx)\n\t"\
				"vmovaps	%%ymm0		,    (%%rcx)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rcx)\n\t"\
			/* Block 3: */\
			"movq	%[__o2]	,%%rax	\n\t"\
			"movq	%[__o5]	,%%rbx	\n\t"\
			"movq	%[__o8]	,%%rcx	\n\t"\
				"vmovaps	    (%%rbx)	,%%ymm2\n\t"\
				"vmovaps	0x20(%%rbx)	,%%ymm3\n\t"\
				"vmovaps	0x40(%%rdx)	,%%ymm6\n\t"/* c2 */\
				"vmovaps	0x60(%%rdx)	,%%ymm7\n\t"\
				"vmovaps	%%ymm2		,%%ymm0\n\t"\
				"vmovaps	%%ymm3		,%%ymm1\n\t"\
				"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
				"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
				"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
				"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t"\
				"vsubpd	%%ymm1,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm0,%%ymm3,%%ymm3\n\t"\
				"vmovaps	    (%%rcx)	,%%ymm4\n\t"\
				"vmovaps	0x20(%%rcx)	,%%ymm5\n\t"\
				"vmovaps	0xc0(%%rdx)	,%%ymm6\n\t"/* c4 */\
				"vmovaps	0xe0(%%rdx)	,%%ymm7\n\t"\
				"vmovaps	%%ymm4		,%%ymm0\n\t"\
				"vmovaps	%%ymm5		,%%ymm1\n\t"\
				"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
				"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t"\
				"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
				"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t"\
				"vsubpd	%%ymm1,%%ymm4,%%ymm4\n\t"\
				"vaddpd	%%ymm0,%%ymm5,%%ymm5\n\t"\
				"vmovaps	0x80(%%rdx)	,%%ymm6\n\t"/* c3m1 */\
				"vmovaps	0xa0(%%rdx)	,%%ymm7\n\t"\
				"vmovaps	    (%%rax)	,%%ymm0\n\t"\
				"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
				"\n\t"\
				"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
				"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t"\
				"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
				"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm4,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm5,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm0		,    (%%rax)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rax)\n\t"\
				"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
				"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t"\
				"vmulpd	%%ymm7,%%ymm2,%%ymm2\n\t"\
				"vmulpd	%%ymm7,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm0,%%ymm4,%%ymm4\n\t"\
				"vaddpd	%%ymm1,%%ymm5,%%ymm5\n\t"\
				"\n\t"\
				"vmovaps	%%ymm4		,%%ymm0\n\t"\
				"vmovaps	%%ymm5		,%%ymm1\n\t"\
				"vsubpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
				"vaddpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm3,%%ymm0,%%ymm0\n\t"\
				"vsubpd	%%ymm2,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm4		,    (%%rbx)\n\t"\
				"vmovaps	%%ymm5		,0x20(%%rbx)\n\t"\
				"vmovaps	%%ymm0		,    (%%rcx)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rcx)\n\t"\
				:					/* outputs: none */\
				: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
				 ,[__i1] "m" (Xi1)\
				 ,[__i2] "m" (Xi2)\
				 ,[__i3] "m" (Xi3)\
				 ,[__i4] "m" (Xi4)\
				 ,[__i5] "m" (Xi5)\
				 ,[__i6] "m" (Xi6)\
				 ,[__i7] "m" (Xi7)\
				 ,[__i8] "m" (Xi8)\
				 ,[__c1] "m" (Xc1)\
				 ,[__o0] "m" (Xo0)\
				 ,[__o1] "m" (Xo1)\
				 ,[__o2] "m" (Xo2)\
				 ,[__o3] "m" (Xo3)\
				 ,[__o4] "m" (Xo4)\
				 ,[__o5] "m" (Xo5)\
				 ,[__o6] "m" (Xo6)\
				 ,[__o7] "m" (Xo7)\
				 ,[__o8] "m" (Xo8)\
				: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
				);\
			}
			
			/*...Radix-9 DIT: Ins in memory locations __i0-8.\
				Outs at memory locations __o0-8, assumed disjoint with inputs:\
			*/\
			#define SSE2_RADIX_09_DIT(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8)\
			{\
			__asm__ volatile (\
			"movq	%[__i0]	,%%rax\n\t"\
			"movq	%[__i1]	,%%rbx\n\t"\
			"movq	%[__i2]	,%%rcx\n\t"\
			"movq	%[__c1]	,%%rdx 	\n\t"/* edx stores trig addresses throughout */\
			/* Block 1: */\
				"vmovaps	    (%%rbx)	,%%ymm2\n\t"\
				"vmovaps	0x20(%%rbx)	,%%ymm3\n\t"\
				"vmovaps	    (%%rax)	,%%ymm0\n\t"\
				"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
				"vmovaps	    (%%rcx)	,%%ymm6\n\t"\
				"vmovaps	0x20(%%rcx)	,%%ymm7\n\t"\
				"vmovaps	%%ymm2		,%%ymm4\n\t"\
				"vmovaps	%%ymm3		,%%ymm5\n\t"\
				"vaddpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm7,%%ymm3,%%ymm3\n\t"\
				"vsubpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
				"vsubpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t"\
				"vmovaps	0x80(%%rdx)	,%%ymm6\n\t"/* c3m1 */\
				"vmovaps	0xa0(%%rdx)	,%%ymm7\n\t"\
				"vmovaps	%%ymm0		,    (%%rax)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rax)\n\t"\
				"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
				"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
				"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t"\
				"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t"\
				"vmovaps	%%ymm2		,%%ymm0\n\t"\
				"vmovaps	%%ymm3		,%%ymm1\n\t"\
				"vaddpd	%%ymm5,%%ymm2,%%ymm2\n\t"\
				"vsubpd	%%ymm4,%%ymm3,%%ymm3\n\t"\
				"vsubpd	%%ymm5,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm4,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm2		,    (%%rbx)\n\t"\
				"vmovaps	%%ymm3		,0x20(%%rbx)\n\t"\
				"vmovaps	%%ymm0		,    (%%rcx)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rcx)\n\t"\
			/* Block 2: */\
			"movq	%[__i3]	,%%rax\n\t"\
			"movq	%[__i4]	,%%rbx\n\t"\
			"movq	%[__i5]	,%%rcx\n\t"\
				"vmovaps	    (%%rbx)	,%%ymm4\n\t"\
				"vmovaps	0x20(%%rbx)	,%%ymm5\n\t"\
				"vmovaps	    (%%rcx)	,%%ymm2\n\t"\
				"vmovaps	0x20(%%rcx)	,%%ymm3\n\t"\
				"vmovaps	    (%%rax)	,%%ymm0\n\t"\
				"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
				"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
				"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm0		,    (%%rax)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rax)\n\t"\
				"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
				"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
				"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t"\
				"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t"\
				"vmovaps	%%ymm2		,%%ymm0\n\t"\
				"vmovaps	%%ymm3		,%%ymm1\n\t"\
				"vaddpd	%%ymm5,%%ymm2,%%ymm2\n\t"\
				"vsubpd	%%ymm4,%%ymm3,%%ymm3\n\t"\
				"vsubpd	%%ymm5,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm4,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm2		,    (%%rbx)\n\t"\
				"vmovaps	%%ymm3		,0x20(%%rbx)\n\t"\
				"vmovaps	%%ymm0		,    (%%rcx)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rcx)\n\t"\
			/* Block 3: */\
			"movq	%[__i6]	,%%rax\n\t"\
			"movq	%[__i7]	,%%rbx\n\t"\
			"movq	%[__i8]	,%%rcx\n\t"\
				"vmovaps	    (%%rbx)	,%%ymm4\n\t"\
				"vmovaps	0x20(%%rbx)	,%%ymm5\n\t"\
				"vmovaps	    (%%rcx)	,%%ymm2\n\t"\
				"vmovaps	0x20(%%rcx)	,%%ymm3\n\t"\
				"vmovaps	    (%%rax)	,%%ymm0\n\t"\
				"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
				"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
				"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm0		,    (%%rax)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rax)\n\t"\
				"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
				"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
				"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t"\
				"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t"\
				"vmovaps	%%ymm2		,%%ymm0\n\t"\
				"vmovaps	%%ymm3		,%%ymm1\n\t"\
				"vaddpd	%%ymm5,%%ymm2,%%ymm2\n\t"\
				"vsubpd	%%ymm4,%%ymm3,%%ymm3\n\t"\
				"vsubpd	%%ymm5,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm4,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm2		,    (%%rbx)\n\t"\
				"vmovaps	%%ymm3		,0x20(%%rbx)\n\t"\
				"vmovaps	%%ymm0		,    (%%rcx)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rcx)\n\t"\
			/*****************************************/\
			/* now do three more radix-3 transforms: */\
			/*****************************************/\
			/* Block 1: */\
			"movq	%[__o0]	,%%rsi 	\n\t"/* __o0-8: esi,edi store output addresses throughout */\
			"movq	%[__i0]	,%%rax\n\t"\
			"movq	%[__i3]	,%%rbx\n\t"\
			"movq	%[__i6]	,%%rcx\n\t"\
				"vmovaps	    (%%rbx)	,%%ymm4\n\t"\
				"vmovaps	0x20(%%rbx)	,%%ymm5\n\t"\
				"vmovaps	    (%%rcx)	,%%ymm2\n\t"\
				"vmovaps	0x20(%%rcx)	,%%ymm3\n\t"\
				"vmovaps	    (%%rax)	,%%ymm0\n\t"\
				"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
				"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
				"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm0		,    (%%rsi)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rsi)\n\t"\
				"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
				"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
				"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t"\
				"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t"\
				"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t"\
				"vmovaps	%%ymm2		,%%ymm0\n\t"\
				"vmovaps	%%ymm3		,%%ymm1\n\t"\
			"movq	%[__o3]	,%%rdi 	\n\t"\
			"movq	%[__o6]	,%%rsi 	\n\t"\
				"vaddpd	%%ymm5,%%ymm2,%%ymm2\n\t"\
				"vsubpd	%%ymm4,%%ymm3,%%ymm3\n\t"\
				"vsubpd	%%ymm5,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm4,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm2		,    (%%rdi)\n\t"\
				"vmovaps	%%ymm3		,0x20(%%rdi)\n\t"\
				"vmovaps	%%ymm0		,    (%%rsi)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rsi)\n\t"\
			/* Block 2: */\
			"movq	%[__i1]	,%%rax\n\t"\
			"movq	%[__i4]	,%%rbx\n\t"\
			"movq	%[__i7]	,%%rcx\n\t"\
				"vmovaps	    (%%rbx)	,%%ymm2\n\t"\
				"vmovaps	0x20(%%rbx)	,%%ymm3\n\t"\
				"vmovaps	    (%%rdx)	,%%ymm6\n\t"/* c1 */\
				"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
				"vmovaps	%%ymm2		,%%ymm0\n\t"\
				"vmovaps	%%ymm3		,%%ymm1\n\t"\
				"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
				"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
				"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
				"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t"\
				"vaddpd	%%ymm1,%%ymm2,%%ymm2\n\t"\
				"vsubpd	%%ymm0,%%ymm3,%%ymm3\n\t"\
				"vmovaps	    (%%rcx)	,%%ymm4\n\t"\
				"vmovaps	0x20(%%rcx)	,%%ymm5\n\t"\
				"vmovaps	0x40(%%rdx)	,%%ymm6\n\t"/* c2 */\
				"vmovaps	0x60(%%rdx)	,%%ymm7\n\t"\
				"vmovaps	%%ymm4		,%%ymm0\n\t"\
				"vmovaps	%%ymm5		,%%ymm1\n\t"\
				"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
				"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t"\
				"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
				"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t"\
				"vaddpd	%%ymm1,%%ymm4,%%ymm4\n\t"\
				"vsubpd	%%ymm0,%%ymm5,%%ymm5\n\t"\
				"vmovaps	0x80(%%rdx)	,%%ymm6\n\t"/* c3m1 */\
				"vmovaps	0xa0(%%rdx)	,%%ymm7\n\t"\
				"vmovaps	    (%%rax)	,%%ymm0\n\t"\
				"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
			"movq	%[__o7]	,%%rsi 	\n\t"\
				"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
				"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t"\
				"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
				"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm4,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm5,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm0		,    (%%rsi)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rsi)\n\t"\
				"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
				"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t"\
				"vmulpd	%%ymm7,%%ymm2,%%ymm2\n\t"\
				"vmulpd	%%ymm7,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm0,%%ymm4,%%ymm4\n\t"\
				"vaddpd	%%ymm1,%%ymm5,%%ymm5\n\t"\
			"movq	%[__o1]	,%%rdi 	\n\t"\
			"movq	%[__o4]	,%%rsi 	\n\t"\
				"vmovaps	%%ymm4		,%%ymm0\n\t"\
				"vmovaps	%%ymm5		,%%ymm1\n\t"\
				"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
				"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
				"vsubpd	%%ymm3,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm2,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm4		,    (%%rdi)\n\t"\
				"vmovaps	%%ymm5		,0x20(%%rdi)\n\t"\
				"vmovaps	%%ymm0		,    (%%rsi)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rsi)\n\t"\
			/* Block 3: */\
			"movq	%[__i2]	,%%rax\n\t"\
			"movq	%[__i5]	,%%rbx\n\t"\
			"movq	%[__i8]	,%%rcx\n\t"\
				"vmovaps	    (%%rbx)	,%%ymm2\n\t"\
				"vmovaps	0x20(%%rbx)	,%%ymm3\n\t"\
				"vmovaps	0x40(%%rdx)	,%%ymm6\n\t"/* c2 */\
				"vmovaps	0x60(%%rdx)	,%%ymm7\n\t"\
				"vmovaps	%%ymm2		,%%ymm0\n\t"\
				"vmovaps	%%ymm3		,%%ymm1\n\t"\
				"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t"\
				"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t"\
				"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
				"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t"\
				"vaddpd	%%ymm1,%%ymm2,%%ymm2\n\t"\
				"vsubpd	%%ymm0,%%ymm3,%%ymm3\n\t"\
				"vmovaps	    (%%rcx)	,%%ymm4\n\t"\
				"vmovaps	0x20(%%rcx)	,%%ymm5\n\t"\
				"vmovaps	0xc0(%%rdx)	,%%ymm6\n\t"/* c4 */\
				"vmovaps	0xe0(%%rdx)	,%%ymm7\n\t"\
				"vmovaps	%%ymm4		,%%ymm0\n\t"\
				"vmovaps	%%ymm5		,%%ymm1\n\t"\
				"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
				"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t"\
				"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t"\
				"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t"\
				"vaddpd	%%ymm1,%%ymm4,%%ymm4\n\t"\
				"vsubpd	%%ymm0,%%ymm5,%%ymm5\n\t"\
				"vmovaps	0x80(%%rdx)	,%%ymm6\n\t"\
				"vmovaps	0xa0(%%rdx)	,%%ymm7\n\t"\
				"vmovaps	    (%%rax)	,%%ymm0\n\t"\
				"vmovaps	0x20(%%rax)	,%%ymm1\n\t"\
			"movq	%[__o5]	,%%rsi 	\n\t"\
				"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t"\
				"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t"\
				"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t"\
				"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t"\
				"vaddpd	%%ymm4,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm5,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm0		,    (%%rsi)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rsi)\n\t"\
				"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t"\
				"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t"\
				"vmulpd	%%ymm7,%%ymm2,%%ymm2\n\t"\
				"vmulpd	%%ymm7,%%ymm3,%%ymm3\n\t"\
				"vaddpd	%%ymm0,%%ymm4,%%ymm4\n\t"\
				"vaddpd	%%ymm1,%%ymm5,%%ymm5\n\t"\
			"movq	%[__o8]	,%%rdi 	\n\t"\
			"movq	%[__o2]	,%%rsi 	\n\t"\
				"vmovaps	%%ymm4		,%%ymm0\n\t"\
				"vmovaps	%%ymm5		,%%ymm1\n\t"\
				"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t"\
				"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t"\
				"vsubpd	%%ymm3,%%ymm0,%%ymm0\n\t"\
				"vaddpd	%%ymm2,%%ymm1,%%ymm1\n\t"\
				"vmovaps	%%ymm4		,    (%%rdi)\n\t"\
				"vmovaps	%%ymm5		,0x20(%%rdi)\n\t"\
				"vmovaps	%%ymm0		,    (%%rsi)\n\t"\
				"vmovaps	%%ymm1		,0x20(%%rsi)\n\t"\
				:					/* outputs: none */\
				: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
				 ,[__i1] "m" (Xi1)\
				 ,[__i2] "m" (Xi2)\
				 ,[__i3] "m" (Xi3)\
				 ,[__i4] "m" (Xi4)\
				 ,[__i5] "m" (Xi5)\
				 ,[__i6] "m" (Xi6)\
				 ,[__i7] "m" (Xi7)\
				 ,[__i8] "m" (Xi8)\
				 ,[__c1] "m" (Xc1)\
				 ,[__o0] "m" (Xo0)\
				 ,[__o1] "m" (Xo1)\
				 ,[__o2] "m" (Xo2)\
				 ,[__o3] "m" (Xo3)\
				 ,[__o4] "m" (Xo4)\
				 ,[__o5] "m" (Xo5)\
				 ,[__o6] "m" (Xo6)\
				 ,[__o7] "m" (Xo7)\
				 ,[__o8] "m" (Xo8)\
				: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
				);\
			}

		#else	// 64-bit SSE2:

			// Simple 64-bit-ified version of the same-named 32-bit ASM macros, using just xmm0-7:

			/*...Radix-9 DIF: Ins in memory locations __i0-8.\
				Outs at memory locations __o0-8, assumed disjoint with inputs:\
			*/\
			#define SSE2_RADIX_09_DIF(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8)\
			{\
			__asm__ volatile (\
			"movq	%[__i0]	,%%rax 	\n\t"/* __i0-8; e[abc]x store input addresses */\
			"movq	%[__o0]	,%%rsi	\n\t"\
			"movq	%[__c1]	,%%rdx 	\n\t"/* edx stores trig addresses throughout */\
			/* Block 1: */\
				"movq	%%rax		,%%rbx\n\t"\
				"movq	%%rax		,%%rcx\n\t"\
				"movq	%[__i3]	,%%rbx 	\n\t"\
				"movq	%[__i6]	,%%rcx 	\n\t"\
				"movaps	    (%%rbx)	,%%xmm2\n\t"\
				"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
				"movaps	    (%%rax)	,%%xmm0\n\t"\
				"movaps	0x10(%%rax)	,%%xmm1\n\t"\
				"movaps	    (%%rcx)	,%%xmm6\n\t"\
				"movaps	0x10(%%rcx)	,%%xmm7\n\t"\
				"movaps	%%xmm2		,%%xmm4\n\t"\
				"movaps	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm6		,%%xmm2\n\t"\
				"addpd	%%xmm7		,%%xmm3\n\t"\
				"subpd	%%xmm6		,%%xmm4\n\t"\
				"subpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm0\n\t"\
				"addpd	%%xmm3		,%%xmm1\n\t"\
				"movaps	0x40(%%rdx)	,%%xmm6\n\t"/* c3m1 */\
				"movaps	0x50(%%rdx)	,%%xmm7\n\t"\
				"movaps	%%xmm0		,    (%%rsi)	/* <- t00 */\n\t"\
				"movaps	%%xmm1		,0x10(%%rsi)	/* <- t01 */\n\t"\
			"movq	%[__o1]	,%%rdi	\n\t"\
			"movq	%[__o2]	,%%rsi	\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm4\n\t"\
				"mulpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm0		,%%xmm2\n\t"\
				"addpd	%%xmm1		,%%xmm3\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"subpd	%%xmm5		,%%xmm2\n\t"\
				"addpd	%%xmm4		,%%xmm3\n\t"\
				"addpd	%%xmm5		,%%xmm0\n\t"\
				"subpd	%%xmm4		,%%xmm1\n\t"\
				"movaps	%%xmm2		,    (%%rdi)	/* <- t02 */\n\t"\
				"movaps	%%xmm3		,0x10(%%rdi)	/* <- t03 */\n\t"\
				"movaps	%%xmm0		,    (%%rsi)	/* <- t04 */\n\t"\
				"movaps	%%xmm1		,0x10(%%rsi)	/* <- t05 */\n\t"\
			/* Block 2: */\
				"movq	%[__i1]	,%%rax 	\n\t"\
				"movq	%[__i4]	,%%rbx 	\n\t"\
				"movq	%[__i7]	,%%rcx 	\n\t"\
				"movaps	    (%%rbx)	,%%xmm4\n\t"\
				"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
				"movaps	    (%%rax)	,%%xmm0\n\t"\
				"movaps	0x10(%%rax)	,%%xmm1\n\t"\
				"movaps	    (%%rcx)	,%%xmm2\n\t"\
				"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movq	%[__o3]	,%%rdi	\n\t"\
				"subpd	%%xmm2		,%%xmm4\n\t"\
				"subpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm2\n\t"\
				"addpd	%%xmm3		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm2\n\t"\
				"addpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm2		,%%xmm0\n\t"\
				"addpd	%%xmm3		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%rdi)	/* <- t06 */\n\t"\
				"movaps	%%xmm1		,0x10(%%rdi)	/* <- t07 */\n\t"\
			"movq	%[__o4]	,%%rsi	\n\t"\
			"movq	%[__o5]	,%%rdi	\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm4\n\t"\
				"mulpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm0		,%%xmm2\n\t"\
				"addpd	%%xmm1		,%%xmm3\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"subpd	%%xmm5		,%%xmm2\n\t"\
				"addpd	%%xmm4		,%%xmm3\n\t"\
				"addpd	%%xmm5		,%%xmm0\n\t"\
				"subpd	%%xmm4		,%%xmm1\n\t"\
				"movaps	%%xmm2		,    (%%rsi)	/* <- t08 */\n\t"\
				"movaps	%%xmm3		,0x10(%%rsi)	/* <- t09 */\n\t"\
				"movaps	%%xmm0		,    (%%rdi)	/* <- t0a */\n\t"\
				"movaps	%%xmm1		,0x10(%%rdi)	/* <- t0b */\n\t"\
			/* Block 3: */\
				"movq	%[__i2]	,%%rax 	\n\t"\
				"movq	%[__i5]	,%%rbx 	\n\t"\
				"movq	%[__i8]	,%%rcx 	\n\t"\
				"movaps	    (%%rbx)	,%%xmm4\n\t"\
				"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
				"movaps	    (%%rax)	,%%xmm0\n\t"\
				"movaps	0x10(%%rax)	,%%xmm1\n\t"\
				"movaps	    (%%rcx)	,%%xmm2\n\t"\
				"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
			"movq	%[__o6]	,%%rsi	\n\t"\
				"subpd	%%xmm2		,%%xmm4\n\t"\
				"subpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm2\n\t"\
				"addpd	%%xmm3		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm2\n\t"\
				"addpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm2		,%%xmm0\n\t"\
				"addpd	%%xmm3		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%rsi)	/* <- t0c */\n\t"\
				"movaps	%%xmm1		,0x10(%%rsi)	/* <- t0d */\n\t"\
			"movq	%[__o7]	,%%rdi	\n\t"\
			"movq	%[__o8]	,%%rsi	\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm4\n\t"\
				"mulpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm0		,%%xmm2\n\t"\
				"addpd	%%xmm1		,%%xmm3\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"subpd	%%xmm5		,%%xmm2\n\t"\
				"addpd	%%xmm4		,%%xmm3\n\t"\
				"addpd	%%xmm5		,%%xmm0\n\t"\
				"subpd	%%xmm4		,%%xmm1\n\t"\
				"movaps	%%xmm2		,    (%%rdi)	/* <- t0e */\n\t"\
				"movaps	%%xmm3		,0x10(%%rdi)	/* <- t0f */\n\t"\
				"movaps	%%xmm0		,    (%%rsi)	/* <- t0g */\n\t"\
				"movaps	%%xmm1		,0x10(%%rsi)	/* <- t0h */\n\t"\
			/*****************************************/\
			/* now do three more radix-3 transforms: */\
			/*****************************************/\
			/* Block 1: */\
			"movq	%[__o0]	,%%rax	\n\t"\
			"movq	%[__o3]	,%%rbx	\n\t"\
			"movq	%[__o6]	,%%rcx	\n\t"\
				"movaps	    (%%rbx)	,%%xmm4\n\t"\
				"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
				"movaps	    (%%rcx)	,%%xmm2\n\t"\
				"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
				"movaps	    (%%rax)	,%%xmm0\n\t"\
				"movaps	0x10(%%rax)	,%%xmm1\n\t"\
				"subpd	%%xmm2		,%%xmm4\n\t"\
				"subpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm2\n\t"\
				"addpd	%%xmm3		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm2\n\t"\
				"addpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm2		,%%xmm0\n\t"\
				"addpd	%%xmm3		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%rax)\n\t"\
				"movaps	%%xmm1		,0x10(%%rax)\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm4\n\t"\
				"mulpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm0		,%%xmm2\n\t"\
				"addpd	%%xmm1		,%%xmm3\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"subpd	%%xmm5		,%%xmm2\n\t"\
				"addpd	%%xmm4		,%%xmm3\n\t"\
				"addpd	%%xmm5		,%%xmm0\n\t"\
				"subpd	%%xmm4		,%%xmm1\n\t"\
				"movaps	%%xmm2		,    (%%rbx)\n\t"\
				"movaps	%%xmm3		,0x10(%%rbx)\n\t"\
				"movaps	%%xmm0		,    (%%rcx)\n\t"\
				"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			/* Block 2: */\
			"movq	%[__o1]	,%%rax	\n\t"\
			"movq	%[__o4]	,%%rbx	\n\t"\
			"movq	%[__o7]	,%%rcx	\n\t"\
				"movaps	    (%%rbx)	,%%xmm2\n\t"\
				"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
				"movaps	    (%%rdx)	,%%xmm6\n\t"/* c1 */\
				"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm0\n\t"\
				"mulpd	%%xmm7		,%%xmm1\n\t"\
				"subpd	%%xmm1		,%%xmm2\n\t"\
				"addpd	%%xmm0		,%%xmm3\n\t"\
				"movaps	    (%%rcx)	,%%xmm4\n\t"\
				"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
				"movaps	0x20(%%rdx)	,%%xmm6\n\t"/* c2 */\
				"movaps	0x30(%%rdx)	,%%xmm7\n\t"\
				"movaps	%%xmm4		,%%xmm0\n\t"\
				"movaps	%%xmm5		,%%xmm1\n\t"\
				"mulpd	%%xmm6		,%%xmm4\n\t"\
				"mulpd	%%xmm6		,%%xmm5\n\t"\
				"mulpd	%%xmm7		,%%xmm0\n\t"\
				"mulpd	%%xmm7		,%%xmm1\n\t"\
				"subpd	%%xmm1		,%%xmm4\n\t"\
				"addpd	%%xmm0		,%%xmm5\n\t"\
				"movaps	0x40(%%rdx)	,%%xmm6\n\t"/* c3m1 */\
				"movaps	0x50(%%rdx)	,%%xmm7\n\t"\
				"movaps	    (%%rax)	,%%xmm0\n\t"\
				"movaps	0x10(%%rax)	,%%xmm1\n\t"\
				"subpd	%%xmm4		,%%xmm2\n\t"\
				"subpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm4\n\t"\
				"addpd	%%xmm5		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm4\n\t"\
				"addpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm4		,%%xmm0\n\t"\
				"addpd	%%xmm5		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%rax)\n\t"\
				"movaps	%%xmm1		,0x10(%%rax)\n\t"\
				"mulpd	%%xmm6		,%%xmm4\n\t"\
				"mulpd	%%xmm6		,%%xmm5\n\t"\
				"mulpd	%%xmm7		,%%xmm2\n\t"\
				"mulpd	%%xmm7		,%%xmm3\n\t"\
				"addpd	%%xmm0		,%%xmm4\n\t"\
				"addpd	%%xmm1		,%%xmm5\n\t"\
				"movaps	%%xmm4		,%%xmm0\n\t"\
				"movaps	%%xmm5		,%%xmm1\n\t"\
				"subpd	%%xmm3		,%%xmm4\n\t"\
				"addpd	%%xmm2		,%%xmm5\n\t"\
				"addpd	%%xmm3		,%%xmm0\n\t"\
				"subpd	%%xmm2		,%%xmm1\n\t"\
				"movaps	%%xmm4		,    (%%rbx)\n\t"\
				"movaps	%%xmm5		,0x10(%%rbx)\n\t"\
				"movaps	%%xmm0		,    (%%rcx)\n\t"\
				"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			/* Block 3: */\
			"movq	%[__o2]	,%%rax	\n\t"\
			"movq	%[__o5]	,%%rbx	\n\t"\
			"movq	%[__o8]	,%%rcx	\n\t"\
				"movaps	    (%%rbx)	,%%xmm2\n\t"\
				"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
				"movaps	0x20(%%rdx)	,%%xmm6\n\t"/* c2 */\
				"movaps	0x30(%%rdx)	,%%xmm7\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm0\n\t"\
				"mulpd	%%xmm7		,%%xmm1\n\t"\
				"subpd	%%xmm1		,%%xmm2\n\t"\
				"addpd	%%xmm0		,%%xmm3\n\t"\
				"movaps	    (%%rcx)	,%%xmm4\n\t"\
				"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
				"movaps	0x60(%%rdx)	,%%xmm6\n\t"/* c4 */\
				"movaps	0x70(%%rdx)	,%%xmm7\n\t"\
				"movaps	%%xmm4		,%%xmm0\n\t"\
				"movaps	%%xmm5		,%%xmm1\n\t"\
				"mulpd	%%xmm6		,%%xmm4\n\t"\
				"mulpd	%%xmm6		,%%xmm5\n\t"\
				"mulpd	%%xmm7		,%%xmm0\n\t"\
				"mulpd	%%xmm7		,%%xmm1\n\t"\
				"subpd	%%xmm1		,%%xmm4\n\t"\
				"addpd	%%xmm0		,%%xmm5\n\t"\
				"movaps	0x40(%%rdx)	,%%xmm6\n\t"/* c3m1 */\
				"movaps	0x50(%%rdx)	,%%xmm7\n\t"\
				"movaps	    (%%rax)	,%%xmm0\n\t"\
				"movaps	0x10(%%rax)	,%%xmm1\n\t"\
				"\n\t"\
				"subpd	%%xmm4		,%%xmm2\n\t"\
				"subpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm4\n\t"\
				"addpd	%%xmm5		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm4\n\t"\
				"addpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm4		,%%xmm0\n\t"\
				"addpd	%%xmm5		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%rax)\n\t"\
				"movaps	%%xmm1		,0x10(%%rax)\n\t"\
				"mulpd	%%xmm6		,%%xmm4\n\t"\
				"mulpd	%%xmm6		,%%xmm5\n\t"\
				"mulpd	%%xmm7		,%%xmm2\n\t"\
				"mulpd	%%xmm7		,%%xmm3\n\t"\
				"addpd	%%xmm0		,%%xmm4\n\t"\
				"addpd	%%xmm1		,%%xmm5\n\t"\
				"\n\t"\
				"movaps	%%xmm4		,%%xmm0\n\t"\
				"movaps	%%xmm5		,%%xmm1\n\t"\
				"subpd	%%xmm3		,%%xmm4\n\t"\
				"addpd	%%xmm2		,%%xmm5\n\t"\
				"addpd	%%xmm3		,%%xmm0\n\t"\
				"subpd	%%xmm2		,%%xmm1\n\t"\
				"movaps	%%xmm4		,    (%%rbx)\n\t"\
				"movaps	%%xmm5		,0x10(%%rbx)\n\t"\
				"movaps	%%xmm0		,    (%%rcx)\n\t"\
				"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
				:					/* outputs: none */\
				: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
				 ,[__i1] "m" (Xi1)\
				 ,[__i2] "m" (Xi2)\
				 ,[__i3] "m" (Xi3)\
				 ,[__i4] "m" (Xi4)\
				 ,[__i5] "m" (Xi5)\
				 ,[__i6] "m" (Xi6)\
				 ,[__i7] "m" (Xi7)\
				 ,[__i8] "m" (Xi8)\
				 ,[__c1] "m" (Xc1)\
				 ,[__o0] "m" (Xo0)\
				 ,[__o1] "m" (Xo1)\
				 ,[__o2] "m" (Xo2)\
				 ,[__o3] "m" (Xo3)\
				 ,[__o4] "m" (Xo4)\
				 ,[__o5] "m" (Xo5)\
				 ,[__o6] "m" (Xo6)\
				 ,[__o7] "m" (Xo7)\
				 ,[__o8] "m" (Xo8)\
				: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
				);\
			}

			/*...Radix-9 DIT: Ins in memory locations __i0-8.\
				Outs at memory locations __o0-8, assumed disjoint with inputs:\
			*/\
			#define SSE2_RADIX_09_DIT(Xi0,Xi1,Xi2,Xi3,Xi4,Xi5,Xi6,Xi7,Xi8, Xc1, Xo0,Xo1,Xo2,Xo3,Xo4,Xo5,Xo6,Xo7,Xo8)\
			{\
			__asm__ volatile (\
			"movq	%[__i0]	,%%rax\n\t"\
			"movq	%[__i1]	,%%rbx\n\t"\
			"movq	%[__i2]	,%%rcx\n\t"\
			"movq	%[__c1]	,%%rdx 	\n\t"/* edx stores trig addresses throughout */\
			/* Block 1: */\
				"movaps	    (%%rbx)	,%%xmm2\n\t"\
				"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
				"movaps	    (%%rax)	,%%xmm0\n\t"\
				"movaps	0x10(%%rax)	,%%xmm1\n\t"\
				"movaps	    (%%rcx)	,%%xmm6\n\t"\
				"movaps	0x10(%%rcx)	,%%xmm7\n\t"\
				"movaps	%%xmm2		,%%xmm4\n\t"\
				"movaps	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm6		,%%xmm2\n\t"\
				"addpd	%%xmm7		,%%xmm3\n\t"\
				"subpd	%%xmm6		,%%xmm4\n\t"\
				"subpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm0\n\t"\
				"addpd	%%xmm3		,%%xmm1\n\t"\
				"movaps	0x40(%%rdx)	,%%xmm6\n\t"/* c3m1 */\
				"movaps	0x50(%%rdx)	,%%xmm7\n\t"\
				"movaps	%%xmm0		,    (%%rax)\n\t"\
				"movaps	%%xmm1		,0x10(%%rax)\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm4\n\t"\
				"mulpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm0		,%%xmm2\n\t"\
				"addpd	%%xmm1		,%%xmm3\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"addpd	%%xmm5		,%%xmm2\n\t"\
				"subpd	%%xmm4		,%%xmm3\n\t"\
				"subpd	%%xmm5		,%%xmm0\n\t"\
				"addpd	%%xmm4		,%%xmm1\n\t"\
				"movaps	%%xmm2		,    (%%rbx)\n\t"\
				"movaps	%%xmm3		,0x10(%%rbx)\n\t"\
				"movaps	%%xmm0		,    (%%rcx)\n\t"\
				"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			/* Block 2: */\
			"movq	%[__i3]	,%%rax\n\t"\
			"movq	%[__i4]	,%%rbx\n\t"\
			"movq	%[__i5]	,%%rcx\n\t"\
				"movaps	    (%%rbx)	,%%xmm4\n\t"\
				"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
				"movaps	    (%%rcx)	,%%xmm2\n\t"\
				"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
				"movaps	    (%%rax)	,%%xmm0\n\t"\
				"movaps	0x10(%%rax)	,%%xmm1\n\t"\
				"subpd	%%xmm2		,%%xmm4\n\t"\
				"subpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm2\n\t"\
				"addpd	%%xmm3		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm2\n\t"\
				"addpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm2		,%%xmm0\n\t"\
				"addpd	%%xmm3		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%rax)\n\t"\
				"movaps	%%xmm1		,0x10(%%rax)\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm4\n\t"\
				"mulpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm0		,%%xmm2\n\t"\
				"addpd	%%xmm1		,%%xmm3\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"addpd	%%xmm5		,%%xmm2\n\t"\
				"subpd	%%xmm4		,%%xmm3\n\t"\
				"subpd	%%xmm5		,%%xmm0\n\t"\
				"addpd	%%xmm4		,%%xmm1\n\t"\
				"movaps	%%xmm2		,    (%%rbx)\n\t"\
				"movaps	%%xmm3		,0x10(%%rbx)\n\t"\
				"movaps	%%xmm0		,    (%%rcx)\n\t"\
				"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			/* Block 3: */\
			"movq	%[__i6]	,%%rax\n\t"\
			"movq	%[__i7]	,%%rbx\n\t"\
			"movq	%[__i8]	,%%rcx\n\t"\
				"movaps	    (%%rbx)	,%%xmm4\n\t"\
				"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
				"movaps	    (%%rcx)	,%%xmm2\n\t"\
				"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
				"movaps	    (%%rax)	,%%xmm0\n\t"\
				"movaps	0x10(%%rax)	,%%xmm1\n\t"\
				"subpd	%%xmm2		,%%xmm4\n\t"\
				"subpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm2\n\t"\
				"addpd	%%xmm3		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm2\n\t"\
				"addpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm2		,%%xmm0\n\t"\
				"addpd	%%xmm3		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%rax)\n\t"\
				"movaps	%%xmm1		,0x10(%%rax)\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm4\n\t"\
				"mulpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm0		,%%xmm2\n\t"\
				"addpd	%%xmm1		,%%xmm3\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"addpd	%%xmm5		,%%xmm2\n\t"\
				"subpd	%%xmm4		,%%xmm3\n\t"\
				"subpd	%%xmm5		,%%xmm0\n\t"\
				"addpd	%%xmm4		,%%xmm1\n\t"\
				"movaps	%%xmm2		,    (%%rbx)\n\t"\
				"movaps	%%xmm3		,0x10(%%rbx)\n\t"\
				"movaps	%%xmm0		,    (%%rcx)\n\t"\
				"movaps	%%xmm1		,0x10(%%rcx)\n\t"\
			/*****************************************/\
			/* now do three more radix-3 transforms: */\
			/*****************************************/\
			/* Block 1: */\
			"movq	%[__o0]	,%%rsi 	\n\t"/* __o0-8: esi,edi store output addresses throughout */\
			"movq	%[__i0]	,%%rax\n\t"\
			"movq	%[__i3]	,%%rbx\n\t"\
			"movq	%[__i6]	,%%rcx\n\t"\
				"movaps	    (%%rbx)	,%%xmm4\n\t"\
				"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
				"movaps	    (%%rcx)	,%%xmm2\n\t"\
				"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
				"movaps	    (%%rax)	,%%xmm0\n\t"\
				"movaps	0x10(%%rax)	,%%xmm1\n\t"\
				"subpd	%%xmm2		,%%xmm4\n\t"\
				"subpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm2\n\t"\
				"addpd	%%xmm3		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm2\n\t"\
				"addpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm2		,%%xmm0\n\t"\
				"addpd	%%xmm3		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%rsi)\n\t"\
				"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm4\n\t"\
				"mulpd	%%xmm7		,%%xmm5\n\t"\
				"addpd	%%xmm0		,%%xmm2\n\t"\
				"addpd	%%xmm1		,%%xmm3\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
			"movq	%[__o3]	,%%rdi 	\n\t"\
			"movq	%[__o6]	,%%rsi 	\n\t"\
				"addpd	%%xmm5		,%%xmm2\n\t"\
				"subpd	%%xmm4		,%%xmm3\n\t"\
				"subpd	%%xmm5		,%%xmm0\n\t"\
				"addpd	%%xmm4		,%%xmm1\n\t"\
				"movaps	%%xmm2		,    (%%rdi)\n\t"\
				"movaps	%%xmm3		,0x10(%%rdi)\n\t"\
				"movaps	%%xmm0		,    (%%rsi)\n\t"\
				"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			/* Block 2: */\
			"movq	%[__i1]	,%%rax\n\t"\
			"movq	%[__i4]	,%%rbx\n\t"\
			"movq	%[__i7]	,%%rcx\n\t"\
				"movaps	    (%%rbx)	,%%xmm2\n\t"\
				"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
				"movaps	    (%%rdx)	,%%xmm6\n\t"/* c1 */\
				"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm0\n\t"\
				"mulpd	%%xmm7		,%%xmm1\n\t"\
				"addpd	%%xmm1		,%%xmm2\n\t"\
				"subpd	%%xmm0		,%%xmm3\n\t"\
				"movaps	    (%%rcx)	,%%xmm4\n\t"\
				"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
				"movaps	0x20(%%rdx)	,%%xmm6\n\t"/* c2 */\
				"movaps	0x30(%%rdx)	,%%xmm7\n\t"\
				"movaps	%%xmm4		,%%xmm0\n\t"\
				"movaps	%%xmm5		,%%xmm1\n\t"\
				"mulpd	%%xmm6		,%%xmm4\n\t"\
				"mulpd	%%xmm6		,%%xmm5\n\t"\
				"mulpd	%%xmm7		,%%xmm0\n\t"\
				"mulpd	%%xmm7		,%%xmm1\n\t"\
				"addpd	%%xmm1		,%%xmm4\n\t"\
				"subpd	%%xmm0		,%%xmm5\n\t"\
				"movaps	0x40(%%rdx)	,%%xmm6\n\t"/* c3m1 */\
				"movaps	0x50(%%rdx)	,%%xmm7\n\t"\
				"movaps	    (%%rax)	,%%xmm0\n\t"\
				"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movq	%[__o7]	,%%rsi 	\n\t"\
				"subpd	%%xmm4		,%%xmm2\n\t"\
				"subpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm4\n\t"\
				"addpd	%%xmm5		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm4\n\t"\
				"addpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm4		,%%xmm0\n\t"\
				"addpd	%%xmm5		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%rsi)\n\t"\
				"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
				"mulpd	%%xmm6		,%%xmm4\n\t"\
				"mulpd	%%xmm6		,%%xmm5\n\t"\
				"mulpd	%%xmm7		,%%xmm2\n\t"\
				"mulpd	%%xmm7		,%%xmm3\n\t"\
				"addpd	%%xmm0		,%%xmm4\n\t"\
				"addpd	%%xmm1		,%%xmm5\n\t"\
			"movq	%[__o1]	,%%rdi 	\n\t"\
			"movq	%[__o4]	,%%rsi 	\n\t"\
				"movaps	%%xmm4		,%%xmm0\n\t"\
				"movaps	%%xmm5		,%%xmm1\n\t"\
				"addpd	%%xmm3		,%%xmm4\n\t"\
				"subpd	%%xmm2		,%%xmm5\n\t"\
				"subpd	%%xmm3		,%%xmm0\n\t"\
				"addpd	%%xmm2		,%%xmm1\n\t"\
				"movaps	%%xmm4		,    (%%rdi)\n\t"\
				"movaps	%%xmm5		,0x10(%%rdi)\n\t"\
				"movaps	%%xmm0		,    (%%rsi)\n\t"\
				"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			/* Block 3: */\
			"movq	%[__i2]	,%%rax\n\t"\
			"movq	%[__i5]	,%%rbx\n\t"\
			"movq	%[__i8]	,%%rcx\n\t"\
				"movaps	    (%%rbx)	,%%xmm2\n\t"\
				"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
				"movaps	0x20(%%rdx)	,%%xmm6\n\t"/* c2 */\
				"movaps	0x30(%%rdx)	,%%xmm7\n\t"\
				"movaps	%%xmm2		,%%xmm0\n\t"\
				"movaps	%%xmm3		,%%xmm1\n\t"\
				"mulpd	%%xmm6		,%%xmm2\n\t"\
				"mulpd	%%xmm6		,%%xmm3\n\t"\
				"mulpd	%%xmm7		,%%xmm0\n\t"\
				"mulpd	%%xmm7		,%%xmm1\n\t"\
				"addpd	%%xmm1		,%%xmm2\n\t"\
				"subpd	%%xmm0		,%%xmm3\n\t"\
				"movaps	    (%%rcx)	,%%xmm4\n\t"\
				"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
				"movaps	0x60(%%rdx)	,%%xmm6\n\t"/* c4 */\
				"movaps	0x70(%%rdx)	,%%xmm7\n\t"\
				"movaps	%%xmm4		,%%xmm0\n\t"\
				"movaps	%%xmm5		,%%xmm1\n\t"\
				"mulpd	%%xmm6		,%%xmm4\n\t"\
				"mulpd	%%xmm6		,%%xmm5\n\t"\
				"mulpd	%%xmm7		,%%xmm0\n\t"\
				"mulpd	%%xmm7		,%%xmm1\n\t"\
				"addpd	%%xmm1		,%%xmm4\n\t"\
				"subpd	%%xmm0		,%%xmm5\n\t"\
				"movaps	0x40(%%rdx)	,%%xmm6\n\t"\
				"movaps	0x50(%%rdx)	,%%xmm7\n\t"\
				"movaps	    (%%rax)	,%%xmm0\n\t"\
				"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movq	%[__o5]	,%%rsi 	\n\t"\
				"subpd	%%xmm4		,%%xmm2\n\t"\
				"subpd	%%xmm5		,%%xmm3\n\t"\
				"addpd	%%xmm4		,%%xmm4\n\t"\
				"addpd	%%xmm5		,%%xmm5\n\t"\
				"addpd	%%xmm2		,%%xmm4\n\t"\
				"addpd	%%xmm3		,%%xmm5\n\t"\
				"addpd	%%xmm4		,%%xmm0\n\t"\
				"addpd	%%xmm5		,%%xmm1\n\t"\
				"movaps	%%xmm0		,    (%%rsi)\n\t"\
				"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
				"mulpd	%%xmm6		,%%xmm4\n\t"\
				"mulpd	%%xmm6		,%%xmm5\n\t"\
				"mulpd	%%xmm7		,%%xmm2\n\t"\
				"mulpd	%%xmm7		,%%xmm3\n\t"\
				"addpd	%%xmm0		,%%xmm4\n\t"\
				"addpd	%%xmm1		,%%xmm5\n\t"\
			"movq	%[__o8]	,%%rdi 	\n\t"\
			"movq	%[__o2]	,%%rsi 	\n\t"\
				"movaps	%%xmm4		,%%xmm0\n\t"\
				"movaps	%%xmm5		,%%xmm1\n\t"\
				"addpd	%%xmm3		,%%xmm4\n\t"\
				"subpd	%%xmm2		,%%xmm5\n\t"\
				"subpd	%%xmm3		,%%xmm0\n\t"\
				"addpd	%%xmm2		,%%xmm1\n\t"\
				"movaps	%%xmm4		,    (%%rdi)\n\t"\
				"movaps	%%xmm5		,0x10(%%rdi)\n\t"\
				"movaps	%%xmm0		,    (%%rsi)\n\t"\
				"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
				:					/* outputs: none */\
				: [__i0] "m" (Xi0)	/* All inputs from memory addresses here */\
				 ,[__i1] "m" (Xi1)\
				 ,[__i2] "m" (Xi2)\
				 ,[__i3] "m" (Xi3)\
				 ,[__i4] "m" (Xi4)\
				 ,[__i5] "m" (Xi5)\
				 ,[__i6] "m" (Xi6)\
				 ,[__i7] "m" (Xi7)\
				 ,[__i8] "m" (Xi8)\
				 ,[__c1] "m" (Xc1)\
				 ,[__o0] "m" (Xo0)\
				 ,[__o1] "m" (Xo1)\
				 ,[__o2] "m" (Xo2)\
				 ,[__o3] "m" (Xo3)\
				 ,[__o4] "m" (Xo4)\
				 ,[__o5] "m" (Xo5)\
				 ,[__o6] "m" (Xo6)\
				 ,[__o7] "m" (Xo7)\
				 ,[__o8] "m" (Xo8)\
				: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
				);\
			}

		#endif	// AVX and 64-bit SSE2

	  #endif	// IF(GCC), USE 32/64-BIT ASM STYLE

	#endif	// MSVC / GCC

#endif	/* radix09_sse_macro_h_included */

