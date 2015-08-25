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

/*******************************************************************************
   We now include this header file if it was not included before.
*******************************************************************************/
#ifndef radix32_ditN_cy_dif1_win_h_included
#define radix32_ditN_cy_dif1_win_h_included

	#define SSE2_RADIX32_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp10,Xp18,Xr00,Xisrt2,Xcc0)\
	{\
		/* SSE2_RADIX4_DIF_IN_PLACE(r00,r20,r10,r30): */\
		\
		__asm	mov	edi,Xisrt2\
		__asm	mov	eax,Xr00\
		__asm	mov	ebx,0x200\
		__asm	mov	ecx,0x100\
		__asm	mov	edx,0x300\
		__asm	add	ebx,eax\
		__asm	add	ecx,eax\
		__asm	add	edx,eax\
		\
		__asm	movaps	xmm0,[eax      ]\
		__asm	movaps	xmm1,[eax+0x010]\
		__asm	movaps	xmm2,[eax      ]\
		__asm	movaps	xmm3,[eax+0x010]\
		__asm	addpd	xmm0,[ebx      ]\
		__asm	addpd	xmm1,[ebx+0x010]\
		__asm	subpd	xmm2,[ebx      ]\
		__asm	subpd	xmm3,[ebx+0x010]\
		__asm	movaps	xmm4,[ecx      ]\
		__asm	movaps	xmm5,[ecx+0x010]\
		__asm	movaps	xmm6,[ecx      ]\
		__asm	movaps	xmm7,[ecx+0x010]\
		__asm	addpd	xmm4,[edx      ]\
		__asm	addpd	xmm5,[edx+0x010]\
		__asm	subpd	xmm6,[edx      ]\
		__asm	subpd	xmm7,[edx+0x010]\
		/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		__asm	subpd	xmm0,xmm4\
		__asm	subpd	xmm1,xmm5\
		__asm	movaps	[ebx      ],xmm0\
		__asm	movaps	[ebx+0x010],xmm1\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm5,xmm1\
		__asm	movaps	[eax      ],xmm4\
		__asm	movaps	[eax+0x010],xmm5\
		\
		__asm	subpd	xmm2,xmm7\
		__asm	subpd	xmm3,xmm6\
		__asm	movaps	[ecx      ],xmm2\
		__asm	movaps	[edx+0x010],xmm3\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm7,xmm2\
		__asm	addpd	xmm6,xmm3\
		__asm	movaps	[edx      ],xmm7\
		__asm	movaps	[ecx+0x010],xmm6\
		\
		/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r08,r28,r18,r38): */\
		__asm	add	eax,0x080\
		__asm	add	ebx,0x080\
		__asm	add	ecx,0x080\
		__asm	add	edx,0x080\
		__asm	movaps	xmm0,[eax      ]\
		__asm	movaps	xmm1,[eax+0x010]\
		__asm	movaps	xmm2,[eax      ]\
		__asm	movaps	xmm3,[eax+0x010]\
		__asm	addpd	xmm0,[ebx      ]\
		__asm	addpd	xmm1,[ebx+0x010]\
		__asm	subpd	xmm2,[ebx      ]\
		__asm	subpd	xmm3,[ebx+0x010]\
		__asm	movaps	xmm4,[ecx      ]\
		__asm	movaps	xmm5,[ecx+0x010]\
		__asm	movaps	xmm6,[ecx      ]\
		__asm	movaps	xmm7,[ecx+0x010]\
		__asm	addpd	xmm6,[edx      ]\
		__asm	addpd	xmm7,[edx+0x010]\
		__asm	subpd	xmm4,[edx      ]\
		__asm	subpd	xmm5,[edx+0x010]\
		/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		__asm	subpd	xmm0,xmm6\
		__asm	subpd	xmm1,xmm7\
		__asm	subpd	xmm2,xmm5\
		__asm	movaps	[ebx      ],xmm0\
		__asm	movaps	[ebx+0x010],xmm1\
		__asm	subpd	xmm3,xmm4\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm6,xmm0\
		__asm	addpd	xmm5,xmm2\
		__asm	addpd	xmm7,xmm1\
		__asm	addpd	xmm4,xmm3\
		__asm	movaps	[eax      ],xmm6\
		__asm	movaps	[eax+0x010],xmm7\
		\
		__asm	movaps	xmm0,[edi]/* isrt2 */\
		__asm	movaps	xmm6,xmm2\
		__asm	movaps	xmm7,xmm5\
		__asm	subpd	xmm2,xmm4\
		__asm	subpd	xmm5,xmm3\
		__asm	addpd	xmm6,xmm4\
		__asm	addpd	xmm7,xmm3\
		__asm	mulpd	xmm2,xmm0\
		__asm	mulpd	xmm5,xmm0\
		__asm	mulpd	xmm6,xmm0\
		__asm	mulpd	xmm7,xmm0\
		__asm	movaps	[ecx      ],xmm2\
		__asm	movaps	[edx      ],xmm5\
		__asm	movaps	[ecx+0x010],xmm6\
		__asm	movaps	[edx+0x010],xmm7\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r00,r10,r20,r30,r08,r18,r28,r38): */\
		__asm	sub	eax,0x080	/* r00 */\
		__asm	sub	ebx,0x200	/* r08 */\
		__asm	add	ecx,0x080	/* r20 */\
		__asm	sub	edx,0x100	/* r28 */\
		__asm	movaps	xmm0,[eax     ]\
		__asm	movaps	xmm4,[ecx     ]\
		__asm	movaps	xmm1,[eax+0x10]\
		__asm	movaps	xmm5,[ecx+0x10]\
		__asm	movaps	xmm2,[ebx     ]\
		__asm	movaps	xmm7,[edx+0x10]\
		__asm	movaps	xmm3,[ebx+0x10]\
		__asm	movaps	xmm6,[edx     ]\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm4,xmm7\
		__asm	subpd	xmm1,xmm3\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm7,xmm4\
		__asm	addpd	xmm3,xmm1\
		__asm	addpd	xmm6,xmm5\
		\
		__asm	movaps	[ebx     ],xmm0\
		__asm	movaps	[ecx     ],xmm4\
		__asm	movaps	[ebx+0x10],xmm1\
		__asm	movaps	[edx+0x10],xmm5\
		__asm	movaps	[eax     ],xmm2\
		__asm	movaps	[edx     ],xmm7\
		__asm	movaps	[eax+0x10],xmm3\
		__asm	movaps	[ecx+0x10],xmm6\
		\
		__asm	add	eax,0x100	/* r10 */\
		__asm	add	ebx,0x100	/* r18 */\
		__asm	add	ecx,0x100	/* r30 */\
		__asm	add	edx,0x100	/* r38 */\
		\
		__asm	movaps	xmm0,[eax     ]\
		__asm	movaps	xmm4,[ecx     ]\
		__asm	movaps	xmm1,[eax+0x10]\
		__asm	movaps	xmm5,[ecx+0x10]\
		__asm	movaps	xmm2,[ebx     ]\
		__asm	movaps	xmm7,[edx+0x10]\
		__asm	movaps	xmm3,[ebx+0x10]\
		__asm	movaps	xmm6,[edx     ]\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm4,xmm7\
		__asm	subpd	xmm1,xmm3\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm7,xmm4\
		__asm	addpd	xmm3,xmm1\
		__asm	addpd	xmm6,xmm5\
		\
		__asm	movaps	[ebx     ],xmm0\
		__asm	movaps	[ecx     ],xmm4\
		__asm	movaps	[ebx+0x10],xmm1\
		__asm	movaps	[edx+0x10],xmm5\
		__asm	movaps	[eax     ],xmm2\
		__asm	movaps	[edx     ],xmm7\
		__asm	movaps	[eax+0x10],xmm3\
		__asm	movaps	[ecx+0x10],xmm6\
		\
		/* SSE2_RADIX4_DIF_IN_PLACE(r04,r24,r14,r34): */\
		\
		__asm	sub	eax,0x0C0	/* r04 */\
		__asm	add	ebx,0x0C0	/* r24 */\
		__asm	sub	ecx,0x1C0	/* r14 */\
		__asm	sub	edx,0x040	/* r34 */\
		\
		__asm	movaps	xmm0,[eax      ]\
		__asm	movaps	xmm1,[eax+0x010]\
		__asm	movaps	xmm2,[eax      ]\
		__asm	movaps	xmm3,[eax+0x010]\
		__asm	addpd	xmm0,[ebx      ]\
		__asm	addpd	xmm1,[ebx+0x010]\
		__asm	subpd	xmm2,[ebx      ]\
		__asm	subpd	xmm3,[ebx+0x010]\
		__asm	movaps	xmm4,[ecx      ]\
		__asm	movaps	xmm5,[ecx+0x010]\
		__asm	movaps	xmm6,[ecx      ]\
		__asm	movaps	xmm7,[ecx+0x010]\
		__asm	addpd	xmm4,[edx      ]\
		__asm	addpd	xmm5,[edx+0x010]\
		__asm	subpd	xmm6,[edx      ]\
		__asm	subpd	xmm7,[edx+0x010]\
		/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		__asm	subpd	xmm0,xmm4\
		__asm	subpd	xmm1,xmm5\
		__asm	movaps	[ebx      ],xmm0\
		__asm	movaps	[ebx+0x010],xmm1\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm5,xmm1\
		__asm	movaps	[eax      ],xmm4\
		__asm	movaps	[eax+0x010],xmm5\
		\
		__asm	subpd	xmm2,xmm7\
		__asm	subpd	xmm3,xmm6\
		__asm	movaps	[ecx      ],xmm2\
		__asm	movaps	[edx+0x010],xmm3\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm7,xmm2\
		__asm	addpd	xmm6,xmm3\
		__asm	movaps	[edx      ],xmm7\
		__asm	movaps	[ecx+0x010],xmm6\
		\
		/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0C,r2C,r1C,r3C): */\
		__asm	add	eax,0x080\
		__asm	add	ebx,0x080\
		__asm	add	ecx,0x080\
		__asm	add	edx,0x080\
		__asm	movaps	xmm0,[eax      ]\
		__asm	movaps	xmm1,[eax+0x010]\
		__asm	movaps	xmm2,[eax      ]\
		__asm	movaps	xmm3,[eax+0x010]\
		__asm	addpd	xmm0,[ebx      ]\
		__asm	addpd	xmm1,[ebx+0x010]\
		__asm	subpd	xmm2,[ebx      ]\
		__asm	subpd	xmm3,[ebx+0x010]\
		__asm	movaps	xmm4,[ecx      ]\
		__asm	movaps	xmm5,[ecx+0x010]\
		__asm	movaps	xmm6,[ecx      ]\
		__asm	movaps	xmm7,[ecx+0x010]\
		__asm	addpd	xmm6,[edx      ]\
		__asm	addpd	xmm7,[edx+0x010]\
		__asm	subpd	xmm4,[edx      ]\
		__asm	subpd	xmm5,[edx+0x010]\
		/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		__asm	subpd	xmm0,xmm6\
		__asm	subpd	xmm1,xmm7\
		__asm	subpd	xmm2,xmm5\
		__asm	movaps	[ebx      ],xmm0\
		__asm	movaps	[ebx+0x010],xmm1\
		__asm	subpd	xmm3,xmm4\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm6,xmm0\
		__asm	addpd	xmm5,xmm2\
		__asm	addpd	xmm7,xmm1\
		__asm	addpd	xmm4,xmm3\
		__asm	movaps	[eax      ],xmm6\
		__asm	movaps	[eax+0x010],xmm7\
		\
		__asm	movaps	xmm0,[edi]/* isrt2 */\
		__asm	movaps	xmm6,xmm2\
		__asm	movaps	xmm7,xmm5\
		__asm	subpd	xmm2,xmm4\
		__asm	subpd	xmm5,xmm3\
		__asm	addpd	xmm6,xmm4\
		__asm	addpd	xmm7,xmm3\
		__asm	mulpd	xmm2,xmm0\
		__asm	mulpd	xmm5,xmm0\
		__asm	mulpd	xmm6,xmm0\
		__asm	mulpd	xmm7,xmm0\
		__asm	movaps	[ecx      ],xmm2\
		__asm	movaps	[edx      ],xmm5\
		__asm	movaps	[ecx+0x010],xmm6\
		__asm	movaps	[edx+0x010],xmm7\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r04,r14,r24,r34,r0C,r1C,r2C,r3C): */\
		__asm	sub	eax,0x080	/* r04 */\
		__asm	sub	ebx,0x200	/* r0C */\
		__asm	add	ecx,0x080	/* r24 */\
		__asm	sub	edx,0x100	/* r2C */\
		__asm	movaps	xmm0,[eax     ]\
		__asm	movaps	xmm4,[ecx     ]\
		__asm	movaps	xmm1,[eax+0x10]\
		__asm	movaps	xmm5,[ecx+0x10]\
		__asm	movaps	xmm2,[ebx     ]\
		__asm	movaps	xmm7,[edx+0x10]\
		__asm	movaps	xmm3,[ebx+0x10]\
		__asm	movaps	xmm6,[edx     ]\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm4,xmm7\
		__asm	subpd	xmm1,xmm3\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm7,xmm4\
		__asm	addpd	xmm3,xmm1\
		__asm	addpd	xmm6,xmm5\
		\
		__asm	movaps	[ebx     ],xmm0\
		__asm	movaps	[ecx     ],xmm4\
		__asm	movaps	[ebx+0x10],xmm1\
		__asm	movaps	[edx+0x10],xmm5\
		__asm	movaps	[eax     ],xmm2\
		__asm	movaps	[edx     ],xmm7\
		__asm	movaps	[eax+0x10],xmm3\
		__asm	movaps	[ecx+0x10],xmm6\
		\
		__asm	add	eax,0x100	/* r14 */\
		__asm	add	ebx,0x100	/* r1C */\
		__asm	add	ecx,0x100	/* r34 */\
		__asm	add	edx,0x100	/* r3C */\
		\
		__asm	movaps	xmm0,[eax     ]\
		__asm	movaps	xmm4,[ecx     ]\
		__asm	movaps	xmm1,[eax+0x10]\
		__asm	movaps	xmm5,[ecx+0x10]\
		__asm	movaps	xmm2,[ebx     ]\
		__asm	movaps	xmm7,[edx+0x10]\
		__asm	movaps	xmm3,[ebx+0x10]\
		__asm	movaps	xmm6,[edx     ]\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm4,xmm7\
		__asm	subpd	xmm1,xmm3\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm7,xmm4\
		__asm	addpd	xmm3,xmm1\
		__asm	addpd	xmm6,xmm5\
		\
		__asm	movaps	[ebx     ],xmm0\
		__asm	movaps	[ecx     ],xmm4\
		__asm	movaps	[ebx+0x10],xmm1\
		__asm	movaps	[edx+0x10],xmm5\
		__asm	movaps	[eax     ],xmm2\
		__asm	movaps	[edx     ],xmm7\
		__asm	movaps	[eax+0x10],xmm3\
		__asm	movaps	[ecx+0x10],xmm6\
		\
		/* SSE2_RADIX4_DIF_IN_PLACE(r02,r22,r12,r32): */\
		\
		__asm	sub	eax,0x120	/* r02 */\
		__asm	add	ebx,0x060	/* r22 */\
		__asm	sub	ecx,0x220	/* r12 */\
		__asm	sub	edx,0x0a0	/* r32 */\
		\
		__asm	movaps	xmm0,[eax      ]\
		__asm	movaps	xmm1,[eax+0x010]\
		__asm	movaps	xmm2,[eax      ]\
		__asm	movaps	xmm3,[eax+0x010]\
		__asm	addpd	xmm0,[ebx      ]\
		__asm	addpd	xmm1,[ebx+0x010]\
		__asm	subpd	xmm2,[ebx      ]\
		__asm	subpd	xmm3,[ebx+0x010]\
		__asm	movaps	xmm4,[ecx      ]\
		__asm	movaps	xmm5,[ecx+0x010]\
		__asm	movaps	xmm6,[ecx      ]\
		__asm	movaps	xmm7,[ecx+0x010]\
		__asm	addpd	xmm4,[edx      ]\
		__asm	addpd	xmm5,[edx+0x010]\
		__asm	subpd	xmm6,[edx      ]\
		__asm	subpd	xmm7,[edx+0x010]\
		/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		__asm	subpd	xmm0,xmm4\
		__asm	subpd	xmm1,xmm5\
		__asm	movaps	[ebx      ],xmm0\
		__asm	movaps	[ebx+0x010],xmm1\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm5,xmm1\
		__asm	movaps	[eax      ],xmm4\
		__asm	movaps	[eax+0x010],xmm5\
		\
		__asm	subpd	xmm2,xmm7\
		__asm	subpd	xmm3,xmm6\
		__asm	movaps	[ecx      ],xmm2\
		__asm	movaps	[edx+0x010],xmm3\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm7,xmm2\
		__asm	addpd	xmm6,xmm3\
		__asm	movaps	[edx      ],xmm7\
		__asm	movaps	[ecx+0x010],xmm6\
		\
		/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0A,r2A,r1A,r3A): */\
		__asm	add	eax,0x080\
		__asm	add	ebx,0x080\
		__asm	add	ecx,0x080\
		__asm	add	edx,0x080\
		__asm	movaps	xmm0,[eax      ]\
		__asm	movaps	xmm1,[eax+0x010]\
		__asm	movaps	xmm2,[eax      ]\
		__asm	movaps	xmm3,[eax+0x010]\
		__asm	addpd	xmm0,[ebx      ]\
		__asm	addpd	xmm1,[ebx+0x010]\
		__asm	subpd	xmm2,[ebx      ]\
		__asm	subpd	xmm3,[ebx+0x010]\
		__asm	movaps	xmm4,[ecx      ]\
		__asm	movaps	xmm5,[ecx+0x010]\
		__asm	movaps	xmm6,[ecx      ]\
		__asm	movaps	xmm7,[ecx+0x010]\
		__asm	addpd	xmm6,[edx      ]\
		__asm	addpd	xmm7,[edx+0x010]\
		__asm	subpd	xmm4,[edx      ]\
		__asm	subpd	xmm5,[edx+0x010]\
		/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		__asm	subpd	xmm0,xmm6\
		__asm	subpd	xmm1,xmm7\
		__asm	subpd	xmm2,xmm5\
		__asm	movaps	[ebx      ],xmm0\
		__asm	movaps	[ebx+0x010],xmm1\
		__asm	subpd	xmm3,xmm4\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm6,xmm0\
		__asm	addpd	xmm5,xmm2\
		__asm	addpd	xmm7,xmm1\
		__asm	addpd	xmm4,xmm3\
		__asm	movaps	[eax      ],xmm6\
		__asm	movaps	[eax+0x010],xmm7\
		\
		__asm	movaps	xmm0,[edi]/* isrt2 */\
		__asm	movaps	xmm6,xmm2\
		__asm	movaps	xmm7,xmm5\
		__asm	subpd	xmm2,xmm4\
		__asm	subpd	xmm5,xmm3\
		__asm	addpd	xmm6,xmm4\
		__asm	addpd	xmm7,xmm3\
		__asm	mulpd	xmm2,xmm0\
		__asm	mulpd	xmm5,xmm0\
		__asm	mulpd	xmm6,xmm0\
		__asm	mulpd	xmm7,xmm0\
		__asm	movaps	[ecx      ],xmm2\
		__asm	movaps	[edx      ],xmm5\
		__asm	movaps	[ecx+0x010],xmm6\
		__asm	movaps	[edx+0x010],xmm7\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r02,r12,r22,r32,r0A,r1A,r2A,r3A): */\
		__asm	sub	eax,0x080	/* r02 */\
		__asm	sub	ebx,0x200	/* r0A */\
		__asm	add	ecx,0x080	/* r22 */\
		__asm	sub	edx,0x100	/* r2A */\
		__asm	movaps	xmm0,[eax     ]\
		__asm	movaps	xmm4,[ecx     ]\
		__asm	movaps	xmm1,[eax+0x10]\
		__asm	movaps	xmm5,[ecx+0x10]\
		__asm	movaps	xmm2,[ebx     ]\
		__asm	movaps	xmm7,[edx+0x10]\
		__asm	movaps	xmm3,[ebx+0x10]\
		__asm	movaps	xmm6,[edx     ]\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm4,xmm7\
		__asm	subpd	xmm1,xmm3\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm7,xmm4\
		__asm	addpd	xmm3,xmm1\
		__asm	addpd	xmm6,xmm5\
		\
		__asm	movaps	[ebx     ],xmm0\
		__asm	movaps	[ecx     ],xmm4\
		__asm	movaps	[ebx+0x10],xmm1\
		__asm	movaps	[edx+0x10],xmm5\
		__asm	movaps	[eax     ],xmm2\
		__asm	movaps	[edx     ],xmm7\
		__asm	movaps	[eax+0x10],xmm3\
		__asm	movaps	[ecx+0x10],xmm6\
		\
		__asm	add	eax,0x100	/* r12 */\
		__asm	add	ebx,0x100	/* r1A */\
		__asm	add	ecx,0x100	/* r32 */\
		__asm	add	edx,0x100	/* r3A */\
		\
		__asm	movaps	xmm0,[eax     ]\
		__asm	movaps	xmm4,[ecx     ]\
		__asm	movaps	xmm1,[eax+0x10]\
		__asm	movaps	xmm5,[ecx+0x10]\
		__asm	movaps	xmm2,[ebx     ]\
		__asm	movaps	xmm7,[edx+0x10]\
		__asm	movaps	xmm3,[ebx+0x10]\
		__asm	movaps	xmm6,[edx     ]\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm4,xmm7\
		__asm	subpd	xmm1,xmm3\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm7,xmm4\
		__asm	addpd	xmm3,xmm1\
		__asm	addpd	xmm6,xmm5\
		\
		__asm	movaps	[ebx     ],xmm0\
		__asm	movaps	[ecx     ],xmm4\
		__asm	movaps	[ebx+0x10],xmm1\
		__asm	movaps	[edx+0x10],xmm5\
		__asm	movaps	[eax     ],xmm2\
		__asm	movaps	[edx     ],xmm7\
		__asm	movaps	[eax+0x10],xmm3\
		__asm	movaps	[ecx+0x10],xmm6\
		\
		/* SSE2_RADIX4_DIF_IN_PLACE(r06,r26,r16,r36): */\
		\
		__asm	sub	eax,0x0C0	/* r06 */\
		__asm	add	ebx,0x0C0	/* r26 */\
		__asm	sub	ecx,0x1C0	/* r16 */\
		__asm	sub	edx,0x040	/* r36 */\
		\
		__asm	movaps	xmm0,[eax      ]\
		__asm	movaps	xmm1,[eax+0x010]\
		__asm	movaps	xmm2,[eax      ]\
		__asm	movaps	xmm3,[eax+0x010]\
		__asm	addpd	xmm0,[ebx      ]\
		__asm	addpd	xmm1,[ebx+0x010]\
		__asm	subpd	xmm2,[ebx      ]\
		__asm	subpd	xmm3,[ebx+0x010]\
		__asm	movaps	xmm4,[ecx      ]\
		__asm	movaps	xmm5,[ecx+0x010]\
		__asm	movaps	xmm6,[ecx      ]\
		__asm	movaps	xmm7,[ecx+0x010]\
		__asm	addpd	xmm4,[edx      ]\
		__asm	addpd	xmm5,[edx+0x010]\
		__asm	subpd	xmm6,[edx      ]\
		__asm	subpd	xmm7,[edx+0x010]\
		/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		__asm	subpd	xmm0,xmm4\
		__asm	subpd	xmm1,xmm5\
		__asm	movaps	[ebx      ],xmm0\
		__asm	movaps	[ebx+0x010],xmm1\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm5,xmm1\
		__asm	movaps	[eax      ],xmm4\
		__asm	movaps	[eax+0x010],xmm5\
		\
		__asm	subpd	xmm2,xmm7\
		__asm	subpd	xmm3,xmm6\
		__asm	movaps	[ecx      ],xmm2\
		__asm	movaps	[edx+0x010],xmm3\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm7,xmm2\
		__asm	addpd	xmm6,xmm3\
		__asm	movaps	[edx      ],xmm7\
		__asm	movaps	[ecx+0x010],xmm6\
		\
		/* SSE2_RADIX4_DIF_IN_PLACE_2NDOFTWO(r0E,r2E,r1E,r3E): */\
		__asm	add	eax,0x080\
		__asm	add	ebx,0x080\
		__asm	add	ecx,0x080\
		__asm	add	edx,0x080\
		__asm	movaps	xmm0,[eax      ]\
		__asm	movaps	xmm1,[eax+0x010]\
		__asm	movaps	xmm2,[eax      ]\
		__asm	movaps	xmm3,[eax+0x010]\
		__asm	addpd	xmm0,[ebx      ]\
		__asm	addpd	xmm1,[ebx+0x010]\
		__asm	subpd	xmm2,[ebx      ]\
		__asm	subpd	xmm3,[ebx+0x010]\
		__asm	movaps	xmm4,[ecx      ]\
		__asm	movaps	xmm5,[ecx+0x010]\
		__asm	movaps	xmm6,[ecx      ]\
		__asm	movaps	xmm7,[ecx+0x010]\
		__asm	addpd	xmm6,[edx      ]\
		__asm	addpd	xmm7,[edx+0x010]\
		__asm	subpd	xmm4,[edx      ]\
		__asm	subpd	xmm5,[edx+0x010]\
		/* Finish radix-4 butterfly and store results into temporary-array slots: */\
		__asm	subpd	xmm0,xmm6\
		__asm	subpd	xmm1,xmm7\
		__asm	subpd	xmm2,xmm5\
		__asm	movaps	[ebx      ],xmm0\
		__asm	movaps	[ebx+0x010],xmm1\
		__asm	subpd	xmm3,xmm4\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm6,xmm0\
		__asm	addpd	xmm5,xmm2\
		__asm	addpd	xmm7,xmm1\
		__asm	addpd	xmm4,xmm3\
		__asm	movaps	[eax      ],xmm6\
		__asm	movaps	[eax+0x010],xmm7\
		\
		__asm	movaps	xmm0,[edi]/* isrt2 */\
		__asm	movaps	xmm6,xmm2\
		__asm	movaps	xmm7,xmm5\
		__asm	subpd	xmm2,xmm4\
		__asm	subpd	xmm5,xmm3\
		__asm	addpd	xmm6,xmm4\
		__asm	addpd	xmm7,xmm3\
		__asm	mulpd	xmm2,xmm0\
		__asm	mulpd	xmm5,xmm0\
		__asm	mulpd	xmm6,xmm0\
		__asm	mulpd	xmm7,xmm0\
		__asm	movaps	[ecx      ],xmm2\
		__asm	movaps	[edx      ],xmm5\
		__asm	movaps	[ecx+0x010],xmm6\
		__asm	movaps	[edx+0x010],xmm7\
		/* SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r06,r16,r26,r36,r0E,r1E,r2E,r3E): */\
		__asm	sub	eax,0x080	/* r02 */\
		__asm	sub	ebx,0x200	/* r0A */\
		__asm	add	ecx,0x080	/* r22 */\
		__asm	sub	edx,0x100	/* r2A */\
		__asm	movaps	xmm0,[eax     ]\
		__asm	movaps	xmm4,[ecx     ]\
		__asm	movaps	xmm1,[eax+0x10]\
		__asm	movaps	xmm5,[ecx+0x10]\
		__asm	movaps	xmm2,[ebx     ]\
		__asm	movaps	xmm7,[edx+0x10]\
		__asm	movaps	xmm3,[ebx+0x10]\
		__asm	movaps	xmm6,[edx     ]\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm4,xmm7\
		__asm	subpd	xmm1,xmm3\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm7,xmm4\
		__asm	addpd	xmm3,xmm1\
		__asm	addpd	xmm6,xmm5\
		\
		__asm	movaps	[ebx     ],xmm0\
		__asm	movaps	[ecx     ],xmm4\
		__asm	movaps	[ebx+0x10],xmm1\
		__asm	movaps	[edx+0x10],xmm5\
		__asm	movaps	[eax     ],xmm2\
		__asm	movaps	[edx     ],xmm7\
		__asm	movaps	[eax+0x10],xmm3\
		__asm	movaps	[ecx+0x10],xmm6\
		\
		__asm	add	eax,0x100	/* r12 */\
		__asm	add	ebx,0x100	/* r1A */\
		__asm	add	ecx,0x100	/* r32 */\
		__asm	add	edx,0x100	/* r3A */\
		\
		__asm	movaps	xmm0,[eax     ]\
		__asm	movaps	xmm4,[ecx     ]\
		__asm	movaps	xmm1,[eax+0x10]\
		__asm	movaps	xmm5,[ecx+0x10]\
		__asm	movaps	xmm2,[ebx     ]\
		__asm	movaps	xmm7,[edx+0x10]\
		__asm	movaps	xmm3,[ebx+0x10]\
		__asm	movaps	xmm6,[edx     ]\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm4,xmm7\
		__asm	subpd	xmm1,xmm3\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm7,xmm4\
		__asm	addpd	xmm3,xmm1\
		__asm	addpd	xmm6,xmm5\
		\
		__asm	movaps	[ebx     ],xmm0\
		__asm	movaps	[ecx     ],xmm4\
		__asm	movaps	[ebx+0x10],xmm1\
		__asm	movaps	[edx+0x10],xmm5\
		__asm	movaps	[eax     ],xmm2\
		__asm	movaps	[edx     ],xmm7\
		__asm	movaps	[eax+0x10],xmm3\
		__asm	movaps	[ecx+0x10],xmm6\
		\
		/*...and now do eight radix-4 transforms, including the internal twiddle factors: */\
		/*...Block 1: t00,t10,t20,t30 in r00,04,02,06 - note swapped middle 2 indices! */\
		\
		__asm	mov	esi,Xr00\
		__asm	mov	eax,Xadd\
		__asm	mov	ebx,Xp01\
		__asm	mov	ecx,Xp02\
		__asm	mov	edx,Xp03\
		__asm	shl	ebx,3\
		__asm	shl	ecx,3\
		__asm	shl	edx,3\
		__asm	add	ebx,eax\
		__asm	add	ecx,eax\
		__asm	add	edx,eax\
		__asm	movaps	xmm0,[esi     ]	/* t00 */\
		__asm	movaps	xmm4,[esi+0x20]	/* t20 */\
		__asm	movaps	xmm1,[esi+0x10]	/* t01 */\
		__asm	movaps	xmm5,[esi+0x30]	/* t21 */\
		__asm	movaps	xmm2,[esi+0x40]	/* t10 */\
		__asm	movaps	xmm6,[esi+0x60]	/* t30 */\
		__asm	movaps	xmm3,[esi+0x50]	/* t11 */\
		__asm	movaps	xmm7,[esi+0x70]	/* t31 */\
		__asm	subpd	xmm0,[esi+0x40]	/* t10=t00-rt */\
		__asm	subpd	xmm4,[esi+0x60]	/* t30=t20-rt */\
		__asm	subpd	xmm1,[esi+0x50]	/* t11=t01-it */\
		__asm	subpd	xmm5,[esi+0x70]	/* t31=t21-it */\
		__asm	addpd	xmm2,[esi     ]	/* t00=t00+rt */\
		__asm	addpd	xmm6,[esi+0x20]	/* t20=t20+rt */\
		__asm	addpd	xmm3,[esi+0x10]	/* t01=t01+it */\
		__asm	addpd	xmm7,[esi+0x30]	/* t21=t21+it */\
		__asm	subpd	xmm2,xmm6	/* t00 <- t00-t20 */\
		__asm	subpd	xmm0,xmm5	/* t10 <- t10-t31 */\
		__asm	subpd	xmm3,xmm7	/* t01 <- t01-t21 */\
		__asm	subpd	xmm1,xmm4	/* t11 <- t11-t30 */\
		__asm	addpd	xmm6,xmm6	/*          2*t20 */\
		__asm	addpd	xmm5,xmm5	/*          2*t31 */\
		__asm	addpd	xmm7,xmm7	/*          2*t21 */\
		__asm	addpd	xmm4,xmm4	/*          2*t30 */\
		__asm	movaps	[ebx     ],xmm2	/* a(jt+p1 ) */\
		__asm	movaps	[ecx     ],xmm0	/* a(jt+p2 ) */\
		__asm	movaps	[ebx+0x10],xmm3	/* a(jp+p1 ) */\
		__asm	movaps	[edx+0x10],xmm1	/* a(jp+p3 ) */\
		__asm	addpd	xmm6,xmm2	/* t20 <- t00+t20 */\
		__asm	addpd	xmm5,xmm0	/* t31 <- t10+t31 */\
		__asm	addpd	xmm7,xmm3	/* t21 <- t01+t21 */\
		__asm	addpd	xmm4,xmm1	/* t30 <- t11+t30 */\
		__asm	movaps	[eax     ],xmm6	/* a(jt+p0 ) */\
		__asm	movaps	[edx     ],xmm5	/* a(jt+p3 ) */\
		__asm	movaps	[eax+0x10],xmm7	/* a(jp+p0 ) */\
		__asm	movaps	[ecx+0x10],xmm4	/* a(jp+p2 ) */\
		/*...Block 5: t08,t18,t28,t38	*/\
		__asm	add	esi,0x80	/* r08 */\
		__asm	mov	edi,Xp04\
		__asm	shl	edi,3\
		__asm	add	eax,edi	/* ap04 */\
		__asm	add	ebx,edi\
		__asm	add	ecx,edi\
		__asm	add	edx,edi\
		__asm	mov	edi,Xisrt2\
		\
		__asm	movaps	xmm4,[esi+0x20]	/* t28 */\
		__asm	movaps	xmm5,[esi+0x30]	/* t29 */\
		__asm	movaps	xmm6,[esi+0x60]	/* t38 */\
		__asm	movaps	xmm7,[esi+0x70]	/* t39 */\
		__asm	movaps	xmm0,[esi     ]	/* t08 */\
		__asm	subpd	xmm4,[esi+0x30]	/* t28-t29 */\
		__asm	movaps	xmm1,[esi+0x10]	/* t09 */\
		__asm	addpd	xmm5,[esi+0x20]	/* t29+t28 */\
		__asm	movaps	xmm2,[esi+0x40]	/* t18 */\
		__asm	mulpd	xmm4,[edi     ]	/* t28 = (t28-t29)*ISRT2 */\
		__asm	movaps	xmm3,[esi+0x50]	/* t19 */\
		__asm	mulpd	xmm5,[edi     ]	/* t29 = (t29+t28)*ISRT2 */\
		__asm	subpd	xmm0,[esi+0x50]	/* t08=t08-t19*/\
		__asm	addpd	xmm6,[esi+0x70]	/* t38+t39 */\
		__asm	subpd	xmm1,[esi+0x40]	/* t19=t09-t18*/\
		__asm	subpd	xmm7,[esi+0x60]	/* t39-t38 */\
		__asm	addpd	xmm2,[esi+0x10]	/* t09=t18+t09*/\
		__asm	mulpd	xmm6,[edi     ]	/*  rt = (t38+t39)*ISRT2 */\
		__asm	addpd	xmm3,[esi     ]	/* t18=t19+t08*/\
		__asm	mulpd	xmm7,[edi     ]	/*  it = (t39-t38)*ISRT2 */\
		\
		__asm	subpd	xmm4,xmm6	/* t28=t28-rt */\
		__asm	subpd	xmm5,xmm7	/* t29=t29-it */\
		__asm	addpd	xmm6,xmm6	/*      2* rt */\
		__asm	addpd	xmm7,xmm7	/*      2* it */\
		__asm	addpd	xmm6,xmm4	/* t38=t28+rt */\
		__asm	addpd	xmm7,xmm5	/* t39=t29+it */\
		__asm	subpd	xmm0,xmm4	/* t08-t28 */\
		__asm	subpd	xmm2,xmm5	/* t09-t29 */\
		__asm	addpd	xmm4,xmm4	/*   2*t28 */\
		__asm	addpd	xmm5,xmm5	/*   2*t29 */\
		\
		__asm	subpd	xmm3,xmm7	/* t18-t39 */\
		__asm	subpd	xmm1,xmm6	/* t19-t38 */\
		__asm	addpd	xmm7,xmm7	/*   2*t39 */\
		__asm	addpd	xmm6,xmm6	/*   2*t38 */\
		__asm	movaps	[ebx     ],xmm0	/* a(jt+p1 ) */\
		__asm	movaps	[ebx+0x10],xmm2	/* a(jp+p1 ) */\
		__asm	movaps	[ecx     ],xmm3	/* a(jt+p2 ) */\
		__asm	movaps	[edx+0x10],xmm1	/* a(jp+p3 ) */\
		__asm	addpd	xmm4,xmm0	/* t08+t28 */\
		__asm	addpd	xmm5,xmm2	/* t09+t29 */\
		__asm	addpd	xmm7,xmm3	/* t18+t39 */\
		__asm	addpd	xmm6,xmm1	/* t19+t38 */\
		__asm	movaps	[eax     ],xmm4	/* a(jt+p0 ) */\
		__asm	movaps	[eax+0x10],xmm5	/* a(jp+p0 ) */\
		__asm	movaps	[edx     ],xmm7	/* a(jt+p3 ) */\
		__asm	movaps	[ecx+0x10],xmm6	/* a(jp+p2 ) */\
		/*...Block 3: t04,t14,t24,t34	*/\
		__asm	add	esi,0x180/* r20 */\
		__asm	mov	edi,Xp04\
		__asm	shl	edi,3\
		__asm	add	eax,edi	/* ap08 */\
		__asm	add	ebx,edi\
		__asm	add	ecx,edi\
		__asm	add	edx,edi\
		__asm	mov	edi,Xcc0\
		\
		__asm	movaps	xmm4,[esi+0x20]	/* t24 */\
		__asm	movaps	xmm6,[esi+0x60]	/* t34 */\
		__asm	movaps	xmm5,[esi+0x30]	/* t25 */\
		__asm	movaps	xmm7,[esi+0x70]	/* t35 */\
		__asm	movaps	xmm0,[esi+0x20]	/* copy t24 */\
		__asm	movaps	xmm2,[esi+0x60]	/* copy t34 */\
		__asm	movaps	xmm1,[esi+0x30]	/* copy t25 */\
		__asm	movaps	xmm3,[esi+0x70]	/* copy t35 */\
		__asm	mulpd	xmm4,[edi     ]	/* t24*c */\
		__asm	mulpd	xmm6,[edi+0x10]	/* t34*s */\
		__asm	mulpd	xmm1,[edi+0x10]	/* t25*s */\
		__asm	mulpd	xmm3,[edi     ]	/* t35*c */\
		__asm	mulpd	xmm5,[edi     ]	/* t25*c */\
		__asm	mulpd	xmm7,[edi+0x10]	/* t35*s */\
		__asm	mulpd	xmm0,[edi+0x10]	/* t24*s */\
		__asm	mulpd	xmm2,[edi     ]	/* t34*c */\
		__asm	subpd	xmm4,xmm1	/* ~t24 */\
		__asm	subpd	xmm6,xmm3	/* rt */\
		__asm	addpd	xmm5,xmm0	/* ~t25 */\
		__asm	addpd	xmm7,xmm2	/* it */\
		__asm	sub	edi,0x10	/* isrt2 */\
		\
		__asm	movaps	xmm2,[esi+0x40]	/* t14 */\
		__asm	subpd	xmm4,xmm6/*~t34=t24-rt */\
		__asm	movaps	xmm3,[esi+0x50]	/* t15 */\
		__asm	subpd	xmm5,xmm7/*~t35=t25-it */\
		__asm	subpd	xmm2,[esi+0x50]	/* t14-t15 */\
		__asm	addpd	xmm6,xmm6/*      2* rt */\
		__asm	addpd	xmm3,[esi+0x40]	/* t15+t14 */\
		__asm	addpd	xmm7,xmm7/*      2* it */\
		__asm	mulpd	xmm2,[edi     ]	/* rt = (t14-t15)*ISRT2 */\
		__asm	addpd	xmm6,xmm4/*~t24=t24+rt */\
		__asm	mulpd	xmm3,[edi     ]	/* it = (t15+t14)*ISRT2 */\
		__asm	addpd	xmm7,xmm5/*~t25=t25+it */\
		__asm	movaps	xmm0,[esi     ]	/* t04 */\
		__asm	movaps	xmm1,[esi+0x10]	/* t05 */\
		__asm	subpd	xmm0,xmm2/*~t14=t04-rt */\
		__asm	subpd	xmm1,xmm3/*~t15=t05-it */\
		__asm	addpd	xmm2,[esi     ]	/*~t04=rt +t04*/\
		__asm	addpd	xmm3,[esi+0x10]	/*~t05=it +t05*/\
		__asm	subpd	xmm2,xmm6	/* t04-t24 */\
		__asm	subpd	xmm0,xmm5	/* t14-t35 */\
		__asm	subpd	xmm3,xmm7	/* t05-t25 */\
		__asm	subpd	xmm1,xmm4	/* t15-t34 */\
		__asm	addpd	xmm6,xmm6	/*   2*t24 */\
		__asm	addpd	xmm5,xmm5	/*          2*t35 */\
		__asm	addpd	xmm7,xmm7	/*   2*t25 */\
		__asm	addpd	xmm4,xmm4	/*          2*t34 */\
		__asm	movaps	[ebx     ],xmm2	/* a(jt+p1 ) */\
		__asm	movaps	[ecx     ],xmm0	/* a(jt+p2 ) */\
		__asm	movaps	[ebx+0x10],xmm3	/* a(jp+p1 ) */\
		__asm	movaps	[edx+0x10],xmm1	/* a(jp+p3 ) */\
		__asm	addpd	xmm6,xmm2	/* t04+t24 */\
		__asm	addpd	xmm5,xmm0	/* t14+t35 */\
		__asm	addpd	xmm7,xmm3	/* t05+t25 */\
		__asm	addpd	xmm4,xmm1	/* t15+t34 */\
		__asm	movaps	[eax     ],xmm6	/* a(jt+p0 ) */\
		__asm	movaps	[edx     ],xmm5	/* a(jt+p3 ) */\
		__asm	movaps	[eax+0x10],xmm7	/* a(jp+p0 ) */\
		__asm	movaps	[ecx+0x10],xmm4	/* a(jp+p2 ) */\
		/*...Block 7: t0C,t1C,t2C,t3C	*/\
		__asm	add	esi,0x80	/* r28 */\
		__asm	mov	edi,Xp04\
		__asm	shl	edi,3\
		__asm	add	eax,edi	/* ap0C */\
		__asm	add	ebx,edi\
		__asm	add	ecx,edi\
		__asm	add	edx,edi\
		__asm	mov	edi,Xcc0\
		\
		__asm	movaps	xmm4,[esi+0x20]	/* t2C */\
		__asm	movaps	xmm6,[esi+0x60]	/* t3C */\
		__asm	movaps	xmm5,[esi+0x30]	/* t2D */\
		__asm	movaps	xmm7,[esi+0x70]	/* t3D */\
		__asm	movaps	xmm0,[esi+0x20]	/* copy t2C */\
		__asm	movaps	xmm2,[esi+0x60]	/* copy t3C */\
		__asm	movaps	xmm1,[esi+0x30]	/* copy t2D */\
		__asm	movaps	xmm3,[esi+0x70]	/* copy t3D */\
		__asm	mulpd	xmm4,[edi+0x10]	/* t2C*s */\
		__asm	mulpd	xmm6,[edi     ]	/* t3C*c */\
		__asm	mulpd	xmm1,[edi     ]	/* t2D*c */\
		__asm	mulpd	xmm3,[edi+0x10]	/* t3D*s */\
		__asm	mulpd	xmm5,[edi+0x10]	/* t2D*s */\
		__asm	mulpd	xmm7,[edi     ]	/* t3D*c */\
		__asm	mulpd	xmm0,[edi     ]	/* t2C*c */\
		__asm	mulpd	xmm2,[edi+0x10]	/* t3C*s */\
		__asm	subpd	xmm4,xmm1	/* ~t24 */\
		__asm	subpd	xmm6,xmm3	/* rt */\
		__asm	addpd	xmm5,xmm0	/* ~t25 */\
		__asm	addpd	xmm7,xmm2	/* it */\
		__asm	sub	edi,0x10	/* isrt2 */\
		\
		__asm	movaps	xmm2,[esi+0x40]	/* t14 */\
		__asm	subpd	xmm4,xmm6	/*~t2C=t2C-rt */\
		__asm	movaps	xmm3,[esi+0x50]	/* t1D */\
		__asm	subpd	xmm5,xmm7	/*~t2D=t2D-it */\
		__asm	addpd	xmm2,[esi+0x50]	/* t1C+t1D */\
		__asm	addpd	xmm6,xmm6	/*      2* rt */\
		__asm	subpd	xmm3,[esi+0x40]	/* t1D-t1C */\
		__asm	addpd	xmm7,xmm7	/*      2* it */\
		__asm	mulpd	xmm2,[edi     ]	/* rt = (t1C+t1D)*ISRT2 */\
		__asm	addpd	xmm6,xmm4	/*~t3C=t2C+rt */\
		__asm	mulpd	xmm3,[edi     ]	/* it = (t1D-t1C)*ISRT2 */\
		__asm	addpd	xmm7,xmm5	/*~t3D=t2D+it */\
		__asm	movaps	xmm0,[esi     ]	/* t0C */\
		__asm	movaps	xmm1,[esi+0x10]	/* t0D */\
		__asm	subpd	xmm0,xmm2	/*~t0C=t0C-rt */\
		__asm	subpd	xmm1,xmm3	/*~t0D=t0D-it */\
		__asm	addpd	xmm2,[esi     ]	/*~t1C=rt +t0C*/\
		__asm	addpd	xmm3,[esi+0x10]	/*~t1D=it +t0D*/\
		__asm	subpd	xmm0,xmm4	/* t0C-t2C */\
		__asm	subpd	xmm2,xmm7	/* t1C-t3D */\
		__asm	subpd	xmm1,xmm5	/* t0D-t2D */\
		__asm	subpd	xmm3,xmm6	/* t1D-t3C */\
		__asm	addpd	xmm4,xmm4	/*   2*t2C */\
		__asm	addpd	xmm7,xmm7	/*   2*t3D */\
		__asm	addpd	xmm5,xmm5	/*   2*t2D */\
		__asm	addpd	xmm6,xmm6	/*   2*t3C */\
		__asm	movaps	[ebx     ],xmm0	/* a(jt+p1 ) */\
		__asm	movaps	[ecx     ],xmm2	/* a(jt+p2 ) */\
		__asm	movaps	[ebx+0x10],xmm1	/* a(jp+p1 ) */\
		__asm	movaps	[edx+0x10],xmm3	/* a(jp+p3 ) */\
		__asm	addpd	xmm4,xmm0	/* t0C+t2C */\
		__asm	addpd	xmm7,xmm2	/* t1C+t3D */\
		__asm	addpd	xmm5,xmm1	/* t0D+t2D */\
		__asm	addpd	xmm6,xmm3	/* t1D+t3C */\
		__asm	movaps	[eax     ],xmm4	/* a(jt+p0 ) */\
		__asm	movaps	[edx     ],xmm7	/* a(jt+p3 ) */\
		__asm	movaps	[eax+0x10],xmm5	/* a(jp+p0 ) */\
		__asm	movaps	[ecx+0x10],xmm6	/* a(jp+p2 ) */\
		/*...Block 2: t02,t12,t22,t32	*/\
		__asm	sub	esi,0x180/* r10 */\
		__asm	mov	edi,Xp04\
		__asm	shl	edi,3\
		__asm	add	eax,edi	/* ap10 */\
		__asm	add	ebx,edi\
		__asm	add	ecx,edi\
		__asm	add	edx,edi\
		__asm	mov	edi,Xcc0\
		\
		__asm	movaps	xmm4,[esi+0x20]	/* t22 */\
		__asm	movaps	xmm6,[esi+0x60]	/* t32 */\
		__asm	movaps	xmm5,[esi+0x30]	/* t23 */\
		__asm	movaps	xmm7,[esi+0x70]	/* t33 */\
		__asm	movaps	xmm0,[esi+0x20]	/* copy t22 */\
		__asm	movaps	xmm2,[esi+0x60]	/* copy t32 */\
		__asm	movaps	xmm1,[esi+0x30]	/* copy t23 */\
		__asm	movaps	xmm3,[esi+0x70]	/* copy t33 */\
		__asm	mulpd	xmm4,[edi+0x20]	/* t22*c32_1 */\
		__asm	mulpd	xmm6,[edi+0x40]	/* t32*c32_3 */\
		__asm	mulpd	xmm1,[edi+0x30]	/* t23*s32_1 */\
		__asm	mulpd	xmm3,[edi+0x50]	/* t33*s32_3 */\
		__asm	mulpd	xmm5,[edi+0x20]	/* t23*c32_1 */\
		__asm	mulpd	xmm7,[edi+0x40]	/* t33*c32_3 */\
		__asm	mulpd	xmm0,[edi+0x30]	/* t22*s32_1 */\
		__asm	mulpd	xmm2,[edi+0x50]	/* t32*s32_3 */\
		__asm	subpd	xmm4,xmm1	/* ~t22 */\
		__asm	subpd	xmm6,xmm3	/* rt */\
		__asm	addpd	xmm5,xmm0	/* ~t23 */\
		__asm	addpd	xmm7,xmm2	/* it */\
		\
		__asm	movaps	xmm2,[esi+0x40]	/* t12 */\
		__asm	movaps	xmm0,[esi+0x50]	/* t13 */\
		__asm	movaps	xmm1,[esi+0x40]	/* copy t12 */\
		__asm	movaps	xmm3,[esi+0x50]	/* copy t13 */\
		__asm	subpd	xmm4,xmm6	/*~t32=t22-rt */\
		__asm	mulpd	xmm2,[edi     ]	/* t12*c */\
		__asm	subpd	xmm5,xmm7	/*~t33=t23-it */\
		__asm	mulpd	xmm0,[edi+0x10]	/* t13*s */\
		__asm	addpd	xmm6,xmm6	/*      2* rt */\
		__asm	mulpd	xmm3,[edi     ]	/* t13*c */\
		__asm	addpd	xmm7,xmm7	/*      2* it */\
		__asm	mulpd	xmm1,[edi+0x10]	/* t12*s */\
		__asm	addpd	xmm6,xmm4	/*~t22=t22+rt */\
		__asm	subpd	xmm2,xmm0	/* rt */\
		__asm	addpd	xmm7,xmm5	/*~t23=t23+it */\
		__asm	addpd	xmm3,xmm1	/* it */\
		\
		__asm	movaps	xmm0,[esi     ]	/* t02 */\
		__asm	movaps	xmm1,[esi+0x10]	/* t03 */\
		__asm	subpd	xmm0,xmm2	/*~t12=t02-rt */\
		__asm	subpd	xmm1,xmm3	/*~t13=t03-it */\
		__asm	addpd	xmm2,[esi     ]	/*~t02=rt+t02 */\
		__asm	addpd	xmm3,[esi+0x10]	/*~t03=it+t03 */\
		__asm	subpd	xmm2,xmm6	/* t02-t22 */\
		__asm	subpd	xmm0,xmm5	/* t12-t33 */\
		__asm	subpd	xmm3,xmm7	/* t03-t23 */\
		__asm	subpd	xmm1,xmm4	/* t13-t32 */\
		__asm	addpd	xmm6,xmm6	/*   2*t22 */\
		__asm	addpd	xmm5,xmm5	/*   2*t33 */\
		__asm	addpd	xmm7,xmm7	/*   2*t23 */\
		__asm	addpd	xmm4,xmm4	/*   2*t32 */\
		__asm	movaps	[ebx     ],xmm2	/* a(jt+p1 ) */\
		__asm	movaps	[ecx     ],xmm0	/* a(jt+p2 ) */\
		__asm	movaps	[ebx+0x10],xmm3	/* a(jp+p1 ) */\
		__asm	movaps	[edx+0x10],xmm1	/* a(jp+p3 ) */\
		__asm	addpd	xmm6,xmm2	/* t02+t22 */\
		__asm	addpd	xmm5,xmm0	/* t12+t33 */\
		__asm	addpd	xmm7,xmm3	/* t03+t23 */\
		__asm	addpd	xmm4,xmm1	/* t13+t32 */\
		__asm	movaps	[eax     ],xmm6	/* a(jt+p0 ) */\
		__asm	movaps	[edx     ],xmm5	/* a(jt+p3 ) */\
		__asm	movaps	[eax+0x10],xmm7	/* a(jp+p0 ) */\
		__asm	movaps	[ecx+0x10],xmm4	/* a(jp+p2 ) */\
		/*...Block 6: t0A,t1A,t2A,t3A	*/\
		__asm	add	esi,0x80	/* r18 */\
		__asm	mov	edi,Xp04\
		__asm	shl	edi,3\
		__asm	add	eax,edi	/* ap14 */\
		__asm	add	ebx,edi\
		__asm	add	ecx,edi\
		__asm	add	edx,edi\
		__asm	mov	edi,Xcc0\
		\
		__asm	movaps	xmm4,[esi+0x20]	/* t2A */\
		__asm	movaps	xmm6,[esi+0x60]	/* t3A */\
		__asm	movaps	xmm5,[esi+0x30]	/* t2B */\
		__asm	movaps	xmm7,[esi+0x70]	/* t3B */\
		__asm	movaps	xmm0,[esi+0x20]	/* copy t2A */\
		__asm	movaps	xmm2,[esi+0x60]	/* copy t3A */\
		__asm	movaps	xmm1,[esi+0x30]	/* copy t2B */\
		__asm	movaps	xmm3,[esi+0x70]	/* copy t3B */\
		__asm	mulpd	xmm4,[edi+0x50]	/* t2A*s32_3 */\
		__asm	mulpd	xmm6,[edi+0x20]	/* t3A*c32_1 */\
		__asm	mulpd	xmm1,[edi+0x40]	/* t2B*c32_3 */\
		__asm	mulpd	xmm3,[edi+0x30]	/* t3B*s32_1 */\
		__asm	mulpd	xmm5,[edi+0x50]	/* t2B*s32_3 */\
		__asm	mulpd	xmm7,[edi+0x20]	/* t3B*c32_1 */\
		__asm	mulpd	xmm0,[edi+0x40]	/* t2A*c32_3 */\
		__asm	mulpd	xmm2,[edi+0x30]	/* t3A*s32_1 */\
		__asm	subpd	xmm4,xmm1	/* ~t2A */\
		__asm	addpd	xmm6,xmm3	/* rt */\
		__asm	addpd	xmm5,xmm0	/* ~t2B */\
		__asm	subpd	xmm7,xmm2	/* it */\
		\
		__asm	movaps	xmm2,[esi+0x40]	/* t1A */\
		__asm	movaps	xmm0,[esi+0x50]	/* t1B */\
		__asm	movaps	xmm1,[esi+0x40]	/* copy t1A */\
		__asm	movaps	xmm3,[esi+0x50]	/* copy t1B */\
		__asm	subpd	xmm4,xmm6	/*~t2A=t2A-rt */\
		__asm	mulpd	xmm2,[edi+0x10]	/* t1A*s */\
		__asm	subpd	xmm5,xmm7	/*~t2B=t2B-it */\
		__asm	mulpd	xmm0,[edi     ]	/* t1B*c */\
		__asm	addpd	xmm6,xmm6	/*      2* rt */\
		__asm	mulpd	xmm3,[edi+0x10]	/* t1B*s */\
		__asm	addpd	xmm7,xmm7	/*      2* it */\
		__asm	mulpd	xmm1,[edi     ]	/* t1A*c */\
		__asm	addpd	xmm6,xmm4	/*~t3A=t2A+rt */\
		__asm	addpd	xmm2,xmm0	/* rt */\
		__asm	addpd	xmm7,xmm5	/*~t3B=t2B+it */\
		__asm	subpd	xmm3,xmm1	/* it */\
		\
		__asm	movaps	xmm0,[esi     ]	/* t0A */\
		__asm	movaps	xmm1,[esi+0x10]	/* t0B */\
		__asm	subpd	xmm0,xmm2	/*~t0A=t0A-rt */\
		__asm	subpd	xmm1,xmm3	/*~t0B=t0B-it */\
		__asm	addpd	xmm2,[esi     ]	/*~t1A=rt+t0A */\
		__asm	addpd	xmm3,[esi+0x10]	/*~t1B=it+t0B */\
		__asm	subpd	xmm0,xmm4	/* t0A-t2A */\
		__asm	subpd	xmm2,xmm7	/* t1A-t3B */\
		__asm	subpd	xmm1,xmm5	/* t0B-t2B */\
		__asm	subpd	xmm3,xmm6	/* t1B-t3A */\
		__asm	addpd	xmm4,xmm4	/*   2*t2A */\
		__asm	addpd	xmm7,xmm7	/*   2*t3B */\
		__asm	addpd	xmm5,xmm5	/*   2*t2B */\
		__asm	addpd	xmm6,xmm6	/*   2*t3A */\
		__asm	movaps	[ebx     ],xmm0	/* a(jt+p1 ) */\
		__asm	movaps	[ecx     ],xmm2	/* a(jt+p2 ) */\
		__asm	movaps	[ebx+0x10],xmm1	/* a(jp+p1 ) */\
		__asm	movaps	[edx+0x10],xmm3	/* a(jp+p3 ) */\
		__asm	addpd	xmm4,xmm0	/* t0A+t2A */\
		__asm	addpd	xmm7,xmm2	/* t1A+t3B */\
		__asm	addpd	xmm5,xmm1	/* t0B+t2B */\
		__asm	addpd	xmm6,xmm3	/* t1B+t3A */\
		__asm	movaps	[eax     ],xmm4	/* a(jt+p0 ) */\
		__asm	movaps	[edx     ],xmm7	/* a(jt+p3 ) */\
		__asm	movaps	[eax+0x10],xmm5	/* a(jp+p0 ) */\
		__asm	movaps	[ecx+0x10],xmm6	/* a(jp+p2 ) */\
		/*...Block 4: t06,t16,t26,t36	*/\
		__asm	add	esi,0x180/* r30 */\
		__asm	mov	edi,Xp04\
		__asm	shl	edi,3\
		__asm	add	eax,edi	/* ap18 */\
		__asm	add	ebx,edi\
		__asm	add	ecx,edi\
		__asm	add	edx,edi\
		__asm	mov	edi,Xcc0\
		\
		__asm	movaps	xmm4,[esi+0x20]	/* t26 */\
		__asm	movaps	xmm6,[esi+0x60]	/* t36 */\
		__asm	movaps	xmm5,[esi+0x30]	/* t27 */\
		__asm	movaps	xmm7,[esi+0x70]	/* t37 */\
		__asm	movaps	xmm0,[esi+0x20]	/* copy t26 */\
		__asm	movaps	xmm2,[esi+0x60]	/* copy t36 */\
		__asm	movaps	xmm1,[esi+0x30]	/* copy t27 */\
		__asm	movaps	xmm3,[esi+0x70]	/* copy t37 */\
		__asm	mulpd	xmm4,[edi+0x40]	/* t26*s32_3 */\
		__asm	mulpd	xmm6,[edi+0x30]	/* t36*s32_1 */\
		__asm	mulpd	xmm1,[edi+0x50]	/* t27*s32_3 */\
		__asm	mulpd	xmm3,[edi+0x20]	/* t37*c32_1 */\
		__asm	mulpd	xmm5,[edi+0x40]	/* t27*c32_3 */\
		__asm	mulpd	xmm7,[edi+0x30]	/* t37*s32_1 */\
		__asm	mulpd	xmm0,[edi+0x50]	/* t26*s32_3 */\
		__asm	mulpd	xmm2,[edi+0x20]	/* t36*c32_1 */\
		__asm	subpd	xmm4,xmm1	/* ~t26 */\
		__asm	addpd	xmm6,xmm3	/* rt */\
		__asm	addpd	xmm5,xmm0	/* ~t27 */\
		__asm	subpd	xmm7,xmm2	/* it */\
		\
		__asm	movaps	xmm2,[esi+0x40]	/* t16 */\
		__asm	movaps	xmm0,[esi+0x50]	/* t17 */\
		__asm	movaps	xmm1,[esi+0x40]	/* copy t16 */\
		__asm	movaps	xmm3,[esi+0x50]	/* copy t17 */\
		__asm	subpd	xmm4,xmm6	/*~t26=t26-rt */\
		__asm	mulpd	xmm2,[edi+0x10]	/* t16*s */\
		__asm	subpd	xmm5,xmm7	/*~t27=t27-it */\
		__asm	mulpd	xmm0,[edi     ]	/* t17*c */\
		__asm	addpd	xmm6,xmm6	/*      2* rt */\
		__asm	mulpd	xmm3,[edi+0x10]	/* t17*s */\
		__asm	addpd	xmm7,xmm7	/*      2* it */\
		__asm	mulpd	xmm1,[edi     ]	/* t16*c */\
		__asm	addpd	xmm6,xmm4	/*~t36=t26+rt */\
		__asm	subpd	xmm2,xmm0	/* rt */\
		__asm	addpd	xmm7,xmm5	/*~t37=t27+it */\
		__asm	addpd	xmm3,xmm1	/* it */\
		\
		__asm	movaps	xmm0,[esi     ]	/* t06 */\
		__asm	movaps	xmm1,[esi+0x10]	/* t07 */\
		__asm	subpd	xmm0,xmm2	/*~t16=t06-rt */\
		__asm	subpd	xmm1,xmm3	/*~t17=t07-it */\
		__asm	addpd	xmm2,[esi     ]	/*~t06=rt+t06 */\
		__asm	addpd	xmm3,[esi+0x10]	/*~t07=it+t07 */\
		__asm	subpd	xmm2,xmm4	/* t06-t26 */\
		__asm	subpd	xmm0,xmm7	/* t16-t37 */\
		__asm	subpd	xmm3,xmm5	/* t07-t27 */\
		__asm	subpd	xmm1,xmm6	/* t17-t36 */\
		__asm	addpd	xmm4,xmm4	/*   2*t26 */\
		__asm	addpd	xmm7,xmm7	/*   2*t37 */\
		__asm	addpd	xmm5,xmm5	/*   2*t27 */\
		__asm	addpd	xmm6,xmm6	/*   2*t36 */\
		__asm	movaps	[ebx     ],xmm2	/* a(jt+p1 ) */\
		__asm	movaps	[ecx     ],xmm0	/* a(jt+p2 ) */\
		__asm	movaps	[ebx+0x10],xmm3	/* a(jp+p1 ) */\
		__asm	movaps	[edx+0x10],xmm1	/* a(jp+p3 ) */\
		__asm	addpd	xmm4,xmm2	/* t06+t26 */\
		__asm	addpd	xmm7,xmm0	/* t16+t37 */\
		__asm	addpd	xmm5,xmm3	/* t07+t27 */\
		__asm	addpd	xmm6,xmm1	/* t17+t36 */\
		__asm	movaps	[eax     ],xmm4	/* a(jt+p0 ) */\
		__asm	movaps	[edx     ],xmm7	/* a(jt+p3 ) */\
		__asm	movaps	[eax+0x10],xmm5	/* a(jp+p0 ) */\
		__asm	movaps	[ecx+0x10],xmm6	/* a(jp+p2 ) */\
		/*...Block 8: t0E,t1E,t2E,t3E	*/\
		__asm	add	esi,0x80	/* r38 */\
		__asm	mov	edi,Xp04\
		__asm	shl	edi,3\
		__asm	add	eax,edi	/* ap1C */\
		__asm	add	ebx,edi\
		__asm	add	ecx,edi\
		__asm	add	edx,edi\
		__asm	mov	edi,Xcc0\
		\
		__asm	movaps	xmm4,[esi+0x20]	/* t2E */\
		__asm	movaps	xmm6,[esi+0x60]	/* t3E */\
		__asm	movaps	xmm5,[esi+0x30]	/* t2F */\
		__asm	movaps	xmm7,[esi+0x70]	/* t3F */\
		__asm	movaps	xmm0,[esi+0x20]	/* copy t2E */\
		__asm	movaps	xmm2,[esi+0x60]	/* copy t3E */\
		__asm	movaps	xmm1,[esi+0x30]	/* copy t2F */\
		__asm	movaps	xmm3,[esi+0x70]	/* copy t3F */\
		__asm	mulpd	xmm4,[edi+0x30]	/* t2E*s32_1 */\
		__asm	mulpd	xmm6,[edi+0x50]	/* t3E*c32_3 */\
		__asm	mulpd	xmm1,[edi+0x20]	/* t2F*c32_1 */\
		__asm	mulpd	xmm3,[edi+0x40]	/* t3F*s32_3 */\
		__asm	mulpd	xmm5,[edi+0x30]	/* t2F*s32_1 */\
		__asm	mulpd	xmm7,[edi+0x50]	/* t3F*c32_3 */\
		__asm	mulpd	xmm0,[edi+0x20]	/* t2E*c32_1 */\
		__asm	mulpd	xmm2,[edi+0x40]	/* t3E*s32_3 */\
		__asm	subpd	xmm4,xmm1	/* ~t2E */\
		__asm	subpd	xmm6,xmm3	/* rt */\
		__asm	addpd	xmm5,xmm0	/* ~t2F */\
		__asm	addpd	xmm7,xmm2	/* it */\
		\
		__asm	movaps	xmm2,[esi+0x40]	/* t1E */\
		__asm	movaps	xmm0,[esi+0x50]	/* t1F */\
		__asm	movaps	xmm1,[esi+0x40]	/* copy t1E */\
		__asm	movaps	xmm3,[esi+0x50]	/* copy t1F */\
		__asm	subpd	xmm4,xmm6	/*~t2E=t2E-rt */\
		__asm	mulpd	xmm2,[edi     ]	/* t1E*c */\
		__asm	subpd	xmm5,xmm7	/*~t2F=t2F-it */\
		__asm	mulpd	xmm0,[edi+0x10]	/* t1F*s */\
		__asm	addpd	xmm6,xmm6	/*      2* rt */\
		__asm	mulpd	xmm3,[edi     ]	/* t1F*c */\
		__asm	addpd	xmm7,xmm7	/*      2* it */\
		__asm	mulpd	xmm1,[edi+0x10]	/* t1E*s */\
		__asm	addpd	xmm6,xmm4	/*~t3E=t2E+rt */\
		__asm	addpd	xmm2,xmm0	/* rt */\
		__asm	addpd	xmm7,xmm5	/*~t3F=t2F+it */\
		__asm	subpd	xmm3,xmm1	/* it */\
		\
		__asm	movaps	xmm0,[esi     ]	/* t0E */\
		__asm	movaps	xmm1,[esi+0x10]	/* t0F */\
		__asm	subpd	xmm0,xmm2	/*~t0E=t0E-rt */\
		__asm	subpd	xmm1,xmm3	/*~t0F=t0F-it */\
		__asm	addpd	xmm2,[esi     ]	/*~t1E=rt+t0E */\
		__asm	addpd	xmm3,[esi+0x10]	/*~t1F=it+t0F */\
		__asm	subpd	xmm0,xmm4	/* t0E-t2E */\
		__asm	subpd	xmm2,xmm7	/* t1E-t3F */\
		__asm	subpd	xmm1,xmm5	/* t0F-t2F */\
		__asm	subpd	xmm3,xmm6	/* t1F-t3E */\
		__asm	addpd	xmm4,xmm4	/*   2*t2E */\
		__asm	addpd	xmm7,xmm7	/*   2*t3F */\
		__asm	addpd	xmm5,xmm5	/*   2*t2F */\
		__asm	addpd	xmm6,xmm6	/*   2*t3E */\
		__asm	movaps	[ebx     ],xmm0	/* a(jt+p1 ) */\
		__asm	movaps	[ecx     ],xmm2	/* a(jt+p2 ) */\
		__asm	movaps	[ebx+0x10],xmm1	/* a(jp+p1 ) */\
		__asm	movaps	[edx+0x10],xmm3	/* a(jp+p3 ) */\
		__asm	addpd	xmm4,xmm0	/* t0E+t2E */\
		__asm	addpd	xmm7,xmm2	/* t1E+t3F */\
		__asm	addpd	xmm5,xmm1	/* t0F+t2F */\
		__asm	addpd	xmm6,xmm3	/* t1F+t3E */\
		__asm	movaps	[eax     ],xmm4	/* a(jt+p0 ) */\
		__asm	movaps	[edx     ],xmm7	/* a(jt+p3 ) */\
		__asm	movaps	[eax+0x10],xmm5	/* a(jp+p0 ) */\
		__asm	movaps	[ecx+0x10],xmm6	/* a(jp+p2 ) */\
	}

	#define SSE2_RADIX32_DIT_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp10,Xp18,Xr00,Xisrt2,Xcc0)\
	{\
		/*...Block 1: */\
		__asm	mov	esi,Xr00\
		__asm	mov	eax,Xadd\
		__asm	mov	ebx,Xp01\
		__asm	mov	ecx,Xp02\
		__asm	mov	edx,Xp03\
		__asm	shl	ebx,3\
		__asm	shl	ecx,3\
		__asm	shl	edx,3\
		__asm	add	ebx,eax\
		__asm	add	ecx,eax\
		__asm	add	edx,eax\
		/* SSE2_RADIX4_DIT_0TWIDDLE[add0, add1, add2, add3, r00): */\
		__asm	movaps	xmm2,[eax     ]\
		__asm	movaps	xmm6,[ecx     ]\
		__asm	movaps	xmm3,[eax+0x10]\
		__asm	movaps	xmm7,[ecx+0x10]\
		__asm	movaps	xmm0,[ebx     ]\
		__asm	movaps	xmm4,[edx     ]\
		__asm	movaps	xmm1,[ebx+0x10]\
		__asm	movaps	xmm5,[edx+0x10]\
		__asm	subpd	xmm2,xmm0\
		__asm	subpd	xmm6,xmm4\
		__asm	subpd	xmm3,xmm1\
		__asm	subpd	xmm7,xmm5\
		__asm	addpd	xmm0,xmm0\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm1,xmm1\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm0,xmm2\
		__asm	addpd	xmm4,xmm6\
		__asm	addpd	xmm1,xmm3\
		__asm	addpd	xmm5,xmm7\
		\
		__asm	subpd	xmm0,xmm4\
		__asm	subpd	xmm2,xmm7\
		__asm	subpd	xmm1,xmm5\
		__asm	subpd	xmm3,xmm6\
		__asm	movaps	[esi+0x040],xmm0\
		__asm	movaps	[esi+0x060],xmm2\
		__asm	movaps	[esi+0x050],xmm1\
		__asm	movaps	[esi+0x030],xmm3\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm7,xmm2\
		__asm	addpd	xmm5,xmm1\
		__asm	addpd	xmm6,xmm3\
		__asm	movaps	[esi      ],xmm4\
		__asm	movaps	[esi+0x020],xmm7\
		__asm	movaps	[esi+0x010],xmm5\
		__asm	movaps	[esi+0x070],xmm6\
		/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r08): */\
		__asm	add	esi,0x80\
		__asm	mov	edi,Xp04\
		__asm	shl	edi,3\
		__asm	add	eax,edi\
		__asm	add	ebx,edi\
		__asm	add	ecx,edi\
		__asm	add	edx,edi\
		\
		__asm	movaps	xmm2,[eax     ]\
		__asm	movaps	xmm6,[ecx     ]\
		__asm	movaps	xmm3,[eax+0x10]\
		__asm	movaps	xmm7,[ecx+0x10]\
		__asm	movaps	xmm0,[ebx     ]\
		__asm	movaps	xmm4,[edx     ]\
		__asm	movaps	xmm1,[ebx+0x10]\
		__asm	movaps	xmm5,[edx+0x10]\
		\
		__asm	mov	edi,Xisrt2\
		__asm	subpd	xmm2,xmm0\
		__asm	subpd	xmm6,xmm4\
		__asm	subpd	xmm3,xmm1\
		__asm	subpd	xmm7,xmm5\
		__asm	addpd	xmm0,xmm0\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm1,xmm1\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm0,xmm2\
		__asm	addpd	xmm4,xmm6\
		__asm	addpd	xmm1,xmm3\
		__asm	addpd	xmm5,xmm7\
		\
		__asm	subpd	xmm0,xmm4\
		__asm	subpd	xmm1,xmm5\
		__asm	movaps	[esi+0x040],xmm0\
		__asm	movaps	[esi+0x050],xmm1\
		__asm	subpd	xmm2,xmm7\
		__asm	addpd	xmm4,xmm4\
		__asm	subpd	xmm3,xmm6\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm5,xmm1\
		__asm	addpd	xmm7,xmm2\
		__asm	movaps	[esi      ],xmm4\
		__asm	addpd	xmm6,xmm3\
		__asm	movaps	[esi+0x010],xmm5\
		\
		__asm	movaps	xmm5,[edi]		/* ISRT2 */\
		__asm	movaps	xmm0,xmm3\
		__asm	movaps	xmm1,xmm6\
		__asm	subpd	xmm3,xmm7\
		__asm	subpd	xmm6,xmm2\
		__asm	addpd	xmm0,xmm7\
		__asm	addpd	xmm1,xmm2\
		__asm	mulpd	xmm3,xmm5\
		__asm	mulpd	xmm6,xmm5\
		__asm	mulpd	xmm0,xmm5\
		__asm	mulpd	xmm1,xmm5\
		__asm	movaps	[esi+0x030],xmm3\
		__asm	movaps	[esi+0x070],xmm6\
		__asm	movaps	[esi+0x020],xmm0\
		__asm	movaps	[esi+0x060],xmm1\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r02,r04,r06,r08,r0A,r0C,r0E): */\
		__asm	movaps	xmm0,[esi-0x080]		/* r00 */\
		__asm	movaps	xmm4,[esi-0x040]		/* r04 */\
		__asm	movaps	xmm1,[esi-0x070]		/* r01 */\
		__asm	movaps	xmm5,[esi-0x030]		/* r05 */\
		__asm	movaps	xmm2,[esi      ]		/* r08 */\
		__asm	movaps	xmm7,[esi+0x050]		/* r0D */\
		__asm	movaps	xmm3,[esi+0x010]		/* r09 */\
		__asm	movaps	xmm6,[esi+0x040]		/* r0C */\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm4,xmm7\
		__asm	subpd	xmm1,xmm3\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm7,xmm4\
		__asm	addpd	xmm3,xmm1\
		__asm	addpd	xmm6,xmm5\
		__asm	movaps	[esi      ],xmm0\
		__asm	movaps	[esi+0x040],xmm4\
		__asm	movaps	[esi+0x010],xmm1\
		__asm	movaps	[esi-0x030],xmm5\
		__asm	movaps	[esi-0x080],xmm2\
		__asm	movaps	[esi-0x040],xmm7\
		__asm	movaps	[esi-0x070],xmm3\
		__asm	movaps	[esi+0x050],xmm6\
		\
		__asm	movaps	xmm0,[esi-0x060]		/* r02 */\
		__asm	movaps	xmm4,[esi-0x020]		/* r06 */\
		__asm	movaps	xmm1,[esi-0x050]		/* r03 */\
		__asm	movaps	xmm5,[esi-0x010]		/* r07 */\
		__asm	movaps	xmm2,[esi+0x020]		/* r0A */\
		__asm	movaps	xmm7,[esi+0x070]		/* r0F */\
		__asm	movaps	xmm3,[esi+0x030]		/* r0B */\
		__asm	movaps	xmm6,[esi+0x060]		/* r0E */\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm4,xmm7\
		__asm	subpd	xmm1,xmm3\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm7,xmm4\
		__asm	addpd	xmm3,xmm1\
		__asm	addpd	xmm6,xmm5\
		__asm	movaps	[esi+0x020],xmm0\
		__asm	movaps	[esi+0x060],xmm4\
		__asm	movaps	[esi+0x030],xmm1\
		__asm	movaps	[esi-0x010],xmm5\
		__asm	movaps	[esi-0x060],xmm2\
		__asm	movaps	[esi-0x020],xmm7\
		__asm	movaps	[esi-0x050],xmm3\
		__asm	movaps	[esi+0x070],xmm6\
		\
		/*...Block 2: */\
		__asm	add	esi,0x80\
		__asm	mov	edi,Xp04\
		__asm	shl	edi,3\
		__asm	add	eax,edi\
		__asm	add	ebx,edi\
		__asm	add	ecx,edi\
		__asm	add	edx,edi\
		\
		/* SSE2_RADIX4_DIT_0TWIDDLE[add0, add1, add2, add3, r10): */\
		__asm	movaps	xmm2,[eax     ]\
		__asm	movaps	xmm6,[ecx     ]\
		__asm	movaps	xmm3,[eax+0x10]\
		__asm	movaps	xmm7,[ecx+0x10]\
		__asm	movaps	xmm0,[ebx     ]\
		__asm	movaps	xmm4,[edx     ]\
		__asm	movaps	xmm1,[ebx+0x10]\
		__asm	movaps	xmm5,[edx+0x10]\
		__asm	subpd	xmm2,xmm0\
		__asm	subpd	xmm6,xmm4\
		__asm	subpd	xmm3,xmm1\
		__asm	subpd	xmm7,xmm5\
		__asm	addpd	xmm0,xmm0\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm1,xmm1\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm0,xmm2\
		__asm	addpd	xmm4,xmm6\
		__asm	addpd	xmm1,xmm3\
		__asm	addpd	xmm5,xmm7\
		\
		__asm	subpd	xmm0,xmm4\
		__asm	subpd	xmm2,xmm7\
		__asm	subpd	xmm1,xmm5\
		__asm	subpd	xmm3,xmm6\
		__asm	movaps	[esi+0x040],xmm0\
		__asm	movaps	[esi+0x060],xmm2\
		__asm	movaps	[esi+0x050],xmm1\
		__asm	movaps	[esi+0x030],xmm3\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm7,xmm2\
		__asm	addpd	xmm5,xmm1\
		__asm	addpd	xmm6,xmm3\
		__asm	movaps	[esi      ],xmm4\
		__asm	movaps	[esi+0x020],xmm7\
		__asm	movaps	[esi+0x010],xmm5\
		__asm	movaps	[esi+0x070],xmm6\
		/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r18): */\
		__asm	add	esi,0x80\
		__asm	mov	edi,Xp04\
		__asm	shl	edi,3\
		__asm	add	eax,edi\
		__asm	add	ebx,edi\
		__asm	add	ecx,edi\
		__asm	add	edx,edi\
		\
		__asm	movaps	xmm2,[eax     ]\
		__asm	movaps	xmm6,[ecx     ]\
		__asm	movaps	xmm3,[eax+0x10]\
		__asm	movaps	xmm7,[ecx+0x10]\
		__asm	movaps	xmm0,[ebx     ]\
		__asm	movaps	xmm4,[edx     ]\
		__asm	movaps	xmm1,[ebx+0x10]\
		__asm	movaps	xmm5,[edx+0x10]\
		\
		__asm	mov	edi,Xisrt2\
		__asm	subpd	xmm2,xmm0\
		__asm	subpd	xmm6,xmm4\
		__asm	subpd	xmm3,xmm1\
		__asm	subpd	xmm7,xmm5\
		__asm	addpd	xmm0,xmm0\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm1,xmm1\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm0,xmm2\
		__asm	addpd	xmm4,xmm6\
		__asm	addpd	xmm1,xmm3\
		__asm	addpd	xmm5,xmm7\
		\
		__asm	subpd	xmm0,xmm4\
		__asm	subpd	xmm1,xmm5\
		__asm	movaps	[esi+0x040],xmm0\
		__asm	movaps	[esi+0x050],xmm1\
		__asm	subpd	xmm2,xmm7\
		__asm	addpd	xmm4,xmm4\
		__asm	subpd	xmm3,xmm6\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm5,xmm1\
		__asm	addpd	xmm7,xmm2\
		__asm	movaps	[esi      ],xmm4\
		__asm	addpd	xmm6,xmm3\
		__asm	movaps	[esi+0x010],xmm5\
		\
		__asm	movaps	xmm5,[edi]		/* ISRT2 */\
		__asm	movaps	xmm0,xmm3\
		__asm	movaps	xmm1,xmm6\
		__asm	subpd	xmm3,xmm7\
		__asm	subpd	xmm6,xmm2\
		__asm	addpd	xmm0,xmm7\
		__asm	addpd	xmm1,xmm2\
		__asm	mulpd	xmm3,xmm5\
		__asm	mulpd	xmm6,xmm5\
		__asm	mulpd	xmm0,xmm5\
		__asm	mulpd	xmm1,xmm5\
		__asm	movaps	[esi+0x030],xmm3\
		__asm	movaps	[esi+0x070],xmm6\
		__asm	movaps	[esi+0x020],xmm0\
		__asm	movaps	[esi+0x060],xmm1\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r10,r12,r14,r16,r18,r1A,r1C,r1E): */\
		__asm	movaps	xmm0,[esi-0x080]		/* r10 */\
		__asm	movaps	xmm4,[esi-0x040]		/* r14 */\
		__asm	movaps	xmm1,[esi-0x070]		/* r11 */\
		__asm	movaps	xmm5,[esi-0x030]		/* r15 */\
		__asm	movaps	xmm2,[esi      ]		/* r18 */\
		__asm	movaps	xmm7,[esi+0x050]		/* r1D */\
		__asm	movaps	xmm3,[esi+0x010]		/* r19 */\
		__asm	movaps	xmm6,[esi+0x040]		/* r1C */\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm4,xmm7\
		__asm	subpd	xmm1,xmm3\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm7,xmm4\
		__asm	addpd	xmm3,xmm1\
		__asm	addpd	xmm6,xmm5\
		__asm	movaps	[esi      ],xmm0\
		__asm	movaps	[esi+0x040],xmm4\
		__asm	movaps	[esi+0x010],xmm1\
		__asm	movaps	[esi-0x030],xmm5\
		__asm	movaps	[esi-0x080],xmm2\
		__asm	movaps	[esi-0x040],xmm7\
		__asm	movaps	[esi-0x070],xmm3\
		__asm	movaps	[esi+0x050],xmm6\
		\
		__asm	movaps	xmm0,[esi-0x060]		/* r12 */\
		__asm	movaps	xmm4,[esi-0x020]		/* r16 */\
		__asm	movaps	xmm1,[esi-0x050]		/* r13 */\
		__asm	movaps	xmm5,[esi-0x010]		/* r17 */\
		__asm	movaps	xmm2,[esi+0x020]		/* r1A */\
		__asm	movaps	xmm7,[esi+0x070]		/* r1F */\
		__asm	movaps	xmm3,[esi+0x030]		/* r1B */\
		__asm	movaps	xmm6,[esi+0x060]		/* r1E */\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm4,xmm7\
		__asm	subpd	xmm1,xmm3\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm7,xmm4\
		__asm	addpd	xmm3,xmm1\
		__asm	addpd	xmm6,xmm5\
		__asm	movaps	[esi+0x020],xmm0\
		__asm	movaps	[esi+0x060],xmm4\
		__asm	movaps	[esi+0x030],xmm1\
		__asm	movaps	[esi-0x010],xmm5\
		__asm	movaps	[esi-0x060],xmm2\
		__asm	movaps	[esi-0x020],xmm7\
		__asm	movaps	[esi-0x050],xmm3\
		__asm	movaps	[esi+0x070],xmm6\
		\
		/*...Block 3: */\
		__asm	add	esi,0x80\
		__asm	mov	edi,Xp04\
		__asm	shl	edi,3\
		__asm	add	eax,edi\
		__asm	add	ebx,edi\
		__asm	add	ecx,edi\
		__asm	add	edx,edi\
		\
		/* SSE2_RADIX4_DIT_0TWIDDLE[add0, add1, add2, add3, r20): */\
		__asm	movaps	xmm2,[eax     ]\
		__asm	movaps	xmm6,[ecx     ]\
		__asm	movaps	xmm3,[eax+0x10]\
		__asm	movaps	xmm7,[ecx+0x10]\
		__asm	movaps	xmm0,[ebx     ]\
		__asm	movaps	xmm4,[edx     ]\
		__asm	movaps	xmm1,[ebx+0x10]\
		__asm	movaps	xmm5,[edx+0x10]\
		__asm	subpd	xmm2,xmm0\
		__asm	subpd	xmm6,xmm4\
		__asm	subpd	xmm3,xmm1\
		__asm	subpd	xmm7,xmm5\
		__asm	addpd	xmm0,xmm0\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm1,xmm1\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm0,xmm2\
		__asm	addpd	xmm4,xmm6\
		__asm	addpd	xmm1,xmm3\
		__asm	addpd	xmm5,xmm7\
		\
		__asm	subpd	xmm0,xmm4\
		__asm	subpd	xmm2,xmm7\
		__asm	subpd	xmm1,xmm5\
		__asm	subpd	xmm3,xmm6\
		__asm	movaps	[esi+0x040],xmm0\
		__asm	movaps	[esi+0x060],xmm2\
		__asm	movaps	[esi+0x050],xmm1\
		__asm	movaps	[esi+0x030],xmm3\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm7,xmm2\
		__asm	addpd	xmm5,xmm1\
		__asm	addpd	xmm6,xmm3\
		__asm	movaps	[esi      ],xmm4\
		__asm	movaps	[esi+0x020],xmm7\
		__asm	movaps	[esi+0x010],xmm5\
		__asm	movaps	[esi+0x070],xmm6\
		/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r28): */\
		__asm	add	esi,0x80\
		__asm	mov	edi,Xp04\
		__asm	shl	edi,3\
		__asm	add	eax,edi\
		__asm	add	ebx,edi\
		__asm	add	ecx,edi\
		__asm	add	edx,edi\
		\
		__asm	movaps	xmm2,[eax     ]\
		__asm	movaps	xmm6,[ecx     ]\
		__asm	movaps	xmm3,[eax+0x10]\
		__asm	movaps	xmm7,[ecx+0x10]\
		__asm	movaps	xmm0,[ebx     ]\
		__asm	movaps	xmm4,[edx     ]\
		__asm	movaps	xmm1,[ebx+0x10]\
		__asm	movaps	xmm5,[edx+0x10]\
		\
		__asm	mov	edi,Xisrt2\
		__asm	subpd	xmm2,xmm0\
		__asm	subpd	xmm6,xmm4\
		__asm	subpd	xmm3,xmm1\
		__asm	subpd	xmm7,xmm5\
		__asm	addpd	xmm0,xmm0\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm1,xmm1\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm0,xmm2\
		__asm	addpd	xmm4,xmm6\
		__asm	addpd	xmm1,xmm3\
		__asm	addpd	xmm5,xmm7\
		\
		__asm	subpd	xmm0,xmm4\
		__asm	subpd	xmm1,xmm5\
		__asm	movaps	[esi+0x040],xmm0\
		__asm	movaps	[esi+0x050],xmm1\
		__asm	subpd	xmm2,xmm7\
		__asm	addpd	xmm4,xmm4\
		__asm	subpd	xmm3,xmm6\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm5,xmm1\
		__asm	addpd	xmm7,xmm2\
		__asm	movaps	[esi      ],xmm4\
		__asm	addpd	xmm6,xmm3\
		__asm	movaps	[esi+0x010],xmm5\
		\
		__asm	movaps	xmm5,[edi]		/* ISRT2 */\
		__asm	movaps	xmm0,xmm3\
		__asm	movaps	xmm1,xmm6\
		__asm	subpd	xmm3,xmm7\
		__asm	subpd	xmm6,xmm2\
		__asm	addpd	xmm0,xmm7\
		__asm	addpd	xmm1,xmm2\
		__asm	mulpd	xmm3,xmm5\
		__asm	mulpd	xmm6,xmm5\
		__asm	mulpd	xmm0,xmm5\
		__asm	mulpd	xmm1,xmm5\
		__asm	movaps	[esi+0x030],xmm3\
		__asm	movaps	[esi+0x070],xmm6\
		__asm	movaps	[esi+0x020],xmm0\
		__asm	movaps	[esi+0x060],xmm1\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r20,r22,r24,r26,r28,r2A,r2C,r2E): */\
		__asm	movaps	xmm0,[esi-0x080]		/* r20 */\
		__asm	movaps	xmm4,[esi-0x040]		/* r24 */\
		__asm	movaps	xmm1,[esi-0x070]		/* r21 */\
		__asm	movaps	xmm5,[esi-0x030]		/* r25 */\
		__asm	movaps	xmm2,[esi      ]		/* r28 */\
		__asm	movaps	xmm7,[esi+0x050]		/* r2D */\
		__asm	movaps	xmm3,[esi+0x010]		/* r29 */\
		__asm	movaps	xmm6,[esi+0x040]		/* r2C */\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm4,xmm7\
		__asm	subpd	xmm1,xmm3\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm7,xmm4\
		__asm	addpd	xmm3,xmm1\
		__asm	addpd	xmm6,xmm5\
		__asm	movaps	[esi      ],xmm0\
		__asm	movaps	[esi+0x040],xmm4\
		__asm	movaps	[esi+0x010],xmm1\
		__asm	movaps	[esi-0x030],xmm5\
		__asm	movaps	[esi-0x080],xmm2\
		__asm	movaps	[esi-0x040],xmm7\
		__asm	movaps	[esi-0x070],xmm3\
		__asm	movaps	[esi+0x050],xmm6\
		\
		__asm	movaps	xmm0,[esi-0x060]		/* r22 */\
		__asm	movaps	xmm4,[esi-0x020]		/* r26 */\
		__asm	movaps	xmm1,[esi-0x050]		/* r23 */\
		__asm	movaps	xmm5,[esi-0x010]		/* r27 */\
		__asm	movaps	xmm2,[esi+0x020]		/* r2A */\
		__asm	movaps	xmm7,[esi+0x070]		/* r2F */\
		__asm	movaps	xmm3,[esi+0x030]		/* r2B */\
		__asm	movaps	xmm6,[esi+0x060]		/* r2E */\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm4,xmm7\
		__asm	subpd	xmm1,xmm3\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm7,xmm4\
		__asm	addpd	xmm3,xmm1\
		__asm	addpd	xmm6,xmm5\
		__asm	movaps	[esi+0x020],xmm0\
		__asm	movaps	[esi+0x060],xmm4\
		__asm	movaps	[esi+0x030],xmm1\
		__asm	movaps	[esi-0x010],xmm5\
		__asm	movaps	[esi-0x060],xmm2\
		__asm	movaps	[esi-0x020],xmm7\
		__asm	movaps	[esi-0x050],xmm3\
		__asm	movaps	[esi+0x070],xmm6\
		\
		/*...Block 4: */\
		__asm	add	esi,0x80\
		__asm	mov	edi,Xp04\
		__asm	shl	edi,3\
		__asm	add	eax,edi\
		__asm	add	ebx,edi\
		__asm	add	ecx,edi\
		__asm	add	edx,edi\
		\
		/* SSE2_RADIX4_DIT_0TWIDDLE[add0, add1, add2, add3, r30): */\
		__asm	movaps	xmm2,[eax     ]\
		__asm	movaps	xmm6,[ecx     ]\
		__asm	movaps	xmm3,[eax+0x10]\
		__asm	movaps	xmm7,[ecx+0x10]\
		__asm	movaps	xmm0,[ebx     ]\
		__asm	movaps	xmm4,[edx     ]\
		__asm	movaps	xmm1,[ebx+0x10]\
		__asm	movaps	xmm5,[edx+0x10]\
		__asm	subpd	xmm2,xmm0\
		__asm	subpd	xmm6,xmm4\
		__asm	subpd	xmm3,xmm1\
		__asm	subpd	xmm7,xmm5\
		__asm	addpd	xmm0,xmm0\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm1,xmm1\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm0,xmm2\
		__asm	addpd	xmm4,xmm6\
		__asm	addpd	xmm1,xmm3\
		__asm	addpd	xmm5,xmm7\
		\
		__asm	subpd	xmm0,xmm4\
		__asm	subpd	xmm2,xmm7\
		__asm	subpd	xmm1,xmm5\
		__asm	subpd	xmm3,xmm6\
		__asm	movaps	[esi+0x040],xmm0\
		__asm	movaps	[esi+0x060],xmm2\
		__asm	movaps	[esi+0x050],xmm1\
		__asm	movaps	[esi+0x030],xmm3\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm7,xmm2\
		__asm	addpd	xmm5,xmm1\
		__asm	addpd	xmm6,xmm3\
		__asm	movaps	[esi      ],xmm4\
		__asm	movaps	[esi+0x020],xmm7\
		__asm	movaps	[esi+0x010],xmm5\
		__asm	movaps	[esi+0x070],xmm6\
		/* SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r38): */\
		__asm	add	esi,0x80\
		__asm	mov	edi,Xp04\
		__asm	shl	edi,3\
		__asm	add	eax,edi\
		__asm	add	ebx,edi\
		__asm	add	ecx,edi\
		__asm	add	edx,edi\
		\
		__asm	movaps	xmm2,[eax     ]\
		__asm	movaps	xmm6,[ecx     ]\
		__asm	movaps	xmm3,[eax+0x10]\
		__asm	movaps	xmm7,[ecx+0x10]\
		__asm	movaps	xmm0,[ebx     ]\
		__asm	movaps	xmm4,[edx     ]\
		__asm	movaps	xmm1,[ebx+0x10]\
		__asm	movaps	xmm5,[edx+0x10]\
		\
		__asm	mov	edi,Xisrt2\
		__asm	subpd	xmm2,xmm0\
		__asm	subpd	xmm6,xmm4\
		__asm	subpd	xmm3,xmm1\
		__asm	subpd	xmm7,xmm5\
		__asm	addpd	xmm0,xmm0\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm1,xmm1\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm0,xmm2\
		__asm	addpd	xmm4,xmm6\
		__asm	addpd	xmm1,xmm3\
		__asm	addpd	xmm5,xmm7\
		\
		__asm	subpd	xmm0,xmm4\
		__asm	subpd	xmm1,xmm5\
		__asm	movaps	[esi+0x040],xmm0\
		__asm	movaps	[esi+0x050],xmm1\
		__asm	subpd	xmm2,xmm7\
		__asm	addpd	xmm4,xmm4\
		__asm	subpd	xmm3,xmm6\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm5,xmm1\
		__asm	addpd	xmm7,xmm2\
		__asm	movaps	[esi      ],xmm4\
		__asm	addpd	xmm6,xmm3\
		__asm	movaps	[esi+0x010],xmm5\
		\
		__asm	movaps	xmm5,[edi]		/* ISRT2 */\
		__asm	movaps	xmm0,xmm3\
		__asm	movaps	xmm1,xmm6\
		__asm	subpd	xmm3,xmm7\
		__asm	subpd	xmm6,xmm2\
		__asm	addpd	xmm0,xmm7\
		__asm	addpd	xmm1,xmm2\
		__asm	mulpd	xmm3,xmm5\
		__asm	mulpd	xmm6,xmm5\
		__asm	mulpd	xmm0,xmm5\
		__asm	mulpd	xmm1,xmm5\
		__asm	movaps	[esi+0x030],xmm3\
		__asm	movaps	[esi+0x070],xmm6\
		__asm	movaps	[esi+0x020],xmm0\
		__asm	movaps	[esi+0x060],xmm1\
		/* SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r30,r32,r34,r36,r38,r3A,r3C,r3E): */\
		__asm	movaps	xmm0,[esi-0x080]	/* r30 */\
		__asm	movaps	xmm4,[esi-0x040]	/* r34 */\
		__asm	movaps	xmm1,[esi-0x070]	/* r31 */\
		__asm	movaps	xmm5,[esi-0x030]	/* r35 */\
		__asm	movaps	xmm2,[esi      ]	/* r38 */\
		__asm	movaps	xmm7,[esi+0x050]	/* r3D */\
		__asm	movaps	xmm3,[esi+0x010]	/* r39 */\
		__asm	movaps	xmm6,[esi+0x040]	/* r3C */\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm4,xmm7\
		__asm	subpd	xmm1,xmm3\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm7,xmm4\
		__asm	addpd	xmm3,xmm1\
		__asm	addpd	xmm6,xmm5\
		__asm	movaps	[esi      ],xmm0\
		__asm	movaps	[esi+0x040],xmm4\
		__asm	movaps	[esi+0x010],xmm1\
		__asm	movaps	[esi-0x030],xmm5\
		__asm	movaps	[esi-0x080],xmm2\
		__asm	movaps	[esi-0x040],xmm7\
		__asm	movaps	[esi-0x070],xmm3\
		__asm	movaps	[esi+0x050],xmm6\
		\
		__asm	movaps	xmm0,[esi-0x060]	/* r32 */\
		__asm	movaps	xmm4,[esi-0x020]	/* r36 */\
		__asm	movaps	xmm1,[esi-0x050]	/* r33 */\
		__asm	movaps	xmm5,[esi-0x010]	/* r37 */\
		__asm	movaps	xmm2,[esi+0x020]	/* r3A */\
		__asm	movaps	xmm7,[esi+0x070]	/* r3F */\
		__asm	movaps	xmm3,[esi+0x030]	/* r3B */\
		__asm	movaps	xmm6,[esi+0x060]	/* r3E */\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm4,xmm7\
		__asm	subpd	xmm1,xmm3\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm7,xmm4\
		__asm	addpd	xmm3,xmm1\
		__asm	addpd	xmm6,xmm5\
		__asm	movaps	[esi+0x020],xmm0\
		__asm	movaps	[esi+0x060],xmm4\
		__asm	movaps	[esi+0x030],xmm1\
		__asm	movaps	[esi-0x010],xmm5\
		__asm	movaps	[esi-0x060],xmm2\
		__asm	movaps	[esi-0x020],xmm7\
		__asm	movaps	[esi-0x050],xmm3\
		__asm	movaps	[esi+0x070],xmm6\
		\
		/*...and now do eight radix-4 transforms, including the internal twiddle factors:,*/\
		/*...Block 1: r00,r10,r20,r30	*/\
		\
		__asm	mov	edi,Xisrt2\
		__asm	mov	esi,Xcc0\
		__asm	mov	eax,Xr00\
		__asm	mov	ebx,eax\
		__asm	mov	ecx,eax\
		__asm	mov	edx,eax\
		__asm	add	ebx,0x100\
		__asm	add	ecx,0x200\
		__asm	add	edx,0x300\
		\
		__asm	movaps	xmm0,[eax     ]\
		__asm	movaps	xmm1,[eax+0x10]\
		__asm	movaps	xmm2,[ebx     ]\
		__asm	movaps	xmm3,[ebx+0x10]\
		\
		__asm	subpd	xmm0,[ebx     ]\
		__asm	subpd	xmm1,[ebx+0x10]\
		__asm	addpd	xmm2,[eax     ]\
		__asm	addpd	xmm3,[eax+0x10]\
		\
		__asm	movaps	xmm4,[ecx     ]\
		__asm	movaps	xmm5,[ecx+0x10]\
		__asm	movaps	xmm6,[edx     ]\
		__asm	movaps	xmm7,[edx+0x10]\
		\
		__asm	subpd	xmm4,[edx     ]\
		__asm	subpd	xmm5,[edx+0x10]\
		__asm	addpd	xmm6,[ecx     ]\
		__asm	addpd	xmm7,[ecx+0x10]\
		\
		/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm6,xmm7,xmm4,xmm5): swap xmm[01<->[23, xmm[45<->[67 */\
		__asm	addpd	xmm2,xmm6\
		__asm	addpd	xmm3,xmm7\
		__asm	movaps	[eax     ],xmm2\
		__asm	movaps	[eax+0x10],xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm7,xmm7\
		__asm	subpd	xmm2,xmm6\
		__asm	subpd	xmm3,xmm7\
		__asm	movaps	[ecx     ],xmm2\
		__asm	movaps	[ecx+0x10],xmm3\
		\
		__asm	addpd	xmm0,xmm5\
		__asm	subpd	xmm1,xmm4\
		__asm	movaps	[ebx     ],xmm0\
		__asm	movaps	[ebx+0x10],xmm1\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm4,xmm4\
		__asm	subpd	xmm0,xmm5\
		__asm	addpd	xmm1,xmm4\
		__asm	movaps	[edx     ],xmm0\
		__asm	movaps	[edx+0x10],xmm1\
		\
		/*...Block 5: r08,r18,r28,r38	*/\
		\
		__asm	add	eax,0x80		/* r08 */\
		__asm	add	ebx,0x80\
		__asm	add	ecx,0x80\
		__asm	add	edx,0x80\
		__asm	movaps	xmm2,[edi]		/* isrt2 */\
		\
		__asm	movaps	xmm4,[ecx     ]\
		__asm	movaps	xmm5,[ecx+0x10]\
		__asm	movaps	xmm0,[edx     ]\
		__asm	movaps	xmm1,[edx+0x10]\
		\
		__asm	addpd	xmm4,[ecx+0x10]\
		__asm	subpd	xmm5,[ecx     ]\
		__asm	subpd	xmm0,[edx+0x10]\
		__asm	addpd	xmm1,[edx     ]\
		__asm	mulpd	xmm4,xmm2\
		__asm	mulpd	xmm5,xmm2\
		__asm	mulpd	xmm0,xmm2\
		__asm	mulpd	xmm1,xmm2\
		__asm	movaps	xmm6,xmm4\
		__asm	movaps	xmm7,xmm5\
		\
		__asm	subpd	xmm4,xmm0\
		__asm	subpd	xmm5,xmm1\
		__asm	addpd	xmm6,xmm0\
		__asm	addpd	xmm7,xmm1\
		\
		__asm	movaps	xmm0,[eax     ]\
		__asm	movaps	xmm1,[eax+0x10]\
		__asm	movaps	xmm2,[ebx     ]\
		__asm	movaps	xmm3,[ebx+0x10]\
		\
		__asm	subpd	xmm0,[ebx+0x10]\
		__asm	subpd	xmm1,[ebx     ]\
		__asm	addpd	xmm3,[eax     ]\
		__asm	addpd	xmm2,[eax+0x10]\
		\
		/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm3,xmm1,xmm0,xmm2,xmm4,xmm5,xmm6,xmm7): swap xmm0123<->3102 */\
		__asm	addpd	xmm3,xmm4\
		__asm	addpd	xmm1,xmm5\
		__asm	movaps	[eax     ],xmm3\
		__asm	movaps	[eax+0x10],xmm1\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm5,xmm5\
		__asm	subpd	xmm3,xmm4\
		__asm	subpd	xmm1,xmm5\
		__asm	movaps	[ecx     ],xmm3\
		__asm	movaps	[ecx+0x10],xmm1\
		\
		__asm	addpd	xmm0,xmm7\
		__asm	subpd	xmm2,xmm6\
		__asm	movaps	[ebx     ],xmm0\
		__asm	movaps	[ebx+0x10],xmm2\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm6,xmm6\
		__asm	subpd	xmm0,xmm7\
		__asm	addpd	xmm2,xmm6\
		__asm	movaps	[edx     ],xmm0\
		__asm	movaps	[edx+0x10],xmm2\
		\
		/*...Block 3: r04,r14,r24,r34	*/\
		\
		__asm	sub	eax ,0x40		/* r04 */\
		__asm	sub	ebx,0x40\
		__asm	sub	ecx,0x40\
		__asm	sub	edx,0x40\
		\
		__asm	movaps	xmm4,[ecx     ]\
		__asm	movaps	xmm0,[edx     ]\
		__asm	movaps	xmm5,[ecx+0x10]\
		__asm	movaps	xmm1,[edx+0x10]\
		__asm	movaps	xmm6,[ecx     ]\
		__asm	movaps	xmm2,[edx     ]\
		__asm	movaps	xmm7,[ecx+0x10]\
		__asm	movaps	xmm3,[edx+0x10]\
		\
		__asm	mulpd	xmm4,[esi     ]\
		__asm	mulpd	xmm0,[esi+0x10]\
		__asm	mulpd	xmm5,[esi     ]\
		__asm	mulpd	xmm1,[esi+0x10]\
		__asm	mulpd	xmm6,[esi+0x10]\
		__asm	mulpd	xmm2,[esi     ]\
		__asm	mulpd	xmm7,[esi+0x10]\
		__asm	mulpd	xmm3,[esi     ]\
		__asm	subpd	xmm5,xmm6\
		__asm	subpd	xmm1,xmm2\
		__asm	addpd	xmm4,xmm7\
		__asm	addpd	xmm0,xmm3\
		__asm	movaps	xmm7,xmm5\
		__asm	movaps	xmm6,xmm4\
		\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm5,xmm1\
		__asm	subpd	xmm6,xmm0\
		__asm	subpd	xmm7,xmm1\
		\
		__asm	movaps	xmm2,[ebx     ]\
		__asm	movaps	xmm3,[ebx+0x10]\
		__asm	movaps	xmm0,[eax     ]\
		__asm	movaps	xmm1,[eax+0x10]\
		__asm	addpd	xmm2,[ebx+0x10]\
		__asm	subpd	xmm3,[ebx     ]\
		__asm	mulpd	xmm2,[edi]\
		__asm	mulpd	xmm3,[edi]\
		\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm1,xmm3\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm3,xmm1\
		\
		/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm4,xmm5,xmm6,xmm7): swap xmm[01<->[23 */\
		__asm	addpd	xmm2,xmm4\
		__asm	addpd	xmm3,xmm5\
		__asm	movaps	[eax     ],xmm2\
		__asm	movaps	[eax+0x10],xmm3\
		__asm	addpd	xmm4,xmm4\
		__asm	addpd	xmm5,xmm5\
		__asm	subpd	xmm2,xmm4\
		__asm	subpd	xmm3,xmm5\
		__asm	movaps	[ecx     ],xmm2\
		__asm	movaps	[ecx+0x10],xmm3\
		\
		__asm	addpd	xmm0,xmm7\
		__asm	subpd	xmm1,xmm6\
		__asm	movaps	[ebx     ],xmm0\
		__asm	movaps	[ebx+0x10],xmm1\
		__asm	addpd	xmm7,xmm7\
		__asm	addpd	xmm6,xmm6\
		__asm	subpd	xmm0,xmm7\
		__asm	addpd	xmm1,xmm6\
		__asm	movaps	[edx     ],xmm0\
		__asm	movaps	[edx+0x10],xmm1\
		\
		/*...Block 7: r0C,r1C,r2C,r3C	*/\
		\
		__asm	add	eax ,0x80		/* r0C */\
		__asm	add	ebx,0x80\
		__asm	add	ecx,0x80\
		__asm	add	edx,0x80\
		\
		__asm	movaps	xmm4,[ecx     ]\
		__asm	movaps	xmm0,[edx     ]\
		__asm	movaps	xmm5,[ecx+0x10]\
		__asm	movaps	xmm1,[edx+0x10]\
		__asm	movaps	xmm6,[ecx     ]\
		__asm	movaps	xmm2,[edx     ]\
		__asm	movaps	xmm7,[ecx+0x10]\
		__asm	movaps	xmm3,[edx+0x10]\
		\
		__asm	mulpd	xmm4,[esi+0x10]\
		__asm	mulpd	xmm0,[esi     ]\
		__asm	mulpd	xmm5,[esi+0x10]\
		__asm	mulpd	xmm1,[esi     ]\
		__asm	mulpd	xmm6,[esi     ]\
		__asm	mulpd	xmm2,[esi+0x10]\
		__asm	mulpd	xmm7,[esi     ]\
		__asm	mulpd	xmm3,[esi+0x10]\
		__asm	subpd	xmm5,xmm6\
		__asm	subpd	xmm1,xmm2\
		__asm	addpd	xmm4,xmm7\
		__asm	addpd	xmm0,xmm3\
		__asm	movaps	xmm7,xmm5\
		__asm	movaps	xmm6,xmm4\
		\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm5,xmm1\
		__asm	subpd	xmm6,xmm0\
		__asm	subpd	xmm7,xmm1\
		\
		__asm	movaps	xmm2,[ebx     ]\
		__asm	movaps	xmm3,[ebx+0x10]\
		__asm	movaps	xmm0,[eax     ]\
		__asm	movaps	xmm1,[eax+0x10]\
		__asm	subpd	xmm2,[ebx+0x10]\
		__asm	addpd	xmm3,[ebx     ]\
		__asm	mulpd	xmm2,[edi]\
		__asm	mulpd	xmm3,[edi]\
		\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm1,xmm3\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm3,xmm1\
		\
		/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm0,xmm1,xmm2,xmm3,xmm6,xmm7,xmm4,xmm5): swap xmm[45<->[67 */\
		__asm	addpd	xmm0,xmm6\
		__asm	addpd	xmm1,xmm7\
		__asm	movaps	[eax     ],xmm0\
		__asm	movaps	[eax+0x10],xmm1\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm7,xmm7\
		__asm	subpd	xmm0,xmm6\
		__asm	subpd	xmm1,xmm7\
		__asm	movaps	[ecx     ],xmm0\
		__asm	movaps	[ecx+0x10],xmm1\
		\
		__asm	addpd	xmm2,xmm5\
		__asm	subpd	xmm3,xmm4\
		__asm	movaps	[ebx     ],xmm2\
		__asm	movaps	[ebx+0x10],xmm3\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm4,xmm4\
		__asm	subpd	xmm2,xmm5\
		__asm	addpd	xmm3,xmm4\
		__asm	movaps	[edx     ],xmm2\
		__asm	movaps	[edx+0x10],xmm3\
		\
		/*...Block 2: r02,r12,r22,r32	*/\
		\
		__asm	sub	eax ,0xa0		/* r02 */\
		__asm	sub	ebx,0xa0\
		__asm	sub	ecx,0xa0\
		__asm	sub	edx,0xa0\
		__asm	add	edi ,0x30		/* cc1 */\
		__asm	add	esi ,0x40		/* cc3 */\
		\
		__asm	movaps	xmm4,[ecx     ]\
		__asm	movaps	xmm0,[edx     ]\
		__asm	movaps	xmm5,[ecx+0x10]\
		__asm	movaps	xmm1,[edx+0x10]\
		__asm	movaps	xmm6,[ecx     ]\
		__asm	movaps	xmm2,[edx     ]\
		__asm	movaps	xmm7,[ecx+0x10]\
		__asm	movaps	xmm3,[edx+0x10]\
		\
		__asm	mulpd	xmm4,[edi     ]\
		__asm	mulpd	xmm0,[esi     ]\
		__asm	mulpd	xmm5,[edi     ]\
		__asm	mulpd	xmm1,[esi     ]\
		__asm	mulpd	xmm6,[edi+0x10]\
		__asm	mulpd	xmm2,[esi+0x10]\
		__asm	mulpd	xmm7,[edi+0x10]\
		__asm	mulpd	xmm3,[esi+0x10]\
		__asm	subpd	xmm5,xmm6\
		__asm	subpd	xmm1,xmm2\
		__asm	addpd	xmm4,xmm7\
		__asm	addpd	xmm0,xmm3\
		__asm	movaps	xmm7,xmm5\
		__asm	movaps	xmm6,xmm4\
		\
		__asm	subpd	xmm4,xmm0\
		__asm	subpd	xmm5,xmm1\
		__asm	addpd	xmm6,xmm0\
		__asm	addpd	xmm7,xmm1\
		\
		__asm	sub	esi ,0x40		/* cc0 */\
		__asm	movaps	xmm2,[ebx     ]\
		__asm	movaps	xmm3,[ebx+0x10]\
		__asm	movaps	xmm0,[ebx     ]\
		__asm	movaps	xmm1,[ebx+0x10]\
		\
		__asm	mulpd	xmm2,[esi     ]\
		__asm	mulpd	xmm1,[esi+0x10]\
		__asm	mulpd	xmm3,[esi     ]\
		__asm	mulpd	xmm0,[esi+0x10]\
		__asm	addpd	xmm2,xmm1\
		__asm	subpd	xmm3,xmm0\
		\
		__asm	movaps	xmm0,[eax     ]\
		__asm	movaps	xmm1,[eax+0x10]\
		\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm1,xmm3\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm3,xmm1\
		\
		/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm6,xmm7,xmm4,xmm5): swap xmm[01<->[23, xmm[45<->[67 */\
		__asm	addpd	xmm2,xmm6\
		__asm	addpd	xmm3,xmm7\
		__asm	movaps	[eax     ],xmm2\
		__asm	movaps	[eax+0x10],xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm7,xmm7\
		__asm	subpd	xmm2,xmm6\
		__asm	subpd	xmm3,xmm7\
		__asm	movaps	[ecx     ],xmm2\
		__asm	movaps	[ecx+0x10],xmm3\
		\
		__asm	addpd	xmm0,xmm5\
		__asm	subpd	xmm1,xmm4\
		__asm	movaps	[ebx     ],xmm0\
		__asm	movaps	[ebx+0x10],xmm1\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm4,xmm4\
		__asm	subpd	xmm0,xmm5\
		__asm	addpd	xmm1,xmm4\
		__asm	movaps	[edx     ],xmm0\
		__asm	movaps	[edx+0x10],xmm1\
		\
		/*...Block 6: r0A,r1A,r2A,r3A	*/\
		\
		__asm	add	eax ,0x80		/* r0A */\
		__asm	add	ebx,0x80\
		__asm	add	ecx,0x80\
		__asm	add	edx,0x80\
		__asm	add	esi ,0x40		/* cc3 */\
		\
		__asm	movaps	xmm4,[ecx     ]\
		__asm	movaps	xmm0,[edx     ]\
		__asm	movaps	xmm5,[ecx+0x10]\
		__asm	movaps	xmm1,[edx+0x10]\
		__asm	movaps	xmm6,[ecx     ]\
		__asm	movaps	xmm2,[edx     ]\
		__asm	movaps	xmm7,[ecx+0x10]\
		__asm	movaps	xmm3,[edx+0x10]\
		\
		__asm	mulpd	xmm4,[esi+0x10]\
		__asm	mulpd	xmm0,[edi     ]\
		__asm	mulpd	xmm5,[esi+0x10]\
		__asm	mulpd	xmm1,[edi     ]\
		__asm	mulpd	xmm6,[esi     ]\
		__asm	mulpd	xmm2,[edi+0x10]\
		__asm	mulpd	xmm7,[esi     ]\
		__asm	mulpd	xmm3,[edi+0x10]\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm1,xmm2\
		__asm	addpd	xmm4,xmm7\
		__asm	subpd	xmm0,xmm3\
		__asm	movaps	xmm7,xmm5\
		__asm	movaps	xmm6,xmm4\
		\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm5,xmm1\
		__asm	subpd	xmm6,xmm0\
		__asm	subpd	xmm7,xmm1\
		\
		__asm	sub	esi ,0x40		/* cc0 */\
		__asm	movaps	xmm2,[ebx     ]\
		__asm	movaps	xmm3,[ebx+0x10]\
		__asm	movaps	xmm0,[ebx     ]\
		__asm	movaps	xmm1,[ebx+0x10]\
		\
		__asm	mulpd	xmm2,[esi+0x10]\
		__asm	mulpd	xmm1,[esi     ]\
		__asm	mulpd	xmm3,[esi+0x10]\
		__asm	mulpd	xmm0,[esi     ]\
		__asm	subpd	xmm2,xmm1\
		__asm	addpd	xmm3,xmm0\
		\
		__asm	movaps	xmm0,[eax     ]\
		__asm	movaps	xmm1,[eax+0x10]\
		\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm1,xmm3\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm3,xmm1\
		\
		/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm0,xmm1,xmm2,xmm3,xmm6,xmm7,xmm4,xmm5): swap xmm[45<->[67 */\
		__asm	addpd	xmm0,xmm6\
		__asm	addpd	xmm1,xmm7\
		__asm	movaps	[eax     ],xmm0\
		__asm	movaps	[eax+0x10],xmm1\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm7,xmm7\
		__asm	subpd	xmm0,xmm6\
		__asm	subpd	xmm1,xmm7\
		__asm	movaps	[ecx     ],xmm0\
		__asm	movaps	[ecx+0x10],xmm1\
		\
		__asm	addpd	xmm2,xmm5\
		__asm	subpd	xmm3,xmm4\
		__asm	movaps	[ebx     ],xmm2\
		__asm	movaps	[ebx+0x10],xmm3\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm4,xmm4\
		__asm	subpd	xmm2,xmm5\
		__asm	addpd	xmm3,xmm4\
		__asm	movaps	[edx     ],xmm2\
		__asm	movaps	[edx+0x10],xmm3\
		\
		/*...Block 4: r06,r16,r26,r36	*/\
		\
		__asm	sub	eax ,0x40		/* r06 */\
		__asm	sub	ebx,0x40\
		__asm	sub	ecx,0x40\
		__asm	sub	edx,0x40\
		__asm	add	esi ,0x40		/* cc3 */\
		\
		__asm	movaps	xmm4,[ecx     ]\
		__asm	movaps	xmm0,[edx     ]\
		__asm	movaps	xmm5,[ecx+0x10]\
		__asm	movaps	xmm1,[edx+0x10]\
		__asm	movaps	xmm6,[ecx     ]\
		__asm	movaps	xmm2,[edx     ]\
		__asm	movaps	xmm7,[ecx+0x10]\
		__asm	movaps	xmm3,[edx+0x10]\
		\
		__asm	mulpd	xmm4,[esi     ]\
		__asm	mulpd	xmm0,[edi+0x10]\
		__asm	mulpd	xmm5,[esi     ]\
		__asm	mulpd	xmm1,[edi+0x10]\
		__asm	mulpd	xmm6,[esi+0x10]\
		__asm	mulpd	xmm2,[edi     ]\
		__asm	mulpd	xmm7,[esi+0x10]\
		__asm	mulpd	xmm3,[edi     ]\
		__asm	subpd	xmm5,xmm6\
		__asm	addpd	xmm1,xmm2\
		__asm	addpd	xmm4,xmm7\
		__asm	subpd	xmm0,xmm3\
		__asm	movaps	xmm7,xmm5\
		__asm	movaps	xmm6,xmm4\
		\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm5,xmm1\
		__asm	subpd	xmm6,xmm0\
		__asm	subpd	xmm7,xmm1\
		\
		__asm	sub	esi ,0x40		/* cc0 */\
		__asm	movaps	xmm2,[ebx     ]\
		__asm	movaps	xmm3,[ebx+0x10]\
		__asm	movaps	xmm0,[ebx     ]\
		__asm	movaps	xmm1,[ebx+0x10]\
		\
		__asm	mulpd	xmm2,[esi+0x10]\
		__asm	mulpd	xmm1,[esi     ]\
		__asm	mulpd	xmm3,[esi+0x10]\
		__asm	mulpd	xmm0,[esi     ]\
		__asm	addpd	xmm2,xmm1\
		__asm	subpd	xmm3,xmm0\
		\
		__asm	movaps	xmm0,[eax     ]\
		__asm	movaps	xmm1,[eax+0x10]\
		\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm1,xmm3\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm3,xmm1\
		\
		/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm2,xmm3,xmm0,xmm1,xmm6,xmm7,xmm4,xmm5): swap xmm[01<->[23, xmm[45<->[67 */\
		__asm	addpd	xmm2,xmm6\
		__asm	addpd	xmm3,xmm7\
		__asm	movaps	[eax     ],xmm2\
		__asm	movaps	[eax+0x10],xmm3\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm7,xmm7\
		__asm	subpd	xmm2,xmm6\
		__asm	subpd	xmm3,xmm7\
		__asm	movaps	[ecx     ],xmm2\
		__asm	movaps	[ecx+0x10],xmm3\
		\
		__asm	addpd	xmm0,xmm5\
		__asm	subpd	xmm1,xmm4\
		__asm	movaps	[ebx     ],xmm0\
		__asm	movaps	[ebx+0x10],xmm1\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm4,xmm4\
		__asm	subpd	xmm0,xmm5\
		__asm	addpd	xmm1,xmm4\
		__asm	movaps	[edx     ],xmm0\
		__asm	movaps	[edx+0x10],xmm1\
		\
		/*...Block 8: r0E,r1E,r2E,r3E	*/\
		\
		__asm	add	eax ,0x80		/* r0E */\
		__asm	add	ebx,0x80\
		__asm	add	ecx,0x80\
		__asm	add	edx,0x80\
		__asm	add	esi ,0x40		/* cc3 */\
		\
		__asm	movaps	xmm4,[ecx     ]\
		__asm	movaps	xmm0,[edx     ]\
		__asm	movaps	xmm5,[ecx+0x10]\
		__asm	movaps	xmm1,[edx+0x10]\
		__asm	movaps	xmm6,[ecx     ]\
		__asm	movaps	xmm2,[edx     ]\
		__asm	movaps	xmm7,[ecx+0x10]\
		__asm	movaps	xmm3,[edx+0x10]\
		\
		__asm	mulpd	xmm4,[edi+0x10]\
		__asm	mulpd	xmm0,[esi+0x10]\
		__asm	mulpd	xmm5,[edi+0x10]\
		__asm	mulpd	xmm1,[esi+0x10]\
		__asm	mulpd	xmm6,[edi     ]\
		__asm	mulpd	xmm2,[esi     ]\
		__asm	mulpd	xmm7,[edi     ]\
		__asm	mulpd	xmm3,[esi     ]\
		__asm	subpd	xmm5,xmm6\
		__asm	subpd	xmm1,xmm2\
		__asm	addpd	xmm4,xmm7\
		__asm	addpd	xmm0,xmm3\
		__asm	movaps	xmm7,xmm5\
		__asm	movaps	xmm6,xmm4\
		\
		__asm	addpd	xmm4,xmm0\
		__asm	addpd	xmm5,xmm1\
		__asm	subpd	xmm6,xmm0\
		__asm	subpd	xmm7,xmm1\
		\
		__asm	sub	esi ,0x40		/* cc0 */\
		__asm	movaps	xmm2,[ebx     ]\
		__asm	movaps	xmm3,[ebx+0x10]\
		__asm	movaps	xmm0,[ebx     ]\
		__asm	movaps	xmm1,[ebx+0x10]\
		\
		__asm	mulpd	xmm2,[esi     ]\
		__asm	mulpd	xmm1,[esi+0x10]\
		__asm	mulpd	xmm3,[esi     ]\
		__asm	mulpd	xmm0,[esi+0x10]\
		__asm	subpd	xmm2,xmm1\
		__asm	addpd	xmm3,xmm0\
		\
		__asm	movaps	xmm0,[eax     ]\
		__asm	movaps	xmm1,[eax+0x10]\
		\
		__asm	subpd	xmm0,xmm2\
		__asm	subpd	xmm1,xmm3\
		__asm	addpd	xmm2,xmm2\
		__asm	addpd	xmm3,xmm3\
		__asm	addpd	xmm2,xmm0\
		__asm	addpd	xmm3,xmm1\
		\
		/* SSE2_RADIX4_DIT_LASTPAIR_IN_PLACE_STRIDE(xmm0,xmm1,xmm2,xmm3,xmm6,xmm7,xmm4,xmm5): swap xmm[45<->[67 */\
		__asm	addpd	xmm0,xmm6\
		__asm	addpd	xmm1,xmm7\
		__asm	movaps	[eax     ],xmm0\
		__asm	movaps	[eax+0x10],xmm1\
		__asm	addpd	xmm6,xmm6\
		__asm	addpd	xmm7,xmm7\
		__asm	subpd	xmm0,xmm6\
		__asm	subpd	xmm1,xmm7\
		__asm	movaps	[ecx     ],xmm0\
		__asm	movaps	[ecx+0x10],xmm1\
		\
		__asm	addpd	xmm2,xmm5\
		__asm	subpd	xmm3,xmm4\
		__asm	movaps	[ebx     ],xmm2\
		__asm	movaps	[ebx+0x10],xmm3\
		__asm	addpd	xmm5,xmm5\
		__asm	addpd	xmm4,xmm4\
		__asm	subpd	xmm2,xmm5\
		__asm	addpd	xmm3,xmm4\
		__asm	movaps	[edx     ],xmm2\
		__asm	movaps	[edx+0x10],xmm3\
	}

#endif	/* radix32_ditN_cy_dif1_win_h_included */

