/*******************************************************************************
*                                                                             *
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

/*******************************************************************************
   We now include this header file if it was not included before.
*******************************************************************************/
#ifndef radix16_wrapper_square_gcc_h_included
#define radix16_wrapper_square_gcc_h_included

	#define SSE2_RADIX16_WRAPPER_DIF(Xadd0,Xadd1,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15)\
	{\
	__asm__ volatile (\
	"/*************************************************************/\n\t"\
	"/*                  1st set of inputs:                       */\n\t"\
	"/*************************************************************/\n\t"\
		"/*...Block 1: */									\n\t	/*...Block 2: */					\n\t"\
		"movq		%[__add0],%%rax							\n\t	movq		%[__r9] ,%%r10			\n\t"\
		"movq		%[__add1],%%rbx							\n\t	movq		%[__c2] ,%%r11			\n\t"\
		"movq		%[__r1] ,%%rcx							\n\t	movaps		0x20(%%rax),%%xmm14		\n\t"\
		"movq		%[__c4] ,%%rdx							\n\t	movaps		0x20(%%rax),%%xmm8		\n\t"\
		"movaps		0x40(%%rax),%%xmm6						\n\t	unpckhpd	0x20(%%rbx),%%xmm14		\n\t"\
		"movaps		0x40(%%rbx),%%xmm2						\n\t	unpcklpd	0x20(%%rbx),%%xmm8		\n\t"\
		"movaps		%%xmm6,%%xmm0							\n\t	movaps			%%xmm14,0x100(%%r10)\n\t"\
		"movaps		%%xmm2,%%xmm3							\n\t	movaps		0x30(%%rax),%%xmm15		\n\t"\
		"unpckhpd	%%xmm2,%%xmm6							\n\t	movaps		0x30(%%rax),%%xmm9		\n\t"\
		"unpcklpd	%%xmm3,%%xmm0							\n\t	unpckhpd	0x30(%%rbx),%%xmm15		\n\t"\
		"movaps		%%xmm6,0x140(%%rcx)						\n\t	unpcklpd	0x30(%%rbx),%%xmm9		\n\t"\
		"movaps		0x50(%%rax),%%xmm7						\n\t	movaps		%%xmm15,0x110(%%r10)	\n\t"\
		"movaps		0x50(%%rbx),%%xmm4						\n\t	movaps		%%xmm8,%%xmm10			\n\t"\
		"movaps		%%xmm7,%%xmm1							\n\t	movaps		%%xmm9,%%xmm11			\n\t"\
		"movaps		%%xmm4,%%xmm5							\n\t	mulpd		    (%%r11),%%xmm8		\n\t"\
		"unpckhpd	%%xmm4,%%xmm7							\n\t	mulpd		    (%%r11),%%xmm9		\n\t"\
		"unpcklpd	%%xmm5,%%xmm1							\n\t	mulpd		0x10(%%r11),%%xmm10		\n\t"\
		"movaps		%%xmm7,0x150(%%rcx)						\n\t	mulpd		0x10(%%r11),%%xmm11		\n\t"\
		"movaps		%%xmm0,%%xmm2							\n\t	addpd		%%xmm10,%%xmm9			\n\t"\
		"movaps		%%xmm1,%%xmm3							\n\t	subpd		%%xmm11,%%xmm8			\n\t"\
		"/* Rest identical to code in radix16_dif_pass: */	\n\t	movaps		%%xmm9,%%xmm11			\n\t"\
		"mulpd		    (%%rdx),%%xmm0						\n\t	movaps		%%xmm8,%%xmm10			\n\t"\
		"mulpd		    (%%rdx),%%xmm1						\n\t	movq		%[__c10],%%r11			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm2						\n\t	movaps		0xa0(%%rax),%%xmm14		\n\t"\
		"mulpd		0x10(%%rdx),%%xmm3						\n\t	movaps		0xa0(%%rax),%%xmm12		\n\t"\
		"addpd		%%xmm2,%%xmm1							\n\t	unpckhpd	0xa0(%%rbx),%%xmm14		\n\t"\
		"subpd		%%xmm3,%%xmm0							\n\t	unpcklpd	0xa0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm1,%%xmm3							\n\t	movaps			%%xmm14,0x120(%%r10)\n\t"\
		"movaps		%%xmm0,%%xmm2							\n\t	movaps		0xb0(%%rax),%%xmm15		\n\t"\
		"movq		%[__c12],%%rdx							\n\t	movaps		0xb0(%%rax),%%xmm13		\n\t"\
		"movaps		0xc0(%%rax),%%xmm6						\n\t	unpckhpd	0xb0(%%rbx),%%xmm15		\n\t"\
		"movaps		0xc0(%%rax),%%xmm4						\n\t	unpcklpd	0xb0(%%rbx),%%xmm13		\n\t"\
		"unpckhpd	0xc0(%%rbx),%%xmm6						\n\t	movaps			%%xmm15,0x130(%%r10)\n\t"\
		"unpcklpd	0xc0(%%rbx),%%xmm4						\n\t	movaps			%%xmm12,%%xmm14		\n\t"\
		"movaps		%%xmm6,0x160(%%rcx)						\n\t	movaps			%%xmm13,%%xmm15		\n\t"\
		"movaps		0xd0(%%rax),%%xmm7						\n\t	mulpd		    (%%r11),%%xmm12		\n\t"\
		"movaps		0xd0(%%rax),%%xmm5						\n\t	mulpd		    (%%r11),%%xmm13		\n\t"\
		"unpckhpd	0xd0(%%rbx),%%xmm7						\n\t	mulpd		0x10(%%r11),%%xmm14		\n\t"\
		"unpcklpd	0xd0(%%rbx),%%xmm5						\n\t	mulpd		0x10(%%r11),%%xmm15		\n\t"\
		"movaps		%%xmm7,0x170(%%rcx)						\n\t	addpd		%%xmm14,%%xmm13			\n\t"\
		"movaps		%%xmm4,%%xmm6							\n\t	subpd		%%xmm15,%%xmm12			\n\t"\
		"movaps		%%xmm5,%%xmm7							\n\t	addpd		%%xmm12,%%xmm8			\n\t"\
		"mulpd		    (%%rdx),%%xmm4						\n\t	addpd		%%xmm13,%%xmm9			\n\t"\
		"mulpd		    (%%rdx),%%xmm5						\n\t	subpd		%%xmm12,%%xmm10			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6						\n\t	subpd		%%xmm13,%%xmm11			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7						\n\t	movq		%[__c14],%%r11			\n\t"\
		"addpd		%%xmm6,%%xmm5							\n\t	movaps		0xe0(%%rax),%%xmm14		\n\t"\
		"subpd		%%xmm7,%%xmm4							\n\t	movaps		0xe0(%%rax),%%xmm12		\n\t"\
		"addpd		%%xmm4,%%xmm0							\n\t	unpckhpd	0xe0(%%rbx),%%xmm14		\n\t"\
		"addpd		%%xmm5,%%xmm1							\n\t	unpcklpd	0xe0(%%rbx),%%xmm12		\n\t"\
		"subpd		%%xmm4,%%xmm2							\n\t	movaps			%%xmm14,0x160(%%r10)\n\t"\
		"subpd		%%xmm5,%%xmm3							\n\t	movaps		0xf0(%%rax),%%xmm15		\n\t"\
		"movq		%[__c8] ,%%rdx							\n\t	movaps		0xf0(%%rax),%%xmm13		\n\t"\
		"movaps		0x80(%%rax),%%xmm6						\n\t	unpckhpd	0xf0(%%rbx),%%xmm15		\n\t"\
		"movaps		0x80(%%rax),%%xmm4						\n\t	unpcklpd	0xf0(%%rbx),%%xmm13		\n\t"\
		"unpckhpd	0x80(%%rbx),%%xmm6						\n\t	movaps			%%xmm15,0x170(%%r10)\n\t"\
		"unpcklpd	0x80(%%rbx),%%xmm4						\n\t	movaps			%%xmm12,%%xmm14		\n\t"\
		"movaps			%%xmm6,0x120(%%rcx)					\n\t	movaps			%%xmm13,%%xmm15		\n\t"\
		"movaps		0x90(%%rax),%%xmm7						\n\t	mulpd		    (%%r11),%%xmm12		\n\t"\
		"movaps		0x90(%%rax),%%xmm5						\n\t	mulpd		    (%%r11),%%xmm13		\n\t"\
		"unpckhpd	0x90(%%rbx),%%xmm7						\n\t	mulpd		0x10(%%r11),%%xmm14		\n\t"\
		"unpcklpd	0x90(%%rbx),%%xmm5						\n\t	mulpd		0x10(%%r11),%%xmm15		\n\t"\
		"movaps		%%xmm7,0x130(%%rcx)						\n\t	addpd		%%xmm14,%%xmm13			\n\t"\
		"movaps		%%xmm4,%%xmm6							\n\t	subpd		%%xmm15,%%xmm12			\n\t"\
		"movaps		%%xmm5,%%xmm7							\n\t	movaps		%%xmm13,0x070(%%r10)	\n\t"\
		"mulpd		    (%%rdx),%%xmm4						\n\t	movaps		%%xmm12,0x060(%%r10)	\n\t"\
		"mulpd		    (%%rdx),%%xmm5						\n\t	movq		%[__c6] ,%%r11			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6						\n\t	movaps		0x60(%%rax),%%xmm14		\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7						\n\t	movaps		0x60(%%rax),%%xmm12		\n\t"\
		"addpd		%%xmm6,%%xmm5							\n\t	unpckhpd	0x60(%%rbx),%%xmm14		\n\t"\
		"subpd		%%xmm7,%%xmm4							\n\t	unpcklpd	0x60(%%rbx),%%xmm12		\n\t"\
		"movaps		    (%%rax),%%xmm6						\n\t	movaps			%%xmm14,0x140(%%r10)\n\t"\
		"movaps		    (%%rax),%%xmm7						\n\t	movaps		0x70(%%rax),%%xmm15		\n\t"\
		"unpckhpd	    (%%rbx),%%xmm6						\n\t	movaps		0x70(%%rax),%%xmm13		\n\t"\
		"unpcklpd	    (%%rbx),%%xmm7						\n\t	unpckhpd	0x70(%%rbx),%%xmm15		\n\t"\
		"movaps		%%xmm6,0x100(%%rcx)						\n\t	unpcklpd	0x70(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7,     (%%rcx)						\n\t	movaps			%%xmm15,0x150(%%r10)\n\t"\
		"movaps		0x10(%%rax),%%xmm6						\n\t	movaps			%%xmm12,%%xmm14		\n\t"\
		"movaps		0x10(%%rax),%%xmm7						\n\t	movaps			%%xmm13,%%xmm15		\n\t"\
		"unpckhpd	0x10(%%rbx),%%xmm6						\n\t	mulpd		    (%%r11),%%xmm12		\n\t"\
		"unpcklpd	0x10(%%rbx),%%xmm7						\n\t	mulpd		    (%%r11),%%xmm13		\n\t"\
		"movaps		%%xmm6,0x110(%%rcx)						\n\t	mulpd		0x10(%%r11),%%xmm14		\n\t"\
		"movaps		    (%%rcx),%%xmm6						\n\t	mulpd		0x10(%%r11),%%xmm15		\n\t"\
		"subpd		%%xmm4,%%xmm6							\n\t	addpd		%%xmm14,%%xmm13			\n\t"\
		"subpd		%%xmm5,%%xmm7							\n\t	subpd		%%xmm15,%%xmm12			\n\t"\
		"addpd		%%xmm4,%%xmm4							\n\t	movaps		%%xmm13,%%xmm15			\n\t"\
		"addpd		%%xmm5,%%xmm5							\n\t	movaps		%%xmm12,%%xmm14			\n\t"\
		"addpd		%%xmm6,%%xmm4							\n\t	subpd		0x060(%%r10),%%xmm12	\n\t"\
		"addpd		%%xmm7,%%xmm5							\n\t	subpd		0x070(%%r10),%%xmm13	\n\t"\
		"/* Finish radix-4 butterfly: */					\n\t	addpd		0x060(%%r10),%%xmm14	\n\t"\
		"subpd		%%xmm0,%%xmm4							\n\t	addpd		0x070(%%r10),%%xmm15	\n\t"\
		"subpd		%%xmm1,%%xmm5							\n\t	/* Finish radix-4 butterfly: */		\n\t"\
		"movaps		%%xmm4,0x040(%%rcx)						\n\t	subpd		%%xmm14,%%xmm8			\n\t"\
		"movaps		%%xmm5,0x050(%%rcx)						\n\t	subpd		%%xmm15,%%xmm9			\n\t"\
		"addpd		%%xmm0,%%xmm0							\n\t	movaps		%%xmm8,0x040(%%r10)		\n\t"\
		"addpd		%%xmm1,%%xmm1							\n\t	movaps		%%xmm9,0x050(%%r10)		\n\t"\
		"addpd		%%xmm4,%%xmm0							\n\t	addpd		%%xmm14,%%xmm14			\n\t"\
		"addpd		%%xmm5,%%xmm1							\n\t	addpd		%%xmm15,%%xmm15			\n\t"\
		"movaps		%%xmm0,     (%%rcx)						\n\t	addpd		%%xmm8,%%xmm14			\n\t"\
		"movaps		%%xmm1,0x010(%%rcx)						\n\t	addpd		%%xmm9,%%xmm15			\n\t"\
		"subpd		%%xmm3,%%xmm6							\n\t	movaps		%%xmm14,     (%%r10)	\n\t"\
		"subpd		%%xmm2,%%xmm7							\n\t	movaps		%%xmm15,0x010(%%r10)	\n\t"\
		"movaps		%%xmm6,0x020(%%rcx)						\n\t	subpd		%%xmm13,%%xmm10			\n\t"\
		"movaps		%%xmm7,0x070(%%rcx)						\n\t	subpd		%%xmm12,%%xmm11			\n\t"\
		"addpd		%%xmm3,%%xmm3							\n\t	movaps		%%xmm10,0x020(%%r10)	\n\t"\
		"addpd		%%xmm2,%%xmm2							\n\t	movaps		%%xmm11,0x070(%%r10)	\n\t"\
		"addpd		%%xmm6,%%xmm3							\n\t	addpd		%%xmm13,%%xmm13			\n\t"\
		"addpd		%%xmm7,%%xmm2							\n\t	addpd		%%xmm12,%%xmm12			\n\t"\
		"movaps		%%xmm3,0x060(%%rcx)						\n\t	addpd		%%xmm10,%%xmm13			\n\t"\
		"movaps		%%xmm2,0x030(%%rcx)						\n\t	addpd		%%xmm11,%%xmm12			\n\t"\
		"													\n\t	movaps		%%xmm13,0x060(%%r10)	\n\t"\
		"													\n\t	movaps		%%xmm12,0x030(%%r10)	\n\t"\
	"/****************************************************************************************************/\n\t"\
	"/* Each of the above 2 radix-4 blocks is unique, but the ensuing 2 blocks [operating on odd-indexed */\n\t"\
	"/* elements from the unpck*pd commands which were stored to temporaries] can use a common macro:    */\n\t"\
	"/****************************************************************************************************/\n\t"\
		"/*...Block 3: */											\n\t"\
		"/*	SSE2_RADIX4_DIF_4TWIDDLE(r17,r21,r19,r23, r17,c1): */	\n\t"\
		"/* Do the p0,p8 combo: */									\n\t"\
		"movq		%[__r17],%%rax									\n\t"\
		"movq		%[__c1] ,%%rbx									\n\t"\
		"movq		%%rax   ,%%rcx									\n\t"\
		"addq		$0x20   ,%%rcx									\n\t"\
		"movaps		    (%%rax),%%xmm0								\n\t"\
		"movaps	        (%%rcx),%%xmm4								\n\t"\
		"movaps		0x10(%%rax),%%xmm1								\n\t"\
		"movaps		0x10(%%rcx),%%xmm5								\n\t"\
		"movaps		    (%%rbx),%%xmm6								\n\t"\
		"movaps		0x10(%%rbx),%%xmm7								\n\t"\
		"movaps		%%xmm0,%%xmm2									\n\t"\
		"movaps		%%xmm1,%%xmm3									\n\t"\
		"mulpd		%%xmm6,%%xmm0									\n\t"\
		"mulpd		%%xmm6,%%xmm1									\n\t"\
		"mulpd		%%xmm7,%%xmm2									\n\t"\
		"mulpd		%%xmm7,%%xmm3									\n\t"\
		"movaps		%%xmm4,%%xmm6									\n\t"\
		"addpd		%%xmm2,%%xmm1									\n\t"\
		"movaps		%%xmm5,%%xmm7									\n\t"\
		"mulpd		0x20(%%rbx),%%xmm4								\n\t"\
		"subpd		%%xmm3,%%xmm0									\n\t"\
		"mulpd		0x20(%%rbx),%%xmm5								\n\t"\
		"mulpd		0x30(%%rbx),%%xmm6								\n\t"\
		"movaps		%%xmm0,%%xmm2									\n\t"\
		"mulpd		0x30(%%rbx),%%xmm7								\n\t"\
		"addpd		%%xmm6,%%xmm5									\n\t"\
		"movaps		%%xmm1,%%xmm3									\n\t"\
		"subpd		%%xmm7,%%xmm4									\n\t"\
		"addq		$0x40,%%rcx										\n\t"\
		"addq		$0x60,%%rbx										\n\t"\
		"movaps			(%%rcx),%%xmm6								\n\t"\
		"movaps		0x10(%%rcx),%%xmm7								\n\t"\
		"addpd		%%xmm4,%%xmm0									\n\t"\
		"addpd		%%xmm5,%%xmm1									\n\t"\
		"subpd		%%xmm4,%%xmm2									\n\t"\
		"subpd		%%xmm5,%%xmm3									\n\t"\
		"/* Do the p4,12 combo: */									\n\t"\
		"movaps		%%xmm6,%%xmm4									\n\t"\
		"movaps		%%xmm7,%%xmm5									\n\t"\
		"mulpd		    (%%rbx),%%xmm4								\n\t"\
		"mulpd		    (%%rbx),%%xmm5								\n\t"\
		"mulpd		0x10(%%rbx),%%xmm6								\n\t"\
		"mulpd		0x10(%%rbx),%%xmm7								\n\t"\
		"movq		%%rax,%%rdx										\n\t"\
		"addpd		%%xmm6,%%xmm5									\n\t"\
		"subpd		%%xmm7,%%xmm4									\n\t"\
		"movaps		%%xmm5,0x010(%%rdx)								\n\t"\
		"movaps		%%xmm4,     (%%rdx)								\n\t"\
		"addq		$0x40 ,%%rax									\n\t"\
		"subq		$0x20 ,%%rbx									\n\t"\
		"movaps		    (%%rax),%%xmm4								\n\t"\
		"movaps		0x10(%%rax),%%xmm5								\n\t"\
		"movaps			%%xmm4,%%xmm6								\n\t"\
		"movaps			%%xmm5,%%xmm7								\n\t"\
		"mulpd		    (%%rbx),%%xmm4								\n\t"\
		"mulpd		    (%%rbx),%%xmm5								\n\t"\
		"mulpd		0x10(%%rbx),%%xmm6								\n\t"\
		"mulpd		0x10(%%rbx),%%xmm7								\n\t"\
		"addpd		%%xmm6,%%xmm5									\n\t"\
		"subpd		%%xmm7,%%xmm4									\n\t"\
		"movaps		%%xmm5,%%xmm7									\n\t"\
		"movaps		%%xmm4,%%xmm6									\n\t"\
		"subpd		     (%%rdx),%%xmm4								\n\t"\
		"subpd		0x010(%%rdx),%%xmm5								\n\t"\
		"addpd		     (%%rdx),%%xmm6								\n\t"\
		"addpd		0x010(%%rdx),%%xmm7								\n\t"\
		"/* Finish radix-4 butterfly: */							\n\t"\
		"subpd		%%xmm6,%%xmm0									\n\t"\
		"subpd		%%xmm5,%%xmm2									\n\t"\
		"subpd		%%xmm7,%%xmm1									\n\t"\
		"subpd		%%xmm4,%%xmm3									\n\t"\
		"movaps		%%xmm0,0x040(%%rdx)								\n\t"\
		"movaps		%%xmm2,0x020(%%rdx)								\n\t"\
		"movaps		%%xmm1,0x050(%%rdx)								\n\t"\
		"movaps		%%xmm3,0x070(%%rdx)								\n\t"\
		"addpd		%%xmm6,%%xmm6									\n\t"\
		"addpd		%%xmm5,%%xmm5									\n\t"\
		"addpd		%%xmm7,%%xmm7									\n\t"\
		"addpd		%%xmm4,%%xmm4									\n\t"\
		"addpd		%%xmm0,%%xmm6									\n\t"\
		"addpd		%%xmm2,%%xmm5									\n\t"\
		"addpd		%%xmm1,%%xmm7									\n\t"\
		"addpd		%%xmm3,%%xmm4									\n\t"\
		"movaps		%%xmm6,     (%%rdx)								\n\t"\
		"movaps		%%xmm5,0x060(%%rdx)								\n\t"\
		"movaps		%%xmm7,0x010(%%rdx)								\n\t"\
		"movaps		%%xmm4,0x030(%%rdx)								\n\t"\
		"/*...Block 4: */											\n\t"\
		"/* SSE2_RADIX4_DIF_4TWIDDLE(r25,r29,r27,r31, r25,c3): */	\n\t"\
		"/* Do the p0,p8 combo: */									\n\t"\
		"movq		%[__r25],%%rax									\n\t"\
		"movq		%[__c3] ,%%rbx									\n\t"\
		"movq		%%rax   ,%%rcx									\n\t"\
		"addq		$0x20   ,%%rcx									\n\t"\
		"movaps		    (%%rax),%%xmm0								\n\t"\
		"movaps	        (%%rcx),%%xmm4								\n\t"\
		"movaps		0x10(%%rax),%%xmm1								\n\t"\
		"movaps		0x10(%%rcx),%%xmm5								\n\t"\
		"movaps		    (%%rbx),%%xmm6								\n\t"\
		"movaps		0x10(%%rbx),%%xmm7								\n\t"\
		"movaps		%%xmm0,%%xmm2									\n\t"\
		"movaps		%%xmm1,%%xmm3									\n\t"\
		"mulpd   	%%xmm6,%%xmm0									\n\t"\
		"mulpd   	%%xmm6,%%xmm1									\n\t"\
		"mulpd   	%%xmm7,%%xmm2									\n\t"\
		"mulpd   	%%xmm7,%%xmm3									\n\t"\
		"movaps		%%xmm4,%%xmm6									\n\t"\
		"addpd   	%%xmm2,%%xmm1									\n\t"\
		"movaps		%%xmm5,%%xmm7									\n\t"\
		"mulpd		0x20(%%rbx),%%xmm4								\n\t"\
		"subpd		%%xmm3,%%xmm0									\n\t"\
		"mulpd		0x20(%%rbx),%%xmm5								\n\t"\
		"mulpd		0x30(%%rbx),%%xmm6								\n\t"\
		"movaps		%%xmm0,%%xmm2									\n\t"\
		"mulpd		0x30(%%rbx),%%xmm7								\n\t"\
		"addpd		%%xmm6			  ,%%xmm5						\n\t"\
		"movaps		%%xmm1,%%xmm3									\n\t"\
		"subpd		%%xmm7,%%xmm4									\n\t"\
		"addq		$0x40	,%%rcx									\n\t"\
		"addq		$0x60	,%%rbx									\n\t"\
		"movaps		    (%%rcx),%%xmm6								\n\t"\
		"movaps		0x10(%%rcx),%%xmm7								\n\t"\
		"addpd		%%xmm4,%%xmm0									\n\t"\
		"addpd		%%xmm5,%%xmm1									\n\t"\
		"subpd		%%xmm4,%%xmm2									\n\t"\
		"subpd		%%xmm5,%%xmm3									\n\t"\
		"/* Do the p4,12 combo: */									\n\t"\
		"movaps		%%xmm6,%%xmm4									\n\t"\
		"movaps		%%xmm7,%%xmm5									\n\t"\
		"mulpd		    (%%rbx),%%xmm4								\n\t"\
		"mulpd		    (%%rbx),%%xmm5								\n\t"\
		"mulpd		0x10(%%rbx),%%xmm6								\n\t"\
		"mulpd		0x10(%%rbx),%%xmm7								\n\t"\
		"movq		%%rax,%%rdx										\n\t"\
		"addpd		%%xmm6,%%xmm5									\n\t"\
		"subpd		%%xmm7,%%xmm4									\n\t"\
		"movaps		%%xmm5,0x010(%%rdx)								\n\t"\
		"movaps		%%xmm4,     (%%rdx)								\n\t"\
		"addq		$0x40 ,%%rax									\n\t"\
		"subq		$0x20 ,%%rbx									\n\t"\
		"movaps		    (%%rax),%%xmm4								\n\t"\
		"movaps		0x10(%%rax),%%xmm5								\n\t"\
		"movaps			%%xmm4,%%xmm6								\n\t"\
		"movaps			%%xmm5,%%xmm7								\n\t"\
		"mulpd		    (%%rbx),%%xmm4								\n\t"\
		"mulpd		    (%%rbx),%%xmm5								\n\t"\
		"mulpd		0x10(%%rbx),%%xmm6								\n\t"\
		"mulpd		0x10(%%rbx),%%xmm7								\n\t"\
		"addpd		%%xmm6,%%xmm5									\n\t"\
		"subpd		%%xmm7,%%xmm4									\n\t"\
		"movaps		%%xmm5,%%xmm7									\n\t"\
		"movaps		%%xmm4,%%xmm6									\n\t"\
		"subpd		     (%%rdx),%%xmm4								\n\t"\
		"subpd		0x010(%%rdx),%%xmm5								\n\t"\
		"addpd		     (%%rdx),%%xmm6								\n\t"\
		"addpd		0x010(%%rdx),%%xmm7								\n\t"\
		"/* Finish radix-4 butterfly: */							\n\t"\
		"subpd		%%xmm6,%%xmm0									\n\t"\
		"subpd		%%xmm5,%%xmm2									\n\t"\
		"subpd		%%xmm7,%%xmm1									\n\t"\
		"subpd		%%xmm4,%%xmm3									\n\t"\
		"movaps		%%xmm0,0x040(%%rdx)								\n\t"\
		"movaps		%%xmm2,0x020(%%rdx)								\n\t"\
		"movaps		%%xmm1,0x050(%%rdx)								\n\t"\
		"movaps		%%xmm3,0x070(%%rdx)								\n\t"\
		"addpd		%%xmm6,%%xmm6									\n\t"\
		"addpd		%%xmm5,%%xmm5									\n\t"\
		"addpd		%%xmm7,%%xmm7									\n\t"\
		"addpd		%%xmm4,%%xmm4									\n\t"\
		"addpd		%%xmm0,%%xmm6									\n\t"\
		"addpd		%%xmm2,%%xmm5									\n\t"\
		"addpd		%%xmm1,%%xmm7									\n\t"\
		"addpd		%%xmm3,%%xmm4									\n\t"\
		"movaps		%%xmm6,     (%%rdx)								\n\t"\
		"movaps		%%xmm5,0x060(%%rdx)								\n\t"\
		"movaps		%%xmm7,0x010(%%rdx)								\n\t"\
		"movaps		%%xmm4,0x030(%%rdx)								\n\t"\
	"/**************************************************************************************/\n\t"\
	"/*...and now do four more radix-4 transforms, including the internal twiddle factors: */\n\t"\
	"/**************************************************************************************/\n\t"\
		"/*...Block 1: t1,9,17,25 */		\n\t"\
		"movq		%[__r1] ,%%rax			\n\t"\
		"movq		%[__r9] ,%%rbx			\n\t"\
		"movq		%[__r17],%%rcx			\n\t"\
		"movq		%[__r25],%%rdx			\n\t"\
		"movaps		     (%%rax),%%xmm0		\n\t"\
		"movaps		0x010(%%rax),%%xmm1		\n\t"\
		"movaps		     (%%rbx),%%xmm2		\n\t"\
		"movaps		0x010(%%rbx),%%xmm3		\n\t"\
		"subpd		     (%%rbx),%%xmm0		\n\t"\
		"subpd		0x010(%%rbx),%%xmm1		\n\t"\
		"addpd		     (%%rax),%%xmm2		\n\t"\
		"addpd		0x010(%%rax),%%xmm3		\n\t"\
		"movaps		     (%%rcx),%%xmm4		\n\t"\
		"movaps		0x010(%%rcx),%%xmm5		\n\t"\
		"movaps		     (%%rdx),%%xmm6		\n\t"\
		"movaps		0x010(%%rdx),%%xmm7		\n\t"\
		"subpd		     (%%rdx),%%xmm4		\n\t"\
		"subpd		0x010(%%rdx),%%xmm5		\n\t"\
		"addpd		     (%%rcx),%%xmm6		\n\t"\
		"addpd		0x010(%%rcx),%%xmm7		\n\t"\
		"subpd		%%xmm6,%%xmm2			\n\t"\
		"subpd		%%xmm7,%%xmm3			\n\t"\
		"addpd		%%xmm6,%%xmm6			\n\t"\
		"addpd		%%xmm7,%%xmm7			\n\t"\
		"movaps		%%xmm2,    (%%rcx)		\n\t"\
		"movaps		%%xmm3,0x10(%%rcx)		\n\t"\
		"addpd		%%xmm2,%%xmm6			\n\t"\
		"addpd		%%xmm3,%%xmm7			\n\t"\
		"movaps		%%xmm6,    (%%rax)		\n\t"\
		"movaps		%%xmm7,0x10(%%rax)		\n\t"\
		"subpd		%%xmm5,%%xmm0			\n\t"\
		"subpd		%%xmm4,%%xmm1			\n\t"\
		"addpd		%%xmm5,%%xmm5			\n\t"\
		"addpd		%%xmm4,%%xmm4			\n\t"\
		"movaps		%%xmm0,    (%%rbx)		\n\t"\
		"movaps		%%xmm1,0x10(%%rdx)		\n\t"\
		"addpd		%%xmm0,%%xmm5			\n\t"\
		"addpd		%%xmm1,%%xmm4			\n\t"\
		"movaps		%%xmm5,    (%%rdx)		\n\t"\
		"movaps		%%xmm4,0x10(%%rbx)		\n\t"\
		"/*...Block 3: t5,13,21,29 */		\n\t"\
		"addq		$0x40,%%rax				\n\t"\
		"addq		$0x40,%%rbx				\n\t"\
		"addq		$0x40,%%rcx				\n\t"\
		"addq		$0x40,%%rdx				\n\t"\
		"movaps		     (%%rax),%%xmm0		\n\t"\
		"movaps		0x010(%%rax),%%xmm1		\n\t"\
		"movaps		0x080(%%rax),%%xmm2		\n\t"\
		"movaps		0x090(%%rax),%%xmm3		\n\t"\
		"subpd		0x090(%%rax),%%xmm0		\n\t"\
		"subpd		0x080(%%rax),%%xmm1		\n\t"\
		"addpd		0x010(%%rax),%%xmm2		\n\t"\
		"addpd		     (%%rax),%%xmm3		\n\t"\
		"movq		%[__isrt2],%%rsi		\n\t"\
		"movaps		0x100(%%rax),%%xmm4		\n\t"\
		"movaps		0x110(%%rax),%%xmm5		\n\t"\
		"movaps		0x180(%%rax),%%xmm6		\n\t"\
		"movaps		0x190(%%rax),%%xmm7		\n\t"\
		"subpd		0x110(%%rax),%%xmm4		\n\t"\
		"addpd		0x100(%%rax),%%xmm5		\n\t"\
		"mulpd		(%%rsi),%%xmm4			\n\t"\
		"mulpd		(%%rsi),%%xmm5			\n\t"\
		"addpd		0x190(%%rax),%%xmm6		\n\t"\
		"subpd		0x180(%%rax),%%xmm7		\n\t"\
		"mulpd		(%%rsi),%%xmm6			\n\t"\
		"mulpd		(%%rsi),%%xmm7			\n\t"\
		"subpd		%%xmm6,%%xmm4			\n\t"\
		"subpd		%%xmm7,%%xmm5			\n\t"\
		"addpd		%%xmm6,%%xmm6			\n\t"\
		"addpd		%%xmm7,%%xmm7			\n\t"\
		"addpd		%%xmm4,%%xmm6			\n\t"\
		"addpd		%%xmm5,%%xmm7			\n\t"\
		"subpd		%%xmm4,%%xmm0			\n\t"\
		"subpd		%%xmm5,%%xmm2			\n\t"\
		"addpd		%%xmm4,%%xmm4			\n\t"\
		"addpd		%%xmm5,%%xmm5			\n\t"\
		"movaps		%%xmm0,    (%%rcx)		\n\t"\
		"movaps		%%xmm2,0x10(%%rcx)		\n\t"\
		"addpd		%%xmm0,%%xmm4			\n\t"\
		"addpd		%%xmm2,%%xmm5			\n\t"\
		"movaps		%%xmm4,    (%%rax)		\n\t"\
		"movaps		%%xmm5,0x10(%%rax)		\n\t"\
		"subpd		%%xmm7,%%xmm3			\n\t"\
		"subpd		%%xmm6,%%xmm1			\n\t"\
		"addpd		%%xmm7,%%xmm7			\n\t"\
		"addpd		%%xmm6,%%xmm6			\n\t"\
		"movaps		%%xmm3,    (%%rbx)		\n\t"\
		"movaps		%%xmm1,0x10(%%rdx)		\n\t"\
		"addpd		%%xmm3,%%xmm7			\n\t"\
		"addpd		%%xmm1,%%xmm6			\n\t"\
		"movaps		%%xmm7,    (%%rdx)		\n\t"\
		"movaps		%%xmm6,0x10(%%rbx)		\n\t"\
		"/*...Block 2: t3,11,19,27 */		\n\t"\
		"subq		$0x20,%%rax				\n\t"\
		"subq		$0x20,%%rbx				\n\t"\
		"subq		$0x20,%%rcx				\n\t"\
		"subq		$0x20,%%rdx				\n\t"\
		"movq		%[__cc0],%%rdi			\n\t"\
		"movaps		0x100(%%rax),%%xmm4		\n\t"\
		"movaps		0x180(%%rax),%%xmm6		\n\t"\
		"movaps		0x110(%%rax),%%xmm5		\n\t"\
		"movaps		0x190(%%rax),%%xmm7		\n\t"\
		"movaps		0x100(%%rax),%%xmm0		\n\t"\
		"movaps		0x180(%%rax),%%xmm2		\n\t"\
		"movaps		0x110(%%rax),%%xmm1		\n\t"\
		"movaps		0x190(%%rax),%%xmm3		\n\t"\
		"mulpd		    (%%rdi),%%xmm4		\n\t"\
		"mulpd		0x10(%%rdi),%%xmm6		\n\t"\
		"mulpd		0x10(%%rdi),%%xmm1		\n\t"\
		"mulpd		    (%%rdi),%%xmm3		\n\t"\
		"mulpd		    (%%rdi),%%xmm5		\n\t"\
		"mulpd		0x10(%%rdi),%%xmm7		\n\t"\
		"mulpd		0x10(%%rdi),%%xmm0		\n\t"\
		"mulpd		(%%rdi),%%xmm2			\n\t"\
		"subpd		%%xmm1,%%xmm4			\n\t"\
		"subpd		%%xmm3,%%xmm6			\n\t"\
		"addpd		%%xmm0,%%xmm5			\n\t"\
		"addpd		%%xmm2,%%xmm7			\n\t"\
		"subpd		%%xmm6,%%xmm4			\n\t"\
		"subpd		%%xmm7,%%xmm5			\n\t"\
		"addpd		%%xmm6,%%xmm6			\n\t"\
		"addpd		%%xmm7,%%xmm7			\n\t"\
		"addpd		%%xmm4,%%xmm6			\n\t"\
		"addpd		%%xmm5,%%xmm7			\n\t"\
		"movaps		0x080(%%rax),%%xmm2		\n\t"\
		"movaps		0x090(%%rax),%%xmm3		\n\t"\
		"subpd		0x090(%%rax),%%xmm2		\n\t"\
		"addpd		0x080(%%rax),%%xmm3		\n\t"\
		"mulpd		(%%rsi),%%xmm2			\n\t"\
		"mulpd		(%%rsi),%%xmm3			\n\t"\
		"movaps		     (%%rax),%%xmm0		\n\t"\
		"movaps		0x010(%%rax),%%xmm1		\n\t"\
		"subpd		%%xmm2,%%xmm0			\n\t"\
		"subpd		%%xmm3,%%xmm1			\n\t"\
		"addpd		     (%%rax),%%xmm2		\n\t"\
		"addpd		0x010(%%rax),%%xmm3		\n\t"\
		"subpd		%%xmm6,%%xmm2			\n\t"\
		"subpd		%%xmm7,%%xmm3			\n\t"\
		"addpd		%%xmm6,%%xmm6			\n\t"\
		"addpd		%%xmm7,%%xmm7			\n\t"\
		"movaps		%%xmm2,    (%%rcx)		\n\t"\
		"movaps		%%xmm3,0x10(%%rcx)		\n\t"\
		"addpd		%%xmm2,%%xmm6			\n\t"\
		"addpd		%%xmm3,%%xmm7			\n\t"\
		"movaps		%%xmm6,    (%%rax)		\n\t"\
		"movaps		%%xmm7,0x10(%%rax)		\n\t"\
		"subpd		%%xmm5,%%xmm0			\n\t"\
		"subpd		%%xmm4,%%xmm1			\n\t"\
		"addpd		%%xmm5,%%xmm5			\n\t"\
		"addpd		%%xmm4,%%xmm4			\n\t"\
		"movaps		%%xmm0,    (%%rbx)		\n\t"\
		"movaps		%%xmm1,0x10(%%rdx)		\n\t"\
		"addpd		%%xmm0,		%%xmm5		\n\t"\
		"addpd		%%xmm1,		%%xmm4		\n\t"\
		"movaps		%%xmm5,    (%%rdx)		\n\t"\
		"movaps		%%xmm4,0x10(%%rbx)		\n\t"\
		"/*...Block 4: t7,15,23,31 */		\n\t"\
		"addq		$0x40,%%rax				\n\t"\
		"addq		$0x40,%%rbx				\n\t"\
		"addq		$0x40,%%rcx				\n\t"\
		"addq		$0x40,%%rdx				\n\t"\
		"movaps		0x100(%%rax),%%xmm4		\n\t"\
		"movaps		0x180(%%rax),%%xmm6		\n\t"\
		"movaps		0x110(%%rax),%%xmm5		\n\t"\
		"movaps		0x190(%%rax),%%xmm7		\n\t"\
		"movaps		0x100(%%rax),%%xmm0		\n\t"\
		"movaps		0x180(%%rax),%%xmm2		\n\t"\
		"movaps		0x110(%%rax),%%xmm1		\n\t"\
		"movaps		0x190(%%rax),%%xmm3		\n\t"\
		"mulpd		0x10(%%rdi),%%xmm4		\n\t"\
		"mulpd		    (%%rdi),%%xmm6		\n\t"\
		"mulpd		    (%%rdi),%%xmm1		\n\t"\
		"mulpd		0x10(%%rdi),%%xmm3		\n\t"\
		"mulpd		0x10(%%rdi),%%xmm5		\n\t"\
		"mulpd		    (%%rdi),%%xmm7		\n\t"\
		"mulpd		    (%%rdi),%%xmm0		\n\t"\
		"mulpd		0x10(%%rdi),%%xmm2		\n\t"\
		"subpd		%%xmm1,%%xmm4			\n\t"\
		"subpd		%%xmm3,%%xmm6			\n\t"\
		"addpd		%%xmm0,%%xmm5			\n\t"\
		"addpd		%%xmm2,%%xmm7			\n\t"\
		"subpd		%%xmm6,%%xmm4			\n\t"\
		"subpd		%%xmm7,%%xmm5			\n\t"\
		"addpd		%%xmm6,%%xmm6			\n\t"\
		"addpd		%%xmm7,%%xmm7			\n\t"\
		"addpd		%%xmm4,%%xmm6			\n\t"\
		"addpd		%%xmm5,%%xmm7			\n\t"\
		"movaps		0x080(%%rax),%%xmm2		\n\t"\
		"movaps		0x090(%%rax),%%xmm3		\n\t"\
		"addpd		0x090(%%rax),%%xmm2		\n\t"\
		"subpd		0x080(%%rax),%%xmm3		\n\t"\
		"mulpd		(%%rsi),%%xmm2			\n\t"\
		"mulpd		(%%rsi),%%xmm3			\n\t"\
		"movaps		     (%%rax),%%xmm0		\n\t"\
		"movaps		0x010(%%rax),%%xmm1		\n\t"\
		"subpd		%%xmm2,%%xmm0			\n\t"\
		"subpd		%%xmm3,%%xmm1			\n\t"\
		"addpd		     (%%rax),%%xmm2		\n\t"\
		"addpd		0x010(%%rax),%%xmm3		\n\t"\
		"subpd		%%xmm4,%%xmm0			\n\t"\
		"subpd		%%xmm5,%%xmm1			\n\t"\
		"addpd		%%xmm4,%%xmm4			\n\t"\
		"addpd		%%xmm5,%%xmm5			\n\t"\
		"movaps		%%xmm0,    (%%rcx)		\n\t"\
		"movaps		%%xmm1,0x10(%%rcx)		\n\t"\
		"addpd		%%xmm0,%%xmm4			\n\t"\
		"addpd		%%xmm1,%%xmm5			\n\t"\
		"movaps		%%xmm4,    (%%rax)		\n\t"\
		"movaps		%%xmm5,0x10(%%rax)		\n\t"\
		"subpd		%%xmm7,%%xmm2			\n\t"\
		"subpd		%%xmm6,%%xmm3			\n\t"\
		"addpd		%%xmm7,%%xmm7			\n\t"\
		"addpd		%%xmm6,%%xmm6			\n\t"\
		"movaps		%%xmm2,    (%%rbx)		\n\t"\
		"movaps		%%xmm3,0x10(%%rdx)		\n\t"\
		"addpd		%%xmm2,%%xmm7			\n\t"\
		"addpd		%%xmm3,%%xmm6			\n\t"\
		"movaps		%%xmm7,    (%%rdx)		\n\t"\
		"movaps		%%xmm6,0x10(%%rbx)		\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
		: "rax","rbx","rcx","rdx","rdi","rsi","r10","r11","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}


	#define SSE2_RADIX16_WRAPPER_DIT(Xadd0,Xadd1,Xr1,Xr9,Xr17,Xr25,Xisrt2,Xcc0,Xc1,Xc2,Xc3,Xc4,Xc5,Xc6,Xc7,Xc8,Xc9,Xc10,Xc11,Xc12,Xc13,Xc14,Xc15)\
	{\
	__asm__ volatile (\
		"/*...Block 1: */						\n\t"\
		"movq	%[__r1],%%rax					\n\t"\
		"movq	%%rax,%%rbx						\n\t"\
		"movq	%%rax,%%rcx						\n\t"\
		"movq	%%rax,%%rdx						\n\t"\
		"addq	$0x100,%%rbx					\n\t"\
		"addq	$0x080,%%rcx					\n\t"\
		"addq	$0x180,%%rdx					\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"movaps	    (%%rax),%%xmm0				\n\t"\
		"movaps	0x10(%%rax),%%xmm1				\n\t"\
		"movaps	    (%%rax),%%xmm2				\n\t"\
		"movaps	0x10(%%rax),%%xmm3				\n\t"\
		"addpd	    (%%rbx),%%xmm0				\n\t"\
		"addpd	0x10(%%rbx),%%xmm1				\n\t"\
		"subpd	    (%%rbx),%%xmm2				\n\t"\
		"subpd	0x10(%%rbx),%%xmm3				\n\t"\
		"movaps	    (%%rcx),%%xmm4				\n\t"\
		"movaps	0x10(%%rcx),%%xmm5				\n\t"\
		"movaps	    (%%rcx),%%xmm6				\n\t"\
		"movaps	0x10(%%rcx),%%xmm7				\n\t"\
		"addpd	    (%%rdx),%%xmm4				\n\t"\
		"addpd	0x10(%%rdx),%%xmm5				\n\t"\
		"subpd	    (%%rdx),%%xmm6				\n\t"\
		"subpd	0x10(%%rdx),%%xmm7				\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t"\
		"movaps	%%xmm0,     (%%rbx)				\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)				\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t"\
		"movaps	%%xmm4,     (%%rax)				\n\t"\
		"movaps	%%xmm5,0x010(%%rax)				\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t"\
		"movaps	%%xmm2,     (%%rdx)				\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)				\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t"\
		"movaps	%%xmm7,     (%%rcx)				\n\t"\
		"movaps	%%xmm6,0x010(%%rdx)				\n\t"\
		"/*...Block 2: */						\n\t"\
		"addq	$0x040,%%rax					\n\t"\
		"addq	$0x040,%%rbx					\n\t"\
		"addq	$0x040,%%rcx					\n\t"\
		"addq	$0x040,%%rdx					\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"movaps	    (%%rax),%%xmm0				\n\t"\
		"movaps	0x10(%%rax),%%xmm1				\n\t"\
		"movaps	    (%%rax),%%xmm2				\n\t"\
		"movaps	0x10(%%rax),%%xmm3				\n\t"\
		"addpd	    (%%rbx),%%xmm0				\n\t"\
		"addpd	0x10(%%rbx),%%xmm1				\n\t"\
		"subpd	    (%%rbx),%%xmm2				\n\t"\
		"subpd	0x10(%%rbx),%%xmm3				\n\t"\
		"movaps	    (%%rcx),%%xmm4				\n\t"\
		"movaps	0x10(%%rcx),%%xmm5				\n\t"\
		"movaps	    (%%rcx),%%xmm6				\n\t"\
		"movaps	0x10(%%rcx),%%xmm7				\n\t"\
		"addpd	    (%%rdx),%%xmm4				\n\t"\
		"addpd	0x10(%%rdx),%%xmm5				\n\t"\
		"subpd	    (%%rdx),%%xmm6				\n\t"\
		"subpd	0x10(%%rdx),%%xmm7				\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t"\
		"movaps	%%xmm0,     (%%rbx)				\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)				\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t"\
		"movaps	%%xmm4,     (%%rax)				\n\t"\
		"movaps	%%xmm5,0x010(%%rax)				\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t"\
		"movaps	%%xmm2,     (%%rdx)				\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)				\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t"\
		"movaps	%%xmm7,     (%%rcx)				\n\t"\
		"movaps	%%xmm6,0x010(%%rdx)				\n\t"\
		"/*...Block 3: */						\n\t"\
		"subq	$0x020,%%rax					\n\t"\
		"subq	$0x020,%%rbx					\n\t"\
		"subq	$0x020,%%rcx					\n\t"\
		"subq	$0x020,%%rdx					\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"movaps	    (%%rax),%%xmm0				\n\t"\
		"movaps	0x10(%%rax),%%xmm1				\n\t"\
		"movaps	    (%%rax),%%xmm2				\n\t"\
		"movaps	0x10(%%rax),%%xmm3				\n\t"\
		"addpd	    (%%rbx),%%xmm0				\n\t"\
		"addpd	0x10(%%rbx),%%xmm1				\n\t"\
		"subpd	    (%%rbx),%%xmm2				\n\t"\
		"subpd	0x10(%%rbx),%%xmm3				\n\t"\
		"movaps	    (%%rcx),%%xmm4				\n\t"\
		"movaps	0x10(%%rcx),%%xmm5				\n\t"\
		"movaps	    (%%rcx),%%xmm6				\n\t"\
		"movaps	0x10(%%rcx),%%xmm7				\n\t"\
		"addpd	    (%%rdx),%%xmm4				\n\t"\
		"addpd	0x10(%%rdx),%%xmm5				\n\t"\
		"subpd	    (%%rdx),%%xmm6				\n\t"\
		"subpd	0x10(%%rdx),%%xmm7				\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t"\
		"movaps	%%xmm0,     (%%rbx)				\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)				\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t"\
		"movaps	%%xmm4,     (%%rax)				\n\t"\
		"movaps	%%xmm5,0x010(%%rax)				\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t"\
		"movaps	%%xmm2,     (%%rdx)				\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)				\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t"\
		"movaps	%%xmm7,     (%%rcx)				\n\t"\
		"movaps	%%xmm6,0x010(%%rdx)				\n\t"\
		"/*...Block 4: */						\n\t"\
		"addq	$0x040,%%rax					\n\t"\
		"addq	$0x040,%%rbx					\n\t"\
		"addq	$0x040,%%rcx					\n\t"\
		"addq	$0x040,%%rdx					\n\t"\
		"/* SSE2_RADIX4_DIT_IN_PLACE(): */		\n\t"\
		"movaps	    (%%rax),%%xmm0				\n\t"\
		"movaps	0x10(%%rax),%%xmm1				\n\t"\
		"movaps	    (%%rax),%%xmm2				\n\t"\
		"movaps	0x10(%%rax),%%xmm3				\n\t"\
		"addpd	    (%%rbx),%%xmm0				\n\t"\
		"addpd	0x10(%%rbx),%%xmm1				\n\t"\
		"subpd	    (%%rbx),%%xmm2				\n\t"\
		"subpd	0x10(%%rbx),%%xmm3				\n\t"\
		"movaps	    (%%rcx),%%xmm4				\n\t"\
		"movaps	0x10(%%rcx),%%xmm5				\n\t"\
		"movaps	    (%%rcx),%%xmm6				\n\t"\
		"movaps	0x10(%%rcx),%%xmm7				\n\t"\
		"addpd	    (%%rdx),%%xmm4				\n\t"\
		"addpd	0x10(%%rdx),%%xmm5				\n\t"\
		"subpd	    (%%rdx),%%xmm6				\n\t"\
		"subpd	0x10(%%rdx),%%xmm7				\n\t"\
		"subpd	%%xmm4,%%xmm0					\n\t"\
		"subpd	%%xmm5,%%xmm1					\n\t"\
		"movaps	%%xmm0,     (%%rbx)				\n\t"\
		"movaps	%%xmm1,0x010(%%rbx)				\n\t"\
		"addpd	%%xmm4,%%xmm4					\n\t"\
		"addpd	%%xmm5,%%xmm5					\n\t"\
		"addpd	%%xmm0,%%xmm4					\n\t"\
		"addpd	%%xmm1,%%xmm5					\n\t"\
		"movaps	%%xmm4,     (%%rax)				\n\t"\
		"movaps	%%xmm5,0x010(%%rax)				\n\t"\
		"subpd	%%xmm7,%%xmm2					\n\t"\
		"subpd	%%xmm6,%%xmm3					\n\t"\
		"movaps	%%xmm2,     (%%rdx)				\n\t"\
		"movaps	%%xmm3,0x010(%%rcx)				\n\t"\
		"addpd	%%xmm7,%%xmm7					\n\t"\
		"addpd	%%xmm6,%%xmm6					\n\t"\
		"addpd	%%xmm2,%%xmm7					\n\t"\
		"addpd	%%xmm3,%%xmm6					\n\t"\
		"movaps	%%xmm7,     (%%rcx)				\n\t"\
		"movaps	%%xmm6,0x010(%%rdx)				\n\t"\
		"/***************************************************************************************************/\n\t"\
		"/*..and now do four more radix-4 transforms, including the internal and external twiddle factors.  */\n\t"\
		"/*  Write even-index 16-byte output pairs to a(j1), odd-index to a(j2), unpack same as on inputs.  */\n\t"\
		"/*  We do the last 2 radix-4 blocks first, to make the unpack-interleaving of outputs smoother.    */\n\t"\
		"/***************************************************************************************************/\n\t"\
		"/* Main-array addresses still in add0,1, no need to re-init: */									  \n\t"\
		"/*...Block 3: t3,11,19,27 -> r9,13,11,15: */	\n\t"\
		"movq		%[__r9],%%rax				\n\t"\
		"movq		%[__isrt2],%%rbx			\n\t"\
		"movq		%[__cc0],%%rcx				\n\t"\
		"movaps		0x020(%%rax),%%xmm4			\n\t"\
		"movaps		0x060(%%rax),%%xmm0			\n\t"\
		"movaps		0x030(%%rax),%%xmm5			\n\t"\
		"movaps		0x070(%%rax),%%xmm1			\n\t"\
		"movaps		0x020(%%rax),%%xmm6			\n\t"\
		"movaps		0x060(%%rax),%%xmm2			\n\t"\
		"movaps		0x030(%%rax),%%xmm7			\n\t"\
		"movaps		0x070(%%rax),%%xmm3			\n\t"\
		"mulpd		    (%%rcx),%%xmm4			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm0			\n\t"\
		"mulpd		    (%%rcx),%%xmm5			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm1			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm6			\n\t"\
		"mulpd		    (%%rcx),%%xmm2			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm7			\n\t"\
		"mulpd		    (%%rcx),%%xmm3			\n\t"\
		"subpd		%%xmm6,%%xmm5				\n\t"\
		"subpd		%%xmm2,%%xmm1				\n\t"\
		"addpd		%%xmm7,%%xmm4				\n\t"\
		"addpd		%%xmm3,%%xmm0				\n\t"\
		"movaps		%%xmm5,%%xmm7				\n\t"\
		"movaps		%%xmm4,%%xmm6				\n\t"\
		"addpd		%%xmm0,%%xmm4				\n\t"\
		"addpd		%%xmm1,%%xmm5				\n\t"\
		"subpd		%%xmm0,%%xmm6				\n\t"\
		"subpd		%%xmm1,%%xmm7				\n\t"\
		"movaps		0x040(%%rax),%%xmm2			\n\t"\
		"movaps		0x050(%%rax),%%xmm3			\n\t"\
		"movaps		     (%%rax),%%xmm0			\n\t"\
		"movaps		0x010(%%rax),%%xmm1			\n\t"\
		"addpd		0x050(%%rax),%%xmm2			\n\t"\
		"subpd		0x040(%%rax),%%xmm3			\n\t"\
		"mulpd		(%%rbx),%%xmm2				\n\t"\
		"mulpd		(%%rbx),%%xmm3				\n\t"\
		"subpd		%%xmm2,%%xmm0				\n\t"\
		"subpd		%%xmm3,%%xmm1				\n\t"\
		"addpd		%%xmm2,%%xmm2				\n\t"\
		"addpd		%%xmm3,%%xmm3				\n\t"\
		"addpd		%%xmm0,%%xmm2				\n\t"\
		"addpd		%%xmm1,%%xmm3				\n\t"\
		"movq		%[__add1],%%rbx				\n\t"\
		"movq		%[__c1],%%rcx				\n\t"\
		"movq		%[__c9],%%rdx				\n\t"\
		"subpd		%%xmm4,%%xmm2				\n\t"\
		"subpd		%%xmm5,%%xmm3				\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t"\
		"addpd		%%xmm2,%%xmm4				\n\t"\
		"addpd		%%xmm3,%%xmm5				\n\t"\
		"movaps		%%xmm2,     (%%rax)			\n\t"\
		"movaps		%%xmm3,0x010(%%rax)			\n\t"\
		"movaps		%%xmm4,%%xmm2				\n\t"\
		"movaps		%%xmm5,%%xmm3				\n\t"\
		"mulpd		    (%%rcx),%%xmm4			\n\t"\
		"mulpd		    (%%rcx),%%xmm5			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm2			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm3			\n\t"\
		"subpd		%%xmm2,%%xmm5				\n\t"\
		"addpd		%%xmm3,%%xmm4				\n\t"\
		"movaps		%%xmm5,0x10(%%rbx)			\n\t"\
		"movaps		%%xmm4,    (%%rbx)			\n\t"\
		"movaps		     (%%rax),%%xmm4			\n\t"\
		"movaps		0x010(%%rax),%%xmm5			\n\t"\
		"movaps		%%xmm4,%%xmm2				\n\t"\
		"movaps		%%xmm5,%%xmm3				\n\t"\
		"mulpd		    (%%rdx),%%xmm4			\n\t"\
		"mulpd		    (%%rdx),%%xmm5			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm2			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm3			\n\t"\
		"subpd		%%xmm2,%%xmm5				\n\t"\
		"addpd		%%xmm3,%%xmm4				\n\t"\
		"movaps		%%xmm5,0x90(%%rbx)			\n\t"\
		"movaps		%%xmm4,0x80(%%rbx)			\n\t"\
		"movq		%[__c5],%%rcx				\n\t"\
		"movq		%[__c13],%%rdx				\n\t"\
		"subpd		%%xmm7,%%xmm0				\n\t"\
		"subpd		%%xmm6,%%xmm1				\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t"\
		"addpd		%%xmm0,%%xmm7				\n\t"\
		"addpd		%%xmm1,%%xmm6				\n\t"\
		"movaps		%%xmm7,%%xmm4				\n\t"\
		"movaps		%%xmm1,%%xmm5				\n\t"\
		"mulpd		    (%%rcx),%%xmm7			\n\t"\
		"mulpd		    (%%rcx),%%xmm1			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm4			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm5			\n\t"\
		"subpd		%%xmm4,%%xmm1				\n\t"\
		"addpd		%%xmm5,%%xmm7				\n\t"\
		"movaps		%%xmm1,0x50(%%rbx)			\n\t"\
		"movaps		%%xmm7,0x40(%%rbx)			\n\t"\
		"movaps		%%xmm0,%%xmm4				\n\t"\
		"movaps		%%xmm6,%%xmm5				\n\t"\
		"mulpd		    (%%rdx),%%xmm0			\n\t"\
		"mulpd		    (%%rdx),%%xmm6			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm4			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm5			\n\t"\
		"subpd		%%xmm4,%%xmm6				\n\t"\
		"addpd		%%xmm5,%%xmm0				\n\t"\
		"movaps		%%xmm6,0xd0(%%rbx)			\n\t"\
		"movaps		%%xmm0,0xc0(%%rbx)			\n\t"\
		"/*...Block 4: t7,15,23,31 -> r25,29,27,31: */	\n\t"\
		"movq		%[__r25],%%rax				\n\t"\
		"movq		%[__isrt2],%%rbx			\n\t"\
		"movq		%[__cc0],%%rcx				\n\t"\
		"movaps		0x020(%%rax),%%xmm4			\n\t"\
		"movaps		0x060(%%rax),%%xmm0			\n\t"\
		"movaps		0x030(%%rax),%%xmm5			\n\t"\
		"movaps		0x070(%%rax),%%xmm1			\n\t"\
		"movaps		0x020(%%rax),%%xmm6			\n\t"\
		"movaps		0x060(%%rax),%%xmm2			\n\t"\
		"movaps		0x030(%%rax),%%xmm7			\n\t"\
		"movaps		0x070(%%rax),%%xmm3			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm4			\n\t"\
		"mulpd		    (%%rcx),%%xmm0			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm5			\n\t"\
		"mulpd		    (%%rcx),%%xmm1			\n\t"\
		"mulpd		    (%%rcx),%%xmm6			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm2			\n\t"\
		"mulpd		    (%%rcx),%%xmm7			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm3			\n\t"\
		"subpd		%%xmm6,%%xmm5				\n\t"\
		"subpd		%%xmm2,%%xmm1				\n\t"\
		"addpd		%%xmm7,%%xmm4				\n\t"\
		"addpd		%%xmm3,%%xmm0				\n\t"\
		"movaps		%%xmm5,%%xmm7				\n\t"\
		"movaps		%%xmm4,%%xmm6				\n\t"\
		"addpd		%%xmm0,%%xmm4				\n\t"\
		"addpd		%%xmm1,%%xmm5				\n\t"\
		"subpd		%%xmm0,%%xmm6				\n\t"\
		"subpd		%%xmm1,%%xmm7				\n\t"\
		"movaps		0x040(%%rax),%%xmm2			\n\t"\
		"movaps		0x050(%%rax),%%xmm3			\n\t"\
		"movaps		     (%%rax),%%xmm0			\n\t"\
		"movaps		0x010(%%rax),%%xmm1			\n\t"\
		"subpd		0x050(%%rax),%%xmm2			\n\t"\
		"addpd		0x040(%%rax),%%xmm3			\n\t"\
		"mulpd		(%%rbx),%%xmm2				\n\t"\
		"mulpd		(%%rbx),%%xmm3				\n\t"\
		"subpd		%%xmm2,%%xmm0				\n\t"\
		"subpd		%%xmm3,%%xmm1				\n\t"\
		"addpd		%%xmm2,%%xmm2				\n\t"\
		"addpd		%%xmm3,%%xmm3				\n\t"\
		"addpd		%%xmm0,%%xmm2				\n\t"\
		"addpd		%%xmm1,%%xmm3				\n\t"\
		"movq		%[__add1],%%rbx				\n\t"\
		"movq		%[__c3],%%rcx				\n\t"\
		"movq		%[__c11],%%rdx				\n\t"\
		"subpd		%%xmm6,%%xmm0				\n\t"\
		"subpd		%%xmm7,%%xmm1				\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t"\
		"addpd		%%xmm0,%%xmm6				\n\t"\
		"addpd		%%xmm1,%%xmm7				\n\t"\
		"movaps		%%xmm0,     (%%rax)			\n\t"\
		"movaps		%%xmm1,0x010(%%rax)			\n\t"\
		"movaps		%%xmm6,%%xmm0				\n\t"\
		"movaps		%%xmm7,%%xmm1				\n\t"\
		"mulpd		    (%%rcx),%%xmm6			\n\t"\
		"mulpd		    (%%rcx),%%xmm7			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm0			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm1			\n\t"\
		"subpd		%%xmm0,%%xmm7				\n\t"\
		"addpd		%%xmm1,%%xmm6				\n\t"\
		"movaps		%%xmm7,0x30(%%rbx)			\n\t"\
		"movaps		%%xmm6,0x20(%%rbx)			\n\t"\
		"movaps		     (%%rax),%%xmm6			\n\t"\
		"movaps		0x010(%%rax),%%xmm7			\n\t"\
		"movaps		%%xmm6,%%xmm0				\n\t"\
		"movaps		%%xmm7,%%xmm1				\n\t"\
		"mulpd		    (%%rdx),%%xmm6			\n\t"\
		"mulpd		    (%%rdx),%%xmm7			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm0			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm1			\n\t"\
		"subpd		%%xmm0,%%xmm7				\n\t"\
		"addpd		%%xmm1,%%xmm6				\n\t"\
		"movaps		%%xmm7,0xb0(%%rbx)			\n\t"\
		"movaps		%%xmm6,0xa0(%%rbx)			\n\t"\
		"movq		%[__c7],%%rcx				\n\t"\
		"movq		%[__c15],%%rdx				\n\t"\
		"subpd		%%xmm5,%%xmm2				\n\t"\
		"subpd		%%xmm4,%%xmm3				\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t"\
		"addpd		%%xmm2,%%xmm5				\n\t"\
		"addpd		%%xmm3,%%xmm4				\n\t"\
		"movaps		%%xmm5,%%xmm0				\n\t"\
		"movaps		%%xmm3,%%xmm1				\n\t"\
		"mulpd		    (%%rcx),%%xmm5			\n\t"\
		"mulpd		    (%%rcx),%%xmm3			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm0			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm1			\n\t"\
		"subpd		%%xmm0,%%xmm3				\n\t"\
		"addpd		%%xmm1,%%xmm5				\n\t"\
		"movaps		%%xmm3,0x70(%%rbx)			\n\t"\
		"movaps		%%xmm5,0x60(%%rbx)			\n\t"\
		"movaps		%%xmm2,%%xmm0				\n\t"\
		"movaps		%%xmm4,%%xmm1				\n\t"\
		"mulpd		    (%%rdx),%%xmm2			\n\t"\
		"mulpd		    (%%rdx),%%xmm4			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm0			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm1			\n\t"\
		"subpd		%%xmm0,%%xmm4				\n\t"\
		"addpd		%%xmm1,%%xmm2				\n\t"\
		"movaps		%%xmm4,0xf0(%%rbx)			\n\t"\
		"movaps		%%xmm2,0xe0(%%rbx)			\n\t"\
		"/*...Block 1: t1,9,17,25 -> r1,5,3,7: */		\n\t"\
		"movq		%[__r1],%%rax				\n\t"\
		"movaps		     (%%rax),%%xmm0			\n\t"\
		"movaps		0x010(%%rax),%%xmm1			\n\t"\
		"movaps		0x040(%%rax),%%xmm2			\n\t"\
		"movaps		0x050(%%rax),%%xmm3			\n\t"\
		"subpd		0x040(%%rax),%%xmm0			\n\t"\
		"subpd		0x050(%%rax),%%xmm1			\n\t"\
		"addpd		     (%%rax),%%xmm2			\n\t"\
		"addpd		0x010(%%rax),%%xmm3			\n\t"\
		"movaps		0x020(%%rax),%%xmm4			\n\t"\
		"movaps		0x030(%%rax),%%xmm5			\n\t"\
		"movaps		0x060(%%rax),%%xmm6			\n\t"\
		"movaps		0x070(%%rax),%%xmm7			\n\t"\
		"subpd		0x060(%%rax),%%xmm4			\n\t"\
		"subpd		0x070(%%rax),%%xmm5			\n\t"\
		"addpd		0x020(%%rax),%%xmm6			\n\t"\
		"addpd		0x030(%%rax),%%xmm7			\n\t"\
		"movq		%[__add0],%%rax				\n\t"\
		"movq		%[__c8],%%rdx				\n\t"\
		"addpd		%%xmm6,%%xmm2				\n\t"\
		"addpd		%%xmm7,%%xmm3				\n\t"\
		"movaps		%%xmm2,    (%%rax)			\n\t"\
		"movaps		%%xmm3,0x10(%%rax)			\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t"\
		"subpd		%%xmm6,%%xmm2				\n\t"\
		"subpd		%%xmm7,%%xmm3				\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t"\
		"mulpd		    (%%rdx),%%xmm2			\n\t"\
		"mulpd		    (%%rdx),%%xmm3			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7			\n\t"\
		"subpd		%%xmm6,%%xmm3				\n\t"\
		"addpd		%%xmm7,%%xmm2				\n\t"\
		"movq		%[__add1],%%rcx				\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t"\
		"unpckhpd	0x90(%%rcx),%%xmm7			\n\t"\
		"unpcklpd	0x90(%%rcx),%%xmm3			\n\t"\
		"movaps		%%xmm7,0x90(%%rcx)			\n\t"\
		"unpckhpd	0x80(%%rcx),%%xmm6			\n\t"\
		"unpcklpd	0x80(%%rcx),%%xmm2			\n\t"\
		"movaps		%%xmm6,0x80(%%rcx)			\n\t"\
		"movaps		%%xmm3,0x90(%%rax)			\n\t"\
		"movaps		%%xmm2,0x80(%%rax)			\n\t"\
		"movaps		0x10(%%rax),%%xmm3			\n\t"\
		"movaps		(%%rax),%%xmm2				\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t"\
		"unpckhpd	0x10(%%rcx),%%xmm7			\n\t"\
		"unpcklpd	0x10(%%rcx),%%xmm3			\n\t"\
		"movaps		%%xmm7,0x10(%%rcx)			\n\t"\
		"unpckhpd	(%%rcx),%%xmm6				\n\t"\
		"unpcklpd	(%%rcx),%%xmm2				\n\t"\
		"movaps		%%xmm6,    (%%rcx)			\n\t"\
		"movaps		%%xmm3,0x10(%%rax)			\n\t"\
		"movaps		%%xmm2,    (%%rax)			\n\t"\
		"movq		%[__c4],%%rcx				\n\t"\
		"movq		%[__c12],%%rdx				\n\t"\
		"addpd		%%xmm5,%%xmm0				\n\t"\
		"subpd		%%xmm4,%%xmm1				\n\t"\
		"movaps		%%xmm0,%%xmm2				\n\t"\
		"movaps		%%xmm1,%%xmm3				\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t"\
		"movaps		%%xmm0,%%xmm6				\n\t"\
		"movaps		%%xmm1,%%xmm7				\n\t"\
		"mulpd		    (%%rcx),%%xmm2			\n\t"\
		"mulpd		    (%%rcx),%%xmm3			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm6			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm7			\n\t"\
		"subpd		%%xmm6,%%xmm3				\n\t"\
		"addpd		%%xmm7,%%xmm2				\n\t"\
		"movq		%[__add1],%%rcx				\n\t"\
		"movaps		%%xmm3,%%xmm7				\n\t"\
		"movaps		%%xmm2,%%xmm6				\n\t"\
		"unpckhpd	0x50(%%rcx),%%xmm7			\n\t"\
		"unpcklpd	0x50(%%rcx),%%xmm3			\n\t"\
		"movaps		%%xmm7,0x50(%%rcx)			\n\t"\
		"unpckhpd	0x40(%%rcx),%%xmm6			\n\t"\
		"unpcklpd	0x40(%%rcx),%%xmm2			\n\t"\
		"movaps		%%xmm6,0x40(%%rcx)			\n\t"\
		"movaps		%%xmm3,0x50(%%rax)			\n\t"\
		"movaps		%%xmm2,0x40(%%rax)			\n\t"\
		"subpd		%%xmm5,%%xmm0				\n\t"\
		"addpd		%%xmm4,%%xmm1				\n\t"\
		"movaps		%%xmm0,%%xmm6				\n\t"\
		"movaps		%%xmm1,%%xmm7				\n\t"\
		"mulpd		    (%%rdx),%%xmm0			\n\t"\
		"mulpd		    (%%rdx),%%xmm1			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm6			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm7			\n\t"\
		"subpd		%%xmm6,%%xmm1				\n\t"\
		"addpd		%%xmm7,%%xmm0				\n\t"\
		"movq		%[__add1],%%rcx				\n\t"\
		"movaps		%%xmm1,%%xmm7				\n\t"\
		"movaps		%%xmm0,%%xmm6				\n\t"\
		"unpckhpd	0xd0(%%rcx),%%xmm7			\n\t"\
		"unpcklpd	0xd0(%%rcx),%%xmm1			\n\t"\
		"movaps		%%xmm7,0xd0(%%rcx)			\n\t"\
		"unpckhpd	0xc0(%%rcx),%%xmm6			\n\t"\
		"unpcklpd	0xc0(%%rcx),%%xmm0			\n\t"\
		"movaps		%%xmm6,0xc0(%%rcx)			\n\t"\
		"movaps		%%xmm1,0xd0(%%rax)			\n\t"\
		"movaps		%%xmm0,0xc0(%%rax)			\n\t"\
		"/*...Block 3: t5,13,21,29 -> r17,21,19,23: */	\n\t"\
		"movq		%[__r17],%%rax				\n\t"\
		"movq		%[__isrt2],%%rbx			\n\t"\
		"movaps		(%%rbx),%%xmm2				\n\t"\
		"movaps		0x020(%%rax),%%xmm4			\n\t"\
		"movaps		0x030(%%rax),%%xmm5			\n\t"\
		"movaps		0x060(%%rax),%%xmm0			\n\t"\
		"movaps		0x070(%%rax),%%xmm1			\n\t"\
		"addpd		0x030(%%rax),%%xmm4			\n\t"\
		"subpd		0x020(%%rax),%%xmm5			\n\t"\
		"subpd		0x070(%%rax),%%xmm0			\n\t"\
		"addpd		0x060(%%rax),%%xmm1			\n\t"\
		"mulpd		%%xmm2,%%xmm4				\n\t"\
		"mulpd		%%xmm2,%%xmm5				\n\t"\
		"mulpd		%%xmm2,%%xmm0				\n\t"\
		"mulpd		%%xmm2,%%xmm1				\n\t"\
		"movaps		%%xmm4,%%xmm6				\n\t"\
		"movaps		%%xmm5,%%xmm7				\n\t"\
		"subpd		%%xmm0,%%xmm4				\n\t"\
		"subpd		%%xmm1,%%xmm5				\n\t"\
		"addpd		%%xmm0,%%xmm6				\n\t"\
		"addpd		%%xmm1,%%xmm7				\n\t"\
		"movaps		     (%%rax),%%xmm0			\n\t"\
		"movaps		0x010(%%rax),%%xmm1			\n\t"\
		"movaps		0x040(%%rax),%%xmm2			\n\t"\
		"movaps		0x050(%%rax),%%xmm3			\n\t"\
		"subpd		0x050(%%rax),%%xmm0			\n\t"\
		"subpd		0x040(%%rax),%%xmm1			\n\t"\
		"addpd		     (%%rax),%%xmm3			\n\t"\
		"addpd		0x010(%%rax),%%xmm2			\n\t"\
		"movq		%[__add0],%%rbx				\n\t"\
		"movq		%[__c2],%%rcx				\n\t"\
		"movq		%[__c10],%%rdx				\n\t"\
		"subpd		%%xmm4,%%xmm3				\n\t"\
		"subpd		%%xmm5,%%xmm1				\n\t"\
		"addpd		%%xmm4,%%xmm4				\n\t"\
		"addpd		%%xmm5,%%xmm5				\n\t"\
		"addpd		%%xmm3,%%xmm4				\n\t"\
		"addpd		%%xmm1,%%xmm5				\n\t"\
		"movaps		%%xmm3,     (%%rax)			\n\t"\
		"movaps		%%xmm1,0x010(%%rax)			\n\t"\
		"movaps		%%xmm4,%%xmm3				\n\t"\
		"movaps		%%xmm5,%%xmm1				\n\t"\
		"mulpd		    (%%rcx),%%xmm4			\n\t"\
		"mulpd		    (%%rcx),%%xmm5			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm3			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm1			\n\t"\
		"subpd		%%xmm3,%%xmm5				\n\t"\
		"addpd		%%xmm1,%%xmm4				\n\t"\
		"movq		%[__add1],%%rcx				\n\t"\
		"movaps		%%xmm5,%%xmm3				\n\t"\
		"movaps		%%xmm4,%%xmm1				\n\t"\
		"unpckhpd	0x30(%%rcx),%%xmm3			\n\t"\
		"unpcklpd	0x30(%%rcx),%%xmm5			\n\t"\
		"movaps		%%xmm3,0x30(%%rcx)			\n\t"\
		"unpckhpd	0x20(%%rcx),%%xmm1			\n\t"\
		"unpcklpd	0x20(%%rcx),%%xmm4			\n\t"\
		"movaps		%%xmm1,0x20(%%rcx)			\n\t"\
		"movaps		%%xmm5,0x30(%%rbx)			\n\t"\
		"movaps		%%xmm4,0x20(%%rbx)			\n\t"\
		"movaps		     (%%rax),%%xmm4			\n\t"\
		"movaps		0x010(%%rax),%%xmm5			\n\t"\
		"movaps		%%xmm4,%%xmm3				\n\t"\
		"movaps		%%xmm5,%%xmm1				\n\t"\
		"mulpd		    (%%rdx),%%xmm4			\n\t"\
		"mulpd		    (%%rdx),%%xmm5			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm3			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm1			\n\t"\
		"subpd		%%xmm3,%%xmm5				\n\t"\
		"addpd		%%xmm1,%%xmm4				\n\t"\
		"movaps		%%xmm5,%%xmm3				\n\t"\
		"movaps		%%xmm4,%%xmm1				\n\t"\
		"unpckhpd	0xb0(%%rcx),%%xmm3			\n\t"\
		"unpcklpd	0xb0(%%rcx),%%xmm5			\n\t"\
		"movaps		%%xmm3,0xb0(%%rcx)			\n\t"\
		"unpckhpd	0xa0(%%rcx),%%xmm1			\n\t"\
		"unpcklpd	0xa0(%%rcx),%%xmm4			\n\t"\
		"movaps		%%xmm1,0xa0(%%rcx)			\n\t"\
		"movaps		%%xmm5,0xb0(%%rbx)			\n\t"\
		"movaps		%%xmm4,0xa0(%%rbx)			\n\t"\
		"movq		%[__c6],%%rcx				\n\t"\
		"movq		%[__c14],%%rdx				\n\t"\
		"subpd		%%xmm7,%%xmm0				\n\t"\
		"subpd		%%xmm6,%%xmm2				\n\t"\
		"addpd		%%xmm7,%%xmm7				\n\t"\
		"addpd		%%xmm6,%%xmm6				\n\t"\
		"addpd		%%xmm0,%%xmm7				\n\t"\
		"addpd		%%xmm2,%%xmm6				\n\t"\
		"movaps		%%xmm7,%%xmm4				\n\t"\
		"movaps		%%xmm2,%%xmm5				\n\t"\
		"mulpd		    (%%rcx),%%xmm7			\n\t"\
		"mulpd		    (%%rcx),%%xmm2			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm4			\n\t"\
		"mulpd		0x10(%%rcx),%%xmm5			\n\t"\
		"subpd		%%xmm4,%%xmm2				\n\t"\
		"addpd		%%xmm5,%%xmm7				\n\t"\
		"movq		%[__add1],%%rcx				\n\t"\
		"movaps		%%xmm2,%%xmm5				\n\t"\
		"movaps		%%xmm7,%%xmm4				\n\t"\
		"unpckhpd	0x70(%%rcx),%%xmm5			\n\t"\
		"unpcklpd	0x70(%%rcx),%%xmm2			\n\t"\
		"movaps		%%xmm5,0x70(%%rcx)			\n\t"\
		"unpckhpd	0x60(%%rcx),%%xmm4			\n\t"\
		"unpcklpd	0x60(%%rcx),%%xmm7			\n\t"\
		"movaps		%%xmm4,0x60(%%rcx)			\n\t"\
		"movaps		%%xmm2,0x70(%%rbx)			\n\t"\
		"movaps		%%xmm7,0x60(%%rbx)			\n\t"\
		"movaps		%%xmm0,%%xmm4				\n\t"\
		"movaps		%%xmm6,%%xmm5				\n\t"\
		"mulpd		    (%%rdx),%%xmm0			\n\t"\
		"mulpd		    (%%rdx),%%xmm6			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm4			\n\t"\
		"mulpd		0x10(%%rdx),%%xmm5			\n\t"\
		"subpd		%%xmm4,%%xmm6				\n\t"\
		"addpd		%%xmm5,%%xmm0				\n\t"\
		"movaps		%%xmm6,%%xmm5				\n\t"\
		"movaps		%%xmm0,%%xmm4				\n\t"\
		"unpckhpd	0xf0(%%rcx),%%xmm5			\n\t"\
		"unpcklpd	0xf0(%%rcx),%%xmm6			\n\t"\
		"movaps		%%xmm5,0xf0(%%rcx)			\n\t"\
		"unpckhpd	0xe0(%%rcx),%%xmm4			\n\t"\
		"unpcklpd	0xe0(%%rcx),%%xmm0			\n\t"\
		"movaps		%%xmm4,0xe0(%%rcx)			\n\t"\
		"movaps		%%xmm6,0xf0(%%rbx)			\n\t"\
		"movaps		%%xmm0,0xe0(%%rbx)			\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r1] "m" (Xr1)\
		 ,[__r9] "m" (Xr9)\
		 ,[__r17] "m" (Xr17)\
		 ,[__r25] "m" (Xr25)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c1] "m" (Xc1)\
		 ,[__c2] "m" (Xc2)\
		 ,[__c3] "m" (Xc3)\
		 ,[__c4] "m" (Xc4)\
		 ,[__c5] "m" (Xc5)\
		 ,[__c6] "m" (Xc6)\
		 ,[__c7] "m" (Xc7)\
		 ,[__c8] "m" (Xc8)\
		 ,[__c9] "m" (Xc9)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c11] "m" (Xc11)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c13] "m" (Xc13)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c15] "m" (Xc15)\
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#endif	/* radix16_wrapper_square_gcc_h_included */

