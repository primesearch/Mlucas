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

/*******************************************************************************
   We now include this header file if it was not included before.
*******************************************************************************/
#ifndef radix32_wrapper_square_gcc_h_included
#define radix32_wrapper_square_gcc_h_included

	#define SSE2_RADIX32_WRAPPER_DIF(Xadd0,Xadd1,Xr00,Xr10,Xr20,Xr30,Xisrt2,Xcc0,Xc00,Xc01,Xc02,Xc03,Xc05,Xc07)\
	{\
	__asm__ volatile (\
		"/************************************************************************/\n\t"\
		"/* Forward DIF radix-32 pass on the interleaved block1 and block2 data: */\n\t"\
		"/************************************************************************/\n\t"\
		"/*...Block 0:	*/									\n\t		movq	%[__isrt2],%%rsi		\n\t"\
		"movq		%[__add0]	,%%rax						\n\t		movq		%%rsi,%%rdi			\n\t"\
		"movq		%[__add1]	,%%rbx						\n\t		addq	$0x480,%%rdi	/* two */	\n\t"\
		"/*****	SSE2_RADIX4_DIF_4WRAPPER(c00,c08,c10,c18,r00) ****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c04,c0C,c14,c1C,r08)	*****/\n\t"\
		"movq		%[__r00]	,%%rcx						\n\t		/*addq		$0x080,%%rcx // __r08 */	\n\t"\
		"movq		%[__c00]	,%%rdx						\n\t		/*addq		$0x080,%%rdx // __c04 */	\n\t"\
		"movaps		     (%%rax),%%xmm6						\n\t		movaps		0x40 (%%rax),%%xmm14		\n\t"\
		"movaps		     (%%rax),%%xmm0						\n\t		movaps		0x40 (%%rax),%%xmm8			\n\t"\
		"unpckhpd	     (%%rbx),%%xmm6						\n\t		unpckhpd	0x40 (%%rbx),%%xmm14		\n\t"\
		"unpcklpd	     (%%rbx),%%xmm0						\n\t		unpcklpd	0x40 (%%rbx),%%xmm8			\n\t"\
		"movaps		%%xmm6		,0x200(%%rcx)				\n\t		movaps		%%xmm14,0x280(%%rcx)		\n\t"\
		"movaps		0x010(%%rax),%%xmm7						\n\t		movaps		0x50 (%%rax),%%xmm15		\n\t"\
		"movaps		0x010(%%rax),%%xmm1						\n\t		movaps		0x50 (%%rax),%%xmm9			\n\t"\
		"unpckhpd	0x010(%%rbx),%%xmm7						\n\t		unpckhpd	0x50 (%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x010(%%rbx),%%xmm1						\n\t		unpcklpd	0x50 (%%rbx),%%xmm9			\n\t"\
		"movaps		%%xmm7		,0x210(%%rcx)				\n\t		movaps		%%xmm15,0x290(%%rcx)		\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		movaps		%%xmm8,%%xmm10				\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		movaps		%%xmm9,%%xmm11				\n\t"\
		"mulpd		     (%%rdx),%%xmm0						\n\t		mulpd		0x080(%%rdx),%%xmm8			\n\t"\
		"mulpd		     (%%rdx),%%xmm1						\n\t		mulpd		0x080(%%rdx),%%xmm9			\n\t"\
		"mulpd		0x010(%%rdx),%%xmm2						\n\t		mulpd		0x090(%%rdx),%%xmm10		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm3						\n\t		mulpd		0x090(%%rdx),%%xmm11		\n\t"\
		"addpd		%%xmm2		,%%xmm1						\n\t		addpd		%%xmm10,%%xmm9				\n\t"\
		"subpd		%%xmm3		,%%xmm0						\n\t		subpd		%%xmm11,%%xmm8				\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		movaps		%%xmm8,%%xmm10				\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		movaps		%%xmm9,%%xmm11				\n\t"\
		"addq		$0x020		,%%rdx																			\n\t"\
		"movaps		0x100(%%rax),%%xmm6						\n\t		movaps		0x140(%%rax),%%xmm14		\n\t"\
		"movaps		0x100(%%rax),%%xmm4						\n\t		movaps		0x140(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x100(%%rbx),%%xmm6						\n\t		unpckhpd	0x140(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x100(%%rbx),%%xmm4						\n\t		unpcklpd	0x140(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x220(%%rcx)				\n\t		movaps		%%xmm14,0x2a0(%%rcx)		\n\t"\
		"movaps		0x110(%%rax),%%xmm7						\n\t		movaps		0x150(%%rax),%%xmm15		\n\t"\
		"movaps		0x110(%%rax),%%xmm5						\n\t		movaps		0x150(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x110(%%rbx),%%xmm7						\n\t		unpckhpd	0x150(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x110(%%rbx),%%xmm5						\n\t		unpcklpd	0x150(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x230(%%rcx)				\n\t		movaps		%%xmm15,0x2b0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd		     (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd		     (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"addpd		%%xmm4		,%%xmm0						\n\t		addpd		%%xmm12,%%xmm8				\n\t"\
		"addpd		%%xmm5		,%%xmm1						\n\t		addpd		%%xmm13,%%xmm9				\n\t"\
		"subpd		%%xmm4		,%%xmm2						\n\t		subpd		%%xmm12,%%xmm10				\n\t"\
		"subpd		%%xmm5		,%%xmm3						\n\t		subpd		%%xmm13,%%xmm11				\n\t"\
		"addq		$0x040		,%%rdx																			\n\t"\
		"movaps		0x180(%%rax),%%xmm6						\n\t		movaps		0x1c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x180(%%rax),%%xmm4						\n\t		movaps		0x1c0(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x180(%%rbx),%%xmm6						\n\t		unpckhpd	0x1c0(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x180(%%rbx),%%xmm4						\n\t		unpcklpd	0x1c0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x260(%%rcx)				\n\t		movaps		%%xmm14,0x2e0(%%rcx)		\n\t"\
		"movaps		0x190(%%rax),%%xmm7						\n\t		movaps		0x1d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x190(%%rax),%%xmm5						\n\t		movaps		0x1d0(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x190(%%rbx),%%xmm7						\n\t		unpckhpd	0x1d0(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x190(%%rbx),%%xmm5						\n\t		unpcklpd	0x1d0(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x270(%%rcx)				\n\t		movaps		%%xmm15,0x2f0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd		     (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd		     (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"movaps		%%xmm5		,0x010(%%rcx)				\n\t		movaps		%%xmm13,0x090(%%rcx)		\n\t"\
		"movaps		%%xmm4		,     (%%rcx)				\n\t		movaps		%%xmm12,0x080(%%rcx)		\n\t"\
		"subq		$0x020		,%%rdx																			\n\t"\
		"movaps		0x080(%%rax),%%xmm6						\n\t		movaps		0x0c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x080(%%rax),%%xmm4						\n\t		movaps		0x0c0(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x080(%%rbx),%%xmm6						\n\t		unpckhpd	0x0c0(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x080(%%rbx),%%xmm4						\n\t		unpcklpd	0x0c0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x240(%%rcx)				\n\t		movaps		%%xmm14,0x2c0(%%rcx)		\n\t"\
		"movaps		0x090(%%rax),%%xmm7						\n\t		movaps		0x0d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x090(%%rax),%%xmm5						\n\t		movaps		0x0d0(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x090(%%rbx),%%xmm7						\n\t		unpckhpd	0x0d0(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x090(%%rbx),%%xmm5						\n\t		unpcklpd	0x0d0(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x250(%%rcx)				\n\t		movaps		%%xmm15,0x2d0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd		     (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd		     (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"subpd		     (%%rcx),%%xmm4						\n\t		subpd		0x080(%%rcx),%%xmm12		\n\t"\
		"subpd		0x010(%%rcx),%%xmm5						\n\t		subpd		0x090(%%rcx),%%xmm13		\n\t"\
		"addpd		     (%%rcx),%%xmm6						\n\t		addpd		0x080(%%rcx),%%xmm14		\n\t"\
		"addpd		0x010(%%rcx),%%xmm7						\n\t		addpd		0x090(%%rcx),%%xmm15		\n\t"\
		"subpd		%%xmm6		,%%xmm0						\n\t		subpd		%%xmm14,%%xmm8				\n\t"\
		"subpd		%%xmm7		,%%xmm1						\n\t		subpd		%%xmm15,%%xmm9				\n\t"\
		"/*movaps	%%xmm0		,0x040(%%rcx)	*/			\n\t	/*	movaps		%%xmm8,0x0c0(%%rcx)		*/	\n\t"\
		"/*movaps	%%xmm1		,0x050(%%rcx)	*/			\n\t	/*	movaps		%%xmm9,0x0d0(%%rcx)		*/	\n\t"\
		"mulpd		(%%rdi)		,%%xmm6						\n\t		mulpd		(%%rdi),%%xmm14				\n\t"\
		"mulpd		(%%rdi)		,%%xmm7						\n\t		mulpd		(%%rdi),%%xmm15				\n\t"\
		"addpd		%%xmm0		,%%xmm6						\n\t		addpd		%%xmm8,%%xmm14				\n\t"\
		"addpd		%%xmm1		,%%xmm7						\n\t		addpd		%%xmm9,%%xmm15				\n\t"\
		"/*movaps	%%xmm6		,     (%%rcx)	*/			\n\t		movaps		%%xmm14,0x080(%%rcx)		\n\t"\
		"/*movaps	%%xmm7		,0x010(%%rcx)	*/			\n\t		movaps		%%xmm15,0x090(%%rcx)		\n\t"\
		"subpd		%%xmm5		,%%xmm2						\n\t		subpd		%%xmm13,%%xmm10				\n\t"\
		"subpd		%%xmm4		,%%xmm3						\n\t		subpd		%%xmm12,%%xmm11				\n\t"\
		"/*movaps	%%xmm2		,0x020(%%rcx)	*/			\n\t		mulpd		(%%rdi),%%xmm13				\n\t"\
		"/*movaps	%%xmm3		,0x070(%%rcx)	*/			\n\t		mulpd		(%%rdi),%%xmm12				\n\t"\
		"mulpd		(%%rdi)		,%%xmm5						\n\t		addpd		%%xmm10,%%xmm13				\n\t"\
		"mulpd		(%%rdi)		,%%xmm4						\n\t		addpd		%%xmm11,%%xmm12				\n\t"\
		"addpd		%%xmm2		,%%xmm5						\n\t		movaps		%%xmm10,%%xmm14				\n\t"\
		"addpd		%%xmm3		,%%xmm4						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"/*movaps	%%xmm5		,0x060(%%rcx)	*/			\n\t		subpd		%%xmm12,%%xmm10				\n\t"\
		"/*movaps	%%xmm4		,0x030(%%rcx)	*/						subpd		%%xmm11,%%xmm13				\n\t"\
		"																addpd		%%xmm12,%%xmm14				\n\t"\
		"																addpd		%%xmm11,%%xmm15				\n\t"\
		"																mulpd		(%%rsi),%%xmm10				\n\t"\
		"																mulpd		(%%rsi),%%xmm13				\n\t"\
		"																mulpd		(%%rsi),%%xmm14				\n\t"\
		"																mulpd		(%%rsi),%%xmm15				\n\t"\
		"															/*	movaps		%%xmm10,0x0a0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm13,0x0e0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm14,0x0b0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm15,0x0f0(%%rcx)	*/	\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00)	*****/\n\t"\
		"																movaps		0x080(%%rcx),%%xmm11		\n\t"\
		"																movaps		0x090(%%rcx),%%xmm12		\n\t"\
		"subpd		%%xmm11		,%%xmm6						\n\t		subpd		%%xmm10		,%%xmm2			\n\t"\
		"subpd		%%xmm9		,%%xmm0						\n\t		subpd		%%xmm15		,%%xmm5			\n\t"\
		"subpd		%%xmm12		,%%xmm7						\n\t		subpd		%%xmm14		,%%xmm4			\n\t"\
		"subpd		%%xmm8		,%%xmm1						\n\t		subpd		%%xmm13		,%%xmm3			\n\t"\
		"mulpd		(%%rdi)		,%%xmm11					\n\t		mulpd		(%%rdi)		,%%xmm10		\n\t"\
		"mulpd		(%%rdi)		,%%xmm9						\n\t		mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		(%%rdi)		,%%xmm12					\n\t		mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"mulpd		(%%rdi)		,%%xmm8						\n\t		mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"addpd		%%xmm6		,%%xmm11					\n\t		addpd		%%xmm2		,%%xmm10		\n\t"\
		"addpd		%%xmm0		,%%xmm9						\n\t		addpd		%%xmm5		,%%xmm15		\n\t"\
		"addpd		%%xmm7		,%%xmm12					\n\t		addpd		%%xmm4		,%%xmm14		\n\t"\
		"addpd		%%xmm1		,%%xmm8						\n\t		addpd		%%xmm3		,%%xmm13		\n\t"\
		"movaps		%%xmm6		,0x080(%%rcx)				\n\t		movaps		%%xmm2		,0x0a0(%%rcx)	\n\t"\
		"movaps		%%xmm0		,0x040(%%rcx)				\n\t		movaps		%%xmm5		,0x060(%%rcx)	\n\t"\
		"movaps		%%xmm7		,0x090(%%rcx)				\n\t		movaps		%%xmm4		,0x0b0(%%rcx)	\n\t"\
		"movaps		%%xmm1		,0x0d0(%%rcx)				\n\t		movaps		%%xmm3		,0x0f0(%%rcx)	\n\t"\
		"movaps		%%xmm11		,     (%%rcx)				\n\t		movaps		%%xmm10		,0x020(%%rcx)	\n\t"\
		"movaps		%%xmm9		,0x0c0(%%rcx)				\n\t		movaps		%%xmm15		,0x0e0(%%rcx)	\n\t"\
		"movaps		%%xmm12		,0x010(%%rcx)				\n\t		movaps		%%xmm14		,0x030(%%rcx)	\n\t"\
		"movaps		%%xmm8		,0x050(%%rcx)				\n\t		movaps		%%xmm13		,0x070(%%rcx)	\n\t"\
		"\n\t"\
		"/*...Block 2:	*/\n\t"\
		"addq		$0x20		,%%rax\n\t"\
		"addq		$0x20		,%%rbx\n\t"\
		"/*****	SSE2_RADIX4_DIF_4WRAPPER(c02,c0A,c12,c1A,r10) ****	SSE2_RADIX4_DIF_4WRAPPER_2NDOFTWO(c06,c0E,c16,c1E,r18)	*****/\n\t"\
		"movq		%[__r10]	,%%rcx						\n\t		/*addq		$0x080,%%rcx // __r18 */	\n\t"\
		"movq		%[__c02]	,%%rdx						\n\t		/*addq		$0x080,%%rdx // __c06 */	\n\t"\
		"movaps		     (%%rax),%%xmm6						\n\t		movaps		0x40 (%%rax),%%xmm14		\n\t"\
		"movaps		     (%%rax),%%xmm0						\n\t		movaps		0x40 (%%rax),%%xmm8			\n\t"\
		"unpckhpd	     (%%rbx),%%xmm6						\n\t		unpckhpd	0x40 (%%rbx),%%xmm14		\n\t"\
		"unpcklpd	     (%%rbx),%%xmm0						\n\t		unpcklpd	0x40 (%%rbx),%%xmm8			\n\t"\
		"movaps		%%xmm6		,0x200(%%rcx)				\n\t		movaps		%%xmm14,0x280(%%rcx)		\n\t"\
		"movaps		0x010(%%rax),%%xmm7						\n\t		movaps		0x50 (%%rax),%%xmm15		\n\t"\
		"movaps		0x010(%%rax),%%xmm1						\n\t		movaps		0x50 (%%rax),%%xmm9			\n\t"\
		"unpckhpd	0x010(%%rbx),%%xmm7						\n\t		unpckhpd	0x50 (%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x010(%%rbx),%%xmm1						\n\t		unpcklpd	0x50 (%%rbx),%%xmm9			\n\t"\
		"movaps		%%xmm7		,0x210(%%rcx)				\n\t		movaps		%%xmm15,0x290(%%rcx)		\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		movaps		%%xmm8,%%xmm10				\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		movaps		%%xmm9,%%xmm11				\n\t"\
		"mulpd		     (%%rdx),%%xmm0						\n\t		mulpd		0x080(%%rdx),%%xmm8			\n\t"\
		"mulpd		     (%%rdx),%%xmm1						\n\t		mulpd		0x080(%%rdx),%%xmm9			\n\t"\
		"mulpd		0x010(%%rdx),%%xmm2						\n\t		mulpd		0x090(%%rdx),%%xmm10		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm3						\n\t		mulpd		0x090(%%rdx),%%xmm11		\n\t"\
		"addpd		%%xmm2		,%%xmm1						\n\t		addpd		%%xmm10,%%xmm9				\n\t"\
		"subpd		%%xmm3		,%%xmm0						\n\t		subpd		%%xmm11,%%xmm8				\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		movaps		%%xmm8,%%xmm10				\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		movaps		%%xmm9,%%xmm11				\n\t"\
		"addq		$0x020		,%%rdx																			\n\t"\
		"movaps		0x100(%%rax),%%xmm6						\n\t		movaps		0x140(%%rax),%%xmm14		\n\t"\
		"movaps		0x100(%%rax),%%xmm4						\n\t		movaps		0x140(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x100(%%rbx),%%xmm6						\n\t		unpckhpd	0x140(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x100(%%rbx),%%xmm4						\n\t		unpcklpd	0x140(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x220(%%rcx)				\n\t		movaps		%%xmm14,0x2a0(%%rcx)		\n\t"\
		"movaps		0x110(%%rax),%%xmm7						\n\t		movaps		0x150(%%rax),%%xmm15		\n\t"\
		"movaps		0x110(%%rax),%%xmm5						\n\t		movaps		0x150(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x110(%%rbx),%%xmm7						\n\t		unpckhpd	0x150(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x110(%%rbx),%%xmm5						\n\t		unpcklpd	0x150(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x230(%%rcx)				\n\t		movaps		%%xmm15,0x2b0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd		     (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd		     (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"addpd		%%xmm4		,%%xmm0						\n\t		addpd		%%xmm12,%%xmm8				\n\t"\
		"addpd		%%xmm5		,%%xmm1						\n\t		addpd		%%xmm13,%%xmm9				\n\t"\
		"subpd		%%xmm4		,%%xmm2						\n\t		subpd		%%xmm12,%%xmm10				\n\t"\
		"subpd		%%xmm5		,%%xmm3						\n\t		subpd		%%xmm13,%%xmm11				\n\t"\
		"addq		$0x040		,%%rdx																			\n\t"\
		"movaps		0x180(%%rax),%%xmm6						\n\t		movaps		0x1c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x180(%%rax),%%xmm4						\n\t		movaps		0x1c0(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x180(%%rbx),%%xmm6						\n\t		unpckhpd	0x1c0(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x180(%%rbx),%%xmm4						\n\t		unpcklpd	0x1c0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x260(%%rcx)				\n\t		movaps		%%xmm14,0x2e0(%%rcx)		\n\t"\
		"movaps		0x190(%%rax),%%xmm7						\n\t		movaps		0x1d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x190(%%rax),%%xmm5						\n\t		movaps		0x1d0(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x190(%%rbx),%%xmm7						\n\t		unpckhpd	0x1d0(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x190(%%rbx),%%xmm5						\n\t		unpcklpd	0x1d0(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x270(%%rcx)				\n\t		movaps		%%xmm15,0x2f0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd		     (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd		     (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"movaps		%%xmm5		,0x010(%%rcx)				\n\t		movaps		%%xmm13,0x090(%%rcx)		\n\t"\
		"movaps		%%xmm4		,     (%%rcx)				\n\t		movaps		%%xmm12,0x080(%%rcx)		\n\t"\
		"subq		$0x020		,%%rdx																			\n\t"\
		"movaps		0x080(%%rax),%%xmm6						\n\t		movaps		0x0c0(%%rax),%%xmm14		\n\t"\
		"movaps		0x080(%%rax),%%xmm4						\n\t		movaps		0x0c0(%%rax),%%xmm12		\n\t"\
		"unpckhpd	0x080(%%rbx),%%xmm6						\n\t		unpckhpd	0x0c0(%%rbx),%%xmm14		\n\t"\
		"unpcklpd	0x080(%%rbx),%%xmm4						\n\t		unpcklpd	0x0c0(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm6		,0x240(%%rcx)				\n\t		movaps		%%xmm14,0x2c0(%%rcx)		\n\t"\
		"movaps		0x090(%%rax),%%xmm7						\n\t		movaps		0x0d0(%%rax),%%xmm15		\n\t"\
		"movaps		0x090(%%rax),%%xmm5						\n\t		movaps		0x0d0(%%rax),%%xmm13		\n\t"\
		"unpckhpd	0x090(%%rbx),%%xmm7						\n\t		unpckhpd	0x0d0(%%rbx),%%xmm15		\n\t"\
		"unpcklpd	0x090(%%rbx),%%xmm5						\n\t		unpcklpd	0x0d0(%%rbx),%%xmm13		\n\t"\
		"movaps		%%xmm7		,0x250(%%rcx)				\n\t		movaps		%%xmm15,0x2d0(%%rcx)		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"mulpd		     (%%rdx),%%xmm4						\n\t		mulpd		0x080(%%rdx),%%xmm12		\n\t"\
		"mulpd		     (%%rdx),%%xmm5						\n\t		mulpd		0x080(%%rdx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm6						\n\t		mulpd		0x090(%%rdx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rdx),%%xmm7						\n\t		mulpd		0x090(%%rdx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14,%%xmm13				\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15,%%xmm12				\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12,%%xmm14				\n\t"\
		"subpd		     (%%rcx),%%xmm4						\n\t		subpd		0x080(%%rcx),%%xmm12		\n\t"\
		"subpd		0x010(%%rcx),%%xmm5						\n\t		subpd		0x090(%%rcx),%%xmm13		\n\t"\
		"addpd		     (%%rcx),%%xmm6						\n\t		addpd		0x080(%%rcx),%%xmm14		\n\t"\
		"addpd		0x010(%%rcx),%%xmm7						\n\t		addpd		0x090(%%rcx),%%xmm15		\n\t"\
		"subpd		%%xmm6		,%%xmm0						\n\t		subpd		%%xmm14,%%xmm8				\n\t"\
		"subpd		%%xmm7		,%%xmm1						\n\t		subpd		%%xmm15,%%xmm9				\n\t"\
		"/*movaps	%%xmm0		,0x040(%%rcx)	*/			\n\t	/*	movaps		%%xmm8,0x0c0(%%rcx)		*/	\n\t"\
		"/*movaps	%%xmm1		,0x050(%%rcx)	*/			\n\t	/*	movaps		%%xmm9,0x0d0(%%rcx)		*/	\n\t"\
		"mulpd		(%%rdi)		,%%xmm6						\n\t		mulpd		(%%rdi),%%xmm14				\n\t"\
		"mulpd		(%%rdi)		,%%xmm7						\n\t		mulpd		(%%rdi),%%xmm15				\n\t"\
		"addpd		%%xmm0		,%%xmm6						\n\t		addpd		%%xmm8,%%xmm14				\n\t"\
		"addpd		%%xmm1		,%%xmm7						\n\t		addpd		%%xmm9,%%xmm15				\n\t"\
		"/*movaps	%%xmm6		,     (%%rcx)	*/			\n\t		movaps		%%xmm14,0x080(%%rcx)		\n\t"\
		"/*movaps	%%xmm7		,0x010(%%rcx)	*/			\n\t		movaps		%%xmm15,0x090(%%rcx)		\n\t"\
		"subpd		%%xmm5		,%%xmm2						\n\t		subpd		%%xmm13,%%xmm10				\n\t"\
		"subpd		%%xmm4		,%%xmm3						\n\t		subpd		%%xmm12,%%xmm11				\n\t"\
		"/*movaps	%%xmm2		,0x020(%%rcx)	*/			\n\t		mulpd		(%%rdi),%%xmm13				\n\t"\
		"/*movaps	%%xmm3		,0x070(%%rcx)	*/			\n\t		mulpd		(%%rdi),%%xmm12				\n\t"\
		"mulpd		(%%rdi)		,%%xmm5						\n\t		addpd		%%xmm10,%%xmm13				\n\t"\
		"mulpd		(%%rdi)		,%%xmm4						\n\t		addpd		%%xmm11,%%xmm12				\n\t"\
		"addpd		%%xmm2		,%%xmm5						\n\t		movaps		%%xmm10,%%xmm14				\n\t"\
		"addpd		%%xmm3		,%%xmm4						\n\t		movaps		%%xmm13,%%xmm15				\n\t"\
		"/*movaps	%%xmm5		,0x060(%%rcx)	*/			\n\t		subpd		%%xmm12,%%xmm10				\n\t"\
		"/*movaps	%%xmm4		,0x030(%%rcx)	*/			\n\t		subpd		%%xmm11,%%xmm13				\n\t"\
		"																addpd		%%xmm12,%%xmm14				\n\t"\
		"																addpd		%%xmm11,%%xmm15				\n\t"\
		"																mulpd		(%%rsi),%%xmm10	/* isrt2 */	\n\t"\
		"																mulpd		(%%rsi),%%xmm13				\n\t"\
		"																mulpd		(%%rsi),%%xmm14				\n\t"\
		"																mulpd		(%%rsi),%%xmm15				\n\t"\
		"															/*	movaps		%%xmm10,0x0a0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm13,0x0e0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm14,0x0b0(%%rcx)	*/	\n\t"\
		"															/*	movaps		%%xmm15,0x0f0(%%rcx)	*/	\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10)	*****/\n\t"\
		"																movaps		0x080(%%rcx),%%xmm11		\n\t"\
		"																movaps		0x090(%%rcx),%%xmm12		\n\t"\
		"subpd		%%xmm11		,%%xmm6						\n\t		subpd		%%xmm10		,%%xmm2			\n\t"\
		"subpd		%%xmm9		,%%xmm0						\n\t		subpd		%%xmm15		,%%xmm5			\n\t"\
		"subpd		%%xmm12		,%%xmm7						\n\t		subpd		%%xmm14		,%%xmm4			\n\t"\
		"subpd		%%xmm8		,%%xmm1						\n\t		subpd		%%xmm13		,%%xmm3			\n\t"\
		"mulpd		(%%rdi)	,%%xmm11						\n\t		mulpd		(%%rdi)		,%%xmm10		\n\t"\
		"mulpd		(%%rdi)	,%%xmm9							\n\t		mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		(%%rdi)	,%%xmm12						\n\t		mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"mulpd		(%%rdi)	,%%xmm8							\n\t		mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"addpd		%%xmm6		,%%xmm11					\n\t		addpd		%%xmm2		,%%xmm10		\n\t"\
		"addpd		%%xmm0		,%%xmm9						\n\t		addpd		%%xmm5		,%%xmm15		\n\t"\
		"addpd		%%xmm7		,%%xmm12					\n\t		addpd		%%xmm4		,%%xmm14		\n\t"\
		"addpd		%%xmm1		,%%xmm8						\n\t		addpd		%%xmm3		,%%xmm13		\n\t"\
		"movaps		%%xmm6		,0x080(%%rcx)				\n\t		movaps		%%xmm2		,0x0a0(%%rcx)	\n\t"\
		"movaps		%%xmm0		,0x040(%%rcx)				\n\t		movaps		%%xmm5		,0x060(%%rcx)	\n\t"\
		"movaps		%%xmm7		,0x090(%%rcx)				\n\t		movaps		%%xmm4		,0x0b0(%%rcx)	\n\t"\
		"movaps		%%xmm1		,0x0d0(%%rcx)				\n\t		movaps		%%xmm3		,0x0f0(%%rcx)	\n\t"\
		"movaps		%%xmm11		,     (%%rcx)				\n\t		movaps		%%xmm10		,0x020(%%rcx)	\n\t"\
		"movaps		%%xmm9		,0x0c0(%%rcx)				\n\t		movaps		%%xmm15		,0x0e0(%%rcx)	\n\t"\
		"movaps		%%xmm12		,0x010(%%rcx)				\n\t		movaps		%%xmm14		,0x030(%%rcx)	\n\t"\
		"movaps		%%xmm8		,0x050(%%rcx)				\n\t		movaps		%%xmm13		,0x070(%%rcx)	\n\t"\
		"/********************************************************************************************************\n\t"\
		" Next 2 blocks operate on odd-indexed elements from the unpck*pd commands which we stored to temporaries:\n\t"\
		"********************************************************************************************************/\n\t"\
		"/*...Block 3:	*/\n\t"\
		"/*	SSE2_RADIX4_DIF_4TWIDDLE         (r20,r24,r22,r26,r20,c01) */\n\t"\
		"addq		$0x100		,%%rcx		/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r28,r2C,r2A,r2E,r28,c05)	*****/\n\t"\
		"movq		%[__c01]	,%%rbx						\n\t	/*	movq		%[__c05]	,%%rbx	*/		\n\t"\
		"movq		%%rcx		,%%rax						\n\t	/*	addq		$0x040		,%%rax	*/		\n\t"\
		"addq		$0x020		,%%rcx						\n\t	/*	addq		$0x040		,%%rcx	*/		\n\t"\
		"movaps		     (%%rax),%%xmm0	\n\t	movq %%rax,%%rdx \n\t	movaps		0x080(%%rax),%%xmm8			\n\t"\
		"movaps		     (%%rcx),%%xmm4						\n\t		movaps		0x080(%%rcx),%%xmm12		\n\t"\
		"movaps		0x010(%%rax),%%xmm1						\n\t		movaps		0x090(%%rax),%%xmm9			\n\t"\
		"movaps		0x010(%%rcx),%%xmm5						\n\t		movaps		0x090(%%rcx),%%xmm13		\n\t"\
		"movaps		     (%%rbx),%%xmm6						\n\t		movaps		0x080(%%rbx),%%xmm14		\n\t"\
		"movaps		0x010(%%rbx),%%xmm7						\n\t		movaps		0x090(%%rbx),%%xmm15		\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		movaps		%%xmm8		,%%xmm10		\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		movaps		%%xmm9		,%%xmm11		\n\t"\
		"mulpd		%%xmm6		,%%xmm0						\n\t		mulpd		%%xmm14		,%%xmm8			\n\t"\
		"mulpd		%%xmm6		,%%xmm1						\n\t		mulpd		%%xmm14		,%%xmm9			\n\t"\
		"mulpd		%%xmm7		,%%xmm2						\n\t		mulpd		%%xmm15		,%%xmm10		\n\t"\
		"mulpd		%%xmm7		,%%xmm3						\n\t		mulpd		%%xmm15		,%%xmm11		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12		,%%xmm14		\n\t"\
		"addpd		%%xmm2		,%%xmm1						\n\t		addpd		%%xmm10		,%%xmm9			\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13		,%%xmm15		\n\t"\
		"mulpd		0x20 (%%rbx),%%xmm4						\n\t		mulpd		0xa0 (%%rbx),%%xmm12		\n\t"\
		"subpd		%%xmm3		,%%xmm0						\n\t		subpd		%%xmm11		,%%xmm8			\n\t"\
		"mulpd		0x20 (%%rbx),%%xmm5						\n\t		mulpd		0xa0 (%%rbx),%%xmm13		\n\t"\
		"mulpd		0x30 (%%rbx),%%xmm6						\n\t		mulpd		0xb0 (%%rbx),%%xmm14		\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		movaps		%%xmm8		,%%xmm10		\n\t"\
		"mulpd		0x30 (%%rbx),%%xmm7						\n\t		mulpd		0xb0 (%%rbx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14		,%%xmm13		\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		movaps		%%xmm9		,%%xmm11		\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15		,%%xmm12		\n\t"\
		"addq		$0x040		,%%rcx						\n\t		addpd		%%xmm12		,%%xmm8			\n\t"\
		"addq		$0x60		,%%rbx						\n\t		addpd		%%xmm13		,%%xmm9			\n\t"\
		"movaps		     (%%rcx),%%xmm6						\n\t		subpd		%%xmm12		,%%xmm10		\n\t"\
		"movaps		0x010(%%rcx),%%xmm7						\n\t		subpd		%%xmm13		,%%xmm11		\n\t"\
		"addpd		%%xmm4		,%%xmm0						\n\t		movaps		0x080(%%rcx),%%xmm12		\n\t"\
		"addpd		%%xmm5		,%%xmm1						\n\t		movaps		0x090(%%rcx),%%xmm13		\n\t"\
		"subpd		%%xmm4		,%%xmm2						\n\t		movaps		0x080(%%rcx),%%xmm14		\n\t"\
		"subpd		%%xmm5		,%%xmm3						\n\t		movaps		0x090(%%rcx),%%xmm15		\n\t"\
		"movaps		%%xmm6		,%%xmm4						\n\t		mulpd		0x080(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm7		,%%xmm5						\n\t		mulpd		0x080(%%rbx),%%xmm13		\n\t"\
		"mulpd		     (%%rbx),%%xmm4						\n\t		mulpd		0x090(%%rbx),%%xmm14		\n\t"\
		"mulpd		     (%%rbx),%%xmm5						\n\t		mulpd		0x090(%%rbx),%%xmm15		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm6						\n\t		addpd		%%xmm14		,%%xmm13		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm7						\n\t		subpd		%%xmm15		,%%xmm12		\n\t"\
		"/*movq		%%rax		,%%rdx		*/				\n\t		movaps		%%xmm13		,0x090(%%rdx)	\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		movaps		%%xmm12		,0x080(%%rdx)	\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t	/*	subq	$0x20	,%%rbx */				\n\t"\
		"movaps		%%xmm5		,0x010(%%rdx)	\n\t	addq	$0x040		,%%rax								\n\t"\
		"movaps		%%xmm4		,     (%%rdx)	\n\t	subq	$0x20		,%%rbx								\n\t"\
		"movaps		     (%%rax),%%xmm4						\n\t		movaps		0x080(%%rax),%%xmm12		\n\t"\
		"movaps		0x010(%%rax),%%xmm5						\n\t		movaps		0x090(%%rax),%%xmm13		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		0x080(%%rax),%%xmm14		\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		0x090(%%rax),%%xmm15		\n\t"\
		"mulpd		     (%%rbx),%%xmm4						\n\t		mulpd		0x080(%%rbx),%%xmm12		\n\t"\
		"mulpd		     (%%rbx),%%xmm5						\n\t		mulpd		0x080(%%rbx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm6						\n\t		mulpd		0x090(%%rbx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm7						\n\t		mulpd		0x090(%%rbx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14		,%%xmm13		\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15		,%%xmm12		\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13		,%%xmm15		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12		,%%xmm14		\n\t"\
		"subpd		     (%%rdx),%%xmm4						\n\t		subpd		0x080(%%rdx),%%xmm12		\n\t"\
		"subpd		0x010(%%rdx),%%xmm5						\n\t		subpd		0x090(%%rdx),%%xmm13		\n\t"\
		"addpd		     (%%rdx),%%xmm6						\n\t		addpd		0x080(%%rdx),%%xmm14		\n\t"\
		"addpd		0x010(%%rdx),%%xmm7						\n\t		addpd		0x090(%%rdx),%%xmm15		\n\t"\
		"subpd		%%xmm6		,%%xmm0						\n\t		subpd		%%xmm14		,%%xmm8			\n\t"\
		"subpd		%%xmm5		,%%xmm2						\n\t		subpd		%%xmm15		,%%xmm9			\n\t"\
		"subpd		%%xmm7		,%%xmm1						\n\t	/*	movaps		%%xmm8		,0x0c0(%%rdx)*/	\n\t"\
		"subpd		%%xmm4		,%%xmm3						\n\t		subpd		%%xmm13		,%%xmm10		\n\t"\
		"/*movaps		%%xmm0		,0x040(%%rdx)	*/		\n\t	/*	movaps		%%xmm9		,0x0d0(%%rdx)*/	\n\t"\
		"/*movaps		%%xmm2		,0x020(%%rdx)	*/		\n\t		subpd		%%xmm12		,%%xmm11		\n\t"\
		"/*movaps		%%xmm1		,0x050(%%rdx)	*/		\n\t		mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"/*movaps		%%xmm3		,0x070(%%rdx)	*/		\n\t		mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"mulpd		(%%rdi)		,%%xmm6						\n\t		mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		(%%rdi)		,%%xmm5						\n\t		mulpd		(%%rdi)		,%%xmm12		\n\t"\
		"mulpd		(%%rdi)		,%%xmm7						\n\t		addpd		%%xmm8		,%%xmm14		\n\t"\
		"mulpd		(%%rdi)		,%%xmm4						\n\t		addpd		%%xmm10		,%%xmm13		\n\t"\
		"addpd		%%xmm0		,%%xmm6						\n\t		addpd		%%xmm9		,%%xmm15		\n\t"\
		"addpd		%%xmm2		,%%xmm5						\n\t		addpd		%%xmm11		,%%xmm12		\n\t"\
		"addpd		%%xmm1		,%%xmm7						\n\t		movaps		%%xmm14		,0x080(%%rdx)	\n\t"\
		"addpd		%%xmm3		,%%xmm4						\n\t		movaps		%%xmm15		,0x090(%%rdx)	\n\t"\
		"/*movaps		%%xmm6		,     (%%rdx)	*/		\n\t		movaps		%%xmm10		,%%xmm14		\n\t"\
		"/*movaps		%%xmm5		,0x060(%%rdx)	*/		\n\t		movaps		%%xmm13		,%%xmm15		\n\t"\
		"/*movaps		%%xmm7		,0x010(%%rdx)	*/		\n\t		subpd		%%xmm12		,%%xmm10		\n\t"\
		"/*movaps		%%xmm4		,0x030(%%rdx)	*/		\n\t		subpd		%%xmm11		,%%xmm13		\n\t"\
		"													\n\t		addpd		%%xmm12		,%%xmm14		\n\t"\
		"													\n\t		addpd		%%xmm11		,%%xmm15		\n\t"\
		"													\n\t		mulpd		(%%rsi)		,%%xmm10		\n\t"\
		"													\n\t		mulpd		(%%rsi)		,%%xmm13		\n\t"\
		"													\n\t		mulpd		(%%rsi)		,%%xmm14		\n\t"\
		"													\n\t		mulpd		(%%rsi)		,%%xmm15		\n\t"\
		"													\n\t	/*	movaps		%%xmm10		,0x0a0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm13		,0x0e0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm14		,0x0b0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm15		,0x0f0(%%rdx)*/	\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20)	*****/\n\t"\
		"																movaps		0x080(%%rdx),%%xmm11		\n\t"\
		"																movaps		0x090(%%rdx),%%xmm12		\n\t"\
		"subpd		%%xmm11		,%%xmm6						\n\t		subpd		%%xmm10		,%%xmm2			\n\t"\
		"subpd		%%xmm9		,%%xmm0						\n\t		subpd		%%xmm15		,%%xmm5			\n\t"\
		"subpd		%%xmm12		,%%xmm7						\n\t		subpd		%%xmm14		,%%xmm4			\n\t"\
		"subpd		%%xmm8		,%%xmm1						\n\t		subpd		%%xmm13		,%%xmm3			\n\t"\
		"mulpd		(%%rdi)		,%%xmm11					\n\t		mulpd		(%%rdi)		,%%xmm10		\n\t"\
		"mulpd		(%%rdi)		,%%xmm9						\n\t		mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		(%%rdi)		,%%xmm12					\n\t		mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"mulpd		(%%rdi)		,%%xmm8						\n\t		mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"addpd		%%xmm6		,%%xmm11					\n\t		addpd		%%xmm2		,%%xmm10		\n\t"\
		"addpd		%%xmm0		,%%xmm9						\n\t		addpd		%%xmm5		,%%xmm15		\n\t"\
		"addpd		%%xmm7		,%%xmm12					\n\t		addpd		%%xmm4		,%%xmm14		\n\t"\
		"addpd		%%xmm1		,%%xmm8						\n\t		addpd		%%xmm3		,%%xmm13		\n\t"\
		"movaps		%%xmm6		,0x080(%%rdx)				\n\t		movaps		%%xmm2		,0x0a0(%%rdx)	\n\t"\
		"movaps		%%xmm0		,0x040(%%rdx)				\n\t		movaps		%%xmm5		,0x060(%%rdx)	\n\t"\
		"movaps		%%xmm7		,0x090(%%rdx)				\n\t		movaps		%%xmm4		,0x0b0(%%rdx)	\n\t"\
		"movaps		%%xmm1		,0x0d0(%%rdx)				\n\t		movaps		%%xmm3		,0x0f0(%%rdx)	\n\t"\
		"movaps		%%xmm11		,     (%%rdx)				\n\t		movaps		%%xmm10		,0x020(%%rdx)	\n\t"\
		"movaps		%%xmm9		,0x0c0(%%rdx)				\n\t		movaps		%%xmm15		,0x0e0(%%rdx)	\n\t"\
		"movaps		%%xmm12		,0x010(%%rdx)				\n\t		movaps		%%xmm14		,0x030(%%rdx)	\n\t"\
		"movaps		%%xmm8		,0x050(%%rdx)				\n\t		movaps		%%xmm13		,0x070(%%rdx)	\n\t"\
		"/*...Block 4:	*/\n\t"\
		"/*****	SSE2_RADIX4_DIF_4TWIDDLE         (r30,r34,r32,r36,r30,c03)	*****/\n\t"\
		"movq		%[__c03]	,%%rbx					\n\t"\
		"movq		%[__r30]	,%%rax					\n\t"\
		"movq		%%rax		,%%rcx		/*****	SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(r38,r3C,r3A,r3E,r38,c07)	*****/\n\t"\
		"addq		$0x020		,%%rcx					\n\t"\
		"movaps		     (%%rax),%%xmm0	\n\t	movq %%rax,%%rdx \n\t	movaps		0x080(%%rax),%%xmm8			\n\t"\
		"movaps		     (%%rcx),%%xmm4						\n\t		movaps		0x080(%%rcx),%%xmm12		\n\t"\
		"movaps		0x010(%%rax),%%xmm1						\n\t		movaps		0x090(%%rax),%%xmm9			\n\t"\
		"movaps		0x010(%%rcx),%%xmm5						\n\t		movaps		0x090(%%rcx),%%xmm13		\n\t"\
		"movaps		     (%%rbx),%%xmm6						\n\t		movaps		0x080(%%rbx),%%xmm14		\n\t"\
		"movaps		0x010(%%rbx),%%xmm7						\n\t		movaps		0x090(%%rbx),%%xmm15		\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		movaps		%%xmm8		,%%xmm10		\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		movaps		%%xmm9		,%%xmm11		\n\t"\
		"mulpd		%%xmm6		,%%xmm0						\n\t		mulpd		%%xmm14		,%%xmm8			\n\t"\
		"mulpd		%%xmm6		,%%xmm1						\n\t		mulpd		%%xmm14		,%%xmm9			\n\t"\
		"mulpd		%%xmm7		,%%xmm2						\n\t		mulpd		%%xmm15		,%%xmm10		\n\t"\
		"mulpd		%%xmm7		,%%xmm3						\n\t		mulpd		%%xmm15		,%%xmm11		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12		,%%xmm14		\n\t"\
		"addpd		%%xmm2		,%%xmm1						\n\t		addpd		%%xmm10		,%%xmm9			\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13		,%%xmm15		\n\t"\
		"mulpd		0x20 (%%rbx),%%xmm4						\n\t		mulpd		0xa0 (%%rbx),%%xmm12		\n\t"\
		"subpd		%%xmm3		,%%xmm0						\n\t		subpd		%%xmm11		,%%xmm8			\n\t"\
		"mulpd		0x20 (%%rbx),%%xmm5						\n\t		mulpd		0xa0 (%%rbx),%%xmm13		\n\t"\
		"mulpd		0x30 (%%rbx),%%xmm6						\n\t		mulpd		0xb0 (%%rbx),%%xmm14		\n\t"\
		"movaps		%%xmm0		,%%xmm2						\n\t		movaps		%%xmm8		,%%xmm10		\n\t"\
		"mulpd		0x30 (%%rbx),%%xmm7						\n\t		mulpd		0xb0 (%%rbx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14		,%%xmm13		\n\t"\
		"movaps		%%xmm1		,%%xmm3						\n\t		movaps		%%xmm9		,%%xmm11		\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15		,%%xmm12		\n\t"\
		"addq		$0x040		,%%rcx						\n\t		addpd		%%xmm12		,%%xmm8			\n\t"\
		"addq		$0x60		,%%rbx						\n\t		addpd		%%xmm13		,%%xmm9			\n\t"\
		"movaps		     (%%rcx),%%xmm6						\n\t		subpd		%%xmm12		,%%xmm10		\n\t"\
		"movaps		0x010(%%rcx),%%xmm7						\n\t		subpd		%%xmm13		,%%xmm11		\n\t"\
		"addpd		%%xmm4		,%%xmm0						\n\t		movaps		0x080(%%rcx),%%xmm12		\n\t"\
		"addpd		%%xmm5		,%%xmm1						\n\t		movaps		0x090(%%rcx),%%xmm13		\n\t"\
		"subpd		%%xmm4		,%%xmm2						\n\t		movaps		0x080(%%rcx),%%xmm14		\n\t"\
		"subpd		%%xmm5		,%%xmm3						\n\t		movaps		0x090(%%rcx),%%xmm15		\n\t"\
		"movaps		%%xmm6		,%%xmm4						\n\t		mulpd		0x080(%%rbx),%%xmm12		\n\t"\
		"movaps		%%xmm7		,%%xmm5						\n\t		mulpd		0x080(%%rbx),%%xmm13		\n\t"\
		"mulpd		     (%%rbx),%%xmm4						\n\t		mulpd		0x090(%%rbx),%%xmm14		\n\t"\
		"mulpd		     (%%rbx),%%xmm5						\n\t		mulpd		0x090(%%rbx),%%xmm15		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm6						\n\t		addpd		%%xmm14		,%%xmm13		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm7						\n\t		subpd		%%xmm15		,%%xmm12		\n\t"\
		"/*movq		%%rax		,%%rdx		*/				\n\t		movaps		%%xmm13		,0x090(%%rdx)	\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		movaps		%%xmm12		,0x080(%%rdx)	\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t	/*	subq	$0x20	,%%rbx */				\n\t"\
		"movaps		%%xmm5		,0x010(%%rdx)	\n\t	addq	$0x040		,%%rax								\n\t"\
		"movaps		%%xmm4		,     (%%rdx)	\n\t	subq	$0x20		,%%rbx								\n\t"\
		"movaps		     (%%rax),%%xmm4						\n\t		movaps		0x080(%%rax),%%xmm12		\n\t"\
		"movaps		0x010(%%rax),%%xmm5						\n\t		movaps		0x090(%%rax),%%xmm13		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		0x080(%%rax),%%xmm14		\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		0x090(%%rax),%%xmm15		\n\t"\
		"mulpd		     (%%rbx),%%xmm4						\n\t		mulpd		0x080(%%rbx),%%xmm12		\n\t"\
		"mulpd		     (%%rbx),%%xmm5						\n\t		mulpd		0x080(%%rbx),%%xmm13		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm6						\n\t		mulpd		0x090(%%rbx),%%xmm14		\n\t"\
		"mulpd		0x010(%%rbx),%%xmm7						\n\t		mulpd		0x090(%%rbx),%%xmm15		\n\t"\
		"addpd		%%xmm6		,%%xmm5						\n\t		addpd		%%xmm14		,%%xmm13		\n\t"\
		"subpd		%%xmm7		,%%xmm4						\n\t		subpd		%%xmm15		,%%xmm12		\n\t"\
		"movaps		%%xmm5		,%%xmm7						\n\t		movaps		%%xmm13		,%%xmm15		\n\t"\
		"movaps		%%xmm4		,%%xmm6						\n\t		movaps		%%xmm12		,%%xmm14		\n\t"\
		"subpd		     (%%rdx),%%xmm4						\n\t		subpd		0x080(%%rdx),%%xmm12		\n\t"\
		"subpd		0x010(%%rdx),%%xmm5						\n\t		subpd		0x090(%%rdx),%%xmm13		\n\t"\
		"addpd		     (%%rdx),%%xmm6						\n\t		addpd		0x080(%%rdx),%%xmm14		\n\t"\
		"addpd		0x010(%%rdx),%%xmm7						\n\t		addpd		0x090(%%rdx),%%xmm15		\n\t"\
		"subpd		%%xmm6		,%%xmm0						\n\t		subpd		%%xmm14		,%%xmm8			\n\t"\
		"subpd		%%xmm5		,%%xmm2						\n\t		subpd		%%xmm15		,%%xmm9			\n\t"\
		"subpd		%%xmm7		,%%xmm1						\n\t	/*	movaps		%%xmm8		,0x0c0(%%rdx)*/	\n\t"\
		"subpd		%%xmm4		,%%xmm3						\n\t		subpd		%%xmm13		,%%xmm10		\n\t"\
		"/*movaps		%%xmm0		,0x040(%%rdx)	*/		\n\t	/*	movaps		%%xmm9		,0x0d0(%%rdx)*/	\n\t"\
		"/*movaps		%%xmm2		,0x020(%%rdx)	*/		\n\t		subpd		%%xmm12		,%%xmm11		\n\t"\
		"/*movaps		%%xmm1		,0x050(%%rdx)	*/		\n\t		mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"/*movaps		%%xmm3		,0x070(%%rdx)	*/		\n\t		mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"mulpd		(%%rdi)		,%%xmm6						\n\t		mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		(%%rdi)		,%%xmm5						\n\t		mulpd		(%%rdi)		,%%xmm12		\n\t"\
		"mulpd		(%%rdi)		,%%xmm7						\n\t		addpd		%%xmm8		,%%xmm14		\n\t"\
		"mulpd		(%%rdi)		,%%xmm4						\n\t		addpd		%%xmm10		,%%xmm13		\n\t"\
		"addpd		%%xmm0		,%%xmm6						\n\t		addpd		%%xmm9		,%%xmm15		\n\t"\
		"addpd		%%xmm2		,%%xmm5						\n\t		addpd		%%xmm11		,%%xmm12		\n\t"\
		"addpd		%%xmm1		,%%xmm7						\n\t		movaps		%%xmm14		,0x080(%%rdx)	\n\t"\
		"addpd		%%xmm3		,%%xmm4						\n\t		movaps		%%xmm15		,0x090(%%rdx)	\n\t"\
		"/*movaps		%%xmm6		,     (%%rdx)	*/		\n\t		movaps		%%xmm10		,%%xmm14		\n\t"\
		"/*movaps		%%xmm5		,0x060(%%rdx)	*/		\n\t		movaps		%%xmm13		,%%xmm15		\n\t"\
		"/*movaps		%%xmm7		,0x010(%%rdx)	*/		\n\t		subpd		%%xmm12		,%%xmm10		\n\t"\
		"/*movaps		%%xmm4		,0x030(%%rdx)	*/		\n\t		subpd		%%xmm11		,%%xmm13		\n\t"\
		"													\n\t		addpd		%%xmm12		,%%xmm14		\n\t"\
		"													\n\t		addpd		%%xmm11		,%%xmm15		\n\t"\
		"													\n\t		mulpd		(%%rsi)		,%%xmm10		\n\t"\
		"													\n\t		mulpd		(%%rsi)		,%%xmm13		\n\t"\
		"													\n\t		mulpd		(%%rsi)		,%%xmm14		\n\t"\
		"													\n\t		mulpd		(%%rsi)		,%%xmm15		\n\t"\
		"													\n\t	/*	movaps		%%xmm10		,0x0a0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm13		,0x0e0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm14		,0x0b0(%%rdx)*/	\n\t"\
		"													\n\t	/*	movaps		%%xmm15		,0x0f0(%%rdx)*/	\n\t"\
		"/*****	SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20)	*****/\n\t"\
		"																movaps		0x080(%%rdx),%%xmm11		\n\t"\
		"																movaps		0x090(%%rdx),%%xmm12		\n\t"\
		"subpd		%%xmm11		,%%xmm6						\n\t		subpd		%%xmm10		,%%xmm2			\n\t"\
		"subpd		%%xmm9		,%%xmm0						\n\t		subpd		%%xmm15		,%%xmm5			\n\t"\
		"subpd		%%xmm12		,%%xmm7						\n\t		subpd		%%xmm14		,%%xmm4			\n\t"\
		"subpd		%%xmm8		,%%xmm1						\n\t		subpd		%%xmm13		,%%xmm3			\n\t"\
		"mulpd		(%%rdi)		,%%xmm11					\n\t		mulpd		(%%rdi)		,%%xmm10		\n\t"\
		"mulpd		(%%rdi)		,%%xmm9						\n\t		mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		(%%rdi)		,%%xmm12					\n\t		mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"mulpd		(%%rdi)		,%%xmm8						\n\t		mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"addpd		%%xmm6		,%%xmm11					\n\t		addpd		%%xmm2		,%%xmm10		\n\t"\
		"addpd		%%xmm0		,%%xmm9						\n\t		addpd		%%xmm5		,%%xmm15		\n\t"\
		"addpd		%%xmm7		,%%xmm12					\n\t		addpd		%%xmm4		,%%xmm14		\n\t"\
		"addpd		%%xmm1		,%%xmm8						\n\t		addpd		%%xmm3		,%%xmm13		\n\t"\
		"movaps		%%xmm6		,0x080(%%rdx)				\n\t		movaps		%%xmm2		,0x0a0(%%rdx)	\n\t"\
		"movaps		%%xmm0		,0x040(%%rdx)				\n\t		movaps		%%xmm5		,0x060(%%rdx)	\n\t"\
		"movaps		%%xmm7		,0x090(%%rdx)				\n\t		movaps		%%xmm4		,0x0b0(%%rdx)	\n\t"\
		"movaps		%%xmm1		,0x0d0(%%rdx)				\n\t		movaps		%%xmm3		,0x0f0(%%rdx)	\n\t"\
		"movaps		%%xmm11		,     (%%rdx)				\n\t		movaps		%%xmm10		,0x020(%%rdx)	\n\t"\
		"movaps		%%xmm9		,0x0c0(%%rdx)				\n\t		movaps		%%xmm15		,0x0e0(%%rdx)	\n\t"\
		"movaps		%%xmm12		,0x010(%%rdx)				\n\t		movaps		%%xmm14		,0x030(%%rdx)	\n\t"\
		"movaps		%%xmm8		,0x050(%%rdx)				\n\t		movaps		%%xmm13		,0x070(%%rdx)	\n\t"\
		"/**********************************************************************************/\n\t"\
		"/*...and now do eight radix-4 transforms, including the internal twiddle factots: */\n\t"\
		"/**********************************************************************************/\n\t"\
		"movq		%[__isrt2]		,%%rsi		\n\t"\
		"/*...Block 1: t00,t10,t20,t30	*/			\n\t	/*...Block 5: t08,t18,t28,t38	*/		\n\t"\
		"movq		%[__r00]		,%%rax			\n\t	movq		$0x080		,%%r10			\n\t"\
		"movq		%[__r10]		,%%rbx			\n\t	movq		$0x080		,%%r11			\n\t"\
		"movq		%[__r20]		,%%rcx			\n\t	movq		$0x080		,%%r12			\n\t"\
		"movq		%[__r30]		,%%rdx			\n\t	movq		$0x080		,%%r13			\n\t"\
		"movaps		      (%%rax)	,%%xmm0			\n\t	addq		%%rax		,%%r10			\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1			\n\t	addq		%%rbx		,%%r11			\n\t"\
		"movaps		      (%%rbx)	,%%xmm2			\n\t	addq		%%rcx		,%%r12			\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm3			\n\t	addq		%%rdx		,%%r13			\n\t"\
		"addpd		       %%xmm0	,%%xmm2			\n\t	movaps		      (%%r10)	,%%xmm8		\n\t"\
		"addpd		       %%xmm1	,%%xmm3			\n\t	movaps		 0x010(%%r10)	,%%xmm9		\n\t"\
		"subpd		      (%%rbx)	,%%xmm0			\n\t	movaps		      (%%r11)	,%%xmm10	\n\t"\
		"subpd		 0x010(%%rbx)	,%%xmm1			\n\t	movaps		 0x010(%%r11)	,%%xmm11	\n\t"\
		"movaps		      (%%rcx)	,%%xmm4			\n\t	addpd		       %%xmm9	,%%xmm10	\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5			\n\t	addpd		       %%xmm8	,%%xmm11	\n\t"\
		"movaps		      (%%rdx)	,%%xmm6			\n\t	subpd		 0x010(%%r11)	,%%xmm8		\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm7			\n\t	subpd		      (%%r11)	,%%xmm9		\n\t"\
		"addpd		       %%xmm4	,%%xmm6			\n\t	movaps		      (%%r12)	,%%xmm12	\n\t"\
		"addpd		       %%xmm5	,%%xmm7			\n\t	movaps		 0x010(%%r12)	,%%xmm13	\n\t"\
		"subpd		      (%%rdx)	,%%xmm4			\n\t	movaps		      (%%r13)	,%%xmm14	\n\t"\
		"subpd		 0x010(%%rdx)	,%%xmm5			\n\t	movaps		 0x010(%%r13)	,%%xmm15	\n\t"\
		"subpd		%%xmm6		,%%xmm2				\n\t	subpd		       %%xmm13	,%%xmm12	\n\t"\
		"subpd		%%xmm7		,%%xmm3				\n\t	addpd		      (%%r12)	,%%xmm13	\n\t"\
		"mulpd		(%%rdi)		,%%xmm6				\n\t	addpd		       %%xmm15	,%%xmm14	\n\t"\
		"mulpd		(%%rdi)		,%%xmm7				\n\t	subpd		      (%%r13)	,%%xmm15	\n\t"\
		"movaps		%%xmm2		,      (%%rcx)		\n\t	mulpd		(%%rsi)		,%%xmm12		\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)		\n\t	mulpd		(%%rsi)		,%%xmm13		\n\t"\
		"addpd		%%xmm2		,%%xmm6				\n\t	mulpd		(%%rsi)		,%%xmm14		\n\t"\
		"addpd		%%xmm3		,%%xmm7				\n\t	mulpd		(%%rsi)		,%%xmm15		\n\t"\
		"movaps		%%xmm6		,      (%%rax)		\n\t	subpd		%%xmm14		,%%xmm12		\n\t"\
		"movaps		%%xmm7		, 0x010(%%rax)		\n\t	subpd		%%xmm15		,%%xmm13		\n\t"\
		"subpd		%%xmm5		,%%xmm0				\n\t	mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"subpd		%%xmm4		,%%xmm1				\n\t	mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		(%%rdi)		,%%xmm5				\n\t	addpd		%%xmm12		,%%xmm14		\n\t"\
		"mulpd		(%%rdi)		,%%xmm4				\n\t	addpd		%%xmm13		,%%xmm15		\n\t"\
		"movaps		%%xmm0		,      (%%rbx)		\n\t	subpd		%%xmm12		,%%xmm8			\n\t"\
		"movaps		%%xmm1		, 0x010(%%rdx)		\n\t	subpd		%%xmm13		,%%xmm10		\n\t"\
		"addpd		%%xmm0		,%%xmm5				\n\t	mulpd		(%%rdi)		,%%xmm12		\n\t"\
		"addpd		%%xmm1		,%%xmm4				\n\t	mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"movaps		%%xmm5		,      (%%rdx)		\n\t	movaps		%%xmm8		,      (%%r12)	\n\t"\
		"movaps		%%xmm4		, 0x010(%%rbx)		\n\t	movaps		%%xmm10		, 0x010(%%r12)	\n\t"\
		"/*...Block 3: t04,t14,t24,t34	*/			\n\t"\
		"addq		$0x040		,%%rax				\n\t	addpd		%%xmm8		,%%xmm12		\n\t"\
		"addq		$0x040		,%%rbx				\n\t	addpd		%%xmm10		,%%xmm13		\n\t"\
		"addq		$0x040		,%%rcx				\n\t	movaps		%%xmm12		,      (%%r10)	\n\t"\
		"addq		$0x040		,%%rdx				\n\t	movaps		%%xmm13		, 0x010(%%r10)	\n\t"\
		"movaps		0x10(%%rsi)	,%%xmm8		/* c */	\n\t	subpd		%%xmm15		,%%xmm11		\n\t"\
		"movaps		0x20(%%rsi)	,%%xmm10	/* s */	\n\t	subpd		%%xmm14		,%%xmm9			\n\t"\
		"movaps		      (%%rcx)	,%%xmm4			\n\t	mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"movaps		      (%%rdx)	,%%xmm6			\n\t	mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5			\n\t	movaps		%%xmm11		,      (%%r11)	\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm7			\n\t	movaps		%%xmm9		, 0x010(%%r13)	\n\t"\
		"movaps		       %%xmm4	,%%xmm0			\n\t	addpd		%%xmm11		,%%xmm15		\n\t"\
		"movaps		       %%xmm6	,%%xmm2			\n\t	addpd		%%xmm9		,%%xmm14		\n\t"\
		"movaps		       %%xmm5	,%%xmm1			\n\t	movaps		%%xmm15		,      (%%r13)	\n\t"\
		"movaps		       %%xmm7	,%%xmm3			\n\t	movaps		%%xmm14		, 0x010(%%r11)	\n\t"\
		"													/*...Block 7: t0C,t1C,t2C,t3C	*/		\n\t"\
		"mulpd		 %%xmm8	,%%xmm4					\n\t	addq		$0x040		,%%r10			\n\t"\
		"mulpd		 %%xmm10,%%xmm6					\n\t	addq		$0x040		,%%r11			\n\t"\
		"mulpd		 %%xmm10,%%xmm1					\n\t	addq		$0x040		,%%r12			\n\t"\
		"mulpd		 %%xmm8	,%%xmm3					\n\t	addq		$0x040		,%%r13			\n\t"\
		"mulpd		 %%xmm8	,%%xmm5					\n\t	movaps		      (%%r12)	,%%xmm12	\n\t"\
		"mulpd		 %%xmm10,%%xmm7					\n\t	movaps		      (%%r13)	,%%xmm14	\n\t"\
		"mulpd		 %%xmm10,%%xmm0					\n\t	movaps		 0x010(%%r12)	,%%xmm13	\n\t"\
		"mulpd		 %%xmm8	,%%xmm2					\n\t	movaps		 0x010(%%r13)	,%%xmm15	\n\t"\
		"subpd		%%xmm1		,%%xmm4				\n\t	movaps		       %%xmm13	,%%xmm9		\n\t"\
		"subpd		%%xmm3		,%%xmm6				\n\t	movaps		       %%xmm15	,%%xmm11	\n\t"\
		"addpd		%%xmm0		,%%xmm5				\n\t	mulpd		 %%xmm10,%%xmm12			\n\t"\
		"addpd		%%xmm2		,%%xmm7				\n\t	mulpd		 %%xmm8	,%%xmm14			\n\t"\
		"subpd		%%xmm6		,%%xmm4				\n\t	mulpd		 %%xmm8	,%%xmm9				\n\t"\
		"subpd		%%xmm7		,%%xmm5				\n\t	mulpd		 %%xmm10,%%xmm11			\n\t"\
		"mulpd		(%%rdi)		,%%xmm6				\n\t	mulpd		 %%xmm10,%%xmm13			\n\t"\
		"mulpd		(%%rdi)		,%%xmm7				\n\t	mulpd		 %%xmm8	,%%xmm15			\n\t"\
		"addpd		%%xmm4		,%%xmm6				\n\t	mulpd		(%%r12)	,%%xmm8				\n\t"\
		"addpd		%%xmm5		,%%xmm7				\n\t	mulpd		(%%r13)	,%%xmm10			\n\t"\
		"movaps		      (%%rbx)	,%%xmm2			\n\t	subpd		%%xmm9		,%%xmm12		\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm3			\n\t	subpd		%%xmm11		,%%xmm14		\n\t"\
		"subpd		 0x010(%%rbx)	,%%xmm2			\n\t	addpd		%%xmm8		,%%xmm13		\n\t"\
		"addpd		      (%%rbx)	,%%xmm3			\n\t	addpd		%%xmm10		,%%xmm15		\n\t"\
		"mulpd		      (%%rsi)	,%%xmm2			\n\t	subpd		%%xmm14		,%%xmm12		\n\t"\
		"mulpd		      (%%rsi)	,%%xmm3			\n\t	subpd		%%xmm15		,%%xmm13		\n\t"\
		"movaps		      (%%rax)	,%%xmm0			\n\t	mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1			\n\t	mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"subpd		%%xmm2		,%%xmm0				\n\t	addpd		%%xmm12		,%%xmm14		\n\t"\
		"subpd		%%xmm3		,%%xmm1				\n\t	addpd		%%xmm13		,%%xmm15		\n\t"\
		"addpd		      (%%rax)	,%%xmm2			\n\t	movaps		      (%%r11)	,%%xmm10	\n\t"\
		"addpd		 0x010(%%rax)	,%%xmm3			\n\t	movaps		 0x010(%%r11)	,%%xmm11	\n\t"\
		"subpd		%%xmm6		,%%xmm2				\n\t	addpd		 0x010(%%r11)	,%%xmm10	\n\t"\
		"subpd		%%xmm7		,%%xmm3				\n\t	subpd		      (%%r11)	,%%xmm11	\n\t"\
		"mulpd		(%%rdi)		,%%xmm6				\n\t	mulpd		      (%%rsi)	,%%xmm10	\n\t"\
		"mulpd		(%%rdi)		,%%xmm7				\n\t	mulpd		      (%%rsi)	,%%xmm11	\n\t"\
		"movaps		%%xmm2		,      (%%rcx)		\n\t	movaps		      (%%r10)	,%%xmm8		\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)		\n\t	movaps		 0x010(%%r10)	,%%xmm9		\n\t"\
		"addpd		%%xmm2		,%%xmm6				\n\t	subpd		%%xmm10		,%%xmm8			\n\t"\
		"addpd		%%xmm3		,%%xmm7				\n\t	subpd		%%xmm11		,%%xmm9			\n\t"\
		"movaps		%%xmm6		,      (%%rax)		\n\t	addpd		      (%%r10)	,%%xmm10	\n\t"\
		"movaps		%%xmm7		, 0x010(%%rax)		\n\t	addpd		 0x010(%%r10)	,%%xmm11	\n\t"\
		"subpd		%%xmm5		,%%xmm0				\n\t	subpd		%%xmm12		,%%xmm8			\n\t"\
		"subpd		%%xmm4		,%%xmm1				\n\t	subpd		%%xmm13		,%%xmm9			\n\t"\
		"mulpd		(%%rdi)		,%%xmm5				\n\t	mulpd		(%%rdi)		,%%xmm12		\n\t"\
		"mulpd		(%%rdi)		,%%xmm4				\n\t	mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"movaps		%%xmm0		,      (%%rbx)		\n\t	movaps		%%xmm8		,      (%%r12)	\n\t"\
		"movaps		%%xmm1		, 0x010(%%rdx)		\n\t	movaps		%%xmm9		, 0x010(%%r12)	\n\t"\
		"addpd		%%xmm0		,%%xmm5				\n\t	addpd		%%xmm8		,%%xmm12		\n\t"\
		"addpd		%%xmm1		,%%xmm4				\n\t	addpd		%%xmm9		,%%xmm13		\n\t"\
		"movaps		%%xmm5		,      (%%rdx)		\n\t	movaps		%%xmm12		,      (%%r10)	\n\t"\
		"movaps		%%xmm4		, 0x010(%%rbx)		\n\t	movaps		%%xmm13		, 0x010(%%r10)	\n\t"\
		"/*...Block 2: t02,t12,t22,t32	*/\n\t"\
		"subq		$0x020		,%%rax				\n\t	subpd		%%xmm15		,%%xmm10		\n\t"\
		"subq		$0x020		,%%rbx				\n\t	subpd		%%xmm14		,%%xmm11		\n\t"\
		"subq		$0x020		,%%rcx				\n\t	mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"subq		$0x020		,%%rdx				\n\t	mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"addq		$0x30		,%%rsi	/* cc1 */	\n\t	movaps		%%xmm10		,      (%%r11)	\n\t"\
		"movaps		      (%%rcx)	,%%xmm4			\n\t	movaps		%%xmm11		, 0x010(%%r13)	\n\t"\
		"movaps		      (%%rdx)	,%%xmm6			\n\t	addpd		%%xmm10		,%%xmm15		\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5			\n\t	addpd		%%xmm11		,%%xmm14		\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm7			\n\t	movaps		%%xmm15		,      (%%r13)	\n\t"\
		"movaps		      (%%rcx)	,%%xmm0			\n\t	movaps		%%xmm14		, 0x010(%%r11)	\n\t"\
		"													/*...Block 6: t0A,t1A,t2A,t3A	*/\n\t"\
		"movaps		      (%%rdx)	,%%xmm2			\n\t	subq		$0x020		,%%r10			\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm1			\n\t	subq		$0x020		,%%r11			\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm3			\n\t	subq		$0x020		,%%r12			\n\t"\
		"mulpd		      (%%rsi)	,%%xmm4			\n\t	subq		$0x020		,%%r13			\n\t"\
		"mulpd		 0x20 (%%rsi)	,%%xmm6			\n\t	movaps		      (%%r12)	,%%xmm12	\n\t"\
		"mulpd		 0x010(%%rsi)	,%%xmm1			\n\t	movaps		      (%%r13)	,%%xmm14	\n\t"\
		"mulpd		 0x30 (%%rsi)	,%%xmm3			\n\t	movaps		 0x010(%%r12)	,%%xmm13	\n\t"\
		"mulpd		      (%%rsi)	,%%xmm5			\n\t	movaps		 0x010(%%r13)	,%%xmm15	\n\t"\
		"mulpd		 0x20 (%%rsi)	,%%xmm7			\n\t	movaps		      (%%r12)	,%%xmm8		\n\t"\
		"mulpd		 0x010(%%rsi)	,%%xmm0			\n\t	movaps		      (%%r13)	,%%xmm10	\n\t"\
		"mulpd		 0x30 (%%rsi)	,%%xmm2			\n\t	movaps		 0x010(%%r12)	,%%xmm9		\n\t"\
		"subpd		%%xmm1		,%%xmm4				\n\t	movaps		 0x010(%%r13)	,%%xmm11	\n\t"\
		"subpd		%%xmm3		,%%xmm6				\n\t	mulpd		 0x30 (%%rsi)	,%%xmm12	\n\t"\
		"addpd		%%xmm0		,%%xmm5				\n\t	mulpd		      (%%rsi)	,%%xmm14	\n\t"\
		"addpd		%%xmm2		,%%xmm7				\n\t	mulpd		 0x20 (%%rsi)	,%%xmm9		\n\t"\
		"subpd		%%xmm6		,%%xmm4				\n\t	mulpd		 0x010(%%rsi)	,%%xmm11	\n\t"\
		"subpd		%%xmm7		,%%xmm5				\n\t	mulpd		 0x30 (%%rsi)	,%%xmm13	\n\t"\
		"mulpd		(%%rdi)		,%%xmm6				\n\t	mulpd		      (%%rsi)	,%%xmm15	\n\t"\
		"mulpd		(%%rdi)		,%%xmm7				\n\t	mulpd		 0x20 (%%rsi)	,%%xmm8		\n\t"\
		"addpd		%%xmm4		,%%xmm6				\n\t	mulpd		 0x010(%%rsi)	,%%xmm10	\n\t"\
		"addpd		%%xmm5		,%%xmm7				\n\t	subpd		%%xmm9		,%%xmm12		\n\t"\
		"movaps		      (%%rbx)	,%%xmm2			\n\t	addpd		%%xmm11		,%%xmm14		\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm0			\n\t	addpd		%%xmm8		,%%xmm13		\n\t"\
		"movaps		      (%%rbx)	,%%xmm1			\n\t	subpd		%%xmm10		,%%xmm15		\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm3			\n\t	subpd		%%xmm14		,%%xmm12		\n\t"\
		"mulpd		-0x020(%%rsi)	,%%xmm2			\n\t	subpd		%%xmm15		,%%xmm13		\n\t"\
		"mulpd		-0x010(%%rsi)	,%%xmm0			\n\t	mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"mulpd		-0x020(%%rsi)	,%%xmm3			\n\t	mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		-0x010(%%rsi)	,%%xmm1			\n\t	addpd		%%xmm12		,%%xmm14		\n\t"\
		"subpd		%%xmm0		,%%xmm2				\n\t	addpd		%%xmm13		,%%xmm15		\n\t"\
		"addpd		%%xmm1		,%%xmm3				\n\t	movaps		      (%%r11)	,%%xmm10	\n\t"\
		"movaps		      (%%rax)	,%%xmm0			\n\t	movaps		 0x010(%%r11)	,%%xmm8		\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1			\n\t	movaps		      (%%r11)	,%%xmm9		\n\t"\
		"subpd		%%xmm2		,%%xmm0				\n\t	movaps		 0x010(%%r11)	,%%xmm11	\n\t"\
		"subpd		%%xmm3		,%%xmm1				\n\t	mulpd		-0x010(%%rsi)	,%%xmm10	\n\t"\
		"addpd		      (%%rax)	,%%xmm2			\n\t	mulpd		-0x020(%%rsi)	,%%xmm8		\n\t"\
		"addpd		 0x010(%%rax)	,%%xmm3			\n\t	mulpd		-0x010(%%rsi)	,%%xmm11	\n\t"\
		"subpd		%%xmm6		,%%xmm2				\n\t	mulpd		-0x020(%%rsi)	,%%xmm9		\n\t"\
		"subpd		%%xmm7		,%%xmm3				\n\t	addpd		%%xmm8		,%%xmm10		\n\t"\
		"mulpd		(%%rdi)		,%%xmm6				\n\t	subpd		%%xmm9		,%%xmm11		\n\t"\
		"mulpd		(%%rdi)		,%%xmm7				\n\t	movaps		      (%%r10)	,%%xmm8		\n\t"\
		"movaps		%%xmm2		,      (%%rcx)		\n\t	movaps		 0x010(%%r10)	,%%xmm9		\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)		\n\t	subpd		%%xmm10		,%%xmm8			\n\t"\
		"addpd		%%xmm2		,%%xmm6				\n\t	subpd		%%xmm11		,%%xmm9			\n\t"\
		"addpd		%%xmm3		,%%xmm7				\n\t	addpd		      (%%r10)	,%%xmm10	\n\t"\
		"movaps		%%xmm6		,      (%%rax)		\n\t	addpd		 0x010(%%r10)	,%%xmm11	\n\t"\
		"movaps		%%xmm7		, 0x010(%%rax)		\n\t	subpd		%%xmm12		,%%xmm8			\n\t"\
		"subpd		%%xmm5		,%%xmm0				\n\t	subpd		%%xmm13		,%%xmm9			\n\t"\
		"subpd		%%xmm4		,%%xmm1				\n\t	mulpd		(%%rdi)		,%%xmm12		\n\t"\
		"mulpd		(%%rdi)		,%%xmm5				\n\t	mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"mulpd		(%%rdi)		,%%xmm4				\n\t	movaps		%%xmm8		,      (%%r12)	\n\t"\
		"movaps		%%xmm0		,      (%%rbx)		\n\t	movaps		%%xmm9		, 0x010(%%r12)	\n\t"\
		"movaps		%%xmm1		, 0x010(%%rdx)		\n\t	addpd		%%xmm8		,%%xmm12		\n\t"\
		"addpd		%%xmm0		,%%xmm5				\n\t	addpd		%%xmm9		,%%xmm13		\n\t"\
		"addpd		%%xmm1		,%%xmm4				\n\t	movaps		%%xmm12		,      (%%r10)	\n\t"\
		"movaps		%%xmm5		,      (%%rdx)		\n\t	movaps		%%xmm13		, 0x010(%%r10)	\n\t"\
		"movaps		%%xmm4		, 0x010(%%rbx)		\n\t	subpd		%%xmm15		,%%xmm10		\n\t"\
		"/*...Block 4: t06,t16,t26,t36	*/\n\t"\
		"addq		$0x040		,%%rax				\n\t	subpd		%%xmm14		,%%xmm11		\n\t"\
		"addq		$0x040		,%%rbx				\n\t	addpd		%%xmm15		,%%xmm15		\n\t"\
		"addq		$0x040		,%%rcx				\n\t	addpd		%%xmm14		,%%xmm14		\n\t"\
		"addq		$0x040		,%%rdx				\n\t	movaps		%%xmm10		,      (%%r11)	\n\t"\
		"movaps		      (%%rcx)	,%%xmm4			\n\t	movaps		%%xmm11		, 0x010(%%r13)	\n\t"\
		"movaps		      (%%rdx)	,%%xmm6			\n\t	addpd		%%xmm10		,%%xmm15		\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5			\n\t	addpd		%%xmm11		,%%xmm14		\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm7			\n\t	movaps		%%xmm15		,      (%%r13)	\n\t"\
		"movaps		      (%%rcx)	,%%xmm0			\n\t	movaps		%%xmm14		, 0x010(%%r11)	\n\t"\
		"													/*...Block 8: t0E,t1E,t2E,t3E	*/		\n\t"\
		"movaps		      (%%rdx)	,%%xmm2			\n\t	addq		$0x040		,%%r10			\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm1			\n\t	addq		$0x040		,%%r11			\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm3			\n\t	addq		$0x040		,%%r12			\n\t"\
		"mulpd		 0x20 (%%rsi)	,%%xmm4			\n\t	addq		$0x040		,%%r13			\n\t"\
		"mulpd		 0x010(%%rsi)	,%%xmm6			\n\t	movaps		      (%%r12)	,%%xmm12	\n\t"\
		"mulpd		 0x30 (%%rsi)	,%%xmm1			\n\t	movaps		      (%%r13)	,%%xmm14	\n\t"\
		"mulpd		      (%%rsi)	,%%xmm3			\n\t	movaps		 0x010(%%r12)	,%%xmm13	\n\t"\
		"mulpd		 0x20 (%%rsi)	,%%xmm5			\n\t	movaps		 0x010(%%r13)	,%%xmm15	\n\t"\
		"mulpd		 0x010(%%rsi)	,%%xmm7			\n\t	movaps		      (%%r12)	,%%xmm8		\n\t"\
		"mulpd		 0x30 (%%rsi)	,%%xmm0			\n\t	movaps		      (%%r13)	,%%xmm10	\n\t"\
		"mulpd		      (%%rsi)	,%%xmm2			\n\t	movaps		 0x010(%%r12)	,%%xmm9		\n\t"\
		"subpd		%%xmm1		,%%xmm4				\n\t	movaps		 0x010(%%r13)	,%%xmm11	\n\t"\
		"addpd		%%xmm3		,%%xmm6				\n\t	mulpd		 0x010(%%rsi)	,%%xmm12	\n\t"\
		"addpd		%%xmm0		,%%xmm5				\n\t	mulpd		 0x30 (%%rsi)	,%%xmm14	\n\t"\
		"subpd		%%xmm2		,%%xmm7				\n\t	mulpd		      (%%rsi)	,%%xmm9		\n\t"\
		"subpd		%%xmm6		,%%xmm4				\n\t	mulpd		 0x20 (%%rsi)	,%%xmm11	\n\t"\
		"subpd		%%xmm7		,%%xmm5				\n\t	mulpd		 0x010(%%rsi)	,%%xmm13	\n\t"\
		"mulpd		(%%rdi)		,%%xmm6				\n\t	mulpd		 0x30 (%%rsi)	,%%xmm15	\n\t"\
		"mulpd		(%%rdi)		,%%xmm7				\n\t	mulpd		      (%%rsi)	,%%xmm8		\n\t"\
		"addpd		%%xmm4		,%%xmm6				\n\t	mulpd		 0x20 (%%rsi)	,%%xmm10	\n\t"\
		"addpd		%%xmm5		,%%xmm7				\n\t	subpd		%%xmm9		,%%xmm12		\n\t"\
		"movaps		      (%%rbx)	,%%xmm2			\n\t	subpd		%%xmm11		,%%xmm14		\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm0			\n\t	addpd		%%xmm8		,%%xmm13		\n\t"\
		"movaps		      (%%rbx)	,%%xmm1			\n\t	addpd		%%xmm10		,%%xmm15		\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm3			\n\t	subpd		%%xmm14		,%%xmm12		\n\t"\
		"mulpd		-0x010(%%rsi)	,%%xmm2			\n\t	subpd		%%xmm15		,%%xmm13		\n\t"\
		"mulpd		-0x020(%%rsi)	,%%xmm0			\n\t	mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"mulpd		-0x010(%%rsi)	,%%xmm3			\n\t	mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"mulpd		-0x020(%%rsi)	,%%xmm1			\n\t	addpd		%%xmm12		,%%xmm14		\n\t"\
		"subpd		%%xmm0		,%%xmm2				\n\t	addpd		%%xmm13		,%%xmm15		\n\t"\
		"addpd		%%xmm1		,%%xmm3				\n\t	movaps		      (%%r11)	,%%xmm10	\n\t"\
		"movaps		      (%%rax)	,%%xmm0			\n\t	movaps		 0x010(%%r11)	,%%xmm8		\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1			\n\t	movaps		      (%%r11)	,%%xmm9		\n\t"\
		"subpd		%%xmm2		,%%xmm0				\n\t	movaps		 0x010(%%r11)	,%%xmm11	\n\t"\
		"subpd		%%xmm3		,%%xmm1				\n\t	mulpd		-0x020(%%rsi)	,%%xmm10	\n\t"\
		"addpd		      (%%rax)	,%%xmm2			\n\t	mulpd		-0x010(%%rsi)	,%%xmm8		\n\t"\
		"addpd		 0x010(%%rax)	,%%xmm3			\n\t	mulpd		-0x020(%%rsi)	,%%xmm11	\n\t"\
		"subpd		%%xmm4		,%%xmm2				\n\t	mulpd		-0x010(%%rsi)	,%%xmm9		\n\t"\
		"subpd		%%xmm5		,%%xmm3				\n\t	addpd		%%xmm8		,%%xmm10		\n\t"\
		"mulpd		(%%rdi)		,%%xmm4				\n\t	subpd		%%xmm9		,%%xmm11		\n\t"\
		"mulpd		(%%rdi)		,%%xmm5				\n\t	movaps		      (%%r10)	,%%xmm8		\n\t"\
		"movaps		%%xmm2		,      (%%rcx)		\n\t	movaps		 0x010(%%r10)	,%%xmm9		\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)		\n\t	subpd		%%xmm10		,%%xmm8			\n\t"\
		"addpd		%%xmm2		,%%xmm4				\n\t	subpd		%%xmm11		,%%xmm9			\n\t"\
		"addpd		%%xmm3		,%%xmm5				\n\t	addpd		      (%%r10)	,%%xmm10	\n\t"\
		"movaps		%%xmm4		,      (%%rax)		\n\t	addpd		 0x010(%%r10)	,%%xmm11	\n\t"\
		"movaps		%%xmm5		, 0x010(%%rax)		\n\t	subpd		%%xmm12		,%%xmm8			\n\t"\
		"subpd		%%xmm7		,%%xmm0				\n\t	subpd		%%xmm13		,%%xmm9			\n\t"\
		"subpd		%%xmm6		,%%xmm1				\n\t	mulpd		(%%rdi)		,%%xmm12		\n\t"\
		"mulpd		(%%rdi)		,%%xmm7				\n\t	mulpd		(%%rdi)		,%%xmm13		\n\t"\
		"mulpd		(%%rdi)		,%%xmm6				\n\t	movaps		%%xmm8		,      (%%r12)	\n\t"\
		"movaps		%%xmm0		,      (%%rbx)		\n\t	movaps		%%xmm9		, 0x010(%%r12)	\n\t"\
		"movaps		%%xmm1		, 0x010(%%rdx)		\n\t	addpd		%%xmm8		,%%xmm12		\n\t"\
		"addpd		%%xmm0		,%%xmm7				\n\t	addpd		%%xmm9		,%%xmm13		\n\t"\
		"addpd		%%xmm1		,%%xmm6				\n\t	movaps		%%xmm12		,      (%%r10)	\n\t"\
		"movaps		%%xmm7		,      (%%rdx)		\n\t	movaps		%%xmm13		, 0x010(%%r10)	\n\t"\
		"movaps		%%xmm6		, 0x010(%%rbx)		\n\t	subpd		%%xmm15		,%%xmm10		\n\t"\
		"													subpd		%%xmm14		,%%xmm11		\n\t"\
		"													mulpd		(%%rdi)		,%%xmm15		\n\t"\
		"													mulpd		(%%rdi)		,%%xmm14		\n\t"\
		"													movaps		%%xmm10		,      (%%r11)	\n\t"\
		"													movaps		%%xmm11		, 0x010(%%r13)	\n\t"\
		"													addpd		%%xmm10		,%%xmm15		\n\t"\
		"													addpd		%%xmm11		,%%xmm14		\n\t"\
		"													movaps		%%xmm15		,      (%%r13)	\n\t"\
		"													movaps		%%xmm14		, 0x010(%%r11)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r30] "m" (Xr30)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__cc0] "m" (Xcc0)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c07] "m" (Xc07)\
		: "rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","r13","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"		/* Clobbered registers */\
	);\
	}

	#define SSE2_RADIX32_WRAPPER_DIT(Xadd0,Xadd1,Xadd2,Xisrt2,Xr00,Xr08,Xr10,Xr20,Xr28,Xr30,Xc00,Xc01,Xc02,Xc03,Xc04,Xc05,Xc06,Xc07,Xc08,Xc0A,Xc0C,Xc0E,Xc10,Xc12,Xc14,Xc16,Xc18,Xc1A,Xc1C,Xc1E)\
	{\
	__asm__ volatile (\
		"/************************************************************************/\n\t"\
		"/*...And do an inverse DIT radix-32 pass on the squared-data blocks.	*/\n\t"\
		"/************************************************************************/\n\t"\
		"/*...Block 1: */\n\t"\
		"movq		%[__isrt2]		,%%rsi\n\t"\
		"movq		%[__add2]		,%%rax\n\t"\
		"movq		%%rax		,%%rbx\n\t"\
		"movq		%%rax		,%%rcx\n\t"\
		"movq		%%rax		,%%rdx\n\t"\
		"addq		$0x200		,%%rbx\n\t"\
		"addq		$0x100		,%%rcx\n\t"\
		"addq		$0x300		,%%rdx\n\t"\
		"/*****	SSE2_RADIX4_DIT_IN_PLACE()	*****/\n\t"\
		"movaps		      (%%rax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1\n\t"\
		"movaps		      (%%rax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3\n\t"\
		"addpd		      (%%rbx)	,%%xmm0\n\t"\
		"addpd		 0x010(%%rbx)	,%%xmm1\n\t"\
		"subpd		      (%%rbx)	,%%xmm2\n\t"\
		"subpd		 0x010(%%rbx)	,%%xmm3\n\t"\
		"movaps		      (%%rcx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5\n\t"\
		"movaps		      (%%rcx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm7\n\t"\
		"addpd		      (%%rdx)	,%%xmm4\n\t"\
		"addpd		 0x010(%%rdx)	,%%xmm5\n\t"\
		"subpd		      (%%rdx)	,%%xmm6\n\t"\
		"subpd		 0x010(%%rdx)	,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"movaps		%%xmm0		,      (%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rbx)\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%rax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%rax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"movaps		%%xmm2		,      (%%rdx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		,      (%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)\n\t"\
		"addq		$0x080		,%%rax\n\t"\
		"addq		$0x080		,%%rbx\n\t"\
		"addq		$0x080		,%%rcx\n\t"\
		"addq		$0x080		,%%rdx\n\t"\
		"movaps		      (%%rax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1\n\t"\
		"movaps		      (%%rax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3\n\t"\
		"addpd		      (%%rbx)	,%%xmm0\n\t"\
		"addpd		 0x010(%%rbx)	,%%xmm1\n\t"\
		"subpd		      (%%rbx)	,%%xmm2\n\t"\
		"subpd		 0x010(%%rbx)	,%%xmm3\n\t"\
		"movaps		      (%%rcx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5\n\t"\
		"movaps		      (%%rcx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm7\n\t"\
		"addpd		      (%%rdx)	,%%xmm4\n\t"\
		"addpd		 0x010(%%rdx)	,%%xmm5\n\t"\
		"subpd		      (%%rdx)	,%%xmm6\n\t"\
		"subpd		 0x010(%%rdx)	,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"movaps		%%xmm0		,      (%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rbx)\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%rax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%rax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm6		,%%xmm1\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm0\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"mulpd		      (%%rsi)	,%%xmm3\n\t"\
		"mulpd		      (%%rsi)	,%%xmm6\n\t"\
		"mulpd		      (%%rsi)	,%%xmm0\n\t"\
		"mulpd		      (%%rsi)	,%%xmm1\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)\n\t"\
		"movaps		%%xmm0		,      (%%rcx)\n\t"\
		"movaps		%%xmm1		,      (%%rdx)\n\t"\
		"/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r10,r20,r30,r08,r18,r28,r38)	*****/\n\t"\
		"movaps		-0x080(%%rax)	,%%xmm0\n\t"\
		"movaps		-0x080(%%rbx)	,%%xmm4\n\t"\
		"movaps		-0x070(%%rax)	,%%xmm1\n\t"\
		"movaps		-0x070(%%rbx)	,%%xmm5\n\t"\
		"movaps		      (%%rax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm7\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3\n\t"\
		"movaps		      (%%rbx)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%rax)\n\t"\
		"movaps		%%xmm4		,      (%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rax)\n\t"\
		"movaps		%%xmm5		,-0x070(%%rbx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%rax)\n\t"\
		"movaps		%%xmm7		,-0x080(%%rbx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%rax)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rbx)\n\t"\
		"movaps		-0x080(%%rcx)	,%%xmm0\n\t"\
		"movaps		-0x080(%%rdx)	,%%xmm4\n\t"\
		"movaps		-0x070(%%rcx)	,%%xmm1\n\t"\
		"movaps		-0x070(%%rdx)	,%%xmm5\n\t"\
		"movaps		      (%%rcx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm7\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm3\n\t"\
		"movaps		      (%%rdx)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%rcx)\n\t"\
		"movaps		%%xmm4		,      (%%rdx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rcx)\n\t"\
		"movaps		%%xmm5		,-0x070(%%rdx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%rcx)\n\t"\
		"movaps		%%xmm7		,-0x080(%%rdx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)\n\t"\
		"/*...Block 2:	*/\n\t"\
		"subq		$0x040		,%%rax\n\t"\
		"subq		$0x040		,%%rbx\n\t"\
		"subq		$0x040		,%%rcx\n\t"\
		"subq		$0x040		,%%rdx\n\t"\
		"movaps		      (%%rax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1\n\t"\
		"movaps		      (%%rax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3\n\t"\
		"addpd		      (%%rbx)	,%%xmm0\n\t"\
		"addpd		 0x010(%%rbx)	,%%xmm1\n\t"\
		"subpd		      (%%rbx)	,%%xmm2\n\t"\
		"subpd		 0x010(%%rbx)	,%%xmm3\n\t"\
		"movaps		      (%%rcx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5\n\t"\
		"movaps		      (%%rcx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm7\n\t"\
		"addpd		      (%%rdx)	,%%xmm4\n\t"\
		"addpd		 0x010(%%rdx)	,%%xmm5\n\t"\
		"subpd		      (%%rdx)	,%%xmm6\n\t"\
		"subpd		 0x010(%%rdx)	,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"movaps		%%xmm0		,      (%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rbx)\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%rax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%rax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"movaps		%%xmm2		,      (%%rdx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		,      (%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)\n\t"\
		"addq		$0x080		,%%rax\n\t"\
		"addq		$0x080		,%%rbx\n\t"\
		"addq		$0x080		,%%rcx\n\t"\
		"addq		$0x080		,%%rdx\n\t"\
		"movaps		      (%%rax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1\n\t"\
		"movaps		      (%%rax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3\n\t"\
		"addpd		      (%%rbx)	,%%xmm0\n\t"\
		"addpd		 0x010(%%rbx)	,%%xmm1\n\t"\
		"subpd		      (%%rbx)	,%%xmm2\n\t"\
		"subpd		 0x010(%%rbx)	,%%xmm3\n\t"\
		"movaps		      (%%rcx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5\n\t"\
		"movaps		      (%%rcx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm7\n\t"\
		"addpd		      (%%rdx)	,%%xmm4\n\t"\
		"addpd		 0x010(%%rdx)	,%%xmm5\n\t"\
		"subpd		      (%%rdx)	,%%xmm6\n\t"\
		"subpd		 0x010(%%rdx)	,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"movaps		%%xmm0		,      (%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rbx)\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%rax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%rax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm6		,%%xmm1\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm0\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"mulpd		      (%%rsi)	,%%xmm3\n\t"\
		"mulpd		      (%%rsi)	,%%xmm6\n\t"\
		"mulpd		      (%%rsi)	,%%xmm0\n\t"\
		"mulpd		      (%%rsi)	,%%xmm1\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)\n\t"\
		"movaps		%%xmm0		,      (%%rcx)\n\t"\
		"movaps		%%xmm1		,      (%%rdx)\n\t"\
		"/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r04,r14,r24,r34,r0C,r1C,r2C,r3C)	*****/\n\t"\
		"movaps		-0x080(%%rax)	,%%xmm0\n\t"\
		"movaps		-0x080(%%rbx)	,%%xmm4\n\t"\
		"movaps		-0x070(%%rax)	,%%xmm1\n\t"\
		"movaps		-0x070(%%rbx)	,%%xmm5\n\t"\
		"movaps		      (%%rax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm7\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3\n\t"\
		"movaps		      (%%rbx)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%rax)\n\t"\
		"movaps		%%xmm4		,      (%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rax)\n\t"\
		"movaps		%%xmm5		,-0x070(%%rbx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%rax)\n\t"\
		"movaps		%%xmm7		,-0x080(%%rbx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%rax)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rbx)\n\t"\
		"movaps		-0x080(%%rcx)	,%%xmm0\n\t"\
		"movaps		-0x080(%%rdx)	,%%xmm4\n\t"\
		"movaps		-0x070(%%rcx)	,%%xmm1\n\t"\
		"movaps		-0x070(%%rdx)	,%%xmm5\n\t"\
		"movaps		      (%%rcx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm7\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm3\n\t"\
		"movaps		      (%%rdx)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%rcx)\n\t"\
		"movaps		%%xmm4		,      (%%rdx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rcx)\n\t"\
		"movaps		%%xmm5		,-0x070(%%rdx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%rcx)\n\t"\
		"movaps		%%xmm7		,-0x080(%%rdx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)\n\t"\
		"/*...Block 3:	*/\n\t"\
		"subq		$0x0a0		,%%rax\n\t"\
		"subq		$0x0a0		,%%rbx\n\t"\
		"subq		$0x0a0		,%%rcx\n\t"\
		"subq		$0x0a0		,%%rdx\n\t"\
		"movaps		      (%%rax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1\n\t"\
		"movaps		      (%%rax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3\n\t"\
		"addpd		      (%%rbx)	,%%xmm0\n\t"\
		"addpd		 0x010(%%rbx)	,%%xmm1\n\t"\
		"subpd		      (%%rbx)	,%%xmm2\n\t"\
		"subpd		 0x010(%%rbx)	,%%xmm3\n\t"\
		"movaps		      (%%rcx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5\n\t"\
		"movaps		      (%%rcx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm7\n\t"\
		"addpd		      (%%rdx)	,%%xmm4\n\t"\
		"addpd		 0x010(%%rdx)	,%%xmm5\n\t"\
		"subpd		      (%%rdx)	,%%xmm6\n\t"\
		"subpd		 0x010(%%rdx)	,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"movaps		%%xmm0		,      (%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rbx)\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%rax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%rax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"movaps		%%xmm2		,      (%%rdx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		,      (%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)\n\t"\
		"addq		$0x080		,%%rax\n\t"\
		"addq		$0x080		,%%rbx\n\t"\
		"addq		$0x080		,%%rcx\n\t"\
		"addq		$0x080		,%%rdx\n\t"\
		"movaps		      (%%rax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1\n\t"\
		"movaps		      (%%rax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3\n\t"\
		"addpd		      (%%rbx)	,%%xmm0\n\t"\
		"addpd		 0x010(%%rbx)	,%%xmm1\n\t"\
		"subpd		      (%%rbx)	,%%xmm2\n\t"\
		"subpd		 0x010(%%rbx)	,%%xmm3\n\t"\
		"movaps		      (%%rcx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5\n\t"\
		"movaps		      (%%rcx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm7\n\t"\
		"addpd		      (%%rdx)	,%%xmm4\n\t"\
		"addpd		 0x010(%%rdx)	,%%xmm5\n\t"\
		"subpd		      (%%rdx)	,%%xmm6\n\t"\
		"subpd		 0x010(%%rdx)	,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"movaps		%%xmm0		,      (%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rbx)\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%rax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%rax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm6		,%%xmm1\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm0\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"mulpd		      (%%rsi)	,%%xmm3\n\t"\
		"mulpd		      (%%rsi)	,%%xmm6\n\t"\
		"mulpd		      (%%rsi)	,%%xmm0\n\t"\
		"mulpd		      (%%rsi)	,%%xmm1\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)\n\t"\
		"movaps		%%xmm0		,      (%%rcx)\n\t"\
		"movaps		%%xmm1		,      (%%rdx)\n\t"\
		"/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r02,r12,r22,r32,r0A,r1A,r2A,r3A)	*****/\n\t"\
		"movaps		-0x080(%%rax)	,%%xmm0\n\t"\
		"movaps		-0x080(%%rbx)	,%%xmm4\n\t"\
		"movaps		-0x070(%%rax)	,%%xmm1\n\t"\
		"movaps		-0x070(%%rbx)	,%%xmm5\n\t"\
		"movaps		      (%%rax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm7\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3\n\t"\
		"movaps		      (%%rbx)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%rax)\n\t"\
		"movaps		%%xmm4		,      (%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rax)\n\t"\
		"movaps		%%xmm5		,-0x070(%%rbx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%rax)\n\t"\
		"movaps		%%xmm7		,-0x080(%%rbx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%rax)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rbx)\n\t"\
		"movaps		-0x080(%%rcx)	,%%xmm0\n\t"\
		"movaps		-0x080(%%rdx)	,%%xmm4\n\t"\
		"movaps		-0x070(%%rcx)	,%%xmm1\n\t"\
		"movaps		-0x070(%%rdx)	,%%xmm5\n\t"\
		"movaps		      (%%rcx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm7\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm3\n\t"\
		"movaps		      (%%rdx)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%rcx)\n\t"\
		"movaps		%%xmm4		,      (%%rdx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rcx)\n\t"\
		"movaps		%%xmm5		,-0x070(%%rdx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%rcx)\n\t"\
		"movaps		%%xmm7		,-0x080(%%rdx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)\n\t"\
		"/*...Block 4:	*/\n\t"\
		"subq		$0x040		,%%rax\n\t"\
		"subq		$0x040		,%%rbx\n\t"\
		"subq		$0x040		,%%rcx\n\t"\
		"subq		$0x040		,%%rdx\n\t"\
		"movaps		      (%%rax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1\n\t"\
		"movaps		      (%%rax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3\n\t"\
		"addpd		      (%%rbx)	,%%xmm0\n\t"\
		"addpd		 0x010(%%rbx)	,%%xmm1\n\t"\
		"subpd		      (%%rbx)	,%%xmm2\n\t"\
		"subpd		 0x010(%%rbx)	,%%xmm3\n\t"\
		"movaps		      (%%rcx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5\n\t"\
		"movaps		      (%%rcx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm7\n\t"\
		"addpd		      (%%rdx)	,%%xmm4\n\t"\
		"addpd		 0x010(%%rdx)	,%%xmm5\n\t"\
		"subpd		      (%%rdx)	,%%xmm6\n\t"\
		"subpd		 0x010(%%rdx)	,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"movaps		%%xmm0		,      (%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rbx)\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%rax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%rax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"movaps		%%xmm2		,      (%%rdx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		,      (%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)\n\t"\
		"addq		$0x080		,%%rax\n\t"\
		"addq		$0x080		,%%rbx\n\t"\
		"addq		$0x080		,%%rcx\n\t"\
		"addq		$0x080		,%%rdx\n\t"\
		"movaps		      (%%rax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1\n\t"\
		"movaps		      (%%rax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3\n\t"\
		"addpd		      (%%rbx)	,%%xmm0\n\t"\
		"addpd		 0x010(%%rbx)	,%%xmm1\n\t"\
		"subpd		      (%%rbx)	,%%xmm2\n\t"\
		"subpd		 0x010(%%rbx)	,%%xmm3\n\t"\
		"movaps		      (%%rcx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm5\n\t"\
		"movaps		      (%%rcx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm7\n\t"\
		"addpd		      (%%rdx)	,%%xmm4\n\t"\
		"addpd		 0x010(%%rdx)	,%%xmm5\n\t"\
		"subpd		      (%%rdx)	,%%xmm6\n\t"\
		"subpd		 0x010(%%rdx)	,%%xmm7\n\t"\
		"subpd		%%xmm4		,%%xmm0\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"movaps		%%xmm0		,      (%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rbx)\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm4		,      (%%rax)\n\t"\
		"movaps		%%xmm5		, 0x010(%%rax)\n\t"\
		"subpd		%%xmm7		,%%xmm2\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm6		,%%xmm1\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm0\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"mulpd		      (%%rsi)	,%%xmm3\n\t"\
		"mulpd		      (%%rsi)	,%%xmm6\n\t"\
		"mulpd		      (%%rsi)	,%%xmm0\n\t"\
		"mulpd		      (%%rsi)	,%%xmm1\n\t"\
		"movaps		%%xmm3		, 0x010(%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)\n\t"\
		"movaps		%%xmm0		,      (%%rcx)\n\t"\
		"movaps		%%xmm1		,      (%%rdx)\n\t"\
		"/*****	SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r06,r16,r26,r36,r0E,r1E,r2E,r3E)	*****/\n\t"\
		"movaps		-0x080(%%rax)	,%%xmm0\n\t"\
		"movaps		-0x080(%%rbx)	,%%xmm4\n\t"\
		"movaps		-0x070(%%rax)	,%%xmm1\n\t"\
		"movaps		-0x070(%%rbx)	,%%xmm5\n\t"\
		"movaps		      (%%rax)	,%%xmm2\n\t"\
		"movaps		 0x010(%%rbx)	,%%xmm7\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm3\n\t"\
		"movaps		      (%%rbx)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%rax)\n\t"\
		"movaps		%%xmm4		,      (%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rax)\n\t"\
		"movaps		%%xmm5		,-0x070(%%rbx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%rax)\n\t"\
		"movaps		%%xmm7		,-0x080(%%rbx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%rax)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rbx)\n\t"\
		"movaps		-0x080(%%rcx)	,%%xmm0\n\t"\
		"movaps		-0x080(%%rdx)	,%%xmm4\n\t"\
		"movaps		-0x070(%%rcx)	,%%xmm1\n\t"\
		"movaps		-0x070(%%rdx)	,%%xmm5\n\t"\
		"movaps		      (%%rcx)	,%%xmm2\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm7\n\t"\
		"movaps		 0x010(%%rcx)	,%%xmm3\n\t"\
		"movaps		      (%%rdx)	,%%xmm6\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm4		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm6\n\t"\
		"movaps		%%xmm0		,      (%%rcx)\n\t"\
		"movaps		%%xmm4		,      (%%rdx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rcx)\n\t"\
		"movaps		%%xmm5		,-0x070(%%rdx)\n\t"\
		"movaps		%%xmm2		,-0x080(%%rcx)\n\t"\
		"movaps		%%xmm7		,-0x080(%%rdx)\n\t"\
		"movaps		%%xmm3		,-0x070(%%rcx)\n\t"\
		"movaps		%%xmm6		, 0x010(%%rdx)\n\t"\
		"/***************************************************************************************/\n\t"\
		"/* Now do eight more radix-4 transforms, including the internal and external twiddles: */\n\t"\
		"/***************************************************************************************/\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 2: t02,t12,t22,t32 -> r10,14,12,16:	*/\n\t"\
		"/************************************************/\n\t"\
		"movq		%[__isrt2]		,%%rsi\n\t"\
		"movq		%[__r10]		,%%rax\n\t"\
		"movq		%[__add1]		,%%rbx\n\t"\
		"movq		%%rsi			,%%rcx\n\t"\
		"movq		%%rsi			,%%rdx\n\t"\
		"movq		%[__c01]		,%%rdi\n\t"\
		"addq		$0x010			,%%rcx	/* cc0 */\n\t"\
		"addq		$0x030			,%%rdx	/* cc1 */\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm4\n\t"\
		"movaps		 0x060(%%rax)	,%%xmm0\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm5\n\t"\
		"movaps		 0x070(%%rax)	,%%xmm1\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm6\n\t"\
		"movaps		 0x060(%%rax)	,%%xmm2\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm7\n\t"\
		"movaps		 0x070(%%rax)	,%%xmm3\n\t"\
		"mulpd		      (%%rdx)	,%%xmm4\n\t"\
		"mulpd		 0x20 (%%rdx)	,%%xmm0\n\t"\
		"mulpd		      (%%rdx)	,%%xmm5\n\t"\
		"mulpd		 0x20 (%%rdx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%rdx)	,%%xmm6\n\t"\
		"mulpd		 0x30 (%%rdx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%rdx)	,%%xmm7\n\t"\
		"mulpd		 0x30 (%%rdx)	,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"subpd		%%xmm0		,%%xmm4\n\t"\
		"subpd		%%xmm1		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		 0x040(%%rax)	,%%xmm2\n\t"\
		"movaps		 0x050(%%rax)	,%%xmm3\n\t"\
		"movaps		 0x040(%%rax)	,%%xmm0\n\t"\
		"movaps		 0x050(%%rax)	,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm0\n\t"\
		"addpd		%%xmm1		,%%xmm2\n\t"\
		"subpd		%%xmm0		,%%xmm3\n\t"\
		"movaps		      (%%rax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm2\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm2		,%%xmm6\n\t"\
		"addpd		%%xmm3		,%%xmm7\n\t"\
		"movaps		%%xmm2		, 0x020(%%rax)\n\t"\
		"movaps		%%xmm3		, 0x030(%%rax)\n\t"\
		"movaps		%%xmm6		,%%xmm2\n\t"\
		"movaps		%%xmm7		,%%xmm3\n\t"\
		"mulpd		      (%%rdi)	,%%xmm6\n\t"\
		"mulpd		      (%%rdi)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		, 0x010(%%rbx)\n\t"\
		"movaps		%%xmm6		,      (%%rbx)\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm6\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm7\n\t"\
		"movaps		%%xmm6		,%%xmm2\n\t"\
		"movaps		%%xmm7		,%%xmm3\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm6\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm7\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm2\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		, 0x110(%%rbx)\n\t"\
		"movaps		%%xmm6		, 0x100(%%rbx)\n\t"\
		"addq		$0x040		,%%rdi\n\t"\
		"subpd		%%xmm5		,%%xmm0\n\t"\
		"subpd		%%xmm4		,%%xmm1\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm0		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"mulpd		      (%%rdi)	,%%xmm5\n\t"\
		"mulpd		      (%%rdi)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm3		,%%xmm5\n\t"\
		"movaps		%%xmm1		, 0x090(%%rbx)\n\t"\
		"movaps		%%xmm5		, 0x080(%%rbx)\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm4		,%%xmm3\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm0\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm4\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm2\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm4		, 0x190(%%rbx)\n\t"\
		"movaps		%%xmm0		, 0x180(%%rbx)\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 6: t0A,t1A,t2A,t3A -> r18,1C,1A,1E:	*/\n\t"\
		"/************************************************/\n\t"\
		"addq		$0x080		,%%rax\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm4\n\t"\
		"movaps		 0x060(%%rax)	,%%xmm0\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm5\n\t"\
		"movaps		 0x070(%%rax)	,%%xmm1\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm6\n\t"\
		"movaps		 0x060(%%rax)	,%%xmm2\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm7\n\t"\
		"movaps		 0x070(%%rax)	,%%xmm3\n\t"\
		"mulpd		 0x30 (%%rdx)	,%%xmm4\n\t"\
		"mulpd		      (%%rdx)	,%%xmm0\n\t"\
		"mulpd		 0x30 (%%rdx)	,%%xmm5\n\t"\
		"mulpd		      (%%rdx)	,%%xmm1\n\t"\
		"mulpd		 0x20 (%%rdx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%rdx)	,%%xmm2\n\t"\
		"mulpd		 0x20 (%%rdx)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%rdx)	,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"subpd		%%xmm0		,%%xmm6\n\t"\
		"subpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		 0x040(%%rax)	,%%xmm2\n\t"\
		"movaps		 0x050(%%rax)	,%%xmm3\n\t"\
		"movaps		 0x040(%%rax)	,%%xmm0\n\t"\
		"movaps		 0x050(%%rax)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm2\n\t"\
		"mulpd		      (%%rcx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm3\n\t"\
		"mulpd		      (%%rcx)	,%%xmm0\n\t"\
		"subpd		%%xmm1		,%%xmm2\n\t"\
		"addpd		%%xmm0		,%%xmm3\n\t"\
		"movaps		      (%%rax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addq		$0x040		,%%rdi\n\t"\
		"subpd		%%xmm6		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm1\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		%%xmm0		, 0x020(%%rax)\n\t"\
		"movaps		%%xmm1		, 0x030(%%rax)\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"mulpd		      (%%rdi)	,%%xmm6\n\t"\
		"mulpd		      (%%rdi)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movaps		%%xmm7		, 0x050(%%rbx)\n\t"\
		"movaps		%%xmm6		, 0x040(%%rbx)\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm6\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm7\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm6\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm7\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm0\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movaps		%%xmm7		, 0x150(%%rbx)\n\t"\
		"movaps		%%xmm6		, 0x140(%%rbx)\n\t"\
		"addq		$0x040		,%%rdi\n\t"\
		"subpd		%%xmm5		,%%xmm2\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm0\n\t"\
		"movaps		%%xmm3		,%%xmm1\n\t"\
		"mulpd		      (%%rdi)	,%%xmm5\n\t"\
		"mulpd		      (%%rdi)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm3\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm3		, 0x0d0(%%rbx)\n\t"\
		"movaps		%%xmm5		, 0x0c0(%%rbx)\n\t"\
		"movaps		%%xmm2		,%%xmm0\n\t"\
		"movaps		%%xmm4		,%%xmm1\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm2\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm4\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm0\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm2\n\t"\
		"movaps		%%xmm4		, 0x1d0(%%rbx)\n\t"\
		"movaps		%%xmm2		, 0x1c0(%%rbx)\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 4: t06,t16,t26,t36 -> r30,34,32,36:	*/\n\t"\
		"/************************************************/\n\t"\
		"addq		$0x180		,%%rax\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm4\n\t"\
		"movaps		 0x060(%%rax)	,%%xmm0\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm5\n\t"\
		"movaps		 0x070(%%rax)	,%%xmm1\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm6\n\t"\
		"movaps		 0x060(%%rax)	,%%xmm2\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm7\n\t"\
		"movaps		 0x070(%%rax)	,%%xmm3\n\t"\
		"mulpd		 0x20 (%%rdx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%rdx)	,%%xmm0\n\t"\
		"mulpd		 0x20 (%%rdx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%rdx)	,%%xmm1\n\t"\
		"mulpd		 0x30 (%%rdx)	,%%xmm6\n\t"\
		"mulpd		      (%%rdx)	,%%xmm2\n\t"\
		"mulpd		 0x30 (%%rdx)	,%%xmm7\n\t"\
		"mulpd		      (%%rdx)	,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm4\n\t"\
		"subpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"subpd		%%xmm0		,%%xmm6\n\t"\
		"subpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		 0x040(%%rax)	,%%xmm2\n\t"\
		"movaps		 0x050(%%rax)	,%%xmm3\n\t"\
		"movaps		 0x040(%%rax)	,%%xmm0\n\t"\
		"movaps		 0x050(%%rax)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm2\n\t"\
		"mulpd		      (%%rcx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm3\n\t"\
		"mulpd		      (%%rcx)	,%%xmm0\n\t"\
		"addpd		%%xmm1		,%%xmm2\n\t"\
		"subpd		%%xmm0		,%%xmm3\n\t"\
		"movaps		      (%%rax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addq		$0x040		,%%rdi\n\t"\
		"subpd		%%xmm6		,%%xmm2\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm2		,%%xmm6\n\t"\
		"addpd		%%xmm3		,%%xmm7\n\t"\
		"movaps		%%xmm2		, 0x020(%%rax)\n\t"\
		"movaps		%%xmm3		, 0x030(%%rax)\n\t"\
		"movaps		%%xmm6		,%%xmm2\n\t"\
		"movaps		%%xmm7		,%%xmm3\n\t"\
		"mulpd		      (%%rdi)	,%%xmm6\n\t"\
		"mulpd		      (%%rdi)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		, 0x030(%%rbx)\n\t"\
		"movaps		%%xmm6		, 0x020(%%rbx)\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm6\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm7\n\t"\
		"movaps		%%xmm6		,%%xmm2\n\t"\
		"movaps		%%xmm7		,%%xmm3\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm6\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm7\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm2\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm7\n\t"\
		"addpd		%%xmm3		,%%xmm6\n\t"\
		"movaps		%%xmm7		, 0x130(%%rbx)\n\t"\
		"movaps		%%xmm6		, 0x120(%%rbx)\n\t"\
		"addq		$0x040		,%%rdi\n\t"\
		"subpd		%%xmm5		,%%xmm0\n\t"\
		"subpd		%%xmm4		,%%xmm1\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm0		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"mulpd		      (%%rdi)	,%%xmm5\n\t"\
		"mulpd		      (%%rdi)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm3		,%%xmm5\n\t"\
		"movaps		%%xmm1		, 0x0b0(%%rbx)\n\t"\
		"movaps		%%xmm5		, 0x0a0(%%rbx)\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm4		,%%xmm3\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm0\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm4\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm2\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm4		, 0x1b0(%%rbx)\n\t"\
		"movaps		%%xmm0		, 0x1a0(%%rbx)\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 8: t0E,t1E,t2E,t3E -> r38,3C,3A,3E:	*/\n\t"\
		"/************************************************/\n\t"\
		"addq		$0x080		,%%rax\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm4\n\t"\
		"movaps		 0x060(%%rax)	,%%xmm0\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm5\n\t"\
		"movaps		 0x070(%%rax)	,%%xmm1\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm6\n\t"\
		"movaps		 0x060(%%rax)	,%%xmm2\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm7\n\t"\
		"movaps		 0x070(%%rax)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%rdx)	,%%xmm4\n\t"\
		"mulpd		 0x30 (%%rdx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%rdx)	,%%xmm5\n\t"\
		"mulpd		 0x30 (%%rdx)	,%%xmm1\n\t"\
		"mulpd		      (%%rdx)	,%%xmm6\n\t"\
		"mulpd		 0x20 (%%rdx)	,%%xmm2\n\t"\
		"mulpd		      (%%rdx)	,%%xmm7\n\t"\
		"mulpd		 0x20 (%%rdx)	,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"subpd		%%xmm0		,%%xmm6\n\t"\
		"subpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		 0x040(%%rax)	,%%xmm2\n\t"\
		"movaps		 0x050(%%rax)	,%%xmm3\n\t"\
		"movaps		 0x040(%%rax)	,%%xmm0\n\t"\
		"movaps		 0x050(%%rax)	,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm0\n\t"\
		"subpd		%%xmm1		,%%xmm2\n\t"\
		"addpd		%%xmm0		,%%xmm3\n\t"\
		"movaps		      (%%rax)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rax)	,%%xmm1\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"addq		$0x040		,%%rdi\n\t"\
		"subpd		%%xmm6		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm1\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		%%xmm0		, 0x020(%%rax)\n\t"\
		"movaps		%%xmm1		, 0x030(%%rax)\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"mulpd		      (%%rdi)	,%%xmm6\n\t"\
		"mulpd		      (%%rdi)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movaps		%%xmm7		, 0x070(%%rbx)\n\t"\
		"movaps		%%xmm6		, 0x060(%%rbx)\n\t"\
		"movaps		 0x020(%%rax)	,%%xmm6\n\t"\
		"movaps		 0x030(%%rax)	,%%xmm7\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm6\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm7\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm0\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movaps		%%xmm7		, 0x170(%%rbx)\n\t"\
		"movaps		%%xmm6		, 0x160(%%rbx)\n\t"\
		"addq		$0x040		,%%rdi\n\t"\
		"subpd		%%xmm5		,%%xmm2\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm0\n\t"\
		"movaps		%%xmm3		,%%xmm1\n\t"\
		"mulpd		      (%%rdi)	,%%xmm5\n\t"\
		"mulpd		      (%%rdi)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%rdi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm3\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm3		, 0x0f0(%%rbx)\n\t"\
		"movaps		%%xmm5		, 0x0e0(%%rbx)\n\t"\
		"movaps		%%xmm2		,%%xmm0\n\t"\
		"movaps		%%xmm4		,%%xmm1\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm2\n\t"\
		"mulpd		 0x20 (%%rdi)	,%%xmm4\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm0\n\t"\
		"mulpd		 0x30 (%%rdi)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm2\n\t"\
		"movaps		%%xmm4		, 0x1f0(%%rbx)\n\t"\
		"movaps		%%xmm2		, 0x1e0(%%rbx)\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 1: t00,t10,t20,t30 -> r00,04,02,06:	*/\n\t"\
		"/************************************************/\n\t"\
		"movq		%[__r00]		,%%rdx\n\t"\
		"movaps		      (%%rdx)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm1\n\t"\
		"movaps		 0x040(%%rdx)	,%%xmm2\n\t"\
		"movaps		 0x050(%%rdx)	,%%xmm3\n\t"\
		"subpd		 0x040(%%rdx)	,%%xmm0\n\t"\
		"subpd		 0x050(%%rdx)	,%%xmm1\n\t"\
		"addpd		      (%%rdx)	,%%xmm2\n\t"\
		"addpd		 0x010(%%rdx)	,%%xmm3\n\t"\
		"movaps		 0x020(%%rdx)	,%%xmm4\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm5\n\t"\
		"movaps		 0x060(%%rdx)	,%%xmm6\n\t"\
		"movaps		 0x070(%%rdx)	,%%xmm7\n\t"\
		"subpd		 0x060(%%rdx)	,%%xmm4\n\t"\
		"subpd		 0x070(%%rdx)	,%%xmm5\n\t"\
		"addpd		 0x020(%%rdx)	,%%xmm6\n\t"\
		"addpd		 0x030(%%rdx)	,%%xmm7\n\t"\
		"movq		%[__add0]		,%%rax\n\t"\
		"movq		%[__c10]		,%%rcx\n\t"\
		"addpd		%%xmm6		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm3\n\t"\
		"movaps		%%xmm2		,      (%%rdx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%rdx)\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm2\n\t"\
		"subpd		%%xmm7		,%%xmm3\n\t"\
		"movaps		%%xmm2		,%%xmm6\n\t"\
		"movaps		%%xmm3		,%%xmm7\n\t"\
		"mulpd		      (%%rcx)	,%%xmm2\n\t"\
		"mulpd		      (%%rcx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"addpd		%%xmm7		,%%xmm2\n\t"\
		"movq		%[__add1]		,%%rbx\n\t"\
		"movaps		%%xmm3		,%%xmm7\n\t"\
		"movaps		%%xmm2		,%%xmm6\n\t"\
		"unpckhpd	 0x110(%%rbx)	,%%xmm7\n\t"\
		"unpcklpd	 0x110(%%rbx)	,%%xmm3\n\t"\
		"movaps		%%xmm7		, 0x110(%%rbx)\n\t"\
		"unpckhpd	 0x100(%%rbx)	,%%xmm6\n\t"\
		"unpcklpd	 0x100(%%rbx)	,%%xmm2\n\t"\
		"movaps		%%xmm6		, 0x100(%%rbx)\n\t"\
		"movaps		%%xmm3		, 0x110(%%rax)\n\t"\
		"movaps		%%xmm2		, 0x100(%%rax)\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm3\n\t"\
		"movaps		      (%%rdx)	,%%xmm2\n\t"\
		"movaps		%%xmm3		,%%xmm7\n\t"\
		"movaps		%%xmm2		,%%xmm6\n\t"\
		"unpckhpd	 0x010(%%rbx)	,%%xmm7\n\t"\
		"unpcklpd	 0x010(%%rbx)	,%%xmm3\n\t"\
		"movaps		%%xmm7		, 0x010(%%rbx)\n\t"\
		"unpckhpd	      (%%rbx)	,%%xmm6\n\t"\
		"unpcklpd	      (%%rbx)	,%%xmm2\n\t"\
		"movaps		%%xmm6		,      (%%rbx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%rax)\n\t"\
		"movaps		%%xmm2		,      (%%rax)\n\t"\
		"movq		%[__c08]		,%%rcx\n\t"\
		"addpd		%%xmm5		,%%xmm0\n\t"\
		"subpd		%%xmm4		,%%xmm1\n\t"\
		"movaps		%%xmm0		,%%xmm2\n\t"\
		"movaps		%%xmm1		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"movaps		%%xmm0		,%%xmm6\n\t"\
		"movaps		%%xmm1		,%%xmm7\n\t"\
		"mulpd		      (%%rcx)	,%%xmm2\n\t"\
		"mulpd		      (%%rcx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm3\n\t"\
		"addpd		%%xmm7		,%%xmm2\n\t"\
		"movaps		%%xmm3		,%%xmm7\n\t"\
		"movaps		%%xmm2		,%%xmm6\n\t"\
		"unpckhpd	 0x090(%%rbx)	,%%xmm7\n\t"\
		"unpcklpd	 0x090(%%rbx)	,%%xmm3\n\t"\
		"movaps		%%xmm7		, 0x090(%%rbx)\n\t"\
		"unpckhpd	 0x080(%%rbx)	,%%xmm6\n\t"\
		"unpcklpd	 0x080(%%rbx)	,%%xmm2\n\t"\
		"movaps		%%xmm6		, 0x080(%%rbx)\n\t"\
		"movaps		%%xmm3		, 0x090(%%rax)\n\t"\
		"movaps		%%xmm2		, 0x080(%%rax)\n\t"\
		"movq		%[__c18]		,%%rcx\n\t"\
		"subpd		%%xmm5		,%%xmm0\n\t"\
		"addpd		%%xmm4		,%%xmm1\n\t"\
		"movaps		%%xmm0		,%%xmm6\n\t"\
		"movaps		%%xmm1		,%%xmm7\n\t"\
		"mulpd		      (%%rcx)	,%%xmm0\n\t"\
		"mulpd		      (%%rcx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm7\n\t"\
		"subpd		%%xmm6		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm0\n\t"\
		"movaps		%%xmm1		,%%xmm7\n\t"\
		"movaps		%%xmm0		,%%xmm6\n\t"\
		"unpckhpd	 0x190(%%rbx)	,%%xmm7\n\t"\
		"unpcklpd	 0x190(%%rbx)	,%%xmm1\n\t"\
		"movaps		%%xmm7		, 0x190(%%rbx)\n\t"\
		"unpckhpd	 0x180(%%rbx)	,%%xmm6\n\t"\
		"unpcklpd	 0x180(%%rbx)	,%%xmm0\n\t"\
		"movaps		%%xmm6		, 0x180(%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x190(%%rax)\n\t"\
		"movaps		%%xmm0		, 0x180(%%rax)\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 5: t08,t18,t28,t38 -> r08,0C,0A,0E	*/\n\t"\
		"/************************************************/\n\t"\
		"movq		%[__r08]		,%%rdx\n\t"\
		"movaps		 (%%rsi)	,%%xmm2\n\t	/* isrt2 */"\
		"movaps		 0x020(%%rdx)	,%%xmm4\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm5\n\t"\
		"movaps		 0x060(%%rdx)	,%%xmm0\n\t"\
		"movaps		 0x070(%%rdx)	,%%xmm1\n\t"\
		"addpd		 0x030(%%rdx)	,%%xmm4\n\t"\
		"subpd		 0x020(%%rdx)	,%%xmm5\n\t"\
		"subpd		 0x070(%%rdx)	,%%xmm0\n\t"\
		"addpd		 0x060(%%rdx)	,%%xmm1\n\t"\
		"mulpd		%%xmm2		,%%xmm4\n\t"\
		"mulpd		%%xmm2		,%%xmm5\n\t"\
		"mulpd		%%xmm2		,%%xmm0\n\t"\
		"mulpd		%%xmm2		,%%xmm1\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"subpd		%%xmm0		,%%xmm4\n\t"\
		"subpd		%%xmm1		,%%xmm5\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		      (%%rdx)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm1\n\t"\
		"movaps		 0x040(%%rdx)	,%%xmm2\n\t"\
		"movaps		 0x050(%%rdx)	,%%xmm3\n\t"\
		"subpd		 0x050(%%rdx)	,%%xmm0\n\t"\
		"subpd		 0x040(%%rdx)	,%%xmm1\n\t"\
		"addpd		      (%%rdx)	,%%xmm3\n\t"\
		"addpd		 0x010(%%rdx)	,%%xmm2\n\t"\
		"movq		%[__add0]		,%%rax\n\t"\
		"movq		%[__c04]		,%%rcx\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"subpd		%%xmm5		,%%xmm1\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm3		,      (%%rdx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rdx)\n\t"\
		"movaps		%%xmm4		,%%xmm3\n\t"\
		"movaps		%%xmm5		,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm4\n\t"\
		"mulpd		      (%%rcx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1\n\t"\
		"subpd		%%xmm3		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm4\n\t"\
		"movq		%[__add1]		,%%rbx\n\t"\
		"movaps		%%xmm5		,%%xmm3\n\t"\
		"movaps		%%xmm4		,%%xmm1\n\t"\
		"unpckhpd	 0x050(%%rbx)	,%%xmm3\n\t"\
		"unpcklpd	 0x050(%%rbx)	,%%xmm5\n\t"\
		"movaps		%%xmm3		, 0x050(%%rbx)\n\t"\
		"unpckhpd	 0x040(%%rbx)	,%%xmm1\n\t"\
		"unpcklpd	 0x040(%%rbx)	,%%xmm4\n\t"\
		"movaps		%%xmm1		, 0x040(%%rbx)\n\t"\
		"movaps		%%xmm5		, 0x050(%%rax)\n\t"\
		"movaps		%%xmm4		, 0x040(%%rax)\n\t"\
		"movq		%[__c14]		,%%rcx\n\t"\
		"movaps		      (%%rdx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm5\n\t"\
		"movaps		%%xmm4		,%%xmm3\n\t"\
		"movaps		%%xmm5		,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm4\n\t"\
		"mulpd		      (%%rcx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1\n\t"\
		"subpd		%%xmm3		,%%xmm5\n\t"\
		"addpd		%%xmm1		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm3\n\t"\
		"movaps		%%xmm4		,%%xmm1\n\t"\
		"unpckhpd	 0x150(%%rbx)	,%%xmm3\n\t"\
		"unpcklpd	 0x150(%%rbx)	,%%xmm5\n\t"\
		"movaps		%%xmm3		, 0x150(%%rbx)\n\t"\
		"unpckhpd	 0x140(%%rbx)	,%%xmm1\n\t"\
		"unpcklpd	 0x140(%%rbx)	,%%xmm4\n\t"\
		"movaps		%%xmm1		, 0x140(%%rbx)\n\t"\
		"movaps		%%xmm5		, 0x150(%%rax)\n\t"\
		"movaps		%%xmm4		, 0x140(%%rax)\n\t"\
		"movq		%[__c0C]		,%%rcx\n\t"\
		"subpd		%%xmm7		,%%xmm0\n\t"\
		"subpd		%%xmm6		,%%xmm2\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm2		,%%xmm6\n\t"\
		"movaps		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm2		,%%xmm5\n\t"\
		"mulpd		      (%%rcx)	,%%xmm7\n\t"\
		"mulpd		      (%%rcx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm5\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"addpd		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm2		,%%xmm5\n\t"\
		"movaps		%%xmm7		,%%xmm4\n\t"\
		"unpckhpd	 0x0d0(%%rbx)	,%%xmm5\n\t"\
		"unpcklpd	 0x0d0(%%rbx)	,%%xmm2\n\t"\
		"movaps		%%xmm5		, 0x0d0(%%rbx)\n\t"\
		"unpckhpd	 0x0c0(%%rbx)	,%%xmm4\n\t"\
		"unpcklpd	 0x0c0(%%rbx)	,%%xmm7\n\t"\
		"movaps		%%xmm4		, 0x0c0(%%rbx)\n\t"\
		"movaps		%%xmm2		, 0x0d0(%%rax)\n\t"\
		"movaps		%%xmm7		, 0x0c0(%%rax)\n\t"\
		"movq		%[__c1C]		,%%rcx\n\t"\
		"movaps		%%xmm0		,%%xmm4\n\t"\
		"movaps		%%xmm6		,%%xmm5\n\t"\
		"mulpd		      (%%rcx)	,%%xmm0\n\t"\
		"mulpd		      (%%rcx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm5\n\t"\
		"subpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm0\n\t"\
		"movaps		%%xmm6		,%%xmm5\n\t"\
		"movaps		%%xmm0		,%%xmm4\n\t"\
		"unpckhpd	 0x1d0(%%rbx)	,%%xmm5\n\t"\
		"unpcklpd	 0x1d0(%%rbx)	,%%xmm6\n\t"\
		"movaps		%%xmm5		, 0x1d0(%%rbx)\n\t"\
		"unpckhpd	 0x1c0(%%rbx)	,%%xmm4\n\t"\
		"unpcklpd	 0x1c0(%%rbx)	,%%xmm0\n\t"\
		"movaps		%%xmm4		, 0x1c0(%%rbx)\n\t"\
		"movaps		%%xmm6		, 0x1d0(%%rax)\n\t"\
		"movaps		%%xmm0		, 0x1c0(%%rax)\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 3: t04,t14,t24,t34 -> r20,24,22,26	*/\n\t"\
		"/************************************************/\n\t"\
		"movq		%[__r20]		,%%rdx\n\t"\
		"movq		%%rsi			,%%rcx\n\t"\
		"addq		$0x010			,%%rcx	/* cc0 */\n\t"\
		"movaps		 0x020(%%rdx)	,%%xmm4\n\t"\
		"movaps		 0x060(%%rdx)	,%%xmm0\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm5\n\t"\
		"movaps		 0x070(%%rdx)	,%%xmm1\n\t"\
		"movaps		 0x020(%%rdx)	,%%xmm6\n\t"\
		"movaps		 0x060(%%rdx)	,%%xmm2\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm7\n\t"\
		"movaps		 0x070(%%rdx)	,%%xmm3\n\t"\
		"mulpd		      (%%rcx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm0\n\t"\
		"mulpd		      (%%rcx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm6\n\t"\
		"mulpd		      (%%rcx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm7\n\t"\
		"mulpd		      (%%rcx)	,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"subpd		%%xmm0		,%%xmm6\n\t"\
		"subpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		 0x040(%%rdx)	,%%xmm2\n\t"\
		"movaps		 0x050(%%rdx)	,%%xmm3\n\t"\
		"movaps		      (%%rdx)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm1\n\t"\
		"addpd		 0x050(%%rdx)	,%%xmm2\n\t"\
		"subpd		 0x040(%%rdx)	,%%xmm3\n\t"\
		"mulpd		      (%%rsi)	,%%xmm2\n\t"\
		"mulpd		      (%%rsi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"movq		%[__add0]		,%%rax\n\t"\
		"movq		%[__c02]		,%%rcx\n\t"\
		"subpd		%%xmm4		,%%xmm2\n\t"\
		"subpd		%%xmm5		,%%xmm3\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm2		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm5\n\t"\
		"movaps		%%xmm2		,      (%%rdx)\n\t"\
		"movaps		%%xmm3		, 0x010(%%rdx)\n\t"\
		"movaps		%%xmm4		,%%xmm2\n\t"\
		"movaps		%%xmm5		,%%xmm3\n\t"\
		"mulpd		      (%%rcx)	,%%xmm4\n\t"\
		"mulpd		      (%%rcx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movq		%[__add1]		,%%rbx\n\t"\
		"movaps		%%xmm5		,%%xmm3\n\t"\
		"movaps		%%xmm4		,%%xmm2\n\t"\
		"unpckhpd	 0x030(%%rbx)	,%%xmm3\n\t"\
		"unpcklpd	 0x030(%%rbx)	,%%xmm5\n\t"\
		"movaps		%%xmm3		, 0x030(%%rbx)\n\t"\
		"unpckhpd	 0x020(%%rbx)	,%%xmm2\n\t"\
		"unpcklpd	 0x020(%%rbx)	,%%xmm4\n\t"\
		"movaps		%%xmm2		, 0x020(%%rbx)\n\t"\
		"movaps		%%xmm5		, 0x030(%%rax)\n\t"\
		"movaps		%%xmm4		, 0x020(%%rax)\n\t"\
		"movq		%[__c12]		,%%rcx\n\t"\
		"movaps		      (%%rdx)	,%%xmm4\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm5\n\t"\
		"movaps		%%xmm4		,%%xmm2\n\t"\
		"movaps		%%xmm5		,%%xmm3\n\t"\
		"mulpd		      (%%rcx)	,%%xmm4\n\t"\
		"mulpd		      (%%rcx)	,%%xmm5\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm2\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm3\n\t"\
		"movaps		%%xmm4		,%%xmm2\n\t"\
		"unpckhpd	 0x130(%%rbx)	,%%xmm3\n\t"\
		"unpcklpd	 0x130(%%rbx)	,%%xmm5\n\t"\
		"movaps		%%xmm3		, 0x130(%%rbx)\n\t"\
		"unpckhpd	 0x120(%%rbx)	,%%xmm2\n\t"\
		"unpcklpd	 0x120(%%rbx)	,%%xmm4\n\t"\
		"movaps		%%xmm2		, 0x120(%%rbx)\n\t"\
		"movaps		%%xmm5		, 0x130(%%rax)\n\t"\
		"movaps		%%xmm4		, 0x120(%%rax)\n\t"\
		"movq		%[__c0A]		,%%rcx\n\t"\
		"subpd		%%xmm7		,%%xmm0\n\t"\
		"subpd		%%xmm6		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movaps		%%xmm7		,%%xmm4\n\t"\
		"movaps		%%xmm1		,%%xmm5\n\t"\
		"mulpd		      (%%rcx)	,%%xmm7\n\t"\
		"mulpd		      (%%rcx)	,%%xmm1\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm5\n\t"\
		"subpd		%%xmm4		,%%xmm1\n\t"\
		"addpd		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm7		,%%xmm4\n\t"\
		"unpckhpd	 0x0b0(%%rbx)	,%%xmm5\n\t"\
		"unpcklpd	 0x0b0(%%rbx)	,%%xmm1\n\t"\
		"movaps		%%xmm5		, 0x0b0(%%rbx)\n\t"\
		"unpckhpd	 0x0a0(%%rbx)	,%%xmm4\n\t"\
		"unpcklpd	 0x0a0(%%rbx)	,%%xmm7\n\t"\
		"movaps		%%xmm4		, 0x0a0(%%rbx)\n\t"\
		"movaps		%%xmm1		, 0x0b0(%%rax)\n\t"\
		"movaps		%%xmm7		, 0x0a0(%%rax)\n\t"\
		"movq		%[__c1A]		,%%rcx\n\t"\
		"movaps		%%xmm0		,%%xmm4\n\t"\
		"movaps		%%xmm6		,%%xmm5\n\t"\
		"mulpd		      (%%rcx)	,%%xmm0\n\t"\
		"mulpd		      (%%rcx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm5\n\t"\
		"subpd		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm5		,%%xmm0\n\t"\
		"movaps		%%xmm6		,%%xmm5\n\t"\
		"movaps		%%xmm0		,%%xmm4\n\t"\
		"unpckhpd	 0x1b0(%%rbx)	,%%xmm5\n\t"\
		"unpcklpd	 0x1b0(%%rbx)	,%%xmm6\n\t"\
		"movaps		%%xmm5		, 0x1b0(%%rbx)\n\t"\
		"unpckhpd	 0x1a0(%%rbx)	,%%xmm4\n\t"\
		"unpcklpd	 0x1a0(%%rbx)	,%%xmm0\n\t"\
		"movaps		%%xmm4		, 0x1a0(%%rbx)\n\t"\
		"movaps		%%xmm6		, 0x1b0(%%rax)\n\t"\
		"movaps		%%xmm0		, 0x1a0(%%rax)\n\t"\
		"/************************************************/\n\t"\
		"/*...Block 7: t0C,t1C,t2C,t3C -> r28,2C,2A,2E	*/\n\t"\
		"/************************************************/\n\t"\
		"movq		%[__r28]		,%%rdx\n\t"\
		"movq		%%rsi			,%%rcx\n\t"\
		"addq		$0x010			,%%rcx	/* cc0 */\n\t"\
		"movaps		 0x020(%%rdx)	,%%xmm4\n\t"\
		"movaps		 0x060(%%rdx)	,%%xmm0\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm5\n\t"\
		"movaps		 0x070(%%rdx)	,%%xmm1\n\t"\
		"movaps		 0x020(%%rdx)	,%%xmm6\n\t"\
		"movaps		 0x060(%%rdx)	,%%xmm2\n\t"\
		"movaps		 0x030(%%rdx)	,%%xmm7\n\t"\
		"movaps		 0x070(%%rdx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm4\n\t"\
		"mulpd		      (%%rcx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm5\n\t"\
		"mulpd		      (%%rcx)	,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm6\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm2\n\t"\
		"mulpd		      (%%rcx)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm3\n\t"\
		"subpd		%%xmm6		,%%xmm5\n\t"\
		"subpd		%%xmm2		,%%xmm1\n\t"\
		"addpd		%%xmm7		,%%xmm4\n\t"\
		"addpd		%%xmm3		,%%xmm0\n\t"\
		"movaps		%%xmm5		,%%xmm7\n\t"\
		"movaps		%%xmm4		,%%xmm6\n\t"\
		"addpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"subpd		%%xmm0		,%%xmm6\n\t"\
		"subpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		 0x040(%%rdx)	,%%xmm2\n\t"\
		"movaps		 0x050(%%rdx)	,%%xmm3\n\t"\
		"movaps		      (%%rdx)	,%%xmm0\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm1\n\t"\
		"subpd		 0x050(%%rdx)	,%%xmm2\n\t"\
		"addpd		 0x040(%%rdx)	,%%xmm3\n\t"\
		"mulpd		      (%%rsi)	,%%xmm2\n\t"\
		"mulpd		      (%%rsi)	,%%xmm3\n\t"\
		"subpd		%%xmm2		,%%xmm0\n\t"\
		"subpd		%%xmm3		,%%xmm1\n\t"\
		"addpd		%%xmm2		,%%xmm2\n\t"\
		"addpd		%%xmm3		,%%xmm3\n\t"\
		"addpd		%%xmm0		,%%xmm2\n\t"\
		"addpd		%%xmm1		,%%xmm3\n\t"\
		"movq		%[__add0]		,%%rax\n\t"\
		"movq		%[__c06]		,%%rcx\n\t"\
		"subpd		%%xmm6		,%%xmm0\n\t"\
		"subpd		%%xmm7		,%%xmm1\n\t"\
		"addpd		%%xmm6		,%%xmm6\n\t"\
		"addpd		%%xmm7		,%%xmm7\n\t"\
		"addpd		%%xmm0		,%%xmm6\n\t"\
		"addpd		%%xmm1		,%%xmm7\n\t"\
		"movaps		%%xmm0		,      (%%rdx)\n\t"\
		"movaps		%%xmm1		, 0x010(%%rdx)\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm6\n\t"\
		"mulpd		      (%%rcx)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movq		%[__add1]		,%%rbx\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"unpckhpd	 0x070(%%rbx)	,%%xmm1\n\t"\
		"unpcklpd	 0x070(%%rbx)	,%%xmm7\n\t"\
		"movaps		%%xmm1		, 0x070(%%rbx)\n\t"\
		"unpckhpd	 0x060(%%rbx)	,%%xmm0\n\t"\
		"unpcklpd	 0x060(%%rbx)	,%%xmm6\n\t"\
		"movaps		%%xmm0		, 0x060(%%rbx)\n\t"\
		"movaps		%%xmm7		, 0x070(%%rax)\n\t"\
		"movaps		%%xmm6		, 0x060(%%rax)\n\t"\
		"movq		%[__c16]		,%%rcx\n\t"\
		"movaps		      (%%rdx)	,%%xmm6\n\t"\
		"movaps		 0x010(%%rdx)	,%%xmm7\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm6\n\t"\
		"mulpd		      (%%rcx)	,%%xmm7\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm7\n\t"\
		"addpd		%%xmm1		,%%xmm6\n\t"\
		"movaps		%%xmm7		,%%xmm1\n\t"\
		"movaps		%%xmm6		,%%xmm0\n\t"\
		"unpckhpd	 0x170(%%rbx)	,%%xmm1\n\t"\
		"unpcklpd	 0x170(%%rbx)	,%%xmm7\n\t"\
		"movaps		%%xmm1		, 0x170(%%rbx)\n\t"\
		"unpckhpd	 0x160(%%rbx)	,%%xmm0\n\t"\
		"unpcklpd	 0x160(%%rbx)	,%%xmm6\n\t"\
		"movaps		%%xmm0		, 0x160(%%rbx)\n\t"\
		"movaps		%%xmm7		, 0x170(%%rax)\n\t"\
		"movaps		%%xmm6		, 0x160(%%rax)\n\t"\
		"movq		%[__c0E]		,%%rcx\n\t"\
		"subpd		%%xmm5		,%%xmm2\n\t"\
		"subpd		%%xmm4		,%%xmm3\n\t"\
		"addpd		%%xmm5		,%%xmm5\n\t"\
		"addpd		%%xmm4		,%%xmm4\n\t"\
		"addpd		%%xmm2		,%%xmm5\n\t"\
		"addpd		%%xmm3		,%%xmm4\n\t"\
		"movaps		%%xmm5		,%%xmm0\n\t"\
		"movaps		%%xmm3		,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm5\n\t"\
		"mulpd		      (%%rcx)	,%%xmm3\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm3\n\t"\
		"addpd		%%xmm1		,%%xmm5\n\t"\
		"movaps		%%xmm3		,%%xmm1\n\t"\
		"movaps		%%xmm5		,%%xmm0\n\t"\
		"unpckhpd	 0x0f0(%%rbx)	,%%xmm1\n\t"\
		"unpcklpd	 0x0f0(%%rbx)	,%%xmm3\n\t"\
		"movaps		%%xmm1		, 0x0f0(%%rbx)\n\t"\
		"unpckhpd	 0x0e0(%%rbx)	,%%xmm0\n\t"\
		"unpcklpd	 0x0e0(%%rbx)	,%%xmm5\n\t"\
		"movaps		%%xmm0		, 0x0e0(%%rbx)\n\t"\
		"movaps		%%xmm3		, 0x0f0(%%rax)\n\t"\
		"movaps		%%xmm5		, 0x0e0(%%rax)\n\t"\
		"movq		%[__c1E]		,%%rcx\n\t"\
		"movaps		%%xmm2		,%%xmm0\n\t"\
		"movaps		%%xmm4		,%%xmm1\n\t"\
		"mulpd		      (%%rcx)	,%%xmm2\n\t"\
		"mulpd		      (%%rcx)	,%%xmm4\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm0\n\t"\
		"mulpd		 0x010(%%rcx)	,%%xmm1\n\t"\
		"subpd		%%xmm0		,%%xmm4\n\t"\
		"addpd		%%xmm1		,%%xmm2\n\t"\
		"movaps		%%xmm4		,%%xmm1\n\t"\
		"movaps		%%xmm2		,%%xmm0\n\t"\
		"unpckhpd	 0x1f0(%%rbx)	,%%xmm1\n\t"\
		"unpcklpd	 0x1f0(%%rbx)	,%%xmm4\n\t"\
		"movaps		%%xmm1		, 0x1f0(%%rbx)\n\t"\
		"unpckhpd	 0x1e0(%%rbx)	,%%xmm0\n\t"\
		"unpcklpd	 0x1e0(%%rbx)	,%%xmm2\n\t"\
		"movaps		%%xmm0		, 0x1e0(%%rbx)\n\t"\
		"movaps		%%xmm4		, 0x1f0(%%rax)\n\t"\
		"movaps		%%xmm2		, 0x1e0(%%rax)\n\t"\
		:					/* outputs: none */\
		: [__add0 ] "m" (Xadd0)	/* All inputs from memory addresses here */\
		 ,[__add1] "m" (Xadd1)\
		 ,[__add2] "m" (Xadd2)\
		 ,[__isrt2] "m" (Xisrt2)\
		 ,[__r00] "m" (Xr00)\
		 ,[__r08] "m" (Xr08)\
		 ,[__r10] "m" (Xr10)\
		 ,[__r20] "m" (Xr20)\
		 ,[__r28] "m" (Xr28)\
		 ,[__r30] "m" (Xr30)\
		 ,[__c00] "m" (Xc00)\
		 ,[__c01] "m" (Xc01)\
		 ,[__c02] "m" (Xc02)\
		 ,[__c03] "m" (Xc03)\
		 ,[__c04] "m" (Xc04)\
		 ,[__c05] "m" (Xc05)\
		 ,[__c06] "m" (Xc06)\
		 ,[__c07] "m" (Xc07)\
		 ,[__c08] "m" (Xc08)\
		 ,[__c0A] "m" (Xc0A)\
		 ,[__c0C] "m" (Xc0C)\
		 ,[__c0E] "m" (Xc0E)\
		 ,[__c10] "m" (Xc10)\
		 ,[__c12] "m" (Xc12)\
		 ,[__c14] "m" (Xc14)\
		 ,[__c16] "m" (Xc16)\
		 ,[__c18] "m" (Xc18)\
		 ,[__c1A] "m" (Xc1A)\
		 ,[__c1C] "m" (Xc1C)\
		 ,[__c1E] "m" (Xc1E)\
		: "rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
	);\
	}

#endif	/* radix32_wrapper_square_gcc_h_included */

