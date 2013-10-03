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
#ifndef radix36_ditN_cy_dif1_gcc_h_included
#define radix36_ditN_cy_dif1_gcc_h_included

 #ifdef USE_AVX

	#define	SSE2_RADIX36_DIT_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xp28,Xp32,Xr00,Xout,Xcc1)\
	{\
	__asm__ volatile (\
		/*	add0,1,3,2 = a+p00+p[0,1,2,3]: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x240, 0x480, r00: */\
			"movq	%[__r00],%%rsi\n\t"\
			"movq	%[__add],%%rax\n\t"\
			"movslq	%[__p01],%%rbx\n\t"\
			"movslq	%[__p02],%%rcx\n\t"\
			"movslq	%[__p03],%%rdx\n\t"\
			"shlq	$3,%%rbx		/* Pointer offset for floating doubles */\n\t"\
			"shlq	$3,%%rcx\n\t"\
			"shlq	$3,%%rdx\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* ecx <-> edx */\n\t"\
			"vmovaps	     (%%rax),%%ymm2	\n\t"\
			"vmovaps	     (%%rdx),%%ymm6	\n\t"\
			"vmovaps	0x020(%%rax),%%ymm3	\n\t"\
			"vmovaps	0x020(%%rdx),%%ymm7	\n\t"\
			"vmovaps	     (%%rbx),%%ymm0	\n\t"\
			"vmovaps	     (%%rcx),%%ymm4	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t"\
			"vmovaps	0x020(%%rcx),%%ymm5	\n\t"\
			"vsubpd	%%ymm0,%%ymm2,%%ymm2	\n\t"\
			"vsubpd	%%ymm4,%%ymm6,%%ymm6	\n\t"\
			"vsubpd	%%ymm1,%%ymm3,%%ymm3	\n\t"\
			"vsubpd	%%ymm5,%%ymm7,%%ymm7	\n\t"\
			"vaddpd	%%ymm0,%%ymm0,%%ymm0	\n\t"		/*	add2,3,0,1 = a+p04+p[0,1,2,3]: cdab, r04: */\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t			movslq	%[__p04],%%rdi\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1	\n\t			shlq	$3,%%rdi \n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t			addq	%%rdi,%%rax	/* add04 */\n\t"\
			"vaddpd	%%ymm2,%%ymm0,%%ymm0	\n\t			addq	%%rdi,%%rbx\n\t"\
			"vaddpd	%%ymm6,%%ymm4,%%ymm4	\n\t			addq	%%rdi,%%rcx\n\t"\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1	\n\t			addq	%%rdi,%%rdx\n\t"\
			"vaddpd	%%ymm7,%%ymm5,%%ymm5	\n\t			/* e[ab]x <-> e[cd]x */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0	\n\t			vmovaps	     (%%rcx),%%ymm10	\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2	\n\t			vmovaps	     (%%rax),%%ymm14	\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1	\n\t			vmovaps	0x020(%%rcx),%%ymm11	\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3	\n\t			vmovaps	0x020(%%rax),%%ymm15	\n\t"\
			"vmovaps	%%ymm0,0x480(%%rsi)	\n\t			vmovaps	     (%%rdx),%%ymm8 	\n\t"\
			"vmovaps	%%ymm2,0x6c0(%%rsi)	\n\t			vmovaps	     (%%rbx),%%ymm12	\n\t"\
			"vmovaps	%%ymm1,0x4a0(%%rsi)	\n\t			vmovaps	0x020(%%rdx),%%ymm9 	\n\t"\
			"vmovaps	%%ymm3,0x260(%%rsi)	\n\t			vmovaps	0x020(%%rbx),%%ymm13	\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t			vsubpd	%%ymm8 ,%%ymm10,%%ymm10	\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t			vsubpd	%%ymm12,%%ymm14,%%ymm14	\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t			vsubpd	%%ymm9 ,%%ymm11,%%ymm11	\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t			vsubpd	%%ymm13,%%ymm15,%%ymm15	\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4	\n\t			vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 	\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7	\n\t			vaddpd	%%ymm12,%%ymm12,%%ymm12	\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5	\n\t			vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 	\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6	\n\t			vaddpd	%%ymm13,%%ymm13,%%ymm13	\n\t"\
			"vmovaps	%%ymm4,     (%%rsi)	\n\t			vaddpd	%%ymm10,%%ymm8 ,%%ymm8 	\n\t"\
			"vmovaps	%%ymm7,0x240(%%rsi)	\n\t			vaddpd	%%ymm14,%%ymm12,%%ymm12	\n\t"\
			"vmovaps	%%ymm5,0x020(%%rsi)	\n\t			vaddpd	%%ymm11,%%ymm9 ,%%ymm9 	\n\t"\
			"vmovaps	%%ymm6,0x6e0(%%rsi)	\n\t			vaddpd	%%ymm15,%%ymm13,%%ymm13	\n\t"\
		/*	add1,0,2,3 = a+p08+p[0,1,2,3]: bacd, r02: */"	addq	$0x080,%%rsi	/* r04 */\n\t"\
			"addq	%%rdi,%%rax	/* add08 */\n\t				vsubpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
			"addq	%%rdi,%%rbx\n\t							vsubpd	%%ymm15,%%ymm10,%%ymm10	\n\t"\
			"addq	%%rdi,%%rcx\n\t							vsubpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
			"addq	%%rdi,%%rdx\n\t							vsubpd	%%ymm14,%%ymm11,%%ymm11	\n\t"\
			"/* eax <-> ebx */\n\t							vmovaps	%%ymm8 ,0x480(%%rsi)	\n\t"\
			"vmovaps	     (%%rbx),%%ymm2	\n\t			vmovaps	%%ymm10,0x6c0(%%rsi)	\n\t"\
			"vmovaps	     (%%rcx),%%ymm6	\n\t			vmovaps	%%ymm9 ,0x4a0(%%rsi)	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm3	\n\t			vmovaps	%%ymm11,0x260(%%rsi)	\n\t"\
			"vmovaps	0x020(%%rcx),%%ymm7	\n\t			vaddpd	%%ymm12,%%ymm12,%%ymm12	\n\t"\
			"vmovaps	     (%%rax),%%ymm0	\n\t			vaddpd	%%ymm15,%%ymm15,%%ymm15	\n\t"\
			"vmovaps	     (%%rdx),%%ymm4	\n\t			vaddpd	%%ymm13,%%ymm13,%%ymm13	\n\t"\
			"vmovaps	0x020(%%rax),%%ymm1	\n\t			vaddpd	%%ymm14,%%ymm14,%%ymm14	\n\t"\
			"vmovaps	0x020(%%rdx),%%ymm5	\n\t			vaddpd	%%ymm8 ,%%ymm12,%%ymm12	\n\t"\
			"vsubpd	%%ymm0,%%ymm2,%%ymm2	\n\t			vaddpd	%%ymm10,%%ymm15,%%ymm15	\n\t"\
			"vsubpd	%%ymm4,%%ymm6,%%ymm6	\n\t			vaddpd	%%ymm9 ,%%ymm13,%%ymm13	\n\t"\
			"vsubpd	%%ymm1,%%ymm3,%%ymm3	\n\t			vaddpd	%%ymm11,%%ymm14,%%ymm14	\n\t"\
			"vsubpd	%%ymm5,%%ymm7,%%ymm7	\n\t			vmovaps	%%ymm12,     (%%rsi)	\n\t"\
			"vaddpd	%%ymm0,%%ymm0,%%ymm0	\n\t			vmovaps	%%ymm15,0x240(%%rsi)	\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t			vmovaps	%%ymm13,0x020(%%rsi)	\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1	\n\t			vmovaps	%%ymm14,0x6e0(%%rsi)	\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"		/*	add3,2,1,0 = a+p12+p[0,1,2,3]: dcba, r0c: */\
			"vaddpd	%%ymm2,%%ymm0,%%ymm0	\n\t			addq	%%rdi,%%rax	/* add12 */\n\t"\
			"vaddpd	%%ymm6,%%ymm4,%%ymm4	\n\t			addq	%%rdi,%%rbx\n\t"\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1	\n\t			addq	%%rdi,%%rcx\n\t"\
			"vaddpd	%%ymm7,%%ymm5,%%ymm5	\n\t			addq	%%rdi,%%rdx\n\t"\
			"subq	$0x040,%%rsi	/* r02 */\n\t			/* eax <-> edx, ebx <-> ecx */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0	\n\t			vmovaps	     (%%rdx),%%ymm10	\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2	\n\t			vmovaps	     (%%rbx),%%ymm14	\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1	\n\t			vmovaps	0x020(%%rdx),%%ymm11	\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3	\n\t			vmovaps	0x020(%%rbx),%%ymm15	\n\t"\
			"vmovaps	%%ymm0,0x480(%%rsi)	\n\t			vmovaps	     (%%rcx),%%ymm8 	\n\t"\
			"vmovaps	%%ymm2,0x6c0(%%rsi)	\n\t			vmovaps	     (%%rax),%%ymm12	\n\t"\
			"vmovaps	%%ymm1,0x4a0(%%rsi)	\n\t			vmovaps	0x020(%%rcx),%%ymm9 	\n\t"\
			"vmovaps	%%ymm3,0x260(%%rsi)	\n\t			vmovaps	0x020(%%rax),%%ymm13	\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t			vsubpd	%%ymm8 ,%%ymm10,%%ymm10	\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t			vsubpd	%%ymm12,%%ymm14,%%ymm14	\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t			vsubpd	%%ymm9 ,%%ymm11,%%ymm11	\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t			vsubpd	%%ymm13,%%ymm15,%%ymm15	\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4	\n\t			vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 	\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7	\n\t			vaddpd	%%ymm12,%%ymm12,%%ymm12	\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5	\n\t			vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 	\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6	\n\t			vaddpd	%%ymm13,%%ymm13,%%ymm13	\n\t"\
			"vmovaps	%%ymm4,     (%%rsi)	\n\t			vaddpd	%%ymm10,%%ymm8 ,%%ymm8 	\n\t"\
			"vmovaps	%%ymm7,0x240(%%rsi)	\n\t			vaddpd	%%ymm14,%%ymm12,%%ymm12	\n\t"\
			"vmovaps	%%ymm5,0x020(%%rsi)	\n\t			vaddpd	%%ymm11,%%ymm9 ,%%ymm9 	\n\t"\
			"vmovaps	%%ymm6,0x6e0(%%rsi)	\n\t			vaddpd	%%ymm15,%%ymm13,%%ymm13	\n\t"\
		/*	add0,1,3,2 = a+p16+p[0,1,2,3]: abdc, r0g: */"	addq	$0x140,%%rsi	/* r0c */\n\t"\
			"addq	%%rdi,%%rax	/* add16 */\n\t				vsubpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
			"addq	%%rdi,%%rbx\n\t							vsubpd	%%ymm15,%%ymm10,%%ymm10	\n\t"\
			"addq	%%rdi,%%rcx\n\t							vsubpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
			"addq	%%rdi,%%rdx\n\t							vsubpd	%%ymm14,%%ymm11,%%ymm11	\n\t"\
			"/* ecx <-> edx */\n\t							vmovaps	%%ymm8 ,0x480(%%rsi)	\n\t"\
			"vmovaps	     (%%rax),%%ymm2	\n\t			vmovaps	%%ymm10,0x6c0(%%rsi)	\n\t"\
			"vmovaps	     (%%rdx),%%ymm6	\n\t			vmovaps	%%ymm9 ,0x4a0(%%rsi)	\n\t"\
			"vmovaps	0x020(%%rax),%%ymm3	\n\t			vmovaps	%%ymm11,0x260(%%rsi)	\n\t"\
			"vmovaps	0x020(%%rdx),%%ymm7	\n\t			vaddpd	%%ymm12,%%ymm12,%%ymm12	\n\t"\
			"vmovaps	     (%%rbx),%%ymm0	\n\t			vaddpd	%%ymm15,%%ymm15,%%ymm15	\n\t"\
			"vmovaps	     (%%rcx),%%ymm4	\n\t			vaddpd	%%ymm13,%%ymm13,%%ymm13	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t			vaddpd	%%ymm14,%%ymm14,%%ymm14	\n\t"\
			"vmovaps	0x020(%%rcx),%%ymm5	\n\t			vaddpd	%%ymm8 ,%%ymm12,%%ymm12	\n\t"\
			"vsubpd	%%ymm0,%%ymm2,%%ymm2	\n\t			vaddpd	%%ymm10,%%ymm15,%%ymm15	\n\t"\
			"vsubpd	%%ymm4,%%ymm6,%%ymm6	\n\t			vaddpd	%%ymm9 ,%%ymm13,%%ymm13	\n\t"\
			"vsubpd	%%ymm1,%%ymm3,%%ymm3	\n\t			vaddpd	%%ymm11,%%ymm14,%%ymm14	\n\t"\
			"vsubpd	%%ymm5,%%ymm7,%%ymm7	\n\t			vmovaps	%%ymm12,     (%%rsi)	\n\t"\
			"vaddpd	%%ymm0,%%ymm0,%%ymm0	\n\t			vmovaps	%%ymm15,0x240(%%rsi)	\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t			vmovaps	%%ymm13,0x020(%%rsi)	\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1	\n\t			vmovaps	%%ymm14,0x6e0(%%rsi)	\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"		/*	add2,3,0,1 = a+p20+p[0,1,2,3]: cdab, r0e: */\
			"vaddpd	%%ymm2,%%ymm0,%%ymm0	\n\t			addq	%%rdi,%%rax	/* add20 */\n\t"\
			"vaddpd	%%ymm6,%%ymm4,%%ymm4	\n\t			addq	%%rdi,%%rbx\n\t"\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1	\n\t			addq	%%rdi,%%rcx\n\t"\
			"vaddpd	%%ymm7,%%ymm5,%%ymm5	\n\t			addq	%%rdi,%%rdx\n\t"\
			"addq	$0x080,%%rsi	/* r0g */\n\t			/* e[ab]x <-> e[cd]x */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0	\n\t			vmovaps	     (%%rcx),%%ymm10	\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2	\n\t			vmovaps	     (%%rax),%%ymm14	\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1	\n\t			vmovaps	0x020(%%rcx),%%ymm11	\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3	\n\t			vmovaps	0x020(%%rax),%%ymm15	\n\t"\
			"vmovaps	%%ymm0,0x480(%%rsi)	\n\t			vmovaps	     (%%rdx),%%ymm8 	\n\t"\
			"vmovaps	%%ymm2,0x6c0(%%rsi)	\n\t			vmovaps	     (%%rbx),%%ymm12	\n\t"\
			"vmovaps	%%ymm1,0x4a0(%%rsi)	\n\t			vmovaps	0x020(%%rdx),%%ymm9 	\n\t"\
			"vmovaps	%%ymm3,0x260(%%rsi)	\n\t			vmovaps	0x020(%%rbx),%%ymm13	\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t			vsubpd	%%ymm8 ,%%ymm10,%%ymm10	\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t			vsubpd	%%ymm12,%%ymm14,%%ymm14	\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t			vsubpd	%%ymm9 ,%%ymm11,%%ymm11	\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t			vsubpd	%%ymm13,%%ymm15,%%ymm15	\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4	\n\t			vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 	\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7	\n\t			vaddpd	%%ymm12,%%ymm12,%%ymm12	\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5	\n\t			vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 	\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6	\n\t			vaddpd	%%ymm13,%%ymm13,%%ymm13	\n\t"\
			"vmovaps	%%ymm4,     (%%rsi)	\n\t			vaddpd	%%ymm10,%%ymm8 ,%%ymm8 	\n\t"\
			"vmovaps	%%ymm7,0x240(%%rsi)	\n\t			vaddpd	%%ymm14,%%ymm12,%%ymm12	\n\t"\
			"vmovaps	%%ymm5,0x020(%%rsi)	\n\t			vaddpd	%%ymm11,%%ymm9 ,%%ymm9 	\n\t"\
			"vmovaps	%%ymm6,0x6e0(%%rsi)	\n\t			vaddpd	%%ymm15,%%ymm13,%%ymm13	\n\t"\
		/*	add1,0,2,3 = a+p24+p[0,1,2,3]: bacd, r08: */"	subq	$0x040,%%rsi	/* r0e */\n\t"\
			"addq	%%rdi,%%rax	/* add24 */\n\t				vsubpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
			"addq	%%rdi,%%rbx\n\t							vsubpd	%%ymm15,%%ymm10,%%ymm10	\n\t"\
			"addq	%%rdi,%%rcx\n\t							vsubpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
			"addq	%%rdi,%%rdx\n\t							vsubpd	%%ymm14,%%ymm11,%%ymm11	\n\t"\
			"/* eax <-> ebx */\n\t							vmovaps	%%ymm8 ,0x480(%%rsi)	\n\t"\
			"vmovaps	     (%%rbx),%%ymm2	\n\t			vmovaps	%%ymm10,0x6c0(%%rsi)	\n\t"\
			"vmovaps	     (%%rcx),%%ymm6	\n\t			vmovaps	%%ymm9 ,0x4a0(%%rsi)	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm3	\n\t			vmovaps	%%ymm11,0x260(%%rsi)	\n\t"\
			"vmovaps	0x020(%%rcx),%%ymm7	\n\t			vaddpd	%%ymm12,%%ymm12,%%ymm12	\n\t"\
			"vmovaps	     (%%rax),%%ymm0	\n\t			vaddpd	%%ymm15,%%ymm15,%%ymm15	\n\t"\
			"vmovaps	     (%%rdx),%%ymm4	\n\t			vaddpd	%%ymm13,%%ymm13,%%ymm13	\n\t"\
			"vmovaps	0x020(%%rax),%%ymm1	\n\t			vaddpd	%%ymm14,%%ymm14,%%ymm14	\n\t"\
			"vmovaps	0x020(%%rdx),%%ymm5	\n\t			vaddpd	%%ymm8 ,%%ymm12,%%ymm12	\n\t"\
			"vsubpd	%%ymm0,%%ymm2,%%ymm2	\n\t			vaddpd	%%ymm10,%%ymm15,%%ymm15	\n\t"\
			"vsubpd	%%ymm4,%%ymm6,%%ymm6	\n\t			vaddpd	%%ymm9 ,%%ymm13,%%ymm13	\n\t"\
			"vsubpd	%%ymm1,%%ymm3,%%ymm3	\n\t			vaddpd	%%ymm11,%%ymm14,%%ymm14	\n\t"\
			"vsubpd	%%ymm5,%%ymm7,%%ymm7	\n\t			vmovaps	%%ymm12,     (%%rsi)	\n\t"\
			"vaddpd	%%ymm0,%%ymm0,%%ymm0	\n\t			vmovaps	%%ymm15,0x240(%%rsi)	\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t			vmovaps	%%ymm13,0x020(%%rsi)	\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1	\n\t			vmovaps	%%ymm14,0x6e0(%%rsi)	\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"		/*	add3,2,1,0 = a+p28+p[0,1,2,3]: dcba, r06: */\
			"vaddpd	%%ymm2,%%ymm0,%%ymm0	\n\t			addq	%%rdi,%%rax	/* add28 */\n\t"\
			"vaddpd	%%ymm6,%%ymm4,%%ymm4	\n\t			addq	%%rdi,%%rbx\n\t"\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1	\n\t			addq	%%rdi,%%rcx\n\t"\
			"vaddpd	%%ymm7,%%ymm5,%%ymm5	\n\t			addq	%%rdi,%%rdx\n\t"\
			"subq	$0x0c0,%%rsi	/* r08 */\n\t			/* eax <-> edx, ebx <-> ecx */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0	\n\t			vmovaps	     (%%rdx),%%ymm10	\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2	\n\t			vmovaps	     (%%rbx),%%ymm14	\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1	\n\t			vmovaps	0x020(%%rdx),%%ymm11	\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3	\n\t			vmovaps	0x020(%%rbx),%%ymm15	\n\t"\
			"vmovaps	%%ymm0,0x480(%%rsi)	\n\t			vmovaps	     (%%rcx),%%ymm8 	\n\t"\
			"vmovaps	%%ymm2,0x6c0(%%rsi)	\n\t			vmovaps	     (%%rax),%%ymm12	\n\t"\
			"vmovaps	%%ymm1,0x4a0(%%rsi)	\n\t			vmovaps	0x020(%%rcx),%%ymm9 	\n\t"\
			"vmovaps	%%ymm3,0x260(%%rsi)	\n\t			vmovaps	0x020(%%rax),%%ymm13	\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t			vsubpd	%%ymm8 ,%%ymm10,%%ymm10	\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t			vsubpd	%%ymm12,%%ymm14,%%ymm14	\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t			vsubpd	%%ymm9 ,%%ymm11,%%ymm11	\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t			vsubpd	%%ymm13,%%ymm15,%%ymm15	\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4	\n\t			vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 	\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7	\n\t			vaddpd	%%ymm12,%%ymm12,%%ymm12	\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5	\n\t			vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 	\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6	\n\t			vaddpd	%%ymm13,%%ymm13,%%ymm13	\n\t"\
			"vmovaps	%%ymm4,     (%%rsi)	\n\t			vaddpd	%%ymm10,%%ymm8 ,%%ymm8 	\n\t"\
			"vmovaps	%%ymm7,0x240(%%rsi)	\n\t			vaddpd	%%ymm14,%%ymm12,%%ymm12	\n\t"\
			"vmovaps	%%ymm5,0x020(%%rsi)	\n\t			vaddpd	%%ymm11,%%ymm9 ,%%ymm9 	\n\t"\
			"vmovaps	%%ymm6,0x6e0(%%rsi)	\n\t			vaddpd	%%ymm15,%%ymm13,%%ymm13	\n\t"\
		/*	add0,1,3,2 = a+p32+p[0,1,2,3]: abdc, r0a: */"	subq	$0x040,%%rsi	/* r04 */\n\t"\
			"addq	%%rdi,%%rax	/* add32 */\n\t				vsubpd	%%ymm12,%%ymm8 ,%%ymm8 	\n\t"\
			"addq	%%rdi,%%rbx\n\t							vsubpd	%%ymm15,%%ymm10,%%ymm10	\n\t"\
			"addq	%%rdi,%%rcx\n\t							vsubpd	%%ymm13,%%ymm9 ,%%ymm9 	\n\t"\
			"addq	%%rdi,%%rdx\n\t							vsubpd	%%ymm14,%%ymm11,%%ymm11	\n\t"\
			"/* ecx <-> edx */\n\t							vmovaps	%%ymm8 ,0x480(%%rsi)	\n\t"\
			"vmovaps	     (%%rax),%%ymm2	\n\t			vmovaps	%%ymm10,0x6c0(%%rsi)	\n\t"\
			"vmovaps	     (%%rdx),%%ymm6	\n\t			vmovaps	%%ymm9 ,0x4a0(%%rsi)	\n\t"\
			"vmovaps	0x020(%%rax),%%ymm3	\n\t			vmovaps	%%ymm11,0x260(%%rsi)	\n\t"\
			"vmovaps	0x020(%%rdx),%%ymm7	\n\t			vaddpd	%%ymm12,%%ymm12,%%ymm12	\n\t"\
			"vmovaps	     (%%rbx),%%ymm0	\n\t			vaddpd	%%ymm15,%%ymm15,%%ymm15	\n\t"\
			"vmovaps	     (%%rcx),%%ymm4	\n\t			vaddpd	%%ymm13,%%ymm13,%%ymm13	\n\t"\
			"vmovaps	0x020(%%rbx),%%ymm1	\n\t			vaddpd	%%ymm14,%%ymm14,%%ymm14	\n\t"\
			"vmovaps	0x020(%%rcx),%%ymm5	\n\t			vaddpd	%%ymm8 ,%%ymm12,%%ymm12	\n\t"\
			"vsubpd	%%ymm0,%%ymm2,%%ymm2	\n\t			vaddpd	%%ymm10,%%ymm15,%%ymm15	\n\t"\
			"vsubpd	%%ymm4,%%ymm6,%%ymm6	\n\t			vaddpd	%%ymm9 ,%%ymm13,%%ymm13	\n\t"\
			"vsubpd	%%ymm1,%%ymm3,%%ymm3	\n\t			vaddpd	%%ymm11,%%ymm14,%%ymm14	\n\t"\
			"vsubpd	%%ymm5,%%ymm7,%%ymm7	\n\t			vmovaps	%%ymm12,     (%%rsi)	\n\t"\
			"vaddpd	%%ymm0,%%ymm0,%%ymm0	\n\t			vmovaps	%%ymm15,0x240(%%rsi)	\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t			vmovaps	%%ymm13,0x020(%%rsi)	\n\t"\
			"vaddpd	%%ymm1,%%ymm1,%%ymm1	\n\t			vmovaps	%%ymm14,0x6e0(%%rsi)	\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
			"vaddpd	%%ymm6,%%ymm4,%%ymm4	\n\t"\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
			"vaddpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
			"addq	$0x080,%%rsi	/* r04 */\n\t"\
			"vsubpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
			"vsubpd	%%ymm7,%%ymm2,%%ymm2	\n\t"\
			"vsubpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
			"vsubpd	%%ymm6,%%ymm3,%%ymm3	\n\t"\
			"vmovaps	%%ymm0,0x480(%%rsi)	\n\t"\
			"vmovaps	%%ymm2,0x6c0(%%rsi)	\n\t"\
			"vmovaps	%%ymm1,0x4a0(%%rsi)	\n\t"\
			"vmovaps	%%ymm3,0x260(%%rsi)	\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4	\n\t"\
			"vaddpd	%%ymm7,%%ymm7,%%ymm7	\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm6,%%ymm6,%%ymm6	\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4	\n\t"\
			"vaddpd	%%ymm2,%%ymm7,%%ymm7	\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5	\n\t"\
			"vaddpd	%%ymm3,%%ymm6,%%ymm6	\n\t"\
			"vmovaps	%%ymm4,     (%%rsi)	\n\t"\
			"vmovaps	%%ymm7,0x240(%%rsi)	\n\t"\
			"vmovaps	%%ymm5,0x020(%%rsi)	\n\t"\
			"vmovaps	%%ymm6,0x6e0(%%rsi)	\n\t"\
			"\n\t"\
			"/***************************************/\n\t"\
			"/*...and now do 4 radix-9 transforms...*/\n\t"\
			"/***************************************/\n\t"\
			"\n\t"\
			/* SSE2_RADIX_09_DIT_0TWIDDLE(r00,s1p[00,32,28,24,20,16,12,08,04]r): */\
			"movq	%[__out],%%rsi 	/* __o0-8: rsi,di store o-addresses */\n\t"\
			"movq	%%rsi		,%%rdi\n\t"\
			"movq	%[__r00]	,%%rax\n\t"/* rcol locals based on r10 = r00+0x240 */\
			"movq	%[__cc1],%%rdx 	/* rdx stores trig addresses throughout */\n\t"\
			"movq	%%rax		,%%rbx\n\t"						/* SSE2_RADIX_09_DIT_0TWIDDLE(r10,s1p[27,23,19,15,11,07,03,35,31]r): */\
			"movq	%%rax		,%%rcx\n\t						/* rcol: r10,r11 store o-addresses */\n\t"\
			"addq	$0x40		,%%rbx\n\t						leaq	0x6c0(%%rsi),%%r10 		/* s1p27 */\n\t"\
			"addq	$0x80		,%%rcx\n\t						movq	%%r10		,%%r11\n\t"\
			"addq	$0x80		,%%rdx 		/* c3m1 */\n\t"\
			"vmovaps	    (%%rbx)	,%%ymm2	\n\t				vmovaps	0x240(%%rbx)	,%%ymm10\n\t"\
			"vmovaps	0x20(%%rbx)	,%%ymm3	\n\t				vmovaps	0x260(%%rbx)	,%%ymm11\n\t"\
			"vmovaps	    (%%rax)	,%%ymm0	\n\t				vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
			"vmovaps	0x20(%%rax)	,%%ymm1	\n\t				vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
			"vmovaps	    (%%rcx)	,%%ymm6	\n\t				vmovaps	0x240(%%rcx)	,%%ymm14\n\t"\
			"vmovaps	0x20(%%rcx)	,%%ymm7	\n\t				vmovaps	0x260(%%rcx)	,%%ymm15\n\t"\
			"vmovaps	%%ymm2		,%%ymm4	\n\t				vmovaps	%%ymm10		,%%ymm12\n\t"\
			"vmovaps	%%ymm3		,%%ymm5	\n\t				vmovaps	%%ymm11		,%%ymm13\n\t"\
			"vaddpd	%%ymm6,%%ymm2,%%ymm2	\n\t				vaddpd	%%ymm14,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm7,%%ymm3,%%ymm3	\n\t				vaddpd	%%ymm15,%%ymm11,%%ymm11\n\t"\
			"vsubpd	%%ymm6,%%ymm4,%%ymm4	\n\t				vsubpd	%%ymm14,%%ymm12,%%ymm12\n\t"\
			"vsubpd	%%ymm7,%%ymm5,%%ymm5	\n\t				vsubpd	%%ymm15,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm2,%%ymm0,%%ymm0	\n\t				vaddpd	%%ymm10,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1	\n\t				vaddpd	%%ymm11,%%ymm9 ,%%ymm9 \n\t"\
			/* These are shared between lcol & rcol: */\
			"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
			"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
			"vmovaps	%%ymm0		,    (%%rax)\n\t			vmovaps	%%ymm8 		,0x240(%%rax)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rax)\n\t			vmovaps	%%ymm9 		,0x260(%%rax)\n\t"\
			"vmulpd	%%ymm6,%%ymm2,%%ymm2	\n\t				vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
			"vmulpd	%%ymm6,%%ymm3,%%ymm3	\n\t				vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
			"vmulpd	%%ymm7,%%ymm4,%%ymm4	\n\t				vmulpd	%%ymm7 ,%%ymm12,%%ymm12\n\t"\
			"vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t				vmulpd	%%ymm7 ,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm0,%%ymm2,%%ymm2	\n\t				vaddpd	%%ymm8 ,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm1,%%ymm3,%%ymm3	\n\t				vaddpd	%%ymm9 ,%%ymm11,%%ymm11\n\t"\
			"vmovaps	%%ymm2		,%%ymm0\n\t					vmovaps	%%ymm10		,%%ymm8 \n\t"\
			"vmovaps	%%ymm3		,%%ymm1\n\t					vmovaps	%%ymm11		,%%ymm9 \n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2	\n\t				vaddpd	%%ymm13,%%ymm10,%%ymm10\n\t"\
			"vsubpd	%%ymm4,%%ymm3,%%ymm3	\n\t				vsubpd	%%ymm12,%%ymm11,%%ymm11\n\t"\
			"vsubpd	%%ymm5,%%ymm0,%%ymm0	\n\t				vsubpd	%%ymm13,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm4,%%ymm1,%%ymm1	\n\t				vaddpd	%%ymm12,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm2		,    (%%rbx)\n\t			vmovaps	%%ymm10		,0x240(%%rbx)\n\t"\
			"vmovaps	%%ymm3		,0x20(%%rbx)\n\t			vmovaps	%%ymm11		,0x260(%%rbx)\n\t"\
			"vmovaps	%%ymm0		,    (%%rcx)\n\t			vmovaps	%%ymm8 		,0x240(%%rcx)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rcx)\n\t			vmovaps	%%ymm9 		,0x260(%%rcx)\n\t"\
			"addq	$0xc0		,%%rax\n\t"\
			"addq	$0xc0		,%%rbx\n\t"\
			"addq	$0xc0		,%%rcx\n\t"\
			"vmovaps	    (%%rbx)	,%%ymm4\n\t					vmovaps	0x240(%%rbx)	,%%ymm12\n\t"\
			"vmovaps	0x20(%%rbx)	,%%ymm5\n\t					vmovaps	0x260(%%rbx)	,%%ymm13\n\t"\
			"vmovaps	    (%%rcx)	,%%ymm2\n\t					vmovaps	0x240(%%rcx)	,%%ymm10\n\t"\
			"vmovaps	0x20(%%rcx)	,%%ymm3\n\t					vmovaps	0x260(%%rcx)	,%%ymm11\n\t"\
			"vmovaps	    (%%rax)	,%%ymm0\n\t					vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
			"vmovaps	0x20(%%rax)	,%%ymm1\n\t					vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
			"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t					vsubpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
			"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t					vsubpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm10,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm11,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t					vaddpd	%%ymm10,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm11,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm0		,    (%%rax)\n\t			vmovaps	%%ymm8 		,0x240(%%rax)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rax)\n\t			vmovaps	%%ymm9 		,0x260(%%rax)\n\t"\
			"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t					vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
			"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t					vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
			"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t					vmulpd	%%ymm7 ,%%ymm12,%%ymm12\n\t"\
			"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t					vmulpd	%%ymm7 ,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm8 ,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm9 ,%%ymm11,%%ymm11\n\t"\
			"vmovaps	%%ymm2		,%%ymm0\n\t					vmovaps	%%ymm10		,%%ymm8 \n\t"\
			"vmovaps	%%ymm3		,%%ymm1\n\t					vmovaps	%%ymm11		,%%ymm9 \n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm13,%%ymm10,%%ymm10\n\t"\
			"vsubpd	%%ymm4,%%ymm3,%%ymm3\n\t					vsubpd	%%ymm12,%%ymm11,%%ymm11\n\t"\
			"vsubpd	%%ymm5,%%ymm0,%%ymm0\n\t					vsubpd	%%ymm13,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm4,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm12,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm2		,    (%%rbx)\n\t			vmovaps	%%ymm10		,0x240(%%rbx)\n\t"\
			"vmovaps	%%ymm3		,0x20(%%rbx)\n\t			vmovaps	%%ymm11		,0x260(%%rbx)\n\t"\
			"vmovaps	%%ymm0		,    (%%rcx)\n\t			vmovaps	%%ymm8 		,0x240(%%rcx)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rcx)\n\t			vmovaps	%%ymm9 		,0x260(%%rcx)\n\t"\
			"addq	$0xc0		,%%rax\n\t"\
			"addq	$0xc0		,%%rbx\n\t"\
			"addq	$0xc0		,%%rcx\n\t"\
			"vmovaps	    (%%rbx)	,%%ymm4\n\t					vmovaps	0x240(%%rbx)	,%%ymm12\n\t"\
			"vmovaps	0x20(%%rbx)	,%%ymm5\n\t					vmovaps	0x260(%%rbx)	,%%ymm13\n\t"\
			"vmovaps	    (%%rcx)	,%%ymm2\n\t					vmovaps	0x240(%%rcx)	,%%ymm10\n\t"\
			"vmovaps	0x20(%%rcx)	,%%ymm3\n\t					vmovaps	0x260(%%rcx)	,%%ymm11\n\t"\
			"vmovaps	    (%%rax)	,%%ymm0\n\t					vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
			"vmovaps	0x20(%%rax)	,%%ymm1\n\t					vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
			"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t					vsubpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
			"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t					vsubpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm10,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm11,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t					vaddpd	%%ymm10,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm11,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm0		,    (%%rax)\n\t			vmovaps	%%ymm8 		,0x240(%%rax)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rax)\n\t			vmovaps	%%ymm9 		,0x260(%%rax)\n\t"\
			"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t					vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
			"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t					vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
			"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t					vmulpd	%%ymm7 ,%%ymm12,%%ymm12\n\t"\
			"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t					vmulpd	%%ymm7 ,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm8 ,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm9 ,%%ymm11,%%ymm11\n\t"\
			"vmovaps	%%ymm2		,%%ymm0\n\t					vmovaps	%%ymm10		,%%ymm8 \n\t"\
			"vmovaps	%%ymm3		,%%ymm1\n\t					vmovaps	%%ymm11		,%%ymm9 \n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm13,%%ymm10,%%ymm10\n\t"\
			"vsubpd	%%ymm4,%%ymm3,%%ymm3\n\t					vsubpd	%%ymm12,%%ymm11,%%ymm11\n\t"\
			"vsubpd	%%ymm5,%%ymm0,%%ymm0\n\t					vsubpd	%%ymm13,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm4,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm12,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm2		,    (%%rbx)\n\t			vmovaps	%%ymm10		,0x240(%%rbx)\n\t"\
			"vmovaps	%%ymm3		,0x20(%%rbx)\n\t			vmovaps	%%ymm11		,0x260(%%rbx)\n\t"\
			"vmovaps	%%ymm0		,    (%%rcx)\n\t			vmovaps	%%ymm8 		,0x240(%%rcx)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rcx)\n\t			vmovaps	%%ymm9 		,0x260(%%rcx)\n\t"\
			"/*****************************************/\n\t"\
			"/* now do three more radix-3 transforms: */\n\t"\
			"/*****************************************/\n\t"\
			"subq	$0x180		,%%rax 		/* __r00 */	\n\t"\
			"subq	$0x100,%%rbx	/* rbx stores r02+0x180 = r00+0x1c0; need r06 = r00+0xc0  = rbx-0x100 */\n\t"\
			"subq	$0x80,%%rcx		/* rcx stores r04+0x180            ; need r0c = r00+0x180 = rcx- 0x80 */\n\t"\
			"vmovaps	    (%%rbx)	,%%ymm4\n\t					vmovaps	0x240(%%rbx)	,%%ymm12\n\t"\
			"vmovaps	0x20(%%rbx)	,%%ymm5\n\t					vmovaps	0x260(%%rbx)	,%%ymm13\n\t"\
			"vmovaps	    (%%rcx)	,%%ymm2\n\t					vmovaps	0x240(%%rcx)	,%%ymm10\n\t"\
			"vmovaps	0x20(%%rcx)	,%%ymm3\n\t					vmovaps	0x260(%%rcx)	,%%ymm11\n\t"\
			"vmovaps	    (%%rax)	,%%ymm0\n\t					vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
			"vmovaps	0x20(%%rax)	,%%ymm1\n\t					vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
			"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t					vsubpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
			"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t					vsubpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm10,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm11,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t					vaddpd	%%ymm10,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm11,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm0		,    (%%rsi)\n\t			vmovaps	%%ymm8 		,    (%%r10)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rsi)\n\t			vmovaps	%%ymm9 		,0x20(%%r10)\n\t"\
			"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t					vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
			"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t					vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
			"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t					vmulpd	%%ymm7 ,%%ymm12,%%ymm12\n\t"\
			"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t					vmulpd	%%ymm7 ,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm8 ,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm9 ,%%ymm11,%%ymm11\n\t"\
			"vmovaps	%%ymm2		,%%ymm0\n\t					vmovaps	%%ymm10		,%%ymm8 \n\t"\
			"vmovaps	%%ymm3		,%%ymm1\n\t					vmovaps	%%ymm11		,%%ymm9 \n\t"\
			"addq	$0x600		,%%rdi 	/* __o3 = s1p24r */\n\t	subq	$0x300		,%%r11 		/* __o3 = s1p15r */\n\t"\
			"addq	$0x300		,%%rsi 	/* __o6 = s1p12r */\n\t	subq	$0x600		,%%r10 		/* __o6 = s1p03r */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm13,%%ymm10,%%ymm10\n\t"\
			"vsubpd	%%ymm4,%%ymm3,%%ymm3\n\t					vsubpd	%%ymm12,%%ymm11,%%ymm11\n\t"\
			"vsubpd	%%ymm5,%%ymm0,%%ymm0\n\t					vsubpd	%%ymm13,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm4,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm12,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm2		,    (%%rdi)\n\t			vmovaps	%%ymm10		,    (%%r11)\n\t"\
			"vmovaps	%%ymm3		,0x20(%%rdi)\n\t			vmovaps	%%ymm11		,0x20(%%r11)\n\t"\
			"vmovaps	%%ymm0		,    (%%rsi)\n\t			vmovaps	%%ymm8 		,    (%%r10)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rsi)\n\t			vmovaps	%%ymm9 		,0x20(%%r10)\n\t"\
			"addq	$0x40		,%%rax\n\t"\
			"addq	$0x40		,%%rbx\n\t"\
			"addq	$0x40		,%%rcx\n\t"\
			"subq	$0x80		,%%rdx 		/* c1 */\n\t"\
			"vmovaps	    (%%rbx)	,%%ymm2\n\t					vmovaps	0x240(%%rbx)	,%%ymm10\n\t"\
			"vmovaps	0x20(%%rbx)	,%%ymm3\n\t					vmovaps	0x260(%%rbx)	,%%ymm11\n\t"\
			"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
			"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
			"vmovaps	%%ymm2		,%%ymm0\n\t					vmovaps	%%ymm10		,%%ymm8 \n\t"\
			"vmovaps	%%ymm3		,%%ymm1\n\t					vmovaps	%%ymm11		,%%ymm9 \n\t"\
			"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t					vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
			"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t					vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
			"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t					vmulpd	%%ymm7 ,%%ymm8 ,%%ymm8 \n\t"\
			"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t					vmulpd	%%ymm7 ,%%ymm9 ,%%ymm9 \n\t"\
			"vaddpd	%%ymm1,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm9 ,%%ymm10,%%ymm10\n\t"\
			"vsubpd	%%ymm0,%%ymm3,%%ymm3\n\t					vsubpd	%%ymm8 ,%%ymm11,%%ymm11\n\t"\
			"addq	$0x40		,%%rdx 		/* c2 */\n\t"\
			"vmovaps	    (%%rcx)	,%%ymm4\n\t					vmovaps	0x240(%%rcx)	,%%ymm12\n\t"\
			"vmovaps	0x20(%%rcx)	,%%ymm5\n\t					vmovaps	0x260(%%rcx)	,%%ymm13\n\t"\
			"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
			"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
			"vmovaps	%%ymm4		,%%ymm0\n\t					vmovaps	%%ymm12		,%%ymm8 \n\t"\
			"vmovaps	%%ymm5		,%%ymm1\n\t					vmovaps	%%ymm13		,%%ymm9 \n\t"\
			"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t					vmulpd	%%ymm6 ,%%ymm12,%%ymm12\n\t"\
			"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t					vmulpd	%%ymm6 ,%%ymm13,%%ymm13\n\t"\
			"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t					vmulpd	%%ymm7 ,%%ymm8 ,%%ymm8 \n\t"\
			"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t					vmulpd	%%ymm7 ,%%ymm9 ,%%ymm9 \n\t"\
			"vaddpd	%%ymm1,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm9 ,%%ymm12,%%ymm12\n\t"\
			"vsubpd	%%ymm0,%%ymm5,%%ymm5\n\t					vsubpd	%%ymm8 ,%%ymm13,%%ymm13\n\t"\
			"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
			"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
			"vmovaps	    (%%rax)	,%%ymm0\n\t					vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
			"vmovaps	0x20(%%rax)	,%%ymm1\n\t					vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
			"subq	$0x100		,%%rsi 	/* __o7 = s1p08r */\n\t	addq	$0x800		,%%r10 		/* __o7 = s1p35r */\n\t"\
			"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t					vsubpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
			"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t					vsubpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm12,%%ymm12,%%ymm12\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t					vaddpd	%%ymm13,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
			"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t					vaddpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm4,%%ymm0,%%ymm0\n\t					vaddpd	%%ymm12,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm5,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm13,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm0		,    (%%rsi)\n\t			vmovaps	%%ymm8 		,    (%%r10)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rsi)\n\t			vmovaps	%%ymm9 		,0x20(%%r10)\n\t"\
			"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t					vmulpd	%%ymm6 ,%%ymm12,%%ymm12\n\t"\
			"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t					vmulpd	%%ymm6 ,%%ymm13,%%ymm13\n\t"\
			"vmulpd	%%ymm7,%%ymm2,%%ymm2\n\t					vmulpd	%%ymm7 ,%%ymm10,%%ymm10\n\t"\
			"vmulpd	%%ymm7,%%ymm3,%%ymm3\n\t					vmulpd	%%ymm7 ,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm8 ,%%ymm12,%%ymm12\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5\n\t					vaddpd	%%ymm9 ,%%ymm13,%%ymm13\n\t"\
			"addq	$0x200		,%%rdi 	/* __o1 = s1p32r */\n\t	addq	$0x200		,%%r11 		/* __o1 = s1p23r */\n\t"\
			"addq	$0x300		,%%rsi 	/* __o4 = s1p20r */\n\t	subq	$0x600		,%%r10 		/* __o4 = s1p11r */\n\t"\
			"vmovaps	%%ymm4		,%%ymm0\n\t					vmovaps	%%ymm12		,%%ymm8 \n\t"\
			"vmovaps	%%ymm5		,%%ymm1\n\t					vmovaps	%%ymm13		,%%ymm9 \n\t"\
			"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm11,%%ymm12,%%ymm12\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t					vsubpd	%%ymm10,%%ymm13,%%ymm13\n\t"\
			"vsubpd	%%ymm3,%%ymm0,%%ymm0\n\t					vsubpd	%%ymm11,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm2,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm10,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm4		,    (%%rdi)\n\t			vmovaps	%%ymm12		,    (%%r11)\n\t"\
			"vmovaps	%%ymm5		,0x20(%%rdi)\n\t			vmovaps	%%ymm13		,0x20(%%r11)\n\t"\
			"vmovaps	%%ymm0		,    (%%rsi)\n\t			vmovaps	%%ymm8 		,    (%%r10)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rsi)\n\t			vmovaps	%%ymm9 		,0x20(%%r10)\n\t"\
			"addq	$0x40		,%%rax\n\t"\
			"addq	$0x40		,%%rbx\n\t"\
			"addq	$0x40		,%%rcx\n\t"\
			"subq	$0x40		,%%rdx 		/* c2 */\n\t"\
			"vmovaps	    (%%rbx)	,%%ymm2\n\t					vmovaps	0x240(%%rbx)	,%%ymm10\n\t"\
			"vmovaps	0x20(%%rbx)	,%%ymm3\n\t					vmovaps	0x260(%%rbx)	,%%ymm11\n\t"\
			"vmovaps	    (%%rdx)	,%%ymm6\n\t	"\
			"vmovaps	0x20(%%rdx)	,%%ymm7\n\t	"\
			"vmovaps	%%ymm2		,%%ymm0\n\t					vmovaps	%%ymm10		,%%ymm8 \n\t"\
			"vmovaps	%%ymm3		,%%ymm1\n\t					vmovaps	%%ymm11		,%%ymm9 \n\t"\
			"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t					vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
			"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t					vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
			"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t					vmulpd	%%ymm7 ,%%ymm8 ,%%ymm8 \n\t"\
			"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t					vmulpd	%%ymm7 ,%%ymm9 ,%%ymm9 \n\t"\
			"vaddpd	%%ymm1,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm9 ,%%ymm10,%%ymm10\n\t"\
			"vsubpd	%%ymm0,%%ymm3,%%ymm3\n\t					vsubpd	%%ymm8 ,%%ymm11,%%ymm11\n\t"\
			"addq	$0x80		,%%rdx 		/* c4 */\n\t"\
			"vmovaps	    (%%rcx)	,%%ymm4\n\t					vmovaps	0x240(%%rcx)	,%%ymm12\n\t"\
			"vmovaps	0x20(%%rcx)	,%%ymm5\n\t					vmovaps	0x260(%%rcx)	,%%ymm13\n\t"\
			"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
			"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
			"subq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"vmovaps	%%ymm4		,%%ymm0\n\t					vmovaps	%%ymm12		,%%ymm8 \n\t"\
			"vmovaps	%%ymm5		,%%ymm1\n\t					vmovaps	%%ymm13		,%%ymm9 \n\t"\
			"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t					vmulpd	%%ymm6 ,%%ymm12,%%ymm12\n\t"\
			"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t					vmulpd	%%ymm6 ,%%ymm13,%%ymm13\n\t"\
			"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t					vmulpd	%%ymm7 ,%%ymm8 ,%%ymm8 \n\t"\
			"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t					vmulpd	%%ymm7 ,%%ymm9 ,%%ymm9 \n\t"\
			"vaddpd	%%ymm1,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm9 ,%%ymm12,%%ymm12\n\t"\
			"vsubpd	%%ymm0,%%ymm5,%%ymm5\n\t					vsubpd	%%ymm8 ,%%ymm13,%%ymm13\n\t"\
			"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
			"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
			"vmovaps	    (%%rax)	,%%ymm0\n\t					vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
			"vmovaps	0x20(%%rax)	,%%ymm1\n\t					vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
			"subq	$0x100		,%%rsi 	/* __o5 = s1p16r */\n\t	subq	$0x100		,%%r10 		/* __o5 = s1p07r */\n\t"\
			"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t					vsubpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
			"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t					vsubpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm12,%%ymm12,%%ymm12\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t					vaddpd	%%ymm13,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
			"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t					vaddpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm4,%%ymm0,%%ymm0\n\t					vaddpd	%%ymm12,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm5,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm13,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm0		,    (%%rsi)\n\t			vmovaps	%%ymm8 		,    (%%r10)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rsi)\n\t			vmovaps	%%ymm9 		,0x20(%%r10)\n\t"\
			"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t					vmulpd	%%ymm6 ,%%ymm12,%%ymm12\n\t"\
			"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t					vmulpd	%%ymm6 ,%%ymm13,%%ymm13\n\t"\
			"vmulpd	%%ymm7,%%ymm2,%%ymm2\n\t					vmulpd	%%ymm7 ,%%ymm10,%%ymm10\n\t"\
			"vmulpd	%%ymm7,%%ymm3,%%ymm3\n\t					vmulpd	%%ymm7 ,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm8 ,%%ymm12,%%ymm12\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5\n\t					vaddpd	%%ymm9 ,%%ymm13,%%ymm13\n\t"\
			"subq	$0x700		,%%rdi 	/* __o8 = s1p04r */\n\t	addq	$0x200		,%%r11 		/* __o8 = s1p31r */\n\t"\
			"addq	$0x300		,%%rsi 	/* __o2 = s1p28r */\n\t	addq	$0x300		,%%r10 		/* __o2 = s1p19r */\n\t"\
			"vmovaps	%%ymm4		,%%ymm0\n\t					vmovaps	%%ymm12		,%%ymm8 \n\t"\
			"vmovaps	%%ymm5		,%%ymm1\n\t					vmovaps	%%ymm13		,%%ymm9 \n\t"\
			"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm11,%%ymm12,%%ymm12\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t					vsubpd	%%ymm10,%%ymm13,%%ymm13\n\t"\
			"vsubpd	%%ymm3,%%ymm0,%%ymm0\n\t					vsubpd	%%ymm11,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm2,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm10,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm4		,    (%%rdi)\n\t			vmovaps	%%ymm12		,    (%%r11)\n\t"\
			"vmovaps	%%ymm5		,0x20(%%rdi)\n\t			vmovaps	%%ymm13		,0x20(%%r11)\n\t"\
			"vmovaps	%%ymm0		,    (%%rsi)\n\t			vmovaps	%%ymm8 		,    (%%r10)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rsi)\n\t			vmovaps	%%ymm9 		,0x20(%%r10)\n\t"\
			/* SSE2_RADIX_09_DIT_0TWIDDLE(r20,s1p[18,14,10,06,02,34,30,26,22]r): */\
			"movq	%[__out]	,%%rsi\n\t						movq	%[__r00]	,%%rax\n\t"/* rcol locals based on r10 = r00+0x240 */\
			"addq	$0x480		,%%rsi 	/* s1p18 */\n\t			addq	$0x480		,%%rax\n\t	/* __r20 */\n\t"/* rcol locals based on r30 = r20+0x120 */\
			"movq	%%rsi		,%%rdi\n\t						movq	%[__cc1]	,%%rdx	\n\t"\
			"movq	%%rax		,%%rbx\n\t"\
			"movq	%%rax		,%%rcx\n\t"						/* SSE2_RADIX_09_DIT_0TWIDDLE(r30,s1p[09,05,01,33,29,25,21,17,13]r): */\
			"addq	$0x40		,%%rbx\n\t						/* rcol: r10,r11 store o-addresses */\n\t"\
			"addq	$0x80		,%%rcx\n\t						leaq	-0x240(%%rsi),%%r10 		/* s1p09 */\n\t"\
			"addq	$0x80		,%%rdx 		/* c3m1 */\n\t		movq	%%r10		,%%r11\n\t"\
			"vmovaps	    (%%rbx)	,%%ymm2\n\t					vmovaps	0x240(%%rbx)	,%%ymm10\n\t"\
			"vmovaps	0x20(%%rbx)	,%%ymm3\n\t					vmovaps	0x260(%%rbx)	,%%ymm11\n\t"\
			"vmovaps	    (%%rax)	,%%ymm0\n\t					vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
			"vmovaps	0x20(%%rax)	,%%ymm1\n\t					vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
			"vmovaps	    (%%rcx)	,%%ymm6\n\t					vmovaps	0x240(%%rcx)	,%%ymm14\n\t"\
			"vmovaps	0x20(%%rcx)	,%%ymm7\n\t					vmovaps	0x260(%%rcx)	,%%ymm15\n\t"\
			"vmovaps	%%ymm2		,%%ymm4\n\t					vmovaps	%%ymm10		,%%ymm12\n\t"\
			"vmovaps	%%ymm3		,%%ymm5\n\t					vmovaps	%%ymm11		,%%ymm13\n\t"\
			"vaddpd	%%ymm6,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm14,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm7,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm15,%%ymm11,%%ymm11\n\t"\
			"vsubpd	%%ymm6,%%ymm4,%%ymm4\n\t					vsubpd	%%ymm14,%%ymm12,%%ymm12\n\t"\
			"vsubpd	%%ymm7,%%ymm5,%%ymm5\n\t					vsubpd	%%ymm15,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t					vaddpd	%%ymm10,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm11,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
			"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
			"vmovaps	%%ymm0		,    (%%rax)\n\t			vmovaps	%%ymm8 		,0x240(%%rax)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rax)\n\t			vmovaps	%%ymm9 		,0x260(%%rax)\n\t"\
			"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t					vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
			"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t					vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
			"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t					vmulpd	%%ymm7 ,%%ymm12,%%ymm12\n\t"\
			"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t					vmulpd	%%ymm7 ,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm8 ,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm9 ,%%ymm11,%%ymm11\n\t"\
			"vmovaps	%%ymm2		,%%ymm0\n\t					vmovaps	%%ymm10		,%%ymm8 \n\t"\
			"vmovaps	%%ymm3		,%%ymm1\n\t					vmovaps	%%ymm11		,%%ymm9 \n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm13,%%ymm10,%%ymm10\n\t"\
			"vsubpd	%%ymm4,%%ymm3,%%ymm3\n\t					vsubpd	%%ymm12,%%ymm11,%%ymm11\n\t"\
			"vsubpd	%%ymm5,%%ymm0,%%ymm0\n\t					vsubpd	%%ymm13,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm4,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm12,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm2		,    (%%rbx)\n\t			vmovaps	%%ymm10		,0x240(%%rbx)\n\t"\
			"vmovaps	%%ymm3		,0x20(%%rbx)\n\t			vmovaps	%%ymm11		,0x260(%%rbx)\n\t"\
			"vmovaps	%%ymm0		,    (%%rcx)\n\t			vmovaps	%%ymm8 		,0x240(%%rcx)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rcx)\n\t			vmovaps	%%ymm9 		,0x260(%%rcx)\n\t"\
			"addq	$0xc0		,%%rax\n\t"\
			"addq	$0xc0		,%%rbx\n\t"\
			"addq	$0xc0		,%%rcx\n\t"\
			"vmovaps	    (%%rbx)	,%%ymm4\n\t					vmovaps	0x240(%%rbx)	,%%ymm12\n\t"\
			"vmovaps	0x20(%%rbx)	,%%ymm5\n\t					vmovaps	0x260(%%rbx)	,%%ymm13\n\t"\
			"vmovaps	    (%%rcx)	,%%ymm2\n\t					vmovaps	0x240(%%rcx)	,%%ymm10\n\t"\
			"vmovaps	0x20(%%rcx)	,%%ymm3\n\t					vmovaps	0x260(%%rcx)	,%%ymm11\n\t"\
			"vmovaps	    (%%rax)	,%%ymm0\n\t					vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
			"vmovaps	0x20(%%rax)	,%%ymm1\n\t					vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
			"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t					vsubpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
			"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t					vsubpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm10,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm11,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t					vaddpd	%%ymm10,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm11,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm0		,    (%%rax)\n\t			vmovaps	%%ymm8 		,0x240(%%rax)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rax)\n\t			vmovaps	%%ymm9 		,0x260(%%rax)\n\t"\
			"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t					vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
			"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t					vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
			"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t					vmulpd	%%ymm7 ,%%ymm12,%%ymm12\n\t"\
			"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t					vmulpd	%%ymm7 ,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm8 ,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm9 ,%%ymm11,%%ymm11\n\t"\
			"vmovaps	%%ymm2		,%%ymm0\n\t					vmovaps	%%ymm10		,%%ymm8 \n\t"\
			"vmovaps	%%ymm3		,%%ymm1\n\t					vmovaps	%%ymm11		,%%ymm9 \n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm13,%%ymm10,%%ymm10\n\t"\
			"vsubpd	%%ymm4,%%ymm3,%%ymm3\n\t					vsubpd	%%ymm12,%%ymm11,%%ymm11\n\t"\
			"vsubpd	%%ymm5,%%ymm0,%%ymm0\n\t					vsubpd	%%ymm13,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm4,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm12,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm2		,    (%%rbx)\n\t			vmovaps	%%ymm10		,0x240(%%rbx)\n\t"\
			"vmovaps	%%ymm3		,0x20(%%rbx)\n\t			vmovaps	%%ymm11		,0x260(%%rbx)\n\t"\
			"vmovaps	%%ymm0		,    (%%rcx)\n\t			vmovaps	%%ymm8 		,0x240(%%rcx)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rcx)\n\t			vmovaps	%%ymm9 		,0x260(%%rcx)\n\t"\
			"addq	$0xc0		,%%rax\n\t"\
			"addq	$0xc0		,%%rbx\n\t"\
			"addq	$0xc0		,%%rcx\n\t"\
			"vmovaps	    (%%rbx)	,%%ymm4\n\t					vmovaps	0x240(%%rbx)	,%%ymm12\n\t"\
			"vmovaps	0x20(%%rbx)	,%%ymm5\n\t					vmovaps	0x260(%%rbx)	,%%ymm13\n\t"\
			"vmovaps	    (%%rcx)	,%%ymm2\n\t					vmovaps	0x240(%%rcx)	,%%ymm10\n\t"\
			"vmovaps	0x20(%%rcx)	,%%ymm3\n\t					vmovaps	0x260(%%rcx)	,%%ymm11\n\t"\
			"vmovaps	    (%%rax)	,%%ymm0\n\t					vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
			"vmovaps	0x20(%%rax)	,%%ymm1\n\t					vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
			"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t					vsubpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
			"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t					vsubpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm10,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm11,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t					vaddpd	%%ymm10,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm11,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm0		,    (%%rax)\n\t			vmovaps	%%ymm8 		,0x240(%%rax)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rax)\n\t			vmovaps	%%ymm9 		,0x260(%%rax)\n\t"\
			"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t					vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
			"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t					vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
			"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t					vmulpd	%%ymm7 ,%%ymm12,%%ymm12\n\t"\
			"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t					vmulpd	%%ymm7 ,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm8 ,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm9 ,%%ymm11,%%ymm11\n\t"\
			"vmovaps	%%ymm2		,%%ymm0\n\t					vmovaps	%%ymm10		,%%ymm8 \n\t"\
			"vmovaps	%%ymm3		,%%ymm1\n\t					vmovaps	%%ymm11		,%%ymm9 \n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm13,%%ymm10,%%ymm10\n\t"\
			"vsubpd	%%ymm4,%%ymm3,%%ymm3\n\t					vsubpd	%%ymm12,%%ymm11,%%ymm11\n\t"\
			"vsubpd	%%ymm5,%%ymm0,%%ymm0\n\t					vsubpd	%%ymm13,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm4,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm12,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm2		,    (%%rbx)\n\t			vmovaps	%%ymm10		,0x240(%%rbx)\n\t"\
			"vmovaps	%%ymm3		,0x20(%%rbx)\n\t			vmovaps	%%ymm11		,0x260(%%rbx)\n\t"\
			"vmovaps	%%ymm0		,    (%%rcx)\n\t			vmovaps	%%ymm8 		,0x240(%%rcx)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rcx)\n\t			vmovaps	%%ymm9 		,0x260(%%rcx)\n\t"\
			"/*****************************************/\n\t"\
			"/* now do three more radix-3 transforms: */\n\t"\
			"/*****************************************/\n\t"\
			"subq	$0x180		,%%rax\n\t"\
			"subq	$0x100		,%%rbx\n\t"\
			"subq	$0x80		,%%rcx\n\t"\
			"vmovaps	    (%%rbx)	,%%ymm4\n\t					vmovaps	0x240(%%rbx)	,%%ymm12\n\t"\
			"vmovaps	0x20(%%rbx)	,%%ymm5\n\t					vmovaps	0x260(%%rbx)	,%%ymm13\n\t"\
			"vmovaps	    (%%rcx)	,%%ymm2\n\t					vmovaps	0x240(%%rcx)	,%%ymm10\n\t"\
			"vmovaps	0x20(%%rcx)	,%%ymm3\n\t					vmovaps	0x260(%%rcx)	,%%ymm11\n\t"\
			"vmovaps	    (%%rax)	,%%ymm0\n\t					vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
			"vmovaps	0x20(%%rax)	,%%ymm1\n\t					vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
			"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t					vsubpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
			"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t					vsubpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm10,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm11,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t					vaddpd	%%ymm10,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm11,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm0		,    (%%rsi)\n\t			vmovaps	%%ymm8 		,    (%%r10)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rsi)\n\t			vmovaps	%%ymm9 		,0x20(%%r10)\n\t"\
			"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t					vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
			"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t					vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
			"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t					vmulpd	%%ymm7 ,%%ymm12,%%ymm12\n\t"\
			"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t					vmulpd	%%ymm7 ,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm8 ,%%ymm10,%%ymm10\n\t"\
			"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t					vaddpd	%%ymm9 ,%%ymm11,%%ymm11\n\t"\
			"vmovaps	%%ymm2		,%%ymm0\n\t					vmovaps	%%ymm10		,%%ymm8 \n\t"\
			"vmovaps	%%ymm3		,%%ymm1\n\t					vmovaps	%%ymm11		,%%ymm9 \n\t"\
			"subq	$0x300		,%%rdi 	/* __o3 = s1p06r */\n\t	addq	$0x600		,%%r11 		/* __o3 = s1p33r */\n\t"\
			"addq	$0x300		,%%rsi 	/* __o6 = s1p30r */\n\t	addq	$0x300		,%%r10 		/* __o6 = s1p21r */\n\t"\
			"vaddpd	%%ymm5,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm13,%%ymm10,%%ymm10\n\t"\
			"vsubpd	%%ymm4,%%ymm3,%%ymm3\n\t					vsubpd	%%ymm12,%%ymm11,%%ymm11\n\t"\
			"vsubpd	%%ymm5,%%ymm0,%%ymm0\n\t					vsubpd	%%ymm13,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm4,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm12,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm2		,    (%%rdi)\n\t			vmovaps	%%ymm10		,    (%%r11)\n\t"\
			"vmovaps	%%ymm3		,0x20(%%rdi)\n\t			vmovaps	%%ymm11		,0x20(%%r11)\n\t"\
			"vmovaps	%%ymm0		,    (%%rsi)\n\t			vmovaps	%%ymm8 		,    (%%r10)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rsi)\n\t			vmovaps	%%ymm9 		,0x20(%%r10)\n\t"\
			"addq	$0x40		,%%rax\n\t"\
			"addq	$0x40		,%%rbx\n\t"\
			"addq	$0x40		,%%rcx\n\t"\
			"subq	$0x80		,%%rdx 		/* c1 */\n\t"\
			"vmovaps	    (%%rbx)	,%%ymm2\n\t					vmovaps	0x240(%%rbx)	,%%ymm10\n\t"\
			"vmovaps	0x20(%%rbx)	,%%ymm3\n\t					vmovaps	0x260(%%rbx)	,%%ymm11\n\t"\
			"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
			"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
			"vmovaps	%%ymm2		,%%ymm0\n\t					vmovaps	%%ymm10		,%%ymm8 \n\t"\
			"vmovaps	%%ymm3		,%%ymm1\n\t					vmovaps	%%ymm11		,%%ymm9 \n\t"\
			"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t					vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
			"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t					vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
			"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t					vmulpd	%%ymm7 ,%%ymm8 ,%%ymm8 \n\t"\
			"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t					vmulpd	%%ymm7 ,%%ymm9 ,%%ymm9 \n\t"\
			"vaddpd	%%ymm1,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm9 ,%%ymm10,%%ymm10\n\t"\
			"vsubpd	%%ymm0,%%ymm3,%%ymm3\n\t					vsubpd	%%ymm8 ,%%ymm11,%%ymm11\n\t"\
			"addq	$0x40		,%%rdx 		/* c2 */\n\t"\
			"vmovaps	    (%%rcx)	,%%ymm4\n\t					vmovaps	0x240(%%rcx)	,%%ymm12\n\t"\
			"vmovaps	0x20(%%rcx)	,%%ymm5\n\t					vmovaps	0x260(%%rcx)	,%%ymm13\n\t"\
			"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
			"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
			"vmovaps	%%ymm4		,%%ymm0\n\t					vmovaps	%%ymm12		,%%ymm8 \n\t"\
			"vmovaps	%%ymm5		,%%ymm1\n\t					vmovaps	%%ymm13		,%%ymm9 \n\t"\
			"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t					vmulpd	%%ymm6 ,%%ymm12,%%ymm12\n\t"\
			"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t					vmulpd	%%ymm6 ,%%ymm13,%%ymm13\n\t"\
			"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t					vmulpd	%%ymm7 ,%%ymm8 ,%%ymm8 \n\t"\
			"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t					vmulpd	%%ymm7 ,%%ymm9 ,%%ymm9 \n\t"\
			"vaddpd	%%ymm1,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm9 ,%%ymm12,%%ymm12\n\t"\
			"vsubpd	%%ymm0,%%ymm5,%%ymm5\n\t					vsubpd	%%ymm8 ,%%ymm13,%%ymm13\n\t"\
			"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
			"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
			"vmovaps	    (%%rax)	,%%ymm0\n\t					vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
			"vmovaps	0x20(%%rax)	,%%ymm1\n\t					vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
			"subq	$0x100		,%%rsi 	/* __o7 = s1p26r */\n\t	subq	$0x100		,%%r10 		/* __o7 = s1p17r */\n\t"\
			"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t					vsubpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
			"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t					vsubpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm12,%%ymm12,%%ymm12\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t					vaddpd	%%ymm13,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
			"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t					vaddpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm4,%%ymm0,%%ymm0\n\t					vaddpd	%%ymm12,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm5,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm13,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm0		,    (%%rsi)\n\t			vmovaps	%%ymm8 		,    (%%r10)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rsi)\n\t			vmovaps	%%ymm9 		,0x20(%%r10)\n\t"\
			"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t					vmulpd	%%ymm6 ,%%ymm12,%%ymm12\n\t"\
			"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t					vmulpd	%%ymm6 ,%%ymm13,%%ymm13\n\t"\
			"vmulpd	%%ymm7,%%ymm2,%%ymm2\n\t					vmulpd	%%ymm7 ,%%ymm10,%%ymm10\n\t"\
			"vmulpd	%%ymm7,%%ymm3,%%ymm3\n\t					vmulpd	%%ymm7 ,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm8 ,%%ymm12,%%ymm12\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5\n\t					vaddpd	%%ymm9 ,%%ymm13,%%ymm13\n\t"\
			"addq	$0x200		,%%rdi 	/* __o1 = s1p14r */\n\t	subq	$0x700		,%%r11 		/* __o1 = s1p05r */\n\t"\
			"subq	$0x600		,%%rsi 	/* __o4 = s1p02r */\n\t	addq	$0x300		,%%r10 		/* __o4 = s1p29r */\n\t"\
			"vmovaps	%%ymm4		,%%ymm0\n\t					vmovaps	%%ymm12		,%%ymm8 \n\t"\
			"vmovaps	%%ymm5		,%%ymm1\n\t					vmovaps	%%ymm13		,%%ymm9 \n\t"\
			"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm11,%%ymm12,%%ymm12\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t					vsubpd	%%ymm10,%%ymm13,%%ymm13\n\t"\
			"vsubpd	%%ymm3,%%ymm0,%%ymm0\n\t					vsubpd	%%ymm11,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm2,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm10,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm4		,    (%%rdi)\n\t			vmovaps	%%ymm12		,    (%%r11)\n\t"\
			"vmovaps	%%ymm5		,0x20(%%rdi)\n\t			vmovaps	%%ymm13		,0x20(%%r11)\n\t"\
			"vmovaps	%%ymm0		,    (%%rsi)\n\t			vmovaps	%%ymm8 		,    (%%r10)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rsi)\n\t			vmovaps	%%ymm9 		,0x20(%%r10)\n\t"\
			"addq	$0x40		,%%rax\n\t"\
			"addq	$0x40		,%%rbx\n\t"\
			"addq	$0x40		,%%rcx\n\t"\
			"subq	$0x40		,%%rdx 		/* c2 */\n\t"\
			"vmovaps	    (%%rbx)	,%%ymm2\n\t					vmovaps	0x240(%%rbx)	,%%ymm10\n\t"\
			"vmovaps	0x20(%%rbx)	,%%ymm3\n\t					vmovaps	0x260(%%rbx)	,%%ymm11\n\t"\
			"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
			"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
			"vmovaps	%%ymm2		,%%ymm0\n\t					vmovaps	%%ymm10		,%%ymm8 \n\t"\
			"vmovaps	%%ymm3		,%%ymm1\n\t					vmovaps	%%ymm11		,%%ymm9 \n\t"\
			"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t					vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
			"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t					vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
			"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t					vmulpd	%%ymm7 ,%%ymm8 ,%%ymm8 \n\t"\
			"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t					vmulpd	%%ymm7 ,%%ymm9 ,%%ymm9 \n\t"\
			"vaddpd	%%ymm1,%%ymm2,%%ymm2\n\t					vaddpd	%%ymm9 ,%%ymm10,%%ymm10\n\t"\
			"vsubpd	%%ymm0,%%ymm3,%%ymm3\n\t					vsubpd	%%ymm8 ,%%ymm11,%%ymm11\n\t"\
			"addq	$0x80		,%%rdx 		/* c4 */\n\t"\
			"vmovaps	    (%%rcx)	,%%ymm4\n\t					vmovaps	0x240(%%rcx)	,%%ymm12\n\t"\
			"vmovaps	0x20(%%rcx)	,%%ymm5\n\t					vmovaps	0x260(%%rcx)	,%%ymm13\n\t"\
			"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
			"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
			"subq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"vmovaps	%%ymm4		,%%ymm0\n\t					vmovaps	%%ymm12		,%%ymm8 \n\t"\
			"vmovaps	%%ymm5		,%%ymm1\n\t					vmovaps	%%ymm13		,%%ymm9 \n\t"\
			"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t					vmulpd	%%ymm6 ,%%ymm12,%%ymm12\n\t"\
			"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t					vmulpd	%%ymm6 ,%%ymm13,%%ymm13\n\t"\
			"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t					vmulpd	%%ymm7 ,%%ymm8 ,%%ymm8 \n\t"\
			"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t					vmulpd	%%ymm7 ,%%ymm9 ,%%ymm9 \n\t"\
			"vaddpd	%%ymm1,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm9 ,%%ymm12,%%ymm12\n\t"\
			"vsubpd	%%ymm0,%%ymm5,%%ymm5\n\t					vsubpd	%%ymm8 ,%%ymm13,%%ymm13\n\t"\
			"vmovaps	    (%%rdx)	,%%ymm6\n\t	"\
			"vmovaps	0x20(%%rdx)	,%%ymm7\n\t	"\
			"vmovaps	    (%%rax)	,%%ymm0\n\t					vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
			"vmovaps	0x20(%%rax)	,%%ymm1\n\t					vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
			"addq	$0x800		,%%rsi 	/* __o5 = s1p34r */\n\t	subq	$0x100		,%%r10 		/* __o5 = s1p25r */\n\t"\
			"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t					vsubpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
			"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t					vsubpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm12,%%ymm12,%%ymm12\n\t"\
			"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t					vaddpd	%%ymm13,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
			"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t					vaddpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
			"vaddpd	%%ymm4,%%ymm0,%%ymm0\n\t					vaddpd	%%ymm12,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm5,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm13,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm0		,    (%%rsi)\n\t			vmovaps	%%ymm8 		,    (%%r10)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rsi)\n\t			vmovaps	%%ymm9 		,0x20(%%r10)\n\t"\
			"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t					vmulpd	%%ymm6 ,%%ymm12,%%ymm12\n\t"\
			"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t					vmulpd	%%ymm6 ,%%ymm13,%%ymm13\n\t"\
			"vmulpd	%%ymm7,%%ymm2,%%ymm2\n\t					vmulpd	%%ymm7 ,%%ymm10,%%ymm10\n\t"\
			"vmulpd	%%ymm7,%%ymm3,%%ymm3\n\t					vmulpd	%%ymm7 ,%%ymm11,%%ymm11\n\t"\
			"vaddpd	%%ymm0,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm8 ,%%ymm12,%%ymm12\n\t"\
			"vaddpd	%%ymm1,%%ymm5,%%ymm5\n\t					vaddpd	%%ymm9 ,%%ymm13,%%ymm13\n\t"\
			"addq	$0x200		,%%rdi 	/* __o8 = s1p22r */\n\t	addq	$0x200		,%%r11 		/* __o8 = s1p13r */\n\t"\
			"subq	$0x600		,%%rsi 	/* __o2 = s1p10r */\n\t	subq	$0x600		,%%r10 		/* __o2 = s1p01r */\n\t"\
			"vmovaps	%%ymm4		,%%ymm0\n\t					vmovaps	%%ymm12		,%%ymm8 \n\t"\
			"vmovaps	%%ymm5		,%%ymm1\n\t					vmovaps	%%ymm13		,%%ymm9 \n\t"\
			"vaddpd	%%ymm3,%%ymm4,%%ymm4\n\t					vaddpd	%%ymm11,%%ymm12,%%ymm12\n\t"\
			"vsubpd	%%ymm2,%%ymm5,%%ymm5\n\t					vsubpd	%%ymm10,%%ymm13,%%ymm13\n\t"\
			"vsubpd	%%ymm3,%%ymm0,%%ymm0\n\t					vsubpd	%%ymm11,%%ymm8 ,%%ymm8 \n\t"\
			"vaddpd	%%ymm2,%%ymm1,%%ymm1\n\t					vaddpd	%%ymm10,%%ymm9 ,%%ymm9 \n\t"\
			"vmovaps	%%ymm4		,    (%%rdi)\n\t			vmovaps	%%ymm12		,    (%%r11)\n\t"\
			"vmovaps	%%ymm5		,0x20(%%rdi)\n\t			vmovaps	%%ymm13		,0x20(%%r11)\n\t"\
			"vmovaps	%%ymm0		,    (%%rsi)\n\t			vmovaps	%%ymm8 		,    (%%r10)\n\t"\
			"vmovaps	%%ymm1		,0x20(%%rsi)\n\t			vmovaps	%%ymm9 		,0x20(%%r10)\n\t"\
			"\n\t"\
			:					/* outputs: none */\
			: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
			 ,[__p01] "m" (Xp01)\
			 ,[__p02] "m" (Xp02)\
			 ,[__p03] "m" (Xp03)\
			 ,[__p04] "m" (Xp04)\
			 ,[__p08] "m" (Xp08)\
			 ,[__p12] "m" (Xp12)\
			 ,[__p16] "m" (Xp16)\
			 ,[__p20] "m" (Xp20)\
			 ,[__p24] "m" (Xp24)\
			 ,[__p28] "m" (Xp28)\
			 ,[__p32] "m" (Xp32)\
			 ,[__r00] "m" (Xr00)\
			 ,[__out] "m" (Xout)\
			 ,[__cc1] "m" (Xcc1)\
			: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
	}


	#define	SSE2_RADIX36_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xp28,Xp32,Xr00,Xin0,Xcc1)\
	{\
	__asm__ volatile (\
		/* SSE2_RADIX_09_DIF_0TWIDDLE(r00,s1p[00,32,28,24,20,16,12,08,04]r: */\
		"movq	%[__in0]	,%%rax 		/* __i0-8; e[abc]x store input addresses */\n\t"\
		"movq	%[__cc1]	,%%rdx 		/* edx stores trig addresses throughout */\n\t"\
		"movq	%[__r00]	,%%rsi\n\t"								/* SSE2_RADIX_09_DIF_0TWIDDLE(r10,s1p[27,23,19,15,11,07,03,35,31]r: */\
		"movq	%%rax		,%%rbx\n\t								leaq	0x6c0(%%rax),%%r10 		/* __i0 = s1p27r */\n\t"\
		"movq	%%rax		,%%rcx\n\t								movq	%%r10		,%%r11\n\t"\
		"addq	$0x600		,%%rbx 	/* __i3 = s1p24r */\n\t			movq	%%r10		,%%r12\n\t"\
		"addq	$0x300		,%%rcx 	/* __i6 = s1p12r */\n\t			subq	$0x300		,%%r11 		/* __i3 = s1p15r */\n\t"\
		"addq	$0x80		,%%rdx 		/* c3m1 */\n\t				subq	$0x600		,%%r12 		/* __i6 = s1p03r */\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2\n\t							vmovaps	    (%%r11)	,%%ymm10\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3\n\t							vmovaps	0x20(%%r11)	,%%ymm11\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t							vmovaps	    (%%r10)	,%%ymm8 \n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t							vmovaps	0x20(%%r10)	,%%ymm9 \n\t"\
		"vmovaps	    (%%rcx)	,%%ymm6\n\t							vmovaps	    (%%r12)	,%%ymm14\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm7\n\t							vmovaps	0x20(%%r12)	,%%ymm15\n\t"\
		"vmovaps	%%ymm2		,%%ymm4\n\t							vmovaps	%%ymm10		,%%ymm12\n\t"\
		"vmovaps	%%ymm3		,%%ymm5\n\t							vmovaps	%%ymm11		,%%ymm13\n\t"\
		"movq	%%rsi		,%%rdi\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm14,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm15,%%ymm11,%%ymm11\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4\n\t							vsubpd	%%ymm14,%%ymm12,%%ymm12\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5\n\t							vsubpd	%%ymm15,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm10,%%ymm8 ,%%ymm8 \n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t							vaddpd	%%ymm11,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)	/* <- t00 */\n\t	vmovaps	%%ymm8 		,0x240(%%rsi)	/* <- t00 */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t01 */\n\t	vmovaps	%%ymm9 		,0x260(%%rsi)	/* <- t01 */\n\t"\
		"addq	$0x80		,%%rsi\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t							vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t							vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t							vmulpd	%%ymm7 ,%%ymm12,%%ymm12\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t							vmulpd	%%ymm7 ,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm8 ,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm9 ,%%ymm11,%%ymm11\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t							vmovaps	%%ymm10		,%%ymm8 \n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t							vmovaps	%%ymm11		,%%ymm9 \n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2\n\t							vsubpd	%%ymm13,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm12,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm13,%%ymm8 ,%%ymm8 \n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1\n\t							vsubpd	%%ymm12,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm2		,    (%%rdi)	/* <- t02 */\n\t	vmovaps	%%ymm10		,0x240(%%rdi)	/* <- t02 */\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rdi)	/* <- t03 */\n\t	vmovaps	%%ymm11		,0x260(%%rdi)	/* <- t03 */\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)	/* <- t04 */\n\t	vmovaps	%%ymm8 		,0x240(%%rsi)	/* <- t04 */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t05 */\n\t	vmovaps	%%ymm9 		,0x260(%%rsi)	/* <- t05 */\n\t"\
		"addq	$0x800		,%%rax 	/* __i1 = s1p32r */\n\t			subq	$0x100		,%%r10 		/* __i1 = s1p23r */\n\t"\
		"subq	$0x100		,%%rbx 	/* __i4 = s1p20r */\n\t			subq	$0x100		,%%r11 		/* __i4 = s1p11r */\n\t"\
		"subq	$0x100		,%%rcx 	/* __i7 = s1p08r */\n\t			addq	$0x800		,%%r12 		/* __i7 = s1p35r */\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4\n\t							vmovaps	    (%%r11)	,%%ymm12\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5\n\t							vmovaps	0x20(%%r11)	,%%ymm13\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t							vmovaps	    (%%r10)	,%%ymm8 \n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t							vmovaps	0x20(%%r10)	,%%ymm9 \n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2\n\t							vmovaps	    (%%r12)	,%%ymm10\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3\n\t							vmovaps	0x20(%%r12)	,%%ymm11\n\t"\
		"addq	$0x80		,%%rdi\n\t"\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t							vsubpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t							vsubpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm10,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm11,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm10,%%ymm8 ,%%ymm8 \n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t							vaddpd	%%ymm11,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm0		,    (%%rdi)	/* <- t06 */\n\t	vmovaps	%%ymm8 		,0x240(%%rdi)	/* <- t06 */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rdi)	/* <- t07 */\n\t	vmovaps	%%ymm9 		,0x260(%%rdi)	/* <- t07 */\n\t"\
		"addq	$0x80		,%%rsi\n\t"\
		"addq	$0x80		,%%rdi\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t							vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t							vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t							vmulpd	%%ymm7 ,%%ymm12,%%ymm12\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t							vmulpd	%%ymm7 ,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm8 ,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm9 ,%%ymm11,%%ymm11\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t							vmovaps	%%ymm10		,%%ymm8 \n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t							vmovaps	%%ymm11		,%%ymm9 \n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2\n\t							vsubpd	%%ymm13,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm12,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm13,%%ymm8 ,%%ymm8 \n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1\n\t							vsubpd	%%ymm12,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm2		,    (%%rsi)	/* <- t08 */\n\t	vmovaps	%%ymm10		,0x240(%%rsi)	/* <- t08 */\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rsi)	/* <- t09 */\n\t	vmovaps	%%ymm11		,0x260(%%rsi)	/* <- t09 */\n\t"\
		"vmovaps	%%ymm0		,    (%%rdi)	/* <- t0a */\n\t	vmovaps	%%ymm8 		,0x240(%%rdi)	/* <- t0a */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rdi)	/* <- t0b */\n\t	vmovaps	%%ymm9 		,0x260(%%rdi)	/* <- t0b */\n\t"\
		"addq	$0x80		,%%rsi\n\t"\
		"subq	$0x100		,%%rax 	/* __i2 = s1p28r */\n\t			subq	$0x100		,%%r10 		/* __i2 = s1p19r */\n\t"\
		"subq	$0x100		,%%rbx 	/* __i5 = s1p16r */\n\t			subq	$0x100		,%%r11 		/* __i5 = s1p07r */\n\t"\
		"subq	$0x100		,%%rcx 	/* __i8 = s1p04r */\n\t			subq	$0x100		,%%r12 		/* __i8 = s1p31r */\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4\n\t							vmovaps	    (%%r11)	,%%ymm12\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5\n\t							vmovaps	0x20(%%r11)	,%%ymm13\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t							vmovaps	    (%%r10)	,%%ymm8 \n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t							vmovaps	0x20(%%r10)	,%%ymm9 \n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2\n\t							vmovaps	    (%%r12)	,%%ymm10\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3\n\t							vmovaps	0x20(%%r12)	,%%ymm11\n\t"\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t							vsubpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t							vsubpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm10,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm11,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm10,%%ymm8 ,%%ymm8 \n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t							vaddpd	%%ymm11,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)	/* <- t0c */\n\t	vmovaps	%%ymm8 		,0x240(%%rsi)	/* <- t0c */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t0d */\n\t	vmovaps	%%ymm9 		,0x260(%%rsi)	/* <- t0d */\n\t"\
		"addq	$0x80		,%%rdi\n\t"\
		"addq	$0x80		,%%rsi\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t							vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t							vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t							vmulpd	%%ymm7 ,%%ymm12,%%ymm12\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t							vmulpd	%%ymm7 ,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm8 ,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm9 ,%%ymm11,%%ymm11\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t							vmovaps	%%ymm10		,%%ymm8 \n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t							vmovaps	%%ymm11		,%%ymm9 \n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2\n\t							vsubpd	%%ymm13,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm12,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm13,%%ymm8 ,%%ymm8 \n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1\n\t							vsubpd	%%ymm12,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm2		,    (%%rdi)	/* <- t0e */\n\t	vmovaps	%%ymm10		,0x240(%%rdi)	/* <- t0e */\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rdi)	/* <- t0f */\n\t	vmovaps	%%ymm11		,0x260(%%rdi)	/* <- t0f */\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)	/* <- t0g */\n\t	vmovaps	%%ymm8 		,0x240(%%rsi)	/* <- t0g */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t0h */\n\t	vmovaps	%%ymm9 		,0x260(%%rsi)	/* <- t0h */\n\t"\
		"/*****************************************/\n\t"\
		"/* now do three more radix-3 transforms: */\n\t"\
		"/*****************************************/\n\t"\
		"movq	%[__r00]	,%%rax 		/* __r00 */	\n\t"\
		"movq	%%rax		,%%rbx					\n\t"\
		"movq	%%rax		,%%rcx					\n\t"\
		"addq	$0xc0		,%%rbx 		/* __r06 */	\n\t"\
		"addq	$0x180		,%%rcx 		/* __r0c */	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4\n\t							vmovaps	0x240(%%rbx)	,%%ymm12\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5\n\t							vmovaps	0x260(%%rbx)	,%%ymm13\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2\n\t							vmovaps	0x240(%%rcx)	,%%ymm10\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3\n\t							vmovaps	0x260(%%rcx)	,%%ymm11\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t							vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t							vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t							vsubpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t							vsubpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm10,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm11,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm10,%%ymm8 ,%%ymm8 \n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t							vaddpd	%%ymm11,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm0		,    (%%rax)\n\t					vmovaps	%%ymm8 		,0x240(%%rax)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)\n\t					vmovaps	%%ymm9 		,0x260(%%rax)\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t							vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t							vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t							vmulpd	%%ymm7 ,%%ymm12,%%ymm12\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t							vmulpd	%%ymm7 ,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm8 ,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm9 ,%%ymm11,%%ymm11\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t							vmovaps	%%ymm10		,%%ymm8 \n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t							vmovaps	%%ymm11		,%%ymm9 \n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2\n\t							vsubpd	%%ymm13,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm12,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm13,%%ymm8 ,%%ymm8 \n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1\n\t							vsubpd	%%ymm12,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm2		,    (%%rbx)\n\t					vmovaps	%%ymm10		,0x240(%%rbx)\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rbx)\n\t					vmovaps	%%ymm11		,0x260(%%rbx)\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)\n\t					vmovaps	%%ymm8 		,0x240(%%rcx)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)\n\t					vmovaps	%%ymm9 		,0x260(%%rcx)\n\t"\
		"subq	$0x80		,%%rdx 		/* c1 */\n\t"\
		"addq	$0x40		,%%rax 		/* __r02 */\n\t"\
		"addq	$0x40		,%%rbx 		/* __r08 */\n\t"\
		"addq	$0x40		,%%rcx 		/* __r0e */\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2\n\t							vmovaps	0x240(%%rbx)	,%%ymm10\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3\n\t							vmovaps	0x260(%%rbx)	,%%ymm11\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
		"addq	$0x40		,%%rdx 		/* c2 */\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t							vmovaps	%%ymm10		,%%ymm8 \n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t							vmovaps	%%ymm11		,%%ymm9 \n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t							vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t							vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
		"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t							vmulpd	%%ymm7 ,%%ymm8 ,%%ymm8 \n\t"\
		"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t							vmulpd	%%ymm7 ,%%ymm9 ,%%ymm9 \n\t"\
		"vsubpd	%%ymm1,%%ymm2,%%ymm2\n\t							vsubpd	%%ymm9 ,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm0,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm8 ,%%ymm11,%%ymm11\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm4\n\t							vmovaps	0x240(%%rcx)	,%%ymm12\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm5\n\t							vmovaps	0x260(%%rcx)	,%%ymm13\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
		"vmovaps	%%ymm4		,%%ymm0\n\t							vmovaps	%%ymm12		,%%ymm8 \n\t"\
		"vmovaps	%%ymm5		,%%ymm1\n\t							vmovaps	%%ymm13		,%%ymm9 \n\t"\
		"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t							vmulpd	%%ymm6 ,%%ymm12,%%ymm12\n\t"\
		"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t							vmulpd	%%ymm6 ,%%ymm13,%%ymm13\n\t"\
		"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t							vmulpd	%%ymm7 ,%%ymm8 ,%%ymm8 \n\t"\
		"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t							vmulpd	%%ymm7 ,%%ymm9 ,%%ymm9 \n\t"\
		"vsubpd	%%ymm1,%%ymm4,%%ymm4\n\t							vsubpd	%%ymm9 ,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm8 ,%%ymm13,%%ymm13\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t							vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t							vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t							vsubpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t							vsubpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t							vaddpd	%%ymm12,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm13,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t							vaddpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm12,%%ymm8 ,%%ymm8 \n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1\n\t							vaddpd	%%ymm13,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm0		,    (%%rax)\n\t					vmovaps	%%ymm8 		,0x240(%%rax)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)\n\t					vmovaps	%%ymm9 		,0x260(%%rax)\n\t"\
		"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t							vmulpd	%%ymm6 ,%%ymm12,%%ymm12\n\t"\
		"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t							vmulpd	%%ymm6 ,%%ymm13,%%ymm13\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2\n\t							vmulpd	%%ymm7 ,%%ymm10,%%ymm10\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3\n\t							vmulpd	%%ymm7 ,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4\n\t							vaddpd	%%ymm8 ,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm9 ,%%ymm13,%%ymm13\n\t"\
		"vmovaps	%%ymm4		,%%ymm0\n\t							vmovaps	%%ymm12		,%%ymm8 \n\t"\
		"vmovaps	%%ymm5		,%%ymm1\n\t							vmovaps	%%ymm13		,%%ymm9 \n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4\n\t							vsubpd	%%ymm11,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm10,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm11,%%ymm8 ,%%ymm8 \n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1\n\t							vsubpd	%%ymm10,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)\n\t					vmovaps	%%ymm12		,0x240(%%rbx)\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rbx)\n\t					vmovaps	%%ymm13		,0x260(%%rbx)\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)\n\t					vmovaps	%%ymm8 		,0x240(%%rcx)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)\n\t					vmovaps	%%ymm9 		,0x260(%%rcx)\n\t"\
		"subq	$0x40		,%%rdx 		/* c2 */\n\t"\
		"addq	$0x40		,%%rax 		/* __r04 */\n\t"\
		"addq	$0x40		,%%rbx 		/* __r0a */\n\t"\
		"addq	$0x40		,%%rcx 		/* __r0g */\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2\n\t							vmovaps	0x240(%%rbx)	,%%ymm10\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3\n\t							vmovaps	0x260(%%rbx)	,%%ymm11\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
		"addq	$0x80		,%%rdx 		/* c4 */\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t							vmovaps	%%ymm10		,%%ymm8 \n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t							vmovaps	%%ymm11		,%%ymm9 \n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t							vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t							vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
		"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t							vmulpd	%%ymm7 ,%%ymm8 ,%%ymm8 \n\t"\
		"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t							vmulpd	%%ymm7 ,%%ymm9 ,%%ymm9 \n\t"\
		"vsubpd	%%ymm1,%%ymm2,%%ymm2\n\t							vsubpd	%%ymm9 ,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm0,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm8 ,%%ymm11,%%ymm11\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm4\n\t							vmovaps	0x240(%%rcx)	,%%ymm12\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm5\n\t							vmovaps	0x260(%%rcx)	,%%ymm13\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
		"subq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
		"vmovaps	%%ymm4		,%%ymm0\n\t							vmovaps	%%ymm12		,%%ymm8 \n\t"\
		"vmovaps	%%ymm5		,%%ymm1\n\t							vmovaps	%%ymm13		,%%ymm9 \n\t"\
		"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t							vmulpd	%%ymm6 ,%%ymm12,%%ymm12\n\t"\
		"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t							vmulpd	%%ymm6 ,%%ymm13,%%ymm13\n\t"\
		"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t							vmulpd	%%ymm7 ,%%ymm8 ,%%ymm8 \n\t"\
		"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t							vmulpd	%%ymm7 ,%%ymm9 ,%%ymm9 \n\t"\
		"vsubpd	%%ymm1,%%ymm4,%%ymm4\n\t							vsubpd	%%ymm9 ,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm8 ,%%ymm13,%%ymm13\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t							vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t							vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t							vsubpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t							vsubpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t							vaddpd	%%ymm12,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm13,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t							vaddpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm12,%%ymm8 ,%%ymm8 \n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1\n\t							vaddpd	%%ymm13,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm0		,    (%%rax)\n\t					vmovaps	%%ymm8 		,0x240(%%rax)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)\n\t					vmovaps	%%ymm9 		,0x260(%%rax)\n\t"\
		"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t							vmulpd	%%ymm6 ,%%ymm12,%%ymm12\n\t"\
		"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t							vmulpd	%%ymm6 ,%%ymm13,%%ymm13\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2\n\t							vmulpd	%%ymm7 ,%%ymm10,%%ymm10\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3\n\t							vmulpd	%%ymm7 ,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4\n\t							vaddpd	%%ymm8 ,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm9 ,%%ymm13,%%ymm13\n\t"\
		"vmovaps	%%ymm4		,%%ymm0\n\t							vmovaps	%%ymm12		,%%ymm8 \n\t"\
		"vmovaps	%%ymm5		,%%ymm1\n\t							vmovaps	%%ymm13		,%%ymm9 \n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4\n\t							vsubpd	%%ymm11,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm10,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm11,%%ymm8 ,%%ymm8 \n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1\n\t							vsubpd	%%ymm10,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)\n\t					vmovaps	%%ymm12		,0x240(%%rbx)\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rbx)\n\t					vmovaps	%%ymm13		,0x260(%%rbx)\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)\n\t					vmovaps	%%ymm8 		,0x240(%%rcx)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)\n\t					vmovaps	%%ymm9 		,0x260(%%rcx)\n\t"\
		/* SSE2_RADIX_09_DIF_0TWIDDLE(r20,s1p[18,14,10,06,02,34,30,26,22]r: */\
		"movq	%[__in0]	,%%rax\n\t"\
		"movq	%[__cc1]	,%%rdx 	\n\t"\
		"movq	%[__r00]	,%%rsi\n\t"								/* SSE2_RADIX_09_DIF_0TWIDDLE(r30,s1p[09,05,01,33,29,25,21,17,13]r: */\
		"addq	$0x480		,%%rsi 		/* r20 */\n\t				leaq	0x240(%%rax),%%r10 		/* __i0 = s1p09r */\n\t"\
		"addq	$0x480		,%%rax 	/* __i0 = s1p18r */\n\t"\
		"movq	%%rax		,%%rbx\n\t								movq	%%r10		,%%r11\n\t"\
		"movq	%%rax		,%%rcx\n\t								movq	%%r10		,%%r12\n\t"\
		"subq	$0x300		,%%rbx 	/* __i3 = s1p06r */\n\t			addq	$0x600		,%%r11 		/* __i3 = s1p33r */\n\t"\
		"addq	$0x300		,%%rcx 	/* __i6 = s1p30r */\n\t			addq	$0x300		,%%r12 		/* __i6 = s1p21r */\n\t"\
		"addq	$0x80		,%%rdx 		/* c3m1 */\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2\n\t							vmovaps	    (%%r11)	,%%ymm10\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3\n\t							vmovaps	0x20(%%r11)	,%%ymm11\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t							vmovaps	    (%%r10)	,%%ymm8 \n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t							vmovaps	0x20(%%r10)	,%%ymm9 \n\t"\
		"vmovaps	    (%%rcx)	,%%ymm6\n\t							vmovaps	    (%%r12)	,%%ymm14\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm7\n\t							vmovaps	0x20(%%r12)	,%%ymm15\n\t"\
		"vmovaps	%%ymm2		,%%ymm4\n\t							vmovaps	%%ymm10		,%%ymm12\n\t"\
		"vmovaps	%%ymm3		,%%ymm5\n\t							vmovaps	%%ymm11		,%%ymm13\n\t"\
		"movq	%%rsi		,%%rdi\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm14,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm15,%%ymm11,%%ymm11\n\t"\
		"vsubpd	%%ymm6,%%ymm4,%%ymm4\n\t							vsubpd	%%ymm14,%%ymm12,%%ymm12\n\t"\
		"vsubpd	%%ymm7,%%ymm5,%%ymm5\n\t							vsubpd	%%ymm15,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm10,%%ymm8 ,%%ymm8 \n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t							vaddpd	%%ymm11,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)	/* <- t00 */\n\t	vmovaps	%%ymm8 		,0x240(%%rsi)	/* <- t00 */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t01 */\n\t	vmovaps	%%ymm9 		,0x260(%%rsi)	/* <- t01 */\n\t"\
		"addq	$0x80		,%%rsi\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t							vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t							vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t							vmulpd	%%ymm7 ,%%ymm12,%%ymm12\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t							vmulpd	%%ymm7 ,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm8 ,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm9 ,%%ymm11,%%ymm11\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t							vmovaps	%%ymm10		,%%ymm8 \n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t							vmovaps	%%ymm11		,%%ymm9 \n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2\n\t							vsubpd	%%ymm13,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm12,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm13,%%ymm8 ,%%ymm8 \n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1\n\t							vsubpd	%%ymm12,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm2		,    (%%rdi)	/* <- t02 */\n\t	vmovaps	%%ymm10		,0x240(%%rdi)	/* <- t02 */\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rdi)	/* <- t03 */\n\t	vmovaps	%%ymm11		,0x260(%%rdi)	/* <- t03 */\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)	/* <- t04 */\n\t	vmovaps	%%ymm8 		,0x240(%%rsi)	/* <- t04 */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t05 */\n\t	vmovaps	%%ymm9 		,0x260(%%rsi)	/* <- t05 */\n\t"\
		"subq	$0x100		,%%rax 	/* __i1 = s1p14r */\n\t			subq	$0x100		,%%r10 		/* __i1 = s1p05r */\n\t"\
		"subq	$0x100		,%%rbx 	/* __i4 = s1p02r */\n\t			subq	$0x100		,%%r11 		/* __i4 = s1p29r */\n\t"\
		"subq	$0x100		,%%rcx 	/* __i7 = s1p26r */\n\t			subq	$0x100		,%%r12 		/* __i7 = s1p17r */\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4\n\t							vmovaps	    (%%r11)	,%%ymm12\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5\n\t							vmovaps	0x20(%%r11)	,%%ymm13\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t							vmovaps	    (%%r10)	,%%ymm8 \n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t							vmovaps	0x20(%%r10)	,%%ymm9 \n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2\n\t							vmovaps	    (%%r12)	,%%ymm10\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3\n\t							vmovaps	0x20(%%r12)	,%%ymm11\n\t"\
		"addq	$0x80		,%%rdi\n\t"\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t							vsubpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t							vsubpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm10,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm11,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm10,%%ymm8 ,%%ymm8 \n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t							vaddpd	%%ymm11,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm0		,    (%%rdi)	/* <- t06 */\n\t	vmovaps	%%ymm8 		,0x240(%%rdi)	/* <- t06 */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rdi)	/* <- t07 */\n\t	vmovaps	%%ymm9 		,0x260(%%rdi)	/* <- t07 */\n\t"\
		"addq	$0x80		,%%rsi\n\t"\
		"addq	$0x80		,%%rdi\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t							vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t							vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t							vmulpd	%%ymm7 ,%%ymm12,%%ymm12\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t							vmulpd	%%ymm7 ,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm8 ,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm9 ,%%ymm11,%%ymm11\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t							vmovaps	%%ymm10		,%%ymm8 \n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t							vmovaps	%%ymm11		,%%ymm9 \n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2\n\t							vsubpd	%%ymm13,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm12,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm13,%%ymm8 ,%%ymm8 \n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1\n\t							vsubpd	%%ymm12,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm2		,    (%%rsi)	/* <- t08 */\n\t	vmovaps	%%ymm10		,0x240(%%rsi)	/* <- t08 */\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rsi)	/* <- t09 */\n\t	vmovaps	%%ymm11		,0x260(%%rsi)	/* <- t09 */\n\t"\
		"vmovaps	%%ymm0		,    (%%rdi)	/* <- t0a */\n\t	vmovaps	%%ymm8 		,0x240(%%rdi)	/* <- t0a */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rdi)	/* <- t0b */\n\t	vmovaps	%%ymm9 		,0x260(%%rdi)	/* <- t0b */\n\t"\
		"addq	$0x80		,%%rsi\n\t"\
		"subq	$0x100		,%%rax 	/* __i2 = s1p10r */\n\t			subq	$0x100		,%%r10 		/* __i2 = s1p01r */\n\t"\
		"addq	$0x800		,%%rbx 	/* __i5 = s1p34r */\n\t			subq	$0x100		,%%r11 		/* __i5 = s1p25r */\n\t"\
		"subq	$0x100		,%%rcx 	/* __i8 = s1p22r */\n\t			subq	$0x100		,%%r12 		/* __i8 = s1p13r */\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4\n\t							vmovaps	    (%%r11)	,%%ymm12\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5\n\t							vmovaps	0x20(%%r11)	,%%ymm13\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t							vmovaps	    (%%r10)	,%%ymm8 \n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t							vmovaps	0x20(%%r10)	,%%ymm9 \n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2\n\t							vmovaps	    (%%r12)	,%%ymm10\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3\n\t							vmovaps	0x20(%%r12)	,%%ymm11\n\t"\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t							vsubpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t							vsubpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm10,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm11,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm10,%%ymm8 ,%%ymm8 \n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t							vaddpd	%%ymm11,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)	/* <- t0c */\n\t	vmovaps	%%ymm8 		,0x240(%%rsi)	/* <- t0c */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t0d */\n\t	vmovaps	%%ymm9 		,0x260(%%rsi)	/* <- t0d */\n\t"\
		"addq	$0x80		,%%rdi\n\t"\
		"addq	$0x80		,%%rsi\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t							vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t							vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t							vmulpd	%%ymm7 ,%%ymm12,%%ymm12\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t							vmulpd	%%ymm7 ,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm8 ,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm9 ,%%ymm11,%%ymm11\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t							vmovaps	%%ymm10		,%%ymm8 \n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t							vmovaps	%%ymm11		,%%ymm9 \n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2\n\t							vsubpd	%%ymm13,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm12,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm13,%%ymm8 ,%%ymm8 \n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1\n\t							vsubpd	%%ymm12,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm2		,    (%%rdi)	/* <- t0e */\n\t	vmovaps	%%ymm10		,0x240(%%rdi)	/* <- t0e */\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rdi)	/* <- t0f */\n\t	vmovaps	%%ymm11		,0x260(%%rdi)	/* <- t0f */\n\t"\
		"vmovaps	%%ymm0		,    (%%rsi)	/* <- t0g */\n\t	vmovaps	%%ymm8 		,0x240(%%rsi)	/* <- t0g */\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rsi)	/* <- t0h */\n\t	vmovaps	%%ymm9 		,0x260(%%rsi)	/* <- t0h */\n\t"\
		"/*****************************************/\n\t"\
		"/* now do three more radix-3 transforms: */\n\t"\
		"/*****************************************/\n\t"\
		"movq	%[__r00]	,%%rax 		/* __r00 */	\n\t"\
		"addq	$0x480		,%%rax 		/* r20 */	\n\t"\
		"movq	%%rax		,%%rbx					\n\t"\
		"movq	%%rax		,%%rcx					\n\t"\
		"addq	$0xc0		,%%rbx 		/* __r26 */	\n\t"\
		"addq	$0x180		,%%rcx 		/* __r2c */	\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm4\n\t							vmovaps	0x240(%%rbx)	,%%ymm12\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm5\n\t							vmovaps	0x260(%%rbx)	,%%ymm13\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm2\n\t							vmovaps	0x240(%%rcx)	,%%ymm10\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm3\n\t							vmovaps	0x260(%%rcx)	,%%ymm11\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t							vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t							vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
		"vsubpd	%%ymm2,%%ymm4,%%ymm4\n\t							vsubpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
		"vsubpd	%%ymm3,%%ymm5,%%ymm5\n\t							vsubpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm10,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm11,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm4,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm5,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm2,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm10,%%ymm8 ,%%ymm8 \n\t"\
		"vaddpd	%%ymm3,%%ymm1,%%ymm1\n\t							vaddpd	%%ymm11,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm0		,    (%%rax)\n\t					vmovaps	%%ymm8 		,0x240(%%rax)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)\n\t					vmovaps	%%ymm9 		,0x260(%%rax)\n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t							vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t							vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
		"vmulpd	%%ymm7,%%ymm4,%%ymm4\n\t							vmulpd	%%ymm7 ,%%ymm12,%%ymm12\n\t"\
		"vmulpd	%%ymm7,%%ymm5,%%ymm5\n\t							vmulpd	%%ymm7 ,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2\n\t							vaddpd	%%ymm8 ,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm9 ,%%ymm11,%%ymm11\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t							vmovaps	%%ymm10		,%%ymm8 \n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t							vmovaps	%%ymm11		,%%ymm9 \n\t"\
		"vsubpd	%%ymm5,%%ymm2,%%ymm2\n\t							vsubpd	%%ymm13,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm4,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm12,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm5,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm13,%%ymm8 ,%%ymm8 \n\t"\
		"vsubpd	%%ymm4,%%ymm1,%%ymm1\n\t							vsubpd	%%ymm12,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm2		,    (%%rbx)\n\t					vmovaps	%%ymm10		,0x240(%%rbx)\n\t"\
		"vmovaps	%%ymm3		,0x20(%%rbx)\n\t					vmovaps	%%ymm11		,0x260(%%rbx)\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)\n\t					vmovaps	%%ymm8 		,0x240(%%rcx)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)\n\t					vmovaps	%%ymm9 		,0x260(%%rcx)\n\t"\
		"subq	$0x80		,%%rdx 		/* c1 */\n\t"\
		"addq	$0x40		,%%rax 		/* __r22 */\n\t"\
		"addq	$0x40		,%%rbx 		/* __r28 */\n\t"\
		"addq	$0x40		,%%rcx 		/* __r2e */\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2\n\t							vmovaps	0x240(%%rbx)	,%%ymm10\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3\n\t							vmovaps	0x260(%%rbx)	,%%ymm11\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6\n\t							vmovaps	    (%%rdx)	,%%ymm14\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm7\n\t							vmovaps	0x20(%%rdx)	,%%ymm15\n\t"\
		"addq	$0x40		,%%rdx 		/* c2 */\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t							vmovaps	%%ymm10		,%%ymm8 \n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t							vmovaps	%%ymm11		,%%ymm9 \n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t							vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t							vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
		"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t							vmulpd	%%ymm7 ,%%ymm8 ,%%ymm8 \n\t"\
		"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t							vmulpd	%%ymm7 ,%%ymm9 ,%%ymm9 \n\t"\
		"vsubpd	%%ymm1,%%ymm2,%%ymm2\n\t							vsubpd	%%ymm9 ,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm0,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm8 ,%%ymm11,%%ymm11\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm4\n\t							vmovaps	0x240(%%rcx)	,%%ymm12\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm5\n\t							vmovaps	0x260(%%rcx)	,%%ymm13\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
		"vmovaps	%%ymm4		,%%ymm0\n\t							vmovaps	%%ymm12		,%%ymm8 \n\t"\
		"vmovaps	%%ymm5		,%%ymm1\n\t							vmovaps	%%ymm13		,%%ymm9 \n\t"\
		"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t							vmulpd	%%ymm6 ,%%ymm12,%%ymm12\n\t"\
		"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t							vmulpd	%%ymm6 ,%%ymm13,%%ymm13\n\t"\
		"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t							vmulpd	%%ymm7 ,%%ymm8 ,%%ymm8 \n\t"\
		"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t							vmulpd	%%ymm7 ,%%ymm9 ,%%ymm9 \n\t"\
		"vsubpd	%%ymm1,%%ymm4,%%ymm4\n\t							vsubpd	%%ymm9 ,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm8 ,%%ymm13,%%ymm13\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t							vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t							vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t							vsubpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t							vsubpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t							vaddpd	%%ymm12,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm13,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t							vaddpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm12,%%ymm8 ,%%ymm8 \n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1\n\t							vaddpd	%%ymm13,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm0		,    (%%rax)\n\t					vmovaps	%%ymm8 		,0x240(%%rax)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)\n\t					vmovaps	%%ymm9 		,0x260(%%rax)\n\t"\
		"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t							vmulpd	%%ymm6 ,%%ymm12,%%ymm12\n\t"\
		"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t							vmulpd	%%ymm6 ,%%ymm13,%%ymm13\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2\n\t							vmulpd	%%ymm7 ,%%ymm10,%%ymm10\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3\n\t							vmulpd	%%ymm7 ,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4\n\t							vaddpd	%%ymm8 ,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm9 ,%%ymm13,%%ymm13\n\t"\
		"vmovaps	%%ymm4		,%%ymm0\n\t							vmovaps	%%ymm12		,%%ymm8 \n\t"\
		"vmovaps	%%ymm5		,%%ymm1\n\t							vmovaps	%%ymm13		,%%ymm9 \n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4\n\t							vsubpd	%%ymm11,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm10,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm11,%%ymm8 ,%%ymm8 \n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1\n\t							vsubpd	%%ymm10,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)\n\t					vmovaps	%%ymm12		,0x240(%%rbx)\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rbx)\n\t					vmovaps	%%ymm13		,0x260(%%rbx)\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)\n\t					vmovaps	%%ymm8 		,0x240(%%rcx)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)\n\t					vmovaps	%%ymm9 		,0x260(%%rcx)\n\t"\
		"subq	$0x40		,%%rdx 		/* c2 */\n\t"\
		"addq	$0x40		,%%rax 		/* __r24 */\n\t"\
		"addq	$0x40		,%%rbx 		/* __r2a */\n\t"\
		"addq	$0x40		,%%rcx 		/* __r2g */\n\t"\
		"vmovaps	    (%%rbx)	,%%ymm2\n\t							vmovaps	0x240(%%rbx)	,%%ymm10\n\t"\
		"vmovaps	0x20(%%rbx)	,%%ymm3\n\t							vmovaps	0x260(%%rbx)	,%%ymm11\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
		"addq	$0x80		,%%rdx 		/* c4 */\n\t"\
		"vmovaps	%%ymm2		,%%ymm0\n\t							vmovaps	%%ymm10		,%%ymm8 \n\t"\
		"vmovaps	%%ymm3		,%%ymm1\n\t							vmovaps	%%ymm11		,%%ymm9 \n\t"\
		"vmulpd	%%ymm6,%%ymm2,%%ymm2\n\t							vmulpd	%%ymm6 ,%%ymm10,%%ymm10\n\t"\
		"vmulpd	%%ymm6,%%ymm3,%%ymm3\n\t							vmulpd	%%ymm6 ,%%ymm11,%%ymm11\n\t"\
		"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t							vmulpd	%%ymm7 ,%%ymm8 ,%%ymm8 \n\t"\
		"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t							vmulpd	%%ymm7 ,%%ymm9 ,%%ymm9 \n\t"\
		"vsubpd	%%ymm1,%%ymm2,%%ymm2\n\t							vsubpd	%%ymm9 ,%%ymm10,%%ymm10\n\t"\
		"vaddpd	%%ymm0,%%ymm3,%%ymm3\n\t							vaddpd	%%ymm8 ,%%ymm11,%%ymm11\n\t"\
		"vmovaps	    (%%rcx)	,%%ymm4\n\t							vmovaps	0x240(%%rcx)	,%%ymm12\n\t"\
		"vmovaps	0x20(%%rcx)	,%%ymm5\n\t							vmovaps	0x260(%%rcx)	,%%ymm13\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
		"subq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
		"vmovaps	%%ymm4		,%%ymm0\n\t							vmovaps	%%ymm12		,%%ymm8 \n\t"\
		"vmovaps	%%ymm5		,%%ymm1\n\t							vmovaps	%%ymm13		,%%ymm9 \n\t"\
		"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t							vmulpd	%%ymm6 ,%%ymm12,%%ymm12\n\t"\
		"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t							vmulpd	%%ymm6 ,%%ymm13,%%ymm13\n\t"\
		"vmulpd	%%ymm7,%%ymm0,%%ymm0\n\t							vmulpd	%%ymm7 ,%%ymm8 ,%%ymm8 \n\t"\
		"vmulpd	%%ymm7,%%ymm1,%%ymm1\n\t							vmulpd	%%ymm7 ,%%ymm9 ,%%ymm9 \n\t"\
		"vsubpd	%%ymm1,%%ymm4,%%ymm4\n\t							vsubpd	%%ymm9 ,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm0,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm8 ,%%ymm13,%%ymm13\n\t"\
		"vmovaps	    (%%rdx)	,%%ymm6\n\t"\
		"vmovaps	0x20(%%rdx)	,%%ymm7\n\t"\
		"vmovaps	    (%%rax)	,%%ymm0\n\t							vmovaps	0x240(%%rax)	,%%ymm8 \n\t"\
		"vmovaps	0x20(%%rax)	,%%ymm1\n\t							vmovaps	0x260(%%rax)	,%%ymm9 \n\t"\
		"vsubpd	%%ymm4,%%ymm2,%%ymm2\n\t							vsubpd	%%ymm12,%%ymm10,%%ymm10\n\t"\
		"vsubpd	%%ymm5,%%ymm3,%%ymm3\n\t							vsubpd	%%ymm13,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm4,%%ymm4,%%ymm4\n\t							vaddpd	%%ymm12,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm5,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm13,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm2,%%ymm4,%%ymm4\n\t							vaddpd	%%ymm10,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm3,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm11,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm12,%%ymm8 ,%%ymm8 \n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1\n\t							vaddpd	%%ymm13,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm0		,    (%%rax)\n\t					vmovaps	%%ymm8 		,0x240(%%rax)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rax)\n\t					vmovaps	%%ymm9 		,0x260(%%rax)\n\t"\
		"vmulpd	%%ymm6,%%ymm4,%%ymm4\n\t							vmulpd	%%ymm6 ,%%ymm12,%%ymm12\n\t"\
		"vmulpd	%%ymm6,%%ymm5,%%ymm5\n\t							vmulpd	%%ymm6 ,%%ymm13,%%ymm13\n\t"\
		"vmulpd	%%ymm7,%%ymm2,%%ymm2\n\t							vmulpd	%%ymm7 ,%%ymm10,%%ymm10\n\t"\
		"vmulpd	%%ymm7,%%ymm3,%%ymm3\n\t							vmulpd	%%ymm7 ,%%ymm11,%%ymm11\n\t"\
		"vaddpd	%%ymm0,%%ymm4,%%ymm4\n\t							vaddpd	%%ymm8 ,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm1,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm9 ,%%ymm13,%%ymm13\n\t"\
		"vmovaps	%%ymm4		,%%ymm0\n\t							vmovaps	%%ymm12		,%%ymm8 \n\t"\
		"vmovaps	%%ymm5		,%%ymm1\n\t							vmovaps	%%ymm13		,%%ymm9 \n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4\n\t							vsubpd	%%ymm11,%%ymm12,%%ymm12\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5\n\t							vaddpd	%%ymm10,%%ymm13,%%ymm13\n\t"\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0\n\t							vaddpd	%%ymm11,%%ymm8 ,%%ymm8 \n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1\n\t							vsubpd	%%ymm10,%%ymm9 ,%%ymm9 \n\t"\
		"vmovaps	%%ymm4		,    (%%rbx)\n\t					vmovaps	%%ymm12		,0x240(%%rbx)\n\t"\
		"vmovaps	%%ymm5		,0x20(%%rbx)\n\t					vmovaps	%%ymm13		,0x260(%%rbx)\n\t"\
		"vmovaps	%%ymm0		,    (%%rcx)\n\t					vmovaps	%%ymm8 		,0x240(%%rcx)\n\t"\
		"vmovaps	%%ymm1		,0x20(%%rcx)\n\t					vmovaps	%%ymm9 		,0x260(%%rcx)\n\t"\
		"\n\t"\
		"/**********************************/"\
		"/*** And now do 9 radix-4 DFTs: ***/"\
		"/**********************************/"\
		"\n\t"\
		"movq	%[__add],%%rax	/* Must use eax as base address throughout, since that is preserved in SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_B */\n\t"\
		"movq	%[__add],%%rdx						\n\t"\
		/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r00, 0x240, 0x480, eax,ebx,edx,ecx) */\
		"movq	%[__r00],%%rdi 						\n\t"\
		"vmovaps	     (%%rdi),%%ymm4				\n\t"\
		"vmovaps	0x240(%%rdi),%%ymm6				\n\t"\
		"vmovaps	0x020(%%rdi),%%ymm5				\n\t"\
		"vmovaps	0x260(%%rdi),%%ymm7				\n\t"\
		"vmovaps	0x480(%%rdi),%%ymm0				\n\t"\
		"vmovaps	0x6c0(%%rdi),%%ymm2				\n\t"\
		"vmovaps	0x4a0(%%rdi),%%ymm1				\n\t"\
		"vmovaps	0x6e0(%%rdi),%%ymm3				\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4				\n\t"\
		"vsubpd	%%ymm2,%%ymm6,%%ymm6				\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5				\n\t"\
		"vsubpd	%%ymm3,%%ymm7,%%ymm7				\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0				\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2				\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1				\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3				\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0				\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2				\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1				\n\t"\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3				\n\t"\
		"/* Swap output addresses d <-> c: */		\n\t"						/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r02, abdc) */\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0				\n\t						addq	$0x40	,%%rdi 	/* r02 */	\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4				\n\t						vmovaps	     (%%rdi),%%ymm12		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1				\n\t						vmovaps	0x240(%%rdi),%%ymm14		\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5				\n\t						vmovaps	0x020(%%rdi),%%ymm13		\n\t"\
		/* add0,1,3,2 = a+p00+p[0,1,2,3] */			"\n\t						vmovaps	0x260(%%rdi),%%ymm15		\n\t"\
		"movslq	%[__p01],%%rsi	/* rsi stores offset-multiples of p04  */\n\t	vmovaps	0x480(%%rdi),%%ymm8 		\n\t"\
		"shlq	$3,%%rsi		/* Pointer offset for floating doubles */\n\t	vmovaps	0x6c0(%%rdi),%%ymm10		\n\t"\
		"addq	%%rsi,%%rdx							\n\t						vmovaps	0x4a0(%%rdi),%%ymm9 		\n\t"\
		"movq	%%rdx,%%rbx		/* add0+p01 */		\n\t						vmovaps	0x6e0(%%rdi),%%ymm11		\n\t"\
		"addq	%%rsi,%%rdx							\n\t						vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"movq	%%rdx,%%rcx		/* add0+p02 */		\n\t						vsubpd	%%ymm10,%%ymm14,%%ymm14		\n\t"\
		"addq	%%rsi,%%rdx		/* add0+p03 */		\n\t						vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,    (%%rbx)				\n\t						vsubpd	%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,    (%%rdx)				\n\t						vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm1,0x20(%%rbx)				\n\t						vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	%%ymm5,0x20(%%rcx)				\n\t						vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2				\n\t						vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7				\n\t						vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3				\n\t						vaddpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6				\n\t						vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2				\n\t						vaddpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7				\n\t						/* Swap output addresses d <-> c: */\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3				\n\t						vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6				\n\t						vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm2,    (%%rax)				\n\t						vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm7,    (%%rcx)				\n\t						vsubpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm3,0x20(%%rax)				\n\t"						/* add0,1,3,2 = a+p32+p[0,1,2,3] */\
		"vmovaps	%%ymm6,0x20(%%rdx)				\n\t						movslq	%[__p32],%%rsi				\n\t"\
		/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r04, cdab) */"\n\t					shlq	$3,%%rsi					\n\t"\
		"addq	$0x40	,%%rdi 	/* r04 */			\n\t						addq	%%rsi,%%rax		/* a+p32] */\n\t"\
		"vmovaps	     (%%rdi),%%ymm4				\n\t						addq	%%rsi,%%rbx					\n\t"\
		"vmovaps	0x240(%%rdi),%%ymm6				\n\t						addq	%%rsi,%%rcx					\n\t"\
		"vmovaps	0x020(%%rdi),%%ymm5				\n\t						addq	%%rsi,%%rdx					\n\t"\
		"vmovaps	0x260(%%rdi),%%ymm7				\n\t						vmovaps	%%ymm8 ,    (%%rbx)			\n\t"\
		"vmovaps	0x480(%%rdi),%%ymm0				\n\t						vmovaps	%%ymm12,    (%%rdx)			\n\t"\
		"vmovaps	0x6c0(%%rdi),%%ymm2				\n\t						vmovaps	%%ymm9 ,0x20(%%rbx)			\n\t"\
		"vmovaps	0x4a0(%%rdi),%%ymm1				\n\t						vmovaps	%%ymm13,0x20(%%rcx)			\n\t"\
		"vmovaps	0x6e0(%%rdi),%%ymm3				\n\t						vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4				\n\t						vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm2,%%ymm6,%%ymm6				\n\t						vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5				\n\t						vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm3,%%ymm7,%%ymm7				\n\t						vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0				\n\t						vaddpd	%%ymm12,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2				\n\t						vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1				\n\t						vaddpd	%%ymm13,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3				\n\t						vmovaps	%%ymm10,    (%%rax)			\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0				\n\t						vmovaps	%%ymm15,    (%%rcx)			\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2				\n\t						vmovaps	%%ymm11,0x20(%%rax)			\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1				\n\t						vmovaps	%%ymm14,0x20(%%rdx)			\n\t"\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3				\n\t"						/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r06, bacd) */\
		"/* Swap output addresses a <-> c ,b <-> d: */\n\t						addq	$0x40	,%%rdi 	/* r06 */	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0				\n\t						vmovaps	     (%%rdi),%%ymm12		\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4				\n\t						vmovaps	0x240(%%rdi),%%ymm14		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1				\n\t						vmovaps	0x020(%%rdi),%%ymm13		\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5				\n\t						vmovaps	0x260(%%rdi),%%ymm15		\n\t"\
		/* add2,3,0,1 = a+p20+p[0,1,2,3] */			"\n\t						vmovaps	0x480(%%rdi),%%ymm8 		\n\t"\
		"movslq	%[__p12],%%rsi						\n\t						vmovaps	0x6c0(%%rdi),%%ymm10		\n\t"\
		"shlq	$3,%%rsi							\n\t						vmovaps	0x4a0(%%rdi),%%ymm9 		\n\t"\
		"subq	%%rsi,%%rax		/* a+p20] */		\n\t						vmovaps	0x6e0(%%rdi),%%ymm11		\n\t"\
		"subq	%%rsi,%%rbx							\n\t						vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"subq	%%rsi,%%rcx							\n\t						vsubpd	%%ymm10,%%ymm14,%%ymm14		\n\t"\
		"subq	%%rsi,%%rdx							\n\t						vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,    (%%rdx)				\n\t						vsubpd	%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,    (%%rax)				\n\t						vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm1,0x20(%%rdx)				\n\t						vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	%%ymm5,0x20(%%rbx)				\n\t						vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2				\n\t						vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7				\n\t						vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3				\n\t						vaddpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6				\n\t						vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2				\n\t						vaddpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7				\n\t						/* Swap output addresses a <-> b: */\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3				\n\t						vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6				\n\t						vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm2,    (%%rcx)				\n\t						vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm7,    (%%rbx)				\n\t						vsubpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm3,0x20(%%rcx)				\n\t"						/* add1,0,2,3 = a+p08+p[0,1,2,3] */\
		"vmovaps	%%ymm6,0x20(%%rax)				\n\t						movslq	%[__p12],%%rsi				\n\t"\
		/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r08, dcba) */"\n\t					shlq	$3,%%rsi					\n\t"\
		"addq	$0x40	,%%rdi 	/* r08 */			\n\t						subq	%%rsi,%%rax		/* a+p08] */\n\t"\
		"vmovaps	     (%%rdi),%%ymm4				\n\t						subq	%%rsi,%%rbx					\n\t"\
		"vmovaps	0x240(%%rdi),%%ymm6				\n\t						subq	%%rsi,%%rcx					\n\t"\
		"vmovaps	0x020(%%rdi),%%ymm5				\n\t						subq	%%rsi,%%rdx					\n\t"\
		"vmovaps	0x260(%%rdi),%%ymm7				\n\t						vmovaps	%%ymm8 ,    (%%rax)			\n\t"\
		"vmovaps	0x480(%%rdi),%%ymm0				\n\t						vmovaps	%%ymm12,    (%%rcx)			\n\t"\
		"vmovaps	0x6c0(%%rdi),%%ymm2				\n\t						vmovaps	%%ymm9 ,0x20(%%rax)			\n\t"\
		"vmovaps	0x4a0(%%rdi),%%ymm1				\n\t						vmovaps	%%ymm13,0x20(%%rdx)			\n\t"\
		"vmovaps	0x6e0(%%rdi),%%ymm3				\n\t						vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4				\n\t						vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm2,%%ymm6,%%ymm6				\n\t						vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5				\n\t						vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm3,%%ymm7,%%ymm7				\n\t						vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0				\n\t						vaddpd	%%ymm12,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2				\n\t						vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1				\n\t						vaddpd	%%ymm13,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3				\n\t						vmovaps	%%ymm10,    (%%rbx)			\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0				\n\t						vmovaps	%%ymm15,    (%%rdx)			\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2				\n\t						vmovaps	%%ymm11,0x20(%%rbx)			\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1				\n\t						vmovaps	%%ymm14,0x20(%%rcx)			\n\t"\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3				\n\t"						/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0a, abdc) */\
		"/* Swap output addresses a <-> d ,b <-> c: */\n\t						addq	$0x40	,%%rdi 	/* r0a */	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0				\n\t						vmovaps	     (%%rdi),%%ymm12		\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4				\n\t						vmovaps	0x240(%%rdi),%%ymm14		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1				\n\t						vmovaps	0x020(%%rdi),%%ymm13		\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5				\n\t						vmovaps	0x260(%%rdi),%%ymm15		\n\t"\
		/* add3,2,1,0 = a+p28+p[0,1,2,3] */			"\n\t						vmovaps	0x480(%%rdi),%%ymm8 		\n\t"\
		"movslq	%[__p20],%%rsi						\n\t						vmovaps	0x6c0(%%rdi),%%ymm10		\n\t"\
		"shlq	$3,%%rsi							\n\t						vmovaps	0x4a0(%%rdi),%%ymm9 		\n\t"\
		"addq	%%rsi,%%rax		/* a+p28] */		\n\t						vmovaps	0x6e0(%%rdi),%%ymm11		\n\t"\
		"addq	%%rsi,%%rbx							\n\t						vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"addq	%%rsi,%%rcx							\n\t						vsubpd	%%ymm10,%%ymm14,%%ymm14		\n\t"\
		"addq	%%rsi,%%rdx							\n\t						vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,    (%%rcx)				\n\t						vsubpd	%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,    (%%rbx)				\n\t						vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm1,0x20(%%rcx)				\n\t						vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	%%ymm5,0x20(%%rax)				\n\t						vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2				\n\t						vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7				\n\t						vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3				\n\t						vaddpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6				\n\t						vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2				\n\t						vaddpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7				\n\t						/* Swap output addresses d <-> c: */\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3				\n\t						vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6				\n\t						vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm2,    (%%rdx)				\n\t						vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm7,    (%%rax)				\n\t						vsubpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm3,0x20(%%rdx)				\n\t"						/* add0,1,3,2 = a+p16+p[0,1,2,3] */\
		"vmovaps	%%ymm6,0x20(%%rbx)				\n\t						movslq	%[__p12],%%rsi				\n\t"\
		/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0c, cdab) */"\n\t					shlq	$3,%%rsi					\n\t"\
		"addq	$0x40	,%%rdi 	/* r0c */			\n\t						subq	%%rsi,%%rax		/* a+p16] */\n\t"\
		"vmovaps	     (%%rdi),%%ymm4				\n\t						subq	%%rsi,%%rbx					\n\t"\
		"vmovaps	0x240(%%rdi),%%ymm6				\n\t						subq	%%rsi,%%rcx					\n\t"\
		"vmovaps	0x020(%%rdi),%%ymm5				\n\t						subq	%%rsi,%%rdx					\n\t"\
		"vmovaps	0x260(%%rdi),%%ymm7				\n\t						vmovaps	%%ymm8 ,    (%%rbx)			\n\t"\
		"vmovaps	0x480(%%rdi),%%ymm0				\n\t						vmovaps	%%ymm12,    (%%rdx)			\n\t"\
		"vmovaps	0x6c0(%%rdi),%%ymm2				\n\t						vmovaps	%%ymm9 ,0x20(%%rbx)			\n\t"\
		"vmovaps	0x4a0(%%rdi),%%ymm1				\n\t						vmovaps	%%ymm13,0x20(%%rcx)			\n\t"\
		"vmovaps	0x6e0(%%rdi),%%ymm3				\n\t						vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4				\n\t						vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm2,%%ymm6,%%ymm6				\n\t						vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5				\n\t						vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm3,%%ymm7,%%ymm7				\n\t						vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0				\n\t						vaddpd	%%ymm12,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2				\n\t						vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1				\n\t						vaddpd	%%ymm13,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3				\n\t						vmovaps	%%ymm10,    (%%rax)			\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0				\n\t						vmovaps	%%ymm15,    (%%rcx)			\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2				\n\t						vmovaps	%%ymm11,0x20(%%rax)			\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1				\n\t						vmovaps	%%ymm14,0x20(%%rdx)			\n\t"\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3				\n\t"						/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0e, bacd) */\
		"/* Swap output addresses a <-> c ,b <-> d: */\n\t						addq	$0x40	,%%rdi 	/* r0e */	\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0				\n\t						vmovaps	     (%%rdi),%%ymm12		\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4				\n\t						vmovaps	0x240(%%rdi),%%ymm14		\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1				\n\t						vmovaps	0x020(%%rdi),%%ymm13		\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5				\n\t						vmovaps	0x260(%%rdi),%%ymm15		\n\t"\
		/* add2,3,0,1 = a+p04+p[0,1,2,3] */			"\n\t						vmovaps	0x480(%%rdi),%%ymm8 		\n\t"\
		"movslq	%[__p12],%%rsi						\n\t						vmovaps	0x6c0(%%rdi),%%ymm10		\n\t"\
		"shlq	$3,%%rsi							\n\t						vmovaps	0x4a0(%%rdi),%%ymm9 		\n\t"\
		"subq	%%rsi,%%rax		/* a+p04] */		\n\t						vmovaps	0x6e0(%%rdi),%%ymm11		\n\t"\
		"subq	%%rsi,%%rbx							\n\t						vsubpd	%%ymm8 ,%%ymm12,%%ymm12		\n\t"\
		"subq	%%rsi,%%rcx							\n\t						vsubpd	%%ymm10,%%ymm14,%%ymm14		\n\t"\
		"subq	%%rsi,%%rdx							\n\t						vsubpd	%%ymm9 ,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm0,    (%%rdx)				\n\t						vsubpd	%%ymm11,%%ymm15,%%ymm15		\n\t"\
		"vmovaps	%%ymm4,    (%%rax)				\n\t						vaddpd	%%ymm8 ,%%ymm8 ,%%ymm8 		\n\t"\
		"vmovaps	%%ymm1,0x20(%%rdx)				\n\t						vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vmovaps	%%ymm5,0x20(%%rbx)				\n\t						vaddpd	%%ymm9 ,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2				\n\t						vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7				\n\t						vaddpd	%%ymm12,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3				\n\t						vaddpd	%%ymm14,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6				\n\t						vaddpd	%%ymm13,%%ymm9 ,%%ymm9 		\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2				\n\t						vaddpd	%%ymm15,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7				\n\t						/* Swap output addresses a <-> b: */\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3				\n\t						vsubpd	%%ymm10,%%ymm8 ,%%ymm8 		\n\t"\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6				\n\t						vsubpd	%%ymm15,%%ymm12,%%ymm12		\n\t"\
		"vmovaps	%%ymm2,    (%%rcx)				\n\t						vsubpd	%%ymm11,%%ymm9 ,%%ymm9 		\n\t"\
		"vmovaps	%%ymm7,    (%%rbx)				\n\t						vsubpd	%%ymm14,%%ymm13,%%ymm13		\n\t"\
		"vmovaps	%%ymm3,0x20(%%rcx)				\n\t"						/* add1,0,2,3 = a+p24+p[0,1,2,3] */\
		"vmovaps	%%ymm6,0x20(%%rax)				\n\t						movslq	%[__p20],%%rsi				\n\t"\
		/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0g, dcba) */"\n\t					shlq	$3,%%rsi					\n\t"\
		"addq	$0x40	,%%rdi 	/* r0g */			\n\t						addq	%%rsi,%%rax		/* a+p24] */\n\t"\
		"vmovaps	     (%%rdi),%%ymm4				\n\t						addq	%%rsi,%%rbx					\n\t"\
		"vmovaps	0x240(%%rdi),%%ymm6				\n\t						addq	%%rsi,%%rcx					\n\t"\
		"vmovaps	0x020(%%rdi),%%ymm5				\n\t						addq	%%rsi,%%rdx					\n\t"\
		"vmovaps	0x260(%%rdi),%%ymm7				\n\t						vmovaps	%%ymm8 ,    (%%rax)			\n\t"\
		"vmovaps	0x480(%%rdi),%%ymm0				\n\t						vmovaps	%%ymm12,    (%%rcx)			\n\t"\
		"vmovaps	0x6c0(%%rdi),%%ymm2				\n\t						vmovaps	%%ymm9 ,0x20(%%rax)			\n\t"\
		"vmovaps	0x4a0(%%rdi),%%ymm1				\n\t						vmovaps	%%ymm13,0x20(%%rdx)			\n\t"\
		"vmovaps	0x6e0(%%rdi),%%ymm3				\n\t						vaddpd	%%ymm10,%%ymm10,%%ymm10		\n\t"\
		"vsubpd	%%ymm0,%%ymm4,%%ymm4				\n\t						vaddpd	%%ymm15,%%ymm15,%%ymm15		\n\t"\
		"vsubpd	%%ymm2,%%ymm6,%%ymm6				\n\t						vaddpd	%%ymm11,%%ymm11,%%ymm11		\n\t"\
		"vsubpd	%%ymm1,%%ymm5,%%ymm5				\n\t						vaddpd	%%ymm14,%%ymm14,%%ymm14		\n\t"\
		"vsubpd	%%ymm3,%%ymm7,%%ymm7				\n\t						vaddpd	%%ymm8 ,%%ymm10,%%ymm10		\n\t"\
		"vaddpd	%%ymm0,%%ymm0,%%ymm0				\n\t						vaddpd	%%ymm12,%%ymm15,%%ymm15		\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2				\n\t						vaddpd	%%ymm9 ,%%ymm11,%%ymm11		\n\t"\
		"vaddpd	%%ymm1,%%ymm1,%%ymm1				\n\t						vaddpd	%%ymm13,%%ymm14,%%ymm14		\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3				\n\t						vmovaps	%%ymm10,    (%%rbx)			\n\t"\
		"vaddpd	%%ymm4,%%ymm0,%%ymm0				\n\t						vmovaps	%%ymm15,    (%%rdx)			\n\t"\
		"vaddpd	%%ymm6,%%ymm2,%%ymm2				\n\t						vmovaps	%%ymm11,0x20(%%rbx)			\n\t"\
		"vaddpd	%%ymm5,%%ymm1,%%ymm1				\n\t						vmovaps	%%ymm14,0x20(%%rcx)			\n\t"\
		"vaddpd	%%ymm7,%%ymm3,%%ymm3				\n\t"\
		"/* Swap output addresses a <-> d ,b <-> c: */\n\t"\
		"vsubpd	%%ymm2,%%ymm0,%%ymm0				\n\t"\
		"vsubpd	%%ymm7,%%ymm4,%%ymm4				\n\t"\
		"vsubpd	%%ymm3,%%ymm1,%%ymm1				\n\t"\
		"vsubpd	%%ymm6,%%ymm5,%%ymm5				\n\t"\
		/* add,2,1,03 = a+p12+p[0,1,2,3] */\
		"movslq	%[__p12],%%rsi						\n\t"\
		"shlq	$3,%%rsi							\n\t"\
		"subq	%%rsi,%%rax		/* a+p12] */		\n\t"\
		"subq	%%rsi,%%rbx							\n\t"\
		"subq	%%rsi,%%rcx							\n\t"\
		"subq	%%rsi,%%rdx							\n\t"\
		"vmovaps	%%ymm0,    (%%rcx)				\n\t"\
		"vmovaps	%%ymm4,    (%%rbx)				\n\t"\
		"vmovaps	%%ymm1,0x20(%%rcx)				\n\t"\
		"vmovaps	%%ymm5,0x20(%%rax)				\n\t"\
		"vaddpd	%%ymm2,%%ymm2,%%ymm2				\n\t"\
		"vaddpd	%%ymm7,%%ymm7,%%ymm7				\n\t"\
		"vaddpd	%%ymm3,%%ymm3,%%ymm3				\n\t"\
		"vaddpd	%%ymm6,%%ymm6,%%ymm6				\n\t"\
		"vaddpd	%%ymm0,%%ymm2,%%ymm2				\n\t"\
		"vaddpd	%%ymm4,%%ymm7,%%ymm7				\n\t"\
		"vaddpd	%%ymm1,%%ymm3,%%ymm3				\n\t"\
		"vaddpd	%%ymm5,%%ymm6,%%ymm6				\n\t"\
		"vmovaps	%%ymm2,    (%%rdx)				\n\t"\
		"vmovaps	%%ymm7,    (%%rax)				\n\t"\
		"vmovaps	%%ymm3,0x20(%%rdx)				\n\t"\
		"vmovaps	%%ymm6,0x20(%%rbx)				\n\t"\
		"\n\t"\
			:					/* outputs: none */\
			: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
			 ,[__p01] "m" (Xp01)\
			 ,[__p02] "m" (Xp02)\
			 ,[__p03] "m" (Xp03)\
			 ,[__p04] "m" (Xp04)\
			 ,[__p08] "m" (Xp08)\
			 ,[__p12] "m" (Xp12)\
			 ,[__p16] "m" (Xp16)\
			 ,[__p20] "m" (Xp20)\
			 ,[__p24] "m" (Xp24)\
			 ,[__p28] "m" (Xp28)\
			 ,[__p32] "m" (Xp32)\
			 ,[__r00] "m" (Xr00)\
			 ,[__in0] "m" (Xin0)\
			 ,[__cc1] "m" (Xcc1)\
			: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
	}

 #else

  #define USE_64BIT_ASM_STYLE	1	// My x86 timings indicate fancier versions below using ymm0-15 for the radix-16 DFT is faster.

  #if USE_64BIT_ASM_STYLE	// True: Deeper 64-bit-ified version of the original 32-bit ASM macros, using all of xmm0-15

	#define	SSE2_RADIX36_DIT_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xp28,Xp32,Xr00,Xout,Xcc1)\
	{\
	__asm__ volatile (\
		/*	add0,1,3,2 = a+p00+p[0,1,2,3]: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x120, 0x240, r00: */\
			"movq	%[__r00],%%rsi\n\t"\
			"movq	%[__add],%%rax\n\t"\
			"movslq	%[__p01],%%rbx\n\t"\
			"movslq	%[__p02],%%rcx\n\t"\
			"movslq	%[__p03],%%rdx\n\t"\
			"shlq	$3,%%rbx		/* Pointer offset for floating doubles */\n\t"\
			"shlq	$3,%%rcx\n\t"\
			"shlq	$3,%%rdx\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* ecx <-> edx */\n\t"\
			"movaps	     (%%rax),%%xmm2	\n\t"\
			"movaps	     (%%rdx),%%xmm6	\n\t"\
			"movaps	0x010(%%rax),%%xmm3	\n\t"\
			"movaps	0x010(%%rdx),%%xmm7	\n\t"\
			"movaps	     (%%rbx),%%xmm0	\n\t"\
			"movaps	     (%%rcx),%%xmm4	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t"\
			"movaps	0x010(%%rcx),%%xmm5	\n\t"\
			"subpd	%%xmm0,%%xmm2		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t"\
			"subpd	%%xmm1,%%xmm3		\n\t"\
			"subpd	%%xmm5,%%xmm7		\n\t"\
			"addpd	%%xmm0,%%xmm0		\n\t"				/*	add2,3,0,1 = a+p04+p[0,1,2,3]: cdab, r04: */\
			"addpd	%%xmm4,%%xmm4		\n\t					movslq	%[__p04],%%rdi\n\t"\
			"addpd	%%xmm1,%%xmm1		\n\t					shlq	$3,%%rdi \n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t					addq	%%rdi,%%rax	/* add04 */\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t					addq	%%rdi,%%rbx\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t					addq	%%rdi,%%rcx\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t					addq	%%rdi,%%rdx\n\t"\
			"addpd	%%xmm7,%%xmm5		\n\t					/* e[ab]x <-> e[cd]x */\n\t"\
			"subpd	%%xmm4,%%xmm0		\n\t					movaps	     (%%rcx),%%xmm10	\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t					movaps	     (%%rax),%%xmm14	\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t					movaps	0x010(%%rcx),%%xmm11	\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t					movaps	0x010(%%rax),%%xmm15	\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t					movaps	     (%%rdx),%%xmm8 	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t					movaps	     (%%rbx),%%xmm12	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t					movaps	0x010(%%rdx),%%xmm9 	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t					movaps	0x010(%%rbx),%%xmm13	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t					subpd	%%xmm8 ,%%xmm10		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t					subpd	%%xmm12,%%xmm14		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t					subpd	%%xmm9 ,%%xmm11		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t					subpd	%%xmm13,%%xmm15		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t					addpd	%%xmm8 ,%%xmm8 		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t					addpd	%%xmm12,%%xmm12		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t					addpd	%%xmm9 ,%%xmm9 		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t					addpd	%%xmm13,%%xmm13		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t					addpd	%%xmm10,%%xmm8 		\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t					addpd	%%xmm14,%%xmm12		\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t					addpd	%%xmm11,%%xmm9 		\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t					addpd	%%xmm15,%%xmm13		\n\t"\
		/*	add1,0,2,3 = a+p08+p[0,1,2,3]: bacd, r02: */"		addq	$0x040,%%rsi	/* r04 */\n\t"\
			"addq	%%rdi,%%rax	/* add08 */\n\t					subpd	%%xmm12,%%xmm8 		\n\t"\
			"addq	%%rdi,%%rbx\n\t								subpd	%%xmm15,%%xmm10		\n\t"\
			"addq	%%rdi,%%rcx\n\t								subpd	%%xmm13,%%xmm9 		\n\t"\
			"addq	%%rdi,%%rdx\n\t								subpd	%%xmm14,%%xmm11		\n\t"\
			"/* eax <-> ebx */\n\t								movaps	%%xmm8 ,0x240(%%rsi)	\n\t"\
			"movaps	     (%%rbx),%%xmm2	\n\t					movaps	%%xmm10,0x360(%%rsi)	\n\t"\
			"movaps	     (%%rcx),%%xmm6	\n\t					movaps	%%xmm9 ,0x250(%%rsi)	\n\t"\
			"movaps	0x010(%%rbx),%%xmm3	\n\t					movaps	%%xmm11,0x130(%%rsi)	\n\t"\
			"movaps	0x010(%%rcx),%%xmm7	\n\t					addpd	%%xmm12,%%xmm12		\n\t"\
			"movaps	     (%%rax),%%xmm0	\n\t					addpd	%%xmm15,%%xmm15		\n\t"\
			"movaps	     (%%rdx),%%xmm4	\n\t					addpd	%%xmm13,%%xmm13		\n\t"\
			"movaps	0x010(%%rax),%%xmm1	\n\t					addpd	%%xmm14,%%xmm14		\n\t"\
			"movaps	0x010(%%rdx),%%xmm5	\n\t					addpd	%%xmm8 ,%%xmm12		\n\t"\
			"subpd	%%xmm0,%%xmm2		\n\t					addpd	%%xmm10,%%xmm15		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t					addpd	%%xmm9 ,%%xmm13		\n\t"\
			"subpd	%%xmm1,%%xmm3		\n\t					addpd	%%xmm11,%%xmm14		\n\t"\
			"subpd	%%xmm5,%%xmm7		\n\t					movaps	%%xmm12,     (%%rsi)	\n\t"\
			"addpd	%%xmm0,%%xmm0		\n\t					movaps	%%xmm15,0x120(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t					movaps	%%xmm13,0x010(%%rsi)	\n\t"\
			"addpd	%%xmm1,%%xmm1		\n\t					movaps	%%xmm14,0x370(%%rsi)	\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"				/*	add3,2,1,0 = a+p12+p[0,1,2,3]: dcba, r0c: */\
			"addpd	%%xmm2,%%xmm0		\n\t					addq	%%rdi,%%rax	/* add12 */\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t					addq	%%rdi,%%rbx\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t					addq	%%rdi,%%rcx\n\t"\
			"addpd	%%xmm7,%%xmm5		\n\t					addq	%%rdi,%%rdx\n\t"\
			"subq	$0x020,%%rsi	/* r02 */\n\t				/* eax <-> edx, ebx <-> ecx */\n\t"\
			"subpd	%%xmm4,%%xmm0		\n\t					movaps	     (%%rdx),%%xmm10	\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t					movaps	     (%%rbx),%%xmm14	\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t					movaps	0x010(%%rdx),%%xmm11	\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t					movaps	0x010(%%rbx),%%xmm15	\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t					movaps	     (%%rcx),%%xmm8 	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t					movaps	     (%%rax),%%xmm12	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t					movaps	0x010(%%rcx),%%xmm9 	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t					movaps	0x010(%%rax),%%xmm13	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t					subpd	%%xmm8 ,%%xmm10		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t					subpd	%%xmm12,%%xmm14		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t					subpd	%%xmm9 ,%%xmm11		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t					subpd	%%xmm13,%%xmm15		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t					addpd	%%xmm8 ,%%xmm8 		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t					addpd	%%xmm12,%%xmm12		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t					addpd	%%xmm9 ,%%xmm9 		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t					addpd	%%xmm13,%%xmm13		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t					addpd	%%xmm10,%%xmm8 		\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t					addpd	%%xmm14,%%xmm12		\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t					addpd	%%xmm11,%%xmm9 		\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t					addpd	%%xmm15,%%xmm13		\n\t"\
		/*	add0,1,3,2 = a+p16+p[0,1,2,3]: abdc, r0g: */"		addq	$0x0a0,%%rsi	/* r0c */\n\t"\
			"addq	%%rdi,%%rax	/* add16 */\n\t					subpd	%%xmm12,%%xmm8 		\n\t"\
			"addq	%%rdi,%%rbx\n\t								subpd	%%xmm15,%%xmm10		\n\t"\
			"addq	%%rdi,%%rcx\n\t								subpd	%%xmm13,%%xmm9 		\n\t"\
			"addq	%%rdi,%%rdx\n\t								subpd	%%xmm14,%%xmm11		\n\t"\
			"/* ecx <-> edx */\n\t								movaps	%%xmm8 ,0x240(%%rsi)	\n\t"\
			"movaps	     (%%rax),%%xmm2	\n\t					movaps	%%xmm10,0x360(%%rsi)	\n\t"\
			"movaps	     (%%rdx),%%xmm6	\n\t					movaps	%%xmm9 ,0x250(%%rsi)	\n\t"\
			"movaps	0x010(%%rax),%%xmm3	\n\t					movaps	%%xmm11,0x130(%%rsi)	\n\t"\
			"movaps	0x010(%%rdx),%%xmm7	\n\t					addpd	%%xmm12,%%xmm12		\n\t"\
			"movaps	     (%%rbx),%%xmm0	\n\t					addpd	%%xmm15,%%xmm15		\n\t"\
			"movaps	     (%%rcx),%%xmm4	\n\t					addpd	%%xmm13,%%xmm13		\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t					addpd	%%xmm14,%%xmm14		\n\t"\
			"movaps	0x010(%%rcx),%%xmm5	\n\t					addpd	%%xmm8 ,%%xmm12		\n\t"\
			"subpd	%%xmm0,%%xmm2		\n\t					addpd	%%xmm10,%%xmm15		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t					addpd	%%xmm9 ,%%xmm13		\n\t"\
			"subpd	%%xmm1,%%xmm3		\n\t					addpd	%%xmm11,%%xmm14		\n\t"\
			"subpd	%%xmm5,%%xmm7		\n\t					movaps	%%xmm12,     (%%rsi)	\n\t"\
			"addpd	%%xmm0,%%xmm0		\n\t					movaps	%%xmm15,0x120(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t					movaps	%%xmm13,0x010(%%rsi)	\n\t"\
			"addpd	%%xmm1,%%xmm1		\n\t					movaps	%%xmm14,0x370(%%rsi)	\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"				/*	add2,3,0,1 = a+p20+p[0,1,2,3]: cdab, r0e: */\
			"addpd	%%xmm2,%%xmm0		\n\t					addq	%%rdi,%%rax	/* add20 */\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t					addq	%%rdi,%%rbx\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t					addq	%%rdi,%%rcx\n\t"\
			"addpd	%%xmm7,%%xmm5		\n\t					addq	%%rdi,%%rdx\n\t"\
			"addq	$0x040,%%rsi	/* r0g */\n\t				/* e[ab]x <-> e[cd]x */\n\t"\
			"subpd	%%xmm4,%%xmm0		\n\t					movaps	     (%%rcx),%%xmm10	\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t					movaps	     (%%rax),%%xmm14	\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t					movaps	0x010(%%rcx),%%xmm11	\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t					movaps	0x010(%%rax),%%xmm15	\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t					movaps	     (%%rdx),%%xmm8 	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t					movaps	     (%%rbx),%%xmm12	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t					movaps	0x010(%%rdx),%%xmm9 	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t					movaps	0x010(%%rbx),%%xmm13	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t					subpd	%%xmm8 ,%%xmm10		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t					subpd	%%xmm12,%%xmm14		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t					subpd	%%xmm9 ,%%xmm11		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t					subpd	%%xmm13,%%xmm15		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t					addpd	%%xmm8 ,%%xmm8 		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t					addpd	%%xmm12,%%xmm12		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t					addpd	%%xmm9 ,%%xmm9 		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t					addpd	%%xmm13,%%xmm13		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t					addpd	%%xmm10,%%xmm8 		\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t					addpd	%%xmm14,%%xmm12		\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t					addpd	%%xmm11,%%xmm9 		\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t					addpd	%%xmm15,%%xmm13		\n\t"\
		/*	add1,0,2,3 = a+p24+p[0,1,2,3]: bacd, r08: */"		subq	$0x020,%%rsi	/* r0e */\n\t"\
			"addq	%%rdi,%%rax	/* add24 */\n\t					subpd	%%xmm12,%%xmm8 		\n\t"\
			"addq	%%rdi,%%rbx\n\t								subpd	%%xmm15,%%xmm10		\n\t"\
			"addq	%%rdi,%%rcx\n\t								subpd	%%xmm13,%%xmm9 		\n\t"\
			"addq	%%rdi,%%rdx\n\t								subpd	%%xmm14,%%xmm11		\n\t"\
			"/* eax <-> ebx */\n\t								movaps	%%xmm8 ,0x240(%%rsi)	\n\t"\
			"movaps	     (%%rbx),%%xmm2	\n\t					movaps	%%xmm10,0x360(%%rsi)	\n\t"\
			"movaps	     (%%rcx),%%xmm6	\n\t					movaps	%%xmm9 ,0x250(%%rsi)	\n\t"\
			"movaps	0x010(%%rbx),%%xmm3	\n\t					movaps	%%xmm11,0x130(%%rsi)	\n\t"\
			"movaps	0x010(%%rcx),%%xmm7	\n\t					addpd	%%xmm12,%%xmm12		\n\t"\
			"movaps	     (%%rax),%%xmm0	\n\t					addpd	%%xmm15,%%xmm15		\n\t"\
			"movaps	     (%%rdx),%%xmm4	\n\t					addpd	%%xmm13,%%xmm13		\n\t"\
			"movaps	0x010(%%rax),%%xmm1	\n\t					addpd	%%xmm14,%%xmm14		\n\t"\
			"movaps	0x010(%%rdx),%%xmm5	\n\t					addpd	%%xmm8 ,%%xmm12		\n\t"\
			"subpd	%%xmm0,%%xmm2		\n\t					addpd	%%xmm10,%%xmm15		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t					addpd	%%xmm9 ,%%xmm13		\n\t"\
			"subpd	%%xmm1,%%xmm3		\n\t					addpd	%%xmm11,%%xmm14		\n\t"\
			"subpd	%%xmm5,%%xmm7		\n\t					movaps	%%xmm12,     (%%rsi)	\n\t"\
			"addpd	%%xmm0,%%xmm0		\n\t					movaps	%%xmm15,0x120(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t					movaps	%%xmm13,0x010(%%rsi)	\n\t"\
			"addpd	%%xmm1,%%xmm1		\n\t					movaps	%%xmm14,0x370(%%rsi)	\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"				/*	add3,2,1,0 = a+p28+p[0,1,2,3]: dcba, r06: */\
			"addpd	%%xmm2,%%xmm0		\n\t					addq	%%rdi,%%rax	/* add28 */\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t					addq	%%rdi,%%rbx\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t					addq	%%rdi,%%rcx\n\t"\
			"addpd	%%xmm7,%%xmm5		\n\t					addq	%%rdi,%%rdx\n\t"\
			"subq	$0x060,%%rsi	/* r08 */\n\t				/* eax <-> edx, ebx <-> ecx */\n\t"\
			"subpd	%%xmm4,%%xmm0		\n\t					movaps	     (%%rdx),%%xmm10	\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t					movaps	     (%%rbx),%%xmm14	\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t					movaps	0x010(%%rdx),%%xmm11	\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t					movaps	0x010(%%rbx),%%xmm15	\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t					movaps	     (%%rcx),%%xmm8 	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t					movaps	     (%%rax),%%xmm12	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t					movaps	0x010(%%rcx),%%xmm9 	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t					movaps	0x010(%%rax),%%xmm13	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t					subpd	%%xmm8 ,%%xmm10		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t					subpd	%%xmm12,%%xmm14		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t					subpd	%%xmm9 ,%%xmm11		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t					subpd	%%xmm13,%%xmm15		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t					addpd	%%xmm8 ,%%xmm8 		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t					addpd	%%xmm12,%%xmm12		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t					addpd	%%xmm9 ,%%xmm9 		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t					addpd	%%xmm13,%%xmm13		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t					addpd	%%xmm10,%%xmm8 		\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t					addpd	%%xmm14,%%xmm12		\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t					addpd	%%xmm11,%%xmm9 		\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t					addpd	%%xmm15,%%xmm13		\n\t"\
		/*	add0,1,3,2 = a+p32+p[0,1,2,3]: abdc, r0a: */"		subq	$0x020,%%rsi	/* r04 */\n\t"\
			"addq	%%rdi,%%rax	/* add32 */\n\t					subpd	%%xmm12,%%xmm8 		\n\t"\
			"addq	%%rdi,%%rbx\n\t								subpd	%%xmm15,%%xmm10		\n\t"\
			"addq	%%rdi,%%rcx\n\t								subpd	%%xmm13,%%xmm9 		\n\t"\
			"addq	%%rdi,%%rdx\n\t								subpd	%%xmm14,%%xmm11		\n\t"\
			"/* ecx <-> edx */\n\t								movaps	%%xmm8 ,0x240(%%rsi)	\n\t"\
			"movaps	     (%%rax),%%xmm2	\n\t					movaps	%%xmm10,0x360(%%rsi)	\n\t"\
			"movaps	     (%%rdx),%%xmm6	\n\t					movaps	%%xmm9 ,0x250(%%rsi)	\n\t"\
			"movaps	0x010(%%rax),%%xmm3	\n\t					movaps	%%xmm11,0x130(%%rsi)	\n\t"\
			"movaps	0x010(%%rdx),%%xmm7	\n\t					addpd	%%xmm12,%%xmm12		\n\t"\
			"movaps	     (%%rbx),%%xmm0	\n\t					addpd	%%xmm15,%%xmm15		\n\t"\
			"movaps	     (%%rcx),%%xmm4	\n\t					addpd	%%xmm13,%%xmm13		\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t					addpd	%%xmm14,%%xmm14		\n\t"\
			"movaps	0x010(%%rcx),%%xmm5	\n\t					addpd	%%xmm8 ,%%xmm12		\n\t"\
			"subpd	%%xmm0,%%xmm2		\n\t					addpd	%%xmm10,%%xmm15		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t					addpd	%%xmm9 ,%%xmm13		\n\t"\
			"subpd	%%xmm1,%%xmm3		\n\t					addpd	%%xmm11,%%xmm14		\n\t"\
			"subpd	%%xmm5,%%xmm7		\n\t					movaps	%%xmm12,     (%%rsi)	\n\t"\
			"addpd	%%xmm0,%%xmm0		\n\t					movaps	%%xmm15,0x120(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t					movaps	%%xmm13,0x010(%%rsi)	\n\t"\
			"addpd	%%xmm1,%%xmm1		\n\t					movaps	%%xmm14,0x370(%%rsi)	\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t"\
			"addpd	%%xmm7,%%xmm5		\n\t"\
			"addq	$0x040,%%rsi	/* r04 */\n\t"\
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
			"\n\t"\
			"/***************************************/\n\t"\
			"/*...and now do 4 radix-9 transforms...*/\n\t"\
			"/***************************************/\n\t"\
			"\n\t"\
			/* SSE2_RADIX_09_DIT_0TWIDDLE(r00,s1p[00,32,28,24,20,16,12,08,04]r): */\
			"movq	%[__out],%%rsi 	/* __o0-8: rsi,di store o-addresses */\n\t"\
			"movq	%%rsi		,%%rdi\n\t"\
			"movq	%[__r00]	,%%rax\n\t"/* rcol locals based on r10 = r00+0x120 */\
		"movq	%[__cc1],%%rdx 	/* rdx stores trig addresses throughout */\n\t"\
			"movq	%%rax		,%%rbx\n\t"						/* SSE2_RADIX_09_DIT_0TWIDDLE(r10,s1p[27,23,19,15,11,07,03,35,31]r): */\
			"movq	%%rax		,%%rcx\n\t						/* rcol: r10,r11 store o-addresses */\n\t"\
			"addq	$0x20		,%%rbx\n\t						leaq	0x360(%%rsi),%%r10 		/* s1p27 */\n\t"\
			"addq	$0x40		,%%rcx\n\t						movq	%%r10		,%%r11\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t						movaps	0x120(%%rbx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t						movaps	0x130(%%rbx)	,%%xmm11\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t						movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t						movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"movaps	    (%%rcx)	,%%xmm6\n\t						movaps	0x120(%%rcx)	,%%xmm14\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm7\n\t						movaps	0x130(%%rcx)	,%%xmm15\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t						movaps	%%xmm10		,%%xmm12\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t						movaps	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t						addpd	%%xmm14		,%%xmm10\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t						addpd	%%xmm15		,%%xmm11\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t						subpd	%%xmm14		,%%xmm12\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t						subpd	%%xmm15		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t						addpd	%%xmm10		,%%xmm8 \n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t						addpd	%%xmm11		,%%xmm9 \n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"	/* These are shared between lcol & rcol */\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t				movaps	%%xmm8 		,0x120(%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t				movaps	%%xmm9 		,0x130(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t						mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t						mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t						mulpd	%%xmm7		,%%xmm12\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t						mulpd	%%xmm7		,%%xmm13\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t						addpd	%%xmm8 		,%%xmm10\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t						addpd	%%xmm9 		,%%xmm11\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t						movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t						movaps	%%xmm11		,%%xmm9 \n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t						addpd	%%xmm13		,%%xmm10\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t						subpd	%%xmm12		,%%xmm11\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t						subpd	%%xmm13		,%%xmm8 \n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t						addpd	%%xmm12		,%%xmm9 \n\t"\
			"movaps	%%xmm2		,    (%%rbx)\n\t				movaps	%%xmm10		,0x120(%%rbx)\n\t"\
			"movaps	%%xmm3		,0x10(%%rbx)\n\t				movaps	%%xmm11		,0x130(%%rbx)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t				movaps	%%xmm8 		,0x120(%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t				movaps	%%xmm9 		,0x130(%%rcx)\n\t"\
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%rbx\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t						movaps	0x120(%%rbx)	,%%xmm12\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t						movaps	0x130(%%rbx)	,%%xmm13\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t						movaps	0x120(%%rcx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t						movaps	0x130(%%rcx)	,%%xmm11\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t						movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t						movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t						subpd	%%xmm10		,%%xmm12\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t						subpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t						addpd	%%xmm10		,%%xmm10\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t						addpd	%%xmm11		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t						addpd	%%xmm12		,%%xmm10\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t						addpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t						addpd	%%xmm10		,%%xmm8 \n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t						addpd	%%xmm11		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t				movaps	%%xmm8 		,0x120(%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t				movaps	%%xmm9 		,0x130(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t						mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t						mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t						mulpd	%%xmm7		,%%xmm12\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t						mulpd	%%xmm7		,%%xmm13\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t						addpd	%%xmm8 		,%%xmm10\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t						addpd	%%xmm9 		,%%xmm11\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t						movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t						movaps	%%xmm11		,%%xmm9 \n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t						addpd	%%xmm13		,%%xmm10\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t						subpd	%%xmm12		,%%xmm11\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t						subpd	%%xmm13		,%%xmm8 \n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t						addpd	%%xmm12		,%%xmm9 \n\t"\
			"movaps	%%xmm2		,    (%%rbx)\n\t				movaps	%%xmm10		,0x120(%%rbx)\n\t"\
			"movaps	%%xmm3		,0x10(%%rbx)\n\t				movaps	%%xmm11		,0x130(%%rbx)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t				movaps	%%xmm8 		,0x120(%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t				movaps	%%xmm9 		,0x130(%%rcx)\n\t"\
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%rbx\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t						movaps	0x120(%%rbx)	,%%xmm12\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t						movaps	0x130(%%rbx)	,%%xmm13\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t						movaps	0x120(%%rcx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t						movaps	0x130(%%rcx)	,%%xmm11\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t						movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t						movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t						subpd	%%xmm10		,%%xmm12\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t						subpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t						addpd	%%xmm10		,%%xmm10\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t						addpd	%%xmm11		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t						addpd	%%xmm12		,%%xmm10\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t						addpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t						addpd	%%xmm10		,%%xmm8 \n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t						addpd	%%xmm11		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t				movaps	%%xmm8 		,0x120(%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t				movaps	%%xmm9 		,0x130(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t						mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t						mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t						mulpd	%%xmm7		,%%xmm12\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t						mulpd	%%xmm7		,%%xmm13\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t						addpd	%%xmm8 		,%%xmm10\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t						addpd	%%xmm9 		,%%xmm11\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t						movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t						movaps	%%xmm11		,%%xmm9 \n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t						addpd	%%xmm13		,%%xmm10\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t						subpd	%%xmm12		,%%xmm11\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t						subpd	%%xmm13		,%%xmm8 \n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t						addpd	%%xmm12		,%%xmm9 \n\t"\
			"movaps	%%xmm2		,    (%%rbx)\n\t				movaps	%%xmm10		,0x120(%%rbx)\n\t"\
			"movaps	%%xmm3		,0x10(%%rbx)\n\t				movaps	%%xmm11		,0x130(%%rbx)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t				movaps	%%xmm8 		,0x120(%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t				movaps	%%xmm9 		,0x130(%%rcx)\n\t"\
			"/*****************************************/\n\t"\
			"/* now do three more radix-3 transforms: */\n\t"\
			"/*****************************************/\n\t"\
			"subq	$0xc0		,%%rax 		/* __r00 */	\n\t"\
			"subq	$0x80,%%rbx /* ebx stores r02+0xc0 = r00+0xe0; need r06 = r00+0x60 =bx-0x80 */\n\t"\
			"subq	$0x40,%%rcx /* ecx stores r04+0xc0           ; need r0c = r00+0xc0 =bx-0x40 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t						movaps	0x120(%%rbx)	,%%xmm12\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t						movaps	0x130(%%rbx)	,%%xmm13\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t						movaps	0x120(%%rcx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t						movaps	0x130(%%rcx)	,%%xmm11\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t						movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t						movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t						subpd	%%xmm10		,%%xmm12\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t						subpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t						addpd	%%xmm10		,%%xmm10\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t						addpd	%%xmm11		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t						addpd	%%xmm12		,%%xmm10\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t						addpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t						addpd	%%xmm10		,%%xmm8 \n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t						addpd	%%xmm11		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t				movaps	%%xmm8 		,    (%%r10)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t				movaps	%%xmm9 		,0x10(%%r10)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t						mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t						mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t						mulpd	%%xmm7		,%%xmm12\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t						mulpd	%%xmm7		,%%xmm13\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t						addpd	%%xmm8 		,%%xmm10\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t						addpd	%%xmm9 		,%%xmm11\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t						movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t						movaps	%%xmm11		,%%xmm9 \n\t"\
			"addq	$0x300		,%%rdi 	/* __o3 = s1p24r */\n\t	subq	$0x180		,%%r11 		/* __o3 = s1p15r */\n\t"\
			"addq	$0x180		,%%rsi 	/* __o6 = s1p12r */\n\t	subq	$0x300		,%%r10 		/* __o6 = s1p03r */\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t						addpd	%%xmm13		,%%xmm10\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t						subpd	%%xmm12		,%%xmm11\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t						subpd	%%xmm13		,%%xmm8 \n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t						addpd	%%xmm12		,%%xmm9 \n\t"\
			"movaps	%%xmm2		,    (%%rdi)\n\t				movaps	%%xmm10		,    (%%r11)\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)\n\t				movaps	%%xmm11		,0x10(%%r11)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t				movaps	%%xmm8 		,    (%%r10)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t				movaps	%%xmm9 		,0x10(%%r10)\n\t"\
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%rbx\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t						movaps	0x120(%%rbx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t						movaps	0x130(%%rbx)	,%%xmm11\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t						movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t						movaps	%%xmm11		,%%xmm9 \n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t						mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t						mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t						mulpd	%%xmm7		,%%xmm8 \n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t						mulpd	%%xmm7		,%%xmm9 \n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t						addpd	%%xmm9 		,%%xmm10\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t						subpd	%%xmm8 		,%%xmm11\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t						movaps	0x120(%%rcx)	,%%xmm12\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t						movaps	0x130(%%rcx)	,%%xmm13\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t						movaps	%%xmm12		,%%xmm8 \n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t						movaps	%%xmm13		,%%xmm9 \n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t						mulpd	%%xmm6		,%%xmm12\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t						mulpd	%%xmm6		,%%xmm13\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t						mulpd	%%xmm7		,%%xmm8 \n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t						mulpd	%%xmm7		,%%xmm9 \n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t						addpd	%%xmm9 		,%%xmm12\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t						subpd	%%xmm8 		,%%xmm13\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t						movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t						movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"subq	$0x080		,%%rsi 	/* __o7 = s1p08r */\n\t	addq	$0x400		,%%r10 		/* __o7 = s1p35r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t						subpd	%%xmm12		,%%xmm10\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t						subpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t						addpd	%%xmm12		,%%xmm12\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t						addpd	%%xmm13		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t						addpd	%%xmm10		,%%xmm12\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t						addpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t						addpd	%%xmm12		,%%xmm8 \n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t						addpd	%%xmm13		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t				movaps	%%xmm8 		,    (%%r10)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t				movaps	%%xmm9 		,0x10(%%r10)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t						mulpd	%%xmm6		,%%xmm12\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t						mulpd	%%xmm6		,%%xmm13\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t						mulpd	%%xmm7		,%%xmm10\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t						mulpd	%%xmm7		,%%xmm11\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t						addpd	%%xmm8 		,%%xmm12\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t						addpd	%%xmm9 		,%%xmm13\n\t"\
			"addq	$0x100		,%%rdi 	/* __o1 = s1p32r */\n\t	addq	$0x100		,%%r11 		/* __o1 = s1p23r */\n\t"\
			"addq	$0x180		,%%rsi 	/* __o4 = s1p20r */\n\t	subq	$0x300		,%%r10 		/* __o4 = s1p11r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t						movaps	%%xmm12		,%%xmm8 \n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t						movaps	%%xmm13		,%%xmm9 \n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t						addpd	%%xmm11		,%%xmm12\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t						subpd	%%xmm10		,%%xmm13\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t						subpd	%%xmm11		,%%xmm8 \n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t						addpd	%%xmm10		,%%xmm9 \n\t"\
			"movaps	%%xmm4		,    (%%rdi)\n\t				movaps	%%xmm12		,    (%%r11)\n\t"\
			"movaps	%%xmm5		,0x10(%%rdi)\n\t				movaps	%%xmm13		,0x10(%%r11)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t				movaps	%%xmm8 		,    (%%r10)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t				movaps	%%xmm9 		,0x10(%%r10)\n\t"\
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%rbx\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t						movaps	0x120(%%rbx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t						movaps	0x130(%%rbx)	,%%xmm11\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t	"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t	"\
			"movaps	%%xmm2		,%%xmm0\n\t						movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t						movaps	%%xmm11		,%%xmm9 \n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t						mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t						mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t						mulpd	%%xmm7		,%%xmm8 \n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t						mulpd	%%xmm7		,%%xmm9 \n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t						addpd	%%xmm9 		,%%xmm10\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t						subpd	%%xmm8 		,%%xmm11\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t						movaps	0x120(%%rcx)	,%%xmm12\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t						movaps	0x130(%%rcx)	,%%xmm13\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t						movaps	%%xmm12		,%%xmm8 \n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t						movaps	%%xmm13		,%%xmm9 \n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t						mulpd	%%xmm6		,%%xmm12\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t						mulpd	%%xmm6		,%%xmm13\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t						mulpd	%%xmm7		,%%xmm8 \n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t						mulpd	%%xmm7		,%%xmm9 \n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t						addpd	%%xmm9 		,%%xmm12\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t						subpd	%%xmm8 		,%%xmm13\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t						movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t						movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"subq	$0x080		,%%rsi 	/* __o5 = s1p16r */\n\t	subq	$0x080		,%%r10 		/* __o5 = s1p07r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t						subpd	%%xmm12		,%%xmm10\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t						subpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t						addpd	%%xmm12		,%%xmm12\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t						addpd	%%xmm13		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t						addpd	%%xmm10		,%%xmm12\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t						addpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t						addpd	%%xmm12		,%%xmm8 \n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t						addpd	%%xmm13		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t				movaps	%%xmm8 		,    (%%r10)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t				movaps	%%xmm9 		,0x10(%%r10)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t						mulpd	%%xmm6		,%%xmm12\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t						mulpd	%%xmm6		,%%xmm13\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t						mulpd	%%xmm7		,%%xmm10\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t						mulpd	%%xmm7		,%%xmm11\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t						addpd	%%xmm8 		,%%xmm12\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t						addpd	%%xmm9 		,%%xmm13\n\t"\
			"subq	$0x380		,%%rdi 	/* __o8 = s1p04r */\n\t	addq	$0x100		,%%r11 		/* __o8 = s1p31r */\n\t"\
			"addq	$0x180		,%%rsi 	/* __o2 = s1p28r */\n\t	addq	$0x180		,%%r10 		/* __o2 = s1p19r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t						movaps	%%xmm12		,%%xmm8 \n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t						movaps	%%xmm13		,%%xmm9 \n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t						addpd	%%xmm11		,%%xmm12\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t						subpd	%%xmm10		,%%xmm13\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t						subpd	%%xmm11		,%%xmm8 \n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t						addpd	%%xmm10		,%%xmm9 \n\t"\
			"movaps	%%xmm4		,    (%%rdi)\n\t				movaps	%%xmm12		,    (%%r11)\n\t"\
			"movaps	%%xmm5		,0x10(%%rdi)\n\t				movaps	%%xmm13		,0x10(%%r11)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t				movaps	%%xmm8 		,    (%%r10)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t				movaps	%%xmm9 		,0x10(%%r10)\n\t"\
			/* SSE2_RADIX_09_DIT_0TWIDDLE(r20,s1p[18,14,10,06,02,34,30,26,22]r): */\
			"movq	%[__out]	,%%rsi\n\t						movq	%[__r00]	,%%rax\n\t"/* rcol locals based on r10 = r00+0x120 */\
			"addq	$0x240		,%%rsi 	/* s1p18 */\n\t			addq	$0x240		,%%rax\n\t	/* __r20 */\n\t"/* rcol locals based on r30 = r20+0x120 */\
			"movq	%%rsi		,%%rdi\n\t						movq	%[__cc1]	,%%rdx	\n\t"\
			"movq	%%rax		,%%rbx\n\t"\
			"movq	%%rax		,%%rcx\n\t"						/* SSE2_RADIX_09_DIT_0TWIDDLE(r30,s1p[09,05,01,33,29,25,21,17,13]r): */\
			"addq	$0x20		,%%rbx\n\t						/* rcol: r10,r11 store o-addresses */\n\t"\
			"addq	$0x40		,%%rcx\n\t						leaq	-0x120(%%rsi),%%r10 		/* s1p09 */\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t			movq	%%r10		,%%r11\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t						movaps	0x120(%%rbx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t						movaps	0x130(%%rbx)	,%%xmm11\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t						movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t						movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"movaps	    (%%rcx)	,%%xmm6\n\t						movaps	0x120(%%rcx)	,%%xmm14\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm7\n\t						movaps	0x130(%%rcx)	,%%xmm15\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t						movaps	%%xmm10		,%%xmm12\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t						movaps	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t						addpd	%%xmm14		,%%xmm10\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t						addpd	%%xmm15		,%%xmm11\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t						subpd	%%xmm14		,%%xmm12\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t						subpd	%%xmm15		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t						addpd	%%xmm10		,%%xmm8 \n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t						addpd	%%xmm11		,%%xmm9 \n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t				movaps	%%xmm8 		,0x120(%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t				movaps	%%xmm9 		,0x130(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t						mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t						mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t						mulpd	%%xmm7		,%%xmm12\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t						mulpd	%%xmm7		,%%xmm13\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t						addpd	%%xmm8 		,%%xmm10\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t						addpd	%%xmm9 		,%%xmm11\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t						movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t						movaps	%%xmm11		,%%xmm9 \n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t						addpd	%%xmm13		,%%xmm10\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t						subpd	%%xmm12		,%%xmm11\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t						subpd	%%xmm13		,%%xmm8 \n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t						addpd	%%xmm12		,%%xmm9 \n\t"\
			"movaps	%%xmm2		,    (%%rbx)\n\t				movaps	%%xmm10		,0x120(%%rbx)\n\t"\
			"movaps	%%xmm3		,0x10(%%rbx)\n\t				movaps	%%xmm11		,0x130(%%rbx)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t				movaps	%%xmm8 		,0x120(%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t				movaps	%%xmm9 		,0x130(%%rcx)\n\t"\
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%rbx\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t						movaps	0x120(%%rbx)	,%%xmm12\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t						movaps	0x130(%%rbx)	,%%xmm13\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t						movaps	0x120(%%rcx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t						movaps	0x130(%%rcx)	,%%xmm11\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t						movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t						movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t						subpd	%%xmm10		,%%xmm12\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t						subpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t						addpd	%%xmm10		,%%xmm10\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t						addpd	%%xmm11		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t						addpd	%%xmm12		,%%xmm10\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t						addpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t						addpd	%%xmm10		,%%xmm8 \n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t						addpd	%%xmm11		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t				movaps	%%xmm8 		,0x120(%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t				movaps	%%xmm9 		,0x130(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t						mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t						mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t						mulpd	%%xmm7		,%%xmm12\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t						mulpd	%%xmm7		,%%xmm13\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t						addpd	%%xmm8 		,%%xmm10\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t						addpd	%%xmm9 		,%%xmm11\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t						movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t						movaps	%%xmm11		,%%xmm9 \n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t						addpd	%%xmm13		,%%xmm10\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t						subpd	%%xmm12		,%%xmm11\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t						subpd	%%xmm13		,%%xmm8 \n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t						addpd	%%xmm12		,%%xmm9 \n\t"\
			"movaps	%%xmm2		,    (%%rbx)\n\t				movaps	%%xmm10		,0x120(%%rbx)\n\t"\
			"movaps	%%xmm3		,0x10(%%rbx)\n\t				movaps	%%xmm11		,0x130(%%rbx)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t				movaps	%%xmm8 		,0x120(%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t				movaps	%%xmm9 		,0x130(%%rcx)\n\t"\
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%rbx\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t						movaps	0x120(%%rbx)	,%%xmm12\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t						movaps	0x130(%%rbx)	,%%xmm13\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t						movaps	0x120(%%rcx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t						movaps	0x130(%%rcx)	,%%xmm11\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t						movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t						movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t						subpd	%%xmm10		,%%xmm12\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t						subpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t						addpd	%%xmm10		,%%xmm10\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t						addpd	%%xmm11		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t						addpd	%%xmm12		,%%xmm10\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t						addpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t						addpd	%%xmm10		,%%xmm8 \n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t						addpd	%%xmm11		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t				movaps	%%xmm8 		,0x120(%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t				movaps	%%xmm9 		,0x130(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t						mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t						mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t						mulpd	%%xmm7		,%%xmm12\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t						mulpd	%%xmm7		,%%xmm13\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t						addpd	%%xmm8 		,%%xmm10\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t						addpd	%%xmm9 		,%%xmm11\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t						movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t						movaps	%%xmm11		,%%xmm9 \n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t						addpd	%%xmm13		,%%xmm10\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t						subpd	%%xmm12		,%%xmm11\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t						subpd	%%xmm13		,%%xmm8 \n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t						addpd	%%xmm12		,%%xmm9 \n\t"\
			"movaps	%%xmm2		,    (%%rbx)\n\t				movaps	%%xmm10		,0x120(%%rbx)\n\t"\
			"movaps	%%xmm3		,0x10(%%rbx)\n\t				movaps	%%xmm11		,0x130(%%rbx)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t				movaps	%%xmm8 		,0x120(%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t				movaps	%%xmm9 		,0x130(%%rcx)\n\t"\
			"/*****************************************/\n\t"\
			"/* now do three more radix-3 transforms: */\n\t"\
			"/*****************************************/\n\t"\
			"subq	$0xc0		,%%rax\n\t"\
			"subq	$0x80		,%%rbx\n\t"\
			"subq	$0x40		,%%rcx\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t						movaps	0x120(%%rbx)	,%%xmm12\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t						movaps	0x130(%%rbx)	,%%xmm13\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t						movaps	0x120(%%rcx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t						movaps	0x130(%%rcx)	,%%xmm11\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t						movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t						movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t						subpd	%%xmm10		,%%xmm12\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t						subpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t						addpd	%%xmm10		,%%xmm10\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t						addpd	%%xmm11		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t						addpd	%%xmm12		,%%xmm10\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t						addpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t						addpd	%%xmm10		,%%xmm8 \n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t						addpd	%%xmm11		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t				movaps	%%xmm8 		,    (%%r10)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t				movaps	%%xmm9 		,0x10(%%r10)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t						mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t						mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t						mulpd	%%xmm7		,%%xmm12\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t						mulpd	%%xmm7		,%%xmm13\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t						addpd	%%xmm8 		,%%xmm10\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t						addpd	%%xmm9 		,%%xmm11\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t						movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t						movaps	%%xmm11		,%%xmm9 \n\t"\
			"subq	$0x180		,%%rdi 	/* __o3 = s1p06r */\n\t	addq	$0x300		,%%r11 		/* __o3 = s1p33r */\n\t"\
			"addq	$0x180		,%%rsi 	/* __o6 = s1p30r */\n\t	addq	$0x180		,%%r10 		/* __o6 = s1p21r */\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t						addpd	%%xmm13		,%%xmm10\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t						subpd	%%xmm12		,%%xmm11\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t						subpd	%%xmm13		,%%xmm8 \n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t						addpd	%%xmm12		,%%xmm9 \n\t"\
			"movaps	%%xmm2		,    (%%rdi)\n\t				movaps	%%xmm10		,    (%%r11)\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)\n\t				movaps	%%xmm11		,0x10(%%r11)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t				movaps	%%xmm8 		,    (%%r10)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t				movaps	%%xmm9 		,0x10(%%r10)\n\t"\
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%rbx\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t						movaps	0x120(%%rbx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t						movaps	0x130(%%rbx)	,%%xmm11\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t						movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t						movaps	%%xmm11		,%%xmm9 \n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t						mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t						mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t						mulpd	%%xmm7		,%%xmm8 \n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t						mulpd	%%xmm7		,%%xmm9 \n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t						addpd	%%xmm9 		,%%xmm10\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t						subpd	%%xmm8 		,%%xmm11\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t						movaps	0x120(%%rcx)	,%%xmm12\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t						movaps	0x130(%%rcx)	,%%xmm13\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t						movaps	%%xmm12		,%%xmm8 \n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t						movaps	%%xmm13		,%%xmm9 \n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t						mulpd	%%xmm6		,%%xmm12\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t						mulpd	%%xmm6		,%%xmm13\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t						mulpd	%%xmm7		,%%xmm8 \n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t						mulpd	%%xmm7		,%%xmm9 \n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t						addpd	%%xmm9 		,%%xmm12\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t						subpd	%%xmm8 		,%%xmm13\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t						movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t						movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"subq	$0x080		,%%rsi 	/* __o7 = s1p26r */\n\t	subq	$0x080		,%%r10 		/* __o7 = s1p17r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t						subpd	%%xmm12		,%%xmm10\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t						subpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t						addpd	%%xmm12		,%%xmm12\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t						addpd	%%xmm13		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t						addpd	%%xmm10		,%%xmm12\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t						addpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t						addpd	%%xmm12		,%%xmm8 \n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t						addpd	%%xmm13		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t				movaps	%%xmm8 		,    (%%r10)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t				movaps	%%xmm9 		,0x10(%%r10)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t						mulpd	%%xmm6		,%%xmm12\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t						mulpd	%%xmm6		,%%xmm13\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t						mulpd	%%xmm7		,%%xmm10\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t						mulpd	%%xmm7		,%%xmm11\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t						addpd	%%xmm8 		,%%xmm12\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t						addpd	%%xmm9 		,%%xmm13\n\t"\
			"addq	$0x100		,%%rdi 	/* __o1 = s1p14r */\n\t	subq	$0x380		,%%r11 		/* __o1 = s1p05r */\n\t"\
			"subq	$0x300		,%%rsi 	/* __o4 = s1p02r */\n\t	addq	$0x180		,%%r10 		/* __o4 = s1p29r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t						movaps	%%xmm12		,%%xmm8 \n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t						movaps	%%xmm13		,%%xmm9 \n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t						addpd	%%xmm11		,%%xmm12\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t						subpd	%%xmm10		,%%xmm13\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t						subpd	%%xmm11		,%%xmm8 \n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t						addpd	%%xmm10		,%%xmm9 \n\t"\
			"movaps	%%xmm4		,    (%%rdi)\n\t				movaps	%%xmm12		,    (%%r11)\n\t"\
			"movaps	%%xmm5		,0x10(%%rdi)\n\t				movaps	%%xmm13		,0x10(%%r11)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t				movaps	%%xmm8 		,    (%%r10)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t				movaps	%%xmm9 		,0x10(%%r10)\n\t"\
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%rbx\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t						movaps	0x120(%%rbx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t						movaps	0x130(%%rbx)	,%%xmm11\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t						movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t						movaps	%%xmm11		,%%xmm9 \n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t						mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t						mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t						mulpd	%%xmm7		,%%xmm8 \n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t						mulpd	%%xmm7		,%%xmm9 \n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t						addpd	%%xmm9 		,%%xmm10\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t						subpd	%%xmm8 		,%%xmm11\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t						movaps	0x120(%%rcx)	,%%xmm12\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t						movaps	0x130(%%rcx)	,%%xmm13\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t						movaps	%%xmm12		,%%xmm8 \n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t						movaps	%%xmm13		,%%xmm9 \n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t						mulpd	%%xmm6		,%%xmm12\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t						mulpd	%%xmm6		,%%xmm13\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t						mulpd	%%xmm7		,%%xmm8 \n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t						mulpd	%%xmm7		,%%xmm9 \n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t						addpd	%%xmm9 		,%%xmm12\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t						subpd	%%xmm8 		,%%xmm13\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t	"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t	"\
			"movaps	    (%%rax)	,%%xmm0\n\t						movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t						movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"addq	$0x400		,%%rsi 	/* __o5 = s1p34r */\n\t	subq	$0x080		,%%r10 		/* __o5 = s1p25r */\n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t						subpd	%%xmm12		,%%xmm10\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t						subpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t						addpd	%%xmm12		,%%xmm12\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t						addpd	%%xmm13		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t						addpd	%%xmm10		,%%xmm12\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t						addpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t						addpd	%%xmm12		,%%xmm8 \n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t						addpd	%%xmm13		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t				movaps	%%xmm8 		,    (%%r10)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t				movaps	%%xmm9 		,0x10(%%r10)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t						mulpd	%%xmm6		,%%xmm12\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t						mulpd	%%xmm6		,%%xmm13\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t						mulpd	%%xmm7		,%%xmm10\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t						mulpd	%%xmm7		,%%xmm11\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t						addpd	%%xmm8 		,%%xmm12\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t						addpd	%%xmm9 		,%%xmm13\n\t"\
			"addq	$0x100		,%%rdi 	/* __o8 = s1p22r */\n\t	addq	$0x100		,%%r11 		/* __o8 = s1p13r */\n\t"\
			"subq	$0x300		,%%rsi 	/* __o2 = s1p10r */\n\t	subq	$0x300		,%%r10 		/* __o2 = s1p01r */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t						movaps	%%xmm12		,%%xmm8 \n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t						movaps	%%xmm13		,%%xmm9 \n\t"\
			"addpd	%%xmm3		,%%xmm4\n\t						addpd	%%xmm11		,%%xmm12\n\t"\
			"subpd	%%xmm2		,%%xmm5\n\t						subpd	%%xmm10		,%%xmm13\n\t"\
			"subpd	%%xmm3		,%%xmm0\n\t						subpd	%%xmm11		,%%xmm8 \n\t"\
			"addpd	%%xmm2		,%%xmm1\n\t						addpd	%%xmm10		,%%xmm9 \n\t"\
			"movaps	%%xmm4		,    (%%rdi)\n\t				movaps	%%xmm12		,    (%%r11)\n\t"\
			"movaps	%%xmm5		,0x10(%%rdi)\n\t				movaps	%%xmm13		,0x10(%%r11)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t				movaps	%%xmm8 		,    (%%r10)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t				movaps	%%xmm9 		,0x10(%%r10)\n\t"\
			"\n\t"\
			:					/* outputs: none */\
			: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
			 ,[__p01] "m" (Xp01)\
			 ,[__p02] "m" (Xp02)\
			 ,[__p03] "m" (Xp03)\
			 ,[__p04] "m" (Xp04)\
			 ,[__p08] "m" (Xp08)\
			 ,[__p12] "m" (Xp12)\
			 ,[__p16] "m" (Xp16)\
			 ,[__p20] "m" (Xp20)\
			 ,[__p24] "m" (Xp24)\
			 ,[__p28] "m" (Xp28)\
			 ,[__p32] "m" (Xp32)\
			 ,[__r00] "m" (Xr00)\
			 ,[__out] "m" (Xout)\
			 ,[__cc1] "m" (Xcc1)\
			: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
	}


	#define	SSE2_RADIX36_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xp28,Xp32,Xr00,Xin0,Xcc1)\
	{\
	__asm__ volatile (\
			/* SSE2_RADIX_09_DIF_0TWIDDLE(r00,s1p[00,32,28,24,20,16,12,08,04]r: */\
			"movq	%[__in0]	,%%rax 		/* __i0-8; e[abc]x store input addresses */\n\t"\
		"movq	%[__cc1]	,%%rdx 		/* edx stores trig addresses throughout */\n\t"\
			"movq	%[__r00]	,%%rsi\n\t"							/* SSE2_RADIX_09_DIF_0TWIDDLE(r10,s1p[27,23,19,15,11,07,03,35,31]r: */\
			"movq	%%rax		,%%rbx\n\t							leaq	0x360(%%rax),%%r10 		/* __i0 = s1p27r */\n\t"\
			"movq	%%rax		,%%rcx\n\t							movq	%%r10		,%%r11\n\t"\
			"addq	$0x300		,%%rbx 	/* __i3 = s1p24r */\n\t		movq	%%r10		,%%r12\n\t"\
			"addq	$0x180		,%%rcx 	/* __i6 = s1p12r */\n\t		subq	$0x180		,%%r11 		/* __i3 = s1p15r */\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t				subq	$0x300		,%%r12 		/* __i6 = s1p03r */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t							movaps	    (%%r11)	,%%xmm10\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t							movaps	0x10(%%r11)	,%%xmm11\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t							movaps	    (%%r10)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t							movaps	0x10(%%r10)	,%%xmm9 \n\t"\
			"movaps	    (%%rcx)	,%%xmm6\n\t							movaps	    (%%r12)	,%%xmm14\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm7\n\t							movaps	0x10(%%r12)	,%%xmm15\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t							movaps	%%xmm10		,%%xmm12\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t							movaps	%%xmm11		,%%xmm13\n\t"\
		"movq	%%rsi		,%%rdi\n\t"\
		"addq	$0x20		,%%rdi\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t							addpd	%%xmm14		,%%xmm10\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t							addpd	%%xmm15		,%%xmm11\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t							subpd	%%xmm14		,%%xmm12\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t							subpd	%%xmm15		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t							addpd	%%xmm10		,%%xmm8 \n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t							addpd	%%xmm11		,%%xmm9 \n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t00 */\n\t	movaps	%%xmm8 		,0x120(%%rsi)	/* <- t00 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t01 */\n\t	movaps	%%xmm9 		,0x130(%%rsi)	/* <- t01 */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t							mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t							mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t							mulpd	%%xmm7		,%%xmm12\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t							mulpd	%%xmm7		,%%xmm13\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t							addpd	%%xmm8 		,%%xmm10\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t							addpd	%%xmm9 		,%%xmm11\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t							movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t							movaps	%%xmm11		,%%xmm9 \n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t							subpd	%%xmm13		,%%xmm10\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t							addpd	%%xmm12		,%%xmm11\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t							addpd	%%xmm13		,%%xmm8 \n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t							subpd	%%xmm12		,%%xmm9 \n\t"\
			"movaps	%%xmm2		,    (%%rdi)	/* <- t02 */\n\t	movaps	%%xmm10		,0x120(%%rdi)	/* <- t02 */\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)	/* <- t03 */\n\t	movaps	%%xmm11		,0x130(%%rdi)	/* <- t03 */\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t04 */\n\t	movaps	%%xmm8 		,0x120(%%rsi)	/* <- t04 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t05 */\n\t	movaps	%%xmm9 		,0x130(%%rsi)	/* <- t05 */\n\t"\
			"addq	$0x400		,%%rax 	/* __i1 = s1p32r */\n\t		subq	$0x080		,%%r10 		/* __i1 = s1p23r */\n\t"\
			"subq	$0x080		,%%rbx 	/* __i4 = s1p20r */\n\t		subq	$0x080		,%%r11 		/* __i4 = s1p11r */\n\t"\
			"subq	$0x080		,%%rcx 	/* __i7 = s1p08r */\n\t		addq	$0x400		,%%r12 		/* __i7 = s1p35r */\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t							movaps	    (%%r11)	,%%xmm12\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t							movaps	0x10(%%r11)	,%%xmm13\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t							movaps	    (%%r10)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t							movaps	0x10(%%r10)	,%%xmm9 \n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t							movaps	    (%%r12)	,%%xmm10\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t							movaps	0x10(%%r12)	,%%xmm11\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t							subpd	%%xmm10		,%%xmm12\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t							subpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t							addpd	%%xmm10		,%%xmm10\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t							addpd	%%xmm11		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t							addpd	%%xmm12		,%%xmm10\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t							addpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t							addpd	%%xmm10		,%%xmm8 \n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t							addpd	%%xmm11		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rdi)	/* <- t06 */\n\t	movaps	%%xmm8 		,0x120(%%rdi)	/* <- t06 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rdi)	/* <- t07 */\n\t	movaps	%%xmm9 		,0x130(%%rdi)	/* <- t07 */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t							mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t							mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t							mulpd	%%xmm7		,%%xmm12\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t							mulpd	%%xmm7		,%%xmm13\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t							addpd	%%xmm8 		,%%xmm10\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t							addpd	%%xmm9 		,%%xmm11\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t							movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t							movaps	%%xmm11		,%%xmm9 \n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t							subpd	%%xmm13		,%%xmm10\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t							addpd	%%xmm12		,%%xmm11\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t							addpd	%%xmm13		,%%xmm8 \n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t							subpd	%%xmm12		,%%xmm9 \n\t"\
			"movaps	%%xmm2		,    (%%rsi)	/* <- t08 */\n\t	movaps	%%xmm10		,0x120(%%rsi)	/* <- t08 */\n\t"\
			"movaps	%%xmm3		,0x10(%%rsi)	/* <- t09 */\n\t	movaps	%%xmm11		,0x130(%%rsi)	/* <- t09 */\n\t"\
			"movaps	%%xmm0		,    (%%rdi)	/* <- t0a */\n\t	movaps	%%xmm8 		,0x120(%%rdi)	/* <- t0a */\n\t"\
			"movaps	%%xmm1		,0x10(%%rdi)	/* <- t0b */\n\t	movaps	%%xmm9 		,0x130(%%rdi)	/* <- t0b */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"subq	$0x080		,%%rax 	/* __i2 = s1p28r */\n\t		subq	$0x080		,%%r10 		/* __i2 = s1p19r */\n\t"\
			"subq	$0x080		,%%rbx 	/* __i5 = s1p16r */\n\t		subq	$0x080		,%%r11 		/* __i5 = s1p07r */\n\t"\
			"subq	$0x080		,%%rcx 	/* __i8 = s1p04r */\n\t		subq	$0x080		,%%r12 		/* __i8 = s1p31r */\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t							movaps	    (%%r11)	,%%xmm12\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t							movaps	0x10(%%r11)	,%%xmm13\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t							movaps	    (%%r10)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t							movaps	0x10(%%r10)	,%%xmm9 \n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t							movaps	    (%%r12)	,%%xmm10\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t							movaps	0x10(%%r12)	,%%xmm11\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t							subpd	%%xmm10		,%%xmm12\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t							subpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t							addpd	%%xmm10		,%%xmm10\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t							addpd	%%xmm11		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t							addpd	%%xmm12		,%%xmm10\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t							addpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t							addpd	%%xmm10		,%%xmm8 \n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t							addpd	%%xmm11		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t0c */\n\t	movaps	%%xmm8 		,0x120(%%rsi)	/* <- t0c */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t0d */\n\t	movaps	%%xmm9 		,0x130(%%rsi)	/* <- t0d */\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t							mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t							mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t							mulpd	%%xmm7		,%%xmm12\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t							mulpd	%%xmm7		,%%xmm13\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t							addpd	%%xmm8 		,%%xmm10\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t							addpd	%%xmm9 		,%%xmm11\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t							movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t							movaps	%%xmm11		,%%xmm9 \n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t							subpd	%%xmm13		,%%xmm10\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t							addpd	%%xmm12		,%%xmm11\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t							addpd	%%xmm13		,%%xmm8 \n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t							subpd	%%xmm12		,%%xmm9 \n\t"\
			"movaps	%%xmm2		,    (%%rdi)	/* <- t0e */\n\t	movaps	%%xmm10		,0x120(%%rdi)	/* <- t0e */\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)	/* <- t0f */\n\t	movaps	%%xmm11		,0x130(%%rdi)	/* <- t0f */\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t0g */\n\t	movaps	%%xmm8 		,0x120(%%rsi)	/* <- t0g */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t0h */\n\t	movaps	%%xmm9 		,0x130(%%rsi)	/* <- t0h */\n\t"\
			"/*****************************************/\n\t"\
			"/* now do three more radix-3 transforms: */\n\t"\
			"/*****************************************/\n\t"\
			"movq	%[__r00]	,%%rax 		/* __r00 */	\n\t"\
			"movq	%%rax		,%%rbx					\n\t"\
			"movq	%%rax		,%%rcx					\n\t"\
			"addq	$0x60		,%%rbx 		/* __r06 */	\n\t"\
			"addq	$0xc0		,%%rcx 		/* __r0c */	\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t							movaps	0x120(%%rbx)	,%%xmm12\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t							movaps	0x130(%%rbx)	,%%xmm13\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t							movaps	0x120(%%rcx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t							movaps	0x130(%%rcx)	,%%xmm11\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t							movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t							movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t							subpd	%%xmm10		,%%xmm12\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t							subpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t							addpd	%%xmm10		,%%xmm10\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t							addpd	%%xmm11		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t							addpd	%%xmm12		,%%xmm10\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t							addpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t							addpd	%%xmm10		,%%xmm8 \n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t							addpd	%%xmm11		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t					movaps	%%xmm8 		,0x120(%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t					movaps	%%xmm9 		,0x130(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t							mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t							mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t							mulpd	%%xmm7		,%%xmm12\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t							mulpd	%%xmm7		,%%xmm13\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t							addpd	%%xmm8 		,%%xmm10\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t							addpd	%%xmm9 		,%%xmm11\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t							movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t							movaps	%%xmm11		,%%xmm9 \n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t							subpd	%%xmm13		,%%xmm10\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t							addpd	%%xmm12		,%%xmm11\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t							addpd	%%xmm13		,%%xmm8 \n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t							subpd	%%xmm12		,%%xmm9 \n\t"\
			"movaps	%%xmm2		,    (%%rbx)\n\t					movaps	%%xmm10		,0x120(%%rbx)\n\t"\
			"movaps	%%xmm3		,0x10(%%rbx)\n\t					movaps	%%xmm11		,0x130(%%rbx)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t					movaps	%%xmm8 		,0x120(%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t					movaps	%%xmm9 		,0x130(%%rcx)\n\t"\
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r02 */\n\t"\
			"addq	$0x20		,%%rbx 		/* __r08 */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r0e */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t							movaps	0x120(%%rbx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t							movaps	0x130(%%rbx)	,%%xmm11\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t							movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t							movaps	%%xmm11		,%%xmm9 \n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t							mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t							mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t							mulpd	%%xmm7		,%%xmm8 \n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t							mulpd	%%xmm7		,%%xmm9 \n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t							subpd	%%xmm9 		,%%xmm10\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t							addpd	%%xmm8 		,%%xmm11\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t							movaps	0x120(%%rcx)	,%%xmm12\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t							movaps	0x130(%%rcx)	,%%xmm13\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t							movaps	%%xmm12		,%%xmm8 \n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t							movaps	%%xmm13		,%%xmm9 \n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t							mulpd	%%xmm6		,%%xmm12\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t							mulpd	%%xmm6		,%%xmm13\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t							mulpd	%%xmm7		,%%xmm8 \n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t							mulpd	%%xmm7		,%%xmm9 \n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t							subpd	%%xmm9 		,%%xmm12\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t							addpd	%%xmm8 		,%%xmm13\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t							movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t							movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t							subpd	%%xmm12		,%%xmm10\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t							subpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t							addpd	%%xmm12		,%%xmm12\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t							addpd	%%xmm13		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t							addpd	%%xmm10		,%%xmm12\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t							addpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t							addpd	%%xmm12		,%%xmm8 \n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t							addpd	%%xmm13		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t					movaps	%%xmm8 		,0x120(%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t					movaps	%%xmm9 		,0x130(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t							mulpd	%%xmm6		,%%xmm12\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t							mulpd	%%xmm6		,%%xmm13\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t							mulpd	%%xmm7		,%%xmm10\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t							mulpd	%%xmm7		,%%xmm11\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t							addpd	%%xmm8 		,%%xmm12\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t							addpd	%%xmm9 		,%%xmm13\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t							movaps	%%xmm12		,%%xmm8 \n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t							movaps	%%xmm13		,%%xmm9 \n\t"\
			"subpd	%%xmm3		,%%xmm4\n\t							subpd	%%xmm11		,%%xmm12\n\t"\
			"addpd	%%xmm2		,%%xmm5\n\t							addpd	%%xmm10		,%%xmm13\n\t"\
			"addpd	%%xmm3		,%%xmm0\n\t							addpd	%%xmm11		,%%xmm8 \n\t"\
			"subpd	%%xmm2		,%%xmm1\n\t							subpd	%%xmm10		,%%xmm9 \n\t"\
			"movaps	%%xmm4		,    (%%rbx)\n\t					movaps	%%xmm12		,0x120(%%rbx)\n\t"\
			"movaps	%%xmm5		,0x10(%%rbx)\n\t					movaps	%%xmm13		,0x130(%%rbx)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t					movaps	%%xmm8 		,0x120(%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t					movaps	%%xmm9 		,0x130(%%rcx)\n\t"\
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r04 */\n\t"\
			"addq	$0x20		,%%rbx 		/* __r0a */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r0g */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t							movaps	0x120(%%rbx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t							movaps	0x130(%%rbx)	,%%xmm11\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t							movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t							movaps	%%xmm11		,%%xmm9 \n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t							mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t							mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t							mulpd	%%xmm7		,%%xmm8 \n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t							mulpd	%%xmm7		,%%xmm9 \n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t							subpd	%%xmm9 		,%%xmm10\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t							addpd	%%xmm8 		,%%xmm11\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t							movaps	0x120(%%rcx)	,%%xmm12\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t							movaps	0x130(%%rcx)	,%%xmm13\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t							movaps	%%xmm12		,%%xmm8 \n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t							movaps	%%xmm13		,%%xmm9 \n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t							mulpd	%%xmm6		,%%xmm12\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t							mulpd	%%xmm6		,%%xmm13\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t							mulpd	%%xmm7		,%%xmm8 \n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t							mulpd	%%xmm7		,%%xmm9 \n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t							subpd	%%xmm9 		,%%xmm12\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t							addpd	%%xmm8 		,%%xmm13\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t							movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t							movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t							subpd	%%xmm12		,%%xmm10\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t							subpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t							addpd	%%xmm12		,%%xmm12\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t							addpd	%%xmm13		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t							addpd	%%xmm10		,%%xmm12\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t							addpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t							addpd	%%xmm12		,%%xmm8 \n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t							addpd	%%xmm13		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t					movaps	%%xmm8 		,0x120(%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t					movaps	%%xmm9 		,0x130(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t							mulpd	%%xmm6		,%%xmm12\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t							mulpd	%%xmm6		,%%xmm13\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t							mulpd	%%xmm7		,%%xmm10\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t							mulpd	%%xmm7		,%%xmm11\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t							addpd	%%xmm8 		,%%xmm12\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t							addpd	%%xmm9 		,%%xmm13\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t							movaps	%%xmm12		,%%xmm8 \n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t							movaps	%%xmm13		,%%xmm9 \n\t"\
			"subpd	%%xmm3		,%%xmm4\n\t							subpd	%%xmm11		,%%xmm12\n\t"\
			"addpd	%%xmm2		,%%xmm5\n\t							addpd	%%xmm10		,%%xmm13\n\t"\
			"addpd	%%xmm3		,%%xmm0\n\t							addpd	%%xmm11		,%%xmm8 \n\t"\
			"subpd	%%xmm2		,%%xmm1\n\t							subpd	%%xmm10		,%%xmm9 \n\t"\
			"movaps	%%xmm4		,    (%%rbx)\n\t					movaps	%%xmm12		,0x120(%%rbx)\n\t"\
			"movaps	%%xmm5		,0x10(%%rbx)\n\t					movaps	%%xmm13		,0x130(%%rbx)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t					movaps	%%xmm8 		,0x120(%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t					movaps	%%xmm9 		,0x130(%%rcx)\n\t"\
			/* SSE2_RADIX_09_DIF_0TWIDDLE(r20,s1p[18,14,10,06,02,34,30,26,22]r: */\
			"movq	%[__in0]	,%%rax\n\t"\
		"movq	%[__cc1]	,%%rdx 	\n\t"\
			"movq	%[__r00]	,%%rsi\n\t"							/* SSE2_RADIX_09_DIF_0TWIDDLE(r30,s1p[09,05,01,33,29,25,21,17,13]r: */\
			"addq	$0x240		,%%rsi 		/* r20 */\n\t			leaq	0x120(%%rax),%%r10 		/* __i0 = s1p09r */\n\t"\
			"addq	$0x240		,%%rax 	/* __i0 = s1p18r */\n\t"\
			"movq	%%rax		,%%rbx\n\t							movq	%%r10		,%%r11\n\t"\
			"movq	%%rax		,%%rcx\n\t							movq	%%r10		,%%r12\n\t"\
			"subq	$0x180		,%%rbx 	/* __i3 = s1p06r */\n\t		addq	$0x300		,%%r11 		/* __i3 = s1p33r */\n\t"\
			"addq	$0x180		,%%rcx 	/* __i6 = s1p30r */\n\t		addq	$0x180		,%%r12 		/* __i6 = s1p21r */\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t							movaps	    (%%r11)	,%%xmm10\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t							movaps	0x10(%%r11)	,%%xmm11\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t							movaps	    (%%r10)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t							movaps	0x10(%%r10)	,%%xmm9 \n\t"\
			"movaps	    (%%rcx)	,%%xmm6\n\t							movaps	    (%%r12)	,%%xmm14\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm7\n\t							movaps	0x10(%%r12)	,%%xmm15\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t							movaps	%%xmm10		,%%xmm12\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t							movaps	%%xmm11		,%%xmm13\n\t"\
		"movq	%%rsi		,%%rdi\n\t"\
		"addq	$0x20		,%%rdi\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t							addpd	%%xmm14		,%%xmm10\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t							addpd	%%xmm15		,%%xmm11\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t							subpd	%%xmm14		,%%xmm12\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t							subpd	%%xmm15		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t							addpd	%%xmm10		,%%xmm8 \n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t							addpd	%%xmm11		,%%xmm9 \n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t00 */\n\t	movaps	%%xmm8 		,0x120(%%rsi)	/* <- t00 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t01 */\n\t	movaps	%%xmm9 		,0x130(%%rsi)	/* <- t01 */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t							mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t							mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t							mulpd	%%xmm7		,%%xmm12\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t							mulpd	%%xmm7		,%%xmm13\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t							addpd	%%xmm8 		,%%xmm10\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t							addpd	%%xmm9 		,%%xmm11\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t							movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t							movaps	%%xmm11		,%%xmm9 \n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t							subpd	%%xmm13		,%%xmm10\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t							addpd	%%xmm12		,%%xmm11\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t							addpd	%%xmm13		,%%xmm8 \n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t							subpd	%%xmm12		,%%xmm9 \n\t"\
			"movaps	%%xmm2		,    (%%rdi)	/* <- t02 */\n\t	movaps	%%xmm10		,0x120(%%rdi)	/* <- t02 */\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)	/* <- t03 */\n\t	movaps	%%xmm11		,0x130(%%rdi)	/* <- t03 */\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t04 */\n\t	movaps	%%xmm8 		,0x120(%%rsi)	/* <- t04 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t05 */\n\t	movaps	%%xmm9 		,0x130(%%rsi)	/* <- t05 */\n\t"\
			"subq	$0x080		,%%rax 	/* __i1 = s1p14r */\n\t		subq	$0x080		,%%r10 		/* __i1 = s1p05r */\n\t"\
			"subq	$0x080		,%%rbx 	/* __i4 = s1p02r */\n\t		subq	$0x080		,%%r11 		/* __i4 = s1p29r */\n\t"\
			"subq	$0x080		,%%rcx 	/* __i7 = s1p26r */\n\t		subq	$0x080		,%%r12 		/* __i7 = s1p17r */\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t							movaps	    (%%r11)	,%%xmm12\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t							movaps	0x10(%%r11)	,%%xmm13\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t							movaps	    (%%r10)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t							movaps	0x10(%%r10)	,%%xmm9 \n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t							movaps	    (%%r12)	,%%xmm10\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t							movaps	0x10(%%r12)	,%%xmm11\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t							subpd	%%xmm10		,%%xmm12\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t							subpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t							addpd	%%xmm10		,%%xmm10\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t							addpd	%%xmm11		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t							addpd	%%xmm12		,%%xmm10\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t							addpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t							addpd	%%xmm10		,%%xmm8 \n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t							addpd	%%xmm11		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rdi)	/* <- t06 */\n\t	movaps	%%xmm8 		,0x120(%%rdi)	/* <- t06 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rdi)	/* <- t07 */\n\t	movaps	%%xmm9 		,0x130(%%rdi)	/* <- t07 */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t							mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t							mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t							mulpd	%%xmm7		,%%xmm12\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t							mulpd	%%xmm7		,%%xmm13\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t							addpd	%%xmm8 		,%%xmm10\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t							addpd	%%xmm9 		,%%xmm11\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t							movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t							movaps	%%xmm11		,%%xmm9 \n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t							subpd	%%xmm13		,%%xmm10\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t							addpd	%%xmm12		,%%xmm11\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t							addpd	%%xmm13		,%%xmm8 \n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t							subpd	%%xmm12		,%%xmm9 \n\t"\
			"movaps	%%xmm2		,    (%%rsi)	/* <- t08 */\n\t	movaps	%%xmm10		,0x120(%%rsi)	/* <- t08 */\n\t"\
			"movaps	%%xmm3		,0x10(%%rsi)	/* <- t09 */\n\t	movaps	%%xmm11		,0x130(%%rsi)	/* <- t09 */\n\t"\
			"movaps	%%xmm0		,    (%%rdi)	/* <- t0a */\n\t	movaps	%%xmm8 		,0x120(%%rdi)	/* <- t0a */\n\t"\
			"movaps	%%xmm1		,0x10(%%rdi)	/* <- t0b */\n\t	movaps	%%xmm9 		,0x130(%%rdi)	/* <- t0b */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"subq	$0x080		,%%rax 	/* __i2 = s1p10r */\n\t		subq	$0x080		,%%r10 		/* __i2 = s1p01r */\n\t"\
			"addq	$0x400		,%%rbx 	/* __i5 = s1p34r */\n\t		subq	$0x080		,%%r11 		/* __i5 = s1p25r */\n\t"\
			"subq	$0x080		,%%rcx 	/* __i8 = s1p22r */\n\t		subq	$0x080		,%%r12 		/* __i8 = s1p13r */\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t							movaps	    (%%r11)	,%%xmm12\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t							movaps	0x10(%%r11)	,%%xmm13\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t							movaps	    (%%r10)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t							movaps	0x10(%%r10)	,%%xmm9 \n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t							movaps	    (%%r12)	,%%xmm10\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t							movaps	0x10(%%r12)	,%%xmm11\n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t							subpd	%%xmm10		,%%xmm12\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t							subpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t							addpd	%%xmm10		,%%xmm10\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t							addpd	%%xmm11		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t							addpd	%%xmm12		,%%xmm10\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t							addpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t							addpd	%%xmm10		,%%xmm8 \n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t							addpd	%%xmm11		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t0c */\n\t	movaps	%%xmm8 		,0x120(%%rsi)	/* <- t0c */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t0d */\n\t	movaps	%%xmm9 		,0x130(%%rsi)	/* <- t0d */\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t							mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t							mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t							mulpd	%%xmm7		,%%xmm12\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t							mulpd	%%xmm7		,%%xmm13\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t							addpd	%%xmm8 		,%%xmm10\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t							addpd	%%xmm9 		,%%xmm11\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t							movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t							movaps	%%xmm11		,%%xmm9 \n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t							subpd	%%xmm13		,%%xmm10\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t							addpd	%%xmm12		,%%xmm11\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t							addpd	%%xmm13		,%%xmm8 \n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t							subpd	%%xmm12		,%%xmm9 \n\t"\
			"movaps	%%xmm2		,    (%%rdi)	/* <- t0e */\n\t	movaps	%%xmm10		,0x120(%%rdi)	/* <- t0e */\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)	/* <- t0f */\n\t	movaps	%%xmm11		,0x130(%%rdi)	/* <- t0f */\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t0g */\n\t	movaps	%%xmm8 		,0x120(%%rsi)	/* <- t0g */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t0h */\n\t	movaps	%%xmm9 		,0x130(%%rsi)	/* <- t0h */\n\t"\
			"/*****************************************/\n\t"\
			"/* now do three more radix-3 transforms: */\n\t"\
			"/*****************************************/\n\t"\
			"movq	%[__r00]	,%%rax 		/* __r00 */	\n\t"\
			"addq	$0x240		,%%rax 		/* r20 */	\n\t"\
			"movq	%%rax		,%%rbx					\n\t"\
			"movq	%%rax		,%%rcx					\n\t"\
			"addq	$0x60		,%%rbx 		/* __r26 */	\n\t"\
			"addq	$0xc0		,%%rcx 		/* __r2c */	\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t							movaps	0x120(%%rbx)	,%%xmm12\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t							movaps	0x130(%%rbx)	,%%xmm13\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t							movaps	0x120(%%rcx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t							movaps	0x130(%%rcx)	,%%xmm11\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t							movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t							movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"subpd	%%xmm2		,%%xmm4\n\t							subpd	%%xmm10		,%%xmm12\n\t"\
			"subpd	%%xmm3		,%%xmm5\n\t							subpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm2\n\t							addpd	%%xmm10		,%%xmm10\n\t"\
			"addpd	%%xmm3		,%%xmm3\n\t							addpd	%%xmm11		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm2\n\t							addpd	%%xmm12		,%%xmm10\n\t"\
			"addpd	%%xmm5		,%%xmm3\n\t							addpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t							addpd	%%xmm10		,%%xmm8 \n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t							addpd	%%xmm11		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t					movaps	%%xmm8 		,0x120(%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t					movaps	%%xmm9 		,0x130(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t							mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t							mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm4\n\t							mulpd	%%xmm7		,%%xmm12\n\t"\
			"mulpd	%%xmm7		,%%xmm5\n\t							mulpd	%%xmm7		,%%xmm13\n\t"\
			"addpd	%%xmm0		,%%xmm2\n\t							addpd	%%xmm8 		,%%xmm10\n\t"\
			"addpd	%%xmm1		,%%xmm3\n\t							addpd	%%xmm9 		,%%xmm11\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t							movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t							movaps	%%xmm11		,%%xmm9 \n\t"\
			"subpd	%%xmm5		,%%xmm2\n\t							subpd	%%xmm13		,%%xmm10\n\t"\
			"addpd	%%xmm4		,%%xmm3\n\t							addpd	%%xmm12		,%%xmm11\n\t"\
			"addpd	%%xmm5		,%%xmm0\n\t							addpd	%%xmm13		,%%xmm8 \n\t"\
			"subpd	%%xmm4		,%%xmm1\n\t							subpd	%%xmm12		,%%xmm9 \n\t"\
			"movaps	%%xmm2		,    (%%rbx)\n\t					movaps	%%xmm10		,0x120(%%rbx)\n\t"\
			"movaps	%%xmm3		,0x10(%%rbx)\n\t					movaps	%%xmm11		,0x130(%%rbx)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t					movaps	%%xmm8 		,0x120(%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t					movaps	%%xmm9 		,0x130(%%rcx)\n\t"\
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r22 */\n\t"\
			"addq	$0x20		,%%rbx 		/* __r28 */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r2e */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t							movaps	0x120(%%rbx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t							movaps	0x130(%%rbx)	,%%xmm11\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t							movaps	    (%%rdx)	,%%xmm14\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t							movaps	0x10(%%rdx)	,%%xmm15\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t							movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t							movaps	%%xmm11		,%%xmm9 \n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t							mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t							mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t							mulpd	%%xmm7		,%%xmm8 \n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t							mulpd	%%xmm7		,%%xmm9 \n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t							subpd	%%xmm9 		,%%xmm10\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t							addpd	%%xmm8 		,%%xmm11\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t							movaps	0x120(%%rcx)	,%%xmm12\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t							movaps	0x130(%%rcx)	,%%xmm13\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t							movaps	%%xmm12		,%%xmm8 \n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t							movaps	%%xmm13		,%%xmm9 \n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t							mulpd	%%xmm6		,%%xmm12\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t							mulpd	%%xmm6		,%%xmm13\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t							mulpd	%%xmm7		,%%xmm8 \n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t							mulpd	%%xmm7		,%%xmm9 \n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t							subpd	%%xmm9 		,%%xmm12\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t							addpd	%%xmm8 		,%%xmm13\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t							movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t							movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t							subpd	%%xmm12		,%%xmm10\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t							subpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t							addpd	%%xmm12		,%%xmm12\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t							addpd	%%xmm13		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t							addpd	%%xmm10		,%%xmm12\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t							addpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t							addpd	%%xmm12		,%%xmm8 \n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t							addpd	%%xmm13		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t					movaps	%%xmm8 		,0x120(%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t					movaps	%%xmm9 		,0x130(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t							mulpd	%%xmm6		,%%xmm12\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t							mulpd	%%xmm6		,%%xmm13\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t							mulpd	%%xmm7		,%%xmm10\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t							mulpd	%%xmm7		,%%xmm11\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t							addpd	%%xmm8 		,%%xmm12\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t							addpd	%%xmm9 		,%%xmm13\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t							movaps	%%xmm12		,%%xmm8 \n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t							movaps	%%xmm13		,%%xmm9 \n\t"\
			"subpd	%%xmm3		,%%xmm4\n\t							subpd	%%xmm11		,%%xmm12\n\t"\
			"addpd	%%xmm2		,%%xmm5\n\t							addpd	%%xmm10		,%%xmm13\n\t"\
			"addpd	%%xmm3		,%%xmm0\n\t							addpd	%%xmm11		,%%xmm8 \n\t"\
			"subpd	%%xmm2		,%%xmm1\n\t							subpd	%%xmm10		,%%xmm9 \n\t"\
			"movaps	%%xmm4		,    (%%rbx)\n\t					movaps	%%xmm12		,0x120(%%rbx)\n\t"\
			"movaps	%%xmm5		,0x10(%%rbx)\n\t					movaps	%%xmm13		,0x130(%%rbx)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t					movaps	%%xmm8 		,0x120(%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t					movaps	%%xmm9 		,0x130(%%rcx)\n\t"\
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r24 */\n\t"\
			"addq	$0x20		,%%rbx 		/* __r2a */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r2g */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t							movaps	0x120(%%rbx)	,%%xmm10\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t							movaps	0x130(%%rbx)	,%%xmm11\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t							movaps	%%xmm10		,%%xmm8 \n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t							movaps	%%xmm11		,%%xmm9 \n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t							mulpd	%%xmm6		,%%xmm10\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t							mulpd	%%xmm6		,%%xmm11\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t							mulpd	%%xmm7		,%%xmm8 \n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t							mulpd	%%xmm7		,%%xmm9 \n\t"\
			"subpd	%%xmm1		,%%xmm2\n\t							subpd	%%xmm9 		,%%xmm10\n\t"\
			"addpd	%%xmm0		,%%xmm3\n\t							addpd	%%xmm8 		,%%xmm11\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t							movaps	0x120(%%rcx)	,%%xmm12\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t							movaps	0x130(%%rcx)	,%%xmm13\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t							movaps	%%xmm12		,%%xmm8 \n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t							movaps	%%xmm13		,%%xmm9 \n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t							mulpd	%%xmm6		,%%xmm12\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t							mulpd	%%xmm6		,%%xmm13\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t							mulpd	%%xmm7		,%%xmm8 \n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t							mulpd	%%xmm7		,%%xmm9 \n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t							subpd	%%xmm9 		,%%xmm12\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t							addpd	%%xmm8 		,%%xmm13\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t							movaps	0x120(%%rax)	,%%xmm8 \n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t							movaps	0x130(%%rax)	,%%xmm9 \n\t"\
			"subpd	%%xmm4		,%%xmm2\n\t							subpd	%%xmm12		,%%xmm10\n\t"\
			"subpd	%%xmm5		,%%xmm3\n\t							subpd	%%xmm13		,%%xmm11\n\t"\
			"addpd	%%xmm4		,%%xmm4\n\t							addpd	%%xmm12		,%%xmm12\n\t"\
			"addpd	%%xmm5		,%%xmm5\n\t							addpd	%%xmm13		,%%xmm13\n\t"\
			"addpd	%%xmm2		,%%xmm4\n\t							addpd	%%xmm10		,%%xmm12\n\t"\
			"addpd	%%xmm3		,%%xmm5\n\t							addpd	%%xmm11		,%%xmm13\n\t"\
			"addpd	%%xmm4		,%%xmm0\n\t							addpd	%%xmm12		,%%xmm8 \n\t"\
			"addpd	%%xmm5		,%%xmm1\n\t							addpd	%%xmm13		,%%xmm9 \n\t"\
			"movaps	%%xmm0		,    (%%rax)\n\t					movaps	%%xmm8 		,0x120(%%rax)\n\t"\
			"movaps	%%xmm1		,0x10(%%rax)\n\t					movaps	%%xmm9 		,0x130(%%rax)\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t							mulpd	%%xmm6		,%%xmm12\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t							mulpd	%%xmm6		,%%xmm13\n\t"\
			"mulpd	%%xmm7		,%%xmm2\n\t							mulpd	%%xmm7		,%%xmm10\n\t"\
			"mulpd	%%xmm7		,%%xmm3\n\t							mulpd	%%xmm7		,%%xmm11\n\t"\
			"addpd	%%xmm0		,%%xmm4\n\t							addpd	%%xmm8 		,%%xmm12\n\t"\
			"addpd	%%xmm1		,%%xmm5\n\t							addpd	%%xmm9 		,%%xmm13\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t							movaps	%%xmm12		,%%xmm8 \n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t							movaps	%%xmm13		,%%xmm9 \n\t"\
			"subpd	%%xmm3		,%%xmm4\n\t							subpd	%%xmm11		,%%xmm12\n\t"\
			"addpd	%%xmm2		,%%xmm5\n\t							addpd	%%xmm10		,%%xmm13\n\t"\
			"addpd	%%xmm3		,%%xmm0\n\t							addpd	%%xmm11		,%%xmm8 \n\t"\
			"subpd	%%xmm2		,%%xmm1\n\t							subpd	%%xmm10		,%%xmm9 \n\t"\
			"movaps	%%xmm4		,    (%%rbx)\n\t					movaps	%%xmm12		,0x120(%%rbx)\n\t"\
			"movaps	%%xmm5		,0x10(%%rbx)\n\t					movaps	%%xmm13		,0x130(%%rbx)\n\t"\
			"movaps	%%xmm0		,    (%%rcx)\n\t					movaps	%%xmm8 		,0x120(%%rcx)\n\t"\
			"movaps	%%xmm1		,0x10(%%rcx)\n\t					movaps	%%xmm9 		,0x130(%%rcx)\n\t"\
			"\n\t"\
		"/**********************************/"\
		"/*** And now do 9 radix-4 DFTs: ***/"\
		"/**********************************/"\
			"\n\t"\
			"movq	%[__add],%%rax	/* Must use eax as base address throughout, since that is preserved in SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_B */\n\t"\
			"movq	%[__add],%%rdx						\n\t"\
			/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r00, 0x120, 0x240, eax,ebx,edx,ecx) */\
			"movq	%[__r00],%%rdi 						\n\t"\
			"movaps	     (%%rdi),%%xmm4					\n\t"\
			"movaps	0x120(%%rdi),%%xmm6					\n\t"\
			"movaps	0x010(%%rdi),%%xmm5					\n\t"\
			"movaps	0x130(%%rdi),%%xmm7					\n\t"\
			"movaps	0x240(%%rdi),%%xmm0					\n\t"\
			"movaps	0x360(%%rdi),%%xmm2					\n\t"\
			"movaps	0x250(%%rdi),%%xmm1					\n\t"\
			"movaps	0x370(%%rdi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses d <-> c: */		\n\t"							/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r02, abdc) */\
			"subpd	%%xmm2,%%xmm0						\n\t							addq	$0x20	,%%rdi 	/* r02 */	\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t							movaps	     (%%rdi),%%xmm12				\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t							movaps	0x120(%%rdi),%%xmm14				\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t							movaps	0x010(%%rdi),%%xmm13				\n\t"\
			/* add0,1,3,2 = a+p00+p[0,1,2,3] */			"\n\t							movaps	0x130(%%rdi),%%xmm15				\n\t"\
			"movslq	%[__p01],%%rsi	/* rsi stores offset-multiples of p04 below */\n\t	movaps	0x240(%%rdi),%%xmm8 				\n\t"\
			"shlq	$3,%%rsi		/* Pointer offset for floating doubles */\n\t		movaps	0x360(%%rdi),%%xmm10				\n\t"\
			"addq	%%rsi,%%rdx							\n\t							movaps	0x250(%%rdi),%%xmm9 				\n\t"\
			"movq	%%rdx,%%rbx		/* add0+p01 */		\n\t							movaps	0x370(%%rdi),%%xmm11				\n\t"\
			"addq	%%rsi,%%rdx							\n\t							subpd	%%xmm8 ,%%xmm12						\n\t"\
			"movq	%%rdx,%%rcx		/* add0+p02 */		\n\t							subpd	%%xmm10,%%xmm14						\n\t"\
			"addq	%%rsi,%%rdx		/* add0+p03 */		\n\t							subpd	%%xmm9 ,%%xmm13						\n\t"\
			"movaps	%%xmm0,    (%%rbx)					\n\t							subpd	%%xmm11,%%xmm15						\n\t"\
			"movaps	%%xmm4,    (%%rdx)					\n\t							addpd	%%xmm8 ,%%xmm8 						\n\t"\
			"movaps	%%xmm1,0x10(%%rbx)					\n\t							addpd	%%xmm10,%%xmm10						\n\t"\
			"movaps	%%xmm5,0x10(%%rcx)					\n\t							addpd	%%xmm9 ,%%xmm9 						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t							addpd	%%xmm11,%%xmm11						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t							addpd	%%xmm12,%%xmm8 						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t							addpd	%%xmm14,%%xmm10						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t							addpd	%%xmm13,%%xmm9 						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t							addpd	%%xmm15,%%xmm11						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t							/* Swap output addresses d <-> c: */		\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t							subpd	%%xmm10,%%xmm8 						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t							subpd	%%xmm15,%%xmm12						\n\t"\
			"movaps	%%xmm2,    (%%rax)					\n\t							subpd	%%xmm11,%%xmm9 						\n\t"\
			"movaps	%%xmm7,    (%%rcx)					\n\t							subpd	%%xmm14,%%xmm13						\n\t"\
			"movaps	%%xmm3,0x10(%%rax)					\n\t"							/* add0,1,3,2 = a+p32+p[0,1,2,3] */\
			"movaps	%%xmm6,0x10(%%rdx)					\n\t							movslq	%[__p32],%%rsi						\n\t"\
			/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r04, cdab) */"\n\t						shlq	$3,%%rsi							\n\t"\
			"addq	$0x20	,%%rdi 	/* r04 */	\n\t									addq	%%rsi,%%rax		/* a+p32] */	\n\t"\
			"movaps	     (%%rdi),%%xmm4					\n\t							addq	%%rsi,%%rbx							\n\t"\
			"movaps	0x120(%%rdi),%%xmm6					\n\t							addq	%%rsi,%%rcx							\n\t"\
			"movaps	0x010(%%rdi),%%xmm5					\n\t							addq	%%rsi,%%rdx							\n\t"\
			"movaps	0x130(%%rdi),%%xmm7					\n\t							movaps	%%xmm8 ,    (%%rbx)					\n\t"\
			"movaps	0x240(%%rdi),%%xmm0					\n\t							movaps	%%xmm12,    (%%rdx)					\n\t"\
			"movaps	0x360(%%rdi),%%xmm2					\n\t							movaps	%%xmm9 ,0x10(%%rbx)					\n\t"\
			"movaps	0x250(%%rdi),%%xmm1					\n\t							movaps	%%xmm13,0x10(%%rcx)					\n\t"\
			"movaps	0x370(%%rdi),%%xmm3					\n\t							addpd	%%xmm10,%%xmm10						\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t							addpd	%%xmm15,%%xmm15						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t							addpd	%%xmm11,%%xmm11						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t							addpd	%%xmm14,%%xmm14						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t							addpd	%%xmm8 ,%%xmm10						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t							addpd	%%xmm12,%%xmm15						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t							addpd	%%xmm9 ,%%xmm11						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t							addpd	%%xmm13,%%xmm14						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t							movaps	%%xmm10,    (%%rax)					\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t							movaps	%%xmm15,    (%%rcx)					\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t							movaps	%%xmm11,0x10(%%rax)					\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t							movaps	%%xmm14,0x10(%%rdx)					\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"							/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r06, bacd) */\
			"/* Swap output addresses a <-> c ,b <-> d: */\n\t							addq	$0x20	,%%rdi 	/* r06 */	\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t							movaps	     (%%rdi),%%xmm12					\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t							movaps	0x120(%%rdi),%%xmm14					\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t							movaps	0x010(%%rdi),%%xmm13					\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t							movaps	0x130(%%rdi),%%xmm15					\n\t"\
			/* add2,3,0,1 = a+p20+p[0,1,2,3] */			"\n\t							movaps	0x240(%%rdi),%%xmm8 					\n\t"\
			"movslq	%[__p12],%%rsi						\n\t							movaps	0x360(%%rdi),%%xmm10					\n\t"\
			"shlq	$3,%%rsi							\n\t							movaps	0x250(%%rdi),%%xmm9 					\n\t"\
			"subq	%%rsi,%%rax		/* a+p20] */		\n\t							movaps	0x370(%%rdi),%%xmm11					\n\t"\
			"subq	%%rsi,%%rbx							\n\t							subpd	%%xmm8 ,%%xmm12						\n\t"\
			"subq	%%rsi,%%rcx							\n\t							subpd	%%xmm10,%%xmm14						\n\t"\
			"subq	%%rsi,%%rdx							\n\t							subpd	%%xmm9 ,%%xmm13						\n\t"\
			"movaps	%%xmm0,    (%%rdx)					\n\t							subpd	%%xmm11,%%xmm15						\n\t"\
			"movaps	%%xmm4,    (%%rax)					\n\t							addpd	%%xmm8 ,%%xmm8 						\n\t"\
			"movaps	%%xmm1,0x10(%%rdx)					\n\t							addpd	%%xmm10,%%xmm10						\n\t"\
			"movaps	%%xmm5,0x10(%%rbx)					\n\t							addpd	%%xmm9 ,%%xmm9 						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t							addpd	%%xmm11,%%xmm11						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t							addpd	%%xmm12,%%xmm8 						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t							addpd	%%xmm14,%%xmm10						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t							addpd	%%xmm13,%%xmm9 						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t							addpd	%%xmm15,%%xmm11						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t							/* Swap output addresses a <-> b: */		\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t							subpd	%%xmm10,%%xmm8 						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t							subpd	%%xmm15,%%xmm12						\n\t"\
			"movaps	%%xmm2,    (%%rcx)					\n\t							subpd	%%xmm11,%%xmm9 						\n\t"\
			"movaps	%%xmm7,    (%%rbx)					\n\t							subpd	%%xmm14,%%xmm13						\n\t"\
			"movaps	%%xmm3,0x10(%%rcx)					\n\t"							/* add1,0,2,3 = a+p08+p[0,1,2,3] */\
			"movaps	%%xmm6,0x10(%%rax)					\n\t							movslq	%[__p12],%%rsi						\n\t"\
			/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r08, dcba) */"\n\t						shlq	$3,%%rsi							\n\t"\
			"addq	$0x20	,%%rdi 	/* r08 */	\n\t									subq	%%rsi,%%rax		/* a+p08] */	\n\t"\
			"movaps	     (%%rdi),%%xmm4					\n\t							subq	%%rsi,%%rbx							\n\t"\
			"movaps	0x120(%%rdi),%%xmm6					\n\t							subq	%%rsi,%%rcx							\n\t"\
			"movaps	0x010(%%rdi),%%xmm5					\n\t							subq	%%rsi,%%rdx							\n\t"\
			"movaps	0x130(%%rdi),%%xmm7					\n\t							movaps	%%xmm8 ,    (%%rax)					\n\t"\
			"movaps	0x240(%%rdi),%%xmm0					\n\t							movaps	%%xmm12,    (%%rcx)					\n\t"\
			"movaps	0x360(%%rdi),%%xmm2					\n\t							movaps	%%xmm9 ,0x10(%%rax)					\n\t"\
			"movaps	0x250(%%rdi),%%xmm1					\n\t							movaps	%%xmm13,0x10(%%rdx)					\n\t"\
			"movaps	0x370(%%rdi),%%xmm3					\n\t							addpd	%%xmm10,%%xmm10						\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t							addpd	%%xmm15,%%xmm15						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t							addpd	%%xmm11,%%xmm11						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t							addpd	%%xmm14,%%xmm14						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t							addpd	%%xmm8 ,%%xmm10						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t							addpd	%%xmm12,%%xmm15						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t							addpd	%%xmm9 ,%%xmm11						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t							addpd	%%xmm13,%%xmm14						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t							movaps	%%xmm10,    (%%rbx)					\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t							movaps	%%xmm15,    (%%rdx)					\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t							movaps	%%xmm11,0x10(%%rbx)					\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t							movaps	%%xmm14,0x10(%%rcx)					\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"							/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0a, abdc) */\
			"/* Swap output addresses a <-> d ,b <-> c: */\n\t							addq	$0x20	,%%rdi 	/* r0a */	\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t							movaps	     (%%rdi),%%xmm12					\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t							movaps	0x120(%%rdi),%%xmm14					\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t							movaps	0x010(%%rdi),%%xmm13					\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t							movaps	0x130(%%rdi),%%xmm15					\n\t"\
			/* add3,2,1,0 = a+p28+p[0,1,2,3] */			"\n\t							movaps	0x240(%%rdi),%%xmm8 					\n\t"\
			"movslq	%[__p20],%%rsi						\n\t							movaps	0x360(%%rdi),%%xmm10					\n\t"\
			"shlq	$3,%%rsi							\n\t							movaps	0x250(%%rdi),%%xmm9 					\n\t"\
			"addq	%%rsi,%%rax		/* a+p28] */	\n\t								movaps	0x370(%%rdi),%%xmm11					\n\t"\
			"addq	%%rsi,%%rbx							\n\t							subpd	%%xmm8 ,%%xmm12						\n\t"\
			"addq	%%rsi,%%rcx							\n\t							subpd	%%xmm10,%%xmm14						\n\t"\
			"addq	%%rsi,%%rdx							\n\t							subpd	%%xmm9 ,%%xmm13						\n\t"\
			"movaps	%%xmm0,    (%%rcx)					\n\t							subpd	%%xmm11,%%xmm15						\n\t"\
			"movaps	%%xmm4,    (%%rbx)					\n\t							addpd	%%xmm8 ,%%xmm8 						\n\t"\
			"movaps	%%xmm1,0x10(%%rcx)					\n\t							addpd	%%xmm10,%%xmm10						\n\t"\
			"movaps	%%xmm5,0x10(%%rax)					\n\t							addpd	%%xmm9 ,%%xmm9 						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t							addpd	%%xmm11,%%xmm11						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t							addpd	%%xmm12,%%xmm8 						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t							addpd	%%xmm14,%%xmm10						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t							addpd	%%xmm13,%%xmm9 						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t							addpd	%%xmm15,%%xmm11						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t							/* Swap output addresses d <-> c: */		\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t							subpd	%%xmm10,%%xmm8 						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t							subpd	%%xmm15,%%xmm12						\n\t"\
			"movaps	%%xmm2,    (%%rdx)					\n\t							subpd	%%xmm11,%%xmm9 						\n\t"\
			"movaps	%%xmm7,    (%%rax)					\n\t							subpd	%%xmm14,%%xmm13						\n\t"\
			"movaps	%%xmm3,0x10(%%rdx)					\n\t"							/* add0,1,3,2 = a+p16+p[0,1,2,3] */\
			"movaps	%%xmm6,0x10(%%rbx)					\n\t							movslq	%[__p12],%%rsi						\n\t"\
			/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0c, cdab) */"\n\t						shlq	$3,%%rsi							\n\t"\
			"addq	$0x20	,%%rdi 	/* r0c */	\n\t									subq	%%rsi,%%rax		/* a+p16] */	\n\t"\
			"movaps	     (%%rdi),%%xmm4					\n\t							subq	%%rsi,%%rbx							\n\t"\
			"movaps	0x120(%%rdi),%%xmm6					\n\t							subq	%%rsi,%%rcx							\n\t"\
			"movaps	0x010(%%rdi),%%xmm5					\n\t							subq	%%rsi,%%rdx							\n\t"\
			"movaps	0x130(%%rdi),%%xmm7					\n\t							movaps	%%xmm8 ,    (%%rbx)					\n\t"\
			"movaps	0x240(%%rdi),%%xmm0					\n\t							movaps	%%xmm12,    (%%rdx)					\n\t"\
			"movaps	0x360(%%rdi),%%xmm2					\n\t							movaps	%%xmm9 ,0x10(%%rbx)					\n\t"\
			"movaps	0x250(%%rdi),%%xmm1					\n\t							movaps	%%xmm13,0x10(%%rcx)					\n\t"\
			"movaps	0x370(%%rdi),%%xmm3					\n\t							addpd	%%xmm10,%%xmm10						\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t							addpd	%%xmm15,%%xmm15						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t							addpd	%%xmm11,%%xmm11						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t							addpd	%%xmm14,%%xmm14						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t							addpd	%%xmm8 ,%%xmm10						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t							addpd	%%xmm12,%%xmm15						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t							addpd	%%xmm9 ,%%xmm11						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t							addpd	%%xmm13,%%xmm14						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t							movaps	%%xmm10,    (%%rax)					\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t							movaps	%%xmm15,    (%%rcx)					\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t							movaps	%%xmm11,0x10(%%rax)					\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t							movaps	%%xmm14,0x10(%%rdx)					\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"							/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0e, bacd) */\
			"/* Swap output addresses a <-> c ,b <-> d: */\n\t							addq	$0x20	,%%rdi 	/* r0e */	\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t							movaps	     (%%rdi),%%xmm12					\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t							movaps	0x120(%%rdi),%%xmm14					\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t							movaps	0x010(%%rdi),%%xmm13					\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t							movaps	0x130(%%rdi),%%xmm15					\n\t"\
			/* add2,3,0,1 = a+p04+p[0,1,2,3] */			"\n\t							movaps	0x240(%%rdi),%%xmm8 					\n\t"\
			"movslq	%[__p12],%%rsi						\n\t							movaps	0x360(%%rdi),%%xmm10					\n\t"\
			"shlq	$3,%%rsi							\n\t							movaps	0x250(%%rdi),%%xmm9 					\n\t"\
			"subq	%%rsi,%%rax		/* a+p04] */	\n\t								movaps	0x370(%%rdi),%%xmm11					\n\t"\
			"subq	%%rsi,%%rbx							\n\t							subpd	%%xmm8 ,%%xmm12						\n\t"\
			"subq	%%rsi,%%rcx							\n\t							subpd	%%xmm10,%%xmm14						\n\t"\
			"subq	%%rsi,%%rdx							\n\t							subpd	%%xmm9 ,%%xmm13						\n\t"\
			"movaps	%%xmm0,    (%%rdx)					\n\t							subpd	%%xmm11,%%xmm15						\n\t"\
			"movaps	%%xmm4,    (%%rax)					\n\t							addpd	%%xmm8 ,%%xmm8 						\n\t"\
			"movaps	%%xmm1,0x10(%%rdx)					\n\t							addpd	%%xmm10,%%xmm10						\n\t"\
			"movaps	%%xmm5,0x10(%%rbx)					\n\t							addpd	%%xmm9 ,%%xmm9 						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t							addpd	%%xmm11,%%xmm11						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t							addpd	%%xmm12,%%xmm8 						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t							addpd	%%xmm14,%%xmm10						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t							addpd	%%xmm13,%%xmm9 						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t							addpd	%%xmm15,%%xmm11						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t							/* Swap output addresses a <-> b: */		\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t							subpd	%%xmm10,%%xmm8 						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t							subpd	%%xmm15,%%xmm12						\n\t"\
			"movaps	%%xmm2,    (%%rcx)					\n\t							subpd	%%xmm11,%%xmm9 						\n\t"\
			"movaps	%%xmm7,    (%%rbx)					\n\t							subpd	%%xmm14,%%xmm13						\n\t"\
			"movaps	%%xmm3,0x10(%%rcx)					\n\t"							/* add1,0,2,3 = a+p24+p[0,1,2,3] */\
			"movaps	%%xmm6,0x10(%%rax)					\n\t							movslq	%[__p20],%%rsi						\n\t"\
			/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0g, dcba) */"\n\t						shlq	$3,%%rsi							\n\t"\
			"addq	$0x20	,%%rdi 	/* r0g */	\n\t									addq	%%rsi,%%rax		/* a+p24] */	\n\t"\
			"movaps	     (%%rdi),%%xmm4					\n\t							addq	%%rsi,%%rbx							\n\t"\
			"movaps	0x120(%%rdi),%%xmm6					\n\t							addq	%%rsi,%%rcx							\n\t"\
			"movaps	0x010(%%rdi),%%xmm5					\n\t							addq	%%rsi,%%rdx							\n\t"\
			"movaps	0x130(%%rdi),%%xmm7					\n\t							movaps	%%xmm8 ,    (%%rax)					\n\t"\
			"movaps	0x240(%%rdi),%%xmm0					\n\t							movaps	%%xmm12,    (%%rcx)					\n\t"\
			"movaps	0x360(%%rdi),%%xmm2					\n\t							movaps	%%xmm9 ,0x10(%%rax)					\n\t"\
			"movaps	0x250(%%rdi),%%xmm1					\n\t							movaps	%%xmm13,0x10(%%rdx)					\n\t"\
			"movaps	0x370(%%rdi),%%xmm3					\n\t							addpd	%%xmm10,%%xmm10						\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t							addpd	%%xmm15,%%xmm15						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t							addpd	%%xmm11,%%xmm11						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t							addpd	%%xmm14,%%xmm14						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t							addpd	%%xmm8 ,%%xmm10						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t							addpd	%%xmm12,%%xmm15						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t							addpd	%%xmm9 ,%%xmm11						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t							addpd	%%xmm13,%%xmm14						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t							movaps	%%xmm10,    (%%rbx)					\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t							movaps	%%xmm15,    (%%rdx)					\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t							movaps	%%xmm11,0x10(%%rbx)					\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t							movaps	%%xmm14,0x10(%%rcx)					\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses a <-> d ,b <-> c: */\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			/* add,2,1,03 = a+p12+p[0,1,2,3] */\
			"movslq	%[__p12],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"subq	%%rsi,%%rax		/* a+p12] */		\n\t"\
			"subq	%%rsi,%%rbx							\n\t"\
			"subq	%%rsi,%%rcx							\n\t"\
			"subq	%%rsi,%%rdx							\n\t"\
			"movaps	%%xmm0,    (%%rcx)					\n\t"\
			"movaps	%%xmm4,    (%%rbx)					\n\t"\
			"movaps	%%xmm1,0x10(%%rcx)					\n\t"\
			"movaps	%%xmm5,0x10(%%rax)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rdx)					\n\t"\
			"movaps	%%xmm7,    (%%rax)					\n\t"\
			"movaps	%%xmm3,0x10(%%rdx)					\n\t"\
			"movaps	%%xmm6,0x10(%%rbx)					\n\t"\
			"\n\t"\
			:					/* outputs: none */\
			: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
			 ,[__p01] "m" (Xp01)\
			 ,[__p02] "m" (Xp02)\
			 ,[__p03] "m" (Xp03)\
			 ,[__p04] "m" (Xp04)\
			 ,[__p08] "m" (Xp08)\
			 ,[__p12] "m" (Xp12)\
			 ,[__p16] "m" (Xp16)\
			 ,[__p20] "m" (Xp20)\
			 ,[__p24] "m" (Xp24)\
			 ,[__p28] "m" (Xp28)\
			 ,[__p32] "m" (Xp32)\
			 ,[__r00] "m" (Xr00)\
			 ,[__in0] "m" (Xin0)\
			 ,[__cc1] "m" (Xcc1)\
			: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","r10","r11","r12","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
	}

  #else // USE_64BIT_ASM_STYLE = False: Use simple 64-bit-ified version of the analogous 32-bit ASM macros, i.e. using just xmm0-7.

	#define	SSE2_RADIX36_DIT_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xp28,Xp32,Xr00,Xout,Xcc1)\
	{\
	__asm__ volatile (\
		/*	add0,1,3,2 = a+p00+p[0,1,2,3]: SSE2_RADIX4_DIT_0TWIDDLE_STRIDE_C(eax,ebx,edx,ecx, 0x120, 0x240, r00: */\
			"movq	%[__r00],%%rsi\n\t"\
			"movq	%[__add],%%rax\n\t"\
			"movslq	%[__p01],%%rbx\n\t"\
			"movslq	%[__p02],%%rcx\n\t"\
			"movslq	%[__p03],%%rdx\n\t"\
			"shlq	$3,%%rbx		/* Pointer offset for floating doubles */\n\t"\
			"shlq	$3,%%rcx\n\t"\
			"shlq	$3,%%rdx\n\t"\
			"addq	%%rax,%%rbx\n\t"\
			"addq	%%rax,%%rcx\n\t"\
			"addq	%%rax,%%rdx\n\t"\
			"/* ecx <-> edx */\n\t"\
			"movaps	     (%%rax),%%xmm2	\n\t"\
			"movaps	     (%%rdx),%%xmm6	\n\t"\
			"movaps	0x010(%%rax),%%xmm3	\n\t"\
			"movaps	0x010(%%rdx),%%xmm7	\n\t"\
			"movaps	     (%%rbx),%%xmm0	\n\t"\
			"movaps	     (%%rcx),%%xmm4	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t"\
			"movaps	0x010(%%rcx),%%xmm5	\n\t"\
			"subpd	%%xmm0,%%xmm2		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t"\
			"subpd	%%xmm1,%%xmm3		\n\t"\
			"subpd	%%xmm5,%%xmm7		\n\t"\
			"addpd	%%xmm0,%%xmm0		\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm1,%%xmm1		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t"\
			"addpd	%%xmm7,%%xmm5		\n\t"\
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
		/*	add2,3,0,1 = a+p04+p[0,1,2,3]: cdab, r04: */\
			"addq	$0x040,%%rsi	/* r04 */\n\t"\
			"movslq	%[__p04],%%rdi\n\t"\
			"shlq	$3,%%rdi \n\t"\
			"addq	%%rdi,%%rax	/* add04 */\n\t"\
			"addq	%%rdi,%%rbx\n\t"\
			"addq	%%rdi,%%rcx\n\t"\
			"addq	%%rdi,%%rdx\n\t"\
			"/* e[ab]x <-> e[cd]x */\n\t"\
			"movaps	     (%%rcx),%%xmm2	\n\t"\
			"movaps	     (%%rax),%%xmm6	\n\t"\
			"movaps	0x010(%%rcx),%%xmm3	\n\t"\
			"movaps	0x010(%%rax),%%xmm7	\n\t"\
			"movaps	     (%%rdx),%%xmm0	\n\t"\
			"movaps	     (%%rbx),%%xmm4	\n\t"\
			"movaps	0x010(%%rdx),%%xmm1	\n\t"\
			"movaps	0x010(%%rbx),%%xmm5	\n\t"\
			"subpd	%%xmm0,%%xmm2		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t"\
			"subpd	%%xmm1,%%xmm3		\n\t"\
			"subpd	%%xmm5,%%xmm7		\n\t"\
			"addpd	%%xmm0,%%xmm0		\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm1,%%xmm1		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t"\
			"addpd	%%xmm7,%%xmm5		\n\t"\
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
		/*	add1,0,2,3 = a+p08+p[0,1,2,3]: bacd, r02: */\
			"subq	$0x020,%%rsi	/* r02 */\n\t"\
			"addq	%%rdi,%%rax	/* add08 */\n\t"\
			"addq	%%rdi,%%rbx\n\t"\
			"addq	%%rdi,%%rcx\n\t"\
			"addq	%%rdi,%%rdx\n\t"\
			"/* eax <-> ebx */\n\t"\
			"movaps	     (%%rbx),%%xmm2	\n\t"\
			"movaps	     (%%rcx),%%xmm6	\n\t"\
			"movaps	0x010(%%rbx),%%xmm3	\n\t"\
			"movaps	0x010(%%rcx),%%xmm7	\n\t"\
			"movaps	     (%%rax),%%xmm0	\n\t"\
			"movaps	     (%%rdx),%%xmm4	\n\t"\
			"movaps	0x010(%%rax),%%xmm1	\n\t"\
			"movaps	0x010(%%rdx),%%xmm5	\n\t"\
			"subpd	%%xmm0,%%xmm2		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t"\
			"subpd	%%xmm1,%%xmm3		\n\t"\
			"subpd	%%xmm5,%%xmm7		\n\t"\
			"addpd	%%xmm0,%%xmm0		\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm1,%%xmm1		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t"\
			"addpd	%%xmm7,%%xmm5		\n\t"\
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
		/*	add3,2,1,0 = a+p12+p[0,1,2,3]: dcba, r0c: */\
			"addq	$0x0a0,%%rsi	/* r0c */\n\t"\
			"addq	%%rdi,%%rax	/* add12 */\n\t"\
			"addq	%%rdi,%%rbx\n\t"\
			"addq	%%rdi,%%rcx\n\t"\
			"addq	%%rdi,%%rdx\n\t"\
			"/* eax <-> edx, ebx <-> ecx */\n\t"\
			"movaps	     (%%rdx),%%xmm2	\n\t"\
			"movaps	     (%%rbx),%%xmm6	\n\t"\
			"movaps	0x010(%%rdx),%%xmm3	\n\t"\
			"movaps	0x010(%%rbx),%%xmm7	\n\t"\
			"movaps	     (%%rcx),%%xmm0	\n\t"\
			"movaps	     (%%rax),%%xmm4	\n\t"\
			"movaps	0x010(%%rcx),%%xmm1	\n\t"\
			"movaps	0x010(%%rax),%%xmm5	\n\t"\
			"subpd	%%xmm0,%%xmm2		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t"\
			"subpd	%%xmm1,%%xmm3		\n\t"\
			"subpd	%%xmm5,%%xmm7		\n\t"\
			"addpd	%%xmm0,%%xmm0		\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm1,%%xmm1		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t"\
			"addpd	%%xmm7,%%xmm5		\n\t"\
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
		/*	add0,1,3,2 = a+p16+p[0,1,2,3]: abdc, r0g: */\
			"addq	$0x040,%%rsi	/* r0g */\n\t"\
			"addq	%%rdi,%%rax	/* add16 */\n\t"\
			"addq	%%rdi,%%rbx\n\t"\
			"addq	%%rdi,%%rcx\n\t"\
			"addq	%%rdi,%%rdx\n\t"\
			"/* ecx <-> edx */\n\t"\
			"movaps	     (%%rax),%%xmm2	\n\t"\
			"movaps	     (%%rdx),%%xmm6	\n\t"\
			"movaps	0x010(%%rax),%%xmm3	\n\t"\
			"movaps	0x010(%%rdx),%%xmm7	\n\t"\
			"movaps	     (%%rbx),%%xmm0	\n\t"\
			"movaps	     (%%rcx),%%xmm4	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t"\
			"movaps	0x010(%%rcx),%%xmm5	\n\t"\
			"subpd	%%xmm0,%%xmm2		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t"\
			"subpd	%%xmm1,%%xmm3		\n\t"\
			"subpd	%%xmm5,%%xmm7		\n\t"\
			"addpd	%%xmm0,%%xmm0		\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm1,%%xmm1		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t"\
			"addpd	%%xmm7,%%xmm5		\n\t"\
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
		/*	add2,3,0,1 = a+p20+p[0,1,2,3]: cdab, r0e: */\
			"subq	$0x020,%%rsi	/* r0e */\n\t"\
			"addq	%%rdi,%%rax	/* add20 */\n\t"\
			"addq	%%rdi,%%rbx\n\t"\
			"addq	%%rdi,%%rcx\n\t"\
			"addq	%%rdi,%%rdx\n\t"\
			"/* e[ab]x <-> e[cd]x */\n\t"\
			"movaps	     (%%rcx),%%xmm2	\n\t"\
			"movaps	     (%%rax),%%xmm6	\n\t"\
			"movaps	0x010(%%rcx),%%xmm3	\n\t"\
			"movaps	0x010(%%rax),%%xmm7	\n\t"\
			"movaps	     (%%rdx),%%xmm0	\n\t"\
			"movaps	     (%%rbx),%%xmm4	\n\t"\
			"movaps	0x010(%%rdx),%%xmm1	\n\t"\
			"movaps	0x010(%%rbx),%%xmm5	\n\t"\
			"subpd	%%xmm0,%%xmm2		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t"\
			"subpd	%%xmm1,%%xmm3		\n\t"\
			"subpd	%%xmm5,%%xmm7		\n\t"\
			"addpd	%%xmm0,%%xmm0		\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm1,%%xmm1		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t"\
			"addpd	%%xmm7,%%xmm5		\n\t"\
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
		/*	add1,0,2,3 = a+p24+p[0,1,2,3]: bacd, r08: */\
			"subq	$0x060,%%rsi	/* r08 */\n\t"\
			"addq	%%rdi,%%rax	/* add24 */\n\t"\
			"addq	%%rdi,%%rbx\n\t"\
			"addq	%%rdi,%%rcx\n\t"\
			"addq	%%rdi,%%rdx\n\t"\
			"/* eax <-> ebx */\n\t"\
			"movaps	     (%%rbx),%%xmm2	\n\t"\
			"movaps	     (%%rcx),%%xmm6	\n\t"\
			"movaps	0x010(%%rbx),%%xmm3	\n\t"\
			"movaps	0x010(%%rcx),%%xmm7	\n\t"\
			"movaps	     (%%rax),%%xmm0	\n\t"\
			"movaps	     (%%rdx),%%xmm4	\n\t"\
			"movaps	0x010(%%rax),%%xmm1	\n\t"\
			"movaps	0x010(%%rdx),%%xmm5	\n\t"\
			"subpd	%%xmm0,%%xmm2		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t"\
			"subpd	%%xmm1,%%xmm3		\n\t"\
			"subpd	%%xmm5,%%xmm7		\n\t"\
			"addpd	%%xmm0,%%xmm0		\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm1,%%xmm1		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t"\
			"addpd	%%xmm7,%%xmm5		\n\t"\
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
		/*	add3,2,1,0 = a+p28+p[0,1,2,3]: dcba, r06: */\
			"subq	$0x020,%%rsi	/* r04 */\n\t"\
			"addq	%%rdi,%%rax	/* add28 */\n\t"\
			"addq	%%rdi,%%rbx\n\t"\
			"addq	%%rdi,%%rcx\n\t"\
			"addq	%%rdi,%%rdx\n\t"\
			"/* eax <-> edx, ebx <-> ecx */\n\t"\
			"movaps	     (%%rdx),%%xmm2	\n\t"\
			"movaps	     (%%rbx),%%xmm6	\n\t"\
			"movaps	0x010(%%rdx),%%xmm3	\n\t"\
			"movaps	0x010(%%rbx),%%xmm7	\n\t"\
			"movaps	     (%%rcx),%%xmm0	\n\t"\
			"movaps	     (%%rax),%%xmm4	\n\t"\
			"movaps	0x010(%%rcx),%%xmm1	\n\t"\
			"movaps	0x010(%%rax),%%xmm5	\n\t"\
			"subpd	%%xmm0,%%xmm2		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t"\
			"subpd	%%xmm1,%%xmm3		\n\t"\
			"subpd	%%xmm5,%%xmm7		\n\t"\
			"addpd	%%xmm0,%%xmm0		\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm1,%%xmm1		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t"\
			"addpd	%%xmm7,%%xmm5		\n\t"\
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
		/*	add0,1,3,2 = a+p32+p[0,1,2,3]: abdc, r0a: */\
			"addq	$0x040,%%rsi	/* r04 */\n\t"\
			"addq	%%rdi,%%rax	/* add32 */\n\t"\
			"addq	%%rdi,%%rbx\n\t"\
			"addq	%%rdi,%%rcx\n\t"\
			"addq	%%rdi,%%rdx\n\t"\
			"/* ecx <-> edx */\n\t"\
			"movaps	     (%%rax),%%xmm2	\n\t"\
			"movaps	     (%%rdx),%%xmm6	\n\t"\
			"movaps	0x010(%%rax),%%xmm3	\n\t"\
			"movaps	0x010(%%rdx),%%xmm7	\n\t"\
			"movaps	     (%%rbx),%%xmm0	\n\t"\
			"movaps	     (%%rcx),%%xmm4	\n\t"\
			"movaps	0x010(%%rbx),%%xmm1	\n\t"\
			"movaps	0x010(%%rcx),%%xmm5	\n\t"\
			"subpd	%%xmm0,%%xmm2		\n\t"\
			"subpd	%%xmm4,%%xmm6		\n\t"\
			"subpd	%%xmm1,%%xmm3		\n\t"\
			"subpd	%%xmm5,%%xmm7		\n\t"\
			"addpd	%%xmm0,%%xmm0		\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm1,%%xmm1		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm2,%%xmm0		\n\t"\
			"addpd	%%xmm6,%%xmm4		\n\t"\
			"addpd	%%xmm3,%%xmm1		\n\t"\
			"addpd	%%xmm7,%%xmm5		\n\t"\
			"subpd	%%xmm4,%%xmm0		\n\t"\
			"subpd	%%xmm7,%%xmm2		\n\t"\
			"subpd	%%xmm5,%%xmm1		\n\t"\
			"subpd	%%xmm6,%%xmm3		\n\t"\
			"movaps	%%xmm0,0x240(%%rsi)	\n\t"\
			"movaps	%%xmm2,0x360(%%rsi)	\n\t"\
			"movaps	%%xmm1,0x250(%%rsi)	\n\t"\
			"movaps	%%xmm3,0x130(%%rsi)	\n\t"\
			"addpd	%%xmm4,%%xmm4		\n\t"\
			"addpd	%%xmm7,%%xmm7		\n\t"\
			"addpd	%%xmm5,%%xmm5		\n\t"\
			"addpd	%%xmm6,%%xmm6		\n\t"\
			"addpd	%%xmm0,%%xmm4		\n\t"\
			"addpd	%%xmm2,%%xmm7		\n\t"\
			"addpd	%%xmm1,%%xmm5		\n\t"\
			"addpd	%%xmm3,%%xmm6		\n\t"\
			"movaps	%%xmm4,     (%%rsi)	\n\t"\
			"movaps	%%xmm7,0x120(%%rsi)	\n\t"\
			"movaps	%%xmm5,0x010(%%rsi)	\n\t"\
			"movaps	%%xmm6,0x370(%%rsi)	\n\t"\
			"\n\t"\
			"/***************************************/\n\t"\
			"/*...and now do 4 radix-9 transforms...*/\n\t"\
			"/***************************************/\n\t"\
			"\n\t"\
			/* SSE2_RADIX_09_DIT_0TWIDDLE(r00,s1p[00,32,28,24,20,16,12,08,04]r): */\
			"movq	%[__out]	,%%rsi 		/* __o0-8: esi,edi store output addresses throughout */\n\t"\
			"movq	%%rsi		,%%rdi\n\t"\
			"movq	%[__r00]	,%%rax\n\t"\
		"movq	%[__cc1]	,%%rdx 		/* edx stores trig addresses throughout */\n\t"\
			"movq	%%rax		,%%rbx\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x20		,%%rbx\n\t"\
			"addq	$0x40		,%%rcx\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
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
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
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
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%rbx\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
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
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%rbx\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
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
			"/*****************************************/\n\t"\
			"/* now do three more radix-3 transforms: */\n\t"\
			"/*****************************************/\n\t"\
			"subq	$0xc0		,%%rax 		/* __r00 */\n\t"\
			"subq	$0x80,%%rbx /* ebx stores r02+0xc0 = r00+0xe0; need r06 = r00+0x60 = r02+0x40 = ebx-0x80 */\n\t"\
			"subq	$0x40,%%rcx /* ecx stores r04+0xc0           ; need r0c = r00+0xc0 = r04+0x80 = ebx-0x40 */\n\t"\
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
			"addq	$0x300		,%%rdi 		/* __o3 = s1p24r */\n\t"\
			"addq	$0x180		,%%rsi 		/* __o6 = s1p12r */\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rdi)\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%rbx\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subq	$0x080		,%%rsi 		/* __o7 = s1p08r */\n\t"\
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
			"addq	$0x100		,%%rdi 		/* __o1 = s1p32r */\n\t"\
			"addq	$0x180		,%%rsi 		/* __o4 = s1p20r */\n\t"\
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
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%rbx\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subq	$0x080		,%%rsi 		/* __o5 = s1p16r */\n\t"\
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
			"subq	$0x380		,%%rdi 		/* __o8 = s1p04r */\n\t"\
			"addq	$0x180		,%%rsi 		/* __o2 = s1p28r */\n\t"\
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
			/* SSE2_RADIX_09_DIT_0TWIDDLE(r10,s1p[27,23,19,15,11,07,03,35,31]r): */\
			"movq	%[__out]	,%%rsi\n\t"\
			"addq	$0x360		,%%rsi 		/* s1p27 */\n\t"\
			"movq	%%rsi		,%%rdi\n\t"\
			"movq	%[__r00]	,%%rax\n\t"\
			"addq	$0x120		,%%rax 		/* __r10 */\n\t"\
		"movq	%[__cc1]	,%%rdx	\n\t"\
			"movq	%%rax		,%%rbx\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x20		,%%rbx\n\t"\
			"addq	$0x40		,%%rcx\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
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
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
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
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%rbx\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
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
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%rbx\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
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
			"/*****************************************/\n\t"\
			"/* now do three more radix-3 transforms: */\n\t"\
			"/*****************************************/\n\t"\
			"subq	$0xc0		,%%rax\n\t"\
			"subq	$0x80		,%%rbx\n\t"\
			"subq	$0x40		,%%rcx\n\t"\
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
			"subq	$0x180		,%%rdi 		/* __o3 = s1p15r */\n\t"\
			"subq	$0x300		,%%rsi 		/* __o6 = s1p03r */\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rdi)\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%rbx\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"addq	$0x400		,%%rsi 		/* __o7 = s1p35r */\n\t"\
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
			"addq	$0x100		,%%rdi 		/* __o1 = s1p23r */\n\t"\
			"subq	$0x300		,%%rsi 		/* __o4 = s1p11r */\n\t"\
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
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%rbx\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subq	$0x080		,%%rsi 		/* __o5 = s1p07r */\n\t"\
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
			"addq	$0x100		,%%rdi 		/* __o8 = s1p31r */\n\t"\
			"addq	$0x180		,%%rsi 		/* __o2 = s1p19r */\n\t"\
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
			/* SSE2_RADIX_09_DIT_0TWIDDLE(r20,s1p[18,14,10,06,02,34,30,26,22]r): */\
			"movq	%[__out]	,%%rsi\n\t"\
			"addq	$0x240		,%%rsi 		/* s1p18 */\n\t"\
			"movq	%%rsi		,%%rdi\n\t"\
			"movq	%[__r00]	,%%rax\n\t"\
			"addq	$0x240		,%%rax 		/* __r20 */\n\t"\
		"movq	%[__cc1]	,%%rdx	\n\t"\
			"movq	%%rax		,%%rbx\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x20		,%%rbx\n\t"\
			"addq	$0x40		,%%rcx\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
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
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
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
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%rbx\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
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
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%rbx\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
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
			"/*****************************************/\n\t"\
			"/* now do three more radix-3 transforms: */\n\t"\
			"/*****************************************/\n\t"\
			"subq	$0xc0		,%%rax\n\t"\
			"subq	$0x80		,%%rbx\n\t"\
			"subq	$0x40		,%%rcx\n\t"\
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
			"subq	$0x180		,%%rdi 		/* __o3 = s1p06r */\n\t"\
			"addq	$0x180		,%%rsi 		/* __o6 = s1p30r */\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rdi)\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%rbx\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subq	$0x080		,%%rsi 		/* __o7 = s1p26r */\n\t"\
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
			"addq	$0x100		,%%rdi 		/* __o1 = s1p14r */\n\t"\
			"subq	$0x300		,%%rsi 		/* __o4 = s1p02r */\n\t"\
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
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%rbx\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"addq	$0x400		,%%rsi 		/* __o5 = s1p34r */\n\t"\
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
			"addq	$0x100		,%%rdi 		/* __o8 = s1p22r */\n\t"\
			"subq	$0x300		,%%rsi 		/* __o2 = s1p10r */\n\t"\
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
			/* SSE2_RADIX_09_DIT_0TWIDDLE(r30,s1p[09,05,01,33,29,25,21,17,13]r): */\
			"movq	%[__out]	,%%rsi\n\t"\
			"addq	$0x120		,%%rsi 		/* s1p09 */\n\t"\
			"movq	%%rsi		,%%rdi\n\t"\
			"movq	%[__r00]	,%%rax\n\t"\
			"addq	$0x360		,%%rax 		/* __r20 */\n\t"\
		"movq	%[__cc1]	,%%rdx	\n\t"\
			"movq	%%rax		,%%rbx\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x20		,%%rbx\n\t"\
			"addq	$0x40		,%%rcx\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
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
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
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
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%rbx\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
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
			"addq	$0x60		,%%rax\n\t"\
			"addq	$0x60		,%%rbx\n\t"\
			"addq	$0x60		,%%rcx\n\t"\
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
			"/*****************************************/\n\t"\
			"/* now do three more radix-3 transforms: */\n\t"\
			"/*****************************************/\n\t"\
			"subq	$0xc0		,%%rax\n\t"\
			"subq	$0x80		,%%rbx\n\t"\
			"subq	$0x40		,%%rcx\n\t"\
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
			"addq	$0x300		,%%rdi 		/* __o3 = s1p33r */\n\t"\
			"addq	$0x180		,%%rsi 		/* __o6 = s1p21r */\n\t"\
			"addpd	%%xmm5		,%%xmm2\n\t"\
			"subpd	%%xmm4		,%%xmm3\n\t"\
			"subpd	%%xmm5		,%%xmm0\n\t"\
			"addpd	%%xmm4		,%%xmm1\n\t"\
			"movaps	%%xmm2		,    (%%rdi)\n\t"\
			"movaps	%%xmm3		,0x10(%%rdi)\n\t"\
			"movaps	%%xmm0		,    (%%rsi)\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)\n\t"\
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%rbx\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subq	$0x080		,%%rsi 		/* __o7 = s1p17r */\n\t"\
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
			"subq	$0x380		,%%rdi 		/* __o1 = s1p05r */\n\t"\
			"addq	$0x180		,%%rsi 		/* __o4 = s1p29r */\n\t"\
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
			"addq	$0x20		,%%rax\n\t"\
			"addq	$0x20		,%%rbx\n\t"\
			"addq	$0x20		,%%rcx\n\t"\
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm0\n\t"\
			"movaps	%%xmm3		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm2\n\t"\
			"mulpd	%%xmm6		,%%xmm3\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm2\n\t"\
			"subpd	%%xmm0		,%%xmm3\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
			"movaps	    (%%rcx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"addpd	%%xmm1		,%%xmm4\n\t"\
			"subpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"subq	$0x080		,%%rsi 		/* __o5 = s1p25r */\n\t"\
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
			"addq	$0x100		,%%rdi 		/* __o8 = s1p13r */\n\t"\
			"subq	$0x300		,%%rsi 		/* __o2 = s1p01r */\n\t"\
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
			"\n\t"\
			:					/* outputs: none */\
			: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
			 ,[__p01] "m" (Xp01)\
			 ,[__p02] "m" (Xp02)\
			 ,[__p03] "m" (Xp03)\
			 ,[__p04] "m" (Xp04)\
			 ,[__p08] "m" (Xp08)\
			 ,[__p12] "m" (Xp12)\
			 ,[__p16] "m" (Xp16)\
			 ,[__p20] "m" (Xp20)\
			 ,[__p24] "m" (Xp24)\
			 ,[__p28] "m" (Xp28)\
			 ,[__p32] "m" (Xp32)\
			 ,[__r00] "m" (Xr00)\
			 ,[__out] "m" (Xout)\
			 ,[__cc1] "m" (Xcc1)\
			: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
	}


	#define	SSE2_RADIX36_DIF_NOTWIDDLE(Xadd,Xp01,Xp02,Xp03,Xp04,Xp08,Xp12,Xp16,Xp20,Xp24,Xp28,Xp32,Xr00,Xin0,Xcc1)\
	{\
	__asm__ volatile (\
			/* SSE2_RADIX_09_DIF_0TWIDDLE(r00,s1p[00,32,28,24,20,16,12,08,04]r: */\
			"movq	%[__in0]	,%%rax 		/* __i0-8; e[abc]x store input addresses */\n\t"\
			"movq	%[__r00]	,%%rsi\n\t"\
		"movq	%[__cc1]	,%%rdx 		/* edx stores trig addresses throughout */\n\t"\
			"movq	%%rax		,%%rbx\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x300		,%%rbx 		/* __i3 = s1p24r */\n\t"\
			"addq	$0x180		,%%rcx 		/* __i6 = s1p12r */\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
		"movq	%%rsi		,%%rdi\n\t"\
		"addq	$0x20		,%%rdi\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t00 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t01 */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
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
			"addq	$0x400		,%%rax 		/* __i1 = s1p32r */\n\t"\
			"subq	$0x080		,%%rbx 		/* __i4 = s1p20r */\n\t"\
			"subq	$0x080		,%%rcx 		/* __i7 = s1p08r */\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
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
		"addq	$0x40		,%%rsi\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
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
		"addq	$0x40		,%%rsi\n\t"\
			"subq	$0x080		,%%rax 		/* __i2 = s1p28r */\n\t"\
			"subq	$0x080		,%%rbx 		/* __i5 = s1p16r */\n\t"\
			"subq	$0x080		,%%rcx 		/* __i8 = s1p04r */\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
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
		"addq	$0x40		,%%rdi\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
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
			"/*****************************************/\n\t"\
			"/* now do three more radix-3 transforms: */\n\t"\
			"/*****************************************/\n\t"\
			"movq	%[__r00]	,%%rax 		/* __r00 */\n\t"\
			"movq	%%rax		,%%rbx\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x60		,%%rbx 		/* __r06 */\n\t"\
			"addq	$0xc0		,%%rcx 		/* __r0c */\n\t"\
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
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r02 */\n\t"\
			"addq	$0x20		,%%rbx 		/* __r08 */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r0e */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
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
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
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
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r04 */\n\t"\
			"addq	$0x20		,%%rbx 		/* __r0a */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r0g */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
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
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
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
			/* SSE2_RADIX_09_DIF_0TWIDDLE(r10,s1p[27,23,19,15,11,07,03,35,31]r: */\
			"movq	%[__in0]	,%%rax\n\t"\
			"movq	%[__r00]	,%%rsi\n\t"\
			"addq	$0x360		,%%rax 		/* __i0 = s1p27r */\n\t"\
			"addq	$0x120		,%%rsi 		/* r10 */\n\t"\
		"movq	%[__cc1]	,%%rdx 	\n\t"\
			"movq	%%rax		,%%rbx\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"subq	$0x180		,%%rbx 		/* __i3 = s1p15r */\n\t"\
			"subq	$0x300		,%%rcx 		/* __i6 = s1p03r */\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
		"movq	%%rsi		,%%rdi\n\t"\
		"addq	$0x20		,%%rdi\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t00 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t01 */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
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
			"subq	$0x080		,%%rax 		/* __i1 = s1p23r */\n\t"\
			"subq	$0x080		,%%rbx 		/* __i4 = s1p11r */\n\t"\
			"addq	$0x400		,%%rcx 		/* __i7 = s1p35r */\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
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
		"addq	$0x40		,%%rsi\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
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
		"addq	$0x40		,%%rsi\n\t"\
			"subq	$0x080		,%%rax 		/* __i2 = s1p19r */\n\t"\
			"subq	$0x080		,%%rbx 		/* __i5 = s1p07r */\n\t"\
			"subq	$0x080		,%%rcx 		/* __i8 = s1p31r */\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
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
		"addq	$0x40		,%%rdi\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
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
			"/*****************************************/\n\t"\
			"/* now do three more radix-3 transforms: */\n\t"\
			"/*****************************************/\n\t"\
			"movq	%[__r00]	,%%rax 		/* __r00 */\n\t"\
			"addq	$0x120		,%%rax 		/* r10 */\n\t"\
			"movq	%%rax		,%%rbx\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x60		,%%rbx 		/* __r16 */\n\t"\
			"addq	$0xc0		,%%rcx 		/* __r1c */\n\t"\
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
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r12 */\n\t"\
			"addq	$0x20		,%%rbx 		/* __r18 */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r1e */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
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
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
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
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r14 */\n\t"\
			"addq	$0x20		,%%rbx 		/* __r1a */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r1g */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
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
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
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
			/* SSE2_RADIX_09_DIF_0TWIDDLE(r20,s1p[18,14,10,06,02,34,30,26,22]r: */\
			"movq	%[__in0]	,%%rax\n\t"\
			"movq	%[__r00]	,%%rsi\n\t"\
			"addq	$0x240		,%%rax 		/* __i0 = s1p18r */\n\t"\
			"addq	$0x240		,%%rsi 		/* r20 */\n\t"\
		"movq	%[__cc1]	,%%rdx 	\n\t"\
			"movq	%%rax		,%%rbx\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"subq	$0x180		,%%rbx 		/* __i3 = s1p06r */\n\t"\
			"addq	$0x180		,%%rcx 		/* __i6 = s1p30r */\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
		"movq	%%rsi		,%%rdi\n\t"\
		"addq	$0x20		,%%rdi\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t00 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t01 */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
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
			"subq	$0x080		,%%rax 		/* __i1 = s1p14r */\n\t"\
			"subq	$0x080		,%%rbx 		/* __i4 = s1p02r */\n\t"\
			"subq	$0x080		,%%rcx 		/* __i7 = s1p26r */\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
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
		"addq	$0x40		,%%rsi\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
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
		"addq	$0x40		,%%rsi\n\t"\
			"subq	$0x080		,%%rax 		/* __i2 = s1p10r */\n\t"\
			"addq	$0x400		,%%rbx 		/* __i5 = s1p34r */\n\t"\
			"subq	$0x080		,%%rcx 		/* __i8 = s1p22r */\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
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
		"addq	$0x40		,%%rdi\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
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
			"/*****************************************/\n\t"\
			"/* now do three more radix-3 transforms: */\n\t"\
			"/*****************************************/\n\t"\
			"movq	%[__r00]	,%%rax 		/* __r00 */\n\t"\
			"addq	$0x240		,%%rax 		/* r20 */\n\t"\
			"movq	%%rax		,%%rbx\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x60		,%%rbx 		/* __r26 */\n\t"\
			"addq	$0xc0		,%%rcx 		/* __r2c */\n\t"\
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
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r22 */\n\t"\
			"addq	$0x20		,%%rbx 		/* __r28 */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r2e */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
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
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
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
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r24 */\n\t"\
			"addq	$0x20		,%%rbx 		/* __r2a */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r2g */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
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
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
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
			/* SSE2_RADIX_09_DIF_0TWIDDLE(r30,s1p[09,05,01,33,29,25,21,17,13]r: */\
			"movq	%[__in0]	,%%rax\n\t"\
			"movq	%[__r00]	,%%rsi\n\t"\
			"addq	$0x120		,%%rax 		/* __i0 = s1p09r */\n\t"\
			"addq	$0x360		,%%rsi 		/* r30 */\n\t"\
		"movq	%[__cc1]	,%%rdx 	\n\t"\
			"movq	%%rax		,%%rbx\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x300		,%%rbx 		/* __i3 = s1p33r */\n\t"\
			"addq	$0x180		,%%rcx 		/* __i6 = s1p21r */\n\t"\
		"addq	$0x40		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm7\n\t"\
			"movaps	%%xmm2		,%%xmm4\n\t"\
			"movaps	%%xmm3		,%%xmm5\n\t"\
		"movq	%%rsi		,%%rdi\n\t"\
		"addq	$0x20		,%%rdi\n\t"\
			"addpd	%%xmm6		,%%xmm2\n\t"\
			"addpd	%%xmm7		,%%xmm3\n\t"\
			"subpd	%%xmm6		,%%xmm4\n\t"\
			"subpd	%%xmm7		,%%xmm5\n\t"\
			"addpd	%%xmm2		,%%xmm0\n\t"\
			"addpd	%%xmm3		,%%xmm1\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
			"movaps	%%xmm0		,    (%%rsi)	/* <- t00 */\n\t"\
			"movaps	%%xmm1		,0x10(%%rsi)	/* <- t01 */\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
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
			"subq	$0x080		,%%rax 		/* __i1 = s1p05r */\n\t"\
			"subq	$0x080		,%%rbx 		/* __i4 = s1p29r */\n\t"\
			"subq	$0x080		,%%rcx 		/* __i7 = s1p17r */\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
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
		"addq	$0x40		,%%rsi\n\t"\
		"addq	$0x40		,%%rdi\n\t"\
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
		"addq	$0x40		,%%rsi\n\t"\
			"subq	$0x080		,%%rax 		/* __i2 = s1p01r */\n\t"\
			"subq	$0x080		,%%rbx 		/* __i5 = s1p25r */\n\t"\
			"subq	$0x080		,%%rcx 		/* __i8 = s1p13r */\n\t"\
			"movaps	    (%%rbx)	,%%xmm4\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm5\n\t"\
			"movaps	    (%%rax)	,%%xmm0\n\t"\
			"movaps	0x10(%%rax)	,%%xmm1\n\t"\
			"movaps	    (%%rcx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rcx)	,%%xmm3\n\t"\
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
		"addq	$0x40		,%%rdi\n\t"\
		"addq	$0x40		,%%rsi\n\t"\
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
			"/*****************************************/\n\t"\
			"/* now do three more radix-3 transforms: */\n\t"\
			"/*****************************************/\n\t"\
			"movq	%[__r00]	,%%rax 		/* __r00 */\n\t"\
			"addq	$0x360		,%%rax 		/* r30 */\n\t"\
			"movq	%%rax		,%%rbx\n\t"\
			"movq	%%rax		,%%rcx\n\t"\
			"addq	$0x60		,%%rbx 		/* __r36 */\n\t"\
			"addq	$0xc0		,%%rcx 		/* __r3c */\n\t"\
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
		"subq	$0x40		,%%rdx 		/* c1 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r32 */\n\t"\
			"addq	$0x20		,%%rbx 		/* __r38 */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r3e */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c2 */\n\t"\
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
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
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
		"subq	$0x20		,%%rdx 		/* c2 */\n\t"\
			"addq	$0x20		,%%rax 		/* __r34 */\n\t"\
			"addq	$0x20		,%%rbx 		/* __r3a */\n\t"\
			"addq	$0x20		,%%rcx 		/* __r3g */\n\t"\
			"movaps	    (%%rbx)	,%%xmm2\n\t"\
			"movaps	0x10(%%rbx)	,%%xmm3\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"addq	$0x40		,%%rdx 		/* c4 */\n\t"\
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
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
		"subq	$0x20		,%%rdx 		/* c3m1 */\n\t"\
			"movaps	%%xmm4		,%%xmm0\n\t"\
			"movaps	%%xmm5		,%%xmm1\n\t"\
			"mulpd	%%xmm6		,%%xmm4\n\t"\
			"mulpd	%%xmm6		,%%xmm5\n\t"\
			"mulpd	%%xmm7		,%%xmm0\n\t"\
			"mulpd	%%xmm7		,%%xmm1\n\t"\
			"subpd	%%xmm1		,%%xmm4\n\t"\
			"addpd	%%xmm0		,%%xmm5\n\t"\
			"movaps	    (%%rdx)	,%%xmm6\n\t"\
			"movaps	0x10(%%rdx)	,%%xmm7\n\t"\
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
			"\n\t"\
		"/**********************************/"\
		"/*** And now do 9 radix-4 DFTs: ***/"\
		"/**********************************/"\
			"\n\t"\
			/* add0,1,3,2 = a+p00+p[0,1,2,3] */\
			"movq	%[__add],%%rax	/* Must use eax as base address throughout, since that is preserved in SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_B */\n\t"\
			"movq	%[__add],%%rdx						\n\t"\
			"movslq	%[__p01],%%rsi	/* esi will store power-of-2 multiples of p01 throughout */\n\t"\
			"shlq	$3,%%rsi		/* Pointer offset for floating doubles */\n\t"\
			"addq	%%rsi,%%rdx							\n\t"\
			"movq	%%rdx,%%rbx		/* add0+p01 */		\n\t"\
			"addq	%%rsi,%%rdx							\n\t"\
			"movq	%%rdx,%%rcx		/* add0+p02 */		\n\t"\
			"addq	%%rsi,%%rdx		/* add0+p03 */		\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r00, 0x120, 0x240, eax,ebx,edx,ecx) */\n\t"\
			"movq	%[__r00],%%rdi 	/* edi = __r00 */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r00 + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses d <-> c: */		\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%rbx)					\n\t"\
			"movaps	%%xmm4,    (%%rdx)					\n\t"\
			"movaps	%%xmm1,0x10(%%rbx)					\n\t"\
			"movaps	%%xmm5,0x10(%%rcx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rax)					\n\t"\
			"movaps	%%xmm7,    (%%rcx)					\n\t"\
			"movaps	%%xmm3,0x10(%%rax)					\n\t"\
			"movaps	%%xmm6,0x10(%%rdx)					\n\t"\
			"\n\t"\
			/* add0,1,3,2 = a+p32+p[0,1,2,3] */\
			"movslq	%[__p32],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"addq	%%rsi,%%rax		/* a+p32] */	\n\t"\
			"addq	%%rsi,%%rbx							\n\t"\
			"addq	%%rsi,%%rcx							\n\t"\
			"addq	%%rsi,%%rdx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r02, abdc) */\n\t"\
			"subq	$0x220	,%%rdi 	/* edi = __r02 */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r02 + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses d <-> c: */		\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%rbx)					\n\t"\
			"movaps	%%xmm4,    (%%rdx)					\n\t"\
			"movaps	%%xmm1,0x10(%%rbx)					\n\t"\
			"movaps	%%xmm5,0x10(%%rcx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rax)					\n\t"\
			"movaps	%%xmm7,    (%%rcx)					\n\t"\
			"movaps	%%xmm3,0x10(%%rax)					\n\t"\
			"movaps	%%xmm6,0x10(%%rdx)					\n\t"\
			"\n\t"\
			/* add2,3,0,1 = a+p20+p[0,1,2,3] */\
			"movslq	%[__p12],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"subq	%%rsi,%%rax		/* a+p20] */	\n\t"\
			"subq	%%rsi,%%rbx							\n\t"\
			"subq	%%rsi,%%rcx							\n\t"\
			"subq	%%rsi,%%rdx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r04, cdab) */\n\t"\
			"subq	$0x220	,%%rdi 	/* edi = __r04 */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r04 + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses a <-> c ,b <-> d: */\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%rdx)					\n\t"\
			"movaps	%%xmm4,    (%%rax)					\n\t"\
			"movaps	%%xmm1,0x10(%%rdx)					\n\t"\
			"movaps	%%xmm5,0x10(%%rbx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rcx)					\n\t"\
			"movaps	%%xmm7,    (%%rbx)					\n\t"\
			"movaps	%%xmm3,0x10(%%rcx)					\n\t"\
			"movaps	%%xmm6,0x10(%%rax)					\n\t"\
			"\n\t"\
			/* add1,0,2,3 = a+p08+p[0,1,2,3] */\
			"movslq	%[__p12],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"subq	%%rsi,%%rax		/* a+p08] */	\n\t"\
			"subq	%%rsi,%%rbx							\n\t"\
			"subq	%%rsi,%%rcx							\n\t"\
			"subq	%%rsi,%%rdx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r06, bacd) */\n\t"\
			"subq	$0x220	,%%rdi 	/* edi = __r06 */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r06 + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses a <-> b: */		\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%rax)					\n\t"\
			"movaps	%%xmm4,    (%%rcx)					\n\t"\
			"movaps	%%xmm1,0x10(%%rax)					\n\t"\
			"movaps	%%xmm5,0x10(%%rdx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rbx)					\n\t"\
			"movaps	%%xmm7,    (%%rdx)					\n\t"\
			"movaps	%%xmm3,0x10(%%rbx)					\n\t"\
			"movaps	%%xmm6,0x10(%%rcx)					\n\t"\
			"\n\t"\
			/* add3,2,1,0 = a+p28+p[0,1,2,3] */\
			"movslq	%[__p20],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"addq	%%rsi,%%rax		/* a+p28] */	\n\t"\
			"addq	%%rsi,%%rbx							\n\t"\
			"addq	%%rsi,%%rcx							\n\t"\
			"addq	%%rsi,%%rdx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r08, dcba) */\n\t"\
			"subq	$0x220	,%%rdi 	/* edi = __r08 */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r08 + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses a <-> d ,b <-> c: */\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%rcx)					\n\t"\
			"movaps	%%xmm4,    (%%rbx)					\n\t"\
			"movaps	%%xmm1,0x10(%%rcx)					\n\t"\
			"movaps	%%xmm5,0x10(%%rax)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rdx)					\n\t"\
			"movaps	%%xmm7,    (%%rax)					\n\t"\
			"movaps	%%xmm3,0x10(%%rdx)					\n\t"\
			"movaps	%%xmm6,0x10(%%rbx)					\n\t"\
			"\n\t"\
			/* add0,1,3,2 = a+p16+p[0,1,2,3] */\
			"movslq	%[__p12],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"subq	%%rsi,%%rax		/* a+p16] */	\n\t"\
			"subq	%%rsi,%%rbx							\n\t"\
			"subq	%%rsi,%%rcx							\n\t"\
			"subq	%%rsi,%%rdx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0a, abdc) */\n\t"\
			"subq	$0x220	,%%rdi 	/* edi = __r0a */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r0a + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses d <-> c: */		\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%rbx)					\n\t"\
			"movaps	%%xmm4,    (%%rdx)					\n\t"\
			"movaps	%%xmm1,0x10(%%rbx)					\n\t"\
			"movaps	%%xmm5,0x10(%%rcx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rax)					\n\t"\
			"movaps	%%xmm7,    (%%rcx)					\n\t"\
			"movaps	%%xmm3,0x10(%%rax)					\n\t"\
			"movaps	%%xmm6,0x10(%%rdx)					\n\t"\
			"\n\t"\
			/* add2,3,0,1 = a+p04+p[0,1,2,3] */\
			"movslq	%[__p12],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"subq	%%rsi,%%rax		/* a+p04] */	\n\t"\
			"subq	%%rsi,%%rbx							\n\t"\
			"subq	%%rsi,%%rcx							\n\t"\
			"subq	%%rsi,%%rdx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0c, cdab) */\n\t"\
			"subq	$0x220	,%%rdi 	/* edi = __r0c */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r0c + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses a <-> c ,b <-> d: */\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%rdx)					\n\t"\
			"movaps	%%xmm4,    (%%rax)					\n\t"\
			"movaps	%%xmm1,0x10(%%rdx)					\n\t"\
			"movaps	%%xmm5,0x10(%%rbx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rcx)					\n\t"\
			"movaps	%%xmm7,    (%%rbx)					\n\t"\
			"movaps	%%xmm3,0x10(%%rcx)					\n\t"\
			"movaps	%%xmm6,0x10(%%rax)					\n\t"\
			"\n\t"\
			/* add1,0,2,3 = a+p24+p[0,1,2,3] */\
			"movslq	%[__p20],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"addq	%%rsi,%%rax		/* a+p24] */	\n\t"\
			"addq	%%rsi,%%rbx							\n\t"\
			"addq	%%rsi,%%rcx							\n\t"\
			"addq	%%rsi,%%rdx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0e, bacd) */\n\t"\
			"subq	$0x220	,%%rdi 	/* edi = __r0e */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r0e + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses a <-> b: */		\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%rax)					\n\t"\
			"movaps	%%xmm4,    (%%rcx)					\n\t"\
			"movaps	%%xmm1,0x10(%%rax)					\n\t"\
			"movaps	%%xmm5,0x10(%%rdx)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rbx)					\n\t"\
			"movaps	%%xmm7,    (%%rdx)					\n\t"\
			"movaps	%%xmm3,0x10(%%rbx)					\n\t"\
			"movaps	%%xmm6,0x10(%%rcx)					\n\t"\
			"\n\t"\
			/* add,2,1,03 = a+p12+p[0,1,2,3] */\
			"movslq	%[__p12],%%rsi						\n\t"\
			"shlq	$3,%%rsi							\n\t"\
			"subq	%%rsi,%%rax		/* a+p12] */	\n\t"\
			"subq	%%rsi,%%rbx							\n\t"\
			"subq	%%rsi,%%rcx							\n\t"\
			"subq	%%rsi,%%rdx							\n\t"\
			"/* SSE2_RADIX4_DIF_0TWIDDLE_STRIDE_C(r0g, dcba) */\n\t"\
			"subq	$0x220	,%%rdi 	/* edi = __r0g */	\n\t"\
			"movq	$0x120	,%%rsi 						\n\t"\
			"addq	%%rdi	,%%rsi 	/* esi = __r0g + 0x120 */\n\t"\
			"movaps	    (%%rdi),%%xmm4					\n\t"\
			"movaps	    (%%rsi),%%xmm6					\n\t"\
			"movaps	0x10(%%rdi),%%xmm5					\n\t"\
			"movaps	0x10(%%rsi),%%xmm7					\n\t"\
			"addq	$0x240,%%rdi 						\n\t"\
			"addq	$0x240,%%rsi 						\n\t"\
			"movaps	    (%%rdi),%%xmm0					\n\t"\
			"movaps	    (%%rsi),%%xmm2					\n\t"\
			"movaps	0x10(%%rdi),%%xmm1					\n\t"\
			"movaps	0x10(%%rsi),%%xmm3					\n\t"\
			"subpd	%%xmm0,%%xmm4						\n\t"\
			"subpd	%%xmm2,%%xmm6						\n\t"\
			"subpd	%%xmm1,%%xmm5						\n\t"\
			"subpd	%%xmm3,%%xmm7						\n\t"\
			"addpd	%%xmm0,%%xmm0						\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm1,%%xmm1						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm4,%%xmm0						\n\t"\
			"addpd	%%xmm6,%%xmm2						\n\t"\
			"addpd	%%xmm5,%%xmm1						\n\t"\
			"addpd	%%xmm7,%%xmm3						\n\t"\
			"/* Swap output addresses a <-> d ,b <-> c: */\n\t"\
			"subpd	%%xmm2,%%xmm0						\n\t"\
			"subpd	%%xmm7,%%xmm4						\n\t"\
			"subpd	%%xmm3,%%xmm1						\n\t"\
			"subpd	%%xmm6,%%xmm5						\n\t"\
			"movaps	%%xmm0,    (%%rcx)					\n\t"\
			"movaps	%%xmm4,    (%%rbx)					\n\t"\
			"movaps	%%xmm1,0x10(%%rcx)					\n\t"\
			"movaps	%%xmm5,0x10(%%rax)					\n\t"\
			"addpd	%%xmm2,%%xmm2						\n\t"\
			"addpd	%%xmm7,%%xmm7						\n\t"\
			"addpd	%%xmm3,%%xmm3						\n\t"\
			"addpd	%%xmm6,%%xmm6						\n\t"\
			"addpd	%%xmm0,%%xmm2						\n\t"\
			"addpd	%%xmm4,%%xmm7						\n\t"\
			"addpd	%%xmm1,%%xmm3						\n\t"\
			"addpd	%%xmm5,%%xmm6						\n\t"\
			"movaps	%%xmm2,    (%%rdx)					\n\t"\
			"movaps	%%xmm7,    (%%rax)					\n\t"\
			"movaps	%%xmm3,0x10(%%rdx)					\n\t"\
			"movaps	%%xmm6,0x10(%%rbx)					\n\t"\
			"\n\t"\
			:					/* outputs: none */\
			: [__add] "m" (Xadd)	/* All inputs from memory addresses here */\
			 ,[__p01] "m" (Xp01)\
			 ,[__p02] "m" (Xp02)\
			 ,[__p03] "m" (Xp03)\
			 ,[__p04] "m" (Xp04)\
			 ,[__p08] "m" (Xp08)\
			 ,[__p12] "m" (Xp12)\
			 ,[__p16] "m" (Xp16)\
			 ,[__p20] "m" (Xp20)\
			 ,[__p24] "m" (Xp24)\
			 ,[__p28] "m" (Xp28)\
			 ,[__p32] "m" (Xp32)\
			 ,[__r00] "m" (Xr00)\
			 ,[__in0] "m" (Xin0)\
			 ,[__cc1] "m" (Xcc1)\
			: "cc","memory","rax","rbx","rcx","rdx","rdi","rsi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"		/* Clobbered registers */\
		);\
	}

  #endif // USE_64BIT_ASM_STYLE

 #endif	// AVX/SSE2 toggle

#endif	/* radix36_ditN_cy_dif1_gcc_h_included */

