/*******************************************************************************
*                                                                              *
*   (C) 1997-2017 by Ernst W. Mayer.                                           *
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
#ifndef radix16_utils_asm_h_included
#define radix16_utils_asm_h_included

#ifdef USE_AVX512	// AVX512 implements a 512-bit-register version of the the AVX2 ALL_FMA-macro

	// SIMD-gen 15 nontrivial twiddles; 1,2,4,8,13 via 2-table-mul, remaining 10 in 5 pair-cmuls from those:
	#define SSE2_RADIX16_CALC_TWIDDLES_LOACC(Xcc0,Xk0,Xk1,Xrt0,Xrt1)\
	{\
		__asm__ volatile (\
		/* Note bytewidth-multiplier = 1 in the vector-addressing since already have incorporated into k1,2_arr data */\
		"movq	%[__k0]	,%%rax	\n\t"\
		"movq	%[__k1]	,%%rbx	\n\t"\
		"movq	%[__rt0],%%rcx	\n\t"\
		"movq	%[__rt1],%%rdx	\n\t"\
		"movq	%[__cc0],%%r8	\n\t"\
	/* In AVX-512 version, do 8x64-bit gather-load of zmm based on 8-element chunks of k0,1_arr[]: */\
		/* Mask-reg zmm9 = 11...11 - the gather-load opmask-reg is stupidly zeroed each time we do gather-load, so need to reinit: */\
		"movl	$-1,%%esi	\n\t"/* Use to init opmask k1 (Only need the low byte) */\
	/* c1,s1: */\
		"vmovaps	0x00(%%rax),%%ymm4		\n\t"/* k0_arr[0- 7] */\
		"vmovaps	0x00(%%rbx),%%ymm5		\n\t"/* k1_arr[0- 7] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x0(%%rcx,%%ymm4),%%zmm0%{%%k1%}	\n\t"/* m0 = [re0.A-H] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x8(%%rcx,%%ymm4),%%zmm1%{%%k1%}	\n\t"/* m1 = [im0.A-H] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x0(%%rdx,%%ymm5),%%zmm2%{%%k1%}	\n\t"/* m2 = [re1.A-H] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x8(%%rdx,%%ymm5),%%zmm3%{%%k1%}	\n\t"/* m3 = [im1.A-H] */\
		"vmovaps	%%zmm0,%%zmm4			\n\t"/* cpy re0 */\
		"vmovaps	%%zmm1,%%zmm5			\n\t"/* cpy im0 */\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t"/* re1*re0, overwrites re0 */\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t"/* im1*im0, overwrites im0 */\
		"vmulpd	%%zmm5,%%zmm2,%%zmm2	\n\t"/* im0*re1, overwrites re1 */\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t"/* re0*im1, overwrites im1 */\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t"/* Re */\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t"/* Im */\
		"vmovaps	%%zmm0,0x480(%%r8)		\n\t"/* cc0 + 0x12 */\
		"vmovaps	%%zmm2,0x4c0(%%r8)		\n\t"\
		"vmovaps	%%zmm0,%%zmm6 	\n\t	vmovaps	%%zmm2,%%zmm7 	\n\t"/* Stash copies of c1,s1 in m6,7 */\
	/* c2,s2: */\
		"vmovaps	0x20(%%rax),%%ymm4		\n\t"/* k0_arr[ 8-15] */\
		"vmovaps	0x20(%%rbx),%%ymm5		\n\t"/* k1_arr[ 8-15] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x0(%%rcx,%%ymm4),%%zmm0%{%%k1%}	\n\t"/* m0 = [re0.A-H] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x8(%%rcx,%%ymm4),%%zmm1%{%%k1%}	\n\t"/* m1 = [im0.A-H] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x0(%%rdx,%%ymm5),%%zmm2%{%%k1%}	\n\t"/* m2 = [re1.A-H] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x8(%%rdx,%%ymm5),%%zmm3%{%%k1%}	\n\t"/* m3 = [im1.A-H] */\
		"vmovaps	%%zmm0,%%zmm4			\n\t"\
		"vmovaps	%%zmm1,%%zmm5			\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t"\
		"vmulpd	%%zmm5,%%zmm2,%%zmm2	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t"\
		"vmovaps	%%zmm0,0x280(%%r8)		\n\t"/* cc0 + 0xa */\
		"vmovaps	%%zmm2,0x2c0(%%r8)		\n\t"\
		"vmovaps	%%zmm0,%%zmm8 	\n\t	vmovaps	%%zmm2,%%zmm9 	\n\t"/* stash c2,s2 in m8,9 */\
	/* c4,s4: */\
		"vmovaps	0x40(%%rax),%%ymm4		\n\t"/* k0_arr[16-23] */\
		"vmovaps	0x40(%%rbx),%%ymm5		\n\t"/* k1_arr[16-23] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x0(%%rcx,%%ymm4),%%zmm0%{%%k1%}	\n\t"/* m0 = [re0.A-H] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x8(%%rcx,%%ymm4),%%zmm1%{%%k1%}	\n\t"/* m1 = [im0.A-H] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x0(%%rdx,%%ymm5),%%zmm2%{%%k1%}	\n\t"/* m2 = [re1.A-H] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x8(%%rdx,%%ymm5),%%zmm3%{%%k1%}	\n\t"/* m3 = [im1.A-H] */\
		"vmovaps	%%zmm0,%%zmm4			\n\t"\
		"vmovaps	%%zmm1,%%zmm5			\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t"\
		"vmulpd	%%zmm5,%%zmm2,%%zmm2	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t"\
		"vmovaps	%%zmm0,0x180(%%r8)		\n\t"/* cc0 + 0x6 */\
		"vmovaps	%%zmm2,0x1c0(%%r8)		\n\t"\
		"vmovaps	%%zmm0,%%zmm10	\n\t	vmovaps	%%zmm2,%%zmm11	\n\t"/* stash c4,s4 in m10,11 */\
	/* c8,s8: */\
		"vmovaps	0x60(%%rax),%%ymm4		\n\t"/* k0_arr[24-31] */\
		"vmovaps	0x60(%%rbx),%%ymm5		\n\t"/* k1_arr[24-31] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x0(%%rcx,%%ymm4),%%zmm0%{%%k1%}	\n\t"/* m0 = [re0.A-H] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x8(%%rcx,%%ymm4),%%zmm1%{%%k1%}	\n\t"/* m1 = [im0.A-H] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x0(%%rdx,%%ymm5),%%zmm2%{%%k1%}	\n\t"/* m2 = [re1.A-H] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x8(%%rdx,%%ymm5),%%zmm3%{%%k1%}	\n\t"/* m3 = [im1.A-H] */\
		"vmovaps	%%zmm0,%%zmm4			\n\t"\
		"vmovaps	%%zmm1,%%zmm5			\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t"\
		"vmulpd	%%zmm5,%%zmm2,%%zmm2	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t"\
		"vmovaps	%%zmm0,0x100(%%r8)		\n\t"/* cc0 + 0x4 */\
		"vmovaps	%%zmm2,0x140(%%r8)		\n\t"\
		"vmovaps	%%zmm0,%%zmm12	\n\t	vmovaps	%%zmm2,%%zmm13	\n\t"/* stash c8,s8 in m12,13 */\
	/* c13,s13: */\
		"vmovaps	0x80(%%rax),%%ymm4		\n\t"/* k0_arr[32-39] */\
		"vmovaps	0x80(%%rbx),%%ymm5		\n\t"/* k1_arr[32-39] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x0(%%rcx,%%ymm4),%%zmm0%{%%k1%}	\n\t"/* m0 = [re0.A-H] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x8(%%rcx,%%ymm4),%%zmm1%{%%k1%}	\n\t"/* m1 = [im0.A-H] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x0(%%rdx,%%ymm5),%%zmm2%{%%k1%}	\n\t"/* m2 = [re1.A-H] */\
		"kmovw	%%esi,%%k1	\n\t	vgatherdpd 0x8(%%rdx,%%ymm5),%%zmm3%{%%k1%}	\n\t"/* m3 = [im1.A-H] */\
		"vmovaps	%%zmm0,%%zmm4			\n\t"\
		"vmovaps	%%zmm1,%%zmm5			\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t"\
		"vmulpd	%%zmm5,%%zmm2,%%zmm2	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t"\
		"vmovaps	%%zmm0,0x600(%%r8)		\n\t"/* cc0 + 0x18 */\
		"vmovaps	%%zmm2,0x640(%%r8)		\n\t"\
		"vmovaps	%%zmm0,%%zmm14	\n\t	vmovaps	%%zmm2,%%zmm15	\n\t"/* stash c13,s13 in m14,15 */\
	/* Mem-mapping: (c,s)[0-15] :  0  1 2  3 4  5 6  7 8  9 a  b c  d  e  f  */\
	/*				(cc0,ss0) + 0x[2,12,a,1a,6,16,e,1e,4,14,c,1c,8,18,10,20] */\
	/* SSE2_CMUL_EXPO(c1,c4 ,c3 ,c5 ): */\
		"vmovaps	%%zmm6 ,%%zmm0	\n\t	vmovaps	%%zmm7 ,%%zmm2\n\t"/* c1,s1 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c1,s1 */\
		"vmovaps	%%zmm10,%%zmm4	\n\t	vmovaps	%%zmm11,%%zmm5\n\t"/* c4,s4 */\
		"vmulpd	%%zmm4,%%zmm0,%%zmm0	\n\t"\
		"vmulpd	%%zmm5,%%zmm1,%%zmm1	\n\t"\
		"vmulpd	%%zmm4,%%zmm2,%%zmm2	\n\t"\
		"vmulpd	%%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmovaps	%%zmm0,%%zmm4	\n\t"\
		"vmovaps	%%zmm1,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm0,%%zmm0	\n\t"\
		"vsubpd	%%zmm2,%%zmm1,%%zmm1	\n\t"\
		"vsubpd	%%zmm3,%%zmm4,%%zmm4	\n\t"\
		"vaddpd	%%zmm2,%%zmm5,%%zmm5	\n\t"\
		"vmovaps	%%zmm0,0x680(%%r8)	\n\t	vmovaps	%%zmm1,0x6c0(%%r8)	\n\t"/* c3 = cc0 + 0x1a */\
		"vmovaps	%%zmm4,0x580(%%r8)	\n\t	vmovaps	%%zmm5,0x5c0(%%r8)	\n\t"/* c5 = cc0 + 0x16 */\
	/* SSE2_CMUL_EXPO(c1,c8 ,c7 ,c9 ): */\
		"vmovaps	%%zmm6 ,%%zmm0	\n\t	vmovaps	%%zmm7 ,%%zmm2\n\t"/* c1,s1 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c1,s1 */\
		"vmovaps	%%zmm12,%%zmm4	\n\t	vmovaps	%%zmm13,%%zmm5\n\t"/* c8,s8 */\
		"vmulpd	%%zmm4,%%zmm0,%%zmm0	\n\t"\
		"vmulpd	%%zmm5,%%zmm1,%%zmm1	\n\t"\
		"vmulpd	%%zmm4,%%zmm2,%%zmm2	\n\t"\
		"vmulpd	%%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmovaps	%%zmm0,%%zmm4	\n\t"\
		"vmovaps	%%zmm1,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm0,%%zmm0	\n\t"\
		"vsubpd	%%zmm2,%%zmm1,%%zmm1	\n\t"\
		"vsubpd	%%zmm3,%%zmm4,%%zmm4	\n\t"\
		"vaddpd	%%zmm2,%%zmm5,%%zmm5	\n\t"\
		"vmovaps	%%zmm0,0x780(%%r8)	\n\t	vmovaps	%%zmm1,0x7c0(%%r8)	\n\t"/* c7 = cc0 + 0x1e */\
		"vmovaps	%%zmm4,0x500(%%r8)	\n\t	vmovaps	%%zmm5,0x540(%%r8)	\n\t"/* c9 = cc0 + 0x14 */\
	/* SSE2_CMUL_EXPO(c2,c8 ,c6 ,c10): */\
		"vmovaps	%%zmm8 ,%%zmm0	\n\t	vmovaps	%%zmm9 ,%%zmm2\n\t"/* c2,s2 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c2,s2 */\
		"vmovaps	%%zmm12,%%zmm4	\n\t	vmovaps	%%zmm13,%%zmm5\n\t"/* c8,s8 */\
		"vmulpd	%%zmm4,%%zmm0,%%zmm0	\n\t"\
		"vmulpd	%%zmm5,%%zmm1,%%zmm1	\n\t"\
		"vmulpd	%%zmm4,%%zmm2,%%zmm2	\n\t"\
		"vmulpd	%%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmovaps	%%zmm0,%%zmm4	\n\t"\
		"vmovaps	%%zmm1,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm0,%%zmm0	\n\t"\
		"vsubpd	%%zmm2,%%zmm1,%%zmm1	\n\t"\
		"vsubpd	%%zmm3,%%zmm4,%%zmm4	\n\t"\
		"vaddpd	%%zmm2,%%zmm5,%%zmm5	\n\t"\
		"vmovaps	%%zmm0,0x380(%%r8)	\n\t	vmovaps	%%zmm1,0x3c0(%%r8)	\n\t"/* c6  = cc0 + 0xe */\
		"vmovaps	%%zmm4,0x300(%%r8)	\n\t	vmovaps	%%zmm5,0x340(%%r8)	\n\t"/* c10 = cc0 + 0xc */\
	/* SSE2_CMUL_EXPO(c1,c13,c12,c14): */\
		"vmovaps	%%zmm6 ,%%zmm0	\n\t	vmovaps	%%zmm7 ,%%zmm2\n\t"/* c1,s1 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c1,s1 */\
		"vmovaps	%%zmm14,%%zmm4	\n\t	vmovaps	%%zmm15,%%zmm5\n\t"/* c13,s13 */\
		"vmulpd	%%zmm4,%%zmm0,%%zmm0	\n\t"\
		"vmulpd	%%zmm5,%%zmm1,%%zmm1	\n\t"\
		"vmulpd	%%zmm4,%%zmm2,%%zmm2	\n\t"\
		"vmulpd	%%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmovaps	%%zmm0,%%zmm4	\n\t"\
		"vmovaps	%%zmm1,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm0,%%zmm0	\n\t"\
		"vsubpd	%%zmm2,%%zmm1,%%zmm1	\n\t"\
		"vsubpd	%%zmm3,%%zmm4,%%zmm4	\n\t"\
		"vaddpd	%%zmm2,%%zmm5,%%zmm5	\n\t"\
		"vmovaps	%%zmm0,0x200(%%r8)	\n\t	vmovaps	%%zmm1,0x240(%%r8)	\n\t"/* c12 = cc0 + 0x08 */\
		"vmovaps	%%zmm4,0x400(%%r8)	\n\t	vmovaps	%%zmm5,0x440(%%r8)	\n\t"/* c14 = cc0 + 0x10 */\
	/* SSE2_CMUL_EXPO(c2,c13,c11,c15): */\
		"vmovaps	%%zmm8 ,%%zmm0	\n\t	vmovaps	%%zmm9 ,%%zmm2\n\t"/* c2,s2 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c2,s2 */\
		"vmovaps	%%zmm14,%%zmm4	\n\t	vmovaps	%%zmm15,%%zmm5\n\t"/* c13,s13 */\
		"vmulpd	%%zmm4,%%zmm0,%%zmm0	\n\t"\
		"vmulpd	%%zmm5,%%zmm1,%%zmm1	\n\t"\
		"vmulpd	%%zmm4,%%zmm2,%%zmm2	\n\t"\
		"vmulpd	%%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmovaps	%%zmm0,%%zmm4	\n\t"\
		"vmovaps	%%zmm1,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm0,%%zmm0	\n\t"\
		"vsubpd	%%zmm2,%%zmm1,%%zmm1	\n\t"\
		"vsubpd	%%zmm3,%%zmm4,%%zmm4	\n\t"\
		"vaddpd	%%zmm2,%%zmm5,%%zmm5	\n\t"\
		"vmovaps	%%zmm0,0x700(%%r8)	\n\t	vmovaps	%%zmm1,0x740(%%r8)	\n\t"/* c11 = cc0 + 0x1c */\
		"vmovaps	%%zmm4,0x800(%%r8)	\n\t	vmovaps	%%zmm5,0x840(%%r8)	\n\t"/* c15 = cc0 + 0x20 */\
		:					/* outputs: none */\
		: [__cc0] "m" (Xcc0)	/* All inputs from memory addresses here */\
		 ,[__k0]  "m" (Xk0)\
		 ,[__k1]  "m" (Xk1)\
		 ,[__rt0] "m" (Xrt0)\
		 ,[__rt1] "m" (Xrt1)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","r8","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
	}

#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers

	// SIMD-gen 15 nontrivial twiddles; 1,2,4,8,13 via 2-table-mul, remaining 10 in 5 pair-cmuls from those:
	#define SSE2_RADIX16_CALC_TWIDDLES_LOACC(Xcc0,Xk0,Xk1,Xrt0,Xrt1)\
	{\
		__asm__ volatile (\
		"movq	%[__k0]	,%%rax\n\t"\
		"movq	%[__k1]	,%%rbx\n\t"\
		"movq	%[__rt0],%%rcx\n\t"\
		"movq	%[__rt1],%%rdx\n\t"\
		"movq	%[__cc0],%%r10\n\t"\
		/* In AVX version, do 2-part loads of lo/hi ymm-halves - lo-half using k0,1_arr[j], hi-half using k0,1_arr[j+10]: */\
		/* c1,s1: */\
		"movl	0x0(%%rax),%%esi		\n\t"/* k0_arr[0] */\
		"movl	0x4(%%rax),%%edi		\n\t"/* k0_arr[1] */\
		"vmovaps	(%%rcx,%%rsi),%%xmm0	\n\t"/* [re0,im0].A; note 128-bit lo-half loads here! */\
		"vmovaps	(%%rcx,%%rdi),%%xmm2	\n\t"/* [re0,im0].B */\
		"movl	0x8(%%rax),%%esi		\n\t"/* k0_arr[2] */\
		"movl	0xc(%%rax),%%edi		\n\t"/* k0_arr[3] */\
	"vinsertf128 $1,(%%rcx,%%rsi),%%ymm0,%%ymm0	\n\t"/* [re0,im0].C into hi-half of m0 */\
	"vinsertf128 $1,(%%rcx,%%rdi),%%ymm2,%%ymm2	\n\t"/* [re0,im0].D into hi-half of m2 */\
		"vmovaps		%%ymm0 ,%%ymm1		\n\t"/* cpy ymm0 */\
		"vshufpd $0x0,%%ymm2,%%ymm0,%%ymm0	\n\t"/* m0 = [re0.A,re0.B,re0.C,re0.D] */\
		"vshufpd $0xf,%%ymm2,%%ymm1,%%ymm1	\n\t"/* m1 = [im0.A,im0.B,im0.C,im0.D] */\
		"movl	0x0(%%rbx),%%esi		\n\t"/* k1_arr[0] */\
		"movl	0x4(%%rbx),%%edi		\n\t"/* k1_arr[1] */\
		"vmovaps	(%%rdx,%%rsi),%%xmm2	\n\t"/* [re1,im1].A */\
		"vmovaps	(%%rdx,%%rdi),%%xmm4	\n\t"/* [re1,im1].B */\
		"movl	0x8(%%rbx),%%esi		\n\t"/* k1_arr[2] */\
		"movl	0xc(%%rbx),%%edi		\n\t"/* k1_arr[3] */\
	"vinsertf128 $1,(%%rdx,%%rsi),%%ymm2,%%ymm2	\n\t"/* [re1,im1].C into hi-half of m2 */\
	"vinsertf128 $1,(%%rdx,%%rdi),%%ymm4,%%ymm4	\n\t"/* [re1,im1].D into hi-half of m4 */\
		"vmovaps		%%ymm2 ,%%ymm3		\n\t"/* cpy ymm2 */\
		"vshufpd $0x0,%%ymm4,%%ymm2,%%ymm2	\n\t"/* m2 = [re1.A,re1.B,re1.C,re1.D] */\
		"vshufpd $0xf,%%ymm4,%%ymm3,%%ymm3	\n\t"/* m3 = [im1.A,im1.B,im1.C,im1.D] */\
		"vmovaps	%%ymm0,%%ymm4			\n\t"/* cpy re0 */\
		"vmovaps	%%ymm1,%%ymm5			\n\t"/* cpy im0 */\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t"/* re1*re0, overwrites re0 */\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t"/* im1*im0, overwrites im0 */\
		"vmulpd	%%ymm5,%%ymm2,%%ymm2	\n\t"/* im0*re1, overwrites re1 */\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t"/* re0*im1, overwrites im1 */\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t"/* Re */\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t"/* Im */\
		"vmovaps	%%ymm0,0x240(%%r10)		\n\t"/* cc0 + 0x12 */\
		"vmovaps	%%ymm2,0x260(%%r10)		\n\t"\
		"vmovaps	%%ymm0,%%ymm6 	\n\t	vmovaps	%%ymm2,%%ymm7 	\n\t"/* Register copies are free; stash c1,s1 in m6,7 */\
		/* c2,s2: */\
		"movl	0x10(%%rax),%%esi		\n\t"/* k0_arr[4] */\
		"movl	0x14(%%rax),%%edi		\n\t"/* k0_arr[5] */\
		"vmovaps	(%%rcx,%%rsi),%%xmm0	\n\t"\
		"vmovaps	(%%rcx,%%rdi),%%xmm2	\n\t"\
		"movl	0x18(%%rax),%%esi		\n\t"/* k0_arr[6] */\
		"movl	0x1c(%%rax),%%edi		\n\t"/* k0_arr[7] */\
	"vinsertf128 $1,(%%rcx,%%rsi),%%ymm0,%%ymm0	\n\t"\
	"vinsertf128 $1,(%%rcx,%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vmovaps		%%ymm0 ,%%ymm1		\n\t"\
		"vshufpd $0x0,%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vshufpd $0xf,%%ymm2,%%ymm1,%%ymm1	\n\t"\
		"movl	0x10(%%rbx),%%esi		\n\t"/* k1_arr[4] */\
		"movl	0x14(%%rbx),%%edi		\n\t"/* k1_arr[5] */\
		"vmovaps	(%%rdx,%%rsi),%%xmm2	\n\t"\
		"vmovaps	(%%rdx,%%rdi),%%xmm4	\n\t"\
		"movl	0x18(%%rbx),%%esi		\n\t"/* k1_arr[6] */\
		"movl	0x1c(%%rbx),%%edi		\n\t"/* k1_arr[7] */\
	"vinsertf128 $1,(%%rdx,%%rsi),%%ymm2,%%ymm2	\n\t"\
	"vinsertf128 $1,(%%rdx,%%rdi),%%ymm4,%%ymm4	\n\t"\
		"vmovaps		%%ymm2 ,%%ymm3		\n\t"\
		"vshufpd $0x0,%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vshufpd $0xf,%%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm0,%%ymm4			\n\t"\
		"vmovaps	%%ymm1,%%ymm5			\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm5,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t"\
		"vmovaps	%%ymm0,0x140(%%r10)		\n\t"/* cc0 + 0xa */\
		"vmovaps	%%ymm2,0x160(%%r10)		\n\t"\
		"vmovaps	%%ymm0,%%ymm8 	\n\t	vmovaps	%%ymm2,%%ymm9 	\n\t"/* stash c2,s2 in m8,9 */\
		/* c4,s4: */\
		"movl	0x20(%%rax),%%esi		\n\t"/* k0_arr[8] */\
		"movl	0x24(%%rax),%%edi		\n\t"/* k0_arr[9] */\
		"vmovaps	(%%rcx,%%rsi),%%xmm0	\n\t"\
		"vmovaps	(%%rcx,%%rdi),%%xmm2	\n\t"\
		"movl	0x28(%%rax),%%esi		\n\t"/* k0_arr[10] */\
		"movl	0x2c(%%rax),%%edi		\n\t"/* k0_arr[11] */\
	"vinsertf128 $1,(%%rcx,%%rsi),%%ymm0,%%ymm0	\n\t"\
	"vinsertf128 $1,(%%rcx,%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vmovaps		%%ymm0 ,%%ymm1		\n\t"\
		"vshufpd $0x0,%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vshufpd $0xf,%%ymm2,%%ymm1,%%ymm1	\n\t"\
		"movl	0x20(%%rbx),%%esi		\n\t"/* k1_arr[8] */\
		"movl	0x24(%%rbx),%%edi		\n\t"/* k1_arr[9] */\
		"vmovaps	(%%rdx,%%rsi),%%xmm2	\n\t"\
		"vmovaps	(%%rdx,%%rdi),%%xmm4	\n\t"\
		"movl	0x28(%%rbx),%%esi		\n\t"/* k1_arr[10] */\
		"movl	0x2c(%%rbx),%%edi		\n\t"/* k1_arr[11] */\
	"vinsertf128 $1,(%%rdx,%%rsi),%%ymm2,%%ymm2	\n\t"\
	"vinsertf128 $1,(%%rdx,%%rdi),%%ymm4,%%ymm4	\n\t"\
		"vmovaps		%%ymm2 ,%%ymm3		\n\t"\
		"vshufpd $0x0,%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vshufpd $0xf,%%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm0,%%ymm4			\n\t"\
		"vmovaps	%%ymm1,%%ymm5			\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm5,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t"\
		"vmovaps	%%ymm0,0x0c0(%%r10)		\n\t"/* cc0 + 0x6 */\
		"vmovaps	%%ymm2,0x0e0(%%r10)		\n\t"\
		"vmovaps	%%ymm0,%%ymm10	\n\t	vmovaps	%%ymm2,%%ymm11	\n\t"/* stash c4,s4 in m10,11 */\
		/* c8,s8: */\
		"movl	0x30(%%rax),%%esi		\n\t"/* k0_arr[12] */\
		"movl	0x34(%%rax),%%edi		\n\t"/* k0_arr[13] */\
		"vmovaps	(%%rcx,%%rsi),%%xmm0	\n\t"\
		"vmovaps	(%%rcx,%%rdi),%%xmm2	\n\t"\
		"movl	0x38(%%rax),%%esi		\n\t"/* k0_arr[14] */\
		"movl	0x3c(%%rax),%%edi		\n\t"/* k0_arr[15] */\
	"vinsertf128 $1,(%%rcx,%%rsi),%%ymm0,%%ymm0	\n\t"\
	"vinsertf128 $1,(%%rcx,%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vmovaps		%%ymm0 ,%%ymm1		\n\t"\
		"vshufpd $0x0,%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vshufpd $0xf,%%ymm2,%%ymm1,%%ymm1	\n\t"\
		"movl	0x30(%%rbx),%%esi		\n\t"/* k1_arr[12] */\
		"movl	0x34(%%rbx),%%edi		\n\t"/* k1_arr[13] */\
		"vmovaps	(%%rdx,%%rsi),%%xmm2	\n\t"\
		"vmovaps	(%%rdx,%%rdi),%%xmm4	\n\t"\
		"movl	0x38(%%rbx),%%esi		\n\t"/* k1_arr[14] */\
		"movl	0x3c(%%rbx),%%edi		\n\t"/* k1_arr[15] */\
	"vinsertf128 $1,(%%rdx,%%rsi),%%ymm2,%%ymm2	\n\t"\
	"vinsertf128 $1,(%%rdx,%%rdi),%%ymm4,%%ymm4	\n\t"\
		"vmovaps		%%ymm2 ,%%ymm3		\n\t"\
		"vshufpd $0x0,%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vshufpd $0xf,%%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm0,%%ymm4			\n\t"\
		"vmovaps	%%ymm1,%%ymm5			\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm5,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t"\
		"vmovaps	%%ymm0,0x080(%%r10)		\n\t"/* cc0 + 0x4 */\
		"vmovaps	%%ymm2,0x0a0(%%r10)		\n\t"\
		"vmovaps	%%ymm0,%%ymm12	\n\t	vmovaps	%%ymm2,%%ymm13	\n\t"/* stash c8,s8 in m12,13 */\
		/* c13,s13: */\
		"movl	0x40(%%rax),%%esi		\n\t"/* k0_arr[16] */\
		"movl	0x44(%%rax),%%edi		\n\t"/* k0_arr[17] */\
		"vmovaps	(%%rcx,%%rsi),%%xmm0	\n\t"\
		"vmovaps	(%%rcx,%%rdi),%%xmm2	\n\t"\
		"movl	0x48(%%rax),%%esi		\n\t"/* k0_arr[18] */\
		"movl	0x4c(%%rax),%%edi		\n\t"/* k0_arr[19] */\
	"vinsertf128 $1,(%%rcx,%%rsi),%%ymm0,%%ymm0	\n\t"\
	"vinsertf128 $1,(%%rcx,%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vmovaps		%%ymm0 ,%%ymm1		\n\t"\
		"vshufpd $0x0,%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vshufpd $0xf,%%ymm2,%%ymm1,%%ymm1	\n\t"\
		"movl	0x40(%%rbx),%%esi		\n\t"/* k1_arr[16] */\
		"movl	0x44(%%rbx),%%edi		\n\t"/* k1_arr[17] */\
		"vmovaps	(%%rdx,%%rsi),%%xmm2	\n\t"\
		"vmovaps	(%%rdx,%%rdi),%%xmm4	\n\t"\
		"movl	0x48(%%rbx),%%esi		\n\t"/* k1_arr[18] */\
		"movl	0x4c(%%rbx),%%edi		\n\t"/* k1_arr[19] */\
	"vinsertf128 $1,(%%rdx,%%rsi),%%ymm2,%%ymm2	\n\t"\
	"vinsertf128 $1,(%%rdx,%%rdi),%%ymm4,%%ymm4	\n\t"\
		"vmovaps		%%ymm2 ,%%ymm3		\n\t"\
		"vshufpd $0x0,%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vshufpd $0xf,%%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm0,%%ymm4			\n\t"\
		"vmovaps	%%ymm1,%%ymm5			\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm5,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t"\
		"vmovaps	%%ymm0,0x300(%%r10)		\n\t"/* cc0 + 0x18 */\
		"vmovaps	%%ymm2,0x320(%%r10)		\n\t"\
		"vmovaps	%%ymm0,%%ymm14	\n\t	vmovaps	%%ymm2,%%ymm15	\n\t"/* stash c13,s13 in m14,15 */\
	/* Mem-mapping: (c,s)[0-15] :  0  1 2  3 4  5 6  7 8  9 a  b c  d  e  f  */\
	/*				(cc0,ss0) + 0x[2,12,a,1a,6,16,e,1e,4,14,c,1c,8,18,10,20] */\
	/* SSE2_CMUL_EXPO(c1,c4 ,c3 ,c5 ): */\
		"vmovaps	%%ymm6 ,%%ymm0	\n\t	vmovaps	%%ymm7 ,%%ymm2\n\t"/* c1,s1 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c1,s1 */\
		"vmovaps	%%ymm10,%%ymm4	\n\t	vmovaps	%%ymm11,%%ymm5\n\t"/* c4,s4 */\
		"vmulpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm0,%%ymm4	\n\t"\
		"vmovaps	%%ymm1,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm0,0x340(%%r10)	\n\t	vmovaps	%%ymm1,0x360(%%r10)	\n\t"/* c3 = cc0 + 0x1a */\
		"vmovaps	%%ymm4,0x2c0(%%r10)	\n\t	vmovaps	%%ymm5,0x2e0(%%r10)	\n\t"/* c5 = cc0 + 0x16 */\
	/* SSE2_CMUL_EXPO(c1,c8 ,c7 ,c9 ): */\
		"vmovaps	%%ymm6 ,%%ymm0	\n\t	vmovaps	%%ymm7 ,%%ymm2\n\t"/* c1,s1 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c1,s1 */\
		"vmovaps	%%ymm12,%%ymm4	\n\t	vmovaps	%%ymm13,%%ymm5\n\t"/* c8,s8 */\
		"vmulpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm0,%%ymm4	\n\t"\
		"vmovaps	%%ymm1,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm0,0x3c0(%%r10)	\n\t	vmovaps	%%ymm1,0x3e0(%%r10)	\n\t"/* c7 = cc0 + 0x1e */\
		"vmovaps	%%ymm4,0x280(%%r10)	\n\t	vmovaps	%%ymm5,0x2a0(%%r10)	\n\t"/* c9 = cc0 + 0x14 */\
	/* SSE2_CMUL_EXPO(c2,c8 ,c6 ,c10): */\
		"vmovaps	%%ymm8 ,%%ymm0	\n\t	vmovaps	%%ymm9 ,%%ymm2\n\t"/* c2,s2 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c2,s2 */\
		"vmovaps	%%ymm12,%%ymm4	\n\t	vmovaps	%%ymm13,%%ymm5\n\t"/* c8,s8 */\
		"vmulpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm0,%%ymm4	\n\t"\
		"vmovaps	%%ymm1,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm0,0x1c0(%%r10)	\n\t	vmovaps	%%ymm1,0x1e0(%%r10)	\n\t"/* c6  = cc0 + 0xe */\
		"vmovaps	%%ymm4,0x180(%%r10)	\n\t	vmovaps	%%ymm5,0x1a0(%%r10)	\n\t"/* c10 = cc0 + 0xc */\
	/* SSE2_CMUL_EXPO(c1,c13,c12,c14): */\
		"vmovaps	%%ymm6 ,%%ymm0	\n\t	vmovaps	%%ymm7 ,%%ymm2\n\t"/* c1,s1 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c1,s1 */\
		"vmovaps	%%ymm14,%%ymm4	\n\t	vmovaps	%%ymm15,%%ymm5\n\t"/* c13,s13 */\
		"vmulpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm0,%%ymm4	\n\t"\
		"vmovaps	%%ymm1,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm0,0x100(%%r10)	\n\t	vmovaps	%%ymm1,0x120(%%r10)	\n\t"/* c12 = cc0 + 0x08 */\
		"vmovaps	%%ymm4,0x200(%%r10)	\n\t	vmovaps	%%ymm5,0x220(%%r10)	\n\t"/* c14 = cc0 + 0x10 */\
	/* SSE2_CMUL_EXPO(c2,c13,c11,c15): */\
		"vmovaps	%%ymm8 ,%%ymm0	\n\t	vmovaps	%%ymm9 ,%%ymm2\n\t"/* c2,s2 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c2,s2 */\
		"vmovaps	%%ymm14,%%ymm4	\n\t	vmovaps	%%ymm15,%%ymm5\n\t"/* c13,s13 */\
		"vmulpd	%%ymm4,%%ymm0,%%ymm0	\n\t"\
		"vmulpd	%%ymm5,%%ymm1,%%ymm1	\n\t"\
		"vmulpd	%%ymm4,%%ymm2,%%ymm2	\n\t"\
		"vmulpd	%%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmovaps	%%ymm0,%%ymm4	\n\t"\
		"vmovaps	%%ymm1,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm0,%%ymm0	\n\t"\
		"vsubpd	%%ymm2,%%ymm1,%%ymm1	\n\t"\
		"vsubpd	%%ymm3,%%ymm4,%%ymm4	\n\t"\
		"vaddpd	%%ymm2,%%ymm5,%%ymm5	\n\t"\
		"vmovaps	%%ymm0,0x380(%%r10)	\n\t	vmovaps	%%ymm1,0x3a0(%%r10)	\n\t"/* c11 = cc0 + 0x1c */\
		"vmovaps	%%ymm4,0x400(%%r10)	\n\t	vmovaps	%%ymm5,0x420(%%r10)	\n\t"/* c15 = cc0 + 0x20 */\
		:					/* outputs: none */\
		: [__cc0] "m" (Xcc0)	/* All inputs from memory addresses here */\
		 ,[__k0]  "m" (Xk0)\
		 ,[__k1]  "m" (Xk1)\
		 ,[__rt0] "m" (Xrt0)\
		 ,[__rt1] "m" (Xrt1)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
	}

#elif OS_BITS == 64	// 64-bit SSE2

  #ifdef USE_ARM_V8_SIMD

	// SIMD-gen 15 nontrivial twiddles; 1,2,4,8,13 via 2-table-mul, remaining 10 in 5 pair-cmuls from those:
	#define SSE2_RADIX16_CALC_TWIDDLES_LOACC(Xcc0,Xk0,Xk1,Xrt0,Xrt1)\
	{\
		__asm__ volatile (\
		"ldr	w0,%[__k0]		\n\t"\
		"ldr	w1,%[__k1]			\n\t"\
		"ldr	x2,%[__rt0]			\n\t"\
		"ldr	x4,%[__rt1]			\n\t"\
		"ldr	x6,%[__cc0]			\n\t"\
		/* c1,s1: */\
		"ldp	w7,w8,[x0]			\n\t"/* k0_arr[0,1] */\
		"add	x3, x2,x8			\n\t"/* rt0 + k0_arr[1], contig-data there = [re0.A,im0.A] */\
		"add	x2, x2,x7			\n\t"/* rt0 + k0_arr[0], contig-data there = [re0.B,im0.B] */\
		/* vector-data-to-be-interleaved are not nec. contiguous, so can't simply use 'ld2 {v0.2d,v1.2d},[base address]' to effect 2x2 transpose: */\
		"ld2	{v0.d,v1.d}[0],[x2]	\n\t"/* m0 = [re0.A,re0.B] */\
		"ld2	{v0.d,v1.d}[1],[x3]	\n\t"/* m1 = [im0.A,im0.B] */\
		"ldp	w7,w8,[x1]			\n\t"/* k1_arr[0,1] */\
		"add	x5, x4,x8			\n\t"/* rt1 + k1_arr[1], contig-data there = [re1.A,im1.A] */\
		"add	x4, x4,x7			\n\t"/* rt1 + k1_arr[0], contig-data there = [re1.B,im1.B] */\
		"ld2	{v2.d,v3.d}[0],[x4]	\n\t"/* m2 = [re1.A,re1.B] */\
		"ld2	{v2.d,v3.d}[1],[x5]	\n\t"/* m3 = [im1.A,im1.B] */\
		"fmul	v6.2d,v0.2d,v2.2d	\n\t"/* re0.re1 */\
		"fmul	v7.2d,v1.2d,v2.2d	\n\t"/* im0.re1 */\
		"fmls	v6.2d,v1.2d,v3.2d	\n\t"/* Re = re0.re1 - im0.im1 */\
		"fmla	v7.2d,v0.2d,v3.2d	\n\t"/* Im = im0.re1 + re0.im1 */\
		"stp	q6,q7,[x6,#0x120]	\n\t"/* cc0 + 0x12; keep persistent copies of c1,s1 in v6,7 */\
		/* c2,s2: */\
		"ldr	x2,%[__rt0]			\n\t"\
		"ldr	x4,%[__rt1]			\n\t"\
		"ldp	w7,w8,[x0,#0x08]	\n\t"/* k0_arr[2,3] */\
		"add	x3, x2,x8			\n\t"\
		"add	x2, x2,x7			\n\t"\
		"ld2	{v0.d,v1.d}[0],[x2]	\n\t"\
		"ld2	{v0.d,v1.d}[1],[x3]	\n\t"\
		"ldp	w7,w8,[x1,#0x08]	\n\t"/* k1_arr[2,3] */\
		"add	x5, x4,x8			\n\t"\
		"add	x4, x4,x7			\n\t"\
		"ld2	{v2.d,v3.d}[0],[x4]	\n\t"\
		"ld2	{v2.d,v3.d}[1],[x5]	\n\t"\
		"fmul	v8.2d,v0.2d,v2.2d	\n\t"\
		"fmul	v9.2d,v1.2d,v2.2d	\n\t"\
		"fmls	v8.2d,v1.2d,v3.2d	\n\t"\
		"fmla	v9.2d,v0.2d,v3.2d	\n\t"\
		"stp	q8,q9,[x6,#0x0a0]	\n\t"/* cc0 + 0xa; keep persistent copies of c2,s2 in v8,9 */\
		/* c4,s4: */\
		"ldr	x2,%[__rt0]			\n\t"\
		"ldr	x4,%[__rt1]			\n\t"\
		"ldp	w7,w8,[x0,#0x10]	\n\t"/* k0_arr[4,5] */\
		"add	x3, x2,x8			\n\t"\
		"add	x2, x2,x7			\n\t"\
		"ld2	{v0.d,v1.d}[0],[x2]	\n\t"\
		"ld2	{v0.d,v1.d}[1],[x3]	\n\t"\
		"ldp	w7,w8,[x1,#0x10]	\n\t"/* k1_arr[4,5] */\
		"add	x5, x4,x8			\n\t"\
		"add	x4, x4,x7			\n\t"\
		"ld2	{v2.d,v3.d}[0],[x4]	\n\t"\
		"ld2	{v2.d,v3.d}[1],[x5]	\n\t"\
		"fmul	v10.2d,v0.2d,v2.2d	\n\t"\
		"fmul	v11.2d,v1.2d,v2.2d	\n\t"\
		"fmls	v10.2d,v1.2d,v3.2d	\n\t"\
		"fmla	v11.2d,v0.2d,v3.2d	\n\t"\
		"stp	q10,q11,[x6,#0x060]	\n\t"/* cc0 + 0x6; keep persistent copies of c4,s4 in v10,11 */\
		/* c8,s8: */\
		"ldr	x2,%[__rt0]			\n\t"\
		"ldr	x4,%[__rt1]			\n\t"\
		"ldp	w7,w8,[x0,#0x18]	\n\t"/* k0_arr[6,7] */\
		"add	x3, x2,x8			\n\t"\
		"add	x2, x2,x7			\n\t"\
		"ld2	{v0.d,v1.d}[0],[x2]	\n\t"\
		"ld2	{v0.d,v1.d}[1],[x3]	\n\t"\
		"ldp	w7,w8,[x1,#0x18]	\n\t"/* k1_arr[6,7] */\
		"add	x5, x4,x8			\n\t"\
		"add	x4, x4,x7			\n\t"\
		"ld2	{v2.d,v3.d}[0],[x4]	\n\t"\
		"ld2	{v2.d,v3.d}[1],[x5]	\n\t"\
		"fmul	v12.2d,v0.2d,v2.2d	\n\t"\
		"fmul	v13.2d,v1.2d,v2.2d	\n\t"\
		"fmls	v12.2d,v1.2d,v3.2d	\n\t"\
		"fmla	v13.2d,v0.2d,v3.2d	\n\t"\
		"stp	q12,q13,[x6,#0x040]	\n\t"/* cc0 + 0x4; keep persistent copies of c8,s8 in v12,13 */\
		/* c13,s13: */\
		"ldr	x2,%[__rt0]			\n\t"\
		"ldr	x4,%[__rt1]			\n\t"\
		"ldp	w7,w8,[x0,#0x20]	\n\t"/* k0_arr[8,9] */\
		"add	x3, x2,x8			\n\t"\
		"add	x2, x2,x7			\n\t"\
		"ld2	{v0.d,v1.d}[0],[x2]	\n\t"\
		"ld2	{v0.d,v1.d}[1],[x3]	\n\t"\
		"ldp	w7,w8,[x1,#0x20]	\n\t"/* k1_arr[8,9] */\
		"add	x5, x4,x8			\n\t"\
		"add	x4, x4,x7			\n\t"\
		"ld2	{v2.d,v3.d}[0],[x4]	\n\t"\
		"ld2	{v2.d,v3.d}[1],[x5]	\n\t"\
		"fmul	v14.2d,v0.2d,v2.2d	\n\t"\
		"fmul	v15.2d,v1.2d,v2.2d	\n\t"\
		"fmls	v14.2d,v1.2d,v3.2d	\n\t"\
		"fmla	v15.2d,v0.2d,v3.2d	\n\t"\
		"stp	q14,q15,[x6,#0x180]	\n\t"/* cc0 + 0x18; keep persistent copies of c13,s13 in v14,15 */\
	/* Mem-mapping: (c,s)[0-15] :  0  1 2  3 4  5 6  7 8  9 a  b c  d  e  f  */\
	/*				(cc0,ss0) + 0x[2,12,a,1a,6,16,e,1e,4,14,c,1c,8,18,10,20] */\
	/* SSE2_CMUL_EXPO(c1,c4 ,c3 ,c5 ): */\
		"fmul	v0.2d,v6.2d,v10.2d	\n\t"/* c1.c4 */\
		"fmul	v1.2d,v6.2d,v11.2d	\n\t"/* c1.s4 */\
		"fmul	v2.2d,v7.2d,v10.2d	\n\t"/* s1.c4 */\
		"fmul	v3.2d,v7.2d,v11.2d	\n\t"/* s1.s4 */\
		"fsub	v4.2d,v0.2d,v3.2d	\n\t"/* c1.c4 - s1.s4 */\
		"fadd	v5.2d,v1.2d,v2.2d	\n\t"/* c1.s4 + s1.c4 */\
		"fadd	v0.2d,v0.2d,v3.2d	\n\t"/* c1.c4 + s1.s4 */\
		"fsub	v1.2d,v1.2d,v2.2d	\n\t"/* c1.s4 - s1.c4 */\
		"stp	q4,q5,[x6,#0x160]	\n\t"/* c5 = cc0 + 0x16 */\
		"stp	q0,q1,[x6,#0x1a0]	\n\t"/* c3 = cc0 + 0x1a */\
	/* SSE2_CMUL_EXPO(c1,c8 ,c7 ,c9 ): */\
		"fmul	v0.2d,v6.2d,v12.2d	\n\t"\
		"fmul	v1.2d,v6.2d,v13.2d	\n\t"\
		"fmul	v2.2d,v7.2d,v12.2d	\n\t"\
		"fmul	v3.2d,v7.2d,v13.2d	\n\t"\
		"fsub	v4.2d,v0.2d,v3.2d	\n\t"\
		"fadd	v5.2d,v1.2d,v2.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v3.2d	\n\t"\
		"fsub	v1.2d,v1.2d,v2.2d	\n\t"\
		"stp	q4,q5,[x6,#0x140]	\n\t"/* c9 = cc0 + 0x14 */\
		"stp	q0,q1,[x6,#0x1e0]	\n\t"/* c7 = cc0 + 0x1e */\
	/* SSE2_CMUL_EXPO(c2,c8 ,c6 ,c10): */\
		"fmul	v0.2d,v8.2d,v12.2d	\n\t"\
		"fmul	v1.2d,v8.2d,v13.2d	\n\t"\
		"fmul	v2.2d,v9.2d,v12.2d	\n\t"\
		"fmul	v3.2d,v9.2d,v13.2d	\n\t"\
		"fsub	v4.2d,v0.2d,v3.2d	\n\t"\
		"fadd	v5.2d,v1.2d,v2.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v3.2d	\n\t"\
		"fsub	v1.2d,v1.2d,v2.2d	\n\t"\
		"stp	q4,q5,[x6,#0x0c0]	\n\t"/* c10 = cc0 + 0xc */\
		"stp	q0,q1,[x6,#0x0e0]	\n\t"/* c6  = cc0 + 0xe */\
	/* SSE2_CMUL_EXPO(c1,c13,c12,c14): */\
		"fmul	v0.2d,v6.2d,v14.2d	\n\t"\
		"fmul	v1.2d,v6.2d,v15.2d	\n\t"\
		"fmul	v2.2d,v7.2d,v14.2d	\n\t"\
		"fmul	v3.2d,v7.2d,v15.2d	\n\t"\
		"fsub	v4.2d,v0.2d,v3.2d	\n\t"\
		"fadd	v5.2d,v1.2d,v2.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v3.2d	\n\t"\
		"fsub	v1.2d,v1.2d,v2.2d	\n\t"\
		"stp	q4,q5,[x6,#0x100]	\n\t"/* c14 = cc0 + 0x10 */\
		"stp	q0,q1,[x6,#0x080]	\n\t"/* c12 = cc0 + 0x08 */\
	/* SSE2_CMUL_EXPO(c2,c13,c11,c15): */\
		"fmul	v0.2d,v8.2d,v14.2d	\n\t"\
		"fmul	v1.2d,v8.2d,v15.2d	\n\t"\
		"fmul	v2.2d,v9.2d,v14.2d	\n\t"\
		"fmul	v3.2d,v9.2d,v15.2d	\n\t"\
		"fsub	v4.2d,v0.2d,v3.2d	\n\t"\
		"fadd	v5.2d,v1.2d,v2.2d	\n\t"\
		"fadd	v0.2d,v0.2d,v3.2d	\n\t"\
		"fsub	v1.2d,v1.2d,v2.2d	\n\t"\
		"stp	q4,q5,[x6,#0x200]	\n\t"/* c15 = cc0 + 0x20 */\
		"stp	q0,q1,[x6,#0x1c0]	\n\t"/* c11 = cc0 + 0x1c */\
		:					/* outputs: none */\
		: [__cc0] "m" (Xcc0)	/* All inputs from memory addresses here */\
		 ,[__k0]  "m" (Xk0)\
		 ,[__k1]  "m" (Xk1)\
		 ,[__rt0] "m" (Xrt0)\
		 ,[__rt1] "m" (Xrt1)\
		: "cc","memory","x0","x1","x2","x3","x4","x5","x6","x7","x8","v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11","v12","v13","v14","v15"	/* Clobbered registers */\
		);\
	}

  #else	// x86 SSE2:

	// SIMD-gen 15 nontrivial twiddles; 1,2,4,8,13 via 2-table-mul, remaining 10 in 5 pair-cmuls from those:
	#define SSE2_RADIX16_CALC_TWIDDLES_LOACC(Xcc0,Xk0,Xk1,Xrt0,Xrt1)\
	{\
		__asm__ volatile (\
		"movq	%[__k0]	,%%rax\n\t"\
		"movq	%[__k1]	,%%rbx\n\t"\
		"movq	%[__rt0],%%rcx\n\t"\
		"movq	%[__rt1],%%rdx\n\t"\
		"movq	%[__cc0],%%r10\n\t"\
		/* c1,s1: */\
		"movl	0x0(%%rax),%%esi		\n\t"/* k0_arr[0] */\
		"movl	0x4(%%rax),%%edi		\n\t"/* k0_arr[1] */\
		"movaps	(%%rcx,%%rsi),%%xmm0	\n\t"/* [re0,im0].A */\
		"movaps	(%%rcx,%%rdi),%%xmm2	\n\t"/* [re0,im0].B */\
		"movaps		%%xmm0 ,%%xmm1		\n\t"/* cpy xmm0 */\
		"shufpd	$0,%%xmm2,%%xmm0		\n\t"/* m0 = [re0.A,re0.B] */\
		"shufpd	$3,%%xmm2,%%xmm1		\n\t"/* m1 = [im0.A,im0.B] */\
		"movl	0x0(%%rbx),%%esi		\n\t"/* k1_arr[0] */\
		"movl	0x4(%%rbx),%%edi		\n\t"/* k1_arr[1] */\
		"movaps	(%%rdx,%%rsi),%%xmm2	\n\t"/* [re1,im1].A */\
		"movaps	(%%rdx,%%rdi),%%xmm4	\n\t"/* [re1,im1].B */\
		"movaps		%%xmm2 ,%%xmm3		\n\t"/* cpy xmm2 */\
		"shufpd	$0,%%xmm4,%%xmm2		\n\t"/* m2 = [re1.A,re1.B] */\
		"shufpd	$3,%%xmm4,%%xmm3		\n\t"/* m3 = [im1.A,im1.B] */\
		"movaps	%%xmm0,%%xmm4			\n\t"/* cpy re0 */\
		"movaps	%%xmm1,%%xmm5			\n\t"/* cpy im0 */\
		"mulpd	%%xmm2,%%xmm0			\n\t"/* re1*re0, overwrites re0 */\
		"mulpd	%%xmm3,%%xmm1			\n\t"/* im1*im0, overwrites im0 */\
		"mulpd	%%xmm5,%%xmm2			\n\t"/* im0*re1, overwrites re1 */\
		"mulpd	%%xmm4,%%xmm3			\n\t"/* re0*im1, overwrites im1 */\
		"subpd	%%xmm1,%%xmm0			\n\t"/* Re */\
		"addpd	%%xmm3,%%xmm2			\n\t"/* Im */\
		"movaps	%%xmm0,0x120(%%r10)		\n\t"/* cc0 + 0x12 */\
		"movaps	%%xmm2,0x130(%%r10)		\n\t"\
		"movaps	%%xmm0,%%xmm6 	\n\t	movaps	%%xmm2,%%xmm7 	\n\t"/* Register copies are free; stash c1,s1 in m6,7 */\
		/* c2,s2: */\
		"movl	0x8(%%rax),%%esi		\n\t"/* k0_arr[2] */\
		"movl	0xc(%%rax),%%edi		\n\t"/* k0_arr[3] */\
		"movaps	(%%rcx,%%rsi),%%xmm0	\n\t"\
		"movaps	(%%rcx,%%rdi),%%xmm2	\n\t"\
		"movaps		%%xmm0 ,%%xmm1		\n\t"\
		"shufpd	$0,%%xmm2,%%xmm0		\n\t"\
		"shufpd	$3,%%xmm2,%%xmm1		\n\t"\
		"movl	0x8(%%rbx),%%esi		\n\t"/* k1_arr[2] */\
		"movl	0xc(%%rbx),%%edi		\n\t"/* k1_arr[3] */\
		"movaps	(%%rdx,%%rsi),%%xmm2	\n\t"\
		"movaps	(%%rdx,%%rdi),%%xmm4	\n\t"\
		"movaps		%%xmm2 ,%%xmm3		\n\t"\
		"shufpd	$0,%%xmm4,%%xmm2		\n\t"\
		"shufpd	$3,%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,%%xmm4			\n\t"\
		"movaps	%%xmm1,%%xmm5			\n\t"\
		"mulpd	%%xmm2,%%xmm0			\n\t"\
		"mulpd	%%xmm3,%%xmm1			\n\t"\
		"mulpd	%%xmm5,%%xmm2			\n\t"\
		"mulpd	%%xmm4,%%xmm3			\n\t"\
		"subpd	%%xmm1,%%xmm0			\n\t"\
		"addpd	%%xmm3,%%xmm2			\n\t"\
		"movaps	%%xmm0,0x0a0(%%r10)		\n\t"/* cc0 + 0xa */\
		"movaps	%%xmm2,0x0b0(%%r10)		\n\t"\
		"movaps	%%xmm0,%%xmm8 	\n\t	movaps	%%xmm2,%%xmm9 	\n\t"/* stash c2,s2 in m8,9 */\
		/* c4,s4: */\
		"movl	0x10(%%rax),%%esi		\n\t"/* k0_arr[4] */\
		"movl	0x14(%%rax),%%edi		\n\t"/* k0_arr[5] */\
		"movaps	(%%rcx,%%rsi),%%xmm0	\n\t"\
		"movaps	(%%rcx,%%rdi),%%xmm2	\n\t"\
		"movaps		%%xmm0 ,%%xmm1		\n\t"\
		"shufpd	$0,%%xmm2,%%xmm0		\n\t"\
		"shufpd	$3,%%xmm2,%%xmm1		\n\t"\
		"movl	0x10(%%rbx),%%esi		\n\t"/* k1_arr[4] */\
		"movl	0x14(%%rbx),%%edi		\n\t"/* k1_arr[5] */\
		"movaps	(%%rdx,%%rsi),%%xmm2	\n\t"\
		"movaps	(%%rdx,%%rdi),%%xmm4	\n\t"\
		"movaps		%%xmm2 ,%%xmm3		\n\t"\
		"shufpd	$0,%%xmm4,%%xmm2		\n\t"\
		"shufpd	$3,%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,%%xmm4			\n\t"\
		"movaps	%%xmm1,%%xmm5			\n\t"\
		"mulpd	%%xmm2,%%xmm0			\n\t"\
		"mulpd	%%xmm3,%%xmm1			\n\t"\
		"mulpd	%%xmm5,%%xmm2			\n\t"\
		"mulpd	%%xmm4,%%xmm3			\n\t"\
		"subpd	%%xmm1,%%xmm0			\n\t"\
		"addpd	%%xmm3,%%xmm2			\n\t"\
		"movaps	%%xmm0,0x060(%%r10)		\n\t"/* cc0 + 0x6 */\
		"movaps	%%xmm2,0x070(%%r10)		\n\t"\
		"movaps	%%xmm0,%%xmm10	\n\t	movaps	%%xmm2,%%xmm11	\n\t"/* stash c4,s4 in m10,11 */\
		/* c8,s8: */\
		"movl	0x18(%%rax),%%esi		\n\t"/* k0_arr[6] */\
		"movl	0x1c(%%rax),%%edi		\n\t"/* k0_arr[7] */\
		"movaps	(%%rcx,%%rsi),%%xmm0	\n\t"\
		"movaps	(%%rcx,%%rdi),%%xmm2	\n\t"\
		"movaps		%%xmm0 ,%%xmm1		\n\t"\
		"shufpd	$0,%%xmm2,%%xmm0		\n\t"\
		"shufpd	$3,%%xmm2,%%xmm1		\n\t"\
		"movl	0x18(%%rbx),%%esi		\n\t"/* k1_arr[6] */\
		"movl	0x1c(%%rbx),%%edi		\n\t"/* k1_arr[7] */\
		"movaps	(%%rdx,%%rsi),%%xmm2	\n\t"\
		"movaps	(%%rdx,%%rdi),%%xmm4	\n\t"\
		"movaps		%%xmm2 ,%%xmm3		\n\t"\
		"shufpd	$0,%%xmm4,%%xmm2		\n\t"\
		"shufpd	$3,%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,%%xmm4			\n\t"\
		"movaps	%%xmm1,%%xmm5			\n\t"\
		"mulpd	%%xmm2,%%xmm0			\n\t"\
		"mulpd	%%xmm3,%%xmm1			\n\t"\
		"mulpd	%%xmm5,%%xmm2			\n\t"\
		"mulpd	%%xmm4,%%xmm3			\n\t"\
		"subpd	%%xmm1,%%xmm0			\n\t"\
		"addpd	%%xmm3,%%xmm2			\n\t"\
		"movaps	%%xmm0,0x040(%%r10)		\n\t"/* cc0 + 0x4 */\
		"movaps	%%xmm2,0x050(%%r10)		\n\t"\
		"movaps	%%xmm0,%%xmm12	\n\t	movaps	%%xmm2,%%xmm13	\n\t"/* stash c8,s8 in m12,13 */\
		/* c13,s13: */\
		"movl	0x20(%%rax),%%esi		\n\t"/* k0_arr[8] */\
		"movl	0x24(%%rax),%%edi		\n\t"/* k0_arr[9] */\
		"movaps	(%%rcx,%%rsi),%%xmm0	\n\t"\
		"movaps	(%%rcx,%%rdi),%%xmm2	\n\t"\
		"movaps		%%xmm0 ,%%xmm1		\n\t"\
		"shufpd	$0,%%xmm2,%%xmm0		\n\t"\
		"shufpd	$3,%%xmm2,%%xmm1		\n\t"\
		"movl	0x20(%%rbx),%%esi		\n\t"/* k1_arr[8] */\
		"movl	0x24(%%rbx),%%edi		\n\t"/* k1_arr[9] */\
		"movaps	(%%rdx,%%rsi),%%xmm2	\n\t"\
		"movaps	(%%rdx,%%rdi),%%xmm4	\n\t"\
		"movaps		%%xmm2 ,%%xmm3		\n\t"\
		"shufpd	$0,%%xmm4,%%xmm2		\n\t"\
		"shufpd	$3,%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,%%xmm4			\n\t"\
		"movaps	%%xmm1,%%xmm5			\n\t"\
		"mulpd	%%xmm2,%%xmm0			\n\t"\
		"mulpd	%%xmm3,%%xmm1			\n\t"\
		"mulpd	%%xmm5,%%xmm2			\n\t"\
		"mulpd	%%xmm4,%%xmm3			\n\t"\
		"subpd	%%xmm1,%%xmm0			\n\t"\
		"addpd	%%xmm3,%%xmm2			\n\t"\
		"movaps	%%xmm0,0x180(%%r10)		\n\t"/* cc0 + 0x18 */\
		"movaps	%%xmm2,0x190(%%r10)		\n\t"\
		"movaps	%%xmm0,%%xmm14	\n\t	movaps	%%xmm2,%%xmm15	\n\t"/* stash c13,s13 in m14,15 */\
	/* Mem-mapping: (c,s)[0-15] :  0  1 2  3 4  5 6  7 8  9 a  b c  d  e  f  */\
	/*				(cc0,ss0) + 0x[2,12,a,1a,6,16,e,1e,4,14,c,1c,8,18,10,20] */\
	/* SSE2_CMUL_EXPO(c1,c4 ,c3 ,c5 ): */\
		"movaps	%%xmm6 ,%%xmm0	\n\t	movaps	%%xmm7 ,%%xmm2\n\t"/* c1,s1 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c1,s1 */\
		"movaps	%%xmm10,%%xmm4	\n\t	movaps	%%xmm11,%%xmm5\n\t"/* c4,s4 */\
		"mulpd	%%xmm4,%%xmm0\n\t"\
		"mulpd	%%xmm5,%%xmm1\n\t"\
		"mulpd	%%xmm4,%%xmm2\n\t"\
		"mulpd	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm0,%%xmm4\n\t"\
		"movaps	%%xmm1,%%xmm5\n\t"\
		"addpd	%%xmm3,%%xmm0\n\t"\
		"subpd	%%xmm2,%%xmm1\n\t"\
		"subpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm2,%%xmm5\n\t"\
		"movaps	%%xmm0,0x1a0(%%r10)	\n\t	movaps	%%xmm1,0x1b0(%%r10)	\n\t"/* c3 = cc0 + 0x1a */\
		"movaps	%%xmm4,0x160(%%r10)	\n\t	movaps	%%xmm5,0x170(%%r10)	\n\t"/* c5 = cc0 + 0x16 */\
	/* SSE2_CMUL_EXPO(c1,c8 ,c7 ,c9 ): */\
		"movaps	%%xmm6 ,%%xmm0	\n\t	movaps	%%xmm7 ,%%xmm2\n\t"/* c1,s1 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c1,s1 */\
		"movaps	%%xmm12,%%xmm4	\n\t	movaps	%%xmm13,%%xmm5\n\t"/* c8,s8 */\
		"mulpd	%%xmm4,%%xmm0\n\t"\
		"mulpd	%%xmm5,%%xmm1\n\t"\
		"mulpd	%%xmm4,%%xmm2\n\t"\
		"mulpd	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm0,%%xmm4\n\t"\
		"movaps	%%xmm1,%%xmm5\n\t"\
		"addpd	%%xmm3,%%xmm0\n\t"\
		"subpd	%%xmm2,%%xmm1\n\t"\
		"subpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm2,%%xmm5\n\t"\
		"movaps	%%xmm0,0x1e0(%%r10)	\n\t	movaps	%%xmm1,0x1f0(%%r10)	\n\t"/* c7 = cc0 + 0x1e */\
		"movaps	%%xmm4,0x140(%%r10)	\n\t	movaps	%%xmm5,0x150(%%r10)	\n\t"/* c9 = cc0 + 0x14 */\
	/* SSE2_CMUL_EXPO(c2,c8 ,c6 ,c10): */\
		"movaps	%%xmm8 ,%%xmm0	\n\t	movaps	%%xmm9 ,%%xmm2\n\t"/* c2,s2 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c2,s2 */\
		"movaps	%%xmm12,%%xmm4	\n\t	movaps	%%xmm13,%%xmm5\n\t"/* c8,s8 */\
		"mulpd	%%xmm4,%%xmm0\n\t"\
		"mulpd	%%xmm5,%%xmm1\n\t"\
		"mulpd	%%xmm4,%%xmm2\n\t"\
		"mulpd	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm0,%%xmm4\n\t"\
		"movaps	%%xmm1,%%xmm5\n\t"\
		"addpd	%%xmm3,%%xmm0\n\t"\
		"subpd	%%xmm2,%%xmm1\n\t"\
		"subpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm2,%%xmm5\n\t"\
		"movaps	%%xmm0,0x0e0(%%r10)	\n\t	movaps	%%xmm1,0x0f0(%%r10)	\n\t"/* c6  = cc0 + 0xe */\
		"movaps	%%xmm4,0x0c0(%%r10)	\n\t	movaps	%%xmm5,0x0d0(%%r10)	\n\t"/* c10 = cc0 + 0xc */\
	/* SSE2_CMUL_EXPO(c1,c13,c12,c14): */\
		"movaps	%%xmm6 ,%%xmm0	\n\t	movaps	%%xmm7 ,%%xmm2\n\t"/* c1,s1 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c1,s1 */\
		"movaps	%%xmm14,%%xmm4	\n\t	movaps	%%xmm15,%%xmm5\n\t"/* c13,s13 */\
		"mulpd	%%xmm4,%%xmm0\n\t"\
		"mulpd	%%xmm5,%%xmm1\n\t"\
		"mulpd	%%xmm4,%%xmm2\n\t"\
		"mulpd	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm0,%%xmm4\n\t"\
		"movaps	%%xmm1,%%xmm5\n\t"\
		"addpd	%%xmm3,%%xmm0\n\t"\
		"subpd	%%xmm2,%%xmm1\n\t"\
		"subpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm2,%%xmm5\n\t"\
		"movaps	%%xmm0,0x080(%%r10)	\n\t	movaps	%%xmm1,0x090(%%r10)	\n\t"/* c12 = cc0 + 0x08 */\
		"movaps	%%xmm4,0x100(%%r10)	\n\t	movaps	%%xmm5,0x110(%%r10)	\n\t"/* c14 = cc0 + 0x10 */\
	/* SSE2_CMUL_EXPO(c2,c13,c11,c15): */\
		"movaps	%%xmm8 ,%%xmm0	\n\t	movaps	%%xmm9 ,%%xmm2\n\t"/* c2,s2 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c2,s2 */\
		"movaps	%%xmm14,%%xmm4	\n\t	movaps	%%xmm15,%%xmm5\n\t"/* c13,s13 */\
		"mulpd	%%xmm4,%%xmm0\n\t"\
		"mulpd	%%xmm5,%%xmm1\n\t"\
		"mulpd	%%xmm4,%%xmm2\n\t"\
		"mulpd	%%xmm5,%%xmm3\n\t"\
		"movaps	%%xmm0,%%xmm4\n\t"\
		"movaps	%%xmm1,%%xmm5\n\t"\
		"addpd	%%xmm3,%%xmm0\n\t"\
		"subpd	%%xmm2,%%xmm1\n\t"\
		"subpd	%%xmm3,%%xmm4\n\t"\
		"addpd	%%xmm2,%%xmm5\n\t"\
		"movaps	%%xmm0,0x1c0(%%r10)	\n\t	movaps	%%xmm1,0x1d0(%%r10)	\n\t"/* c11 = cc0 + 0x1c */\
		"movaps	%%xmm4,0x200(%%r10)	\n\t	movaps	%%xmm5,0x210(%%r10)	\n\t"/* c15 = cc0 + 0x20 */\
		:					/* outputs: none */\
		: [__cc0] "m" (Xcc0)	/* All inputs from memory addresses here */\
		 ,[__k0]  "m" (Xk0)\
		 ,[__k1]  "m" (Xk1)\
		 ,[__rt0] "m" (Xrt0)\
		 ,[__rt1] "m" (Xrt1)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
	}

  #endif	// ARMv8 or x86_64?

#elif OS_BITS == 32	// 32-bit SSE2

	// SIMD-gen 15 nontrivial twiddles; 1,2,4,8,13 via 2-table-mul, remaining 10 in 5 pair-cmuls from those:
	#define SSE2_RADIX16_CALC_TWIDDLES_1_2_4_8_13(Xcc0,Xk0,Xk1,Xrt0,Xrt1)\
	{\
		__asm__ volatile (\
		"movl	%[__rt0],%%ecx\n\t"\
		"movl	%[__rt1],%%edx\n\t"\
		/* c1,s1: */\
		"movl	%[__k0]	,%%eax\n\t"\
		"movl	0x0(%%eax),%%esi		\n\t"/* k0_arr[0] */\
		"movl	0x4(%%eax),%%edi		\n\t"/* k0_arr[1] */\
		"movaps	(%%ecx,%%esi),%%xmm0	\n\t"/* [re0,im0].A */\
		"movaps	(%%ecx,%%edi),%%xmm2	\n\t"/* [re0,im0].B */\
		"movaps		%%xmm0 ,%%xmm1		\n\t"/* cpy xmm0 */\
		"shufpd	$0,%%xmm2,%%xmm0		\n\t"/* m0 = [re0.A,re0.B] */\
		"shufpd	$3,%%xmm2,%%xmm1		\n\t"/* m1 = [im0.A,im0.B] */\
		"movl	%[__k1]	,%%eax\n\t"\
		"movl	0x0(%%eax),%%esi		\n\t"/* k1_arr[0] */\
		"movl	0x4(%%eax),%%edi		\n\t"/* k1_arr[1] */\
		"movaps	(%%edx,%%esi),%%xmm2	\n\t"/* [re1,im1].A */\
		"movaps	(%%edx,%%edi),%%xmm4	\n\t"/* [re1,im1].B */\
		"movaps		%%xmm2 ,%%xmm3		\n\t"/* cpy xmm2 */\
		"shufpd	$0,%%xmm4,%%xmm2		\n\t"/* m2 = [re1.A,re1.B] */\
		"shufpd	$3,%%xmm4,%%xmm3		\n\t"/* m3 = [im1.A,im1.B] */\
		"movaps	%%xmm0,%%xmm4			\n\t"/* cpy re0 */\
		"movaps	%%xmm1,%%xmm5			\n\t"/* cpy im0 */\
		"mulpd	%%xmm2,%%xmm0			\n\t"/* re1*re0, overwrites re0 */\
		"mulpd	%%xmm3,%%xmm1			\n\t"/* im1*im0, overwrites im0 */\
		"mulpd	%%xmm5,%%xmm2			\n\t"/* im0*re1, overwrites re1 */\
		"mulpd	%%xmm4,%%xmm3			\n\t"/* re0*im1, overwrites im1 */\
		"subpd	%%xmm1,%%xmm0			\n\t"/* Re */\
		"addpd	%%xmm3,%%xmm2			\n\t"/* Im */\
		"movl	%[__cc0],%%esi\n\t"\
		"movaps	%%xmm0,0x120(%%esi)		\n\t"/* cc0 + 0x12 */\
		"movaps	%%xmm2,0x130(%%esi)		\n\t"\
		/* c2,s2: */\
		"movl	%[__k0]	,%%eax\n\t"\
		"movl	0x8(%%eax),%%esi		\n\t"/* k0_arr[2] */\
		"movl	0xc(%%eax),%%edi		\n\t"/* k0_arr[3] */\
		"movaps	(%%ecx,%%esi),%%xmm0	\n\t"\
		"movaps	(%%ecx,%%edi),%%xmm2	\n\t"\
		"movaps		%%xmm0 ,%%xmm1		\n\t"\
		"shufpd	$0,%%xmm2,%%xmm0		\n\t"\
		"shufpd	$3,%%xmm2,%%xmm1		\n\t"\
		"movl	%[__k1]	,%%eax\n\t"\
		"movl	0x8(%%eax),%%esi		\n\t"/* k1_arr[2] */\
		"movl	0xc(%%eax),%%edi		\n\t"/* k1_arr[3] */\
		"movaps	(%%edx,%%esi),%%xmm2	\n\t"\
		"movaps	(%%edx,%%edi),%%xmm4	\n\t"\
		"movaps		%%xmm2 ,%%xmm3		\n\t"\
		"shufpd	$0,%%xmm4,%%xmm2		\n\t"\
		"shufpd	$3,%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,%%xmm4			\n\t"\
		"movaps	%%xmm1,%%xmm5			\n\t"\
		"mulpd	%%xmm2,%%xmm0			\n\t"\
		"mulpd	%%xmm3,%%xmm1			\n\t"\
		"mulpd	%%xmm5,%%xmm2			\n\t"\
		"mulpd	%%xmm4,%%xmm3			\n\t"\
		"subpd	%%xmm1,%%xmm0			\n\t"\
		"addpd	%%xmm3,%%xmm2			\n\t"\
		"movl	%[__cc0],%%esi\n\t"\
		"movaps	%%xmm0,0x0a0(%%esi)		\n\t"/* cc0 + 0xa */\
		"movaps	%%xmm2,0x0b0(%%esi)		\n\t"\
		/* c4,s4: */\
		"movl	%[__k0]	,%%eax\n\t"\
		"movl	0x10(%%eax),%%esi		\n\t"/* k0_arr[4] */\
		"movl	0x14(%%eax),%%edi		\n\t"/* k0_arr[5] */\
		"movaps	(%%ecx,%%esi),%%xmm0	\n\t"\
		"movaps	(%%ecx,%%edi),%%xmm2	\n\t"\
		"movaps		%%xmm0 ,%%xmm1		\n\t"\
		"shufpd	$0,%%xmm2,%%xmm0		\n\t"\
		"shufpd	$3,%%xmm2,%%xmm1		\n\t"\
		"movl	%[__k1]	,%%eax\n\t"\
		"movl	0x10(%%eax),%%esi		\n\t"/* k1_arr[4] */\
		"movl	0x14(%%eax),%%edi		\n\t"/* k1_arr[5] */\
		"movaps	(%%edx,%%esi),%%xmm2	\n\t"\
		"movaps	(%%edx,%%edi),%%xmm4	\n\t"\
		"movaps		%%xmm2 ,%%xmm3		\n\t"\
		"shufpd	$0,%%xmm4,%%xmm2		\n\t"\
		"shufpd	$3,%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,%%xmm4			\n\t"\
		"movaps	%%xmm1,%%xmm5			\n\t"\
		"mulpd	%%xmm2,%%xmm0			\n\t"\
		"mulpd	%%xmm3,%%xmm1			\n\t"\
		"mulpd	%%xmm5,%%xmm2			\n\t"\
		"mulpd	%%xmm4,%%xmm3			\n\t"\
		"subpd	%%xmm1,%%xmm0			\n\t"\
		"addpd	%%xmm3,%%xmm2			\n\t"\
		"movl	%[__cc0],%%esi\n\t"\
		"movaps	%%xmm0,0x060(%%esi)		\n\t"/* cc0 + 0x6 */\
		"movaps	%%xmm2,0x070(%%esi)		\n\t"\
		/* c8,s8: */\
		"movl	%[__k0]	,%%eax\n\t"\
		"movl	0x18(%%eax),%%esi		\n\t"/* k0_arr[6] */\
		"movl	0x1c(%%eax),%%edi		\n\t"/* k0_arr[7] */\
		"movaps	(%%ecx,%%esi),%%xmm0	\n\t"\
		"movaps	(%%ecx,%%edi),%%xmm2	\n\t"\
		"movaps		%%xmm0 ,%%xmm1		\n\t"\
		"shufpd	$0,%%xmm2,%%xmm0		\n\t"\
		"shufpd	$3,%%xmm2,%%xmm1		\n\t"\
		"movl	%[__k1]	,%%eax\n\t"\
		"movl	0x18(%%eax),%%esi		\n\t"/* k1_arr[6] */\
		"movl	0x1c(%%eax),%%edi		\n\t"/* k1_arr[7] */\
		"movaps	(%%edx,%%esi),%%xmm2	\n\t"\
		"movaps	(%%edx,%%edi),%%xmm4	\n\t"\
		"movaps		%%xmm2 ,%%xmm3		\n\t"\
		"shufpd	$0,%%xmm4,%%xmm2		\n\t"\
		"shufpd	$3,%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,%%xmm4			\n\t"\
		"movaps	%%xmm1,%%xmm5			\n\t"\
		"mulpd	%%xmm2,%%xmm0			\n\t"\
		"mulpd	%%xmm3,%%xmm1			\n\t"\
		"mulpd	%%xmm5,%%xmm2			\n\t"\
		"mulpd	%%xmm4,%%xmm3			\n\t"\
		"subpd	%%xmm1,%%xmm0			\n\t"\
		"addpd	%%xmm3,%%xmm2			\n\t"\
		"movl	%[__cc0],%%esi\n\t"\
		"movaps	%%xmm0,0x040(%%esi)		\n\t"/* cc0 + 0x4 */\
		"movaps	%%xmm2,0x050(%%esi)		\n\t"\
		/* c13,s13: */\
		"movl	%[__k0]	,%%eax\n\t"\
		"movl	0x20(%%eax),%%esi		\n\t"/* k0_arr[8] */\
		"movl	0x24(%%eax),%%edi		\n\t"/* k0_arr[9] */\
		"movaps	(%%ecx,%%esi),%%xmm0	\n\t"\
		"movaps	(%%ecx,%%edi),%%xmm2	\n\t"\
		"movaps		%%xmm0 ,%%xmm1		\n\t"\
		"shufpd	$0,%%xmm2,%%xmm0		\n\t"\
		"shufpd	$3,%%xmm2,%%xmm1		\n\t"\
		"movl	%[__k1]	,%%eax\n\t"\
		"movl	0x20(%%eax),%%esi		\n\t"/* k1_arr[8] */\
		"movl	0x24(%%eax),%%edi		\n\t"/* k1_arr[9] */\
		"movaps	(%%edx,%%esi),%%xmm2	\n\t"\
		"movaps	(%%edx,%%edi),%%xmm4	\n\t"\
		"movaps		%%xmm2 ,%%xmm3		\n\t"\
		"shufpd	$0,%%xmm4,%%xmm2		\n\t"\
		"shufpd	$3,%%xmm4,%%xmm3		\n\t"\
		"movaps	%%xmm0,%%xmm4			\n\t"\
		"movaps	%%xmm1,%%xmm5			\n\t"\
		"mulpd	%%xmm2,%%xmm0			\n\t"\
		"mulpd	%%xmm3,%%xmm1			\n\t"\
		"mulpd	%%xmm5,%%xmm2			\n\t"\
		"mulpd	%%xmm4,%%xmm3			\n\t"\
		"subpd	%%xmm1,%%xmm0			\n\t"\
		"addpd	%%xmm3,%%xmm2			\n\t"\
		"movl	%[__cc0],%%esi\n\t"\
		"movaps	%%xmm0,0x180(%%esi)		\n\t"/* cc0 + 0x18 */\
		"movaps	%%xmm2,0x190(%%esi)		\n\t"\
		:					/* outputs: none */\
		: [__cc0] "m" (Xcc0)	/* All inputs from memory addresses here */\
		 ,[__k0]  "m" (Xk0)\
		 ,[__k1]  "m" (Xk1)\
		 ,[__rt0] "m" (Xrt0)\
		 ,[__rt1] "m" (Xrt1)\
		: "cc","memory","eax","ecx","edx","esi","edi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5"		/* Clobbered registers */\
		);\
	}

#else

	#error Unhandled combination of #defs!

#endif	// SSE2 or AVX?

#endif	/* radix16_utils_asm_h_included */

