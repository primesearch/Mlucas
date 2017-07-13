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
#ifndef radix32_utils_asm_h_included
#define radix32_utils_asm_h_included

#ifdef USE_AVX512	// AVX512 implements a 512-bit-register version of the the AVX2 ALL_FMA-macro

	// SIMD-gen 31 nontrivial twiddles; 1,2,3,7,14,21,28, via 2-table-mul, remaining 24 in 12 pair-cmuls from those:
	#define SSE2_RADIX32_CALC_TWIDDLES_LOACC(Xcc0,Xk0,Xk1,Xrt0,Xrt1)\
	{\
		__asm__ volatile (\
		/* Note bytewidth-multiplier = 1 in the vector-addressing since already have incorporated into k1,2_arr data */\
		"movq	%[__k0]	,%%rax	\n\t"\
		"movq	%[__k1]	,%%rbx	\n\t"\
		"movq	%[__rt0],%%rcx	\n\t"\
		"movq	%[__rt1],%%rdx	\n\t"\
		"movq	%[__cc0],%%r10	\n\t"\
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
		"vmovaps	%%zmm0,0x980(%%r10)		\n\t"/* cc0 + 0x26 */\
		"vmovaps	%%zmm2,0x9c0(%%r10)		\n\t"\
		"vmovaps	%%zmm0,%%zmm6 	\n\t	vmovaps	%%zmm2,%%zmm7 	\n\t"/* stash c1,s1 in m6,7 */\
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
		"vmovaps	%%zmm0,0x580(%%r10)		\n\t"/* cc0 + 0x16 */\
		"vmovaps	%%zmm2,0x5c0(%%r10)		\n\t"\
		"vmovaps	%%zmm0,%%zmm8 	\n\t	vmovaps	%%zmm2,%%zmm9 	\n\t"/* stash c2,s2 in m8,9 */\
	/* c3,s3: */\
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
		"vmovaps	%%zmm0,0xd80(%%r10)		\n\t"/* cc0 + 0x36 */\
		"vmovaps	%%zmm2,0xdc0(%%r10)		\n\t"\
		"vmovaps	%%zmm0,%%zmm10	\n\t	vmovaps	%%zmm2,%%zmm11	\n\t"/* stash c3,s3 in m10,11 */\
	/* c7,s7: */\
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
		"vmovaps	%%zmm0,0xf80(%%r10)		\n\t"/* cc0 + 0x3e */\
		"vmovaps	%%zmm2,0xfc0(%%r10)		\n\t"\
		"vmovaps	%%zmm0,%%zmm12	\n\t	vmovaps	%%zmm2,%%zmm13	\n\t"/* stash c7,s7 in m12,13 */\
	/* c0E,s0E: */\
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
		"vmovaps	%%zmm0,0x880(%%r10)		\n\t"/* cc0 + 0x22 */\
		"vmovaps	%%zmm2,0x8c0(%%r10)		\n\t"\
		"vmovaps	%%zmm0,%%zmm14	\n\t	vmovaps	%%zmm2,%%zmm15	\n\t"/* stash c0E,s0E in m14,15 */\
	/* c15,s15 (hex 15, i.e. 0x15 = 21): */\
		"vmovaps	0xa0(%%rax),%%ymm4		\n\t"/* k0_arr[40-47] */\
		"vmovaps	0xa0(%%rbx),%%ymm5		\n\t"/* k1_arr[40-47] */\
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
		"vmovaps	%%zmm0,0xc00(%%r10)		\n\t"/* cc0 + 0x30 */\
		"vmovaps	%%zmm2,0xc40(%%r10)		\n\t"\
	/* c1C,s1C (0x1C = 28): */\
		"vmovaps	0xc0(%%rax),%%ymm4		\n\t"/* k0_arr[48-55] */\
		"vmovaps	0xc0(%%rbx),%%ymm5		\n\t"/* k1_arr[48-66] */\
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
		"vmovaps	%%zmm0,0x500(%%r10)		\n\t"/* cc0 + 0x14 */\
		"vmovaps	%%zmm2,0x540(%%r10)		\n\t"\
	/* Mem-mapping:   ** ** **          **                   **                    *                    * [** = precomputed, in-reg; * = precomputed, in-mem]
	cs1 in xmm6,7; cs2 in xmm8,9; cs3 in xmm10,11; cs7 in xmm12,13; csE in xmm14,15; cs15 in cc0+0x30; cs1C in cc0+0x14
	(c,s)[0-31]:cc[00 01 02 03 04 05 06 07 08 09 0A 0B 0C 0D 0E 0F 10 11 12 13 14 15 16 17 18 19 1A 1B 1C 1D 1E 1F]
	(cc0,ss0) + 0x[06,26,16,36,0e,2e,1e,3e,0a,2a,1a,3a,12,32,22,42,08,28,18,38,10,30,20,40,0c,2c,1c,3c,14,34,24,44].
	*/\
	/* SSE2_CMUL_EXPO(c01,c07,c06,c08): */\
		"vmovaps	%%zmm6 ,%%zmm0	\n\t	vmovaps	%%zmm7 ,%%zmm2\n\t"/* c1,s1 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c1,s1 */\
		"vmovaps	%%zmm12,%%zmm4	\n\t	vmovaps	%%zmm13,%%zmm5\n\t"/* c7,s7 */\
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
		"vmovaps	%%zmm0,0x780(%%r10)	\n\t	vmovaps	%%zmm1,0x7c0(%%r10)	\n\t"/* c6 = cc0 + 0x1e */\
		"vmovaps	%%zmm4,0x280(%%r10)	\n\t	vmovaps	%%zmm5,0x2c0(%%r10)	\n\t"/* c8 = cc0 + 0x0a */\
	/* SSE2_CMUL_EXPO(c02,c07,c05,c09): */\
		"vmovaps	%%zmm8 ,%%zmm0	\n\t	vmovaps	%%zmm9 ,%%zmm2\n\t"/* c2,s2 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c2,s2 */\
		"vmovaps	%%zmm12,%%zmm4	\n\t	vmovaps	%%zmm13,%%zmm5\n\t"/* c7,s7 */\
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
		"vmovaps	%%zmm0,0xb80(%%r10)	\n\t	vmovaps	%%zmm1,0xbc0(%%r10)	\n\t"/* c5 = cc0 + 0x2e */\
		"vmovaps	%%zmm4,0xa80(%%r10)	\n\t	vmovaps	%%zmm5,0xac0(%%r10)	\n\t"/* c9 = cc0 + 0x2a */\
	/* SSE2_CMUL_EXPO(c03,c07,c04,c0A): */\
		"vmovaps	%%zmm10,%%zmm0	\n\t	vmovaps	%%zmm11,%%zmm2\n\t"/* c3,s3 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c3,s3 */\
		"vmovaps	%%zmm12,%%zmm4	\n\t	vmovaps	%%zmm13,%%zmm5\n\t"/* c7,s7 */\
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
		"vmovaps	%%zmm0,0x380(%%r10)	\n\t	vmovaps	%%zmm1,0x3c0(%%r10)	\n\t"/* c4 = cc0 + 0x0e */\
		"vmovaps	%%zmm4,0x680(%%r10)	\n\t	vmovaps	%%zmm5,0x6c0(%%r10)	\n\t"/* cA = cc0 + 0x1a */\
	/* SSE2_CMUL_EXPO(c01,c0E,c0D,c0F): */\
		"vmovaps	%%zmm6 ,%%zmm0	\n\t	vmovaps	%%zmm7 ,%%zmm2\n\t"/* c1,s1 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c1,s1 */\
		"vmovaps	%%zmm14,%%zmm4	\n\t	vmovaps	%%zmm15,%%zmm5\n\t"/* cE,sE */\
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
		"vmovaps	%%zmm0,0xc80(%%r10)	\n\t	vmovaps	%%zmm1,0xcc0(%%r10)	\n\t"/* cD = cc0 + 0x32 */\
		"vmovaps	%%zmm4,0x1080(%%r10)\n\t	vmovaps	%%zmm5,0x10c0(%%r10)\n\t"/* cF = cc0 + 0x42 */\
	/* SSE2_CMUL_EXPO(c02,c0E,c0C,c10): */\
		"vmovaps	%%zmm8 ,%%zmm0	\n\t	vmovaps	%%zmm9 ,%%zmm2\n\t"/* c2,s2 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c2,s2 */\
		"vmovaps	%%zmm14,%%zmm4	\n\t	vmovaps	%%zmm15,%%zmm5\n\t"/* cE,sE */\
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
		"vmovaps	%%zmm0,0x480(%%r10)	\n\t	vmovaps	%%zmm1,0x4c0(%%r10)	\n\t"/* c5 = cc0 + 0x12 */\
		"vmovaps	%%zmm4,0x200(%%r10)	\n\t	vmovaps	%%zmm5,0x240(%%r10)	\n\t"/* c9 = cc0 + 0x08 */\
	/* SSE2_CMUL_EXPO(c03,c0E,c0B,c11): */\
		"vmovaps	%%zmm10,%%zmm0	\n\t	vmovaps	%%zmm11,%%zmm2\n\t"/* c3,s3 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c3,s3 */\
		"vmovaps	%%zmm14,%%zmm4	\n\t	vmovaps	%%zmm15,%%zmm5\n\t"/* cE,sE */\
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
		"vmovaps	%%zmm0,0xe80(%%r10)	\n\t	vmovaps	%%zmm1,0xec0(%%r10)	\n\t"/* c0B = cc0 + 0x3a */\
		"vmovaps	%%zmm4,0xa00(%%r10)	\n\t	vmovaps	%%zmm5,0xa40(%%r10)	\n\t"/* c11 = cc0 + 0x28 */\
	/* Done with cE,sE; move c15,s15 into zmm14,15: */\
	"vmovaps	0xc00(%%r10),%%zmm14	\n\t	vmovaps	0xc40(%%r10),%%zmm15\n\t"\
	/* SSE2_CMUL_EXPO(c01,c15,c14,c16): */\
		"vmovaps	%%zmm6 ,%%zmm0	\n\t	vmovaps	%%zmm7 ,%%zmm2\n\t"/* c1,s1 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c1,s1 */\
		"vmovaps	%%zmm14,%%zmm4	\n\t	vmovaps	%%zmm15,%%zmm5\n\t"/* c15,s15 */\
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
		"vmovaps	%%zmm0,0x400(%%r10)	\n\t	vmovaps	%%zmm1,0x440(%%r10)	\n\t"/* c14 = cc0 + 0x10 */\
		"vmovaps	%%zmm4,0x800(%%r10)	\n\t	vmovaps	%%zmm5,0x840(%%r10)	\n\t"/* c16 = cc0 + 0x20 */\
	/* SSE2_CMUL_EXPO(c02,c15,c13,c17): */\
		"vmovaps	%%zmm8 ,%%zmm0	\n\t	vmovaps	%%zmm9 ,%%zmm2\n\t"/* c2,s2 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c2,s2 */\
		"vmovaps	%%zmm14,%%zmm4	\n\t	vmovaps	%%zmm15,%%zmm5\n\t"/* c15,s15 */\
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
		"vmovaps	%%zmm0,0xe00(%%r10)	\n\t	vmovaps	%%zmm1,0xe40(%%r10)	\n\t"/* c13 = cc0 + 0x38 */\
		"vmovaps	%%zmm4,0x1000(%%r10)\n\t	vmovaps	%%zmm5,0x1040(%%r10)\n\t"/* c17 = cc0 + 0x40 */\
	/* SSE2_CMUL_EXPO(c03,c15,c12,c18): */\
		"vmovaps	%%zmm10,%%zmm0	\n\t	vmovaps	%%zmm11,%%zmm2\n\t"/* c3,s3 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c3,s3 */\
		"vmovaps	%%zmm14,%%zmm4	\n\t	vmovaps	%%zmm15,%%zmm5\n\t"/* c15,s15 */\
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
		"vmovaps	%%zmm0,0x600(%%r10)	\n\t	vmovaps	%%zmm1,0x640(%%r10)	\n\t"/* c12 = cc0 + 0x18 */\
		"vmovaps	%%zmm4,0x300(%%r10)	\n\t	vmovaps	%%zmm5,0x340(%%r10)	\n\t"/* c18 = cc0 + 0x0c */\
	/* Done with c15,s15; move c1C,s1C into zmm14,15: */\
	"vmovaps	0x500(%%r10),%%zmm14	\n\t	vmovaps	0x540(%%r10),%%zmm15\n\t"\
	/* SSE2_CMUL_EXPO(c01,c1C,c1B,c1D): */\
		"vmovaps	%%zmm6 ,%%zmm0	\n\t	vmovaps	%%zmm7 ,%%zmm2\n\t"/* c1,s1 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c1,s1 */\
		"vmovaps	%%zmm14,%%zmm4	\n\t	vmovaps	%%zmm15,%%zmm5\n\t"/* c1C,s1C */\
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
		"vmovaps	%%zmm0,0xf00(%%r10)	\n\t	vmovaps	%%zmm1,0xf40(%%r10)	\n\t"/* c1B = cc0 + 0x3c */\
		"vmovaps	%%zmm4,0xd00(%%r10)	\n\t	vmovaps	%%zmm5,0xd40(%%r10)	\n\t"/* c1D = cc0 + 0x34 */\
	/* SSE2_CMUL_EXPO(c02,c1C,c1A,c1E): */\
		"vmovaps	%%zmm8 ,%%zmm0	\n\t	vmovaps	%%zmm9 ,%%zmm2\n\t"/* c2,s2 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c2,s2 */\
		"vmovaps	%%zmm14,%%zmm4	\n\t	vmovaps	%%zmm15,%%zmm5\n\t"/* c1C,s1C */\
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
		"vmovaps	%%zmm0,0x700(%%r10)	\n\t	vmovaps	%%zmm1,0x740(%%r10)	\n\t"/* c1A = cc0 + 0x1c */\
		"vmovaps	%%zmm4,0x900(%%r10)	\n\t	vmovaps	%%zmm5,0x940(%%r10)	\n\t"/* c1E = cc0 + 0x24 */\
	/* SSE2_CMUL_EXPO(c03,c1C,c19,c1F): */\
		"vmovaps	%%zmm10,%%zmm0	\n\t	vmovaps	%%zmm11,%%zmm2\n\t"/* c3,s3 */\
		"vmovaps	%%zmm0 ,%%zmm1	\n\t	vmovaps	%%zmm2 ,%%zmm3\n\t"/* copy c3,s3 */\
		"vmovaps	%%zmm14,%%zmm4	\n\t	vmovaps	%%zmm15,%%zmm5\n\t"/* c1C,s1C */\
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
		"vmovaps	%%zmm0,0xb00(%%r10)	\n\t	vmovaps	%%zmm1,0xb40(%%r10)	\n\t"/* c19 = cc0 + 0x2c */\
		"vmovaps	%%zmm4,0x1100(%%r10)\n\t	vmovaps	%%zmm5,0x1140(%%r10)\n\t"/* c1F = cc0 + 0x44 */\
		:					/* outputs: none */\
		: [__cc0] "m" (Xcc0)	/* All inputs from memory addresses here */\
		 ,[__k0]  "m" (Xk0)\
		 ,[__k1]  "m" (Xk1)\
		 ,[__rt0] "m" (Xrt0)\
		 ,[__rt1] "m" (Xrt1)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
	}

#elif defined(USE_AVX)	// AVX and AVX2 both use 256-bit registers

	// SIMD-gen 31 nontrivial twiddles; 1,2,3,7,14,21,28, via 2-table-mul, remaining 24 in 12 pair-cmuls from those:
	#define SSE2_RADIX32_CALC_TWIDDLES_LOACC(Xcc0,Xk0,Xk1,Xrt0,Xrt1)\
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
		"vmovaps	%%ymm0,0x4c0(%%r10)		\n\t"/* cc0 + 0x26 */\
		"vmovaps	%%ymm2,0x4e0(%%r10)		\n\t"\
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
		"vmovaps	%%ymm0,0x2c0(%%r10)		\n\t"/* cc0 + 0x16 */\
		"vmovaps	%%ymm2,0x2e0(%%r10)		\n\t"\
		"vmovaps	%%ymm0,%%ymm8 	\n\t	vmovaps	%%ymm2,%%ymm9 	\n\t"/* stash c2,s2 in m8,9 */\
		/* c3,s3: */\
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
		"vmovaps	%%ymm0,0x6c0(%%r10)		\n\t"/* cc0 + 0x36 */\
		"vmovaps	%%ymm2,0x6e0(%%r10)		\n\t"\
		"vmovaps	%%ymm0,%%ymm10	\n\t	vmovaps	%%ymm2,%%ymm11	\n\t"/* stash c3,s3 in m10,11 */\
		/* c7,s7: */\
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
		"vmovaps	%%ymm0,0x7c0(%%r10)		\n\t"/* cc0 + 0x3e */\
		"vmovaps	%%ymm2,0x7e0(%%r10)		\n\t"\
		"vmovaps	%%ymm0,%%ymm12	\n\t	vmovaps	%%ymm2,%%ymm13	\n\t"/* stash c7,s7 in m12,13 */\
		/* c0E,s0E: */\
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
		"vmovaps	%%ymm0,0x440(%%r10)		\n\t"/* cc0 + 0x22 */\
		"vmovaps	%%ymm2,0x460(%%r10)		\n\t"\
		"vmovaps	%%ymm0,%%ymm14	\n\t	vmovaps	%%ymm2,%%ymm15	\n\t"/* stash c0E,s0E in m14,15 */\
		/* c15,s15 (hex 15, i.e. 0x15 = 21): */\
		"movl	0x50(%%rax),%%esi		\n\t"/* k0_arr[20] */\
		"movl	0x54(%%rax),%%edi		\n\t"/* k0_arr[21] */\
		"vmovaps	(%%rcx,%%rsi),%%xmm0	\n\t"\
		"vmovaps	(%%rcx,%%rdi),%%xmm2	\n\t"\
		"movl	0x58(%%rax),%%esi		\n\t"/* k0_arr[22] */\
		"movl	0x5c(%%rax),%%edi		\n\t"/* k0_arr[23] */\
	"vinsertf128 $1,(%%rcx,%%rsi),%%ymm0,%%ymm0	\n\t"\
	"vinsertf128 $1,(%%rcx,%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vmovaps		%%ymm0 ,%%ymm1		\n\t"\
		"vshufpd $0x0,%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vshufpd $0xf,%%ymm2,%%ymm1,%%ymm1	\n\t"\
		"movl	0x50(%%rbx),%%esi		\n\t"/* k1_arr[20] */\
		"movl	0x54(%%rbx),%%edi		\n\t"/* k1_arr[21] */\
		"vmovaps	(%%rdx,%%rsi),%%xmm2	\n\t"\
		"vmovaps	(%%rdx,%%rdi),%%xmm4	\n\t"\
		"movl	0x58(%%rbx),%%esi		\n\t"/* k1_arr[23] */\
		"movl	0x5c(%%rbx),%%edi		\n\t"/* k1_arr[23] */\
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
		"vmovaps	%%ymm0,0x600(%%r10)		\n\t"/* cc0 + 0x30 */\
		"vmovaps	%%ymm2,0x620(%%r10)		\n\t"\
		/* c1C,s1C (0x1C = 28): */\
		"movl	0x60(%%rax),%%esi		\n\t"/* k0_arr[24] */\
		"movl	0x64(%%rax),%%edi		\n\t"/* k0_arr[25] */\
		"vmovaps	(%%rcx,%%rsi),%%xmm0	\n\t"\
		"vmovaps	(%%rcx,%%rdi),%%xmm2	\n\t"\
		"movl	0x68(%%rax),%%esi		\n\t"/* k0_arr[26] */\
		"movl	0x6c(%%rax),%%edi		\n\t"/* k0_arr[27] */\
	"vinsertf128 $1,(%%rcx,%%rsi),%%ymm0,%%ymm0	\n\t"\
	"vinsertf128 $1,(%%rcx,%%rdi),%%ymm2,%%ymm2	\n\t"\
		"vmovaps		%%ymm0 ,%%ymm1		\n\t"\
		"vshufpd $0x0,%%ymm2,%%ymm0,%%ymm0	\n\t"\
		"vshufpd $0xf,%%ymm2,%%ymm1,%%ymm1	\n\t"\
		"movl	0x60(%%rbx),%%esi		\n\t"/* k1_arr[24] */\
		"movl	0x64(%%rbx),%%edi		\n\t"/* k1_arr[25] */\
		"vmovaps	(%%rdx,%%rsi),%%xmm2	\n\t"\
		"vmovaps	(%%rdx,%%rdi),%%xmm4	\n\t"\
		"movl	0x68(%%rbx),%%esi		\n\t"/* k1_arr[26] */\
		"movl	0x6c(%%rbx),%%edi		\n\t"/* k1_arr[27] */\
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
		"vmovaps	%%ymm0,0x280(%%r10)		\n\t"/* cc0 + 0x14 */\
		"vmovaps	%%ymm2,0x2a0(%%r10)		\n\t"\
	/* Mem-mapping:   ** ** **          **                   **                    *                    * [** = precomputed, in-reg; * = precomputed, in-mem]
	cs1 in xmm6,7; cs2 in xmm8,9; cs3 in xmm10,11; cs7 in xmm12,13; csE in xmm14,15; cs15 in cc0+0x30; cs1C in cc0+0x14
	(c,s)[0-31]:cc[00 01 02 03 04 05 06 07 08 09 0A 0B 0C 0D 0E 0F 10 11 12 13 14 15 16 17 18 19 1A 1B 1C 1D 1E 1F]
	(cc0,ss0) + 0x[06,26,16,36,0e,2e,1e,3e,0a,2a,1a,3a,12,32,22,42,08,28,18,38,10,30,20,40,0c,2c,1c,3c,14,34,24,44].
	*/\
	/* SSE2_CMUL_EXPO(c01,c07,c06,c08): */\
		"vmovaps	%%ymm6 ,%%ymm0	\n\t	vmovaps	%%ymm7 ,%%ymm2\n\t"/* c1,s1 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c1,s1 */\
		"vmovaps	%%ymm12,%%ymm4	\n\t	vmovaps	%%ymm13,%%ymm5\n\t"/* c7,s7 */\
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
		"vmovaps	%%ymm0,0x3c0(%%r10)	\n\t	vmovaps	%%ymm1,0x3e0(%%r10)	\n\t"/* c6 = cc0 + 0x1e */\
		"vmovaps	%%ymm4,0x140(%%r10)	\n\t	vmovaps	%%ymm5,0x160(%%r10)	\n\t"/* c8 = cc0 + 0x0a */\
	/* SSE2_CMUL_EXPO(c02,c07,c05,c09): */\
		"vmovaps	%%ymm8 ,%%ymm0	\n\t	vmovaps	%%ymm9 ,%%ymm2\n\t"/* c2,s2 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c2,s2 */\
		"vmovaps	%%ymm12,%%ymm4	\n\t	vmovaps	%%ymm13,%%ymm5\n\t"/* c7,s7 */\
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
		"vmovaps	%%ymm0,0x5c0(%%r10)	\n\t	vmovaps	%%ymm1,0x5e0(%%r10)	\n\t"/* c5 = cc0 + 0x2e */\
		"vmovaps	%%ymm4,0x540(%%r10)	\n\t	vmovaps	%%ymm5,0x560(%%r10)	\n\t"/* c9 = cc0 + 0x2a */\
	/* SSE2_CMUL_EXPO(c03,c07,c04,c0A): */\
		"vmovaps	%%ymm10,%%ymm0	\n\t	vmovaps	%%ymm11,%%ymm2\n\t"/* c3,s3 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c3,s3 */\
		"vmovaps	%%ymm12,%%ymm4	\n\t	vmovaps	%%ymm13,%%ymm5\n\t"/* c7,s7 */\
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
		"vmovaps	%%ymm0,0x1c0(%%r10)	\n\t	vmovaps	%%ymm1,0x1e0(%%r10)	\n\t"/* c4 = cc0 + 0x0e */\
		"vmovaps	%%ymm4,0x340(%%r10)	\n\t	vmovaps	%%ymm5,0x360(%%r10)	\n\t"/* cA = cc0 + 0x1a */\
	/* SSE2_CMUL_EXPO(c01,c0E,c0D,c0F): */\
		"vmovaps	%%ymm6 ,%%ymm0	\n\t	vmovaps	%%ymm7 ,%%ymm2\n\t"/* c1,s1 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c1,s1 */\
		"vmovaps	%%ymm14,%%ymm4	\n\t	vmovaps	%%ymm15,%%ymm5\n\t"/* cE,sE */\
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
		"vmovaps	%%ymm0,0x640(%%r10)	\n\t	vmovaps	%%ymm1,0x660(%%r10)	\n\t"/* cD = cc0 + 0x32 */\
		"vmovaps	%%ymm4,0x840(%%r10)	\n\t	vmovaps	%%ymm5,0x860(%%r10)	\n\t"/* cF = cc0 + 0x42 */\
	/* SSE2_CMUL_EXPO(c02,c0E,c0C,c10): */\
		"vmovaps	%%ymm8 ,%%ymm0	\n\t	vmovaps	%%ymm9 ,%%ymm2\n\t"/* c2,s2 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c2,s2 */\
		"vmovaps	%%ymm14,%%ymm4	\n\t	vmovaps	%%ymm15,%%ymm5\n\t"/* cE,sE */\
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
		"vmovaps	%%ymm0,0x240(%%r10)	\n\t	vmovaps	%%ymm1,0x260(%%r10)	\n\t"/* c5 = cc0 + 0x12 */\
		"vmovaps	%%ymm4,0x100(%%r10)	\n\t	vmovaps	%%ymm5,0x120(%%r10)	\n\t"/* c9 = cc0 + 0x08 */\
	/* SSE2_CMUL_EXPO(c03,c0E,c0B,c11): */\
		"vmovaps	%%ymm10,%%ymm0	\n\t	vmovaps	%%ymm11,%%ymm2\n\t"/* c3,s3 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c3,s3 */\
		"vmovaps	%%ymm14,%%ymm4	\n\t	vmovaps	%%ymm15,%%ymm5\n\t"/* cE,sE */\
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
		"vmovaps	%%ymm0,0x740(%%r10)	\n\t	vmovaps	%%ymm1,0x760(%%r10)	\n\t"/* c0B = cc0 + 0x3a */\
		"vmovaps	%%ymm4,0x500(%%r10)	\n\t	vmovaps	%%ymm5,0x520(%%r10)	\n\t"/* c11 = cc0 + 0x28 */\
	/* Done with cE,sE; move c15,s15 into ymm14,15: */\
	"vmovaps	0x600(%%r10),%%ymm14	\n\t	vmovaps	0x620(%%r10),%%ymm15\n\t"\
	/* SSE2_CMUL_EXPO(c01,c15,c14,c16): */\
		"vmovaps	%%ymm6 ,%%ymm0	\n\t	vmovaps	%%ymm7 ,%%ymm2\n\t"/* c1,s1 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c1,s1 */\
		"vmovaps	%%ymm14,%%ymm4	\n\t	vmovaps	%%ymm15,%%ymm5\n\t"/* c15,s15 */\
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
		"vmovaps	%%ymm0,0x200(%%r10)	\n\t	vmovaps	%%ymm1,0x220(%%r10)	\n\t"/* c14 = cc0 + 0x10 */\
		"vmovaps	%%ymm4,0x400(%%r10)	\n\t	vmovaps	%%ymm5,0x420(%%r10)	\n\t"/* c16 = cc0 + 0x20 */\
	/* SSE2_CMUL_EXPO(c02,c15,c13,c17): */\
		"vmovaps	%%ymm8 ,%%ymm0	\n\t	vmovaps	%%ymm9 ,%%ymm2\n\t"/* c2,s2 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c2,s2 */\
		"vmovaps	%%ymm14,%%ymm4	\n\t	vmovaps	%%ymm15,%%ymm5\n\t"/* c15,s15 */\
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
		"vmovaps	%%ymm0,0x700(%%r10)	\n\t	vmovaps	%%ymm1,0x720(%%r10)	\n\t"/* c13 = cc0 + 0x38 */\
		"vmovaps	%%ymm4,0x800(%%r10)	\n\t	vmovaps	%%ymm5,0x820(%%r10)	\n\t"/* c17 = cc0 + 0x40 */\
	/* SSE2_CMUL_EXPO(c03,c15,c12,c18): */\
		"vmovaps	%%ymm10,%%ymm0	\n\t	vmovaps	%%ymm11,%%ymm2\n\t"/* c3,s3 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c3,s3 */\
		"vmovaps	%%ymm14,%%ymm4	\n\t	vmovaps	%%ymm15,%%ymm5\n\t"/* c15,s15 */\
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
		"vmovaps	%%ymm0,0x300(%%r10)	\n\t	vmovaps	%%ymm1,0x320(%%r10)	\n\t"/* c12 = cc0 + 0x18 */\
		"vmovaps	%%ymm4,0x180(%%r10)	\n\t	vmovaps	%%ymm5,0x1a0(%%r10)	\n\t"/* c18 = cc0 + 0x0c */\
	/* Done with c15,s15; move c1C,s1C into ymm14,15: */\
	"vmovaps	0x280(%%r10),%%ymm14	\n\t	vmovaps	0x2a0(%%r10),%%ymm15\n\t"\
	/* SSE2_CMUL_EXPO(c01,c1C,c1B,c1D): */\
		"vmovaps	%%ymm6 ,%%ymm0	\n\t	vmovaps	%%ymm7 ,%%ymm2\n\t"/* c1,s1 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c1,s1 */\
		"vmovaps	%%ymm14,%%ymm4	\n\t	vmovaps	%%ymm15,%%ymm5\n\t"/* c1C,s1C */\
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
		"vmovaps	%%ymm0,0x780(%%r10)	\n\t	vmovaps	%%ymm1,0x7a0(%%r10)	\n\t"/* c1B = cc0 + 0x3c */\
		"vmovaps	%%ymm4,0x680(%%r10)	\n\t	vmovaps	%%ymm5,0x6a0(%%r10)	\n\t"/* c1D = cc0 + 0x34 */\
	/* SSE2_CMUL_EXPO(c02,c1C,c1A,c1E): */\
		"vmovaps	%%ymm8 ,%%ymm0	\n\t	vmovaps	%%ymm9 ,%%ymm2\n\t"/* c2,s2 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c2,s2 */\
		"vmovaps	%%ymm14,%%ymm4	\n\t	vmovaps	%%ymm15,%%ymm5\n\t"/* c1C,s1C */\
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
		"vmovaps	%%ymm0,0x380(%%r10)	\n\t	vmovaps	%%ymm1,0x3a0(%%r10)	\n\t"/* c1A = cc0 + 0x1c */\
		"vmovaps	%%ymm4,0x480(%%r10)	\n\t	vmovaps	%%ymm5,0x4a0(%%r10)	\n\t"/* c1E = cc0 + 0x24 */\
	/* SSE2_CMUL_EXPO(c03,c1C,c19,c1F): */\
		"vmovaps	%%ymm10,%%ymm0	\n\t	vmovaps	%%ymm11,%%ymm2\n\t"/* c3,s3 */\
		"vmovaps	%%ymm0,%%ymm1	\n\t	vmovaps	%%ymm2,%%ymm3\n\t"/* copy c3,s3 */\
		"vmovaps	%%ymm14,%%ymm4	\n\t	vmovaps	%%ymm15,%%ymm5\n\t"/* c1C,s1C */\
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
		"vmovaps	%%ymm0,0x580(%%r10)	\n\t	vmovaps	%%ymm1,0x5a0(%%r10)	\n\t"/* c19 = cc0 + 0x2c */\
		"vmovaps	%%ymm4,0x880(%%r10)	\n\t	vmovaps	%%ymm5,0x8a0(%%r10)	\n\t"/* c1F = cc0 + 0x44 */\
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

	// SIMD-gen 31 nontrivial twiddles; 1,2,3,7,14,21,28, via 2-table-mul, remaining 24 in 12 pair-cmuls from those:
	#define SSE2_RADIX32_CALC_TWIDDLES_LOACC(Xcc0,Xk0,Xk1,Xrt0,Xrt1)\
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
		"movaps	%%xmm0,0x260(%%r10)		\n\t"/* cc0 + 0x26 */\
		"movaps	%%xmm2,0x270(%%r10)		\n\t"\
		"movaps	%%xmm0,%%xmm6 	\n\t	movaps	%%xmm2,%%xmm7 	\n\t"/* stash c1,s1 in m6,7 */\
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
		"movaps	%%xmm0,0x160(%%r10)		\n\t"/* cc0 + 0x16 */\
		"movaps	%%xmm2,0x170(%%r10)		\n\t"\
		"movaps	%%xmm0,%%xmm8 	\n\t	movaps	%%xmm2,%%xmm9 	\n\t"/* stash c2,s2 in m8,9 */\
		/* c3,s3: */\
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
		"movaps	%%xmm0,0x360(%%r10)		\n\t"/* cc0 + 0x36 */\
		"movaps	%%xmm2,0x370(%%r10)		\n\t"\
		"movaps	%%xmm0,%%xmm10	\n\t	movaps	%%xmm2,%%xmm11	\n\t"/* stash c3,s3 in m10,11 */\
		/* c7,s7: */\
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
		"movaps	%%xmm0,0x3e0(%%r10)		\n\t"/* cc0 + 0x3e */\
		"movaps	%%xmm2,0x3f0(%%r10)		\n\t"\
		"movaps	%%xmm0,%%xmm12	\n\t	movaps	%%xmm2,%%xmm13	\n\t"/* stash c7,s7 in m12,13 */\
		/* c0E,s0E: */\
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
		"movaps	%%xmm0,0x220(%%r10)		\n\t"/* cc0 + 0x22 */\
		"movaps	%%xmm2,0x230(%%r10)		\n\t"\
		"movaps	%%xmm0,%%xmm14	\n\t	movaps	%%xmm2,%%xmm15	\n\t"/* stash c0E,s0E in m14,15 */\
		/* c15,s15 (hex 15, i.e. 0x15 = 21): */\
		"movl	0x28(%%rax),%%esi		\n\t"/* k0_arr[10] */\
		"movl	0x2c(%%rax),%%edi		\n\t"/* k0_arr[11] */\
		"movaps	(%%rcx,%%rsi),%%xmm0	\n\t"\
		"movaps	(%%rcx,%%rdi),%%xmm2	\n\t"\
		"movaps		%%xmm0 ,%%xmm1		\n\t"\
		"shufpd	$0,%%xmm2,%%xmm0		\n\t"\
		"shufpd	$3,%%xmm2,%%xmm1		\n\t"\
		"movl	0x28(%%rbx),%%esi		\n\t"/* k1_arr[10] */\
		"movl	0x2c(%%rbx),%%edi		\n\t"/* k1_arr[11] */\
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
		"movaps	%%xmm0,0x300(%%r10)		\n\t"/* cc0 + 0x30 */\
		"movaps	%%xmm2,0x310(%%r10)		\n\t"\
		/* c1C,s1C (0x1C = 28): */\
		"movl	0x30(%%rax),%%esi		\n\t"/* k0_arr[12] */\
		"movl	0x34(%%rax),%%edi		\n\t"/* k0_arr[13] */\
		"movaps	(%%rcx,%%rsi),%%xmm0	\n\t"\
		"movaps	(%%rcx,%%rdi),%%xmm2	\n\t"\
		"movaps		%%xmm0 ,%%xmm1		\n\t"\
		"shufpd	$0,%%xmm2,%%xmm0		\n\t"\
		"shufpd	$3,%%xmm2,%%xmm1		\n\t"\
		"movl	0x30(%%rbx),%%esi		\n\t"/* k1_arr[12] */\
		"movl	0x34(%%rbx),%%edi		\n\t"/* k1_arr[13] */\
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
		"movaps	%%xmm0,0x140(%%r10)		\n\t"/* cc0 + 0x14 */\
		"movaps	%%xmm2,0x150(%%r10)		\n\t"\
	/* Mem-mapping:   ** ** **          **                   **                    *                    * [** = precomputed, in-reg; * = precomputed, in-mem]
	cs1 in xmm6,7; cs2 in xmm8,9; cs3 in xmm10,11; cs7 in xmm12,13; csE in xmm14,15; cs15 in cc0+0x30; cs1C in cc0+0x14
	(c,s)[0-31]:cc[00 01 02 03 04 05 06 07 08 09 0A 0B 0C 0D 0E 0F 10 11 12 13 14 15 16 17 18 19 1A 1B 1C 1D 1E 1F]
	(cc0,ss0) + 0x[06,26,16,36,0e,2e,1e,3e,0a,2a,1a,3a,12,32,22,42,08,28,18,38,10,30,20,40,0c,2c,1c,3c,14,34,24,44].
	*/\
	/* SSE2_CMUL_EXPO(c01,c07,c06,c08): */\
		"movaps	%%xmm6 ,%%xmm0	\n\t	movaps	%%xmm7 ,%%xmm2\n\t"/* c1,s1 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c1,s1 */\
		"movaps	%%xmm12,%%xmm4	\n\t	movaps	%%xmm13,%%xmm5\n\t"/* c7,s7 */\
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
		"movaps	%%xmm0,0x1e0(%%r10)	\n\t	movaps	%%xmm1,0x1f0(%%r10)	\n\t"/* c6 = cc0 + 0x1e */\
		"movaps	%%xmm4,0x0a0(%%r10)	\n\t	movaps	%%xmm5,0x0b0(%%r10)	\n\t"/* c8 = cc0 + 0x0a */\
	/* SSE2_CMUL_EXPO(c02,c07,c05,c09): */\
		"movaps	%%xmm8 ,%%xmm0	\n\t	movaps	%%xmm9 ,%%xmm2\n\t"/* c2,s2 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c2,s2 */\
		"movaps	%%xmm12,%%xmm4	\n\t	movaps	%%xmm13,%%xmm5\n\t"/* c7,s7 */\
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
		"movaps	%%xmm0,0x2e0(%%r10)	\n\t	movaps	%%xmm1,0x2f0(%%r10)	\n\t"/* c5 = cc0 + 0x2e */\
		"movaps	%%xmm4,0x2a0(%%r10)	\n\t	movaps	%%xmm5,0x2b0(%%r10)	\n\t"/* c9 = cc0 + 0x2a */\
	/* SSE2_CMUL_EXPO(c03,c07,c04,c0A): */\
		"movaps	%%xmm10,%%xmm0	\n\t	movaps	%%xmm11,%%xmm2\n\t"/* c3,s3 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c3,s3 */\
		"movaps	%%xmm12,%%xmm4	\n\t	movaps	%%xmm13,%%xmm5\n\t"/* c7,s7 */\
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
		"movaps	%%xmm0,0x0e0(%%r10)	\n\t	movaps	%%xmm1,0x0f0(%%r10)	\n\t"/* c4 = cc0 + 0x0e */\
		"movaps	%%xmm4,0x1a0(%%r10)	\n\t	movaps	%%xmm5,0x1b0(%%r10)	\n\t"/* cA = cc0 + 0x1a */\
	/* SSE2_CMUL_EXPO(c01,c0E,c0D,c0F): */\
		"movaps	%%xmm6 ,%%xmm0	\n\t	movaps	%%xmm7 ,%%xmm2\n\t"/* c1,s1 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c1,s1 */\
		"movaps	%%xmm14,%%xmm4	\n\t	movaps	%%xmm15,%%xmm5\n\t"/* cE,sE */\
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
		"movaps	%%xmm0,0x320(%%r10)	\n\t	movaps	%%xmm1,0x330(%%r10)	\n\t"/* cD = cc0 + 0x32 */\
		"movaps	%%xmm4,0x420(%%r10)	\n\t	movaps	%%xmm5,0x430(%%r10)	\n\t"/* cF = cc0 + 0x42 */\
	/* SSE2_CMUL_EXPO(c02,c0E,c0C,c10): */\
		"movaps	%%xmm8 ,%%xmm0	\n\t	movaps	%%xmm9 ,%%xmm2\n\t"/* c2,s2 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c2,s2 */\
		"movaps	%%xmm14,%%xmm4	\n\t	movaps	%%xmm15,%%xmm5\n\t"/* cE,sE */\
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
		"movaps	%%xmm0,0x120(%%r10)	\n\t	movaps	%%xmm1,0x130(%%r10)	\n\t"/* c5 = cc0 + 0x12 */\
		"movaps	%%xmm4,0x080(%%r10)	\n\t	movaps	%%xmm5,0x090(%%r10)	\n\t"/* c9 = cc0 + 0x08 */\
	/* SSE2_CMUL_EXPO(c03,c0E,c0B,c11): */\
		"movaps	%%xmm10,%%xmm0	\n\t	movaps	%%xmm11,%%xmm2\n\t"/* c3,s3 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c3,s3 */\
		"movaps	%%xmm14,%%xmm4	\n\t	movaps	%%xmm15,%%xmm5\n\t"/* cE,sE */\
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
		"movaps	%%xmm0,0x3a0(%%r10)	\n\t	movaps	%%xmm1,0x3b0(%%r10)	\n\t"/* c0B = cc0 + 0x3a */\
		"movaps	%%xmm4,0x280(%%r10)	\n\t	movaps	%%xmm5,0x290(%%r10)	\n\t"/* c11 = cc0 + 0x28 */\
	/* Done with cE,sE; move c15,s15 into xmm14,15: */\
	"movaps	0x300(%%r10),%%xmm14	\n\t	movaps	0x310(%%r10),%%xmm15\n\t"\
	/* SSE2_CMUL_EXPO(c01,c15,c14,c16): */\
		"movaps	%%xmm6 ,%%xmm0	\n\t	movaps	%%xmm7 ,%%xmm2\n\t"/* c1,s1 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c1,s1 */\
		"movaps	%%xmm14,%%xmm4	\n\t	movaps	%%xmm15,%%xmm5\n\t"/* c15,s15 */\
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
		"movaps	%%xmm0,0x100(%%r10)	\n\t	movaps	%%xmm1,0x110(%%r10)	\n\t"/* c14 = cc0 + 0x10 */\
		"movaps	%%xmm4,0x200(%%r10)	\n\t	movaps	%%xmm5,0x210(%%r10)	\n\t"/* c16 = cc0 + 0x20 */\
	/* SSE2_CMUL_EXPO(c02,c15,c13,c17): */\
		"movaps	%%xmm8 ,%%xmm0	\n\t	movaps	%%xmm9 ,%%xmm2\n\t"/* c2,s2 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c2,s2 */\
		"movaps	%%xmm14,%%xmm4	\n\t	movaps	%%xmm15,%%xmm5\n\t"/* c15,s15 */\
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
		"movaps	%%xmm0,0x380(%%r10)	\n\t	movaps	%%xmm1,0x390(%%r10)	\n\t"/* c13 = cc0 + 0x38 */\
		"movaps	%%xmm4,0x400(%%r10)	\n\t	movaps	%%xmm5,0x410(%%r10)	\n\t"/* c17 = cc0 + 0x40 */\
	/* SSE2_CMUL_EXPO(c03,c15,c12,c18): */\
		"movaps	%%xmm10,%%xmm0	\n\t	movaps	%%xmm11,%%xmm2\n\t"/* c3,s3 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c3,s3 */\
		"movaps	%%xmm14,%%xmm4	\n\t	movaps	%%xmm15,%%xmm5\n\t"/* c15,s15 */\
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
		"movaps	%%xmm0,0x180(%%r10)	\n\t	movaps	%%xmm1,0x190(%%r10)	\n\t"/* c12 = cc0 + 0x18 */\
		"movaps	%%xmm4,0x0c0(%%r10)	\n\t	movaps	%%xmm5,0x0d0(%%r10)	\n\t"/* c18 = cc0 + 0x0c */\
	/* Done with c15,s15; move c1C,s1C into xmm14,15: */\
	"movaps	0x140(%%r10),%%xmm14	\n\t	movaps	0x150(%%r10),%%xmm15\n\t"\
	/* SSE2_CMUL_EXPO(c01,c1C,c1B,c1D): */\
		"movaps	%%xmm6 ,%%xmm0	\n\t	movaps	%%xmm7 ,%%xmm2\n\t"/* c1,s1 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c1,s1 */\
		"movaps	%%xmm14,%%xmm4	\n\t	movaps	%%xmm15,%%xmm5\n\t"/* c1C,s1C */\
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
		"movaps	%%xmm0,0x3c0(%%r10)	\n\t	movaps	%%xmm1,0x3d0(%%r10)	\n\t"/* c1B = cc0 + 0x3c */\
		"movaps	%%xmm4,0x340(%%r10)	\n\t	movaps	%%xmm5,0x350(%%r10)	\n\t"/* c1D = cc0 + 0x34 */\
	/* SSE2_CMUL_EXPO(c02,c1C,c1A,c1E): */\
		"movaps	%%xmm8 ,%%xmm0	\n\t	movaps	%%xmm9 ,%%xmm2\n\t"/* c2,s2 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c2,s2 */\
		"movaps	%%xmm14,%%xmm4	\n\t	movaps	%%xmm15,%%xmm5\n\t"/* c1C,s1C */\
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
		"movaps	%%xmm0,0x1c0(%%r10)	\n\t	movaps	%%xmm1,0x1d0(%%r10)	\n\t"/* c1A = cc0 + 0x1c */\
		"movaps	%%xmm4,0x240(%%r10)	\n\t	movaps	%%xmm5,0x250(%%r10)	\n\t"/* c1E = cc0 + 0x24 */\
	/* SSE2_CMUL_EXPO(c03,c1C,c19,c1F): */\
		"movaps	%%xmm10,%%xmm0	\n\t	movaps	%%xmm11,%%xmm2\n\t"/* c3,s3 */\
		"movaps	%%xmm0,%%xmm1	\n\t	movaps	%%xmm2,%%xmm3\n\t"/* copy c3,s3 */\
		"movaps	%%xmm14,%%xmm4	\n\t	movaps	%%xmm15,%%xmm5\n\t"/* c1C,s1C */\
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
		"movaps	%%xmm0,0x2c0(%%r10)	\n\t	movaps	%%xmm1,0x2d0(%%r10)	\n\t"/* c19 = cc0 + 0x2c */\
		"movaps	%%xmm4,0x440(%%r10)	\n\t	movaps	%%xmm5,0x450(%%r10)	\n\t"/* c1F = cc0 + 0x44 */\
		:					/* outputs: none */\
		: [__cc0] "m" (Xcc0)	/* All inputs from memory addresses here */\
		 ,[__k0]  "m" (Xk0)\
		 ,[__k1]  "m" (Xk1)\
		 ,[__rt0] "m" (Xrt0)\
		 ,[__rt1] "m" (Xrt1)\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r10","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9","xmm10","xmm11","xmm12","xmm13","xmm14","xmm15"	/* Clobbered registers */\
		);\
	}

#elif OS_BITS == 32	// 32-bit SSE2

	// SIMD-gen 31 nontrivial twiddles; 1,2,3,7,14,21,28, via 2-table-mul, remaining 24 in 12 pair-cmuls from those:
	#define SSE2_RADIX32_CALC_TWIDDLES_1_2_3_7_14_21_28(Xcc0,Xk0,Xk1,Xrt0,Xrt1)\
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
		"movaps	%%xmm0,0x260(%%esi)		\n\t"/* cc0 + 0x26 */\
		"movaps	%%xmm2,0x270(%%esi)		\n\t"\
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
		"movaps	%%xmm0,0x160(%%esi)		\n\t"/* cc0 + 0x16 */\
		"movaps	%%xmm2,0x170(%%esi)		\n\t"\
		/* c3,s3: */\
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
		"movaps	%%xmm0,0x360(%%esi)		\n\t"/* cc0 + 0x36 */\
		"movaps	%%xmm2,0x370(%%esi)		\n\t"\
		/* c7,s7: */\
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
		"movaps	%%xmm0,0x3e0(%%esi)		\n\t"/* cc0 + 0x3e */\
		"movaps	%%xmm2,0x3f0(%%esi)		\n\t"\
		/* c0E,s0E: */\
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
		"movaps	%%xmm0,0x220(%%esi)		\n\t"/* cc0 + 0x22 */\
		"movaps	%%xmm2,0x230(%%esi)		\n\t"\
		/* c15,s15 (hex 15, i.e. 0x15 = 21): */\
		"movl	%[__k0]	,%%eax\n\t"\
		"movl	0x28(%%eax),%%esi		\n\t"/* k0_arr[10] */\
		"movl	0x2c(%%eax),%%edi		\n\t"/* k0_arr[11] */\
		"movaps	(%%ecx,%%esi),%%xmm0	\n\t"\
		"movaps	(%%ecx,%%edi),%%xmm2	\n\t"\
		"movaps		%%xmm0 ,%%xmm1		\n\t"\
		"shufpd	$0,%%xmm2,%%xmm0		\n\t"\
		"shufpd	$3,%%xmm2,%%xmm1		\n\t"\
		"movl	%[__k1]	,%%eax\n\t"\
		"movl	0x28(%%eax),%%esi		\n\t"/* k1_arr[10] */\
		"movl	0x2c(%%eax),%%edi		\n\t"/* k1_arr[11] */\
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
		"movaps	%%xmm0,0x300(%%esi)		\n\t"/* cc0 + 0x30 */\
		"movaps	%%xmm2,0x310(%%esi)		\n\t"\
		/* c1C,s1C (0x1C = 28): */\
		"movl	%[__k0]	,%%eax\n\t"\
		"movl	0x30(%%eax),%%esi		\n\t"/* k0_arr[12] */\
		"movl	0x34(%%eax),%%edi		\n\t"/* k0_arr[13] */\
		"movaps	(%%ecx,%%esi),%%xmm0	\n\t"\
		"movaps	(%%ecx,%%edi),%%xmm2	\n\t"\
		"movaps		%%xmm0 ,%%xmm1		\n\t"\
		"shufpd	$0,%%xmm2,%%xmm0		\n\t"\
		"shufpd	$3,%%xmm2,%%xmm1		\n\t"\
		"movl	%[__k1]	,%%eax\n\t"\
		"movl	0x30(%%eax),%%esi		\n\t"/* k1_arr[12] */\
		"movl	0x34(%%eax),%%edi		\n\t"/* k1_arr[13] */\
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
		"movaps	%%xmm0,0x140(%%esi)		\n\t"/* cc0 + 0x14 */\
		"movaps	%%xmm2,0x150(%%esi)		\n\t"\
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

#endif	/* radix32_utils_asm_h_included */

