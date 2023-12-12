/*******************************************************************************
*                                                                              *
*   (C) 1997-2021 by Ernst W. Mayer.                                           *
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

/************************************************************************/
/* SIMD macros for dyadic square, mul, submul a*(b-c) and mulsub a*b-c: */
/************************************************************************/

/*********** Mulsub a*b - c: ***********/

#ifdef USE_ARM_V8_SIMD

	// No Fermat-mod support on ARMv8, just supply a stub macro:
	#define SIMD_MULSUB(Xa,Xb,Xc)\
	{\
	  __asm__ volatile (\
		"ldr x0,%[__add0]	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xa)	/* All inputs from memory addresses here */\
		 ,[__bdd0] "m" (Xb)	/* B-array base-address */\
		 ,[__cdd0] "m" (Xc)	/* C-array base-address */\
		: "cc","memory","x0"	/* Clobbered registers */\
	  );\
	}

#elif defined(USE_AVX512)	// AVX512 implements a 512-bit-register version of the the AVX2 ALL_FMA-macro

	/*** Started with SUBMUL macro, simply moved vsubpd line-pair to just before output-writes and changed SIMD register indices
	in the vsubpd to match those of output-registers. Also reordered reads to 0213 to match order in which inputs are used: ***/
	#define SIMD_MULSUB(Xa,Xb,Xc)\
	{\
	  __asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movq	%[__bdd0],%%rbx		\n\t"\
		"movq	%[__cdd0],%%rcx		\n\t"\
		/* x0.y0: */							/* x10.y10: */\
		"vmovaps		 (%%rax),%%zmm0	\n\t	vmovaps	0x800(%%rax),%%zmm5	\n\t"/* x.re */\
		"vmovaps		 (%%rbx),%%zmm2	\n\t	vmovaps	0x800(%%rbx),%%zmm7	\n\t"/* y.re */\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t	vmovaps	0x840(%%rax),%%zmm6	\n\t"/* x.im */\
		"vmovaps	0x040(%%rbx),%%zmm3	\n\t	vmovaps	0x840(%%rbx),%%zmm8	\n\t"/* y.im */\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"/* copy x.re */\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"/* x.re *= y.re */\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"/* y.re *= x.im */\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"/* x.im *= y.im */\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"/* y.im *= x.re[copy] */\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"/* xy.re = x.re*y.re - x.im*y.im */\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"/* xy.im = y.re*x.im + y.im*x.re */\
	"vsubpd		 (%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x800(%%rcx),%%zmm5,%%zmm5	\n\t"/* re = xy.re - z.re */\
	"vsubpd	0x040(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0x840(%%rcx),%%zmm7,%%zmm7	\n\t"/* im = xy.im - z.im */\
		"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm5,0x800(%%rax)	\n\t"/* write z.re */\
		"vmovaps	%%zmm2,0x040(%%rax)	\n\t	vmovaps	%%zmm7,0x840(%%rax)	\n\t"/* write z.im */\
		/* x1.y1: */							/* x11.y11: */\
		"vmovaps	0x080(%%rax),%%zmm0	\n\t	vmovaps	0x880(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x080(%%rbx),%%zmm2	\n\t	vmovaps	0x880(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x0c0(%%rax),%%zmm1	\n\t	vmovaps	0x8c0(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x0c0(%%rbx),%%zmm3	\n\t	vmovaps	0x8c0(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x080(%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x880(%%rcx),%%zmm5,%%zmm5	\n\t"\
	"vsubpd	0x0c0(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0x8c0(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x080(%%rax)	\n\t	vmovaps	%%zmm5,0x880(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x0c0(%%rax)	\n\t	vmovaps	%%zmm7,0x8c0(%%rax)	\n\t"\
		/* x2.y2: */							/* x12.y12: */\
		"vmovaps	0x100(%%rax),%%zmm0	\n\t	vmovaps	0x900(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x100(%%rbx),%%zmm2	\n\t	vmovaps	0x900(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x140(%%rax),%%zmm1	\n\t	vmovaps	0x940(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x140(%%rbx),%%zmm3	\n\t	vmovaps	0x940(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x100(%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x900(%%rcx),%%zmm5,%%zmm5	\n\t"\
	"vsubpd	0x140(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0x940(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x100(%%rax)	\n\t	vmovaps	%%zmm5,0x900(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x140(%%rax)	\n\t	vmovaps	%%zmm7,0x940(%%rax)	\n\t"\
		/* x3.y3: */							/* x13.y13: */\
		"vmovaps	0x180(%%rax),%%zmm0	\n\t	vmovaps	0x980(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x180(%%rbx),%%zmm2	\n\t	vmovaps	0x980(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x1c0(%%rax),%%zmm1	\n\t	vmovaps	0x9c0(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x1c0(%%rbx),%%zmm3	\n\t	vmovaps	0x9c0(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x180(%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0x980(%%rcx),%%zmm5,%%zmm5	\n\t"\
	"vsubpd	0x1c0(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0x9c0(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x180(%%rax)	\n\t	vmovaps	%%zmm5,0x980(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x1c0(%%rax)	\n\t	vmovaps	%%zmm7,0x9c0(%%rax)	\n\t"\
		/* x4.y4: */							/* x14.y14: */\
		"vmovaps	0x200(%%rax),%%zmm0	\n\t	vmovaps	0xa00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x200(%%rbx),%%zmm2	\n\t	vmovaps	0xa00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x240(%%rax),%%zmm1	\n\t	vmovaps	0xa40(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x240(%%rbx),%%zmm3	\n\t	vmovaps	0xa40(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x200(%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0xa00(%%rcx),%%zmm5,%%zmm5	\n\t"\
	"vsubpd	0x240(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xa40(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x200(%%rax)	\n\t	vmovaps	%%zmm5,0xa00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x240(%%rax)	\n\t	vmovaps	%%zmm7,0xa40(%%rax)	\n\t"\
		/* x5.y5: */							/* x15.y15: */\
		"vmovaps	0x280(%%rax),%%zmm0	\n\t	vmovaps	0xa80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x280(%%rbx),%%zmm2	\n\t	vmovaps	0xa80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x2c0(%%rax),%%zmm1	\n\t	vmovaps	0xac0(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x2c0(%%rbx),%%zmm3	\n\t	vmovaps	0xac0(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x280(%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0xa80(%%rcx),%%zmm5,%%zmm5	\n\t"\
	"vsubpd	0x2c0(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xac0(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x280(%%rax)	\n\t	vmovaps	%%zmm5,0xa80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x2c0(%%rax)	\n\t	vmovaps	%%zmm7,0xac0(%%rax)	\n\t"\
		/* x6.y6: */							/* x16.y16: */\
		"vmovaps	0x300(%%rax),%%zmm0	\n\t	vmovaps	0xb00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x300(%%rbx),%%zmm2	\n\t	vmovaps	0xb00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x340(%%rax),%%zmm1	\n\t	vmovaps	0xb40(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x340(%%rbx),%%zmm3	\n\t	vmovaps	0xb40(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x300(%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0xb00(%%rcx),%%zmm5,%%zmm5	\n\t"\
	"vsubpd	0x340(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xb40(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x300(%%rax)	\n\t	vmovaps	%%zmm5,0xb00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x340(%%rax)	\n\t	vmovaps	%%zmm7,0xb40(%%rax)	\n\t"\
		/* x7.y7: */							/* x17.y17: */\
		"vmovaps	0x380(%%rax),%%zmm0	\n\t	vmovaps	0xb80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x380(%%rbx),%%zmm2	\n\t	vmovaps	0xb80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x3c0(%%rax),%%zmm1	\n\t	vmovaps	0xbc0(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x3c0(%%rbx),%%zmm3	\n\t	vmovaps	0xbc0(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x380(%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0xb80(%%rcx),%%zmm5,%%zmm5	\n\t"\
	"vsubpd	0x3c0(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xbc0(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x380(%%rax)	\n\t	vmovaps	%%zmm5,0xb80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x3c0(%%rax)	\n\t	vmovaps	%%zmm7,0xbc0(%%rax)	\n\t"\
		/* x8.y8: */							/* x18.y18: */\
		"vmovaps	0x400(%%rax),%%zmm0	\n\t	vmovaps	0xc00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x400(%%rbx),%%zmm2	\n\t	vmovaps	0xc00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x440(%%rax),%%zmm1	\n\t	vmovaps	0xc40(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x440(%%rbx),%%zmm3	\n\t	vmovaps	0xc40(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x400(%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0xc00(%%rcx),%%zmm5,%%zmm5	\n\t"\
	"vsubpd	0x440(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xc40(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x400(%%rax)	\n\t	vmovaps	%%zmm5,0xc00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x440(%%rax)	\n\t	vmovaps	%%zmm7,0xc40(%%rax)	\n\t"\
		/* x9.y9: */							/* x19.y19: */\
		"vmovaps	0x480(%%rax),%%zmm0	\n\t	vmovaps	0xc80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x480(%%rbx),%%zmm2	\n\t	vmovaps	0xc80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x4c0(%%rax),%%zmm1	\n\t	vmovaps	0xcc0(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x4c0(%%rbx),%%zmm3	\n\t	vmovaps	0xcc0(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x480(%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0xc80(%%rcx),%%zmm5,%%zmm5	\n\t"\
	"vsubpd	0x4c0(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xcc0(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x480(%%rax)	\n\t	vmovaps	%%zmm5,0xc80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x4c0(%%rax)	\n\t	vmovaps	%%zmm7,0xcc0(%%rax)	\n\t"\
		/* xA.yA: */							/* x1A.y1A: */\
		"vmovaps	0x500(%%rax),%%zmm0	\n\t	vmovaps	0xd00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x500(%%rbx),%%zmm2	\n\t	vmovaps	0xd00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x540(%%rax),%%zmm1	\n\t	vmovaps	0xd40(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x540(%%rbx),%%zmm3	\n\t	vmovaps	0xd40(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x500(%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0xd00(%%rcx),%%zmm5,%%zmm5	\n\t"\
	"vsubpd	0x540(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xd40(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x500(%%rax)	\n\t	vmovaps	%%zmm5,0xd00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x540(%%rax)	\n\t	vmovaps	%%zmm7,0xd40(%%rax)	\n\t"\
		/* xB.yB: */							/* x1B.y1B: */\
		"vmovaps	0x580(%%rax),%%zmm0	\n\t	vmovaps	0xd80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x580(%%rbx),%%zmm2	\n\t	vmovaps	0xd80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x5c0(%%rax),%%zmm1	\n\t	vmovaps	0xdc0(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x5c0(%%rbx),%%zmm3	\n\t	vmovaps	0xdc0(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x580(%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0xd80(%%rcx),%%zmm5,%%zmm5	\n\t"\
	"vsubpd	0x5c0(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xdc0(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x580(%%rax)	\n\t	vmovaps	%%zmm5,0xd80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x5c0(%%rax)	\n\t	vmovaps	%%zmm7,0xdc0(%%rax)	\n\t"\
		/* xC.yC: */							/* x1C.y1C: */\
		"vmovaps	0x600(%%rax),%%zmm0	\n\t	vmovaps	0xe00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x600(%%rbx),%%zmm2	\n\t	vmovaps	0xe00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x640(%%rax),%%zmm1	\n\t	vmovaps	0xe40(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x640(%%rbx),%%zmm3	\n\t	vmovaps	0xe40(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x600(%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0xe00(%%rcx),%%zmm5,%%zmm5	\n\t"\
	"vsubpd	0x640(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xe40(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x600(%%rax)	\n\t	vmovaps	%%zmm5,0xe00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x640(%%rax)	\n\t	vmovaps	%%zmm7,0xe40(%%rax)	\n\t"\
		/* xD.yD: */							/* x1D.y1D: */\
		"vmovaps	0x680(%%rax),%%zmm0	\n\t	vmovaps	0xe80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x680(%%rbx),%%zmm2	\n\t	vmovaps	0xe80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x6c0(%%rax),%%zmm1	\n\t	vmovaps	0xec0(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x6c0(%%rbx),%%zmm3	\n\t	vmovaps	0xec0(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x680(%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0xe80(%%rcx),%%zmm5,%%zmm5	\n\t"\
	"vsubpd	0x6c0(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xec0(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x680(%%rax)	\n\t	vmovaps	%%zmm5,0xe80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x6c0(%%rax)	\n\t	vmovaps	%%zmm7,0xec0(%%rax)	\n\t"\
		/* xE.yE: */							/* x1E.y1E: */\
		"vmovaps	0x700(%%rax),%%zmm0	\n\t	vmovaps	0xf00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x700(%%rbx),%%zmm2	\n\t	vmovaps	0xf00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x740(%%rax),%%zmm1	\n\t	vmovaps	0xf40(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x740(%%rbx),%%zmm3	\n\t	vmovaps	0xf40(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x700(%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0xf00(%%rcx),%%zmm5,%%zmm5	\n\t"\
	"vsubpd	0x740(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xf40(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x700(%%rax)	\n\t	vmovaps	%%zmm5,0xf00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x740(%%rax)	\n\t	vmovaps	%%zmm7,0xf40(%%rax)	\n\t"\
		/* xF.yF: */							/* x1F.y1F: */\
		"vmovaps	0x780(%%rax),%%zmm0	\n\t	vmovaps	0xf80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x780(%%rbx),%%zmm2	\n\t	vmovaps	0xf80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x7c0(%%rax),%%zmm1	\n\t	vmovaps	0xfc0(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x7c0(%%rbx),%%zmm3	\n\t	vmovaps	0xfc0(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x780(%%rcx),%%zmm0,%%zmm0	\n\t	vsubpd	0xf80(%%rcx),%%zmm5,%%zmm5	\n\t"\
	"vsubpd	0x7c0(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xfc0(%%rcx),%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x780(%%rax)	\n\t	vmovaps	%%zmm5,0xf80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x7c0(%%rax)	\n\t	vmovaps	%%zmm7,0xfc0(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xa)	/* All inputs from memory addresses here */\
		 ,[__bdd0] "m" (Xb)	/* B-array base-address */\
		 ,[__cdd0] "m" (Xc)	/* C-array base-address */\
		: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"	/* Clobbered registers */\
	  );\
	}

#elif defined(USE_AVX)

	#define SIMD_MULSUB(Xa,Xb,Xc)\
	{\
	  __asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movq	%[__bdd0],%%rbx		\n\t"\
		"movq	%[__cdd0],%%rcx		\n\t"\
		/* x0.y0: */							/* x10.y10: */\
		"vmovaps		 (%%rax),%%ymm0	\n\t	vmovaps	0x400(%%rax),%%ymm5	\n\t"/* x.re */\
		"vmovaps		 (%%rbx),%%ymm2	\n\t	vmovaps	0x400(%%rbx),%%ymm7	\n\t"/* y.re */\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t	vmovaps	0x420(%%rax),%%ymm6	\n\t"/* x.im */\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t	vmovaps	0x420(%%rbx),%%ymm8	\n\t"/* y.im */\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"/* copy x.re */\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"/* x.re *= y.re */\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"/* y.re *= x.im */\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"/* x.im *= y.im */\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"/* y.im *= x.re[copy] */\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"/* xy.re = x.re*y.re - x.im*y.im */\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"/* xy.im = y.re*x.im + y.im*x.re */\
	"vsubpd		 (%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x400(%%rcx),%%ymm7,%%ymm7	\n\t"/* re = xy.re - z.re */\
	"vsubpd	0x020(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x420(%%rcx),%%ymm8,%%ymm8	\n\t"/* im = xy.im - z.im */\
		"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm5,0x400(%%rax)	\n\t"/* write z.re */\
		"vmovaps	%%ymm2,0x020(%%rax)	\n\t	vmovaps	%%ymm7,0x420(%%rax)	\n\t"/* write z.im */\
		/* x1.y1: */							/* x11.y11: */\
		"vmovaps	0x040(%%rax),%%ymm0	\n\t	vmovaps	0x440(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x040(%%rbx),%%ymm2	\n\t	vmovaps	0x440(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x060(%%rax),%%ymm1	\n\t	vmovaps	0x460(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x060(%%rbx),%%ymm3	\n\t	vmovaps	0x460(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x440(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x460(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	%%ymm0,0x040(%%rax)	\n\t	vmovaps	%%ymm5,0x440(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x060(%%rax)	\n\t	vmovaps	%%ymm7,0x460(%%rax)	\n\t"\
		/* x2.y2: */							/* x12.y12: */\
		"vmovaps	0x080(%%rax),%%ymm0	\n\t	vmovaps	0x480(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x080(%%rbx),%%ymm2	\n\t	vmovaps	0x480(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x0a0(%%rax),%%ymm1	\n\t	vmovaps	0x4a0(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x0a0(%%rbx),%%ymm3	\n\t	vmovaps	0x4a0(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x080(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x480(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x0a0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x4a0(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	%%ymm0,0x080(%%rax)	\n\t	vmovaps	%%ymm5,0x480(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x0a0(%%rax)	\n\t	vmovaps	%%ymm7,0x4a0(%%rax)	\n\t"\
		/* x3.y3: */							/* x13.y13: */\
		"vmovaps	0x0c0(%%rax),%%ymm0	\n\t	vmovaps	0x4c0(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x0c0(%%rbx),%%ymm2	\n\t	vmovaps	0x4c0(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x0e0(%%rax),%%ymm1	\n\t	vmovaps	0x4e0(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x0e0(%%rbx),%%ymm3	\n\t	vmovaps	0x4e0(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x0c0(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x4c0(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x0e0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x4e0(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	%%ymm0,0x0c0(%%rax)	\n\t	vmovaps	%%ymm5,0x4c0(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x0e0(%%rax)	\n\t	vmovaps	%%ymm7,0x4e0(%%rax)	\n\t"\
		/* x4.y4: */							/* x14.y14: */\
		"vmovaps	0x100(%%rax),%%ymm0	\n\t	vmovaps	0x500(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x100(%%rbx),%%ymm2	\n\t	vmovaps	0x500(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x120(%%rax),%%ymm1	\n\t	vmovaps	0x520(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x120(%%rbx),%%ymm3	\n\t	vmovaps	0x520(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x100(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x500(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x120(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x520(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	%%ymm0,0x100(%%rax)	\n\t	vmovaps	%%ymm5,0x500(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x120(%%rax)	\n\t	vmovaps	%%ymm7,0x520(%%rax)	\n\t"\
		/* x5.y5: */							/* x15.y15: */\
		"vmovaps	0x140(%%rax),%%ymm0	\n\t	vmovaps	0x540(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x140(%%rbx),%%ymm2	\n\t	vmovaps	0x540(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x160(%%rax),%%ymm1	\n\t	vmovaps	0x560(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x160(%%rbx),%%ymm3	\n\t	vmovaps	0x560(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x140(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x540(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x160(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x560(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	%%ymm0,0x140(%%rax)	\n\t	vmovaps	%%ymm5,0x540(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x160(%%rax)	\n\t	vmovaps	%%ymm7,0x560(%%rax)	\n\t"\
		/* x6.y6: */							/* x16.y16: */\
		"vmovaps	0x180(%%rax),%%ymm0	\n\t	vmovaps	0x580(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x180(%%rbx),%%ymm2	\n\t	vmovaps	0x580(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x1a0(%%rax),%%ymm1	\n\t	vmovaps	0x5a0(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x1a0(%%rbx),%%ymm3	\n\t	vmovaps	0x5a0(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x180(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x580(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x1a0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x5a0(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	%%ymm0,0x180(%%rax)	\n\t	vmovaps	%%ymm5,0x580(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x1a0(%%rax)	\n\t	vmovaps	%%ymm7,0x5a0(%%rax)	\n\t"\
		/* x7.y7: */							/* x17.y17: */\
		"vmovaps	0x1c0(%%rax),%%ymm0	\n\t	vmovaps	0x5c0(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x1c0(%%rbx),%%ymm2	\n\t	vmovaps	0x5c0(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x1e0(%%rax),%%ymm1	\n\t	vmovaps	0x5e0(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x1e0(%%rbx),%%ymm3	\n\t	vmovaps	0x5e0(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x1c0(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x5c0(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x1e0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x5e0(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	%%ymm0,0x1c0(%%rax)	\n\t	vmovaps	%%ymm5,0x5c0(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x1e0(%%rax)	\n\t	vmovaps	%%ymm7,0x5e0(%%rax)	\n\t"\
		/* x8.y8: */							/* x18.y18: */\
		"vmovaps	0x200(%%rax),%%ymm0	\n\t	vmovaps	0x600(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x200(%%rbx),%%ymm2	\n\t	vmovaps	0x600(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x220(%%rax),%%ymm1	\n\t	vmovaps	0x620(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x220(%%rbx),%%ymm3	\n\t	vmovaps	0x620(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x200(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x600(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x220(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x620(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	%%ymm0,0x200(%%rax)	\n\t	vmovaps	%%ymm5,0x600(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x220(%%rax)	\n\t	vmovaps	%%ymm7,0x620(%%rax)	\n\t"\
		/* x9.y9: */							/* x19.y19: */\
		"vmovaps	0x240(%%rax),%%ymm0	\n\t	vmovaps	0x640(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x240(%%rbx),%%ymm2	\n\t	vmovaps	0x640(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x260(%%rax),%%ymm1	\n\t	vmovaps	0x660(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x260(%%rbx),%%ymm3	\n\t	vmovaps	0x660(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x240(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x640(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x260(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x660(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	%%ymm0,0x240(%%rax)	\n\t	vmovaps	%%ymm5,0x640(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x260(%%rax)	\n\t	vmovaps	%%ymm7,0x660(%%rax)	\n\t"\
		/* xA.yA: */							/* x1A.y1A: */\
		"vmovaps	0x280(%%rax),%%ymm0	\n\t	vmovaps	0x680(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x280(%%rbx),%%ymm2	\n\t	vmovaps	0x680(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x2a0(%%rax),%%ymm1	\n\t	vmovaps	0x6a0(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x2a0(%%rbx),%%ymm3	\n\t	vmovaps	0x6a0(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x280(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x680(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x2a0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x6a0(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	%%ymm0,0x280(%%rax)	\n\t	vmovaps	%%ymm5,0x680(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x2a0(%%rax)	\n\t	vmovaps	%%ymm7,0x6a0(%%rax)	\n\t"\
		/* xB.yB: */							/* x1B.y1B: */\
		"vmovaps	0x2c0(%%rax),%%ymm0	\n\t	vmovaps	0x6c0(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x2c0(%%rbx),%%ymm2	\n\t	vmovaps	0x6c0(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x2e0(%%rax),%%ymm1	\n\t	vmovaps	0x6e0(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x2e0(%%rbx),%%ymm3	\n\t	vmovaps	0x6e0(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x2c0(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x6c0(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x2e0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x6e0(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	%%ymm0,0x2c0(%%rax)	\n\t	vmovaps	%%ymm5,0x6c0(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x2e0(%%rax)	\n\t	vmovaps	%%ymm7,0x6e0(%%rax)	\n\t"\
		/* xC.yC: */							/* x1C.y1C: */\
		"vmovaps	0x300(%%rax),%%ymm0	\n\t	vmovaps	0x700(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x300(%%rbx),%%ymm2	\n\t	vmovaps	0x700(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x320(%%rax),%%ymm1	\n\t	vmovaps	0x720(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x320(%%rbx),%%ymm3	\n\t	vmovaps	0x720(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x300(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x700(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x320(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x720(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	%%ymm0,0x300(%%rax)	\n\t	vmovaps	%%ymm5,0x700(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x320(%%rax)	\n\t	vmovaps	%%ymm7,0x720(%%rax)	\n\t"\
		/* xD.yD: */							/* x1D.y1D: */\
		"vmovaps	0x340(%%rax),%%ymm0	\n\t	vmovaps	0x740(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x340(%%rbx),%%ymm2	\n\t	vmovaps	0x740(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x360(%%rax),%%ymm1	\n\t	vmovaps	0x760(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x360(%%rbx),%%ymm3	\n\t	vmovaps	0x760(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x340(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x740(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x360(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x760(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	%%ymm0,0x340(%%rax)	\n\t	vmovaps	%%ymm5,0x740(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x360(%%rax)	\n\t	vmovaps	%%ymm7,0x760(%%rax)	\n\t"\
		/* xE.yE: */							/* x1E.y1E: */\
		"vmovaps	0x380(%%rax),%%ymm0	\n\t	vmovaps	0x780(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x380(%%rbx),%%ymm2	\n\t	vmovaps	0x780(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x3a0(%%rax),%%ymm1	\n\t	vmovaps	0x7a0(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x3a0(%%rbx),%%ymm3	\n\t	vmovaps	0x7a0(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x380(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x780(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x3a0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x7a0(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	%%ymm0,0x380(%%rax)	\n\t	vmovaps	%%ymm5,0x780(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x3a0(%%rax)	\n\t	vmovaps	%%ymm7,0x7a0(%%rax)	\n\t"\
		/* xF.yF: */							/* x1F.y1F: */\
		"vmovaps	0x3c0(%%rax),%%ymm0	\n\t	vmovaps	0x7c0(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x3c0(%%rbx),%%ymm2	\n\t	vmovaps	0x7c0(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x3e0(%%rax),%%ymm1	\n\t	vmovaps	0x7e0(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x3e0(%%rbx),%%ymm3	\n\t	vmovaps	0x7e0(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x3c0(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x7c0(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x3e0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x7e0(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	%%ymm0,0x3c0(%%rax)	\n\t	vmovaps	%%ymm5,0x7c0(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x3e0(%%rax)	\n\t	vmovaps	%%ymm7,0x7e0(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xa)	/* All inputs from memory addresses here */\
		 ,[__bdd0] "m" (Xb)	/* B-array base-address */\
		 ,[__cdd0] "m" (Xc)	/* C-array base-address */\
		: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"	/* Clobbered registers */\
	  );\
	}

#else	// 64-bit SSE2:

	// Inputs are X = [x.re,x.im], Y = [y.re,y.im]
	// Output is Z = [x.re*y.re - x.im*y.im, x.re*y.im + x.im*y.re], overwriting X
	#define SIMD_MULSUB(Xa,Xb,Xc)\
	{\
	  __asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movq	%[__bdd0],%%rbx		\n\t"\
		"movq	%[__cdd0],%%rcx		\n\t"\
		/* x0.y0: */						/* x10.y10: */\
		"movaps		 (%%rax),%%xmm0	\n\t	movaps	0x200(%%rax),%%xmm5	\n\t"/* x.re */\
		"movaps		 (%%rbx),%%xmm2	\n\t	movaps	0x200(%%rbx),%%xmm7	\n\t"/* y.re */\
		"movaps	0x010(%%rax),%%xmm1	\n\t	movaps	0x210(%%rax),%%xmm6	\n\t"/* x.im */\
		"movaps	0x010(%%rbx),%%xmm3	\n\t	movaps	0x210(%%rbx),%%xmm8	\n\t"/* y.im */\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"/* copy x.re */\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"/* x.re *= y.re */\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"/* y.re *= x.im */\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"/* x.im *= y.im */\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"/* y.im *= x.re[copy] */\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"/* xy.re = x.re*y.re - x.im*y.im */\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"/* xy.im = y.re*x.im + y.im*x.re */\
		"subpd	     (%%rcx),%%xmm2	\n\t	subpd	0x200(%%rcx),%%xmm7	\n\t"/* re = xy.re - z.re */\
		"subpd	0x010(%%rcx),%%xmm3	\n\t	subpd	0x210(%%rcx),%%xmm8	\n\t"/* im = xy.im - z.im */\
		"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm5,0x200(%%rax)	\n\t"/* write z.re */\
		"movaps	%%xmm2,0x010(%%rax)	\n\t	movaps	%%xmm7,0x210(%%rax)	\n\t"/* write z.im */\
		/* x1.y1: */						/* x11.y11: */\
		"movaps	0x020(%%rax),%%xmm0	\n\t	movaps	0x220(%%rax),%%xmm5	\n\t"\
		"movaps	0x020(%%rbx),%%xmm2	\n\t	movaps	0x220(%%rbx),%%xmm7	\n\t"\
		"movaps	0x030(%%rax),%%xmm1	\n\t	movaps	0x230(%%rax),%%xmm6	\n\t"\
		"movaps	0x030(%%rbx),%%xmm3	\n\t	movaps	0x230(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"subpd	0x020(%%rcx),%%xmm2	\n\t	subpd	0x220(%%rcx),%%xmm7	\n\t"\
		"subpd	0x030(%%rcx),%%xmm3	\n\t	subpd	0x230(%%rcx),%%xmm8	\n\t"\
		"movaps	%%xmm0,0x020(%%rax)	\n\t	movaps	%%xmm5,0x220(%%rax)	\n\t"\
		"movaps	%%xmm2,0x030(%%rax)	\n\t	movaps	%%xmm7,0x230(%%rax)	\n\t"\
		/* x2.y2: */						/* x12.y12: */\
		"movaps	0x040(%%rax),%%xmm0	\n\t	movaps	0x240(%%rax),%%xmm5	\n\t"\
		"movaps	0x040(%%rbx),%%xmm2	\n\t	movaps	0x240(%%rbx),%%xmm7	\n\t"\
		"movaps	0x050(%%rax),%%xmm1	\n\t	movaps	0x250(%%rax),%%xmm6	\n\t"\
		"movaps	0x050(%%rbx),%%xmm3	\n\t	movaps	0x250(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"subpd	0x040(%%rcx),%%xmm2	\n\t	subpd	0x240(%%rcx),%%xmm7	\n\t"\
		"subpd	0x050(%%rcx),%%xmm3	\n\t	subpd	0x250(%%rcx),%%xmm8	\n\t"\
		"movaps	%%xmm0,0x040(%%rax)	\n\t	movaps	%%xmm5,0x240(%%rax)	\n\t"\
		"movaps	%%xmm2,0x050(%%rax)	\n\t	movaps	%%xmm7,0x250(%%rax)	\n\t"\
		/* x3.y3: */						/* x13.y13: */\
		"movaps	0x060(%%rax),%%xmm0	\n\t	movaps	0x260(%%rax),%%xmm5	\n\t"\
		"movaps	0x060(%%rbx),%%xmm2	\n\t	movaps	0x260(%%rbx),%%xmm7	\n\t"\
		"movaps	0x070(%%rax),%%xmm1	\n\t	movaps	0x270(%%rax),%%xmm6	\n\t"\
		"movaps	0x070(%%rbx),%%xmm3	\n\t	movaps	0x270(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"subpd	0x060(%%rcx),%%xmm2	\n\t	subpd	0x260(%%rcx),%%xmm7	\n\t"\
		"subpd	0x070(%%rcx),%%xmm3	\n\t	subpd	0x270(%%rcx),%%xmm8	\n\t"\
		"movaps	%%xmm0,0x060(%%rax)	\n\t	movaps	%%xmm5,0x260(%%rax)	\n\t"\
		"movaps	%%xmm2,0x070(%%rax)	\n\t	movaps	%%xmm7,0x270(%%rax)	\n\t"\
		/* x4.y4: */						/* x14.y14: */\
		"movaps	0x080(%%rax),%%xmm0	\n\t	movaps	0x280(%%rax),%%xmm5	\n\t"\
		"movaps	0x080(%%rbx),%%xmm2	\n\t	movaps	0x280(%%rbx),%%xmm7	\n\t"\
		"movaps	0x090(%%rax),%%xmm1	\n\t	movaps	0x290(%%rax),%%xmm6	\n\t"\
		"movaps	0x090(%%rbx),%%xmm3	\n\t	movaps	0x290(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"subpd	0x080(%%rcx),%%xmm2	\n\t	subpd	0x280(%%rcx),%%xmm7	\n\t"\
		"subpd	0x090(%%rcx),%%xmm3	\n\t	subpd	0x290(%%rcx),%%xmm8	\n\t"\
		"movaps	%%xmm0,0x080(%%rax)	\n\t	movaps	%%xmm5,0x280(%%rax)	\n\t"\
		"movaps	%%xmm2,0x090(%%rax)	\n\t	movaps	%%xmm7,0x290(%%rax)	\n\t"\
		/* x5.y5: */						/* x15.y15: */\
		"movaps	0x0a0(%%rax),%%xmm0	\n\t	movaps	0x2a0(%%rax),%%xmm5	\n\t"\
		"movaps	0x0a0(%%rbx),%%xmm2	\n\t	movaps	0x2a0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x0b0(%%rax),%%xmm1	\n\t	movaps	0x2b0(%%rax),%%xmm6	\n\t"\
		"movaps	0x0b0(%%rbx),%%xmm3	\n\t	movaps	0x2b0(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"subpd	0x0a0(%%rcx),%%xmm2	\n\t	subpd	0x2a0(%%rcx),%%xmm7	\n\t"\
		"subpd	0x0b0(%%rcx),%%xmm3	\n\t	subpd	0x2b0(%%rcx),%%xmm8	\n\t"\
		"movaps	%%xmm0,0x0a0(%%rax)	\n\t	movaps	%%xmm5,0x2a0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x0b0(%%rax)	\n\t	movaps	%%xmm7,0x2b0(%%rax)	\n\t"\
		/* x6.y6: */						/* x16.y16: */\
		"movaps	0x0c0(%%rax),%%xmm0	\n\t	movaps	0x2c0(%%rax),%%xmm5	\n\t"\
		"movaps	0x0c0(%%rbx),%%xmm2	\n\t	movaps	0x2c0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x0d0(%%rax),%%xmm1	\n\t	movaps	0x2d0(%%rax),%%xmm6	\n\t"\
		"movaps	0x0d0(%%rbx),%%xmm3	\n\t	movaps	0x2d0(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"subpd	0x0c0(%%rcx),%%xmm2	\n\t	subpd	0x2c0(%%rcx),%%xmm7	\n\t"\
		"subpd	0x0d0(%%rcx),%%xmm3	\n\t	subpd	0x2d0(%%rcx),%%xmm8	\n\t"\
		"movaps	%%xmm0,0x0c0(%%rax)	\n\t	movaps	%%xmm5,0x2c0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x0d0(%%rax)	\n\t	movaps	%%xmm7,0x2d0(%%rax)	\n\t"\
		/* x7.y7: */						/* x17.y17: */\
		"movaps	0x0e0(%%rax),%%xmm0	\n\t	movaps	0x2e0(%%rax),%%xmm5	\n\t"\
		"movaps	0x0e0(%%rbx),%%xmm2	\n\t	movaps	0x2e0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x0f0(%%rax),%%xmm1	\n\t	movaps	0x2f0(%%rax),%%xmm6	\n\t"\
		"movaps	0x0f0(%%rbx),%%xmm3	\n\t	movaps	0x2f0(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"subpd	0x0e0(%%rcx),%%xmm2	\n\t	subpd	0x2e0(%%rcx),%%xmm7	\n\t"\
		"subpd	0x0f0(%%rcx),%%xmm3	\n\t	subpd	0x2f0(%%rcx),%%xmm8	\n\t"\
		"movaps	%%xmm0,0x0e0(%%rax)	\n\t	movaps	%%xmm5,0x2e0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x0f0(%%rax)	\n\t	movaps	%%xmm7,0x2f0(%%rax)	\n\t"\
		/* x8.y8: */						/* x18.y18: */\
		"movaps	0x100(%%rax),%%xmm0	\n\t	movaps	0x300(%%rax),%%xmm5	\n\t"\
		"movaps	0x100(%%rbx),%%xmm2	\n\t	movaps	0x300(%%rbx),%%xmm7	\n\t"\
		"movaps	0x110(%%rax),%%xmm1	\n\t	movaps	0x310(%%rax),%%xmm6	\n\t"\
		"movaps	0x110(%%rbx),%%xmm3	\n\t	movaps	0x310(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"subpd	0x100(%%rcx),%%xmm2	\n\t	subpd	0x300(%%rcx),%%xmm7	\n\t"\
		"subpd	0x110(%%rcx),%%xmm3	\n\t	subpd	0x310(%%rcx),%%xmm8	\n\t"\
		"movaps	%%xmm0,0x100(%%rax)	\n\t	movaps	%%xmm5,0x300(%%rax)	\n\t"\
		"movaps	%%xmm2,0x110(%%rax)	\n\t	movaps	%%xmm7,0x310(%%rax)	\n\t"\
		/* x9.y9: */						/* x19.y19: */\
		"movaps	0x120(%%rax),%%xmm0	\n\t	movaps	0x320(%%rax),%%xmm5	\n\t"\
		"movaps	0x120(%%rbx),%%xmm2	\n\t	movaps	0x320(%%rbx),%%xmm7	\n\t"\
		"movaps	0x130(%%rax),%%xmm1	\n\t	movaps	0x330(%%rax),%%xmm6	\n\t"\
		"movaps	0x130(%%rbx),%%xmm3	\n\t	movaps	0x330(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"subpd	0x120(%%rcx),%%xmm2	\n\t	subpd	0x320(%%rcx),%%xmm7	\n\t"\
		"subpd	0x130(%%rcx),%%xmm3	\n\t	subpd	0x330(%%rcx),%%xmm8	\n\t"\
		"movaps	%%xmm0,0x120(%%rax)	\n\t	movaps	%%xmm5,0x320(%%rax)	\n\t"\
		"movaps	%%xmm2,0x130(%%rax)	\n\t	movaps	%%xmm7,0x330(%%rax)	\n\t"\
		/* xA.yA: */						/* x1A.y1A: */\
		"movaps	0x140(%%rax),%%xmm0	\n\t	movaps	0x340(%%rax),%%xmm5	\n\t"\
		"movaps	0x140(%%rbx),%%xmm2	\n\t	movaps	0x340(%%rbx),%%xmm7	\n\t"\
		"movaps	0x150(%%rax),%%xmm1	\n\t	movaps	0x350(%%rax),%%xmm6	\n\t"\
		"movaps	0x150(%%rbx),%%xmm3	\n\t	movaps	0x350(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"subpd	0x140(%%rcx),%%xmm2	\n\t	subpd	0x340(%%rcx),%%xmm7	\n\t"\
		"subpd	0x150(%%rcx),%%xmm3	\n\t	subpd	0x350(%%rcx),%%xmm8	\n\t"\
		"movaps	%%xmm0,0x140(%%rax)	\n\t	movaps	%%xmm5,0x340(%%rax)	\n\t"\
		"movaps	%%xmm2,0x150(%%rax)	\n\t	movaps	%%xmm7,0x350(%%rax)	\n\t"\
		/* xB.yB: */						/* x1B.y1B: */\
		"movaps	0x160(%%rax),%%xmm0	\n\t	movaps	0x360(%%rax),%%xmm5	\n\t"\
		"movaps	0x160(%%rbx),%%xmm2	\n\t	movaps	0x360(%%rbx),%%xmm7	\n\t"\
		"movaps	0x170(%%rax),%%xmm1	\n\t	movaps	0x370(%%rax),%%xmm6	\n\t"\
		"movaps	0x170(%%rbx),%%xmm3	\n\t	movaps	0x370(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"subpd	0x160(%%rcx),%%xmm2	\n\t	subpd	0x360(%%rcx),%%xmm7	\n\t"\
		"subpd	0x170(%%rcx),%%xmm3	\n\t	subpd	0x370(%%rcx),%%xmm8	\n\t"\
		"movaps	%%xmm0,0x160(%%rax)	\n\t	movaps	%%xmm5,0x360(%%rax)	\n\t"\
		"movaps	%%xmm2,0x170(%%rax)	\n\t	movaps	%%xmm7,0x370(%%rax)	\n\t"\
		/* xC.yC: */						/* x1C.y1C: */\
		"movaps	0x180(%%rax),%%xmm0	\n\t	movaps	0x380(%%rax),%%xmm5	\n\t"\
		"movaps	0x180(%%rbx),%%xmm2	\n\t	movaps	0x380(%%rbx),%%xmm7	\n\t"\
		"movaps	0x190(%%rax),%%xmm1	\n\t	movaps	0x390(%%rax),%%xmm6	\n\t"\
		"movaps	0x190(%%rbx),%%xmm3	\n\t	movaps	0x390(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"subpd	0x180(%%rcx),%%xmm2	\n\t	subpd	0x380(%%rcx),%%xmm7	\n\t"\
		"subpd	0x190(%%rcx),%%xmm3	\n\t	subpd	0x390(%%rcx),%%xmm8	\n\t"\
		"movaps	%%xmm0,0x180(%%rax)	\n\t	movaps	%%xmm5,0x380(%%rax)	\n\t"\
		"movaps	%%xmm2,0x190(%%rax)	\n\t	movaps	%%xmm7,0x390(%%rax)	\n\t"\
		/* xD.yD: */						/* x1D.y1D: */\
		"movaps	0x1a0(%%rax),%%xmm0	\n\t	movaps	0x3a0(%%rax),%%xmm5	\n\t"\
		"movaps	0x1a0(%%rbx),%%xmm2	\n\t	movaps	0x3a0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x1b0(%%rax),%%xmm1	\n\t	movaps	0x3b0(%%rax),%%xmm6	\n\t"\
		"movaps	0x1b0(%%rbx),%%xmm3	\n\t	movaps	0x3b0(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"subpd	0x1a0(%%rcx),%%xmm2	\n\t	subpd	0x3a0(%%rcx),%%xmm7	\n\t"\
		"subpd	0x1b0(%%rcx),%%xmm3	\n\t	subpd	0x3b0(%%rcx),%%xmm8	\n\t"\
		"movaps	%%xmm0,0x1a0(%%rax)	\n\t	movaps	%%xmm5,0x3a0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x1b0(%%rax)	\n\t	movaps	%%xmm7,0x3b0(%%rax)	\n\t"\
		/* xE.yE: */						/* x1E.y1E: */\
		"movaps	0x1c0(%%rax),%%xmm0	\n\t	movaps	0x3c0(%%rax),%%xmm5	\n\t"\
		"movaps	0x1c0(%%rbx),%%xmm2	\n\t	movaps	0x3c0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x1d0(%%rax),%%xmm1	\n\t	movaps	0x3d0(%%rax),%%xmm6	\n\t"\
		"movaps	0x1d0(%%rbx),%%xmm3	\n\t	movaps	0x3d0(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"subpd	0x1c0(%%rcx),%%xmm2	\n\t	subpd	0x3c0(%%rcx),%%xmm7	\n\t"\
		"subpd	0x1d0(%%rcx),%%xmm3	\n\t	subpd	0x3d0(%%rcx),%%xmm8	\n\t"\
		"movaps	%%xmm0,0x1c0(%%rax)	\n\t	movaps	%%xmm5,0x3c0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x1d0(%%rax)	\n\t	movaps	%%xmm7,0x3d0(%%rax)	\n\t"\
		/* xF.yF: */						/* x1F.y1F: */\
		"movaps	0x1e0(%%rax),%%xmm0	\n\t	movaps	0x3e0(%%rax),%%xmm5	\n\t"\
		"movaps	0x1e0(%%rbx),%%xmm2	\n\t	movaps	0x3e0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x1f0(%%rax),%%xmm1	\n\t	movaps	0x3f0(%%rax),%%xmm6	\n\t"\
		"movaps	0x1f0(%%rbx),%%xmm3	\n\t	movaps	0x3f0(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"subpd	0x1e0(%%rcx),%%xmm2	\n\t	subpd	0x3e0(%%rcx),%%xmm7	\n\t"\
		"subpd	0x1f0(%%rcx),%%xmm3	\n\t	subpd	0x3f0(%%rcx),%%xmm8	\n\t"\
		"movaps	%%xmm0,0x1e0(%%rax)	\n\t	movaps	%%xmm5,0x3e0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x1f0(%%rax)	\n\t	movaps	%%xmm7,0x3f0(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xa)	/* All inputs from memory addresses here */\
		 ,[__bdd0] "m" (Xb)	/* B-array base-address */\
		 ,[__cdd0] "m" (Xc)	/* C-array base-address */\
		: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"	/* Clobbered registers */\
	  );\
	}

#endif	// AVX or SSE2?

/*********** Submul a*(b-c): ***********/

#ifdef USE_ARM_V8_SIMD

	// No Fermat-mod support on ARMv8, just supply a stub macro:
	#define SIMD_SUBMUL(Xa,Xb,Xc)\
	{\
	  __asm__ volatile (\
		"ldr x0,%[__add0]	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xa)	/* All inputs from memory addresses here */\
		 ,[__bdd0] "m" (Xb)	/* B-array base-address */\
		 ,[__cdd0] "m" (Xc)	/* C-array base-address */\
		: "cc","memory","x0"	/* Clobbered registers */\
	  );\
	}

#elif defined(USE_AVX512)	// AVX512 implements a 512-bit-register version of the the AVX2 ALL_FMA-macro

	#define SIMD_SUBMUL(Xa,Xb,Xc)\
	{\
	  __asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movq	%[__bdd0],%%rbx		\n\t"\
		"movq	%[__cdd0],%%rcx		\n\t"\
		/* x0.y0: */							/* x10.y10: */\
		"vmovaps		 (%%rbx),%%zmm2	\n\t	vmovaps	0x800(%%rbx),%%zmm7	\n\t"/* y.re */\
		"vmovaps	0x040(%%rbx),%%zmm3	\n\t	vmovaps	0x840(%%rbx),%%zmm8	\n\t"/* y.im */\
	"vsubpd		 (%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0x800(%%rcx),%%zmm7,%%zmm7	\n\t"/* y.re -= z.re */\
	"vsubpd	0x040(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0x840(%%rcx),%%zmm8,%%zmm8	\n\t"/* y.im -= z.im */\
		"vmovaps		 (%%rax),%%zmm0	\n\t	vmovaps	0x800(%%rax),%%zmm5	\n\t"/* x.re */\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t	vmovaps	0x840(%%rax),%%zmm6	\n\t"/* x.im */\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"/* copy x.re */\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"/* x.re *= y.re */\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"/* y.re *= x.im */\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"/* x.im *= y.im */\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"/* y.im *= x.re[copy] */\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"/* z.re = x.re*y.re - x.im*y.im */\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"/* z.im = y.re*x.im + y.im*x.re */\
		"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm5,0x800(%%rax)	\n\t"/* write z.re */\
		"vmovaps	%%zmm2,0x040(%%rax)	\n\t	vmovaps	%%zmm7,0x840(%%rax)	\n\t"/* write z.im */\
		/* x1.y1: */							/* x11.y11: */\
		"vmovaps	0x080(%%rbx),%%zmm2	\n\t	vmovaps	0x880(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x0c0(%%rbx),%%zmm3	\n\t	vmovaps	0x8c0(%%rbx),%%zmm8	\n\t"\
	"vsubpd	0x080(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0x880(%%rcx),%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x0c0(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0x8c0(%%rcx),%%zmm8,%%zmm8	\n\t"\
		"vmovaps	0x080(%%rax),%%zmm0	\n\t	vmovaps	0x880(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x0c0(%%rax),%%zmm1	\n\t	vmovaps	0x8c0(%%rax),%%zmm6	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x080(%%rax)	\n\t	vmovaps	%%zmm5,0x880(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x0c0(%%rax)	\n\t	vmovaps	%%zmm7,0x8c0(%%rax)	\n\t"\
		/* x2.y2: */							/* x12.y12: */\
		"vmovaps	0x100(%%rbx),%%zmm2	\n\t	vmovaps	0x900(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x140(%%rbx),%%zmm3	\n\t	vmovaps	0x940(%%rbx),%%zmm8	\n\t"\
	"vsubpd	0x100(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0x900(%%rcx),%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x140(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0x940(%%rcx),%%zmm8,%%zmm8	\n\t"\
		"vmovaps	0x100(%%rax),%%zmm0	\n\t	vmovaps	0x900(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x140(%%rax),%%zmm1	\n\t	vmovaps	0x940(%%rax),%%zmm6	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x100(%%rax)	\n\t	vmovaps	%%zmm5,0x900(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x140(%%rax)	\n\t	vmovaps	%%zmm7,0x940(%%rax)	\n\t"\
		/* x3.y3: */							/* x13.y13: */\
		"vmovaps	0x180(%%rbx),%%zmm2	\n\t	vmovaps	0x980(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x1c0(%%rbx),%%zmm3	\n\t	vmovaps	0x9c0(%%rbx),%%zmm8	\n\t"\
	"vsubpd	0x180(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0x980(%%rcx),%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x1c0(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0x9c0(%%rcx),%%zmm8,%%zmm8	\n\t"\
		"vmovaps	0x180(%%rax),%%zmm0	\n\t	vmovaps	0x980(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x1c0(%%rax),%%zmm1	\n\t	vmovaps	0x9c0(%%rax),%%zmm6	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x180(%%rax)	\n\t	vmovaps	%%zmm5,0x980(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x1c0(%%rax)	\n\t	vmovaps	%%zmm7,0x9c0(%%rax)	\n\t"\
		/* x4.y4: */							/* x14.y14: */\
		"vmovaps	0x200(%%rbx),%%zmm2	\n\t	vmovaps	0xa00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x240(%%rbx),%%zmm3	\n\t	vmovaps	0xa40(%%rbx),%%zmm8	\n\t"\
	"vsubpd	0x200(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xa00(%%rcx),%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x240(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0xa40(%%rcx),%%zmm8,%%zmm8	\n\t"\
		"vmovaps	0x200(%%rax),%%zmm0	\n\t	vmovaps	0xa00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x240(%%rax),%%zmm1	\n\t	vmovaps	0xa40(%%rax),%%zmm6	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x200(%%rax)	\n\t	vmovaps	%%zmm5,0xa00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x240(%%rax)	\n\t	vmovaps	%%zmm7,0xa40(%%rax)	\n\t"\
		/* x5.y5: */							/* x15.y15: */\
		"vmovaps	0x280(%%rbx),%%zmm2	\n\t	vmovaps	0xa80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x2c0(%%rbx),%%zmm3	\n\t	vmovaps	0xac0(%%rbx),%%zmm8	\n\t"\
	"vsubpd	0x280(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xa80(%%rcx),%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x2c0(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0xac0(%%rcx),%%zmm8,%%zmm8	\n\t"\
		"vmovaps	0x280(%%rax),%%zmm0	\n\t	vmovaps	0xa80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x2c0(%%rax),%%zmm1	\n\t	vmovaps	0xac0(%%rax),%%zmm6	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x280(%%rax)	\n\t	vmovaps	%%zmm5,0xa80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x2c0(%%rax)	\n\t	vmovaps	%%zmm7,0xac0(%%rax)	\n\t"\
		/* x6.y6: */							/* x16.y16: */\
		"vmovaps	0x300(%%rbx),%%zmm2	\n\t	vmovaps	0xb00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x340(%%rbx),%%zmm3	\n\t	vmovaps	0xb40(%%rbx),%%zmm8	\n\t"\
	"vsubpd	0x300(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xb00(%%rcx),%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x340(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0xb40(%%rcx),%%zmm8,%%zmm8	\n\t"\
		"vmovaps	0x300(%%rax),%%zmm0	\n\t	vmovaps	0xb00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x340(%%rax),%%zmm1	\n\t	vmovaps	0xb40(%%rax),%%zmm6	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x300(%%rax)	\n\t	vmovaps	%%zmm5,0xb00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x340(%%rax)	\n\t	vmovaps	%%zmm7,0xb40(%%rax)	\n\t"\
		/* x7.y7: */							/* x17.y17: */\
		"vmovaps	0x380(%%rbx),%%zmm2	\n\t	vmovaps	0xb80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x3c0(%%rbx),%%zmm3	\n\t	vmovaps	0xbc0(%%rbx),%%zmm8	\n\t"\
	"vsubpd	0x380(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xb80(%%rcx),%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x3c0(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0xbc0(%%rcx),%%zmm8,%%zmm8	\n\t"\
		"vmovaps	0x380(%%rax),%%zmm0	\n\t	vmovaps	0xb80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x3c0(%%rax),%%zmm1	\n\t	vmovaps	0xbc0(%%rax),%%zmm6	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x380(%%rax)	\n\t	vmovaps	%%zmm5,0xb80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x3c0(%%rax)	\n\t	vmovaps	%%zmm7,0xbc0(%%rax)	\n\t"\
		/* x8.y8: */							/* x18.y18: */\
		"vmovaps	0x400(%%rbx),%%zmm2	\n\t	vmovaps	0xc00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x440(%%rbx),%%zmm3	\n\t	vmovaps	0xc40(%%rbx),%%zmm8	\n\t"\
	"vsubpd	0x400(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xc00(%%rcx),%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x440(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0xc40(%%rcx),%%zmm8,%%zmm8	\n\t"\
		"vmovaps	0x400(%%rax),%%zmm0	\n\t	vmovaps	0xc00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x440(%%rax),%%zmm1	\n\t	vmovaps	0xc40(%%rax),%%zmm6	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x400(%%rax)	\n\t	vmovaps	%%zmm5,0xc00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x440(%%rax)	\n\t	vmovaps	%%zmm7,0xc40(%%rax)	\n\t"\
		/* x9.y9: */							/* x19.y19: */\
		"vmovaps	0x480(%%rbx),%%zmm2	\n\t	vmovaps	0xc80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x4c0(%%rbx),%%zmm3	\n\t	vmovaps	0xcc0(%%rbx),%%zmm8	\n\t"\
	"vsubpd	0x480(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xc80(%%rcx),%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x4c0(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0xcc0(%%rcx),%%zmm8,%%zmm8	\n\t"\
		"vmovaps	0x480(%%rax),%%zmm0	\n\t	vmovaps	0xc80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x4c0(%%rax),%%zmm1	\n\t	vmovaps	0xcc0(%%rax),%%zmm6	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x480(%%rax)	\n\t	vmovaps	%%zmm5,0xc80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x4c0(%%rax)	\n\t	vmovaps	%%zmm7,0xcc0(%%rax)	\n\t"\
		/* xA.yA: */							/* x1A.y1A: */\
		"vmovaps	0x500(%%rbx),%%zmm2	\n\t	vmovaps	0xd00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x540(%%rbx),%%zmm3	\n\t	vmovaps	0xd40(%%rbx),%%zmm8	\n\t"\
	"vsubpd	0x500(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xd00(%%rcx),%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x540(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0xd40(%%rcx),%%zmm8,%%zmm8	\n\t"\
		"vmovaps	0x500(%%rax),%%zmm0	\n\t	vmovaps	0xd00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x540(%%rax),%%zmm1	\n\t	vmovaps	0xd40(%%rax),%%zmm6	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x500(%%rax)	\n\t	vmovaps	%%zmm5,0xd00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x540(%%rax)	\n\t	vmovaps	%%zmm7,0xd40(%%rax)	\n\t"\
		/* xB.yB: */							/* x1B.y1B: */\
		"vmovaps	0x580(%%rbx),%%zmm2	\n\t	vmovaps	0xd80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x5c0(%%rbx),%%zmm3	\n\t	vmovaps	0xdc0(%%rbx),%%zmm8	\n\t"\
	"vsubpd	0x580(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xd80(%%rcx),%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x5c0(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0xdc0(%%rcx),%%zmm8,%%zmm8	\n\t"\
		"vmovaps	0x580(%%rax),%%zmm0	\n\t	vmovaps	0xd80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x5c0(%%rax),%%zmm1	\n\t	vmovaps	0xdc0(%%rax),%%zmm6	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x580(%%rax)	\n\t	vmovaps	%%zmm5,0xd80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x5c0(%%rax)	\n\t	vmovaps	%%zmm7,0xdc0(%%rax)	\n\t"\
		/* xC.yC: */							/* x1C.y1C: */\
		"vmovaps	0x600(%%rbx),%%zmm2	\n\t	vmovaps	0xe00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x640(%%rbx),%%zmm3	\n\t	vmovaps	0xe40(%%rbx),%%zmm8	\n\t"\
	"vsubpd	0x600(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xe00(%%rcx),%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x640(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0xe40(%%rcx),%%zmm8,%%zmm8	\n\t"\
		"vmovaps	0x600(%%rax),%%zmm0	\n\t	vmovaps	0xe00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x640(%%rax),%%zmm1	\n\t	vmovaps	0xe40(%%rax),%%zmm6	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x600(%%rax)	\n\t	vmovaps	%%zmm5,0xe00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x640(%%rax)	\n\t	vmovaps	%%zmm7,0xe40(%%rax)	\n\t"\
		/* xD.yD: */							/* x1D.y1D: */\
		"vmovaps	0x680(%%rbx),%%zmm2	\n\t	vmovaps	0xe80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x6c0(%%rbx),%%zmm3	\n\t	vmovaps	0xec0(%%rbx),%%zmm8	\n\t"\
	"vsubpd	0x680(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xe80(%%rcx),%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x6c0(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0xec0(%%rcx),%%zmm8,%%zmm8	\n\t"\
		"vmovaps	0x680(%%rax),%%zmm0	\n\t	vmovaps	0xe80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x6c0(%%rax),%%zmm1	\n\t	vmovaps	0xec0(%%rax),%%zmm6	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x680(%%rax)	\n\t	vmovaps	%%zmm5,0xe80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x6c0(%%rax)	\n\t	vmovaps	%%zmm7,0xec0(%%rax)	\n\t"\
		/* xE.yE: */							/* x1E.y1E: */\
		"vmovaps	0x700(%%rbx),%%zmm2	\n\t	vmovaps	0xf00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x740(%%rbx),%%zmm3	\n\t	vmovaps	0xf40(%%rbx),%%zmm8	\n\t"\
	"vsubpd	0x700(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xf00(%%rcx),%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x740(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0xf40(%%rcx),%%zmm8,%%zmm8	\n\t"\
		"vmovaps	0x700(%%rax),%%zmm0	\n\t	vmovaps	0xf00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x740(%%rax),%%zmm1	\n\t	vmovaps	0xf40(%%rax),%%zmm6	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x700(%%rax)	\n\t	vmovaps	%%zmm5,0xf00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x740(%%rax)	\n\t	vmovaps	%%zmm7,0xf40(%%rax)	\n\t"\
		/* xF.yF: */							/* x1F.y1F: */\
		"vmovaps	0x780(%%rbx),%%zmm2	\n\t	vmovaps	0xf80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x7c0(%%rbx),%%zmm3	\n\t	vmovaps	0xfc0(%%rbx),%%zmm8	\n\t"\
	"vsubpd	0x780(%%rcx),%%zmm2,%%zmm2	\n\t	vsubpd	0xf80(%%rcx),%%zmm7,%%zmm7	\n\t"\
	"vsubpd	0x7c0(%%rcx),%%zmm3,%%zmm3	\n\t	vsubpd	0xfc0(%%rcx),%%zmm8,%%zmm8	\n\t"\
		"vmovaps	0x780(%%rax),%%zmm0	\n\t	vmovaps	0xf80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x7c0(%%rax),%%zmm1	\n\t	vmovaps	0xfc0(%%rax),%%zmm6	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x780(%%rax)	\n\t	vmovaps	%%zmm5,0xf80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x7c0(%%rax)	\n\t	vmovaps	%%zmm7,0xfc0(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xa)	/* All inputs from memory addresses here */\
		 ,[__bdd0] "m" (Xb)	/* B-array base-address */\
		 ,[__cdd0] "m" (Xc)	/* C-array base-address */\
		: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"	/* Clobbered registers */\
	  );\
	}

#elif defined(USE_AVX)

	#define SIMD_SUBMUL(Xa,Xb,Xc)\
	{\
	  __asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movq	%[__bdd0],%%rbx		\n\t"\
		"movq	%[__cdd0],%%rcx		\n\t"\
		/* x0.y0: */							/* x10.y10: */\
		"vmovaps		 (%%rbx),%%ymm2	\n\t	vmovaps	0x400(%%rbx),%%ymm7	\n\t"/* y.re */\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t	vmovaps	0x420(%%rbx),%%ymm8	\n\t"/* y.im */\
	"vsubpd	     (%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x400(%%rcx),%%ymm7,%%ymm7	\n\t"/* y.re -= z.re */\
	"vsubpd	0x020(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x420(%%rcx),%%ymm8,%%ymm8	\n\t"/* y.im -= z.im */\
		"vmovaps		 (%%rax),%%ymm0	\n\t	vmovaps	0x400(%%rax),%%ymm5	\n\t"/* x.re */\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t	vmovaps	0x420(%%rax),%%ymm6	\n\t"/* x.im */\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"/* copy x.re */\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"/* x.re *= y.re */\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"/* y.re *= x.im */\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"/* x.im *= y.im */\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"/* y.im *= x.re[copy] */\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"/* z.re = x.re*y.re - x.im*y.im */\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"/* z.im = y.re*x.im + y.im*x.re */\
		"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm5,0x400(%%rax)	\n\t"/* write z.re */\
		"vmovaps	%%ymm2,0x020(%%rax)	\n\t	vmovaps	%%ymm7,0x420(%%rax)	\n\t"/* write z.im */\
		/* x1.y1: */							/* x11.y11: */\
		"vmovaps	0x040(%%rbx),%%ymm2	\n\t	vmovaps	0x440(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x060(%%rbx),%%ymm3	\n\t	vmovaps	0x460(%%rbx),%%ymm8	\n\t"\
	"vsubpd	0x040(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x440(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x060(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x460(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	0x040(%%rax),%%ymm0	\n\t	vmovaps	0x440(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x060(%%rax),%%ymm1	\n\t	vmovaps	0x460(%%rax),%%ymm6	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x040(%%rax)	\n\t	vmovaps	%%ymm5,0x440(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x060(%%rax)	\n\t	vmovaps	%%ymm7,0x460(%%rax)	\n\t"\
		/* x2.y2: */							/* x12.y12: */\
		"vmovaps	0x080(%%rbx),%%ymm2	\n\t	vmovaps	0x480(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x0a0(%%rbx),%%ymm3	\n\t	vmovaps	0x4a0(%%rbx),%%ymm8	\n\t"\
	"vsubpd	0x080(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x480(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x0a0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x4a0(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	0x080(%%rax),%%ymm0	\n\t	vmovaps	0x480(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x0a0(%%rax),%%ymm1	\n\t	vmovaps	0x4a0(%%rax),%%ymm6	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x080(%%rax)	\n\t	vmovaps	%%ymm5,0x480(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x0a0(%%rax)	\n\t	vmovaps	%%ymm7,0x4a0(%%rax)	\n\t"\
		/* x3.y3: */							/* x13.y13: */\
		"vmovaps	0x0c0(%%rbx),%%ymm2	\n\t	vmovaps	0x4c0(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x0e0(%%rbx),%%ymm3	\n\t	vmovaps	0x4e0(%%rbx),%%ymm8	\n\t"\
	"vsubpd	0x0c0(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x4c0(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x0e0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x4e0(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	0x0c0(%%rax),%%ymm0	\n\t	vmovaps	0x4c0(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x0e0(%%rax),%%ymm1	\n\t	vmovaps	0x4e0(%%rax),%%ymm6	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x0c0(%%rax)	\n\t	vmovaps	%%ymm5,0x4c0(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x0e0(%%rax)	\n\t	vmovaps	%%ymm7,0x4e0(%%rax)	\n\t"\
		/* x4.y4: */							/* x14.y14: */\
		"vmovaps	0x100(%%rbx),%%ymm2	\n\t	vmovaps	0x500(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x120(%%rbx),%%ymm3	\n\t	vmovaps	0x520(%%rbx),%%ymm8	\n\t"\
	"vsubpd	0x100(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x500(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x120(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x520(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	0x100(%%rax),%%ymm0	\n\t	vmovaps	0x500(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x120(%%rax),%%ymm1	\n\t	vmovaps	0x520(%%rax),%%ymm6	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x100(%%rax)	\n\t	vmovaps	%%ymm5,0x500(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x120(%%rax)	\n\t	vmovaps	%%ymm7,0x520(%%rax)	\n\t"\
		/* x5.y5: */							/* x15.y15: */\
		"vmovaps	0x140(%%rbx),%%ymm2	\n\t	vmovaps	0x540(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x160(%%rbx),%%ymm3	\n\t	vmovaps	0x560(%%rbx),%%ymm8	\n\t"\
	"vsubpd	0x140(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x540(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x160(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x560(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	0x140(%%rax),%%ymm0	\n\t	vmovaps	0x540(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x160(%%rax),%%ymm1	\n\t	vmovaps	0x560(%%rax),%%ymm6	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x140(%%rax)	\n\t	vmovaps	%%ymm5,0x540(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x160(%%rax)	\n\t	vmovaps	%%ymm7,0x560(%%rax)	\n\t"\
		/* x6.y6: */							/* x16.y16: */\
		"vmovaps	0x180(%%rbx),%%ymm2	\n\t	vmovaps	0x580(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x1a0(%%rbx),%%ymm3	\n\t	vmovaps	0x5a0(%%rbx),%%ymm8	\n\t"\
	"vsubpd	0x180(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x580(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x1a0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x5a0(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	0x180(%%rax),%%ymm0	\n\t	vmovaps	0x580(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x1a0(%%rax),%%ymm1	\n\t	vmovaps	0x5a0(%%rax),%%ymm6	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x180(%%rax)	\n\t	vmovaps	%%ymm5,0x580(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x1a0(%%rax)	\n\t	vmovaps	%%ymm7,0x5a0(%%rax)	\n\t"\
		/* x7.y7: */							/* x17.y17: */\
		"vmovaps	0x1c0(%%rbx),%%ymm2	\n\t	vmovaps	0x5c0(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x1e0(%%rbx),%%ymm3	\n\t	vmovaps	0x5e0(%%rbx),%%ymm8	\n\t"\
	"vsubpd	0x1c0(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x5c0(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x1e0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x5e0(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	0x1c0(%%rax),%%ymm0	\n\t	vmovaps	0x5c0(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x1e0(%%rax),%%ymm1	\n\t	vmovaps	0x5e0(%%rax),%%ymm6	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x1c0(%%rax)	\n\t	vmovaps	%%ymm5,0x5c0(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x1e0(%%rax)	\n\t	vmovaps	%%ymm7,0x5e0(%%rax)	\n\t"\
		/* x8.y8: */							/* x18.y18: */\
		"vmovaps	0x200(%%rbx),%%ymm2	\n\t	vmovaps	0x600(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x220(%%rbx),%%ymm3	\n\t	vmovaps	0x620(%%rbx),%%ymm8	\n\t"\
	"vsubpd	0x200(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x600(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x220(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x620(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	0x200(%%rax),%%ymm0	\n\t	vmovaps	0x600(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x220(%%rax),%%ymm1	\n\t	vmovaps	0x620(%%rax),%%ymm6	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x200(%%rax)	\n\t	vmovaps	%%ymm5,0x600(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x220(%%rax)	\n\t	vmovaps	%%ymm7,0x620(%%rax)	\n\t"\
		/* x9.y9: */							/* x19.y19: */\
		"vmovaps	0x240(%%rbx),%%ymm2	\n\t	vmovaps	0x640(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x260(%%rbx),%%ymm3	\n\t	vmovaps	0x660(%%rbx),%%ymm8	\n\t"\
	"vsubpd	0x240(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x640(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x260(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x660(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	0x240(%%rax),%%ymm0	\n\t	vmovaps	0x640(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x260(%%rax),%%ymm1	\n\t	vmovaps	0x660(%%rax),%%ymm6	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x240(%%rax)	\n\t	vmovaps	%%ymm5,0x640(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x260(%%rax)	\n\t	vmovaps	%%ymm7,0x660(%%rax)	\n\t"\
		/* xA.yA: */							/* x1A.y1A: */\
		"vmovaps	0x280(%%rbx),%%ymm2	\n\t	vmovaps	0x680(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x2a0(%%rbx),%%ymm3	\n\t	vmovaps	0x6a0(%%rbx),%%ymm8	\n\t"\
	"vsubpd	0x280(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x680(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x2a0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x6a0(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	0x280(%%rax),%%ymm0	\n\t	vmovaps	0x680(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x2a0(%%rax),%%ymm1	\n\t	vmovaps	0x6a0(%%rax),%%ymm6	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x280(%%rax)	\n\t	vmovaps	%%ymm5,0x680(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x2a0(%%rax)	\n\t	vmovaps	%%ymm7,0x6a0(%%rax)	\n\t"\
		/* xB.yB: */							/* x1B.y1B: */\
		"vmovaps	0x2c0(%%rbx),%%ymm2	\n\t	vmovaps	0x6c0(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x2e0(%%rbx),%%ymm3	\n\t	vmovaps	0x6e0(%%rbx),%%ymm8	\n\t"\
	"vsubpd	0x2c0(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x6c0(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x2e0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x6e0(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	0x2c0(%%rax),%%ymm0	\n\t	vmovaps	0x6c0(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x2e0(%%rax),%%ymm1	\n\t	vmovaps	0x6e0(%%rax),%%ymm6	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x2c0(%%rax)	\n\t	vmovaps	%%ymm5,0x6c0(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x2e0(%%rax)	\n\t	vmovaps	%%ymm7,0x6e0(%%rax)	\n\t"\
		/* xC.yC: */							/* x1C.y1C: */\
		"vmovaps	0x300(%%rbx),%%ymm2	\n\t	vmovaps	0x700(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x320(%%rbx),%%ymm3	\n\t	vmovaps	0x720(%%rbx),%%ymm8	\n\t"\
	"vsubpd	0x300(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x700(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x320(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x720(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	0x300(%%rax),%%ymm0	\n\t	vmovaps	0x700(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x320(%%rax),%%ymm1	\n\t	vmovaps	0x720(%%rax),%%ymm6	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x300(%%rax)	\n\t	vmovaps	%%ymm5,0x700(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x320(%%rax)	\n\t	vmovaps	%%ymm7,0x720(%%rax)	\n\t"\
		/* xD.yD: */							/* x1D.y1D: */\
		"vmovaps	0x340(%%rbx),%%ymm2	\n\t	vmovaps	0x740(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x360(%%rbx),%%ymm3	\n\t	vmovaps	0x760(%%rbx),%%ymm8	\n\t"\
	"vsubpd	0x340(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x740(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x360(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x760(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	0x340(%%rax),%%ymm0	\n\t	vmovaps	0x740(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x360(%%rax),%%ymm1	\n\t	vmovaps	0x760(%%rax),%%ymm6	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x340(%%rax)	\n\t	vmovaps	%%ymm5,0x740(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x360(%%rax)	\n\t	vmovaps	%%ymm7,0x760(%%rax)	\n\t"\
		/* xE.yE: */							/* x1E.y1E: */\
		"vmovaps	0x380(%%rbx),%%ymm2	\n\t	vmovaps	0x780(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x3a0(%%rbx),%%ymm3	\n\t	vmovaps	0x7a0(%%rbx),%%ymm8	\n\t"\
	"vsubpd	0x380(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x780(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x3a0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x7a0(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	0x380(%%rax),%%ymm0	\n\t	vmovaps	0x780(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x3a0(%%rax),%%ymm1	\n\t	vmovaps	0x7a0(%%rax),%%ymm6	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x380(%%rax)	\n\t	vmovaps	%%ymm5,0x780(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x3a0(%%rax)	\n\t	vmovaps	%%ymm7,0x7a0(%%rax)	\n\t"\
		/* xF.yF: */							/* x1F.y1F: */\
		"vmovaps	0x3c0(%%rbx),%%ymm2	\n\t	vmovaps	0x7c0(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x3e0(%%rbx),%%ymm3	\n\t	vmovaps	0x7e0(%%rbx),%%ymm8	\n\t"\
	"vsubpd	0x3c0(%%rcx),%%ymm2,%%ymm2	\n\t	vsubpd	0x7c0(%%rcx),%%ymm7,%%ymm7	\n\t"\
	"vsubpd	0x3e0(%%rcx),%%ymm3,%%ymm3	\n\t	vsubpd	0x7e0(%%rcx),%%ymm8,%%ymm8	\n\t"\
		"vmovaps	0x3c0(%%rax),%%ymm0	\n\t	vmovaps	0x7c0(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x3e0(%%rax),%%ymm1	\n\t	vmovaps	0x7e0(%%rax),%%ymm6	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x3c0(%%rax)	\n\t	vmovaps	%%ymm5,0x7c0(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x3e0(%%rax)	\n\t	vmovaps	%%ymm7,0x7e0(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xa)	/* All inputs from memory addresses here */\
		 ,[__bdd0] "m" (Xb)	/* B-array base-address */\
		 ,[__cdd0] "m" (Xc)	/* C-array base-address */\
		: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"	/* Clobbered registers */\
	  );\
	}

#else	// 64-bit SSE2:

	// Inputs are X = [x.re,x.im], Y = [y.re,y.im]
	// Output is Z = [x.re*y.re - x.im*y.im, x.re*y.im + x.im*y.re], overwriting X
	#define SIMD_SUBMUL(Xa,Xb,Xc)\
	{\
	  __asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movq	%[__bdd0],%%rbx		\n\t"\
		"movq	%[__cdd0],%%rcx		\n\t"\
		/* x0.y0: */						/* x10.y10: */\
		"movaps		 (%%rbx),%%xmm2	\n\t	movaps	0x200(%%rbx),%%xmm7	\n\t"/* y.re */\
		"movaps	0x010(%%rbx),%%xmm3	\n\t	movaps	0x210(%%rbx),%%xmm8	\n\t"/* y.im */\
		"subpd	     (%%rcx),%%xmm2	\n\t	subpd	0x200(%%rcx),%%xmm7	\n\t"/* y.re -= z.re */\
		"subpd	0x010(%%rcx),%%xmm3	\n\t	subpd	0x210(%%rcx),%%xmm8	\n\t"/* y.im -= z.im */\
		"movaps		 (%%rax),%%xmm0	\n\t	movaps	0x200(%%rax),%%xmm5	\n\t"/* x.re */\
		"movaps	0x010(%%rax),%%xmm1	\n\t	movaps	0x210(%%rax),%%xmm6	\n\t"/* x.im */\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"/* copy x.re */\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"/* x.re *= y.re */\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"/* y.re *= x.im */\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"/* x.im *= y.im */\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"/* y.im *= x.re[copy] */\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"/* z.re = x.re*y.re - x.im*y.im */\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"/* z.im = y.re*x.im + y.im*x.re */\
		"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm5,0x200(%%rax)	\n\t"/* write z.re */\
		"movaps	%%xmm2,0x010(%%rax)	\n\t	movaps	%%xmm7,0x210(%%rax)	\n\t"/* write z.im */\
		/* x1.y1: */						/* x11.y11: */\
		"movaps	0x020(%%rbx),%%xmm2	\n\t	movaps	0x220(%%rbx),%%xmm7	\n\t"\
		"movaps	0x030(%%rbx),%%xmm3	\n\t	movaps	0x230(%%rbx),%%xmm8	\n\t"\
		"subpd	0x020(%%rcx),%%xmm2	\n\t	subpd	0x220(%%rcx),%%xmm7	\n\t"\
		"subpd	0x030(%%rcx),%%xmm3	\n\t	subpd	0x230(%%rcx),%%xmm8	\n\t"\
		"movaps	0x020(%%rax),%%xmm0	\n\t	movaps	0x220(%%rax),%%xmm5	\n\t"\
		"movaps	0x030(%%rax),%%xmm1	\n\t	movaps	0x230(%%rax),%%xmm6	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x020(%%rax)	\n\t	movaps	%%xmm5,0x220(%%rax)	\n\t"\
		"movaps	%%xmm2,0x030(%%rax)	\n\t	movaps	%%xmm7,0x230(%%rax)	\n\t"\
		/* x2.y2: */						/* x12.y12: */\
		"movaps	0x040(%%rbx),%%xmm2	\n\t	movaps	0x240(%%rbx),%%xmm7	\n\t"\
		"movaps	0x050(%%rbx),%%xmm3	\n\t	movaps	0x250(%%rbx),%%xmm8	\n\t"\
		"subpd	0x040(%%rcx),%%xmm2	\n\t	subpd	0x240(%%rcx),%%xmm7	\n\t"\
		"subpd	0x050(%%rcx),%%xmm3	\n\t	subpd	0x250(%%rcx),%%xmm8	\n\t"\
		"movaps	0x040(%%rax),%%xmm0	\n\t	movaps	0x240(%%rax),%%xmm5	\n\t"\
		"movaps	0x050(%%rax),%%xmm1	\n\t	movaps	0x250(%%rax),%%xmm6	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x040(%%rax)	\n\t	movaps	%%xmm5,0x240(%%rax)	\n\t"\
		"movaps	%%xmm2,0x050(%%rax)	\n\t	movaps	%%xmm7,0x250(%%rax)	\n\t"\
		/* x3.y3: */						/* x13.y13: */\
		"movaps	0x060(%%rbx),%%xmm2	\n\t	movaps	0x260(%%rbx),%%xmm7	\n\t"\
		"movaps	0x070(%%rbx),%%xmm3	\n\t	movaps	0x270(%%rbx),%%xmm8	\n\t"\
		"subpd	0x060(%%rcx),%%xmm2	\n\t	subpd	0x260(%%rcx),%%xmm7	\n\t"\
		"subpd	0x070(%%rcx),%%xmm3	\n\t	subpd	0x270(%%rcx),%%xmm8	\n\t"\
		"movaps	0x060(%%rax),%%xmm0	\n\t	movaps	0x260(%%rax),%%xmm5	\n\t"\
		"movaps	0x070(%%rax),%%xmm1	\n\t	movaps	0x270(%%rax),%%xmm6	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x060(%%rax)	\n\t	movaps	%%xmm5,0x260(%%rax)	\n\t"\
		"movaps	%%xmm2,0x070(%%rax)	\n\t	movaps	%%xmm7,0x270(%%rax)	\n\t"\
		/* x4.y4: */						/* x14.y14: */\
		"movaps	0x080(%%rbx),%%xmm2	\n\t	movaps	0x280(%%rbx),%%xmm7	\n\t"\
		"movaps	0x090(%%rbx),%%xmm3	\n\t	movaps	0x290(%%rbx),%%xmm8	\n\t"\
		"subpd	0x080(%%rcx),%%xmm2	\n\t	subpd	0x280(%%rcx),%%xmm7	\n\t"\
		"subpd	0x090(%%rcx),%%xmm3	\n\t	subpd	0x290(%%rcx),%%xmm8	\n\t"\
		"movaps	0x080(%%rax),%%xmm0	\n\t	movaps	0x280(%%rax),%%xmm5	\n\t"\
		"movaps	0x090(%%rax),%%xmm1	\n\t	movaps	0x290(%%rax),%%xmm6	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x080(%%rax)	\n\t	movaps	%%xmm5,0x280(%%rax)	\n\t"\
		"movaps	%%xmm2,0x090(%%rax)	\n\t	movaps	%%xmm7,0x290(%%rax)	\n\t"\
		/* x5.y5: */						/* x15.y15: */\
		"movaps	0x0a0(%%rbx),%%xmm2	\n\t	movaps	0x2a0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x0b0(%%rbx),%%xmm3	\n\t	movaps	0x2b0(%%rbx),%%xmm8	\n\t"\
		"subpd	0x0a0(%%rcx),%%xmm2	\n\t	subpd	0x2a0(%%rcx),%%xmm7	\n\t"\
		"subpd	0x0b0(%%rcx),%%xmm3	\n\t	subpd	0x2b0(%%rcx),%%xmm8	\n\t"\
		"movaps	0x0a0(%%rax),%%xmm0	\n\t	movaps	0x2a0(%%rax),%%xmm5	\n\t"\
		"movaps	0x0b0(%%rax),%%xmm1	\n\t	movaps	0x2b0(%%rax),%%xmm6	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x0a0(%%rax)	\n\t	movaps	%%xmm5,0x2a0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x0b0(%%rax)	\n\t	movaps	%%xmm7,0x2b0(%%rax)	\n\t"\
		/* x6.y6: */						/* x16.y16: */\
		"movaps	0x0c0(%%rbx),%%xmm2	\n\t	movaps	0x2c0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x0d0(%%rbx),%%xmm3	\n\t	movaps	0x2d0(%%rbx),%%xmm8	\n\t"\
		"subpd	0x0c0(%%rcx),%%xmm2	\n\t	subpd	0x2c0(%%rcx),%%xmm7	\n\t"\
		"subpd	0x0d0(%%rcx),%%xmm3	\n\t	subpd	0x2d0(%%rcx),%%xmm8	\n\t"\
		"movaps	0x0c0(%%rax),%%xmm0	\n\t	movaps	0x2c0(%%rax),%%xmm5	\n\t"\
		"movaps	0x0d0(%%rax),%%xmm1	\n\t	movaps	0x2d0(%%rax),%%xmm6	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x0c0(%%rax)	\n\t	movaps	%%xmm5,0x2c0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x0d0(%%rax)	\n\t	movaps	%%xmm7,0x2d0(%%rax)	\n\t"\
		/* x7.y7: */						/* x17.y17: */\
		"movaps	0x0e0(%%rbx),%%xmm2	\n\t	movaps	0x2e0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x0f0(%%rbx),%%xmm3	\n\t	movaps	0x2f0(%%rbx),%%xmm8	\n\t"\
		"subpd	0x0e0(%%rcx),%%xmm2	\n\t	subpd	0x2e0(%%rcx),%%xmm7	\n\t"\
		"subpd	0x0f0(%%rcx),%%xmm3	\n\t	subpd	0x2f0(%%rcx),%%xmm8	\n\t"\
		"movaps	0x0e0(%%rax),%%xmm0	\n\t	movaps	0x2e0(%%rax),%%xmm5	\n\t"\
		"movaps	0x0f0(%%rax),%%xmm1	\n\t	movaps	0x2f0(%%rax),%%xmm6	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x0e0(%%rax)	\n\t	movaps	%%xmm5,0x2e0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x0f0(%%rax)	\n\t	movaps	%%xmm7,0x2f0(%%rax)	\n\t"\
		/* x8.y8: */						/* x18.y18: */\
		"movaps	0x100(%%rbx),%%xmm2	\n\t	movaps	0x300(%%rbx),%%xmm7	\n\t"\
		"movaps	0x110(%%rbx),%%xmm3	\n\t	movaps	0x310(%%rbx),%%xmm8	\n\t"\
		"subpd	0x100(%%rcx),%%xmm2	\n\t	subpd	0x300(%%rcx),%%xmm7	\n\t"\
		"subpd	0x110(%%rcx),%%xmm3	\n\t	subpd	0x310(%%rcx),%%xmm8	\n\t"\
		"movaps	0x100(%%rax),%%xmm0	\n\t	movaps	0x300(%%rax),%%xmm5	\n\t"\
		"movaps	0x110(%%rax),%%xmm1	\n\t	movaps	0x310(%%rax),%%xmm6	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x100(%%rax)	\n\t	movaps	%%xmm5,0x300(%%rax)	\n\t"\
		"movaps	%%xmm2,0x110(%%rax)	\n\t	movaps	%%xmm7,0x310(%%rax)	\n\t"\
		/* x9.y9: */						/* x19.y19: */\
		"movaps	0x120(%%rbx),%%xmm2	\n\t	movaps	0x320(%%rbx),%%xmm7	\n\t"\
		"movaps	0x130(%%rbx),%%xmm3	\n\t	movaps	0x330(%%rbx),%%xmm8	\n\t"\
		"subpd	0x120(%%rcx),%%xmm2	\n\t	subpd	0x320(%%rcx),%%xmm7	\n\t"\
		"subpd	0x130(%%rcx),%%xmm3	\n\t	subpd	0x330(%%rcx),%%xmm8	\n\t"\
		"movaps	0x120(%%rax),%%xmm0	\n\t	movaps	0x320(%%rax),%%xmm5	\n\t"\
		"movaps	0x130(%%rax),%%xmm1	\n\t	movaps	0x330(%%rax),%%xmm6	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x120(%%rax)	\n\t	movaps	%%xmm5,0x320(%%rax)	\n\t"\
		"movaps	%%xmm2,0x130(%%rax)	\n\t	movaps	%%xmm7,0x330(%%rax)	\n\t"\
		/* xA.yA: */						/* x1A.y1A: */\
		"movaps	0x140(%%rbx),%%xmm2	\n\t	movaps	0x340(%%rbx),%%xmm7	\n\t"\
		"movaps	0x150(%%rbx),%%xmm3	\n\t	movaps	0x350(%%rbx),%%xmm8	\n\t"\
		"subpd	0x140(%%rcx),%%xmm2	\n\t	subpd	0x340(%%rcx),%%xmm7	\n\t"\
		"subpd	0x150(%%rcx),%%xmm3	\n\t	subpd	0x350(%%rcx),%%xmm8	\n\t"\
		"movaps	0x140(%%rax),%%xmm0	\n\t	movaps	0x340(%%rax),%%xmm5	\n\t"\
		"movaps	0x150(%%rax),%%xmm1	\n\t	movaps	0x350(%%rax),%%xmm6	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x140(%%rax)	\n\t	movaps	%%xmm5,0x340(%%rax)	\n\t"\
		"movaps	%%xmm2,0x150(%%rax)	\n\t	movaps	%%xmm7,0x350(%%rax)	\n\t"\
		/* xB.yB: */						/* x1B.y1B: */\
		"movaps	0x160(%%rbx),%%xmm2	\n\t	movaps	0x360(%%rbx),%%xmm7	\n\t"\
		"movaps	0x170(%%rbx),%%xmm3	\n\t	movaps	0x370(%%rbx),%%xmm8	\n\t"\
		"subpd	0x160(%%rcx),%%xmm2	\n\t	subpd	0x360(%%rcx),%%xmm7	\n\t"\
		"subpd	0x170(%%rcx),%%xmm3	\n\t	subpd	0x370(%%rcx),%%xmm8	\n\t"\
		"movaps	0x160(%%rax),%%xmm0	\n\t	movaps	0x360(%%rax),%%xmm5	\n\t"\
		"movaps	0x170(%%rax),%%xmm1	\n\t	movaps	0x370(%%rax),%%xmm6	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x160(%%rax)	\n\t	movaps	%%xmm5,0x360(%%rax)	\n\t"\
		"movaps	%%xmm2,0x170(%%rax)	\n\t	movaps	%%xmm7,0x370(%%rax)	\n\t"\
		/* xC.yC: */						/* x1C.y1C: */\
		"movaps	0x180(%%rbx),%%xmm2	\n\t	movaps	0x380(%%rbx),%%xmm7	\n\t"\
		"movaps	0x190(%%rbx),%%xmm3	\n\t	movaps	0x390(%%rbx),%%xmm8	\n\t"\
		"subpd	0x180(%%rcx),%%xmm2	\n\t	subpd	0x380(%%rcx),%%xmm7	\n\t"\
		"subpd	0x190(%%rcx),%%xmm3	\n\t	subpd	0x390(%%rcx),%%xmm8	\n\t"\
		"movaps	0x180(%%rax),%%xmm0	\n\t	movaps	0x380(%%rax),%%xmm5	\n\t"\
		"movaps	0x190(%%rax),%%xmm1	\n\t	movaps	0x390(%%rax),%%xmm6	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x180(%%rax)	\n\t	movaps	%%xmm5,0x380(%%rax)	\n\t"\
		"movaps	%%xmm2,0x190(%%rax)	\n\t	movaps	%%xmm7,0x390(%%rax)	\n\t"\
		/* xD.yD: */						/* x1D.y1D: */\
		"movaps	0x1a0(%%rbx),%%xmm2	\n\t	movaps	0x3a0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x1b0(%%rbx),%%xmm3	\n\t	movaps	0x3b0(%%rbx),%%xmm8	\n\t"\
		"subpd	0x1a0(%%rcx),%%xmm2	\n\t	subpd	0x3a0(%%rcx),%%xmm7	\n\t"\
		"subpd	0x1b0(%%rcx),%%xmm3	\n\t	subpd	0x3b0(%%rcx),%%xmm8	\n\t"\
		"movaps	0x1a0(%%rax),%%xmm0	\n\t	movaps	0x3a0(%%rax),%%xmm5	\n\t"\
		"movaps	0x1b0(%%rax),%%xmm1	\n\t	movaps	0x3b0(%%rax),%%xmm6	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x1a0(%%rax)	\n\t	movaps	%%xmm5,0x3a0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x1b0(%%rax)	\n\t	movaps	%%xmm7,0x3b0(%%rax)	\n\t"\
		/* xE.yE: */						/* x1E.y1E: */\
		"movaps	0x1c0(%%rbx),%%xmm2	\n\t	movaps	0x3c0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x1d0(%%rbx),%%xmm3	\n\t	movaps	0x3d0(%%rbx),%%xmm8	\n\t"\
		"subpd	0x1c0(%%rcx),%%xmm2	\n\t	subpd	0x3c0(%%rcx),%%xmm7	\n\t"\
		"subpd	0x1d0(%%rcx),%%xmm3	\n\t	subpd	0x3d0(%%rcx),%%xmm8	\n\t"\
		"movaps	0x1c0(%%rax),%%xmm0	\n\t	movaps	0x3c0(%%rax),%%xmm5	\n\t"\
		"movaps	0x1d0(%%rax),%%xmm1	\n\t	movaps	0x3d0(%%rax),%%xmm6	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x1c0(%%rax)	\n\t	movaps	%%xmm5,0x3c0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x1d0(%%rax)	\n\t	movaps	%%xmm7,0x3d0(%%rax)	\n\t"\
		/* xF.yF: */						/* x1F.y1F: */\
		"movaps	0x1e0(%%rbx),%%xmm2	\n\t	movaps	0x3e0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x1f0(%%rbx),%%xmm3	\n\t	movaps	0x3f0(%%rbx),%%xmm8	\n\t"\
		"subpd	0x1e0(%%rcx),%%xmm2	\n\t	subpd	0x3e0(%%rcx),%%xmm7	\n\t"\
		"subpd	0x1f0(%%rcx),%%xmm3	\n\t	subpd	0x3f0(%%rcx),%%xmm8	\n\t"\
		"movaps	0x1e0(%%rax),%%xmm0	\n\t	movaps	0x3e0(%%rax),%%xmm5	\n\t"\
		"movaps	0x1f0(%%rax),%%xmm1	\n\t	movaps	0x3f0(%%rax),%%xmm6	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x1e0(%%rax)	\n\t	movaps	%%xmm5,0x3e0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x1f0(%%rax)	\n\t	movaps	%%xmm7,0x3f0(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xa)	/* All inputs from memory addresses here */\
		 ,[__bdd0] "m" (Xb)	/* B-array base-address */\
		 ,[__cdd0] "m" (Xc)	/* C-array base-address */\
		: "cc","memory","rax","rbx","rcx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"	/* Clobbered registers */\
	  );\
	}

#endif	// AVX or SSE2?

/*********** Mul a*b: ***********/

#ifdef USE_ARM_V8_SIMD

	// No Fermat-mod support on ARMv8, just supply a stub macro:
	#define SIMD_MUL(Xa,Xb)\
	{\
	  __asm__ volatile (\
		"ldr x0,%[__add0]	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xa)	/* All inputs from memory addresses here */\
		 ,[__bdd0] "m" (Xb)	/* B-array base-address */\
		: "cc","memory","x0"	/* Clobbered registers */\
	  );\
	}

#elif defined(USE_AVX512)	// AVX512 implements a 512-bit-register version of the the AVX2 ALL_FMA-macro

	#define SIMD_MUL(Xa,Xb)\
	{\
	  __asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movq	%[__bdd0],%%rbx		\n\t"\
		/* x0.y0: */							/* x10.y10: */\
		"vmovaps		 (%%rax),%%zmm0	\n\t	vmovaps	0x800(%%rax),%%zmm5	\n\t"/* x.re */\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t	vmovaps	0x840(%%rax),%%zmm6	\n\t"/* x.im */\
		"vmovaps		 (%%rbx),%%zmm2	\n\t	vmovaps	0x800(%%rbx),%%zmm7	\n\t"/* y.re */\
		"vmovaps	0x040(%%rbx),%%zmm3	\n\t	vmovaps	0x840(%%rbx),%%zmm8	\n\t"/* y.im */\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"/* copy x.re */\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"/* x.re *= y.re */\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"/* y.re *= x.im */\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"/* x.im *= y.im */\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"/* y.im *= x.re[copy] */\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"/* z.re = x.re*y.re - x.im*y.im */\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"/* z.im = y.re*x.im + y.im*x.re */\
		"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm5,0x800(%%rax)	\n\t"/* write z.re */\
		"vmovaps	%%zmm2,0x040(%%rax)	\n\t	vmovaps	%%zmm7,0x840(%%rax)	\n\t"/* write z.im */\
		/* x1.y1: */							/* x11.y11: */\
		"vmovaps	0x080(%%rax),%%zmm0	\n\t	vmovaps	0x880(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x0c0(%%rax),%%zmm1	\n\t	vmovaps	0x8c0(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x080(%%rbx),%%zmm2	\n\t	vmovaps	0x880(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x0c0(%%rbx),%%zmm3	\n\t	vmovaps	0x8c0(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x080(%%rax)	\n\t	vmovaps	%%zmm5,0x880(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x0c0(%%rax)	\n\t	vmovaps	%%zmm7,0x8c0(%%rax)	\n\t"\
		/* x2.y2: */							/* x12.y12: */\
		"vmovaps	0x100(%%rax),%%zmm0	\n\t	vmovaps	0x900(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x140(%%rax),%%zmm1	\n\t	vmovaps	0x940(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x100(%%rbx),%%zmm2	\n\t	vmovaps	0x900(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x140(%%rbx),%%zmm3	\n\t	vmovaps	0x940(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x100(%%rax)	\n\t	vmovaps	%%zmm5,0x900(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x140(%%rax)	\n\t	vmovaps	%%zmm7,0x940(%%rax)	\n\t"\
		/* x3.y3: */							/* x13.y13: */\
		"vmovaps	0x180(%%rax),%%zmm0	\n\t	vmovaps	0x980(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x1c0(%%rax),%%zmm1	\n\t	vmovaps	0x9c0(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x180(%%rbx),%%zmm2	\n\t	vmovaps	0x980(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x1c0(%%rbx),%%zmm3	\n\t	vmovaps	0x9c0(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x180(%%rax)	\n\t	vmovaps	%%zmm5,0x980(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x1c0(%%rax)	\n\t	vmovaps	%%zmm7,0x9c0(%%rax)	\n\t"\
		/* x4.y4: */							/* x14.y14: */\
		"vmovaps	0x200(%%rax),%%zmm0	\n\t	vmovaps	0xa00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x240(%%rax),%%zmm1	\n\t	vmovaps	0xa40(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x200(%%rbx),%%zmm2	\n\t	vmovaps	0xa00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x240(%%rbx),%%zmm3	\n\t	vmovaps	0xa40(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x200(%%rax)	\n\t	vmovaps	%%zmm5,0xa00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x240(%%rax)	\n\t	vmovaps	%%zmm7,0xa40(%%rax)	\n\t"\
		/* x5.y5: */							/* x15.y15: */\
		"vmovaps	0x280(%%rax),%%zmm0	\n\t	vmovaps	0xa80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x2c0(%%rax),%%zmm1	\n\t	vmovaps	0xac0(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x280(%%rbx),%%zmm2	\n\t	vmovaps	0xa80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x2c0(%%rbx),%%zmm3	\n\t	vmovaps	0xac0(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x280(%%rax)	\n\t	vmovaps	%%zmm5,0xa80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x2c0(%%rax)	\n\t	vmovaps	%%zmm7,0xac0(%%rax)	\n\t"\
		/* x6.y6: */							/* x16.y16: */\
		"vmovaps	0x300(%%rax),%%zmm0	\n\t	vmovaps	0xb00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x340(%%rax),%%zmm1	\n\t	vmovaps	0xb40(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x300(%%rbx),%%zmm2	\n\t	vmovaps	0xb00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x340(%%rbx),%%zmm3	\n\t	vmovaps	0xb40(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x300(%%rax)	\n\t	vmovaps	%%zmm5,0xb00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x340(%%rax)	\n\t	vmovaps	%%zmm7,0xb40(%%rax)	\n\t"\
		/* x7.y7: */							/* x17.y17: */\
		"vmovaps	0x380(%%rax),%%zmm0	\n\t	vmovaps	0xb80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x3c0(%%rax),%%zmm1	\n\t	vmovaps	0xbc0(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x380(%%rbx),%%zmm2	\n\t	vmovaps	0xb80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x3c0(%%rbx),%%zmm3	\n\t	vmovaps	0xbc0(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x380(%%rax)	\n\t	vmovaps	%%zmm5,0xb80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x3c0(%%rax)	\n\t	vmovaps	%%zmm7,0xbc0(%%rax)	\n\t"\
		/* x8.y8: */							/* x18.y18: */\
		"vmovaps	0x400(%%rax),%%zmm0	\n\t	vmovaps	0xc00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x440(%%rax),%%zmm1	\n\t	vmovaps	0xc40(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x400(%%rbx),%%zmm2	\n\t	vmovaps	0xc00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x440(%%rbx),%%zmm3	\n\t	vmovaps	0xc40(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x400(%%rax)	\n\t	vmovaps	%%zmm5,0xc00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x440(%%rax)	\n\t	vmovaps	%%zmm7,0xc40(%%rax)	\n\t"\
		/* x9.y9: */							/* x19.y19: */\
		"vmovaps	0x480(%%rax),%%zmm0	\n\t	vmovaps	0xc80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x4c0(%%rax),%%zmm1	\n\t	vmovaps	0xcc0(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x480(%%rbx),%%zmm2	\n\t	vmovaps	0xc80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x4c0(%%rbx),%%zmm3	\n\t	vmovaps	0xcc0(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x480(%%rax)	\n\t	vmovaps	%%zmm5,0xc80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x4c0(%%rax)	\n\t	vmovaps	%%zmm7,0xcc0(%%rax)	\n\t"\
		/* xA.yA: */							/* x1A.y1A: */\
		"vmovaps	0x500(%%rax),%%zmm0	\n\t	vmovaps	0xd00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x540(%%rax),%%zmm1	\n\t	vmovaps	0xd40(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x500(%%rbx),%%zmm2	\n\t	vmovaps	0xd00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x540(%%rbx),%%zmm3	\n\t	vmovaps	0xd40(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x500(%%rax)	\n\t	vmovaps	%%zmm5,0xd00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x540(%%rax)	\n\t	vmovaps	%%zmm7,0xd40(%%rax)	\n\t"\
		/* xB.yB: */							/* x1B.y1B: */\
		"vmovaps	0x580(%%rax),%%zmm0	\n\t	vmovaps	0xd80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x5c0(%%rax),%%zmm1	\n\t	vmovaps	0xdc0(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x580(%%rbx),%%zmm2	\n\t	vmovaps	0xd80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x5c0(%%rbx),%%zmm3	\n\t	vmovaps	0xdc0(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x580(%%rax)	\n\t	vmovaps	%%zmm5,0xd80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x5c0(%%rax)	\n\t	vmovaps	%%zmm7,0xdc0(%%rax)	\n\t"\
		/* xC.yC: */							/* x1C.y1C: */\
		"vmovaps	0x600(%%rax),%%zmm0	\n\t	vmovaps	0xe00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x640(%%rax),%%zmm1	\n\t	vmovaps	0xe40(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x600(%%rbx),%%zmm2	\n\t	vmovaps	0xe00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x640(%%rbx),%%zmm3	\n\t	vmovaps	0xe40(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x600(%%rax)	\n\t	vmovaps	%%zmm5,0xe00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x640(%%rax)	\n\t	vmovaps	%%zmm7,0xe40(%%rax)	\n\t"\
		/* xD.yD: */							/* x1D.y1D: */\
		"vmovaps	0x680(%%rax),%%zmm0	\n\t	vmovaps	0xe80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x6c0(%%rax),%%zmm1	\n\t	vmovaps	0xec0(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x680(%%rbx),%%zmm2	\n\t	vmovaps	0xe80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x6c0(%%rbx),%%zmm3	\n\t	vmovaps	0xec0(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x680(%%rax)	\n\t	vmovaps	%%zmm5,0xe80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x6c0(%%rax)	\n\t	vmovaps	%%zmm7,0xec0(%%rax)	\n\t"\
		/* xE.yE: */							/* x1E.y1E: */\
		"vmovaps	0x700(%%rax),%%zmm0	\n\t	vmovaps	0xf00(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x740(%%rax),%%zmm1	\n\t	vmovaps	0xf40(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x700(%%rbx),%%zmm2	\n\t	vmovaps	0xf00(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x740(%%rbx),%%zmm3	\n\t	vmovaps	0xf40(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x700(%%rax)	\n\t	vmovaps	%%zmm5,0xf00(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x740(%%rax)	\n\t	vmovaps	%%zmm7,0xf40(%%rax)	\n\t"\
		/* xF.yF: */							/* x1F.y1F: */\
		"vmovaps	0x780(%%rax),%%zmm0	\n\t	vmovaps	0xf80(%%rax),%%zmm5	\n\t"\
		"vmovaps	0x7c0(%%rax),%%zmm1	\n\t	vmovaps	0xfc0(%%rax),%%zmm6	\n\t"\
		"vmovaps	0x780(%%rbx),%%zmm2	\n\t	vmovaps	0xf80(%%rbx),%%zmm7	\n\t"\
		"vmovaps	0x7c0(%%rbx),%%zmm3	\n\t	vmovaps	0xfc0(%%rbx),%%zmm8	\n\t"\
		"vmovaps		  %%zmm0,%%zmm4	\n\t	vmovaps		  %%zmm5,%%zmm9	\n\t"\
		"vmulpd	%%zmm2,%%zmm0,%%zmm0	\n\t	vmulpd	%%zmm7,%%zmm5,%%zmm5	\n\t"\
		"vmulpd	%%zmm1,%%zmm2,%%zmm2	\n\t	vmulpd	%%zmm6,%%zmm7,%%zmm7	\n\t"\
		"vmulpd	%%zmm3,%%zmm1,%%zmm1	\n\t	vmulpd	%%zmm8,%%zmm6,%%zmm6	\n\t"\
		"vmulpd	%%zmm4,%%zmm3,%%zmm3	\n\t	vmulpd	%%zmm9,%%zmm8,%%zmm8	\n\t"\
		"vsubpd	%%zmm1,%%zmm0,%%zmm0	\n\t	vsubpd	%%zmm6,%%zmm5,%%zmm5	\n\t"\
		"vaddpd	%%zmm3,%%zmm2,%%zmm2	\n\t	vaddpd	%%zmm8,%%zmm7,%%zmm7	\n\t"\
		"vmovaps	%%zmm0,0x780(%%rax)	\n\t	vmovaps	%%zmm5,0xf80(%%rax)	\n\t"\
		"vmovaps	%%zmm2,0x7c0(%%rax)	\n\t	vmovaps	%%zmm7,0xfc0(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xa)	/* All inputs from memory addresses here */\
		 ,[__bdd0] "m" (Xb)	/* B-array base-address */\
		: "cc","memory","rax","rbx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"	/* Clobbered registers */\
	  );\
	}

#elif defined(USE_AVX)

	#define SIMD_MUL(Xa,Xb)\
	{\
	  __asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movq	%[__bdd0],%%rbx		\n\t"\
		/* x0.y0: */							/* x10.y10: */\
		"vmovaps		 (%%rax),%%ymm0	\n\t	vmovaps	0x400(%%rax),%%ymm5	\n\t"/* x.re */\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t	vmovaps	0x420(%%rax),%%ymm6	\n\t"/* x.im */\
		"vmovaps		 (%%rbx),%%ymm2	\n\t	vmovaps	0x400(%%rbx),%%ymm7	\n\t"/* y.re */\
		"vmovaps	0x020(%%rbx),%%ymm3	\n\t	vmovaps	0x420(%%rbx),%%ymm8	\n\t"/* y.im */\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"/* copy x.re */\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"/* x.re *= y.re */\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"/* y.re *= x.im */\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"/* x.im *= y.im */\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"/* y.im *= x.re[copy] */\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"/* z.re = x.re*y.re - x.im*y.im */\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"/* z.im = y.re*x.im + y.im*x.re */\
		"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm5,0x400(%%rax)	\n\t"/* write z.re */\
		"vmovaps	%%ymm2,0x020(%%rax)	\n\t	vmovaps	%%ymm7,0x420(%%rax)	\n\t"/* write z.im */\
		/* x1.y1: */							/* x11.y11: */\
		"vmovaps	0x040(%%rax),%%ymm0	\n\t	vmovaps	0x440(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x060(%%rax),%%ymm1	\n\t	vmovaps	0x460(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x040(%%rbx),%%ymm2	\n\t	vmovaps	0x440(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x060(%%rbx),%%ymm3	\n\t	vmovaps	0x460(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x040(%%rax)	\n\t	vmovaps	%%ymm5,0x440(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x060(%%rax)	\n\t	vmovaps	%%ymm7,0x460(%%rax)	\n\t"\
		/* x2.y2: */							/* x12.y12: */\
		"vmovaps	0x080(%%rax),%%ymm0	\n\t	vmovaps	0x480(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x0a0(%%rax),%%ymm1	\n\t	vmovaps	0x4a0(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x080(%%rbx),%%ymm2	\n\t	vmovaps	0x480(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x0a0(%%rbx),%%ymm3	\n\t	vmovaps	0x4a0(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x080(%%rax)	\n\t	vmovaps	%%ymm5,0x480(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x0a0(%%rax)	\n\t	vmovaps	%%ymm7,0x4a0(%%rax)	\n\t"\
		/* x3.y3: */							/* x13.y13: */\
		"vmovaps	0x0c0(%%rax),%%ymm0	\n\t	vmovaps	0x4c0(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x0e0(%%rax),%%ymm1	\n\t	vmovaps	0x4e0(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x0c0(%%rbx),%%ymm2	\n\t	vmovaps	0x4c0(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x0e0(%%rbx),%%ymm3	\n\t	vmovaps	0x4e0(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x0c0(%%rax)	\n\t	vmovaps	%%ymm5,0x4c0(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x0e0(%%rax)	\n\t	vmovaps	%%ymm7,0x4e0(%%rax)	\n\t"\
		/* x4.y4: */							/* x14.y14: */\
		"vmovaps	0x100(%%rax),%%ymm0	\n\t	vmovaps	0x500(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x120(%%rax),%%ymm1	\n\t	vmovaps	0x520(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x100(%%rbx),%%ymm2	\n\t	vmovaps	0x500(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x120(%%rbx),%%ymm3	\n\t	vmovaps	0x520(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x100(%%rax)	\n\t	vmovaps	%%ymm5,0x500(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x120(%%rax)	\n\t	vmovaps	%%ymm7,0x520(%%rax)	\n\t"\
		/* x5.y5: */							/* x15.y15: */\
		"vmovaps	0x140(%%rax),%%ymm0	\n\t	vmovaps	0x540(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x160(%%rax),%%ymm1	\n\t	vmovaps	0x560(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x140(%%rbx),%%ymm2	\n\t	vmovaps	0x540(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x160(%%rbx),%%ymm3	\n\t	vmovaps	0x560(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x140(%%rax)	\n\t	vmovaps	%%ymm5,0x540(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x160(%%rax)	\n\t	vmovaps	%%ymm7,0x560(%%rax)	\n\t"\
		/* x6.y6: */							/* x16.y16: */\
		"vmovaps	0x180(%%rax),%%ymm0	\n\t	vmovaps	0x580(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x1a0(%%rax),%%ymm1	\n\t	vmovaps	0x5a0(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x180(%%rbx),%%ymm2	\n\t	vmovaps	0x580(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x1a0(%%rbx),%%ymm3	\n\t	vmovaps	0x5a0(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x180(%%rax)	\n\t	vmovaps	%%ymm5,0x580(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x1a0(%%rax)	\n\t	vmovaps	%%ymm7,0x5a0(%%rax)	\n\t"\
		/* x7.y7: */							/* x17.y17: */\
		"vmovaps	0x1c0(%%rax),%%ymm0	\n\t	vmovaps	0x5c0(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x1e0(%%rax),%%ymm1	\n\t	vmovaps	0x5e0(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x1c0(%%rbx),%%ymm2	\n\t	vmovaps	0x5c0(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x1e0(%%rbx),%%ymm3	\n\t	vmovaps	0x5e0(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x1c0(%%rax)	\n\t	vmovaps	%%ymm5,0x5c0(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x1e0(%%rax)	\n\t	vmovaps	%%ymm7,0x5e0(%%rax)	\n\t"\
		/* x8.y8: */							/* x18.y18: */\
		"vmovaps	0x200(%%rax),%%ymm0	\n\t	vmovaps	0x600(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x220(%%rax),%%ymm1	\n\t	vmovaps	0x620(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x200(%%rbx),%%ymm2	\n\t	vmovaps	0x600(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x220(%%rbx),%%ymm3	\n\t	vmovaps	0x620(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x200(%%rax)	\n\t	vmovaps	%%ymm5,0x600(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x220(%%rax)	\n\t	vmovaps	%%ymm7,0x620(%%rax)	\n\t"\
		/* x9.y9: */							/* x19.y19: */\
		"vmovaps	0x240(%%rax),%%ymm0	\n\t	vmovaps	0x640(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x260(%%rax),%%ymm1	\n\t	vmovaps	0x660(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x240(%%rbx),%%ymm2	\n\t	vmovaps	0x640(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x260(%%rbx),%%ymm3	\n\t	vmovaps	0x660(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x240(%%rax)	\n\t	vmovaps	%%ymm5,0x640(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x260(%%rax)	\n\t	vmovaps	%%ymm7,0x660(%%rax)	\n\t"\
		/* xA.yA: */							/* x1A.y1A: */\
		"vmovaps	0x280(%%rax),%%ymm0	\n\t	vmovaps	0x680(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x2a0(%%rax),%%ymm1	\n\t	vmovaps	0x6a0(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x280(%%rbx),%%ymm2	\n\t	vmovaps	0x680(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x2a0(%%rbx),%%ymm3	\n\t	vmovaps	0x6a0(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x280(%%rax)	\n\t	vmovaps	%%ymm5,0x680(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x2a0(%%rax)	\n\t	vmovaps	%%ymm7,0x6a0(%%rax)	\n\t"\
		/* xB.yB: */							/* x1B.y1B: */\
		"vmovaps	0x2c0(%%rax),%%ymm0	\n\t	vmovaps	0x6c0(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x2e0(%%rax),%%ymm1	\n\t	vmovaps	0x6e0(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x2c0(%%rbx),%%ymm2	\n\t	vmovaps	0x6c0(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x2e0(%%rbx),%%ymm3	\n\t	vmovaps	0x6e0(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x2c0(%%rax)	\n\t	vmovaps	%%ymm5,0x6c0(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x2e0(%%rax)	\n\t	vmovaps	%%ymm7,0x6e0(%%rax)	\n\t"\
		/* xC.yC: */							/* x1C.y1C: */\
		"vmovaps	0x300(%%rax),%%ymm0	\n\t	vmovaps	0x700(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x320(%%rax),%%ymm1	\n\t	vmovaps	0x720(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x300(%%rbx),%%ymm2	\n\t	vmovaps	0x700(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x320(%%rbx),%%ymm3	\n\t	vmovaps	0x720(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x300(%%rax)	\n\t	vmovaps	%%ymm5,0x700(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x320(%%rax)	\n\t	vmovaps	%%ymm7,0x720(%%rax)	\n\t"\
		/* xD.yD: */							/* x1D.y1D: */\
		"vmovaps	0x340(%%rax),%%ymm0	\n\t	vmovaps	0x740(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x360(%%rax),%%ymm1	\n\t	vmovaps	0x760(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x340(%%rbx),%%ymm2	\n\t	vmovaps	0x740(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x360(%%rbx),%%ymm3	\n\t	vmovaps	0x760(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x340(%%rax)	\n\t	vmovaps	%%ymm5,0x740(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x360(%%rax)	\n\t	vmovaps	%%ymm7,0x760(%%rax)	\n\t"\
		/* xE.yE: */							/* x1E.y1E: */\
		"vmovaps	0x380(%%rax),%%ymm0	\n\t	vmovaps	0x780(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x3a0(%%rax),%%ymm1	\n\t	vmovaps	0x7a0(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x380(%%rbx),%%ymm2	\n\t	vmovaps	0x780(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x3a0(%%rbx),%%ymm3	\n\t	vmovaps	0x7a0(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x380(%%rax)	\n\t	vmovaps	%%ymm5,0x780(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x3a0(%%rax)	\n\t	vmovaps	%%ymm7,0x7a0(%%rax)	\n\t"\
		/* xF.yF: */							/* x1F.y1F: */\
		"vmovaps	0x3c0(%%rax),%%ymm0	\n\t	vmovaps	0x7c0(%%rax),%%ymm5	\n\t"\
		"vmovaps	0x3e0(%%rax),%%ymm1	\n\t	vmovaps	0x7e0(%%rax),%%ymm6	\n\t"\
		"vmovaps	0x3c0(%%rbx),%%ymm2	\n\t	vmovaps	0x7c0(%%rbx),%%ymm7	\n\t"\
		"vmovaps	0x3e0(%%rbx),%%ymm3	\n\t	vmovaps	0x7e0(%%rbx),%%ymm8	\n\t"\
		"vmovaps		  %%ymm0,%%ymm4	\n\t	vmovaps		  %%ymm5,%%ymm9	\n\t"\
		"vmulpd	%%ymm2,%%ymm0,%%ymm0	\n\t	vmulpd	%%ymm7,%%ymm5,%%ymm5	\n\t"\
		"vmulpd	%%ymm1,%%ymm2,%%ymm2	\n\t	vmulpd	%%ymm6,%%ymm7,%%ymm7	\n\t"\
		"vmulpd	%%ymm3,%%ymm1,%%ymm1	\n\t	vmulpd	%%ymm8,%%ymm6,%%ymm6	\n\t"\
		"vmulpd	%%ymm4,%%ymm3,%%ymm3	\n\t	vmulpd	%%ymm9,%%ymm8,%%ymm8	\n\t"\
		"vsubpd	%%ymm1,%%ymm0,%%ymm0	\n\t	vsubpd	%%ymm6,%%ymm5,%%ymm5	\n\t"\
		"vaddpd	%%ymm3,%%ymm2,%%ymm2	\n\t	vaddpd	%%ymm8,%%ymm7,%%ymm7	\n\t"\
		"vmovaps	%%ymm0,0x3c0(%%rax)	\n\t	vmovaps	%%ymm5,0x7c0(%%rax)	\n\t"\
		"vmovaps	%%ymm2,0x3e0(%%rax)	\n\t	vmovaps	%%ymm7,0x7e0(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xa)	/* All inputs from memory addresses here */\
		 ,[__bdd0] "m" (Xb)	/* B-array base-address */\
		: "cc","memory","rax","rbx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"	/* Clobbered registers */\
	  );\
	}

#else	// 64-bit SSE2:

	// Inputs are X = [x.re,x.im], Y = [y.re,y.im]
	// Output is Z = [x.re*y.re - x.im*y.im, x.re*y.im + x.im*y.re], overwriting X
	#define SIMD_MUL(Xa,Xb)\
	{\
	  __asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		"movq	%[__bdd0],%%rbx		\n\t"\
		/* x0.y0: */						/* x10.y10: */\
		"movaps		 (%%rax),%%xmm0	\n\t	movaps	0x200(%%rax),%%xmm5	\n\t"/* x.re */\
		"movaps	0x010(%%rax),%%xmm1	\n\t	movaps	0x210(%%rax),%%xmm6	\n\t"/* x.im */\
		"movaps		 (%%rbx),%%xmm2	\n\t	movaps	0x200(%%rbx),%%xmm7	\n\t"/* y.re */\
		"movaps	0x010(%%rbx),%%xmm3	\n\t	movaps	0x210(%%rbx),%%xmm8	\n\t"/* y.im */\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"/* copy x.re */\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"/* x.re *= y.re */\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"/* y.re *= x.im */\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"/* x.im *= y.im */\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"/* y.im *= x.re[copy] */\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"/* z.re = x.re*y.re - x.im*y.im */\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"/* z.im = y.re*x.im + y.im*x.re */\
		"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm5,0x200(%%rax)	\n\t"/* write z.re */\
		"movaps	%%xmm2,0x010(%%rax)	\n\t	movaps	%%xmm7,0x210(%%rax)	\n\t"/* write z.im */\
		/* x1.y1: */						/* x11.y11: */\
		"movaps	0x020(%%rax),%%xmm0	\n\t	movaps	0x220(%%rax),%%xmm5	\n\t"\
		"movaps	0x030(%%rax),%%xmm1	\n\t	movaps	0x230(%%rax),%%xmm6	\n\t"\
		"movaps	0x020(%%rbx),%%xmm2	\n\t	movaps	0x220(%%rbx),%%xmm7	\n\t"\
		"movaps	0x030(%%rbx),%%xmm3	\n\t	movaps	0x230(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x020(%%rax)	\n\t	movaps	%%xmm5,0x220(%%rax)	\n\t"\
		"movaps	%%xmm2,0x030(%%rax)	\n\t	movaps	%%xmm7,0x230(%%rax)	\n\t"\
		/* x2.y2: */						/* x12.y12: */\
		"movaps	0x040(%%rax),%%xmm0	\n\t	movaps	0x240(%%rax),%%xmm5	\n\t"\
		"movaps	0x050(%%rax),%%xmm1	\n\t	movaps	0x250(%%rax),%%xmm6	\n\t"\
		"movaps	0x040(%%rbx),%%xmm2	\n\t	movaps	0x240(%%rbx),%%xmm7	\n\t"\
		"movaps	0x050(%%rbx),%%xmm3	\n\t	movaps	0x250(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x040(%%rax)	\n\t	movaps	%%xmm5,0x240(%%rax)	\n\t"\
		"movaps	%%xmm2,0x050(%%rax)	\n\t	movaps	%%xmm7,0x250(%%rax)	\n\t"\
		/* x3.y3: */						/* x13.y13: */\
		"movaps	0x060(%%rax),%%xmm0	\n\t	movaps	0x260(%%rax),%%xmm5	\n\t"\
		"movaps	0x070(%%rax),%%xmm1	\n\t	movaps	0x270(%%rax),%%xmm6	\n\t"\
		"movaps	0x060(%%rbx),%%xmm2	\n\t	movaps	0x260(%%rbx),%%xmm7	\n\t"\
		"movaps	0x070(%%rbx),%%xmm3	\n\t	movaps	0x270(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x060(%%rax)	\n\t	movaps	%%xmm5,0x260(%%rax)	\n\t"\
		"movaps	%%xmm2,0x070(%%rax)	\n\t	movaps	%%xmm7,0x270(%%rax)	\n\t"\
		/* x4.y4: */						/* x14.y14: */\
		"movaps	0x080(%%rax),%%xmm0	\n\t	movaps	0x280(%%rax),%%xmm5	\n\t"\
		"movaps	0x090(%%rax),%%xmm1	\n\t	movaps	0x290(%%rax),%%xmm6	\n\t"\
		"movaps	0x080(%%rbx),%%xmm2	\n\t	movaps	0x280(%%rbx),%%xmm7	\n\t"\
		"movaps	0x090(%%rbx),%%xmm3	\n\t	movaps	0x290(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x080(%%rax)	\n\t	movaps	%%xmm5,0x280(%%rax)	\n\t"\
		"movaps	%%xmm2,0x090(%%rax)	\n\t	movaps	%%xmm7,0x290(%%rax)	\n\t"\
		/* x5.y5: */						/* x15.y15: */\
		"movaps	0x0a0(%%rax),%%xmm0	\n\t	movaps	0x2a0(%%rax),%%xmm5	\n\t"\
		"movaps	0x0b0(%%rax),%%xmm1	\n\t	movaps	0x2b0(%%rax),%%xmm6	\n\t"\
		"movaps	0x0a0(%%rbx),%%xmm2	\n\t	movaps	0x2a0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x0b0(%%rbx),%%xmm3	\n\t	movaps	0x2b0(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x0a0(%%rax)	\n\t	movaps	%%xmm5,0x2a0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x0b0(%%rax)	\n\t	movaps	%%xmm7,0x2b0(%%rax)	\n\t"\
		/* x6.y6: */						/* x16.y16: */\
		"movaps	0x0c0(%%rax),%%xmm0	\n\t	movaps	0x2c0(%%rax),%%xmm5	\n\t"\
		"movaps	0x0d0(%%rax),%%xmm1	\n\t	movaps	0x2d0(%%rax),%%xmm6	\n\t"\
		"movaps	0x0c0(%%rbx),%%xmm2	\n\t	movaps	0x2c0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x0d0(%%rbx),%%xmm3	\n\t	movaps	0x2d0(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x0c0(%%rax)	\n\t	movaps	%%xmm5,0x2c0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x0d0(%%rax)	\n\t	movaps	%%xmm7,0x2d0(%%rax)	\n\t"\
		/* x7.y7: */						/* x17.y17: */\
		"movaps	0x0e0(%%rax),%%xmm0	\n\t	movaps	0x2e0(%%rax),%%xmm5	\n\t"\
		"movaps	0x0f0(%%rax),%%xmm1	\n\t	movaps	0x2f0(%%rax),%%xmm6	\n\t"\
		"movaps	0x0e0(%%rbx),%%xmm2	\n\t	movaps	0x2e0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x0f0(%%rbx),%%xmm3	\n\t	movaps	0x2f0(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x0e0(%%rax)	\n\t	movaps	%%xmm5,0x2e0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x0f0(%%rax)	\n\t	movaps	%%xmm7,0x2f0(%%rax)	\n\t"\
		/* x8.y8: */						/* x18.y18: */\
		"movaps	0x100(%%rax),%%xmm0	\n\t	movaps	0x300(%%rax),%%xmm5	\n\t"\
		"movaps	0x110(%%rax),%%xmm1	\n\t	movaps	0x310(%%rax),%%xmm6	\n\t"\
		"movaps	0x100(%%rbx),%%xmm2	\n\t	movaps	0x300(%%rbx),%%xmm7	\n\t"\
		"movaps	0x110(%%rbx),%%xmm3	\n\t	movaps	0x310(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x100(%%rax)	\n\t	movaps	%%xmm5,0x300(%%rax)	\n\t"\
		"movaps	%%xmm2,0x110(%%rax)	\n\t	movaps	%%xmm7,0x310(%%rax)	\n\t"\
		/* x9.y9: */						/* x19.y19: */\
		"movaps	0x120(%%rax),%%xmm0	\n\t	movaps	0x320(%%rax),%%xmm5	\n\t"\
		"movaps	0x130(%%rax),%%xmm1	\n\t	movaps	0x330(%%rax),%%xmm6	\n\t"\
		"movaps	0x120(%%rbx),%%xmm2	\n\t	movaps	0x320(%%rbx),%%xmm7	\n\t"\
		"movaps	0x130(%%rbx),%%xmm3	\n\t	movaps	0x330(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x120(%%rax)	\n\t	movaps	%%xmm5,0x320(%%rax)	\n\t"\
		"movaps	%%xmm2,0x130(%%rax)	\n\t	movaps	%%xmm7,0x330(%%rax)	\n\t"\
		/* xA.yA: */						/* x1A.y1A: */\
		"movaps	0x140(%%rax),%%xmm0	\n\t	movaps	0x340(%%rax),%%xmm5	\n\t"\
		"movaps	0x150(%%rax),%%xmm1	\n\t	movaps	0x350(%%rax),%%xmm6	\n\t"\
		"movaps	0x140(%%rbx),%%xmm2	\n\t	movaps	0x340(%%rbx),%%xmm7	\n\t"\
		"movaps	0x150(%%rbx),%%xmm3	\n\t	movaps	0x350(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x140(%%rax)	\n\t	movaps	%%xmm5,0x340(%%rax)	\n\t"\
		"movaps	%%xmm2,0x150(%%rax)	\n\t	movaps	%%xmm7,0x350(%%rax)	\n\t"\
		/* xB.yB: */						/* x1B.y1B: */\
		"movaps	0x160(%%rax),%%xmm0	\n\t	movaps	0x360(%%rax),%%xmm5	\n\t"\
		"movaps	0x170(%%rax),%%xmm1	\n\t	movaps	0x370(%%rax),%%xmm6	\n\t"\
		"movaps	0x160(%%rbx),%%xmm2	\n\t	movaps	0x360(%%rbx),%%xmm7	\n\t"\
		"movaps	0x170(%%rbx),%%xmm3	\n\t	movaps	0x370(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x160(%%rax)	\n\t	movaps	%%xmm5,0x360(%%rax)	\n\t"\
		"movaps	%%xmm2,0x170(%%rax)	\n\t	movaps	%%xmm7,0x370(%%rax)	\n\t"\
		/* xC.yC: */						/* x1C.y1C: */\
		"movaps	0x180(%%rax),%%xmm0	\n\t	movaps	0x380(%%rax),%%xmm5	\n\t"\
		"movaps	0x190(%%rax),%%xmm1	\n\t	movaps	0x390(%%rax),%%xmm6	\n\t"\
		"movaps	0x180(%%rbx),%%xmm2	\n\t	movaps	0x380(%%rbx),%%xmm7	\n\t"\
		"movaps	0x190(%%rbx),%%xmm3	\n\t	movaps	0x390(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x180(%%rax)	\n\t	movaps	%%xmm5,0x380(%%rax)	\n\t"\
		"movaps	%%xmm2,0x190(%%rax)	\n\t	movaps	%%xmm7,0x390(%%rax)	\n\t"\
		/* xD.yD: */						/* x1D.y1D: */\
		"movaps	0x1a0(%%rax),%%xmm0	\n\t	movaps	0x3a0(%%rax),%%xmm5	\n\t"\
		"movaps	0x1b0(%%rax),%%xmm1	\n\t	movaps	0x3b0(%%rax),%%xmm6	\n\t"\
		"movaps	0x1a0(%%rbx),%%xmm2	\n\t	movaps	0x3a0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x1b0(%%rbx),%%xmm3	\n\t	movaps	0x3b0(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x1a0(%%rax)	\n\t	movaps	%%xmm5,0x3a0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x1b0(%%rax)	\n\t	movaps	%%xmm7,0x3b0(%%rax)	\n\t"\
		/* xE.yE: */						/* x1E.y1E: */\
		"movaps	0x1c0(%%rax),%%xmm0	\n\t	movaps	0x3c0(%%rax),%%xmm5	\n\t"\
		"movaps	0x1d0(%%rax),%%xmm1	\n\t	movaps	0x3d0(%%rax),%%xmm6	\n\t"\
		"movaps	0x1c0(%%rbx),%%xmm2	\n\t	movaps	0x3c0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x1d0(%%rbx),%%xmm3	\n\t	movaps	0x3d0(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x1c0(%%rax)	\n\t	movaps	%%xmm5,0x3c0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x1d0(%%rax)	\n\t	movaps	%%xmm7,0x3d0(%%rax)	\n\t"\
		/* xF.yF: */						/* x1F.y1F: */\
		"movaps	0x1e0(%%rax),%%xmm0	\n\t	movaps	0x3e0(%%rax),%%xmm5	\n\t"\
		"movaps	0x1f0(%%rax),%%xmm1	\n\t	movaps	0x3f0(%%rax),%%xmm6	\n\t"\
		"movaps	0x1e0(%%rbx),%%xmm2	\n\t	movaps	0x3e0(%%rbx),%%xmm7	\n\t"\
		"movaps	0x1f0(%%rbx),%%xmm3	\n\t	movaps	0x3f0(%%rbx),%%xmm8	\n\t"\
		"movaps		  %%xmm0,%%xmm4	\n\t	movaps		  %%xmm5,%%xmm9	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm7,%%xmm5	\n\t"\
		"mulpd		  %%xmm1,%%xmm2	\n\t	mulpd		  %%xmm6,%%xmm7	\n\t"\
		"mulpd		  %%xmm3,%%xmm1	\n\t	mulpd		  %%xmm8,%%xmm6	\n\t"\
		"mulpd		  %%xmm4,%%xmm3	\n\t	mulpd		  %%xmm9,%%xmm8	\n\t"\
		"subpd		  %%xmm1,%%xmm0	\n\t	subpd		  %%xmm6,%%xmm5	\n\t"\
		"addpd		  %%xmm3,%%xmm2	\n\t	addpd		  %%xmm8,%%xmm7	\n\t"\
		"movaps	%%xmm0,0x1e0(%%rax)	\n\t	movaps	%%xmm5,0x3e0(%%rax)	\n\t"\
		"movaps	%%xmm2,0x1f0(%%rax)	\n\t	movaps	%%xmm7,0x3f0(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xa)	/* All inputs from memory addresses here */\
		 ,[__bdd0] "m" (Xb)	/* B-array base-address */\
		: "cc","memory","rax","rbx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7","xmm8","xmm9"	/* Clobbered registers */\
	  );\
	}

#endif	// AVX or SSE2?

/*********** Square a^2: ***********/

#ifdef USE_ARM_V8_SIMD

	// No Fermat-mod support on ARMv8, just supply a stub macro:
	#define SIMD_SQR(Xa)\
	{\
	__asm__ volatile (\
		"ldr x0,%[__add0]	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xa)	/* All inputs from memory addresses here */\
		: "cc","memory","x0"	/* Clobbered registers */\
	  );\
	}

#elif defined(USE_AVX512)	// AVX512 implements a 512-bit-register version of the the AVX2 ALL_FMA-macro

	#define SIMD_SQR(Xa)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax			\n\t"\
		/* z0^2: */								/*z10^2: */\
		"vmovaps		 (%%rax),%%zmm0	\n\t	vmovaps	0x800(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x040(%%rax),%%zmm1	\n\t	vmovaps	0x840(%%rax),%%zmm4	\n\t"\
		"vmovaps		  %%zmm0,%%zmm2	\n\t	vmovaps		  %%zmm3,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm0,%%zmm0		\n\t	vaddpd		  %%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd		  %%zmm1,%%zmm2,%%zmm2		\n\t	vsubpd		  %%zmm4,%%zmm5,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd		  %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmulpd		  %%zmm2,%%zmm0,%%zmm0		\n\t	vmulpd		  %%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmulpd		 (%%rax),%%zmm1,%%zmm1		\n\t	vmulpd	0x800(%%rax),%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm0,     (%%rax)	\n\t	vmovaps	%%zmm3,0x800(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x040(%%rax)	\n\t	vmovaps	%%zmm4,0x840(%%rax)	\n\t"\
		/* z1^2: */								/*z11^2: */\
		"vmovaps	0x080(%%rax),%%zmm0	\n\t	vmovaps	0x880(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x0c0(%%rax),%%zmm1	\n\t	vmovaps	0x8c0(%%rax),%%zmm4	\n\t"\
		"vmovaps		  %%zmm0,%%zmm2	\n\t	vmovaps		  %%zmm3,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm0,%%zmm0		\n\t	vaddpd		  %%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd		  %%zmm1,%%zmm2,%%zmm2		\n\t	vsubpd		  %%zmm4,%%zmm5,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd		  %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmulpd		  %%zmm2,%%zmm0,%%zmm0		\n\t	vmulpd		  %%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x080(%%rax),%%zmm1,%%zmm1		\n\t	vmulpd	0x880(%%rax),%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm0,0x080(%%rax)	\n\t	vmovaps	%%zmm3,0x880(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x0c0(%%rax)	\n\t	vmovaps	%%zmm4,0x8c0(%%rax)	\n\t"\
		/* z2^2: */								/*z12^2: */\
		"vmovaps	0x100(%%rax),%%zmm0	\n\t	vmovaps	0x900(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x140(%%rax),%%zmm1	\n\t	vmovaps	0x940(%%rax),%%zmm4	\n\t"\
		"vmovaps		  %%zmm0,%%zmm2	\n\t	vmovaps		  %%zmm3,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm0,%%zmm0		\n\t	vaddpd		  %%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd		  %%zmm1,%%zmm2,%%zmm2		\n\t	vsubpd		  %%zmm4,%%zmm5,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd		  %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmulpd		  %%zmm2,%%zmm0,%%zmm0		\n\t	vmulpd		  %%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x100(%%rax),%%zmm1,%%zmm1		\n\t	vmulpd	0x900(%%rax),%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm0,0x100(%%rax)	\n\t	vmovaps	%%zmm3,0x900(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x140(%%rax)	\n\t	vmovaps	%%zmm4,0x940(%%rax)	\n\t"\
		/* z3^2: */								/*z13^2: */\
		"vmovaps	0x180(%%rax),%%zmm0	\n\t	vmovaps	0x980(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x1c0(%%rax),%%zmm1	\n\t	vmovaps	0x9c0(%%rax),%%zmm4	\n\t"\
		"vmovaps		  %%zmm0,%%zmm2	\n\t	vmovaps		  %%zmm3,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm0,%%zmm0		\n\t	vaddpd		  %%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd		  %%zmm1,%%zmm2,%%zmm2		\n\t	vsubpd		  %%zmm4,%%zmm5,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd		  %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmulpd		  %%zmm2,%%zmm0,%%zmm0		\n\t	vmulpd		  %%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x180(%%rax),%%zmm1,%%zmm1		\n\t	vmulpd	0x980(%%rax),%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm0,0x180(%%rax)	\n\t	vmovaps	%%zmm3,0x980(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x1c0(%%rax)	\n\t	vmovaps	%%zmm4,0x9c0(%%rax)	\n\t"\
		/* z4^2: */								/*z14^2: */\
		"vmovaps	0x200(%%rax),%%zmm0	\n\t	vmovaps	0xa00(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x240(%%rax),%%zmm1	\n\t	vmovaps	0xa40(%%rax),%%zmm4	\n\t"\
		"vmovaps		  %%zmm0,%%zmm2	\n\t	vmovaps		  %%zmm3,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm0,%%zmm0		\n\t	vaddpd		  %%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd		  %%zmm1,%%zmm2,%%zmm2		\n\t	vsubpd		  %%zmm4,%%zmm5,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd		  %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmulpd		  %%zmm2,%%zmm0,%%zmm0		\n\t	vmulpd		  %%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x200(%%rax),%%zmm1,%%zmm1		\n\t	vmulpd	0xa00(%%rax),%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm0,0x200(%%rax)	\n\t	vmovaps	%%zmm3,0xa00(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x240(%%rax)	\n\t	vmovaps	%%zmm4,0xa40(%%rax)	\n\t"\
		/* z5^2: */								/*z15^2: */\
		"vmovaps	0x280(%%rax),%%zmm0	\n\t	vmovaps	0xa80(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x2c0(%%rax),%%zmm1	\n\t	vmovaps	0xac0(%%rax),%%zmm4	\n\t"\
		"vmovaps		  %%zmm0,%%zmm2	\n\t	vmovaps		  %%zmm3,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm0,%%zmm0		\n\t	vaddpd		  %%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd		  %%zmm1,%%zmm2,%%zmm2		\n\t	vsubpd		  %%zmm4,%%zmm5,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd		  %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmulpd		  %%zmm2,%%zmm0,%%zmm0		\n\t	vmulpd		  %%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x280(%%rax),%%zmm1,%%zmm1		\n\t	vmulpd	0xa80(%%rax),%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm0,0x280(%%rax)	\n\t	vmovaps	%%zmm3,0xa80(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x2c0(%%rax)	\n\t	vmovaps	%%zmm4,0xac0(%%rax)	\n\t"\
		/* z6^2: */								/*z16^2: */\
		"vmovaps	0x300(%%rax),%%zmm0	\n\t	vmovaps	0xb00(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x340(%%rax),%%zmm1	\n\t	vmovaps	0xb40(%%rax),%%zmm4	\n\t"\
		"vmovaps		  %%zmm0,%%zmm2	\n\t	vmovaps		  %%zmm3,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm0,%%zmm0		\n\t	vaddpd		  %%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd		  %%zmm1,%%zmm2,%%zmm2		\n\t	vsubpd		  %%zmm4,%%zmm5,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd		  %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmulpd		  %%zmm2,%%zmm0,%%zmm0		\n\t	vmulpd		  %%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x300(%%rax),%%zmm1,%%zmm1		\n\t	vmulpd	0xb00(%%rax),%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm0,0x300(%%rax)	\n\t	vmovaps	%%zmm3,0xb00(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x340(%%rax)	\n\t	vmovaps	%%zmm4,0xb40(%%rax)	\n\t"\
		/* z7^2: */								/*z17^2: */\
		"vmovaps	0x380(%%rax),%%zmm0	\n\t	vmovaps	0xb80(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x3c0(%%rax),%%zmm1	\n\t	vmovaps	0xbc0(%%rax),%%zmm4	\n\t"\
		"vmovaps		  %%zmm0,%%zmm2	\n\t	vmovaps		  %%zmm3,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm0,%%zmm0		\n\t	vaddpd		  %%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd		  %%zmm1,%%zmm2,%%zmm2		\n\t	vsubpd		  %%zmm4,%%zmm5,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd		  %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmulpd		  %%zmm2,%%zmm0,%%zmm0		\n\t	vmulpd		  %%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x380(%%rax),%%zmm1,%%zmm1		\n\t	vmulpd	0xb80(%%rax),%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm0,0x380(%%rax)	\n\t	vmovaps	%%zmm3,0xb80(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x3c0(%%rax)	\n\t	vmovaps	%%zmm4,0xbc0(%%rax)	\n\t"\
		/* z8^2: */								/*z18^2: */\
		"vmovaps	0x400(%%rax),%%zmm0	\n\t	vmovaps	0xc00(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x440(%%rax),%%zmm1	\n\t	vmovaps	0xc40(%%rax),%%zmm4	\n\t"\
		"vmovaps		  %%zmm0,%%zmm2	\n\t	vmovaps		  %%zmm3,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm0,%%zmm0		\n\t	vaddpd		  %%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd		  %%zmm1,%%zmm2,%%zmm2		\n\t	vsubpd		  %%zmm4,%%zmm5,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd		  %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmulpd		  %%zmm2,%%zmm0,%%zmm0		\n\t	vmulpd		  %%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x400(%%rax),%%zmm1,%%zmm1		\n\t	vmulpd	0xc00(%%rax),%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm0,0x400(%%rax)	\n\t	vmovaps	%%zmm3,0xc00(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x440(%%rax)	\n\t	vmovaps	%%zmm4,0xc40(%%rax)	\n\t"\
		/* z9^2: */								/*z19^2: */\
		"vmovaps	0x480(%%rax),%%zmm0	\n\t	vmovaps	0xc80(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x4c0(%%rax),%%zmm1	\n\t	vmovaps	0xcc0(%%rax),%%zmm4	\n\t"\
		"vmovaps		  %%zmm0,%%zmm2	\n\t	vmovaps		  %%zmm3,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm0,%%zmm0		\n\t	vaddpd		  %%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd		  %%zmm1,%%zmm2,%%zmm2		\n\t	vsubpd		  %%zmm4,%%zmm5,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd		  %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmulpd		  %%zmm2,%%zmm0,%%zmm0		\n\t	vmulpd		  %%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x480(%%rax),%%zmm1,%%zmm1		\n\t	vmulpd	0xc80(%%rax),%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm0,0x480(%%rax)	\n\t	vmovaps	%%zmm3,0xc80(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x4c0(%%rax)	\n\t	vmovaps	%%zmm4,0xcc0(%%rax)	\n\t"\
		/* zA^2: */								/*z1A^2: */\
		"vmovaps	0x500(%%rax),%%zmm0	\n\t	vmovaps	0xd00(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x540(%%rax),%%zmm1	\n\t	vmovaps	0xd40(%%rax),%%zmm4	\n\t"\
		"vmovaps		  %%zmm0,%%zmm2	\n\t	vmovaps		  %%zmm3,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm0,%%zmm0		\n\t	vaddpd		  %%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd		  %%zmm1,%%zmm2,%%zmm2		\n\t	vsubpd		  %%zmm4,%%zmm5,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd		  %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmulpd		  %%zmm2,%%zmm0,%%zmm0		\n\t	vmulpd		  %%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x500(%%rax),%%zmm1,%%zmm1		\n\t	vmulpd	0xd00(%%rax),%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm0,0x500(%%rax)	\n\t	vmovaps	%%zmm3,0xd00(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x540(%%rax)	\n\t	vmovaps	%%zmm4,0xd40(%%rax)	\n\t"\
		/* zB^2: */								/*z1B^2: */\
		"vmovaps	0x580(%%rax),%%zmm0	\n\t	vmovaps	0xd80(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x5c0(%%rax),%%zmm1	\n\t	vmovaps	0xdc0(%%rax),%%zmm4	\n\t"\
		"vmovaps		  %%zmm0,%%zmm2	\n\t	vmovaps		  %%zmm3,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm0,%%zmm0		\n\t	vaddpd		  %%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd		  %%zmm1,%%zmm2,%%zmm2		\n\t	vsubpd		  %%zmm4,%%zmm5,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd		  %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmulpd		  %%zmm2,%%zmm0,%%zmm0		\n\t	vmulpd		  %%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x580(%%rax),%%zmm1,%%zmm1		\n\t	vmulpd	0xd80(%%rax),%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm0,0x580(%%rax)	\n\t	vmovaps	%%zmm3,0xd80(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x5c0(%%rax)	\n\t	vmovaps	%%zmm4,0xdc0(%%rax)	\n\t"\
		/* zC^2: */								/*z1C^2: */\
		"vmovaps	0x600(%%rax),%%zmm0	\n\t	vmovaps	0xe00(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x640(%%rax),%%zmm1	\n\t	vmovaps	0xe40(%%rax),%%zmm4	\n\t"\
		"vmovaps		  %%zmm0,%%zmm2	\n\t	vmovaps		  %%zmm3,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm0,%%zmm0		\n\t	vaddpd		  %%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd		  %%zmm1,%%zmm2,%%zmm2		\n\t	vsubpd		  %%zmm4,%%zmm5,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd		  %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmulpd		  %%zmm2,%%zmm0,%%zmm0		\n\t	vmulpd		  %%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x600(%%rax),%%zmm1,%%zmm1		\n\t	vmulpd	0xe00(%%rax),%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm0,0x600(%%rax)	\n\t	vmovaps	%%zmm3,0xe00(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x640(%%rax)	\n\t	vmovaps	%%zmm4,0xe40(%%rax)	\n\t"\
		/* zD^2: */								/*z1D^2: */\
		"vmovaps	0x680(%%rax),%%zmm0	\n\t	vmovaps	0xe80(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x6c0(%%rax),%%zmm1	\n\t	vmovaps	0xec0(%%rax),%%zmm4	\n\t"\
		"vmovaps		  %%zmm0,%%zmm2	\n\t	vmovaps		  %%zmm3,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm0,%%zmm0		\n\t	vaddpd		  %%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd		  %%zmm1,%%zmm2,%%zmm2		\n\t	vsubpd		  %%zmm4,%%zmm5,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd		  %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmulpd		  %%zmm2,%%zmm0,%%zmm0		\n\t	vmulpd		  %%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x680(%%rax),%%zmm1,%%zmm1		\n\t	vmulpd	0xe80(%%rax),%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm0,0x680(%%rax)	\n\t	vmovaps	%%zmm3,0xe80(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x6c0(%%rax)	\n\t	vmovaps	%%zmm4,0xec0(%%rax)	\n\t"\
		/* zE^2: */								/*z1E^2: */\
		"vmovaps	0x700(%%rax),%%zmm0	\n\t	vmovaps	0xf00(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x740(%%rax),%%zmm1	\n\t	vmovaps	0xf40(%%rax),%%zmm4	\n\t"\
		"vmovaps		  %%zmm0,%%zmm2	\n\t	vmovaps		  %%zmm3,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm0,%%zmm0		\n\t	vaddpd		  %%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd		  %%zmm1,%%zmm2,%%zmm2		\n\t	vsubpd		  %%zmm4,%%zmm5,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd		  %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmulpd		  %%zmm2,%%zmm0,%%zmm0		\n\t	vmulpd		  %%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x700(%%rax),%%zmm1,%%zmm1		\n\t	vmulpd	0xf00(%%rax),%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm0,0x700(%%rax)	\n\t	vmovaps	%%zmm3,0xf00(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x740(%%rax)	\n\t	vmovaps	%%zmm4,0xf40(%%rax)	\n\t"\
		/* zF^2: */								/*z1F^2: */\
		"vmovaps	0x780(%%rax),%%zmm0	\n\t	vmovaps	0xf80(%%rax),%%zmm3	\n\t"\
		"vmovaps	0x7c0(%%rax),%%zmm1	\n\t	vmovaps	0xfc0(%%rax),%%zmm4	\n\t"\
		"vmovaps		  %%zmm0,%%zmm2	\n\t	vmovaps		  %%zmm3,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm0,%%zmm0		\n\t	vaddpd		  %%zmm4,%%zmm3,%%zmm3	\n\t"\
		"vsubpd		  %%zmm1,%%zmm2,%%zmm2		\n\t	vsubpd		  %%zmm4,%%zmm5,%%zmm5	\n\t"\
		"vaddpd		  %%zmm1,%%zmm1,%%zmm1		\n\t	vaddpd		  %%zmm4,%%zmm4,%%zmm4	\n\t"\
		"vmulpd		  %%zmm2,%%zmm0,%%zmm0		\n\t	vmulpd		  %%zmm5,%%zmm3,%%zmm3	\n\t"\
		"vmulpd	0x780(%%rax),%%zmm1,%%zmm1		\n\t	vmulpd	0xf80(%%rax),%%zmm4,%%zmm4		\n\t"\
		"vmovaps	%%zmm0,0x780(%%rax)	\n\t	vmovaps	%%zmm3,0xf80(%%rax)	\n\t"\
		"vmovaps	%%zmm1,0x7c0(%%rax)	\n\t	vmovaps	%%zmm4,0xfc0(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xa)	/* All inputs from memory addresses here */\
		: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5"	/* Clobbered registers */\
	  );\
	}

#elif defined(USE_AVX)

	#define SIMD_SQR(Xa)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax			\n\t"\
		/* z0^2: */								/*z10^2: */\
		"vmovaps		 (%%rax),%%ymm0	\n\t	vmovaps	0x400(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x020(%%rax),%%ymm1	\n\t	vmovaps	0x420(%%rax),%%ymm4	\n\t"\
		"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmulpd		 (%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x400(%%rax),%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,     (%%rax)	\n\t	vmovaps	%%ymm3,0x400(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x020(%%rax)	\n\t	vmovaps	%%ymm4,0x420(%%rax)	\n\t"\
		/* z1^2: */								/*z11^2: */\
		"vmovaps	0x040(%%rax),%%ymm0	\n\t	vmovaps	0x440(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x060(%%rax),%%ymm1	\n\t	vmovaps	0x460(%%rax),%%ymm4	\n\t"\
		"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x040(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x440(%%rax),%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,0x040(%%rax)	\n\t	vmovaps	%%ymm3,0x440(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x060(%%rax)	\n\t	vmovaps	%%ymm4,0x460(%%rax)	\n\t"\
		/* z2^2: */								/*z12^2: */\
		"vmovaps	0x080(%%rax),%%ymm0	\n\t	vmovaps	0x480(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x0a0(%%rax),%%ymm1	\n\t	vmovaps	0x4a0(%%rax),%%ymm4	\n\t"\
		"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x080(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x480(%%rax),%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,0x080(%%rax)	\n\t	vmovaps	%%ymm3,0x480(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x0a0(%%rax)	\n\t	vmovaps	%%ymm4,0x4a0(%%rax)	\n\t"\
		/* z3^2: */								/*z13^2: */\
		"vmovaps	0x0c0(%%rax),%%ymm0	\n\t	vmovaps	0x4c0(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x0e0(%%rax),%%ymm1	\n\t	vmovaps	0x4e0(%%rax),%%ymm4	\n\t"\
		"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x0c0(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x4c0(%%rax),%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,0x0c0(%%rax)	\n\t	vmovaps	%%ymm3,0x4c0(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x0e0(%%rax)	\n\t	vmovaps	%%ymm4,0x4e0(%%rax)	\n\t"\
		/* z4^2: */								/*z14^2: */\
		"vmovaps	0x100(%%rax),%%ymm0	\n\t	vmovaps	0x500(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x120(%%rax),%%ymm1	\n\t	vmovaps	0x520(%%rax),%%ymm4	\n\t"\
		"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x100(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x500(%%rax),%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,0x100(%%rax)	\n\t	vmovaps	%%ymm3,0x500(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x120(%%rax)	\n\t	vmovaps	%%ymm4,0x520(%%rax)	\n\t"\
		/* z5^2: */								/*z15^2: */\
		"vmovaps	0x140(%%rax),%%ymm0	\n\t	vmovaps	0x540(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x160(%%rax),%%ymm1	\n\t	vmovaps	0x560(%%rax),%%ymm4	\n\t"\
		"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x140(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x540(%%rax),%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,0x140(%%rax)	\n\t	vmovaps	%%ymm3,0x540(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x160(%%rax)	\n\t	vmovaps	%%ymm4,0x560(%%rax)	\n\t"\
		/* z6^2: */								/*z16^2: */\
		"vmovaps	0x180(%%rax),%%ymm0	\n\t	vmovaps	0x580(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x1a0(%%rax),%%ymm1	\n\t	vmovaps	0x5a0(%%rax),%%ymm4	\n\t"\
		"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x180(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x580(%%rax),%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,0x180(%%rax)	\n\t	vmovaps	%%ymm3,0x580(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x1a0(%%rax)	\n\t	vmovaps	%%ymm4,0x5a0(%%rax)	\n\t"\
		/* z7^2: */								/*z17^2: */\
		"vmovaps	0x1c0(%%rax),%%ymm0	\n\t	vmovaps	0x5c0(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x1e0(%%rax),%%ymm1	\n\t	vmovaps	0x5e0(%%rax),%%ymm4	\n\t"\
		"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x1c0(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x5c0(%%rax),%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,0x1c0(%%rax)	\n\t	vmovaps	%%ymm3,0x5c0(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x1e0(%%rax)	\n\t	vmovaps	%%ymm4,0x5e0(%%rax)	\n\t"\
		/* z8^2: */								/*z18^2: */\
		"vmovaps	0x200(%%rax),%%ymm0	\n\t	vmovaps	0x600(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x220(%%rax),%%ymm1	\n\t	vmovaps	0x620(%%rax),%%ymm4	\n\t"\
		"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x200(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x600(%%rax),%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,0x200(%%rax)	\n\t	vmovaps	%%ymm3,0x600(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x220(%%rax)	\n\t	vmovaps	%%ymm4,0x620(%%rax)	\n\t"\
		/* z9^2: */								/*z19^2: */\
		"vmovaps	0x240(%%rax),%%ymm0	\n\t	vmovaps	0x640(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x260(%%rax),%%ymm1	\n\t	vmovaps	0x660(%%rax),%%ymm4	\n\t"\
		"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x240(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x640(%%rax),%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,0x240(%%rax)	\n\t	vmovaps	%%ymm3,0x640(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x260(%%rax)	\n\t	vmovaps	%%ymm4,0x660(%%rax)	\n\t"\
		/* zA^2: */								/*z1A^2: */\
		"vmovaps	0x280(%%rax),%%ymm0	\n\t	vmovaps	0x680(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x2a0(%%rax),%%ymm1	\n\t	vmovaps	0x6a0(%%rax),%%ymm4	\n\t"\
		"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x280(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x680(%%rax),%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,0x280(%%rax)	\n\t	vmovaps	%%ymm3,0x680(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x2a0(%%rax)	\n\t	vmovaps	%%ymm4,0x6a0(%%rax)	\n\t"\
		/* zB^2: */								/*z1B^2: */\
		"vmovaps	0x2c0(%%rax),%%ymm0	\n\t	vmovaps	0x6c0(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x2e0(%%rax),%%ymm1	\n\t	vmovaps	0x6e0(%%rax),%%ymm4	\n\t"\
		"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x2c0(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x6c0(%%rax),%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,0x2c0(%%rax)	\n\t	vmovaps	%%ymm3,0x6c0(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x2e0(%%rax)	\n\t	vmovaps	%%ymm4,0x6e0(%%rax)	\n\t"\
		/* zC^2: */								/*z1C^2: */\
		"vmovaps	0x300(%%rax),%%ymm0	\n\t	vmovaps	0x700(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x320(%%rax),%%ymm1	\n\t	vmovaps	0x720(%%rax),%%ymm4	\n\t"\
		"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x300(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x700(%%rax),%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,0x300(%%rax)	\n\t	vmovaps	%%ymm3,0x700(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x320(%%rax)	\n\t	vmovaps	%%ymm4,0x720(%%rax)	\n\t"\
		/* zD^2: */								/*z1D^2: */\
		"vmovaps	0x340(%%rax),%%ymm0	\n\t	vmovaps	0x740(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x360(%%rax),%%ymm1	\n\t	vmovaps	0x760(%%rax),%%ymm4	\n\t"\
		"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x340(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x740(%%rax),%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,0x340(%%rax)	\n\t	vmovaps	%%ymm3,0x740(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x360(%%rax)	\n\t	vmovaps	%%ymm4,0x760(%%rax)	\n\t"\
		/* zE^2: */								/*z1E^2: */\
		"vmovaps	0x380(%%rax),%%ymm0	\n\t	vmovaps	0x780(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x3a0(%%rax),%%ymm1	\n\t	vmovaps	0x7a0(%%rax),%%ymm4	\n\t"\
		"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x380(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x780(%%rax),%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,0x380(%%rax)	\n\t	vmovaps	%%ymm3,0x780(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x3a0(%%rax)	\n\t	vmovaps	%%ymm4,0x7a0(%%rax)	\n\t"\
		/* zF^2: */								/*z1F^2: */\
		"vmovaps	0x3c0(%%rax),%%ymm0	\n\t	vmovaps	0x7c0(%%rax),%%ymm3	\n\t"\
		"vmovaps	0x3e0(%%rax),%%ymm1	\n\t	vmovaps	0x7e0(%%rax),%%ymm4	\n\t"\
		"vmovaps		  %%ymm0,%%ymm2	\n\t	vmovaps		  %%ymm3,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm0,%%ymm0		\n\t	vaddpd		  %%ymm4,%%ymm3,%%ymm3	\n\t"\
		"vsubpd		  %%ymm1,%%ymm2,%%ymm2		\n\t	vsubpd		  %%ymm4,%%ymm5,%%ymm5	\n\t"\
		"vaddpd		  %%ymm1,%%ymm1,%%ymm1		\n\t	vaddpd		  %%ymm4,%%ymm4,%%ymm4	\n\t"\
		"vmulpd		  %%ymm2,%%ymm0,%%ymm0		\n\t	vmulpd		  %%ymm5,%%ymm3,%%ymm3	\n\t"\
		"vmulpd	0x3c0(%%rax),%%ymm1,%%ymm1		\n\t	vmulpd	0x7c0(%%rax),%%ymm4,%%ymm4		\n\t"\
		"vmovaps	%%ymm0,0x3c0(%%rax)	\n\t	vmovaps	%%ymm3,0x7c0(%%rax)	\n\t"\
		"vmovaps	%%ymm1,0x3e0(%%rax)	\n\t	vmovaps	%%ymm4,0x7e0(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xa)	/* All inputs from memory addresses here */\
		: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5"	/* Clobbered registers */\
	  );\
	}

#else	// 64-bit SSE2:

	#define SIMD_SQR(Xa)\
	{\
	__asm__ volatile (\
		"movq	%[__add0],%%rax		\n\t"\
		/* z0^2: */							/*z10^2: */\
		"movaps		 (%%rax),%%xmm0	\n\t	movaps	0x200(%%rax),%%xmm3	\n\t"\
		"movaps	0x010(%%rax),%%xmm1	\n\t	movaps	0x210(%%rax),%%xmm4	\n\t"\
		"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
		"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
		"mulpd		 (%%rax),%%xmm1	\n\t	mulpd	0x200(%%rax),%%xmm4	\n\t"\
		"movaps	%%xmm0,     (%%rax)	\n\t	movaps	%%xmm3,0x200(%%rax)	\n\t"\
		"movaps	%%xmm1,0x010(%%rax)	\n\t	movaps	%%xmm4,0x210(%%rax)	\n\t"\
		/* z1^2: */							/*z11^2: */\
		"movaps	0x020(%%rax),%%xmm0	\n\t	movaps	0x220(%%rax),%%xmm3	\n\t"\
		"movaps	0x030(%%rax),%%xmm1	\n\t	movaps	0x230(%%rax),%%xmm4	\n\t"\
		"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
		"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
		"mulpd	0x020(%%rax),%%xmm1	\n\t	mulpd	0x220(%%rax),%%xmm4	\n\t"\
		"movaps	%%xmm0,0x020(%%rax)	\n\t	movaps	%%xmm3,0x220(%%rax)	\n\t"\
		"movaps	%%xmm1,0x030(%%rax)	\n\t	movaps	%%xmm4,0x230(%%rax)	\n\t"\
		/* z2^2: */							/*z12^2: */\
		"movaps	0x040(%%rax),%%xmm0	\n\t	movaps	0x240(%%rax),%%xmm3	\n\t"\
		"movaps	0x050(%%rax),%%xmm1	\n\t	movaps	0x250(%%rax),%%xmm4	\n\t"\
		"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
		"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
		"mulpd	0x040(%%rax),%%xmm1	\n\t	mulpd	0x240(%%rax),%%xmm4	\n\t"\
		"movaps	%%xmm0,0x040(%%rax)	\n\t	movaps	%%xmm3,0x240(%%rax)	\n\t"\
		"movaps	%%xmm1,0x050(%%rax)	\n\t	movaps	%%xmm4,0x250(%%rax)	\n\t"\
		/* z3^2: */							/*z13^2: */\
		"movaps	0x060(%%rax),%%xmm0	\n\t	movaps	0x260(%%rax),%%xmm3	\n\t"\
		"movaps	0x070(%%rax),%%xmm1	\n\t	movaps	0x270(%%rax),%%xmm4	\n\t"\
		"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
		"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
		"mulpd	0x060(%%rax),%%xmm1	\n\t	mulpd	0x260(%%rax),%%xmm4	\n\t"\
		"movaps	%%xmm0,0x060(%%rax)	\n\t	movaps	%%xmm3,0x260(%%rax)	\n\t"\
		"movaps	%%xmm1,0x070(%%rax)	\n\t	movaps	%%xmm4,0x270(%%rax)	\n\t"\
		/* z4^2: */							/*z14^2: */\
		"movaps	0x080(%%rax),%%xmm0	\n\t	movaps	0x280(%%rax),%%xmm3	\n\t"\
		"movaps	0x090(%%rax),%%xmm1	\n\t	movaps	0x290(%%rax),%%xmm4	\n\t"\
		"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
		"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
		"mulpd	0x080(%%rax),%%xmm1	\n\t	mulpd	0x280(%%rax),%%xmm4	\n\t"\
		"movaps	%%xmm0,0x080(%%rax)	\n\t	movaps	%%xmm3,0x280(%%rax)	\n\t"\
		"movaps	%%xmm1,0x090(%%rax)	\n\t	movaps	%%xmm4,0x290(%%rax)	\n\t"\
		/* z5^2: */							/*z15^2: */\
		"movaps	0x0a0(%%rax),%%xmm0	\n\t	movaps	0x2a0(%%rax),%%xmm3	\n\t"\
		"movaps	0x0b0(%%rax),%%xmm1	\n\t	movaps	0x2b0(%%rax),%%xmm4	\n\t"\
		"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
		"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
		"mulpd	0x0a0(%%rax),%%xmm1	\n\t	mulpd	0x2a0(%%rax),%%xmm4	\n\t"\
		"movaps	%%xmm0,0x0a0(%%rax)	\n\t	movaps	%%xmm3,0x2a0(%%rax)	\n\t"\
		"movaps	%%xmm1,0x0b0(%%rax)	\n\t	movaps	%%xmm4,0x2b0(%%rax)	\n\t"\
		/* z6^2: */							/*z16^2: */\
		"movaps	0x0c0(%%rax),%%xmm0	\n\t	movaps	0x2c0(%%rax),%%xmm3	\n\t"\
		"movaps	0x0d0(%%rax),%%xmm1	\n\t	movaps	0x2d0(%%rax),%%xmm4	\n\t"\
		"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
		"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
		"mulpd	0x0c0(%%rax),%%xmm1	\n\t	mulpd	0x2c0(%%rax),%%xmm4	\n\t"\
		"movaps	%%xmm0,0x0c0(%%rax)	\n\t	movaps	%%xmm3,0x2c0(%%rax)	\n\t"\
		"movaps	%%xmm1,0x0d0(%%rax)	\n\t	movaps	%%xmm4,0x2d0(%%rax)	\n\t"\
		/* z7^2: */							/*z17^2: */\
		"movaps	0x0e0(%%rax),%%xmm0	\n\t	movaps	0x2e0(%%rax),%%xmm3	\n\t"\
		"movaps	0x0f0(%%rax),%%xmm1	\n\t	movaps	0x2f0(%%rax),%%xmm4	\n\t"\
		"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
		"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
		"mulpd	0x0e0(%%rax),%%xmm1	\n\t	mulpd	0x2e0(%%rax),%%xmm4	\n\t"\
		"movaps	%%xmm0,0x0e0(%%rax)	\n\t	movaps	%%xmm3,0x2e0(%%rax)	\n\t"\
		"movaps	%%xmm1,0x0f0(%%rax)	\n\t	movaps	%%xmm4,0x2f0(%%rax)	\n\t"\
		/* z8^2: */							/*z18^2: */\
		"movaps	0x100(%%rax),%%xmm0	\n\t	movaps	0x300(%%rax),%%xmm3	\n\t"\
		"movaps	0x110(%%rax),%%xmm1	\n\t	movaps	0x310(%%rax),%%xmm4	\n\t"\
		"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
		"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
		"mulpd	0x100(%%rax),%%xmm1	\n\t	mulpd	0x300(%%rax),%%xmm4	\n\t"\
		"movaps	%%xmm0,0x100(%%rax)	\n\t	movaps	%%xmm3,0x300(%%rax)	\n\t"\
		"movaps	%%xmm1,0x110(%%rax)	\n\t	movaps	%%xmm4,0x310(%%rax)	\n\t"\
		/* z9^2: */							/*z19^2: */\
		"movaps	0x120(%%rax),%%xmm0	\n\t	movaps	0x320(%%rax),%%xmm3	\n\t"\
		"movaps	0x130(%%rax),%%xmm1	\n\t	movaps	0x330(%%rax),%%xmm4	\n\t"\
		"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
		"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
		"mulpd	0x120(%%rax),%%xmm1	\n\t	mulpd	0x320(%%rax),%%xmm4	\n\t"\
		"movaps	%%xmm0,0x120(%%rax)	\n\t	movaps	%%xmm3,0x320(%%rax)	\n\t"\
		"movaps	%%xmm1,0x130(%%rax)	\n\t	movaps	%%xmm4,0x330(%%rax)	\n\t"\
		/* zA^2: */							/*z1A^2: */\
		"movaps	0x140(%%rax),%%xmm0	\n\t	movaps	0x340(%%rax),%%xmm3	\n\t"\
		"movaps	0x150(%%rax),%%xmm1	\n\t	movaps	0x350(%%rax),%%xmm4	\n\t"\
		"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
		"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
		"mulpd	0x140(%%rax),%%xmm1	\n\t	mulpd	0x340(%%rax),%%xmm4	\n\t"\
		"movaps	%%xmm0,0x140(%%rax)	\n\t	movaps	%%xmm3,0x340(%%rax)	\n\t"\
		"movaps	%%xmm1,0x150(%%rax)	\n\t	movaps	%%xmm4,0x350(%%rax)	\n\t"\
		/* zB^2: */							/*z1B^2: */\
		"movaps	0x160(%%rax),%%xmm0	\n\t	movaps	0x360(%%rax),%%xmm3	\n\t"\
		"movaps	0x170(%%rax),%%xmm1	\n\t	movaps	0x370(%%rax),%%xmm4	\n\t"\
		"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
		"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
		"mulpd	0x160(%%rax),%%xmm1	\n\t	mulpd	0x360(%%rax),%%xmm4	\n\t"\
		"movaps	%%xmm0,0x160(%%rax)	\n\t	movaps	%%xmm3,0x360(%%rax)	\n\t"\
		"movaps	%%xmm1,0x170(%%rax)	\n\t	movaps	%%xmm4,0x370(%%rax)	\n\t"\
		/* zC^2: */							/*z1C^2: */\
		"movaps	0x180(%%rax),%%xmm0	\n\t	movaps	0x380(%%rax),%%xmm3	\n\t"\
		"movaps	0x190(%%rax),%%xmm1	\n\t	movaps	0x390(%%rax),%%xmm4	\n\t"\
		"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
		"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
		"mulpd	0x180(%%rax),%%xmm1	\n\t	mulpd	0x380(%%rax),%%xmm4	\n\t"\
		"movaps	%%xmm0,0x180(%%rax)	\n\t	movaps	%%xmm3,0x380(%%rax)	\n\t"\
		"movaps	%%xmm1,0x190(%%rax)	\n\t	movaps	%%xmm4,0x390(%%rax)	\n\t"\
		/* zD^2: */							/*z1D^2: */\
		"movaps	0x1a0(%%rax),%%xmm0	\n\t	movaps	0x3a0(%%rax),%%xmm3	\n\t"\
		"movaps	0x1b0(%%rax),%%xmm1	\n\t	movaps	0x3b0(%%rax),%%xmm4	\n\t"\
		"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
		"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
		"mulpd	0x1a0(%%rax),%%xmm1	\n\t	mulpd	0x3a0(%%rax),%%xmm4	\n\t"\
		"movaps	%%xmm0,0x1a0(%%rax)	\n\t	movaps	%%xmm3,0x3a0(%%rax)	\n\t"\
		"movaps	%%xmm1,0x1b0(%%rax)	\n\t	movaps	%%xmm4,0x3b0(%%rax)	\n\t"\
		/* zE^2: */							/*z1E^2: */\
		"movaps	0x1c0(%%rax),%%xmm0	\n\t	movaps	0x3c0(%%rax),%%xmm3	\n\t"\
		"movaps	0x1d0(%%rax),%%xmm1	\n\t	movaps	0x3d0(%%rax),%%xmm4	\n\t"\
		"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
		"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
		"mulpd	0x1c0(%%rax),%%xmm1	\n\t	mulpd	0x3c0(%%rax),%%xmm4	\n\t"\
		"movaps	%%xmm0,0x1c0(%%rax)	\n\t	movaps	%%xmm3,0x3c0(%%rax)	\n\t"\
		"movaps	%%xmm1,0x1d0(%%rax)	\n\t	movaps	%%xmm4,0x3d0(%%rax)	\n\t"\
		/* zF^2: */							/*z1F^2: */\
		"movaps	0x1e0(%%rax),%%xmm0	\n\t	movaps	0x3e0(%%rax),%%xmm3	\n\t"\
		"movaps	0x1f0(%%rax),%%xmm1	\n\t	movaps	0x3f0(%%rax),%%xmm4	\n\t"\
		"movaps		  %%xmm0,%%xmm2	\n\t	movaps		  %%xmm3,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm0	\n\t	addpd		  %%xmm4,%%xmm3	\n\t"\
		"subpd		  %%xmm1,%%xmm2	\n\t	subpd		  %%xmm4,%%xmm5	\n\t"\
		"addpd		  %%xmm1,%%xmm1	\n\t	addpd		  %%xmm4,%%xmm4	\n\t"\
		"mulpd		  %%xmm2,%%xmm0	\n\t	mulpd		  %%xmm5,%%xmm3	\n\t"\
		"mulpd	0x1e0(%%rax),%%xmm1	\n\t	mulpd	0x3e0(%%rax),%%xmm4	\n\t"\
		"movaps	%%xmm0,0x1e0(%%rax)	\n\t	movaps	%%xmm3,0x3e0(%%rax)	\n\t"\
		"movaps	%%xmm1,0x1f0(%%rax)	\n\t	movaps	%%xmm4,0x3f0(%%rax)	\n\t"\
		:					/* outputs: none */\
		: [__add0] "m" (Xa)	/* All inputs from memory addresses here */\
		: "cc","memory","rax","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5"	/* Clobbered registers */\
	  );\
	}

#endif	// AVX or SSE2?
