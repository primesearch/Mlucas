/*******************************************************************************
*                                                                              *
*   (C) 1997-2015 by Ernst W. Mayer.                                           *
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
#ifndef carry_gcc_h_included
#define carry_gcc_h_included

	/************** See the Visual-studio-style 32-bit analogs of these in carry.h for commented versions: **********/

	/*************************************************************/
	/**************** FERMAT  -MOD CARRY MACROS ******************/
	/*************************************************************/

#if 1	// Experimental low-reg version with just 1-push/pop and 2 compiler-visible GPR clobbers (max allowed by the fermat-mod
		// carry code in compact-obj-code loop-wrapped form, and even 2-reg is too mch for GCC, so this is clang--buildable-only):

	#define SSE2_fermat_carry_norm_pow2_errcheck(Xdata,Xcy,Xnrt_bits,Xnrtm1,Xidx_offset,Xidx_incr,Xhalf_arr,Xsign_mask,Xadd1,Xadd2, Xadd0)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"movl		%[__add0],%%ecx		\n\t"/* base address for prefetch-from-main-data-array */\
	"prefetcht0	(%%ecx)	\n\t"\
		"movl		%[__nrt_bits],%%ecx	\n\t"\
		"movl	%[__idx_offset],%%eax	\n\t"\
		"shrl		$1,%%eax			\n\t"\
		"movl		%%eax,%%ebx			\n\t"\
		"andl		%[__nrtm1],%%eax	\n\t"\
		"shrl		%%cl,%%ebx			\n\t"\
		"shll		$4,%%eax			\n\t"\
		"shll		$4,%%ebx			\n\t"\
		"addl		%[__add1],%%eax		\n\t"\
		"addl		%[__add2],%%ebx		\n\t"\
		"movaps		(%%eax),%%xmm0		\n\t"\
		"movaps		(%%ebx),%%xmm1		\n\t"\
		"movl	%[__idx_offset],%%eax	\n\t"\
		"movaps		%%xmm1,%%xmm2		\n\t"\
		"shufpd	$1,	%%xmm2,%%xmm2		\n\t"\
		"mulpd		%%xmm0,%%xmm1		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"addl		$2,%%eax			\n\t"\
		"shrl		$1,%%eax			\n\t"\
		"movl		%%eax,%%ebx			\n\t"\
		"andl		%[__nrtm1],%%eax	\n\t"\
		"shrl		%%cl,%%ebx			\n\t"\
		"shll		$4,%%eax			\n\t"\
		"shll		$4,%%ebx			\n\t"\
		"addl		%[__add1],%%eax		\n\t"\
		"addl		%[__add2],%%ebx		\n\t"\
		"movaps		(%%eax),%%xmm0		\n\t"\
		"movaps		(%%ebx),%%xmm3		\n\t"\
		"movl	%[__idx_offset],%%eax	\n\t"\
		"movaps		%%xmm3,%%xmm4		\n\t"\
		"shufpd	$1,	%%xmm4,%%xmm4		\n\t"\
		"mulpd		%%xmm0,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"movaps		%%xmm1,%%xmm0		\n\t"\
		"unpcklpd	%%xmm3,%%xmm0		\n\t"\
		"unpckhpd	%%xmm3,%%xmm1		\n\t"\
		"subpd		%%xmm1,%%xmm0		\n\t"\
		"movaps		%%xmm2,%%xmm1		\n\t"\
		"unpcklpd	%%xmm4,%%xmm1		\n\t"\
		"unpckhpd	%%xmm4,%%xmm2		\n\t"\
		"addpd		%%xmm2,%%xmm1		\n\t"\
		"movl		%[__half_arr],%%ecx	\n\t"\
		"movl		%[__data],%%edx		\n\t"\
		"movaps		     (%%edx),%%xmm4	\n\t"\
		"movaps		 0x10(%%edx),%%xmm2	\n\t"\
		"movaps		0x020(%%ecx),%%xmm5	\n\t"\
		"mulpd		%%xmm5,%%xmm4		\n\t"\
		"mulpd		%%xmm5,%%xmm2		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t"\
		"addpd		%%xmm3,%%xmm4		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t"\
		"movl		%[__cy],%%ebx		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"shufpd	$0,	%%xmm2,%%xmm4		\n\t"\
		"shufpd	$3,	%%xmm2,%%xmm5		\n\t"\
		"addpd		     (%%ebx),%%xmm4	\n\t"\
		"movaps		-0x20(%%ecx),%%xmm6	\n\t"\
		"movaps		-0x10(%%ecx),%%xmm7	\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"addpd		%%xmm7,%%xmm4		\n\t"\
		"subpd		%%xmm7,%%xmm4		\n\t"\
		"movl		%[__sign_mask],%%eax\n\t"\
		"subpd		%%xmm4,%%xmm2		\n\t"\
		"andpd		     (%%eax),%%xmm2	\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"mulpd		0x10(%%ecx),%%xmm2	\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t"\
		"subpd		%%xmm7,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		    (%%ecx),%%xmm3	\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t"\
		"addpd		%%xmm7,%%xmm5		\n\t"\
		"subpd		%%xmm7,%%xmm5		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t"\
		"andpd		     (%%eax),%%xmm2	\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t"\
		"mulpd		 0x10(%%ecx),%%xmm2	\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t"\
		"subpd		%%xmm7,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		     (%%ecx),%%xmm3	\n\t"\
		"subpd		%%xmm3,%%xmm5		\n\t"\
		"movaps		%%xmm2,(%%ebx)		\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"shufpd	$0,	%%xmm5,%%xmm4		\n\t"\
		"shufpd	$3,	%%xmm5,%%xmm2		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"movaps		%%xmm6,-0x20(%%ecx)	\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t"\
		"movaps		%%xmm4,    (%%edx)	\n\t"\
		"movaps		%%xmm5,0x10(%%edx)	\n\t"\
		"movl	%[__idx_offset],%%ebx	\n\t"\
		"addl	%[__idx_incr],%%ebx		\n\t"\
		"movl	%%ebx, %[__idx_offset]	\n\t"/* Store incremented idx_offset */\
	"popl %%ebx	\n\t"\
	:						/* outputs: none */\
	:	[__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
	,	[__cy]			"m" (Xcy)\
	,	[__nrt_bits]	"m" (Xnrt_bits)\
	,	[__nrtm1]		"m" (Xnrtm1)\
	,	[__idx_offset]	"m" (Xidx_offset)\
	,	[__idx_incr]	"m" (Xidx_incr)\
	,	[__half_arr]	"m" (Xhalf_arr)\
	,	[__sign_mask]	"m" (Xsign_mask)\
	,	[__add1]		"m" (Xadd1)\
	,	[__add2]		"m" (Xadd2)\
	/* Prefetch address */\
	,	[__add0] "m" (Xadd0)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

#else

	#define SSE2_fermat_carry_norm_pow2_errcheck(Xdata,Xcy,Xnrt_bits,Xnrtm1,Xidx_offset,Xidx_incr,Xhalf_arr,Xsign_mask,Xadd1,Xadd2, Xadd0)\
	{\
	__asm__ volatile (\
	"movl		%[__add0],%%ecx		\n\t"/* base address for prefetch-from-main-data-array */\
	"prefetcht0	(%%ecx)	\n\t"\
		"movl		%[__nrt_bits],%%ecx	\n\t"\
		"movl	%[__idx_offset],%%eax	\n\t"\
		"shrl		$1,%%eax			\n\t"\
		"movl		%%eax,%%esi			\n\t"\
		"andl		%[__nrtm1],%%eax	\n\t"\
		"shrl		%%cl,%%esi			\n\t"\
		"shll		$4,%%eax			\n\t"\
		"shll		$4,%%esi			\n\t"\
		"addl		%[__add1],%%eax		\n\t"\
		"addl		%[__add2],%%esi		\n\t"\
		"movaps		(%%eax),%%xmm0		\n\t"\
		"movaps		(%%esi),%%xmm1		\n\t"\
		"movl	%[__idx_offset],%%eax	\n\t"\
		"movaps		%%xmm1,%%xmm2		\n\t"\
		"shufpd	$1,	%%xmm2,%%xmm2		\n\t"\
		"mulpd		%%xmm0,%%xmm1		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"addl		$2,%%eax			\n\t"\
		"shrl		$1,%%eax			\n\t"\
		"movl		%%eax,%%esi			\n\t"\
		"andl		%[__nrtm1],%%eax	\n\t"\
		"shrl		%%cl,%%esi			\n\t"\
		"shll		$4,%%eax			\n\t"\
		"shll		$4,%%esi			\n\t"\
		"addl		%[__add1],%%eax		\n\t"\
		"addl		%[__add2],%%esi		\n\t"\
		"movaps		(%%eax),%%xmm0		\n\t"\
		"movaps		(%%esi),%%xmm3		\n\t"\
		"movl	%[__idx_offset],%%eax	\n\t"\
		"movaps		%%xmm3,%%xmm4		\n\t"\
		"shufpd	$1,	%%xmm4,%%xmm4		\n\t"\
		"mulpd		%%xmm0,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"movaps		%%xmm1,%%xmm0		\n\t"\
		"unpcklpd	%%xmm3,%%xmm0		\n\t"\
		"unpckhpd	%%xmm3,%%xmm1		\n\t"\
		"subpd		%%xmm1,%%xmm0		\n\t"\
		"movaps		%%xmm2,%%xmm1		\n\t"\
		"unpcklpd	%%xmm4,%%xmm1		\n\t"\
		"unpckhpd	%%xmm4,%%xmm2		\n\t"\
		"addpd		%%xmm2,%%xmm1		\n\t"\
		"movl		%[__half_arr],%%ecx	\n\t"\
		"movl		%[__data],%%edx		\n\t"\
		"movaps		     (%%edx),%%xmm4	\n\t"\
		"movaps		 0x10(%%edx),%%xmm2	\n\t"\
		"movaps		0x020(%%ecx),%%xmm5	\n\t"\
		"mulpd		%%xmm5,%%xmm4		\n\t"\
		"mulpd		%%xmm5,%%xmm2		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t"\
		"addpd		%%xmm3,%%xmm4		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t"\
		"movl		%[__cy],%%esi		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"shufpd	$0,	%%xmm2,%%xmm4		\n\t"\
		"shufpd	$3,	%%xmm2,%%xmm5		\n\t"\
		"addpd		     (%%esi),%%xmm4	\n\t"\
		"movaps		-0x20(%%ecx),%%xmm6	\n\t"\
		"movaps		-0x10(%%ecx),%%xmm7	\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"addpd		%%xmm7,%%xmm4		\n\t"\
		"subpd		%%xmm7,%%xmm4		\n\t"\
		"movl		%[__sign_mask],%%eax\n\t"\
		"subpd		%%xmm4,%%xmm2		\n\t"\
		"andpd		     (%%eax),%%xmm2	\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"mulpd		0x10(%%ecx),%%xmm2	\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t"\
		"subpd		%%xmm7,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		    (%%ecx),%%xmm3	\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t"\
		"addpd		%%xmm7,%%xmm5		\n\t"\
		"subpd		%%xmm7,%%xmm5		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t"\
		"andpd		     (%%eax),%%xmm2	\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t"\
		"mulpd		 0x10(%%ecx),%%xmm2	\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t"\
		"subpd		%%xmm7,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		     (%%ecx),%%xmm3	\n\t"\
		"subpd		%%xmm3,%%xmm5		\n\t"\
		"movaps		%%xmm2,(%%esi)		\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"shufpd	$0,	%%xmm5,%%xmm4		\n\t"\
		"shufpd	$3,	%%xmm5,%%xmm2		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"movaps		%%xmm6,-0x20(%%ecx)	\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t"\
		"movaps		%%xmm4,    (%%edx)	\n\t"\
		"movaps		%%xmm5,0x10(%%edx)	\n\t"\
		"movl	%[__idx_offset],%%esi	\n\t"\
		"addl	%[__idx_incr],%%esi		\n\t"\
		"movl	%%esi, %[__idx_offset]	\n\t"/* Store incremented idx_offset */\
	:						/* outputs: none */\
	:	[__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
	,	[__cy]			"m" (Xcy)\
	,	[__nrt_bits]	"m" (Xnrt_bits)\
	,	[__nrtm1]		"m" (Xnrtm1)\
	,	[__idx_offset]	"m" (Xidx_offset)\
	,	[__idx_incr]	"m" (Xidx_incr)\
	,	[__half_arr]	"m" (Xhalf_arr)\
	,	[__sign_mask]	"m" (Xsign_mask)\
	,	[__add1]		"m" (Xadd1)\
	,	[__add2]		"m" (Xadd2)\
	/* Prefetch address */\
	,	[__add0] "m" (Xadd0)\
		: "cc","memory","eax","esi","ecx","edx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

#endif

#if 1	// Experimental low-reg version with just 1-push/pop and 2 compiler-visible GPR clobbers (max allowed by the fermat-mod
		// carry code in compact-obj-code loop-wrapped form, and even 2-reg is too mch for GCC, so this is clang--buildable-only):

	/* 01 Sep 2014: Unable to cut GPR-use to < 3 here, and only ebx allows push/pop ... this worked fie in single-testcase mode
	but gave "bus error" in radix-240 call of the macro when doing all-radices-test of -f 24 -fftlen 960 ... debug-digging
	(unable to sue gcc due to reg-constraints, had to use clang -Os -g3 -ggdb) showed 0-addr error, but all macro-args looked
	OK. Attempting to cut one more reg by replacing the %%cl shift-count in the 2 shrl using it with $%[__nrt_bits] gave this error:

		../radix240_main_carry_loop.h:689:4: error: unknown token in expression
							SSE2_fermat...
							^
	In file included from ../radix240_ditN_cy_dif1.c:23:
	In file included from ../Mlucas.h:30:
	In file included from ../carry.h:5060:
	../carry_gcc32.h:170:2: note: instantiated from:
			"pushl %%ebx    \n\t"## Explicit save/restore of PIC register ##\
			^
	<inline asm>:22:9: note: instantiated into assembly here
			shrl    $(%ebx),%ebx                    
					  ^

	i.e. problem is that compiler is storing __nrt_bits in ebx on input, collides with above syntax. (I believe this lack of
	user cntrol over what-compiler-uses-ebx-for is responsible for the bus error that started this debug-digging, as well.)

	After trying several unsuccessful workaround, eventually reverted to original "bus error" version of macro, tried the
	same self-tests once more ... and it worked! Difference proved to be the addition of the -g3 -ggdb flags to the -Os opt-level
	in the compile args. *** Need to see if this "accdental workaround" is more generally usefu;. ***
	*/
	#define SSE2_fermat_carry_norm_errcheck(Xdata,Xcy,Xnrt_bits,Xnrtm1,Xidx_offset,Xidx_incr,Xodd_radix,Xhalf_arr,Xsign_mask,Xadd1,Xadd2,Xoffset0,Xoffset1, Xadd0)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"movl		%[__add0],%%ecx		\n\t"/* base address for prefetch-from-main-data-array */\
	"prefetcht0	(%%ecx)	\n\t"\
		"movl	%[__nrt_bits],%%ecx		\n\t"/* __odd_radix was in edx, no more in 3-reg version using just e[abc]x */\
		"movl	%[__idx_offset],%%eax	\n\t"\
		"shrl		$1,%%eax			\n\t"\
		"movl		%%eax,%%ebx			\n\t"\
		"andl		%[__nrtm1],%%eax	\n\t"\
		"shrl		%%cl,%%ebx			\n\t"\
		"shll		$4,%%eax			\n\t"\
		"shll		$4,%%ebx			\n\t"\
		"addl		%[__add1],%%eax		\n\t"\
		"addl		%[__add2],%%ebx		\n\t"\
		"movaps		(%%eax),%%xmm0		\n\t"\
		"movaps		(%%ebx),%%xmm1		\n\t"\
		"movl	%[__idx_offset],%%eax	\n\t"\
		"movaps		%%xmm1,%%xmm2		\n\t"\
		"shufpd	$1,	%%xmm2,%%xmm2		\n\t"\
		"mulpd		%%xmm0,%%xmm1		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"addl		$2,%%eax			\n\t"\
		"shrl		$1,%%eax			\n\t"\
		"movl		%%eax,%%ebx			\n\t"\
		"andl		%[__nrtm1],%%eax	\n\t"\
		"shrl		%%cl,%%ebx			\n\t"\
		"shll		$4,%%eax			\n\t"\
		"shll		$4,%%ebx			\n\t"\
		"addl		%[__add1],%%eax		\n\t"\
		"addl		%[__add2],%%ebx		\n\t"\
		"movaps		(%%eax),%%xmm0		\n\t"/* eax free, since store sign_mask into xmm-reg the 2x we use it below */\
		"movaps		(%%ebx),%%xmm3		\n\t"\
		"movaps		%%xmm3,%%xmm4		\n\t"/* __idx_offset was in eax, no more */\
		"shufpd	$1,	%%xmm4,%%xmm4		\n\t"\
		"mulpd		%%xmm0,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"movaps		%%xmm1,%%xmm0		\n\t"\
		"movl		%[__data],%%ecx		\n\t"\
		"unpcklpd	%%xmm3,%%xmm0		\n\t"\
		"unpckhpd	%%xmm3,%%xmm1		\n\t"\
		"subpd		%%xmm1,%%xmm0		\n\t"\
		"movaps		%%xmm2,%%xmm1		\n\t"\
		"unpcklpd	%%xmm4,%%xmm1		\n\t"\
		"unpckhpd	%%xmm4,%%xmm2		\n\t"\
		"addpd		%%xmm2,%%xmm1		\n\t"\
		"movaps		     (%%ecx),%%xmm4	\n\t"\
		"movaps		 0x10(%%ecx),%%xmm2	\n\t"\
		"movl		%[__half_arr],%%ecx	\n\t"\
		"movl	%[__odd_radix],%%eax	\n\t"/* In this 32-bit lo-reg version assume __odd_radix already lshifted << 4 on input */\
		"addl		%[__offset0],%%ecx	\n\t"\
		"movaps	(%%ecx,%%eax),%%xmm5	\n\t"\
		"subl		%[__offset0],%%ecx	\n\t"\
		"mulpd		%%xmm5,%%xmm4		\n\t"\
		"mulpd		%%xmm5,%%xmm2		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t"\
		"addpd		%%xmm3,%%xmm4		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t"\
		"movl		%[__cy],%%ebx		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"shufpd	$0,	%%xmm2,%%xmm4		\n\t"\
		"shufpd	$3,	%%xmm2,%%xmm5		\n\t"\
		"addpd		     (%%ebx),%%xmm4	\n\t"\
		"movaps		-0x20(%%ecx),%%xmm6	\n\t"\
		"movaps		-0x10(%%ecx),%%xmm7	\n\t"\
		"addl	   %[__offset0],%%ecx	\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
	/*========= use xmm1 to store (%[__sign_mask]) =========*/\
	"movl	%[__data],%%eax			\n\t"\
	"movaps		%%xmm1,(%%eax)		\n\t"/* spill xmm1 into in-data slot (which is otherwise used only at start and end) */\
	"movl	%[__sign_mask],%%eax	\n\t"\
	"movaps		(%%eax),%%xmm1		\n\t"\
	"movl	%[__odd_radix],%%eax	\n\t"/* restore odd_radix to eax */\
		"shll	   $1,%%eax				\n\t"\
		"addpd		%%xmm7,%%xmm4		\n\t"\
		"subpd		%%xmm7,%%xmm4		\n\t"\
		"subpd		%%xmm4,%%xmm2		\n\t"\
		"addl		%%eax,%%ecx			\n\t"\
		"andpd		%%xmm1,%%xmm2		\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t"\
		"shrl		$1,%%eax			\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"mulpd  (%%ecx,%%eax),%%xmm2	\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t"\
		"subpd		%%xmm7,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		(%%ecx),%%xmm3		\n\t"\
		"subl		%[__offset0],%%ecx	\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t"\
		"addpd		%%xmm7,%%xmm5		\n\t"\
		"subpd		%%xmm7,%%xmm5		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t"\
		"andpd		%%xmm1,%%xmm2		\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t"\
		"addl		%[__offset1],%%ecx	\n\t"\
		"mulpd  (%%ecx,%%eax),%%xmm2	\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t"\
		"subpd		%%xmm7,%%xmm2		\n\t"\
		"shll		$1,%%eax			\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		(%%ecx),%%xmm3		\n\t"\
		"subl		%[__offset1],%%ecx	\n\t"\
		"subpd		%%xmm3,%%xmm5		\n\t"\
		"movaps		%%xmm2,(%%ebx)		\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"shufpd	$0,	%%xmm5,%%xmm4		\n\t"\
		"shufpd	$3,	%%xmm5,%%xmm2		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"subl		%%eax,%%ecx			\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
	"movl	%[__data],%%eax			\n\t"\
	"movaps		(%%eax),%%xmm1		\n\t"/* restore earlier spill of xmm1 */\
		"mulpd		%%xmm1,%%xmm3		\n\t"\
		"movaps	 	%%xmm6,-0x20(%%ecx)	\n\t"\
		"addl		%[__offset0],%%ecx	\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t"\
		"movaps		(%%ecx),%%xmm0		\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm0,%%xmm5		\n\t"\
		"movaps		%%xmm4,    (%%eax)	\n\t"/* Write outputs to __data */\
		"movaps		%%xmm5,0x10(%%eax)	\n\t"\
		"movl	%[__idx_offset],%%ecx	\n\t"\
		"addl	%[__idx_incr],%%ecx		\n\t"\
		"movl	%%ecx, %[__idx_offset]	\n\t"/* Store incremented idx_offset */\
	"popl %%ebx	\n\t"\
	:						/* outputs: none */\
	:	[__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
	,	[__cy]			"m" (Xcy)\
	,	[__nrt_bits]	"m" (Xnrt_bits)\
	,	[__nrtm1]		"m" (Xnrtm1)\
	,	[__idx_offset]	"m" (Xidx_offset)\
	,	[__idx_incr]	"m" (Xidx_incr)\
	,	[__odd_radix]   "m" (Xodd_radix)\
	,	[__half_arr]	"m" (Xhalf_arr)\
	,	[__sign_mask]	"m" (Xsign_mask)\
	,	[__add1]		"m" (Xadd1)\
	,	[__add2]		"m" (Xadd2)\
	,	[__offset0]		"m" (Xoffset0)\
	,	[__offset1]		"m" (Xoffset1)\
	/* Prefetch address */\
	,	[__add0] "m" (Xadd0)\
		: "cc","memory","eax",/*"ebx",*/"ecx","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

#else

	#define SSE2_fermat_carry_norm_errcheck(Xdata,Xcy,Xnrt_bits,Xnrtm1,Xidx_offset,Xidx_incr,Xodd_radix,Xhalf_arr,Xsign_mask,Xadd1,Xadd2,Xoffset0,Xoffset1, Xadd0)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"movl		%[__add0],%%ecx		\n\t"/* base address for prefetch-from-main-data-array */\
	"prefetcht0	(%%ecx)	\n\t"\
		"movl	%[__idx_offset],%%esi	\n\t"\
		"movl %[__odd_radix],%%edi		\n\t"\
		"movl	%[__nrt_bits],%%ecx		\n\t"\
		"movl		%%esi,%%eax			\n\t"\
		"shrl		$1,%%eax			\n\t"\
		"movl		%%eax,%%ebx			\n\t"\
		"andl		%[__nrtm1],%%eax	\n\t"\
		"shrl		%%cl,%%ebx			\n\t"\
		"shll		$4,%%eax			\n\t"\
		"shll		$4,%%ebx			\n\t"\
		"shll		$4,%%edi			\n\t"\
		"addl		%[__add1],%%eax		\n\t"\
		"addl		%[__add2],%%ebx		\n\t"\
		"movaps		(%%eax),%%xmm0		\n\t"\
		"movaps		(%%ebx),%%xmm1		\n\t"\
		"movl		%%esi,%%eax			\n\t"\
		"movaps		%%xmm1,%%xmm2		\n\t"\
		"shufpd	$1,	%%xmm2,%%xmm2		\n\t"\
		"mulpd		%%xmm0,%%xmm1		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"addl		$2,%%eax			\n\t"\
		"shrl		$1,%%eax			\n\t"\
		"movl		%%eax,%%ebx			\n\t"\
		"andl		%[__nrtm1],%%eax	\n\t"\
		"shrl		%%cl,%%ebx			\n\t"\
		"shll		$4,%%eax			\n\t"\
		"shll		$4,%%ebx			\n\t"\
		"addl		%[__add1],%%eax		\n\t"\
		"addl		%[__add2],%%ebx		\n\t"\
		"movaps		(%%eax),%%xmm0		\n\t"\
		"movaps		(%%ebx),%%xmm3		\n\t"\
		"movl		%%esi,%%eax			\n\t"\
		"movaps		%%xmm3,%%xmm4		\n\t"\
		"shufpd	$1,	%%xmm4,%%xmm4		\n\t"\
		"mulpd		%%xmm0,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"movaps		%%xmm1,%%xmm0		\n\t"\
		"unpcklpd	%%xmm3,%%xmm0		\n\t"\
		"unpckhpd	%%xmm3,%%xmm1		\n\t"\
		"subpd		%%xmm1,%%xmm0		\n\t"\
		"movaps		%%xmm2,%%xmm1		\n\t"\
		"unpcklpd	%%xmm4,%%xmm1		\n\t"\
		"unpckhpd	%%xmm4,%%xmm2		\n\t"\
		"addpd		%%xmm2,%%xmm1		\n\t"\
		"movl		%[__half_arr],%%ecx	\n\t"\
		"movl		%[__data],%%edx		\n\t"\
		"movaps		     (%%edx),%%xmm4	\n\t"\
		"movaps		 0x10(%%edx),%%xmm2	\n\t"\
		"addl		%[__offset0],%%ecx	\n\t"\
		"movaps	(%%ecx,%%edi),%%xmm5	\n\t"\
		"subl		%[__offset0],%%ecx	\n\t"\
		"mulpd		%%xmm5,%%xmm4		\n\t"\
		"mulpd		%%xmm5,%%xmm2		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t"\
		"addpd		%%xmm3,%%xmm4		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t"\
		"movl		%[__cy],%%ebx		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"shufpd	$0,	%%xmm2,%%xmm4		\n\t"\
		"shufpd	$3,	%%xmm2,%%xmm5		\n\t"\
		"addpd		     (%%ebx),%%xmm4	\n\t"\
		"movaps		-0x20(%%ecx),%%xmm6	\n\t"\
		"movaps		-0x10(%%ecx),%%xmm7	\n\t"\
		"addl	   %[__offset0],%%ecx	\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"shll	   $1,%%edi				\n\t"\
		"addpd		%%xmm7,%%xmm4		\n\t"\
		"subpd		%%xmm7,%%xmm4		\n\t"\
		"movl		%[__sign_mask],%%eax\n\t"\
		"subpd		%%xmm4,%%xmm2		\n\t"\
		"andpd		(%%eax),%%xmm2		\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t"\
		"addl		%%edi,%%ecx			\n\t"\
		"shrl		$1,%%edi			\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"mulpd  (%%ecx,%%edi),%%xmm2	\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t"\
		"subpd		%%xmm7,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		(%%ecx),%%xmm3		\n\t"\
		"subl		%[__offset0],%%ecx	\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t"\
		"addpd		%%xmm7,%%xmm5		\n\t"\
		"subpd		%%xmm7,%%xmm5		\n\t"\
		"subpd		%%xmm5,%%xmm2		\n\t"\
		"andpd		(%%eax),%%xmm2		\n\t"\
		"maxpd		%%xmm6,%%xmm2		\n\t"\
		"movaps		%%xmm2,%%xmm6		\n\t"\
		"movaps		%%xmm5,%%xmm2		\n\t"\
		"addl		%[__offset1],%%ecx	\n\t"\
		"mulpd  (%%ecx,%%edi),%%xmm2	\n\t"\
		"addpd		%%xmm7,%%xmm2		\n\t"\
		"subpd		%%xmm7,%%xmm2		\n\t"\
		"shll		$1,%%edi			\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"mulpd		(%%ecx),%%xmm3		\n\t"\
		"subl		%[__offset1],%%ecx	\n\t"\
		"subpd		%%xmm3,%%xmm5		\n\t"\
		"movaps		%%xmm2,(%%ebx)		\n\t"\
		"movaps		%%xmm4,%%xmm2		\n\t"\
		"shufpd	$0,	%%xmm5,%%xmm4		\n\t"\
		"shufpd	$3,	%%xmm5,%%xmm2		\n\t"\
		"movaps		%%xmm4,%%xmm5		\n\t"\
		"movaps		%%xmm2,%%xmm3		\n\t"\
		"subl		%%edi,%%ecx			\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm1,%%xmm3		\n\t"\
		"movaps	 	%%xmm6,-0x20(%%ecx)	\n\t"\
		"addl		%[__offset0],%%ecx	\n\t"\
		"mulpd		%%xmm0,%%xmm2		\n\t"\
		"mulpd		%%xmm1,%%xmm5		\n\t"\
		"movaps		(%%ecx),%%xmm0		\n\t"\
		"subpd		%%xmm3,%%xmm4		\n\t"\
		"addpd		%%xmm2,%%xmm5		\n\t"\
		"mulpd		%%xmm0,%%xmm4		\n\t"\
		"mulpd		%%xmm0,%%xmm5		\n\t"\
		"movaps		%%xmm4,    (%%edx)	\n\t"\
		"movaps		%%xmm5,0x10(%%edx)	\n\t"\
		"addl	%[__idx_incr],%%esi		\n\t"\
		/* Store incremented idx offsetL:  */\
		"movl	%%esi, %[__idx_offset]	\n\t"\
	"popl %%ebx	\n\t"\
	:						/* outputs: none */\
	:	[__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
	,	[__cy]			"m" (Xcy)\
	,	[__nrt_bits]	"m" (Xnrt_bits)\
	,	[__nrtm1]		"m" (Xnrtm1)\
	,	[__idx_offset]	"m" (Xidx_offset)\
	,	[__idx_incr]	"m" (Xidx_incr)\
	,	[__odd_radix]   "m" (Xodd_radix)\
	,	[__half_arr]	"m" (Xhalf_arr)\
	,	[__sign_mask]	"m" (Xsign_mask)\
	,	[__add1]		"m" (Xadd1)\
	,	[__add2]		"m" (Xadd2)\
	,	[__offset0]		"m" (Xoffset0)\
	,	[__offset1]		"m" (Xoffset1)\
	/* Prefetch address */\
	,	[__add0] "m" (Xadd0)\
		: "cc","memory","eax",/*"ebx",*/"ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

#endif

	/*************************************************************/
	/**************** MERSENNE-MOD CARRY MACROS ******************/
	/*************************************************************/

	/********* Packed 32-bit-int version of SSE2_cmplx_carry_norm_pow2_errcheck0_2x:***********/
	#define SSE2_cmplx_carry_norm_pow2_errcheck1_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xi,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw, Xadd0,Xp1)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"movl		%[__add0],%%ecx		\n\t"/* base address for prefetch-from-main-data-array */\
	"prefetcht0	(%%ecx)	\n\t"\
	"/***************Unpack the data:*************************/\n\t"\
		"movl	%[__data]	,%%eax	\n\t"\
		"movaps		    (%%eax)	,%%xmm1	\n\t	movaps		0x40(%%eax)	,%%xmm5	\n\t"\
		"movaps		    (%%eax)	,%%xmm2	\n\t	movaps		0x40(%%eax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%eax)	,%%xmm1	\n\t	unpcklpd	0x60(%%eax)	,%%xmm5	\n\t"\
		"unpckhpd	0x20(%%eax)	,%%xmm2	\n\t	unpckhpd	0x60(%%eax)	,%%xmm6	\n\t"\
		"movaps		%%xmm2, 0x20(%%eax)	\n\t	movaps		%%xmm6, 0x60(%%eax)	\n\t"\
		"\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2	\n\t	movaps		0x50(%%eax)	,%%xmm6	\n\t"\
		"movaps		0x10(%%eax)	,%%xmm3	\n\t	movaps		0x50(%%eax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%eax)	,%%xmm2	\n\t	unpcklpd	0x70(%%eax)	,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%eax)	,%%xmm3	\n\t	unpckhpd	0x70(%%eax)	,%%xmm7	\n\t"\
		"movaps		%%xmm2, 0x10(%%eax)	\n\t	movaps		%%xmm6, 0x50(%%eax)	\n\t"\
		"movaps		%%xmm3, 0x30(%%eax)	\n\t	movaps		%%xmm7, 0x70(%%eax)	\n\t"\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__bjmod_0],	%%eax	\n\t"\
		"movaps		(%%eax)	,	%%xmm0	\n\t"\
		"movl	%[__sse_sw]	,	%%ebx	\n\t"\
		"movaps		(%%ebx)	,	%%xmm7	\n\t"\
		"psubd		%%xmm0	,	%%xmm7	\n\t"\
		"movmskps	%%xmm7	,	%%esi	\n\t"\
	"movl	%[__i]	,%%ecx			\n\t"/* I == 1 if it's the 0-word, in which case we force-bigword-ness by XORing esi (whose */\
	"xorl	%%ecx	,%%esi			\n\t"/* low bit will == 0 on input in this case) with I. Otherwise I == 0, thus XOR = no-op. */\
		"shll	$24		,%%esi			\n\t"\
		"movaps		%%xmm0	,%%xmm7		\n\t"\
		"movl	%[__n_minus_sil],%%ecx	\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16		,%%ecx			\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwt]	,%%edx		\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm7		\n\t"\
		"movmskps	%%xmm7	,%%edx		\n\t"\
		"shll	$8		,%%edx			\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax		\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"movaps		-0x10(%%eax),%%xmm4	\n\t"\
		"movaps		     (%%ebx),%%xmm2	\n\t	movaps		 0x10(%%ebx),%%xmm6	\n\t"\
		"movhpd		     (%%ecx),%%xmm3	\n\t	movhpd		-0x10(%%ecx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%ecx),%%xmm3	\n\t	movlpd		-0x08(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ecx	,%[__wtB]		\n\t"\
		"\n\t"\
		"shld	$20,	%%esi	,%%edi	\n\t	shld	$18,	%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%eax),%%xmm2	\n\t	mulpd	 	0x100(%%eax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%eax),%%xmm3	\n\t	mulpd	 	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	 	     (%%edi),%%xmm2	\n\t	mulpd	 	     (%%ebx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%edx),%%xmm3	\n\t	mulpd	 	0x040(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%ecx)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%eax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"movaps		%%xmm3	,(%%ecx)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%eax)\n\t	movaps		%%xmm5	,0x40(%%eax)\n\t"\
		"\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"movl	%[__sse_nm1]	,%%ebx	\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"pand		(%%ebx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
	"movl		%[__add0],%%ecx		\n\t"/* base address for prefetch-from-main-data-array */\
	"movl		%[__p1],%%esi	\n\t"	\
	"prefetcht0	(%%ecx,%%esi,8)	\n\t"\
		"movl	%[__sse_sw]	,%%edx		\n\t"\
		"movaps	(%%edx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%esi		\n\t"\
		"shll	$24	,%%esi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movl	%[__n_minus_silp1],%%ecx\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16	,%%ecx				\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwtm1]	,%%edx	\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8	,%%edx				\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"movaps	 0x10(%%eax)	,%%xmm1	\n\t"\
		"movaps	 0x50(%%eax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtC]	,%%ecx		\n\t"\
		"\n\t"\
		"movaps	     (%%ebx)	,%%xmm2	\n\t	movaps	 0x10(%%ebx)	,%%xmm6	\n\t"\
		"movhpd	     (%%ecx)	,%%xmm3	\n\t	movhpd	-0x10(%%ecx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%ecx)	,%%xmm3	\n\t	movlpd	-0x08(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addl	$0x20	,%%ebx			\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ebx	,%[__wtA]		\n\t"\
		"movl	%%ecx	,%[__wtC]		\n\t"\
		"shld	$20		,%%esi	,%%edi	\n\t	shld	$18		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%eax)	,%%xmm2	\n\t	mulpd	 0x120(%%eax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%eax)	,%%xmm3	\n\t	mulpd	 0x130(%%eax)	,%%xmm7	\n\t"\
		"mulpd	      (%%edi)	,%%xmm2	\n\t	mulpd	      (%%ebx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%edx)	,%%xmm3	\n\t	mulpd	 0x040(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%ecx)	,%%xmm1			\n\t	addpd	(%%edx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"andpd	(%%ebx)	,%%xmm1			\n\t	andpd	(%%ebx)	,%%xmm5			\n\t"\
		"maxpd	%%xmm5	,%%xmm1			\n\t"\
		"maxpd	-0x20(%%eax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"movaps	%%xmm3	,(%%ecx)		\n\t	movaps	%%xmm7	,(%%edx)		\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%eax)	\n\t	movaps	%%xmm5	, 0x50(%%eax)	\n\t"\
		"\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"movl	%[__sse_nm1]	,%%ebx	\n\t"\
		"paddd	(%%eax)	,%%xmm0			\n\t"\
		"pand	(%%ebx)	,%%xmm0			\n\t"\
		"movl	%[__bjmod_0],%%ecx		\n\t"\
		"movaps	%%xmm0,(%%ecx)			\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__wtC]		"m" (XwtC)		\
		, [__cyA]		"m" (XcyA)		\
		, [__cyB]		"m" (XcyB)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__i]			"m" (Xi)			\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sign_mask]	"m" (Xsign_mask)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]		"m" (Xsinwtm1)		\
		, [__sse_bw]		"m" (Xsse_bw)		\
		, [__sse_nm1]		"m" (Xsse_nm1)		\
		, [__sse_sw]		"m" (Xsse_sw)		\
	/* Prefetch: base address and 1 index offsets */\
	,	[__add0] "m" (Xadd0)\
	,	[__p1] "m" (Xp1)\
		: "cc","memory","eax","ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	/********* Packed 32-bit-int version of SSE2_cmplx_carry_norm_pow2_errcheck2_2x:***********/
	#define SSE2_cmplx_carry_norm_pow2_errcheck2_2B(Xdata,XwtA,XwtB,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_nm1,Xsse_sw, Xadd0,Xp2,Xp3)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"movl		%[__add0],%%ecx		\n\t"/* base address for prefetch-from-main-data-array */\
	"movl		%[__p2],%%esi	\n\t"	\
	"prefetcht0	(%%ecx,%%esi,8)	\n\t"\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__bjmod_0],	%%eax	\n\t"\
		"movaps		(%%eax)	,	%%xmm0	\n\t"\
		"movl	%[__sse_sw]	,	%%ebx	\n\t"\
		"movaps		(%%ebx)	,	%%xmm1	\n\t"\
		"psubd		%%xmm0	,	%%xmm1	\n\t"\
		"movmskps	%%xmm1	,	%%esi	\n\t"\
		"\n\t"\
		"shll	$24		,%%esi			\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movl	%[__n_minus_sil],%%ecx	\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16		,%%ecx			\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwt]	,%%edx		\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm1		\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8		,%%edx			\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"\n\t"\
		"movl	%[__data],%%eax			\n\t"\
		"movaps	 0x20(%%eax),%%xmm1		\n\t	movaps		 0x60(%%eax),%%xmm5\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"movaps		-0x10(%%eax),%%xmm4	\n\t"\
		"movaps		     (%%ebx),%%xmm2	\n\t	movaps		 0x10(%%ebx),%%xmm6	\n\t"\
		"movhpd		     (%%ecx),%%xmm3	\n\t	movhpd		-0x10(%%ecx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%ecx),%%xmm3	\n\t	movlpd		-0x08(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"shld	$20,	%%esi	,%%edi	\n\t	shld	$18,	%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 	0x110(%%eax),%%xmm3	\n\t	mulpd	 	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	 	0x100(%%eax),%%xmm2	\n\t	mulpd	 	0x100(%%eax),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%edx),%%xmm3	\n\t	mulpd	 	0x040(%%ecx),%%xmm7	\n\t"\
		"mulpd	 	     (%%edi),%%xmm2	\n\t	mulpd	 	     (%%ebx),%%xmm6	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%ecx)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%eax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"movaps		%%xmm3	,(%%ecx)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,0x20(%%eax)\n\t	movaps		%%xmm5	,0x60(%%eax)\n\t"\
		"\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"movl	%[__sse_nm1]	,%%ebx	\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"pand		(%%ebx)	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
	"movl		%[__add0],%%ecx		\n\t"/* base address for prefetch-from-main-data-array */\
	"movl		%[__p3],%%esi	\n\t"	\
	"prefetcht0	(%%ecx,%%esi,8)	\n\t"\
		"movl	%[__sse_sw]	,%%edx		\n\t"\
		"movaps	(%%edx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%esi		\n\t"\
		"shll	$24	,%%esi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movl	%[__n_minus_silp1],%%ecx\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16	,%%ecx				\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwtm1]	,%%edx	\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8	,%%edx				\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"movaps	 0x30(%%eax)	,%%xmm1	\n\t"\
		"movaps	 0x70(%%eax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"\n\t"\
		"movaps	     (%%ebx)	,%%xmm2	\n\t	movaps	 0x10(%%ebx)	,%%xmm6	\n\t"\
		"movhpd	     (%%ecx)	,%%xmm3	\n\t	movhpd	-0x10(%%ecx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%ecx)	,%%xmm3	\n\t	movlpd	-0x08(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addl	$0x20	,%%ebx			\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ebx	,%[__wtA]		\n\t"\
		"movl	%%ecx	,%[__wtB]		\n\t"\
		"shld	$20		,%%esi	,%%edi	\n\t	shld	$18		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 0x130(%%eax)	,%%xmm3	\n\t	mulpd	 0x130(%%eax)	,%%xmm7	\n\t"\
		"mulpd	 0x120(%%eax)	,%%xmm2	\n\t	mulpd	 0x120(%%eax)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%edx)	,%%xmm3	\n\t	mulpd	 0x040(%%ecx)	,%%xmm7	\n\t"\
		"mulpd	      (%%edi)	,%%xmm2	\n\t	mulpd	      (%%ebx)	,%%xmm6	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%ecx)	,%%xmm1			\n\t	addpd	(%%edx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"andpd	(%%ebx)	,%%xmm1			\n\t	andpd	(%%ebx)	,%%xmm5			\n\t"\
		"maxpd	%%xmm5	,%%xmm1			\n\t"\
		"maxpd	-0x20(%%eax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"movaps	%%xmm3	,(%%ecx)		\n\t	movaps	%%xmm7	,(%%edx)		\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x30(%%eax)	\n\t	movaps	%%xmm5	, 0x70(%%eax)	\n\t"\
		"\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"movl	%[__sse_nm1]	,%%ebx	\n\t"\
		"paddd	(%%eax)	,%%xmm0			\n\t"\
		"pand	(%%ebx)	,%%xmm0			\n\t"\
		"movl	%[__bjmod_0],%%ecx		\n\t"\
		"movaps	%%xmm0,(%%ecx)			\n\t"\
	"/**********************************************/\n\t"\
	"/*          Repack the data:              */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__data],%%eax			\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1	\n\t	movaps		0x50(%%eax)	,%%xmm5	\n\t"\
		"movaps		    (%%eax)	,%%xmm0	\n\t	movaps		0x40(%%eax)	,%%xmm4	\n\t"\
		"movaps		%%xmm1		,%%xmm3	\n\t	movaps		%%xmm5		,%%xmm7	\n\t"\
		"movaps		%%xmm0		,%%xmm2	\n\t	movaps		%%xmm4		,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%eax)	,%%xmm3	\n\t	unpckhpd	0x70(%%eax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%eax)	,%%xmm1	\n\t	unpcklpd	0x70(%%eax)	,%%xmm5	\n\t"\
		"movaps		%%xmm3,0x30(%%eax)	\n\t	movaps		%%xmm7,0x70(%%eax)	\n\t"\
		"unpckhpd	0x20(%%eax)	,%%xmm2	\n\t	unpckhpd	0x60(%%eax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%eax)	,%%xmm0	\n\t	unpcklpd	0x60(%%eax)	,%%xmm4	\n\t"\
		"movaps		%%xmm2,0x20(%%eax)	\n\t	movaps		%%xmm6,0x60(%%eax)	\n\t"\
		"movaps		%%xmm1,0x10(%%eax)	\n\t	movaps		%%xmm5,0x50(%%eax)	\n\t"\
		"movaps		%%xmm0,    (%%eax)	\n\t	movaps		%%xmm4,0x40(%%eax)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__cyA]		"m" (XcyA)		\
		, [__cyB]		"m" (XcyB)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sign_mask]	"m" (Xsign_mask)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]		"m" (Xsinwtm1)		\
		, [__sse_bw]		"m" (Xsse_bw)		\
		, [__sse_nm1]		"m" (Xsse_nm1)		\
		, [__sse_sw]		"m" (Xsse_sw)		\
	/* Prefetch: base address and 2 index offsets */\
	,	[__add0] "m" (Xadd0)\
	,	[__p2] "m" (Xp2)\
	,	[__p3] "m" (Xp3)\
		: "cc","memory","eax","ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	/********* Fused fast-LOACC-algo [cf. carry.h::cmplx_carry_fast_pow2_errcheck] version of above 2 macros: ***********/
	// "Fused" means that - like the HIACC macros - we process 4 carry chains, one from each separate array section corr.
	// to each wide-strided final-iFFT-pass output, at a time, but fuse the [j,j+2] linear-index-within-each-array-section
	// processing (done separately in the HIACC case by the 1_2B and 2_2B SSE2 carry macros) into a single macro. This
	// fusion is eased by the fact that the LOACC chained-weights-computation needs no weights-reinit-from-scalar-data
	// step for the [j+2] data.
	//
	#define SSE2_cmplx_carry_fast_pow2_errcheck(Xdata,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xi,Xsign_mask,Xsse_bw,Xsse_nm1,Xsse_sw, Xadd0,Xp1,Xp2,Xp3)\
	{\
	__asm__ volatile (\
		"movl	%[__add0],%%esi	\n\t"/* base address for 4 prefetches-from-main-data-array spread through this macro */\
		"prefetcht0	(%%esi)		\n\t"\
		/***************Unpack the data:*************************/\
			"movl	%[__data]	,%%eax	\n\t"\
			"movaps		    (%%eax)	,%%xmm1	\n\t	movaps		0x40(%%eax)	,%%xmm5	\n\t"/* r1, this is the active  xmm register */\
			"movaps		    %%xmm1	,%%xmm2	\n\t	movaps			%%xmm5	,%%xmm6	\n\t"/* r1, this is the scratch xmm register */\
			"movaps		0x20(%%eax)	,%%xmm0	\n\t	movaps		0x60(%%eax)	,%%xmm4	\n\t"\
			"unpcklpd		%%xmm0	,%%xmm1	\n\t	unpcklpd		%%xmm4	,%%xmm5	\n\t"/* r1 -x- r3 (lo halves) ==> R0~ */\
			"unpckhpd		%%xmm0	,%%xmm2	\n\t	unpckhpd		%%xmm4	,%%xmm6	\n\t"/* r1 -x- r3 (hi halves) ==> R1~ */\
			"movaps		%%xmm2, 0x20(%%eax)	\n\t	movaps		%%xmm6, 0x60(%%eax)	\n\t"/* Tmp store R1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
			"\n\t"\
			"movaps		0x10(%%eax)	,%%xmm2	\n\t	movaps		0x50(%%eax)	,%%xmm6	\n\t"\
			"movaps			%%xmm2	,%%xmm3	\n\t	movaps			%%xmm6	,%%xmm7	\n\t"\
			"movaps		0x30(%%eax)	,%%xmm0	\n\t	movaps		0x70(%%eax)	,%%xmm4	\n\t"\
			"unpcklpd		%%xmm0	,%%xmm2	\n\t	unpcklpd		%%xmm4	,%%xmm6	\n\t"/* r2 -x- r4 (lo halves) ==> I0~ */\
			"unpckhpd		%%xmm0	,%%xmm3	\n\t	unpckhpd		%%xmm4	,%%xmm7	\n\t"/* r2 -x- r4 (hi halves) ==> I1~ */\
			"movaps		%%xmm2, 0x10(%%eax)	\n\t	movaps		%%xmm6, 0x50(%%eax)	\n\t"/* Tmp store I0~ until needed by imaginary-part-processing section */\
			"movaps		%%xmm3, 0x30(%%eax)	\n\t	movaps		%%xmm7, 0x70(%%eax)	\n\t"/* Tmp store I1~ until needed on 2nd set of SSE2_cmplx_carry.calls */\
		/* Active data in xmm1,5 here - avoid using those registers in index computation. */\
		/**********************************************/\
		/*          Real      parts                   */\
		/**********************************************/\
			"movl	%[__bjmod_0],	%%eax	\n\t"/* Pointer to bjmodn data */\
			"movaps		(%%eax)	,	%%xmm0	\n\t"/* bjmod[0:3] */\
			"movl	%[__sse_sw]	,	%%ebx	\n\t"\
			"movaps		(%%ebx)	,	%%xmm7	\n\t"/* sw[0:3] */\
			"psubd		%%xmm0	,	%%xmm7	\n\t"/* sw[0:3] - bjmod[0:3] */\
			"movmskps	%%xmm7	,	%%esi	\n\t"/* Extract sign bits into 4-bit signmask */\
			"movl	%[__i]	,%%edi			\n\t"/* I == 1 if it's the 0-word, in which case we force-bigword-ness by XORing esi (whose */\
			"xorl	%%edi	,%%esi			\n\t"/* low bit will == 0 on input in this case) with I. Otherwise I == 0, thus XOR = no-op. */\
			"movl	%[__half_arr]	,%%eax	\n\t"\
			"movaps	0x100(%%eax),%%xmm2		\n\t	movaps	 0x110(%%eax),%%xmm6	\n\t"/* wt_re for our 2 independent carry-chain pairs */\
			"movaps	-0x10(%%eax),%%xmm4		\n\t"/* sse2_rnd */\
			"movl	%[__cyA]	,%%edi		\n\t	movl	%[__cyB]	,%%edx		\n\t"/* cy_in */\
			"mulpd	0x120(%%eax),%%xmm1		\n\t	mulpd	0x130(%%eax),%%xmm5		\n\t"/* x = x*wi_re */\
			"addpd		(%%edi)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"/* x = x*wi_re + cy */\
			"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"/* temp = x */\
			"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
			"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"/* temp = DNINT(x) */\
			"movl	%[__sign_mask],%%ebx	\n\t"\
			"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"/* x - temp */\
			"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"/* frac = fabs(x-temp) */\
			"maxpd			%%xmm5	,%%xmm1	\n\t	maxpd		-0x20(%%eax),%%xmm1	\n\t"/* if(frac > maxerr) maxerr=frac */\
			"movaps		%%xmm1,-0x20(%%eax)	\n\t"/* Note serialization here! */\
			"movl	%%esi,	%%ebx													\n\t"/* byte offsets for bits <0:1> */\
			"shll	 $4,	%%ebx			\n\t	shll	 $2,	%%esi			\n\t"/* and <2:3> go into rbx, rsi */\
			"andl	$0x30,	%%ebx			\n\t	andl	$0x30,	%%esi			\n\t"/* respectively. */\
			"addl	%%eax,	%%ebx			\n\t	addl	%%eax,	%%esi			\n\t"\
			"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"/* cpy temp */\
			"mulpd	 	0xc0(%%ebx),%%xmm3	\n\t	mulpd	 	0xc0(%%esi),%%xmm7	\n\t"/* temp*baseinv[i] */\
			"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
			"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"/* cy_out */\
			"movaps		%%xmm3	,(%%edi)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"/* store cy_out */\
			"movl	%[__data]	,%%eax		\n\t"\
			"mulpd	 0x80(%%ebx),	%%xmm3	\n\t	mulpd	 0x80(%%esi),	%%xmm7	\n\t"/* cy*base[i] */\
			"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"/* x = (temp-cy*base[i]) */\
			"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"/* x*= wt_re */\
			"movaps		%%xmm1	,    (%%eax)\n\t	movaps		%%xmm5	,0x40(%%eax)\n\t"/* store x */\
		/* Update and store weights: */\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movaps	0x120(%%eax),%%xmm3		\n\t	movaps	 0x130(%%eax),%%xmm7	\n\t"/* wi_re, inverse-wt_re */\
		"movaps	 0x70(%%eax),%%xmm1		\n\t	movaps		%%xmm1,%%xmm5		\n\t"/* inv_mult[1], 2 copies */\
		/* Do compare as (inv_mult[1] < wt_re), result overwrites copy of inv_mult[1]: */\
		"cmppd	$1,%%xmm2,%%xmm1		\n\t	cmppd	$1,%%xmm6,%%xmm5		\n\t"/* i = (wt_re >= inv_mult[1]) */\
		"movmskpd	%%xmm1,%%edi		\n\t	movmskpd	%%xmm5,%%edx		\n\t"/* Extract cmp-results into pair of 2-bit signmasks */\
		"shll			$4,%%edi		\n\t	shll			$4,%%edx		\n\t"/* ...and mpy by xmm bytewidth. */\
		"leal	(%%eax,%%edi),%%edi		\n\t	leal	(%%eax,%%edx),%%edx		\n\t"/* address = half_arr + i */\
		"mulpd	     (%%edi),%%xmm2		\n\t	mulpd	     (%%edx),%%xmm6		\n\t"/* wt_re *= wts_mult[i] */\
		"mulpd	 0x40(%%edi),%%xmm3		\n\t	mulpd	 0x40(%%edx),%%xmm7		\n\t"/* wi_re *= inv_mult[i] */\
		"movaps	%%xmm2,0x100(%%eax)		\n\t	movaps	 %%xmm6,0x110(%%eax)	\n\t"/* Store wt_re */\
		"movaps	%%xmm3,0x120(%%eax)		\n\t	movaps	 %%xmm7,0x130(%%eax)	\n\t"/* Store wi_re */\
			/* Get ready for next set [IM0~] : */\
			"movl	%[__sse_bw]	,%%eax		\n\t"\
			"movl	%[__sse_nm1],%%ebx		\n\t"\
			"paddd		(%%eax)	,%%xmm0		\n\t"/* bjmod[0:3] += bw  */\
			"pand		(%%ebx)	,%%xmm0		\n\t"/* bjmod[0:3] &= nm1 */\
		/**********************************************/\
		/*          Imaginary parts                   */\
		/**********************************************/\
		"movl		%[__add0],%%esi	\n\t"\
		"movl		%[__p1],%%eax	\n\t"	\
		"prefetcht0	(%%esi,%%eax,8)	\n\t"\
			"movl	%[__sse_sw]	,	%%ebx	\n\t"\
			"movaps		(%%ebx)	,	%%xmm7	\n\t"/* sw[0:3] */\
			"psubd		%%xmm0	,	%%xmm7	\n\t"/* sw[0:3] - bjmod[0:3] */\
			"movmskps	%%xmm7	,	%%esi	\n\t"/* Extract sign bits into 4-bit signmask */\
			"movl	%[__data]	,%%eax		\n\t"\
			"movaps	 0x10(%%eax)	,%%xmm1	\n\t	movaps	 0x50(%%eax)	,%%xmm5	\n\t"/* I0~ */\
			"movl	%[__half_arr]	,%%eax	\n\t"\
			"movaps	0x140(%%eax),%%xmm2		\n\t	movaps	 0x150(%%eax),%%xmm6	\n\t"/* wt_re for our 2 independent carry-chain pairs */\
		/*	"movaps	-0x10(%%eax),%%xmm4		\n\t"// sse2_rnd already in xmm4 */\
			"movl	%[__cyA]	,%%edi		\n\t	movl	%[__cyB]	,%%edx		\n\t"/* cy_in */\
			"mulpd	0x160(%%eax),%%xmm1		\n\t	mulpd	0x170(%%eax),%%xmm5		\n\t"/* y = y*wi_im */\
			"addpd		(%%edi)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"/* y = y*wi_im + cy */\
			"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"/* temp = y */\
			"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
			"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"/* temp = DNINT(y) */\
			"movl	%[__sign_mask],%%ebx	\n\t"\
			"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"/* y - temp */\
			"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"/* frac = fabs(y-temp) */\
			"maxpd			%%xmm5	,%%xmm1	\n\t	maxpd		-0x20(%%eax),%%xmm1	\n\t"/* if(frac > maxerr) maxerr=frac */\
			"movaps		%%xmm1,-0x20(%%eax)	\n\t"/* Note serialization here! */\
			"movl	%%esi,	%%ebx													\n\t"/* byte offsets for bits <0:1> */\
			"shll	 $4,	%%ebx			\n\t	shll	 $2,	%%esi			\n\t"/* and <2:3> go into rbx, rsi */\
			"andl	$0x30,	%%ebx			\n\t	andl	$0x30,	%%esi			\n\t"/* respectively. */\
			"addl	%%eax,	%%ebx			\n\t	addl	%%eax,	%%esi			\n\t"\
			"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"/* cpy temp */\
			"mulpd	 	0xc0(%%ebx),%%xmm3	\n\t	mulpd	 	0xc0(%%esi),%%xmm7	\n\t"/* temp*baseinv[i] */\
			"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
			"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"/* cy_out */\
			"movaps		%%xmm3	,(%%edi)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"/* store cy_out */\
			"movl	%[__data]	,%%eax		\n\t"\
			"mulpd	 0x80(%%ebx),	%%xmm3	\n\t	mulpd	 0x80(%%esi),	%%xmm7	\n\t"/* cy*base[i] */\
			"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"/* y = (temp-cy*base[i]) */\
			"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"/* y*= wt_im */\
			"movaps		%%xmm1	,0x10(%%eax)\n\t	movaps		%%xmm5	,0x50(%%eax)\n\t"/* store y */\
		/* Update and store weights: */\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movaps	0x160(%%eax),%%xmm3		\n\t	movaps	 0x170(%%eax),%%xmm7	\n\t"/* wi_im, inverse-wt_im */\
		"movaps	 0x70(%%eax),%%xmm1		\n\t	movaps		%%xmm1,%%xmm5		\n\t"/* inv_mult[1], 2 copies */\
		/* Do compare as (inv_mult[1] < wt_im), result overwrites copy of inv_mult[1]: */\
		"cmppd	$1,%%xmm2,%%xmm1		\n\t	cmppd	$1,%%xmm6,%%xmm5		\n\t"/* i = (wt_im >= inv_mult[1]) */\
		"movmskpd	%%xmm1,%%edi		\n\t	movmskpd	%%xmm5,%%edx		\n\t"/* Extract cmp-results into pair of 2-bit signmasks */\
		"shll			$4,%%edi		\n\t	shll			$4,%%edx		\n\t"/* ...and mpy by xmm bytewidth. */\
		"leal	(%%eax,%%edi),%%edi		\n\t	leal	(%%eax,%%edx),%%edx		\n\t"/* address = half_arr + i */\
		"mulpd	     (%%edi),%%xmm2		\n\t	mulpd	     (%%edx),%%xmm6		\n\t"/* wt_im *= wts_mult[i] */\
		"mulpd	 0x40(%%edi),%%xmm3		\n\t	mulpd	 0x40(%%edx),%%xmm7		\n\t"/* wi_im *= inv_mult[i] */\
		"movaps	%%xmm2,0x140(%%eax)		\n\t	movaps	 %%xmm6,0x150(%%eax)	\n\t"/* Store wt_im */\
		"movaps	%%xmm3,0x160(%%eax)		\n\t	movaps	 %%xmm7,0x170(%%eax)	\n\t"/* Store wi_im */\
			/* Get ready for next set [RE1~] : */\
			"movl	%[__sse_bw]	,%%eax		\n\t"\
			"movl	%[__sse_nm1],%%ebx		\n\t"\
			"paddd		(%%eax)	,%%xmm0		\n\t"/* bjmod[0:3] += bw  */\
			"pand		(%%ebx)	,%%xmm0		\n\t"/* bjmod[0:3] &= nm1 */\
	/**********************************************/\
	/*          Now do the (j+2) data:            */\
	/**********************************************/\
		"movl		%[__add0],%%esi	\n\t"\
		"movl		%[__p2],%%eax	\n\t"	\
		"prefetcht0	(%%esi,%%eax,8)	\n\t"\
		/**********************************************/\
		/*          Real      parts                   */\
		/**********************************************/\
			"movl	%[__sse_sw]	,	%%ebx	\n\t"\
			"movaps		(%%ebx)	,	%%xmm7	\n\t"/* sw[0:3] */\
			"psubd		%%xmm0	,	%%xmm7	\n\t"/* sw[0:3] - bjmod[0:3] */\
			"movmskps	%%xmm7	,	%%esi	\n\t"/* Extract sign bits into 4-bit signmask */\
			"movl	%[__data]	,%%eax		\n\t"\
			"movaps	 0x20(%%eax)	,%%xmm1	\n\t	movaps	 0x60(%%eax)	,%%xmm5	\n\t"/* R1~ */\
			"movl	%[__half_arr]	,%%eax	\n\t"\
			"movaps	0x180(%%eax),%%xmm2		\n\t	movaps	 0x190(%%eax),%%xmm6	\n\t"/* wt_re for our 2 independent carry-chain pairs */\
		/*	"movaps	-0x10(%%eax),%%xmm4		\n\t"// sse2_rnd already in xmm4 */\
			"movl	%[__cyA]	,%%edi		\n\t	movl	%[__cyB]	,%%edx		\n\t"/* cy_in */\
			"mulpd	0x1a0(%%eax),%%xmm1		\n\t	mulpd	0x1b0(%%eax),%%xmm5		\n\t"/* x = x*wi_re */\
			"addpd		(%%edi)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"/* x = x*wi_re + cy */\
			"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"/* temp = x */\
			"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
			"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"/* temp = DNINT(x) */\
			"movl	%[__sign_mask],%%ebx	\n\t"\
			"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"/* x - temp */\
			"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"/* frac = fabs(x-temp) */\
			"maxpd			%%xmm5	,%%xmm1	\n\t	maxpd		-0x20(%%eax),%%xmm1	\n\t"/* if(frac > maxerr) maxerr=frac */\
			"movaps		%%xmm1,-0x20(%%eax)	\n\t"/* Note serialization here! */\
			"movl	%%esi,	%%ebx													\n\t"/* byte offsets for bits <0:1> */\
			"shll	 $4,	%%ebx			\n\t	shll	 $2,	%%esi			\n\t"/* and <2:3> go into rbx, rsi */\
			"andl	$0x30,	%%ebx			\n\t	andl	$0x30,	%%esi			\n\t"/* respectively. */\
			"addl	%%eax,	%%ebx			\n\t	addl	%%eax,	%%esi			\n\t"\
			"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"/* cpy temp */\
			"mulpd	 	0xc0(%%ebx),%%xmm3	\n\t	mulpd	 	0xc0(%%esi),%%xmm7	\n\t"/* temp*baseinv[i] */\
			"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
			"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"/* cy_out */\
			"movaps		%%xmm3	,(%%edi)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"/* store cy_out */\
			"movl	%[__data]	,%%eax		\n\t"\
			"mulpd	 0x80(%%ebx),	%%xmm3	\n\t	mulpd	 0x80(%%esi),	%%xmm7	\n\t"/* cy*base[i] */\
			"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"/* x = (temp-cy*base[i]) */\
			"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"/* x*= wt_re */\
			"movaps		%%xmm1	,0x20(%%eax)\n\t	movaps		%%xmm5	,0x60(%%eax)\n\t"/* store x */\
		/* Update and store weights ... the (j+2) wts are in the 8 slots (half_arr+[0x180-0x1f0]) above the (j) ones: */\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movaps	0x1a0(%%eax),%%xmm3		\n\t	movaps	 0x1b0(%%eax),%%xmm7	\n\t"/* wi_re, inverse-wt_re */\
		"movaps	 0x70(%%eax),%%xmm1		\n\t	movaps		%%xmm1,%%xmm5		\n\t"/* inv_mult[1], 2 copies */\
		/* Do compare as (inv_mult[1] < wt_re), result overwrites copy of inv_mult[1]: */\
		"cmppd	$1,%%xmm2,%%xmm1		\n\t	cmppd	$1,%%xmm6,%%xmm5		\n\t"/* i = (wt_re >= inv_mult[1]) */\
		"movmskpd	%%xmm1,%%edi		\n\t	movmskpd	%%xmm5,%%edx		\n\t"/* Extract cmp-results into pair of 2-bit signmasks */\
		"shll			$4,%%edi		\n\t	shll			$4,%%edx		\n\t"/* ...and mpy by xmm bytewidth. */\
		"leal	(%%eax,%%edi),%%edi		\n\t	leal	(%%eax,%%edx),%%edx		\n\t"/* address = half_arr + i */\
		"mulpd	     (%%edi),%%xmm2		\n\t	mulpd	     (%%edx),%%xmm6		\n\t"/* wt_re *= wts_mult[i] */\
		"mulpd	 0x40(%%edi),%%xmm3		\n\t	mulpd	 0x40(%%edx),%%xmm7		\n\t"/* wi_re *= inv_mult[i] */\
		"movaps	%%xmm2,0x180(%%eax)		\n\t	movaps	 %%xmm6,0x190(%%eax)	\n\t"/* Store wt_re */\
		"movaps	%%xmm3,0x1a0(%%eax)		\n\t	movaps	 %%xmm7,0x1b0(%%eax)	\n\t"/* Store wi_re */\
			/* Get ready for next set [IM0~] : */\
			"movl	%[__sse_bw]	,%%eax		\n\t"\
			"movl	%[__sse_nm1],%%ebx		\n\t"\
			"paddd		(%%eax)	,%%xmm0		\n\t"/* bjmod[0:3] += bw  */\
			"pand		(%%ebx)	,%%xmm0		\n\t"/* bjmod[0:3] &= nm1 */\
		/**********************************************/\
		/*          Imaginary parts                   */\
		/**********************************************/\
		"movl		%[__add0],%%esi	\n\t"\
		"movl		%[__p3],%%eax	\n\t"	\
		"prefetcht0	(%%esi,%%eax,8)	\n\t"\
			"movl	%[__sse_sw]	,	%%ebx	\n\t"\
			"movaps		(%%ebx)	,	%%xmm7	\n\t"/* sw[0:3] */\
			"psubd		%%xmm0	,	%%xmm7	\n\t"/* sw[0:3] - bjmod[0:3] */\
			"movmskps	%%xmm7	,	%%esi	\n\t"/* Extract sign bits into 4-bit signmask */\
			"movl	%[__data]	,%%eax		\n\t"\
			"movaps	 0x30(%%eax)	,%%xmm1	\n\t	movaps	 0x70(%%eax)	,%%xmm5	\n\t"/* I1~ */\
			"movl	%[__half_arr]	,%%eax	\n\t"\
			"movaps	0x1c0(%%eax),%%xmm2		\n\t	movaps	 0x1d0(%%eax),%%xmm6	\n\t"/* wt_im for our 2 independent carry-chain pairs */\
		/*	"movaps	-0x10(%%eax),%%xmm4		\n\t"// sse2_rnd already in xmm4 */\
			"movl	%[__cyA]	,%%edi		\n\t	movl	%[__cyB]	,%%edx		\n\t"/* cy_in */\
			"mulpd	0x1e0(%%eax),%%xmm1		\n\t	mulpd	0x1f0(%%eax),%%xmm5		\n\t"/* y = y*wi_im */\
			"addpd		(%%edi)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"/* y = y*wi_im + cy */\
			"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"/* temp = y */\
			"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
			"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"/* temp = DNINT(y) */\
			"movl	%[__sign_mask],%%ebx	\n\t"\
			"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"/* y - temp */\
			"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"/* frac = fabs(y-temp) */\
			"maxpd			%%xmm5	,%%xmm1	\n\t	maxpd		-0x20(%%eax),%%xmm1	\n\t"/* if(frac > maxerr) maxerr=frac */\
			"movaps		%%xmm1,-0x20(%%eax)	\n\t"/* Note serialization here! */\
			"movl	%%esi,	%%ebx													\n\t"/* byte offsets for bits <0:1> */\
			"shll	 $4,	%%ebx			\n\t	shll	 $2,	%%esi			\n\t"/* and <2:3> go into rbx, rsi */\
			"andl	$0x30,	%%ebx			\n\t	andl	$0x30,	%%esi			\n\t"/* respectively. */\
			"addl	%%eax,	%%ebx			\n\t	addl	%%eax,	%%esi			\n\t"\
			"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"/* cpy temp */\
			"mulpd	 	0xc0(%%ebx),%%xmm3	\n\t	mulpd	 	0xc0(%%esi),%%xmm7	\n\t"/* temp*baseinv[i] */\
			"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
			"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"/* cy_out */\
			"movaps		%%xmm3	,(%%edi)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"/* store cy_out */\
			"movl	%[__data]	,%%eax		\n\t"\
			"mulpd	 0x80(%%ebx),	%%xmm3	\n\t	mulpd	 0x80(%%esi),	%%xmm7	\n\t"/* cy*base[i] */\
			"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"/* y = (temp-cy*base[i]) */\
			"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"/* y*= wt_im */\
			"movaps		%%xmm1	,0x30(%%eax)\n\t	movaps		%%xmm5	,0x70(%%eax)\n\t"/* store y */\
		/* Update and store weights ... the (j+2) wts are in the 8 slots (half_arr+[0x180-0x1f0]) above the (j) ones: */\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movaps	0x1e0(%%eax),%%xmm3		\n\t	movaps	 0x1f0(%%eax),%%xmm7	\n\t"/* wi_im, inverse-wt_im */\
		"movaps	 0x70(%%eax),%%xmm1		\n\t	movaps		%%xmm1,%%xmm5		\n\t"/* inv_mult[1], 2 copies */\
		/* Do compare as (inv_mult[1] < wt_im), result overwrites copy of inv_mult[1]: */\
		"cmppd	$1,%%xmm2,%%xmm1		\n\t	cmppd	$1,%%xmm6,%%xmm5		\n\t"/* i = (wt_im >= inv_mult[1]) */\
		"movmskpd	%%xmm1,%%edi		\n\t	movmskpd	%%xmm5,%%edx		\n\t"/* Extract cmp-results into pair of 2-bit signmasks */\
		"shll			$4,%%edi		\n\t	shll			$4,%%edx		\n\t"/* ...and mpy by xmm bytewidth. */\
		"leal	(%%eax,%%edi),%%edi		\n\t	leal	(%%eax,%%edx),%%edx		\n\t"/* address = half_arr + i */\
		"mulpd	     (%%edi),%%xmm2		\n\t	mulpd	     (%%edx),%%xmm6		\n\t"/* wt_im *= wts_mult[i] */\
		"mulpd	 0x40(%%edi),%%xmm3		\n\t	mulpd	 0x40(%%edx),%%xmm7		\n\t"/* wi_im *= inv_mult[i] */\
		"movaps	%%xmm2,0x1c0(%%eax)		\n\t	movaps	 %%xmm6,0x1d0(%%eax)	\n\t"/* Store wt_im */\
		"movaps	%%xmm3,0x1e0(%%eax)		\n\t	movaps	 %%xmm7,0x1f0(%%eax)	\n\t"/* Store wi_im */\
			/* Get ready for next set [RE1~] : */\
			"movl	%[__sse_bw]	,%%eax		\n\t"\
			"movl	%[__sse_nm1],%%ebx		\n\t"\
			"paddd		(%%eax)	,%%xmm0		\n\t"/* bjmod[0:3] += bw  */\
			"pand		(%%ebx)	,%%xmm0		\n\t"/* bjmod[0:3] &= nm1 */\
			"movl	%[__bjmod_0],%%edi		\n\t"\
			"movaps	%%xmm0,(%%edi)			\n\t"/* Write bjmod[0:3] */\
		/**********************************************/\
		/*              Repack the data:              */\
		/**********************************************/\
			"movl	%[__data],%%eax			\n\t"\
			"movaps		0x10(%%eax)	,%%xmm1	\n\t	movaps		0x50(%%eax)	,%%xmm5	\n\t"/* reload a[jp+p0 ] */\
			"movaps		0x30(%%eax)	,%%xmm0	\n\t	movaps		0x70(%%eax)	,%%xmm4	\n\t"\
			"movaps		%%xmm1		,%%xmm3	\n\t	movaps		%%xmm5		,%%xmm7	\n\t"/* cpy a[jp    ] */\
			"unpcklpd		%%xmm0	,%%xmm1	\n\t	unpcklpd		%%xmm4	,%%xmm5	\n\t"\
			"unpckhpd		%%xmm0	,%%xmm3	\n\t	unpckhpd		%%xmm4	,%%xmm7	\n\t"\
			"movaps		%%xmm1,0x10(%%eax)	\n\t	movaps		%%xmm5,0x50(%%eax)	\n\t"/* store a[jp+p0 ] */\
			"movaps		%%xmm3,0x30(%%eax)	\n\t	movaps		%%xmm7,0x70(%%eax)	\n\t"/* Store hi imag in aj2 */\
			"movaps		    (%%eax)	,%%xmm0	\n\t	movaps		0x40(%%eax)	,%%xmm4	\n\t"/* reload a[jt+p0 ] */\
			"movaps		0x20(%%eax)	,%%xmm1	\n\t	movaps		0x60(%%eax)	,%%xmm5	\n\t"\
			"movaps		%%xmm0		,%%xmm2	\n\t	movaps		%%xmm4		,%%xmm6	\n\t"/* cpy a[jt    ] */\
			"unpckhpd		%%xmm1	,%%xmm2	\n\t	unpckhpd		%%xmm5	,%%xmm6	\n\t"\
			"unpcklpd		%%xmm1	,%%xmm0	\n\t	unpcklpd		%%xmm5	,%%xmm4	\n\t"\
			"movaps		%%xmm2,0x20(%%eax)	\n\t	movaps		%%xmm6,0x60(%%eax)	\n\t"/* Store hi real in aj2 */\
			"movaps		%%xmm0,    (%%eax)	\n\t	movaps		%%xmm4,0x40(%%eax)	\n\t"/* store a[jt+p0 ] */\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__cyA]		"m" (XcyA)		\
		, [__cyB]		"m" (XcyB)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__i]			"m" (Xi)			\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__sign_mask]	"m" (Xsign_mask)	\
		, [__sse_bw]		"m" (Xsse_bw)		\
		, [__sse_nm1]		"m" (Xsse_nm1)		\
		, [__sse_sw]		"m" (Xsse_sw)		\
		/* Prefetch: base address and 3 index offsets */\
		,	[__add0] "m" (Xadd0)\
		,	[__p1] "m" (Xp1)\
		,	[__p2] "m" (Xp2)\
		,	[__p3] "m" (Xp3)\
		: "cc","memory","eax","ebx","ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	/***************************************************************************************************************************************************/
	/********* Non-power-of-2-FFT versions of SSE2_cmplx_carry_norm_pow2_errcheck0_2B,1_2B,2_2B (only give sans-error-check version of latter 2: *******/
	/***************************************************************************************************************************************************/

	#define SSE2_cmplx_carry_norm_errcheck1_2B(Xdata,XwtA,XwtB,XwtC,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xi,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw, Xadd0,Xp1)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"movl		%[__add0],%%ecx		\n\t"/* base address for prefetch-from-main-data-array */\
	"prefetcht0	(%%ecx)	\n\t"\
	/***************Unpack the data:*************************/\
		"movl	%[__data]	,%%eax	\n\t"\
		"movaps		    (%%eax)	,%%xmm1	\n\t	movaps		0x40(%%eax)	,%%xmm5	\n\t"\
		"movaps		    (%%eax)	,%%xmm2	\n\t	movaps		0x40(%%eax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%eax)	,%%xmm1	\n\t	unpcklpd	0x60(%%eax)	,%%xmm5	\n\t"\
		"unpckhpd	0x20(%%eax)	,%%xmm2	\n\t	unpckhpd	0x60(%%eax)	,%%xmm6	\n\t"\
		"movaps		%%xmm2, 0x20(%%eax)	\n\t	movaps		%%xmm6, 0x60(%%eax)	\n\t"\
		"\n\t"\
		"movaps		0x10(%%eax)	,%%xmm2	\n\t	movaps		0x50(%%eax)	,%%xmm6	\n\t"\
		"movaps		0x10(%%eax)	,%%xmm3	\n\t	movaps		0x50(%%eax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%eax)	,%%xmm2	\n\t	unpcklpd	0x70(%%eax)	,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%eax)	,%%xmm3	\n\t	unpckhpd	0x70(%%eax)	,%%xmm7	\n\t"\
		"movaps		%%xmm2, 0x10(%%eax)	\n\t	movaps		%%xmm6, 0x50(%%eax)	\n\t"\
		"movaps		%%xmm3, 0x30(%%eax)	\n\t	movaps		%%xmm7, 0x70(%%eax)	\n\t"\
	/**********************************************/\
	/*          Real      parts                   */\
	/**********************************************/\
		"movl	%[__bjmod_0],	%%eax	\n\t"\
		"movaps		(%%eax)	,	%%xmm0	\n\t"\
		"movl	%[__sse_sw]	,	%%ebx	\n\t"\
		"movaps		(%%ebx)	,	%%xmm7	\n\t"\
		"psubd		%%xmm0	,	%%xmm7	\n\t"\
		"movmskps	%%xmm7	,	%%esi	\n\t"\
	"movl	%[__i]	,%%ecx			\n\t"/* I == 1 if it's the 0-word, in which case we force-bigword-ness by XORing esi (whose */\
	"xorl	%%ecx	,%%esi			\n\t"/* low bit will == 0 on input in this case) with I. Otherwise I == 0, thus XOR = no-op. */\
		"shll	$24		,%%esi			\n\t"\
		"movaps		%%xmm0	,%%xmm7		\n\t"\
		"movl	%[__n_minus_sil],%%ecx	\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16		,%%ecx			\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwt]	,%%edx		\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm7		\n\t"\
		"movmskps	%%xmm7	,%%edx		\n\t"\
		"shll	$8		,%%edx			\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"movaps		-0x10(%%eax),%%xmm4	\n\t"\
		"movaps		     (%%ebx),%%xmm2	\n\t	movaps		 0x10(%%ebx),%%xmm6	\n\t"\
		"movhpd		     (%%ecx),%%xmm3	\n\t	movhpd		-0x10(%%ecx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%ecx),%%xmm3	\n\t	movlpd		-0x08(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ecx	,%[__wtB]		\n\t"\
		"\n\t"\
		"shld	$20,	%%esi	,%%edi	\n\t	shld	$18,	%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%eax),%%xmm2	\n\t	mulpd	 	0x100(%%eax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%eax),%%xmm3	\n\t	mulpd	 	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	 	     (%%edi),%%xmm2	\n\t	mulpd	 	     (%%ebx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%edx),%%xmm3	\n\t	mulpd	 	0x040(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%ecx)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%eax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"movaps		%%xmm3	,(%%ecx)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,    (%%eax)\n\t	movaps		%%xmm5	,0x40(%%eax)\n\t"\
		"\n\t"\
		"movl	%[__sse_n]	,%%ebx		\n\t"\
		"movaps		(%%ebx)	,%%xmm2		\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%ebx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	/**********************************************/\
	/*          Imaginary parts                   */\
	/**********************************************/\
	"movl		%[__add0],%%ecx		\n\t"/* base address for prefetch-from-main-data-array */\
	"movl		%[__p1],%%esi	\n\t"	\
	"prefetcht0	(%%ecx,%%esi,8)	\n\t"\
		"movl	%[__sse_sw]	,%%edx		\n\t"\
		"movaps	(%%edx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%esi		\n\t"\
		"shll	$24	,%%esi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movl	%[__n_minus_silp1],%%ecx\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16	,%%ecx				\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwtm1]	,%%edx	\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8	,%%edx				\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"movaps	 0x10(%%eax)	,%%xmm1	\n\t"\
		"movaps	 0x50(%%eax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtC]	,%%ecx		\n\t"\
		"\n\t"\
		"movaps	     (%%ebx)	,%%xmm2	\n\t	movaps	 0x10(%%ebx)	,%%xmm6	\n\t"\
		"movhpd	     (%%ecx)	,%%xmm3	\n\t	movhpd	-0x10(%%ecx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%ecx)	,%%xmm3	\n\t	movlpd	-0x08(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addl	$0x20	,%%ebx			\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ebx	,%[__wtA]		\n\t"\
		"movl	%%ecx	,%[__wtC]		\n\t"\
		"shld	$20		,%%esi	,%%edi	\n\t	shld	$18		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%eax)	,%%xmm2	\n\t	mulpd	 0x120(%%eax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%eax)	,%%xmm3	\n\t	mulpd	 0x130(%%eax)	,%%xmm7	\n\t"\
		"mulpd	      (%%edi)	,%%xmm2	\n\t	mulpd	      (%%ebx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%edx)	,%%xmm3	\n\t	mulpd	 0x040(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%ecx)	,%%xmm1			\n\t	addpd	(%%edx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"andpd	(%%ebx)	,%%xmm1			\n\t	andpd	(%%ebx)	,%%xmm5			\n\t"\
		"maxpd	%%xmm5	,%%xmm1			\n\t"\
		"maxpd	-0x20(%%eax)	,%%xmm1	\n\t"\
		"movaps	%%xmm1	,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"movaps	%%xmm3	,(%%ecx)		\n\t	movaps	%%xmm7	,(%%edx)		\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x10(%%eax)	\n\t	movaps	%%xmm5	, 0x50(%%eax)	\n\t"\
		"\n\t"\
		"movl	%[__sse_n]	,%%ebx		\n\t"\
		"movaps		(%%ebx)	,%%xmm2		\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%ebx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
		"movl	%[__bjmod_0],%%ecx		\n\t"\
		"movaps	%%xmm0,(%%ecx)			\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__wtC]		"m" (XwtC)		\
		, [__cyA]		"m" (XcyA)		\
		, [__cyB]		"m" (XcyB)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__i]			"m" (Xi)			\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sign_mask]	"m" (Xsign_mask)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]		"m" (Xsinwtm1)		\
		, [__sse_bw]		"m" (Xsse_bw)		\
		, [__sse_n]		"m" (Xsse_n)		\
		, [__sse_sw]		"m" (Xsse_sw)		\
	/* Prefetch: base address and 1 index offsets */\
	,	[__add0] "m" (Xadd0)\
	,	[__p1] "m" (Xp1)\
		: "cc","memory","eax","ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

	#define SSE2_cmplx_carry_norm_errcheck2_2B(Xdata,XwtA,XwtB,XcyA,XcyB,Xbjmod_0,Xhalf_arr,Xn_minus_silp1,Xn_minus_sil,Xsign_mask,Xsinwt,Xsinwtm1,Xsse_bw,Xsse_n,Xsse_sw, Xadd0,Xp2,Xp3)\
	{\
	__asm__ volatile (\
	"pushl %%ebx	\n\t"/* Explicit save/restore of PIC register */\
	"movl		%[__add0],%%ecx		\n\t"/* base address for prefetch-from-main-data-array */\
	"movl		%[__p2],%%esi	\n\t"	\
	"prefetcht0	(%%ecx,%%esi,8)	\n\t"\
	"/**********************************************/\n\t"\
	"/*          Real      parts                   */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__bjmod_0],	%%eax	\n\t"\
		"movaps		(%%eax)	,	%%xmm0	\n\t"\
		"movl	%[__sse_sw]	,	%%ebx	\n\t"\
		"movaps		(%%ebx)	,	%%xmm1	\n\t"\
		"psubd		%%xmm0	,	%%xmm1	\n\t"\
		"movmskps	%%xmm1	,	%%esi	\n\t"\
		"\n\t"\
		"shll	$24		,%%esi			\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"movl	%[__n_minus_sil],%%ecx	\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0,	%%xmm2	,%%xmm2		\n\t"\
		"psubd		%%xmm0	,%%xmm2		\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16		,%%ecx			\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwt]	,%%edx		\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0,	%%xmm3	,%%xmm3		\n\t"\
		"psubd		%%xmm3	,%%xmm1		\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8		,%%edx			\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"\n\t"\
		"movl	%[__data],%%eax			\n\t"\
		"movaps	 0x20(%%eax),%%xmm1		\n\t	movaps		 0x60(%%eax),%%xmm5\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"movaps		-0x10(%%eax),%%xmm4	\n\t"\
		"movaps		     (%%ebx),%%xmm2	\n\t	movaps		 0x10(%%ebx),%%xmm6	\n\t"\
		"movhpd		     (%%ecx),%%xmm3	\n\t	movhpd		-0x10(%%ecx),%%xmm7	\n\t"\
		"movlpd		 0x08(%%ecx),%%xmm3	\n\t	movlpd		-0x08(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"shld	$20,	%%esi	,%%edi	\n\t	shld	$18,	%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 	0x100(%%eax),%%xmm2	\n\t	mulpd	 	0x100(%%eax),%%xmm6	\n\t"\
		"mulpd	 	0x110(%%eax),%%xmm3	\n\t	mulpd	 	0x110(%%eax),%%xmm7	\n\t"\
		"mulpd	 	     (%%edi),%%xmm2	\n\t	mulpd	 	     (%%ebx),%%xmm6	\n\t"\
		"mulpd	 	0x040(%%edx),%%xmm3	\n\t	mulpd	 	0x040(%%ecx),%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd		%%xmm3	,%%xmm1		\n\t	mulpd		%%xmm7	,%%xmm5		\n\t"\
		"addpd		(%%ecx)	,%%xmm1		\n\t	addpd		(%%edx)	,%%xmm5		\n\t"\
		"movaps		%%xmm1	,%%xmm3		\n\t	movaps		%%xmm5	,%%xmm7		\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%eax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps		%%xmm3	,%%xmm1		\n\t	movaps		%%xmm7	,%%xmm5		\n\t"\
		"mulpd	 	0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 	0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd		%%xmm4	,%%xmm3		\n\t	addpd		%%xmm4	,%%xmm7		\n\t"\
		"subpd		%%xmm4	,%%xmm3		\n\t	subpd		%%xmm4	,%%xmm7		\n\t"\
		"movaps		%%xmm3	,(%%ecx)	\n\t	movaps		%%xmm7	,(%%edx)	\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd		%%xmm3	,	%%xmm1	\n\t	subpd		%%xmm7	,	%%xmm5	\n\t"\
		"mulpd		%%xmm2	,	%%xmm1	\n\t	mulpd		%%xmm6	,	%%xmm5	\n\t"\
		"movaps		%%xmm1	,0x20(%%eax)\n\t	movaps		%%xmm5	,0x60(%%eax)\n\t"\
		"\n\t"\
		"movl	%[__sse_n]	,%%ebx		\n\t"\
		"movaps		(%%ebx)	,%%xmm2		\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%ebx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
	"/**********************************************/\n\t"\
	"/*          Imaginary parts               */\n\t"\
	"/**********************************************/\n\t"\
	"movl		%[__add0],%%ecx		\n\t"/* base address for prefetch-from-main-data-array */\
	"movl		%[__p3],%%esi	\n\t"	\
	"prefetcht0	(%%ecx,%%esi,8)	\n\t"\
		"movl	%[__sse_sw]	,%%edx		\n\t"\
		"movaps	(%%edx)	,%%xmm1			\n\t"\
		"psubd	%%xmm0	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%esi		\n\t"\
		"shll	$24	,%%esi				\n\t"\
		"movaps	%%xmm0	,%%xmm1			\n\t"\
		"movl	%[__n_minus_silp1],%%ecx\n\t"\
		"movd	%%ecx	,%%xmm2			\n\t"\
		"pshufd	$0	,%%xmm2	,%%xmm2		\n\t"\
		"psubd	%%xmm0	,%%xmm2			\n\t"\
		"movmskps	%%xmm2	,%%ecx		\n\t"\
		"shll	$16	,%%ecx				\n\t"\
		"addl	%%ecx	,%%esi			\n\t"\
		"movl	%[__sinwtm1]	,%%edx	\n\t"\
		"movd	%%edx	,%%xmm3			\n\t"\
		"pshufd	$0	,%%xmm3	,%%xmm3		\n\t"\
		"psubd	%%xmm3	,%%xmm1			\n\t"\
		"movmskps	%%xmm1	,%%edx		\n\t"\
		"shll	$8	,%%edx				\n\t"\
		"addl	%%edx	,%%esi			\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"movaps	 0x30(%%eax)	,%%xmm1	\n\t"\
		"movaps	 0x70(%%eax)	,%%xmm5	\n\t"\
		"\n\t"\
		"movl	%[__half_arr]	,%%eax	\n\t"\
		"movl	%[__wtA]	,%%ebx		\n\t"\
		"movl	%[__wtB]	,%%ecx		\n\t"\
		"\n\t"\
		"movaps	     (%%ebx)	,%%xmm2	\n\t	movaps	 0x10(%%ebx)	,%%xmm6	\n\t"\
		"movhpd	     (%%ecx)	,%%xmm3	\n\t	movhpd	-0x10(%%ecx)	,%%xmm7	\n\t"\
		"movlpd	 0x08(%%ecx)	,%%xmm3	\n\t	movlpd	-0x08(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"addl	$0x20	,%%ebx			\n\t"\
		"subl	$0x20	,%%ecx			\n\t"\
		"movl	%%ebx	,%[__wtA]		\n\t"\
		"movl	%%ecx	,%[__wtB]		\n\t"\
		"shld	$20		,%%esi	,%%edi	\n\t	shld	$18		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"shld	$28,	%%esi,	%%edx	\n\t	shld	$26,	%%esi,	%%ecx	\n\t"\
		"andl	$0x00000030	,%%edx		\n\t	andl	$0x00000030	,%%ecx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"addl	%%eax	,%%edx			\n\t	addl	%%eax	,%%ecx			\n\t"\
		"\n\t"\
		"mulpd	 0x120(%%eax)	,%%xmm2	\n\t	mulpd	 0x120(%%eax)	,%%xmm6	\n\t"\
		"mulpd	 0x130(%%eax)	,%%xmm3	\n\t	mulpd	 0x130(%%eax)	,%%xmm7	\n\t"\
		"mulpd	      (%%edi)	,%%xmm2	\n\t	mulpd	      (%%ebx)	,%%xmm6	\n\t"\
		"mulpd	 0x040(%%edx)	,%%xmm3	\n\t	mulpd	 0x040(%%ecx)	,%%xmm7	\n\t"\
		"\n\t"\
		"movl	%[__cyA]	,%%ecx		\n\t	movl	%[__cyB]	,%%edx		\n\t"\
		"mulpd	%%xmm3	,%%xmm1			\n\t	mulpd	%%xmm7	,%%xmm5			\n\t"\
		"addpd	(%%ecx)	,%%xmm1			\n\t	addpd	(%%edx)	,%%xmm5			\n\t"\
		"movaps	%%xmm1	,%%xmm3			\n\t	movaps	%%xmm5	,%%xmm7			\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"\n\t"\
		"movl	%[__sign_mask],%%ebx	\n\t"\
		"subpd			%%xmm3	,%%xmm1	\n\t	subpd			%%xmm7	,%%xmm5	\n\t"\
		"andpd			(%%ebx)	,%%xmm1	\n\t	andpd			(%%ebx)	,%%xmm5	\n\t"\
		"maxpd			%%xmm5	,%%xmm1	\n\t"\
		"maxpd		-0x20(%%eax),%%xmm1	\n\t"\
		"movaps		%%xmm1,-0x20(%%eax)	\n\t"\
		"\n\t"\
		"shld	$12		,%%esi	,%%edi	\n\t	shld	$10		,%%esi	,%%ebx	\n\t"\
		"andl	$0x00000030	,%%edi		\n\t	andl	$0x00000030	,%%ebx		\n\t"\
		"addl	%%eax	,%%edi			\n\t	addl	%%eax	,%%ebx			\n\t"\
		"movaps	%%xmm3	,%%xmm1			\n\t	movaps	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	 0xc0(%%edi)	,%%xmm3	\n\t	mulpd	 0xc0(%%ebx)	,%%xmm7	\n\t"\
		"addpd	%%xmm4	,%%xmm3			\n\t	addpd	%%xmm4	,%%xmm7			\n\t"\
		"subpd	%%xmm4	,%%xmm3			\n\t	subpd	%%xmm4	,%%xmm7			\n\t"\
		"movaps	%%xmm3	,(%%ecx)		\n\t	movaps	%%xmm7	,(%%edx)		\n\t"\
		"\n\t"\
		"movl	%[__data]	,%%eax		\n\t"\
		"mulpd	 0x80(%%edi)	,%%xmm3	\n\t	mulpd	 0x80(%%ebx)	,%%xmm7	\n\t"\
		"subpd	%%xmm3	,%%xmm1			\n\t	subpd	%%xmm7	,%%xmm5			\n\t"\
		"mulpd	%%xmm2	,%%xmm1			\n\t	mulpd	%%xmm6	,%%xmm5			\n\t"\
		"movaps	%%xmm1	, 0x30(%%eax)	\n\t	movaps	%%xmm5	, 0x70(%%eax)	\n\t"\
		"\n\t"\
		"movl	%[__sse_n]	,%%ebx		\n\t"\
		"movaps		(%%ebx)	,%%xmm2		\n\t"\
		"movl	%[__sse_bw]	,%%eax		\n\t"\
		"paddd		(%%eax)	,%%xmm0		\n\t"\
		"movaps		%%xmm0	,%%xmm1		\n\t"\
		"pcmpgtd	%%xmm2	,%%xmm1		\n\t"\
		"pand		(%%ebx)	,%%xmm1		\n\t"\
		"psubd		%%xmm1	,%%xmm0		\n\t"\
		"movl	%[__bjmod_0],%%ecx		\n\t"\
		"movaps	%%xmm0,(%%ecx)			\n\t"\
	"/**********************************************/\n\t"\
	"/*          Repack the data:              */\n\t"\
	"/**********************************************/\n\t"\
		"movl	%[__data],%%eax			\n\t"\
		"movaps		0x10(%%eax)	,%%xmm1	\n\t	movaps		0x50(%%eax)	,%%xmm5	\n\t"\
		"movaps		    (%%eax)	,%%xmm0	\n\t	movaps		0x40(%%eax)	,%%xmm4	\n\t"\
		"movaps		%%xmm1		,%%xmm3	\n\t	movaps		%%xmm5		,%%xmm7	\n\t"\
		"movaps		%%xmm0		,%%xmm2	\n\t	movaps		%%xmm4		,%%xmm6	\n\t"\
		"unpckhpd	0x30(%%eax)	,%%xmm3	\n\t	unpckhpd	0x70(%%eax)	,%%xmm7	\n\t"\
		"unpcklpd	0x30(%%eax)	,%%xmm1	\n\t	unpcklpd	0x70(%%eax)	,%%xmm5	\n\t"\
		"movaps		%%xmm3,0x30(%%eax)	\n\t	movaps		%%xmm7,0x70(%%eax)	\n\t"\
		"unpckhpd	0x20(%%eax)	,%%xmm2	\n\t	unpckhpd	0x60(%%eax)	,%%xmm6	\n\t"\
		"unpcklpd	0x20(%%eax)	,%%xmm0	\n\t	unpcklpd	0x60(%%eax)	,%%xmm4	\n\t"\
		"movaps		%%xmm2,0x20(%%eax)	\n\t	movaps		%%xmm6,0x60(%%eax)	\n\t"\
		"movaps		%%xmm1,0x10(%%eax)	\n\t	movaps		%%xmm5,0x50(%%eax)	\n\t"\
		"movaps		%%xmm0,    (%%eax)	\n\t	movaps		%%xmm4,0x40(%%eax)	\n\t"\
	"popl %%ebx	\n\t"\
		:					/* outputs: none */\
		: [__data]		"m" (Xdata)	/* All inputs from memory addresses here */\
		, [__wtA]		"m" (XwtA)		\
		, [__wtB]		"m" (XwtB)		\
		, [__cyA]		"m" (XcyA)		\
		, [__cyB]		"m" (XcyB)		\
		, [__bjmod_0]	"m" (Xbjmod_0)		\
		, [__half_arr]	"m" (Xhalf_arr)		\
		, [__n_minus_silp1] "m" (Xn_minus_silp1)\
		, [__n_minus_sil]	"m" (Xn_minus_sil)	\
		, [__sign_mask]	"m" (Xsign_mask)	\
		, [__sinwt]		"m" (Xsinwt)		\
		, [__sinwtm1]		"m" (Xsinwtm1)		\
		, [__sse_bw]		"m" (Xsse_bw)		\
		, [__sse_n]		"m" (Xsse_n)		\
		, [__sse_sw]		"m" (Xsse_sw)		\
	/* Prefetch: base address and 2 index offsets */\
	,	[__add0] "m" (Xadd0)\
	,	[__p2] "m" (Xp2)\
	,	[__p3] "m" (Xp3)\
		: "cc","memory","eax","ecx","edx","edi","esi","xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"	/* Clobbered registers */\
	);\
	}

#endif	/* carry_gcc_h_included */

