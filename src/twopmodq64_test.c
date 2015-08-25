/*******************************************************************************
*                                                                              *
*   (C) 1997-2012 by Ernst W. Mayer.                                           *
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

#include "factor.h"	/* Need this for include of platform-specific defs used in YES_ASM test below */

#undef YES_ASM
#if(defined(CPU_IS_X86_64) && defined(COMPILER_TYPE_GCC) && (OS_BITS == 64))
	#define YES_ASM
#endif

#ifdef YES_ASM
char str0[64];
/*** 4-trial-factor version, optimized x86_64 asm for main powering loop ***/
uint64 twopmodq64_q4(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3)
{
	int32 j;
	uint64 q0 = 1+(p<<1)*k0, q1 = 1+(p<<1)*k1, q2 = 1+(p<<1)*k2, q3 = 1+(p<<1)*k3
	, qinv0, qinv1, qinv2, qinv3
	, x0, x1, x2, x3
	, y0, y1, y2, y3
	, lo0, lo1, lo2, lo3, r;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;
	
	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 64;
		start_index = 64-leadz64(pshift)-6;
		
		zshift = 63-ibits64(pshift,start_index,6);	/* Since 6 bits with leftmost bit = 1 is guaranteed
		 to be in [32,63], the shift count here is in [0, 31].
		 That means that zstart < 2^32. Together with the fact that
		 squaring a power of two gives another power of two, we can
		 simplify the modmul code sequence for the first iteration.
		 Every little bit counts (literally in this case :), right? */
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */
		
		pshift = ~pshift;
	}
	*checksum1 += q0 + q1 + q2 + q3;

	/* q must be odd for Montgomery-style modmul to work: */
#if FAC_DEBUG
	ASSERT(HERE, (q0 & (uint64)1) == 1, "(q0 & (uint64)1) != 1!");
	ASSERT(HERE, (q1 & (uint64)1) == 1, "(q1 & (uint64)1) != 1!");
	ASSERT(HERE, (q2 & (uint64)1) == 1, "(q2 & (uint64)1) != 1!");
	ASSERT(HERE, (q3 & (uint64)1) == 1, "(q3 & (uint64)1) != 1!");
#endif
	qinv0 = (q0+q0+q0) ^ (uint64)2;
	qinv1 = (q1+q1+q1) ^ (uint64)2;
	qinv2 = (q2+q2+q2) ^ (uint64)2;
	qinv3 = (q3+q3+q3) ^ (uint64)2;
	
	for(j = 0; j < 4; j++)
	{
		qinv0 = qinv0*((uint64)2 - q0*qinv0);
		qinv1 = qinv1*((uint64)2 - q1*qinv1);
		qinv2 = qinv2*((uint64)2 - q2*qinv2);
		qinv3 = qinv3*((uint64)2 - q3*qinv3);
	}
	
	/* Since zstart is a power of two < 2^128, use a
	 streamlined code sequence for the first iteration: */
	j = start_index-1;
	
	/* MULL64(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	lo0 = qinv0 << zshift;
	lo1 = qinv1 << zshift;
	lo2 = qinv2 << zshift;
	lo3 = qinv3 << zshift;
	MULH64(q0,lo0,lo0);
	MULH64(q1,lo1,lo1);
	MULH64(q2,lo2,lo2);
	MULH64(q3,lo3,lo3);

	/* hi = 0 in this instance, which simplifies things. */
	y0 = q0 - lo0;
	y1 = q1 - lo1;
	y2 = q2 - lo2;
	y3 = q3 - lo3;
	
	if((pshift >> j) & (uint64)1)
	{
		x0 = y0 + y0;	x0 -= q0 & -(x0 >= q0 || x0 < y0);
		x1 = y1 + y1;	x1 -= q1 & -(x1 >= q1 || x1 < y1);
		x2 = y2 + y2;	x2 -= q2 & -(x2 >= q2 || x2 < y2);
		x3 = y3 + y3;	x3 -= q3 & -(x3 >= q3 || x3 < y3);
	}
	else
	{
		x0 = y0;
		x1 = y1;
		x2 = y2;
		x3 = y3;
	}

//printf("twopmodq64_q4 : x1 = %s\n", &str0[convert_uint64_base10_char(str0, x1)] );
//for(j = start_index-2; j >= 0; j--) {
	__asm__ volatile (\
	"/* Load the q0|1|2 into rbx|rsi|rdi, keeping rcx free for loop index and rax:rdx for double-width MULs: */	\n\t"\
		"movq	%[__q0],%%rbx	\n\t"\
		"movq	%[__q1],%%rsi	\n\t"\
		"movq	%[__q2],%%rdi	\n\t"\
		"/* Must load q3 as-needed due to lack of user register to store it */\n\t"\
	"/* Load the x's into r12-15: */	\n\t"\
		"movq	%[__x0],%%r12	\n\t"\
		"movq	%[__x1],%%r13	\n\t"\
		"movq	%[__x2],%%r14	\n\t"\
		"movq	%[__x3],%%r15	\n\t"\
	"/* Pure-ASM loop control: for(j = start_index-2; j >= 0; j--) */\n\t"\
		"movslq	%[__start_index], %%rcx		\n\t"\
		"subq $2,%%rcx						\n\t"\
		"test %%rcx, %%rcx					\n\t"\
		"jl LoopEnd4		/* Skip if n < 0 */	\n\t"\
	"LoopBeg4:								\n\t"\
	"/* SQR_LOHI_q4(x*, lo*, hi*): */	\n\t"\
		"movq	%%r12,%%rax	/* x0-3 in r8-11. */\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r8 	/* lo0 */	\n\t"\
		"movq	%%rdx,%%r12	/* hi0 */	\n\t"\
		"movq	%%r13,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r9 	/* lo1 */	\n\t"\
		"movq	%%rdx,%%r13	/* hi1 */	\n\t"\
		"movq	%%r14,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r10	/* lo2 */	\n\t"\
		"movq	%%rdx,%%r14	/* hi2 */	\n\t"\
		"movq	%%r15,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r11	/* lo3 */	\n\t"\
		"movq	%%rdx,%%r15	/* hi3 */	\n\t"\
		"\n\t"\
		"\n\t"\
	"/* MULL_q4((qinv*, lo*, lo*): */	\n\t"\
		"\n\t"\
		"imulq	%[__qinv0],%%r8 	\n\t"\
		"imulq	%[__qinv1],%%r9 	\n\t"\
		"imulq	%[__qinv2],%%r10	\n\t"\
		"imulq	%[__qinv3],%%r11	\n\t"\
		"\n\t"\
	"/* UMULH_q4((q*, lo*, lo*): lo0-3 in r8-11: */	\n\t"\
		"\n\t"\
		"movq	%%r8 ,%%rax	/* lo0 */\n\t"\
		"mulq	%%rbx	/* q0 in rbx; q0*lo0 in rax:rdx */\n\t"\
		"movq	%%rdx,%%r8 	/* Discard low 64 bits [rax] */\n\t"\
		"\n\t"\
		"movq	%%r9 ,%%rax	/* lo1 */\n\t"\
		"mulq	%%rsi	/* q1 in rsi; q1*lo1 in rax:rdx */\n\t"\
		"movq	%%rdx,%%r9 	/* Discard low 64 bits [rax] */\n\t"\
		"\n\t"\
		"movq	%%r10,%%rax	/* lo2 */\n\t"\
		"mulq	%%rdi	/* q2 in rdi; q2*lo2 in rax:rdx */\n\t"\
		"movq	%%rdx,%%r10	/* Discard low 64 bits [rax] */\n\t"\
		"\n\t"\
		"movq	%%r11,%%rax	/* lo3 */\n\t"\
		"mulq	%[__q3]		\n\t"\
		"movq	%%rdx,%%r11	/* Discard low 64 bits [rax] */\n\t"\
		"\n\t"\
	"/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\n\t"\
		"subq	%%r8 ,%%r12			\n\t" /* r12 = (h-l) */\
		"leaq (%%r12,%%rbx),%%r8 	\n\t" /* r8  = (h-l)+q */\
		"cmovcq %%r8 ,%%r12	\n\t" /* if carry = true (i.e. h > l), copy source (r8 = h-l+q) into dest (r12), else leave dest = h-l. */\
		"\n\t"\
		"subq	%%r9 ,%%r13			\n\t"\
		"leaq (%%r13,%%rsi),%%r9 	\n\t"\
		"cmovcq %%r9 ,%%r13	\n\t"\
		"\n\t"\
		"subq	%%r10,%%r14			\n\t"\
		"leaq (%%r14,%%rdi),%%r10	\n\t"\
		"cmovcq %%r10,%%r14	\n\t"\
		"\n\t"\
		"movq	%[__q3],%%rdx	/* Re-use q3 several times below, so store in free register */\n\t"\
		"subq	%%r11,%%r15			\n\t"\
		"leaq (%%r15,%%rdx),%%r11	\n\t"\
		"cmovcq %%r11,%%r15	\n\t"\
		"\n\t"\
"/* If current bit of pshift == 1, double each output modulo q: */	\n\t"\
		"/* if((pshift >> j) & (uint64)1) { */	\n\t"\
		"movl	%[__pshift],%%eax	/* Need to follow this with load-j-into-ecx if use HLL loop control in debug mode */\n\t"\
		"shrq	%%cl,%%rax				\n\t"\
		"andq	$0x1,%%rax				\n\t"\
	"je	twopmodq64_q4_pshiftjmp			\n\t"\
		"\n\t"\
		"movq	%%r12,%%r8 	/* r8  <- Copy of x */\n\t"\
		"movq	%%r12,%%rax	/* rax <- Copy of x */\n\t"\
		"leaq (%%r12,%%r12),%%r12	/* x+x */\n\t"\
		"subq	%%rbx,%%r8	/* r8  <- x-q */\n\t"\
		"addq	%%rax,%%r8 	/* r8 <- 2x-q */\n\t"\
		"cmovcq %%r8 ,%%r12	\n\t"/* if carry (i.e. x+x needed modding), copy source (r8 = x+x-q) into dest (r12), else leave dest = x+x. */\
		"\n\t"\
		"movq	%%r13,%%r9 	\n\t"\
		"movq	%%r13,%%rax	\n\t"\
		"leaq (%%r13,%%r13),%%r13	\n\t"\
		"subq	%%rsi,%%r9	\n\t"\
		"addq	%%rax,%%r9 	\n\t"\
		"cmovcq %%r9 ,%%r13	\n\t"\
		"\n\t"\
		"movq	%%r14,%%r10	\n\t"\
		"movq	%%r14,%%rax	\n\t"\
		"leaq (%%r14,%%r14),%%r14	\n\t"\
		"subq	%%rdi,%%r10	\n\t"\
		"addq	%%rax,%%r10	\n\t"\
		"cmovcq %%r10,%%r14	\n\t"\
		"\n\t"\
		"movq	%%r15,%%r11	\n\t"\
		"movq	%%r15,%%rax	\n\t"\
		"leaq (%%r15,%%r15),%%r15	\n\t"\
		"subq	%[__q3],%%r11	/* Weird...using rdx in place of __q3 here made timings 1-2 percent  worse. */\n\t"\
		"addq	%%rax,%%r11	\n\t"\
		"cmovcq %%r11,%%r15	\n\t"\
		"\n\t"\
	"twopmodq64_q4_pshiftjmp:					\n\t"\
		"/* } endif((pshift >> j) & (uint64)1) */						\n\t"\
		"subq	$1,%%rcx	/* j-- */		\n\t"\
		"cmpq	$0,%%rcx	/* compare j vs 0 */\n\t"\
		"jge	LoopBeg4	/* if (j >= 0), Loop */	\n\t"\
	"LoopEnd4:							\n\t"\
		"movq	%%r12,%[__x0]	\n\t"\
		"movq	%%r13,%[__x1]	\n\t"\
		"movq	%%r14,%[__x2]	\n\t"\
		"movq	%%r15,%[__x3]	\n\t"\
		:	/* outputs: none */\
		: [__q0] "m" (q0)	/* All inputs from memory addresses here */\
		 ,[__q1] "m" (q1)	\
		 ,[__q2] "m" (q2)	\
		 ,[__q3] "m" (q3)	\
		 ,[__qinv0] "m" (qinv0)	\
		 ,[__qinv1] "m" (qinv1)	\
		 ,[__qinv2] "m" (qinv2)	\
		 ,[__qinv3] "m" (qinv3)	\
		 ,[__x0] "m" (x0)	\
		 ,[__x1] "m" (x1)	\
		 ,[__x2] "m" (x2)	\
		 ,[__x3] "m" (x3)	\
		 ,[__pshift] "m" (pshift)	\
		 ,[__j] "m" (j)	/* Only need this if debug and explicit loop enabled, but try with/sans for each version of the asm, pivk faster one. */\
		 ,[__start_index] "m" (start_index)	\
		: "cc","memory","rax","rbx","rcx","rdx","rsi","rdi","r8","r9","r10","r11","r12","r13","r14","r15"		/* Clobbered registers */\
		);
//printf("twopmodq64_q4 : x1 = %s\n", &str0[convert_uint64_base10_char(str0, x1)] );
//	}
//printf("twopmodq64_q4 : x1 = %s\n", &str0[convert_uint64_base10_char(str0, x1+x1-q1)] );
//exit(0);
	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2. */

	*checksum2 += x0 + x1 + x2 + x3;

	r = 0;
	if(x0+x0-q0 == 1) r +=  1;
	if(x1+x1-q1 == 1) r +=  2;
	if(x2+x2-q2 == 1) r +=  4;
	if(x3+x3-q3 == 1) r +=  8;
	return(r);
}


/*** 8-trial-factor version ***/
uint64 twopmodq64_q8(uint64*checksum1, uint64*checksum2, uint64 p, uint64 k0, uint64 k1, uint64 k2, uint64 k3, uint64 k4, uint64 k5, uint64 k6, uint64 k7)
{
	 int32 j;
	uint64 q0 = 1+(p<<1)*k0, q1 = 1+(p<<1)*k1, q2 = 1+(p<<1)*k2, q3 = 1+(p<<1)*k3, q4 = 1+(p<<1)*k4, q5 = 1+(p<<1)*k5, q6 = 1+(p<<1)*k6, q7 = 1+(p<<1)*k7
		, qinv0, qinv1, qinv2, qinv3, qinv4, qinv5, qinv6, qinv7
		, x0, x1, x2, x3, x4, x5, x6, x7
		, y0, y1, y2, y3, y4, y5, y6, y7
		, lo0, lo1, lo2, lo3, lo4, lo5, lo6, lo7, r;
	static uint64 psave = 0, pshift;
	static uint32 start_index, zshift, first_entry = TRUE;

	if(first_entry || p != psave)
	{
		first_entry = FALSE;
		psave  = p;
		pshift = p + 64;
		start_index = 64-leadz64(pshift)-6;

		zshift = 63-ibits64(pshift,start_index,6);	/* Since 6 bits with leftmost bit = 1 is guaranteed
								to be in [32,63], the shift count here is in [0, 31].
								That means that zstart < 2^32. Together with the fact that
								squaring a power of two gives another power of two, we can
								simplify the modmul code sequence for the first iteration.
								Every little bit counts (literally in this case :), right? */
		zshift <<= 1;				/* Doubling the shift count here takes cares of the first SQR_LOHI */

		pshift = ~pshift;
	}
	*checksum1 += q0 + q1 + q2 + q3 + q4 + q5 + q6 + q7;

	/* q must be odd for Montgomery-style modmul to work: */
	qinv0 = (q0+q0+q0) ^ (uint64)2;
	qinv1 = (q1+q1+q1) ^ (uint64)2;
	qinv2 = (q2+q2+q2) ^ (uint64)2;
	qinv3 = (q3+q3+q3) ^ (uint64)2;
	qinv4 = (q4+q4+q4) ^ (uint64)2;
	qinv5 = (q5+q5+q5) ^ (uint64)2;
	qinv6 = (q6+q6+q6) ^ (uint64)2;
	qinv7 = (q7+q7+q7) ^ (uint64)2;

	for(j = 0; j < 4; j++)
	{
		qinv0 = qinv0*((uint64)2 - q0*qinv0);
		qinv1 = qinv1*((uint64)2 - q1*qinv1);
		qinv2 = qinv2*((uint64)2 - q2*qinv2);
		qinv3 = qinv3*((uint64)2 - q3*qinv3);
		qinv4 = qinv4*((uint64)2 - q4*qinv4);
		qinv5 = qinv5*((uint64)2 - q5*qinv5);
		qinv6 = qinv6*((uint64)2 - q6*qinv6);
		qinv7 = qinv7*((uint64)2 - q7*qinv7);
	}

	/* Since zstart is a power of two < 2^128, use a
	streamlined code sequence for the first iteration: */
	j = start_index-1;

	/* MULL64(zstart,qinv,lo) simply amounts to a left-shift of the bits of qinv: */
	lo0 = qinv0 << zshift;
	lo1 = qinv1 << zshift;
	lo2 = qinv2 << zshift;
	lo3 = qinv3 << zshift;
	lo4 = qinv4 << zshift;
	lo5 = qinv5 << zshift;
	lo6 = qinv6 << zshift;
	lo7 = qinv7 << zshift;

	MULH64(q0,lo0,lo0);
	MULH64(q1,lo1,lo1);
	MULH64(q2,lo2,lo2);
	MULH64(q3,lo3,lo3);
	MULH64(q4,lo4,lo4);
	MULH64(q5,lo5,lo5);
	MULH64(q6,lo6,lo6);
	MULH64(q7,lo7,lo7);

	/* hi = 0 in this instance, which simplifies things. */
	y0 = q0 - lo0;
	y1 = q1 - lo1;
	y2 = q2 - lo2;
	y3 = q3 - lo3;
	y4 = q4 - lo4;
	y5 = q5 - lo5;
	y6 = q6 - lo6;
	y7 = q7 - lo7;

	if((pshift >> j) & (uint64)1)
	{
		x0 = y0 + y0;	x0 -= q0 & -(x0 >= q0 || x0 < y0);
		x1 = y1 + y1;	x1 -= q1 & -(x1 >= q1 || x1 < y1);
		x2 = y2 + y2;	x2 -= q2 & -(x2 >= q2 || x2 < y2);
		x3 = y3 + y3;	x3 -= q3 & -(x3 >= q3 || x3 < y3);
		x4 = y4 + y4;	x4 -= q4 & -(x4 >= q4 || x4 < y4);
		x5 = y5 + y5;	x5 -= q5 & -(x5 >= q5 || x5 < y5);
		x6 = y6 + y6;	x6 -= q6 & -(x6 >= q6 || x6 < y6);
		x7 = y7 + y7;	x7 -= q7 & -(x7 >= q7 || x7 < y7);
	}
	else
	{
		x0 = y0;
		x1 = y1;
		x2 = y2;
		x3 = y3;
		x4 = y4;
		x5 = y5;
		x6 = y6;
		x7 = y7;
	}

	__asm__ volatile (\
	"/* Load the x's into r8-15: */	\n\t"\
		"movq	%[__x0],%%r8 	\n\t"\
		"movq	%[__x1],%%r9 	\n\t"\
		"movq	%[__x2],%%r10	\n\t"\
		"movq	%[__x3],%%r11	\n\t"\
		"movq	%[__x4],%%r12	\n\t"\
		"movq	%[__x5],%%r13	\n\t"\
		"movq	%[__x6],%%r14	\n\t"\
		"movq	%[__x7],%%r15	\n\t"\
	"/* Pure-ASM loop control: for(j = start_index-2; j >= 0; j--) */\n\t"\
		"movslq	%[__start_index], %%rcx		\n\t"\
		"subq $2,%%rcx						\n\t"\
		"test %%rcx, %%rcx					\n\t"\
		"jl LoopEnd8		/* Skip if n < 0 */	\n\t"\
	"LoopBeg8:								\n\t"\
	"/* SQR_LOHI_q4(x*, lo*, hi*): */	\n\t"\
		"movq	%%r8 ,%%rax		\n\t"\
		"mulq	%%rax			\n\t"\
		"movq	%%rax,%%r8		/* lo */	\n\t"\
		"movq	%%rdx,%[__x0]	/* Move hi back into x */	\n\t"\
		"movq	%%r9 ,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r9  	\n\t"\
		"movq	%%rdx,%[__x1]	\n\t"\
		"movq	%%r10,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r10 	\n\t"\
		"movq	%%rdx,%[__x2]	\n\t"\
		"movq	%%r11,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r11 	\n\t"\
		"movq	%%rdx,%[__x3]	\n\t"\
		"movq	%%r12,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r12 	\n\t"\
		"movq	%%rdx,%[__x4]	\n\t"\
		"movq	%%r13,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r13 	\n\t"\
		"movq	%%rdx,%[__x5]	\n\t"\
		"movq	%%r14,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r14		\n\t"\
		"movq	%%rdx,%[__x6]	\n\t"\
		"movq	%%r15,%%rax	\n\t"\
		"mulq	%%rax		\n\t"\
		"movq	%%rax,%%r15		\n\t"\
		"movq	%%rdx,%[__x7]	\n\t"\
		"\n\t"\
		"\n\t"\
	"/* MULL_q4((qinv*, lo*, lo*): */	\n\t"\
		"\n\t"\
		"imulq	%[__qinv0],%%r8 	\n\t"\
		"imulq	%[__qinv1],%%r9 	\n\t"\
		"imulq	%[__qinv2],%%r10	\n\t"\
		"imulq	%[__qinv3],%%r11	\n\t"\
		"imulq	%[__qinv4],%%r12	\n\t"\
		"imulq	%[__qinv5],%%r13	\n\t"\
		"imulq	%[__qinv6],%%r14	\n\t"\
		"imulq	%[__qinv7],%%r15	\n\t"\
		"\n\t"\
	"/* UMULH_q4((q*, lo*, lo*): lo0-3 in r8-11: */	\n\t"\
		"\n\t"\
		"movq	%%r8 ,%%rax	/* lo0 */\n\t"\
		"mulq	%[__q0]		/* q0*lo0 in rax:rdx */\n\t"\
		"movq	%%rdx,%%r8 	/* Discard low 64 bits [rax] */\n\t"\
		"\n\t"\
		"movq	%%r9 ,%%rax	\n\t"\
		"mulq	%[__q1]		\n\t"\
		"movq	%%rdx,%%r9 	\n\t"\
		"\n\t"\
		"movq	%%r10,%%rax	\n\t"\
		"mulq	%[__q2]		\n\t"\
		"movq	%%rdx,%%r10	\n\t"\
		"\n\t"\
		"movq	%%r11,%%rax	\n\t"\
		"mulq	%[__q3]		\n\t"\
		"movq	%%rdx,%%r11	\n\t"\
		"\n\t"\
		"movq	%%r12,%%rax	\n\t"\
		"mulq	%[__q4]		\n\t"\
		"movq	%%rdx,%%r12	\n\t"\
		"\n\t"\
		"movq	%%r13,%%rax	\n\t"\
		"mulq	%[__q5]		\n\t"\
		"movq	%%rdx,%%r13	\n\t"\
		"\n\t"\
		"movq	%%r14,%%rax	\n\t"\
		"mulq	%[__q6]		\n\t"\
		"movq	%%rdx,%%r14	\n\t"\
		"\n\t"\
		"movq	%%r15,%%rax	\n\t"\
		"mulq	%[__q7]		\n\t"\
		"movq	%%rdx,%%r15	\n\t"\
		"\n\t"\
	"/* If h < l, then calculate h-l+q < q; otherwise calculate h-l. */\n\t"\
		"movq	%[__x0],%%rax	\n\t" /* h */\
		"movq	%[__q0],%%rbx	\n\t" /* h */\
		"subq	%%r8 ,%%rax		\n\t" /* rax = (h-l), sets cf used by cmov */\
		"leaq (%%rax,%%rbx),%%r8 	\n\t" /* r8  = (h-l)+q, no cf set */\
		"cmovncq %%rax,%%r8 	\n\t" /* if carry = false (i.e. h >= l), copy source (rax = h-l) into dest (r8), else leave dest = (h-l+q. */\
		"\n\t"\
		"movq	%[__x1],%%rax	\n\t"\
		"movq	%[__q1],%%rbx	\n\t"\
		"subq	%%r9 ,%%rax		\n\t"\
		"leaq (%%rax,%%rbx),%%r9 	\n\t"\
		"cmovncq %%rax,%%r9 	\n\t"\
		"\n\t"\
		"movq	%[__x2],%%rax	\n\t"\
		"movq	%[__q2],%%rbx	\n\t"\
		"subq	%%r10,%%rax		\n\t"\
		"leaq (%%rax,%%rbx),%%r10	\n\t"\
		"cmovncq %%rax,%%r10	\n\t"\
		"\n\t"\
		"movq	%[__x3],%%rax	\n\t"\
		"movq	%[__q3],%%rbx	\n\t"\
		"subq	%%r11,%%rax		\n\t"\
		"leaq (%%rax,%%rbx),%%r11	\n\t"\
		"cmovncq %%rax,%%r11	\n\t"\
		"\n\t"\
		"movq	%[__x4],%%rax	\n\t"\
		"movq	%[__q4],%%rbx	\n\t"\
		"subq	%%r12,%%rax		\n\t"\
		"leaq (%%rax,%%rbx),%%r12	\n\t"\
		"cmovncq %%rax,%%r12	\n\t"\
		"\n\t"\
		"movq	%[__x5],%%rax	\n\t"\
		"movq	%[__q5],%%rbx	\n\t"\
		"subq	%%r13,%%rax		\n\t"\
		"leaq (%%rax,%%rbx),%%r13	\n\t"\
		"cmovncq %%rax,%%r13	\n\t"\
		"\n\t"\
		"movq	%[__x6],%%rax	\n\t"\
		"movq	%[__q6],%%rbx	\n\t"\
		"subq	%%r14,%%rax		\n\t"\
		"leaq (%%rax,%%rbx),%%r14	\n\t"\
		"cmovncq %%rax,%%r14	\n\t"\
		"\n\t"\
		"movq	%[__x7],%%rax	\n\t"\
		"movq	%[__q7],%%rbx	\n\t"\
		"subq	%%r15,%%rax		\n\t"\
		"leaq (%%rax,%%rbx),%%r15	\n\t"\
		"cmovncq %%rax,%%r15	\n\t"\
		"\n\t"\
"/* If current bit of pshift == 1, double each output modulo q: */	\n\t"\
		"/* if((pshift >> j) & (uint64)1) { */	\n\t"\
		"movl	%[__pshift],%%eax	/* Need to follow this with load-j-into-ecx if use HLL loop control in debug mode */\n\t"\
		"shrq	%%cl,%%rax				\n\t"\
		"andq	$0x1,%%rax				\n\t"\
	"je	twopmodq64_q8_pshiftjmp			\n\t"\
		"\n\t"\
		"movq	%[__q0],%%rax  	\n\t"/* rax <- Copy of q */\
		"leaq (%%r8 ,%%r8 ),%%rbx	\n\t"/* rbx <- x+x, no CF set, leave r8 (x) intact */\
		"subq	%%r8 ,%%rax		\n\t"/* rax <- q-x, sets CF */\
		"subq	%%rax,%%r8 		\n\t"/* r8  <- x-(q-x) = x-q+x */\
		"cmovcq %%rbx,%%r8		\n\t" /* if borrow (i.e. x+x did not need modding), copy source (rax = x+x) into dest (r8 ), else leave dest = x+x-q. */\
		"\n\t"\
		"movq	%[__q1],%%rax  	\n\t"\
		"leaq (%%r9 ,%%r9 ),%%rbx	\n\t"\
		"subq	%%r9 ,%%rax		\n\t"\
		"subq	%%rax,%%r9 		\n\t"\
		"cmovcq %%rbx,%%r9 		\n\t"\
		"\n\t"\
		"movq	%[__q2],%%rax  	\n\t"\
		"leaq (%%r10,%%r10),%%rbx	\n\t"\
		"subq	%%r10,%%rax		\n\t"\
		"subq	%%rax,%%r10		\n\t"\
		"cmovcq %%rbx,%%r10 	\n\t"\
		"\n\t"\
		"movq	%[__q3],%%rax  	\n\t"\
		"leaq (%%r11,%%r11),%%rbx	\n\t"\
		"subq	%%r11,%%rax		\n\t"\
		"subq	%%rax,%%r11		\n\t"\
		"cmovcq %%rbx,%%r11		\n\t"\
		"\n\t"\
		"movq	%[__q4],%%rax  	\n\t"\
		"leaq (%%r12,%%r12),%%rbx	\n\t"\
		"subq	%%r12,%%rax		\n\t"\
		"subq	%%rax,%%r12		\n\t"\
		"cmovcq %%rbx,%%r12		\n\t"\
		"\n\t"\
		"movq	%[__q5],%%rax  	\n\t"\
		"leaq (%%r13,%%r13),%%rbx	\n\t"\
		"subq	%%r13,%%rax		\n\t"\
		"subq	%%rax,%%r13		\n\t"\
		"cmovcq %%rbx,%%r13		\n\t"\
		"\n\t"\
		"movq	%[__q6],%%rax  	\n\t"\
		"leaq (%%r14,%%r14),%%rbx	\n\t"\
		"subq	%%r14,%%rax		\n\t"\
		"subq	%%rax,%%r14		\n\t"\
		"cmovcq %%rbx,%%r14		\n\t"\
		"\n\t"\
		"movq	%[__q7],%%rax  	\n\t"\
		"leaq (%%r15,%%r15),%%rbx	\n\t"\
		"subq	%%r15,%%rax		\n\t"\
		"subq	%%rax,%%r15		\n\t"\
		"cmovcq %%rbx,%%r15		\n\t"\
		"\n\t"\
	"twopmodq64_q8_pshiftjmp:					\n\t"\
		"/* } endif((pshift >> j) & (uint64)1) */	\n\t"\
		"subq	$1,%%rcx	/* j-- */		\n\t"\
		"cmpq	$0,%%rcx	/* compare j vs 0 */\n\t"\
		"jge	LoopBeg8	/* if (j >= 0), Loop */	\n\t"\
	"LoopEnd8:							\n\t"\
		"movq	%%r8 ,%[__x0]	\n\t"\
		"movq	%%r9 ,%[__x1]	\n\t"\
		"movq	%%r10,%[__x2]	\n\t"\
		"movq	%%r11,%[__x3]	\n\t"\
		"movq	%%r12,%[__x4]	\n\t"\
		"movq	%%r13,%[__x5]	\n\t"\
		"movq	%%r14,%[__x6]	\n\t"\
		"movq	%%r15,%[__x7]	\n\t"\
		:	/* outputs: none */\
		: [__q0] "m" (q0)	/* All inputs from memory addresses here */\
		 ,[__q1] "m" (q1)	\
		 ,[__q2] "m" (q2)	\
		 ,[__q3] "m" (q3)	\
		 ,[__q4] "m" (q4)	\
		 ,[__q5] "m" (q5)	\
		 ,[__q6] "m" (q6)	\
		 ,[__q7] "m" (q7)	\
		 ,[__qinv0] "m" (qinv0)	\
		 ,[__qinv1] "m" (qinv1)	\
		 ,[__qinv2] "m" (qinv2)	\
		 ,[__qinv3] "m" (qinv3)	\
		 ,[__qinv4] "m" (qinv4)	\
		 ,[__qinv5] "m" (qinv5)	\
		 ,[__qinv6] "m" (qinv6)	\
		 ,[__qinv7] "m" (qinv7)	\
		 ,[__x0] "m" (x0)	\
		 ,[__x1] "m" (x1)	\
		 ,[__x2] "m" (x2)	\
		 ,[__x3] "m" (x3)	\
		 ,[__x4] "m" (x4)	\
		 ,[__x5] "m" (x5)	\
		 ,[__x6] "m" (x6)	\
		 ,[__x7] "m" (x7)	\
		 ,[__pshift] "m" (pshift)	\
		 ,[__j] "m" (j)	/* Only need this if debug and explicit loop enabled, but try with/sans for each version of the asm, pivk faster one. */\
		 ,[__start_index] "m" (start_index)	\
		: "cc","memory","rax","rbx","rcx","rdx","r8","r9","r10","r11","r12","r13","r14","r15"		/* Clobbered registers */\
		);

	/*...Double and return.	These are specialized for the case where 2^p == 1 mod q implies divisibility, in which case x = (q+1)/2.
	*/

	*checksum2 += x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7;

	r = 0;
	if(x0+x0-q0 == 1) r +=  1;
	if(x1+x1-q1 == 1) r +=  2;
	if(x2+x2-q2 == 1) r +=  4;
	if(x3+x3-q3 == 1) r +=  8;
	if(x4+x4-q4 == 1) r += 16;
	if(x5+x5-q5 == 1) r += 32;
	if(x6+x6-q6 == 1) r += 64;
	if(x7+x7-q7 == 1) r +=128;
	return(r);
}

#endif /* YES_ASM */

/*
Testing 64-bit factors...
twopmodq64_q4 : x0 = 2296955013172626669
twopmodq64_q4 : x0 = 10409682588550727410
twopmodq64_q4 : x0 = 880670618513687225
twopmodq64_q4 : x0 = 2866340383748859111
twopmodq64_q4 : x0 = 5073482287602503363
twopmodq64_q4 : x0 = 6951418682011146289
twopmodq64_q4 : x0 = 6509011776062795213
twopmodq64_q4 : x0 = 272552355776550471
twopmodq64_q4 : x0 = 4838859802115059419
twopmodq64_q4 : x0 = 1529997424122302037
twopmodq64_q4 : x0 = 512414035628167992
twopmodq64_q4 : x0 = 6146322359422961128
twopmodq64_q4 : x0 = 8278789907500224298
twopmodq64_q4 : x0 = 4711846519350180584
twopmodq64_q4 : x0 = 7430191046458322498
twopmodq64_q4 : x0 = 4067943528132455686
twopmodq64_q4 : x0 = 6760717124811562117
twopmodq64_q4 : x0 = 9342645731491044701
twopmodq64_q4 : x0 = 5250182887251718101
twopmodq64_q4 : x0 = 1
rm: t17922587: No such file or directory
*/

