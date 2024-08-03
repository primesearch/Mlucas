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

/*******************************************************************************
   We now include this header file if it was not included before.
*******************************************************************************/
#ifndef f2psp_h_included
#define f2psp_h_included

#ifdef __cplusplus
extern "C" {
#endif

	#define MI64_IS_DIV_BY_SCALAR32P_X8_SSE2(\
		array_64x8inputs,\
		q,		\
		qinv,	\
		retval	\
	)\
	{\
		DBG_ASSERT(qinv == qinv*((uint32)2 - q*qinv), "mi64_is_div_by_scalar32p: bad qinv!");\
		DBG_ASSERT(((uint32)&a[0] & 0x3f) == 0, "A-array not 64-byte aligned!");\
		__asm	mov	eax, array_64x8inputs	/* Assumes inputs a,b,c,d,... are 64-bit separated and &a[0} is 64-byte aligned */\
		__asm	lea	ebx, q\
		__asm	lea	ecx, qinv\
		__asm	movaps	xmm0,[eax     ]	/* ab: d3210 = [bhi|blo|ahi|alo] */\
		__asm	movaps	xmm1,[eax+0x10]	/* cd: d3210 = [dhi|dlo|chi|clo] */\
		__asm	movaps	xmm2,[eax+0x20]	/* ef: d3210 = [fhi|flo|ehi|elo] */\
		__asm	movaps	xmm3,[eax+0x30]	/* gh: d3210 = [hhi|hlo|ghi|glo] */\
		__asm	movaps	xmm6,xmm0	/* Circularly-permute [4,6,7] -> [6,7,4] here so the 2 packed outputs end up in xmm6,7 */\
		__asm	movaps	xmm5,xmm1\
		__asm	movaps	xmm7,xmm2\
		__asm	movaps	xmm4,xmm3\
		__asm	psrlq	xmm6, 32		/* d3210 = [  0|bhi|  0|ahi] */\
		__asm	psrlq	xmm5, 32		/* d3210 = [  0|dhi|  0|chi] */\
		__asm	psrlq	xmm7, 32		/* d3210 = [  0|fhi|  0|ehi] */\
		__asm	psrlq	xmm4, 32		/* d3210 = [  0|hhi|  0|ghi] */\
		__asm	psllq	xmm5, 32		/* d3210 = [dhi|  0|chi|  0] */\
		__asm	psllq	xmm4, 32		/* d3210 = [hhi|  0|ghi|  0] */\
		__asm	paddd	xmm6,xmm5		/* d3210 = [dhi|bhi|chi|ahi], xmm5 FREE */\
		__asm	paddd	xmm7,xmm4		/* d3210 = [hhi|fhi|ghi|ehi], xmm4 FREE */\
		__asm	movd	xmm4,[ebx]\
		__asm	movd	xmm5,[ecx]\
		__asm	pshufd	xmm4,xmm4,0x44	/* Broadcast q    to slots 0,2 of xmm4 */\
		__asm	pshufd	xmm5,xmm5,0x44	/* Broadcast qinv to slots 0,2 of xmm5 */\
		/* (a-h)[0]*qinv; Alas SSE2 has no 32-bit low-half packed MUL, so use 32x32->64 -bit and discard high halves */\
		__asm	pmuludq	xmm0,xmm5\
		__asm	pmuludq	xmm1,xmm5\
		__asm	pmuludq	xmm2,xmm5\
		__asm	pmuludq	xmm3,xmm5\
		/* cy[0-7] = MULH32(tmp[0-7]*q) - high halves of above MULQs automatically get overwritten: */\
		__asm	pmuludq	xmm0,xmm4\
		__asm	pmuludq	xmm1,xmm4\
		__asm	pmuludq	xmm2,xmm4\
		__asm	pmuludq	xmm3,xmm4\
		__asm	psrlq	xmm0, 32		/* d3210 = [  0|cy1|  0|cy0] */\
		__asm	psrlq	xmm1, 32		/* d3210 = [  0|cy3|  0|cy2] */\
		__asm	psrlq	xmm2, 32		/* d3210 = [  0|cy5|  0|cy4] */\
		__asm	psrlq	xmm3, 32		/* d3210 = [  0|cy7|  0|cy6] */\
		__asm	psllq	xmm1, 32		/* d3210 = [cy3|  0|cy2|  0] */\
		__asm	psllq	xmm3, 32		/* d3210 = [cy7|  0|cy6|  0] */\
		__asm	paddd	xmm0,xmm1		/* d3210 = [cy3|cy1|cy2|cy0], xmm1 FREE */\
		__asm	paddd	xmm2,xmm3		/* d3210 = [cy7|cy5|cy6|cy4], xmm3 FREE */\
		__asm	movaps	xmm3,xmm6		/* Copy of acbd[1] */\
		__asm	movaps	xmm1,xmm7		/* Copy of efgh[1] */\
		__asm	psubd	xmm6,xmm0		/* acbd[1] - cy0213, xmm0 FREE */\
		__asm	psubd	xmm7,xmm2		/* egfh[1] - cy4657, xmm2 FREE */\
		__asm	movaps	xmm2,xmm6		/* Copy of acbd[1] - cy0213 */\
		__asm	movaps	xmm0,xmm7		/* Copy of efgh[1] - cy4657 */\
		/* Had a borrow? Frickin' SSE2 only gives us signed packed-integer compares,\
		so need to emulate unsigned (x > y) via signed (x ^ 0x80000000) < (y ^ 0x80000000): */\
		__asm	pcmpeqd	xmm4,xmm4		/* All 1s  - will need to restore q to this register later */\
		__asm	pslld	xmm4, 31		/* 4-way 0x80000000 */\
		__asm	pxor	xmm6,xmm4		/* (acbd[1]-cy0213) ^ 0x80000000 */\
		__asm	pxor	xmm7,xmm4		/* (egfh[1]-cy4657) ^ 0x80000000 */\
		__asm	pxor	xmm3,xmm4		/* (acbd[1]) ^ 0x80000000 */\
		__asm	pxor	xmm1,xmm4		/* (egfh[1]) ^ 0x80000000 */\
		__asm	pcmpgtd	xmm6,xmm3		/* cy0213 = (acbd[1]-cy0213) > abcd[1], xmm3 FREE */\
		__asm	pcmpgtd	xmm7,xmm1		/* cy4657 = (egfh[1]-cy4657) > efgh[1], xmm1 FREE */\
		__asm	pshufd	xmm3,xmm2,0x31	/* xmm2 = [----|tmp1|----|tmp0], xmm3 = [----|tmp3|----|tmp2], don't care what's in ---- slots */\
		__asm	pshufd	xmm1,xmm0,0x31	/* xmm0 = [----|tmp5|----|tmp4], xmm1 = [----|tmp7|----|tmp6], don't care what's in ---- slots */\
		__asm	movd	xmm4,[ebx]		/* Restore q to xmm4 */\
		__asm	pshufd	xmm4,xmm4,0x44	/* Broadcast q    to slots 0,2 of xmm4 */\
		/* tmp[0-7]*qinv; Alas SSE2 has no 32-bit low-half packed MUL, so use 32x32->64 -bit and discard high halves */\
		__asm	pmuludq	xmm3,xmm5\
		__asm	pmuludq	xmm1,xmm5\
		__asm	pmuludq	xmm2,xmm5\
		__asm	pmuludq	xmm0,xmm5\
		/* Add carries 01/45, scatter carries 23/67 into slots of 01/45, add those...Since SSE2 compare result is ~()ed, add really means sub: */\
		__asm	psubd	xmm2,xmm6		/* xmm6 = [----|tmp1|----|tmp0], don't care what's in ---- slots */\
		__asm	psubd	xmm0,xmm7		/* xmm7 = [----|tmp5|----|tmp4], don't care what's in ---- slots */\
		__asm	pshufd	xmm6,xmm6,0x31\
		__asm	pshufd	xmm7,xmm7,0x31\
		__asm	psubd	xmm3,xmm6		/* xmm3 = [----|tmp3|----|tmp2], don't care what's in ---- slots */\
		__asm	psubd	xmm1,xmm7		/* xmm1 = [----|tmp7|----|tmp6], don't care what's in ---- slots */\
		/* cy[0-7] = MULH32(tmp[0-7]*q) - high halves of above MULQs automatically get overwritten: */\
		__asm	pmuludq	xmm2,xmm4\
		__asm	pmuludq	xmm0,xmm4\
		__asm	pmuludq	xmm3,xmm4\
		__asm	pmuludq	xmm1,xmm4\
		__asm	psrlq	xmm2, 32		/* d3210 = [  0|cy1|  0|cy0] */\
		__asm	psrlq	xmm0, 32		/* d3210 = [  0|cy5|  0|cy4] */\
		__asm	psrlq	xmm3, 32		/* d3210 = [  0|cy3|  0|cy2] */\
		__asm	psrlq	xmm1, 32		/* d3210 = [  0|cy7|  0|cy6] */\
		__asm	pshufd	xmm2,xmm2,0x58	/* [  0|  0|cy1|cy0] */\
		__asm	pshufd	xmm0,xmm0,0x58	/* [  0|  0|cy5|cy4] */\
		__asm	pshufd	xmm3,xmm3,0x85	/* [cy3|cy2|  0|  0] */\
		__asm	pshufd	xmm1,xmm1,0x85	/* [cy7|cy6|  0|  0] */\
		__asm	paddd	xmm2,xmm3		/* d3210 = [cy3|cy1|cy2|cy0] */\
		__asm	paddd	xmm0,xmm1		/* d3210 = [cy7|cy5|cy6|cy4] */\
		__asm	pcmpgtd	xmm7,xmm7		/* All 0s */\
		__asm	pcmpeqd	xmm2,xmm7		/* retval[0-3] */\
		__asm	pcmpeqd	xmm0,xmm7		/* retval[4-7] */\
		__asm	movmskps eax,xmm2		/* retval[0-3] */\
		__asm	movmskps ebx,xmm0		/* retval[4-7] */\
		__asm	shl		 ebx, 4		/* retval[4-7] << 4 */\
		__asm	add		 eax,ebx	/* retval[0-7] */\
		__asm	mov	retval,  eax	\
	}

#ifdef __cplusplus
}
#endif

#endif	/* f2psp_h_included */

