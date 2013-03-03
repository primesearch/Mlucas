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

#include "Mlucas.h"

/* Use for testing higher-accuracy version of the twiddles computation */
#define HIACC 1

#ifdef USE_SSE2

	#ifdef COMPILER_TYPE_MSVC
		#include "sse2_macro.h"
	#endif

	#undef DEBUG_SSE2
//	#define DEBUG_SSE2

	#if defined(COMPILER_TYPE_MSVC)

		/* Full-inline-asm version of SSE2_RADIX8_DIT_0TWIDDLE_B. Assumes base addresses add4-7 in e[a-d]x,
		base offset in edi, these must remain unchanged, so esi is available for temporary storage. */
		#define SSE2_RADIX8_DIT_0TWIDDLE_B(__tmp)\
		{\
			/*** 2nd of the 2 length-4 subtransforms gets done first, due to e.g. t1-+t9 combos in final step: ***/\
			/*\
			t9 =a[j1+p4];	t10=a[j2+p4];\
			t11=a[j1+p5];	t12=a[j2+p5];\
			~t11=t9 -t11;	~t12=t10-t12;\
			~t9 =t9 +t11;	~t10=t10+t12;\
			*/\
			__asm	mov	esi,__tmp\
			__asm	movaps	xmm0,[eax]		/* xmm0 <- a[j1+p4] = t9 */\
			__asm	movaps	xmm1,[eax+0x10]	/* xmm1 <- a[j2+p4] = t10*/\
			__asm	movaps	xmm2,xmm0		/* xmm4 <- copy of t9    */\
			__asm	movaps	xmm3,xmm1		/* xmm5 <- copy of t10   */\
			__asm	addpd	xmm2,[ebx]		/* xmm2 <- t9  */\
			__asm	addpd	xmm3,[ebx+0x10]	/* xmm3 <- t10 */\
			__asm	subpd	xmm0,[ebx]		/* xmm0 <- t11 */\
			__asm	subpd	xmm1,[ebx+0x10]	/* xmm1 <- t12 */\
			/*\
			t13=a[j1+p6];	t14=a[j2+p6];\
			rt =a[j1+p7];	it =a[j2+p7];\
			~t15=t13-t15	~t16=t14-t16\
			~t13=t13+t15	~t14=t14+t16\
			*/\
			__asm	movaps	xmm4,[ecx]		/* xmm4 <- a[j1+p6] = t13*/\
			__asm	movaps	xmm5,[ecx+0x10]	/* xmm5 <- a[j2+p6] = t14*/\
			__asm	movaps	xmm6,xmm4		/* xmm6 <- copy of t13   */\
			__asm	movaps	xmm7,xmm5		/* xmm7 <- copy of t14   */\
			__asm	addpd	xmm6,[edx]		/* xmm6 <- t13 */\
			__asm	addpd	xmm7,[edx+0x10]	/* xmm7 <- t14 */\
			__asm	subpd	xmm4,[edx]		/* xmm4 <- t15 */\
			__asm	subpd	xmm5,[edx+0x10]	/* xmm5 <- t16 */\
			/* Copy t13,14 into temp-array slot add6 */\
			__asm	movaps	[esi+0xc0],xmm6\
			__asm	movaps	[esi+0xd0],xmm7\
			\
			/** GPRs: ***** SSE Regs: ***** Temp array: ********\\
			*	eax, add4	xmm0 <- t11		add0 <- unused		*\
			*	ebx, add5	xmm1 <- t12		add1 <- unused 		*\
			*	ecx, add6	xmm2 <- t9 		add2 <- unused		*\
			*	edx, add7	xmm3 <- t10		add3 <- unused 		*\
			*				xmm4 <- t15		add4 <- unused		*\
			*				xmm5 <- t16		add5 <- unused		*\
			*				xmm6 <- t13		add6 <- t13,14		*\
			*				xmm7 <- t14		add7 <- t15,16		*\
			\***************************************************/\
			\
			/*\
			rt =t13;	t13=t9 -rt ;	t9 =t9 +rt ;	copies of t13 in add6     , xmm6\
			it =t14;	t14=t10-it ;	t10=t10+it ;	copies of t14 in add6+0x10, xmm7\
			\
			rt =t15;	t15=t11-t16;	t11=t11+t16;	copies of t15 in add7     , xmm4\
						t16=t12+rt ;	t12=t12-rt ;	copies of t16 in add7+0x10, xmm5\
			*/\
			/* Move outputs t11,12 into a[j1,j2+p5], first doing the addsub and mul by ISRT2: */\
			/* Move outputs t15,16 into a[j1,j2+p7], first doing the addsub and mul by ISRT2: */\
			\
			__asm	addpd	xmm6,xmm2		/* xmm6 <- ~t9  */\
			__asm	addpd	xmm7,xmm3		/* xmm7 <- ~t10 */\
			__asm	subpd	xmm2,[esi+0xc0]	/* xmm2 <- ~t13 */\
			__asm	subpd	xmm3,[esi+0xd0]	/* xmm3 <- ~t14 */\
			/* Move t13,14 into a[j1,j2+p6] */\
			__asm	movaps	[esi+0xc0],xmm2	/* add6r <- ~t13 */\
			__asm	movaps	[esi+0xd0],xmm3	/* add6i <- ~t14 */\
			\
			__asm	movaps	xmm2,xmm4	/* xmm2 <- copy of t15 */\
			__asm	movaps	xmm3,xmm5	/* xmm3 <- copy of t16 */\
			__asm	addpd	xmm5,xmm0	/* xmm5 <-~t11 */\
			__asm	subpd	xmm0,xmm3	/* xmm0 <-~t15 */\
			__asm	addpd	xmm4,xmm1	/* xmm4 <-~t16 */\
			__asm	subpd	xmm1,xmm2	/* xmm1 <-~t12 */\
			\
			__asm	mov	esi,isrt2\
			__asm	movaps	xmm2,xmm5	/* xmm2 <- copy of~t11 */\
			__asm	movaps	xmm3,xmm1	/* xmm3 <- copy of~t12 */\
			__asm	addpd	xmm5,xmm1	/* xmm5 <-~(t11+t12), xmm1 FREE */\
			__asm	movaps	xmm1,[esi]	/* xmm1 <- ISRT2 */\
			__asm	subpd	xmm2,xmm3	/* xmm2 <-~(t11-t12), xmm3 FREE */\
			__asm	mov	esi,__tmp\
			__asm	mulpd	xmm5,xmm1	/* xmm5 <- (t11+t12)*ISRT2 */\
			__asm	mulpd	xmm2,xmm1	/* xmm2 <- (t11-t12)*ISRT2 */\
			__asm	movaps	xmm3,xmm0	/* xmm3 <- copy of~t15 */\
			\
			__asm	movaps	[esi+0xa0],xmm5	/* add5r<- (t11+t12)*ISRT2, xmm5 FREE */\
			__asm	movaps	xmm5,xmm4	/* xmm5 <- copy of~t16 */\
			__asm	addpd	xmm0,xmm4	/* xmm0 <-~(t15+t16) */\
			__asm	movaps	[esi+0xb0],xmm2	/* add5i<- (t11-t12)*ISRT2 */\
			__asm	subpd	xmm3,xmm5	/* xmm3 <-~(t15-t16) */\
			__asm	mulpd	xmm0,xmm1	/* xmm0 <- (t15+t16)*ISRT2 */\
			__asm	mulpd	xmm3,xmm1	/* xmm3 <- (t15-t16)*ISRT2 */\
			__asm	movaps	[esi+0xe0],xmm0	/* add7r<- (t15+t16)*ISRT2 */\
			__asm	movaps	[esi+0xf0],xmm3	/* add7i<- (t15-t16)*ISRT2 */\
			\
			/** GPRs: ***** SSE Regs: ***** Temp array: ******************\\
			*    eax, add4   xmm0 unused        add0 <- unused            *\
			*    ebx, add5   xmm1 unused        add1 <- unused            *\
			*    ecx, add6   xmm2 unused        add2 <- unused            *\
			*    edx, add7   xmm3 unused        add3 <- unused            *\
			*                xmm4 unused        add4 <- unused            *\
			*                xmm5 unused        add5 <- (t11+-t12)*ISRT2  *\
			*                xmm6 t9            add6 <-  t13,14           *\
			*                xmm7 t10           add7 <- (t15+-t16)*ISRT2  *\
			\*************************************************************/\
			\
			/**************** 1st of the 2 length-4 subtransforms... **************/\
			/*\
			t1 =a[j1   ];	t2 =a[j2   ];\
			rt =a[j1+p1];	it =a[j2+p1];\
			t3 =t1 -rt;		t4 =t2 -it;\
			t1 =t1 +rt;		t2 =t2 +it;\
			*/\
			__asm	sub	eax, edi\
			__asm	sub	ebx, edi\
			__asm	sub	ecx, edi\
			__asm	sub	edx, edi\
			\
			__asm	movaps	xmm0,[eax]		/* xmm0 <- a[j1   ] = t1 */\
			__asm	movaps	xmm1,[eax+0x10]	/* xmm1 <- a[j2   ] = t2 */\
			__asm	movaps	xmm2,xmm0		/* xmm2 <- copy of t1    */\
			__asm	movaps	xmm3,xmm1		/* xmm3 <- copy of t2    */\
			__asm	addpd	xmm2,[ebx]		/*~xmm2 <- t1  */\
			__asm	addpd	xmm3,[ebx+0x10]	/*~xmm3 <- t2  */\
			__asm	subpd	xmm0,[ebx]		/*~xmm0 <- t3  */\
			__asm	subpd	xmm1,[ebx+0x10]	/*~xmm1 <- t4  */\
			\
			/* Move t9,10 into temp[j1+p0] in anticipation of final outputs t1+-t9, t2+-t10 which will go there: */\
			__asm	movaps	[esi     ],xmm6\
			__asm	movaps	[esi+0x10],xmm7	/* add0 <-  t9,t10 */\
			__asm	addpd	xmm6,xmm6		/* xmm6 <- 2*t9  */\
			__asm	addpd	xmm7,xmm7		/* xmm7 <- 2*t10 */\
			/* Move 2*t9,10 into temp[j1+p4] in anticipation of final outputs t1+-t9, t2+-t10 which will go there: */\
			__asm	movaps	[esi+0x80],xmm6	/* add4 <- 2*t9,2*t10  */\
			__asm	movaps	[esi+0x90],xmm7	/* xmm4-7 FREE */\
			/*\
			t5 =a[j1+p2];	t6 =a[j2+p2];\
			rt =a[j1+p3];	it =a[j2+p3];\
			t7 =t5 -rt;		t8 =t6 -it;\
			t5 =t5 +rt;		t6 =t6 +it;\
			*/\
			__asm	movaps	xmm4,[ecx]		/* xmm4 <- a[j1+p2] = t5 */\
			__asm	movaps	xmm5,[ecx+0x10]	/* xmm5 <- a[j2+p2] = t6 */\
			__asm	movaps	xmm6,xmm4		/* xmm6 <- copy of t5    */\
			__asm	movaps	xmm7,xmm5		/* xmm7 <- copy of t6    */\
			__asm	addpd	xmm6,[edx]		/*~xmm6 <- t5  */\
			__asm	addpd	xmm7,[edx+0x10]	/*~xmm7 <- t6  */\
			__asm	subpd	xmm4,[edx]		/*~xmm4 <- t7  */\
			__asm	subpd	xmm5,[edx+0x10]	/*~xmm5 <- t8  */\
			/* Copy t5,6 into temp-array slots a[j1,j2+p2] */\
			__asm	movaps	[esi+0x40],xmm6\
			__asm	movaps	[esi+0x50],xmm7	/* add2 <-  t5,t6 */\
			\
			/** GPRs: ***** SSE Regs: ***** Temp array: ******************\\
			*    eax, add0   xmm0 t3            add0 <-  t9,t10           *\
			*    ebx, add1   xmm1 t4            add1 <- unused            *\
			*    ecx, add2   xmm2 t1            add2 <-  t5,t6            *\
			*    edx, add3   xmm3 t2            add3 <- unused            *\
			*                xmm4 t7            add4 <- 2*t9,2*t10        *\
			*                xmm5 t8            add5 <- (t11+-t12)*ISRT2  *\
			*                xmm6 t5            add6 <-  t13,14           *\
			*                xmm7 t6            add7 <- (t15+-t16)*ISRT2  *\
			\*************************************************************/\
			\
			/*\
			rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;	copies of t5 in add2     , xmm6\
			it =t6;	t6 =t2 -it;	t2 =t2 +it;	copies of t6 in add2+0x10, xmm7\
			\
			rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;\
					t8 =t4 +rt;	t4 =t4 -rt;\
			*/\
			__asm	addpd	xmm6,xmm2		/* xmm6 <- ~t1 */\
			__asm	addpd	xmm7,xmm3		/* xmm7 <- ~t2 */\
			__asm	subpd	xmm2,[esi+0x40]	/* xmm2 <- ~t5 */\
			__asm	subpd	xmm3,[esi+0x50]	/* xmm3 <- ~t6 */\
			\
			/* Compute and dump first 2 outputs now, in order to free up 2 registers: */\
			__asm	addpd	xmm6,[esi     ]	/* t1+t9 */\
			__asm	addpd	xmm7,[esi+0x10]	/* t2+t10*/\
			__asm	movaps	[esi     ],xmm6	/* a[j1   ], DONE. */\
			__asm	movaps	[esi+0x10],xmm7	/* a[j2   ], DONE. */\
			\
			__asm	subpd	xmm6,[esi+0x80]	/* t1-t9  = [t1+t9 ] - 2*t9  */\
			__asm	subpd	xmm7,[esi+0x90]	/* t2-t10 = [t2+t10] - 2*t10 */\
			__asm	movaps	[esi+0x80],xmm6	/* a[j1+p4], DONE. */\
			__asm	movaps	[esi+0x90],xmm7	/* a[j2+p4], DONE. */\
			\
			__asm	movaps	xmm6,xmm4		/* xmm6 <- copy of t7 */\
			__asm	movaps	xmm7,xmm5		/* xmm7 <- copy of t8 */\
			__asm	addpd	xmm5,xmm0	/* xmm5 <- ~t3 */\
			__asm	subpd	xmm0,xmm7	/* xmm0 <- ~t7 */\
			__asm	addpd	xmm4,xmm1	/* xmm4 <- ~t8 */\
			__asm	subpd	xmm1,xmm6	/* xmm1 <- ~t4 */\
			 \
			/** GPRs: ***** SSE Regs: ***** Temp array: ******************\\
			*    eax, add0   xmm0 t7            add0 <- DONE              *\
			*    ebx, add1   xmm1 t4            add1 <- unused            *\
			*    ecx, add2   xmm2 t5            add2 <- unused            *\
			*    edx, add3   xmm3 t6            add3 <- unused            *\
			*    edi, p4     xmm4 t8            add4 <- DONE              *\
			*    esi, add4   xmm5 t3            add5 <- (t11+-t12)*ISRT2  *\
			*                xmm6 unused        add6 <-  t13,14           *\
			*                xmm7 unused        add7 <- (t15+-t16)*ISRT2  *\
			\*************************************************************/\
			\
			/*\
			rt=(t11+t12)*ISRT2;			it=(t11-t12)*ISRT2;	precomputed, in add5\
			a[j1+p1]=t3+rt;				a[j2+p1]=t4-it;\
			a[j1+p5]=t3-rt;				a[j2+p5]=t4+it;\
			*/\
			__asm	movaps	xmm6,xmm5		/* xmm6 <- copy of t3 */\
			__asm	movaps	xmm7,xmm1		/* xmm7 <- copy of t4 */\
			__asm	addpd	xmm5,[esi+0xa0]	/* t3+rt */\
			__asm	subpd	xmm1,[esi+0xb0]	/* t4-it */\
			__asm	subpd	xmm6,[esi+0xa0]	/* t3-rt */\
			__asm	addpd	xmm7,[esi+0xb0]	/* t4+it */\
			__asm	movaps	[esi+0x20],xmm5	/* a[j1+p1] */\
			__asm	movaps	[esi+0x30],xmm1	/* a[j2+p1] */\
			__asm	movaps	[esi+0xa0],xmm6	/* a[j1+p5] */\
			__asm	movaps	[esi+0xb0],xmm7	/* a[j2+p5] */\
			\
			/*\
			a[j1+p2]=t5+t14;			a[j1+p2]=t6-t13;\
			a[j1+p6]=t5-t14;			a[j1+p6]=t6+t13;\
			*/\
			__asm	movaps	xmm6,xmm2		/* xmm6 <- copy of t5 */\
			__asm	movaps	xmm7,xmm3		/* xmm7 <- copy of t6 */\
			__asm	addpd	xmm2,[esi+0xd0]	/* t5+t14*/\
			__asm	subpd	xmm3,[esi+0xc0]	/* t6-t13*/\
			__asm	subpd	xmm6,[esi+0xd0]	/* t5-t14*/\
			__asm	addpd	xmm7,[esi+0xc0]	/* t6+t13*/\
			__asm	movaps	[esi+0x40],xmm2	/* a[j1+p2] */\
			__asm	movaps	[esi+0x50],xmm3	/* a[j2+p2] */\
			__asm	movaps	[esi+0xc0],xmm6	/* a[j1+p6] */\
			__asm	movaps	[esi+0xd0],xmm7	/* a[j2+p6] */\
			\
			/*\
			rt=(t15-t16)*ISRT2;			it=(t15+t16)*ISRT2;	precomputed, it,rt] in add7; NOTE reversed order!\
			a[j1+p3]=t7-rt;				a[j1+p3]  =t8-it;\
			a[j1+p7]=t7+rt;				a[j1+p7]  =t8+it;\
			*/\
			__asm	movaps	xmm6,xmm0		/* xmm6 <- copy of t7 */\
			__asm	movaps	xmm7,xmm4		/* xmm7 <- copy of t8 */\
			__asm	subpd	xmm0,[esi+0xf0]	/* t7-rt */\
			__asm	subpd	xmm4,[esi+0xe0]	/* t8-it */\
			__asm	addpd	xmm6,[esi+0xf0]	/* t7+rt */\
			__asm	addpd	xmm7,[esi+0xe0]	/* t8+it */\
			__asm	movaps	[esi+0x60],xmm0	/* a[j1+p3] */\
			__asm	movaps	[esi+0x70],xmm4	/* a[j2+p3] */\
			__asm	movaps	[esi+0xe0],xmm6	/* a[j1+p7] */\
			__asm	movaps	[esi+0xf0],xmm7	/* a[j2+p7] */\
													/* Totals: 97 load/store [61 movaps, 36 implied], 54 add/subpd,  4 mulpd, 59 address-compute */\
		}

	#else	/* GCC-style inline ASM: */

		#if OS_BITS == 32

			#include "radix32_dif_dit_pass_gcc32.h"

		#else

			#include "radix32_dif_dit_pass_gcc64.h"

		#endif

	#endif

#endif

/***************/

/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform a radix-32 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details.
*/
void radix32_dif_pass(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr, int init_sse2, int thr_id)
{
	static int max_threads = 0;
	int i,j,j1,j2,jlo,jhi,m,iroot_prim,iroot,k1,k2;
	int p01,p02,p03,p04,p08,p0C,p10,p14,p18,p1C;
	const double c = 0.92387953251128675613, s     = 0.38268343236508977173	/* exp[  i*(twopi/16)]	*/
			,c32_1 = 0.98078528040323044912, s32_1 = 0.19509032201612826784	/* exp(  i*twopi/32), the radix-32 fundamental sincos datum	*/
			,c32_3 = 0.83146961230254523708, s32_3 = 0.55557023301960222473;/* exp(3*i*twopi/32)	*/
	double rt,it,re0,im0,re1,im1;

#ifdef USE_SSE2

  #if !(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))
	#error SSE2 code not supported for this compiler!
  #endif

	static struct complex *sc_arr = 0x0, *sc_ptr;
	double *add0;	/* Addresses into array sections */
#ifdef COMPILER_TYPE_MSVC
	double *add1, *add2, *add3, *add4, *add5, *add6, *add7;
#endif
	struct complex *c_tmp,*s_tmp;

  #ifdef MULTITHREAD
	static struct complex *__r0;					// Base address for discrete per-thread local stores
	// In || mode, only above base-pointer (shared by all threads) is static:
	struct complex *isrt2, *cc0, *ss0, *cc1, *ss1, *cc3, *ss3, *two, *r00;
  #elif defined(COMPILER_TYPE_GCC)
	static struct complex *isrt2, *cc0, *ss0, *cc1, *ss1, *cc3, *ss3, *two, *r00;
  #else
	static struct complex *isrt2, *two, *cc0, *ss0, *cc1, *ss1, *cc3, *ss3
		,*c00,*c01,*c02,*c03,*c04,*c05,*c06,*c07,*c08,*c09,*c0A,*c0B,*c0C,*c0D,*c0E,*c0F
		,*c10,*c11,*c12,*c13,*c14,*c15,*c16,*c17,*c18,*c19,*c1A,*c1B,*c1C,*c1D,*c1E,*c1F
		,*s00,*s01,*s02,*s03,*s04,*s05,*s06,*s07,*s08,*s09,*s0A,*s0B,*s0C,*s0D,*s0E,*s0F
		,*s10,*s11,*s12,*s13,*s14,*s15,*s16,*s17,*s18,*s19,*s1A,*s1B,*s1C,*s1D,*s1E,*s1F
		,*r00,*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09,*r0A,*r0B,*r0C,*r0D,*r0E,*r0F
		,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r1A,*r1B,*r1C,*r1D,*r1E,*r1F
		,*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r2A,*r2B,*r2C,*r2D,*r2E,*r2F
		,*r30,*r31,*r32,*r33,*r34,*r35,*r36,*r37,*r38,*r39,*r3A,*r3B,*r3C,*r3D,*r3E,*r3F;
  #endif

#else

	int jp,jt;
	double *addr, *addp;
	int prefetch_offset;
	double	 c01,c02,c03,c04,c05,c06,c07,c08,c09,c0A,c0B,c0C,c0D,c0E,c0F
		,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c1A,c1B,c1C,c1D,c1E,c1F
		    ,s01,s02,s03,s04,s05,s06,s07,s08,s09,s0A,s0B,s0C,s0D,s0E,s0F
		,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s1A,s1B,s1C,s1D,s1E,s1F
		,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F
		,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F
		,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F;

#endif

#ifdef USE_SSE2

	/* In order to eschew complex thread-block-and-sync logic related to the local-store-init step in multithread mode,
	switch to a special-init-mode-call paradigm, in which the function is inited once (for as many threads as needed)
	prior to being executed:
	*/
	if(init_sse2 && !sc_ptr)	// Just check (init_sse2 != 0) here, to allow the *value* of init_sse2 to store #threads
	{
		max_threads = init_sse2;
	#ifndef COMPILER_TYPE_GCC
		ASSERT(HERE, NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
	#endif
		ASSERT(HERE, max_threads >= NTHREADS, "Multithreading requires max_threads >= NTHREADS!");
		ASSERT(HERE, sc_arr == 0x0, "Init-mode call conflicts with already-malloc'ed local storage!");
		ASSERT(HERE, thr_id == -1, "Init-mode call must be outside of any multithreading!");
		sc_arr = ALLOC_COMPLEX(sc_arr, 0x90*max_threads);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 64 16-byte slots of sc_arr for temporaries, next 7 for the nontrivial complex 32nd roots,
	last 64 for the doubled sincos twiddles, plus at least 3 more slots to allow for 64-byte alignment of the array.
	*/
		#ifdef MULTITHREAD
	//	if(max_threads > 1) {
			__r0  = sc_ptr;
			isrt2 = sc_ptr + 0x40;
			cc0	  = sc_ptr + 0x41;
			ss0	  = sc_ptr + 0x42;
			cc1	  = sc_ptr + 0x43;
			ss1	  = sc_ptr + 0x44;
			cc3	  = sc_ptr + 0x45;
			ss3	  = sc_ptr + 0x46;
			two   = sc_ptr + 0x87;
			for(i = 0; i < max_threads; ++i) {
				/* These remain fixed within each per-thread local store: */
				isrt2->re = ISRT2;	isrt2->im = ISRT2;
				cc0  ->re = c	;	cc0  ->im = c	;		ss0  ->re = s	;	ss0  ->im = s	;
				cc1  ->re = c32_1;	cc1  ->im = c32_1;		ss1  ->re = s32_1;	ss1  ->im = s32_1;
				cc3  ->re = c32_3;	cc3  ->im = c32_3;		ss3  ->re = s32_3;	ss3  ->im = s32_3;
				two  ->re = 2.0;	two  ->im = 2.0;
				/* Move on to next thread's local store */
				isrt2 += 0x90;
				cc0   += 0x90;
				ss0   += 0x90;
				cc1   += 0x90;
				ss1   += 0x90;
				cc3   += 0x90;
				ss3   += 0x90;
				two   += 0x90;
			}
		#elif defined(COMPILER_TYPE_GCC)
			r00   = sc_ptr;
			isrt2 = sc_ptr + 0x40;
			cc0	  = sc_ptr + 0x41;
			ss0	  = sc_ptr + 0x42;
			cc1	  = sc_ptr + 0x43;
			ss1	  = sc_ptr + 0x44;
			cc3	  = sc_ptr + 0x45;
			ss3	  = sc_ptr + 0x46;
			two   = sc_ptr + 0x87;
			/* These remain fixed: */
			isrt2->re = ISRT2;	isrt2->im = ISRT2;		two->re   = 2.0;	two->im   = 2.0;
			cc0  ->re = c	;	cc0  ->im = c	;		ss0  ->re = s	;	ss0  ->im = s	;
			cc1  ->re = c32_1;	cc1  ->im = c32_1;		ss1  ->re = s32_1;	ss1  ->im = s32_1;
			cc3  ->re = c32_3;	cc3  ->im = c32_3;		ss3  ->re = s32_3;	ss3  ->im = s32_3;
		#else
	//	} else {
			r00		= sc_ptr + 0x00;	c00		= sc_ptr + 0x47;
			r01		= sc_ptr + 0x01;	s00		= sc_ptr + 0x48;
			r02		= sc_ptr + 0x02;	c10		= sc_ptr + 0x49;
			r03		= sc_ptr + 0x03;	s10		= sc_ptr + 0x4a;
			r04		= sc_ptr + 0x04;	c08		= sc_ptr + 0x4b;
			r05		= sc_ptr + 0x05;	s08		= sc_ptr + 0x4c;
			r06		= sc_ptr + 0x06;	c18		= sc_ptr + 0x4d;
			r07		= sc_ptr + 0x07;	s18		= sc_ptr + 0x4e;
			r08		= sc_ptr + 0x08;	c04		= sc_ptr + 0x4f;
			r09		= sc_ptr + 0x09;	s04		= sc_ptr + 0x50;
			r0A		= sc_ptr + 0x0a;	c14		= sc_ptr + 0x51;
			r0B		= sc_ptr + 0x0b;	s14		= sc_ptr + 0x52;
			r0C		= sc_ptr + 0x0c;	c0C		= sc_ptr + 0x53;
			r0D		= sc_ptr + 0x0d;	s0C		= sc_ptr + 0x54;
			r0E		= sc_ptr + 0x0e;	c1C		= sc_ptr + 0x55;
			r0F		= sc_ptr + 0x0f;	s1C		= sc_ptr + 0x56;
			r10		= sc_ptr + 0x10;	c02		= sc_ptr + 0x57;
			r11		= sc_ptr + 0x11;	s02		= sc_ptr + 0x58;
			r12		= sc_ptr + 0x12;	c12		= sc_ptr + 0x59;
			r13		= sc_ptr + 0x13;	s12		= sc_ptr + 0x5a;
			r14		= sc_ptr + 0x14;	c0A		= sc_ptr + 0x5b;
			r15		= sc_ptr + 0x15;	s0A		= sc_ptr + 0x5c;
			r16		= sc_ptr + 0x16;	c1A		= sc_ptr + 0x5d;
			r17		= sc_ptr + 0x17;	s1A		= sc_ptr + 0x5e;
			r18		= sc_ptr + 0x18;	c06		= sc_ptr + 0x5f;
			r19		= sc_ptr + 0x19;	s06		= sc_ptr + 0x60;
			r1A		= sc_ptr + 0x1a;	c16		= sc_ptr + 0x61;
			r1B		= sc_ptr + 0x1b;	s16		= sc_ptr + 0x62;
			r1C		= sc_ptr + 0x1c;	c0E		= sc_ptr + 0x63;
			r1D		= sc_ptr + 0x1d;	s0E		= sc_ptr + 0x64;
			r1E		= sc_ptr + 0x1e;	c1E		= sc_ptr + 0x65;
			r1F		= sc_ptr + 0x1f;	s1E		= sc_ptr + 0x66;
			r20		= sc_ptr + 0x20;	c01		= sc_ptr + 0x67;
			r21		= sc_ptr + 0x21;	s01		= sc_ptr + 0x68;
			r22		= sc_ptr + 0x22;	c11		= sc_ptr + 0x69;
			r23		= sc_ptr + 0x23;	s11		= sc_ptr + 0x6a;
			r24		= sc_ptr + 0x24;	c09		= sc_ptr + 0x6b;
			r25		= sc_ptr + 0x25;	s09		= sc_ptr + 0x6c;
			r26		= sc_ptr + 0x26;	c19		= sc_ptr + 0x6d;
			r27		= sc_ptr + 0x27;	s19		= sc_ptr + 0x6e;
			r28		= sc_ptr + 0x28;	c05		= sc_ptr + 0x6f;
			r29		= sc_ptr + 0x29;	s05		= sc_ptr + 0x70;
			r2A		= sc_ptr + 0x2a;	c15		= sc_ptr + 0x71;
			r2B		= sc_ptr + 0x2b;	s15		= sc_ptr + 0x72;
			r2C		= sc_ptr + 0x2c;	c0D		= sc_ptr + 0x73;
			r2D		= sc_ptr + 0x2d;	s0D		= sc_ptr + 0x74;
			r2E		= sc_ptr + 0x2e;	c1D		= sc_ptr + 0x75;
			r2F		= sc_ptr + 0x2f;	s1D		= sc_ptr + 0x76;
			r30		= sc_ptr + 0x30;	c03		= sc_ptr + 0x77;
			r31		= sc_ptr + 0x31;	s03		= sc_ptr + 0x78;
			r32		= sc_ptr + 0x32;	c13		= sc_ptr + 0x79;
			r33		= sc_ptr + 0x33;	s13		= sc_ptr + 0x7a;
			r34		= sc_ptr + 0x34;	c0B		= sc_ptr + 0x7b;
			r35		= sc_ptr + 0x35;	s0B		= sc_ptr + 0x7c;
			r36		= sc_ptr + 0x36;	c1B		= sc_ptr + 0x7d;
			r37		= sc_ptr + 0x37;	s1B		= sc_ptr + 0x7e;
			r38		= sc_ptr + 0x38;	c07		= sc_ptr + 0x7f;
			r39		= sc_ptr + 0x39;	s07		= sc_ptr + 0x80;
			r3A		= sc_ptr + 0x3a;	c17		= sc_ptr + 0x81;
			r3B		= sc_ptr + 0x3b;	s17		= sc_ptr + 0x82;
			r3C		= sc_ptr + 0x3c;	c0F		= sc_ptr + 0x83;
			r3D		= sc_ptr + 0x3d;	s0F		= sc_ptr + 0x84;
			r3E		= sc_ptr + 0x3e;	c1F		= sc_ptr + 0x85;
			r3F		= sc_ptr + 0x3f;	s1F		= sc_ptr + 0x86;
			isrt2	= sc_ptr + 0x40;	two		= sc_ptr + 0x87;
			cc0		= sc_ptr + 0x41;
			ss0		= sc_ptr + 0x42;
			cc1		= sc_ptr + 0x43;
			ss1		= sc_ptr + 0x44;
			cc3		= sc_ptr + 0x45;
			ss3		= sc_ptr + 0x46;
			/* These remain fixed: */
			isrt2->re = ISRT2;	isrt2->im = ISRT2;		two->re   = 2.0;	two->im   = 2.0;
			cc0  ->re = c	;	cc0  ->im = c	;		ss0  ->re = s	;	ss0  ->im = s	;
			cc1  ->re = c32_1;	cc1  ->im = c32_1;		ss1  ->re = s32_1;	ss1  ->im = s32_1;
			cc3  ->re = c32_3;	cc3  ->im = c32_3;		ss3  ->re = s32_3;	ss3  ->im = s32_3;
		#endif
	//	}
		return;
	}	/* end of inits */

	/* If multithreaded, set the local-store pointers needed for the current thread; */
	#ifdef MULTITHREAD
		ASSERT(HERE, (uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
		r00 = __r0 + thr_id*0x90;
		cc0	= r00 + 0x41;
	#endif

#endif

	p01 = incr >> 5;
	p02 = p01 +p01;
	p03 = p02 +p01;
	p04 = p03 +p01;
	p08 = p04 +p04;
	p0C = p08 +p04;
	p10 = p0C +p04;
	p14 = p10 +p04;
	p18 = p14 +p04;
	p1C = p18 +p04;

	p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
	p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
	p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
	p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
	p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
	p0C = p0C + ( (p0C >> DAT_BITS) << PAD_BITS );
	p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
	p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
	p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
	p1C = p1C + ( (p1C >> DAT_BITS) << PAD_BITS );
	ASSERT(HERE, p04+p04 == p08, "p04+p04 != p08");
	ASSERT(HERE, p04+p08 == p0C, "p04+p08 != p0C");

/*...The radix-32 pass is here.	*/

	iroot_prim=(incr >> 6);		/* (incr/2)/radix_now		*/

	for(m=0; m < nloops; m++)	/* NLOOPS may range from 1 (if first pass radix = 16) to P*N/32 (last pass radix = 16).x	*/
	{				/* NLOOPS satisfies the identity NLOOPS * INCR = P*N, which says that we process the entire
					   array each time this subroutine is executed (since P*N = vector length, sans paddring.)	*/

/*	here are the needed sincos data - these are processed below in bit-reversed order.	*/
	  iroot = index[m]*iroot_prim;
	  i = iroot;

/* In the DIF pass we may be able to afford a more-accurate twiddles computation,
   since we only compute each set of twiddles once and re-use it for many data blocks. */
#if HIACC
	#ifdef USE_SSE2
		/* Due to roots-locality considerations, roots (c,s)[0-31] are offset w.r.to the thread-local ptr pair as
					cc[00 01 02 03 04 05 06 07 08 09 0A 0B 0C 0D 0E 0F 10 11 12 13 14 15 16 17 18 19 1A 1B 1C 1D 1E 1F]
		(cc0,ss0) + 0x[06,26,16,36|0e,2e,1e,3e|0a,2a,1a,3a|12,32,22,42|08,28,18,38|10,30,20,40|0c,2c,1c,3c|14,34,24,44].
		*/
		c_tmp = cc0 + 0x06; s_tmp = c_tmp+1;	/* c0,s0 */
		c_tmp->re = c_tmp->im = 1.0;
		s_tmp->re = s_tmp->im = 0.0;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x26; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c01=rt;		s01=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x16; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c02=rt;		s02=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x36; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c03=rt;		s03=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x0e; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c04=rt;		s04=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x2e; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c05=rt;		s05=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x1e; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c06=rt;		s06=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x3e; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c07=rt;		s07=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x0a; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c08=rt;		s08=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x2a; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c09=rt;		s09=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x1a; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c0A=rt;		s0A=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x3a; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c0B=rt;		s0B=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x12; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c0C=rt;		s0C=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x32; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c0D=rt;		s0D=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x22; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c0E=rt;		s0E=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x42; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c0F=rt;		s0F=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x08; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c10=rt;		s10=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x28; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c11=rt;		s11=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x18; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c12=rt;		s12=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x38; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c13=rt;		s13=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x10; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c14=rt;		s14=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x30; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c15=rt;		s15=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x20; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c16=rt;		s16=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x40; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c17=rt;		s17=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x0c; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c18=rt;		s18=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x2c; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c19=rt;		s19=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x1c; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c1A=rt;		s1A=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x3c; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c1B=rt;		s1B=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x14; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c1C=rt;		s1C=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x34; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c1D=rt;		s1D=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x24; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c1E=rt;		s1E=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x44; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c1F=rt;		s1F=it;
	#endif

#else
	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;			/* 2*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c01=t00*rt -t01*it;	s01=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;				/* 3*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c02=t00*rt -t01*it;	s02=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += 4*iroot;				/* 7*iroot	*/
	    iroot = i;				/* 7*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c03=t00*rt -t01*it;	s03=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;				/* 14*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c07=t00*rt -t01*it;	s07=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;				/* 21*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c0E=t00*rt -t01*it;	s0E=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;				/* 28*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c15=t00*rt -t01*it;	s15=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c1C=t00*rt -t01*it;	s1C=t00*it +t01*rt;

/*(c,s)4-10 (decimal indices), 4-0A (hexadecimal indices):	*/
	    t00=c01*c07; t01=c01*s07; t02=s01*c07; t03=s01*s07;
	    c06=t00+t03; s06=t01-t02; c08=t00-t03; s08=t01+t02;

	    t00=c02*c07; t01=c02*s07; t02=s02*c07; t03=s02*s07;
	    c05=t00+t03; s05=t01-t02; c09=t00-t03; s09=t01+t02;

	    t00=c03*c07; t01=c03*s07; t02=s03*c07; t03=s03*s07;
	    c04=t00+t03; s04=t01-t02; c0A=t00-t03; s0A=t01+t02;

/*(c,s)11-17 (decimal indices), 0B-11 (hexadecimal indices):;	*/
	    t00=c01*c0E; t01=c01*s0E; t02=s01*c0E; t03=s01*s0E;
	    c0D=t00+t03; s0D=t01-t02; c0F=t00-t03; s0F=t01+t02;

	    t00=c02*c0E; t01=c02*s0E; t02=s02*c0E; t03=s02*s0E;
	    c0C=t00+t03; s0C=t01-t02; c10=t00-t03; s10=t01+t02;

	    t00=c03*c0E; t01=c03*s0E; t02=s03*c0E; t03=s03*s0E;
	    c0B=t00+t03; s0B=t01-t02; c11=t00-t03; s11=t01+t02;

/*(c,s)18-24 (decimal indices), 12-18 (hexadecimal indices):	*/
	    t00=c01*c15; t01=c01*s15; t02=s01*c15; t03=s01*s15;
	    c14=t00+t03; s14=t01-t02; c16=t00-t03; s16=t01+t02;

	    t00=c02*c15; t01=c02*s15; t02=s02*c15; t03=s02*s15;
	    c13=t00+t03; s13=t01-t02; c17=t00-t03; s17=t01+t02;

	    t00=c03*c15; t01=c03*s15; t02=s03*c15; t03=s03*s15;
	    c12=t00+t03; s12=t01-t02; c18=t00-t03; s18=t01+t02;

/*(c,s)25-31 (decimal indices), 19-1F (hexadecimal indices):	*/
	    t00=c01*c1C; t01=c01*s1C; t02=s01*c1C; t03=s01*s1C;
	    c1B=t00+t03; s1B=t01-t02; c1D=t00-t03; s1D=t01+t02;

	    t00=c02*c1C; t01=c02*s1C; t02=s02*c1C; t03=s02*s1C;
	    c1A=t00+t03; s1A=t01-t02; c1E=t00-t03; s1E=t01+t02;

	    t00=c03*c1C; t01=c03*s1C; t02=s03*c1C; t03=s03*s1C;
	    c19=t00+t03; s19=t01-t02; c1F=t00-t03; s1F=t01+t02;
#endif

	/* Define the inner-loop parameters in terms of the outer-loop ones to make OpenMP's job easier: */
	  jlo = m*incr;
	  jhi = jlo+(incr >> 5);

	#ifdef USE_SSE2
	  for(j=jlo; j < jhi; j += 4)
	  {
		/* In SSE2 mode, data are arranged in [re0,re1,im0,im1] quartets, not the usual [re0,im0],[re1,im1] pairs.
		Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
		but to point to the proper real datum, we need to bit-reverse bits <0:1> of j, i.e. [0,1,2,3] ==> [0,2,1,3].
		*/
		j1 = (j & mask01) + br4[j&3];
	#elif defined(USE_SSE2)	/* This allows us to use #if 0 above and disable sse2-based *computation*, while still using sse2-style data layout */
	  for(j=jlo; j < jhi; j += 2)
	  {
		j1 = (j & mask01) + br4[j&3];
	#else
	  for(j=jlo; j < jhi; j += 2)
	  {
		j1 =  j;
	#endif
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;

#ifdef USE_SSE2

	/*       gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...
			 We process the sincos data in bit-reversed order.	*/

	#if 0//DEBUG_SSE2
		if(s01->re != 0.0)
		{
			rng_isaac_init(TRUE);
			jt = j1;		jp = j2;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p04;	jp = j2 + p04;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p08;	jp = j2 + p08;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p0C;	jp = j2 + p0C;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p10;	jp = j2 + p10;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p14;	jp = j2 + p14;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p18;	jp = j2 + p18;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p1C;	jp = j2 + p1C;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();

			jt = j1; jp = j2;
			fprintf(stderr, "radix32_dif_pass: a:00 = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dif_pass: a:02 = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dif_pass: a:04 = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dif_pass: a:06 = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = j1 + p04;	jp = j2 + p04;
			fprintf(stderr, "radix32_dif_pass: a:08 = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dif_pass: a:0A = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dif_pass: a:0C = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dif_pass: a:0E = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = j1 + p08;	jp = j2 + p08;

			fprintf(stderr, "radix32_dif_pass: a:10 = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dif_pass: a:12 = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dif_pass: a:14 = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dif_pass: a:16 = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = j1 + p0C;	jp = j2 + p0C;
			fprintf(stderr, "radix32_dif_pass: a:18 = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dif_pass: a:1A = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dif_pass: a:1C = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dif_pass: a:1E = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = j1 + p10;	jp = j2 + p10;

			fprintf(stderr, "radix32_dif_pass: a:20 = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dif_pass: a:22 = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dif_pass: a:24 = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dif_pass: a:26 = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = j1 + p14;	jp = j2 + p14;
			fprintf(stderr, "radix32_dif_pass: a:28 = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dif_pass: a:2A = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dif_pass: a:2C = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dif_pass: a:2E = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = j1 + p18;	jp = j2 + p18;

			fprintf(stderr, "radix32_dif_pass: a:30 = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dif_pass: a:32 = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dif_pass: a:34 = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dif_pass: a:36 = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = j1 + p1C;	jp = j2 + p1C;
			fprintf(stderr, "radix32_dif_pass: a:38 = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dif_pass: a:3A = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dif_pass: a:3C = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dif_pass: a:0E = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);
		}
	#endif

	#ifdef COMPILER_TYPE_MSVC

	/*...Block 1: */
	#if 1	// if(1) - test out pure-asm version
		add0 = &a[j1];
		__asm	mov	eax, add0
		__asm	mov	ebx, p08		// Can't get these via simple load-one-and-shift-as-needed due to array padding scheme
		__asm	mov	ecx, p10
		__asm	mov	edx, p18
		__asm	shl	ebx,  3
		__asm	shl	ecx,  3
		__asm	shl	edx,  3
		__asm	add	ebx, eax
		__asm	add	ecx, eax
		__asm	add	edx, eax
		SSE2_RADIX4_DIF_4TWIDDLE_A         (r00,c10)	// Saves one CMUL compared to _B version [since first twiddle is unity]
//		SSE2_RADIX4_DIF_4TWIDDLE_B         (r00,c00)
		__asm	mov	edi, p04
		__asm	shl	edi,  3		// 8-bytes for array-of-doubles
		__asm	add	eax, edi	// &a[j1+p4];
		__asm	add	ebx, edi
		__asm	add	ecx, edi
		__asm	add	edx, edi
		SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r08,c04)
		/* Combine the 2 radix-4 subtransforms: */
		SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r00)
	#else
		add0 = &a[j1];
		add1 = add0+p04;
		add2 = add0+p08;
		add3 = add0+p0C;
		add4 = add0+p10;
		add5 = add0+p14;
		add6 = add0+p18;
		add7 = add0+p1C;

		SSE2_RADIX4_DIF_4TWIDDLE         (add0,add2,add4,add6,r00,c00)
		SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(add1,add3,add5,add7,r08,c04)
		/* Combine the 2 radix-4 subtransforms: */
		SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r00,r02,r04,r06,r08,r0A,r0C,r0E)
	#endif

	/*...Block 2:	*/
	#if 1	// if(1) - test out pure-asm version
		__asm	mov	edi, p02	// Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p+p1 not guaranteed.
		__asm	shl	edi,  3		// 8-bytes for array-of-doubles
		__asm	sub	eax, edi	// &a[j1+p2]
		__asm	sub	ebx, edi
		__asm	sub	ecx, edi
		__asm	sub	edx, edi
		SSE2_RADIX4_DIF_4TWIDDLE_B         (r10,c02)
		__asm	mov	edi, p04
		__asm	shl	edi,  3		// 8-bytes for array-of-doubles
		__asm	add	eax, edi	// &a[j1+p6];
		__asm	add	ebx, edi
		__asm	add	ecx, edi
		__asm	add	edx, edi
		SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r18,c06)
		/* Combine the 2 radix-4 subtransforms: */
		SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r10)
	#else
		add0 = &a[j1+p02];
		add1 = add0+p04;
		add2 = add0+p08;
		add3 = add0+p0C;
		add4 = add0+p10;
		add5 = add0+p14;
		add6 = add0+p18;
		add7 = add0+p1C;

		SSE2_RADIX4_DIF_4TWIDDLE         (add0,add2,add4,add6,r10,c02)
		SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(add1,add3,add5,add7,r18,c06)
		/* Combine the 2 radix-4 subtransforms: */
		SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r10,r12,r14,r16,r18,r1A,r1C,r1E)
	#endif

	/*...Block 3:	*/
	#if 1	// if(1) - test out pure-asm version
		__asm	mov	eax, add0
		__asm	mov	edi, p01	// Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p+p1 not guaranteed.
		__asm	mov	ebx, p08		// Can't get these via simple load-one-and-shift-as-needed due to array padding scheme
		__asm	mov	ecx, p10
		__asm	mov	edx, p18
		__asm	shl	edi,  3		// 8-bytes for array-of-doubles
		__asm	shl	ebx,  3
		__asm	shl	ecx,  3
		__asm	shl	edx,  3
		__asm	add	eax, edi	// &a[j1+p1];
		__asm	add	ebx, eax
		__asm	add	ecx, eax
		__asm	add	edx, eax
		SSE2_RADIX4_DIF_4TWIDDLE_B         (r20,c01)
		__asm	mov	edi, p04
		__asm	shl	edi,  3		// 8-bytes for array-of-doubles
		__asm	add	eax, edi	// &a[j1+p5];
		__asm	add	ebx, edi
		__asm	add	ecx, edi
		__asm	add	edx, edi
		SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r28,c05)
		/* Combine the 2 radix-4 subtransforms: */
		SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r20)
	#else
		add0 = &a[j1+p01];
		add1 = add0+p04;
		add2 = add0+p08;
		add3 = add0+p0C;
		add4 = add0+p10;
		add5 = add0+p14;
		add6 = add0+p18;
		add7 = add0+p1C;

		SSE2_RADIX4_DIF_4TWIDDLE         (add0,add2,add4,add6,r20,c01)
		SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(add1,add3,add5,add7,r28,c05)
		/* Combine the 2 radix-4 subtransforms: */
		SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r20,r22,r24,r26,r28,r2A,r2C,r2E)
	#endif

	/*...Block 4:	*/
	#if 1	// if(1) - test out pure-asm version
		__asm	mov	edi, p02	// Do things this way [rather than repeatedly adding p1] since array-padding scheme means p2 == p+p1 not guaranteed.
		__asm	shl	edi,  3		// 8-bytes for array-of-doubles
		__asm	sub	eax, edi	// &a[j1+p3]
		__asm	sub	ebx, edi
		__asm	sub	ecx, edi
		__asm	sub	edx, edi
		SSE2_RADIX4_DIF_4TWIDDLE_B         (r30,c03)
		__asm	mov	edi, p04
		__asm	shl	edi,  3		// 8-bytes for array-of-doubles
		__asm	add	eax, edi	// &a[j1+p7];
		__asm	add	ebx, edi
		__asm	add	ecx, edi
		__asm	add	edx, edi
		SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO_B(r38,c07)
		/* Combine the 2 radix-4 subtransforms: */
		SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS_B(r30)
	#else
		add0 = &a[jt];
		add1 = add0+p04;
		add2 = add0+p08;
		add3 = add0+p0C;
		add4 = add0+p10;
		add5 = add0+p14;
		add6 = add0+p18;
		add7 = add0+p1C;

		SSE2_RADIX4_DIF_4TWIDDLE         (add0,add2,add4,add6,r30,c03)
		SSE2_RADIX4_DIF_4TWIDDLE_2NDOFTWO(add1,add3,add5,add7,r38,c07)
		/* Combine the 2 radix-4 subtransforms: */
		SSE2_RADIX8_DIF_COMBINE_RAD4_SUBS(r30,r32,r34,r36,r38,r3A,r3C,r3E)
	#endif

	/*...and now do eight radix-4 transforms, including the internal twiddle factors:
		1, exp(i* 1*twopi/32) =       ( c32_1, s32_1), exp(i* 2*twopi/32) =       ( c    , s    ), exp(i* 3*twopi/32) =       ( c32_3, s32_3) (for inputs to transform block 2),
		1, exp(i* 2*twopi/32) =       ( c    , s    ), exp(i* 4*twopi/32) = ISRT2*( 1    , 1    ), exp(i* 6*twopi/32) =       ( s    , c    ) (for inputs to transform block 3),
		1, exp(i* 3*twopi/32) =       ( c32_3, s32_3), exp(i* 6*twopi/32) =       ( s    , c    ), exp(i* 9*twopi/32) =       (-s32_1, c32_1) (for inputs to transform block 4),
		1, exp(i* 4*twopi/32) = ISRT2*( 1    , 1    ), exp(i* 8*twopi/32) =       ( 0    , 1    ), exp(i*12*twopi/32) = ISRT2*(-1    , 1    ) (for inputs to transform block 5),
		1, exp(i* 5*twopi/32) =       ( s32_3, c32_3), exp(i*10*twopi/32) =       (-s    , c    ), exp(i*15*twopi/32) =       (-c32_1, s32_1) (for inputs to transform block 6),
		1, exp(i* 6*twopi/32) =       ( s    , c    ), exp(i*12*twopi/32) = ISRT2*(-1    , 1    ), exp(i*18*twopi/32) =       (-c    ,-s    ) (for inputs to transform block 7),
		1, exp(i* 7*twopi/32) =       ( s32_1, c32_1), exp(i*14*twopi/32) =       (-c    , s    ), exp(i*21*twopi/32) =       (-s32_3,-c32_3) (for inputs to transform block 8),
		 and only the last 3 inputs to each of the radix-4 transforms 2 through 8 are multiplied by non-unity twiddles.	*/

	/*...Block 1: t00,t10,t20,t30	*/
	#if 1
		__asm	mov eax, add0	// &a[j1]
		__asm	mov ebx, p01
		__asm	mov ecx, p02
		__asm	mov edx, p03
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	shl	ecx, 3
		__asm	shl	edx, 3
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
	#else
		add0 = &a[j1];
		add1 = add0+p01;
		add2 = add0+p02;
		add3 = add0+p03;
		__asm	mov	eax, add0	/* restore main-array indices */
		__asm	mov	ebx, add1
		__asm	mov	ecx, add2
		__asm	mov	edx, add3
	#endif

		__asm	mov	edi, r00

		__asm	movaps	xmm0,[edi      ]	/* t00 */				__asm	movaps	xmm4,[edi+0x200]	/* t20 */
		__asm	movaps	xmm1,[edi+0x010]	/* t01 */				__asm	movaps	xmm5,[edi+0x210]	/* t21 */
		__asm	movaps	xmm2,[edi+0x100]	/* t10 */				__asm	movaps	xmm6,[edi+0x300]	/* t30 */
		__asm	movaps	xmm3,[edi+0x110]	/* t11 */				__asm	movaps	xmm7,[edi+0x310]	/* t31 */

		__asm	subpd	xmm0,[edi+0x100]	/* t10=t00-rt */		__asm	subpd	xmm4,[edi+0x300]	/* t30=t20-rt */
		__asm	subpd	xmm1,[edi+0x110]	/* t11=t01-it */		__asm	subpd	xmm5,[edi+0x310]	/* t31=t21-it */
		__asm	addpd	xmm2,[edi      ]	/* t00=t00+rt */		__asm	addpd	xmm6,[edi+0x200]	/* t20=t20+rt */
		__asm	addpd	xmm3,[edi+0x010]	/* t01=t01+it */		__asm	addpd	xmm7,[edi+0x210]	/* t21=t21+it */

		__asm	subpd	xmm2,xmm6		/* t00 <- t00-t20 */		__asm	subpd	xmm0,xmm5		/* t10 <- t10-t31 */
		__asm	subpd	xmm3,xmm7		/* t01 <- t01-t21 */		__asm	subpd	xmm1,xmm4		/* t11 <- t11-t30 */
		__asm	addpd	xmm6,xmm6		/*          2*t20 */		__asm	addpd	xmm5,xmm5		/*          2*t31 */
		__asm	addpd	xmm7,xmm7		/*          2*t21 */		__asm	addpd	xmm4,xmm4		/*          2*t30 */
		__asm	movaps	[ebx     ],xmm2	/* a[jt+p1 ] */				__asm	movaps	[ecx     ],xmm0	/* a[jt+p2 ] */
		__asm	movaps	[ebx+0x10],xmm3	/* a[jp+p1 ] */				__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
		__asm	addpd	xmm6,xmm2		/* t20 <- t00+t20 */		__asm	addpd	xmm5,xmm0		/* t31 <- t10+t31 */
		__asm	addpd	xmm7,xmm3		/* t21 <- t01+t21 */		__asm	addpd	xmm4,xmm1		/* t30 <- t11+t30 */
		__asm	movaps	[eax     ],xmm6	/* a[jt+p0 ] */				__asm	movaps	[edx     ],xmm5	/* a[jt+p3 ] */
		__asm	movaps	[eax+0x10],xmm7	/* a[jp+p0 ] */				__asm	movaps	[ecx+0x10],xmm4	/* a[jp+p2 ] */

	/*...Block 5: t08,t18,t28,t38	*/
	#if 1
		__asm	sub ebx, eax
		__asm	sub ecx, eax
		__asm	sub edx, eax
		__asm	mov	edi, p04
		__asm	shl	edi, 3
		__asm	add eax, edi	// &a[j1+p04]
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
	#else
		add0 = &a[j1+p04];
		add1 = add0+p01;
		add2 = add0+p02;
		add3 = add0+p03;
	#endif
		__asm	mov	edi, r08
		__asm	mov	esi, isrt2
		__asm	movaps	xmm3,[esi ]	/* isrt2 */
																__asm	movaps	xmm4,[edi+0x200]	/* t28 */
																__asm	movaps	xmm5,[edi+0x210]	/* t29 */
																__asm	movaps	xmm6,[edi+0x300]	/* t38 */
																__asm	movaps	xmm7,[edi+0x310]	/* t39 */
																__asm	mulpd	xmm4,xmm3	/* t28 *ISRT2 */
		__asm	movaps	xmm0,[edi      ]	/* t08 */			__asm	mulpd	xmm5,xmm3	/* t29 *ISRT2 */
		__asm	movaps	xmm1,[edi+0x010]	/* t09 */			__asm	mulpd	xmm6,xmm3	/* t38 *ISRT2 */
		__asm	movaps	xmm2,[edi+0x100]	/* t18 */			__asm	mulpd	xmm7,xmm3	/* t39 *ISRT2 */
		__asm	movaps	xmm3,[edi+0x110]	/* t19; this must execute after the last mul-by-ISRT2 above */

		__asm	subpd	xmm0,xmm3		/* ~t08= t08-t19*/		__asm	subpd	xmm4,xmm5		/* ~t28=t28-t29*/
		__asm	subpd	xmm1,xmm2		/* ~t19= t09-t18*/		__asm	subpd	xmm7,xmm6		/*  it =t39-t38*/
		__asm	addpd	xmm3,xmm3		/*         2*t19*/		__asm	addpd	xmm5,xmm5		/*        2*t29*/
		__asm	addpd	xmm2,xmm2		/*         2*t18*/		__asm	addpd	xmm6,xmm6		/*        2*t38*/
		__asm	addpd	xmm3,xmm0		/* ~t18= t19+t08*/		__asm	addpd	xmm5,xmm4		/* ~t29=t29+t28*/
		__asm	addpd	xmm2,xmm1		/* ~t09= t18+t09*/		__asm	addpd	xmm6,xmm7		/*  rt =t38+t39*/

		__asm	subpd	xmm4,xmm6		/* t28=t28-rt */
		__asm	subpd	xmm5,xmm7		/* t29=t29-it */
		__asm	addpd	xmm6,xmm6		/*      2* rt */
		__asm	addpd	xmm7,xmm7		/*      2* it */
		__asm	addpd	xmm6,xmm4		/* t38=t28+rt */
		__asm	addpd	xmm7,xmm5		/* t39=t29+it */

		__asm	subpd	xmm0,xmm4		/* t08-t28 */			__asm	subpd	xmm3,xmm7		/* t18-t39 */
		__asm	subpd	xmm2,xmm5		/* t09-t29 */			__asm	subpd	xmm1,xmm6		/* t19-t38 */
		__asm	addpd	xmm4,xmm4		/*   2*t28 */			__asm	addpd	xmm7,xmm7		/*   2*t39 */
		__asm	addpd	xmm5,xmm5		/*   2*t29 */			__asm	addpd	xmm6,xmm6		/*   2*t38 */
		__asm	movaps	[ebx     ],xmm0	/* a[jt+p1 ] */			__asm	movaps	[ecx     ],xmm3	/* a[jt+p2 ] */
		__asm	movaps	[ebx+0x10],xmm2	/* a[jp+p1 ] */			__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
		__asm	addpd	xmm4,xmm0		/* t08+t28 */			__asm	addpd	xmm7,xmm3		/* t18+t39 */
		__asm	addpd	xmm5,xmm2		/* t09+t29 */			__asm	addpd	xmm6,xmm1		/* t19+t38 */
		__asm	movaps	[eax     ],xmm4	/* a[jt+p0 ] */			__asm	movaps	[edx     ],xmm7	/* a[jt+p3 ] */
		__asm	movaps	[eax+0x10],xmm5	/* a[jp+p0 ] */			__asm	movaps	[ecx+0x10],xmm6	/* a[jp+p2 ] */

	/*...Block 3: t04,t14,t24,t34	*/
	#if 1
		__asm	sub ebx, eax
		__asm	sub ecx, eax
		__asm	sub edx, eax
		__asm	mov	edi, p04
		__asm	shl	edi, 3
		__asm	sub eax, edi	// &a[j1]
		__asm	mov	edi, p08
		__asm	shl	edi, 3
		__asm	add eax, edi	// &a[j1+p08]
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
	#else
		add0 = &a[j1+p08];
		add1 = add0+p01;
		add2 = add0+p02;
		add3 = add0+p03;
		__asm	mov	eax, add0	/* restore main-array indices */
		__asm	mov	ebx, add1
		__asm	mov	ecx, add2
		__asm	mov	edx, add3
	#endif
		__asm	mov	edi, r04
		__asm	mov	esi, cc0
		__asm	movaps	xmm4,[edi+0x200]	/* t24 */
		__asm	movaps	xmm5,[edi+0x210]	/* t25 */
		__asm	movaps	xmm3,[esi      ]	/* c */
		__asm	movaps	xmm2,[esi+0x010]	/* s */
		__asm	movaps	xmm6,xmm4		/* copy t24 */
		__asm	movaps	xmm7,xmm5		/* copy t25 */

		__asm	mulpd	xmm4,xmm3		/* t24*c */
		__asm	mulpd	xmm5,xmm3		/* t25*c */
		__asm	mulpd	xmm6,xmm2		/* t24*s */				__asm	movaps	xmm0,[edi+0x300]	/* t34 */
		__asm	mulpd	xmm7,xmm2		/* t25*s */				__asm	movaps	xmm1,[edi+0x310]	/* t35 */
		__asm	addpd	xmm5,xmm6	/* ~t25 */					__asm	movaps	xmm6,xmm0		/* copy t34 */
		__asm	subpd	xmm4,xmm7	/* ~t24 */					__asm	movaps	xmm7,xmm1		/* copy t35 */

																__asm	mulpd	xmm6,xmm2		/* t34*s */
																__asm	mulpd	xmm7,xmm2		/* t35*s */
																__asm	mulpd	xmm0,xmm3		/* t34*c */
																__asm	mulpd	xmm1,xmm3		/* t35*c */
																__asm	addpd	xmm7,xmm0	/* it */
																__asm	subpd	xmm6,xmm1	/* rt */

		__asm	movaps	xmm2,xmm4		/* copy t24 */
		__asm	movaps	xmm3,xmm5		/* copy t25 */
		__asm	subpd	xmm4,xmm6		/*~t34=t24-rt */
		__asm	subpd	xmm5,xmm7		/*~t35=t25-it */
		__asm	addpd	xmm6,xmm2		/*~t24=t24+rt */
		__asm	addpd	xmm7,xmm3		/*~t25=t25+it */

		__asm	mov	esi, isrt2
		__asm	movaps	xmm2,[edi+0x100]	/* t14 */
		__asm	movaps	xmm3,[edi+0x110]	/* t15 */
		__asm	movaps	xmm1,[esi]	/* isrt2 */
		__asm	movaps	xmm0,xmm2	/* cpy t14 */
		__asm	subpd	xmm2,xmm3	/*~t14=t14-t15 */
		__asm	addpd	xmm3,xmm0	/*~t15=t15+t14 */
		__asm	mulpd	xmm2,xmm1	/* rt = (t14-t15)*ISRT2 */
		__asm	mulpd	xmm3,xmm1	/* it = (t15+t14)*ISRT2 */

		__asm	movaps	xmm0,[edi      ]	/* t04 */
		__asm	movaps	xmm1,[edi+0x010]	/* t05 */

		__asm	subpd	xmm0,xmm2			/*~t14=t04-rt */
		__asm	subpd	xmm1,xmm3			/*~t15=t05-it */
		__asm	addpd	xmm2,xmm2			/*      2* rt */
		__asm	addpd	xmm3,xmm3			/*      2* it */
		__asm	addpd	xmm2,xmm0			/*~t04=rt +t04*/
		__asm	addpd	xmm3,xmm1			/*~t05=it +t05*/

		__asm	subpd	xmm2,xmm6		/* t04-t24 */			__asm	subpd	xmm0,xmm5		/* t14-t35 */
		__asm	subpd	xmm3,xmm7		/* t05-t25 */			__asm	subpd	xmm1,xmm4		/* t15-t34 */
		__asm	addpd	xmm6,xmm6		/*   2*t24 */			__asm	addpd	xmm5,xmm5		/*          2*t35 */
		__asm	addpd	xmm7,xmm7		/*   2*t25 */			__asm	addpd	xmm4,xmm4		/*          2*t34 */
		__asm	movaps	[ebx     ],xmm2	/* a[jt+p1 ] */			__asm	movaps	[ecx     ],xmm0	/* a[jt+p2 ] */
		__asm	movaps	[ebx+0x10],xmm3	/* a[jp+p1 ] */			__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
		__asm	addpd	xmm6,xmm2		/* t04+t24 */			__asm	addpd	xmm5,xmm0		/* t14+t35 */
		__asm	addpd	xmm7,xmm3		/* t05+t25 */			__asm	addpd	xmm4,xmm1		/* t15+t34 */
		__asm	movaps	[eax     ],xmm6	/* a[jt+p0 ] */			__asm	movaps	[edx     ],xmm5	/* a[jt+p3 ] */
		__asm	movaps	[eax+0x10],xmm7	/* a[jp+p0 ] */			__asm	movaps	[ecx+0x10],xmm4	/* a[jp+p2 ] */

	/*...Block 7: t0C,t1C,t2C,t3C	*/
	#if 1
		__asm	sub ebx, eax
		__asm	sub ecx, eax
		__asm	sub edx, eax
		__asm	mov	edi, p08
		__asm	shl	edi, 3
		__asm	sub eax, edi	// &a[j1]
		__asm	mov	edi, p0C
		__asm	shl	edi, 3
		__asm	add eax, edi	// &a[j1+p0C]
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
	#else
		add0 = &a[j1+p0C];
		add1 = add0+p01;
		add2 = add0+p02;
		add3 = add0+p03;
		__asm	mov	eax, add0	/* restore main-array indices */
		__asm	mov	ebx, add1
		__asm	mov	ecx, add2
		__asm	mov	edx, add3
	#endif
		__asm	mov	edi, r0C
		__asm	mov	esi, cc0
		__asm	movaps	xmm4,[edi+0x200]	/* t2C */
		__asm	movaps	xmm5,[edi+0x210]	/* t2D */
		__asm	movaps	xmm2,[esi      ]	/* c */
		__asm	movaps	xmm3,[esi+0x010]	/* s */
		__asm	movaps	xmm6,xmm4		/* copy t2C */
		__asm	movaps	xmm7,xmm5		/* copy t2D */

		__asm	mulpd	xmm4,xmm3		/* t2C*s */
		__asm	mulpd	xmm5,xmm3		/* t2D*s */
		__asm	mulpd	xmm6,xmm2		/* t2C*c */				__asm	movaps	xmm0,[edi+0x300]	/* t3C */
		__asm	mulpd	xmm7,xmm2		/* t2D*c */				__asm	movaps	xmm1,[edi+0x310]	/* t3D */
		__asm	addpd	xmm5,xmm6	/* ~t2D */					__asm	movaps	xmm6,xmm0		/* copy t3C */
		__asm	subpd	xmm4,xmm7	/* ~t2C */					__asm	movaps	xmm7,xmm1		/* copy t3D */

																__asm	mulpd	xmm6,xmm2		/* t3C*c */
																__asm	mulpd	xmm7,xmm2		/* t3D*c */
																__asm	mulpd	xmm0,xmm3		/* t3C*s */
																__asm	mulpd	xmm1,xmm3		/* t3D*s */
																__asm	addpd	xmm7,xmm0	/* it */
																__asm	subpd	xmm6,xmm1	/* rt */

		__asm	movaps	xmm2,xmm4		/* copy t2C */
		__asm	movaps	xmm3,xmm5		/* copy t2D */
		__asm	subpd	xmm4,xmm6		/*~t2C=t2C-rt */
		__asm	subpd	xmm5,xmm7		/*~t2D=t2D-it */
		__asm	addpd	xmm6,xmm2		/*~t3C=t2C+rt */
		__asm	addpd	xmm7,xmm3		/*~t3D=t2D+it */

		__asm	mov	esi, isrt2
		__asm	movaps	xmm2,[edi+0x100]	/* t1C */
		__asm	movaps	xmm3,[edi+0x110]	/* t1D */
		__asm	movaps	xmm1,[esi]	/* isrt2 */
		__asm	movaps	xmm0,xmm2	/* cpy t1C */
		__asm	addpd	xmm2,xmm3	/*~t1C=t1C+t1D */
		__asm	subpd	xmm3,xmm0	/*~t1D=t1D-t1C */
		__asm	mulpd	xmm2,xmm1	/* rt = (t1C+t1D)*ISRT2 */
		__asm	mulpd	xmm3,xmm1	/* it = (t1D-t1C)*ISRT2 */

		__asm	movaps	xmm0,[edi      ]	/* t0C */
		__asm	movaps	xmm1,[edi+0x010]	/* t0D */

		__asm	subpd	xmm0,xmm2			/*~t0C=t0C-rt */
		__asm	subpd	xmm1,xmm3			/*~t0D=t0D-it */
		__asm	addpd	xmm2,xmm2			/*      2* rt */
		__asm	addpd	xmm3,xmm3			/*      2* it */
		__asm	addpd	xmm2,xmm0			/*~t1C=rt +t0C*/
		__asm	addpd	xmm3,xmm1			/*~t1D=it +t0D*/

		__asm	subpd	xmm0,xmm4		/* t0C-t2C */			__asm	subpd	xmm2,xmm7		/* t1C-t3D */
		__asm	subpd	xmm1,xmm5		/* t0D-t2D */			__asm	subpd	xmm3,xmm6		/* t1D-t3C */
		__asm	addpd	xmm4,xmm4		/*   2*t2C */			__asm	addpd	xmm7,xmm7		/*   2*t3D */
		__asm	addpd	xmm5,xmm5		/*   2*t2D */			__asm	addpd	xmm6,xmm6		/*   2*t3C */
		__asm	movaps	[ebx     ],xmm0	/* a[jt+p1 ] */			__asm	movaps	[ecx     ],xmm2	/* a[jt+p2 ] */
		__asm	movaps	[ebx+0x10],xmm1	/* a[jp+p1 ] */			__asm	movaps	[edx+0x10],xmm3	/* a[jp+p3 ] */
		__asm	addpd	xmm4,xmm0		/* t0C+t2C */			__asm	addpd	xmm7,xmm2		/* t1C+t3D */
		__asm	addpd	xmm5,xmm1		/* t0D+t2D */			__asm	addpd	xmm6,xmm3		/* t1D+t3C */
		__asm	movaps	[eax     ],xmm4	/* a[jt+p0 ] */			__asm	movaps	[edx     ],xmm7	/* a[jt+p3 ] */
		__asm	movaps	[eax+0x10],xmm5	/* a[jp+p0 ] */			__asm	movaps	[ecx+0x10],xmm6	/* a[jp+p2 ] */

	/*...Block 2: t02,t12,t22,t32	*/
		__asm	sub ebx, eax
		__asm	sub ecx, eax
		__asm	sub edx, eax
		__asm	mov	edi, p0C
		__asm	shl	edi, 3
		__asm	sub eax, edi	// &a[j1]
		__asm	mov	edi, p10
		__asm	shl	edi, 3
		__asm	add eax, edi	// &a[j1+p10]
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax

		__asm	mov	edi, r02
		__asm	mov	esi, cc1
		__asm	movaps	xmm4,[edi+0x200]	/* t22 */
		__asm	movaps	xmm5,[edi+0x210]	/* t23 */
		__asm	movaps	xmm2,[esi      ]	/* c32_1 */
		__asm	movaps	xmm3,[esi+0x010]	/* s32_1 */
		__asm	movaps	xmm6,xmm4		/* copy t22 */
		__asm	movaps	xmm7,xmm5		/* copy t23 */

																__asm	add	esi, 0x20	/* cc3 */
		__asm	mulpd	xmm4,xmm2		/* t22*c32_1 */
		__asm	mulpd	xmm5,xmm2		/* t23*c32_1 */
		__asm	mulpd	xmm6,xmm3		/* t22*s32_1 */			__asm	movaps	xmm0,[edi+0x300]	/* t32 */
		__asm	mulpd	xmm7,xmm3		/* t23*s32_1 */			__asm	movaps	xmm1,[edi+0x310]	/* t33 */
																__asm	movaps	xmm2,[esi      ]	/* c32_3 */
																__asm	movaps	xmm3,[esi+0x010]	/* s32_3 */
		__asm	addpd	xmm5,xmm6	/* ~t23 */					__asm	movaps	xmm6,xmm0		/* copy t32 */
		__asm	subpd	xmm4,xmm7	/* ~t22 */					__asm	movaps	xmm7,xmm1		/* copy t33 */

																__asm	mulpd	xmm6,xmm2		/* t32*c32_3 */
																__asm	mulpd	xmm7,xmm2		/* t33*c32_3 */
																__asm	mulpd	xmm0,xmm3		/* t32*s32_3 */
																__asm	mulpd	xmm1,xmm3		/* t33*s32_3 */
																__asm	addpd	xmm7,xmm0	/* it */
																__asm	subpd	xmm6,xmm1	/* rt */

		__asm	movaps	xmm2,xmm4		/* copy t22 */
		__asm	movaps	xmm3,xmm5		/* copy t23 */
		__asm	subpd	xmm4,xmm6		/*~t32=t22-rt */
		__asm	subpd	xmm5,xmm7		/*~t33=t23-it */
		__asm	addpd	xmm6,xmm2		/*~t22=t22+rt */
		__asm	addpd	xmm7,xmm3		/*~t23=t23+it */

		__asm	sub	esi, 0x40	/* cc0 */
		__asm	movaps	xmm1,[edi+0x100]	/* t12 */
		__asm	movaps	xmm3,[edi+0x110]	/* t13 */
		__asm	movaps	xmm0,[esi+0x010]	/* s */
		__asm	movaps	xmm2,xmm1			/* cpy t12 */
		__asm	mulpd	xmm1,xmm0		/* t12*s */
		__asm	mulpd	xmm0,xmm3		/* s*t13 */
		__asm	mulpd	xmm2,[esi]		/* t12*c */
		__asm	mulpd	xmm3,[esi]		/* t13*c */
		__asm	subpd	xmm2,xmm0	/* rt =t12*c - t13*s */
		__asm	addpd	xmm3,xmm1	/* it =t13*c + t12*s */

		__asm	movaps	xmm0,[edi      ]	/* t02 */
		__asm	movaps	xmm1,[edi+0x010]	/* t03 */

		__asm	subpd	xmm0,xmm2			/*~t12=t02-rt */
		__asm	subpd	xmm1,xmm3			/*~t13=t03-it */
		__asm	addpd	xmm2,xmm2			/*      2* rt */
		__asm	addpd	xmm3,xmm3			/*      2* it */
		__asm	addpd	xmm2,xmm0			/*~t02=rt +t02*/
		__asm	addpd	xmm3,xmm1			/*~t03=it +t03*/

		__asm	subpd	xmm2,xmm6		/* t02-t22 */			__asm	subpd	xmm0,xmm5		/* t12-t33 */
		__asm	subpd	xmm3,xmm7		/* t03-t23 */			__asm	subpd	xmm1,xmm4		/* t13-t32 */
		__asm	addpd	xmm6,xmm6		/*   2*t22 */			__asm	addpd	xmm5,xmm5		/*   2*t33 */
		__asm	addpd	xmm7,xmm7		/*   2*t23 */			__asm	addpd	xmm4,xmm4		/*   2*t32 */
		__asm	movaps	[ebx     ],xmm2	/* a[jt+p1 ] */			__asm	movaps	[ecx     ],xmm0	/* a[jt+p2 ] */
		__asm	movaps	[ebx+0x10],xmm3	/* a[jp+p1 ] */			__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
		__asm	addpd	xmm6,xmm2		/* t02+t22 */			__asm	addpd	xmm5,xmm0		/* t12+t33 */
		__asm	addpd	xmm7,xmm3		/* t03+t23 */			__asm	addpd	xmm4,xmm1		/* t13+t32 */
		__asm	movaps	[eax     ],xmm6	/* a[jt+p0 ] */			__asm	movaps	[edx     ],xmm5	/* a[jt+p3 ] */
		__asm	movaps	[eax+0x10],xmm7	/* a[jp+p0 ] */			__asm	movaps	[ecx+0x10],xmm4	/* a[jp+p2 ] */

	/*...Block 6: t0A,t1A,t2A,t3A	*/
		__asm	sub ebx, eax
		__asm	sub ecx, eax
		__asm	sub edx, eax
		__asm	mov	edi, p10
		__asm	shl	edi, 3
		__asm	sub eax, edi	// &a[j1]
		__asm	mov	edi, p14
		__asm	shl	edi, 3
		__asm	add eax, edi	// &a[j1+p14]
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax

		__asm	mov	edi, r0A
		__asm	mov	esi, cc3
		__asm	movaps	xmm4,[edi+0x200]	/* t22 */
		__asm	movaps	xmm5,[edi+0x210]	/* t23 */
		__asm	movaps	xmm3,[esi      ]	/* c32_3 */
		__asm	movaps	xmm2,[esi+0x010]	/* s32_3 */
		__asm	movaps	xmm6,xmm4		/* copy t22 */
		__asm	movaps	xmm7,xmm5		/* copy t23 */

																__asm	sub	esi, 0x20	/* cc1 */
		__asm	mulpd	xmm4,xmm2		/* t22*s32_3 */
		__asm	mulpd	xmm5,xmm2		/* t23*s32_3 */
		__asm	mulpd	xmm6,xmm3		/* t22*c32_3 */			__asm	movaps	xmm0,[edi+0x300]	/* t32 */
		__asm	mulpd	xmm7,xmm3		/* t23*c32_3 */			__asm	movaps	xmm1,[edi+0x310]	/* t33 */
																__asm	movaps	xmm2,[esi      ]	/* c32_1 */
																__asm	movaps	xmm3,[esi+0x010]	/* s32_1 */
		__asm	addpd	xmm5,xmm6	/* ~t23 */					__asm	movaps	xmm6,xmm0		/* copy t32 */
		__asm	subpd	xmm4,xmm7	/* ~t22 */					__asm	movaps	xmm7,xmm1		/* copy t33 */

																__asm	mulpd	xmm6,xmm2		/* t32*c32_1 */
																__asm	mulpd	xmm7,xmm2		/* t33*c32_1 */
																__asm	mulpd	xmm0,xmm3		/* t32*s32_1 */
																__asm	mulpd	xmm1,xmm3		/* t33*s32_1 */
																__asm	subpd	xmm7,xmm0	/* it */
																__asm	addpd	xmm6,xmm1	/* rt */

		__asm	movaps	xmm2,xmm4		/* copy t22 */
		__asm	movaps	xmm3,xmm5		/* copy t23 */
		__asm	subpd	xmm4,xmm6		/*~t22=t22-rt */
		__asm	subpd	xmm5,xmm7		/*~t23=t23-it */
		__asm	addpd	xmm6,xmm2		/*~t32=t22+rt */
		__asm	addpd	xmm7,xmm3		/*~t33=t23+it */

		__asm	sub	esi, 0x20	/* cc0 */
		__asm	movaps	xmm1,[edi+0x100]	/* t12 */
		__asm	movaps	xmm3,[edi+0x110]	/* t13 */
		__asm	movaps	xmm0,[esi]		/* c */
		__asm	movaps	xmm2,xmm1			/* cpy t12 */
		__asm	mulpd	xmm1,xmm0		/* t12*c */
		__asm	mulpd	xmm0,xmm3		/* c*t13 */
		__asm	mulpd	xmm2,[esi+0x010]/* t12*s */
		__asm	mulpd	xmm3,[esi+0x010]/* t13*s */
		__asm	addpd	xmm2,xmm0	/* rt =t12*s + t13*c */
		__asm	subpd	xmm3,xmm1	/* it =t13*s - t12*c */

		__asm	movaps	xmm0,[edi      ]	/* t02 */
		__asm	movaps	xmm1,[edi+0x010]	/* t03 */

		__asm	subpd	xmm0,xmm2			/*~t02=t02-rt */
		__asm	subpd	xmm1,xmm3			/*~t03=t03-it */
		__asm	addpd	xmm2,xmm2			/*      2* rt */
		__asm	addpd	xmm3,xmm3			/*      2* it */
		__asm	addpd	xmm2,xmm0			/*~t12=rt +t02*/
		__asm	addpd	xmm3,xmm1			/*~t13=it +t03*/

		__asm	subpd	xmm0,xmm4		/* t02-t22 */			__asm	subpd	xmm2,xmm7		/* t12-t33 */
		__asm	subpd	xmm1,xmm5		/* t03-t23 */			__asm	subpd	xmm3,xmm6		/* t13-t32 */
		__asm	addpd	xmm4,xmm4		/*   2*t22 */			__asm	addpd	xmm7,xmm7		/*   2*t33 */
		__asm	addpd	xmm5,xmm5		/*   2*t23 */			__asm	addpd	xmm6,xmm6		/*   2*t32 */
		__asm	movaps	[ebx     ],xmm0	/* a[jt+p1 ] */			__asm	movaps	[ecx     ],xmm2	/* a[jt+p2 ] */
		__asm	movaps	[ebx+0x10],xmm1	/* a[jp+p1 ] */			__asm	movaps	[edx+0x10],xmm3	/* a[jp+p3 ] */
		__asm	addpd	xmm4,xmm0		/* t02+t22 */			__asm	addpd	xmm7,xmm2		/* t12+t33 */
		__asm	addpd	xmm5,xmm1		/* t03+t23 */			__asm	addpd	xmm6,xmm3		/* t13+t32 */
		__asm	movaps	[eax     ],xmm4	/* a[jt+p0 ] */			__asm	movaps	[edx     ],xmm7	/* a[jt+p3 ] */
		__asm	movaps	[eax+0x10],xmm5	/* a[jp+p0 ] */			__asm	movaps	[ecx+0x10],xmm6	/* a[jp+p2 ] */

	/*...Block 4: t06,t16,t26,t36	*/
		__asm	sub ebx, eax
		__asm	sub ecx, eax
		__asm	sub edx, eax
		__asm	mov	edi, p14
		__asm	shl	edi, 3
		__asm	sub eax, edi	// &a[j1]
		__asm	mov	edi, p18
		__asm	shl	edi, 3
		__asm	add eax, edi	// &a[j1+p18]
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax

		__asm	mov	edi, r06
		__asm	mov	esi, cc3
		__asm	movaps	xmm4,[edi+0x200]	/* t22 */
		__asm	movaps	xmm5,[edi+0x210]	/* t23 */
		__asm	movaps	xmm2,[esi      ]	/* c32_3 */
		__asm	movaps	xmm3,[esi+0x010]	/* s32_3 */
		__asm	movaps	xmm6,xmm4		/* copy t22 */
		__asm	movaps	xmm7,xmm5		/* copy t23 */

																__asm	sub	esi, 0x20	/* cc1 */
		__asm	mulpd	xmm4,xmm2		/* t22*c32_3 */
		__asm	mulpd	xmm5,xmm2		/* t23*c32_3 */
		__asm	mulpd	xmm6,xmm3		/* t22*s32_3 */			__asm	movaps	xmm0,[edi+0x300]	/* t32 */
		__asm	mulpd	xmm7,xmm3		/* t23*s32_3 */			__asm	movaps	xmm1,[edi+0x310]	/* t33 */
																__asm	movaps	xmm3,[esi      ]	/* c32_1 */
																__asm	movaps	xmm2,[esi+0x010]	/* s32_1 */
		__asm	addpd	xmm5,xmm6	/* ~t23 */					__asm	movaps	xmm6,xmm0		/* copy t32 */
		__asm	subpd	xmm4,xmm7	/* ~t22 */					__asm	movaps	xmm7,xmm1		/* copy t33 */

																__asm	mulpd	xmm6,xmm2		/* t32*s32_1 */
																__asm	mulpd	xmm7,xmm2		/* t33*s32_1 */
																__asm	mulpd	xmm0,xmm3		/* t32*c32_1 */
																__asm	mulpd	xmm1,xmm3		/* t33*c32_1 */
																__asm	subpd	xmm7,xmm0	/* it */
																__asm	addpd	xmm6,xmm1	/* rt */

		__asm	movaps	xmm2,xmm4		/* copy t22 */
		__asm	movaps	xmm3,xmm5		/* copy t23 */
		__asm	subpd	xmm4,xmm6		/*~t22=t22-rt */
		__asm	subpd	xmm5,xmm7		/*~t23=t23-it */
		__asm	addpd	xmm6,xmm2		/*~t32=t22+rt */
		__asm	addpd	xmm7,xmm3		/*~t33=t23+it */

		__asm	sub	esi, 0x20	/* cc0 */
		__asm	movaps	xmm1,[edi+0x100]	/* t12 */
		__asm	movaps	xmm3,[edi+0x110]	/* t13 */
		__asm	movaps	xmm0,[esi]		/* c */
		__asm	movaps	xmm2,xmm1			/* cpy t12 */
		__asm	mulpd	xmm1,xmm0		/* t12*c */
		__asm	mulpd	xmm0,xmm3		/* c*t13 */
		__asm	mulpd	xmm2,[esi+0x010]/* t12*s */
		__asm	mulpd	xmm3,[esi+0x010]/* t13*s */
		__asm	subpd	xmm2,xmm0	/* rt =t12*s - t13*c */
		__asm	addpd	xmm3,xmm1	/* it =t13*s + t12*c */

		__asm	movaps	xmm0,[edi      ]	/* t02 */
		__asm	movaps	xmm1,[edi+0x010]	/* t03 */

		__asm	subpd	xmm0,xmm2			/*~t12=t02-rt */
		__asm	subpd	xmm1,xmm3			/*~t13=t03-it */
		__asm	addpd	xmm2,xmm2			/*      2* rt */
		__asm	addpd	xmm3,xmm3			/*      2* it */
		__asm	addpd	xmm2,xmm0			/*~t02=rt +t02*/
		__asm	addpd	xmm3,xmm1			/*~t03=it +t03*/

		__asm	subpd	xmm2,xmm4		/* t02-t22 */			__asm	subpd	xmm0,xmm7		/* t12-t33 */
		__asm	subpd	xmm3,xmm5		/* t03-t23 */			__asm	subpd	xmm1,xmm6		/* t13-t32 */
		__asm	addpd	xmm4,xmm4		/*   2*t22 */			__asm	addpd	xmm7,xmm7		/*   2*t33 */
		__asm	addpd	xmm5,xmm5		/*   2*t23 */			__asm	addpd	xmm6,xmm6		/*   2*t32 */
		__asm	movaps	[ebx     ],xmm2	/* a[jt+p1 ] */			__asm	movaps	[ecx     ],xmm0	/* a[jt+p2 ] */
		__asm	movaps	[ebx+0x10],xmm3	/* a[jp+p1 ] */			__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
		__asm	addpd	xmm4,xmm2		/* t02+t22 */			__asm	addpd	xmm7,xmm0		/* t12+t33 */
		__asm	addpd	xmm5,xmm3		/* t03+t23 */			__asm	addpd	xmm6,xmm1		/* t13+t32 */
		__asm	movaps	[eax     ],xmm4	/* a[jt+p0 ] */			__asm	movaps	[edx     ],xmm7	/* a[jt+p3 ] */
		__asm	movaps	[eax+0x10],xmm5	/* a[jp+p0 ] */			__asm	movaps	[ecx+0x10],xmm6	/* a[jp+p2 ] */

	/*...Block 8: t0E,t1E,t2E,t3E	*/
		__asm	sub ebx, eax
		__asm	sub ecx, eax
		__asm	sub edx, eax
		__asm	mov	edi, p18
		__asm	shl	edi, 3
		__asm	sub eax, edi	// &a[j1]
		__asm	mov	edi, p1C
		__asm	shl	edi, 3
		__asm	add eax, edi	// &a[j1+p1C]
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax

		__asm	mov	edi, r0E
		__asm	mov	esi, cc1
		__asm	movaps	xmm4,[edi+0x200]	/* t22 */
		__asm	movaps	xmm5,[edi+0x210]	/* t23 */
		__asm	movaps	xmm3,[esi      ]	/* c32_1 */
		__asm	movaps	xmm2,[esi+0x010]	/* s32_1 */
		__asm	movaps	xmm6,xmm4		/* copy t22 */
		__asm	movaps	xmm7,xmm5		/* copy t23 */

																__asm	add	esi, 0x20	/* cc3 */
		__asm	mulpd	xmm4,xmm2		/* t22*s32_1 */
		__asm	mulpd	xmm5,xmm2		/* t23*s32_1 */
		__asm	mulpd	xmm6,xmm3		/* t22*c32_1 */			__asm	movaps	xmm0,[edi+0x300]	/* t32 */
		__asm	mulpd	xmm7,xmm3		/* t23*c32_1 */			__asm	movaps	xmm1,[edi+0x310]	/* t33 */
																__asm	movaps	xmm3,[esi      ]	/* c32_3 */
																__asm	movaps	xmm2,[esi+0x010]	/* s32_3 */
		__asm	addpd	xmm5,xmm6	/* ~t23 */					__asm	movaps	xmm6,xmm0		/* copy t32 */
		__asm	subpd	xmm4,xmm7	/* ~t22 */					__asm	movaps	xmm7,xmm1		/* copy t33 */

																__asm	mulpd	xmm6,xmm2		/* t32*s32_3 */
																__asm	mulpd	xmm7,xmm2		/* t33*s32_3 */
																__asm	mulpd	xmm0,xmm3		/* t32*c32_3 */
																__asm	mulpd	xmm1,xmm3		/* t33*c32_3 */
																__asm	addpd	xmm7,xmm0	/* it */
																__asm	subpd	xmm6,xmm1	/* rt */

		__asm	movaps	xmm2,xmm4		/* copy t22 */
		__asm	movaps	xmm3,xmm5		/* copy t23 */
		__asm	subpd	xmm4,xmm6		/*~t22=t22-rt */
		__asm	subpd	xmm5,xmm7		/*~t23=t23-it */
		__asm	addpd	xmm6,xmm2		/*~t32=t22+rt */
		__asm	addpd	xmm7,xmm3		/*~t33=t23+it */

		__asm	sub	esi, 0x40	/* cc0 */
		__asm	movaps	xmm1,[edi+0x100]	/* t12 */
		__asm	movaps	xmm3,[edi+0x110]	/* t13 */
		__asm	movaps	xmm0,[esi+0x010]	/* s */
		__asm	movaps	xmm2,xmm1			/* cpy t12 */
		__asm	mulpd	xmm1,xmm0		/* t12*s */
		__asm	mulpd	xmm0,xmm3		/* s*t13 */
		__asm	mulpd	xmm2,[esi]		/* t12*c */
		__asm	mulpd	xmm3,[esi]		/* t13*c */
		__asm	addpd	xmm2,xmm0	/* rt =t12*c + t13*s */
		__asm	subpd	xmm3,xmm1	/* it =t13*c - t12*s */

		__asm	movaps	xmm0,[edi      ]	/* t02 */
		__asm	movaps	xmm1,[edi+0x010]	/* t03 */

		__asm	subpd	xmm0,xmm2			/*~t02=t02-rt */
		__asm	subpd	xmm1,xmm3			/*~t03=t03-it */
		__asm	addpd	xmm2,xmm2			/*      2* rt */
		__asm	addpd	xmm3,xmm3			/*      2* it */
		__asm	addpd	xmm2,xmm0			/*~t12=rt +t02*/
		__asm	addpd	xmm3,xmm1			/*~t13=it +t03*/

		__asm	subpd	xmm0,xmm4		/* t02-t22 */			__asm	subpd	xmm2,xmm7		/* t12-t33 */
		__asm	subpd	xmm1,xmm5		/* t03-t23 */			__asm	subpd	xmm3,xmm6		/* t13-t32 */
		__asm	addpd	xmm4,xmm4		/*   2*t22 */			__asm	addpd	xmm7,xmm7		/*   2*t33 */
		__asm	addpd	xmm5,xmm5		/*   2*t23 */			__asm	addpd	xmm6,xmm6		/*   2*t32 */
		__asm	movaps	[ebx     ],xmm0	/* a[jt+p1 ] */			__asm	movaps	[ecx     ],xmm2	/* a[jt+p2 ] */
		__asm	movaps	[ebx+0x10],xmm1	/* a[jp+p1 ] */			__asm	movaps	[edx+0x10],xmm3	/* a[jp+p3 ] */
		__asm	addpd	xmm4,xmm0		/* t02+t22 */			__asm	addpd	xmm7,xmm2		/* t12+t33 */
		__asm	addpd	xmm5,xmm1		/* t03+t23 */			__asm	addpd	xmm6,xmm3		/* t13+t32 */
		__asm	movaps	[eax     ],xmm4	/* a[jt+p0 ] */			__asm	movaps	[edx     ],xmm7	/* a[jt+p3 ] */
		__asm	movaps	[eax+0x10],xmm5	/* a[jp+p0 ] */			__asm	movaps	[ecx+0x10],xmm6	/* a[jp+p2 ] */

	#else	/* GCC-style inline ASM: */

		add0 = &a[j1];
		SSE2_RADIX32_DIF_TWIDDLE(add0,p01,p02,p03,p04,p08,p0C,p10,p18,r00)

	#endif

	#if 0//DEBUG_SSE2
		if(s01->re != 0.0)
		{
			fprintf(stderr, "radix32_dif_pass: r00 = %20.5f, %20.5f\n", r00->re, r01->re);
			fprintf(stderr, "radix32_dif_pass: r02 = %20.5f, %20.5f\n", r02->re, r03->re);
			fprintf(stderr, "radix32_dif_pass: r04 = %20.5f, %20.5f\n", r04->re, r05->re);
			fprintf(stderr, "radix32_dif_pass: r06 = %20.5f, %20.5f\n", r06->re, r07->re);
			fprintf(stderr, "radix32_dif_pass: r08 = %20.5f, %20.5f\n", r08->re, r09->re);
			fprintf(stderr, "radix32_dif_pass: r0A = %20.5f, %20.5f\n", r0A->re, r0B->re);
			fprintf(stderr, "radix32_dif_pass: r0C = %20.5f, %20.5f\n", r0C->re, r0D->re);
			fprintf(stderr, "radix32_dif_pass: r0E = %20.5f, %20.5f\n", r0E->re, r0F->re);

			fprintf(stderr, "radix32_dif_pass: r10 = %20.5f, %20.5f\n", r10->re, r11->re);
			fprintf(stderr, "radix32_dif_pass: r12 = %20.5f, %20.5f\n", r12->re, r13->re);
			fprintf(stderr, "radix32_dif_pass: r14 = %20.5f, %20.5f\n", r14->re, r15->re);
			fprintf(stderr, "radix32_dif_pass: r16 = %20.5f, %20.5f\n", r16->re, r17->re);
			fprintf(stderr, "radix32_dif_pass: r18 = %20.5f, %20.5f\n", r18->re, r19->re);
			fprintf(stderr, "radix32_dif_pass: r1A = %20.5f, %20.5f\n", r1A->re, r1B->re);
			fprintf(stderr, "radix32_dif_pass: r1C = %20.5f, %20.5f\n", r1C->re, r1D->re);
			fprintf(stderr, "radix32_dif_pass: r1E = %20.5f, %20.5f\n", r1E->re, r1F->re);

			fprintf(stderr, "radix32_dif_pass: r20 = %20.5f, %20.5f\n", r20->re, r21->re);
			fprintf(stderr, "radix32_dif_pass: r22 = %20.5f, %20.5f\n", r22->re, r23->re);
			fprintf(stderr, "radix32_dif_pass: r24 = %20.5f, %20.5f\n", r24->re, r25->re);
			fprintf(stderr, "radix32_dif_pass: r26 = %20.5f, %20.5f\n", r26->re, r27->re);
			fprintf(stderr, "radix32_dif_pass: r28 = %20.5f, %20.5f\n", r28->re, r29->re);
			fprintf(stderr, "radix32_dif_pass: r2A = %20.5f, %20.5f\n", r2A->re, r2B->re);
			fprintf(stderr, "radix32_dif_pass: r2C = %20.5f, %20.5f\n", r2C->re, r2D->re);
			fprintf(stderr, "radix32_dif_pass: r2E = %20.5f, %20.5f\n", r2E->re, r2F->re);

			fprintf(stderr, "radix32_dif_pass: r30 = %20.5f, %20.5f\n", r30->re, r31->re);
			fprintf(stderr, "radix32_dif_pass: r32 = %20.5f, %20.5f\n", r32->re, r33->re);
			fprintf(stderr, "radix32_dif_pass: r34 = %20.5f, %20.5f\n", r34->re, r35->re);
			fprintf(stderr, "radix32_dif_pass: r36 = %20.5f, %20.5f\n", r36->re, r37->re);
			fprintf(stderr, "radix32_dif_pass: r38 = %20.5f, %20.5f\n", r38->re, r39->re);
			fprintf(stderr, "radix32_dif_pass: r3A = %20.5f, %20.5f\n", r3A->re, r3B->re);
			fprintf(stderr, "radix32_dif_pass: r3C = %20.5f, %20.5f\n", r3C->re, r3D->re);
			fprintf(stderr, "radix32_dif_pass: r0E = %20.5f, %20.5f\n", r0E->re, r0F->re);
			exit(0);
		}
	#endif

#else

	/*       gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...
			 We process the sincos data in bit-reversed order.	*/
	#ifdef PFETCH_AGGRESSIVE
		addr = &a[j1];
	#elif PFETCH
		prefetch_offset = ((j >> 1) & 7)*p04 + 4;	/* Cycle among p00, p04, p08, p0C, p10, p14, p18 and p1C. */
	#endif
	/*...Block 1: */

		jt = j1;
		jp = j2;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr;
		prefetch_p_doubles(addp);
	#elif PFETCH
		addr = &a[jt];
		addp = addr+prefetch_offset;
		prefetch_p_doubles(addp);
	#endif
			t00=a[jt    ];									t01=a[jp    ];
			rt =a[jt+p10]*c10 - a[jp+p10]*s10;	it =a[jp+p10]*c10 + a[jt+p10]*s10;
			t02=t00-rt;		t03=t01-it;
			t00=t00+rt;		t01=t01+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t04=a[jt+p08]*c08 - a[jp+p08]*s08;	t05=a[jp+p08]*c08 + a[jt+p08]*s08;
			rt =a[jt+p18]*c18 - a[jp+p18]*s18;	it =a[jp+p18]*c18 + a[jt+p18]*s18;
			t06=t04-rt;		t07=t05-it;
			t04=t04+rt;		t05=t05+it;

			rt =t04;		it =t05;
			t04=t00-rt;		t05=t01-it;
			t00=t00+rt;		t01=t01+it;

			rt =t06;		it =t07;
			t06=t02+it;		t07=t03-rt;
			t02=t02-it;		t03=t03+rt;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t08=a[jt+p04]*c04 - a[jp+p04]*s04;	t09=a[jp+p04]*c04 + a[jt+p04]*s04;
			rt =a[jt+p14]*c14 - a[jp+p14]*s14;	it =a[jp+p14]*c14 + a[jt+p14]*s14;
			t0A=t08-rt;		t0B=t09-it;
			t08=t08+rt;		t09=t09+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t0C=a[jt+p0C]*c0C - a[jp+p0C]*s0C;	t0D=a[jp+p0C]*c0C + a[jt+p0C]*s0C;
			rt =a[jt+p1C]*c1C - a[jp+p1C]*s1C;	it =a[jp+p1C]*c1C + a[jt+p1C]*s1C;
			t0E=t0C-rt;		t0F=t0D-it;
			t0C=t0C+rt;		t0D=t0D+it;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p04;
		prefetch_p_doubles(addp);
	#endif
			rt =t0C;		it =t0D;
			t0C=t08-rt;		t0D=t09-it;
			t08=t08+rt;		t09=t09+it;

			rt =t0E;		it =t0F;
			t0E=t0A+it;		t0F=t0B-rt;
			t0A=t0A-it;		t0B=t0B+rt;

			rt =t08;		it =t09;
			t08=t00-rt;		t09=t01-it;
			t00=t00+rt;		t01=t01+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t0C;		it =t0D;
			t0C=t04+it;		t0D=t05-rt;
			t04=t04-it;		t05=t05+rt;

			rt =(t0A-t0B)*ISRT2;it =(t0A+t0B)*ISRT2;
			t0A=t02-rt;		t0B=t03-it;
			t02=t02+rt;		t03=t03+it;

			rt =(t0E+t0F)*ISRT2;it =(t0F-t0E)*ISRT2;
			t0E=t06+rt;		t0F=t07+it;
			t06=t06-rt;		t07=t07-it;

	/*...Block 2:	*/

		jt = j1 + p02;
		jp = j2 + p02;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#elif PFETCH
		addr = &a[jt];
		addp = addr+prefetch_offset;
		prefetch_p_doubles(addp);
	#endif
			t10=a[jt    ]*c02 - a[jp    ]*s02;	t11=a[jp    ]*c02 + a[jt    ]*s02;
			rt =a[jt+p10]*c12 - a[jp+p10]*s12;	it =a[jp+p10]*c12 + a[jt+p10]*s12;
			t12=t10-rt;		t13=t11-it;
			t10=t10+rt;		t11=t11+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t14=a[jt+p08]*c0A - a[jp+p08]*s0A;	t15=a[jp+p08]*c0A + a[jt+p08]*s0A;
			rt =a[jt+p18]*c1A - a[jp+p18]*s1A;	it =a[jp+p18]*c1A + a[jt+p18]*s1A;
			t16=t14-rt;		t17=t15-it;
			t14=t14+rt;		t15=t15+it;

			rt =t14;		it =t15;
			t14=t10-rt;		t15=t11-it;
			t10=t10+rt;		t11=t11+it;

			rt =t16;		it =t17;
			t16=t12+it;		t17=t13-rt;
			t12=t12-it;		t13=t13+rt;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p08;
		prefetch_p_doubles(addp);
	#endif
			t18=a[jt+p04]*c06 - a[jp+p04]*s06;	t19=a[jp+p04]*c06 + a[jt+p04]*s06;
			rt =a[jt+p14]*c16 - a[jp+p14]*s16;	it =a[jp+p14]*c16 + a[jt+p14]*s16;
			t1A=t18-rt;		t1B=t19-it;
			t18=t18+rt;		t19=t19+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t1C=a[jt+p0C]*c0E - a[jp+p0C]*s0E;	t1D=a[jp+p0C]*c0E + a[jt+p0C]*s0E;
			rt =a[jt+p1C]*c1E - a[jp+p1C]*s1E;	it =a[jp+p1C]*c1E + a[jt+p1C]*s1E;
			t1E=t1C-rt;		t1F=t1D-it;
			t1C=t1C+rt;		t1D=t1D+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t1C;		it =t1D;
			t1C=t18-rt;		t1D=t19-it;
			t18=t18+rt;		t19=t19+it;

			rt =t1E;		it =t1F;
			t1E=t1A+it;		t1F=t1B-rt;
			t1A=t1A-it;		t1B=t1B+rt;

			rt =t18;		it =t19;
			t18=t10-rt;		t19=t11-it;
			t10=t10+rt;		t11=t11+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t1C;		it =t1D;
			t1C=t14+it;		t1D=t15-rt;
			t14=t14-it;		t15=t15+rt;

			rt =(t1A-t1B)*ISRT2;it =(t1A+t1B)*ISRT2;
			t1A=t12-rt;		t1B=t13-it;
			t12=t12+rt;		t13=t13+it;

			rt =(t1E+t1F)*ISRT2;it =(t1F-t1E)*ISRT2;
			t1E=t16+rt;		t1F=t17+it;
			t16=t16-rt;		t17=t17-it;

	/*...Block 3:	*/

		jt = j1 + p01;
		jp = j2 + p01;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p0C;
		prefetch_p_doubles(addp);
	#elif PFETCH
		addr = &a[jt];
		addp = addr+prefetch_offset;
		prefetch_p_doubles(addp);
	#endif
			t20=a[jt    ]*c01 - a[jp    ]*s01;	t21=a[jp    ]*c01 + a[jt    ]*s01;
			rt =a[jt+p10]*c11 - a[jp+p10]*s11;	it =a[jp+p10]*c11 + a[jt+p10]*s11;
			t22=t20-rt;		t23=t21-it;
			t20=t20+rt;		t21=t21+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t24=a[jt+p08]*c09 - a[jp+p08]*s09;	t25=a[jp+p08]*c09 + a[jt+p08]*s09;
			rt =a[jt+p18]*c19 - a[jp+p18]*s19;	it =a[jp+p18]*c19 + a[jt+p18]*s19;
			t26=t24-rt;		t27=t25-it;
			t24=t24+rt;		t25=t25+it;

			rt =t24;		it =t25;
			t24=t20-rt;		t25=t21-it;
			t20=t20+rt;		t21=t21+it;

			rt =t26;		it =t27;
			t26=t22+it;		t27=t23-rt;
			t22=t22-it;		t23=t23+rt;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t28=a[jt+p04]*c05 - a[jp+p04]*s05;	t29=a[jp+p04]*c05 + a[jt+p04]*s05;
			rt =a[jt+p14]*c15 - a[jp+p14]*s15;	it =a[jp+p14]*c15 + a[jt+p14]*s15;
			t2A=t28-rt;		t2B=t29-it;
			t28=t28+rt;		t29=t29+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t2C=a[jt+p0C]*c0D - a[jp+p0C]*s0D;	t2D=a[jp+p0C]*c0D + a[jt+p0C]*s0D;
			rt =a[jt+p1C]*c1D - a[jp+p1C]*s1D;	it =a[jp+p1C]*c1D + a[jt+p1C]*s1D;
			t2E=t2C-rt;		t2F=t2D-it;
			t2C=t2C+rt;		t2D=t2D+it;
	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p10;
		prefetch_p_doubles(addp);
	#endif

			rt =t2C;		it =t2D;
			t2C=t28-rt;		t2D=t29-it;
			t28=t28+rt;		t29=t29+it;

			rt =t2E;		it =t2F;
			t2E=t2A+it;		t2F=t2B-rt;
			t2A=t2A-it;		t2B=t2B+rt;

			rt =t28;		it =t29;
			t28=t20-rt;		t29=t21-it;
			t20=t20+rt;		t21=t21+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t2C;		it =t2D;
			t2C=t24+it;		t2D=t25-rt;
			t24=t24-it;		t25=t25+rt;

			rt =(t2A-t2B)*ISRT2;it =(t2A+t2B)*ISRT2;
			t2A=t22-rt;		t2B=t23-it;
			t22=t22+rt;		t23=t23+it;

			rt =(t2E+t2F)*ISRT2;it =(t2F-t2E)*ISRT2;
			t2E=t26+rt;		t2F=t27+it;
			t26=t26-rt;		t27=t27-it;

	/*...Block 4:	*/

		jt = j1 + p03;
		jp = j2 + p03;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#elif PFETCH
		addr = &a[jt];
		addp = addr+prefetch_offset;
		prefetch_p_doubles(addp);
	#endif
			t30=a[jt    ]*c03 - a[jp    ]*s03;	t31=a[jp    ]*c03 + a[jt    ]*s03;
			rt =a[jt+p10]*c13 - a[jp+p10]*s13;	it =a[jp+p10]*c13 + a[jt+p10]*s13;
			t32=t30-rt;		t33=t31-it;
			t30=t30+rt;		t31=t31+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t34=a[jt+p08]*c0B - a[jp+p08]*s0B;	t35=a[jp+p08]*c0B + a[jt+p08]*s0B;
			rt =a[jt+p18]*c1B - a[jp+p18]*s1B;	it =a[jp+p18]*c1B + a[jt+p18]*s1B;
			t36=t34-rt;		t37=t35-it;
			t34=t34+rt;		t35=t35+it;

			rt =t34;		it =t35;
			t34=t30-rt;		t35=t31-it;
			t30=t30+rt;		t31=t31+it;

			rt =t36;		it =t37;
			t36=t32+it;		t37=t33-rt;
			t32=t32-it;		t33=t33+rt;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p14;
		prefetch_p_doubles(addp);
	#endif
			t38=a[jt+p04]*c07 - a[jp+p04]*s07;	t39=a[jp+p04]*c07 + a[jt+p04]*s07;
			rt =a[jt+p14]*c17 - a[jp+p14]*s17;	it =a[jp+p14]*c17 + a[jt+p14]*s17;
			t3A=t38-rt;		t3B=t39-it;
			t38=t38+rt;		t39=t39+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			t3C=a[jt+p0C]*c0F - a[jp+p0C]*s0F;	t3D=a[jp+p0C]*c0F + a[jt+p0C]*s0F;
			rt =a[jt+p1C]*c1F - a[jp+p1C]*s1F;	it =a[jp+p1C]*c1F + a[jt+p1C]*s1F;
			t3E=t3C-rt;		t3F=t3D-it;
			t3C=t3C+rt;		t3D=t3D+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t3C;		it =t3D;
			t3C=t38-rt;		t3D=t39-it;
			t38=t38+rt;		t39=t39+it;

			rt =t3E;		it =t3F;
			t3E=t3A+it;		t3F=t3B-rt;
			t3A=t3A-it;		t3B=t3B+rt;

			rt =t38;		it =t39;
			t38=t30-rt;		t39=t31-it;
			t30=t30+rt;		t31=t31+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t3C;		it =t3D;
			t3C=t34+it;		t3D=t35-rt;
			t34=t34-it;		t35=t35+rt;

			rt =(t3A-t3B)*ISRT2;it =(t3A+t3B)*ISRT2;
			t3A=t32-rt;		t3B=t33-it;
			t32=t32+rt;		t33=t33+it;

			rt =(t3E+t3F)*ISRT2;it =(t3F-t3E)*ISRT2;
			t3E=t36+rt;		t3F=t37+it;
			t36=t36-rt;		t37=t37-it;

	/*...and now do eight radix-4 transforms, including the internal twiddle factors:
		1, exp(i* 1*twopi/32) =       ( c32_1, s32_1), exp(i* 2*twopi/32) =       ( c    , s    ), exp(i* 3*twopi/32) =       ( c32_3, s32_3) (for inputs to transform block 2),
		1, exp(i* 2*twopi/32) =       ( c    , s    ), exp(i* 4*twopi/32) = ISRT2*( 1    , 1    ), exp(i* 6*twopi/32) =       ( s    , c    ) (for inputs to transform block 3),
		1, exp(i* 3*twopi/32) =       ( c32_3, s32_3), exp(i* 6*twopi/32) =       ( s    , c    ), exp(i* 9*twopi/32) =       (-s32_1, c32_1) (for inputs to transform block 4),
		1, exp(i* 4*twopi/32) = ISRT2*( 1    , 1    ), exp(i* 8*twopi/32) =       ( 0    , 1    ), exp(i*12*twopi/32) = ISRT2*(-1    , 1    ) (for inputs to transform block 5),
		1, exp(i* 5*twopi/32) =       ( s32_3, c32_3), exp(i*10*twopi/32) =       (-s    , c    ), exp(i*15*twopi/32) =       (-c32_1, s32_1) (for inputs to transform block 6),
		1, exp(i* 6*twopi/32) =       ( s    , c    ), exp(i*12*twopi/32) = ISRT2*(-1    , 1    ), exp(i*18*twopi/32) =       (-c    ,-s    ) (for inputs to transform block 7),
		1, exp(i* 7*twopi/32) =       ( s32_1, c32_1), exp(i*14*twopi/32) =       (-c    , s    ), exp(i*21*twopi/32) =       (-s32_3,-c32_3) (for inputs to transform block 8),
		 and only the last 3 inputs to each of the radix-4 transforms 2 through 8 are multiplied by non-unity twiddles.	*/

	/*...Block 1: t00,t10,t20,t30	*/

		jt = j1;
		jp = j2;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p18;
		prefetch_p_doubles(addp);
	#endif
			rt =t10;	t10=t00-rt;	t00=t00+rt;
			it =t11;	t11=t01-it;	t01=t01+it;

			rt =t30;	t30=t20-rt;	t20=t20+rt;
			it =t31;	t31=t21-it;	t21=t21+it;

			a[jt    ]=t00+t20;		a[jp    ]=t01+t21;
			a[jt+p01]=t00-t20;		a[jp+p01]=t01-t21;

			a[jt+p02]=t10-t31;		a[jp+p02]=t11+t30;	/* mpy by E^4=i is inlined here...	*/
			a[jt+p03]=t10+t31;		a[jp+p03]=t11-t30;

	/*...Block 5: t08,t18,t28,t38	*/

		jt = j1 + p04;
		jp = j2 + p04;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t18;	t18=t08+t19;	t08=t08-t19;		/* twiddle mpy by E^4 = I	*/
				t19=t09-rt;	t09=t09+rt;

			rt =(t28-t29)*ISRT2;	t29=(t28+t29)*ISRT2;		t28=rt;	/* twiddle mpy by E^8	*/
			rt =(t39+t38)*ISRT2;	it =(t39-t38)*ISRT2;			/* twiddle mpy by -E^12 is here...	*/
			t38=t28+rt;			t28=t28-rt;				/* ...and get E^12 by flipping signs here.	*/
			t39=t29+it;			t29=t29-it;

			a[jt    ]=t08+t28;		a[jp    ]=t09+t29;
			a[jt+p01]=t08-t28;		a[jp+p01]=t09-t29;

			a[jt+p02]=t18-t39;		a[jp+p02]=t19+t38;	/* mpy by E^4=i is inlined here...	*/
			a[jt+p03]=t18+t39;		a[jp+p03]=t19-t38;

	/*...Block 3: t04,t14,t24,t34	*/

		jt = j1 + p08;
		jp = j2 + p08;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =(t14-t15)*ISRT2;	it =(t14+t15)*ISRT2;			/* twiddle mpy by E^4	*/
			t14=t04-rt;			t04=t04+rt;
			t15=t05-it;			t05=t05+it;

			rt =t24*c - t25*s;		t25=t25*c + t24*s;		t24=rt;	/* twiddle mpy by E^2	*/
			rt =t34*s - t35*c;		it =t35*s + t34*c;			/* twiddle mpy by E^6	*/
			t34=t24-rt;			t24=t24+rt;
			t35=t25-it;			t25=t25+it;

			a[jt    ]=t04+t24;		a[jp    ]=t05+t25;
			a[jt+p01]=t04-t24;		a[jp+p01]=t05-t25;

			a[jt+p02]=t14-t35;		a[jp+p02]=t15+t34;	/* mpy by E^4=i is inlined here...	*/
			a[jt+p03]=t14+t35;		a[jp+p03]=t15-t34;

	/*...Block 7: t0C,t1C,t2C,t3C	*/

		jt = j1 + p0C;
		jp = j2 + p0C;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =(t1D+t1C)*ISRT2;	it =(t1D-t1C)*ISRT2;			/* twiddle mpy by -E^12 is here...	*/
			t1C=t0C+rt;			t0C=t0C-rt;				/* ...and get E^12 by flipping signs here.	*/
			t1D=t0D+it;			t0D=t0D-it;

			rt =t2C*s - t2D*c;		t2D=t2D*s + t2C*c;		t2C=rt;	/* twiddle mpy by E^6	*/
			rt =t3C*c - t3D*s;		it =t3D*c + t3C*s;			/* twiddle mpy by E^18 is here...	*/
			t3C=t2C+rt;			t2C=t2C-rt;				/* ...and get E^18 by flipping signs here.	*/
			t3D=t2D+it;			t2D=t2D-it;

			a[jt    ]=t0C+t2C;		a[jp    ]=t0D+t2D;
			a[jt+p01]=t0C-t2C;		a[jp+p01]=t0D-t2D;

			a[jt+p02]=t1C-t3D;		a[jp+p02]=t1D+t3C;	/* mpy by E^4=i is inlined here...	*/
			a[jt+p03]=t1C+t3D;		a[jp+p03]=t1D-t3C;

	/*...Block 2: t02,t12,t22,t32	*/

		jt = j1 + p10;
		jp = j2 + p10;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p1C;
		prefetch_p_doubles(addp);
	#endif
			rt =t12*c - t13*s;		it =t13*c + t12*s;			/* twiddle mpy by E^2	*/
			t12=t02-rt;			t02=t02+rt;
			t13=t03-it;			t03=t03+it;

			rt =t22*c32_1 - t23*s32_1;	t23=t23*c32_1 + t22*s32_1;	t22=rt;	/* twiddle mpy by E^1	*/
			rt =t32*c32_3 - t33*s32_3;	it =t33*c32_3 + t32*s32_3;		/* twiddle mpy by E^3	*/
			t32=t22-rt;			t22=t22+rt;
			t33=t23-it;			t23=t23+it;

			a[jt    ]=t02+t22;		a[jp    ]=t03+t23;
			a[jt+p01]=t02-t22;		a[jp+p01]=t03-t23;

			a[jt+p02]=t12-t33;		a[jp+p02]=t13+t32;	/* mpy by E^4=i is inlined here...;	*/
			a[jt+p03]=t12+t33;		a[jp+p03]=t13-t32;

	/*...Block 6: t0A,t1A,t2A,t3A	*/

		jt = j1 + p14;
		jp = j2 + p14;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t1A*s + t1B*c;		it =t1B*s - t1A*c;			/* twiddle mpy by -E^10 is here...	*/
			t1A=t0A+rt;			t0A =t0A-rt;				/* ...and get E^10 by flipping signs here.	*/
			t1B=t0B+it;			t0B =t0B-it;

			rt =t2A*s32_3 - t2B*c32_3;	t2B=t2B*s32_3 + t2A*c32_3;	t2A=rt;	/* twiddle mpy by E^5	*/
			rt =t3A*c32_1 + t3B*s32_1;	it =t3B*c32_1 - t3A*s32_1;	 	/* twiddle mpy by -E^15 is here...	*/
			t3A=t2A+rt;			t2A=t2A-rt;				/* ...and get E^15 by flipping signs here.	*/
			t3B=t2B+it;			t2B=t2B-it;

			a[jt    ]=t0A+t2A;		a[jp    ]=t0B+t2B;
			a[jt+p01]=t0A-t2A;		a[jp+p01]=t0B-t2B;

			a[jt+p02]=t1A-t3B;		a[jp+p02]=t1B+t3A;	/* mpy by E^4=i is inlined here...	*/
			a[jt+p03]=t1A+t3B;		a[jp+p03]=t1B-t3A;

	/*...Block 4: t06,t16,t26,t36	*/

		jt = j1 + p18;
		jp = j2 + p18;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t16*s - t17*c;		it =t17*s + t16*c;			/* twiddle mpy by E^6	*/
			t16=t06-rt;			t06 =t06+rt;
			t17=t07-it;			t07 =t07+it;

			rt =t26*c32_3 - t27*s32_3;	t27=t27*c32_3 + t26*s32_3;	t26=rt;	/* twiddle mpy by E^3	*/
			rt =t36*s32_1 + t37*c32_1;	it =t37*s32_1 - t36*c32_1;		/* twiddle mpy by -E^9 is here...	*/
			t36=t26+rt;			t26=t26-rt;				/* ...and get E^9 by flipping signs here.	*/
			t37=t27+it;			t27=t27-it;

			a[jt    ]=t06+t26;		a[jp    ]=t07+t27;
			a[jt+p01]=t06-t26;		a[jp+p01]=t07-t27;

			a[jt+p02]=t16-t37;		a[jp+p02]=t17+t36;	/* mpy by E^4=i is inlined here...	*/
			a[jt+p03]=t16+t37;		a[jp+p03]=t17-t36;

	/*...Block 8: t0E,t1E,t2E,t3E	*/

		jt = j1 + p1C;
		jp = j2 + p1C;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t1E*c + t1F*s;		it =t1F*c - t1E*s;			/* twiddle mpy by -E^14 is here...	*/
			t1E=t0E+rt;			t0E =t0E-rt;				/* ...and get E^14 by flipping signs here.	*/
			t1F=t0F+it;			t0F =t0F-it;

			rt =t2E*s32_1 - t2F*c32_1;	t2F=t2F*s32_1 + t2E*c32_1;	t2E=rt;	/* twiddle mpy by E^7	*/
			rt =t3E*s32_3 - t3F*c32_3;	it =t3F*s32_3 + t3E*c32_3;		/* twiddle mpy by -E^21 is here...	*/
			t3E=t2E+rt;			t2E=t2E-rt;				/* ...and get E^21 by flipping signs here.	*/
			t3F=t2F+it;			t2F=t2F-it;

			a[jt    ]=t0E+t2E;		a[jp    ]=t0F+t2F;
			a[jt+p01]=t0E-t2E;		a[jp+p01]=t0F-t2F;

			a[jt+p02]=t1E-t3F;		a[jp+p02]=t1F+t3E;	/* mpy by E^4=i is inlined here...	*/
			a[jt+p03]=t1E+t3F;		a[jp+p03]=t1F-t3E;

#endif	/* USE_SSE2 */

	  }	/* endfor(j=jlo; j < jhi; j += 4) */

	  /*jlo=jlo+incr; jhi=jhi+incr; See my note about OpenMP above. */

	}	/* endfor(m=0; m < nloops; m++) */

}

/***************/

/*
!...Post-twiddles implementation of radix32_dit_pass.
*/
void radix32_dit_pass(double a[], int n, struct complex rt0[], struct complex rt1[], int index[], int nloops, int incr, int init_sse2, int thr_id)
{
	static int max_threads = 0;
	int i,j,j1,j2,jlo,jhi,m,iroot_prim,iroot,k1,k2;
	int p01,p02,p03,p04,p05,p06,p07,p08,p10,p18;
	const double c = 0.92387953251128675613, s     = 0.38268343236508977173	/* exp[  i*(twopi/16)]	*/
			,c32_1 = 0.98078528040323044912, s32_1 = 0.19509032201612826784	/* exp(  i*twopi/32), the radix-32 fundamental sincos datum	*/
			,c32_3 = 0.83146961230254523708, s32_3 = 0.55557023301960222473;/* exp(3*i*twopi/32)	*/
	double rt,it,re0,im0,re1,im1;

#ifdef USE_SSE2

  #if !(defined(COMPILER_TYPE_MSVC) || defined(COMPILER_TYPE_GCC) || defined(COMPILER_TYPE_SUNC))
	#error SSE2 code not supported for this compiler!
  #endif

	static struct complex *sc_arr = 0x0, *sc_ptr;
	double *add0;	/* Addresses into array sections */
#ifdef COMPILER_TYPE_MSVC
	double *add1, *add2, *add3, *add4, *add5, *add6, *add7;
#endif
	struct complex *c_tmp,*s_tmp;

  #ifdef MULTITHREAD
	static struct complex *__r0;	/* Base address for discrete per-thread local stores */
	// In || mode, only above base-pointer (shared by all threads) is static:
	struct complex *isrt2, *cc0, *ss0, *cc1, *ss1, *cc3, *ss3, *two
		,*r00,*r02,*r04,*r06,*r08,*r0A,*r0C,*r0E
		,*r10,*r12,*r14,*r16,*r18,*r1A,*r1C,*r1E,*r20,*r30;
  #else
	static struct complex *isrt2, *two, *cc0, *ss0, *cc1, *ss1, *cc3, *ss3
		,*c00,*c01,*c02,*c03,*c04,*c05,*c06,*c07,*c08,*c09,*c0A,*c0B,*c0C,*c0D,*c0E,*c0F
		,*c10,*c11,*c12,*c13,*c14,*c15,*c16,*c17,*c18,*c19,*c1A,*c1B,*c1C,*c1D,*c1E,*c1F
		,*s00,*s01,*s02,*s03,*s04,*s05,*s06,*s07,*s08,*s09,*s0A,*s0B,*s0C,*s0D,*s0E,*s0F
		,*s10,*s11,*s12,*s13,*s14,*s15,*s16,*s17,*s18,*s19,*s1A,*s1B,*s1C,*s1D,*s1E,*s1F
		,*r00,*r01,*r02,*r03,*r04,*r05,*r06,*r07,*r08,*r09,*r0A,*r0B,*r0C,*r0D,*r0E,*r0F
		,*r10,*r11,*r12,*r13,*r14,*r15,*r16,*r17,*r18,*r19,*r1A,*r1B,*r1C,*r1D,*r1E,*r1F
		,*r20,*r21,*r22,*r23,*r24,*r25,*r26,*r27,*r28,*r29,*r2A,*r2B,*r2C,*r2D,*r2E,*r2F
		,*r30,*r31,*r32,*r33,*r34,*r35,*r36,*r37,*r38,*r39,*r3A,*r3B,*r3C,*r3D,*r3E,*r3F;
  #endif

#else

	int jp,jt;
	double *addr, *addp;
	int prefetch_offset;
	double	 c01,c02,c03,c04,c05,c06,c07,c08,c09,c0A,c0B,c0C,c0D,c0E,c0F
		,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c1A,c1B,c1C,c1D,c1E,c1F
		    ,s01,s02,s03,s04,s05,s06,s07,s08,s09,s0A,s0B,s0C,s0D,s0E,s0F
		,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,s1A,s1B,s1C,s1D,s1E,s1F
		,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t0A,t0B,t0C,t0D,t0E,t0F
		,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t1A,t1B,t1C,t1D,t1E,t1F
		,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t2A,t2B,t2C,t2D,t2E,t2F
		,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t3A,t3B,t3C,t3D,t3E,t3F;

#endif

#ifdef USE_SSE2

	/* In order to eschew complex thread-block-and-sync logic related to the local-store-init step in multithread mode,
	switch to a special-init-mode-call paradigm, in which the function is inited once (for as many threads as needed)
	prior to being executed:
	*/
	if(init_sse2 && !sc_ptr)	// Just check (init_sse2 != 0) here, to allow the *value* of init_sse2 to store #threads
	{
		max_threads = init_sse2;
	#ifndef COMPILER_TYPE_GCC
		ASSERT(HERE, NTHREADS == 1, "Multithreading currently only supported for GCC builds!");
	#endif
		ASSERT(HERE, max_threads >= NTHREADS, "Multithreading requires max_threads >= NTHREADS!");
		ASSERT(HERE, sc_arr == 0x0, "Init-mode call conflicts with already-malloc'ed local storage!");
		ASSERT(HERE, thr_id == -1, "Init-mode call must be outside of any multithreading!");
		sc_arr = ALLOC_COMPLEX(sc_arr, 0x90*max_threads);	if(!sc_arr){ sprintf(cbuf, "FATAL: unable to allocate sc_arr!.\n"); fprintf(stderr,"%s", cbuf);	ASSERT(HERE, 0,cbuf); }
		sc_ptr = ALIGN_COMPLEX(sc_arr);
		ASSERT(HERE, ((uint32)sc_ptr & 0x3f) == 0, "sc_ptr not 64-byte aligned!");

	/* Use low 64 16-byte slots of sc_arr for temporaries, next 7 for the nontrivial complex 32nd roots,
	last 64 for the doubled sincos twiddles, plus at least 3 more slots to allow for 64-byte alignment of the array.
	*/
		#ifdef MULTITHREAD
	//	if(max_threads > 1) {
			__r0  = sc_ptr;
			isrt2 = sc_ptr + 0x40;
			cc0	  = sc_ptr + 0x41;
			ss0	  = sc_ptr + 0x42;
			cc1	  = sc_ptr + 0x43;
			ss1	  = sc_ptr + 0x44;
			cc3	  = sc_ptr + 0x45;
			ss3	  = sc_ptr + 0x46;
			two   = sc_ptr + 0x87;
			for(i = 0; i < max_threads; ++i) {
				/* These remain fixed within each per-thread local store: */
				isrt2->re = ISRT2;	isrt2->im = ISRT2;
				cc0  ->re = c	;	cc0  ->im = c	;		ss0  ->re = s	;	ss0  ->im = s	;
				cc1  ->re = c32_1;	cc1  ->im = c32_1;		ss1  ->re = s32_1;	ss1  ->im = s32_1;
				cc3  ->re = c32_3;	cc3  ->im = c32_3;		ss3  ->re = s32_3;	ss3  ->im = s32_3;
				two  ->re = 2.0;	two  ->im = 2.0;
				/* Move on to next thread's local store */
				isrt2 += 0x90;
				cc0   += 0x90;
				ss0   += 0x90;
				cc1   += 0x90;
				ss1   += 0x90;
				cc3   += 0x90;
				ss3   += 0x90;
				two   += 0x90;
			}
		#else
	//	} else {
			r00		= sc_ptr + 0x00;	c00		= sc_ptr + 0x47;
			r01		= sc_ptr + 0x01;	s00		= sc_ptr + 0x48;
			r02		= sc_ptr + 0x02;	c08		= sc_ptr + 0x49;
			r03		= sc_ptr + 0x03;	s08		= sc_ptr + 0x4a;
			r04		= sc_ptr + 0x04;	c10		= sc_ptr + 0x4b;
			r05		= sc_ptr + 0x05;	s10		= sc_ptr + 0x4c;
			r06		= sc_ptr + 0x06;	c18		= sc_ptr + 0x4d;
			r07		= sc_ptr + 0x07;	s18		= sc_ptr + 0x4e;
			r08		= sc_ptr + 0x08;	c04		= sc_ptr + 0x4f;
			r09		= sc_ptr + 0x09;	s04		= sc_ptr + 0x50;
			r0A		= sc_ptr + 0x0a;	c0C		= sc_ptr + 0x51;
			r0B		= sc_ptr + 0x0b;	s0C		= sc_ptr + 0x52;
			r0C		= sc_ptr + 0x0c;	c14		= sc_ptr + 0x53;
			r0D		= sc_ptr + 0x0d;	s14		= sc_ptr + 0x54;
			r0E		= sc_ptr + 0x0e;	c1C		= sc_ptr + 0x55;
			r0F		= sc_ptr + 0x0f;	s1C		= sc_ptr + 0x56;
			r10		= sc_ptr + 0x10;	c02		= sc_ptr + 0x57;
			r11		= sc_ptr + 0x11;	s02		= sc_ptr + 0x58;
			r12		= sc_ptr + 0x12;	c0A		= sc_ptr + 0x59;
			r13		= sc_ptr + 0x13;	s0A		= sc_ptr + 0x5a;
			r14		= sc_ptr + 0x14;	c12		= sc_ptr + 0x5b;
			r15		= sc_ptr + 0x15;	s12		= sc_ptr + 0x5c;
			r16		= sc_ptr + 0x16;	c1A		= sc_ptr + 0x5d;
			r17		= sc_ptr + 0x17;	s1A		= sc_ptr + 0x5e;
			r18		= sc_ptr + 0x18;	c06		= sc_ptr + 0x5f;
			r19		= sc_ptr + 0x19;	s06		= sc_ptr + 0x60;
			r1A		= sc_ptr + 0x1a;	c0E		= sc_ptr + 0x61;
			r1B		= sc_ptr + 0x1b;	s0E		= sc_ptr + 0x62;
			r1C		= sc_ptr + 0x1c;	c16		= sc_ptr + 0x63;
			r1D		= sc_ptr + 0x1d;	s16		= sc_ptr + 0x64;
			r1E		= sc_ptr + 0x1e;	c1E		= sc_ptr + 0x65;
			r1F		= sc_ptr + 0x1f;	s1E		= sc_ptr + 0x66;
			r20		= sc_ptr + 0x20;	c01		= sc_ptr + 0x67;
			r21		= sc_ptr + 0x21;	s01		= sc_ptr + 0x68;
			r22		= sc_ptr + 0x22;	c09		= sc_ptr + 0x69;
			r23		= sc_ptr + 0x23;	s09		= sc_ptr + 0x6a;
			r24		= sc_ptr + 0x24;	c11		= sc_ptr + 0x6b;
			r25		= sc_ptr + 0x25;	s11		= sc_ptr + 0x6c;
			r26		= sc_ptr + 0x26;	c19		= sc_ptr + 0x6d;
			r27		= sc_ptr + 0x27;	s19		= sc_ptr + 0x6e;
			r28		= sc_ptr + 0x28;	c05		= sc_ptr + 0x6f;
			r29		= sc_ptr + 0x29;	s05		= sc_ptr + 0x70;
			r2A		= sc_ptr + 0x2a;	c0D		= sc_ptr + 0x71;
			r2B		= sc_ptr + 0x2b;	s0D		= sc_ptr + 0x72;
			r2C		= sc_ptr + 0x2c;	c15		= sc_ptr + 0x73;
			r2D		= sc_ptr + 0x2d;	s15		= sc_ptr + 0x74;
			r2E		= sc_ptr + 0x2e;	c1D		= sc_ptr + 0x75;
			r2F		= sc_ptr + 0x2f;	s1D		= sc_ptr + 0x76;
			r30		= sc_ptr + 0x30;	c03		= sc_ptr + 0x77;
			r31		= sc_ptr + 0x31;	s03		= sc_ptr + 0x78;
			r32		= sc_ptr + 0x32;	c0B		= sc_ptr + 0x79;
			r33		= sc_ptr + 0x33;	s0B		= sc_ptr + 0x7a;
			r34		= sc_ptr + 0x34;	c13		= sc_ptr + 0x7b;
			r35		= sc_ptr + 0x35;	s13		= sc_ptr + 0x7c;
			r36		= sc_ptr + 0x36;	c1B		= sc_ptr + 0x7d;
			r37		= sc_ptr + 0x37;	s1B		= sc_ptr + 0x7e;
			r38		= sc_ptr + 0x38;	c07		= sc_ptr + 0x7f;
			r39		= sc_ptr + 0x39;	s07		= sc_ptr + 0x80;
			r3A		= sc_ptr + 0x3a;	c0F		= sc_ptr + 0x81;
			r3B		= sc_ptr + 0x3b;	s0F		= sc_ptr + 0x82;
			r3C		= sc_ptr + 0x3c;	c17		= sc_ptr + 0x83;
			r3D		= sc_ptr + 0x3d;	s17		= sc_ptr + 0x84;
			r3E		= sc_ptr + 0x3e;	c1F		= sc_ptr + 0x85;
			r3F		= sc_ptr + 0x3f;	s1F		= sc_ptr + 0x86;
			isrt2	= sc_ptr + 0x40;	two     = sc_ptr + 0x87;
			cc0		= sc_ptr + 0x41;
			ss0		= sc_ptr + 0x42;
			cc1		= sc_ptr + 0x43;
			ss1		= sc_ptr + 0x44;
			cc3		= sc_ptr + 0x45;
			ss3		= sc_ptr + 0x46;

			/* These remain fixed: */
			isrt2->re = ISRT2;	isrt2->im = ISRT2;		two->re   = 2.0;	two->im   = 2.0;
			cc0  ->re = c	;	cc0  ->im = c	;		ss0  ->re = s	;	ss0  ->im = s	;
			cc1  ->re = c32_1;	cc1  ->im = c32_1;		ss1  ->re = s32_1;	ss1  ->im = s32_1;
			cc3  ->re = c32_3;	cc3  ->im = c32_3;		ss3  ->re = s32_3;	ss3  ->im = s32_3;
		#endif
	//	}
		return;
	}	/* end of inits */

	/* If multithreaded, set the local-store pointers needed for the current thread; */
	#ifdef MULTITHREAD
		ASSERT(HERE, (uint32)thr_id < (uint32)max_threads, "Bad thread ID!");
		r00 = __r0 + thr_id*0x90;
		r02 = r00 + 0x02;
		r04 = r00 + 0x04;
		r06 = r00 + 0x06;
		r08 = r00 + 0x08;
		r0A = r00 + 0x0A;
		r0C = r00 + 0x0C;
		r0E = r00 + 0x0E;
		r10 = r00 + 0x10;
		r12 = r00 + 0x12;
		r14 = r00 + 0x14;
		r16 = r00 + 0x16;
		r18 = r00 + 0x18;
		r1A = r00 + 0x1A;
		r1C = r00 + 0x1C;
		r1E = r00 + 0x1E;
		r20 = r00 + 0x20;
		r30 = r00 + 0x30;
		isrt2 = r00 + 0x40;
		cc0	= r00 + 0x41;
	#endif

#endif

	p01 = incr >> 5;
	p02 = p01 +p01;
	p03 = p02 +p01;
	p04 = p03 +p01;
	p05 = p04 +p01;
	p06 = p05 +p01;
	p07 = p06 +p01;
	p08 = p07 +p01;
	p10 = p08 +p08;
	p18 = p10 +p08;

	p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
	p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
	p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
	p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
	p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
	p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
	p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
	p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
	p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
	p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );

/*...The radix-32 pass is here.	*/

	iroot_prim=(incr >> 6);		/* (incr/2)/radix_now		*/

	for(m=0; m < nloops; m++)	/* NLOOPS may range from 1 (if first pass radix = 16) to P*N/32 (last pass radix = 16).x	*/
	{				/* NLOOPS satisfies the identity NLOOPS * INCR = P*N, which says that we process the entire
					   array each time this subroutine is executed (since P*N = vector length, sans padding.)	*/

/*	here are the needed sincos data - these are processed below in bit-reversed order.	*/
	  iroot = index[m]*iroot_prim;
	  i = iroot;

/* In the DIT pass we may be able to afford a more-accurate twiddles computation,
   since we only compute each set of twiddles once and re-use it for many data blocks. */
#if HIACC
	#ifdef USE_SSE2
		/* Due to roots-locality considerations, roots (c,s)[0-31] are offset w.r.to the thread-local ptr pair as
					cc[00 01 02 03 04 05 06 07 08 09 0A 0B 0C 0D 0E 0F 10 11 12 13 14 15 16 17 18 19 1A 1B 1C 1D 1E 1F]
		(cc0,ss0) + 0x[06,26,16,36|0e,2e,1e,3e|08,28,18,38|10,30,20,40|0a,2a,1a,3a|12,32,22,42|0c,2c,1c,3c|14,34,24,44].
	*** NOTE: This is the same pattern as for DIF version, but with the middle 2 roots octets [08-40] and [0a-42] swapped ***
		*/
		c_tmp = cc0 + 0x06; s_tmp = c_tmp+1;	/* c0,s0 */
		c_tmp->re = c_tmp->im = 1.0;
		s_tmp->re = s_tmp->im = 0.0;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x26; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c01=rt;		s01=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x16; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c02=rt;		s02=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x36; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c03=rt;		s03=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x0e; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c04=rt;		s04=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x2e; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c05=rt;		s05=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x1e; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c06=rt;		s06=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x3e; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c07=rt;		s07=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x08; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c08=rt;		s08=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x28; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c09=rt;		s09=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x18; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c0A=rt;		s0A=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x38; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c0B=rt;		s0B=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x10; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c0C=rt;		s0C=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x30; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c0D=rt;		s0D=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x20; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c0E=rt;		s0E=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x40; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c0F=rt;		s0F=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x0a; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c10=rt;		s10=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x2a; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c11=rt;		s11=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x1a; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c12=rt;		s12=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x3a; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c13=rt;		s13=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x12; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c14=rt;		s14=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x32; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c15=rt;		s15=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x22; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c16=rt;		s16=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x42; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c17=rt;		s17=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x0c; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c18=rt;		s18=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x2c; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c19=rt;		s19=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x1c; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c1A=rt;		s1A=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x3c; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c1B=rt;		s1B=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x14; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c1C=rt;		s1C=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x34; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c1D=rt;		s1D=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x24; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c1E=rt;		s1E=it;
	#endif

		k1=(i & NRTM1);
		k2=(i >> NRT_BITS);
		i += iroot;
		re0=rt0[k1].re;	im0=rt0[k1].im;
		re1=rt1[k2].re;	im1=rt1[k2].im;
		rt=re0*re1-im0*im1;	it=re0*im1+im0*re1;
	#ifdef USE_SSE2
		c_tmp = cc0 + 0x44; s_tmp = c_tmp+1;
		c_tmp->re = c_tmp->im = rt;
		s_tmp->re = s_tmp->im = it;
	#else
		c1F=rt;		s1F=it;
	#endif

#else

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;			/* 2*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c01=t00*rt -t01*it;	s01=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;				/* 3*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c02=t00*rt -t01*it;	s02=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += 4*iroot;				/* 7*iroot	*/
	    iroot = i;				/* 7*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c03=t00*rt -t01*it;	s03=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;				/* 14*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c07=t00*rt -t01*it;	s07=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;				/* 21*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c0E=t00*rt -t01*it;	s0E=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    i += iroot;				/* 28*iroot	*/
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c15=t00*rt -t01*it;	s15=t00*it +t01*rt;

	    k1=(i & NRTM1);
	    k2=(i >> NRT_BITS);
	    t00=rt0[k1].re;	t01=rt0[k1].im;
	    rt =rt1[k2].re;	it =rt1[k2].im;
	    c1C=t00*rt -t01*it;	s1C=t00*it +t01*rt;

/*(c,s)4-10 (decimal indices), 4-0A (hexadecimal indices):	*/
	    t00=c01*c07; t01=c01*s07; t02=s01*c07; t03=s01*s07;
	    c06=t00+t03; s06=t01-t02; c08=t00-t03; s08=t01+t02;

	    t00=c02*c07; t01=c02*s07; t02=s02*c07; t03=s02*s07;
	    c05=t00+t03; s05=t01-t02; c09=t00-t03; s09=t01+t02;

	    t00=c03*c07; t01=c03*s07; t02=s03*c07; t03=s03*s07;
	    c04=t00+t03; s04=t01-t02; c0A=t00-t03; s0A=t01+t02;

/*(c,s)11-17 (decimal indices), 0B-11 (hexadecimal indices):;	*/
	    t00=c01*c0E; t01=c01*s0E; t02=s01*c0E; t03=s01*s0E;
	    c0D=t00+t03; s0D=t01-t02; c0F=t00-t03; s0F=t01+t02;

	    t00=c02*c0E; t01=c02*s0E; t02=s02*c0E; t03=s02*s0E;
	    c0C=t00+t03; s0C=t01-t02; c10=t00-t03; s10=t01+t02;

	    t00=c03*c0E; t01=c03*s0E; t02=s03*c0E; t03=s03*s0E;
	    c0B=t00+t03; s0B=t01-t02; c11=t00-t03; s11=t01+t02;

/*(c,s)18-24 (decimal indices), 12-18 (hexadecimal indices):	*/
	    t00=c01*c15; t01=c01*s15; t02=s01*c15; t03=s01*s15;
	    c14=t00+t03; s14=t01-t02; c16=t00-t03; s16=t01+t02;

	    t00=c02*c15; t01=c02*s15; t02=s02*c15; t03=s02*s15;
	    c13=t00+t03; s13=t01-t02; c17=t00-t03; s17=t01+t02;

	    t00=c03*c15; t01=c03*s15; t02=s03*c15; t03=s03*s15;
	    c12=t00+t03; s12=t01-t02; c18=t00-t03; s18=t01+t02;

/*(c,s)25-31 (decimal indices), 19-1F (hexadecimal indices):	*/
	    t00=c01*c1C; t01=c01*s1C; t02=s01*c1C; t03=s01*s1C;
	    c1B=t00+t03; s1B=t01-t02; c1D=t00-t03; s1D=t01+t02;

	    t00=c02*c1C; t01=c02*s1C; t02=s02*c1C; t03=s02*s1C;
	    c1A=t00+t03; s1A=t01-t02; c1E=t00-t03; s1E=t01+t02;

	    t00=c03*c1C; t01=c03*s1C; t02=s03*c1C; t03=s03*s1C;
	    c19=t00+t03; s19=t01-t02; c1F=t00-t03; s1F=t01+t02;
#endif

	/* Define the inner-loop parameters in terms of the outer-loop ones to make OpenMP's job easier: */
	  jlo = m*incr;
	  jhi = jlo+(incr >> 5);

	#ifdef USE_SSE2
	  for(j=jlo; j < jhi; j += 4)
	  {
		/* In SSE2 mode, data are arranged in [re0,re1,im0,im1] quartets, not the usual [re0,im0],[re1,im1] pairs.
		Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
		but to point to the proper real datum, we need to bit-reverse bits <0:1> of j, i.e. [0,1,2,3] ==> [0,2,1,3].
		*/
		j1 = (j & mask01) + br4[j&3];
	#elif defined(USE_SSE2)	/* This allows us to use #if 0 above and disable sse2-based *computation*, while still using sse2-style data layout */
	  for(j=jlo; j < jhi; j += 2)
	  {
		j1 = (j & mask01) + br4[j&3];
	#else
	  for(j=jlo; j < jhi; j += 2)
	  {
		j1 =  j;
	#endif
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1+RE_IM_STRIDE;

#ifdef USE_SSE2

	#ifdef DEBUG_SSE2
		if(s01->re != 0.0)
		{
			rng_isaac_init(TRUE);
			jt = j1;		jp = j2;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = jt + p04;	jp = jp + p04;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p08;	jp = j2 + p08;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = jt + p04;	jp = jp + p04;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p10;	jp = j2 + p10;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = jt + p04;	jp = jp + p04;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = j1 + p18;	jp = j2 + p18;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	jt = jt + p04;	jp = jp + p04;
			a[jt    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp    ] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p01] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p02] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();
			a[jt+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();	a[jp+p03] = 1024.0*1024.0*rng_isaac_rand_double_norm_pm1();

			jt = j1; jp = j2;
			fprintf(stderr, "radix32_dit_pass: a_in[00] = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dit_pass: a_in[02] = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dit_pass: a_in[04] = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dit_pass: a_in[06] = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = jt + p04;	jp = jp + p04;
			fprintf(stderr, "radix32_dit_pass: a_in[08] = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dit_pass: a_in[0A] = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dit_pass: a_in[0C] = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dit_pass: a_in[0E] = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = j1 + p08;	jp = j2 + p08;
			fprintf(stderr, "radix32_dit_pass: a_in[10] = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dit_pass: a_in[12] = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dit_pass: a_in[14] = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dit_pass: a_in[16] = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = jt + p04;	jp = jp + p04;
			fprintf(stderr, "radix32_dit_pass: a_in[18] = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dit_pass: a_in[1A] = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dit_pass: a_in[1C] = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dit_pass: a_in[1E] = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = j1 + p10;	jp = j2 + p10;
			fprintf(stderr, "radix32_dit_pass: a_in[20] = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dit_pass: a_in[22] = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dit_pass: a_in[24] = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dit_pass: a_in[26] = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = jt + p04;	jp = jp + p04;
			fprintf(stderr, "radix32_dit_pass: a_in[28] = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dit_pass: a_in[2A] = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dit_pass: a_in[2C] = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dit_pass: a_in[2E] = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = j1 + p18;	jp = j2 + p18;
			fprintf(stderr, "radix32_dit_pass: a_in[30] = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dit_pass: a_in[32] = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dit_pass: a_in[34] = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dit_pass: a_in[36] = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = jt + p04;	jp = jp + p04;
			fprintf(stderr, "radix32_dit_pass: a_in[38] = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dit_pass: a_in[3A] = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dit_pass: a_in[3C] = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dit_pass: a_in[0E] = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);
		}
	#endif

	/*       gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...
			 We process the sincos data in bit-reversed order.	*/

	#ifdef COMPILER_TYPE_MSVC

	/*...Block 1: */
	#if 0
		add0 = &a[j1];
		add1 = add0+p08;
		add2 = add0+p10;
		add3 = add0+p18;
		__asm	mov eax, add0
		__asm	mov ebx, p05
		__asm	mov ecx, p06
		__asm	mov edx, p07
		__asm	mov	edi, p04		/* edi will store copy of p4 throughout */
		__asm	shl	ebx, 3			/* Pointer offset for floating doubles */
		__asm	shl	ecx, 3
		__asm	shl	edx, 3
		__asm	shl	edi, 3
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
		__asm	add eax, edi
		SSE2_RADIX8_DIT_0TWIDDLE_B(r00)
		/* Totals: 97 load/store [61 movaps, 36 implied], 54 add/subpd,  4 mulpd, 59 address-compute */
	#elif 1
		add0 = &a[j1];
		add1 = add0+p08;
		add2 = add0+p10;
		add3 = add0+p18;
		__asm	mov eax, add0
		__asm	mov ebx, p01
		__asm	mov ecx, p02
		__asm	mov edx, p03
		__asm	mov	edi, p04		/* edi will store copy of p4 throughout */
		__asm	shl	ebx, 3			/* Pointer offset for floating doubles */
		__asm	shl	ecx, 3
		__asm	shl	edx, 3
		__asm	shl	edi, 3
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
		SSE2_RADIX4_DIT_0TWIDDLE_B(r00)
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
		SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r08)
		SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r02,r04,r06,r08,r0A,r0C,r0E)
		/* Totals: 67 load/store [67 movaps,  0 implied], 76 add/subpd,  4 mulpd, 22 address-compute */
	#else
		add0 = &a[jt];
		add1 = add0+p01;
		add2 = add0+p02;
		add3 = add0+p03;
		add4 = add0+p04;
		add5 = add0+p05;
		add6 = add0+p06;
		add7 = add0+p07;

		SSE2_RADIX4_DIT_0TWIDDLE         (add0, add1, add2, add3, r00)
		SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r08)
		SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r00,r02,r04,r06,r08,r0A,r0C,r0E)
	#endif

	/*...Block 2:;	*/
	#if 0
		__asm	mov eax, add1
		__asm	mov ebx, p05
		__asm	mov ecx, p06
		__asm	mov edx, p07
		__asm	mov	edi, p04		/* edi will store copy of p4 throughout */
		__asm	shl	ebx, 3			/* Pointer offset for floating doubles */
		__asm	shl	ecx, 3
		__asm	shl	edx, 3
		__asm	shl	edi, 3
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
		__asm	add eax, edi
		SSE2_RADIX8_DIT_0TWIDDLE_B(r10)
	#elif 1
		__asm	mov eax, add1
		__asm	mov ebx, p01
		__asm	mov ecx, p02
		__asm	mov edx, p03
		__asm	mov	edi, p04		/* edi will store copy of p4 throughout */
		__asm	shl	ebx, 3			/* Pointer offset for floating doubles */
		__asm	shl	ecx, 3
		__asm	shl	edx, 3
		__asm	shl	edi, 3
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
		SSE2_RADIX4_DIT_0TWIDDLE_B(r10)
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
		SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r18)
		SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r10,r12,r14,r16,r18,r1A,r1C,r1E)
	#else
		add0 = &a[jt];
		add1 = add0+p01;
		add2 = add0+p02;
		add3 = add0+p03;
		add4 = add0+p04;
		add5 = add0+p05;
		add6 = add0+p06;
		add7 = add0+p07;

		SSE2_RADIX4_DIT_0TWIDDLE         (add0, add1, add2, add3, r10)
		SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r18)
		SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r10,r12,r14,r16,r18,r1A,r1C,r1E)
	#endif

	/*...Block 3:	*/
	#if 0
		__asm	mov eax, add2
		__asm	mov ebx, p05
		__asm	mov ecx, p06
		__asm	mov edx, p07
		__asm	mov	edi, p04		/* edi will store copy of p4 throughout */
		__asm	shl	ebx, 3			/* Pointer offset for floating doubles */
		__asm	shl	ecx, 3
		__asm	shl	edx, 3
		__asm	shl	edi, 3
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
		__asm	add eax, edi
		SSE2_RADIX8_DIT_0TWIDDLE_B(r20)
	#elif 1
		__asm	mov eax, add2
		__asm	mov ebx, p01
		__asm	mov ecx, p02
		__asm	mov edx, p03
		__asm	mov	edi, p04		/* edi will store copy of p4 throughout */
		__asm	shl	ebx, 3			/* Pointer offset for floating doubles */
		__asm	shl	ecx, 3
		__asm	shl	edx, 3
		__asm	shl	edi, 3
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
		SSE2_RADIX4_DIT_0TWIDDLE_B(r20)
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
		SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r28)
		SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r20,r22,r24,r26,r28,r2A,r2C,r2E)
	#else
		add0 = &a[jt];
		add1 = add0+p01;
		add2 = add0+p02;
		add3 = add0+p03;
		add4 = add0+p04;
		add5 = add0+p05;
		add6 = add0+p06;
		add7 = add0+p07;

		SSE2_RADIX4_DIT_0TWIDDLE         (add0, add1, add2, add3, r20)
		SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r28)
		SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r20,r22,r24,r26,r28,r2A,r2C,r2E)
	#endif

	/*...Block 4:	*/
	#if 0
		__asm	mov eax, add3
		__asm	mov ebx, p05
		__asm	mov ecx, p06
		__asm	mov edx, p07
		__asm	mov	edi, p04		/* edi will store copy of p4 throughout */
		__asm	shl	ebx, 3			/* Pointer offset for floating doubles */
		__asm	shl	ecx, 3
		__asm	shl	edx, 3
		__asm	shl	edi, 3
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
		__asm	add eax, edi
		SSE2_RADIX8_DIT_0TWIDDLE_B(r30)
	#elif 1
		__asm	mov eax, add3
		__asm	mov ebx, p01
		__asm	mov ecx, p02
		__asm	mov edx, p03
		__asm	mov	edi, p04		/* edi will store copy of p4 throughout */
		__asm	shl	ebx, 3			/* Pointer offset for floating doubles */
		__asm	shl	ecx, 3
		__asm	shl	edx, 3
		__asm	shl	edi, 3
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
		SSE2_RADIX4_DIT_0TWIDDLE_B(r30)
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
		SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO_B(r38)
		SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r30,r32,r34,r36,r38,r3A,r3C,r3E)
	#else
		add0 = &a[jt];
		add1 = add0+p01;
		add2 = add0+p02;
		add3 = add0+p03;
		add4 = add0+p04;
		add5 = add0+p05;
		add6 = add0+p06;
		add7 = add0+p07;

		SSE2_RADIX4_DIT_0TWIDDLE         (add0, add1, add2, add3, r30)
		SSE2_RADIX4_DIT_0TWIDDLE_2NDOFTWO(add4, add5, add6, add7, r38)
		SSE2_RADIX8_DIT_COMBINE_RAD4_SUBS(r30,r32,r34,r36,r38,r3A,r3C,r3E)
	#endif

	/*...and now do eight radix-4 transforms, including the internal twiddle factors:
		1, exp(-i* 1*twopi/32) =       ( c32_1,-s32_1), exp(-i* 2*twopi/32) =       ( c    ,-s    ), exp(-i* 3*twopi/32) =       ( c32_3,-s32_3) (for inputs to transform block 2),
		1, exp(-i* 2*twopi/32) =       ( c    ,-s    ), exp(-i* 4*twopi/32) = ISRT2*( 1    ,-1    ), exp(-i* 6*twopi/32) =       ( s    ,-c    ) (for inputs to transform block 3),
		1, exp(-i* 3*twopi/32) =       ( c32_3,-s32_3), exp(-i* 6*twopi/32) =       ( s    ,-c    ), exp(-i* 9*twopi/32) =       (-s32_1,-c32_1) (for inputs to transform block 4),
		1, exp(-i* 4*twopi/32) = ISRT2*( 1    ,-1    ), exp(-i* 8*twopi/32) =       ( 0    ,-1    ), exp(-i*12*twopi/32) = ISRT2*(-1    ,-1    ) (for inputs to transform block 5),
		1, exp(-i* 5*twopi/32) =       ( s32_3,-c32_3), exp(-i*10*twopi/32) =       (-s    ,-c    ), exp(-i*15*twopi/32) =       (-c32_1,-s32_1) (for inputs to transform block 6),
		1, exp(-i* 6*twopi/32) =       ( s    ,-c    ), exp(-i*12*twopi/32) = ISRT2*(-1    ,-1    ), exp(-i*18*twopi/32) =       (-c    , s    ) (for inputs to transform block 7),
		1, exp(-i* 7*twopi/32) =       ( s32_1,-c32_1), exp(-i*14*twopi/32) =       (-c    ,-s    ), exp(-i*21*twopi/32) =       (-s32_3, c32_3) (for inputs to transform block 8),
		 and only the last 3 inputs to each of the radix-4 transforms 2 through 8 are multiplied by non-unity twiddles.	*/

	/*...Block 1: t00,t10,t20,t30	*/
	#if 0
		add0 = &a[j1];
		add1 = add0+p08;
		add2 = add0+p10;
		add3 = add0+p18;

		__asm	mov	eax, r00
		__asm	mov	edx, eax
		__asm	add	edx, 0x100	/* r10 */

		__asm	movaps	xmm2,[edx      ]	/* t10 */				__asm	movaps	xmm4,[edx+0x200]	/* t30 */
		__asm	movaps	xmm3,[edx+0x010]	/* t11 */				__asm	movaps	xmm5,[edx+0x210]	/* t31 */
		__asm	movaps	xmm0,[eax      ]	/* t00 */				__asm	movaps	xmm6,[eax+0x200]	/* t20 */
		__asm	movaps	xmm1,[eax+0x010]	/* t01 */				__asm	movaps	xmm7,[eax+0x210]	/* t21 */

		__asm	subpd	xmm0,xmm2			/*~t10=t00-t10*/		__asm	subpd	xmm6,xmm4			/*~t30=t20-t30*/
		__asm	subpd	xmm1,xmm3			/*~t11=t01-t11*/		__asm	subpd	xmm7,xmm5			/*~t31=t21-t31*/
		__asm	addpd	xmm2,xmm2			/*       2*t10*/		__asm	addpd	xmm4,xmm4			/*       2*t30*/
		__asm	addpd	xmm3,xmm3			/*       2*t11*/		__asm	addpd	xmm5,xmm5			/*       2*t31*/
		__asm	addpd	xmm2,xmm0			/*~t00=t00+t10*/		__asm	addpd	xmm4,xmm6			/*~t20=t20+t30*/
		__asm	addpd	xmm3,xmm1			/*~t01=t01+t11*/		__asm	addpd	xmm5,xmm7			/*~t21=t21+t31*/

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF(add0, add1, add2, add3, c00)
	#else
		__asm	mov eax, add0
		__asm	mov ebx, p08
		__asm	mov ecx, r00
		__asm	mov edx, r10
		__asm	mov edi, p10	/* edi will store copy of p10 throughout */
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	shl	edi, 3
		__asm	add ebx, eax	/* add1 = add0+p8 */

		__asm	movaps	xmm2,[edx      ]	/* t10 */				__asm	movaps	xmm4,[edx+0x200]	/* t30 */
		__asm	movaps	xmm3,[edx+0x010]	/* t11 */				__asm	movaps	xmm5,[edx+0x210]	/* t31 */
		__asm	movaps	xmm0,[ecx      ]	/* t00 */				__asm	movaps	xmm6,[ecx+0x200]	/* t20 */
		__asm	movaps	xmm1,[ecx+0x010]	/* t01 */				__asm	movaps	xmm7,[ecx+0x210]	/* t21 */

		__asm	subpd	xmm0,xmm2			/*~t10=t00-t10*/		__asm	subpd	xmm6,xmm4			/*~t30=t20-t30*/
		__asm	subpd	xmm1,xmm3			/*~t11=t01-t11*/		__asm	subpd	xmm7,xmm5			/*~t31=t21-t31*/
		__asm	addpd	xmm2,xmm2			/*       2*t10*/		__asm	addpd	xmm4,xmm4			/*       2*t30*/
		__asm	addpd	xmm3,xmm3			/*       2*t11*/		__asm	addpd	xmm5,xmm5			/*       2*t31*/
		__asm	addpd	xmm2,xmm0			/*~t00=t00+t10*/		__asm	addpd	xmm4,xmm6			/*~t20=t20+t30*/
		__asm	addpd	xmm3,xmm1			/*~t01=t01+t11*/		__asm	addpd	xmm5,xmm7			/*~t21=t21+t31*/

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c00)
	#endif

	/*...Block 2: t02,t12,t22,t32	*/
		__asm	mov eax, add0
		__asm	mov esi, p01
		__asm	mov ebx, p08
		__asm	shl	esi, 3
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	add eax, esi	/* add0 = &a[j1+p1] */
		__asm	add ebx, eax	/* add1 = add0+p8 */
		__asm	mov ecx, r02
		__asm	mov edx, r12
		__asm	add	ecx, 0x200
		__asm	add	edx, 0x200
		__asm	mov	esi, cc1

		__asm	movaps	xmm4,[ecx      ]	/* t22 */
		__asm	movaps	xmm5,[ecx+0x010]	/* t23 */
		__asm	movaps	xmm2,[esi      ]	/* c32_1 */
		__asm	movaps	xmm3,[esi+0x010]	/* s32_1 */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy t22 */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy t23 */

		__asm	mulpd	xmm4,xmm2		/* t22*c32_1 */
		__asm	mulpd	xmm5,xmm2		/* t23*c32_1 */
		__asm	mulpd	xmm6,xmm3		/* t22*s32_1 */	__asm	movaps	xmm0,[edx      ]	/* t32 */
		__asm	mulpd	xmm7,xmm3		/* t23*s32_1 */	__asm	movaps	xmm1,[edx+0x010]	/* t33 */
														__asm	add	esi, 0x20
														__asm	movaps	xmm2,[esi      ]	/* c32_3 */
														__asm	movaps	xmm3,[esi+0x010]	/* s32_3 */
		__asm	subpd	xmm5,xmm6	/* xmm5 <-~t23 */	__asm	movaps	xmm6,xmm0	/* xmm6 <- cpy t32 */
		__asm	addpd	xmm4,xmm7	/* xmm4 <-~t22 */	__asm	movaps	xmm7,xmm1	/* xmm7 <- cpy t33 */

														__asm	mulpd	xmm0,xmm2		/* t32*c32_3 */
														__asm	mulpd	xmm1,xmm2		/* t33*c32_3 */
														__asm	mulpd	xmm6,xmm3		/* t32*s32_3 */
														__asm	mulpd	xmm7,xmm3		/* t33*s32_3 */
														__asm	subpd	xmm1,xmm6	/* xmm1 <- it */
														__asm	addpd	xmm0,xmm7	/* xmm0 <- rt */

		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy~t23*/
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy~t22*/

		__asm	addpd	xmm4,xmm0	/* ~t22 <- t22+rt */
		__asm	addpd	xmm5,xmm1	/* ~t23 <- t23+it */
		__asm	subpd	xmm6,xmm0	/* ~t32 <- t22-rt */
		__asm	subpd	xmm7,xmm1	/* ~t33 <- t23-it */
	/*
		rt =t12*c + t13*s;			it =t13*c - t12*s;
		t12=t02-rt;					t02=t02+rt;
		t13=t03-it;					t03=t03+it;
	*/
		__asm	sub	ecx, 0x200
		__asm	sub	edx, 0x200
		__asm	mov	esi, cc0
		__asm	movaps	xmm1,[edx      ]	/* t12 */
		__asm	movaps	xmm3,[edx+0x010]	/* t13 */
		__asm	movaps	xmm2,[esi+0x010]	/* s */
		__asm	movaps	xmm0,xmm1			/* cpy t12 */
		__asm	mulpd	xmm1,xmm2		/* t12*s */
		__asm	mulpd	xmm2,xmm3		/* t13*s */
		__asm	mulpd	xmm0,[esi]		/* t12*c */
		__asm	mulpd	xmm3,[esi]		/* t13*c */

		__asm	addpd	xmm2,xmm0	/* rt =t12*c + t13*s */
		__asm	subpd	xmm3,xmm1	/* it =t13*c - t12*s */

		__asm	movaps	xmm0,[ecx      ]	/* t02 */
		__asm	movaps	xmm1,[ecx+0x010]	/* t03 */
		__asm	subpd	xmm0,xmm2	/*~t12 <- t02- rt */
		__asm	subpd	xmm1,xmm3	/*~t13 <- t03- it */
		__asm	addpd	xmm2,xmm2	/*          2* rt */
		__asm	addpd	xmm3,xmm3	/*          2* it */
		__asm	addpd	xmm2,xmm0	/*~t02 <- t02+ rt */
		__asm	addpd	xmm3,xmm1	/*~t03 <- t03+ it */

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c01)

	/*...Block 3: t04,t14,t24,t34	*/
		__asm	mov eax, add0
		__asm	mov esi, p02
		__asm	mov ebx, p08
		__asm	shl	esi, 3
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	add eax, esi	/* add0 = &a[j1+p2] */
		__asm	add ebx, eax	/* add1 = add0+p8 */
		__asm	mov ecx, r04
		__asm	mov edx, r14
		__asm	add	ecx, 0x200
		__asm	add	edx, 0x200
		__asm	mov	esi, cc0

		__asm	movaps	xmm4,[ecx      ]	/* t24 */
		__asm	movaps	xmm5,[ecx+0x010]	/* t25 */
		__asm	movaps	xmm2,[esi      ]	/* c */
		__asm	movaps	xmm3,[esi+0x010]	/* s */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy t24 */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy t25 */

		__asm	mulpd	xmm4,xmm2		/* t24*c */
		__asm	mulpd	xmm5,xmm2		/* t25*c */
		__asm	mulpd	xmm6,xmm3		/* t24*s */		__asm	movaps	xmm0,[edx      ]	/* t34 */
		__asm	mulpd	xmm7,xmm3		/* t25*s */		__asm	movaps	xmm1,[edx+0x010]	/* t35 */
		__asm	subpd	xmm5,xmm6	/* xmm1 <-~t25 */	__asm	movaps	xmm6,xmm0	/* xmm6 <- cpy t34 */
		__asm	addpd	xmm4,xmm7	/* xmm0 <-~t24 */	__asm	movaps	xmm7,xmm1	/* xmm7 <- cpy t35 */

														__asm	mulpd	xmm0,xmm3		/* t34*s */
														__asm	mulpd	xmm1,xmm3		/* t35*s */
														__asm	mulpd	xmm6,xmm2		/* t34*c */
														__asm	mulpd	xmm7,xmm2		/* t35*c */
														__asm	subpd	xmm1,xmm6	/* xmm5 <- it */
														__asm	addpd	xmm0,xmm7	/* xmm4 <- rt */

		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy~t25*/
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy~t24*/

		__asm	addpd	xmm4,xmm0	/* ~t24 <- t24+rt */
		__asm	addpd	xmm5,xmm1	/* ~t25 <- t25+it */
		__asm	subpd	xmm6,xmm0	/* ~t34 <- t24-rt */
		__asm	subpd	xmm7,xmm1	/* ~t35 <- t25-it */

		__asm	sub	ecx, 0x200
		__asm	sub	edx, 0x200
		__asm	mov	esi, isrt2
		__asm	movaps	xmm2,[edx      ]	/* t14 */
		__asm	movaps	xmm3,[edx+0x010]	/* t15 */
		__asm	movaps	xmm1,[esi]	/* isrt2 */
		__asm	movaps	xmm0,xmm3	/* cpy t15 */
		__asm	subpd	xmm3,xmm2	/*~t15=t15-t14 */
		__asm	addpd	xmm2,xmm0	/*~t14=t14+t15 */
		__asm	mulpd	xmm2,xmm1	/* rt */
		__asm	mulpd	xmm3,xmm1	/* it */

		__asm	movaps	xmm0,[ecx      ]	/* t04 */
		__asm	movaps	xmm1,[ecx+0x010]	/* t05 */
		__asm	subpd	xmm0,xmm2	/*~t14 <- t04- rt */
		__asm	subpd	xmm1,xmm3	/*~t15 <- t05- it */
		__asm	addpd	xmm2,xmm2	/*          2* rt */
		__asm	addpd	xmm3,xmm3	/*          2* it */
		__asm	addpd	xmm2,xmm0	/*~t04 <- t04+ rt */
		__asm	addpd	xmm3,xmm1	/*~t05 <- t05+ it */

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c02)

	/*...Block 4: t06,t16,t26,t36	*/
		__asm	mov eax, add0
		__asm	mov esi, p03
		__asm	mov ebx, p08
		__asm	shl	esi, 3
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	add eax, esi	/* add0 = &a[j1+p3] */
		__asm	add ebx, eax	/* add1 = add0+p8 */
		__asm	mov ecx, r06
		__asm	mov edx, r16
		__asm	add	ecx, 0x200
		__asm	add	edx, 0x200
		__asm	mov	esi, cc3

		__asm	movaps	xmm4,[ecx      ]	/* t26 */
		__asm	movaps	xmm5,[ecx+0x010]	/* t27 */
		__asm	movaps	xmm2,[esi      ]	/* c32_3 */
		__asm	movaps	xmm3,[esi+0x010]	/* s32_3 */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy t26 */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy t27 */

		__asm	mulpd	xmm4,xmm2		/* t26*c32_3 */
		__asm	mulpd	xmm5,xmm2		/* t27*c32_3 */
		__asm	mulpd	xmm6,xmm3		/* t26*s32_3 */	__asm	movaps	xmm0,[edx      ]	/* t36 */
		__asm	mulpd	xmm7,xmm3		/* t27*s32_3 */	__asm	movaps	xmm1,[edx+0x010]	/* t37 */
														__asm	sub	esi, 0x20
														__asm	movaps	xmm3,[esi      ]	/* c32_1 */
														__asm	movaps	xmm2,[esi+0x010]	/* s32_1 */
		__asm	subpd	xmm5,xmm6	/* xmm5 <-~t27 */	__asm	movaps	xmm6,xmm0	/* xmm6 <- cpy t36 */
		__asm	addpd	xmm4,xmm7	/* xmm4 <-~t26 */	__asm	movaps	xmm7,xmm1	/* xmm7 <- cpy t37 */

														__asm	mulpd	xmm0,xmm2		/* t36*s32_1 */
														__asm	mulpd	xmm1,xmm2		/* t37*s32_1 */
														__asm	mulpd	xmm6,xmm3		/* t36*c32_1 */
														__asm	mulpd	xmm7,xmm3		/* t37*c32_1 */
														__asm	addpd	xmm1,xmm6	/* xmm1 <- it */
														__asm	subpd	xmm0,xmm7	/* xmm0 <- rt */

		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy~t27*/
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy~t26*/

		__asm	addpd	xmm6,xmm0	/* ~t36 <- t26+rt */
		__asm	addpd	xmm7,xmm1	/* ~t37 <- t27+it */
		__asm	subpd	xmm4,xmm0	/* ~t26 <- t26-rt */
		__asm	subpd	xmm5,xmm1	/* ~t27 <- t27-it */
	/*
		rt =t16*s + t17*c;			it =t17*s - t16*c;
		t16=t06-rt;					t06 =t06+rt;
		t17=t07-it;					t07 =t07+it;
	*/
		__asm	sub	ecx, 0x200
		__asm	sub	edx, 0x200
		__asm	mov	esi, cc0
		__asm	movaps	xmm2,[edx      ]	/* t16 */
		__asm	movaps	xmm0,[edx+0x010]	/* t17 */
		__asm	movaps	xmm3,[esi+0x010]	/* s */
		__asm	movaps	xmm1,xmm2			/* cpy t16 */
		__asm	mulpd	xmm2,xmm3		/* t16*s */
		__asm	mulpd	xmm3,xmm0		/* s*t17 */
		__asm	mulpd	xmm1,[esi]		/* t16*c */
		__asm	mulpd	xmm0,[esi]		/* t17*c */

		__asm	addpd	xmm2,xmm0	/* rt =t16*s - t17*c */
		__asm	subpd	xmm3,xmm1	/* it =t17*s + t16*c */

		__asm	movaps	xmm0,[ecx      ]	/* t06 */
		__asm	movaps	xmm1,[ecx+0x010]	/* t07 */
		__asm	subpd	xmm0,xmm2	/*~t16 <- t06- rt */
		__asm	subpd	xmm1,xmm3	/*~t17 <- t07- it */
		__asm	addpd	xmm2,xmm2	/*          2* rt */
		__asm	addpd	xmm3,xmm3	/*          2* it */
		__asm	addpd	xmm2,xmm0	/*~t06 <- t06+ rt */
		__asm	addpd	xmm3,xmm1	/*~t07 <- t07+ it */

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c03)

	/*...Block 5: t08,t18,t28,t38	*/
		__asm	mov eax, add0
		__asm	mov esi, p04
		__asm	mov ebx, p08
		__asm	shl	esi, 3
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	add eax, esi	/* add0 = &a[j1+p4] */
		__asm	add ebx, eax	/* add1 = add0+p8 */
		__asm	mov ecx, r08
		__asm	mov edx, r18
		__asm	add	ecx, 0x200
		__asm	add	edx, 0x200
		__asm	mov	esi, isrt2
		__asm	movaps	xmm2,[esi]	/* isrt2 */

		__asm	movaps	xmm4,[ecx      ]	/* t28 */
		__asm	movaps	xmm5,[ecx+0x010]	/* t29 */
		__asm	movaps	xmm6,[edx      ]	/* t38 */
		__asm	movaps	xmm7,[edx+0x010]	/* t39 */
		__asm	sub	ecx, 0x200
		__asm	sub	edx, 0x200
		__asm	mulpd	xmm4,xmm2
		__asm	mulpd	xmm5,xmm2
		__asm	mulpd	xmm6,xmm2
		__asm	mulpd	xmm7,xmm2

		__asm	subpd	xmm5,xmm4			/*~t29=t29-t28*/		__asm	movaps	xmm0,[ecx      ]	/* t08 */
		__asm	subpd	xmm6,xmm7			/* rt =t38-t39*/		__asm	movaps	xmm2,[edx+0x010]	/* t19 */
		__asm	addpd	xmm4,xmm4			/*       2*t28*/		__asm	movaps	xmm3,[ecx+0x010]	/* t09 */
		__asm	addpd	xmm7,xmm7			/*       2*t39*/		__asm	movaps	xmm1,[edx      ]	/* t18 */
		__asm	addpd	xmm4,xmm5			/*~t28=t28+t29*/
		__asm	addpd	xmm7,xmm6			/* it =t39+t38*/

		__asm	subpd	xmm4,xmm6			/*~t28=t28-rt */		__asm	subpd	xmm0,xmm2			/*~t18=t08-t19*/
		__asm	subpd	xmm5,xmm7			/*~t29=t29-it */		__asm	subpd	xmm3,xmm1			/*~t09=t09-t18*/
		__asm	addpd	xmm6,xmm6			/*       2*rt */		__asm	addpd	xmm2,xmm2			/*       2*t08*/
		__asm	addpd	xmm7,xmm7			/*       2*it */		__asm	addpd	xmm1,xmm1			/*       2*t09*/
		__asm	addpd	xmm6,xmm4			/*~t38=t28+rt */		__asm	addpd	xmm2,xmm0			/*~t08=t19+t08*/
		__asm	addpd	xmm7,xmm5			/*~t39=t29+it */		__asm	addpd	xmm1,xmm3			/*~t19=t18+t09*/

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c04)

	/*...Block 6: t0A,t1A,t2A,t3A	*/
		__asm	mov eax, add0
		__asm	mov esi, p05
		__asm	mov ebx, p08
		__asm	shl	esi, 3
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	add eax, esi	/* add0 = &a[j1+p5] */
		__asm	add ebx, eax	/* add1 = add0+p8 */
		__asm	mov ecx, r0A
		__asm	mov edx, r1A
		__asm	add	ecx, 0x200
		__asm	add	edx, 0x200
		__asm	mov	esi, cc3

		__asm	movaps	xmm4,[ecx      ]	/* t2A */
		__asm	movaps	xmm5,[ecx+0x010]	/* t2B */
		__asm	movaps	xmm3,[esi      ]	/* c32_3 */
		__asm	movaps	xmm2,[esi+0x010]	/* s32_3 */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy t2A */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy t2B */

		__asm	mulpd	xmm4,xmm2		/* t2A*s32_3 */
		__asm	mulpd	xmm5,xmm2		/* t2B*s32_3 */
		__asm	mulpd	xmm6,xmm3		/* t2A*c32_3 */	__asm	movaps	xmm0,[edx      ]	/* t3A */
		__asm	mulpd	xmm7,xmm3		/* t2B*c32_3 */	__asm	movaps	xmm1,[edx+0x010]	/* t3B */
														__asm	sub	esi, 0x20
														__asm	movaps	xmm2,[esi      ]	/* c32_1 */
														__asm	movaps	xmm3,[esi+0x010]	/* s32_1 */
		__asm	subpd	xmm5,xmm6	/* xmm5 <-~t2B */	__asm	movaps	xmm6,xmm0	/* xmm6 <- cpy t3A */
		__asm	addpd	xmm4,xmm7	/* xmm4 <-~t2A */	__asm	movaps	xmm7,xmm1	/* xmm7 <- cpy t3B */

														__asm	mulpd	xmm0,xmm2		/* t3A*c32_1 */
														__asm	mulpd	xmm1,xmm2		/* t3B*c32_1 */
														__asm	mulpd	xmm6,xmm3		/* t3A*s32_1 */
														__asm	mulpd	xmm7,xmm3		/* t3B*s32_1 */
														__asm	addpd	xmm1,xmm6	/* xmm1 <- it */
														__asm	subpd	xmm0,xmm7	/* xmm0 <- rt */

		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy~t2B*/
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy~t2A*/

		__asm	addpd	xmm6,xmm0	/* ~t3A <- t2A+rt */
		__asm	addpd	xmm7,xmm1	/* ~t3B <- t2B+it */
		__asm	subpd	xmm4,xmm0	/* ~t2A <- t2A-rt */
		__asm	subpd	xmm5,xmm1	/* ~t2B <- t2B-it */
	/*
		rt =t1A*s - t1B*c;		it =t1B*s + t1A*c;
		t1A=t0A+rt;				t0A =t0A-rt;
		t1B=t0B+it;				t0B =t0B-it;
	*/
		__asm	sub	ecx, 0x200
		__asm	sub	edx, 0x200
		__asm	mov	esi, cc0
		__asm	movaps	xmm0,[edx      ]	/* t1A */
		__asm	movaps	xmm2,[edx+0x010]	/* t1B */
		__asm	movaps	xmm1,[esi+0x010]	/* s */
		__asm	movaps	xmm3,xmm0			/* cpy t1A */
		__asm	mulpd	xmm0,xmm1		/* t1A*s */
		__asm	mulpd	xmm1,xmm2		/* s*t1B */
		__asm	mulpd	xmm3,[esi]		/* t1A*c */
		__asm	mulpd	xmm2,[esi]		/* t1B*c */

		__asm	subpd	xmm0,xmm2	/* rt =t1A*s - t1B*c */
		__asm	addpd	xmm1,xmm3	/* it =t1B*s + t1A*c */

		__asm	movaps	xmm2,[ecx      ]	/* t0A */
		__asm	movaps	xmm3,[ecx+0x010]	/* t0B */
		__asm	subpd	xmm2,xmm0	/*~t0A <- t0A- rt */
		__asm	subpd	xmm3,xmm1	/*~t0B <- t0B- it */
		__asm	addpd	xmm0,xmm0	/*          2* rt */
		__asm	addpd	xmm1,xmm1	/*          2* it */
		__asm	addpd	xmm0,xmm2	/*~t1A <- t0A+ rt */
		__asm	addpd	xmm1,xmm3	/*~t1B <- t0B+ it */

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c05)

	/*...Block 7: t0C,t1C,t2C,t3C	*/
		__asm	mov eax, add0
		__asm	mov esi, p06
		__asm	mov ebx, p08
		__asm	shl	esi, 3
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	add eax, esi	/* add0 = &a[j1+p5] */
		__asm	add ebx, eax	/* add1 = add0+p8 */
		__asm	mov ecx, r0C
		__asm	mov edx, r1C
		__asm	add	ecx, 0x200
		__asm	add	edx, 0x200
		__asm	mov	esi, cc0

		__asm	movaps	xmm4,[ecx      ]	/* t2C */
		__asm	movaps	xmm5,[ecx+0x010]	/* t2D */
		__asm	movaps	xmm3,[esi      ]	/* c */
		__asm	movaps	xmm2,[esi+0x010]	/* s */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy t2C */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy t2D */

		__asm	mulpd	xmm4,xmm2		/* t2C*s */
		__asm	mulpd	xmm5,xmm2		/* t2D*s */
		__asm	mulpd	xmm6,xmm3		/* t2C*c */		__asm	movaps	xmm0,[edx      ]	/* t3C */
		__asm	mulpd	xmm7,xmm3		/* t2D*c */		__asm	movaps	xmm1,[edx+0x010]	/* t3D */
		__asm	subpd	xmm5,xmm6	/* xmm5 <-~t2D */	__asm	movaps	xmm6,xmm0	/* xmm6 <- cpy t3C */
		__asm	addpd	xmm4,xmm7	/* xmm4 <-~t2C */	__asm	movaps	xmm7,xmm1	/* xmm7 <- cpy t3D */

														__asm	mulpd	xmm0,xmm3		/* t3C*c */
														__asm	mulpd	xmm1,xmm3		/* t3D*c */
														__asm	mulpd	xmm6,xmm2		/* t3C*s */
														__asm	mulpd	xmm7,xmm2		/* t3D*s */
														__asm	subpd	xmm1,xmm6	/* xmm1 <- it */
														__asm	addpd	xmm0,xmm7	/* xmm0 <- rt */

		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy~t2D*/
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy~t2C*/

		__asm	addpd	xmm6,xmm0	/* ~t3C <- t2C+rt */
		__asm	addpd	xmm7,xmm1	/* ~t3D <- t2D+it */
		__asm	subpd	xmm4,xmm0	/* ~t2C <- t2C-rt */
		__asm	subpd	xmm5,xmm1	/* ~t2D <- t2D-it */

		__asm	sub	ecx, 0x200
		__asm	sub	edx, 0x200
		__asm	mov	esi, isrt2
		__asm	movaps	xmm0,[edx      ]	/* t1C */
		__asm	movaps	xmm1,[edx+0x010]	/* t1D */
		__asm	movaps	xmm3,[esi]	/* isrt2 */
		__asm	movaps	xmm2,xmm0	/* cpy t1C */
		__asm	subpd	xmm0,xmm1	/*~t1C=t1C-t1D */
		__asm	addpd	xmm1,xmm2	/*~t1D=t1D+t1C */
		__asm	mulpd	xmm0,xmm3	/* it */
		__asm	mulpd	xmm1,xmm3	/* rt */

		__asm	movaps	xmm2,[ecx      ]	/* t0C */
		__asm	movaps	xmm3,[ecx+0x010]	/* t0D */
		__asm	subpd	xmm2,xmm0	/*~t0C <- t0C- rt */
		__asm	subpd	xmm3,xmm1	/*~t0D <- t0D- it */
		__asm	addpd	xmm0,xmm0	/*          2* rt */
		__asm	addpd	xmm1,xmm1	/*          2* it */
		__asm	addpd	xmm0,xmm2	/*~t1C <- t0C+ rt */
		__asm	addpd	xmm1,xmm3	/*~t1D <- t0D+ it */

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c06)

	/*...Block 8: t0E,t1E,t2E,t3E	*/
		__asm	mov eax, add0
		__asm	mov esi, p07
		__asm	mov ebx, p08
		__asm	shl	esi, 3
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	add eax, esi	/* add0 = &a[j1+p5] */
		__asm	add ebx, eax	/* add1 = add0+p8 */
		__asm	mov ecx, r0E
		__asm	mov edx, r1E
		__asm	add	ecx, 0x200
		__asm	add	edx, 0x200
		__asm	mov	esi, cc1

		__asm	movaps	xmm4,[ecx      ]	/* t2E */
		__asm	movaps	xmm5,[ecx+0x010]	/* t2F */
		__asm	movaps	xmm3,[esi      ]	/* c32_1 */
		__asm	movaps	xmm2,[esi+0x010]	/* s32_1 */
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy t2E */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy t2F */

		__asm	mulpd	xmm4,xmm2		/* t2E*s32_1 */
		__asm	mulpd	xmm5,xmm2		/* t2F*s32_1 */
		__asm	mulpd	xmm6,xmm3		/* t2E*c32_1 */	__asm	movaps	xmm0,[edx      ]	/* t3E */
		__asm	mulpd	xmm7,xmm3		/* t2F*c32_1 */	__asm	movaps	xmm1,[edx+0x010]	/* t3F */
														__asm	add	esi, 0x20
														__asm	movaps	xmm3,[esi      ]	/* c32_3 */
														__asm	movaps	xmm2,[esi+0x010]	/* s32_3 */
		__asm	subpd	xmm5,xmm6	/* xmm5 <-~t2F */	__asm	movaps	xmm6,xmm0	/* xmm6 <- cpy t3E */
		__asm	addpd	xmm4,xmm7	/* xmm4 <-~t2E */	__asm	movaps	xmm7,xmm1	/* xmm7 <- cpy t3F */

														__asm	mulpd	xmm0,xmm2		/* t3E*s32_3 */
														__asm	mulpd	xmm1,xmm2		/* t3F*s32_3 */
														__asm	mulpd	xmm6,xmm3		/* t3E*c32_3 */
														__asm	mulpd	xmm7,xmm3		/* t3F*c32_3 */
														__asm	subpd	xmm1,xmm6	/* xmm1 <- it */
														__asm	addpd	xmm0,xmm7	/* xmm0 <- rt */

		__asm	movaps	xmm7,xmm5	/* xmm7 <- cpy~t2F*/
		__asm	movaps	xmm6,xmm4	/* xmm6 <- cpy~t2E*/

		__asm	subpd	xmm4,xmm0	/* ~t2E <- t2E-rt */
		__asm	subpd	xmm5,xmm1	/* ~t2F <- t2F-it */
		__asm	addpd	xmm6,xmm0	/* ~t3E <- t2E+rt */
		__asm	addpd	xmm7,xmm1	/* ~t3F <- t2F+it */
	/*
		rt =t1E*c - t1F*s;		it =t1F*c + t1E*s;
		t1E=t0E+rt;				t0E =t0E-rt;
		t1F=t0F+it;				t0F =t0F-it;
	*/
		__asm	sub	ecx, 0x200
		__asm	sub	edx, 0x200
		__asm	mov	esi, cc0
		__asm	movaps	xmm3,[edx      ]	/* t1E */
		__asm	movaps	xmm1,[edx+0x010]	/* t1F */
		__asm	movaps	xmm2,[esi+0x010]	/* s */
		__asm	movaps	xmm0,xmm3			/* cpy t1E */
		__asm	mulpd	xmm3,xmm2		/* t1E*s */
		__asm	mulpd	xmm2,xmm1		/* t1F*s */
		__asm	mulpd	xmm0,[esi]		/* t1E*c */
		__asm	mulpd	xmm1,[esi]		/* t1F*c */

		__asm	subpd	xmm0,xmm2	/* rt =t1E*c - t1F*s */
		__asm	addpd	xmm1,xmm3	/* it =t1F*c + t1E*s */

		__asm	movaps	xmm2,[ecx      ]	/* t0E */
		__asm	movaps	xmm3,[ecx+0x010]	/* t0F */
		__asm	subpd	xmm2,xmm0	/*~t0E <- t0E- rt */
		__asm	subpd	xmm3,xmm1	/*~t0F <- t0F- it */
		__asm	addpd	xmm0,xmm0	/*          2* rt */
		__asm	addpd	xmm1,xmm1	/*          2* it */
		__asm	addpd	xmm0,xmm2	/*~t1E <- t0E+ rt */
		__asm	addpd	xmm1,xmm3	/*~t1F <- t0F+ it */

		SSE2_RADIX4_DIT_4TWIDDLE_2ND_HALF_B(c07)

	#else	/* GCC-style inline ASM: */

		add0 = &a[j1];
	  #if (OS_BITS == 32) || !USE_64BIT_ASM_STYLE	// 2nd of these is a toggle defined and set in radix32_dif_dit_pass_gcc64.h
		SSE2_RADIX32_DIT_TWIDDLE(add0,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,r00,r02,r04,r06,r08,r0A,r0C,r0E,r10,r12,r14,r16,r18,r1A,r1C,r1E,r20,r30,isrt2)
	  #else
		SSE2_RADIX32_DIT_TWIDDLE(add0,p01,p02,p03,p04,p05,p06,p07,p08,p10,p18,r00,r10,r20,r30,isrt2)
	  #endif

	#endif

	#ifdef DEBUG_SSE2
		if(s01->re != 0.0)
		{
			jt = j1; jp = j2;
			fprintf(stderr, "radix32_dit_pass: a_out[00] = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dit_pass: a_out[02] = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dit_pass: a_out[04] = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dit_pass: a_out[06] = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = jt + p04;	jp = jp + p04;
			fprintf(stderr, "radix32_dit_pass: a_out[08] = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dit_pass: a_out[0A] = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dit_pass: a_out[0C] = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dit_pass: a_out[0E] = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = j1 + p08;	jp = j2 + p08;
			fprintf(stderr, "radix32_dit_pass: a_out[10] = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dit_pass: a_out[12] = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dit_pass: a_out[14] = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dit_pass: a_out[16] = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = jt + p04;	jp = jp + p04;
			fprintf(stderr, "radix32_dit_pass: a_out[18] = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dit_pass: a_out[1A] = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dit_pass: a_out[1C] = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dit_pass: a_out[1E] = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = j1 + p10;	jp = j2 + p10;
			fprintf(stderr, "radix32_dit_pass: a_out[20] = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dit_pass: a_out[22] = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dit_pass: a_out[24] = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dit_pass: a_out[26] = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = jt + p04;	jp = jp + p04;
			fprintf(stderr, "radix32_dit_pass: a_out[28] = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dit_pass: a_out[2A] = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dit_pass: a_out[2C] = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dit_pass: a_out[2E] = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = j1 + p18;	jp = j2 + p18;
			fprintf(stderr, "radix32_dit_pass: a_out[30] = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dit_pass: a_out[32] = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dit_pass: a_out[34] = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dit_pass: a_out[36] = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);	jt = jt + p04;	jp = jp + p04;
			fprintf(stderr, "radix32_dit_pass: a_out[38] = %20.5f, %20.5f\n", a[jt    ], a[jp    ]);
			fprintf(stderr, "radix32_dit_pass: a_out[3A] = %20.5f, %20.5f\n", a[jt+p01], a[jp+p01]);
			fprintf(stderr, "radix32_dit_pass: a_out[3C] = %20.5f, %20.5f\n", a[jt+p02], a[jp+p02]);
			fprintf(stderr, "radix32_dit_pass: a_out[0E] = %20.5f, %20.5f\n", a[jt+p03], a[jp+p03]);
			exit(0);
		}
	#endif

#else	/* USE_SSE2 */

	/*       gather the needed data (32 64-bit complex, i.e. 64 64-bit reals) and do the first set of four length-8 transforms...
			 We process the sincos data in bit-reversed order.	*/
	#ifdef PFETCH_AGGRESSIVE
		addr = &a[j1];
	#elif PFETCH
		prefetch_offset = ((j >> 1) & 7)*p04 + 4;	/* Cycle among p00, p04, p08, p0C, p10, p14, p18 and p1C. */
	#endif

	/*...Block 1: */

		jt = j1;
		jp = j2;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr;
		prefetch_p_doubles(addr);
	#elif PFETCH
		addr = &a[jt];
		addp = addr+prefetch_offset;
		prefetch_p_doubles(addp);
	#endif
			t00=a[jt    ];	t01=a[jp    ];
			rt =a[jt+p01];	it =a[jp+p01];
			t02=t00-rt;		t03=t01-it;
			t00=t00+rt;		t01=t01+it;

			t04=a[jt+p02];	t05=a[jp+p02];
			rt =a[jt+p03];	it =a[jp+p03];
			t06=t04-rt;		t07=t05-it;
			t04=t04+rt;		t05=t05+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t04;		it =t05;
			t04=t00-rt;		t05=t01-it;
			t00=t00+rt;		t01=t01+it;

			rt =t06;		it =t07;
			t06=t02-it;		t07=t03+rt;
			t02=t02+it;		t03=t03-rt;

			t08=a[jt+p04];	t09=a[jp+p04];
			rt =a[jt+p05];	it =a[jp+p05];
			t0A=t08-rt;		t0B=t09-it;
			t08=t08+rt;		t09=t09+it;

			t0C=a[jt+p06];	t0D=a[jp+p06];
			rt =a[jt+p07];	it =a[jp+p07];
			t0E=t0C-rt;		t0F=t0D-it;
			t0C=t0C+rt;		t0D=t0D+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t0C;		it =t0D;
			t0C=t08-rt;		t0D=t09-it;
			t08=t08+rt;		t09=t09+it;

			rt =t0E;		it =t0F;
			t0E=t0A-it;		t0F=t0B+rt;
			t0A=t0A+it;		t0B=t0B-rt;

			rt =t08;		it =t09;
			t08=t00-rt;		t09=t01-it;
			t00=t00+rt;		t01=t01+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t0C;		it =t0D;
			t0C=t04-it;		t0D=t05+rt;
			t04=t04+it;		t05=t05-rt;

			rt =(t0A+t0B)*ISRT2;it =(t0A-t0B)*ISRT2;
			t0A=t02-rt;		t0B=t03+it;
			t02=t02+rt;		t03=t03-it;

			rt =(t0E-t0F)*ISRT2;it =(t0F+t0E)*ISRT2;
			t0E=t06+rt;		t0F=t07+it;
			t06=t06-rt;		t07=t07-it;

	/*...Block 2:;	*/

		jt = j1 + p08;
		jp = j2 + p08;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#elif PFETCH
		addr = &a[jt];
		addp = addr+prefetch_offset;
		prefetch_p_doubles(addp);
	#endif
			t10=a[jt    ];	t11=a[jp    ];
			rt =a[jt+p01];	it =a[jp+p01];
			t12=t10-rt;		t13=t11-it;
			t10=t10+rt;		t11=t11+it;

			t14=a[jt+p02];	t15=a[jp+p02];
			rt =a[jt+p03];	it =a[jp+p03];
			t16=t14-rt;		t17=t15-it;
			t14=t14+rt;		t15=t15+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t14;		it =t15;
			t14=t10-rt;		t15=t11-it;
			t10=t10+rt;		t11=t11+it;

			rt =t16;		it =t17;
			t16=t12-it;		t17=t13+rt;
			t12=t12+it;		t13=t13-rt;

			t18=a[jt+p04];	t19=a[jp+p04];
			rt =a[jt+p05];	it =a[jp+p05];
			t1A=t18-rt;		t1B=t19-it;
			t18=t18+rt;		t19=t19+it;

			t1C=a[jt+p06];	t1D=a[jp+p06];
			rt =a[jt+p07];	it =a[jp+p07];
			t1E=t1C-rt;		t1F=t1D-it;
			t1C=t1C+rt;		t1D=t1D+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t1C;		it =t1D;
			t1C=t18-rt;		t1D=t19-it;
			t18=t18+rt;		t19=t19+it;

			rt =t1E;		it =t1F;
			t1E=t1A-it;		t1F=t1B+rt;
			t1A=t1A+it;		t1B=t1B-rt;

			rt =t18;		it =t19;
			t18=t10-rt;		t19=t11-it;
			t10=t10+rt;		t11=t11+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t1C;		it =t1D;
			t1C=t14-it;		t1D=t15+rt;
			t14=t14+it;		t15=t15-rt;

			rt =(t1A+t1B)*ISRT2;it =(t1A-t1B)*ISRT2;
			t1A=t12-rt;		t1B=t13+it;
			t12=t12+rt;		t13=t13-it;

			rt =(t1E-t1F)*ISRT2;it =(t1F+t1E)*ISRT2;
			t1E=t16+rt;		t1F=t17+it;
			t16=t16-rt;		t17=t17-it;

	/*...Block 3:	*/

		jt = j1 + p10;
		jp = j2 + p10;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#elif PFETCH
		addr = &a[jt];
		addp = addr+prefetch_offset;
		prefetch_p_doubles(addp);
	#endif
			t20=a[jt    ];	t21=a[jp    ];
			rt =a[jt+p01];	it =a[jp+p01];
			t22=t20-rt;		t23=t21-it;
			t20=t20+rt;		t21=t21+it;

			t24=a[jt+p02];	t25=a[jp+p02];
			rt =a[jt+p03];	it =a[jp+p03];
			t26=t24-rt;		t27=t25-it;
			t24=t24+rt;		t25=t25+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t24;		it =t25;
			t24=t20-rt;		t25=t21-it;
			t20=t20+rt;		t21=t21+it;

			rt =t26;		it =t27;
			t26=t22-it;		t27=t23+rt;
			t22=t22+it;		t23=t23-rt;

			t28=a[jt+p04];	t29=a[jp+p04];
			rt =a[jt+p05];	it =a[jp+p05];
			t2A=t28-rt;		t2B=t29-it;
			t28=t28+rt;		t29=t29+it;

			t2C=a[jt+p06];	t2D=a[jp+p06];
			rt =a[jt+p07];	it =a[jp+p07];
			t2E=t2C-rt;		t2F=t2D-it;
			t2C=t2C+rt;		t2D=t2D+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t2C;		it =t2D;
			t2C=t28-rt;		t2D=t29-it;
			t28=t28+rt;		t29=t29+it;

			rt =t2E;		it =t2F;
			t2E=t2A-it;		t2F=t2B+rt;
			t2A=t2A+it;		t2B=t2B-rt;

			rt =t28;		it =t29;
			t28=t20-rt;		t29=t21-it;
			t20=t20+rt;		t21=t21+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t2C;		it =t2D;
			t2C=t24-it;		t2D=t25+rt;
			t24=t24+it;		t25=t25-rt;

			rt =(t2A+t2B)*ISRT2;it =(t2A-t2B)*ISRT2;
			t2A=t22-rt;		t2B=t23+it;
			t22=t22+rt;		t23=t23-it;

			rt =(t2E-t2F)*ISRT2;it =(t2F+t2E)*ISRT2;
			t2E=t26+rt;		t2F=t27+it;
			t26=t26-rt;		t27=t27-it;

	/*...Block 4:	*/

		jt = j1 + p18;
		jp = j2 + p18;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#elif PFETCH
		addr = &a[jt];
		addp = addr+prefetch_offset;
		prefetch_p_doubles(addp);
	#endif
			t30=a[jt    ];	t31=a[jp    ];
			rt =a[jt+p01];	it =a[jp+p01];
			t32=t30-rt;		t33=t31-it;
			t30=t30+rt;		t31=t31+it;

			t34=a[jt+p02];	t35=a[jp+p02];
			rt =a[jt+p03];	it =a[jp+p03];
			t36=t34-rt;		t37=t35-it;
			t34=t34+rt;		t35=t35+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t34;		it =t35;
			t34=t30-rt;		t35=t31-it;
			t30=t30+rt;		t31=t31+it;

			rt =t36;		it =t37;
			t36=t32-it;		t37=t33+rt;
			t32=t32+it;		t33=t33-rt;

			t38=a[jt+p04];	t39=a[jp+p04];
			rt =a[jt+p05];	it =a[jp+p05];
			t3A=t38-rt;		t3B=t39-it;
			t38=t38+rt;		t39=t39+it;

			t3C=a[jt+p06];	t3D=a[jp+p06];
			rt =a[jt+p07];	it =a[jp+p07];
			t3E=t3C-rt;		t3F=t3D-it;
			t3C=t3C+rt;		t3D=t3D+it;
	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif

			rt =t3C;		it =t3D;
			t3C=t38-rt;		t3D=t39-it;
			t38=t38+rt;		t39=t39+it;

			rt =t3E;		it =t3F;
			t3E=t3A-it;		t3F=t3B+rt;
			t3A=t3A+it;		t3B=t3B-rt;

			rt =t38;		it =t39;
			t38=t30-rt;		t39=t31-it;
			t30=t30+rt;		t31=t31+it;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t3C;		it =t3D;
			t3C=t34-it;		t3D=t35+rt;
			t34=t34+it;		t35=t35-rt;

			rt =(t3A+t3B)*ISRT2;it =(t3A-t3B)*ISRT2;
			t3A=t32-rt;		t3B=t33+it;
			t32=t32+rt;		t33=t33-it;

			rt =(t3E-t3F)*ISRT2;it =(t3F+t3E)*ISRT2;
			t3E=t36+rt;		t3F=t37+it;
			t36=t36-rt;		t37=t37-it;

	/*...and now do eight radix-4 transforms, including the internal twiddle factors:
		1, exp(-i* 1*twopi/32) =       ( c32_1,-s32_1), exp(-i* 2*twopi/32) =       ( c    ,-s    ), exp(-i* 3*twopi/32) =       ( c32_3,-s32_3) (for inputs to transform block 2),
		1, exp(-i* 2*twopi/32) =       ( c    ,-s    ), exp(-i* 4*twopi/32) = ISRT2*( 1    ,-1    ), exp(-i* 6*twopi/32) =       ( s    ,-c    ) (for inputs to transform block 3),
		1, exp(-i* 3*twopi/32) =       ( c32_3,-s32_3), exp(-i* 6*twopi/32) =       ( s    ,-c    ), exp(-i* 9*twopi/32) =       (-s32_1,-c32_1) (for inputs to transform block 4),
		1, exp(-i* 4*twopi/32) = ISRT2*( 1    ,-1    ), exp(-i* 8*twopi/32) =       ( 0    ,-1    ), exp(-i*12*twopi/32) = ISRT2*(-1    ,-1    ) (for inputs to transform block 5),
		1, exp(-i* 5*twopi/32) =       ( s32_3,-c32_3), exp(-i*10*twopi/32) =       (-s    ,-c    ), exp(-i*15*twopi/32) =       (-c32_1,-s32_1) (for inputs to transform block 6),
		1, exp(-i* 6*twopi/32) =       ( s    ,-c    ), exp(-i*12*twopi/32) = ISRT2*(-1    ,-1    ), exp(-i*18*twopi/32) =       (-c    , s    ) (for inputs to transform block 7),
		1, exp(-i* 7*twopi/32) =       ( s32_1,-c32_1), exp(-i*14*twopi/32) =       (-c    ,-s    ), exp(-i*21*twopi/32) =       (-s32_3, c32_3) (for inputs to transform block 8),
		 and only the last 3 inputs to each of the radix-4 transforms 2 through 8 are multiplied by non-unity twiddles.	*/

	/*...Block 1: t00,t10,t20,t30	*/

		jt = j1;
		jp = j2;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p10;
		prefetch_p_doubles(addp);
	#endif
			rt =t10;	t10=t00-rt;	t00=t00+rt;
			it =t11;	t11=t01-it;	t01=t01+it;

			rt =t30;	t30=t20-rt;	t20=t20+rt;
			it =t31;	t31=t21-it;	t21=t21+it;

			a[jt    ]=t00+t20;		a[jp    ]=t01+t21;
			t00      =t00-t20;		t01        =t01-t21;
			a[jt+p10]=t00*c10+t01*s10;	a[jp+p10]=t01*c10-t00*s10;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt       =t10+t31;		it         =t11-t30;	/* mpy by E^-4 = -I is inlined here...	*/
			t10      =t10-t31;		t11        =t11+t30;
			a[jt+p08]=rt *c08+it *s08;	a[jp+p08]=it *c08-rt *s08;
			a[jt+p18]=t10*c18+t11*s18;	a[jp+p18]=t11*c18-t10*s18;

	/*...Block 5: t08,t18,t28,t38	*/

		jt = j1 + p04;
		jp = j2 + p04;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t18;	t18=t08-t19;	t08=t08+t19;		/* twiddle mpy by E^8 =-I	*/
				t19=t09+rt;	t09=t09-rt;

			rt =(t29+t28)*ISRT2;	t29=(t29-t28)*ISRT2;		t28=rt;	/* twiddle mpy by E^-4	*/
			rt =(t38-t39)*ISRT2;	it =(t38+t39)*ISRT2;			/* twiddle mpy by E^4 = -E^-12 is here...	*/
			t38=t28+rt;			t28=t28-rt;				/* ...and get E^-12 by flipping signs here.	*/
			t39=t29+it;			t29=t29-it;

			rt       =t08+t28;		it         =t09+t29;
			t08      =t08-t28;		t09        =t09-t29;
			a[jt    ]=rt *c04+it *s04;	a[jp    ]=it *c04-rt *s04;
			a[jt+p10]=t08*c14+t09*s14;	a[jp+p10]=t09*c14-t08*s14;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt       =t18+t39;		it         =t19-t38;	/* mpy by E^-4 = -I is inlined here...	*/
			t18      =t18-t39;		t19        =t19+t38;
			a[jt+p08]=rt *c0C+it *s0C;	a[jp+p08]=it *c0C-rt *s0C;
			a[jt+p18]=t18*c1C+t19*s1C;	a[jp+p18]=t19*c1C-t18*s1C;

	/*...Block 3: t04,t14,t24,t34	*/

		jt = j1 + p02;
		jp = j2 + p02;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =(t15+t14)*ISRT2;	it =(t15-t14)*ISRT2;			/* twiddle mpy by E^-4	*/
			t14=t04-rt;			t04=t04+rt;
			t15=t05-it;			t05=t05+it;

			rt =t24*c + t25*s;		t25=t25*c - t24*s;		t24=rt;	/* twiddle mpy by E^-2	*/
			rt =t34*s + t35*c;		it =t35*s - t34*c;			/* twiddle mpy by E^-6	*/
			t34=t24-rt;			t24=t24+rt;
			t35=t25-it;			t25=t25+it;

			rt       =t04+t24;		it         =t05+t25;
			t04      =t04-t24;		t05        =t05-t25;
			a[jt    ]=rt *c02+it *s02;	a[jp    ]=it *c02-rt *s02;
			a[jt+p10]=t04*c12+t05*s12;	a[jp+p10]=t05*c12-t04*s12;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt       =t14+t35;		it         =t15-t34;	/* mpy by E^-4 = -I is inlined here...	*/
			t14      =t14-t35;		t15        =t15+t34;
			a[jt+p08]=rt *c0A+it *s0A;	a[jp+p08]=it *c0A-rt *s0A;
			a[jt+p18]=t14*c1A+t15*s1A;	a[jp+p18]=t15*c1A-t14*s1A;

	/*...Block 7: t0C,t1C,t2C,t3C	*/

		jt = j1 + p06;
		jp = j2 + p06;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =(t1C-t1D)*ISRT2;	it =(t1C+t1D)*ISRT2;			/* twiddle mpy by E^4 = -E^-12 is here...	*/
			t1C=t0C+rt;			t0C=t0C-rt;				/* ...and get E^-12 by flipping signs here.	*/
			t1D=t0D+it;			t0D=t0D-it;

			rt =t2C*s + t2D*c;		t2D=t2D*s - t2C*c;		t2C=rt;	/* twiddle mpy by E^-6	*/
			rt =t3C*c + t3D*s;		it =t3D*c - t3C*s;			/* twiddle mpy by E^-18 is here...	*/
			t3C=t2C+rt;			t2C=t2C-rt;				/* ...and get E^-18 by flipping signs here.	*/
			t3D=t2D+it;			t2D=t2D-it;

			rt       =t0C+t2C;		it         =t0D+t2D;
			t0C      =t0C-t2C;		t0D        =t0D-t2D;
			a[jt    ]=rt *c06+it *s06;	a[jp    ]=it *c06-rt *s06;
			a[jt+p10]=t0C*c16+t0D*s16;	a[jp+p10]=t0D*c16-t0C*s16;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt       =t1C+t3D;		it         =t1D-t3C;	/* mpy by E^-4 = -I is inlined here...	*/
			t1C      =t1C-t3D;		t1D        =t1D+t3C;
			a[jt+p08]=rt *c0E+it *s0E;	a[jp+p08]=it *c0E-rt *s0E;
			a[jt+p18]=t1C*c1E+t1D*s1E;	a[jp+p18]=t1D*c1E-t1C*s1E;

	/*...Block 2: t02,t12,t22,t32	*/

		jt = j1 + p01;
		jp = j2 + p01;

	#ifdef PFETCH_AGGRESSIVE
		addp = addr + p18;
		prefetch_p_doubles(addp);
	#endif
			rt =t12*c + t13*s;		it =t13*c - t12*s;			/* twiddle mpy by E^-2	*/
			t12=t02-rt;			t02=t02+rt;
			t13=t03-it;			t03=t03+it;

			rt =t22*c32_1 + t23*s32_1;	t23=t23*c32_1 - t22*s32_1;	t22=rt;	/* twiddle mpy by E^-1	*/
			rt =t32*c32_3 + t33*s32_3;	it =t33*c32_3 - t32*s32_3;		/* twiddle mpy by E^-3	*/
			t32=t22-rt;			t22=t22+rt;
			t33=t23-it;			t23=t23+it;

			rt       =t02+t22;		it         =t03+t23;
			t02      =t02-t22;		t03        =t03-t23;
			a[jt    ]=rt *c01+it *s01;	a[jp    ]=it *c01-rt *s01;
			a[jt+p10]=t02*c11+t03*s11;	a[jp+p10]=t03*c11-t02*s11;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt       =t12+t33;		it         =t13-t32;	/* mpy by E^-4 = -I is inlined here...	*/
			t12      =t12-t33;		t13        =t13+t32;
			a[jt+p08]=rt *c09+it *s09;	a[jp+p08]=it *c09-rt *s09;
			a[jt+p18]=t12*c19+t13*s19;	a[jp+p18]=t13*c19-t12*s19;

	/*...Block 6: t0A,t1A,t2A,t3A	*/

		jt = j1 + p05;
		jp = j2 + p05;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t1A*s - t1B*c;		it =t1B*s + t1A*c;			/* twiddle mpy by -E^-10 is here...	*/
			t1A=t0A+rt;			t0A =t0A-rt;				/* ...and get E^-10 by flipping signs here.	*/
			t1B=t0B+it;			t0B =t0B-it;

			rt =t2A*s32_3 + t2B*c32_3;	t2B=t2B*s32_3 - t2A*c32_3;	t2A=rt;	/* twiddle mpy by E^-5	*/
			rt =t3A*c32_1 - t3B*s32_1;	it =t3B*c32_1 + t3A*s32_1;		/* twiddle mpy by -E^-15 is here...	*/
			t3A=t2A+rt;			t2A=t2A-rt;				/* ...and get E^-15 by flipping signs here.	*/
			t3B=t2B+it;			t2B=t2B-it;

			rt       =t0A+t2A;		it         =t0B+t2B;
			t0A      =t0A-t2A;		t0B        =t0B-t2B;
			a[jt    ]=rt *c05+it *s05;	a[jp    ]=it *c05-rt *s05;
			a[jt+p10]=t0A*c15+t0B*s15;	a[jp+p10]=t0B*c15-t0A*s15;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt       =t1A+t3B;		it         =t1B-t3A;	/* mpy by E^-4 = -I is inlined here...	*/
			t1A      =t1A-t3B;		t1B        =t1B+t3A;
			a[jt+p08]=rt *c0D+it *s0D;	a[jp+p08]=it *c0D-rt *s0D;
			a[jt+p18]=t1A*c1D+t1B*s1D;	a[jp+p18]=t1B*c1D-t1A*s1D;

	/*...Block 4: t06,t16,t26,t36	*/

		jt = j1 + p03;
		jp = j2 + p03;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t16*s + t17*c;		it =t17*s - t16*c;			/* twiddle mpy by E^-6	*/
			t16=t06-rt;			t06 =t06+rt;
			t17=t07-it;			t07 =t07+it;

			rt =t26*c32_3 + t27*s32_3;	t27=t27*c32_3 - t26*s32_3;	t26=rt;	/* twiddle mpy by E^-3	*/
			rt =t36*s32_1 - t37*c32_1;	it =t37*s32_1 + t36*c32_1;		/* twiddle mpy by -E^-9 is here...	*/
			t36=t26+rt;			t26=t26-rt;				/* ...and get E^-9 by flipping signs here.	*/
			t37=t27+it;			t27=t27-it;

			rt       =t06+t26;		it         =t07+t27;
			t06      =t06-t26;		t07        =t07-t27;
			a[jt    ]=rt *c03+it *s03;	a[jp    ]=it *c03-rt *s03;
			a[jt+p10]=t06*c13+t07*s13;	a[jp+p10]=t07*c13-t06*s13;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt       =t16+t37;		it         =t17-t36;	/* mpy by E^-4 = -I is inlined here...	*/
			t16      =t16-t37;		t17        =t17+t36;
			a[jt+p08]=rt *c0B+it *s0B;	a[jp+p08]=it *c0B-rt *s0B;
			a[jt+p18]=t16*c1B+t17*s1B;	a[jp+p18]=t17*c1B-t16*s1B;

	/*...Block 8: t0E,t1E,t2E,t3E	*/

		jt = j1 + p07;
		jp = j2 + p07;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt =t1E*c - t1F*s;		it =t1F*c + t1E*s;			/* twiddle mpy by -E^-14 is here...	*/
			t1E=t0E+rt;			t0E =t0E-rt;				/* ...and get E^-14 by flipping signs here.	*/
			t1F=t0F+it;			t0F =t0F-it;

			rt =t2E*s32_1 + t2F*c32_1;	t2F=t2F*s32_1 - t2E*c32_1;	t2E=rt;	/* twiddle mpy by E^-7	*/
			rt =t3E*s32_3 + t3F*c32_3;	it =t3F*s32_3 - t3E*c32_3;		/* twiddle mpy by -E^-21 is here...	*/
			t3E=t2E+rt;			t2E=t2E-rt;				/* ...and get E^-21 by flipping signs here.	*/
			t3F=t2F+it;			t2F=t2F-it;

			rt       =t0E+t2E;		it         =t0F+t2F;
			t0E      =t0E-t2E;		t0F        =t0F-t2F;
			a[jt    ]=rt *c07+it *s07;	a[jp    ]=it *c07-rt *s07;
			a[jt+p10]=t0E*c17+t0F*s17;	a[jp+p10]=t0F*c17-t0E*s17;

	#ifdef PFETCH_AGGRESSIVE
		addp += p01;
		prefetch_p_doubles(addp);
	#endif
			rt       =t1E+t3F;		it         =t1F-t3E;	/* mpy by E^-4 = -I is inlined here...	*/
			t1E      =t1E-t3F;		t1F        =t1F+t3E;
			a[jt+p08]=rt *c0F+it *s0F;	a[jp+p08]=it *c0F-rt *s0F;
			a[jt+p18]=t1E*c1F+t1F*s1F;	a[jp+p18]=t1F*c1F-t1E*s1F;

#endif	/* USE_SSE2 */

	  }	/* endfor(j=jlo; j < jhi; j += 4) */

	  /*jlo=jlo+incr; jhi=jhi+incr; See my note about OpenMP above. */

	}	/* endfor(m=0; m < nloops; m++) */

}

