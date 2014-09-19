/*******************************************************************************
*                                                                              *
*   (C) 1997-2014 by Ernst W. Mayer.                                           *
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

// This main loop is same for un-and-multithreaded, so stick into a header file
// (can't use a macro because of the #if-enclosed stuff).

for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
{
	/* In SIMD mode, data are arranged in [re_0,...,re_n-1,im_0,...,im_n-1] groups, not the usual [re_0,im_0],...,[re_n-1,im_n-1] pairs.
	Thus we can still increment the j-index as if stepping through the residue array-of-doubles in strides of 2,
	but to point to the proper real datum, we need to index-map e.g. [0,1,2,3] ==> [0,2,1,3] in 2-way SIMD mode.
	(But only ever need to explicitly do this in debug mode).
	*/
	for(j = jstart; j < jhi; j += stride)
	{
		j1 =  j;
		j1 = j1 + ( (j1 >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		j2 = j1 + RE_IM_STRIDE;

/*...The radix-32 DIT pass is here:	*/

#ifdef USE_SSE2

  #if defined(COMPILER_TYPE_MSVC)

	/*...Block 1: */
	#if 1
		add0 = &a[j1];
		__asm	mov eax, add0
		__asm	mov ebx, p1
		__asm	mov ecx, p2
		__asm	mov edx, p3
		__asm	mov	edi, p4		/* edi will store copy of p4 throughout */
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	shl	ecx, 3
		__asm	shl	edx, 3
		__asm	shl	edi, 3
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
		SSE2_RADIX4_DIT_0TWIDDLE_B(r00)
	#else
		add0 = &a[j1];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r00)
	#endif

	/*...Block 2: */
	#if 1
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
		SSE2_RADIX4_DIT_0TWIDDLE_B(r08)
	#else
		add0 = &a[j1+p4];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r08)
	#endif

	/*...Block 3: */
	#if 1
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
		SSE2_RADIX4_DIT_0TWIDDLE_B(r10)
	#else
		add0 = &a[j1+p8];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r10)
	#endif

	/*...Block 4: */
	#if 1
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
		SSE2_RADIX4_DIT_0TWIDDLE_B(r18)
	#else
		add0 = &a[j1+p12];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;
		/* DIT radix-4 subconvolution, sans twiddles.	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		SSE2_RADIX4_DIT_0TWIDDLE(add0, add1, add2, add3, r18)
	#endif

	/****************************************************************************************************
	!...and now do four more radix-4 transforms, including the internal [no external]twiddle factors:   !
	****************************************************************************************************/

	/*...Block 1: Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
		__asm	mov	eax, r00

		__asm	movaps	xmm0,[eax      ]	/* t1  */
		__asm	movaps	xmm1,[eax+0x010]	/* t2  */
		__asm	movaps	xmm2,[eax+0x080]	/* t9  */
		__asm	movaps	xmm3,[eax+0x090]	/* t10 */

		__asm	subpd	xmm0,[eax+0x080]	/*~t9 =t1 -t9 */
		__asm	subpd	xmm1,[eax+0x090]	/*~t10=t2 -t10*/
		__asm	addpd	xmm2,[eax      ]	/*~t1 =t9 +t1 */
		__asm	addpd	xmm3,[eax+0x010]	/*~t2 =t10+t2 */

		__asm	movaps	xmm4,[eax+0x100]	/* t17 */
		__asm	movaps	xmm5,[eax+0x110]	/* t18 */
		__asm	movaps	xmm6,[eax+0x180]	/* t25 */
		__asm	movaps	xmm7,[eax+0x190]	/* t26 */

		__asm	subpd	xmm4,[eax+0x180]	/*~t25=t17-t25*/
		__asm	subpd	xmm5,[eax+0x190]	/*~t26=t18-t26*/
		__asm	addpd	xmm6,[eax+0x100]	/*~t17=t25+t17*/
		__asm	addpd	xmm7,[eax+0x110]	/*~t18=t26+t18*/

	/*
		t1       =t1+t17;				t2       =t2+t18;
		t17     *=     2;				t18     *=     2;
		t17      =t1-t17;				t18      =t2-t18;
	*/
		__asm	subpd	xmm2,xmm6		/* t1 -t17 */
		__asm	subpd	xmm3,xmm7		/* t2 -t18 */
		__asm	movaps	[eax+0x100],xmm2	/* a[jt+p8 ]/t17 */
		__asm	movaps	[eax+0x110],xmm3	/* a[jp+p8 ]/t18 */
		__asm	addpd	xmm6,xmm6		/*   2*t17 */
		__asm	addpd	xmm7,xmm7		/*   2*t18 */
		__asm	addpd	xmm6,xmm2		/*~t17 <- t1 +t17 */
		__asm	addpd	xmm7,xmm3		/*~t18 <- t2 +t18 */
		__asm	movaps	[eax      ],xmm6	/* a[jt+p0 ]/t1  */
		__asm	movaps	[eax+0x010],xmm7	/* a[jp+p0 ]/t2  */

	/*
		t9       =t9 +t26;				t10      =t10-t25;	// mpy by E^-4 = -I is inlined here...
		t26     *=     2;				t25     *=     2;
		t26      =t9 -t26;				t25      =t10+t25;
	*/
		__asm	subpd	xmm0,xmm5		/* t9 -t26 */
		__asm	subpd	xmm1,xmm4		/* t10-t25 */
		__asm	movaps	[eax+0x180],xmm0	/* a[jt+p12]/t26 */
		__asm	movaps	[eax+0x090],xmm1	/* a[jp+p4 ]/t10 */
		__asm	addpd	xmm5,xmm5		/*   2*t26 */
		__asm	addpd	xmm4,xmm4		/*   2*t25 */
		__asm	addpd	xmm5,xmm0		/* t9 +t26 */
		__asm	addpd	xmm4,xmm1		/* t10+t25 */
		__asm	movaps	[eax+0x080],xmm5	/* a[jt+p4 ]/t9  */
		__asm	movaps	[eax+0x190],xmm4	/* a[jp+p12]/t25 */

	/*...Block 2: Cost: 19 MOVapd, 24 ADD/SUBpd,  4 MULpd */
		__asm	mov	eax, r04
		__asm	mov	ebx, isrt2
		__asm	movaps	xmm2,[ebx]	/* isrt2 */

		__asm	movaps	xmm4,[eax+0x100]	/* t21 */
		__asm	movaps	xmm5,[eax+0x110]	/* t22 */
		__asm	movaps	xmm0,[eax+0x180]	/* t29 */
		__asm	movaps	xmm1,[eax+0x190]	/* t30 */

		__asm	addpd	xmm4,[eax+0x110]	/*~t21=t21+t22*/
		__asm	subpd	xmm5,[eax+0x100]	/*~t22=t22-t21*/
		__asm	subpd	xmm0,[eax+0x190]	/* rt =t29-t30*/
		__asm	addpd	xmm1,[eax+0x180]	/* it =t30+t29*/
		__asm	mulpd	xmm4,xmm2
		__asm	mulpd	xmm5,xmm2
		__asm	mulpd	xmm0,xmm2
		__asm	mulpd	xmm1,xmm2
		__asm	movaps	xmm6,xmm4			/* t21 copy */
		__asm	movaps	xmm7,xmm5			/* t22 copy */

		__asm	subpd	xmm4,xmm0			/*~t21=t21-rt */
		__asm	subpd	xmm5,xmm1			/*~t22=t22-it */
		__asm	addpd	xmm6,xmm0			/*~t29=t21+rt */
		__asm	addpd	xmm7,xmm1			/*~t30=t22+it */

		__asm	movaps	xmm0,[eax      ]	/* t5  */
		__asm	movaps	xmm1,[eax+0x010]	/* t6  */
		__asm	movaps	xmm2,[eax+0x080]	/* t13 */
		__asm	movaps	xmm3,[eax+0x090]	/* t14 */

		__asm	subpd	xmm0,[eax+0x090]	/*~t13=t5 -t14*/
		__asm	subpd	xmm1,[eax+0x080]	/*~t6 =t6 -t13*/
		__asm	addpd	xmm3,[eax      ]	/*~t5 =t14+t5 */
		__asm	addpd	xmm2,[eax+0x010]	/*~t14=t13+t6 */

	/*
		t5       =t5 +t21;			t6       =t6 +t22;
		t21     *=     2;			t22     *=     2;
		t21      =t5 -t21;			t22      =t6 -t22;
	*/
		__asm	subpd	xmm3,xmm4		/*~t21 <- t5 -t21 */
		__asm	subpd	xmm1,xmm5		/*~t22 <- t6 -t22 */
		__asm	movaps	[eax+0x100],xmm3	/* a[jt+p8 ]/t21 */
		__asm	movaps	[eax+0x110],xmm1	/* a[jp+p8 ]/t22 */
		__asm	addpd	xmm4,xmm4		/*          2*t21 */
		__asm	addpd	xmm5,xmm5		/*          2*t22 */
		__asm	addpd	xmm4,xmm3		/*~t5  <- t5 +t21 */
		__asm	addpd	xmm5,xmm1		/*~t6  <- t6 +t22 */
		__asm	movaps	[eax      ],xmm4	/* a[jt+p0 ]/t5  */
		__asm	movaps	[eax+0x010],xmm5	/* a[jp+p0 ]/t6  */
	/*
		t13      =t13+t30;			t14      =t14-t29;	// mpy by E^-4 = -I is inlined here...
		t30     *=     2;			t29     *=     2;
		t30      =t13-t30;			t29      =t14+t29;
	*/
		__asm	subpd	xmm0,xmm7		/*~t30 <- t13-t30 */
		__asm	subpd	xmm2,xmm6		/*~t14 <- t14-t29 */
		__asm	movaps	[eax+0x180],xmm0	/* a[jt+p12]/t30 */
		__asm	movaps	[eax+0x090],xmm2	/* a[jp+p4 ]/t14 */
		__asm	addpd	xmm7,xmm7		/*          2*t30 */
		__asm	addpd	xmm6,xmm6		/*          2*t29 */
		__asm	addpd	xmm7,xmm0		/*~t13 <- t13+t30 */
		__asm	addpd	xmm6,xmm2		/*~t29 <- t14+t29 */
		__asm	movaps	[eax+0x080],xmm7	/* a[jt+p4 ]/t13 */
		__asm	movaps	[eax+0x190],xmm6	/* a[jp+p12]/t29 */

	/*...Block 3: Cost: 22 MOVapd, 28 ADD/SUBpd, 10 MULpd */
		__asm	mov	eax, r02
		__asm	mov	ebx, isrt2
		__asm	mov	ecx, cc0

		__asm	movaps	xmm4,[eax+0x100]	/* t19 */				__asm	movaps	xmm0,[eax+0x180]	/* t27 */
		__asm	movaps	xmm5,[eax+0x110]	/* t20 */				__asm	movaps	xmm1,[eax+0x190]	/* t28 */
		__asm	movaps	xmm6,[eax+0x100]	/* xmm2 <- cpy t19 */	__asm	movaps	xmm2,[eax+0x180]	/* xmm6 <- cpy t27 */
		__asm	movaps	xmm7,[eax+0x110]	/* xmm3 <- cpy t20 */	__asm	movaps	xmm3,[eax+0x190]	/* xmm7 <- cpy t28 */

		__asm	mulpd	xmm4,[ecx     ]	/* t19*c */					__asm	mulpd	xmm0,[ecx+0x10]	/* t27*s */
		__asm	mulpd	xmm5,[ecx     ]	/* t20*c */					__asm	mulpd	xmm1,[ecx+0x10]	/* t28*s */
		__asm	mulpd	xmm6,[ecx+0x10]	/* t19*s */					__asm	mulpd	xmm2,[ecx     ]	/* t27*c */
		__asm	mulpd	xmm7,[ecx+0x10]	/* t20*s */					__asm	mulpd	xmm3,[ecx     ]	/* t28*c */
		__asm	subpd	xmm5,xmm6	/* xmm1 <-~t20*/				__asm	subpd	xmm1,xmm2	/* xmm5 <- it */
		__asm	addpd	xmm4,xmm7	/* xmm0 <-~t19*/				__asm	addpd	xmm0,xmm3	/* xmm4 <- rt */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy~t20*/
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy~t19*/

		__asm	addpd	xmm4,xmm0	/* ~t19 <- t19+rt */
		__asm	addpd	xmm5,xmm1	/* ~t20 <- t20+it */
		__asm	subpd	xmm6,xmm0	/* ~t27 <- t19-rt */
		__asm	subpd	xmm7,xmm1	/* ~t28 <- t20-it */

		__asm	movaps	xmm2,[eax+0x080]	/* t11 */
		__asm	movaps	xmm3,[eax+0x090]	/* t12 */
		__asm	movaps	xmm0,[eax      ]	/* t3  */
		__asm	movaps	xmm1,[eax+0x010]	/* t4  */
		__asm	addpd	xmm2,[eax+0x090]	/*~t11=t11+t12*/
		__asm	subpd	xmm3,[eax+0x080]	/*~t12=t12-t11*/
		__asm	mulpd	xmm2,[ebx]	/* rt */
		__asm	mulpd	xmm3,[ebx]	/* it */

		__asm	subpd	xmm0,xmm2	/*~t11 <- t3 - rt */
		__asm	subpd	xmm1,xmm3	/*~t12 <- t4 - it */
		__asm	addpd	xmm2,xmm2	/*          2* rt */
		__asm	addpd	xmm3,xmm3	/*          2* it */
		__asm	addpd	xmm2,xmm0	/*~t3  <- t3 + rt */
		__asm	addpd	xmm3,xmm1	/*~t4  <- t4 + it */

	/*
		t3       =t3 +t19;			t4       =t4 +t20;
		t19     *=     2;			t20     *=     2;
		t19      =t3 -t19;			t20      =t4 -t20;
	*/
		__asm	subpd	xmm2,xmm4		/*~t19 <- t3 -t19 */
		__asm	subpd	xmm3,xmm5		/*~t20 <- t4 -t20 */
		__asm	movaps	[eax+0x100],xmm2	/* a[jt+p8 ]/t19 */
		__asm	movaps	[eax+0x110],xmm3	/* a[jp+p8 ]/t20 */
		__asm	addpd	xmm4,xmm4		/*          2*t19 */
		__asm	addpd	xmm5,xmm5		/*          2*t20 */
		__asm	addpd	xmm4,xmm2		/* rt  <- t3 +t19 */
		__asm	addpd	xmm5,xmm3		/* it  <- t4 +t20 */
		__asm	movaps	[eax      ],xmm4	/* a[jt    ]/t3  */
		__asm	movaps	[eax+0x010],xmm5	/* a[jp    ]/t4  */

	/*
		t11      =t11+t28;			t12      =t12-t27;	// mpy by E^-4 = -I is inlined here...
		t28     *=     2;			t27     *=     2;
		t28      =t11-t28;			t27      =t12+t27;
	*/
		__asm	subpd	xmm0,xmm7		/*~t28 <- t11-t28 */
		__asm	subpd	xmm1,xmm6		/*~t12 <- t12-t27 */
		__asm	movaps	[eax+0x180],xmm0	/* a[jt+p12]/t28 */
		__asm	movaps	[eax+0x090],xmm1	/* a[jp+p4 ]/t12 */
		__asm	addpd	xmm7,xmm7		/*          2*t28 */
		__asm	addpd	xmm6,xmm6		/*          2*t27 */
		__asm	addpd	xmm7,xmm0		/*~t11 <- t11+t28 */
		__asm	addpd	xmm6,xmm1		/*~t27 <- t12+t27 */
		__asm	movaps	[eax+0x080],xmm7	/* a[jt+p4 ]/t11 */
		__asm	movaps	[eax+0x190],xmm6	/* a[jp+p12]/t27 */

	/*...Block 4: Cost: 22 MOVapd, 28 ADD/SUBpd, 10 MULpd */
		__asm	mov	eax, r06
		__asm	mov	ebx, isrt2
		__asm	mov	ecx, cc0

		__asm	movaps	xmm4,[eax+0x100]	/* t23 */				__asm	movaps	xmm0,[eax+0x180]	/* t31 */
		__asm	movaps	xmm5,[eax+0x110]	/* t24 */				__asm	movaps	xmm1,[eax+0x190]	/* t32 */
		__asm	movaps	xmm6,[eax+0x100]	/* xmm2 <- cpy t23 */	__asm	movaps	xmm2,[eax+0x180]	/* xmm6 <- cpy t31 */
		__asm	movaps	xmm7,[eax+0x110]	/* xmm3 <- cpy t24 */	__asm	movaps	xmm3,[eax+0x190]	/* xmm7 <- cpy t32 */

		__asm	mulpd	xmm4,[ecx+0x10]	/* t23*s */					__asm	mulpd	xmm0,[ecx     ]	/* t31*c */
		__asm	mulpd	xmm5,[ecx+0x10]	/* t24*s */					__asm	mulpd	xmm1,[ecx     ]	/* t32*c */
		__asm	mulpd	xmm6,[ecx     ]	/* t23*c */					__asm	mulpd	xmm2,[ecx+0x10]	/* t31*s */
		__asm	mulpd	xmm7,[ecx     ]	/* t24*c */					__asm	mulpd	xmm3,[ecx+0x10]	/* t32*s */
		__asm	subpd	xmm5,xmm6	/* xmm1 <-~t24*/				__asm	subpd	xmm1,xmm2	/* xmm5 <- it */
		__asm	addpd	xmm4,xmm7	/* xmm0 <-~t23*/				__asm	addpd	xmm0,xmm3	/* xmm4 <- rt */
		__asm	movaps	xmm7,xmm5	/* xmm3 <- cpy~t24*/
		__asm	movaps	xmm6,xmm4	/* xmm2 <- cpy~t23*/

		__asm	addpd	xmm4,xmm0	/* ~t31 <- t23+rt */
		__asm	addpd	xmm5,xmm1	/* ~t32 <- t24+it */
		__asm	subpd	xmm6,xmm0	/* ~t23 <- t23-rt */
		__asm	subpd	xmm7,xmm1	/* ~t24 <- t24-it */

		__asm	movaps	xmm2,[eax+0x080]	/* t15 */
		__asm	movaps	xmm3,[eax+0x090]	/* t16 */
		__asm	movaps	xmm0,[eax      ]	/* t7  */
		__asm	movaps	xmm1,[eax+0x010]	/* t8  */
		__asm	subpd	xmm2,[eax+0x090]	/*~t15=t15-t16*/
		__asm	addpd	xmm3,[eax+0x080]	/*~t16=t16+t15*/
		__asm	mulpd	xmm2,[ebx]	/* rt */
		__asm	mulpd	xmm3,[ebx]	/* it */

		__asm	subpd	xmm0,xmm2	/*~t7  <- t7 - rt */
		__asm	subpd	xmm1,xmm3	/*~t8  <- t8 - it */
		__asm	addpd	xmm2,xmm2	/*          2* rt */
		__asm	addpd	xmm3,xmm3	/*          2* it */
		__asm	addpd	xmm2,xmm0	/*~t15 <- t7 + rt */
		__asm	addpd	xmm3,xmm1	/*~t16 <- t8 + it */

	/*
		t7       =t7 +t23;			t8       =t8 +t24;
		t23     *=     2;			t24     *=     2;
		t23      =t7 -t23;			t24      =t8 -t24;
	*/
		__asm	subpd	xmm0,xmm6		/*~t23 <- t7 -t23 */
		__asm	subpd	xmm1,xmm7		/*~t24 <- t8 -t24 */
		__asm	movaps	[eax+0x100],xmm0	/* a[jt+p8 ]/t23 */
		__asm	movaps	[eax+0x110],xmm1	/* a[jp+p8 ]/t24 */
		__asm	addpd	xmm6,xmm6		/*          2*t23 */
		__asm	addpd	xmm7,xmm7		/*          2*t24 */
		__asm	addpd	xmm6,xmm0		/* rt  <- t7 +t23 */
		__asm	addpd	xmm7,xmm1		/* it  <- t8 +t24 */
		__asm	movaps	[eax      ],xmm6	/* a[jt+p0 ]/t7  */
		__asm	movaps	[eax+0x010],xmm7	/* a[jp+p0 ]/t8  */

	/*
		t15      =t15+t32;			t16      =t16-t31;	// mpy by E^-4 = -I is inlined here...
		t32     *=     2;			t31     *=     2;
		t32      =t15-t32;			t31      =t16+t31;
	*/
		__asm	subpd	xmm2,xmm5		/*~t32 <- t15-t32 */
		__asm	subpd	xmm3,xmm4		/*~t16 <- t16-t31 */
		__asm	movaps	[eax+0x180],xmm2	/* a[jt+p12]/t19 */
		__asm	movaps	[eax+0x090],xmm3	/* a[jp+p4 ]/t20 */
		__asm	addpd	xmm5,xmm5		/*          2*t32 */
		__asm	addpd	xmm4,xmm4		/*          2*t31 */
		__asm	addpd	xmm5,xmm2		/*~t15 <- t15+t32 */
		__asm	addpd	xmm4,xmm3		/*~t31 <- t16+t31 */
		__asm	movaps	[eax+0x080],xmm5	/* a[jt+p4 ]/t19 */
		__asm	movaps	[eax+0x190],xmm4	/* a[jp+p12]/t20 */

		/***************************************************/
		/* DIT Totals: 143 MOVapd, 180 ADD/SUBpd, 24 MULpd */
		/***************************************************/

  #else	/* GCC-style inline ASM: */

		add0 = &a[j1];
		SSE2_RADIX16_DIT_NOTWIDDLE(add0,p1,p2,p3,p4,r00,r02,r04,r06,r08,r0a,r10,r18,isrt2,cc0);

  #endif

#else	/* !USE_SSE2 */

	#if USE_SCALAR_DFT_MACRO
		RADIX_16_DIT(a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
					,a1p0r,a1p0i,a1p1r,a1p1i,a1p2r,a1p2i,a1p3r,a1p3i,a1p4r,a1p4i,a1p5r,a1p5i,a1p6r,a1p6i,a1p7r,a1p7i,a1p8r,a1p8i,a1p9r,a1p9i,a1pAr,a1pAi,a1pBr,a1pBi,a1pCr,a1pCi,a1pDr,a1pDi,a1pEr,a1pEi,a1pFr,a1pFi
					,c,s)
	#else
		/*...Block 1:	*/
		t1 =a[j1    ];	t2 =a[j2    ];
		rt =a[j1+p1 ];	it =a[j2+p1 ];
		t3 =t1 -rt;		t1 =t1 +rt;
		t4 =t2 -it;		t2 =t2 +it;

		t5 =a[j1+p2 ];	t6 =a[j2+p2 ];
		rt =a[j1+p3 ];	it =a[j2+p3 ];
		t7 =t5 -rt;  	t5 =t5 +rt;
		t8 =t6 -it;  	t6 =t6 +it;

		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;

		rt =t7;	t7 =t3 -t8;	t3 =t3 +t8;
		t8 =t4 +rt;	t4 =t4 -rt;

		/*...Block 2:	*/
		t9 =a[j1+p4 ];	t10=a[j2+p4 ];
		rt =a[j1+p5 ];	it =a[j2+p5 ];
		t11=t9 -rt;		t9 =t9 +rt;
		t12=t10-it;		t10=t10+it;

		t13=a[j1+p6 ];	t14=a[j2+p6 ];
		rt =a[j1+p7 ];	it =a[j2+p7 ];
		t15=t13-rt;  	t13=t13+rt;
		t16=t14-it;		t14=t14+it;

		rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
		it =t14;	t14=t10-it;	t10=t10+it;

		rt =t15;	t15=t11-t16;	t11=t11+t16;
		t16=t12+rt;	t12=t12-rt;

		/*...Block 3:	*/
		t17=a[j1+p8 ];	t18=a[j2+p8 ];
		rt =a[j1+p9 ];	it =a[j2+p9 ];
		t19=t17-rt;		t17=t17+rt;
		t20=t18-it;		t18=t18+it;

		t21=a[j1+p10];	t22=a[j2+p10];
		rt =a[j1+p11];	it =a[j2+p11];
		t23=t21-rt;  	t21=t21+rt;
		t24=t22-it;		t22=t22+it;

		rt =t21;	t21=t17-rt;	t17=t17+rt;
		it =t22;	t22=t18-it;	t18=t18+it;

		rt =t23;	t23=t19-t24;	t19=t19+t24;
		t24=t20+rt;	t20=t20-rt;

		/*...Block 4:	*/
		t25=a[j1+p12];	t26=a[j2+p12];
		rt =a[j1+p13];	it =a[j2+p13];
		t27=t25-rt;		t25=t25+rt;
		t28=t26-it;		t26=t26+it;

		t29=a[j1+p14];	t30=a[j2+p14];
		rt =a[j1+p15];	it =a[j2+p15];
		t31=t29-rt;  	t29=t29+rt;
		t32=t30-it;		t30=t30+it;

		rt =t29;	t29=t25-rt;	t25=t25+rt;
		it =t30;	t30=t26-it;	t26=t26+it;

		rt =t31;	t31=t27-t32;	t27=t27+t32;
		t32=t28+rt;	t28=t28-rt;

		/*...and now do four more radix-4 transforms, including the internal twiddle factors:
		1, exp(-i* 1*twopi/16) =       ( c,-s), exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 3*twopi/16) =       ( s,-c) (for inputs to transform block 2)
		1, exp(-i* 2*twopi/16) = ISRT2*( 1,-1), exp(-i* 4*twopi/16) =       ( 0,-1), exp(-i* 6*twopi/16) = ISRT2*(-1,-1) (for inputs to transform block 3)
		1, exp(-i* 3*twopi/16) =       ( s,-c), exp(-i* 6*twopi/16) = ISRT2*(-1,-1), exp(-i* 9*twopi/16) =       (-c, s) (for inputs to transform block 4).
		(This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
		I.e. do similar as above, except inputs a(j1  +p0:15:1) are replaced by t0:30:2,
									 a(j2+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
		and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.	*/

		/*...Block 1: t1,9,17,25	*/
		rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;

		a1p0r =t1+t17;	a1p0i =t2+t18;
		a1p8r =t1-t17;	a1p8i =t2-t18;

		a1p4r =t9 +t26;	a1p4i =t10-t25;	/* mpy by E^-4 = -I is inlined here...	*/
		a1pCr=t9 -t26;	a1pCi=t10+t25;

		/*...Block 3: t5,13,21,29	*/
		rt =t13;	t13=t5 -t14;	t5 =t5 +t14;		/* twiddle mpy by E^4 =-I	*/
		t14=t6 +rt;	t6 =t6 -rt;

		rt =(t22+t21)*ISRT2;t22=(t22-t21)*ISRT2;	t21=rt;	/* twiddle mpy by E^-2	*/
		rt =(t29-t30)*ISRT2;it =(t29+t30)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here...	*/
		t29=t21+rt;		t21=t21-rt;			/* ...and get E^-6 by flipping signs here.	*/
		t30=t22+it;		t22=t22-it;

		a1p2r =t5+t21;	a1p2i =t6+t22;
		a1pAr=t5-t21;	a1pAi=t6-t22;

		a1p6r =t13+t30;	a1p6i =t14-t29;	/* mpy by E^-4 = -I is inlined here...	*/
		a1pEr=t13-t30;	a1pEi=t14+t29;

		/*...Block 2: t3,11,19,27	*/
		rt =(t12+t11)*ISRT2;it =(t12-t11)*ISRT2;		/* twiddle mpy by E^-2	*/
		t11=t3 -rt;		t3 =t3 +rt;
		t12=t4 -it;		t4 =t4 +it;

		rt =t19*c + t20*s;	t20=t20*c - t19*s;	t19=rt;	/* twiddle mpy by E^-1	*/
		rt =t27*s + t28*c;	it =t28*s - t27*c;		/* twiddle mpy by E^-3	*/
		t27=t19-rt;		t19=t19+rt;
		t28=t20-it;		t20=t20+it;

		a1p1r =t3+t19;	a1p1i =t4+t20;
		a1p9r =t3-t19;	a1p9i =t4-t20;

		a1p5r =t11+t28;	a1p5i =t12-t27;	/* mpy by E^-4 = -I is inlined here...	*/
		a1pDr=t11-t28;	a1pDi=t12+t27;

		/*...Block 4: t7,15,23,31	*/
		rt =(t15-t16)*ISRT2;it =(t15+t16)*ISRT2;		/* twiddle mpy by E^2 = -E^-6 is here...	*/
		t15=t7 +rt;		t7 =t7 -rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here.	*/
		t16=t8 +it;		t8 =t8 -it;

		rt =t23*s + t24*c;	t24=t24*s - t23*c;	t23=rt;	/* twiddle mpy by E^-3	*/
		rt =t31*c + t32*s;	it =t32*c - t31*s;		/* twiddle mpy by E^-1 = -E^-9...	*/
		t31=t23+rt;		t23=t23-rt;			/* ...and get E^9 by flipping signs here.	*/
		t32=t24+it;		t24=t24-it;

		a1p3r =t7+t23;	a1p3i =t8+t24;
		a1pBr=t7-t23;	a1pBi=t8-t24;

		a1p7r =t15+t32;	a1p7i =t16-t31;	/* mpy by E^-4 = -I is inlined here...	*/
		a1pFr=t15-t32;	a1pFi=t16+t31;
	#endif	// USE_SCALAR_DFT_MACRO ?

#endif	/* USE_SSE2 */

	/*...and combine those to complete the radix-16 transform and do the carries. Since the outputs would
	normally be getting dispatched to 16 separate blocks of the A-array, we need 16 separate carries.	*/

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
	#ifdef USE_AVX

		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		l= j & (nwt-1);						tmp = half_arr + 64;	/* ptr to local storage for the doubled wtl,wtn terms: */
		n_minus_sil  ->d0 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d0 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d0 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d0 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+2) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
		n_minus_sil  ->d1 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d1 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d1 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d1 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+4) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
		n_minus_sil  ->d2 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d2 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d2 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d2 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+6) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
		n_minus_sil  ->d3 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d3 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d3 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d3 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		AVX_cmplx_carry_norm_pow2_errcheck0_X4(r00,add1,add2,add3,cy_r0,bjmodn0,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
		AVX_cmplx_carry_norm_pow2_errcheck1_X4(r08,add1,add2,add3,cy_r4,bjmodn4,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
		AVX_cmplx_carry_norm_pow2_errcheck1_X4(r10,add1,add2,add3,cy_r8,bjmodn8,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
		AVX_cmplx_carry_norm_pow2_errcheck1_X4(r18,add1,add2,add3,cy_rC,bjmodnC,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		i =((uint32)(sw - *bjmodn0) >> 31);	/* get ready for the next set...	*/

	#elif defined(USE_SSE2)

		/* In SSE2 mode, the data are arranged in memory like so, where we view things in 16-byte chunks:

			R0:	a0.re,b0.re		I0:	a0.im,b0.im
			R1:	a1.re,b1.re		I1:	a1.im,b1.im
			R2:	a2.re,b2.re		I2:	a2.im,b2.im
			R3:	a3.re,b3.re		I3:	a3.im,b3.im
			R4:	a4.re,b4.re		I4:	a4.im,b4.im
			R5:	a5.re,b5.re		I5:	a5.im,b5.im
			R6:	a6.re,b6.re		I6:	a6.im,b6.im
			R7:	a7.re,b7.re		I7:	a7.im,b7.im
			R8:	a8.re,b8.re		I8:	a8.im,b8.im
			R9:	a9.re,b9.re		I9:	a9.im,b9.im
			Ra:	aA.re,bA.re		Ia:	aA.im,bA.im
			Rb:	aB.re,bB.re		Ib:	aB.im,bB.im
			Rc:	aC.re,bC.re		Ic:	aC.im,bC.im
			Rd:	aD.re,bD.re		Id:	aD.im,bD.im
			Re:	aE.re,bE.re		Ie:	aE.im,bE.im
			Rf:	aF.re,bF.re		If:	aF.im,bF.im

		Where the R's and I's map to the local temps as follows: R0:f ==> r00:1E:2, I0:f ==> r01:1F:2 , and the
		a's and b's of each pair represent doubles which in non-SSE2 mode would be getting processed in the same relative
		position on subsequent passes through the main-array for-loop, i.e. among which the carry propagation proceeds as

			a0.re -> a0.im -> b0.re -> b0.im .

		Because of the undesirable intra-xmm-register data dependency this leads to, we instead require data arranged as

			R0/r00:	a0.re,a1.re		I0/r01:	a0.im,a1.im
			R1/r02:	b0.re,b1.re		I1/r03:	b0.im,b1.im

		We need to interleave these pairwise so as to swap the high word of each even-indexed R-and-I-pair
		with the low word of the subsequent odd-indexed pair, e.g. for R0/r00 and R1/r02:

				low		high	low		high
			R0	[a0.re,b0.re]	[a1.re,b1.re]	R1
				   |      \       /      |
				   |        \   /        |
				   |          x          |
				   |        /   \        |
				   V      /       \      V
			R0~	[a0.re,a1.re]	[b0.re,b1.re]	R1~, and analogously for I0/r01 and I1/r03.

		This is the same butterfly swap pattern as is used in the wrapper_square routines. The other nice things about this:

			1) Even though e.g. a0 and a1 appear adjacent, they are actually n/16 memory locations apart, i.e. there
			   is no carry propagation between them;

			2) Processing a[j] and a[j+1] together means we access the following elements of the wt1[] array paiwise in the carry step:

				xmm.lo:			xmm.hi:
				wt1[col+j]		wt1[col+(j+1)]
				wt1[co2-j]		wt1[co2-(j+1)]
				wt1[co3-j]		wt1[co3-(j+1)]

			Thus these wt-array elements are also adjacent in memory and can be loaded pairwise into an XMM register
			[With an unaligned movupd load and a shufpd-based lo/hi-word swap needed on the latter two.]
		*/
		/* These indices remain constant throughout the carry block below - only change when loop index j does: */

		l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
		n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
		n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
		sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

	  #if defined(COMPILER_TYPE_MSVC)

	   #ifdef ERR_CHECK_ALL
		SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r00,add1,add2,add3,cy_r0,cy_r2,bjmodn0);
		SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r08,add1,add2,add3,cy_r4,cy_r6,bjmodn4);
		SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r10,add1,add2,add3,cy_r8,cy_rA,bjmodn8);
		SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r18,add1,add2,add3,cy_rC,cy_rE,bjmodnC);
	   #else
		SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r00,add1,add2,add3,cy_r0,cy_r2,bjmodn0);
		SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r08,add1,add2,add3,cy_r4,cy_r6,bjmodn4);
		SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r10,add1,add2,add3,cy_r8,cy_rA,bjmodn8);
		SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r18,add1,add2,add3,cy_rC,cy_rE,bjmodnC);
	   #endif

	  #else	/* GCC-style inline ASM: */

	   #ifdef ERR_CHECK_ALL
		SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r00,add1,add2,add3,cy_r0,cy_r2,bjmodn0,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
		SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r08,add1,add2,add3,cy_r4,cy_r6,bjmodn4,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
		SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r10,add1,add2,add3,cy_r8,cy_rA,bjmodn8,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
		SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r18,add1,add2,add3,cy_rC,cy_rE,bjmodnC,half_arr,  n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
	   #else
		SSE2_cmplx_carry_norm_pow2_errcheck0_2B(r00,add1,add2,add3,cy_r0,cy_r2,bjmodn0,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
		SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r08,add1,add2,add3,cy_r4,cy_r6,bjmodn4,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
		SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r10,add1,add2,add3,cy_r8,cy_rA,bjmodn8,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
		SSE2_cmplx_carry_norm_pow2_nocheck1_2B (r18,add1,add2,add3,cy_rC,cy_rE,bjmodnC,half_arr,  n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
	   #endif

	  #endif

		l= (j+2) & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
		n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
		n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
		sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		ctmp = (struct complex *)half_arr + 16;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/
		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

	  #if defined(COMPILER_TYPE_MSVC)

	   #ifdef ERR_CHECK_ALL
		SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r00,add1,add2,cy_r0,cy_r2,bjmodn0);
		SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r08,add1,add2,cy_r4,cy_r6,bjmodn4);
		SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r10,add1,add2,cy_r8,cy_rA,bjmodn8);
		SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r18,add1,add2,cy_rC,cy_rE,bjmodnC);
	   #else
		SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r00,add1,add2,cy_r0,cy_r2,bjmodn0);
		SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r08,add1,add2,cy_r4,cy_r6,bjmodn4);
		SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r10,add1,add2,cy_r8,cy_rA,bjmodn8);
		SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r18,add1,add2,cy_rC,cy_rE,bjmodnC);
	   #endif

	  #else	/* GCC-style inline ASM: */

	   #ifdef ERR_CHECK_ALL
		SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r00,add1,add2,cy_r0,cy_r2,bjmodn0,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
		SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r08,add1,add2,cy_r4,cy_r6,bjmodn4,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
		SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r10,add1,add2,cy_r8,cy_rA,bjmodn8,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
		SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r18,add1,add2,cy_rC,cy_rE,bjmodnC,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
	   #else
		SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r00,add1,add2,cy_r0,cy_r2,bjmodn0,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
		SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r08,add1,add2,cy_r4,cy_r6,bjmodn4,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
		SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r10,add1,add2,cy_r8,cy_rA,bjmodn8,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
		SSE2_cmplx_carry_norm_pow2_nocheck2_2B (r18,add1,add2,cy_rC,cy_rE,bjmodnC,half_arr,n_minus_silp1,n_minus_sil,          sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw);
	   #endif

	  #endif

		i =((uint32)(sw - *bjmodn0) >> 31);	/* get ready for the next set...	*/

	#else	// Scalar-double mode:

		l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
		n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
		n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
		sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

	/*...set0 is slightly different from others:	*/
	   cmplx_carry_norm_pow2_errcheck0(a1p0r,a1p0i,cy_r0,bjmodn0    );
		cmplx_carry_norm_pow2_errcheck(a1p1r,a1p1i,cy_r1,bjmodn1,0x1);
		cmplx_carry_norm_pow2_errcheck(a1p2r,a1p2i,cy_r2,bjmodn2,0x2);
		cmplx_carry_norm_pow2_errcheck(a1p3r,a1p3i,cy_r3,bjmodn3,0x3);
		cmplx_carry_norm_pow2_errcheck(a1p4r,a1p4i,cy_r4,bjmodn4,0x4);
		cmplx_carry_norm_pow2_errcheck(a1p5r,a1p5i,cy_r5,bjmodn5,0x5);
		cmplx_carry_norm_pow2_errcheck(a1p6r,a1p6i,cy_r6,bjmodn6,0x6);
		cmplx_carry_norm_pow2_errcheck(a1p7r,a1p7i,cy_r7,bjmodn7,0x7);
		cmplx_carry_norm_pow2_errcheck(a1p8r,a1p8i,cy_r8,bjmodn8,0x8);
		cmplx_carry_norm_pow2_errcheck(a1p9r,a1p9i,cy_r9,bjmodn9,0x9);
		cmplx_carry_norm_pow2_errcheck(a1pAr,a1pAi,cy_rA,bjmodnA,0xA);
		cmplx_carry_norm_pow2_errcheck(a1pBr,a1pBi,cy_rB,bjmodnB,0xB);
		cmplx_carry_norm_pow2_errcheck(a1pCr,a1pCi,cy_rC,bjmodnC,0xC);
		cmplx_carry_norm_pow2_errcheck(a1pDr,a1pDi,cy_rD,bjmodnD,0xD);
		cmplx_carry_norm_pow2_errcheck(a1pEr,a1pEi,cy_rE,bjmodnE,0xE);
		cmplx_carry_norm_pow2_errcheck(a1pFr,a1pFi,cy_rF,bjmodnF,0xF);

		i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?

	}
	else	/* MODULUS_TYPE_FERMAT */
	{

	#ifdef USE_AVX

		// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.
		tmp = half_arr+2;
		VEC_DBL_INIT(tmp, scale);
		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		idx_offset = j;
		idx_incr = NDIVR;

		tmp = base_negacyclic_root;	tm2 = tmp+1;

		// Hi-accuracy version needs 4 copies of each base root:
		l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
		VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
		VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
		VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
		VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
		VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
		VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
		VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
		VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
		VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
		VEC_DBL_INIT(tmp+  8,wt_re);	VEC_DBL_INIT(tm2+  8,wt_im);
		VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
		VEC_DBL_INIT(tmp+ 24,wt_re);	VEC_DBL_INIT(tm2+ 24,wt_im);

		// AVX-custom 4-way carry macro - each contains 4 of the RADIX stride-n/RADIX-separated carries
		// (processed independently in parallel), and steps through sequential-data indicies j,j+2,j+4,j+6:

		// The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*2*sizeof(vec_dbl)
		// which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
		// precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
		// uses the same 0x40-byte up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls:
		l = 0x400;
		tmp = base_negacyclic_root+  0;	SSE2_fermat_carry_norm_pow2_errcheck_X4(r00,tmp,l,cy_r0,cy_i0,half_arr,sign_mask); l -= 0xc0;
		tmp = base_negacyclic_root+  8;	SSE2_fermat_carry_norm_pow2_errcheck_X4(r08,tmp,l,cy_r4,cy_i4,half_arr,sign_mask); l -= 0xc0;
		tmp = base_negacyclic_root+ 16;	SSE2_fermat_carry_norm_pow2_errcheck_X4(r10,tmp,l,cy_r8,cy_i8,half_arr,sign_mask); l -= 0xc0;
		tmp = base_negacyclic_root+ 24;	SSE2_fermat_carry_norm_pow2_errcheck_X4(r18,tmp,l,cy_rC,cy_iC,half_arr,sign_mask);

	#elif defined(USE_SSE2)

		/* In SSE2 mode, carry propagation proceeds as

			a0.re -> b0.re;		a0.im -> b0.im, where these imaginary parts really represent elements
												a0.im = a[n/2] and b0.im = a[n/2+1] of the right-angle transform.

		This data layout is ideal for the negacyclic unweighting/reweighting step bracketing the carry step, but in the latter,
		because of the undesirable intra-xmm-register data dependency this leads to, we instead require data arranged as

			R0/r00:	a0.re,a0.im		I0/r01:	b0.re,b0.im, i.e. the non-SSE2 data layout works best in the carry step!

		We need to interleave these pairwise so as to swap the high word of each R-element
		with the low word of the corresponding I-element, e.g. for R0/r00 and I0/r01:

				low		high	low		high
			R0	[a0.re,b0.re]	[a0.im,b0.im]	I0
				   |      \       /      |
				   |        \   /        |
				   |          x          |
				   |        /   \        |
				   V      /       \      V
			R0~	[a0.re,a0.im]	[b0.re,b0.im]	I0~.

		Note that even though e.g. a0 and a1 appear adjacent in terms of their a-subscripts, they are actually
		n/16 memory locations apart, i.e. there is no carry propagation between them.
		*/

		tmp = half_arr+2;
		VEC_DBL_INIT(tmp, scale);
		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		idx_offset = j;
		idx_incr = NDIVR;

	  #if defined(COMPILER_TYPE_MSVC)
		/* The cy_[r|i]_idx[A|B] names here are not meaningful, each simple stores one [re,im] carry pair,
		e.g. cy_r0 stores the carries our of [a0.re,a0.im], cy_r2 stores the carries our of [a1.re,a1.im], etc.
		Here is the actual mapping between these SSE2-mode 2-vector carry pairs and the scalar carries:
											  2-vector                               Scalar
											 ----------                            ----------- */
		SSE2_fermat_carry_norm_pow2_errcheck(r00,cy_r0,idx_offset,idx_incr);	/* cy_r0,cy_i0 */
		SSE2_fermat_carry_norm_pow2_errcheck(r02,cy_r2,idx_offset,idx_incr);	/* cy_r1,cy_i1 */
		SSE2_fermat_carry_norm_pow2_errcheck(r04,cy_r4,idx_offset,idx_incr);	/* cy_r2,cy_i2 */
		SSE2_fermat_carry_norm_pow2_errcheck(r06,cy_r6,idx_offset,idx_incr);	/* cy_r3,cy_i3 */
		SSE2_fermat_carry_norm_pow2_errcheck(r08,cy_r8,idx_offset,idx_incr);	/* cy_r4,cy_i4 */
		SSE2_fermat_carry_norm_pow2_errcheck(r0A,cy_rA,idx_offset,idx_incr);	/* cy_r5,cy_i5 */
		SSE2_fermat_carry_norm_pow2_errcheck(r0C,cy_rC,idx_offset,idx_incr);	/* cy_r6,cy_i6 */
		SSE2_fermat_carry_norm_pow2_errcheck(r0E,cy_rE,idx_offset,idx_incr);	/* cy_r7,cy_i7 */
		SSE2_fermat_carry_norm_pow2_errcheck(r10,cy_i0,idx_offset,idx_incr);	/* cy_r8,cy_i8 */
		SSE2_fermat_carry_norm_pow2_errcheck(r12,cy_i2,idx_offset,idx_incr);	/* cy_r9,cy_i9 */
		SSE2_fermat_carry_norm_pow2_errcheck(r14,cy_i4,idx_offset,idx_incr);	/* cy_rA,cy_iA */
		SSE2_fermat_carry_norm_pow2_errcheck(r16,cy_i6,idx_offset,idx_incr);	/* cy_rB,cy_iB */
		SSE2_fermat_carry_norm_pow2_errcheck(r18,cy_i8,idx_offset,idx_incr);	/* cy_rC,cy_iC */
		SSE2_fermat_carry_norm_pow2_errcheck(r1A,cy_iA,idx_offset,idx_incr);	/* cy_rD,cy_iD */
		SSE2_fermat_carry_norm_pow2_errcheck(r1C,cy_iC,idx_offset,idx_incr);	/* cy_rE,cy_iE */
		SSE2_fermat_carry_norm_pow2_errcheck(r1E,cy_iE,idx_offset,idx_incr);	/* cy_rF,cy_iF */

	  #elif (OS_BITS == 32)

		SSE2_fermat_carry_norm_pow2_errcheck(r00,cy_r0,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck(r02,cy_r2,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck(r04,cy_r4,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck(r06,cy_r6,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck(r08,cy_r8,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck(r0A,cy_rA,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck(r0C,cy_rC,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck(r0E,cy_rE,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck(r10,cy_i0,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck(r12,cy_i2,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck(r14,cy_i4,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck(r16,cy_i6,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck(r18,cy_i8,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck(r1A,cy_iA,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck(r1C,cy_iC,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck(r1E,cy_iE,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);

	  #else	// 64-bit SSE2

		SSE2_fermat_carry_norm_pow2_errcheck_X2(r00,cy_r0,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck_X2(r04,cy_r4,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck_X2(r08,cy_r8,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck_X2(r0C,cy_rC,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck_X2(r10,cy_i0,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck_X2(r14,cy_i4,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck_X2(r18,cy_i8,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);
		SSE2_fermat_carry_norm_pow2_errcheck_X2(r1C,cy_iC,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2);

	  #endif

	#else	// Scalar-double mode:

		ntmp = 0;
		fermat_carry_norm_pow2_errcheck(a1p0r,a1p0i,cy_r0,cy_i0,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a1p1r,a1p1i,cy_r1,cy_i1,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a1p2r,a1p2i,cy_r2,cy_i2,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a1p3r,a1p3i,cy_r3,cy_i3,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a1p4r,a1p4i,cy_r4,cy_i4,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a1p5r,a1p5i,cy_r5,cy_i5,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a1p6r,a1p6i,cy_r6,cy_i6,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a1p7r,a1p7i,cy_r7,cy_i7,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a1p8r,a1p8i,cy_r8,cy_i8,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a1p9r,a1p9i,cy_r9,cy_i9,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a1pAr,a1pAi,cy_rA,cy_iA,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a1pBr,a1pBi,cy_rB,cy_iB,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a1pCr,a1pCi,cy_rC,cy_iC,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a1pDr,a1pDi,cy_rD,cy_iD,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a1pEr,a1pEi,cy_rE,cy_iE,ntmp,NRTM1,NRT_BITS);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a1pFr,a1pFi,cy_rF,cy_iF,ntmp,NRTM1,NRT_BITS);

	#endif	/* #ifdef USE_SSE2 */

	}	/* if(MODULUS_TYPE == ...) */

/*...The radix-16 DIF pass is here:	*/

/* Four DIF radix-4 subconvolution, sans twiddles.	Cost each: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */

#ifdef USE_SSE2

  #if defined(COMPILER_TYPE_MSVC)

		SSE2_RADIX4_DIF_IN_PLACE(r00, r10, r08, r18)
		SSE2_RADIX4_DIF_IN_PLACE(r04, r14, r0C, r1C)
		SSE2_RADIX4_DIF_IN_PLACE(r02, r12, r0A, r1A)
		SSE2_RADIX4_DIF_IN_PLACE(r06, r16, r0E, r1E)

	/****************************************************************************************
	!...and now do four more radix-4 transforms, including the internal twiddle factors.	!
	!																						!
	!	This is identical to latter half of radix16 DIF, except for the r-vector indexing,	!
	!	which permutes as follows:															!
	!																						!
	!			t1	t3	t5	t7	t9	t11	t13	t15	t17	t19	t21	t23	t25	t27	t29	t31				!
	!		==>	r00	r08	r10	r18	r04	r0C	r14	r1C	r02	r0A	r12	r1A	r06	r0E	r16	r1E				!
	!																						!
	****************************************************************************************/

	/*...Block 1: t1,9,17,25 ==> r00,04,02,06	Cost: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */
	#if 1	// if(1) - test out pure-asm version
	//	add0 = &a[j1];
		__asm	mov eax, add0
		__asm	mov ebx, p1
		__asm	mov ecx, p2
		__asm	mov edx, p3
		__asm	mov	edi, p4		/* edi will store copy of p4 throughout */
		__asm	shl	ebx, 3		/* Pointer offset for floating doubles */
		__asm	shl	ecx, 3
		__asm	shl	edx, 3
		__asm	shl	edi, 3
		__asm	add ebx, eax
		__asm	add ecx, eax
		__asm	add edx, eax
	#else
		add0 = &a[j1];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;

		__asm	mov	eax, add0	/* restore main-array indices */
		__asm	mov	ebx, add1
		__asm	mov	ecx, add2
		__asm	mov	edx, add3
	#endif
		__asm	mov	esi, r00

		__asm	movaps	xmm0,[esi      ]	/* t1  */
		__asm	movaps	xmm1,[esi+0x010]	/* t2  */
		__asm	movaps	xmm2,[esi+0x040]	/* t9  */
		__asm	movaps	xmm3,[esi+0x050]	/* t10 */

		__asm	subpd	xmm0,[esi+0x040]	/* t9 =t1 -rt */
		__asm	subpd	xmm1,[esi+0x050]	/* t10=t2 -it */
		__asm	addpd	xmm2,[esi      ]	/* t1 =t1 +rt */
		__asm	addpd	xmm3,[esi+0x010]	/* t2 =t2 +it */

		__asm	movaps	xmm4,[esi+0x020]	/* t17 */
		__asm	movaps	xmm5,[esi+0x030]	/* t18 */
		__asm	movaps	xmm6,[esi+0x060]	/* t25 */
		__asm	movaps	xmm7,[esi+0x070]	/* t26 */

		__asm	subpd	xmm4,[esi+0x060]	/* t25=t17-rt */
		__asm	subpd	xmm5,[esi+0x070]	/* t26=t18-it */
		__asm	addpd	xmm6,[esi+0x020]	/* t17=t17+rt */
		__asm	addpd	xmm7,[esi+0x030]	/* t18=t18+it */

		__asm	subpd	xmm2,xmm6		/* t1  <- t1 -t17 */
		__asm	subpd	xmm3,xmm7		/* t2  <- t2 -t18 */
		__asm	addpd	xmm6,xmm6		/*          2*t17 */
		__asm	addpd	xmm7,xmm7		/*          2*t18 */
		__asm	movaps	[ebx     ],xmm2	/* a[jt+p1 ] */
		__asm	movaps	[ebx+0x10],xmm3	/* a[jp+p1 ] */
		__asm	addpd	xmm6,xmm2		/* t17 <- t1 +t17 */
		__asm	addpd	xmm7,xmm3		/* t18 <- t2 +t18 */
		__asm	movaps	[eax     ],xmm6	/* a[jt+p0 ] */
		__asm	movaps	[eax+0x10],xmm7	/* a[jp+p0 ] */

		__asm	subpd	xmm0,xmm5		/* t9  <- t9 -t26 */
		__asm	subpd	xmm1,xmm4		/* t10 <- t10-t25 */
		__asm	addpd	xmm5,xmm5		/*          2*t26 */
		__asm	addpd	xmm4,xmm4		/*          2*t25 */
		__asm	movaps	[ecx     ],xmm0	/* a[jt+p2 ] */
		__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
		__asm	addpd	xmm5,xmm0		/* t26 <- t9 +t26 */
		__asm	addpd	xmm4,xmm1		/* t25 <- t10+t25 */
		__asm	movaps	[edx     ],xmm5	/* a[jt+p3 ] */
		__asm	movaps	[ecx+0x10],xmm4	/* a[jp+p2 ] */

	/*...Block 3: t5,13,21,29 ==> r10,14,12,16	Cost: 16 MOVapd, 26 ADD/SUBpd,  4 MULpd */
	#if 1
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
	#else
		add0 = &a[j1+p4];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;

		__asm	mov	eax, add0	/* restore main-array indices */
		__asm	mov	ebx, add1
		__asm	mov	ecx, add2
		__asm	mov	edx, add3
	#endif
		__asm	mov	esi, r10

		__asm	movaps	xmm0,[esi      ]	/* t5  */
		__asm	movaps	xmm1,[esi+0x010]	/* t6  */
		__asm	movaps	xmm2,[esi+0x040]	/* t13 */
		__asm	movaps	xmm3,[esi+0x050]	/* t14 */

		__asm	subpd	xmm0,[esi+0x050]	/* t5 =t5 -t14*/
		__asm	subpd	xmm1,[esi+0x040]	/* t14=t6 -t13*/
		__asm	addpd	xmm2,[esi+0x010]	/* t6 =t13+t6 */
		__asm	addpd	xmm3,[esi      ]	/* t13=t14+t5 */

		__asm	movaps	xmm4,[esi+0x020]	/* t21 */
		__asm	movaps	xmm5,[esi+0x030]	/* t22 */
		__asm	movaps	xmm6,[esi+0x060]	/* t29 */
		__asm	movaps	xmm7,[esi+0x070]	/* t30 */

		__asm	subpd	xmm4,[esi+0x030]	/* t21-t22 */
		__asm	addpd	xmm5,[esi+0x020]	/* t22+t21 */
		__asm	addpd	xmm6,[esi+0x070]	/* t29+t30 */
		__asm	subpd	xmm7,[esi+0x060]	/* t30-t29 */

		__asm	mov	esi, isrt2
		__asm	mulpd	xmm4,[esi]	/* t21 = (t21-t22)*ISRT2 */
		__asm	mulpd	xmm5,[esi]	/* t22 = (t22+t21)*ISRT2 */
		__asm	mulpd	xmm6,[esi]	/*  rt = (t29+t30)*ISRT2 */
		__asm	mulpd	xmm7,[esi]	/*  it = (t30-t29)*ISRT2 */

		__asm	subpd	xmm4,xmm6		/* t21=t21-rt */
		__asm	subpd	xmm5,xmm7		/* t22=t22-it */
		__asm	addpd	xmm6,xmm6		/*      2* rt */
		__asm	addpd	xmm7,xmm7		/*      2* it */
		__asm	addpd	xmm6,xmm4		/* t29=t21+rt */
		__asm	addpd	xmm7,xmm5		/* t30=t22+it */

		__asm	subpd	xmm0,xmm4		/* t5 -t21 */
		__asm	subpd	xmm2,xmm5		/* t6 -t22 */
		__asm	addpd	xmm4,xmm4		/*   2*t21 */
		__asm	addpd	xmm5,xmm5		/*   2*t22 */

		__asm	movaps	[ebx     ],xmm0	/* a[jt+p1 ] */
		__asm	movaps	[ebx+0x10],xmm2	/* a[jp+p1 ] */
		__asm	addpd	xmm4,xmm0		/* t5 +t21 */
		__asm	addpd	xmm5,xmm2		/* t6 +t22 */
		__asm	movaps	[eax     ],xmm4	/* a[jt+p0 ] */
		__asm	movaps	[eax+0x10],xmm5	/* a[jp+p0 ] */

		__asm	subpd	xmm3,xmm7		/* t13-t30 */
		__asm	subpd	xmm1,xmm6		/* t14-t29 */
		__asm	addpd	xmm7,xmm7		/*   2*t30 */
		__asm	addpd	xmm6,xmm6		/*   2*t29 */
		__asm	movaps	[ecx     ],xmm3	/* a[jt+p2 ] */
		__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
		__asm	addpd	xmm7,xmm3		/* t13+t30 */
		__asm	addpd	xmm6,xmm1		/* t14+t29 */
		__asm	movaps	[edx     ],xmm7	/* a[jt+p3 ] */
		__asm	movaps	[ecx+0x10],xmm6	/* a[jp+p2 ] */

	/*...Block 2: t3,11,19,27 ==> r08,0C,0A,0E	Cost: 18 MOVapd, 28 ADD/SUBpd, 10 MULpd */
	#if 1
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
	#else
		add0 = &a[j1+p8];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;

		__asm	mov	eax, add0	/* restore main-array indices */
		__asm	mov	ebx, add1
		__asm	mov	ecx, add2
		__asm	mov	edx, add3
	#endif
		__asm	mov	esi, r08
		__asm	mov	edi, cc0
		__asm	movaps	xmm4,[esi+0x020]	/* t19 */		__asm	movaps	xmm6,[esi+0x060]	/* t27 */
		__asm	movaps	xmm5,[esi+0x030]	/* t20 */		__asm	movaps	xmm7,[esi+0x070]	/* t28 */
		__asm	movaps	xmm0,xmm4			/* copy t19 */	__asm	movaps	xmm2,xmm6			/* copy t27 */
		__asm	movaps	xmm1,xmm5			/* copy t20 */	__asm	movaps	xmm3,xmm7			/* copy t28 */

		__asm	mulpd	xmm4,[edi     ]	/* t19*c */			__asm	mulpd	xmm6,[edi+0x10]	/* t27*s */
		__asm	mulpd	xmm1,[edi+0x10]	/* t20*s */			__asm	mulpd	xmm3,[edi     ]	/* t28*c */
		__asm	mulpd	xmm5,[edi     ]	/* t20*c */			__asm	mulpd	xmm7,[edi+0x10]	/* t28*s */
		__asm	mulpd	xmm0,[edi+0x10]	/* t19*s */			__asm	mulpd	xmm2,[edi     ]	/* t27*c */
		__asm	subpd	xmm4,xmm1	/* ~t19 */				__asm	subpd	xmm6,xmm3	/* rt */
		__asm	addpd	xmm5,xmm0	/* ~t20 */				__asm	addpd	xmm7,xmm2	/* it */

		__asm	subpd	xmm4,xmm6		/*~t27=t19-rt */
		__asm	subpd	xmm5,xmm7		/*~t28=t20-it */
		__asm	addpd	xmm6,xmm6		/*      2* rt */
		__asm	addpd	xmm7,xmm7		/*      2* it */
		__asm	addpd	xmm6,xmm4		/*~t19=t19+rt */
		__asm	addpd	xmm7,xmm5		/*~t20=t20+it */

		__asm	mov	edi, isrt2
		__asm	movaps	xmm2,[esi+0x040]	/* t11 */
		__asm	movaps	xmm3,[esi+0x050]	/* t12 */
		__asm	subpd	xmm2,[esi+0x050]	/* t11-t12 */
		__asm	addpd	xmm3,[esi+0x040]	/* t12+t11 */
		__asm	mulpd	xmm2,[edi]	/* rt = (t11-t12)*ISRT2 */
		__asm	mulpd	xmm3,[edi]	/* it = (t12+t11)*ISRT2 */

		__asm	movaps	xmm0,[esi      ]	/* t3  */
		__asm	movaps	xmm1,[esi+0x010]	/* t4  */

		__asm	mov	edi, p4		/* restore p4-based value of edi */
		__asm	shl	edi, 3

		__asm	subpd	xmm0,xmm2			/*~t11=t3 -rt */
		__asm	subpd	xmm1,xmm3			/*~t12=t4 -it */
		__asm	addpd	xmm2,[esi      ]	/*~t3 =rt +t3 */
		__asm	addpd	xmm3,[esi+0x010]	/*~t4 =it +t4 */

		__asm	subpd	xmm2,xmm6		/* t3 -t19 */
		__asm	subpd	xmm3,xmm7		/* t4 -t20 */
		__asm	addpd	xmm6,xmm6		/*   2*t19 */
		__asm	addpd	xmm7,xmm7		/*   2*t20 */
		__asm	movaps	[ebx     ],xmm2	/* a[jt+p1 ] */
		__asm	movaps	[ebx+0x10],xmm3	/* a[jp+p1 ] */
		__asm	addpd	xmm6,xmm2		/* t3 +t19 */
		__asm	addpd	xmm7,xmm3		/* t4 +t20 */
		__asm	movaps	[eax     ],xmm6	/* a[jt+p0 ] */
		__asm	movaps	[eax+0x10],xmm7	/* a[jp+p0 ] */

		__asm	subpd	xmm0,xmm5		/* t11-t28 */
		__asm	subpd	xmm1,xmm4		/* t12-t27 */
		__asm	addpd	xmm5,xmm5		/*          2*t28 */
		__asm	addpd	xmm4,xmm4		/*          2*t27 */
		__asm	movaps	[ecx     ],xmm0	/* a[jt+p2 ] */
		__asm	movaps	[edx+0x10],xmm1	/* a[jp+p3 ] */
		__asm	addpd	xmm5,xmm0		/* t11+t28 */
		__asm	addpd	xmm4,xmm1		/* t12+t27 */
		__asm	movaps	[edx     ],xmm5	/* a[jt+p3 ] */
		__asm	movaps	[ecx+0x10],xmm4	/* a[jp+p2 ] */

	/*...Block 4: t7,15,23,31 ==> r18,1C,1A,1E	Cost: 18 MOVapd, 28 ADD/SUBpd, 10 MULpd */
	#if 1
		__asm	add eax, edi
		__asm	add ebx, edi
		__asm	add ecx, edi
		__asm	add edx, edi
	#else
		add0 = &a[j1+p12];
		add1 = add0+p1;
		add2 = add0+p2;
		add3 = add0+p3;

		__asm	mov	eax, add0	/* restore main-array indices */
		__asm	mov	ebx, add1
		__asm	mov	ecx, add2
		__asm	mov	edx, add3
	#endif

		__asm	mov	esi, r18
		__asm	mov	edi, cc0
		__asm	movaps	xmm4,[esi+0x020]	/* t23 */		__asm	movaps	xmm6,[esi+0x060]	/* t31 */
		__asm	movaps	xmm5,[esi+0x030]	/* t24 */		__asm	movaps	xmm7,[esi+0x070]	/* t32 */
		__asm	movaps	xmm0,xmm4			/* copy t23 */	__asm	movaps	xmm2,xmm6			/* copy t31 */
		__asm	movaps	xmm1,xmm5			/* copy t24 */	__asm	movaps	xmm3,xmm7			/* copy t32 */

		__asm	mulpd	xmm4,[edi+0x10]	/* t23*s */			__asm	mulpd	xmm6,[edi     ]	/* t31*c */
		__asm	mulpd	xmm1,[edi     ]	/* t24*c */			__asm	mulpd	xmm3,[edi+0x10]	/* t32*s */
		__asm	mulpd	xmm5,[edi+0x10]	/* t24*s */			__asm	mulpd	xmm7,[edi     ]	/* t32*c */
		__asm	mulpd	xmm0,[edi     ]	/* t23*c */			__asm	mulpd	xmm2,[edi+0x10]	/* t31*s */
		__asm	subpd	xmm4,xmm1	/* ~t23 */				__asm	subpd	xmm6,xmm3	/* rt */
		__asm	addpd	xmm5,xmm0	/* ~t24 */				__asm	addpd	xmm7,xmm2	/* it */

		__asm	subpd	xmm4,xmm6		/*~t23=t23-rt */
		__asm	subpd	xmm5,xmm7		/*~t24=t24-it */
		__asm	addpd	xmm6,xmm6		/*      2* rt */
		__asm	addpd	xmm7,xmm7		/*      2* it */
		__asm	addpd	xmm6,xmm4		/*~t31=t23+rt */
		__asm	addpd	xmm7,xmm5		/*~t32=t24+it */

		__asm	mov	edi, isrt2
		__asm	movaps	xmm2,[esi+0x040]	/* t15 */
		__asm	movaps	xmm3,[esi+0x050]	/* t16 */
		__asm	addpd	xmm2,[esi+0x050]	/* t15+t16 */
		__asm	subpd	xmm3,[esi+0x040]	/* t16-t15 */
		__asm	mulpd	xmm2,[edi]	/* rt = (t15+t16)*ISRT2 */
		__asm	mulpd	xmm3,[edi]	/* it = (t16-t15)*ISRT2 */

		__asm	movaps	xmm0,[esi      ]	/* t7  */
		__asm	movaps	xmm1,[esi+0x010]	/* t8  */

		__asm	subpd	xmm0,xmm2			/*~t7 =t7 -rt */
		__asm	subpd	xmm1,xmm3			/*~t8 =t8 -it */
		__asm	addpd	xmm2,[esi      ]	/*~t15=rt +t7 */
		__asm	addpd	xmm3,[esi+0x010]	/*~t16=it +t8 */

		__asm	subpd	xmm0,xmm4		/* t7 -t23 */
		__asm	subpd	xmm1,xmm5		/* t8 -t24 */
		__asm	addpd	xmm4,xmm4		/*   2*t23 */
		__asm	addpd	xmm5,xmm5		/*   2*t24 */
		__asm	movaps	[ebx     ],xmm0	/* a[jt+p1 ] */
		__asm	movaps	[ebx+0x10],xmm1	/* a[jp+p1 ] */
		__asm	addpd	xmm4,xmm0		/* t7 +t23 */
		__asm	addpd	xmm5,xmm1		/* t8 +t24 */
		__asm	movaps	[eax     ],xmm4	/* a[jt+p0 ] */
		__asm	movaps	[eax+0x10],xmm5	/* a[jp+p0 ] */

		__asm	subpd	xmm2,xmm7		/* t15-t32 */
		__asm	subpd	xmm3,xmm6		/* t16-t31 */
		__asm	addpd	xmm7,xmm7		/*   2*t32 */
		__asm	addpd	xmm6,xmm6		/*   2*t31 */
		__asm	movaps	[ecx     ],xmm2	/* a[jt+p2 ] */
		__asm	movaps	[edx+0x10],xmm3	/* a[jp+p3 ] */
		__asm	addpd	xmm7,xmm2		/* t15+t32 */
		__asm	addpd	xmm6,xmm3		/* t16+t31 */
		__asm	movaps	[edx     ],xmm7	/* a[jt+p3 ] */
		__asm	movaps	[ecx+0x10],xmm6	/* a[jp+p2 ] */

		/***************************************************/
		/* DIF Totals: 132 MOVapd, 182 ADD/SUBpd, 24 MULpd */
		/***************************************************/

  #else	/* GCC-style inline ASM: */

	SSE2_RADIX16_DIF_NOTWIDDLE(add0,p1,p2,p3,p4,r00,r02,r04,r06,r08,r0A,r0C,r0E,r10,r12,r14,r16,r18,r1A,r1C,r1E,isrt2,cc0);

  #endif

#else	/* !USE_SSE2 */

	#if USE_SCALAR_DFT_MACRO
		RADIX_16_DIF(a1p0r,a1p0i,a1p1r,a1p1i,a1p2r,a1p2i,a1p3r,a1p3i,a1p4r,a1p4i,a1p5r,a1p5i,a1p6r,a1p6i,a1p7r,a1p7i,a1p8r,a1p8i,a1p9r,a1p9i,a1pAr,a1pAi,a1pBr,a1pBi,a1pCr,a1pCi,a1pDr,a1pDi,a1pEr,a1pEi,a1pFr,a1pFi
					,a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
					,c,s)
	#else
		/*       gather the needed data (16 64-bit complex, i.e. 32 64-bit reals) and do the first set of four length-4 transforms.	*/
	  #if PFETCH
		addr = &a[j1];
		prefetch_p_doubles(addr);
	  #endif
		/*...Block 1:	*/
		t3 =a1p0r -a1p8r;  	t1 =a1p0r +a1p8r;
		t4 =a1p0i -a1p8i;	t2 =a1p0i +a1p8i;

		t7 =a1p4r -a1pCr;	t5 =a1p4r +a1pCr;
		t8 =a1p4i -a1pCi;	t6 =a1p4i +a1pCi;
	  #if PFETCH
		addp = addr+p1;
		prefetch_p_doubles(addp);
	  #endif
		rt =t5;	t5 =t1 -rt;	t1 =t1 +rt;
		it =t6;	t6 =t2 -it;	t2 =t2 +it;

		rt =t7;	t7 =t3 +t8;	t3 =t3 -t8;
				t8 =t4 -rt;	t4 =t4 +rt;
	  #if PFETCH
		addp = addr+p2;
		prefetch_p_doubles(addp);
	  #endif
		/*...Block 2:	*/
		t11=a1p2r -a1pAr;	t9 =a1p2r +a1pAr;
		t12=a1p2i -a1pAi;	t10=a1p2i +a1pAi;

		t15=a1p6r -a1pEr;	t13=a1p6r +a1pEr;
		t16=a1p6i -a1pEi;	t14=a1p6i +a1pEi;
	  #if PFETCH
		addp = addr+p3;
		prefetch_p_doubles(addp);
	  #endif
		rt =t13;	t13=t9 -rt;	t9 =t9 +rt;
		it =t14;	t14=t10-it;	t10=t10+it;

		rt =t15;	t15=t11+t16;t11=t11-t16;
					t16=t12-rt;	t12=t12+rt;
	  #if PFETCH
		addp = addr+p4;
		prefetch_p_doubles(addp);
	  #endif
		/*...Block 3:	*/
		t19=a1p1r -a1p9r;  	t17=a1p1r +a1p9r;
		t20=a1p1i -a1p9i;	t18=a1p1i +a1p9i;

		t23=a1p5r -a1pDr;	t21=a1p5r +a1pDr;
		t24=a1p5i -a1pDi;	t22=a1p5i +a1pDi;
	  #if PFETCH
		addp = addr+p5;
		prefetch_p_doubles(addp);
	  #endif
		rt =t21;	t21=t17-rt;	t17=t17+rt;
		it =t22;	t22=t18-it;	t18=t18+it;

		rt =t23;	t23=t19+t24;	t19=t19-t24;
					t24=t20-rt;	t20=t20+rt;
	  #if PFETCH
		addp = addr+p6;
		prefetch_p_doubles(addp);
	  #endif
		/*...Block 4:	*/
		t27=a1p3r -a1pBr;	t25=a1p3r +a1pBr;
		t28=a1p3i -a1pBi;	t26=a1p3i +a1pBi;

		t31=a1p7r -a1pFr;	t29=a1p7r +a1pFr;
		t32=a1p7i -a1pFi;	t30=a1p7i +a1pFi;
	  #if PFETCH
		addp = addr+p7;
		prefetch_p_doubles(addp);
	  #endif
		rt =t29;	t29=t25-rt;	t25=t25+rt;
		it =t30;	t30=t26-it;	t26=t26+it;

		rt =t31;	t31=t27+t32;t27=t27-t32;
					t32=t28-rt;	t28=t28+rt;
	  #if PFETCH
		addp = addr+p8;
		prefetch_p_doubles(addp);
	  #endif

		/*...and now do four more radix-4 transforms, including the internal twiddle factors:
		1, exp(i* 1*twopi/16) =       ( c, s), exp(i* 2*twopi/16) = ISRT2*( 1, 1), exp(i* 3*twopi/16) =       ( s, c) (for inputs to transform block 2)
		1, exp(i* 2*twopi/16) = ISRT2*( 1, 1), exp(i* 4*twopi/16) =       ( 0, 1), exp(i* 6*twopi/16) = ISRT2*(-1, 1) (for inputs to transform block 3)
		1, exp(i* 3*twopi/16) =       ( s, c), exp(i* 6*twopi/16) = ISRT2*(-1, 1), exp(i* 9*twopi/16) =       (-c,-s) (for inputs to transform block 4).
		(This is 4 real*complex and 4 complex*complex multiplies (= 24 FMUL), compared to 6 and 4, respectively (= 28 FMUL) for my old scheme.)
		I.e. do similar as above, except inputs a(j1  +p0:15:1) are replaced by t0:30:2,
									 a(j2+p0:15:1) are replaced by t1:31:2, and v.v. for outputs,
		and only the last 3 inputs to each of the radix-4 transforms 2 through 4 are multiplied by non-unity twiddles.	*/

		/*...Block 1: t1,9,17,25	*/
		rt =t9 ;	t9 =t1 -rt;	t1 =t1 +rt;
		it =t10;	t10=t2 -it;	t2 =t2 +it;

		rt =t25;	t25=t17-rt;	t17=t17+rt;
		it =t26;	t26=t18-it;	t18=t18+it;
	  #if PFETCH
		addp = addr+p9;
		prefetch_p_doubles(addp);
	  #endif
		a[j1    ]=t1+t17;	a[j2    ]=t2+t18;
		a[j1+p1 ]=t1-t17;	a[j2+p1 ]=t2-t18;

		a[j1+p2 ]=t9 -t26;	a[j2+p2 ]=t10+t25;	/* mpy by E^4=i is inlined here...	*/
		a[j1+p3 ]=t9 +t26;	a[j2+p3 ]=t10-t25;
	  #if PFETCH
		addp = addr+p10;
		prefetch_p_doubles(addp);
	  #endif
		/*...Block 3: t5,13,21,29	*/
		rt =t13;	t13=t5 +t14;	t5 =t5 -t14;		/* twiddle mpy by E^4 = I	*/
		t14=t6 -rt;	t6 =t6 +rt;

		rt =(t21-t22)*ISRT2;t22=(t21+t22)*ISRT2;	t21=rt;	/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t21,22,29,30	*/
		rt =(t30+t29)*ISRT2;it =(t30-t29)*ISRT2;		/* twiddle mpy by -E^6 is here...	*/
		t29=t21+rt;		t21=t21-rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here.	*/
		t30=t22+it;		t22=t22-it;
	  #if PFETCH
		addp = addr+p11;
		prefetch_p_doubles(addp);
	  #endif
		a[j1+p4 ]=t5+t21;	a[j2+p4 ]=t6+t22;
		a[j1+p5 ]=t5-t21;	a[j2+p5 ]=t6-t22;

		a[j1+p6 ]=t13-t30;	a[j2+p6 ]=t14+t29;	/* mpy by E^4=i is inlined here...	*/
		a[j1+p7 ]=t13+t30;	a[j2+p7 ]=t14-t29;
	  #if PFETCH
		addp = addr+p12;
		prefetch_p_doubles(addp);
	  #endif
		/*...Block 2: t3,11,19,27	*/
		rt =(t11-t12)*ISRT2;it =(t11+t12)*ISRT2;		/* twiddle mpy by E^2; could also do mpy by ISRT2 directly on t11,12	*/
		t11=t3 -rt;		t3 =t3 +rt;
		t12=t4 -it;		t4 =t4 +it;

		rt =t19*c - t20*s;	t20=t20*c + t19*s;	t19=rt;	/* twiddle mpy by E^1	*/
		rt =t27*s - t28*c;	it =t28*s + t27*c;		/* twiddle mpy by E^3	*/
		t27=t19-rt;		t19=t19+rt;
		t28=t20-it;		t20=t20+it;
	  #if PFETCH
		addp = addr+p13;
		prefetch_p_doubles(addp);
	  #endif
		a[j1+p8 ]=t3+t19;	a[j2+p8 ]=t4+t20;
		a[j1+p9 ]=t3-t19;	a[j2+p9 ]=t4-t20;

		a[j1+p10]=t11-t28;	a[j2+p10]=t12+t27;	/* mpy by E^4=i is inlined here...	*/
		a[j1+p11]=t11+t28;	a[j2+p11]=t12-t27;
	  #if PFETCH
		addp = addr+p14;
		prefetch_p_doubles(addp);
	  #endif
		/*...Block 4: t7,15,23,31	*/
		rt =(t16+t15)*ISRT2;it =(t16-t15)*ISRT2;		/* twiddle mpy by -E^6 is here...	*/
		t15=t7 +rt;		t7 =t7 -rt;			/* ...and get E^6=(i-1)/sqrt by flipping signs here.	*/
		t16=t8 +it;		t8 =t8 -it;

		rt =t23*s - t24*c;	t24=t24*s + t23*c;	t23=rt;	/* twiddle mpy by E^3	*/
		rt =t31*c - t32*s;	it =t32*c + t31*s;		/* twiddle mpy by E^1 = -E^9...	*/
		t31=t23+rt;		t23=t23-rt;			/* ...and get E^9 by flipping signs here.	*/
		t32=t24+it;		t24=t24-it;			/* Note: t23+rt = t23*(s+1)	*/
	  #if PFETCH
		addp = addr+p15;
		prefetch_p_doubles(addp);
	  #endif
		a[j1+p12]=t7+t23;	a[j2+p12]=t8+t24;
		a[j1+p13]=t7-t23;	a[j2+p13]=t8-t24;

		a[j1+p14]=t15-t32;	a[j2+p14]=t16+t31;	/* mpy by E^4=i is inlined here...	*/
		a[j1+p15]=t15+t32;	a[j2+p15]=t16-t31;
	#endif	// USE_SCALAR_DFT_MACRO ?

#endif	/* USE_SSE2 */

	}	/* end for(j=_jstart[ithread]; j < _jhi[ithread]; j += 2) */

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart += nwt;
		jhi    += nwt;

		col += 16;
		co3 -= 16;
	}
}	/* end for(k=1; k <= khi; k++) */

