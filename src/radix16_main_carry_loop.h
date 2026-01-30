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

// This main loop is same for un-and-multithreaded, so stick into a header file
// (can't use a macro because of the #if-enclosed stuff).

for(int k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
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
#ifndef USE_SSE2
		j2 = j1 + RE_IM_STRIDE;
#endif

/*...The radix-32 DIT pass is here:	*/

#ifdef USE_SSE2

	add0 = &a[j1];
	SSE2_RADIX16_DIT_NOTWIDDLE(add0,p1,p2,p3,p4,p8,r00,r02,r04,r06,r08,r0a,r10,r18,isrt2,cc0);

#else	/* !USE_SSE2 */

  #if USE_SCALAR_DFT_MACRO

		// v20: Make in-place - need outputs back in array to support shifted-residue carry injection:
		RADIX_16_DIT(a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
		//			,a1p0r,a1p0i,a1p1r,a1p1i,a1p2r,a1p2i,a1p3r,a1p3i,a1p4r,a1p4i,a1p5r,a1p5i,a1p6r,a1p6i,a1p7r,a1p7i,a1p8r,a1p8i,a1p9r,a1p9i,a1pAr,a1pAi,a1pBr,a1pBi,a1pCr,a1pCi,a1pDr,a1pDi,a1pEr,a1pEi,a1pFr,a1pFi
					,a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
					,c,s)
  #else

	#ifdef USE_FGT61

	/*...Block 1: */
		jt = j1;		jp = j2;
		t1 =a[jt ];		t2 =a[jp ];					m1 =b[jt ];		m$ =b[jp ];
		rt =a[jt+p1 ];	it =a[jp+p1 ];				rm =b[jt+p1 ];	im =b[jp+p1 ];
		t3 =t1 -rt;		t1 =t1 +rt;					m3 =m1 -rm;		m1 =m1 +rm;		// 1,2 in 0,2b
		t4 =t2 -it;		t2 =t2 +it;					m4 =m$ -im;		m$ =m$ +im;		// 3,4 in -b,b

		t5 =a[jt+p2 ];	t6 =a[jp+p2 ];				m5 =b[jt+p2 ];	m6 =b[jp+p2 ];
		rt =a[jt+p3 ];	it =a[jp+p3 ];				rm =b[jt+p3 ];	im =b[jp+p3 ];
		t7 =t5 -rt;		t5 =t5 +rt;					m7 =m5 -rm;		m5 =m5 +rm;		// 5,6 in 0,2b
		t8 =t6 -it;		t6 =t6 +it;					m8 =m6 -im;		m6 =m6 +im;		// 7,8 in -b,b

		rt =t5;	t5 =t1 -rt ;	t1 =t1 +rt;			rm =m5;	m5 =m1 -rm ;	m1 =m1 +rm;	// 1,2 in   0,4b
		it =t6;	t6 =t2 -it ;	t2 =t2 +it;			im =m6;	m6 =m$ -im ;	m$ =m$ +im;	// 5,6 in -2b,2b

		rt =t7;	t7 =t3 -t8 ;	t3 =t3 +t8;			rm =m7;	m7 =m3 -m8 ;	m3 =m3 +m8;	// 3,4,7,8 all in -2b,2b
				t8 =t4 +rt ;	t4 =t4 -rt;					m8 =m4 +rm ;	m4 =m4 -rm;	// rest of 4-DFTs same

	/*...Block 2: */
		jt = j1 + p4;	jp = j2 + p4;
		t9 =a[jt    ];	t10=a[jp    ];				m9 =b[jt    ];	m10=b[jp    ];
		rt =a[jt+p1 ];	it =a[jp+p1 ];				rm =b[jt+p1 ];	im =b[jp+p1 ];
		t11=t9 -rt;		t9 =t9 +rt;					m11=m9 -rm;		m9 =m9 +rm;		//  9,10 in 0,2b
		t12=t10-it;		t10=t10+it;					m12=m10-im;		m10=m10+im;		// 11,12 in -b,b

		t13=a[jt+p2 ];	t14=a[jp+p2 ];				m13=b[jt+p2 ];	m14=b[jp+p2 ];
		rt =a[jt+p3 ];	it =a[jp+p3 ];				rm =b[jt+p3 ];	im =b[jp+p3 ];
		t15=t13-rt;		t13=t13+rt;					m15=m13-rm;		m13=m13+rm;		// 13,14 in 0,2b
		t16=t14-it;		t14=t14+it;					m16=m14-im;		m14=m14+im;		// 15,16 in -b,b

		rt =t13;	t13=t9 -rt ;	t9 =t9 +rt;		rm =m13;	m13=m9 -rm ;	m9 =m9 +rm;
		it =t14;	t14=t10-it ;	t10=t10+it;		im =m14;	m14=m10-im ;	m10=m10+im;

		rt =t15;	t15=t11-t16;	t11=t11+t16;	rm =m15;	m15=m11-m16;	m11=m11+m16;
					t16=t12+rt ;	t12=t12-rt;					m16=m12+rm ;	m12=m12-rm;

	/*...Block 3: */
		jt = j1 + p8;	jp = j2 + p8;
		t17=a[jt    ];	t18=a[jp    ];				m17=b[jt    ];	m18=b[jp    ];
		rt =a[jt+p1 ];	it =a[jp+p1 ];				rm =b[jt+p1 ];	im =b[jp+p1 ];
		t19=t17-rt;		t17=t17+rt;					m19=m17-rm;		m17=m17+rm;		// 17,18 in 0,2b
		t20=t18-it;		t18=t18+it;					m20=m18-im;		m18=m18+im;		// 19,20 in -b,b

		t21=a[jt+p2 ];	t22=a[jp+p2 ];				m21=b[jt+p2 ];	m22=b[jp+p2 ];
		rt =a[jt+p3 ];	it =a[jp+p3 ];				rm =b[jt+p3 ];	im =b[jp+p3 ];
		t23=t21-rt;		t21=t21+rt;					m23=m21-rm;		m21=m21+rm;	// 21,22 in 0,2b
		t24=t22-it;		t22=t22+it;					m24=m22-im;		m22=m22+im;	// 23,24 in -b,b

		rt =t21;	t21=t17-rt ;	t17=t17+rt;		rm =m21;	m21=m17-rm ;	m17=m17+rm;	// 17,18 in   0,4b
		it =t22;	t22=t18-it ;	t18=t18+it;		im =m22;	m22=m18-im ;	m18=m18+im;	// 21,22 in -2b,2b
													// Prior to qreduce(...+q4), all 4 outs here in -2b,2b:
		rt =t23;	t23=t19-t24;	t19=t19+t24;	rm =m23;	m23=qreduce(m19-m24+q4);	m19=qreduce(m19+m24+q4);
					t24=t20+rt ;	t20=t20-rt;					m24=qreduce(m20+rm +q4);	m20=qreduce(m20-rm +q4);
													// m19,20,23,24 all needed for CMUL, so reduce.
	/*...Block 4: */
		jt = j1 + p12;	jp = j2 + p12;
		t25=a[jt    ];	t26=a[jp    ];				m25=b[jt    ];	m26=b[jp    ];
		rt =a[jt+p1 ];	it =a[jp+p1 ];				rm =b[jt+p1 ];	im =b[jp+p1 ];
		t27=t25-rt;		t25=t25+rt;					m27=m25-rm;		m25=m25+rm;	// 25,26 in 0,2b
		t28=t26-it;		t26=t26+it;					m28=m26-im;		m26=m26+im;	// 27,28 in -b,b

		t29=a[jt+p2 ];	t30=a[jp+p2 ];				m29=b[jt+p2 ];	m30=b[jp+p2 ];
		rt =a[jt+p3 ];	it =a[jp+p3 ];				rm =b[jt+p3 ];	im =b[jp+p3 ];
		t31=t29-rt;		t29=t29+rt;					m31=m29-rm;		m29=m29+rm;	// 29,30 in 0,2b
		t32=t30-it;		t30=t30+it;					m32=m30-im;		m30=m30+im;	// 31,32 in -b,b

		rt =t29;	t29=t25-rt ;	t25=t25+rt;		rm =m29;	m29=m25-rm ;	m25=m25+rm;
		it =t30;	t30=t26-it ;	t26=t26+it;		im =m30;	m30=m26-im ;	m26=m26+im;
													// Prior to qreduce(...+q4), all 4 outs here in -2b,2b:
		rt =t31;	t31=t27-t32;	t27=t27+t32;	rm =m31;	m31=qreduce(m27-m32+q4);	m27=qreduce(m27+m32+q4);
					t32=t28+rt ;	t28=t28-rt;					m32=qreduce(m28+rm +q4);	m28=qreduce(m28-rm +q4);
													// m27,28,31,32 all needed for CMUL, so reduce.
	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	*/
													/*===============
													Input bounds on modular terms:
														-2b,2b: m3-8, 11-16, 21,22, 29,30
														  0,4b: m1,2, 9 ,10, 17,18, 25,26
													These are all reduced (in 0,b) because they are inputs to CMUL:
																m19,20,23,24,27,28,31,32
													===============*/
		// And since radix16_dit_pass1 is twiddleless version, use it as template for phase 2:
	/*...Block 1: t1,9,17,25	*/
		rt =t9 ;	t9 =t1 -rt;	t1 =t1 +rt;			rm =m9 ;	m9 =qreduce(m1 -rm+q4);	m1 =qreduce(m1 +rm);	// +:   0,8b -> 0,b
		it =t10;	t10=t2 -it;	t2 =t2 +it;			im =m10;	m10=qreduce(m$ -im+q4);	m$ =qreduce(m$ +im);	// -: -4b,4b -> 0,b

		rt =t25;	t25=t17-rt;	t17=t17+rt;			rm =m25;	m25=qreduce(m17-rm+q4);	m17=qreduce(m17+rm);	// same as above quartet
		it =t26;	t26=t18-it;	t18=t18+it;			im =m26;	m26=qreduce(m18-im+q4);	m18=qreduce(m18+im);

		a1p0r = t1 +t17;	a1p0i = t2 +t18;		b1p0r = qreduce(m1 +m17  );	b1p0i = qreduce(m$ +m18  );	// treat 0,b as = 0,q here
		a1p8r = t1 -t17;	a1p8i = t2 -t18;		b1p8r = qreduce(m1 -m17+q);	b1p8i = qreduce(m$ -m18+q);	// ~(1/2^58 odds of failure)

		a1p4r = t9 +t26;	a1p4i = t10-t25;		b1p4r = qreduce(m9 +m26  );	b1p4i = qreduce(m10-m25+q);
		a1pCr = t9 -t26;	a1pCi = t10+t25;		b1pCr = qreduce(m9 -m26+q);	b1pCi = qreduce(m10+m25  );

	/*...Block 3: t5,13,21,29	*/
		rt =t13;	t13=t5 -t14;	t5 =t5 +t14;	rm =m13;m13=qreduce(m5-m14+q4);	m5 =qreduce(m5 +m14+q4);	// all 4 outs in -4b,4b;
					t14=t6 +rt;		t6 =t6 -rt;				m14=qreduce(m6+rm +q4);	m6 =qreduce(m6 -rm +q4);	// reduce all 4 to 0,b.

		rt =(t22+t21)*ISRT2;t22=(t22-t21)*ISRT2;	rm = mul_i2(m22+m21+q4);	m22= mul_i2(m22-m21+q4);	m21=rm;
t21=rt;	rt =(t29-t30)*ISRT2;it =(t29+t30)*ISRT2;	rm = mul_i2(m29-m30+q4);	im = mul_i2(m29+m30+q4);	// m21,22,rm,im in 0,b30
		t29=t21+rt;		t21=t21-rt;					m29 = m21+rm;		m21 = m21-rm;	// m21,22 in -b30,b30
		t30=t22+it;		t22=t22-it;					m30 = m22+im;		m22 = m22-im;	// m29,30 in 0,2*b30

		a1p2r = t5 +t21;	a1p2i = t6 +t22;		b1p2r = qreduce(m5 +m21+q2);	b1p2i = qreduce(m6 +m22+q2);	// + in -b30,b+b30
		a1pAr = t5 -t21;	a1pAi = t6 -t22;		b1pAr = qreduce(m5 -m21+q2);	b1pAi = qreduce(m6 -m22+q2);	// - same; reduce all to 0,b

		a1p6r = t13+t30;	a1p6i = t14-t29;		b1p6r = qreduce(m13+m30   );	b1p6i = qreduce(m14-m29+q3);	// + in 0,b+2*b30
		a1pEr = t13-t30;	a1pEi = t14+t29;		b1pEr = qreduce(m13-m30+q3);	b1pEi = qreduce(m14+m29   );	// - in -2*b30, b

	/*...Block 2: t3,11,19,27	*/
		rt =(t12+t11)*ISRT2;it =(t12-t11)*ISRT2;	rm = mul_i2(m12+m11+q4);im = mul_i2(m12-m11+q4);	// 0,b30
		t11 = t3 -rt;		t3 = t3 +rt;			m11 = m3 -rm;		m3 = m3 +rm;	//  3, 4 in -2b,2b+b30
		t12 = t4 -it;		t4 = t4 +it;			m12 = m4 -im;		m4 = m4 +im;	// 11,12 in -2b-b30,2b
											// Remember, mod-root [cm,sm] premultiplied by 8, thus negate by subbing from q*8:
		rt =t19*c + t20*s;	t20=t20*c - t19*s;		cmul_modq8(m19,m20, cm,q8-sm, &m19,&m20);
t19=rt;	rt =t27*s + t28*c;	it =t28*s - t27*c;		cmul_modq8(m27,m28, sm,q8-cm, &rm ,&im );
		t27 = t19-rt;		t19 = t19+rt;			m27 =qreduce(m19-rm+q4);	m19 =qreduce(m19+rm);	// +:   0,8q ==> 0,b
		t28 = t20-it;		t20 = t20+it;			m28 =qreduce(m20-im+q4);	m20 =qreduce(m20+im);	// -: -4q,4q ==> 0,b
						// NB: +=3q for rm,im below too large, since odds of 4b+b30+3q overflowing 64-bits much larger than those of +2q not being enough to make result >= 0; analogously for m3,4)
		a1p1r = t3 +t19;	a1p1i = t4 +t20;		b1p1r = qreduce(m3 +m19+q3);	b1p1i = qreduce(m4 +m20+q3);	// + in [-2b,2b+b30] + [0,b] = [-2b,3b+b30], add 3q and reduce
		a1p9r = t3 -t19;	a1p9i = t4 -t20;		b1p9r = qreduce(m3 -m19+q4);	b1p9i = qreduce(m4 -m20+q4);	// - in [-2b,2b+b30] - [0,b] = [-3b,2b+b30], add 4q and reduce

		a1p5r = t11+t28;	a1p5i = t12-t27;		b1p5r = qreduce(m11+m28+q4);	b1p5i = qreduce(m12-m27+q4);	// +- in [-2b-b30,2b] +- [0,b] = -2b-b30,3b;
		a1pDr = t11-t28;	a1pDi = t12+t27;		b1pDr = qreduce(m11-m28+q4);	b1pDi = qreduce(m12+m27+q4);	// add 4q and reduce.

	/*...Block 4: t7,15,23,31	*/
		rt =(t15-t16)*ISRT2;it =(t15+t16)*ISRT2;	rm = mul_i2(m15-m16+q4);im = mul_i2(m15+m16+q4);	// 0,b30
		t15 = t7 +rt;			t7 = t7 -rt;		m15 = m7 +rm;			m7 = m7 -rm;	//  7, 8 in -2b-b30,2b
		t16 = t8 +it;			t8 = t8 -it;		m16 = m8 +im;			m8 = m8 -im;	// 15,16 in -2b,2b+b30

		rt =t23*s + t24*c;	t24=t24*s - t23*c;		cmul_modq8(m23,m24, sm,q8-cm, &m23,&m24);
t23=rt;	rt =t31*c + t32*s;	it =t32*c - t31*s;		cmul_modq8(m31,m32, cm,q8-sm, &rm ,&im );
		t31 = t23+rt;			t23 = t23-rt;		m31 =qreduce(m23+rm);	m23 =qreduce(m23-rm+q4);// +:   0,8q ==> 0,b
		t32 = t24+it;			t24 = t24-it;		m32 =qreduce(m24+im);	m24 =qreduce(m24-im+q4);// -: -4q,4q ==> 0,b

		a1p3r = t7 +t23;	a1p3i = t8 +t24;		b1p3r = qreduce(m7 +m23+q4);	b1p3i = qreduce(m8 +m24+q4);	// + in [-2b-b30,2b] + [0,b] = -2b-b30,3b, add 4q and reduce
		a1pBr = t7 -t23;	a1pBi = t8 -t24;		b1pBr = qreduce(m7 -m23+q5);	b1pBi = qreduce(m8 -m24+q5);	// - in [-2b-b30,2b] - [0,b] = -3b-b30,2b, add 5q and reduce

		a1p7r = t15+t32;	a1p7i = t16-t31;		b1p7r = qreduce(m15+m32+q3);	b1p7i = qreduce(m16-m31+q4);	// + in [-2b,2b+b30] + [0,b] = -2b,3b+b30, add 3q and reduce
		a1pFr = t15-t32;	a1pFi = t16+t31;		b1pFr = qreduce(m15-m32+q4);	b1pFi = qreduce(m16+m31+q3);	// - in [-2b,2b+b30] - [0,b] = -3b,2b+b30, add 4q and reduce

			/**********************************************/
	#else	// USE_FGT61 = False; Basic scalar-double mode:
			/**********************************************/
		#error 16-DFT outputs-into-scalar-temps no longer supported!
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

	#endif	// USE_FGT61 ?

  #endif	// USE_SCALAR_DFT_MACRO ?

#endif	/* USE_SSE2 */

	/*...and combine those to complete the radix-16 transform and do the carries. Since the outputs would
	normally be getting dispatched to 16 separate blocks of the A-array, we need 16 separate carries.	*/

	if(MODULUS_TYPE == MODULUS_TYPE_GENFFTMUL)
	{	//			Indices in rightmost col are debug-usage only: vvv
	#ifdef USE_SSE2
		#warning Experimental GENFFTMUL carry-macros only available in non-SIMD build mode.
	#else
		genfftmul_carry_norm_pow2_errcheck(a[j1       ],a[j2       ],cy_r0,cy_i0,0x0);
		genfftmul_carry_norm_pow2_errcheck(a[j1+p1    ],a[j2+p1    ],cy_r1,cy_i1,0x1);
		genfftmul_carry_norm_pow2_errcheck(a[j1+p2    ],a[j2+p2    ],cy_r2,cy_i2,0x2);
		genfftmul_carry_norm_pow2_errcheck(a[j1+p3    ],a[j2+p3    ],cy_r3,cy_i3,0x3);
		genfftmul_carry_norm_pow2_errcheck(a[j1   +p4 ],a[j2   +p4 ],cy_r4,cy_i4,0x4);
		genfftmul_carry_norm_pow2_errcheck(a[j1+p1+p4 ],a[j2+p1+p4 ],cy_r5,cy_i5,0x5);
		genfftmul_carry_norm_pow2_errcheck(a[j1+p2+p4 ],a[j2+p2+p4 ],cy_r6,cy_i6,0x6);
		genfftmul_carry_norm_pow2_errcheck(a[j1+p3+p4 ],a[j2+p3+p4 ],cy_r7,cy_i7,0x7);
		genfftmul_carry_norm_pow2_errcheck(a[j1   +p8 ],a[j2   +p8 ],cy_r8,cy_i8,0x8);
		genfftmul_carry_norm_pow2_errcheck(a[j1+p1+p8 ],a[j2+p1+p8 ],cy_r9,cy_i9,0x9);
		genfftmul_carry_norm_pow2_errcheck(a[j1+p2+p8 ],a[j2+p2+p8 ],cy_rA,cy_iA,0xA);
		genfftmul_carry_norm_pow2_errcheck(a[j1+p3+p8 ],a[j2+p3+p8 ],cy_rB,cy_iB,0xB);
		genfftmul_carry_norm_pow2_errcheck(a[j1   +p12],a[j2   +p12],cy_rC,cy_iC,0xC);
		genfftmul_carry_norm_pow2_errcheck(a[j1+p1+p12],a[j2+p1+p12],cy_rD,cy_iD,0xD);
		genfftmul_carry_norm_pow2_errcheck(a[j1+p2+p12],a[j2+p2+p12],cy_rE,cy_iE,0xE);
		genfftmul_carry_norm_pow2_errcheck(a[j1+p3+p12],a[j2+p3+p12],cy_rF,cy_iF,0xF);
	#endif
	}
	else if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		// Check if current index-interval contains the target index for rotated-residue carry injection.
		// In data-init we set target_idx = -1 on wraparound-carry mini-pass, so if() only taken on full pass:
		if(target_idx == j) {
		#ifdef USE_SSE2
			addr = (double *)r00 + target_set;
			*addr += target_cy*(n>>1);	// target_cy = [-2 << within-word-shift]*[DWT weight]*n/2, i.e. includes fwd DWT weight and n/2 factor
		#else
			// target_set in [0,2*RADIX); tidx_mod_stride [even|odd] means shifted-carry goes into [Re|Im] part of the complex FFT datum:
			l = target_set&1;	target_set >>= 1;
			a[j1+poff[target_set>>2]+p0123[target_set&3]+l] += target_cy*(n>>1);
		#endif
			target_idx = -1;
		}

	#ifdef USE_AVX
		// For AVX512-and-beyond we support only the fast Mers-carry macros.
		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];
		/* ptr to local storage for the doubled wtl,wtn terms: */
	  #ifdef USE_AVX512
		tmp = half_arr +  64;	// No lookup-tables used in avx-512; instead use opmasked conditional-doubling; 1st 32 slots hold outputs of wtsinit call
	  #else
		tmp = half_arr + 128;	// 1st 64 slots are basic-4 LUTs, next 32 are the additional 2 LOACC LUTs, next 32 hold outputs of wtsinit call
	  #endif
		l= j & (nwt-1);						// These rcol wts-terms are for individual-double-broadcast-to-full-vector-width,
		n_minus_sil  ->d0 = n-si[l  ];		tmp->d0 = wt0[    l  ];	// hence the mixing of fwd/inv wts, which is normally taboo.
		n_minus_silp1->d0 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d0 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d0 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+2) & (nwt-1);					++tmp;	/* Get ready for next 4 weights-related doubles... */
		n_minus_sil  ->d1 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d1 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d1 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d1 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+4) & (nwt-1);					++tmp;
		n_minus_sil  ->d2 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d2 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d2 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d2 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;

		l= (j+6) & (nwt-1);					++tmp;
		n_minus_sil  ->d3 = n-si[l  ];		tmp->d0 = wt0[    l  ];
		n_minus_silp1->d3 = n-si[l+1];		tmp->d1 = wt0[nwt-l  ]*scale;
		sinwt        ->d3 = si[nwt-l  ];	tmp->d2 = wt0[    l+1];
		sinwtm1      ->d3 = si[nwt-l-1];	tmp->d3 = wt0[nwt-l-1]*scale;
	  #ifdef USE_AVX512
		l= (j+8) & (nwt-1);					tmp -= 3;	// Reset to same tmp-startval as above, now copy data into d4-7 slots of vec_dbl
		n_minus_sil  ->d4 = n-si[l  ];		tmp->d4 = wt0[    l  ];
		n_minus_silp1->d4 = n-si[l+1];		tmp->d5 = wt0[nwt-l  ]*scale;
		sinwt        ->d4 = si[nwt-l  ];	tmp->d6 = wt0[    l+1];
		sinwtm1      ->d4 = si[nwt-l-1];	tmp->d7 = wt0[nwt-l-1]*scale;

		l= (j+10) & (nwt-1);				++tmp;
		n_minus_sil  ->d5 = n-si[l  ];		tmp->d4 = wt0[    l  ];
		n_minus_silp1->d5 = n-si[l+1];		tmp->d5 = wt0[nwt-l  ]*scale;
		sinwt        ->d5 = si[nwt-l  ];	tmp->d6 = wt0[    l+1];
		sinwtm1      ->d5 = si[nwt-l-1];	tmp->d7 = wt0[nwt-l-1]*scale;

		l= (j+12) & (nwt-1);				++tmp;
		n_minus_sil  ->d6 = n-si[l  ];		tmp->d4 = wt0[    l  ];
		n_minus_silp1->d6 = n-si[l+1];		tmp->d5 = wt0[nwt-l  ]*scale;
		sinwt        ->d6 = si[nwt-l  ];	tmp->d6 = wt0[    l+1];
		sinwtm1      ->d6 = si[nwt-l-1];	tmp->d7 = wt0[nwt-l-1]*scale;

		l= (j+14) & (nwt-1);				++tmp;
		n_minus_sil  ->d7 = n-si[l  ];		tmp->d4 = wt0[    l  ];
		n_minus_silp1->d7 = n-si[l+1];		tmp->d5 = wt0[nwt-l  ]*scale;
		sinwt        ->d7 = si[nwt-l  ];	tmp->d6 = wt0[    l+1];
		sinwtm1      ->d7 = si[nwt-l-1];	tmp->d7 = wt0[nwt-l-1]*scale;
	  #endif
	 #ifndef USE_AVX512	// Power-of-2 radix means LOACC-mode only for AVX-512
	  if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {	// LOACC with tunable DWT-weights chaining
	 #endif
		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
		// This is based on the AVX 4-way init macro:
	   #ifdef CARRY_8_WAY
		AVX_cmplx_carry_fast_pow2_wtsinit_X8(add1,add2,add3, bjmodn0, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_nm1)
	   #else
		AVX_cmplx_carry_fast_pow2_wtsinit_X4(add1,add2,add3, bjmodn0, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, sse_bw,sse_nm1)
	   #endif

		i = (!j);
		addr = &prp_mult;

	  #ifdef CARRY_8_WAY

	   #ifdef USE_AVX512

		// Each carry macro call also processes 8 prefetches of main-array data:
		add0 = a + j1 + pfetch_dist;
		AVX_cmplx_carry_fast_pow2_errcheck_X8(r00,cy_r0,bjmodn0,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3,p4, addr); i = 0;
		add0 = a + j1 + pfetch_dist + p8;
		AVX_cmplx_carry_fast_pow2_errcheck_X8(r10,cy_r8,bjmodn8,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3,p4, addr);

	   #elif defined(USE_AVX2) && defined(GCC_5PLUS)
			// gcc 4.x !support the needed AVX2 8-way int instructions (while still being fine for for the floating-FMA
			// used for the FFT), so require an added compile-time define to enable this version of 8-way. Based on my
			// [Has,Broad]well timings this is no faster than the version using half-width 128-bit arithmetic for the
			// integer math attendant to Mersenne-mod IBDWT, so leave these calls around just in this one radix-16
			// carry routine which I used to do my comparative timing tests:

		// Each carry macro call also processes 8 prefetches of main-array data:
		add0 = a + j1 + pfetch_dist;
		AVX2_cmplx_carry_fast_pow2_errcheck_X8(r00,cy_r0,cy_r4,bjmodn0,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3,p4, addr); i = 0;
		add0 = a + j1 + pfetch_dist + p8 ;	// poff[] = p0,4,8,...
		AVX2_cmplx_carry_fast_pow2_errcheck_X8(r10,cy_r8,cy_rC,bjmodn8,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3,p4, addr);

	   #else	// Version targeting AVX-only and AVX2-under-gcc 4.x, which has FMA3 support, but not AVX2 8-way-int support

		// Each carry macro call also processes 8 prefetches of main-array data:
		add0 = a + j1 + pfetch_dist;
		AVX_cmplx_carry_fast_pow2_errcheck_X8(r00,cy_r0,cy_r4,bjmodn0,bjmodn4,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3,p4, addr); i = 0;
		add0 = a + j1 + pfetch_dist + p8 ;	// poff[] = p0,4,8,...
		AVX_cmplx_carry_fast_pow2_errcheck_X8(r10,cy_r8,cy_rC,bjmodn8,bjmodnC,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3,p4, addr);

	   #endif	// ifdef GCC_5PLUS ?

	  #else	// USE_AVX, LOACC, 4-way carry:

		// Each carry macro call also processes 4 prefetches of main-array data:
		add0 = a + j1 + pfetch_dist;
		AVX_cmplx_carry_fast_pow2_errcheck_X4(r00,cy_r0,bjmodn0,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3, addr); i = 0;
		add0 = a + j1 + pfetch_dist + p4 ;	// poff[] = p0,4,8,...
		AVX_cmplx_carry_fast_pow2_errcheck_X4(r08,cy_r4,bjmodn4,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3, addr);
		add0 = a + j1 + pfetch_dist + p8 ;	// poff[] = p0,4,8,...
		AVX_cmplx_carry_fast_pow2_errcheck_X4(r10,cy_r8,bjmodn8,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3, addr);
		add0 = a + j1 + pfetch_dist + p12;	// poff[] = p0,4,8,...
		AVX_cmplx_carry_fast_pow2_errcheck_X4(r18,cy_rC,bjmodnC,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3, addr);

	  #endif	// USE_AVX, LOACC, (8-way or 4-way carry) ?

	 #ifndef USE_AVX512
	  } else {	// HiACC:

		// Each AVX carry macro call also processes 4 prefetches of main-array data
		i = (!j);
		addr = &prp_mult;
		add0 = a + j1 + pfetch_dist;
		AVX_cmplx_carry_norm_pow2_errcheck_X4(r00,add1,add2,add3,cy_r0,bjmodn0,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3, addr); i = 0;
		add0 = a + j1 + pfetch_dist + p4 ;	// poff[] = p0,4,8,...
		AVX_cmplx_carry_norm_pow2_errcheck_X4(r08,add1,add2,add3,cy_r4,bjmodn4,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3, addr);
		add0 = a + j1 + pfetch_dist + p8 ;	// poff[] = p0,4,8,...
		AVX_cmplx_carry_norm_pow2_errcheck_X4(r10,add1,add2,add3,cy_r8,bjmodn8,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3, addr);
		add0 = a + j1 + pfetch_dist + p12;	// poff[] = p0,4,8,...
		AVX_cmplx_carry_norm_pow2_errcheck_X4(r18,add1,add2,add3,cy_rC,bjmodnC,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3, addr);

		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

	  }	// USE_AVX: (8-way or 4-way LOACC) or (4-way HIACC) ?
	 #endif

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

	 #ifdef USE_ARM_V8_SIMD
	  if(1) {	// No HIACC mode for ARMv8
	 #else
	  if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {	// LOACC with tunable DWT-weights chaining
	 #endif
		/************ wt_re,wi_re,wt_im,wi_im inits: **************/
	   #if 1	// SSE2 inlines-asm-ized version of the scalar-double macro-call sequence in the #else below:

		uint32 k0,k1,k2,k3;
		l = j & (nwt-1);
		n_minus_sil   = n-si[l  ];
		n_minus_silp1 = n-si[l+1];
		sinwt   = si[nwt-l  ];
		sinwtm1 = si[nwt-l-1];
		wtl     = wt0[    l  ];
		wtn     = wt0[nwt-l  ]*scale;
		wtlp1   = wt0[    l+1];
		wtnm1   = wt0[nwt-l-1]*scale;

		ctmp = (struct complex *)half_arr + 24;	/* ptr to local storage for the doubled wtl,wtn terms: */
		// (j)-data occupy the 8 xmm-sized slots above the 16 ones used by the fixed auxiliary-data, and overwrite these inits:
		ctmp->re = ctmp->im = wtl;		ctmp += 2;
		ctmp->re = ctmp->im = wtn;		ctmp += 2;
		ctmp->re = ctmp->im = wtlp1;	ctmp += 2;
		ctmp->re = ctmp->im = wtnm1;

		l = (j+2) & (nwt-1);
		k0 = n-si[l  ];
		k1 = n-si[l+1];
		k2 = si[nwt-l  ];
		k3 = si[nwt-l-1];
		wtl     = wt0[    l  ];
		wtn     = wt0[nwt-l  ]*scale;
		wtlp1   = wt0[    l+1];
		wtnm1   = wt0[nwt-l-1]*scale;

		ctmp = (struct complex *)half_arr + 32;	// (j+2) data start at ctmp + 8
		ctmp->re = ctmp->im = wtl;		ctmp += 2;
		ctmp->re = ctmp->im = wtn;		ctmp += 2;
		ctmp->re = ctmp->im = wtlp1;	ctmp += 2;
		ctmp->re = ctmp->im = wtnm1;

		add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).
		SSE2_cmplx_carry_fast_pow2_wtsinit(add1,add2,add3, bjmodn0, half_arr,sign_mask, n_minus_sil,n_minus_silp1,sinwt,sinwtm1, k0,k1,k2,k3, sse_bw,sse_nm1)

	   #else

		uint32 i0,i1,i2,i3;
		double wtA,wtB,wtC;
		const double one_half[3] = {1.0, 0.5, 0.25};

		// First quartet of roots corr. to calling cmplx_carry_norm_pow2_errcheck with loop index j and set = 0-3.

		l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 32 carries, so precompute	*/
		n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
		n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
		sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		ctmp = (struct complex *)half_arr + 24;	/* ptr to local storage for the doubled wtl,wtn terms: */
			// (j)-data occupy the 8 xmm-sized slots above the 16 ones used by the fixed auxiliary-data.

		// Since the wtsinit macro calls here are chained and linked via bjmodn[0-3]-updates, use i to store those:
		i0 = *bjmodn0;	cmplx_carry_fast_pow2_wtsinit((ctmp+0)->re,(ctmp+2)->re,(ctmp+4)->re,(ctmp+6)->re,i0,0)	// set 0
		i1 = *bjmodn1;	cmplx_carry_fast_pow2_wtsinit((ctmp+0)->im,(ctmp+2)->im,(ctmp+4)->im,(ctmp+6)->im,i1,1)	// set 1
		i2 = *bjmodn2;	cmplx_carry_fast_pow2_wtsinit((ctmp+1)->re,(ctmp+3)->re,(ctmp+5)->re,(ctmp+7)->re,i2,2)	// set 2
		i3 = *bjmodn3;	cmplx_carry_fast_pow2_wtsinit((ctmp+1)->im,(ctmp+3)->im,(ctmp+5)->im,(ctmp+7)->im,i3,3)	// set 3

		// 2nd quartet of roots corr. to calling cmplx_carry_norm_pow2_errcheck with loop index j+2 and set = 0-3.

		l= (j+2) & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 32 carries, so precompute	*/
		n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
		n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
		sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		// Since use wt1-array in the wtsinit macro, need to fiddle this here:
		co2 = co3;	// For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					// (and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).

		ctmp = (struct complex *)half_arr + 32;	/* (j+2) data occupy the 8 xmm-sized slots above the (j) ones */

		// Since the wtsinit macro calls here are chained and linked via bjmodn[0-3]-updates, use i to store those:
						cmplx_carry_fast_pow2_wtsinit((ctmp+0)->re,(ctmp+2)->re,(ctmp+4)->re,(ctmp+6)->re,i0,0)	// set 0
						cmplx_carry_fast_pow2_wtsinit((ctmp+0)->im,(ctmp+2)->im,(ctmp+4)->im,(ctmp+6)->im,i1,1)	// set 1
						cmplx_carry_fast_pow2_wtsinit((ctmp+1)->re,(ctmp+3)->re,(ctmp+5)->re,(ctmp+7)->re,i2,2)	// set 2
						cmplx_carry_fast_pow2_wtsinit((ctmp+1)->im,(ctmp+3)->im,(ctmp+5)->im,(ctmp+7)->im,i3,3)	// set 3

	   #endif	// SSE2 inline-asm or above scalar-double macro-call sequence

		i = (!j);
		addr = &prp_mult;
		// Each SSE2 LOACC carry macro call also processes 4 prefetches of main-array data:
		add0 = a + j1 + pfetch_dist;
		SSE2_cmplx_carry_fast_pow2_errcheck(r00,cy_r0,cy_r2,bjmodn0,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3, addr); i = 0;
		add0 += p4;
		SSE2_cmplx_carry_fast_pow2_errcheck(r08,cy_r4,cy_r6,bjmodn4,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3, addr);
		add0 = a + j1 + pfetch_dist + p8;
		SSE2_cmplx_carry_fast_pow2_errcheck(r10,cy_r8,cy_rA,bjmodn8,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3, addr);
		add0 += p4;
		SSE2_cmplx_carry_fast_pow2_errcheck(r18,cy_rC,cy_rE,bjmodnC,half_arr,i,sign_mask,sse_bw,sse_nm1,sse_sw, add0,p1,p2,p3, addr);

	  } else {	// HiACC:

		l= j & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
		n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
		n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
		sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		ctmp = (struct complex *)half_arr + 24;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		add1 = &wt1[col  ];	/* Don't use add0 here, to avoid need to reload main-array address */
		add2 = &wt1[co2-1];
		add3 = &wt1[co3-1];

		// Each SSE2 carry macro call also processes 2 prefetches of main-array data
		i = (!j);
		addr = &prp_mult;
		add0 = a + j1 + pfetch_dist;
		SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r00,add1,add2,add3,cy_r0,cy_r2,bjmodn0,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p1, addr); i = 0;
		add0 += p4;
		SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r08,add1,add2,add3,cy_r4,cy_r6,bjmodn4,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p1, addr);
		add0 = a + j1 + pfetch_dist + p8;
		SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r10,add1,add2,add3,cy_r8,cy_rA,bjmodn8,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p1, addr);
		add0 += p4;
		SSE2_cmplx_carry_norm_pow2_errcheck1_2B(r18,add1,add2,add3,cy_rC,cy_rE,bjmodnC,half_arr,i,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p1, addr);

		l= (j+2) & (nwt-1);			/* We want (S*J mod N) - SI(L) for all 16 carries, so precompute	*/
		n_minus_sil   = n-si[l  ];		/* N - SI(L) and for each J, find N - (B*J mod N) - SI(L)		*/
		n_minus_silp1 = n-si[l+1];		/* For the inverse weight, want (S*(N - J) mod N) - SI(NWT - L) =	*/
		sinwt   = si[nwt-l  ];		/*	= N - (S*J mod N) - SI(NWT - L) = (B*J mod N) - SI(NWT - L).	*/
		sinwtm1 = si[nwt-l-1];

		wtl     =wt0[    l  ];
		wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
		wtlp1   =wt0[    l+1];
		wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

		ctmp = (struct complex *)half_arr + 24;	/* ptr to local storage for the doubled wtl,wtn terms: */
		ctmp->re = wtl;		ctmp->im = wtl;	++ctmp;
		ctmp->re = wtn;		ctmp->im = wtn;	++ctmp;
		ctmp->re = wtlp1;	ctmp->im = wtlp1;++ctmp;
		ctmp->re = wtnm1;	ctmp->im = wtnm1;

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/
		add1 = &wt1[col  ];
		add2 = &wt1[co2-1];

		// Each SSE2 carry macro call also processes 2 prefetches of main-array data
		add0 = a + j1 + pfetch_dist;
		SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r00,add1,add2,cy_r0,cy_r2,bjmodn0,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p2,p3, addr);
		add0 += p4;
		SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r08,add1,add2,cy_r4,cy_r6,bjmodn4,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p2,p3, addr);
		add0 = a + j1 + pfetch_dist + p8;
		SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r10,add1,add2,cy_r8,cy_rA,bjmodn8,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p2,p3, addr);
		add0 += p4;
		SSE2_cmplx_carry_norm_pow2_errcheck2_2B(r18,add1,add2,cy_rC,cy_rE,bjmodnC,half_arr,n_minus_silp1,n_minus_sil,sign_mask,sinwt,sinwtm1,sse_bw,sse_nm1,sse_sw, add0,p2,p3, addr);

	  }	// LOACC or HIACC?

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

	#ifdef USE_FGT61

if(!j) {
	printf("J = 0, carry-step INputs:\n");
	printf("a1p0r,a1p0i, b1p0r,b1p0i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p0r,a1p0i, b1p0r,b1p0i);
	printf("a1p1r,a1p1i, b1p1r,b1p1i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p1r,a1p1i, b1p1r,b1p1i);
	printf("a1p2r,a1p2i, b1p2r,b1p2i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p2r,a1p2i, b1p2r,b1p2i);
	printf("a1p3r,a1p3i, b1p3r,b1p3i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p3r,a1p3i, b1p3r,b1p3i);
	printf("a1p4r,a1p4i, b1p4r,b1p4i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p4r,a1p4i, b1p4r,b1p4i);
	printf("a1p5r,a1p5i, b1p5r,b1p5i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p5r,a1p5i, b1p5r,b1p5i);
	printf("a1p6r,a1p6i, b1p6r,b1p6i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p6r,a1p6i, b1p6r,b1p6i);
	printf("a1p7r,a1p7i, b1p7r,b1p7i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p7r,a1p7i, b1p7r,b1p7i);
	printf("a1p8r,a1p8i, b1p8r,b1p8i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p8r,a1p8i, b1p8r,b1p8i);
	printf("a1p9r,a1p9i, b1p9r,b1p9i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p9r,a1p9i, b1p9r,b1p9i);
	printf("a1pAr,a1pAi, b1pAr,b1pAi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1pAr,a1pAi, b1pAr,b1pAi);
	printf("a1pBr,a1pBi, b1pBr,b1pBi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1pBr,a1pBi, b1pBr,b1pBi);
	printf("a1pCr,a1pCi, b1pCr,b1pCi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1pCr,a1pCi, b1pCr,b1pCi);
	printf("a1pDr,a1pDi, b1pDr,b1pDi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1pDr,a1pDi, b1pDr,b1pDi);
	printf("a1pEr,a1pEi, b1pEr,b1pEi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1pEr,a1pEi, b1pEr,b1pEi);
	printf("a1pFr,a1pFi, b1pFr,b1pFi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1pFr,a1pFi, b1pFr,b1pFi);
}
if(!j) {
	if(full_pass)printf("\n");
	else
		printf("Wrapper!\n");
	printf("iter %2u [full-pass = %u]: a01 IN : %20.10e, %20.10e, cy = %20.10e\n",iter,full_pass,a1p0r,a1p0i,cy_r0);
}
// function form for step-thru debugging:
		// Mod-power-of-2 here used for *= 1/n2 if full_pass, *= 1 otherwise
		int shift = 0, ipow = (-full_pass) & l2_n2;
if(!j)
	printf("IPOW = %d, [l2_n2 = %d]\n",ipow,l2_n2);
		mixed_carry(sw,bw,nm1,bits,j,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,base,baseinv,wt1,wtl,wtlp1,wtn,wtnm1,&maxerr, &a1p0r,&a1p0i,&b1p0r,&b1p0i,&cy_r0,bjmodn0,0x0, shift, ipow);
		mixed_carry(sw,bw,nm1,bits,j,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,base,baseinv,wt1,wtl,wtlp1,wtn,wtnm1,&maxerr, &a1p1r,&a1p1i,&b1p1r,&b1p1i,&cy_r1,bjmodn1,0x1, shift, ipow);
		mixed_carry(sw,bw,nm1,bits,j,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,base,baseinv,wt1,wtl,wtlp1,wtn,wtnm1,&maxerr, &a1p2r,&a1p2i,&b1p2r,&b1p2i,&cy_r2,bjmodn2,0x2, shift, ipow);
		mixed_carry(sw,bw,nm1,bits,j,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,base,baseinv,wt1,wtl,wtlp1,wtn,wtnm1,&maxerr, &a1p3r,&a1p3i,&b1p3r,&b1p3i,&cy_r3,bjmodn3,0x3, shift, ipow);
		mixed_carry(sw,bw,nm1,bits,j,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,base,baseinv,wt1,wtl,wtlp1,wtn,wtnm1,&maxerr, &a1p4r,&a1p4i,&b1p4r,&b1p4i,&cy_r4,bjmodn4,0x4, shift, ipow);
		mixed_carry(sw,bw,nm1,bits,j,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,base,baseinv,wt1,wtl,wtlp1,wtn,wtnm1,&maxerr, &a1p5r,&a1p5i,&b1p5r,&b1p5i,&cy_r5,bjmodn5,0x5, shift, ipow);
		mixed_carry(sw,bw,nm1,bits,j,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,base,baseinv,wt1,wtl,wtlp1,wtn,wtnm1,&maxerr, &a1p6r,&a1p6i,&b1p6r,&b1p6i,&cy_r6,bjmodn6,0x6, shift, ipow);
		mixed_carry(sw,bw,nm1,bits,j,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,base,baseinv,wt1,wtl,wtlp1,wtn,wtnm1,&maxerr, &a1p7r,&a1p7i,&b1p7r,&b1p7i,&cy_r7,bjmodn7,0x7, shift, ipow);
		mixed_carry(sw,bw,nm1,bits,j,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,base,baseinv,wt1,wtl,wtlp1,wtn,wtnm1,&maxerr, &a1p8r,&a1p8i,&b1p8r,&b1p8i,&cy_r8,bjmodn8,0x8, shift, ipow);
		mixed_carry(sw,bw,nm1,bits,j,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,base,baseinv,wt1,wtl,wtlp1,wtn,wtnm1,&maxerr, &a1p9r,&a1p9i,&b1p9r,&b1p9i,&cy_r9,bjmodn9,0x9, shift, ipow);
		mixed_carry(sw,bw,nm1,bits,j,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,base,baseinv,wt1,wtl,wtlp1,wtn,wtnm1,&maxerr, &a1pAr,&a1pAi,&b1pAr,&b1pAi,&cy_rA,bjmodnA,0xA, shift, ipow);
		mixed_carry(sw,bw,nm1,bits,j,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,base,baseinv,wt1,wtl,wtlp1,wtn,wtnm1,&maxerr, &a1pBr,&a1pBi,&b1pBr,&b1pBi,&cy_rB,bjmodnB,0xB, shift, ipow);
		mixed_carry(sw,bw,nm1,bits,j,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,base,baseinv,wt1,wtl,wtlp1,wtn,wtnm1,&maxerr, &a1pCr,&a1pCi,&b1pCr,&b1pCi,&cy_rC,bjmodnC,0xC, shift, ipow);
		mixed_carry(sw,bw,nm1,bits,j,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,base,baseinv,wt1,wtl,wtlp1,wtn,wtnm1,&maxerr, &a1pDr,&a1pDi,&b1pDr,&b1pDi,&cy_rD,bjmodnD,0xD, shift, ipow);
		mixed_carry(sw,bw,nm1,bits,j,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,base,baseinv,wt1,wtl,wtlp1,wtn,wtnm1,&maxerr, &a1pEr,&a1pEi,&b1pEr,&b1pEi,&cy_rE,bjmodnE,0xE, shift, ipow);
		mixed_carry(sw,bw,nm1,bits,j,col,co2,co3,n_minus_sil,n_minus_silp1,sinwt,sinwtm1,base,baseinv,wt1,wtl,wtlp1,wtn,wtnm1,&maxerr, &a1pFr,&a1pFi,&b1pFr,&b1pFi,&cy_rF,bjmodnF,0xF, shift, ipow);
/* macro form, once have dbugged code:
		mixed_carry(a1p0r,a1p0i,b1p0r,b1p0i,cy_r0,bjmodn0,0x0);
		mixed_carry(a1p1r,a1p1i,b1p1r,b1p1i,cy_r1,bjmodn1,0x1);
		mixed_carry(a1p2r,a1p2i,b1p2r,b1p2i,cy_r2,bjmodn2,0x2);
		mixed_carry(a1p3r,a1p3i,b1p3r,b1p3i,cy_r3,bjmodn3,0x3);
		mixed_carry(a1p4r,a1p4i,b1p4r,b1p4i,cy_r4,bjmodn4,0x4);
		mixed_carry(a1p5r,a1p5i,b1p5r,b1p5i,cy_r5,bjmodn5,0x5);
		mixed_carry(a1p6r,a1p6i,b1p6r,b1p6i,cy_r6,bjmodn6,0x6);
		mixed_carry(a1p7r,a1p7i,b1p7r,b1p7i,cy_r7,bjmodn7,0x7);
		mixed_carry(a1p8r,a1p8i,b1p8r,b1p8i,cy_r8,bjmodn8,0x8);
		mixed_carry(a1p9r,a1p9i,b1p9r,b1p9i,cy_r9,bjmodn9,0x9);
		mixed_carry(a1pAr,a1pAi,b1pAr,b1pAi,cy_rA,bjmodnA,0xA);
		mixed_carry(a1pBr,a1pBi,b1pBr,b1pBi,cy_rB,bjmodnB,0xB);
		mixed_carry(a1pCr,a1pCi,b1pCr,b1pCi,cy_rC,bjmodnC,0xC);
		mixed_carry(a1pDr,a1pDi,b1pDr,b1pDi,cy_rD,bjmodnD,0xD);
		mixed_carry(a1pEr,a1pEi,b1pEr,b1pEi,cy_rE,bjmodnE,0xE);
		mixed_carry(a1pFr,a1pFi,b1pFr,b1pFi,cy_rF,bjmodnF,0xF);
*/
if(!j) {
	printf("J = 0, carry-step OUTputs:\n");
	printf("a1p0r,a1p0i, b1p0r,b1p0i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p0r,a1p0i, b1p0r,b1p0i);
	printf("a1p1r,a1p1i, b1p1r,b1p1i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p1r,a1p1i, b1p1r,b1p1i);
	printf("a1p2r,a1p2i, b1p2r,b1p2i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p2r,a1p2i, b1p2r,b1p2i);
	printf("a1p3r,a1p3i, b1p3r,b1p3i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p3r,a1p3i, b1p3r,b1p3i);
	printf("a1p4r,a1p4i, b1p4r,b1p4i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p4r,a1p4i, b1p4r,b1p4i);
	printf("a1p5r,a1p5i, b1p5r,b1p5i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p5r,a1p5i, b1p5r,b1p5i);
	printf("a1p6r,a1p6i, b1p6r,b1p6i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p6r,a1p6i, b1p6r,b1p6i);
	printf("a1p7r,a1p7i, b1p7r,b1p7i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p7r,a1p7i, b1p7r,b1p7i);
	printf("a1p8r,a1p8i, b1p8r,b1p8i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p8r,a1p8i, b1p8r,b1p8i);
	printf("a1p9r,a1p9i, b1p9r,b1p9i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1p9r,a1p9i, b1p9r,b1p9i);
	printf("a1pAr,a1pAi, b1pAr,b1pAi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1pAr,a1pAi, b1pAr,b1pAi);
	printf("a1pBr,a1pBi, b1pBr,b1pBi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1pBr,a1pBi, b1pBr,b1pBi);
	printf("a1pCr,a1pCi, b1pCr,b1pCi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1pCr,a1pCi, b1pCr,b1pCi);
	printf("a1pDr,a1pDi, b1pDr,b1pDi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1pDr,a1pDi, b1pDr,b1pDi);
	printf("a1pEr,a1pEi, b1pEr,b1pEi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1pEr,a1pEi, b1pEr,b1pEi);
	printf("a1pFr,a1pFi, b1pFr,b1pFi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a1pFr,a1pFi, b1pFr,b1pFi);

	printf("\niter %2u [full-pass = %u]: a01 OUT: %20.10e, %20.10e, cy = %20.10e\n",iter,full_pass,a1p0r,a1p0i,cy_r0);
}

	#else

	  if(USE_SHORT_CY_CHAIN < USE_SHORT_CY_CHAIN_MAX) {	// LOACC with tunable DWT-weights chaining

		// Re-init weights every 4th macro invocation to keep errors under control:
		jt = j1; jp = j2;
		cmplx_carry_norm_pow2_errcheck0(a[jt   ],a[jp   ],cy_r0,bjmodn0,0x0,prp_mult);
		cmplx_carry_fast_pow2_errcheck (a[jt+p1],a[jp+p1],cy_r1,bjmodn1,0x1,prp_mult);
		cmplx_carry_fast_pow2_errcheck (a[jt+p2],a[jp+p2],cy_r2,bjmodn2,0x2,prp_mult);
		cmplx_carry_fast_pow2_errcheck (a[jt+p3],a[jp+p3],cy_r3,bjmodn3,0x3,prp_mult);
		jt = j1+p4; jp = j2+p4;
		cmplx_carry_norm_pow2_errcheck0(a[jt   ],a[jp   ],cy_r4,bjmodn4,0x4,prp_mult);
		cmplx_carry_fast_pow2_errcheck (a[jt+p1],a[jp+p1],cy_r5,bjmodn5,0x5,prp_mult);
		cmplx_carry_fast_pow2_errcheck (a[jt+p2],a[jp+p2],cy_r6,bjmodn6,0x6,prp_mult);
		cmplx_carry_fast_pow2_errcheck (a[jt+p3],a[jp+p3],cy_r7,bjmodn7,0x7,prp_mult);
		jt = j1+p8; jp = j2+p8;
		cmplx_carry_norm_pow2_errcheck0(a[jt   ],a[jp   ],cy_r8,bjmodn8,0x8,prp_mult);
		cmplx_carry_fast_pow2_errcheck (a[jt+p1],a[jp+p1],cy_r9,bjmodn9,0x9,prp_mult);
		cmplx_carry_fast_pow2_errcheck (a[jt+p2],a[jp+p2],cy_rA,bjmodnA,0xA,prp_mult);
		cmplx_carry_fast_pow2_errcheck (a[jt+p3],a[jp+p3],cy_rB,bjmodnB,0xB,prp_mult);
		jt = j1+p12; jp = j2+p12;
		cmplx_carry_norm_pow2_errcheck0(a[jt   ],a[jp   ],cy_rC,bjmodnC,0xC,prp_mult);
		cmplx_carry_fast_pow2_errcheck (a[jt+p1],a[jp+p1],cy_rD,bjmodnD,0xD,prp_mult);
		cmplx_carry_fast_pow2_errcheck (a[jt+p2],a[jp+p2],cy_rE,bjmodnE,0xE,prp_mult);
		cmplx_carry_fast_pow2_errcheck (a[jt+p3],a[jp+p3],cy_rF,bjmodnF,0xF,prp_mult);

	  } else {	// HiACC:

	/*...set0 is slightly different from others:	*/
		jt = j1; jp = j2;
	   cmplx_carry_norm_pow2_errcheck0(a[jt   ],a[jp   ],cy_r0,bjmodn0,0x0,prp_mult);
		cmplx_carry_norm_pow2_errcheck(a[jt+p1],a[jp+p1],cy_r1,bjmodn1,0x1,prp_mult);
		cmplx_carry_norm_pow2_errcheck(a[jt+p2],a[jp+p2],cy_r2,bjmodn2,0x2,prp_mult);
		cmplx_carry_norm_pow2_errcheck(a[jt+p3],a[jp+p3],cy_r3,bjmodn3,0x3,prp_mult);
		jt = j1+p4; jp = j2+p4;
		cmplx_carry_norm_pow2_errcheck(a[jt   ],a[jp   ],cy_r4,bjmodn4,0x4,prp_mult);
		cmplx_carry_norm_pow2_errcheck(a[jt+p1],a[jp+p1],cy_r5,bjmodn5,0x5,prp_mult);
		cmplx_carry_norm_pow2_errcheck(a[jt+p2],a[jp+p2],cy_r6,bjmodn6,0x6,prp_mult);
		cmplx_carry_norm_pow2_errcheck(a[jt+p3],a[jp+p3],cy_r7,bjmodn7,0x7,prp_mult);
		jt = j1+p8; jp = j2+p8;
		cmplx_carry_norm_pow2_errcheck(a[jt   ],a[jp   ],cy_r8,bjmodn8,0x8,prp_mult);
		cmplx_carry_norm_pow2_errcheck(a[jt+p1],a[jp+p1],cy_r9,bjmodn9,0x9,prp_mult);
		cmplx_carry_norm_pow2_errcheck(a[jt+p2],a[jp+p2],cy_rA,bjmodnA,0xA,prp_mult);
		cmplx_carry_norm_pow2_errcheck(a[jt+p3],a[jp+p3],cy_rB,bjmodnB,0xB,prp_mult);
		jt = j1+p12; jp = j2+p12;
		cmplx_carry_norm_pow2_errcheck(a[jt   ],a[jp   ],cy_rC,bjmodnC,0xC,prp_mult);
		cmplx_carry_norm_pow2_errcheck(a[jt+p1],a[jp+p1],cy_rD,bjmodnD,0xD,prp_mult);
		cmplx_carry_norm_pow2_errcheck(a[jt+p2],a[jp+p2],cy_rE,bjmodnE,0xE,prp_mult);
		cmplx_carry_norm_pow2_errcheck(a[jt+p3],a[jp+p3],cy_rF,bjmodnF,0xF,prp_mult);

	  }	// LOACC or HIACC?

	#endif

		i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/

		co2 = co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					(and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

	#endif	// USE_AVX?

	}
	else	/* MODULUS_TYPE_FERMAT */
	{
	#ifdef USE_SSE2
		add0 = &prp_mult;
	#endif
	#ifdef USE_AVX512

		// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.
		tmp = half_arr+2;	// In Fermat-mod mode, use first 2 vector slots above half_arr for base and 1/base
		VEC_DBL_INIT(tmp, scale);
		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		tmp = base_negacyclic_root;	tm2 = tmp+1;

		// Hi-accuracy version needs 2 copies of each base root, one for each invocation of the SSE2_fermat_carry_norm_pow2 carry macri:
		l = (j >> 1);	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
		VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
		VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
		VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
		VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
		VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
		VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
		VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);
		tmp += 2;	tm2 += 2;
		l += 1;	k1=(l & NRTM1);	k2=(l >> NRT_BITS);
		dtmp=rn0[k1].re;			wt_im=rn0[k1].im;
		rt  =rn1[k2].re;			it   =rn1[k2].im;
		wt_re =dtmp*rt-wt_im*it;	wt_im =dtmp*it+wt_im*rt;
		VEC_DBL_INIT(tmp+  0,wt_re);	VEC_DBL_INIT(tm2+  0,wt_im);
		VEC_DBL_INIT(tmp+ 16,wt_re);	VEC_DBL_INIT(tm2+ 16,wt_im);

		// AVX-custom 8-way carry macro - each contains 8 of the [RADIX] stride-n/RADIX-separated carries
		// (processed independently in parallel), and steps through sequential-data indicies j,j+2,j+4,j+6,j+8,j+10,j+12,j+14:

		// The starting value of the literal pointer offsets following 'tmp' in these macro calls = RADIX*4*sizeof(vec_dbl)
		// which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
		// precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x200-byte quartet of base roots
		// uses the same 0x80-byte complex up-multiplier, so the literal offsets advance (+0x200-0x80) = -0x180 bytes between macro calls:
		l = 0x800;
		// Each AVX carry macro call also processes 4 prefetches of main-array data
		addr = a + j1 + pfetch_dist;	// tmp+l = base_negacyclic_root + 2*RADIX; after each macro-invoc incr by 16 and decr by 14 (via l -= 0x380) for net incr += 2, thus covering the [rbase * (j*I*Pi)/(2*RADIX)] const-data:
		tmp = base_negacyclic_root+  0;	SSE2_fermat_carry_norm_pow2_errcheck_X8(r00,tmp,l,cy_r0,cy_i0,half_arr,sign_mask, addr,p1,p2,p3,p4, add0); l -= 0x380;
		addr = a + j1 + pfetch_dist + p8;
		tmp = base_negacyclic_root+ 16;	SSE2_fermat_carry_norm_pow2_errcheck_X8(r10,tmp,l,cy_r8,cy_i8,half_arr,sign_mask, addr,p1,p2,p3,p4, add0);

	#elif defined(USE_AVX)

		// For a description of the data movement in AVX mode, see radix28_ditN_cy_dif1.
		tmp = half_arr+2;	// In Fermat-mod mode, use first 2 vector slots above half_arr for base and 1/base
		VEC_DBL_INIT(tmp, scale);
		/* Get the needed Nth root of -1: */
		add1 = (double *)&rn0[0];
		add2 = (double *)&rn1[0];

		tmp = base_negacyclic_root;	tm2 = tmp+1;

		// Hi-accuracy version needs 4 copies of each base root, one for each invocation of the SSE2_fermat_carry_norm_pow2 carry macro:
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

		// Startvalue of lit-pointer offsets following 'tmp' in these macro calls = #VEC_DBL_INIT calls in above section (RADIX*2) * sizeof(vec_dbl)
		// which is the byte offset between the 'active' negacyclic weights [pointed to by base_negacyclic_root] and the
		// precomputed multipliers in the HIACC-wrapped section of the SIMD data initializations. Each 0x100-byte quartet of base roots
		// uses the same 0x40-byte complex up-multiplier, so the literal offsets advance (+0x100-0x40) = -0xc0 bytes between macro calls:
		l = 0x400;
		// Each AVX carry macro call also processes 4 prefetches of main-array data
		addr = a + j1 + pfetch_dist;	// tmp+l = base_negacyclic_root + 2*RADIX; after each macro-invoc incr by 8 and decr by 6 (via l -= 0xc0) for net incr += 2, thus covering the [rbase * (j*I*Pi)/(2*RADIX)] const-data:
		tmp = base_negacyclic_root+  0;	SSE2_fermat_carry_norm_pow2_errcheck_X4(r00,tmp,l,cy_r0,cy_i0,half_arr,sign_mask, addr,p1,p2,p3, add0); l -= 0xc0;
		addr = a + j1 + pfetch_dist + p4 ;	// poff[] = p0,4,8,...
		tmp = base_negacyclic_root+  8;	SSE2_fermat_carry_norm_pow2_errcheck_X4(r08,tmp,l,cy_r4,cy_i4,half_arr,sign_mask, addr,p1,p2,p3, add0); l -= 0xc0;
		addr = a + j1 + pfetch_dist + p8 ;	// poff[] = p0,4,8,...
		tmp = base_negacyclic_root+ 16;	SSE2_fermat_carry_norm_pow2_errcheck_X4(r10,tmp,l,cy_r8,cy_i8,half_arr,sign_mask, addr,p1,p2,p3, add0); l -= 0xc0;
		addr = a + j1 + pfetch_dist + p12;	// poff[] = p0,4,8,...
		tmp = base_negacyclic_root+ 24;	SSE2_fermat_carry_norm_pow2_errcheck_X4(r18,tmp,l,cy_rC,cy_iC,half_arr,sign_mask, addr,p1,p2,p3, add0);

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

		// Each SSE2 carry macro call also processes 2 prefetches of main-array data
		addr = a + j1 + pfetch_dist;
		SSE2_fermat_carry_norm_pow2_errcheck_X2(r00,cy_r0,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2, addr,p1, add0);
		addr += p2;
		SSE2_fermat_carry_norm_pow2_errcheck_X2(r04,cy_r4,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2, addr,p1, add0);
		addr = a + j1 + pfetch_dist + p4;
		SSE2_fermat_carry_norm_pow2_errcheck_X2(r08,cy_r8,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2, addr,p1, add0);
		addr += p2;
		SSE2_fermat_carry_norm_pow2_errcheck_X2(r0C,cy_rC,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2, addr,p1, add0);
		addr = a + j1 + pfetch_dist + p8;
		SSE2_fermat_carry_norm_pow2_errcheck_X2(r10,cy_i0,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2, addr,p1, add0);
		addr += p2;
		SSE2_fermat_carry_norm_pow2_errcheck_X2(r14,cy_i4,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2, addr,p1, add0);
		addr = a + j1 + pfetch_dist + p12;
		SSE2_fermat_carry_norm_pow2_errcheck_X2(r18,cy_i8,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2, addr,p1, add0);
		addr += p2;
		SSE2_fermat_carry_norm_pow2_errcheck_X2(r1C,cy_iC,NRT_BITS,NRTM1,idx_offset,idx_incr,half_arr,sign_mask,add1,add2, addr,p1, add0);

	#else	// Scalar-double mode:

		ntmp = 0;	jt = j1; jp = j2;
		fermat_carry_norm_pow2_errcheck(a[jt   ],a[jp   ],cy_r0,cy_i0,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a[jt+p1],a[jp+p1],cy_r1,cy_i1,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a[jt+p2],a[jp+p2],cy_r2,cy_i2,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a[jt+p3],a[jp+p3],cy_r3,cy_i3,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR;
		jt = j1+p4; jp = j2+p4;
		fermat_carry_norm_pow2_errcheck(a[jt   ],a[jp   ],cy_r4,cy_i4,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a[jt+p1],a[jp+p1],cy_r5,cy_i5,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a[jt+p2],a[jp+p2],cy_r6,cy_i6,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a[jt+p3],a[jp+p3],cy_r7,cy_i7,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR;
		jt = j1+p8; jp = j2+p8;
		fermat_carry_norm_pow2_errcheck(a[jt   ],a[jp   ],cy_r8,cy_i8,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a[jt+p1],a[jp+p1],cy_r9,cy_i9,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a[jt+p2],a[jp+p2],cy_rA,cy_iA,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a[jt+p3],a[jp+p3],cy_rB,cy_iB,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR;
		jt = j1+p12; jp = j2+p12;
		fermat_carry_norm_pow2_errcheck(a[jt   ],a[jp   ],cy_rC,cy_iC,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a[jt+p1],a[jp+p1],cy_rD,cy_iD,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a[jt+p2],a[jp+p2],cy_rE,cy_iE,ntmp,NRTM1,NRT_BITS,prp_mult);	ntmp += NDIVR;
		fermat_carry_norm_pow2_errcheck(a[jt+p3],a[jp+p3],cy_rF,cy_iF,ntmp,NRTM1,NRT_BITS,prp_mult);

	#endif	/* #ifdef USE_SSE2 */

	}	/* if(MODULUS_TYPE == ...) */

/*...The radix-16 DIF pass is here:	*/

/* Four DIF radix-4 subconvolution, sans twiddles.	Cost each: 16 MOVapd, 20 ADD/SUBpd,  0 MULpd */

#ifdef USE_SSE2

	add0 = &a[j1    ];	// re-init this, because ptr used as a prefetch address in carry step above

	SSE2_RADIX16_DIF_NOTWIDDLE(add0,p1,p2,p3,p4,p8,r00,isrt2,cc0);

#else	/* !USE_SSE2 */

	#if USE_SCALAR_DFT_MACRO

		RADIX_16_DIF(a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
				//	a1p0r,a1p0i,a1p1r,a1p1i,a1p2r,a1p2i,a1p3r,a1p3i,a1p4r,a1p4i,a1p5r,a1p5i,a1p6r,a1p6i,a1p7r,a1p7i,a1p8r,a1p8i,a1p9r,a1p9i,a1pAr,a1pAi,a1pBr,a1pBi,a1pCr,a1pCi,a1pDr,a1pDi,a1pEr,a1pEi,a1pFr,a1pFi
					,a[j1    ],a[j2    ],a[j1+p1 ],a[j2+p1 ],a[j1+p2 ],a[j2+p2 ],a[j1+p3 ],a[j2+p3 ],a[j1+p4 ],a[j2+p4 ],a[j1+p5 ],a[j2+p5 ],a[j1+p6 ],a[j2+p6 ],a[j1+p7 ],a[j2+p7 ],a[j1+p8 ],a[j2+p8 ],a[j1+p9 ],a[j2+p9 ],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],a[j1+p15],a[j2+p15]
					,c,s)

	#else

	#ifdef USE_FGT61

	  #if PFETCH
		addr = &a[j1];
		prefetch_p_doubles(addr);
	  #endif
		/*...Block 1:	*/
		t3 =a1p0r -a1p8r;  	t1 =a1p0r +a1p8r;			m3 =b1p0r -b1p8r;  	m1 =b1p0r +b1p8r;
		t4 =a1p0i -a1p8i;	t2 =a1p0i +a1p8i;			m4 =b1p0i -b1p8i;	m$ =b1p0i +b1p8i;

		t7 =a1p4r -a1pCr;	t5 =a1p4r +a1pCr;			m7 =b1p4r -b1pCr;	m5 =b1p4r +b1pCr;
		t8 =a1p4i -a1pCi;	t6 =a1p4i +a1pCi;			m8 =b1p4i -b1pCi;	m6 =b1p4i +b1pCi;
	  #if PFETCH
		addp = addr+p1;
		prefetch_p_doubles(addp);
	  #endif
														rm =m5;			im =m6		;		// 0,2b
		rt =t5;	t5 =t1 -rt;			t1 =t1 +rt;			m5 =m1 -rm;		m6 =m$ -im	;		// m5,6: -2b,2b
		it =t6;	t6 =t2 -it;			t2 =t2 +it;			m1 =m1 +rm;		m$ =m$ +im	;		// m1,2: 0,4b
														rm =m7;			im =m8		;
		rt =t7;	t7 =t3 +t8;			t3 =t3 -t8;			m7 =m3 +im;		m8 =m4 -rm	;		// m3,4,7,8 all in -2b,2b
				t8 =t4 -rt;			t4 =t4 +rt;			m3 =m3 -im;		m4 =m4 +rm	;		// (Blocks 2-4 similar)
	  #if PFETCH
		addp = addr+p2;
		prefetch_p_doubles(addp);
	  #endif
		/*...Block 2:	*/
		t11=a1p2r -a1pAr;	t9 =a1p2r +a1pAr;			m11=b1p2r -b1pAr;	m9 =b1p2r +b1pAr;
		t12=a1p2i -a1pAi;	t10=a1p2i +a1pAi;			m12=b1p2i -b1pAi;	m10=b1p2i +b1pAi;

		t15=a1p6r -a1pEr;	t13=a1p6r +a1pEr;			m15=b1p6r -b1pEr;	m13=b1p6r +b1pEr;
		t16=a1p6i -a1pEi;	t14=a1p6i +a1pEi;			m16=b1p6i -b1pEi;	m14=b1p6i +b1pEi;
	  #if PFETCH
		addp = addr+p3;
		prefetch_p_doubles(addp);
	  #endif
														rm =m13;		im =m14		;
		rt =t13;	t13=t9 -rt;		t9 =t9 +rt;			m13=m9 -rm;		m14=m10-im	;
		it =t14;	t14=t10-it;		t10=t10+it;			m9 =m9 +rm;		m10=m10+im	;
														rm =m15;		im =m16		;
		rt =t15;	t15=t11+t16;	t11=t11-t16;		m15=m11+im;		m16=m12-rm	;
					t16=t12-rt;		t12=t12+rt;			m11=m11-im;		m12=m12+rm	;
	  #if PFETCH
		addp = addr+p4;
		prefetch_p_doubles(addp);
	  #endif
		/*...Block 3:	*/
		t19=a1p1r -a1p9r;  	t17=a1p1r +a1p9r;			m19=b1p1r -b1p9r;  	m17=b1p1r +b1p9r;
		t20=a1p1i -a1p9i;	t18=a1p1i +a1p9i;			m20=b1p1i -b1p9i;	m18=b1p1i +b1p9i;

		t23=a1p5r -a1pDr;	t21=a1p5r +a1pDr;			m23=b1p5r -b1pDr;	m21=b1p5r +b1pDr;
		t24=a1p5i -a1pDi;	t22=a1p5i +a1pDi;			m24=b1p5i -b1pDi;	m22=b1p5i +b1pDi;
	  #if PFETCH
		addp = addr+p5;
		prefetch_p_doubles(addp);
	  #endif
														rm =m21;					im =m22		;
		rt =t21;	t21=t17-rt;		t17=t17+rt;			m21=m17-rm;					m22=m18-im	;
		it =t22;	t22=t18-it;		t18=t18+it;			m17=m17+rm;					m18=m18+im	;
														rm =m23;					im =m24		;
		rt =t23;	t23=t19+t24;	t19=t19-t24;		m23=qreduce(m19+im+q4);		m24=qreduce(m20-rm+q4);	// all in -2b,2b
					t24=t20-rt;		t20=t20+rt;			m19=qreduce(m19-im+q4);		m20=qreduce(m20+rm+q4);	// prior to reduction
														// m19,20,23,24 are needed for CMUL, so reduce.
	  #if PFETCH
		addp = addr+p6;
		prefetch_p_doubles(addp);
	  #endif
		/*...Block 4:	*/
		t27=a1p3r -a1pBr;	t25=a1p3r +a1pBr;			m27=b1p3r -b1pBr;	m25=b1p3r +b1pBr;
		t28=a1p3i -a1pBi;	t26=a1p3i +a1pBi;			m28=b1p3i -b1pBi;	m26=b1p3i +b1pBi;

		t31=a1p7r -a1pFr;	t29=a1p7r +a1pFr;			m31=b1p7r -b1pFr;	m29=b1p7r +b1pFr;
		t32=a1p7i -a1pFi;	t30=a1p7i +a1pFi;			m32=b1p7i -b1pFi;	m30=b1p7i +b1pFi;
	  #if PFETCH
		addp = addr+p7;
		prefetch_p_doubles(addp);
	  #endif
														rm =m29;					im =m30		;
		rt =t29;	t29=t25-rt;		t25=t25+rt;			m29=m25-rm;					m30=m26-im	;
		it =t30;	t30=t26-it;		t26=t26+it;			m25=m25+rm;					m26=m26+im	;
														rm =m31;					im =m32		;
		rt =t31;	t31=t27+t32;	t27=t27-t32;		m31=qreduce(m27+im+q4);		m32=qreduce(m28-rm+q4);	// all in -2b,2b
					t32=t28-rt;		t28=t28+rt;			m27=qreduce(m27-im+q4);		m28=qreduce(m28+rm+q4);	// prior to reduction
														// m27,28,31,32 are needed for CMUL, so reduce.
	  #if PFETCH
		addp = addr+p8;
		prefetch_p_doubles(addp);
	  #endif

	/*
	!...and now do four more radix-4 transforms, including the internal twiddle factors:
	*/
													/*===============
													Input bounds on modular terms:
														-2b,2b: m3-8, 11-16, 21,22, 29,30
														  0,4b: m1,2, 9 ,10, 17,18, 25,26
													These are all reduced (in 0,b) because they are inputs to CMUL:
																m19,20,23,24,27,28,31,32
													===============*/
	/*...Block 1: t1,9,17,25 */
		jt = j1;		jp = j2;
		/* Debug: check for overflow of + terms: */	ASSERT(m1+m9 >= m1 && m$+m10 >= m$,"Overflow of [0,8b] term!");
		rt =t9;	t9 =t1 -rt;	t1 =t1 +rt;				rm =m9;	m9 =qreduce(m1 -rm+q4);	m1 =qreduce(m1 +rm   );	// +:   0,8b -> 0,b
		it =t10;t10=t2 -it;	t2 =t2 +it;				im =m10;m10=qreduce(m$ -im+q4);	m$ =qreduce(m$ +im   );	// -: -4b,4b -> 0,b

		rt =t25;t25=t17-rt;	t17=t17+rt;				rm =m25;m25=qreduce(m17-rm+q4);	m17=qreduce(m17+rm   );	// +:   0,8b -> 0,b
		it =t26;t26=t18-it;	t18=t18+it;				im =m26;m26=qreduce(m18-im+q4);	m18=qreduce(m18+im   );	// -: -4b,4b -> 0,b
													// + terms in 0,2b; - tems in -b,b:
		a[jt    ]= t1+t17;	a[jp    ]= t2+t18;		b[jt    ]=qreduce(m1 +m17   );	b[jp    ]=qreduce(m$ +m18   );
		a[jt+p1 ]= t1-t17;	a[jp+p1 ]= t2-t18;		b[jt+p1 ]=qreduce(m1 -m17+q2);	b[jp+p1 ]=qreduce(m$ -m18+q2);
		// mpy by E^4=i is inlined here:
		a[jt+p2 ]=t9 -t26;	a[jp+p2 ]=t10+t25;		b[jt+p2 ]=qreduce(m9 -m26+q2);	b[jp+p2 ]=qreduce(m10+m25+q2);
		a[jt+p3 ]=t9 +t26;	a[jp+p3 ]=t10-t25;		b[jt+p3 ]=qreduce(m9 +m26+q2);	b[jp+p3 ]=qreduce(m10-m25+q2);

	/*...Block 3: t5,13,21,29 */
		jt = j1 + p4;		jp = j2 + p4;			// All inputs in -2b,2b:
		rt =t13;t13=t5 +t14;t5 =t5 -t14;			rm =m13;m13=qreduce(m5 +m14+q4);m5 =qreduce(m5 -m14+q4);	// All 4 outputs in
				t14=t6 -rt;	t6 =t6 +rt;						m14=qreduce(m6 -rm +q4);m6 =qreduce(m6 +rm +q4);	// -4b,4b -> 0,b
	// twiddle mpy by E^2:
		rt =(t21-t22)*ISRT2;t22=(t21+t22)*ISRT2;			rm = mul_i2(m21-m22+q4);m22= mul_i2(m21+m22+q4);	// All 4 outputs in
t21=rt;	rt =(t30+t29)*ISRT2;it =(t30-t29)*ISRT2;	m21=rm;	rm = mul_i2(m30+m29+q4);im = mul_i2(m30-m29+q4);	// -4b,4b -> 0,b30
		t29=t21+rt;			t21=t21-rt;						m29=m21+rm;			m21=m21-rm;		// m21,22 in [-b30,b30]
		t30=t22+it;			t22=t22-it;						m30=m22+im;			m22=m22-im;		// m29,30 in [0,2*b30]

		a[jt    ]= t5+t21;	a[jp    ]= t6+t22;		b[jt    ]=qreduce( m5+m21+q4);	b[jp    ]=qreduce( m6+m22+q4);	// + and - in [0,b] + [-b30,b30] = [-b30,b+b30]
		a[jt+p1 ]= t5-t21;	a[jp+p1 ]= t6-t22;		b[jt+p1 ]=qreduce( m5-m21+q4);	b[jp+p1 ]=qreduce( m6-m22+q4);
		// mpy by E^4=i is inlined here:
		a[jt+p2 ]=t13-t30;	a[jp+p2 ]=t14+t29;		b[jt+p2 ]=qreduce(m13-m30+q3);	b[jp+p2 ]=qreduce(m14+m29   );	// + in [0,b] + [0,2*b30] = [0,b+2*b30)]
		a[jt+p3 ]=t13+t30;	a[jp+p3 ]=t14-t29;		b[jt+p3 ]=qreduce(m13+m30   );	b[jp+p3 ]=qreduce(m14-m29+q3);	// - in [0,b] - [0,2*b30] = [-2*b30,b]

	/*...Block 2: t3,11,19,27 */
		jt = j1 + p8;		jp = j2 + p8;
	// twiddle mpy by E^2							// m3,4,11,12 all in -2b,2b -> m11+-m12 in -4b,4b -> rm,im in 0,b30:
		rt =(t11-t12)*ISRT2;it =(t11+t12)*ISRT2;	rm =mul_i2(m11-m12+q4);	im =mul_i2(m11+m12+q4);
		t11=t3 -rt;			t3 =t3 +rt;				m11=m3 -rm;				m12=m4 -im;	// + in [-2b,2b] + [0,b30] = -2b,2b+b30
		t12=t4 -it;			t4 =t4 +it;				m3 =m3 +rm;				m4 =m4 +im;	// - in [-2b,2b] - [0,b30] = -2b-b30,2b

		rt =t19*c - t20*s;	t20=t20*c + t19*s;		cmul_modq8(m19,m20, cm,sm, &m19,&m20);	// 0,4q
t19=rt;	rt =t27*s - t28*c;	it =t28*s + t27*c;		cmul_modq8(m27,m28, sm,cm,  &rm, &im);	// 0,4q
		t27=t19-rt;			t19=t19+rt;				m27=qreduce(m19-rm+q4);	m19=qreduce(m19+rm);	// +:   0,8q ==> 0,b
		t28=t20-it;			t20=t20+it;				m28=qreduce(m20-im+q4);	m20=qreduce(m20+im);	// -: -4q,4q ==> 0,b
													// m3,4 in -2b,2b+b30;  m19,20 in 0,b:
		a[jt    ]= t3+t19;	a[jp    ]= t4+t20;		b[jt    ]=qreduce(m3+m19+q3);	b[jp    ]=qreduce(m4+m20+q3);	// + in [-2b,2b+b30] + [0,b] = -2b,3b+b30
		a[jt+p1 ]= t3-t19;	a[jp+p1 ]= t4-t20;		b[jt+p1 ]=qreduce(m3-m19+q4);	b[jp+p1 ]=qreduce(m4-m20+q4);	// - in [-2b,2b+b30] - [0,b] = -3b,2b+b30
		// mpy by E^4=i is inlined here:			// m11,12 in -2b-b30,2b;  m27,28 in 0,b:
		a[jt+p2 ]=t11-t28;	a[jp+p2 ]=t12+t27;		b[jt+p2 ]=qreduce(m11-m28+q5);	b[jp+p2 ]=qreduce(m12+m27+q4);	// + in [-2b-b30,2b] + [0,b] = -2b-b30,3b
		a[jt+p3 ]=t11+t28;	a[jp+p3 ]=t12-t27;		b[jt+p3 ]=qreduce(m11+m28+q4);	b[jp+p3 ]=qreduce(m12-m27+q5);	// - in [-2b-b30,2b] - [0,b] = -3b-b30,2b

	/*...Block 4: t7,15,23,31 */
		jt = j1 + p12;		jp = j2 + p12;
		/* twiddle mpy by -E^6 is here... */		// m7,8,15,16 all in -2b,2b -> m16+-m15 in -4b,4b -> rm,im in 0,b30:
		rt =(t16+t15)*ISRT2;it =(t16-t15)*ISRT2;	rm =mul_i2(m16+m15+q4);	im =mul_i2(m16-m15+q4);
		t15=t7 +rt;			t7 =t7 -rt;				m15=m7 +rm;				m16=m8 +im;	// + in [-2b,2b] + [0,b30] = -2b,2b+b30
		t16=t8 +it;			t8 =t8 -it;				m7 =m7 -rm;				m8 =m8 -im;	// - in [-2b,2b] - [0,b30] = -2b-b30,2b
		// E^3 = s,c; E^9 = -c,-s = -E^1, do latter via twiddle-mul by E^1, then flip signs on ensuing +- [re,im]:
		rt =t23*s - t24*c;	t24=t24*s + t23*c;		cmul_modq8(m23,m24, sm,cm, &m23,&m24);	// 0,4q
t23=rt;	rt =t31*c - t32*s;	it =t32*c + t31*s;		cmul_modq8(m31,m32, cm,sm,  &rm, &im);	// 0,4q
		t31=t23+rt;			t23=t23-rt;				m31=qreduce(m23+rm);	m23=qreduce(m23-rm+q4);	// +:   0,8q ==> 0,b
		t32=t24+it;			t24=t24-it;				m32=qreduce(m24+im);	m24=qreduce(m24-im+q4);	// -: -4q,4q ==> 0,b
													// m7,8 in -2b-b30,2b;  m23,24 in 0,b:
		a[jt    ]= t7+t23;	a[jp    ]= t8+t24;		b[jt    ]=qreduce( m7+m23+q4);	b[jp    ]=qreduce( m8+m24+q4);	// + in [-2b-b30,2b] + [0,b] = -2b-b30,3b
		a[jt+p1 ]= t7-t23;	a[jp+p1 ]= t8-t24;		b[jt+p1 ]=qreduce( m7-m23+q5);	b[jp+p1 ]=qreduce( m8-m24+q5);	// - in [-2b-b30,2b] - [0,b] = -3b-b30,2b
		// mpy by E^4=i is inlined here:			// m15,16 in -2b,2b+b30;  m31,32 in 0,b:
		a[jt+p2 ]=t15-t32;	a[jp+p2 ]=t16+t31;		b[jt+p2 ]=qreduce(m15-m32+q4);	b[jp+p2 ]=qreduce(m16+m31+q3);	// + in [-2b,2b+b30] + [0,b] = -2b,3b+b30
		a[jt+p3 ]=t15+t32;	a[jp+p3 ]=t16-t31;		b[jt+p3 ]=qreduce(m15+m32+q3);	b[jp+p3 ]=qreduce(m16-m31+q4);	// - in [-2b,2b+b30] - [0,b] = -3b,2b+b30

if(!j) {
	printf("J = 0, DIF1 OUTputs:\n");
	printf("a1p0r,a1p0i, b1p0r,b1p0i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j1    ],a[j1    +1], b[j1    ],b[j1    +1]);
	printf("a1p1r,a1p1i, b1p1r,b1p1i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j1+p1 ],a[j1+p1 +1], b[j1+p1 ],b[j1+p1 +1]);
	printf("a1p2r,a1p2i, b1p2r,b1p2i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j1+p2 ],a[j1+p2 +1], b[j1+p2 ],b[j1+p2 +1]);
	printf("a1p3r,a1p3i, b1p3r,b1p3i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j1+p3 ],a[j1+p3 +1], b[j1+p3 ],b[j1+p3 +1]);
	printf("a1p4r,a1p4i, b1p4r,b1p4i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j1+p4 ],a[j1+p4 +1], b[j1+p4 ],b[j1+p4 +1]);
	printf("a1p5r,a1p5i, b1p5r,b1p5i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j1+p5 ],a[j1+p5 +1], b[j1+p5 ],b[j1+p5 +1]);
	printf("a1p6r,a1p6i, b1p6r,b1p6i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j1+p6 ],a[j1+p6 +1], b[j1+p6 ],b[j1+p6 +1]);
	printf("a1p7r,a1p7i, b1p7r,b1p7i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j1+p7 ],a[j1+p7 +1], b[j1+p7 ],b[j1+p7 +1]);
	printf("a1p8r,a1p8i, b1p8r,b1p8i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j1+p8 ],a[j1+p8 +1], b[j1+p8 ],b[j1+p8 +1]);
	printf("a1p9r,a1p9i, b1p9r,b1p9i = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j1+p9 ],a[j1+p9 +1], b[j1+p9 ],b[j1+p9 +1]);
	printf("a1pAr,a1pAi, b1pAr,b1pAi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j1+p10],a[j1+p10+1], b[j1+p10],b[j1+p10+1]);
	printf("a1pBr,a1pBi, b1pBr,b1pBi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j1+p11],a[j1+p11+1], b[j1+p11],b[j1+p11+1]);
	printf("a1pCr,a1pCi, b1pCr,b1pCi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j1+p12],a[j1+p12+1], b[j1+p12],b[j1+p12+1]);
	printf("a1pDr,a1pDi, b1pDr,b1pDi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j1+p13],a[j1+p13+1], b[j1+p13],b[j1+p13+1]);
	printf("a1pEr,a1pEi, b1pEr,b1pEi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j1+p14],a[j1+p14+1], b[j1+p14],b[j1+p14+1]);
	printf("a1pFr,a1pFi, b1pFr,b1pFi = %20.10e, %20.10e, %20" PRIu64 ", %20" PRIu64 "\n",a[j1+p15],a[j1+p15+1], b[j1+p15],b[j1+p15+1]);
}
			/**********************************************/
	#else	// USE_FGT61 = False; Basic scalar-double mode:
			/**********************************************/
		#error 16-DFT outputs-into-scalar-temps no longer supported!
	// DIF1: Gather needed data (16 64-bit complex, i.e. 32 64-bit reals) and do first set of four length-4 transforms:
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

	#endif	// USE_FGT61 ?

  #endif	// USE_SCALAR_DFT_MACRO ?

#endif	/* USE_SSE2 */

	}	/* end for(j=_jstart[ithread]; j < _jhi[ithread]; j += 2) */

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart += nwt;
		jhi    += nwt;

		col += RADIX;
		co3 -= RADIX;
	}
}	/* end for(int k=1; k <= khi; k++) */

