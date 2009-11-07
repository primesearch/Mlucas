/*******************************************************************************
*                                                                              *
*   (C) 1997-2009 by Ernst W. Mayer.                                           *
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

/***************/

int radix13_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-13 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-14 complex DIF pass on the data in
!   the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n13,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodn10,bjmodn11,bjmodn12
		,i,j,j1,j2,jstart,jhi,iroot,root_incr,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
/*...Fast length-6 cyclic & acyclic convolution schemes need the following: */
	static double
		bp0  = -.08333333333333333333,	/* -1/12 */
		bp1  =  .30046260628866577442,
		bp20 =  .34691580671669062393,
		bp30 = -.09194269158406167643,
		bp21 = -.17731083280989152506,
		bp31 = -.33569135524251108482,
		bp2B =  .16960497390679909887,
		bp3B = -.24374866365844940839,

		bq00 = -.17413860115213590500,
		bq01 =  .57514072947400312138,
		bq0B =  .40100212832186721638,
		bq10 = -.46657980490972778969,
		bq11 = -.19455606720203287579,
		bq12 =  .64257503892782293874,
		bq13 =  .05328706921762039739,
		bq1A =  .17599523401809514904,
		bq1B = -.14126899798441247840;

	static double radix_inv, n2inv;
	double r1,i1,r2,i2,r3,i3,r4,i4
		,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25
		,aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r,aj1p6r,aj1p7r,aj1p8r,aj1p9r,aj1p10r,aj1p11r,aj1p12r
		,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i,aj1p6i,aj1p7i,aj1p8i,aj1p9i,aj1p10i,aj1p11i,aj1p12i
		,cy0,cy1,cy2,cy3,cy4,cy5,cy6,cy7,cy8,cy9,cy10,cy11,cy12,temp,frac,scale;
	double maxerr = 0.0;
#if PFETCH
	double *add0, *addr;
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;

/*...change n13 and n_div_wt to non-static to work around a gcc compiler bug. */
	n13   = n/13;
	n_div_nwt = n13 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n13)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/13 in radix13_ditN_cy_dif1.\n",iter);
		if(INTERACT)fprintf(stderr,"%s",cbuf);
		fp = fopen(   OFILE,"a");
		fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;
		err=ERR_CARRY;
		return(err);
	}

	if(p != psave)
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)13));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

/*   constant index offsets for array load/stores are here.	*/

		p1 = n13;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;
		p8 = p7 +p1;
		p9 = p8 +p1;
		p10= p9 +p1;
		p11= p10+p1;
		p12= p11+p1;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 = p9 + ( (p9 >> DAT_BITS) << PAD_BITS );
		p10= p10+ ( (p10>> DAT_BITS) << PAD_BITS );
		p11= p11+ ( (p11>> DAT_BITS) << PAD_BITS );
		p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );

		bjmodnini=0;
		for(j=0; j < n13; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
	}

/*...The radix-13 final DIT pass is here.	*/

	cy0 = 0;		/* init carries	*/
	cy1 = 0;
	cy2 = 0;
	cy3 = 0;
	cy4 = 0;
	cy5 = 0;
	cy6 = 0;
	cy7 = 0;
	cy8 = 0;
	cy9 = 0;
	cy10= 0;
	cy11= 0;
	cy12= 0;

	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		cy0 = -2;
	}
	else
	{
		ASSERT(HERE,0,"Radix-13 currently only supports LL test mode!");
	}

	*fracmax=0;	/* init max. fractional error	*/

	iroot = 0;	/* init sincos array index	*/
	root_incr = 1;	/* init sincos array index increment (set = 1 for normal carry pass, = 0 for wrapper pass)	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

	jstart = 0;
	jhi = jstart+nwt-1;
	khi = n_div_nwt;

for(outer=0; outer <= 1; outer++)
{
	i = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (i = 1).	*/

	bjmodn0 = 0;
	bjmodn1 = bjmodnini;
	bjmodn2 = bjmodn1 +bjmodnini-n; bjmodn2 = bjmodn2 + ( (-(int)((uint32)bjmodn2 >> 31)) & n);
	bjmodn3 = bjmodn2 +bjmodnini-n; bjmodn3 = bjmodn3 + ( (-(int)((uint32)bjmodn3 >> 31)) & n);
	bjmodn4 = bjmodn3 +bjmodnini-n; bjmodn4 = bjmodn4 + ( (-(int)((uint32)bjmodn4 >> 31)) & n);
	bjmodn5 = bjmodn4 +bjmodnini-n; bjmodn5 = bjmodn5 + ( (-(int)((uint32)bjmodn5 >> 31)) & n);
	bjmodn6 = bjmodn5 +bjmodnini-n; bjmodn6 = bjmodn6 + ( (-(int)((uint32)bjmodn6 >> 31)) & n);
	bjmodn7 = bjmodn6 +bjmodnini-n; bjmodn7 = bjmodn7 + ( (-(int)((uint32)bjmodn7 >> 31)) & n);
	bjmodn8 = bjmodn7 +bjmodnini-n; bjmodn8 = bjmodn8 + ( (-(int)((uint32)bjmodn8 >> 31)) & n);
	bjmodn9 = bjmodn8 +bjmodnini-n; bjmodn9 = bjmodn9 + ( (-(int)((uint32)bjmodn9 >> 31)) & n);
	bjmodn10= bjmodn9 +bjmodnini-n; bjmodn10= bjmodn10+ ( (-(int)((uint32)bjmodn10>> 31)) & n);
	bjmodn11= bjmodn10+bjmodnini-n; bjmodn11= bjmodn11+ ( (-(int)((uint32)bjmodn11>> 31)) & n);
	bjmodn12= bjmodn11+bjmodnini-n; bjmodn12= bjmodn12+ ( (-(int)((uint32)bjmodn12>> 31)) & n);

	col=0;
	co2=(n >> nwt_bits)-1+13;
	co3=co2-13;		/* At the start of each new j-loop, co3=co2-radix(1)	*/

	for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
	{
		for(j=jstart; j<jhi; j += 2)	/* Each inner loop execution processes (radix(1)*nwt) array data.	*/
		{
		#ifdef USE_SSE2
			j1 = (j & mask01) + br4[j&3];
			j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
		#else
			j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
		#endif
			j2 = j1+RE_IM_STRIDE;

			aj1p0r = a[j1  ];		aj1p0i = a[j2];

			t2  =a[j1+p1 ];			t3  =a[j2+p1 ];
			r1  =a[j1+p12];			i1  =a[j2+p12];
			t14 =t2  -r1;			t15 =t3  -i1;
			t2  =t2  +r1;			t3  =t3  +i1;

			t4  =a[j1+p2 ];			t5  =a[j2+p2 ];
			r1  =a[j1+p11];			i1  =a[j2+p11];
			t16 =t4  -r1;			t17 =t5  -i1;
			t4  =t4  +r1;			t5  =t5  +i1;

			t10 =a[j1+p3 ];			t11 =a[j2+p3 ];
			r1  =a[j1+p10];			i1  =a[j2+p10];
			t22 =t10 -r1;			t23 =t11 -i1;
			t10 =t10 +r1;			t11 =t11 +i1;

			t6  =a[j1+p4 ];			t7  =a[j2+p4 ];
			r1  =a[j1+p9 ];			i1  =a[j2+p9 ];
			t18 =t6  -r1;			t19 =t7  -i1;
			t6  =t6  +r1;			t7  =t7  +i1;

			t8  =a[j1+p5 ];			t9  =a[j2+p5 ];
			r1  =a[j1+p8 ];			i1  =a[j2+p8 ];
			t20 =r1  -t8;			t21 =i1  -t9;		/* flip signs on as3 */
			t8  =t8  +r1;			t9  =t9  +i1;

			t12 =a[j1+p6 ];			t13 =a[j2+p6 ];
			r1  =a[j1+p7 ];			i1  =a[j2+p7 ];
			t24 =t12 -r1;			t25 =t13 -i1;
			t12 =t12 +r1;			t13 =t13 +i1;

		/*** COSINE TERMS ***/

			t0  = t2 +t6 +t10;		t1  = t3 +t7 +t11;
			t2  = t2 -t6     ;		t3  = t3 -t7     ;
			t6  =     t6 -t10;		t7  =     t7 -t11;

			r1  = t4 +t8 +t12;		i1  = t5 +t9 +t13;
			t4  = t4     -t12;		t5  = t5     -t13;
			t8  =     t8 -t12;		t9  =     t9 -t13;

			t10 = t0+r1;			t11 = t1+i1;
			t0  = t0-r1;			t1  = t1-i1;
			t12 = t2-t8;			t13 = t3-t9;
			t2  = t2+t8;			t3  = t3+t9;
			t8  = t4+t6;			t9  = t5+t7;
			t6  = t4-t6;			t7  = t5-t7;

			/* polypro(A,B) mod P0 + x0 */
			r1 = t10*bp0 + aj1p0r;	i1 = t11*bp0 + aj1p0i;

			/* polypro(A,B) mod P1      */
			t0 = t0 *bp1;			t1 = t1 *bp1;

			aj1p0r += t10;			aj1p0i += t11;

			/* polypro(A,B) mod P2 = x^2-x+1: */
			t4  =(t12+t8 )*bp2B;	t5  =(t13+t9 )*bp2B;
			t12 = t12*bp20;			t13 = t13*bp20;
			t8  = t8 *bp21;			t9  = t9 *bp21;
			t8  = t12 - t8 ;		t9  = t13 - t9 ;
			t4  = t4  - t12;		t5  = t5  - t13;

			/* polypro(A,B) mod P3 = x^2+x+1: */
			t10 =(t2 -t6 )*bp3B;	t11 =(t3 -t7 )*bp3B;
			t2  = t2 *bp30;			t3  = t3 *bp30;
			t6  = t6 *bp31;			t7  = t7 *bp31;
			t6  = t2  - t6 ;		t7  = t3  - t7 ;
			t10 = t10 + t2 ;		t11 = t11 + t3 ;

			/* Combine the 4 modular polynomial product outputs according to CRT. */
			t2 = r1 +t0;			t3 = i1 +t1;	/* a */
			t0 = r1 -t0;			t1 = i1 -t1;	/* b */
			r1 = t8 +t6;			i1 = t9 +t7;	/* c */
			t8 = t8 -t6;			t9 = t9 -t7;	/* d */
			t6 = t4 +t10;			t7 = t5 +t11;	/* e */
			t4 = t4 -t10;			t5 = t5 -t11;	/* f */
	#if(LO_ADD)
			t2 = t2 - r1 + t4;		t3 = t3 -i1 + t5;	/* C3 = x^2 = a - c + f		2 add        */
			t12= t2 + 3*r1;			t13= t3 + 3*i1;		/* C1 = x^0 = x^2 + 3.c		1 add, 1 mul */
			t4 = t2 - 3*t4;			t5 = t3 - 3*t5;		/* C4 = x^4 = x^2 - 3.f		1 add, 1 mul */

			t10= t0 + t8 - t6;		t11= t1 + t9 - t7;	/* C2 = x^5 = b + d - e		2 add        */
			t6 = t10+ 3*t6;			t7 = t11+ 3*t7;		/* C6 = x^1 = x^5 + 3.e		1 add, 1 mul */
			t8 = t10- 3*t8;			t9 = t11- 3*t9;		/* C5 = x^3 = x^5 - 3.d		1 add, 1 mul */
	#else
			t10 = t2 -r1;			t11 = t3 -i1;
			t12 = t2 +r1 +r1 +t4;	t13 = t3 +i1 +i1 +t5;	/* C1 */
			t2  = t10        +t4;	t3  = t11        +t5;	/* C3 */
			t4  = t10    -t4 -t4;	t5  = t11    -t5 -t5;	/* C4 */

			r1  = t0 +t8;			i1  = t1 +t9;
			t10 = t0 -t8 -t8 -t6;	t11 = t1 -t9 -t9 -t7;	/* C5 */
			t8  = r1     +t6 +t6;	t9  = i1     +t7 +t7;	/* C6 */
			t6  = r1         -t6;	t7  = i1         -t7;	/* C2 */
	#endif
		/*** SINE TERMS ***/

			r1 = t14-t18+t22;		i1 = t15-t19+t23;
			t14= t14-t22;			t15= t15-t23;
			t18= t18+t22;			t19= t19+t23;

			t0 = t16-t20+t24;		t1 = t17-t21+t25;
			t16= t16-t24;			t17= t17-t25;
			t20= t20+t24;			t21= t21+t25;

			/* polypro(A,B) mod Q0 = x^2+1: */
			t22 =(r1 +t0 )*bq0B;	t23 =(i1 +t1 )*bq0B;
			r1  = r1 *bq00;			i1  = i1 *bq00;
			t0  = t0 *bq01;			t1  = t1 *bq01;
			t24 = r1  - t0;			t25 = i1  - t1;
			t22 = t22 - r1  - t0;	t23 = t23 - i1  - t1;	/* t0,1 free */

			/* polypro(A,B) mod Q1 = x^4-x^2+1: */
			r1  = t14+t18;			i1  = t15+t19;
			t0  = t14*bq10;			t1  = t15*bq10;
			r2  = t18*bq12;			i2  = t19*bq12;
			r2  = t0  - r2 ;		i2  = t1  - i2 ;
			t0  = r1 *bq1A - t0 ;	t1  = i1 *bq1A - t1 ;

			t14 = t14*bq11;			t15 = t15*bq11;
			t18 = t18*bq13;			t19 = t19*bq13;
			t18 = t14 - t18;		t19 = t15 - t19;
			t14 = r1 *bq1B - t14;	t15 = i1 *bq1B - t15;

			r1  = t16+t20;			i1  = t17+t21;
			r3  = t16*bq10;			i3  = t17*bq10;
			r4  = t20*bq12;			i4  = t21*bq12;
			r4  = r3  - r4;			i4  = i3  - i4;
			r3  = r1 *bq1A - r3 ;	i3  = i1 *bq1A - i3 ;

			t16 = t16*bq11;			t17 = t17*bq11;
			t20 = t20*bq13;			t21 = t21*bq13;
			t20 = t16 - t20;		t21 = t17 - t21;
			t16 = r1 *bq1B - t16;	t17 = i1 *bq1B - t17;

			t20 = t20+t0 +t16;		t21 = t21+t1 +t17;
			t16 = r2 -t16;			t17 = i2 -t17;
			t18 = t18+r4;			t19 = t19+i4;
			t14 = t14+r3 ;			t15 = t15+i3 ;

			/* Combine the 2 modular polynomial product outputs according to CRT. */
	#if(LO_ADD)
			t24 = t24 - t16 + t20;	t25 = t25 - t17 + t21;	/* S4 = x^4 = q00 - q10 + q12	2 add        */
			t16 = t24 + 3*t16;		t17 = t25 + 3*t17;	/* S1 = x^0 = x^4 + 3.q10	1 add, 1 mul */
			t20 = 3*t20 - t24;		t21 = 3*t21 - t25;	/* S3 = x^2 = 3.q12 - x^4	1 add, 1 mul */

			t22 = t22 - t18 + t14;	t23 = t23 - t19 + t15;	/* S2 = x^5 = q01 - q11 + q13	2 add        */
			t18 = t22 + 3*t18;		t19 = t23 + 3*t19;	/* S6 = x^1 = x^5 + 3.q11	1 add, 1 mul */
			t14 = 3*t14 - t22;		t15 = 3*t15 - t23;	/* S5 = x^3 = 3.q13 - x^5	1 add, 1 mul */

		/*...Inline multiply of sine parts by +-I into finishing phase, remembering that the acyclic
				 convo algorithm used above returns -s1 and -s5. 	*/

			aj1p1r =t12-t17;	aj1p1i =t13+t16;	/* C1 - S1 */
			aj1p2r =t10+t23;	aj1p2i =t11-t22;	/* C2 - S2 */
			aj1p3r =t2 +t21;	aj1p3i =t3 -t20;	/* C3 - S3 */
			aj1p4r =t4 +t25;	aj1p4i =t5 -t24;	/* C4 - S4 */
			aj1p5r =t8 -t15;	aj1p5i =t9 +t14;	/* C5 - S5 */
			aj1p6r =t6 +t19;	aj1p6i =t7 -t18;	/* C6 - S6 */
			aj1p7r =t6 -t19;	aj1p7i =t7 +t18;	/* C6 + S6 */
			aj1p8r =t8 +t15;	aj1p8i =t9 -t14;	/* C5 + S5 */
			aj1p9r =t4 -t25;	aj1p9i =t5 +t24;	/* C4 + S4 */
			aj1p10r=t2 -t21;	aj1p10i=t3 +t20;	/* C3 + S3 */
			aj1p11r=t10-t23;	aj1p11i=t11+t22;	/* C2 + S2 */
			aj1p12r=t12+t17;	aj1p12i=t13-t16;	/* C1 + S1 */
	#else
			r1  = t24+t16;		i1  = t25+t17;
			t24 = t24-t16;		t25 = t25-t17;
			t16 = t16+t20;		t17 = t17+t21;
			t0  = t22+t18;		t1  = t23+t19;
			t22 = t22-t18;		t23 = t23-t19;
			t18 = t18+t14;		t19 = t19+t15;

			t24 = t24+t20;		t25 = t25+t21;		/* S4 */
			t22 = t22+t14;		t23 = t23+t15;		/* S2 */
			t14 = r1 +t16;		t15 = i1 +t17;		/* S1 */
			t16 = t16+t16-r1 ;	t17 = t17+t17-i1 ;	/* S3 */
			t20 = t0 +t18;		t21 = t1 +t19;		/* S6 */
			t18 = t18+t18-t0 ;	t19 = t19+t19-t1 ;	/* S5 */

		/*...Inline multiply of sine parts by +-I into finishing phase, remembering that the acyclic
				 convo algorithm used above returns -s1 and -s5. 	*/

			aj1p1r =t12-t15;	aj1p1i =t13+t14;	/* C1 - S1 */
			aj1p2r =t6 +t23;	aj1p2i =t7 -t22;	/* C2 - S2 */
			aj1p3r =t2 +t17;	aj1p3i =t3 -t16;	/* C3 - S3 */
			aj1p4r =t4 +t25;	aj1p4i =t5 -t24;	/* C4 - S4 */
			aj1p5r =t10-t19;	aj1p5i =t11+t18;	/* C5 - S5 */
			aj1p6r =t8 +t21;	aj1p6i =t9 -t20;	/* C6 - S6 */
			aj1p7r =t8 -t21;	aj1p7i =t9 +t20;	/* C6 + S6 */
			aj1p8r =t10+t19;	aj1p8i =t11-t18;	/* C5 + S5 */
			aj1p9r =t4 -t25;	aj1p9i =t5 +t24;	/* C4 + S4 */
			aj1p10r=t2 -t17;	aj1p10i=t3 +t16;	/* C3 + S3 */
			aj1p11r=t6 -t23;	aj1p11i=t7 +t22;	/* C2 + S2 */
			aj1p12r=t12+t15;	aj1p12i=t13-t14;	/* C1 + S1 */
	#endif

/*..Now do the carries. Since the outputs would
    normally be getting dispatched to 13 separate blocks of the A-array, we need 13 separate carries.	*/

			l= j & (nwt-1);
			n_minus_sil   = n-si[l  ];
			n_minus_silp1 = n-si[l+1];
			sinwt   = si[nwt-l  ];
			sinwtm1 = si[nwt-l-1];

			wtl     =wt0[    l  ];
			wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
			wtlp1   =wt0[    l+1];
			wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

/*...set0 is slightly different from others:	*/
		 cmplx_carry_norm_errcheck0(aj1p0r ,aj1p0i ,cy0 ,bjmodn0 );
			cmplx_carry_norm_errcheck(aj1p1r ,aj1p1i ,cy1 ,bjmodn1 ,1 );
			cmplx_carry_norm_errcheck(aj1p2r ,aj1p2i ,cy2 ,bjmodn2 ,2 );
			cmplx_carry_norm_errcheck(aj1p3r ,aj1p3i ,cy3 ,bjmodn3 ,3 );
			cmplx_carry_norm_errcheck(aj1p4r ,aj1p4i ,cy4 ,bjmodn4 ,4 );
			cmplx_carry_norm_errcheck(aj1p5r ,aj1p5i ,cy5 ,bjmodn5 ,5 );
			cmplx_carry_norm_errcheck(aj1p6r ,aj1p6i ,cy6 ,bjmodn6 ,6 );
			cmplx_carry_norm_errcheck(aj1p7r ,aj1p7i ,cy7 ,bjmodn7 ,7 );
			cmplx_carry_norm_errcheck(aj1p8r ,aj1p8i ,cy8 ,bjmodn8 ,8 );
			cmplx_carry_norm_errcheck(aj1p9r ,aj1p9i ,cy9 ,bjmodn9 ,9 );
			cmplx_carry_norm_errcheck(aj1p10r,aj1p10i,cy10,bjmodn10,10);
			cmplx_carry_norm_errcheck(aj1p11r,aj1p11i,cy11,bjmodn11,11);
			cmplx_carry_norm_errcheck(aj1p12r,aj1p12i,cy12,bjmodn12,12);

			i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
				 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

/*...The radix-13 DIF pass is here:	*/
#if PFETCH
add0 = & a[j2];
prefetch_p_doubles(add0);
#endif
	/*    gather the needed data (13 64-bit complex, i.e. 26 64-bit reals) and get needed xi+-xj combos:	*/
			t2 =aj1p1r+aj1p12r;	t3 =aj1p1i+aj1p12i;
			t14=aj1p1r-aj1p12r;	t15=aj1p1i-aj1p12i;

			t4 =aj1p2r+aj1p11r;	t5 =aj1p2i+aj1p11i;
			t16=aj1p2r-aj1p11r;	t17=aj1p2i-aj1p11i;

			t10=aj1p3r+aj1p10r;	t11=aj1p3i+aj1p10i;
			t22=aj1p3r-aj1p10r;	t23=aj1p3i-aj1p10i;

			t6 =aj1p4r+aj1p9r ;	t7 =aj1p4i+aj1p9i ;
			t18=aj1p4r-aj1p9r ;	t19=aj1p4i-aj1p9i ;

			t8 =aj1p5r+aj1p8r ;	t9 =aj1p5i+aj1p8i ;
			t20=aj1p8r-aj1p5r ;	t21=aj1p8i-aj1p5i ;

			t12=aj1p6r+aj1p7r ;	t13=aj1p6i+aj1p7i ;
			t24=aj1p6r-aj1p7r ;	t25=aj1p6i-aj1p7i ;
#if PFETCH
addr = add0+p1;
prefetch_p_doubles(addr);
#endif
		/*** COSINE TERMS ***/

			t0  = t2 +t6 +t10;		t1  = t3 +t7 +t11;
			t2  = t2 -t6     ;		t3  = t3 -t7     ;
			t6  =     t6 -t10;		t7  =     t7 -t11;

			r1  = t4 +t8 +t12;		i1  = t5 +t9 +t13;
			t4  = t4     -t12;		t5  = t5     -t13;
			t8  =     t8 -t12;		t9  =     t9 -t13;

			t10 = t0+r1;			t11 = t1+i1;
			t0  = t0-r1;			t1  = t1-i1;
			t12 = t2-t8;			t13 = t3-t9;
			t2  = t2+t8;			t3  = t3+t9;
			t8  = t4+t6;			t9  = t5+t7;
			t6  = t4-t6;			t7  = t5-t7;
#if PFETCH
addr = add0+p2;
prefetch_p_doubles(addr);
#endif
			/* polypro(A,B) mod P0 + x0 */
			r1 = t10*bp0 + aj1p0r;		i1 = t11*bp0 + aj1p0i;

			/* polypro(A,B) mod P1      */
			t0 = t0 *bp1;			t1 = t1 *bp1;

			a[j1  ] = aj1p0r+t10;		a[j2] = aj1p0i+t11;

			/* polypro(A,B) mod P2 = x^2-x+1: */
			t4  =(t12+t8 )*bp2B;	t5  =(t13+t9 )*bp2B;
			t12 = t12*bp20;		t13 = t13*bp20;
			t8  = t8 *bp21;		t9  = t9 *bp21;
			t8  = t12 - t8 ;	t9  = t13 - t9 ;
			t4  = t4  - t12;	t5  = t5  - t13;
#if PFETCH
addr = add0+p3;
prefetch_p_doubles(addr);
#endif
			/* polypro(A,B) mod P3 = x^2+x+1: */
			t10 =(t2 -t6 )*bp3B;	t11 =(t3 -t7 )*bp3B;
			t2  = t2 *bp30;		t3  = t3 *bp30;
			t6  = t6 *bp31;		t7  = t7 *bp31;
			t6  = t2  - t6 ;	t7  = t3  - t7 ;
			t10 = t10 + t2 ;	t11 = t11 + t3 ;
#if PFETCH
addr = add0+p4;
prefetch_p_doubles(addr);
#endif
			/* Combine the 4 modular polynomial product outputs according to CRT. */
			t2 = r1 +t0;		t3 = i1 +t1;
			t0 = r1 -t0;		t1 = i1 -t1;
			r1 = t8 +t6;		i1 = t9 +t7;
			t8 = t8 -t6;		t9 = t9 -t7;
			t6 = t4 +t10;		t7 = t5 +t11;
			t4 = t4 -t10;		t5 = t5 -t11;
#if PFETCH
addr = add0+p5;
prefetch_p_doubles(addr);
#endif
	#if(LO_ADD)
			t2 = t2 - r1 + t4;	t3 = t3 -i1 + t5;	/* C3 = x^2 = a - c + f		2 add        */
			t12= t2 + 3*r1;		t13= t3 + 3*i1;		/* C1 = x^0 = x^2 + 3.c		1 add, 1 mul */
			t4 = t2 - 3*t4;		t5 = t3 - 3*t5;		/* C4 = x^4 = x^2 - 3.f		1 add, 1 mul */

			t10= t0 + t8 - t6;	t11= t1 + t9 - t7;	/* C2 = x^5 = b + d - e		2 add        */
			t6 = t10+ 3*t6;		t7 = t11+ 3*t7;		/* C6 = x^1 = x^5 + 3.e		1 add, 1 mul */
			t8 = t10- 3*t8;		t9 = t11- 3*t9;		/* C5 = x^3 = x^5 - 3.d		1 add, 1 mul */
	#else
			t10 = t2 -r1;		t11 = t3 -i1;
			t12 = t2 +r1 +r1 +t4;	t13 = t3 +i1 +i1 +t5;	/* C1 */
			t2  = t10        +t4;	t3  = t11        +t5;	/* C3 */
			t4  = t10    -t4 -t4;	t5  = t11    -t5 -t5;	/* C4 */

			r1  = t0 +t8;		i1  = t1 +t9;
			t10 = t0 -t8 -t8 -t6;	t11 = t1 -t9 -t9 -t7;	/* C5 */
			t8  = r1     +t6 +t6;	t9  = i1     +t7 +t7;	/* C6 */
			t6  = r1         -t6;	t7  = i1         -t7;	/* C2 */
	#endif
#if PFETCH
addr = add0+p6;
prefetch_p_doubles(addr);
#endif

		/*** SINE TERMS ***/

			r1 = t14-t18+t22;	i1 = t15-t19+t23;
			t14= t14-t22;		t15= t15-t23;
			t18= t18+t22;		t19= t19+t23;

			t0 = t16-t20+t24;	t1 = t17-t21+t25;
			t16= t16-t24;		t17= t17-t25;
			t20= t20+t24;		t21= t21+t25;
#if PFETCH
addr = add0+p7;
prefetch_p_doubles(addr);
#endif
			/* polypro(A,B) mod Q0 = x^2+1: */
			t22 =(r1 +t0 )*bq0B;	t23 =(i1 +t1 )*bq0B;
			r1  = r1 *bq00;		i1  = i1 *bq00;
			t0  = t0 *bq01;		t1  = t1 *bq01;
			t24 = r1  - t0;		t25 = i1  - t1;
			t22 = t22 - r1  - t0;	t23 = t23 - i1  - t1;	/* t0,1 free */

			/* polypro(A,B) mod Q1 = x^4-x^2+1: */
			r1  = t14+t18;		i1  = t15+t19;
			t0  = t14*bq10;		t1  = t15*bq10;
			r2  = t18*bq12;		i2  = t19*bq12;
			r2  = t0  - r2 ;	i2  = t1  - i2 ;
			t0  = r1 *bq1A - t0 ;	t1  = i1 *bq1A - t1 ;
#if PFETCH
addr = add0+p8;
prefetch_p_doubles(addr);
#endif
			t14 = t14*bq11;		t15 = t15*bq11;
			t18 = t18*bq13;		t19 = t19*bq13;
			t18 = t14 - t18;	t19 = t15 - t19;
			t14 = r1 *bq1B - t14;	t15 = i1 *bq1B - t15;

			r1  = t16+t20;		i1  = t17+t21;
			r3  = t16*bq10;		i3  = t17*bq10;
			r4  = t20*bq12;		i4  = t21*bq12;
			r4  = r3  - r4;		i4  = i3  - i4;
			r3  = r1 *bq1A - r3 ;	i3  = i1 *bq1A - i3 ;
#if PFETCH
addr = add0+p9;
prefetch_p_doubles(addr);
#endif
			t16 = t16*bq11;		t17 = t17*bq11;
			t20 = t20*bq13;		t21 = t21*bq13;
			t20 = t16 - t20;	t21 = t17 - t21;
			t16 = r1 *bq1B - t16;	t17 = i1 *bq1B - t17;

			t20 = t20+t0 +t16;	t21 = t21+t1 +t17;
			t16 = r2 -t16;		t17 = i2 -t17;
			t18 = t18+r4;		t19 = t19+i4;
			t14 = t14+r3 ;		t15 = t15+i3 ;
#if PFETCH
addr = add0+p10;
prefetch_p_doubles(addr);
#endif
		/* Combine the 2 modular polynomial product outputs according to CRT. */
	#if(LO_ADD)
			t24 = t24 - t16 + t20;	t25 = t25 - t17 + t21;	/* S4 = x^4 = q00 - q10 + q12	2 add        */
			t16 = t24 + 3*t16;	t17 = t25 + 3*t17;	/* S1 = x^0 = x^4 + 3.q10	1 add, 1 mul */
			t20 = 3*t20 - t24;	t21 = 3*t21 - t25;	/* S3 = x^2 = 3.q12 - x^4	1 add, 1 mul */

			t22 = t22 - t18 + t14;	t23 = t23 - t19 + t15;	/* S2 = x^5 = q01 - q11 + q13	2 add        */
			t18 = t22 + 3*t18;	t19 = t23 + 3*t19;	/* S6 = x^1 = x^5 + 3.q11	1 add, 1 mul */
			t14 = 3*t14 - t22;	t15 = 3*t15 - t23;	/* S5 = x^3 = 3.q13 - x^5	1 add, 1 mul */

		/*...Inline multiply of sine parts by +-I into finishing phase, remembering that the acyclic
				 convo algorithm used above returns -s1 and -s5. 	*/
		#if PFETCH
		addr = add0+p11;
		prefetch_p_doubles(addr);
		#endif
			a[j1+p1 ]=t12+t17;	a[j2+p1 ]=t13-t16;	/* C1 - S1 */
			a[j1+p2 ]=t10-t23;	a[j2+p2 ]=t11+t22;	/* C2 - S2 */
			a[j1+p3 ]=t2 -t21;	a[j2+p3 ]=t3 +t20;	/* C3 - S3 */
			a[j1+p4 ]=t4 -t25;	a[j2+p4 ]=t5 +t24;	/* C4 - S4 */
			a[j1+p5 ]=t8 +t15;	a[j2+p5 ]=t9 -t14;	/* C5 - S5 */
			a[j1+p6 ]=t6 -t19;	a[j2+p6 ]=t7 +t18;	/* C6 - S6 */
			a[j1+p7 ]=t6 +t19;	a[j2+p7 ]=t7 -t18;	/* C6 + S6 */
			a[j1+p8 ]=t8 -t15;	a[j2+p8 ]=t9 +t14;	/* C5 + S5 */
			a[j1+p9 ]=t4 +t25;	a[j2+p9 ]=t5 -t24;	/* C4 + S4 */
			a[j1+p10]=t2 +t21;	a[j2+p10]=t3 -t20;	/* C3 + S3 */
			a[j1+p11]=t10+t23;	a[j2+p11]=t11-t22;	/* C2 + S2 */
			a[j1+p12]=t12-t17;	a[j2+p12]=t13+t16;	/* C1 + S1 */
	#else
			r1  = t24+t16;		i1  = t25+t17;
			t24 = t24-t16;		t25 = t25-t17;
			t16 = t16+t20;		t17 = t17+t21;
			t0  = t22+t18;		t1  = t23+t19;
			t22 = t22-t18;		t23 = t23-t19;
			t18 = t18+t14;		t19 = t19+t15;

			t24 = t24+t20;		t25 = t25+t21;		/* S4 */
			t22 = t22+t14;		t23 = t23+t15;		/* S2 */
			t14 = r1 +t16;		t15 = i1 +t17;		/* S1 */
			t16 = t16+t16-r1 ;	t17 = t17+t17-i1 ;	/* S3 */
			t20 = t0 +t18;		t21 = t1 +t19;		/* S6 */
			t18 = t18+t18-t0 ;	t19 = t19+t19-t1 ;	/* S5 */

		/*...Inline multiply of sine parts by +-I into finishing phase, remembering that the acyclic
				 convo algorithm used above returns -s1 and -s5. 	*/
		#if PFETCH
		addr = add0+p11;
		prefetch_p_doubles(addr);
		#endif
			a[j1+p1 ]=t12+t15;	a[j2+p1 ]=t13-t14;	/* C1 - S1 */
			a[j1+p2 ]=t6 -t23;	a[j2+p2 ]=t7 +t22;	/* C2 - S2 */
			a[j1+p3 ]=t2 -t17;	a[j2+p3 ]=t3 +t16;	/* C3 - S3 */
			a[j1+p4 ]=t4 -t25;	a[j2+p4 ]=t5 +t24;	/* C4 - S4 */
			a[j1+p5 ]=t10+t19;	a[j2+p5 ]=t11-t18;	/* C5 - S5 */
			a[j1+p6 ]=t8 -t21;	a[j2+p6 ]=t9 +t20;	/* C6 - S6 */
			a[j1+p7 ]=t8 +t21;	a[j2+p7 ]=t9 -t20;	/* C6 + S6 */
			a[j1+p8 ]=t10-t19;	a[j2+p8 ]=t11+t18;	/* C5 + S5 */
			a[j1+p9 ]=t4 +t25;	a[j2+p9 ]=t5 -t24;	/* C4 + S4 */
			a[j1+p10]=t2 +t17;	a[j2+p10]=t3 -t16;	/* C3 + S3 */
			a[j1+p11]=t6 +t23;	a[j2+p11]=t7 -t22;	/* C2 + S2 */
			a[j1+p12]=t12-t15;	a[j2+p12]=t13+t14;	/* C1 + S1 */
	#endif
#if PFETCH
addr = add0+p12;
prefetch_p_doubles(addr);
#endif
		iroot += root_incr;		/* increment sincos index.	*/

		}

		jstart += nwt;
		jhi    += nwt;
		col += 13;
		co3 -= 13;

	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-13 forward DIF FFT of the first block of 13 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 13 outputs of (1);
!   (3) Reweight and perform a radix-13 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 13 elements and repeat (1-4).
*/

	t1  = cy12;
	cy12= cy11;
	cy11= cy10;
	cy10= cy9;
	cy9 = cy8;
	cy8 = cy7;
	cy7 = cy6;
	cy6 = cy5;
	cy5 = cy4;
	cy4 = cy3;
	cy3 = cy2;
	cy2 = cy1;
	cy1 = cy0;
	cy0 = t1;

	iroot = 0;
	root_incr = 0;
	scale = 1;

	jstart = 0;
	jhi = 7;
	khi = 1;

	for(j=0; j<=jhi; j++)
	{
		a[j    ] *= radix_inv;
		a[j+p1 ] *= radix_inv;
		a[j+p2 ] *= radix_inv;
		a[j+p3 ] *= radix_inv;
		a[j+p4 ] *= radix_inv;
		a[j+p5 ] *= radix_inv;
		a[j+p6 ] *= radix_inv;
		a[j+p7 ] *= radix_inv;
		a[j+p8 ] *= radix_inv;
		a[j+p9 ] *= radix_inv;
		a[j+p10] *= radix_inv;
		a[j+p11] *= radix_inv;
		a[j+p12] *= radix_inv;
	}
}

	if(fabs(cy0)+fabs(cy1)+fabs(cy2)+fabs(cy3)+fabs(cy4)+fabs(cy5)+fabs(cy6)+fabs(cy7)+fabs(cy8)+fabs(cy9)+fabs(cy10)+fabs(cy11)+fabs(cy12) != 0.0)
	{
			sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix13_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
			if(INTERACT)fprintf(stderr,"%s",cbuf);
			fp = fopen(   OFILE,"a");
			fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
			fclose(fp);	fp = 0x0;
			fclose(fq);	fq = 0x0;
			err=ERR_CARRY;
			return(err);
	}
	*fracmax = maxerr;
	return(0);
}


/***************/

int radix13_ditN_cy_dif1_nochk(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter                 , uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-13 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-14 complex DIF pass on the data in
!   the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n13,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodn10,bjmodn11,bjmodn12
		,i,j,j1,j2,jstart,jhi,iroot,root_incr,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
/*...Fast length-6 cyclic & acyclic convolution schemes need the following: */
	static double
		bp0  = -.08333333333333333333,	/* -1/12 */
		bp1  =  .30046260628866577442,
		bp20 =  .34691580671669062393,
		bp30 = -.09194269158406167643,
		bp21 = -.17731083280989152506,
		bp31 = -.33569135524251108482,
		bp2B =  .16960497390679909887,
		bp3B = -.24374866365844940839,

		bq00 = -.17413860115213590500,
		bq01 =  .57514072947400312138,
		bq0B =  .40100212832186721638,
		bq10 = -.46657980490972778969,
		bq11 = -.19455606720203287579,
		bq12 =  .64257503892782293874,
		bq13 =  .05328706921762039739,
		bq1A =  .17599523401809514904,
		bq1B = -.14126899798441247840;

	static double radix_inv, n2inv;
	double r1,i1,r2,i2,r3,i3,r4,i4
		,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25
		,aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r,aj1p6r,aj1p7r,aj1p8r,aj1p9r,aj1p10r,aj1p11r,aj1p12r
		,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i,aj1p6i,aj1p7i,aj1p8i,aj1p9i,aj1p10i,aj1p11i,aj1p12i
		,cy0,cy1,cy2,cy3,cy4,cy5,cy6,cy7,cy8,cy9,cy10,cy11,cy12,temp,scale;
#if PFETCH
	double *add0, *addr;
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;

/*...change n13 and n_div_wt to non-static to work around a gcc compiler bug. */
	n13   = n/13;
	n_div_nwt = n13 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n13)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/13 in radix13_ditN_cy_dif1.\n",iter);
		if(INTERACT)fprintf(stderr,"%s",cbuf);
		fp = fopen(   OFILE,"a");
		fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
		fclose(fp);	fp = 0x0;
		fclose(fq);	fq = 0x0;
		err=ERR_CARRY;
		return(err);
	}

	if(p != psave)
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		psave = p;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)13));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

/*   constant index offsets for array load/stores are here.	*/

		p1 = n13;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;
		p8 = p7 +p1;
		p9 = p8 +p1;
		p10= p9 +p1;
		p11= p10+p1;
		p12= p11+p1;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 = p9 + ( (p9 >> DAT_BITS) << PAD_BITS );
		p10= p10+ ( (p10>> DAT_BITS) << PAD_BITS );
		p11= p11+ ( (p11>> DAT_BITS) << PAD_BITS );
		p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );

		bjmodnini=0;
		for(j=0; j < n13; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
	}

/*...The radix-13 final DIT pass is here.	*/

	cy0 = 0;		/* init carries	*/
	cy1 = 0;
	cy2 = 0;
	cy3 = 0;
	cy4 = 0;
	cy5 = 0;
	cy6 = 0;
	cy7 = 0;
	cy8 = 0;
	cy9 = 0;
	cy10= 0;
	cy11= 0;
	cy12= 0;

	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		cy0 = -2;
	}
	else
	{
		ASSERT(HERE,0,"Radix-13 currently only supports LL test mode!");
	}

	iroot = 0;	/* init sincos array index	*/
	root_incr = 1;	/* init sincos array index increment (set = 1 for normal carry pass, = 0 for wrapper pass)	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

	jstart = 0;
	jhi = jstart+nwt-1;
	khi = n_div_nwt;

for(outer=0; outer <= 1; outer++)
{
	i = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (i = 1).	*/

	bjmodn0 = 0;
	bjmodn1 = bjmodnini;
	bjmodn2 = bjmodn1 +bjmodnini-n; bjmodn2 = bjmodn2 + ( (-(int)((uint32)bjmodn2 >> 31)) & n);
	bjmodn3 = bjmodn2 +bjmodnini-n; bjmodn3 = bjmodn3 + ( (-(int)((uint32)bjmodn3 >> 31)) & n);
	bjmodn4 = bjmodn3 +bjmodnini-n; bjmodn4 = bjmodn4 + ( (-(int)((uint32)bjmodn4 >> 31)) & n);
	bjmodn5 = bjmodn4 +bjmodnini-n; bjmodn5 = bjmodn5 + ( (-(int)((uint32)bjmodn5 >> 31)) & n);
	bjmodn6 = bjmodn5 +bjmodnini-n; bjmodn6 = bjmodn6 + ( (-(int)((uint32)bjmodn6 >> 31)) & n);
	bjmodn7 = bjmodn6 +bjmodnini-n; bjmodn7 = bjmodn7 + ( (-(int)((uint32)bjmodn7 >> 31)) & n);
	bjmodn8 = bjmodn7 +bjmodnini-n; bjmodn8 = bjmodn8 + ( (-(int)((uint32)bjmodn8 >> 31)) & n);
	bjmodn9 = bjmodn8 +bjmodnini-n; bjmodn9 = bjmodn9 + ( (-(int)((uint32)bjmodn9 >> 31)) & n);
	bjmodn10= bjmodn9 +bjmodnini-n; bjmodn10= bjmodn10+ ( (-(int)((uint32)bjmodn10>> 31)) & n);
	bjmodn11= bjmodn10+bjmodnini-n; bjmodn11= bjmodn11+ ( (-(int)((uint32)bjmodn11>> 31)) & n);
	bjmodn12= bjmodn11+bjmodnini-n; bjmodn12= bjmodn12+ ( (-(int)((uint32)bjmodn12>> 31)) & n);

	col=0;
	co2=(n >> nwt_bits)-1+13;
	co3=co2-13;		/* At the start of each new j-loop, co3=co2-radix(1)	*/

	for(k=1; k <= khi; k++)	/* Do n/(radix(1)*nwt) outer loop executions...	*/
	{
		for(j=jstart; j<jhi; j += 2)	/* Each inner loop execution processes (radix(1)*nwt) array data.	*/
		{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
			j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
			j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

			aj1p0r = a[j1  ];	aj1p0i = a[j2];

		t2  =a[j1+p1 ];	t3  =a[j2+p1 ];
		r1  =a[j1+p12];	i1  =a[j2+p12];
		t14 =t2  -r1;	t15 =t3  -i1;
		t2  =t2  +r1;	t3  =t3  +i1;

		t4  =a[j1+p2 ];	t5  =a[j2+p2 ];
		r1  =a[j1+p11];	i1  =a[j2+p11];
		t16 =t4  -r1;	t17 =t5  -i1;
		t4  =t4  +r1;	t5  =t5  +i1;

		t10 =a[j1+p3 ];	t11 =a[j2+p3 ];
		r1  =a[j1+p10];	i1  =a[j2+p10];
		t22 =t10 -r1;	t23 =t11 -i1;
		t10 =t10 +r1;	t11 =t11 +i1;

		t6  =a[j1+p4 ];	t7  =a[j2+p4 ];
		r1  =a[j1+p9 ];	i1  =a[j2+p9 ];
		t18 =t6  -r1;	t19 =t7  -i1;
		t6  =t6  +r1;	t7  =t7  +i1;

		t8  =a[j1+p5 ];	t9  =a[j2+p5 ];
		r1  =a[j1+p8 ];	i1  =a[j2+p8 ];
		t20 =r1  -t8;	t21 =i1  -t9;		/* flip signs on as3 */
		t8  =t8  +r1;	t9  =t9  +i1;

		t12 =a[j1+p6 ];	t13 =a[j2+p6 ];
		r1  =a[j1+p7 ];	i1  =a[j2+p7 ];
		t24 =t12 -r1;	t25 =t13 -i1;
		t12 =t12 +r1;	t13 =t13 +i1;

	/*** COSINE TERMS ***/

		t0  = t2 +t6 +t10;		t1  = t3 +t7 +t11;
		t2  = t2 -t6     ;		t3  = t3 -t7     ;
		t6  =     t6 -t10;		t7  =     t7 -t11;

		r1  = t4 +t8 +t12;		i1  = t5 +t9 +t13;
		t4  = t4     -t12;		t5  = t5     -t13;
		t8  =     t8 -t12;		t9  =     t9 -t13;

		t10 = t0+r1;			t11 = t1+i1;
		t0  = t0-r1;			t1  = t1-i1;
		t12 = t2-t8;			t13 = t3-t9;
		t2  = t2+t8;			t3  = t3+t9;
		t8  = t4+t6;			t9  = t5+t7;
		t6  = t4-t6;			t7  = t5-t7;

		/* polypro(A,B) mod P0 + x0 */
		r1 = t10*bp0 + aj1p0r;		i1 = t11*bp0 + aj1p0i;

		/* polypro(A,B) mod P1      */
		t0 = t0 *bp1;			t1 = t1 *bp1;

		aj1p0r += t10;			aj1p0i += t11;

		/* polypro(A,B) mod P2 = x^2-x+1: */
		t4  =(t12+t8 )*bp2B;	t5  =(t13+t9 )*bp2B;
		t12 = t12*bp20;		t13 = t13*bp20;
		t8  = t8 *bp21;		t9  = t9 *bp21;
		t8  = t12 - t8 ;	t9  = t13 - t9 ;
		t4  = t4  - t12;	t5  = t5  - t13;

		/* polypro(A,B) mod P3 = x^2+x+1: */
		t10 =(t2 -t6 )*bp3B;	t11 =(t3 -t7 )*bp3B;
		t2  = t2 *bp30;		t3  = t3 *bp30;
		t6  = t6 *bp31;		t7  = t7 *bp31;
		t6  = t2  - t6 ;	t7  = t3  - t7 ;
		t10 = t10 + t2 ;	t11 = t11 + t3 ;

		/* Combine the 4 modular polynomial product outputs according to CRT. */
		t2 = r1 +t0;		t3 = i1 +t1;
		t0 = r1 -t0;		t1 = i1 -t1;
		r1 = t8 +t6;		i1 = t9 +t7;
		t8 = t8 -t6;		t9 = t9 -t7;
		t6 = t4 +t10;		t7 = t5 +t11;
		t4 = t4 -t10;		t5 = t5 -t11;

#if(LO_ADD)
		t2 = t2 - r1 + t4;	t3 = t3 -i1 + t5;	/* C3 = x^2 = a - c + f		2 add        */
		t12= t2 + 3*r1;		t13= t3 + 3*i1;		/* C1 = x^0 = x^2 + 3.c		1 add, 1 mul */
		t4 = t2 - 3*t4;		t5 = t3 - 3*t5;		/* C4 = x^4 = x^2 - 3.f		1 add, 1 mul */

		t10= t0 + t8 - t6;	t11= t1 + t9 - t7;	/* C2 = x^5 = b + d - e		2 add        */
		t6 = t10+ 3*t6;		t7 = t11+ 3*t7;		/* C6 = x^1 = x^5 + 3.e		1 add, 1 mul */
		t8 = t10- 3*t8;		t9 = t11- 3*t9;		/* C5 = x^3 = x^5 - 3.d		1 add, 1 mul */
#else
		t10 = t2 -r1;		t11 = t3 -i1;
		t12 = t2 +r1 +r1 +t4;	t13 = t3 +i1 +i1 +t5;	/* C1 */
		t2  = t10        +t4;	t3  = t11        +t5;	/* C3 */
		t4  = t10    -t4 -t4;	t5  = t11    -t5 -t5;	/* C4 */

		r1  = t0 +t8;		i1  = t1 +t9;
		t10 = t0 -t8 -t8 -t6;	t11 = t1 -t9 -t9 -t7;	/* C5 */
		t8  = r1     +t6 +t6;	t9  = i1     +t7 +t7;	/* C6 */
		t6  = r1         -t6;	t7  = i1         -t7;	/* C2 */
#endif
	/*** SINE TERMS ***/

		r1 = t14-t18+t22;	i1 = t15-t19+t23;
		t14= t14-t22;		t15= t15-t23;
		t18= t18+t22;		t19= t19+t23;

		t0 = t16-t20+t24;	t1 = t17-t21+t25;
		t16= t16-t24;		t17= t17-t25;
		t20= t20+t24;		t21= t21+t25;

		/* polypro(A,B) mod Q0 = x^2+1: */
		t22 =(r1 +t0 )*bq0B;	t23 =(i1 +t1 )*bq0B;
		r1  = r1 *bq00;		i1  = i1 *bq00;
		t0  = t0 *bq01;		t1  = t1 *bq01;
		t24 = r1  - t0;		t25 = i1  - t1;
		t22 = t22 - r1  - t0;	t23 = t23 - i1  - t1;	/* t0,1 free */

		/* polypro(A,B) mod Q1 = x^4-x^2+1: */
		r1  = t14+t18;		i1  = t15+t19;
		t0  = t14*bq10;		t1  = t15*bq10;
		r2  = t18*bq12;		i2  = t19*bq12;
		r2  = t0  - r2 ;	i2  = t1  - i2 ;
		t0  = r1 *bq1A - t0 ;	t1  = i1 *bq1A - t1 ;

		t14 = t14*bq11;		t15 = t15*bq11;
		t18 = t18*bq13;		t19 = t19*bq13;
		t18 = t14 - t18;	t19 = t15 - t19;
		t14 = r1 *bq1B - t14;	t15 = i1 *bq1B - t15;

		r1  = t16+t20;		i1  = t17+t21;
		r3  = t16*bq10;		i3  = t17*bq10;
		r4  = t20*bq12;		i4  = t21*bq12;
		r4  = r3  - r4;		i4  = i3  - i4;
		r3  = r1 *bq1A - r3 ;	i3  = i1 *bq1A - i3 ;

		t16 = t16*bq11;		t17 = t17*bq11;
		t20 = t20*bq13;		t21 = t21*bq13;
		t20 = t16 - t20;	t21 = t17 - t21;
		t16 = r1 *bq1B - t16;	t17 = i1 *bq1B - t17;

		t20 = t20+t0 +t16;	t21 = t21+t1 +t17;
		t16 = r2 -t16;		t17 = i2 -t17;
		t18 = t18+r4;		t19 = t19+i4;
		t14 = t14+r3 ;		t15 = t15+i3 ;

		/* Combine the 2 modular polynomial product outputs according to CRT. */
#if(LO_ADD)
		t24 = t24 - t16 + t20;	t25 = t25 - t17 + t21;	/* S4 = x^4 = q00 - q10 + q12	2 add        */
		t16 = t24 + 3*t16;	t17 = t25 + 3*t17;	/* S1 = x^0 = x^4 + 3.q10	1 add, 1 mul */
		t20 = 3*t20 - t24;	t21 = 3*t21 - t25;	/* S3 = x^2 = 3.q12 - x^4	1 add, 1 mul */

		t22 = t22 - t18 + t14;	t23 = t23 - t19 + t15;	/* S2 = x^5 = q01 - q11 + q13	2 add        */
		t18 = t22 + 3*t18;	t19 = t23 + 3*t19;	/* S6 = x^1 = x^5 + 3.q11	1 add, 1 mul */
		t14 = 3*t14 - t22;	t15 = 3*t15 - t23;	/* S5 = x^3 = 3.q13 - x^5	1 add, 1 mul */

	/*...Inline multiply of sine parts by +-I into finishing phase, remembering that the acyclic
			 convo algorithm used above returns -s1 and -s5. 	*/

		aj1p1r =t12-t17;	aj1p1i =t13+t16;	/* C1 - S1 */
		aj1p2r =t10+t23;	aj1p2i =t11-t22;	/* C2 - S2 */
		aj1p3r =t2 +t21;	aj1p3i =t3 -t20;	/* C3 - S3 */
		aj1p4r =t4 +t25;	aj1p4i =t5 -t24;	/* C4 - S4 */
		aj1p5r =t8 -t15;	aj1p5i =t9 +t14;	/* C5 - S5 */
		aj1p6r =t6 +t19;	aj1p6i =t7 -t18;	/* C6 - S6 */
		aj1p7r =t6 -t19;	aj1p7i =t7 +t18;	/* C6 + S6 */
		aj1p8r =t8 +t15;	aj1p8i =t9 -t14;	/* C5 + S5 */
		aj1p9r =t4 -t25;	aj1p9i =t5 +t24;	/* C4 + S4 */
		aj1p10r=t2 -t21;	aj1p10i=t3 +t20;	/* C3 + S3 */
		aj1p11r=t10-t23;	aj1p11i=t11+t22;	/* C2 + S2 */
		aj1p12r=t12+t17;	aj1p12i=t13-t16;	/* C1 + S1 */
#else
		r1  = t24+t16;		i1  = t25+t17;
		t24 = t24-t16;		t25 = t25-t17;
		t16 = t16+t20;		t17 = t17+t21;
		t0  = t22+t18;		t1  = t23+t19;
		t22 = t22-t18;		t23 = t23-t19;
		t18 = t18+t14;		t19 = t19+t15;

		t24 = t24+t20;		t25 = t25+t21;		/* S4 */
		t22 = t22+t14;		t23 = t23+t15;		/* S2 */
		t14 = r1 +t16;		t15 = i1 +t17;		/* S1 */
		t16 = t16+t16-r1 ;	t17 = t17+t17-i1 ;	/* S3 */
		t20 = t0 +t18;		t21 = t1 +t19;		/* S6 */
		t18 = t18+t18-t0 ;	t19 = t19+t19-t1 ;	/* S5 */

	/*...Inline multiply of sine parts by +-I into finishing phase, remembering that the acyclic
			 convo algorithm used above returns -s1 and -s5. 	*/

		aj1p1r =t12-t15;	aj1p1i =t13+t14;	/* C1 - S1 */
		aj1p2r =t6 +t23;	aj1p2i =t7 -t22;	/* C2 - S2 */
		aj1p3r =t2 +t17;	aj1p3i =t3 -t16;	/* C3 - S3 */
		aj1p4r =t4 +t25;	aj1p4i =t5 -t24;	/* C4 - S4 */
		aj1p5r =t10-t19;	aj1p5i =t11+t18;	/* C5 - S5 */
		aj1p6r =t8 +t21;	aj1p6i =t9 -t20;	/* C6 - S6 */
		aj1p7r =t8 -t21;	aj1p7i =t9 +t20;	/* C6 + S6 */
		aj1p8r =t10+t19;	aj1p8i =t11-t18;	/* C5 + S5 */
		aj1p9r =t4 -t25;	aj1p9i =t5 +t24;	/* C4 + S4 */
		aj1p10r=t2 -t17;	aj1p10i=t3 +t16;	/* C3 + S3 */
		aj1p11r=t6 -t23;	aj1p11i=t7 +t22;	/* C2 + S2 */
		aj1p12r=t12+t15;	aj1p12i=t13-t14;	/* C1 + S1 */
#endif

/*..Now do the carries. Since the outputs would
    normally be getting dispatched to 13 separate blocks of the A-array, we need 13 separate carries.	*/

			l= j & (nwt-1);
			n_minus_sil   = n-si[l  ];
			n_minus_silp1 = n-si[l+1];
			sinwt   = si[nwt-l  ];
			sinwtm1 = si[nwt-l-1];

			wtl     =wt0[    l  ];
			wtn     =wt0[nwt-l  ]*scale;	/* Include 1/(n/2) scale factor of inverse transform here...	*/
			wtlp1   =wt0[    l+1];
			wtnm1   =wt0[nwt-l-1]*scale;	/* ...and here.	*/

/*...set0 is slightly different from others:	*/
		 cmplx_carry_norm_nocheck0(aj1p0r ,aj1p0i ,cy0 ,bjmodn0 );
			cmplx_carry_norm_nocheck(aj1p1r ,aj1p1i ,cy1 ,bjmodn1 ,1 );
			cmplx_carry_norm_nocheck(aj1p2r ,aj1p2i ,cy2 ,bjmodn2 ,2 );
			cmplx_carry_norm_nocheck(aj1p3r ,aj1p3i ,cy3 ,bjmodn3 ,3 );
			cmplx_carry_norm_nocheck(aj1p4r ,aj1p4i ,cy4 ,bjmodn4 ,4 );
			cmplx_carry_norm_nocheck(aj1p5r ,aj1p5i ,cy5 ,bjmodn5 ,5 );
			cmplx_carry_norm_nocheck(aj1p6r ,aj1p6i ,cy6 ,bjmodn6 ,6 );
			cmplx_carry_norm_nocheck(aj1p7r ,aj1p7i ,cy7 ,bjmodn7 ,7 );
			cmplx_carry_norm_nocheck(aj1p8r ,aj1p8i ,cy8 ,bjmodn8 ,8 );
			cmplx_carry_norm_nocheck(aj1p9r ,aj1p9i ,cy9 ,bjmodn9 ,9 );
			cmplx_carry_norm_nocheck(aj1p10r,aj1p10i,cy10,bjmodn10,10);
			cmplx_carry_norm_nocheck(aj1p11r,aj1p11i,cy11,bjmodn11,11);
			cmplx_carry_norm_nocheck(aj1p12r,aj1p12i,cy12,bjmodn12,12);

			i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
				 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

/*...The radix-13 DIF pass is here:	*/
#if PFETCH
add0 = & a[j2];
prefetch_p_doubles(add0);
#endif
	/*    gather the needed data (13 64-bit complex, i.e. 26 64-bit reals) and get needed xi+-xj combos:	*/
		t2 =aj1p1r+aj1p12r;	t3 =aj1p1i+aj1p12i;	t14=aj1p1r-aj1p12r;	t15=aj1p1i-aj1p12i;
		t4 =aj1p2r+aj1p11r;	t5 =aj1p2i+aj1p11i;	t16=aj1p2r-aj1p11r;	t17=aj1p2i-aj1p11i;
		t10=aj1p3r+aj1p10r;	t11=aj1p3i+aj1p10i;	t22=aj1p3r-aj1p10r;	t23=aj1p3i-aj1p10i;
		t6 =aj1p4r+aj1p9r ;	t7 =aj1p4i+aj1p9i ;	t18=aj1p4r-aj1p9r ;	t19=aj1p4i-aj1p9i ;
		t8 =aj1p5r+aj1p8r ;	t9 =aj1p5i+aj1p8i ;	t20=aj1p8r-aj1p5r ;	t21=aj1p8i-aj1p5i ;
		t12=aj1p6r+aj1p7r ;	t13=aj1p6i+aj1p7i ;	t24=aj1p6r-aj1p7r ;	t25=aj1p6i-aj1p7i ;
#if PFETCH
addr = add0+p1;
prefetch_p_doubles(addr);
#endif
	/*** COSINE TERMS ***/

		t0  = t2 +t6 +t10;		t1  = t3 +t7 +t11;
		t2  = t2 -t6     ;		t3  = t3 -t7     ;
		t6  =     t6 -t10;		t7  =     t7 -t11;

		r1  = t4 +t8 +t12;		i1  = t5 +t9 +t13;
		t4  = t4     -t12;		t5  = t5     -t13;
		t8  =     t8 -t12;		t9  =     t9 -t13;

		t10 = t0+r1;			t11 = t1+i1;
		t0  = t0-r1;			t1  = t1-i1;
		t12 = t2-t8;			t13 = t3-t9;
		t2  = t2+t8;			t3  = t3+t9;
		t8  = t4+t6;			t9  = t5+t7;
		t6  = t4-t6;			t7  = t5-t7;
#if PFETCH
addr = add0+p2;
prefetch_p_doubles(addr);
#endif
		/* polypro(A,B) mod P0 + x0 */
		r1 = t10*bp0 + aj1p0r;		i1 = t11*bp0 + aj1p0i;

		/* polypro(A,B) mod P1      */
		t0 = t0 *bp1;			t1 = t1 *bp1;

		a[j1  ] = aj1p0r+t10;		a[j2] = aj1p0i+t11;

		/* polypro(A,B) mod P2 = x^2-x+1: */
		t4  =(t12+t8 )*bp2B;	t5  =(t13+t9 )*bp2B;
		t12 = t12*bp20;		t13 = t13*bp20;
		t8  = t8 *bp21;		t9  = t9 *bp21;
		t8  = t12 - t8 ;	t9  = t13 - t9 ;
		t4  = t4  - t12;	t5  = t5  - t13;
#if PFETCH
addr = add0+p3;
prefetch_p_doubles(addr);
#endif
		/* polypro(A,B) mod P3 = x^2+x+1: */
		t10 =(t2 -t6 )*bp3B;	t11 =(t3 -t7 )*bp3B;
		t2  = t2 *bp30;		t3  = t3 *bp30;
		t6  = t6 *bp31;		t7  = t7 *bp31;
		t6  = t2  - t6 ;	t7  = t3  - t7 ;
		t10 = t10 + t2 ;	t11 = t11 + t3 ;
#if PFETCH
addr = add0+p4;
prefetch_p_doubles(addr);
#endif
		/* Combine the 4 modular polynomial product outputs according to CRT. */
		t2 = r1 +t0;		t3 = i1 +t1;
		t0 = r1 -t0;		t1 = i1 -t1;
		r1 = t8 +t6;		i1 = t9 +t7;
		t8 = t8 -t6;		t9 = t9 -t7;
		t6 = t4 +t10;		t7 = t5 +t11;
		t4 = t4 -t10;		t5 = t5 -t11;
#if PFETCH
addr = add0+p5;
prefetch_p_doubles(addr);
#endif
#if(LO_ADD)
		t2 = t2 - r1 + t4;	t3 = t3 -i1 + t5;	/* C3 = x^2 = a - c + f		2 add        */
		t12= t2 + 3*r1;		t13= t3 + 3*i1;		/* C1 = x^0 = x^2 + 3.c		1 add, 1 mul */
		t4 = t2 - 3*t4;		t5 = t3 - 3*t5;		/* C4 = x^4 = x^2 - 3.f		1 add, 1 mul */

		t10= t0 + t8 - t6;	t11= t1 + t9 - t7;	/* C2 = x^5 = b + d - e		2 add        */
		t6 = t10+ 3*t6;		t7 = t11+ 3*t7;		/* C6 = x^1 = x^5 + 3.e		1 add, 1 mul */
		t8 = t10- 3*t8;		t9 = t11- 3*t9;		/* C5 = x^3 = x^5 - 3.d		1 add, 1 mul */
#else
		t10 = t2 -r1;		t11 = t3 -i1;
		t12 = t2 +r1 +r1 +t4;	t13 = t3 +i1 +i1 +t5;	/* C1 */
		t2  = t10        +t4;	t3  = t11        +t5;	/* C3 */
		t4  = t10    -t4 -t4;	t5  = t11    -t5 -t5;	/* C4 */

		r1  = t0 +t8;		i1  = t1 +t9;
		t10 = t0 -t8 -t8 -t6;	t11 = t1 -t9 -t9 -t7;	/* C5 */
		t8  = r1     +t6 +t6;	t9  = i1     +t7 +t7;	/* C6 */
		t6  = r1         -t6;	t7  = i1         -t7;	/* C2 */
#endif
#if PFETCH
addr = add0+p6;
prefetch_p_doubles(addr);
#endif

	/*** SINE TERMS ***/

		r1 = t14-t18+t22;	i1 = t15-t19+t23;
		t14= t14-t22;		t15= t15-t23;
		t18= t18+t22;		t19= t19+t23;

		t0 = t16-t20+t24;	t1 = t17-t21+t25;
		t16= t16-t24;		t17= t17-t25;
		t20= t20+t24;		t21= t21+t25;
#if PFETCH
addr = add0+p7;
prefetch_p_doubles(addr);
#endif
		/* polypro(A,B) mod Q0 = x^2+1: */
		t22 =(r1 +t0 )*bq0B;	t23 =(i1 +t1 )*bq0B;
		r1  = r1 *bq00;		i1  = i1 *bq00;
		t0  = t0 *bq01;		t1  = t1 *bq01;
		t24 = r1  - t0;		t25 = i1  - t1;
		t22 = t22 - r1  - t0;	t23 = t23 - i1  - t1;	/* t0,1 free */

		/* polypro(A,B) mod Q1 = x^4-x^2+1: */
		r1  = t14+t18;		i1  = t15+t19;
		t0  = t14*bq10;		t1  = t15*bq10;
		r2  = t18*bq12;		i2  = t19*bq12;
		r2  = t0  - r2 ;	i2  = t1  - i2 ;
		t0  = r1 *bq1A - t0 ;	t1  = i1 *bq1A - t1 ;
#if PFETCH
addr = add0+p8;
prefetch_p_doubles(addr);
#endif
		t14 = t14*bq11;		t15 = t15*bq11;
		t18 = t18*bq13;		t19 = t19*bq13;
		t18 = t14 - t18;	t19 = t15 - t19;
		t14 = r1 *bq1B - t14;	t15 = i1 *bq1B - t15;

		r1  = t16+t20;		i1  = t17+t21;
		r3  = t16*bq10;		i3  = t17*bq10;
		r4  = t20*bq12;		i4  = t21*bq12;
		r4  = r3  - r4;		i4  = i3  - i4;
		r3  = r1 *bq1A - r3 ;	i3  = i1 *bq1A - i3 ;
#if PFETCH
addr = add0+p9;
prefetch_p_doubles(addr);
#endif
		t16 = t16*bq11;		t17 = t17*bq11;
		t20 = t20*bq13;		t21 = t21*bq13;
		t20 = t16 - t20;	t21 = t17 - t21;
		t16 = r1 *bq1B - t16;	t17 = i1 *bq1B - t17;

		t20 = t20+t0 +t16;	t21 = t21+t1 +t17;
		t16 = r2 -t16;		t17 = i2 -t17;
		t18 = t18+r4;		t19 = t19+i4;
		t14 = t14+r3 ;		t15 = t15+i3 ;
#if PFETCH
addr = add0+p10;
prefetch_p_doubles(addr);
#endif
		/* Combine the 2 modular polynomial product outputs according to CRT. */
#if(LO_ADD)
		t24 = t24 - t16 + t20;	t25 = t25 - t17 + t21;	/* S4 = x^4 = q00 - q10 + q12	2 add        */
		t16 = t24 + 3*t16;	t17 = t25 + 3*t17;	/* S1 = x^0 = x^4 + 3.q10	1 add, 1 mul */
		t20 = 3*t20 - t24;	t21 = 3*t21 - t25;	/* S3 = x^2 = 3.q12 - x^4	1 add, 1 mul */

		t22 = t22 - t18 + t14;	t23 = t23 - t19 + t15;	/* S2 = x^5 = q01 - q11 + q13	2 add        */
		t18 = t22 + 3*t18;	t19 = t23 + 3*t19;	/* S6 = x^1 = x^5 + 3.q11	1 add, 1 mul */
		t14 = 3*t14 - t22;	t15 = 3*t15 - t23;	/* S5 = x^3 = 3.q13 - x^5	1 add, 1 mul */

	/*...Inline multiply of sine parts by +-I into finishing phase, remembering that the acyclic
			 convo algorithm used above returns -s1 and -s5. 	*/
#if PFETCH
addr = add0+p11;
prefetch_p_doubles(addr);
#endif
		a[j1+p1 ]=t12+t17;	a[j2+p1 ]=t13-t16;	/* C1 - S1 */
		a[j1+p2 ]=t10-t23;	a[j2+p2 ]=t11+t22;	/* C2 - S2 */
		a[j1+p3 ]=t2 -t21;	a[j2+p3 ]=t3 +t20;	/* C3 - S3 */
		a[j1+p4 ]=t4 -t25;	a[j2+p4 ]=t5 +t24;	/* C4 - S4 */
		a[j1+p5 ]=t8 +t15;	a[j2+p5 ]=t9 -t14;	/* C5 - S5 */
		a[j1+p6 ]=t6 -t19;	a[j2+p6 ]=t7 +t18;	/* C6 - S6 */
		a[j1+p7 ]=t6 +t19;	a[j2+p7 ]=t7 -t18;	/* C6 + S6 */
		a[j1+p8 ]=t8 -t15;	a[j2+p8 ]=t9 +t14;	/* C5 + S5 */
		a[j1+p9 ]=t4 +t25;	a[j2+p9 ]=t5 -t24;	/* C4 + S4 */
		a[j1+p10]=t2 +t21;	a[j2+p10]=t3 -t20;	/* C3 + S3 */
		a[j1+p11]=t10+t23;	a[j2+p11]=t11-t22;	/* C2 + S2 */
		a[j1+p12]=t12-t17;	a[j2+p12]=t13+t16;	/* C1 + S1 */
#else
		r1  = t24+t16;		i1  = t25+t17;
		t24 = t24-t16;		t25 = t25-t17;
		t16 = t16+t20;		t17 = t17+t21;
		t0  = t22+t18;		t1  = t23+t19;
		t22 = t22-t18;		t23 = t23-t19;
		t18 = t18+t14;		t19 = t19+t15;

		t24 = t24+t20;		t25 = t25+t21;		/* S4 */
		t22 = t22+t14;		t23 = t23+t15;		/* S2 */
		t14 = r1 +t16;		t15 = i1 +t17;		/* S1 */
		t16 = t16+t16-r1 ;	t17 = t17+t17-i1 ;	/* S3 */
		t20 = t0 +t18;		t21 = t1 +t19;		/* S6 */
		t18 = t18+t18-t0 ;	t19 = t19+t19-t1 ;	/* S5 */

	/*...Inline multiply of sine parts by +-I into finishing phase, remembering that the acyclic
			 convo algorithm used above returns -s1 and -s5. 	*/
#if PFETCH
addr = add0+p11;
prefetch_p_doubles(addr);
#endif
		a[j1+p1 ]=t12+t15;	a[j2+p1 ]=t13-t14;	/* C1 - S1 */
		a[j1+p2 ]=t6 -t23;	a[j2+p2 ]=t7 +t22;	/* C2 - S2 */
		a[j1+p3 ]=t2 -t17;	a[j2+p3 ]=t3 +t16;	/* C3 - S3 */
		a[j1+p4 ]=t4 -t25;	a[j2+p4 ]=t5 +t24;	/* C4 - S4 */
		a[j1+p5 ]=t10+t19;	a[j2+p5 ]=t11-t18;	/* C5 - S5 */
		a[j1+p6 ]=t8 -t21;	a[j2+p6 ]=t9 +t20;	/* C6 - S6 */
		a[j1+p7 ]=t8 +t21;	a[j2+p7 ]=t9 -t20;	/* C6 + S6 */
		a[j1+p8 ]=t10-t19;	a[j2+p8 ]=t11+t18;	/* C5 + S5 */
		a[j1+p9 ]=t4 +t25;	a[j2+p9 ]=t5 -t24;	/* C4 + S4 */
		a[j1+p10]=t2 +t17;	a[j2+p10]=t3 -t16;	/* C3 + S3 */
		a[j1+p11]=t6 +t23;	a[j2+p11]=t7 -t22;	/* C2 + S2 */
		a[j1+p12]=t12-t15;	a[j2+p12]=t13+t14;	/* C1 + S1 */
#endif
#if PFETCH
addr = add0+p12;
prefetch_p_doubles(addr);
#endif
		iroot += root_incr;		/* increment sincos index.	*/

		}

		jstart += nwt;
		jhi    += nwt;
		col += 13;
		co3 -= 13;

	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-13 forward DIF FFT of the first block of 13 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 13 outputs of (1);
!   (3) Reweight and perform a radix-13 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 13 elements and repeat (1-4).
*/

	t1  = cy12;
	cy12= cy11;
	cy11= cy10;
	cy10= cy9;
	cy9 = cy8;
	cy8 = cy7;
	cy7 = cy6;
	cy6 = cy5;
	cy5 = cy4;
	cy4 = cy3;
	cy3 = cy2;
	cy2 = cy1;
	cy1 = cy0;
	cy0 = t1;

	iroot = 0;
	root_incr = 0;
	scale = 1;

	jstart = 0;
	jhi = 7;
	khi = 1;

	for(j=0; j<=jhi; j++)
	{
		a[j    ] *= radix_inv;
		a[j+p1 ] *= radix_inv;
		a[j+p2 ] *= radix_inv;
		a[j+p3 ] *= radix_inv;
		a[j+p4 ] *= radix_inv;
		a[j+p5 ] *= radix_inv;
		a[j+p6 ] *= radix_inv;
		a[j+p7 ] *= radix_inv;
		a[j+p8 ] *= radix_inv;
		a[j+p9 ] *= radix_inv;
		a[j+p10] *= radix_inv;
		a[j+p11] *= radix_inv;
		a[j+p12] *= radix_inv;
	}
}

	if(fabs(cy0)+fabs(cy1)+fabs(cy2)+fabs(cy3)+fabs(cy4)+fabs(cy5)+fabs(cy6)+fabs(cy7)+fabs(cy8)+fabs(cy9)+fabs(cy10)+fabs(cy11)+fabs(cy12) != 0.0)
	{
			sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix13_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
			if(INTERACT)fprintf(stderr,"%s",cbuf);
			fp = fopen(   OFILE,"a");
			fq = fopen(STATFILE,"a");
		fprintf(fp,"%s",cbuf);
		fprintf(fq,"%s",cbuf);
			fclose(fp);	fp = 0x0;
			fclose(fq);	fq = 0x0;
			err=ERR_CARRY;
			return(err);
	}
	return(0);
}

/***************/

void radix13_dif_pass1(double a[], int n)
{
/*
!
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-13 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.

================================================================================================

Radix-13 DFT using CRT for 6x6 subconvos:

The cosine terms of the output need the following combos:

c1*(x1+x12)+c2*(x2+x11)+c3*(x3+x10)+c4*(x4+x9)+c5*(x5+x8)+c6*(x6+x7)
c2*(x1+x12)+c4*(x2+x11)+c6*(x3+x10)+c5*(x4+x9)+c3*(x5+x8)+c1*(x6+x7)
c3*(x1+x12)+c6*(x2+x11)+c4*(x3+x10)+c1*(x4+x9)+c2*(x5+x8)+c5*(x6+x7)
c4*(x1+x12)+c5*(x2+x11)+c1*(x3+x10)+c3*(x4+x9)+c6*(x5+x8)+c2*(x6+x7)
c5*(x1+x12)+c3*(x2+x11)+c2*(x3+x10)+c6*(x4+x9)+c1*(x5+x8)+c4*(x6+x7)
c6*(x1+x12)+c1*(x2+x11)+c5*(x3+x10)+c2*(x4+x9)+c4*(x5+x8)+c3*(x6+x7)

where cX := cos(X*2*pi/13). These give no hint as to what permutation of input and/or trig. term indices we'll
need in order to put this in the form of a convolution. For the sine terms things seem even more dicy, since there
we not only have to put the index pattern in the form of a convolution (and the same perutation that accomplishes
this for the cosine terms also does the same for the sine terms), but we also need to get the + and - signs
in the form of an acayclic convo. This is because (unlike radix-11) there doesn't appear to be a pattern of
sign flips on the inputs ((xi-xj) and sine multipliers) which makes the sine part look like a cyclic convo.
But first let's deal just with the index permutation: my program permute.c gives the permutation (and perhaps
it's not unique, since the program stops as soon as it finds the first successful permutation) (0,5,2,4,3,1)
as one which yields a convolution index pattern. Let's apply it to the sine terms, defined by

+s1*(x1-x12)+s2*(x2-x11)+s3*(x3-x10)+s4*(x4-x9)+s5*(x5-x8)+s6*(x6-x7) = S1
+s2*(x1-x12)+s4*(x2-x11)+s6*(x3-x10)-s5*(x4-x9)-s3*(x5-x8)-s1*(x6-x7) = S2
+s3*(x1-x12)+s6*(x2-x11)-s4*(x3-x10)-s1*(x4-x9)+s2*(x5-x8)+s5*(x6-x7) = S3
+s4*(x1-x12)-s5*(x2-x11)-s1*(x3-x10)+s3*(x4-x9)-s6*(x5-x8)-s2*(x6-x7) = S4
+s5*(x1-x12)-s3*(x2-x11)+s2*(x3-x10)-s6*(x4-x9)-s1*(x5-x8)+s4*(x6-x7) = S5
+s6*(x1-x12)-s1*(x2-x11)+s5*(x3-x10)-s2*(x4-x9)+s4*(x5-x8)-s3*(x6-x7) = S6 .

Letting s(1,2,3,4,5,6)
      = b(0,5,2,4,3,1) gives

+b0*(x1-x12)+b5*(x2-x11)+b2*(x3-x10)+b4*(x4-x9)+b3*(x5-x8)+b1*(x6-x7) = S1
+b5*(x1-x12)+b4*(x2-x11)+b1*(x3-x10)-b3*(x4-x9)-b2*(x5-x8)-b0*(x6-x7) = S2
+b2*(x1-x12)+b1*(x2-x11)-b4*(x3-x10)-b0*(x4-x9)+b5*(x5-x8)+b3*(x6-x7) = S3
+b4*(x1-x12)-b3*(x2-x11)-b0*(x3-x10)+b2*(x4-x9)-b1*(x5-x8)-b5*(x6-x7) = S4
+b3*(x1-x12)-b2*(x2-x11)+b5*(x3-x10)-b1*(x4-x9)-b0*(x5-x8)+b4*(x6-x7) = S5
+b1*(x1-x12)-b0*(x2-x11)+b3*(x3-x10)-b5*(x4-x9)+b4*(x5-x8)-b2*(x6-x7) = S6 ,

and rearranging the rows to reorder the b-terms of the leftmost column yields

+b0*(x1-x12)+b5*(x2-x11)+b2*(x3-x10)+b4*(x4-x9)+b3*(x5-x8)+b1*(x6-x7) = S1
+b1*(x1-x12)-b0*(x2-x11)+b3*(x3-x10)-b5*(x4-x9)+b4*(x5-x8)-b2*(x6-x7) = S6
+b2*(x1-x12)+b1*(x2-x11)-b4*(x3-x10)-b0*(x4-x9)+b5*(x5-x8)+b3*(x6-x7) = S3
+b3*(x1-x12)-b2*(x2-x11)+b5*(x3-x10)-b1*(x4-x9)-b0*(x5-x8)+b4*(x6-x7) = S5
+b4*(x1-x12)-b3*(x2-x11)-b0*(x3-x10)+b2*(x4-x9)-b1*(x5-x8)-b5*(x6-x7) = S4
+b5*(x1-x12)+b4*(x2-x11)+b1*(x3-x10)-b3*(x4-x9)-b2*(x5-x8)-b0*(x6-x7) = S2

We can see that the columns have convolution index patterns, i.e. the six distinct
circular shifts of the vector (0,1,2,3,4,5). Next, letting

a0 = (x1-x12), a1 = (x2-x11), a4 = (x3-x10), a2 = (x4-x9), a3 = (x5-x8), a5 = (x6-x7),
(and doing analogously for the cosine terms, except there we don't need to worry about signs) we get

+a0.b0+a1.b5+a2.b4+a3.b3+a4.b2+a5.b1 = S1
+a0.b1-a1.b0-a2.b5+a3.b4+a4.b3-a5.b2 = S6
+a0.b2+a1.b1-a2.b0+a3.b5-a4.b4+a5.b3 = S3
+a0.b3-a1.b2-a2.b1-a3.b0+a4.b5+a5.b4 = S5
+a0.b4-a1.b3+a2.b2-a3.b1-a4.b0-a5.b5 = S4
+a0.b5+a1.b4-a2.b3-a3.b2+a4.b1-a5.b0 = S2 .

Now we need to generate the proper sign pattern for an acyclic.

Flipping the sign of b0 (i.e. s1) and then flipping the sign of the entire first row gives:

+a0.b0-a1.b5-a2.b4-a3.b3-a4.b2-a5.b1 =-S1;	a0, b0 frozen
+a0.b1+a1.b0-a2.b5+a3.b4+a4.b3-a5.b2 = S6;	need to flip a3.b4,a4.b3
+a0.b2+a1.b1+a2.b0+a3.b5-a4.b4+a5.b3 = S3;	need to flip a5.b3,a3.b5
+a0.b3-a1.b2-a2.b1+a3.b0+a4.b5+a5.b4 = S5;	need to flip a0.b3,a3.b0, then entire row
+a0.b4-a1.b3+a2.b2-a3.b1+a4.b0-a5.b5 = S4;	need to flip a1.b3,a3.b1
+a0.b5+a1.b4-a2.b3-a3.b2+a4.b1+a5.b0 = S2;	need to flip a2.b3,a3.b2

Note that every one of the needed 2-term sign flips involves a3 in one term and b3 in the other, and we can flip both
without affecting the already-finished first row, so flipping signs of a3 and b3 (i.e. s5) and then of row 4:

+a0.b0-a1.b5-a2.b4-a3.b3-a4.b2-a5.b1 =-S1
+a0.b1+a1.b0-a2.b5-a3.b4-a4.b3-a5.b2 = S6
+a0.b2+a1.b1+a2.b0-a3.b5-a4.b4-a5.b3 = S3
+a0.b3+a1.b2+a2.b1+a3.b0-a4.b5-a5.b4 =-S5
+a0.b4+a1.b3+a2.b2+a3.b1+a4.b0-a5.b5 = S4
+a0.b5+a1.b4+a2.b3+a3.b2+a4.b1+a5.b0 = S2 .

The LHS is now in the form of a length-6 acyclic convolution. We next turn to the issue of how to efficiently perform the
two needed convolutions (length-6 cyclic for cosine terms, length-6 acyclic for sine terms), which involves factoring the
two underlying order-6 polynomials into their irreducible components and then applying the Chinese Remainder Theorem (CRT).


Given variable A = [a0,a1,a2,a3,a4,a5], constant B = [b0,b1,b2,b3,b4,b5],
compute convo(A,B) modulo P := x^6 - 1 (cyclic) and modulo Q := x^6 + 1 (acyclic).
As our trial vectors for debugging, we use A = [3,-1,4,-1,-5,-9] and B = [2,-7,-1,8,2,8],
and (for checksum purposes) assume these digits are with respect to a base-10 expansion,
i.e.
	A = 3+10*(-1+10*( 4+10*(-1+10*(-5+10*(-9))))) = -950607
	B = 2+10*(-7+10*(-1+10*( 8+10*( 2+10*( 8))))) = +827832 .

********* CYCLIC CONVO: **********

Outputs are, in terms of coefficients of various powers of x:

x^0: a0.b0+a1.b5+a2.b4+a3.b3+a4.b2+a5.b1 =  66
x^1: a0.b1+a1.b0+a2.b5+a3.b4+a4.b3+a5.b2 = -24
x^2: a0.b2+a1.b1+a2.b0+a3.b5+a4.b4+a5.b3 = -78
x^3: a0.b3+a1.b2+a2.b1+a3.b0+a4.b5+a5.b4 = -63
x^4: a0.b4+a1.b3+a2.b2+a3.b1+a4.b0+a5.b5 = -81
x^5: a0.b5+a1.b4+a2.b3+a3.b2+a4.b1+a5.b0 =  72


P factors into x^6-1 = (x^3-1)*(x^3+1) = (x-1)*(x+1)*(x^2-x+1)*(x^2+x+1) := P0*P1*P2*P3,
so first reduce both input polys modulo the factors, using that:

modulo P0, x^n == +1 ;
modulo P1, x^n == +1 if n even, -1 if n odd ;
modulo P2, x^2 == +x-1, x^3 == -1, x^4 == -x, x^5 == -x^2 == -x+1 ;
modulo P3, x^2 == -x-1, x^3 == +1, x^4 == +x, x^5 == +x^2 == -x-1 .

											COST:
A mod P0 :=     x-1 = a0+a1+a2+a3+a4+a5			-9	5 ADD (form 0+2+4, 1+3+5, add)
      P1 :=     x+1 = a0-a1+a2-a3+a4-a5			13	1 ADD (subtract the two 3-term sums found above)
      P2 := x^2-x+1 = [a0-a2-a3+a5] + [a1+a2-a4-a5]*x	-9+17x	8 ADD (form d=0-2, e=3-5, f=1-5, g=2-4,
      P3 := x^2+x+1 = [a0-a2+a3-a5] + [a1-a2+a4-a5]*x	 7-  x		 then A mod P2 = (d-e)+(f+g)*x
												A mod P3 = (d+e)+(f-g)*x .

similar for B, but we can precompute the latter, so they're basically free.	12, -6, 3-18x, 3-12x

Check: for x := 10,
A%(x-1) = A%9 = 0 == 9, A%(x+1) = A%11 = 2 == 13, A%(x^2-x+1) = A%91 = -21 == -9+17*10, A%(x^2+x+1) = A%111 = -3 = 7-10, ok.
B%(x-1) = B%9 = 3 ==12, B%(x+1) = B%11 = 5 == -6, B%(x^2-x+1) = B%91 =   5 ==  3-18*10, B%(x^2+x+1) = B%111 =105 = 3-12*10, ok.

				COST:
polypro(A,B) mod P0:	0 ADD, 1 MUL, output = p00		-108 ==   0 mod   9, ok
                 P1:	0 ADD, 1 MUL, output = p10		-78  ==  10 mod  11, ok
                 P2:	3 ADD, 3 MUL, output = p20 + x*p21	(-9.3) + (9.18+17.3)*x + (-17.18).x^2 = -27 + 213.x - 306.x^2
                                                             == -27 + 213.x - 306.(x-1) modulo P2
                                                              =  279 - 93.x == -14 mod  91, ok
                 P3:	3 ADD, 3 MUL, output = p30 + x*p31 .	( 7.3) + (-7.12-1.3)*x + (  1.12).x^2 =  21 -  87.x +  12.x^2
                                                             ==  21 -  87.x -  12.(x+1) modulo P3
                                                              =   9 -  99.x ==  18 mod 111, ok
LENGTH-2 SUBCONVOS: full-length polypro output is

(a+b.x)*(c+d.x) = a.c + (b.c+a.d).x + b.d.x^2, cyclic convo is modulo x^2-1,
 i.e. folds the x^2 coeff. in with the x^0, giving (a.c+b.d) + (b.c+a.d).x .

Now,

(a+b).(c+d) = a.c+a.d+b.c+b.d
(a-b).(c-d) = a.c-a.d-b.c+b.d

so to get convo using just 2 muls, precompute e = (c+d)/2, f = (c-d)/2 (free), then calculate

y0 = (a+b).e
y1 = (a-b).f
z0 = y0+y1	(x^0 output)
z1 = y0-y1	(x^1 output).

Now, modulo P2 = x^2-x+1, our convo output is (a.c-b.d) + (b.c+a.d + b.d).x,
i.e. x^0 has 2.b.d subtracted w.r.to convo output, x^1 has extra b.d added in.
In fact this looks like a complex multiply (x = I) with the extra b.d. added to the imaginary output.
So we use a modified Karatsuba: precompute (c+d), then calculate
					a=-9, b=17, c=3, d=-18
m0 = a.c				-27		(d-e).(m-n)
m1 = b.d				-306		(f+g).(o+p)	a+b = d-e+f+g = (d+f) + (g-e)
y0 = (a+b).(c+d) = a.c+a.d+b.c+b.d	8.-15 = -120	NB: a+b = [a0-a2-a3+a5] + [a1+a2-a4-a5] = a0+a1-a3-a4 = (a0-a4) + (a1-a3) = u + v
							similarly, c+d = x + y, so the product = u.x+v.x+u.y+v.y
z0 = m0-m1	(x^0 output)		279
z1 = y0-m	(x^1 output).		-93	COST: 3 mul, 3 add	(same # total ops as length-2 cyclic convo, but equal # of mul and add.)

Modulo P3 = x^2+x+1, our convo output is (a.c-b.d) + (b.c+a.d - b.d).x,
i.e. x^0 has 2.b.d subtracted w.r.to convo output, x^1 has b.d subtracted.
To get this, we precompute (d-c), then calculate
					a=7, b=-1, c=3, d=-12
m0 = a.c				21		(d+e).(m+n)
m1 = b.d				12		(g-f).(p-o)	a-b = d+e+f-g = (d+f) - (g-e)
y0 = (a-b).(d-c) =-a.c+a.d+b.c-b.d	8.-15 = -120	NB: a-b = [a0-a2+a3-a5] - [a1-a2+a4-a5] = a0-a1+a3-a4 = (a0-a4) - (a1-a3) = u - v
							similarly, d-c = y - x, so the product =-u.x+v.x+u.y-v.y
z0 = m0-m1	(x^0 output)		9
z1 = y0+m0	(x^1 output).		-99	COST: 3 mul, 3 add


Lastly, need to combine the 4 modular polynomial product outputs according to CRT.
First, find the four inverse polys, t0-3:

t0*P1*P2*P3 == 1 modulo P0.
Now:
  P1*P2*P3 = (x+1)*(x^4+x^2+1) = x^5+x^4+x^3+x^2+x+1 == +6 modulo x-1, so t0 = +1/6.

t1*P0*P2*P3 == 1 modulo P1.
Now:
  P0*P2*P3 = (x-1)*(x^4+x^2+1) = x^5-x^4+x^3-x^2+x-1 == -6 modulo x+1, so t1 = -1/6.

t2*P0*P1*P3 == 1 modulo P2.
Now:
  P0*P1*P3 = (x^2-1)*(x^2+x+1) =     x^4+x^3    -x-1 == -2x-2  modulo x^2-x+1. Try multiplying by x:
x*P0*P1*P3                     = x^5+x^4    -x^2-x   == (-x+1)-x-(x-1)-x = -4x+2,
so (x-2)*P0*P1*P3 == -4x+2+4x+4 = 6, so t2 = (x-2)/6.

t3*P0*P1*P2 == 1 modulo P3.
Now:
  P0*P1*P2 = (x^2-1)*(x^2-x+1) =     x^4-x^3    +x-1 == +2x-2  modulo x^2+x+1. Try multiplying by x:
x*P0*P1*P2                     = x^5-x^4    +x^2-x   == -4x-2,
so (x+2)*P0*P1*P3 == -6, so t3 = -(x+2)/6.


Next, find the s-terms:

s0 := t0*P1*P2*P3 =  (x^5+x^4+x^3+x^2+x+1)/6 ,
s1 := t1*P0*P2*P3 = -(x^5-x^4+x^3-x^2+x-1)/6 ,
s2 := t2*P0*P1*P3 = (x-2)*(x^4+x^3    -x-1)/6 = (x^5-x^4-2.x^3-x^2+x+2)/6 ,
s3 := t3*P0*P1*P2 =-(x+2)*(x^4-x^3    +x-1)/6 =-(x^5+x^4-2.x^3+x^2+x-2)/6 ,

Then, the polynomial product mod (x^6 - 1) = p00*s0 + p10*s1 + (p20 + x*p21)*s2 + (p30 + x*p31)*s3,
where the signature of each p-term times the respective s-polynomial is

p00: s0*6    =  1  +x+x^2  +x^3  +x^4+x^5
p10: s1*6    =  1  -x+x^2  -x^3  +x^4-x^5
p20: s2*6    =  2  +x-x^2-2.x^3  -x^4+x^5
p21: s2*6*x ==  1+2.x+x^2  -x^3-2.x^4-x^5 mod (x^6 - 1)
p30: s3*6    =  2  -x-x^2+2.x^3  -x^4-x^5
p31: s3*6*x == -1+2.x-x^2  -x^3+2.x^4-x^5 mod (x^6 - 1) .

So, crunching all the numbers and dividing by 6 (the divide can be absorbed into the b-constants,
so is free) gives the following output coefficients, where we group p00&p10, p20&p30, p21&p31 by letting
a,b = (p00+-p10)/6, c,d = (p20+-p30)/6, e,f = (p21+-p31)/6                   p00 = -108, p10 = -78		COST: 6 ADD
                                                                             p20 =  279, p30 =   9
                                                                             p21 = - 93, p31 = -99, so
											 a = -31, b= -5, c = +48, d = +45, e = -32, f = +1

x^0: a + 2.c +   f					66
x^1: b +   d + 2.e				 -24
x^2: a -   c +   f = (x^0 coeff) - 3.(c)	 -78
x^3: b - 2.d -   e = (x^5 coeff) - 3.(d)	 -63
x^4: a -   c - 2.f = (x^2 coeff) - 3.(f)	 -81
x^5: b +   d -   e = (x^1 coeff) - 3.(e)		72, looks good!

so the total cost of the reconstruction phase = 16 add, and we can save one add by calculating (a-c)
and using it in both the x^2 and x^4 terms, and can save another add by calculating (b+d)
and using it in both the x^1 and x^5 terms, or using (b-e) in both the x^3 and x^5 terms.			COST: 14 ADD
We can also use the following sequence, which trades 4 multiplies for 6 adds, and thus
might be preferable for a floating-point implementation, where the muls could be done
alongside the adds:

x^2 = a - c + f		2 add
x^0 = x^2 + 3.c		1 add, 1 mul
x^4 = x^2 - 3.f		1 add, 1 mul
x^5 = b + d - e		2 add
x^1 = x^5 + 3.e		1 add, 1 mul
x^3 = x^5 - 3.d		1 add, 1 mul

GRAND TOTAL: 40 add, 8 mul (or 34 add, 12 mul), compared to 30 add, 36 mul for the naive scheme.
We can save some further adds by using that the X0 output of the radix-13 DFT = x0 + (A mod P0), i.e. needs just 1 more add,
and by rolling the x0 term into the above calculation (i.e. it needs just
1 add to include the real part of x0 in the p00 output, which then propagates it to all 6 x^k coeffs),
and lastly by using the fact that P1 also appears in the X0 output computation (X0 = x0 + A mod P1).


********* ACYCLIC CONVO: **********

For A = [3,-1,4,-1,-5,-9] and B = [2,-7,-1,8,2,8], outputs are, in terms of coefficients of various powers of x:

x^0: a0.b0-a1.b5-a2.b4-a3.b3-a4.b2-a5.b1 = -54
x^1: a0.b1+a1.b0-a2.b5-a3.b4-a4.b3-a5.b2 = -22
x^2: a0.b2+a1.b1+a2.b0-a3.b5-a4.b4-a5.b3 = 102
x^3: a0.b3+a1.b2+a2.b1+a3.b0-a4.b5-a5.b4 =  53
x^4: a0.b4+a1.b3+a2.b2+a3.b1+a4.b0-a5.b5 =  63
x^5: a0.b5+a1.b4+a2.b3+a3.b2+a4.b1+a5.b0 =  72


Q factors into x^6+1 = (x^2+1)*(x^4-x^2+1) := Q0*Q1,
so first reduce both input polys modulo the factors, using that:

modulo Q0, x^2 == -1, x^3 == -x, x^4 == +1, x^5 == +x ;
modulo Q1, x^4 == x^2-1, x^5 == x^3-x .

Check: for x := 10,
A%(x^2+1) = A%101 =  5, A%(x^4-x^2+1) = A%9901 =  -111
B%(x^2+1) = B%101 = 36, B%(x^4-x^2+1) = B%9901 = +6049
A*B%(x^2+1) = 5.36%101 = 79, A*B%(x^4-x^2+1) = 1829.
													COST:
A mod Q0 :=     x^2+1 = [a0-a2+a4] + [a1-a3+a5]*x				-6 -9.x			4 ADD
      Q1 := x^4-x^2+1 = [a0-a4] + [a1-a5]*x + [a2+a4]*x^2 + [a3+a5]*x^3		8+8.x-x^2-10.x^3	4 ADD

similar for B, but we can precompute the latter, so they're basically free.		B mod Q0 = 5-7.x, B mod Q1 = 0-15.x+x^2+16.x^3

				COST:
convo(A,B) mod Q0:	 3 MUL, 4 ADD, output = q00 + x*q01
               Q1:	12 MUL,15 ADD, output = q10 + x*q11 + x^2*q21 + x^3*q31.

LENGTH-2 SUBCONVO: full-length polypro output is

(a+b.x)*(c+d.x) = a.c + (b.c+a.d).x + b.d.x^2, acyclic convo is modulo x^2+1,
 i.e. folds the -x^2 coeff. in with the x^0, giving (a.c-b.d) + (b.c+a.d).x , which is equivalent to complex multiply with constant b,c.
So we use Karatsuba: precompute (c+d), then calculate
						a=-6, b=-9, c=5, d=-7
m0 = a.c					-30
m1 = b.d					63
y0 = (a+b).(c+d) = a.c+a.d+b.c+b.d		30
z0 = m0-m1	(x^0 output)			q00 = -93
z1 = y0-m0-m1	(x^1 output).			COST: 3 mul, 4 add	q01 =  -3, q00+q01.10 = -93-30 == 79 mod 101, ok.

LENGTH-4 SUBCONVO: full-length polypro output is

(a+b.x+c.x^2+d.x^3)*(e+f.x+g.x^2+h.x^3) = a.e + (a.f+b.e).x + (a.g+b.f+c.e).x^2 + (a.h+b.g+c.f+d.e).x^3 + (b.h+c.g+d.f).x^4 + (c.h+d.g).x^5 + d.h.x^6,

Length-4 cyclic convo is modulo x^4-1, i.e. folds this in such a fashion as to give the following output coefficients:

x^0: a.e+b.h+c.g+d.f
x^1: a.f+b.e+c.h+d.g
x^2: a.g+b.f+c.e+d.h
x^3: a.h+b.g+c.f+d.e

Length-4 acyclic convo is modulo x^4+1, i.e. folds this in such a fashion as to give the following output coefficients:

x^0: a.e-b.h-c.g-d.f
x^1: a.f+b.e-c.h-d.g
x^2: a.g+b.f+c.e-d.h
x^3: a.h+b.g+c.f+d.e

modulo Q1 = x^4-x^2+1, the polypro outputs are (using that x^4 == x^2-1, x^5 == x^3-x, x^6 == x^4-x^2 == -1 mod Q1)
						a,b,c,d = 8,8,-1,-10	e,f,g,h = 0,-15,+1,16
x^0: a.e-b.h-c.g-d.f - d.h			8.  0 - 8. 16 - -1. +1 - -10.-15 - (-10.16) =    0-128+  1-150 + 160 = -117
x^1: a.f+b.e-c.h-d.g				8.-15 + 8.  0 - -1. 16 - -10. +1            = -120+  0+ 16+ 10       =  -94
x^2: a.g+b.f+c.e + (b.h+c.g+d.f)		8. +1 + 8.-15 + -1.  0 + (8.16-1.1+10.15)   =   +8-120+  0 + 279     = +165
x^3: a.h+b.g+c.f+d.e + (c.h+d.g)		8. 16 + 8. +1 + -1.-15 + -10.  0 + (-16-10) =  128+  8+ 15+  0 -  26 = +125

(a+b*x+c*x^2+d*x^3)*(e+f*x+g*x^2+h*x^3)%9901 = 1829, ok.

****************
08/26/02: Try the following sequence:

1) Do length-4 cyclic convo (5 mul, 15 add)
2) Via side calculation, obtain x = (d.h), y = (c.h+d.g), z = (b.h+c.g+d.f); (naive opcount is 6 mul, 3 add)
3) Obtain polypro modulo Q1 by modifying the cyclic convo outputs as follows:

	(x^0 term) - 2.z - x;
	(x^1 term) - 2.y
	(x^2 term) - 2.x
	(x^3 term)       + y, which needs an additional (3 mul, 5 add) or (8 add).

Thus, via this route the mod-Q1 polypro costs (14 mul, 23 add) or (11 mul, 26 add). Both are much add-heavier than we'd like.

****************

...doesn't look very promising. Unless can find an algorithm in Nuss. for polypro modulo Q1,
perhaps should try separating even and odd-order terms and then doing 4 subconvos modulo x^2-x+1:

									intermediates needed:			total ops needed:
(a+c.x)*(e+g.x) == (a.e-c.g) + (c.e+a.g + c.g).x mod x^2-x+1 = m0+n0.x	a.e, c.g, a+c plus 1 mul, 2 add		3 mul, 3 add		+1   +7.x
(a+c.x)*(f+h.x) == (a.f-c.h) + (c.f+a.h + c.h).x mod x^2-x+1 = m1+n1.x	a.f, c.h      plus 1 mul, 2 add		3 mul, 2 add	-104 +127.x
(b+d.x)*(e+g.x) == (b.e-d.g) + (d.e+b.g + d.g).x mod x^2-x+1 = m2+n2.x	b.e, d.g, b+d plus 1 mul, 2 add		3 mul, 3 add	 +10   -2.x
(b+d.x)*(f+h.x) == (b.f-d.h) + (d.f+b.h + d.h).x mod x^2-x+1 = m3+n3.x	b.f, d.h,     plus 1 mul, 2 add		3 mul, 2 add	 +40 +118.x

NB: a+c = a0+a2, b+d = a1+a3, more accurate to get directly from the a-inputs this way.

%%%%%%%%%%%%%% STUFF to TRY %%%%%%%%%%%%%%%%

Rewriting the mod-Q1 convo in terms of indexed vectors (a,b,c,d) = a0-3 and (e,f,g,h) = b0-3
may make it easier to look for structure in the computation:

x^0: a0.b0-a1.b3-a2.b2-a3.b1 - a3.b3
x^1: a0.b1+a1.b0-a2.b3-a3.b2
x^2: a0.b2+a1.b1+a2.b0       + (a1.b3+a2.b2+a3.b1)
x^3: a0.b3+a1.b2+a2.b1+a3.b0 + (a2.b3+a3.b2)

Note: all elts except b1 appear either 4 or 6 times; b1 appears 5 times.

Now, regular length-4 acyclic (i.e. convo modulo x^4+1) gives

x^0: a0.b0-a1.b3-a2.b2-a3.b1		need to sub a3.b3
x^1: a0.b1+a1.b0-a2.b3-a3.b2		no change needed
x^2: a0.b2+a1.b1+a2.b0-a3.b3		need to add a3.b3 + (a1.b3+a2.b2+a3.b1), which = -[(x^0 term) - a0.b0 - a3.b3]
x^3: a0.b3+a1.b2+a2.b1+a3.b0		need to add (a2.b3+a3.b2)
					TOTAL: 4 mul, 6 add to effect corrections

If replace inputs b0,b1,b2,b3 with b0+b2,b1+b3,b2,b3, then get

x^0: a0.b0-a1.b3-a2.b2-a3.b1 + a0.b2-a3.b3	need to sub a0.b2
x^1: a0.b1+a1.b0-a2.b3-a3.b2 + a0.b3+a1.b2	need to sub a0.b3 and a1.b2
x^2: a0.b2+a1.b1+a2.b0-a3.b3 + a1.b3+a2.b2	need to add a3.b1
x^3: a0.b3+a1.b2+a2.b1+a3.b0 + a2.b3+a3.b2	no change needed
					TOTAL: 4 mul, 4 add to effect corrections. Can we do even better?

replace b1 input with b1+b3, a2 input with a0+a2:

x^0: a0.b0-a1.b3-a2.b2-a3.b1 - a3.b3-a0.b2	need to add a0.b2
x^1: a0.b1+a1.b0-a2.b3-a3.b2			no change needed
x^2: a0.b2+a1.b1+a2.b0-a3.b3 + a1.b3+a0.b0	need to add a3.b3-a0.b0+a2.b2+a3.b1
x^3: a0.b3+a1.b2+a2.b1+a3.b0 + a0.b1+a2.b3	need to add (a2.b3+a3.b2)

Try using sign flips to make mod-Q1 convo look more like a length-4 cyclic:

x^0: a0.b0+a1.b3+a2.b2+a3.b1 + a3.b3		chg sign of a0, then of entire row; a0, b0 frozen
x^1:-a0.b1+a1.b0-a2.b3-a3.b2
x^2:-a0.b2+a1.b1+a2.b0       + (a1.b3+a2.b2+a3.b1)
x^3:-a0.b3+a1.b2+a2.b1+a3.b0 + (a2.b3+a3.b2)


Aside:
Nussbaumer's length-3 cyclic convo algorithm needs 4 mul, 11 add; his length-3 acyclic, OTOH, needs 5 mul and a whopping 20 add.
This is insane, since we can very easily convert a length-3 acyclic into a length-3 cyclic, as follows:

length-3 acyclic:

x^0: a0.b0-a1.b2-a2.b1
x^1: a0.b1+a1.b0-a2.b2
x^2: a0.b2+a1.b1+a2.b0

x^0: a0.b0+a1.b2+a2.b1		flip sign of a0, and then of entire row; a0, b0 frozen
x^1: a0.b1+a1.b0+a2.b2		flip sign of a2 and b1 (which leaves first row unchanged)
x^2: a0.b2+a1.b1+a2.b0		flip sign of entire row, done!

For length-2 acyclic, this fails, since we can only flip an even number of signs this way, and the acyclic has an odd # of minuses.

x^0: a0.b0-a1.b1
x^1: a0.b1+a1.b0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

We then get the final x^0-3 coeffs via

x^0: m0-n3		1 add		-117
x^1: m1+m2		1 add		 -94
x^2: m3+n0+n3		2 add		+165
x^3: n1+n2		1 add		+125

So the total cost for the polypro mod Q1 is 12 MUL, 15 ADD, compared to 5 mul, 15 add for simple length-4 cyclic.
(Multiply count seems high - can reduce to 9 MUL via length-4 cyclic convo with certain input sign changes + clever
combinations of correction terms, but only at the cost of ~10 more adds.)

Lastly, need to combine the 2 modular polynomial product outputs according to CRT.

First, find the two inverse polys, t0 and t1:

t0*Q1 == 1 modulo Q0.
Now:
  1*Q1 = x^4-x^2+1 == 3, so t0 = 1/3.

t1*Q0 == 1 modulo Q1.
Now:
  1*Q0 = x^2+  1,
  x*Q0 = x^3+  x,
x^2*Q0 = x^4+x^2 == 2.x^2-1 modulo Q1, so (x^2-2  )*Q0 == -3 modulo Q1, so t1 = (2-x^2)/3 .

Next, find the s-terms:

s0 := t0*Q1 mod Q = Q1/3              = (1-x^2+x^4)/3 ;
s1 := t1*Q0 mod Q = (2-x^2)*(1+x^2)/3 = (2+x^2-x^4)/3 .

product mod (x^6 + 1) = (q00 + x*q01)*s0 + (q10 + x*q11 + x^2*q12 + x^3*q13)*s1,
where the signature of each q-term times the respective t-polynomial is

q00: s0*3      =  1      -x^2      +x^4
q01: s0*3*x   ==      x        -x^3    +x^5
q10: s1*3      =  2      +x^2      -x^4
q11: s1*3*x   ==    2.x        +x^3    -x^5
q12: s1*3*x^2 ==  1    +2.x^2      +x^4     mod (x^6 + 1)
q13: s1*3*x^3 ==      x      +2.x^3    +x^5 mod (x^6 + 1) .

So, crunching all the numbers and dividing by 3 (the divide can be absorbed into the b-constants,
so is free) gives the following output coefficients, where we define a,b = q00+-q10, c = q10+q12, d,e = q01+-q11, f = q11+q13:
						q00,01 = -93,-3	q10,11,12,13 = -117,-94,+165,+125	a,b,c = -210,+24,+48	d,e,f = -97,+91,+31
x^0:   q00      +2.q10        +q12       = a+c					-162/3 = -54, ok
x^1:         q01      +2.q11        +q13 = d+f					 +66/3 = +22, ok
x^2:  -q00        +q10      +2.q12       = c+c-a				+306/3 = 102, ok
x^3:        -q01        +q11      +2.q13 = f+f-d				+159/3 = +53, ok
x^4:   q00        -q10        +q12       = b+q12 (don't really need b)	+189/3 = +63, ok
x^5:         q01        -q11        +q13 = e+q13 (don't really need e)	+216/3 = +72, ok. It works!

So the total cost of the reconstruction = 14 add.
We can also use the following sequence, which trades 4 multiplies for 6 adds, and thus
might be preferable for a floating-point implementation, where the muls could be done
alongside the adds:

x^4 = q00 - q10 + q12	2 add
x^0 = x^4 + 3.q10	1 add, 1 mul
x^2 = 3.q12 - x^4	1 add, 1 mul
x^5 = q01 - q11 + q13	2 add
x^1 = x^5 + 3.q11	1 add, 1 mul
x^3 = 3.q13 - x^5	1 add, 1 mul


GRAND TOTAL: 41 add, 15 mul (or 35 add, 19 mul), compared to 30 add, 36 mul for the naive scheme.


To do the complete radix-13 DFT, we then need to do as follows:

!   We refer to the terms C1,2,3,4,5,6 (which do not explicitly involving the imaginary constant I)
!   as the "cosine part" of the output, and S1,2,3,4,5,6 (those multiplied by I) as the "sine part."
!                                                                                   opcount for general odd-prime radix R:
!   Form      (x1+-x12),  (x2+-x11),  (x3+-x10),  (x4+-x9),  (x5+-x8),  (x6+-x7) :  0 FMUL, 24 FADD
!   Form X0                                                                      :           2 FADD
!   Form x0+c1*(x1+x12)+c2*(x2+x11)+c3*(x3+x10)+c4*(x4+x9)+c5*(x5+x8)+c6*(x6+x7) :
!   Form x0+c2*(x1+x12)+c4*(x2+x11)+c6*(x3+x10)+c5*(x4+x9)+c3*(x5+x8)+c1*(x6+x7) :
!   Form x0+c3*(x1+x12)+c6*(x2+x11)+c4*(x3+x10)+c1*(x4+x9)+c2*(x5+x8)+c5*(x6+x7) :
!   Form x0+c4*(x1+x12)+c5*(x2+x11)+c1*(x3+x10)+c3*(x4+x9)+c6*(x5+x8)+c2*(x6+x7) :
!   Form x0+c5*(x1+x12)+c3*(x2+x11)+c2*(x3+x10)+c6*(x4+x9)+c1*(x5+x8)+c4*(x6+x7) :
!   Form x0+c6*(x1+x12)+c1*(x2+x11)+c5*(x3+x10)+c2*(x4+x9)+c4*(x5+x8)+c3*(x6+x7) : 16 FMUL, 82 FADD (2 real length-6 cyclics, plus 2 adds for x0+...)
!
!   Form    s1*(x1-x12)+s2*(x2-x11)+s3*(x3-x10)+s4*(x4-x9)+s5*(x5-x8)+s6*(x6-x7) :
!   Form    s2*(x1-x12)+s4*(x2-x11)+s6*(x3-x10)-s5*(x4-x9)-s3*(x5-x8)-s1*(x6-x7) :
!   Form    s3*(x1-x12)+s6*(x2-x11)-s4*(x3-x10)-s1*(x4-x9)+s2*(x5-x8)+s5*(x6-x7) :
!   Form    s4*(x1-x12)-s5*(x2-x11)-s1*(x3-x10)+s3*(x4-x9)-s6*(x5-x8)-s2*(x6-x7) :
!   Form    s5*(x1-x12)-s3*(x2-x11)+s2*(x3-x10)-s6*(x4-x9)-s1*(x5-x8)+s4*(x6-x7) :
!   Form    s6*(x1-x12)-s1*(x2-x11)+s5*(x3-x10)-s2*(x4-x9)+s4*(x5-x8)-s3*(x6-x7) : 30 FMUL, 82 FADD (2 real length-6 acyclics)
!   Form X1,2,3,4,5,6,7,8,9,10,11,12                                             :  0 FMUL, 24 FADD
!
!   Totals :                                                                       46 FMUL, 214 FADD (62 FMUL, 190 FADD for LO_ADD version)
!                                                                     compared to 144 FMUL, 192 FADD for naive scheme,
!                                                                        and just  16 FMUL,  98 FADD for radix-12,
!                                                                        and just  32 FMUL, 160 FADD for radix-14.

UPSHOT: Definitely better than naive radix-13; cut multiply count by more than a factor of 3, with ~10% more adds.
        Still relatively less efficient than radix-12, but not terribly worse than the next-higher composite radix,
        R = 14 (1.55x as many muls per point, 1.43x as many adds per point.)

==========================================================================================================

*/
	int j,j1,j2;
	static int n13,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, first_entry=TRUE;
	static double
					/*cc1 =  .88545602565320989590,	* Real part of exp(i*2*pi/13), the radix-13 fundamental sincos datum	*/
					/*ss1 =  .46472317204376854565,	* Imag part of exp(i*2*pi/13).	*/
					/*cc2 =  .56806474673115580252,	* cos(2u)	*/
					/*ss2 =  .82298386589365639457,	* sin(2u)	*/
					/*cc3 =  .12053668025532305336,	* cos(3u)	*/
					/*ss3 =  .99270887409805399279,	* sin(3u)	*/
					/*cc4 = -.35460488704253562594,	* cos(4u)	*/
					/*ss4 =  .93501624268541482344,	* sin(4u)	*/
					/*cc5 = -.74851074817110109861,	* cos(5u)	*/
					/*ss5 =  .66312265824079520239,	* sin(5u)	*/
					/*cc6 = -.97094181742605202714,	* cos(6u)	*/
					/*ss6 =  .23931566428755776718;	* sin(6u)	*/

		/*sr  =  .10856463647766622055,	* (cc1+cc3+cc4)/6 */
		/*dr  =  .12748655756631447375,	* (cc1-cc3    )/6 */
		/*gr  =  .07919026121630977988,	* (    cc3-cc4)/6 */

		/*tr  = -.19189796981099955387,	* (cc6+cc5+cc2)/6 */
		/*fr  = -.25650109402620130494,	* (cc6    -cc2)/6 */
		/*er  = -.21942924915037615018,	* (    cc5-cc2)/6 */

		bp0  = -.08333333333333333333,	/* sr+tr */
		bp1  =  .30046260628866577442,	/* sr-tr */
		bp20 =  .34691580671669062393,	/* dr-er */
		bp30 = -.09194269158406167643,	/* dr+er */
		bp21 = -.17731083280989152506,	/* fr+gr */
		bp31 = -.33569135524251108482,	/* fr-gr */
		bp2B =  .16960497390679909887,	/* bp20+bp21 = fr+gr+dr-er = [ (cc1-cc4)+(cc6-cc5)]/6 */
		bp3B = -.24374866365844940839,	/* bp31-bp30 = fr-gr-dr-er = [-(cc1-cc4)+(cc6-cc5)]/6 */

		bq00 = -.17413860115213590500,	/* (-ss1-ss3+ss4)/3 */
		bq01 =  .57514072947400312138,	/* ( ss6+ss5+ss2)/3 */
		bq0B =  .40100212832186721638,	/*  bq00+bq01 = (-ss1+ss2-ss3+ss4+ss5+ss6)/3 */
		bq10 = -.46657980490972778969,	/* (-ss1-ss4)/3 */
		bq11 = -.19455606720203287579,	/* ( ss6-ss2)/3 */
		bq12 =  .64257503892782293874,	/* ( ss3+ss4)/3 */
		bq13 =  .05328706921762039739,	/* (-ss5+ss2)/3 */
		bq1A =  .17599523401809514904,	/*  bq10+bq12 = (ss3-ss1)/3 */
		bq1B = -.14126899798441247840;	/*  bq11+bq13 = (ss6-ss5)/3 */

	double y0r,y0i,m0r,m0i,m1r,m1i
	,ar,br,cr,dr,er,fr,gr,sr,tr
	,ai,bi,ci,di,ei,fi,gi,si,ti
	,u0r,u1r,u2r,u3r
	,u0i,u1i,u2i,u3i
	,v0r,v1r,v2r,v3r
	,v0i,v1i,v2i,v3i
	,ac0r,ac1r,ac2r,ac3r,ac4r,ac5r
	,ac0i,ac1i,ac2i,ac3i,ac4i,ac5i
	,as0r,as1r,as2r,as3r,as4r,as5r
	,as0i,as1i,as2i,as3i,as4i,as5i
	,ap0r,ap1r,ap20r,ap21r,ap30r,ap31r
	,ap0i,ap1i,ap20i,ap21i,ap30i,ap31i
	,aq00r,aq01r,aq10r,aq11r,aq12r,aq13r
	,aq00i,aq01i,aq10i,aq11i,aq12i,aq13i
	,p00r,p10r,p20r,p21r,p30r,p31r
	,p00i,p10i,p20i,p21i,p30i,p31i
	,q00r,q01r,q10r,q11r,q12r,q13r
	,q00i,q01i,q10i,q11i,q12i,q13i
	,c1r,c2r,c3r,c4r,c5r,c6r
	,c1i,c2i,c3i,c4i,c5i,c6i
	,s1r,s2r,s3r,s4r,s5r,s6r
	,s1i,s2i,s3i,s4i,s5i,s6i;

	if(!first_entry && (n/13) != n13)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n13=n/13;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n13;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;
		p8 = p7 +p1;
		p9 = p8 +p1;
		p10= p9 +p1;
		p11= p10+p1;
		p12= p11+p1;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 = p9 + ( (p9 >> DAT_BITS) << PAD_BITS );
		p10= p10+ ( (p10>> DAT_BITS) << PAD_BITS );
		p11= p11+ ( (p11>> DAT_BITS) << PAD_BITS );
		p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );
	}

/*...The radix-13 pass is here.	*/

	for(j=0; j < n13; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*    gather the needed data (13 64-bit complex, i.e. 26 64-bit reals) and get needed xi+-xj combos:	*/
		ac0r=a[j1+p1 ]+a[j1+p12];	ac0i=a[j2+p1 ]+a[j2+p12];	as0r=a[j1+p1 ]-a[j1+p12];	as0i=a[j2+p1 ]-a[j2+p12];
		ac1r=a[j1+p2 ]+a[j1+p11];	ac1i=a[j2+p2 ]+a[j2+p11];	as1r=a[j1+p2 ]-a[j1+p11];	as1i=a[j2+p2 ]-a[j2+p11];
		ac4r=a[j1+p3 ]+a[j1+p10];	ac4i=a[j2+p3 ]+a[j2+p10];	as4r=a[j1+p3 ]-a[j1+p10];	as4i=a[j2+p3 ]-a[j2+p10];
		ac2r=a[j1+p4 ]+a[j1+p9 ];	ac2i=a[j2+p4 ]+a[j2+p9 ];	as2r=a[j1+p4 ]-a[j1+p9 ];	as2i=a[j2+p4 ]-a[j2+p9 ];
		ac3r=a[j1+p5 ]+a[j1+p8 ];	ac3i=a[j2+p5 ]+a[j2+p8 ];	as3r=a[j1+p8 ]-a[j1+p5 ];	as3i=a[j2+p8 ]-a[j2+p5 ];
		ac5r=a[j1+p6 ]+a[j1+p7 ];	ac5i=a[j2+p6 ]+a[j2+p7 ];	as5r=a[j1+p6 ]-a[j1+p7 ];	as5i=a[j2+p6 ]-a[j2+p7 ];

	/*** COSINE TERMS ***/

		sr  = ac0r+ac2r+ac4r;		si  = ac0i+ac2i+ac4i;
		dr  = ac0r-ac2r;			di  = ac0i-ac2i;
		gr  =      ac2r-ac4r;		gi  =      ac2i-ac4i;
	/*	ar  = ac0r     -ac4r;		ai  = ac0i     -ac4i;	*/

		tr  = ac1r+ac3r+ac5r;		ti  = ac1i+ac3i+ac5i;
		fr  = ac1r     -ac5r;		fi  = ac1i     -ac5i;
		er  =      ac3r-ac5r;		ei  =      ac3i-ac5i;
	/*	br  = ac1r-ac3r;			bi  = ac1i-ac3i;	*/

		ap0r = sr+tr;				ap0i = si+ti;
		ap1r = sr-tr;				ap1i = si-ti;
		ap20r= dr-er;				ap20i= di-ei;
		ap30r= dr+er;				ap30i= di+ei;
		ap21r= fr+gr;				ap21i= fi+gi;
		ap31r= fr-gr;				ap31i= fi-gi;

		/* polypro(A,B) mod P0 + x0 */
		p00r = ap0r*bp0 + a[j1  ];	p00i = ap0i*bp0 + a[j2];

		/* polypro(A,B) mod P1      */
		p10r = ap1r*bp1;			p10i = ap1i*bp1;

		a[j1  ] += ap0r;			a[j2] += ap0i;

		/* polypro(A,B) mod P2 = x^2-x+1: */
		y0r =(ap20r+ap21r)*bp2B;	y0i =(ap20i+ap21i)*bp2B;/* y0 = (a+b).(c+d) = a.c+a.d+b.c+b.d */
	/*	y0r =(ar+br)*bp2B;			y0i =(ai+bi)*bp2B;	 * use that (ap20+ap21) = (d+g)+(f-e) = (ac0-ac4)+(ac1-ac3) = b+a */
														/* (This is more accurate than above combo, but costs more adds). */
		m0r = ap20r*bp20;			m0i = ap20i*bp20;	/* m0 = a.c	*/
		m1r = ap21r*bp21;			m1i = ap21i*bp21;	/* m1 = b.d	*/
		p20r = m0r - m1r;			p20i = m0i - m1i;	/* z0 = m0-m1	(x^0 output)	*/
		p21r = y0r - m0r;			p21i = y0i - m0i;	/* z1 = y0-m0	(x^1 output).	*/

		/* polypro(A,B) mod P3 = x^2+x+1: */
		y0r =(ap30r-ap31r)*bp3B;	y0i =(ap30i-ap31i)*bp3B;/* y0 = (a-b).(d-c) =-a.c+a.d+b.c-b.d */
	/*	y0r =(ar-br)*bp3B;			y0i =(ai-bi)*bp3B;	 * use that (ap30-ap31) = (d+g)-(f-e) = (ac0-ac4)-(ac1-ac3) = b-a */
														/* (This is more accurate than above combo, but costs more adds). */
		m0r = ap30r*bp30;			m0i = ap30i*bp30;	/* m0 = a.c	*/
		m1r = ap31r*bp31;			m1i = ap31i*bp31;	/* m1 = b.d	*/
		p30r = m0r - m1r;			p30i = m0i - m1i;	/* z0 = m0-m1	(x^0 output)	*/
		p31r = y0r + m0r;			p31i = y0i + m0i;	/* z1 = y0+m0	(x^1 output).	*/

		/* Combine the 4 modular polynomial product outputs according to CRT. */
		ar = p00r+p10r;				ai = p00i+p10i;	/* a,b = (p00+-p10), c,d = (p20+-p30), e,f = (p21+-p31), g = b+d */
		br = p00r-p10r;				bi = p00i-p10i;
		cr = p20r+p30r;				ci = p20i+p30i;
		dr = p20r-p30r;				di = p20i-p30i;
		er = p21r+p31r;				ei = p21i+p31i;
		fr = p21r-p31r;				fi = p21i-p31i;

#if LO_ADD
		c3r = ar-cr+fr;				c3i = ai-ci+fi;		/* x^2 = a - c + f	*/
		c1r = c3r + 3.0*cr;			c1i = c3i + 3.0*ci;	/* x^0 = x^2 + 3.c	*/
		c4r = c3r - 3.0*fr;			c4i = c3i - 3.0*fi;	/* x^4 = x^2 - 3.f	*/
		c2r = br+dr-er;				c2i = bi+di-ei;		/* x^5 = b + d - e	*/
		c6r = c2r + 3.0*er;			c6i = c2i + 3.0*ei;	/* x^1 = x^5 + 3.e	*/
		c5r = c2r - 3.0*dr;			c5i = c2i - 3.0*di;	/* x^3 = x^5 - 3.d	*/
#else
		c1r = ar+cr+cr+fr;			c1i = ai+ci+ci+fi;	/* x^0: a + 2.c +   f				*/
		ar  = ar-cr;				ai  = ai-ci;
		c3r = ar      +fr;			c3i = ai      +fi;	/* x^2: a -   c +   f = (x^0 coeff) - 3.(c)	*/
		c4r = ar   -fr-fr;			c4i = ai   -fi-fi;	/* x^4: a -   c - 2.f = (x^2 coeff) - 3.(f)	*/
		c5r = br-dr-dr-er;			c5i = bi-di-di-ei;	/* x^3: b - 2.d -   e = (x^5 coeff) - 3.(d)	*/
		br  = br+dr;				bi  = bi+di;
		c6r = br   +er+er;			c6i = bi   +ei+ei;	/* x^1: b +   d + 2.e				*/
		c2r = br      -er;			c2i = bi      -ei;	/* x^5: b +   d -   e = (x^1 coeff) - 3.(e)	*/
#endif

	/*** SINE TERMS ***/

		aq00r= as0r-as2r+as4r;		aq00i= as0i-as2i+as4i;
		aq01r= as1r-as3r+as5r;		aq01i= as1i-as3i+as5i;
		aq10r= as0r-as4r;			aq10i= as0i-as4i;
		aq11r= as1r-as5r;			aq11i= as1i-as5i;
		aq12r= as2r+as4r;			aq12i= as2i+as4i;
		aq13r= as3r+as5r;			aq13i= as3i+as5i;

		/* polypro(A,B) mod Q0 = x^2+1: */
		y0r =(aq00r+aq01r)*bq0B;	y0i =(aq00i+aq01i)*bq0B;/* y0 = (a+b).(c+d) = a.c+a.d+b.c+b.d */
		m0r = aq00r*bq00;			m0i = aq00i*bq00;		/* m0 = a.c	*/
		m1r = aq01r*bq01;			m1i = aq01i*bq01;		/* m1 = b.d	*/
		q00r = m0r - m1r;			q00i = m0i - m1i;		/* z0 = m0-m1	(x^0 output)	*/
		q01r = y0r - m0r - m1r;		q01i = y0i - m0i - m1i;	/* z1 = y0-m0-m1(x^1 output).	*/

		/* polypro(A,B) mod Q1 = x^4-x^2+1: */
	/*	y0r = aq10r+aq12r;			y0i = aq10i+aq12i;	 * y0 = (a+c), below, multiply this by (e+g) = a.e+a.g+c.e+c.g */
		y0r = as0r+as2r;			y0i = as0i+as2i;
		m0r = aq10r*bq10;			m0i = aq10i*bq10;	/* m0 = a.e	*/
		m1r = aq12r*bq12;			m1i = aq12i*bq12;	/* m1 = c.g	*/
		u0r = m0r - m1r;			u0i = m0i - m1i;	/* s0 = m0-m1	(x^0 output)	*/
		v0r = y0r*bq1A - m0r;		v0i = y0i*bq1A - m0i;	/* t0 = y0-m0	(x^1 output).	*/

		m0r = aq10r*bq11;			m0i = aq10i*bq11;	/* m0 = a.f	*/
		m1r = aq12r*bq13;			m1i = aq12i*bq13;	/* m1 = c.h	*/
		u1r = m0r - m1r;			u1i = m0i - m1i;	/* s1 = m0-m1	(x^0 output)	*/
		v1r = y0r*bq1B - m0r;		v1i = y0i*bq1B - m0i;	/* t1 = y0-m0	(x^1 output).	*/

	/*	y0r = aq11r+aq13r;			y0i = aq11i+aq13i;	 * y0 = (b+d), below, multiply this by (e+g) = b.e+b.g+d.e+d.g */
		y0r = as1r+as3r;			y0i = as1i+as3i;	/* y0 = (b+d), below, multiply this by (e+g) = b.e+b.g+d.e+d.g */
		m0r = aq11r*bq10;			m0i = aq11i*bq10;	/* m0 = a.e	*/
		m1r = aq13r*bq12;			m1i = aq13i*bq12;	/* m1 = c.g	*/
		u2r = m0r - m1r;			u2i = m0i - m1i;	/* s2 = m0-m1	(x^0 output)	*/
		v2r = y0r*bq1A - m0r;		v2i = y0i*bq1A - m0i;	/* t2 = y0-m0	(x^1 output).	*/

		m0r = aq11r*bq11;			m0i = aq11i*bq11;	/* m0 = b.f	*/
		m1r = aq13r*bq13;			m1i = aq13i*bq13;	/* m1 = d.h	*/
		u3r = m0r - m1r;			u3i = m0i - m1i;	/* s3 = m0-m1	(x^0 output)	*/
		v3r = y0r*bq1B - m0r;		v3i = y0i*bq1B - m0i;	/* t3 = y0-m0	(x^1 output).	*/

		q10r= u0r-v3r;				q10i= u0i-v3i;
		q11r= u1r+u2r;				q11i= u1i+u2i;
		q12r= u3r+v0r+v3r;			q12i= u3i+v0i+v3i;
		q13r= v1r+v2r;				q13i= v1i+v2i;

		/* Combine the 2 modular polynomial product outputs according to CRT. */
#if LO_ADD
		s4r = q00r-q10r+q12r;		s4i = q00i-q10i+q12i;	/* x^4 = q00 - q10 + q12*/
		s1r = s4r + 3.0*q10r;		s1i = s4i + 3.0*q10i;	/* x^0 = x^4 + 3.q10	*/
		s3r = 3.0*q12r - s4r;		s3i = 3.0*q12i - s4i;	/* x^2 = 3.q12 - x^4	*/
		s2r = q01r-q11r+q13r;		s2i = q01i-q11i+q13i;	/* x^5 = q01 - q11 + q13*/
		s6r = s2r + 3.0*q11r;		s6i = s2i + 3.0*q11i;	/* x^1 = x^5 + 3.q11	*/
		s5r = 3.0*q13r - s2r;		s5i = 3.0*q13i - s2i;	/* x^3 = 3.q13 - x^5	*/
#else
		ar = q00r+q10r;		ai = q00i+q10i;	/* a,b = q00+-q10, c = q10+q12, d,e = q01+-q11, f = q11+q13 */
		br = q00r-q10r;		bi = q00i-q10i;
		cr = q10r+q12r;		ci = q10i+q12i;
		dr = q01r+q11r;		di = q01i+q11i;
		er = q01r-q11r;		ei = q01i-q11i;
		fr = q11r+q13r;		fi = q11i+q13i;

		s1r = ar+cr;		s1i = ai+ci;	/* x^0: a + c	*/
		s6r = dr+fr;		s6i = di+fi;	/* x^1: d + f	*/
		s3r = cr+cr-ar;		s3i = ci+ci-ai;	/* x^2: 2.c - a	*/
		s5r = fr+fr-dr;		s5i = fi+fi-di;	/* x^3: 2.f - d	*/
		s4r = br+q12r;		s4i = bi+q12i;	/* x^4: b + q12	*/
		s2r = er+q13r;		s2i = ei+q13i;	/* x^5: e + q13 */
#endif
	/*...Inline multiply of sine parts by +-I into finishing phase, remembering that the acyclic
			 convo algorithm used above returns -s1 and -s5. 	*/

		a[j1+p1 ]=c1r+s1i;	a[j2+p1 ]=c1i-s1r;
		a[j1+p2 ]=c2r-s2i;	a[j2+p2 ]=c2i+s2r;
		a[j1+p3 ]=c3r-s3i;	a[j2+p3 ]=c3i+s3r;
		a[j1+p4 ]=c4r-s4i;	a[j2+p4 ]=c4i+s4r;
		a[j1+p5 ]=c5r+s5i;	a[j2+p5 ]=c5i-s5r;
		a[j1+p6 ]=c6r-s6i;	a[j2+p6 ]=c6i+s6r;
		a[j1+p7 ]=c6r+s6i;	a[j2+p7 ]=c6i-s6r;
		a[j1+p8 ]=c5r-s5i;	a[j2+p8 ]=c5i+s5r;
		a[j1+p9 ]=c4r+s4i;	a[j2+p9 ]=c4i-s4r;
		a[j1+p10]=c3r+s3i;	a[j2+p10]=c3i-s3r;
		a[j1+p11]=c2r+s2i;	a[j2+p11]=c2i-s2r;
		a[j1+p12]=c1r-s1i;	a[j2+p12]=c1i+s1r;
								/* Totals: 46 FMUL, 214 FADD (62 FMUL, 190 FADD for USE_MUL3 version.) */
	}
}

/*************************** Several different versions of the radix-13 DIT to contrast *************************/

void radix13_dit_pass1A(double a[], int n)
{
/*
!
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-13 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2;
	static int n13,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, first_entry=TRUE;
	static double
		bp0  = -.08333333333333333333,
		bp1  =  .30046260628866577442,
		bp20 =  .34691580671669062393,
		bp30 = -.09194269158406167643,
		bp21 = -.17731083280989152506,
		bp31 = -.33569135524251108482,
		bp2B =  .16960497390679909887,
		bp3B = -.24374866365844940839,

		bq00 = -.17413860115213590500,
		bq01 =  .57514072947400312138,
		bq0B =  .40100212832186721638,
		bq10 = -.46657980490972778969,
		bq11 = -.19455606720203287579,
		bq12 =  .64257503892782293874,
		bq13 =  .05328706921762039739,
		bq1A =  .17599523401809514904,
		bq1B = -.14126899798441247840;

	double y0r,y0i,m0r,m0i,m1r,m1i
	,ar,br,cr,dr,er,fr,gr,sr,tr
	,ai,bi,ci,di,ei,fi,gi,si,ti
	,u0r,u1r,u2r,u3r
	,u0i,u1i,u2i,u3i
	,v0r,v1r,v2r,v3r
	,v0i,v1i,v2i,v3i
	,ac0r,ac1r,ac2r,ac3r,ac4r,ac5r
	,ac0i,ac1i,ac2i,ac3i,ac4i,ac5i
	,as0r,as1r,as2r,as3r,as4r,as5r
	,as0i,as1i,as2i,as3i,as4i,as5i
	,ap0r,ap1r,ap20r,ap21r,ap30r,ap31r
	,ap0i,ap1i,ap20i,ap21i,ap30i,ap31i
	,aq00r,aq01r,aq10r,aq11r,aq12r,aq13r
	,aq00i,aq01i,aq10i,aq11i,aq12i,aq13i
	,p00r,p10r,p20r,p21r,p30r,p31r
	,p00i,p10i,p20i,p21i,p30i,p31i
	,q00r,q01r,q10r,q11r,q12r,q13r
	,q00i,q01i,q10i,q11i,q12i,q13i
	,c1r,c2r,c3r,c4r,c5r,c6r
	,c1i,c2i,c3i,c4i,c5i,c6i
	,s1r,s2r,s3r,s4r,s5r,s6r
	,s1i,s2i,s3i,s4i,s5i,s6i;

	if(!first_entry && (n/13) != n13)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n13=n/13;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n13;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;
		p8 = p7 +p1;
		p9 = p8 +p1;
		p10= p9 +p1;
		p11= p10+p1;
		p12= p11+p1;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 = p9 + ( (p9 >> DAT_BITS) << PAD_BITS );
		p10= p10+ ( (p10>> DAT_BITS) << PAD_BITS );
		p11= p11+ ( (p11>> DAT_BITS) << PAD_BITS );
		p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );
	}

/*...The radix-13 pass is here.	*/

	for(j=0; j < n13; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*    gather the needed data (13 64-bit complex, i.e. 26 64-bit reals) and get needed xi+-xj combos:	*/
		ac0r=a[j1+p1 ]+a[j1+p12];	ac0i=a[j2+p1 ]+a[j2+p12];	as0r=a[j1+p1 ]-a[j1+p12];	as0i=a[j2+p1 ]-a[j2+p12];
		ac1r=a[j1+p2 ]+a[j1+p11];	ac1i=a[j2+p2 ]+a[j2+p11];	as1r=a[j1+p2 ]-a[j1+p11];	as1i=a[j2+p2 ]-a[j2+p11];
		ac4r=a[j1+p3 ]+a[j1+p10];	ac4i=a[j2+p3 ]+a[j2+p10];	as4r=a[j1+p3 ]-a[j1+p10];	as4i=a[j2+p3 ]-a[j2+p10];
		ac2r=a[j1+p4 ]+a[j1+p9 ];	ac2i=a[j2+p4 ]+a[j2+p9 ];	as2r=a[j1+p4 ]-a[j1+p9 ];	as2i=a[j2+p4 ]-a[j2+p9 ];
		ac3r=a[j1+p5 ]+a[j1+p8 ];	ac3i=a[j2+p5 ]+a[j2+p8 ];	as3r=a[j1+p8 ]-a[j1+p5 ];	as3i=a[j2+p8 ]-a[j2+p5 ];
		ac5r=a[j1+p6 ]+a[j1+p7 ];	ac5i=a[j2+p6 ]+a[j2+p7 ];	as5r=a[j1+p6 ]-a[j1+p7 ];	as5i=a[j2+p6 ]-a[j2+p7 ];

	/*** COSINE TERMS ***/

		sr  = ac0r+ac2r+ac4r;		si  = ac0i+ac2i+ac4i;
		dr  = ac0r-ac2r;			di  = ac0i-ac2i;
		gr  =      ac2r-ac4r;		gi  =      ac2i-ac4i;
	/*	ar  = ac0r     -ac4r;		ai  = ac0i     -ac4i;	*/

		tr  = ac1r+ac3r+ac5r;		ti  = ac1i+ac3i+ac5i;
		fr  = ac1r     -ac5r;		fi  = ac1i     -ac5i;
		er  =      ac3r-ac5r;		ei  =      ac3i-ac5i;
	/*	br  = ac1r-ac3r;			bi  = ac1i-ac3i;	*/

		ap0r = sr+tr;				ap0i = si+ti;
		ap1r = sr-tr;				ap1i = si-ti;
		ap20r= dr-er;				ap20i= di-ei;
		ap30r= dr+er;				ap30i= di+ei;
		ap21r= fr+gr;				ap21i= fi+gi;
		ap31r= fr-gr;				ap31i= fi-gi;

		/* polypro(A,B) mod P0 + x0 */
		p00r = ap0r*bp0 + a[j1  ];	p00i = ap0i*bp0 + a[j2];

		/* polypro(A,B) mod P1      */
		p10r = ap1r*bp1;			p10i = ap1i*bp1;

		a[j1  ] += ap0r;			a[j2] += ap0i;

		/* polypro(A,B) mod P2 = x^2-x+1: */
		y0r =(ap20r+ap21r)*bp2B;	y0i =(ap20i+ap21i)*bp2B;
	/*	y0r =(ar+br)*bp2B;			y0i =(ai+bi)*bp2B;	*/

		m0r = ap20r*bp20;			m0i = ap20i*bp20;
		m1r = ap21r*bp21;			m1i = ap21i*bp21;
		p20r = m0r - m1r;			p20i = m0i - m1i;
		p21r = y0r - m0r;			p21i = y0i - m0i;

		/* polypro(A,B) mod P3 = x^2+x+1: */
		y0r =(ap30r-ap31r)*bp3B;	y0i =(ap30i-ap31i)*bp3B;
	/*	y0r =(ar-br)*bp3B;			y0i =(ai-bi)*bp3B;	*/
		m0r = ap30r*bp30;			m0i = ap30i*bp30;
		m1r = ap31r*bp31;			m1i = ap31i*bp31;
		p30r = m0r - m1r;			p30i = m0i - m1i;
		p31r = y0r + m0r;			p31i = y0i + m0i;

		/* Combine the 4 modular polynomial product outputs according to CRT. */
		ar = p00r+p10r;				ai = p00i+p10i;
		br = p00r-p10r;				bi = p00i-p10i;
		cr = p20r+p30r;				ci = p20i+p30i;
		dr = p20r-p30r;				di = p20i-p30i;
		er = p21r+p31r;				ei = p21i+p31i;
		fr = p21r-p31r;				fi = p21i-p31i;

		c1r = ar+cr+cr+fr;			c1i = ai+ci+ci+fi;
		ar  = ar-cr;				ai  = ai-ci;
		c3r = ar      +fr;			c3i = ai      +fi;
		c4r = ar   -fr-fr;			c4i = ai   -fi-fi;
		c5r = br-dr-dr-er;			c5i = bi-di-di-ei;
		br  = br+dr;				bi  = bi+di;
		c6r = br   +er+er;			c6i = bi   +ei+ei;
		c2r = br      -er;			c2i = bi      -ei;

	/*** SINE TERMS ***/

		aq00r= as0r-as2r+as4r;		aq00i= as0i-as2i+as4i;
		aq01r= as1r-as3r+as5r;		aq01i= as1i-as3i+as5i;
		aq10r= as0r-as4r;			aq10i= as0i-as4i;
		aq11r= as1r-as5r;			aq11i= as1i-as5i;
		aq12r= as2r+as4r;			aq12i= as2i+as4i;
		aq13r= as3r+as5r;			aq13i= as3i+as5i;

		/* polypro(A,B) mod Q0 = x^2+1: */
		y0r =(aq00r+aq01r)*bq0B;	y0i =(aq00i+aq01i)*bq0B;
		m0r = aq00r*bq00;			m0i = aq00i*bq00;
		m1r = aq01r*bq01;			m1i = aq01i*bq01;
		q00r = m0r - m1r;			q00i = m0i - m1i;
		q01r = y0r - m0r - m1r;		q01i = y0i - m0i - m1i;

		/* polypro(A,B) mod Q1 = x^4-x^2+1: */
		y0r = as0r+as2r;			y0i = as0i+as2i;
		m0r = aq10r*bq10;			m0i = aq10i*bq10;
		m1r = aq12r*bq12;			m1i = aq12i*bq12;
		u0r = m0r - m1r;			u0i = m0i - m1i;
		v0r = y0r*bq1A - m0r;		v0i = y0i*bq1A - m0i;

		m0r = aq10r*bq11;			m0i = aq10i*bq11;
		m1r = aq12r*bq13;			m1i = aq12i*bq13;
		u1r = m0r - m1r;			u1i = m0i - m1i;
		v1r = y0r*bq1B - m0r;		v1i = y0i*bq1B - m0i;

		y0r = as1r+as3r;			y0i = as1i+as3i;
		m0r = aq11r*bq10;			m0i = aq11i*bq10;
		m1r = aq13r*bq12;			m1i = aq13i*bq12;
		u2r = m0r - m1r;			u2i = m0i - m1i;
		v2r = y0r*bq1A - m0r;		v2i = y0i*bq1A - m0i;

		m0r = aq11r*bq11;			m0i = aq11i*bq11;
		m1r = aq13r*bq13;			m1i = aq13i*bq13;
		u3r = m0r - m1r;			u3i = m0i - m1i;
		v3r = y0r*bq1B - m0r;		v3i = y0i*bq1B - m0i;

		q10r= u0r-v3r;				q10i= u0i-v3i;
		q11r= u1r+u2r;				q11i= u1i+u2i;
		q12r= u3r+v0r+v3r;			q12i= u3i+v0i+v3i;
		q13r= v1r+v2r;				q13i= v1i+v2i;

		/* Combine the 2 modular polynomial product outputs according to CRT. */
		ar = q00r+q10r;		ai = q00i+q10i;
		br = q00r-q10r;		bi = q00i-q10i;
		cr = q10r+q12r;		ci = q10i+q12i;
		dr = q01r+q11r;		di = q01i+q11i;
		er = q01r-q11r;		ei = q01i-q11i;
		fr = q11r+q13r;		fi = q11i+q13i;

		s1r = ar+cr;		s1i = ai+ci;
		s6r = dr+fr;		s6i = di+fi;
		s3r = cr+cr-ar;		s3i = ci+ci-ai;
		s5r = fr+fr-dr;		s5i = fi+fi-di;
		s4r = br+q12r;		s4i = bi+q12i;
		s2r = er+q13r;		s2i = ei+q13i;

	/*...Inline multiply of sine parts by +-I into finishing phase, remembering that the acyclic
			 convo algorithm used above returns -s1 and -s5. 	*/

		a[j1+p1 ]=c1r-s1i;	a[j2+p1 ]=c1i+s1r;
		a[j1+p2 ]=c2r+s2i;	a[j2+p2 ]=c2i-s2r;
		a[j1+p3 ]=c3r+s3i;	a[j2+p3 ]=c3i-s3r;
		a[j1+p4 ]=c4r+s4i;	a[j2+p4 ]=c4i-s4r;
		a[j1+p5 ]=c5r-s5i;	a[j2+p5 ]=c5i+s5r;
		a[j1+p6 ]=c6r+s6i;	a[j2+p6 ]=c6i-s6r;
		a[j1+p7 ]=c6r-s6i;	a[j2+p7 ]=c6i+s6r;
		a[j1+p8 ]=c5r+s5i;	a[j2+p8 ]=c5i-s5r;
		a[j1+p9 ]=c4r-s4i;	a[j2+p9 ]=c4i+s4r;
		a[j1+p10]=c3r-s3i;	a[j2+p10]=c3i+s3r;
		a[j1+p11]=c2r-s2i;	a[j2+p11]=c2i+s2r;
		a[j1+p12]=c1r+s1i;	a[j2+p12]=c1i-s1r;
								/* Totals: 46 FMUL, 214 FADD. */
	}
}


/****** TRY TO MINIMIZE THE NUMBER OF TEMPORARIES *********/

void radix13_dit_pass1B(double a[], int n)
{
/*
!
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-13 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2;
	static int n13,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, first_entry=TRUE;
	static double
		bp0  = -.08333333333333333333,
		bp1  =  .30046260628866577442,
		bp20 =  .34691580671669062393,
		bp30 = -.09194269158406167643,
		bp21 = -.17731083280989152506,
		bp31 = -.33569135524251108482,
		bp2B =  .16960497390679909887,
		bp3B = -.24374866365844940839,

		bq00 = -.17413860115213590500,
		bq01 =  .57514072947400312138,
		bq0B =  .40100212832186721638,
		bq10 = -.46657980490972778969,
		bq11 = -.19455606720203287579,
		bq12 =  .64257503892782293874,
		bq13 =  .05328706921762039739,
		bq1A =  .17599523401809514904,
		bq1B = -.14126899798441247840;

	double r1,i1,r2,i2,r3,i3,r4,i4
	,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25;

	if(!first_entry && (n/13) != n13)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n13=n/13;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n13;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;
		p8 = p7 +p1;
		p9 = p8 +p1;
		p10= p9 +p1;
		p11= p10+p1;
		p12= p11+p1;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 = p9 + ( (p9 >> DAT_BITS) << PAD_BITS );
		p10= p10+ ( (p10>> DAT_BITS) << PAD_BITS );
		p11= p11+ ( (p11>> DAT_BITS) << PAD_BITS );
		p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );
	}

/*...The radix-13 pass is here.	*/

	for(j=0; j < n13; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*    gather the needed data (13 64-bit complex, i.e. 26 64-bit reals) and get needed xi+-xj combos:	*/
		t2 =a[j1+p1 ]+a[j1+p12];	t3 =a[j2+p1 ]+a[j2+p12];	t14=a[j1+p1 ]-a[j1+p12];	t15=a[j2+p1 ]-a[j2+p12];
		t4 =a[j1+p2 ]+a[j1+p11];	t5 =a[j2+p2 ]+a[j2+p11];	t16=a[j1+p2 ]-a[j1+p11];	t17=a[j2+p2 ]-a[j2+p11];
		t10=a[j1+p3 ]+a[j1+p10];	t11=a[j2+p3 ]+a[j2+p10];	t22=a[j1+p3 ]-a[j1+p10];	t23=a[j2+p3 ]-a[j2+p10];
		t6 =a[j1+p4 ]+a[j1+p9 ];	t7 =a[j2+p4 ]+a[j2+p9 ];	t18=a[j1+p4 ]-a[j1+p9 ];	t19=a[j2+p4 ]-a[j2+p9 ];
		t8 =a[j1+p5 ]+a[j1+p8 ];	t9 =a[j2+p5 ]+a[j2+p8 ];	t20=a[j1+p8 ]-a[j1+p5 ];	t21=a[j2+p8 ]-a[j2+p5 ];
		t12=a[j1+p6 ]+a[j1+p7 ];	t13=a[j2+p6 ]+a[j2+p7 ];	t24=a[j1+p6 ]-a[j1+p7 ];	t25=a[j2+p6 ]-a[j2+p7 ];

	/*** COSINE TERMS ***/

		t0  = t2 +t6 +t10;		t1  = t3 +t7 +t11;
		t2  = t2 -t6     ;		t3  = t3 -t7     ;
		t6  =     t6 -t10;		t7  =     t7 -t11;

		r1  = t4 +t8 +t12;		i1  = t5 +t9 +t13;
		t4  = t4     -t12;		t5  = t5     -t13;
		t8  =     t8 -t12;		t9  =     t9 -t13;

		t10 = t0+r1;			t11 = t1+i1;
		t0  = t0-r1;			t1  = t1-i1;
		t12 = t2-t8;			t13 = t3-t9;
		t2  = t2+t8;			t3  = t3+t9;
		t8  = t4+t6;			t9  = t5+t7;
		t6  = t4-t6;			t7  = t5-t7;

		/* polypro(A,B) mod P0 + x0 */
		r1 = t10*bp0 + a[j1  ];		i1 = t11*bp0 + a[j2];

		/* polypro(A,B) mod P1      */
		t0 = t0 *bp1;			t1 = t1 *bp1;

		a[j1  ] += t10;			a[j2] += t11;
		/* polypro(A,B) mod P2 = x^2-x+1: */

		t4  =(t12+t8 )*bp2B;	t5  =(t13+t9 )*bp2B;
		t12 = t12*bp20;		t13 = t13*bp20;
		t8  = t8 *bp21;		t9  = t9 *bp21;
		t8  = t12 - t8 ;	t9  = t13 - t9 ;
		t4  = t4  - t12;	t5  = t5  - t13;

		/* polypro(A,B) mod P3 = x^2+x+1: */
		t10 =(t2 -t6 )*bp3B;	t11 =(t3 -t7 )*bp3B;
		t2  = t2 *bp30;		t3  = t3 *bp30;
		t6  = t6 *bp31;		t7  = t7 *bp31;
		t6  = t2  - t6 ;	t7  = t3  - t7 ;
		t10 = t10 + t2 ;	t11 = t11 + t3 ;

		/* Combine the 4 modular polynomial product outputs according to CRT. */
		t2 = r1 +t0;		t3 = i1 +t1;
		t0 = r1 -t0;		t1 = i1 -t1;
		r1 = t8 +t6;		i1 = t9 +t7;
		t8 = t8 -t6;		t9 = t9 -t7;
		t6 = t4 +t10;		t7 = t5 +t11;
		t4 = t4 -t10;		t5 = t5 -t11;

		t10 = t2 -r1;		t11 = t3 -i1;
		t12 = t2 +r1 +r1 +t4;	t13 = t3 +i1 +i1 +t5;	/* C1 */
		t2  = t10        +t4;	t3  = t11        +t5;	/* C3 */
		t4  = t10    -t4 -t4;	t5  = t11    -t5 -t5;	/* C4 */

		r1  = t0 +t8;		i1  = t1 +t9;
		t10 = t0 -t8 -t8 -t6;	t11 = t1 -t9 -t9 -t7;	/* C5 */
		t8  = r1     +t6 +t6;	t9  = i1     +t7 +t7;	/* C6 */
		t6  = r1         -t6;	t7  = i1         -t7;	/* C2 */

	/*** SINE TERMS ***/

		r1 = t14-t18+t22;	i1 = t15-t19+t23;
		t14= t14-t22;		t15= t15-t23;
		t18= t18+t22;		t19= t19+t23;

		t0 = t16-t20+t24;	t1 = t17-t21+t25;
		t16= t16-t24;		t17= t17-t25;
		t20= t20+t24;		t21= t21+t25;

		/* polypro(A,B) mod Q0 = x^2+1: */
		t22 =(r1 +t0 )*bq0B;	t23 =(i1 +t1 )*bq0B;
		r1  = r1 *bq00;		i1  = i1 *bq00;
		t0  = t0 *bq01;		t1  = t1 *bq01;
		t24 = r1  - t0;		t25 = i1  - t1;
		t22 = t22 - r1  - t0;	t23 = t23 - i1  - t1;	/* t0,1 free */

		/* polypro(A,B) mod Q1 = x^4-x^2+1: */
		r1  = t14+t18;		i1  = t15+t19;
		t0  = t14*bq10;		t1  = t15*bq10;
		r2  = t18*bq12;		i2  = t19*bq12;
		r2  = t0  - r2 ;	i2  = t1  - i2 ;
		t0  = r1 *bq1A - t0 ;	t1  = i1 *bq1A - t1 ;

		t14 = t14*bq11;		t15 = t15*bq11;
		t18 = t18*bq13;		t19 = t19*bq13;
		t18 = t14 - t18;	t19 = t15 - t19;
		t14 = r1 *bq1B - t14;	t15 = i1 *bq1B - t15;

		r1  = t16+t20;		i1  = t17+t21;
		r3  = t16*bq10;		i3  = t17*bq10;
		r4  = t20*bq12;		i4  = t21*bq12;
		r4  = r3  - r4;		i4  = i3  - i4;
		r3  = r1 *bq1A - r3 ;	i3  = i1 *bq1A - i3 ;

		t16 = t16*bq11;		t17 = t17*bq11;
		t20 = t20*bq13;		t21 = t21*bq13;
		t20 = t16 - t20;	t21 = t17 - t21;
		t16 = r1 *bq1B - t16;	t17 = i1 *bq1B - t17;

		t20 = t20+t0 +t16;	t21 = t21+t1 +t17;
		t16 = r2 -t16;		t17 = i2 -t17;
		t18 = t18+r4;		t19 = t19+i4;
		t14 = t14+r3 ;		t15 = t15+i3 ;

		/* Combine the 2 modular polynomial product outputs according to CRT. */
		r1  = t24+t16;		i1  = t25+t17;
		t24 = t24-t16;		t25 = t25-t17;
		t16 = t16+t20;		t17 = t17+t21;
		t0  = t22+t18;		t1  = t23+t19;
		t22 = t22-t18;		t23 = t23-t19;
		t18 = t18+t14;		t19 = t19+t15;

		t24 = t24+t20;		t25 = t25+t21;		/* S4 */
		t22 = t22+t14;		t23 = t23+t15;		/* S2 */
		t14 = r1 +t16;		t15 = i1 +t17;		/* S1 */
		t16 = t16+t16-r1 ;	t17 = t17+t17-i1 ;	/* S3 */
		t20 = t0 +t18;		t21 = t1 +t19;		/* S6 */
		t18 = t18+t18-t0 ;	t19 = t19+t19-t1 ;	/* S5 */

	/*...Inline multiply of sine parts by +-I into finishing phase, remembering that the acyclic
			 convo algorithm used above returns -s1 and -s5. 	*/

		a[j1+p1 ]=t12-t15;	a[j2+p1 ]=t13+t14;	/* C1 - S1 */
		a[j1+p2 ]=t6 +t23;	a[j2+p2 ]=t7 -t22;	/* C2 - S2 */
		a[j1+p3 ]=t2 +t17;	a[j2+p3 ]=t3 -t16;	/* C3 - S3 */
		a[j1+p4 ]=t4 +t25;	a[j2+p4 ]=t5 -t24;	/* C4 - S4 */
		a[j1+p5 ]=t10-t19;	a[j2+p5 ]=t11+t18;	/* C5 - S5 */
		a[j1+p6 ]=t8 +t21;	a[j2+p6 ]=t9 -t20;	/* C6 - S6 */
		a[j1+p7 ]=t8 -t21;	a[j2+p7 ]=t9 +t20;	/* C6 + S6 */
		a[j1+p8 ]=t10+t19;	a[j2+p8 ]=t11-t18;	/* C5 + S5 */
		a[j1+p9 ]=t4 -t25;	a[j2+p9 ]=t5 +t24;	/* C4 + S4 */
		a[j1+p10]=t2 -t17;	a[j2+p10]=t3 +t16;	/* C3 + S3 */
		a[j1+p11]=t6 -t23;	a[j2+p11]=t7 +t22;	/* C2 + S2 */
		a[j1+p12]=t12+t15;	a[j2+p12]=t13-t14;	/* C1 + S1 */
								/* Totals: 46 FMUL, 214 FADD. */
	}
}


/****** TRY TO MINIMIZE THE NUMBER OF ADDS *********/

void radix13_dit_pass1(double a[], int n)
{
/*
!
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-13 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2;
	static int n13,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, first_entry=TRUE;
	static double
		bp0  = -.08333333333333333333,
		bp1  =  .30046260628866577442,
		bp20 =  .34691580671669062393,
		bp30 = -.09194269158406167643,
		bp21 = -.17731083280989152506,
		bp31 = -.33569135524251108482,
		bp2B =  .16960497390679909887,
		bp3B = -.24374866365844940839,

		bq00 = -.17413860115213590500,
		bq01 =  .57514072947400312138,
		bq0B =  .40100212832186721638,
		bq10 = -.46657980490972778969,
		bq11 = -.19455606720203287579,
		bq12 =  .64257503892782293874,
		bq13 =  .05328706921762039739,
		bq1A =  .17599523401809514904,
		bq1B = -.14126899798441247840;

	double r1,i1,r2,i2,r3,i3,r4,i4
	,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25;

	if(!first_entry && (n/13) != n13)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n13=n/13;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n13;
		p2 = p1 +p1;
		p3 = p2 +p1;
		p4 = p3 +p1;
		p5 = p4 +p1;
		p6 = p5 +p1;
		p7 = p6 +p1;
		p8 = p7 +p1;
		p9 = p8 +p1;
		p10= p9 +p1;
		p11= p10+p1;
		p12= p11+p1;

		p1 = p1 + ( (p1 >> DAT_BITS) << PAD_BITS );
		p2 = p2 + ( (p2 >> DAT_BITS) << PAD_BITS );
		p3 = p3 + ( (p3 >> DAT_BITS) << PAD_BITS );
		p4 = p4 + ( (p4 >> DAT_BITS) << PAD_BITS );
		p5 = p5 + ( (p5 >> DAT_BITS) << PAD_BITS );
		p6 = p6 + ( (p6 >> DAT_BITS) << PAD_BITS );
		p7 = p7 + ( (p7 >> DAT_BITS) << PAD_BITS );
		p8 = p8 + ( (p8 >> DAT_BITS) << PAD_BITS );
		p9 = p9 + ( (p9 >> DAT_BITS) << PAD_BITS );
		p10= p10+ ( (p10>> DAT_BITS) << PAD_BITS );
		p11= p11+ ( (p11>> DAT_BITS) << PAD_BITS );
		p12= p12+ ( (p12>> DAT_BITS) << PAD_BITS );
	}

/*...The radix-13 pass is here.	*/

	for(j=0; j < n13; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*    gather the needed data (13 64-bit complex, i.e. 26 64-bit reals) and get needed xi+-xj combos:	*/
		t2 =a[j1+p1 ]+a[j1+p12];	t3 =a[j2+p1 ]+a[j2+p12];	t14=a[j1+p1 ]-a[j1+p12];	t15=a[j2+p1 ]-a[j2+p12];
		t4 =a[j1+p2 ]+a[j1+p11];	t5 =a[j2+p2 ]+a[j2+p11];	t16=a[j1+p2 ]-a[j1+p11];	t17=a[j2+p2 ]-a[j2+p11];
		t10=a[j1+p3 ]+a[j1+p10];	t11=a[j2+p3 ]+a[j2+p10];	t22=a[j1+p3 ]-a[j1+p10];	t23=a[j2+p3 ]-a[j2+p10];
		t6 =a[j1+p4 ]+a[j1+p9 ];	t7 =a[j2+p4 ]+a[j2+p9 ];	t18=a[j1+p4 ]-a[j1+p9 ];	t19=a[j2+p4 ]-a[j2+p9 ];
		t8 =a[j1+p5 ]+a[j1+p8 ];	t9 =a[j2+p5 ]+a[j2+p8 ];	t20=a[j1+p8 ]-a[j1+p5 ];	t21=a[j2+p8 ]-a[j2+p5 ];
		t12=a[j1+p6 ]+a[j1+p7 ];	t13=a[j2+p6 ]+a[j2+p7 ];	t24=a[j1+p6 ]-a[j1+p7 ];	t25=a[j2+p6 ]-a[j2+p7 ];

	/*** COSINE TERMS ***/
		t0  = t2 +t6 +t10;		t1  = t3 +t7 +t11;
		t2  = t2 -t6     ;		t3  = t3 -t7     ;
		t6  =     t6 -t10;		t7  =     t7 -t11;

		r1  = t4 +t8 +t12;		i1  = t5 +t9 +t13;
		t4  = t4     -t12;		t5  = t5     -t13;
		t8  =     t8 -t12;		t9  =     t9 -t13;

		t10 = t0+r1;			t11 = t1+i1;
		t0  = t0-r1;			t1  = t1-i1;
		t12 = t2-t8;			t13 = t3-t9;
		t2  = t2+t8;			t3  = t3+t9;
		t8  = t4+t6;			t9  = t5+t7;
		t6  = t4-t6;			t7  = t5-t7;

		/* polypro(A,B) mod P0 + x0 */
		r1 = t10*bp0 + a[j1  ];		i1 = t11*bp0 + a[j2];

		/* polypro(A,B) mod P1      */
		t0 = t0 *bp1;			t1 = t1 *bp1;

		a[j1  ] += t10;			a[j2] += t11;
		/* polypro(A,B) mod P2 = x^2-x+1: */

		t4  =(t12+t8 )*bp2B;	t5  =(t13+t9 )*bp2B;
		t12 = t12*bp20;		t13 = t13*bp20;
		t8  = t8 *bp21;		t9  = t9 *bp21;
		t8  = t12 - t8 ;	t9  = t13 - t9 ;
		t4  = t4  - t12;	t5  = t5  - t13;

		/* polypro(A,B) mod P3 = x^2+x+1: */
		t10 =(t2 -t6 )*bp3B;	t11 =(t3 -t7 )*bp3B;
		t2  = t2 *bp30;		t3  = t3 *bp30;
		t6  = t6 *bp31;		t7  = t7 *bp31;
		t6  = t2  - t6 ;	t7  = t3  - t7 ;
		t10 = t10 + t2 ;	t11 = t11 + t3 ;

		/* Combine the 4 modular polynomial product outputs according to CRT. */
		t2 = r1 +t0;		t3 = i1 +t1;
		t0 = r1 -t0;		t1 = i1 -t1;
		r1 = t8 +t6;		i1 = t9 +t7;
		t8 = t8 -t6;		t9 = t9 -t7;
		t6 = t4 +t10;		t7 = t5 +t11;
		t4 = t4 -t10;		t5 = t5 -t11;

		t2 = t2 - r1 + t4;	t3 = t3 -i1 + t5;	/* C3 = x^2 = a - c + f		2 add        */
		t12= t2 + 3*r1;		t13= t3 + 3*i1;		/* C1 = x^0 = x^2 + 3.c		1 add, 1 mul */
		t4 = t2 - 3*t4;		t5 = t3 - 3*t5;		/* C4 = x^4 = x^2 - 3.f		1 add, 1 mul */

		t10= t0 + t8 - t6;	t11= t1 + t9 - t7;	/* C2 = x^5 = b + d - e		2 add        */
		t6 = t10+ 3*t6;		t7 = t11+ 3*t7;		/* C6 = x^1 = x^5 + 3.e		1 add, 1 mul */
		t8 = t10- 3*t8;		t9 = t11- 3*t9;		/* C5 = x^3 = x^5 - 3.d		1 add, 1 mul */

	/*** SINE TERMS ***/

		r1 = t14-t18+t22;	i1 = t15-t19+t23;
		t14= t14-t22;		t15= t15-t23;
		t18= t18+t22;		t19= t19+t23;

		t0 = t16-t20+t24;	t1 = t17-t21+t25;
		t16= t16-t24;		t17= t17-t25;
		t20= t20+t24;		t21= t21+t25;

		/* polypro(A,B) mod Q0 = x^2+1: */
		t22 =(r1 +t0 )*bq0B;	t23 =(i1 +t1 )*bq0B;
		r1  = r1 *bq00;		i1  = i1 *bq00;
		t0  = t0 *bq01;		t1  = t1 *bq01;
		t24 = r1  - t0;		t25 = i1  - t1;		/* q00 */
		t22 = t22 - r1  - t0;	t23 = t23 - i1  - t1;	/* q01 */

		/* polypro(A,B) mod Q1 = x^4-x^2+1: */
		r1  = t14+t18;		i1  = t15+t19;
		t0  = t14*bq10;		t1  = t15*bq10;
		r2  = t18*bq12;		i2  = t19*bq12;
		r2  = t0  - r2 ;	i2  = t1  - i2 ;
		t0  = r1 *bq1A - t0 ;	t1  = i1 *bq1A - t1 ;

		t14 = t14*bq11;		t15 = t15*bq11;
		t18 = t18*bq13;		t19 = t19*bq13;
		t18 = t14 - t18;	t19 = t15 - t19;
		t14 = r1 *bq1B - t14;	t15 = i1 *bq1B - t15;

		r1  = t16+t20;		i1  = t17+t21;
		r3  = t16*bq10;		i3  = t17*bq10;
		r4  = t20*bq12;		i4  = t21*bq12;
		r4  = r3  - r4;		i4  = i3  - i4;
		r3  = r1 *bq1A - r3 ;	i3  = i1 *bq1A - i3 ;

		t16 = t16*bq11;		t17 = t17*bq11;
		t20 = t20*bq13;		t21 = t21*bq13;
		t20 = t16 - t20;	t21 = t17 - t21;
		t16 = r1 *bq1B - t16;	t17 = i1 *bq1B - t17;

		t20 = t20+t0 +t16;	t21 = t21+t1 +t17;	/* q12 */
		t16 = r2 -t16;		t17 = i2 -t17;		/* q10 */
		t18 = t18+r4;		t19 = t19+i4;		/* q11 */
		t14 = t14+r3;		t15 = t15+i3;		/* q13 */

		/* Combine the 2 modular polynomial product outputs according to CRT. */

		t24 = t24 - t16 + t20;	t25 = t25 - t17 + t21;	/* S4 = x^4 = q00 - q10 + q12	2 add        */
		t16 = t24 + 3*t16;	t17 = t25 + 3*t17;	/* S1 = x^0 = x^4 + 3.q10	1 add, 1 mul */
		t20 = 3*t20 - t24;	t21 = 3*t21 - t25;	/* S3 = x^2 = 3.q12 - x^4	1 add, 1 mul */

		t22 = t22 - t18 + t14;	t23 = t23 - t19 + t15;	/* S2 = x^5 = q01 - q11 + q13	2 add        */
		t18 = t22 + 3*t18;	t19 = t23 + 3*t19;	/* S6 = x^1 = x^5 + 3.q11	1 add, 1 mul */
		t14 = 3*t14 - t22;	t15 = 3*t15 - t23;	/* S5 = x^3 = 3.q13 - x^5	1 add, 1 mul */

	/*...Inline multiply of sine parts by +-I into finishing phase, remembering that the acyclic
			 convo algorithm used above returns -s1 and -s5. 	*/

		a[j1+p1 ]=t12-t17;	a[j2+p1 ]=t13+t16;	/* C1 - S1 */
		a[j1+p2 ]=t10+t23;	a[j2+p2 ]=t11-t22;	/* C2 - S2 */
		a[j1+p3 ]=t2 +t21;	a[j2+p3 ]=t3 -t20;	/* C3 - S3 */
		a[j1+p4 ]=t4 +t25;	a[j2+p4 ]=t5 -t24;	/* C4 - S4 */
		a[j1+p5 ]=t8 -t15;	a[j2+p5 ]=t9 +t14;	/* C5 - S5 */
		a[j1+p6 ]=t6 +t19;	a[j2+p6 ]=t7 -t18;	/* C6 - S6 */
		a[j1+p7 ]=t6 -t19;	a[j2+p7 ]=t7 +t18;	/* C6 + S6 */
		a[j1+p8 ]=t8 +t15;	a[j2+p8 ]=t9 -t14;	/* C5 + S5 */
		a[j1+p9 ]=t4 -t25;	a[j2+p9 ]=t5 +t24;	/* C4 + S4 */
		a[j1+p10]=t2 -t21;	a[j2+p10]=t3 +t20;	/* C3 + S3 */
		a[j1+p11]=t10-t23;	a[j2+p11]=t11+t22;	/* C2 + S2 */
		a[j1+p12]=t12+t17;	a[j2+p12]=t13-t16;	/* C1 + S1 */
								/* Totals: 62 FMUL, 190 FADD. */
	}
}

