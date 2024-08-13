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

#include "Mlucas.h"

/***************/

int radix15_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-15 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-15 complex DIF pass on the data in
!   the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix7/8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n15,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14
		,i,j,j1,j2,jstart,jhi,full_pass,k1,k2,k,khi,l,outer;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14, nsave = 0;
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/

					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292,	/* [sin(u)-sin(2u)] */
					radix_inv, n2inv;
	double rt,it
	  #if 0
		,t0,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
	  #endif
		,t1,t2
		,aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r,aj1p6r,aj1p7r,aj1p8r,aj1p9r,aj1p10r,aj1p11r,aj1p12r,aj1p13r,aj1p14r
		,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i,aj1p6i,aj1p7i,aj1p8i,aj1p9i,aj1p10i,aj1p11i,aj1p12i,aj1p13i,aj1p14i
		,cy_r0,cy_r1,cy_r2,cy_r3,cy_r4,cy_r5,cy_r6,cy_r7,cy_r8,cy_r9,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14
		,cy_i0,cy_i1,cy_i2,cy_i3,cy_i4,cy_i5,cy_i6,cy_i7,cy_i8,cy_i9,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14
		,temp,frac,scale;
	double maxerr = 0.0;
#if 0
#if PFETCH
	double *add0, *addr;
#endif
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	int ii0,ii1,ii2,ii3,ii4,ii5,ii6,ii7,ii8,ii9,ii10,ii11,ii12,ii13,ii14;	/* indices into weights arrays (mod NWT) */
	double wt_re,wt_im, wi_re,wi_im;					/* Fermat/LOACC weights stuff */

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=ii0=ii1=ii2=ii3=ii4=ii5=ii6=ii7=ii8=ii9=ii10=ii11=ii12=ii13=ii14=-1;

	if(RES_SHIFT) {
		WARN(HERE, "CY routines with radix < 16 do not support shifted residues!", "", 1);
		return(ERR_ASSERT);
	}

	// Jan 2018: To support PRP-testing, read the LR-modpow-scalar-multiply-needed bit for the current iteration from the global array:
	double prp_mult = 1.0;
	if((TEST_TYPE & 0xfffffffe) == TEST_TYPE_PRP) {	// Mask off low bit to lump together PRP and PRP-C tests
		i = (iter-1) % ITERS_BETWEEN_CHECKPOINTS;	// Bit we need to read...iter-counter is unit-offset w.r.to iter-interval, hence the -1
		if((BASE_MULTIPLIER_BITS[i>>6] >> (i&63)) & 1)
			prp_mult = PRP_BASE;
	}

/*...change n15 and n_div_wt to non-static to work around a gcc compiler bug. */
	n15   = n/15;
	n_div_nwt = n15 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n15)
	{
		sprintf(cbuf,"ERROR: iter = %10d; NWT_BITS does not divide N/15 in radix15_ditN_cy_dif1.\n",iter);
		mlucas_fprint(cbuf,INTERACT);
		err=ERR_CARRY;
		return(err);
	}

	if(p != psave || n != nsave) {	/* Exponent or array length change triggers re-init */
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		psave = p;	nsave = n;
		first_entry=FALSE;
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)15));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

		/*   constant index offsets for array load/stores are here.	*/

		p1 = n15;
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
		p13= p12+p1;
		p14= p13+p1;

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
		p13= p13+ ( (p13>> DAT_BITS) << PAD_BITS );
		p14= p14+ ( (p14>> DAT_BITS) << PAD_BITS );

		/* For Fermat-mod, since 'adjacent' words are actually stride-2 separated
		in terms of the floating residue array, block boundaries have half the i-index
		(e.g. as in sw*i and bw*i) value they do in the Mersenne-mod case:
		*/
		khi = n15 >> (MODULUS_TYPE != MODULUS_TYPE_MERSENNE);	// loop up to n15 for Mers, n15/2 for Fermat-mod
		bjmodnini=0;
		for(j=0; j < n15; j++) {
			bjmodnini -= sw; bjmodnini += ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
	}

	/* init carries	*/
	cy_r0 = 0;	cy_i0 = 0;
	cy_r1 = 0;	cy_i1 = 0;
	cy_r2 = 0;	cy_i2 = 0;
	cy_r3 = 0;	cy_i3 = 0;
	cy_r4 = 0;	cy_i4 = 0;
	cy_r5 = 0;	cy_i5 = 0;
	cy_r6 = 0;	cy_i6 = 0;
	cy_r7 = 0;	cy_i7 = 0;
	cy_r8 = 0;	cy_i8 = 0;
	cy_r9 = 0;	cy_i9 = 0;
	cy_r10= 0;	cy_i10= 0;
	cy_r11= 0;	cy_i11= 0;
	cy_r12= 0;	cy_i12= 0;
	cy_r13= 0;	cy_i13= 0;
	cy_r14= 0;	cy_i14= 0;

	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		cy_r0 = -2;
	}

	*fracmax=0;	/* init max. fractional error	*/
	full_pass = 1;	/* set = 1 for normal carry pass, = 0 for wrapper pass	*/
	scale = n2inv;	// init inverse-weight scale factor = 2/n for normal carry pass, 1 for wrapper pass

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart = 0;
		jhi = jstart+nwt-1;
		khi = n_div_nwt;
	} else {
		jstart = 0;
		jhi = n_div_nwt;
		khi = 1;
	}

for(outer=0; outer <= 1; outer++)
{
	// i = Index into the BASE and BASEINV arrays:
	i = (bw > 0);	/* If bw > 0 (i.e. n does not divide p), lowest-order digit always a bigword:	*/

	bjmodn0 = 0;	l = bjmodnini-n;
	bjmodn1 = bjmodnini;
	bjmodn2 = bjmodn1 +l; bjmodn2  += ( (-(int)((uint32)bjmodn2 >> 31)) & n);
	bjmodn3 = bjmodn2 +l; bjmodn3  += ( (-(int)((uint32)bjmodn3 >> 31)) & n);
	bjmodn4 = bjmodn3 +l; bjmodn4  += ( (-(int)((uint32)bjmodn4 >> 31)) & n);
	bjmodn5 = bjmodn4 +l; bjmodn5  += ( (-(int)((uint32)bjmodn5 >> 31)) & n);
	bjmodn6 = bjmodn5 +l; bjmodn6  += ( (-(int)((uint32)bjmodn6 >> 31)) & n);
	bjmodn7 = bjmodn6 +l; bjmodn7  += ( (-(int)((uint32)bjmodn7 >> 31)) & n);
	bjmodn8 = bjmodn7 +l; bjmodn8  += ( (-(int)((uint32)bjmodn8 >> 31)) & n);
	bjmodn9 = bjmodn8 +l; bjmodn9  += ( (-(int)((uint32)bjmodn9 >> 31)) & n);
	bjmodn10= bjmodn9 +l; bjmodn10 += ( (-(int)((uint32)bjmodn10>> 31)) & n);
	bjmodn11= bjmodn10+l; bjmodn11 += ( (-(int)((uint32)bjmodn11>> 31)) & n);
	bjmodn12= bjmodn11+l; bjmodn12 += ( (-(int)((uint32)bjmodn12>> 31)) & n);
	bjmodn13= bjmodn12+l; bjmodn13 += ( (-(int)((uint32)bjmodn13>> 31)) & n);
	bjmodn14= bjmodn13+l; bjmodn14 += ( (-(int)((uint32)bjmodn14>> 31)) & n);

	/* For Fermat-mod, IBDWT access patterns repeat with period NWT = {odd part of radix0},
	so for even radix0's only really need that many bjmodn and ii's, but that would require
	specialized carry macros tha`t don't update ii and bjmodn - not worth the trouble.
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		col=0;
		co2=(n >> nwt_bits)-1+15;
		co3=co2-15;		/* At the start of each new j-loop, co3=co2-radix(1)	*/
	}
	else
	{
		/* indices into IBDWT weights arrays (mod NWT) is here: */
		ii0 = 0;
		ii1 = (SW_DIV_N*n15/2) % nwt;
		ii2 = (ii1 + ii1) % nwt;
		ii3 = (ii2 + ii1) % nwt;
		ii4 = (ii3 + ii1) % nwt;
		ii5 = (ii4 + ii1) % nwt;
		ii6 = (ii5 + ii1) % nwt;
		ii7 = (ii6 + ii1) % nwt;
		ii8 = (ii7 + ii1) % nwt;
		ii9 = (ii8 + ii1) % nwt;
		ii10= (ii9 + ii1) % nwt;
		ii11= (ii10+ ii1) % nwt;
		ii12= (ii11+ ii1) % nwt;
		ii13= (ii12+ ii1) % nwt;
		ii14= (ii13+ ii1) % nwt;

		/* Start this value off at N in Fermat-mod case, so (bjmodn >= sw) check in
		fermat_carry_norm_errcheck (cf. carry.h) yields a bigword (i == 1) for j= 0:
		*/
		bjmodn0 = n;
	}

	/* Move this cleanup-pass-specific khi setting here, since need regular-pass khi value for above inits: */
	if(!full_pass)
	{
		khi = 1;
	}

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

	// Gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do five radix-3 transforms...
		#if 1
		RADIX_15_DIT(
			a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],
			aj1p0r,aj1p0i,aj1p1r,aj1p1i,aj1p2r,aj1p2i,aj1p3r,aj1p3i,aj1p4r,aj1p4i,aj1p5r,aj1p5i,aj1p6r,aj1p6i,aj1p7r,aj1p7i,aj1p8r,aj1p8i,aj1p9r,aj1p9i,aj1p10r,aj1p10i,aj1p11r,aj1p11i,aj1p12r,aj1p12i,aj1p13r,aj1p13i,aj1p14r,aj1p14i
		)
		#else
		/* ...Block 1:	*/
			t0 =a[j1    ];	t1 =a[j2    ];
			t2 =a[j1+p2 ];	t3 =a[j2+p2 ];
			rt =a[j1+p1 ];	it =a[j2+p1 ];
			t4 =t2 -rt;				t5 =t3 -it;
			t2 =t2 +rt;				t3 =t3 +it;
			t0 =t0+t2;				t1 =t1+t3;
			t2 =t0+c3m1*t2;			t3 =t1+c3m1*t3;
			rt =s*t4;				it =s*t5;
			t4 =t2-it;				t5 =t3+rt;
			t2 =t2+it;				t3 =t3-rt;
		/* ...Block 2:	*/
			t6 =a[j1+p8 ];	t7 =a[j2+p8 ];
			t8 =a[j1+p7 ];	t9 =a[j2+p7 ];
			rt =a[j1+p6 ];	it =a[j2+p6 ];
			t10=t8 -rt;				t11=t9 -it;
			t8 =t8 +rt;				t9 =t9 +it;
			t6 =t6+t8;				t7 =t7+t9;
			t8 =t6+c3m1*t8;			t9 =t7+c3m1*t9;
			rt =s*t10;				it =s*t11;
			t10=t8-it;				t11=t9+rt;
			t8 =t8+it;				t9 =t9-rt;
		/* ...Block 3:	*/
			t12=a[j1+p13];	t13=a[j2+p13];
			t14=a[j1+p12];	t15=a[j2+p12];
			rt =a[j1+p14];	it =a[j2+p14];
			t16=t14 -rt;			t17=t15 -it;
			t14=t14 +rt;			t15=t15 +it;
			t12=t12+t14;			t13=t13+t15;
			t14=t12+c3m1*t14;			t15=t13+c3m1*t15;
			rt =s*t16;				it =s*t17;
			t16=t14-it;				t17=t15+rt;
			t14=t14+it;				t15=t15-rt;
		/* ...Block 4:	*/
			t18=a[j1+p4 ];	t19=a[j2+p4 ];
			t20=a[j1+p3 ];	t21=a[j2+p3 ];
			rt =a[j1+p5 ];	it =a[j2+p5 ];
			t22=t20 -rt;			t23=t21 -it;
			t20=t20 +rt;			t21=t21 +it;
			t18=t18+t20;			t19=t19+t21;
			t20=t18+c3m1*t20;			t21=t19+c3m1*t21;
			rt =s*t22;				it =s*t23;
			t22=t20-it;				t23=t21+rt;
			t20=t20+it;				t21=t21-rt;
		/* ...Block 5:	*/
			t24=a[j1+p9 ];	t25=a[j2+p9 ];
			t26=a[j1+p11];	t27=a[j2+p11];
			rt =a[j1+p10];	it =a[j2+p10];
			t28=t26 -rt;			t29=t27 -it;
			t26=t26 +rt;			t27=t27 +it;
			t24=t24+t26;			t25=t25+t27;
			t26=t24+c3m1*t26;			t27=t25+c3m1*t27;
			rt =s*t28;				it =s*t29;
			t28=t26-it;				t29=t27+rt;
			t26=t26+it;				t27=t27-rt;

	// ...and now do three radix-5 transforms:
		/* ...Block 1:	*/
			rt =t24;			it =t25;
			t24=t6-rt;			t25=t7-it;
			t6 =t6+rt;			t7 =t7+it;
			rt =t18;			it =t19;
			t18=t12-rt;			t19=t13-it;
			t12=t12+rt;			t13=t13+it;

			rt = t6+t12;		it = t7+t13;
			t0 = t0+rt;			t1 = t1+it;
			rt = t0+cn1*rt;		it = t1+cn1*it;
			t12= cn2*(t6-t12);		t13= cn2*(t7-t13);
			t6 = rt+t12;		t7 = it+t13;
			t12= rt-t12;		t13= it-t13;
			rt = ss3*(t18-t24);		it = ss3*(t19-t25);
			t18= rt-sn1*t18;		t19= it-sn1*t19;
			t24= rt+sn2*t24;		t25= it+sn2*t25;

			aj1p0r =t0;			aj1p0i =t1;
			aj1p9r =t6-t19;		aj1p9i =t7+t18;
			aj1p3r =t12-t25;		aj1p3i =t13+t24;
			aj1p12r=t12+t25;		aj1p12i=t13-t24;
			aj1p6r =t6+t19;		aj1p6i =t7-t18;

		/* ...Block 2:	*/
			rt =t26;			it =t27;
			t26=t8-rt;			t27=t9-it;
			t8 =t8+rt;			t9 =t9+it;
			rt =t20;			it =t21;
			t20=t14-rt;			t21=t15-it;
			t14=t14+rt;			t15=t15+it;

			rt = t8+t14;		it = t9+t15;
			t2 = t2+rt;			t3 = t3+it;
			rt = t2+cn1*rt;		it = t3+cn1*it;
			t14= cn2*(t8-t14);		t15= cn2*(t9-t15);
			t8 = rt+t14;		t9 = it+t15;
			t14= rt-t14;		t15= it-t15;
			rt = ss3*(t20-t26);		it = ss3*(t21-t27);
			t20= rt-sn1*t20;		t21= it-sn1*t21;
			t26= rt+sn2*t26;		t27= it+sn2*t27;

			aj1p5r =t2;			aj1p5i =t3;
			aj1p14r=t8-t21;		aj1p14i=t9+t20;
			aj1p8r =t14-t27;		aj1p8i =t15+t26;
			aj1p2r =t14+t27;		aj1p2i =t15-t26;
			aj1p11r=t8+t21;		aj1p11i=t9-t20;

		/* ...Block 3:	*/
			rt =t28;			it =t29;
			t28=t10-rt;			t29=t11-it;
			t10=t10+rt;			t11=t11+it;
			rt =t22;			it =t23;
			t22=t16-rt;			t23=t17-it;
			t16=t16+rt;			t17=t17+it;

			rt = t10+t16;		it = t11+t17;
			t4 = t4+rt;			t5 = t5+it;
			rt = t4+cn1*rt;		it = t5+cn1*it;
			t16= cn2*(t10-t16);		t17= cn2*(t11-t17);
			t10= rt+t16;		t11= it+t17;
			t16= rt-t16;		t17= it-t17;
			rt = ss3*(t22-t28);		it = ss3*(t23-t29);
			t22= rt-sn1*t22;		t23= it-sn1*t23;
			t28= rt+sn2*t28;		t29= it+sn2*t29;

			aj1p10r =t4;		aj1p10i=t5;
			aj1p4r =t10-t23;		aj1p4i =t11+t22;
			aj1p13r=t16-t29;		aj1p13i=t17+t28;
			aj1p7r =t16+t29;		aj1p7i =t17-t28;
			aj1p1r =t10+t23;		aj1p1i =t11-t22;
		#endif

	/*...Now do the carries. Since the outputs would
		normally be getting dispatched to 15 separate blocks of the A-array, we need 15 separate carries.	*/

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
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
		   cmplx_carry_norm_errcheck0(aj1p0r ,aj1p0i ,cy_r0 ,bjmodn0 ,0 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p1r ,aj1p1i ,cy_r1 ,bjmodn1 ,1 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p2r ,aj1p2i ,cy_r2 ,bjmodn2 ,2 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p3r ,aj1p3i ,cy_r3 ,bjmodn3 ,3 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p4r ,aj1p4i ,cy_r4 ,bjmodn4 ,4 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p5r ,aj1p5i ,cy_r5 ,bjmodn5 ,5 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p6r ,aj1p6i ,cy_r6 ,bjmodn6 ,6 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p7r ,aj1p7i ,cy_r7 ,bjmodn7 ,7 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p8r ,aj1p8i ,cy_r8 ,bjmodn8 ,8 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p9r ,aj1p9i ,cy_r9 ,bjmodn9 ,9 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p10r,aj1p10i,cy_r10,bjmodn10,10,prp_mult);
			cmplx_carry_norm_errcheck(aj1p11r,aj1p11i,cy_r11,bjmodn11,11,prp_mult);
			cmplx_carry_norm_errcheck(aj1p12r,aj1p12i,cy_r12,bjmodn12,12,prp_mult);
			cmplx_carry_norm_errcheck(aj1p13r,aj1p13i,cy_r13,bjmodn13,13,prp_mult);
			cmplx_carry_norm_errcheck(aj1p14r,aj1p14i,cy_r14,bjmodn14,14,prp_mult);

			i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
						and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/
		}
		else
		{
			fermat_carry_norm_errcheck(aj1p0r ,aj1p0i ,cy_r0 ,cy_i0 ,ii0 ,bjmodn0 ,0 *n15,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p1r ,aj1p1i ,cy_r1 ,cy_i1 ,ii1 ,bjmodn1 ,1 *n15,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p2r ,aj1p2i ,cy_r2 ,cy_i2 ,ii2 ,bjmodn2 ,2 *n15,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p3r ,aj1p3i ,cy_r3 ,cy_i3 ,ii3 ,bjmodn3 ,3 *n15,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p4r ,aj1p4i ,cy_r4 ,cy_i4 ,ii4 ,bjmodn4 ,4 *n15,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p5r ,aj1p5i ,cy_r5 ,cy_i5 ,ii5 ,bjmodn5 ,5 *n15,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p6r ,aj1p6i ,cy_r6 ,cy_i6 ,ii6 ,bjmodn6 ,6 *n15,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p7r ,aj1p7i ,cy_r7 ,cy_i7 ,ii7 ,bjmodn7 ,7 *n15,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p8r ,aj1p8i ,cy_r8 ,cy_i8 ,ii8 ,bjmodn8 ,8 *n15,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p9r ,aj1p9i ,cy_r9 ,cy_i9 ,ii9 ,bjmodn9 ,9 *n15,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p10r,aj1p10i,cy_r10,cy_i10,ii10,bjmodn10,10*n15,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p11r,aj1p11i,cy_r11,cy_i11,ii11,bjmodn11,11*n15,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p12r,aj1p12i,cy_r12,cy_i12,ii12,bjmodn12,12*n15,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p13r,aj1p13i,cy_r13,cy_i13,ii13,bjmodn13,13*n15,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p14r,aj1p14i,cy_r14,cy_i14,ii14,bjmodn14,14*n15,NRTM1,NRT_BITS,prp_mult);
		}

	/* The radix-15 DIF pass is here: */

		// Gather the needed data (5 64-bit complex, i.e. 10 64-bit reals) and do the first set of three length-5 transforms...
		#if 1
		RADIX_15_DIF(
			aj1p0r,aj1p0i,aj1p1r,aj1p1i,aj1p2r,aj1p2i,aj1p3r,aj1p3i,aj1p4r,aj1p4i,aj1p5r,aj1p5i,aj1p6r,aj1p6i,aj1p7r,aj1p7i,aj1p8r,aj1p8i,aj1p9r,aj1p9i,aj1p10r,aj1p10i,aj1p11r,aj1p11i,aj1p12r,aj1p12i,aj1p13r,aj1p13i,aj1p14r,aj1p14i,
			a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14]
		)
		#else
		  #if PFETCH
			add0 = &a[j1]; prefetch_p_doubles(add0);
		  #endif
		/* ...Block 1:	*/
			t0 =aj1p0r;				t1 =aj1p0i;
			t2 =aj1p12r;			t3 =aj1p12i;
			rt =aj1p3r; 			it =aj1p3i ;
			t6 =t2 -rt;				t7 =t3 -it;
			t2 =t2 +rt;				t3 =t3 +it;
			t4 =aj1p9r;				t5 =aj1p9i;
			rt =aj1p6r;				it =aj1p6i;
			t8 =t4 -rt;				t9 =t5 -it;
			t4 =t4 +rt;				t5 =t5 +it;
		  #if PFETCH
			addr = add0+p1; prefetch_p_doubles(addr);
		  #endif
			rt = t2+t4;				it = t3+t5;
			t0 = t0+rt;				t1 = t1+it;
			rt = t0+cn1*rt;			it = t1+cn1*it;
			t4 = cn2*(t2-t4);			t5 = cn2*(t3-t5);
			t2 = rt+t4;				t3 = it+t5;
			t4 = rt-t4;				t5 = it-t5;
			rt = ss3*(t6-t8);			it = ss3*(t7-t9);
		  #if PFETCH
			addr = add0+p2; prefetch_p_doubles(addr);
		  #endif
			t8 = rt+sn1*t8;			t9 = it+sn1*t9;
			t6 = rt-sn2*t6;			t7 = it-sn2*t7;
			rt=t8;				it=t9;
			t8 =t2+it;				t9 =t3-rt;	/* <==prefer these to be stored in t8,9	*/
			t2 =t2-it;				t3 =t3+rt;
			rt=t6;				it=t7;
			t6 =t4+it;				t7 =t5-rt;	/* <==prefer these to be stored in t6,7	*/
			t4 =t4-it;				t5 =t5+rt;
		  #if PFETCH
			addr = add0+p3; prefetch_p_doubles(addr);
		  #endif
		/* ...Block 2:	*/
			t10=aj1p10r;			t11=aj1p10i;
			t12=aj1p7r;				t13=aj1p7i;
			rt =aj1p13r;			it =aj1p13i;
			t16=t12 -rt;			t17=t13 -it;
			t12=t12 +rt;			t13=t13 +it;
			t14=aj1p4r;				t15=aj1p4i;
			rt =aj1p1r;				it =aj1p1i;
			t18=t14 -rt;			t19=t15 -it;
			t14=t14 +rt;			t15=t15 +it;
		  #if PFETCH
			addr = add0+p4; prefetch_p_doubles(addr);
		  #endif
			rt = t12+t14;			it = t13+t15;
			t10= t10+rt;			t11= t11+it;
			rt = t10+cn1*rt;			it = t11+cn1*it;
			t14= cn2*(t12-t14);			t15= cn2*(t13-t15);
			t12= rt+t14;			t13= it+t15;
			t14= rt-t14;			t15= it-t15;
			rt = ss3*(t16-t18);			it = ss3*(t17-t19);
		  #if PFETCH
			addr = add0+p5; prefetch_p_doubles(addr);
		  #endif
			t18= rt+sn1*t18;			t19= it+sn1*t19;
			t16= rt-sn2*t16;			t17= it-sn2*t17;
			rt =t18;				it =t19;
			t18=t12+it;				t19=t13-rt;
			t12=t12-it;				t13=t13+rt;
			rt =t16;				it =t17;
			t16=t14+it;				t17=t15-rt;
			t14=t14-it;				t15=t15+rt;
		  #if PFETCH
			addr = add0+p6; prefetch_p_doubles(addr);
		  #endif
		/* ...Block 3:	*/
			t20=aj1p5r;				t21=aj1p5i;
			t22=aj1p2r;				t23=aj1p2i;
			rt =aj1p8r;				it =aj1p8i;
			t26=t22 -rt;			t27=t23 -it;
			t22=t22 +rt;			t23=t23 +it;
			t24=aj1p14r;			t25=aj1p14i;
			rt =aj1p11r;			it =aj1p11i;
			t28=t24 -rt;			t29=t25 -it;
			t24=t24 +rt;			t25=t25 +it;
		  #if PFETCH
			addr = add0+p7; prefetch_p_doubles(addr);
		  #endif
			rt = t22+t24;			it = t23+t25;
			t20= t20+rt;			t21= t21+it;
			rt = t20+cn1*rt;			it = t21+cn1*it;
			t24= cn2*(t22-t24);			t25= cn2*(t23-t25);
			t22= rt+t24;			t23= it+t25;
			t24= rt-t24;			t25= it-t25;
			rt = ss3*(t26-t28);			it = ss3*(t27-t29);
		  #if PFETCH
			addr = add0+p8; prefetch_p_doubles(addr);
		  #endif
			t28= rt+sn1*t28;			t29= it+sn1*t29;
			t26= rt-sn2*t26;			t27= it-sn2*t27;
			rt =t28;				it =t29;
			t28=t22+it;				t29=t23-rt;
			t22=t22-it;				t23=t23+rt;
			rt =t26;				it =t27;
			t26=t24+it;				t27=t25-rt;
			t24=t24-it;				t25=t25+rt;
		  #if PFETCH
			addr = add0+p9; prefetch_p_doubles(addr);
		  #endif
	// ...and now do five radix-3 transforms:
		/* ...Block 1:	*/
			rt =t20;			it =t21;
			t20=t10-rt;			t21=t11-it;
			t10=t10+rt;			t11=t11+it;
			t0 =t0+t10;			t1 =t1+t11;
			a[j1   ]=t0;			a[j2   ]=t1;
			t10=t0+c3m1*t10;		t11=t1+c3m1*t11;
			rt =s*t20;			it =s*t21;
			a[j1+p1]=t10-it;		a[j2+p1]=t11+rt;
			a[j1+p2]=t10+it;		a[j2+p2]=t11-rt;
		  #if PFETCH
			addr = add0+p10; prefetch_p_doubles(addr);
		  #endif
		/* ...Block 2:	*/
			rt =t22;			it =t23;
			t22=t12-rt;			t23=t13-it;
			t12=t12+rt;			t13=t13+it;
			t2 =t2+t12;			t3 =t3+t13;
			a[j1+p13]=t2;			a[j2+p13]=t3;
			t12=t2+c3m1*t12;		t13=t3+c3m1*t13;
			rt =s*t22;			it =s*t23;
			a[j1+p14]=t12-it;		a[j2+p14]=t13+rt;
			a[j1+p12]=t12+it;		a[j2+p12]=t13-rt;
		  #if PFETCH
			addr = add0+p11; prefetch_p_doubles(addr);
		  #endif
		/* ...Block 3:	*/
			rt =t24;			it =t25;
			t24=t14-rt;			t25=t15-it;
			t14=t14+rt;			t15=t15+it;
			t4 =t4+t14;			t5 =t5+t15;
			a[j1+p9 ]=t4;			a[j2+p9 ]=t5;
			t14=t4+c3m1*t14;		t15=t5+c3m1*t15;
			rt =s*t24;			it =s*t25;
			a[j1+p10]=t14-it;		a[j2+p10]=t15+rt;
			a[j1+p11]=t14+it;		a[j2+p11]=t15-rt;
		  #if PFETCH
			addr = add0+p12; prefetch_p_doubles(addr);
		  #endif
		/* ...Block 4:	*/
			rt =t26;			it =t27;
			t26=t16-rt;			t27=t17-it;
			t16=t16+rt;			t17=t17+it;
			t6 =t6+t16;			t7 =t7+t17;
			a[j1+p8 ]=t6;			a[j2+p8 ]=t7;
			t16=t6+c3m1*t16;		t17=t7+c3m1*t17;
			rt =s*t26;			it =s*t27;
			a[j1+p6 ]=t16-it;		a[j2+p6 ]=t17+rt;
			a[j1+p7 ]=t16+it;		a[j2+p7 ]=t17-rt;
		  #if PFETCH
			addr = add0+p13; prefetch_p_doubles(addr);
		  #endif
		/* ...Block 5:	*/
			rt =t28;			it =t29;
			t28=t18-rt;			t29=t19-it;
			t18=t18+rt;			t19=t19+it;
			t8 =t8+t18;			t9 =t9+t19;
			a[j1+p4 ]=t8;			a[j2+p4 ]=t9;
			t18=t8+c3m1*t18;		t19=t9+c3m1*t19;
			rt =s*t28;			it =s*t29;
			a[j1+p5 ]=t18-it;		a[j2+p5 ]=t19+rt;
			a[j1+p3 ]=t18+it;		a[j2+p3 ]=t19-rt;
		  #if PFETCH
			addr = add0+p14; prefetch_p_doubles(addr);
		  #endif
		#endif
		}

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE) {
			jstart += nwt;
			jhi    += nwt;
			col += 15;
			co3 -= 15;
		}
	}

	if(full_pass) {
	//	printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
	} else {
		break;
	}

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-15 forward DIF FFT of the first block of 15 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 15 outputs of (1);
!   (3) Reweight and perform a radix-15 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 15 elements and repeat (1-4).
*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t1    = cy_r14;
		cy_r14= cy_r13;
		cy_r13= cy_r12;
		cy_r12= cy_r11;
		cy_r11= cy_r10;
		cy_r10= cy_r9 ;
		cy_r9 = cy_r8 ;
		cy_r8 = cy_r7 ;
		cy_r7 = cy_r6 ;
		cy_r6 = cy_r5 ;
		cy_r5 = cy_r4 ;
		cy_r4 = cy_r3 ;
		cy_r3 = cy_r2 ;
		cy_r2 = cy_r1 ;
		cy_r1 = cy_r0 ;
		cy_r0 =    t1 ;
	}
	else
	{
		/* ...The 2 Mo"bius carries are here: */
		t1    = cy_r14;	t2    = cy_i14;
		cy_r14= cy_r13;	cy_i14= cy_i13;
		cy_r13= cy_r12;	cy_i13= cy_i12;
		cy_r12= cy_r11;	cy_i12= cy_i11;
		cy_r11= cy_r10;	cy_i11= cy_i10;
		cy_r10= cy_r9 ;	cy_i10= cy_i9 ;
		cy_r9 = cy_r8 ;	cy_i9 = cy_i8 ;
		cy_r8 = cy_r7 ;	cy_i8 = cy_i7 ;
		cy_r7 = cy_r6 ;	cy_i7 = cy_i6 ;
		cy_r6 = cy_r5 ;	cy_i6 = cy_i5 ;
		cy_r5 = cy_r4 ;	cy_i5 = cy_i4 ;
		cy_r4 = cy_r3 ;	cy_i4 = cy_i3 ;
		cy_r3 = cy_r2 ;	cy_i3 = cy_i2 ;
		cy_r2 = cy_r1 ;	cy_i2 = cy_i1 ;
		cy_r1 = cy_r0 ;	cy_i1 = cy_i0 ;
		cy_r0 =   -t2 ;	cy_i0 =   +t1 ;
	}

	full_pass = 0;
	scale = prp_mult = 1;

	/*
	For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
	*/
	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		jhi =15;
	} else {
		jhi = 7;
	}

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
		a[j+p13] *= radix_inv;
		a[j+p14] *= radix_inv;
	}
}

	if(fabs(cy_r0)+fabs(cy_r1)+fabs(cy_r2)+fabs(cy_r3)+fabs(cy_r4)+fabs(cy_r5)+fabs(cy_r6)+fabs(cy_r7)+fabs(cy_r8)+fabs(cy_r9)+fabs(cy_r10)+fabs(cy_r11)+fabs(cy_r12)+fabs(cy_r13)+fabs(cy_r14)
		+fabs(cy_i0)+fabs(cy_i1)+fabs(cy_i2)+fabs(cy_i3)+fabs(cy_i4)+fabs(cy_i5)+fabs(cy_i6)+fabs(cy_i7)+fabs(cy_i8)+fabs(cy_i9)+fabs(cy_i10)+fabs(cy_i11)+fabs(cy_i12)+fabs(cy_i13)+fabs(cy_i14) != 0.0)
	{
		sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in radix15_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
		mlucas_fprint(cbuf,INTERACT);
		err=ERR_CARRY;
		return(err);
	}
	*fracmax = maxerr;
	return(0);
}

/***************/

void radix15_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-15 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing,
!   and radix3,5,15_dif_pass1 for details on the algorithm.
*/
	int j,j1,j2;
	static int n15,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14, first_entry=TRUE;
	static double	c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/

					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */

	if(!first_entry && (n/15) != n15)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n15=n/15;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n15;
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
		p13= p12+p1;
		p14= p13+p1;

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
		p13= p13+ ( (p13>> DAT_BITS) << PAD_BITS );
		p14= p14+ ( (p14>> DAT_BITS) << PAD_BITS );
	}

/*...The radix-15 pass is here.	*/

	for(j=0; j < n15; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*       gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do three length-5 transforms...	*/
	/*
	Twiddleless version requires us to swap inputs as follows:
	indices  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14
		  -> 0, -5,-10, 12,  7,  2,  9,  4, -1,  6,  1, -4,  3, -2, -7
		  == 0, 10,  5, 12,  7,  2,  9,  4, 14,  6,  1, 11,  3, 13,  8 modulo 15.
	I.e. start out with first quintet of indices {0,3,6,9,12}, permute those according to
	{0,3,6,9,12}*14%15 = {0,12,9,6,3}, then each is head of a length-3 list of indices with decrement 5.

	Upshot: swap x3 <-> x12, x6 <-> x9, x1 <-> x10, x4 <-> x7, x2 <-> x5, x8 <-> x14.

	In other words, the input permutation is

		0     0
		1     10
		2     5
		3     12
		4     7
		5     2
		6     9
		7  => 4
		8     14
		9     6
		10    1
		11    11
		12    3
		13    13
		14    8

	The output permutation is

		0     0
		1     1
		2     2
		3     13
		4     14
		5     12
		6     9
		7  => 10
		8     11
		9     8
		10    6
		11    7
		12    4
		13    5
		14    3
	*/

	#if 1
		 /*  inputs, intermediates, outputs */
		RADIX_15_DIF(
			a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],
			a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14]
		);

	#else

	/*...Block 1:	*/
		t0 =a[j1    ];			t1 =a[j2    ];
		t2 =a[j1+p12];			t3 =a[j2+p12];
		rt =a[j1+p3 ];			it =a[j2+p3 ];
		t6 =t2 -rt;				t7 =t3 -it;
		t2 =t2 +rt;				t3 =t3 +it;
		t4 =a[j1+p9 ];			t5 =a[j2+p9 ];
		rt =a[j1+p6 ];			it =a[j2+p6 ];
		t8 =t4 -rt;				t9 =t5 -it;
		t4 =t4 +rt;				t5 =t5 +it;

		rt = t2+t4;				it = t3+t5;
		t0 = t0+rt;				t1 = t1+it;
		rt = t0+cn1*rt;			it = t1+cn1*it;
		t4 = cn2*(t2-t4);		t5 = cn2*(t3-t5);
		t2 = rt+t4;				t3 = it+t5;
		t4 = rt-t4;				t5 = it-t5;
		rt = ss3*(t6-t8);		it = ss3*(t7-t9);
		t8 = rt+sn1*t8;			t9 = it+sn1*t9;
		t6 = rt-sn2*t6;			t7 = it-sn2*t7;
		rt=t8;					it=t9;
		t8 =t2+it;				t9 =t3-rt;	/*<==prefer these to be stored in t8,9	*/
		t2 =t2-it;				t3 =t3+rt;
		rt=t6;					it=t7;
		t6 =t4+it;				t7 =t5-rt;	/*<==prefer these to be stored in t6,7	*/
		t4 =t4-it;				t5 =t5+rt;

	/*...Block 2:	*/
		t10=a[j1+p10];			t11=a[j2+p10];
		t12=a[j1+p7 ];			t13=a[j2+p7 ];
		rt =a[j1+p13];			it =a[j2+p13];
		t16=t12 -rt;			t17=t13 -it;
		t12=t12 +rt;			t13=t13 +it;
		t14=a[j1+p4 ];			t15=a[j2+p4 ];
		rt =a[j1+p1 ];			it =a[j2+p1 ];
		t18=t14 -rt;			t19=t15 -it;
		t14=t14 +rt;			t15=t15 +it;

		rt = t12+t14;			it = t13+t15;
		t10= t10+rt;			t11= t11+it;
		rt = t10+cn1*rt;		it = t11+cn1*it;
		t14= cn2*(t12-t14);		t15= cn2*(t13-t15);
		t12= rt+t14;			t13= it+t15;
		t14= rt-t14;			t15= it-t15;
		rt = ss3*(t16-t18);		it = ss3*(t17-t19);
		t18= rt+sn1*t18;		t19= it+sn1*t19;
		t16= rt-sn2*t16;		t17= it-sn2*t17;
		rt =t18;				it =t19;
		t18=t12+it;				t19=t13-rt;
		t12=t12-it;				t13=t13+rt;
		rt =t16;				it =t17;
		t16=t14+it;				t17=t15-rt;
		t14=t14-it;				t15=t15+rt;

	/*...Block 3:	*/
		t20=a[j1+p5 ];			t21=a[j2+p5 ];
		t22=a[j1+p2 ];			t23=a[j2+p2 ];
		rt =a[j1+p8 ];			it =a[j2+p8 ];
		t26=t22 -rt;			t27=t23 -it;
		t22=t22 +rt;			t23=t23 +it;
		t24=a[j1+p14];			t25=a[j2+p14];
		rt =a[j1+p11];			it =a[j2+p11];
		t28=t24 -rt;			t29=t25 -it;
		t24=t24 +rt;			t25=t25 +it;

		rt = t22+t24;			it = t23+t25;
		t20= t20+rt;			t21= t21+it;
		rt = t20+cn1*rt;		it = t21+cn1*it;
		t24= cn2*(t22-t24);		t25= cn2*(t23-t25);
		t22= rt+t24;			t23= it+t25;
		t24= rt-t24;			t25= it-t25;
		rt = ss3*(t26-t28);		it = ss3*(t27-t29);
		t28= rt+sn1*t28;		t29= it+sn1*t29;
		t26= rt-sn2*t26;		t27= it-sn2*t27;
		rt =t28;				it =t29;
		t28=t22+it;				t29=t23-rt;
		t22=t22-it;				t23=t23+rt;
		rt =t26;				it =t27;
		t26=t24+it;				t27=t25-rt;
		t24=t24-it;				t25=t25+rt;

	/*...and now do five radix-3 transforms:	*/
	/*...Block 1:	*/
		rt =t20;				it =t21;
		t20=t10-rt;				t21=t11-it;
		t10=t10+rt;				t11=t11+it;
		t0 =t0+t10;				t1 =t1+t11;
		a[j1    ]=t0;			a[j2    ]=t1;
		t10=t0+c3m1*t10;		t11=t1+c3m1*t11;
		rt =s*t20;				it =s*t21;
		a[j1+p1 ]=t10-it;		a[j2+p1 ]=t11+rt;
		a[j1+p2 ]=t10+it;		a[j2+p2 ]=t11-rt;

	/*...Block 2:	*/
		rt =t22;				it =t23;
		t22=t12-rt;				t23=t13-it;
		t12=t12+rt;				t13=t13+it;
		t2 =t2+t12;				t3 =t3+t13;
		a[j1+p13]=t2;			a[j2+p13]=t3;
		t12=t2+c3m1*t12;		t13=t3+c3m1*t13;
		rt =s*t22;				it =s*t23;
		a[j1+p14]=t12-it;		a[j2+p14]=t13+rt;
		a[j1+p12]=t12+it;		a[j2+p12]=t13-rt;

	/*...Block 3:	*/
		rt =t24;				it =t25;
		t24=t14-rt;				t25=t15-it;
		t14=t14+rt;				t15=t15+it;
		t4 =t4+t14;				t5 =t5+t15;
		a[j1+p9 ]=t4;			a[j2+p9 ]=t5;
		t14=t4+c3m1*t14;		t15=t5+c3m1*t15;
		rt =s*t24;				it =s*t25;
		a[j1+p10]=t14-it;		a[j2+p10]=t15+rt;
		a[j1+p11]=t14+it;		a[j2+p11]=t15-rt;

	/*...Block 4:	*/
		rt =t26;				it =t27;
		t26=t16-rt;				t27=t17-it;
		t16=t16+rt;				t17=t17+it;
		t6 =t6+t16;				t7 =t7+t17;
		a[j1+p8 ]=t6;			a[j2+p8 ]=t7;
		t16=t6+c3m1*t16;		t17=t7+c3m1*t17;
		rt =s*t26;				it =s*t27;
		a[j1+p6 ]=t16-it;		a[j2+p6 ]=t17+rt;
		a[j1+p7 ]=t16+it;		a[j2+p7 ]=t17-rt;

	/*...Block 5:	*/
		rt =t28;				it =t29;
		t28=t18-rt;				t29=t19-it;
		t18=t18+rt;				t19=t19+it;
		t8 =t8+t18;				t9 =t9+t19;
		a[j1+p4 ]=t8;			a[j2+p4 ]=t9;
		t18=t8+c3m1*t18;		t19=t9+c3m1*t19;
		rt =s*t28;				it =s*t29;
		a[j1+p5 ]=t18-it;		a[j2+p5 ]=t19+rt;
		a[j1+p3 ]=t18+it;		a[j2+p3 ]=t19-rt;
								/* Totals: 17*6+6*10 = 162 FADD, 5*6+2*10 = 50 FMUL.	*/
#endif
	}
}

/***************/

void radix15_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-15 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix15_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
!
!   See the documentation in radix16_dif_pass for further details on the array indexing,
!   and radix3,5,15_dif_pass1 for details on the algorithm.
*/
	int j,j1,j2;
	static int n15,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14, first_entry=TRUE;
	static double c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/

					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */

	if(!first_entry && (n/15) != n15)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n15=n/15;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n15;
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
		p13= p12+p1;
		p14= p13+p1;

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
		p13= p13+ ( (p13>> DAT_BITS) << PAD_BITS );
		p14= p14+ ( (p14>> DAT_BITS) << PAD_BITS );
	}

/*...The radix-15 pass is here.	*/

	for(j=0; j < n15; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*       gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do five radix-3 transforms...	*/
	/*
	Twiddleless version requires us to swap inputs as follows:
	indices  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14
		  -> 0, -3, -6, -9,-12, 10,  7,  4,  1, -2,  5,  2, -1, -4, -7
		  == 0, 12,  9,  6,  3, 10,  7,  4,  1, 13,  5,  2, 14, 11,  8 modulo 15.
	I.e. start out with first triplet of indices {0,5,10}, permute those according to
	{0,5,10}*14%15 = {0,10,5}, then each is head of a length-5 list of indices with decrement 3.

	Remember, inputs to DIT are bit-reversed, so a[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14]
										 contain x[0, 5,10, 1, 6,11, 2, 7,12, 3, 8,13, 4, 9,14], so
	x1 ->x12	means a3 -> a8
	x2 ->x9 	means a6 -> a13
	x3 ->x6 	means a9 -> a4
	x4 ->x3 	means a12-> a9
	x5 ->x10	means a1 -> a2
	x6 ->x7 	means a4 -> a7
	x7 ->x4 	means a7 -> a12
	x8 ->x1 	means a10-> a3
	x9 ->x13	means a13-> a11
	x10->x5 	means a2 -> a1
	x11->x2 	means a5 -> a6
	x12->x14	means a8 -> a14
	x13->x11	means a11-> a5
	x14->x8 	means a14-> a10

	In other words, the input permutation is

		0     0
		1     2
		2     1
		3     8
		4     7
		5     6
		6     13
		7  => 12
		8     14
		9     4
		10    3
		11    5
		12    9
		13    11
		14    10

	The output permutation is

		0     0
		1     5
		2     10
		3     9
		4     14
		5     4
		6     3
		7  => 8
		8     13
		9     12
		10    2
		11    7
		12    6
		13    11
		14    1
	*/

	#if 1

		RADIX_15_DIT(
			a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14],
			a[j1],a[j2],a[j1+p1],a[j2+p1],a[j1+p2],a[j2+p2],a[j1+p3],a[j2+p3],a[j1+p4],a[j2+p4],a[j1+p5],a[j2+p5],a[j1+p6],a[j2+p6],a[j1+p7],a[j2+p7],a[j1+p8],a[j2+p8],a[j1+p9],a[j2+p9],a[j1+p10],a[j2+p10],a[j1+p11],a[j2+p11],a[j1+p12],a[j2+p12],a[j1+p13],a[j2+p13],a[j1+p14],a[j2+p14]
		);

	#else

	/*...Block 1:	*/
		t0 =a[j1    ];			t1 =a[j2    ];
		t2 =a[j1+p2 ];			t3 =a[j2+p2 ];
		rt =a[j1+p1 ];			it =a[j2+p1 ];
		t4 =t2 -rt;				t5 =t3 -it;
		t2 =t2 +rt;				t3 =t3 +it;
		t0 =t0+t2;				t1 =t1+t3;
		t2 =t0+c3m1*t2;			t3 =t1+c3m1*t3;
		rt =s*t4;				it =s*t5;
		t4 =t2-it;				t5 =t3+rt;
		t2 =t2+it;				t3 =t3-rt;

	/*...Block 2:	*/
		t6 =a[j1+p8 ];			t7 =a[j2+p8 ];
		t8 =a[j1+p7 ];			t9 =a[j2+p7 ];
		rt =a[j1+p6 ];			it =a[j2+p6 ];
		t10=t8 -rt;				t11=t9 -it;
		t8 =t8 +rt;				t9 =t9 +it;
		t6 =t6+t8;				t7 =t7+t9;
		t8 =t6+c3m1*t8;			t9 =t7+c3m1*t9;
		rt =s*t10;				it =s*t11;
		t10=t8-it;				t11=t9+rt;
		t8 =t8+it;				t9 =t9-rt;

	/*...Block 3:	*/
		t12=a[j1+p13];			t13=a[j2+p13];
		t14=a[j1+p12];			t15=a[j2+p12];
		rt =a[j1+p14];			it =a[j2+p14];
		t16=t14 -rt;			t17=t15 -it;
		t14=t14 +rt;			t15=t15 +it;
		t12=t12+t14;			t13=t13+t15;
		t14=t12+c3m1*t14;			t15=t13+c3m1*t15;
		rt =s*t16;				it =s*t17;
		t16=t14-it;				t17=t15+rt;
		t14=t14+it;				t15=t15-rt;

	/*...Block 4:	*/
		t18=a[j1+p4 ];			t19=a[j2+p4 ];
		t20=a[j1+p3 ];			t21=a[j2+p3 ];
		rt =a[j1+p5 ];			it =a[j2+p5 ];
		t22=t20 -rt;			t23=t21 -it;
		t20=t20 +rt;			t21=t21 +it;
		t18=t18+t20;			t19=t19+t21;
		t20=t18+c3m1*t20;			t21=t19+c3m1*t21;
		rt =s*t22;				it =s*t23;
		t22=t20-it;				t23=t21+rt;
		t20=t20+it;				t21=t21-rt;

	/*...Block 5:	*/
		t24=a[j1+p9 ];			t25=a[j2+p9 ];
		t26=a[j1+p11];			t27=a[j2+p11];
		rt =a[j1+p10];			it =a[j2+p10];
		t28=t26 -rt;			t29=t27 -it;
		t26=t26 +rt;			t27=t27 +it;
		t24=t24+t26;			t25=t25+t27;
		t26=t24+c3m1*t26;			t27=t25+c3m1*t27;
		rt =s*t28;				it =s*t29;
		t28=t26-it;				t29=t27+rt;
		t26=t26+it;				t27=t27-rt;

	/*...and now do three radix-5 transforms:	*/
	/*...Block 1:	*/
		rt =t24;			it =t25;
		t24=t6-rt;			t25=t7-it;
		t6 =t6+rt;			t7 =t7+it;
		rt =t18;			it =t19;
		t18=t12-rt;			t19=t13-it;
		t12=t12+rt;			t13=t13+it;

		rt = t6+t12;		it = t7+t13;
		t0 = t0+rt;			t1 = t1+it;
		rt = t0+cn1*rt;		it = t1+cn1*it;
		t12= cn2*(t6 -t12);		t13= cn2*(t7 -t13);
		t6 = rt+t12;		t7 = it+t13;
		t12= rt-t12;		t13= it-t13;
		rt = ss3*(t18-t24);		it = ss3*(t19-t25);
		t18= rt-sn1*t18;		t19= it-sn1*t19;
		t24= rt+sn2*t24;		t25= it+sn2*t25;

		a[j1    ]=t0;		a[j2    ]=t1;
		a[j1+p9 ]=t6-t19;		a[j2+p9 ]=t7+t18;
		a[j1+p3 ]=t12-t25;		a[j2+p3 ]=t13+t24;
		a[j1+p12]=t12+t25;		a[j2+p12]=t13-t24;
		a[j1+p6 ]=t6+t19;		a[j2+p6 ]=t7-t18;

	/*...Block 2:	*/
		rt =t26;			it =t27;
		t26=t8-rt;			t27=t9-it;
		t8 =t8+rt;			t9 =t9+it;
		rt =t20;			it =t21;
		t20=t14-rt;			t21=t15-it;
		t14=t14+rt;			t15=t15+it;

		rt = t8+t14;		it = t9+t15;
		t2 = t2+rt;			t3 = t3+it;
		rt = t2+cn1*rt;		it = t3+cn1*it;
		t14= cn2*(t8-t14);		t15= cn2*(t9-t15);
		t8 = rt+t14;		t9 = it+t15;
		t14= rt-t14;		t15= it-t15;
		rt = ss3*(t20-t26);		it = ss3*(t21-t27);
		t20= rt-sn1*t20;		t21= it-sn1*t21;
		t26= rt+sn2*t26;		t27= it+sn2*t27;

		a[j1+p5 ]=t2;		a[j2+p5 ]=t3;
		a[j1+p14]=t8-t21;		a[j2+p14]=t9+t20;
		a[j1+p8 ]=t14-t27;		a[j2+p8 ]=t15+t26;
		a[j1+p2 ]=t14+t27;		a[j2+p2 ]=t15-t26;
		a[j1+p11]=t8+t21;		a[j2+p11]=t9-t20;

	/*...Block 3:	*/
		rt =t28;			it =t29;
		t28=t10-rt;			t29=t11-it;
		t10=t10+rt;			t11=t11+it;
		rt =t22;			it =t23;
		t22=t16-rt;			t23=t17-it;
		t16=t16+rt;			t17=t17+it;

		rt = t10+t16;		it = t11+t17;
		t4 = t4+rt;			t5 = t5+it;
		rt = t4+cn1*rt;		it = t5+cn1*it;
		t16= cn2*(t10-t16);		t17= cn2*(t11-t17);
		t10= rt+t16;		t11= it+t17;
		t16= rt-t16;		t17= it-t17;
		rt = ss3*(t22-t28);		it = ss3*(t23-t29);
		t22= rt-sn1*t22;		t23= it-sn1*t23;
		t28= rt+sn2*t28;		t29= it+sn2*t29;

		a[j1+p10]=t4;		a[j2+p10]=t5;
		a[j1+p4 ]=t10-t23;		a[j2+p4 ]=t11+t22;
		a[j1+p13]=t16-t29;		a[j2+p13]=t17+t28;
		a[j1+p7 ]=t16+t29;		a[j2+p7 ]=t17-t28;
		a[j1+p1 ]=t10+t23;		a[j2+p1 ]=t11-t22;
	#endif
	}
}

