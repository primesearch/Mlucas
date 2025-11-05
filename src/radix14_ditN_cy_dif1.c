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

int radix14_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-14 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-14 complex DIF pass on the data in
!   the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix7/8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n14,bjmodn0,bjmodn1,bjmodn2,bjmodn3,bjmodn4,bjmodn5,bjmodn6,bjmodn7,bjmodn8,bjmodn9,bjmodn10
		,bjmodn11,bjmodn12,bjmodn13,i,j,j1,j2,jstart,jhi,root_incr,k1,k2,k,khi,l,outer;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13, nsave = 0;
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
#if LO_ADD
	static double cc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					ss1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					cc2 =-.22252093395631440426,	 /* cos(2u)	*/
					ss2 = .97492791218182360702,	 /* sin(2u)	*/
					cc3 =-.90096886790241912622,	 /* cos(3u)	*/
					ss3 = .43388373911755812050;	 /* sin(3u)	*/
					/*
					Radix-14-specific roots are all symmetric with Radix-7-specific ones.
					cs  = .90096886790241912623 = -cc3	Real part of exp(i*2*pi/14), the radix-14 fundamental sincos datum.
					sn  = .43388373911755812047 = +ss3	Imag part of exp(i*2*pi/14).
					cs3 = .22252093395631440430 = -cc2
					sn3 = .97492791218182360701 = +ss2
					*/
	double re,im;
#else
	static double cc0  =-0.16666666666666666667,	/* (cc1+cc2+cc3)/3	*/
					cc1  = 1.52445866976115265675, 	/*  cc1-cc3		*/
					cc2  = 0.67844793394610472196, 	/*  cc2-cc3		*/
					cc3  = 0.73430220123575245957,	/* (cc1+cc2-2*cc3)/3	*/
		/* Switch the sign of ss3 in these: */
					ss0  = 0.44095855184409843174,	/* (ss1+ss2-ss3)/3	*/
					ss1  = 1.21571522158558792920, 	/*  ss1+ss3		*/
					ss2  = 1.40881165129938172752, 	/*  ss2+ss3		*/
					ss3  = 0.87484229096165655224;	/* (ss1+ss2+2*ss3)/3	*/
#endif
	static double radix_inv, n2inv;
	double rt,it
		,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28
		,aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r,aj1p6r,aj1p7r,aj1p8r,aj1p9r,aj1p10r,aj1p11r,aj1p12r,aj1p13r
		,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i,aj1p6i,aj1p7i,aj1p8i,aj1p9i,aj1p10i,aj1p11i,aj1p12i,aj1p13i
		,cy_r0,cy_r1,cy_r2,cy_r3,cy_r4,cy_r5,cy_r6,cy_r7,cy_r8,cy_r9,cy_r10,cy_r11,cy_r12,cy_r13
		,cy_i0,cy_i1,cy_i2,cy_i3,cy_i4,cy_i5,cy_i6,cy_i7,cy_i8,cy_i9,cy_i10,cy_i11,cy_i12,cy_i13
		,temp,frac,scale;
	double maxerr = 0.0;
#if PFETCH
	double *add0, *addr;
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	int ii0,ii1,ii2,ii3,ii4,ii5,ii6,ii7,ii8,ii9,ii10,ii11,ii12,ii13;	/* indices into weights arrays (mod NWT) */
	double wt_re,wt_im, wi_re,wi_im;					/* Fermat/LOACC weights stuff */

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=ii0=ii1=ii2=ii3=ii4=ii5=ii6=ii7=ii8=ii9=ii10=ii11=ii12=ii13=-1;

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

/*...change n14 and n_div_wt to non-static to work around a gcc compiler bug. */
	n14   = n/14;
	n_div_nwt = n14 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n14)
	{
		sprintf(cbuf,"ERROR: iter = %10d; NWT_BITS does not divide N/14 in radix14_ditN_cy_dif1.\n",iter);
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
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)14));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

		/*   constant index offsets for load/stores are here.	*/

		p1 = n14;
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

		/* For Fermat-mod, since 'adjacent' words are actually stride-2 separated
		in terms of the floating residue array, block boundaries have half the i-index
		(e.g. as in sw*i and bw*i) value they do in the Mersenne-mod case:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			bjmodnini=0;
			for(j=0; j < n14; j++)
			{
				bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
			}
		}
		else
		{
			bjmodnini=0;
			for(j=0; j < n14/2; j++)
			{
				bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
			}
		}
	}

/*...The radix-14 final DIT pass is here.	*/

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

	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		cy_r0 = -2;
	}

	*fracmax=0;	/* init max. fractional error	*/

	root_incr = 1;	/* init sincos array index increment (set = 1 for normal carry pass, = 0 for wrapper pass)	*/

	scale = n2inv;	// init inverse-weight scale factor = 2/n for normal carry pass, 1 for wrapper pass

	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		jstart = 0;
		jhi = jstart+nwt-1;
		khi = n_div_nwt;
	}
	else
	{
		jstart = 0;
		jhi = n_div_nwt;
		khi = 1;
	}

for(outer=0; outer <= 1; outer++)
{
	i = 0;		/* Index into the BASE and BASEINV arrays. */
	/* If bw > 0 (i.e. n does not divide p), lowest-order digit always a bigword:	*/
	if(bw > 0)
		i = 1;

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
	bjmodn13= bjmodn12+bjmodnini-n; bjmodn13= bjmodn13+ ( (-(int)((uint32)bjmodn13>> 31)) & n);

	/* For Fermat-mod, IBDWT access patterns repeat with period NWT = {odd part of radix0},
	so for even radix0's only really need that many bjmodn and ii's, but that would require
	specialized carry macros tha`t don't update ii and bjmodn - not worth the trouble.
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		col=0;
		co2=(n >> nwt_bits)-1+14;
		co3=co2-14;		/* At the start of each new j-loop, co3=co2-radix(1)	*/
	}
	else
	{
		/* indices into IBDWT weights arrays (mod NWT) is here: */
		ii0 = 0;
		ii1 = (SW_DIV_N*n14/2) % nwt;
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

		/* Start this value off at N in Fermat-mod case, so (bjmodn >= sw) check in
		fermat_carry_norm_errcheck (cf. carry.h) yields a bigword (i == 1) for j= 0:
		*/
		bjmodn0 = n;
		bjmodn7 = n;
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

			t1 =a[j1    ];	t2 =a[j2    ];
			rt =a[j1+p1 ];	it =a[j2+p1 ];
			t3 =t1 -rt;		t4 =t2 -it;
			t1 =t1 +rt;		t2 =t2 +it;

			t5 =a[j1+p11];	t6 =a[j2+p11];
			rt =a[j1+p10];	it =a[j2+p10];
			t7 =t5 -rt;		t8 =t6 -it;
			t5 =t5 +rt;		t6 =t6 +it;

			t9 =a[j1+p7 ];	t10=a[j2+p7 ];
			rt =a[j1+p6 ];	it =a[j2+p6 ];
			t11=t9 -rt;		t12=t10-it;
			t9 =t9 +rt;		t10=t10+it;

			t13=a[j1+p3 ];	t14=a[j2+p3 ];
			rt =a[j1+p2 ];	it =a[j2+p2 ];
			t15=t13-rt;		t16=t14-it;
			t13=t13+rt;		t14=t14+it;

			t17=a[j1+p12];	t18=a[j2+p12];
			rt =a[j1+p13];	it =a[j2+p13];
			t19=t17-rt;		t20=t18-it;
			t17=t17+rt;		t18=t18+it;

			t21=a[j1+p8 ];	t22=a[j2+p8 ];
			rt =a[j1+p9 ];	it =a[j2+p9 ];
			t23=t21-rt;		t24=t22-it;
			t21=t21+rt;		t22=t22+it;

			t25=a[j1+p4 ];	t26=a[j2+p4 ];
			rt =a[j1+p5 ];	it =a[j2+p5 ];
			t27=t25-rt;		t28=t26-it;
			t25=t25+rt;		t26=t26+it;

/*
!...Now do the two radix-7 transforms.
!   normally be getting dispatched to 14 separate blocks of the A-array, we need 14 separate carries.
*/
/*...aj1p[0:12:2]r use t[1:25:4]; aj1p[0:12:2]i use t[2:26:4]	*/

			rt =t25;						it =t26;
			t25=t5 -rt;						t26=t6 -it;
			t5 =t5 +rt;						t6 =t6 +it;
			rt =t21;						it =t22;
			t21=t9 -rt;						t22=t10 -it;
			t9 =t9 +rt;						t10=t10 +it;
			rt =t17;						it =t18;
			t17=t13 -rt;					t18=t14 -it;
			t13=t13 +rt;					t14=t14 +it;

	/* Version 1: direct computation of 3x3 matrix-vector product. */
	#if LO_ADD
			aj1p0r = t1+t5+t9+t13;			aj1p0i = t2+t6+t10+t14;
			rt = t1+cc1*t5+cc2*t9+cc3*t13;	it = t2+cc1*t6+cc2*t10+cc3*t14;
			re = t1+cc2*t5+cc3*t9+cc1*t13;	im = t2+cc2*t6+cc3*t10+cc1*t14;
			t1 = t1+cc3*t5+cc1*t9+cc2*t13;	t2 = t2+cc3*t6+cc1*t10+cc2*t14;

			t5 = ss1*t25+ss2*t21+ss3*t17;	t6 = ss1*t26+ss2*t22+ss3*t18;
			t9 = ss2*t25-ss3*t21-ss1*t17;	t10= ss2*t26-ss3*t22-ss1*t18;
			t13= ss3*t25-ss1*t21+ss2*t17;	t14= ss3*t26-ss1*t22+ss2*t18;
		/* Output permutation causes signs to get flipped here: */
			aj1p6r =rt+t6;					aj1p6i =it-t5;
			aj1p12r=re+t10;					aj1p12i=im-t9;
			aj1p4r =t1+t14;					aj1p4i =t2-t13;
			aj1p10r=t1-t14;					aj1p10i=t2+t13;
			aj1p2r =re-t10;					aj1p2i =im+t9;
			aj1p8r =rt-t6;					aj1p8i =it+t5;

		/* Version 2: */
	#else
		/* Get cyclic convolution of (cc1,cc2,cc3) with (t5,t13,t9) and with (t6,t14,t10),
			 (note switch of last 2 t-elements in each case), add t1 and t2 to respective outputs: */

			rt = t5+t13+t9; aj1p0r=rt+t1;	it = t6+t14+t10; aj1p0i=it+t2;
			t1 = rt*cc0+t1;					t2 = it*cc0+t2;
			t5 = t5-t9;						t6 = t6 -t10;
			t9 = t13-t9;					t10= t14-t10;
			t13=(t5+t9)*cc3;				t14=(t6+t10)*cc3;
			t5 = t5*cc1;					t6 = t6 *cc1;
			t9 = t9*cc2;					t10= t10*cc2;
			rt = t5-t13;					it = t6 -t14;
			t9 = t9-t13;					t10= t10-t14;

			t5 = t1-rt-t9;					t6 = t2-it-t10;	/* C2	*/
			t9 = t1+t9;						t10= t2+t10;	/* C3	*/
			t1 = t1+rt;						t2 = t2+it;	/* C1	*/

		/* Get cyclic convolution of (ss1,ss2,-ss3) with (t25,-t17,t21) and with (t26,-t18,t22),
			 (note permutation of t-elements and sign flips on ss3, t17 and t18): */

			t13=(t25-t17 +t21)*ss0;			t14=(t26-t18+t22)*ss0;
			t25= t25-t21;					t26= t26-t22;
			t21= t17+t21;					t22= t18+t22;		/* Flip sign here... */
			t17=(t21-t25)*ss3;				t18=(t22-t26)*ss3;	/* ...and here... */
			t25= t25*ss1;					t26= t26*ss1;
			t21= t21*ss2;					t22= t22*ss2;
			t25= t17+t25;					t26= t18+t26;		/* Sign flips get undone in these next 2. */
			t21= t17-t21;					t22= t18-t22;

			t17= t13-t25-t21;				t18= t14-t26-t22;	/* S2 */
			t21= t13+t21;					t22= t14+t22;		/*-S3 */
			t13= t13+t25;					t14= t14+t26;		/* S1 */

			aj1p6r =t1+t14;					aj1p6i =t2-t13;
			aj1p12r=t5+t18;					aj1p12i=t6-t17;
			aj1p4r =t9-t22;					aj1p4i =t10+t21;
			aj1p10r=t9+t22;					aj1p10i=t10-t21;
			aj1p2r =t5-t18;					aj1p2i =t6+t17;
			aj1p8r =t1-t14;					aj1p8i =t2+t13;
	#endif
	/*...aj1p[1:13:2]r use t[3:27:4]; aj1p[1:13:2]i use t[4:28:4]	*/

			rt =t27;						it =t28;
			t27=t7 -rt;						t28=t8 -it;
			t7 =t7 +rt;						t8 =t8 +it;
			rt =t23;						it =t24;
			t23=t11 -rt;					t24=t12 -it;
			t11=t11 +rt;					t12=t12 +it;
			rt =t19;						it =t20;
			t19=t15 -rt;					t20=t16 -it;
			t15=t15 +rt;					t16=t16 +it;

		/* Version 1: direct computation of 3x3 matrix-vector product. */
	#if LO_ADD
			aj1p7r = t3+t7+t11+t15;			aj1p7i = t4+t8+t12+t16;
			rt = t3+cc1*t7+cc2*t11+cc3*t15;	it = t4+cc1*t8+cc2*t12+cc3*t16;
			re = t3+cc2*t7+cc3*t11+cc1*t15;	im = t4+cc2*t8+cc3*t12+cc1*t16;
			t3 = t3+cc3*t7+cc1*t11+cc2*t15;	t4 = t4+cc3*t8+cc1*t12+cc2*t16;

			t7 = ss1*t27+ss2*t23+ss3*t19;	t8 = ss1*t28+ss2*t24+ss3*t20;
			t11= ss2*t27-ss3*t23-ss1*t19;	t12= ss2*t28-ss3*t24-ss1*t20;
			t15= ss3*t27-ss1*t23+ss2*t19;	t16= ss3*t28-ss1*t24+ss2*t20;
		/* Output permutation causes signs to get flipped here: */
			aj1p13r=rt+t8 ;					aj1p13i=it-t7 ;
			aj1p5r =re+t12;					aj1p5i =im-t11;
			aj1p11r=t3+t16;					aj1p11i=t4-t15;
			aj1p3r =t3-t16;					aj1p3i =t4+t15;
			aj1p9r =re-t12;					aj1p9i =im+t11;
			aj1p1r =rt-t8 ;					aj1p1i =it+t7 ;

		/* Version 2: */
	#else
		/* Get cyclic convolution of (cc1,cc2,cc3) with (t7,t15,t11) and with (t8,t16,t12),
			 (note switch of last 2 t-elements in each case), add t1 and t2 to respective outputs: */

			rt = t7+t15+t11;aj1p7r=rt+t3;	it = t8+t16+t12; aj1p7i=it+t4;
			t3 = rt*cc0+t3;					t4 = it*cc0+t4;
			t7 = t7-t11;					t8 = t8 -t12;
			t11= t15-t11;					t12= t16-t12;
			t15=(t7+t11)*cc3;				t16=(t8+t12)*cc3;
			t7 = t7*cc1;					t8 = t8 *cc1;
			t11= t11*cc2;					t12= t12*cc2;
			rt = t7 -t15;					it = t8 -t16;
			t11= t11-t15;					t12= t12-t16;

			t7 = t3-rt-t11;					t8 = t4-it-t12;	/* C2	*/
			t11= t3+t11;					t12= t4+t12;	/* C3	*/
			t3 = t3+rt;						t4 = t4+it;	/* C1	*/

		/* Get cyclic convolution of (ss1,ss2,-ss3) with (t27,-t19,t23) and with (t28,-t20,t24),
			 (note permutation of t-elements and sign flips on ss3, t19 and t20): */

			t15=(t27-t19 +t23)*ss0;			t16=(t28-t20+t24)*ss0;
			t27= t27-t23;					t28= t28-t24;
			t23= t19+t23;					t24= t20+t24;		/* Flip sign here... */
			t19=(t23-t27)*ss3;				t20=(t24-t28)*ss3;	/* ...and here... */
			t27= t27*ss1;					t28= t28*ss1;
			t23= t23*ss2;					t24= t24*ss2;
			t27= t19+t27;					t28= t20+t28;		/* Sign flips get undone in these next 2. */
			t23= t19-t23;					t24= t20-t24;

			t19= t15-t27-t23;				t20= t16-t28-t24;	/* S2 */
			t23= t15+t23;					t24= t16+t24;		/*-S3 */
			t15= t15+t27;					t16= t16+t28;		/* S1 */

			aj1p13r=t3+t16;					aj1p13i=t4-t15;
			aj1p5r =t7+t20;					aj1p5i =t8-t19;
			aj1p11r=t11-t24;				aj1p11i=t12+t23;
			aj1p3r =t11+t24;				aj1p3i =t12-t23;
			aj1p9r =t7-t20;					aj1p9i =t8+t19;
			aj1p1r =t3-t16;					aj1p1i =t4+t15;
	#endif
/*...Now do the carries. Since the outputs would
    normally be getting dispatched to 14 separate blocks of the A-array, we need 14 separate carries.	*/

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

			i =((uint32)(sw - bjmodn0) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
						and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/
		}
		else
		{
			fermat_carry_norm_errcheck(aj1p0r ,aj1p0i ,cy_r0 ,cy_i0 ,ii0 ,bjmodn0 ,0 *n14,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p1r ,aj1p1i ,cy_r1 ,cy_i1 ,ii1 ,bjmodn1 ,1 *n14,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p2r ,aj1p2i ,cy_r2 ,cy_i2 ,ii2 ,bjmodn2 ,2 *n14,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p3r ,aj1p3i ,cy_r3 ,cy_i3 ,ii3 ,bjmodn3 ,3 *n14,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p4r ,aj1p4i ,cy_r4 ,cy_i4 ,ii4 ,bjmodn4 ,4 *n14,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p5r ,aj1p5i ,cy_r5 ,cy_i5 ,ii5 ,bjmodn5 ,5 *n14,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p6r ,aj1p6i ,cy_r6 ,cy_i6 ,ii6 ,bjmodn6 ,6 *n14,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p7r ,aj1p7i ,cy_r7 ,cy_i7 ,ii7 ,bjmodn7 ,7 *n14,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p8r ,aj1p8i ,cy_r8 ,cy_i8 ,ii8 ,bjmodn8 ,8 *n14,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p9r ,aj1p9i ,cy_r9 ,cy_i9 ,ii9 ,bjmodn9 ,9 *n14,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p10r,aj1p10i,cy_r10,cy_i10,ii10,bjmodn10,10*n14,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p11r,aj1p11i,cy_r11,cy_i11,ii11,bjmodn11,11*n14,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p12r,aj1p12i,cy_r12,cy_i12,ii12,bjmodn12,12*n14,NRTM1,NRT_BITS,prp_mult);
			fermat_carry_norm_errcheck(aj1p13r,aj1p13i,cy_r13,cy_i13,ii13,bjmodn13,13*n14,NRTM1,NRT_BITS,prp_mult);
		}

	/*
	!...The radix-14 DIF pass is here:
	!
	!   Twiddleless version requires us to swap inputs as follows:
	!   indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13
	!         -> 0,-7,12, 5,10, 3, 8, 1, 6,-1, 4,-3, 2,-5
	!         == 0, 7,12, 5,10, 3, 8, 1, 6,13, 4,11, 2, 9 modulo 14.
	!   I.e. start out with first septet of indices {0,2,4,6,8,10,12}, permute those according to
	!   {0,2,4,6,8,10,12}*13%14 = {0,12,10,8,6,4,2,}, then each is head of a length-2 list of indices with decrement 7.
	*/
	/*...First radix-7 block uses aj1p[0:12:2] as inputs:	*/
	#if PFETCH
	add0 = &a[j1];
	prefetch_p_doubles(add0);
	#endif
			t1 =aj1p0r;						t2 =aj1p0i ;
			t3 =aj1p12r+aj1p2r ;			t4 =aj1p12i+aj1p2i ;	/* x1+x6 */
			t5 =aj1p10r+aj1p4r ;			t6 =aj1p10i+aj1p4i ;	/* x2+x5 */
			t7 =aj1p8r +aj1p6r ;			t8 =aj1p8i +aj1p6i ;	/* x3+x4 */
	#if PFETCH
	addr = add0+p1;
	prefetch_p_doubles(addr);
	#endif
			t9 =aj1p8r -aj1p6r ;			t10=aj1p8i -aj1p6i ;	/* x3-x4 */
			t11=aj1p10r-aj1p4r ;			t12=aj1p10i-aj1p4i ;	/* x2-x5 */
			t13=aj1p12r-aj1p2r ;			t14=aj1p12i-aj1p2i ;	/* x1-x6 */
	#if PFETCH
	addr = add0+p2;
	prefetch_p_doubles(addr);
	#endif
		/* Version 1: direct computation of 3x3 matrix-vector product. */
	#if LO_ADD
			aj1p0r = t1+t3+t5+t7;			aj1p0i = t2+t4+t6+t8	;
			rt = t1+cc1*t3+cc2*t5+cc3*t7;	it = t2+cc1*t4+cc2*t6+cc3*t8;
			re = t1+cc2*t3+cc3*t5+cc1*t7;	im = t2+cc2*t4+cc3*t6+cc1*t8;
			t1 = t1+cc3*t3+cc1*t5+cc2*t7;	t2 = t2+cc3*t4+cc1*t6+cc2*t8;
	#if PFETCH
	addr = add0+p3;
	prefetch_p_doubles(addr);
	#endif
			t3 = ss1*t13+ss2*t11+ss3*t9;	t4 = ss1*t14+ss2*t12+ss3*t10;
			t5 = ss2*t13-ss3*t11-ss1*t9;	t6 = ss2*t14-ss3*t12-ss1*t10;
			t7 = ss3*t13-ss1*t11+ss2*t9;	t8 = ss3*t14-ss1*t12+ss2*t10;
	#if PFETCH
	addr = add0+p4;
	prefetch_p_doubles(addr);
	#endif
		/* Output permutation causes signs to get flipped here: */
			aj1p2r =rt-t4;					aj1p2i =it+t3;
			aj1p4r =re-t6;					aj1p4i =im+t5;
			aj1p6r =t1-t8;					aj1p6i =t2+t7;
	#if PFETCH
	addr = add0+p5;
	prefetch_p_doubles(addr);
	#endif
			aj1p8r =t1+t8;					aj1p8i =t2-t7;
			aj1p10r=re+t6;					aj1p10i=im-t5;
			aj1p12r=rt+t4;					aj1p12i=it-t3;

		/* Version 2: get cyclic convolutionof (cc1,cc2,cc3) with (t3,t7,t5) and with (t4,t8,t6), (note switch of
		last 2 t-elements in each case), addt1 and t2 to respective outputs. */
	#else
			rt = t3+t7+t5; aj1p0r=rt+t1;	it = t4+t8+t6; aj1p0i=it+t2;	/* X0 = w+x+y+z	*/
			t1 = rt*cc0+t1;					t2 = it*cc0+t2;	/* w + (x+y+z)*(a+b+c)/3	*/
			t3 = t3-t5;						t4 = t4-t6;	/* x-z	*/
			t5 = t7-t5;						t6 = t8-t6;	/* y-z	*/
			t7 =(t3+t5)*cc3;				t8 =(t4+t6)*cc3;	/* (x+y-2z)*(a+b-2c)/3	*/
	#if PFETCH
	addr = add0+p3;
	prefetch_p_doubles(addr);
	#endif
			t3 = t3*cc1;					t4 = t4*cc1;	/* (x-z)*(a-c)	*/
			t5 = t5*cc2;					t6 = t6*cc2;	/* (y-z)*(b-c)	*/
			rt = t3-t7;						it = t4-t8;	/* x*(2a-b-c)/3 + y*(-a-b+2c)/3 + z*(-a+2b-c)/3	*/
			t5 = t5-t7;						t6 = t6-t8;	/* x*(-a-b+2c)/3+ y*(-a+2b-c)/3 + z*(2a-b-c)/3	*/
											/* v0+v1 = x*(a-2b+c)/3 + y*(-2a+b+c)/3 + z*(a+b-2c)/3	*/
			t3 = t1-rt-t5;					t4 = t2-it-t6;	/* C2 = x*b + y*a + z*c	*/
			t5 = t1+t5;						t6 = t2+t6;	/* C3 = x*c + y*b + z*a	*/
			t1 = t1+rt;						t2 = t2+it;	/* C1 = x*a + y*c + z*b	*/

	/*  Get cyclic convolution of (ss1,ss2,-ss3) with (t13,-t9,t11) and with (t14,-t10,t12),
	!   (note permutation of t-elements and sign flips on ss3, t9 and t10): */
			t7 =(t13-t9 +t11)*ss0;			t8 =(t14-t10+t12)*ss0;
			t13= t13-t11;					t14= t14-t12;
			t11= t9 +t11;					t12= t10+t12;		/* Flip sign here... */
			t9 =(t11-t13)*ss3;				t10=(t12-t14)*ss3;	/* ...and here... */
	#if PFETCH
	addr = add0+p4;
	prefetch_p_doubles(addr);
	#endif
			t13= t13*ss1;					t14= t14*ss1;
			t11= t11*ss2;					t12= t12*ss2;
			t13= t9+t13;					t14= t10+t14;		/* Sign flips get undone in these next 2. */
			t11= t9-t11;					t12= t10-t12;

			t9 = t7-t13-t11;				t10= t8-t14-t12;	/* S2 */
			t11= t7+t11;					t12= t8+t12;		/*-S3 */
			t7 = t7+t13;					t8 = t8+t14;		/* S1 */
	#if PFETCH
	addr = add0+p5;
	prefetch_p_doubles(addr);
	#endif
			aj1p2r =t1-t8;					aj1p2i =t2+t7;		/* X1 = C1 + I*S1	*/
			aj1p4r =t3-t10;					aj1p4i =t4+t9;		/* X2 = C2 + I*S2	*/
			aj1p6r =t5+t12;					aj1p6i =t6-t11;		/* X3 = C3 + I*S3	*/
			aj1p8r =t5-t12;					aj1p8i =t6+t11;		/* X4 =	C3 - I*S3	*/
			aj1p10r=t3+t10;					aj1p10i=t4-t9;		/* X5 =	C2 - I*S2	*/
			aj1p12r=t1+t8;					aj1p12i=t2-t7;		/* X6 =	C1 - I*S1	*/
	#endif
	/*...Second radix-7 block uses aj1p[1:13:2] as inputs.	*/
	#if PFETCH
	addr = add0+p6;
	prefetch_p_doubles(addr);
	#endif
			t1 =aj1p7r ;					t2 =aj1p7i ;
			t3 =aj1p5r +aj1p9r ;			t4 =aj1p5i +aj1p9i ;
			t5 =aj1p3r +aj1p11r;			t6 =aj1p3i +aj1p11i;
			t7 =aj1p1r +aj1p13r;			t8 =aj1p1i +aj1p13i;
	#if PFETCH
	addr = add0+p7;
	prefetch_p_doubles(addr);
	#endif
			t9 =aj1p1r -aj1p13r;			t10=aj1p1i -aj1p13i;
			t11=aj1p3r -aj1p11r;			t12=aj1p3i -aj1p11i;
			t13=aj1p5r -aj1p9r ;			t14=aj1p5i -aj1p9i ;
	#if PFETCH
	addr = add0+p8;
	prefetch_p_doubles(addr);
	#endif
		/* Version 1: direct computation of 3x3 matrix-vector product. Length-3 convo cost: 36 ADD, 36 MUL */
	#if LO_ADD
			aj1p1r = t1+t3+t5+t7;			aj1p1i = t2+t4+t6+t8	;
			rt = t1+cc1*t3+cc2*t5+cc3*t7;	it = t2+cc1*t4+cc2*t6+cc3*t8;
			re = t1+cc2*t3+cc3*t5+cc1*t7;	im = t2+cc2*t4+cc3*t6+cc1*t8;
			t1 = t1+cc3*t3+cc1*t5+cc2*t7;	t2 = t2+cc3*t4+cc1*t6+cc2*t8;
	#if PFETCH
	addr = add0+p9;
	prefetch_p_doubles(addr);
	#endif
			t3 = ss1*t13+ss2*t11+ss3*t9;	t4 = ss1*t14+ss2*t12+ss3*t10;
			t5 = ss2*t13-ss3*t11-ss1*t9;	t6 = ss2*t14-ss3*t12-ss1*t10;
			t7 = ss3*t13-ss1*t11+ss2*t9;	t8 = ss3*t14-ss1*t12+ss2*t10;
	#if PFETCH
	addr = add0+p10;
	prefetch_p_doubles(addr);
	#endif
		/* Output permutation causes signs to get flipped here: */
			aj1p3r =rt-t4;					aj1p3i =it+t3;
			aj1p5r =re-t6;					aj1p5i =im+t5;
			aj1p7r =t1-t8;					aj1p7i =t2+t7;
	#if PFETCH
	addr = add0+p11;
	prefetch_p_doubles(addr);
	#endif
			aj1p9r =t1+t8;					aj1p9i =t2-t7;
			aj1p11r=re+t6;					aj1p11i=im-t5;
			aj1p13r=rt+t4;					aj1p13i=it-t3;

		/* Version 2: get cyclic convolution of (cc1,cc2,cc3) with (t3,t7,t5) and with (t4,t8,t6), (note switch of
		last 2 t-elements in each case), add t1 and t2 to respective outputs. Length-3 convo cost: 48 ADD, 16 MUL */
	#else
			rt = t3+t7+t5; aj1p1r=rt+t1;	it = t4+t8+t6; aj1p1i=it+t2;
			t1 = rt*cc0+t1;					t2 = it*cc0+t2;
			t3 = t3-t5;						t4 = t4-t6;
			t5 = t7-t5;						t6 = t8-t6;
			t7 =(t3+t5)*cc3;				t8 =(t4+t6)*cc3;
	#if PFETCH
	addr = add0+p9;
	prefetch_p_doubles(addr);
	#endif
			t3 = t3*cc1;					t4 = t4*cc1;
			t5 = t5*cc2;					t6 = t6*cc2;
			rt = t3-t7;						it = t4-t8;
			t5 = t5-t7;						t6 = t6-t8;

			t3 = t1-rt-t5;					t4 = t2-it-t6;		/* C2	*/
			t5 = t1+t5;						t6 = t2+t6;			/* C3	*/
			t1 = t1+rt;						t2 = t2+it;			/* C1	*/

	/*  Get cyclic convolution of (ss1,ss2,-ss3) with (t13,-t9,t11) and with (t14,-t10,t12),
	!   (note permutation of t-elements and sign flips on ss3, t9 and t10): */
			t7 =(t13-t9 +t11)*ss0;			t8 =(t14-t10+t12)*ss0;
			t13= t13-t11;					t14= t14-t12;
			t11= t9 +t11;					t12= t10+t12;		/* Flip sign here... */
			t9 =(t11-t13)*ss3;				t10=(t12-t14)*ss3;	/* ...and here... */
	#if PFETCH
	addr = add0+p10;
	prefetch_p_doubles(addr);
	#endif
			t13= t13*ss1;					t14= t14*ss1;
			t11= t11*ss2;					t12= t12*ss2;
			t13= t9+t13;					t14= t10+t14;		/* Sign flips get undone in these next 2. */
			t11= t9-t11;					t12= t10-t12;

			t9 = t7-t13-t11;				t10= t8-t14-t12;	/* S2 */
			t11= t7+t11;					t12= t8+t12;		/*-S3 */
			t7 = t7+t13;					t8 = t8+t14;		/* S1 */
	#if PFETCH
	addr = add0+p11;
	prefetch_p_doubles(addr);
	#endif
			aj1p3r =t1-t8;					aj1p3i =t2+t7;		/* Y1 = C1 + I*S1	*/
			aj1p5r =t3-t10;					aj1p5i =t4+t9;		/* Y2 = C2 + I*S2	*/
			aj1p7r =t5+t12;					aj1p7i =t6-t11;		/* Y3 = C3 + I*S3	*/
			aj1p9r =t5-t12;					aj1p9i =t6+t11;		/* Y4 =	C3 - I*S3	*/
			aj1p11r=t3+t10;					aj1p11i=t4-t9;		/* Y5 =	C2 - I*S2	*/
			aj1p13r=t1+t8;					aj1p13i=t2-t7;		/* Y6 =	C1 - I*S1	*/
	#if PFETCH
	addr = add0+p11;
	prefetch_p_doubles(addr);
	#endif
	#endif
	/*...and now do seven radix-2 transforms:	*/

			a[j1    ]=aj1p0r +aj1p1r;		a[j2    ]=aj1p0i +aj1p1i;
			a[j1+p1 ]=aj1p0r -aj1p1r;		a[j2+p1 ]=aj1p0i -aj1p1i;

			a[j1+p12]=aj1p2r +aj1p3r;		a[j2+p12]=aj1p2i +aj1p3i;
			a[j1+p13]=aj1p2r -aj1p3r;		a[j2+p13]=aj1p2i -aj1p3i;

			a[j1+p11]=aj1p4r +aj1p5r;		a[j2+p11]=aj1p4i +aj1p5i;
			a[j1+p10]=aj1p4r -aj1p5r;		a[j2+p10]=aj1p4i -aj1p5i;
	#if PFETCH
	addr = add0+p12;
	prefetch_p_doubles(addr);
	#endif
			a[j1+p8 ]=aj1p6r +aj1p7r;		a[j2+p8 ]=aj1p6i +aj1p7i;
			a[j1+p9 ]=aj1p6r -aj1p7r;		a[j2+p9 ]=aj1p6i -aj1p7i;

			a[j1+p7 ]=aj1p8r +aj1p9r;		a[j2+p7 ]=aj1p8i +aj1p9i;
			a[j1+p6 ]=aj1p8r -aj1p9r;		a[j2+p6 ]=aj1p8i -aj1p9i;

			a[j1+p4 ]=aj1p10r+aj1p11r;		a[j2+p4 ]=aj1p10i+aj1p11i;
			a[j1+p5 ]=aj1p10r-aj1p11r;		a[j2+p5 ]=aj1p10i-aj1p11i;

			a[j1+p3 ]=aj1p12r+aj1p13r;		a[j2+p3 ]=aj1p12i+aj1p13i;
			a[j1+p2 ]=aj1p12r-aj1p13r;		a[j2+p2 ]=aj1p12i-aj1p13i;
	#if PFETCH
	addr = add0+p13;
	prefetch_p_doubles(addr);
	#endif
		}

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			jstart += nwt;
			jhi    += nwt;

			col += 14;
			co3 -= 14;
		}
	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-14 forward DIF FFT of the first block of 14 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 14 outputs of (1);
!   (3) Reweight and perform a radix-14 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 14 elements and repeat (1-4).
*/
/*
printf("iter = %10d; maxerr = %20.15f\n",iter,maxerr);
printf("carries = %10d %10d %10d %10d %10d %10d %10d\n",(int)cy_r0,(int)cy_r1,(int)cy_r2,(int)cy_r3,(int)cy_r4,(int)cy_r5,(int)cy_r6);
*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t1    = cy_r13;
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
		t1    = cy_r13;	t2    = cy_i13;
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

	root_incr = 0;
	scale = prp_mult = 1;

	jstart = 0;
	/*
	For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
	*/
	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		jhi = 7;
	}
	else
	{
		jhi = 3;
	}
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
		a[j+p13] *= radix_inv;
	}
}

	if(fabs(cy_r0)+fabs(cy_r1)+fabs(cy_r2)+fabs(cy_r3)+fabs(cy_r4)+fabs(cy_r5)+fabs(cy_r6)+fabs(cy_r7)+fabs(cy_r8)+fabs(cy_r9)+fabs(cy_r10)+fabs(cy_r11)+fabs(cy_r12)+fabs(cy_r13)
		+fabs(cy_i0)+fabs(cy_i1)+fabs(cy_i2)+fabs(cy_i3)+fabs(cy_i4)+fabs(cy_i5)+fabs(cy_i6)+fabs(cy_i7)+fabs(cy_i8)+fabs(cy_i9)+fabs(cy_i10)+fabs(cy_i11)+fabs(cy_i12)+fabs(cy_i13) != 0.0)
	{
		sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in radix14_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
		mlucas_fprint(cbuf,INTERACT);
		err=ERR_CARRY;
		return(err);
	}
	*fracmax = maxerr;
	return(0);
}

/***************/

void radix14_dif_pass1(double a[], int n)
{
/*
!
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-14 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
*/
	int j,j1,j2;
	static int n14,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13, first_entry=TRUE;
#if LO_ADD
	static double cc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					ss1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					cc2 =-.22252093395631440426,	 /* cos(2u)	*/
					ss2 = .97492791218182360702,	 /* sin(2u)	*/
					cc3 =-.90096886790241912622,	 /* cos(3u)	*/
					ss3 = .43388373911755812050;	 /* sin(3u)	*/
					/*
					Radix-14-specific roots are all symmetric with Radix-7-specific ones.
					cs  = .90096886790241912623 = -cc3	Real part of exp(i*2*pi/14), the radix-14 fundamental sincos datum.
					sn  = .43388373911755812047 = +ss3	Imag part of exp(i*2*pi/14).
					cs3 = .22252093395631440430 = -cc2
					sn3 = .97492791218182360701 = +ss2
					*/
	double re,im;
#else
	static double c0  =-0.16666666666666666667,	/* (cc1+cc2+cc3)/3	*/
					c1  = 1.52445866976115265675, 	/*  cc1-cc3		*/
					c2  = 0.67844793394610472196, 	/*  cc2-cc3		*/
					c3  = 0.73430220123575245957,	/* (cc1+cc2-2*cc3)/3	*/
		/* Switch the sign of ss3 in these: */
					s0  = 0.44095855184409843174,	/* (ss1+ss2-ss3)/3	*/
					s1  = 1.21571522158558792920, 	/*  ss1+ss3		*/
					s2  = 1.40881165129938172752, 	/*  ss2+ss3		*/
					s3  = 0.87484229096165655224;	/* (ss1+ss2+2*ss3)/3	*/
#endif
	double rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14
	,aj1p0r,aj1p1r,aj1p2r,aj1p3r,aj1p4r,aj1p5r,aj1p6r,aj1p7r,aj1p8r,aj1p9r,aj1p10r,aj1p11r,aj1p12r,aj1p13r
	,aj1p0i,aj1p1i,aj1p2i,aj1p3i,aj1p4i,aj1p5i,aj1p6i,aj1p7i,aj1p8i,aj1p9i,aj1p10i,aj1p11i,aj1p12i,aj1p13i;

	if(!first_entry && (n/14) != n14)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n14=n/14;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n14;
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
	}

/*...The radix-14 pass is here.	*/

      for(j=0; j < n14; j += 2)
      {
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

/*       gather the needed data (14 64-bit complex, i.e. 28 64-bit reals) and do two radix-7 transforms...	*/

/*
Twiddleless version requires us to swap inputs as follows:
indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13
      -> 0,-7,12, 5,10, 3, 8, 1, 6,-1, 4,-3, 2,-5
      == 0, 7,12, 5,10, 3, 8, 1, 6,13, 4,11, 2, 9 modulo 14.
I.e. start out with first septet of indices {0,2,4,6,8,10,12}, permute those according to
{0,2,4,6,8,10,12}*13%14 = {0,12,10,8,6,4,2,}, then each is head of a length-2 list of indices with decrement 7.
*/

/*...First radix-7 block uses aj1p[0:12:2] as inputs:	*/

		t1 =a[j1    ];					t2 =a[j2    ];
		t3 =a[j1+p12]+a[j1+p2 ];		t4 =a[j2+p12]+a[j2+p2 ];	/* x1+x6 */
		t5 =a[j1+p10]+a[j1+p4 ];		t6 =a[j2+p10]+a[j2+p4 ];	/* x2+x5 */
		t7 =a[j1+p8 ]+a[j1+p6 ];		t8 =a[j2+p8 ]+a[j2+p6 ];	/* x3+x4 */
		t9 =a[j1+p8 ]-a[j1+p6 ];		t10=a[j2+p8 ]-a[j2+p6 ];	/* x3-x4 */
		t11=a[j1+p10]-a[j1+p4 ];		t12=a[j2+p10]-a[j2+p4 ];	/* x2-x5 */
		t13=a[j1+p12]-a[j1+p2 ];		t14=a[j2+p12]-a[j2+p2 ];	/* x1-x6 */

	/* Version 1: direct computation of 3x3 matrix-vector product. */
#if LO_ADD
		aj1p0r = t1+t3+t5+t7;			aj1p0i = t2+t4+t6+t8	;
		rt = t1+cc1*t3+cc2*t5+cc3*t7;	it = t2+cc1*t4+cc2*t6+cc3*t8;
		re = t1+cc2*t3+cc3*t5+cc1*t7;	im = t2+cc2*t4+cc3*t6+cc1*t8;
		t1 = t1+cc3*t3+cc1*t5+cc2*t7;	t2 = t2+cc3*t4+cc1*t6+cc2*t8;

		t3 = ss1*t13+ss2*t11+ss3*t9;	t4 = ss1*t14+ss2*t12+ss3*t10;
		t5 = ss2*t13-ss3*t11-ss1*t9;	t6 = ss2*t14-ss3*t12-ss1*t10;
		t7 = ss3*t13-ss1*t11+ss2*t9;	t8 = ss3*t14-ss1*t12+ss2*t10;
	/* Output permutation causes signs to get flipped here: */
		aj1p2r =rt-t4;					aj1p2i =it+t3;
		aj1p4r =re-t6;					aj1p4i =im+t5;
		aj1p6r =t1-t8;					aj1p6i =t2+t7;
		aj1p8r =t1+t8;					aj1p8i =t2-t7;
		aj1p10r=re+t6;					aj1p10i=im-t5;
		aj1p12r=rt+t4;					aj1p12i=it-t3;

	/* Version 2: get cyclic convolution of (cc1,cc2,cc3) with (t3,t7,t5) and with (t4,t8,t6), (note switch of
	last 2 t-elements in each case), add t1 and t2 to respective outputs. */
#else
		rt = t3+t7+t5; aj1p0r=rt+t1;	it = t4+t8+t6; aj1p0i=it+t2;
		t1 = rt*c0+t1;					t2 = it*c0+t2;
		t3 = t3-t5;						t4 = t4-t6;
		t5 = t7-t5;						t6 = t8-t6;
		t7 =(t3+t5)*c3;					t8 =(t4+t6)*c3;
		t3 = t3*c1;						t4 = t4*c1;
		t5 = t5*c2;						t6 = t6*c2;
		rt = t3-t7;						it = t4-t8;
		t5 = t5-t7;						t6 = t6-t8;

		t3 = t1-rt-t5;					t4 = t2-it-t6;		/* C2	*/
		t5 = t1+t5;						t6 = t2+t6;		/* C3	*/
		t1 = t1+rt;						t2 = t2+it;		/* C1	*/

	/* Get cyclic convolution of (ss1,ss2,-ss3) with (t13,-t9,t11) and with (t14,-t10,t12),
		 (note permutation of t-elements and sign flips on ss3, t9 and t10): */

		t7 =(t13-t9 +t11)*s0;			t8 =(t14-t10+t12)*s0;
		t13= t13-t11;					t14= t14-t12;
		t11= t9 +t11;					t12= t10+t12;		/* Flip sign here... */
		t9 =(t11-t13)*s3;				t10=(t12-t14)*s3;	/* ...and here... */
		t13= t13*s1;					t14= t14*s1;
		t11= t11*s2;					t12= t12*s2;
		t13= t9+t13;					t14= t10+t14;		/* Sign flips get undone in these next 2. */
		t11= t9-t11;					t12= t10-t12;

		t9 = t7-t13-t11;				t10= t8-t14-t12;	/* S2 */
		t11= t7+t11;					t12= t8+t12;		/*-S3 */
		t7 = t7+t13;					t8 = t8+t14;		/* S1 */
	/* Output permutation causes signs to get flipped here: */
		aj1p2r =t1-t8;					aj1p2i =t2+t7;		/* X1 = C1 + I*S1	*/
		aj1p4r =t3-t10;					aj1p4i =t4+t9;		/* X2 = C2 + I*S2	*/
		aj1p6r =t5+t12;					aj1p6i =t6-t11;		/* X3 = C3 + I*S3	*/
		aj1p8r =t5-t12;					aj1p8i =t6+t11;		/* X4 =	C3 - I*S3	*/
		aj1p10r=t3+t10;					aj1p10i=t4-t9;		/* X5 =	C2 - I*S2	*/
		aj1p12r=t1+t8;					aj1p12i=t2-t7;		/* X6 =	C1 - I*S1	*/
#endif
/*
printf("%20.15f  %20.15f\n", aj1p0r , aj1p0i );
printf("%20.15f  %20.15f\n", aj1p2r , aj1p2i );
printf("%20.15f  %20.15f\n", aj1p4r , aj1p4i );
printf("%20.15f  %20.15f\n", aj1p6r , aj1p6i );
printf("%20.15f  %20.15f\n", aj1p8r , aj1p8i );
printf("%20.15f  %20.15f\n", aj1p10r, aj1p10i);
printf("%20.15f  %20.15f\n", aj1p12r, aj1p12i);
*/
/*...Second radix-7 block uses aj1p[1:13:2] as inputs.	*/

		t1 =a[j1+p7 ];					t2 =a[j2+p7 ];
		t3 =a[j1+p5 ]+a[j1+p9 ];		t4 =a[j2+p5 ]+a[j2+p9 ];
		t5 =a[j1+p3 ]+a[j1+p11];		t6 =a[j2+p3 ]+a[j2+p11];
		t7 =a[j1+p1 ]+a[j1+p13];		t8 =a[j2+p1 ]+a[j2+p13];
		t9 =a[j1+p1 ]-a[j1+p13];		t10=a[j2+p1 ]-a[j2+p13];
		t11=a[j1+p3 ]-a[j1+p11];		t12=a[j2+p3 ]-a[j2+p11];
		t13=a[j1+p5 ]-a[j1+p9 ];		t14=a[j2+p5 ]-a[j2+p9 ];

	/* Version 1: direct computation of 3x3 matrix-vector product. Length-3 convo cost: 36 ADD, 36 MUL */
#if LO_ADD
		aj1p1r = t1+t3+t5+t7;			aj1p1i = t2+t4+t6+t8	;
		rt = t1+cc1*t3+cc2*t5+cc3*t7;	it = t2+cc1*t4+cc2*t6+cc3*t8;
		re = t1+cc2*t3+cc3*t5+cc1*t7;	im = t2+cc2*t4+cc3*t6+cc1*t8;
		t1 = t1+cc3*t3+cc1*t5+cc2*t7;	t2 = t2+cc3*t4+cc1*t6+cc2*t8;

		t3 = ss1*t13+ss2*t11+ss3*t9;	t4 = ss1*t14+ss2*t12+ss3*t10;
		t5 = ss2*t13-ss3*t11-ss1*t9;	t6 = ss2*t14-ss3*t12-ss1*t10;
		t7 = ss3*t13-ss1*t11+ss2*t9;	t8 = ss3*t14-ss1*t12+ss2*t10;
	/* Output permutation causes signs to get flipped here: */
		aj1p3r =rt-t4;					aj1p3i =it+t3;
		aj1p5r =re-t6;					aj1p5i =im+t5;
		aj1p7r =t1-t8;					aj1p7i =t2+t7;
		aj1p9r =t1+t8;					aj1p9i =t2-t7;
		aj1p11r=re+t6;					aj1p11i=im-t5;
		aj1p13r=rt+t4;					aj1p13i=it-t3;

	/* Version 2: get cyclic convolution of (cc1,cc2,cc3) with (t3,t7,t5) and with (t4,t8,t6), (note switch of
	last 2 t-elements in each case), add t1 and t2 to respective outputs. Length-3 convo cost: 48 ADD, 16 MUL */
#else
		rt = t3+t7+t5; aj1p1r=rt+t1;	it = t4+t8+t6; aj1p1i=it+t2;
		t1 = rt*c0+t1;					t2 = it*c0+t2;
		t3 = t3-t5;						t4 = t4-t6;
		t5 = t7-t5;						t6 = t8-t6;
		t7 =(t3+t5)*c3;					t8 =(t4+t6)*c3;
		t3 = t3*c1;						t4 = t4*c1;
		t5 = t5*c2;						t6 = t6*c2;
		rt = t3-t7;						it = t4-t8;
		t5 = t5-t7;						t6 = t6-t8;

		t3 = t1-rt-t5;					t4 = t2-it-t6;		/* C2	*/
		t5 = t1+t5;						t6 = t2+t6;		/* C3	*/
		t1 = t1+rt;						t2 = t2+it;		/* C1	*/

	/* Get cyclic convolution of (ss1,ss2,-ss3) with (t13,-t9,t11) and with (t14,-t10,t12),
		 (note permutation of t-elements and sign flips on ss3, t9 and t10): */

		t7 =(t13-t9 +t11)*s0;			t8 =(t14-t10+t12)*s0;
		t13= t13-t11;					t14= t14-t12;
		t11= t9 +t11;					t12= t10+t12;		/* Flip sign here... */
		t9 =(t11-t13)*s3;				t10=(t12-t14)*s3;	/* ...and here... */
		t13= t13*s1;					t14= t14*s1;
		t11= t11*s2;					t12= t12*s2;
		t13= t9+t13;					t14= t10+t14;		/* Sign flips get undone in these next 2. */
		t11= t9-t11;					t12= t10-t12;

		t9 = t7-t13-t11;				t10= t8-t14-t12;	/* S2 */
		t11= t7+t11;					t12= t8+t12;		/*-S3 */
		t7 = t7+t13;					t8 = t8+t14;		/* S1 */
	/* Output permutation causes signs to get flipped here: */
		aj1p3r =t1-t8;					aj1p3i =t2+t7;		/* Y1 = C1 + I*S1	*/
		aj1p5r =t3-t10;					aj1p5i =t4+t9;		/* Y2 = C2 + I*S2	*/
		aj1p7r =t5+t12;					aj1p7i =t6-t11;		/* Y3 = C3 + I*S3	*/
		aj1p9r =t5-t12;					aj1p9i =t6+t11;		/* Y4 =	C3 - I*S3	*/
		aj1p11r=t3+t10;					aj1p11i=t4-t9;		/* Y5 =	C2 - I*S2	*/
		aj1p13r=t1+t8;					aj1p13i=t2-t7;		/* Y6 =	C1 - I*S1	*/
#endif
/*
printf("%20.15f  %20.15f\n", aj1p1r , aj1p1i );
printf("%20.15f  %20.15f\n", aj1p3r , aj1p3i );
printf("%20.15f  %20.15f\n", aj1p5r , aj1p5i );
printf("%20.15f  %20.15f\n", aj1p7r , aj1p7i );
printf("%20.15f  %20.15f\n", aj1p9r , aj1p9i );
printf("%20.15f  %20.15f\n", aj1p11r, aj1p11i);
printf("%20.15f  %20.15f\n", aj1p13r, aj1p13i);
*/
/*...and now do seven radix-2 transforms:	*/

		a[j1    ]=aj1p0r +aj1p1r;		a[j2    ]=aj1p0i +aj1p1i;
		a[j1+p1 ]=aj1p0r -aj1p1r;		a[j2+p1 ]=aj1p0i -aj1p1i;

		a[j1+p12]=aj1p2r +aj1p3r;		a[j2+p12]=aj1p2i +aj1p3i;
		a[j1+p13]=aj1p2r -aj1p3r;		a[j2+p13]=aj1p2i -aj1p3i;

		a[j1+p11]=aj1p4r +aj1p5r;		a[j2+p11]=aj1p4i +aj1p5i;
		a[j1+p10]=aj1p4r -aj1p5r;		a[j2+p10]=aj1p4i -aj1p5i;

		a[j1+p8 ]=aj1p6r +aj1p7r;		a[j2+p8 ]=aj1p6i +aj1p7i;
		a[j1+p9 ]=aj1p6r -aj1p7r;		a[j2+p9 ]=aj1p6i -aj1p7i;

		a[j1+p7 ]=aj1p8r +aj1p9r;		a[j2+p7 ]=aj1p8i +aj1p9i;
		a[j1+p6 ]=aj1p8r -aj1p9r;		a[j2+p6 ]=aj1p8i -aj1p9i;

		a[j1+p4 ]=aj1p10r+aj1p11r;		a[j2+p4 ]=aj1p10i+aj1p11i;
		a[j1+p5 ]=aj1p10r-aj1p11r;		a[j2+p5 ]=aj1p10i-aj1p11i;

		a[j1+p3 ]=aj1p12r+aj1p13r;		a[j2+p3 ]=aj1p12i+aj1p13i;
		a[j1+p2 ]=aj1p12r-aj1p13r;		a[j2+p2 ]=aj1p12i-aj1p13i;
								/* Version 1 Totals: 120+7*4 = 148 FADD, 72 FMUL,	*/
								/* 24 fewer adds and 3/4 the multiplies of the naive version. */
								/* Version 2 Totals: 144+7*4 = 172 FADD, 32 FMUL,	*/
								/* same #adds, just 1/3 the multiplies of the naive version. */
	}
}

/***************/

void radix14_dit_pass1(double a[], int n)
{
/*
!
!...Acronym: DIT = Decimation In Time
!
!...Subroutine to perform an initial radix-14 complex DIT FFT pass on the data in the length-N real vector A.
!
!   This routine is designed exclusively to undo the effects of radix14_dif_pass1,
!   i.e. to reobtain the raw all-integer residue vector at the end of an iteration cycle.
!
!   See the documentation in radix16_dif_pass for further details on the array indexing,
!   and radix14_dif_pass1 for details on the algorithm.
!
!   See the documentation in radix7_dif_pass1 for notes on the algorithm.
*/
	int j,j1,j2;
	static int n14,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13, first_entry=TRUE;
#if LO_ADD
	static double cc1 = .62348980185873353053,	 /* cos(u) = Real part of exp(i*2*pi/7), the radix-7 fundamental sincos datum	*/
					ss1 = .78183148246802980870,	 /* sin(u) = Imag part of exp(i*2*pi/7).	*/
					cc2 =-.22252093395631440426,	 /* cos(2u)	*/
					ss2 = .97492791218182360702,	 /* sin(2u)	*/
					cc3 =-.90096886790241912622,	 /* cos(3u)	*/
					ss3 = .43388373911755812050;	 /* sin(3u)	*/
					/*
					Radix-14-specific roots are all symmetric with Radix-7-specific ones.
					cs  = .90096886790241912623 = -cc3	Real part of exp(i*2*pi/14), the radix-14 fundamental sincos datum.
					sn  = .43388373911755812047 = +ss3	Imag part of exp(i*2*pi/14).
					cs3 = .22252093395631440430 = -cc2
					sn3 = .97492791218182360701 = +ss2
					*/
	double re,im;
#else
	static double c0  =-0.16666666666666666667,	/* (cc1+cc2+cc3)/3	*/
					c1  = 1.52445866976115265675, 	/*  cc1-cc3		*/
					c2  = 0.67844793394610472196, 	/*  cc2-cc3		*/
					c3  = 0.73430220123575245957,	/* (cc1+cc2-2*cc3)/3	*/
		/* Switch the sign of ss3 in these: */
					s0  = 0.44095855184409843174,	/* (ss1+ss2-ss3)/3	*/
					s1  = 1.21571522158558792920, 	/*  ss1+ss3		*/
					s2  = 1.40881165129938172752, 	/*  ss2+ss3		*/
					s3  = 0.87484229096165655224;	/* (ss1+ss2+2*ss3)/3	*/
#endif
	double rt,it
	,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28;

	if(!first_entry && (n/14) != n14)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n14=n/14;

/*   constant index offsets for array load/stores are here.	*/

		p1 = n14;
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
	}

/*...The radix-14 pass is here.	*/

	for(j=0; j < n14; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*...gather the needed data (14 64-bit complex, i.e. 28 64-bit reals) and do seven radix-2 transforms,	*/
	/*
	Twiddleless version requires us to swap inputs as follows:
	indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13
		  -> 0,-2,-4,-6,-8,-A,-C, 7, 5, 3, 1,-1,-3,-5
		  == 0,12,10, 8, 6, 4, 2, 7, 5, 3, 1,13,11, 9 modulo 14.
	I.e. start out with first pair of indices {0,7}, permute those according to
	{0,7}*13%14 = {0,7}, then each is head of a length-7 list of indices with decrement 2.

	Remember, inputs to DIT are bit-reversed, so
	a[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13] contain
	x[0, 7, 1, 8, 2, 9, 3,10, 4,11, 5,12, 6,13], which get swapped to
	x[0, 7,12, 5,10, 3, 8, 1, 6,13, 4,11, 2, 9], which means the a-indices get swapped as
	a[0, 1,11,10, 7, 6, 3, 2,12,13, 8, 9, 4, 5].
	*/
		t1 =a[j1    ];					t2 =a[j2    ];
		rt =a[j1+p1 ];					it =a[j2+p1 ];
		t3 =t1 -rt;						t4 =t2 -it;
		t1 =t1 +rt;						t2 =t2 +it;

		t5 =a[j1+p11];					t6 =a[j2+p11];
		rt =a[j1+p10];					it =a[j2+p10];
		t7 =t5 -rt;						t8 =t6 -it;
		t5 =t5 +rt;						t6 =t6 +it;

		t9 =a[j1+p7 ];					t10=a[j2+p7 ];
		rt =a[j1+p6 ];					it =a[j2+p6 ];
		t11=t9 -rt;						t12=t10-it;
		t9 =t9 +rt;						t10=t10+it;

		t13=a[j1+p3 ];					t14=a[j2+p3 ];
		rt =a[j1+p2 ];					it =a[j2+p2 ];
		t15=t13-rt;						t16=t14-it;
		t13=t13+rt;						t14=t14+it;

		t17=a[j1+p12];					t18=a[j2+p12];
		rt =a[j1+p13];					it =a[j2+p13];
		t19=t17-rt;						t20=t18-it;
		t17=t17+rt;						t18=t18+it;

		t21=a[j1+p8 ];					t22=a[j2+p8 ];
		rt =a[j1+p9 ];					it =a[j2+p9 ];
		t23=t21-rt;						t24=t22-it;
		t21=t21+rt;						t22=t22+it;

		t25=a[j1+p4 ];					t26=a[j2+p4 ];
		rt =a[j1+p5 ];					it =a[j2+p5 ];
		t27=t25-rt;						t28=t26-it;
		t25=t25+rt;						t26=t26+it;

/*       ...and now do two radix-7 transforms.	*/

/*...aj1p[0:12:2]r use t[1:25:4]; aj1p[0:12:2]i use t[2:26:4]	*/

		rt =t25;						it =t26;
		t25=t5 -rt;						t26=t6 -it;
		t5 =t5 +rt;						t6 =t6 +it;
		rt =t21;						it =t22;
		t21=t9 -rt;						t22=t10 -it;
		t9 =t9 +rt;						t10=t10 +it;
		rt =t17;						it =t18;
		t17=t13 -rt;					t18=t14 -it;
		t13=t13 +rt;					t14=t14 +it;

	/* Version 1: direct computation of 3x3 matrix-vector product. */
#if LO_ADD
		a[j1] = t1+t5+t9+t13;			a[j2] = t2+t6+t10+t14;
		rt = t1+cc1*t5+cc2*t9+cc3*t13;	it = t2+cc1*t6+cc2*t10+cc3*t14;
		re = t1+cc2*t5+cc3*t9+cc1*t13;	im = t2+cc2*t6+cc3*t10+cc1*t14;
		t1 = t1+cc3*t5+cc1*t9+cc2*t13;	t2 = t2+cc3*t6+cc1*t10+cc2*t14;

		t5 = ss1*t25+ss2*t21+ss3*t17;	t6 = ss1*t26+ss2*t22+ss3*t18;
		t9 = ss2*t25-ss3*t21-ss1*t17;	t10= ss2*t26-ss3*t22-ss1*t18;
		t13= ss3*t25-ss1*t21+ss2*t17;	t14= ss3*t26-ss1*t22+ss2*t18;
	/* Output permutation causes signs to get flipped here: */
		a[j1+p6 ]=rt+t6;				a[j2+p6 ]=it-t5;
		a[j1+p12]=re+t10;				a[j2+p12]=im-t9;
		a[j1+p4 ]=t1+t14;				a[j2+p4 ]=t2-t13;
		a[j1+p10]=t1-t14;				a[j2+p10]=t2+t13;
		a[j1+p2 ]=re-t10;				a[j2+p2 ]=im+t9;
		a[j1+p8 ]=rt-t6;				a[j2+p8 ]=it+t5;

	/* Version 2: */
#else
	/* Get cyclic convolution of (cc1,cc2,cc3) with (t5,t13,t9) and with (t6,t14,t10),
		 (note switch of last 2 t-elements in each case), add t1 and t2 to respective outputs: */

		rt = t5+t13+t9; a[j1  ]=rt+t1;	it = t6+t14+t10; a[j2]=it+t2;
		t1 = rt*c0+t1;					t2 = it*c0+t2;
		t5 = t5-t9;						t6 = t6 -t10;
		t9 = t13-t9;					t10= t14-t10;
		t13=(t5+t9)*c3;					t14=(t6+t10)*c3;
		t5 = t5*c1;						t6 = t6 *c1;
		t9 = t9*c2;						t10= t10*c2;
		rt = t5-t13;					it = t6 -t14;
		t9 = t9-t13;					t10= t10-t14;

		t5 = t1-rt-t9;					t6 = t2-it-t10;	/* C2	*/
		t9 = t1+t9;						t10= t2+t10;	/* C3	*/
		t1 = t1+rt;						t2 = t2+it;	/* C1	*/

	/* Get cyclic convolution of (ss1,ss2,-ss3) with (t25,-t17,t21) and with (t26,-t18,t22),
		 (note permutation of t-elements and sign flips on ss3, t17 and t18): */

		t13=(t25-t17 +t21)*s0;			t14=(t26-t18+t22)*s0;
		t25= t25-t21;					t26= t26-t22;
		t21= t17+t21;					t22= t18+t22;		/* Flip sign here... */
		t17=(t21-t25)*s3;				t18=(t22-t26)*s3;	/* ...and here... */
		t25= t25*s1;					t26= t26*s1;
		t21= t21*s2;					t22= t22*s2;
		t25= t17+t25;					t26= t18+t26;		/* Sign flips get undone in these next 2. */
		t21= t17-t21;					t22= t18-t22;

		t17= t13-t25-t21;				t18= t14-t26-t22;	/* S2 */
		t21= t13+t21;					t22= t14+t22;		/*-S3 */
		t13= t13+t25;					t14= t14+t26;		/* S1 */

		a[j1+p6 ]=t1+t14;				a[j2+p6 ]=t2-t13;
		a[j1+p12]=t5+t18;				a[j2+p12]=t6-t17;
		a[j1+p4 ]=t9-t22;				a[j2+p4 ]=t10+t21;
		a[j1+p10]=t9+t22;				a[j2+p10]=t10-t21;
		a[j1+p2 ]=t5-t18;				a[j2+p2 ]=t6+t17;
		a[j1+p8 ]=t1-t14;				a[j2+p8 ]=t2+t13;
#endif
/*
printf("%20.15f  %20.15f\n", a[j1    ], a[j2    ]);
printf("%20.15f  %20.15f\n", a[j1+p2 ], a[j2+p2 ]);
printf("%20.15f  %20.15f\n", a[j1+p4 ], a[j2+p4 ]);
printf("%20.15f  %20.15f\n", a[j1+p6 ], a[j2+p6 ]);
printf("%20.15f  %20.15f\n", a[j1+p8 ], a[j2+p8 ]);
printf("%20.15f  %20.15f\n", a[j1+p10], a[j2+p10]);
printf("%20.15f  %20.15f\n", a[j1+p12], a[j2+p12]);
*/
/*...aj1p[1:13:2]r use t[3:27:4]; aj1p[1:13:2]i use t[4:28:4]	*/

		rt =t27;						it =t28;
		t27=t7 -rt;						t28=t8 -it;
		t7 =t7 +rt;						t8 =t8 +it;
		rt =t23;						it =t24;
		t23=t11 -rt;					t24=t12 -it;
		t11=t11 +rt;					t12=t12 +it;
		rt =t19;						it =t20;
		t19=t15 -rt;					t20=t16 -it;
		t15=t15 +rt;					t16=t16 +it;

	/* Version 1: direct computation of 3x3 matrix-vector product. */
#if LO_ADD
		a[j1+p7] = t3+t7+t11+t15;		a[j2+p7] = t4+t8+t12+t16;
		rt = t3+cc1*t7+cc2*t11+cc3*t15;	it = t4+cc1*t8+cc2*t12+cc3*t16;
		re = t3+cc2*t7+cc3*t11+cc1*t15;	im = t4+cc2*t8+cc3*t12+cc1*t16;
		t3 = t3+cc3*t7+cc1*t11+cc2*t15;	t4 = t4+cc3*t8+cc1*t12+cc2*t16;

		t7 = ss1*t27+ss2*t23+ss3*t19;	t8 = ss1*t28+ss2*t24+ss3*t20;
		t11= ss2*t27-ss3*t23-ss1*t19;	t12= ss2*t28-ss3*t24-ss1*t20;
		t15= ss3*t27-ss1*t23+ss2*t19;	t16= ss3*t28-ss1*t24+ss2*t20;
	/* Output permutation causes signs to get flipped here: */
		a[j1+p13]=rt+t8 ;				a[j2+p13]=it-t7 ;
		a[j1+p5 ]=re+t12;				a[j2+p5 ]=im-t11;
		a[j1+p11]=t3+t16;				a[j2+p11]=t4-t15;
		a[j1+p3 ]=t3-t16;				a[j2+p3 ]=t4+t15;
		a[j1+p9 ]=re-t12;				a[j2+p9 ]=im+t11;
		a[j1+p1 ]=rt-t8 ;				a[j2+p1 ]=it+t7 ;

	/* Version 2: */
#else
	/* Get cyclic convolution of (cc1,cc2,cc3) with (t7,t15,t11) and with (t8,t16,t12),
		 (note switch of last 2 t-elements in each case), add t1 and t2 to respective outputs: */

		rt = t7+t15+t11;a[j1+p7]=rt+t3;	it = t8+t16+t12; a[j2+p7]=it+t4;
		t3 = rt*c0+t3;					t4 = it*c0+t4;
		t7 = t7-t11;					t8 = t8 -t12;
		t11= t15-t11;					t12= t16-t12;
		t15=(t7+t11)*c3;				t16=(t8+t12)*c3;
		t7 = t7*c1;						t8 = t8 *c1;
		t11= t11*c2;					t12= t12*c2;
		rt = t7 -t15;					it = t8 -t16;
		t11= t11-t15;					t12= t12-t16;

		t7 = t3-rt-t11;					t8 = t4-it-t12;	/* C2	*/
		t11= t3+t11;					t12= t4+t12;	/* C3	*/
		t3 = t3+rt;						t4 = t4+it;	/* C1	*/

	/* Get cyclic convolution of (ss1,ss2,-ss3) with (t27,-t19,t23) and with (t28,-t20,t24),
		 (note permutation of t-elements and sign flips on ss3, t19 and t20): */

		t15=(t27-t19 +t23)*s0;			t16=(t28-t20+t24)*s0;
		t27= t27-t23;					t28= t28-t24;
		t23= t19+t23;					t24= t20+t24;		/* Flip sign here... */
		t19=(t23-t27)*s3;				t20=(t24-t28)*s3;	/* ...and here... */
		t27= t27*s1;					t28= t28*s1;
		t23= t23*s2;					t24= t24*s2;
		t27= t19+t27;					t28= t20+t28;		/* Sign flips get undone in these next 2. */
		t23= t19-t23;					t24= t20-t24;

		t19= t15-t27-t23;				t20= t16-t28-t24;	/* S2 */
		t23= t15+t23;					t24= t16+t24;		/*-S3 */
		t15= t15+t27;					t16= t16+t28;		/* S1 */

		a[j1+p13]=t3+t16;				a[j2+p13]=t4-t15;
		a[j1+p5 ]=t7+t20;				a[j2+p5 ]=t8-t19;
		a[j1+p11]=t11-t24;				a[j2+p11]=t12+t23;
		a[j1+p3 ]=t11+t24;				a[j2+p3 ]=t12-t23;
		a[j1+p9 ]=t7-t20;				a[j2+p9 ]=t8+t19;
		a[j1+p1 ]=t3-t16;				a[j2+p1 ]=t4+t15;
#endif
/*
printf("%20.15f  %20.15f\n", a[j1+p1 ], a[j2+p1 ]);
printf("%20.15f  %20.15f\n", a[j1+p3 ], a[j2+p3 ]);
printf("%20.15f  %20.15f\n", a[j1+p5 ], a[j2+p5 ]);
printf("%20.15f  %20.15f\n", a[j1+p7 ], a[j2+p7 ]);
printf("%20.15f  %20.15f\n", a[j1+p9 ], a[j2+p9 ]);
printf("%20.15f  %20.15f\n", a[j1+p11], a[j2+p11]);
printf("%20.15f  %20.15f\n", a[j1+p13], a[j2+p13]);
*/
	}
}

