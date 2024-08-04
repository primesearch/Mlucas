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

/**************/

int radix22_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-22 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-22 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n22, bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21
		,i,j,j1,j2,jstart,jhi,root_incr,k,khi,l,outer;
	static uint64 psave = 0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21, nsave = 0;
	const double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double radix_inv, n2inv;
	static double	cc1 =  0.84125353283118116886,	/* Real part of exp(i*2*pi/11), the radix-11 fundamental sincos datum	*/
			ss1 =  0.54064081745559758210,	/* Imag part of exp(i*2*pi/11).	*/
			cc2 =  0.41541501300188642553,	/* cos(2u)	*/
			ss2 =  0.90963199535451837140,	/* sin(2u)	*/
			cc3 = -0.14231483827328514043,	/* cos(3u)	*/
			ss3 =  0.98982144188093273237,	/* sin(3u)	*/
			cc4 = -0.65486073394528506404,	/* cos(4u)	*/
			ss4 =  0.75574957435425828378,	/* sin(4u)	*/
			cc5 = -0.95949297361449738988,	/* cos(5u)	*/
			ss5 =  0.28173255684142969773;	/* sin(5u)	*/
	double rt,it
		,cr1,cr2,cr3,cr4,cr5,ci1,ci2,ci3,ci4,ci5
		,sr1,sr2,sr3,sr4,sr5,si1,si2,si3,si4,si5
		,u01,u02,u03,u04,u05,u06,u07,u08,u09,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,u31,u32,u33,u34,u35,u36,u37,u38,u39,u40,u41,u42,u43,u44
		,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22
	,aj1p00r,aj1p01r,aj1p02r,aj1p03r,aj1p04r,aj1p05r,aj1p06r,aj1p07r,aj1p08r,aj1p09r,aj1p10r,aj1p11r,aj1p12r,aj1p13r,aj1p14r,aj1p15r,aj1p16r,aj1p17r,aj1p18r,aj1p19r,aj1p20r,aj1p21r
	,aj1p00i,aj1p01i,aj1p02i,aj1p03i,aj1p04i,aj1p05i,aj1p06i,aj1p07i,aj1p08i,aj1p09i,aj1p10i,aj1p11i,aj1p12i,aj1p13i,aj1p14i,aj1p15i,aj1p16i,aj1p17i,aj1p18i,aj1p19i,aj1p20i,aj1p21i
		,cy00,cy01,cy02,cy03,cy04,cy05,cy06,cy07,cy08,cy09,cy10,cy11,cy12,cy13,cy14,cy15,cy16,cy17,cy18,cy19,cy20,cy21,temp,frac,scale;
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
	double wt_re,wt_im, wi_re,wi_im;					/* Fermat/LOACC weights stuff */

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=-1;
	// Jan 2018: To support PRP-testing, read the LR-modpow-scalar-multiply-needed bit for the current iteration from the global array:
	double prp_mult = 1.0;
	if((TEST_TYPE & 0xfffffffe) == TEST_TYPE_PRP) {	// Mask off low bit to lump together PRP and PRP-C tests
		i = (iter-1) % ITERS_BETWEEN_CHECKPOINTS;	// Bit we need to read...iter-counter is unit-offset w.r.to iter-interval, hence the -1
		if((BASE_MULTIPLIER_BITS[i>>6] >> (i&63)) & 1)
			prp_mult = PRP_BASE;
	}

/*...change n22 and n_div_wt to non-static to work around a gcc compiler bug. */
	n22   = n/22;
	n_div_nwt = n22 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n22)
	{
		sprintf(cbuf,"ERROR: iter = %10d; NWT_BITS does not divide N/22 in radix22_ditN_cy_dif1.\n",iter);
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
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)22));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

/*   constant index offsets for array load/stores are here.	*/

		p01 = n22;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p10 = p09 + p01;
		p11 = p10 + p01;
		p12 = p11 + p01;
		p13 = p12 + p01;
		p14 = p13 + p01;
		p15 = p14 + p01;
		p16 = p15 + p01;
		p17 = p16 + p01;
		p18 = p17 + p01;
		p19 = p18 + p01;
		p20 = p19 + p01;
		p21 = p20 + p01;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p17 = p17 + ( (p17 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p19 = p19 + ( (p19 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p21 = p21 + ( (p21 >> DAT_BITS) << PAD_BITS );

		bjmodnini=0;
		for(j=0; j < n22; j++)
		{
			bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
		}
	}

/*...The radix-22 final DIT pass is here.	*/

	cy00= 0;		/* init carries	*/
	cy01= 0;
	cy02= 0;
	cy03= 0;
	cy04= 0;
	cy05= 0;
	cy06= 0;
	cy07= 0;
	cy08= 0;
	cy09= 0;
	cy10= 0;
	cy11= 0;
	cy12= 0;
	cy13= 0;
	cy14= 0;
	cy15= 0;
	cy16= 0;
	cy17= 0;
	cy18= 0;
	cy19= 0;
	cy20= 0;
	cy21= 0;

	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		cy00= -2;
	}

	*fracmax=0;	/* init max. fractional error	*/

	root_incr = 1;	/* init sincos array index increment (set = 1 for normal carry pass, = 0 for wrapper pass)	*/
	scale = n2inv;	// init inverse-weight scale factor = 2/n for normal carry pass, 1 for wrapper pass

	jstart = 0;
	jhi = jstart+nwt-1;
	khi = n_div_nwt;

for(outer=0; outer <= 1; outer++)
{
	i = 1;		/* Pointer to the BASE and BASEINV arrays. lowest-order digit is always a bigword (i = 1).	*/

	bjmodn00= 0;
	bjmodn01= bjmodnini;
	bjmodn02= bjmodn01+bjmodnini-n; bjmodn02= bjmodn02+ ( (-(int)((uint32)bjmodn02>> 31)) & n);
	bjmodn03= bjmodn02+bjmodnini-n; bjmodn03= bjmodn03+ ( (-(int)((uint32)bjmodn03>> 31)) & n);
	bjmodn04= bjmodn03+bjmodnini-n; bjmodn04= bjmodn04+ ( (-(int)((uint32)bjmodn04>> 31)) & n);
	bjmodn05= bjmodn04+bjmodnini-n; bjmodn05= bjmodn05+ ( (-(int)((uint32)bjmodn05>> 31)) & n);
	bjmodn06= bjmodn05+bjmodnini-n; bjmodn06= bjmodn06+ ( (-(int)((uint32)bjmodn06>> 31)) & n);
	bjmodn07= bjmodn06+bjmodnini-n; bjmodn07= bjmodn07+ ( (-(int)((uint32)bjmodn07>> 31)) & n);
	bjmodn08= bjmodn07+bjmodnini-n; bjmodn08= bjmodn08+ ( (-(int)((uint32)bjmodn08>> 31)) & n);
	bjmodn09= bjmodn08+bjmodnini-n; bjmodn09= bjmodn09+ ( (-(int)((uint32)bjmodn09>> 31)) & n);
	bjmodn10= bjmodn09+bjmodnini-n; bjmodn10= bjmodn10+ ( (-(int)((uint32)bjmodn10>> 31)) & n);
	bjmodn11= bjmodn10+bjmodnini-n; bjmodn11= bjmodn11+ ( (-(int)((uint32)bjmodn11>> 31)) & n);
	bjmodn12= bjmodn11+bjmodnini-n; bjmodn12= bjmodn12+ ( (-(int)((uint32)bjmodn12>> 31)) & n);
	bjmodn13= bjmodn12+bjmodnini-n; bjmodn13= bjmodn13+ ( (-(int)((uint32)bjmodn13>> 31)) & n);
	bjmodn14= bjmodn13+bjmodnini-n; bjmodn14= bjmodn14+ ( (-(int)((uint32)bjmodn14>> 31)) & n);
	bjmodn15= bjmodn14+bjmodnini-n; bjmodn15= bjmodn15+ ( (-(int)((uint32)bjmodn15>> 31)) & n);
	bjmodn16= bjmodn15+bjmodnini-n; bjmodn16= bjmodn16+ ( (-(int)((uint32)bjmodn16>> 31)) & n);
	bjmodn17= bjmodn16+bjmodnini-n; bjmodn17= bjmodn17+ ( (-(int)((uint32)bjmodn17>> 31)) & n);
	bjmodn18= bjmodn17+bjmodnini-n; bjmodn18= bjmodn18+ ( (-(int)((uint32)bjmodn18>> 31)) & n);
	bjmodn19= bjmodn18+bjmodnini-n; bjmodn19= bjmodn19+ ( (-(int)((uint32)bjmodn19>> 31)) & n);
	bjmodn20= bjmodn19+bjmodnini-n; bjmodn20= bjmodn20+ ( (-(int)((uint32)bjmodn20>> 31)) & n);
	bjmodn21= bjmodn20+bjmodnini-n; bjmodn21= bjmodn21+ ( (-(int)((uint32)bjmodn21>> 31)) & n);

	col=0;
	co2=(n >> nwt_bits)-1+22;
	co3=co2-22;		/* At the start of each new j-loop, co3=co2-radix(1)	*/

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
	/*
	!...gather the needed data (22 64-bit complex, i.e. 44 64-bit reals) and do 11 radix-2 transforms...
	*/
			u01=a[j1      ];		u02=a[j2    ];	/* x0a,b */
			rt =a[j1  +p01];		it =a[j2+p01];
			u03=u01-rt;				u04=u02-it;
			u01=u01+rt;				u02=u02+it;

			u05=a[j1  +p19];		u06=a[j2+p19];	/* x1a,b */
			rt =a[j1  +p18];		it =a[j2+p18];
			u07=u05-rt;				u08=u06-it;
			u05=u05+rt;				u06=u06+it;

			u09=a[j1  +p15];		u10=a[j2+p15];	/* x2a,b */
			rt =a[j1  +p14];		it =a[j2+p14];
			u11=u09-rt;				u12=u10-it;
			u09=u09+rt;				u10=u10+it;

			u13=a[j1  +p11];		u14=a[j2+p11];	/* x3a,b */
			rt =a[j1  +p10];		it =a[j2+p10];
			u15=u13-rt;				u16=u14-it;
			u13=u13+rt;				u14=u14+it;

			u17=a[j1  +p07];		u18=a[j2+p07];	/* x4a,b */
			rt =a[j1  +p06];		it =a[j2+p06];
			u19=u17-rt;				u20=u18-it;
			u17=u17+rt;				u18=u18+it;

			u21=a[j1  +p03];		u22=a[j2+p03];	/* x5a,b */
			rt =a[j1  +p02];		it =a[j2+p02];
			u23=u21-rt;				u24=u22-it;
			u21=u21+rt;				u22=u22+it;

			u25=a[j1  +p20];		u26=a[j2+p20];	/* x6a,b */
			rt =a[j1  +p21];		it =a[j2+p21];
			u27=u25-rt;				u28=u26-it;
			u25=u25+rt;				u26=u26+it;

			u29=a[j1  +p16];		u30=a[j2+p16];	/* x7a,b */
			rt =a[j1  +p17];		it =a[j2+p17];
			u31=u29-rt;				u32=u30-it;
			u29=u29+rt;				u30=u30+it;

			u33=a[j1  +p12];		u34=a[j2+p12];	/* x8a,b */
			rt =a[j1  +p13];		it =a[j2+p13];
			u35=u33-rt;				u36=u34-it;
			u33=u33+rt;				u34=u34+it;

			u37=a[j1  +p08];		u38=a[j2+p08];	/* x9a,b */
			rt =a[j1  +p09];		it =a[j2+p09];
			u39=u37-rt;				u40=u38-it;
			u37=u37+rt;				u38=u38+it;

			u41=a[j1  +p04];		u42=a[j2+p04];	/* x10a,b */
			rt =a[j1  +p05];		it =a[j2+p05];
			u43=u41-rt;				u44=u42-it;
			u41=u41+rt;				u42=u42+it;

	/*       ...and now do two radix-11 transforms.	*/

	/*...aj1p[0:20:2]r use u[1:41:4]; aj1p[0:20:2]i use u[2:42:4]	*/

			t01= u01;				t02= u02;		/* x0		*/

			t03= u05 + u41;			t04= u06 + u42;		/* x1 + x10	*/
			t21= u05 - u41;			t22= u06 - u42;		/* x1 - x10	*/

			t05= u09 + u37;			t06= u10 + u38;		/* x2 + x9	*/
			t19= u09 - u37;			t20= u10 - u38;		/* x2 - x9	*/

			t07= u13 + u33;			t08= u14 + u34;		/* x3 + x8	*/
			t17= u13 - u33;			t18= u14 - u34;		/* x3 - x8	*/

			t09= u17 + u29;			t10= u18 + u30;		/* x4 + x7	*/
			t15= u17 - u29;			t16= u18 - u30;		/* x4 - x7	*/

			t11= u21 + u25;			t12= u22 + u26;		/* x5 + x6	*/
			t13= u21 - u25;			t14= u22 - u26;		/* x5 - x6	*/

			aj1p00r = t01+t03+t05+t07+t09+t11;	aj1p00i = t02+t04+t06+t08+t10+t12;	/* X0	*/

			cr1= t01+cc1*t03+cc2*t05+cc3*t07+cc4*t09+cc5*t11;	ci1= t02+cc1*t04+cc2*t06+cc3*t08+cc4*t10+cc5*t12;	/* C1	*/
			cr2= t01+cc2*t03+cc4*t05+cc5*t07+cc3*t09+cc1*t11;	ci2= t02+cc2*t04+cc4*t06+cc5*t08+cc3*t10+cc1*t12;	/* C2	*/
			cr3= t01+cc3*t03+cc5*t05+cc2*t07+cc1*t09+cc4*t11;	ci3= t02+cc3*t04+cc5*t06+cc2*t08+cc1*t10+cc4*t12;	/* C3	*/
			cr4= t01+cc4*t03+cc3*t05+cc1*t07+cc5*t09+cc2*t11;	ci4= t02+cc4*t04+cc3*t06+cc1*t08+cc5*t10+cc2*t12;	/* C4	*/
			cr5= t01+cc5*t03+cc1*t05+cc4*t07+cc2*t09+cc3*t11;	ci5= t02+cc5*t04+cc1*t06+cc4*t08+cc2*t10+cc3*t12;	/* C5	*/

			sr1=     ss1*t21+ss2*t19+ss3*t17+ss4*t15+ss5*t13;	si1=     ss1*t22+ss2*t20+ss3*t18+ss4*t16+ss5*t14;	/* S1	*/
			sr2=     ss2*t21+ss4*t19-ss5*t17-ss3*t15-ss1*t13;	si2=     ss2*t22+ss4*t20-ss5*t18-ss3*t16-ss1*t14;	/* S2	*/
			sr3=     ss3*t21-ss5*t19-ss2*t17+ss1*t15+ss4*t13;	si3=     ss3*t22-ss5*t20-ss2*t18+ss1*t16+ss4*t14;	/* S3	*/
			sr4=     ss4*t21-ss3*t19+ss1*t17+ss5*t15-ss2*t13;	si4=     ss4*t22-ss3*t20+ss1*t18+ss5*t16-ss2*t14;	/* S4	*/
			sr5=     ss5*t21-ss1*t19+ss4*t17-ss2*t15+ss3*t13;	si5=     ss5*t22-ss1*t20+ss4*t18-ss2*t16+ss3*t14;	/* S5	*/

	/*...Inline multiply of sine parts by +-I into finishing phase...	*/

			aj1p10r=cr1+si1;		aj1p10i=ci1-sr1;	/* X1 = C1 + I*S1	*/
			aj1p20r=cr2+si2;		aj1p20i=ci2-sr2;	/* X2 = C2 + I*S2	*/
			aj1p08r=cr3+si3;		aj1p08i=ci3-sr3;	/* X3 = C3 + I*S3	*/
			aj1p18r=cr4+si4;		aj1p18i=ci4-sr4;	/* X4 = C4 + I*S4	*/
			aj1p06r=cr5+si5;		aj1p06i=ci5-sr5;	/* X5 = C5 + I*S5	*/
			aj1p16r=cr5-si5;		aj1p16i=ci5+sr5;	/* X6 =	C5 - I*S5	*/
			aj1p04r=cr4-si4;		aj1p04i=ci4+sr4;	/* X7 =	C4 - I*S4	*/
			aj1p14r=cr3-si3;		aj1p14i=ci3+sr3;	/* X8 =	C3 - I*S3	*/
			aj1p02r=cr2-si2;		aj1p02i=ci2+sr2;	/* X9 =	C2 - I*S2	*/
			aj1p12r=cr1-si1;		aj1p12i=ci1+sr1;	/* X10=	C1 - I*S1	*/

	/*...aj1p[1:21:2]r use u[3:43:4]; aj1p[1:21:2]i use u[4:44:4]	*/

			t01= u03;				t02= u04;		/* x0		*/

			t03= u07 + u43;			t04= u08 + u44;		/* x1 + x10	*/
			t21= u07 - u43;			t22= u08 - u44;		/* x1 - x10	*/

			t05= u11 + u39;			t06= u12 + u40;		/* x2 + x9	*/
			t19= u11 - u39;			t20= u12 - u40;		/* x2 - x9	*/

			t07= u15 + u35;			t08= u16 + u36;		/* x3 + x8	*/
			t17= u15 - u35;			t18= u16 - u36;		/* x3 - x8	*/

			t09= u19 + u31;			t10= u20 + u32;		/* x4 + x7	*/
			t15= u19 - u31;			t16= u20 - u32;		/* x4 - x7	*/

			t11= u23 + u27;			t12= u24 + u28;		/* x5 + x6	*/
			t13= u23 - u27;			t14= u24 - u28;		/* x5 - x6	*/

			aj1p11r = t01+t03+t05+t07+t09+t11;	aj1p11i = t02+t04+t06+t08+t10+t12;	/* X0	*/

			cr1= t01+cc1*t03+cc2*t05+cc3*t07+cc4*t09+cc5*t11;	ci1= t02+cc1*t04+cc2*t06+cc3*t08+cc4*t10+cc5*t12;	/* C1	*/
			cr2= t01+cc2*t03+cc4*t05+cc5*t07+cc3*t09+cc1*t11;	ci2= t02+cc2*t04+cc4*t06+cc5*t08+cc3*t10+cc1*t12;	/* C2	*/
			cr3= t01+cc3*t03+cc5*t05+cc2*t07+cc1*t09+cc4*t11;	ci3= t02+cc3*t04+cc5*t06+cc2*t08+cc1*t10+cc4*t12;	/* C3	*/
			cr4= t01+cc4*t03+cc3*t05+cc1*t07+cc5*t09+cc2*t11;	ci4= t02+cc4*t04+cc3*t06+cc1*t08+cc5*t10+cc2*t12;	/* C4	*/
			cr5= t01+cc5*t03+cc1*t05+cc4*t07+cc2*t09+cc3*t11;	ci5= t02+cc5*t04+cc1*t06+cc4*t08+cc2*t10+cc3*t12;	/* C5	*/

			sr1=     ss1*t21+ss2*t19+ss3*t17+ss4*t15+ss5*t13;	si1=     ss1*t22+ss2*t20+ss3*t18+ss4*t16+ss5*t14;	/* S1	*/
			sr2=     ss2*t21+ss4*t19-ss5*t17-ss3*t15-ss1*t13;	si2=     ss2*t22+ss4*t20-ss5*t18-ss3*t16-ss1*t14;	/* S2	*/
			sr3=     ss3*t21-ss5*t19-ss2*t17+ss1*t15+ss4*t13;	si3=     ss3*t22-ss5*t20-ss2*t18+ss1*t16+ss4*t14;	/* S3	*/
			sr4=     ss4*t21-ss3*t19+ss1*t17+ss5*t15-ss2*t13;	si4=     ss4*t22-ss3*t20+ss1*t18+ss5*t16-ss2*t14;	/* S4	*/
			sr5=     ss5*t21-ss1*t19+ss4*t17-ss2*t15+ss3*t13;	si5=     ss5*t22-ss1*t20+ss4*t18-ss2*t16+ss3*t14;	/* S5	*/

	/*...Inline multiply of sine parts by +-I into finishing phase...	*/

			aj1p21r=cr1+si1;		aj1p21i=ci1-sr1;	/* X1 = C1 + I*S1	*/
			aj1p09r=cr2+si2;		aj1p09i=ci2-sr2;	/* X2 = C2 + I*S2	*/
			aj1p19r=cr3+si3;		aj1p19i=ci3-sr3;	/* X3 = C3 + I*S3	*/
			aj1p07r=cr4+si4;		aj1p07i=ci4-sr4;	/* X4 = C4 + I*S4	*/
			aj1p17r=cr5+si5;		aj1p17i=ci5-sr5;	/* X5 = C5 + I*S5	*/
			aj1p05r=cr5-si5;		aj1p05i=ci5+sr5;	/* X6 =	C5 - I*S5	*/
			aj1p15r=cr4-si4;		aj1p15i=ci4+sr4;	/* X7 =	C4 - I*S4	*/
			aj1p03r=cr3-si3;		aj1p03i=ci3+sr3;	/* X8 =	C3 - I*S3	*/
			aj1p13r=cr2-si2;		aj1p13i=ci2+sr2;	/* X9 =	C2 - I*S2	*/
			aj1p01r=cr1-si1;		aj1p01i=ci1+sr1;	/* X10=	C1 - I*S1	*/

/*...Now do the carries. Since the outputs would
    normally be getting dispatched to 22 separate blocks of the A-array, we need 22 separate carries.	*/

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

		   cmplx_carry_norm_errcheck0(aj1p00r,aj1p00i,cy00,bjmodn00,0 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p01r,aj1p01i,cy01,bjmodn01,1 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p02r,aj1p02i,cy02,bjmodn02,2 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p03r,aj1p03i,cy03,bjmodn03,3 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p04r,aj1p04i,cy04,bjmodn04,4 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p05r,aj1p05i,cy05,bjmodn05,5 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p06r,aj1p06i,cy06,bjmodn06,6 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p07r,aj1p07i,cy07,bjmodn07,7 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p08r,aj1p08i,cy08,bjmodn08,8 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p09r,aj1p09i,cy09,bjmodn09,9 ,prp_mult);
			cmplx_carry_norm_errcheck(aj1p10r,aj1p10i,cy10,bjmodn10,10,prp_mult);
			cmplx_carry_norm_errcheck(aj1p11r,aj1p11i,cy11,bjmodn11,11,prp_mult);
			cmplx_carry_norm_errcheck(aj1p12r,aj1p12i,cy12,bjmodn12,12,prp_mult);
			cmplx_carry_norm_errcheck(aj1p13r,aj1p13i,cy13,bjmodn13,13,prp_mult);
			cmplx_carry_norm_errcheck(aj1p14r,aj1p14i,cy14,bjmodn14,14,prp_mult);
			cmplx_carry_norm_errcheck(aj1p15r,aj1p15i,cy15,bjmodn15,15,prp_mult);
			cmplx_carry_norm_errcheck(aj1p16r,aj1p16i,cy16,bjmodn16,16,prp_mult);
			cmplx_carry_norm_errcheck(aj1p17r,aj1p17i,cy17,bjmodn17,17,prp_mult);
			cmplx_carry_norm_errcheck(aj1p18r,aj1p18i,cy18,bjmodn18,18,prp_mult);
			cmplx_carry_norm_errcheck(aj1p19r,aj1p19i,cy19,bjmodn19,19,prp_mult);
			cmplx_carry_norm_errcheck(aj1p20r,aj1p20i,cy20,bjmodn20,20,prp_mult);
			cmplx_carry_norm_errcheck(aj1p21r,aj1p21i,cy21,bjmodn21,21,prp_mult);

			i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
				 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/

/*...The radix-22 DIF pass is here:	*/
	#if PFETCH
	add0 = &a[j1];
	prefetch_p_doubles(add0);
	#endif

	/*...First radix-11 block uses aj1p[0:20:2] as inputs:	*/

			t01 = aj1p00r;					t02 = aj1p00i;		/* x0		*/
			t03 = aj1p20r+aj1p02r;				t04 = aj1p20i+aj1p02i;	/* x1 + x10	*/
			t05 = aj1p18r+aj1p04r;				t06 = aj1p18i+aj1p04i;	/* x2 + x9	*/
			t07 = aj1p16r+aj1p06r;				t08 = aj1p16i+aj1p06i;	/* x3 + x8	*/
			t09 = aj1p14r+aj1p08r;				t10 = aj1p14i+aj1p08i;	/* x4 + x7	*/
			t11 = aj1p12r+aj1p10r;				t12 = aj1p12i+aj1p10i;	/* x5 + x6	*/
	#if PFETCH
	addr = add0+p01;
	prefetch_p_doubles(addr);
	#endif
			t13 = aj1p12r-aj1p10r;				t14 = aj1p12i-aj1p10i;	/* x5 - x6	*/
			t15 = aj1p14r-aj1p08r;				t16 = aj1p14i-aj1p08i;	/* x4 - x7	*/
			t17 = aj1p16r-aj1p06r;				t18 = aj1p16i-aj1p06i;	/* x3 - x8	*/
			t19 = aj1p18r-aj1p04r;				t20 = aj1p18i-aj1p04i;	/* x2 - x9	*/
			t21 = aj1p20r-aj1p02r;				t22 = aj1p20i-aj1p02i;	/* x1 - x10	*/

	#if PFETCH
	addr = add0+p02;
	prefetch_p_doubles(addr);
	#endif
			aj1p00r = t01+t03+t05+t07+t09+t11;			aj1p00i = t02+t04+t06+t08+t10+t12;	/* X0	*/

			cr1= t01+cc1*t03+cc2*t05+cc3*t07+cc4*t09+cc5*t11;	ci1= t02+cc1*t04+cc2*t06+cc3*t08+cc4*t10+cc5*t12;	/* C1	*/
			cr2= t01+cc2*t03+cc4*t05+cc5*t07+cc3*t09+cc1*t11;	ci2= t02+cc2*t04+cc4*t06+cc5*t08+cc3*t10+cc1*t12;	/* C2	*/
	#if PFETCH
	addr = add0+p03;
	prefetch_p_doubles(addr);
	#endif
			cr3= t01+cc3*t03+cc5*t05+cc2*t07+cc1*t09+cc4*t11;	ci3= t02+cc3*t04+cc5*t06+cc2*t08+cc1*t10+cc4*t12;	/* C3	*/
			cr4= t01+cc4*t03+cc3*t05+cc1*t07+cc5*t09+cc2*t11;	ci4= t02+cc4*t04+cc3*t06+cc1*t08+cc5*t10+cc2*t12;	/* C4	*/
			cr5= t01+cc5*t03+cc1*t05+cc4*t07+cc2*t09+cc3*t11;	ci5= t02+cc5*t04+cc1*t06+cc4*t08+cc2*t10+cc3*t12;	/* C5	*/
	#if PFETCH
	addr = add0+p04;
	prefetch_p_doubles(addr);
	#endif

			sr1=     ss1*t21+ss2*t19+ss3*t17+ss4*t15+ss5*t13;	si1=     ss1*t22+ss2*t20+ss3*t18+ss4*t16+ss5*t14;	/* S1	*/
			sr2=     ss2*t21+ss4*t19-ss5*t17-ss3*t15-ss1*t13;	si2=     ss2*t22+ss4*t20-ss5*t18-ss3*t16-ss1*t14;	/* S2	*/
			sr3=     ss3*t21-ss5*t19-ss2*t17+ss1*t15+ss4*t13;	si3=     ss3*t22-ss5*t20-ss2*t18+ss1*t16+ss4*t14;	/* S3	*/
	#if PFETCH
	addr = add0+p05;
	prefetch_p_doubles(addr);
	#endif
			sr4=     ss4*t21-ss3*t19+ss1*t17+ss5*t15-ss2*t13;	si4=     ss4*t22-ss3*t20+ss1*t18+ss5*t16-ss2*t14;	/* S4	*/
			sr5=     ss5*t21-ss1*t19+ss4*t17-ss2*t15+ss3*t13;	si5=     ss5*t22-ss1*t20+ss4*t18-ss2*t16+ss3*t14;	/* S5	*/

	/*...Inline multiply of sine parts by +-I into finishing phase...	*/

	#if PFETCH
	addr = add0+p06;
	prefetch_p_doubles(addr);
	#endif
			aj1p02r = cr1-si1;		aj1p02i = ci1+sr1;	/* X1 = C1 + I*S1	*/
			aj1p04r = cr2-si2;		aj1p04i = ci2+sr2;	/* X2 = C2 + I*S2	*/
			aj1p06r = cr3-si3;		aj1p06i = ci3+sr3;	/* X3 = C3 + I*S3	*/
			aj1p08r = cr4-si4;		aj1p08i = ci4+sr4;	/* X4 = C4 + I*S4	*/
			aj1p10r = cr5-si5;		aj1p10i = ci5+sr5;	/* X5 = C5 + I*S5	*/
	#if PFETCH
	addr = add0+p07;
	prefetch_p_doubles(addr);
	#endif
			aj1p12r = cr5+si5;		aj1p12i = ci5-sr5;	/* X6 =	C5 - I*S5	*/
			aj1p14r = cr4+si4;		aj1p14i = ci4-sr4;	/* X7 =	C4 - I*S4	*/
			aj1p16r = cr3+si3;		aj1p16i = ci3-sr3;	/* X8 =	C3 - I*S3	*/
			aj1p18r = cr2+si2;		aj1p18i = ci2-sr2;	/* X9 =	C2 - I*S2	*/
			aj1p20r = cr1+si1;		aj1p20i = ci1-sr1;	/* X10=	C1 - I*S1	*/

	/*...Second radix-11 block uses aj1p[1:21:2] as inputs:	*/

	#if PFETCH
	addr = add0+p08;
	prefetch_p_doubles(addr);
	#endif
			t01 = aj1p11r;					t02 = aj1p11i;		/* x0		*/
			t03 = aj1p09r+aj1p13r;				t04 = aj1p09i+aj1p13i;	/* x1 + x10	*/
			t05 = aj1p07r+aj1p15r;				t06 = aj1p07i+aj1p15i;	/* x2 + x9	*/
			t07 = aj1p05r+aj1p17r;				t08 = aj1p05i+aj1p17i;	/* x3 + x8	*/
			t09 = aj1p03r+aj1p19r;				t10 = aj1p03i+aj1p19i;	/* x4 + x7	*/
			t11 = aj1p01r+aj1p21r;				t12 = aj1p01i+aj1p21i;	/* x5 + x6	*/
	#if PFETCH
	addr = add0+p09;
	prefetch_p_doubles(addr);
	#endif
			t13 = aj1p01r-aj1p21r;				t14 = aj1p01i-aj1p21i;	/* x5 - x6	*/
			t15 = aj1p03r-aj1p19r;				t16 = aj1p03i-aj1p19i;	/* x4 - x7	*/
			t17 = aj1p05r-aj1p17r;				t18 = aj1p05i-aj1p17i;	/* x3 - x8	*/
			t19 = aj1p07r-aj1p15r;				t20 = aj1p07i-aj1p15i;	/* x2 - x9	*/
			t21 = aj1p09r-aj1p13r;				t22 = aj1p09i-aj1p13i;	/* x1 - x10	*/

	#if PFETCH
	addr = add0+p10;
	prefetch_p_doubles(addr);
	#endif
			aj1p01r = t01+t03+t05+t07+t09+t11;			aj1p01i = t02+t04+t06+t08+t10+t12;	/* X0	*/

			cr1= t01+cc1*t03+cc2*t05+cc3*t07+cc4*t09+cc5*t11;	ci1= t02+cc1*t04+cc2*t06+cc3*t08+cc4*t10+cc5*t12;	/* C1	*/
			cr2= t01+cc2*t03+cc4*t05+cc5*t07+cc3*t09+cc1*t11;	ci2= t02+cc2*t04+cc4*t06+cc5*t08+cc3*t10+cc1*t12;	/* C2	*/
	#if PFETCH
	addr = add0+p11;
	prefetch_p_doubles(addr);
	#endif
			cr3= t01+cc3*t03+cc5*t05+cc2*t07+cc1*t09+cc4*t11;	ci3= t02+cc3*t04+cc5*t06+cc2*t08+cc1*t10+cc4*t12;	/* C3	*/
			cr4= t01+cc4*t03+cc3*t05+cc1*t07+cc5*t09+cc2*t11;	ci4= t02+cc4*t04+cc3*t06+cc1*t08+cc5*t10+cc2*t12;	/* C4	*/
			cr5= t01+cc5*t03+cc1*t05+cc4*t07+cc2*t09+cc3*t11;	ci5= t02+cc5*t04+cc1*t06+cc4*t08+cc2*t10+cc3*t12;	/* C5	*/

	#if PFETCH
	addr = add0+p12;
	prefetch_p_doubles(addr);
	#endif
			sr1=     ss1*t21+ss2*t19+ss3*t17+ss4*t15+ss5*t13;	si1=     ss1*t22+ss2*t20+ss3*t18+ss4*t16+ss5*t14;	/* S1	*/
			sr2=     ss2*t21+ss4*t19-ss5*t17-ss3*t15-ss1*t13;	si2=     ss2*t22+ss4*t20-ss5*t18-ss3*t16-ss1*t14;	/* S2	*/
			sr3=     ss3*t21-ss5*t19-ss2*t17+ss1*t15+ss4*t13;	si3=     ss3*t22-ss5*t20-ss2*t18+ss1*t16+ss4*t14;	/* S3	*/
	#if PFETCH
	addr = add0+p13;
	prefetch_p_doubles(addr);
	#endif
			sr4=     ss4*t21-ss3*t19+ss1*t17+ss5*t15-ss2*t13;	si4=     ss4*t22-ss3*t20+ss1*t18+ss5*t16-ss2*t14;	/* S4	*/
			sr5=     ss5*t21-ss1*t19+ss4*t17-ss2*t15+ss3*t13;	si5=     ss5*t22-ss1*t20+ss4*t18-ss2*t16+ss3*t14;	/* S5	*/

	/*...Inline multiply of sine parts by +-I into finishing phase...	*/

	#if PFETCH
	addr = add0+p14;
	prefetch_p_doubles(addr);
	#endif
			aj1p03r = cr1-si1;		aj1p03i = ci1+sr1;	/* X1 = C1 + I*S1	*/
			aj1p05r = cr2-si2;		aj1p05i = ci2+sr2;	/* X2 = C2 + I*S2	*/
			aj1p07r = cr3-si3;		aj1p07i = ci3+sr3;	/* X3 = C3 + I*S3	*/
			aj1p09r = cr4-si4;		aj1p09i = ci4+sr4;	/* X4 = C4 + I*S4	*/
			aj1p11r = cr5-si5;		aj1p11i = ci5+sr5;	/* X5 = C5 + I*S5	*/
	#if PFETCH
	addr = add0+p15;
	prefetch_p_doubles(addr);
	#endif
			aj1p13r = cr5+si5;		aj1p13i = ci5-sr5;	/* X6 =	C5 - I*S5	*/
			aj1p15r = cr4+si4;		aj1p15i = ci4-sr4;	/* X7 =	C4 - I*S4	*/
			aj1p17r = cr3+si3;		aj1p17i = ci3-sr3;	/* X8 =	C3 - I*S3	*/
			aj1p19r = cr2+si2;		aj1p19i = ci2-sr2;	/* X9 =	C2 - I*S2	*/
			aj1p21r = cr1+si1;		aj1p21i = ci1-sr1;	/* X10=	C1 - I*S1	*/

	/*...and now do 11 radix-2 transforms:	*/

	#if PFETCH
	addr = add0+p16;
	prefetch_p_doubles(addr);
	#endif
			a[j1    ]=aj1p00r+aj1p01r;	a[j2    ]=aj1p00i+aj1p01i;
			a[j1+p01]=aj1p00r-aj1p01r;	a[j2+p01]=aj1p00i-aj1p01i;

			a[j1+p20]=aj1p02r+aj1p03r;	a[j2+p20]=aj1p02i+aj1p03i;
			a[j1+p21]=aj1p02r-aj1p03r;	a[j2+p21]=aj1p02i-aj1p03i;

	#if PFETCH
	addr = add0+p17;
	prefetch_p_doubles(addr);
	#endif
			a[j1+p19]=aj1p04r+aj1p05r;	a[j2+p19]=aj1p04i+aj1p05i;
			a[j1+p18]=aj1p04r-aj1p05r;	a[j2+p18]=aj1p04i-aj1p05i;

			a[j1+p16]=aj1p06r+aj1p07r;	a[j2+p16]=aj1p06i+aj1p07i;
			a[j1+p17]=aj1p06r-aj1p07r;	a[j2+p17]=aj1p06i-aj1p07i;

	#if PFETCH
	addr = add0+p18;
	prefetch_p_doubles(addr);
	#endif
			a[j1+p15]=aj1p08r+aj1p09r;	a[j2+p15]=aj1p08i+aj1p09i;
			a[j1+p14]=aj1p08r-aj1p09r;	a[j2+p14]=aj1p08i-aj1p09i;

			a[j1+p12]=aj1p10r+aj1p11r;	a[j2+p12]=aj1p10i+aj1p11i;
			a[j1+p13]=aj1p10r-aj1p11r;	a[j2+p13]=aj1p10i-aj1p11i;

	#if PFETCH
	addr = add0+p19;
	prefetch_p_doubles(addr);
	#endif
			a[j1+p11]=aj1p12r+aj1p13r;	a[j2+p11]=aj1p12i+aj1p13i;
			a[j1+p10]=aj1p12r-aj1p13r;	a[j2+p10]=aj1p12i-aj1p13i;

			a[j1+p08]=aj1p14r+aj1p15r;	a[j2+p08]=aj1p14i+aj1p15i;
			a[j1+p09]=aj1p14r-aj1p15r;	a[j2+p09]=aj1p14i-aj1p15i;

	#if PFETCH
	addr = add0+p20;
	prefetch_p_doubles(addr);
	#endif
			a[j1+p07]=aj1p16r+aj1p17r;	a[j2+p07]=aj1p16i+aj1p17i;
			a[j1+p06]=aj1p16r-aj1p17r;	a[j2+p06]=aj1p16i-aj1p17i;

			a[j1+p04]=aj1p18r+aj1p19r;	a[j2+p04]=aj1p18i+aj1p19i;
			a[j1+p05]=aj1p18r-aj1p19r;	a[j2+p05]=aj1p18i-aj1p19i;

	#if PFETCH
	addr = add0+p21;
	prefetch_p_doubles(addr);
	#endif
			a[j1+p03]=aj1p20r+aj1p21r;	a[j2+p03]=aj1p20i+aj1p21i;
			a[j1+p02]=aj1p20r-aj1p21r;	a[j2+p02]=aj1p20i-aj1p21i;
		}

		jstart += nwt;
		jhi    += nwt;
		col += 22;
		co3 -= 22;

	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-22 forward DIF FFT of the first block of 22 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 22 outputs of (1);
!   (3) Reweight and perform a radix-22 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 22 elements and repeat (1-4).
*/
	t01 = cy21;
	cy21= cy20;
	cy20= cy19;
	cy19= cy18;
	cy18= cy17;
	cy17= cy16;
	cy16= cy15;
	cy15= cy14;
	cy14= cy13;
	cy13= cy12;
	cy12= cy11;
	cy11= cy10;
	cy10= cy09;
	cy09= cy08;
	cy08= cy07;
	cy07= cy06;
	cy06= cy05;
	cy05= cy04;
	cy04= cy03;
	cy03= cy02;
	cy02= cy01;
	cy01= cy00;
	cy00= t01 ;

	root_incr = 0;
	scale = prp_mult = 1;

	jstart = 0;
	jhi = 7;
	khi = 1;

	for(j=0; j<=jhi; j++)
	{
		a[j    ] *= radix_inv;
		a[j+p01] *= radix_inv;
		a[j+p02] *= radix_inv;
		a[j+p03] *= radix_inv;
		a[j+p04] *= radix_inv;
		a[j+p05] *= radix_inv;
		a[j+p06] *= radix_inv;
		a[j+p07] *= radix_inv;
		a[j+p08] *= radix_inv;
		a[j+p09] *= radix_inv;
		a[j+p10] *= radix_inv;
		a[j+p11] *= radix_inv;
		a[j+p12] *= radix_inv;
		a[j+p13] *= radix_inv;
		a[j+p14] *= radix_inv;
		a[j+p15] *= radix_inv;
		a[j+p16] *= radix_inv;
		a[j+p17] *= radix_inv;
		a[j+p18] *= radix_inv;
		a[j+p19] *= radix_inv;
		a[j+p20] *= radix_inv;
		a[j+p21] *= radix_inv;
	}
}

	if(fabs(cy00)+fabs(cy01)+fabs(cy02)+fabs(cy03)+fabs(cy04)+fabs(cy05)+fabs(cy06)+fabs(cy07)+fabs(cy08)+fabs(cy09)+fabs(cy10)+fabs(cy11)+fabs(cy12)+fabs(cy13)+fabs(cy14)+fabs(cy15)+fabs(cy16)+fabs(cy17)+fabs(cy18)+fabs(cy19)+fabs(cy20)+fabs(cy21) != 0.0)
	{
		sprintf(cbuf,"ERROR: iter = %10d; nonzero exit carry in radix22_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
		mlucas_fprint(cbuf,INTERACT);
		err=ERR_CARRY;
		return(err);
	}
	*fracmax = maxerr;
	return(0);
}

/***************/

void radix22_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-22 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix11_dif_pass for details on the radix-11 subtransforms.
*/
	int j,j1,j2;
	static int n22,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21, first_entry=TRUE;
	static double	cc1 =  0.84125353283118116886,	/* Real part of exp(i*2*pi/11), the radix-11 fundamental sincos datum	*/
			ss1 =  0.54064081745559758210,	/* Imag part of exp(i*2*pi/11).	*/
			cc2 =  0.41541501300188642553,	/* cos(2u)	*/
			ss2 =  0.90963199535451837140,	/* sin(2u)	*/
			cc3 = -0.14231483827328514043,	/* cos(3u)	*/
			ss3 =  0.98982144188093273237,	/* sin(3u)	*/
			cc4 = -0.65486073394528506404,	/* cos(4u)	*/
			ss4 =  0.75574957435425828378,	/* sin(4u)	*/
			cc5 = -0.95949297361449738988,	/* cos(5u)	*/
			ss5 =  0.28173255684142969773;	/* sin(5u)	*/
	double   cr1,cr2,cr3,cr4,cr5,ci1,ci2,ci3,ci4,ci5
		,sr1,sr2,sr3,sr4,sr5,si1,si2,si3,si4,si5
		,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22
	,aj1p00r,aj1p01r,aj1p02r,aj1p03r,aj1p04r,aj1p05r,aj1p06r,aj1p07r,aj1p08r,aj1p09r,aj1p10r,aj1p11r,aj1p12r,aj1p13r,aj1p14r,aj1p15r,aj1p16r,aj1p17r,aj1p18r,aj1p19r,aj1p20r,aj1p21r
	,aj1p00i,aj1p01i,aj1p02i,aj1p03i,aj1p04i,aj1p05i,aj1p06i,aj1p07i,aj1p08i,aj1p09i,aj1p10i,aj1p11i,aj1p12i,aj1p13i,aj1p14i,aj1p15i,aj1p16i,aj1p17i,aj1p18i,aj1p19i,aj1p20i,aj1p21i;

	if(!first_entry && (n/22) != n22)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n22=n/22;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n22;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p10 = p09 + p01;
		p11 = p10 + p01;
		p12 = p11 + p01;
		p13 = p12 + p01;
		p14 = p13 + p01;
		p15 = p14 + p01;
		p16 = p15 + p01;
		p17 = p16 + p01;
		p18 = p17 + p01;
		p19 = p18 + p01;
		p20 = p19 + p01;
		p21 = p20 + p01;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p17 = p17 + ( (p17 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p19 = p19 + ( (p19 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p21 = p21 + ( (p21 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-22 pass is here.	*/

	for(j=0; j < n22; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*...gather the needed data (22 64-bit complex, i.e. 44 64-bit reals)...*/

	/*
		Twiddleless version requires us to swap inputs as follows:
		indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21
			  => 0,11,20, 9,18, 7,16, 5,14, 3,12, 1,10,21, 8,19, 6,17, 4,15, 2,13 modulo 22.
		I.e. start out with first 11-vector of indices {0,2,4,6,8,10,12,14,16,18,20},
		permute those according to {0, 2, 4, 6, 8,10,12,14,16,18,20}*21%22
								=> {0,20,18,16,14,12,10, 8, 6, 4, 2},
		then each is head of a length-2 list of indices with decrement 11.
	*/

	/*...First radix-11 block uses aj1p[0:20:2] as inputs:	*/

		t01 = a[j1    ];							t02 = a[j2    ];		/* x0		*/
		t03 = a[j1+p20]+a[j1+p02];				t04 = a[j2+p20]+a[j2+p02];	/* x1 + x10	*/
		t05 = a[j1+p18]+a[j1+p04];				t06 = a[j2+p18]+a[j2+p04];	/* x2 + x9	*/
		t07 = a[j1+p16]+a[j1+p06];				t08 = a[j2+p16]+a[j2+p06];	/* x3 + x8	*/
		t09 = a[j1+p14]+a[j1+p08];				t10 = a[j2+p14]+a[j2+p08];	/* x4 + x7	*/
		t11 = a[j1+p12]+a[j1+p10];				t12 = a[j2+p12]+a[j2+p10];	/* x5 + x6	*/
		t13 = a[j1+p12]-a[j1+p10];				t14 = a[j2+p12]-a[j2+p10];	/* x5 - x6	*/
		t15 = a[j1+p14]-a[j1+p08];				t16 = a[j2+p14]-a[j2+p08];	/* x4 - x7	*/
		t17 = a[j1+p16]-a[j1+p06];				t18 = a[j2+p16]-a[j2+p06];	/* x3 - x8	*/
		t19 = a[j1+p18]-a[j1+p04];				t20 = a[j2+p18]-a[j2+p04];	/* x2 - x9	*/
		t21 = a[j1+p20]-a[j1+p02];				t22 = a[j2+p20]-a[j2+p02];	/* x1 - x10	*/

		aj1p00r = t01+t03+t05+t07+t09+t11;			aj1p00i = t02+t04+t06+t08+t10+t12;	/* X0	*/

		cr1= t01+cc1*t03+cc2*t05+cc3*t07+cc4*t09+cc5*t11;	ci1= t02+cc1*t04+cc2*t06+cc3*t08+cc4*t10+cc5*t12;	/* C1	*/
		cr2= t01+cc2*t03+cc4*t05+cc5*t07+cc3*t09+cc1*t11;	ci2= t02+cc2*t04+cc4*t06+cc5*t08+cc3*t10+cc1*t12;	/* C2	*/
		cr3= t01+cc3*t03+cc5*t05+cc2*t07+cc1*t09+cc4*t11;	ci3= t02+cc3*t04+cc5*t06+cc2*t08+cc1*t10+cc4*t12;	/* C3	*/
		cr4= t01+cc4*t03+cc3*t05+cc1*t07+cc5*t09+cc2*t11;	ci4= t02+cc4*t04+cc3*t06+cc1*t08+cc5*t10+cc2*t12;	/* C4	*/
		cr5= t01+cc5*t03+cc1*t05+cc4*t07+cc2*t09+cc3*t11;	ci5= t02+cc5*t04+cc1*t06+cc4*t08+cc2*t10+cc3*t12;	/* C5	*/

		sr1=     ss1*t21+ss2*t19+ss3*t17+ss4*t15+ss5*t13;	si1=     ss1*t22+ss2*t20+ss3*t18+ss4*t16+ss5*t14;	/* S1	*/
		sr2=     ss2*t21+ss4*t19-ss5*t17-ss3*t15-ss1*t13;	si2=     ss2*t22+ss4*t20-ss5*t18-ss3*t16-ss1*t14;	/* S2	*/
		sr3=     ss3*t21-ss5*t19-ss2*t17+ss1*t15+ss4*t13;	si3=     ss3*t22-ss5*t20-ss2*t18+ss1*t16+ss4*t14;	/* S3	*/
		sr4=     ss4*t21-ss3*t19+ss1*t17+ss5*t15-ss2*t13;	si4=     ss4*t22-ss3*t20+ss1*t18+ss5*t16-ss2*t14;	/* S4	*/
		sr5=     ss5*t21-ss1*t19+ss4*t17-ss2*t15+ss3*t13;	si5=     ss5*t22-ss1*t20+ss4*t18-ss2*t16+ss3*t14;	/* S5	*/

/*...Inline multiply of sine parts by +-I into finishing phase...	*/

		aj1p02r = cr1-si1;		aj1p02i = ci1+sr1;	/* X1 = C1 + I*S1	*/
		aj1p04r = cr2-si2;		aj1p04i = ci2+sr2;	/* X2 = C2 + I*S2	*/
		aj1p06r = cr3-si3;		aj1p06i = ci3+sr3;	/* X3 = C3 + I*S3	*/
		aj1p08r = cr4-si4;		aj1p08i = ci4+sr4;	/* X4 = C4 + I*S4	*/
		aj1p10r = cr5-si5;		aj1p10i = ci5+sr5;	/* X5 = C5 + I*S5	*/
		aj1p12r = cr5+si5;		aj1p12i = ci5-sr5;	/* X6 =	C5 - I*S5	*/
		aj1p14r = cr4+si4;		aj1p14i = ci4-sr4;	/* X7 =	C4 - I*S4	*/
		aj1p16r = cr3+si3;		aj1p16i = ci3-sr3;	/* X8 =	C3 - I*S3	*/
		aj1p18r = cr2+si2;		aj1p18i = ci2-sr2;	/* X9 =	C2 - I*S2	*/
		aj1p20r = cr1+si1;		aj1p20i = ci1-sr1;	/* X10=	C1 - I*S1	*/

/*...Second radix-11 block uses aj1p[1:21:2] as inputs:	*/

		t01 = a[j1+p11];							t02 = a[j2+p11];		/* x0		*/
		t03 = a[j1+p09]+a[j1+p13];				t04 = a[j2+p09]+a[j2+p13];	/* x1 + x10	*/
		t05 = a[j1+p07]+a[j1+p15];				t06 = a[j2+p07]+a[j2+p15];	/* x2 + x9	*/
		t07 = a[j1+p05]+a[j1+p17];				t08 = a[j2+p05]+a[j2+p17];	/* x3 + x8	*/
		t09 = a[j1+p03]+a[j1+p19];				t10 = a[j2+p03]+a[j2+p19];	/* x4 + x7	*/
		t11 = a[j1+p01]+a[j1+p21];				t12 = a[j2+p01]+a[j2+p21];	/* x5 + x6	*/
		t13 = a[j1+p01]-a[j1+p21];				t14 = a[j2+p01]-a[j2+p21];	/* x5 - x6	*/
		t15 = a[j1+p03]-a[j1+p19];				t16 = a[j2+p03]-a[j2+p19];	/* x4 - x7	*/
		t17 = a[j1+p05]-a[j1+p17];				t18 = a[j2+p05]-a[j2+p17];	/* x3 - x8	*/
		t19 = a[j1+p07]-a[j1+p15];				t20 = a[j2+p07]-a[j2+p15];	/* x2 - x9	*/
		t21 = a[j1+p09]-a[j1+p13];				t22 = a[j2+p09]-a[j2+p13];	/* x1 - x10	*/

		aj1p01r = t01+t03+t05+t07+t09+t11;			aj1p01i = t02+t04+t06+t08+t10+t12;	/* X0	*/

		cr1= t01+cc1*t03+cc2*t05+cc3*t07+cc4*t09+cc5*t11;	ci1= t02+cc1*t04+cc2*t06+cc3*t08+cc4*t10+cc5*t12;	/* C1	*/
		cr2= t01+cc2*t03+cc4*t05+cc5*t07+cc3*t09+cc1*t11;	ci2= t02+cc2*t04+cc4*t06+cc5*t08+cc3*t10+cc1*t12;	/* C2	*/
		cr3= t01+cc3*t03+cc5*t05+cc2*t07+cc1*t09+cc4*t11;	ci3= t02+cc3*t04+cc5*t06+cc2*t08+cc1*t10+cc4*t12;	/* C3	*/
		cr4= t01+cc4*t03+cc3*t05+cc1*t07+cc5*t09+cc2*t11;	ci4= t02+cc4*t04+cc3*t06+cc1*t08+cc5*t10+cc2*t12;	/* C4	*/
		cr5= t01+cc5*t03+cc1*t05+cc4*t07+cc2*t09+cc3*t11;	ci5= t02+cc5*t04+cc1*t06+cc4*t08+cc2*t10+cc3*t12;	/* C5	*/

		sr1=     ss1*t21+ss2*t19+ss3*t17+ss4*t15+ss5*t13;	si1=     ss1*t22+ss2*t20+ss3*t18+ss4*t16+ss5*t14;	/* S1	*/
		sr2=     ss2*t21+ss4*t19-ss5*t17-ss3*t15-ss1*t13;	si2=     ss2*t22+ss4*t20-ss5*t18-ss3*t16-ss1*t14;	/* S2	*/
		sr3=     ss3*t21-ss5*t19-ss2*t17+ss1*t15+ss4*t13;	si3=     ss3*t22-ss5*t20-ss2*t18+ss1*t16+ss4*t14;	/* S3	*/
		sr4=     ss4*t21-ss3*t19+ss1*t17+ss5*t15-ss2*t13;	si4=     ss4*t22-ss3*t20+ss1*t18+ss5*t16-ss2*t14;	/* S4	*/
		sr5=     ss5*t21-ss1*t19+ss4*t17-ss2*t15+ss3*t13;	si5=     ss5*t22-ss1*t20+ss4*t18-ss2*t16+ss3*t14;	/* S5	*/

/*...Inline multiply of sine parts by +-I into finishing phase...	*/

		aj1p03r = cr1-si1;		aj1p03i = ci1+sr1;	/* X1 = C1 + I*S1	*/
		aj1p05r = cr2-si2;		aj1p05i = ci2+sr2;	/* X2 = C2 + I*S2	*/
		aj1p07r = cr3-si3;		aj1p07i = ci3+sr3;	/* X3 = C3 + I*S3	*/
		aj1p09r = cr4-si4;		aj1p09i = ci4+sr4;	/* X4 = C4 + I*S4	*/
		aj1p11r = cr5-si5;		aj1p11i = ci5+sr5;	/* X5 = C5 + I*S5	*/
		aj1p13r = cr5+si5;		aj1p13i = ci5-sr5;	/* X6 =	C5 - I*S5	*/
		aj1p15r = cr4+si4;		aj1p15i = ci4-sr4;	/* X7 =	C4 - I*S4	*/
		aj1p17r = cr3+si3;		aj1p17i = ci3-sr3;	/* X8 =	C3 - I*S3	*/
		aj1p19r = cr2+si2;		aj1p19i = ci2-sr2;	/* X9 =	C2 - I*S2	*/
		aj1p21r = cr1+si1;		aj1p21i = ci1-sr1;	/* X10=	C1 - I*S1	*/

/*...and now do 11 radix-2 transforms:	*/

		a[j1    ]=aj1p00r+aj1p01r;	a[j2    ]=aj1p00i+aj1p01i;
		a[j1+p01]=aj1p00r-aj1p01r;	a[j2+p01]=aj1p00i-aj1p01i;

		a[j1+p20]=aj1p02r+aj1p03r;	a[j2+p20]=aj1p02i+aj1p03i;
		a[j1+p21]=aj1p02r-aj1p03r;	a[j2+p21]=aj1p02i-aj1p03i;

		a[j1+p19]=aj1p04r+aj1p05r;	a[j2+p19]=aj1p04i+aj1p05i;
		a[j1+p18]=aj1p04r-aj1p05r;	a[j2+p18]=aj1p04i-aj1p05i;

		a[j1+p16]=aj1p06r+aj1p07r;	a[j2+p16]=aj1p06i+aj1p07i;
		a[j1+p17]=aj1p06r-aj1p07r;	a[j2+p17]=aj1p06i-aj1p07i;

		a[j1+p15]=aj1p08r+aj1p09r;	a[j2+p15]=aj1p08i+aj1p09i;
		a[j1+p14]=aj1p08r-aj1p09r;	a[j2+p14]=aj1p08i-aj1p09i;

		a[j1+p12]=aj1p10r+aj1p11r;	a[j2+p12]=aj1p10i+aj1p11i;
		a[j1+p13]=aj1p10r-aj1p11r;	a[j2+p13]=aj1p10i-aj1p11i;

		a[j1+p11]=aj1p12r+aj1p13r;	a[j2+p11]=aj1p12i+aj1p13i;
		a[j1+p10]=aj1p12r-aj1p13r;	a[j2+p10]=aj1p12i-aj1p13i;

		a[j1+p08]=aj1p14r+aj1p15r;	a[j2+p08]=aj1p14i+aj1p15i;
		a[j1+p09]=aj1p14r-aj1p15r;	a[j2+p09]=aj1p14i-aj1p15i;

		a[j1+p07]=aj1p16r+aj1p17r;	a[j2+p07]=aj1p16i+aj1p17i;
		a[j1+p06]=aj1p16r-aj1p17r;	a[j2+p06]=aj1p16i-aj1p17i;

		a[j1+p04]=aj1p18r+aj1p19r;	a[j2+p04]=aj1p18i+aj1p19i;
		a[j1+p05]=aj1p18r-aj1p19r;	a[j2+p05]=aj1p18i-aj1p19i;

		a[j1+p03]=aj1p20r+aj1p21r;	a[j2+p03]=aj1p20i+aj1p21i;
		a[j1+p02]=aj1p20r-aj1p21r;	a[j2+p02]=aj1p20i-aj1p21i;

	}
}

/*
   0 (   0)  117.0000000000    96.0000000000   117.0000000000    96.0000000000
   1 (  11)  -37.0000000000    -8.0000000000   -37.0000000000    -8.0000000000
   2 (   1)   -8.6477389971     7.2270730891    -8.6477389971     7.2270730891
   3 (  12)   -3.0223004569    23.7316133212    -3.0223004569    23.7316133212
   4 (   2)  -10.1405228139   -23.6488721056   -10.1405228139   -23.6488721056
   5 (  13)   -8.5921367524    21.6413258456    -8.5921367524    21.6413258456
   6 (   3)    7.6869203188   -13.4632968675     7.6869203188   -13.4632968675
   7 (  14)   -8.1246193995    16.7880694327    -8.1246193995    16.7880694327
   8 (   4)    5.5042888214    -5.1970405614     5.5042888214    -5.1970405614
   9 (  15)   -0.5272180707    -5.1403494614    -0.5272180707    -5.1403494614
  10 (   5)   -7.1799286117    -1.3369247201    -7.1799286117    -1.3369247201
  11 (  16)   -6.6000045704    -2.9051777499    -6.6000045704    -2.9051777499
  12 (   6)    5.1967297167   -25.3180414893     5.1967297167   -25.3180414893
  13 (  17)    5.7888188888   -20.3540605997     5.7888188888   -20.3540605997
  14 (   7)   13.5656531024   -15.3886651675    13.5656531024   -15.3886651675
  15 (  18)    2.5054644265   -11.6289177615     2.5054644265   -11.6289177615
  16 (   8)  -17.2478597963    11.6031662658   -17.2478597963    11.6031662658
  17 (  19)   -3.4317337204     0.6891503972    -3.4317337204     0.6891503972
  18 (   9)    9.0096529655     0.2134979599     9.0096529655     0.2134979599
  19 (  20)    1.5843232103    -3.1680925873     1.5843232103    -3.1680925873
  20 (  10)   12.3445008621   -21.2567067648    12.3445008621   -21.2567067648
  21 (  21)   -3.6722891232     0.9122495245    -3.6722891232     0.9122495245
*/

/***************/

void radix22_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-22 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix11_dif_pass for details on the radix-11 subtransforms.
*/
	int j,j1,j2;
	static int n22,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21, first_entry=TRUE;
	static double	cc1 =  0.84125353283118116886,	/* Real part of exp(i*2*pi/11), the radix-11 fundamental sincos datum	*/
			ss1 =  0.54064081745559758210,	/* Imag part of exp(i*2*pi/11).	*/
			cc2 =  0.41541501300188642553,	/* cos(2u)	*/
			ss2 =  0.90963199535451837140,	/* sin(2u)	*/
			cc3 = -0.14231483827328514043,	/* cos(3u)	*/
			ss3 =  0.98982144188093273237,	/* sin(3u)	*/
			cc4 = -0.65486073394528506404,	/* cos(4u)	*/
			ss4 =  0.75574957435425828378,	/* sin(4u)	*/
			cc5 = -0.95949297361449738988,	/* cos(5u)	*/
			ss5 =  0.28173255684142969773;	/* sin(5u)	*/
	double  rt,it
		,cr1,cr2,cr3,cr4,cr5,ci1,ci2,ci3,ci4,ci5
		,sr1,sr2,sr3,sr4,sr5,si1,si2,si3,si4,si5
		,u01,u02,u03,u04,u05,u06,u07,u08,u09,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,u31,u32,u33,u34,u35,u36,u37,u38,u39,u40,u41,u42,u43,u44
		,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22;

	if(!first_entry && (n/22) != n22)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n22=n/22;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n22;
		p02 = p01 + p01;
		p03 = p02 + p01;
		p04 = p03 + p01;
		p05 = p04 + p01;
		p06 = p05 + p01;
		p07 = p06 + p01;
		p08 = p07 + p01;
		p09 = p08 + p01;
		p10 = p09 + p01;
		p11 = p10 + p01;
		p12 = p11 + p01;
		p13 = p12 + p01;
		p14 = p13 + p01;
		p15 = p14 + p01;
		p16 = p15 + p01;
		p17 = p16 + p01;
		p18 = p17 + p01;
		p19 = p18 + p01;
		p20 = p19 + p01;
		p21 = p20 + p01;

		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p17 = p17 + ( (p17 >> DAT_BITS) << PAD_BITS );
		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p19 = p19 + ( (p19 >> DAT_BITS) << PAD_BITS );
		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p21 = p21 + ( (p21 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-22 pass is here.	*/

	for(j=0; j < n22; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*...gather the needed data (22 64-bit complex, i.e. 44 64-bit reals) and do 11 radix-2 transforms,	*/
	/*
		Twiddleless version requires us to swap inputs as follows:
		indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21
			  => 0,20,18,16,14,12,10, 8, 6, 4, 2,11, 9, 7, 5, 3, 1,21,19,17,15,13 modulo 22.
		I.e. start out with first pair of indices {0,11}, permute those according to
		{0,11}*21%22 = {0,11}, then each is head of a length-11 list of indices with decrement 2.

		Remember, inputs to DIT are bit-reversed, so
		a[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21] contain
		x[0,11, 1,12, 2,13, 3,14, 4,15, 5,16, 6,17, 7,18, 8,19, 9,20,10,21], which get swapped to
		x[0,11,20, 9,18, 7,16, 5,14, 3,12, 1,10,21, 8,19, 6,17, 4,15, 2,13], which means the a-indices get swapped as
		a[0, 1,19,18,15,14,11,10, 7, 6, 3, 2,20,21,16,17,12,13, 8, 9, 4, 5].
	*/
		u01=a[j1    ];			u02=a[j2    ];	/* x0a,b */
		rt =a[j1+p01];			it =a[j2+p01];
		u03=u01-rt;				u04=u02-it;
		u01=u01+rt;				u02=u02+it;

		u05=a[j1+p19];			u06=a[j2+p19];	/* x1a,b */
		rt =a[j1+p18];			it =a[j2+p18];
		u07=u05-rt;				u08=u06-it;
		u05=u05+rt;				u06=u06+it;

		u09=a[j1+p15];			u10=a[j2+p15];	/* x2a,b */
		rt =a[j1+p14];			it =a[j2+p14];
		u11=u09-rt;				u12=u10-it;
		u09=u09+rt;				u10=u10+it;

		u13=a[j1+p11];			u14=a[j2+p11];	/* x3a,b */
		rt =a[j1+p10];			it =a[j2+p10];
		u15=u13-rt;				u16=u14-it;
		u13=u13+rt;				u14=u14+it;

		u17=a[j1+p07];			u18=a[j2+p07];	/* x4a,b */
		rt =a[j1+p06];			it =a[j2+p06];
		u19=u17-rt;				u20=u18-it;
		u17=u17+rt;				u18=u18+it;

		u21=a[j1+p03];			u22=a[j2+p03];	/* x5a,b */
		rt =a[j1+p02];			it =a[j2+p02];
		u23=u21-rt;				u24=u22-it;
		u21=u21+rt;				u22=u22+it;

		u25=a[j1+p20];			u26=a[j2+p20];	/* x6a,b */
		rt =a[j1+p21];			it =a[j2+p21];
		u27=u25-rt;				u28=u26-it;
		u25=u25+rt;				u26=u26+it;

		u29=a[j1+p16];			u30=a[j2+p16];	/* x7a,b */
		rt =a[j1+p17];			it =a[j2+p17];
		u31=u29-rt;				u32=u30-it;
		u29=u29+rt;				u30=u30+it;

		u33=a[j1+p12];			u34=a[j2+p12];	/* x8a,b */
		rt =a[j1+p13];			it =a[j2+p13];
		u35=u33-rt;				u36=u34-it;
		u33=u33+rt;				u34=u34+it;

		u37=a[j1+p08];			u38=a[j2+p08];	/* x9a,b */
		rt =a[j1+p09];			it =a[j2+p09];
		u39=u37-rt;				u40=u38-it;
		u37=u37+rt;				u38=u38+it;

		u41=a[j1+p04];			u42=a[j2+p04];	/* x10a,b */
		rt =a[j1+p05];			it =a[j2+p05];
		u43=u41-rt;				u44=u42-it;
		u41=u41+rt;				u42=u42+it;

/*       ...and now do two radix-11 transforms.	*/

/*...aj1p[0:20:2]r use u[1:41:4]; aj1p[0:12:2]i use u[2:42:4]	*/

		t01= u01;				t02= u02;		/* x0		*/

		t03= u05 + u41;			t04= u06 + u42;		/* x1 + x10	*/
		t21= u05 - u41;			t22= u06 - u42;		/* x1 - x10	*/

		t05= u09 + u37;			t06= u10 + u38;		/* x2 + x9	*/
		t19= u09 - u37;			t20= u10 - u38;		/* x2 - x9	*/

		t07= u13 + u33;			t08= u14 + u34;		/* x3 + x8	*/
		t17= u13 - u33;			t18= u14 - u34;		/* x3 - x8	*/

		t09= u17 + u29;			t10= u18 + u30;		/* x4 + x7	*/
		t15= u17 - u29;			t16= u18 - u30;		/* x4 - x7	*/

		t11= u21 + u25;			t12= u22 + u26;		/* x5 + x6	*/
		t13= u21 - u25;			t14= u22 - u26;		/* x5 - x6	*/

		a[j1   ] = t01+t03+t05+t07+t09+t11;	a[j2] = t02+t04+t06+t08+t10+t12;	/* X0	*/

		cr1= t01+cc1*t03+cc2*t05+cc3*t07+cc4*t09+cc5*t11;	ci1= t02+cc1*t04+cc2*t06+cc3*t08+cc4*t10+cc5*t12;	/* C1	*/
		cr2= t01+cc2*t03+cc4*t05+cc5*t07+cc3*t09+cc1*t11;	ci2= t02+cc2*t04+cc4*t06+cc5*t08+cc3*t10+cc1*t12;	/* C2	*/
		cr3= t01+cc3*t03+cc5*t05+cc2*t07+cc1*t09+cc4*t11;	ci3= t02+cc3*t04+cc5*t06+cc2*t08+cc1*t10+cc4*t12;	/* C3	*/
		cr4= t01+cc4*t03+cc3*t05+cc1*t07+cc5*t09+cc2*t11;	ci4= t02+cc4*t04+cc3*t06+cc1*t08+cc5*t10+cc2*t12;	/* C4	*/
		cr5= t01+cc5*t03+cc1*t05+cc4*t07+cc2*t09+cc3*t11;	ci5= t02+cc5*t04+cc1*t06+cc4*t08+cc2*t10+cc3*t12;	/* C5	*/

		sr1=     ss1*t21+ss2*t19+ss3*t17+ss4*t15+ss5*t13;	si1=     ss1*t22+ss2*t20+ss3*t18+ss4*t16+ss5*t14;	/* S1	*/
		sr2=     ss2*t21+ss4*t19-ss5*t17-ss3*t15-ss1*t13;	si2=     ss2*t22+ss4*t20-ss5*t18-ss3*t16-ss1*t14;	/* S2	*/
		sr3=     ss3*t21-ss5*t19-ss2*t17+ss1*t15+ss4*t13;	si3=     ss3*t22-ss5*t20-ss2*t18+ss1*t16+ss4*t14;	/* S3	*/
		sr4=     ss4*t21-ss3*t19+ss1*t17+ss5*t15-ss2*t13;	si4=     ss4*t22-ss3*t20+ss1*t18+ss5*t16-ss2*t14;	/* S4	*/
		sr5=     ss5*t21-ss1*t19+ss4*t17-ss2*t15+ss3*t13;	si5=     ss5*t22-ss1*t20+ss4*t18-ss2*t16+ss3*t14;	/* S5	*/

/*...Inline multiply of sine parts by +-I into finishing phase...	*/

		a[j1+p10]=cr1+si1;		a[j2+p10]=ci1-sr1;	/* X1 = C1 + I*S1	*/
		a[j1+p20]=cr2+si2;		a[j2+p20]=ci2-sr2;	/* X2 = C2 + I*S2	*/
		a[j1+p08]=cr3+si3;		a[j2+p08]=ci3-sr3;	/* X3 = C3 + I*S3	*/
		a[j1+p18]=cr4+si4;		a[j2+p18]=ci4-sr4;	/* X4 = C4 + I*S4	*/
		a[j1+p06]=cr5+si5;		a[j2+p06]=ci5-sr5;	/* X5 = C5 + I*S5	*/
		a[j1+p16]=cr5-si5;		a[j2+p16]=ci5+sr5;	/* X6 =	C5 - I*S5	*/
		a[j1+p04]=cr4-si4;		a[j2+p04]=ci4+sr4;	/* X7 =	C4 - I*S4	*/
		a[j1+p14]=cr3-si3;		a[j2+p14]=ci3+sr3;	/* X8 =	C3 - I*S3	*/
		a[j1+p02]=cr2-si2;		a[j2+p02]=ci2+sr2;	/* X9 =	C2 - I*S2	*/
		a[j1+p12]=cr1-si1;		a[j2+p12]=ci1+sr1;	/* X10=	C1 - I*S1	*/

/*...aj1p[1:21:2]r use u[3:43:4]; aj1p[1:13:2]i use u[4:44:4]	*/

		t01= u03;				t02= u04;		/* x0		*/

		t03= u07 + u43;			t04= u08 + u44;		/* x1 + x10	*/
		t21= u07 - u43;			t22= u08 - u44;		/* x1 - x10	*/

		t05= u11 + u39;			t06= u12 + u40;		/* x2 + x9	*/
		t19= u11 - u39;			t20= u12 - u40;		/* x2 - x9	*/

		t07= u15 + u35;			t08= u16 + u36;		/* x3 + x8	*/
		t17= u15 - u35;			t18= u16 - u36;		/* x3 - x8	*/

		t09= u19 + u31;			t10= u20 + u32;		/* x4 + x7	*/
		t15= u19 - u31;			t16= u20 - u32;		/* x4 - x7	*/

		t11= u23 + u27;			t12= u24 + u28;		/* x5 + x6	*/
		t13= u23 - u27;			t14= u24 - u28;		/* x5 - x6	*/

		a[j1+p11] = t01+t03+t05+t07+t09+t11;	a[j2+p11] = t02+t04+t06+t08+t10+t12;	/* X0	*/

		cr1= t01+cc1*t03+cc2*t05+cc3*t07+cc4*t09+cc5*t11;	ci1= t02+cc1*t04+cc2*t06+cc3*t08+cc4*t10+cc5*t12;	/* C1	*/
		cr2= t01+cc2*t03+cc4*t05+cc5*t07+cc3*t09+cc1*t11;	ci2= t02+cc2*t04+cc4*t06+cc5*t08+cc3*t10+cc1*t12;	/* C2	*/
		cr3= t01+cc3*t03+cc5*t05+cc2*t07+cc1*t09+cc4*t11;	ci3= t02+cc3*t04+cc5*t06+cc2*t08+cc1*t10+cc4*t12;	/* C3	*/
		cr4= t01+cc4*t03+cc3*t05+cc1*t07+cc5*t09+cc2*t11;	ci4= t02+cc4*t04+cc3*t06+cc1*t08+cc5*t10+cc2*t12;	/* C4	*/
		cr5= t01+cc5*t03+cc1*t05+cc4*t07+cc2*t09+cc3*t11;	ci5= t02+cc5*t04+cc1*t06+cc4*t08+cc2*t10+cc3*t12;	/* C5	*/

		sr1=     ss1*t21+ss2*t19+ss3*t17+ss4*t15+ss5*t13;	si1=     ss1*t22+ss2*t20+ss3*t18+ss4*t16+ss5*t14;	/* S1	*/
		sr2=     ss2*t21+ss4*t19-ss5*t17-ss3*t15-ss1*t13;	si2=     ss2*t22+ss4*t20-ss5*t18-ss3*t16-ss1*t14;	/* S2	*/
		sr3=     ss3*t21-ss5*t19-ss2*t17+ss1*t15+ss4*t13;	si3=     ss3*t22-ss5*t20-ss2*t18+ss1*t16+ss4*t14;	/* S3	*/
		sr4=     ss4*t21-ss3*t19+ss1*t17+ss5*t15-ss2*t13;	si4=     ss4*t22-ss3*t20+ss1*t18+ss5*t16-ss2*t14;	/* S4	*/
		sr5=     ss5*t21-ss1*t19+ss4*t17-ss2*t15+ss3*t13;	si5=     ss5*t22-ss1*t20+ss4*t18-ss2*t16+ss3*t14;	/* S5	*/

/*...Inline multiply of sine parts by +-I into finishing phase...	*/

		a[j1+p21]=cr1+si1;		a[j2+p21]=ci1-sr1;	/* X1 = C1 + I*S1	*/
		a[j1+p09]=cr2+si2;		a[j2+p09]=ci2-sr2;	/* X2 = C2 + I*S2	*/
		a[j1+p19]=cr3+si3;		a[j2+p19]=ci3-sr3;	/* X3 = C3 + I*S3	*/
		a[j1+p07]=cr4+si4;		a[j2+p07]=ci4-sr4;	/* X4 = C4 + I*S4	*/
		a[j1+p17]=cr5+si5;		a[j2+p17]=ci5-sr5;	/* X5 = C5 + I*S5	*/
		a[j1+p05]=cr5-si5;		a[j2+p05]=ci5+sr5;	/* X6 =	C5 - I*S5	*/
		a[j1+p15]=cr4-si4;		a[j2+p15]=ci4+sr4;	/* X7 =	C4 - I*S4	*/
		a[j1+p03]=cr3-si3;		a[j2+p03]=ci3+sr3;	/* X8 =	C3 - I*S3	*/
		a[j1+p13]=cr2-si2;		a[j2+p13]=ci2+sr2;	/* X9 =	C2 - I*S2	*/
		a[j1+p01]=cr1-si1;		a[j2+p01]=ci1+sr1;	/* X10=	C1 - I*S1	*/

	}
}

