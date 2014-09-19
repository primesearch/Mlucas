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

#define RADIX 30	// Use #define rather than const int to ensure it's really a compile-time const in the C sense
#define ODD_RADIX 15	// ODD_RADIX = [radix >> trailz(radix)]

/**************/

int radix30_ditN_cy_dif1(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter, double *fracmax, uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-30 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-30 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix7/8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int NDIVR,i,j,j1,j2,jt,jp,jstart,jhi,iroot,root_incr,k1,k2,k,khi,l,ntmp,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29;
	static int poff[(RADIX+2)>>2];
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double radix_inv, n2inv;
	const double c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
				s   =  0.86602540378443864675,	/* sin(twopi/3)		*/
				cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
				cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
				ss3 =  0.95105651629515357210,	/*  sin(u) */
				sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
				sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
	,u01,u02,u03,u04,u05,u06,u07,u08,u09,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,u31,u32,u33,u34,u35,u36,u37,u38,u39,u40,u41,u42,u43,u44,u45,u46,u47,u48,u49,u50,u51,u52,u53,u54,u55,u56,u57,u58,u59,u60
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
	,t00r,t01r,t02r,t03r,t04r,t05r,t06r,t07r,t08r,t09r,t10r,t11r,t12r,t13r,t14r,t15r,t16r,t17r,t18r,t19r,t20r,t21r,t22r,t23r,t24r,t25r,t26r,t27r,t28r,t29r
	,t00i,t01i,t02i,t03i,t04i,t05i,t06i,t07i,t08i,t09i,t10i,t11i,t12i,t13i,t14i,t15i,t16i,t17i,t18i,t19i,t20i,t21i,t22i,t23i,t24i,t25i,t26i,t27i,t28i,t29i;
	int bjmodn[RADIX];
	double temp,frac,
		cy_r[RADIX],cy_i[RADIX];
	double scale,dtmp, maxerr = 0.0;
	double *addr,*addi;
	int *itmp;	// Pointer into the bjmodn array
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	/* indices into weights arrays (mod NWT) */
	int ic,ii[RADIX] = {
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
		-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
	};
	double wt_re,wt_im;		/* Fermat-mod weights stuff */

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=-1;

/*...change NDIVR and n_div_wt to non-static to work around a gcc compiler bug. */
	NDIVR   = n/RADIX;
	n_div_nwt = NDIVR >> nwt_bits;

	if((n_div_nwt << nwt_bits) != NDIVR)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/30 in radix30_ditN_cy_dif1.\n",iter);
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
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)RADIX));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

		/*   constant index offsets for array load/stores are here.	*/

		p01 = NDIVR;
		p02 = p01 + NDIVR;		p01 += ( (p01 >> DAT_BITS) << PAD_BITS );
		p03 = p02 + NDIVR;		p02 += ( (p02 >> DAT_BITS) << PAD_BITS );
		p04 = p03 + NDIVR;		p03 += ( (p03 >> DAT_BITS) << PAD_BITS );
		p05 = p04 + NDIVR;		p04 += ( (p04 >> DAT_BITS) << PAD_BITS );
		p06 = p05 + NDIVR;		p05 += ( (p05 >> DAT_BITS) << PAD_BITS );
		p07 = p06 + NDIVR;		p06 += ( (p06 >> DAT_BITS) << PAD_BITS );
		p08 = p07 + NDIVR;		p07 += ( (p07 >> DAT_BITS) << PAD_BITS );
		p09 = p08 + NDIVR;		p08 += ( (p08 >> DAT_BITS) << PAD_BITS );
		p10 = p09 + NDIVR;		p09 += ( (p09 >> DAT_BITS) << PAD_BITS );
		p11 = p10 + NDIVR;		p10 += ( (p10 >> DAT_BITS) << PAD_BITS );
		p12 = p11 + NDIVR;		p11 += ( (p11 >> DAT_BITS) << PAD_BITS );
		p13 = p12 + NDIVR;		p12 += ( (p12 >> DAT_BITS) << PAD_BITS );
		p14 = p13 + NDIVR;		p13 += ( (p13 >> DAT_BITS) << PAD_BITS );
		p15 = p14 + NDIVR;		p14 += ( (p14 >> DAT_BITS) << PAD_BITS );
		p16 = p15 + NDIVR;		p15 += ( (p15 >> DAT_BITS) << PAD_BITS );
		p17 = p16 + NDIVR;		p16 += ( (p16 >> DAT_BITS) << PAD_BITS );
		p18 = p17 + NDIVR;		p17 += ( (p17 >> DAT_BITS) << PAD_BITS );
		p19 = p18 + NDIVR;		p18 += ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 = p19 + NDIVR;		p19 += ( (p19 >> DAT_BITS) << PAD_BITS );
		p21 = p20 + NDIVR;		p20 += ( (p20 >> DAT_BITS) << PAD_BITS );
		p22 = p21 + NDIVR;		p21 += ( (p21 >> DAT_BITS) << PAD_BITS );
		p23 = p22 + NDIVR;		p22 += ( (p22 >> DAT_BITS) << PAD_BITS );
		p24 = p23 + NDIVR;		p23 += ( (p23 >> DAT_BITS) << PAD_BITS );
		p25 = p24 + NDIVR;		p24 += ( (p24 >> DAT_BITS) << PAD_BITS );
		p26 = p25 + NDIVR;		p25 += ( (p25 >> DAT_BITS) << PAD_BITS );
		p27 = p26 + NDIVR;		p26 += ( (p26 >> DAT_BITS) << PAD_BITS );
		p28 = p27 + NDIVR;		p27 += ( (p27 >> DAT_BITS) << PAD_BITS );
		p29 = p28 + NDIVR;		p28 += ( (p28 >> DAT_BITS) << PAD_BITS );
								p29 += ( (p29 >> DAT_BITS) << PAD_BITS );

		poff[0x0] =   0; poff[0x1] = p04; poff[0x2] = p08; poff[0x3] = p12;
		poff[0x4] = p16; poff[0x5] = p20; poff[0x6] = p24; poff[0x7] = p28;

		/* For Fermat-mod, since 'adjacent' words are actually stride-2 separated
		in terms of the floating residue array, block boundaries have half the i-index
		(e.g. as in sw*i and bw*i) value they do in the Mersenne-mod case:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			bjmodnini=0;
			for(j=0; j < NDIVR; j++)
			{
				bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
			}
		}
		else
		{
			bjmodnini=0;
			for(j=0; j < NDIVR/2; j++)
			{
				bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
			}
		}
	}

/*...The radix-30 final DIT pass is here.	*/

	for(l = 0; l < RADIX; l++) {
		cy_r[l] = cy_i[l] = 0;
	}
	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		cy_r[0] = -2;
	}

	*fracmax=0;	/* init max. fractional error	*/

	iroot = 0;	/* init sincos array index	*/
	root_incr = 1;	/* init sincos array index increment (set = 1 for normal carry pass, = 0 for wrapper pass)	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

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
int full_pass = (root_incr!=0);

	i = 0;		/* Index into the BASE and BASEINV arrays. */
	/* If bw > 0 (i.e. n does not divide p), lowest-order digit always a bigword:	*/
	if(bw > 0)
		i = 1;

	bjmodn[0] = 0;
	bjmodn[1] = bjmodnini;
	j = bjmodnini-n;	// Const addend
	for(l = 2; l < RADIX; l++) {
		MOD_ADD32(bjmodn[l-1], j, n, bjmodn[l]);
	}

	/* For Fermat-mod, IBDWT access patterns repeat with period NWT = {odd part of radix0},
	so for even radix0's only really need that many bjmodn and ii's, but that would require
	specialized carry macros tha`t don't update ii and bjmodn - not worth the trouble.
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		col=0;
		co2=(n >> nwt_bits)-1+RADIX;
		co3=co2-RADIX;		/* At the start of each new j-loop, co3=co2-radix(1)	*/
	}
	else
	{
		/* indices into IBDWT weights arrays (mod NWT) is here: */
		ii[0]= 0;
		ii[1]= (SW_DIV_N*(NDIVR >> 1)) % nwt;	// nwt *not* a power of 2, must use library-mod!
		for(l = 2; l < RADIX; l++) {
			MOD_ADD32(ii[l-1], ii[1], nwt, ii[l]);
		}

		/* Start this value off at N in Fermat-mod case, so (bjmodn >= sw) check in
		fermat_carry_norm_errcheck (cf. carry.h) yields a bigword (i == 1) for j= 0:
		*/
		bjmodn[ 0] = n;
		bjmodn[15] = n;
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
		/*
		!...gather the needed data (30 64-bit complex, i.e. 60 64-bit reals) and do 15 radix-2 transforms...
		*/
			u01=a[j1      ];		u02=a[j2    ];	/* x0a,b */
			rt =a[j1  +p01];		it =a[j2+p01];
			u03=u01-rt;				u04=u02-it;
			u01=u01+rt;				u02=u02+it;

			u05=a[j1  +p03];		u06=a[j2+p03];	/* x1a,b */
			rt =a[j1  +p02];		it =a[j2+p02];
			u07=u05-rt;				u08=u06-it;
			u05=u05+rt;				u06=u06+it;

			u09=a[j1  +p04];		u10=a[j2+p04];	/* x2a,b */
			rt =a[j1  +p05];		it =a[j2+p05];
			u11=u09-rt;				u12=u10-it;
			u09=u09+rt;				u10=u10+it;

			u13=a[j1  +p23];		u14=a[j2+p23];	/* x3a,b */
			rt =a[j1  +p22];		it =a[j2+p22];
			u15=u13-rt;				u16=u14-it;
			u13=u13+rt;				u14=u14+it;

			u17=a[j1  +p19];		u18=a[j2+p19];	/* x4a,b */
			rt =a[j1  +p18];		it =a[j2+p18];
			u19=u17-rt;				u20=u18-it;
			u17=u17+rt;				u18=u18+it;

			u21=a[j1  +p20];		u22=a[j2+p20];	/* x5a,b */
			rt =a[j1  +p21];		it =a[j2+p21];
			u23=u21-rt;				u24=u22-it;
			u21=u21+rt;				u22=u22+it;

			u25=a[j1  +p11];		u26=a[j2+p11];	/* x6a,b */
			rt =a[j1  +p10];		it =a[j2+p10];
			u27=u25-rt;				u28=u26-it;
			u25=u25+rt;				u26=u26+it;

			u29=a[j1  +p07];		u30=a[j2+p07];	/* x7a,b */
			rt =a[j1  +p06];		it =a[j2+p06];
			u31=u29-rt;				u32=u30-it;
			u29=u29+rt;				u30=u30+it;

			u33=a[j1  +p08];		u34=a[j2+p08];	/* x8a,b */
			rt =a[j1  +p09];		it =a[j2+p09];
			u35=u33-rt;				u36=u34-it;
			u33=u33+rt;				u34=u34+it;

			u37=a[j1  +p27];		u38=a[j2+p27];	/* x9a,b */
			rt =a[j1  +p26];		it =a[j2+p26];
			u39=u37-rt;				u40=u38-it;
			u37=u37+rt;				u38=u38+it;

			u41=a[j1  +p28];		u42=a[j2+p28];	/* x10a,b */
			rt =a[j1  +p29];		it =a[j2+p29];
			u43=u41-rt;				u44=u42-it;
			u41=u41+rt;				u42=u42+it;

			u45=a[j1  +p24];		u46=a[j2+p24];	/* x11a,b */
			rt =a[j1  +p25];		it =a[j2+p25];
			u47=u45-rt;				u48=u46-it;
			u45=u45+rt;				u46=u46+it;

			u49=a[j1  +p15];		u50=a[j2+p15];	/* x12a,b */
			rt =a[j1  +p14];		it =a[j2+p14];
			u51=u49-rt;				u52=u50-it;
			u49=u49+rt;				u50=u50+it;

			u53=a[j1  +p16];		u54=a[j2+p16];	/* x13a,b */
			rt =a[j1  +p17];		it =a[j2+p17];
			u55=u53-rt;				u56=u54-it;
			u53=u53+rt;				u54=u54+it;

			u57=a[j1  +p12];		u58=a[j2+p12];	/* x14a,b */
			rt =a[j1  +p13];		it =a[j2+p13];
			u59=u57-rt;				u60=u58-it;
			u57=u57+rt;				u58=u58+it;

		//...and now do two radix-15 transforms:

	/*...a1p[0:28:2]r use u[1:57:4]; a1p[0:28:2]i use u[2:58:4]	*/

	/*   Gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do five radix-3 transforms...	*/
		/*...Block 1:	*/
			t00=u01;			t01=u02;
			t02=u09;			t03=u10;
			rt =u05;			it =u06;
			t04=t02-rt;			t05=t03-it;
			t02=t02+rt;			t03=t03+it;
			t00=t00+t02;		t01=t01+t03;
			t02=t00+c3m1*t02;		t03=t01+c3m1*t03;
			rt =s*t04;			it =s*t05;
			t04=t02-it;			t05=t03+rt;
			t02=t02+it;			t03=t03-rt;

		/*...Block 2:	*/
			t06=u33;			t07=u34;
			t08=u29;			t09=u30;
			rt =u25;			it =u26;
			t10=t08-rt;			t11=t09-it;
			t08=t08+rt;			t09=t09+it;
			t06=t06+t08;		t07=t07+t09;
			t08=t06+c3m1*t08;		t09=t07+c3m1*t09;
			rt =s*t10;			it =s*t11;
			t10=t08-it;			t11=t09+rt;
			t08=t08+it;			t09=t09-rt;

		/*...Block 3:	*/
			t12=u53;			t13=u54;
			t14=u49;			t15=u50;
			rt =u57;			it =u58;
			t16=t14 -rt;		t17=t15 -it;
			t14=t14 +rt;		t15=t15 +it;
			t12=t12+t14;		t13=t13+t15;
			t14=t12+c3m1*t14;		t15=t13+c3m1*t15;
			rt =s*t16;			it =s*t17;
			t16=t14-it;			t17=t15+rt;
			t14=t14+it;			t15=t15-rt;

		/*...Block 4:	*/
			t18=u17;			t19=u18;
			t20=u13;			t21=u14;
			rt =u21;			it =u22;
			t22=t20 -rt;		t23=t21 -it;
			t20=t20 +rt;		t21=t21 +it;
			t18=t18+t20;		t19=t19+t21;
			t20=t18+c3m1*t20;		t21=t19+c3m1*t21;
			rt =s*t22;			it =s*t23;
			t22=t20-it;			t23=t21+rt;
			t20=t20+it;			t21=t21-rt;

		/*...Block 5:	*/
			t24=u37;			t25=u38;
			t26=u45;			t27=u46;
			rt =u41;			it =u42;
			t28=t26 -rt;		t29=t27 -it;
			t26=t26 +rt;		t27=t27 +it;
			t24=t24+t26;		t25=t25+t27;
			t26=t24+c3m1*t26;		t27=t25+c3m1*t27;
			rt =s*t28;			it =s*t29;
			t28=t26-it;			t29=t27+rt;
			t26=t26+it;			t27=t27-rt;

		/*...and now do three radix-5 transforms:	*/
		/*...Block 1:	*/
			rt = t24;			it = t25;
			t24= t06-rt;		t25= t07-it;
			t06= t06+rt;		t07= t07+it;
			rt = t18;			it = t19;
			t18= t12-rt;		t19= t13-it;
			t12= t12+rt;		t13= t13+it;

			rt = t06+t12;		it = t07+t13;
			t00= t00+rt;		t01= t01+it;
			rt = t00+cn1*rt;		it = t01+cn1*it;
			t12= cn2*(t06-t12);		t13= cn2*(t07-t13);
			t06= rt+t12;		t07= it+t13;
			t12= rt-t12;		t13= it-t13;
			rt = ss3*(t18-t24);		it = ss3*(t19-t25);
			t18= rt-sn1*t18;		t19= it-sn1*t19;
			t24= rt+sn2*t24;		t25= it+sn2*t25;

			a[j1    ]=t00;		a[j2    ]=t01;
			a[j1+p06]=t06-t19;		a[j2+p06]=t07+t18;
			a[j1+p12]=t12-t25;		a[j2+p12]=t13+t24;
			a[j1+p18]=t12+t25;		a[j2+p18]=t13-t24;
			a[j1+p24]=t06+t19;		a[j2+p24]=t07-t18;

		/*...Block 2:	*/
			rt = t26;			it = t27;
			t26= t08-rt;		t27= t09-it;
			t08= t08+rt;		t09= t09+it;
			rt = t20;			it = t21;
			t20= t14-rt;		t21= t15-it;
			t14= t14+rt;		t15= t15+it;

			rt = t08+t14;		it = t09+t15;
			t02= t02+rt;		t03= t03+it;
			rt = t02+cn1*rt;		it = t03+cn1*it;
			t14= cn2*(t08-t14);		t15= cn2*(t09-t15);
			t08= rt+t14;		t09= it+t15;
			t14= rt-t14;		t15= it-t15;
			rt = ss3*(t20-t26);		it = ss3*(t21-t27);
			t20= rt-sn1*t20;		t21= it-sn1*t21;
			t26= rt+sn2*t26;		t27= it+sn2*t27;

			a[j1+p10]=t02;		a[j2+p10]=t03;
			a[j1+p16]=t08-t21;		a[j2+p16]=t09+t20;
			a[j1+p22]=t14-t27;		a[j2+p22]=t15+t26;
			a[j1+p28]=t14+t27;		a[j2+p28]=t15-t26;
			a[j1+p04]=t08+t21;		a[j2+p04]=t09-t20;

		/*...Block 3:	*/
			rt = t28;			it = t29;
			t28= t10-rt;		t29= t11-it;
			t10= t10+rt;		t11= t11+it;
			rt = t22;			it = t23;
			t22= t16-rt;		t23= t17-it;
			t16= t16+rt;		t17= t17+it;

			rt = t10+t16;		it = t11+t17;
			t04= t04+rt;		t05= t05+it;
			rt = t04+cn1*rt;		it = t05+cn1*it;
			t16= cn2*(t10-t16);		t17= cn2*(t11-t17);
			t10= rt+t16;		t11= it+t17;
			t16= rt-t16;		t17= it-t17;
			rt = ss3*(t22-t28);		it = ss3*(t23-t29);
			t22= rt-sn1*t22;		t23= it-sn1*t23;
			t28= rt+sn2*t28;		t29= it+sn2*t29;

			a[j1+p20]=t04;		a[j2+p20]=t05;
			a[j1+p26]=t10-t23;		a[j2+p26]=t11+t22;
			a[j1+p02]=t16-t29;		a[j2+p02]=t17+t28;
			a[j1+p08]=t16+t29;		a[j2+p08]=t17-t28;
			a[j1+p14]=t10+t23;		a[j2+p14]=t11-t22;

		/*...a1p[1:29:2]r use u[3:59:4]; a1p[1:29:2]i use u[4:60:4]	*/

		/*   Gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do five radix-3 transforms...	*/
		/*...Block 1:	*/
			t00=u03;			t01=u04;
			t02=u11;			t03=u12;
			rt =u07;			it =u08;
			t04=t02-rt;			t05=t03-it;
			t02=t02+rt;			t03=t03+it;
			t00=t00+t02;		t01=t01+t03;
			t02=t00+c3m1*t02;		t03=t01+c3m1*t03;
			rt =s*t04;			it =s*t05;
			t04=t02-it;			t05=t03+rt;
			t02=t02+it;			t03=t03-rt;

		/*...Block 2:	*/
			t06=u35;			t07=u36;
			t08=u31;			t09=u32;
			rt =u27;			it =u28;
			t10=t08-rt;			t11=t09-it;
			t08=t08+rt;			t09=t09+it;
			t06=t06+t08;		t07=t07+t09;
			t08=t06+c3m1*t08;		t09=t07+c3m1*t09;
			rt =s*t10;			it =s*t11;
			t10=t08-it;			t11=t09+rt;
			t08=t08+it;			t09=t09-rt;

		/*...Block 3:	*/
			t12=u55;			t13=u56;
			t14=u51;			t15=u52;
			rt =u59;			it =u60;
			t16=t14 -rt;		t17=t15 -it;
			t14=t14 +rt;		t15=t15 +it;
			t12=t12+t14;		t13=t13+t15;
			t14=t12+c3m1*t14;		t15=t13+c3m1*t15;
			rt =s*t16;			it =s*t17;
			t16=t14-it;			t17=t15+rt;
			t14=t14+it;			t15=t15-rt;

		/*...Block 4:	*/
			t18=u19;			t19=u20;
			t20=u15;			t21=u16;
			rt =u23;			it =u24;
			t22=t20 -rt;		t23=t21 -it;
			t20=t20 +rt;		t21=t21 +it;
			t18=t18+t20;		t19=t19+t21;
			t20=t18+c3m1*t20;		t21=t19+c3m1*t21;
			rt =s*t22;			it =s*t23;
			t22=t20-it;			t23=t21+rt;
			t20=t20+it;			t21=t21-rt;

		/*...Block 5:	*/
			t24=u39;			t25=u40;
			t26=u47;			t27=u48;
			rt =u43;			it =u44;
			t28=t26 -rt;		t29=t27 -it;
			t26=t26 +rt;		t27=t27 +it;
			t24=t24+t26;		t25=t25+t27;
			t26=t24+c3m1*t26;		t27=t25+c3m1*t27;
			rt =s*t28;			it =s*t29;
			t28=t26-it;			t29=t27+rt;
			t26=t26+it;			t27=t27-rt;

		/*...and now do three radix-5 transforms:	*/
		/*...Block 1:	*/
			rt = t24;			it = t25;
			t24= t06-rt;		t25= t07-it;
			t06= t06+rt;		t07= t07+it;
			rt = t18;			it = t19;
			t18= t12-rt;		t19= t13-it;
			t12= t12+rt;		t13= t13+it;

			rt = t06+t12;		it = t07+t13;
			t00= t00+rt;		t01= t01+it;
			rt = t00+cn1*rt;		it = t01+cn1*it;
			t12= cn2*(t06-t12);		t13= cn2*(t07-t13);
			t06= rt+t12;		t07= it+t13;
			t12= rt-t12;		t13= it-t13;
			rt = ss3*(t18-t24);		it = ss3*(t19-t25);
			t18= rt-sn1*t18;		t19= it-sn1*t19;
			t24= rt+sn2*t24;		t25= it+sn2*t25;

			a[j1+p15]=t00;		a[j2+p15]=t01;
			a[j1+p21]=t06-t19;		a[j2+p21]=t07+t18;
			a[j1+p27]=t12-t25;		a[j2+p27]=t13+t24;
			a[j1+p03]=t12+t25;		a[j2+p03]=t13-t24;
			a[j1+p09]=t06+t19;		a[j2+p09]=t07-t18;

		/*...Block 2:	*/
			rt = t26;			it = t27;
			t26= t08-rt;		t27= t09-it;
			t08= t08+rt;		t09= t09+it;
			rt = t20;			it = t21;
			t20= t14-rt;		t21= t15-it;
			t14= t14+rt;		t15= t15+it;

			rt = t08+t14;		it = t09+t15;
			t02= t02+rt;		t03= t03+it;
			rt = t02+cn1*rt;		it = t03+cn1*it;
			t14= cn2*(t08-t14);		t15= cn2*(t09-t15);
			t08= rt+t14;		t09= it+t15;
			t14= rt-t14;		t15= it-t15;
			rt = ss3*(t20-t26);		it = ss3*(t21-t27);
			t20= rt-sn1*t20;		t21= it-sn1*t21;
			t26= rt+sn2*t26;		t27= it+sn2*t27;

			a[j1+p25]=t02;		a[j2+p25]=t03;
			a[j1+p01]=t08-t21;		a[j2+p01]=t09+t20;
			a[j1+p07]=t14-t27;		a[j2+p07]=t15+t26;
			a[j1+p13]=t14+t27;		a[j2+p13]=t15-t26;
			a[j1+p19]=t08+t21;		a[j2+p19]=t09-t20;

		/*...Block 3:	*/
			rt = t28;			it = t29;
			t28= t10-rt;		t29= t11-it;
			t10= t10+rt;		t11= t11+it;
			rt = t22;			it = t23;
			t22= t16-rt;		t23= t17-it;
			t16= t16+rt;		t17= t17+it;

			rt = t10+t16;		it = t11+t17;
			t04= t04+rt;		t05= t05+it;
			rt = t04+cn1*rt;		it = t05+cn1*it;
			t16= cn2*(t10-t16);		t17= cn2*(t11-t17);
			t10= rt+t16;		t11= it+t17;
			t16= rt-t16;		t17= it-t17;
			rt = ss3*(t22-t28);		it = ss3*(t23-t29);
			t22= rt-sn1*t22;		t23= it-sn1*t23;
			t28= rt+sn2*t28;		t29= it+sn2*t29;

			a[j1+p05]=t04;			a[j2+p05]=t05;
			a[j1+p11]=t10-t23;		a[j2+p11]=t11+t22;
			a[j1+p17]=t16-t29;		a[j2+p17]=t17+t28;
			a[j1+p23]=t16+t29;		a[j2+p23]=t17-t28;
			a[j1+p29]=t10+t23;		a[j2+p29]=t11-t22;

/*...Now do the carries. Since the outputs would
    normally be getting dispatched to 30 separate blocks of the A-array, we need 30 separate carries.	*/

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

			/*...set0 is slightly different from others; divide work into blocks of 4 macro calls, 1st set of which gets pulled out of loop: */
			l = 0; addr = cy_r; itmp = bjmodn;
//if(!j)printf("iter %d, full = %d, a0in = %20.15f\n",iter,full_pass,a[0]/(n>>1));
		   cmplx_carry_norm_errcheck0(a[j1    ],a[j2    ],*addr,*itmp  ); ++l; ++addr; ++itmp;
			// Middle 7 quartets of macro calls done in loop:
			jt = j1; jp = j2;
			for(ntmp = 1; ntmp < 8; ntmp++) {
				cmplx_carry_norm_errcheck(a[jt+p01],a[jp+p01],*addr,*itmp,l); ++l; ++addr; ++itmp;
				cmplx_carry_norm_errcheck(a[jt+p02],a[jp+p02],*addr,*itmp,l); ++l; ++addr; ++itmp;
				cmplx_carry_norm_errcheck(a[jt+p03],a[jp+p03],*addr,*itmp,l); ++l; ++addr; ++itmp;
				jt = j1 + poff[ntmp]; jp = j2 + poff[ntmp];	// poff[] = p04,p08,...,p56
				cmplx_carry_norm_errcheck(a[jt    ],a[jp    ],*addr,*itmp,l); ++l; ++addr; ++itmp;
			}
			cmplx_carry_norm_errcheck(a[j1+p29],a[j2+p29],*addr,*itmp,l);
//if(!j)printf("iter %d, full = %d, a0out = %20.15f\n",iter,full_pass,a[0]);
			i =((uint32)(sw - bjmodn[0]) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/
		}
		else
		{
			// Can't use l as loop index here, since it gets used in the Fermat-mod carry macro (as are k1,k2):
			ntmp = 0; addr = cy_r; addi = cy_i; ic = 0;	// ic = idx into ii-array
			for(m = 0; m < 7; m++) {
				jt = j1 + poff[m]; jp = j2 + poff[m];	// poff[] = p04,p08,...
				fermat_carry_norm_errcheck(a[jt    ],a[jp    ],*addr,*addi,ii[ic],bjmodn[ic],ntmp,NRTM1,NRT_BITS);	ic++; ntmp += NDIVR; ++addr; ++addi;
				fermat_carry_norm_errcheck(a[jt+p01],a[jp+p01],*addr,*addi,ii[ic],bjmodn[ic],ntmp,NRTM1,NRT_BITS);	ic++; ntmp += NDIVR; ++addr; ++addi;
				fermat_carry_norm_errcheck(a[jt+p02],a[jp+p02],*addr,*addi,ii[ic],bjmodn[ic],ntmp,NRTM1,NRT_BITS);	ic++; ntmp += NDIVR; ++addr; ++addi;
				fermat_carry_norm_errcheck(a[jt+p03],a[jp+p03],*addr,*addi,ii[ic],bjmodn[ic],ntmp,NRTM1,NRT_BITS);	ic++; ntmp += NDIVR; ++addr; ++addi;
			}
			fermat_carry_norm_errcheck(a[j1+p28],a[j2+p28],*addr,*addi,ii[ic],bjmodn[ic],ntmp,NRTM1,NRT_BITS);	ic++; ntmp += NDIVR; ++addr; ++addi;
			fermat_carry_norm_errcheck(a[j1+p29],a[j2+p29],*addr,*addi,ii[ic],bjmodn[ic],ntmp,NRTM1,NRT_BITS);
		}

	/*...The radix-30 DIF pass is here:	*/
		/*...gather the needed data (30 64-bit complex, i.e. 60 64-bit reals) and do 2 radix-15 DIF transforms...*/
			RADIX_15_DIF(
				a[j1    ],a[j2    ],a[j1+p28],a[j2+p28],a[j1+p26],a[j2+p26],a[j1+p24],a[j2+p24],a[j1+p22],a[j2+p22],a[j1+p20],a[j2+p20],a[j1+p18],a[j2+p18],a[j1+p16],a[j2+p16],a[j1+p14],a[j2+p14],a[j1+p12],a[j2+p12],a[j1+p10],a[j2+p10],a[j1+p08],a[j2+p08],a[j1+p06],a[j2+p06],a[j1+p04],a[j2+p04],a[j1+p02],a[j2+p02],
				t00r,t00i,t02r,t02i,t04r,t04i,t26r,t26i,t28r,t28i,t24r,t24i,t18r,t18i,t20r,t20i,t22r,t22i,t16r,t16i,t12r,t12i,t14r,t14i,t08r,t08i,t10r,t10i,t06r,t06i
			);
			RADIX_15_DIF(
				a[j1+p15],a[j2+p15],a[j1+p13],a[j2+p13],a[j1+p11],a[j2+p11],a[j1+p09],a[j2+p09],a[j1+p07],a[j2+p07],a[j1+p05],a[j2+p05],a[j1+p03],a[j2+p03],a[j1+p01],a[j2+p01],a[j1+p29],a[j2+p29],a[j1+p27],a[j2+p27],a[j1+p25],a[j2+p25],a[j1+p23],a[j2+p23],a[j1+p21],a[j2+p21],a[j1+p19],a[j2+p19],a[j1+p17],a[j2+p17],
				t01r,t01i,t03r,t03i,t05r,t05i,t27r,t27i,t29r,t29i,t25r,t25i,t19r,t19i,t21r,t21i,t23r,t23i,t17r,t17i,t13r,t13i,t15r,t15i,t09r,t09i,t11r,t11i,t07r,t07i
			);
		/*...and now do 15 radix-2 transforms:	*/
			a[j1    ]=t00r+t01r;	a[j2    ]=t00i+t01i;
			a[j1+p01]=t00r-t01r;	a[j2+p01]=t00i-t01i;
			a[j1+p04]=t02r+t03r;	a[j2+p04]=t02i+t03i;
			a[j1+p05]=t02r-t03r;	a[j2+p05]=t02i-t03i;
			a[j1+p03]=t04r+t05r;	a[j2+p03]=t04i+t05i;
			a[j1+p02]=t04r-t05r;	a[j2+p02]=t04i-t05i;
			a[j1+p07]=t06r+t07r;	a[j2+p07]=t06i+t07i;
			a[j1+p06]=t06r-t07r;	a[j2+p06]=t06i-t07i;
			a[j1+p11]=t08r+t09r;	a[j2+p11]=t08i+t09i;
			a[j1+p10]=t08r-t09r;	a[j2+p10]=t08i-t09i;
			a[j1+p08]=t10r+t11r;	a[j2+p08]=t10i+t11i;
			a[j1+p09]=t10r-t11r;	a[j2+p09]=t10i-t11i;
			a[j1+p15]=t12r+t13r;	a[j2+p15]=t12i+t13i;
			a[j1+p14]=t12r-t13r;	a[j2+p14]=t12i-t13i;
			a[j1+p12]=t14r+t15r;	a[j2+p12]=t14i+t15i;
			a[j1+p13]=t14r-t15r;	a[j2+p13]=t14i-t15i;
			a[j1+p16]=t16r+t17r;	a[j2+p16]=t16i+t17i;
			a[j1+p17]=t16r-t17r;	a[j2+p17]=t16i-t17i;
			a[j1+p23]=t18r+t19r;	a[j2+p23]=t18i+t19i;
			a[j1+p22]=t18r-t19r;	a[j2+p22]=t18i-t19i;
			a[j1+p20]=t20r+t21r;	a[j2+p20]=t20i+t21i;
			a[j1+p21]=t20r-t21r;	a[j2+p21]=t20i-t21i;
			a[j1+p19]=t22r+t23r;	a[j2+p19]=t22i+t23i;
			a[j1+p18]=t22r-t23r;	a[j2+p18]=t22i-t23i;
			a[j1+p24]=t24r+t25r;	a[j2+p24]=t24i+t25i;
			a[j1+p25]=t24r-t25r;	a[j2+p25]=t24i-t25i;
			a[j1+p28]=t26r+t27r;	a[j2+p28]=t26i+t27i;
			a[j1+p29]=t26r-t27r;	a[j2+p29]=t26i-t27i;
			a[j1+p27]=t28r+t29r;	a[j2+p27]=t28i+t29i;
			a[j1+p26]=t28r-t29r;	a[j2+p26]=t28i-t29i;

			iroot += root_incr;		/* increment sincos index.	*/
		}

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			jstart += nwt;
			jhi    += nwt;

			col += RADIX;
			co3 -= RADIX;
		}
	}

	if(full_pass) {
		printf("Iter = %d, maxerr = %20.15f\n",iter,maxerr);
		printf("Carries = ");
		for(l = 0; l < RADIX; l++) {
			printf("%12.5e, ",cy_r[l]);
		}
		printf("\n");
	} else {
		break;
	}

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-30 forward DIF FFT of the first block of 30 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 30 outputs of (1);
!   (3) Reweight and perform a radix-30 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 30 elements and repeat (1-4).
*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t00 = cy_r[RADIX-1];
		for(l = RADIX-1; l > 0; l--) {
			cy_r[l] = cy_r[l-1];
		}
		cy_r[0] = t00;
	} else {
		/* ...The 2 Mo"bius carries are here: */
		t00 = cy_r[RADIX-1]; t01 = cy_i[RADIX-1];
		for(l = RADIX-1; l > 0; l--) {
			cy_r[l] = cy_r[l-1]; cy_i[l] = cy_i[l-1];
		}
		cy_r[0] = -t01; cy_i[0] = +t00;
	}

	iroot = 0;
	root_incr = 0;
	scale = 1;

	jstart = 0;
	/*
	For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
	*/
	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		jhi =15;
	}
	else
	{
		jhi = 7;
	}
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
		a[j+p22] *= radix_inv;
		a[j+p23] *= radix_inv;
		a[j+p24] *= radix_inv;
		a[j+p25] *= radix_inv;
		a[j+p26] *= radix_inv;
		a[j+p27] *= radix_inv;
		a[j+p28] *= radix_inv;
		a[j+p29] *= radix_inv;
	}
}

	dtmp = 0;
	for(l = 0; l < RADIX; l++) {
		dtmp += fabs(cy_r[l]) + fabs(cy_i[l]);
	}
	if(dtmp != 0.0) {
			sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix30_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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

int radix30_ditN_cy_dif1_nochk(double a[], int n, int nwt, int nwt_bits, double wt0[], double wt1[], int si[], struct complex rn0[], struct complex rn1[], double base[], double baseinv[], int iter                 , uint64 p)
{
/*
!...Acronym: DWT = Discrete Weighted Transform, DIT = Decimation In Time, DIF = Decimation In Frequency
!
!...Performs a final radix-30 complex DIT pass, an inverse DWT weighting, a carry propagation,
!   a forward DWT weighting, and an initial radix-30 complex DIF pass on the data in the length-N real vector A.
!
!   Data enter and are returned in the A-array.
!
!   See the documentation in mers_mod_square and radix16_dif_pass for further details on the array
!   storage scheme, and radix7/8_ditN_cy_dif1 for details on the reduced-length weights array scheme.
*/
	int n30, bjmodn00,bjmodn01,bjmodn02,bjmodn03,bjmodn04,bjmodn05,bjmodn06,bjmodn07,bjmodn08,bjmodn09,bjmodn10,bjmodn11,bjmodn12,bjmodn13,bjmodn14,bjmodn15,bjmodn16,bjmodn17,bjmodn18,bjmodn19,bjmodn20,bjmodn21,bjmodn22,bjmodn23,bjmodn24,bjmodn25,bjmodn26,bjmodn27,bjmodn28,bjmodn29
		,i,j,j1,j2,jstart,jhi,iroot,root_incr,k1,k2,k,khi,l,outer;
	static uint64 psave=0;
	static uint32 bw,sw,bjmodnini,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29;
	const  double one_half[3] = {1.0, 0.5, 0.25};	/* Needed for small-weights-tables scheme */
	static double radix_inv, n2inv;
	static double c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/

					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
	,u01,u02,u03,u04,u05,u06,u07,u08,u09,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,u31,u32,u33,u34,u35,u36,u37,u38,u39,u40,u41,u42,u43,u44,u45,u46,u47,u48,u49,u50,u51,u52,u53,u54,u55,u56,u57,u58,u59,u60
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29
	,a1p00r,a1p01r,a1p02r,a1p03r,a1p04r,a1p05r,a1p06r,a1p07r,a1p08r,a1p09r,a1p10r,a1p11r,a1p12r,a1p13r,a1p14r,a1p15r,a1p16r,a1p17r,a1p18r,a1p19r,a1p20r,a1p21r,a1p22r,a1p23r,a1p24r,a1p25r,a1p26r,a1p27r,a1p28r,a1p29r
	,a1p00i,a1p01i,a1p02i,a1p03i,a1p04i,a1p05i,a1p06i,a1p07i,a1p08i,a1p09i,a1p10i,a1p11i,a1p12i,a1p13i,a1p14i,a1p15i,a1p16i,a1p17i,a1p18i,a1p19i,a1p20i,a1p21i,a1p22i,a1p23i,a1p24i,a1p25i,a1p26i,a1p27i,a1p28i,a1p29i
	,cy_r00,cy_r01,cy_r02,cy_r03,cy_r04,cy_r05,cy_r06,cy_r07,cy_r08,cy_r09,cy_r10,cy_r11,cy_r12,cy_r13,cy_r14,cy_r15,cy_r16,cy_r17,cy_r18,cy_r19,cy_r20,cy_r21,cy_r22,cy_r23,cy_r24,cy_r25,cy_r26,cy_r27,cy_r28,cy_r29
	,cy_i00,cy_i01,cy_i02,cy_i03,cy_i04,cy_i05,cy_i06,cy_i07,cy_i08,cy_i09,cy_i10,cy_i11,cy_i12,cy_i13,cy_i14,cy_i15,cy_i16,cy_i17,cy_i18,cy_i19,cy_i20,cy_i21,cy_i22,cy_i23,cy_i24,cy_i25,cy_i26,cy_i27,cy_i28,cy_i29
	,temp,scale;

#if PFETCH
	double *add0, *addr;
#endif
	int err;
	static int first_entry=TRUE;

/*...stuff for the reduced-length DWT weights array is here:	*/
	int n_div_nwt;
	int col,co2,co3,m,m2,n_minus_sil,n_minus_silp1,sinwt,sinwtm1;
	double wt,wtinv,wtl,wtlp1,wtn,wtnm1,wtA,wtB,wtC;	/* Mersenne-mod weights stuff */
	int ii00,ii01,ii02,ii03,ii04,ii05,ii06,ii07,ii08,ii09,ii10,ii11,ii12,ii13,ii14,ii15,ii16,ii17,ii18,ii19,ii20,ii21,ii22,ii23,ii24,ii25,ii26,ii27,ii28,ii29;	/* indices into weights arrays (mod NWT) */
	double wt_re,wt_im;									/* Fermat-mod weights stuff */

	// Init these to get rid of GCC "may be used uninitialized in this function" warnings:
	col=co2=co3=ii00=ii01=ii02=ii03=ii04=ii05=ii06=ii07=ii08=ii09=ii10=ii11=ii12=ii13=ii14=ii15=ii16=ii17=ii18=ii19=ii20=ii21=ii22=ii23=ii24=ii25=ii26=ii27=ii28=ii29=-1;

/*...change n30 and n_div_wt to non-static to work around a gcc compiler bug. */
	n30   = n/30;
	n_div_nwt = n30 >> nwt_bits;

	if((n_div_nwt << nwt_bits) != n30)
	{
		sprintf(cbuf,"FATAL: iter = %10d; NWT_BITS does not divide N/30 in radix30_ditN_cy_dif1.\n",iter);
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
		radix_inv = qfdbl(qf_rational_quotient((int64)1, (int64)30));
		n2inv     = qfdbl(qf_rational_quotient((int64)1, (int64)(n/2)));

		bw    = p%n;		/* Number of bigwords in the Crandall/Fagin mixed-radix representation = (Mersenne exponent) mod (vector length).	*/
		sw    = n - bw;	/* Number of smallwords.	*/

		/*   constant index offsets for array load/stores are here.	*/

		p01 = n30;
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
		p22 = p21 + p01;
		p23 = p22 + p01;
		p24 = p23 + p01;
		p25 = p24 + p01;
		p26 = p25 + p01;
		p27 = p26 + p01;
		p28 = p27 + p01;
		p29 = p28 + p01;

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
		p22 = p22 + ( (p22 >> DAT_BITS) << PAD_BITS );
		p23 = p23 + ( (p23 >> DAT_BITS) << PAD_BITS );
		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p25 = p25 + ( (p25 >> DAT_BITS) << PAD_BITS );
		p26 = p26 + ( (p26 >> DAT_BITS) << PAD_BITS );
		p27 = p27 + ( (p27 >> DAT_BITS) << PAD_BITS );
		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
		p29 = p29 + ( (p29 >> DAT_BITS) << PAD_BITS );

		/* For Fermat-mod, since 'adjacent' words are actually stride-2 separated
		in terms of the floating residue array, block boundaries have half the i-index
		(e.g. as in sw*i and bw*i) value they do in the Mersenne-mod case:
		*/
		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			bjmodnini=0;
			for(j=0; j < n30; j++)
			{
				bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
			}
		}
		else
		{
			bjmodnini=0;
			for(j=0; j < n30/2; j++)
			{
				bjmodnini -= sw; bjmodnini = bjmodnini + ( (-(int)((uint32)bjmodnini >> 31)) & n);
			}
		}
	}

/*...The radix-30 final DIT pass is here.	*/

	cy_r00= 0;	cy_i00= 0;
	cy_r01= 0;	cy_i01= 0;
	cy_r02= 0;	cy_i02= 0;
	cy_r03= 0;	cy_i03= 0;
	cy_r04= 0;	cy_i04= 0;
	cy_r05= 0;	cy_i05= 0;
	cy_r06= 0;	cy_i06= 0;
	cy_r07= 0;	cy_i07= 0;
	cy_r08= 0;	cy_i08= 0;
	cy_r09= 0;	cy_i09= 0;
	cy_r10= 0;	cy_i10= 0;
	cy_r11= 0;	cy_i11= 0;
	cy_r12= 0;	cy_i12= 0;
	cy_r13= 0;	cy_i13= 0;
	cy_r14= 0;	cy_i14= 0;
	cy_r15= 0;	cy_i15= 0;
	cy_r16= 0;	cy_i16= 0;
	cy_r17= 0;	cy_i17= 0;
	cy_r18= 0;	cy_i18= 0;
	cy_r19= 0;	cy_i19= 0;
	cy_r20= 0;	cy_i20= 0;
	cy_r21= 0;	cy_i21= 0;
	cy_r22= 0;	cy_i22= 0;
	cy_r23= 0;	cy_i23= 0;
	cy_r24= 0;	cy_i24= 0;
	cy_r25= 0;	cy_i25= 0;
	cy_r26= 0;	cy_i26= 0;
	cy_r27= 0;	cy_i27= 0;
	cy_r28= 0;	cy_i28= 0;
	cy_r29= 0;	cy_i29= 0;

	/* If an LL test, init the subtract-2: */
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE && TEST_TYPE == TEST_TYPE_PRIMALITY)
	{
		cy_r00 = -2;
	}

	iroot = 0;	/* init sincos array index	*/
	root_incr = 1;	/* init sincos array index increment (set = 1 for normal carry pass, = 0 for wrapper pass)	*/
	scale = n2inv;	/* init inverse-weight scale factor  (set = 2/n for normal carry pass, = 1 for wrapper pass)	*/

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
	bjmodn22= bjmodn21+bjmodnini-n; bjmodn22= bjmodn22+ ( (-(int)((uint32)bjmodn22>> 31)) & n);
	bjmodn23= bjmodn22+bjmodnini-n; bjmodn23= bjmodn23+ ( (-(int)((uint32)bjmodn23>> 31)) & n);
	bjmodn24= bjmodn23+bjmodnini-n; bjmodn24= bjmodn24+ ( (-(int)((uint32)bjmodn24>> 31)) & n);
	bjmodn25= bjmodn24+bjmodnini-n; bjmodn25= bjmodn25+ ( (-(int)((uint32)bjmodn25>> 31)) & n);
	bjmodn26= bjmodn25+bjmodnini-n; bjmodn26= bjmodn26+ ( (-(int)((uint32)bjmodn26>> 31)) & n);
	bjmodn27= bjmodn26+bjmodnini-n; bjmodn27= bjmodn27+ ( (-(int)((uint32)bjmodn27>> 31)) & n);
	bjmodn28= bjmodn27+bjmodnini-n; bjmodn28= bjmodn28+ ( (-(int)((uint32)bjmodn28>> 31)) & n);
	bjmodn29= bjmodn28+bjmodnini-n; bjmodn29= bjmodn29+ ( (-(int)((uint32)bjmodn29>> 31)) & n);

	/* For Fermat-mod, IBDWT access patterns repeat with period NWT = {odd part of radix0},
	so for even radix0's only really need that many bjmodn and ii's, but that would require
	specialized carry macros tha`t don't update ii and bjmodn - not worth the trouble.
	*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		col=0;
		co2=(n >> nwt_bits)-1+30;
		co3=co2-30;		/* At the start of each new j-loop, co3=co2-radix(1)	*/
	}
	else
	{
		/* indices into IBDWT weights arrays (mod NWT) is here: */
		ii00= 0;
		ii01= (SW_DIV_N*n30/2) % nwt;
		ii02= (ii01+ ii01) % nwt;
		ii03= (ii02+ ii01) % nwt;
		ii04= (ii03+ ii01) % nwt;
		ii05= (ii04+ ii01) % nwt;
		ii06= (ii05+ ii01) % nwt;
		ii07= (ii06+ ii01) % nwt;
		ii08= (ii07+ ii01) % nwt;
		ii09= (ii08+ ii01) % nwt;
		ii10= (ii09+ ii01) % nwt;
		ii11= (ii10+ ii01) % nwt;
		ii12= (ii11+ ii01) % nwt;
		ii13= (ii12+ ii01) % nwt;
		ii14= (ii13+ ii01) % nwt;
		ii15= (ii14+ ii01) % nwt;
		ii16= (ii15+ ii01) % nwt;
		ii17= (ii16+ ii01) % nwt;
		ii18= (ii17+ ii01) % nwt;
		ii19= (ii18+ ii01) % nwt;
		ii20= (ii19+ ii01) % nwt;
		ii21= (ii20+ ii01) % nwt;
		ii22= (ii21+ ii01) % nwt;
		ii23= (ii22+ ii01) % nwt;
		ii24= (ii23+ ii01) % nwt;
		ii25= (ii24+ ii01) % nwt;
		ii26= (ii25+ ii01) % nwt;
		ii27= (ii26+ ii01) % nwt;
		ii28= (ii27+ ii01) % nwt;
		ii29= (ii28+ ii01) % nwt;

		/* Start this value off at N in Fermat-mod case, so (bjmodn >= sw) check in
		fermat_carry_norm_errcheck (cf. carry.h) yields a bigword (i == 1) for j= 0:
		*/
		bjmodn00= n;
		bjmodn15= n;
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
	/*
	!...gather the needed data (30 64-bit complex, i.e. 60 64-bit reals) and do 15 radix-2 transforms...
	*/
			u01=a[j1      ];			u02=a[j2    ];	/* x0a,b */
			rt =a[j1  +p01];			it =a[j2+p01];
			u03=u01-rt;				u04=u02-it;
			u01=u01+rt;				u02=u02+it;

			u05=a[j1  +p03];			u06=a[j2+p03];	/* x1a,b */
			rt =a[j1  +p02];			it =a[j2+p02];
			u07=u05-rt;				u08=u06-it;
			u05=u05+rt;				u06=u06+it;

			u09=a[j1  +p04];			u10=a[j2+p04];	/* x2a,b */
			rt =a[j1  +p05];			it =a[j2+p05];
			u11=u09-rt;				u12=u10-it;
			u09=u09+rt;				u10=u10+it;

			u13=a[j1  +p23];			u14=a[j2+p23];	/* x3a,b */
			rt =a[j1  +p22];			it =a[j2+p22];
			u15=u13-rt;				u16=u14-it;
			u13=u13+rt;				u14=u14+it;

			u17=a[j1  +p19];			u18=a[j2+p19];	/* x4a,b */
			rt =a[j1  +p18];			it =a[j2+p18];
			u19=u17-rt;				u20=u18-it;
			u17=u17+rt;				u18=u18+it;

			u21=a[j1  +p20];			u22=a[j2+p20];	/* x5a,b */
			rt =a[j1  +p21];			it =a[j2+p21];
			u23=u21-rt;				u24=u22-it;
			u21=u21+rt;				u22=u22+it;

			u25=a[j1  +p11];			u26=a[j2+p11];	/* x6a,b */
			rt =a[j1  +p10];			it =a[j2+p10];
			u27=u25-rt;				u28=u26-it;
			u25=u25+rt;				u26=u26+it;

			u29=a[j1  +p07];			u30=a[j2+p07];	/* x7a,b */
			rt =a[j1  +p06];			it =a[j2+p06];
			u31=u29-rt;				u32=u30-it;
			u29=u29+rt;				u30=u30+it;

			u33=a[j1  +p08];			u34=a[j2+p08];	/* x8a,b */
			rt =a[j1  +p09];			it =a[j2+p09];
			u35=u33-rt;				u36=u34-it;
			u33=u33+rt;				u34=u34+it;

			u37=a[j1  +p27];			u38=a[j2+p27];	/* x9a,b */
			rt =a[j1  +p26];			it =a[j2+p26];
			u39=u37-rt;				u40=u38-it;
			u37=u37+rt;				u38=u38+it;

			u41=a[j1  +p28];			u42=a[j2+p28];	/* x10a,b */
			rt =a[j1  +p29];			it =a[j2+p29];
			u43=u41-rt;				u44=u42-it;
			u41=u41+rt;				u42=u42+it;

			u45=a[j1  +p24];			u46=a[j2+p24];	/* x11a,b */
			rt =a[j1  +p25];			it =a[j2+p25];
			u47=u45-rt;				u48=u46-it;
			u45=u45+rt;				u46=u46+it;

			u49=a[j1  +p15];			u50=a[j2+p15];	/* x12a,b */
			rt =a[j1  +p14];			it =a[j2+p14];
			u51=u49-rt;				u52=u50-it;
			u49=u49+rt;				u50=u50+it;

			u53=a[j1  +p16];			u54=a[j2+p16];	/* x13a,b */
			rt =a[j1  +p17];			it =a[j2+p17];
			u55=u53-rt;				u56=u54-it;
			u53=u53+rt;				u54=u54+it;

			u57=a[j1  +p12];			u58=a[j2+p12];	/* x14a,b */
			rt =a[j1  +p13];			it =a[j2+p13];
			u59=u57-rt;				u60=u58-it;
			u57=u57+rt;				u58=u58+it;

/*       ...and now do two radix-15 transforms.	*/

/*...a1p[0:28:2]r use u[1:57:4]; a1p[0:28:2]i use u[2:58:4]	*/

/*   Gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do five radix-3 transforms...	*/
/*...Block 1:	*/
			t00=u01;			t01=u02;
			t02=u09;			t03=u10;
			rt =u05;			it =u06;
			t04=t02-rt;			t05=t03-it;
			t02=t02+rt;			t03=t03+it;
			t00=t00+t02;		t01=t01+t03;
			t02=t00+c3m1*t02;		t03=t01+c3m1*t03;
			rt =s*t04;			it =s*t05;
			t04=t02-it;			t05=t03+rt;
			t02=t02+it;			t03=t03-rt;

/*...Block 2:	*/
			t06=u33;			t07=u34;
			t08=u29;			t09=u30;
			rt =u25;			it =u26;
			t10=t08-rt;			t11=t09-it;
			t08=t08+rt;			t09=t09+it;
			t06=t06+t08;		t07=t07+t09;
			t08=t06+c3m1*t08;		t09=t07+c3m1*t09;
			rt =s*t10;			it =s*t11;
			t10=t08-it;			t11=t09+rt;
			t08=t08+it;			t09=t09-rt;

/*...Block 3:	*/
			t12=u53;			t13=u54;
			t14=u49;			t15=u50;
			rt =u57;			it =u58;
			t16=t14 -rt;		t17=t15 -it;
			t14=t14 +rt;		t15=t15 +it;
			t12=t12+t14;		t13=t13+t15;
			t14=t12+c3m1*t14;		t15=t13+c3m1*t15;
			rt =s*t16;			it =s*t17;
			t16=t14-it;			t17=t15+rt;
			t14=t14+it;			t15=t15-rt;

/*...Block 4:	*/
			t18=u17;			t19=u18;
			t20=u13;			t21=u14;
			rt =u21;			it =u22;
			t22=t20 -rt;		t23=t21 -it;
			t20=t20 +rt;		t21=t21 +it;
			t18=t18+t20;		t19=t19+t21;
			t20=t18+c3m1*t20;		t21=t19+c3m1*t21;
			rt =s*t22;			it =s*t23;
			t22=t20-it;			t23=t21+rt;
			t20=t20+it;			t21=t21-rt;

/*...Block 5:	*/
			t24=u37;			t25=u38;
			t26=u45;			t27=u46;
			rt =u41;			it =u42;
			t28=t26 -rt;		t29=t27 -it;
			t26=t26 +rt;		t27=t27 +it;
			t24=t24+t26;		t25=t25+t27;
			t26=t24+c3m1*t26;		t27=t25+c3m1*t27;
			rt =s*t28;			it =s*t29;
			t28=t26-it;			t29=t27+rt;
			t26=t26+it;			t27=t27-rt;

/*...and now do three radix-5 transforms:	*/
/*...Block 1:	*/
			rt = t24;			it = t25;
			t24= t06-rt;		t25= t07-it;
			t06= t06+rt;		t07= t07+it;
			rt = t18;			it = t19;
			t18= t12-rt;		t19= t13-it;
			t12= t12+rt;		t13= t13+it;

			rt = t06+t12;		it = t07+t13;
			t00= t00+rt;		t01= t01+it;
			rt = t00+cn1*rt;		it = t01+cn1*it;
			t12= cn2*(t06-t12);		t13= cn2*(t07-t13);
			t06= rt+t12;		t07= it+t13;
			t12= rt-t12;		t13= it-t13;
			rt = ss3*(t18-t24);		it = ss3*(t19-t25);
			t18= rt-sn1*t18;		t19= it-sn1*t19;
			t24= rt+sn2*t24;		t25= it+sn2*t25;

			a1p00r=t00;		a1p00i=t01;
			a1p06r=t06-t19;		a1p06i=t07+t18;
			a1p12r=t12-t25;		a1p12i=t13+t24;
			a1p18r=t12+t25;		a1p18i=t13-t24;
			a1p24r=t06+t19;		a1p24i=t07-t18;

/*...Block 2:	*/
			rt = t26;			it = t27;
			t26= t08-rt;		t27= t09-it;
			t08= t08+rt;		t09= t09+it;
			rt = t20;			it = t21;
			t20= t14-rt;		t21= t15-it;
			t14= t14+rt;		t15= t15+it;

			rt = t08+t14;		it = t09+t15;
			t02= t02+rt;		t03= t03+it;
			rt = t02+cn1*rt;		it = t03+cn1*it;
			t14= cn2*(t08-t14);		t15= cn2*(t09-t15);
			t08= rt+t14;		t09= it+t15;
			t14= rt-t14;		t15= it-t15;
			rt = ss3*(t20-t26);		it = ss3*(t21-t27);
			t20= rt-sn1*t20;		t21= it-sn1*t21;
			t26= rt+sn2*t26;		t27= it+sn2*t27;

			a1p10r=t02;		a1p10i=t03;
			a1p16r=t08-t21;		a1p16i=t09+t20;
			a1p22r=t14-t27;		a1p22i=t15+t26;
			a1p28r=t14+t27;		a1p28i=t15-t26;
			a1p04r=t08+t21;		a1p04i=t09-t20;

/*...Block 3:	*/
			rt = t28;			it = t29;
			t28= t10-rt;		t29= t11-it;
			t10= t10+rt;		t11= t11+it;
			rt = t22;			it = t23;
			t22= t16-rt;		t23= t17-it;
			t16= t16+rt;		t17= t17+it;

			rt = t10+t16;		it = t11+t17;
			t04= t04+rt;		t05= t05+it;
			rt = t04+cn1*rt;		it = t05+cn1*it;
			t16= cn2*(t10-t16);		t17= cn2*(t11-t17);
			t10= rt+t16;		t11= it+t17;
			t16= rt-t16;		t17= it-t17;
			rt = ss3*(t22-t28);		it = ss3*(t23-t29);
			t22= rt-sn1*t22;		t23= it-sn1*t23;
			t28= rt+sn2*t28;		t29= it+sn2*t29;

			a1p20r=t04;		a1p20i=t05;
			a1p26r=t10-t23;		a1p26i=t11+t22;
			a1p02r=t16-t29;		a1p02i=t17+t28;
			a1p08r=t16+t29;		a1p08i=t17-t28;
			a1p14r=t10+t23;		a1p14i=t11-t22;

/*...a1p[1:29:2]r use u[3:59:4]; a1p[1:29:2ii use u[4:60:4]	*/

/*   Gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do five radix-3 transforms...	*/
/*...Block 1:	*/
			t00=u03;			t01=u04;
			t02=u11;			t03=u12;
			rt =u07;			it =u08;
			t04=t02-rt;			t05=t03-it;
			t02=t02+rt;			t03=t03+it;
			t00=t00+t02;		t01=t01+t03;
			t02=t00+c3m1*t02;		t03=t01+c3m1*t03;
			rt =s*t04;			it =s*t05;
			t04=t02-it;			t05=t03+rt;
			t02=t02+it;			t03=t03-rt;

/*...Block 2:	*/
			t06=u35;			t07=u36;
			t08=u31;			t09=u32;
			rt =u27;			it =u28;
			t10=t08-rt;			t11=t09-it;
			t08=t08+rt;			t09=t09+it;
			t06=t06+t08;		t07=t07+t09;
			t08=t06+c3m1*t08;		t09=t07+c3m1*t09;
			rt =s*t10;			it =s*t11;
			t10=t08-it;			t11=t09+rt;
			t08=t08+it;			t09=t09-rt;

/*...Block 3:	*/
			t12=u55;			t13=u56;
			t14=u51;			t15=u52;
			rt =u59;			it =u60;
			t16=t14 -rt;		t17=t15 -it;
			t14=t14 +rt;		t15=t15 +it;
			t12=t12+t14;		t13=t13+t15;
			t14=t12+c3m1*t14;		t15=t13+c3m1*t15;
			rt =s*t16;			it =s*t17;
			t16=t14-it;			t17=t15+rt;
			t14=t14+it;			t15=t15-rt;

/*...Block 4:	*/
			t18=u19;			t19=u20;
			t20=u15;			t21=u16;
			rt =u23;			it =u24;
			t22=t20 -rt;		t23=t21 -it;
			t20=t20 +rt;		t21=t21 +it;
			t18=t18+t20;		t19=t19+t21;
			t20=t18+c3m1*t20;		t21=t19+c3m1*t21;
			rt =s*t22;			it =s*t23;
			t22=t20-it;			t23=t21+rt;
			t20=t20+it;			t21=t21-rt;

/*...Block 5:	*/
			t24=u39;			t25=u40;
			t26=u47;			t27=u48;
			rt =u43;			it =u44;
			t28=t26 -rt;		t29=t27 -it;
			t26=t26 +rt;		t27=t27 +it;
			t24=t24+t26;		t25=t25+t27;
			t26=t24+c3m1*t26;		t27=t25+c3m1*t27;
			rt =s*t28;			it =s*t29;
			t28=t26-it;			t29=t27+rt;
			t26=t26+it;			t27=t27-rt;

/*...and now do three radix-5 transforms:	*/
/*...Block 1:	*/
			rt = t24;			it = t25;
			t24= t06-rt;		t25= t07-it;
			t06= t06+rt;		t07= t07+it;
			rt = t18;			it = t19;
			t18= t12-rt;		t19= t13-it;
			t12= t12+rt;		t13= t13+it;

			rt = t06+t12;		it = t07+t13;
			t00= t00+rt;		t01= t01+it;
			rt = t00+cn1*rt;		it = t01+cn1*it;
			t12= cn2*(t06-t12);		t13= cn2*(t07-t13);
			t06= rt+t12;		t07= it+t13;
			t12= rt-t12;		t13= it-t13;
			rt = ss3*(t18-t24);		it = ss3*(t19-t25);
			t18= rt-sn1*t18;		t19= it-sn1*t19;
			t24= rt+sn2*t24;		t25= it+sn2*t25;

			a1p15r=t00;		a1p15i=t01;
			a1p21r=t06-t19;		a1p21i=t07+t18;
			a1p27r=t12-t25;		a1p27i=t13+t24;
			a1p03r=t12+t25;		a1p03i=t13-t24;
			a1p09r=t06+t19;		a1p09i=t07-t18;

/*...Block 2:	*/
			rt = t26;			it = t27;
			t26= t08-rt;		t27= t09-it;
			t08= t08+rt;		t09= t09+it;
			rt = t20;			it = t21;
			t20= t14-rt;		t21= t15-it;
			t14= t14+rt;		t15= t15+it;

			rt = t08+t14;		it = t09+t15;
			t02= t02+rt;		t03= t03+it;
			rt = t02+cn1*rt;		it = t03+cn1*it;
			t14= cn2*(t08-t14);		t15= cn2*(t09-t15);
			t08= rt+t14;		t09= it+t15;
			t14= rt-t14;		t15= it-t15;
			rt = ss3*(t20-t26);		it = ss3*(t21-t27);
			t20= rt-sn1*t20;		t21= it-sn1*t21;
			t26= rt+sn2*t26;		t27= it+sn2*t27;

			a1p25r=t02;		a1p25i=t03;
			a1p01r=t08-t21;		a1p01i=t09+t20;
			a1p07r=t14-t27;		a1p07i=t15+t26;
			a1p13r=t14+t27;		a1p13i=t15-t26;
			a1p19r=t08+t21;		a1p19i=t09-t20;

/*...Block 3:	*/
			rt = t28;			it = t29;
			t28= t10-rt;		t29= t11-it;
			t10= t10+rt;		t11= t11+it;
			rt = t22;			it = t23;
			t22= t16-rt;		t23= t17-it;
			t16= t16+rt;		t17= t17+it;

			rt = t10+t16;		it = t11+t17;
			t04= t04+rt;		t05= t05+it;
			rt = t04+cn1*rt;		it = t05+cn1*it;
			t16= cn2*(t10-t16);		t17= cn2*(t11-t17);
			t10= rt+t16;		t11= it+t17;
			t16= rt-t16;		t17= it-t17;
			rt = ss3*(t22-t28);		it = ss3*(t23-t29);
			t22= rt-sn1*t22;		t23= it-sn1*t23;
			t28= rt+sn2*t28;		t29= it+sn2*t29;

			a1p05r=t04;		a1p05i=t05;
			a1p11r=t10-t23;		a1p11i=t11+t22;
			a1p17r=t16-t29;		a1p17i=t17+t28;
			a1p23r=t16+t29;		a1p23i=t17-t28;
			a1p29r=t10+t23;		a1p29i=t11-t22;

/*...Now do the carries. Since the outputs would
    normally be getting dispatched to 30 separate blocks of the A-array, we need 30 separate carries.	*/

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
			 cmplx_carry_norm_nocheck0(a1p00r,a1p00i,cy_r00,bjmodn00   );
			cmplx_carry_norm_nocheck(a1p01r,a1p01i,cy_r01,bjmodn01,1 );
			cmplx_carry_norm_nocheck(a1p02r,a1p02i,cy_r02,bjmodn02,2 );
			cmplx_carry_norm_nocheck(a1p03r,a1p03i,cy_r03,bjmodn03,3 );
			cmplx_carry_norm_nocheck(a1p04r,a1p04i,cy_r04,bjmodn04,4 );
			cmplx_carry_norm_nocheck(a1p05r,a1p05i,cy_r05,bjmodn05,5 );
			cmplx_carry_norm_nocheck(a1p06r,a1p06i,cy_r06,bjmodn06,6 );
			cmplx_carry_norm_nocheck(a1p07r,a1p07i,cy_r07,bjmodn07,7 );
			cmplx_carry_norm_nocheck(a1p08r,a1p08i,cy_r08,bjmodn08,8 );
			cmplx_carry_norm_nocheck(a1p09r,a1p09i,cy_r09,bjmodn09,9 );
			cmplx_carry_norm_nocheck(a1p10r,a1p10i,cy_r10,bjmodn10,10);
			cmplx_carry_norm_nocheck(a1p11r,a1p11i,cy_r11,bjmodn11,11);
			cmplx_carry_norm_nocheck(a1p12r,a1p12i,cy_r12,bjmodn12,12);
			cmplx_carry_norm_nocheck(a1p13r,a1p13i,cy_r13,bjmodn13,13);
			cmplx_carry_norm_nocheck(a1p14r,a1p14i,cy_r14,bjmodn14,14);
			cmplx_carry_norm_nocheck(a1p15r,a1p15i,cy_r15,bjmodn15,15);
			cmplx_carry_norm_nocheck(a1p16r,a1p16i,cy_r16,bjmodn16,16);
			cmplx_carry_norm_nocheck(a1p17r,a1p17i,cy_r17,bjmodn17,17);
			cmplx_carry_norm_nocheck(a1p18r,a1p18i,cy_r18,bjmodn18,18);
			cmplx_carry_norm_nocheck(a1p19r,a1p19i,cy_r19,bjmodn19,19);
			cmplx_carry_norm_nocheck(a1p20r,a1p20i,cy_r20,bjmodn20,20);
			cmplx_carry_norm_nocheck(a1p21r,a1p21i,cy_r21,bjmodn21,21);
			cmplx_carry_norm_nocheck(a1p22r,a1p22i,cy_r22,bjmodn22,22);
			cmplx_carry_norm_nocheck(a1p23r,a1p23i,cy_r23,bjmodn23,23);
			cmplx_carry_norm_nocheck(a1p24r,a1p24i,cy_r24,bjmodn24,24);
			cmplx_carry_norm_nocheck(a1p25r,a1p25i,cy_r25,bjmodn25,25);
			cmplx_carry_norm_nocheck(a1p26r,a1p26i,cy_r26,bjmodn26,26);
			cmplx_carry_norm_nocheck(a1p27r,a1p27i,cy_r27,bjmodn27,27);
			cmplx_carry_norm_nocheck(a1p28r,a1p28i,cy_r28,bjmodn28,28);
			cmplx_carry_norm_nocheck(a1p29r,a1p29i,cy_r29,bjmodn29,29);

			i =((uint32)(sw - bjmodn00) >> 31);	/* get ready for the next set...	*/
			co2=co3;	/* For all data but the first set in each j-block, co2=co3. Thus, after the first block of data is done
					 and only then: for all subsequent blocks it's superfluous), this assignment decrements co2 by radix(1).	*/
		}
		else
		{
			fermat_carry_norm_nocheck(a1p00r,a1p00i,cy_r00,cy_i00,ii00,bjmodn00,0 *n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p01r,a1p01i,cy_r01,cy_i01,ii01,bjmodn01,1 *n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p02r,a1p02i,cy_r02,cy_i02,ii02,bjmodn02,2 *n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p03r,a1p03i,cy_r03,cy_i03,ii03,bjmodn03,3 *n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p04r,a1p04i,cy_r04,cy_i04,ii04,bjmodn04,4 *n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p05r,a1p05i,cy_r05,cy_i05,ii05,bjmodn05,5 *n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p06r,a1p06i,cy_r06,cy_i06,ii06,bjmodn06,6 *n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p07r,a1p07i,cy_r07,cy_i07,ii07,bjmodn07,7 *n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p08r,a1p08i,cy_r08,cy_i08,ii08,bjmodn08,8 *n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p09r,a1p09i,cy_r09,cy_i09,ii09,bjmodn09,9 *n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p10r,a1p10i,cy_r10,cy_i10,ii10,bjmodn10,10*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p11r,a1p11i,cy_r11,cy_i11,ii11,bjmodn11,11*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p12r,a1p12i,cy_r12,cy_i12,ii12,bjmodn12,12*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p13r,a1p13i,cy_r13,cy_i13,ii13,bjmodn13,13*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p14r,a1p14i,cy_r14,cy_i14,ii14,bjmodn14,14*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p15r,a1p15i,cy_r15,cy_i15,ii15,bjmodn15,15*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p16r,a1p16i,cy_r16,cy_i16,ii16,bjmodn16,16*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p17r,a1p17i,cy_r17,cy_i17,ii17,bjmodn17,17*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p18r,a1p18i,cy_r18,cy_i18,ii18,bjmodn18,18*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p19r,a1p19i,cy_r19,cy_i19,ii19,bjmodn19,19*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p20r,a1p20i,cy_r20,cy_i20,ii20,bjmodn20,20*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p21r,a1p21i,cy_r21,cy_i21,ii21,bjmodn21,21*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p22r,a1p22i,cy_r22,cy_i22,ii22,bjmodn22,22*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p23r,a1p23i,cy_r23,cy_i23,ii23,bjmodn23,23*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p24r,a1p24i,cy_r24,cy_i24,ii24,bjmodn24,24*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p25r,a1p25i,cy_r25,cy_i25,ii25,bjmodn25,25*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p26r,a1p26i,cy_r26,cy_i26,ii26,bjmodn26,26*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p27r,a1p27i,cy_r27,cy_i27,ii27,bjmodn27,27*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p28r,a1p28i,cy_r28,cy_i28,ii28,bjmodn28,28*n30,NRTM1,NRT_BITS);
			fermat_carry_norm_nocheck(a1p29r,a1p29i,cy_r29,cy_i29,ii29,bjmodn29,29*n30,NRTM1,NRT_BITS);
		}

/*...The radix-30 DIF pass is here:	*/
#if PFETCH
add0 = &a[j1];
prefetch_p_doubles(add0);
#endif

/*...First radix-15 block uses a1p[0:28:2] as inputs:	*/
/*   Gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do three length-5 transforms...	*/
/*...Block 1:	*/
		t00=a1p00r;			t01=a1p00i;	/* permuted versions of p00,24,06,18,12 */
		t02=a1p06r;			t03=a1p06i;
		rt =a1p24r;			it =a1p24i;
		t06=t02-rt;			t07=t03-it;
		t02=t02+rt;			t03=t03+it;
		t04=a1p12r;			t05=a1p12i;
		rt =a1p18r;			it =a1p18i;
		t08=t04-rt;			t09=t05-it;
		t04=t04+rt;			t05=t05+it;
#if PFETCH
addr = add0+p01;
prefetch_p_doubles(addr);
#endif
		rt = t02+t04;			it = t03+t05;
		t00= t00+rt;			t01= t01+it;
		rt = t00+cn1*rt;		it = t01+cn1*it;
		t04= cn2*(t02-t04);		t05= cn2*(t03-t05);
		t02= rt+t04;			t03= it+t05;
		t04= rt-t04;			t05= it-t05;
		rt = ss3*(t06-t08);		it = ss3*(t07-t09);
#if PFETCH
addr = add0+p02;
prefetch_p_doubles(addr);
#endif
		t08= rt+sn1*t08;		t09= it+sn1*t09;
		t06= rt-sn2*t06;		t07= it-sn2*t07;
		rt = t08;			it = t09;
		t08= t02+it;			t09= t03-rt;	/*<==prefer these to be stored in t08,9 	*/
		t02= t02-it;			t03= t03+rt;
		rt = t06;			it = t07;
		t06= t04+it;			t07= t05-rt;	/*<==prefer these to be stored in t06,7 	*/
		t04= t04-it;			t05= t05+rt;
#if PFETCH
addr = add0+p03;
prefetch_p_doubles(addr);
#endif

/*...Block 2:	*/
		t10=a1p10r;			t11=a1p10i;	/* permuted versions of p20,14,26,08,02 */
		t12=a1p16r;			t13=a1p16i;
		rt =a1p04r;			it =a1p04i;
		t16=t12-rt;			t17=t13-it;
		t12=t12+rt;			t13=t13+it;
		t14=a1p22r;			t15=a1p22i;
		rt =a1p28r;			it =a1p28i;
		t18=t14-rt;			t19=t15-it;
		t14=t14+rt;			t15=t15+it;
#if PFETCH
addr = add0+p04;
prefetch_p_doubles(addr);
#endif
		rt = t12+t14;			it = t13+t15;
		t10= t10+rt;			t11= t11+it;
		rt = t10+cn1*rt;		it = t11+cn1*it;
		t14= cn2*(t12-t14);		t15= cn2*(t13-t15);
		t12= rt+t14;			t13= it+t15;
		t14= rt-t14;			t15= it-t15;
		rt = ss3*(t16-t18);		it = ss3*(t17-t19);
#if PFETCH
addr = add0+p05;
prefetch_p_doubles(addr);
#endif
		t18= rt+sn1*t18;		t19= it+sn1*t19;
		t16= rt-sn2*t16;		t17= it-sn2*t17;
		rt = t18;			it = t19;
		t18= t12+it;			t19= t13-rt;
		t12= t12-it;			t13= t13+rt;
		rt = t16;			it = t17;
		t16= t14+it;			t17= t15-rt;
		t14= t14-it;			t15= t15+rt;
#if PFETCH
addr = add0+p06;
prefetch_p_doubles(addr);
#endif

/*...Block 3:	*/
		t20=a1p20r;			t21=a1p20i;	/* permuted versions of p10,04,16,28,22 */
		t22=a1p26r;			t23=a1p26i;
		rt =a1p14r;			it =a1p14i;
		t26=t22-rt;			t27=t23-it;
		t22=t22+rt;			t23=t23+it;
		t24=a1p02r;			t25=a1p02i;
		rt =a1p08r;			it =a1p08i;
		t28=t24-rt;			t29=t25-it;
		t24=t24+rt;			t25=t25+it;
#if PFETCH
addr = add0+p07;
prefetch_p_doubles(addr);
#endif
		rt = t22+t24;			it = t23+t25;
		t20= t20+rt;			t21= t21+it;
		rt = t20+cn1*rt;		it = t21+cn1*it;
		t24= cn2*(t22-t24);		t25= cn2*(t23-t25);
		t22= rt+t24;			t23= it+t25;
		t24= rt-t24;			t25= it-t25;
		rt = ss3*(t26-t28);		it = ss3*(t27-t29);
#if PFETCH
addr = add0+p08;
prefetch_p_doubles(addr);
#endif
		t28= rt+sn1*t28;		t29= it+sn1*t29;
		t26= rt-sn2*t26;		t27= it-sn2*t27;
		rt = t28;			it = t29;
		t28= t22+it;			t29= t23-rt;
		t22= t22-it;			t23= t23+rt;
		rt = t26;			it = t27;
		t26= t24+it;			t27= t25-rt;
		t24= t24-it;			t25= t25+rt;
#if PFETCH
addr = add0+p09;
prefetch_p_doubles(addr);
#endif

/*...and now do five radix-3 transforms:	*/
/*...Block 1:	*/
		rt =t20;			it =t21;
		t20=t10-rt;			t21=t11-it;
		t10=t10+rt;			t11=t11+it;
		t00=t00+t10;			t01=t01+t11;
		a1p00r=t00;			a1p00i=t01;
		t10=t00+c3m1*t10;		t11=t01+c3m1*t11;
		rt =s*t20;			it =s*t21;
		a1p02r=t10-it;		a1p02i=t11+rt;
		a1p04r=t10+it;		a1p04i=t11-rt;
#if PFETCH
addr = add0+p10;
prefetch_p_doubles(addr);
#endif

/*...Block 2:	*/
		rt =t22;			it =t23;
		t22=t12-rt;			t23=t13-it;
		t12=t12+rt;			t13=t13+it;
		t02=t02+t12;			t03=t03+t13;
		a1p26r=t02;			a1p26i=t03;
		t12=t02+c3m1*t12;		t13=t03+c3m1*t13;
		rt =s*t22;			it =s*t23;
		a1p28r=t12-it;		a1p28i=t13+rt;
		a1p24r=t12+it;		a1p24i=t13-rt;
#if PFETCH
addr = add0+p11;
prefetch_p_doubles(addr);
#endif

/*...Block 3:	*/
		rt =t24;			it =t25;
		t24=t14-rt;			t25=t15-it;
		t14=t14+rt;			t15=t15+it;
		t04=t04+t14;			t05=t05+t15;
		a1p18r=t04;			a1p18i=t05;
		t14=t04+c3m1*t14;		t15=t05+c3m1*t15;
		rt =s*t24;			it =s*t25;
		a1p20r=t14-it;		a1p20i=t15+rt;
		a1p22r=t14+it;		a1p22i=t15-rt;
#if PFETCH
addr = add0+p12;
prefetch_p_doubles(addr);
#endif

/*...Block 4:	*/
		rt =t26;			it =t27;
		t26=t16-rt;			t27=t17-it;
		t16=t16+rt;			t17=t17+it;
		t06=t06+t16;			t07=t07+t17;
		a1p16r=t06;			a1p16i=t07;
		t16=t06+c3m1*t16;		t17=t07+c3m1*t17;
		rt =s*t26;			it =s*t27;
		a1p12r=t16-it;		a1p12i=t17+rt;
		a1p14r=t16+it;		a1p14i=t17-rt;
#if PFETCH
addr = add0+p13;
prefetch_p_doubles(addr);
#endif

/*...Block 5:	*/
		rt =t28;			it =t29;
		t28=t18-rt;			t29=t19-it;
		t18=t18+rt;			t19=t19+it;
		t08=t08+t18;			t09=t09+t19;
		a1p08r=t08;			a1p08i=t09;
		t18=t08+c3m1*t18;		t19=t09+c3m1*t19;
		rt =s*t28;			it =s*t29;
		a1p10r=t18-it;		a1p10i=t19+rt;
		a1p06r=t18+it;		a1p06i=t19-rt;
#if PFETCH
addr = add0+p14;
prefetch_p_doubles(addr);
#endif

/*...Second radix-15 block uses a1p[1:29:2] as inputs:	*/
/*   Gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do three length-5 transforms...	*/
/*...Block 1:	*/
		t00=a1p15r;			t01=a1p15i;	/* permuted versions of p01,25,07,19,13 */
		t02=a1p21r;			t03=a1p21i;
		rt =a1p09r;			it =a1p09i;
		t06=t02-rt;			t07=t03-it;
		t02=t02+rt;			t03=t03+it;
		t04=a1p27r;			t05=a1p27i;
		rt =a1p03r;			it =a1p03i;
		t08=t04-rt;			t09=t05-it;
		t04=t04+rt;			t05=t05+it;
#if PFETCH
addr = add0+p15;
prefetch_p_doubles(addr);
#endif
		rt = t02+t04;			it = t03+t05;
		t00= t00+rt;			t01= t01+it;
		rt = t00+cn1*rt;		it = t01+cn1*it;
		t04= cn2*(t02-t04);		t05= cn2*(t03-t05);
		t02= rt+t04;			t03= it+t05;
		t04= rt-t04;			t05= it-t05;
		rt = ss3*(t06-t08);		it = ss3*(t07-t09);
#if PFETCH
addr = add0+p16;
prefetch_p_doubles(addr);
#endif
		t08= rt+sn1*t08;		t09= it+sn1*t09;
		t06= rt-sn2*t06;		t07= it-sn2*t07;
		rt = t08;			it = t09;
		t08= t02+it;			t09= t03-rt;	/*<==prefer these to be stored in t08,9 	*/
		t02= t02-it;			t03= t03+rt;
		rt = t06;			it = t07;
		t06= t04+it;			t07= t05-rt;	/*<==prefer these to be stored in t06,7 	*/
		t04= t04-it;			t05= t05+rt;
#if PFETCH
addr = add0+p17;
prefetch_p_doubles(addr);
#endif

/*...Block 2:	*/
		t10=a1p25r;			t11=a1p25i;	/* permuted versions of p21,15,27,09,03 */
		t12=a1p01r;			t13=a1p01i;
		rt =a1p19r;			it =a1p19i;
		t16=t12-rt;			t17=t13-it;
		t12=t12+rt;			t13=t13+it;
		t14=a1p07r;			t15=a1p07i;
		rt =a1p13r;			it =a1p13i;
		t18=t14-rt;			t19=t15-it;
		t14=t14+rt;			t15=t15+it;
#if PFETCH
addr = add0+p18;
prefetch_p_doubles(addr);
#endif
		rt = t12+t14;			it = t13+t15;
		t10= t10+rt;			t11= t11+it;
		rt = t10+cn1*rt;		it = t11+cn1*it;
		t14= cn2*(t12-t14);		t15= cn2*(t13-t15);
		t12= rt+t14;			t13= it+t15;
		t14= rt-t14;			t15= it-t15;
		rt = ss3*(t16-t18);		it = ss3*(t17-t19);
#if PFETCH
addr = add0+p19;
prefetch_p_doubles(addr);
#endif
		t18= rt+sn1*t18;		t19= it+sn1*t19;
		t16= rt-sn2*t16;		t17= it-sn2*t17;
		rt = t18;			it = t19;
		t18= t12+it;			t19= t13-rt;
		t12= t12-it;			t13= t13+rt;
		rt = t16;			it = t17;
		t16= t14+it;			t17= t15-rt;
		t14= t14-it;			t15= t15+rt;
#if PFETCH
addr = add0+p20;
prefetch_p_doubles(addr);
#endif

/*...Block 3:	*/
		t20=a1p05r;			t21=a1p05i;	/* permuted versions of p11,05,17,29,23 */
		t22=a1p11r;			t23=a1p11i;
		rt =a1p29r;			it =a1p29i;
		t26=t22-rt;			t27=t23-it;
		t22=t22+rt;			t23=t23+it;
		t24=a1p17r;			t25=a1p17i;
		rt =a1p23r;			it =a1p23i;
		t28=t24-rt;			t29=t25-it;
		t24=t24+rt;			t25=t25+it;
#if PFETCH
addr = add0+p21;
prefetch_p_doubles(addr);
#endif
		rt = t22+t24;			it = t23+t25;
		t20= t20+rt;			t21= t21+it;
		rt = t20+cn1*rt;		it = t21+cn1*it;
		t24= cn2*(t22-t24);		t25= cn2*(t23-t25);
		t22= rt+t24;			t23= it+t25;
		t24= rt-t24;			t25= it-t25;
		rt = ss3*(t26-t28);		it = ss3*(t27-t29);
#if PFETCH
addr = add0+p22;
prefetch_p_doubles(addr);
#endif
		t28= rt+sn1*t28;		t29= it+sn1*t29;
		t26= rt-sn2*t26;		t27= it-sn2*t27;
		rt = t28;			it = t29;
		t28= t22+it;			t29= t23-rt;
		t22= t22-it;			t23= t23+rt;
		rt = t26;			it = t27;
		t26= t24+it;			t27= t25-rt;
		t24= t24-it;			t25= t25+rt;
#if PFETCH
addr = add0+p23;
prefetch_p_doubles(addr);
#endif

/*...and now do five radix-3 transforms:	*/
/*...Block 1:	*/
		rt =t20;			it =t21;
		t20=t10-rt;			t21=t11-it;
		t10=t10+rt;			t11=t11+it;
		t00=t00+t10;			t01=t01+t11;
		a1p01r=t00;			a1p01i=t01;
		t10=t00+c3m1*t10;		t11=t01+c3m1*t11;
		rt =s*t20;			it =s*t21;
		a1p03r=t10-it;		a1p03i=t11+rt;
		a1p05r=t10+it;		a1p05i=t11-rt;
#if PFETCH
addr = add0+p24;
prefetch_p_doubles(addr);
#endif

/*...Block 2:	*/
		rt =t22;			it =t23;
		t22=t12-rt;			t23=t13-it;
		t12=t12+rt;			t13=t13+it;
		t02=t02+t12;			t03=t03+t13;
		a1p27r=t02;			a1p27i=t03;
		t12=t02+c3m1*t12;		t13=t03+c3m1*t13;
		rt =s*t22;			it =s*t23;
		a1p29r=t12-it;		a1p29i=t13+rt;
		a1p25r=t12+it;		a1p25i=t13-rt;
#if PFETCH
addr = add0+p25;
prefetch_p_doubles(addr);
#endif

/*...Block 3:	*/
		rt =t24;			it =t25;
		t24=t14-rt;			t25=t15-it;
		t14=t14+rt;			t15=t15+it;
		t04=t04+t14;			t05=t05+t15;
		a1p19r=t04;			a1p19i=t05;
		t14=t04+c3m1*t14;		t15=t05+c3m1*t15;
		rt =s*t24;			it =s*t25;
		a1p21r=t14-it;		a1p21i=t15+rt;
		a1p23r=t14+it;		a1p23i=t15-rt;
#if PFETCH
addr = add0+p26;
prefetch_p_doubles(addr);
#endif

/*...Block 4:	*/
		rt =t26;			it =t27;
		t26=t16-rt;			t27=t17-it;
		t16=t16+rt;			t17=t17+it;
		t06=t06+t16;			t07=t07+t17;
		a1p17r=t06;			a1p17i=t07;
		t16=t06+c3m1*t16;		t17=t07+c3m1*t17;
		rt =s*t26;			it =s*t27;
		a1p13r=t16-it;		a1p13i=t17+rt;
		a1p15r=t16+it;		a1p15i=t17-rt;
#if PFETCH
addr = add0+p27;
prefetch_p_doubles(addr);
#endif

/*...Block 5:	*/
		rt =t28;			it =t29;
		t28=t18-rt;			t29=t19-it;
		t18=t18+rt;			t19=t19+it;
		t08=t08+t18;			t09=t09+t19;
		a1p09r=t08;			a1p09i=t09;
		t18=t08+c3m1*t18;		t19=t09+c3m1*t19;
		rt =s*t28;			it =s*t29;
		a1p11r=t18-it;		a1p11i=t19+rt;
		a1p07r=t18+it;		a1p07i=t19-rt;
#if PFETCH
addr = add0+p28;
prefetch_p_doubles(addr);
#endif

/*...and now do 15 radix-2 transforms:	*/

		a[j1    ]=a1p00r+a1p01r;	a[j2    ]=a1p00i+a1p01i;
		a[j1+p01]=a1p00r-a1p01r;	a[j2+p01]=a1p00i-a1p01i;

		a[j1+p04]=a1p02r+a1p03r;	a[j2+p04]=a1p02i+a1p03i;
		a[j1+p05]=a1p02r-a1p03r;	a[j2+p05]=a1p02i-a1p03i;

		a[j1+p03]=a1p04r+a1p05r;	a[j2+p03]=a1p04i+a1p05i;
		a[j1+p02]=a1p04r-a1p05r;	a[j2+p02]=a1p04i-a1p05i;

		a[j1+p28]=a1p06r+a1p07r;	a[j2+p28]=a1p06i+a1p07i;
		a[j1+p29]=a1p06r-a1p07r;	a[j2+p29]=a1p06i-a1p07i;

		a[j1+p27]=a1p08r+a1p09r;	a[j2+p27]=a1p08i+a1p09i;
		a[j1+p26]=a1p08r-a1p09r;	a[j2+p26]=a1p08i-a1p09i;

		a[j1+p24]=a1p10r+a1p11r;	a[j2+p24]=a1p10i+a1p11i;
		a[j1+p25]=a1p10r-a1p11r;	a[j2+p25]=a1p10i-a1p11i;

		a[j1+p23]=a1p12r+a1p13r;	a[j2+p23]=a1p12i+a1p13i;
		a[j1+p22]=a1p12r-a1p13r;	a[j2+p22]=a1p12i-a1p13i;
#if PFETCH
addr = add0+p29;
prefetch_p_doubles(addr);
#endif
		a[j1+p20]=a1p14r+a1p15r;	a[j2+p20]=a1p14i+a1p15i;
		a[j1+p21]=a1p14r-a1p15r;	a[j2+p21]=a1p14i-a1p15i;

		a[j1+p19]=a1p16r+a1p17r;	a[j2+p19]=a1p16i+a1p17i;
		a[j1+p18]=a1p16r-a1p17r;	a[j2+p18]=a1p16i-a1p17i;

		a[j1+p16]=a1p18r+a1p19r;	a[j2+p16]=a1p18i+a1p19i;
		a[j1+p17]=a1p18r-a1p19r;	a[j2+p17]=a1p18i-a1p19i;

		a[j1+p15]=a1p20r+a1p21r;	a[j2+p15]=a1p20i+a1p21i;
		a[j1+p14]=a1p20r-a1p21r;	a[j2+p14]=a1p20i-a1p21i;

		a[j1+p12]=a1p22r+a1p23r;	a[j2+p12]=a1p22i+a1p23i;
		a[j1+p13]=a1p22r-a1p23r;	a[j2+p13]=a1p22i-a1p23i;

		a[j1+p11]=a1p24r+a1p25r;	a[j2+p11]=a1p24i+a1p25i;
		a[j1+p10]=a1p24r-a1p25r;	a[j2+p10]=a1p24i-a1p25i;

		a[j1+p08]=a1p26r+a1p27r;	a[j2+p08]=a1p26i+a1p27i;
		a[j1+p09]=a1p26r-a1p27r;	a[j2+p09]=a1p26i-a1p27i;

		a[j1+p07]=a1p28r+a1p29r;	a[j2+p07]=a1p28i+a1p29i;
		a[j1+p06]=a1p28r-a1p29r;	a[j2+p06]=a1p28i-a1p29i;

		iroot += root_incr;		/* increment sincos index.	*/
		}

		if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
		{
			jstart += nwt;
			jhi    += nwt;

			col += 30;
			co3 -= 30;
		}
	}

	if(root_incr==0)break;

/*   Wraparound carry cleanup loop is here: ***
!
!   (1) Invert the radix-30 forward DIF FFT of the first block of 30 complex elements in A and unweight;
!   (2) Propagate cleanup carries among the real and imaginary parts of the 30 outputs of (1);
!   (3) Reweight and perform a radix-30 forward DIF FFT on the result of (2);
!   (4) If any of the exit carries from (2) are nonzero, advance to the next 30 elements and repeat (1-4).
*/
	if(MODULUS_TYPE == MODULUS_TYPE_MERSENNE)
	{
		t00 = cy_r29;
		cy_r29= cy_r28;
		cy_r28= cy_r27;
		cy_r27= cy_r26;
		cy_r26= cy_r25;
		cy_r25= cy_r24;
		cy_r24= cy_r23;
		cy_r23= cy_r22;
		cy_r22= cy_r21;
		cy_r21= cy_r20;
		cy_r20= cy_r19;
		cy_r19= cy_r18;
		cy_r18= cy_r17;
		cy_r17= cy_r16;
		cy_r16= cy_r15;
		cy_r15= cy_r14;
		cy_r14= cy_r13;
		cy_r13= cy_r12;
		cy_r12= cy_r11;
		cy_r11= cy_r10;
		cy_r10= cy_r09;
		cy_r09= cy_r08;
		cy_r08= cy_r07;
		cy_r07= cy_r06;
		cy_r06= cy_r05;
		cy_r05= cy_r04;
		cy_r04= cy_r03;
		cy_r03= cy_r02;
		cy_r02= cy_r01;
		cy_r01= cy_r00;
		cy_r00= t00 ;
	}
	else
	{
		/* ...The 2 Mo"bius carries are here: */
		t00   = cy_r29;	t01   = cy_i29;
		cy_r29= cy_r28;	cy_i29= cy_i28;
		cy_r28= cy_r27;	cy_i28= cy_i27;
		cy_r27= cy_r26;	cy_i27= cy_i26;
		cy_r26= cy_r25;	cy_i26= cy_i25;
		cy_r25= cy_r24;	cy_i25= cy_i24;
		cy_r24= cy_r23;	cy_i24= cy_i23;
		cy_r23= cy_r22;	cy_i23= cy_i22;
		cy_r22= cy_r21;	cy_i22= cy_i21;
		cy_r21= cy_r20;	cy_i21= cy_i20;
		cy_r20= cy_r19;	cy_i20= cy_i19;
		cy_r19= cy_r18;	cy_i19= cy_i18;
		cy_r18= cy_r17;	cy_i18= cy_i17;
		cy_r17= cy_r16;	cy_i17= cy_i16;
		cy_r16= cy_r15;	cy_i16= cy_i15;
		cy_r15= cy_r14;	cy_i15= cy_i14;
		cy_r14= cy_r13;	cy_i14= cy_i13;
		cy_r13= cy_r12;	cy_i13= cy_i12;
		cy_r12= cy_r11;	cy_i12= cy_i11;
		cy_r11= cy_r10;	cy_i11= cy_i10;
		cy_r10= cy_r09;	cy_i10= cy_i09;
		cy_r09= cy_r08;	cy_i09= cy_i08;
		cy_r08= cy_r07;	cy_i08= cy_i07;
		cy_r07= cy_r06;	cy_i07= cy_i06;
		cy_r06= cy_r05;	cy_i06= cy_i05;
		cy_r05= cy_r04;	cy_i05= cy_i04;
		cy_r04= cy_r03;	cy_i04= cy_i03;
		cy_r03= cy_r02;	cy_i03= cy_i02;
		cy_r02= cy_r01;	cy_i02= cy_i01;
		cy_r01= cy_r00;	cy_i01= cy_i00;
		cy_r00= -t01 ;	cy_i00= +t00 ;
	}

	iroot = 0;
	root_incr = 0;
	scale = 1;

	jstart = 0;
	/*
	For right-angle transform need *complex* elements for wraparound, so jhi needs to be twice as large
	*/
	if(TRANSFORM_TYPE == RIGHT_ANGLE)
	{
		jhi =15;
	}
	else
	{
		jhi = 7;
	}
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
		a[j+p22] *= radix_inv;
		a[j+p23] *= radix_inv;
		a[j+p24] *= radix_inv;
		a[j+p25] *= radix_inv;
		a[j+p26] *= radix_inv;
		a[j+p27] *= radix_inv;
		a[j+p28] *= radix_inv;
		a[j+p29] *= radix_inv;
	}
}

	if(fabs(cy_r00)+fabs(cy_r01)+fabs(cy_r02)+fabs(cy_r03)+fabs(cy_r04)+fabs(cy_r05)+fabs(cy_r06)+fabs(cy_r07)+fabs(cy_r08)+fabs(cy_r09)+fabs(cy_r10)+fabs(cy_r11)+fabs(cy_r12)+fabs(cy_r13)+fabs(cy_r14)+fabs(cy_r15)+fabs(cy_r16)+fabs(cy_r17)+fabs(cy_r18)+fabs(cy_r19)+fabs(cy_r20)+fabs(cy_r21)+fabs(cy_r22)+fabs(cy_r23)+fabs(cy_r24)+fabs(cy_r25)+fabs(cy_r26)+fabs(cy_r27)+fabs(cy_r28)+fabs(cy_r29)
		+fabs(cy_i00)+fabs(cy_i01)+fabs(cy_i02)+fabs(cy_i03)+fabs(cy_i04)+fabs(cy_i05)+fabs(cy_i06)+fabs(cy_i07)+fabs(cy_i08)+fabs(cy_i09)+fabs(cy_i10)+fabs(cy_i11)+fabs(cy_i12)+fabs(cy_i13)+fabs(cy_i14)+fabs(cy_i15)+fabs(cy_i16)+fabs(cy_i17)+fabs(cy_i18)+fabs(cy_i19)+fabs(cy_i20)+fabs(cy_i21)+fabs(cy_i22)+fabs(cy_i23)+fabs(cy_i24)+fabs(cy_i25)+fabs(cy_i26)+fabs(cy_i27)+fabs(cy_i28)+fabs(cy_i29) != 0.0)
	{
			sprintf(cbuf,"FATAL: iter = %10d; nonzero exit carry in radix30_ditN_cy_dif1 - input wordsize may be too small.\n",iter);
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

void radix30_dif_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-30 complex DIF FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix15_dif_pass for details on the radix-15 subtransforms.
*/
	int j,j1,j2;
	static int n30,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29, first_entry=TRUE;
	static double c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/

					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double
	 t00r,t01r,t02r,t03r,t04r,t05r,t06r,t07r,t08r,t09r,t10r,t11r,t12r,t13r,t14r,t15r,t16r,t17r,t18r,t19r,t20r,t21r,t22r,t23r,t24r,t25r,t26r,t27r,t28r,t29r
	,t00i,t01i,t02i,t03i,t04i,t05i,t06i,t07i,t08i,t09i,t10i,t11i,t12i,t13i,t14i,t15i,t16i,t17i,t18i,t19i,t20i,t21i,t22i,t23i,t24i,t25i,t26i,t27i,t28i,t29i;

	if(!first_entry && (n/30) != n30)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n30=n/30;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n30;
		p02 = p01 + n30;		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p03 = p02 + n30;		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p04 = p03 + n30;		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p05 = p04 + n30;		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p06 = p05 + n30;		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p07 = p06 + n30;		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p08 = p07 + n30;		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p09 = p08 + n30;		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p10 = p09 + n30;		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p11 = p10 + n30;		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p12 = p11 + n30;		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p13 = p12 + n30;		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p14 = p13 + n30;		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p15 = p14 + n30;		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p16 = p15 + n30;		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
		p17 = p16 + n30;		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p18 = p17 + n30;		p17 = p17 + ( (p17 >> DAT_BITS) << PAD_BITS );
		p19 = p18 + n30;		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 = p19 + n30;		p19 = p19 + ( (p19 >> DAT_BITS) << PAD_BITS );
		p21 = p20 + n30;		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p22 = p21 + n30;		p21 = p21 + ( (p21 >> DAT_BITS) << PAD_BITS );
		p23 = p22 + n30;		p22 = p22 + ( (p22 >> DAT_BITS) << PAD_BITS );
		p24 = p23 + n30;		p23 = p23 + ( (p23 >> DAT_BITS) << PAD_BITS );
		p25 = p24 + n30;		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p26 = p25 + n30;		p25 = p25 + ( (p25 >> DAT_BITS) << PAD_BITS );
		p27 = p26 + n30;		p26 = p26 + ( (p26 >> DAT_BITS) << PAD_BITS );
		p28 = p27 + n30;		p27 = p27 + ( (p27 >> DAT_BITS) << PAD_BITS );
		p29 = p28 + n30;		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
								p29 = p29 + ( (p29 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-30 pass is here.	*/

	for(j=0; j < n30; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

/*...gather the needed data (30 64-bit complex, i.e. 60 64-bit reals)...*/

	/*
	Twiddleless version requires us to swap inputs as follows:
	indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29
		  => 0,15,28,13,26,11,24, 9,22, 7,20, 5,18, 3,16, 1,14,29,12,27,10,25, 8,23, 6,21, 4,19, 2,17 modulo 30.
	I.e. start out with first 15-vector of indices {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28},
	permute those according to {0, 2, 4, 6, 8,10,12,14,16,18,20,22,24,26,28}*29%30 => {0,28,26,24,22,20,18,16,14,12,10, 8, 6, 4, 2},
	then each is head of a length-2 list of indices with decrement 15:

		EVEN	[00,28,26,24,22,20,18,16,14,12,10,08,06,04,02]
		ODD		[15,13,11,09,07,05,03,01,29,27,25,23,21,19,17]

	The output permutations are inlined below.
	*/
	#if 1
	/*...gather the needed data (30 64-bit complex, i.e. 60 64-bit reals) and do 2 radix-15 DIF transforms...*/
		RADIX_15_DIF(
			a[j1    ],a[j2    ],a[j1+p28],a[j2+p28],a[j1+p26],a[j2+p26],a[j1+p24],a[j2+p24],a[j1+p22],a[j2+p22],a[j1+p20],a[j2+p20],a[j1+p18],a[j2+p18],a[j1+p16],a[j2+p16],a[j1+p14],a[j2+p14],a[j1+p12],a[j2+p12],a[j1+p10],a[j2+p10],a[j1+p08],a[j2+p08],a[j1+p06],a[j2+p06],a[j1+p04],a[j2+p04],a[j1+p02],a[j2+p02],
			t00r,t00i,t02r,t02i,t04r,t04i,t26r,t26i,t28r,t28i,t24r,t24i,t18r,t18i,t20r,t20i,t22r,t22i,t16r,t16i,t12r,t12i,t14r,t14i,t08r,t08i,t10r,t10i,t06r,t06i
		);
		RADIX_15_DIF(
			a[j1+p15],a[j2+p15],a[j1+p13],a[j2+p13],a[j1+p11],a[j2+p11],a[j1+p09],a[j2+p09],a[j1+p07],a[j2+p07],a[j1+p05],a[j2+p05],a[j1+p03],a[j2+p03],a[j1+p01],a[j2+p01],a[j1+p29],a[j2+p29],a[j1+p27],a[j2+p27],a[j1+p25],a[j2+p25],a[j1+p23],a[j2+p23],a[j1+p21],a[j2+p21],a[j1+p19],a[j2+p19],a[j1+p17],a[j2+p17],
			t01r,t01i,t03r,t03i,t05r,t05i,t27r,t27i,t29r,t29i,t25r,t25i,t19r,t19i,t21r,t21i,t23r,t23i,t17r,t17i,t13r,t13i,t15r,t15i,t09r,t09i,t11r,t11i,t07r,t07i
		);
	/*...and now do 15 radix-2 transforms:	*/
		a[j1    ]=t00r+t01r;	a[j2    ]=t00i+t01i;
		a[j1+p01]=t00r-t01r;	a[j2+p01]=t00i-t01i;
		a[j1+p04]=t02r+t03r;	a[j2+p04]=t02i+t03i;
		a[j1+p05]=t02r-t03r;	a[j2+p05]=t02i-t03i;
		a[j1+p03]=t04r+t05r;	a[j2+p03]=t04i+t05i;
		a[j1+p02]=t04r-t05r;	a[j2+p02]=t04i-t05i;
		a[j1+p07]=t06r+t07r;	a[j2+p07]=t06i+t07i;
		a[j1+p06]=t06r-t07r;	a[j2+p06]=t06i-t07i;
		a[j1+p11]=t08r+t09r;	a[j2+p11]=t08i+t09i;
		a[j1+p10]=t08r-t09r;	a[j2+p10]=t08i-t09i;
		a[j1+p08]=t10r+t11r;	a[j2+p08]=t10i+t11i;
		a[j1+p09]=t10r-t11r;	a[j2+p09]=t10i-t11i;
		a[j1+p15]=t12r+t13r;	a[j2+p15]=t12i+t13i;
		a[j1+p14]=t12r-t13r;	a[j2+p14]=t12i-t13i;
		a[j1+p12]=t14r+t15r;	a[j2+p12]=t14i+t15i;
		a[j1+p13]=t14r-t15r;	a[j2+p13]=t14i-t15i;
		a[j1+p16]=t16r+t17r;	a[j2+p16]=t16i+t17i;
		a[j1+p17]=t16r-t17r;	a[j2+p17]=t16i-t17i;
		a[j1+p23]=t18r+t19r;	a[j2+p23]=t18i+t19i;
		a[j1+p22]=t18r-t19r;	a[j2+p22]=t18i-t19i;
		a[j1+p20]=t20r+t21r;	a[j2+p20]=t20i+t21i;
		a[j1+p21]=t20r-t21r;	a[j2+p21]=t20i-t21i;
		a[j1+p19]=t22r+t23r;	a[j2+p19]=t22i+t23i;
		a[j1+p18]=t22r-t23r;	a[j2+p18]=t22i-t23i;
		a[j1+p24]=t24r+t25r;	a[j2+p24]=t24i+t25i;
		a[j1+p25]=t24r-t25r;	a[j2+p25]=t24i-t25i;
		a[j1+p28]=t26r+t27r;	a[j2+p28]=t26i+t27i;
		a[j1+p29]=t26r-t27r;	a[j2+p29]=t26i-t27i;
		a[j1+p27]=t28r+t29r;	a[j2+p27]=t28i+t29i;
		a[j1+p26]=t28r-t29r;	a[j2+p26]=t28i-t29i;

	#else

/*...First radix-15 block uses a1p[0:28:2] as inputs:	*/
/*   Gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do three length-5 transforms...	*/
/*...Block 1:	*/
		t00=a[j1    ];		t01=a[j2    ];	/* permuted versions of p00,24,06,18,12 */
		t02=a[j1+p06];		t03=a[j2+p06];
		rt =a[j1+p24];		it =a[j2+p24];
		t06=t02-rt;			t07=t03-it;
		t02=t02+rt;			t03=t03+it;
		t04=a[j1+p12];		t05=a[j2+p12];
		rt =a[j1+p18];		it =a[j2+p18];
		t08=t04-rt;			t09=t05-it;
		t04=t04+rt;			t05=t05+it;

		rt = t02+t04;			it = t03+t05;
		t00= t00+rt;			t01= t01+it;
		rt = t00+cn1*rt;		it = t01+cn1*it;
		t04= cn2*(t02-t04);		t05= cn2*(t03-t05);
		t02= rt+t04;			t03= it+t05;
		t04= rt-t04;			t05= it-t05;
		rt = ss3*(t06-t08);		it = ss3*(t07-t09);
		t08= rt+sn1*t08;		t09= it+sn1*t09;
		t06= rt-sn2*t06;		t07= it-sn2*t07;
		rt = t08;			it = t09;
		t08= t02+it;			t09= t03-rt;	/*<==prefer these to be stored in t08,9 	*/
		t02= t02-it;			t03= t03+rt;
		rt = t06;			it = t07;
		t06= t04+it;			t07= t05-rt;	/*<==prefer these to be stored in t06,7 	*/
		t04= t04-it;			t05= t05+rt;

/*...Block 2:	*/
		t10=a[j1+p10];		t11=a[j2+p10];	/* permuted versions of p20,14,26,08,02 */
		t12=a[j1+p16];		t13=a[j2+p16];
		rt =a[j1+p04];		it =a[j2+p04];
		t16=t12-rt;			t17=t13-it;
		t12=t12+rt;			t13=t13+it;
		t14=a[j1+p22];		t15=a[j2+p22];
		rt =a[j1+p28];		it =a[j2+p28];
		t18=t14-rt;			t19=t15-it;
		t14=t14+rt;			t15=t15+it;

		rt = t12+t14;			it = t13+t15;
		t10= t10+rt;			t11= t11+it;
		rt = t10+cn1*rt;		it = t11+cn1*it;
		t14= cn2*(t12-t14);		t15= cn2*(t13-t15);
		t12= rt+t14;			t13= it+t15;
		t14= rt-t14;			t15= it-t15;
		rt = ss3*(t16-t18);		it = ss3*(t17-t19);
		t18= rt+sn1*t18;		t19= it+sn1*t19;
		t16= rt-sn2*t16;		t17= it-sn2*t17;
		rt = t18;			it = t19;
		t18= t12+it;			t19= t13-rt;
		t12= t12-it;			t13= t13+rt;
		rt = t16;			it = t17;
		t16= t14+it;			t17= t15-rt;
		t14= t14-it;			t15= t15+rt;

/*...Block 3:	*/
		t20=a[j1+p20];		t21=a[j2+p20];	/* permuted versions of p10,04,16,28,22 */
		t22=a[j1+p26];		t23=a[j2+p26];
		rt =a[j1+p14];		it =a[j2+p14];
		t26=t22-rt;			t27=t23-it;
		t22=t22+rt;			t23=t23+it;
		t24=a[j1+p02];		t25=a[j2+p02];
		rt =a[j1+p08];		it =a[j2+p08];
		t28=t24-rt;			t29=t25-it;
		t24=t24+rt;			t25=t25+it;

		rt = t22+t24;			it = t23+t25;
		t20= t20+rt;			t21= t21+it;
		rt = t20+cn1*rt;		it = t21+cn1*it;
		t24= cn2*(t22-t24);		t25= cn2*(t23-t25);
		t22= rt+t24;			t23= it+t25;
		t24= rt-t24;			t25= it-t25;
		rt = ss3*(t26-t28);		it = ss3*(t27-t29);
		t28= rt+sn1*t28;		t29= it+sn1*t29;
		t26= rt-sn2*t26;		t27= it-sn2*t27;
		rt = t28;			it = t29;
		t28= t22+it;			t29= t23-rt;
		t22= t22-it;			t23= t23+rt;
		rt = t26;			it = t27;
		t26= t24+it;			t27= t25-rt;
		t24= t24-it;			t25= t25+rt;

/*...and now do five radix-3 transforms:	*/
/*...Block 1:	*/
		rt =t20;			it =t21;
		t20=t10-rt;			t21=t11-it;
		t10=t10+rt;			t11=t11+it;
		t00=t00+t10;			t01=t01+t11;
		a1p00r=t00;			a1p00i=t01;
		t10=t00+c3m1*t10;		t11=t01+c3m1*t11;
		rt =s*t20;			it =s*t21;
		a1p02r=t10-it;		a1p02i=t11+rt;
		a1p04r=t10+it;		a1p04i=t11-rt;

/*...Block 2:	*/
		rt =t22;			it =t23;
		t22=t12-rt;			t23=t13-it;
		t12=t12+rt;			t13=t13+it;
		t02=t02+t12;			t03=t03+t13;
		a1p26r=t02;			a1p26i=t03;
		t12=t02+c3m1*t12;		t13=t03+c3m1*t13;
		rt =s*t22;			it =s*t23;
		a1p28r=t12-it;		a1p28i=t13+rt;
		a1p24r=t12+it;		a1p24i=t13-rt;

/*...Block 3:	*/
		rt =t24;			it =t25;
		t24=t14-rt;			t25=t15-it;
		t14=t14+rt;			t15=t15+it;
		t04=t04+t14;			t05=t05+t15;
		a1p18r=t04;			a1p18i=t05;
		t14=t04+c3m1*t14;		t15=t05+c3m1*t15;
		rt =s*t24;			it =s*t25;
		a1p20r=t14-it;		a1p20i=t15+rt;
		a1p22r=t14+it;		a1p22i=t15-rt;

/*...Block 4:	*/
		rt =t26;			it =t27;
		t26=t16-rt;			t27=t17-it;
		t16=t16+rt;			t17=t17+it;
		t06=t06+t16;			t07=t07+t17;
		a1p16r=t06;			a1p16i=t07;
		t16=t06+c3m1*t16;		t17=t07+c3m1*t17;
		rt =s*t26;			it =s*t27;
		a1p12r=t16-it;		a1p12i=t17+rt;
		a1p14r=t16+it;		a1p14i=t17-rt;

/*...Block 5:	*/
		rt =t28;			it =t29;
		t28=t18-rt;			t29=t19-it;
		t18=t18+rt;			t19=t19+it;
		t08=t08+t18;			t09=t09+t19;
		a1p08r=t08;			a1p08i=t09;
		t18=t08+c3m1*t18;		t19=t09+c3m1*t19;
		rt =s*t28;			it =s*t29;
		a1p10r=t18-it;		a1p10i=t19+rt;
		a1p06r=t18+it;		a1p06i=t19-rt;

/*...Second radix-15 block uses a1p[1:29:2] as inputs:	*/
/*   Gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do three length-5 transforms...	*/
/*...Block 1:	*/
		t00=a[j1+p15];		t01=a[j2+p15];	/* permuted versions of p01,25,07,19,13 */
		t02=a[j1+p21];		t03=a[j2+p21];
		rt =a[j1+p09];		it =a[j2+p09];
		t06=t02-rt;			t07=t03-it;
		t02=t02+rt;			t03=t03+it;
		t04=a[j1+p27];		t05=a[j2+p27];
		rt =a[j1+p03];		it =a[j2+p03];
		t08=t04-rt;			t09=t05-it;
		t04=t04+rt;			t05=t05+it;

		rt = t02+t04;			it = t03+t05;
		t00= t00+rt;			t01= t01+it;
		rt = t00+cn1*rt;		it = t01+cn1*it;
		t04= cn2*(t02-t04);		t05= cn2*(t03-t05);
		t02= rt+t04;			t03= it+t05;
		t04= rt-t04;			t05= it-t05;
		rt = ss3*(t06-t08);		it = ss3*(t07-t09);
		t08= rt+sn1*t08;		t09= it+sn1*t09;
		t06= rt-sn2*t06;		t07= it-sn2*t07;
		rt = t08;			it = t09;
		t08= t02+it;			t09= t03-rt;	/*<==prefer these to be stored in t08,9 	*/
		t02= t02-it;			t03= t03+rt;
		rt = t06;			it = t07;
		t06= t04+it;			t07= t05-rt;	/*<==prefer these to be stored in t06,7 	*/
		t04= t04-it;			t05= t05+rt;

/*...Block 2:	*/
		t10=a[j1+p25];		t11=a[j2+p25];	/* permuted versions of p21,15,27,09,03 */
		t12=a[j1+p01];		t13=a[j2+p01];
		rt =a[j1+p19];		it =a[j2+p19];
		t16=t12-rt;			t17=t13-it;
		t12=t12+rt;			t13=t13+it;
		t14=a[j1+p07];		t15=a[j2+p07];
		rt =a[j1+p13];		it =a[j2+p13];
		t18=t14-rt;			t19=t15-it;
		t14=t14+rt;			t15=t15+it;

		rt = t12+t14;			it = t13+t15;
		t10= t10+rt;			t11= t11+it;
		rt = t10+cn1*rt;		it = t11+cn1*it;
		t14= cn2*(t12-t14);		t15= cn2*(t13-t15);
		t12= rt+t14;			t13= it+t15;
		t14= rt-t14;			t15= it-t15;
		rt = ss3*(t16-t18);		it = ss3*(t17-t19);
		t18= rt+sn1*t18;		t19= it+sn1*t19;
		t16= rt-sn2*t16;		t17= it-sn2*t17;
		rt = t18;			it = t19;
		t18= t12+it;			t19= t13-rt;
		t12= t12-it;			t13= t13+rt;
		rt = t16;			it = t17;
		t16= t14+it;			t17= t15-rt;
		t14= t14-it;			t15= t15+rt;

/*...Block 3:	*/
		t20=a[j1+p05];		t21=a[j2+p05];	/* permuted versions of p11,05,17,29,23 */
		t22=a[j1+p11];		t23=a[j2+p11];
		rt =a[j1+p29];		it =a[j2+p29];
		t26=t22-rt;			t27=t23-it;
		t22=t22+rt;			t23=t23+it;
		t24=a[j1+p17];		t25=a[j2+p17];
		rt =a[j1+p23];		it =a[j2+p23];
		t28=t24-rt;			t29=t25-it;
		t24=t24+rt;			t25=t25+it;

		rt = t22+t24;			it = t23+t25;
		t20= t20+rt;			t21= t21+it;
		rt = t20+cn1*rt;		it = t21+cn1*it;
		t24= cn2*(t22-t24);		t25= cn2*(t23-t25);
		t22= rt+t24;			t23= it+t25;
		t24= rt-t24;			t25= it-t25;
		rt = ss3*(t26-t28);		it = ss3*(t27-t29);
		t28= rt+sn1*t28;		t29= it+sn1*t29;
		t26= rt-sn2*t26;		t27= it-sn2*t27;
		rt = t28;			it = t29;
		t28= t22+it;			t29= t23-rt;
		t22= t22-it;			t23= t23+rt;
		rt = t26;			it = t27;
		t26= t24+it;			t27= t25-rt;
		t24= t24-it;			t25= t25+rt;

/*...and now do five radix-3 transforms:	*/
/*...Block 1:	*/
		rt =t20;			it =t21;
		t20=t10-rt;			t21=t11-it;
		t10=t10+rt;			t11=t11+it;
		t00=t00+t10;			t01=t01+t11;
		a1p01r=t00;			a1p01i=t01;
		t10=t00+c3m1*t10;		t11=t01+c3m1*t11;
		rt =s*t20;			it =s*t21;
		a1p03r=t10-it;		a1p03i=t11+rt;
		a1p05r=t10+it;		a1p05i=t11-rt;

/*...Block 2:	*/
		rt =t22;			it =t23;
		t22=t12-rt;			t23=t13-it;
		t12=t12+rt;			t13=t13+it;
		t02=t02+t12;			t03=t03+t13;
		a1p27r=t02;			a1p27i=t03;
		t12=t02+c3m1*t12;		t13=t03+c3m1*t13;
		rt =s*t22;			it =s*t23;
		a1p29r=t12-it;		a1p29i=t13+rt;
		a1p25r=t12+it;		a1p25i=t13-rt;

/*...Block 3:	*/
		rt =t24;			it =t25;
		t24=t14-rt;			t25=t15-it;
		t14=t14+rt;			t15=t15+it;
		t04=t04+t14;			t05=t05+t15;
		a1p19r=t04;			a1p19i=t05;
		t14=t04+c3m1*t14;		t15=t05+c3m1*t15;
		rt =s*t24;			it =s*t25;
		a1p21r=t14-it;		a1p21i=t15+rt;
		a1p23r=t14+it;		a1p23i=t15-rt;

/*...Block 4:	*/
		rt =t26;			it =t27;
		t26=t16-rt;			t27=t17-it;
		t16=t16+rt;			t17=t17+it;
		t06=t06+t16;			t07=t07+t17;
		a1p17r=t06;			a1p17i=t07;
		t16=t06+c3m1*t16;		t17=t07+c3m1*t17;
		rt =s*t26;			it =s*t27;
		a1p13r=t16-it;		a1p13i=t17+rt;
		a1p15r=t16+it;		a1p15i=t17-rt;

/*...Block 5:	*/
		rt =t28;			it =t29;
		t28=t18-rt;			t29=t19-it;
		t18=t18+rt;			t19=t19+it;
		t08=t08+t18;			t09=t09+t19;
		a1p09r=t08;			a1p09i=t09;
		t18=t08+c3m1*t18;		t19=t09+c3m1*t19;
		rt =s*t28;			it =s*t29;
		a1p11r=t18-it;		a1p11i=t19+rt;
		a1p07r=t18+it;		a1p07i=t19-rt;

/*...and now do 15 radix-2 transforms:	*/

		a[j1    ]=a1p00r+a1p01r;	a[j2    ]=a1p00i+a1p01i;
		a[j1+p01]=a1p00r-a1p01r;	a[j2+p01]=a1p00i-a1p01i;

		a[j1+p04]=a1p02r+a1p03r;	a[j2+p04]=a1p02i+a1p03i;
		a[j1+p05]=a1p02r-a1p03r;	a[j2+p05]=a1p02i-a1p03i;

		a[j1+p03]=a1p04r+a1p05r;	a[j2+p03]=a1p04i+a1p05i;
		a[j1+p02]=a1p04r-a1p05r;	a[j2+p02]=a1p04i-a1p05i;

		a[j1+p28]=a1p06r+a1p07r;	a[j2+p28]=a1p06i+a1p07i;
		a[j1+p29]=a1p06r-a1p07r;	a[j2+p29]=a1p06i-a1p07i;

		a[j1+p27]=a1p08r+a1p09r;	a[j2+p27]=a1p08i+a1p09i;
		a[j1+p26]=a1p08r-a1p09r;	a[j2+p26]=a1p08i-a1p09i;

		a[j1+p24]=a1p10r+a1p11r;	a[j2+p24]=a1p10i+a1p11i;
		a[j1+p25]=a1p10r-a1p11r;	a[j2+p25]=a1p10i-a1p11i;

		a[j1+p23]=a1p12r+a1p13r;	a[j2+p23]=a1p12i+a1p13i;
		a[j1+p22]=a1p12r-a1p13r;	a[j2+p22]=a1p12i-a1p13i;

		a[j1+p20]=a1p14r+a1p15r;	a[j2+p20]=a1p14i+a1p15i;
		a[j1+p21]=a1p14r-a1p15r;	a[j2+p21]=a1p14i-a1p15i;

		a[j1+p19]=a1p16r+a1p17r;	a[j2+p19]=a1p16i+a1p17i;
		a[j1+p18]=a1p16r-a1p17r;	a[j2+p18]=a1p16i-a1p17i;

		a[j1+p16]=a1p18r+a1p19r;	a[j2+p16]=a1p18i+a1p19i;
		a[j1+p17]=a1p18r-a1p19r;	a[j2+p17]=a1p18i-a1p19i;

		a[j1+p15]=a1p20r+a1p21r;	a[j2+p15]=a1p20i+a1p21i;
		a[j1+p14]=a1p20r-a1p21r;	a[j2+p14]=a1p20i-a1p21i;

		a[j1+p12]=a1p22r+a1p23r;	a[j2+p12]=a1p22i+a1p23i;
		a[j1+p13]=a1p22r-a1p23r;	a[j2+p13]=a1p22i-a1p23i;

		a[j1+p11]=a1p24r+a1p25r;	a[j2+p11]=a1p24i+a1p25i;
		a[j1+p10]=a1p24r-a1p25r;	a[j2+p10]=a1p24i-a1p25i;

		a[j1+p08]=a1p26r+a1p27r;	a[j2+p08]=a1p26i+a1p27i;
		a[j1+p09]=a1p26r-a1p27r;	a[j2+p09]=a1p26i-a1p27i;

		a[j1+p07]=a1p28r+a1p29r;	a[j2+p07]=a1p28i+a1p29i;
		a[j1+p06]=a1p28r-a1p29r;	a[j2+p06]=a1p28i-a1p29i;
#endif
	}
}

/***************/

void radix30_dit_pass1(double a[], int n)
{
/*
!...Acronym: DIF = Decimation In Frequency
!
!...Subroutine to perform an initial radix-30 complex DIT FFT pass on the data in the length-N real vector A.
!
!   See the documentation in radix16_dif_pass for further details on storage and indexing.
!
!   See the documentation in radix15_dif_pass for details on the radix-15 subtransforms.
*/
	int j,j1,j2;
	static int n30,p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29, first_entry=TRUE;
	static double c3m1= -1.50000000000000000000,	/* cos(twopi/3)-1	*/
					s   =  0.86602540378443864675,	/* sin(twopi/3)		*/

					cn1 = -1.25000000000000000000,	/* [cos(u)+cos(2u)]/2-1 = -5/4 */
					cn2 =  0.55901699437494742409,	/* [cos(u)-cos(2u)]/2 */
					ss3 =  0.95105651629515357210,	/*  sin(u) */
					sn1 =  1.53884176858762670130,	/* [sin(u)+sin(2u)] */
					sn2 =  0.36327126400268044292;	/* [sin(u)-sin(2u)] */
	double rt,it
	,u01,u02,u03,u04,u05,u06,u07,u08,u09,u10,u11,u12,u13,u14,u15,u16,u17,u18,u19,u20,u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,u31,u32,u33,u34,u35,u36,u37,u38,u39,u40,u41,u42,u43,u44,u45,u46,u47,u48,u49,u50,u51,u52,u53,u54,u55,u56,u57,u58,u59,u60
	,t00,t01,t02,t03,t04,t05,t06,t07,t08,t09,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29;

	if(!first_entry && (n/30) != n30)	/* New runlength?	*/
	{
		first_entry=TRUE;
	}

/*...initialize things upon first entry	*/

	if(first_entry)
	{
		first_entry=FALSE;
		n30=n/30;

/*   constant index offsets for array load/stores are here.	*/

		p01 = n30;
		p02 = p01 + n30;		p01 = p01 + ( (p01 >> DAT_BITS) << PAD_BITS );
		p03 = p02 + n30;		p02 = p02 + ( (p02 >> DAT_BITS) << PAD_BITS );
		p04 = p03 + n30;		p03 = p03 + ( (p03 >> DAT_BITS) << PAD_BITS );
		p05 = p04 + n30;		p04 = p04 + ( (p04 >> DAT_BITS) << PAD_BITS );
		p06 = p05 + n30;		p05 = p05 + ( (p05 >> DAT_BITS) << PAD_BITS );
		p07 = p06 + n30;		p06 = p06 + ( (p06 >> DAT_BITS) << PAD_BITS );
		p08 = p07 + n30;		p07 = p07 + ( (p07 >> DAT_BITS) << PAD_BITS );
		p09 = p08 + n30;		p08 = p08 + ( (p08 >> DAT_BITS) << PAD_BITS );
		p10 = p09 + n30;		p09 = p09 + ( (p09 >> DAT_BITS) << PAD_BITS );
		p11 = p10 + n30;		p10 = p10 + ( (p10 >> DAT_BITS) << PAD_BITS );
		p12 = p11 + n30;		p11 = p11 + ( (p11 >> DAT_BITS) << PAD_BITS );
		p13 = p12 + n30;		p12 = p12 + ( (p12 >> DAT_BITS) << PAD_BITS );
		p14 = p13 + n30;		p13 = p13 + ( (p13 >> DAT_BITS) << PAD_BITS );
		p15 = p14 + n30;		p14 = p14 + ( (p14 >> DAT_BITS) << PAD_BITS );
		p16 = p15 + n30;		p15 = p15 + ( (p15 >> DAT_BITS) << PAD_BITS );
		p17 = p16 + n30;		p16 = p16 + ( (p16 >> DAT_BITS) << PAD_BITS );
		p18 = p17 + n30;		p17 = p17 + ( (p17 >> DAT_BITS) << PAD_BITS );
		p19 = p18 + n30;		p18 = p18 + ( (p18 >> DAT_BITS) << PAD_BITS );
		p20 = p19 + n30;		p19 = p19 + ( (p19 >> DAT_BITS) << PAD_BITS );
		p21 = p20 + n30;		p20 = p20 + ( (p20 >> DAT_BITS) << PAD_BITS );
		p22 = p21 + n30;		p21 = p21 + ( (p21 >> DAT_BITS) << PAD_BITS );
		p23 = p22 + n30;		p22 = p22 + ( (p22 >> DAT_BITS) << PAD_BITS );
		p24 = p23 + n30;		p23 = p23 + ( (p23 >> DAT_BITS) << PAD_BITS );
		p25 = p24 + n30;		p24 = p24 + ( (p24 >> DAT_BITS) << PAD_BITS );
		p26 = p25 + n30;		p25 = p25 + ( (p25 >> DAT_BITS) << PAD_BITS );
		p27 = p26 + n30;		p26 = p26 + ( (p26 >> DAT_BITS) << PAD_BITS );
		p28 = p27 + n30;		p27 = p27 + ( (p27 >> DAT_BITS) << PAD_BITS );
		p29 = p28 + n30;		p28 = p28 + ( (p28 >> DAT_BITS) << PAD_BITS );
								p29 = p29 + ( (p29 >> DAT_BITS) << PAD_BITS );
	}

/*...The radix-30 pass is here.	*/

	for(j=0; j < n30; j += 2)
	{
	#ifdef USE_SSE2
		j1 = (j & mask01) + br4[j&3];
		j1 =j1 + ( (j1>> DAT_BITS) << PAD_BITS );
	#else
		j1 = j + ( (j >> DAT_BITS) << PAD_BITS );	/* padded-array fetch index is here */
	#endif
		j2 = j1+RE_IM_STRIDE;

	/*...gather the needed data (30 64-bit complex, i.e. 60 64-bit reals) and do 15 radix-2 transforms,	*/
	/*
		Twiddleless version requires us to swap inputs as follows:
		indices  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29
			  => 0,28,26,24,22,20,18,16,14,12,10, 8, 6, 4, 2,15,13,11, 9, 7, 5, 3, 1,29,27,25,23,21,19,17 modulo 30.
		I.e. start out with first pair of indices {0,15}, permute those according to
		{0,15}*29%30 = {0,15}, then each is head of a length-15 list of indices with decrement 2.

		Remember, inputs to DIT are bit-reversed, so
		a[0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29] contain
		x[0,15, 5,20,10,25, 1,16, 6,21,11,26, 2,17, 7,22,12,27, 3,18, 8,23,13,28, 4,19, 9,24,14,29], which get swapped to
		x[0,15,20, 5,10,25,28,13,18, 3, 8,23,26,11,16, 1, 6,21,24, 9,14,29, 4,19,22, 7,12,27, 2,17], which means the a-indices get swapped as
		a[0, 1, 3, 2, 4, 5,23,22,19,18,20,21,11,10, 7, 6, 8, 9,27,26,28,29,24,25,15,14,16,17,12,13].
	*/
		u01=a[j1      ];		u02=a[j2    ];	/* x0a,b */
		rt =a[j1  +p01];		it =a[j2+p01];
		u03=u01-rt;				u04=u02-it;
		u01=u01+rt;				u02=u02+it;

		u05=a[j1  +p03];		u06=a[j2+p03];	/* x1a,b */
		rt =a[j1  +p02];		it =a[j2+p02];
		u07=u05-rt;				u08=u06-it;
		u05=u05+rt;				u06=u06+it;

		u09=a[j1  +p04];		u10=a[j2+p04];	/* x2a,b */
		rt =a[j1  +p05];		it =a[j2+p05];
		u11=u09-rt;				u12=u10-it;
		u09=u09+rt;				u10=u10+it;

		u13=a[j1  +p23];		u14=a[j2+p23];	/* x3a,b */
		rt =a[j1  +p22];		it =a[j2+p22];
		u15=u13-rt;				u16=u14-it;
		u13=u13+rt;				u14=u14+it;

		u17=a[j1  +p19];		u18=a[j2+p19];	/* x4a,b */
		rt =a[j1  +p18];		it =a[j2+p18];
		u19=u17-rt;				u20=u18-it;
		u17=u17+rt;				u18=u18+it;

		u21=a[j1  +p20];		u22=a[j2+p20];	/* x5a,b */
		rt =a[j1  +p21];		it =a[j2+p21];
		u23=u21-rt;				u24=u22-it;
		u21=u21+rt;				u22=u22+it;

		u25=a[j1  +p11];		u26=a[j2+p11];	/* x6a,b */
		rt =a[j1  +p10];		it =a[j2+p10];
		u27=u25-rt;				u28=u26-it;
		u25=u25+rt;				u26=u26+it;

		u29=a[j1  +p07];		u30=a[j2+p07];	/* x7a,b */
		rt =a[j1  +p06];		it =a[j2+p06];
		u31=u29-rt;				u32=u30-it;
		u29=u29+rt;				u30=u30+it;

		u33=a[j1  +p08];		u34=a[j2+p08];	/* x8a,b */
		rt =a[j1  +p09];		it =a[j2+p09];
		u35=u33-rt;				u36=u34-it;
		u33=u33+rt;				u34=u34+it;

		u37=a[j1  +p27];		u38=a[j2+p27];	/* x9a,b */
		rt =a[j1  +p26];		it =a[j2+p26];
		u39=u37-rt;				u40=u38-it;
		u37=u37+rt;				u38=u38+it;

		u41=a[j1  +p28];		u42=a[j2+p28];	/* x10a,b */
		rt =a[j1  +p29];		it =a[j2+p29];
		u43=u41-rt;				u44=u42-it;
		u41=u41+rt;				u42=u42+it;

		u45=a[j1  +p24];		u46=a[j2+p24];	/* x11a,b */
		rt =a[j1  +p25];		it =a[j2+p25];
		u47=u45-rt;				u48=u46-it;
		u45=u45+rt;				u46=u46+it;

		u49=a[j1  +p15];		u50=a[j2+p15];	/* x12a,b */
		rt =a[j1  +p14];		it =a[j2+p14];
		u51=u49-rt;				u52=u50-it;
		u49=u49+rt;				u50=u50+it;

		u53=a[j1  +p16];		u54=a[j2+p16];	/* x13a,b */
		rt =a[j1  +p17];		it =a[j2+p17];
		u55=u53-rt;				u56=u54-it;
		u53=u53+rt;				u54=u54+it;

		u57=a[j1  +p12];		u58=a[j2+p12];	/* x14a,b */
		rt =a[j1  +p13];		it =a[j2+p13];
		u59=u57-rt;				u60=u58-it;
		u57=u57+rt;				u58=u58+it;

/*       ...and now do two radix-15 transforms.	*/

	/*...a1p[0:28:2]r use u[1:57:4]; a1p[0:28:2]i use u[2:58:4]	*/

	/*   Gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do five radix-3 transforms...	*/
	/*...Block 1:	*/
		t00=u01;			t01=u02;
		t02=u09;			t03=u10;
		rt =u05;			it =u06;
		t04=t02-rt;			t05=t03-it;
		t02=t02+rt;			t03=t03+it;
		t00=t00+t02;		t01=t01+t03;
		t02=t00+c3m1*t02;		t03=t01+c3m1*t03;
		rt =s*t04;			it =s*t05;
		t04=t02-it;			t05=t03+rt;
		t02=t02+it;			t03=t03-rt;

	/*...Block 2:	*/
		t06=u33;			t07=u34;
		t08=u29;			t09=u30;
		rt =u25;			it =u26;
		t10=t08-rt;			t11=t09-it;
		t08=t08+rt;			t09=t09+it;
		t06=t06+t08;		t07=t07+t09;
		t08=t06+c3m1*t08;		t09=t07+c3m1*t09;
		rt =s*t10;			it =s*t11;
		t10=t08-it;			t11=t09+rt;
		t08=t08+it;			t09=t09-rt;

	/*...Block 3:	*/
		t12=u53;			t13=u54;
		t14=u49;			t15=u50;
		rt =u57;			it =u58;
		t16=t14 -rt;		t17=t15 -it;
		t14=t14 +rt;		t15=t15 +it;
		t12=t12+t14;		t13=t13+t15;
		t14=t12+c3m1*t14;		t15=t13+c3m1*t15;
		rt =s*t16;			it =s*t17;
		t16=t14-it;			t17=t15+rt;
		t14=t14+it;			t15=t15-rt;

	/*...Block 4:	*/
		t18=u17;			t19=u18;
		t20=u13;			t21=u14;
		rt =u21;			it =u22;
		t22=t20 -rt;		t23=t21 -it;
		t20=t20 +rt;		t21=t21 +it;
		t18=t18+t20;		t19=t19+t21;
		t20=t18+c3m1*t20;		t21=t19+c3m1*t21;
		rt =s*t22;			it =s*t23;
		t22=t20-it;			t23=t21+rt;
		t20=t20+it;			t21=t21-rt;

	/*...Block 5:	*/
		t24=u37;			t25=u38;
		t26=u45;			t27=u46;
		rt =u41;			it =u42;
		t28=t26 -rt;		t29=t27 -it;
		t26=t26 +rt;		t27=t27 +it;
		t24=t24+t26;		t25=t25+t27;
		t26=t24+c3m1*t26;		t27=t25+c3m1*t27;
		rt =s*t28;			it =s*t29;
		t28=t26-it;			t29=t27+rt;
		t26=t26+it;			t27=t27-rt;

	/*...and now do three radix-5 transforms:	*/
	/*...Block 1:	*/
		rt = t24;			it = t25;
		t24= t06-rt;		t25= t07-it;
		t06= t06+rt;		t07= t07+it;
		rt = t18;			it = t19;
		t18= t12-rt;		t19= t13-it;
		t12= t12+rt;		t13= t13+it;

		rt = t06+t12;		it = t07+t13;
		t00= t00+rt;		t01= t01+it;
		rt = t00+cn1*rt;		it = t01+cn1*it;
		t12= cn2*(t06-t12);		t13= cn2*(t07-t13);
		t06= rt+t12;		t07= it+t13;
		t12= rt-t12;		t13= it-t13;
		rt = ss3*(t18-t24);		it = ss3*(t19-t25);
		t18= rt-sn1*t18;		t19= it-sn1*t19;
		t24= rt+sn2*t24;		t25= it+sn2*t25;

		a[j1    ]=t00;		a[j2    ]=t01;
		a[j1+p06]=t06-t19;		a[j2+p06]=t07+t18;
		a[j1+p12]=t12-t25;		a[j2+p12]=t13+t24;
		a[j1+p18]=t12+t25;		a[j2+p18]=t13-t24;
		a[j1+p24]=t06+t19;		a[j2+p24]=t07-t18;

	/*...Block 2:	*/
		rt = t26;			it = t27;
		t26= t08-rt;		t27= t09-it;
		t08= t08+rt;		t09= t09+it;
		rt = t20;			it = t21;
		t20= t14-rt;		t21= t15-it;
		t14= t14+rt;		t15= t15+it;

		rt = t08+t14;		it = t09+t15;
		t02= t02+rt;		t03= t03+it;
		rt = t02+cn1*rt;		it = t03+cn1*it;
		t14= cn2*(t08-t14);		t15= cn2*(t09-t15);
		t08= rt+t14;		t09= it+t15;
		t14= rt-t14;		t15= it-t15;
		rt = ss3*(t20-t26);		it = ss3*(t21-t27);
		t20= rt-sn1*t20;		t21= it-sn1*t21;
		t26= rt+sn2*t26;		t27= it+sn2*t27;

		a[j1+p10]=t02;		a[j2+p10]=t03;
		a[j1+p16]=t08-t21;		a[j2+p16]=t09+t20;
		a[j1+p22]=t14-t27;		a[j2+p22]=t15+t26;
		a[j1+p28]=t14+t27;		a[j2+p28]=t15-t26;
		a[j1+p04]=t08+t21;		a[j2+p04]=t09-t20;

	/*...Block 3:	*/
		rt = t28;			it = t29;
		t28= t10-rt;		t29= t11-it;
		t10= t10+rt;		t11= t11+it;
		rt = t22;			it = t23;
		t22= t16-rt;		t23= t17-it;
		t16= t16+rt;		t17= t17+it;

		rt = t10+t16;		it = t11+t17;
		t04= t04+rt;		t05= t05+it;
		rt = t04+cn1*rt;		it = t05+cn1*it;
		t16= cn2*(t10-t16);		t17= cn2*(t11-t17);
		t10= rt+t16;		t11= it+t17;
		t16= rt-t16;		t17= it-t17;
		rt = ss3*(t22-t28);		it = ss3*(t23-t29);
		t22= rt-sn1*t22;		t23= it-sn1*t23;
		t28= rt+sn2*t28;		t29= it+sn2*t29;

		a[j1+p20]=t04;		a[j2+p20]=t05;
		a[j1+p26]=t10-t23;		a[j2+p26]=t11+t22;
		a[j1+p02]=t16-t29;		a[j2+p02]=t17+t28;
		a[j1+p08]=t16+t29;		a[j2+p08]=t17-t28;
		a[j1+p14]=t10+t23;		a[j2+p14]=t11-t22;

	/*...a1p[1:29:2]r use u[3:59:4]; a1p[1:29:2]i use u[4:60:4]	*/

	/*   Gather the needed data (15 64-bit complex, i.e. 30 64-bit reals) and do five radix-3 transforms...	*/
	/*...Block 1:	*/
		t00=u03;			t01=u04;
		t02=u11;			t03=u12;
		rt =u07;			it =u08;
		t04=t02-rt;			t05=t03-it;
		t02=t02+rt;			t03=t03+it;
		t00=t00+t02;		t01=t01+t03;
		t02=t00+c3m1*t02;		t03=t01+c3m1*t03;
		rt =s*t04;			it =s*t05;
		t04=t02-it;			t05=t03+rt;
		t02=t02+it;			t03=t03-rt;

	/*...Block 2:	*/
		t06=u35;			t07=u36;
		t08=u31;			t09=u32;
		rt =u27;			it =u28;
		t10=t08-rt;			t11=t09-it;
		t08=t08+rt;			t09=t09+it;
		t06=t06+t08;		t07=t07+t09;
		t08=t06+c3m1*t08;		t09=t07+c3m1*t09;
		rt =s*t10;			it =s*t11;
		t10=t08-it;			t11=t09+rt;
		t08=t08+it;			t09=t09-rt;

	/*...Block 3:	*/
		t12=u55;			t13=u56;
		t14=u51;			t15=u52;
		rt =u59;			it =u60;
		t16=t14 -rt;		t17=t15 -it;
		t14=t14 +rt;		t15=t15 +it;
		t12=t12+t14;		t13=t13+t15;
		t14=t12+c3m1*t14;		t15=t13+c3m1*t15;
		rt =s*t16;			it =s*t17;
		t16=t14-it;			t17=t15+rt;
		t14=t14+it;			t15=t15-rt;

	/*...Block 4:	*/
		t18=u19;			t19=u20;
		t20=u15;			t21=u16;
		rt =u23;			it =u24;
		t22=t20 -rt;		t23=t21 -it;
		t20=t20 +rt;		t21=t21 +it;
		t18=t18+t20;		t19=t19+t21;
		t20=t18+c3m1*t20;		t21=t19+c3m1*t21;
		rt =s*t22;			it =s*t23;
		t22=t20-it;			t23=t21+rt;
		t20=t20+it;			t21=t21-rt;

	/*...Block 5:	*/
		t24=u39;			t25=u40;
		t26=u47;			t27=u48;
		rt =u43;			it =u44;
		t28=t26 -rt;		t29=t27 -it;
		t26=t26 +rt;		t27=t27 +it;
		t24=t24+t26;		t25=t25+t27;
		t26=t24+c3m1*t26;		t27=t25+c3m1*t27;
		rt =s*t28;			it =s*t29;
		t28=t26-it;			t29=t27+rt;
		t26=t26+it;			t27=t27-rt;

	/*...and now do three radix-5 transforms:	*/
	/*...Block 1:	*/
		rt = t24;			it = t25;
		t24= t06-rt;		t25= t07-it;
		t06= t06+rt;		t07= t07+it;
		rt = t18;			it = t19;
		t18= t12-rt;		t19= t13-it;
		t12= t12+rt;		t13= t13+it;

		rt = t06+t12;		it = t07+t13;
		t00= t00+rt;		t01= t01+it;
		rt = t00+cn1*rt;		it = t01+cn1*it;
		t12= cn2*(t06-t12);		t13= cn2*(t07-t13);
		t06= rt+t12;		t07= it+t13;
		t12= rt-t12;		t13= it-t13;
		rt = ss3*(t18-t24);		it = ss3*(t19-t25);
		t18= rt-sn1*t18;		t19= it-sn1*t19;
		t24= rt+sn2*t24;		t25= it+sn2*t25;

		a[j1+p15]=t00;		a[j2+p15]=t01;
		a[j1+p21]=t06-t19;		a[j2+p21]=t07+t18;
		a[j1+p27]=t12-t25;		a[j2+p27]=t13+t24;
		a[j1+p03]=t12+t25;		a[j2+p03]=t13-t24;
		a[j1+p09]=t06+t19;		a[j2+p09]=t07-t18;

	/*...Block 2:	*/
		rt = t26;			it = t27;
		t26= t08-rt;		t27= t09-it;
		t08= t08+rt;		t09= t09+it;
		rt = t20;			it = t21;
		t20= t14-rt;		t21= t15-it;
		t14= t14+rt;		t15= t15+it;

		rt = t08+t14;		it = t09+t15;
		t02= t02+rt;		t03= t03+it;
		rt = t02+cn1*rt;		it = t03+cn1*it;
		t14= cn2*(t08-t14);		t15= cn2*(t09-t15);
		t08= rt+t14;		t09= it+t15;
		t14= rt-t14;		t15= it-t15;
		rt = ss3*(t20-t26);		it = ss3*(t21-t27);
		t20= rt-sn1*t20;		t21= it-sn1*t21;
		t26= rt+sn2*t26;		t27= it+sn2*t27;

		a[j1+p25]=t02;		a[j2+p25]=t03;
		a[j1+p01]=t08-t21;		a[j2+p01]=t09+t20;
		a[j1+p07]=t14-t27;		a[j2+p07]=t15+t26;
		a[j1+p13]=t14+t27;		a[j2+p13]=t15-t26;
		a[j1+p19]=t08+t21;		a[j2+p19]=t09-t20;

	/*...Block 3:	*/
		rt = t28;			it = t29;
		t28= t10-rt;		t29= t11-it;
		t10= t10+rt;		t11= t11+it;
		rt = t22;			it = t23;
		t22= t16-rt;		t23= t17-it;
		t16= t16+rt;		t17= t17+it;

		rt = t10+t16;		it = t11+t17;
		t04= t04+rt;		t05= t05+it;
		rt = t04+cn1*rt;		it = t05+cn1*it;
		t16= cn2*(t10-t16);		t17= cn2*(t11-t17);
		t10= rt+t16;		t11= it+t17;
		t16= rt-t16;		t17= it-t17;
		rt = ss3*(t22-t28);		it = ss3*(t23-t29);
		t22= rt-sn1*t22;		t23= it-sn1*t23;
		t28= rt+sn2*t28;		t29= it+sn2*t29;

		a[j1+p05]=t04;		a[j2+p05]=t05;
		a[j1+p11]=t10-t23;		a[j2+p11]=t11+t22;
		a[j1+p17]=t16-t29;		a[j2+p17]=t17+t28;
		a[j1+p23]=t16+t29;		a[j2+p23]=t17-t28;
		a[j1+p29]=t10+t23;		a[j2+p29]=t11-t22;

	}
}

